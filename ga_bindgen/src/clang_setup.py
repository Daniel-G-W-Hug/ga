"""Locate libclang and assemble compile flags for parsing the `ga/` headers.

The pip-installed `libclang` package ships a `.dylib`/`.so`/`.dll` and a
minimal clang resource directory. On macOS and Linux the bundled resource dir
is not complete enough to parse std headers (stdarg.h, __type_traits, …), so
we probe for a system LLVM install that carries a fuller resource dir. On
Windows the pip bundle is sufficient — no separate LLVM install is needed.
"""
from __future__ import annotations

import platform
import shutil
import subprocess
import sys
from pathlib import Path

import clang.cindex as cx

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent  # repo root
GA_DIR = PROJECT_ROOT / "ga"


def find_libclang() -> Path:
    """Return path to a libclang shared library.

    macOS / Linux: prefer a system or homebrew LLVM whose resource dir sits
    next to the dylib/so. Falls back to whatever `clang.cindex` defaults to.

    Windows: the pip `libclang` package bundles a usable `libclang.dll` inside
    `clang/native/` — no separate LLVM install required. A system LLVM install
    at `C:/Program Files/LLVM/` is tried first, falling back to the bundled DLL.
    """
    candidates: list[Path] = []
    system = platform.system()

    if system == "Darwin":
        candidates += [
            Path("/opt/homebrew/opt/llvm/lib/libclang.dylib"),
            Path("/usr/local/opt/llvm/lib/libclang.dylib"),
        ]
    elif system == "Linux":
        for v in ("18", "17", "16", "15", "14"):
            candidates += [
                Path(f"/usr/lib/x86_64-linux-gnu/libclang-{v}.so.1"),
                Path(f"/usr/lib/llvm-{v}/lib/libclang.so.1"),
            ]
        candidates += [Path("/usr/lib/x86_64-linux-gnu/libclang.so")]
    elif system == "Windows":
        # Standard LLVM installer locations (optional, better resource dir)
        candidates += [
            Path("C:/Program Files/LLVM/lib/libclang.dll"),
            Path("C:/Program Files (x86)/LLVM/lib/libclang.dll"),
        ]
        # Bundled DLL from the pip `libclang` package — works out of the box
        try:
            import clang.native as _cn
            bundled = Path(_cn.__file__).parent / "libclang.dll"
            candidates.append(bundled)
        except (ImportError, AttributeError):
            pass

    for c in candidates:
        if c.exists():
            return c

    raise RuntimeError(
        "Could not locate a usable libclang.\n"
        "  macOS:  brew install llvm\n"
        "  Linux:  apt install libclang-dev  (or distro equivalent)\n"
        "  Windows: pip install libclang  (inside the ga_bindgen venv)\n"
        "See ga_bindgen/README.md §Platform notes for the full candidate list."
    )


def discover_include_paths(libclang_path: Path) -> list[str]:
    """Use the matching `clang` driver to discover include paths.

    Runs `clang -E -v` to get the same header search list the driver uses.
    Returns an empty list if the clang binary cannot be found.

    Not called on Windows — use discover_msvc_include_paths() instead.
    """
    clang_bin: Path | None = None
    for cand in [
        libclang_path.parent.parent / "bin" / "clang",
        libclang_path.parent.parent / "bin" / "clang++",
        Path(shutil.which("clang") or ""),
    ]:
        if cand and cand.exists():
            clang_bin = cand
            break

    if clang_bin is None:
        return []

    null_dev = "NUL" if platform.system() == "Windows" else "/dev/null"
    proc = subprocess.run(
        [str(clang_bin), "-E", "-v", "-x", "c++", "-std=c++23", null_dev],
        capture_output=True, text=True,
    )
    out = proc.stderr
    in_block = False
    paths: list[str] = []
    for line in out.splitlines():
        if "#include <...> search starts here:" in line:
            in_block = True
            continue
        if "End of search list" in line:
            break
        if in_block:
            stripped = line.strip()
            if stripped and "(framework directory)" not in stripped:
                paths.append(stripped)
    return paths


def discover_msvc_include_paths() -> list[str]:
    """Find MSVC and Windows SDK include directories on Windows.

    Uses vswhere.exe to locate the latest Visual Studio installation, then
    resolves the MSVC toolchain and Windows SDK include directories.
    Returns an empty list if Visual Studio is not installed.
    """
    paths: list[str] = []
    vswhere = Path("C:/Program Files (x86)/Microsoft Visual Studio/Installer/vswhere.exe")
    if not vswhere.exists():
        return paths

    proc = subprocess.run(
        [str(vswhere), "-latest", "-property", "installationPath"],
        capture_output=True, text=True,
    )
    vs_path = proc.stdout.strip()
    if not vs_path:
        return paths

    # MSVC STL headers
    msvc_base = Path(vs_path) / "VC" / "Tools" / "MSVC"
    if msvc_base.exists():
        msvc_versions = sorted(msvc_base.iterdir())
        if msvc_versions:
            msvc_inc = msvc_versions[-1] / "include"
            if msvc_inc.exists():
                paths.append(str(msvc_inc))

    # Windows SDK (ucrt, shared, um)
    try:
        import winreg
        key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                             r"SOFTWARE\Microsoft\Windows Kits\Installed Roots")
        kits_root, _ = winreg.QueryValueEx(key, "KitsRoot10")
        winreg.CloseKey(key)
    except (ImportError, OSError):
        kits_root = None

    if kits_root:
        inc_root = Path(kits_root) / "Include"
        if inc_root.exists():
            sdk_versions = sorted(inc_root.iterdir())
            if sdk_versions:
                sdk_inc = sdk_versions[-1]
                for sub in ("ucrt", "shared", "um"):
                    p = sdk_inc / sub
                    if p.exists():
                        paths.append(str(p))

    return paths


def compile_args() -> list[str]:
    """Compile flags suitable for parsing `ga/*.hpp` with libclang."""
    libclang_path = find_libclang()
    cx.Config.set_library_file(str(libclang_path))

    args: list[str] = [
        "-x", "c++",
        "-std=c++23",
        "-DFMT_HEADER_ONLY=1",
        "-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO",
        f"-I{PROJECT_ROOT}",
    ]

    system = platform.system()
    if system == "Windows":
        # Allow using the pip-bundled libclang 18 with MSVC STL 2022+ which
        # requires Clang 19+. The define suppresses the compile-time assertion.
        args.append("-D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH")
        for p in discover_msvc_include_paths():
            args += ["-isystem", p]
        # fmt headers via vcpkg (the standard location in this project)
        vcpkg_inc = Path("C:/vcpkg/installed/x64-windows/include")
        if vcpkg_inc.exists():
            args += ["-isystem", str(vcpkg_inc)]
    else:
        args += ["-isystem", "/opt/homebrew/include"]  # fmt on macOS
        for p in discover_include_paths(libclang_path):
            args += ["-isystem", p]

    if system == "Darwin":
        for sdk_root in (
            "/Library/Developer/CommandLineTools/SDKs/MacOSX26.sdk",
            "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
        ):
            if Path(sdk_root).exists():
                args += ["-isysroot", sdk_root]
                break

    return args


def in_ga_tree(location) -> bool:
    """Filter: True if cursor location is inside ga/."""
    if location.file is None:
        return False
    try:
        Path(str(location.file)).relative_to(GA_DIR)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    lc = find_libclang()
    print(f"libclang: {lc}")
    args = compile_args()
    print("compile args:")
    for a in args:
        print(f"  {a}")
