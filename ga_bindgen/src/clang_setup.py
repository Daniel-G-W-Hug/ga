"""Locate libclang and assemble compile flags for parsing the `ga/` headers.

The pip-installed `libclang` package ships a `.dylib`/`.so` but no clang
resource directory, which silently breaks parsing of std headers (stdarg.h,
__type_traits, …). We therefore probe for a usable system / homebrew libclang
and assemble compile flags that mirror what `clang++ -E -v` reports.
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
    """Return path to a libclang shared library whose resource dir is reachable.

    macOS preference: homebrew LLVM (its resource dir sits next to the dylib).
    Linux preference: system /usr/lib*/libclang*.so.
    Falls back to whatever `clang.cindex` defaults to.
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
        # On Windows the LLVM installer typically puts libclang.dll on PATH.
        for v in ("18", "17", "16"):
            candidates += [
                Path(f"C:/Program Files/LLVM/lib/libclang.dll"),
                Path(f"C:/Program Files (x86)/LLVM/lib/libclang.dll"),
            ]

    for c in candidates:
        if c.exists():
            return c

    # Last resort: hope the bundled libclang's resource dir is reachable.
    raise RuntimeError(
        "Could not locate a usable libclang. Install LLVM (e.g. `brew install llvm` "
        "on macOS, `apt install libclang-dev` on Debian/Ubuntu) and re-run."
    )


def discover_include_paths(libclang_path: Path) -> list[str]:
    """Use the matching `clang` driver (sibling of libclang) to discover the
    libc++ + builtin + SDK include paths the way the driver itself would."""
    # Find the clang executable next to libclang, e.g. /opt/homebrew/opt/llvm/bin/clang
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

    proc = subprocess.run(
        [str(clang_bin), "-E", "-v", "-x", "c++", "-std=c++23", "/dev/null"],
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
        "-isystem", "/opt/homebrew/include",  # fmt on macOS
    ]
    for p in discover_include_paths(libclang_path):
        args += ["-isystem", p]

    # Add sysroot on macOS (if using CommandLineTools SDK).
    if platform.system() == "Darwin":
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
    return str(location.file).startswith(str(GA_DIR) + "/")


if __name__ == "__main__":
    lc = find_libclang()
    print(f"libclang: {lc}")
    args = compile_args()
    print("compile args:")
    for a in args:
        print(f"  {a}")
