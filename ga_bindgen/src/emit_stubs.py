#!/usr/bin/env python3
"""Generate type stubs (.pyi) for ga_py from the compiled extension.

The pipeline runs `nanobind.stubgen` against the loaded `_ga_py` module
to capture the bindings' Python surface, then post-processes the output
in two ways before writing it into the wrapper package
(`ga_py/python/ga_py/`):

  1. Inject `__format__(self, format_spec: str = "") -> str` into every
     bound class. nanobind.stubgen filters this dunder out of its output,
     but we define it explicitly on every type via `emit_format_lambda`
     in emit_nanobind.py and the hand-written physics bindings — so it
     belongs in the stub.

  2. Append signatures for the runtime forwarders in
     `ga_py/python/ga_py/__init__.py` (pure-Python wrappers injected
     into the `pga` submodule at import time, e.g. `pga.ortho_proj3dp`,
     `pga.dist2dp`, the `*_contract*` / `*_expand*` families). nanobind
     can't see these because they're added at runtime, not at the
     binding site.

After post-processing, the script writes:
  * ga_py/python/ga_py/__init__.pyi    (top-level + exports)
  * ga_py/python/ga_py/ega.pyi          (Euclidean GA surface)
  * ga_py/python/ga_py/pga.pyi          (Projective GA surface + forwarders)
  * ga_py/python/ga_py/py.typed         (PEP 561 marker)

Forwarders are typed loosely as `(a: Any, b: Any) -> Any` since they are
generic compositions of bound C++ functions whose return type depends on
the (a, b) input combination — the same reason they live in
`__init__.py` rather than the bindgen output. Tightening the signatures
is a future refinement; the current shape unblocks IDE completion and
prevents `pga.<forwarder>` from triggering "missing attribute" errors.

Run via the CMake target `regenerate_python_stubs` (see ga_py/CMakeLists.txt),
or directly from the project root:

    PYTHONPATH=build/ga_py:ga_py/python \\
        ga_py/.venv/bin/python ga_bindgen/src/emit_stubs.py
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths.
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TARGET_DIR = PROJECT_ROOT / "ga_py" / "python" / "ga_py"
INIT_PY = TARGET_DIR / "__init__.py"

STUB_FILES = ("__init__.pyi", "ega.pyi", "pga.pyi")


# ---------------------------------------------------------------------------
# Run nanobind.stubgen against the loaded `_ga_py` module.
# ---------------------------------------------------------------------------

def run_stubgen(out_dir: Path) -> None:
    """Invoke `python -m nanobind.stubgen` as a subprocess targeting the
    compiled extension. Subprocess isolation keeps any sys.path mutations
    inside stubgen from leaking into this script."""
    cmd = [sys.executable, "-m", "nanobind.stubgen",
           "-m", "_ga_py", "-r", "-O", str(out_dir)]
    subprocess.run(cmd, check=True, cwd=str(PROJECT_ROOT))


# ---------------------------------------------------------------------------
# Post-pass 1 — inject __format__ into every class block.
# ---------------------------------------------------------------------------

# All bound classes declare both __repr__ and __str__ (emit_format_lambda
# in emit_nanobind.py and the hand-written bindings_*.cpp). We anchor the
# injection on `def __str__(self) -> str: ...` because it sits at a fixed
# indentation (4 spaces) and is unique per class.
_STR_LINE_RE = re.compile(r"(\n    def __str__\(self\) -> str: \.\.\.\n)",
                          re.MULTILINE)
_FORMAT_LINE = ('\n    def __format__(self, format_spec: str = "") -> str: '
                "...\n")


def inject_format(text: str) -> str:
    return _STR_LINE_RE.sub(r"\1" + _FORMAT_LINE, text)


# ---------------------------------------------------------------------------
# Post-pass 2 — discover forwarders in ga_py/__init__.py and append
# matching `def NAME(a, b) -> Any: ...` lines to the relevant submodule
# stub.
# ---------------------------------------------------------------------------

# Matches lines like `pga.ortho_proj2dp = _pga2dp_ortho_proj2dp` — the
# wrapper assignments in ga_py/python/ga_py/__init__.py that inject
# pure-Python helpers into the bound submodule at import time.
_FORWARDER_RE = re.compile(r"^([a-z]+)\.(\w+)\s*=\s*\w+\s*$", re.MULTILINE)


def extract_forwarders() -> dict[str, list[str]]:
    """Return {submodule: [func_name, ...]} for every `<submod>.NAME = ...`
    assignment in __init__.py. Order of appearance is preserved; duplicates
    (last-write-wins) are deduplicated."""
    text = INIT_PY.read_text()
    found: dict[str, list[str]] = {}
    for submod, func in _FORWARDER_RE.findall(text):
        bucket = found.setdefault(submod, [])
        if func not in bucket:
            bucket.append(func)
    return found


def append_forwarder_stubs(text: str, names: list[str]) -> str:
    if not names:
        return text
    lines = ["", "",
             "# Pure-Python forwarders injected at import time from "
             "ga_py/__init__.py",
             "# (typed loosely — return type depends on the (a, b) "
             "combination at the call site)",
             ""]
    for name in names:
        lines.append(f"def {name}(a: Any, b: Any) -> Any: ...")
    return text.rstrip() + "\n" + "\n".join(lines) + "\n"


def ensure_any_import(text: str) -> str:
    """Make sure `from typing import Any` is present in the stub. Stubgen
    typically imports `overload` from `typing`; we widen that line if
    needed."""
    if re.search(r"^from typing import .*\bAny\b", text, re.MULTILINE):
        return text
    # Try to extend an existing `from typing import ...` line.
    pat = re.compile(r"^(from typing import )(.+)$", re.MULTILINE)
    m = pat.search(text)
    if m:
        existing = [s.strip() for s in m.group(2).split(",")]
        if "Any" not in existing:
            existing.append("Any")
        return pat.sub(m.group(1) + ", ".join(sorted(existing)), text, count=1)
    # No typing import yet — add one near the top, after the module
    # docstring and any other `from ... import ...` lines.
    return "from typing import Any\n\n" + text


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="Generate type stubs (.pyi) for ga_py.")
    p.add_argument("--out-dir", default=str(TARGET_DIR),
                   help="Output directory (default: %(default)s).")
    args = p.parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        run_stubgen(Path(td))
        src_pkg = Path(td) / "_ga_py"  # stubgen writes the C-ext name
        forwarders = extract_forwarders()
        for stub_name in STUB_FILES:
            text = (src_pkg / stub_name).read_text()
            text = inject_format(text)
            submod = stub_name.removesuffix(".pyi")  # __init__ / ega / pga
            if submod in forwarders:
                text = ensure_any_import(text)
                text = append_forwarder_stubs(text, forwarders[submod])
            (out_dir / stub_name).write_text(text)

    # PEP 561 marker so mypy / pyright honour the inline stubs.
    (out_dir / "py.typed").write_text("")

    print(f"Stubs written to {out_dir}/")
    print("  - " + ", ".join(STUB_FILES) + ", py.typed")
    if forwarders:
        for submod, names in forwarders.items():
            print(f"  - {submod}.pyi: appended {len(names)} forwarder(s)")


if __name__ == "__main__":
    main()
