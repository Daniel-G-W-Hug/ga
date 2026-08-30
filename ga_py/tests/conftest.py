"""Shared test fixtures and helpers for the ga_py test suite.

This file is loaded automatically by pytest. It provides:

  * Discovery of the built extension + the Python package, so no PYTHONPATH is needed.
  * A clear error if the test session was launched from the wrong venv.
  * Per-algebra / per-operation tolerance fixtures.
  * An `approx_eq()` helper that mirrors the C++ `is_congruent` relative-
    tolerance pattern (see CLAUDE.md "Critical Numerical Precision Lessons").

Run from the wrapper venv, from the project root:
    source ga_py/.venv/bin/activate
    pytest ga_py/tests/

Setting PYTHONPATH is no longer required --- the bootstrap below locates the build
tree and prepends it, which also means it takes precedence over any PYTHONPATH
already pointing at a different checkout.

Do NOT use ga_bindgen/.venv for tests --- that one is for ga_bindgen.
See ga_py/README.md > "Two virtual environments".
"""

from __future__ import annotations

import importlib.machinery
import sys
from pathlib import Path

import pytest

# --------------------------------------------------------------------------- #
# Import bootstrap
# --------------------------------------------------------------------------- #
# The suite needs two directories on sys.path: the one holding the compiled
# `_ga_py` extension, and `ga_py/python/` (the pure-Python package). Discover
# them relative to this file instead of requiring PYTHONPATH, so that a plain
# `pytest ga_py/tests/` works in a fresh checkout with no environment setup.
#
# The binary directory is conventionally `build/` beside the source root. When
# this repository is configured as part of a larger enclosing build (added via
# add_subdirectory), CMake mirrors the source layout underneath that build's
# binary directory --- the artifacts then sit at `<ancestor>/build/<this root's
# path relative to that ancestor>/ga_py`. Walking a few ancestors and deriving
# that relative path covers both shapes without hardcoding either.
#
# Multi-config generators (Visual Studio, Xcode) additionally nest the artifact
# in a per-configuration subdirectory; single-config ones (Ninja, Makefiles)
# leave it in place. Both are probed, optimized configurations first.

_ROOT = Path(__file__).resolve().parents[2]  # <root>/ga_py/tests -> <root>
_ENCLOSING_LEVELS = 4  # how far up to look for an enclosing build tree

# Match the module itself, not build-system scaffolding that shares its stem --- an
# MSBuild tree keeps `_ga_py.dir/` and `_ga_py.vcxproj` next to the configuration
# subdirectories, and a bare `_ga_py.*` glob happily matches those.
_EXT_SUFFIXES = tuple(importlib.machinery.EXTENSION_SUFFIXES)  # .pyd / .so variants


def _holds_extension(d: Path) -> bool:
    return d.is_dir() and any(
        f.is_file() and f.name.endswith(_EXT_SUFFIXES) for f in d.glob("_ga_py*")
    )


def _extension_dirs() -> list[Path]:
    """Directories that actually contain a built `_ga_py`, in preference order."""
    found = []
    for ancestor in (_ROOT, *list(_ROOT.parents)[:_ENCLOSING_LEVELS]):
        base = ancestor / "build" / _ROOT.relative_to(ancestor) / "ga_py"
        for d in (base, base / "Release", base / "RelWithDebInfo", base / "Debug"):
            if _holds_extension(d):
                found.append(d)
    return found


# Prepended in reverse so the preference order above survives into sys.path.
for _p in reversed([_ROOT / "ga_py" / "python", *_extension_dirs()]):
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


# --------------------------------------------------------------------------- #
# Venv sanity check
# --------------------------------------------------------------------------- #

_MISSING_DEP_HINT = (
    "\n"
    "ga_py tests require the wrapper venv (ga_py/.venv), not the bindgen venv.\n"
    "Activate it and re-run:\n"
    "    source ga_py/.venv/bin/activate\n"
    "    PYTHONPATH='build/ga_py:ga_py/python' pytest ga_py/tests/\n"
    "If ga_py/.venv does not exist, create it once:\n"
    "    python3 -m venv ga_py/.venv\n"
    "    ga_py/.venv/bin/pip install pytest hypothesis numpy\n"
    "See ga_py/README.md > 'Two virtual environments'.\n"
)

try:
    import hypothesis  # noqa: F401  (import-only check)
except ImportError as exc:  # pragma: no cover --- environment guard
    raise ImportError(_MISSING_DEP_HINT) from exc

try:
    import ga_py  # noqa: F401  (import-only check)
except ImportError as exc:  # pragma: no cover --- environment guard
    raise ImportError(
        "ga_py is not importable. The bootstrap in conftest.py did not find a built\n"
        "extension --- build it first:\n"
        "    cmake --build build --target _ga_py\n"
        "(configure with -D_GA_BUILD_PYTHON=ON if you have not already; on a\n"
        "multi-config generator add --config Release).\n"
        f"searched: {[str(d) for d in _extension_dirs()] or 'nothing found'}\n"
        f"sys.path = {sys.path}\n"
    ) from exc


# --------------------------------------------------------------------------- #
# Tolerance constants
# --------------------------------------------------------------------------- #
# Mirrors the C++ side's is_congruent relative-tolerance approach.
# Default to 1e-12 for elementary ops; loosen for ops that accumulate FP
# error (wedge-then-divide, motor compositions, etc.).

EPS_TIGHT = 1e-12   # elementary ops on basis values
EPS_DEFAULT = 1e-9  # most algebraic identities
EPS_LOOSE = 1e-7    # accumulator-heavy ops (motor sandwich, normalize chains)


@pytest.fixture
def eps_tight() -> float:
    return EPS_TIGHT


@pytest.fixture
def eps() -> float:
    return EPS_DEFAULT


@pytest.fixture
def eps_loose() -> float:
    return EPS_LOOSE


# --------------------------------------------------------------------------- #
# Numerical comparison helpers
# --------------------------------------------------------------------------- #

def approx_eq(a: float, b: float, tol: float = EPS_DEFAULT) -> bool:
    """Relative-tolerance scalar equality: |a-b| <= tol * max(|a|, |b|, 1).

    Mirrors the relative-tolerance pattern documented in CLAUDE.md.
    Falls back to absolute tolerance for values near zero (the max(..., 1)).
    """
    return abs(a - b) <= tol * max(abs(a), abs(b), 1.0)


# --------------------------------------------------------------------------- #
# Component extraction --- canonical field order for every bound user type
# --------------------------------------------------------------------------- #
# Used by test_constants.py and test_cross_check.py to compare bound objects
# against expected component tuples. Order MUST match the C++ side and the
# ga_test/python_utilities/export_python_cases.cpp emitter.

FIELD_ORDER = {
    "vec2d": ("x", "y"),
    "vec3d": ("x", "y", "z"),
    "vec2dp": ("x", "y", "z"),
    "vec3dp": ("x", "y", "z", "w"),
    "bivec3d": ("x", "y", "z"),
    "bivec2dp": ("x", "y", "z"),
    "bivec3dp": ("vx", "vy", "vz", "mx", "my", "mz"),
    "trivec3dp": ("x", "y", "z", "w"),
    "vec4ds": ("x", "y", "z", "w"),
    "bivec4ds": ("vx", "vy", "vz", "mx", "my", "mz"),
    "trivec4ds": ("x", "y", "z", "w"),
    "dualnum2dp": ("c0", "c1"),
    "dualnum3dp": ("c0", "c1"),
    # cga: the field names are positional reuses of the shared templates --
    # Vec4_t for vec2dc/trivec2dc, BVec6_t for bivec2dc, Vec5_t for
    # vec3dc/quadvec3dc, BVec10_t for bivec3dc/trivec3dc (see
    # ga/detail/type_t/ga_type{2,3}dc.hpp for what each slot means).
    "vec2dc": ("x", "y", "z", "w"),
    "bivec2dc": ("vx", "vy", "vz", "mx", "my", "mz"),
    "trivec2dc": ("x", "y", "z", "w"),
    "vec3dc": ("x", "y", "z", "w", "u"),
    "bivec3dc": ("vx", "vy", "vz", "mx", "my", "mz", "px", "py", "pz", "pw"),
    "trivec3dc": ("vx", "vy", "vz", "mx", "my", "mz", "px", "py", "pz", "pw"),
    "quadvec3dc": ("x", "y", "z", "w", "u"),
}


def components(obj) -> tuple[float, ...]:
    """Return the numeric components of obj as a tuple, in canonical order."""
    cls = type(obj).__name__
    if cls.startswith("mvec"):
        n = sum(1 for a in dir(obj) if a.startswith("c") and a[1:].isdigit())
        return tuple(getattr(obj, f"c{i}") for i in range(n))
    if cls.startswith("pscalar") or cls.startswith("scalar"):
        return (float(obj),)
    if cls in ("float", "int"):
        return (float(obj),)
    return tuple(getattr(obj, f) for f in FIELD_ORDER[cls])
