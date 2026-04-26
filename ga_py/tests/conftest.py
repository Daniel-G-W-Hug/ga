"""Shared test fixtures and helpers for the ga_py test suite.

This file is loaded automatically by pytest. It provides:

  * A clear error if the test session was launched from the wrong venv.
  * Per-algebra / per-operation tolerance fixtures.
  * An `approx_eq()` helper that mirrors the C++ `is_congruent` relative-
    tolerance pattern (see CLAUDE.md "Critical Numerical Precision Lessons").

Run from the wrapper venv:
    source ga_py/.venv/bin/activate
    PYTHONPATH="build/ga_py:ga_py/python" pytest ga_py/tests/

Do NOT use build/spike_libclang/.venv for tests --- that one is for
ga_bindgen. See ga_py/README.md > "Two virtual environments".
"""

from __future__ import annotations

import sys

import pytest

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
        "ga_py is not importable. Did you forget to set PYTHONPATH?\n"
        "    PYTHONPATH='build/ga_py:ga_py/python' pytest ga_py/tests/\n"
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
# ga_test/src/export_python_cases.cpp emitter.

FIELD_ORDER = {
    "vec2d": ("x", "y"),
    "vec3d": ("x", "y", "z"),
    "vec4d": ("x", "y", "z", "w"),
    "vec2dp": ("x", "y", "z"),
    "vec3dp": ("x", "y", "z", "w"),
    "bivec3d": ("x", "y", "z"),
    "bivec4d": ("vx", "vy", "vz", "mx", "my", "mz"),
    "bivec2dp": ("x", "y", "z"),
    "bivec3dp": ("vx", "vy", "vz", "mx", "my", "mz"),
    "trivec3dp": ("x", "y", "z", "w"),
    "trivec4d": ("x", "y", "z", "w"),
    "dualnum2dp": ("c0", "c1"),
    "dualnum3dp": ("c0", "c1"),
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
