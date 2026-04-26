"""T5 --- Cross-check the Python wrapper against C++ reference output.

Loads `data/ga_test_cases.json` (emitted by ga_test/src/export_python_cases.cpp)
and replays each case through the Python bindings, asserting the result
matches the captured C++ output within tolerance.

This is the most direct equivalence test: any drift between the C++
implementation and the Python wrapper --- whether from a binding bug,
a wrong overload selection, a sign flip, or a coercion mishap ---
surfaces here.

Workflow:
  1. Build the exporter:   cmake --build build --target ga_export_python_cases
  2. Run it once:          build/ga_test/ga_export_python_cases ga_py/tests/data/ga_test_cases.json
  3. Run pytest:           pytest ga_py/tests/test_cross_check.py
The JSON snapshot is committed to git, so step 2 is only needed when the
case set is updated (new ops, new types, new representative inputs).

Plan reference: TODO/considerations_python_wrapper.md sec. 7.3, step T5.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import ga_py
from conftest import components


# --------------------------------------------------------------------------- #
# Locate and load the JSON snapshot
# --------------------------------------------------------------------------- #

_CASES_PATH = Path(__file__).parent / "data" / "ga_test_cases.json"


def _load_cases() -> list[dict]:
    if not _CASES_PATH.exists():
        pytest.skip(
            f"{_CASES_PATH.relative_to(Path(__file__).parent.parent.parent)} "
            f"is missing. Generate it with:\n"
            f"  cmake --build build --target ga_export_python_cases\n"
            f"  build/ga_test/ga_export_python_cases "
            f"{_CASES_PATH}\n"
        )
    with _CASES_PATH.open() as f:
        doc = json.load(f)
    assert doc["version"] == 1, f"unsupported schema version: {doc['version']}"
    return doc["cases"]


_CASES = _load_cases()


# --------------------------------------------------------------------------- #
# Type reconstruction --- {"type": "...", "fields": [...]}  ->  Python object
# --------------------------------------------------------------------------- #

# Where each type lives in the Python module hierarchy.
_TYPE_LOCATIONS: dict[str, str] = {}
for n in dir(ga_py.ega):
    if isinstance(getattr(ga_py.ega, n), type):
        _TYPE_LOCATIONS[n] = "ega"
for n in dir(ga_py.pga):
    if isinstance(getattr(ga_py.pga, n), type):
        _TYPE_LOCATIONS[n] = "pga"


def _reconstruct(value: dict):
    """Reconstruct a {"type": ..., "fields": [...]} entry as a Python object."""
    type_name = value["type"]
    fields = value["fields"]
    if type_name == "double":
        # Plain numeric scalar (e.g. nrm_sq result).
        assert len(fields) == 1, f"'double' must have one field, got {fields}"
        return float(fields[0])
    submod_name = _TYPE_LOCATIONS.get(type_name)
    if submod_name is None:
        raise ValueError(f"Unknown bound type: {type_name!r}. "
                         f"Add to _TYPE_LOCATIONS or register in ga_bindgen.")
    cls = getattr(getattr(ga_py, submod_name), type_name)
    return cls(*fields)


# --------------------------------------------------------------------------- #
# Operation dispatch --- name (+ submodule)  ->  callable
# --------------------------------------------------------------------------- #

def _dispatch(submodule: str, op: str):
    """Resolve op name to a callable. Handles operator symbols and free functions."""
    if op == "+":
        return lambda a, b: a + b
    if op == "-":
        return lambda a, b: a - b
    if op == "*":
        return lambda a, b: a * b
    if op == "neg":
        return lambda a: -a
    submod = getattr(ga_py, submodule) if submodule != "top" else ga_py
    return getattr(submod, op)


# --------------------------------------------------------------------------- #
# Comparison
# --------------------------------------------------------------------------- #

def _close_enough(actual, expected, tolerance: float) -> tuple[bool, str]:
    """Component-wise comparison with relative tolerance.

    Returns (ok, diagnostic_message).
    """
    a_comps = (float(actual),) if isinstance(actual, (int, float)) else components(actual)
    e_comps = (float(expected),) if isinstance(expected, (int, float)) else components(expected)
    if len(a_comps) != len(e_comps):
        return False, (f"component count mismatch: actual {len(a_comps)} "
                       f"({type(actual).__name__}) vs expected {len(e_comps)} "
                       f"({type(expected).__name__})")
    for i, (a, e) in enumerate(zip(a_comps, e_comps)):
        # Relative tolerance --- mirrors the C++ is_congruent pattern.
        bound = tolerance * max(abs(a), abs(e), 1.0)
        if abs(a - e) > bound:
            return False, (f"component {i} differs: {a!r} vs {e!r} "
                           f"(diff {abs(a - e):.3g}, allowed {bound:.3g})")
    return True, ""


# --------------------------------------------------------------------------- #
# The actual cross-check test
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("case", _CASES, ids=lambda c: c["id"])
def test_cross_check(case: dict):
    """Replay one C++ case through the Python bindings and compare.

    Each case fixes (op, args, expected, tolerance). The Python result must
    match the captured C++ result within the per-case tolerance.
    """
    submodule = case["submodule"]
    op = case["op"]
    tolerance = case["tolerance"]

    # Reconstruct args.
    args = [_reconstruct(a) for a in case["args"]]

    # Dispatch op.
    fn = _dispatch(submodule, op)
    actual = fn(*args)

    # Reconstruct expected and verify it matches the actual return type.
    expected = _reconstruct(case["expected"])

    # Type identity check, with one allowed coercion: C++ "double" maps to
    # Python `float`, so an actual `int`/`float` should compare against an
    # expected `float`.
    expected_type = case["expected"]["type"]
    actual_type = (
        "double" if isinstance(actual, (int, float))
        else type(actual).__name__
    )
    assert actual_type == expected_type, (
        f"return type mismatch for case {case['id']!r}: "
        f"got {actual_type}, expected {expected_type}"
    )

    ok, msg = _close_enough(actual, expected, tolerance)
    assert ok, f"value mismatch for case {case['id']!r}: {msg}\n" \
               f"  actual:   {actual!r}\n" \
               f"  expected: {expected!r}"


# --------------------------------------------------------------------------- #
# Coverage sanity check --- not a value test, just a count
# --------------------------------------------------------------------------- #

def test_case_set_is_nonempty():
    """At least 50 cases should be present (the JSON exporter currently emits ~80).
       Catches an empty / corrupted JSON without exploding the parametrize loop."""
    assert len(_CASES) >= 50, (
        f"Only {len(_CASES)} cross-check cases loaded from "
        f"{_CASES_PATH.name}. Regenerate via "
        f"`cmake --build build --target ga_export_python_cases && "
        f"build/ga_test/ga_export_python_cases ga_py/tests/data/ga_test_cases.json`"
    )
