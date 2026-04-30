"""Smoke test for the type stubs in ga_py/python/ga_py/{__init__,ega,pga}.pyi.

The stubs are produced by `ga_bindgen/src/emit_stubs.py` (a thin wrapper
around `nanobind.stubgen` plus two post-passes — see that script's
docstring for what it does).

Goals here:
  * stubs exist alongside the package, plus the PEP 561 marker;
  * stubs parse as valid Python (catches mistakes in the post-passes);
  * key symbols a user would touch are present (smoke for stubgen
    coverage regressions);
  * post-pass invariants hold — `__format__` is on every bound class
    and the runtime forwarders from __init__.py appear in pga.pyi.

If the stubs drift out of sync with the bindings, the right fix is
`cmake --build build --target regenerate_python_stubs`, not editing
the stubs by hand.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest


# Resolve the stubs once. Doing it at import time keeps each test fast.
PKG_DIR = Path(__file__).resolve().parents[1] / "python" / "ga_py"
STUBS = {name: PKG_DIR / f"{name}.pyi"
         for name in ("__init__", "ega", "pga")}


# ---------------------------------------------------------------------------
# File presence + valid Python
# ---------------------------------------------------------------------------

def test_py_typed_marker_exists():
    """PEP 561: without this file mypy / pyright will not look at the stubs."""
    assert (PKG_DIR / "py.typed").is_file()


@pytest.mark.parametrize("name", list(STUBS))
def test_stub_file_exists(name: str):
    assert STUBS[name].is_file(), \
        f"missing {name}.pyi — run `cmake --build build " \
        "--target regenerate_python_stubs`"


@pytest.mark.parametrize("name", list(STUBS))
def test_stub_parses_as_python(name: str):
    """Catches syntax errors introduced by the post-passes (broken regex,
    indentation, etc.)."""
    text = STUBS[name].read_text()
    ast.parse(text)


# ---------------------------------------------------------------------------
# Coverage of the key public surface
# ---------------------------------------------------------------------------

def test_init_pyi_contains_top_level_symbols():
    text = STUBS["__init__"].read_text()
    for sym in ("def Hz2radps", "def deg2rad", "def rk4_get_time",
                "def rk4_step", "def sign"):
        assert sym in text, f"{sym} missing from __init__.pyi"


def test_ega_pyi_contains_core_types():
    text = STUBS["ega"].read_text()
    for cls in ("class scalar3d:", "class vec3d:", "class bivec3d:",
                "class pscalar3d:", "class mvec3d:", "class mvec3d_e:",
                "class mvec3d_u:"):
        assert cls in text


def test_pga_pyi_contains_core_types():
    text = STUBS["pga"].read_text()
    for cls in ("class vec2dp:", "class bivec2dp:", "class vec3dp:",
                "class bivec3dp:", "class trivec3dp:",
                "class mvec2dp:", "class mvec3dp:",
                # Hand-written physics surface.
                "class inertia2dp:", "class inertia3dp:"):
        assert cls in text


def test_pga_pyi_contains_physics_factories():
    text = STUBS["pga"].read_text()
    for fn in ("def get_point_inertia", "def get_plate_inertia",
               "def get_cuboid_inertia", "def get_inertia_inverse",
               "def compute_omega_dot"):
        assert fn in text


# ---------------------------------------------------------------------------
# Post-pass invariants
# ---------------------------------------------------------------------------

def _bound_classes(text: str) -> list[str]:
    """All `class NAME:` (top-level / no inheritance) and
    `class NAME(BASE):` declarations in a stub."""
    return re.findall(r"^class (\w+)(?:\([^)]*\))?:\s*$", text, re.MULTILINE)


@pytest.mark.parametrize("name", ("ega", "pga"))
def test_every_bound_class_has_format(name: str):
    """__format__ post-pass: nanobind.stubgen filters this dunder out;
    emit_stubs.py re-injects it into every class with a __str__ since
    the bindings always emit it (see emit_format_lambda in
    emit_nanobind.py and the hand-written bindings_*.cpp)."""
    text = STUBS[name].read_text()
    classes = _bound_classes(text)
    assert classes, f"no classes found in {name}.pyi"

    # Slice the file at each class boundary and verify __format__ is
    # present in every class block. Slicing keeps the assertion message
    # informative if a single class is missing the method.
    boundaries = [m.start() for m in
                  re.finditer(r"^class \w+", text, re.MULTILINE)]
    boundaries.append(len(text))
    blocks = {classes[i]: text[boundaries[i]:boundaries[i + 1]]
              for i in range(len(classes))}
    missing = [c for c, body in blocks.items()
               if "def __format__" not in body]
    assert not missing, \
        f"__format__ missing from {len(missing)} class(es): {missing}"


def test_pga_pyi_contains_runtime_forwarders():
    """Forwarder post-pass: pure-Python wrappers in __init__.py are
    appended as module-level decls in pga.pyi so static type checkers
    don't flag them as missing."""
    text = STUBS["pga"].read_text()
    expected = ("dist2dp", "dist3dp",
                "ortho_proj2dp", "ortho_proj3dp",
                "central_proj2dp", "central_proj3dp",
                "ortho_antiproj2dp", "ortho_antiproj3dp",
                "l_bulk_contract2dp", "r_bulk_contract2dp",
                "l_weight_contract3dp", "r_weight_expand3dp")
    for name in expected:
        # The forwarders are emitted as `def NAME(a: Any, b: Any) -> Any: ...`.
        assert re.search(rf"^def {name}\(", text, re.MULTILINE), \
            f"forwarder {name} missing from pga.pyi"


def test_pga_pyi_imports_Any_for_forwarders():
    """The forwarder block uses Any, so the stub must import it."""
    text = STUBS["pga"].read_text()
    assert re.search(r"^from typing import .*\bAny\b", text, re.MULTILINE)
