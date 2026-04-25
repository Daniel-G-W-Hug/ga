"""Coverage smoke tests across the 27 Vec/MVec-shape user types.

These are light tests: confirm each type is exposed under the right submodule,
can be constructed with N doubles, prints cleanly, and supports the basic
operator set. Heavier numerical / identity tests come with the cross-check
harness described in TODO/considerations_python_wrapper.md.
"""
from __future__ import annotations

import pytest

import ga_py


EGA_TYPES = [
    ("vec2d", 2), ("vec3d", 3), ("vec4d", 4),
    ("bivec3d", 3), ("bivec4d", 6),
    ("trivec4d", 4),
    ("mvec2d", 4), ("mvec2d_e", 2),
    ("mvec3d", 8), ("mvec3d_e", 4), ("mvec3d_u", 4),
    ("mvec4d", 16), ("mvec4d_e", 8), ("mvec4d_u", 8),
]

PGA_TYPES = [
    ("vec2dp", 3), ("vec3dp", 4),
    ("bivec2dp", 3), ("bivec3dp", 6),
    ("trivec3dp", 4),
    ("mvec2dp", 8), ("mvec2dp_e", 4), ("mvec2dp_u", 4),
    ("mvec3dp", 16), ("mvec3dp_e", 8), ("mvec3dp_u", 8),
    ("dualnum2dp", 2), ("dualnum3dp", 2),
]


@pytest.mark.parametrize("name,arity", EGA_TYPES)
def test_ega_type_exposed(name, arity):
    cls = getattr(ga_py.ega, name)
    obj = cls(*[float(i) for i in range(arity)])
    s = repr(obj)
    assert s.startswith(f"{name}(")
    assert str(obj)  # non-empty


@pytest.mark.parametrize("name,arity", PGA_TYPES)
def test_pga_type_exposed(name, arity):
    cls = getattr(ga_py.pga, name)
    obj = cls(*[float(i) for i in range(arity)])
    s = repr(obj)
    assert s.startswith(f"{name}(")
    assert str(obj)


def test_ega_and_pga_are_disjoint():
    """A PGA-only type must not appear in the EGA submodule, and vice versa."""
    assert hasattr(ga_py.pga, "vec3dp")
    assert not hasattr(ga_py.ega, "vec3dp")
    assert hasattr(ga_py.ega, "vec3d")
    assert not hasattr(ga_py.pga, "vec3d")


def test_vec2d_arithmetic():
    a = ga_py.ega.vec2d(1.0, 2.0)
    b = ga_py.ega.vec2d(10.0, 20.0)
    s = a + b
    assert s.x == 11.0 and s.y == 22.0
    d = b - a
    assert d.x == 9.0 and d.y == 18.0
    assert (-a).x == -1.0 and (-a).y == -2.0
    assert (2.5 * a).x == 2.5


def test_bivec3d_arithmetic():
    a = ga_py.ega.bivec3d(1.0, 2.0, 3.0)
    b = ga_py.ega.bivec3d(0.5, 0.5, 0.5)
    s = a + b
    assert (s.x, s.y, s.z) == (1.5, 2.5, 3.5)


def test_mvec3d_e_arithmetic():
    """Even multivector in 3D has 4 components (c0..c3)."""
    a = ga_py.ega.mvec3d_e(1.0, 2.0, 3.0, 4.0)
    b = ga_py.ega.mvec3d_e(10.0, 20.0, 30.0, 40.0)
    s = a + b
    assert (s.c0, s.c1, s.c2, s.c3) == (11.0, 22.0, 33.0, 44.0)
    # scalar multiplication
    d = s / 11.0
    assert d.c0 == pytest.approx(1.0) and d.c3 == pytest.approx(4.0)


def test_dualnum2dp_arithmetic():
    """Dual number in PGA 2dp has 2 components."""
    a = ga_py.pga.dualnum2dp(1.0, 2.0)
    b = ga_py.pga.dualnum2dp(3.0, 4.0)
    s = a + b
    assert (s.c0, s.c1) == (4.0, 6.0)


def test_trivec3dp_has_four_components():
    t = ga_py.pga.trivec3dp(1.0, 2.0, 3.0, 4.0)
    assert (t.x, t.y, t.z, t.w) == (1.0, 2.0, 3.0, 4.0)


def test_mvec_fields_named_cN():
    """Multivectors expose components as c0..c{N-1}, not x/y/z."""
    m = ga_py.ega.mvec3d(1, 2, 3, 4, 5, 6, 7, 8)
    for i in range(8):
        assert getattr(m, f"c{i}") == float(i + 1)


def test_equality_is_type_aware():
    a = ga_py.ega.vec3d(1, 2, 3)
    b = ga_py.ega.vec3d(1, 2, 3)
    c = ga_py.ega.vec3d(1, 2, 4)
    assert a == b
    assert a != c
