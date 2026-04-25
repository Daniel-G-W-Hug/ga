"""Tests for PGA geometric primitives — point*, vector*, line*, plane3d.

These types are C++ subclasses of vec/bivec/trivec bases (e.g.,
`Point2d : public Vec2d<T>`). In Python they are bound with nanobind
inheritance, so `isinstance(point2d_instance, vec2d) == True`.
"""
from __future__ import annotations

import pytest

import ga_py


# --- Module placement ------------------------------------------------------

PGA_PRIMITIVES = ["vector2d", "vector3d",
                  "point2d", "point3d", "point2dp", "point3dp",
                  "line2d", "line3d", "plane3d"]


@pytest.mark.parametrize("name", PGA_PRIMITIVES)
def test_lives_in_pga_submodule(name):
    """All PGA geometric primitives are routed to the pga submodule, even
    when the user-typedef name doesn't end in 'dp' (e.g. point2d, line2d).
    The submodule routing must use canonical_underlying, not just name suffix."""
    assert hasattr(ga_py.pga, name), f"{name} missing from ga_py.pga"
    assert not hasattr(ga_py.ega, name), f"{name} unexpectedly in ga_py.ega"


# --- Inheritance preserved -------------------------------------------------

def test_point2d_is_subclass_of_vec2d():
    p = ga_py.pga.point2d(3.0, 4.0)
    assert isinstance(p, ga_py.ega.vec2d)


def test_vector2d_is_subclass_of_vec2d():
    v = ga_py.pga.vector2d(1.0, 0.0)
    assert isinstance(v, ga_py.ega.vec2d)


def test_point2dp_is_subclass_of_vec2dp():
    p = ga_py.pga.point2dp(1.0, 2.0, 1.0)
    assert isinstance(p, ga_py.pga.vec2dp)


def test_line2d_is_subclass_of_bivec2dp():
    ell = ga_py.pga.line2d(1.0, 2.0, 3.0)
    assert isinstance(ell, ga_py.pga.bivec2dp)


def test_plane3d_is_subclass_of_trivec3dp():
    """Three coplanar points define a plane3d (TriVec3dp under the hood)."""
    pl = ga_py.pga.plane3d(ga_py.pga.point3d(1, 0, 0),
                           ga_py.pga.point3d(0, 1, 0),
                           ga_py.pga.point3d(0, 0, 1))
    assert isinstance(pl, ga_py.pga.trivec3dp)


# --- Inherited fields and arithmetic ---------------------------------------

def test_point2d_inherits_xy_fields():
    p = ga_py.pga.point2d(3.0, 4.0)
    assert p.x == 3.0
    assert p.y == 4.0


def test_inherited_field_assignment():
    p = ga_py.pga.point2d(1.0, 2.0)
    p.x = 99.0
    assert p.x == 99.0


def test_point2d_inherits_arithmetic_from_vec2d():
    """Adding two points slices to vec2d (this matches C++ behavior — no
    operator+ defined directly on point2d, so the inherited Vec2d operator
    runs and returns Vec2d)."""
    p = ga_py.pga.point2d(1.0, 2.0)
    q = ga_py.pga.point2d(3.0, 4.0)
    s = p + q
    # Result is the base vec2d, not point2d (C++ slicing).
    assert isinstance(s, ga_py.ega.vec2d)
    assert s.x == 4.0 and s.y == 6.0


# --- Convenience constructors ---------------------------------------------

def test_line2d_from_two_points():
    """Line2d(p, q) constructs the line through p and q.
    Per the C++ ctor in ga_type2dp.hpp:
        BiVec2dp(q.x - p.x, q.y - p.y, p.x * q.y - p.y * q.x)
    """
    p = ga_py.pga.point2d(0.0, 0.0)
    q = ga_py.pga.point2d(1.0, 1.0)
    ell = ga_py.pga.line2d(p, q)
    assert isinstance(ell, ga_py.pga.bivec2dp)
    # Expected from formula: (1-0, 1-0, 0*1 - 0*1) = (1, 1, 0)
    assert ell.x == 1.0 and ell.y == 1.0 and ell.z == 0.0


def test_line2d_from_point_and_direction():
    """Line2d(p, v) builds a line at point p with direction v."""
    p = ga_py.pga.point2d(2.0, 0.0)
    v = ga_py.ega.vec2d(0.0, 1.0)
    ell = ga_py.pga.line2d(p, v)
    assert isinstance(ell, ga_py.pga.bivec2dp)


def test_point2dp_from_point2d_lifts_to_homogeneous():
    """Point2dp(Point2d) sets the homogeneous z component to 1."""
    p2 = ga_py.pga.point2d(3.0, 4.0)
    p2dp = ga_py.pga.point2dp(p2)
    assert p2dp.x == 3.0
    assert p2dp.y == 4.0
    assert p2dp.z == 1.0


def test_repr_uses_derived_class_name():
    """`repr(point2d(...))` says `point2d(...)`, not `vec2d(...)` even
    though the fields come from the inherited base."""
    p = ga_py.pga.point2d(1.0, 2.0)
    r = repr(p)
    assert r.startswith("point2d(")


def test_str_via_inherited_fmt_formatter():
    p = ga_py.pga.point2d(1.0, 2.0)
    assert str(p)  # any non-empty string from the inherited fmt formatter
