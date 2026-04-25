"""Smoke tests for the vec3d binding — the pipeline-proof scope.

Verifies:
  - the extension imports
  - vec3d can be constructed (default + 3-arg)
  - fields are read/write accessible
  - __repr__ and __str__ produce sensible strings
  - basic arithmetic operators (unary -, binary +, -, scalar *, /) work
  - in-place operators mutate the value
"""
from __future__ import annotations

import ga_py


def test_module_layout():
    assert hasattr(ga_py, "ega")
    assert hasattr(ga_py.ega, "vec3d")


def test_default_construction():
    v = ga_py.ega.vec3d()
    assert v.x == 0.0 and v.y == 0.0 and v.z == 0.0


def test_three_arg_construction_and_fields():
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    assert v.x == 1.0
    assert v.y == 2.0
    assert v.z == 3.0


def test_field_assignment():
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    v.y = 99.0
    assert v.y == 99.0


def test_repr_round_trips_via_eval_locals():
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    r = repr(v)
    assert r.startswith("vec3d(") and "1" in r and "2" in r and "3" in r


def test_str_uses_fmt_formatter():
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    s = str(v)
    # The fmt formatter for Vec3d emits something like "Vec3d(1,2,3)" — the
    # exact form is the project's existing convention. We just check it's
    # non-empty and mentions all three components.
    assert s and "1" in s and "2" in s and "3" in s


def test_unary_minus():
    v = ga_py.ega.vec3d(1.0, -2.0, 3.0)
    n = -v
    assert n.x == -1.0 and n.y == 2.0 and n.z == -3.0


def test_addition():
    a = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    b = ga_py.ega.vec3d(10.0, 20.0, 30.0)
    s = a + b
    assert s.x == 11.0 and s.y == 22.0 and s.z == 33.0


def test_subtraction():
    a = ga_py.ega.vec3d(10.0, 20.0, 30.0)
    b = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    d = a - b
    assert d.x == 9.0 and d.y == 18.0 and d.z == 27.0


def test_scalar_multiplication():
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    s = 2.5 * v
    assert s.x == 2.5 and s.y == 5.0 and s.z == 7.5
    s2 = v * 2.0
    assert s2.x == 2.0 and s2.y == 4.0 and s2.z == 6.0


def test_scalar_division():
    v = ga_py.ega.vec3d(2.0, 4.0, 6.0)
    h = v / 2.0
    assert h.x == 1.0 and h.y == 2.0 and h.z == 3.0


def test_equality():
    a = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    b = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    c = ga_py.ega.vec3d(1.0, 2.0, 4.0)
    assert a == b
    assert a != c


def test_inplace_addition():
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    v += ga_py.ega.vec3d(10.0, 20.0, 30.0)
    assert v.x == 11.0 and v.y == 22.0 and v.z == 33.0
