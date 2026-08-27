"""body2dp / body3dp: the inertial record of a frame, and the make_*_body builders.

The two structs are pure data (inertia map, its cached inverse, the mass) whose
field type inertia{2,3}dp is HAND-bound; the generator resolves them through an
injected map entry (emit_nanobind.py). Before that entry every make_*_body
builder was bound but failed at the return value.
"""
from __future__ import annotations

import pytest

from ga_py import pga


def test_body3dp_default_and_fields():
    b = pga.body3dp()
    assert b.mass == 0.0
    assert isinstance(b.I, pga.inertia3dp)
    assert "body3dp" in repr(b) and "mass" in repr(b)


def test_make_cuboid_body_returns_a_body3dp():
    b = pga.make_cuboid_body(2.0, 0.3, 0.2, 0.5)
    assert isinstance(b, pga.body3dp)
    assert b.mass == pytest.approx(2.0)
    # the inertia map sends a pure translation twist to m times it (linear momentum)
    v = pga.bivec3dp(0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    p = b.I(v)
    comps = (p.vx, p.vy, p.vz, p.mx, p.my, p.mz)
    assert sum(c * c for c in comps) == pytest.approx(4.0)  # |m v| = 2
    # a massless carrier: the builder accepts m = 0 (the inverse is left zero)
    assert pga.make_cuboid_body(0.0, 0.05, 0.05, 0.05).mass == 0.0


def test_make_point_body3dp_and_2d_twin():
    b3 = pga.make_point_body3dp(1.5)
    assert isinstance(b3, pga.body3dp) and b3.mass == pytest.approx(1.5)
    b2 = pga.make_point_body2dp(1.5)
    assert isinstance(b2, pga.body2dp) and b2.mass == pytest.approx(1.5)


def test_body2dp_by_field_ctor_round_trips():
    plate = pga.make_plate_body(3.0, 0.8, 0.08)
    b = pga.body2dp(plate.I, pga.get_inertia_inverse(plate.I), 3.0)
    assert b.mass == pytest.approx(3.0)
    assert b.I[0, 0] == pytest.approx(plate.I[0, 0])
    assert "body2dp" in str(b)
