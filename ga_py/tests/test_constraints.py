"""Tests for the closed-loop constraint bindings and the scoped-enum support
that unlocks them.

The C++ closed-loop layer (ga/ga_pga2dp_ops_constraints.hpp) is built around a
pure-data descriptor, loop_constraint2dp, whose `type` field is the scoped enum
constraint2dp. Binding a struct with an enum field required teaching ga_bindgen
to bind scoped enums (nb::enum_); that same support also exposes the joint
enums joint2dp / joint3dp and the joint_state2dp / joint_state3dp descriptors.

The stateful closed_loop_system2dp class itself is NOT bound (like
dynamic_system2dp it is reconstructed in Python from primitives when needed);
these tests cover the bound value descriptors and enums.
"""
from __future__ import annotations

import pytest

import ga_py
from ga_py import pga


# --- scoped enums --------------------------------------------------------

def test_constraint2dp_enum_members():
    # the loop-joint kinds: pin, rod, weld (D1)
    assert list(pga.constraint2dp.__members__) == ["coincidence", "distance", "frame"]


def test_constraint3dp_enum_members():
    assert list(pga.constraint3dp.__members__) == ["coincidence", "distance", "frame"]


def test_joint_enum_members():
    # the 1-dof screw joints plus, in 3D, the motor joints (F): the enum ORDER is part
    # of the surface (it is the bound value)
    assert list(pga.joint2dp.__members__) == ["free", "revolute", "prismatic"]
    assert list(pga.joint3dp.__members__) == [
        "free", "revolute", "prismatic", "helical", "cylindrical", "spherical", "planar",
    ]


def test_enum_values_distinct_and_equal_to_self():
    assert pga.joint2dp.revolute == pga.joint2dp.revolute
    assert pga.joint2dp.revolute != pga.joint2dp.prismatic


# --- loop_constraint2dp (the closed-loop descriptor) ---------------------

def _make_coincidence(fa=0, fb=1):
    # the bound ctor is all-fields positional: (..., type, length, active)
    return pga.loop_constraint2dp(
        fa, pga.vec2dp(1.0, 0.0, 1.0),
        fb, pga.vec2dp(2.0, 0.0, 1.0),
        pga.constraint2dp.coincidence, 0.0, True,
    )


def test_loop_constraint2dp_kinds_and_state_round_trip():
    # the three constraint kinds, the rod length and the active flag (D1)
    rod = pga.loop_constraint2dp(0, pga.vec2dp(1.0, 0.0, 1.0),
                                 1, pga.vec2dp(2.0, 0.0, 1.0),
                                 pga.constraint2dp.distance, 0.75, False)
    assert rod.type == pga.constraint2dp.distance
    assert rod.length == 0.75
    assert rod.active is False
    rod.active = True
    rod.type = pga.constraint2dp.frame
    assert rod.active is True and rod.type == pga.constraint2dp.frame
    assert "frame" in str(rod) and "active" in str(rod)


def test_loop_constraint2dp_field_access():
    c = _make_coincidence()
    assert c.frame_a == 0
    assert c.frame_b == 1
    assert c.anchor_a == pga.vec2dp(1.0, 0.0, 1.0)
    assert c.anchor_b == pga.vec2dp(2.0, 0.0, 1.0)
    assert c.type == pga.constraint2dp.coincidence


def test_loop_constraint2dp_default_ctor():
    c = pga.loop_constraint2dp()
    # default-initialised: zeroed frames/anchors, coincidence type (the in-class
    # default member initialiser)
    assert c.frame_a == 0
    assert c.frame_b == 0
    assert c.type == pga.constraint2dp.coincidence


def test_loop_constraint2dp_fields_mutable():
    c = _make_coincidence()
    c.frame_a = 5
    c.anchor_a = pga.vec2dp(3.0, 4.0, 1.0)
    assert c.frame_a == 5
    assert c.anchor_a == pga.vec2dp(3.0, 4.0, 1.0)


def test_loop_constraint2dp_equality():
    assert _make_coincidence() == _make_coincidence()
    assert _make_coincidence(fa=0) != _make_coincidence(fa=9)


def test_loop_constraint2dp_repr_and_str():
    c = _make_coincidence()
    for s in (repr(c), str(c)):
        assert "loop_constraint2dp" in s
        assert "coincidence" in s
        assert "frame_a = 0" in s


def test_loop_constraint2dp_format_spec():
    c = _make_coincidence()
    # a precision spec flows into the vec2dp anchors' scalars
    out = format(c, ".2f")
    assert "1.00" in out and "2.00" in out
    # an empty spec equals the default str()
    assert format(c, "") == str(c)


def test_loop_constraint2dp_bad_format_spec_raises():
    c = _make_coincidence()
    with pytest.raises(ValueError):
        format(c, "not-a-real-spec")


# --- loop_constraint3dp (the 3D closed-loop descriptor) ------------------

def _make_coincidence3d(fa=0, fb=1):
    return pga.loop_constraint3dp(
        fa, pga.vec3dp(1.0, 0.0, 0.0, 1.0),
        fb, pga.vec3dp(2.0, 0.0, 0.0, 1.0),
        pga.constraint3dp.coincidence, 0.0, True,
    )


def test_loop_constraint3dp_kinds_and_state_round_trip():
    rod = pga.loop_constraint3dp(0, pga.vec3dp(1.0, 0.0, 0.0, 1.0),
                                 1, pga.vec3dp(2.0, 0.0, 0.0, 1.0),
                                 pga.constraint3dp.distance, 0.75, False)
    assert rod.type == pga.constraint3dp.distance
    assert rod.length == 0.75
    assert rod.active is False
    assert "distance" in str(rod)


def test_loop_constraint3dp_field_access():
    c = _make_coincidence3d()
    assert c.frame_a == 0
    assert c.frame_b == 1
    assert c.anchor_a == pga.vec3dp(1.0, 0.0, 0.0, 1.0)
    assert c.anchor_b == pga.vec3dp(2.0, 0.0, 0.0, 1.0)
    assert c.type == pga.constraint3dp.coincidence


def test_loop_constraint3dp_equality_and_repr():
    assert _make_coincidence3d() == _make_coincidence3d()
    assert _make_coincidence3d(fa=0) != _make_coincidence3d(fa=7)
    s = repr(_make_coincidence3d())
    assert "loop_constraint3dp" in s and "coincidence" in s


def test_loop_constraint3dp_format_spec():
    out = format(_make_coincidence3d(), ".2f")
    assert "1.00" in out and "2.00" in out


# --- joint_state2dp / joint_state3dp (unlocked by the joint enums) -------

def test_joint_state2dp_roundtrips():
    # trailing args = the spring/damper force elements (stiffness, damping, q_rest)
    j = pga.joint_state2dp(
        pga.joint2dp.revolute, pga.vec2dp(0.0, 0.0, 1.0),
        pga.mvec2dp_u(), 0.5, -1.0, 8.0, 0.3, 0.1, [], [], pga.mvec2dp_u(),
    )
    assert j.type == pga.joint2dp.revolute
    assert j.screw_b == pga.vec2dp(0.0, 0.0, 1.0)
    assert j.phi == pytest.approx(0.5)
    assert j.omega == pytest.approx(-1.0)
    assert j.stiffness == pytest.approx(8.0)
    assert j.damping == pytest.approx(0.3)
    assert j.q_rest == pytest.approx(0.1)
    assert "joint_state2dp" in repr(j)
    assert "revolute" in repr(j)


def test_joint_state3dp_roundtrips():
    # trailing args = the spring/damper force elements (stiffness, damping, q_rest)
    j = pga.joint_state3dp(
        pga.joint3dp.prismatic, pga.bivec3dp(0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
        pga.mvec3dp_e(), 2.0, 3.0, 8.0, 0.3, 0.1, [], [], pga.mvec3dp_e(),
    )
    assert j.type == pga.joint3dp.prismatic
    assert j.phi == pytest.approx(2.0)
    assert j.omega == pytest.approx(3.0)
    assert j.stiffness == pytest.approx(8.0)
    assert j.damping == pytest.approx(0.3)
    assert j.q_rest == pytest.approx(0.1)
    assert "joint_state3dp" in repr(j)
    assert "prismatic" in repr(j)


# --- module surface ------------------------------------------------------

def test_constraint_symbols_exported_under_pga():
    for name in ("constraint2dp", "constraint3dp", "joint2dp", "joint3dp",
                 "loop_constraint2dp", "loop_constraint3dp",
                 "joint_state2dp", "joint_state3dp"):
        assert hasattr(ga_py.pga, name), f"ga_py.pga.{name} missing"
