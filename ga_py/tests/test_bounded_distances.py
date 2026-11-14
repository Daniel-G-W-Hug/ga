"""Distances of bounded objects and the turn of a motor about a direction.

`segment_distance`, `box_faces` and `turn_about` are bound by the generator;
`convex_distance` is hand-bound (its faces are a std::span in C++), with the segment
form's out-parameters returned as a tuple (distance, t, face).
"""
from __future__ import annotations

import math

import pytest

from ga_py import pga


def p3(x, y, z):
    return pga.vec3dp(x, y, z, 1.0)


def p2(x, y):
    return pga.vec2dp(x, y, 1.0)


def test_segment_distance_3d():
    # skew, beyond an end (the unbounded lines would read 0.3), parallel
    assert pga.segment_distance(p3(-1, 0, 0), p3(1, 0, 0), p3(0, -1, 0.3),
                                p3(0, 1, 0.3)) == pytest.approx(0.3, abs=1e-14)
    assert pga.segment_distance(p3(-1, 0, 0), p3(1, 0, 0), p3(2, -1, 0.3),
                                p3(2, 1, 0.3)) == pytest.approx(math.hypot(1.0, 0.3))
    assert pga.segment_distance(p3(0, 0, 0), p3(1, 0, 0), p3(3, 0.5, 0),
                                p3(4, 0.5, 0)) == pytest.approx(math.hypot(2.0, 0.5))
    # the points' weights do not matter
    assert pga.segment_distance(pga.vec3dp(-2, 0, 0, 2), pga.vec3dp(-3, 0, 0, -3),
                                pga.vec3dp(0, -0.5, 0.15, 0.5),
                                pga.vec3dp(0, 4, 1.2, 4)) == pytest.approx(0.3)


def test_segment_distance_2d():
    assert pga.segment_distance(p2(-1, 0), p2(1, 0), p2(0, -1), p2(0, 1)) == \
        pytest.approx(0.0, abs=1e-15)
    assert pga.segment_distance(p2(-1, 0), p2(1, 0), p2(2, 0.3), p2(2, 1)) == \
        pytest.approx(math.hypot(1.0, 0.3))


def test_convex_distance_box_3d():
    faces = pga.box_faces(0.2, 0.3, 0.5, p3(1, 2, 3))
    assert len(faces) == 6
    assert pga.convex_distance(faces, p3(1.35, 2.0, 3.0)) == pytest.approx(0.25)
    assert pga.convex_distance(faces, p3(1.0, 2.0, 3.0)) == pytest.approx(-0.1)
    # the point form is the largest wedge with a face
    q = p3(1.13, 2.19, 3.0)
    assert pga.convex_distance(faces, q) == \
        pytest.approx(max(float(pga.wdg(q, f)) for f in faces))
    # the segment form: straight through the box along x, least at its middle
    d, t, face = pga.convex_distance(faces, p3(0.0, 2.0, 3.0), p3(2.0, 2.0, 3.0))
    assert d == pytest.approx(-0.1)
    assert t == pytest.approx(0.5)
    assert face in (0, 1)
    # past a corner: the least is where two faces' lines cross
    d, t, face = pga.convex_distance(faces, p3(1.3, 2.0, 3.0), p3(1.0, 2.4, 3.0))
    sampled = min(
        pga.convex_distance(faces, p3(1.3 - 0.3 * k / 1000, 2.0 + 0.4 * k / 1000, 3.0))
        for k in range(1001))
    assert d <= sampled + 1e-13
    assert sampled - d < 1e-3


def test_convex_distance_rectangle_2d():
    faces = pga.box_faces(0.2, 0.3, p2(1, 2))
    assert len(faces) == 4
    assert pga.convex_distance(faces, p2(1.35, 2.0)) == pytest.approx(0.25)
    assert pga.convex_distance(faces, p2(1.0, 2.0)) == pytest.approx(-0.1)
    d, t, face = pga.convex_distance(faces, p2(0.0, 2.0), p2(2.0, 2.0))
    assert d == pytest.approx(-0.1)
    assert t == pytest.approx(0.5)


def test_turn_about():
    z = pga.vec3dp(0.0, 0.0, 1.0, 0.0)
    x = pga.vec3dp(1.0, 0.0, 0.0, 0.0)
    origin = pga.vec3dp(0.0, 0.0, 0.0, 1.0)
    z_line = pga.unitize(pga.wdg(origin, p3(0, 0, 1)))
    x_line = pga.unitize(pga.wdg(origin, p3(1, 0, 0)))
    twist = pga.get_motor(z_line, 0.4)
    swing = pga.get_motor(x_line, 1.1)
    m = pga.rgpr(swing, twist)
    assert pga.turn_about(m, z) == pytest.approx(0.4)
    assert pga.turn_about(m, pga.vec3dp(0.0, 0.0, -2.0, 0.0)) == pytest.approx(-0.4)
    # a rotation about a perpendicular axis has no turn about z
    assert pga.turn_about(swing, z) == pytest.approx(0.0, abs=1e-15)
    assert pga.turn_about(swing, x) == pytest.approx(1.1)
