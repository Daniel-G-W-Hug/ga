"""Conformal geometric algebra through the Python bindings.

The C++ suites pin the algebra itself (ga_cga_test); what this file checks is
that the SAME questions can be asked from Python and give the same answers --
the objects, the properties read off them, and the conformal transformations.
The numbers here are the ones the C++ tests and ga_docu use, so a divergence
between the two surfaces shows up as a failure rather than as two truths.
"""

from __future__ import annotations

import math

import pytest

from ga_py import cga

# the bindings expose no default tolerance -- eps_congruent (C++ default) spelled out
TOL = 4.5e3 * 2.220446049250313e-16


# --------------------------------------------------------------------------- #
# objects: joins, meets, incidence
# --------------------------------------------------------------------------- #


def test_null_embedding_and_distance():
    """A radius-zero round point squares to zero, and the dot product of two of
    them is -1/2 the squared Euclidean distance."""
    a = cga.round_point3dc(0.0, 0.0, 0.0, 0.0)
    b = cga.round_point3dc(3.0, 4.0, 0.0, 0.0)
    assert float((a * a).gr0()) == pytest.approx(0.0, abs=1e-12)
    assert float(cga.dot(a, b)) == pytest.approx(-12.5)  # -0.5 * 5^2


def test_round_objects_are_joins_of_their_points():
    """dipole = 2 points, circle = 3, sphere = 4 (cga3dc)."""
    a = cga.round_point3dc(1.0, 0.0, 0.0, 0.0)
    b = cga.round_point3dc(0.0, 1.0, 0.0, 0.0)
    c = cga.round_point3dc(0.0, 0.0, 1.0, 0.0)
    d = cga.round_point3dc(-1.0, 0.0, 0.0, 0.0)

    unit_sphere = cga.sphere3dc(0.0, 0.0, 0.0, 1.0)
    assert cga.is_congruent(cga.wdg(cga.wdg(cga.wdg(a, b), c), d), unit_sphere, TOL)

    unit_circle = cga.circle3dc(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0)
    assert cga.is_congruent(cga.wdg(cga.wdg(a, b), d), unit_circle, TOL)


def test_flats_are_rounds_through_infinity():
    """A flat is the same join with the point at infinity appended."""
    a = cga.round_point3dc(1.0, 0.0, 0.0, 0.0)
    assert cga.is_congruent(cga.wdg(a, cga.e5_3dc), cga.flat_point3dc(1.0, 0.0, 0.0), TOL)


def test_incidence_is_the_vanishing_wedge():
    s = cga.sphere3dc(0.0, 0.0, 0.0, 1.0)
    on = cga.round_point3dc(0.0, 1.0, 0.0, 0.0)
    off = cga.round_point3dc(0.0, 2.0, 0.0, 0.0)
    assert float(cga.wdg(on, s)) == pytest.approx(0.0, abs=1e-12)
    assert float(cga.wdg(off, s)) != pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# center, radius and the properties
# --------------------------------------------------------------------------- #


def test_radius_is_signed_and_center_is_a_round_point():
    s = cga.unitize(cga.sphere3dc(3.0, 4.0, 0.0, 2.0))
    assert cga.radius_sq(s) == pytest.approx(4.0)
    assert cga.center_nrm(s) == pytest.approx(5.0)  # |(3,4,0)|
    assert float(cga.dot(s, s)) == pytest.approx(-4.0)  # dot square = -r^2
    # the center carries position AND radius, so it is the same round point
    assert cga.is_congruent(cga.cen(s), cga.round_point3dc(3.0, 4.0, 0.0, 2.0), TOL)


def test_properties_car_con_par():
    a = cga.round_point3dc(1.0, 2.0, 3.0, 2.0)
    s = cga.sphere3dc(1.0, 2.0, 3.0, 2.0)
    c = cga.circle3dc(1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 1.0)
    # container of a round point is the sphere of the same center and radius
    assert cga.is_congruent(cga.con(a), s, TOL)
    # carrier of a circle is its plane; carrier and container factor the object
    assert cga.is_congruent(cga.rwdg(cga.car(c), cga.con(c)), c, TOL)
    # the partner negates r^2 and is orthogonal to the object
    assert float(cga.dot(s, cga.par(s))) == pytest.approx(0.0, abs=1e-12)
    assert cga.radius_sq(cga.unitize(cga.par(s))) == pytest.approx(-4.0)


# --------------------------------------------------------------------------- #
# intersections: one chain, no case distinction (ga_docu 1.4.4)
# --------------------------------------------------------------------------- #


def _chord(sphere, line):
    """meet -> radius sign -> unitize -> cen +/- sqrt(r^2) att."""
    d = cga.rwdg(sphere, line)
    r_sq = cga.radius_sq(d)
    if r_sq < 0.0:
        return None
    du = cga.unitize(d)
    ctr = cga.unitize(cga.cen(du))
    dir_ = cga.att(du)
    r = math.sqrt(cga.radius_sq(du))
    return (
        (ctr.x + r * dir_.x, ctr.y + r * dir_.y, ctr.z + r * dir_.z),
        (ctr.x - r * dir_.x, ctr.y - r * dir_.y, ctr.z - r * dir_.z),
    )


@pytest.mark.parametrize(
    "y,r_sq",
    [(1.0, 3.0), (2.0, 0.0), (3.0, -5.0)],  # secant (off-center), tangent, missing
)
def test_line_sphere_intersection_off_center(y, r_sq):
    s = cga.sphere3dc(0.0, 0.0, 0.0, 2.0)
    line = cga.line3dc(0.0, y, 0.0, 1.0, 0.0, 0.0)
    assert cga.radius_sq(cga.rwdg(s, line)) == pytest.approx(r_sq)

    pts = _chord(s, line)
    if r_sq < 0.0:
        assert pts is None
    else:
        x = math.sqrt(4.0 - y * y)
        got = sorted(p[0] for p in pts)
        assert got == pytest.approx([-x, x])
        assert all(p[1] == pytest.approx(y) for p in pts)


def test_plane_cut_gives_a_smaller_circle():
    """An off-center plane cut yields a circle of squared radius r^2 - d^2,
    centered on the cutting plane."""
    s = cga.sphere3dc(0.0, 0.0, 0.0, 2.0)
    c = cga.rwdg(s, cga.plane3dc(0.0, 0.0, 1.0, 1.0))
    assert cga.radius_sq(c) == pytest.approx(3.0)
    assert cga.is_congruent(
        cga.cen(cga.unitize(c)), cga.round_point3dc(0.0, 0.0, 1.0, math.sqrt(3.0)), TOL
    )


def test_two_circles_in_2dc():
    """The 2d showcase: intersecting, tangent and disjoint circles."""
    c1 = cga.circle2dc(0.0, 0.0, 2.0)
    assert cga.radius_sq(cga.rwdg(c1, cga.circle2dc(3.0, 0.0, 2.0))) == pytest.approx(
        1.75
    )
    assert cga.radius_sq(
        cga.rwdg(cga.circle2dc(0.0, 0.0, 1.0), cga.circle2dc(2.0, 0.0, 1.0))
    ) == pytest.approx(0.0, abs=1e-12)
    assert cga.radius_sq(
        cga.rwdg(cga.circle2dc(0.0, 0.0, 1.0), cga.circle2dc(5.0, 0.0, 1.0))
    ) == pytest.approx(-5.25)


# --------------------------------------------------------------------------- #
# conformal transformations
# --------------------------------------------------------------------------- #


def test_translation_moves_the_center_and_keeps_the_radius():
    s = cga.sphere3dc(0.0, 0.0, 0.0, 2.0)
    moved = cga.unitize(cga.transform(s, cga.get_translation(1.0, 2.0, 3.0)))
    assert cga.radius_sq(moved) == pytest.approx(4.0)
    assert cga.is_congruent(cga.cen(moved), cga.round_point3dc(1.0, 2.0, 3.0, 2.0), TOL)


def test_dilation_scales_the_radius_by_sigma_squared():
    s = cga.sphere3dc(0.0, 0.0, 0.0, 2.0)
    scaled = cga.unitize(cga.transform(s, cga.get_dilation(0.0, 0.0, 0.0, 3.0)))
    assert cga.radius_sq(scaled) == pytest.approx(4.0 * 9.0)


def test_versors_double_cover_the_transformations():
    """M and -M act identically, which is what is_same_transform decides."""
    M = cga.get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.7)
    assert cga.is_same_transform(M, -M, TOL)
    assert not cga.is_same_transform(
        M, cga.get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.9), TOL
    )


def test_rexp_rlog_round_trip():
    M = cga.get_rotation(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.5)
    assert cga.is_same_transform(cga.rexp(cga.rlog(M)), M, TOL)
    # the square root applied twice returns the motion
    h = cga.rsqrt(M)
    assert cga.is_same_transform(cga.rgpr(h, h), M, TOL)


def test_inversion_in_a_sphere_is_its_own_inverse():
    s = cga.sphere3dc(0.0, 0.0, 0.0, 1.0)
    p = cga.round_point3dc(2.0, 0.0, 0.0, 0.0)
    once = cga.unitize(cga.invert_on(p, s))
    assert once.x == pytest.approx(0.5)  # r^2 / x
    twice = cga.unitize(cga.invert_on(cga.invert_on(p, s), s))
    assert twice.x == pytest.approx(2.0)
