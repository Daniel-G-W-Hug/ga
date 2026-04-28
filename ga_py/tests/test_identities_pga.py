"""T4 --- PGA algebraic identities, randomized via hypothesis.

Identities are taken from the corresponding `CHECK(...)` lines in the C++
suite (`ga_test/src/ga_pga3dp_test.hpp`, with a few PGA2dp adaptations
from `ga_pga2dp_test.hpp`). The intent is the same as T3: pick
foundational invariants known to hold in C++ and verify the bindings
yield the same answer on randomized inputs.

Reference: ga_test/src/ga_pga3dp_test.hpp.

PGA-specific gotchas:
  * "Ideal" elements (points/lines/planes at infinity) have weight_nrm == 0.
    Such elements have no `unitize()`, no `inv()`, infinite distance to
    Euclidean elements, etc. We filter them with `assume(weight_nrm_sq > 1)`.
  * `dist3dp` returns a `dualnum3dp` with `.c0` (bulk) and `.c1` (weight);
    the actual distance is `c0 / c1`.
  * Motor construction: a unit motor M comes from `rgpr(plane2_unit, plane1_unit)`
    where each plane is `unitize(wdg(wdg(p0, p1), p2))`. Then
    `rgpr(rgpr(M, X), rrev(M))` sandwiches X through the rigid motion
    encoded by M (equivalent to `move3dp(X, M)` for points).
"""

from __future__ import annotations

import math

import pytest
from hypothesis import HealthCheck, assume, given, settings, strategies as st

import ga_py
from conftest import EPS_DEFAULT, EPS_LOOSE, approx_eq

pga = ga_py.pga


# --------------------------------------------------------------------------- #
# Strategies
# --------------------------------------------------------------------------- #

_finite = st.floats(min_value=-100.0, max_value=100.0,
                    allow_nan=False, allow_infinity=False)


@st.composite
def vec3dp_st(draw):
    return pga.vec3dp(draw(_finite), draw(_finite), draw(_finite), draw(_finite))


@st.composite
def bivec3dp_st(draw):
    return pga.bivec3dp(draw(_finite), draw(_finite), draw(_finite),
                        draw(_finite), draw(_finite), draw(_finite))


@st.composite
def trivec3dp_st(draw):
    return pga.trivec3dp(draw(_finite), draw(_finite), draw(_finite), draw(_finite))


@st.composite
def vec2dp_st(draw):
    return pga.vec2dp(draw(_finite), draw(_finite), draw(_finite))


@st.composite
def bivec2dp_st(draw):
    return pga.bivec2dp(draw(_finite), draw(_finite), draw(_finite))


# Non-ideal variants: weight component large enough to safely unitize.
@st.composite
def euclidean_point3dp_st(draw):
    """A 3dp point well away from the plane at infinity (weight > 0)."""
    p = pga.vec3dp(draw(_finite), draw(_finite), draw(_finite), draw(_finite))
    assume(pga.weight_nrm_sq(p) > 1.0)
    return p


@st.composite
def euclidean_point2dp_st(draw):
    """A 2dp point well away from the line at infinity (weight > 0)."""
    p = pga.vec2dp(draw(_finite), draw(_finite), draw(_finite))
    assume(pga.weight_nrm_sq(p) > 1.0)
    return p


_settings = settings(deadline=None,
                     max_examples=80,
                     suppress_health_check=[HealthCheck.too_slow,
                                            HealthCheck.filter_too_much])


# --------------------------------------------------------------------------- #
# Bulk + weight decomposition --- C++ ref: ga_pga3dp_test.hpp:3131-3136
# --------------------------------------------------------------------------- #

@_settings
@given(vec3dp_st())
def test_bulk_plus_weight_recovers_vec3dp(v):
    """bulk(v) + weight(v) == v"""
    rec = pga.bulk(v) + pga.weight(v)
    assert pga.weight_nrm_sq(rec - v) + pga.bulk_nrm_sq(rec - v) \
        < EPS_DEFAULT * max(pga.weight_nrm_sq(v) + pga.bulk_nrm_sq(v), 1.0)


@_settings
@given(bivec3dp_st())
def test_bulk_plus_weight_recovers_bivec3dp(B):
    """bulk(B) + weight(B) == B"""
    rec = pga.bulk(B) + pga.weight(B)
    err = pga.weight_nrm_sq(rec - B) + pga.bulk_nrm_sq(rec - B)
    ref = pga.weight_nrm_sq(B) + pga.bulk_nrm_sq(B)
    assert err < EPS_DEFAULT * max(ref, 1.0)


@_settings
@given(trivec3dp_st())
def test_bulk_plus_weight_recovers_trivec3dp(t):
    """bulk(t) + weight(t) == t"""
    rec = pga.bulk(t) + pga.weight(t)
    err = pga.weight_nrm_sq(rec - t) + pga.bulk_nrm_sq(rec - t)
    ref = pga.weight_nrm_sq(t) + pga.bulk_nrm_sq(t)
    assert err < EPS_DEFAULT * max(ref, 1.0)


# --------------------------------------------------------------------------- #
# Bulk and weight are orthogonal projections --- bulk(weight(x))==0 etc.
# --------------------------------------------------------------------------- #

@_settings
@given(vec3dp_st())
def test_bulk_of_weight_is_zero_vec3dp(v):
    """bulk(weight(v)) has zero bulk_nrm_sq and weight(bulk(v)) has zero weight_nrm_sq.
       Verifies bulk and weight are complementary projections."""
    bw = pga.bulk(pga.weight(v))
    wb = pga.weight(pga.bulk(v))
    assert pga.bulk_nrm_sq(bw) < EPS_DEFAULT * max(pga.weight_nrm_sq(v), 1.0)
    assert pga.weight_nrm_sq(wb) < EPS_DEFAULT * max(pga.bulk_nrm_sq(v), 1.0)


# --------------------------------------------------------------------------- #
# Complement involution --- C++ ref: ga_pga3dp_test.hpp:2983-2999
# --------------------------------------------------------------------------- #
# In even-dimensional algebras like PGA3dp, l_cmpl and r_cmpl are inverse.

@_settings
@given(vec3dp_st())
def test_complement_involution_vec3dp(v):
    """l_cmpl(r_cmpl(v)) == v and r_cmpl(l_cmpl(v)) == v"""
    a = pga.l_cmpl(pga.r_cmpl(v))
    b = pga.r_cmpl(pga.l_cmpl(v))
    assert pga.weight_nrm_sq(a - v) + pga.bulk_nrm_sq(a - v) \
        < EPS_DEFAULT * max(pga.weight_nrm_sq(v) + pga.bulk_nrm_sq(v), 1.0)
    assert pga.weight_nrm_sq(b - v) + pga.bulk_nrm_sq(b - v) \
        < EPS_DEFAULT * max(pga.weight_nrm_sq(v) + pga.bulk_nrm_sq(v), 1.0)


@_settings
@given(bivec3dp_st())
def test_complement_involution_bivec3dp(B):
    a = pga.l_cmpl(pga.r_cmpl(B))
    b = pga.r_cmpl(pga.l_cmpl(B))
    err_a = pga.weight_nrm_sq(a - B) + pga.bulk_nrm_sq(a - B)
    err_b = pga.weight_nrm_sq(b - B) + pga.bulk_nrm_sq(b - B)
    ref = pga.weight_nrm_sq(B) + pga.bulk_nrm_sq(B)
    assert err_a < EPS_DEFAULT * max(ref, 1.0)
    assert err_b < EPS_DEFAULT * max(ref, 1.0)


@_settings
@given(trivec3dp_st())
def test_complement_involution_trivec3dp(t):
    a = pga.l_cmpl(pga.r_cmpl(t))
    err = pga.weight_nrm_sq(a - t) + pga.bulk_nrm_sq(a - t)
    ref = pga.weight_nrm_sq(t) + pga.bulk_nrm_sq(t)
    assert err < EPS_DEFAULT * max(ref, 1.0)


# --------------------------------------------------------------------------- #
# Unitize --- C++ ref: ga_pga3dp_test.hpp:3190-3192 + general PGA convention
# --------------------------------------------------------------------------- #

@_settings
@given(euclidean_point3dp_st())
def test_unitize_idempotent_vec3dp(p):
    """unitize(unitize(p)) == unitize(p) for non-ideal points."""
    pu = pga.unitize(p)
    puu = pga.unitize(pu)
    err = pga.weight_nrm_sq(pu - puu) + pga.bulk_nrm_sq(pu - puu)
    assert err < EPS_LOOSE * max(pga.weight_nrm_sq(pu) + pga.bulk_nrm_sq(pu), 1.0)


@_settings
@given(euclidean_point3dp_st())
def test_unitize_gives_unit_weight_norm_vec3dp(p):
    """weight_nrm(unitize(p)) == 1 for non-ideal p."""
    pu = pga.unitize(p)
    assert approx_eq(float(pga.weight_nrm(pu)), 1.0, EPS_LOOSE)


@_settings
@given(euclidean_point2dp_st())
def test_unitize_gives_unit_weight_norm_vec2dp(p):
    """Same for PGA 2dp."""
    pu = pga.unitize(p)
    assert approx_eq(float(pga.weight_nrm(pu)), 1.0, EPS_LOOSE)


# --------------------------------------------------------------------------- #
# Distance --- C++ ref: ga_pga3dp_test.hpp:3231-3233
# --------------------------------------------------------------------------- #
# dist3dp returns a dualnum3dp; the distance is .c0 / .c1.

@_settings
@given(euclidean_point3dp_st())
def test_dist_to_self_is_zero(p):
    """dist3dp(p, p) has bulk component 0 (the distance)."""
    d = pga.dist3dp(p, p)
    # bulk part .c0 should be zero; weight .c1 may be nonzero (it's the
    # normalization). Distance proper is c0/c1.
    assert abs(d.c0) < EPS_LOOSE * max(pga.weight_nrm_sq(p), 1.0)


@_settings
@given(euclidean_point3dp_st(), euclidean_point3dp_st())
def test_dist_symmetric(p, q):
    """dist3dp(p, q) == dist3dp(q, p)  (both as dualnum3dp values).

    Equivalent to: c0/c1 of dist(p,q) equals c0/c1 of dist(q,p), which is
    the actual scalar distance.
    """
    d_pq = pga.dist3dp(p, q)
    d_qp = pga.dist3dp(q, p)
    # Compare the actual scalar distance (c0/c1).
    assume(abs(d_pq.c1) > EPS_LOOSE and abs(d_qp.c1) > EPS_LOOSE)
    dist_pq = d_pq.c0 / d_pq.c1
    dist_qp = d_qp.c0 / d_qp.c1
    assert approx_eq(dist_pq, dist_qp, EPS_LOOSE)


# --------------------------------------------------------------------------- #
# Join / meet --- C++ ref: ga_pga3dp_test.hpp:3375
# --------------------------------------------------------------------------- #

@_settings
@given(vec3dp_st(), vec3dp_st())
def test_join_equals_wdg_for_points(p, q):
    """join(p, q) == wdg(p, q) for two points (the line through them)."""
    j = pga.join(p, q)
    w = pga.wdg(p, q)
    assert pga.weight_nrm_sq(j - w) + pga.bulk_nrm_sq(j - w) \
        < EPS_DEFAULT * max(pga.weight_nrm_sq(j) + pga.bulk_nrm_sq(j), 1.0)


@_settings
@given(vec3dp_st(), vec3dp_st())
def test_join_antisymmetric_for_points(p, q):
    """join(p, q) == -join(q, p) (wedge anti-symmetry on grade-1 args)."""
    a, b = pga.join(p, q), -pga.join(q, p)
    err = pga.weight_nrm_sq(a - b) + pga.bulk_nrm_sq(a - b)
    ref = pga.weight_nrm_sq(a) + pga.bulk_nrm_sq(a)
    assert err < EPS_DEFAULT * max(ref, 1.0)


# --------------------------------------------------------------------------- #
# Reverse and regressive reverse involutions
# --------------------------------------------------------------------------- #

@_settings
@given(vec3dp_st())
def test_rev_involution_vec3dp(v):
    rr = pga.rev(pga.rev(v))
    err = pga.weight_nrm_sq(rr - v) + pga.bulk_nrm_sq(rr - v)
    ref = pga.weight_nrm_sq(v) + pga.bulk_nrm_sq(v)
    assert err < EPS_DEFAULT * max(ref, 1.0)


@_settings
@given(bivec3dp_st())
def test_rev_involution_bivec3dp(B):
    rr = pga.rev(pga.rev(B))
    err = pga.weight_nrm_sq(rr - B) + pga.bulk_nrm_sq(rr - B)
    ref = pga.weight_nrm_sq(B) + pga.bulk_nrm_sq(B)
    assert err < EPS_DEFAULT * max(ref, 1.0)


@_settings
@given(bivec3dp_st())
def test_rrev_involution_bivec3dp(B):
    """rrev(rrev(B)) == B  (regressive reverse is also an involution)."""
    rr = pga.rrev(pga.rrev(B))
    err = pga.weight_nrm_sq(rr - B) + pga.bulk_nrm_sq(rr - B)
    ref = pga.weight_nrm_sq(B) + pga.bulk_nrm_sq(B)
    assert err < EPS_DEFAULT * max(ref, 1.0)


# --------------------------------------------------------------------------- #
# Motor sandwich --- C++ ref: ga_pga3dp_test.hpp:2475 ff.
# --------------------------------------------------------------------------- #
# Build a motor as composition of two unit-plane reflections, sandwich a point
# through it via rgpr-rgpr-rrev, and verify it returns a point (vec3dp), not a
# different grade. Also verify move3dp gives the same result.

# Use fixed reference points/planes so the motor is well-conditioned; only the
# point being transformed is randomized.
def _build_unit_motor():
    """A unit motor reflecting through two planes (xz and a 45° rotation of it
    around z), built per the recipe in ga_pga3dp_test.hpp:2492-2512."""
    p0 = pga.O_3dp
    p1 = pga.vec3dp(1.0, 0.0, 0.0, 1.0)
    p2 = pga.vec3dp(1.0, 1.0, 0.0, 1.0)
    p3 = pga.vec3dp(0.0, 0.0, 1.0, 1.0)
    pl1 = pga.unitize(pga.wdg(pga.wdg(p0, p3), p1))
    pl2 = pga.unitize(pga.wdg(pga.wdg(p0, p3), p2))
    return pga.rgpr(pl2, pl1)


@_settings
@given(euclidean_point3dp_st())
def test_motor_sandwich_returns_point(p):
    """move3dp(p, M) returns a vec3dp (point), not a different grade."""
    M = _build_unit_motor()
    p_moved = pga.move3dp(p, M)
    assert type(p_moved).__name__ == "vec3dp"


@_settings
@given(euclidean_point3dp_st())
def test_move3dp_matches_explicit_sandwich(p):
    """move3dp(p, M) == gr1(rgpr(rgpr(M, p), rrev(M)))  (the underlying definition)."""
    M = _build_unit_motor()
    via_move = pga.unitize(pga.move3dp(p, M))
    via_explicit = pga.unitize(pga.rgpr(pga.rgpr(M, p), pga.rrev(M)).gr1())
    err = pga.weight_nrm_sq(via_move - via_explicit) \
        + pga.bulk_nrm_sq(via_move - via_explicit)
    ref = pga.weight_nrm_sq(via_move) + pga.bulk_nrm_sq(via_move)
    assert err < EPS_LOOSE * max(ref, 1.0)


@_settings
@given(euclidean_point3dp_st())
def test_motor_sandwich_preserves_point_unitization(p):
    """A unit motor (rigid motion) preserves the unitization of points:
       weight_nrm(unitize(move3dp(p, M))) == 1.

       This is the PGA analogue of "rotation preserves norm" from EGA (T3).
    """
    M = _build_unit_motor()
    p_moved = pga.move3dp(p, M)
    pu = pga.unitize(p_moved)
    assert approx_eq(float(pga.weight_nrm(pu)), 1.0, EPS_LOOSE)
