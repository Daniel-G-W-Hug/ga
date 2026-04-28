"""T3 --- EGA algebraic identities, randomized via hypothesis.

The identities chosen below all appear as `CHECK(...)` lines in the C++
suite (`ga_test/src/ga_ega3d_test.hpp` and ga_ega2d_test.hpp): inner
product symmetry/linearity, wedge anti-symmetry, the geometric-product
split `v1*v2 = dot(v1,v2) + wdg(v1,v2)`, inverse, reverse involution,
and rotor unitarity. The C++ tests verify each on a handful of fixed
inputs; here we randomize the inputs via hypothesis so the same
algebraic invariants get hit on hundreds of cases per run.

Reference: ga_test/src/ga_ega3d_test.hpp.

Tolerance policy: relative tolerance via conftest.approx_eq, with
a slightly looser default for chains that involve wedge + inverse +
division (where rounding accumulates). Bounds on inputs are kept to
[-100, 100] so squared norms stay comfortably in double precision.
"""

from __future__ import annotations

import math

import pytest
from hypothesis import HealthCheck, assume, given, settings, strategies as st

import ga_py
from conftest import EPS_DEFAULT, EPS_LOOSE, approx_eq

ega = ga_py.ega


# --------------------------------------------------------------------------- #
# Strategies
# --------------------------------------------------------------------------- #
# Bound floats to a moderate range so squared norms stay well within double
# precision and so degenerate-input filtering (`assume nrm_sq > eps`) only
# rarely shrinks too aggressively.

_finite = st.floats(min_value=-100.0, max_value=100.0,
                    allow_nan=False, allow_infinity=False)
_nonzero = st.floats(min_value=-100.0, max_value=100.0,
                     allow_nan=False, allow_infinity=False) \
            .filter(lambda x: abs(x) > 1e-3)


@st.composite
def vec2d_st(draw):
    return ega.vec2d(draw(_finite), draw(_finite))


@st.composite
def vec3d_st(draw):
    return ega.vec3d(draw(_finite), draw(_finite), draw(_finite))


@st.composite
def bivec3d_st(draw):
    return ega.bivec3d(draw(_finite), draw(_finite), draw(_finite))


@st.composite
def nonnull_vec3d_st(draw):
    """Vector with squared norm comfortably above zero --- safe for inv()."""
    v = ega.vec3d(draw(_finite), draw(_finite), draw(_finite))
    assume(ega.nrm_sq(v) > 1.0)  # guards inv() and the inverse identities
    return v


@st.composite
def nonnull_bivec3d_st(draw):
    b = ega.bivec3d(draw(_finite), draw(_finite), draw(_finite))
    assume(ega.nrm_sq(b) > 1.0)
    return b


# Hypothesis settings: silence the "function is too slow" health check for
# the chains that go through inv(); they are still well under 1 ms / case.
_settings = settings(deadline=None,
                     max_examples=80,
                     suppress_health_check=[HealthCheck.too_slow])


# --------------------------------------------------------------------------- #
# Inner product properties --- C++ ref: ga_ega3d_test.hpp:346
# --------------------------------------------------------------------------- #

@_settings
@given(vec3d_st(), vec3d_st())
def test_dot_symmetry_vec3d(u, v):
    """dot(u, v) == dot(v, u)"""
    assert approx_eq(float(ega.dot(u, v)), float(ega.dot(v, u)), EPS_DEFAULT)


@_settings
@given(_finite, vec3d_st(), vec3d_st())
def test_dot_linearity_vec3d(a, u, v):
    """dot(a*u, v) == a * dot(u, v)"""
    lhs = float(ega.dot(a * u, v))
    rhs = a * float(ega.dot(u, v))
    assert approx_eq(lhs, rhs, EPS_DEFAULT)


@_settings
@given(vec3d_st(), vec3d_st(), vec3d_st())
def test_dot_distributive_vec3d(u, v, w):
    """dot(u + v, w) == dot(u, w) + dot(v, w)"""
    lhs = float(ega.dot(u + v, w))
    rhs = float(ega.dot(u, w)) + float(ega.dot(v, w))
    assert approx_eq(lhs, rhs, EPS_DEFAULT)


# --------------------------------------------------------------------------- #
# Wedge product properties --- C++ ref: ga_ega3d_test.hpp:541, :1220
# --------------------------------------------------------------------------- #

@_settings
@given(vec3d_st(), vec3d_st())
def test_wdg_antisymmetric_vec3d(u, v):
    """wdg(u, v) == -wdg(v, u)"""
    a, b = ega.wdg(u, v), -ega.wdg(v, u)
    # bivec3d equality compares all components; tolerance via component diff
    diff = a - b
    assert ega.nrm_sq(diff) < EPS_DEFAULT * max(ega.nrm_sq(a), 1.0)


@_settings
@given(vec3d_st())
def test_wdg_self_is_zero_vec3d(v):
    """wdg(v, v) == 0  (a vector wedged with itself spans no plane)"""
    w = ega.wdg(v, v)
    assert ega.nrm_sq(w) < EPS_DEFAULT * max(ega.nrm_sq(v) ** 2, 1.0)


@_settings
@given(vec2d_st(), vec2d_st())
def test_wdg_antisymmetric_vec2d(u, v):
    """Same anti-symmetry in 2D --- result is a pseudoscalar instead of bivector."""
    lhs = float(ega.wdg(u, v))      # pscalar2d in 2D
    rhs = -float(ega.wdg(v, u))
    assert approx_eq(lhs, rhs, EPS_DEFAULT)


# --------------------------------------------------------------------------- #
# Geometric-product split --- C++ ref: ga_ega3d_test.hpp:1715 (vector case)
# --------------------------------------------------------------------------- #

@_settings
@given(vec3d_st(), vec3d_st())
def test_gpr_split_vec_vec(u, v):
    """For vectors: u*v = dot(u,v) + wdg(u,v).
       Since gpr(vec, vec) -> mvec3d_e, check the grade-0 and grade-2 parts."""
    prod = u * v
    s = float(prod.gr0())
    B = prod.gr2()
    # Grade-0 part equals the dot product
    assert approx_eq(s, float(ega.dot(u, v)), EPS_DEFAULT)
    # Grade-2 part equals the wedge product
    diff = B - ega.wdg(u, v)
    assert ega.nrm_sq(diff) < EPS_DEFAULT * max(ega.nrm_sq(B), 1.0)


# --------------------------------------------------------------------------- #
# Reverse involution --- rev(rev(x)) == x for any grade
# --------------------------------------------------------------------------- #

@_settings
@given(vec3d_st())
def test_rev_involution_vec3d(v):
    rr = ega.rev(ega.rev(v))
    assert ega.nrm_sq(rr - v) < EPS_DEFAULT * max(ega.nrm_sq(v), 1.0)


@_settings
@given(bivec3d_st())
def test_rev_involution_bivec3d(B):
    rr = ega.rev(ega.rev(B))
    assert ega.nrm_sq(rr - B) < EPS_DEFAULT * max(ega.nrm_sq(B), 1.0)


# --------------------------------------------------------------------------- #
# Inverse --- C++ ref: ga_ega3d_test.hpp:401-409
# --------------------------------------------------------------------------- #
# For non-null v: v * inv(v) is a pure scalar 1 (its grade-0 part is 1,
# grade-2 part is 0).

@_settings
@given(nonnull_vec3d_st())
def test_inv_vec3d_left(v):
    """gr0(v * inv(v)) == 1, gr2(v * inv(v)) == 0"""
    prod = v * ega.inv(v)
    assert approx_eq(float(prod.gr0()), 1.0, EPS_LOOSE)
    assert ega.nrm_sq(prod.gr2()) < EPS_LOOSE


@_settings
@given(nonnull_vec3d_st())
def test_inv_equals_rev_over_nrm_sq_vec3d(v):
    """inv(v) == rev(v) / nrm_sq(v)  (closed form for vectors in EGA3D)"""
    expected = ega.rev(v) / ega.nrm_sq(v)
    diff = ega.inv(v) - expected
    assert ega.nrm_sq(diff) < EPS_LOOSE * max(ega.nrm_sq(ega.inv(v)), 1.0)


@_settings
@given(nonnull_bivec3d_st())
def test_inv_bivec3d_left(B):
    """gr0(B * inv(B)) == 1, gr2(B * inv(B)) == 0"""
    prod = B * ega.inv(B)
    assert approx_eq(float(prod.gr0()), 1.0, EPS_LOOSE)
    assert ega.nrm_sq(prod.gr2()) < EPS_LOOSE


# --------------------------------------------------------------------------- #
# Rotor properties --- C++ ref: ga_ega3d_test.hpp:727 (rotations)
# --------------------------------------------------------------------------- #
# A unit rotor R = exp(-B*phi/2) for unit bivector B and angle phi satisfies
# R * rev(R) == 1 and rotate(v, R) preserves nrm_sq(v).

# Use a fixed unit bivector (e12) and a randomized angle to keep the input
# space simple --- rotor construction with random bivectors invites numerical
# corner cases that aren't really testing the identity.
_angle_st = st.floats(min_value=-math.pi, max_value=math.pi,
                      allow_nan=False, allow_infinity=False)


@_settings
@given(_angle_st)
def test_unit_rotor_R_times_rev_R_is_one(phi):
    """R = get_rotor(e12, phi) is a unit rotor: R * rev(R) == 1."""
    R = ega.get_rotor(ega.e12_3d, phi)
    prod = R * ega.rev(R)
    # gr0 should be 1, gr2 should be 0
    assert approx_eq(float(prod.gr0()), 1.0, EPS_LOOSE)
    assert ega.nrm_sq(prod.gr2()) < EPS_LOOSE


@_settings
@given(_angle_st, vec3d_st())
def test_rotation_preserves_nrm_sq(phi, v):
    """rotate(v, R) preserves |v|^2 when R is a unit rotor."""
    R = ega.get_rotor(ega.e12_3d, phi)
    v_rot = ega.rotate(v, R)
    assert approx_eq(ega.nrm_sq(v_rot), ega.nrm_sq(v),
                     EPS_LOOSE * max(ega.nrm_sq(v), 1.0))


@_settings
@given(_angle_st, vec3d_st())
def test_rotate_matches_rotate_opt(phi, v):
    """rotate() and rotate_opt() must agree (both are bound; user might pick either)."""
    R = ega.get_rotor(ega.e12_3d, phi)
    a, b = ega.rotate(v, R), ega.rotate_opt(v, R)
    assert ega.nrm_sq(a - b) < EPS_LOOSE * max(ega.nrm_sq(a), 1.0)
