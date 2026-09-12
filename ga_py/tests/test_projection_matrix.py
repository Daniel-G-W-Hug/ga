"""The projection / rejection overloads completed by the 2026-09 audit.

The four pairs added there -- ega3d bivector onto bivector, sta4ds vector and
bivector onto a hyperplane (trivector), pga3dp line rejected from a plane --
are checked here through the Python bindings, on the same properties the C++
suite pins (`ga_test/src/ga_{ega3d,sta4ds,pga3dp}_test.hpp`, the
"projection contract" cases) but with randomized inputs:

  * idempotence -- a blade lying in the target comes back unchanged. This is
    the property that catches a reversion sign, which containment and both
    scaling contracts do not see.
  * the target's scale must not reach the result.
  * a projection is an outermorphism: P(a ^ b) == P(a) ^ P(b). This ties the
    grade-1 and grade-2 overloads together; note it is blind to a sign error
    in the grade-1 one, where the flip cancels.
  * the rejection stays a BLADE -- in 4d a bivector is a 2-plane / line only
    if wdg(B, B) == 0. It holds because every target here is a hyperplane; for
    a bivector target in 4d the same subtraction leaves a non-simple bivector,
    which is why no such overload exists.
"""

from __future__ import annotations

import pytest
from hypothesis import HealthCheck, assume, given, settings, strategies as st

import ga_py
from conftest import EPS_LOOSE

ega = ga_py.ega
pga = ga_py.pga
sta = ga_py.sta

_settings = settings(max_examples=200, deadline=None,
                     suppress_health_check=[HealthCheck.filter_too_much])

_f = st.floats(min_value=-20.0, max_value=20.0, allow_nan=False,
               allow_infinity=False)


def _close(a, b, nrm_sq, scale=1.0):
    """relative comparison on GA objects, via the algebra's own squared norm"""
    return abs(nrm_sq(a - b)) <= EPS_LOOSE * max(abs(nrm_sq(a)), abs(nrm_sq(b)),
                                                 scale)


@st.composite
def _bivec3d(draw):
    return ega.bivec3d(draw(_f), draw(_f), draw(_f))


@st.composite
def _vec4ds(draw):
    return sta.vec4ds(draw(_f), draw(_f), draw(_f), draw(_f))


@st.composite
def _trivec4ds(draw):
    return sta.trivec4ds(draw(_f), draw(_f), draw(_f), draw(_f))


# --------------------------------------------------------------------------- #
# ega3d: bivector onto bivector
# --------------------------------------------------------------------------- #

@_settings
@given(_bivec3d(), _bivec3d())
def test_ega3d_bivec_projection(A, B):
    assume(ega.nrm_sq(B) > 1.0)
    par = ega.project_onto(A, B)
    rej = ega.reject_from(A, B)
    # the split is complete and orthogonal
    assert _close(par + rej, A, ega.nrm_sq, ega.nrm_sq(A))
    # the bound needs a floor: hypothesis happily draws a subnormal A (1e-217),
    # where nrm(A) underflows to zero and a purely relative bound is unsatisfiable
    assert abs(float(ega.dot(rej, B))) <= EPS_LOOSE * max(ega.nrm(A) * ega.nrm(B),
                                                          1.0)
    # idempotence -- the sign gate
    assert _close(ega.project_onto(par, B), par, ega.nrm_sq, ega.nrm_sq(A))
    assert _close(ega.project_onto(B, B), B, ega.nrm_sq, ega.nrm_sq(B))
    # the target's scale stays out
    assert _close(ega.project_onto(A, 3.0 * B), par, ega.nrm_sq, ega.nrm_sq(A))


# --------------------------------------------------------------------------- #
# sta4ds: onto a hyperplane, both source grades
# --------------------------------------------------------------------------- #

@_settings
@given(_vec4ds(), _vec4ds(), _trivec4ds())
def test_sta4ds_hyperplane_projection(a, b, t):
    assume(abs(sta.nrm_sq(t)) > 1.0)
    pa = sta.project_onto(a, t)
    ra = sta.reject_from(a, t)
    assert _close(pa + ra, a, sta.nrm_sq, abs(sta.nrm_sq(a)))
    # idempotence, and the target's scale stays out
    assert _close(sta.project_onto(pa, t), pa, sta.nrm_sq, abs(sta.nrm_sq(a)))
    assert _close(sta.project_onto(a, 3.0 * t), pa, sta.nrm_sq,
                  abs(sta.nrm_sq(a)))
    # the grades agree: a projection distributes over the wedge
    B = sta.wdg(a, b)
    lhs = sta.wdg(pa, sta.project_onto(b, t))
    rhs = sta.project_onto(B, t)
    assert _close(lhs, rhs, sta.nrm_sq, abs(sta.nrm_sq(lhs)) + 1.0)
    # the rejection of the bivector is still a 2-plane (the target is a
    # hyperplane, so the leftover is simple)
    rB = sta.reject_from(B, t)
    scale = abs(sta.nrm_sq(rB)) + 1.0
    assert abs(float(sta.wdg(rB, rB))) <= EPS_LOOSE * scale


# --------------------------------------------------------------------------- #
# pga3dp: a line rejected from a plane
# --------------------------------------------------------------------------- #

def test_pga3dp_line_rejected_from_plane():
    """the rejection is the perpendicular through the piercing point"""
    T = pga.wdg(pga.wdg(pga.vec3dp(-1, -1, 0, 1), pga.vec3dp(1, -1, 1, 1)),
                pga.vec3dp(-1, 1, 1, 1))
    L = pga.wdg(pga.vec3dp(-1, 0.5, -1.5, 1), pga.vec3dp(1, 0.5, 2.5, 1))
    par = pga.project_onto(L, T)
    rej = pga.reject_from(L, T)
    assert par + rej == L
    # a line again: simple, and with a weight
    assert abs(float(pga.wdg(rej, rej))) < 1e-12
    assert pga.weight_nrm(rej) > 0.1
    # through the point where L pierces T, and perpendicular to T
    X = pga.unitize(pga.rwdg(L, T))
    # the bound is_congruent takes its tolerance explicitly (no default through the
    # binding): eps_congruent, the same 4.5e3 machine epsilons the C++ default uses
    tol = 4.5e3 * 2.220446049250313e-16
    assert pga.is_congruent(rej, pga.r_weight_expand3dp(X, T), tol)
    assert pga.is_congruent(pga.rwdg(rej, T), pga.rwdg(L, T), tol)
    # the target's scale stays out, the source's does not
    assert pga.is_close(pga.reject_from(L, 3.0 * T), rej, tol)
    assert pga.is_close(pga.reject_from(2.0 * L, T), 2.0 * rej, tol)
