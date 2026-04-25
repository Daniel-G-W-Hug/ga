"""Tests for multivector grade-composition constructors and grade extraction.

Validates the round trip: build an mvec from grade components, then extract
each grade — the extracted parts should equal the originals.
"""
from __future__ import annotations

import pytest

import ga_py


# --- Grade-composition constructors ----------------------------------------

def test_mvec3d_full_composition_from_four_grades():
    s = ga_py.ega.scalar3d(1.0)
    v = ga_py.ega.vec3d(2.0, 3.0, 4.0)
    B = ga_py.ega.bivec3d(5.0, 6.0, 7.0)
    p = ga_py.ega.pscalar3d(8.0)
    m = ga_py.ega.mvec3d(s, v, B, p)
    # Components stored in order: c0=s, c1..c3=v.x..z, c4..c6=B.x..z, c7=ps
    assert m.c0 == 1.0
    assert (m.c1, m.c2, m.c3) == (2.0, 3.0, 4.0)
    assert (m.c4, m.c5, m.c6) == (5.0, 6.0, 7.0)
    assert m.c7 == 8.0


def test_mvec3d_from_single_grade_zeros_others():
    """mvec3d(v) should populate the vector slots only; other grades = 0."""
    v = ga_py.ega.vec3d(2.0, 3.0, 4.0)
    m = ga_py.ega.mvec3d(v)
    assert m.c0 == 0.0
    assert (m.c1, m.c2, m.c3) == (2.0, 3.0, 4.0)
    assert (m.c4, m.c5, m.c6) == (0.0, 0.0, 0.0)
    assert m.c7 == 0.0


def test_mvec3d_e_from_scalar_and_bivector():
    s = ga_py.ega.scalar3d(1.0)
    B = ga_py.ega.bivec3d(2.0, 3.0, 4.0)
    me = ga_py.ega.mvec3d_e(s, B)
    assert me.c0 == 1.0
    assert (me.c1, me.c2, me.c3) == (2.0, 3.0, 4.0)


def test_mvec3d_from_even_subalgebra():
    """mvec3d(MVec3d_E) should embed even-grade components in the right slots."""
    me = ga_py.ega.mvec3d_e(ga_py.ega.scalar3d(1.0),
                            ga_py.ega.bivec3d(2.0, 3.0, 4.0))
    m = ga_py.ega.mvec3d(me)
    assert m.c0 == 1.0                     # scalar part
    assert (m.c1, m.c2, m.c3) == (0, 0, 0) # vector slots empty
    assert (m.c4, m.c5, m.c6) == (2.0, 3.0, 4.0)
    assert m.c7 == 0.0                     # pseudoscalar empty


def test_dualnum_from_scalar_and_pscalar():
    """DualNum is grade 0 + highest grade. dualnum2dp has scalar + pscalar parts."""
    dn = ga_py.pga.dualnum2dp(ga_py.pga.scalar2dp(7.0),
                              ga_py.pga.pscalar2dp(11.0))
    assert dn.c0 == 7.0
    assert dn.c1 == 11.0


def test_mvec2dp_full_composition():
    s = ga_py.pga.scalar2dp(1.0)
    v = ga_py.pga.vec2dp(2.0, 3.0, 4.0)
    B = ga_py.pga.bivec2dp(5.0, 6.0, 7.0)
    p = ga_py.pga.pscalar2dp(8.0)
    m = ga_py.pga.mvec2dp(s, v, B, p)
    assert m.c0 == 1.0
    assert (m.c1, m.c2, m.c3) == (2.0, 3.0, 4.0)
    assert (m.c4, m.c5, m.c6) == (5.0, 6.0, 7.0)
    assert m.c7 == 8.0


# --- Grade extraction round trip -------------------------------------------

def test_mvec3d_grade_extraction_round_trip():
    """Build mvec from grades; extract each grade; equality checks."""
    s_in = ga_py.ega.scalar3d(1.0)
    v_in = ga_py.ega.vec3d(2.0, 3.0, 4.0)
    B_in = ga_py.ega.bivec3d(5.0, 6.0, 7.0)
    p_in = ga_py.ega.pscalar3d(8.0)
    m = ga_py.ega.mvec3d(s_in, v_in, B_in, p_in)

    assert isinstance(m.gr0(), ga_py.ega.scalar3d)
    assert isinstance(m.gr1(), ga_py.ega.vec3d)
    assert isinstance(m.gr2(), ga_py.ega.bivec3d)
    assert isinstance(m.gr3(), ga_py.ega.pscalar3d)

    assert m.gr0() == s_in
    assert m.gr1() == v_in
    assert m.gr2() == B_in
    assert m.gr3() == p_in


def test_mvec3d_e_only_has_even_grades():
    me = ga_py.ega.mvec3d_e()
    # gr0 and gr2 must exist; gr1 and gr3 must NOT exist for even multivectors.
    assert hasattr(me, "gr0")
    assert hasattr(me, "gr2")
    assert not hasattr(me, "gr1")
    assert not hasattr(me, "gr3")


def test_mvec3d_u_only_has_odd_grades():
    mu = ga_py.ega.mvec3d_u()
    assert not hasattr(mu, "gr0")
    assert hasattr(mu, "gr1")
    assert not hasattr(mu, "gr2")
    assert hasattr(mu, "gr3")


def test_dualnum2dp_grades_are_0_and_3():
    """In PGA 2dp the dualnum carries scalar (grade 0) + pscalar (grade 3)."""
    dn = ga_py.pga.dualnum2dp(ga_py.pga.scalar2dp(2.0),
                              ga_py.pga.pscalar2dp(5.0))
    assert hasattr(dn, "gr0")
    assert hasattr(dn, "gr3")
    assert not hasattr(dn, "gr1")
    assert not hasattr(dn, "gr2")
    assert isinstance(dn.gr0(), ga_py.pga.scalar2dp)
    assert isinstance(dn.gr3(), ga_py.pga.pscalar2dp)
    assert float(dn.gr0()) == 2.0
    assert float(dn.gr3()) == 5.0


def test_grade_extraction_returns_correct_python_type():
    """Cross-algebra check: gr2 of an mvec2d returns pscalar2d (grade 2 in 2D
    is the pseudoscalar), not bivec — there is no bivec in 2D."""
    m = ga_py.ega.mvec2d(1.0, 2.0, 3.0, 4.0)
    g2 = m.gr2()
    assert isinstance(g2, ga_py.ega.pscalar2d)


def test_compose_then_decompose_then_recompose():
    """Build → extract grades → rebuild → equality."""
    m1 = ga_py.ega.mvec3d(ga_py.ega.scalar3d(1.0),
                          ga_py.ega.vec3d(2.0, 3.0, 4.0),
                          ga_py.ega.bivec3d(5.0, 6.0, 7.0),
                          ga_py.ega.pscalar3d(8.0))
    m2 = ga_py.ega.mvec3d(m1.gr0(), m1.gr1(), m1.gr2(), m1.gr3())
    assert m1 == m2
