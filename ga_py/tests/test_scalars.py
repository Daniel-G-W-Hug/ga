"""Tests for scalar and pseudoscalar types.

Key concern: because each tagged C++ specialization binds to a distinct
Python class, the C++ type-safety (no mixing scalar2d with scalar3d) is
preserved automatically — no explicit tag mechanism needed in Python.
"""
from __future__ import annotations

import pytest

import ga_py


ALL_SCALAR_TYPES = [
    ("ega", "scalar2d"), ("ega", "scalar3d"), ("ega", "scalar4d"),
    ("pga", "scalar2dp"), ("pga", "scalar3dp"),
    ("ega", "pscalar2d"), ("ega", "pscalar3d"), ("ega", "pscalar4d"),
    ("pga", "pscalar2dp"), ("pga", "pscalar3dp"),
]


@pytest.mark.parametrize("submod,name", ALL_SCALAR_TYPES)
def test_scalar_exposed(submod, name):
    cls = getattr(getattr(ga_py, submod), name)
    s = cls(3.5)
    assert float(s) == 3.5
    assert s.value == 3.5
    assert repr(s).startswith(f"{name}(")


@pytest.mark.parametrize("submod,name", ALL_SCALAR_TYPES)
def test_default_construction_is_zero(submod, name):
    cls = getattr(getattr(ga_py, submod), name)
    s = cls()
    assert float(s) == 0.0


def test_value_is_writable():
    s = ga_py.ega.scalar3d(1.0)
    s.value = 42.0
    assert float(s) == 42.0


def test_explicit_float_conversion():
    """Python __float__ mirrors C++ `explicit operator T&()`."""
    s = ga_py.ega.scalar3d(2.0)
    assert float(s) == 2.0
    # No implicit conversion: scalar is NOT a float in Python either.
    assert not isinstance(s, float)


def test_distinct_scalars_are_different_python_classes():
    """C++ uses type tags; Python preserves the distinction via class identity."""
    s2 = ga_py.ega.scalar2d(1.0)
    s3 = ga_py.ega.scalar3d(1.0)
    assert type(s2) is not type(s3)
    assert not isinstance(s2, type(s3))


def test_mixing_distinct_scalars_raises_type_error():
    """scalar2d + scalar3d must fail — they represent different algebras."""
    s2 = ga_py.ega.scalar2d(1.0)
    s3 = ga_py.ega.scalar3d(1.0)
    with pytest.raises(TypeError):
        s2 + s3


def test_scalar_and_pscalar_are_different_classes_in_same_algebra():
    """Even within a single algebra, scalar and pseudoscalar are distinct types.

    Cross-grade addition between them IS defined in C++ — it produces a
    multivector. For ega3d, scalar3d + pscalar3d yields mvec3d_e (the even
    multivector containing both grade 0 and grade 3 — wait, in ega3d the
    pseudoscalar is grade 3, which is odd, so the result is mvec3d, not
    mvec3d_e). For pga2dp/3dp, scalar*dp + pscalar*dp specifically yields
    DualNum*dp (see ga_dualnum*dp.hpp); other algebras compose into mvec.
    """
    s3 = ga_py.ega.scalar3d(1.0)
    p3 = ga_py.ega.pscalar3d(1.0)
    assert type(s3) is not type(p3)  # type-tag distinction preserved
    result = s3 + p3
    # In ega3d the pscalar is grade 3 (odd), so the result is mvec3d
    # (mixed even+odd). Not mvec3d_e.
    assert isinstance(result, ga_py.ega.mvec3d)


def test_pga_scalar_plus_pscalar_is_dualnum():
    """In pga2dp / pga3dp, scalar + pscalar specifically lands in DualNum*dp
    (the disambiguated canonical return type chosen for this combination)."""
    s = ga_py.pga.scalar2dp(1.0)
    p = ga_py.pga.pscalar2dp(2.0)
    assert isinstance(s + p, ga_py.pga.dualnum2dp)
    assert isinstance(p + s, ga_py.pga.dualnum2dp)
    assert isinstance(s - p, ga_py.pga.dualnum2dp)
    s3 = ga_py.pga.scalar3dp(1.0)
    p3 = ga_py.pga.pscalar3dp(2.0)
    assert isinstance(s3 + p3, ga_py.pga.dualnum3dp)
    assert isinstance(p3 - s3, ga_py.pga.dualnum3dp)


def test_same_type_addition():
    a = ga_py.ega.scalar3d(3.0)
    b = ga_py.ega.scalar3d(4.0)
    s = a + b
    assert isinstance(s, ga_py.ega.scalar3d)
    assert float(s) == 7.0


def test_same_type_subtraction():
    # pscalar2dp lives in pga (ends in 'dp'), not ega — partitioning check.
    a = ga_py.pga.pscalar2dp(10.0)
    b = ga_py.pga.pscalar2dp(3.0)
    d = a - b
    assert float(d) == 7.0


def test_unary_minus():
    s = ga_py.ega.scalar3d(5.0)
    n = -s
    assert isinstance(n, ga_py.ega.scalar3d)
    assert float(n) == -5.0


def test_scalar_times_double():
    s = ga_py.ega.scalar3d(2.0)
    r = s * 3.0
    assert float(r) == 6.0
    r2 = 4.0 * s
    assert float(r2) == 8.0


def test_scalar_divided_by_double():
    s = ga_py.ega.scalar3d(6.0)
    r = s / 2.0
    assert isinstance(r, ga_py.ega.scalar3d)
    assert float(r) == 3.0


def test_pscalar_arithmetic():
    """Pseudoscalar arithmetic works the same as scalar arithmetic."""
    p = ga_py.pga.pscalar3dp(2.0)
    q = ga_py.pga.pscalar3dp(3.0)
    assert float(p + q) == 5.0
    assert float(-p) == -2.0
    assert float(p * 2.5) == 5.0
