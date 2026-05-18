"""T1 --- Constants verification.

Every named constant exposed by ga_py.ega / ga_py.pga must hold the value
declared in ga/ga_usr_consts.hpp. This file enumerates all 120 currently
bound constants and verifies their components, plus the alias and cross-
relation invariants spelled out in the C++ source.

Source of truth: ga/ga_usr_consts.hpp.
"""

from __future__ import annotations

import pytest

import ga_py
from conftest import components

# --------------------------------------------------------------------------- #
# Expected values --- one entry per bound constant
# --------------------------------------------------------------------------- #
# Format: (submodule, name, type_name, expected_components_tuple)
# mvec basis orderings (from ga/detail/type_t/ga_mvec*.hpp):
#   mvec2d     (4): 1, e1, e2, e12
#   mvec2d_e   (2): 1, e12
#   mvec3d     (8): 1, e1, e2, e3, e23, e31, e12, e123
#   mvec3d_e   (4): 1, e23, e31, e12
#   mvec3d_u   (4): e1, e2, e3, e123
#   mvec2dp    (8): 1, e1, e2, e3, e31, e32, e12, e321
#   mvec2dp_e  (4): 1, e31, e32, e12
#   mvec2dp_u  (4): e1, e2, e3, e321
#   mvec3dp   (16): 1, e1, e2, e3, e4, e41, e42, e43, e23, e31, e12,
#                   e423, e431, e412, e321, e1234
#   mvec3dp_e  (8): 1, e41, e42, e43, e23, e31, e12, e1234
#   mvec3dp_u  (8): e1, e2, e3, e4, e423, e431, e412, e321
#
# Aliases collapse to identical components (verified separately below).

# Format used below: (submod, name, type_name, expected_components)
EXPECTED: list[tuple[str, str, str, tuple[float, ...]]] = [
    # ---------- EGA 2D ----------
    ("ega", "e1_2d", "vec2d", (1.0, 0.0)),
    ("ega", "e2_2d", "vec2d", (0.0, 1.0)),
    ("ega", "x_dir_2d", "vec2d", (1.0, 0.0)),
    ("ega", "y_dir_2d", "vec2d", (0.0, 1.0)),
    ("ega", "O_2d", "vec2d", (0.0, 0.0)),
    ("ega", "e1_2d_mv", "mvec2d", (0.0, 1.0, 0.0, 0.0)),
    ("ega", "e2_2d_mv", "mvec2d", (0.0, 0.0, 1.0, 0.0)),
    ("ega", "I_2d", "pscalar2d", (1.0,)),
    ("ega", "e12_2d", "pscalar2d", (1.0,)),
    ("ega", "I_2d_mv", "mvec2d", (0.0, 0.0, 0.0, 1.0)),
    ("ega", "I_2d_mv_e", "mvec2d_e", (0.0, 1.0)),
    ("ega", "one_2d", "scalar2d", (1.0,)),
    ("ega", "one_2d_mv", "mvec2d", (1.0, 0.0, 0.0, 0.0)),
    ("ega", "one_2d_mv_e", "mvec2d_e", (1.0, 0.0)),
    # ---------- EGA 3D ----------
    ("ega", "e1_3d", "vec3d", (1.0, 0.0, 0.0)),
    ("ega", "e2_3d", "vec3d", (0.0, 1.0, 0.0)),
    ("ega", "e3_3d", "vec3d", (0.0, 0.0, 1.0)),
    ("ega", "x_dir_3d", "vec3d", (1.0, 0.0, 0.0)),
    ("ega", "y_dir_3d", "vec3d", (0.0, 1.0, 0.0)),
    ("ega", "z_dir_3d", "vec3d", (0.0, 0.0, 1.0)),
    ("ega", "O_3d", "vec3d", (0.0, 0.0, 0.0)),
    ("ega", "e1_3d_mv", "mvec3d", (0, 1, 0, 0, 0, 0, 0, 0)),
    ("ega", "e2_3d_mv", "mvec3d", (0, 0, 1, 0, 0, 0, 0, 0)),
    ("ega", "e3_3d_mv", "mvec3d", (0, 0, 0, 1, 0, 0, 0, 0)),
    ("ega", "e23_3d", "bivec3d", (1.0, 0.0, 0.0)),
    ("ega", "e31_3d", "bivec3d", (0.0, 1.0, 0.0)),
    ("ega", "e12_3d", "bivec3d", (0.0, 0.0, 1.0)),
    ("ega", "yz_3d", "bivec3d", (1.0, 0.0, 0.0)),
    ("ega", "zx_3d", "bivec3d", (0.0, 1.0, 0.0)),
    ("ega", "xy_3d", "bivec3d", (0.0, 0.0, 1.0)),
    ("ega", "e23_3d_mv", "mvec3d", (0, 0, 0, 0, 1, 0, 0, 0)),
    ("ega", "e31_3d_mv", "mvec3d", (0, 0, 0, 0, 0, 1, 0, 0)),
    ("ega", "e12_3d_mv", "mvec3d", (0, 0, 0, 0, 0, 0, 1, 0)),
    ("ega", "e23_3d_mv_e", "mvec3d_e", (0, 1, 0, 0)),
    ("ega", "e31_3d_mv_e", "mvec3d_e", (0, 0, 1, 0)),
    ("ega", "e12_3d_mv_e", "mvec3d_e", (0, 0, 0, 1)),
    ("ega", "I_3d", "pscalar3d", (1.0,)),
    ("ega", "e123_3d", "pscalar3d", (1.0,)),
    ("ega", "I_3d_mv", "mvec3d", (0, 0, 0, 0, 0, 0, 0, 1)),
    ("ega", "I_3d_mv_u", "mvec3d_u", (0, 0, 0, 1)),
    ("ega", "one_3d", "scalar3d", (1.0,)),
    ("ega", "one_3d_mv", "mvec3d", (1, 0, 0, 0, 0, 0, 0, 0)),
    ("ega", "one_3d_mv_e", "mvec3d_e", (1, 0, 0, 0)),
    # ---------- PGA 2dp ----------
    ("pga", "e1_2dp", "vec2dp", (1.0, 0.0, 0.0)),
    ("pga", "e2_2dp", "vec2dp", (0.0, 1.0, 0.0)),
    ("pga", "e3_2dp", "vec2dp", (0.0, 0.0, 1.0)),
    ("pga", "O_2dp", "vec2dp", (0.0, 0.0, 1.0)),
    ("pga", "x_dir_2dp", "vec2dp", (1.0, 0.0, 0.0)),
    ("pga", "y_dir_2dp", "vec2dp", (0.0, 1.0, 0.0)),
    ("pga", "e1_2dp_mv", "mvec2dp", (0, 1, 0, 0, 0, 0, 0, 0)),
    ("pga", "e2_2dp_mv", "mvec2dp", (0, 0, 1, 0, 0, 0, 0, 0)),
    ("pga", "e3_2dp_mv", "mvec2dp", (0, 0, 0, 1, 0, 0, 0, 0)),
    ("pga", "e31_2dp", "bivec2dp", (1.0, 0.0, 0.0)),
    ("pga", "e32_2dp", "bivec2dp", (0.0, 1.0, 0.0)),
    ("pga", "e12_2dp", "bivec2dp", (0.0, 0.0, 1.0)),
    ("pga", "x_axis_2dp", "bivec2dp", (1.0, 0.0, 0.0)),
    ("pga", "y_axis_2dp", "bivec2dp", (0.0, 1.0, 0.0)),
    # H_2dp is initialized as `-e12_2dp` -> components are (-0.0, -0.0, -1.0);
    # tuple equality treats -0.0 == 0.0, so we compare against (0, 0, -1).
    ("pga", "H_2dp", "bivec2dp", (0.0, 0.0, -1.0)),
    ("pga", "e31_2dp_mv", "mvec2dp", (0, 0, 0, 0, 1, 0, 0, 0)),
    ("pga", "e32_2dp_mv", "mvec2dp", (0, 0, 0, 0, 0, 1, 0, 0)),
    ("pga", "e12_2dp_mv", "mvec2dp", (0, 0, 0, 0, 0, 0, 1, 0)),
    ("pga", "e31_2dp_mv_e", "mvec2dp_e", (0, 1, 0, 0)),
    ("pga", "e32_2dp_mv_e", "mvec2dp_e", (0, 0, 1, 0)),
    ("pga", "e12_2dp_mv_e", "mvec2dp_e", (0, 0, 0, 1)),
    ("pga", "I_2dp", "pscalar2dp", (1.0,)),
    ("pga", "e321_2dp", "pscalar2dp", (1.0,)),
    ("pga", "I_2dp_mv", "mvec2dp", (0, 0, 0, 0, 0, 0, 0, 1)),
    ("pga", "I_2dp_mv_u", "mvec2dp_u", (0, 0, 0, 1)),
    ("pga", "one_2dp", "scalar2dp", (1.0,)),
    ("pga", "one_2dp_mv", "mvec2dp", (1, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "one_2dp_mv_e", "mvec2dp_e", (1, 0, 0, 0)),
    # ---------- PGA 3dp ----------
    ("pga", "e1_3dp", "vec3dp", (1.0, 0.0, 0.0, 0.0)),
    ("pga", "e2_3dp", "vec3dp", (0.0, 1.0, 0.0, 0.0)),
    ("pga", "e3_3dp", "vec3dp", (0.0, 0.0, 1.0, 0.0)),
    ("pga", "e4_3dp", "vec3dp", (0.0, 0.0, 0.0, 1.0)),
    ("pga", "O_3dp", "vec3dp", (0.0, 0.0, 0.0, 1.0)),
    ("pga", "x_dir_3dp", "vec3dp", (1.0, 0.0, 0.0, 0.0)),
    ("pga", "y_dir_3dp", "vec3dp", (0.0, 1.0, 0.0, 0.0)),
    ("pga", "z_dir_3dp", "vec3dp", (0.0, 0.0, 1.0, 0.0)),
    ("pga", "e1_3dp_mv", "mvec3dp", (0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e2_3dp_mv", "mvec3dp", (0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e3_3dp_mv", "mvec3dp", (0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e4_3dp_mv", "mvec3dp", (0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e41_3dp", "bivec3dp", (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
    ("pga", "e42_3dp", "bivec3dp", (0.0, 1.0, 0.0, 0.0, 0.0, 0.0)),
    ("pga", "e43_3dp", "bivec3dp", (0.0, 0.0, 1.0, 0.0, 0.0, 0.0)),
    ("pga", "e23_3dp", "bivec3dp", (0.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
    ("pga", "e31_3dp", "bivec3dp", (0.0, 0.0, 0.0, 0.0, 1.0, 0.0)),
    ("pga", "e12_3dp", "bivec3dp", (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
    ("pga", "x_axis_3dp", "bivec3dp", (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
    ("pga", "y_axis_3dp", "bivec3dp", (0.0, 1.0, 0.0, 0.0, 0.0, 0.0)),
    ("pga", "z_axis_3dp", "bivec3dp", (0.0, 0.0, 1.0, 0.0, 0.0, 0.0)),
    ("pga", "e41_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e42_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e43_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e23_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "e31_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)),
    ("pga", "e12_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)),
    ("pga", "e23_3dp_mv_e", "mvec3dp_e", (0, 0, 0, 0, 1, 0, 0, 0)),
    ("pga", "e31_3dp_mv_e", "mvec3dp_e", (0, 0, 0, 0, 0, 1, 0, 0)),
    ("pga", "e12_3dp_mv_e", "mvec3dp_e", (0, 0, 0, 0, 0, 0, 1, 0)),
    ("pga", "e423_3dp", "trivec3dp", (1.0, 0.0, 0.0, 0.0)),
    ("pga", "e431_3dp", "trivec3dp", (0.0, 1.0, 0.0, 0.0)),
    ("pga", "e412_3dp", "trivec3dp", (0.0, 0.0, 1.0, 0.0)),
    ("pga", "e321_3dp", "trivec3dp", (0.0, 0.0, 0.0, 1.0)),
    ("pga", "yz_3dp", "trivec3dp", (1.0, 0.0, 0.0, 0.0)),
    ("pga", "zx_3dp", "trivec3dp", (0.0, 1.0, 0.0, 0.0)),
    ("pga", "xy_3dp", "trivec3dp", (0.0, 0.0, 1.0, 0.0)),
    ("pga", "H_3dp", "trivec3dp", (0.0, 0.0, 0.0, 1.0)),
    ("pga", "e423_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)),
    ("pga", "e431_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)),
    ("pga", "e412_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)),
    ("pga", "e321_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)),
    ("pga", "I_3dp", "pscalar3dp", (1.0,)),
    ("pga", "e1234_3dp", "pscalar3dp", (1.0,)),
    ("pga", "I_3dp_mv", "mvec3dp", (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)),
    ("pga", "I_3dp_mv_e", "mvec3dp_e", (0, 0, 0, 0, 0, 0, 0, 1)),
    ("pga", "one_3dp", "scalar3dp", (1.0,)),
    ("pga", "one_3dp_mv", "mvec3dp", (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    ("pga", "one_3dp_mv_e", "mvec3dp_e", (1, 0, 0, 0, 0, 0, 0, 0)),
]


# --------------------------------------------------------------------------- #
# Per-constant value verification (one test per constant)
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "submod,name,type_name,expected",
    EXPECTED,
    ids=[f"{s}.{n}" for s, n, _, _ in EXPECTED],
)
def test_constant_value(submod, name, type_name, expected):
    """Each bound constant has the components declared in ga_usr_consts.hpp."""
    obj = getattr(getattr(ga_py, submod), name)
    assert (
        type(obj).__name__ == type_name
    ), f"expected type {type_name}, got {type(obj).__name__}"
    actual = components(obj)
    assert actual == pytest.approx(
        expected, abs=0.0
    ), f"{submod}.{name}: expected {expected}, got {actual}"


# --------------------------------------------------------------------------- #
# Alias relationships --- these constants are declared as `auto const x = y;`
# in ga_usr_consts.hpp and must therefore be equal as Python objects too.
# --------------------------------------------------------------------------- #

ALIASES = [
    # EGA 2D
    ("ega", "x_dir_2d", "e1_2d"),
    ("ega", "y_dir_2d", "e2_2d"),
    ("ega", "e12_2d", "I_2d"),
    # EGA 3D
    ("ega", "x_dir_3d", "e1_3d"),
    ("ega", "y_dir_3d", "e2_3d"),
    ("ega", "z_dir_3d", "e3_3d"),
    ("ega", "yz_3d", "e23_3d"),
    ("ega", "zx_3d", "e31_3d"),
    ("ega", "xy_3d", "e12_3d"),
    ("ega", "e123_3d", "I_3d"),
    # PGA 2dp
    ("pga", "O_2dp", "e3_2dp"),
    ("pga", "x_dir_2dp", "e1_2dp"),
    ("pga", "y_dir_2dp", "e2_2dp"),
    ("pga", "x_axis_2dp", "e31_2dp"),
    ("pga", "y_axis_2dp", "e32_2dp"),
    ("pga", "e321_2dp", "I_2dp"),
    # PGA 3dp
    ("pga", "O_3dp", "e4_3dp"),
    ("pga", "x_dir_3dp", "e1_3dp"),
    ("pga", "y_dir_3dp", "e2_3dp"),
    ("pga", "z_dir_3dp", "e3_3dp"),
    ("pga", "x_axis_3dp", "e41_3dp"),
    ("pga", "y_axis_3dp", "e42_3dp"),
    ("pga", "z_axis_3dp", "e43_3dp"),
    ("pga", "yz_3dp", "e423_3dp"),
    ("pga", "zx_3dp", "e431_3dp"),
    ("pga", "xy_3dp", "e412_3dp"),
    ("pga", "H_3dp", "e321_3dp"),
    ("pga", "e1234_3dp", "I_3dp"),
]


@pytest.mark.parametrize(
    "submod,alias,target",
    ALIASES,
    ids=[f"{s}.{a}=={s}.{t}" for s, a, t in ALIASES],
)
def test_alias_equals_target(submod, alias, target):
    """Aliases declared via `auto const a = b;` in C++ must compare equal in Python."""
    sub = getattr(ga_py, submod)
    a, b = getattr(sub, alias), getattr(sub, target)
    assert a == b, f"{submod}.{alias} should equal {submod}.{target}"


# --------------------------------------------------------------------------- #
# A handful of mathematical sanity checks --- catches dispatch regressions
# on basis values themselves
# --------------------------------------------------------------------------- #


def test_ega3d_basis_orthonormality():
    """e_i . e_j == delta_ij for the EGA3D orthonormal basis."""
    basis = [ga_py.ega.e1_3d, ga_py.ega.e2_3d, ga_py.ega.e3_3d]
    for i, ei in enumerate(basis):
        for j, ej in enumerate(basis):
            expected = 1.0 if i == j else 0.0
            assert (
                float(ga_py.ega.dot(ei, ej)) == expected
            ), f"dot(e{i+1}_3d, e{j+1}_3d) should be {expected}"


def test_ega3d_wedge_of_basis_vectors_yields_basis_bivectors():
    """e_i ^ e_j matches the named EGA3D basis bivectors."""
    e1, e2, e3 = ga_py.ega.e1_3d, ga_py.ega.e2_3d, ga_py.ega.e3_3d
    assert (e2 ^ e3) == ga_py.ega.e23_3d
    assert (e3 ^ e1) == ga_py.ega.e31_3d
    assert (e1 ^ e2) == ga_py.ega.e12_3d


def test_ega2d_wedge_of_basis_vectors_yields_pseudoscalar():
    """e1_2d ^ e2_2d == I_2d == e12_2d."""
    e1, e2 = ga_py.ega.e1_2d, ga_py.ega.e2_2d
    assert (e1 ^ e2) == ga_py.ega.I_2d


def test_pga2dp_basis_bivector_definitions():
    """e31 = e3 ^ e1, e32 = e3 ^ e2, e12 = e1 ^ e2 in PGA2dp."""
    e1, e2, e3 = ga_py.pga.e1_2dp, ga_py.pga.e2_2dp, ga_py.pga.e3_2dp
    assert (e3 ^ e1) == ga_py.pga.e31_2dp
    assert (e3 ^ e2) == ga_py.pga.e32_2dp
    assert (e1 ^ e2) == ga_py.pga.e12_2dp


def test_pga3dp_pseudoscalar_factorization():
    """I_3dp = e1 ^ e2 ^ e3 ^ e4 = e1234_3dp."""
    e1 = ga_py.pga.e1_3dp
    e2 = ga_py.pga.e2_3dp
    e3 = ga_py.pga.e3_3dp
    e4 = ga_py.pga.e4_3dp
    # wedge into bivec, then trivec, then pscalar
    b = e1 ^ e2
    t = b ^ e3
    ps = t ^ e4
    assert ps == ga_py.pga.I_3dp


# --------------------------------------------------------------------------- #
# one_* multiplicative-identity checks --- one_* is the geometric-product unit
# (`*` is the geometric product). Covered for all three forms: the scalar
# one_*, the full multivector one_*_mv, and the even-grade one_*_mv_e.
# --------------------------------------------------------------------------- #


def test_ega2d_one_is_geometric_product_identity():
    """one_2d / one_2d_mv / one_2d_mv_e are the EGA2d geometric-product unit."""
    sub = ga_py.ega
    v, mv, mv_e = sub.e1_2d, sub.e1_2d_mv, sub.I_2d_mv_e
    assert sub.one_2d * v == v and v * sub.one_2d == v
    assert sub.one_2d_mv * mv == mv and mv * sub.one_2d_mv == mv
    assert sub.one_2d_mv_e * mv_e == mv_e and mv_e * sub.one_2d_mv_e == mv_e


def test_ega3d_one_is_geometric_product_identity():
    """one_3d / one_3d_mv / one_3d_mv_e are the EGA3d geometric-product unit."""
    sub = ga_py.ega
    v, mv, mv_e = sub.e1_3d, sub.e1_3d_mv, sub.e23_3d_mv_e
    assert sub.one_3d * v == v and v * sub.one_3d == v
    assert sub.one_3d_mv * mv == mv and mv * sub.one_3d_mv == mv
    assert sub.one_3d_mv_e * mv_e == mv_e and mv_e * sub.one_3d_mv_e == mv_e


def test_pga2dp_one_is_geometric_product_identity():
    """one_2dp / one_2dp_mv / one_2dp_mv_e are the PGA2dp geometric-product unit."""
    sub = ga_py.pga
    v, mv, mv_e = sub.e1_2dp, sub.e1_2dp_mv, sub.e31_2dp_mv_e
    assert sub.one_2dp * v == v and v * sub.one_2dp == v
    assert sub.one_2dp_mv * mv == mv and mv * sub.one_2dp_mv == mv
    assert sub.one_2dp_mv_e * mv_e == mv_e and mv_e * sub.one_2dp_mv_e == mv_e


def test_pga3dp_one_is_geometric_product_identity():
    """one_3dp / one_3dp_mv / one_3dp_mv_e are the PGA3dp geometric-product unit."""
    sub = ga_py.pga
    v, mv, mv_e = sub.e1_3dp, sub.e1_3dp_mv, sub.e23_3dp_mv_e
    assert sub.one_3dp * v == v and v * sub.one_3dp == v
    assert sub.one_3dp_mv * mv == mv and mv * sub.one_3dp_mv == mv
    assert sub.one_3dp_mv_e * mv_e == mv_e and mv_e * sub.one_3dp_mv_e == mv_e


# --------------------------------------------------------------------------- #
# I_* regressive-identity checks --- in PGA the pseudoscalar I_* is the unit
# of the regressive geometric product rgpr() (dual to one_* being the gpr
# unit). Covered for all three forms: the pseudoscalar I_*, the full
# multivector I_*_mv, and the parity-grade form (I_2dp_mv_u, since I_2dp is
# odd grade; I_3dp_mv_e, since I_3dp is even grade).
# --------------------------------------------------------------------------- #


def test_pga2dp_pseudoscalar_is_regressive_product_identity():
    """I_2dp / I_2dp_mv / I_2dp_mv_u are the PGA2dp regressive-product unit."""
    p = ga_py.pga
    v = p.e1_2dp
    mv = p.mvec2dp(
        p.scalar2dp(2.0), p.vec2dp(1, 2, 3), p.bivec2dp(4, 5, 6), p.pscalar2dp(7.0)
    )
    mv_u = p.mvec2dp_u(p.vec2dp(1, 2, 3), p.pscalar2dp(7.0))
    assert p.rgpr(p.I_2dp, v) == v and p.rgpr(v, p.I_2dp) == v
    assert p.rgpr(p.I_2dp_mv, mv) == mv and p.rgpr(mv, p.I_2dp_mv) == mv
    assert p.rgpr(p.I_2dp_mv_u, mv_u) == mv_u and p.rgpr(mv_u, p.I_2dp_mv_u) == mv_u


def test_pga3dp_pseudoscalar_is_regressive_product_identity():
    """I_3dp / I_3dp_mv / I_3dp_mv_e are the PGA3dp regressive-product unit."""
    p = ga_py.pga
    v = p.e1_3dp
    mv = p.mvec3dp(
        p.scalar3dp(2.0),
        p.vec3dp(1, 2, 3, 4),
        p.bivec3dp(1, 2, 3, 4, 5, 6),
        p.trivec3dp(1, 2, 3, 4),
        p.pscalar3dp(7.0),
    )
    mv_e = p.mvec3dp_e(
        p.scalar3dp(2.0), p.bivec3dp(1, 2, 3, 4, 5, 6), p.pscalar3dp(7.0)
    )
    assert p.rgpr(p.I_3dp, v) == v and p.rgpr(v, p.I_3dp) == v
    assert p.rgpr(p.I_3dp_mv, mv) == mv and p.rgpr(mv, p.I_3dp_mv) == mv
    assert p.rgpr(p.I_3dp_mv_e, mv_e) == mv_e and p.rgpr(mv_e, p.I_3dp_mv_e) == mv_e


# --------------------------------------------------------------------------- #
# Coverage / drift detection --- proves we test every bound constant
# --------------------------------------------------------------------------- #


def _all_bound_constants() -> set[tuple[str, str]]:
    """Enumerate every const-like attribute on ga_py.ega and ga_py.pga."""
    found: set[tuple[str, str]] = set()
    for sub_name in ("ega", "pga"):
        sub = getattr(ga_py, sub_name)
        for n in dir(sub):
            if n.startswith("_"):
                continue
            attr = getattr(sub, n)
            if isinstance(attr, type) or callable(attr):
                continue
            found.add((sub_name, n))
    return found


def test_expected_table_covers_every_bound_constant():
    """If a constant is added to C++ and bindings, it must be added here too."""
    bound = _all_bound_constants()
    covered = {(s, n) for s, n, _, _ in EXPECTED}
    missing = bound - covered
    extra = covered - bound
    assert not missing, f"bound but not tested: {sorted(missing)}"
    assert not extra, f"tested but no longer bound: {sorted(extra)}"


def test_total_bound_constants_is_120():
    """Stable count --- changes here flag intentional or accidental drift.

    120 = the 108 previously bound + 12 `one_*` scalar/multivector identities
    (one_2d{,_mv,_mv_e}, one_3d{,_mv,_mv_e}, one_2dp{,_mv,_mv_e},
    one_3dp{,_mv,_mv_e}) added in ga/ga_usr_consts.hpp.
    """
    assert len(_all_bound_constants()) == 120


def test_h_2dp_value_matches_negated_e12_2dp():
    """H_2dp is declared as `-e12_2dp` in ga_usr_consts.hpp; verify."""
    sub = ga_py.pga
    assert sub.H_2dp.x == 0.0
    assert sub.H_2dp.y == 0.0
    assert sub.H_2dp.z == -1.0
    assert sub.e12_2dp.z == 1.0  # invariant: H_2dp.z == -e12_2dp.z
