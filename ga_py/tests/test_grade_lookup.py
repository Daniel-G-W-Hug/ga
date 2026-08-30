"""T2 --- Grade lookup verification.

`gr()` returns the grade (0..n) of a graded GA element. `rgr()` returns
the regressive (co-)grade, or anti-grade, `algebra_dim - grade`.

By design, neither is defined for every type:
  * `gr()` is defined only for graded types (scalar, vector, bivector,
    trivector, pscalar). Multivectors are multi-graded; use their
    per-grade extractor methods `.gr0()`, `.gr1()`, ... `.grN()` instead.
  * `rgr()` is defined for PGA and CGA only. The criterion is not the
    metric -- anti-grade is just `n - grade` and needs none -- but the
    ANTI-PRODUCT family: pga2dp/3dp and cga2dc/3dc carry rgpr, rdot and
    rcmt, and the anti-grade is their grade bookkeeping. ega2d/3d and
    sta4ds carry only rwdg, and define no rgr. (sta4ds did until
    2026-08-30, with no caller and no test in any language; it was
    removed to make the rule uniform.)

This file tests only the cases where each function is defined.
"""

from __future__ import annotations

import pytest

import ga_py


# --------------------------------------------------------------------------- #
# (submodule, type_name, ctor_args, expected_grade) for every graded type
# --------------------------------------------------------------------------- #
# Algebra dimensions (used for rgr expectations below):
#   ega2d=2, ega3d=3, sta4ds=4, pga2dp=3, pga3dp=4
# (PGA's "dimension" here is the highest grade, which equals the dimension
#  of its underlying linear space: e321 for 2dp = grade 3, e1234 for 3dp = grade 4.)

GRADED_TYPES: list[tuple[str, str, tuple, int]] = [
    # EGA 2D --- algebra dim 2
    ("ega", "scalar2d",   (1.0,),                  0),
    ("ega", "vec2d",      (1.0, 2.0),              1),
    ("ega", "pscalar2d",  (1.0,),                  2),
    # EGA 3D --- algebra dim 3
    ("ega", "scalar3d",   (1.0,),                  0),
    ("ega", "vec3d",      (1.0, 2.0, 3.0),         1),
    ("ega", "bivec3d",    (1.0, 2.0, 3.0),         2),
    ("ega", "pscalar3d",  (1.0,),                  3),
    # STA 4ds --- G(1,3,0), algebra dim 4 (g1234 is the pseudoscalar)
    ("sta", "scalar4ds",  (1.0,),                  0),
    ("sta", "vec4ds",     (1.0, 2.0, 3.0, 4.0),    1),
    ("sta", "bivec4ds",   (1.0, 2.0, 3.0, 4.0, 5.0, 6.0), 2),
    ("sta", "trivec4ds",  (1.0, 2.0, 3.0, 4.0),    3),
    ("sta", "pscalar4ds", (1.0,),                  4),
    # PGA 2dp --- algebra dim 3 (e321 is the pseudoscalar)
    ("pga", "scalar2dp",  (1.0,),                  0),
    ("pga", "vec2dp",     (1.0, 2.0, 3.0),         1),
    ("pga", "bivec2dp",   (1.0, 2.0, 3.0),         2),
    ("pga", "pscalar2dp", (1.0,),                  3),
    # PGA 3dp --- algebra dim 4 (e1234 is the pseudoscalar)
    ("pga", "scalar3dp",  (1.0,),                  0),
    ("pga", "vec3dp",     (1.0, 2.0, 3.0, 4.0),    1),
    ("pga", "bivec3dp",   (1.0, 2.0, 3.0, 4.0, 5.0, 6.0), 2),
    ("pga", "trivec3dp",  (1.0, 2.0, 3.0, 4.0),    3),
    ("pga", "pscalar3dp", (1.0,),                  4),
]


def _make(submod: str, type_name: str, ctor_args: tuple):
    return getattr(getattr(ga_py, submod), type_name)(*ctor_args)


# --------------------------------------------------------------------------- #
# gr() --- must return the natural grade of every bound graded type
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize(
    "submod,type_name,ctor_args,expected_grade", GRADED_TYPES,
    ids=[f"{s}.{t}" for s, t, _, _ in GRADED_TYPES],
)
def test_gr_returns_correct_grade(submod, type_name, ctor_args, expected_grade):
    obj = _make(submod, type_name, ctor_args)
    assert ga_py.gr(obj) == expected_grade


# --------------------------------------------------------------------------- #
# rgr() --- must return algebra_dim - grade. Parametrized over PGA here; the CGA
# overloads are pinned on the C++ side instead (ga_cga2dc_test / ga_cga3dc_test).
# --------------------------------------------------------------------------- #
# pga2dp algebra_dim = 3, pga3dp algebra_dim = 4.

_ALGEBRA_DIM = {"2dp": 3, "3dp": 4}

PGA_GRADED = [
    (s, t, ca, g) for (s, t, ca, g) in GRADED_TYPES if s == "pga"
]


@pytest.mark.parametrize(
    "submod,type_name,ctor_args,grade", PGA_GRADED,
    ids=[f"{s}.{t}" for s, t, _, _ in PGA_GRADED],
)
def test_rgr_returns_algebra_dim_minus_grade(submod, type_name, ctor_args, grade):
    obj = _make(submod, type_name, ctor_args)
    suffix = type_name.removeprefix("scalar").removeprefix("vec").removeprefix("bivec") \
                      .removeprefix("trivec").removeprefix("pscalar")
    expected = _ALGEBRA_DIM[suffix] - grade
    assert ga_py.rgr(obj) == expected
