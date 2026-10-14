// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_arith_check.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// STA contracts: the arithmetic of every value type, and the operations of
// ga_sta4ds_ops*.hpp that the topical tests did not reach
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("STA contracts")
{
    TEST_CASE("sta4ds: arithmetic contract of every value type")
    {
        arith_check::check<scalar4ds>("scalar4ds");
        arith_check::check<vec4ds>("vec4ds");
        arith_check::check<bivec4ds>("bivec4ds");
        arith_check::check<trivec4ds>("trivec4ds");
        arith_check::check<pscalar4ds>("pscalar4ds");
        arith_check::check<mvec4ds_e>("mvec4ds_e");
        arith_check::check<mvec4ds_u>("mvec4ds_u");
        arith_check::check<mvec4ds>("mvec4ds");
    }

    TEST_CASE("sta4ds: the grade of every value type")
    {
        CHECK(gr(scalar4ds{}) == 0);
        CHECK(gr(vec4ds{}) == 1);
        CHECK(gr(bivec4ds{}) == 2);
        CHECK(gr(trivec4ds{}) == 3);
        CHECK(gr(pscalar4ds{}) == 4);
    }

    TEST_CASE("sta4ds: normalize scales a blade to nrm_sq = +/-1 by its causal kind")
    {
        // timelike -> +1, spacelike -> -1, lightlike (null) is returned unchanged
        vec4ds const t(0.0, 0.0, 0.0, 2.0);
        vec4ds const s(3.0, 0.0, 0.0, 0.0);
        vec4ds const n(1.0, 0.0, 0.0, 1.0);
        CHECK(normalize(t) == vec4ds(0.0, 0.0, 0.0, 1.0));
        CHECK(value_t(nrm_sq(normalize(t))) == doctest::Approx(1.0));
        CHECK(normalize(s) == vec4ds(1.0, 0.0, 0.0, 0.0));
        CHECK(value_t(nrm_sq(normalize(s))) == doctest::Approx(-1.0));
        CHECK(normalize(n) == n);

        // a bivector and a trivector: unit |nrm_sq|, same direction
        bivec4ds const B(0.0, 0.0, 0.0, 0.0, 0.0, 4.0); // 4 g12
        CHECK(normalize(B) == bivec4ds(0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
        CHECK(std::abs(value_t(nrm_sq(normalize(
                  bivec4ds(1.0, 2.0, 0.5, -1.0, 0.0, 3.0))))) == doctest::Approx(1.0));
        trivec4ds const T(0.0, 0.0, 0.0, 5.0); // 5 g123
        CHECK(normalize(T) == trivec4ds(0.0, 0.0, 0.0, 1.0));

        // a scalar normalizes to its sign
        CHECK(normalize(scalar4ds(-3.0)) == scalar4ds(-1.0));
        CHECK(normalize(scalar4ds(0.25)) == scalar4ds(1.0));
    }

    TEST_CASE("sta4ds: transform of a full multivector transforms each grade")
    {
        auto const R = get_rotor(g12_4ds, deg2rad(30.0)) * get_boost(g14_4ds, 0.4);
        vec4ds const v(1.0, -2.0, 0.5, 3.0);
        bivec4ds const B(1.0, 2.0, 0.5, -1.0, 0.0, 3.0);
        trivec4ds const T(2.0, 0.0, -1.0, 0.5);
        mvec4ds const X = mvec4ds(scalar4ds(2.0)) + mvec4ds(v) + mvec4ds(B) + mvec4ds(T) +
                          mvec4ds(pscalar4ds(-1.5));
        // a unit rotor leaves the scalar and the pseudoscalar alone
        CHECK(transform(X, R) == mvec4ds(scalar4ds(2.0)) + mvec4ds(transform(v, R)) +
                                     mvec4ds(transform(B, R)) + mvec4ds(transform(T, R)) +
                                     mvec4ds(pscalar4ds(-1.5)));
    }

    TEST_CASE("sta4ds: is_congruent asks for the same subspace, is_close for the value")
    {
        // any two non-zero scalars (pseudoscalars) span the same subspace
        CHECK(is_congruent(scalar4ds(2.0), scalar4ds(-7.0)));
        CHECK(is_congruent(pscalar4ds(2.0), pscalar4ds(-7.0)));

        vec4ds const v(1.0, -2.0, 0.5, 3.0);
        CHECK(is_congruent(v, -4.0 * v));
        CHECK_FALSE(is_congruent(v, vec4ds(1.0, -2.0, 0.5, 3.1)));
        bivec4ds const B(1.0, 2.0, 0.5, -1.0, 0.0, 3.0);
        CHECK(is_congruent(B, 0.5 * B));
        CHECK_FALSE(is_congruent(B, bivec4ds(1.0, 2.0, 0.5, -1.0, 0.1, 3.0)));
        trivec4ds const T(2.0, 0.0, -1.0, 0.5);
        CHECK(is_congruent(T, -T));
        CHECK_FALSE(is_congruent(T, trivec4ds(2.0, 0.0, -1.0, 0.6)));

        CHECK(is_close(scalar4ds(3.0), scalar4ds(3.0 + 3.0e-14)));
        CHECK_FALSE(is_close(scalar4ds(3.0), scalar4ds(-3.0)));
        CHECK(is_close(pscalar4ds(0.5), pscalar4ds(0.5 + 1.0e-15)));
        CHECK_FALSE(is_close(pscalar4ds(0.5), pscalar4ds(0.6)));
    }
}
