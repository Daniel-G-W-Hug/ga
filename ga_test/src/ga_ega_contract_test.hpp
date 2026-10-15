// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_arith_check.hpp"

#include <cmath>
#include <type_traits>
#include <vector>

/////////////////////////////////////////////////////////////////////////////////////////
// EGA contracts: the arithmetic of every value type, and the operations of
// ga_ega{2,3}d_ops*.hpp that the topical tests did not reach, each pinned against an
// independently known answer
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA contracts")
{
    TEST_CASE("ega2d: arithmetic contract of every value type")
    {
        arith_check::check<scalar2d>("scalar2d");
        arith_check::check<vec2d>("vec2d");
        arith_check::check<pscalar2d>("pscalar2d");
        arith_check::check<mvec2d_e>("mvec2d_e");
        arith_check::check<mvec2d>("mvec2d");
    }

    TEST_CASE("ega2d: an even multivector is built from (scalar, pseudoscalar) only")
    {
        static_assert(std::is_constructible_v<mvec2d_e, scalar2d, pscalar2d>);
        CHECK(mvec2d_e(scalar2d(1.0), pscalar2d(2.0)) == mvec2d_e(1.0, 2.0));
        // swapped parts, or parts of another algebra, do not compile
        static_assert(!std::is_constructible_v<mvec2d_e, pscalar2d, scalar2d>);
        static_assert(!std::is_constructible_v<mvec2d_e, scalar3d, pscalar3d>);
    }

    TEST_CASE("ega3d: arithmetic contract of every value type")
    {
        arith_check::check<scalar3d>("scalar3d");
        arith_check::check<vec3d>("vec3d");
        arith_check::check<bivec3d>("bivec3d");
        arith_check::check<pscalar3d>("pscalar3d");
        arith_check::check<mvec3d_e>("mvec3d_e");
        arith_check::check<mvec3d_u>("mvec3d_u");
        arith_check::check<mvec3d>("mvec3d");

        // the even multivector from a scalar and a pseudoscalar alone
        CHECK(mvec3d(scalar3d(2.0), pscalar3d(-3.0)) ==
              mvec3d(scalar3d(2.0)) + mvec3d(pscalar3d(-3.0)));
    }

    TEST_CASE("ega3d: angle between planes, and between a vector and a plane")
    {
        // the e12 plane and the e23 plane meet at a right angle; e12 against a plane
        // tilted by 30 deg about e1 encloses 30 deg
        value_t const a = deg2rad(30.0);
        bivec3d const tilted = wdg(e1_3d, vec3d(0.0, std::cos(a), std::sin(a)));
        CHECK(angle(e12_3d, e23_3d) == doctest::Approx(pi / 2));
        CHECK(angle(e12_3d, tilted) == doctest::Approx(a));
        CHECK(angle(e12_3d, e12_3d) == doctest::Approx(0.0));

        // a vector against a plane: the SIGNED angle between the line and the plane, the
        // same in both argument orders -- 0 in the plane, +pi/2 along its normal
        // cmpl(e12) = e3, -pi/2 against it, so the range is [-pi/2, pi/2]; reversing
        // the vector or the plane flips the sign
        vec3d const v(std::cos(a), 0.0, std::sin(a)); // 30 deg out of the e12 plane
        CHECK(angle(v, e12_3d) == doctest::Approx(a));
        CHECK(angle(e12_3d, v) == doctest::Approx(a));
        CHECK(angle(-v, e12_3d) == doctest::Approx(-a));
        CHECK(angle(e12_3d, -v) == doctest::Approx(-a));
        CHECK(angle(v, -e12_3d) == doctest::Approx(-a));
        CHECK(angle(e3_3d, e12_3d) == doctest::Approx(pi / 2));
        CHECK(angle(-e3_3d, e12_3d) == doctest::Approx(-pi / 2));
        CHECK(angle(e1_3d, e12_3d) == doctest::Approx(0.0));
        // scale-free, and pi/2 minus the angle to the plane's normal for every pair
        CHECK(angle(5.0 * v, 0.25 * e12_3d) == doctest::Approx(a));
        vec3d const w(0.3, -1.7, -0.6);
        bivec3d const B(-0.8, 0.4, 1.1);
        CHECK(angle(w, B) == doctest::Approx(pi / 2 - angle(w, vec3d(cmpl(B)))));
        CHECK(angle(-w, B) == doctest::Approx(pi / 2 - angle(-w, vec3d(cmpl(B)))));
    }

    TEST_CASE("ega3d: rotate_opt on a batch of bivectors is rotate, element-wise")
    {
        auto const R = get_rotor(e12_3d, deg2rad(40.0));
        std::vector<bivec3d> const B{e23_3d, e31_3d, e12_3d, bivec3d(1.0, -2.0, 0.5)};
        auto const rotated = rotate_opt(B, R);
        REQUIRE(rotated.size() == B.size());
        for (size_t i = 0; i < B.size(); ++i)
            CHECK(rotated[i] == rotate(B[i], R));
    }

    TEST_CASE("ega: is_close on scalars, pseudoscalars and even multivectors")
    {
        // a relative tolerance: a difference of 1e-14 is close at any scale, 1e-6 is not
        CHECK(is_close(scalar2d(3.0), scalar2d(3.0 + 3.0e-14)));
        CHECK_FALSE(is_close(scalar2d(3.0), scalar2d(3.0 + 3.0e-6)));
        CHECK(is_close(pscalar2d(-2.0e6), pscalar2d(-2.0e6 * (1.0 + 1.0e-14))));
        CHECK_FALSE(is_close(pscalar2d(1.0), pscalar2d(-1.0)));
        CHECK(is_close(mvec2d_e(1.0, -2.0), mvec2d_e(1.0, -2.0 + 1.0e-14)));
        CHECK_FALSE(is_close(mvec2d_e(1.0, -2.0), mvec2d_e(1.0, 2.0)));
        // a tolerance can be widened explicitly
        CHECK(is_close(mvec2d_e(1.0, -2.0), mvec2d_e(1.0, -2.001), 1.0e-3));

        CHECK(is_close(scalar3d(3.0), scalar3d(3.0 + 3.0e-14)));
        CHECK_FALSE(is_close(scalar3d(3.0), scalar3d(3.1)));
        CHECK(is_close(pscalar3d(0.5), pscalar3d(0.5 + 1.0e-15)));
        CHECK_FALSE(is_close(pscalar3d(0.5), pscalar3d(0.6)));
        mvec3d_e const R = get_rotor(e12_3d, 0.3);
        CHECK(is_close(R, R * (1.0 + 1.0e-14)));
        CHECK_FALSE(is_close(R, -R)); // the same rotation, but not the same value
    }

    TEST_CASE("ega2d: exp of a bivector angle is the unit complex number")
    {
        // exp(I theta) = cos(theta) + I sin(theta); its square is the double angle
        value_t const theta = 0.7;
        auto const z = exp(I_2d, theta);
        CHECK(z == mvec2d_e(std::cos(theta), std::sin(theta)));
        CHECK(z * z == exp(I_2d, 2 * theta));
        CHECK(nrm(z) == doctest::Approx(1.0));
    }

    TEST_CASE("ega2d calculus: time derivative and gradient of a travelling wave")
    {
        // u(q, t) = cos(k.q - w t), |k| = w: du/dt = w sin(ph), grad u = -sin(ph) k,
        // d2u/dt2 = -w^2 cos(ph)
        value_t const w = 1.3;
        vec2d const k(0.6 * w, 0.8 * w);
        auto wave = [&](vec2d const& q, value_t t) {
            return scalar2d(std::cos(value_t(dot(k, q)) - w * t));
        };
        vec2d const r0(0.37, -0.62);
        value_t const t0 = 0.29;
        value_t const ph = value_t(dot(k, r0)) - w * t0;

        CHECK(value_t(gr0(d_dt(wave, r0, t0, 1.0e-5))) ==
              doctest::Approx(w * std::sin(ph)));
        CHECK(value_t(gr0(d_dt(wave, r0, t0, central_scheme(2, 2), 1.0e-3))) ==
              doctest::Approx(-w * w * std::cos(ph)).epsilon(1.0e-5));
        auto const g = nabla(wave, r0, t0, 1.0e-5);
        CHECK(nrm(mvec2d(g) - mvec2d(-std::sin(ph) * k)) < 1.0e-8);
    }

    TEST_CASE("ega: exp, log and sqrt at the identity")
    {
        // a zero angle or a zero plane is no rotation; log of no rotation is zero; and
        // the square root of the 2 pi rotor (-1, the identity rotation) is the identity
        mvec3d_e const id3(scalar3d(1.0), bivec3d(0.0, 0.0, 0.0));
        CHECK(exp(e12_3d, 0.0) == id3);
        CHECK(exp(bivec3d(0.0, 0.0, 0.0), 1.2) == id3);
        CHECK(log(id3) == bivec3d(0.0, 0.0, 0.0));
        CHECK(sqrt(mvec3d_e(-1.0, 0.0, 0.0, 0.0)) == id3);
        CHECK(sqrt(mvec2d_e(-1.0, 0.0)) == mvec2d_e(1.0, 0.0));
    }
}
