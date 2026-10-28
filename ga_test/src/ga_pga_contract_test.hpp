// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_arith_check.hpp"
#include "ga_test_equivariance.hpp"

#include <array>
#include <bit>
#include <type_traits>
#include <vector>

/////////////////////////////////////////////////////////////////////////////////////////
// PGA contracts: the arithmetic of every value type (dual numbers included), the scalar
// type's operators, and the operations of ga_pga{2,3}dp_ops*.hpp that the topical tests
// did not reach -- each pinned against an independently known answer or a documented
// identity
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA contracts")
{
    TEST_CASE("pga2dp / pga3dp: arithmetic contract of every value type")
    {
        arith_check::check<scalar2dp>("scalar2dp");
        arith_check::check<vec2dp>("vec2dp");
        arith_check::check<bivec2dp>("bivec2dp");
        arith_check::check<pscalar2dp>("pscalar2dp");
        arith_check::check<mvec2dp_e>("mvec2dp_e");
        arith_check::check<mvec2dp_u>("mvec2dp_u");
        arith_check::check<mvec2dp>("mvec2dp");
        arith_check::check<dualnum2dp>("dualnum2dp");

        arith_check::check<scalar3dp>("scalar3dp");
        arith_check::check<vec3dp>("vec3dp");
        arith_check::check<bivec3dp>("bivec3dp");
        arith_check::check<trivec3dp>("trivec3dp");
        arith_check::check<pscalar3dp>("pscalar3dp");
        arith_check::check<mvec3dp_e>("mvec3dp_e");
        arith_check::check<mvec3dp_u>("mvec3dp_u");
        arith_check::check<mvec3dp>("mvec3dp");
        arith_check::check<dualnum3dp>("dualnum3dp");
    }

    TEST_CASE("Scalar_t: comparisons and arithmetic with raw numbers")
    {
        scalar2dp const s(2.5);

        // a Scalar_t orders against a raw number from either side
        CHECK(s < 3.0);
        CHECK(s <= 2.5);
        CHECK(s > 2.0);
        CHECK(s >= 2.5);
        CHECK_FALSE(s < 2.5);
        CHECK(2.0 < s);
        CHECK(2.5 <= s);
        CHECK(3.0 > s);
        CHECK(2.5 >= s);
        CHECK(3.0 != s);
        CHECK_FALSE(2.5 != s);

        // ... and between two Scalar_t
        CHECK(scalar2dp(1.0) < s);
        CHECK(scalar2dp(2.5) <= s);
        CHECK(scalar2dp(3.0) > s);
        CHECK(scalar2dp(2.5) >= s);

        // raw numbers add and subtract from either side, keeping the scalar's type
        CHECK(s + 1.0 == scalar2dp(3.5));
        CHECK(1.0 + s == scalar2dp(3.5));
        CHECK(s - 1.0 == scalar2dp(1.5));
        CHECK(1.0 - s == scalar2dp(-1.5));

        // abs and sqrt accept the raw type and the Scalar_t alike
        CHECK(hd::ga::abs(-2.5) == 2.5);
        CHECK(hd::ga::abs(scalar2dp(-2.5)) == 2.5);
        CHECK(hd::ga::sqrt(6.25) == 2.5);
        CHECK(hd::ga::sqrt(scalar2dp(6.25)) == 2.5);
    }

    TEST_CASE("pga2dp / pga3dp: a dual number sums grade by grade")
    {
        // every +/- involving a dual number is the sum of the operands' grade parts --
        // pinned by lifting both sides into the full multivector (a dual number lifts
        // through its two grades)
        auto const lift2 = [](dualnum2dp const& d) {
            return mvec2dp(scalar2dp(gr0(d))) + mvec2dp(pscalar2dp(gr3(d)));
        };
        scalar2dp const s2(1.5);
        pscalar2dp const p2(-2.0);
        dualnum2dp const d2(0.5, 3.0);
        vec2dp const v2(1.0, -1.0, 2.0);
        bivec2dp const b2(0.25, 0.5, -0.75);
        mvec2dp_e const e2(1.0, 2.0, 3.0, 4.0);
        mvec2dp_u const u2(-1.0, 0.5, 2.0, -3.0);
        mvec2dp const m2(1, 2, 3, 4, 5, 6, 7, 8);

        CHECK(s2 + p2 == dualnum2dp(1.5, -2.0));
        CHECK(p2 + s2 == dualnum2dp(1.5, -2.0));
        CHECK(s2 - p2 == dualnum2dp(1.5, 2.0));
        CHECK(p2 - s2 == dualnum2dp(-1.5, -2.0));
        CHECK(s2 + d2 == dualnum2dp(2.0, 3.0));
        CHECK(d2 + s2 == dualnum2dp(2.0, 3.0));
        CHECK(s2 - d2 == dualnum2dp(1.0, -3.0));
        CHECK(d2 - s2 == dualnum2dp(-1.0, 3.0));
        CHECK(p2 + d2 == dualnum2dp(0.5, 1.0));
        CHECK(d2 + p2 == dualnum2dp(0.5, 1.0));
        CHECK(p2 - d2 == dualnum2dp(-0.5, -5.0));
        CHECK(d2 - p2 == dualnum2dp(0.5, 5.0));

        CHECK(d2 + v2 == lift2(d2) + mvec2dp(v2));
        CHECK(v2 + d2 == lift2(d2) + mvec2dp(v2));
        CHECK(d2 - v2 == lift2(d2) - mvec2dp(v2));
        CHECK(v2 - d2 == mvec2dp(v2) - lift2(d2));
        CHECK(d2 + b2 == lift2(d2) + mvec2dp(b2));
        CHECK(b2 + d2 == lift2(d2) + mvec2dp(b2));
        CHECK(d2 - b2 == lift2(d2) - mvec2dp(b2));
        CHECK(b2 - d2 == mvec2dp(b2) - lift2(d2));
        CHECK(d2 + e2 == lift2(d2) + mvec2dp(e2));
        CHECK(e2 + d2 == lift2(d2) + mvec2dp(e2));
        CHECK(d2 - e2 == lift2(d2) - mvec2dp(e2));
        CHECK(e2 - d2 == mvec2dp(e2) - lift2(d2));
        CHECK(d2 + u2 == lift2(d2) + mvec2dp(u2));
        CHECK(u2 + d2 == lift2(d2) + mvec2dp(u2));
        CHECK(d2 - u2 == lift2(d2) - mvec2dp(u2));
        CHECK(u2 - d2 == mvec2dp(u2) - lift2(d2));
        CHECK(d2 + m2 == lift2(d2) + m2);
        CHECK(m2 + d2 == lift2(d2) + m2);
        CHECK(d2 - m2 == lift2(d2) - m2);
        CHECK(m2 - d2 == m2 - lift2(d2));

        auto const lift3 = [](dualnum3dp const& d) {
            return mvec3dp(scalar3dp(gr0(d))) + mvec3dp(pscalar3dp(gr4(d)));
        };
        scalar3dp const s3(1.5);
        pscalar3dp const p3(-2.0);
        dualnum3dp const d3(0.5, 3.0);
        vec3dp const v3(1.0, -1.0, 2.0, 0.5);
        bivec3dp const b3(0.25, 0.5, -0.75, 1.0, 2.0, -1.0);
        trivec3dp const t3(2.0, 0.0, -1.0, 3.0);
        mvec3dp_e const e3(1, 2, 3, 4, 5, 6, 7, 8);
        mvec3dp const m3(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

        CHECK(s3 + p3 == dualnum3dp(1.5, -2.0));
        CHECK(p3 + s3 == dualnum3dp(1.5, -2.0));
        CHECK(s3 - p3 == dualnum3dp(1.5, 2.0));
        CHECK(p3 - s3 == dualnum3dp(-1.5, -2.0));
        CHECK(s3 + d3 == dualnum3dp(2.0, 3.0));
        CHECK(d3 + s3 == dualnum3dp(2.0, 3.0));
        CHECK(s3 - d3 == dualnum3dp(1.0, -3.0));
        CHECK(d3 - s3 == dualnum3dp(-1.0, 3.0));
        CHECK(p3 + d3 == dualnum3dp(0.5, 1.0));
        CHECK(d3 + p3 == dualnum3dp(0.5, 1.0));
        CHECK(p3 - d3 == dualnum3dp(-0.5, -5.0));
        CHECK(d3 - p3 == dualnum3dp(0.5, 5.0));

        // in 3dp a dual number is even: with a bivector or an even multivector it stays
        // even, with an odd part it needs the full multivector
        CHECK(mvec3dp(d3 + b3) == lift3(d3) + mvec3dp(b3));
        CHECK(mvec3dp(b3 + d3) == lift3(d3) + mvec3dp(b3));
        CHECK(mvec3dp(d3 - b3) == lift3(d3) - mvec3dp(b3));
        CHECK(mvec3dp(b3 - d3) == mvec3dp(b3) - lift3(d3));
        CHECK(mvec3dp(d3 + e3) == lift3(d3) + mvec3dp(e3));
        CHECK(mvec3dp(e3 + d3) == lift3(d3) + mvec3dp(e3));
        CHECK(mvec3dp(d3 - e3) == lift3(d3) - mvec3dp(e3));
        CHECK(mvec3dp(e3 - d3) == mvec3dp(e3) - lift3(d3));
        CHECK(d3 + v3 == lift3(d3) + mvec3dp(v3));
        CHECK(v3 + d3 == lift3(d3) + mvec3dp(v3));
        CHECK(d3 - v3 == lift3(d3) - mvec3dp(v3));
        CHECK(v3 - d3 == mvec3dp(v3) - lift3(d3));
        CHECK(d3 + t3 == lift3(d3) + mvec3dp(t3));
        CHECK(t3 + d3 == lift3(d3) + mvec3dp(t3));
        CHECK(d3 - t3 == lift3(d3) - mvec3dp(t3));
        CHECK(t3 - d3 == mvec3dp(t3) - lift3(d3));
        mvec3dp_u const u3(1, -2, 3, -4, 5, -6, 7, -8);
        CHECK(d3 + u3 == lift3(d3) + mvec3dp(u3));
        CHECK(u3 + d3 == lift3(d3) + mvec3dp(u3));
        CHECK(d3 - u3 == lift3(d3) - mvec3dp(u3));
        CHECK(u3 - d3 == mvec3dp(u3) - lift3(d3));
        CHECK(d3 + m3 == lift3(d3) + m3);
        CHECK(m3 + d3 == lift3(d3) + m3);
        CHECK(d3 - m3 == lift3(d3) - m3);
        CHECK(m3 - d3 == m3 - lift3(d3));
    }

    TEST_CASE("pga2dp / pga3dp: a dual number is built from (scalar, pseudoscalar) only")
    {
        // the right parts in the right order construct it ...
        static_assert(std::is_constructible_v<dualnum2dp, scalar2dp, pscalar2dp>);
        static_assert(std::is_constructible_v<dualnum3dp, scalar3dp, pscalar3dp>);
        CHECK(dualnum2dp(scalar2dp(1.0), pscalar2dp(2.0)) == dualnum2dp(1.0, 2.0));
        CHECK(dualnum3dp(scalar3dp(1.0), pscalar3dp(2.0)) == dualnum3dp(1.0, 2.0));
        // ... while swapped parts, or parts of another algebra, do not compile (they
        // used to, silently: the swapped call put the pseudoscalar into the scalar slot)
        static_assert(!std::is_constructible_v<dualnum2dp, pscalar2dp, scalar2dp>);
        static_assert(!std::is_constructible_v<dualnum3dp, pscalar3dp, scalar3dp>);
        static_assert(!std::is_constructible_v<dualnum2dp, scalar3dp, pscalar3dp>);
        static_assert(!std::is_constructible_v<dualnum3dp, scalar2dp, pscalar2dp>);
    }

    TEST_CASE("pga2dp / pga3dp: the angle between two vectors lies in [0, pi]")
    {
        // an acos: never negative, the same in both argument orders, pi when antiparallel
        CHECK(angle(vec2dp(1.0, 0.0, 0.0), vec2dp(-1.0, 0.0, 0.0)) ==
              doctest::Approx(pi));
        CHECK(angle(vec2dp(0.0, 1.0, 0.0), vec2dp(1.0, 0.0, 0.0)) ==
              doctest::Approx(pi / 2));
        CHECK(angle(vec2dp(1.0, 0.0, 0.0), vec2dp(0.0, 1.0, 0.0)) ==
              doctest::Approx(pi / 2));
        CHECK(angle(vec3dp(1.0, 0.0, 0.0, 0.0), vec3dp(-1.0, 0.0, 0.0, 0.0)) ==
              doctest::Approx(pi));
        CHECK(angle(vec3dp(0.0, 1.0, 0.0, 0.0), vec3dp(1.0, 0.0, 0.0, 0.0)) ==
              doctest::Approx(pi / 2));
    }

    TEST_CASE("pga3dp: the angle between a line and a plane is signed, in [-pi/2, pi/2]")
    {
        // 0 in the plane, +pi/2 along its normal, -pi/2 against it; the same in both
        // argument orders; reversing the line or the plane flips the sign
        value_t const a = deg2rad(30.0);
        trivec3dp const pl(0.0, 0.0, 1.0, -3.0); // z == 3, normal +z
        vec3dp const u(std::cos(a), 0.0, std::sin(a), 0.0);
        bivec3dp const L = wdg(vec3dp(1.0, 2.0, 5.0, 1.0), u); // above the plane
        CHECK(angle(L, pl) == doctest::Approx(a));
        CHECK(angle(pl, L) == doctest::Approx(a));
        CHECK(angle(-L, pl) == doctest::Approx(-a));
        CHECK(angle(pl, -L) == doctest::Approx(-a));
        CHECK(angle(L, -pl) == doctest::Approx(-a));
        CHECK(angle(wdg(O_3dp, e3_3dp), pl) == doctest::Approx(pi / 2));
        CHECK(angle(wdg(O_3dp, -e3_3dp), pl) == doctest::Approx(-pi / 2));
        CHECK(angle(wdg(O_3dp, e1_3dp), pl) == doctest::Approx(0.0));

        // only directions enter: the same direction below the plane gives the same
        // angle, and the weights of line and plane divide out
        CHECK(angle(wdg(vec3dp(-4.0, 0.5, -2.0, 1.0), u), pl) == doctest::Approx(a));
        CHECK(angle(7.0 * L, 0.5 * pl) == doctest::Approx(a));

        // pi/2 minus the angle between the line's direction and the plane's normal
        bivec3dp const L2 =
            wdg(vec3dp(0.3, -1.2, 0.7, 1.0), vec3dp(-0.4, 1.3, -0.9, 0.0));
        trivec3dp const pl2(0.6, -0.2, 0.8, 1.5);
        // (both lifted back as ideal vectors, i.e. directions)
        vec3dp const n2(to_vec3d(pl2), 0.0);
        CHECK(angle(L2, pl2) ==
              doctest::Approx(pi / 2 - angle(vec3dp(to_vec3d(L2), 0.0), n2)));
        CHECK(angle(-L2, pl2) ==
              doctest::Approx(pi / 2 - angle(vec3dp(to_vec3d(-L2), 0.0), n2)));
    }

    TEST_CASE("pga2dp / pga3dp: the regressive grade complements the grade")
    {
        // rgr(x) = n - gr(x): the grade of the complement
        CHECK(rgr(scalar2dp{}) == 3);
        CHECK(rgr(vec2dp{}) == 2);
        CHECK(rgr(bivec2dp{}) == 1);
        CHECK(rgr(pscalar2dp{}) == 0);
        CHECK(rgr(scalar3dp{}) == 4);
        CHECK(rgr(vec3dp{}) == 3);
        CHECK(rgr(bivec3dp{}) == 2);
        CHECK(rgr(trivec3dp{}) == 1);
        CHECK(rgr(pscalar3dp{}) == 0);
        CHECK(gr(bivec3dp{}) + rgr(bivec3dp{}) == 4);
    }

    TEST_CASE("pga2dp / pga3dp: the exported metric arrays are the live dot and rdot")
    {
        // G[i, j] == dot(e_i, e_j) and Gr[i, j] == rdot(e_i, e_j) over every pair of
        // basis blades (in the order the arrays document), each blade the full
        // multivector with a single unit coefficient
        auto check = [](auto zero, auto G, auto Gr) {
            using M = decltype(zero);
            constexpr size_t n = sizeof(M) / sizeof(value_t);
            auto blade = [](size_t i) {
                std::array<value_t, n> c{};
                c[i] = 1.0;
                return std::bit_cast<M>(c);
            };
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    CHECK(value_t(dot(blade(i), blade(j))) == G[i, j]);
                    CHECK(value_t(rdot(blade(i), blade(j))) == Gr[i, j]);
                }
            }
        };
        check(mvec2dp{}, pga2dp_metric_view(), pga2dp_rmetric_view());
        check(mvec3dp{}, pga3dp_metric_view(), pga3dp_rmetric_view());
    }

    TEST_CASE("pga2dp: a Euclidean vector lifts with its third component")
    {
        // the EGA -> PGA lift: vec2dp(v, 1) is the point at v, vec2dp(v, 0) the direction
        CHECK(vec2dp(vec2d(1.5, -2.0), 1.0) == vec2dp(1.5, -2.0, 1.0));
        CHECK(vec2dp(vec2d(1.5, -2.0), 0.0) == vec2dp(1.5, -2.0, 0.0));
    }

    TEST_CASE("pga3dp: the geometric objects' constructors")
    {
        // a Euclidean point lifts with w = 1, from a Point3d or straight from a vec3d
        CHECK(point3d(vec3d(1.0, 2.0, 3.0)) == point3d(1.0, 2.0, 3.0));
        CHECK(point3dp(point3d(1.0, 2.0, 3.0)) == vec3dp(1.0, 2.0, 3.0, 1.0));
        CHECK(point3dp(vec3d(1.0, 2.0, 3.0)) == vec3dp(1.0, 2.0, 3.0, 1.0));

        // a line from its direction and moment is the bivector of those parts, and it is
        // a line: the moment is perpendicular to the direction
        vec3d const dir(0.0, 0.0, 1.0);
        bivec3d const mom = hd::ga::ega::wdg(vec3d(1.0, 2.0, 0.0), dir); // via (1, 2, *)
        line3d const L(dir, mom);
        CHECK(L == bivec3dp(dir, mom));
        CHECK(is_congruent(L,
                           join(vec3dp(1.0, 2.0, 0.0, 1.0), vec3dp(1.0, 2.0, 5.0, 1.0))));

        // a plane from a line and a point off it is their join
        plane3d const P(L, point3d(4.0, 2.0, 0.0));
        CHECK(P == join(L, vec3dp(4.0, 2.0, 0.0, 1.0)));
        CHECK(is_congruent(P, trivec3dp(0.0, 1.0, 0.0, -2.0))); // the plane y = 2
    }

    TEST_CASE("pga2dp / pga3dp: att is the meet with the horizon, for every overload")
    {
        // documented: att(u) = rwdg(u, H) with H the horizon
        vec2dp const v2(1.0, 2.0, 3.0);
        CHECK(att(v2) == rwdg(v2, H_2dp));
        CHECK(att(point2dp(v2)) == rwdg(v2, H_2dp));
        CHECK(att(bivec2dp(1.0, 2.0, 3.0)) == rwdg(bivec2dp(1.0, 2.0, 3.0), H_2dp));
        CHECK(att(line2d(1.0, 2.0, 3.0)) == rwdg(bivec2dp(1.0, 2.0, 3.0), H_2dp));
        CHECK(att(pscalar2dp(2.0)) == rwdg(pscalar2dp(2.0), H_2dp));

        vec3dp const v3(1.0, 2.0, 3.0, 4.0);
        CHECK(att(v3) == rwdg(v3, H_3dp));
        CHECK(att(point3dp(v3)) == rwdg(v3, H_3dp));
        bivec3dp const b3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        CHECK(att(b3) == rwdg(b3, H_3dp));
        trivec3dp const t3(1.0, 2.0, 3.0, 4.0);
        CHECK(att(t3) == rwdg(t3, H_3dp));
        CHECK(att(plane3d(t3)) == rwdg(t3, H_3dp));
        CHECK(att(pscalar3dp(2.0)) == rwdg(pscalar3dp(2.0), H_3dp));
    }

    TEST_CASE("pga2dp / pga3dp: a translation motor from a Vector moves a point by it")
    {
        auto const M2 = get_motor(vector2d(vec2d(1.5, -2.0)));
        CHECK(move2dp(vec2dp(1.0, 2.0, 1.0), M2) == vec2dp(2.5, 0.0, 1.0));
        // a direction (ideal point) is not translated
        CHECK(move2dp(vec2dp(1.0, 2.0, 0.0), M2) == vec2dp(1.0, 2.0, 0.0));

        auto const M3 = get_motor(vector3d(vec3d(1.5, -2.0, 0.5)));
        CHECK(move3dp(vec3dp(1.0, 2.0, 3.0, 1.0), M3) == vec3dp(2.5, 0.0, 3.5, 1.0));
        CHECK(move3dp(vec3dp(1.0, 2.0, 3.0, 0.0), M3) == vec3dp(1.0, 2.0, 3.0, 0.0));
    }

    TEST_CASE("pga2dp: moving a motor is the change of frame")
    {
        // M acts on the motor M0 so that the moved motor does in the moved frame what M0
        // did in the old one:  move(move(P, M0), M) == move(move(P, M), move(M0, M))
        auto const M0 = get_motor(vec2dp(1.0, 0.5, 1.0), deg2rad(40.0));
        auto const M = rgpr(get_motor(vec2dp(-1.0, 2.0, 1.0), deg2rad(-75.0)),
                            get_motor(vector2d(vec2d(0.3, 0.7))));
        for (vec2dp const P :
             {vec2dp(0.0, 0.0, 1.0), vec2dp(2.0, -1.0, 1.0), vec2dp(-0.5, 3.0, 1.0)}) {
            CHECK(move2dp(move2dp(P, M0), M) == move2dp(move2dp(P, M), move2dp(M0, M)));
        }
    }

    TEST_CASE("pga2dp / pga3dp: move*_opt and the batch overloads agree with move*")
    {
        auto const M2 = rgpr(get_motor(vec2dp(1.0, 0.5, 1.0), deg2rad(40.0)),
                             get_motor(vector2d(vec2d(0.3, 0.7))));
        std::vector<bivec2dp> const B2{bivec2dp(1.0, 2.0, 3.0), bivec2dp(-1.0, 0.5, 0.0),
                                       e12_2dp};
        auto const moved2 = move2dp(B2, M2);
        REQUIRE(moved2.size() == B2.size());
        for (size_t i = 0; i < B2.size(); ++i) {
            CHECK(move2dp_opt(B2[i], M2) == move2dp(B2[i], M2));
            CHECK(moved2[i] == move2dp(B2[i], M2));
        }

        auto const M3 =
            rgpr(get_motor(join(vec3dp(1.0, 0.0, 0.0, 1.0), vec3dp(1.0, 1.0, 1.0, 1.0)),
                           deg2rad(35.0)),
                 get_motor(vector3d(vec3d(0.3, -0.7, 0.2))));
        std::vector<bivec3dp> const B3{bivec3dp(1, 2, 3, 4, 5, 6),
                                       bivec3dp(0, 0, 1, 1, 0, 0)};
        std::vector<trivec3dp> const T3{trivec3dp(1, 2, 3, 4), trivec3dp(0, 0, 1, -3)};
        auto const movedB = move3dp(B3, M3);
        auto const movedT = move3dp(T3, M3);
        REQUIRE(movedB.size() == B3.size());
        REQUIRE(movedT.size() == T3.size());
        for (size_t i = 0; i < B3.size(); ++i) {
            CHECK(move3dp_opt(B3[i], M3) == move3dp(B3[i], M3));
            CHECK(movedB[i] == move3dp(B3[i], M3));
        }
        for (size_t i = 0; i < T3.size(); ++i) {
            CHECK(move3dp_opt(T3[i], M3) == move3dp(T3[i], M3));
            CHECK(movedT[i] == move3dp(T3[i], M3));
        }

        // a moved plane keeps the moved point: incidence survives the motion
        vec3dp const on_plane(1.0, 2.0, 3.0, 1.0); // on z = 3
        trivec3dp const plane(0.0, 0.0, 1.0, -3.0);
        REQUIRE(value_t(wdg(on_plane, plane)) == doctest::Approx(0.0));
        CHECK(value_t(wdg(move3dp(on_plane, M3), move3dp(plane, M3))) ==
              doctest::Approx(0.0).epsilon(1.0e-12));

        // a motion of a full multivector moves each grade (a unit motor leaves the scalar
        // and the pseudoscalar parts unchanged)
        mvec3dp const X(mvec3dp(scalar3dp(2.0)) + mvec3dp(B3[0]) + mvec3dp(T3[0]) +
                        mvec3dp(vec3dp(1.0, -1.0, 0.5, 1.0)) + mvec3dp(pscalar3dp(3.0)));
        CHECK(move3dp(X, M3) == mvec3dp(scalar3dp(2.0)) + mvec3dp(move3dp(B3[0], M3)) +
                                    mvec3dp(move3dp(T3[0], M3)) +
                                    mvec3dp(move3dp(vec3dp(1.0, -1.0, 0.5, 1.0), M3)) +
                                    mvec3dp(pscalar3dp(3.0)));
    }

    TEST_CASE("pga3dp: expand builds the perpendicular object through the given one")
    {
        // the line through a point, perpendicular to a plane
        auto const L =
            expand(point3d(1.0, 2.0, 3.0), plane3d(trivec3dp(0.0, 0.0, 1.0, 0.0)));
        CHECK(is_congruent(L,
                           join(vec3dp(1.0, 2.0, 3.0, 1.0), vec3dp(1.0, 2.0, 4.0, 1.0))));

        // the plane through a point, perpendicular to a line: the z axis -> z = 3
        auto const z_axis = line3d(join(O_3dp, vec3dp(0.0, 0.0, 1.0, 1.0)));
        CHECK(is_congruent(expand(point3d(1.0, 2.0, 3.0), z_axis),
                           trivec3dp(0.0, 0.0, 1.0, -3.0)));

        // the plane through a line, perpendicular to a plane: the line along x at z = 1,
        // against z = 0 -> the plane y = 0
        auto const x_line =
            line3d(join(vec3dp(0.0, 0.0, 1.0, 1.0), vec3dp(1.0, 0.0, 1.0, 1.0)));
        CHECK(is_congruent(expand(x_line, plane3d(trivec3dp(0.0, 0.0, 1.0, 0.0))),
                           trivec3dp(0.0, 1.0, 0.0, 0.0)));
    }

    TEST_CASE("pga2dp / pga3dp: is_close on scalars, pseudoscalars and motors")
    {
        CHECK(is_close(scalar2dp(3.0), scalar2dp(3.0 + 3.0e-14)));
        CHECK_FALSE(is_close(scalar2dp(3.0), scalar2dp(3.1)));
        CHECK(is_close(pscalar2dp(-2.0e6), pscalar2dp(-2.0e6 * (1.0 + 1.0e-14))));
        CHECK_FALSE(is_close(pscalar2dp(1.0), pscalar2dp(-1.0)));
        auto const M2 = get_motor(vec2dp(1.0, 0.5, 1.0), deg2rad(40.0));
        CHECK(is_close(M2, M2 * (1.0 + 1.0e-14)));
        CHECK_FALSE(is_close(M2, -M2)); // the same motion, not the same value

        CHECK(is_close(scalar3dp(3.0), scalar3dp(3.0 + 3.0e-14)));
        CHECK_FALSE(is_close(scalar3dp(3.0), scalar3dp(3.1)));
        CHECK(is_close(pscalar3dp(0.5), pscalar3dp(0.5 + 1.0e-15)));
        CHECK_FALSE(is_close(pscalar3dp(0.5), pscalar3dp(0.6)));
        auto const M3 = get_motor(join(O_3dp, vec3dp(0.0, 0.0, 1.0, 1.0)), deg2rad(40.0));
        CHECK(is_close(M3, M3 * (1.0 + 1.0e-14)));
        CHECK_FALSE(is_close(M3, -M3));
    }

    TEST_CASE("pga2dp / pga3dp: unitize divides a single-coefficient weight out signed")
    {
        // where the weight is ONE coefficient (a point, a pseudoscalar, a dual number)
        // unitize() divides by it signed, so the result has weight +1 whatever the sign
        // it came with. A multivector's weight has several coefficients in general, so
        // its unitize() divides by the unsigned weight norm and keeps a negative weight's
        // sign: the same point unitizes differently held as a multivector.
        vec2dp const p2(2.0, 4.0, -2.0);
        CHECK(unitize(p2) == vec2dp(-1.0, -2.0, 1.0));
        CHECK(unitize(point2dp(p2)) == vec2dp(-1.0, -2.0, 1.0));
        CHECK(unitize(mvec2dp(p2)) == mvec2dp(vec2dp(1.0, 2.0, -1.0)));
        // with a positive weight all three agree
        CHECK(mvec2dp(unitize(-p2)) == unitize(mvec2dp(-p2)));

        vec3dp const p3(2.0, 4.0, 6.0, -2.0);
        CHECK(unitize(p3) == vec3dp(-1.0, -2.0, -3.0, 1.0));
        CHECK(unitize(point3dp(p3)) == vec3dp(-1.0, -2.0, -3.0, 1.0));
        CHECK(unitize(mvec3dp(p3)) == mvec3dp(vec3dp(1.0, 2.0, 3.0, -1.0)));

        // a pseudoscalar's weight is the pseudoscalar itself
        CHECK(unitize(pscalar2dp(-3.0)) == pscalar2dp(1.0));
        CHECK(unitize(pscalar2dp(0.25)) == pscalar2dp(1.0));
        CHECK(unitize(pscalar3dp(-3.0)) == pscalar3dp(1.0));
        CHECK(unitize(pscalar3dp(0.25)) == pscalar3dp(1.0));
        CHECK(unitize(mvec2dp(pscalar2dp(-3.0))) == mvec2dp(pscalar2dp(-1.0)));

        // a dual number's weight is its pseudoscalar part, and it is divided out signed
        CHECK(unitize(dualnum2dp(2.0, 4.0)) == dualnum2dp(0.5, 1.0));
        CHECK(unitize(dualnum2dp(1.0, -2.0)) == dualnum2dp(-0.5, 1.0));
        CHECK(unitize(dualnum3dp(3.0, -1.5)) == dualnum3dp(-2.0, 1.0));
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        CHECK_THROWS_AS((void)unitize(dualnum2dp(1.0, 0.0)), std::runtime_error);
        CHECK_THROWS_AS((void)unitize(dualnum3dp(1.0, 0.0)), std::runtime_error);
#endif
    }

    TEST_CASE("pga2dp / pga3dp: try_unitize unitizes what has weight, reports the rest")
    {
        bool unitized = false;
        CHECK(try_unitize(vec2dp(2.0, 4.0, 2.0), &unitized) == vec2dp(1.0, 2.0, 1.0));
        CHECK(unitized);
        // an ideal point has no weight: returned unchanged and reported
        CHECK(try_unitize(vec2dp(2.0, 4.0, 0.0), &unitized) == vec2dp(2.0, 4.0, 0.0));
        CHECK_FALSE(unitized);
        CHECK(try_unitize(point2dp(vec2dp(2.0, 4.0, 2.0)), &unitized) ==
              vec2dp(1.0, 2.0, 1.0));
        CHECK(unitized);
        // a weight of -1 is not yet +1: it is unitized too (the test is on the signed
        // weight, not on its square)
        CHECK(try_unitize(vec2dp(2.0, 4.0, -1.0), &unitized) == vec2dp(-2.0, -4.0, 1.0));
        CHECK(unitized);
        CHECK(try_unitize(point2dp(vec2dp(2.0, 4.0, -1.0))) == vec2dp(-2.0, -4.0, 1.0));
        CHECK(try_unitize(vec3dp(1.0, 2.0, 3.0, -1.0), &unitized) ==
              vec3dp(-1.0, -2.0, -3.0, 1.0));
        CHECK(unitized);
        CHECK(try_unitize(point3dp(vec3dp(1.0, 2.0, 3.0, -1.0))) ==
              vec3dp(-1.0, -2.0, -3.0, 1.0));

        bivec2dp const B(3.0, 4.0, 10.0); // weight norm 5
        CHECK(try_unitize(B, &unitized) == unitize(B));
        CHECK(unitized);
        CHECK(try_unitize(line2d(B), &unitized) == unitize(B));
        CHECK(unitized);
        CHECK(try_unitize(bivec2dp(0.0, 0.0, 3.0), &unitized) == bivec2dp(0.0, 0.0, 3.0));
        CHECK_FALSE(unitized);
        // the out-parameter is optional
        CHECK(try_unitize(B) == unitize(B));

        trivec3dp const t(0.0, 3.0, 4.0, 10.0);
        CHECK(try_unitize(t, &unitized) == unitize(t));
        CHECK(unitized);
        CHECK(try_unitize(trivec3dp(0.0, 0.0, 0.0, 2.0), &unitized) ==
              trivec3dp(0.0, 0.0, 0.0, 2.0));
        CHECK_FALSE(unitized);
    }

    TEST_CASE(
        "pga2dp / pga3dp: try_unitize leaves an ideal point, line or plane as it is")
    {
        bool unitized = true;
        CHECK(try_unitize(point2dp(vec2dp(2.0, 4.0, 0.0)), &unitized) ==
              vec2dp(2.0, 4.0, 0.0));
        CHECK_FALSE(unitized);
        CHECK(try_unitize(line2d(0.0, 0.0, 3.0), &unitized) == bivec2dp(0.0, 0.0, 3.0));
        CHECK_FALSE(unitized);
        CHECK(try_unitize(point3dp(vec3dp(1.0, 2.0, 3.0, 0.0)), &unitized) ==
              vec3dp(1.0, 2.0, 3.0, 0.0));
        CHECK_FALSE(unitized);
        CHECK(try_unitize(line3d(0.0, 0.0, 0.0, 1.0, 2.0, 3.0), &unitized) ==
              bivec3dp(0.0, 0.0, 0.0, 1.0, 2.0, 3.0));
        CHECK_FALSE(unitized);
        CHECK(try_unitize(plane3d(trivec3dp(0.0, 0.0, 0.0, 2.0)), &unitized) ==
              trivec3dp(0.0, 0.0, 0.0, 2.0));
        CHECK_FALSE(unitized);
    }

    TEST_CASE("pga3dp: get_motor_from_lines carries l1 onto l2, whatever their scale")
    {
        // documented: the motor that moves line l1 to line l2. It used to unitize copies
        // of its arguments and then use the originals, so only unit lines worked.
        bivec3dp const l1 =
            unitize(join(vec3dp(0.0, 0.0, 0.0, 1.0), vec3dp(1.0, 0.0, 0.0, 1.0)));
        bivec3dp const l2 =
            unitize(join(vec3dp(0.0, 1.0, 1.0, 1.0), vec3dp(0.0, 1.0, 2.0, 1.0)));
        auto const M = get_motor_from_lines(l1, l2);
        CHECK(is_close(unitize(move3dp(l1, M)), l2));
        for (value_t const k1 : {2.0, 0.25})
            for (value_t const k2 : {3.0, 0.5})
                CHECK(is_same_motion(get_motor_from_lines(k1 * l1, k2 * l2), M));
    }

    TEST_CASE("pga2dp: get_motor_from_lines is the product of the two reflections")
    {
        // documented: reflect across B1, then across B2 -- a rotation about their
        // intersection by TWICE the angle between them. The x axis and the 45 deg line
        // through (0, 1) meet at (-1, 0): a point on the x axis turns by 90 deg there.
        bivec2dp const m1 = unitize(join(vec2dp(0.0, 0.0, 1.0), vec2dp(1.0, 0.0, 1.0)));
        bivec2dp const m2 = unitize(join(vec2dp(0.0, 1.0, 1.0), vec2dp(1.0, 2.0, 1.0)));
        auto const M = get_motor_from_lines(m1, m2);
        CHECK(unitize(move2dp(vec2dp(0.5, 0.0, 1.0), M)) == vec2dp(-1.0, 1.5, 1.0));
        CHECK(unitize(move2dp(vec2dp(-1.0, 0.0, 1.0), M)) ==
              vec2dp(-1.0, 0.0, 1.0)); // fixed
        // scaled lines give the same motion
        CHECK(is_same_motion(get_motor_from_lines(2.0 * m1, 0.5 * m2), M));
    }

    TEST_CASE("pga2dp: rlog on the far sheet -- two half turns make a translation")
    {
        // Two half turns, about (0, 0) and then about (1, 0), are a translation by
        // (2, 0), and their product lies on the FAR sheet of the double cover
        // (pseudoscalar -1). Until 2026-10-27 the 2D rlog returned that translation
        // negated: rexp(rlog(M)) moved the origin to (-2, 0), and screw_axis said so too.
        auto const M = rgpr(get_motor(vec2dp{1.0, 0.0, 1.0}, pi),
                            get_motor(vec2dp{0.0, 0.0, 1.0}, pi));
        REQUIRE(M.c3 < -0.99);
        CHECK(is_close(unitize(move2dp(O_2dp, M)), vec2dp{2.0, 0.0, 1.0}));
        CHECK(is_same_motion(rexp(rlog(M)), M));
        CHECK(is_close(unitize(move2dp(O_2dp, rexp(rlog(M)))), vec2dp{2.0, 0.0, 1.0}));
        auto const sa = screw_axis(M);
        CHECK(std::abs(sa.translation.x - 2.0) < 1.0e-12);
        CHECK(std::abs(sa.translation.y) < 1.0e-12);
    }

    TEST_CASE("pga2dp / pga3dp: the motor builders unitize a SMALL weight, not just a "
              "large one")
    {
        // "Carries weight" is detail::has_weight, the predicate unitize() uses. Until
        // 2026-10-27 the motor builders tested the squared weight against eps, so an
        // input of weight below ~3e-8 -- the meet of two nearly parallel lines, a line
        // through two close points -- was used as it stood: a 90 deg turn of (2, 2)
        // about the point (1, 2) at weight 1e-9 returned (2, 2), and get_motor(L, theta)
        // refused a line unitize() accepts.
        vec2dp const P{1.0e-9, 2.0e-9, 1.0e-9};
        CHECK(is_same_motion(get_motor(P, pi / 2), get_motor(unitize(P), pi / 2)));
        CHECK(is_close(unitize(move2dp(vec2dp{2.0, 2.0, 1.0}, get_motor(P, pi / 2))),
                       vec2dp{1.0, 3.0, 1.0}));

        bivec3dp const shortx =
            wdg(vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{1.0e-8, 0.0, 0.0, 1.0});
        bivec3dp const y_axis =
            wdg(vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{0.0, 1.0, 0.0, 1.0});
        CHECK(is_same_motion(get_motor_from_lines(shortx, y_axis),
                             get_motor_from_lines(unitize(shortx), y_axis)));
        CHECK(is_same_motion(get_motor(shortx, 0.7), get_motor(unitize(shortx), 0.7)));
        CHECK(is_same_motion(get_motor(shortx, 0.7, 0.2),
                             get_motor(unitize(shortx), 0.7, 0.2)));

        // ... and an IDEAL line has no axis: refused by both forms of get_motor
        bivec3dp const ideal{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        CHECK_THROWS_AS(get_motor(ideal, 0.7), std::invalid_argument);
        CHECK_THROWS_AS(get_motor(ideal, 0.7, 0.2), std::invalid_argument);
    }

    TEST_CASE("pga2dp / pga3dp: rlog round-trips on both sheets, through both roots")
    {
        // A versor M and -M are ONE transformation (the double cover). exp(log(M)) must
        // perform M's transformation at every angle -- across pi and 2 pi, where the
        // sine of the half angle passes 0 -- and on BOTH sheets; it is checked by the
        // transformation, not the representation. A systematic check of every
        // algebra's logarithm, made 2026-10-27 after a far-sheet defect in one of them:
        // it found the far sheet wrong in the pga2dp rlog, in sta4ds for boosts and null
        // rotors, and in both cga rlogs, and a screw near a full turn 1e-7 off in pga3dp.
        bivec3dp const L =
            unitize(wdg(vec3dp{1.0, 2.0, 0.5, 1.0}, vec3dp{1.3, 1.6, 1.4, 1.0}));
        for (double const s : {1.0, -1.0}) {
            CAPTURE(s);
            for (double const a :
                 {0.3, pi - 1.0e-9, pi, pi + 1.0e-9, 2.0 * pi - 1.0e-9, 2.0 * pi,
                  2.0 * pi + 1.0e-9, 2.0 * pi + 0.3, 3.0 * pi, 4.0 * pi - 1.0e-9}) {
                CAPTURE(a);
                auto const M2 = s * get_motor(vec2dp{1.0, 2.0, 1.0}, a);
                CHECK(is_same_motion(rexp(rlog(M2)), M2));
                auto const M3 = s * get_motor(L, a);
                CHECK(is_same_motion(rexp(rlog(M3)), M3));
                auto const S3 = s * get_motor(L, a, 0.7); // a screw
                CHECK(is_same_motion(rexp(rlog(S3)), S3));
            }
            auto const T2 = s * get_motor(vec2dp{0.4, -1.1, 0.0});
            CHECK(is_same_motion(rexp(rlog(T2)), T2));
            auto const T3 = s * get_motor(vec3dp{0.4, -1.1, 0.3, 0.0});
            CHECK(is_same_motion(rexp(rlog(T3)), T3));
            // two half turns about parallel axes: a translation on the far sheet
            auto const H2 = s * rgpr(get_motor(vec2dp{1.0, 0.0, 1.0}, pi),
                                     get_motor(vec2dp{0.0, 0.0, 1.0}, pi));
            CHECK(is_same_motion(rexp(rlog(H2)), H2));
            auto const L1 =
                unitize(wdg(vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 1.0}));
            auto const L2 =
                unitize(wdg(vec3dp{1.0, 0.0, 0.0, 1.0}, vec3dp{1.0, 0.0, 1.0, 1.0}));
            auto const H3 = s * rgpr(get_motor(L2, pi), get_motor(L1, pi));
            CHECK(is_same_motion(rexp(rlog(H3)), H3));
        }
    }

    TEST_CASE("pga2dp / pga3dp: a small object inverts and is projected onto")
    {
        // A SMALL object is not a degenerate one (2026-10-27). The storage types'
        // operator/ refused any divisor below 5.5e-16, absolutely, and overruled the
        // relative guards in front of it; the projections, angles and reciprocal frames
        // carried their own absolute floors. Each answer here is the unit-size one,
        // scaled as the contract says, and a ZERO object is still refused.
        bivec3dp const B{0.3, 0.1, -0.2, 0.5, 0.7, -0.4};
        auto const Bs = bivec3dp(1.0e-4 * B);
        auto const P = Bs * inv(Bs);
        CHECK(std::abs(value_t(gr0(P)) - 1.0) < 1.0e-12);
        CHECK(is_close(inv(Bs), bivec3dp(1.0e4 * inv(B))));
        // a projection onto a plane scaled by 1e-8 is the projection onto the plane
        trivec3dp const pl =
            wdg(wdg(vec3dp{0.0, 0.0, 1.0, 1.0}, vec3dp{1.0, 0.0, 1.0, 1.0}),
                vec3dp{0.0, 1.0, 1.2, 1.0});
        vec3dp const X{0.4, -0.3, 2.0, 1.0};
        CHECK(is_close(unitize(project_onto(X, trivec3dp(1.0e-8 * pl))),
                       unitize(project_onto(X, pl))));
        // angle between small directions
        CHECK(
            std::abs(angle(vec3dp{1.0e-9, 0.0, 0.0, 0.0}, vec3dp{0.0, 1.0e-9, 0.0, 0.0}) -
                     0.5 * pi) < 1.0e-12);
        CHECK(std::abs(angle(vec2dp{1.0e-9, 0.0, 0.0}, vec2dp{0.0, 1.0e-9, 0.0}) -
                       0.5 * pi) < 1.0e-12);
        CHECK_THROWS(angle(vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{1.0, 0.0, 0.0, 0.0}));
    }
    TEST_CASE("pga2dp / pga3dp: reading a motion off a motor -- the shortest motion")
    {
        // rlog() reads the full angle in (-2 pi, 2 pi] and -M a full turn off;
        // rlog_shortest() reads the motion's shortest representative, [-pi, pi], the
        // same for M and -M
        auto yaw = [](value_t a) {
            return rexp(0.5 * bivec3dp{0.0, 0.0, a, 0.0, 0.0, 0.0});
        };
        value_t const deg = pi / 180.0;
        for (auto const [a, e] :
             {std::pair{90.0, 90.0}, std::pair{270.0, -90.0}, std::pair{359.0, -1.0},
              std::pair{361.0, 1.0}, std::pair{450.0, 90.0}, std::pair{-200.0, 160.0}}) {
            CAPTURE(a);
            auto const Y = yaw(a * deg);
            CHECK(std::abs(2.0 * rlog_shortest(Y).vz - e * deg) < 1.0e-12);
            CHECK(std::abs(2.0 * rlog_shortest(mvec3dp_e(-Y)).vz - e * deg) < 1.0e-12);
            auto const R = get_motor(vec2dp{0.0, 0.0, 1.0}, a * deg);
            CHECK(std::abs(2.0 * rlog_shortest(R).z - e * deg) < 1.0e-12);
            CHECK(std::abs(2.0 * rlog_shortest(mvec2dp_u(-R)).z - e * deg) < 1.0e-12);
        }
        // the shortest generator still generates the motion, a screw included
        bivec3dp const L =
            unitize(wdg(vec3dp{1.0, 2.0, 0.5, 1.0}, vec3dp{1.3, 1.6, 1.4, 1.0}));
        for (value_t const a : {0.4, 2.0, 3.5, 6.0, 9.0}) {
            auto const S = get_motor(L, a, 0.7);
            CHECK(is_same_motion(rexp(rlog_shortest(S)), S));
            auto const P = get_motor(vec2dp{1.0, 2.0, 1.0}, a);
            CHECK(is_same_motion(rexp(rlog_shortest(P)), P));
        }
    }

    TEST_CASE("pga2dp / pga3dp: reading a motion off a motor -- the error between two")
    {
        fmt::println("pga2dp / pga3dp: motion_error_body / motion_error_world");

        // A roll of 0.01 under a heading psi, M = Y(psi) R(0.01): the error from Y in
        // body axes is the roll, at every heading. The difference of the ABSOLUTE
        // logarithms is not (the falsification): the log of a product is not the sum
        // of the logs, and the roll reads turned by psi / 2 and scaled by
        // (psi / 2) / sin(psi / 2) -- 311 times too large at 359 degrees.
        value_t const deg = pi / 180.0, eps_roll = 0.01;
        auto yaw = [](value_t a) {
            return rexp(0.5 * bivec3dp{0.0, 0.0, a, 0.0, 0.0, 0.0});
        };
        auto const R = rexp(0.5 * bivec3dp{eps_roll, 0.0, 0.0, 0.0, 0.0, 0.0});
        value_t worst_rel = 0.0, worst_abs = 0.0;
        for (value_t const psi : {0.0, 90.0, 180.0, 270.0, 350.0, 359.0, 365.0}) {
            auto const Y = yaw(psi * deg);
            auto const M = rgpr(Y, R);
            bivec3dp const e = motion_error_body(Y, M);
            worst_rel = std::max(
                worst_rel,
                equivariance::max_diff(e, bivec3dp{eps_roll, 0.0, 0.0, 0.0, 0.0, 0.0}));
            bivec3dp const d = 2.0 * rlog(M) - 2.0 * rlog(Y);
            worst_abs = std::max(worst_abs, std::hypot(d.vx, d.vy) / eps_roll);
            // what the errors mean: M = M_ref (x) rexp(e_b / 2) = rexp(e_w / 2) (x) M_ref
            CHECK(is_same_motion(rgpr(Y, rexp(0.5 * e)), M));
            CHECK(is_same_motion(rgpr(rexp(0.5 * motion_error_world(Y, M)), Y), M));
        }
        CHECK(worst_rel < 1.0e-14);
        CHECK(worst_abs > 100.0);
        fmt::println(
            "    the roll from the relative motor: off by {:.1e}; from the absolute "
            "logarithms: up to {:.0f} times the roll",
            worst_rel, worst_abs);

        // THE EQUIVARIANCE PATTERN (ga_test_equivariance.hpp): the whole scene turned by
        // 0, 90, 180, 270, 365 and 719 degrees about a vertical axis off the origin. The
        // body error does not change, the world error turns with the scene; the
        // difference of absolute logarithms does neither.
        bivec3dp const L1 =
            unitize(wdg(vec3dp{1.0, 2.0, 0.5, 1.0}, vec3dp{1.3, 1.6, 1.4, 1.0}));
        bivec3dp const L2 =
            unitize(wdg(vec3dp{-0.5, 0.2, 0.0, 1.0}, vec3dp{0.1, 0.9, 0.3, 1.0}));
        std::vector<mvec3dp_e> const scene3{get_motor(L1, 0.8, 0.3),
                                            get_motor(L2, 2.6, -0.4)};
        auto body3 = [](std::vector<mvec3dp_e> const& m) {
            return motion_error_body(m[0], m[1]);
        };
        auto world3 = [](std::vector<mvec3dp_e> const& m) {
            return motion_error_world(m[0], m[1]);
        };
        auto absolute3 = [](std::vector<mvec3dp_e> const& m) {
            return bivec3dp(2.0 * rlog(m[1]) - 2.0 * rlog(m[0]));
        };
        CHECK(equivariance::body_invariance3dp(body3, scene3) < 1.0e-12);
        CHECK(equivariance::world_covariance3dp(world3, scene3) < 1.0e-12);
        CHECK(equivariance::body_invariance3dp(absolute3, scene3) > 0.1);
        CHECK(equivariance::world_covariance3dp(absolute3, scene3) > 0.1);
        std::vector<mvec2dp_u> const scene2{get_motor(vec2dp{1.0, 2.0, 1.0}, 0.8),
                                            get_motor(vec2dp{-0.5, 0.3, 1.0}, 2.6)};
        auto body2 = [](std::vector<mvec2dp_u> const& m) {
            return motion_error_body(m[0], m[1]);
        };
        auto world2 = [](std::vector<mvec2dp_u> const& m) {
            return motion_error_world(m[0], m[1]);
        };
        auto absolute2 = [](std::vector<mvec2dp_u> const& m) {
            return vec2dp(2.0 * rlog(m[1]) - 2.0 * rlog(m[0]));
        };
        CHECK(equivariance::body_invariance2dp(body2, scene2) < 1.0e-12);
        CHECK(equivariance::world_covariance2dp(world2, scene2) < 1.0e-12);
        CHECK(equivariance::body_invariance2dp(absolute2, scene2) > 0.1);
        CHECK(equivariance::world_covariance2dp(absolute2, scene2) > 0.1);
        fmt::println(
            "    turned scenes: body error invariant to {:.1e}, world error covariant "
            "to {:.1e}; absolute logs off by {:.2f}",
            equivariance::body_invariance3dp(body3, scene3),
            equivariance::world_covariance3dp(world3, scene3),
            equivariance::body_invariance3dp(absolute3, scene3));
        fmt::println("");
    }

    TEST_CASE("pga2dp / pga3dp: reading a motion off a motor -- following the turns")
    {
        fmt::println(
            "pga2dp / pga3dp: rlog_continued follows a heading through four turns");

        // A heading growing by 37 degrees per read, about 4 turns, and one falling by 23
        // per read: rlog_continued, fed its own previous value, reads the accumulated
        // angle; rlog alone wraps (the falsification). A screw keeps its axis and pitch.
        value_t const deg = pi / 180.0;
        bivec3dp const Lz =
            unitize(wdg(vec3dp{0.7, -0.4, 0.0, 1.0}, vec3dp{0.7, -0.4, 1.0, 1.0}));
        for (value_t const step : {37.0, -23.0}) {
            CAPTURE(step);
            bivec3dp L3{};
            vec2dp L2{};
            value_t worst = 0.0, worst_plain = 0.0;
            for (int k = 1; k <= 40; ++k) {
                value_t const psi = k * step * deg;
                auto const S = get_motor(Lz, psi, 0.25); // a screw about a vertical axis
                L3 = rlog_continued(S, L3);
                CHECK(is_same_motion(rexp(L3), S));
                auto const P = get_motor(vec2dp{1.0, 2.0, 1.0}, psi);
                L2 = rlog_continued(P, L2);
                CHECK(is_same_motion(rexp(L2), P));
                worst = std::max(
                    {worst, std::abs(2.0 * L3.vz - psi), std::abs(2.0 * L2.z - psi)});
                worst_plain = std::max(worst_plain, std::abs(2.0 * rlog(S).vz - psi));
            }
            CHECK(worst < 1.0e-12);
            CHECK(worst_plain > 6.0); // whole turns lost
            fmt::println(
                "    {:+.0f} deg per read, 40 reads: the continued angle off by {:.1e} "
                "rad, rlog alone by up to {:.1f}",
                step, worst, worst_plain);
        }
        // exactly at whole turns (a read every 90 degrees) the motor is +-identity and
        // has no axis; the turns continue along the previous axis. Exact for a rotation
        // and for a screw, whose translation runs along the axis
        bivec3dp L3{};
        vec2dp L2{};
        for (int k = 1; k <= 41; ++k) {
            auto const S = get_motor(Lz, k * 90.0 * deg, 0.25);
            L3 = rlog_continued(S, L3);
            CHECK(is_same_motion(rexp(L3), S));
            auto const P = get_motor(vec2dp{1.0, 2.0, 1.0}, k * 90.0 * deg);
            L2 = rlog_continued(P, L2);
            CHECK(is_same_motion(rexp(L2), P));
        }
        CHECK(std::abs(2.0 * L3.vz - 41 * 90.0 * deg) < 1.0e-12);
        CHECK(std::abs(2.0 * L2.z - 41 * 90.0 * deg) < 1.0e-12);
        // a translation, with no turns behind it, is rlog's
        auto const T = get_motor(vec3dp{0.4, -1.1, 0.3, 0.0});
        CHECK(rlog_continued(T, bivec3dp{}) == rlog(T));
        auto const T2 = get_motor(vec2dp{0.4, -1.1, 0.0});
        CHECK(rlog_continued(T2, vec2dp{}) == rlog(T2));
        fmt::println("");
    }
}
