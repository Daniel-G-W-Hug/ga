// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include <variant>
#include <vector>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_sta.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::sta; // use specific operations of STA (Space-Time Algebra)


/////////////////////////////////////////////////////////////////////////////////////////
// Space-Time Algebra: STA4ds (G<1,3,0>)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("STA 3D Tests")
{

    TEST_CASE("G<1,3,0>: sta4ds")
    {
        fmt::println("G<1,3,0>: sta4ds");
        // Space-Time Algebra with signature (-,-,-,+)
        const algebra<1, 3, 0> alg;
        CHECK(alg.p() == 1);
        CHECK(alg.n() == 3);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   sta4ds: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   sta4ds: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    TEST_CASE("Scalar4ds: comparison")
    {
        fmt::println("Scalar4ds: comparison");
        scalar4ds s1{3.0};
        scalar4ds s2{3.0};
        pscalar4ds ps1{3.0};
        pscalar4ds ps2{4.0};

        // same tag comparisons should work
        CHECK(s1 == s2);
        CHECK(s1 != scalar4ds{4.0});
        CHECK(s1 < scalar4ds{4.0});
        CHECK(s1 <= scalar4ds{3.0});
        CHECK(scalar4ds{4.0} > s1);
        CHECK(scalar4ds{3.0} >= s1);
        CHECK(ps1 != ps2);

        // different floating point types with same tag should work
        Scalar_t<float, scalar4ds_tag> sf{3.0f};
        CHECK(s1 == sf); // double vs float with same tag
    }

    ////////////////////////////////////////////////////////////////////////////////
    // complement mathematical properties tests
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: cmpl composition properties")
    {
        fmt::println("G<1,3,0>: cmpl composition properties");

        // For STA4ds (even-dimensional):
        // r_cmpl(l_cmpl(u)) = u
        // l_cmpl(r_cmpl(u)) = u

        scalar4ds s{3.0};
        vec4ds v{2.0, 5.0, 7.0, 9.0};
        bivec4ds B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds t{1.0, 2.0, 3.0, 4.0};
        pscalar4ds ps{11.0};
        mvec4ds_u M_u{vec4ds{1.0, 2.0, 3.0, 4.0}, trivec4ds{11.0, 12.0, 13.0, 14.0}};
        mvec4ds_e M_e{scalar4ds{1.0}, bivec4ds{2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
                      pscalar4ds{8.0}};
        mvec4ds M{scalar4ds{1.0}, vec4ds{2.0, 3.0, 4.0, 5.0},
                  bivec4ds{6.0, 7.0, 8.0, 9.0, 10.0, 11.0},
                  trivec4ds{12.0, 13.0, 14.0, 15.0}, pscalar4ds{16.0}};

        // l_cmpl(r_cmpl(u)) == u
        CHECK(l_cmpl(r_cmpl(s)) == s);
        CHECK(l_cmpl(r_cmpl(v)) == v);
        CHECK(l_cmpl(r_cmpl(B)) == B);
        CHECK(l_cmpl(r_cmpl(t)) == t);
        CHECK(l_cmpl(r_cmpl(ps)) == ps);
        CHECK(l_cmpl(r_cmpl(M_u)) == M_u);
        CHECK(l_cmpl(r_cmpl(M_e)) == M_e);
        CHECK(l_cmpl(r_cmpl(M)) == M);

        // r_cmpl(l_cmpl(u)) == u
        CHECK(r_cmpl(l_cmpl(s)) == s);
        CHECK(r_cmpl(l_cmpl(v)) == v);
        CHECK(r_cmpl(l_cmpl(B)) == B);
        CHECK(r_cmpl(l_cmpl(t)) == t);
        CHECK(r_cmpl(l_cmpl(ps)) == ps);
        CHECK(r_cmpl(l_cmpl(M_u)) == M_u);
        CHECK(r_cmpl(l_cmpl(M_e)) == M_e);
        CHECK(r_cmpl(l_cmpl(M)) == M);

        fmt::println("l_cmpl(r_cmpl(u)) == r_cmpl(l_cmpl(u)) == u");

        // l_cmpl(u) ^ u == (Euclidean coordinate sum of squares) * I_4ds, and
        // u ^ r_cmpl(u) == (Euclidean coordinate sum of squares) * I_4ds.
        //
        // NOTE: The test is only valid in this metric-independent form. The
        // complement is a purely combinatorial (non-metric) operation: it maps
        // each basis blade to its complementary blade with a sign chosen so that
        // l_cmpl(e_i) ^ e_i == I_4ds for every basis vector. Hence the resulting
        // coefficient is always the *Euclidean* sum of the coordinates squared,
        // regardless of the STA signature. The familiar EGA form
        // "l_cmpl(u) ^ u == nrm_sq(u) * I" holds only when every basis vector
        // squares to +1; in the mixed STA metric nrm_sq(u) carries the signature
        // and differs from the Euclidean sum, so that form no longer applies.
        value_t const sq_s = value_t(s) * value_t(s);
        value_t const sq_v = v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
        value_t const sq_B = B.vx * B.vx + B.vy * B.vy + B.vz * B.vz + B.mx * B.mx +
                             B.my * B.my + B.mz * B.mz;
        value_t const sq_t = t.x * t.x + t.y * t.y + t.z * t.z + t.w * t.w;
        value_t const sq_ps = value_t(ps) * value_t(ps);

        CHECK(wdg(l_cmpl(s), s) == sq_s * I_4ds);
        CHECK(wdg(l_cmpl(v), v) == sq_v * I_4ds);
        CHECK(wdg(l_cmpl(B), B) == sq_B * I_4ds);
        CHECK(wdg(l_cmpl(t), t) == sq_t * I_4ds);
        CHECK(wdg(l_cmpl(ps), ps) == sq_ps * I_4ds);

        CHECK(wdg(s, r_cmpl(s)) == sq_s * I_4ds);
        CHECK(wdg(v, r_cmpl(v)) == sq_v * I_4ds);
        CHECK(wdg(B, r_cmpl(B)) == sq_B * I_4ds);
        CHECK(wdg(t, r_cmpl(t)) == sq_t * I_4ds);
        CHECK(wdg(ps, r_cmpl(ps)) == sq_ps * I_4ds);

        fmt::println("l_cmpl(u) ^ u == u ^ r_cmpl(u) == (sum of squares) * I_4ds");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // metric signature tests (dot / gpr / nrm_sq)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: metric signature (dot / gpr / nrm_sq)")
    {
        fmt::println("G<1,3,0>: metric signature (dot / gpr / nrm_sq)");

        // Unlike the complement (a purely combinatorial, metric-blind operation),
        // dot(), the geometric product and nrm_sq() all carry the metric. These
        // tests pin the signature (-,-,-,+) and its exomorphism extension to the
        // higher grades, and cross-check the three independent code paths against
        // each other.

        // (1) Signature pinned via nrm_sq() on each unit basis blade. Under the
        //     P-unify metric nrm_sq is the reverse-norm <X ~X>_0 = the metric
        //     exomorphism P(e_S) = prod_{i in S} g_i (g1^2 = g2^2 = g3^2 = -1; g4^2 = +1;
        //     higher grades via G(a^b) = G(a)^G(b)). NOTE: this is the METRIC, not the
        //     geometric square -- causal character (timelike/spacelike) reads gr0(X*X) =
        //     sigma(k)*P, which differs in sign at grades 2,3 (see the geom_sq guard).
        // vectors: g1,g2,g3 -> -1; g4 -> +1
        CHECK(nrm_sq(g1_4ds) == -1.0);
        CHECK(nrm_sq(g2_4ds) == -1.0);
        CHECK(nrm_sq(g3_4ds) == -1.0);
        CHECK(nrm_sq(g4_4ds) == 1.0);
        // bivectors: P = product of the two vector metrics (with g4 -> -1; without -> +1)
        CHECK(nrm_sq(g14_4ds) == -1.0);
        CHECK(nrm_sq(g24_4ds) == -1.0);
        CHECK(nrm_sq(g34_4ds) == -1.0);
        CHECK(nrm_sq(g23_4ds) == 1.0);
        CHECK(nrm_sq(g31_4ds) == 1.0);
        CHECK(nrm_sq(g12_4ds) == 1.0);
        // trivectors: P = product of the three vector metrics (with g4 -> +1; g123 -> -1)
        CHECK(nrm_sq(g234_4ds) == 1.0);
        CHECK(nrm_sq(g314_4ds) == 1.0);
        CHECK(nrm_sq(g124_4ds) == 1.0);
        CHECK(nrm_sq(g123_4ds) == -1.0);
        // pseudoscalar: det(g) = -1
        CHECK(nrm_sq(I_4ds) == -1.0);

        // (2) Defining metric law via the geometric product: g_i (.) g_i = g_ii.
        //     Routes the signature through gpr(), an independent code path.
        CHECK(value_t(gr0(g1_4ds * g1_4ds)) == -1.0);
        CHECK(value_t(gr0(g2_4ds * g2_4ds)) == -1.0);
        CHECK(value_t(gr0(g3_4ds * g3_4ds)) == -1.0);
        CHECK(value_t(gr0(g4_4ds * g4_4ds)) == 1.0);

        // (3) Cross-check the three metric-bearing paths against each other:
        //     codegen dot(), hand-written nrm_sq() and (for vectors) the
        //     contraction axiom v (.) v = v . v carried by gpr().
        scalar4ds s{3.0};
        vec4ds v{2.0, 5.0, 7.0, 9.0};
        bivec4ds B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds t{1.0, 2.0, 3.0, 4.0};
        pscalar4ds ps{11.0};

        CHECK(value_t(dot(s, s)) == nrm_sq(s));
        CHECK(value_t(dot(v, v)) == nrm_sq(v));
        CHECK(value_t(dot(B, B)) == nrm_sq(B));
        CHECK(value_t(dot(t, t)) == nrm_sq(t));
        CHECK(value_t(dot(ps, ps)) == nrm_sq(ps));

        // vector contraction axiom: v (.) v == v . v (scalar), v ^ v == 0
        CHECK(value_t(gr0(v * v)) == nrm_sq(v));
        CHECK(gr2(v * v) == bivec4ds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

        fmt::println("signature (-,-,-,+) verified: dot == nrm_sq, g_i (.) g_i == g_ii");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // geometric-square helper guard: detail::sta4ds_geom_sq(X) == gr0(X*X)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: detail::sta4ds_geom_sq matches gr0(X*X)")
    {
        fmt::println("G<1,3,0>: detail::sta4ds_geom_sq matches gr0(X*X)");

        // detail::sta4ds_geom_sq is a hand-written component transcription of the
        // GEOMETRIC square B^2 = gr0(X*X). The causal-character predicates
        // (is_timelike / is_spacelike / is_lightlike) read it INSTEAD of nrm_sq -- which
        // under the P-unify metric is the reverse-norm = sigma(k)*B^2, sign-flipped at
        // grades 2,3. This guard pins the transcription to the actual geometric product,
        // so any future signature / basis-ordering / gpr change that desyncs them is
        // caught here (integer-valued components keep the == comparison exact).
        auto chk = [](auto const& X) {
            CHECK(detail::sta4ds_geom_sq(X) == value_t(gr0(X * X)));
        };

        // grade 1 (vectors)
        chk(g1_4ds);
        chk(g2_4ds);
        chk(g3_4ds);
        chk(g4_4ds);
        // grade 2 (bivectors)
        chk(g14_4ds);
        chk(g24_4ds);
        chk(g34_4ds);
        chk(g23_4ds);
        chk(g31_4ds);
        chk(g12_4ds);
        // grade 3 (trivectors)
        chk(g234_4ds);
        chk(g314_4ds);
        chk(g124_4ds);
        chk(g123_4ds);

        // general (non-basis) blades: the full quadratic form must agree too
        auto const v = vec4ds{1.0, -2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, -3.0, 4.0, -5.0, 6.0};
        auto const t = trivec4ds{-2.0, 1.0, 4.0, -3.0};
        chk(v);
        chk(B);
        chk(t);

        // P-unify relation to nrm_sq: geom_sq = sigma(k) * nrm_sq -- equal for vectors
        // (sigma(1) = +1), sign-flipped for bi-/trivectors (sigma(2) = sigma(3) = -1)
        CHECK(detail::sta4ds_geom_sq(v) == nrm_sq(v));
        CHECK(detail::sta4ds_geom_sq(B) == -nrm_sq(B));
        CHECK(detail::sta4ds_geom_sq(t) == -nrm_sq(t));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // causal character, nrm and normalization of blades
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: causal character, nrm and normalization of blades")
    {
        fmt::println("G<1,3,0>: causal character, nrm and normalization of blades");

        // ---- causal classification (grades 1..n-1: Vec, BiVec, TriVec) ----
        // sourced from the geometric square B^2 = gr0(X*X) (NOT nrm_sq):
        // timelike: B^2 > 0 ; spacelike: B^2 < 0 ; lightlike: B^2 == 0
        // vectors: g1,g2,g3 -> spacelike ; g4 -> timelike
        CHECK(is_spacelike(g1_4ds));
        CHECK(is_spacelike(g2_4ds));
        CHECK(is_spacelike(g3_4ds));
        CHECK(is_timelike(g4_4ds));
        // null vector g1 + g4 -> B^2 = -1 + 1 = 0
        vec4ds vl{1.0, 0.0, 0.0, 1.0};
        CHECK(is_lightlike(vl));
        CHECK(!is_timelike(vl)); // mutually exclusive
        CHECK(!is_spacelike(vl));

        // bivectors: g14,g24,g34 -> timelike ; g23,g31,g12 -> spacelike
        CHECK(is_timelike(g14_4ds));
        CHECK(is_timelike(g24_4ds));
        CHECK(is_timelike(g34_4ds));
        CHECK(is_spacelike(g23_4ds));
        CHECK(is_spacelike(g31_4ds));
        CHECK(is_spacelike(g12_4ds));
        // null bivector g14 + g23 -> B^2 = +1 - 1 = 0
        bivec4ds Bl{1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        CHECK(is_lightlike(Bl));

        // trivectors: g234,g314,g124 -> spacelike ; g123 -> timelike
        CHECK(is_spacelike(g234_4ds));
        CHECK(is_spacelike(g314_4ds));
        CHECK(is_spacelike(g124_4ds));
        CHECK(is_timelike(g123_4ds));
        // null trivector g234 + g123 -> B^2 = -1 + 1 = 0
        trivec4ds tl{1.0, 0.0, 0.0, 1.0};
        CHECK(is_lightlike(tl));

        // a mixed k-vector takes the character of its dominant contribution:
        // (timelike weight = sum of timelike-blade squares, likewise spacelike)
        CHECK(is_timelike(vec4ds{1.0, 0.0, 0.0, 2.0}));  // time 2^2 > space 1^2
        CHECK(is_spacelike(vec4ds{2.0, 0.0, 0.0, 1.0})); // space 2^2 > time 1^2
        // bivector: g14 timelike (+1), g23 spacelike (-1)
        CHECK(is_spacelike(bivec4ds{1.0, 0.0, 0.0, 2.0, 0.0, 0.0})); // g23 dominates
        CHECK(is_timelike(bivec4ds{2.0, 0.0, 0.0, 1.0, 0.0, 0.0}));  // g14 dominates
        // trivector: g234 spacelike (-1), g123 timelike (+1)
        CHECK(is_spacelike(trivec4ds{2.0, 0.0, 0.0, 1.0})); // g234 dominates
        CHECK(is_timelike(trivec4ds{1.0, 0.0, 0.0, 2.0}));  // g123 dominates

        // ---- nrm(u) == sqrt(|nrm_sq(u)|), always >= 0, zero for null ----
        scalar4ds s{3.0};
        vec4ds v{2.0, 5.0, 7.0, 9.0}; // nrm_sq = 4+25+49-81 = -3
        bivec4ds B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds t{1.0, 2.0, 3.0, 4.0};
        pscalar4ds ps{11.0};

        CHECK(value_t(nrm(s)) == doctest::Approx(std::sqrt(std::abs(nrm_sq(s)))));
        CHECK(value_t(nrm(v)) == doctest::Approx(std::sqrt(std::abs(nrm_sq(v)))));
        CHECK(value_t(nrm(B)) == doctest::Approx(std::sqrt(std::abs(nrm_sq(B)))));
        CHECK(value_t(nrm(t)) == doctest::Approx(std::sqrt(std::abs(nrm_sq(t)))));
        CHECK(value_t(nrm(ps)) == doctest::Approx(std::sqrt(std::abs(nrm_sq(ps)))));
        // scalar / pseudoscalar magnitude is just |coeff|
        CHECK(value_t(nrm(s)) == doctest::Approx(3.0));
        CHECK(value_t(nrm(ps)) == doctest::Approx(11.0));
        // null blades have zero magnitude
        CHECK(value_t(nrm(vl)) == doctest::Approx(0.0));
        CHECK(value_t(nrm(Bl)) == doctest::Approx(0.0));
        CHECK(value_t(nrm(tl)) == doctest::Approx(0.0));

        // ---- normalize(): scales to nrm_sq == sign(nrm_sq) == +/-1, lightlike unchanged
        // direction is preserved (scale = 1/nrm > 0, never flips sign). The resulting
        // sign is that of the reverse-norm nrm_sq = P (NOT the causal character B^2):
        // they agree for vectors but are opposite for bi-/trivectors.
        // timelike vector 3*g4: nrm_sq = 9 -> +1 (P agrees with B^2 for vectors)
        vec4ds vt{0.0, 0.0, 0.0, 3.0};
        CHECK(nrm_sq(normalize(vt)) == doctest::Approx(1.0));
        CHECK(normalize(vt).w > 0.0);
        // spacelike vector 3*g1: nrm_sq = -9 -> -1
        vec4ds vs{3.0, 0.0, 0.0, 0.0};
        CHECK(nrm_sq(normalize(vs)) == doctest::Approx(-1.0));
        CHECK(normalize(vs).x > 0.0);
        // lightlike vector: returned unchanged (exact, no scaling applied)
        CHECK(normalize(vl) == vl);

        // bivector: 2*g14 (timelike, B^2>0) has nrm_sq = P = -1 ; 2*g23 (spacelike,
        // B^2<0) has nrm_sq = P = +1 ; null unchanged  (signs opposite to causal char)
        CHECK(nrm_sq(normalize(bivec4ds{2.0, 0.0, 0.0, 0.0, 0.0, 0.0})) ==
              doctest::Approx(-1.0));
        CHECK(nrm_sq(normalize(bivec4ds{0.0, 0.0, 0.0, 2.0, 0.0, 0.0})) ==
              doctest::Approx(1.0));
        CHECK(normalize(Bl) == Bl);

        // trivector: 2*g123 (timelike, B^2>0) has nrm_sq = P = -1 ; 2*g234 (spacelike,
        // B^2<0) has nrm_sq = P = +1 ; null unchanged
        CHECK(nrm_sq(normalize(trivec4ds{0.0, 0.0, 0.0, 2.0})) == doctest::Approx(-1.0));
        CHECK(nrm_sq(normalize(trivec4ds{2.0, 0.0, 0.0, 0.0})) == doctest::Approx(1.0));
        CHECK(normalize(tl) == tl);

        // scalar normalize still works (no causal logic): scales to unit sign
        CHECK(value_t(normalize(scalar4ds{3.0})) == doctest::Approx(1.0));
        CHECK(value_t(normalize(scalar4ds{-5.0})) == doctest::Approx(-1.0));

        fmt::println("causal predicates, nrm = sqrt(|nrm_sq|), normalize -> nrm_sq == "
                     "+/-1 (null unchanged)");
    }

    TEST_CASE("G<1,3,0>: inverses (inv w.r.t. geometric product)")
    {
        fmt::println("G<1,3,0>: inverses (inv w.r.t. geometric product)");

        // STA has a non-degenerate (Lorentzian) metric, so every non-null blade and
        // every non-zero-divisor multivector has a geometric inverse with
        //     u * inv(u) == inv(u) * u == 1   (the scalar, identity of gpr)
        // Unlike the degenerate PGA pseudoscalars, the STA pseudoscalar is also
        // invertible (I_4ds^2 = -1, so inv(I_4ds) = -I_4ds).

        scalar4ds s1{3.2};
        vec4ds v1{2.0, 1.0, 2.0, -2.0}; // nrm_sq = -4-1-4+4 = -5 (non-null)
        bivec4ds b1{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds t1{1.0, 2.0, 3.0, 4.0};
        pscalar4ds ps1{-4.7};
        mvec4ds_e mve1{s1, b1, ps1};
        mvec4ds_u mvu1{v1, t1};
        mvec4ds mv1{s1, v1, b1, t1, ps1};

        // scalar
        CHECK(value_t(s1 * inv(s1)) == doctest::Approx(1.0));
        CHECK(value_t(inv(s1) - rev(s1) / nrm_sq(s1)) == doctest::Approx(0.0));

        // vector  (v * inv(v) -> even-grade multivector)
        CHECK(value_t(gr0(v1 * inv(v1))) == doctest::Approx(1.0));
        CHECK(value_t(nrm(gr2(v1 * inv(v1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr4(v1 * inv(v1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(inv(v1) - rev(v1) / nrm_sq(v1))) == doctest::Approx(0.0));

        // bivector
        CHECK(value_t(gr0(b1 * inv(b1))) == doctest::Approx(1.0));
        CHECK(value_t(nrm(gr2(b1 * inv(b1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr4(b1 * inv(b1)))) == doctest::Approx(0.0));

        // trivector
        CHECK(value_t(gr0(t1 * inv(t1))) == doctest::Approx(1.0));
        CHECK(value_t(nrm(gr2(t1 * inv(t1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr4(t1 * inv(t1)))) == doctest::Approx(0.0));

        // pseudoscalar (invertible: non-degenerate metric, I_4ds^2 = -1)
        CHECK(value_t(ps1 * inv(ps1)) == doctest::Approx(1.0));
        CHECK(inv(pscalar4ds{1.0}) == -I_4ds); // inv(I) = -I

        // even-grade multivector
        CHECK(value_t(gr0(mve1 * inv(mve1))) == doctest::Approx(1.0));
        CHECK(value_t(nrm(gr2(mve1 * inv(mve1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr4(mve1 * inv(mve1)))) == doctest::Approx(0.0));

        // odd-grade multivector (u * inv(u) -> even-grade multivector)
        CHECK(value_t(gr0(mvu1 * inv(mvu1))) == doctest::Approx(1.0));
        CHECK(value_t(nrm(gr2(mvu1 * inv(mvu1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr4(mvu1 * inv(mvu1)))) == doctest::Approx(0.0));

        // full multivector (left and right inverse coincide)
        CHECK(value_t(gr0(mv1 * inv(mv1))) == doctest::Approx(1.0));
        CHECK(value_t(nrm(gr1(mv1 * inv(mv1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr2(mv1 * inv(mv1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr3(mv1 * inv(mv1)))) == doctest::Approx(0.0));
        CHECK(value_t(nrm(gr4(mv1 * inv(mv1)))) == doctest::Approx(0.0));
        CHECK(value_t(gr0(inv(mv1) * mv1)) == doctest::Approx(1.0)); // left inverse
    }

    ////////////////////////////////////////////////////////////////////////////////
    // sanity checks for the contraction / expansion / rwdg / cmt products
    // (identities mirrored from the EGA3D test suite, adapted to G(1,3,0))
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: contraction / rwdg / cmt / expansion sanity checks")
    {
        fmt::println("G<1,3,0>: contraction / rwdg / cmt / expansion sanity checks");

        vec4ds v1{1.0, 2.0, 3.0, 4.0};
        vec4ds v2{5.0, 6.0, 7.0, 8.0};
        bivec4ds B1{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        bivec4ds B2{2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
        trivec4ds t{1.0, 2.0, 3.0, 4.0};
        vec4ds v3{2.0, 0.0, 1.0, 3.0};
        bivec4ds B3{1.0, 0.0, 2.0, 0.0, 3.0, 1.0};
        trivec4ds t2{4.0, 3.0, 2.0, 1.0};
        scalar4ds sc{3.0};
        mvec4ds M1{scalar4ds{1.0}, v1, B1, t, pscalar4ds{2.0}};
        mvec4ds M2{scalar4ds{2.0}, v2, B2, trivec4ds{4.0, 3.0, 2.0, 1.0},
                   pscalar4ds{5.0}};

        // ---- contractions ----
        // for equal grades the contraction reduces to the metric inner product
        // (holds at every grade, not just vectors)
        CHECK((v1 << v2) == dot(v1, v2));
        CHECK((v1 >> v2) == dot(v1, v2));
        CHECK((B1 << B2) == dot(B1, B2));
        CHECK((B1 >> B2) == dot(B1, B2));
        CHECK((t << t2) == dot(t, t2));
        CHECK((t >> t2) == dot(t, t2));
        // the contraction carries the metric: timelike g4 -> +1, spacelike g1 -> -1
        CHECK((g4_4ds << g4_4ds) == scalar4ds{1.0});
        CHECK((g1_4ds << g1_4ds) == scalar4ds{-1.0});
        // scalar contraction is just scalar multiplication
        CHECK((sc << B1) == value_t(sc) * B1);
        CHECK((B1 >> sc) == value_t(sc) * B1);
        // geometric product decomposes into contraction + wedge (grade split)
        CHECK(v1 * v2 == (v2 >> v1) + wdg(v1, v2));
        CHECK(B1 * v1 == (v1 << B1) + wdg(B1, v1));
        CHECK(v1 * B1 == (B1 >> v1) + wdg(v1, B1));
        // Lengyel's defining identity: contraction == rwdg(dual(a), b), exact at every
        // grade. grade(a)=2 (B1<<B2) is the case that pinned the bivector dual sign.
        CHECK((v1 << B1) == rwdg(l_dual(v1), B1));
        CHECK((B1 << B2) == rwdg(l_dual(B1), B2));
        CHECK((t << t) == rwdg(l_dual(t), t));
        CHECK((B1 >> B2) == rwdg(B1, r_dual(B2)));

        // ---- regressive wedge (meet) ----
        // defining identity: rwdg(a,b) = l_cmpl(wdg(r_cmpl(a), r_cmpl(b)))
        CHECK(rwdg(B1, B2) == l_cmpl(wdg(r_cmpl(B1), r_cmpl(B2))));
        // the pseudoscalar is the neutral element of the meet
        CHECK(rwdg(I_4ds, sc) == sc);
        CHECK(rwdg(sc, I_4ds) == sc);
        // complement turns the meet into the join (De Morgan duality)
        CHECK(r_cmpl(rwdg(B1, B2)) == wdg(r_cmpl(B1), r_cmpl(B2)));

        // ---- commutator (asymmetric part of the geometric product) ----
        CHECK(cmt(v1, v2) == wdg(v1, v2));               // = wedge for two vectors
        CHECK(cmt(M1, M2) == 0.5 * (M1 * M2 - M2 * M1)); // definition, full mv
        CHECK(cmt(v1, v2) == -cmt(v2, v1));              // antisymmetric
        CHECK(cmt(B1, B2) == -cmt(B2, B1));
        // the bivectors close under the commutator -> Lorentz Lie algebra so(1,3)
        CHECK(cmt(B1, B2) == gr2(0.5 * (B1 * B2 - B2 * B1)));
        CHECK(cmt(cmt(B1, B2), B3) + cmt(cmt(B2, B3), B1) + cmt(cmt(B3, B1), B2) ==
              bivec4ds{}); // Jacobi identity

        // ---- expansions (wdg-duals of the contractions) ----
        // ground-truth values independently computed (NOT the self-referential wrapper
        // form l_expand4ds(a,b)==wdg(l_dual(a),b), which can't catch a dual sign bug).
        // Under P-unify these negate vs the old Q-dual values, since the bivector dual
        // flips sign at grade 2 (verified by hand against the named-blade signs).
        CHECK(l_expand4ds(B1, v1) == trivec4ds{-7.0, -2.0, -15.0, -14.0});
        CHECK(r_expand4ds(v1, B1) == trivec4ds{-7.0, -2.0, -15.0, -14.0});
        // meet with the pseudoscalar dual: l_dual(I_4ds) is a unit scalar, so
        // expanding a blade by I_4ds returns it scaled by det(g) = -1
        CHECK(l_expand4ds(I_4ds, t) == -t);
        CHECK(r_expand4ds(t, I_4ds) == -t);

        // ---- bilinearity (regression guards against gross wiring breakage) ----
        CHECK(((v1 + v2) << B1) == (v1 << B1) + (v2 << B1));
        CHECK(rwdg(B1, B2 + B3) == rwdg(B1, B2) + rwdg(B1, B3));
        CHECK(cmt(v1, v2 + v3) == cmt(v1, v2) + cmt(v1, v3));
        CHECK(l_expand4ds(B1 + B2, v1) == l_expand4ds(B1, v1) + l_expand4ds(B2, v1));

        fmt::println("contraction == dot (equal grades), gpr = contraction + wedge, "
                     "rwdg/cmt/expansion identities hold");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // transcription gate: dual == metric * complement, per basis blade
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: dual == metric * complement (transcription gate)")
    {
        fmt::println("G<1,3,0>: dual == metric * complement (transcription gate)");

        // The dual is "complement after multiplication with the extended metric":
        //   l_dual(e) = l_cmpl(G*e),   r_dual(e) = r_cmpl(G*e).
        // For a unit basis blade e, G*e = nrm_sq(e)*e, so the relation must hold
        // element-wise. This is a direct guard against sign-transcription errors
        // when the dual/complement tables are hand-coded from the output of
        // ga_prdxpr_rule_generator_test -- it catches exactly the kind of bivector
        // dual sign slip that previously went unnoticed (the rule generator is the
        // source of truth; this pins the hand-written code to it).
        auto gate = [](auto const& e) {
            CHECK(l_dual(e) == nrm_sq(e) * l_cmpl(e));
            CHECK(r_dual(e) == nrm_sq(e) * r_cmpl(e));
        };
        gate(one_4ds);
        gate(g1_4ds);
        gate(g2_4ds);
        gate(g3_4ds);
        gate(g4_4ds);
        gate(g14_4ds);
        gate(g24_4ds);
        gate(g34_4ds);
        gate(g23_4ds);
        gate(g31_4ds);
        gate(g12_4ds);
        gate(g234_4ds);
        gate(g314_4ds);
        gate(g124_4ds);
        gate(g123_4ds);
        gate(I_4ds);

        fmt::println("l_dual(e) == nrm_sq(e)*l_cmpl(e) and r_dual analogue verified");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 1: rotor exponential (rotations & Lorentz boosts)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: rotor exponential exp() (rotation & boost)")
    {
        fmt::println("G<1,3,0>: rotor exponential exp() (rotation & boost)");

        value_t const a = 0.7;

        // exp(0) == identity rotor (scalar 1, no bivector/pseudoscalar)
        CHECK(exp(bivec4ds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) ==
              mvec4ds_e{scalar4ds{1.0}, bivec4ds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}});

        // spatial plane g12 (B^2 < 0): exp(a*g12) = cos(a) + sin(a) g12  (circular)
        auto const Rs = exp(a * g12_4ds);
        CHECK(value_t(gr0(Rs)) == doctest::Approx(std::cos(a)));
        CHECK(gr2(Rs).mz == doctest::Approx(std::sin(a))); // g12 component

        // boost plane g14 (B^2 > 0): exp(a*g14) = cosh(a) + sinh(a) g14  (hyperbolic)
        auto const Rb = exp(a * g14_4ds);
        CHECK(value_t(gr0(Rb)) == doctest::Approx(std::cosh(a)));
        CHECK(gr2(Rb).vx == doctest::Approx(std::sinh(a))); // g14 component

        // builders carry the half-angle for the sandwich (validated semantically in
        // step 2); here check their form and that they are unit rotors rev(R)*R == 1.
        auto const Rr = get_rotor(g12_4ds, a); // cos(a/2) + sin(a/2) g12
        CHECK(value_t(gr0(Rr)) == doctest::Approx(std::cos(a / 2)));
        CHECK(gr2(Rr).mz == doctest::Approx(std::sin(a / 2)));
        CHECK(value_t(gr0(rev(Rr) * Rr)) == doctest::Approx(1.0));
        CHECK(nrm_sq(gr2(rev(Rr) * Rr)) == doctest::Approx(0.0)); // no leftover bivector

        auto const Rbo = get_boost(g14_4ds, a); // cosh(a/2) + sinh(a/2) g14
        CHECK(value_t(gr0(Rbo)) == doctest::Approx(std::cosh(a / 2)));
        CHECK(gr2(Rbo).vx == doctest::Approx(std::sinh(a / 2)));
        CHECK(value_t(gr0(rev(Rbo) * Rbo)) == doctest::Approx(1.0));
        CHECK(nrm_sq(gr2(rev(Rbo) * Rbo)) == doctest::Approx(0.0));

        fmt::println("exp: spatial -> cos/sin, boost -> cosh/sinh; rotor/boost are unit");
    }

    TEST_CASE("G<1,3,0>: rotor logarithm log() (inverse of exp, simple rotors)")
    {
        fmt::println("G<1,3,0>: rotor logarithm log() (inverse of exp, simple rotors)");

        // log() is the inverse of exp() for simple rotors: log(exp(B)) == B and
        // exp(log(R)) == R, for both a spatial rotation (spacelike plane, via atan2) and
        // a Lorentz boost (timelike plane, via acosh).

        // a) rotations: spatial plane g12 (B^2 < 0), angles in (0, pi)
        for (value_t const a : {0.3, 0.7, 1.2, 2.0, 3.0}) {
            auto const B = a * g12_4ds;
            auto const R = exp(B);
            auto const B_back = log(R);
            CHECK(value_t(nrm(B_back - B)) == doctest::Approx(0.0)); // recovers generator
            auto const R_back = exp(B_back);
            CHECK(value_t(gr0(R_back)) == doctest::Approx(value_t(gr0(R))));
            CHECK(value_t(nrm(gr2(R_back) - gr2(R))) == doctest::Approx(0.0));
        }

        // b) boosts: timelike plane g14 (B^2 > 0), rapidity a
        for (value_t const a : {0.3, 0.7, 1.2, 2.0}) {
            auto const B = a * g14_4ds;
            auto const R = exp(B);
            auto const B_back = log(R);
            CHECK(value_t(nrm(B_back - B)) == doctest::Approx(0.0));
            auto const R_back = exp(B_back);
            CHECK(value_t(gr0(R_back)) == doctest::Approx(value_t(gr0(R))));
            CHECK(value_t(nrm(gr2(R_back) - gr2(R))) == doctest::Approx(0.0));
        }

        // c) a tilted (multi-component) spatial plane, and the identity rotor
        {
            auto const B = 0.9 * normalize(bivec4ds{0.0, 0.0, 0.0, 1.0, -2.0, 2.0});
            CHECK(value_t(nrm(log(exp(B)) - B)) == doctest::Approx(0.0));
            auto const I = exp(bivec4ds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}); // identity rotor
            CHECK(value_t(nrm(log(I))) == doctest::Approx(0.0)); // -> zero generator
        }

        fmt::println("log: rotation via atan2, boost via acosh; exp(log(R)) == R");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 2: transform() sandwich R*X*rev(R) (rotations & Lorentz boosts)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: transform() sandwich (rotation / boost / invariance)")
    {
        fmt::println("G<1,3,0>: transform() sandwich (rotation / boost / invariance)");

        value_t const th = 0.6;  // rotation angle
        value_t const phi = 0.5; // boost rapidity
        auto const Rr = get_rotor(g12_4ds, th);
        auto const Rb = get_boost(g14_4ds, phi);

        // ---- interval invariance: Lorentz transforms preserve nrm_sq (all grades) ----
        vec4ds x{2.0, 3.0, 5.0, 7.0};
        bivec4ds B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds t{1.0, 2.0, 3.0, 4.0};
        CHECK(nrm_sq(transform(x, Rr)) == doctest::Approx(nrm_sq(x)));
        CHECK(nrm_sq(transform(x, Rb)) == doctest::Approx(nrm_sq(x)));
        CHECK(nrm_sq(transform(B, Rb)) == doctest::Approx(nrm_sq(B)));
        CHECK(nrm_sq(transform(t, Rb)) == doctest::Approx(nrm_sq(t)));

        // ---- spatial rotation in the g1-g2 plane: g1 -> cos(th) g1 + sin(th) g2 ----
        // (identical to a 3D rotation; the time direction g4 is untouched)
        auto const vr = transform(g1_4ds, Rr);
        CHECK(vr.x == doctest::Approx(std::cos(th)));
        CHECK(vr.y == doctest::Approx(std::sin(th)));
        CHECK(vr.z == doctest::Approx(0.0));
        CHECK(vr.w == doctest::Approx(0.0));

        // ---- Lorentz boost of the rest 4-velocity g4 in the g14 plane ----
        // g4 -> cosh(phi) g4 + sinh(phi) g1 ; gamma = cosh(phi), beta = tanh(phi)
        auto const u = transform(g4_4ds, Rb);
        value_t const beta = std::tanh(phi);
        CHECK(u.w == doctest::Approx(std::cosh(phi))); // gamma (time comp.)
        CHECK(u.x == doctest::Approx(std::sinh(phi))); // beta*gamma (space comp.)
        CHECK(u.x / u.w == doctest::Approx(beta));     // beta = tanh(phi)
        CHECK(u.w ==
              doctest::Approx(1.0 / std::sqrt(1.0 - beta * beta))); // gamma identity

        // ---- collinear boosts compose: rapidities add (velocities via tanh) ----
        value_t const p1 = 0.3, p2 = 0.4;
        auto const Rc = get_boost(g14_4ds, p2) * get_boost(g14_4ds, p1); // p1 then p2
        CHECK(transform(g4_4ds, Rc).w == doctest::Approx(std::cosh(p1 + p2)));
        CHECK(nrm_sq(transform(g4_4ds, Rc) -
                     transform(g4_4ds, get_boost(g14_4ds, p1 + p2))) ==
              doctest::Approx(0.0));

        // ---- optimized closed-form transform_opt() matches the direct sandwich ----
        // (validated for a general rotor = boost composed with a spatial rotation,
        //  exercising all 8 even-grade rotor coefficients, on a generic vector)
        auto const Rgen = Rb * Rr; // mixes boost + rotation -> all c0..c7 populated
        CHECK(nrm_sq(transform_opt(x, Rgen) - transform(x, Rgen)) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(transform_opt(g1_4ds, Rgen) - transform(g1_4ds, Rgen)) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(transform_opt(g4_4ds, Rgen) - transform(g4_4ds, Rgen)) ==
              doctest::Approx(0.0));
        // also preserves the interval, like the direct form
        CHECK(nrm_sq(transform_opt(x, Rgen)) == doctest::Approx(nrm_sq(x)));

        // ---- transform_opt for bivectors and trivectors matches the direct form ----
        bivec4ds Bgen{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds Tgen{1.0, 2.0, 3.0, 4.0};
        CHECK(nrm_sq(transform_opt(Bgen, Rgen) - transform(Bgen, Rgen)) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(transform_opt(Tgen, Rgen) - transform(Tgen, Rgen)) ==
              doctest::Approx(0.0));
        // each basis bivector / trivector too (exercises every matrix column)
        CHECK(nrm_sq(transform_opt(g14_4ds, Rgen) - transform(g14_4ds, Rgen)) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(transform_opt(g23_4ds, Rgen) - transform(g23_4ds, Rgen)) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(transform_opt(g123_4ds, Rgen) - transform(g123_4ds, Rgen)) ==
              doctest::Approx(0.0));
        // batch overloads agree with the scalar form
        std::vector<vec4ds> vs{x, g1_4ds, g4_4ds};
        std::vector<bivec4ds> Bs{Bgen, g14_4ds, g23_4ds};
        std::vector<trivec4ds> Ts{Tgen, g123_4ds};
        auto vo = transform_opt(vs, Rgen);
        auto Bo = transform_opt(Bs, Rgen);
        auto To = transform_opt(Ts, Rgen);
        CHECK(nrm_sq(vo[0] - transform(x, Rgen)) == doctest::Approx(0.0));
        CHECK(nrm_sq(Bo[0] - transform(Bgen, Rgen)) == doctest::Approx(0.0));
        CHECK(nrm_sq(To[0] - transform(Tgen, Rgen)) == doctest::Approx(0.0));

        fmt::println("transform: interval-invariant; rotation==3D; boost gamma=cosh, "
                     "beta=tanh; collinear boosts add rapidity; transform_opt==transform "
                     "(vec/bivec/trivec, scalar + batch)");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 3: spacetime split of a vector (time + relative space)
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // convention gate: a positive parameter must produce a positive-sense transform,
    // and it must mean the SAME thing as in ega3d. Both builders use the POSITIVE
    // half-argument, exp(+p/2 B_hat) -- the sign flip that ega3d's get_rotor gets from
    // its explicit -theta/2 is supplied here by the metric instead (spatial vectors
    // square to -1), so copying ega3d's literal sign would reverse the rotation.
    // Guarded because nothing else in the suite compares the two algebras' senses:
    // angle() and rapidity() are magnitude-only and cannot see a reversed rotation.
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: rotation and boost sense (convention gate vs. ega3d)")
    {
        fmt::println("G<1,3,0>: rotation and boost sense (convention gate vs. ega3d)");

        value_t const th = 0.6;  // rotation angle
        value_t const phi = 0.5; // boost rapidity

        // ---- both builders are the POSITIVE half-argument exponential ----
        // (is_close compares component-wise: nrm_sq would be a broken witness here,
        //  since a NULL difference also has nrm_sq == 0)
        CHECK(is_close(get_rotor(g12_4ds, th), exp(0.5 * th * g12_4ds)));
        CHECK(is_close(get_boost(g14_4ds, phi), exp(0.5 * phi * g14_4ds)));

        // ---- rotation turns g1 TOWARD g2, as ega3d turns e1 toward e2 ----
        // (ga_ega3d_test.hpp pins rotate(e1_3d, get_rotor(e12_3d, phi)) == (cos, sin, 0))
        auto const Rr = get_rotor(g12_4ds, th);
        auto const g1r = transform(g1_4ds, Rr);
        CHECK(g1r.x == doctest::Approx(std::cos(th)));
        CHECK(g1r.y == doctest::Approx(std::sin(th)));
        // and g2 goes to -sin(th) g1 + cos(th) g2 -- the same sense, one plane
        auto const g2r = transform(g2_4ds, Rr);
        CHECK(g2r.x == doctest::Approx(-std::sin(th)));
        CHECK(g2r.y == doctest::Approx(std::cos(th)));
        // the orthogonal directions are untouched
        CHECK(is_close(transform(g3_4ds, Rr), g3_4ds));
        CHECK(is_close(transform(g4_4ds, Rr), g4_4ds));

        // ---- boost carries g4 TOWARD +g1, the first factor of the plane g14 ----
        auto const u = transform(g4_4ds, get_boost(g14_4ds, phi));
        CHECK(u.x == doctest::Approx(std::sinh(phi)));
        CHECK(u.w == doctest::Approx(std::cosh(phi)));

        // ---- a negative parameter reverses each, exactly ----
        CHECK(is_close(transform(g1_4ds, get_rotor(g12_4ds, -th)),
                       vec4ds{std::cos(th), -std::sin(th), 0.0, 0.0}));
        CHECK(transform(g4_4ds, get_boost(g14_4ds, -phi)).x ==
              doctest::Approx(-std::sinh(phi)));
    }

    TEST_CASE("G<1,3,0>: spacetime split of a vector")
    {
        fmt::println("G<1,3,0>: spacetime split of a vector");

        // ---- standard frame u = g4: time = g4 component, space = g_k4 part ----
        vec4ds x{2.0, 3.0, 5.0, 7.0};
        CHECK(value_t(time_split(x, g4_4ds)) == doctest::Approx(x.w)); // time = g4 comp.
        auto const s = space_split(x, g4_4ds);                         // x ^ g4
        CHECK(s.vx == doctest::Approx(x.x));                           // g14
        CHECK(s.vy == doctest::Approx(x.y));                           // g24
        CHECK(s.vz == doctest::Approx(x.z));                           // g34
        CHECK(s.mx == doctest::Approx(0.0)); // no g23/g31/g12 part in the g4 frame
        CHECK(s.my == doctest::Approx(0.0));
        CHECK(s.mz == doctest::Approx(0.0));
        // reconstruction: x = gr1( (time + space) * g4 )
        CHECK(nrm_sq(
                  gr1(mvec4ds_e{time_split(x, g4_4ds), space_split(x, g4_4ds)} * g4_4ds) -
                  x) == doctest::Approx(0.0));

        // ---- moving observer: time dilation between frames ----
        // observer 4-velocity u = g4 boosted by rapidity phi (relative speed beta=tanh)
        value_t const phi = 0.5;
        auto const u = transform(g4_4ds, get_boost(g14_4ds, phi));
        // the rest 4-velocity g4 has time component gamma = cosh(phi) in u's frame
        CHECK(value_t(time_split(g4_4ds, u)) == doctest::Approx(std::cosh(phi)));
        // and is seen to move at speed beta = |relative vector| / time = tanh(phi)
        CHECK(value_t(nrm(space_split(g4_4ds, u))) / value_t(time_split(g4_4ds, u)) ==
              doctest::Approx(std::tanh(phi)));

        fmt::println("spacetime split: g4-frame reconstruct; moving observer -> "
                     "gamma=cosh(phi), beta=tanh(phi)");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 3b: spacetime split of a bivector (relative E / B parts)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: spacetime split of a bivector (relative E / B)")
    {
        fmt::println("G<1,3,0>: spacetime split of a bivector (relative E / B)");

        // ---- g4 frame: g_k4 -> relative-vector (E) part ; g_jk -> rel-bivector (B) ----
        bivec4ds F{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // (g14,g24,g34, g23,g31,g12)
        auto const E = rel_vec_split(F, g4_4ds);
        auto const Bp = rel_bivec_split(F, g4_4ds);
        CHECK(E == bivec4ds{1.0, 2.0, 3.0, 0.0, 0.0, 0.0});  // g_k4 part only
        CHECK(Bp == bivec4ds{0.0, 0.0, 0.0, 4.0, 5.0, 6.0}); // g_jk part only
        CHECK(E + Bp == F);                                  // parts recombine

        // ---- EM field transformation: boost a pure magnetic field B_z = g12 ----
        // along x (the g14 plane). Expect an induced electric field E_y = -gamma*beta
        // B_z and B_z' = gamma B_z (gamma = cosh phi, gamma*beta = sinh phi).
        value_t const phi = 0.5;
        auto const Fb = transform(g12_4ds, get_boost(g14_4ds, phi));
        auto const Eb = rel_vec_split(Fb, g4_4ds);
        auto const Bb = rel_bivec_split(Fb, g4_4ds);
        CHECK(Eb.vy == doctest::Approx(-std::sinh(phi))); // induced E_y (g24 component)
        CHECK(Bb.mz == doctest::Approx(std::cosh(phi)));  // B_z' = gamma (g12 component)
        CHECK(nrm_sq(Eb + Bb - Fb) == doctest::Approx(0.0)); // recombine after boost

        fmt::println("bivector split: g4-frame E/B; boost mixes B_z -> E_y "
                     "(EM field transformation)");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // the Faraday bivector and the relativistic particle (electrodynamics tier)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: the Faraday bivector -- E/B convention pinned by the force law")
    {
        fmt::println("G<1,3,0>: the Faraday bivector -- E/B convention pinned by the "
                     "force law");

        // THE CONVENTION, stated once and pinned here -- nothing else in the library
        // fixes it, and it cannot be read off the transformation law (see step 4):
        //
        //     F = Ex g14 + Ey g24 + Ez g34 - (Bx g23 + By g31 + Bz g12)
        //
        // The electric field enters the time-space slots with a PLUS, the magnetic
        // field the space-space slots with a MINUS. Its companion is that
        // transform(X, get_boost(B, phi)) is an ACTIVE boost (it carries g4 -> +g1 for
        // g14, see the header), so a PASSIVE change to a frame moving at +v uses -phi.
        //
        // The PAIR is what matters: flipping either half alone breaks the physics, and
        // the field transformation law alone cannot see the difference. Natural units,
        // c = 1. Griffiths, Introduction to Electrodynamics 4th ed., Ch. 12.3.
        auto faraday = [](double Ex, double Ey, double Ez, double Bx, double By,
                          double Bz) { return bivec4ds{Ex, Ey, Ez, -Bx, -By, -Bz}; };
        auto four_vel = [](double vx, double vy, double vz) {
            double const g = 1.0 / std::sqrt(1.0 - (vx * vx + vy * vy + vz * vz));
            return vec4ds{g * vx, g * vy, g * vz, g};
        };

        // ---- 1. the boost sense the convention is matched to: ACTIVE ---------------
        value_t const phi = 0.55;
        auto const ub = transform(g4_4ds, get_boost(g14_4ds, phi));
        CHECK(ub.x == doctest::Approx(std::sinh(phi)));
        CHECK(ub.w == doctest::Approx(std::cosh(phi)));

        // ---- 2. the hand-checkable case, worked on paper before the library --------
        // B along +z, motion along +x:  q v x B = q (v g1) x (Bz g3) = -q v Bz g2
        {
            double const Bz = 1.3, v = 0.4;
            double const g = 1.0 / std::sqrt(1.0 - v * v);
            auto const K = gr1(faraday(0.0, 0.0, 0.0, 0.0, 0.0, Bz) * four_vel(v, 0, 0));
            CHECK(K.x == doctest::Approx(0.0).epsilon(1.0e-12));
            CHECK(K.y == doctest::Approx(-g * v * Bz));
            CHECK(K.z == doctest::Approx(0.0).epsilon(1.0e-12));
            // magnetic forces do no work (Griffiths 5.11): the time component vanishes
            CHECK(K.w == doctest::Approx(0.0).epsilon(1.0e-12));
        }

        // ---- 3. the general Lorentz force, K = q gr1(F u)  (Griffiths 12.128) ------
        double const Ex = 0.3, Ey = -1.2, Ez = 0.8, Bx = 1.1, By = 0.4, Bz = -0.6;
        double const vx = 0.2, vy = -0.35, vz = 0.1;
        double const gv = 1.0 / std::sqrt(1.0 - (vx * vx + vy * vy + vz * vz));
        auto const F = faraday(Ex, Ey, Ez, Bx, By, Bz);
        auto const u = four_vel(vx, vy, vz);
        auto const K = gr1(F * u);
        // classical reference v x B, written out independently
        double const cx = vy * Bz - vz * By;
        double const cy = vz * Bx - vx * Bz;
        double const cz = vx * By - vy * Bx;
        CHECK(K.x == doctest::Approx(gv * (Ex + cx))); // Griffiths (12.129)
        CHECK(K.y == doctest::Approx(gv * (Ey + cy)));
        CHECK(K.z == doctest::Approx(gv * (Ez + cz)));
        CHECK(K.w == doctest::Approx(gv * (Ex * vx + Ey * vy + Ez * vz))); // the power
        // F . u is orthogonal to u, so u . u = 1 is preserved by the motion
        CHECK(value_t(dot(K, u)) == doctest::Approx(0.0).epsilon(1.0e-12));

        // ---- 4. all SIX field transformation rules, from ONE boost (12.109) --------
        double const v = 0.6;
        double const gam = 1.0 / std::sqrt(1.0 - v * v), ph = std::atanh(v);
        auto const Fb = transform(F, get_boost(g14_4ds, -ph)); // passive, frame at +v
        CHECK(Fb.vx == doctest::Approx(Ex));                   // Ex unchanged
        CHECK(Fb.vy == doctest::Approx(gam * (Ey - v * Bz)));
        CHECK(Fb.vz == doctest::Approx(gam * (Ez + v * By)));
        CHECK(-Fb.mx == doctest::Approx(Bx)); // Bx unchanged
        CHECK(-Fb.my == doctest::Approx(gam * (By + v * Ez)));
        CHECK(-Fb.mz == doctest::Approx(gam * (Bz - v * Ey)));

        // ---- 5. THE TRAP: the transformation law cannot see the magnetic sign ------
        // the WRONG embedding (+B) combined with the opposite boost satisfies the same
        // six rules -- the two sign errors cancel. Only the force law separates them,
        // which is why both are asserted in one test case.
        auto const Fw = bivec4ds{Ex, Ey, Ez, Bx, By, Bz}; // +B: wrong
        auto const Fwb = transform(Fw, get_boost(g14_4ds, ph));
        CHECK(Fwb.vy == doctest::Approx(gam * (Ey - v * Bz))); // still passes
        CHECK(Fwb.mz == doctest::Approx(gam * (Bz - v * Ey))); // still passes
        auto const Kw = gr1(Fw * u);
        CHECK(Kw.y != doctest::Approx(gv * (Ey + cy))); // only this fires ...
        CHECK(Kw.y == doctest::Approx(gv * (Ey - cy))); // ... it is E - v x B

        // ---- 6. the observer split reads the two halves back out ------------------
        CHECK(rel_vec_split(F, g4_4ds) == bivec4ds{Ex, Ey, Ez, 0.0, 0.0, 0.0});
        CHECK(rel_bivec_split(F, g4_4ds) == bivec4ds{0.0, 0.0, 0.0, -Bx, -By, -Bz});

        fmt::println("Faraday bivector: E on g_k4 (+), B on g_jk (-); the force law "
                     "fixes the sign the boost law cannot see");
    }

    TEST_CASE("G<1,3,0>: the rotor equation of motion (cyclotron, hyperbolic, E x B)")
    {
        fmt::println("G<1,3,0>: the rotor equation of motion (cyclotron, hyperbolic, "
                     "E x B)");

        // The relativistic charge in an electromagnetic field, carried by a rotor:
        //
        //     dR/dtau = (q/2m) F R ,   u = transform(g4, R) ,   dx/dtau = u
        //
        // For constant F this integrates in closed form to R(tau) = exp((q tau/2m) F) R0
        // -- ONE exp of a general bivector (a boost times a rotation). RK4 (the shared
        // rk4_step) is checked against that closed form, and both against the textbook
        // answers. Same Faraday convention and natural units as the test above.
        auto faraday = [](double Ex, double Ey, double Ez, double Bx, double By,
                          double Bz) { return bivec4ds{Ex, Ey, Ez, -Bx, -By, -Bz}; };
        // the VERSOR norm -- NOT nrm_sq, which is the grade-wise Lorentzian sum and is
        // not 1 for a unit rotor (see the sqrt(rotor) note in the ops header)
        auto versor_nrm = [](mvec4ds_e const& R) {
            return std::sqrt(std::abs(value_t(gr0(rev(R) * R))));
        };

        // proper-time RK4 on the 12-component state (rotor: 8, world line: 4)
        auto integrate = [&](bivec4ds const& F, double qo2m, mvec4ds_e R0, vec4ds x0,
                             double dtau, int nsteps) {
            std::array<value_t, 12> u_mem{}, rhs_mem{};
            std::array<value_t, 24> uh_mem{};
            auto pack = [&](mvec4ds_e const& R, vec4ds const& x) {
                u_mem = {R.c0, R.c1, R.c2, R.c3, R.c4, R.c5,
                         R.c6, R.c7, x.x,  x.y,  x.z,  x.w};
            };
            auto unpack_R = [&]() {
                return mvec4ds_e{u_mem[0], u_mem[1], u_mem[2], u_mem[3],
                                 u_mem[4], u_mem[5], u_mem[6], u_mem[7]};
            };
            auto unpack_x = [&]() {
                return vec4ds{u_mem[8], u_mem[9], u_mem[10], u_mem[11]};
            };
            pack(R0, x0);
            auto uv = mdspan<value_t, dextents<size_t, 1>>(u_mem.data(), 12);
            auto uh = mdspan<value_t, dextents<size_t, 2>>(uh_mem.data(), 2, 12);
            auto rhs = mdspan<value_t const, dextents<size_t, 1>>(rhs_mem.data(), 12);
            for (int n = 0; n < nsteps; ++n) {
                for (size_t rk = 1; rk <= 4; ++rk) {
                    auto const R = unpack_R();
                    auto const dR = qo2m * (F * R); // bivector * even -> even
                    auto const uu = transform(g4_4ds, R);
                    rhs_mem = {dR.c0, dR.c1, dR.c2, dR.c3, dR.c4, dR.c5,
                               dR.c6, dR.c7, uu.x,  uu.y,  uu.z,  uu.w};
                    rk4_step(uv, uh, rhs, dtau, rk);
                }
                auto R = unpack_R();
                R = R * (1.0 / versor_nrm(R)); // renormalise the rotor
                pack(R, unpack_x());
            }
            return std::pair{unpack_R(), unpack_x()};
        };

        // nrm_sq is NULL-CAPABLE in a Minkowski metric: it can vanish (or go negative)
        // for a non-zero operand, so every "is this zero?" here is COMPONENTWISE -- the
        // same rule the source-free calculus case states
        auto max_abs = [](mvec4ds const& M) {
            value_t const* c = &M.c0;
            value_t m = 0.0;
            for (int i = 0; i < 16; ++i)
                m = std::max(m, std::abs(c[i]));
            return m;
        };

        auto const R_id = mvec4ds_e{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        auto const x_0 = vec4ds{0.0, 0.0, 0.0, 0.0};

        // ---- gate 1: cyclotron motion in a constant magnetic field ----------------
        {
            double const Bz = 0.8, qo2m = 0.5; // q/m = 1
            double const v0 = 0.5, g0 = 1.0 / std::sqrt(1.0 - v0 * v0);
            double const om_p = 2.0 * qo2m * Bz; // proper-time turn rate q B / m
            auto const F = faraday(0.0, 0.0, 0.0, 0.0, 0.0, Bz);
            auto const R0 = get_boost(g14_4ds, std::atanh(v0)); // start along +x

            double const tau = 0.9;
            auto const [R, x] = integrate(F, qo2m, R0, x_0, tau / 2000, 2000);

            // RK4 reproduces the closed form R(tau) = exp((q tau / 2m) F) R0
            auto const Rc = exp((qo2m * tau) * F) * R0;
            CHECK(value_t(gr0(R)) == doctest::Approx(value_t(gr0(Rc))));
            CHECK(max_abs(mvec4ds(gr2(R) - gr2(Rc))) < 1.0e-8);

            // gamma is constant: a magnetic field does no work (Griffiths 5.11)
            auto const uu = transform(g4_4ds, R);
            CHECK(uu.w == doctest::Approx(g0));
            CHECK(value_t(dot(uu, uu)) == doctest::Approx(1.0));

            // the velocity turns CLOCKWISE in the xy plane (q > 0, B along +z) at the
            // proper rate q B / m -- so in coordinate time at q B / (gamma m). The sign
            // is asserted, not its magnitude: it is the second gate on the B convention
            CHECK(std::atan2(uu.y, uu.x) == doctest::Approx(-om_p * tau));

            // half a proper period turns the velocity around and displaces the particle
            // by one orbit diameter: r = gamma m v / (q B), i.e. Griffiths (5.3) p = QBR
            double const tau_half = pi / om_p;
            auto const [Rh, xh] = integrate(F, qo2m, R0, x_0, tau_half / 4000, 4000);
            double const r_orbit = g0 * v0 / (2.0 * qo2m * Bz);
            CHECK(xh.x == doctest::Approx(0.0).epsilon(1.0e-8));
            CHECK(xh.y == doctest::Approx(-2.0 * r_orbit));
            CHECK(xh.z == doctest::Approx(0.0).epsilon(1.0e-12));
            CHECK(xh.w == doctest::Approx(g0 * tau_half)); // dt = gamma dtau
        }

        // ---- gate 2: hyperbolic motion in a constant electric field ---------------
        {
            double const Ex = 0.7, qo2m = 0.5;    // q/m = 1
            double const alpha = 2.0 * qo2m * Ex; // proper acceleration q E / m
            auto const F = faraday(Ex, 0.0, 0.0, 0.0, 0.0, 0.0);
            double const tau = 1.3;
            auto const [R, x] = integrate(F, qo2m, R_id, x_0, tau / 2000, 2000);
            auto const uu = transform(g4_4ds, R);

            // the world line is a boost with rapidity linear in proper time
            CHECK(uu.w == doctest::Approx(std::cosh(alpha * tau)));        // gamma
            CHECK(uu.x == doctest::Approx(std::sinh(alpha * tau)));        // gamma v
            CHECK(x.w == doctest::Approx(std::sinh(alpha * tau) / alpha)); // t(tau)

            // Griffiths (12.60)/(12.61) in COORDINATE time -- independent statements,
            // since t here came from integrating dx/dtau alongside the rotor
            double const t = x.w, at = alpha * t;
            CHECK(uu.x / uu.w == doctest::Approx(at / std::sqrt(1.0 + at * at)));
            CHECK(x.x == doctest::Approx((std::sqrt(1.0 + at * at) - 1.0) / alpha));
        }

        // ---- gate 3: the E x B drift ---------------------------------------------
        {
            double const Ey = 0.35, Bz = 0.9, qo2m = 0.5;
            double const vd = Ey / Bz; // Griffiths (5.10): the drift speed is E/B
            double const gd = 1.0 / std::sqrt(1.0 - vd * vd);
            auto const F = faraday(0.0, Ey, 0.0, 0.0, 0.0, Bz);

            // (a) boosting into the drift frame kills the electric field entirely and
            // leaves a pure magnetic field, weaker by 1/gamma_d
            auto const Fd = transform(F, get_boost(g14_4ds, -std::atanh(vd)));
            // (the E-part lies entirely in the g_k4 slots, whose metric signs agree,
            // so here nrm_sq IS definite -- checked componentwise anyway)
            CHECK(max_abs(mvec4ds(rel_vec_split(Fd, g4_4ds))) < 1.0e-12);
            CHECK(-rel_bivec_split(Fd, g4_4ds).mz == doctest::Approx(Bz / gd));

            // (b) |E| < |B| makes F rotation-like, so the rotor closes after an exact
            // proper period tau_c = 2 pi m / (q sqrt(B^2 - E^2)). Over one period the
            // cycloid returns to its guiding centre and only the drift remains.
            CHECK(is_simple(F));
            double const tau_c = 2.0 * pi / (2.0 * qo2m * std::sqrt(Bz * Bz - Ey * Ey));
            auto const [R, x] = integrate(F, qo2m, R_id, x_0, tau_c / 4000, 4000);
            CHECK(x.y == doctest::Approx(0.0).epsilon(1.0e-6));
            CHECK(x.z == doctest::Approx(0.0).epsilon(1.0e-12));
            CHECK(x.x / x.w == doctest::Approx(vd).epsilon(1.0e-6)); // Griffiths (5.10)
        }

        fmt::println("rotor EOM: cyclotron orbit, hyperbolic motion and the E x B drift");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 4: projections and rejections (onto vector / bivector)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: projections and rejections")
    {
        fmt::println("G<1,3,0>: projections and rejections");

        vec4ds v{2.0, 3.0, 5.0, 7.0};

        // ---- onto a timelike axis g4: projection = time part, rejection = space part
        // ----
        CHECK(project_onto(v, g4_4ds) == vec4ds{0.0, 0.0, 0.0, 7.0});
        CHECK(reject_from(v, g4_4ds) == vec4ds{2.0, 3.0, 5.0, 0.0});
        CHECK(project_onto(v, g4_4ds) + reject_from(v, g4_4ds) == v);

        // ---- onto a spacelike axis g1 ----
        CHECK(project_onto(v, g1_4ds) == vec4ds{2.0, 0.0, 0.0, 0.0});
        CHECK(project_onto(v, g1_4ds) + reject_from(v, g1_4ds) == v);
        // the projection is parallel to the target (wedge vanishes)
        CHECK(nrm_sq(wdg(project_onto(v, g1_4ds), g1_4ds)) == doctest::Approx(0.0));

        // ---- onto a bivector (the spatial plane g12) ----
        auto const p = project_onto(v, g12_4ds);
        auto const r = reject_from(v, g12_4ds);
        CHECK(nrm_sq(p + r - v) == doctest::Approx(0.0));       // project + reject == v
        CHECK(nrm_sq(wdg(p, g12_4ds)) == doctest::Approx(0.0)); // projection lies in g12
        CHECK(p.x == doctest::Approx(2.0));                     // g1 part
        CHECK(p.y == doctest::Approx(3.0));                     // g2 part
        CHECK(p.z == doctest::Approx(0.0));
        CHECK(p.w == doctest::Approx(0.0));

        fmt::println("project + reject == original; projection lies in / parallel to "
                     "the (non-null) target");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 5: reflections (gpr-based sandwich, Macdonald p.129)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: reflections (involution, nrm_sq, discrete symmetries)")
    {
        fmt::println("G<1,3,0>: reflections (involution, nrm_sq, discrete symmetries)");

        vec4ds v{2.0, 3.0, 5.0, 7.0}; // (g1, g2, g3, g4=time)

        // ---- reflect on a hyperplane via its timelike normal g4 == time reflection T
        // ---- the component along g4 (time) flips, the spatial part is kept
        CHECK(reflect_on(v, g4_4ds) == vec4ds{2.0, 3.0, 5.0, -7.0});
        // reflect on the same hyperplane given directly as the spatial trivector g123:
        // reflection is invariant under scaling/orientation of the mirror, and g123 is
        // dual to the normal g4 -> must give the identical result as the normal form
        CHECK(nrm_sq(reflect_on(v, g123_4ds) - reflect_on(v, g4_4ds)) ==
              doctest::Approx(0.0));

        // ---- reflect onto the time line g4 == parity P (spatial inversion) ----
        // the time component is kept, the three spatial components flip
        CHECK(reflect_on_vec(v, g4_4ds) == vec4ds{-2.0, -3.0, -5.0, 7.0});

        // ---- every reflection preserves the spacetime interval nrm_sq ----
        CHECK(nrm_sq(reflect_on(v, g4_4ds)) == doctest::Approx(nrm_sq(v)));
        CHECK(nrm_sq(reflect_on(v, g1_4ds)) == doctest::Approx(nrm_sq(v)));
        CHECK(nrm_sq(reflect_on(v, g123_4ds)) == doctest::Approx(nrm_sq(v)));
        CHECK(nrm_sq(reflect_on(v, g14_4ds)) == doctest::Approx(nrm_sq(v)));
        CHECK(nrm_sq(reflect_on_vec(v, g4_4ds)) == doctest::Approx(nrm_sq(v)));

        // ---- every reflection is an involution: applying it twice is the identity ----
        CHECK(nrm_sq(reflect_on(reflect_on(v, g4_4ds), g4_4ds) - v) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(reflect_on(reflect_on(v, g14_4ds), g14_4ds) - v) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(reflect_on(reflect_on(v, g123_4ds), g123_4ds) - v) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(reflect_on_vec(reflect_on_vec(v, g4_4ds), g4_4ds) - v) ==
              doctest::Approx(0.0));

        // ---- bivector reflection in a bivector is also an nrm_sq-preserving involution
        bivec4ds B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        CHECK(nrm_sq(reflect_on(B, g14_4ds)) == doctest::Approx(nrm_sq(B)));
        CHECK(nrm_sq(reflect_on(reflect_on(B, g14_4ds), g14_4ds) - B) ==
              doctest::Approx(0.0));

        // ---- two reflections compose to a proper rotation (a rotor sandwich) ----
        // reflecting on g1 then on a vector b at spatial half-angle a rotates by 2a in
        // the g1-g2 plane; here check the composite preserves nrm_sq and equals a rotor
        // transform (the classic "rotation = two reflections" theorem)
        value_t const a = 0.35;
        vec4ds b{std::cos(a), std::sin(a), 0.0, 0.0}; // unit spatial direction at angle a
        auto const twice = reflect_on_vec(reflect_on_vec(v, g1_4ds), b);
        CHECK(nrm_sq(twice) == doctest::Approx(nrm_sq(v)));
        // composite of the two spatial reflections is a rotation by 2a in the g12
        // plane, in the sense from the first mirror toward the second (g1 -> b), which
        // is get_rotor's positive sense (is_close, not nrm_sq: a NULL difference would
        // also have nrm_sq == 0)
        auto const Rrot = get_rotor(g12_4ds, 2.0 * a);
        CHECK(is_close(twice, transform(v, Rrot)));

        fmt::println("reflections: T / P symmetries; interval-preserving involutions; "
                     "two reflections == one rotation");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 6: angle (spacelike) / rapidity (timelike)
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: angle (spacelike) and rapidity (timelike)")
    {
        fmt::println("G<1,3,0>: angle (spacelike) and rapidity (timelike)");

        using std::numbers::pi;

        // ---- Euclidean angle between spacelike vectors ----
        CHECK(angle(g1_4ds, g1_4ds) == doctest::Approx(0.0));      // self -> 0
        CHECK(angle(g1_4ds, g2_4ds) == doctest::Approx(pi / 2.0)); // orthogonal axes
        CHECK(angle(g1_4ds, g3_4ds) == doctest::Approx(pi / 2.0));
        // angle equals the rotation parameter: rotate g1 by theta in the g12 plane
        value_t const th = 0.7;
        auto const g1_rot = transform(g1_4ds, get_rotor(g12_4ds, th));
        CHECK(is_spacelike(g1_rot)); // rotation keeps a spacelike vector spacelike
        CHECK(angle(g1_4ds, g1_rot) == doctest::Approx(th));
        // symmetric and scale-independent
        CHECK(angle(g2_4ds, g1_4ds) == doctest::Approx(pi / 2.0));
        CHECK(angle(2.0 * g1_4ds, 3.0 * g2_4ds) == doctest::Approx(pi / 2.0));

        // ---- relative rapidity between timelike 4-velocities ----
        CHECK(rapidity(g4_4ds, g4_4ds) == doctest::Approx(0.0)); // same rest frame -> 0
        // boost g4 by rapidity phi in the g14 plane; rapidity back to g4 must be phi
        value_t const phi = 0.9;
        auto const u = transform(g4_4ds, get_boost(g14_4ds, phi));
        CHECK(is_timelike(u)); // a boosted 4-velocity stays timelike
        CHECK(rapidity(g4_4ds, u) == doctest::Approx(phi));
        // gamma = cosh(rapidity) = dot(g4, u) for unit 4-velocities
        CHECK(value_t(dot(g4_4ds, u)) == doctest::Approx(std::cosh(phi)));
        // collinear boosts add rapidity (consistency with the step-2 finding)
        value_t const p1 = 0.4, p2 = 0.5;
        auto const u2 = transform(g4_4ds, get_boost(g14_4ds, p1 + p2));
        CHECK(rapidity(g4_4ds, u2) == doctest::Approx(p1 + p2));

        // ---- causal-character guards: wrong domain throws ----
        CHECK_THROWS(angle(g4_4ds, g4_4ds));    // timelike into spacelike angle()
        CHECK_THROWS(angle(g1_4ds, g4_4ds));    // mixed
        CHECK_THROWS(rapidity(g1_4ds, g1_4ds)); // spacelike into rapidity()
        CHECK_THROWS(rapidity(g1_4ds, g4_4ds)); // mixed

        fmt::println(
            "angle: spacelike Euclidean angle == rotation param; "
            "rapidity: timelike hyperbolic angle == boost param; domain-guarded");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp: sqrt(rotor) -- halves the rotation angle / boost rapidity
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: sqrt(rotor) (half-angle / half-rapidity)")
    {
        fmt::println("G<1,3,0>: sqrt(rotor) (half-angle / half-rapidity)");

        // max abs component difference between two even rotors (convention-agnostic)
        auto rotor_diff = [](mvec4ds_e const& a, mvec4ds_e const& b) {
            return std::max({std::abs(a.c0 - b.c0), std::abs(a.c1 - b.c1),
                             std::abs(a.c2 - b.c2), std::abs(a.c3 - b.c3),
                             std::abs(a.c4 - b.c4), std::abs(a.c5 - b.c5),
                             std::abs(a.c6 - b.c6), std::abs(a.c7 - b.c7)});
        };

        value_t const th = 0.8;  // rotation angle
        value_t const phi = 0.9; // boost rapidity

        // ---- defining property: sqrt(R) * sqrt(R) == R  (rotation and boost) ----
        auto const Rr = get_rotor(g12_4ds, th);
        auto const Rb = get_boost(g14_4ds, phi);
        CHECK(rotor_diff(sqrt(Rr) * sqrt(Rr), Rr) == doctest::Approx(0.0));
        CHECK(rotor_diff(sqrt(Rb) * sqrt(Rb), Rb) == doctest::Approx(0.0));

        // ---- sqrt halves the parameter: sqrt(R(x)) == R(x/2) ----
        CHECK(rotor_diff(sqrt(Rr), get_rotor(g12_4ds, th / 2)) == doctest::Approx(0.0));
        CHECK(rotor_diff(sqrt(Rb), get_boost(g14_4ds, phi / 2)) == doctest::Approx(0.0));

        // ---- result is a unit rotor (versor norm gr0(rev(s) s) == 1) ----
        auto const sr = sqrt(Rr);
        auto const sb = sqrt(Rb);
        CHECK(value_t(gr0(rev(sr) * sr)) == doctest::Approx(1.0));
        CHECK(value_t(gr0(rev(sb) * sb)) == doctest::Approx(1.0));

        // ---- geometric meaning: applying sqrt(R) twice == applying R once ----
        vec4ds x{2.0, 3.0, 5.0, 7.0};
        CHECK(nrm_sq(transform(transform(x, sr), sr) - transform(x, Rr)) ==
              doctest::Approx(0.0));
        CHECK(nrm_sq(transform(transform(x, sb), sb) - transform(x, Rb)) ==
              doctest::Approx(0.0));

        // ---- identity rotor: sqrt(1) == 1 ----
        mvec4ds_e id{scalar4ds{1.0}, bivec4ds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                     pscalar4ds{0.0}};
        CHECK(rotor_diff(sqrt(id), id) == doctest::Approx(0.0));

        // ---- degenerate: sqrt(-1) (a 2*pi rotation rotor) -> identity, not a crash ----
        mvec4ds_e neg_id{scalar4ds{-1.0}, bivec4ds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                         pscalar4ds{0.0}};
        CHECK(rotor_diff(sqrt(neg_id), id) == doctest::Approx(0.0));

        fmt::println("sqrt: sqrt(R)^2 == R; halves angle/rapidity; unit rotor; "
                     "two half-transforms == one; degenerate -1 -> identity");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp: exp / log / sqrt for GENERAL (non-simple) rotors
    // -- a Lorentz boost and a spatial rotation combined in orthogonal (dual) planes,
    //    so R carries a g1234 part (gr4(R) != 0). This is the case the earlier
    //    simple-only implementation could not handle; guards against regressions in the
    //    invariant-decomposition (exp/log) and Study-renormalisation (sqrt) code paths.
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: general (non-simple) rotor exp / log / sqrt")
    {
        fmt::println("G<1,3,0>: general (non-simple) rotor exp / log / sqrt");

        auto rotor_diff = [](mvec4ds_e const& a, mvec4ds_e const& b) {
            return std::max({std::abs(a.c0 - b.c0), std::abs(a.c1 - b.c1),
                             std::abs(a.c2 - b.c2), std::abs(a.c3 - b.c3),
                             std::abs(a.c4 - b.c4), std::abs(a.c5 - b.c5),
                             std::abs(a.c6 - b.c6), std::abs(a.c7 - b.c7)});
        };

        // reference exp via truncated Taylor series (ground truth for moderate |B|)
        auto exp_series = [](bivec4ds const& B, int N) {
            mvec4ds_e term{scalar4ds{1.0}, bivec4ds{0, 0, 0, 0, 0, 0}, pscalar4ds{0.0}};
            mvec4ds_e sum = term;
            mvec4ds_e const Bmv{scalar4ds{0.0}, B, pscalar4ds{0.0}};
            for (int k = 1; k <= N; ++k) {
                term = term * Bmv;
                term = term * (1.0 / k);
                sum = sum + term;
            }
            return sum;
        };

        // ---- a concrete non-simple bivector: boost in g14 + rotation in g23 (dual
        //      planes) plus off-axis components so no coordinate is accidentally zero
        //      ----
        {
            bivec4ds const B{0.5, -0.2, 0.1, 0.4, -0.3, 0.25}; // v-part boost, m-part rot
            REQUIRE(std::abs(value_t(gr4(mvec4ds_e{B} * B))) >
                    1.0e-3); // genuinely non-simple

            auto const R = exp(B);
            CHECK(std::abs(value_t(gr4(R))) > 1.0e-3); // rotor has a g1234 part

            // exp matches the series definition
            CHECK(rotor_diff(R, exp_series(B, 30)) == doctest::Approx(0.0));

            // log is the inverse of exp (both directions)
            auto const B_back = log(R);
            CHECK(value_t(nrm(B_back - B)) == doctest::Approx(0.0));
            CHECK(rotor_diff(exp(B_back), R) == doctest::Approx(0.0));

            // sqrt: defining property and unit-rotor / geometric meaning
            auto const sR = sqrt(R);
            CHECK(rotor_diff(sR * sR, R) == doctest::Approx(0.0));
            CHECK(value_t(gr0(rev(sR) * sR)) == doctest::Approx(1.0)); // unit rotor
            vec4ds const x{2.0, 3.0, 5.0, 7.0};
            CHECK(nrm_sq(transform(transform(x, sR), sR) - transform(x, R)) ==
                  doctest::Approx(0.0)); // two half-transforms == one full
        }

        // ---- edge: rotation angle == pi/2, where the boost plane vanishes from gr2(R)
        //      and must be recovered from the pseudoscalar part via the dual ----
        {
            bivec4ds const B = 0.8 * g14_4ds + (pi / 2.0) * g23_4ds;
            auto const R = exp(B);
            CHECK(value_t(nrm(log(R) - B)) == doctest::Approx(0.0));
            CHECK(rotor_diff(sqrt(R) * sqrt(R), R) == doctest::Approx(0.0));
        }

        // ---- larger boost + rotation, still non-simple (outside series range) ----
        {
            bivec4ds const B = 1.3 * g14_4ds + 1.1 * g31_4ds;
            auto const R = exp(B);
            CHECK(value_t(nrm(log(R) - B)) == doctest::Approx(0.0));
            CHECK(rotor_diff(sqrt(R) * sqrt(R), R) == doctest::Approx(0.0));
            CHECK(value_t(gr0(rev(sqrt(R)) * sqrt(R))) == doctest::Approx(1.0));
        }

        // ---- deterministic random sweep of non-simple rotors ----
        {
            std::mt19937 gen(20240617u);
            std::uniform_real_distribution<value_t> d(-1.0, 1.0);
            int tested = 0;
            value_t worst_series = 0.0, worst_logexp = 0.0, worst_explog = 0.0;
            value_t worst_sqrt = 0.0, worst_unit = 0.0;
            for (int i = 0; i < 2000; ++i) {
                bivec4ds B{d(gen), d(gen), d(gen), d(gen), d(gen), d(gen)};
                B = 0.35 * B; // moderate norm so the series converges
                if (std::abs(value_t(gr4(mvec4ds_e{B} * B))) < 1.0e-3) continue; // simple
                ++tested;
                auto const R = exp(B);
                worst_series = std::max(worst_series, rotor_diff(R, exp_series(B, 30)));
                auto const Bb = log(R);
                worst_logexp = std::max(worst_logexp, value_t(nrm(Bb - B)));
                worst_explog = std::max(worst_explog, rotor_diff(exp(Bb), R));
                auto const sR = sqrt(R);
                worst_sqrt = std::max(worst_sqrt, rotor_diff(sR * sR, R));
                worst_unit =
                    std::max(worst_unit, std::abs(value_t(gr0(rev(sR) * sR)) - 1.0));
            }
            REQUIRE(tested > 1500); // most random bivectors are non-simple
            CHECK(worst_series == doctest::Approx(0.0).epsilon(1.0e-9));
            CHECK(worst_logexp == doctest::Approx(0.0).epsilon(1.0e-9));
            CHECK(worst_explog == doctest::Approx(0.0).epsilon(1.0e-9));
            CHECK(worst_sqrt == doctest::Approx(0.0).epsilon(1.0e-9));
            CHECK(worst_unit == doctest::Approx(0.0).epsilon(1.0e-9));
            fmt::println("   random non-simple rotors tested: {}", tested);
            fmt::println("   worst: exp-series {:.2e}  log(exp)-B {:.2e}  exp(log)-R "
                         "{:.2e}  sqrt^2-R {:.2e}  unit {:.2e}",
                         worst_series, worst_logexp, worst_explog, worst_sqrt,
                         worst_unit);
        }

        fmt::println("exp/log/sqrt: non-simple rotors (boost x rotation in dual planes) "
                     "round-trip to machine precision");
    }

    TEST_CASE("G<1,3,0>: left-right complement composition")
    {
        fmt::println("G<1,3,0>: left-right complement composition");

        // For even-dimensional algebras: l_cmpl(r_cmpl(u)) = r_cmpl(l_cmpl(u)) = u
        scalar4ds s{3.0};
        vec4ds v{2.0, 5.0, 7.0, 9.0};
        bivec4ds B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds T{1.0, 2.0, 3.0, 4.0};
        pscalar4ds ps{11.0};
        mvec4ds_e M_e{scalar4ds{1.0}, bivec4ds{2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
                      pscalar4ds{8.0}};
        mvec4ds M{scalar4ds{1.0}, vec4ds{2.0, 3.0, 4.0, 5.0},
                  bivec4ds{6.0, 7.0, 8.0, 9.0, 10.0, 11.0},
                  trivec4ds{12.0, 13.0, 14.0, 15.0}, pscalar4ds{16.0}};

        CHECK(l_cmpl(r_cmpl(s)) == s);
        CHECK(l_cmpl(r_cmpl(v)) == v);
        CHECK(l_cmpl(r_cmpl(B)) == B);
        CHECK(l_cmpl(r_cmpl(T)) == T);
        CHECK(l_cmpl(r_cmpl(ps)) == ps);
        CHECK(l_cmpl(r_cmpl(M_e)) == M_e);
        CHECK(l_cmpl(r_cmpl(M)) == M);

        // Also test r_cmpl(l_cmpl(u)) = u
        CHECK(r_cmpl(l_cmpl(s)) == s);
        CHECK(r_cmpl(l_cmpl(v)) == v);
        CHECK(r_cmpl(l_cmpl(B)) == B);
        CHECK(r_cmpl(l_cmpl(T)) == T);
        CHECK(r_cmpl(l_cmpl(ps)) == ps);
        CHECK(r_cmpl(l_cmpl(M_e)) == M_e);
        CHECK(r_cmpl(l_cmpl(M)) == M);

        fmt::println("l_cmpl(r_cmpl(u)) = r_cmpl(l_cmpl(u)) = u for even-dimensional");
    }

    TEST_CASE("G<1,3,0>: MVec4ds_E and MVec4ds_U formatting tests")
    {
        fmt::println("G<3,0,0>: MVec4ds_E and MVec4ds_U formatting tests");

        // Test MVec4ds_E formatting (MVec8_t with 8 components: c0 through c7)
        MVec4ds_E<double> mvec4ds_e_val{10.1, 20.2, 30.3, 40.4, 50.5, 60.6, 70.7, 80.8};

        // Basic output
        std::string mvec4ds_e_basic = fmt::format("{}", mvec4ds_e_val);
        CHECK(mvec4ds_e_basic == "MVec4ds_E(10.1,20.2,30.3,40.4,50.5,60.6,70.7,80.8)");

        // Two decimal places
        std::string mvec4ds_e_two_decimals = fmt::format("{:.2f}", mvec4ds_e_val);
        CHECK(mvec4ds_e_two_decimals ==
              "MVec4ds_E(10.10,20.20,30.30,40.40,50.50,60.60,70.70,80.80)");

        // Scientific notation
        std::string mvec4ds_e_scientific = fmt::format("{:.2e}", mvec4ds_e_val);
        CHECK(mvec4ds_e_scientific == "MVec4ds_E(1.01e+01,2.02e+01,3.03e+01,4.04e+01,"
                                      "5.05e+01,6.06e+01,7.07e+01,8.08e+01)");

        // Contextual usage
        std::string mvec4ds_e_contextual =
            fmt::format("Even 4ds MV: {:.1f}", mvec4ds_e_val);
        CHECK(mvec4ds_e_contextual ==
              "Even 4ds MV: MVec4ds_E(10.1,20.2,30.3,40.4,50.5,60.6,70.7,80.8)");

        // Test MVec4ds_U formatting (MVec8_t with 8 components: c0 through c7)
        MVec4ds_U<double> mvec4ds_u_val{11.1, 22.2, 33.3, 44.4, 55.5, 66.6, 77.7, 88.8};

        // Basic output
        std::string mvec4ds_u_basic = fmt::format("{}", mvec4ds_u_val);
        CHECK(mvec4ds_u_basic == "MVec4ds_U(11.1,22.2,33.3,44.4,55.5,66.6,77.7,88.8)");

        // Two decimal places
        std::string mvec4ds_u_two_decimals = fmt::format("{:.2f}", mvec4ds_u_val);
        CHECK(mvec4ds_u_two_decimals ==
              "MVec4ds_U(11.10,22.20,33.30,44.40,55.50,66.60,77.70,88.80)");

        // Scientific notation
        std::string mvec4ds_u_scientific = fmt::format("{:.2e}", mvec4ds_u_val);
        CHECK(mvec4ds_u_scientific == "MVec4ds_U(1.11e+01,2.22e+01,3.33e+01,4.44e+01,"
                                      "5.55e+01,6.66e+01,7.77e+01,8.88e+01)");

        // Contextual usage
        std::string mvec4ds_u_contextual =
            fmt::format("Odd 4ds MV: {:.1f}", mvec4ds_u_val);
        CHECK(mvec4ds_u_contextual ==
              "Odd 4ds MV: MVec4ds_U(11.1,22.2,33.3,44.4,55.5,66.6,77.7,88.8)");

        fmt::println("   MVec4ds_E basic: {}", mvec4ds_e_basic);
        fmt::println("   MVec4ds_E 2-decimal: {}", mvec4ds_e_two_decimals);
        fmt::println("   MVec4ds_E scientific: {}", mvec4ds_e_scientific);
        fmt::println("   MVec4ds_E contextual: {}", mvec4ds_e_contextual);

        fmt::println("   MVec4ds_U basic: {}", mvec4ds_u_basic);
        fmt::println("   MVec4ds_U 2-decimal: {}", mvec4ds_u_two_decimals);
        fmt::println("   MVec4ds_U scientific: {}", mvec4ds_u_scientific);
        fmt::println("   MVec4ds_U contextual: {}", mvec4ds_u_contextual);
    }

    TEST_CASE("G<1,3,0>: defining basic types and ctor checks")
    {
        fmt::println("G<1,3,0>: defining basic types and ctor checks");

        auto s1 = scalar4ds{5.0};
        auto v1 = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto b1 = bivec4ds{-1.0, 2.0, -3.0, 4.0, 5.0, -6.0};
        auto t1 = trivec4ds{-1.0, -2.0, -3.0, -4.0};
        auto ps1 = pscalar4ds{-5.0};

        auto mv1 = mvec4ds{s1};
        auto mv2 = mvec4ds{v1};
        auto mv3 = mvec4ds{b1};
        auto mv4 = mvec4ds{t1};
        auto mv5 = mvec4ds{ps1};

        auto mv6a = mvec4ds{s1, b1, ps1};
        auto mv6b = mvec4ds_e{s1};
        auto mv6c = mvec4ds_e{b1};
        auto mv6d = mvec4ds_e{ps1};
        auto mv6e = mvec4ds_e{s1, b1, ps1};
        auto mv6f = mvec4ds{mv6e};

        auto mv7a = mvec4ds{v1, t1};
        auto mv7b = mvec4ds_u{v1};
        auto mv7c = mvec4ds_u{t1};
        auto mv7d = mvec4ds_u{v1, t1};
        auto mv7e = mvec4ds{mv7d};

        auto mv8 = mvec4ds{s1, v1, b1, t1, ps1};

        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv4  = {}", mv4);
        // fmt::println("   mv5  = {}", mv5);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6b = {}", mv6b);
        // fmt::println("   mv6c = {}", mv6c);
        // fmt::println("   mv6d = {}", mv6d);
        // fmt::println("   mv6e = {}", mv6e);
        // fmt::println("   mv6f = {}", mv6f);
        // fmt::println("   mv7a = {}", mv7a);
        // fmt::println("   mv7b = {}", mv7b);
        // fmt::println("   mv7c = {}", mv7c);
        // fmt::println("   mv7d = {}", mv7d);
        // fmt::println("   mv7e = {}", mv7e);
        // fmt::println("   mv8  = {}", mv8);

        CHECK(gr0(mv1) == s1);
        CHECK(gr1(mv2) == v1);
        CHECK(gr2(mv3) == b1);
        CHECK(gr3(mv4) == t1);
        CHECK(gr4(mv5) == ps1);

        CHECK(gr0(mv6a) == s1);
        CHECK(gr1(mv6a) == vec4ds{});
        CHECK(gr2(mv6a) == b1);
        CHECK(gr3(mv6a) == trivec4ds{});
        CHECK(gr4(mv6a) == ps1);

        CHECK(gr0(mv6b) == s1);
        CHECK(gr2(mv6b) == bivec4ds{});
        CHECK(gr4(mv6b) == pscalar4ds{});

        CHECK(gr0(mv6c) == scalar4ds{});
        CHECK(gr2(mv6c) == b1);
        CHECK(gr4(mv6c) == pscalar4ds{});

        CHECK(gr0(mv6d) == scalar4ds{});
        CHECK(gr2(mv6d) == bivec4ds{});
        CHECK(gr4(mv6d) == ps1);

        CHECK(gr0(mv6e) == s1);
        CHECK(gr2(mv6e) == b1);
        CHECK(gr4(mv6e) == ps1);

        CHECK(gr0(mv6f) == s1);
        CHECK(gr1(mv6f) == vec4ds{});
        CHECK(gr2(mv6f) == b1);
        CHECK(gr3(mv6f) == trivec4ds{});
        CHECK(gr4(mv6f) == ps1);

        CHECK(gr0(mv7a) == scalar4ds{});
        CHECK(gr1(mv7a) == v1);
        CHECK(gr2(mv7a) == bivec4ds{});
        CHECK(gr3(mv7a) == t1);
        CHECK(gr4(mv7a) == pscalar4ds{});

        CHECK(gr1(mv7b) == v1);
        CHECK(gr3(mv7b) == trivec4ds{});

        CHECK(gr1(mv7c) == vec4ds{});
        CHECK(gr3(mv7c) == t1);

        CHECK(gr1(mv7d) == v1);
        CHECK(gr3(mv7d) == t1);

        CHECK(gr0(mv7e) == scalar4ds{});
        CHECK(gr1(mv7e) == v1);
        CHECK(gr2(mv7e) == bivec4ds{});
        CHECK(gr3(mv7e) == t1);
        CHECK(gr4(mv7e) == pscalar4ds{});

        CHECK(gr0(mv8) == s1);
        CHECK(gr1(mv8) == v1);
        CHECK(gr2(mv8) == b1);
        CHECK(gr3(mv8) == t1);
        CHECK(gr4(mv8) == ps1);
    }

    TEST_CASE("G<1,3,0>: Scalar4ds and PScalar4ds formatting tests")
    {
        fmt::println("G<1,3,0>: Scalar4ds and PScalar4ds formatting tests");

        // Test Scalar4ds formatting
        Scalar4ds<double> scalar_val{2.23607};

        // Basic output
        std::string scalar_basic = fmt::format("{}", scalar_val);
        CHECK(scalar_basic == "Scalar4ds(2.23607)");

        // Two decimal places
        std::string scalar_two_decimals = fmt::format("{:.2f}", scalar_val);
        CHECK(scalar_two_decimals == "Scalar4ds(2.24)");

        // Scientific notation
        std::string scalar_scientific = fmt::format("{:.2e}", scalar_val);
        CHECK(scalar_scientific == "Scalar4ds(2.24e+00)");

        // Contextual usage
        std::string scalar_contextual = fmt::format("Norm: {:.3f}", scalar_val);
        CHECK(scalar_contextual == "Norm: Scalar4ds(2.236)");

        // Test PScalar4ds formatting
        PScalar4ds<double> pscalar_val{2.44949};

        // Basic output
        std::string pscalar_basic = fmt::format("{}", pscalar_val);
        CHECK(pscalar_basic == "PScalar4ds(2.44949)");

        // Two decimal places
        std::string pscalar_two_decimals = fmt::format("{:.2f}", pscalar_val);
        CHECK(pscalar_two_decimals == "PScalar4ds(2.45)");

        // Scientific notation
        std::string pscalar_scientific = fmt::format("{:.2e}", pscalar_val);
        CHECK(pscalar_scientific == "PScalar4ds(2.45e+00)");

        // Contextual usage
        std::string pscalar_contextual = fmt::format("4ds Volume: {:.3f}", pscalar_val);
        CHECK(pscalar_contextual == "4ds Volume: PScalar4ds(2.449)");

        fmt::println("   Scalar4ds basic: {}", scalar_basic);
        fmt::println("   Scalar4ds 2-decimal: {}", scalar_two_decimals);
        fmt::println("   Scalar4ds scientific: {}", scalar_scientific);
        fmt::println("   Scalar4ds contextual: {}", scalar_contextual);

        fmt::println("   PScalar4ds basic: {}", pscalar_basic);
        fmt::println("   PScalar4ds 2-decimal: {}", pscalar_two_decimals);
        fmt::println("   PScalar4ds scientific: {}", pscalar_scientific);
        fmt::println("   PScalar4ds contextual: {}", pscalar_contextual);
    }

    TEST_CASE("G<1,3,0>: Vec4ds and TriVec4ds formatting tests")
    {
        fmt::println("G<1,3,0>: Vec4ds and TriVec4ds formatting tests");

        // Test Vec4ds formatting
        Vec4ds<double> vec4ds_val{3.14159, 2.71828, 1.41421, 1.73205};

        // Basic output
        std::string basic = fmt::format("{}", vec4ds_val);
        CHECK(basic == "Vec4ds(3.14159,2.71828,1.41421,1.73205)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", vec4ds_val);
        CHECK(two_decimals == "Vec4ds(3.14,2.72,1.41,1.73)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", vec4ds_val);
        CHECK(scientific == "Vec4ds(3.14e+00,2.72e+00,1.41e+00,1.73e+00)");

        // Contextual usage
        std::string contextual = fmt::format("Position: {:.2f}", vec4ds_val);
        CHECK(contextual == "Position: Vec4ds(3.14,2.72,1.41,1.73)");

        // Test TriVec4ds formatting
        TriVec4ds<double> trivec4ds_val{1.2, 2.3, 3.4, 4.5};

        // Basic output
        std::string trivec_basic = fmt::format("{}", trivec4ds_val);
        CHECK(trivec_basic == "TriVec4ds(1.2,2.3,3.4,4.5)");

        // Two decimal places
        std::string trivec_two_decimals = fmt::format("{:.2f}", trivec4ds_val);
        CHECK(trivec_two_decimals == "TriVec4ds(1.20,2.30,3.40,4.50)");

        // Scientific notation
        std::string trivec_scientific = fmt::format("{:.2e}", trivec4ds_val);
        CHECK(trivec_scientific == "TriVec4ds(1.20e+00,2.30e+00,3.40e+00,4.50e+00)");

        // Contextual usage
        std::string trivec_contextual = fmt::format("Volume: {:.2f}", trivec4ds_val);
        CHECK(trivec_contextual == "Volume: TriVec4ds(1.20,2.30,3.40,4.50)");

        fmt::println("   Vec4ds basic: {}", basic);
        fmt::println("   Vec4ds 2-decimal: {}", two_decimals);
        fmt::println("   Vec4ds scientific: {}", scientific);
        fmt::println("   Vec4ds contextual: {}", contextual);

        fmt::println("   TriVec4ds basic: {}", trivec_basic);
        fmt::println("   TriVec4ds 2-decimal: {}", trivec_two_decimals);
        fmt::println("   TriVec4ds scientific: {}", trivec_scientific);
        fmt::println("   TriVec4ds contextual: {}", trivec_contextual);
    }

    TEST_CASE("G<1,3,0>: BiVec4ds formatting tests")
    {
        fmt::println("G<1,3,0>: BiVec4ds formatting tests");

        // Test BiVec4ds formatting (BVec6_t with 6 components: vx, vy, vz, mx, my, mz)
        BiVec4ds<double> bivec4ds_val{1.1, 2.2, 3.3, 4.4, 5.5, 6.6};

        // Basic output
        std::string basic = fmt::format("{}", bivec4ds_val);
        CHECK(basic == "BiVec4ds(1.1,2.2,3.3,4.4,5.5,6.6)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", bivec4ds_val);
        CHECK(two_decimals == "BiVec4ds(1.10,2.20,3.30,4.40,5.50,6.60)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", bivec4ds_val);
        CHECK(scientific ==
              "BiVec4ds(1.10e+00,2.20e+00,3.30e+00,4.40e+00,5.50e+00,6.60e+00)");

        // Contextual usage
        std::string contextual = fmt::format("Plane: {:.2f}", bivec4ds_val);
        CHECK(contextual == "Plane: BiVec4ds(1.10,2.20,3.30,4.40,5.50,6.60)");

        fmt::println("   BiVec4ds basic: {}", basic);
        fmt::println("   BiVec4ds 2-decimal: {}", two_decimals);
        fmt::println("   BiVec4ds scientific: {}", scientific);
        fmt::println("   BiVec4ds contextual: {}", contextual);
    }

    TEST_CASE("G<1,3,0>: MVec4ds formatting tests")
    {
        fmt::println("G<1,3,0>: MVec4ds formatting tests");

        // Test MVec4ds formatting (MVec16_t with 16 components: c0 through c15)
        MVec4ds<double> mvec4ds_val{1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,
                                    9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

        // Basic output
        std::string basic = fmt::format("{}", mvec4ds_val);
        CHECK(basic == "MVec4ds(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", mvec4ds_val);
        CHECK(two_decimals == "MVec4ds(1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,"
                              "9.00,10.00,11.00,12.00,13.00,14.00,15.00,16.00)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", mvec4ds_val);
        CHECK(scientific == "MVec4ds(1.00e+00,2.00e+00,3.00e+00,4.00e+00,5.00e+00,"
                            "6.00e+00,7.00e+00,8.00e+00,9.00e+00,1.00e+01,1.10e+01,"
                            "1.20e+01,1.30e+01,1.40e+01,1.50e+01,1.60e+01)");

        // Contextual usage
        std::string contextual = fmt::format("Full multivector: {:.1f}", mvec4ds_val);
        CHECK(contextual == "Full multivector: MVec4ds(1.0,2.0,3.0,4.0,5.0,6.0,7.0,"
                            "8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0)");

        fmt::println("   MVec4ds basic: {}", basic);
        fmt::println("   MVec4ds 2-decimal: {}", two_decimals);
        fmt::println("   MVec4ds scientific: {}", scientific);
        fmt::println("   MVec4ds contextual: {}", contextual);
    }

    TEST_CASE("MVec4ds: metric / antimetric exomorphisms (G, rG) - comparison table")
    {
        fmt::println(
            "MVec4ds: metric / antimetric exomorphisms (G, rG) - comparison table");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Exomorphism comparison table for sta4ds = G(1,3,0), det = -1 (Minkowski).
        // Signature m = (-1,-1,-1,+1): g1,g2,g3 spacelike, g4 timelike.
        //
        // Same TWO exomorphisms as ega/pga (Lengyel, "PGA Illuminated"):
        //   - the metric      G : a WEDGE exomorphism
        //         G(wdg(a,b))   == wdg(G(a),G(b))
        //   - the antimetric rG : an ANTIWEDGE exomorphism
        //         rG(rwdg(a,b)) == rwdg(rG(a),rG(b))
        // related by  G * rG = det(metric) * I. Non-degenerate (det = -1 != 0), so unlike
        // PGA there is no bulk/weight split; the question is purely whether the STORED
        // metric is multiplicative.
        //
        // P, Q, g_S defined exactly as in the ega/pga tables:
        //   P(e_S) = prod_{i in S} m_i   (pure product of the vector squares)
        //   Q(e_S) = sigma(k) * P        (the blade square e_S . e_S),
        //            sigma(k) = (-1)^(k(k-1)/2)
        // g_S = the value STORED in sta4ds_metric (ga_usr_consts.hpp), read straight off
        // the diagonal. Under P-unify sta4ds stores g_S == P, exactly like ega/pga (the
        // is_minkowski override that used to store Q = sigma(k)*P was removed). P is
        // multiplicative, so BOTH verdicts are YES. (Q -- which differs from P at grades
        // 2,3 -- now lives only in the geometric square gr0(X*X), used for causal
        // character, NOT in the metric.)
        ////////////////////////////////////////////////////////////////////////////////

        int const m[4] = {-1, -1, -1, 1}; // g1^2, g2^2, g3^2, g4^2 (g4 timelike)
        char const* nm[16] = {"1",    "g1",   "g2",   "g3",   "g4",  "g14",
                              "g24",  "g34",  "g23",  "g31",  "g12", "g234",
                              "g314", "g124", "g123", "g1234"};
        int const msk[16] = {0, 1, 2, 4, 8, 9, 10, 12, 6, 5, 3, 14, 13, 11, 7, 15};
        int const full = 0b1111;
        auto const grade = [](int mask) {
            int k = 0;
            for (int i = 0; i < 4; ++i)
                k += (mask >> i) & 1;
            return k;
        };
        auto const Pof = [&](int mask) {
            int p = 1;
            for (int i = 0; i < 4; ++i)
                if (mask & (1 << i)) p *= m[i];
            return p;
        };
        auto const sigma = [](int k) { return ((k * (k - 1) / 2) & 1) ? -1 : 1; };
        auto const idx_of = [&](int mask) {
            for (int i = 0; i < 16; ++i)
                if (msk[i] == mask) return i;
            return 0;
        };
        auto const Gv = sta4ds_metric_view();
        auto const storedM = [&](int mask) {
            int const i = idx_of(mask);
            return Gv[i, i];
        };

        int const det = Pof(full);

        fmt::println("   blade |  k |  P |  Q | g_S | rG_S | g_S*rG_S(=det={:+})", det);
        fmt::println("   ------+----+----+----+-----+------+-------------------");
        for (int i = 0; i < 16; ++i) {
            int const k = grade(msk[i]);
            int const P = Pof(msk[i]);
            int const Q = sigma(k) * P;
            int const g = storedM(msk[i]);
            int const rg = storedM(msk[i] ^ full); // rG_S = g_{S^c}
            fmt::println(
                "   {:>5} | {:>2} | {:>+2} | {:>+2} | {:>+3} | {:>+4} | {:>+12}  {}",
                nm[i], k, P, Q, g, rg, g * rg, (g == P ? "" : "[g_S != P!]"));
        }
        fmt::println("");

        // WEDGE verdict: g multiplicative iff g_{S u T} == g_S * g_T for disjoint S,T.
        bool wedge_exo = true;
        for (int a2 = 0; a2 <= full; ++a2)
            for (int b2 = 0; b2 <= full; ++b2)
                if ((a2 & b2) == 0 && storedM(a2 | b2) != storedM(a2) * storedM(b2))
                    wedge_exo = false;
        // ANTIWEDGE verdict: rwdg(e_S,e_T) != 0 iff S u T == full; result e_{S n T}.
        // rG_S = g_{S^c} is an antiwedge exomorphism iff rG_{S n T} == rG_S * rG_T there.
        auto const rGof = [&](int mask) { return storedM(mask ^ full); };
        bool antiwedge_exo = true;
        for (int a2 = 0; a2 <= full; ++a2)
            for (int b2 = 0; b2 <= full; ++b2)
                if ((a2 | b2) == full && rGof(a2 & b2) != rGof(a2) * rGof(b2))
                    antiwedge_exo = false;
        fmt::println("   stored metric g    multiplicative (WEDGE     exomorphism): {}",
                     wedge_exo ? "YES" : "NO");
        fmt::println("   antimetric rG=g_Sc multiplicative (ANTIWEDGE exomorphism): {}",
                     antiwedge_exo ? "YES" : "NO");
        fmt::println("");

        CHECK(det == -1);
        // sta4ds stores the metric exomorphism P, exactly like ega and pga -- the
        // extended metric is the wedge exomorphism seeded by the vector signature, and
        // it is a different object from the blade square: g_S == P at every blade, so
        // BOTH verdicts are YES -- the stored metric is a wedge exomorphism and its
        // antimetric rG = g_{S^c} is an antiwedge exomorphism. (The blade square Q, which
        // is NOT multiplicative, was moved out of the metric slot into the geometric
        // square gr0(X*X), where it carries causal character / rotor structure.)
        CHECK(wedge_exo);
        CHECK(antiwedge_exo);
        // every diagonal entry now equals P, and g_S * rG_S == det uniformly
        // (P_S * P_{S^c} = det) -- no grade-dependent sigma(k)sigma(n-k) alternation
        // (that alternation was the Q fingerprint).
        for (int i = 0; i < 16; ++i) {
            CHECK(storedM(msk[i]) == Pof(msk[i]));
            CHECK(storedM(msk[i]) * storedM(msk[i] ^ full) == det);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // Library handles (same even-dim non-degenerate form as ega2d):
        //   Gx  = l_cmpl o r_dual         -- the stored metric (reproduces sta4ds_metric)
        //   rGx = r_cmpl o Gx o l_cmpl    -- the antimetric; on a blade e_S this is
        //                                    g_{S^c} * e_S
        ////////////////////////////////////////////////////////////////////////////////
        auto const Gx = [](auto const& X) { return l_cmpl(r_dual(X)); };
        auto const rGx = [&Gx](auto const& X) { return r_cmpl(Gx(l_cmpl(X))); };

        // concrete witnesses (mixed-grade blade pairs): both exomorphism laws now HOLD
        CHECK(Gx(wdg(g1_4ds, g4_4ds)) == wdg(Gx(g1_4ds), Gx(g4_4ds)));
        CHECK(rGx(rwdg(g1_4ds, g234_4ds)) == rwdg(rGx(g1_4ds), rGx(g234_4ds)));
    }

    TEST_CASE("MVec4ds: metric / antimetric exomorphisms (G, rG) and metric-indep. dual")
    {
        fmt::println("MVec4ds: metric / antimetric exomorphisms (G, rG)");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Following Lengyel ("PGA Illuminated"), every metric algebra carries TWO
        // exomorphisms tied to the metric:
        //   - the metric      G : a wedge      exomorphism,  G(wdg(a,b))  ==
        //   wdg(G(a),G(b))
        //   - the antimetric rG : an antiwedge exomorphism, rG(rwdg(a,b)) ==
        //   rwdg(rG(a),rG(b))
        // related by  G * rG = det(metric) * I  (ordinary matrix product, I = identity).
        //
        // In PGA the metric is degenerate (det = 0), so G*rG = 0 -- G and rG annihilate
        // each other, which is exactly the bulk/weight split (disjoint support). In a
        // NON-degenerate algebra like sta4ds there is no null direction, det != 0, so the
        // relation is invertible. With det = -1 and G^2 = I (every basis blade squares to
        // +/-1) it collapses to
        //
        //                          rG = -G          (the antimetric is the negated
        //                          metric)
        //
        // Under P-unify the LIBRARY metric IS the faithful exomorphism (no reversion
        // fix):
        //     G(X) = l_cmpl(r_dual(X))     (since r_dual = r_cmpl o G, l_cmpl o r_cmpl =
        //     id)
        // gives G(e_S) = P(e_S) e_S = (prod_{i in S} g_i) e_S -- the pure product P,
        // which IS multiplicative, so G is a genuine wedge exomorphism at every grade.
        // (Under the old Q-metric, l_cmpl o r_dual gave the blade square sigma(k)*P and
        // one needed rev() to recover P; that reversion was the gr_inv that broke l_dual
        // o r_dual. It is gone now -- the blade square Q lives only in the geometric
        // square gr0(X*X).)
        //
        // Signature duality (the physics remark): G(g1)=-g1 realizes (-,-,-,+) while
        // rG = -G realizes (+,+,+,-). The two signatures are wedge-metric and
        // antiwedge-antimetric of one another -- the same geometry, read off in dual
        // spaces. One school works it out in the space, the other in its dual.
        ////////////////////////////////////////////////////////////////////////////////

        auto const sf = [](auto const& x) { return fmt::format("{}", x); };
        auto const G = [](auto const& X) {
            return l_cmpl(r_dual(X));
        }; // metric P (exomorphism)
        // antimetric: rG(e_S) = P(e_S^c) e_S = -P(e_S) e_S, i.e. exactly rG = -G (det =
        // -1)
        auto const rG = [&G](auto const& X) { return -G(X); };

        auto const s = scalar4ds{2.0};
        auto const v = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto const t = trivec4ds{1.0, 2.0, 3.0, 4.0};
        auto const ps = pscalar4ds{5.0};

        fmt::println("   metric vs antimetric on the basis vectors (signature duality):");
        fmt::println("     G(g1)={:>13}  rG(g1)={:>13}", sf(G(g1_4ds)), sf(rG(g1_4ds)));
        fmt::println("     G(g4)={:>13}  rG(g4)={:>13}", sf(G(g4_4ds)), sf(rG(g4_4ds)));
        fmt::println("");

        // (0) the library metric G is the pure product P on each basis blade
        CHECK(G(g14_4ds) == -g14_4ds); // P(g14) = m1*m4 = -1
        CHECK(G(g23_4ds) == g23_4ds);  // P(g23) = m2*m3 = +1
        // and the antimetric is the negated metric: rG = -G
        CHECK(rG(v) == -G(v));
        CHECK(rG(B) == -G(B));
        CHECK(rG(t) == -G(t));

        // (1) G * rG = det(metric) * I, det = -1  <=>  rG(G(X)) == -X
        CHECK(rG(G(s)) == -s);
        CHECK(rG(G(v)) == -v);
        CHECK(rG(G(B)) == -B);
        CHECK(rG(G(t)) == -t);
        CHECK(rG(G(ps)) == -ps);

        // signature duality: G realizes (-,-,-,+), rG = -G realizes (+,+,+,-)
        CHECK(G(g1_4ds) == -g1_4ds);
        CHECK(rG(g1_4ds) == g1_4ds);
        CHECK(G(g4_4ds) == g4_4ds);
        CHECK(rG(g4_4ds) == -g4_4ds);

        // (2) G is a WEDGE exomorphism at every grade (incl. mixed grades) -- under
        //     P-unify the library metric is multiplicative, so this now holds directly
        CHECK(G(wdg(g1_4ds, g4_4ds)) == wdg(G(g1_4ds), G(g4_4ds)));   // vector ^ vector
        CHECK(G(wdg(g2_4ds, g3_4ds)) == wdg(G(g2_4ds), G(g3_4ds)));   // vector ^ vector
        CHECK(G(wdg(g1_4ds, g23_4ds)) == wdg(G(g1_4ds), G(g23_4ds))); // vector ^ bivector
        CHECK(G(wdg(g4_4ds, g23_4ds)) == wdg(G(g4_4ds), G(g23_4ds))); // vector ^ bivector

        // (3) rG is an ANTIWEDGE exomorphism: rG(rwdg(a,b)) == rwdg(rG(a),rG(b))
        CHECK(rG(rwdg(g234_4ds, g123_4ds)) == rwdg(rG(g234_4ds), rG(g123_4ds)));
        CHECK(rG(rwdg(g234_4ds, g14_4ds)) == rwdg(rG(g234_4ds), rG(g14_4ds)));

        ////////////////////////////////////////////////////////////////////////////////
        // Dual round trip (clean, no gr_inv under P-unify).
        //
        // With the exomorphism metric, l_dual and r_dual are mutual inverses up to the
        // determinant -- a single global sign, NOT the grade-dependent gr_inv that the
        // old Q-metric forced:
        //     l_dual(r_dual(X)) == r_dual(l_dual(X)) == det * X == -X    (det = -1)
        // For a Euclidean det = +1 algebra this is the familiar l_dual(r_dual) == X; the
        // only difference here is the global -1. (The metric-on-the-dual-side r_undual
        // = l_cmpl o r_dual o l_cmpl still inverts r_dual EXACTLY, with no det factor.)
        ////////////////////////////////////////////////////////////////////////////////

        // clean round trip: l_dual o r_dual == det * id == -id (no grade involution)
        CHECK(l_dual(r_dual(s)) == -s);
        CHECK(l_dual(r_dual(v)) == -v);
        CHECK(l_dual(r_dual(B)) == -B);
        CHECK(l_dual(r_dual(t)) == -t);
        CHECK(l_dual(r_dual(ps)) == -ps);
        CHECK(r_dual(l_dual(v)) == -v);
        CHECK(r_dual(l_dual(t)) == -t);

        // exact metric-independent inverse: the library r_undual (metric on the dual
        // side)
        CHECK(r_undual(r_dual(s)) == s);
        CHECK(r_undual(r_dual(v)) == v);
        CHECK(r_undual(r_dual(B)) == B);
        CHECK(r_undual(r_dual(t)) == t);
        CHECK(r_undual(r_dual(ps)) == ps);
        CHECK(r_dual(r_undual(v)) == v);
        CHECK(r_dual(r_undual(B)) == B);
    }

    TEST_CASE("MVec4ds: r_undual / l_undual; the wedge product fixes space vs dual")
    {
        fmt::println(
            "MVec4ds: r_undual / l_undual; the wedge product fixes space vs dual");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Q1 -- there are TWO unduals, one inverse for each handedness of the dual:
        //     r_undual (inverts r_dual):  r_undual(D) = l_cmpl(r_dual(l_cmpl(D))) = G o
        //     l_cmpl l_undual (inverts l_dual):  l_undual(D) = r_cmpl(l_dual(r_cmpl(D)))
        //     = G o r_cmpl
        //   (metric on the dual side; the complement that opens/closes matches the dual's
        //   handedness). On even-grade inputs l_dual == r_dual, so r_undual == l_undual
        //   there; on odd grades they differ by a sign, exactly as l_dual and r_dual do.
        //   Each inverts only its own handedness -- the cross pairing flips the odd
        //   grades.
        //
        // Q2 -- space vs dual space is fixed by WHICH product is the wedge. Each
        // signature
        //   is an exomorphism of EXACTLY ONE product (they are not interchangeable):
        //     G  = (-,-,-,+)  is the WEDGE     exomorphism, NOT the antiwedge one
        //     rG = (+,+,+,-)  is the ANTIWEDGE exomorphism, NOT the wedge one
        //   So, relative to the library's wedge convention, the space is uniquely the
        //   signature whose metric respects the wedge product: (-,-,-,+). Its dual space
        //   is
        //   (+,+,+,-). The structure is symmetric under (wdg <-> rwdg, G <-> rG,
        //   space <-> dual) -- which is exactly why both signatures are physically
        //   equivalent: each is "the space" for one of the two products.
        //
        // Under P-unify the LIBRARY metric G = l_cmpl o r_dual IS the exomorphism P (no
        // reversion fix needed), so the Q2 laws hold directly on the library's own dual.
        ////////////////////////////////////////////////////////////////////////////////

        auto const G = [](auto const& X) {
            return l_cmpl(r_dual(X));
        }; // metric P = (-,-,-,+)
        auto const rG = [&G](auto const& X) { return -G(X); }; // antimetric (+,+,+,-)
        // r_undual / l_undual are the library functions (ga_sta4ds_ops_basics.hpp)

        auto const s = scalar4ds{2.0};
        auto const v = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto const t = trivec4ds{1.0, 2.0, 3.0, 4.0};
        auto const ps = pscalar4ds{5.0};

        // Q1: r_undual inverts r_dual at every grade
        CHECK(r_undual(r_dual(s)) == s);
        CHECK(r_undual(r_dual(v)) == v);
        CHECK(r_undual(r_dual(B)) == B);
        CHECK(r_undual(r_dual(t)) == t);
        CHECK(r_undual(r_dual(ps)) == ps);
        // l_undual inverts l_dual at every grade
        CHECK(l_undual(l_dual(s)) == s);
        CHECK(l_undual(l_dual(v)) == v);
        CHECK(l_undual(l_dual(B)) == B);
        CHECK(l_undual(l_dual(t)) == t);
        CHECK(l_undual(l_dual(ps)) == ps);
        // r_undual and l_undual are genuinely distinct: each inverts only its own
        // handedness, the cross-pairing flips the odd grades (and coincides on the even
        // grades)
        CHECK(l_undual(r_dual(v)) == -v); // l_undual does NOT invert r_dual (odd grade)
        CHECK(r_undual(l_dual(v)) == -v); // r_undual does NOT invert l_dual (odd grade)
        CHECK(l_undual(r_dual(B)) == B);  // ... but they agree on even grades
        CHECK(r_undual(l_dual(B)) == B);

        // Q2: G = (-,-,-,+) is the WEDGE exomorphism, and ONLY the wedge one
        CHECK(G(wdg(g1_4ds, g4_4ds)) == wdg(G(g1_4ds), G(g4_4ds))); // wedge: holds
        CHECK(G(wdg(g2_4ds, g3_4ds)) == wdg(G(g2_4ds), G(g3_4ds)));
        CHECK(G(rwdg(g234_4ds, g14_4ds)) !=
              rwdg(G(g234_4ds), G(g14_4ds))); // antiwedge: fails
        // rG = (+,+,+,-) is the ANTIWEDGE exomorphism, and ONLY the antiwedge one
        CHECK(rG(rwdg(g234_4ds, g14_4ds)) ==
              rwdg(rG(g234_4ds), rG(g14_4ds))); // antiwedge: holds
        CHECK(rG(rwdg(g234_4ds, g123_4ds)) == rwdg(rG(g234_4ds), rG(g123_4ds)));
        CHECK(rG(wdg(g1_4ds, g4_4ds)) != wdg(rG(g1_4ds), rG(g4_4ds))); // wedge: fails

        fmt::println(
            "   Q1: r_undual inverts r_dual, l_undual inverts l_dual (distinct on "
            "odd grades)");
        fmt::println("   Q2: wedge-exomorphism uniquely selects (-,-,-,+) as the space; "
                     "(+,+,+,-) is its dual");
        fmt::println("");
    }

    TEST_CASE("MVec4ds: explicit space -> dual -> space round trip")
    {
        fmt::println("MVec4ds: explicit space -> dual -> space round trip");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Round trip in sta4ds (non-degenerate, Minkowski). The metric must be applied on
        // OPPOSITE sides going out and coming back:
        //
        //   forward  (space -> dual):  D = r_dual(A)                = r_cmpl(G * A)
        //   backward (dual  -> space): A = r_undual(D)              =
        //   l_cmpl(r_dual(l_cmpl(D)))
        //                                                           = (1/det) * rG *
        //                                                           l_cmpl(D)
        //
        // The "naive" backward leg l_dual(D) (metric on the primal side again -- the
        // recipe that works in the Euclidean even algebra ega2d) does NOT recover A here:
        // it returns det*A = -A. Under P-unify that is a CLEAN global sign at every grade
        // (the determinant), not the grade-dependent det*gr_inv the old Q-metric produced
        // -- so dividing by det (a single -1) repairs it; no grade involution is needed.
        ////////////////////////////////////////////////////////////////////////////////

        auto const sf = [](auto const& x) { return fmt::format("{}", x); };
        // r_undual is the library function (ga_sta4ds_ops_basics.hpp)

        auto const s = scalar4ds{2.0};
        auto const v = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto const t = trivec4ds{1.0, 2.0, 3.0, 4.0};
        auto const ps = pscalar4ds{5.0};

        auto const row = [&sf](std::string_view name, auto const& A) {
            auto const D = r_dual(A);
            fmt::println("   {:<3} | {:<24} | {:<24} | {:<24} | {}", name, sf(A), sf(D),
                         sf(r_undual(D)), sf(l_dual(D)));
        };
        fmt::println("   {:<3} | {:<24} | {:<24} | {:<24} | {}", "A", "value (space)",
                     "D = r_dual(A) (dual)", "r_undual(D) [== A]",
                     "l_dual(D) [== det*A]");
        fmt::println("   {:-<108}", "");
        row("s", s);
        row("v", v);
        row("B", B);
        row("t", t);
        row("ps", ps);
        fmt::println("");

        // forward then backward recovers A at EVERY grade (r_undual = correct inverse)
        CHECK(r_undual(r_dual(s)) == s);
        CHECK(r_undual(r_dual(v)) == v);
        CHECK(r_undual(r_dual(B)) == B);
        CHECK(r_undual(r_dual(t)) == t);
        CHECK(r_undual(r_dual(ps)) == ps);

        // full multivector, both orders
        auto const M = mvec4ds{s, v, B, t, ps};
        CHECK(r_undual(r_dual(M)) == M);
        CHECK(r_dual(r_undual(M)) == M);

        // the naive backward leg l_dual yields det*A = -A at EVERY grade (clean global
        // sign, no grade-dependent involution -- nothing "survives by luck")
        CHECK(l_dual(r_dual(s)) == -s);
        CHECK(l_dual(r_dual(v)) == -v);
        CHECK(l_dual(r_dual(B)) == -B);
        CHECK(l_dual(r_dual(t)) == -t);
        CHECK(l_dual(r_dual(ps)) == -ps);

        // the backward leg expressed with the antimetric rG = -G  (uses both G and rG):
        //   r_undual(D) = (1/det) * rG * l_cmpl(D),  with det = I^2 = -1, rG(Y) = -G(Y)
        value_t const det = value_t(I_4ds * I_4ds);                      // = -1
        auto const Gx = [](auto const& X) { return l_cmpl(r_dual(X)); }; // metric
        auto const rG = [&Gx](auto const& X) { return -Gx(X); }; // antimetric = -G
        CHECK((value_t(1.0) / det) * rG(l_cmpl(r_dual(v))) == v);
        CHECK((value_t(1.0) / det) * rG(l_cmpl(r_dual(B))) == B);
    }

    TEST_CASE("MVec4ds: dualization - complement vs. pseudoscalar-multiplication")
    {
        fmt::println("MVec4ds: dualization - complement vs. pseudoscalar-multiplication");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Two schools of thought for dualizing in sta4ds = G(1,3,0). Like ega2d this is
        // an EVEN-dimensional algebra (l_dual != r_dual for odd grades), but unlike the
        // Euclidean algebras it carries a NON-trivial (Minkowski) metric, so the dual is
        //   l_dual(A) = l_cmpl(G * A),   r_dual(A) = r_cmpl(G * A)
        // (complement after multiplying with the extended metric G, det(G) = -1).
        //
        //   (1) THIS LIBRARY: l_dual / r_dual (metric duals), and the non-metric
        //       complements l_cmpl / r_cmpl.
        //
        //   (2) OTHER SCHOOL: to_dual(A) = A * I_4ds,  from_dual(A) = A * inv(I_4ds),
        //       with I_4ds*I_4ds == -1  =>  inv(I_4ds) == -I_4ds.
        //
        // Three-way contrast across the algebras (under P-unify sta4ds joins the EGA
        // pattern -- the Minkowski "surprises" of the old Q-metric are gone):
        //
        //   ega3d (odd, Euclidean):  one dual, it is an involution dual(dual(A))==A;
        //                            A*I differs from dual(A) by the reversion sign.
        //   ega2d (even, Euclidean): l_dual != r_dual; metric = identity so the dual
        //                            round trip works: l_dual(r_dual(A)) == A;
        //                            A*I differs from r_dual by the reversion sign.
        //   sta4ds (even, Minkowski): l_dual != r_dual, det(G) = -1, and the SAME two
        //                            facts as ega2d, up to the determinant:
        //          (a) r_dual(A) == rev(A) * I_4ds at every grade (the metric dual is the
        //              pseudoscalar-mult of the REVERSE). So A*I_4ds == r_dual(A) only
        //              where the reversion sign is +1 (grades 0,1,4); at grades 2,3 it is
        //              A*I_4ds == -r_dual(A). Exactly the EGA "differ by reversion".
        //          (b) the metric-dual round trip is det(G) * A -- a CLEAN global sign:
        //              l_dual(r_dual(A)) == r_dual(l_dual(A)) == det(G) * A == -A at
        //              every grade (no grade involution). Genuine recovery: divide by
        //              det(G), or the metric-free complement pair l_cmpl(r_cmpl(A)) == A,
        //              or the other school's (A * I_4ds) * inv(I_4ds) == A.
        ////////////////////////////////////////////////////////////////////////////////

        auto const s = scalar4ds{2.0};
        auto const v = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto const t = trivec4ds{1.0, 2.0, 3.0, 4.0};
        auto const ps = pscalar4ds{5.0};

        // inv(I_4ds): no graded inv(PScalar4ds) overload exists, but I_4ds*I_4ds == -1
        // makes -I_4ds the multiplicative inverse (verified below).
        auto const I_inv = -I_4ds;

        auto const sf = [](auto const& x) { return fmt::format("{}", x); };
        auto const cmp = [&sf](std::string_view g, auto const& a, auto const& b,
                               auto const& c, std::string_view note) {
            fmt::println("   {:^5} | {:<22} | {:<24} | {:<24} | {}", g, sf(a), sf(b),
                         sf(c), note);
        };
        auto const rule = [] { fmt::println("   {:-<98}", ""); };

        // --- school-1 r_dual vs school-2 A*I_4ds: differ by the reversion sign (as EGA)
        // ---
        fmt::println(
            "   FORWARD: ours r_dual(A) vs theirs A*I_4ds  (differ by rev at gr 2,3):");
        rule();
        fmt::println("   {:^5} | {:<22} | {:<24} | {:<24} | {}", "grade", "primal A",
                     "ours: r_dual(A)", "theirs: A*I_4ds", "relation");
        rule();
        cmp("0", s, r_dual(s), s * I_4ds, "equal");
        cmp("1", v, r_dual(v), v * I_4ds, "equal");
        cmp("2", B, r_dual(B), B * I_4ds, "A*I = -r_dual");
        cmp("3", t, r_dual(t), t * I_4ds, "A*I = -r_dual");
        cmp("4", ps, r_dual(ps), ps * I_4ds, "equal");
        fmt::println("");

        // clean identity at every grade: the metric dual is pseudoscalar-mult of the
        // reverse
        CHECK(r_dual(s) == rev(s) * I_4ds);
        CHECK(r_dual(v) == rev(v) * I_4ds);
        CHECK(r_dual(B) == rev(B) * I_4ds);
        CHECK(r_dual(t) == rev(t) * I_4ds);
        CHECK(r_dual(ps) == rev(ps) * I_4ds);

        // school-2 A*I_4ds coincides with r_dual where rev = +1 (grades 0,1,4) and flips
        // at grades 2,3 (where rev = -1)
        CHECK(s * I_4ds == r_dual(s));
        CHECK(v * I_4ds == r_dual(v));
        CHECK(B * I_4ds == -r_dual(B));
        CHECK(t * I_4ds == -r_dual(t));
        CHECK(ps * I_4ds == r_dual(ps));

        // basis spot-checks across grades
        CHECK(one_4ds * I_4ds == r_dual(one_4ds));
        CHECK(g1_4ds * I_4ds == r_dual(g1_4ds));
        CHECK(g14_4ds * I_4ds == -r_dual(g14_4ds));   // bivector: reversion flips
        CHECK(g234_4ds * I_4ds == -r_dual(g234_4ds)); // trivector: reversion flips
        CHECK(I_4ds * I_4ds == r_dual(I_4ds));

        // --- l_dual vs r_dual: even-dimensional split (odd grades flip sign) ---
        fmt::println("   l_dual(A) vs r_dual(A)  (even dim: odd grades flip):");
        rule();
        fmt::println("   {:^5} | {:<22} | {:<24} | {:<24} | {}", "grade", "primal A",
                     "ours: l_dual(A)", "ours: r_dual(A)", "relation");
        rule();
        cmp("0", s, l_dual(s), r_dual(s), "equal");
        cmp("1", v, l_dual(v), r_dual(v), "l_dual = -r_dual");
        cmp("2", B, l_dual(B), r_dual(B), "equal");
        cmp("3", t, l_dual(t), r_dual(t), "l_dual = -r_dual");
        cmp("4", ps, l_dual(ps), r_dual(ps), "equal");
        fmt::println("");

        CHECK(l_dual(s) == r_dual(s));
        CHECK(l_dual(v) == -r_dual(v));
        CHECK(l_dual(B) == r_dual(B));
        CHECK(l_dual(t) == -r_dual(t));
        CHECK(l_dual(ps) == r_dual(ps));
        // basis: odd-grade vectors / trivectors flip, even-grade bivectors agree
        CHECK(l_dual(g1_4ds) == -r_dual(g1_4ds));
        CHECK(l_dual(g234_4ds) == -r_dual(g234_4ds));
        CHECK(l_dual(g14_4ds) == r_dual(g14_4ds));

        ////////////////////////////////////////////////////////////////////////////////
        // full multivector: r_dual(M) == rev(M) * I_4ds (the reverse absorbs the
        // per-grade reversion sign). Plain M * I_4ds differs from r_dual(M) at the
        // grade-2,3 parts.
        ////////////////////////////////////////////////////////////////////////////////

        auto const M = mvec4ds{s, v, B, t, ps};

        fmt::println("   full multivector  M = {}:", M);
        fmt::println("     {:<22} = {}", "ours : r_dual(M)", sf(r_dual(M)));
        fmt::println("     {:<22} = {}", "them : rev(M) * I_4ds", sf(rev(M) * I_4ds));
        fmt::println("     {:<22} = {}", "them : M * I_4ds", sf(M * I_4ds));
        fmt::println("     {:<22} = {}", "ours : l_dual(M)", sf(l_dual(M)));
        fmt::println(
            "     (r_dual == rev(M)*I_4ds exactly; plain M*I_4ds differs at gr 2,3)");
        fmt::println("");

        CHECK(rev(M) * I_4ds == r_dual(M));

        ////////////////////////////////////////////////////////////////////////////////
        // inverse dualization (no Minkowski surprise under P-unify):
        //
        //   - the two metric-dual compositions AGREE with each other (as in any even
        //     dimension):  l_dual(r_dual(A)) == r_dual(l_dual(A))
        //   - they do NOT recover A: they return det(G) * A == -A, a CLEAN global sign at
        //     EVERY grade (det(G) = -1) -- NOT the grade-dependent det(G)*gr_inv that the
        //     old Q-metric produced. (For a Euclidean even algebra det = +1 and this
        //     composition IS the identity, exactly as in ega2d.)
        //   - genuine recovery: divide by det(G) (a single -1), the metric-free
        //   complement
        //     pair, or the other school's multiplication by inv(I_4ds).
        //
        // METRIC-INDEPENDENCE: det(G) = (-1)^q with q the number of negative directions,
        // so BOTH physically-equivalent STA signatures give -1 ((-,-,-,+) has q=3,
        // (+,+,+,-) has q=1, and (-1)^3 == (-1)^1 == -1). (I_4ds*I_4ds == det(G) since
        // sigma(n) == +1 for n = 4.)
        ////////////////////////////////////////////////////////////////////////////////

        value_t const detG =
            value_t(I_4ds * I_4ds); // == det(G) == -1 for both signatures

        fmt::println("   pseudoscalar / round-trip facts:");
        fmt::println("     {:<26} = {}", "I_4ds * I_4ds (= det(G))", sf(I_4ds * I_4ds));
        fmt::println("     {:<26} = {}", "I_4ds * inv(I_4ds)", sf(I_4ds * I_inv));
        fmt::println("     {:<26} = {}  (= det(G)*s = -s)", "l_dual(r_dual(s))",
                     sf(l_dual(r_dual(s))));
        fmt::println("     {:<26} = {}  (= det(G)*v = -v)", "l_dual(r_dual(v))",
                     sf(l_dual(r_dual(v))));
        fmt::println("");

        CHECK(I_4ds * I_4ds == -one_4ds);
        CHECK(I_4ds * I_inv == one_4ds); // -I_4ds is the multiplicative inverse

        // the two metric-dual compositions agree with each other ...
        CHECK(l_dual(r_dual(s)) == r_dual(l_dual(s)));
        CHECK(l_dual(r_dual(v)) == r_dual(l_dual(v)));
        CHECK(l_dual(r_dual(B)) == r_dual(l_dual(B)));
        CHECK(l_dual(r_dual(t)) == r_dual(l_dual(t)));
        CHECK(l_dual(r_dual(ps)) == r_dual(l_dual(ps)));

        // ... and return det(G) * A == -A at EVERY grade (clean global sign, no gr_inv)
        CHECK(l_dual(r_dual(s)) == detG * s);
        CHECK(l_dual(r_dual(v)) == detG * v);
        CHECK(l_dual(r_dual(B)) == detG * B);
        CHECK(l_dual(r_dual(t)) == detG * t);
        CHECK(l_dual(r_dual(ps)) == detG * ps);

        // healing is now just the determinant (no grade involution): divide by det(G)
        CHECK(detG * l_dual(r_dual(s)) == s);
        CHECK(detG * l_dual(r_dual(v)) == v);
        CHECK(detG * l_dual(r_dual(B)) == B);
        CHECK(detG * l_dual(r_dual(t)) == t);
        CHECK(detG * l_dual(r_dual(ps)) == ps);

        // equivalently: det(G) is a global sign with det(G)^2 = 1, so applying the
        // metric-dual pair TWICE recovers A (no determinant correction needed at all)
        CHECK(l_dual(r_dual(l_dual(r_dual(s)))) == s);
        CHECK(l_dual(r_dual(l_dual(r_dual(v)))) == v);
        CHECK(l_dual(r_dual(l_dual(r_dual(B)))) == B);
        CHECK(l_dual(r_dual(l_dual(r_dual(t)))) == t);
        CHECK(l_dual(r_dual(l_dual(r_dual(ps)))) == ps);

        // genuine recovery (1): the metric-free complement pair round-trips for all
        // grades
        CHECK(r_cmpl(l_cmpl(s)) == s);
        CHECK(l_cmpl(r_cmpl(s)) == s);
        CHECK(r_cmpl(l_cmpl(v)) == v);
        CHECK(l_cmpl(r_cmpl(v)) == v);
        CHECK(r_cmpl(l_cmpl(B)) == B);
        CHECK(l_cmpl(r_cmpl(B)) == B);
        CHECK(r_cmpl(l_cmpl(t)) == t);
        CHECK(l_cmpl(r_cmpl(t)) == t);
        CHECK(r_cmpl(l_cmpl(ps)) == ps);
        CHECK(l_cmpl(r_cmpl(ps)) == ps);

        // genuine recovery (2): the other school multiplies by inv(I_4ds)
        CHECK((s * I_4ds) * I_inv == s);
        CHECK((v * I_4ds) * I_inv == v);
        CHECK((B * I_4ds) * I_inv == B);
        CHECK((t * I_4ds) * I_inv == t);
        CHECK((ps * I_4ds) * I_inv == ps);

        ////////////////////////////////////////////////////////////////////////////////
        // concrete round trips that recover the input, for all grades:
        //   ours : forward = r_cmpl, backward = l_cmpl (the metric-free complement pair)
        //   them : forward = *I_4ds, backward = *inv(I_4ds)
        ////////////////////////////////////////////////////////////////////////////////

        auto const rt_row = [&sf](std::string_view input, std::string_view school,
                                  auto const& fwd, auto const& bwd) {
            fmt::println("   {:<11} | {:<6} | {:<24} | {}", input, school, sf(fwd),
                         sf(bwd));
        };
        fmt::println("   concrete round trips (forward, then backward recovers input):");
        fmt::println("   {:-<76}", "");
        fmt::println("   {:<11} | {:<6} | {:<24} | {}", "input", "school", "forward",
                     "backward (-> input)");
        fmt::println("   {:-<76}", "");
        rt_row("s  (gr0)", "ours", r_cmpl(s), l_cmpl(r_cmpl(s)));
        rt_row("", "them", s * I_4ds, (s * I_4ds) * I_inv);
        rt_row("v  (gr1)", "ours", r_cmpl(v), l_cmpl(r_cmpl(v)));
        rt_row("", "them", v * I_4ds, (v * I_4ds) * I_inv);
        rt_row("B  (gr2)", "ours", r_cmpl(B), l_cmpl(r_cmpl(B)));
        rt_row("", "them", B * I_4ds, (B * I_4ds) * I_inv);
        rt_row("t  (gr3)", "ours", r_cmpl(t), l_cmpl(r_cmpl(t)));
        rt_row("", "them", t * I_4ds, (t * I_4ds) * I_inv);
        rt_row("ps (gr4)", "ours", r_cmpl(ps), l_cmpl(r_cmpl(ps)));
        rt_row("", "them", ps * I_4ds, (ps * I_4ds) * I_inv);
        fmt::println("");
    }


    TEST_CASE("sta4ds: is_close and is_same_transform")
    {
        fmt::println("sta4ds: is_close and is_same_transform");

        auto const big = vec4ds{1.0e6, 2.0e6, 3.0e6, 4.0e6};
        auto const big_1ulp = vec4ds{std::nextafter(big.x, 1e9), big.y, big.z, big.w};
        CHECK(big != big_1ulp);
        CHECK(is_close(big, big_1ulp));
        CHECK(!is_close(big, vec4ds{big.x * (1.0 + 1e-9), big.y, big.z, big.w}));

        // An STA rotor double-covers its transformation just as a Euclidean one does --
        // and it covers BOOSTS as well as spatial rotations, which is why the function
        // is named transform rather than rotation (matching the transform() sandwich).
        auto const R_rot = get_rotor(g12_4ds, 0.7); // spatial rotation
        auto const R_bst = get_boost(g14_4ds, 0.5); // boost

        for (auto const& R : {R_rot, R_bst}) {
            auto const R_neg = mvec4ds_e{-R};
            CHECK(R != R_neg);                  // representations differ ...
            CHECK(is_same_transform(R, R_neg)); // ... the transformation does not
            CHECK(transform(g1_4ds, R) == transform(g1_4ds, R_neg));
            CHECK(transform(g4_4ds, R) == transform(g4_4ds, R_neg));
        }

        // a rotation and a boost of the same parameter are NOT the same transformation
        CHECK(!is_same_transform(R_rot, R_bst));
        CHECK(!is_same_transform(R_bst, get_boost(g14_4ds, 0.6)));

        fmt::println("");
    }


    TEST_CASE("MVec4ds: is_close for the odd, the full multivector and the dual number")
    {
        fmt::println("MVec4ds: is_close for the odd and the full multivector");

        auto A = mvec4ds{1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,
                         9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
        auto B = A;
        CHECK(is_close(A, B));
        B.c0 += 1.0e-6;
        CHECK(!is_close(A, B));

        auto U = mvec4ds_u{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        auto V = U;
        CHECK(is_close(U, V));
        V.c0 += 1.0e-6;
        CHECK(!is_close(U, V));
    }

    TEST_CASE("PScalar4ds: exp of a pseudoscalar (duality rotation on bivectors)")
    {
        fmt::println("PScalar4ds: exp of a pseudoscalar (duality rotation on bivectors)");

        auto const al = 0.63;
        auto E = exp(pscalar4ds(al));

        CHECK(std::abs(value_t(gr0(E)) - std::cos(al)) < eps);
        CHECK(std::abs(value_t(gr4(E)) - std::sin(al)) < eps);
        CHECK(nrm(gr2(E)) < eps);
        CHECK(is_close(exp(pscalar4ds(al)) * exp(pscalar4ds(-0.41)),
                       exp(pscalar4ds(al - 0.41))));

        // I_4ds is NOT central: it commutes with bivectors and anticommutes with the
        // odd grades, the anticommuting part surfacing in grade 3 rather than grade 1
        auto F = bivec4ds{1.0, 0.3, -0.2, 0.5, 1.1, -0.7};
        auto v = vec4ds{0.3, -0.5, 0.2, 1.4};
        CHECK(is_close(bivec4ds(gr2(mvec4ds(E) * mvec4ds(F))),
                       bivec4ds(gr2(mvec4ds(F) * mvec4ds(E)))));
        CHECK(is_close(trivec4ds(gr3(mvec4ds(E) * mvec4ds(v))),
                       trivec4ds(-gr3(mvec4ds(v) * mvec4ds(E)))));

        // and it is not a unit versor: E rev(E) = exp(2 alpha I), not 1
        CHECK(is_close(E * rev(E), exp(pscalar4ds(2.0 * al))));
    }

    TEST_CASE("BiVec4ds: invariant decomposition (is_simple, boost_part, rot_part)")
    {
        fmt::println(
            "BiVec4ds: invariant decomposition (is_simple, boost_part, rot_part)");

        auto F = bivec4ds{1.0, 0.3, -0.2, 0.5, 1.1, -0.7}; // non-simple
        auto S = wdg(vec4ds{0.3, -0.5, 0.2, 1.4}, vec4ds{1.1, 0.4, -0.7, 2.0}); // simple
        auto Tl = g14_4ds;                      // simple, time-like
        auto Sp = g12_4ds;                      // simple, space-like
        auto Nl = wdg(g4_4ds + g1_4ds, g2_4ds); // simple, light-like

        CHECK(!is_simple(F));
        CHECK(is_simple(S));
        CHECK(is_simple(Tl));
        CHECK(is_simple(Sp));
        CHECK(is_simple(Nl));

        // B = boost_part + rot_part, both simple, squares straddling zero
        CHECK(is_close(F, boost_part(F) + rot_part(F)));
        CHECK(is_simple(boost_part(F)));
        CHECK(is_simple(rot_part(F)));
        CHECK(value_t(gr0(boost_part(F) * boost_part(F))) >= 0.0);
        CHECK(value_t(gr0(rot_part(F) * rot_part(F))) <= 0.0);

        // orthogonal and commuting
        CHECK(std::abs(value_t(dot(boost_part(F), rot_part(F)))) < eps);
        CHECK(nrm(cmt(boost_part(F), rot_part(F))) < eps);

        // the two planes commute, so the rotor factors -- this is what exp() relies on
        CHECK(is_close(exp(F), exp(boost_part(F)) * exp(rot_part(F))));
        CHECK(is_close(exp(boost_part(F)) * exp(rot_part(F)),
                       exp(rot_part(F)) * exp(boost_part(F))));

        // a simple input lands in the slot matching its causal character
        CHECK(is_close(boost_part(Tl), Tl));
        CHECK(nrm(rot_part(Tl)) < eps);
        CHECK(is_close(rot_part(Sp), Sp));
        CHECK(nrm(boost_part(Sp)) < eps);
        CHECK(is_close(boost_part(Nl), Nl));
        CHECK(nrm(rot_part(Nl)) < eps);
    }

    TEST_CASE("MVec4ds_E: exp / log / sqrt unchanged by the shared decomposition")
    {
        fmt::println("MVec4ds_E: exp / log / sqrt unchanged by the shared decomposition");

        // regression guard: exp() was refactored onto detail::sta4ds_biv_decompose()
        auto F = bivec4ds{1.0, 0.3, -0.2, 0.5, 1.1, -0.7}; // non-simple
        auto S = wdg(vec4ds{0.3, -0.5, 0.2, 1.4}, vec4ds{1.1, 0.4, -0.7, 2.0});
        auto Tl = 0.6 * g14_4ds;
        auto Sp = 0.4 * g12_4ds;
        auto Nl = wdg(g4_4ds + g1_4ds, g2_4ds);
        auto Z = bivec4ds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        for (auto const& B : {F, S, Tl, Sp, Nl, Z}) {
            auto R = exp(B);
            CHECK(std::abs(value_t(gr0(R * rev(R))) - 1.0) < eps); // unit rotor
        }
        // round trips
        CHECK(is_close(exp(log(exp(Tl))), exp(Tl)));
        CHECK(is_close(exp(log(exp(Sp))), exp(Sp)));
        CHECK(is_close(exp(log(exp(F))), exp(F)));
        CHECK(is_close(sqrt(exp(F)) * sqrt(exp(F)), exp(F)));
    }


    TEST_CASE("Vec4ds: reciprocal_frame (the non-Euclidean case)")
    {
        fmt::println("Vec4ds: reciprocal_frame (the non-Euclidean case)");

        // with the metric {-1,-1,-1,+1} the reciprocals of the BASIS are not the basis:
        // this is what a vector derivative written over e_mu instead of e^mu gets wrong
        auto r = reciprocal_frame(g1_4ds, g2_4ds, g3_4ds, g4_4ds);
        CHECK(is_close(r[0], -g1_4ds));
        CHECK(is_close(r[1], -g2_4ds));
        CHECK(is_close(r[2], -g3_4ds));
        CHECK(is_close(r[3], g4_4ds));

        vec4ds const a[4] = {g1_4ds, g2_4ds, g3_4ds, g4_4ds};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                CHECK(std::abs(value_t(dot(r[i], a[j])) - (i == j ? 1.0 : 0.0)) < eps);
            }
        }
        // the frame is orthogonal, so the shortcut inv(a_i) agrees
        CHECK(is_close(r[0], inv(g1_4ds)));
        CHECK(is_close(r[3], inv(g4_4ds)));

        // a genuinely skewed spacetime frame: boosted, with one axis tilted
        auto R = get_boost(wdg(g1_4ds, g4_4ds), 0.4);
        vec4ds const b[4] = {transform(g1_4ds, R), transform(g2_4ds, R) + 0.3 * g1_4ds,
                             transform(g3_4ds, R), transform(g4_4ds, R)};
        auto rb = reciprocal_frame(b[0], b[1], b[2], b[3]);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                CHECK(std::abs(value_t(dot(rb[i], b[j])) - (i == j ? 1.0 : 0.0)) < eps);
            }
        }
        CHECK_THROWS_AS(reciprocal_frame(g1_4ds, g2_4ds, g3_4ds, g1_4ds + g2_4ds),
                        std::invalid_argument);
    }


    TEST_CASE("sta4ds calculus: the spacetime vector derivative and Maxwell")
    {
        fmt::println("sta4ds calculus: the spacetime vector derivative and Maxwell");

        // nrm() is NULL-CAPABLE in a Minkowski metric: nrm(X) == 0 does NOT mean X == 0,
        // so a result is tested for zero COMPONENTWISE
        auto is_zero = [](mvec4ds const& M, double tol) {
            double const* c = &M.c0;
            for (int i = 0; i < 16; ++i) {
                if (std::abs(c[i]) > tol) return false;
            }
            return true;
        };

        auto x0 = vec4ds{0.31, -0.22, 0.44, 0.57}; // (x, y, z, t)

        // the operator sums over the RECIPROCAL basis, so the space components of a
        // gradient carry the metric sign: grad(2x + 3t) = -2 g1 + 3 g4
        auto phi = [](vec4ds const& x) { return scalar4ds(2.0 * x.x + 3.0 * x.w); };
        CHECK(is_zero(mvec4ds(nabla_wdg(phi, x0) - vec4ds{-2.0, 0.0, 0.0, 3.0}), 1.0e-8));

        // a vacuum plane wave is built on the NULL propagation vector k = g4 + g3
        double const w = 1.1, E0 = 0.7;
        auto F = [&](vec4ds const& x) {
            double const ph = w * (x.z - x.w);
            return E0 * std::cos(ph) * wdg(g4_4ds + g3_4ds, g1_4ds);
        };
        // Maxwell in one equation, and its two grade halves: sources and monopoles
        CHECK(is_zero(nabla(F, x0, 1.0e-5), 1.0e-8));
        CHECK(is_zero(mvec4ds(nabla_dot(F, x0, 1.0e-5)), 1.0e-8));
        CHECK(is_zero(mvec4ds(nabla_wdg(F, x0, 1.0e-5)), 1.0e-8));
        CHECK(is_zero(dalembertian(F, x0, 1.0e-3), 1.0e-5));

        // a field that does not solve Maxwell must be detected -- otherwise the checks
        // above would pass for anything
        auto G = [](vec4ds const& x) { return (x.x * x.w) * g14_4ds + x.y * g23_4ds; };
        CHECK(!is_zero(nabla(G, x0, 1.0e-5), 1.0e-3));

        // the specific trap: this wave has the WRONG magnetic pairing, so nabla F has
        // components of order 0.1 -- yet its Minkowski norm is exactly zero. A test
        // written with nrm() would pass here; a componentwise one does not.
        auto Bad = [&](vec4ds const& x) {
            double const ph = w * (x.z - x.w);
            return E0 * std::cos(ph) * g14_4ds + E0 * std::cos(ph) * g12_4ds;
        };
        auto Nb = nabla(Bad, x0, 1.0e-5);
        CHECK(std::abs(value_t(nrm(Nb))) < 1.0e-12); // the norm sees nothing
        CHECK(!is_zero(Nb, 1.0e-3));                 // the components do

        // the d'Alembertian carries the signature: box = d_t^2 - d_x^2 - d_y^2 - d_z^2,
        // so box(t^2 + 2x^2) = 2 - 4 = -2
        auto q = [](vec4ds const& x) { return scalar4ds(x.w * x.w + 2.0 * x.x * x.x); };
        CHECK(std::abs(value_t(gr0(dalembertian(q, x0, 1.0e-3))) + 2.0) < 1.0e-5);

        // grade extremes vanish, and compact schemes are rejected as elsewhere
        auto ps = [](vec4ds const& x) { return pscalar4ds(x.x * x.w); };
        auto sc = [](vec4ds const& x) { return scalar4ds(x.y); };
        CHECK(std::abs(value_t(nabla_wdg(ps, x0))) < eps);
        CHECK(std::abs(value_t(nabla_dot(sc, x0))) < eps);
        CHECK_THROWS_AS(nabla(sc, x0, compact_scheme(1)), std::invalid_argument);
    }

    TEST_CASE("sta4ds calculus: Maxwell WITH sources -- nabla F = J from a potential")
    {
        fmt::println("sta4ds calculus: Maxwell WITH sources -- nabla F = J from a "
                     "potential");

        // The existing calculus case above solves the SOURCE-FREE equation. This one
        // closes the other half: a four-potential in Lorenz gauge, its Faraday bivector
        // F = nabla ^ A, and the four-current it carries.
        //
        //   Griffiths (12.133)  F^{mu nu} = d^mu A^nu - d^nu A^mu     ->  F = nabla ^ A
        //   Griffiths (12.136)  d_mu A^mu = 0        (Lorenz gauge)  ->  nabla . A = 0
        //   Griffiths (12.127)  d_nu F^{mu nu} = mu0 J^mu            ->  gr1(nabla F) = J
        //                       d_nu G^{mu nu} = 0   (no monopoles)  ->  gr3(nabla F) = 0
        //   Griffiths (12.137)  box A^mu = -mu0 J^mu                 ->  dalembertian(A)
        //
        // Griffiths' box is nabla^2 - (1/c^2) d_t^2; the library's dalembertian is its
        // NEGATIVE (d_t^2 - nabla^2, see the header), which is why the two signs agree
        // here without a minus. Natural units, c = eps0 = mu0 = 1.
        //
        // A = (0, 0, Az, V) with Az = (x^2+y^2)/2 and V = (x^2+y^2+z^2)/2: static, so
        // the gauge condition holds identically, and both fields are QUADRATIC, so the
        // central differences are exact and the nested application below carries only
        // round-off.
        // nrm_sq is NULL-CAPABLE in a Minkowski metric: it can vanish (or go negative)
        // for a non-zero operand, so every "is this zero?" here is COMPONENTWISE -- the
        // same rule the source-free calculus case states
        auto max_abs = [](mvec4ds const& M) {
            value_t const* c = &M.c0;
            value_t m = 0.0;
            for (int i = 0; i < 16; ++i)
                m = std::max(m, std::abs(c[i]));
            return m;
        };

        value_t const h = 1.0e-3;
        auto A = [](vec4ds const& y) {
            return vec4ds{0.0, 0.0, 0.5 * (y.x * y.x + y.y * y.y),
                          0.5 * (y.x * y.x + y.y * y.y + y.z * y.z)};
        };
        auto x0 = vec4ds{0.31, -0.22, 0.44, 0.57}; // (x, y, z, t)

        // ---- the gauge condition, and the four-current from the wave operator -------
        CHECK(std::abs(value_t(nabla_dot(A, x0, h))) < 1.0e-9);
        auto const J = gr1(dalembertian(A, x0, h));
        CHECK(J.x == doctest::Approx(0.0).epsilon(1.0e-9));
        CHECK(J.y == doctest::Approx(0.0).epsilon(1.0e-9));
        CHECK(J.z == doctest::Approx(-2.0)); // -laplacian(Az) = -2
        CHECK(J.w == doctest::Approx(-3.0)); // -laplacian(V)  = -3

        // ---- the Faraday bivector, and the classical fields it carries --------------
        auto Ffield = [&](vec4ds const& y) { return nabla_wdg(A, y, h); };
        auto const F = Ffield(x0);
        // read E and B back out for the standard observer, in the convention of the
        // "Faraday bivector" test case (E on g_k4 with +, B on g_jk with -)
        auto const Ep = rel_vec_split(F, g4_4ds);
        auto const Bp = rel_bivec_split(F, g4_4ds);
        // classical references, computed by hand from A: E = -grad V = -r (static),
        // B = curl A = (y, -x, 0)
        CHECK(Ep.vx == doctest::Approx(-x0.x));
        CHECK(Ep.vy == doctest::Approx(-x0.y));
        CHECK(Ep.vz == doctest::Approx(-x0.z));
        CHECK(-Bp.mx == doctest::Approx(x0.y));  // Bx =  y
        CHECK(-Bp.my == doctest::Approx(-x0.x)); // By = -x
        CHECK(-Bp.mz == doctest::Approx(0.0).epsilon(1.0e-9));

        // ---- Maxwell in ONE equation: nabla F = J ----------------------------------
        auto const N = nabla(Ffield, x0, h);
        CHECK(max_abs(mvec4ds(gr1(N) - J)) < 1.0e-9); // (12.127) first equation
        CHECK(max_abs(mvec4ds(gr3(N))) < 1.0e-9);     // (12.127) second: no monopoles
        // the remaining grades are empty: nabla F has no scalar, bivector or
        // pseudoscalar part for a bivector field
        CHECK(std::abs(value_t(gr0(N))) < 1.0e-9);
        CHECK(max_abs(mvec4ds(gr2(N))) < 1.0e-9);
        CHECK(std::abs(value_t(gr4(N))) < 1.0e-9);

        // the classical reading of the same four numbers, computed by hand:
        //   div E   = -3      -> c rho = J_w = -3           (Gauss, Griffiths 2.14)
        //   curl B  = (0,0,-2) -> mu0 J_z = J_z = -2        (Ampere, Griffiths 5.56)
        CHECK(J.w == doctest::Approx(-3.0));
        CHECK(J.z == doctest::Approx(-2.0));

        // ---- falsify: break the gauge and the source no longer equals box A ---------
        // adding a t-dependent piece to V violates nabla . A = 0, and then
        // nabla F = box A - nabla(nabla . A) picks up the extra term
        auto Ag = [](vec4ds const& y) {
            return vec4ds{0.0, 0.0, 0.5 * (y.x * y.x + y.y * y.y),
                          0.5 * (y.x * y.x + y.y * y.y + y.z * y.z) + 0.7 * y.w * y.z};
        };
        CHECK(std::abs(value_t(nabla_dot(Ag, x0, h))) > 0.1); // gauge broken
        auto Fg = [&](vec4ds const& y) { return nabla_wdg(Ag, y, h); };
        // NOTE the componentwise form: the difference here has nrm_sq = -0.49, so a
        // "> 0.1" test on nrm_sq FAILS on a residual that is plainly non-zero
        CHECK(max_abs(mvec4ds(gr1(nabla(Fg, x0, h)) - gr1(dalembertian(Ag, x0, h)))) >
              0.1);

        fmt::println("nabla F = J from a Lorenz-gauge potential; grade 1 = sources, "
                     "grade 3 = no monopoles");
    }

    TEST_CASE("MVec4ds_E: the null rotor -- the third kind of Lorentz transformation")
    {
        fmt::println(
            "MVec4ds_E: the null rotor -- the third kind of Lorentz transformation");

        // Jancewicz sec.6.3 names three kinds of Lorentz transformation, one per causal
        // character of the generating plane: a BOOST from a time-like volutor, a spatial
        // ROTATION from a space-like one, and "a Lorentz transformation of a new kind"
        // from a LIGHT-LIKE volutor W = k ^ n with k^2 = 0, k.n = 0, n^2 = 1. The third
        // is the null (parabolic) rotation. Because W^2 = 0 the exponential series
        // truncates after two terms -- exp(W) = 1 + W, exactly, no cos/cosh involved.
        //
        // This is the branch of sta4ds_exp_simple() that neither ega nor a Euclidean
        // metric can reach: it exists only because a non-zero plane can square to zero.

        auto k = g4_4ds + g1_4ds; // null vector, k^2 = +1 - 1 = 0
        auto n = g2_4ds;          // unit space-like, orthogonal to k
        auto W = wdg(k, n);

        CHECK(is_lightlike(k));
        CHECK(std::abs(value_t(dot(k, n))) < eps);
        CHECK(is_lightlike(W)); // W^2 = 0 -- neither a boost nor a rotation plane
        CHECK(is_simple(W));

        // the series truncates: exp(W) == 1 + W with no transcendental factor
        auto R = exp(W);
        CHECK(is_close(R, mvec4ds_e(scalar4ds(1.0), W, pscalar4ds(0.0))));

        // it is a genuine unit versor, so transform() is a motion
        CHECK(std::abs(value_t(gr0(rev(R) * R)) - 1.0) < eps);

        // log inverts the null branch as well (bn ~ 0 -> return Bv)
        CHECK(is_close(log(R), W));

        // PARABOLIC: R - 1 is nilpotent. This is what separates the null rotor from the
        // other two kinds -- a boost or a rotation rotor has no nilpotent part.
        auto Rm1 = mvec4ds_e(scalar4ds(0.0), W, pscalar4ds(0.0));
        CHECK(is_close(Rm1 * Rm1,
                       mvec4ds_e(scalar4ds(0.0), bivec4ds(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                                 pscalar4ds(0.0))));

        // the book's defining property: k W = -W k = 0, so the null rotation FIXES its
        // own null direction while moving everything else -- the parabolic signature
        CHECK(nrm(mvec4ds(k) * mvec4ds_e(W)) < eps);
        CHECK(is_close(transform(k, R), k));
        CHECK(!is_close(transform(n, R), n)); // not the identity

        // and it IS a Lorentz transformation: the pseudolength is preserved at every
        // grade, exactly as for a boost or a rotation
        // RELATIVE tolerance: these pass through a sandwich product, so the difference
        // accumulates several roundings and machine eps is not the right yardstick for
        // an O(1)..O(10) quantity (the vec case lands at 1.3e-15, just over it)
        auto preserved = [](value_t before, value_t after) {
            return std::abs(after - before) <=
                   eps_congruent * std::max(std::abs(before), 1.0);
        };
        vec4ds const x{0.3, -0.5, 0.2, 1.4};
        bivec4ds const B{1.0, 0.3, -0.2, 0.5, 1.1, -0.7};
        trivec4ds const T{1.0, 2.0, 3.0, 0.5};
        CHECK(preserved(value_t(gr0(x * x)),
                        value_t(gr0(transform(x, R) * transform(x, R)))));
        CHECK(preserved(value_t(gr0(B * B)),
                        value_t(gr0(transform(B, R) * transform(B, R)))));
        CHECK(preserved(value_t(gr0(T * T)),
                        value_t(gr0(transform(T, R) * transform(T, R)))));
    }

    TEST_CASE("sta4ds: the complement is NOT the dual (ega cannot see this)")
    {
        fmt::println("sta4ds: the complement is NOT the dual (ega cannot see this)");

        // Two maps that ega conflates: the complement is combinatorial (which basis
        // blade completes this one to the pseudoscalar), the dual is that composed with
        // the metric. Under an identity metric they are the same function, so every
        // ega relation reads correctly with either name and a wrong choice there is
        // invisible. This algebra separates them, and pins which one is which.

        vec4ds const v{0.3, -0.5, 0.2, 1.4};
        bivec4ds const B{1.0, 0.3, -0.2, 0.5, 1.1, -0.7};
        trivec4ds const T{1.0, 2.0, 3.0, 0.5};

        // they differ on exactly the blades whose extended-metric entry is -1
        CHECK(!is_close(l_cmpl(v), l_dual(v)));
        CHECK(!is_close(l_cmpl(B), l_dual(B)));
        CHECK(!is_close(l_cmpl(T), l_dual(T)));

        // the pseudoscalar picks out the DUAL, and of the REVERSE -- both halves of
        // that statement are load-bearing, and each fails on its own
        CHECK(is_close(trivec4ds(I_4ds * v), l_dual(rev(v))));
        CHECK(is_close(bivec4ds(I_4ds * B), l_dual(rev(B))));
        CHECK(is_close(vec4ds(I_4ds * T), l_dual(rev(T))));

        CHECK(!is_close(trivec4ds(I_4ds * v), l_cmpl(rev(v)))); // not the complement
        CHECK(!is_close(bivec4ds(I_4ds * B), l_cmpl(rev(B))));
        CHECK(!is_close(vec4ds(I_4ds * T), l_cmpl(rev(T))));

        CHECK(!is_close(bivec4ds(I_4ds * B), l_dual(B))); // not without rev()
        CHECK(!is_close(vec4ds(I_4ds * T), l_dual(T)));

        // Jancewicz Problem 6.14: the Hodge map of a time-like volutor is space-like and
        // vice versa; of a light-like volutor, light-like. The book DEFINES the Hodge map
        // as left multiplication by the unit quadrivector (sec.6.3, "M -> J M") -- which
        // is the identity pinned just above, our l_dual composed with a reversion.
        auto const Tl = bivec4ds(g14_4ds);            // time-like,  B^2 = +1
        auto const Sp = bivec4ds(g12_4ds);            // space-like, B^2 = -1
        auto const Nl = wdg(g4_4ds + g1_4ds, g2_4ds); // light-like, B^2 =  0
        auto csq = [](bivec4ds const& b) { return value_t(gr0(b * b)); };

        CHECK(csq(Tl) > 0.0);
        CHECK(csq(l_dual(Tl)) < 0.0); // time  -> space
        CHECK(csq(Sp) < 0.0);
        CHECK(csq(l_dual(Sp)) > 0.0); // space -> time
        CHECK(std::abs(csq(Nl)) < eps);
        CHECK(std::abs(csq(l_dual(Nl))) < eps); // null  -> null

        // This does NOT discriminate the dual from the complement, though it reads as if
        // it should: the two differ by a per-blade metric sign, and causal character is
        // QUADRATIC in the blade, so the sign cancels and l_cmpl reproduces the table
        // exactly. The discriminating checks are the ones above; do not let this weaker
        // property be substituted for them.
        CHECK(csq(l_cmpl(Tl)) < 0.0); // the complement passes it too
        CHECK(csq(l_cmpl(Sp)) > 0.0);
    }

} // STA 3D Tests
