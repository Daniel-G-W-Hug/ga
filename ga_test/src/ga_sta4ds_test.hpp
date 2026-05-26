// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <chrono>
#include <iostream>
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

        // (1) Signature pinned via nrm_sq() on each unit basis blade. The values
        //     follow the STA extended metric (g1^2 = g2^2 = g3^2 = -1; g4^2 = +1;
        //     higher grades via G(a^b) = G(a)^G(b)).
        // vectors: g1,g2,g3 -> -1 (spacelike); g4 -> +1 (timelike)
        CHECK(nrm_sq(g1_4ds) == -1.0);
        CHECK(nrm_sq(g2_4ds) == -1.0);
        CHECK(nrm_sq(g3_4ds) == -1.0);
        CHECK(nrm_sq(g4_4ds) == 1.0);
        // bivectors: with g4 -> +1 (timelike); without g4 -> -1 (spacelike)
        CHECK(nrm_sq(g14_4ds) == 1.0);
        CHECK(nrm_sq(g24_4ds) == 1.0);
        CHECK(nrm_sq(g34_4ds) == 1.0);
        CHECK(nrm_sq(g23_4ds) == -1.0);
        CHECK(nrm_sq(g31_4ds) == -1.0);
        CHECK(nrm_sq(g12_4ds) == -1.0);
        // trivectors: with g4 -> -1 (spacelike); g123 (without g4) -> 1 (timelike)
        CHECK(nrm_sq(g234_4ds) == -1.0);
        CHECK(nrm_sq(g314_4ds) == -1.0);
        CHECK(nrm_sq(g124_4ds) == -1.0);
        CHECK(nrm_sq(g123_4ds) == 1.0);
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
    // causal character, nrm and normalization of blades
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: causal character, nrm and normalization of blades")
    {
        fmt::println("G<1,3,0>: causal character, nrm and normalization of blades");

        // ---- causal classification (grades 1..n-1: Vec, BiVec, TriVec) ----
        // timelike: nrm_sq > 0 ; spacelike: nrm_sq < 0 ; lightlike: nrm_sq == 0
        // vectors: g1,g2,g3 -> spacelike ; g4 -> timelike
        CHECK(is_spacelike(g1_4ds));
        CHECK(is_spacelike(g2_4ds));
        CHECK(is_spacelike(g3_4ds));
        CHECK(is_timelike(g4_4ds));
        // null vector g1 + g4 -> nrm_sq = +1 - 1 = 0
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
        // null bivector g14 + g23 -> nrm_sq = +1 - 1 = 0
        bivec4ds Bl{1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        CHECK(is_lightlike(Bl));

        // trivectors: g234,g314,g124 -> spacelike ; g123 -> timelike
        CHECK(is_spacelike(g234_4ds));
        CHECK(is_spacelike(g314_4ds));
        CHECK(is_spacelike(g124_4ds));
        CHECK(is_timelike(g123_4ds));
        // null trivector g234 + g123 -> nrm_sq = -1 + 1 = 0
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

        // ---- normalize(): timelike -> +1, spacelike -> -1, lightlike unchanged ----
        // direction is preserved (scale = 1/nrm > 0, never flips sign)
        // timelike vector 3*g4: nrm_sq = 9 -> 1
        vec4ds vt{0.0, 0.0, 0.0, 3.0};
        CHECK(nrm_sq(normalize(vt)) == doctest::Approx(1.0));
        CHECK(normalize(vt).w > 0.0);
        // spacelike vector 3*g1: nrm_sq = -9 -> -1
        vec4ds vs{3.0, 0.0, 0.0, 0.0};
        CHECK(nrm_sq(normalize(vs)) == doctest::Approx(-1.0));
        CHECK(normalize(vs).x > 0.0);
        // lightlike vector: returned unchanged (exact, no scaling applied)
        CHECK(normalize(vl) == vl);

        // bivector: timelike 2*g14 -> +1 ; spacelike 2*g23 -> -1 ; null unchanged
        CHECK(nrm_sq(normalize(bivec4ds{2.0, 0.0, 0.0, 0.0, 0.0, 0.0})) ==
              doctest::Approx(1.0));
        CHECK(nrm_sq(normalize(bivec4ds{0.0, 0.0, 0.0, 2.0, 0.0, 0.0})) ==
              doctest::Approx(-1.0));
        CHECK(normalize(Bl) == Bl);

        // trivector: timelike 2*g123 -> +1 ; spacelike 2*g234 -> -1 ; null unchanged
        CHECK(nrm_sq(normalize(trivec4ds{0.0, 0.0, 0.0, 2.0})) == doctest::Approx(1.0));
        CHECK(nrm_sq(normalize(trivec4ds{2.0, 0.0, 0.0, 0.0})) == doctest::Approx(-1.0));
        CHECK(normalize(tl) == tl);

        // scalar normalize still works (no causal logic): scales to unit sign
        CHECK(value_t(normalize(scalar4ds{3.0})) == doctest::Approx(1.0));
        CHECK(value_t(normalize(scalar4ds{-5.0})) == doctest::Approx(-1.0));

        fmt::println("causal predicates, nrm = sqrt(|nrm_sq|), normalize -> nrm_sq == "
                     "+/-1 (null unchanged)");
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
        // ground-truth values from ga_prdxpr (NOT the self-referential wrapper form
        // l_expand4ds(a,b)==wdg(l_dual(a),b), which can't catch a dual sign bug)
        CHECK(l_expand4ds(B1, v1) == trivec4ds{7.0, 2.0, 15.0, 14.0});
        CHECK(r_expand4ds(v1, B1) == trivec4ds{7.0, 2.0, 15.0, 14.0});
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
        auto const Rr = get_rotor(g12_4ds, a); // cos(a/2) + sin(-a/2) g12
        CHECK(value_t(gr0(Rr)) == doctest::Approx(std::cos(a / 2)));
        CHECK(gr2(Rr).mz == doctest::Approx(std::sin(-a / 2)));
        CHECK(value_t(gr0(rev(Rr) * Rr)) == doctest::Approx(1.0));
        CHECK(nrm_sq(gr2(rev(Rr) * Rr)) == doctest::Approx(0.0)); // no leftover bivector

        auto const Rbo = get_boost(g14_4ds, a); // cosh(a/2) + sinh(a/2) g14
        CHECK(value_t(gr0(Rbo)) == doctest::Approx(std::cosh(a / 2)));
        CHECK(gr2(Rbo).vx == doctest::Approx(std::sinh(a / 2)));
        CHECK(value_t(gr0(rev(Rbo) * Rbo)) == doctest::Approx(1.0));
        CHECK(nrm_sq(gr2(rev(Rbo) * Rbo)) == doctest::Approx(0.0));

        fmt::println("exp: spatial -> cos/sin, boost -> cosh/sinh; rotor/boost are unit");
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

        // ---- spatial rotation in the g1-g2 plane: g1 -> cos(th) g1 - sin(th) g2 ----
        // (identical to a 3D rotation; the time direction g4 is untouched)
        auto const vr = transform(g1_4ds, Rr);
        CHECK(vr.x == doctest::Approx(std::cos(th)));
        CHECK(vr.y == doctest::Approx(-std::sin(th)));
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

        fmt::println(
            "transform: interval-invariant; rotation==3D; boost gamma=cosh, "
            "beta=tanh; collinear boosts add rapidity; transform_opt==transform "
            "(vec/bivec/trivec, scalar + batch)");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ops.hpp step 3: spacetime split of a vector (time + relative space)
    ////////////////////////////////////////////////////////////////////////////////

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
        // plane; the rotor turns g1 toward -g2 for positive angle, so the mirror
        // order (g1 then b) corresponds to angle -2a in get_rotor's convention
        auto const Rrot = get_rotor(g12_4ds, -2.0 * a);
        CHECK(nrm_sq(twice - transform(v, Rrot)) == doctest::Approx(0.0));

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

    TEST_CASE("MVec4ds: metric / antimetric exomorphisms (G, rG) and metric-indep. dual")
    {
        fmt::println("MVec4ds: metric / antimetric exomorphisms (G, rG)");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Following Lengyel ("PGA Illuminated"), every metric algebra carries TWO
        // exomorphisms tied to the metric:
        //   - the metric      G : a wedge      exomorphism,  G(wdg(a,b))  == wdg(G(a),G(b))
        //   - the antimetric rG : an antiwedge exomorphism, rG(rwdg(a,b)) == rwdg(rG(a),rG(b))
        // related by  G * rG = det(metric) * I  (ordinary matrix product, I = identity).
        //
        // In PGA the metric is degenerate (det = 0), so G*rG = 0 -- G and rG annihilate
        // each other, which is exactly the bulk/weight split (disjoint support). In a
        // NON-degenerate algebra like sta4ds there is no null direction, det != 0, so the
        // relation is invertible. With det = -1 and G^2 = I (every basis blade squares to
        // +/-1) it collapses to
        //
        //                          rG = -G          (the antimetric is the negated metric)
        //
        // FAITHFUL METRIC vs LIBRARY METRIC. The faithful exomorphism is the PURE-PRODUCT
        // metric G(e_S) = (prod_{i in S} g_i) e_S. The library's own extended metric
        //     Gx(X) = l_cmpl(r_dual(X))     (since r_dual = r_cmpl o Gx, l_cmpl o r_cmpl = id)
        // is the BLADE-SQUARE form e_S^2 = sigma(k) * prod g_i (it carries the reversion
        // sign sigma(k) = (-1)^(k(k-1)/2)). The two are related by exactly the reversion:
        //     G(X) = rev(Gx(X))            (pure product = blade square with reversion removed)
        // and that reversion content is precisely the gr_inv that makes l_dual o r_dual
        // (which uses Gx on the primal side twice) miss the identity. Using the FAITHFUL
        // G / rG, all of Lengyel's relations hold cleanly at every grade.
        //
        // Signature duality (the physics remark): G(g1)=-g1 realizes (-,-,-,+) while
        // rG = -G realizes (+,+,+,-). The two signatures are wedge-metric and
        // antiwedge-antimetric of one another -- the same geometry, read off in dual
        // spaces. One school works it out in the space, the other in its dual.
        ////////////////////////////////////////////////////////////////////////////////

        auto const sf = [](auto const& x) { return fmt::format("{}", x); };
        auto const Gx = [](auto const& X) { return l_cmpl(r_dual(X)); }; // library (blade square)
        auto const G = [&Gx](auto const& X) { return rev(Gx(X)); };      // faithful wedge metric
        auto const rG = [&G](auto const& X) { return -G(X); };           // antimetric (= -G)

        auto const s = scalar4ds{2.0};
        auto const v = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto const t = trivec4ds{1.0, 2.0, 3.0, 4.0};
        auto const ps = pscalar4ds{5.0};

        fmt::println("   metric vs antimetric on the basis vectors (signature duality):");
        fmt::println("     G(g1)={:>13}  rG(g1)={:>13}", sf(G(g1_4ds)), sf(rG(g1_4ds)));
        fmt::println("     G(g4)={:>13}  rG(g4)={:>13}", sf(G(g4_4ds)), sf(rG(g4_4ds)));
        fmt::println("");

        // (0) faithful metric = library metric with the reversion removed
        CHECK(G(g14_4ds) == rev(Gx(g14_4ds)));
        CHECK(G(g23_4ds) == rev(Gx(g23_4ds)));

        // (1) G * rG = det(metric) * I, det = -1  <=>  rG(G(X)) == -X  (and rG = -G)
        CHECK(rG(G(s)) == -s);
        CHECK(rG(G(v)) == -v);
        CHECK(rG(G(B)) == -B);
        CHECK(rG(G(t)) == -t);
        CHECK(rG(G(ps)) == -ps);
        CHECK(G(v) == -rG(v)); // rG = -G

        // signature duality: G realizes (-,-,-,+), rG = -G realizes (+,+,+,-)
        CHECK(G(g1_4ds) == -g1_4ds);  CHECK(rG(g1_4ds) == g1_4ds);
        CHECK(G(g4_4ds) == g4_4ds);   CHECK(rG(g4_4ds) == -g4_4ds);

        // (2) G is a faithful WEDGE exomorphism at every grade (incl. mixed grades, where
        //     the blade-square library metric Gx would fail by the reversion sign)
        CHECK(G(wdg(g1_4ds, g4_4ds)) == wdg(G(g1_4ds), G(g4_4ds)));   // vector ^ vector
        CHECK(G(wdg(g2_4ds, g3_4ds)) == wdg(G(g2_4ds), G(g3_4ds)));   // vector ^ vector
        CHECK(G(wdg(g1_4ds, g23_4ds)) == wdg(G(g1_4ds), G(g23_4ds))); // vector ^ bivector
        CHECK(G(wdg(g4_4ds, g23_4ds)) == wdg(G(g4_4ds), G(g23_4ds))); // vector ^ bivector
        // the library's blade-square Gx is NOT a faithful exomorphism (off by reversion):
        CHECK(Gx(wdg(g1_4ds, g4_4ds)) != wdg(Gx(g1_4ds), Gx(g4_4ds)));

        // (3) rG is a faithful ANTIWEDGE exomorphism: rG(rwdg(a,b)) == rwdg(rG(a),rG(b))
        CHECK(rG(rwdg(g234_4ds, g123_4ds)) == rwdg(rG(g234_4ds), rG(g123_4ds)));
        CHECK(rG(rwdg(g234_4ds, g14_4ds)) == rwdg(rG(g234_4ds), rG(g14_4ds)));
        // ... and G is NOT an antiwedge exomorphism (the roles are swapped)
        CHECK(G(rwdg(g234_4ds, g14_4ds)) != rwdg(G(g234_4ds), G(g14_4ds)));

        ////////////////////////////////////////////////////////////////////////////////
        // Metric-independent dualization using the metric on the DUAL side.
        //
        // r_dual = r_cmpl o Gx puts the metric on the primal side; its inverse therefore
        // puts the metric on the dual side. Since Gx^2 = I, the inverse is Gx o l_cmpl,
        // i.e. the "antidual":
        //     antidual(D) = l_cmpl(r_dual(l_cmpl(D)))           (= Gx o l_cmpl)
        // This inverts r_dual at EVERY grade and depends on the metric only through Gx^2=I,
        // so it works for any non-degenerate signature -- including the equivalent
        // (+,+,+,-) -- without the det*gr_inv correction that l_dual o r_dual needs.
        // (Cannot be exercised against G(3,1,0) here, as the library only instantiates
        //  G(1,3,0); but the construction uses no signature-specific value beyond Gx^2=I.)
        ////////////////////////////////////////////////////////////////////////////////

        auto const antidual = [](auto const& D) { return l_cmpl(r_dual(l_cmpl(D))); };

        CHECK(antidual(r_dual(s)) == s);
        CHECK(antidual(r_dual(v)) == v);
        CHECK(antidual(r_dual(B)) == B);
        CHECK(antidual(r_dual(t)) == t);
        CHECK(antidual(r_dual(ps)) == ps);
        // symmetric statement: r_dual inverts the antidual too
        CHECK(r_dual(antidual(v)) == v);
        CHECK(r_dual(antidual(B)) == B);
    }

    TEST_CASE("MVec4ds: left/right antidual; the wedge product fixes space vs dual")
    {
        fmt::println("MVec4ds: left/right antidual; the wedge product fixes space vs dual");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Q1 -- there are TWO antiduals, one inverse for each handedness of the dual:
        //     right antidual (inverts r_dual):  r_undual(D) = l_cmpl(r_dual(l_cmpl(D))) = Gx o l_cmpl
        //     left  antidual (inverts l_dual):  l_undual(D) = r_cmpl(l_dual(r_cmpl(D))) = Gx o r_cmpl
        //   (metric on the dual side; the complement that opens/closes matches the dual's
        //   handedness). On even-grade inputs l_dual == r_dual, so r_undual == l_undual there;
        //   on odd grades they differ by a sign, exactly as l_dual and r_dual do. Each
        //   inverts only its own handedness -- the cross pairing flips the odd grades.
        //
        // Q2 -- space vs dual space is fixed by WHICH product is the wedge. Each signature
        //   is an exomorphism of EXACTLY ONE product (they are not interchangeable):
        //     G  = (-,-,-,+)  is the WEDGE     exomorphism, NOT the antiwedge one
        //     rG = (+,+,+,-)  is the ANTIWEDGE exomorphism, NOT the wedge one
        //   So, relative to the library's wedge convention, the space is uniquely the
        //   signature whose metric respects the wedge product: (-,-,-,+). Its dual space is
        //   (+,+,+,-). The structure is symmetric under (wdg <-> rwdg, G <-> rG,
        //   space <-> dual) -- which is exactly why both signatures are physically
        //   equivalent: each is "the space" for one of the two products.
        ////////////////////////////////////////////////////////////////////////////////

        auto const Gx = [](auto const& X) { return l_cmpl(r_dual(X)); };
        auto const G = [&Gx](auto const& X) { return rev(Gx(X)); }; // faithful (-,-,-,+)
        auto const rG = [&G](auto const& X) { return -G(X); };      // antimetric (+,+,+,-)
        auto const r_undual = [](auto const& D) { return l_cmpl(r_dual(l_cmpl(D))); };
        auto const l_undual = [](auto const& D) { return r_cmpl(l_dual(r_cmpl(D))); };

        auto const s = scalar4ds{2.0};
        auto const v = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto const t = trivec4ds{1.0, 2.0, 3.0, 4.0};
        auto const ps = pscalar4ds{5.0};

        // Q1: right antidual inverts r_dual at every grade
        CHECK(r_undual(r_dual(s)) == s);
        CHECK(r_undual(r_dual(v)) == v);
        CHECK(r_undual(r_dual(B)) == B);
        CHECK(r_undual(r_dual(t)) == t);
        CHECK(r_undual(r_dual(ps)) == ps);
        // left antidual inverts l_dual at every grade
        CHECK(l_undual(l_dual(s)) == s);
        CHECK(l_undual(l_dual(v)) == v);
        CHECK(l_undual(l_dual(B)) == B);
        CHECK(l_undual(l_dual(t)) == t);
        CHECK(l_undual(l_dual(ps)) == ps);
        // the two antiduals are genuinely distinct: each inverts only its own handedness,
        // the cross-pairing flips the odd grades (and coincides on the even grades)
        CHECK(l_undual(r_dual(v)) == -v); // left antidual does NOT invert r_dual (odd grade)
        CHECK(r_undual(l_dual(v)) == -v); // right antidual does NOT invert l_dual (odd grade)
        CHECK(l_undual(r_dual(B)) == B);  // ... but they agree on even grades
        CHECK(r_undual(l_dual(B)) == B);

        // Q2: G = (-,-,-,+) is the WEDGE exomorphism, and ONLY the wedge one
        CHECK(G(wdg(g1_4ds, g4_4ds)) == wdg(G(g1_4ds), G(g4_4ds)));    // wedge: holds
        CHECK(G(wdg(g2_4ds, g3_4ds)) == wdg(G(g2_4ds), G(g3_4ds)));
        CHECK(G(rwdg(g234_4ds, g14_4ds)) != rwdg(G(g234_4ds), G(g14_4ds))); // antiwedge: fails
        // rG = (+,+,+,-) is the ANTIWEDGE exomorphism, and ONLY the antiwedge one
        CHECK(rG(rwdg(g234_4ds, g14_4ds)) == rwdg(rG(g234_4ds), rG(g14_4ds)));   // antiwedge: holds
        CHECK(rG(rwdg(g234_4ds, g123_4ds)) == rwdg(rG(g234_4ds), rG(g123_4ds)));
        CHECK(rG(wdg(g1_4ds, g4_4ds)) != wdg(rG(g1_4ds), rG(g4_4ds)));     // wedge: fails

        fmt::println("   Q1: r_undual inverts r_dual, l_undual inverts l_dual (distinct on "
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
        //   backward (dual  -> space): A = antidual(D)              = l_cmpl(r_dual(l_cmpl(D)))
        //                                                           = (1/det) * rG * l_cmpl(D)
        //
        // The "naive" backward leg l_dual(D) (metric on the primal side again -- the recipe
        // that works in the Euclidean even algebra ega2d) does NOT recover A here. It
        // returns det*gr_inv(A) = -gr_inv(A): it recovers the odd grades by luck but flips
        // every even grade.
        ////////////////////////////////////////////////////////////////////////////////

        auto const sf = [](auto const& x) { return fmt::format("{}", x); };
        auto const antidual = [](auto const& D) { return l_cmpl(r_dual(l_cmpl(D))); };

        auto const s = scalar4ds{2.0};
        auto const v = vec4ds{1.0, 2.0, 3.0, 4.0};
        auto const B = bivec4ds{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto const t = trivec4ds{1.0, 2.0, 3.0, 4.0};
        auto const ps = pscalar4ds{5.0};

        auto const row = [&sf, &antidual](std::string_view name, auto const& A) {
            auto const D = r_dual(A);
            fmt::println("   {:<3} | {:<24} | {:<24} | {:<24} | {}", name, sf(A), sf(D),
                         sf(antidual(D)), sf(l_dual(D)));
        };
        fmt::println("   {:<3} | {:<24} | {:<24} | {:<24} | {}", "A", "value (space)",
                     "D = r_dual(A) (dual)", "antidual(D) [== A]", "l_dual(D) [naive]");
        fmt::println("   {:-<108}", "");
        row("s", s);
        row("v", v);
        row("B", B);
        row("t", t);
        row("ps", ps);
        fmt::println("");

        // forward then backward recovers A at EVERY grade (antidual = correct inverse)
        CHECK(antidual(r_dual(s)) == s);
        CHECK(antidual(r_dual(v)) == v);
        CHECK(antidual(r_dual(B)) == B);
        CHECK(antidual(r_dual(t)) == t);
        CHECK(antidual(r_dual(ps)) == ps);

        // full multivector, both orders
        auto const M = mvec4ds{s, v, B, t, ps};
        CHECK(antidual(r_dual(M)) == M);
        CHECK(r_dual(antidual(M)) == M);

        // the naive backward leg l_dual yields det*gr_inv(A): even grades flip, odd survive
        CHECK(l_dual(r_dual(s)) == -s);  // grade 0 (even) -> WRONG
        CHECK(l_dual(r_dual(B)) == -B);  // grade 2 (even) -> WRONG
        CHECK(l_dual(r_dual(ps)) == -ps);// grade 4 (even) -> WRONG
        CHECK(l_dual(r_dual(v)) == v);   // grade 1 (odd)  -> survives by luck
        CHECK(l_dual(r_dual(t)) == t);   // grade 3 (odd)  -> survives by luck

        // the backward leg expressed with the antimetric rG = -G  (uses both G and rG):
        //   antidual(D) = (1/det) * rG * l_cmpl(D),  with det = I^2 = -1, rG(Y) = -G(Y)
        value_t const det = value_t(I_4ds * I_4ds); // = -1
        auto const Gx = [](auto const& X) { return l_cmpl(r_dual(X)); }; // metric
        auto const rG = [&Gx](auto const& X) { return -Gx(X); };         // antimetric = -G
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
        // Three-way contrast across the algebras:
        //
        //   ega3d (odd, Euclidean):  one dual, it is an involution dual(dual(A))==A;
        //                            A*I differs from dual(A) by the reversion sign.
        //   ega2d (even, Euclidean): l_dual != r_dual; metric = identity so the dual
        //                            round trip works: l_dual(r_dual(A)) == A;
        //                            A*I differs from r_dual by the reversion sign.
        //   sta4ds (even, Minkowski): l_dual != r_dual, AND two surprises caused by the
        //                            mixed metric:
        //          (a) A * I_4ds == r_dual(A) EXACTLY at every grade -- the metric
        //              absorbs the reversion sign that makes them differ in EGA.
        //          (b) the metric-dual round trip is NOT the identity:
        //              l_dual(r_dual(A)) == r_dual(l_dual(A)) == det(G)-signed A
        //              (= +A for odd grades, -A for even grades). Genuine recovery
        //              comes from the metric-free complement pair
        //              r_cmpl(l_cmpl(A)) == l_cmpl(r_cmpl(A)) == A, or from the other
        //              school's (A * I_4ds) * inv(I_4ds) == A.
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

        // --- school-1 r_dual vs school-2 A*I_4ds: identical at every grade in sta ---
        fmt::println("   FORWARD: ours r_dual(A) vs theirs A*I_4ds  (coincide in sta!):");
        rule();
        fmt::println("   {:^5} | {:<22} | {:<24} | {:<24} | {}", "grade", "primal A",
                     "ours: r_dual(A)", "theirs: A*I_4ds", "match");
        rule();
        cmp("0", s, r_dual(s), s * I_4ds, "equal");
        cmp("1", v, r_dual(v), v * I_4ds, "equal");
        cmp("2", B, r_dual(B), B * I_4ds, "equal");
        cmp("3", t, r_dual(t), t * I_4ds, "equal");
        cmp("4", ps, r_dual(ps), ps * I_4ds, "equal");
        fmt::println("");

        CHECK(s * I_4ds == r_dual(s));
        CHECK(v * I_4ds == r_dual(v));
        CHECK(B * I_4ds == r_dual(B));
        CHECK(t * I_4ds == r_dual(t));
        CHECK(ps * I_4ds == r_dual(ps));

        // basis spot-checks across grades: A*I_4ds == r_dual(A)
        CHECK(one_4ds * I_4ds == r_dual(one_4ds));
        CHECK(g1_4ds * I_4ds == r_dual(g1_4ds));
        CHECK(g14_4ds * I_4ds == r_dual(g14_4ds));
        CHECK(g234_4ds * I_4ds == r_dual(g234_4ds));
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
        // full multivector: M * I_4ds coincides with r_dual(M) exactly (no grade flips,
        // unlike the Euclidean cases).
        ////////////////////////////////////////////////////////////////////////////////

        auto const M = mvec4ds{s, v, B, t, ps};

        fmt::println("   full multivector  M = {}:", M);
        fmt::println("     {:<20} = {}", "ours : r_dual(M)", sf(r_dual(M)));
        fmt::println("     {:<20} = {}", "them : M * I_4ds", sf(M * I_4ds));
        fmt::println("     {:<20} = {}", "ours : l_dual(M)", sf(l_dual(M)));
        fmt::println("     (r_dual == theirs exactly; l_dual flips the odd-grade-image "
                     "parts)");
        fmt::println("");

        CHECK(M * I_4ds == r_dual(M));

        ////////////////////////////////////////////////////////////////////////////////
        // inverse dualization and the Minkowski-metric surprise:
        //
        //   - the two metric-dual compositions AGREE with each other (as in any even
        //     dimension) ...
        //         l_dual(r_dual(A)) == r_dual(l_dual(A))
        //   - ... but in the mixed metric they do NOT recover A: they return det(G)=-1
        //     signed by grade (= +A for the odd grades 1,3; = -A for the even 0,2,4).
        //   - genuine recovery: the metric-free complement pair, or the other school's
        //     multiplication by inv(I_4ds).
        ////////////////////////////////////////////////////////////////////////////////

        fmt::println("   pseudoscalar / round-trip facts:");
        fmt::println("     {:<26} = {}", "I_4ds * I_4ds", sf(I_4ds * I_4ds));
        fmt::println("     {:<26} = {}", "I_4ds * inv(I_4ds)", sf(I_4ds * I_inv));
        fmt::println("     {:<26} = {}  (= -s, NOT s!)", "l_dual(r_dual(s))",
                     sf(l_dual(r_dual(s))));
        fmt::println("     {:<26} = {}  (= v)", "l_dual(r_dual(v))",
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

        // ... but they return det(G)-signed A: -A for even grades, +A for odd grades
        CHECK(l_dual(r_dual(s)) == -s);
        CHECK(l_dual(r_dual(v)) == v);
        CHECK(l_dual(r_dual(B)) == -B);
        CHECK(l_dual(r_dual(t)) == t);
        CHECK(l_dual(r_dual(ps)) == -ps);

        ////////////////////////////////////////////////////////////////////////////////
        // CAN THE MISMATCH BE HEALED BY INCLUDING THE DETERMINANT?
        //
        // For a grade-k blade the composition evaluates in closed form to
        //     l_dual(r_dual(A)) = sigma(k)*sigma(n-k) * det(G) * A,   sigma(j)=(-1)^(j(j-1)/2)
        // and for n = 4 the combinatorial factor sigma(k)*sigma(n-k) == (-1)^k, so
        //     l_dual(r_dual(A)) == det(G) * gr_inv(A)        with gr_inv(A) = (-1)^k A.
        //
        // => the determinant ALONE cannot heal it: multiplying by det(G) leaves the
        //    grade-involution sign (-1)^k behind (det(G)*comp == gr_inv(A) != A). The
        //    fix needs BOTH the determinant AND the grade involution:
        //        A == det(G) * gr_inv( l_dual(r_dual(A)) ).
        //
        // METRIC-INDEPENDENCE (the real test of a "solution"): det(G) = (-1)^q where q
        // is the number of negative directions, so BOTH physically-equivalent STA
        // signatures give the same value -- (-,-,-,+) has q=3 and (+,+,+,-) has q=1, and
        // (-1)^3 == (-1)^1 == -1. sigma(k)*sigma(n-k) is purely combinatorial (no metric
        // at all). Hence the closed form and the healing formula are identical for both
        // signatures -- this is a genuine metric-independent correction, not a fix tuned
        // to G(1,3,0). (I_4ds*I_4ds == det(G) here because sigma(n)==+1 for n=4.)
        ////////////////////////////////////////////////////////////////////////////////

        value_t const detG = value_t(I_4ds * I_4ds); // == det(G) == -1 for both signatures

        fmt::println("   healing the metric-dual round trip:");
        fmt::println("     det(G) = I_4ds*I_4ds = {}", detG);
        fmt::println("     {:<34} = {}", "l_dual(r_dual(B))", sf(l_dual(r_dual(B))));
        fmt::println("     {:<34} = {}", "det(G) * gr_inv(B)", sf(detG * gr_inv(B)));
        fmt::println("     {:<34} = {}", "det(G)*gr_inv(l_dual(r_dual(B)))",
                     sf(detG * gr_inv(l_dual(r_dual(B)))));
        fmt::println("");

        // closed form: l_dual(r_dual(A)) == det(G) * gr_inv(A) at every grade
        CHECK(l_dual(r_dual(s)) == detG * gr_inv(s));
        CHECK(l_dual(r_dual(v)) == detG * gr_inv(v));
        CHECK(l_dual(r_dual(B)) == detG * gr_inv(B));
        CHECK(l_dual(r_dual(t)) == detG * gr_inv(t));
        CHECK(l_dual(r_dual(ps)) == detG * gr_inv(ps));

        // determinant ALONE is not enough: det(G) * composition == gr_inv(A) != A
        CHECK(detG * l_dual(r_dual(v)) == gr_inv(v));
        CHECK(detG * l_dual(r_dual(v)) != v); // (gr_inv flips the odd-grade vector)

        // the healing: det(G) AND grade involution recover A for every grade
        CHECK(detG * gr_inv(l_dual(r_dual(s))) == s);
        CHECK(detG * gr_inv(l_dual(r_dual(v))) == v);
        CHECK(detG * gr_inv(l_dual(r_dual(B))) == B);
        CHECK(detG * gr_inv(l_dual(r_dual(t))) == t);
        CHECK(detG * gr_inv(l_dual(r_dual(ps))) == ps);

        // equivalently: det(G)*gr_inv is an involution, so applying the metric-dual pair
        // TWICE recovers A (no determinant or gr_inv needed at all)
        CHECK(l_dual(r_dual(l_dual(r_dual(s)))) == s);
        CHECK(l_dual(r_dual(l_dual(r_dual(v)))) == v);
        CHECK(l_dual(r_dual(l_dual(r_dual(B)))) == B);
        CHECK(l_dual(r_dual(l_dual(r_dual(t)))) == t);
        CHECK(l_dual(r_dual(l_dual(r_dual(ps)))) == ps);

        // genuine recovery (1): the metric-free complement pair round-trips for all grades
        CHECK(r_cmpl(l_cmpl(s)) == s);   CHECK(l_cmpl(r_cmpl(s)) == s);
        CHECK(r_cmpl(l_cmpl(v)) == v);   CHECK(l_cmpl(r_cmpl(v)) == v);
        CHECK(r_cmpl(l_cmpl(B)) == B);   CHECK(l_cmpl(r_cmpl(B)) == B);
        CHECK(r_cmpl(l_cmpl(t)) == t);   CHECK(l_cmpl(r_cmpl(t)) == t);
        CHECK(r_cmpl(l_cmpl(ps)) == ps); CHECK(l_cmpl(r_cmpl(ps)) == ps);

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
        fmt::println("   {:<11} | {:<6} | {:<24} | {}", "input", "school",
                     "forward", "backward (-> input)");
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

} // STA 3D Tests
