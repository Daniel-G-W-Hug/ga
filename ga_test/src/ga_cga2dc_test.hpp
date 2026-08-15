// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <array>
#include <random>

#include "fmt/format.h" // formatting

// include functions to be tested
#include "ga/ga_cga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::cga; // use specific operations of CGA (Conformal Algebra)


/////////////////////////////////////////////////////////////////////////////////////////
// Conformal Geometric Algebra of 2d Euclidean space: cga2dc
// (G(3,1,0) after diagonalization of the native null-basis metric)
/////////////////////////////////////////////////////////////////////////////////////////

namespace {

// embed a Euclidean point (x, y) as a normalized round point on the null cone
// (w = 1, u = 0.5*(x^2 + y^2)); q * q = 0 for every embedded point
vec2dc point2dc(value_t x, value_t y) { return vec2dc(x, y, 1.0, 0.5 * (x * x + y * y)); }

// random multivectors with small integer coefficients: all products and
// grade-wise sums stay exact in double, so tests can compare with operator==
std::mt19937 rnd_gen(20260815);

value_t rnd_int()
{
    static std::uniform_int_distribution<int> dist(-3, 3);
    return value_t(dist(rnd_gen));
}

mvec2dc rnd_mv()
{
    return mvec2dc(rnd_int(), rnd_int(), rnd_int(), rnd_int(), rnd_int(), rnd_int(),
                   rnd_int(), rnd_int(), rnd_int(), rnd_int(), rnd_int(), rnd_int(),
                   rnd_int(), rnd_int(), rnd_int(), rnd_int());
}

vec2dc rnd_vec() { return vec2dc(rnd_int(), rnd_int(), rnd_int(), rnd_int()); }

bivec2dc rnd_bivec()
{
    return bivec2dc(rnd_int(), rnd_int(), rnd_int(), rnd_int(), rnd_int(), rnd_int());
}

trivec2dc rnd_trivec() { return trivec2dc(rnd_int(), rnd_int(), rnd_int(), rnd_int()); }

} // anonymous namespace


TEST_SUITE("CGA 2dc Tests")
{

    TEST_CASE("G<3,1,0>: cga2dc")
    {
        fmt::println("G<3,1,0>: cga2dc");
        // Conformal algebra of 2d Euclidean space; the native null-basis
        // metric (e3.e3 = e4.e4 = 0, e3.e4 = -1) diagonalizes to G(3,1,0)
        const algebra<3, 1, 0> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 1);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   cga2dc: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   cga2dc: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    TEST_CASE("cga2dc: type layer round trips")
    {
        fmt::println("cga2dc: type layer round trips");

        auto v = rnd_vec();
        auto B = rnd_bivec();
        auto t = rnd_trivec();
        auto s = scalar2dc(3.0);
        auto ps = pscalar2dc(-2.0);

        // grades of the graded types
        CHECK(gr(s) == 0);
        CHECK(gr(v) == 1);
        CHECK(gr(B) == 2);
        CHECK(gr(t) == 3);
        CHECK(gr(ps) == 4);
        CHECK(rgr(s) == 4);
        CHECK(rgr(v) == 3);
        CHECK(rgr(B) == 2);
        CHECK(rgr(t) == 1);
        CHECK(rgr(ps) == 0);

        // graded-part constructors and grade selection round-trip
        mvec2dc M(s, v, B, t, ps);
        CHECK(gr0(M) == s);
        CHECK(gr1(M) == v);
        CHECK(gr2(M) == B);
        CHECK(gr3(M) == t);
        CHECK(gr4(M) == ps);

        mvec2dc_e Me(s, B, ps);
        CHECK(gr0(Me) == s);
        CHECK(gr2(Me) == B);
        CHECK(gr4(Me) == ps);

        mvec2dc_u Mu(v, t);
        CHECK(gr1(Mu) == v);
        CHECK(gr3(Mu) == t);

        // even/odd embedding into the full multivector
        CHECK(mvec2dc(Me) == mvec2dc(s, vec2dc(), B, trivec2dc(), ps));
        CHECK(mvec2dc(Mu) == mvec2dc(scalar2dc(0.0), v, bivec2dc(), t, pscalar2dc(0.0)));
    }

    TEST_CASE("cga2dc: involutions")
    {
        fmt::println("cga2dc: involutions");

        // grade-wise sign patterns (dimension 4):
        // gr_inv: + - + - +, rev: + + - - +, rrev: + - - + +, conj: + - - + +
        CHECK(gr_inv(scalar2dc(1.0)) == scalar2dc(1.0));
        CHECK(gr_inv(e1_2dc) == -e1_2dc);
        CHECK(gr_inv(e12_2dc) == e12_2dc);
        CHECK(gr_inv(e314_2dc) == -e314_2dc);
        CHECK(gr_inv(I_2dc) == I_2dc);

        CHECK(rev(scalar2dc(1.0)) == scalar2dc(1.0));
        CHECK(rev(e1_2dc) == e1_2dc);
        CHECK(rev(e12_2dc) == -e12_2dc);
        CHECK(rev(e314_2dc) == -e314_2dc);
        CHECK(rev(I_2dc) == I_2dc);

        CHECK(rrev(scalar2dc(1.0)) == scalar2dc(1.0));
        CHECK(rrev(e1_2dc) == -e1_2dc);
        CHECK(rrev(e12_2dc) == -e12_2dc);
        CHECK(rrev(e314_2dc) == e314_2dc);
        CHECK(rrev(I_2dc) == I_2dc);

        CHECK(conj(scalar2dc(1.0)) == scalar2dc(1.0));
        CHECK(conj(e1_2dc) == -e1_2dc);
        CHECK(conj(e12_2dc) == -e12_2dc);
        CHECK(conj(e314_2dc) == e314_2dc);
        CHECK(conj(I_2dc) == I_2dc);

        // structural identities on random multivectors
        for (int i = 0; i < 3; ++i) {
            auto M = rnd_mv();
            CHECK(rev(rev(M)) == M);
            CHECK(gr_inv(gr_inv(M)) == M);
            CHECK(conj(M) == gr_inv(rev(M)));
            CHECK(rrev(M) == l_cmpl(rev(r_cmpl(M))));
        }
    }

    TEST_CASE("cga2dc: complements")
    {
        fmt::println("cga2dc: complements");

        // defining property per unit basis blade: l_cmpl(u) ^ u = u ^ r_cmpl(u) = I
        CHECK(wdg(l_cmpl(scalar2dc(1.0)), scalar2dc(1.0)) == I_2dc);
        CHECK(wdg(scalar2dc(1.0), r_cmpl(scalar2dc(1.0))) == I_2dc);
        CHECK(wdg(l_cmpl(e1_2dc), e1_2dc) == I_2dc);
        CHECK(wdg(e1_2dc, r_cmpl(e1_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e2_2dc), e2_2dc) == I_2dc);
        CHECK(wdg(e2_2dc, r_cmpl(e2_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e3_2dc), e3_2dc) == I_2dc);
        CHECK(wdg(e3_2dc, r_cmpl(e3_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e4_2dc), e4_2dc) == I_2dc);
        CHECK(wdg(e4_2dc, r_cmpl(e4_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e31_2dc), e31_2dc) == I_2dc);
        CHECK(wdg(e31_2dc, r_cmpl(e31_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e32_2dc), e32_2dc) == I_2dc);
        CHECK(wdg(e32_2dc, r_cmpl(e32_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e12_2dc), e12_2dc) == I_2dc);
        CHECK(wdg(e12_2dc, r_cmpl(e12_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e14_2dc), e14_2dc) == I_2dc);
        CHECK(wdg(e14_2dc, r_cmpl(e14_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e24_2dc), e24_2dc) == I_2dc);
        CHECK(wdg(e24_2dc, r_cmpl(e24_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e34_2dc), e34_2dc) == I_2dc);
        CHECK(wdg(e34_2dc, r_cmpl(e34_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e314_2dc), e314_2dc) == I_2dc);
        CHECK(wdg(e314_2dc, r_cmpl(e314_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e324_2dc), e324_2dc) == I_2dc);
        CHECK(wdg(e324_2dc, r_cmpl(e324_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e124_2dc), e124_2dc) == I_2dc);
        CHECK(wdg(e124_2dc, r_cmpl(e124_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(e321_2dc), e321_2dc) == I_2dc);
        CHECK(wdg(e321_2dc, r_cmpl(e321_2dc)) == I_2dc);
        CHECK(wdg(l_cmpl(I_2dc), I_2dc) == I_2dc);
        CHECK(wdg(I_2dc, r_cmpl(I_2dc)) == I_2dc);

        // round trips on random multivectors
        for (int i = 0; i < 3; ++i) {
            auto M = rnd_mv();
            CHECK(l_cmpl(r_cmpl(M)) == M);
            CHECK(r_cmpl(l_cmpl(M)) == M);
        }
    }

    TEST_CASE("cga2dc: dual/antidual transcription gate")
    {
        fmt::println("cga2dc: dual/antidual transcription gate");

        // The extended metric G is a signed permutation (hand-pinned from the
        // literature-reviewed rule reference, reference_output/
        // ga_prdxpr_rulegen_cga2dc.txt):
        //
        //     G(1) = 1, G(e1) = e1, G(e2) = e2, G(e3) = -e4, G(e4) = -e3,
        //     G(e31) = e14, G(e32) = e24, G(e12) = e12,
        //     G(e14) = e31, G(e24) = e32, G(e34) = -e34,
        //     G(e314) = -e314, G(e324) = -e324, G(e124) = e321, G(e321) = e124,
        //     G(e1234) = -e1234
        //
        // gate: dual(e) == cmpl(G(e)) per unit basis blade (both l/r variants)
        CHECK(l_dual(scalar2dc(1.0)) == l_cmpl(scalar2dc(1.0)));
        CHECK(l_dual(e1_2dc) == l_cmpl(e1_2dc));
        CHECK(l_dual(e2_2dc) == l_cmpl(e2_2dc));
        CHECK(l_dual(e3_2dc) == l_cmpl(vec2dc(-e4_2dc)));
        CHECK(l_dual(e4_2dc) == l_cmpl(vec2dc(-e3_2dc)));
        CHECK(l_dual(e31_2dc) == l_cmpl(e14_2dc));
        CHECK(l_dual(e32_2dc) == l_cmpl(e24_2dc));
        CHECK(l_dual(e12_2dc) == l_cmpl(e12_2dc));
        CHECK(l_dual(e14_2dc) == l_cmpl(e31_2dc));
        CHECK(l_dual(e24_2dc) == l_cmpl(e32_2dc));
        CHECK(l_dual(e34_2dc) == l_cmpl(bivec2dc(-e34_2dc)));
        CHECK(l_dual(e314_2dc) == l_cmpl(trivec2dc(-e314_2dc)));
        CHECK(l_dual(e324_2dc) == l_cmpl(trivec2dc(-e324_2dc)));
        CHECK(l_dual(e124_2dc) == l_cmpl(e321_2dc));
        CHECK(l_dual(e321_2dc) == l_cmpl(e124_2dc));
        CHECK(l_dual(I_2dc) == l_cmpl(pscalar2dc(-1.0)));

        CHECK(r_dual(scalar2dc(1.0)) == r_cmpl(scalar2dc(1.0)));
        CHECK(r_dual(e1_2dc) == r_cmpl(e1_2dc));
        CHECK(r_dual(e2_2dc) == r_cmpl(e2_2dc));
        CHECK(r_dual(e3_2dc) == r_cmpl(vec2dc(-e4_2dc)));
        CHECK(r_dual(e4_2dc) == r_cmpl(vec2dc(-e3_2dc)));
        CHECK(r_dual(e31_2dc) == r_cmpl(e14_2dc));
        CHECK(r_dual(e32_2dc) == r_cmpl(e24_2dc));
        CHECK(r_dual(e12_2dc) == r_cmpl(e12_2dc));
        CHECK(r_dual(e14_2dc) == r_cmpl(e31_2dc));
        CHECK(r_dual(e24_2dc) == r_cmpl(e32_2dc));
        CHECK(r_dual(e34_2dc) == r_cmpl(bivec2dc(-e34_2dc)));
        CHECK(r_dual(e314_2dc) == r_cmpl(trivec2dc(-e314_2dc)));
        CHECK(r_dual(e324_2dc) == r_cmpl(trivec2dc(-e324_2dc)));
        CHECK(r_dual(e124_2dc) == r_cmpl(e321_2dc));
        CHECK(r_dual(e321_2dc) == r_cmpl(e124_2dc));
        CHECK(r_dual(I_2dc) == r_cmpl(pscalar2dc(-1.0)));

        // direct value pins against the same reference (transcription gate for
        // the complement side as well)
        CHECK(l_dual(scalar2dc(1.0)) == I_2dc);
        CHECK(l_dual(e1_2dc) == e324_2dc);
        CHECK(l_dual(e2_2dc) == -e314_2dc);
        CHECK(l_dual(e3_2dc) == e321_2dc);
        CHECK(l_dual(e4_2dc) == e124_2dc);
        CHECK(l_dual(e31_2dc) == -e32_2dc);
        CHECK(l_dual(e32_2dc) == e31_2dc);
        CHECK(l_dual(e12_2dc) == e34_2dc);
        CHECK(l_dual(e14_2dc) == e24_2dc);
        CHECK(l_dual(e24_2dc) == -e14_2dc);
        CHECK(l_dual(e34_2dc) == -e12_2dc);
        CHECK(l_dual(e314_2dc) == -e2_2dc);
        CHECK(l_dual(e324_2dc) == e1_2dc);
        CHECK(l_dual(e124_2dc) == e4_2dc);
        CHECK(l_dual(e321_2dc) == e3_2dc);
        CHECK(l_dual(I_2dc) == scalar2dc(-1.0));

        CHECK(r_dual(scalar2dc(1.0)) == I_2dc);
        CHECK(r_dual(e1_2dc) == -e324_2dc);
        CHECK(r_dual(e2_2dc) == e314_2dc);
        CHECK(r_dual(e3_2dc) == -e321_2dc);
        CHECK(r_dual(e4_2dc) == -e124_2dc);
        CHECK(r_dual(e31_2dc) == -e32_2dc);
        CHECK(r_dual(e32_2dc) == e31_2dc);
        CHECK(r_dual(e12_2dc) == e34_2dc);
        CHECK(r_dual(e14_2dc) == e24_2dc);
        CHECK(r_dual(e24_2dc) == -e14_2dc);
        CHECK(r_dual(e34_2dc) == -e12_2dc);
        CHECK(r_dual(e314_2dc) == e2_2dc);
        CHECK(r_dual(e324_2dc) == -e1_2dc);
        CHECK(r_dual(e124_2dc) == -e4_2dc);
        CHECK(r_dual(e321_2dc) == -e3_2dc);
        CHECK(r_dual(I_2dc) == scalar2dc(-1.0));

        // antidual == negated dual (the anti-exomorphism equals -G for cga2dc;
        // the identity is verified at rule-generation time, pinned here per
        // unit basis blade)
        auto const check_antidual = [](auto const& e) {
            CHECK(l_antidual(e) == -l_dual(e));
            CHECK(r_antidual(e) == -r_dual(e));
        };
        check_antidual(scalar2dc(1.0));
        check_antidual(e1_2dc);
        check_antidual(e2_2dc);
        check_antidual(e3_2dc);
        check_antidual(e4_2dc);
        check_antidual(e31_2dc);
        check_antidual(e32_2dc);
        check_antidual(e12_2dc);
        check_antidual(e14_2dc);
        check_antidual(e24_2dc);
        check_antidual(e34_2dc);
        check_antidual(e314_2dc);
        check_antidual(e324_2dc);
        check_antidual(e124_2dc);
        check_antidual(e321_2dc);
        check_antidual(I_2dc);

        // aggregate duals delegate grade-wise
        for (int i = 0; i < 3; ++i) {
            auto M = rnd_mv();
            auto D = l_dual(M);
            CHECK(gr1(D) == l_dual(gr3(M)));
            CHECK(gr3(D) == l_dual(gr1(M)));
            CHECK(gr2(D) == l_dual(gr2(M)));
        }
    }

    TEST_CASE("cga2dc: dot and the extended metric")
    {
        fmt::println("cga2dc: dot and the extended metric");

        // vector level: the native null-pair metric
        CHECK(dot(e1_2dc, e1_2dc) == scalar2dc(1.0));
        CHECK(dot(e2_2dc, e2_2dc) == scalar2dc(1.0));
        CHECK(dot(e3_2dc, e3_2dc) == scalar2dc(0.0)); // null vector
        CHECK(dot(e4_2dc, e4_2dc) == scalar2dc(0.0)); // null vector
        CHECK(dot(e3_2dc, e4_2dc) == scalar2dc(-1.0));
        CHECK(dot(e1_2dc, e2_2dc) == scalar2dc(0.0));

        // blade level: G is a signed permutation -- w-blades pair with u-blades
        CHECK(dot(e31_2dc, e31_2dc) == scalar2dc(0.0));
        CHECK(dot(e14_2dc, e14_2dc) == scalar2dc(0.0));
        CHECK(dot(e31_2dc, e14_2dc) == scalar2dc(1.0));
        CHECK(dot(e32_2dc, e24_2dc) == scalar2dc(1.0));
        CHECK(dot(e12_2dc, e12_2dc) == scalar2dc(1.0));
        CHECK(dot(e34_2dc, e34_2dc) == scalar2dc(-1.0));
        CHECK(dot(e314_2dc, e314_2dc) == scalar2dc(-1.0));
        CHECK(dot(e324_2dc, e324_2dc) == scalar2dc(-1.0));
        CHECK(dot(e124_2dc, e124_2dc) == scalar2dc(0.0));
        CHECK(dot(e124_2dc, e321_2dc) == scalar2dc(1.0));
        CHECK(dot(I_2dc, I_2dc) == scalar2dc(-1.0));

        // multivector dot == sum of the graded dots
        for (int i = 0; i < 3; ++i) {
            auto A = rnd_mv();
            auto B = rnd_mv();
            auto expected =
                scalar2dc(value_t(dot(gr0(A), gr0(B))) + value_t(dot(gr1(A), gr1(B))) +
                          value_t(dot(gr2(A), gr2(B))) + value_t(dot(gr3(A), gr3(B))) +
                          value_t(dot(gr4(A), gr4(B))));
            CHECK(dot(A, B) == expected);
        }
    }

    TEST_CASE("cga2dc: geometric product identities")
    {
        fmt::println("cga2dc: geometric product identities");

        // the null pair makes single basis products multi-term
        auto m34 = e3_2dc * e4_2dc; // -1 + e34
        CHECK(gr0(m34) == scalar2dc(-1.0));
        CHECK(gr2(m34) == e34_2dc);
        CHECK(gr4(m34) == pscalar2dc(0.0));
        CHECK(e34_2dc * e34_2dc == mvec2dc_e(scalar2dc(1.0)));
        CHECK(I_2dc * I_2dc == scalar2dc(-1.0));

        // vector anticommutators reproduce the vector Gram matrix: v_i v_j +
        // v_j v_i = 2 G[i,j]
        std::array<vec2dc, 4> const basis_vec{e1_2dc, e2_2dc, e3_2dc, e4_2dc};
        std::array<std::array<value_t, 4>, 4> const G{{{1.0, 0.0, 0.0, 0.0},
                                                       {0.0, 1.0, 0.0, 0.0},
                                                       {0.0, 0.0, 0.0, -1.0},
                                                       {0.0, 0.0, -1.0, 0.0}}};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                auto acmt = basis_vec[i] * basis_vec[j] + basis_vec[j] * basis_vec[i];
                CHECK(gr0(acmt) == scalar2dc(2.0 * G[i][j]));
                CHECK(gr2(acmt) == bivec2dc());
                CHECK(gr4(acmt) == pscalar2dc(0.0));
            }
        }

        // gpr(v1,v2) = dot(v1,v2) + wdg(v1,v2) for vectors
        for (int i = 0; i < 3; ++i) {
            auto v1 = rnd_vec();
            auto v2 = rnd_vec();
            CHECK(v1 * v2 == mvec2dc_e(dot(v1, v2), wdg(v1, v2), pscalar2dc(0.0)));
        }

        // associativity (exact for integer coefficients)
        for (int i = 0; i < 3; ++i) {
            auto A = rnd_mv();
            auto B = rnd_mv();
            auto C = rnd_mv();
            CHECK((A * B) * C == A * (B * C));
        }

        // commutator product
        CHECK(cmt(e3_2dc, e4_2dc) == e34_2dc);
        for (int i = 0; i < 3; ++i) {
            auto A = rnd_mv();
            auto B = rnd_mv();
            CHECK(cmt(A, B) == 0.5 * (A * B - B * A));
        }
    }

    TEST_CASE("cga2dc: null-point embedding")
    {
        fmt::println("cga2dc: null-point embedding");

        // every embedded point lies on the null cone: q * q = 0
        CHECK(point2dc(3.0, 4.0) * point2dc(3.0, 4.0) == mvec2dc_e(scalar2dc(0.0)));
        CHECK(point2dc(-1.0, 2.0) * point2dc(-1.0, 2.0) == mvec2dc_e(scalar2dc(0.0)));
        CHECK(point2dc(0.0, 0.0) * point2dc(0.0, 0.0) == mvec2dc_e(scalar2dc(0.0)));

        // the conformal distance identity for normalized points:
        // dot(p, q) = -0.5 * |p - q|^2
        auto const check_dist = [](value_t px, value_t py, value_t qx, value_t qy) {
            auto p = point2dc(px, py);
            auto q = point2dc(qx, qy);
            value_t const d_sq = (px - qx) * (px - qx) + (py - qy) * (py - qy);
            CHECK(dot(p, q) == scalar2dc(-0.5 * d_sq));
        };
        check_dist(3.0, 4.0, 1.0, 2.0);
        check_dist(0.0, 0.0, 5.0, 0.0);
        check_dist(-2.0, 1.0, 2.0, -1.0);
    }

    TEST_CASE("cga2dc: regressive product family")
    {
        fmt::println("cga2dc: regressive product family");

        // the pseudoscalar is the identity element of rgpr
        for (int i = 0; i < 3; ++i) {
            auto M = rnd_mv();
            CHECK(rgpr(M, I_2dc) == M);
            CHECK(rgpr(I_2dc, M) == M);
        }

        // rgpr(1,1) = l_cmpl(gpr(r_cmpl(1), r_cmpl(1))) = l_cmpl(I * I) = -e1234
        // (non-degenerate metric: I * I = -1; this case is identically zero in
        // the degenerate pga algebras -- do not transplant the zero here)
        CHECK(rgpr(scalar2dc(1.0), scalar2dc(1.0)) == pscalar2dc(-1.0));
        CHECK(rdot(scalar2dc(1.0), scalar2dc(1.0)) == pscalar2dc(-1.0));

        // complement conjugation identities on random multivectors
        for (int i = 0; i < 3; ++i) {
            auto A = rnd_mv();
            auto B = rnd_mv();
            CHECK(rgpr(A, B) == l_cmpl(r_cmpl(A) * r_cmpl(B)));
            CHECK(rwdg(A, B) == l_cmpl(wdg(r_cmpl(A), r_cmpl(B))));
            CHECK(rdot(A, B) == l_cmpl(dot(r_cmpl(A), r_cmpl(B))));
            CHECK(rcmt(A, B) == 0.5 * (rgpr(A, B) - rgpr(B, A)));
        }
    }

    TEST_CASE("cga2dc: contractions and expansions")
    {
        fmt::println("cga2dc: contractions and expansions");

        // contraction = dual + regressive wedge (single pair, non-degenerate)
        for (int i = 0; i < 3; ++i) {
            auto A = rnd_mv();
            auto B = rnd_mv();
            CHECK((A << B) == rwdg(l_dual(A), B));
            CHECK((A >> B) == rwdg(A, r_dual(B)));
        }

        // on equal grades the contraction reduces to the metric inner product
        for (int i = 0; i < 3; ++i) {
            auto v1 = rnd_vec();
            auto v2 = rnd_vec();
            auto B1 = rnd_bivec();
            auto B2 = rnd_bivec();
            auto t1 = rnd_trivec();
            auto t2 = rnd_trivec();
            CHECK((v1 << v2) == dot(v1, v2));
            CHECK((v1 >> v2) == dot(v1, v2));
            CHECK((B1 << B2) == dot(B1, B2));
            CHECK((B1 >> B2) == dot(B1, B2));
            CHECK((t1 << t2) == dot(t1, t2));
            CHECK((t1 >> t2) == dot(t1, t2));
        }

        // expansions: wdg-based duals of the contractions, with a hand-pinned
        // value each: l_dual(e12) = e34, wdg(e34, e1) = e3^e4^e1 = -e314;
        // r_dual(e12) = e34, wdg(e1, e34) = e1^e3^e4 = -e314
        CHECK(l_expand2dc(e12_2dc, e1_2dc) == -e314_2dc);
        CHECK(r_expand2dc(e1_2dc, e12_2dc) == -e314_2dc);
        for (int i = 0; i < 3; ++i) {
            auto B = rnd_bivec();
            auto v = rnd_vec();
            CHECK(l_expand2dc(B, v) == wdg(l_dual(B), v));
            CHECK(r_expand2dc(v, B) == wdg(v, r_dual(B)));
        }
    }

    TEST_CASE("cga2dc: fmt printing")
    {
        fmt::println("cga2dc: fmt printing");
        vec2dc v(1.0, 2.0, 3.0, 4.0);
        auto const s = fmt::format("{}", v);
        CHECK(s.size() > 0);
        fmt::println("   vec2dc: {}", v);
        fmt::println("");
    }

} // TEST_SUITE("CGA 2dc Tests")
