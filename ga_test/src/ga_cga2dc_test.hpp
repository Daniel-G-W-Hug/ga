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

    TEST_CASE("cga2dc: inv and rinv")
    {
        fmt::println("cga2dc: inv and rinv");

        // pinned values (non-degenerate metric: BOTH pseudoscalar inv and
        // scalar rinv exist, unlike in the degenerate pga algebras)
        CHECK(inv(scalar2dc(2.0)) == scalar2dc(0.5));
        CHECK(inv(I_2dc) == pscalar2dc(-1.0));          // I * I = -1
        CHECK(rinv(scalar2dc(1.0)) == scalar2dc(-1.0)); // rgpr(1, -1) = I
        CHECK(is_close(rgpr(scalar2dc(1.0), rinv(scalar2dc(1.0))), pscalar2dc(1.0)));

        // null elements are not invertible -- including every embedded point
        CHECK_THROWS(inv(e3_2dc));
        CHECK_THROWS(inv(e4_2dc));
        CHECK_THROWS(inv(point2dc(3.0, 4.0)));

        // gate: u * inv(u) = 1 for random (generically non-null) elements
        auto const one_e = mvec2dc_e(scalar2dc(1.0));
        auto const one_m = mvec2dc(scalar2dc(1.0));
        auto const I_e = mvec2dc_e(pscalar2dc(1.0));
        auto const I_m = mvec2dc(pscalar2dc(1.0));
        for (int i = 0; i < 3; ++i) {
            auto v = rnd_vec();
            auto B = rnd_bivec();
            auto t = rnd_trivec();
            auto Me = mvec2dc_e(rnd_bivec()) + mvec2dc_e(scalar2dc(rnd_int()));
            auto Mu = mvec2dc_u(rnd_vec(), rnd_trivec());
            auto M = rnd_mv();
            CHECK(is_close(v * inv(v), one_e));
            CHECK(is_close(inv(v) * v, one_e));
            CHECK(is_close(B * inv(B), one_e));
            CHECK(is_close(inv(B) * B, one_e));
            CHECK(is_close(t * inv(t), one_e));
            CHECK(is_close(inv(t) * t, one_e));
            CHECK(is_close(Me * inv(Me), one_e));
            CHECK(is_close(Mu * inv(Mu), one_e));
            CHECK(is_close(M * inv(M), one_m));
            CHECK(is_close(inv(M) * M, one_m));

            // gate: rgpr(u, rinv(u)) = I (the rgpr identity element)
            CHECK(is_close(rgpr(v, rinv(v)), I_e));
            CHECK(is_close(rgpr(rinv(v), v), I_e));
            CHECK(is_close(rgpr(B, rinv(B)), I_e));
            CHECK(is_close(rgpr(t, rinv(t)), I_e));
            CHECK(is_close(rgpr(M, rinv(M)), I_m));
            CHECK(is_close(rgpr(rinv(M), M), I_m));
        }
    }

    TEST_CASE("cga2dc: congruence and closeness")
    {
        fmt::println("cga2dc: congruence and closeness");

        auto v = vec2dc(1.0, -2.0, 3.0, 0.5);
        auto B = bivec2dc(1.0, 0.0, -2.0, 3.0, 0.0, 1.5);
        auto t = trivec2dc(-1.0, 2.0, 0.0, 4.0);

        // congruent: same subspace up to any non-zero scale (sign included)
        CHECK(is_congruent(v, vec2dc(2.0 * v)));
        CHECK(is_congruent(v, vec2dc(-3.0 * v)));
        CHECK(is_congruent(B, bivec2dc(-0.5 * B)));
        CHECK(is_congruent(t, trivec2dc(7.0 * t)));
        CHECK_FALSE(is_congruent(v, e1_2dc));
        CHECK_FALSE(is_congruent(B, e12_2dc));
        CHECK(is_congruent(scalar2dc(2.0), scalar2dc(-5.0)));
        CHECK(is_congruent(I_2dc, pscalar2dc(-3.0)));

        // embedded points scale-normalize to the same round point
        auto p = point2dc(3.0, 4.0);
        CHECK(is_congruent(p, vec2dc(2.5 * p)));

        // is_close: same value within a relative tolerance
        CHECK(is_close(v, v));
        CHECK_FALSE(is_close(v, vec2dc(2.0 * v)));
        CHECK(is_close(mvec2dc(v, t), mvec2dc(v, t)));
    }

    TEST_CASE("cga2dc: geometric objects, split and norms")
    {
        fmt::println("cga2dc: geometric objects, split and norms");

        // round point (1,2) r=1; null point embedding is the r=0 special case
        auto a = round_point2dc(1.0, 2.0, 1.0);
        CHECK(a == vec2dc(1.0, 2.0, 1.0, 3.0));
        CHECK(round_point2dc(3.0, 4.0, 0.0) == point2dc(3.0, 4.0));

        // circle (1,2) r=1: negative carrier weight, center readable directly
        auto c = circle2dc(1.0, 2.0, 1.0);
        CHECK(c == trivec2dc(2.0, -1.0, -2.0, -1.0));

        // a dipole's two points lie perpendicular to its normal: for center
        // (1,2), r=2, normal (0,1) the points are (-1,2) and (3,2); the join
        // of those embedded points is congruent to the constructed dipole
        auto d = dipole2dc(1.0, 2.0, 2.0, 0.0, 1.0);
        CHECK(is_congruent(d, wdg(point2dc(-1.0, 2.0), point2dc(3.0, 2.0))));

        // the flat point is the dipole with vanishing carrier-line part
        auto fp = flat_point2dc(3.0, 4.0);
        CHECK(round_weight_nrm_sq(fp) == 0.0);
        CHECK(round_bulk_nrm_sq(fp) == 0.0);

        // split: the four parts are disjoint and sum to the object
        auto B = rnd_bivec();
        auto t = rnd_trivec();
        auto v = rnd_vec();
        CHECK(round_bulk(B) + round_weight(B) + flat_bulk(B) + flat_weight(B) == B);
        CHECK(round_bulk(t) + round_weight(t) + flat_bulk(t) + flat_weight(t) == t);
        CHECK(round_bulk(v) + round_weight(v) + flat_bulk(v) + flat_weight(v) == v);

        // norm pins (reference-table values through the dictionaries):
        // round point: rb = |center dist|*w, rw = |w|, fb = |u-comp|, fw = 0
        CHECK(round_bulk_nrm_sq(a) == 5.0);
        CHECK(round_weight_nrm_sq(a) == 1.0);
        CHECK(flat_bulk_nrm_sq(a) == 9.0);
        CHECK(flat_weight_nrm_sq(a) == 0.0);
        CHECK(center_nrm_sq(a) == 5.0); // origin<->center distance squared
        // circle: rw = |w|, fw = center dist * |w|
        CHECK(round_weight_nrm_sq(c) == 1.0);
        CHECK(flat_weight_nrm_sq(c) == 5.0);
        CHECK(center_nrm_sq(c) == 5.0);

        // radius norms: the antidot square r^2 (signed), the dot square -r^2
        CHECK(radius_nrm_sq(a) == 1.0);
        CHECK(radius_nrm_sq(c) == 1.0);
        CHECK(is_close(scalar2dc(radius_sq(d)), scalar2dc(4.0)));
        CHECK(dot(c, c) == scalar2dc(-1.0));
        CHECK(radius_nrm_sq(point2dc(3.0, 4.0)) == 0.0); // null point

        // unitize: round weight norm becomes one, sign preserved
        auto c3 = trivec2dc(3.0 * c.x, 3.0 * c.y, 3.0 * c.z, 3.0 * c.w);
        CHECK(unitize(c3) == c);
        CHECK_THROWS(unitize(fp)); // flat objects have no round weight
    }

    TEST_CASE("cga2dc: carrier, cocarrier, attitude")
    {
        fmt::println("cga2dc: carrier, cocarrier, attitude");

        auto a = round_point2dc(1.0, 2.0, 1.0);
        auto c = circle2dc(1.0, 2.0, 1.0);
        auto d = dipole2dc(1.0, 2.0, 2.0, 0.0, 1.0);

        // carrier: lowest-dimensional flat containing the object. For a round
        // point: the flat point at the same position (congruent to it)
        CHECK(is_congruent(car(a), flat_point2dc(1.0, 2.0)));
        // for a dipole: the carrier line -- every point of the dipole lies on
        // it: join of carrier line with either end point vanishes
        CHECK(wdg(point2dc(-1.0, 2.0), car(d)) == pscalar2dc(0.0));
        CHECK(wdg(point2dc(3.0, 2.0), car(d)) == pscalar2dc(0.0));
        // for a circle: the whole space (pseudoscalar; the trailing-wedge
        // convention gives car(c) = -c.w * I -- odd-grade orientation flip
        // vs. the reference's leading-wedge 2D table)
        CHECK(car(c) == pscalar2dc(1.0));

        // cocarrier = carrier of the antidual; meet of cocarrier and carrier
        // is the object's (flat) center
        CHECK(is_congruent(rwdg(ccr(d), car(d)), flat_point2dc(1.0, 2.0)));

        // attitude: origin removed; unitized round point -> weight 1; dipole:
        // carrier direction (perpendicular to the normal) in the v-fields
        CHECK(att(a) == scalar2dc(1.0));
        auto ad = att(d);
        CHECK(is_congruent(vec2dc(ad.x, ad.y, 0.0, 0.0), vec2dc(2.0, 0.0, 0.0, 0.0)));
    }

    TEST_CASE("cga2dc: center, container, partner")
    {
        fmt::println("cga2dc: center, container, partner");

        auto a = round_point2dc(1.0, 2.0, 1.0);
        auto c = circle2dc(1.0, 2.0, 1.0);
        auto d = dipole2dc(1.0, 2.0, 2.0, 0.0, 1.0);

        // cen: the round point with the object's position AND radius
        CHECK(cen(a) == a); // unitized round point is its own center
        CHECK(cen(c) == a); // exact: pins the ccr orientation convention, not
                            // just the subspace
        CHECK(is_congruent(cen(d), round_point2dc(1.0, 2.0, 2.0)));
        // center coordinates: normalize by the origin component
        auto cd = cen(d);
        CHECK(is_close(vec2dc(cd.x / cd.z, cd.y / cd.z, 1.0, cd.w / cd.z),
                       round_point2dc(1.0, 2.0, 2.0)));

        // con: the smallest circle containing the object
        CHECK(is_congruent(con(a), c)); // container of a round point: same r
        CHECK(is_congruent(con(d), circle2dc(1.0, 2.0, 2.0)));
        CHECK(is_congruent(con(c), c)); // a circle is its own container

        // par: same center and carrier, radius squared negated,
        // and always orthogonal to the object (u . par(u) = 0)
        CHECK(par(a) == vec2dc(1.0, 2.0, 1.0, 2.0)); // exact: pins the par sign
        CHECK(is_close(scalar2dc(radius_sq(par(a))), scalar2dc(-1.0)));
        CHECK(is_close(scalar2dc(radius_sq(par(d))), scalar2dc(-4.0)));
        CHECK(is_close(scalar2dc(radius_sq(par(c))), scalar2dc(-1.0)));
        CHECK(is_congruent(cen(par(d)), round_point2dc(1.0, 2.0, 2.0)) ==
              false); // radius differs -> cen not congruent ...
        CHECK(is_close(vec2dc(cen(par(d)).x / cen(par(d)).z,
                              cen(par(d)).y / cen(par(d)).z, 1.0, 0.0),
                       vec2dc(1.0, 2.0, 1.0, 0.0))); // ... but same position
        CHECK(dot(a, par(a)) == scalar2dc(0.0));
        CHECK(dot(d, par(d)) == scalar2dc(0.0));
        CHECK(dot(c, par(c)) == scalar2dc(0.0));
        // par of par restores the object up to a positive factor
        CHECK(is_congruent(par(par(d)), d));
    }

    TEST_CASE("cga2dc: conformal distance and same-type dot products")
    {
        fmt::println("cga2dc: conformal distance and same-type dot products");

        // unitized round points: a1.a2 = -1/2 (v^2 + r1^2 + r2^2)
        auto a1 = round_point2dc(1.0, 2.0, 1.0);
        auto a2 = round_point2dc(4.0, 6.0, 0.0);
        CHECK(dot(a1, a2) == scalar2dc(-0.5 * (25.0 + 1.0 + 0.0)));

        // unitized circles (the codim-1 rounds): c1.c2 = +1/2 (v^2 - r1^2 - r2^2)
        auto c1 = circle2dc(1.0, 0.0, 2.0);
        auto c2 = circle2dc(3.0, 0.0, 2.0);
        CHECK(dot(c1, c2) == scalar2dc(0.5 * (4.0 - 4.0 - 4.0)));

        // orthogonal circles: r1^2 + r2^2 = v^2 (right-aligned intersection)
        auto o1 = circle2dc(0.0, 0.0, 3.0);
        auto o2 = circle2dc(5.0, 0.0, 4.0);
        CHECK(dot(o1, o2) == scalar2dc(0.0));
    }

    TEST_CASE("cga2dc: circle intersection end to end")
    {
        fmt::println("cga2dc: circle intersection end to end");

        // the showcase chain: meet -> radius sign test -> unitize ->
        // p_1,2 = cen(d) +/- (radius/weight) att(d), x/y coordinates
        auto const intersect = [](trivec2dc const& c1, trivec2dc const& c2, vec2dc& p1,
                                  vec2dc& p2) -> bool {
            auto d = rwdg(c1, c2);
            if (radius_sq(d) < 0.0) return false; // no intersection
            auto du = unitize(d);
            auto ctr = cen(du);
            auto dir = att(du);
            value_t const r = std::sqrt(radius_sq(du));
            // center normalized to weight one; only x/y carry meaning here
            value_t const izw = 1.0 / ctr.z;
            p1 = vec2dc(ctr.x * izw + r * dir.x, ctr.y * izw + r * dir.y, 1.0, 0.0);
            p2 = vec2dc(ctr.x * izw - r * dir.x, ctr.y * izw - r * dir.y, 1.0, 0.0);
            return true;
        };

        vec2dc p1, p2;

        // two r=2 circles at (1,0), (3,0): intersection (2, +/-rsqrt(3))
        CHECK(intersect(circle2dc(1.0, 0.0, 2.0), circle2dc(3.0, 0.0, 2.0), p1, p2));
        value_t const s3 = std::sqrt(3.0);
        bool const order_a = is_close(p1, vec2dc(2.0, s3, 1.0, 0.0)) &&
                             is_close(p2, vec2dc(2.0, -s3, 1.0, 0.0));
        bool const order_b = is_close(p1, vec2dc(2.0, -s3, 1.0, 0.0)) &&
                             is_close(p2, vec2dc(2.0, s3, 1.0, 0.0));
        CHECK((order_a || order_b));

        // externally tangent circles: single intersection point
        CHECK(intersect(circle2dc(0.0, 0.0, 1.0), circle2dc(3.0, 0.0, 2.0), p1, p2));
        CHECK(is_close(p1, vec2dc(1.0, 0.0, 1.0, 0.0)));
        CHECK(is_close(p2, vec2dc(1.0, 0.0, 1.0, 0.0)));

        // disjoint circles: no intersection (imaginary meet radius)
        CHECK_FALSE(
            intersect(circle2dc(0.0, 0.0, 1.0), circle2dc(5.0, 0.0, 1.0), p1, p2));

        // concentric circles: the meet has zero round weight (a flat object),
        // radius_sq must throw rather than give a bogus answer
        auto dcc = rwdg(circle2dc(0.0, 0.0, 1.0), circle2dc(0.0, 0.0, 2.0));
        CHECK_THROWS(radius_sq(dcc));

        // the same chain intersects a circle with a LINE, and the line need not
        // pass through the center: a line is the join of two of its points with
        // the point at infinity (a circle through infinity, see the flat objects
        // above). Circle r=2 about the origin against y = 1 (half the radius),
        // y = 2 (tangent) and y = 3 (missing it):
        auto const line_through = [](value_t px, value_t py, value_t dx,
                                     value_t dy) -> trivec2dc {
            return wdg(
                wdg(round_point2dc(px, py, 0.0), round_point2dc(px + dx, py + dy, 0.0)),
                e4_2dc);
        };

        auto const c_r2 = circle2dc(0.0, 0.0, 2.0);
        CHECK(intersect(c_r2, line_through(0.0, 1.0, 1.0, 0.0), p1, p2));
        bool const chord_a = is_close(p1, vec2dc(s3, 1.0, 1.0, 0.0)) &&
                             is_close(p2, vec2dc(-s3, 1.0, 1.0, 0.0));
        bool const chord_b = is_close(p1, vec2dc(-s3, 1.0, 1.0, 0.0)) &&
                             is_close(p2, vec2dc(s3, 1.0, 1.0, 0.0));
        CHECK((chord_a || chord_b)); // x = +/- sqrt(r^2 - 1^2) = +/- sqrt(3)

        CHECK(intersect(c_r2, line_through(0.0, 2.0, 1.0, 0.0), p1, p2)); // tangent
        CHECK(is_close(p1, vec2dc(0.0, 2.0, 1.0, 0.0)));
        CHECK(is_close(p2, vec2dc(0.0, 2.0, 1.0, 0.0)));

        CHECK_FALSE(intersect(c_r2, line_through(0.0, 3.0, 1.0, 0.0), p1, p2));
        CHECK(radius_sq(rwdg(c_r2, line_through(0.0, 3.0, 1.0, 0.0))) == -5.0);

        // the same chain as the library's dipole_points(), which returns NULL points
        // (round_point2dc(x, y, 0), not the weight-only vectors above): the secant
        // pair, the tangent pair coinciding, the imaginary meet throwing
        {
            auto const [q1, q2] =
                dipole_points(rwdg(circle2dc(1.0, 0.0, 2.0), circle2dc(3.0, 0.0, 2.0)));
            bool const oa = is_close(q1, round_point2dc(2.0, s3, 0.0)) &&
                            is_close(q2, round_point2dc(2.0, -s3, 0.0));
            bool const ob = is_close(q1, round_point2dc(2.0, -s3, 0.0)) &&
                            is_close(q2, round_point2dc(2.0, s3, 0.0));
            CHECK((oa || ob));
            auto const [t1, t2] =
                dipole_points(rwdg(circle2dc(0.0, 0.0, 1.0), circle2dc(3.0, 0.0, 2.0)));
            CHECK(is_close(t1, round_point2dc(1.0, 0.0, 0.0)));
            CHECK(is_close(t2, round_point2dc(1.0, 0.0, 0.0)));
            CHECK_THROWS(
                dipole_points(rwdg(circle2dc(0.0, 0.0, 1.0), circle2dc(5.0, 0.0, 1.0))));
        }

        // two-link inverse kinematics as that meet: base (0, 0), tip (2, 0), links
        // of length 2 -> the joint at (1, +/- sqrt 3), each solution at L1 from the
        // base and L2 from the tip; a tip beyond L1 + L2 is out of reach (throws),
        // and so is one closer than |L1 - L2|
        {
            auto const [j1, j2] = two_link_ik2dc(0.0, 0.0, 2.0, 0.0, 2.0, 2.0);
            bool const ka = is_close(j1, round_point2dc(1.0, s3, 0.0)) &&
                            is_close(j2, round_point2dc(1.0, -s3, 0.0));
            bool const kb = is_close(j1, round_point2dc(1.0, -s3, 0.0)) &&
                            is_close(j2, round_point2dc(1.0, s3, 0.0));
            CHECK((ka || kb));
            // an unequal pair, checked by the link lengths alone
            auto const [k1, k2] = two_link_ik2dc(0.3, -0.2, 1.1, 0.9, 0.8, 0.7);
            for (auto const& k : {k1, k2}) {
                CHECK(std::hypot(k.x - 0.3, k.y + 0.2) == doctest::Approx(0.8));
                CHECK(std::hypot(k.x - 1.1, k.y - 0.9) == doctest::Approx(0.7));
            }
            CHECK_THROWS(two_link_ik2dc(0.0, 0.0, 5.0, 0.0, 2.0, 2.0));
            CHECK_THROWS(two_link_ik2dc(0.0, 0.0, 0.5, 0.0, 2.0, 1.0));
        }
    }

    TEST_CASE("cga2dc: conformal conjugate and containment")
    {
        fmt::println("cga2dc: conformal conjugate and containment");

        auto a = round_point2dc(1.0, 2.0, 1.0);
        auto d = dipole2dc(1.0, 2.0, 2.0, 0.0, 1.0);
        auto c = circle2dc(1.0, 2.0, 1.0);

        // involution
        CHECK(cconj(cconj(a)) == a);
        CHECK(cconj(cconj(d)) == d);
        CHECK(cconj(cconj(c)) == c);

        // alternate center norm: dot(u, cconj(u)) = center_nrm_sq(u)
        // = squared distance origin <-> center for unitized round objects
        CHECK(dot(a, cconj(a)) == scalar2dc(5.0));
        CHECK(dot(c, cconj(c)) == scalar2dc(5.0));
        CHECK(scalar2dc(center_nrm_sq(a)) == dot(a, cconj(a)));
        CHECK(scalar2dc(center_nrm_sq(c)) == dot(c, cconj(c)));
        CHECK(is_close(scalar2dc(center_nrm_sq(d)),
                       scalar2dc(value_t(dot(d, cconj(d))) / round_weight_nrm_sq(d))));

        // the conjugate's center is reflected through the origin, its radius
        // adjusted: for a unitized round point r'^2 = -2*|p|^2 - r^2, for a
        // unitized circle r'^2 = 2*|p|^2 - r^2
        CHECK(is_close(scalar2dc(radius_sq(cconj(a))), scalar2dc(-11.0)));
        CHECK(is_close(scalar2dc(radius_sq(cconj(c))), scalar2dc(9.0)));
        auto cc = cen(cconj(c));
        CHECK(is_close(vec2dc(cc.x / cc.z, cc.y / cc.z, 1.0, 0.0),
                       vec2dc(-1.0, -2.0, 1.0, 0.0)));

        // containment: wdg(a, c) = -1/2 (v^2 + r^2 - R^2) * I for a unitized
        // round point (radius r) and circle (radius R) with center distance v
        auto const contain = [](value_t px, value_t py, value_t r, value_t cx, value_t cy,
                                value_t R) {
            auto val = wdg(round_point2dc(px, py, r), circle2dc(cx, cy, R));
            value_t const v_sq = (px - cx) * (px - cx) + (py - cy) * (py - cy);
            CHECK(val == pscalar2dc(-0.5 * (v_sq + r * r - R * R)));
        };
        contain(0.0, 0.0, 0.0, 1.0, 0.0, 2.0); // off-circle null point
        contain(3.0, 0.0, 0.0, 1.0, 0.0, 2.0); // on-circle null point -> 0
        contain(2.0, 1.0, 1.0, 1.0, 1.0, 2.0); // polar aligned: v^2+r^2 = R^2 -> 0
        contain(2.0, -1.0, 1.5, -1.0, 0.5, 2.0);

        // every round object contains its own round center: u ^ cen(u) = 0
        CHECK(wdg(a, cen(a)) == bivec2dc());
        CHECK(is_close(wdg(d, cen(d)), trivec2dc()));
        CHECK(wdg(c, cen(c)) == pscalar2dc(0.0));

        // dipole surface points p_pm = cen(d) +/- radius_nrm * att(d) (the
        // formula behind the circle-intersection chain), unitized dipole
        auto du = unitize(d);
        auto ctr = cen(du);
        auto dir = att(du);
        value_t const rn = std::sqrt(radius_sq(du));
        auto p1 =
            vec2dc(ctr.x / ctr.z + rn * dir.x, ctr.y / ctr.z + rn * dir.y, 1.0, 0.0);
        auto p2 =
            vec2dc(ctr.x / ctr.z - rn * dir.x, ctr.y / ctr.z - rn * dir.y, 1.0, 0.0);
        // d was built with center (1,2), r=2, normal (0,1): points (-1,2), (3,2)
        bool const ok_a = is_close(p1, vec2dc(3.0, 2.0, 1.0, 0.0)) &&
                          is_close(p2, vec2dc(-1.0, 2.0, 1.0, 0.0));
        bool const ok_b = is_close(p1, vec2dc(-1.0, 2.0, 1.0, 0.0)) &&
                          is_close(p2, vec2dc(3.0, 2.0, 1.0, 0.0));
        CHECK((ok_a || ok_b));
    }

    TEST_CASE("cga2dc: exp, log and sqrt of regressive versors")
    {
        fmt::println("cga2dc: exp, log and sqrt of regressive versors");

        auto const I_e = mvec2dc_e(pscalar2dc(1.0));

        // line through two Euclidean points (a circle through infinity)
        auto const line_through = [](value_t x1, value_t y1, value_t x2, value_t y2) {
            return wdg(wdg(point2dc(x1, y1), point2dc(x2, y2)), e4_2dc);
        };
        // even-versor sandwich, spelled out with rgpr until transform() lands
        auto const sandwich = [](mvec2dc_e const& M, vec2dc const& X) {
            return gr1(rgpr(rgpr(M, X), rrev(M)));
        };

        // the identity: exp of the zero generator is the rgpr identity I
        CHECK(rexp(bivec2dc()) == I_e);

        // group gates on generic generators
        for (int i = 0; i < 3; ++i) {
            auto B = bivec2dc(0.1 * rnd_int(), 0.1 * rnd_int(), 0.1 * rnd_int(),
                              0.1 * rnd_int(), 0.1 * rnd_int(), 0.1 * rnd_int());
            auto U = rexp(B);
            CHECK(is_close(rgpr(U, rrev(U)), I_e)); // unit versor
            CHECK(is_close(rgpr(rexp(0.3 * B), rexp(0.5 * B)), rexp(0.8 * B)));
            CHECK(is_close(rlog(U), B)); // principal
            CHECK(is_close(rgpr(rsqrt(U), rsqrt(U)), U));
            CHECK(is_close(rsqrt(U), rexp(0.5 * B)));
        }

        // TRANSLATION: reflections in the parallel lines x=0 and x=1 compose
        // to the translation by (2,0); the versor is exactly I - e24 and its
        // generator is exactly parabolic (rgpr(B,B) = 0, series terminates)
        auto Mt = rgpr(line_through(1.0, 0.0, 1.0, 1.0),  // x = 1
                       line_through(0.0, 0.0, 0.0, 1.0)); // x = 0
        CHECK(Mt == mvec2dc_e(scalar2dc(0.0), bivec2dc(0.0, 0.0, 0.0, 0.0, -1.0, 0.0),
                              pscalar2dc(1.0)));
        CHECK(rgpr(Mt, rrev(Mt)) == I_e);
        CHECK(sandwich(Mt, point2dc(0.5, 2.0)) == point2dc(2.5, 2.0)); // exact
        auto Bt = rlog(Mt);
        CHECK(rgpr(Bt, Bt) == mvec2dc_e(scalar2dc(0.0))); // parabolic generator
        CHECK(rexp(Bt) == Mt);

        // ROTATION about the origin by pi/2 (line pair at 45 degrees); the
        // raw two-line versor is not unit -- normalize the sandwiched point
        auto Mr =
            rgpr(line_through(0.0, 0.0, 1.0, 1.0), line_through(0.0, 0.0, 1.0, 0.0));
        CHECK(is_close(unitize(sandwich(Mr, point2dc(1.0, 0.0))), point2dc(0.0, 1.0)));

        // a full turn (2 pi) exponentiates to -I, where the principal
        // logarithm is singular by construction
        auto Br = rlog(mvec2dc_e(1.0 / std::sqrt(2.0) * Mr)); // pi/2 generator
        CHECK(is_close(rexp(4.0 * Br), mvec2dc_e(pscalar2dc(-1.0))));
        CHECK_THROWS(rlog(rexp(4.0 * Br)));

        // DILATION: inversion in concentric circles r=1 then r=2 scales by 4
        auto Md = rgpr(circle2dc(0.0, 0.0, 2.0), circle2dc(0.0, 0.0, 1.0));
        CHECK(is_close(unitize(sandwich(Md, point2dc(1.0, 0.0))), point2dc(4.0, 0.0)));
        CHECK(is_close(unitize(sandwich(Md, point2dc(-0.5, 0.25))), point2dc(-2.0, 1.0)));
    }

    TEST_CASE("cga2dc: conformal transformations")
    {
        fmt::println("cga2dc: conformal transformations");

        auto const up = [](vec2dc const& v) { return unitize(v); };

        // TRANSLATION: exact on embedded points; matches the two-line versor
        auto Tr = get_translation(2.0, 0.0);
        CHECK(transform(point2dc(0.5, 2.0), Tr) == point2dc(2.5, 2.0));
        auto line_x0 = wdg(wdg(point2dc(0.0, 0.0), point2dc(0.0, 1.0)), e4_2dc);
        auto line_x1 = wdg(wdg(point2dc(1.0, 0.0), point2dc(1.0, 1.0)), e4_2dc);
        CHECK(is_same_transform(Tr, mvec2dc_e(rgpr(line_x1, line_x0))));
        // a translated circle keeps its radius
        auto c = circle2dc(1.0, 2.0, 1.5);
        CHECK(is_close(scalar2dc(radius_sq(transform(c, Tr))), scalar2dc(radius_sq(c))));

        // ROTATION about a point (counterclockwise for positive angles)
        auto Ro = get_rotation(0.0, 0.0, 0.5 * pi);
        CHECK(is_close(up(transform(point2dc(1.0, 0.0), Ro)), point2dc(0.0, 1.0)));
        auto Rc = get_rotation(1.0, 1.0, 0.5 * pi);
        CHECK(is_close(up(transform(point2dc(2.0, 1.0), Rc)), point2dc(1.0, 2.0)));
        CHECK(is_close(scalar2dc(radius_sq(transform(c, Rc))), scalar2dc(radius_sq(c))));
        // rotation angles add; a full turn acts as the identity
        CHECK(is_same_transform(rgpr(Ro, Ro), get_rotation(0.0, 0.0, pi)));
        CHECK(is_same_transform(get_rotation(0.5, 0.5, 2.0 * pi),
                                mvec2dc_e(pscalar2dc(1.0))));

        // DILATION about a point; radius scales by sigma
        auto Di = get_dilation(0.0, 0.0, 4.0);
        CHECK(is_close(up(transform(point2dc(1.0, 0.0), Di)), point2dc(4.0, 0.0)));
        auto Dm = get_dilation(1.0, 1.0, 2.0);
        CHECK(is_close(up(transform(point2dc(2.0, 1.0), Dm)), point2dc(3.0, 1.0)));
        CHECK(is_close(scalar2dc(radius_sq(transform(c, Dm))),
                       scalar2dc(4.0 * radius_sq(c))));
        // matches the concentric-circle inversion pair
        auto Dc = rgpr(circle2dc(0.0, 0.0, 2.0), circle2dc(0.0, 0.0, 1.0));
        CHECK(is_same_transform(Di, mvec2dc_e(Dc)));
        CHECK_THROWS(get_dilation(0.0, 0.0, 0.0));

        // INVERSION in the unit circle: geometric reciprocal r^2/x, on-circle
        // points fixed, involutive
        auto uc = circle2dc(0.0, 0.0, 1.0);
        CHECK(is_congruent(invert_on(point2dc(2.0, 0.0), uc), point2dc(0.5, 0.0)));
        CHECK(is_congruent(invert_on(point2dc(0.0, -1.0), uc), point2dc(0.0, -1.0)));
        CHECK(is_congruent(invert_on(invert_on(point2dc(0.3, -0.7), uc), uc),
                           point2dc(0.3, -0.7)));
        // a circle through the center of inversion maps to a LINE (round
        // weight zero); its image contains the image of any circle point
        auto c0 = circle2dc(1.0, 0.0, 1.0); // passes through the origin
        auto img = invert_on(c0, uc);
        CHECK(std::abs(round_weight_nrm_sq(img)) < 1.0e-12);
        CHECK(std::abs(value_t(wdg(invert_on(point2dc(2.0, 0.0), uc), img))) < 1.0e-12);
        // a circle NOT through the center maps to a circle
        CHECK(round_weight_nrm_sq(invert_on(circle2dc(3.0, 0.0, 1.0), uc)) > 0.1);

        // versor equality is decided by action: M vs -M and uniform rescaling
        // act identically; genuinely different motors do not
        CHECK(is_same_transform(Ro, mvec2dc_e(-1.0 * Ro)));
        CHECK(is_same_transform(Ro, mvec2dc_e(2.0 * Ro)));
        CHECK_FALSE(is_same_transform(Ro, get_rotation(0.0, 0.0, pi)));
        CHECK_FALSE(is_same_transform(Tr, Di));

        // composition sanity: rotation about m = T(m) o rotation about origin
        // o T(-m)
        auto Rc2 = rgpr(rgpr(get_translation(1.0, 1.0), Ro), get_translation(-1.0, -1.0));
        CHECK(is_same_transform(Rc, mvec2dc_e(Rc2)));
    }

    TEST_CASE("cga2dc: exported extended metric arrays")
    {
        fmt::println("cga2dc: exported extended metric arrays");

        // the regressive extended metric is the anti-exomorphism = -G
        for (size_t i = 0; i < cga2dc_metric.size(); ++i) {
            CHECK(cga2dc_rmetric[i] == -cga2dc_metric[i]);
        }

        // G against the LIVE dot product, per grade block (dot(e_i, e_j) is
        // the extended metric pairing); rdot likewise against the regressive
        // metric
        auto G = cga2dc_metric_view();
        auto Gr = cga2dc_rmetric_view();
        CHECK(value_t(dot(scalar2dc(1.0), scalar2dc(1.0))) == G[0, 0]);
        CHECK(value_t(rdot(scalar2dc(1.0), scalar2dc(1.0))) == Gr[0, 0]);
        std::array<vec2dc, 4> const vs{e1_2dc, e2_2dc, e3_2dc, e4_2dc};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                CHECK(value_t(dot(vs[i], vs[j])) == G[i + 1, j + 1]);
                CHECK(value_t(rdot(vs[i], vs[j])) == Gr[i + 1, j + 1]);
            }
        }
        std::array<bivec2dc, 6> const bs{e31_2dc, e32_2dc, e12_2dc,
                                         e14_2dc, e24_2dc, e34_2dc};
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                CHECK(value_t(dot(bs[i], bs[j])) == G[i + 5, j + 5]);
                CHECK(value_t(rdot(bs[i], bs[j])) == Gr[i + 5, j + 5]);
            }
        }
        std::array<trivec2dc, 4> const ts{e314_2dc, e324_2dc, e124_2dc, e321_2dc};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                CHECK(value_t(dot(ts[i], ts[j])) == G[i + 11, j + 11]);
                CHECK(value_t(rdot(ts[i], ts[j])) == Gr[i + 11, j + 11]);
            }
        }
        CHECK(value_t(dot(I_2dc, I_2dc)) == G[15, 15]);
        CHECK(value_t(rdot(I_2dc, I_2dc)) == Gr[15, 15]);
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


    TEST_CASE("DualNum2dc: is_close")
    {
        fmt::println("DualNum2dc: is_close");

        auto D = dualnum2dc{2.0, -3.0};
        auto E = D;
        CHECK(is_close(D, E));
        E.c0 += 1.0e-6;
        CHECK(!is_close(D, E));
    }

} // TEST_SUITE("CGA 2dc Tests")
