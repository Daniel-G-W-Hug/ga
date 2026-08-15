// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include "fmt/format.h" // formatting

// include functions to be tested
#include "ga/ga_cga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::cga; // use specific operations of CGA (Conformal Algebra)


/////////////////////////////////////////////////////////////////////////////////////////
// Conformal Geometric Algebra of 3d Euclidean space: cga3dc — type layer
// (G(4,1,0) after diagonalization of the native null-basis metric; the
// operations layer follows with the generator round)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("CGA 3dc Tests")
{

    TEST_CASE("G<4,1,0>: cga3dc")
    {
        fmt::println("G<4,1,0>: cga3dc");
        // Conformal algebra of 3d Euclidean space; the native null-basis
        // metric (e4.e4 = e5.e5 = 0, e4.e5 = -1) diagonalizes to G(4,1,0)
        const algebra<4, 1, 0> alg;
        CHECK(alg.p() == 4);
        CHECK(alg.n() == 1);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 5);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 32);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 6); // == dim_space + 1
        CHECK(alg.num_components_grade[2] == 10);
        CHECK(alg.num_components_grade[3] == 10);
        fmt::println("   cga3dc: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   cga3dc: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    TEST_CASE("cga3dc: type layer round trips")
    {
        fmt::println("cga3dc: type layer round trips");

        auto s = scalar3dc(3.0);
        auto v = vec3dc(1.0, 2.0, 3.0, 4.0, 5.0);
        auto B = bivec3dc(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
        auto t = trivec3dc(10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);
        auto Q = quadvec3dc(5.0, 4.0, 3.0, 2.0, 1.0);
        auto ps = pscalar3dc(-2.0);

        // grades and antigrades of the graded types
        CHECK(gr(s) == 0);
        CHECK(gr(v) == 1);
        CHECK(gr(B) == 2);
        CHECK(gr(t) == 3);
        CHECK(gr(Q) == 4);
        CHECK(gr(ps) == 5);
        CHECK(rgr(s) == 5);
        CHECK(rgr(v) == 4);
        CHECK(rgr(B) == 3);
        CHECK(rgr(t) == 2);
        CHECK(rgr(Q) == 1);
        CHECK(rgr(ps) == 0);

        // graded-part constructors and grade selection round-trip
        mvec3dc M(s, v, B, t, Q, ps);
        CHECK(gr0(M) == s);
        CHECK(gr1(M) == v);
        CHECK(gr2(M) == B);
        CHECK(gr3(M) == t);
        CHECK(gr4(M) == Q);
        CHECK(gr5(M) == ps);

        mvec3dc_e Me(s, B, Q);
        CHECK(gr0(Me) == s);
        CHECK(gr2(Me) == B);
        CHECK(gr4(Me) == Q);

        mvec3dc_u Mu(v, t, ps);
        CHECK(gr1(Mu) == v);
        CHECK(gr3(Mu) == t);
        CHECK(gr5(Mu) == ps);

        // even/odd embedding into the full multivector recomposes it
        CHECK(mvec3dc(Me) + mvec3dc(Mu) == M);

        // graded-part combination operators land in the right subalgebra
        CHECK(gr2(s + B) == B);
        CHECK(gr4(B + Q) == Q);
        CHECK(gr0(s + Q) == s);
        CHECK(gr3(v + t) == t);
        CHECK(gr5(v + ps) == ps);
        CHECK(gr1(t + v) == v);
        // mixed parity lands in the full multivector
        CHECK(gr1(s + v) == v);
        CHECK(gr5(B - ps) == -ps);
        CHECK(gr4(t + Q) == Q);

        // basis constants: null pair and pseudoscalar
        CHECK(e4_3dc == vec3dc(0.0, 0.0, 0.0, 1.0, 0.0));
        CHECK(e5_3dc == vec3dc(0.0, 0.0, 0.0, 0.0, 1.0));
        CHECK(I_3dc == pscalar3dc(1.0));
        CHECK(gr2(e45_3dc_mv) == e45_3dc);
        CHECK(gr3(e321_3dc_mv_u) == e321_3dc);
        CHECK(gr4(e1234_3dc_mv_e) == e1234_3dc);

        // scalar multiplication and division on the new templates
        CHECK(2.0 * v / 2.0 == v);
        CHECK(0.5 * (B + B) == B);
        CHECK(-(-M) == M);
    }

    TEST_CASE("cga3dc: exported extended metric arrays")
    {
        fmt::println("cga3dc: exported extended metric arrays");

        // the regressive extended metric is the anti-exomorphism = -G
        for (size_t i = 0; i < cga3dc_metric.size(); ++i) {
            CHECK(cga3dc_rmetric[i] == -cga3dc_metric[i]);
        }

        // pinned rows of the reviewed metric (the operation-level cross-check
        // against dot/rdot follows with the cga3dc ops layer): the null pair,
        // the w<->u bivector pairing, the quadvector partner pair, and the
        // pseudoscalar
        auto G = cga3dc_metric_view();
        CHECK(G[0, 0] == 1);
        CHECK(G[4, 5] == -1); // e4 . e5
        CHECK(G[5, 4] == -1);
        CHECK(G[4, 4] == 0); // null vectors
        CHECK(G[5, 5] == 0);
        CHECK(G[6, 12] == 1);   // e41 <-> e15
        CHECK(G[15, 15] == -1); // e45
        CHECK(G[26, 26] == -1); // e4235 (contains the null pair e45)
        CHECK(G[29, 30] == 1);  // e3215 <-> e1234
        CHECK(G[30, 29] == 1);
        CHECK(G[31, 31] == -1); // e12345
        // every row is a signed permutation: exactly one non-zero entry
        for (size_t i = 0; i < 32; ++i) {
            int nz = 0;
            for (size_t j = 0; j < 32; ++j) {
                if (G[i, j] != 0) ++nz;
            }
            CHECK(nz == 1);
        }
    }

    TEST_CASE("cga3dc: involutions, complement, dual, antidual")
    {
        fmt::println("cga3dc: involutions, complement, dual, antidual");

        // unit basis blade k as a full multivector
        auto blade = [](size_t k) {
            auto c = [k](size_t i) { return i == k ? 1.0 : 0.0; };
            return mvec3dc(c(0), c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9),
                           c(10), c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18),
                           c(19), c(20), c(21), c(22), c(23), c(24), c(25), c(26), c(27),
                           c(28), c(29), c(30), c(31));
        };

        // grade sign patterns (n=5): gr_inv - + at even grades; rev/rrev
        // + + - - + +; conj + - - + + -
        auto v = vec3dc(1.0, 2.0, 3.0, 4.0, 5.0);
        auto B = bivec3dc(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
        auto t = trivec3dc(10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);
        auto Q = quadvec3dc(5.0, 4.0, 3.0, 2.0, 1.0);
        CHECK(gr_inv(v) == -v);
        CHECK(gr_inv(B) == B);
        CHECK(gr_inv(t) == -t);
        CHECK(gr_inv(Q) == Q);
        CHECK(gr_inv(pscalar3dc(1.0)) == pscalar3dc(-1.0));
        CHECK(rev(v) == v);
        CHECK(rev(B) == -B);
        CHECK(rev(t) == -t);
        CHECK(rev(Q) == Q);
        CHECK(rev(pscalar3dc(1.0)) == pscalar3dc(1.0));
        CHECK(conj(v) == -v);
        CHECK(conj(B) == -B);
        CHECK(conj(t) == t);
        CHECK(conj(Q) == Q);
        CHECK(conj(pscalar3dc(1.0)) == pscalar3dc(-1.0));

        mvec3dc const M(scalar3dc(3.0), v, B, t, Q, pscalar3dc(-2.0));
        // rev == rrev in 5d; conj = gr_inv(rev); all involutions square to id
        CHECK(rev(M) == rrev(M));
        CHECK(conj(M) == gr_inv(rev(M)));
        CHECK(gr_inv(gr_inv(M)) == M);
        CHECK(rev(rev(M)) == M);
        CHECK(conj(conj(M)) == M);

        // complement: u ^ cmpl(u) = I_3dc and cmpl(cmpl(u)) = u per unit blade
        for (size_t i = 0; i < 32; ++i) {
            auto const u = blade(i);
            CHECK(wdg(u, cmpl(u)) == blade(31));
            CHECK(cmpl(cmpl(u)) == u);
        }

        // transcription gate: dual(e_i) = cmpl(G e_i) = sum_j G[i,j] cmpl(e_j),
        // antidual(e_i) = sum_j AG[i,j] cmpl(e_j) with AG = -G -- exact,
        // sign-sensitive pins against the exported extended metric arrays
        auto G = cga3dc_metric_view();
        auto Gr = cga3dc_rmetric_view();
        for (size_t i = 0; i < 32; ++i) {
            mvec3dc d_expect, ad_expect;
            for (size_t j = 0; j < 32; ++j) {
                if (G[i, j] != 0) d_expect += value_t(G[i, j]) * cmpl(blade(j));
                if (Gr[i, j] != 0) ad_expect += value_t(Gr[i, j]) * cmpl(blade(j));
            }
            CHECK(dual(blade(i)) == d_expect);
            CHECK(antidual(blade(i)) == ad_expect);
            CHECK(antidual(blade(i)) == -dual(blade(i)));
        }
    }

    TEST_CASE("cga3dc: products - identities and anchors")
    {
        fmt::println("cga3dc: products - identities and anchors");

        auto blade = [](size_t k) {
            auto c = [k](size_t i) { return i == k ? 1.0 : 0.0; };
            return mvec3dc(c(0), c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9),
                           c(10), c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18),
                           c(19), c(20), c(21), c(22), c(23), c(24), c(25), c(26), c(27),
                           c(28), c(29), c(30), c(31));
        };

        // two full multivectors with small integer coefficients (exact arithmetic)
        mvec3dc M1, M2;
        for (size_t i = 0; i < 32; ++i) {
            M1 += value_t(i % 5 + 1) * blade(i);
            M2 += value_t(double(i * 7 % 9) - 4.0) * blade(i);
        }

        // product identity elements
        CHECK(M1 * mvec3dc(scalar3dc(1.0)) == M1);
        CHECK(mvec3dc(scalar3dc(1.0)) * M1 == M1);
        CHECK(rgpr(M1, blade(31)) == M1);
        CHECK(rgpr(blade(31), M1) == M1);

        // contractions are the rwdg-duals of the interior products
        CHECK((M1 << M2) == rwdg(dual(M1), M2));
        CHECK((M1 >> M2) == rwdg(M1, dual(M2)));

        // rwdg is the complement-mapped wedge (odd-dim: singular cmpl)
        CHECK(rwdg(M1, M2) == cmpl(wdg(cmpl(M1), cmpl(M2))));

        // rgpr is the complement-mapped geometric product
        CHECK(rgpr(M1, M2) == cmpl(cmpl(M1) * cmpl(M2)));

        // commutators are the asymmetric parts of gpr/rgpr
        CHECK(cmt(M1, M2) == 0.5 * (M1 * M2 - M2 * M1));
        CHECK(rcmt(M1, M2) == 0.5 * (rgpr(M1, M2) - rgpr(M2, M1)));

        // equal-grade contraction reduces to the metric inner product
        auto v1 = vec3dc(1.0, 2.0, 3.0, 4.0, 5.0);
        auto v2 = vec3dc(-2.0, 1.0, 0.0, 3.0, -1.0);
        CHECK(value_t(v1 << v2) == value_t(dot(v1, v2)));
        CHECK(value_t(v1 >> v2) == value_t(dot(v1, v2)));

        // null-pair anchors (book-reviewed gpr table): e4 * e5 = -1 + e45,
        // e45 * e45 = 1; the null basis vectors square to zero
        CHECK(e4_3dc * e5_3dc == mvec3dc_e(scalar3dc(-1.0), e45_3dc));
        CHECK(e45_3dc * e45_3dc == mvec3dc_e(scalar3dc(1.0)));
        CHECK(e4_3dc * e4_3dc == mvec3dc_e(scalar3dc(0.0)));
        CHECK(e5_3dc * e5_3dc == mvec3dc_e(scalar3dc(0.0)));

        // center of the rgpr algebra: rgpr(1,1) = -I (the scalar acts as the
        // imaginary unit w.r.t. the regressive product, as in cga2dc)
        CHECK(rgpr(scalar3dc(1.0), scalar3dc(1.0)) == pscalar3dc(-1.0));

        // expansions are the wdg-based duals of the contractions
        CHECK(l_expand3dc(v1, e12_3dc) == wdg(dual(v1), e12_3dc));
        CHECK(r_expand3dc(e12_3dc, v1) == wdg(e12_3dc, dual(v1)));
    }

    TEST_CASE("cga3dc: inv and rinv")
    {
        fmt::println("cga3dc: inv and rinv");

        // pinned values (non-degenerate metric: BOTH pseudoscalar inv and
        // scalar rinv exist, unlike in the degenerate pga algebras)
        CHECK(inv(scalar3dc(2.0)) == scalar3dc(0.5));
        CHECK(inv(I_3dc) == pscalar3dc(-1.0));          // I * I = -1
        CHECK(rinv(scalar3dc(1.0)) == scalar3dc(-1.0)); // rgpr(1, -1) = I
        CHECK(is_close(rgpr(scalar3dc(1.0), rinv(scalar3dc(1.0))), pscalar3dc(1.0)));

        // null elements are not invertible -- including the null basis vectors
        CHECK_THROWS(inv(e4_3dc));
        CHECK_THROWS(inv(e5_3dc));

        // gate: u * inv(u) = 1 for generically non-null elements, both
        // directions (H&S n=5: a left and a right inverse coincide)
        auto const one_e = mvec3dc_e(scalar3dc(1.0));
        auto const one_m = mvec3dc(scalar3dc(1.0));
        auto const I_e = mvec3dc(pscalar3dc(1.0)); // rgpr(graded,...) returns mv types
        auto v = vec3dc(1.0, 2.0, -1.0, 0.5, 2.0);
        auto B = bivec3dc(1.0, -2.0, 0.5, 3.0, 0.0, 1.0, -1.5, 2.0, 0.0, 1.0);
        auto t = trivec3dc(0.5, 1.0, -2.0, 0.0, 3.0, -1.0, 2.0, 0.5, 1.0, -0.5);
        auto Q = quadvec3dc(2.0, -1.0, 0.5, 3.0, 1.0);
        auto Me = mvec3dc_e(scalar3dc(2.0), B, Q);
        auto Mu = mvec3dc_u(v, t, pscalar3dc(1.5));
        mvec3dc const M(scalar3dc(1.0), v, B, t, Q, pscalar3dc(-2.0));

        CHECK(is_close(v * inv(v), one_e));
        CHECK(is_close(inv(v) * v, one_e));
        CHECK(is_close(mvec3dc(B * inv(B)), mvec3dc(one_e)));
        CHECK(is_close(mvec3dc(inv(B) * B), mvec3dc(one_e)));
        CHECK(is_close(mvec3dc(t * inv(t)), mvec3dc(one_e)));
        CHECK(is_close(mvec3dc(inv(t) * t), mvec3dc(one_e)));
        CHECK(is_close(mvec3dc(Q * inv(Q)), mvec3dc(one_e)));
        CHECK(is_close(mvec3dc(inv(Q) * Q), mvec3dc(one_e)));
        CHECK(is_close(Me * inv(Me), one_e));
        CHECK(is_close(inv(Me) * Me, one_e));
        CHECK(is_close(mvec3dc(Mu * inv(Mu)), mvec3dc(one_e)));
        CHECK(is_close(mvec3dc(inv(Mu) * Mu), mvec3dc(one_e)));
        CHECK(is_close(M * inv(M), one_m));
        CHECK(is_close(inv(M) * M, one_m));

        // gate: rgpr(u, rinv(u)) = I (the rgpr identity element)
        CHECK(is_close(mvec3dc(rgpr(v, rinv(v))), I_e));
        CHECK(is_close(mvec3dc(rgpr(rinv(v), v)), I_e));
        CHECK(is_close(mvec3dc(rgpr(B, rinv(B))), I_e));
        CHECK(is_close(mvec3dc(rgpr(t, rinv(t))), I_e));
        CHECK(is_close(mvec3dc(rgpr(Q, rinv(Q))), I_e));
        CHECK(is_close(rgpr(M, rinv(M)), I_e));
        CHECK(is_close(rgpr(rinv(M), M), I_e));
    }

    TEST_CASE("cga3dc: congruence and closeness")
    {
        fmt::println("cga3dc: congruence and closeness");

        auto v = vec3dc(1.0, -2.0, 3.0, 0.5, 1.0);
        auto B = bivec3dc(1.0, 0.0, -2.0, 3.0, 0.0, 1.5, 0.5, -1.0, 0.0, 2.0);
        auto t = trivec3dc(-1.0, 2.0, 0.0, 4.0, 1.0, 0.0, -2.0, 0.5, 3.0, 0.0);
        auto Q = quadvec3dc(1.0, -2.0, 0.0, 3.0, 0.5);

        // congruent: same subspace up to any non-zero scale (sign included)
        CHECK(is_congruent(v, vec3dc(2.0 * v)));
        CHECK(is_congruent(v, vec3dc(-3.0 * v)));
        CHECK(is_congruent(B, bivec3dc(-0.5 * B)));
        CHECK(is_congruent(t, trivec3dc(7.0 * t)));
        CHECK(is_congruent(Q, quadvec3dc(-2.5 * Q)));
        CHECK_FALSE(is_congruent(v, e1_3dc));
        CHECK_FALSE(is_congruent(B, e12_3dc));
        CHECK(is_congruent(scalar3dc(2.0), scalar3dc(-5.0)));
        CHECK(is_congruent(I_3dc, pscalar3dc(-3.0)));

        // is_close: same value within a relative tolerance
        CHECK(is_close(v, v));
        CHECK_FALSE(is_close(v, vec3dc(2.0 * v)));
        CHECK(is_close(mvec3dc(v) + mvec3dc(t), mvec3dc(v) + mvec3dc(t)));
        CHECK(is_close(mvec3dc_e(scalar3dc(1.0), B, Q), mvec3dc_e(scalar3dc(1.0), B, Q)));
        CHECK_FALSE(is_close(mvec3dc(v), mvec3dc(t)));
    }

    TEST_CASE("cga3dc: fmt printing")
    {
        fmt::println("cga3dc: fmt printing");
        auto v = vec3dc(1.0, 2.0, 3.0, 4.0, 5.0);
        auto B = bivec3dc(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0);
        auto Q = quadvec3dc(1.0, 0.0, 0.0, 0.0, -1.0);
        CHECK(fmt::format("{}", v).size() > 0);
        fmt::println("   vec3dc:     {}", v);
        fmt::println("   bivec3dc:   {}", B);
        fmt::println("   quadvec3dc: {}", Q);
        fmt::println("   mvec3dc_e:  {}", mvec3dc_e(scalar3dc(1.0), B, Q));
        fmt::println("   mvec3dc:    {}", mvec3dc(v));
        fmt::println("");
    }

} // TEST_SUITE("CGA 3dc Tests")
