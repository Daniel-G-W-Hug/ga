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

    TEST_CASE("cga3dc: geometric objects, split and norms")
    {
        fmt::println("cga3dc: geometric objects, split and norms");

        // round point at (1,2,3), r = 2 (all exact pins probe-verified)
        auto a = round_point3dc(1.0, 2.0, 3.0, 2.0);
        CHECK(a == vec3dc(1.0, 2.0, 3.0, 1.0, 9.0)); // u = (p^2 + r^2)/2 = 9
        CHECK(radius_sq(a) == 4.0);
        CHECK(center_nrm_sq(a) == 14.0); // round bulk^2 (no flat weight at gr 1)
        // null-point embedding: q * q = 0
        auto q = round_point3dc(3.0, 2.0, 3.0, 0.0);
        CHECK(is_close(mvec3dc(q * q), mvec3dc(scalar3dc(0.0))));

        // sphere: matches the container of the round point exactly
        auto s = sphere3dc(1.0, 2.0, 3.0, 2.0);
        CHECK(s == quadvec3dc(1.0, 2.0, 3.0, -5.0, -1.0));
        CHECK(s == con(a));
        CHECK(radius_sq(s) == 4.0);
        // containment: wdg(point, sphere) = -1/2 (v^2 + r^2 - R^2) * I for
        // unitized operands (q lies ON s: v = 2, r = 0, R = 2)
        CHECK(wdg(q, s) == pscalar3dc(0.0));
        CHECK(wdg(round_point3dc(4.0, 2.0, 3.0, 0.0), s) == pscalar3dc(-2.5));
        CHECK(wdg(round_point3dc(4.0, 2.0, 3.0, 1.0), s) == pscalar3dc(-3.0));

        // circle at (1,2,3), r = 2, normal +z; dipole same center/radius, axis +x
        auto c = circle3dc(1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 1.0);
        CHECK(c == trivec3dc(2.0, -1.0, 0.0, 3.0, 6.0, 4.0, 0.0, 0.0, 1.0, -3.0));
        CHECK(radius_sq(c) == 4.0);
        auto d = dipole3dc(1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 0.0);
        CHECK(d == bivec3dc(1.0, 0.0, 0.0, 0.0, 3.0, -2.0, -8.0, 2.0, 3.0, 1.0));
        CHECK(radius_sq(d) == 4.0);

        // flats: zero round weight (radius_sq/unitize throw), exact pins
        auto fp = flat_point3dc(1.0, 2.0, 3.0);
        CHECK(fp == bivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 1.0));
        auto l = line3dc(1.0, 2.0, 3.0, 1.0, 0.0, 0.0);
        CHECK(l == trivec3dc(1.0, 0.0, 0.0, 0.0, 3.0, -2.0, 0.0, 0.0, 0.0, 0.0));
        auto pl = plane3dc(0.0, 0.0, 1.0, 3.0);
        CHECK(pl == quadvec3dc(0.0, 0.0, 1.0, -3.0, 0.0));
        CHECK_THROWS(radius_sq(pl));
        CHECK_THROWS(unitize(l));
        // incidence of flats: q on the line / on the plane wedges to zero
        CHECK(wdg(round_point3dc(5.0, 2.0, 3.0, 0.0), l) ==
              quadvec3dc(0.0, 0.0, 0.0, 0.0, 0.0));
        CHECK(wdg(round_point3dc(5.0, 7.0, 3.0, 0.0), pl) == pscalar3dc(0.0));
        CHECK(wdg(round_point3dc(5.0, 7.0, 4.0, 0.0), pl) != pscalar3dc(0.0));

        // bulk/weight split partition: the four parts sum to the object
        CHECK(round_bulk(d) + round_weight(d) + flat_bulk(d) + flat_weight(d) == d);
        CHECK(round_bulk(c) + round_weight(c) + flat_bulk(c) + flat_weight(c) == c);
        CHECK(round_weight_nrm_sq(s) == 1.0);
        CHECK(round_weight_nrm_sq(pl) == 0.0); // flat

        // conformal conjugate: involution + alternate center norm
        CHECK(cconj(cconj(d)) == d);
        CHECK(value_t(dot(unitize(a), cconj(unitize(a)))) == 14.0);

        // dot products of unitized rounds (reference Table 4.14): points
        // a1.a2 = -1/2 (v^2 + r1^2 + r2^2), spheres s1.s2 = +1/2 (v^2 - r1^2
        // - r2^2), null points give the distance identity -1/2 d^2
        auto a1 = round_point3dc(0.0, 0.0, 0.0, 1.0);
        auto a2 = round_point3dc(3.0, 0.0, 0.0, 2.0);
        CHECK(value_t(dot(a1, a2)) == -7.0);
        auto s1 = sphere3dc(0.0, 0.0, 0.0, 1.0);
        auto s2 = sphere3dc(3.0, 0.0, 0.0, 2.0);
        CHECK(value_t(dot(s1, s2)) == 2.0);
        CHECK(value_t(dot(round_point3dc(0.0, 0.0, 0.0, 0.0),
                          round_point3dc(3.0, 4.0, 0.0, 0.0))) == -12.5);
    }

    TEST_CASE("cga3dc: flat/round classification")
    {
        fmt::println("cga3dc: flat/round classification");

        // the grade does not decide the kind: each grade from 2 to 4 carries a
        // round and a flat object, and is_flat separates them
        auto rp = round_point3dc(1.0, 2.0, 3.0, 0.0);
        auto d = dipole3dc(1.0, 2.0, 3.0, 2.5, 0.0, 0.0, 1.0);
        auto c = circle3dc(1.0, 2.0, 3.0, 2.5, 0.0, 0.0, 1.0);
        auto s = sphere3dc(1.0, 2.0, 3.0, 2.5);
        auto fp = flat_point3dc(1.0, 2.0, 3.0);
        auto l = line3dc(1.0, 2.0, 3.0, 1.0, 0.0, 0.0);
        auto p = plane3dc(0.0, 0.0, 1.0, 4.0);

        CHECK(is_round(rp));
        CHECK(is_round(d)); // grade 2, round
        CHECK(is_flat(fp)); // grade 2, flat
        CHECK(is_round(c)); // grade 3, round
        CHECK(is_flat(l));  // grade 3, flat
        CHECK(is_round(s)); // grade 4, round
        CHECK(is_flat(p));  // grade 4, flat

        // the two are complementary
        CHECK(is_round(d) == !is_flat(d));
        CHECK(is_round(l) == !is_flat(l));

        // a COMPUTED flat is classified too: the carrier of a round point is
        // the flat point through it (4.25), and a transformed line stays flat
        CHECK(is_flat(car(rp)));
        CHECK(is_flat(transform(l, get_translation(1.0, 2.0, 3.0))));
        CHECK(
            is_round(transform(s, get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, pi / 3.0))));

        // the test is RELATIVE: the object's weight must not change the answer
        CHECK(is_flat(1.0e8 * l));
        CHECK(is_flat(1.0e-8 * l));
        CHECK(is_round(1.0e8 * s));

        // and it discriminates at the tolerance, not at exact zero: a line
        // contaminated by a fraction f of a circle is flat below eps_congruent
        // and round above it (falsification of the tolerance itself -- an
        // absolute test would answer differently at the two scales below)
        CHECK(is_flat(l + 1.0e-13 * c));
        CHECK(is_round(l + 1.0e-10 * c));
        CHECK(is_flat(1.0e8 * (l + 1.0e-13 * c)));
        CHECK(is_round(1.0e8 * (l + 1.0e-10 * c)));

        // an explicit tolerance moves the threshold
        CHECK(is_round(l + 1.0e-6 * c));
        CHECK(is_flat(l + 1.0e-6 * c, 1.0e-3));

        // a flat has no radius: radius_sq divides by the vanishing round weight
        CHECK(round_weight_nrm_sq(l) == 0.0);
        CHECK(round_weight_nrm_sq(p) == 0.0);
        CHECK(round_weight_nrm_sq(s) != 0.0);
    }

    TEST_CASE("cga3dc: directional-argument convention and vector overloads")
    {
        fmt::println("cga3dc: directional-argument convention and vector overloads");

        auto const p = vec3d(1.0, 2.0, 3.0);
        double const r = 2.0;

        // DIRECTION objects: the argument is the object's own direction, which
        // is pinned geometrically rather than by reading components -- the two
        // points of a dipole lie at center +/- r * direction
        for (auto const& d :
             {vec3d(1.0, 0.0, 0.0), vec3d(0.0, 1.0, 0.0), vec3d(0.0, 0.0, 1.0)}) {
            auto dp = dipole3dc(p.x, p.y, p.z, r, d.x, d.y, d.z);
            auto pts = dipole_points(dp);
            CHECK(pts.first.x == doctest::Approx(p.x + r * d.x));
            CHECK(pts.first.y == doctest::Approx(p.y + r * d.y));
            CHECK(pts.first.z == doctest::Approx(p.z + r * d.z));
            CHECK(pts.second.x == doctest::Approx(p.x - r * d.x));
            CHECK(pts.second.y == doctest::Approx(p.y - r * d.y));
            CHECK(pts.second.z == doctest::Approx(p.z - r * d.z));

            // a line's direction: points along it are incident with it
            auto l = line3dc(p.x, p.y, p.z, d.x, d.y, d.z);
            auto on =
                round_point3dc(p.x + 3.0 * d.x, p.y + 3.0 * d.y, p.z + 3.0 * d.z, 0.0);
            CHECK(is_close(mvec3dc(wdg(on, l)), mvec3dc(scalar3dc(0.0))));
        }

        // NORMAL objects: the hyperplane datum. A plane contains every point
        // it should, and a circle's plane is the one with the given normal
        auto n = vec3d(0.0, 0.0, 1.0);
        auto pl = plane3dc(n.x, n.y, n.z, 3.0);
        CHECK(is_close(mvec3dc(wdg(round_point3dc(5.0, -7.0, 3.0, 0.0), pl)),
                       mvec3dc(scalar3dc(0.0)))); // z = d = 3 lies in it
        auto ci = circle3dc(p.x, p.y, p.z, r, n.x, n.y, n.z);
        // a point of the circle: centre + r * (unit vector perpendicular to n)
        CHECK(is_close(mvec3dc(wdg(round_point3dc(p.x + r, p.y, p.z, 0.0), ci)),
                       mvec3dc(scalar3dc(0.0))));

        // the vector-taking overloads build exactly the component-form object
        auto d = vec3d(0.0, 1.0, 0.0);
        CHECK(round_point3dc(p, r) == round_point3dc(p.x, p.y, p.z, r));
        CHECK(flat_point3dc(p) == flat_point3dc(p.x, p.y, p.z));
        CHECK(sphere3dc(p, r) == sphere3dc(p.x, p.y, p.z, r));
        CHECK(line3dc(p, d) == line3dc(p.x, p.y, p.z, d.x, d.y, d.z));
        CHECK(plane3dc(n, 3.0) == plane3dc(n.x, n.y, n.z, 3.0));
        CHECK(circle3dc(p, r, n) == circle3dc(p.x, p.y, p.z, r, n.x, n.y, n.z));
        CHECK(dipole3dc(p, r, d) == dipole3dc(p.x, p.y, p.z, r, d.x, d.y, d.z));

        // and the accessors feed straight back into them (the reason for the
        // overloads): position() returns the vec3d a constructor takes
        auto s = sphere3dc(p, r);
        CHECK(sphere3dc(position(s), radius(s)) == s);
    }

    TEST_CASE("cga3dc: direction as one accessor over the object kinds")
    {
        fmt::println("cga3dc: direction as one accessor over the object kinds");

        auto const p = vec3d(1.0, 2.0, 3.0);
        auto const d = vec3d(0.6, 0.0, 0.8); // unit
        double const r = 2.0;

        auto same = [](vec3d const& a, vec3d const& b) {
            return a.x == doctest::Approx(b.x) && a.y == doctest::Approx(b.y) &&
                   a.z == doctest::Approx(b.z);
        };

        // ONE name over four kinds: the axis of a dipole, the direction of a
        // line, the normal of a circle's plane and of a plane -- each the
        // constructor's own directional argument, recovered
        CHECK(same(direction(dipole3dc(p, r, d)), d));
        CHECK(same(direction(line3dc(p, d)), d));
        CHECK(same(direction(circle3dc(p, r, d)), d)); // d is the plane normal
        CHECK(same(direction(plane3dc(d, 3.0)), d));

        // weight-independent: the result is a unit vector whatever the scale
        CHECK(same(direction(1.0e7 * dipole3dc(p, r, d)), d));
        CHECK(same(direction(1.0e-7 * line3dc(p, d)), d));
        CHECK(same(direction(1.0e7 * plane3dc(d, 3.0)), d));

        // both branches of the shared grade 3 are served, which is why is_flat
        // is consulted: a line gives a direction, a circle a normal
        auto l = line3dc(p, d);
        auto c = circle3dc(p, r, d);
        CHECK(is_flat(l));
        CHECK(is_round(c));
        CHECK(same(direction(l), direction(c))); // same vector, different meaning

        // the kinds that have none
        CHECK_THROWS_AS(direction(flat_point3dc(p)), std::runtime_error);
        CHECK_THROWS_AS(direction(sphere3dc(p, r)), std::runtime_error);

        // it survives a rigid motion the way a direction must: a rotation by
        // pi/2 about z carries x to y for every kind at once
        auto R = get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, pi / 2.0);
        auto const ex = vec3d(1.0, 0.0, 0.0), ey = vec3d(0.0, 1.0, 0.0);
        CHECK(same(direction(transform(dipole3dc(p, r, ex), R)), ey));
        CHECK(same(direction(transform(line3dc(p, ex), R)), ey));
        CHECK(same(direction(transform(plane3dc(ex, 3.0), R)), ey));
        // and a translation leaves every direction alone
        auto M = get_translation(4.0, -5.0, 6.0);
        CHECK(same(direction(transform(line3dc(p, d), M)), d));
        CHECK(same(direction(transform(circle3dc(p, r, d), M)), d));
    }

    TEST_CASE("cga3dc: Euclidean accessors - radius and position")
    {
        fmt::println("cga3dc: Euclidean accessors - radius and position");

        double const px = 1.0, py = 2.0, pz = 3.0, r = 2.5;
        auto rp = round_point3dc(px, py, pz, r);
        auto d = dipole3dc(px, py, pz, r, 0.0, 0.0, 1.0);
        auto c = circle3dc(px, py, pz, r, 0.0, 0.0, 1.0);
        auto s = sphere3dc(px, py, pz, r);
        auto fp = flat_point3dc(px, py, pz);

        // round trip against the constructor's own arguments, at every grade
        CHECK(radius(rp) == doctest::Approx(r));
        CHECK(radius(d) == doctest::Approx(r));
        CHECK(radius(c) == doctest::Approx(r));
        CHECK(radius(s) == doctest::Approx(r));

        // position returns an EGA vec3d; ga_cga.hpp carries the EGA types but
        // not the EGA operations, so compare componentwise here
        auto same = [](vec3d const& a, vec3d const& b) {
            return a.x == doctest::Approx(b.x) && a.y == doctest::Approx(b.y) &&
                   a.z == doctest::Approx(b.z);
        };
        auto const expect = vec3d(px, py, pz);
        CHECK(same(position(rp), expect));
        CHECK(same(position(d), expect));
        CHECK(same(position(c), expect));
        CHECK(same(position(s), expect));
        // the flat point has a position too, held in the flat part instead
        CHECK(same(position(fp), expect));

        // both are independent of the object's weight
        CHECK(radius(1.0e6 * s) == doctest::Approx(r));
        CHECK(same(position(1.0e6 * s), expect));
        CHECK(same(position(1.0e6 * fp), expect));

        // radius is the square root of the signed square, and the partner --
        // same centre, r^2 negated -- is exactly the imaginary case that throws
        CHECK(radius_sq(s) == doctest::Approx(r * r));
        CHECK(radius_sq(par(s)) == doctest::Approx(-r * r));
        CHECK_THROWS_AS(radius(par(s)), std::runtime_error);

        // a flat has no radius, and a line or plane has no centre
        auto l = line3dc(px, py, pz, 1.0, 0.0, 0.0);
        auto p = plane3dc(0.0, 0.0, 1.0, 4.0);
        CHECK_THROWS(radius(l));
        CHECK_THROWS(radius(p));
        CHECK_THROWS_AS(position(l), std::runtime_error);
        CHECK_THROWS_AS(position(p), std::runtime_error);

        // an unreachable two-link target makes the meet imaginary: radius_sq
        // is the sign test that decides BEFORE the square root is taken
        auto miss = rwdg(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(10.0, 0.0, 0.0, 1.0));
        CHECK(radius_sq(miss) < 0.0);
        CHECK_THROWS_AS(radius(miss), std::runtime_error);

        // the centre survives a rigid motion the way a position must
        auto M = get_translation(1.0, -2.0, 0.5);
        CHECK(same(position(transform(s, M)), vec3d(px + 1.0, py - 2.0, pz + 0.5)));
        CHECK(radius(transform(s, M)) == doctest::Approx(r));
        // and a dilation scales the radius by sigma while fixing its centre
        auto D = get_dilation(px, py, pz, 2.0);
        CHECK(radius(transform(s, D)) == doctest::Approx(2.0 * r));
        CHECK(same(position(transform(s, D)), expect));
    }

    TEST_CASE("cga3dc: off-center intersections with lines and planes")
    {
        fmt::println("cga3dc: off-center intersections with lines and planes");

        // The meet does not care whether the line or plane passes through the
        // center of the sphere. One chain covers secant, tangent and missing:
        // meet -> sign of radius_sq -> unitize -> cen +/- sqrt(r^2) * att.
        auto const s = sphere3dc(0.0, 0.0, 0.0, 2.0);

        // LINE along e1 at height y, i.e. at distance y from the center
        auto const chord = [&](value_t y, vec3dc& p1, vec3dc& p2) -> bool {
            auto d = rwdg(s, line3dc(0.0, y, 0.0, 1.0, 0.0, 0.0));
            if (radius_sq(d) < 0.0) return false;
            auto du = unitize(d);
            auto ctr = unitize(cen(du));
            auto dir = att(du);
            value_t const r = std::sqrt(radius_sq(du));
            p1 =
                vec3dc(ctr.x + r * dir.x, ctr.y + r * dir.y, ctr.z + r * dir.z, 1.0, 0.0);
            p2 =
                vec3dc(ctr.x - r * dir.x, ctr.y - r * dir.y, ctr.z - r * dir.z, 1.0, 0.0);
            return true;
        };

        vec3dc p1, p2;
        value_t const s3 = std::sqrt(3.0);

        // y = 1 (half the radius): chord ends at x = +/- sqrt(r^2 - y^2)
        CHECK(chord(1.0, p1, p2));
        bool const order_a = is_close(p1, vec3dc(s3, 1.0, 0.0, 1.0, 0.0)) &&
                             is_close(p2, vec3dc(-s3, 1.0, 0.0, 1.0, 0.0));
        bool const order_b = is_close(p1, vec3dc(-s3, 1.0, 0.0, 1.0, 0.0)) &&
                             is_close(p2, vec3dc(s3, 1.0, 0.0, 1.0, 0.0));
        CHECK((order_a || order_b));
        CHECK(radius_sq(rwdg(s, line3dc(0.0, 1.0, 0.0, 1.0, 0.0, 0.0))) == 3.0);

        // y = 2: tangent, the two points coincide and the radius vanishes
        CHECK(chord(2.0, p1, p2));
        CHECK(is_close(p1, vec3dc(0.0, 2.0, 0.0, 1.0, 0.0)));
        CHECK(is_close(p2, vec3dc(0.0, 2.0, 0.0, 1.0, 0.0)));

        // y = 3: the line misses the sphere -- imaginary radius, no real points
        CHECK_FALSE(chord(3.0, p1, p2));
        CHECK(radius_sq(rwdg(s, line3dc(0.0, 3.0, 0.0, 1.0, 0.0, 0.0))) == -5.0);

        // PLANE cutting off-center: the meet is a smaller circle, its radius
        // following r^2 - d^2, and its center sitting on the plane
        auto const cut = [&](value_t z) { return rwdg(s, plane3dc(0.0, 0.0, 1.0, z)); };
        CHECK(radius_sq(cut(1.0)) == 3.0); // r^2 - d^2 = 4 - 1
        CHECK(
            is_close(unitize(cen(unitize(cut(1.0)))), round_point3dc(0.0, 0.0, 1.0, s3)));
        CHECK(radius_sq(cut(2.0)) == 0.0);  // tangent plane
        CHECK(radius_sq(cut(3.0)) == -5.0); // plane misses the sphere

        // read the cut circle back with the Euclidean accessors: its plane is
        // the cutting plane, so its normal is that plane's normal, and its
        // centre is the foot of the perpendicular from the sphere's centre
        auto const same = [](vec3d const& a, vec3d const& b) {
            return a.x == doctest::Approx(b.x) && a.y == doctest::Approx(b.y) &&
                   a.z == doctest::Approx(b.z);
        };
        auto const c1 = cut(1.0);
        CHECK(radius(c1) == doctest::Approx(s3));
        CHECK(same(position(c1), vec3d(0.0, 0.0, 1.0)));
        CHECK(std::abs(direction(c1).z) == doctest::Approx(1.0)); // plane normal
        CHECK(direction(c1).x == doctest::Approx(0.0));
        CHECK(direction(c1).y == doctest::Approx(0.0));

        // the TANGENT plane touches at a single point: a zero-radius circle
        // sitting exactly on the sphere
        CHECK(radius(cut(2.0)) == doctest::Approx(0.0));
        CHECK(same(position(cut(2.0)), vec3d(0.0, 0.0, 2.0)));
        // and the missing plane has no real radius to take
        CHECK_THROWS_AS(radius(cut(3.0)), std::runtime_error);

        // the same for the LINE meets: the chord is a dipole whose axis is the
        // line's own direction and whose centre is the chord's midpoint
        auto const d1 = rwdg(s, line3dc(0.0, 1.0, 0.0, 1.0, 0.0, 0.0));
        CHECK(is_round(d1)); // a dipole, not a flat point
        CHECK(radius(d1) == doctest::Approx(s3));
        CHECK(same(position(d1), vec3d(0.0, 1.0, 0.0)));
        CHECK(std::abs(direction(d1).x) == doctest::Approx(1.0)); // along the line
        CHECK(direction(d1).y == doctest::Approx(0.0));
        CHECK(direction(d1).z == doctest::Approx(0.0));
        // tangent line: the dipole degenerates to the touch point
        auto const d2 = rwdg(s, line3dc(0.0, 2.0, 0.0, 1.0, 0.0, 0.0));
        CHECK(radius(d2) == doctest::Approx(0.0));
        CHECK(same(position(d2), vec3d(0.0, 2.0, 0.0)));
        // missing line: imaginary, so no real points
        CHECK_THROWS_AS(radius(rwdg(s, line3dc(0.0, 3.0, 0.0, 1.0, 0.0, 0.0))),
                        std::runtime_error);
        CHECK_THROWS_AS(dipole_points(rwdg(s, line3dc(0.0, 3.0, 0.0, 1.0, 0.0, 0.0))),
                        std::runtime_error);
    }

    TEST_CASE("cga3dc: sphere-sphere intersection end to end")
    {
        fmt::println("cga3dc: sphere-sphere intersection end to end");

        // The 3d counterpart of the cga2dc circle-circle case. Two spheres meet
        // in a CIRCLE (quadvector v quadvector -> trivector), and one chain
        // covers every configuration: meet -> sign of radius_sq -> the circle.
        // Nothing branches on whether they intersect; the radius decides.
        // (sphere v plane and sphere v line are the neighbouring case,
        // "off-center intersections with lines and planes".)
        auto const meet2 = [](quadvec3dc const& a, quadvec3dc const& b) {
            return rwdg(a, b);
        };
        auto const same = [](vec3d const& a, vec3d const& b) {
            return a.x == doctest::Approx(b.x) && a.y == doctest::Approx(b.y) &&
                   a.z == doctest::Approx(b.z);
        };

        // SECANT: two r=2 spheres at (0,0,0) and (2,0,0) meet in the circle of
        // radius sqrt(3) in the plane x = 1, whose normal is the centre line
        auto sec = meet2(sphere3dc(0.0, 0.0, 0.0, 2.0), sphere3dc(2.0, 0.0, 0.0, 2.0));
        CHECK(radius_sq(sec) == doctest::Approx(3.0)); // r^2 - d^2 = 4 - 1
        CHECK(radius(sec) == doctest::Approx(std::sqrt(3.0)));
        CHECK(same(position(sec), vec3d(1.0, 0.0, 0.0)));
        CHECK(is_round(sec)); // the meet of two spheres is a round circle
        // its plane's normal lies along the line of centres (sign is the
        // meet's orientation, so compare up to it)
        auto n = direction(sec);
        CHECK(std::abs(n.x) == doctest::Approx(1.0));
        CHECK(n.y == doctest::Approx(0.0));
        CHECK(n.z == doctest::Approx(0.0));
        // every point of that circle lies on BOTH spheres
        auto const on = round_point3dc(1.0, std::sqrt(3.0), 0.0, 0.0);
        CHECK(is_close(mvec3dc(wdg(on, sphere3dc(0.0, 0.0, 0.0, 2.0))),
                       mvec3dc(scalar3dc(0.0))));
        CHECK(is_close(mvec3dc(wdg(on, sphere3dc(2.0, 0.0, 0.0, 2.0))),
                       mvec3dc(scalar3dc(0.0))));

        // EXTERNALLY TANGENT: r=1 at the origin and r=2 at (3,0,0) touch at
        // (1,0,0) -- the circle degenerates to a point, radius exactly zero
        auto ext = meet2(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(3.0, 0.0, 0.0, 2.0));
        CHECK(radius_sq(ext) == doctest::Approx(0.0));
        CHECK(radius(ext) == doctest::Approx(0.0));
        CHECK(same(position(ext), vec3d(1.0, 0.0, 0.0)));

        // INTERNALLY TANGENT: r=3 at the origin and r=2 at (1,0,0) touch at
        // (3,0,0) -- also a zero-radius circle, from the inside
        auto in_ = meet2(sphere3dc(0.0, 0.0, 0.0, 3.0), sphere3dc(1.0, 0.0, 0.0, 2.0));
        CHECK(radius_sq(in_) == doctest::Approx(0.0));
        CHECK(same(position(in_), vec3d(3.0, 0.0, 0.0)));

        // DISJOINT: the meet still exists, with an imaginary radius -- this is
        // the sign test a caller makes before taking the square root
        auto dis = meet2(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(5.0, 0.0, 0.0, 1.0));
        CHECK(radius_sq(dis) < 0.0);
        CHECK(radius_sq(dis) == doctest::Approx(-5.25));
        CHECK_THROWS_AS(radius(dis), std::runtime_error);

        // ONE SPHERE INSIDE THE OTHER, not touching: also imaginary, so the
        // same sign test covers the second way of not intersecting
        auto nes = meet2(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(0.5, 0.0, 0.0, 3.0));
        CHECK(radius_sq(nes) < 0.0);
        CHECK_THROWS_AS(radius(nes), std::runtime_error);

        // CONCENTRIC: a degenerate meet. It is NOT flat -- it keeps a round
        // bulk -- yet it has no round weight at all, so it has neither a
        // radius nor a centre, and both accessors must refuse rather than
        // divide by zero (the guard is on the round weight, not on is_flat)
        auto con2 = meet2(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(0.0, 0.0, 0.0, 2.0));
        CHECK_FALSE(is_flat(con2));
        CHECK(round_weight_nrm_sq(con2) == 0.0);
        // the MESSAGE is pinned, not just the throw: this build defines
        // _HD_GA_EXTENDED_TEST_DIV_BY_ZERO, so radius_sq's division check would
        // throw here even with the old, narrower is_flat guard -- and a bare
        // CHECK_THROWS therefore passes whether or not the round-weight guard
        // exists. Naming the message is what actually pins it.
        CHECK_THROWS_WITH_AS(radius_sq(con2), doctest::Contains("no round weight"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS(radius(con2), doctest::Contains("no round weight"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS(position(con2), doctest::Contains("no round weight"),
                             std::runtime_error);
        // is_degenerate names that condition, and is INDEPENDENT of flat/round:
        // the concentric meet is round AND degenerate
        CHECK(is_degenerate(con2));
        CHECK(is_round(con2));
        // every ordinary meet above is non-degenerate, whatever its radius sign
        CHECK_FALSE(is_degenerate(sec));
        CHECK_FALSE(is_degenerate(ext));
        CHECK_FALSE(is_degenerate(dis));
        // a flat is degenerate too (no round weight at all), so the two
        // questions overlap without either implying the other
        CHECK(is_degenerate(plane3dc(0.0, 0.0, 1.0, 1.0)));
        CHECK(is_flat(plane3dc(0.0, 0.0, 1.0, 1.0)));

        // the reason it is a FUNCTION and not `round_weight_nrm_sq(m) == 0.0`
        // at the call site: a meet of NEARLY concentric spheres leaves a tiny
        // but non-zero weight, which an exact test misses. The relative test
        // catches it, and widening rel_tol moves the threshold.
        auto const near_ =
            meet2(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(1.0e-16, 0.0, 0.0, 2.0));
        CHECK(round_weight_nrm_sq(near_) > 0.0); // an exact == 0.0 test FAILS here
        CHECK(is_degenerate(near_));             // the relative one does not
        // further apart, it is an ordinary (here: non-intersecting) pair again
        auto const far_ =
            meet2(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(0.1, 0.0, 0.0, 2.0));
        CHECK_FALSE(is_degenerate(far_));
        CHECK(radius_sq(far_) < 0.0);    // nested, so no real intersection
        CHECK(is_degenerate(far_, 0.5)); // a wide tolerance calls it degenerate

        // the documented three-way ladder, exercised as a caller would write it
        auto const classify = [](quadvec3dc const& a, quadvec3dc const& b) -> int {
            auto m = rwdg(a, b);
            if (is_degenerate(m)) return 0;   // concentric: no centre, no radius
            if (radius_sq(m) < 0.0) return 1; // no real intersection
            return 2;                         // real circle or tangent point
        };
        CHECK(classify(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(0.0, 0.0, 0.0, 2.0)) ==
              0);
        CHECK(classify(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(5.0, 0.0, 0.0, 1.0)) ==
              1);
        CHECK(classify(sphere3dc(0.0, 0.0, 0.0, 2.0), sphere3dc(2.0, 0.0, 0.0, 2.0)) ==
              2);
        CHECK(classify(sphere3dc(0.0, 0.0, 0.0, 1.0), sphere3dc(3.0, 0.0, 0.0, 2.0)) ==
              2); // tangent

        // the whole ladder is monotone in the centre distance: as the second
        // sphere moves away, the meet radius shrinks through zero into the
        // imaginary range, with no branch anywhere in the chain
        value_t prev = 1.0e30;
        for (value_t d : {0.0, 1.0, 2.0, 3.0, 4.0, 5.0}) {
            auto m = meet2(sphere3dc(0.0, 0.0, 0.0, 2.0), sphere3dc(d, 0.0, 0.0, 2.0));
            if (d == 0.0) continue; // concentric: no radius (checked above)
            value_t const rs = radius_sq(m);
            CHECK(rs < prev); // strictly decreasing
            prev = rs;
        }
        CHECK(prev < 0.0); // ends imaginary, i.e. separated
    }

    TEST_CASE("cga3dc: dipole points and two-link inverse kinematics")
    {
        fmt::println("cga3dc: dipole points and two-link inverse kinematics");

        // the two points of a dipole built directly: center (1,2,3), r = 2, axis x
        {
            auto const [p1, p2] =
                dipole_points(dipole3dc(1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 0.0));
            bool const oa = is_close(p1, round_point3dc(3.0, 2.0, 3.0, 0.0)) &&
                            is_close(p2, round_point3dc(-1.0, 2.0, 3.0, 0.0));
            bool const ob = is_close(p1, round_point3dc(-1.0, 2.0, 3.0, 0.0)) &&
                            is_close(p2, round_point3dc(3.0, 2.0, 3.0, 0.0));
            CHECK((oa || ob));
        }

        // two-link IK as sphere v sphere v plane: base at the origin, tip (2,0,0),
        // links of length 2; in the plane z = 0 the joint is (1, +/- sqrt 3, 0), in
        // the plane y = 0 it is (1, 0, +/- sqrt 3) -- the sphere meet is the whole
        // circle of joint positions, the plane picks two of them
        value_t const s3 = std::sqrt(3.0);
        {
            auto const [j1, j2] = two_link_ik3dc(0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0,
                                                 plane3dc(0.0, 0.0, 1.0, 0.0));
            bool const ka = is_close(j1, round_point3dc(1.0, s3, 0.0, 0.0)) &&
                            is_close(j2, round_point3dc(1.0, -s3, 0.0, 0.0));
            bool const kb = is_close(j1, round_point3dc(1.0, -s3, 0.0, 0.0)) &&
                            is_close(j2, round_point3dc(1.0, s3, 0.0, 0.0));
            CHECK((ka || kb));
            auto const [m1, m2] = two_link_ik3dc(0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0,
                                                 plane3dc(0.0, 1.0, 0.0, 0.0));
            for (auto const& m : {m1, m2}) {
                CHECK(m.x / m.w == doctest::Approx(1.0));
                CHECK(m.y / m.w == doctest::Approx(0.0));
                CHECK(std::abs(m.z / m.w) == doctest::Approx(s3));
            }
            // a general case checked by the link lengths and the plane alone
            auto const [k1, k2] = two_link_ik3dc(0.2, -0.1, 0.3, 0.9, 0.7, 0.8, 0.8, 0.6,
                                                 plane3dc(0.0, 0.0, 1.0, 0.5));
            for (auto const& k : {k1, k2}) {
                value_t const x = k.x / k.w, y = k.y / k.w, z = k.z / k.w;
                CHECK(std::hypot(x - 0.2, y + 0.1, z - 0.3) == doctest::Approx(0.8));
                CHECK(std::hypot(x - 0.9, y - 0.7, z - 0.8) == doctest::Approx(0.6));
                CHECK(z == doctest::Approx(0.5));
            }
            // out of reach: the dipole is imaginary
            CHECK_THROWS(two_link_ik3dc(0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 2.0, 2.0,
                                        plane3dc(0.0, 0.0, 1.0, 0.0)));
        }
    }

    TEST_CASE("cga3dc: object properties (car/ccr/cen/con/par/att)")
    {
        fmt::println("cga3dc: object properties (car/ccr/cen/con/par/att)");

        auto a = round_point3dc(1.0, 2.0, 3.0, 2.0);
        auto s = sphere3dc(1.0, 2.0, 3.0, 2.0);
        auto c = circle3dc(1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 1.0);
        auto d = dipole3dc(1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 0.0);

        // centers: every object about (1,2,3) with r = 2 has the SAME round
        // center -- exact pins (orientation-sensitive, lesson: congruence
        // gates are sign-blind)
        CHECK(cen(a) == a);
        CHECK(cen(unitize(s)) == a);
        CHECK(cen(unitize(c)) == a);
        CHECK(cen(unitize(d)) == a);

        // attitudes recover the construction directions with POSITIVE sign
        CHECK(att(a) == scalar3dc(1.0));
        CHECK(att(unitize(d)) == vec3dc(1.0, 0.0, 0.0, 0.0, 1.0));
        CHECK(att(unitize(c)) ==
              bivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, -1.0, 0.0, 0.0));
        CHECK(att(line3dc(1.0, 2.0, 3.0, 1.0, 0.0, 0.0)) ==
              bivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0));
        CHECK(att(plane3dc(0.0, 0.0, 1.0, 3.0)) ==
              trivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0));
        // Table 4.10 rows for the remaining types: flat point att = pw e5,
        // sphere att = su e321 + sx e235 + sy e315 + sz e125
        CHECK(att(flat_point3dc(1.0, 2.0, 3.0)) == vec3dc(0.0, 0.0, 0.0, 0.0, 1.0));
        CHECK(att(s) == trivec3dc(0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, -1.0));

        // carriers: car of a round point IS the flat point at its center
        // (reference 4.25); car of a flat is zero; car of the circle is
        // (congruent to) its plane, car of the dipole its line
        CHECK(car(a) == flat_point3dc(1.0, 2.0, 3.0));
        CHECK(car(flat_point3dc(1.0, 2.0, 3.0)) ==
              trivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        CHECK(car(line3dc(1.0, 2.0, 3.0, 1.0, 0.0, 0.0)) ==
              quadvec3dc(0.0, 0.0, 0.0, 0.0, 0.0));
        CHECK(is_congruent(car(c), plane3dc(0.0, 0.0, 1.0, 3.0)));
        CHECK(is_congruent(car(d), line3dc(1.0, 2.0, 3.0, 1.0, 0.0, 0.0)));

        // containers: smallest containing sphere; the round point's container
        // is pinned exactly above (== sphere3dc); circle/dipole congruent
        CHECK(is_congruent(con(c), s));
        CHECK(is_congruent(con(d), s));

        // factorization u = car(u) v con(u) (reference 4.63), up to scale
        CHECK(is_congruent(rwdg(car(c), con(c)), c));
        CHECK(is_congruent(rwdg(car(d), con(d)), d));

        // partners: same center and carrier, r^2 negated; orthogonal to u;
        // exact pin for the round point
        CHECK(par(a) == vec3dc(1.0, 2.0, 3.0, 1.0, 5.0));
        CHECK(radius_sq(par(a)) == -4.0);
        CHECK(radius_sq(unitize(par(c))) == -4.0);
        CHECK(value_t(dot(a, par(a))) == 0.0);
        CHECK(value_t(dot(c, par(c))) == 0.0);

        // dipole surface points: p+- = cen +- rsqrt(radius_sq) * att (exact:
        // the two round points at center +- r * axis)
        auto du = unitize(d);
        auto p_plus = vec3dc(cen(du) + std::sqrt(radius_sq(du)) * att(du));
        auto p_minus = vec3dc(cen(du) - std::sqrt(radius_sq(du)) * att(du));
        CHECK(p_plus == round_point3dc(3.0, 2.0, 3.0, 0.0));
        CHECK(p_minus == round_point3dc(-1.0, 2.0, 3.0, 0.0));

        // u ^ cen(u) = 0 (the center lies on every centered object's carrier
        // chain -- reference 4.58)
        CHECK(wdg(a, cen(a)) ==
              bivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    TEST_CASE("cga3dc: exp/log/sqrt and conformal transformations")
    {
        fmt::println("cga3dc: exp/log/sqrt and conformal transformations");

        auto const I_u = mvec3dc_u(
            vec3dc(0.0, 0.0, 0.0, 0.0, 0.0),
            trivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), pscalar3dc(1.0));
        auto euclid = [](vec3dc const& q) {
            auto u = unitize(q);
            return vec3dc(u.x, u.y, u.z, 0.0, 0.0);
        };

        // rexp(0) = I; unit versors for all three builders
        CHECK(rexp(trivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)) == I_u);
        auto Tm = get_translation(2.0, 0.0, 0.0);
        auto Rm = get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, pi / 2.0);
        auto Dm = get_dilation(0.0, 0.0, 0.0, 4.0);
        CHECK(is_close(rgpr(Tm, rrev(Tm)), I_u));
        CHECK(is_close(rgpr(Rm, rrev(Rm)), I_u));
        CHECK(is_close(rgpr(Dm, rrev(Dm)), I_u));

        // translation is EXACT (parabolic, series terminates)
        CHECK(transform(round_point3dc(1.0, 2.0, 3.0, 0.0), Tm) ==
              round_point3dc(3.0, 2.0, 3.0, 0.0));
        // matches the two-plane reflection composition by construction; the
        // generator is exactly parabolic
        CHECK(rgpr(gr3(Tm), gr3(Tm)) ==
              mvec3dc_u(vec3dc(0.0, 0.0, 0.0, 0.0, 0.0),
                        trivec3dc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                        pscalar3dc(0.0)));

        // rotation: right-handed about +z; radius invariant under rigid motion
        CHECK(is_close(euclid(transform(round_point3dc(1.0, 0.0, 0.0, 0.0), Rm)),
                       vec3dc(0.0, 1.0, 0.0, 0.0, 0.0)));
        auto s = sphere3dc(1.0, 2.0, 3.0, 2.0);
        CHECK(std::abs(radius_sq(unitize(transform(s, Rm))) - 4.0) < 1e-12);

        // dilation: sigma = 4 about the origin; radius_sq scales by sigma^2
        CHECK(is_close(euclid(transform(round_point3dc(1.0, 0.0, 0.0, 0.0), Dm)),
                       vec3dc(4.0, 0.0, 0.0, 0.0, 0.0)));
        CHECK(std::abs(radius_sq(unitize(transform(sphere3dc(0.0, 0.0, 0.0, 1.0), Dm))) -
                       16.0) < 1e-12);
        CHECK_THROWS(get_dilation(0.0, 0.0, 0.0, 0.0));

        // inversion on the unit sphere: reciprocal distance, fixed sphere
        // points, involution
        auto S1 = sphere3dc(0.0, 0.0, 0.0, 1.0);
        CHECK(is_close(euclid(invert_on(round_point3dc(2.0, 0.0, 0.0, 0.0), S1)),
                       vec3dc(0.5, 0.0, 0.0, 0.0, 0.0)));
        CHECK(is_close(euclid(invert_on(round_point3dc(1.0, 0.0, 0.0, 0.0), S1)),
                       vec3dc(1.0, 0.0, 0.0, 0.0, 0.0)));
        auto q = round_point3dc(1.0, 2.0, 3.0, 0.0);
        CHECK(is_close(unitize(invert_on(invert_on(q, S1), S1)), unitize(q)));

        // log/exp round trip (elliptic + hyperbolic + parabolic) and sqrt
        CHECK(is_close(rexp(rlog(Rm)), Rm));
        CHECK(is_close(rexp(rlog(Dm)), Dm));
        CHECK(is_close(rexp(rlog(Tm)), Tm));
        CHECK(is_close(rgpr(rsqrt(Rm), rsqrt(Rm)), Rm));
        CHECK(is_close(rgpr(rsqrt(Dm), rsqrt(Dm)), Dm));

        // one-parameter group: R(a) v R(b) = R(a+b)
        auto Ra = get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.4);
        auto Rb = get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.9);
        CHECK(is_close(rgpr(Ra, Rb), get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.3)));

        // conjugation identity: T(m) R T(-m) = rotation about the shifted line
        // (compared by ACTION -- versors double-cover, the conjugate may land
        // on the -M sheet)
        auto Tsh = get_translation(1.0, 2.0, 0.0);
        auto Rsh = get_rotation(1.0, 2.0, 0.0, 0.0, 0.0, 1.0, pi / 2.0);
        auto Mc = rgpr(rgpr(Tsh, Rm), rrev(Tsh));
        auto qc = round_point3dc(2.0, 2.0, 0.0, 0.0);
        CHECK(is_close(euclid(transform(qc, Mc)), euclid(transform(qc, Rsh))));
        CHECK(is_close(euclid(transform(qc, Rsh)), vec3dc(1.0, 3.0, 0.0, 0.0, 0.0)));

        // non-simple (screw) generators throw; composed motors still work
        auto lz = gr3(get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0));
        auto tz = gr3(get_translation(0.0, 0.0, 2.0));
        CHECK_THROWS(rexp(trivec3dc(lz + tz)));
        auto screw = rgpr(get_translation(0.0, 0.0, 2.0), Rm); // compose instead
        CHECK(is_close(rgpr(screw, rrev(screw)), I_u));
        CHECK(is_close(euclid(transform(round_point3dc(1.0, 0.0, 0.0, 0.0), screw)),
                       vec3dc(0.0, 1.0, 2.0, 0.0, 0.0)));
        CHECK_THROWS(rlog(screw)); // outside the simple-generator image

        // is_same_transform: action-equality across the double cover and
        // uniform rescaling; distinguishes genuinely different maps
        CHECK(is_same_transform(Rm, -1.0 * Rm));
        CHECK(is_same_transform(Tm, 2.5 * Tm));
        CHECK(is_same_transform(
            Rm, rgpr(get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, pi / 4.0),
                     get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, pi / 4.0))));
        // translation equals its two-plane reflection composition
        CHECK(
            is_same_transform(Tm, rgpr(mvec3dc_u(gr1(rgpr(plane3dc(1.0, 0.0, 0.0, 1.0),
                                                          plane3dc(1.0, 0.0, 0.0, 0.0))),
                                                 gr3(rgpr(plane3dc(1.0, 0.0, 0.0, 1.0),
                                                          plane3dc(1.0, 0.0, 0.0, 0.0))),
                                                 gr5(rgpr(plane3dc(1.0, 0.0, 0.0, 1.0),
                                                          plane3dc(1.0, 0.0, 0.0, 0.0)))),
                                       I_u)));
        CHECK_FALSE(is_same_transform(Rm, Tm));
        CHECK_FALSE(is_same_transform(Rm, Dm));
        // transversion: x' = (x + x^2 t)/(1 + 2 t.x + t^2 x^2); on the t-axis
        // x -> x/(1 + tau x)
        auto V = get_transversion(0.5, 0.0, 0.0);
        CHECK(is_close(rgpr(V, rrev(V)), I_u));
        CHECK(is_close(euclid(transform(round_point3dc(1.0, 0.0, 0.0, 0.0), V)),
                       vec3dc(2.0 / 3.0, 0.0, 0.0, 0.0, 0.0)));
        CHECK(is_close(euclid(transform(round_point3dc(2.0, 0.0, 0.0, 0.0), V)),
                       vec3dc(1.0, 0.0, 0.0, 0.0, 0.0)));

        // loxodromic about the dipole with points (0,0,+-1): fixes both points;
        // the pure hyperbolic flow moves the equator point (1,0,0) to
        // (sech(2 delta), 0, tanh(2 delta)); the pure elliptic part is the
        // rotation about the carrier line; the two parts commute
        auto dz = dipole3dc(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0);
        auto Lx = get_loxodromic(dz, 0.7, 0.3);
        CHECK(is_close(rgpr(Lx, rrev(Lx)), I_u));
        CHECK(is_congruent(transform(round_point3dc(0.0, 0.0, 1.0, 0.0), Lx),
                           round_point3dc(0.0, 0.0, 1.0, 0.0)));
        CHECK(is_congruent(transform(round_point3dc(0.0, 0.0, -1.0, 0.0), Lx),
                           round_point3dc(0.0, 0.0, -1.0, 0.0)));
        auto Lh = get_loxodromic(dz, 0.0, 0.3);
        CHECK(is_close(euclid(transform(round_point3dc(1.0, 0.0, 0.0, 0.0), Lh)),
                       vec3dc(1.0 / std::cosh(0.6), 0.0, std::tanh(0.6), 0.0, 0.0)));
        CHECK(is_same_transform(get_loxodromic(dz, 0.7, 0.0),
                                get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.7)));
        CHECK(is_same_transform(Lx, rgpr(get_loxodromic(dz, 0.7, 0.0), Lh)));
        CHECK(is_same_transform(Lx, rgpr(Lh, get_loxodromic(dz, 0.7, 0.0))));

        // conjugated rotation == rotation about the shifted line (the check
        // the double cover blocked above)
        CHECK(is_same_transform(Mc, Rsh));
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

    TEST_CASE("MVec3dc: +/- closure over every type pair")
    {
        fmt::println("MVec3dc: +/- closure over every type pair");

        // Every pair of cga3dc types must add and subtract, and the result must
        // equal the same sum carried out in the full multivector. The operators
        // themselves keep the NARROWEST type that holds the grade union; widening
        // here is only so one reference expression can check them all.

        auto const Scalar3dc_ = scalar3dc(2.0);
        auto const Vec3dc_ = vec3dc(1.0, 2.0, 3.0, 4.0, 5.0);
        auto const BiVec3dc_ =
            bivec3dc(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
        auto const TriVec3dc_ =
            trivec3dc(2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0);
        auto const QuadVec3dc_ = quadvec3dc(3.0, 4.0, 5.0, 6.0, 7.0);
        auto const PScalar3dc_ = pscalar3dc(-13.0);
        auto const MVec3dc_E_ = mvec3dc_e(
            scalar3dc(1.5), bivec3dc(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0),
            quadvec3dc(2.0, 3.0, 4.0, 5.0, 6.0));
        auto const MVec3dc_U_ =
            mvec3dc_u(vec3dc(1.0, 2.0, 3.0, 4.0, 5.0),
                      trivec3dc(2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0),
                      pscalar3dc(7.5));
        auto const MVec3dc_ =
            mvec3dc(scalar3dc(0.5), vec3dc(1.0, 2.0, 3.0, 4.0, 5.0),
                    bivec3dc(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0),
                    trivec3dc(2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0),
                    quadvec3dc(3.0, 4.0, 5.0, 6.0, 7.0), pscalar3dc(9.0));

        CHECK(mvec3dc(Scalar3dc_ + Scalar3dc_) ==
              mvec3dc(Scalar3dc_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(Scalar3dc_ - Scalar3dc_) ==
              mvec3dc(Scalar3dc_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(Scalar3dc_ + Vec3dc_) == mvec3dc(Scalar3dc_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(Scalar3dc_ - Vec3dc_) == mvec3dc(Scalar3dc_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(Scalar3dc_ + BiVec3dc_) ==
              mvec3dc(Scalar3dc_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(Scalar3dc_ - BiVec3dc_) ==
              mvec3dc(Scalar3dc_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(Scalar3dc_ + TriVec3dc_) ==
              mvec3dc(Scalar3dc_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(Scalar3dc_ - TriVec3dc_) ==
              mvec3dc(Scalar3dc_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(Scalar3dc_ + QuadVec3dc_) ==
              mvec3dc(Scalar3dc_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(Scalar3dc_ - QuadVec3dc_) ==
              mvec3dc(Scalar3dc_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(Scalar3dc_ + PScalar3dc_) ==
              mvec3dc(Scalar3dc_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(Scalar3dc_ - PScalar3dc_) ==
              mvec3dc(Scalar3dc_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(Scalar3dc_ + MVec3dc_E_) ==
              mvec3dc(Scalar3dc_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(Scalar3dc_ - MVec3dc_E_) ==
              mvec3dc(Scalar3dc_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(Scalar3dc_ + MVec3dc_U_) ==
              mvec3dc(Scalar3dc_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(Scalar3dc_ - MVec3dc_U_) ==
              mvec3dc(Scalar3dc_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(Scalar3dc_ + MVec3dc_) == mvec3dc(Scalar3dc_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(Scalar3dc_ - MVec3dc_) == mvec3dc(Scalar3dc_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(Vec3dc_ + Scalar3dc_) == mvec3dc(Vec3dc_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(Vec3dc_ - Scalar3dc_) == mvec3dc(Vec3dc_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(Vec3dc_ + Vec3dc_) == mvec3dc(Vec3dc_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(Vec3dc_ - Vec3dc_) == mvec3dc(Vec3dc_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(Vec3dc_ + BiVec3dc_) == mvec3dc(Vec3dc_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(Vec3dc_ - BiVec3dc_) == mvec3dc(Vec3dc_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(Vec3dc_ + TriVec3dc_) == mvec3dc(Vec3dc_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(Vec3dc_ - TriVec3dc_) == mvec3dc(Vec3dc_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(Vec3dc_ + QuadVec3dc_) == mvec3dc(Vec3dc_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(Vec3dc_ - QuadVec3dc_) == mvec3dc(Vec3dc_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(Vec3dc_ + PScalar3dc_) == mvec3dc(Vec3dc_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(Vec3dc_ - PScalar3dc_) == mvec3dc(Vec3dc_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(Vec3dc_ + MVec3dc_E_) == mvec3dc(Vec3dc_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(Vec3dc_ - MVec3dc_E_) == mvec3dc(Vec3dc_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(Vec3dc_ + MVec3dc_U_) == mvec3dc(Vec3dc_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(Vec3dc_ - MVec3dc_U_) == mvec3dc(Vec3dc_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(Vec3dc_ + MVec3dc_) == mvec3dc(Vec3dc_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(Vec3dc_ - MVec3dc_) == mvec3dc(Vec3dc_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ + Scalar3dc_) ==
              mvec3dc(BiVec3dc_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(BiVec3dc_ - Scalar3dc_) ==
              mvec3dc(BiVec3dc_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(BiVec3dc_ + Vec3dc_) == mvec3dc(BiVec3dc_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(BiVec3dc_ - Vec3dc_) == mvec3dc(BiVec3dc_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(BiVec3dc_ + BiVec3dc_) == mvec3dc(BiVec3dc_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ - BiVec3dc_) == mvec3dc(BiVec3dc_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ + TriVec3dc_) ==
              mvec3dc(BiVec3dc_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ - TriVec3dc_) ==
              mvec3dc(BiVec3dc_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ + QuadVec3dc_) ==
              mvec3dc(BiVec3dc_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ - QuadVec3dc_) ==
              mvec3dc(BiVec3dc_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ + PScalar3dc_) ==
              mvec3dc(BiVec3dc_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(BiVec3dc_ - PScalar3dc_) ==
              mvec3dc(BiVec3dc_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(BiVec3dc_ + MVec3dc_E_) ==
              mvec3dc(BiVec3dc_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(BiVec3dc_ - MVec3dc_E_) ==
              mvec3dc(BiVec3dc_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(BiVec3dc_ + MVec3dc_U_) ==
              mvec3dc(BiVec3dc_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(BiVec3dc_ - MVec3dc_U_) ==
              mvec3dc(BiVec3dc_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(BiVec3dc_ + MVec3dc_) == mvec3dc(BiVec3dc_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(BiVec3dc_ - MVec3dc_) == mvec3dc(BiVec3dc_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ + Scalar3dc_) ==
              mvec3dc(TriVec3dc_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(TriVec3dc_ - Scalar3dc_) ==
              mvec3dc(TriVec3dc_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(TriVec3dc_ + Vec3dc_) == mvec3dc(TriVec3dc_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(TriVec3dc_ - Vec3dc_) == mvec3dc(TriVec3dc_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(TriVec3dc_ + BiVec3dc_) ==
              mvec3dc(TriVec3dc_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ - BiVec3dc_) ==
              mvec3dc(TriVec3dc_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ + TriVec3dc_) ==
              mvec3dc(TriVec3dc_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ - TriVec3dc_) ==
              mvec3dc(TriVec3dc_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ + QuadVec3dc_) ==
              mvec3dc(TriVec3dc_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ - QuadVec3dc_) ==
              mvec3dc(TriVec3dc_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ + PScalar3dc_) ==
              mvec3dc(TriVec3dc_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(TriVec3dc_ - PScalar3dc_) ==
              mvec3dc(TriVec3dc_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(TriVec3dc_ + MVec3dc_E_) ==
              mvec3dc(TriVec3dc_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(TriVec3dc_ - MVec3dc_E_) ==
              mvec3dc(TriVec3dc_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(TriVec3dc_ + MVec3dc_U_) ==
              mvec3dc(TriVec3dc_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(TriVec3dc_ - MVec3dc_U_) ==
              mvec3dc(TriVec3dc_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(TriVec3dc_ + MVec3dc_) == mvec3dc(TriVec3dc_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(TriVec3dc_ - MVec3dc_) == mvec3dc(TriVec3dc_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ + Scalar3dc_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(QuadVec3dc_ - Scalar3dc_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(QuadVec3dc_ + Vec3dc_) == mvec3dc(QuadVec3dc_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ - Vec3dc_) == mvec3dc(QuadVec3dc_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ + BiVec3dc_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ - BiVec3dc_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ + TriVec3dc_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ - TriVec3dc_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ + QuadVec3dc_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ - QuadVec3dc_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ + PScalar3dc_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(QuadVec3dc_ - PScalar3dc_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(QuadVec3dc_ + MVec3dc_E_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(QuadVec3dc_ - MVec3dc_E_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(QuadVec3dc_ + MVec3dc_U_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(QuadVec3dc_ - MVec3dc_U_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(QuadVec3dc_ + MVec3dc_) ==
              mvec3dc(QuadVec3dc_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(QuadVec3dc_ - MVec3dc_) ==
              mvec3dc(QuadVec3dc_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ + Scalar3dc_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(PScalar3dc_ - Scalar3dc_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(PScalar3dc_ + Vec3dc_) == mvec3dc(PScalar3dc_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(PScalar3dc_ - Vec3dc_) == mvec3dc(PScalar3dc_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(PScalar3dc_ + BiVec3dc_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ - BiVec3dc_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ + TriVec3dc_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ - TriVec3dc_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ + QuadVec3dc_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ - QuadVec3dc_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ + PScalar3dc_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(PScalar3dc_ - PScalar3dc_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(PScalar3dc_ + MVec3dc_E_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(PScalar3dc_ - MVec3dc_E_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(PScalar3dc_ + MVec3dc_U_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(PScalar3dc_ - MVec3dc_U_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(PScalar3dc_ + MVec3dc_) ==
              mvec3dc(PScalar3dc_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(PScalar3dc_ - MVec3dc_) ==
              mvec3dc(PScalar3dc_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ + Scalar3dc_) ==
              mvec3dc(MVec3dc_E_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(MVec3dc_E_ - Scalar3dc_) ==
              mvec3dc(MVec3dc_E_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(MVec3dc_E_ + Vec3dc_) == mvec3dc(MVec3dc_E_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ - Vec3dc_) == mvec3dc(MVec3dc_E_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ + BiVec3dc_) ==
              mvec3dc(MVec3dc_E_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ - BiVec3dc_) ==
              mvec3dc(MVec3dc_E_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ + TriVec3dc_) ==
              mvec3dc(MVec3dc_E_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ - TriVec3dc_) ==
              mvec3dc(MVec3dc_E_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ + QuadVec3dc_) ==
              mvec3dc(MVec3dc_E_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ - QuadVec3dc_) ==
              mvec3dc(MVec3dc_E_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ + PScalar3dc_) ==
              mvec3dc(MVec3dc_E_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(MVec3dc_E_ - PScalar3dc_) ==
              mvec3dc(MVec3dc_E_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(MVec3dc_E_ + MVec3dc_E_) ==
              mvec3dc(MVec3dc_E_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(MVec3dc_E_ - MVec3dc_E_) ==
              mvec3dc(MVec3dc_E_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(MVec3dc_E_ + MVec3dc_U_) ==
              mvec3dc(MVec3dc_E_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(MVec3dc_E_ - MVec3dc_U_) ==
              mvec3dc(MVec3dc_E_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(MVec3dc_E_ + MVec3dc_) == mvec3dc(MVec3dc_E_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(MVec3dc_E_ - MVec3dc_) == mvec3dc(MVec3dc_E_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ + Scalar3dc_) ==
              mvec3dc(MVec3dc_U_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(MVec3dc_U_ - Scalar3dc_) ==
              mvec3dc(MVec3dc_U_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(MVec3dc_U_ + Vec3dc_) == mvec3dc(MVec3dc_U_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ - Vec3dc_) == mvec3dc(MVec3dc_U_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ + BiVec3dc_) ==
              mvec3dc(MVec3dc_U_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ - BiVec3dc_) ==
              mvec3dc(MVec3dc_U_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ + TriVec3dc_) ==
              mvec3dc(MVec3dc_U_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ - TriVec3dc_) ==
              mvec3dc(MVec3dc_U_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ + QuadVec3dc_) ==
              mvec3dc(MVec3dc_U_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ - QuadVec3dc_) ==
              mvec3dc(MVec3dc_U_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ + PScalar3dc_) ==
              mvec3dc(MVec3dc_U_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(MVec3dc_U_ - PScalar3dc_) ==
              mvec3dc(MVec3dc_U_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(MVec3dc_U_ + MVec3dc_E_) ==
              mvec3dc(MVec3dc_U_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(MVec3dc_U_ - MVec3dc_E_) ==
              mvec3dc(MVec3dc_U_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(MVec3dc_U_ + MVec3dc_U_) ==
              mvec3dc(MVec3dc_U_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(MVec3dc_U_ - MVec3dc_U_) ==
              mvec3dc(MVec3dc_U_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(MVec3dc_U_ + MVec3dc_) == mvec3dc(MVec3dc_U_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(MVec3dc_U_ - MVec3dc_) == mvec3dc(MVec3dc_U_) - mvec3dc(MVec3dc_));
        CHECK(mvec3dc(MVec3dc_ + Scalar3dc_) == mvec3dc(MVec3dc_) + mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(MVec3dc_ - Scalar3dc_) == mvec3dc(MVec3dc_) - mvec3dc(Scalar3dc_));
        CHECK(mvec3dc(MVec3dc_ + Vec3dc_) == mvec3dc(MVec3dc_) + mvec3dc(Vec3dc_));
        CHECK(mvec3dc(MVec3dc_ - Vec3dc_) == mvec3dc(MVec3dc_) - mvec3dc(Vec3dc_));
        CHECK(mvec3dc(MVec3dc_ + BiVec3dc_) == mvec3dc(MVec3dc_) + mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(MVec3dc_ - BiVec3dc_) == mvec3dc(MVec3dc_) - mvec3dc(BiVec3dc_));
        CHECK(mvec3dc(MVec3dc_ + TriVec3dc_) == mvec3dc(MVec3dc_) + mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(MVec3dc_ - TriVec3dc_) == mvec3dc(MVec3dc_) - mvec3dc(TriVec3dc_));
        CHECK(mvec3dc(MVec3dc_ + QuadVec3dc_) ==
              mvec3dc(MVec3dc_) + mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(MVec3dc_ - QuadVec3dc_) ==
              mvec3dc(MVec3dc_) - mvec3dc(QuadVec3dc_));
        CHECK(mvec3dc(MVec3dc_ + PScalar3dc_) ==
              mvec3dc(MVec3dc_) + mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(MVec3dc_ - PScalar3dc_) ==
              mvec3dc(MVec3dc_) - mvec3dc(PScalar3dc_));
        CHECK(mvec3dc(MVec3dc_ + MVec3dc_E_) == mvec3dc(MVec3dc_) + mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(MVec3dc_ - MVec3dc_E_) == mvec3dc(MVec3dc_) - mvec3dc(MVec3dc_E_));
        CHECK(mvec3dc(MVec3dc_ + MVec3dc_U_) == mvec3dc(MVec3dc_) + mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(MVec3dc_ - MVec3dc_U_) == mvec3dc(MVec3dc_) - mvec3dc(MVec3dc_U_));
        CHECK(mvec3dc(MVec3dc_ + MVec3dc_) == mvec3dc(MVec3dc_) + mvec3dc(MVec3dc_));
        CHECK(mvec3dc(MVec3dc_ - MVec3dc_) == mvec3dc(MVec3dc_) - mvec3dc(MVec3dc_));
    }

} // TEST_SUITE("CGA 3dc Tests")
