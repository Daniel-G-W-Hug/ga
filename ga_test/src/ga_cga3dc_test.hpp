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
