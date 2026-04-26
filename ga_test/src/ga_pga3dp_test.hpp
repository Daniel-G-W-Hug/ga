// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <chrono>
#include <iostream>
#include <tuple>
#include <vector>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_pga.hpp"

#include "ga/ga_ega.hpp" // for hd::ga::ega::cmpl() and hd::ga::ega::dot()

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// projective geometric algebra 3d: pga3dp (embedded in a 4d representational space)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA 3DP Tests")
{
    TEST_CASE("G<3,0,1>: pga3dp")
    {
        fmt::println("");
        fmt::println("G<3,0,1>: pga3dp");
        // 3d projective geometric algebra (3dp) models three-dimensional Euclidean
        // algebra embedding it in a four-dimensional representational space
        const algebra<3, 0, 1> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   pga3dp: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga3dp: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3dp: default init")
    {
        fmt::println("Vec3dp: default init");
        vec3dp v;
        // fmt::println("   v = {}", v);
        CHECK(abs(v.x) < eps);
        CHECK(abs(v.y) < eps);
        CHECK(abs(v.z) < eps);
        CHECK(abs(v.w) < eps);
    }

    TEST_CASE("Vec3dp: with curly braced intializer")
    {
        fmt::println("Vec3dp: with curly braced intializer");
        vec3dp v{0.0, 0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(abs(v.x) < eps);
        CHECK(abs(v.y) < eps);
        CHECK(abs(v.z) < eps);
        CHECK(abs(v.w) < eps);
    }

    TEST_CASE("Vec3dp: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec3dp: cp ctor & cp assign incl. type deduction");
        vec3dp v1{1.0, 2.0, 3.0, 1.0}; // init with double (type deduction)
        vec3dp v2{v1};                 // cp ctor
        vec3dp v3 = v2;                // cp assign
        vec3dp v4 = -v2;               // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(abs(v1.x - 1.0) < eps);
        CHECK(abs(v1.y - 2.0) < eps);
        CHECK(abs(v1.z - 3.0) < eps);
        CHECK(abs(v1.w - 1.0) < eps);
        CHECK(abs(v2.x - 1.0) < eps);
        CHECK(abs(v2.y - 2.0) < eps);
        CHECK(abs(v2.z - 3.0) < eps);
        CHECK(abs(v2.w - 1.0) < eps);
        CHECK(abs(v3.x - 1.0) < eps);
        CHECK(abs(v3.y - 2.0) < eps);
        CHECK(abs(v3.z - 3.0) < eps);
        CHECK(abs(v3.w - 1.0) < eps);
        CHECK(v4 == -v2);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }

    TEST_CASE("Vec3dp: fmt & cout printing")
    {
        fmt::println("Vec3dp: fmt & cout printing");

        vec3dp pf{1.0f, 2.00001f, 3.0f, 1.0f};
        vec3dp pd{1.0, 2.00001, 3.0, 1.0};

        fmt::println("");
        // std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.8f}", pf);

        // std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.8f}", pd);

        std::vector<Vec3dp<double>> vp1{{1.0, 1.0, 1.0, 1.0}, {1.5, 2.0, 3.0, 1.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:e}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(geom_nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("Vec3dp: comparison float")
    {
        fmt::println("Vec3dp: comparison float");

        Vec3dp<float> v1f{1.0f, 2.0f, 3.0f, 1.0f};
        Vec3dp<float> v2f{2.0f, 4.0f, 3.0f, 1.0f};
        Vec3dp<float> v3f{1.0f, 2.0000001f, 3.0f, 1.0f};
        Vec3dp<float> v4f{v1f};

        // fmt::println("   v1f = {}", v1f);
        // fmt::println("   v2f = {}", v2f);
        // fmt::println("   v3f = {}", v3f);
        // fmt::println("   v4f = {}", v4f);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<float>::epsilon());

        CHECK(v1f == v4f);                           // comparison (equality)
        CHECK(v1f != v2f);                           // comparison (inequality)
        CHECK(bulk_nrm_sq(v1f) < bulk_nrm_sq(v2f));  // comparison (less than)
        CHECK(bulk_nrm_sq(v2f) >= bulk_nrm_sq(v1f)); // comparison (greater than or equal)
        CHECK(v3f == v1f);                           // comparison (eqality)
    }

    TEST_CASE("Vec3dp: comparison double")
    {
        fmt::println("Vec3dp: comparison double");

        Vec3dp<double> v1d{1.0, 2.0, 3.0, 1.0};
        Vec3dp<double> v2d{2.0, 4.0, 3.0, 1.0};
        Vec3dp<double> v3d{1.0, 2.0000000000000001, 3.0, 1.0};
        Vec3dp<double> v4d{v1d};

        // fmt::println("   v1d = {}", v1d);
        // fmt::println("   v2d = {}", v2d);
        // fmt::println("   v3d = {}", v3d);
        // fmt::println("   v4d = {}", v4d);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<double>::epsilon());

        CHECK(v1d == v4d);                     // comparison (equality)
        CHECK(v1d != v2d);                     // comparison (inequality)
        CHECK(bulk_nrm(v1d) < bulk_nrm(v2d));  // comparison norm
        CHECK(bulk_nrm(v2d) >= bulk_nrm(v1d)); // comparison norm
        CHECK(v3d == v1d);                     // comparison (eqality)
    }

    TEST_CASE("Vec3dp: vector space and linearity tests")
    {
        fmt::println("Vec3dp: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        vec3dp p0;
        vec3dp p1{1.0, 2.0, 3.0, 1.0};
        vec3dp p2{2.0, 4.0, 6.0, 2.0};
        vec3dp p3{3.0, 6.0, 9.0, 3.0};
        vec3dp p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t. multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t.addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
        CHECK(p4 == -p1);
    }

    TEST_CASE("Vec3dp: inner product properties")
    {
        fmt::println("Vec3dp: inner product properties");

        double a = 2.35;
        vec3dp u{1.0, 2.0, 1.0, 1.0};
        vec3dp v{-0.5, 3.0, 0.5, 1.0};
        vec3dp w{3.0, 6.0, -3.0, 1.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    TEST_CASE("Vec3dp: linking inner and outer product")
    {
        fmt::println("Vec3dp: linking inner and outer product");

        auto s1 = scalar3dp{2.0};
        auto s2 = scalar3dp{3.0};
        auto v1 = vec3dp{1.0, 2.0, 3.0, 1.0};
        auto v2 = vec3dp{-1.0, 4.0, -5.0, 1.0};
        auto B1 = bivec3dp{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        auto B2 = bivec3dp{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0};
        auto t1 = trivec3dp{10.0, 20.0, 30.0, 10.0};
        auto t2 = trivec3dp{-10.0, 40.0, -50.0, 10.0};
        auto ps1 = pscalar3dp{-2.0};
        auto ps2 = pscalar3dp{4.0};

        CHECK(r_bulk_dual(wdg(v1, B1)) == rwdg(r_bulk_dual(v1), r_bulk_dual(B1)));

        fmt::println("");
        fmt::println("wdg(s1, r_bulk_dual(s2)) = {}", wdg(s1, r_bulk_dual(s2)));
        fmt::println("dot(s1, s2) ^ I_3dp = {}", wdg(dot(s1, s2), I_3dp));
        fmt::println("");
        fmt::println("wdg(v1, r_bulk_dual(v2)) = {}", wdg(v1, r_bulk_dual(v2)));
        fmt::println("dot(v1, v2) ^ I_3dp = {}", wdg(dot(v1, v2), I_3dp));
        fmt::println("");
        fmt::println("wdg(B1, r_bulk_dual(B2)) = {}", wdg(B1, r_bulk_dual(B2)));
        fmt::println("dot(B1, B2) ^ I_3dp = {}", wdg(dot(B1, B2), I_3dp));
        fmt::println("");
        fmt::println("wdg(t1, r_bulk_dual(t2)) = {}", wdg(t1, r_bulk_dual(t2)));
        fmt::println("dot(t1, t2) ^ I_3dp = {}", wdg(dot(t1, t2), I_3dp));
        fmt::println("");
        fmt::println("wdg(ps1, r_bulk_dual(ps2)) = {}", wdg(ps1, r_bulk_dual(ps2)));
        fmt::println("dot(ps1, ps2) ^ I_3dp = {}", wdg(dot(ps1, ps2), I_3dp));
        fmt::println("");

        CHECK(wdg(s1, r_bulk_dual(s2)) == wdg(dot(s1, s2), I_3dp));
        //
        CHECK(wdg(v1, r_bulk_dual(v2)) == wdg(dot(v1, v2), I_3dp));
        //
        CHECK(wdg(B1, r_bulk_dual(B2)) == wdg(dot(B1, B2), I_3dp));
        //
        CHECK(wdg(t1, r_bulk_dual(t2)) == wdg(dot(t1, t2), I_3dp));
        //
        CHECK(wdg(ps1, r_bulk_dual(ps2)) == wdg(dot(ps1, ps2), I_3dp));

        CHECK(wdg(v1, r_weight_dual(v2)) == rdot(v1, v2));

        fmt::println("");

        // zero values can always be compared (only true if the value is actually zero!)
        CHECK(to_val(r_weight_dual(rwdg(v1, B1))) ==
              to_val(wdg(r_weight_dual(v1), r_weight_dual(B1))));
        CHECK(to_val(r_weight_dual(rwdg(v1, B1))) == 0.0);
        CHECK(to_val(wdg(r_weight_dual(v1), r_weight_dual(B1))) == 0.0);
        // alternative valid way to test this in one step is to test this via conversion
        // to full multivector:
        CHECK(mvec3dp(r_weight_dual(rwdg(v1, B1))) ==
              mvec3dp(wdg(r_weight_dual(v1), r_weight_dual(B1))));

        CHECK(rwdg(v1, r_bulk_dual(v2)) == dot(v1, v2));

        fmt::println("");
        fmt::println("rwdg(s1, r_weight_dual(s2)) = {}", rwdg(s1, r_weight_dual(s2)));
        fmt::println("rwdg(rdot(s1, s2), I_3dp) = {}",
                     rwdg(rdot(s1, s2), scalar3dp{1.0}));
        fmt::println("");
        fmt::println("rwdg(v1, r_weight_dual(v2)) = {}", rwdg(v1, r_weight_dual(v2)));
        fmt::println("rwdg(rdot(v1, v2), I_3dp) = {}",
                     rwdg(rdot(v1, v2), scalar3dp{1.0}));
        fmt::println("");
        fmt::println("rwdg(B1, r_weight_dual(B2)) = {}", rwdg(B1, r_weight_dual(B2)));
        fmt::println("rwdg(rdot(B1, B2), I_3dp) = {}",
                     rwdg(rdot(B1, B2), scalar3dp{1.0}));
        fmt::println("");
        fmt::println("rwdg(t1, r_weight_dual(t2)) = {}", rwdg(t1, r_weight_dual(t2)));
        fmt::println("rwdg(rdot(t1, t2), I_3dp) = {}",
                     rwdg(rdot(t1, t2), scalar3dp{1.0}));
        fmt::println("");
        fmt::println("rwdg(ps1, r_weight_dual(ps2)) = {}", rwdg(ps1, r_weight_dual(ps2)));
        fmt::println("rwdg(rdot(ps1, ps2), I_3dp) = {}",
                     rwdg(rdot(ps1, ps2), scalar3dp{1.0}));
        fmt::println("");

        CHECK(rwdg(s1, r_weight_dual(s2)) == rwdg(rdot(s1, s2), scalar3dp{1.0}));
        //
        CHECK(rwdg(v1, r_weight_dual(v2)) == rwdg(rdot(v1, v2), scalar3dp{1.0}));
        //
        CHECK(rwdg(B1, r_weight_dual(B2)) == rwdg(rdot(B1, B2), scalar3dp{1.0}));
        //
        CHECK(rwdg(t1, r_weight_dual(t2)) == rwdg(rdot(t1, t2), scalar3dp{1.0}));
        //
        CHECK(rwdg(ps1, r_weight_dual(ps2)) == rwdg(rdot(ps1, ps2), scalar3dp{1.0}));

        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3dp<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3dp: operations - norm, inverse, dot")
    {
        fmt::println("Vec3dp: operations - norm, inverse, dot");

        scalar3dp s1{3.2};

        vec3dp v1{2.0, 1.0, 2.0, -2.0};
        vec3dp v1n{bulk_normalize(v1)};
        vec3dp v1u{unitize(v1)};

        vec3dp v3{2.0, 6.0, -4.0, 2.0};
        vec3dp v4{inv(v3)};

        bivec3dp b1{-2.3, 1.2, 4.5, 4.1, -2.6, 5.2};
        trivec3dp t1{-2.0, 6.2, 4.0, -12.0};
        pscalar3dp ps1{-4.7};
        mvec3dp_e mve1{s1, b1, ps1};
        mvec3dp_u mvu1{v1, t1};
        mvec3dp mv1{s1, v1, b1, t1, ps1};

        // auto v3m = mvec3dp{v3};
        // auto v4m = mvec3dp{v4};
        // fmt::println("v1                  = {: .4f}, nrm(v1)         = {: .4f}",
        // v1,
        //              nrm(v1));
        // fmt::println("v1                  = {: .4f}, bulk_nrm(v1)    = {: .4f}",
        // v1,
        //              bulk_nrm(v1));
        // fmt::println("v1                  = {: .4f}, weight_nrm(v1)  = {: .4f}",
        // v1,
        //              weight_nrm(v1));
        // fmt::println("v1n = normalize(v1) = {: .4f}, nrm(v1n)        = {: .4f}",
        // v1n,
        //              nrm(v1n));
        // fmt::println("v1n                 = {: .4f}, bulk_nrm(v1n)   = {: .4f}",
        // v1n,
        //              bulk_nrm(v1n));
        // fmt::println("v1n                 = {: .4f}, weight_nrm(v1n) = {: .4f}",
        // v1n,
        //              weight_nrm(v1n));
        // fmt::println("v1u = unitize(v1)   = {: .4f}, nrm(v1u)        = {: .4f}",
        // v1u,
        //              nrm(v1u));
        // fmt::println("v1u                 = {: .4f}, bulk_nrm(v1u)   = {: .4f}",
        // v1u,
        //              bulk_nrm(v1u));
        // fmt::println("v1u                 = {: .4f}, weight_nrm(v1u) = {: .4f}",
        // v1u,
        //              weight_nrm(v1u));
        // fmt::println("");
        // fmt::println("v3                  = {: .4f}, nrm(v1) = {: .4f}", v3,
        // nrm(v3)); fmt::println("v4 = inv(v3)        = {: .4f}, "
        //              "nrm(v3) = {: .4f}, nrm(v3)*nrm(v4) = {: .4f}",
        //              v4, nrm(v4), nrm(v3) * nrm(v4));
        // fmt::println("");
        // fmt::println("v3m = {}, nrm(v3m) = {}", v3m, nrm(v3m));
        // fmt::println("v4m = {}, nrm(v4m) = {}", v4m, nrm(v4m));
        // fmt::println("0.5*(v3m*v4m + v4m*v3m) = {}", 0.5 * (v3m * v4m + v4m *
        // v3m)); fmt::println("0.5*(v3m*v4m - v4m*v3m) = {}", 0.5 * (v3m * v4m - v4m
        // * v3m));

        CHECK(abs(bulk_nrm_sq(v1) + weight_nrm_sq(v1) - 13) < eps);
        CHECK(abs(bulk_nrm_sq(v1) - 9) < eps);
        CHECK(abs(weight_nrm_sq(v1) - 4) < eps);

        CHECK(abs(bulk_nrm_sq(v1n) - 1) < eps); // v1n = bulk_normalize(v1)
        CHECK(abs(weight_nrm_sq(v1u) - 1) < eps);

        CHECK(abs(bulk_nrm_sq(v3) - 56) < eps);
        CHECK(abs(to_val(bulk_nrm(v4 * v3)) - 1.0) < eps);
        CHECK(abs(dot(v4, v3) - 1.0) < eps);
        CHECK(abs(to_val(bulk_nrm(wdg(v4, v3)))) < eps);

        // check inverses - scalar
        // fmt::println("");
        // fmt::println("s1 * inv(s1) = {}", s1 * inv(s1)); // s
        CHECK(abs(to_val(bulk_nrm(s1 * inv(s1))) - 1) < eps);
        CHECK(abs(inv(s1) - rev(s1) / bulk_nrm_sq(s1)) < eps);

        // check inverses - vector
        // fmt::println("v1 * inv(v1) = {}", v1 * inv(v1)); // mv_e
        CHECK(abs(to_val(bulk_nrm(gr0(v1 * inv(v1)))) - 1) < eps);
        CHECK(abs(to_val(bulk_nrm(gr2(v1 * inv(v1)))) - 0) < eps);
        CHECK(abs(to_val(bulk_nrm(gr4(v1 * inv(v1))))) < eps);
        CHECK(abs(to_val(bulk_nrm(inv(v1) - rev(v1) / bulk_nrm_sq(v1)))) < eps);

        // check inverses - bivector
        // fmt::println("b1 * inv(b1) = {}", b1 * inv(b1)); // mv_e
        CHECK(abs(to_val(bulk_nrm(gr0(b1 * inv(b1)))) - 1) < eps);
        CHECK(abs(to_val(bulk_nrm(gr2(b1 * inv(b1)))) - 0) < eps);
        CHECK(abs(to_val(bulk_nrm(gr4(b1 * inv(b1))))) < eps);

        // check inverses - trivector
        // fmt::println("t1 * inv(t1) = {}", t1 * inv(t1)); // mv_e
        CHECK(abs(to_val(bulk_nrm(gr0(t1 * inv(t1)))) - 1) < eps);
        CHECK(abs(to_val(bulk_nrm(gr2(t1 * inv(t1)))) - 0) < eps);
        CHECK(abs(to_val(bulk_nrm(gr4(t1 * inv(t1))))) < eps);


        // check inverses - pseudoscalar
        // due to the degenerate metric there is no inverse of the pseudoscalar

        // check inverses - even grade multivector
        // fmt::println("mve1 * inv(mve1) = {}", mve1 * inv(mve1)); // mv_e
        CHECK(abs(to_val(bulk_nrm(gr0(mve1 * inv(mve1)))) - 1) < eps);
        CHECK(abs(to_val(bulk_nrm(gr2(mve1 * inv(mve1)))) - 0) < eps);
        CHECK(abs(to_val(bulk_nrm(gr4(mve1 * inv(mve1))))) < eps);

        // check inverses - odd grade multivector
        // fmt::println("mvu1 * inv(mvu1) = {}", mvu1 * inv(mvu1)); // mv_e
        CHECK(abs(to_val(bulk_nrm(gr0(mvu1 * inv(mvu1)))) - 1) < eps);
        CHECK(abs(to_val(bulk_nrm(gr2(mvu1 * inv(mvu1)))) - 0) < eps);
        CHECK(abs(to_val(bulk_nrm(gr4(mvu1 * inv(mvu1))))) < eps);

        // check inverses - multivector
        // fmt::println("mv1 * inv(mv1) = {}", mv1 * inv(mv1)); // mv
        CHECK(abs(to_val(bulk_nrm(gr0(mv1 * inv(mv1)))) - 1) < eps);
        CHECK(abs(to_val(bulk_nrm(gr1(mv1 * inv(mv1))))) < eps);
        CHECK(abs(to_val(bulk_nrm(gr2(mv1 * inv(mv1))))) < eps);
        CHECK(abs(to_val(bulk_nrm(gr3(mv1 * inv(mv1))))) < eps);
        CHECK(abs(to_val(bulk_nrm(gr4(mv1 * inv(mv1))))) < eps);
        CHECK(abs(to_val(bulk_nrm(gr0(inv(mv1) * mv1))) - 1) < eps); // left and right
                                                                     // inverse are equal
        // fmt::println("");
    }

    TEST_CASE("Vec3dp: operations - bulk_dual, weight_dual")
    {
        fmt::println("Vec3dp: operations - bulk_dual, weight_dual");

        auto s = scalar3dp{5.0};
        auto v = vec3dp{1.0, 2.0, 1.0, 1.0};
        auto B = bivec3dp{-1.0, 2.0, 1.0, -10.0, 20.0, 10.0};
        auto t = trivec3dp{3.0, 6.0, 3.0, 3.0};
        auto ps = pscalar3dp{-5.0};

        auto v2 = vec3dp{3.0, -2.0, 2.0, 5.0};
        auto B2 = bivec3dp{1.0, -2.0, 3.0, -10.0, -20.0, 30.0};

        CHECK(r_bulk_dual(s) == pscalar3dp(to_val(s)));
        CHECK(r_weight_dual(s) == pscalar3dp{0.0});

        CHECK(r_bulk_dual(v) == trivec3dp{1.0, 2.0, 1.0, 0.0});
        CHECK(r_weight_dual(v) == trivec3dp{0.0, 0.0, 0.0, 1.0});

        CHECK(r_bulk_dual(B) == bivec3dp{10.0, -20.0, -10.0, 0.0, 0.0, 0.0});
        CHECK(r_weight_dual(B) == bivec3dp{0.0, 0.0, 0.0, 1.0, -2.0, -1.0});

        CHECK(r_bulk_dual(t) == vec3dp{0.0, 0.0, 0.0, -3.0});
        CHECK(r_weight_dual(t) == vec3dp{-3.0, -6.0, -3.0, 0.0});

        CHECK(r_bulk_dual(ps) == scalar3dp{0.0});
        CHECK(r_weight_dual(ps) == scalar3dp(to_val(ps)));

        // duality of wdg and rwdg based on complements
        CHECK(r_bulk_dual(wdg(v, v2)) == rwdg(r_bulk_dual(v), r_bulk_dual(v2)));
        CHECK(r_bulk_dual(wdg(v, B)) == rwdg(r_bulk_dual(v), r_bulk_dual(B)));
        CHECK(r_bulk_dual(wdg(v, t)) == rwdg(r_bulk_dual(v), r_bulk_dual(t)));
        CHECK(r_bulk_dual(wdg(t, v)) == rwdg(r_bulk_dual(t), r_bulk_dual(v)));
        CHECK(r_bulk_dual(wdg(B, v)) == rwdg(r_bulk_dual(B), r_bulk_dual(v)));
        CHECK(r_bulk_dual(wdg(v, B)) == rwdg(r_bulk_dual(v), r_bulk_dual(B)));
        CHECK(r_bulk_dual(wdg(B, B2)) == rwdg(r_bulk_dual(B), r_bulk_dual(B2)));

        // contractions
        CHECK(l_bulk_contract3dp(v, v2) == rwdg(l_bulk_dual(v), v2));
        CHECK(l_bulk_contract3dp(B, B2) == rwdg(l_bulk_dual(B), B2));
        CHECK(l_bulk_contract3dp(B, v) == rwdg(l_bulk_dual(B), v));

        CHECK(l_weight_contract3dp(v, v2) == rwdg(l_weight_dual(v), v2));
        CHECK(l_weight_contract3dp(B, B2) == rwdg(l_weight_dual(B), B2));
        CHECK(l_weight_contract3dp(B, v) == rwdg(l_weight_dual(B), v));

        CHECK(r_bulk_contract3dp(v, v2) == rwdg(v, r_bulk_dual(v2)));
        CHECK(r_bulk_contract3dp(B, B2) == rwdg(B, r_bulk_dual(B2)));
        CHECK(r_bulk_contract3dp(B, v) == rwdg(B, r_bulk_dual(v)));

        CHECK(r_weight_contract3dp(v, v2) == rwdg(v, r_weight_dual(v2)));
        CHECK(r_weight_contract3dp(B, B2) == rwdg(B, r_weight_dual(B2)));
        CHECK(r_weight_contract3dp(B, v) == rwdg(B, r_weight_dual(v)));

        // expansions
        CHECK(l_weight_expand3dp(v, v2) == wdg(l_weight_dual(v), v2));
        CHECK(l_weight_expand3dp(B, B2) == wdg(l_weight_dual(B), B2));
        CHECK(l_weight_expand3dp(B, v) == wdg(l_weight_dual(B), v));

        CHECK(l_bulk_expand3dp(v, v2) == wdg(l_bulk_dual(v), v2));
        CHECK(l_bulk_expand3dp(B, B2) == wdg(l_bulk_dual(B), B2));
        CHECK(l_bulk_expand3dp(B, v) == wdg(l_bulk_dual(B), v));

        CHECK(r_weight_expand3dp(v, v2) == wdg(v, r_weight_dual(v2)));
        CHECK(r_weight_expand3dp(B, B2) == wdg(B, r_weight_dual(B2)));
        CHECK(r_weight_expand3dp(v, B) == wdg(v, r_weight_dual(B)));

        CHECK(r_bulk_expand3dp(v, v2) == wdg(v, r_bulk_dual(v2)));
        CHECK(r_bulk_expand3dp(B, B2) == wdg(B, r_bulk_dual(B2)));
        CHECK(r_bulk_expand3dp(v, B) == wdg(v, r_bulk_dual(B)));
    }

    TEST_CASE("Vec3dp: operations - angle I")
    {
        fmt::println("Vec3dp: operations - angle I");

        vec3dp v1{1.0, 0.0, 0.0, 0.0};
        vec3dp v2{bulk_normalize(vec3dp(1.0, 1.0, 0.0, 0.0))};
        vec3dp v3{0.0, 1.0, 0.0, 0.0};
        vec3dp v4{bulk_normalize(vec3dp(-1.0, 1.0, 0.0, 0.0))};
        vec3dp v5{-1.0, 0.0, 0.0, 0.0};
        vec3dp v6{bulk_normalize(vec3dp(-1.0, -1.0, 0.0, 0.0))};
        vec3dp v7{0.0, -1.0, 0.0, 0.0};
        vec3dp v8{bulk_normalize(vec3dp(1.0, -1.0, 0.0, 0.0))};

        // fmt::println("v1 = {: .4f}, nrm(v1) = {:.8f}, "
        //              "angle(v1,v1) = {:.8f}, {:.8f}",
        //              v1, nrm(v1), angle(v1, v1), angle(v1, v1) / pi);
        // fmt::println("v2 = {: .4f}, nrm(v2) = {:.8f}, "
        //              "angle(v1,v2) = {:.8f}, {:.8f}",
        //              v2, nrm(v2), angle(v1, v2), angle(v1, v2) / pi);
        // fmt::println("v3 = {: .4f}, nrm(v3) = {:.8f}, "
        //              "angle(v1,v3) = {:.8f}, {:.8f} ",
        //              v3, nrm(v3), angle(v1, v3), angle(v1, v3) / pi);
        // fmt::println("v4 = {: .4f}, nrm(v4) = {:.8f}, "
        //              "angle(v1,v4) = {:.8f}, {:.8f} ",
        //              v4, nrm(v4), angle(v1, v4), angle(v1, v4) / pi);
        // fmt::println("v5 = {: .4f}, nrm(v5) = {:.8f}, "
        //              "angle(v1,v5) = {:.8f}, {:.8f} ",
        //              v5, nrm(v5), angle(v1, v5), angle(v1, v5) / pi);
        // fmt::println("v6 = {: .4f}, nrm(v6) = {:.8f}, "
        //              "angle(v1,v6) = {:.8f}, {:.8f} ",
        //              v6, nrm(v6), angle(v1, v6), angle(v1, v6) / pi);
        // fmt::println("v7 = {: .4f}, nrm(v7) = {:.8f}, "
        //              "angle(v1,v7) = {:.8f}, {:.8f} ",
        //              v7, nrm(v7), angle(v1, v7), angle(v1, v7) / pi);
        // fmt::println("v8 = {: .4f}, nrm(v8) = {:.8f}, "
        //              "angle(v1,v8) = {:.8f}, {:.8f} ",
        //              v8, nrm(v8), angle(v1, v8), angle(v1, v8) / pi);

        CHECK(abs(angle(v1, v1) - 0.0) < eps);
        CHECK(abs(angle(v1, v2) - pi * 0.25) < eps);
        CHECK(abs(angle(v1, v3) - pi * 0.5) < eps);
        CHECK(abs(angle(v1, v4) - pi * 0.75) < eps);
        CHECK(abs(angle(v1, v5) - pi) < eps);

        // just to suppress unused variable warnings
        CHECK(v6 == bulk_normalize(vec3dp(-1.0, -1.0, 0.0, 0.0)));
        CHECK(v7 == bulk_normalize(vec3dp(0.0, -1.0, 0.0, 0.0)));
        CHECK(v8 == bulk_normalize(vec3dp(1.0, -1.0, 0.0, 0.0)));
    }

    TEST_CASE("Vec3dp: operations - angle II")
    {
        fmt::println("Vec3dp: operations - angle II");

        std::vector<std::tuple<double, Vec3dp<double>>> v1;
        std::vector<std::tuple<double, Vec3dp<double>>> v2;
        std::vector<std::tuple<double, Vec3dp<double>>> v3;

        // only positive angles are easy to implement vs. the 2d case

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3dp<double>(std::cos(phi), std::sin(phi), 0.0, 0.0);
            v1.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_3dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c =
                Vec3dp<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2), 0.0, 0.0);
            v2.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e2_3dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c =
                Vec3dp<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4), 0.0, 0.0);
            v3.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_3dp + e2_3dp, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(abs(phi - angle(e1_3dp, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(abs(phi - angle(e2_3dp, c)) < eps);
        }
        auto ref_vec = bulk_normalize(e1_3dp + e2_3dp);
        for (auto const& [phi, c] : v3) {
            CHECK(abs(phi - angle(ref_vec, c)) < eps);
        }
    }


    TEST_CASE("Vec3dp: operations - wedge")
    {
        fmt::println("Vec3dp: operations - wedge");

        vec3dp v1{1.0, 0.0, 0.0, 0.0};
        vec3dp v2{bulk_normalize(vec3dp(1.0, 1.0, 0.0, 0.0))};
        vec3dp v3{0.0, 1.0, 0.0, 0.0};
        vec3dp v4{bulk_normalize(vec3dp(-1.0, 1.0, 0.0, 0.0))};
        vec3dp v5{-1.0, 0.0, 0.0, 0.0};
        vec3dp v6{bulk_normalize(vec3dp(-1.0, -1.0, 0.0, 0.0))};
        vec3dp v7{0.0, -1.0, 0.0, 0.0};
        vec3dp v8{bulk_normalize(vec3dp(1.0, -1.0, 0.0, 0.0))};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3dp{sd};
        auto t = scalar3dp{st};

        // fmt::println("v1 = {: .4f}, wdg(v1,v1) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v1, wdg(v1, v1), angle(v1, v1));
        // fmt::println("v2 = {: .4f}, wdg(v1,v2) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v2, wdg(v1, v2), angle(v1, v2));
        // fmt::println("v3 = {: .4f}, wdg(v1,v3) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v3, wdg(v1, v3), angle(v1, v3));
        // fmt::println("v4 = {: .4f}, wdg(v1,v4) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v4, wdg(v1, v4), angle(v1, v4));
        // fmt::println("v5 = {: .4f}, wdg(v1,v5) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v5, wdg(v1, v5), angle(v1, v5));
        // fmt::println("v6 = {: .4f}, wdg(v1,v6) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v6, wdg(v1, v6), angle(v1, v6));
        // fmt::println("v7 = {: .4f}, wdg(v1,v7) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v7, wdg(v1, v7), angle(v1, v7));
        // fmt::println("v8 = {: .4f}, wdg(v1,v8) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v8, wdg(v1, v8), angle(v1, v8));

        CHECK(wdg(v1, v1) == bivec3dp{});                     // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s));   // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, v1) == wdg(v1, s)); // wdg between scalar and vector
        CHECK(wdg(s, v1) == sd * v1);    // wdg between scalar and vector

        CHECK(abs(bulk_nrm(wdg(v1, v1)) - sin(angle(v1, v1))) < eps);
        CHECK(abs(bulk_nrm(wdg(v1, v2)) - sin(angle(v1, v2))) < eps);
        CHECK(abs(bulk_nrm(wdg(v1, v3)) - sin(angle(v1, v3))) < eps);
        CHECK(abs(bulk_nrm(wdg(v1, v4)) - sin(angle(v1, v4))) < eps);
        CHECK(abs(bulk_nrm(wdg(v1, v5)) - sin(angle(v1, v5))) < eps);
        CHECK(abs(bulk_nrm(wdg(v1, v6)) - sin(angle(v1, v6))) < eps);
        CHECK(abs(bulk_nrm(wdg(v1, v7)) - sin(angle(v1, v7))) < eps);
        CHECK(abs(bulk_nrm(wdg(v1, v8)) - sin(angle(v1, v8))) < eps);
    }

    TEST_CASE("Vec3dp: operations - project / reject / reflect")
    {
        fmt::println("Vec3dp: operations - project / reject / reflect");

        vec3dp v1{5, 1, 1, 0};
        vec3dp v2{2, 2, 1, 0};


        // vectors (=directions) projected and rejected
        vec3dp v3{project_onto(v1, v2)};
        vec3dp v4{reject_from(v1, v2)};
        CHECK(v3 + v4 == v1);
        CHECK(std::abs(to_val(dot(v3, v4))) <= eps);

        fmt::println("v1  = {: .4f}, bulk_nrm(v1) = {: .4f}", v1, bulk_nrm(v1));
        fmt::println("v2  = {: .4f}, bulk_nrm(v2) = {: .4f}", v2, bulk_nrm(v2));
        fmt::println("");
        fmt::println("v3 = project_onto(v1, v2) = {: .4f}, bulk_nrm(v3) = {: .4f}", v3,
                     bulk_nrm(v3));
        fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, bulk_nrm(v4) = {: .4f}", v4,
                     bulk_nrm(v4));
        fmt::println("v3 + v4              = {: .4f}, bulk_nrm(v3 + v4) = {: .4f}",
                     v3 + v4, bulk_nrm(v3 + v4));
        fmt::println("");

        // points
        vec3dp p1{5, 1, 1, 1};
        vec3dp p2{2, 2, 1, 1};

        // points projected onto and rejected from line
        bivec3dp B{e41_3dp};

        vec3dp p3{project_onto(p1, B)};
        vec3dp p3v{project_onto(v1, B)};
        vec3dp v5{reject_from(p1, B)};
        CHECK(p1 == p3 + v5);
        CHECK(dot(p3, v5) == 0);

        fmt::println("p3 = project_onto(p1,B) = {: .4f}, bulk_nrm(p3)  = {: .4f}", p3,
                     bulk_nrm(p3));
        fmt::println("v5 = reject_from(p1,B)  = {: .4f}, bulk_nrm(p3)  = {: .4f}", v5,
                     bulk_nrm(v5));
        fmt::println("p3v  = {: .4f}, bulk_nrm(p3v) = {: .4f}", p3v, bulk_nrm(p3v));
        fmt::println("");

        // project a line on a plane
        auto l = unitize(wdg(vec3dp{1, 1, 1, 1}, vec3dp{1, -2, 3, 1}));

        // plane
        trivec3dp t{e431_3dp};

        auto lp = project_onto(l, t);
        fmt::println("l  = {: .4f}", l);
        fmt::println("t  = {: .4f}", t);
        fmt::println("lp = {: .4f}", lp);
        fmt::println("");


        // point to reflect
        auto v = vec3dp{1, 3, 1, 1};
        auto vr = vec3dp{1, -3, 1, 1};

        // fmt::println("");
        // fmt::println("v = {}", v);
        // fmt::println("reflect_on(v, t) = {}", reflect_on(v, t));
        // fmt::println("unitize(reflect_on(v, t)) = {}", unitize(reflect_on(v, t)));
        // fmt::println("vr = {}", vr);
        // fmt::println("");

        CHECK(reflect_on(v, t) == vr);

        // line to reflect
        p1 = vec3dp{1, 1, 0, 1};
        p2 = vec3dp{3, 2, 0, 1};
        auto p1r = vec3dp{1, -1, 0, 1};
        auto p2r = vec3dp{3, -2, 0, 1};
        auto l12 = wdg(p1, p2);
        auto l12r = wdg(p1r, p2r);
        auto lu12 = unitize(wdg(p1, p2));
        auto lu12r = unitize(wdg(p1r, p2r));

        // fmt::println("");
        // fmt::println("l12 = {}, att(l12) = {}", l12, att(l12));
        // fmt::println("l12r = {}, att(l12r) = {}", l12r, att(l12r));
        // fmt::println("lu12 = {}, att(lu12) = {}", lu12, att(lu12));
        // fmt::println("lu12r = {}, att(lu12r) = {}", lu12r, att(lu12r));
        // fmt::println("");

        CHECK(reflect_on(l12, t) == l12r);
        CHECK(reflect_on(lu12, t) == lu12r);

        // plane to reflect
        p3 = vec3dp{1, 1, 1, 1};
        auto p3r = vec3dp{1, -1, 1, 1};
        auto p123 = unitize(wdg(-l12, p3));
        auto p123r = unitize(wdg(l12r, p3r));

        // fmt::println("");
        // fmt::println("p123 = {}, att(p123) = {}", p123, att(p123));
        // fmt::println("p123r = {}, att(p123r) = {}", p123r, att(p123r));
        // fmt::println("");

        CHECK(reflect_on(p123, t) == p123r);


        // point reflected on a plane
        vec3dp p{4, 1, 1, 1};
        CHECK(reflect_on(p, zx_3dp) == vec3dp{4, -1, 1, 1});

        // coordinate axis reflected on perpendicular base planes yield their negatives
        CHECK(reflect_on(x_axis_3dp, yz_3dp) == -x_axis_3dp);
        CHECK(reflect_on(y_axis_3dp, zx_3dp) == -y_axis_3dp);
        CHECK(reflect_on(z_axis_3dp, xy_3dp) == -z_axis_3dp);

        // lines parallel to coordinate axis after reflexion:
        // remain parallel, have same orientation, but are on other side of axis
        p1 = vec3dp{0, 1, 1, 1};
        p2 = vec3dp{1, 1, 1, 1};
        l12 = unitize(wdg(p1, p2)); // l12 = BiVec3dp(1, 0, 0, 0, 1, -1)
        p1r = reflect_on(p1, zx_3dp);
        p2r = reflect_on(p2, zx_3dp);
        l12r = unitize(wdg(p1r, p2r)); // l12r = BiVec3dp(1, 0, 0, 0, 1, 1)
        auto basept = unitize(rwdg(l12, yz_3dp));
        auto baseptr = unitize(rwdg(l12r, yz_3dp));

        // fmt::println("");
        // fmt::println("p1  = {},  p2 = {}", p1, p2);
        // fmt::println("p1r = {}, p2r = {}", p1r, p2r);
        // fmt::println("");
        // fmt::println("l12 = {}, att(l12) = {}", l12, att(l12));
        // fmt::println("l12r = {}, att(l12r) = {}", l12r, att(l12r));
        // fmt::println("");
        // fmt::println("basept  = {}", basept);
        // fmt::println("baseptr = {}", baseptr);
        // fmt::println("");

        CHECK(p1r == vec3dp{0, -1, 1, 1});
        CHECK(p2r == vec3dp{1, -1, 1, 1});
        CHECK(reflect_on(l12, zx_3dp) == bivec3dp{1, 0, 0, 0, 1, 1});
        CHECK(reflect_on(l12, xy_3dp) == bivec3dp{1, 0, 0, 0, -1, -1});
        CHECK(basept == vec3dp{0, 1, 1, 1});
        CHECK(baseptr == vec3dp{0, -1, 1, 1});

        // a coordinate plane reflected on itself changes orientation, i.e. direction of
        // its normal vector!
        CHECK(reflect_on(yz_3dp, yz_3dp) == -yz_3dp);
        CHECK(reflect_on(zx_3dp, zx_3dp) == -zx_3dp);
        CHECK(reflect_on(xy_3dp, xy_3dp) == -xy_3dp);

        // a plane with a given distance to the origin, reflected on a plane parallel to
        // itself containing the origin, changes orientation and has the same distance to
        // the origin after the reflection
        p3 = {1, 1, 0, 1};
        p123 = unitize(wdg(wdg(p1, p2), p3)); // TriVec3dp(0, 1, 0, -1)
        p123r = reflect_on(p123, zx_3dp);     // TriVec3dp(-0, -1, -0, -1)

        auto s123 = sup(p123);   // point in p123 that is closest to the origin
        auto s123r = sup(p123r); // point in p123r that is closest to the origin

        // fmt::println("");
        // fmt::println("p123  = {}", p123);
        // fmt::println("p123r = {}", p123r);
        // fmt::println("att(p123)  = {}", att(p123));
        // fmt::println("att(p123r) = {}", att(p123r));
        // fmt::println("sup(p123)  = {}", s123);
        // fmt::println("sup(p123r) = {}", s123r);
        // fmt::println("");

        CHECK(p123.x == -p123r.x); // minus sign is the mirrored normal
        CHECK(p123.y == -p123r.y);
        CHECK(p123.z == -p123r.z);
        CHECK(p123.w == p123r.w);  // distance to origin is unchanged
        CHECK(s123.y == -s123r.y); // just the sign has changed
        CHECK(att(p123) == -att(p123r));

        // reflect planes on planes directly
        CHECK(reflect_on(e423_3dp + e412_3dp, e412_3dp) == e423_3dp - e412_3dp);
    }

    TEST_CASE("Vec3dp: operations - reflect in plane (=trivec)")
    {
        fmt::println("Vec3dp: operations - reflect in plane (=trivec)");

        CHECK(reflect_on(x_dir_3dp, yz_3dp) == -x_dir_3dp);
        CHECK(reflect_on(y_dir_3dp, zx_3dp) == -y_dir_3dp);
        CHECK(reflect_on(z_dir_3dp, xy_3dp) == -z_dir_3dp);

        CHECK(reflect_on(vec3dp{1, 0, 0, 1}, yz_3dp) == vec3dp{-1, 0, 0, 1});
        CHECK(reflect_on(vec3dp{0, 1, 0, 1}, zx_3dp) == vec3dp{0, -1, 0, 1});
        CHECK(reflect_on(vec3dp{0, 0, 1, 1}, xy_3dp) == vec3dp{0, 0, -1, 1});

        CHECK(reflect_on(x_axis_3dp, yz_3dp) == -x_axis_3dp);
        CHECK(reflect_on(y_axis_3dp, zx_3dp) == -y_axis_3dp);
        CHECK(reflect_on(z_axis_3dp, xy_3dp) == -z_axis_3dp);

        CHECK(reflect_on(trivec3dp{1, 0, 0, 1}, yz_3dp) == trivec3dp{-1, 0, 0, 1});
        CHECK(reflect_on(trivec3dp{0, 1, 0, 1}, zx_3dp) == trivec3dp{0, -1, 0, 1});
        CHECK(reflect_on(trivec3dp{0, 0, 1, 1}, xy_3dp) == trivec3dp{0, 0, -1, 1});
    }

    TEST_CASE("Vec3dp: operations - rotations")
    {
        fmt::println("Vec3dp: operations - rotations");

        {
            std::vector<std::tuple<double, Vec3dp<double>>> v;

            // fmt::println("");
            for (int i = -12; i <= 12; ++i) {
                double phi = i * pi / 12;
                auto c = vec3dp(std::cos(phi), std::sin(phi), 0.0, 0.0);
                auto d = move3dp(x_dir_3dp, get_motor(z_axis_3dp, phi));
                v.emplace_back(std::make_tuple(phi, c));
                // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
                //              " angle={: .4f}",
                //              i, phi, rad2deg(phi), c, rad2deg(angle(e1_3dp, c)));
                // fmt::println("                                  d={: .3f}", d);
                CHECK(c == d);
            }
        }

        {

            // rotations in a plane (using motors)
            std::vector<std::tuple<double, Vec3dp<double>>> v;

            // fmt::println("");
            for (int i = -12; i <= 12; ++i) {
                double phi = i * pi / 12;
                auto c = vec3dp(std::cos(phi), std::sin(phi), 0.0, 1.0);
                // auto b = bivec3dp(std::cos(phi), std::sin(phi), 0.0, 0.0, 0.0, 0.0);
                // auto t = trivec3dp(std::cos(phi), std::sin(phi), 0.0, 0.0);
                v.emplace_back(std::make_tuple(phi, c));
                // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
                //              " angle={: .4f}, b={: .3f}",
                //              i, phi, rad2deg(phi), c, angle(e1_3dp, c), b);
                CHECK(c == move3dp(e1_3dp + O_3dp, get_motor(e43_3dp, phi)));
                // CHECK(b == move3dp(b, get_motor(e43_3dp, phi)));
                // CHECK(t == move3dp(t, get_motor(e43_3dp, phi)));

                // check rotate optimization (use points on the projective plan)
                CHECK(move3dp(e1_3dp + O_3dp, get_motor(e43_3dp, phi)) ==
                      move3dp_opt(e1_3dp + O_3dp, get_motor(e43_3dp, phi)));

                // rotation of a ega3dp multivector is to rotate the vector and
                // trivector parts; rotation plane bivector,scalar and pscalar should be
                // unchanged

                // auto mv_rev = mvec3dp{scalar3dp{4.0}, e1_3dp + O_3dp, e41_3dp,
                // e423_3dp,
                //                       pscalar3dp{-2.0}};
                // auto mv_rot = mvec3dp{scalar3dp{4.0}, c, b, t, pscalar3dp{-2.0}};
                // CHECK(mv_rot == move3dp(mv_rev, get_motor(e43_3dp, phi)));
            }

            // check same rotation applied to many points
            std::vector<Vec3dp<double>> vec_ref(100);
            std::vector<Vec3dp<double>> vec_rot(100);
            std::vector<Vec3dp<double>> vec_rot_calc(100);

            // prepare for checking vector-based rotation (after the loop)
            double phi = pi / 12;
            for (size_t i = 0; i < 100; ++i) {
                vec_ref.push_back(e1_3dp);
                vec_rot.emplace_back(move3dp(e1_3dp, get_motor(e43_3dp, phi)));
            }

            vec_rot_calc = move3dp(vec_ref, get_motor(e43_3dp, phi));
            for (size_t i = 0; i < 100; ++i) {
                CHECK(vec_rot_calc[i] == vec_rot[i]);
            }
        }

        fmt::println("");
    }


    TEST_CASE("Vec3dp: operations - simple rotation")
    {
        fmt::println("Vec3dp: operations - simple rotation");

        // rotation in e431 plane, normalized rotation axis is showing in e2 direction
        auto const axis = unitize(wdg(vec3dp{1, 0, 1, 1}, vec3dp{1, 1, 1, 1}));
        auto const phi = deg2rad(90);
        auto const X0 = vec3dp{0, 0, 2, 1};
        auto const dx = vec3dp{1, 0, -1, 0};
        auto const nrm_dx = bulk_nrm(dx);

        auto const M_rot = get_motor(axis, phi);
        auto const M_tra = get_motor(dx);
        auto const X_rot = move3dp(X0, M_rot);
        auto const X_tra = move3dp(X0, M_tra);

        fmt::println("axis = {}", axis);
        fmt::println("phi  = {}", phi);
        fmt::println("X0   = {}", X0);
        fmt::println("");
        fmt::println("M_rot    = {}", M_rot);
        fmt::println("X_rot    = {}", X_rot);
        fmt::println("");
        fmt::println("M_tra    = {}", M_tra);
        fmt::println("X_tra    = {}", X_tra);
        fmt::println("");
        fmt::println("phi/2*I_3dp = {}", phi / 2 * I_3dp);
        fmt::println("rgpr(phi/2*I_3dp,axis) = {}", rgpr(phi / 2 * I_3dp, axis));
        fmt::println("rgpr(axis,axis) = {}", rgpr(axis, axis));
        fmt::println("");
        fmt::println("rgpr(nrm_dx/2,axis) = {}", rgpr(nrm_dx / 2.0, axis));
        fmt::println("r_weight_dual(axis) = {}", r_weight_dual(axis));
        fmt::println("");

        CHECK(X_rot == vec3dp{2.0, 0.0, 2.0, 1.0});
        CHECK(X_tra == vec3dp{1.0, 0.0, 1.0, 1.0});
        //
        CHECK(rgpr(axis, rrev(axis)) == mvec3dp_e(pscalar3dp(1.0)));
        CHECK(rgpr(axis, axis) == mvec3dp_e(pscalar3dp(-1.0)));

        fmt::println("");
    }

    TEST_CASE("Vec3dp: operations - vanishing points, plane at infinity")
    {
        fmt::println("Vec3dp: operations - vanishing points, plane at infinity");

        // simple line
        auto const line = unitize(wdg(vec3dp{1, 0, 1, 1}, vec3dp{1, 1, 1, 1}));

        // horizon = plane at infinity (H_3dp = rmpl(O_3dp) = e321_3dp)

        // vanishing point vp = intersection of line with horizon plane
        auto vp = rgpr(line, H_3dp);

        fmt::println("");
        fmt::println("vp      = {}", vp);
        fmt::println("gr1(vp) = {}", gr1(vp)); // this is what we are looking for
        fmt::println("gr3(vp) = {}", gr3(vp));

        CHECK(gr1(vp) == vec3dp{0.0, 1.0, 0.0, 0.0}); // e2 direction
        CHECK(gr3(vp) == trivec3dp{});                // no plane part

        fmt::println("");
    }

    TEST_CASE("Vec3dp: operations - different kinds of lines")
    {
        fmt::println("Vec3dp: operations - different kinds of lines");

        // simple line defined by two points in e1 direction
        auto const l1 = unitize(wdg(vec3dp{0, 1, -1, 1}, vec3dp{1, 1, -1, 1}));

        // line defined by containing the first point and the direction is identical
        // (what Dorst calls a "vanishing line", because it contains a point at infinity)
        CHECK(unitize(wdg(vec3dp{0, 1, -1, 1}, e1_3dp)) == l1);

        // simple line in e2 direction
        auto const l2 = unitize(wdg(vec3dp{1, 0, 1, 1}, vec3dp{1, 1, 1, 1}));
        CHECK(unitize(wdg(vec3dp{1, 0, 1, 1}, e2_3dp)) == l2);

        // line l1 + l2
        auto const l3 = unitize(l1 + l2);

        // let's take an non-blade bivector now
        auto const l4 = wdg(e1_3dp - e3_3dp, e2_3dp) + wdg(e4_3dp, e3_3dp);

        // horizon = plane at infinity (H_3dp = rmpl(O_3dp) = e321_3dp)
        // vanishing point vp1 = intersection of line with horizon plane
        // => results in a direction towards infinity, not a localized point!
        auto vp1 = rgpr(l1, H_3dp);
        auto vp2 = rgpr(l2, H_3dp);
        auto vp3 = rgpr(l3, H_3dp);

        fmt::println("");
        fmt::println("l1  = {: 5.3f}, att(l1) = {: 5.3f}, sup(l1) = {: 5.3f}", l1,
                     att(l1), sup(l1));
        fmt::println("l2  = {: 5.3f}, att(l2) = {: 5.3f}, sup(l2) = {: 5.3f}", l2,
                     att(l2), sup(l2));
        fmt::println("l3  = {: 5.3f}, att(l3) = {: 5.3f}, sup(l3) = {: 5.3f}", l3,
                     att(l3), sup(l3));
        fmt::println("l4  = {: 5.3f}, att(l4) = {: 5.3f}, sup(l4) = {: 5.3f}", l4,
                     att(l4), sup(l4));
        fmt::println("");
        fmt::println("vp1 = {: 5.3f}", vp1);
        fmt::println("vp2 = {: 5.3f}", vp2);
        fmt::println("vp3 = {: 5.3f}", vp3);
        fmt::println("");
        fmt::println("rgpr(l1,l1) = {: 5.3f}", rgpr(l1, l1));
        fmt::println("rgpr(l1,rrev(l1)) = {: 5.3f}", rgpr(l1, rrev(l1)));
        fmt::println("rwdg(l1,l1) = {: 5.3f}", rwdg(l1, l1));
        fmt::println("rcmt(l1,l1) = {: 5.3f}", rcmt(l1, l1));
        fmt::println("rdot(l1,l1) = {: 5.3f}", rdot(l1, l1));
        fmt::println("");
        fmt::println("rgpr(l4,l4) = {: 5.3f}", rgpr(l4, l4));
        fmt::println("rgpr(l4,rrev(l4)) = {: 5.3f}", rgpr(l4, rrev(l4)));
        fmt::println("rwdg(l4,l4) = {: 5.3f}", rwdg(l4, l4));
        fmt::println("rcmt(l4,l4) = {: 5.3f}", rcmt(l4, l4));
        fmt::println("rdot(l4,l4) = {: 5.3f}", rdot(l4, l4));
        fmt::println("");

        CHECK(gr1(vp1) == e1_3dp);
        CHECK(gr1(vp2) == e2_3dp);

        fmt::println("");
    }

    TEST_CASE("Vec3dp: operations - test point lies in plane")
    {
        fmt::println("Vec3dp: operations - test point in plane");

        // plane is e431 plane (bivector of the plane is e31)
        auto const plane1 = e431_3dp;
        // plane is parallel to e431 plane (defined by three points in the new plane)
        auto const plane2 = unitize(
            wdg(wdg(vec3dp{0, 1, 0, 1}, vec3dp{1, 1, 0, 1}), vec3dp{1, 1, -1, 1}));

        // alternative way to shift the plane:
        // shift a plane with its reference point in a direction x,
        // while keeping the same orientation:
        // 1.) pl = P ^ u ^ v
        // 2.) pl_shifted = P ^ u ^ v + x ^ u ^ v = (P + x) ^ (u ^ v)
        //
        // see Browne, "Grassmann Algebra": a vector acts as carrier of points
        //                                  a bivector acts as carrier of lines
        //                                  a trivector acts as carrier of planes
        auto const x = vec3dp{0, 1, 0, 0};
        CHECK(plane2 == plane1 + wdg(x, e31_3dp));

        // point P
        auto const P1 = vec3dp{2.0, 0.0, 1.0, 1.0};
        auto const P2 = vec3dp{2.0, 1.0, 1.0, 1.0};

        // P = O + p is in plane, if wp == 0 (=> p is part of the trivector of the plane)
        auto wp11 = wdg(plane1, P1);
        auto wp12 = wdg(plane1, P2);
        auto wp21 = wdg(plane2, P1);
        auto wp22 = wdg(plane2, P2);

        fmt::println("");
        fmt::println("plane1      = {}", plane1);
        fmt::println("plane2      = {}", plane2);
        fmt::println("");
        fmt::println("P1      = {}", P1);
        fmt::println("P2      = {}", P2);
        fmt::println("");
        fmt::println("wp11      = {}", wp11);
        fmt::println("wp12      = {}", wp12);
        fmt::println("");
        fmt::println("wp21      = {}", wp21);
        fmt::println("wp22      = {}", wp22);

        CHECK(wp11 == 0.0); // == works because scalar_t can be compared to raw values
        CHECK(wp12 != 0.0);
        CHECK(wp21 != 0.0);
        CHECK(wp22 == 0.0);

        fmt::println("");
    }

    TEST_CASE("Vec3dp: operations - rotation + translation combined")
    {
        fmt::println("Vec3dp: operations - rotation + translation combined");

        // three points representing two lines l0=(p0,p1), l1 = (p0,p2)
        auto const p0 = vec3dp{0, 0, -1, 1};
        auto const p1 = vec3dp{1, 0, 1, 1};
        auto const p2 = vec3dp{0, 1, 0, 1};

        auto vp = std::vector<vec3dp>{p0, p1, p2};


        auto phi = deg2rad(15);           // turning angle
        auto B = e42_3dp + 0.5 * e43_3dp; // turning axis (a bivector)
        auto t = vec3dp{-2, 1, 1, 0};     // translation vector

        // stepwise movement: turn around turning axis by 15°, then shift by t
        auto m_rot = get_motor(B, phi);
        auto m_tra = get_motor(t);

        auto m_rt_c = rgpr(m_tra, m_rot); // concatenation: translation after rotation
        auto m_tr_c = rgpr(m_rot, m_tra); // concatenation: rotation after translation
                                          // sequence matters!

        auto vr = move3dp(vp, m_rot);
        auto vt = move3dp(vp, m_tra);
        auto vrt = move3dp(vr, m_tra);
        auto vtr = move3dp(vt, m_rot);
        auto vrt_c = move3dp(vp, m_rt_c);
        auto vtr_c = move3dp(vp, m_tr_c);

        fmt::println("m_rot  = {: 5.3f}", m_rot);
        fmt::println("m_tra  = {: 5.3f}", m_tra);
        fmt::println("m_rt_c = rgpr(m_tra,m_rot)  = {: 5.3f}", m_rt_c);
        fmt::println("m_tr_c = rgpr(m_rot,m_tra)  = {: 5.3f}", m_tr_c);
        fmt::println("");
        fmt::println("vp     = {: 5.3f}", fmt::join(vp, ", "));
        fmt::println("vr     = {: 5.3f}", fmt::join(vr, ", "));
        fmt::println("vt     = {: 5.3f}", fmt::join(vt, ", "));
        fmt::println("");
        fmt::println("vrt    = {: 5.3f}", fmt::join(vrt, ", "));
        fmt::println("vrt_c  = {: 5.3f}", fmt::join(vrt_c, ", "));
        fmt::println("");
        fmt::println("vtr    = {: 5.3f}", fmt::join(vtr, ", "));
        fmt::println("vtr_c  = {: 5.3f}", fmt::join(vtr_c, ", "));
        fmt::println("");

        CHECK(vrt == vrt_c); // same result stepwise or combined
        CHECK(vtr == vtr_c); // same result stepwise or combined

        fmt::println("m_rot                          = {: 5.3f}", m_rot);
        fmt::println("rrev(m_rot)                    = {: 5.3f}", rrev(m_rot));
        fmt::println("rgpr(m_rot, rrev(m_rot))       = {: 5.3f}",
                     rgpr(m_rot, rrev(m_rot)));
        fmt::println("gr4(rgpr(m_rot, rrev(m_rot)))  = {: 5.3f}",
                     rgpr(m_rot, rrev(m_rot)));
        fmt::println("gr4(rgpr(rrev(m_rot),m_rot))  = {: 5.3f}",
                     rgpr(rrev(m_rot), m_rot));
        // rotors in pga multiply to the pseudoscalar as neutral element of rgpr(),
        // not to 1 (which is the neutral element of gpr()) like in ega
        CHECK(gr0(rgpr(m_rot, rrev(m_rot))) == scalar3dp{});
        CHECK(gr2(rgpr(m_rot, rrev(m_rot))) == bivec3dp{});
        CHECK(gr4(rgpr(m_rot, rrev(m_rot))) == pscalar3dp{1.0});

        CHECK(gr0(rgpr(rrev(m_rot), m_rot)) == scalar3dp{});
        CHECK(gr2(rgpr(rrev(m_rot), m_rot)) == bivec3dp{});
        CHECK(gr4(rgpr(rrev(m_rot), m_rot)) == pscalar3dp{1.0});


        fmt::println("");
    }

    TEST_CASE("Vec3dp: operations - rotation + translation combined (screw motion)")
    {
        fmt::println(
            "Vec3dp: operations - rotation + translation combined (screw motion)");

        auto axis = x_axis_3dp + 0.3 * y_axis_3dp; // axis direction
        auto theta = deg2rad(180);                 // = will turn to neg. direction
        auto d = 2.0;                              // = dist d for half a turn
        auto t =
            bulk_normalize(att(axis)) * d; // vector t in direction of axis with length d

        auto vp = vec3dp{0, 1, 0, 1}; // reference point for transformation

        auto m_rot = get_motor(axis, theta);
        auto m_tra = get_motor(t);
        auto m_c = get_motor(axis, theta, d);

        auto m_rt_c = rgpr(m_tra, m_rot); // concatenation: translation after rotation
        auto m_tr_c = rgpr(m_rot, m_tra); // concatenation: rotation after translation
                                          // sequence does not matter for screw motion

        auto vr = move3dp(vp, m_rot);
        auto vt = move3dp(vp, m_tra);
        auto vrt = move3dp(vr, m_tra);
        auto vtr = move3dp(vt, m_rot);
        auto vrt_c = move3dp(vp, m_rt_c);
        auto vtr_c = move3dp(vp, m_tr_c);
        auto v_c = move3dp(vp, m_c);

        fmt::println("m_rot                      = {: 5.3f}", m_rot);
        fmt::println("m_tra                      = {: 5.3f}", m_tra);
        fmt::println("m_rt_c = rgpr(m_tra,m_rot) = {: 5.3f}", m_rt_c);
        fmt::println("m_tr_c = rgpr(m_rot,m_tra) = {: 5.3f}", m_tr_c);
        fmt::println("m_c                        = {: 5.3f}", m_c);
        fmt::println("");
        fmt::println("vp     = {: 5.3f}", vp, ", ");
        fmt::println("vr     = {: 5.3f}", vr, ", ");
        fmt::println("t      = {: 5.3f}", t, ", ");
        fmt::println("vt     = {: 5.3f}", vt, ", ");
        fmt::println("");
        fmt::println("vrt    = {: 5.3f}", vrt, ", ");
        fmt::println("vrt_c  = {: 5.3f}", vrt_c, ", ");
        fmt::println("vtr    = {: 5.3f}", vtr, ", ");
        fmt::println("vrt_c  = {: 5.3f}", vtr_c, ", ");
        fmt::println("");
        fmt::println("v_c    = {: 5.3f}", v_c, ", ");

        CHECK(vrt == vrt_c); // same result stepwise or combined
        CHECK(vtr == vtr_c); // same result stepwise or combined
        CHECK(vtr_c == v_c); // same result stepwise or combined
        CHECK(vrt_c == v_c); // same result stepwise or combined

        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3dp: default init")
    {
        fmt::println("MVec3dp: default init");
        // default initialization
        mvec3dp v;
        // fmt::println("   v = {}", v);
        CHECK(abs(v.c0) < eps);
        CHECK(abs(v.c1) < eps);
        CHECK(abs(v.c2) < eps);
        CHECK(abs(v.c3) < eps);
        CHECK(abs(v.c4) < eps);
        CHECK(abs(v.c5) < eps);
        CHECK(abs(v.c6) < eps);
        CHECK(abs(v.c7) < eps);
        CHECK(abs(v.c8) < eps);
        CHECK(abs(v.c9) < eps);
        CHECK(abs(v.c10) < eps);
        CHECK(abs(v.c11) < eps);
        CHECK(abs(v.c12) < eps);
        CHECK(abs(v.c13) < eps);
        CHECK(abs(v.c14) < eps);
        CHECK(abs(v.c15) < eps);
    }
    TEST_CASE("MVec3dp: with curly braced intializer")
    {
        fmt::println("MVec3dp: with curly braced intializer");
        // default initialization
        mvec3dp v{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0,
                  0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        // fmt::println("   v = {}", v);
        CHECK(abs(v.c0 - 0.0) < eps);
        CHECK(abs(v.c1 - 1.0) < eps);
        CHECK(abs(v.c2 - 2.0) < eps);
        CHECK(abs(v.c3 - 3.0) < eps);
        CHECK(abs(v.c4 - 23.0) < eps);
        CHECK(abs(v.c5 - 31.0) < eps);
        CHECK(abs(v.c6 - 12.0) < eps);
        CHECK(abs(v.c7 - 123.0) < eps);
        CHECK(abs(v.c8 - 0.0) < eps);
        CHECK(abs(v.c9 - 1.0) < eps);
        CHECK(abs(v.c10 - 2.0) < eps);
        CHECK(abs(v.c11 - 3.0) < eps);
        CHECK(abs(v.c12 - 23.0) < eps);
        CHECK(abs(v.c13 - 31.0) < eps);
        CHECK(abs(v.c14 - 12.0) < eps);
        CHECK(abs(v.c15 - 123.0) < eps);
    }

    TEST_CASE("MVec3dp: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec3dp: cp ctor & cp assign incl. type deduction");
        // default initialization
        mvec3dp v1{
            0.0, 1.0, 2.0, 3.0,  23.0, 31.0, 12.0, 123.0, 0.0,
            1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0}; // init with double (type deduction)

        mvec3dp v2{v1};   // cp ctor
        mvec3dp v3 = v2;  // cp assign
        mvec3dp v4 = -v3; // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(abs(v2.c0 - 0.0) < eps);
        CHECK(abs(v2.c1 - 1.0) < eps);
        CHECK(abs(v2.c2 - 2.0) < eps);
        CHECK(abs(v2.c3 - 3.0) < eps);
        CHECK(abs(v2.c4 - 23.0) < eps);
        CHECK(abs(v2.c5 - 31.0) < eps);
        CHECK(abs(v2.c6 - 12.0) < eps);
        CHECK(abs(v2.c7 - 123.0) < eps);
        CHECK(abs(v2.c8 - 0.0) < eps);
        CHECK(abs(v2.c9 - 1.0) < eps);
        CHECK(abs(v2.c10 - 2.0) < eps);
        CHECK(abs(v2.c11 - 3.0) < eps);
        CHECK(abs(v2.c12 - 23.0) < eps);
        CHECK(abs(v2.c13 - 31.0) < eps);
        CHECK(abs(v2.c14 - 12.0) < eps);
        CHECK(abs(v2.c15 - 123.0) < eps);

        CHECK(abs(v3.c0 - 0.0) < eps);
        CHECK(abs(v3.c1 - 1.0) < eps);
        CHECK(abs(v3.c2 - 2.0) < eps);
        CHECK(abs(v3.c3 - 3.0) < eps);
        CHECK(abs(v3.c4 - 23.0) < eps);
        CHECK(abs(v3.c5 - 31.0) < eps);
        CHECK(abs(v3.c6 - 12.0) < eps);
        CHECK(abs(v3.c7 - 123.0) < eps);
        CHECK(abs(v3.c8 - 0.0) < eps);
        CHECK(abs(v3.c9 - 1.0) < eps);
        CHECK(abs(v3.c10 - 2.0) < eps);
        CHECK(abs(v3.c11 - 3.0) < eps);
        CHECK(abs(v3.c12 - 23.0) < eps);
        CHECK(abs(v3.c13 - 31.0) < eps);
        CHECK(abs(v3.c14 - 12.0) < eps);
        CHECK(abs(v3.c15 - 123.0) < eps);

        CHECK(v4 == -v3);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }

    TEST_CASE("MVec3dp: defining basic types and ctor checks")
    {
        fmt::println("MVec3dp: defining basic types and ctor checks");

        auto s = scalar3dp{5.0};
        auto v = vec3dp{1.0, 2.0, 3.0, 1.0};
        auto B = bivec3dp{-1.0, 2.0, 1.0, -10.0, 20.0, 10.0};
        auto t = trivec3dp{3.0, 6.0, 9.0, 3.0};
        auto ps = pscalar3dp{-5.0};

        auto mv0 = mvec3dp{s};
        auto mv1 = mvec3dp{v};
        auto mv2 = mvec3dp{B};
        auto mv3 = mvec3dp{t};
        auto mv4 = mvec3dp{ps};

        auto mv5a = mvec3dp{s, B, ps};
        auto mv5b = mvec3dp_e{s};
        auto mv5c = mvec3dp_e{B};
        auto mv5d = mvec3dp_e{ps};
        auto mv5e = mvec3dp_e{s, B, ps};
        auto mv5f{mv5e};           // cp assign to even grade
        auto mv5g = mvec3dp{mv5e}; // assign to mv

        auto mv5h = mvec3dp_e{s, B};
        auto mv5i = mvec3dp_e{B, ps};
        auto mv5j = mvec3dp_e{s, ps};

        auto mv6a = mvec3dp{v, t};
        auto mv6b = mvec3dp_u{v};
        auto mv6c = mvec3dp_u{t};
        auto mv6d = mvec3dp_u{v, t};
        auto mv6e{mv6d};           // cp assign to odd grade
        auto mv6f = mvec3dp{mv6d}; // assign to mv

        auto mv7 = mvec3dp{s, v, B, t, ps};

        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv4  = {}", mv4);
        // fmt::println("   mv5a = {}", mv5a);
        // fmt::println("   mv5b = {}", mv5b);
        // fmt::println("   mv5c = {}", mv5c);
        // fmt::println("   mv5d = {}", mv5d);
        // fmt::println("   mv5e = {}", mv5e);
        // fmt::println("   mv5f = {}", mv5f);
        // fmt::println("   mv5g = {}", mv5g);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6b = {}", mv6b);
        // fmt::println("   mv6c = {}", mv6c);
        // fmt::println("   mv6d = {}", mv6d);
        // fmt::println("   mv6e = {}", mv6e);
        // fmt::println("   mv6f = {}", mv6f);
        // fmt::println("   mv7  = {}", mv7);

        CHECK(gr0(mv0) == s);
        CHECK(gr1(mv1) == v);
        CHECK(gr2(mv2) == B);
        CHECK(gr3(mv3) == t);
        CHECK(gr4(mv4) == ps);

        CHECK(gr0(mv5a) == s);
        CHECK(gr1(mv5a) == vec3dp{});
        CHECK(gr2(mv5a) == B);
        CHECK(gr3(mv5a) == trivec3dp{});
        CHECK(gr4(mv5a) == ps);

        CHECK(gr0(mv5b) == s);
        CHECK(gr2(mv5b) == bivec3dp{});
        CHECK(gr4(mv5b) == pscalar3dp{});

        CHECK(gr0(mv5c) == scalar3dp{});
        CHECK(gr2(mv5c) == B);
        CHECK(gr4(mv5c) == pscalar3dp{});

        CHECK(gr0(mv5d) == scalar3dp{});
        CHECK(gr2(mv5d) == bivec3dp{});
        CHECK(gr4(mv5d) == ps);

        CHECK(gr0(mv5e) == s);
        CHECK(gr2(mv5e) == B);
        CHECK(gr4(mv5e) == ps);

        CHECK(mv5f == mv5e);

        CHECK(gr0(mv5g) == s);
        CHECK(gr2(mv5g) == B);
        CHECK(gr4(mv5g) == ps);

        CHECK(gr0(mv5h) == s);
        CHECK(gr2(mv5h) == B);
        CHECK(gr4(mv5h) == pscalar3dp{});

        CHECK(gr0(mv5i) == scalar3dp{});
        CHECK(gr2(mv5i) == B);
        CHECK(gr4(mv5i) == ps);

        CHECK(gr0(mv5j) == s);
        CHECK(gr2(mv5j) == bivec3dp{});
        CHECK(gr4(mv5j) == ps);

        CHECK(gr0(mv6a) == scalar3dp{});
        CHECK(gr1(mv6a) == v);
        CHECK(gr2(mv6a) == bivec3dp{});
        CHECK(gr3(mv6a) == t);
        CHECK(gr4(mv6a) == pscalar3dp{});

        CHECK(gr1(mv6b) == v);
        CHECK(gr3(mv6b) == trivec3dp{});

        CHECK(gr1(mv6c) == vec3dp{});
        CHECK(gr3(mv6c) == t);

        CHECK(gr1(mv6d) == v);
        CHECK(gr3(mv6d) == t);

        CHECK(mv6e == mv6d);

        CHECK(gr0(mv6f) == scalar3dp{});
        CHECK(gr1(mv6f) == v);
        CHECK(gr2(mv6f) == bivec3dp{});
        CHECK(gr3(mv6f) == t);
        CHECK(gr4(mv6f) == pscalar3dp{});

        CHECK(gr0(mv7) == s);
        CHECK(gr1(mv7) == v);
        CHECK(gr2(mv7) == B);
        CHECK(gr3(mv7) == t);
        CHECK(gr4(mv7) == ps);

        CHECK(gr(s) == 0);
        CHECK(gr(v) == 1);
        CHECK(gr(B) == 2);
        CHECK(gr(t) == 3);
        CHECK(gr(ps) == 4);
    }

    TEST_CASE("MVec3dp: fmt & cout printing")
    {
        fmt::println("MVec3dp: fmt & cout printing");

        mvec3dp pf{1.0f, 2.00001f, 0.0f, 3.0f, 1.0f, 2.00001f, 0.0f, 3.0f,
                   1.0f, 2.00001f, 1.0f, 3.0f, 1.0f, 2.00001f, 2.0f, 3.0f};
        mvec3dp pd{1.0, 2.00001, 0.0, 3.0, 1.0, 2.00001, 0.0, 3.0,
                   1.0, 2.00001, 1.0, 3.0, 1.0, 2.00001, 2.0, 3.0};

        fmt::println("");
        std::cout << "    cout: pf = " << pf << std::endl;
        fmt::println("    fmt:  pf = {}", pf);
        fmt::println("    fmt:  pf = {:.8f}", pf);

        std::cout << "    cout: pd = " << pd << std::endl;
        fmt::println("    fmt:  pd = {}", pd);
        fmt::println("    fmt:  pd = {:.8f}", pd);

        std::vector<MVec3dp<double>> vp1{{1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 10.0,
                                          10.0, 10.0, 20.0, 10.0, 10.0, 10.0, 20.0},
                                         {0.5, 1.5, 2.0, 2.5, 1.0, 1.0, 1.0, 2.0, 10.5,
                                          10.5, 20.0, 20.5, 10.0, 10.0, 10.0, 20.0}};
        fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("    fmt: vp1 = {:e}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(geom_nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("MVec3dp: vector space and linearity tests")
    {
        fmt::println("MVec3dp: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        mvec3dp p0;
        mvec3dp p1{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0,
                   0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
        mvec3dp p2{0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0,
                   0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0};
        mvec3dp p3{0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0,
                   0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0};
        mvec3dp p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t.multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t. addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
        CHECK(p4 == -p1);
    }

    TEST_CASE("MVec3dp: wedge product - basic properties")
    {
        fmt::println("MVec3dp: wedge product - basic properties");

        vec3dp v1{1.0, 2.0, 1.0, 1.0};
        vec3dp v2{0.5, 3.0, 2.0, 1.0};
        vec3dp v3{-2.0, 6.0, 3.0, 1.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3dp{sd};
        auto t = scalar3dp{st};

        CHECK(wdg(v1, v1) == bivec3dp{});                     // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s));   // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, v1) == wdg(v1, s)); // wdg between scalar and vector
        CHECK(wdg(s, v1) == sd * v1);    // wdg between scalar and vector
    }

    TEST_CASE("MVec3dp: geometric product - basic properties")
    {
        fmt::println("MVec3dp: geometric product - basic properties");

        vec3dp v1{1.0, 2.0, 1.0, 1.0};
        vec3dp v2{0.5, 3.0, 2.0, 1.0};
        vec3dp v3{-2.0, 6.0, 3.0, 1.0};
        bivec3dp B1{-4, 2, 1, 1, -1, 2};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3dp{sd};
        auto t = scalar3dp{st};

        // bulk_nrm_sq(v1) = dot(v1,v1)  for every projective vector
        CHECK(bulk_nrm_sq(v1) == value_t(dot(v1, v1)));

        CHECK(v1 * v2 == dot(v1, v2) + wdg(v1, v2)); // valid for vectors only
        CHECK(v1 * v2 == (v2 >> v1) + wdg(v1, v2));  // contraction = dot for same grades

        // valid also for gr(B1) != gr(v1)
        CHECK(v1 * B1 == (B1 >> v1) + wdg(v1, B1));
        CHECK(B1 * v1 == (v1 << B1) + wdg(B1, v1));

        // only valid for vectors (!):
        // dot = gpr symmetric part
        CHECK(dot(v1, v2) == gr0(0.5 * (v1 * v2 + v2 * v1)));
        // wdg = gpr anti-symmetric part
        CHECK(wdg(v1, v2) == gr2(0.5 * (v1 * v2 - v2 * v1)));

        // mathematical characteristics
        CHECK((v1 * v2) * v3 == v1 * (v2 * v3));    // gpr is associative
        CHECK(v1 * (v2 + v3) == v1 * v2 + v1 * v3); // gpr distributes over addition
        CHECK((v1 + v2) * v3 == v1 * v3 + v2 * v3); // wdg distributes over addition
        CHECK((sd * v1) * v2 == v1 * (sd * v2));    // scalars can be factored out of gpr
        CHECK((sd * v1) * v2 == sd * (v1 * v2));    // scalars can be factored out of gpr
        CHECK(s * t == t * s);    // gpr between scalars equivalent to scalar mult.
        CHECK(s * v1 == v1 * s);  // gpr between scalar and vector
        CHECK(s * v1 == sd * v1); // gpr between scalar and vector
    }

    TEST_CASE("MVec3dp: geometric product - combinatorial tests")
    {
        fmt::println("MVec3dp: geometric product - combinatorial tests");

        scalar3dp s1{3.0};
        scalar3dp s2{-1.5};
        vec3dp v1{1.0, 2.0, 3.0, 4.0};
        vec3dp v2{0.5, 3.0, -2.0, 1.0};
        bivec3dp b1{-1.0, -3.0, -5.0, -10.0, -30.0, -50.0};
        bivec3dp b2{2.0, 0.5, -4.0, 20.0, 5.0, -40.0};
        trivec3dp t1{1.0, 2.0, 3.0, 4.0};
        trivec3dp t2{0.5, 3.0, -2.0, 1.0};
        pscalar3dp ps1{-5.0};
        pscalar3dp ps2{2.0};

        // checks ctor combinations and equivalence of simplified products with
        // geometric product of a fully populated multivector
        //
        // assumption: the geometric product of a fully populated multivector serves as
        //             reference and is correctly implemented
        //             results are relative to full geometric product

        // MVec3dp_E * MVec3dp
        CHECK(mvec3dp{s1, b1, ps1} * mvec3dp{s1, v1, b1, t1, ps1} ==
              mvec3dp{mvec3dp_e{s1, b1, ps1} * mvec3dp{s1, v1, b1, t1, ps1}});

        // MVec3dp * MVec3dp_E
        CHECK(mvec3dp{s1, v1, b1, t1, ps1} * mvec3dp{s1, b1, ps1} ==
              mvec3dp{mvec3dp{s1, v1, b1, t1, ps1} * mvec3dp_e{s1, b1, ps1}});

        // MVec3dp_E * MVec3dp_E
        CHECK(mvec3dp{s1, b1, ps1} * mvec3dp{s2, b2, ps2} ==
              mvec3dp{mvec3dp_e{s1, b1, ps1} * mvec3dp_e{s2, b2, ps2}});

        // MVec3dp_U * MVec3dp_U
        CHECK(mvec3dp{v1, t1} * mvec3dp{v2, t2} ==
              mvec3dp{mvec3dp_u{v1, t1} * mvec3dp_u{v2, t2}});

        // MVec3dp_E * MVec3dp_U
        CHECK(mvec3dp{s1, b1, ps1} * mvec3dp{v2, t2} ==
              mvec3dp{mvec3dp_e{s1, b1, ps1} * mvec3dp_u{v2, t2}});

        // MVec3dp_U * MVec3dp_E
        CHECK(mvec3dp{v1, t1} * mvec3dp{s2, b2, ps2} ==
              mvec3dp{mvec3dp_u{v1, t1} * mvec3dp_e{s2, b2, ps2}});

        // MVec3dp_E * TriVec3dp
        CHECK(mvec3dp{s1, b1, ps1} * mvec3dp{t2} == mvec3dp{mvec3dp_e{s1, b1, ps1} * t2});

        // MVec3dp_E * BiVec3dp
        CHECK(mvec3dp{s1, b1, ps1} * mvec3dp{b2} == mvec3dp{mvec3dp_e{s1, b1, ps1} * b2});

        // BiVec3dp * MVec3dp_E
        CHECK(mvec3dp{b1} * mvec3dp{s2, b2, ps2} ==
              mvec3dp{mvec3dp{b1} * mvec3dp_e{s2, b2, ps2}});

        // MVec3dp_E * Vec3dp
        CHECK(mvec3dp{s1, b1, ps1} * mvec3dp{v2} == mvec3dp{mvec3dp_e{s1, b1, ps1} * v2});

        // Vec3dp * MVec3dp_E
        CHECK(mvec3dp{v1} * mvec3dp{s2, b2, ps2} == mvec3dp{v1 * mvec3dp_e{s2, b2, ps2}});

        // MVec3dp_U * BiVec3dp
        CHECK(mvec3dp{v1, t1} * mvec3dp{b2} == mvec3dp{mvec3dp_u{v1, t1} * b2});

        // BiVec3dp * MVec3dp_U
        CHECK(mvec3dp{b1} * mvec3dp{v2, t2} == mvec3dp{b1 * mvec3dp_u{v2, t2}});

        // TriVec3dp * TriVec3dp
        CHECK(mvec3dp{t1} * mvec3dp{t2} == mvec3dp{t1 * t2});

        // TriVec3dp * BiVec3dp
        CHECK(mvec3dp{t1} * mvec3dp{b1} == mvec3dp{t1 * b1});

        // BiVec3dp * TriVec3dp
        CHECK(mvec3dp{b1} * mvec3dp{t1} == mvec3dp{b1 * t1});

        // TriVec3dp * Vec3dp
        CHECK(mvec3dp{t1} * mvec3dp{v1} == mvec3dp{t1 * v1});

        // Vec3dp * TriVec3dp
        CHECK(mvec3dp{v1} * mvec3dp{t1} == mvec3dp{v1 * t1});

        // BiVec3dp * BiVec3dp
        CHECK(mvec3dp{b1} * mvec3dp{b2} == mvec3dp{b1 * b2});

        // BiVec3dp * Vec3dp
        CHECK(mvec3dp{b1} * mvec3dp{v1} == mvec3dp{b1 * v1});

        // Vec3dp * BiVec3dp
        CHECK(mvec3dp{v1} * mvec3dp{b1} == mvec3dp{v1 * b1});

        // Vec3dp * Vec3dp
        CHECK(mvec3dp{v1} * mvec3dp{v2} == mvec3dp{v1 * v2});

        // PScalar3dp * MVec3dp
        CHECK(mvec3dp{ps1} * mvec3dp{s2, v2, b2, t2, ps2} ==
              ps1 * mvec3dp{s2, v2, b2, t2, ps2});

        // MVec3dp * Pscalar3dp
        CHECK(mvec3dp{s1, v1, b1, t1, ps1} * mvec3dp{ps2} ==
              mvec3dp{mvec3dp{s1, v1, b1, t1, ps1} * ps2});

        // PScalar3dp * MVec3dp_E
        CHECK(mvec3dp{ps1} * mvec3dp{s2, b2, ps2} ==
              mvec3dp{ps1 * mvec3dp_e{s2, b2, ps2}});

        // MVec3dp_E * Pscalar3dp
        CHECK(mvec3dp{s1, b1, ps1} * mvec3dp{ps2} ==
              mvec3dp{mvec3dp_e{s1, b1, ps1} * ps2});

        // PScalar3dp * MVec3dp_U
        CHECK(mvec3dp{ps1} * mvec3dp{v2, t2} == mvec3dp{ps1 * mvec3dp_u{v2, t2}});

        // MVec3dp_U * Pscalar3dp
        CHECK(mvec3dp{v1, t1} * mvec3dp{ps2} == mvec3dp{mvec3dp_u{v1, t1} * ps2});

        // TODO: add geometric products including pseudoscalars

        // // PScalar3dp * TriVec3dp
        // CHECK(mvec3dp{ps1} * mvec3dp{t2} == mvec3dp{ps1 * t2});

        // // TriVec3dp * Pscalar3dp
        // CHECK(mvec3dp{t1} * mvec3dp{ps2} == mvec3dp{t1 * ps2});

        // // PScalar3dp * BiVec3dp
        // CHECK(mvec3dp{ps1} * mvec3dp{b2} == mvec3dp{ps1 * b2});

        // // BiVec3dp * Pscalar3dp
        // CHECK(mvec3dp{b1} * mvec3dp{ps2} == mvec3dp{b1 * ps2});

        // // PScalar3dp * Vec3dp
        // CHECK(mvec3dp{ps1} * mvec3dp{v2} == mvec3dp{ps1 * v2});

        // // Vec3dp * Pscalar3dp
        // CHECK(mvec3dp{v1} * mvec3dp{ps2} == mvec3dp{v1 * ps2});

        // // Pscalar3dp * Pscalar3dp
        // CHECK(mvec3dp{ps1} * mvec3dp{ps2} == mvec3dp{ps1 * ps2});

        // Scalar * Scalar
        CHECK(mvec3dp{s1} * mvec3dp{s2} == mvec3dp{s1 * s2});
    }

    TEST_CASE("MVec3dp: geometric product tests - vec * vec")
    {
        fmt::println("MVec3dp: geometric product tests - vec * vec");

        // ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
        //
        // dot(a,b) = 0.5*(ab + ba)   (symmetric part)
        // wdg(a,b) = 0.5*(ab - ba)   (antisymmetric part)

        vec3dp a{1, 2, 3, 1};
        vec3dp b{0.5, 3, -2, 1};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec3dp mva{a};
        mvec3dp mvb{b};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb + mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb - mvb * mva);

        CHECK(dot_ab == gr0(mvab));
        CHECK(dot_ab == gr0(mvab_sym));
        CHECK(wdg_ab == gr2(mvab));
        CHECK(wdg_ab == gr2(mvab_asym));
    }

    TEST_CASE("MVec3dp: geometric product tests - bivec * vec")
    {
        fmt::println("MVec3dp: geometric product tests - bivec * vec");

        // Ab = (b << A) + wdg(A,b) = gr1(Ab) + gr3(Ab)
        //
        // (b << A) = 0.5*(Ab - Aa)   (antisymmetric part)
        // wdg(A,b) = 0.5*(Ab + Aa)   (symmetric part)

        bivec3dp A{1.0, 2.0, 3.0, 10.0, 20.0, 30.0};
        vec3dp b{0.5, 3.0, -2.0, 1};
        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        mvec3dp mva{A};
        mvec3dp mvb{b};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb + mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb - mvb * mva);

        // fmt::println("");
        // fmt::println("   A = {}", A);
        // fmt::println("   b = {}", b);
        // fmt::println("   (b << A) = {}, gr1(A*b) = {}", b_onto_A, gr1(A * b));
        // fmt::println("   wdg(A,b) = {}, gr3(A*b) = {}", wdg_Ab, gr3(A * b));
        // fmt::println("");
        // fmt::println("   mva  = {}", mva);
        // fmt::println("   mvb  = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym  = 0.5*(mva * mvb + mvb * mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva * mvb - mvb * mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));
        // fmt::println("");

        CHECK(b_onto_A == gr1(mvab));
        CHECK(b_onto_A == gr1(mvab_asym));
        CHECK(wdg_Ab == gr3(mvab));
        CHECK(wdg_Ab == gr3(mvab_sym));
    }

    TEST_CASE("MVec3dp: geometric product tests - vec * bivec")
    {
        fmt::println("MVec3dp: geometric product tests - vec * bivec");

        // a*B = (B >> a) + wdg(a,B) = gr1(aB) + gr3(aB)
        //
        // (B >> a) = 0.5*(aB - Ba)   (antisymmetric part)
        // wdg(a,B) = 0.5*(aB + Ba)   (symmetric part)

        vec3dp a{0.5, 3.0, -2.0, 1};
        bivec3dp B{1.0, 2.0, 3.0, 10.0, 20.0, 30.0};
        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec3dp mva{a};
        mvec3dp mvb{B};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb + mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb - mvb * mva);

        // fmt::println("");
        // fmt::println("   a = {}", a);
        // fmt::println("   B = {}", B);
        // fmt::println("   (B >> a) = {}", B_by_a);
        // fmt::println("   wdg(a,b) = {}", wdg_aB);
        // fmt::println("");
        // fmt::println("   mva  = {}", mva);
        // fmt::println("   mvb  = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym  = 0.5*(mva * mvb + mvb * mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva * mvb - mvb * mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));
        // fmt::println("");

        CHECK(B_by_a == gr1(mvab));
        CHECK(B_by_a == gr1(mvab_asym));
        CHECK(wdg_aB == gr3(mvab));
        CHECK(wdg_aB == gr3(mvab_sym));
    }

    TEST_CASE("MVec3dp: geometric product tests - trivec * vec")
    {
        fmt::println("MVec3dp: geometric product tests - trivec * vec");

        // old: Ab = dot(A,b) + wdg(A,b) = gr2(Ab) + gr4(Ab)
        // new: Ab = (A >> b) + wdg(A,b) = gr2(Ab) + gr4(Ab)
        //
        // dot(A,b) = 0.5*(Ab + Aa)   (symmetric part)
        // wdg(A,b) = 0.5*(Ab - Aa)   (antisymmetric part)

        trivec3dp A{1.0, 2.0, 3.0, 10.0};
        vec3dp b{0.5, 3.0, -2.0, 1};
        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        mvec3dp mva{A};
        mvec3dp mvb{b};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb - mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb + mvb * mva);

        // fmt::println("   A = {}", A);
        // fmt::println("   b = {}", b);
        // fmt::println("   b_onto_A = {}, gr2(A*b) = {}", b_onto_A, gr2(A * b));
        // fmt::println("   wdg(A,b) = {}, gr4(A*b) = {}", wdg_Ab, gr4(A * b));
        // fmt::println("");
        // fmt::println("   mva  = {}", mva);
        // fmt::println("   mvb  = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym  = 0.5*(mva * mvb + mvb * mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva * mvb - mvb * mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));
        // fmt::println("   gr4(mvab) = {}", gr4(mvab));

        CHECK(b_onto_A == gr2(mvab));
        CHECK(b_onto_A == gr2(mvab_asym));
        CHECK(wdg_Ab == gr4(mvab));
        CHECK(wdg_Ab == gr4(mvab_sym));
    }

    TEST_CASE("MVec3dp: geometric product - link to inner and outer products")
    {
        fmt::println("MVec3dp: geometric product - link to inner and outer products");
        // vec3d u{1.0, 2.0, 3.0};    // 3d vector from ega3d case
        vec3dp u{1.0, 2.0, 3.0, 1.0}; // 3d point in pga3dp

        scalar3dp s{3};
        // vec3d v{-3.0, 2.5, -0.5};    // 3d vector from ega3d case
        vec3dp v{-3.0, 2.5, -0.5, 1.0}; // 3d point in pga3dp
        // bivec3dp B{2.5, 3.5, 1.5};
        bivec3dp B{2.5, 3.5, 1.5, -2.5, -3.5, -1.5};
        trivec3dp t{-7.0, 2.0, -1.0, 2.0};
        pscalar3dp ps{4.0};

        // fmt::println("");
        // fmt::println("u = {}", u);
        // fmt::println("s = {}", s);
        // fmt::println("v = {}", v);
        // fmt::println("B = {}", B);
        // fmt::println("t = {}", t);
        // fmt::println("ps = {}", ps);
        // fmt::println("");

        // fmt::println("");
        // fmt::println("scalar case: ");
        CHECK(u * s == gr1((s >> u) + wdg(u, s)));
        CHECK(u * s == gr1(rwdg(s, r_bulk_dual(u)) + wdg(u, s)));

        // fmt::println("");
        // fmt::println("u * s = {}", u * s);
        // fmt::println("");
        // fmt::println("(s >> u) = {}", (s >> u));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("(s >> u) + wdg(u, s) = {}", (s >> u) + wdg(u, s));
        // fmt::println("gr1((s >> u) + wdg(u, s)) = {}", gr1((s >> u) + wdg(u, s)));
        // fmt::println("");
        // fmt::println("r_bulk_dual(u) = {}", r_bulk_dual(u));
        // fmt::println("rwdg(s,r_bulk_dual(u))= {}", rwdg(s, r_bulk_dual(u)));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("rwdg(s,r_bulk_dual(u)) + wdg(u, s) = {}",
        //              rwdg(s, r_bulk_dual(u)) + wdg(u, s));
        // fmt::println("gr1(rwdg(s,r_bulk_dual(u)) + wdg(u, s)) = {}",
        //              gr1(rwdg(s, r_bulk_dual(u)) + wdg(u, s)));
        // fmt::println("");

        CHECK(s * u == gr1((u << s) + wdg(s, u)));
        CHECK(s * u == gr1(rwdg(l_bulk_dual(u), s) + wdg(s, u)));

        // fmt::println("");
        // fmt::println("s * u = {}", s * u);
        // fmt::println("");
        // fmt::println("(u << s) = {}", (u << s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("(u << s) + wdg(s, u) = {}", (u << s) + wdg(s, u));
        // fmt::println("gr1((u << s) + wdg(s, u)) = {}", gr1((u << s) + wdg(s, u)));
        // fmt::println("");
        // fmt::println("l_bulk_dual(u) = {}", l_bulk_dual(u));
        // fmt::println("rwdg(l_bulk_dual(u), s)= {}", rwdg(l_bulk_dual(u), s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("rwdg(l_bulk_dual(u), s) + wdg(s, u) = {}",
        //              rwdg(l_bulk_dual(u), s) + wdg(s, u));
        // fmt::println("gr1(rwdg(l_bulk_dual(u), s) + wdg(s, u)) = {}",
        //              gr1(rwdg(l_bulk_dual(u), s) + wdg(s, u)));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("vector case: ");
        CHECK(u * v == (v >> u) + wdg(u, v));
        CHECK(u * v == rwdg(v, r_bulk_dual(u)) + wdg(u, v));

        // fmt::println("");
        // fmt::println("u * v = {}", u * v);
        // fmt::println("");
        // fmt::println("(v >> u) = {}", (v >> u));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("(v >> u) + wdg(u, v) = {}", (v >> u) + wdg(u, v));
        // fmt::println("");
        // fmt::println("r_bulk_dual(u) = {}", r_bulk_dual(u));
        // fmt::println("rwdg(v,r_bulk_dual(u))= {}", rwdg(v, r_bulk_dual(u)));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("rwdg(v,r_bulk_dual(u)) + wdg(u, v) = {}",
        //              rwdg(v, r_bulk_dual(u)) + wdg(u, v));
        // fmt::println("");

        CHECK(v * u == (u << v) + wdg(v, u));
        CHECK(v * u == rwdg(l_bulk_dual(u), v) + wdg(v, u));

        // fmt::println("");
        // fmt::println("v * u = {}", v * u);
        // fmt::println("");
        // fmt::println("(u << v) = {}", (u << v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("(u << v) + wdg(v, u) = {}", (u << v) + wdg(v, u));
        // fmt::println("");
        // fmt::println("l_bulk_dual(u) = {}", l_bulk_dual(u));
        // fmt::println("rwdg(l_bulk_dual(u), v)= {}", rwdg(l_bulk_dual(u), v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("rwdg(l_bulk_dual(u), v) + wdg(v, u) = {}",
        //              rwdg(l_bulk_dual(u), v) + wdg(v, u));
        // fmt::println("");


        // fmt::println("");
        // fmt::println("bivector case: ");
        CHECK(u * B == (B >> u) + wdg(u, B));
        CHECK(u * B == rwdg(B, r_bulk_dual(u)) + wdg(u, B));

        // fmt::println("");
        // fmt::println("u * B = {}", u * B);
        // fmt::println("");
        // fmt::println("(B >> u) = {}", (B >> u));
        // fmt::println("wdg(u, B) = {}", wdg(u, B));
        // fmt::println("(B >> u) + wdg(u, B) = {}", (B >> u) + wdg(u, B));
        // fmt::println("");
        // fmt::println("r_bulk_dual(u) = {}", r_bulk_dual(u));
        // fmt::println("rwdg(B,r_bulk_dual(u))= {}", rwdg(B, r_bulk_dual(u)));
        // fmt::println("wdg(u, B) = {}", wdg(u, B));
        // fmt::println("rwdg(B,r_bulk_dual(u)) + wdg(u, B) = {}",
        //              rwdg(B, r_bulk_dual(u)) + wdg(u, B));
        // fmt::println("");

        CHECK(B * u == (u << B) + wdg(B, u));
        CHECK(B * u == rwdg(l_bulk_dual(u), B) + wdg(B, u));

        // fmt::println("");
        // fmt::println("B * u = {}", B * u);
        // fmt::println("");
        // fmt::println("(u << B) = {}", (u << B));
        // fmt::println("wdg(B, u) = {}", wdg(B, u));
        // fmt::println("(u << B) + wdg(B, u) = {}", (u << B) + wdg(B, u));
        // fmt::println("");
        // fmt::println("l_bulk_dual(u) = {}", l_bulk_dual(u));
        // fmt::println("rwdg(l_bulk_dual(u), B)= {}", rwdg(l_bulk_dual(u), B));
        // fmt::println("wdg(B, u) = {}", wdg(B, u));
        // fmt::println("rwdg(l_bulk_dual(u), B) + wdg(B, u) = {}",
        //              rwdg(l_bulk_dual(u), B) + wdg(B, u));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("trivector case: ");
        CHECK(u * t == (t >> u) + wdg(u, t));
        CHECK(u * t == rwdg(t, r_bulk_dual(u)) + wdg(u, t));

        // fmt::println("");
        // fmt::println("u * t = {}", u * t);
        // fmt::println("");
        // fmt::println("(t >> u) = {}", (t >> u));
        // fmt::println("wdg(u, t) = {}", wdg(u, t));
        // fmt::println("(t >> u) + wdg(u, t) = {}", (t >> u) + wdg(u, t));
        // fmt::println("");
        // fmt::println("r_bulk_dual(u) = {}", r_bulk_dual(u));
        // fmt::println("rwdg(t,r_bulk_dual(u))= {}", rwdg(t, r_bulk_dual(u)));
        // fmt::println("wdg(u, t) = {}", wdg(u, t));
        // fmt::println("rwdg(t,r_bulk_dual(u)) + wdg(u, t) = {}",
        //              rwdg(t, r_bulk_dual(u)) + wdg(u, t));
        // fmt::println("");

        CHECK(t * u == (u << t) + wdg(t, u));
        CHECK(t * u == rwdg(l_bulk_dual(u), t) + wdg(t, u));

        // fmt::println("");
        // fmt::println("t * u = {}", t * u);
        // fmt::println("");
        // fmt::println("(u << t) = {}", (u << t));
        // fmt::println("wdg(t, u) = {}", wdg(t, u));
        // fmt::println("(u << t) + wdg(t, u) = {}", (u << t) + wdg(t, u));
        // fmt::println("");
        // fmt::println("l_bulk_dual(u) = {}", l_bulk_dual(u));
        // fmt::println("rwdg(l_bulk_dual(u), t)= {}", rwdg(l_bulk_dual(u), t));
        // fmt::println("wdg(t, u) = {}", wdg(t, u));
        // fmt::println("rwdg(l_bulk_dual(u), t) + wdg(t, u) = {}",
        //              rwdg(l_bulk_dual(u), t) + wdg(t, u));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("pscalar case: ");
        CHECK(u * ps == gr3((ps >> u) + wdg(u, ps)));
        CHECK(u * ps == gr3(rwdg(ps, r_bulk_dual(u)) + wdg(u, ps)));

        // fmt::println("");
        // fmt::println("u * ps = {}", u * ps);
        // fmt::println("");
        // fmt::println("(ps >> u) = {}", (ps >> u));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("(ps >> u) + wdg(u, ps) = {}", (ps >> u) + wdg(u, ps));
        // fmt::println("gr3((ps >> u) + wdg(u, ps)) = {}", gr3((ps >> u) + wdg(u, ps)));
        // fmt::println("");
        // fmt::println("r_bulk_dual(u) = {}", r_bulk_dual(u));
        // fmt::println("rwdg(ps,r_bulk_dual(u))= {}", rwdg(ps, r_bulk_dual(u)));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("rwdg(ps,r_bulk_dual(u)) + wdg(u, ps) = {}",
        //              rwdg(ps, r_bulk_dual(u)) + wdg(u, ps));
        // fmt::println("gr3(rwdg(ps,r_bulk_dual(u)) + wdg(u, ps)) = {}",
        //              gr3(rwdg(ps, r_bulk_dual(u)) + wdg(u, ps)));
        // fmt::println("");

        CHECK(ps * u == gr3((u << ps) + wdg(ps, u)));
        CHECK(ps * u == gr3(rwdg(l_bulk_dual(u), ps) + wdg(ps, u)));

        // fmt::println("");
        // fmt::println("ps * u = {}", ps * u);
        // fmt::println("");
        // fmt::println("(u << ps) = {}", (u << ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("(u << ps) + wdg(ps, u) = {}", (u << ps) + wdg(ps, u));
        // fmt::println("gr3((u << ps) + wdg(ps, u)) = {}", gr3((u << ps) + wdg(ps, u)));
        // fmt::println("");
        // fmt::println("l_bulk_dual(u) = {}", l_bulk_dual(u));
        // fmt::println("rwdg(l_bulk_dual(u), ps)= {}", rwdg(l_bulk_dual(u), ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("rwdg(l_bulk_dual(u), ps) + wdg(ps, u) = {}",
        //              rwdg(l_bulk_dual(u), ps) + wdg(ps, u));
        // fmt::println("gr3(rwdg(l_bulk_dual(u), ps) + wdg(ps, u)) = {}",
        //              gr3(rwdg(l_bulk_dual(u), ps) + wdg(ps, u)));
        // fmt::println("");

        // fmt::println("");
    }

    TEST_CASE("MVec3dp: geometric product tests - equivalence tests")
    {
        fmt::println("MVec3dp: geometric product tests - equivalence tests");

        vec3dp a{1.0, 2.0, 3.0, 1.0};
        vec3dp b{0.5, 3.0, -4.0, 1.0};
        mvec3dp mva{a};
        mvec3dp mvb{b};

        bivec3dp A{1.0, 2.0, 3.0, 10.0, 20.0, 30.0};
        bivec3dp B{0.5, 3.0, -4.0, 5.0, 30.0, -40.0};
        mvec3dp mvA{A};
        mvec3dp mvB{B};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec3dp_e ab = a * b;
        mvec3dp abm = mva * mvb;
        mvec3dp abd(dot_ab, wdg_ab, pscalar3dp{});

        mvec3dp_u Ab = A * b;
        mvec3dp Abm = mvA * mvb;
        mvec3dp Abd{b_onto_A, wdg_Ab};

        mvec3dp_u aB = a * B;
        mvec3dp aBm = mva * mvB;
        mvec3dp aBd{B_by_a, wdg_aB};

        // fmt::println("");
        // fmt::println("   a                                 = {}", a);
        // fmt::println("   mva                               = {}", mva);
        // fmt::println("   b                                 = {}", b);
        // fmt::println("   mvb                               = {}", mvb);
        // fmt::println("   ab  = mvec3dp_e(a * b)            = {}", ab);
        // fmt::println("   abm = mva * mvb                   = {}", abm);
        // fmt::println("   abd = mvec3dp(dot(a,b), wdg(a,b)) = {}", abd);
        // fmt::println("");
        // fmt::println("   A                                 = {}", A);
        // fmt::println("   mvA                               = {}", mvA);
        // fmt::println("   b                                 = {}", b);
        // fmt::println("   mvb                               = {}", mvb);
        // fmt::println("   Ab  = mvec3dp_u(A * b)            = {}", Ab);
        // fmt::println("   Abm = mvA * mvb                   = {}", Abm);
        // fmt::println("   Abd = mvec3dp((b << A), wdg(A,b)) = {}", Abd);
        // fmt::println("");
        // fmt::println("   a                                 = {}", a);
        // fmt::println("   mva                               = {}", mva);
        // fmt::println("   B                                 = {}", B);
        // fmt::println("   mvB                               = {}", mvB);
        // fmt::println("   aB  = mvec3dp_u(a * B)            = {}", aB);
        // fmt::println("   aBm = mva * mvB                   = {}", aBm);
        // fmt::println("   aBd = mvec3dp((B >> a), wdg(a,B)) = {}", aBd);
        // fmt::println("");

        CHECK(gr0(ab) == gr0(abm));
        CHECK(gr1(abm) == vec3dp{});
        CHECK(gr2(ab) == gr2(abm));
        CHECK(gr3(abm) == trivec3dp{});
        CHECK(gr4(abm) == pscalar3dp{});

        CHECK(gr0(ab) == gr0(abd));
        CHECK(gr1(abd) == vec3dp{});
        CHECK(gr2(ab) == gr2(abd));
        CHECK(gr3(abd) == trivec3dp{});
        CHECK(gr4(ab) == pscalar3dp{});

        CHECK(gr0(Abm) == 0);
        CHECK(unitize(gr1(Ab)) == unitize(gr1(Abm)));
        CHECK(gr2(Abm) == bivec3dp{});
        CHECK(gr3(Ab) == gr3(Abm));
        CHECK(gr4(Abm) == pscalar3dp{});

        CHECK(gr0(Abd) == 0);
        CHECK(unitize(gr1(Ab)) == unitize(gr1(Abd)));
        CHECK(gr2(Abd) == bivec3dp{});
        CHECK(gr3(Ab) == gr3(Abd));
        CHECK(gr4(Abd) == pscalar3dp{});

        CHECK(gr0(aBm) == 0);
        CHECK(unitize(gr1(aB)) == unitize(gr1(aBm)));
        CHECK(gr2(aBm) == bivec3dp{});
        CHECK(gr3(aB) == gr3(aBm));
        CHECK(gr4(aBm) == pscalar3dp{});

        CHECK(gr0(aBd) == 0);
        CHECK(unitize(gr1(aB)) == unitize(gr1(aBd)));
        CHECK(gr2(aBd) == bivec3dp{});
        CHECK(gr3(aB) == gr3(aBd));
        CHECK(gr4(aBd) == pscalar3dp{});
    }

    TEST_CASE("MVec3dp: assignment tests")
    {
        fmt::println("MVec3dp: assignment tests");

        scalar3dp s1{5.0};
        vec3dp v1{1.0, 2.0, 3.0, 23.0};
        vec3dp v2{0.5, 1.0, 1.5, 11.5};
        vec3dp v3{0.5, 1.0, -4.5, 1.0};
        bivec3dp b1{1.0, 2.0, 3.0, 10.0, 20.0, 30.0};
        trivec3dp t1{1.5, 3.0, 4.5, 7.0};
        pscalar3dp ps1{6.0};

        mvec3dp mv1{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0,
                    0.0, 0.5, 1.0, 1.5, 11.5, 15.5, 6.0,  61.5};
        mvec3dp mv2{0.0, 0.5, 1.0, 1.5, 11.5, 15.5, 6.0,  61.5,
                    0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        mvec3dp mv3{mv1};
        mvec3dp mv4 = mv2;

        mvec3dp mv5(s1);
        mvec3dp mv6(ps1);
        mvec3dp mv7{v1};
        mvec3dp mv8{b1};
        mvec3dp mv9{dot(v1, v3), wdg(v1, v3), ps1};

        mvec3dp mv10{v1, t1};

        // fmt::println("");
        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   mv3 = {}", mv3);
        // fmt::println("   mv4 = {}", mv4);
        // fmt::println("   mv5 = {}", mv5);
        // fmt::println("   mv6 = {}", mv6);
        // fmt::println("");
        // fmt::println("   gr1(mv1) = {}", gr1(mv1));
        // fmt::println("   gr1(mv2) = {}", gr1(mv2));
        // fmt::println("   gr1(mv3) = {}", gr1(mv3));
        // fmt::println("   gr1(mv4) = {}", gr1(mv4));
        // fmt::println("");
        // fmt::println("   v1       = {}", v1);
        // fmt::println("   mv7 = v1 = {}", mv7);
        // fmt::println("   b1       = {}", b1);
        // fmt::println("   mv8 = b1 = {}", mv8);
        // fmt::println("");
        // fmt::println("   mv9  = {}", mv9);
        // fmt::println("   mv10 = {}", mv10);
        // fmt::println("");

        CHECK(gr1(mv1) == v1);
        CHECK(gr1(mv2) == v2);
        CHECK(gr1(mv3) == v1);
        CHECK(gr1(mv4) == v2);
        CHECK(gr0(mv5) == 5.0);
        CHECK(gr4(mv6) == 6.0);
        CHECK(mv1 == mv3);
        CHECK(mv4 == mv2);
        CHECK(gr1(mv7) == v1);
        CHECK(gr2(mv8) == b1);
        CHECK(gr0(mv9) == dot(v1, v3));
        CHECK(gr2(mv9) == wdg(v1, v3));
        CHECK(gr4(mv9) == ps1);
        CHECK(gr1(mv10) == v1);
        CHECK(gr3(mv10) == t1);
    }

    TEST_CASE("MVec3dp: bivector product properties")
    {
        fmt::println("MVec3dp: bivector product properties");

        bivec3dp b1{1.0, 2.0, 3.0, 10.0, 20.0, 30.0};
        mvec3dp mb1{b1};
        bivec3dp b2{-3.0, 1.0, 2.0, -30.0, 10.0, 20.0};
        mvec3dp mb2{b2};

        auto gpr12_m = mb1 * mb2;
        auto gpr21_m = mb2 * mb1;
        auto gpr12_m_sym = 0.5 * (gpr12_m + gpr21_m);
        auto gpr12_m_asym = 0.5 * (gpr12_m - gpr21_m);

        auto gpr12_d = b1 * b2;
        auto gpr21_d = b2 * b1;
        auto gpr12_d_sym = 0.5 * (gpr12_d + gpr21_d);
        auto gpr12_d_asym = 0.5 * (gpr12_d - gpr21_d);

        // fmt::println("");
        // fmt::println("   b1  = {}", b1);
        // fmt::println("   mb1 = {}", mb1);
        // fmt::println("   b2  = {}", b2);
        // fmt::println("   mb2 = {}", mb2);
        // fmt::println("");
        // fmt::println("   dot(b1, b2)         = {}", dot(b1, b2));
        // fmt::println("   cmt(b1, b2)         = {}", cmt(b1, b2));
        // fmt::println("   wdg(b1, b2)         = {}", wdg(b1, b2));
        // fmt::println("   gpr12_m = mb1 * mb2 = {}", gpr12_m);
        // fmt::println("");
        // fmt::println("   dot(b2, b1)         = {}", dot(b2, b1));
        // fmt::println("   cmt(b2, b1)         = {}", cmt(b2, b1));
        // fmt::println("   wdg(b2, b1)         = {}", wdg(b2, b1));
        // fmt::println("   gpr21_m = mb2 * mb1 = {}", gpr21_m);
        // fmt::println("");
        // fmt::println("   gpr12_m_sym  = 0.5*(gpr12_d + gpr21_d) = {}", gpr12_m_sym);
        // fmt::println("   gpr12_m_asym = 0.5*(gpr12_m - gpr21_m) = {}", gpr12_m_asym);
        // fmt::println("");
        // fmt::println("   gpr12_d = b1 * b2                      = {}", gpr12_d);
        // fmt::println("   gpr21_d = b2 * b1)                     = {}", gpr21_d);
        // fmt::println("   gpr12_d_sym  = 0.5*(gpr12_d + gpr21_d) = {}", gpr12_d_sym);
        // fmt::println("   gpr12_d_asym = 0.5*(gpr12_d - gpr21_d) = {}", gpr12_d_asym);
        // fmt::println("");

        CHECK(gr2(mb1) == b1);

        // just to silence unused variable warnings
        CHECK(gpr12_m_sym == 0.5 * (gpr12_m + gpr21_m));
        CHECK(gpr12_m_asym == 0.5 * (gpr12_m - gpr21_m));
        CHECK(gpr12_d_sym == 0.5 * (gpr12_d + gpr21_d));
        CHECK(gpr12_d_asym == 0.5 * (gpr12_d - gpr21_d));
    }


    ////////////////////////////////////////////////////////////////////////////////
    // MVec3dp_E<T> and MVec3dp_U<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3dp_E: product characteristics")
    {
        fmt::println("MVec3dp_E: product characteristics");

        auto s = scalar3dp{5.0};
        auto v = vec3dp{1.0, 2.0, 3.0, 1.0};
        auto B = bivec3dp{-1.0, 2.0, 1.0, -10.0, 20.0, 10.0};
        auto t = trivec3dp{3.0, 6.0, 9.0, 3.0};
        auto ps = pscalar3dp{-5.0};

        auto mv = mvec3dp{s, v, B, t, ps};
        auto mv_e = mvec3dp_e{s, B, ps};
        auto mv_u = mvec3dp{v, t};

        auto s2 = scalar3dp{2.5};
        // auto v2 = vec3dp{0.5, 1.0, 1.5, 2.0};
        auto B2 = bivec3dp{5.0, 10.0, 15.0, 50.0, 100.0, 150.0};
        // auto t2 = trivec3dp{1.5, 3.0, 4.5, 6.0};
        auto ps2 = pscalar3dp{-1.5};

        // scalar and bivector commute
        CHECK(s * B == B * s);

        // pscalar and bivector commute
        CHECK(ps * B == B * ps);

        // vector and bivector anti-commute
        CHECK(v * B == -rev(B * v));
        CHECK(rev(v * B) == rev(B) * rev(v));

        // scalar and pscalar parts commute while bivector anticommutes
        CHECK(B * B2 == rev(B2 * B));
        CHECK(rev(B * B2) == rev(B2) * rev(B));

        auto mv_e1 = mvec3dp_e{s, B, ps};
        auto mv_e2 = mvec3dp_e{s2, B2, ps2};
        auto mv_e3 = get_motor(
            z_axis_3dp, deg2rad(15)); // only scalar and bivector part (as in a rotor)

        CHECK(rev(mv_e1 * mv_e2) == rev(mv_e2) * rev(mv_e1));

        CHECK(rev(mv_e1 * B2) == rev(B2) * rev(mv_e1));

        // an even-grade multivector and and arbitrary bivector do not commute
        CHECK(wdg(B2, B2) != 0); // blade test, B2 is not a blade
        // CHECK(wdg(mv_e1, mv_e1) == mvec3dp_e{}); // blade test
        // fmt::println("mv_e1             = {}", mv_e1);
        // fmt::println("wdg(mv_e1, mv_e1) = {}", wdg(mv_e1, mv_e1));
        CHECK(mv_e1 * B2 != B2 * mv_e1);

        // a bivector and its motor commute
        CHECK(wdg(z_axis_3dp, z_axis_3dp) == 0); // blade test, z_axis_3dp is a blade
        // CHECK(wdg(mv_e3, mv_e3) == mvec3dp_e{}); // blade test
        // fmt::println("mv_e3             = {}", mv_e3);
        // fmt::println("wdg(mv_e3, mv_e3) = {}", wdg(mv_e3, mv_e3));
        CHECK(z_axis_3dp * mv_e3 == mv_e3 * z_axis_3dp);
        CHECK(rgpr(z_axis_3dp, mv_e3) == rgpr(mv_e3, z_axis_3dp));
    }

    TEST_CASE("MVec3dp_E/_U: reflections and motors (rotations, translations)")
    {
        fmt::println("MVec3dp_E/_U: reflections and motors (rotations, translations)");

        {
            /////////////////////////////////////////////////////////////////////////////
            // 1st: start with points use planes through the origin and reflect
            // points
            /////////////////////////////////////////////////////////////////////////////

            // define points and planes
            auto p0 = O_3dp;
            auto p1 = vec3dp{1, 0, 0, 1};
            auto p2 = vec3dp{1, 1, 0, 1};
            auto p3 = vec3dp{0, 0, 1, 1};
            auto p = vec3dp{1, -0.5, 0, 1};

            auto pl1 = unitize(wdg(wdg(p0, p3), p1));
            auto pl2 = unitize(wdg(wdg(p0, p3), p2));

            // reflect p on pl1 -> pr and reflect pr on pl2 -> prr
            auto pr = -gr1(rgpr(rgpr(pl1, p), pl1));

            auto prr = -gr1(rgpr(rgpr(pl2, pr), pl2));

            // fmt::println("");
            // fmt::println("pl1: {}, pl1u: {}, pl2: {}, pl2u: {}", pl1,
            // unitize(pl1), pl2,
            //              unitize(pl2));
            // fmt::println("pr: {}, pru: {}", pr, unitize(pr));
            // fmt::println("prr: {}, prru: {}", prr, unitize(prr));
            // fmt::println("");

            CHECK(unitize(pr) == vec3dp{1, 0.5, 0, 1});
            CHECK(unitize(prr) == vec3dp{0.5, 1, 0, 1});

            // show that prr can be obtained directly from a motor as pm
            auto motor = rgpr(pl2, pl1);

            auto pm = gr1(rgpr(rgpr(motor, p), rrev(motor)));

            // fmt::println("");
            // fmt::println("pm: {}, pmu: {}", pm, unitize(pm));
            // fmt::println("");

            CHECK(unitize(prr) == unitize(pm));

            // TODO: add inverse inv() for an even grade multivector
            //       see paper Hitzer et.al.

            // fmt::println("");
            // fmt::println("motor: {}", motor);
            // fmt::println("rrev(motor): {}", rrev(motor));
            // // fmt::println("inv(motor): {}", inv(motor));
            // fmt::println("rrev(motor): {}", rrev(motor));
            // fmt::println("rgpr(motor, rrev(motor)): {}", rgpr(motor, rrev(motor)));
            // fmt::println("gr4(rgpr(motor, rrev(motor))): {}",
            //              gr4(rgpr(motor, rrev(motor))));
            // fmt::println("");

            // CHECK(abs(double(gr4(rgpr(motor, rrev(motor))) - I_3dp)) < eps);
            // CHECK(rrev(motor) == inv(motor));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 2nd: start with lines through the origin and reflect them in planes
            // through the origin
            /////////////////////////////////////////////////////////////////////////////

            // define points and planes
            auto p0 = O_3dp;
            auto p1 = vec3dp{1, 0, 0, 1};
            auto p2 = vec3dp{1, 1, 0, 1};
            auto p3 = vec3dp{0, 0, 1, 1};
            auto p = vec3dp{1, -0.5, 0, 1};

            auto pl1 = unitize(wdg(wdg(p0, p3), p1));
            auto pl2 = unitize(wdg(wdg(p0, p3), p2));

            auto l = wdg(p0, p);

            // reflect l on pl1 -> lr and reflect lr on pl2 -> lrr
            auto lr = -gr2(rgpr(rgpr(pl1, l), pl1));
            auto lrr = -gr2(rgpr(rgpr(pl2, lr), pl2));

            // fmt::println("");
            // fmt::println("l: {}, lu: {}, att(l): {}", l, unitize(l), att(l));
            // fmt::println("");
            // fmt::println("pl1: {}, pl1u: {}, att(pl1): {}", pl1, unitize(pl1),
            // att(pl1)); fmt::println("pl2: {}, pl2u: {}, att(pl2): {}", pl2,
            // unitize(pl2), att(pl2)); fmt::println(""); fmt::println("lr: {},
            // lru: {}, att(lr): {}", lr, unitize(lr), att(lr));
            // fmt::println("lrr: {}, lrru: {}, att(lrr): {}", lrr, unitize(lrr),
            // att(lrr)); fmt::println("");

            CHECK(unitize(lr) == unitize(join(p0, vec3dp{1, 0.5, 0, 1})));
            CHECK(unitize(lrr) == unitize(join(p0, vec3dp{0.5, 1, 0, 1})));

            // show that prr can be obtained directly from a rotation via a motor
            // as pm
            auto motor = rgpr(pl2, pl1);

            auto lm = gr2(rgpr(rgpr(motor, l), rrev(motor)));

            // fmt::println("");
            // fmt::println("lm: {}, lmu: {}, att(lm): {}", lm, unitize(lm),
            // att(lm)); fmt::println("");

            CHECK(lrr == lm);
            CHECK(unitize(lrr) == unitize(lm));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 3rd: as 1st step (reflection with points), but with the planes for
            // reflection and the center of rotation not in the origin
            /////////////////////////////////////////////////////////////////////////////

            // define points and planes
            // p0 was O_3dp, now shifted to new point
            auto p0 = vec3dp{1, 0.5, 0, 1};
            auto delta = p0 - O_3dp;
            auto p1 = vec3dp{1, 0, 0, 1} + delta;
            auto p2 = vec3dp{1, 1, 0, 1} + delta;
            auto p3 = vec3dp{0, 0, 1, 1} + delta;
            auto p = vec3dp{1, -0.5, 0, 1} + delta;

            auto pl1 = unitize(wdg(wdg(p0, p3), p1));
            auto pl2 = unitize(wdg(wdg(p0, p3), p2));

            // reflect p on pl1 -> pr and reflect pr on pl2 -> prr
            auto pr = -gr1(rgpr(rgpr(pl1, p), pl1));
            auto prr = -gr1(rgpr(rgpr(pl2, pr), pl2));

            // fmt::println("");
            // fmt::println("pl1: {}, pl1u: {}, pl2: {}, pl2u: {}", pl1,
            // unitize(pl1), pl2,
            //              unitize(pl2));
            // fmt::println("pr: {}, pru: {}", pr, unitize(pr));
            // fmt::println("prr: {}, prru: {}", prr, unitize(prr));
            // fmt::println("");

            CHECK(unitize(pr) == vec3dp{1, 0.5, 0, 1} + delta);
            CHECK(unitize(prr) == vec3dp{0.5, 1, 0, 1} + delta);

            // show that prr can be obtained directly from a movement via a motor
            // as pm
            auto motor = rgpr(pl2, pl1);
            auto pm = gr1(rgpr(rgpr(motor, p), rrev(motor)));

            // fmt::println("");
            // fmt::println("pm: {}, pmu: {}", pm, unitize(pm));
            // fmt::println("");

            CHECK(unitize(prr) == unitize(pm));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 4th: same at the 2nd step (reflecion of planes on planes), but with
            // planes for reflection and the center of rotation not in the origin
            /////////////////////////////////////////////////////////////////////////////

            // define points and planes
            // p0 was O_3dp, now shifted to new point
            auto p0 = vec3dp{1, 0.5, 0, 1};
            auto delta = p0 - O_3dp;
            auto p1 = vec3dp{1, 0, 0, 1} + delta;
            auto p2 = vec3dp{1, 1, 0, 1} + delta;
            auto p3 = vec3dp{0, 0, 1, 1} + delta;
            auto p = vec3dp{1, -0.5, 0, 1} + delta;

            auto pl1 = unitize(wdg(wdg(p0, p3), p1));
            auto pl2 = unitize(wdg(wdg(p0, p3), p2));

            auto l = wdg(p0, p);

            // reflect l on pl1 -> lr and reflect lr on pl2 -> lrr
            auto lr = -gr2(rgpr(rgpr(pl1, l), pl1));
            auto lrr = -gr2(rgpr(rgpr(pl2, lr), pl2));

            // fmt::println("");
            // fmt::println("l: {}, lu: {}, att(l): {}", l, unitize(l), att(l));
            // fmt::println("");
            // fmt::println("pl1: {}, pl1u: {}, att(pl1): {}", pl1, unitize(pl1),
            // att(pl1)); fmt::println("pl2: {}, pl2u: {}, att(pl2): {}", pl2,
            // unitize(pl2), att(pl2)); fmt::println(""); fmt::println("lr: {},
            // lru: {},"
            //              "att(lr): {}",
            //              lr, unitize(lr), att(lr));
            // fmt::println("lrr: {}, lrru: {}, "
            //              "att(lrr): {}",
            //              lrr, unitize(lrr), att(lrr));
            // fmt::println("");

            CHECK(unitize(lr) == unitize(join(p0, vec3dp{1, 0.5, 0, 1} + delta)));
            CHECK(unitize(lrr) == unitize(join(p0, vec3dp{0.5, 1, 0, 1} + delta)));

            // show that prr can be obtained directly from a rotation via a motor
            // as pm
            auto motor = rgpr(pl2, pl1);

            auto lm = gr2(rgpr(rgpr(motor, l), rrev(motor)));

            // fmt::println("");
            // fmt::println("lm: {}, lmu: {}, att(lm): {}", lm, unitize(lm),
            // att(lm)); fmt::println("");

            CHECK(lrr == lm);
            CHECK(unitize(lrr) == unitize(lm));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 5th step: do the same as above, but use the optimized motor
            // formulas
            /////////////////////////////////////////////////////////////////////////////

            // define points and planes
            // p0 was O_3dp, now shifted to that new point
            auto p0 = vec3dp{1, 0.5, 0, 1};
            auto delta = p0 - O_3dp;
            auto p1 = vec3dp{1, 0, 0, 1} + delta;
            auto p2 = vec3dp{1, 1, 0, 1} + delta;
            auto p3 = vec3dp{0, 0, 1, 1} + delta;
            auto p = vec3dp{1, -0.5, 0, 1} + delta;

            auto pl1 = unitize(wdg(wdg(p0, p3), p1));
            auto pl2 = unitize(wdg(wdg(p0, p3), p2));

            auto l = wdg(p0, p);

            auto R = get_motor_from_planes(pl1, pl2);
            CHECK(R == unitize(rgpr(pl2, pl1)));

            auto pm_manual = gr1(rgpr(rgpr(R, p), rrev(R)));
            auto pm_orig = move3dp(p, R);
            auto pm = move3dp(p, R);

            CHECK(pm_manual == pm_orig);
            CHECK(pm_manual == pm);

            auto lm_manual = gr2(rgpr(rgpr(R, l), rrev(R)));
            auto lm_orig = move3dp(l, R);
            auto lm = move3dp(l, R);

            // fmt::println("");
            // fmt::println("lm_orig: {}, lmu_orig: {}, att(lm_orig): {}",
            // lm_orig,
            //              unitize(lm_orig), att(lm_orig));
            // fmt::println("");
            // fmt::println("");
            // fmt::println("lm: {}, lmu: {}, att(lm): {}", lm, unitize(lm),
            // att(lm)); fmt::println("");

            CHECK(lm_manual == lm_orig);
            CHECK(lm_manual == lm);
        }

        //     {
        //         /////////////////////////////////////////////////////////////////////////////
        //         // 6th step: optional speed testing for the optimized motor formulas
        //         // -> advantage for debug mode, but not when enabling optimization
        //         /////////////////////////////////////////////////////////////////////////////

        //         // define points and planes
        //         // p0 was O_3dp, now shifted to that new position
        //         auto p0 = vec3dp{1, 0.5, 0, 1};
        //         auto delta = p0 - O_3dp;
        //         auto p1 = vec3dp{1, 0, 0, 1} + delta;
        //         auto p2 = vec3dp{1, 1, 0, 1} + delta;
        //         auto p3 = vec3dp{0, 0, 1, 1} + delta;
        //         auto p = vec3dp{1, -0.5, 0, 1} + delta;

        //         auto pl1 = unitize(wdg(wdg(p0, p3), p1));
        //         auto pl2 = unitize(wdg(wdg(p0, p3), p2));

        //         auto l = wdg(p0, p);
        //         auto pl = wdg(wdg(p0, p), p3);

        //         auto R = get_motor_from_planes(pl1, pl2);
        //         CHECK(R == unitize(rgpr(pl2, pl1)));

        //         // checking time required
        //         constexpr size_t steps = 10'000'000;
        //         vec3dp pm_sum_orig{};
        //         auto start = std::chrono::system_clock::now();
        //         for (size_t i = 0; i < steps; ++i) {
        //             auto pm_orig = move3dp(p, R);
        //             pm_sum_orig += pm_orig; // just to avoid replacement with opt
        //         }
        //         auto end = std::chrono::system_clock::now();
        //         auto elapsed =
        //             std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        //         fmt::println("");
        //         fmt::println("The measurement orig for point took {}", elapsed);
        //         fmt::println("pm_sum_orig = {}", unitize(pm_sum_orig * (1.0 / steps)));

        //         vec3dp pm_sum{};
        //         start = std::chrono::system_clock::now();
        //         for (size_t i = 0; i < steps; ++i) {
        //             auto pm = move3dp_opt(p, R);
        //             pm_sum += pm; // just to avoid full replacement with opt
        //         }
        //         end = std::chrono::system_clock::now();
        //         elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //         start); fmt::println("The measurement opt for point took {}", elapsed);
        //         fmt::println("pm_sum = {}", unitize(pm_sum * (1.0 / steps)));
        //         fmt::println("");

        //         bivec3dp lm_sum_orig{};
        //         start = std::chrono::system_clock::now();
        //         for (size_t i = 0; i < steps; ++i) {
        //             auto lm_orig = move3dp(l, R);
        //             lm_sum_orig += lm_orig; // just to avoid replacement with opt
        //         }
        //         end = std::chrono::system_clock::now();
        //         elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //         start); fmt::println("The measurement orig for line took {}", elapsed);
        //         fmt::println("lm_sum_orig = {}", lm_sum_orig * (1.0 / steps));

        //         bivec3dp lm_sum{};
        //         start = std::chrono::system_clock::now();
        //         for (size_t i = 0; i < steps; ++i) {
        //             auto lm = move3dp(l, R);
        //             lm_sum += lm; // just to avoid replacement with opt
        //         }
        //         end = std::chrono::system_clock::now();
        //         elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //         start); fmt::println("The measurement opt for line took {}", elapsed);
        //         fmt::println("lm_sum = {}", lm_sum * (1.0 / steps));
        //         fmt::println("");

        //         trivec3dp tm_sum_orig{};
        //         start = std::chrono::system_clock::now();
        //         for (size_t i = 0; i < steps; ++i) {
        //             auto tm_orig = move3dp(pl, R);
        //             tm_sum_orig += tm_orig; // just to avoid replacement with opt
        //         }
        //         end = std::chrono::system_clock::now();
        //         elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //         start); fmt::println("The measurement orig for line took {}", elapsed);
        //         fmt::println("tm_sum_orig = {}", tm_sum_orig * (1.0 / steps));

        //         trivec3dp tm_sum{};
        //         start = std::chrono::system_clock::now();
        //         for (size_t i = 0; i < steps; ++i) {
        //             auto tm = move3dp(pl, R);
        //             tm_sum += tm; // just to avoid replacement with opt
        //         }
        //         end = std::chrono::system_clock::now();
        //         elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //         start); fmt::println("The measurement opt for line took {}", elapsed);
        //         fmt::println("tm_sum = {}", tm_sum * (1.0 / steps));
        //         fmt::println("");

        //         CHECK(pm_sum_orig * (1.0 / steps) == pm_sum * (1.0 / steps));
        //         CHECK(lm_sum_orig * (1.0 / steps) == lm_sum * (1.0 / steps));
        //         CHECK(tm_sum_orig * (1.0 / steps) == tm_sum * (1.0 / steps));

        //         // summary: manual optimization brings benefit in debug-mode.
        //         //          In release-mode with -O3 optimization there is no speed
        //         delta
        //         //          => manual optimization does NOT bring benefit        }
        //     }


        {
            /////////////////////////////////////////////////////////////////////////////
            // 7th: create the motors directly
            /////////////////////////////////////////////////////////////////////////////

            // define points and planes
            // p0 was O_3dp, now shifted to that new point
            auto p0 = vec3dp{1, 0.5, 0, 1};
            auto delta = p0 - O_3dp;
            auto p1 = vec3dp{1, 0, 0, 1} + delta;
            auto p2 = vec3dp{1, 1, 0, 1} + delta;
            auto p3 = vec3dp{0, 0, 1, 1} + delta;
            auto p = vec3dp{1, -0.5, 0, 1} + delta;
            auto pt = vec3dp{0.5, 1, 0, 1} + delta; // target point after applying motor

            auto pl1 = unitize(wdg(wdg(p0, p3), p1));
            auto pl2 = unitize(wdg(wdg(p0, p3), p2));

            auto l = wdg(p0, p);

            auto lis = unitize(rwdg(pl2, pl1)); // intersection line of planes
                                                // consistent orientation of line
                                                // w.r.t to rgpr used for motor only when
                                                // using this sequence, otherwise the
                                                // attitude will be in the other direction

            auto R = get_motor_from_planes(pl1, pl2);
            CHECK(R == unitize(rgpr(pl2, pl1)));

            auto pm = gr1(rgpr(rgpr(R, p), rrev(R)));  // transformation
            auto pb = gr1(rgpr(rgpr(rrev(R), pm), R)); // reverse transformation

            // fmt::println("");
            // fmt::println("pl1:        {:.4g}", pl1);
            // fmt::println("bulk:       {:.4g}", bulk(pl1));
            // fmt::println("weight:     {:.4g}", weight(pl1));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(pl1));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(pl1));
            // fmt::println("");
            // fmt::println("pl2:        {:.4g}", pl2);
            // fmt::println("bulk:       {:.4g}", bulk(pl2));
            // fmt::println("weight:     {:.4g}", weight(pl2));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(pl2));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(pl2));
            // fmt::println("");
            // fmt::println("R:          {:.4g}", R);
            // fmt::println("bulk:       {:.4g}", bulk(R));
            // fmt::println("weight:     {:.4g}", weight(R));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(R));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(R));
            // fmt::println("");
            // fmt::println("lis:        {:.4g}", lis);
            // fmt::println("att(lis):   {:.4g}", att(lis));
            // fmt::println("");
            // fmt::println("pt: {:.4g}, ptu: {:.4g} (target value)", pt, unitize(pt));
            // fmt::println("pm: {:.4g}, pmu: {:.4g} (after trafo)", pm, unitize(pm));
            // fmt::println("pb: {:.4g}, pbu: {:.4g} (after backtrafo)", pb, unitize(pb));
            // fmt::println("");

            CHECK(unitize(pm) == unitize(pt));
            CHECK(unitize(pb) == unitize(p));

            // now build the motor by hand
            auto S = get_motor(lis, pi / 2.); // ATTENTION: sign of turning angle depends
                                              //            on attitude of line!
            auto ps = move3dp(p, S);

            // fmt::println("");
            // fmt::println("S:          {:.4g}", S);
            // fmt::println("bulk:       {:.4g}", bulk(S));
            // fmt::println("weight:     {:.4g}", weight(S));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(S));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(S));
            // fmt::println("");

            CHECK(R == S);
            CHECK(ps == pt);
            CHECK(move3dp(l, S) == wdg(p0, pt));

            auto T = get_motor(delta);
            auto pst = move3dp(p, T);
            // fmt::println("");
            // fmt::println("p: {:.4g}, pu: {:.4g}", p, unitize(p));
            // fmt::println("delta: {:.4g}", delta);
            // fmt::println("T: {:.4g}, Tu: {:.4g}", T, unitize(T));
            // fmt::println("pst: {:.4g}, pstu: {:.4g} (after trafo)", pst, unitize(pst));
            // fmt::println("");

            CHECK(pst == p + delta);
        }
    }


    TEST_CASE("MVec3dp: complement operation")
    {
        fmt::println("MVec3dp: complement operation");

        auto s = scalar3dp{5.0};
        auto v = vec3dp{1.0, 2.0, 3.0, 1.0};
        auto B = bivec3dp{-1.0, 2.0, 1.0, -10.0, 20.0, 10.0};
        auto t = trivec3dp{3.0, 6.0, 9.0, 3.0};
        auto ps = pscalar3dp{-5.0};

        auto mv = mvec3dp{s, v, B, t, ps};
        auto mv_e = mvec3dp_e{s, B, ps};
        auto mv_u = mvec3dp{v, t};

        auto s2 = scalar3dp{2.5};
        auto v2 = vec3dp{0.5, 1.0, 1.5, 2.0};
        auto B2 = bivec3dp{5.0, 10.0, 15.0, 50.0, 100.0, 150.0};
        auto t2 = trivec3dp{1.5, 3.0, 4.5, 6.0};
        auto ps2 = pscalar3dp{-1.5};

        fmt::println("s = {}, r_cmpl(s) = {}", s, r_cmpl(s));
        fmt::println("bulk_nrm_sq(s) = {}, wdg(s, r_cmpl(s)) = {}", bulk_nrm_sq(s),
                     wdg(s, r_cmpl(s)));
        fmt::println("");
        fmt::println("v = {}, r_cmpl(v) = {}", v, r_cmpl(v));
        fmt::println("bulk_nrm_sq(v) = {}", bulk_nrm_sq(v));
        fmt::println("weight_nrm_sq(v) = {}, ", weight_nrm_sq(v));
        fmt::println("geom_nrm(v) = {}", geom_nrm(v));
        fmt::println("geom_nrm_sq(v) = {}, ", geom_nrm_sq(v));
        fmt::println("wdg(v, r_cmpl(v)) = {}", wdg(v, r_cmpl(v)));
        fmt::println("");

        // following checks are only valid for blases, not for multivectors
        CHECK(wdg(s, r_cmpl(s)) == geom_nrm_sq(s) * I_3dp);
        CHECK(wdg(v, r_cmpl(v)) == geom_nrm_sq(v) * I_3dp);
        CHECK(wdg(B, r_cmpl(B)) == geom_nrm_sq(B) * I_3dp);
        CHECK(wdg(t, r_cmpl(t)) == geom_nrm_sq(t) * I_3dp);
        CHECK(wdg(ps, r_cmpl(ps)) == geom_nrm_sq(ps) * I_3dp);

        CHECK(wdg(l_cmpl(s), s) == geom_nrm_sq(s) * I_3dp);
        CHECK(wdg(l_cmpl(v), v) == geom_nrm_sq(v) * I_3dp);
        CHECK(wdg(l_cmpl(B), B) == geom_nrm_sq(B) * I_3dp);
        CHECK(wdg(l_cmpl(t), t) == geom_nrm_sq(t) * I_3dp);
        CHECK(wdg(l_cmpl(ps), ps) == geom_nrm_sq(ps) * I_3dp);

        // complement properties (right and left complements are inverse to each
        // other)
        CHECK(l_cmpl(r_cmpl(s)) == s);
        CHECK(l_cmpl(r_cmpl(v)) == v);
        CHECK(l_cmpl(r_cmpl(B)) == B);
        CHECK(l_cmpl(r_cmpl(t)) == t);
        CHECK(l_cmpl(r_cmpl(ps)) == ps);
        CHECK(l_cmpl(r_cmpl(mv)) == mv);
        CHECK(l_cmpl(r_cmpl(mv_e)) == mv_e);
        CHECK(l_cmpl(r_cmpl(mv_u)) == mv_u);

        CHECK(r_cmpl(l_cmpl(s)) == s);
        CHECK(r_cmpl(l_cmpl(v)) == v);
        CHECK(r_cmpl(l_cmpl(B)) == B);
        CHECK(r_cmpl(l_cmpl(t)) == t);
        CHECK(r_cmpl(l_cmpl(ps)) == ps);
        CHECK(r_cmpl(l_cmpl(mv)) == mv);
        CHECK(r_cmpl(l_cmpl(mv_e)) == mv_e);
        CHECK(r_cmpl(l_cmpl(mv_u)) == mv_u);

        // complement values
        CHECK(r_cmpl(scalar3dp(1.0)) == I_3dp);
        CHECK(r_cmpl(e1_3dp) == e423_3dp);
        CHECK(r_cmpl(e2_3dp) == e431_3dp);
        CHECK(r_cmpl(e3_3dp) == e412_3dp);
        CHECK(r_cmpl(e4_3dp) == e321_3dp);
        CHECK(r_cmpl(e41_3dp) == -e23_3dp);
        CHECK(r_cmpl(e42_3dp) == -e31_3dp);
        CHECK(r_cmpl(e43_3dp) == -e12_3dp);
        CHECK(r_cmpl(e23_3dp) == -e41_3dp);
        CHECK(r_cmpl(e31_3dp) == -e42_3dp);
        CHECK(r_cmpl(e12_3dp) == -e43_3dp);
        CHECK(r_cmpl(e423_3dp) == -e1_3dp);
        CHECK(r_cmpl(e431_3dp) == -e2_3dp);
        CHECK(r_cmpl(e412_3dp) == -e3_3dp);
        CHECK(r_cmpl(e321_3dp) == -e4_3dp);
        CHECK(r_cmpl(I_3dp) == scalar3dp(1.0));

        CHECK(l_cmpl(scalar3dp(1.0)) == I_3dp);
        CHECK(l_cmpl(e1_3dp) == -e423_3dp);
        CHECK(l_cmpl(e2_3dp) == -e431_3dp);
        CHECK(l_cmpl(e3_3dp) == -e412_3dp);
        CHECK(l_cmpl(e4_3dp) == -e321_3dp);
        CHECK(l_cmpl(e41_3dp) == -e23_3dp);
        CHECK(l_cmpl(e42_3dp) == -e31_3dp);
        CHECK(l_cmpl(e43_3dp) == -e12_3dp);
        CHECK(l_cmpl(e23_3dp) == -e41_3dp);
        CHECK(l_cmpl(e31_3dp) == -e42_3dp);
        CHECK(l_cmpl(e12_3dp) == -e43_3dp);
        CHECK(l_cmpl(e423_3dp) == e1_3dp);
        CHECK(l_cmpl(e431_3dp) == e2_3dp);
        CHECK(l_cmpl(e412_3dp) == e3_3dp);
        CHECK(l_cmpl(e321_3dp) == e4_3dp);
        CHECK(l_cmpl(I_3dp) == scalar3dp(1.0));


        // complements are defined for basis elements only
        // => magnitude has to be covered separately for non-normalized elements

        // left complements
        CHECK(wdg(l_cmpl(s), s) / geom_nrm_sq(s) == I_3dp);
        CHECK(wdg(l_cmpl(v), v) / geom_nrm_sq(v) == I_3dp);
        CHECK(wdg(l_cmpl(B), B) / geom_nrm_sq(B) == I_3dp);
        CHECK(wdg(l_cmpl(t), t) / geom_nrm_sq(t) == I_3dp);
        CHECK(wdg(l_cmpl(ps), ps) / geom_nrm_sq(ps) == I_3dp);

        // right complements
        CHECK(wdg(s, r_cmpl(s)) / geom_nrm_sq(s) == I_3dp);
        CHECK(wdg(v, r_cmpl(v)) / geom_nrm_sq(v) == I_3dp);
        CHECK(wdg(B, r_cmpl(B)) / geom_nrm_sq(B) == I_3dp);
        CHECK(wdg(t, r_cmpl(t)) / geom_nrm_sq(t) == I_3dp);
        CHECK(wdg(ps, r_cmpl(ps)) / geom_nrm_sq(ps) == I_3dp);

        // correspondence of complements with geometric products:
        // bulk_duals differ from complements in pga3dp (influence of degenerate metric)
        // left and right bulk duals correspond to
        // left and right complements of the bulk in pga3dp
        CHECK((I_3dp * rev(s)) == l_bulk_dual(s));
        CHECK((I_3dp * rev(v)) == l_bulk_dual(v));
        CHECK((I_3dp * rev(B)) == l_bulk_dual(B));
        CHECK((I_3dp * rev(t)) == l_bulk_dual(t));
        CHECK((I_3dp * rev(ps)) == l_bulk_dual(ps));

        CHECK((rev(s) * I_3dp) == r_bulk_dual(s));
        CHECK((rev(v) * I_3dp) == r_bulk_dual(v));
        CHECK((rev(B) * I_3dp) == r_bulk_dual(B));
        CHECK((rev(t) * I_3dp) == r_bulk_dual(t));
        CHECK((rev(ps) * I_3dp) == r_bulk_dual(ps));


        // linearity of the right and left complement operation
        double a = 2.0;
        CHECK(r_cmpl(a * s) == a * r_cmpl(s));
        CHECK(r_cmpl(a * v) == a * r_cmpl(v));
        CHECK(r_cmpl(a * B) == a * r_cmpl(B));
        CHECK(r_cmpl(a * t) == a * r_cmpl(t));
        CHECK(r_cmpl(a * ps) == a * r_cmpl(ps));

        CHECK(r_cmpl(s + s2) == r_cmpl(s) + r_cmpl(s2));
        CHECK(r_cmpl(v + v2) == r_cmpl(v) + r_cmpl(v2));
        CHECK(r_cmpl(B + B2) == r_cmpl(B) + r_cmpl(B2));
        CHECK(r_cmpl(t + t2) == r_cmpl(t) + r_cmpl(t2));
        CHECK(r_cmpl(ps + ps2) == r_cmpl(ps) + r_cmpl(ps2));

        CHECK(l_cmpl(a * s) == a * l_cmpl(s));
        CHECK(l_cmpl(a * v) == a * l_cmpl(v));
        CHECK(l_cmpl(a * B) == a * l_cmpl(B));
        CHECK(l_cmpl(a * t) == a * l_cmpl(t));
        CHECK(l_cmpl(a * ps) == a * l_cmpl(ps));

        CHECK(l_cmpl(s + s2) == l_cmpl(s) + l_cmpl(s2));
        CHECK(l_cmpl(v + v2) == l_cmpl(v) + l_cmpl(v2));
        CHECK(l_cmpl(B + B2) == l_cmpl(B) + l_cmpl(B2));
        CHECK(l_cmpl(t + t2) == l_cmpl(t) + l_cmpl(t2));
        CHECK(l_cmpl(ps + ps2) == l_cmpl(ps) + l_cmpl(ps2));

        // check regressive wedge product expressed in terms of
        // wegde product and complement operation
        CHECK(rwdg(t, t2) == l_cmpl(wdg(r_cmpl(t), r_cmpl(t2))));
        CHECK(rwdg(t, B) == l_cmpl(wdg(r_cmpl(t), r_cmpl(B))));
        CHECK(rwdg(B, t) == l_cmpl(wdg(r_cmpl(B), r_cmpl(t))));
        CHECK(rwdg(t, v) == l_cmpl(wdg(r_cmpl(t), r_cmpl(v))));
        CHECK(rwdg(v, t) == l_cmpl(wdg(r_cmpl(v), r_cmpl(t))));
        CHECK(rwdg(B, B2) == l_cmpl(wdg(r_cmpl(B), r_cmpl(B2))));

        CHECK(l_cmpl(rwdg(B, B2)) == wdg(r_cmpl(B), r_cmpl(B2)));
        CHECK(l_cmpl(wdg(v, v2)) == rwdg(r_cmpl(v), r_cmpl(v2)));

        // origin and horizon
        CHECK(H_3dp == r_cmpl(O_3dp));
        CHECK(O_3dp == l_cmpl(H_3dp));
    }

    TEST_CASE("MVec3dp: bulk and weight")
    {
        fmt::println("MVec3dp: bulk and weight");

        auto p1 = vec3dp{-2.0, 1.0, 1.0, 1.0};
        auto p2 = vec3dp{2.0, 1.0, 2.0, 2.0};
        auto p3 = vec3dp{2.0, 1.0, 1.0, 3.0};

        auto l1 = bivec3dp{2.0, 1.0, -1.0, 20.0, 10.0, -10.0};
        auto l2 = bivec3dp{-2.0, -1.0, 3.0, -20.0, -10.0, 30.0};
        auto l3 = wdg(p1, p3);

        auto t1 = trivec3dp{-6.0, 3.0, 3.0, 3.0};
        auto t2 = trivec3dp{6.0, 3.0, 6.0, 6.0};
        auto t3 = trivec3dp{6.0, 3.0, 3.0, 9.0};

        // bulk and weight
        CHECK(bulk(p1) == vec3dp{p1.x, p1.y, p1.z, 0.0});
        CHECK(weight(p1) == vec3dp{0.0, 0.0, 0.0, p1.w});
        CHECK(bulk(l1) == bivec3dp{0.0, 0.0, 0.0, l1.mx, l1.my, l1.mz});
        CHECK(weight(l1) == bivec3dp{l1.vx, l1.vy, l1.vz, 0.0, 0.0, 0.0});
        CHECK(bulk(t1) == trivec3dp{0.0, 0.0, 0.0, t1.w});
        CHECK(weight(t1) == trivec3dp{t1.x, t1.y, t1.z, 0.0});

        // bulk_nrm_sq and weight_nrm_sq
        CHECK(bulk_nrm_sq(p1) == p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);
        CHECK(weight_nrm_sq(p1) == p1.w * p1.w);
        CHECK(bulk_nrm_sq(l1) == l1.mx * l1.mx + l1.my * l1.my + l1.mz * l1.mz);
        CHECK(weight_nrm_sq(l1) == l1.vx * l1.vx + l1.vy * l1.vy + l1.vz * l1.vz);
        CHECK(bulk_nrm_sq(t1) == t1.w * t1.w);
        CHECK(weight_nrm_sq(t1) == t1.x * t1.x + t1.y * t1.y + t1.z * t1.z);

        // bulk_nrm and weight_nrm
        CHECK(bulk_nrm(p1) == std::sqrt(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z));
        CHECK(weight_nrm(p1) == std::sqrt(p1.w * p1.w));
        CHECK(bulk_nrm(l1) == std::sqrt(l1.mx * l1.mx + l1.my * l1.my + l1.mz * l1.mz));
        CHECK(weight_nrm(l1) == std::sqrt(l1.vx * l1.vx + l1.vy * l1.vy + l1.vz * l1.vz));
        CHECK(bulk_nrm(t1) == std::sqrt(t1.w * t1.w));
        CHECK(weight_nrm(t1) == std::sqrt(t1.x * t1.x + t1.y * t1.y + t1.z * t1.z));


        // geom_nrm
        CHECK(geom_nrm(p1).c0 == bulk_nrm(p1));
        CHECK(geom_nrm(p1).c1 == weight_nrm(p1));
        CHECK(geom_nrm(l1).c0 == bulk_nrm(l1));
        CHECK(geom_nrm(l1).c1 == weight_nrm(l1));
        CHECK(geom_nrm(t1).c0 == bulk_nrm(t1));
        CHECK(geom_nrm(t1).c1 == weight_nrm(t1));
        CHECK(geom_nrm(l1).c0 / geom_nrm(l1).c1 == geom_nrm(l2).c0 / geom_nrm(l2).c1);

        // attitude
        CHECK(att(p1) == p1.w);
        CHECK(att(p1) == rwdg(p1, r_cmpl(e4_3dp)));
        CHECK(att(p2) == p2.w);
        CHECK(att(p2) == rwdg(p2, r_cmpl(e4_3dp)));
        CHECK(att(l1) == vec3dp{l1.vx, l1.vy, l1.vz, 0.0});
        CHECK(att(l1) == rwdg(l1, r_cmpl(e4_3dp)));
        CHECK(att(l2) == rwdg(l2, r_cmpl(e4_3dp)));
        CHECK(att(l3) == rwdg(l3, r_cmpl(e4_3dp)));
        CHECK(att(t1) == bivec3dp{0.0, 0.0, 0.0, t1.x, t1.y, t1.z});
        CHECK(att(t1) == rwdg(t1, r_cmpl(e4_3dp)));
        CHECK(att(t2) == bivec3dp{0.0, 0.0, 0.0, t2.x, t2.y, t2.z});
        CHECK(att(t2) == rwdg(t2, r_cmpl(e4_3dp)));
        CHECK(att(t3) == bivec3dp{0.0, 0.0, 0.0, t3.x, t3.y, t3.z});
        CHECK(att(t3) == rwdg(t3, r_cmpl(e4_3dp)));
    }

    TEST_CASE("MVec3dp: euclidean distance")
    {
        fmt::println("MVec3dp: euclidean distance");

        auto p0 = O_3dp;
        auto p1 = vec3dp{0.0, 3.0, 0.0, 1.0};
        auto p2 = vec3dp{1.0, 3.0, 0.0, 1.0};
        auto p3 = vec3dp{0.0, 0.0, 1.0, 1.0};

        auto p1u = unitize(p1);
        auto p2u = unitize(p2);
        auto delta = p2u - p1u;

        // fmt::println("");
        // fmt::println("p1  = {}, p2  = {}", p1, p2);
        // fmt::println("p1u = {}, p2u = {}", p1u, p2u);
        // fmt::println("delta = {}, nrm(delta) = {}", delta, nrm(delta));
        // fmt::println("");

        auto l1 = join(p1, p2);
        auto l2 = join(p0, p3);

        auto t1 = e431_3dp;


        auto dist_p0_p1 = dist3dp(p0, p1);
        auto dist_p1_p2 = dist3dp(p1, p2);
        auto dist_p1_l1 = dist3dp(p1, l1);
        auto dist_l1_p1 = dist3dp(l1, p1);
        auto dist_p2_t1 = dist3dp(p2, t1);
        auto dist_l1_l2 = dist3dp(l1, l2);

        // fmt::println("");
        // fmt::println("dist_p1_p2 = {}, dist = {}", dist_p1_p2,
        //              dist_p1_p2.c0 / dist_p1_p2.c1);
        // fmt::println("dist_p1_l1 = {}, dist = {}", dist_p1_l1,
        //              dist_p1_l1.c0 / dist_p1_l1.c1);
        // fmt::println("dist_l1_p1 = {}, dist = {}", dist_l1_p1,
        //              dist_l1_p1.c0 / dist_l1_p1.c1);
        // fmt::println("");
        // fmt::println("dist_p1_t1 = {}, dist = {}", dist_p1_t1,
        //              dist_p1_t1.c0 / dist_p1_t1.c1);
        // fmt::println("");
        // fmt::println("dist_l1_l2 = {}, dist = {}", dist_l1_l2,
        //              dist_l1_l2.c0 / dist_l1_l2.c1);
        // fmt::println("");

        CHECK(to_val(bulk_nrm(delta)) == dist_p1_p2.c0 / dist_p1_p2.c1);
        CHECK(bulk_nrm(bulk_normalize(delta)) == dist_p1_p2.c0 / dist_p1_p2.c1);
        CHECK(dist_p1_p2.c0 / dist_p1_p2.c1 == 1);
        CHECK(dist_p1_l1 == dist_l1_p1);
        CHECK(dist_p2_t1.c0 / dist_p2_t1.c1 == 3);
        CHECK(dist_l1_l2 == dist_p0_p1);

        // auto res = rwdg(pscalar3dp(2.5), H_3dp);
        // fmt::println("res = {}", res);
    }

    TEST_CASE("MVec3dp: product tests (dot, rdot)")
    {
        fmt::println("MVec3dp: product tests (dot, rdot)");
        // tests based on "The inner products of geometric algebra", Leo Dorst


        // auto s1 = scalar3dp{2.0};
        auto v1 = vec3dp{1.0, -3.0, 0.0, 1.0};
        auto b1 = bivec3dp{2.0, -4.0, 1.0, -4.0, 8.0, -2.0};
        // auto t1 = trivec3dp{1.0, -5.0, 6.0, 7.0};
        // auto ps1 = pscalar3dp{-2.0};

        // auto s2 = scalar3dp{-1.0};
        // auto v2 = vec3dp{2.0, 1.5, 0.0};
        // // auto v2 = vec3dp{0.0, 1.0, 0.0};
        // auto b2 = bivec3dp{-2.0, 10.5, 3.0};
        // auto ps2 = pscalar3dp{3.0};

        // auto s3 = scalar3dp{-10.0};
        // auto v3 = vec3dp{-2.0, 4.0, -10.0};
        // auto b3 = bivec3dp{1.0, -2.0, 3.0};
        // auto ps3 = pscalar3dp{-12.0};

        // dot product (difference due to different definition of scalar product)
        CHECK(dot(scalar3dp(1.0), scalar3dp(1.0)) == scalar3dp(1.0));
        CHECK(dot(e1_3dp, e1_3dp) == scalar3dp(1.0));
        CHECK(dot(e2_3dp, e2_3dp) == scalar3dp(1.0));
        CHECK(dot(e3_3dp, e3_3dp) == scalar3dp(1.0));
        CHECK(dot(e4_3dp, e4_3dp) == scalar3dp(0.0));

        CHECK(dot(e41_3dp, e41_3dp) == scalar3dp(0.0));
        CHECK(dot(e42_3dp, e42_3dp) == scalar3dp(0.0));
        CHECK(dot(e43_3dp, e43_3dp) == scalar3dp(0.0));
        CHECK(dot(e23_3dp, e23_3dp) == scalar3dp(1.0));
        CHECK(dot(e31_3dp, e31_3dp) == scalar3dp(1.0));
        CHECK(dot(e12_3dp, e12_3dp) == scalar3dp(1.0));

        CHECK(dot(e423_3dp, e423_3dp) == scalar3dp(0.0));
        CHECK(dot(e431_3dp, e431_3dp) == scalar3dp(0.0));
        CHECK(dot(e412_3dp, e412_3dp) == scalar3dp(0.0));
        CHECK(dot(e321_3dp, e321_3dp) == scalar3dp(1.0));
        CHECK(dot(pscalar3dp(1.0), pscalar3dp(1.0)) == scalar3dp(0.0));

        // regressive dot product
        CHECK(rdot(scalar3dp(1.0), scalar3dp(1.0)) == pscalar3dp(0.0));
        CHECK(rdot(e1_3dp, e1_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e2_3dp, e2_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e3_3dp, e3_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e4_3dp, e4_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e41_3dp, e41_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e42_3dp, e42_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e43_3dp, e43_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e23_3dp, e23_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e31_3dp, e31_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e12_3dp, e12_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e423_3dp, e423_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e431_3dp, e431_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e412_3dp, e412_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e321_3dp, e321_3dp) == pscalar3dp(0.0));
        CHECK(rdot(pscalar3dp(1.0), pscalar3dp(1.0)) == pscalar3dp(1.0));

        // commutator and contraction are equal for a specific case with bivectors
        CHECK(cmt(v1, b1) == (b1 >> v1));
        CHECK(cmt(b1, v1) == (v1 << b1));
        CHECK(cmt(v1, b1) == -cmt(b1, v1));
    }

    TEST_CASE("MVec3dp: simple applications, complements, contraction, "
              "expansions")
    {
        fmt::println("MVec3dp: simple applications, complements, contraction, "
                     "expansions");

        auto s1 = scalar3dp{2.0};
        auto v1 = vec3dp{1.0, -3.0, 0.0, 1.0};
        auto b1 = bivec3dp{2.0, -4.0, 1.0, -4.0, 8.0, -2.0};
        auto t1 = trivec3dp{1.0, -5.0, 6.0, 7.0};
        auto ps1 = pscalar3dp{-2.0};
        auto M1 = mvec3dp{s1, v1, b1, t1, ps1};

        auto s2 = scalar3dp{1.0};
        auto v2 = vec3dp{-1.0, 3.0, 2.0, 1.0};
        auto b2 = bivec3dp{-2.0, 6.0, -1.0, 12.0, 6.0, -3.0};
        auto t2 = trivec3dp{-1.0, -3.0, -1.5, 2.0};
        auto ps2 = pscalar3dp{1.0};
        auto M2 = mvec3dp{s2, v2, b2, t2, ps2};

        // connection between inner product and geometric product
        CHECK(dot(M1, M2) == gr0(M1 * rev(M2)));
        CHECK(dot(M1, M2) == gr0(M2 * rev(M1)));

        // fmt::println("");
        // fmt::println("   M1         = {:g}", M1);
        // fmt::println("   M2         = {:g}", M2);
        // fmt::println("   dot(M1,M2) = {:g}", dot(M1, M2));
        // fmt::println("   M1*rev(M2) = {:g}", M1 * rev(M2));
        // fmt::println("   M2*rev(M1) = {:g}", M2 * rev(M1));
        // fmt::println("");

        // connection between contraction and regressive wedge with dualized args
        l_bulk_contract3dp(M1, M2) = rwdg(l_cmpl(M1), M2);
        r_bulk_contract3dp(M1, M2) = rwdg(M1, r_cmpl(M2));
    }

    TEST_CASE("MVec3dp: join and meet (wdg, rwdg)")
    {
        fmt::println("MVec3dp: join and meet (wdg, rwdg)");

        auto p1 = point3d{2, 0, 1};
        auto p2 = point3d{4, 3, 1};
        auto l12 = join(p1, p2); // = wdg(p1, p2)

        auto p3 = point3d{0, 6, 1};
        auto pl123 = unitize(join(l12, p3)); // plane parallel to e12

        auto ln1 = unitize(join(point3d{0, 3, 0}, point3d{0, 3, 3}));
        auto pt12 = unitize(meet(ln1, pl123)); // should provide pt(0,3,1)

        auto pl123_2 = unitize(Plane3d(p1, p2, p3)); // check the ctor

        auto pl456 =
            unitize(Plane3d(point3d{2, 0, 3}, point3d{4, 3, 3}, point3d{0, 6, 3}));

        auto l_inf = meet(pl123, pl456); // gives a line at infinity
                                         // weighted by the distance of the planes

        // fmt::println("");
        // fmt::println("p1 = {}, p2 = {}, p3 = {}", p1, p2, p3);
        // fmt::println("l12      = {}, pl123      = {}", l12, pl123);
        // fmt::println("att(l12) = {}, att(pl123) = {}", att(l12), att(pl123));
        // fmt::println("pt12  = meet(l1, pl123)   = {}", pt12);
        // fmt::println("");
        // fmt::println("pl456 = {}", pl456);
        // fmt::println("l_inf = meet(pl123, pl456) = {}", l_inf);
        // fmt::println("att(l_inf)) = {}", att(l_inf));

        CHECK(join(p1, p2) == wdg(p1, p2));
        CHECK(pt12 == vec3dp{0, 3, 1, 1});
        CHECK(pl123 == pl123_2);
        CHECK(att(l_inf) == vec3dp{}); // line is at infinity
    }

    TEST_CASE("MVec3dp: orthogonal projections")
    {
        fmt::println("MVec3dp: orthogonal projections");

        auto pt0 = O_3dp;
        auto pt1 = vec3dp{1, 0, 0, 1};
        auto pt2 = vec3dp{1, 1, 0, 1};
        // auto pt3 = vec3dp{0, 1, 0, 1};
        // auto pt4 = vec3dp{0, 1, 1, 1};
        auto pt5 = vec3dp{0, 0, 1, 1};
        auto pt6 = vec3dp{1, 0, 1, 1};
        auto pt7 = vec3dp{1, 1, 1, 1};
        // auto pt8 = vec3dp{1, -1, 0, 1};
        auto pt9 = vec3dp{0, -1, 0, 1};

        auto ln01 = join(pt0, pt1); // line containing both points
        auto ln02 = join(pt0, pt2);
        auto ln06 = join(pt0, pt6);
        auto ln07 = join(pt0, pt7);
        auto ln09 = join(pt0, pt9);
        auto ln56 = join(pt5, pt6);

        auto a12 = angle(ln01, ln02);
        auto a21 = angle(ln02, ln01);
        auto a29 = angle(ln02, ln09);
        auto a92 = angle(ln09, ln02);
        auto a67 = angle(ln06, ln07);
        auto a256 = angle(ln02, ln56);

        auto pl = e431_3dp; // plane 31 "at the bottom" of the coordinate system
        auto pl2 = join(join(pt0, pt7), pt5);

        auto apl07 = angle(pl, ln07);
        auto apl70 = angle(ln07, pl);

        auto apl02 = angle(pl, pl2);

        // fmt::println("");
        // fmt::println("pt0  = {}, pt1  = {}, pt2  = {}, pt3  = {}", pt0, pt1,
        // pt2, pt3); fmt::println("pt4  = {}, pt5  = {}, pt6  = {}, pt7  = {}",
        // pt4, pt5, pt6, pt7); fmt::println("ln02 = {}", ln01);
        // fmt::println("ln02 = {}", ln02); fmt::println("ln07 = {}", ln07);
        // fmt::println("ln09 = {}", ln09); fmt::println("ln56 = {}", ln56);
        // fmt::println("pl = {}", pl);
        // fmt::println("pl2 = {}", pl2);
        // fmt::println("");
        // fmt::println("a12 = {}", rad2deg(a12));
        // fmt::println("a21 = {}", rad2deg(a21));
        // fmt::println("a29 = {}", rad2deg(a29));
        // fmt::println("a92 = {}", rad2deg(a92));
        // fmt::println("a67 = {}", rad2deg(a67));
        // fmt::println("a256 = {}", rad2deg(a256));
        // fmt::println("");
        // fmt::println("apl07 = {}", rad2deg(apl07));
        // fmt::println("apl70 = {}", rad2deg(apl70));
        // fmt::println("");
        // fmt::println("apl02 = {}", rad2deg(apl02));

        auto pt1p = ortho_proj3dp(pt7, pl);
        CHECK(pt1p == vec3dp{1, 0, 1, 1});

        auto pt4p = unitize(ortho_proj3dp(pt1, ln02));
        auto pt7p = unitize(ortho_proj3dp(vec3dp{1, 0, 1, 1}, ln07));
        auto s = 0.5;
        auto l = 2. / 3.;

        // fmt::println("pt4p  = {}", pt4p);
        // fmt::println("pt7p  = {}, dist07 = {}", pt7p, dist3dp(pt0, pt7));
        // fmt::println("dist07p = {}", dist3dp(pt0, pt7p));

        CHECK(pt4p == vec3dp{s, s, 0, 1});
        CHECK(pt7p == vec3dp{l, l, l, 1});

        auto ln67p = join(pt6, pt7p);
        auto phi = angle(ln67p, ln07);
        auto phi0607 = angle(ln06, ln07);

        // fmt::println("phi = {}", rad2deg(phi));
        // fmt::println("phi0607 = {}", rad2deg(phi0607));

        CHECK(phi0607 - std::atan(1. / std::sqrt(2.0)) < eps);

        auto ln7p = unitize(ortho_proj3dp(ln07, pl));
        CHECK(ln7p == unitize(ln06));
        CHECK(abs(phi - pi / 2) < eps);
        CHECK(apl07 == apl70);
        CHECK(a21 == a12);
        CHECK(a29 == a92);
        CHECK(a67 - apl07 < eps);
        CHECK(rad2deg(apl02) == 135);
        CHECK(a256 == a12);
    }

    TEST_CASE("MVec3dp: central projections")
    {
        fmt::println("MVec3dp: central projections");

        // central projection onto plane towards origin
        auto p1 = vec3dp{1, 0, 1, 1};
        auto p2 = vec3dp{1, 0, 0, 1};
        auto p3 = vec3dp{0, 1, 0, 1};
        auto p = vec3dp{10, 10, 0, 1};

        auto pl = join(join(p1, p2), p3);
        auto pp = unitize(central_proj3dp(p, pl));

        // fmt::println("");
        // fmt::println("pl = {}", pl);
        // fmt::println("p  = {}, pp  = {}", p, pp);
        CHECK(pp == vec3dp(0.5, 0.5, 0, 1));

        // central projection onto line towards origin
        auto ln = join(p2, p3);
        auto lpt = unitize(central_proj3dp(p, ln));

        // fmt::println("");
        // fmt::println("ln = {}", ln);
        // fmt::println("p  = {}, lpt  = {}", p, lpt);
        CHECK(lpt == vec3dp(0.5, 0.5, 0, 1));

        // central projection of line onto plane towards origin
        auto p4 = vec3dp{2, 0, 2, 1};
        auto p5 = vec3dp{0, 2, 0, 1};
        auto ln13 = unitize(join(p1, p3));
        auto ln2 = unitize(join(p4, p5));
        auto ln2p = unitize(central_proj3dp(ln2, pl));

        // fmt::println("");
        // fmt::println("pl = {}", pl);
        // fmt::println("ln2 = {}", ln2);
        // fmt::println("ln13  = {}", ln13);
        // fmt::println("ln2p  = {}", ln2p);
        CHECK(ln2p == ln13);
    }

    TEST_CASE("MVec3dp: orthogonal antiprojections")
    {
        fmt::println("MVec3dp: orthogonal antiprojections");

        // orthogonal anitprojection of a plane to a point,
        // i.e. create a new plane with the same orientation of the normal that
        // contains the point
        auto pl = e431_3dp + e423_3dp;
        auto p1 = vec3dp{7, 3, -2, 1};
        auto dist_p1_pl = dist3dp(p1, pl).c0 / dist3dp(p1, pl).c1;

        auto plp = ortho_antiproj3dp(pl, p1);
        auto dist_p1_plp = dist3dp(p1, plp).c0 / dist3dp(p1, plp).c1;

        // fmt::println("");
        // fmt::println("p1  = {}", p1);
        // fmt::println("pl  = {}, dist_p1_pl  = {}", pl, dist_p1_pl);
        // fmt::println("plp = {}, dist_p1_plp = {}", plp, dist_p1_plp);
        CHECK(weight_nrm(pl) == weight_nrm(plp)); // planes have same orientation
        CHECK(dist_p1_pl != 0.0);                 // point is not in original plane
        CHECK(dist_p1_plp == 0.0);                // point is in antiprojected plane

        // orthogonal anitprojection of a line onto a point,
        // i.e. create a new line with the same orientation that contains the
        // point
        auto p2 = vec3dp{-9, 1, -3, 1};
        auto p3 = vec3dp{1, 2, -3, 1};
        auto ln = join(p1, p2);
        auto dist_p3_ln = dist3dp(p3, ln).c0 / dist3dp(p3, ln).c1;

        auto lnp = ortho_antiproj3dp(ln, p3);
        auto dist_p3_lnp = dist3dp(p3, lnp).c0 / dist3dp(p3, lnp).c1;

        // fmt::println("");
        // fmt::println("p1  = {}", p1);
        // fmt::println("ln  = {}, dist_p1_ln  = {}", ln, dist_p1_ln);
        // fmt::println("lnp = {}, dist_p1_lnp = {}", lnp, dist_p1_lnp);
        CHECK(weight_nrm(pl) == weight_nrm(plp)); // lines have same orientation
        CHECK(dist_p3_ln != 0.0);
        CHECK(dist_p3_lnp == 0.0); // point is in antiprojected line

        auto plp2 = ortho_antiproj3dp(pl, ln);
        auto dist_p1_plp2 = dist3dp(p1, plp2).c0 / dist3dp(p1, plp2).c1;
        // fmt::println("");
        // fmt::println("pl  = {}", pl);
        // fmt::println("ln  = {}", ln);
        // fmt::println("plp2 = {}", plp2);
        CHECK(dist_p1_plp2 == 0);
    }

    TEST_CASE("G<3,0,1>: Scalar3dp and PScalar3dp formatting tests")
    {
        fmt::println("G<3,0,1>: Scalar3dp and PScalar3dp formatting tests");

        // Test scalar formatting with different format specifiers
        scalar3dp s{3.14159};
        pscalar3dp ps{2.71828};

        // Default formatting
        std::string default_scalar = fmt::format("{}", s);
        std::string default_pscalar = fmt::format("{}", ps);
        CHECK(default_scalar == "Scalar3dp(3.14159)");
        CHECK(default_pscalar == "PScalar3dp(2.71828)");

        // Two decimal places
        std::string two_decimals_scalar = fmt::format("{:.2f}", s);
        std::string two_decimals_pscalar = fmt::format("{:.2f}", ps);
        CHECK(two_decimals_scalar == "Scalar3dp(3.14)");
        CHECK(two_decimals_pscalar == "PScalar3dp(2.72)");

        // Scientific notation
        std::string scientific_scalar = fmt::format("{:.2e}", s);
        std::string scientific_pscalar = fmt::format("{:.2e}", ps);
        CHECK(scientific_scalar == "Scalar3dp(3.14e+00)");
        CHECK(scientific_pscalar == "PScalar3dp(2.72e+00)");

        // Six decimal places
        std::string six_decimals_scalar = fmt::format("{:.6f}", s);
        std::string six_decimals_pscalar = fmt::format("{:.6f}", ps);
        CHECK(six_decimals_scalar == "Scalar3dp(3.141590)");
        CHECK(six_decimals_pscalar == "PScalar3dp(2.718280)");

        fmt::println(
            "   All Scalar3dp and PScalar3dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("Vec3dp: comprehensive format specifier tests with pre-computed values")
    {
        fmt::println(
            "Vec3dp: comprehensive format specifier tests with pre-computed values:\n");

        // Test vector with specific values for predictable output
        vec3dp v{3.14159, 2.71828, 1.41421, 1.73205};

        // Default formatting - using value_t precision
        std::string default_format = fmt::format("{}", v);
        CHECK(default_format == "Vec3dp(3.14159,2.71828,1.41421,1.73205)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", v);
        CHECK(two_decimals == "Vec3dp(3.14,2.72,1.41,1.73)");

        // Scientific notation with 2 decimal places
        std::string scientific = fmt::format("{:.2e}", v);
        CHECK(scientific == "Vec3dp(3.14e+00,2.72e+00,1.41e+00,1.73e+00)");

        // Six decimal places
        std::string six_decimals = fmt::format("{:.6f}", v);
        CHECK(six_decimals == "Vec3dp(3.141590,2.718280,1.414210,1.732050)");

        // Test with different vector values for edge cases
        vec3dp v_zero{0.0, 0.0, 0.0, 0.0};
        vec3dp v_negative{-1.5, -2.5, -3.5, -4.5};

        std::string zero_format = fmt::format("{:.1f}", v_zero);
        std::string negative_format = fmt::format("{:.1f}", v_negative);
        CHECK(zero_format == "Vec3dp(0.0,0.0,0.0,0.0)");
        CHECK(negative_format == "Vec3dp(-1.5,-2.5,-3.5,-4.5)");

        fmt::println("   All Vec3dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: TriVec3dp formatting tests")
    {
        fmt::println("G<3,0,1>: TriVec3dp formatting tests");

        // Test trivector formatting with different format specifiers
        trivec3dp tv{1.234, 5.678, 9.012, 2.468};

        // Default formatting
        std::string default_format = fmt::format("{}", tv);
        CHECK(default_format == "TriVec3dp(1.234,5.678,9.012,2.468)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", tv);
        CHECK(two_decimals == "TriVec3dp(1.23,5.68,9.01,2.47)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", tv);
        CHECK(scientific == "TriVec3dp(1.23e+00,5.68e+00,9.01e+00,2.47e+00)");

        // Four decimal places
        std::string four_decimals = fmt::format("{:.4f}", tv);
        CHECK(four_decimals == "TriVec3dp(1.2340,5.6780,9.0120,2.4680)");

        fmt::println("   All TriVec3dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: BiVec3dp formatting tests")
    {
        fmt::println("G<3,0,1>: BiVec3dp formatting tests");

        // Test bivector formatting with different format specifiers
        bivec3dp bv{1.11, 2.22, 3.33, 4.44, 5.55, 6.66};

        // Default formatting
        std::string default_format = fmt::format("{}", bv);
        CHECK(default_format == "BiVec3dp(1.11,2.22,3.33,4.44,5.55,6.66)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", bv);
        CHECK(two_decimals == "BiVec3dp(1.11,2.22,3.33,4.44,5.55,6.66)");

        // Scientific notation
        std::string scientific = fmt::format("{:.1e}", bv);
        CHECK(scientific == "BiVec3dp(1.1e+00,2.2e+00,3.3e+00,4.4e+00,5.5e+00,6.7e+00)");

        // Three decimal places
        std::string three_decimals = fmt::format("{:.3f}", bv);
        CHECK(three_decimals == "BiVec3dp(1.110,2.220,3.330,4.440,5.550,6.660)");

        fmt::println("   All BiVec3dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - default init")
    {
        fmt::println("G<3,0,1>: DualNum3dp - default init");
        dualnum3dp dn;
        CHECK(abs(dn.c0) < eps);
        CHECK(abs(dn.c1) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - constructor from two components")
    {
        fmt::println("G<3,0,1>: DualNum3dp - constructor from two components");
        dualnum3dp dn{2.5, 3.7};
        CHECK(abs(dn.c0 - 2.5) < eps);
        CHECK(abs(dn.c1 - 3.7) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - constructor from scalar")
    {
        fmt::println("G<3,0,1>: DualNum3dp - constructor from scalar");
        scalar3dp s{5.0};
        dualnum3dp dn{s};
        CHECK(abs(dn.c0 - 5.0) < eps);
        CHECK(abs(dn.c1) < eps);
        CHECK(gr0(dn) == s);
        CHECK(gr4(dn) == pscalar3dp{0.0});
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - constructor from pseudoscalar")
    {
        fmt::println("G<3,0,1>: DualNum3dp - constructor from pseudoscalar");
        pscalar3dp ps{7.5};
        dualnum3dp dn{ps};
        CHECK(abs(dn.c0) < eps);
        CHECK(abs(dn.c1 - 7.5) < eps);
        CHECK(gr0(dn) == scalar3dp{0.0});
        CHECK(gr4(dn) == ps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - constructor from scalar and pseudoscalar")
    {
        fmt::println("G<3,0,1>: DualNum3dp - constructor from scalar and pseudoscalar");
        scalar3dp s{3.0};
        pscalar3dp ps{4.0};
        dualnum3dp dn{s, ps};
        CHECK(abs(dn.c0 - 3.0) < eps);
        CHECK(abs(dn.c1 - 4.0) < eps);
        CHECK(gr0(dn) == s);
        CHECK(gr4(dn) == ps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - copy constructor and assignment")
    {
        fmt::println("G<3,0,1>: DualNum3dp - copy constructor and assignment");
        dualnum3dp dn1{2.5, 3.7};
        dualnum3dp dn2{dn1};   // copy constructor
        dualnum3dp dn3 = dn1;  // copy assignment
        dualnum3dp dn4 = -dn1; // unary minus

        CHECK(abs(dn2.c0 - 2.5) < eps);
        CHECK(abs(dn2.c1 - 3.7) < eps);
        CHECK(abs(dn3.c0 - 2.5) < eps);
        CHECK(abs(dn3.c1 - 3.7) < eps);
        CHECK(abs(dn4.c0 + 2.5) < eps);
        CHECK(abs(dn4.c1 + 3.7) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - grade extraction")
    {
        fmt::println("G<3,0,1>: DualNum3dp - grade extraction");
        dualnum3dp dn{5.0, 7.0};

        scalar3dp s = gr0(dn);
        pscalar3dp ps = gr4(dn);

        CHECK(abs(value_t(s) - 5.0) < eps);
        CHECK(abs(value_t(ps) - 7.0) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - addition operations")
    {
        fmt::println("G<3,0,1>: DualNum3dp - addition operations");

        scalar3dp s{2.0};
        pscalar3dp ps{3.0};
        dualnum3dp dn1{1.0, 2.0};
        dualnum3dp dn2{3.0, 4.0};

        // scalar + pseudoscalar => dual number (explicit construction to avoid ambiguity)
        // MSVC workaround: explicit cast to value_t for cross-platform compatibility
        dualnum3dp result1{s, ps};
        CHECK(abs(result1.c0 - 2.0) < eps);
        CHECK(abs(result1.c1 - 3.0) < eps);

        // scalar + dual number
        dualnum3dp result3 = s + dn1;
        CHECK(abs(result3.c0 - 3.0) < eps);
        CHECK(abs(result3.c1 - 2.0) < eps);

        // dual number + scalar
        dualnum3dp result4 = dn1 + s;
        CHECK(abs(result4.c0 - 3.0) < eps);
        CHECK(abs(result4.c1 - 2.0) < eps);

        // pseudoscalar + dual number
        dualnum3dp result5 = ps + dn1;
        CHECK(abs(result5.c0 - 1.0) < eps);
        CHECK(abs(result5.c1 - 5.0) < eps);

        // dual number + pseudoscalar
        dualnum3dp result6 = dn1 + ps;
        CHECK(abs(result6.c0 - 1.0) < eps);
        CHECK(abs(result6.c1 - 5.0) < eps);

        // dual number + dual number (using compound assignment)
        dualnum3dp dn3 = dn1;
        dn3 += dn2;
        CHECK(abs(dn3.c0 - 4.0) < eps);
        CHECK(abs(dn3.c1 - 6.0) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - subtraction operations")
    {
        fmt::println("G<3,0,1>: DualNum3dp - subtraction operations");

        scalar3dp s{5.0};
        pscalar3dp ps{3.0};
        dualnum3dp dn1{7.0, 4.0};
        dualnum3dp dn2{2.0, 1.0};

        // scalar - dual number
        dualnum3dp result3 = s - dn1;
        CHECK(abs(result3.c0 + 2.0) < eps);
        CHECK(abs(result3.c1 + 4.0) < eps);

        // dual number - scalar
        dualnum3dp result4 = dn1 - s;
        CHECK(abs(result4.c0 - 2.0) < eps);
        CHECK(abs(result4.c1 - 4.0) < eps);

        // pseudoscalar - dual number
        dualnum3dp result5 = ps - dn1;
        CHECK(abs(result5.c0 + 7.0) < eps);
        CHECK(abs(result5.c1 + 1.0) < eps);

        // dual number - pseudoscalar
        dualnum3dp result6 = dn1 - ps;
        CHECK(abs(result6.c0 - 7.0) < eps);
        CHECK(abs(result6.c1 - 1.0) < eps);

        // dual number - dual number (using compound assignment)
        dualnum3dp dn3 = dn1;
        dn3 -= dn2;
        CHECK(abs(dn3.c0 - 5.0) < eps);
        CHECK(abs(dn3.c1 - 3.0) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - scalar multiplication and division")
    {
        fmt::println("G<3,0,1>: DualNum3dp - scalar multiplication and division");

        dualnum3dp dn{6.0, 8.0};

        // scalar multiplication
        dualnum3dp dn2 = dn;
        dn2 *= 2.0;
        CHECK(abs(dn2.c0 - 12.0) < eps);
        CHECK(abs(dn2.c1 - 16.0) < eps);

        // scalar division
        dualnum3dp dn3 = dn;
        dn3 /= 2.0;
        CHECK(abs(dn3.c0 - 3.0) < eps);
        CHECK(abs(dn3.c1 - 4.0) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - compound assignment operators")
    {
        fmt::println("G<3,0,1>: DualNum3dp - compound assignment operators");

        dualnum3dp dn1{10.0, 20.0};
        dualnum3dp dn2{3.0, 5.0};

        // test +=
        dualnum3dp result1 = dn1;
        result1 += dn2;
        CHECK(abs(result1.c0 - 13.0) < eps);
        CHECK(abs(result1.c1 - 25.0) < eps);

        // test -=
        dualnum3dp result2 = dn1;
        result2 -= dn2;
        CHECK(abs(result2.c0 - 7.0) < eps);
        CHECK(abs(result2.c1 - 15.0) < eps);

        // test *=
        dualnum3dp result3 = dn1;
        result3 *= 0.5;
        CHECK(abs(result3.c0 - 5.0) < eps);
        CHECK(abs(result3.c1 - 10.0) < eps);

        // test /=
        dualnum3dp result4 = dn1;
        result4 /= 2.0;
        CHECK(abs(result4.c0 - 5.0) < eps);
        CHECK(abs(result4.c1 - 10.0) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - mixed type operations")
    {
        fmt::println("G<3,0,1>: DualNum3dp - mixed type operations");

        // Test that operations work with different floating point types
        scalar3dp s{2.0};
        pscalar3dp ps{3.0};
        dualnum3dp dn{1.5, 2.5};

        // Test dual number operations with scalar and pseudoscalar
        dualnum3dp result2 = dn + s;
        CHECK(abs(result2.c0 - 3.5) < eps);
        CHECK(abs(result2.c1 - 2.5) < eps);

        dualnum3dp result3 = ps + dn;
        CHECK(abs(result3.c0 - 1.5) < eps);
        CHECK(abs(result3.c1 - 5.5) < eps);

        dualnum3dp result4 = dn - s;
        CHECK(abs(result4.c0 + 0.5) < eps);
        CHECK(abs(result4.c1 - 2.5) < eps);
    }

    TEST_CASE("G<3,0,1>: DualNum3dp - property: I_3dp^2 = 0")
    {
        fmt::println("G<3,0,1>: DualNum3dp - property: I_3dp^2 = 0");

        // The pseudoscalar in PGA3DP squares to zero (degenerate metric)
        // For dual numbers: (c0 + c1*I)(c0 + c1*I) = c0^2 + 2*c0*c1*I + c1^2*I^2
        // Since I^2 = 0: result = c0^2 + 2*c0*c1*I

        // Verify that I_3dp^2 = 0 through geometric product (using operator*)
        auto I_squared = I_3dp * I_3dp;
        CHECK(abs(value_t(I_squared)) < eps);

        fmt::println("I_3dp^2 = 0 (degenerate metric property verified)");
    }

    TEST_CASE("G<3,0,1>: DualNum3dp formatting tests")
    {
        fmt::println("G<3,0,1>: DualNum3dp formatting tests");

        // Test dual number formatting with different format specifiers
        dualnum3dp dn{2.5, 3.7};

        // Default formatting
        std::string default_format = fmt::format("{}", dn);
        CHECK(default_format == "DualNum3dp(2.5,3.7)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", dn);
        CHECK(two_decimals == "DualNum3dp(2.50,3.70)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", dn);
        CHECK(scientific == "DualNum3dp(2.50e+00,3.70e+00)");

        // Three decimal places
        std::string three_decimals = fmt::format("{:.3f}", dn);
        CHECK(three_decimals == "DualNum3dp(2.500,3.700)");

        fmt::println("   All DualNum3dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: MVec3dp_E, MVec3dp_U, and MVec3dp formatting tests")
    {
        fmt::println("G<3,0,1>: MVec3dp_E, MVec3dp_U, and MVec3dp formatting tests");

        // Test multivector even formatting
        mvec3dp_e me{1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8};
        std::string even_default = fmt::format("{}", me);
        std::string even_two_dec = fmt::format("{:.1f}", me);
        CHECK(even_default == "MVec3dp_E(1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8)");
        CHECK(even_two_dec == "MVec3dp_E(1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8)");

        // Test multivector odd formatting
        mvec3dp_u mu{1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
        std::string odd_default = fmt::format("{}", mu);
        std::string odd_two_dec = fmt::format("{:.1f}", mu);
        CHECK(odd_default == "MVec3dp_U(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)");
        CHECK(odd_two_dec == "MVec3dp_U(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)");

        // Test full multivector formatting (16 components)
        mvec3dp m{1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,
                  9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
        std::string full_default = fmt::format("{}", m);
        std::string full_scientific = fmt::format("{:.1e}", m);
        CHECK(full_default == "MVec3dp(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)");
        CHECK(full_scientific == "MVec3dp(1.0e+00,2.0e+00,3.0e+00,4.0e+00,5.0e+00,"
                                 "6.0e+00,7.0e+00,8.0e+00,9.0e+00,1.0e+01,1.1e+01,"
                                 "1.2e+01,1.3e+01,1.4e+01,1.5e+01,1.6e+01)");

        fmt::println("   All MVec3dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("PGA3DP: congruence tests")
    {
        fmt::println("PGA3DP: congruence tests");

        // Test scalars
        scalar3dp s1{5.0};
        scalar3dp s2{-3.0}; // different sign
        scalar3dp s3{2.5};  // same sign as s1
        scalar3dp s_zero{0.0};

        CHECK(is_congruent(s1, s2) == true);         // different signs are congruent
        CHECK(is_congruent(s1, s3) == true);         // same signs are congruent
        CHECK(is_congruent(s1, s_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(s_zero, s_zero) == true); // zero vs zero

        // Test vectors (points in PGA3DP)
        vec3dp v1{1.0, 0.0, 0.0, 1.0};   // point on x-axis
        vec3dp v2{0.0, 1.0, 0.0, 1.0};   // point on y-axis
        vec3dp v3{2.0, 0.0, 0.0, 2.0};   // parallel to v1
        vec3dp v4{-1.0, 0.0, 0.0, -1.0}; // antiparallel to v1
        vec3dp v_zero{0.0, 0.0, 0.0, 0.0};

        CHECK(is_congruent(v1, v2) == false);        // different directions
        CHECK(is_congruent(v1, v3) == true);         // parallel directions
        CHECK(is_congruent(v1, v4) == true);         // antiparallel directions
        CHECK(is_congruent(v1, v_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(v_zero, v_zero) == true); // zero vs zero

        // Test bivectors (lines in PGA3DP) - use wedge product since it exists in 4D
        bivec3dp b1{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};  // line in xy-plane
        bivec3dp b2{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};  // different line orientation
        bivec3dp b3{2.0, 0.0, 0.0, 0.0, 0.0, 0.0};  // parallel to b1
        bivec3dp b4{-1.5, 0.0, 0.0, 0.0, 0.0, 0.0}; // antiparallel to b1
        bivec3dp b_zero{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        CHECK(is_congruent(b1, b2) == false);        // different line orientations
        CHECK(is_congruent(b1, b3) == true);         // same line orientation
        CHECK(is_congruent(b1, b4) == true);         // same line, opposite orientation
        CHECK(is_congruent(b1, b_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(b_zero, b_zero) == true); // zero vs zero

        // Test trivectors (planes in PGA3DP) - use dot product method
        trivec3dp t1{1.0, 0.0, 0.0, 0.0};  // plane
        trivec3dp t2{0.0, 1.0, 0.0, 0.0};  // different plane orientation
        trivec3dp t3{2.0, 0.0, 0.0, 0.0};  // parallel to t1
        trivec3dp t4{-1.5, 0.0, 0.0, 0.0}; // antiparallel to t1
        trivec3dp t_zero{0.0, 0.0, 0.0, 0.0};

        CHECK(is_congruent(t1, t2) == false);        // different plane orientations
        CHECK(is_congruent(t1, t3) == true);         // same plane orientation
        CHECK(is_congruent(t1, t4) == true);         // same plane, opposite orientation
        CHECK(is_congruent(t1, t_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(t_zero, t_zero) == true); // zero vs zero

        // Test pseudoscalars (max grade in PGA3DP)
        pscalar3dp p1{1.0};
        pscalar3dp p2{2.0};  // different magnitude
        pscalar3dp p3{-1.0}; // opposite sign
        pscalar3dp p_zero{0.0};

        CHECK(is_congruent(p1, p2) == true);      // all non-zero pseudoscalars congruent
        CHECK(is_congruent(p1, p3) == true);      // different signs still congruent
        CHECK(is_congruent(p1, p_zero) == false); // zero vs non-zero
        CHECK(is_congruent(p_zero, p_zero) == true); // zero vs zero

        // fmt::println("   All PGA3DP congruence tests passed");
    }

    TEST_CASE("G<3,0,1>: transwedge and regressive transwedge products")
    {
        fmt::println("G<3,0,1>: transwedge and regressive transwedge products");

        auto s = scalar3dp{4.0};
        auto v = vec3dp{1.0, 2.0, 3.0, 1.0};
        auto v2 = vec3dp{-1.0, 3.0, -2.0, 1.0};
        auto B = bivec3dp{10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
        auto B2 = bivec3dp{-20.0, -40.0, -60.0, -12.0, -15.0, -10.0};
        auto t = trivec3dp{100.0, 200.0, 300.0, 400.0};
        auto t2 = trivec3dp{-200.0, -400.0, -600.0, -120.0};
        auto ps = pscalar3dp{-3.0};

        // check identities of transwedge product
        CHECK(twdg1(ps, v) == ps * v);
        CHECK(twdg1(ps, v) == (v << ps));
        CHECK(twdg1(ps, v) ==
              rwdg(l_bulk_dual(v), ps)); // contraction expressed as rwdg prod.

        CHECK(twdg1(ps, v) == -twdg1(v, ps)); // vec and ps do not commute

        CHECK(twdg1(v, ps) == v * ps);
        CHECK(twdg1(v, ps) == (ps >> v));
        CHECK(twdg1(v, ps) ==
              rwdg(ps, r_bulk_dual(v))); // contraction expressed as rwdg prod.


        CHECK(twdg1(t, v) == (v << t));
        CHECK(twdg1(t, v) == (t >> v));

        CHECK(twdg1(t, v) == twdg1(v, t));

        CHECK(twdg1(v, t) == (t >> v));
        CHECK(twdg1(v, t) == (v << t));

        CHECK(twdg1(t, B) == cmt(t, B));
        CHECK(twdg1(B, t) == cmt(B, t));
        CHECK(twdg1(B, B2) == cmt(B, B2));
        CHECK(twdg1(B, v) == (v << B));
        CHECK(twdg1(v, B) == (B >> v));

        CHECK(twdg1(v, v2) == dot(v, v2));
        CHECK(twdg1(v, v2) == dot(v2, v));

        // check identities of regressive transwedge product
        CHECK(rtwdg1(t, t2) == rdot(t, t2));
        CHECK(rtwdg1(t, t2) == rdot(t2, t));

        CHECK(rtwdg1(t, B) == rcmt(t, B));
        CHECK(rtwdg1(B, t) == rcmt(B, t));

        CHECK(rtwdg1(t, v) == r_weight_expand3dp(v, t));
        CHECK(rtwdg1(t, v) == wdg(v, r_weight_dual(t)));
        CHECK(rtwdg1(v, t) == l_weight_expand3dp(t, v));
        CHECK(rtwdg1(v, t) == wdg(l_weight_dual(t), v));

        CHECK(rtwdg1(t, s) == r_weight_expand3dp(s, t));
        CHECK(rtwdg1(t, s) == wdg(s, r_weight_dual(t)));
        CHECK(rtwdg1(s, t) == l_weight_expand3dp(t, s));
        CHECK(rtwdg1(s, t) == wdg(l_weight_dual(t), s));

        CHECK(rtwdg1(B, B2) == rcmt(B, B2));

        CHECK(rtwdg1(B, v) == rcmt(B, v));
        CHECK(rtwdg1(B, v) == -rtwdg1(v, B)); // vector ans bivector anticommute
        CHECK(rtwdg1(v, B) == rcmt(v, B));
    }

    TEST_CASE("G<3,0,1>: adding and reducing bound vectors")
    {
        fmt::println("G<3,0,1>: adding and reducing bound vectors");

        // case taken from: "Grassmann Algebra, Vol. 1 - Foundations", John Browne, p. 177

        auto P1 = vec3dp{1, 3, -4, 1};
        auto P2 = vec3dp{2, -1, -2, 1};
        auto P3 = vec3dp{-5, 3, -6, 1};
        auto P4 = vec3dp{4, 2, -9, 1};

        auto x1 = vec3dp{1, 0, -1, 0};
        auto x2 = vec3dp{1, -1, 1, 0};
        auto x3 = vec3dp{2, 3, 0, 0};
        auto x4 = vec3dp{0, 0, 5, 0};

        auto L1 = wdg(P1, x1);
        auto L2 = wdg(P2, x2);
        auto L3 = wdg(P3, x3);
        auto L4 = wdg(P4, x4);

        auto dir_vec = x1 + x2 + x3 + x4;
        auto L = L1 + L2 + L3 + L4;
        CHECK(dir_vec == vec3dp{4, 2, 5, 0});
        CHECK(L == bivec3dp{4, 2, 5, 22, -39, -25});

        // moving a bound vector: L = O_3dp ^ dir_vec + bivec_part
        // (= first three components: direction vector; last three components: bivec part)
        // now expressed through another reference point:
        // R = O_3dp + x
        // => L = R ^ dir_vec + L - R ^ dir_vec

        fmt::println("dir_vec =   {}", dir_vec);
        fmt::println("L       = {}", L);
        fmt::println("");
    }

    TEST_CASE("MVec3dp: contraction vs. projection")
    {
        fmt::println("MVec3dp: contraction vs. projection");

        // get vector to poject into coordinate planes
        auto v = vec3dp{1, 1, 1, 0};
        auto B = e31_3dp; // bivector of e31 plane (see ga_usr_consts.hpp)

        // from v ⟑ B = B >> v + v ^ B
        // get  v = (v ⟑ B) ⟑ B^(-1) = (B >> v) ⟑ B^(-1) + (v ^ B) ⟑ B^(-1)
        //                           =      v_parallel   + v_perpendicular

        fmt::println("v = {}", v);
        fmt::println("B = {}", B);
        fmt::println("");
        fmt::println("B >> v = {}", B >> v);
        fmt::println("v ^ B = {}", wdg(v, B));
        fmt::println("");
        fmt::println("inv(B) = {}", inv(B));
        fmt::println("");
        fmt::println("v_para = (B >> v) * inv(B) = {}", gr1((B >> v) * inv(B)));
        fmt::println("v_perp = (v ^ B) * inv(B) = {}", gr1(wdg(v, B) * inv(B)));
        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Complement Mathematical Properties Tests
    // Note: Dual involution tests are excluded for PGA due to degenerate metric
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<3,0,1>: left-right complement composition")
    {
        fmt::println("G<3,0,1>: left-right complement composition");

        // For even-dimensional algebras: l_cmpl(r_cmpl(u)) = r_cmpl(l_cmpl(u)) = u
        scalar3dp s{3.0};
        vec3dp v{2.0, 5.0, 7.0, 1.0};
        bivec3dp B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec3dp t{1.0, 2.0, 3.0, 4.0};
        pscalar3dp ps{11.0};
        mvec3dp_e M_e{s, B, ps};
        mvec3dp_u M_u{v, t};
        mvec3dp M{s, v, B, t, ps};

        CHECK(l_cmpl(r_cmpl(s)) == s);
        CHECK(l_cmpl(r_cmpl(v)) == v);
        CHECK(l_cmpl(r_cmpl(B)) == B);
        CHECK(l_cmpl(r_cmpl(t)) == t);
        CHECK(l_cmpl(r_cmpl(ps)) == ps);
        CHECK(l_cmpl(r_cmpl(M_e)) == M_e);
        CHECK(l_cmpl(r_cmpl(M_u)) == M_u);
        CHECK(l_cmpl(r_cmpl(M)) == M);

        // Also test r_cmpl(l_cmpl(u)) = u
        CHECK(r_cmpl(l_cmpl(s)) == s);
        CHECK(r_cmpl(l_cmpl(v)) == v);
        CHECK(r_cmpl(l_cmpl(B)) == B);
        CHECK(r_cmpl(l_cmpl(t)) == t);
        CHECK(r_cmpl(l_cmpl(ps)) == ps);
        CHECK(r_cmpl(l_cmpl(M_e)) == M_e);
        CHECK(r_cmpl(l_cmpl(M)) == M);

        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: axis definitions")
    {
        fmt::println("G<3,0,1>: axis definitions");

        CHECK(wdg(O_3dp, e1_3dp) == x_axis_3dp);
        CHECK(wdg(O_3dp, e2_3dp) == y_axis_3dp);
        CHECK(wdg(O_3dp, e3_3dp) == z_axis_3dp);
        CHECK(wdg(e2_3dp, e3_3dp) == e23_3dp);
        CHECK(wdg(e3_3dp, e1_3dp) == e31_3dp);
        CHECK(wdg(e1_3dp, e2_3dp) == e12_3dp);

        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: equivalence tests")
    {
        fmt::println("G<3,0,1>: equivalence tests");

        // auto const s = scalar3dp(2);
        auto const v = vec3dp(1, 2, 3, 4);
        auto const B = bivec3dp(5, 6, 7, 8, 9, 10);
        auto const t = trivec3dp(11, 12, 13, 14);
        // auto const ps = pscalar3dp(-2);

        // grp(v,B) -> mvec3dp_u
        CHECK(v * B == (B >> v) + wdg(v, B));
        CHECK(v * B == rwdg(B, r_bulk_dual(v)) + wdg(v, B));
        CHECK(v * B == cmt(v, B) + wdg(v, B)); // only valid for vec*bivec
        // grp(B,v) -> mvec3dp_u
        CHECK(B * v == (v << B) + wdg(B, v));
        CHECK(B * v == rwdg(l_bulk_dual(v), B) + wdg(B, v));
        CHECK(B * v == cmt(B, v) + wdg(B, v));
        // rgpr(v,B) -> mvec3dp_u
        CHECK(rgpr(v, B) == -l_weight_expand3dp(B, v) + rcmt(v, B));
        CHECK(rgpr(v, B) == -wdg(l_weight_dual(B), v) + rcmt(v, B));
        //  rgpr(B,v) -> mvec3dp_u
        CHECK(rgpr(B, v) == -r_weight_expand3dp(v, B) + rcmt(B, v));
        CHECK(rgpr(B, v) == -wdg(v, r_weight_dual(B)) + rcmt(B, v));

        // grp(v,t) -> mvec3dp_e
        CHECK(v * t == (t >> v) + wdg(v, t));
        CHECK(v * t == rwdg(t, r_bulk_dual(v)) + wdg(v, t));

        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: extended metric matrix validation")
    {
        fmt::println("G<3,0,1>: extended metric matrix validation");

        // Get the extended metric matrix from stored constants
        auto G = pga3dp_metric_view();

        // PGA3DP has degenerate metric: e1²=1, e2²=1, e3²=1, e4²=0
        // Basis order: 1, e1, e2, e3, e4, e41, e42, e43, e23, e31, e12, e423, e431, e412,
        // e321, e1234

        // Verify expected diagonal values
        CHECK(G[0, 0] == 1);   // scalar: always 1
        CHECK(G[1, 1] == 1);   // e1: +1 (Euclidean)
        CHECK(G[2, 2] == 1);   // e2: +1 (Euclidean)
        CHECK(G[3, 3] == 1);   // e3: +1 (Euclidean)
        CHECK(G[4, 4] == 0);   // e4: 0 (null/degenerate)
        CHECK(G[5, 5] == 0);   // e41: involves e4, so 0
        CHECK(G[6, 6] == 0);   // e42: involves e4, so 0
        CHECK(G[7, 7] == 0);   // e43: involves e4, so 0
        CHECK(G[8, 8] == 1);   // e23: both Euclidean, so +1
        CHECK(G[9, 9] == 1);   // e31: both Euclidean, so +1
        CHECK(G[10, 10] == 1); // e12: both Euclidean, so +1
        CHECK(G[11, 11] == 0); // e423: involves e4, so 0
        CHECK(G[12, 12] == 0); // e431: involves e4, so 0
        CHECK(G[13, 13] == 0); // e412: involves e4, so 0
        CHECK(G[14, 14] == 1); // e321: all Euclidean, so +1
        CHECK(G[15, 15] == 0); // e1234: involves e4, so 0

        // Verify off-diagonal elements are zero (orthogonal basis)
        for (size_t i = 0; i < 16; ++i) {
            for (size_t j = 0; j < 16; ++j) {
                if (i != j) {
                    CHECK(G[i, j] == 0);
                }
            }
        }

        fmt::println("Extended metric diagonal correct for PGA3DP (degenerate e4)");
    }

    TEST_CASE("G<3,0,1>: extended metric recursive extraction via wedge products")
    {
        fmt::println("G<3,0,1>: extended metric recursive extraction via wedge products");

        auto G = pga3dp_metric_view();

        // Level 0: Scalar (always 1)
        CHECK(G[0, 0] == 1);

        // Level 1: Vectors (extract from dot products)
        // Note: e4 is the null vector (origin in PGA3DP)
        value_t g_e1 = value_t(dot(e1_3dp, e1_3dp));
        value_t g_e2 = value_t(dot(e2_3dp, e2_3dp));
        value_t g_e3 = value_t(dot(e3_3dp, e3_3dp));
        value_t g_e4 = value_t(dot(e4_3dp, e4_3dp)); // Should be 0 (null)
        CHECK(abs(g_e1 - value_t(G[1, 1])) < eps);
        CHECK(abs(g_e2 - value_t(G[2, 2])) < eps);
        CHECK(abs(g_e3 - value_t(G[3, 3])) < eps);
        CHECK(abs(g_e4 - value_t(G[4, 4])) < eps);

        // Level 2: Bivectors (extract from wedge products + dot)
        // Basis order: e41, e42, e43, e23, e31, e12
        auto e41_constructed = wdg(e4_3dp, e1_3dp);
        auto e42_constructed = wdg(e4_3dp, e2_3dp);
        auto e43_constructed = wdg(e4_3dp, e3_3dp);
        auto e23_constructed = wdg(e2_3dp, e3_3dp);
        auto e31_constructed = wdg(e3_3dp, e1_3dp);
        auto e12_constructed = wdg(e1_3dp, e2_3dp);

        value_t g_e41 =
            value_t(dot(e41_constructed, e41_constructed)); // Should be 0 (involves e4)
        value_t g_e42 =
            value_t(dot(e42_constructed, e42_constructed)); // Should be 0 (involves e4)
        value_t g_e43 =
            value_t(dot(e43_constructed, e43_constructed)); // Should be 0 (involves e4)
        value_t g_e23 =
            value_t(dot(e23_constructed, e23_constructed)); // Should be 1 (Euclidean)
        value_t g_e31 =
            value_t(dot(e31_constructed, e31_constructed)); // Should be 1 (Euclidean)
        value_t g_e12 =
            value_t(dot(e12_constructed, e12_constructed)); // Should be 1 (Euclidean)

        CHECK(abs(g_e41 - value_t(G[5, 5])) < eps);
        CHECK(abs(g_e42 - value_t(G[6, 6])) < eps);
        CHECK(abs(g_e43 - value_t(G[7, 7])) < eps);
        CHECK(abs(g_e23 - value_t(G[8, 8])) < eps);
        CHECK(abs(g_e31 - value_t(G[9, 9])) < eps);
        CHECK(abs(g_e12 - value_t(G[10, 10])) < eps);

        // Level 3: Trivectors (extract from wedge products + dot)
        auto e423_constructed = wdg(wdg(e4_3dp, e2_3dp), e3_3dp);
        auto e431_constructed = wdg(wdg(e4_3dp, e3_3dp), e1_3dp);
        auto e412_constructed = wdg(wdg(e4_3dp, e1_3dp), e2_3dp);
        auto e321_constructed = wdg(wdg(e3_3dp, e2_3dp), e1_3dp);

        value_t g_e423 =
            value_t(dot(e423_constructed, e423_constructed)); // Should be 0 (involves e4)
        value_t g_e431 =
            value_t(dot(e431_constructed, e431_constructed)); // Should be 0 (involves e4)
        value_t g_e412 =
            value_t(dot(e412_constructed, e412_constructed)); // Should be 0 (involves e4)
        value_t g_e321 =
            value_t(dot(e321_constructed, e321_constructed)); // Should be 1 (Euclidean)

        CHECK(abs(g_e423 - value_t(G[11, 11])) < eps);
        CHECK(abs(g_e431 - value_t(G[12, 12])) < eps);
        CHECK(abs(g_e412 - value_t(G[13, 13])) < eps);
        CHECK(abs(g_e321 - value_t(G[14, 14])) < eps);

        // Level 4: Pseudoscalar (extract from wedge products + dot)
        auto e1234_constructed = wdg(wdg(wdg(e1_3dp, e2_3dp), e3_3dp), e4_3dp);
        value_t g_e1234 = value_t(
            dot(e1234_constructed, e1234_constructed)); // Should be 0 (involves e4)
        CHECK(abs(g_e1234 - value_t(G[15, 15])) < eps);

        fmt::println(
            "Recursive extraction: all levels match (including null dimension e4)");
    }

    TEST_CASE("G<3,0,1>: extended metric vs dot products")
    {
        fmt::println("G<3,0,1>: extended metric vs dot products");

        // The extended metric diagonal should match dot(basis, basis)
        auto G = pga3dp_metric_view();

        // Scalars: always 1
        CHECK(G[0, 0] == 1);

        // Vectors: e1·e1=1, e2·e2=1, e3·e3=1, e4·e4=0
        CHECK(abs(value_t(dot(e1_3dp, e1_3dp)) - value_t(G[1, 1])) < eps);
        CHECK(abs(value_t(dot(e2_3dp, e2_3dp)) - value_t(G[2, 2])) < eps);
        CHECK(abs(value_t(dot(e3_3dp, e3_3dp)) - value_t(G[3, 3])) < eps);
        CHECK(abs(value_t(dot(e4_3dp, e4_3dp)) - value_t(G[4, 4])) < eps);

        // Bivectors: e41·e41=0, e42·e42=0, e43·e43=0 (involve e4)
        //            e23·e23=1, e31·e31=1, e12·e12=1 (all Euclidean)
        CHECK(abs(value_t(dot(e41_3dp, e41_3dp)) - value_t(G[5, 5])) < eps);
        CHECK(abs(value_t(dot(e42_3dp, e42_3dp)) - value_t(G[6, 6])) < eps);
        CHECK(abs(value_t(dot(e43_3dp, e43_3dp)) - value_t(G[7, 7])) < eps);
        CHECK(abs(value_t(dot(e23_3dp, e23_3dp)) - value_t(G[8, 8])) < eps);
        CHECK(abs(value_t(dot(e31_3dp, e31_3dp)) - value_t(G[9, 9])) < eps);
        CHECK(abs(value_t(dot(e12_3dp, e12_3dp)) - value_t(G[10, 10])) < eps);

        // Trivectors: e423·e423=0, e431·e431=0, e412·e412=0 (involve e4)
        //             e321·e321=1 (all Euclidean)
        CHECK(abs(value_t(dot(e423_3dp, e423_3dp)) - value_t(G[11, 11])) < eps);
        CHECK(abs(value_t(dot(e431_3dp, e431_3dp)) - value_t(G[12, 12])) < eps);
        CHECK(abs(value_t(dot(e412_3dp, e412_3dp)) - value_t(G[13, 13])) < eps);
        CHECK(abs(value_t(dot(e321_3dp, e321_3dp)) - value_t(G[14, 14])) < eps);

        // Pseudoscalar: e1234·e1234=0 (involves null dimension e4)
        CHECK(abs(value_t(dot(e1234_3dp, e1234_3dp)) - value_t(G[15, 15])) < eps);

        fmt::println("Extended metric diagonal matches dot products (including null e4)");
    }

    TEST_CASE("G<3,0,1>: bivec3dp with vec3d lv, lm")
    {
        fmt::println("G<3,0,1>: bivec3dp with vec3d lv, lm");

        auto lv_in = vec3d{1, 1, 1};
        auto lm_in = vec3d{-1, 2, -1};

        auto l = bivec3dp(lv_in, lm_in); // ctor bivev3dp from two vec3d

        // get vector component as vec3d
        auto lv = l.lv();
        auto lm = l.lm();
        auto lm_as_bivec = hd::ga::ega::cmpl(lm);

        CHECK(lv_in == lv);
        CHECK(lm_in == lm);

        fmt::println("l           = {}", l);
        fmt::println("lv          = {}", lv);
        fmt::println("lm          = {}", lm);
        fmt::println("lm as bivec = {}", lm_as_bivec);
        fmt::println("dot(lv,lm)  = {}", hd::ga::ega::dot(lv, lm));
        fmt::println("");
    }

    TEST_CASE("G<3,0,1>: exponential function")
    {
        fmt::println("G<3,0,1>: exponential function");

        //////////////////////////
        // a) pure translation
        //////////////////////////
        fmt::println("\npure translation:\n");

        auto X0 = vec3dp{1, 0, 0, 1};
        auto delta = vec3dp{-2, -1, 3, 0};
        auto X = X0 + delta;

        auto B = att(r_bulk_dual(delta)); // bivector representing delta
        auto M = get_motor(delta);        // translation motor

        fmt::println("X0                   = {}", X0);
        fmt::println("delta                = {}", delta);
        fmt::println("");
        fmt::println("B_tra                = {}", B);
        fmt::println("M (from get_motor()) = {}", M);
        fmt::println("M (from exp(0.5 * B)) = {}", exp(0.5 * B));
        fmt::println("");
        fmt::println("X = M ⟇ X0 ⟇ rrev(m) = {}", move3dp(X0, M));
        fmt::println("X (target value)     = {}", X);
        fmt::println("");

        CHECK(X == move3dp(X0, M)); // motor moves to X
        CHECK(M == exp(0.5 * B));   // bivector creates motor via exp()

        //////////////////////////
        // b) pure rotation
        //////////////////////////
        fmt::println("\npure rotation:\n");

        auto Q = vec3dp{1, 0, 1, 1};
        double phi = deg2rad(90);
        auto phi_vec = phi * bulk_normalize(vec3dp{0, 1, 0, 0});
        auto l = wdg(Q, phi_vec);
        auto l_hat = unitize(l);

        X0 = vec3dp{2, 0, 1, 1};
        X = vec3dp{1, 0, 0, 1};

        B = l_hat * phi;
        M = get_motor(l, phi);

        fmt::println("X0                    = {}", X0);
        fmt::println("l                     = {}", l);
        fmt::println("l_hat                 = {}", l_hat);
        fmt::println("phi                   = {}", phi);
        fmt::println("");
        fmt::println("B_rot                 = {}", B);
        fmt::println("M (from get_motor())  = {}", M);
        fmt::println("M (from exp(0.5 * B)) = {}", exp(0.5 * B));
        fmt::println("");
        fmt::println("X = M ⟇ X0 ⟇ rrev(m)  = {}", move3dp(X0, M));
        fmt::println("X (target value)      = {}", X);
        fmt::println("");

        CHECK(X == move3dp(X0, M));
        CHECK(M == exp(0.5 * B)); // bivector creates motor via exp()

        //////////////////////////
        // c) screw motion (rotation + translation in direction of rotation axis)
        //////////////////////////
        fmt::println("\nscrew motion (rotation and translation):\n");

        // rotational part with axis l
        Q = vec3dp{1, 0, 1, 1};
        phi = deg2rad(90);
        auto phi_hat = bulk_normalize(vec3dp{0, 1, 0, 0});
        l = wdg(Q, phi_vec);
        l_hat = unitize(l);
        auto B_rot = l_hat * phi;

        // translational part with translation along l
        double dist = 1.0;
        delta = dist * phi_hat;
        auto B_tra = att(r_bulk_dual(delta)); // bivector representing delta

        // resulting bivector for screw motion
        B = B_rot + B_tra;

        // start and target point
        X0 = vec3dp{2, 0, 1, 1};
        X = vec3dp{1, 1, 0, 1};

        M = get_motor(l, phi, dist);

        fmt::println("X0                    = {}", X0);
        fmt::println("l                     = {}", l);
        fmt::println("l_hat                 = {}", l_hat);
        fmt::println("phi                   = {}", phi);
        fmt::println("dist                  = {}", dist);
        fmt::println("");
        fmt::println("B_rot                 = {}", B_rot);
        fmt::println("B_tra                 = {}", B_tra);
        fmt::println("B = B_rot + B_tra     = {}", B);
        fmt::println("M (from get_motor())  = {}", M);
        fmt::println("M (from exp(0.5 * B)) = {}", exp(0.5 * B));
        fmt::println("");
        fmt::println("X = M ⟇ X0 ⟇ rrev(m)  = {}", move3dp(X0, M));
        fmt::println("X (target value)      = {}", X);
        fmt::println("");

        CHECK(X == move3dp(X0, M));
        CHECK(M == exp(0.5 * B)); // bivector creates motor via exp()


        fmt::println("exp check:");
        fmt::println("");
        fmt::println("e41_3dp * pi = {}, bulk_nrm_sq = {}, weight_nrm_sq = {}",
                     e41_3dp * pi, bulk_nrm_sq(e41_3dp * pi),
                     weight_nrm_sq(e41_3dp * pi));
        fmt::println("e42_3dp * pi = {}, bulk_nrm_sq = {}, weight_nrm_sq = {}",
                     e42_3dp * pi, bulk_nrm_sq(e42_3dp * pi),
                     weight_nrm_sq(e42_3dp * pi));
        fmt::println("e43_3dp * pi = {}, bulk_nrm_sq = {}, weight_nrm_sq = {}",
                     e43_3dp * pi, bulk_nrm_sq(e43_3dp * pi),
                     weight_nrm_sq(e43_3dp * pi));
        fmt::println("e23_3dp * pi = {}, bulk_nrm_sq = {}, weight_nrm_sq = {}",
                     e23_3dp * pi, bulk_nrm_sq(e23_3dp * pi),
                     weight_nrm_sq(e23_3dp * pi));
        fmt::println("e31_3dp * pi = {}, bulk_nrm_sq = {}, weight_nrm_sq = {}",
                     e31_3dp * pi, bulk_nrm_sq(e31_3dp * pi),
                     weight_nrm_sq(e31_3dp * pi));
        fmt::println("e12_3dp * pi = {}, bulk_nrm_sq = {}, weight_nrm_sq = {}",
                     e12_3dp * pi, bulk_nrm_sq(e12_3dp * pi),
                     weight_nrm_sq(e12_3dp * pi));
        fmt::println("");
        fmt::println("The weight carries rotational information. "
                     "The bivector squares to -1 w.r.t. rgpr().");
        fmt::println("rgpr(e41_3dp,e41_3dp) = {}", rgpr(e41_3dp, e41_3dp));
        fmt::println("rgpr(e42_3dp,e42_3dp) = {}", rgpr(e42_3dp, e42_3dp));
        fmt::println("rgpr(e43_3dp,e43_3dp) = {}", rgpr(e43_3dp, e43_3dp));
        fmt::println("exp(e41_3dp * pi) = {}", exp(e41_3dp * pi));
        fmt::println("exp(e42_3dp * pi) = {}", exp(e42_3dp * pi));
        fmt::println("exp(e43_3dp * pi) = {}", exp(e43_3dp * pi));
        fmt::println("");

        CHECK(exp(e41_3dp * pi) == mvec3dp_e{pscalar3dp{-1}});
        CHECK(exp(e42_3dp * pi) == mvec3dp_e{pscalar3dp{-1}});
        CHECK(exp(e43_3dp * pi) == mvec3dp_e{pscalar3dp{-1}});

        fmt::println("The bulk carries translational information. "
                     "The bivector squares to 0 w.r.t. rgpr().");
        fmt::println("rgpr(e23_3dp,e23_3dp) = {}", rgpr(e23_3dp, e23_3dp));
        fmt::println("rgpr(e31_3dp,e31_3dp) = {}", rgpr(e31_3dp, e31_3dp));
        fmt::println("rgpr(e12_3dp,e12_3dp) = {}", rgpr(e12_3dp, e12_3dp));
        fmt::println("exp(e23_3dp * pi) = {}", exp(e23_3dp * pi));
        fmt::println("exp(e31_3dp * pi) = {}", exp(e31_3dp * pi));
        fmt::println("exp(e12_3dp * pi) = {}", exp(e12_3dp * pi));
        fmt::println("");

        CHECK(exp(e23_3dp * pi) == mvec3dp_e{scalar3dp{0}, e23_3dp * pi, pscalar3dp{1}});
        CHECK(exp(e31_3dp * pi) == mvec3dp_e{scalar3dp{0}, e31_3dp * pi, pscalar3dp{1}});
        CHECK(exp(e12_3dp * pi) == mvec3dp_e{scalar3dp{0}, e12_3dp * pi, pscalar3dp{1}});
    }

    TEST_CASE("G<3,0,1>: sqrt(motor) function")
    {
        fmt::println("G<3,0,1>: sqrt(motor) function");

        // how to transfer plane pl1 into plane pl2, such that pl2 = M ⟇ pl1 ⟇ rrev(M)?
        //
        // We know R = pl2 ⟇ pl1 represents two consequtive reflections across pl1 & pl2.
        // So we need "half of the transformation" M = R^{1/2} = sqrt(R)

        // case a) rotation as defined by consequtive mirroring across intersecting planes
        fmt::println("\nG<3,0,1>: case a) rotation:\n");

        auto pl1 =
            unitize(wdg(vec3dp{0, 0, 0, 1}, wdg(vec3dp{0, 0, 1, 1}, vec3dp{1, 2, 0, 1})));
        auto pl2 =
            unitize(wdg(vec3dp{0, 0, 0, 1}, wdg(vec3dp{0, 0, 1, 1}, vec3dp{3, 2, 0, 1})));

        auto M = get_motor_from_planes(pl1, pl2);

        auto P0 = vec3dp{0, 7, 0, 1};

        auto phi = angle(pl1, pl2);
        auto l_fix = unitize(rwdg(pl1, pl2));

        auto R = sqrt(M);

        fmt::println("pl1             = {}", pl1);
        fmt::println("pl1 ⟇ pl1       = {}", rgpr(pl1, pl1));
        fmt::println("pl1 ⟇ rrev(l1)  = {}", rgpr(pl1, rrev(pl1)));
        fmt::println("pl2             = {}", pl2);
        fmt::println("pl2 ⟇ pl2       = {}", rgpr(pl2, pl2));
        fmt::println("pl2 ⟇ rrev(pl2) = {}", rgpr(pl2, rrev(pl2)));
        fmt::println("");
        fmt::println("l_fix = {}, phi = {}°", l_fix, rad2deg(phi));
        fmt::println("l_fix ⟇ l_fix = {}", rgpr(l_fix, l_fix));
        fmt::println("");
        fmt::println("M(pl1,pl2)  = {}", M);
        fmt::println("M ⟇ rrev(M) = {}", rgpr(M, rrev(M)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = M ⟇ P0 ⟇ rrev(M) = {}", move3dp(P0, M));
        fmt::println("");
        fmt::println("R = sqrt(M) = {}", R);
        fmt::println("R ⟇ rrev(R) = {}", rgpr(R, rrev(R)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = R ⟇ P0 ⟇ rrev(R) = {}", move3dp(P0, R));
        fmt::println("");
        fmt::println("pl1                   = {}", pl1);
        fmt::println("att(pl1)              = {}", att(pl1));
        fmt::println("l_weight_dual(pl1) = {}", l_weight_dual(pl1));
        fmt::println("");

        // case b) translation as defined by parallel planes
        fmt::println("\nG<3,0,1>: case b) translation:\n");

        pl1 =
            unitize(wdg(vec3dp{1, 0, 0, 1}, wdg(vec3dp{1, 0, 1, 1}, vec3dp{0, 0, 0, 1})));
        pl2 =
            unitize(wdg(vec3dp{1, 1, 0, 1}, wdg(vec3dp{1, 1, 1, 1}, vec3dp{0, 1, 0, 1})));

        M = get_motor_from_planes(pl1, pl2);

        P0 = vec3dp{0, 7, 0, 1};

        R = sqrt(M);

        fmt::println("pl1             = {}", pl1);
        fmt::println("pl1 ⟇ pl1       = {}", rgpr(pl1, pl1));
        fmt::println("pl1 ⟇ rrev(l1)  = {}", rgpr(pl1, rrev(pl1)));
        fmt::println("pl2             = {}", pl2);
        fmt::println("pl2 ⟇ pl2       = {}", rgpr(pl2, pl2));
        fmt::println("pl2 ⟇ rrev(pl2) = {}", rgpr(pl2, rrev(pl2)));
        fmt::println("");
        fmt::println("M(pl1,pl2)  = {}", M);
        fmt::println("M ⟇ rrev(M) = {}", rgpr(M, rrev(M)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = M ⟇ P0 ⟇ rrev(M) = {}", move3dp(P0, M));
        fmt::println("");
        fmt::println("R = sqrt(M) = {}", R);
        fmt::println("R ⟇ rrev(R) = {}", rgpr(R, rrev(R)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = R ⟇ P0 ⟇ rrev(R) = {}", move3dp(P0, R));
        fmt::println("");
        fmt::println("pl1                   = {}", pl1);
        fmt::println("att(pl1)              = {}", att(pl1));
        fmt::println("l_weight_dual(pl1) = {}", l_weight_dual(pl1));
        fmt::println("");

        // case c) combined rotation and translation as defined by transforming
        //         one line l1 into another line l2, if l1 and l2 are not intersecting
        fmt::println("\nG<3,0,1>: case c) combined rotation and translation:\n");

        auto l1 = x_axis_3dp;
        auto l2 = wdg(vec3dp{0, 1, 0, 1}, vec3dp{0, 1, -1, 1});
        auto l3 = rcmt(l2, rrev(l1)); // rotation axis (shortest dist between l1, l2)

        M = get_motor(l3, deg2rad(90), 1.0);
        auto M_double = rgpr(l2, rrev(l1));
        auto M_calc = sqrt(M_double);
        auto M2 = get_motor_from_lines(l1, l2);

        P0 = vec3dp{1, 0, 0, 1};
        auto P1 = vec3dp{0, 1, -1, 1};

        fmt::println("l1 = {}", l1);
        fmt::println("l2 = {}", l2);
        fmt::println("l3 = {}", l3);
        fmt::println("");
        fmt::println("M = {}", M);
        fmt::println("M_double = {}", M_double);
        fmt::println("M_calc = {}", M_calc);
        fmt::println("M2 = {}", M2);

        CHECK(l3 == y_axis_3dp);
        CHECK(P1 == move3dp(P0, M));
        CHECK(P1 == move3dp(P0, M2));

        fmt::println("l1            = {}", l1);
        fmt::println("l1 ⟇ l1       = {}", rgpr(l1, l1));
        fmt::println("l1 ⟇ rrev(l1) = {}", rgpr(l1, rrev(l1)));
        fmt::println("l2            = {}", l2);
        fmt::println("l2 ⟇ l2       = {}", rgpr(l2, l2));
        fmt::println("l2 ⟇ rrev(l2) = {}", rgpr(l2, rrev(l2)));
        fmt::println("");
        fmt::println("M(l1,l2)    = {}", M);
        fmt::println("M ⟇ rrev(M) = {}", rgpr(M, rrev(M)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = M ⟇ P0 ⟇ rrev(M) = {}", move3dp(P0, M));
        fmt::println("");
        fmt::println("M2                      = {}", M2);
        fmt::println("M2 ⟇ rrev(M2)           = {}", rgpr(M2, rrev(M2)));
        fmt::println("");
        fmt::println("P0                      = {}", P0);
        fmt::println("P  = M2 ⟇ P0 ⟇ rrev(M2) = {}", move3dp(P0, M2));
        fmt::println("");
    }

    TEST_CASE("Vec3dp: operations - inverses II")
    {
        fmt::println("Vec3dp: operations - inverses II");

        vec3dp v0{2, 1, 2, 0};
        vec3dp v1{2, 1, 2, 1};
        vec3dp v2{2, 1, 2, 2};

        fmt::println("v0 = {}, inv(v0) = {}, , rinv(v0) = n.a.", v0, inv(v0));
        fmt::println("v1 = {}, inv(v1) = {}, , rinv(v1) = v1 / (-v1.w * v1.w) = {}", v1,
                     inv(v1), v1 / (-v1.w * v1.w));
        fmt::println("v2 = {}, inv(v2) = {}, , rinv(v2) = v2 / (-v2.w * v2.w) = {}", v2,
                     inv(v2), v2 / (-v2.w * v2.w));
        fmt::println("");
        fmt::println("dot(v0,v0) = {}, rdot(v0,v0) = {}", dot(v0, v0), rdot(v0, v0));
        fmt::println("dot(v1,v1) = {}, rdot(v1,v1) = {}", dot(v1, v1), rdot(v1, v1));
        fmt::println("dot(v2,v2) = {}, rdot(v2,v2) = {}", dot(v2, v2), rdot(v2, v2));
        fmt::println("");
        fmt::println("v0 ⟑ v0 = {}, v0 ⟑ inv(v0) = {}", v0 * v0, v0 * inv(v0));
        fmt::println("v1 ⟑ v1 = {}, v1 ⟑ inv(v1) = {}", v1 * v1, v1 * inv(v1));
        fmt::println("v2 ⟑ v2 = {}, v2 ⟑ inv(v2) = {}", v2 * v2, v2 * inv(v2));
        fmt::println("");
        fmt::println("v0 ⟇ v0 = {}, v0 ⟇ rinv(v0) = n.a.", rgpr(v0, v0));
        fmt::println("v1 ⟇ v1 = {}, v1 ⟇ rinv(v1) = {}", rgpr(v1, v1),
                     rgpr(v1, v1 / (-v1.w * v1.w)));
        fmt::println("v2 ⟇ v2 = {}, v2 ⟇ rinv(v2) = {}", rgpr(v2, v2),
                     rgpr(v2, v2 / (-v2.w * v2.w)));
        fmt::println("");
    }


} // PGA 3DP Tests

// | ⟑ | U+27D1 | (direct Unicode) | Geometric product |
// | ⟇ | U+27C7 | (direct Unicode) | Regressive geometric product |
