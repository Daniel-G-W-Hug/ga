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

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// projective geometric algebra 2d: pga2dp (embedded in a 3d representational space)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA 2DP Tests")
{
    TEST_CASE("G<2,0,1>: pga2dp")
    {
        fmt::println("G<2, 0, 1>: pga2dp");
        // 2d projective geometric algebra (2dp) models two-dimensional Euclidean algebra
        // embedding it in a three-dimensional representational space
        const algebra<2, 0, 1> alg;
        CHECK(alg.p() == 2);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("   pga2dp: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga2dp: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Vec2dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2dp: default init")
    {
        fmt::println("Vec2dp: default init");
        vec2dp v;
        // fmt::println("   v = {}", v);
        CHECK(abs(v.x) < eps);
        CHECK(abs(v.y) < eps);
        CHECK(abs(v.z) < eps);
    }

    TEST_CASE("Vec2dp: with curly braced intializer")
    {
        fmt::println("Vec2dp: with curly braced intializer");
        vec2dp v{0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(abs(v.x) < eps);
        CHECK(abs(v.y) < eps);
        CHECK(abs(v.z) < eps);
    }

    TEST_CASE("Vec2dp: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec2dp: cp ctor & cp assign incl. type deduction");
        vec2dp v1{1.0, 2.0, 3.0}; // init with double (type deduction)
        vec2dp v2{v1};            // cp ctor
        vec2dp v3 = v2;           // cp assign
        vec2dp v4 = -v2;          // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(abs(v1.x - 1.0) < eps);
        CHECK(abs(v1.y - 2.0) < eps);
        CHECK(abs(v1.z - 3.0) < eps);
        CHECK(abs(v2.x - 1.0) < eps);
        CHECK(abs(v2.y - 2.0) < eps);
        CHECK(abs(v2.z - 3.0) < eps);
        CHECK(abs(v3.x - 1.0) < eps);
        CHECK(abs(v3.y - 2.0) < eps);
        CHECK(abs(v3.z - 3.0) < eps);
        CHECK(v4 == -v2);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }


    TEST_CASE("Vec2dp: fmt & cout printing")
    {
        fmt::println("Vec2dp: fmt & cout printing");

        vec2dp pf{1.0f, 2.00001f, 3.0f};
        vec2dp pd{1.0, 2.00001, 3.0};

        fmt::println("");
        // std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.8f}", pf);

        // std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.8f}", pd);

        std::vector<Vec2dp<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:e}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("Vec2dp: comparison float")
    {
        fmt::println("Vec2dp: comparison float");

        Vec2dp<float> v1f{1.0f, 2.0f, 3.0f};
        Vec2dp<float> v2f{2.0f, 4.0f, 3.0f};
        Vec2dp<float> v3f{1.0f, 2.0000001f, 3.0f};
        Vec2dp<float> v4f{v1f};

        // fmt::println("   v1f = {}", v1f);
        // fmt::println("   v2f = {}", v2f);
        // fmt::println("   v3f = {}", v3f);
        // fmt::println("   v4f = {}", v4f);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<float>::epsilon());

        CHECK(v1f == v4f);           // comparison (equality)
        CHECK(v1f != v2f);           // comparison (inequality)
        CHECK(nrm(v1f) < nrm(v2f));  // comparison (less than)
        CHECK(nrm(v2f) >= nrm(v1f)); // comparison (greater than or equal)
        CHECK(v3f == v1f);           // comparison (eqality)
    }

    TEST_CASE("Vec2dp: comparison double")
    {
        fmt::println("Vec2dp: comparison double");

        Vec2dp<double> v1d{1.0, 2.0, 3.0};
        Vec2dp<double> v2d{2.0, 4.0, 3.0};
        Vec2dp<double> v3d{1.0, 2.0000000000000001, 3.0};
        Vec2dp<double> v4d{v1d};

        // fmt::println("   v1d = {}", v1d);
        // fmt::println("   v2d = {}", v2d);
        // fmt::println("   v3d = {}", v3d);
        // fmt::println("   v4d = {}", v4d);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<double>::epsilon());

        CHECK(v1d == v4d);           // comparison (equality)
        CHECK(v1d != v2d);           // comparison (inequality)
        CHECK(nrm(v1d) < nrm(v2d));  // comparison norm
        CHECK(nrm(v2d) >= nrm(v1d)); // comparison norm
        CHECK(v3d == v1d);           // comparison (eqality)
    }

    TEST_CASE("Vec2dp: vector space and linearity tests")
    {
        fmt::println("Vec2dp: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        vec2dp p0;
        vec2dp p1{1.0, 2.0, 3.0};
        vec2dp p2{2.0, 4.0, 6.0};
        vec2dp p3{3.0, 6.0, 9.0};
        vec2dp p4 = -p1; // assignment using unary minus
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

    TEST_CASE("Vec2dp: inner product properties")
    {
        fmt::println("Vec2dp: inner product properties");

        double a = 2.35;
        vec2dp u{1.0, 2.0, 1.0};
        vec2dp v{-0.5, 3.0, 0.5};
        vec2dp w{3.0, 6.0, -3.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    TEST_CASE("Vec2dp: linking inner and outer product")
    {
        fmt::println("Vec2dp: linking inner and outer product");

        auto s1 = scalar2dp{2.0};
        auto s2 = scalar2dp{3.0};
        auto v1 = vec2dp{1.0, 2.0, 3.0};
        auto v2 = vec2dp{-1.0, 4.0, -5.0};
        auto B1 = bivec2dp{4.0, 5.0, 6.0};
        auto B2 = bivec2dp{-4.0, -5.0, -6.0};
        auto ps1 = pscalar2dp{-2.0};
        auto ps2 = pscalar2dp{4.0};

        CHECK(bulk_dual(wdg(v1, B1)) == rwdg(bulk_dual(v1), bulk_dual(B1)));

        fmt::println("");
        fmt::println("wdg(s1, bulk_dual(s2)) = {}", wdg(s1, bulk_dual(s2)));
        fmt::println("dot(s1, s2) ^ I_2dp = {}", wdg(dot(s1, s2), I_2dp));
        fmt::println("");
        fmt::println("wdg(v1, bulk_dual(v2)) = {}", wdg(v1, bulk_dual(v2)));
        fmt::println("dot(v1, v2) ^ I_2dp = {}", wdg(dot(v1, v2), I_2dp));
        fmt::println("");
        fmt::println("wdg(B1, bulk_dual(B2)) = {}", wdg(B1, bulk_dual(B2)));
        fmt::println("dot(B1, B2) ^ I_2dp = {}", wdg(dot(B1, B2), I_2dp));
        fmt::println("");
        fmt::println("wdg(ps1, bulk_dual(ps2)) = {}", wdg(ps1, bulk_dual(ps2)));
        fmt::println("dot(ps1, ps2) ^ I_2dp = {}", wdg(dot(ps1, ps2), I_2dp));
        fmt::println("");

        CHECK(wdg(s1, bulk_dual(s2)) == wdg(dot(s1, s2), I_2dp));
        //
        CHECK(wdg(v1, bulk_dual(v2)) == wdg(dot(v1, v2), I_2dp));
        //
        CHECK(wdg(B1, bulk_dual(B2)) == wdg(dot(B1, B2), I_2dp));
        //
        CHECK(wdg(ps1, bulk_dual(ps2)) == wdg(dot(ps1, ps2), I_2dp));

        CHECK(wdg(v1, weight_dual(v2)) == rdot(v1, v2));

        fmt::println("");


        CHECK(weight_dual(rwdg(v1, B1)) == wdg(weight_dual(v1), weight_dual(B1)));

        CHECK(rwdg(v1, bulk_dual(v2)) == dot(v1, v2));

        fmt::println("");
        fmt::println("rwdg(s1, weight_dual(s2)) = {}", rwdg(s1, weight_dual(s2)));
        fmt::println("rwdg(rdot(s1, s2), I_2dp) = {}",
                     rwdg(rdot(s1, s2), scalar2dp{1.0}));
        fmt::println("");
        fmt::println("rwdg(v1, weight_dual(v2)) = {}", rwdg(v1, weight_dual(v2)));
        fmt::println("rwdg(rdot(v1, v2), I_2dp) = {}",
                     rwdg(rdot(v1, v2), scalar2dp{1.0}));
        fmt::println("");
        fmt::println("rwdg(B1, weight_dual(B2)) = {}", rwdg(B1, weight_dual(B2)));
        fmt::println("rwdg(rdot(B1, B2), I_2dp) = {}",
                     rwdg(rdot(B1, B2), scalar2dp{1.0}));
        fmt::println("");
        fmt::println("rwdg(ps1, weight_dual(ps2)) = {}", rwdg(ps1, weight_dual(ps2)));
        fmt::println("rwdg(rdot(ps1, ps2), I_2dp) = {}",
                     rwdg(rdot(ps1, ps2), scalar2dp{1.0}));
        fmt::println("");

        CHECK(rwdg(s1, weight_dual(s2)) == rwdg(rdot(s1, s2), scalar2dp{1.0}));
        //
        CHECK(rwdg(v1, weight_dual(v2)) == rwdg(rdot(v1, v2), scalar2dp{1.0}));
        //
        CHECK(rwdg(B1, weight_dual(B2)) == rwdg(rdot(B1, B2), scalar2dp{1.0}));
        //
        CHECK(rwdg(ps1, weight_dual(ps2)) == rwdg(rdot(ps1, ps2), scalar2dp{1.0}));

        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2dp<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2dp: operations - norm, inverse, dot")
    {
        fmt::println("Vec2dp: operations - norm, inverse, dot");

        scalar2dp s1{3.2};

        vec2dp v1{2.0, 1.0, 2.0};
        vec2dp v2{normalize(v1)};

        vec2dp v3{2.0, 6.0, -4.0};
        vec2dp v4{inv(v3)};

        bivec2dp b1{-2.3, 1.2, 4.5};
        pscalar2dp ps1{-4.7};
        mvec2dp_e mve1{s1, b1};
        mvec2dp_u mvu1{v1, ps1};
        mvec2dp mv1{s1, v1, b1, ps1};

        // auto v3m = mvec2dp{v3};
        // auto v4m = mvec2dp{v4};
        // fmt::println("v1 = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2 = normalize(v1) = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("v3 = {: .4f}, nrm(v1) = {: .4f}", v3, nrm(v3));
        // fmt::println(
        //     "v4 = inv(v3) = {: .4f}, nrm(v3) = {: .4f}, nrm(v3)*nrm(v4) = {: .4f}",
        //     v4, nrm(v4), nrm(v3) * nrm(v4));
        // fmt::println("");
        // fmt::println("v3m = {}, nrm(v3m) = {}", v3m, nrm(v3m));
        // fmt::println("v4m = {}, nrm(v4m) = {}", v4m, nrm(v4m));
        // fmt::println("0.5*(v3m*v4m + v4m*v3m) = {}", 0.5 * (v3m * v4m + v4m * v3m));
        // fmt::println("0.5*(v3m*v4m - v4m*v3m) = {}", 0.5 * (v3m * v4m - v4m * v3m));

        CHECK(abs(nrm_sq(v1) - 9.0) < eps);
        CHECK(abs(bulk_nrm_sq(v1) - 5.0) < eps);
        CHECK(abs(nrm_sq(v2) - 1.0) < eps);
        CHECK(abs(bulk_nrm_sq(v3) - 40.0) < eps);
        CHECK(abs(dot(v4, v3) - 1.0) < eps);

        // check inverses - scalar
        // fmt::println("");
        // fmt::println("s1 * inv(s1) = {}", s1 * inv(s1)); // s
        CHECK(abs(nrm(s1 * inv(s1)) - 1) < eps);
        CHECK(abs(inv(s1) - rev(s1) / bulk_nrm_sq(s1)) < eps);

        // check inverses - vector
        // fmt::println("v1 * inv(v1) = {}", v1 * inv(v1)); // mv_e
        CHECK(abs(nrm(gr0(v1 * inv(v1))) - 1) < eps);
        CHECK(abs(nrm(gr2(v1 * inv(v1))) - 0) < eps);
        CHECK(abs(nrm(inv(v1) - rev(v1) / bulk_nrm_sq(v1))) < eps);

        // check inverses - bivector
        // fmt::println("b1 * inv(b1) = {}", b1 * inv(b1)); // mv_e
        CHECK(abs(nrm(gr0(b1 * inv(b1))) - 1) < eps);
        CHECK(abs(nrm(gr2(b1 * inv(b1))) - 0) < eps);
        CHECK(abs(nrm(inv(b1) - rev(b1) / bulk_nrm_sq(b1))) < eps);

        // check inverses - pseudoscalar
        // due to the degenerate metric there is no inverse of the pseudoscalar

        // check inverses - even grade multivector
        // fmt::println("mve1 * inv(mve1) = {}", mve1 * inv(mve1)); // mv_e
        CHECK(abs(nrm(gr0(mve1 * inv(mve1))) - 1) < eps);
        CHECK(abs(nrm(gr2(mve1 * inv(mve1))) - 0) < eps);
        CHECK(abs(nrm(inv(mve1) - rev(mve1) / bulk_nrm_sq(mve1))) < eps);

        // check inverses - odd grade multivector
        // fmt::println("mvu1 * inv(mvu1) = {}", mvu1 * inv(mvu1)); // mv_e
        CHECK(abs(nrm(gr0(mvu1 * inv(mvu1))) - 1) < eps);
        CHECK(abs(nrm(gr2(mvu1 * inv(mvu1))) - 0) < eps);
        CHECK(abs(nrm(inv(mvu1) - rev(mvu1) / bulk_nrm_sq(mvu1))) < eps);

        // check inverses - multivector
        // fmt::println("mv1 * inv(mv1) = {}", mv1 * inv(mv1)); // mv
        CHECK(abs(nrm(gr0(mv1 * inv(mv1))) - 1) < eps);
        CHECK(abs(nrm(gr1(mv1 * inv(mv1))) - 0) < eps);
        CHECK(abs(nrm(gr2(mv1 * inv(mv1))) - 0) < eps);
        CHECK(abs(nrm(gr3(mv1 * inv(mv1))) - 0) < eps);
        CHECK(abs(nrm(gr0(inv(mv1) * mv1)) - 1) < eps); // left and right inverse
                                                        // are equal
        // fmt::println("");
    }

    TEST_CASE("Vec2dp: operations - bulk_dual, weight_dual")
    {
        fmt::println("Vec2dp: operations - bulk_dual, weight_dual");

        auto s = scalar2dp{5.0};
        auto v = vec2dp{1.0, 2.0, 1.0};
        auto B = bivec2dp{-1.0, 2.0, 1.0};
        auto ps = pscalar2dp{-5.0};

        auto v2 = vec2dp{3.0, -2.0, 1.0};
        auto B2 = bivec2dp{-3.0, 2.0, 5.0};

        CHECK(bulk_dual(s) == pscalar2dp(to_val(s)));
        CHECK(weight_dual(s) == pscalar2dp{0.0});

        CHECK(bulk_dual(v) == bivec2dp{-v.y, v.x, 0.0});
        CHECK(weight_dual(v) == bivec2dp{0.0, 0.0, -v.z});

        CHECK(bulk_dual(B) == vec2dp{0.0, 0.0, -B.z});
        CHECK(weight_dual(B) == vec2dp{B.y, -B.x, 0.0});

        CHECK(bulk_dual(ps) == scalar2dp{0.0});
        CHECK(weight_dual(ps) == scalar2dp(to_val(ps)));

        // duality of wdg and rwdg based on complements
        CHECK(bulk_dual(wdg(v, v2)) == rwdg(bulk_dual(v), bulk_dual(v2)));
        CHECK(bulk_dual(wdg(v, B)) == rwdg(bulk_dual(v), bulk_dual(B)));
        CHECK(bulk_dual(wdg(B, v)) == rwdg(bulk_dual(B), bulk_dual(v)));

        // contractions
        CHECK(left_bulk_contract2dp(v, v2) == rwdg(bulk_dual(v), v2));
        CHECK(left_bulk_contract2dp(B, B2) == rwdg(bulk_dual(B), B2));
        CHECK(left_bulk_contract2dp(B, v) == rwdg(bulk_dual(B), v));

        CHECK(left_weight_contract2dp(v, v2) == rwdg(weight_dual(v), v2));
        CHECK(left_weight_contract2dp(B, B2) == rwdg(weight_dual(B), B2));
        CHECK(left_weight_contract2dp(B, v) == rwdg(weight_dual(B), v));

        CHECK(right_bulk_contract2dp(v, v2) == rwdg(v, bulk_dual(v2)));
        CHECK(right_bulk_contract2dp(B, B2) == rwdg(B, bulk_dual(B2)));
        CHECK(right_bulk_contract2dp(B, v) ==
              rwdg(B, bulk_dual(v))); // contracts v onto B

        CHECK(right_weight_contract2dp(v, v2) == rwdg(v, weight_dual(v2)));
        CHECK(right_weight_contract2dp(B, B2) == rwdg(B, weight_dual(B2)));
        CHECK(right_weight_contract2dp(B, v) == rwdg(B, weight_dual(v)));

        // expansions
        CHECK(left_weight_expand2dp(v, v2) == wdg(weight_dual(v), v2));
        CHECK(left_weight_expand2dp(B, B2) == wdg(weight_dual(B), B2));
        CHECK(left_weight_expand2dp(v, B) == wdg(weight_dual(v), B));

        CHECK(left_bulk_expand2dp(v, v2) == wdg(bulk_dual(v), v2));
        CHECK(left_bulk_expand2dp(B, B2) == wdg(bulk_dual(B), B2));
        CHECK(left_bulk_expand2dp(v, B) == wdg(bulk_dual(v), B));

        CHECK(right_weight_expand2dp(v, v2) == wdg(v, weight_dual(v2)));
        CHECK(right_weight_expand2dp(B, B2) == wdg(B, weight_dual(B2)));
        CHECK(right_weight_expand2dp(v, B) == wdg(v, weight_dual(B)));

        CHECK(right_bulk_expand2dp(v, v2) == wdg(v, bulk_dual(v2)));
        CHECK(right_bulk_expand2dp(B, B2) == wdg(B, bulk_dual(B2)));
        CHECK(right_bulk_expand2dp(v, B) == wdg(v, bulk_dual(B)));
    }

    TEST_CASE("Vec2dp: operations - angle I")
    {
        fmt::println("Vec2dp: operations - angle I");

        vec2dp v1{1.0, 0.0, 0.0};
        vec2dp v2{normalize(vec2dp(1.0, 1.0, 0.0))};
        vec2dp v3{0.0, 1.0, 0.0};
        vec2dp v4{normalize(vec2dp(-1.0, 1.0, 0.0))};
        vec2dp v5{-1.0, 0.0, 0.0};
        vec2dp v6{normalize(vec2dp(-1.0, -1.0, 0.0))};
        vec2dp v7{0.0, -1.0, 0.0};
        vec2dp v8{normalize(vec2dp(1.0, -1.0, 0.0))};

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
        CHECK(v6 == normalize(vec2dp(-1.0, -1.0, 0.0)));
        CHECK(v7 == normalize(vec2dp(0.0, -1.0, 0.0)));
        CHECK(v8 == normalize(vec2dp(1.0, -1.0, 0.0)));
    }

    TEST_CASE("Vec2dp: operations - angle II")
    {
        fmt::println("Vec2dp: operations - angle II");

        std::vector<std::tuple<double, Vec2dp<double>>> v1;
        std::vector<std::tuple<double, Vec2dp<double>>> v2;
        std::vector<std::tuple<double, Vec2dp<double>>> v3;

        // only positive angles are easy to implement vs. the 2d case

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi), std::sin(phi), 0.0);
            v1.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_2dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2), 0.0);
            v2.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e2_2dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4), 0.0);
            v3.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_2dp + e2_2dp, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(abs(phi - angle(e1_2dp, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(abs(phi - angle(e2_2dp, c)) < eps);
        }
        auto ref_vec = normalize(e1_2dp + e2_2dp);
        for (auto const& [phi, c] : v3) {
            CHECK(abs(phi - angle(ref_vec, c)) < eps);
        }
    }

    // TEST_CASE("Vec2dp: operations - angle III (angle between lines)")
    // {
    //     fmt::println("Vec2dp: operations - angle III (angle between lines)");

    //     auto x_axis = bivec2dp(0, 1, 0);
    //     auto y_axis = bivec2dp(1, 0, 0); // really this is -y_axis_2dp

    //     fmt::println("");

    //     for (int i = 0; i <= 23; ++i) {
    //         double phi_y = i * pi / 12;
    //         double phi_x = phi_y + pi / 2.;
    //         bivec2dp b(cos(phi_x), sin(phi_x),
    //                    1.5); // tangent lines to circle with r = 1.5

    //         auto phi_xc = angle(x_axis, b);
    //         fmt::println("i={:3}, phi_x={: 8.3f}, phi_y={: 8.3f}, phi_xc={: 8.3f}", i,
    //                      rad2deg(phi_x), rad2deg(phi_y), rad2deg(phi_xc));
    //     }

    //     fmt::println("");
    // }

    TEST_CASE("Vec2dp: operations - wedge")
    {
        fmt::println("Vec2dp: operations - wedge");

        vec2dp v1{1.0, 0.0, 0.0};
        vec2dp v2{normalize(vec2dp(1.0, 1.0, 0.0))};
        vec2dp v3{0.0, 1.0, 0.0};
        vec2dp v4{normalize(vec2dp(-1.0, 1.0, 0.0))};
        vec2dp v5{-1.0, 0.0, 0.0};
        vec2dp v6{normalize(vec2dp(-1.0, -1.0, 0.0))};
        vec2dp v7{0.0, -1.0, 0.0};
        vec2dp v8{normalize(vec2dp(1.0, -1.0, 0.0))};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar2dp{sd};
        auto t = scalar2dp{st};

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

        CHECK(wdg(v1, v1) == bivec2dp{});                     // wdg=0 for collin. vectors
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

    TEST_CASE("Vec2dp: operations - project / reject / reflect")
    {
        fmt::println("Vec2dp: operations - project / reject / reflect");

        // vectors (=directions or projective points at infinity)
        vec2dp v1{5, 1, 0};
        vec2dp v2{2, 2, 0};

        // points
        vec2dp p1{5, 1, 1};

        // vectors (=directions) projected and rejected
        vec2dp v3{project_onto(v1, v2)};
        vec2dp v4{reject_from(v1, v2)};
        vec2dp v5{v3 + v4};
        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);

        // points projected and rejected

        // line
        bivec2dp B1{e31_2dp};
        vec2dp p3{project_onto(p1, B1)};
        vec2dp v6{reject_from(p1, B1)};
        vec2dp p8{p3 + v6};
        CHECK(p3 + v6 == p8);
        CHECK(p8 == p1);

        // point reflected on a line
        vec2dp p{4, 1, 1};
        auto B2 = e31_2dp; // B2 = x_axis_2dp, just modifyies y-coordinate
        CHECK(unitize(reflect_on(p, B2)) == vec2dp{4, -1, 1});

        // coordinate axis reflected on perpendicular axis yield their negatives
        CHECK(reflect_on(y_axis_2dp, x_axis_2dp) == -y_axis_2dp);
        CHECK(reflect_on(x_axis_2dp, y_axis_2dp) == -x_axis_2dp);

        // coordinate axis reflected on itself remains itself (identity)
        CHECK(reflect_on(x_axis_2dp, x_axis_2dp) == x_axis_2dp);
        CHECK(reflect_on(y_axis_2dp, y_axis_2dp) == y_axis_2dp);

        // lines parallel to coordinate axis after reflexion:
        // remain parallel, have same orientation, but are on other side of axis
        CHECK(reflect_on(bivec2dp{1, 0, 1}, x_axis_2dp) == bivec2dp{1, 0, -1});
        CHECK(reflect_on(bivec2dp{0, 1, 1}, y_axis_2dp) == bivec2dp{0, 1, -1});
    }

    TEST_CASE("Vec2dp: operations - rotation")
    {
        fmt::println("Vec2dp: operations - rotation");
        {
            std::vector<std::tuple<double, Vec2dp<double>>> v;

            // fmt::println("");
            for (int i = -12; i <= 12; ++i) {
                double phi = i * pi / 12;
                auto c = vec2dp(std::cos(phi), std::sin(phi), 0.0);
                auto d = move2dp(x_dir_2dp, get_motor(O_2dp, phi));
                v.emplace_back(std::make_tuple(phi, c));
                // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
                //              " angle={: .4f}",
                //              i, phi, rad2deg(phi), c,
                //              rad2deg(angle(x_dir_2dp, c)));
                // fmt::println("                                  d={: .3f}", d);
                CHECK(c == d);
            }
        }

        {
            // rotations in a plane (using motors)
            std::vector<std::tuple<double, Vec2dp<double>>> v;

            // fmt::println("");
            for (int i = -12; i <= 12; ++i) {
                double phi = i * pi / 12;
                auto c = vec2dp(std::cos(phi), std::sin(phi), 1.0);
                v.emplace_back(std::make_tuple(phi, c));
                // auto b = cmpl(c); // bivec, representing plane with normal c
                // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
                //              " angle={: .4f}, b={: .3f}",
                //              i, phi, rad2deg(phi), c, angle(e1_2dp, c), b);
                CHECK(c == move2dp(e1_2dp + O_2dp, get_motor(O_2dp, phi)));

                // check rotate optimization (use points on the projective plan)
                CHECK(move2dp(e1_2dp + O_2dp, get_motor(O_2dp, phi)) ==
                      move2dp_opt(e1_2dp + O_2dp, get_motor(O_2dp, phi)));

                // rotation of a ega2dp multivector is to rotate the vector part
                // exclusively
                // the bivector part should remain unchanged for rotations in the e12
                // plane
                auto mv_rev =
                    mvec2dp{scalar2dp{4.0}, e1_2dp + O_2dp, e12_2dp, pscalar2dp{-2.0}};
                auto mv_rot = mvec2dp{scalar2dp{4.0}, c, e12_2dp, pscalar2dp{-2.0}};
                CHECK(mv_rot == move2dp(mv_rev, get_motor(O_2dp, phi)));
            }

            // check same rotation applied to many points
            std::vector<Vec2dp<double>> vec_ref(100);
            std::vector<Vec2dp<double>> vec_rot(100);
            std::vector<Vec2dp<double>> vec_rot_calc(100);

            // prepare for checking vector-based rotation (after the loop)
            double phi = pi / 12;
            for (size_t i = 0; i < 100; ++i) {
                vec_ref.push_back(e1_2dp);
                vec_rot.emplace_back(move2dp(e1_2dp, get_motor(O_2dp, phi)));
            }

            vec_rot_calc = move2dp(vec_ref, get_motor(O_2dp, phi));
            for (size_t i = 0; i < 100; ++i) {
                CHECK(vec_rot_calc[i] == vec_rot[i]);
            }
        }

        fmt::println("");
    }

    TEST_CASE("Vec2dp: operations - translation")
    {
        fmt::println("Vec2dp: operations - translation");

        auto const P = vec2dp{3.5, 0.5, 1}; // point P
        auto const p = vec2dp{3.5, 0.5, 0}; // vector from O to P: p = P-O
        auto const d = bulk_nrm(p);

        auto m_tra = get_motor(p);
        auto Pt = move2dp(P, m_tra);

        fmt::println("m_tra        = {: 5.3f}", m_tra);
        fmt::println("");
        fmt::println("P            = {: 5.3f}", P);
        fmt::println("p            = {: 5.3f}", p);
        fmt::println("Pt           = {: 5.3f}", Pt);
        fmt::println("O_2dp + 2*p  = {: 5.3f}", O_2dp + 2 * p);
        fmt::println("");

        CHECK(Pt == O_2dp + 2.0 * p);

        auto ln = unitize(wdg(O_2dp, P));

        fmt::println("ln                = {: 5.3f}", ln);
        fmt::println("rgpr(ln,ln)       = {: 5.3f}", rgpr(ln, ln));
        fmt::println("rgpr(ln,rrev(ln)) = {: 5.3f}", rgpr(ln, rrev(ln)));
        fmt::println("rgpr(d/2.0*ln,P)  = {: 5.3f}", rgpr(d / 2.0 * ln, P));
        fmt::println("rgpr(P,d/2.0*ln)  = {: 5.3f}", rgpr(P, d / 2.0 * ln));
        fmt::println("rgrp(rgpr(d / 2.0 * ln,P),rrev(d / 2.0 * ln))  = {: 5.3f}",
                     rgpr(rgpr(d / 2.0 * ln, P), rrev(d / 2.0 * ln)));
        fmt::println("unitize(gr1(rgrp(rgpr(d/2.0*ln,P),rrev(d/2.0*ln))))  = {: 5.3f}",
                     unitize(gr1(rgpr(rgpr(d / 2.0 * ln, P), rrev(d / 2.0 * ln)))));

        fmt::println("");
        fmt::println("rgpr(-d/2*ln,d/2*ln)  = {: 5.3f}",
                     rgpr((-d / 2.0) * ln, (d / 2.0) * ln));
        fmt::println("-d*d/4  = {: 5.3f}", -d * d / 4.0);

        auto m_tra_c = rgpr(d / 2.0, ln) + pscalar2dp{1.0};
        fmt::println("m_tra_c  = {: 5.3f}", m_tra_c);

        fmt::println("");
    }


    TEST_CASE("Vec2dp: operations - rotation + translation combined")
    {
        fmt::println("Vec2dp: operations - rotation + translation combined");

        auto pc = vec2dp{-1, -1, 1}; // center or turning point

        // three points representing two lines l0=(p0,p1), l1 = (p0,p2)
        auto const p0 = vec2dp{0, 0, 1};
        auto const p1 = vec2dp{1, 0, 1};
        auto const p2 = vec2dp{0, 1, 1};

        auto vp = std::vector<vec2dp>{p0, p1, p2};

        double phi = deg2rad(15);  // turning angle (turning plane is e12)
        auto t = vec2dp{-2, 1, 0}; // translation vector

        // stepwise movement: turn around pc by 15°, then shift by t
        auto m_rot = get_motor(pc, phi);
        auto m_tra = get_motor(t);

        auto m_rt_c = rgpr(m_tra, m_rot); // concatenation: translation after rotation
        auto m_tr_c = rgpr(m_rot, m_tra); // concatenation: rotation after translation
                                          // sequence matters!

        auto vr = move2dp(vp, m_rot);
        auto vt = move2dp(vp, m_tra);
        auto vrt = move2dp(vr, m_tra);
        auto vtr = move2dp(vt, m_rot);
        auto vrt_c = move2dp(vp, m_rt_c);
        auto vtr_c = move2dp(vp, m_tr_c);

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

        CHECK(vrt == vrt_c); // same result stepwise or combined
        CHECK(vtr == vtr_c); // same result stepwise or combined

        fmt::println("");
    }

    TEST_CASE("Vec2dp: modeling force & torque = forque")
    {

        fmt::println("Vec2dp: modeling force & torque = forque");

        // points and directions are distinguished (different from ega)
        auto O = O_2dp;            // O is the origin (= an arbitrary reference point)
        auto R = vec2dp{1, 3, 1};  // R is the point of action of force f
        auto r = R - O;            // r is the direction from origin towards R
        auto C = vec2dp{1, 1, 1};  // C in an arbitrary point (e.g. of a rigid body B)
        auto c = C - O;            // c is the direction from origin towards C
        auto f = vec2dp{4, -1, 0}; // force f acting on B at point R

        auto lo = wdg(O, f) + wdg(r, f);     // line of action of f w.r.t. origin O
        auto lc = wdg(C, f) + wdg(r - c, f); // line of action of f w.r.t. C

        // fmt::println("");
        // fmt::println("R  = {}", R);
        // fmt::println("r  = {}", r);
        // fmt::println("C  = {}", C);
        // fmt::println("c  = {}", c);
        // fmt::println("f  = {}", f);
        // fmt::println("");
        // fmt::println("wdg(O,f)                       "
        //              "- effect of f on O = {}",
        //              wdg(O, f));
        // fmt::println("wdg(r,f)        "
        //              "- effect of f on lever r w.r.t. O = {}",
        //              wdg(r, f));
        // fmt::println("lo = wdg(O+r,f) = wdg(R,f)      "
        //              "- combined effect = {}",
        //              wdg(R, f));
        // fmt::println("");
        // fmt::println("wdg(C,f) = wdg(O+c,f)          "
        //              "- effect of f on C = {}",
        //              wdg(C, f));
        // fmt::println(" wdg(O,f)                      "
        //              "- effect of f on O = {}",
        //              wdg(O, f));
        // fmt::println(" wdg(c,f)       "
        //              "- effect of f on lever c w.r.t. O = {}",
        //              wdg(c, f));
        // fmt::println("wdg(r-c,f)  "
        //              "- effect of f on lever (r-c) w.r.t. C = {}",
        //              wdg(r - c, f));
        // fmt::println("lc = wdg(C,f) + wdg(r-c,f) = wdg(r,f) "
        //              "- combined effect = {}",
        //              wdg(O + c + r - c, f));
        // fmt::println("");
        // fmt::println("bulk(lo)                 "
        //              "- models torque w.r.t. O = {}",
        //              bulk(lo));
        // fmt::println("att(lo)         "
        //              "- models force felt a every point = {}",
        //              att(lo));
        // fmt::println("");
        // fmt::println("bulk(lc)=                "
        //              "- models torque w.r.t. O = {}",
        //              bulk(lc));
        // fmt::println("lc - wdg(C,f)            "
        //              "- models torque w.r.t. C = {}",
        //              lc - wdg(C, f));
        // fmt::println("wdg(R-C,f)               "
        //              "- models torque w.r.t. C = {}",
        //              wdg(R - C, f));
        // fmt::println("wdg(r-c,f)               "
        //              "- models torque w.r.t. C = {}",
        //              wdg(r - c, f));
        // fmt::println("att(lc)         "
        //              "- models force felt a every point = {}",
        //              att(lc));
        // fmt::println("");

        CHECK(wdg(O, f) + wdg(r, f) == lo);
        CHECK(wdg(C, f) + wdg(r - c, f) == lc);
        CHECK(wdg(R, f) == wdg(O, f) + wdg(r, f));
        CHECK(wdg(C, f) == wdg(O, f) + wdg(c, f));
        CHECK(lo == lc);
        CHECK(bulk(lo) == bulk(lc));
        CHECK(att(lo) == f);
        CHECK(lc - wdg(C, f) == wdg(R - C, f));
        CHECK(wdg(R - C, f) == wdg(r - c, f));
    }


    ////////////////////////////////////////////////////////////////////////////////
    // MVec2dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2dp: default init")
    {
        fmt::println("MVec2dp: default init");
        // default initialization
        mvec2dp v;
        // fmt::println("   v = {}", v);
        CHECK(abs(v.c0) < eps);
        CHECK(abs(v.c1) < eps);
        CHECK(abs(v.c2) < eps);
        CHECK(abs(v.c3) < eps);
        CHECK(abs(v.c4) < eps);
        CHECK(abs(v.c5) < eps);
        CHECK(abs(v.c6) < eps);
        CHECK(abs(v.c7) < eps);
    }
    TEST_CASE("MVec2dp: with curly braced intializer")
    {
        fmt::println("MVec2dp: with curly braced intializer");
        // default initialization
        mvec2dp v{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        // fmt::println("   v = {}", v);
        CHECK(abs(v.c0 - 0.0) < eps);
        CHECK(abs(v.c1 - 1.0) < eps);
        CHECK(abs(v.c2 - 2.0) < eps);
        CHECK(abs(v.c3 - 3.0) < eps);
        CHECK(abs(v.c4 - 23.0) < eps);
        CHECK(abs(v.c5 - 31.0) < eps);
        CHECK(abs(v.c6 - 12.0) < eps);
        CHECK(abs(v.c7 - 123.0) < eps);
    }

    TEST_CASE("MVec2dp: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec2dp: cp ctor & cp assign incl. type deduction");
        // default initialization
        mvec2dp v1{0.0,  1.0,  2.0,  3.0,
                   23.0, 31.0, 12.0, 123.0}; // init with double (type deduction)
        mvec2dp v2{v1};                      // cp ctor
        mvec2dp v3 = v2;                     // cp assign
        mvec2dp v4 = -v3;                    // cp assign with unary minus

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

        CHECK(abs(v3.c0 - 0.0) < eps);
        CHECK(abs(v3.c1 - 1.0) < eps);
        CHECK(abs(v3.c2 - 2.0) < eps);
        CHECK(abs(v3.c3 - 3.0) < eps);
        CHECK(abs(v3.c4 - 23.0) < eps);
        CHECK(abs(v3.c5 - 31.0) < eps);
        CHECK(abs(v3.c6 - 12.0) < eps);
        CHECK(abs(v3.c7 - 123.0) < eps);

        CHECK(v4 == -v3);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }

    TEST_CASE("MVec2dp: defining basic types and ctor checks")
    {
        fmt::println("");
        fmt::println("MVec2dp: defining basic types and ctor checks");

        auto s = scalar2dp{5.0};
        auto v = vec2dp{1.0, 2.0, 1.0};
        auto B = bivec2dp{-1.0, 2.0, 1.0};
        auto ps = pscalar2dp{-5.0};

        auto mv0 = mvec2dp{s};
        auto mv1 = mvec2dp{v};
        auto mv2 = mvec2dp{B};
        auto mv3 = mvec2dp{ps};

        auto mv5a = mvec2dp{s, B};
        auto mv5b = mvec2dp_e{s};
        auto mv5c = mvec2dp_e{B};
        auto mv5d = mvec2dp_e{s, B};
        auto mv5e = mvec2dp{mv5d};

        auto mv6a = mvec2dp{v, ps};
        auto mv6b = mvec2dp_u{v};
        auto mv6c = mvec2dp_u{ps};
        auto mv6d = mvec2dp_u{v, ps};
        auto mv6e = mvec2dp{mv6d};

        auto mv7 = mvec2dp{s, v, B, ps};

        // fmt::println("   mv0  = {}", mv0);
        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv5a = {}", mv5a);
        // fmt::println("   mv5b = {}", mv5b);
        // fmt::println("   mv5c = {}", mv5c);
        // fmt::println("   mv5d = {}", mv5d);
        // fmt::println("   mv5e = {}", mv5e);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6b = {}", mv6b);
        // fmt::println("   mv6c = {}", mv6c);
        // fmt::println("   mv6d = {}", mv6d);
        // fmt::println("   mv6e = {}", mv6e);
        // fmt::println("   mv7  = {}", mv7);

        CHECK(gr0(mv0) == s);
        CHECK(gr1(mv1) == v);
        CHECK(gr2(mv2) == B);
        CHECK(gr3(mv3) == ps);

        CHECK(gr0(mv5a) == s);
        CHECK(gr1(mv5a) == vec2dp{});
        CHECK(gr2(mv5a) == B);
        CHECK(gr3(mv5a) == pscalar2dp{});

        CHECK(gr0(mv5b) == s);
        CHECK(gr2(mv5b) == bivec2dp{});

        CHECK(gr0(mv5c) == scalar2dp{});
        CHECK(gr2(mv5c) == B);

        CHECK(gr0(mv5d) == s);
        CHECK(gr2(mv5d) == B);

        CHECK(gr0(mv5e) == s);
        CHECK(gr1(mv5e) == vec2dp{});
        CHECK(gr2(mv5e) == B);
        CHECK(gr3(mv5e) == pscalar2dp{});

        CHECK(gr0(mv6a) == scalar2dp{});
        CHECK(gr1(mv6a) == v);
        CHECK(gr2(mv6a) == bivec2dp{});
        CHECK(gr3(mv6a) == ps);

        CHECK(gr1(mv6b) == v);
        CHECK(gr3(mv6b) == pscalar2dp{});

        CHECK(gr1(mv6c) == vec2dp{});
        CHECK(gr3(mv6c) == ps);

        CHECK(gr1(mv6d) == v);
        CHECK(gr3(mv6d) == ps);

        CHECK(gr0(mv6e) == scalar2dp{});
        CHECK(gr1(mv6e) == v);
        CHECK(gr2(mv6e) == bivec2dp{});
        CHECK(gr3(mv6e) == ps);

        CHECK(gr0(mv7) == s);
        CHECK(gr1(mv7) == v);
        CHECK(gr2(mv7) == B);
        CHECK(gr3(mv7) == ps);

        CHECK(gr(s) == 0);
        CHECK(gr(v) == 1);
        CHECK(gr(B) == 2);
        CHECK(gr(ps) == 3);
    }

    TEST_CASE("MVec2dp: fmt & cout printing")
    {
        fmt::println("MVec2dp: fmt & cout printing");

        mvec2dp pf{1.0f, 2.00001f, 0.0f, 3.0f, 1.0f, 2.00001f, 0.0f, 3.0f};
        mvec2dp pd{1.0, 2.00001, 0.0, 3.0, 1.0, 2.00001, 0.0, 3.0};

        fmt::println("");
        std::cout << "    cout: pf = " << pf << std::endl;
        fmt::println("    fmt:  pf = {}", pf);
        fmt::println("    fmt:  pf = {:.8f}", pf);

        std::cout << "    cout: pd = " << pd << std::endl;
        fmt::println("    fmt:  pd = {}", pd);
        fmt::println("    fmt:  pd = {:.8f}", pd);

        std::vector<MVec2dp<double>> vp1{{1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0},
                                         {0.5, 1.5, 2.0, 2.5, 1.0, 1.0, 1.0, 2.0}};
        fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("    fmt: vp1 = {:e}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("MVec2dp: vector space and linearity tests")
    {
        fmt::println("MVec2dp: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        mvec2dp p0;
        mvec2dp p1{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
        mvec2dp p2{0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0};
        mvec2dp p3{0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0};
        mvec2dp p4 = -p1; // assignment using unary minus
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

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2dp<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2dp: wedge product - basic properties")
    {
        fmt::println("MVec2dp: wedge product - basic properties");

        vec2dp v1{1.0, 2.0, 1.0};
        vec2dp v2{0.5, 3.0, 2.0};
        vec2dp v3{-2.0, 6.0, 3.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar2dp{sd};
        auto t = scalar2dp{st};

        CHECK(wdg(v1, v1) == bivec2dp{});                     // wdg=0 for collin. vectors
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

    TEST_CASE("MVec2dp: geometric product - basic properties")
    {
        fmt::println("MVec2dp: geometric product - basic properties");

        vec2dp v1{1.0, 2.0, 1.0};
        vec2dp v2{0.5, 3.0, 2.0};
        vec2dp v3{-2.0, 6.0, 3.0};
        bivec2dp B1{-4, 2, 1};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar2dp{sd};
        auto t = scalar2dp{st};

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

    TEST_CASE("MVec2dp: geometric product - combinatorial tests")
    {
        fmt::println("MVec2dp: geometric product - combinatorial tests");

        scalar2dp s1{3.0};
        scalar2dp s2{-1.5};
        vec2dp v1{1.0, 2.0, 3.0};
        vec2dp v2{0.5, 3.0, -2.0};
        bivec2dp b1{-1.0, -3.0, -5.0};
        bivec2dp b2{2.0, 0.5, -4.0};
        pscalar2dp ps1{-5.0};
        pscalar2dp ps2{2.0};

        // checks all ctor combinations and equivalence of simplified products with
        // geometric product of a fully populated multivector
        //
        // assumption: the geometric product of a fully populated multivector serves as
        //             reference and is correctly implemented
        //             results are relative to full geometric product

        // MVec2dp_E * MVec2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{s1, v1, b1, ps1} ==
              mvec2dp{mvec2dp_e{s1, b1} * mvec2dp{s1, v1, b1, ps1}});

        // MVec2dp * MVec2dp_E
        CHECK(mvec2dp{s1, v1, b1, ps1} * mvec2dp{s1, b1} ==
              mvec2dp{mvec2dp{s1, v1, b1, ps1} * mvec2dp_e{s1, b1}});

        // MVec2dp_E * MVec2dp_E
        CHECK(mvec2dp{s1, b1} * mvec2dp{s2, b2} ==
              mvec2dp{mvec2dp_e{s1, b1} * mvec2dp_e{s2, b2}});

        // MVec2dp_U * MVec2dp_U
        CHECK(mvec2dp{v1, ps1} * mvec2dp{v2, ps2} ==
              mvec2dp{mvec2dp_u{v1, ps1} * mvec2dp_u{v2, ps2}});

        // MVec2dp_E * MVec2dp_U
        CHECK(mvec2dp{s1, b1} * mvec2dp{v2, ps2} ==
              mvec2dp{mvec2dp_e{s1, b1} * mvec2dp_u{v2, ps2}});

        // MVec2dp_U * MVec2dp_E
        CHECK(mvec2dp{v1, ps1} * mvec2dp{s2, b2} ==
              mvec2dp{mvec2dp_u{v1, ps1} * mvec2dp_e{s2, b2}});

        // MVec2dp_E * BiVec2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{b2} == mvec2dp{mvec2dp_e{s1, b1} * mvec2dp{b2}});

        // BiVec2dp * MVec2dp_E
        CHECK(mvec2dp{b1} * mvec2dp{s2, b2} == mvec2dp{mvec2dp{b1} * mvec2dp_e{s2, b2}});

        // MVec2dp_E * Vec2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{v2} == mvec2dp{mvec2dp_e{s1, b1} * mvec2dp{v2}});

        // Vec2dp * MVec2dp_E
        CHECK(mvec2dp{v1} * mvec2dp{s2, b2} == mvec2dp{mvec2dp{v1} * mvec2dp_e{s2, b2}});

        // MVec2dp_U * BiVec2dp
        CHECK(mvec2dp{v1, ps1} * mvec2dp{b2} ==
              mvec2dp{mvec2dp_u{v1, ps1} * mvec2dp{b2}});

        // BiVec2dp * MVec2dp_U
        CHECK(mvec2dp{b1} * mvec2dp{v2, ps2} ==
              mvec2dp{mvec2dp{b1} * mvec2dp_u{v2, ps2}});

        // BiVec2dp * BiVec2dp
        CHECK(mvec2dp{b1} * mvec2dp{b2} == mvec2dp{b1 * b2});

        // BiVec2dp * Vec2dp
        CHECK(mvec2dp{b1} * mvec2dp{v1} == mvec2dp{b1 * v1});

        // Vec2dp * BiVec2dp
        CHECK(mvec2dp{v1} * mvec2dp{b1} == mvec2dp{v1 * b1});

        // Vec2dp * Vec2dp
        CHECK(mvec2dp{v1} * mvec2dp{v2} == mvec2dp{v1 * v2});


        // PScalar2dp * MVec2dp
        CHECK(mvec2dp{ps1} * mvec2dp{s2, v2, b2, ps2} ==
              mvec2dp{ps1 * mvec2dp{s2, v2, b2, ps2}});

        // MVec2dp * Pscalar2dp
        CHECK(mvec2dp{s1, v1, b1, ps1} * mvec2dp{ps2} ==
              mvec2dp{mvec2dp{s1, v1, b1, ps1} * ps2});

        // PScalar2dp * MVec2dp_E
        CHECK(mvec2dp{ps1} * mvec2dp{s2, b2} == mvec2dp{ps1 * mvec2dp_e{s2, b2}});

        // MVec2dp_E * Pscalar2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{ps2} == mvec2dp{mvec2dp_e{s1, b1} * ps2});

        // PScalar2dp * MVec2dp_U
        CHECK(mvec2dp{ps1} * mvec2dp{v2, ps2} == mvec2dp{ps1 * mvec2dp_u{v2, ps2}});

        // MVec2dp_U * Pscalar2dp
        CHECK(mvec2dp{v1, ps1} * mvec2dp{ps2} == mvec2dp{mvec2dp_u{v1, ps1} * ps2});

        // PScalar2dp * BiVec2dp
        CHECK(mvec2dp{ps1} * mvec2dp{b2} == mvec2dp{ps1 * b2});

        // BiVec2dp * Pscalar2dp
        CHECK(mvec2dp{b1} * mvec2dp{ps2} == mvec2dp{b1 * ps2});

        // PScalar2dp * Vec2dp
        CHECK(mvec2dp{ps1} * mvec2dp{v2} == mvec2dp{ps1 * v2});

        // Vec2dp * Pscalar2dp
        CHECK(mvec2dp{v1} * mvec2dp{ps2} == mvec2dp{v1 * ps2});

        // Pscalar2dp * Pscalar2dp
        CHECK(mvec2dp{ps1} * mvec2dp{ps2} == mvec2dp{ps1 * ps2});

        // Scalar * Scalar
        CHECK(mvec2dp{s1} * mvec2dp{s2} == mvec2dp{s1 * s2});
    }

    TEST_CASE("MVec2dp: geometric product tests - vec * vec")
    {
        fmt::println("MVec2dp: geometric product tests - vec * vec");

        // ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
        //
        // dot(a,b) = 0.5*(ab + ba)   (symmetric part)
        // wdg(a,b) = 0.5*(ab - ba)   (antisymmetric part)

        vec2dp a{1.0, 2.0, 3.0};
        vec2dp b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec2dp mva{a};
        mvec2dp mvb{b};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb + mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb - mvb * mva);

        // fmt::println("   a = {}", a);
        // fmt::println("   b = {}", b);
        // fmt::println("   dot(a,b) = {}", dot_ab);
        // fmt::println("   wdg(a,b) = {}", wdg_ab);
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

        CHECK(dot_ab == gr0(mvab));
        CHECK(dot_ab == gr0(mvab_sym));
        CHECK(wdg_ab == gr2(mvab));
        CHECK(wdg_ab == gr2(mvab_asym));
    }

    TEST_CASE("MVec2dp: geometric product tests - bivec * vec")
    {
        fmt::println("MVec2dp: geometric product tests - bivec * vec");

        // Ab = (b << A) + wdg(A,b) = gr1(Ab) + gr3(Ab)
        //
        // (b << A) = 0.5*(Ab - Aa)   (antisymmetric part)
        // wdg(A,b) = 0.5*(Ab + Aa)   (symmetric part)

        bivec2dp A{1.0, 2.0, 3.0};
        vec2dp b{0.5, 3.0, -2.0};
        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        mvec2dp mva{A};
        mvec2dp mvb{b};
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

    TEST_CASE("MVec2dp: geometric product tests - vec * bivec")
    {
        fmt::println("MVec2dp: geometric product tests - vec * bivec");

        // a*B = (B >> a) + wdg(a,B) = gr1(aB) + gr3(aB)
        //
        // (B >> a) = 0.5*(aB - Ba)   (antisymmetric part)
        // wdg(a,B) = 0.5*(aB + Ba)   (symmetric part)

        vec2dp a{1.0, 2.0, 3.0};
        bivec2dp B{0.5, 3.0, -2.0};
        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec2dp mva{a};
        mvec2dp mvb{B};
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

    TEST_CASE("MVec2dp: geometric product - link to inner and outer products")
    {
        fmt::println("MVec2dp: geometric product - link to inner and outer products");
        // vec2dp u{1.0, 2.0, 3.0};  // 3d vector from ega3d case
        vec2dp u{1.0, 2.0, 1.0}; //     2d point in pga2dp

        scalar2dp s{3};
        // vec2dp v{-3.0, 2.5, -0.5}; // 3d vector from ega3d case
        vec2dp v{-3.0, 2.5, 1.0}; //     2d point in pga2dp
        bivec2dp B{2.5, 3.5, 1.5};
        pscalar2dp ps{4.0};

        // fmt::println("");
        // fmt::println("u = {}", u);
        // fmt::println("s = {}", s);
        // fmt::println("v = {}", v);
        // fmt::println("B = {}", B);
        // fmt::println("ps = {}", ps);
        // fmt::println("");

        // fmt::println("");
        // fmt::println("scalar case: ");
        CHECK(u * s == gr1((s >> u) + wdg(u, s)));
        CHECK(u * s == gr1(rwdg(s, bulk_dual(u)) + wdg(u, s)));

        // fmt::println("");
        // fmt::println("u * s = {}", u * s);
        // fmt::println("");
        // fmt::println("(s >> u) = {}", (s >> u));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("(s >> u) + wdg(u, s) = {}", (s >> u) + wdg(u, s));
        // fmt::println("gr1((s >> u) + wdg(u, s)) = {}", gr1((s >> u) + wdg(u, s)));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(s,bulk_dual(u))= {}", rwdg(s, bulk_dual(u)));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("rwdg(s,bulk_dual(u)) + wdg(u, s) = {}",
        //              rwdg(s, bulk_dual(u)) + wdg(u, s));
        // fmt::println("gr1(rwdg(s,bulk_dual(u)) + wdg(u, s)) = {}",
        //              gr1(rwdg(s, bulk_dual(u)) + wdg(u, s)));
        // fmt::println("");

        CHECK(s * u == gr1((u << s) + wdg(s, u)));
        CHECK(s * u == gr1(rwdg(bulk_dual(u), s) + wdg(s, u)));

        // fmt::println("");
        // fmt::println("s * u = {}", s * u);
        // fmt::println("");
        // fmt::println("(u << s) = {}", (u << s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("(u << s) + wdg(s, u) = {}", (u << s) + wdg(s, u));
        // fmt::println("gr1((u << s) + wdg(s, u)) = {}", gr1((u << s) + wdg(s, u)));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(bulk_dual(u), s)= {}", rwdg(bulk_dual(u), s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("rwdg(lbulk_dual(u), s) + wdg(s, u) = {}",
        //              rwdg(bulk_dual(u), s) + wdg(s, u));
        // fmt::println("gr1(rwdg(bulk_dual(u), s) + wdg(s, u)) = {}",
        //              gr1(rwdg(bulk_dual(u), s) + wdg(s, u)));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("vector case: ");
        CHECK(u * v == (v >> u) + wdg(u, v));
        CHECK(u * v == rwdg(v, bulk_dual(u)) + wdg(u, v));

        // fmt::println("");
        // fmt::println("u * v = {}", u * v);
        // fmt::println("");
        // fmt::println("(v >> u) = {}", (v >> u));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("(v >> u) + wdg(u, v) = {}", (v >> u) + wdg(u, v));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(v,bulk_dual(u))= {}", rwdg(v, bulk_dual(u)));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("rwdg(v,bulk_dual(u)) + wdg(u, v) = {}",
        //              rwdg(v, bulk_dual(u)) + wdg(u, v));
        // fmt::println("");

        CHECK(v * u == (u << v) + wdg(v, u));
        CHECK(v * u == rwdg(bulk_dual(u), v) + wdg(v, u));

        // fmt::println("");
        // fmt::println("v * u = {}", v * u);
        // fmt::println("");
        // fmt::println("(u << v) = {}", (u << v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("(u << v) + wdg(v, u) = {}", (u << v) + wdg(v, u));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(bulk_dual(u), v)= {}", rwdg(bulk_dual(u), v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("rwdg(bulk_dual(u), v) + wdg(v, u) = {}",
        //              rwdg(bulk_dual(u), v) + wdg(v, u));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("bivector case: ");
        CHECK(u * B == (B >> u) + wdg(u, B));
        CHECK(u * B == rwdg(B, bulk_dual(u)) + wdg(u, B));

        // fmt::println("");
        // fmt::println("u * B = {}", u * B);
        // fmt::println("");
        // fmt::println("(B >> u) = {}", (B >> u));
        // fmt::println("wdg(u, B) = {}", wdg(u, B));
        // fmt::println("(B >> u) + wdg(u, B) = {}", (B >> u) + wdg(u, B));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(B,bulk_dual(u))= {}", rwdg(B, bulk_dual(u)));
        // fmt::println("wdg(u, B) = {}", wdg(u, B));
        // fmt::println("rwdg(B,bulk_dual(u)) + wdg(u, B) = {}",
        //              rwdg(B, bulk_dual(u)) + wdg(u, B));
        // fmt::println("");

        CHECK(B * u == (u << B) + wdg(B, u));
        CHECK(B * u == rwdg(bulk_dual(u), B) + wdg(B, u));

        // fmt::println("");
        // fmt::println("B * u = {}", B * u);
        // fmt::println("");
        // fmt::println("(u << B) = {}", (u << B));
        // fmt::println("wdg(B, u) = {}", wdg(B, u));
        // fmt::println("(u << B) + wdg(B, u) = {}", (u << B) + wdg(B, u));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(bulk_dual(u), B)= {}", rwdg(bulk_dual(u), B));
        // fmt::println("wdg(B, u) = {}", wdg(B, u));
        // fmt::println("rwdg(bulk_dual(u), B) + wdg(B, u) = {}",
        //              rwdg(bulk_dual(u), B) + wdg(B, u));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("pscalar case: ");
        CHECK(u * ps == gr2((ps >> u) + wdg(u, ps)));
        CHECK(u * ps == gr2(rwdg(ps, bulk_dual(u)) + wdg(u, ps)));

        // fmt::println("");
        // fmt::println("u * ps = {}", u * ps);
        // fmt::println("");
        // fmt::println("(ps >> u) = {}", (ps >> u));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("(ps >> u) + wdg(u, ps) = {}", (ps >> u) + wdg(u, ps));
        // fmt::println("gr2((ps >> u) + wdg(u, ps)) = {}", gr2((ps >> u) + wdg(u, ps)));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(ps,bulk_dual(u))= {}", rwdg(ps, bulk_dual(u)));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("rwdg(ps,bulk_dual(u)) + wdg(u, ps) = {}",
        //              rwdg(ps, bulk_dual(u)) + wdg(u, ps));
        // fmt::println("gr2(rwdg(ps,bulk_dual(u)) + wdg(u, ps)) = {}",
        //              gr2(rwdg(ps, bulk_dual(u)) + wdg(u, ps)));
        // fmt::println("");

        CHECK(ps * u == gr2((u << ps) + wdg(ps, u)));
        CHECK(ps * u == gr2(rwdg(bulk_dual(u), ps) + wdg(ps, u)));

        // fmt::println("");
        // fmt::println("ps * u = {}", ps * u);
        // fmt::println("");
        // fmt::println("(u << ps) = {}", (u << ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("(u << ps) + wdg(ps, u) = {}", (u << ps) + wdg(ps, u));
        // fmt::println("gr2((u << ps) + wdg(ps, u)) = {}", gr2((u << ps) + wdg(ps, u)));
        // fmt::println("");
        // fmt::println("bulk_dual(u) = {}", bulk_dual(u));
        // fmt::println("rwdg(bulk_dual(u), ps)= {}", rwdg(bulk_dual(u), ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("rwdg(bulk_dual(u), ps) + wdg(ps, u) = {}",
        //              rwdg(bulk_dual(u), ps) + wdg(ps, u));
        // fmt::println("gr2(rwdg(bulk_dual(u), ps) + wdg(ps, u)) = {}",
        //              gr2(rwdg(bulk_dual(u), ps) + wdg(ps, u)));
        // fmt::println("");

        // fmt::println("");
    }

    TEST_CASE("MVec2dp: geometric product tests - equivalence tests")
    {
        fmt::println("MVec2dp: geometric product tests - equivalence tests");

        vec2dp a{1.0, 2.0, 3.0};
        vec2dp b{0.5, 3.0, -4.0};
        mvec2dp mva{a};
        mvec2dp mvb{b};

        bivec2dp A{1.0, 2.0, 3.0};
        bivec2dp B{0.5, 3.0, -4.0};
        mvec2dp mvA{A};
        mvec2dp mvB{B};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec2dp_e ab = a * b;
        mvec2dp abm = mva * mvb;
        mvec2dp abd{scalar2dp{dot_ab}, wdg_ab};

        mvec2dp_u Ab = A * b;
        mvec2dp Abm = mvA * mvb;
        mvec2dp Abd{b_onto_A, wdg_Ab};

        mvec2dp_u aB = a * B;
        mvec2dp aBm = mva * mvB;
        mvec2dp aBd{B_by_a, wdg_aB};

        // fmt::println("");
        // fmt::println("   a                                 = {}", a);
        // fmt::println("   mva                               = {}", mva);
        // fmt::println("   b                                 = {}", b);
        // fmt::println("   mvb                               = {}", mvb);
        // fmt::println("   ab  = mvec2dp_e(a * b)            = {}", ab);
        // fmt::println("   abm = mva * mvb                   = {}", abm);
        // fmt::println("   abd = mvec2dp(dot(a,b), wdg(a,b)) = {}", abd);
        // fmt::println("");
        // fmt::println("   A                                 = {}", A);
        // fmt::println("   mvA                               = {}", mvA);
        // fmt::println("   b                                 = {}", b);
        // fmt::println("   mvb                               = {}", mvb);
        // fmt::println("   Ab  = mvec2dp_u(A * b)            = {}", Ab);
        // fmt::println("   Abm = mvA * mvb                   = {}", Abm);
        // fmt::println("   Abd = mvec2dp((b << A), wdg(A,b)) = {}", Abd);
        // fmt::println("");
        // fmt::println("   a                                 = {}", a);
        // fmt::println("   mva                               = {}", mva);
        // fmt::println("   B                                 = {}", B);
        // fmt::println("   mvB                               = {}", mvB);
        // fmt::println("   aB  = mvec2dp_u(a * B)            = {}", aB);
        // fmt::println("   aBm = mva * mvB                   = {}", aBm);
        // fmt::println("   aBd = mvec2dp((B >> a), wdg(a,B)) = {}", aBd);
        // fmt::println("");

        CHECK(gr0(ab) == gr0(abm));
        CHECK(gr1(abm) == vec2dp{});
        CHECK(gr2(ab) == gr2(abm));
        CHECK(gr3(abm) == pscalar2dp{0.0});

        CHECK(gr0(ab) == gr0(abd));
        CHECK(gr1(abd) == vec2dp{});
        CHECK(gr2(ab) == gr2(abd));
        CHECK(gr3(abd) == pscalar2dp{0.0});

        CHECK(gr0(Abm) == 0);
        CHECK(gr1(Ab) == gr1(Abm));
        CHECK(gr2(Abm) == bivec2dp{});
        CHECK(gr3(Ab) == gr3(Abm));

        CHECK(gr0(Abd) == 0);
        CHECK(gr1(Ab) == gr1(Abd));
        CHECK(gr2(Abd) == bivec2dp{});
        CHECK(gr3(Ab) == gr3(Abd));

        CHECK(gr0(aBm) == 0);
        CHECK(gr1(aB) == gr1(aBm));
        CHECK(gr2(aBm) == bivec2dp{});
        CHECK(gr3(aB) == gr3(aBm));

        CHECK(gr0(aBd) == 0);
        CHECK(gr1(aB) == gr1(aBd));
        CHECK(gr2(aBd) == bivec2dp{});
        CHECK(gr3(aB) == gr3(aBd));
    }

    TEST_CASE("MVec2dp: assignment tests")
    {
        fmt::println("MVec2dp: assignment tests");

        vec2dp v1{1.0, 2.0, 3.0};
        vec2dp v2{0.5, 1.0, 1.5};
        vec2dp v3{0.5, 1.0, -4.5};
        bivec2dp b1{1.0, 2.0, 3.0};

        mvec2dp mv1{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        mvec2dp mv2{0.0, 0.5, 1.0, 1.5, 11.5, 15.5, 6.0, 61.5};
        mvec2dp mv3{mv1};
        mvec2dp mv4 = mv2;

        mvec2dp mv5(scalar2dp(5.0));
        mvec2dp mv6(pscalar2dp(6.0));
        mvec2dp mv7{v1};
        mvec2dp mv8{b1};
        mvec2dp mv9{scalar2dp{dot(v1, v3)}, wdg(v1, v3)};

        mvec2dp mv10{v1, pscalar2dp(10.0)};

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
        // fmt::println("   gr1(mv3) = {}", gr1(mv4));
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
        CHECK(gr3(mv6) == 6.0);
        CHECK(mv1 == mv3);
        CHECK(mv4 == mv2);
        CHECK(gr1(mv7) == v1);
        CHECK(gr2(mv8) == b1);
        CHECK(gr0(mv9) == dot(v1, v3));
        CHECK(gr2(mv9) == wdg(v1, v3));
        CHECK(gr1(mv10) == v1);
        CHECK(gr3(mv10) == 10.0);
    }

    TEST_CASE("MVec2dp: bivector product properties")
    {
        fmt::println("MVec2dp: bivector product properties");

        bivec2dp b1{1.0, 2.0, 3.0};
        mvec2dp mb1{b1};
        bivec2dp b2{-3.0, 1.0, 2.0};
        mvec2dp mb2{b2};

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
        // fmt::println("   gpr12_m = mb1 * mb2 = {}", gpr12_m);
        // fmt::println("");
        // fmt::println("   dot(b2, b1)         = {}", dot(b2, b1));
        // fmt::println("   cmt(b2, b1)         = {}", cmt(b2, b1));
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
    // MVec2dp_E<T> and MVec2dp_U<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2dp_E/_U: reflections and motors (rotations, translations)")
    {
        fmt::println("MVec2dp_E/_U: reflections and motors (rotations, translations)");

        {
            /////////////////////////////////////////////////////////////////////////////
            // 1st: start with points use lines through the origin and reflect points
            /////////////////////////////////////////////////////////////////////////////

            // define points and lines
            auto p0 = O_2dp;
            auto p1 = vec2dp{1, 0, 1};
            auto p2 = vec2dp{1, 1, 1};
            auto p = vec2dp{1, -0.5, 1};
            // auto p = vec2dp{2, -1, 2};

            auto l1 = wdg(p0, p1);
            auto l2 = wdg(p0, p2);

            // reflect p on l1 -> pr and reflect pr on l2 -> prr
            auto pr = -gr1(rgpr(rgpr(l1, p), l1));

            auto prr = -gr1(rgpr(rgpr(l2, pr), l2));

            fmt::println("test_refl");
            fmt::println("");
            fmt::println("l1: {}, l1u: {}, l2: {}, l2u: {}", l1, unitize(l1), l2,
                         unitize(l2));
            fmt::println("pr: {}, pru: {}", pr, unitize(pr));
            fmt::println("prr: {}, prru: {}", prr, unitize(prr));
            fmt::println("");

            CHECK(unitize(pr) == vec2dp{1, 0.5, 1});
            CHECK(unitize(prr) == vec2dp{0.5, 1, 1});

            // show that prr can be obtained directly from a rotation via a motor as pm
            auto motor = rgpr(l2, l1);
            auto rmotor = rrev(motor);

            auto pm = gr1(rgpr(rgpr(motor, p), rmotor));

            fmt::println("");
            fmt::println("motor: {}, rmotor: {}", motor, rmotor);
            fmt::println("");
            fmt::println("pm: {}, pmu: {}", pm, unitize(pm));
            fmt::println("");
            auto motor_calc = get_motor(O_2dp, deg2rad(90));
            fmt::println("motor_calc: {}, rmotor_calc: {}", motor_calc, rrev(motor_calc));
            fmt::println("p_calc = {}", move2dp(p, motor_calc));


            CHECK(unitize(prr) == unitize(pm));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 2nd: start with line through the origin and reflect lines through the
            //      origin with them
            /////////////////////////////////////////////////////////////////////////////

            // define points and lines
            auto p0 = O_2dp;
            auto p1 = vec2dp{1, 0, 1};
            auto p2 = vec2dp{1, 1, 1};
            auto p = vec2dp{1, -0.5, 1};

            auto l = wdg(p0, p);
            auto l1 = unitize(wdg(p0, p1));
            auto l2 = unitize(wdg(p0, p2));

            // reflect l on l1 -> lr and reflect lr on l2 -> lrr
            auto lr = -gr2(rgpr(rgpr(l1, l), l1));

            auto lrr = -gr2(rgpr(rgpr(l2, lr), l2));

            // fmt::println("");
            // fmt::println("l: {}, lu: {}, att(l): {}", l, unitize(l), att(l));
            // fmt::println("");
            // fmt::println("l1: {}, l1u: {}, att(l1): {}", l1, unitize(l1), att(l1));
            // fmt::println("l2: {}, l2u: {}, att(l2): {}", l2, unitize(l2), att(l2));
            // fmt::println("");
            // fmt::println("lr: {}, lru: {}, att(lr): {}", lr, unitize(lr), att(lr));
            // fmt::println("lrr: {}, lrru: {}, att(lrr): {}", lrr, unitize(lrr),
            // att(lrr)); fmt::println("");

            CHECK(lr == -join(p0, vec2dp{1, 0.5, 1})); // reflection turns the direction
            CHECK(lrr == join(p0, vec2dp{0.5, 1, 1}));

            // show that prr can be obtained directly from a rotation via a motor as pm
            auto motor = rgpr(l2, l1);
            auto rmotor = rrev(motor);

            auto lm = gr2(rgpr(rgpr(motor, l), rmotor));

            // fmt::println("");
            // fmt::println("lm: {}, lmu: {}, att(lm): {}", lm, unitize(lm), att(lm));
            // fmt::println("");

            CHECK(lrr == lm);
            CHECK(unitize(lrr) == unitize(lm));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 3rd: as 1st step (reflection with points), but with the lines for
            // reflection and the center of rotation not in the origin
            /////////////////////////////////////////////////////////////////////////////

            // define points and lines
            auto p0 = vec2dp{1, 0.5, 1}; // was O_2dp, now shifted to that new point
            auto delta = p0 - O_2dp;
            auto p1 = vec2dp{1, 0, 1} + delta;
            auto p2 = vec2dp{1, 1, 1} + delta;
            auto p = vec2dp{1, -0.5, 1} + delta;

            auto l1 = wdg(p0, p1);
            auto l2 = wdg(p0, p2);

            // reflect p on l1 -> pr and reflect pr on l2 -> prr
            auto pr = -gr1(rgpr(rgpr(l1, p), l1));

            auto prr = -gr1(rgpr(rgpr(l2, pr), l2));

            // fmt::println("");
            // fmt::println("l1: {}, l1u: {}, l2: {}, l2u: {}", l1, unitize(l1), l2,
            //              unitize(l2));
            // fmt::println("pr: {}, pru: {}", pr, unitize(pr));
            // fmt::println("prr: {}, prru: {}", prr, unitize(prr));
            // fmt::println("");

            CHECK(unitize(pr) == vec2dp{1, 0.5, 1} + delta);
            CHECK(unitize(prr) == vec2dp{0.5, 1, 1} + delta);

            // show that prr can be obtained directly from a movement via a motor as pm
            auto motor = rgpr(l2, l1);
            auto rmotor = rrev(motor);

            auto pm = gr1(rgpr(rgpr(motor, p), rmotor));

            // fmt::println("");
            // fmt::println("pm: {}, pmu: {}", pm, unitize(pm));
            // fmt::println("");

            CHECK(unitize(prr) == unitize(pm));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 4th: same at the 2nd step (reflecion of lines on lines), but with
            // lines for reflection and the center of rotation not in the origin
            /////////////////////////////////////////////////////////////////////////////

            // define points and lines
            auto p0 = vec2dp{1, 0.5, 1}; // was O_2dp, now shifted to that new point
            auto delta = p0 - O_2dp;
            auto p1 = vec2dp{1, 0, 1} + delta;
            auto p2 = vec2dp{1, 1, 1} + delta;
            auto p = vec2dp{1, -0.5, 1} + delta;

            auto l = wdg(p0, p);
            auto l1 = unitize(wdg(p0, p1));
            auto l2 = unitize(wdg(p0, p2));

            // reflect l on l1 -> lr and reflect lr on l2 -> lrr
            auto lr = -gr2(rgpr(rgpr(l1, l), l1));

            auto lrr = -gr2(rgpr(rgpr(l2, lr), l2));

            // fmt::println("");
            // fmt::println("l: {}, lu: {}, att(l): {}", l, unitize(l), att(l));
            // fmt::println("");
            // fmt::println("l1: {}, l1u: {}, att(l1): {}", l1, unitize(l1), att(l1));
            // fmt::println("l2: {}, l2u: {}, att(l2): {}", l2, unitize(l2), att(l2));
            // fmt::println("");
            // fmt::println("lr: {}, lru: {}, att(l2): {}", lr, unitize(lr), att(lr));
            // fmt::println("lrr: {}, lrru: {}, att(lrr): {}", lrr, unitize(lrr),
            // att(lrr)); fmt::println("");

            CHECK(lr ==
                  -join(p0, vec2dp{1, 0.5, 1} + delta)); // reflection turns the direction
            CHECK(lrr == join(p0, vec2dp{0.5, 1, 1} + delta));

            // show that prr can be obtained directly from a rotation via a motor as pm
            auto motor = rgpr(l2, l1);
            auto rmotor = rrev(motor);

            auto lm = gr2(rgpr(rgpr(motor, l), rmotor));

            // fmt::println("");
            // fmt::println("lm: {}, lmu: {}, att(lm): {}", lm, unitize(lm), att(lm));
            // fmt::println("");

            CHECK(lrr == lm);
            CHECK(unitize(lrr) == unitize(lm));
        }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 5th step: do the same as above, but use the optimized motor formulas
            /////////////////////////////////////////////////////////////////////////////

            // define points and lines
            auto p0 = vec2dp{1, 0.5, 1}; // was O_2dp, now shifted to that new point
            auto delta = p0 - O_2dp;
            auto p1 = vec2dp{1, 0, 1} + delta;
            auto p2 = vec2dp{1, 1, 1} + delta;
            auto p = vec2dp{1, -0.5, 1} + delta;

            auto l = wdg(p0, p);
            auto l1 = unitize(wdg(p0, p1));
            auto l2 = unitize(wdg(p0, p2));

            auto R = get_motor_from_lines(l1, l2);
            CHECK(R == rgpr(l2, l1));

            auto pm_manual = gr1(rgpr(rgpr(R, p), rrev(R)));
            auto pm_orig = move2dp(p, R);
            auto pm = move2dp(p, R);

            CHECK(pm_manual == pm_orig);
            CHECK(pm_manual == unitize(pm));

            auto lm_manual = gr2(rgpr(rgpr(R, l), rrev(R)));
            auto lm_orig = move2dp(l, R);
            auto lm = move2dp(l, R);

            // fmt::println("");
            // fmt::println("lm_orig: {}, lmu_orig: {}, att(lm_orig): {}", lm_orig,
            //              unitize(lm_orig), att(lm_orig));
            // fmt::println("");
            // fmt::println("");
            // fmt::println("lm: {}, lmu: {}, att(lm): {}", lm, unitize(lm), att(lm));
            // fmt::println("");

            CHECK(lm_manual == lm_orig);
            CHECK(lm_manual == lm);
        }

        // {
        //     /////////////////////////////////////////////////////////////////////////////
        //     // 6th step: optional speed testing for the optimized motor formulas
        //     // -> advantage for debug mode, but not when enabling optimization
        //     /////////////////////////////////////////////////////////////////////////////

        //     // define points and lines
        //     auto p0 = vec2dp{1, 0.5, 1}; // was O_2dp,
        //     // now shifted to that new point

        //     auto delta = p0 - O_2dp;
        //     auto p1 = vec2dp{1, 0, 1} + delta;
        //     auto p2 = vec2dp{1, 1, 1} + delta;
        //     auto p = vec2dp{1, -0.5, 1} + delta;

        //     auto l = wdg(p0, p);
        //     auto l1 = unitize(wdg(p0, p1));
        //     auto l2 = unitize(wdg(p0, p2));

        //     auto R = get_motor_from_lines(l1, l2);
        //     CHECK(R == rgpr(l2, l1));

        //     // checking time required
        //     constexpr size_t steps = 10'000'000;
        //     vec2dp pm_sum_orig{};
        //     auto start = std::chrono::system_clock::now();
        //     for (size_t i = 0; i < steps; ++i) {
        //         auto pm_orig = move2dp(p, R);
        //         pm_sum_orig += pm_orig; // just to avoid replacement with opt
        //     }
        //     auto end = std::chrono::system_clock::now();
        //     auto elapsed =
        //         std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        //     fmt::println("");
        //     fmt::println("The measurement orig for point took {}", elapsed);
        //     fmt::println("pm_sum_orig = {}", unitize(pm_sum_orig * (1.0 / steps)));

        //     vec2dp pm_sum{};
        //     start = std::chrono::system_clock::now();
        //     for (size_t i = 0; i < steps; ++i) {
        //         auto pm = move2dp_opt(p, R);
        //         pm_sum += pm; // just to avoid full replacement with opt
        //     }
        //     end = std::chrono::system_clock::now();
        //     elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //     start); fmt::println("The measurement opt for point took {}", elapsed);
        //     fmt::println("pm_sum = {}", unitize(pm_sum * (1.0 / steps)));
        //     fmt::println("");

        //     vec2dp pm_vsum{};
        //     std::vector<vec2dp> pm_vecsum(steps);
        //     for (auto& e : pm_vecsum) {
        //         e = p;
        //     }
        //     start = std::chrono::system_clock::now();
        //     auto pm_vec = move2dp(pm_vecsum, R);
        //     for (size_t i = 0; i < steps; ++i) {
        //         pm_vsum += pm_vec[i]; // just to avoid full replacement with opt
        //     }
        //     end = std::chrono::system_clock::now();
        //     elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //     start); fmt::println("The measurement opt for vec of points took {}",
        //     elapsed); fmt::println("pm_vsum = {}", unitize(pm_vsum * (1.0 / steps)));
        //     fmt::println("");

        //     bivec2dp lm_sum_orig{};
        //     start = std::chrono::system_clock::now();
        //     for (size_t i = 0; i < steps; ++i) {
        //         auto lm_orig = move2dp(l, R);
        //         lm_sum_orig += lm_orig; // just to avoid full replacement with opt
        //     }
        //     end = std::chrono::system_clock::now();
        //     elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //     start); fmt::println("The measurement orig for line took {}", elapsed);
        //     fmt::println("lm_sum_orig = {}", lm_sum_orig * (1.0 / steps));

        //     bivec2dp lm_sum{};
        //     start = std::chrono::system_clock::now();
        //     for (size_t i = 0; i < steps; ++i) {
        //         auto lm = move2dp_opt(l, R);
        //         lm_sum += lm; // just to avoid full replacement with opt
        //     }
        //     end = std::chrono::system_clock::now();
        //     elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        //     start); fmt::println("The measurement opt for line took {}", elapsed);
        //     fmt::println("lm_sum = {}", lm_sum * (1.0 / steps));
        //     fmt::println("");

        //     CHECK(pm_sum_orig * (1.0 / steps) == pm_sum * (1.0 / steps));
        //     CHECK(lm_sum_orig * (1.0 / steps) == lm_sum * (1.0 / steps));

        //     // summary: manual optimization brings benefit in debug-mode.
        //     //          In release-mode with -O3 optimization there is no speed delta
        //     //          => manual optimization does NOT bring benefit
        //     //          => use unoptimized function calls directly (rely on compiler)
        // }

        {
            /////////////////////////////////////////////////////////////////////////////
            // 7th: create the motors directly
            /////////////////////////////////////////////////////////////////////////////

            // define points and lines
            // auto p0 = O_2dp;
            auto p0 = vec2dp{1, 0.5, 1}; // was origin initally

            auto delta = p0 - O_2dp;
            auto p1 = vec2dp{1, 0, 1} + delta;
            auto p2 = vec2dp{1, 1, 1} + delta;
            // auto p2 = vec2dp{1, 0.5, 1} + delta;
            auto p = vec2dp{1, -0.5, 1} + delta; // initial point
            auto pt = vec2dp{0.5, 1, 1} + delta; // target point after applying motor

            auto l1 = unitize(wdg(p0, p1)); // horizontal line
            auto l2 = unitize(wdg(p0, p2)); // line with 45° elevation
            // auto l1 = wdg(p0, p1); // horizontal line
            // auto l2 = wdg(p0, p2); // line with 45° elevation

            auto pis = unitize(rwdg(l1, l2)); // intersection point of lines

            auto R = get_motor_from_lines(l1, l2);
            CHECK(R == unitize(rgpr(l2, l1)));
            auto pfix = unitize(vec2dp{R.c0, R.c1, R.c2});

            auto pm = gr1(rgpr(rgpr(R, p), rrev(R)));  // transformation
            auto pb = gr1(rgpr(rgpr(rrev(R), pm), R)); // reverse transformation

            // fmt::println("");
            // fmt::println("l1:         {:.4g}", l1);
            // fmt::println("bulk:       {:.4g}", bulk(l1));
            // fmt::println("weight:     {:.4g}", weight(l1));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(l1));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(l1));
            // fmt::println("");
            // fmt::println("l2:         {:.4g}", l2);
            // fmt::println("bulk:       {:.4g}", bulk(l2));
            // fmt::println("weight:     {:.4g}", weight(l2));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(l2));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(l2));
            // fmt::println("");
            // fmt::println("R:          {:.4g}", R);
            // fmt::println("bulk:       {:.4g}", bulk(R));
            // fmt::println("weight:     {:.4g}", weight(R));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(R));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(R));
            // fmt::println("");
            // fmt::println("pis:   {:.4g}", pis);
            // fmt::println("pfix:  {:.4g}", pfix);
            // fmt::println("R.c0:  {:.4g}", R.c3);
            // if (R.c3 != 0.0) {
            //     fmt::println("angle: {:.4g}°", rad2deg(std::atan(R.c2 / R.c3)));
            // }
            // else {
            // }
            // fmt::println("");
            // fmt::println("pt: {:.4g}, ptu: {:.4g} (target value)", pt, unitize(pt));
            // fmt::println("pm: {:.4g}, pmu: {:.4g} (after trafo)", pm, unitize(pm));
            // fmt::println("pb: {:.4g}, pbu: {:.4g} (after backtrafo)", pb, unitize(pb));
            // fmt::println("");

            CHECK(unitize(pm) == unitize(pt));
            CHECK(unitize(pb) == unitize(p));
            CHECK(pis == p0); // intersection point is the fix point of the transformation
            CHECK(pfix == p0);

            // now build the motor by hand
            auto S = get_motor(pis, pi / 2.);
            auto ps = move2dp(p, S);

            // fmt::println("");
            // fmt::println("S:          {:.4g}", S);
            // fmt::println("bulk:       {:.4g}", bulk(S));
            // fmt::println("weight:     {:.4g}", weight(S));
            // fmt::println("bulk_nrm:   {:.4g}", bulk_nrm(S));
            // fmt::println("weight_nrm: {:.4g}", weight_nrm(S));
            // fmt::println("");

            CHECK(R == S);
            CHECK(ps == pt);

            auto T = get_motor(delta);
            auto pst = move2dp(p, T);
            fmt::println("");
            fmt::println("p: {:.4g}, pu: {:.4g}", p, unitize(p));
            fmt::println("delta: {:.4g}", delta);
            fmt::println("T: {:.4g}, Tu: {:.4g}", T, unitize(T));
            fmt::println("pst: {:.4g}, pstu: {:.4g} (after trafo)", pst, unitize(pst));
            fmt::println("");

            CHECK(pst == p + delta);
        }
    }

    TEST_CASE("MVec2dp: complement operation")
    {
        fmt::println("MVec2dp: complement operation");

        auto s = scalar2dp{5.0};
        auto v = vec2dp{1.0, 2.0, 3.0};
        auto B = bivec2dp{10.0, 20.0, 30.0};
        auto ps = pscalar2dp{-3.0};

        auto mv = mvec2dp{s, v, B, ps};
        auto mv_e = mvec2dp{s, B};
        auto mv_u = mvec2dp{v, ps};

        auto s2 = scalar2dp{2.5};
        auto v2 = vec2dp{0.5, 1.0, 1.5};
        auto B2 = bivec2dp{5.0, 10.0, 15.0};
        auto ps2 = pscalar2dp{-1.5};

        CHECK(wdg(s, cmpl(s)) == nrm_sq(s) * I_2dp);
        CHECK(wdg(v, cmpl(v)) == nrm_sq(v) * I_2dp);
        CHECK(wdg(B, cmpl(B)) == nrm_sq(B) * I_2dp);
        CHECK(wdg(ps, cmpl(ps)) == nrm_sq(ps) * I_2dp);

        // complement properties
        CHECK(cmpl(cmpl(s)) == s);
        CHECK(cmpl(cmpl(v)) == v);
        CHECK(cmpl(cmpl(B)) == B);
        CHECK(cmpl(cmpl(ps)) == ps);

        // complement values
        CHECK(cmpl(scalar2dp(1.0)) == I_2dp);
        CHECK(cmpl(scalar2dp(1.0)) == rev(scalar2dp(1.0)) * I_2dp);
        CHECK(cmpl(e1_2dp) == e32_2dp);
        CHECK(cmpl(e2_2dp) == -e31_2dp);
        CHECK(cmpl(e3_2dp) == -e12_2dp);
        CHECK(cmpl(e31_2dp) == -e2_2dp);
        CHECK(cmpl(e32_2dp) == e1_2dp);
        CHECK(cmpl(e12_2dp) == -e3_2dp);
        CHECK(cmpl(I_2dp) == scalar2dp(1.0));
        //
        CHECK(cmpl(cmpl(mv)) == mv);
        CHECK(cmpl(cmpl(mv_e)) == mv_e);
        CHECK(cmpl(cmpl(mv_u)) == mv_u);
        //
        CHECK(wdg(s, cmpl(s)) / nrm_sq(s) == I_2dp);
        CHECK(wdg(cmpl(s), s) / nrm_sq(s) == I_2dp);
        CHECK(wdg(v, cmpl(v)) / nrm_sq(v) == I_2dp);
        CHECK(wdg(cmpl(v), v) / nrm_sq(v) == I_2dp);
        CHECK(wdg(B, cmpl(B)) / nrm_sq(B) == I_2dp);
        CHECK(wdg(cmpl(B), B) / nrm_sq(B) == I_2dp);
        CHECK(wdg(ps, cmpl(ps)) / nrm_sq(ps) == I_2dp);
        CHECK(wdg(cmpl(ps), ps) / nrm_sq(ps) == I_2dp);

        // linearity of the complement operation
        double a = 2.0;
        CHECK(cmpl(a * s) == a * cmpl(s));
        CHECK(cmpl(a * v) == a * cmpl(v));
        CHECK(cmpl(a * B) == a * cmpl(B));
        CHECK(cmpl(a * ps) == a * cmpl(ps));

        CHECK(cmpl(s + s2) == cmpl(s) + cmpl(s2));
        CHECK(cmpl(v + v2) == cmpl(v) + cmpl(v2));
        CHECK(cmpl(B + B2) == cmpl(B) + cmpl(B2));
        CHECK(cmpl(ps + ps2) == cmpl(ps) + cmpl(ps2));

        // check regressive wedge product expressed in terms of
        // wegde product and complement operation
        CHECK(rwdg(B, B2) == cmpl(wdg(cmpl(B), cmpl(B2))));
        CHECK(rwdg(B, v) == cmpl(wdg(cmpl(B), cmpl(v))));
        CHECK(rwdg(v, B) == cmpl(wdg(cmpl(v), cmpl(B))));

        CHECK(cmpl(rwdg(B, B2)) == wdg(cmpl(B), cmpl(B2)));
        CHECK(cmpl(wdg(v, v2)) == rwdg(cmpl(v), cmpl(v2)));

        // regressive reverse operation
        CHECK(rrev(s) == cmpl(rev(cmpl(s))));
        CHECK(rrev(v) == cmpl(rev(cmpl(v))));
        CHECK(rrev(B) == cmpl(rev(cmpl(B))));
        CHECK(rrev(ps) == cmpl(rev(cmpl(ps))));

        // origin and horizon
        CHECK(H_2dp == cmpl(O_2dp));
        CHECK(O_2dp == cmpl(H_2dp));
    }

    TEST_CASE("MVec2dp: convenience types")
    {
        fmt::println("MVec2dp: convenience types");

        auto v2d = Vec2d<double>(2.0, 1.5);
        auto p2d = Point2d<double>(2.0, 1.5);
        auto v = Vec2dp<double>(6.0, 4.5, 3.0);
        auto p = Point2dp<double>(6.0, 4.5, 3.0);
        auto p1 = Point2d<double>(1.0, 2.0);
        auto p2 = Point2d<double>(2.0, 3.0);
        auto l1 = wdg(p1, p2);

        // a Point2d is-a Vec2dp, thus all operations defined for Vec2dp should work
        // directly for Point2dp
        CHECK(dot(v, v) == dot(p, p));
        CHECK(p2d == Point2d<double>(v2d));
        CHECK(unitize(p) == Point2dp<double>(p2d));
        CHECK(unitize(p) == Point2dp<double>(v2d));
        CHECK(dot(BiVec2dp<double>(1, 7, 3), BiVec2dp<double>(4, 5, 6)) == 18.0);

        auto vector = Vector2d<double>(1.0, 2.0);

        // fmt::println("vector = {}", vector);

        CHECK(nrm_sq(vector) == 5.0);

        //
        // fmt::println("p1 = {}", p1);
        // fmt::println("p2 = {}", p2);
        // fmt::println("l1 = {}", l1);
        //
        // Line2d ctors
        auto l2 = Line2d<double>(p1, p2);
        auto l3 = Line2d<double>(BiVec2dp<double>{1, 1, -1});
        auto l4 = Line2d<double>(p1, Vec2d<double>{1, 1});
        auto l5 = Line2d<double>(1, 1, -1);
        auto l6 = Line2d<double>();
        CHECK(l1 == l2);
        CHECK(l1 == l3);
        CHECK(l1 == l4);
        CHECK(l1 == l5);
        CHECK(l6 == BiVec2dp<double>(0, 0, 0));
    }

    TEST_CASE("MVec2dp: bulk and weight")
    {
        fmt::println("MVec2dp: bulk and weight");

        auto p1 = vec2dp{-2.0, 1.0, 1.0};
        auto p2 = vec2dp{2.0, 1.0, 2.0};
        auto p3 = vec2dp{2.0, 1.0, 1.0};

        auto l1 = bivec2dp{2.0, 1.0, 1.0};
        auto l2 = bivec2dp{-2.0, 1.0, 3.0};
        auto l3 = wdg(p1, p3);

        // bulk and weight
        CHECK(bulk(p1) == vec2dp{p1.x, p1.y, 0.0});
        CHECK(weight(p1) == vec2dp{0.0, 0.0, p1.z});
        CHECK(bulk(l1) == bivec2dp{0.0, 0.0, l1.z});
        CHECK(weight(l1) == bivec2dp{l1.x, l1.y, 0.0});

        // bulk_nrm_sq and weight_nrm_sq
        CHECK(bulk_nrm_sq(p1) == p1.x * p1.x + p1.y * p1.y);
        CHECK(weight_nrm_sq(p1) == p1.z * p1.z);
        CHECK(bulk_nrm_sq(l1) == l1.z * l1.z);
        CHECK(weight_nrm_sq(l1) == l1.x * l1.x + l1.y * l1.y);

        // bulk_nrm and weight_nrm
        CHECK(bulk_nrm(p1) == std::sqrt(p1.x * p1.x + p1.y * p1.y));
        CHECK(weight_nrm(p1) == std::sqrt(p1.z * p1.z));
        CHECK(bulk_nrm(l1) == std::sqrt(l1.z * l1.z));
        CHECK(weight_nrm(l1) == std::sqrt(l1.x * l1.x + l1.y * l1.y));

        // geom_nrm
        CHECK(geom_nrm(p1).c0 / geom_nrm(p1).c1 ==
              value_t(bulk_nrm(p1)) / value_t(weight_nrm(p1)));
        CHECK(geom_nrm(l1).c0 / geom_nrm(l1).c1 ==
              value_t(bulk_nrm(l1)) / value_t(weight_nrm(l1)));
        CHECK(3.0 * geom_nrm(l1).c0 == geom_nrm(l2).c0);
        CHECK(geom_nrm(l1).c1 == geom_nrm(l2).c1);

        // attitude
        CHECK(att(p1) == p1.z);
        CHECK(att(p1) == rwdg(p1, cmpl(e3_2dp)));
        CHECK(att(p2) == p2.z);
        CHECK(att(p2) == rwdg(p2, cmpl(e3_2dp)));
        CHECK(att(l1) == vec2dp{l1.x, l1.y, 0.0});
        CHECK(att(l1) == rwdg(l1, cmpl(e3_2dp)));
        CHECK(att(l2) == rwdg(l2, cmpl(e3_2dp)));
        CHECK(att(l3) == rwdg(l3, cmpl(e3_2dp)));

        // intersections
        auto l4 = bivec2dp{1.0, 1.0, 3.0};
        auto l5 = bivec2dp{1.0, 1.0, 0.0};
        auto int_sec =
            rwdg(l4, l5); // should intersect at infinity, i.e. att(int_sec) == 0.0
        fmt::println("   att(l4) = {}", att(l4));
        fmt::println("   att(l5) = {}", att(l5));
        fmt::println("   normalize(int_sec(l4,l5)) = {}", normalize(int_sec));
        CHECK(att(int_sec) == 0.0);

        fmt::println("   att(l1) = {}", att(l1));
        fmt::println("   att(l2) = {}", att(l2));
        fmt::println("   sup(l1) = {}", sup(l1));
        fmt::println("   sup(l2) = {}", sup(l2));
        fmt::println("   unitize(rwdg(l1, l2)) = {}", unitize(rwdg(l1, l2)));
        CHECK(unitize(rwdg(l1, l2)) == vec2dp{2.0, 0.5, 1.0});

        auto P1 = vec2dp{1, 1, 1};
        auto P2 = vec2dp{3, 2, 1};
        auto l12 = wdg(P1, P2);

        fmt::println("   P1           = {}", P1);
        fmt::println("   P2           = {}", P2);
        fmt::println("   wdg(P1,P2)   = {}", l12);
        fmt::println("   att(l12)     = {}", att(l12));
        fmt::println("   sup(l12) = {}", sup(l12));
    }

    TEST_CASE("MVec2dp: euclidean distance")
    {
        fmt::println("MVec2dp: euclidean distance");

        auto p1 = vec2dp{-1.0, -1.0, 1.0};
        auto p2 = vec2dp{1.0, -1.0, 1.0};
        auto p3 = vec2dp{1.0, 1.0, 1.0};

        auto l1 = wdg(p1, p2);
        auto l2 = wdg(p2, p3);
        auto l3 = wdg(p1, p3);


        auto dp1p2 = dist2dp(p1, p2);
        auto dp2p3 = dist2dp(p2, p3);
        auto dp1p3 = dist2dp(p1, p3);

        auto dp3l1 = dist2dp(p3, l1);
        auto dp1l2 = dist2dp(p1, l2);
        auto dp2l3 = dist2dp(p2, l3);

        // fmt::println("");
        // fmt::println("dp1p2 = {}, dist = {}", dp1p2, dp1p2.c0 / dp1p2.c1);
        // fmt::println("dp2p3 = {}, dist = {}", dp2p3, dp2p3.c0 / dp2p3.c1);
        // fmt::println("dp1p3 = {}, dist = {}", dp1p3, dp1p3.c0 / dp1p3.c1);
        // fmt::println("");
        // fmt::println("dp3l1 = {}, dist = {}", dp3l1, dp3l1.c0 / dp3l1.c1);
        // fmt::println("dp1l2 = {}, dist = {}", dp1l2, dp1l2.c0 / dp1l2.c1);
        // fmt::println("dp2l3 = {}, dist = {}", dp2l3, dp2l3.c0 / dp2l3.c1);
        // fmt::println("");

        CHECK(dp1p2 == dp2p3);
        CHECK(dp3l1 == dp1l2);
        CHECK(dp1p3.c0 / dp1p3.c1 - 2.0 * std::sqrt(2.0) < eps);
        CHECK(dp2l3.c0 / dp2l3.c1 - std::sqrt(2.0) < eps);

        // auto res = rwdg(pscalar2dp(2.5), H_2dp);
        // fmt::println("res = {}", res);
    }

    TEST_CASE("MVec2dp: product tests")
    {
        fmt::println("MVec2dp: product tests");

        auto s1 = scalar2dp{2.0};
        auto v1 = vec2dp{1.0, -3.0, 0.0};
        // auto v1 = vec2dp{1.0, 0.0, 0.0};
        auto b1 = bivec2dp{2.0, -4.0, 1.0};
        auto ps1 = pscalar2dp{-2.0};

        // auto s2 = scalar2dp{-1.0};
        auto v2 = vec2dp{2.0, 1.5, 0.0};
        // auto v2 = vec2dp{0.0, 1.0, 0.0};
        auto b2 = bivec2dp{-2.0, 10.5, 3.0};
        // auto ps2 = pscalar2dp{3.0};

        // auto s3 = scalar2dp{-10.0};
        auto v3 = vec2dp{-2.0, 4.0, -10.0};
        // auto b3 = bivec2dp{1.0, -2.0, 3.0};
        // auto ps3 = pscalar2dp{-12.0};

        // dot product
        CHECK(dot(scalar2dp(1.0), scalar2dp(1.0)) == scalar2dp(1.0));
        CHECK(dot(e1_2dp, e1_2dp) == scalar2dp(1.0));
        CHECK(dot(e2_2dp, e2_2dp) == scalar2dp(1.0));
        CHECK(dot(e3_2dp, e3_2dp) == scalar2dp(0.0));
        CHECK(dot(e31_2dp, e31_2dp) == scalar2dp(0.0));
        CHECK(dot(e32_2dp, e32_2dp) == scalar2dp(0.0));
        CHECK(dot(e12_2dp, e12_2dp) == scalar2dp(1.0));
        CHECK(dot(pscalar2dp(1.0), pscalar2dp(1.0)) == scalar2dp(0.0));

        // regressive dot product
        CHECK(rdot(scalar2dp(1.0), scalar2dp(1.0)) == pscalar2dp(0.0));
        CHECK(rdot(e1_2dp, e1_2dp) == pscalar2dp(0.0));
        CHECK(rdot(e2_2dp, e2_2dp) == pscalar2dp(0.0));
        CHECK(rdot(e3_2dp, e3_2dp) == pscalar2dp(1.0));
        CHECK(rdot(e31_2dp, e31_2dp) == pscalar2dp(1.0));
        CHECK(rdot(e32_2dp, e32_2dp) == pscalar2dp(1.0));
        CHECK(rdot(e12_2dp, e12_2dp) == pscalar2dp(0.0));
        CHECK(rdot(pscalar2dp(1.0), pscalar2dp(1.0)) == pscalar2dp(1.0));


        // wedge product
        // check full permissible range of arguments, even the ones delivering 0)
        CHECK(wdg(b1, b2) == scalar2dp(0.0));
        CHECK(wdg(v1, ps1) == scalar2dp(0.0));
        CHECK(wdg(ps1, v1) == scalar2dp(0.0));
        CHECK(wdg(b1, ps1) == scalar2dp(0.0));
        CHECK(wdg(ps1, b1) == scalar2dp(0.0));
        CHECK(wdg(ps1, ps1) == scalar2dp(0.0));

        // fmt::println("");
        // fmt::println("   b1           = {}", b1);
        // fmt::println("   b2           = {}", b2);
        // fmt::println("   b1*b2        = {}", b1 * b2);
        // fmt::println("   cmt(b1, b2)  = {}", cmt(b1, b2));
        // fmt::println("   dot(b1,b2)   = {}", dot(b1, b2));
        // fmt::println("");
        CHECK(wdg(b1, b2) == nrm_sq(b1 * b2 - cmt(b1, b2) + dot(b1, b2)));
        CHECK(wdg(v1, ps1) == nrm_sq(0.5 * (v1 * ps1 + rev(ps1) * v1)));
        CHECK(wdg(ps1, v1) == nrm_sq(0.5 * (ps1 * v1 + v1 * rev(ps1))));


        CHECK(rwdg(v1, bulk_dual(wdg(v2, v3))) ==
              rwdg(v1, rwdg(bulk_dual(v2), bulk_dual(v3))));
        CHECK((wdg(v1, v2) << b1) == (v1 << (v2 << b1)));

        // contractions - check full permissible range of arguments, even the ones
        // delivering 0 as a result
        CHECK((ps1 << s1) == 0.0);
        CHECK((ps1 << v1) == 0.0);
        CHECK((ps1 << b1) == 0.0);
        CHECK((b1 << s1) == 0.0);
        CHECK((b1 << v1) == 0.0);
        CHECK((v1 << s1) == 0.0);
        //
        CHECK((s1 >> ps1) == 0.0);
        CHECK((v1 >> ps1) == 0.0);
        CHECK((b1 >> ps1) == 0.0);
        CHECK((s1 >> b1) == 0.0);
        CHECK((v1 >> b1) == 0.0);
        CHECK((s1 >> v1) == 0.0);

        // 2.3.2

        // 3rd with vectors and bivectors directly
        CHECK(v1 * b1 == (b1 >> v1) + wdg(v1, b1));
        CHECK(b1 * v1 == (v1 << b1) + wdg(b1, v1));
        CHECK((b1 >> v1) == -(v1 << gr_inv(b1)));

        // fmt::println("   wdg(v1, b1)     = {}", wdg(v1, b1));
        // fmt::println("   v1 * b1         = {}", v1 * b1);
        // fmt::println("   gr_inv(b1) * v1 = {}", gr_inv(b1) * v1);
        CHECK(wdg(v1, b1) == gr3(0.5 * (v1 * b1 + gr_inv(b1) * v1)));
        CHECK(wdg(b1, v1) == gr3(0.5 * (b1 * v1 + v1 * gr_inv(b1))));

        // cross-check direct implementation of rwdg by comparing with wdg
        CHECK(rwdg(b1, b2) == cmpl(wdg(cmpl(b1), cmpl(b2))));

        // commutator and contraction are equal for a specific case with bivectors
        CHECK(cmt(v1, b1) == (b1 >> v1));
        CHECK(cmt(b1, v1) == (v1 << b1));
        CHECK(cmt(v1, b1) == -cmt(b1, v1));
    }

    TEST_CASE("MVec2dp: simple applications, complements, contraction, expansions")
    {
        fmt::println("MVec2dp: simple applications, complements, contraction, "
                     "expansions");

        auto s1 = scalar2dp{2.0};
        auto v1 = vec2dp{1.0, -3.0, 0.0};
        auto b1 = bivec2dp{2.0, -4.0, 1.0};
        auto ps1 = pscalar2dp{-2.0};
        auto M1 = mvec2dp{s1, v1, b1, ps1};


        // simple projections
        auto v = vec2dp{4.0, 3.5, 0.0};
        auto u = vec2dp{1.0, 2.0, 0.0};
        auto B = e12_2dp;

        auto M2 = mvec2dp{scalar2dp{1.0}, v, B, pscalar2dp{1.0}};

        auto v_in_u = project_onto(v, u);
        auto v_perp_u = reject_from(v, u);

        auto v_in_B = project_onto(v, B);
        auto v_perp_B = reject_from(v, B);

        // fmt::println("   v         = {:g}", v);
        // fmt::println("   u         = {:g}", u);
        // fmt::println("   B         = {:g}", B);
        // fmt::println("   cmpl(v)   = {:g}", cmpl(v));
        // fmt::println("   cmpl(B)   = {:g}", cmpl(B));
        // fmt::println("");
        // fmt::println("   v_in_u    = {:g}", v_in_u);
        // fmt::println("   v_perp_u  = {:g}", v_perp_u);
        // fmt::println("   v << u    = {:g}", v << u);
        // fmt::println("   u >> v    = {:g}", u >> v);
        // fmt::println("");
        // fmt::println("");
        // fmt::println("   v_in_B    = {:g}", v_in_B);
        // fmt::println("   v_perp_B  = {:g}", v_perp_B);
        // fmt::println("   v << B    = {:g}", v << B);
        // fmt::println("   B >> v    = {:g}", B >> v);
        // fmt::println("");

        // x^B = 0 for every point x in B
        CHECK(nrm_sq(wdg(v_in_B, B)) < eps);

        // v_perp_B should be proportional to the normal vector n of B
        // n = att(B)
        // thus, wdg(v_perp_B, att(B)) == 0 is required
        CHECK(nrm_sq(wdg(v_perp_B, att(B))) < eps);

        // v_in_B and v_perp_B should be perpendicular to each other
        CHECK(nrm_sq(dot(v_in_B, v_perp_B)) < eps);

        // v should be the sum of v_in_B and v_perp_B
        CHECK(v == v_in_B + v_perp_B);

        // v should be the sum of v_in_u and v_perp_u
        CHECK(v == v_in_u + v_perp_u);

        // check complements
        CHECK(cmpl(scalar2dp(1.0)) == I_2dp);
        CHECK(cmpl(e1_2dp) == e32_2dp);
        CHECK(cmpl(e2_2dp) == -e31_2dp);
        CHECK(cmpl(e3_2dp) == -e12_2dp);
        CHECK(cmpl(e31_2dp) == -e2_2dp);
        CHECK(cmpl(e32_2dp) == e1_2dp);
        CHECK(cmpl(e12_2dp) == -e3_2dp);
        CHECK(cmpl(I_2dp) == scalar2dp(1.0));
        //
        CHECK(cmpl(cmpl(M1)) == M1);
        CHECK(cmpl(cmpl(mvec2dp_e(s1, b1))) == mvec2dp_e(s1, b1));
        CHECK(cmpl(cmpl(mvec2dp_u(v1, ps1))) == mvec2dp_u(v1, ps1));
        //
        // complements are defined for basis elements only
        // => magnitude has to be covered separately for non-normalized elements

        // left complements = complements in spaces of odd dimension
        CHECK(wdg(cmpl(s1), s1) / nrm_sq(s1) == I_2dp);
        CHECK(wdg(cmpl(v1), v1) / nrm_sq(v1) == I_2dp);
        CHECK(wdg(cmpl(b1), b1) / nrm_sq(b1) == I_2dp);
        CHECK(wdg(cmpl(ps1), ps1) / nrm_sq(ps1) == I_2dp);

        // right complements = complements in spaces of odd dimension
        CHECK(wdg(s1, cmpl(s1)) / nrm_sq(s1) == I_2dp);
        CHECK(wdg(v1, cmpl(v1)) / nrm_sq(v1) == I_2dp);
        CHECK(wdg(b1, cmpl(b1)) / nrm_sq(b1) == I_2dp);
        CHECK(wdg(ps1, cmpl(ps1)) / nrm_sq(ps1) == I_2dp);

        // correspondence of complements with geometric products:
        // bulk_duals differ from complements in pga2dp (influence of degenerate metric)
        // they correspond to complements of the bulk in pga3dp
        CHECK((I_2dp * rev(s1)) == bulk_dual(s1));
        CHECK((I_2dp * rev(v1)) == bulk_dual(v1));
        CHECK((I_2dp * rev(b1)) == bulk_dual(b1));
        CHECK((I_2dp * rev(ps1)) == bulk_dual(ps1));

        CHECK((rev(s1) * I_2dp) == bulk_dual(s1));
        CHECK((rev(v1) * I_2dp) == bulk_dual(v1));
        CHECK((rev(b1) * I_2dp) == bulk_dual(b1));
        CHECK((rev(ps1) * I_2dp) == bulk_dual(ps1));

        // check contractions: <<, >> and rwdg( u, compl(v) )
        // fmt::println("");
        // fmt::println("   v       = {:g}", v);
        // fmt::println("   cmpl(v) = {:g}", cmpl(v));
        // fmt::println("");
        // fmt::println("   v << B  = {:g}", v << B);
        // fmt::println("   B >> v  = {:g}", B >> v);
        // fmt::println("");
        // fmt::println("   B >> v  = rwdg(B, cmpl(v)) = {:g}", rwdg(B, cmpl(v)));
        // fmt::println("   v << B  = rwdg(cmpl(v), B) = {:g}", rwdg(cmpl(v), B));
        // fmt::println("");
        CHECK((B >> v) == rwdg(B, cmpl(v)));
        CHECK((v << B) == rwdg(cmpl(v), B));

        // check expansions: v ^ cmpl(B)  and  cmpl(B) ^ v
        // (create new bivector that contains v and is perpendicular to B)
        //
        // also checks the duality correspondence:
        //      cmpl(v >> B) == v ^ cmpl(B)
        //      cmpl(v >> B) == cmpl(B) ^ v
        // fmt::println("   v      = {:g}", v);
        // fmt::println("   B      = {:g}", B);
        // fmt::println("   v << B = {:g}", v << B);
        // fmt::println("   B >> v = {:g}", B >> v);
        // fmt::println("");
        // fmt::println("   wdg(cmpl(B), v)         = {:g}", wdg(cmpl(B), v));
        // fmt::println("   n=cmpl(wdg(cmpl(B), v)) = {:g}", cmpl(wdg(cmpl(B), v)));
        // fmt::println("");
        // fmt::println("   wdg(v, cmpl(B))         = {:g}", wdg(v, cmpl(B)));
        // fmt::println("   n=cmpl(wdg(v, cmpl(B))) = {:g}", cmpl(wdg(v, cmpl(B))));
        // vector is in plane defined by the expansion
        CHECK(wdg(v, wdg(v, cmpl(B))) == 0.0);
        CHECK(wdg(v, wdg(cmpl(B), v)) == 0.0);
        // duality of the contraction and the wedge product (based on complement)
        CHECK(cmpl(v << B) == wdg(v, cmpl(B)));
        CHECK(cmpl(B >> v) == wdg(cmpl(B), v));

        // check identity with the dot product for same grade vectors
        CHECK((v1 << v) == dot(v1, v));
        CHECK((b1 << B) == dot(b1, B));

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
        left_bulk_contract2dp(M1, M2) = rwdg(cmpl(M1), M2);
        right_bulk_contract2dp(M1, M2) = rwdg(M1, cmpl(M2));
    }

    TEST_CASE("MVec2dp: join and meet (wdg, rwdg)")
    {
        fmt::println("MVec2dp: join and meet (wdg, rwdg)");

        auto p1 = point2d{2, 0};
        auto p2 = point2d{4, 3};
        auto l12 = join(p1, p2); // = wdg(p1, p2)

        auto p3 = point2d{0, 6};
        auto p4 = point2d{2, 6};
        auto l34 = join(p3, p4); // = wdg(p3, p4)

        auto p1p = point2dp{p1}; // all works with projective points as well
        auto p2p = point2dp{p2};
        auto l12p = join(p1p, p2p);
        auto p3p = point2dp{p3};
        auto p4p = point2dp{p4};
        auto l34p = join(p3p, p4p);

        auto p5 = meet(l12, l34);
        auto p5p = meet(l12p, l34p);

        auto p6 = point2d{0, 3}; // the expansion: line perpendicular to l12 through p6
        auto lp6 = expand(p6, l12);

        // fmt::println("");
        // fmt::println("p1 = {}, p2 = {}", p1, p2);
        // fmt::println("l12  = join(p1, p2)   = wdg(p1, p2) = {}", l12);
        // fmt::println("att(l12) = {}", att(l12));
        // fmt::println("");
        // fmt::println("l12p = join(p1p, p2p) = wdg(p1p, p2p) = {}", l12p);
        // fmt::println("att(l12p) = {}", att(l12p));
        // fmt::println("");
        // fmt::println("p5 = meet(l12, l34) = {}", p5);
        // fmt::println("p5p = meet(l12p, l34p) = {}", p5p);
        // fmt::println("");
        // fmt::println("lp6  = {}", lp6);
        // fmt::println("att(lp6) = {}", att(lp6));
        // fmt::println("");

        auto ln1 = unitize(join(point2d{0, 3}, point2d{3, 3}));
        auto ln2 = unitize(join(point2d{0, 6}, point2d{3, 6}));
        auto pt12 = meet(ln1, ln2); // gives direction to intersection at infinity
                                    // weighted by the distance of the lines

        auto ln1p = unitize(join(point2dp{0, 3, 1}, point2dp{3, 3, 1}));
        auto ln2p = unitize(join(point2dp{0, 6, 1}, point2dp{3, 6, 1}));
        auto pt12p = meet(ln1p, ln2p);

        // fmt::println("");
        // fmt::println("ln1      = {}, ln2      = {}", ln1, ln2);
        // fmt::println("att(ln1) = {}, att(ln2) = {}", att(ln1), att(ln2));
        // fmt::println("pt12  = meet(ln1, ln2) = {}, att(pt12) = {}", pt12,
        // att(pt12)); fmt::println(""); fmt::println("ln1p      = {}, ln2p      =
        // {}", ln1p, ln2p); fmt::println("att(ln1p) = {}, att(ln2p) = {}", att(ln1p),
        // att(ln2p)); fmt::println("pt12p  = meet(ln1p, ln2p) = {}", pt12p);
        // fmt::println("");

        CHECK(join(p1, p2) == wdg(p1, p2));
        CHECK(point2dp(p5) == point2dp(p5p));
        CHECK(point2dp(unitize(meet(l12, l34))) == point2dp(unitize(rwdg(l12, l34))));
        CHECK(dot(att(l12), att(lp6)) == 0);

        CHECK(pt12 == pt12p);

        // constructing lines from points and from points and directions
        auto lp = wdg(vec2dp{2, 2, 1}, vec2dp{3, 2, 1}); // line from two points
        auto ld = wdg(vec2dp{2, 2, 1}, vec2dp{1, 0, 0}); // line from a point and a

        // fmt::println("");
        // fmt::println("lp      = {}, ld      = {}", lp, ld);
        // fmt::println("");

        CHECK(lp == ld); // both lines are equivalent
    }

    TEST_CASE("G<2,0,1>: Scalar2dp and PScalar2dp formatting tests")
    {
        fmt::println("G<2,0,1>: Scalar2dp and PScalar2dp formatting tests");

        // Test scalar formatting with different format specifiers
        scalar2dp s{3.14159};
        pscalar2dp ps{2.71828};

        // Default formatting
        std::string default_scalar = fmt::format("{}", s);
        std::string default_pscalar = fmt::format("{}", ps);
        CHECK(default_scalar == "Scalar2dp(3.14159)");
        CHECK(default_pscalar == "PScalar2dp(2.71828)");

        // Two decimal places
        std::string two_decimals_scalar = fmt::format("{:.2f}", s);
        std::string two_decimals_pscalar = fmt::format("{:.2f}", ps);
        CHECK(two_decimals_scalar == "Scalar2dp(3.14)");
        CHECK(two_decimals_pscalar == "PScalar2dp(2.72)");

        // Scientific notation
        std::string scientific_scalar = fmt::format("{:.2e}", s);
        std::string scientific_pscalar = fmt::format("{:.2e}", ps);
        CHECK(scientific_scalar == "Scalar2dp(3.14e+00)");
        CHECK(scientific_pscalar == "PScalar2dp(2.72e+00)");

        // Six decimal places
        std::string six_decimals_scalar = fmt::format("{:.6f}", s);
        std::string six_decimals_pscalar = fmt::format("{:.6f}", ps);
        CHECK(six_decimals_scalar == "Scalar2dp(3.141590)");
        CHECK(six_decimals_pscalar == "PScalar2dp(2.718280)");

        fmt::println(
            "   All Scalar2dp and PScalar2dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("Vec2dp: comprehensive format specifier tests with pre-computed values")
    {
        fmt::println(
            "Vec2dp: comprehensive format specifier tests with pre-computed values:\n");

        // Test vector with specific values for predictable output
        vec2dp v{3.14159, 2.71828, 1.41421};

        // Default formatting - using value_t precision
        std::string default_format = fmt::format("{}", v);
        CHECK(default_format == "Vec2dp(3.14159,2.71828,1.41421)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", v);
        CHECK(two_decimals == "Vec2dp(3.14,2.72,1.41)");

        // Scientific notation with 2 decimal places
        std::string scientific = fmt::format("{:.2e}", v);
        CHECK(scientific == "Vec2dp(3.14e+00,2.72e+00,1.41e+00)");

        // Six decimal places
        std::string six_decimals = fmt::format("{:.6f}", v);
        CHECK(six_decimals == "Vec2dp(3.141590,2.718280,1.414210)");

        // Test with different vector values for edge cases
        vec2dp v_zero{0.0, 0.0, 0.0};
        vec2dp v_negative{-1.5, -2.5, -3.5};

        std::string zero_format = fmt::format("{:.1f}", v_zero);
        std::string negative_format = fmt::format("{:.1f}", v_negative);
        CHECK(zero_format == "Vec2dp(0.0,0.0,0.0)");
        CHECK(negative_format == "Vec2dp(-1.5,-2.5,-3.5)");

        fmt::println("   All Vec2dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("G<2,0,1>: BiVec2dp formatting tests")
    {
        fmt::println("G<2,0,1>: BiVec2dp formatting tests");

        // Test bivector formatting with different format specifiers
        bivec2dp bv{1.234, 5.678, 9.012};

        // Default formatting
        std::string default_format = fmt::format("{}", bv);
        CHECK(default_format == "BiVec2dp(1.234,5.678,9.012)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", bv);
        CHECK(two_decimals == "BiVec2dp(1.23,5.68,9.01)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", bv);
        CHECK(scientific == "BiVec2dp(1.23e+00,5.68e+00,9.01e+00)");

        // Four decimal places
        std::string four_decimals = fmt::format("{:.4f}", bv);
        CHECK(four_decimals == "BiVec2dp(1.2340,5.6780,9.0120)");

        fmt::println("   All BiVec2dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - default init")
    {
        fmt::println("G<2,0,1>: DualNum2dp - default init");
        dualnum2dp dn;
        CHECK(abs(dn.c0) < eps);
        CHECK(abs(dn.c1) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - constructor from two components")
    {
        fmt::println("G<2,0,1>: DualNum2dp - constructor from two components");
        dualnum2dp dn{2.5, 3.7};
        CHECK(abs(dn.c0 - 2.5) < eps);
        CHECK(abs(dn.c1 - 3.7) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - constructor from scalar")
    {
        fmt::println("G<2,0,1>: DualNum2dp - constructor from scalar");
        scalar2dp s{5.0};
        dualnum2dp dn{s};
        CHECK(abs(dn.c0 - 5.0) < eps);
        CHECK(abs(dn.c1) < eps);
        CHECK(gr0(dn) == s);
        CHECK(gr3(dn) == pscalar2dp{0.0});
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - constructor from pseudoscalar")
    {
        fmt::println("G<2,0,1>: DualNum2dp - constructor from pseudoscalar");
        pscalar2dp ps{7.5};
        dualnum2dp dn{ps};
        CHECK(abs(dn.c0) < eps);
        CHECK(abs(dn.c1 - 7.5) < eps);
        CHECK(gr0(dn) == scalar2dp{0.0});
        CHECK(gr3(dn) == ps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - constructor from scalar and pseudoscalar")
    {
        fmt::println("G<2,0,1>: DualNum2dp - constructor from scalar and pseudoscalar");
        scalar2dp s{3.0};
        pscalar2dp ps{4.0};
        dualnum2dp dn{s, ps};
        CHECK(abs(dn.c0 - 3.0) < eps);
        CHECK(abs(dn.c1 - 4.0) < eps);
        CHECK(gr0(dn) == s);
        CHECK(gr3(dn) == ps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - copy constructor and assignment")
    {
        fmt::println("G<2,0,1>: DualNum2dp - copy constructor and assignment");
        dualnum2dp dn1{2.5, 3.7};
        dualnum2dp dn2{dn1};   // copy constructor
        dualnum2dp dn3 = dn1;  // copy assignment
        dualnum2dp dn4 = -dn1; // unary minus

        CHECK(abs(dn2.c0 - 2.5) < eps);
        CHECK(abs(dn2.c1 - 3.7) < eps);
        CHECK(abs(dn3.c0 - 2.5) < eps);
        CHECK(abs(dn3.c1 - 3.7) < eps);
        CHECK(abs(dn4.c0 + 2.5) < eps);
        CHECK(abs(dn4.c1 + 3.7) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - grade extraction")
    {
        fmt::println("G<2,0,1>: DualNum2dp - grade extraction");
        dualnum2dp dn{5.0, 7.0};

        scalar2dp s = gr0(dn);
        pscalar2dp ps = gr3(dn);

        CHECK(abs(value_t(s) - 5.0) < eps);
        CHECK(abs(value_t(ps) - 7.0) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - addition operations")
    {
        fmt::println("G<2,0,1>: DualNum2dp - addition operations");

        scalar2dp s{2.0};
        pscalar2dp ps{3.0};
        dualnum2dp dn1{1.0, 2.0};
        dualnum2dp dn2{3.0, 4.0};

        // scalar + pseudoscalar => dual number (explicit construction to avoid ambiguity)
        // MSVC workaround: explicit cast to value_t for cross-platform compatibility
        dualnum2dp result1{s, ps};
        CHECK(abs(result1.c0 - 2.0) < eps);
        CHECK(abs(result1.c1 - 3.0) < eps);

        // scalar + dual number
        dualnum2dp result3 = s + dn1;
        CHECK(abs(result3.c0 - 3.0) < eps);
        CHECK(abs(result3.c1 - 2.0) < eps);

        // dual number + scalar
        dualnum2dp result4 = dn1 + s;
        CHECK(abs(result4.c0 - 3.0) < eps);
        CHECK(abs(result4.c1 - 2.0) < eps);

        // pseudoscalar + dual number
        dualnum2dp result5 = ps + dn1;
        CHECK(abs(result5.c0 - 1.0) < eps);
        CHECK(abs(result5.c1 - 5.0) < eps);

        // dual number + pseudoscalar
        dualnum2dp result6 = dn1 + ps;
        CHECK(abs(result6.c0 - 1.0) < eps);
        CHECK(abs(result6.c1 - 5.0) < eps);

        // dual number + dual number (using compound assignment)
        dualnum2dp dn3 = dn1;
        dn3 += dn2;
        CHECK(abs(dn3.c0 - 4.0) < eps);
        CHECK(abs(dn3.c1 - 6.0) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - subtraction operations")
    {
        fmt::println("G<2,0,1>: DualNum2dp - subtraction operations");

        scalar2dp s{5.0};
        pscalar2dp ps{3.0};
        dualnum2dp dn1{7.0, 4.0};
        dualnum2dp dn2{2.0, 1.0};

        // scalar - dual number
        dualnum2dp result3 = s - dn1;
        CHECK(abs(result3.c0 + 2.0) < eps);
        CHECK(abs(result3.c1 + 4.0) < eps);

        // dual number - scalar
        dualnum2dp result4 = dn1 - s;
        CHECK(abs(result4.c0 - 2.0) < eps);
        CHECK(abs(result4.c1 - 4.0) < eps);

        // pseudoscalar - dual number
        dualnum2dp result5 = ps - dn1;
        CHECK(abs(result5.c0 + 7.0) < eps);
        CHECK(abs(result5.c1 + 1.0) < eps);

        // dual number - pseudoscalar
        dualnum2dp result6 = dn1 - ps;
        CHECK(abs(result6.c0 - 7.0) < eps);
        CHECK(abs(result6.c1 - 1.0) < eps);

        // dual number - dual number (using compound assignment)
        dualnum2dp dn3 = dn1;
        dn3 -= dn2;
        CHECK(abs(dn3.c0 - 5.0) < eps);
        CHECK(abs(dn3.c1 - 3.0) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - scalar multiplication and division")
    {
        fmt::println("G<2,0,1>: DualNum2dp - scalar multiplication and division");

        dualnum2dp dn{6.0, 8.0};

        // scalar multiplication
        dualnum2dp dn2 = dn;
        dn2 *= 2.0;
        CHECK(abs(dn2.c0 - 12.0) < eps);
        CHECK(abs(dn2.c1 - 16.0) < eps);

        // scalar division
        dualnum2dp dn3 = dn;
        dn3 /= 2.0;
        CHECK(abs(dn3.c0 - 3.0) < eps);
        CHECK(abs(dn3.c1 - 4.0) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - compound assignment operators")
    {
        fmt::println("G<2,0,1>: DualNum2dp - compound assignment operators");

        dualnum2dp dn1{10.0, 20.0};
        dualnum2dp dn2{3.0, 5.0};

        // test +=
        dualnum2dp result1 = dn1;
        result1 += dn2;
        CHECK(abs(result1.c0 - 13.0) < eps);
        CHECK(abs(result1.c1 - 25.0) < eps);

        // test -=
        dualnum2dp result2 = dn1;
        result2 -= dn2;
        CHECK(abs(result2.c0 - 7.0) < eps);
        CHECK(abs(result2.c1 - 15.0) < eps);

        // test *=
        dualnum2dp result3 = dn1;
        result3 *= 0.5;
        CHECK(abs(result3.c0 - 5.0) < eps);
        CHECK(abs(result3.c1 - 10.0) < eps);

        // test /=
        dualnum2dp result4 = dn1;
        result4 /= 2.0;
        CHECK(abs(result4.c0 - 5.0) < eps);
        CHECK(abs(result4.c1 - 10.0) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - mixed type operations")
    {
        fmt::println("G<2,0,1>: DualNum2dp - mixed type operations");

        // Test that operations work with different floating point types
        scalar2dp s{2.0};
        pscalar2dp ps{3.0};
        dualnum2dp dn{1.5, 2.5};

        // Test dual number operations with scalar and pseudoscalar
        dualnum2dp result2 = dn + s;
        CHECK(abs(result2.c0 - 3.5) < eps);
        CHECK(abs(result2.c1 - 2.5) < eps);

        dualnum2dp result3 = ps + dn;
        CHECK(abs(result3.c0 - 1.5) < eps);
        CHECK(abs(result3.c1 - 5.5) < eps);

        dualnum2dp result4 = dn - s;
        CHECK(abs(result4.c0 + 0.5) < eps);
        CHECK(abs(result4.c1 - 2.5) < eps);
    }

    TEST_CASE("G<2,0,1>: DualNum2dp - property: I_2dp^2 = 0")
    {
        fmt::println("G<2,0,1>: DualNum2dp - property: I_2dp^2 = 0");

        // The pseudoscalar in PGA2DP squares to zero (degenerate metric)
        // For dual numbers: (c0 + c1*I)(c0 + c1*I) = c0^2 + 2*c0*c1*I + c1^2*I^2
        // Since I^2 = 0: result = c0^2 + 2*c0*c1*I

        // Verify that I_2dp^2 = 0 through geometric product (using operator*)
        auto I_squared = I_2dp * I_2dp;
        CHECK(abs(value_t(I_squared)) < eps);

        fmt::println("I_2dp^2 = 0 (degenerate metric property verified)");
    }

    TEST_CASE("G<2,0,1>: DualNum2dp formatting tests")
    {
        fmt::println("G<2,0,1>: DualNum2dp formatting tests");

        // Test dual number formatting with different format specifiers
        dualnum2dp dn{2.5, 3.7};

        // Default formatting
        std::string default_format = fmt::format("{}", dn);
        CHECK(default_format == "DualNum2dp(2.5,3.7)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", dn);
        CHECK(two_decimals == "DualNum2dp(2.50,3.70)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", dn);
        CHECK(scientific == "DualNum2dp(2.50e+00,3.70e+00)");

        // Three decimal places
        std::string three_decimals = fmt::format("{:.3f}", dn);
        CHECK(three_decimals == "DualNum2dp(2.500,3.700)");

        fmt::println("   All DualNum2dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("G<2,0,1>: MVec2dp_E, MVec2dp_U, and MVec2dp formatting tests")
    {
        fmt::println("G<2,0,1>: MVec2dp_E, MVec2dp_U, and MVec2dp formatting tests");

        // Test multivector even formatting
        mvec2dp_e me{1.1, 2.2, 3.3, 4.4};
        std::string even_default = fmt::format("{}", me);
        std::string even_two_dec = fmt::format("{:.2f}", me);
        CHECK(even_default == "MVec2dp_E(1.1,2.2,3.3,4.4)");
        CHECK(even_two_dec == "MVec2dp_E(1.10,2.20,3.30,4.40)");

        // Test multivector odd formatting
        mvec2dp_u mu{1.5, 2.5, 3.5, 4.5};
        std::string odd_default = fmt::format("{}", mu);
        std::string odd_two_dec = fmt::format("{:.2f}", mu);
        CHECK(odd_default == "MVec2dp_U(1.5,2.5,3.5,4.5)");
        CHECK(odd_two_dec == "MVec2dp_U(1.50,2.50,3.50,4.50)");

        // Test full multivector formatting
        mvec2dp m{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        std::string full_default = fmt::format("{}", m);
        std::string full_scientific = fmt::format("{:.1e}", m);
        CHECK(full_default == "MVec2dp(1,2,3,4,5,6,7,8)");
        CHECK(full_scientific == "MVec2dp(1.0e+00,2.0e+00,3.0e+00,4.0e+00,5.0e+00,"
                                 "6.0e+00,7.0e+00,8.0e+00)");

        fmt::println("   All MVec2dp format tests passed with expected values!");
        fmt::println("");
    }

    TEST_CASE("PGA2DP: congruence tests")
    {
        fmt::println("PGA2DP: congruence tests");

        // Test scalars
        scalar2dp s1{5.0};
        scalar2dp s2{-3.0}; // different sign
        scalar2dp s3{2.5};  // same sign as s1
        scalar2dp s_zero{0.0};

        CHECK(is_congruent(s1, s2) == true);         // different signs are congruent
        CHECK(is_congruent(s1, s3) == true);         // same signs are congruent
        CHECK(is_congruent(s1, s_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(s_zero, s_zero) == true); // zero vs zero

        // Test vectors (points in PGA2DP)
        vec2dp v1{1.0, 0.0, 1.0};   // point on x-axis
        vec2dp v2{0.0, 1.0, 1.0};   // point on y-axis
        vec2dp v3{2.0, 0.0, 2.0};   // parallel to v1 (same direction from origin)
        vec2dp v4{-1.0, 0.0, -1.0}; // antiparallel to v1
        vec2dp v_zero{0.0, 0.0, 0.0};

        CHECK(is_congruent(v1, v2) == false);        // different directions
        CHECK(is_congruent(v1, v3) == true);         // parallel directions
        CHECK(is_congruent(v1, v4) == true);         // antiparallel directions
        CHECK(is_congruent(v1, v_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(v_zero, v_zero) == true); // zero vs zero

        // Test bivectors (lines in PGA2DP)
        bivec2dp b1{1.0, 0.0, 0.0};  // line parallel to y-axis
        bivec2dp b2{0.0, 1.0, 0.0};  // line parallel to x-axis
        bivec2dp b3{2.0, 0.0, 0.0};  // parallel to b1
        bivec2dp b4{-1.5, 0.0, 0.0}; // antiparallel to b1
        bivec2dp b_zero{0.0, 0.0, 0.0};

        CHECK(is_congruent(b1, b2) == false);        // different line orientations
        CHECK(is_congruent(b1, b3) == true);         // same line orientation
        CHECK(is_congruent(b1, b4) == true);         // same line, opposite orientation
        CHECK(is_congruent(b1, b_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(b_zero, b_zero) == true); // zero vs zero

        // Test pseudoscalars (max grade in PGA2DP)
        pscalar2dp p1{1.0};
        pscalar2dp p2{2.0};  // different magnitude
        pscalar2dp p3{-1.0}; // opposite sign
        pscalar2dp p_zero{0.0};

        CHECK(is_congruent(p1, p2) == true);      // all non-zero pseudoscalars congruent
        CHECK(is_congruent(p1, p3) == true);      // different signs still congruent
        CHECK(is_congruent(p1, p_zero) == false); // zero vs non-zero
        CHECK(is_congruent(p_zero, p_zero) == true); // zero vs zero

        // fmt::println("   All PGA2DP congruence tests passed");
    }

    TEST_CASE("G<2,0,1>: transwedge and regressive transwedge products")
    {
        fmt::println("G<2,0,1>: transwedge and regressive transwedge products");

        auto s = scalar2dp{4.0};
        auto v = vec2dp{1.0, 2.0, 1.0};
        auto v2 = vec2dp{-1.0, 3.0, 1.0};
        auto B = bivec2dp{10.0, 20.0, 30.0};
        auto B2 = bivec2dp{-20.0, -40.0, -60.0};
        auto ps = pscalar2dp{-3.0};

        // check identities of transwedge product
        CHECK(twdg1(ps, v) == ps * v);
        CHECK(twdg1(ps, v) == (v << ps));
        CHECK(twdg1(ps, v) == (ps >> v));
        CHECK(twdg1(ps, v) ==
              rwdg(ps, bulk_dual(v))); // contraction expressed as rwdg prod.

        CHECK(twdg1(ps, v) == twdg1(v, ps));

        CHECK(twdg1(v, ps) == v * ps);
        CHECK(twdg1(v, ps) == (ps >> v));
        CHECK(twdg1(v, ps) == (v << ps));
        CHECK(twdg1(v, ps) ==
              rwdg(bulk_dual(v), ps)); // contraction expressed as rwdg prod.

        CHECK(twdg1(B, B2) == cmt(B, B2));

        CHECK(twdg1(B, v) == (v << B));
        CHECK(twdg1(v, B) == (B >> v));

        CHECK(twdg1(v, v2) == dot(v, v2));
        CHECK(twdg1(v, v2) == dot(v2, v));

        // check identities of regressive transwedge product
        CHECK(rtwdg1(B, B2) == rdot(B, B2));
        CHECK(rtwdg1(B, B2) == rdot(B2, B));

        CHECK(rtwdg1(B, v) == wdg(v, weight_dual(B))); // expansion expressed as wdg prod.
        CHECK(rtwdg1(B, v) ==
              right_weight_expand2dp(v, B)); // expansion expressed as wdg prod.
        CHECK(rtwdg1(B, v) == rcmt(B, v));
        CHECK(rtwdg1(B, v) == -rtwdg1(v, B));          // vector ans bivector anticommute
        CHECK(rtwdg1(v, B) == wdg(weight_dual(B), v)); // expansion expressed as wdg prod.
        CHECK(rtwdg1(v, B) ==
              left_weight_expand2dp(B, v)); // expansion expressed as wdg prod.
        CHECK(rtwdg1(v, B) == rcmt(v, B));

        CHECK(rtwdg1(B, s) == wdg(s, weight_dual(B))); // expansion expressed as wdg prod.
        CHECK(rtwdg1(B, s) ==
              right_weight_expand2dp(s, B));           // expansion expressed as wdg prod.
        CHECK(rtwdg1(s, B) == rtwdg1(B, s));           // bivector and scalar commute
        CHECK(rtwdg1(s, B) == wdg(weight_dual(B), s)); // expansion expressed as wdg prod.
        CHECK(rtwdg1(s, B) ==
              left_weight_expand2dp(B, s)); // expansion expressed as wdg prod.

        CHECK(rtwdg1(v, v2) == rcmt(v, v2)); // identity to regressive commutator product

        fmt::println("");
    }

    TEST_CASE("pga2dp: application tests - force and moment (I/II)")
    {
        fmt::println("pga2dp: application tests - force and moment (I/II)");


        // PGA implementation
        auto P = vec2dp{2.25, 1, 1};
        auto f = vec2dp{0, 0.75, 0};
        auto F = wdg(P, f); // force line resulting from f acting at P (=forque)

        fmt::println("P = {}", P);
        fmt::println("f = {}", f);
        fmt::println("F = P^f = {}, att(f) = {}, M_O = bulk(F) = {}", F, att(F), bulk(F));
        fmt::println("");
        fmt::println("O_2dp^f = {}, att(O_2dp^f) = {} (=force)", wdg(O_2dp, f),
                     att(wdg(O_2dp, f)));
        fmt::println("p^f = {}, bulk(p^f) = {} (=torque)", wdg(P - O_2dp, f),
                     bulk(wdg(P - O_2dp, f)));
        fmt::println("");
        CHECK(sup(F) == vec2dp{2.25, 0, 1});
        CHECK(wdg(P, f) == wdg(O_2dp, f) + wdg(P - O_2dp, f));

        auto R1 = vec2dp{1.5, 2, 1};
        auto R2 = vec2dp{3.0, 2, 1};
        auto R3 = vec2dp{3.0, 0, 1};
        auto R4 = vec2dp{1.5, 0, 1};
        auto F_R1 = wdg(P - R1, f);
        auto F_R2 = wdg(P - R2, f);
        auto F_R3 = wdg(P - R3, f);
        auto F_R4 = wdg(P - R4, f);

        fmt::println("R1 = {}", R1);
        fmt::println("R2 = {}", R2);
        fmt::println("R3 = {}", R3);
        fmt::println("R4 = {}", R4);
        fmt::println("");
        fmt::println("F_R1 = (P - R1)^f = {}, att(f) = {}, M_R1 = bulk(F_R1) = {}", F_R1,
                     att(F_R1), bulk(F_R1));
        fmt::println("F_R2 = (P - R2)^f = {}, att(f) = {}, M_R2 = bulk(F_R2) = {}", F_R2,
                     att(F_R2), bulk(F_R2));
        fmt::println("F_R3 = (P - R3)^f = {}, att(f) = {}, M_R3 = bulk(F_R3) = {}", F_R3,
                     att(F_R3), bulk(F_R3));
        fmt::println("F_R4 = (P - R4)^f = {}, att(f) = {}, M_R4 = bulk(F_R4) = {}", F_R4,
                     att(F_R4), bulk(F_R4));
        fmt::println("");

        CHECK(F_R1 == F - wdg(R1, f));
        CHECK(F_R2 == F - wdg(R2, f));
        CHECK(F_R3 == F - wdg(R3, f));
        CHECK(F_R4 == F - wdg(R4, f));

        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Complement Mathematical Properties Tests
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<2,0,1>: complement involution properties")
    {
        fmt::println("G<2,0,1>: complement involution properties");

        // For odd-dimensional algebras: cmpl(cmpl(u)) = +u
        scalar2dp s{3.0};
        vec2dp v{2.0, 5.0, 1.0};
        bivec2dp B{1.0, 2.0, 3.0};
        pscalar2dp ps{11.0};
        mvec2dp_e M_e{s, B};
        mvec2dp_u M_u{v, ps};
        mvec2dp M{s, v, B, ps};

        CHECK(cmpl(cmpl(s)) == s);
        CHECK(cmpl(cmpl(v)) == v);
        CHECK(cmpl(cmpl(B)) == B);
        CHECK(cmpl(cmpl(ps)) == ps);
        CHECK(cmpl(cmpl(M_e)) == M_e);
        CHECK(cmpl(cmpl(M_u)) == M_u);
        CHECK(cmpl(cmpl(M)) == M);

        fmt::println("complement involution: cmpl(cmpl(u)) = u");
        fmt::println("");
    }

    TEST_CASE("G<2,0,1>: extended metric matrix validation")
    {
        fmt::println("G<2,0,1>: extended metric matrix validation");

        // Get the extended metric matrix from stored constants
        auto G = pga2dp_metric_view();

        // PGA2DP has degenerate metric: e1²=1, e2²=1, e3²=0
        // Basis order: 1, e1, e2, e3, e23, e31, e12, e321

        // Verify expected diagonal values
        CHECK(G[0, 0] == 1); // scalar: always 1
        CHECK(G[1, 1] == 1); // e1: +1 (Euclidean)
        CHECK(G[2, 2] == 1); // e2: +1 (Euclidean)
        CHECK(G[3, 3] == 0); // e3: 0 (null/degenerate)
        CHECK(G[4, 4] == 0); // e23: involves e3, so 0
        CHECK(G[5, 5] == 0); // e31: involves e3, so 0
        CHECK(G[6, 6] == 1); // e12: both Euclidean, so +1
        CHECK(G[7, 7] == 0); // e321: involves e3, so 0

        // Verify off-diagonal elements are zero (orthogonal basis)
        for (size_t i = 0; i < 8; ++i) {
            for (size_t j = 0; j < 8; ++j) {
                if (i != j) {
                    CHECK(G[i, j] == 0);
                }
            }
        }

        fmt::println("Extended metric diagonal correct for PGA2DP (degenerate e3)");
    }

    TEST_CASE("G<2,0,1>: extended metric recursive extraction via wedge products")
    {
        fmt::println("G<2,0,1>: extended metric recursive extraction via wedge products");

        auto G = pga2dp_metric_view();

        // Level 0: Scalar (always 1)
        CHECK(G[0, 0] == 1);

        // Level 1: Vectors (extract from dot products)
        // Note: e3 is the null vector (origin in PGA2DP)
        value_t g_e1 = value_t(dot(e1_2dp, e1_2dp));
        value_t g_e2 = value_t(dot(e2_2dp, e2_2dp));
        value_t g_e3 = value_t(dot(e3_2dp, e3_2dp)); // Should be 0 (null)
        CHECK(abs(g_e1 - value_t(G[1, 1])) < eps);
        CHECK(abs(g_e2 - value_t(G[2, 2])) < eps);
        CHECK(abs(g_e3 - value_t(G[3, 3])) < eps);

        // Level 2: Bivectors (extract from wedge products + dot)
        // Basis order: e31, e32, e12
        auto e31_constructed = wdg(e3_2dp, e1_2dp);
        auto e32_constructed = wdg(e3_2dp, e2_2dp);
        auto e12_constructed = wdg(e1_2dp, e2_2dp);

        value_t g_e31 =
            value_t(dot(e31_constructed, e31_constructed)); // Should be 0 (involves e3)
        value_t g_e32 =
            value_t(dot(e32_constructed, e32_constructed)); // Should be 0 (involves e3)
        value_t g_e12 = value_t(
            dot(e12_constructed, e12_constructed)); // Should be 1 (both Euclidean)

        CHECK(abs(g_e31 - value_t(G[4, 4])) < eps);
        CHECK(abs(g_e32 - value_t(G[5, 5])) < eps);
        CHECK(abs(g_e12 - value_t(G[6, 6])) < eps);

        // Verify constructed bivectors match canonical ones
        CHECK(e31_constructed == e31_2dp);
        CHECK(e32_constructed == e32_2dp);
        CHECK(e12_constructed == e12_2dp);

        // Level 3: Trivector/Pseudoscalar (extract from wedge products + dot)
        auto e321_constructed = wdg(wdg(e3_2dp, e2_2dp), e1_2dp);
        value_t g_e321 =
            value_t(dot(e321_constructed, e321_constructed)); // Should be 0 (involves e3)
        CHECK(abs(g_e321 - value_t(G[7, 7])) < eps);
        CHECK(e321_constructed == e321_2dp);

        fmt::println("Recursive extraction: scalar → vectors → bivectors → trivector "
                     "all match (including null dimension)");
    }

    TEST_CASE("G<2,0,1>: extended metric vs dot products")
    {
        fmt::println("G<2,0,1>: extended metric vs dot products");

        // The extended metric diagonal should match dot(basis, basis)
        auto G = pga2dp_metric_view();

        // Scalars: always 1
        CHECK(G[0, 0] == 1);

        // Vectors: e1·e1=1, e2·e2=1, e3·e3=0
        CHECK(abs(value_t(dot(e1_2dp, e1_2dp)) - value_t(G[1, 1])) < eps);
        CHECK(abs(value_t(dot(e2_2dp, e2_2dp)) - value_t(G[2, 2])) < eps);
        CHECK(abs(value_t(dot(e3_2dp, e3_2dp)) - value_t(G[3, 3])) < eps);

        // Bivectors: e23·e23=0, e31·e31=0, e12·e12=1
        CHECK(abs(value_t(dot(e31_2dp, e31_2dp)) - value_t(G[4, 4])) < eps);
        CHECK(abs(value_t(dot(e32_2dp, e32_2dp)) - value_t(G[5, 5])) < eps);
        CHECK(abs(value_t(dot(e12_2dp, e12_2dp)) - value_t(G[6, 6])) < eps);

        // Pseudoscalar: e321·e321=0 (involves null dimension)
        CHECK(abs(value_t(dot(e321_2dp, e321_2dp)) - value_t(G[7, 7])) < eps);

        fmt::println("Extended metric diagonal matches dot products (including null e3)");
    }

    TEST_CASE("G<2,0,1>: exponential function")
    {
        fmt::println("G<2,0,1>: exponential function");

        auto X0 = vec2dp{1, 0, 1};
        auto delta = vec2dp{-2, -1, 0};
        auto X = X0 + delta;

        auto B = att(bulk_dual(delta)); // vector
        auto M = get_motor(delta);      // translation motor

        fmt::println("B_tra                 = {}", B);
        fmt::println("M (from get_motor())  = {}", M);
        fmt::println("M (from exp(0.5 * B)) = {}", exp(0.5 * B));
        fmt::println("");

        CHECK(X == move2dp(X0, M)); // motor moves to X
        CHECK(M == exp(0.5 * B));   // bivector creates motor via exp()

        auto P_fix = vec2dp{-2, -3, 1};
        auto X_rot = vec2dp{-2, -3 + std::sqrt(18), 1};

        B = P_fix * deg2rad(45);
        M = get_motor(P_fix, deg2rad(45));

        fmt::println("B_rot                 = {}", B);
        fmt::println("M (from get_motor())  = {}", M);
        fmt::println("M (from exp(0.5 * B)) = {}", exp(0.5 * B));
        fmt::println("");

        CHECK(X_rot == move2dp(X0, M));
        CHECK(M == exp(0.5 * B)); // bivector creates motor via exp()

        auto B0 = vec2dp{0, 0, deg2rad(45)};
        auto M0 = get_motor(O_2dp, deg2rad(45));
        fmt::println("B0 = {}, exp(0.5 * B0) = {}", B0, exp(0.5 * B0));
        fmt::println("M0 = {}", M0);
        CHECK(M0 == exp(0.5 * B0));
        CHECK(move2dp(X0, M0) == vec2dp{0.5 * std::sqrt(2.0), 0.5 * std::sqrt(2), 1});
        CHECK(move2dp(O_2dp, M0) == O_2dp); // fixed point
        //
        auto B1 = vec2dp{1, 1, 1} * deg2rad(45);
        auto M1 = get_motor(vec2dp{1, 1, 1}, deg2rad(45));
        fmt::println("B1 = {}, exp(0.5 * B1) = {}", B1, exp(0.5 * B1));
        fmt::println("M1 = {}", M1);
        CHECK(M1 == exp(0.5 * B1));
        CHECK(move2dp(X0, M1) ==
              vec2dp{1.0 + 0.5 * std::sqrt(2.0), 1.0 - 0.5 * std::sqrt(2), 1});
        CHECK(move2dp(vec2dp{1, 1, 1}, M1) == vec2dp{1, 1, 1}); // fixed point
        fmt::println("");

        // now as stepwise approach with concatenated transformations
        auto Ma = get_motor(vec2dp{-1, -1, 1});  // move from fixed point to O_2dp
        auto Mb = get_motor(O_2dp, deg2rad(45)); // rotate around O_2dp
        auto Mc = get_motor(vec2dp{1, 1, 1});    // move back from O_2dp to fixed point
        auto M_res = rgpr(rgpr(Mc, Mb), Ma);

        fmt::println("Ma = {}", Ma);
        fmt::println("Mb = {}", Mb);
        fmt::println("Mc = {}", Mc);
        fmt::println("Mres = rgpr(rgpr(Mc, Mb), Ma) = {}", M_res);
        fmt::println("");

        auto back_and_forth = rgpr(Mc, Ma);
        CHECK(gr1(back_and_forth) == vec2dp{});
        CHECK(gr3(back_and_forth) == pscalar2dp{1.0});
        CHECK(M_res == M1);
    }

    TEST_CASE("G<2,0,1>: sqrt(motor) function")
    {
        fmt::println("G<2,0,1>: sqrt(motor) function");

        // how to transfer line l1 into line l2, such that l2 = M ⟇ l1 ⟇ rrev(M)?
        //
        // We know R = l2 ⟇ l1 represents two consequtive reflections across l1 and l2.
        // So we need "half of the transformation" M = R^{1/2} = sqrt(R)

        // case a) rotation as defined by consequtive mirroring across intersecting lines
        fmt::println("\nG<2,0,1>: case a) rotation:\n");

        auto l1 = unitize(wdg(vec2dp{1, 1, 1}, vec2dp{2, 4, 1}));
        auto l2 = unitize(wdg(vec2dp{2, 1, 1}, vec2dp{4, 3, 1}));

        auto M = get_motor_from_lines(l1, l2);

        auto P0 = vec2dp{0, 7, 1};

        auto phi = angle(l1, l2);
        auto P_fix = unitize(rwdg(l1, l2));

        auto R = sqrt(M);

        fmt::println("l1            = {}", l1);
        fmt::println("l1 ⟇ l1       = {}", rgpr(l1, l1));
        fmt::println("l1 ⟇ rrev(l1) = {}", rgpr(l1, rrev(l1)));
        fmt::println("l2            = {}", l2);
        fmt::println("l2 ⟇ l2       = {}", rgpr(l2, l2));
        fmt::println("l2 ⟇ rrev(l2) = {}", rgpr(l2, rrev(l2)));
        fmt::println("");
        fmt::println("P_fix = {}, phi = {}°", P_fix, rad2deg(phi));
        fmt::println("P_fix ⟇ P_fix = {}", rgpr(P_fix, P_fix));
        fmt::println("");
        fmt::println("M(l1,l2)    = {}", M);
        fmt::println("M ⟇ rrev(M) = {}", rgpr(M, rrev(M)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = M ⟇ P0 ⟇ rrev(M) = {}", move2dp(P0, M));
        fmt::println("");
        fmt::println("R = sqrt(M) = {}", R);
        fmt::println("R ⟇ rrev(R) = {}", rgpr(R, rrev(R)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = R ⟇ P0 ⟇ rrev(R) = {}", move2dp(P0, R));
        fmt::println("");
        fmt::println("l1               = {}", l1);
        fmt::println("att(l1)          = {}", att(l1));
        fmt::println("cmpl(l1)         = {}", cmpl(l1));
        fmt::println("-weight_dual(l1) = {}", -weight_dual(l1));
        fmt::println("");
        fmt::println("gpr(l1,I_2dp)   = {}", l1 * I_2dp);
        fmt::println("rgpr(l1,I_2dp)  = {}", rgpr(l1, I_2dp));
        fmt::println("");

        // case b) translation as defined by parallel lines
        fmt::println("\nG<2,0,1>: case b) translation:\n");

        l1 = unitize(wdg(vec2dp{1, 0, 1}, vec2dp{1, 1, 1}));
        l2 = unitize(wdg(vec2dp{3, 0, 1}, vec2dp{3, 1, 1}));

        M = get_motor_from_lines(l1, l2);

        P0 = vec2dp{0, 7, 1};

        R = sqrt(M);

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
        fmt::println("P  = M ⟇ P0 ⟇ rrev(M) = {}", move2dp(P0, M));
        fmt::println("");
        fmt::println("R = sqrt(M) = {}", R);
        fmt::println("R ⟇ rrev(R) = {}", rgpr(R, rrev(R)));
        fmt::println("");
        fmt::println("P0                    = {}", P0);
        fmt::println("P  = R ⟇ P0 ⟇ rrev(R) = {}", move2dp(P0, R));
        fmt::println("");
        fmt::println("l1               = {}", l1);
        fmt::println("att(l1)          = {}", att(l1));
        fmt::println("cmpl(l1)         = {}", cmpl(l1));
        fmt::println("-weight_dual(l1) = {}", -weight_dual(l1));
        fmt::println("");
        fmt::println("gpr(l1,I_2dp)   = {}", l1 * I_2dp);
        fmt::println("rgpr(l1,I_2dp)  = {}", rgpr(l1, I_2dp));
        fmt::println("");
    }

} // PGA 2DP Tests