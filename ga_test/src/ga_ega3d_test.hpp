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
#include "ga/ga_ega.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using std::abs;              // bring std::abs into scope for ADL with our abs overload


/////////////////////////////////////////////////////////////////////////////////////////
// Euclideon geometric algebra of three-dimensional space: ega3d
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA 3D Tests")
{

    TEST_CASE("G<3,0,0>: ega3d")
    {
        fmt::println("");
        fmt::println("G<3,0,0>: ega3d");
        // 3d euklidean geometric algebra
        const algebra<3> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("   ega3d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   ega3d: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: default init")
    {
        fmt::println("Vec3d: default init");
        vec3d v;
        // fmt::println("   v = {}", v);
        CHECK(abs(v.x) < eps);
        CHECK(abs(v.y) < eps);
        CHECK(abs(v.z) < eps);
    }
    TEST_CASE("Vec3d: with curly braced intializer")
    {
        fmt::println("Vec3d: with curly braced intializer");
        vec3d v{0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(abs(v.x) < eps);
        CHECK(abs(v.y) < eps);
        CHECK(abs(v.z) < eps);
    }
    TEST_CASE("Vec3d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec3d: cp ctor & cp assign incl. type deduction");
        vec3d v1{1.0, 2.0, 3.0}; // init with double (type deduction)
        vec3d v2{v1};            // cp ctor
        vec3d v3 = v2;           // cp assign
        vec3d v4 = -v2;          // cp assign with unary minus

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


    TEST_CASE("Vec3d: fmt & cout printing")
    {
        fmt::println("Vec3d: fmt & cout printing:\n");

        vec3d pf{1.0f, 2.00001f, 3.0f};
        vec3d pd{1.0, 2.00001, 3.0};

        // std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.4f}", pf);
        fmt::println("");

        // std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.4f}", pd);
        fmt::println("");

        std::vector<Vec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:g}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("G<3,0,0>: Scalar3d and PScalar3d formatting tests")
    {
        fmt::println("G<3,0,0>: Scalar3d and PScalar3d formatting tests");

        // Test Scalar3d formatting
        Scalar3d<double> scalar_val{1.41421};

        // Basic output
        std::string scalar_basic = fmt::format("{}", scalar_val);
        CHECK(scalar_basic == "Scalar3d(1.41421)");

        // Two decimal places
        std::string scalar_two_decimals = fmt::format("{:.2f}", scalar_val);
        CHECK(scalar_two_decimals == "Scalar3d(1.41)");

        // Scientific notation
        std::string scalar_scientific = fmt::format("{:.2e}", scalar_val);
        CHECK(scalar_scientific == "Scalar3d(1.41e+00)");

        // Contextual usage
        std::string scalar_contextual = fmt::format("Magnitude: {:.3f}", scalar_val);
        CHECK(scalar_contextual == "Magnitude: Scalar3d(1.414)");

        // Test PScalar3d formatting
        PScalar3d<double> pscalar_val{1.73205};

        // Basic output
        std::string pscalar_basic = fmt::format("{}", pscalar_val);
        CHECK(pscalar_basic == "PScalar3d(1.73205)");

        // Two decimal places
        std::string pscalar_two_decimals = fmt::format("{:.2f}", pscalar_val);
        CHECK(pscalar_two_decimals == "PScalar3d(1.73)");

        // Scientific notation
        std::string pscalar_scientific = fmt::format("{:.2e}", pscalar_val);
        CHECK(pscalar_scientific == "PScalar3d(1.73e+00)");

        // Contextual usage
        std::string pscalar_contextual = fmt::format("Volume: {:.3f}", pscalar_val);
        CHECK(pscalar_contextual == "Volume: PScalar3d(1.732)");

        fmt::println("   Scalar3d basic: {}", scalar_basic);
        fmt::println("   Scalar3d 2-decimal: {}", scalar_two_decimals);
        fmt::println("   Scalar3d scientific: {}", scalar_scientific);
        fmt::println("   Scalar3d contextual: {}", scalar_contextual);

        fmt::println("   PScalar3d basic: {}", pscalar_basic);
        fmt::println("   PScalar3d 2-decimal: {}", pscalar_two_decimals);
        fmt::println("   PScalar3d scientific: {}", pscalar_scientific);
        fmt::println("   PScalar3d contextual: {}", pscalar_contextual);
    }

    TEST_CASE("Vec3d: comprehensive format specifier tests with pre-computed values")
    {
        fmt::println(
            "Vec3d: comprehensive format specifier tests with pre-computed values:\n");

        // Test Vec3d and BiVec3d with precise values for predictable output
        vec3d vd{3.14159, 2.71828, 1.41421};
        Vec3d<float> vf{3.14159f, 2.71828f, 1.41421f};
        bivec3d bv{1.2, 2.3, 3.4};

        // Test default formatting and verify expected output
        std::string default_vd = fmt::format("{}", vd);
        std::string default_vf = fmt::format("{}", vf);
        std::string default_bv = fmt::format("{}", bv);
        fmt::println("   Default Vec3d:   '{}'", default_vd);
        fmt::println("   Default Vec3d<f>: '{}'", default_vf);
        fmt::println("   Default BiVec3d: '{}'", default_bv);

        // Expected patterns for Vec3d
        CHECK(default_vd.find("Vec3d(") == 0);
        CHECK(default_vd.find(",") != std::string::npos);
        CHECK(default_vd.back() == ')');

        // Expected patterns for BiVec3d
        CHECK(default_bv.find("BiVec3d(") == 0);
        CHECK(default_bv.back() == ')');

        // Test precision formatting with expected outputs for Vec3d
        std::string two_decimals = fmt::format("{:.2f}", vd);
        std::string expected_2f = "Vec3d(3.14,2.72,1.41)";
        fmt::println("   Vec3d .2f:       '{}' (expected: '{}')", two_decimals,
                     expected_2f);
        CHECK(two_decimals == expected_2f);

        std::string three_decimals = fmt::format("{:.3f}", vd);
        std::string expected_3f = "Vec3d(3.142,2.718,1.414)";
        fmt::println("   Vec3d .3f:       '{}' (expected: '{}')", three_decimals,
                     expected_3f);
        CHECK(three_decimals == expected_3f);

        // Test BiVec3d formatting
        std::string bivec_2f = fmt::format("{:.2f}", bv);
        std::string expected_bv_2f = "BiVec3d(1.20,2.30,3.40)";
        fmt::println("   BiVec3d .2f:     '{}' (expected: '{}')", bivec_2f,
                     expected_bv_2f);
        CHECK(bivec_2f == expected_bv_2f);

        // Test scientific notation
        std::string scientific = fmt::format("{:.2e}", vd);
        fmt::println("   Scientific .2e:  '{}'", scientific);
        CHECK(scientific.find("Vec3d(") == 0);
        bool has_scientific = (scientific.find("e+") != std::string::npos ||
                               scientific.find("e-") != std::string::npos);
        CHECK(has_scientific);

        // Test with float type
        std::string float_2f = fmt::format("{:.2f}", vf);
        std::string expected_float_2f = "Vec3d(3.14,2.72,1.41)";
        fmt::println("   Float .2f:       '{}' (expected: '{}')", float_2f,
                     expected_float_2f);
        CHECK(float_2f == expected_float_2f);

        // Test contextual usage like in the reference example
        std::string position_3d = fmt::format("Position: {:.2f}", vd);
        std::string expected_pos_3d = "Position: Vec3d(3.14,2.72,1.41)";
        fmt::println("   Context Vec3d:   '{}' (expected: '{}')", position_3d,
                     expected_pos_3d);
        CHECK(position_3d == expected_pos_3d);

        std::string orientation_log = fmt::format("Orientation: {:.2f}", bv);
        std::string expected_orient = "Orientation: BiVec3d(1.20,2.30,3.40)";
        fmt::println("   Context BiVec:   '{}' (expected: '{}')", orientation_log,
                     expected_orient);
        CHECK(orientation_log == expected_orient);

        fmt::println("");
        fmt::println(
            "   All Vec3d and BiVec3d format tests passed with expected values!");
    }

    TEST_CASE("Vec3d: comparison float")
    {
        fmt::println("Vec3d: comparison float");

        Vec3d<float> v1f{1.0f, 2.0f, 3.0f};
        Vec3d<float> v2f{2.0f, 4.0f, 3.0f};
        Vec3d<float> v3f{1.0f, 2.0000001f, 3.0f};
        Vec3d<float> v4f{v1f};

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

    TEST_CASE("Vec3d: comparison double")
    {
        fmt::println("Vec3d: comparison double");

        Vec3d<double> v1d{1.0, 2.0, 3.0};
        Vec3d<double> v2d{2.0, 4.0, 3.0};
        Vec3d<double> v3d{1.0, 2.0000000000000001, 3.0};
        Vec3d<double> v4d{v1d};

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

    TEST_CASE("Vec3d: vector space and linearity tests")
    {
        fmt::println("Vec3d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        vec3d p0;
        vec3d p1{1.0, 2.0, 3.0};
        vec3d p2{2.0, 4.0, 6.0};
        vec3d p3{3.0, 6.0, 9.0};
        vec3d p4 = -p1; // assignment using unary minus
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

    TEST_CASE("Vec3d: inner product properties")
    {
        fmt::println("Vec3d: inner product properties");

        double a = 2.35;
        vec3d u{1.0, 2.0, 1.0};
        vec3d v{-0.5, 3.0, 0.5};
        vec3d w{3.0, 6.0, -3.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: operations - norm, inverse, dot")
    {
        fmt::println("Vec3d: operations - norm, inverse, dot");

        scalar3d s1{3.2};

        vec3d v1{2.0, 1.0, 2.0};
        vec3d v2{normalize(v1)};

        vec3d v3{2.0, 6.0, -4.0};
        vec3d v4{inv(v3)};

        bivec3d b1{-2.3, 1.2, 4.5};
        pscalar3d ps1{-4.7};
        mvec3d_e mve1{s1, b1};
        mvec3d_u mvu1{v1, ps1};
        mvec3d mv1{s1, v1, b1, ps1};

        // fmt::println("v1 = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2 = normalize(v1) = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("v3 = {: .4f}, nrm(v1) = {: .4f}", v3, nrm(v3));
        // fmt::println(
        //     "v4 = inv(v3) = {: .4f}, nrm(v3) = {: .4f}, nrm(v3)*nrm(v4) = {: .4f}", v4,
        //     nrm(v4), nrm(v3) * nrm(v4));

        CHECK(abs(nrm_sq(v1) - 9.0) < eps);
        CHECK(abs(nrm_sq(v2) - 1.0) < eps);
        CHECK(abs(value_t(dot(v4, v3)) - 1.0) < eps);

        // check inverses - scalar
        // fmt::println("");
        // fmt::println("s1 * inv(s1) = {}", s1 * inv(s1)); // s
        CHECK(abs(nrm(s1 * inv(s1)) - 1) < eps);
        CHECK(abs(value_t(inv(s1) - rev(s1) / nrm_sq(s1))) < eps);

        // check inverses - vector
        // fmt::println("v1 * inv(v1) = {}", v1 * inv(v1)); // mv_e
        CHECK(abs(nrm(gr0(v1 * inv(v1))) - 1) < eps);
        CHECK(abs(nrm(gr2(v1 * inv(v1))) - 0) < eps);
        CHECK(abs(nrm(inv(v1) - rev(v1) / nrm_sq(v1))) < eps);

        // check inverses - bivector
        // fmt::println("b1 * inv(b1) = {}", b1 * inv(b1)); // mv_e
        CHECK(abs(nrm(gr0(b1 * inv(b1))) - 1) < eps);
        CHECK(abs(nrm(gr2(b1 * inv(b1))) - 0) < eps);
        CHECK(abs(nrm(inv(b1) - rev(b1) / nrm_sq(b1))) < eps);

        // check inverses - pseudoscalar
        // fmt::println("ps1 * inv(ps1) = {}", ps1 * inv(ps1)); // s
        CHECK(abs(nrm(ps1 * inv(ps1)) - 1) < eps);
        CHECK(abs(nrm(inv(ps1) - rev(ps1) / nrm_sq(ps1))) < eps);

        // check inverses - even-grade multivector
        // fmt::println("mve1 * inv(mve1) = {}", mve1 * inv(mve1)); // mv_e
        CHECK(abs(nrm(gr0(mve1 * inv(mve1))) - 1) < eps);
        CHECK(abs(nrm(gr2(mve1 * inv(mve1))) - 0) < eps);
        CHECK(abs(nrm(inv(mve1) - rev(mve1) / nrm_sq(mve1))) < eps);

        // check inverses - odd-grade multivector
        // fmt::println("mvu1 * inv(mvu1) = {}", mvu1 * inv(mvu1)); // mv_e
        CHECK(abs(nrm(gr0(mvu1 * inv(mvu1))) - 1) < eps);
        CHECK(abs(nrm(gr2(mvu1 * inv(mvu1))) - 0) < eps);
        CHECK(abs(nrm(inv(mvu1) - rev(mvu1) / nrm_sq(mvu1))) < eps);

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

    TEST_CASE("Vec3d: operations - angle I")
    {
        fmt::println("Vec3d: operations - angle");

        vec3d v1{1.0, 0.0, 0.0};
        vec3d v2{normalize(vec3d(1.0, 1.0, 0.0))};
        vec3d v3{0.0, 1.0, 0.0};
        vec3d v4{normalize(vec3d(-1.0, 1.0, 0.0))};
        vec3d v5{-1.0, 0.0, 0.0};
        vec3d v6{normalize(vec3d(-1.0, -1.0, 0.0))};
        vec3d v7{0.0, -1.0, 0.0};
        vec3d v8{normalize(vec3d(1.0, -1.0, 0.0))};

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
        CHECK(v6 == normalize(vec3d(-1.0, -1.0, 0.0)));
        CHECK(v7 == normalize(vec3d(0.0, -1.0, 0.0)));
        CHECK(v8 == normalize(vec3d(1.0, -1.0, 0.0)));
    }

    TEST_CASE("Vec3d: operations - angle II")
    {
        fmt::println("Vec3d: operations - angle II");

        std::vector<std::tuple<double, Vec3d<double>>> v1;
        std::vector<std::tuple<double, Vec3d<double>>> v2;
        std::vector<std::tuple<double, Vec3d<double>>> v3;

        // only positive angles are easy to implement vs. the 2d case

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi), std::sin(phi), 0.0);
            v1.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_3d, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2), 0.0);
            v2.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e2_3d, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4), 0.0);
            v3.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_3d + e2_3d, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(abs(phi - angle(e1_3d, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(abs(phi - angle(e2_3d, c)) < eps);
        }
        auto ref_vec = normalize(e1_3d + e2_3d);
        for (auto const& [phi, c] : v3) {
            CHECK(abs(phi - angle(ref_vec, c)) < eps);
        }
    }

    TEST_CASE("Vec3d: operations - wedge")
    {
        fmt::println("Vec3d: operations - wedge");

        vec3d v1{1.0, 0.0, 0.0};
        vec3d v2{normalize(vec3d(1.0, 1.0, 0.0))};
        vec3d v3{0.0, 1.0, 0.0};
        vec3d v4{normalize(vec3d(-1.0, 1.0, 0.0))};
        vec3d v5{-1.0, 0.0, 0.0};
        vec3d v6{normalize(vec3d(-1.0, -1.0, 0.0))};
        vec3d v7{0.0, -1.0, 0.0};
        vec3d v8{normalize(vec3d(1.0, -1.0, 0.0))};

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

        CHECK(abs(nrm(wdg(v1, v1)) - sin(angle(v1, v1))) < eps);
        CHECK(abs(nrm(wdg(v1, v2)) - sin(angle(v1, v2))) < eps);
        CHECK(abs(nrm(wdg(v1, v3)) - sin(angle(v1, v3))) < eps);
        CHECK(abs(nrm(wdg(v1, v4)) - sin(angle(v1, v4))) < eps);
        CHECK(abs(nrm(wdg(v1, v5)) - sin(angle(v1, v5))) < eps);
        CHECK(abs(nrm(wdg(v1, v6)) - sin(angle(v1, v6))) < eps);
        CHECK(abs(nrm(wdg(v1, v7)) - sin(angle(v1, v7))) < eps);
        CHECK(abs(nrm(wdg(v1, v8)) - sin(angle(v1, v8))) < eps);
    }

    TEST_CASE("Vec3d: operations - project / reject / reflect (vector - vector)")
    {
        fmt::println("Vec3d: operations - project / reject / reflect (vector - vector)");

        vec3d v1{5.0, 1.0, 1.0};
        vec3d v2{2.0, 2.0, 1.0};

        vec3d v3{project_onto(v1, v2)};
        vec3d v4{reject_from(v1, v2)};
        vec3d v5{v3 + v4};


        // fmt::println("");
        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f}, nrm(v3) = {: .4f}", v3,
        //              nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, nrm(v4) = {: .4f}", v4,
        //              nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f}, nrm(v5) = {: .4f}", v5,
        //              nrm(v5));
        // fmt::println("");
        // // this helps to understand, why the dot-product is sufficient
        // auto w = wdg(v1, v2);
        // auto i = inv(v2);
        // fmt::println("wdg(v1,v2)         = {: .4f}", w);
        // fmt::println("inv(v2)            = {: .4f}", i);
        // fmt::println("wdg(v1,v2)*inv(v2) = {: .4f}", w * i);
        // fmt::println("");

        vec3d v{4.0, 1.0, 1.0};
        vec3d b{e2_3d};
        auto B = bivec3d{e12_3d};

        // auto UB = bivec3d{e23_3d + e12_3d};
        // fmt::println("v   = {}", v);
        // fmt::println("b   = {}", b);
        // fmt::println("B   = {}", B);
        // fmt::println("UB  = {}", UB);
        // fmt::println("reflect_on_vec(v,b)  = {}", reflect_on_vec(v, b));
        // fmt::println("reflect_on(v,e3_3d)  = {}", reflect_on(v, e3_3d));
        // fmt::println("reflect_on(v,B)  = {}", reflect_on(v, B));
        // fmt::println("reflect_on(UB,B) = {}", reflect_on(UB, B));
        // fmt::println("");

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);
        CHECK(std::abs(to_val(dot(v3, v4))) <= eps);

        // just to suppress unused variable warnings
        CHECK(b == e2_3d);

        CHECK(reflect_on(v, e3_3d) == reflect_on(v, B));

        // checking time required
        //
        // auto start = std::chrono::system_clock::now();
        // for (size_t i = 0; i < 10000000; ++i) {
        //     vec3d v = reject_from(v1, v2);
        // }
        // auto end = std::chrono::system_clock::now();
        // auto elapsed =
        // std::chrono::duration_cast<std::chrono::milliseconds>(end
        // -
        // start); fmt::println("The measurement took {}", elapsed);

        // point reflected on a vector
        vec3d p{4, 1, 0};
        CHECK(reflect_on_vec(p, x_dir_3d) == vec3d{4, -1, 0});

        // coordinate axis reflected on perpendicular axis yield their negatives
        CHECK(reflect_on_vec(y_dir_3d, x_dir_3d) == -y_dir_3d);
        CHECK(reflect_on_vec(z_dir_3d, x_dir_3d) == -z_dir_3d);
        CHECK(reflect_on_vec(x_dir_3d, y_dir_3d) == -x_dir_3d);
        CHECK(reflect_on_vec(z_dir_3d, y_dir_3d) == -z_dir_3d);
        CHECK(reflect_on_vec(x_dir_3d, z_dir_3d) == -x_dir_3d);
        CHECK(reflect_on_vec(y_dir_3d, z_dir_3d) == -y_dir_3d);

        // coordinate axis reflected on itself remains itself (identity)
        CHECK(reflect_on_vec(x_dir_3d, x_dir_3d) == x_dir_3d);
        CHECK(reflect_on_vec(y_dir_3d, y_dir_3d) == y_dir_3d);
        CHECK(reflect_on_vec(z_dir_3d, z_dir_3d) == z_dir_3d);

        // point reflected on a hyperplane that the vector is a normal to
        // the hyperplane can be created by taking the dual (or the rcmpl) of the normal
        CHECK(reflect_on(p, dual(y_dir_3d)) == vec3d{4, -1, 0});
        // alternatively the plane can be used directly (represented by a bivector)
        CHECK(reflect_on(p, zx_3d) == vec3d{4, -1, 0});

        // coordinate axis reflected on perpendicular base planes yield their negatives
        CHECK(reflect_on(x_dir_3d, yz_3d) == -x_dir_3d);
        CHECK(reflect_on(x_dir_3d, dual(x_dir_3d)) == -x_dir_3d);
        CHECK(reflect_on(y_dir_3d, zx_3d) == -y_dir_3d);
        CHECK(reflect_on(y_dir_3d, dual(y_dir_3d)) == -y_dir_3d);
        CHECK(reflect_on(z_dir_3d, xy_3d) == -z_dir_3d);
        CHECK(reflect_on(z_dir_3d, dual(z_dir_3d)) == -z_dir_3d);

        // a coordinate plane reflected on itself remains itself (identity)
        CHECK(reflect_on(yz_3d, yz_3d) == yz_3d);
        CHECK(reflect_on(zx_3d, zx_3d) == zx_3d);
        CHECK(reflect_on(xy_3d, xy_3d) == xy_3d);

        // reflect planes on planes directly
        CHECK(reflect_on(e23_3d + e12_3d, e12_3d) == -e23_3d + e12_3d);
    }

    TEST_CASE("Vec3d: operations - project / reject / reflect (vector - bivector)")
    {
        fmt::println(
            "Vec3d: operations - project / reject / reflect (vector - bivector)");

        vec3d v1{5.0, 3.0, 1.0};
        bivec3d v2 = wdg(vec3d{0.0, 0.0, 2.0}, vec3d{2.0, 0.0, 0.0});

        vec3d v3{project_onto(v1, v2)};
        vec3d v4{reject_from(v1, v2)};
        vec3d v5{v3 + v4};

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f}, nrm(v3) = {: .4f}", v3,
        //              nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, nrm(v4) = {: .4f}", v4,
        //              nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f}, nrm(v5) = {: .4f}", v5,
        //              nrm(v5));

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);
        CHECK(project_onto(v1, e31_3d) == vec3d{5, 0, 1});
        CHECK(reject_from(v1, e31_3d) == vec3d{0, 3, 0});

        // from Macdonald, "Linear and Geometric Algebra", Exercise 7.14, p. 129
        CHECK(reflect_on(wdg(e1_3d + e3_3d, e2_3d), e12_3d) == wdg(e1_3d - e3_3d, e2_3d));
    }

    TEST_CASE("Vec3d: operations - rotations")
    {
        fmt::println("Vec3d: operations - rotations");

        // fmt::println("");

        double phi = deg2rad(30.0);
        auto c = vec3d(std::cos(phi), std::sin(phi), 0.0);
        auto d = rotate(e1_3d, get_rotor(e12_3d, phi));
        auto cm = vec3d(std::cos(-phi), std::sin(-phi), 0.0);
        auto dm = rotate(e1_3d, get_rotor(e12_3d, -phi));

        // fmt::println("c = {}", c);
        // fmt::println("d = {}", d);
        // fmt::println("cm = {}", cm);
        // fmt::println("dm = {}", dm);

        CHECK(c == d);
        CHECK(cm == dm);

        // rotations in a plane (using rotors)
        std::vector<std::tuple<double, Vec3d<double>>> v;

        // fmt::println("");
        for (int i = -12; i <= 12; ++i) {
            double angle = i * pi / 12;
            auto c_val = vec3d(std::cos(angle), std::sin(angle), 0.0);
            auto b = cmpl(c_val); // bivec, representing plane with normal cos_val
            v.emplace_back(std::make_tuple(angle, c_val));
            // fmt::println("   i={: 3}: angle={: .4f}, angle={: 4.0f}°, c_val={: .3f},"
            //              " angle={: .4f}",
            //              i, angle, rad2deg(angle), c_val, angle(e1_3d, c_val));
            CHECK(c_val == rotate(e1_3d, get_rotor(e12_3d, angle)));

            // check rotate optimization
            CHECK(rotate(e1_3d, get_rotor(e12_3d, angle)) ==
                  rotate_opt(e1_3d, get_rotor(e12_3d, angle)));

            // rotation of a ega3d multivector is to rotate the vector and bivector parts
            // exclusively
            auto mv_rev = mvec3d{scalar3d{4.0}, e1_3d, e23_3d, pscalar3d{-2.0}};
            auto mv_rot = mvec3d{scalar3d{4.0}, c_val, b, pscalar3d{-2.0}};
            CHECK(mv_rot == rotate(mv_rev, get_rotor(e12_3d, angle)));
        }

        // check same rotation applied to many points
        std::vector<Vec3d<double>> vec_ref(100);
        std::vector<Vec3d<double>> vec_rot(100);
        std::vector<Vec3d<double>> vec_rot_calc(100);

        // prepare for checking vector-based rotation (after the loop)
        phi = pi / 12;
        for (size_t i = 0; i < 100; ++i) {
            vec_ref.push_back(e1_3d);
            vec_rot.emplace_back(rotate(e1_3d, get_rotor(e12_3d, phi)));
        }

        vec_rot_calc = rotate_opt(vec_ref, get_rotor(e12_3d, phi));
        for (size_t i = 0; i < 100; ++i) {
            CHECK(vec_rot_calc[i] == vec_rot[i]);
        }

        // fmt::println("");
    }

    TEST_CASE("Vec3d: operations - simple rotation")
    {
        fmt::println("Vec3d: operations - simple rotation");

        // rotation in e31 plane, normalized rotation axis is showing in e2 direction
        auto const axis = e2_3d;
        auto const phi = deg2rad(90);
        auto const X0 = vec3d{0, 0, 2};

        auto const M_rot = get_rotor(cmpl(axis), phi);
        auto const X_rot = rotate(X0, M_rot);

        fmt::println("axis = {}", axis);
        fmt::println("phi  = {}", phi);
        fmt::println("X0   = {}", X0);
        fmt::println("");
        fmt::println("M_rot    = {}", M_rot);
        fmt::println("X_rot    = {}", X_rot);
        fmt::println("");
        auto const Xh = 0.5 * (X_rot + X0);
        fmt::println("Xh = 0.5*(X + X0) = {}", Xh);
        fmt::println("X0*Xh = {}", normalize(X0 * Xh));
        fmt::println("Xh*X0 = {}", normalize(Xh * X0));

        CHECK(X_rot == vec3d{2.0, 0.0, 0.0});

        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> Gram-Schmidt-Orthogonalization
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: Gram-Schmidt-Orthogonalization")
    {
        fmt::println("Vec3d: Gram-Schmidt-Orthogonalization - part 1 (plane)");

        vec3d v1{2.0, 2.0, 0.0};
        vec3d v2{-1.5, 2.5, 1.0};
        vec3d v3{1.5, -2.5, -1.0};

        std::vector<vec3d> og12 = gs_orthogonal(v1, v2);
        std::vector<vec3d> og13 = gs_orthogonal(v1, v3);

        std::vector<vec3d> on12 = gs_orthonormal(v1, v2);
        std::vector<vec3d> on13 = gs_orthonormal(v1, v3);

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}, angle = {: .2f}°", v1,
        // nrm(v1),
        //              rad2deg(angle(e1_3d, v1)));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v2,
        // nrm(v2),
        //              rad2deg(angle(e1_3d, v2)));
        // fmt::println("v3  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v3,
        // nrm(v3),
        //              rad2deg(angle(e1_3d, v3)));
        // fmt::println("");
        // fmt::println("og12[0]  = {: .4f}, nrm(og12[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              og12[0], nrm(og12[0]), rad2deg(angle(e1_3d, og12[0])));
        // fmt::println("og12[1]  = {: .4f}, nrm(og12[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              og12[1], nrm(og12[1]), rad2deg(angle(e1_3d, og12[1])));
        // fmt::println("on12[0]  = {: .4f}, nrm(on12[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              on12[0], nrm(on12[0]), rad2deg(angle(e1_3d, on12[0])));
        // fmt::println("on12[1]  = {: .4f}, nrm(on12[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              on12[1], nrm(on12[1]), rad2deg(angle(e1_3d, on12[1])));
        // fmt::println("");
        // fmt::println("og13[0]  = {: .4f}, nrm(og13[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              og13[0], nrm(og13[0]), rad2deg(angle(e1_3d, og13[0])));
        // fmt::println("og13[1]  = {: .4f}, nrm(og13[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              og13[1], nrm(og13[1]), rad2deg(angle(e1_3d, og13[1])));
        // fmt::println("on13[0]  = {: .4f}, nrm(on13[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              on13[0], nrm(on13[0]), rad2deg(angle(e1_3d, on13[0])));
        // fmt::println("on13[1]  = {: .4f}, nrm(on13[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              on13[1], nrm(on13[1]), rad2deg(angle(e1_3d, on13[1])));
        // fmt::println("");

        CHECK(abs(value_t(dot(og12[0], og12[1]))) < eps);
        CHECK(abs(value_t(dot(og13[0], og13[1]))) < eps);

        CHECK(abs(value_t(dot(on12[0], on12[1]))) < eps);
        CHECK(abs(nrm(on12[0]) - 1.0) < eps);
        CHECK(abs(nrm(on12[1]) - 1.0) < eps);

        CHECK(abs(value_t(dot(on13[0], on13[1]))) < eps);
        CHECK(abs(nrm(on13[0]) - 1.0) < eps);
        CHECK(abs(nrm(on13[1]) - 1.0) < eps);

        fmt::println("Vec3d: Gram-Schmidt-Orthogonalization - part 2 (space)");

        v1 = vec3d{2.0, 2.0, 0.0};
        v2 = vec3d{-1.5, 2.5, 1.0};
        v3 = vec3d{1.5, -2.5, -6.0}; // all three vectors must be linear independent

        std::vector<vec3d> og123 = gs_orthogonal(v1, v2, v3);
        std::vector<vec3d> on123 = gs_orthonormal(v1, v2, v3);

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}, angle = {: .2f}°", v1,
        // nrm(v1),
        //              rad2deg(angle(e1_3d, v1)));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v2,
        // nrm(v2),
        //              rad2deg(angle(e1_3d, v2)));
        // fmt::println("v3  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v3,
        // nrm(v3),
        //              rad2deg(angle(e1_3d, v3)));
        // fmt::println("");
        // fmt::println("og123[0]  = {: .4f}, nrm(og123[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              og123[0], nrm(og123[0]), rad2deg(angle(e1_3d, og123[0])));
        // fmt::println("og123[1]  = {: .4f}, nrm(og123[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              og123[1], nrm(og123[1]), rad2deg(angle(e1_3d, og123[1])));
        // fmt::println("og123[2]  = {: .4f}, nrm(og123[2]) = {: .4f}, angle = {:
        // .2f}°",
        //              og123[2], nrm(og123[2]), rad2deg(angle(e1_3d, og123[2])));
        // fmt::println("");
        // fmt::println("on123[0]  = {: .4f}, nrm(on123[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              on123[0], nrm(on123[0]), rad2deg(angle(e1_3d, on123[0])));
        // fmt::println("on123[1]  = {: .4f}, nrm(on123[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              on123[1], nrm(on123[1]), rad2deg(angle(e1_3d, on123[1])));
        // fmt::println("on123[2]  = {: .4f}, nrm(on123[2]) = {: .4f}, angle = {:
        // .2f}°",
        //              on123[2], nrm(on123[2]), rad2deg(angle(e1_3d, on123[2])));
        // fmt::println("");
        // fmt::println("sign(wdg(on123[0], wdg(on123[1], on123[2]))/I_3d) = {} "
        //              "(+1: right-handed, -1: left-handed)",
        //              sign(wdg(on123[0], wdg(on123[1], on123[2])) / I_3d));
        // fmt::println("");

        CHECK(abs(value_t(dot(og123[0], og123[1]))) < eps);
        CHECK(abs(value_t(dot(og123[1], og123[2]))) < eps);
        CHECK(abs(value_t(dot(og123[0], og123[2]))) < eps);

        CHECK(abs(value_t(dot(on123[0], on123[1]))) < eps);
        CHECK(abs(value_t(dot(on123[1], on123[2]))) < eps);
        CHECK(abs(value_t(dot(on123[0], on123[2]))) < eps);

        CHECK(abs(nrm(on123[0]) - 1.0) < eps);
        CHECK(abs(nrm(on123[1]) - 1.0) < eps);
        CHECK(abs(nrm(on123[2]) - 1.0) < eps);
    }

    TEST_CASE("Vec3d: cross-product")
    {
        fmt::println("Vec3d: cross-product");

        vec3d u{1.0, 1.0, 0.0};
        vec3d v{0.0, 1.0, 1.0};
        vec3d w{1.0, 1.0, 1.0};

        vec3d u_cross_v = cross(u, v);
        bivec3d u_wdg_v = wdg(u, v);

        CHECK(u_cross_v == dual(u_wdg_v));
        CHECK(u_wdg_v == dual(u_cross_v));

        // definitions using the pseudoscalars directly are valid independent of the
        // duality definitions

        CHECK(u_cross_v == u_wdg_v * rev(I_3d));
        CHECK(u_cross_v == rev(u_wdg_v) * I_3d);
        CHECK(u_cross_v == -I_3d * u_wdg_v);
        CHECK(u_cross_v == -u_wdg_v * I_3d);

        CHECK(wdg(u, v) == u_cross_v * I_3d);
        CHECK(wdg(u, v) == I_3d * u_cross_v);
        CHECK(wdg(u, v) == -u_cross_v * rev(I_3d));

        // double cross product identity
        CHECK(cross(u, cross(v, w)) == (u << wdg(v, w)));
        CHECK(cross(u, cross(v, w)) == (dot(u, w) * v - dot(u, v) * w));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d: default init")
    {
        fmt::println("MVec3d: default init");
        // default initialization
        mvec3d v;
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

    TEST_CASE("MVec3d: with curly braced intializer")
    {
        fmt::println("MVec3d: with curly braced intializer");
        // default initialization
        mvec3d v{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
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

    TEST_CASE("MVec3d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec3d: cp ctor & cp assign incl. type deduction");
        // default initialization
        mvec3d v1{0.0,  1.0,  2.0,  3.0,
                  23.0, 31.0, 12.0, 123.0}; // init with double (type deduction)
        mvec3d v2{v1};                      // cp ctor
        mvec3d v3 = v2;                     // cp assign
        mvec3d v4 = -v3;                    // cp assign with unary minus

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

    TEST_CASE("MVec3d: defining basic types and ctor checks")
    {
        fmt::println("MVec3d: defining basic types and ctor checks");

        auto s = scalar3d(5.0);
        auto v = vec3d{1.0, 2.0, 1.0};
        auto B = bivec3d{-1.0, 2.0, 1.0};
        auto ps = pscalar3d(-5.0);

        auto mv0 = mvec3d{s};
        auto mv1 = mvec3d{v};
        auto mv2 = mvec3d{B};
        auto mv3 = mvec3d{ps};

        auto mv5a = mvec3d{s, B};
        auto mv5b = mvec3d_e{s};
        auto mv5c = mvec3d_e{B};
        auto mv5d = mvec3d_e{s, B};
        auto mv5e = mvec3d{mv5d};

        auto mv6a = mvec3d{v, ps};
        auto mv6b = mvec3d_u{v};
        auto mv6c = mvec3d_u{ps};
        auto mv6d = mvec3d_u{v, ps};
        auto mv6e = mvec3d{mv6d};

        auto mv7 = mvec3d{s, v, B, ps};

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
        CHECK(gr1(mv5a) == vec3d{});
        CHECK(gr2(mv5a) == B);
        CHECK(gr3(mv5a) == pscalar3d{});

        CHECK(gr0(mv5b) == s);
        CHECK(gr2(mv5b) == bivec3d{});

        CHECK(gr0(mv5c) == scalar3d{});
        CHECK(gr2(mv5c) == B);

        CHECK(gr0(mv5d) == s);
        CHECK(gr2(mv5d) == B);

        CHECK(gr0(mv5e) == s);
        CHECK(gr1(mv5e) == vec3d{});
        CHECK(gr2(mv5e) == B);
        CHECK(gr3(mv5e) == pscalar3d{});

        CHECK(gr0(mv6a) == scalar3d{});
        CHECK(gr1(mv6a) == v);
        CHECK(gr2(mv6a) == bivec3d{});
        CHECK(gr3(mv6a) == ps);

        CHECK(gr1(mv6b) == v);
        CHECK(gr3(mv6b) == pscalar3d{});

        CHECK(gr1(mv6c) == vec3d{});
        CHECK(gr3(mv6c) == ps);

        CHECK(gr1(mv6d) == v);
        CHECK(gr3(mv6d) == ps);

        CHECK(gr0(mv6e) == scalar3d{});
        CHECK(gr1(mv6e) == v);
        CHECK(gr2(mv6e) == bivec3d{});
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

    TEST_CASE("MVec3d: fmt & cout printing")
    {
        fmt::println("MVec3d: fmt & cout printing:\n");

        mvec3d pf{1.0f, 2.00001f, 0.0f, 3.0f, 1.0f, 2.00001f, 0.0f, 3.0f};
        mvec3d pd{1.0, 2.00001, 0.0, 3.0, 1.0, 2.00001, 0.0, 3.0};

        std::cout << "    cout: pf = " << pf << std::endl;
        fmt::println("    fmt:  pf = {}", pf);
        fmt::println("    fmt:  pf = {:.4f}", pf);
        fmt::println("");

        std::cout << "    cout: pd = " << pd << std::endl;
        fmt::println("    fmt:  pd = {}", pd);
        fmt::println("    fmt:  pd = {:.4f}", pd);
        fmt::println("");

        std::vector<MVec3d<double>> vp1{{1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0},
                                        {0.5, 1.5, 2.0, 2.5, 1.0, 1.0, 1.0, 2.0}};
        fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("    fmt: vp1 = {:g}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("MVec3d: vector space and linearity tests")
    {
        fmt::println("MVec3d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        mvec3d p0;
        mvec3d p1{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
        mvec3d p2{0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0};
        mvec3d p3{0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0};
        mvec3d p4 = -p1; // assignment using unary minus
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
    // MVec3d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d: wedge product - basic properties")
    {
        fmt::println("MVec3d: wedge product - basic properties");

        vec3d v1{1.0, 2.0, 1.0};
        vec3d v2{0.5, 3.0, 2.0};
        vec3d v3{-2.0, 6.0, 3.0};
        bivec3d B{-3.4, 2.1, 0.7};
        pscalar3d ps{1.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3d{sd};
        auto t = scalar3d{st};

        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v1) == bivec3d{});                      // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s)); // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, t) == s * t);
        CHECK(wdg(s, v1) == wdg(v1, s));       // wdg between scalar and vector
        CHECK(wdg(s, B) == wdg(B, s));         // wdg between scalar and bivector
        CHECK(wdg(scalar3d{1.0}, I_3d) == ps); // wdg between scalar and trivector
        CHECK(wdg(I_3d, scalar3d{1.0}) == ps); // (=pseudoscalar); 1.0 is neutral element
                                               // for the wedge product
        CHECK(wdg(s, v1) == sd * v1);          // wdg between scalar and vector
    }

    TEST_CASE("MVec3d: regressive wedge product - basic properties")
    {
        fmt::println("MVec3d: regressive wedge product - basic properties");

        bivec3d B1{1.0, 2.0, 1.0};
        bivec3d B2{0.5, 3.0, 2.0};
        bivec3d B3{-2.0, 6.0, 3.0};
        vec3d v1{-3.4, 2.1, 0.7};
        vec3d v2{1, 2, 3};
        vec3d z = {-2, 3, 7};
        scalar3d sc{2.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = pscalar3d{sd};
        auto t = pscalar3d{st};

        CHECK(rwdg(rwdg(B1, B2), B3) == rwdg(B1, rwdg(B2, B3))); // rwdg is associative
        CHECK(rwdg(B1, B1) == vec3d{});       // wdg=0 for coplanar bivectors
        CHECK(rwdg(B1, B2) == -rwdg(B2, B1)); // rwdg is anticommutative for bivectors
        CHECK(rwdg(B1, B2 + B3) ==
              rwdg(B1, B2) + rwdg(B1, B3)); // rwdg distributes over add.
        CHECK(rwdg(B1 + B2, B3) ==
              rwdg(B1, B3) + rwdg(B2, B3)); // rwdg distributes over add.
        CHECK(rwdg(sd * B1, B2) ==
              rwdg(B1, sd * B2)); // scalars can be factored out of wdg
        CHECK(rwdg(sd * B1, B2) ==
              sd * rwdg(B1, B2));        // scalars can be factored out of wdg
        CHECK(rwdg(s, t) == rwdg(t, s)); // rwdg between pscalars equivalent to smult.
        CHECK(rwdg(s, t) == double(s) * double(t));
        CHECK(rwdg(s, B1) == rwdg(B1, s)); // rwdg between pscalar and bivector
        CHECK(rwdg(s, v1) == rwdg(v1, s)); // rwdg between pscalar and vector
        CHECK(rwdg(sc, I_3d) == sc);       // rwdg between pscalar and scalar
        CHECK(rwdg(I_3d, sc) == sc);       // (=pseudoscalar); I_3d is neutral element
                                           // for the wedge product
        CHECK(rwdg(s, B1) == sd * B1);     // rwdg between pseudoscalar and vector
        // (I_3d is neutral element: s = sd * I_3d)

        // congruence
        CHECK(is_congruent(v1, 2.3 * v1));
        CHECK(is_congruent(v1, -12.3 * v1));
        CHECK(is_congruent(wdg(v1, v2), wdg(v2, v1)));
        CHECK_FALSE(is_congruent(v1, v2));

        // common factor axiom
        auto lhs = rwdg(wdg(v1, z), wdg(v2, z));
        auto rhs = rwdg(wdg(wdg(v1, v2), z), z);
        CHECK(lhs.x == doctest::Approx(rhs.x).epsilon(eps));
        CHECK(lhs.y == doctest::Approx(rhs.y).epsilon(eps));
        CHECK(lhs.z == doctest::Approx(rhs.z).epsilon(eps));
        CHECK(is_congruent(rwdg(wdg(v1, z), wdg(v2, z)), z));
    }

    TEST_CASE("MVec3d: geometric product - basic properties")
    {
        fmt::println("MVec3d: geometric product - basic properties");

        vec3d v1{1.0, 2.0, 1.0};
        vec3d v2{0.5, 3.0, 2.0};
        vec3d v3{-2.0, 6.0, 3.0};
        bivec3d B1{-4, 2, 1};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3d{sd};
        auto t = scalar3d{st};

        // nrm_sq(v1) = dot(v1,v1)  for every vector
        CHECK(nrm_sq(v1) == dot(v1, v1));

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

    TEST_CASE("MVec3d: one_3d as geometric-product identity")
    {
        fmt::println("MVec3d: one_3d as geometric-product identity");

        auto v = vec3d{1.0, 2.0, 1.0};
        auto B = bivec3d{-4.0, 2.0, 1.0};
        auto ps = pscalar3d{-3.0};
        auto mv = mvec3d{scalar3d{4.0}, v, B, ps};
        auto mv_e = mvec3d_e{scalar3d{4.0}, B};

        // scalar one_3d is the geometric-product unit
        CHECK(one_3d * v == v);
        CHECK(v * one_3d == v);

        // one_3d_mv is the unit of the full geometric product
        CHECK(one_3d_mv * mv == mv);
        CHECK(mv * one_3d_mv == mv);

        // one_3d_mv_e is the unit of the even-grade geometric product
        CHECK(one_3d_mv_e * mv_e == mv_e);
        CHECK(mv_e * one_3d_mv_e == mv_e);
    }

    TEST_CASE("MVec3d: geometric product tests - vec * vec")
    {
        fmt::println("MVec3d: geometric product tests - vec * vec");

        // ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
        //
        // dot(a,b) = 0.5*(ab + ba)   (symmetric part)
        // wdg(a,b) = 0.5*(ab - ba)   (antisymmetric part)

        vec3d a{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec3d mva{a};
        mvec3d mvb{b};
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

    TEST_CASE("MVec3d: geometric product tests - bivec * vec")
    {
        fmt::println("MVec3d: geometric product tests - bivec * vec");

        // Ab = (b << A) + wdg(A,b) = gr1(Ab) + gr3(Ab)
        //
        // (b << A) = 0.5*(Ab - Aa)   (antisymmetric part)
        // wdg(A,b) = 0.5*(Ab + Aa)   (symmetric part)

        bivec3d A{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -2.0};
        // auto dot_ab = dot(a, b);
        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        mvec3d mva{A};
        mvec3d mvb{b};
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

    TEST_CASE("MVec3d: geometric product tests - vec * bivec")
    {
        fmt::println("MVec3d: geometric product tests - vec * bivec");

        // a*B = (B >> a) + wdg(a,B) = gr1(aB) + gr3(aB)
        //
        // (B >> a) = 0.5*(aB - Ba)   (antisymmetric part)
        // wdg(a,B) = 0.5*(aB + Ba)   (symmetric part)

        vec3d a{1.0, 2.0, 3.0};
        bivec3d B{0.5, 3.0, -2.0};
        // auto dot_ab = dot(a, b);
        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec3d mva{a};
        mvec3d mvb{B};
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

    TEST_CASE("MVec3d: geometric product - combinatorial tests")
    {
        fmt::println("MVec3d: geometric product - combinatorial tests");

        scalar3d s1{3.0};
        scalar3d s2{-1.5};
        vec3d v1{1.0, 2.0, 3.0};
        vec3d v2{0.5, 3.0, -2.0};
        bivec3d b1{-1.0, -3.0, -5.0};
        bivec3d b2{2.0, 0.5, -4.0};
        pscalar3d ps1{-5.0};
        pscalar3d ps2{2.0};

        // checks all ctor combinations and equivalence of simplified products with
        // geometric product of a fully populated multivector
        //
        // assumption: the geometric product of a fully populated multivector serves
        // as
        //             reference and is correctly implemented
        //             results are relative to full geometric product

        // MVec3d_E * MVec3d
        CHECK(mvec3d{s1, b1} * mvec3d{s1, v1, b1, ps1} ==
              mvec3d{mvec3d_e{s1, b1} * mvec3d{s1, v1, b1, ps1}});

        // MVec3d * MVec3d_E
        CHECK(mvec3d{s1, v1, b1, ps1} * mvec3d{s1, b1} ==
              mvec3d{mvec3d{s1, v1, b1, ps1} * mvec3d_e{s1, b1}});

        // MVec3d_E * MVec3d_E
        CHECK(mvec3d{s1, b1} * mvec3d{s2, b2} ==
              mvec3d{mvec3d_e{s1, b1} * mvec3d_e{s2, b2}});

        // MVec3d_U * MVec3d_U
        CHECK(mvec3d{v1, ps1} * mvec3d{v2, ps2} ==
              mvec3d{mvec3d_u{v1, ps1} * mvec3d_u{v2, ps2}});

        // MVec3d_E * MVec3d_U
        CHECK(mvec3d{s1, b1} * mvec3d{v2, ps2} ==
              mvec3d{mvec3d_e{s1, b1} * mvec3d_u{v2, ps2}});

        // MVec3d_U * MVec3d_E
        CHECK(mvec3d{v1, ps1} * mvec3d{s2, b2} ==
              mvec3d{mvec3d_u{v1, ps1} * mvec3d_e{s2, b2}});

        // MVec3d_E * BiVec3d
        CHECK(mvec3d{s1, b1} * mvec3d{b2} == mvec3d{mvec3d_e{s1, b1} * mvec3d{b2}});

        // BiVec3d * MVec3d_E
        CHECK(mvec3d{b1} * mvec3d{s2, b2} == mvec3d{mvec3d{b1} * mvec3d_e{s2, b2}});

        // MVec3d_E * Vec3d
        CHECK(mvec3d{s1, b1} * mvec3d{v2} == mvec3d{mvec3d_e{s1, b1} * mvec3d{v2}});

        // Vec3d * MVec3d_E
        CHECK(mvec3d{v1} * mvec3d{s2, b2} == mvec3d{mvec3d{v1} * mvec3d_e{s2, b2}});

        // MVec3d_U * BiVec3d
        CHECK(mvec3d{v1, ps1} * mvec3d{b2} == mvec3d{mvec3d_u{v1, ps1} * mvec3d{b2}});

        // BiVec3d * MVec3d_U
        CHECK(mvec3d{b1} * mvec3d{v2, ps2} == mvec3d{mvec3d{b1} * mvec3d_u{v2, ps2}});

        // BiVec3d * BiVec3d
        CHECK(mvec3d{b1} * mvec3d{b2} == mvec3d{b1 * b2});

        // BiVec3d * Vec3d
        CHECK(mvec3d{b1} * mvec3d{v1} == mvec3d{b1 * v1});

        // Vec3d * BiVec3d
        CHECK(mvec3d{v1} * mvec3d{b1} == mvec3d{v1 * b1});

        // Vec3d * Vec3d
        CHECK(mvec3d{v1} * mvec3d{v2} == mvec3d{v1 * v2});

        // PScalar3d * MVec3d
        CHECK(mvec3d{ps1} * mvec3d{s2, v2, b2, ps2} ==
              mvec3d{ps1 * mvec3d{s2, v2, b2, ps2}});

        // MVec3d * Pscalar3d
        CHECK(mvec3d{s1, v1, b1, ps1} * mvec3d{ps2} ==
              mvec3d{mvec3d{s1, v1, b1, ps1} * ps2});

        // PScalar3d * MVec3d_E
        CHECK(mvec3d{ps1} * mvec3d{s2, b2} == mvec3d{ps1 * mvec3d_e{s2, b2}});

        // MVec3d_E * Pscalar3d
        CHECK(mvec3d{s1, b1} * mvec3d{ps2} == mvec3d{mvec3d_e{s1, b1} * ps2});

        // PScalar3d * MVec3d_U
        CHECK(mvec3d{ps1} * mvec3d{v2, ps2} == mvec3d{ps1 * mvec3d_u{v2, ps2}});

        // MVec3d_U * Pscalar3d
        CHECK(mvec3d{v1, ps1} * mvec3d{ps2} == mvec3d{mvec3d_u{v1, ps1} * ps2});

        // PScalar3d * BiVec3d
        CHECK(mvec3d{ps1} * mvec3d{b2} == mvec3d{ps1 * b2});

        // BiVec3d * Pscalar3d
        CHECK(mvec3d{b1} * mvec3d{ps2} == mvec3d{b1 * ps2});

        // PScalar3d * Vec3d
        CHECK(mvec3d{ps1} * mvec3d{v2} == mvec3d{ps1 * v2});

        // Vec3d * Pscalar3d
        CHECK(mvec3d{v1} * mvec3d{ps2} == mvec3d{v1 * ps2});

        // Pscalar3d * Pscalar3d
        CHECK(mvec3d{ps1} * mvec3d{ps2} == mvec3d{ps1 * ps2});

        // Scalar * Scalar
        CHECK(mvec3d{s1} * mvec3d{s2} == mvec3d{s1 * s2});
    }

    TEST_CASE("MVec3d: geometric product tests - recovering vectors from the"
              "geometric product")
    {
        fmt::println("MVec3d: geometric product tests - recovering vectors from the "
                     "geometric product");

        // Two multivectors mv1 and mv2 formed from vectors v1 and v2.
        // (gr0(mv1)==0 && gr1(mv1) != 0 && gr2(mv1)==0 &&
        //  gr0(mv2)==0 && gr1(mv2) != 0 && gr2(mv2)==0 )
        //
        // They are multiplied by the geometric product to form a multivector C
        // C = mv1(v1) * mv2(v2) = mv1 * mv2
        //
        // C contains a scalar part and a bivector part exclusively,
        // the remaining components are zero.
        // (gr0(C) != 0 && gr1(C)==0 && gr2(C) !=0)
        //
        // The scalar part of C represents the parts of v1 and v2
        // that are parallel to each other.
        // The bivector part of C represents the parts of v1 and v2
        // that are perpendicular to each other.
        //
        // multiply C from the right with inv(v2) recovers v1
        // multiply C from the left the the inv(v1) recovers v2

        vec3d a{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -4.0};
        mvec3d mva{a};
        mvec3d mvb{b};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);
        mvec3d C = a * b;
        mvec3d Cm = mva * mvb;
        mvec3d Cd{scalar3d(dot_ab), wdg_ab};

        mvec3d gpr_right = C * mvec3d{inv(b)};
        mvec3d gpr_left = mvec3d{inv(a)} * C;

        // fmt::println("   a                           = {}", a);
        // fmt::println("   b                           = {}", b);
        // fmt::println("   C  = a * b                  = {}", C);
        // fmt::println("   Cm = mva * mvb              = {}", Cm);
        // fmt::println("   Cd = mv(dot(a,b), wdg(a,b)) = {}", Cd);
        // fmt::println("");
        // fmt::println("   C * bi = gpr_right = {}", gpr_right);
        // fmt::println("   ai * C = gpr_left  = {}", gpr_left);
        // fmt::println("   gr1(gpr_right) = a = {}", gr1(gpr_right));
        // fmt::println("   gr1(gpr_left)  = b = {}", gr1(gpr_left));

        CHECK(C == Cm);
        CHECK(C == Cd);
        CHECK(a == gr1(gpr_right));
        CHECK(b == gr1(gpr_left));

        auto m = mvec2d{13.0, -27.0, 3.0, 5.0};
        auto prd = m * inv(m);
        CHECK(abs(gr0(prd) - 1.0) < eps);
        CHECK(nrm(gr1(prd)) < eps);
        CHECK(abs(gr2(prd) - 0.0) < eps);
    }


    TEST_CASE("MVec3d: geometric product - link to inner and outer products")
    {
        fmt::println("MVec3d: geometric product - link to inner and outer products");
        vec3d u{1.0, 2.0, 3.0};

        scalar3d s{3};
        vec3d v{-3.0, 2.5, -0.5};
        bivec3d B{2.5, 3.5, 1.5};
        pscalar3d ps{4.0};

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
        CHECK(u * s == gr1(rwdg(s, dual(u)) + wdg(u, s)));

        // fmt::println("");
        // fmt::println("u * s = {}", u * s);
        // fmt::println("");
        // fmt::println("(s >> u) = {}", (s >> u));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("(s >> u) + wdg(u, s) = {}", (s >> u) + wdg(u, s));
        // fmt::println("gr1((s >> u) + wdg(u, s)) = {}", gr1((s >> u) + wdg(u, s)));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(s,dual(u))= {}", rwdg(s, dual(u)));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("rwdg(s,dual(u)) + wdg(u, s) = {}", rwdg(s, dual(u)) + wdg(u, s));
        // fmt::println("gr1(rwdg(s,dual(u)) + wdg(u, s)) = {}",
        //              gr1(rwdg(s, dual(u)) + wdg(u, s)));
        // fmt::println("");

        CHECK(s * u == gr1((u << s) + wdg(s, u)));
        CHECK(s * u == gr1(rwdg(dual(u), s) + wdg(s, u)));

        // fmt::println("");
        // fmt::println("s * u = {}", s * u);
        // fmt::println("");
        // fmt::println("(u << s) = {}", (u << s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("(u << s) + wdg(s, u) = {}", (u << s) + wdg(s, u));
        // fmt::println("gr1((u << s) + wdg(s, u)) = {}", gr1((u << s) + wdg(s, u)));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(dual(u), s)= {}", rwdg(dual(u), s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("rwdg(ldual(u), s) + wdg(s, u) = {}", rwdg(dual(u), s) + wdg(s,
        // u)); fmt::println("gr1(rwdg(dual(u), s) + wdg(s, u)) = {}",
        //              gr1(rwdg(dual(u), s) + wdg(s, u)));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("vector case: ");
        CHECK(u * v == (v >> u) + wdg(u, v));
        CHECK(u * v == rwdg(v, dual(u)) + wdg(u, v));

        // fmt::println("");
        // fmt::println("u * v = {}", u * v);
        // fmt::println("");
        // fmt::println("(v >> u) = {}", (v >> u));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("(v >> u) + wdg(u, v) = {}", (v >> u) + wdg(u, v));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(v,dual(u))= {}", rwdg(v, dual(u)));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("rwdg(v,dual(u)) + wdg(u, v) = {}", rwdg(v, dual(u)) + wdg(u, v));
        // fmt::println("");

        CHECK(v * u == (u << v) + wdg(v, u));
        CHECK(v * u == rwdg(dual(u), v) + wdg(v, u));

        // fmt::println("");
        // fmt::println("v * u = {}", v * u);
        // fmt::println("");
        // fmt::println("(u << v) = {}", (u << v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("(u << v) + wdg(v, u) = {}", (u << v) + wdg(v, u));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(dual(u), v)= {}", rwdg(dual(u), v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("rwdg(dual(u), v) + wdg(v, u) = {}", rwdg(dual(u), v) + wdg(v,
        // u)); fmt::println("");

        // fmt::println("");
        // fmt::println("bivector case: ");
        CHECK(u * B == (B >> u) + wdg(u, B));
        CHECK(u * B == rwdg(B, dual(u)) + wdg(u, B));

        // fmt::println("");
        // fmt::println("u * B = {}", u * B);
        // fmt::println("");
        // fmt::println("(B >> u) = {}", (B >> u));
        // fmt::println("wdg(u, B) = {}", wdg(u, B));
        // fmt::println("(B >> u) + wdg(u, B) = {}", (B >> u) + wdg(u, B));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(B,dual(u))= {}", rwdg(B, dual(u)));
        // fmt::println("wdg(u, B) = {}", wdg(u, B));
        // fmt::println("rwdg(B,dual(u)) + wdg(u, B) = {}", rwdg(B, dual(u)) + wdg(u, B));
        // fmt::println("");

        CHECK(B * u == (u << B) + wdg(B, u));
        CHECK(B * u == rwdg(dual(u), B) + wdg(B, u));

        // fmt::println("");
        // fmt::println("B * u = {}", B * u);
        // fmt::println("");
        // fmt::println("(u << B) = {}", (u << B));
        // fmt::println("wdg(B, u) = {}", wdg(B, u));
        // fmt::println("(u << B) + wdg(B, u) = {}", (u << B) + wdg(B, u));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(dual(u), B)= {}", rwdg(dual(u), B));
        // fmt::println("wdg(B, u) = {}", wdg(B, u));
        // fmt::println("rwdg(dual(u), B) + wdg(B, u) = {}", rwdg(dual(u), B) + wdg(B,
        // u)); fmt::println("");

        // fmt::println("");
        // fmt::println("pscalar case: ");
        CHECK(u * ps == gr2((ps >> u) + wdg(u, ps)));
        CHECK(u * ps == gr2(rwdg(ps, dual(u)) + wdg(u, ps)));

        // fmt::println("");
        // fmt::println("u * ps = {}", u * ps);
        // fmt::println("");
        // fmt::println("(ps >> u) = {}", (ps >> u));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("(ps >> u) + wdg(u, ps) = {}", (ps >> u) + wdg(u, ps));
        // fmt::println("gr1((ps >> u) + wdg(u, ps)) = {}", gr2((ps >> u) + wdg(u, ps)));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(ps,dual(u))= {}", rwdg(ps, dual(u)));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("rwdg(ps,dual(u)) + wdg(u, ps) = {}",
        //              rwdg(ps, dual(u)) + wdg(u, ps));
        // fmt::println("gr2(rwdg(ps,dual(u)) + wdg(u, ps)) = {}",
        //              gr2(rwdg(ps, dual(u)) + wdg(u, ps)));
        // fmt::println("");

        CHECK(ps * u == gr2((u << ps) + wdg(ps, u)));
        CHECK(ps * u == gr2(rwdg(dual(u), ps) + wdg(ps, u)));

        // fmt::println("");
        // fmt::println("ps * u = {}", ps * u);
        // fmt::println("");
        // fmt::println("(u << ps) = {}", (u << ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("(u << ps) + wdg(ps, u) = {}", (u << ps) + wdg(ps, u));
        // fmt::println("gr1((u << ps) + wdg(ps, u)) = {}", gr2((u << ps) + wdg(ps, u)));
        // fmt::println("");
        // fmt::println("dual(u) = {}", dual(u));
        // fmt::println("rwdg(dual(u), ps)= {}", rwdg(dual(u), ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("rwdg(dual(u), ps) + wdg(ps, u) = {}",
        //              rwdg(dual(u), ps) + wdg(ps, u));
        // fmt::println("gr2(rwdg(dual(u), ps) + wdg(ps, u)) = {}",
        //              gr2(rwdg(dual(u), ps) + wdg(ps, u)));
        // fmt::println("");

        // fmt::println("");
    }

    TEST_CASE("MVec3d: geometric product tests - equivalence tests")
    {
        fmt::println("MVec3d: geometric product tests - equivalence tests");

        vec3d a{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -4.0};
        mvec3d mva{a};
        mvec3d mvb{b};

        bivec3d A{1.0, 2.0, 3.0};
        bivec3d B{0.5, 3.0, -4.0};
        mvec3d mvA{A};
        mvec3d mvB{B};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec3d_e ab = a * b;
        mvec3d abm = mva * mvb;
        mvec3d abd{scalar3d{dot_ab}, wdg_ab};

        mvec3d_u Ab = A * b;
        mvec3d Abm = mvA * mvb;
        mvec3d Abd{b_onto_A, wdg_Ab};

        mvec3d_u aB = a * B;
        mvec3d aBm = mva * mvB;
        mvec3d aBd{B_by_a, wdg_aB};

        // fmt::println("");
        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   ab  = mvec3d_e(a * b)            = {}", ab);
        // fmt::println("   abm = mva * mvb                  = {}", abm);
        // fmt::println("   abd = mvec3d(dot(a,b), wdg(a,b)) = {}", abd);
        // fmt::println("");
        // fmt::println("   A                                = {}", A);
        // fmt::println("   mvA                              = {}", mvA);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   Ab  = mvec3d_u(A * b)            = {}", Ab);
        // fmt::println("   Abm = mvA * mvb                  = {}", Abm);
        // fmt::println("   Abd = mvec3d((b << A), wdg(A,b)) = {}", Abd);
        // fmt::println("");
        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   B                                = {}", B);
        // fmt::println("   mvB                              = {}", mvB);
        // fmt::println("   aB  = mvec3d_u(a * B)            = {}", aB);
        // fmt::println("   aBm = mva * mvB                  = {}", aBm);
        // fmt::println("   aBd = mvec3d((B >> a), wdg(a,B)) = {}", aBd);
        // fmt::println("");

        CHECK(gr0(ab) == gr0(abm));
        CHECK(gr1(abm) == vec3d{});
        CHECK(gr2(ab) == gr2(abm));
        CHECK(gr3(abm) == pscalar3d{0.0});

        CHECK(gr0(ab) == gr0(abd));
        CHECK(gr1(abd) == vec3d{});
        CHECK(gr2(ab) == gr2(abd));
        CHECK(gr3(abd) == pscalar3d{0.0});

        CHECK(gr0(Abm) == 0);
        CHECK(gr1(Ab) == gr1(Abm));
        CHECK(gr2(Abm) == bivec3d{});
        CHECK(gr3(Ab) == gr3(Abm));

        CHECK(gr0(Abd) == 0);
        CHECK(gr1(Ab) == gr1(Abd));
        CHECK(gr2(Abd) == bivec3d{});
        CHECK(gr3(Ab) == gr3(Abd));

        CHECK(gr0(aBm) == 0);
        CHECK(gr1(aB) == gr1(aBm));
        CHECK(gr2(aBm) == bivec3d{});
        CHECK(gr3(aB) == gr3(aBm));

        CHECK(gr0(aBd) == 0);
        CHECK(gr1(aB) == gr1(aBd));
        CHECK(gr2(aBd) == bivec3d{});
        CHECK(gr3(aB) == gr3(aBd));
    }

    TEST_CASE("MVec3d: assignment tests")
    {
        fmt::println("MVec3d: assignment tests");

        vec3d v1{1.0, 2.0, 3.0};
        vec3d v2{0.5, 1.0, 1.5};
        vec3d v3{0.5, 1.0, -4.5};
        bivec3d b1{1.0, 2.0, 3.0};

        mvec3d mv1{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        mvec3d mv2{0.0, 0.5, 1.0, 1.5, 11.5, 15.5, 6.0, 61.5};
        mvec3d mv3{mv1};
        mvec3d mv4 = mv2;

        mvec3d mv5(scalar3d(5.0));
        mvec3d mv6(pscalar3d(6.0));
        mvec3d mv7{v1};
        mvec3d mv8{b1};
        mvec3d mv9{scalar3d{dot(v1, v3)}, wdg(v1, v3)};

        mvec3d mv10{v1, pscalar3d(10.0)};
        // This must not compile! Implict conversion to vec3d possible
        // possible solution: explicitly deleted constructor for MVec3d
        // mvec3d mv11{b1, pscalar3d_t(10.0)};

        // this does not compile (which is fine, a base cannot convert to derived)
        // mvec3d mv12{scalar3d(10.0), v1};

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
        // fmt::println("   v1 = {}", v1);
        // fmt::println("   mv7 = v1 = {}", mv7);
        // fmt::println("   b1 = {}", b1);
        // fmt::println("   mv8 = b1 = {}", mv8);
        // fmt::println("");
        // fmt::println("   mv9 = {}", mv9);
        // fmt::println("   mv10 = {}", mv10);

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

    TEST_CASE("MVec3d: bivector product properties")
    {
        fmt::println("MVec3d: bivector product properties");

        bivec3d b1{1.0, 2.0, 3.0};
        mvec3d mb1{b1};
        bivec3d b2{-3.0, 1.0, 2.0};
        mvec3d mb2{b2};

        auto gpr12_m = mb1 * mb2;
        auto gpr21_m = mb2 * mb1;
        auto gpr12_m_sym = 0.5 * (gpr12_m + gpr21_m);
        auto gpr12_m_asym = 0.5 * (gpr12_m - gpr21_m);

        auto gpr12_d = b1 * b2;
        auto gpr21_d = b2 * b1;
        auto gpr12_d_sym = 0.5 * (gpr12_d + gpr21_d);
        auto gpr12_d_asym = 0.5 * (gpr12_d - gpr21_d);

        // fmt::println("   b1  = {}", b1);
        // fmt::println("   mb1 = {}", mb1);
        // fmt::println("   b2  = {}", b2);
        // fmt::println("   mb2 = {}", mb2);
        // fmt::println("");
        // fmt::println("   dot(b1, b2) = {}", dot(b1, b2));
        // fmt::println("   cmt(b1, b2) = {}", cmt(b1, b2));
        // fmt::println("   dot(b2, b1) = {}", dot(b2, b1));
        // fmt::println("   cmt(b2, b1) = {}", cmt(b2, b1));
        // fmt::println("");
        // fmt::println("   gpr12_m = mb1 * mb2 = {}", gpr12_m);
        // fmt::println("   gpr21_m = mb2 * mb1 = {}", gpr21_m);
        // fmt::println("   gpr12_m_sym  = 0.5*(gpr12_d + gpr21_d) = {}",
        // gpr12_m_sym); fmt::println("   gpr12_m_asym = 0.5*(gpr12_m - gpr21_m) =
        // {}", gpr12_m_asym); fmt::println(""); fmt::println("   gpr12_d = b1 * b2 =
        // {} ", gpr12_d); fmt::println("   gpr21_d = b2 * b1) = {} ", gpr21_d);
        // fmt::println("   gpr12_d_sym  = 0.5*(gpr12_d + gpr21_d) = {}",
        // gpr12_d_sym); fmt::println("   gpr12_d_asym = 0.5*(gpr12_d - gpr21_d) =
        // {}", gpr12_d_asym); fmt::println("");

        CHECK(gr2(mb1) == b1);

        // just to silence unused variable warnings
        CHECK(gpr12_m_sym == 0.5 * (gpr12_m + gpr21_m));
        CHECK(gpr12_m_asym == 0.5 * (gpr12_m - gpr21_m));
        CHECK(gpr12_d_sym == 0.5 * (gpr12_d + gpr21_d));
        CHECK(gpr12_d_asym == 0.5 * (gpr12_d - gpr21_d));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d_E<T> and MVec3d_U<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d_E/_U: modelling even and odd parts of 3d algebra - basics")
    {
        fmt::println("MVec3d_E/_U: modelling even and odd parts of 3d algebra - basics");

        // defining a complex number in all three forms as multivector
        auto u = normalize(vec3d{1.0, 0.0, 0.0});
        auto v =
            normalize(vec3d(std::cos(pi / 12), std::sin(pi / 12), 0.0)); // unit vec +15%
        auto angle_uv = angle(u, v);
        auto B = wdg(u, v); // normalized bivector describing the plane spanned by u and v

        auto my_exp = exp(-B * angle_uv);
        auto my_rot = get_rotor(B, 2.0 * angle_uv);

        // definition of rotor used here: B = u^v
        // => B determines the meaning of the positive sign of the rotation
        //
        auto R_m =
            mvec3d(exp(-B * angle_uv)); // Rotor formed by u and v (normalized bivector)
        auto Rr_m = mvec3d(rev(R_m));   // and its reverse

        auto c = vec3d{1.0, 1.0, 1.0};
        auto c_m = mvec3d{c};

        auto c_tmp_m = R_m * c_m;
        auto c_rot_m = c_tmp_m * Rr_m;

        auto R = exp(-B * angle_uv); // Rotor formed by u and v (normalized bivector)
        auto Rr = rev(R);            // and its reverse

        auto c_tmp_l = R * c;
        auto c_rot_u_l = c_tmp_l * Rr;
        auto c_rot_l = gr1(c_rot_u_l);
        // due to symmetry of R and Rr the gr3(c_rot) part will be zero
        // and thus can be assumed to be zero for further computations

        auto c_tmp_r = c * Rr;
        auto c_rot_u_r = R * c_tmp_r;
        auto c_rot_r = gr1(c_rot_u_r);
        // due to symmetry of R and Rr the gr3(c_rot) part will be zero
        // and thus can be assumed to be zero for further computations

        auto angle_c_c_rot = angle(c, c_rot_l); // not that easy in 3D!
        // (angle in plane of both vectors is not the angle in the plane
        // represented by the bivector!)
        // => requires projection of vectors onto plane and then taking
        // the angle between the projected vectors to be correct (bivector angle!)

        auto c_proj = project_onto(c, B);
        auto c_rot_proj = project_onto(c_rot_l, B);
        auto angle_proj = angle(c_proj, c_rot_proj);

        // fmt::println("   u                     = {: .3}", u);
        // fmt::println("   v                     = {: .3}", v);
        // fmt::println("   B = u^v = wdg(u,v)    = {: .3}", B);
        // fmt::println("   angle(u,v)            = {: .3}°", rad2deg(angle_uv));
        // fmt::println("   sin(angle(u,v))       = {: .3}", std::sin(angle_uv));
        // fmt::println("");
        // fmt::println("   c                     = {: .3}", c);
        // fmt::println("");
        // fmt::println("Implemented as full multivector operation:");
        // fmt::println("   R_m  = mvec3d(exp(-B,angle_uv))  = {: .3}", R_m);
        // fmt::println("   Rr_m = rev(R_m)                  = {: .3}", Rr_m);
        // fmt::println("   Rr_m*R_m                         = {: .3}", Rr_m * R_m);
        // fmt::println("   c_m                              = {: .3}", c_m);
        // fmt::println("   c_tmp_m = R_m*c_m                = {: .3}", c_tmp_m);
        // fmt::println("   c_rot_m = c_tmp_m*Rr_m           = {: .3}", c_rot_m);
        // fmt::println("   gr1(c_rot_m)                     = {: .3}", gr1(c_rot_m));
        // fmt::println("");
        // fmt::println("Implemented as reduced grade multivector operation:");
        // fmt::println("   R  = exp(-B,angle_uv)            = {: .3}", R);
        // fmt::println("   Rr = rev(R)                      = {: .3}", Rr);
        // fmt::println("   my_exp = exp(-B, angle_uv)       = {: .3}", my_exp);
        // fmt::println("   my_rot = get_rotor(B, 2*angle_uv)    = {: .3}", my_rot);
        // fmt::println("");
        // fmt::println("Left multiplication of rotor first:");
        // fmt::println("   c_tmp_l = R*c            = {: .3}", c_tmp_l);
        // fmt::println("   c_rot_u_l = c_tmp_l*Rr   = {: .3}", c_rot_u_l);
        // fmt::println("   c_rot_l = gr1(c_rot_u_l) = {: .3}", c_rot_l);
        // fmt::println("");
        // fmt::println("Right multiplication of rotor first:");
        // fmt::println("   c_tmp_r = c*Rr           = {: .3}", c_tmp_r);
        // fmt::println("   c_rot_u_r = R*c_tmp_r    = {: .3}", c_rot_u_r);
        // fmt::println("   c_rot_r = gr1(c_rot_u_r) = {: .3}", c_rot_r);
        // fmt::println("");
        // fmt::println("   angle(c, c_rot_l) = {: .3}°", rad2deg(angle_c_c_rot));
        // fmt::println("   angle(projected)  = {: .3}°", rad2deg(angle_proj));
        // fmt::println("");
        // fmt::println("direct calclulation:");
        // fmt::println("   c_rot = rotate(c,R)          = {: .3}", rotate(c, R));

        CHECK(std::abs(nrm(rotate(c, R)) - nrm(c)) < eps * nrm(c));
        auto const rot_c = rotate(c, R);
        auto const rot_c_m = gr1(c_rot_m);
        auto const rot_c_opt = rotate_opt(c, R);
        CHECK(std::abs(rot_c_m.x - rot_c.x) < eps * nrm(rot_c));
        CHECK(std::abs(rot_c_m.y - rot_c.y) < eps * nrm(rot_c));
        CHECK(std::abs(rot_c_m.z - rot_c.z) < eps * nrm(rot_c));
        CHECK(std::abs(rot_c_opt.x - rot_c.x) < eps * nrm(rot_c));
        CHECK(std::abs(rot_c_opt.y - rot_c.y) < eps * nrm(rot_c));
        CHECK(std::abs(rot_c_opt.z - rot_c.z) < eps * nrm(rot_c));
        // using I_3d approach:
        auto const rot_x_I =
            rotate(vec3d{1.0, 0.0, 0.0}, get_rotor(e3_3d * I_3d, pi / 4));
        auto const ref_xy = normalize(vec3d{1.0, 1.0, 0.0});
        CHECK(std::abs(rot_x_I.x - ref_xy.x) < eps * nrm(ref_xy));
        CHECK(std::abs(rot_x_I.y - ref_xy.y) < eps * nrm(ref_xy));
        CHECK(std::abs(rot_x_I.z - ref_xy.z) < eps * nrm(ref_xy));
        // using a bivector directly:
        auto const rot_x_bv = rotate(vec3d{1.0, 0.0, 0.0}, get_rotor(e12_3d, pi / 4));
        CHECK(std::abs(rot_x_bv.x - ref_xy.x) < eps * nrm(ref_xy));
        CHECK(std::abs(rot_x_bv.y - ref_xy.y) < eps * nrm(ref_xy));
        CHECK(std::abs(rot_x_bv.z - ref_xy.z) < eps * nrm(ref_xy));

        // direct rotation of a bivector
        auto const rot_bv = rotate(bivec3d{0.0, 0.0, 1.0}, get_rotor(e23_3d, pi / 2));
        CHECK(std::abs(rot_bv.x - (-e31_3d).x) < eps * nrm(-e31_3d));
        CHECK(std::abs(rot_bv.y - (-e31_3d).y) < eps * nrm(-e31_3d));
        CHECK(std::abs(rot_bv.z - (-e31_3d).z) < eps * nrm(-e31_3d));

        // example see Macdonald "Linear and Geometric Algebra", Exercise 7.12, p. 127
        auto Bv =
            wdg(e2_3d, e1_3d + std::sqrt(3.0) * e3_3d); // bivector describing the plane
        CHECK(abs(nrm(Bv) - 2.0) < eps * 2.0);
        auto const rot_Bv = rotate(Bv, get_rotor(e31_3d, pi / 3));
        auto const ref_Bv = -2.0 * e12_3d;
        CHECK(std::abs(rot_Bv.x - ref_Bv.x) < eps * nrm(ref_Bv));
        CHECK(std::abs(rot_Bv.y - ref_Bv.y) < eps * nrm(ref_Bv));
        CHECK(std::abs(rot_Bv.z - ref_Bv.z) < eps * nrm(ref_Bv));
        auto const rot_Bv_opt = rotate_opt(Bv, get_rotor(e31_3d, pi / 3));
        CHECK(std::abs(rot_Bv_opt.x - rot_Bv.x) < eps * nrm(rot_Bv));
        CHECK(std::abs(rot_Bv_opt.y - rot_Bv.y) < eps * nrm(rot_Bv));
        CHECK(std::abs(rot_Bv_opt.z - rot_Bv.z) < eps * nrm(rot_Bv));

        // just to silence unused variable warnings
        CHECK(my_exp == exp(-B * angle_uv));
        CHECK(my_rot == get_rotor(B, 2.0 * angle_uv));
        CHECK(c_rot_r == gr1(c_rot_u_r));
        CHECK(angle_c_c_rot == angle(c, c_rot_l));
        CHECK(angle_proj == angle(c_proj, c_rot_proj));
    }

    TEST_CASE("MVec3d: dualization")
    {
        fmt::println("MVec3d: dualization");

        vec3d v{1.0, 2.0, 3.0};                                   // 3d vector
        bivec3d B{10.0, 20.0, 30.0};                              // 3d bivector
        mvec3d vm{100.0, 1.0, 2.0, 3.0, 10.0, 20.0, 30.0, 300.0}; // full 3d multivector
        mvec3d vm2{-20.0, 3.0,  7.0,    -4.0,
                   -1.0,  40.0, -330.0, -70.0}; // full 3d multivector

        // full 3d multivector - even content
        mvec3d vm_even{100.0, 0.0, 0.0, 0.0, 10.0, 20.0, 30.0, 0.0};
        // even-grade 3d multivector
        mvec3d_e vm_E{100.0, 10.0, 20.0, 30.0};

        // full 3d multivector - odd content
        mvec3d vm_odd{0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 300.0};
        // odd-grade 3d multivector
        mvec3d_u vm_U{1.0, 2.0, 3.0, 300.0};

        ////////////////////////////////////////////////////////////////////////////////
        // duality (as defined in Lengyel, "PGA illuminated"):
        ////////////////////////////////////////////////////////////////////////////////
        //
        // dual(A) = cmpl(A) in spaces of odd dimension
        //         = r_cmpl(A) in spaces of even dimension
        //

        auto vm_dual_manual = rev(vm) * I_3d_mv;
        auto vm_dual = dual(vm);

        auto vm_dual_even_manual = rev(vm_even) * I_3d_mv;
        auto vm_dual_even = dual(vm_even);

        auto vm_dual_odd_manual = rev(vm_odd) * I_3d_mv;
        auto vm_dual_odd = dual(vm_odd);

        // result is odd, naming chosen for consistency
        auto vm_dual_manual_E = rev(vm_E) * I_3d;
        auto vm_dual_E = dual(vm_E);

        // result is even, naming chosen for consistency
        auto vm_dual_manual_U = rev(vm_U) * I_3d_mv_u;
        auto vm_dual_U = dual(vm_U);

        auto v_dual_manual = rev(v) * I_3d;
        auto v_dual = dual(v);

        auto B_dual_manual = rev(B) * I_3d;
        auto B_dual = dual(B);

        // fmt::println("");
        // fmt::println("   I_3d                   = {}", I_3d);
        // fmt::println("   I_3d_mv                  = {}", I_3d_mv);
        // fmt::println("   I_3d_mv_u                = {}", I_3d_mv_u);
        // fmt::println("");
        // fmt::println("   v                      = {}", v);
        // fmt::println("   B                      = {}", B);
        // fmt::println("");
        // fmt::println("   vm                     = {}", vm);
        // fmt::println("   rev(vm) * I_3d_mv        = {}", vm_dual_manual);
        // fmt::println("   dual(vm)               = {}", vm_dual);
        // fmt::println("");
        // fmt::println("   vm_even                = {}", vm_even);
        // fmt::println("   rev(vm_even) * I_3d_mv   = {}", vm_dual_even_manual);
        // fmt::println("   dual(vm_even)          = {}", vm_dual_even);
        // fmt::println("");
        // fmt::println("   vm_E                   = {}", vm_E);
        // fmt::println("   rev(vm_E) * I_3d       = {}", vm_dual_manual_E);
        // fmt::println("   dual(vm_E)             = {}", vm_dual_E);
        // fmt::println("");
        // fmt::println("   vm_odd              = {}", vm_odd);
        // fmt::println("   rev(vm_odd) * I_3d_mv = {}", vm_dual_odd_manual);
        // fmt::println("   dual(vm_odd)        = {}", vm_dual_odd);
        // fmt::println("");
        // fmt::println("   vm_U                   = {}", vm_U);
        // fmt::println("   rev(vm_U) * I_3d_mv_u    = {}", vm_dual_manual_U);
        // fmt::println("   dual(vm_U)             = {}", vm_dual_U);
        // fmt::println("");
        // fmt::println("   v                      = {}", v);
        // fmt::println("   rev(v)*I_3d            = {} - bivec ", v_dual_manual);
        // fmt::println("   dual(v)                = {} - bivec ", v_dual);
        // fmt::println("");
        // fmt::println("   B                      = {}", B);
        // fmt::println("   rev(B)*I_3d            = {} - vec", B_dual_manual);
        // fmt::println("   dual(B)                = {} - vec", B_dual);
        // fmt::println("");

        CHECK(vm_dual == vm_dual_manual);
        CHECK(vm_dual_even == vm_dual_even_manual);
        CHECK(vm_dual_odd == vm_dual_odd_manual);
        CHECK(vm_dual_E == vm_dual_manual_E);
        CHECK(vm_dual_U == vm_dual_manual_U);
        CHECK(dual(v) == bivec3d{1.0, 2.0, 3.0});
        CHECK(dual(B) == vec3d{10.0, 20.0, 30.0});
        CHECK(dual(scalar3d(5)) == pscalar3d(5));
        CHECK(dual(pscalar3d(6)) == scalar3d(6));

        // dual properties
        CHECK(dual(3.0 * vm) == 3.0 * dual(vm));
        CHECK(dual(vm + vm2) == dual(vm) + dual(vm2));
        CHECK(dual(dual(vm)) == vm);
        CHECK(dual(I_3d) == 1);
        CHECK(dual(v) == bivec3d(1.0, 2.0, 3.0));

        // inner and outer products are in G^n are dual (A. Macdonald, p. 111):
        //
        // dual(dot(A,B)) = wdg(A,dual(B)), dual(wdg(A,B)) = dot(A,dual(B))

        auto a = vec3d{2.0, -3.0, 2.0};
        auto b = vec3d{-5.0, 1.0, -7.0};

        // fmt::println("   a                 = {}", a);
        // fmt::println("   b                 = {}", b);
        // fmt::println("   dual(b)           = {}", dual(b));
        // fmt::println("   dot(a, b)         = {}", dot(a, b));
        // fmt::println("   wdg(a, b)         = {}", wdg(a, b));
        // fmt::println("");
        // fmt::println("   dual(dot(a, b))   = {}", dual(dot(a, b)));
        // fmt::println("   wdg(a, dual(b))   = {}", wdg(a, dual(b)));
        // fmt::println("");
        // fmt::println("   dual(wdg(a, b))   = {}", dual(wdg(a, b)));
        // fmt::println("   dot(a, dual(b))   = {}", dot(a, dual(b)));

        CHECK(dual(scalar3d{dot(a, b)}) == wdg(a, dual(b)));
        CHECK(dual(scalar3d(dot(a, b))) == wdg(a, dual(b)));
        // CHECK(dual(wdg(a, b)) == dot(a, dual(b)));
        CHECK(wdg(a, dual(b)) == dot(a, b) * I_3d);

        // just to silence unused variable warning
        CHECK(v_dual_manual == rev(v) * I_3d);
        CHECK(v_dual == dual(v));
        CHECK(B_dual_manual == rev(B) * I_3d);
        CHECK(B_dual == dual(B));
    }

    TEST_CASE("MVec3d: metric / antimetric exomorphisms (G, rG) - comparison table")
    {
        fmt::println(
            "MVec3d: metric / antimetric exomorphisms (G, rG) - comparison table");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Exomorphism comparison table for ega3d = G(3,0,0), det(metric) = +1.
        //
        // A metric algebra carries TWO exomorphisms (Lengyel, "PGA Illuminated"):
        //   - the metric      G : a WEDGE exomorphism
        //         G(wdg(a,b))   == wdg(G(a),G(b))
        //   - the antimetric rG : an ANTIWEDGE exomorphism
        //         rG(rwdg(a,b)) == rwdg(rG(a),rG(b))
        // related by  G * rG = det(metric) * I. For ega3d det = +1 and the metric is the
        // identity, so rG == G == I; both verdicts are trivially YES.
        //
        // For each basis blade e_S (vectors in index set S):
        //   P(e_S) = prod_{i in S} m_i   (pure product of the vector squares)
        //   Q(e_S) = sigma(k) * P        (the blade square e_S . e_S),
        //            sigma(k) = (-1)^(k(k-1)/2)
        // and g_S = the value STORED in ega3d_metric (ga_usr_consts.hpp), read straight
        // off the array diagonal. EGA stores P, so g_S follows P at every grade; Q now
        // diverges from P at BOTH grade 2 (e23,e31,e12) and grade 3 (e123).
        ////////////////////////////////////////////////////////////////////////////////

        int const m[3] = {1, 1, 1}; // e1^2, e2^2, e3^2
        char const* nm[8] = {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"};
        int const msk[8] = {0, 1, 2, 4, 6, 5, 3, 7}; // bitmask of each basis blade
        int const full = 0b111;
        auto const grade = [](int mask) {
            int k = 0;
            for (int i = 0; i < 3; ++i)
                k += (mask >> i) & 1;
            return k;
        };
        auto const Pof = [&](int mask) {
            int p = 1;
            for (int i = 0; i < 3; ++i)
                if (mask & (1 << i)) p *= m[i];
            return p;
        };
        auto const sigma = [](int k) { return ((k * (k - 1) / 2) & 1) ? -1 : 1; };
        auto const idx_of = [&](int mask) {
            for (int i = 0; i < 8; ++i)
                if (msk[i] == mask) return i;
            return 0;
        };
        auto const Gv = ega3d_metric_view();
        auto const storedM = [&](int mask) {
            int const i = idx_of(mask);
            return Gv[i, i];
        };

        int const det = Pof(full);

        fmt::println("   blade |  k |  P |  Q | g_S | rG_S | g_S*rG_S(=det={:+})", det);
        fmt::println("   ------+----+----+----+-----+------+-------------------");
        for (int i = 0; i < 8; ++i) {
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

        CHECK(wedge_exo);
        CHECK(antiwedge_exo);
        CHECK(det == 1);

        ////////////////////////////////////////////////////////////////////////////////
        // Library handles. ega3d is odd-dimensional, so the single complement cmpl and
        // single dual stand in for the l_/r_ pair (l_cmpl == r_cmpl == cmpl here):
        //   Gx  = cmpl o dual          -- the stored metric (reproduces the array)
        //   rGx = cmpl o Gx o cmpl     -- the antimetric; on a blade e_S this is
        //                                 g_{S^c} * e_S (= complement product)
        ////////////////////////////////////////////////////////////////////////////////
        auto const Gx = [](auto const& X) { return cmpl(dual(X)); };
        auto const rGx = [&Gx](auto const& X) { return cmpl(Gx(cmpl(X))); };

        auto const s = scalar3d{2.0};
        auto const v = vec3d{3.0, 5.0, 7.0};
        auto const B = bivec3d{1.0, 2.0, 3.0};
        auto const ps = pscalar3d{4.0};

        // G * rG = det * I, det = +1  =>  rGx(Gx(X)) == X at every grade
        CHECK(rGx(Gx(s)) == s);
        CHECK(rGx(Gx(v)) == v);
        CHECK(rGx(Gx(B)) == B);
        CHECK(rGx(Gx(ps)) == ps);
        CHECK(rGx(v) == Gx(v)); // identity metric => rG == G here

        // (a) Gx is a WEDGE exomorphism
        CHECK(Gx(wdg(e1_3d, e2_3d)) == wdg(Gx(e1_3d), Gx(e2_3d)));   // v ^ v -> B
        CHECK(Gx(wdg(e1_3d, e23_3d)) == wdg(Gx(e1_3d), Gx(e23_3d))); // v ^ B -> ps
        // (b) rGx is an ANTIWEDGE exomorphism
        CHECK(rGx(rwdg(e23_3d, e31_3d)) == rwdg(rGx(e23_3d), rGx(e31_3d))); // B v B -> v
        CHECK(rGx(rwdg(I_3d, e12_3d)) == rwdg(rGx(I_3d), rGx(e12_3d)));     // ps v B -> B
    }

    TEST_CASE("MVec3d: product tests")
    {
        fmt::println("MVec3d: product tests");
        // tests based on "The inner products of geometric algebra", Leo Dorst


        auto s1 = scalar3d{2.0};
        auto v1 = vec3d{1.0, -3.0, 2.0};
        auto b1 = bivec3d{2.0, -4.0, 1.0};
        auto ps1 = pscalar3d{-2.0};
        auto M1 = mvec3d{s1, v1, b1, ps1};

        auto s2 = scalar3d{-1.0};
        auto v2 = vec3d{2.0, 1.5, -1.0};
        auto b2 = bivec3d{-2.0, 10.5, 3.0};
        auto ps2 = pscalar3d{3.0};
        auto M2 = mvec3d{s2, v2, b2, ps2};

        auto s3 = scalar3d{-10.0};
        auto v3 = vec3d{-2.0, 4.0, -10.0};
        auto b3 = bivec3d{1.0, -2.0, 3.0};
        auto ps3 = pscalar3d{-12.0};
        auto M3 = mvec3d{s3, v3, b3, ps3};

        auto V1 = mvec3d{v1};
        auto V2 = mvec3d{v2};
        auto V3 = mvec3d{v3};
        auto B12 = mvec3d{wdg(v1, v2)};

        CHECK(gr1(V1) == v1);
        CHECK(gr1(V2) == v2);
        CHECK(gr1(V3) == v3);
        CHECK(gr2(B12) == wdg(v1, v2));

        CHECK(nrm_sq(M1) == nrm_sq(s1) + nrm_sq(v1) + nrm_sq(b1) + nrm_sq(ps1));
        CHECK(nrm_sq(M2) == nrm_sq(s2) + nrm_sq(v2) + nrm_sq(b2) + nrm_sq(ps2));
        CHECK(nrm_sq(M3) == nrm_sq(s3) + nrm_sq(v3) + nrm_sq(b3) + nrm_sq(ps3));

        CHECK(dot(M1, M2) == dot(s1, s2) + dot(v1, v2) + dot(b1, b2) + dot(ps1, ps2));
        CHECK(dot(M1, M1) == nrm_sq(M1));
        CHECK(dot(M2, M2) == nrm_sq(M2));
        CHECK(dot(M3, M3) == nrm_sq(M3));

        // wedge product
        CHECK(wdg(v1, v2) == gr2(wdg(mvec3d(v1), mvec3d(v2))));
        // check full permissible range of arguments, even the ones delivering 0)
        CHECK(wdg(b1, b2) == scalar3d(0.0));
        CHECK(wdg(v1, ps1) == scalar3d(0.0));
        CHECK(wdg(ps1, v1) == scalar3d(0.0));
        CHECK(wdg(b1, ps1) == scalar3d(0.0));
        CHECK(wdg(ps1, b1) == scalar3d(0.0));
        CHECK(wdg(ps1, ps1) == scalar3d(0.0));
        CHECK(wdg(b1, b2) == nrm_sq(b1 * b2 - cmt(b1, b2) + dot(b1, b2)));
        CHECK(wdg(v1, ps1) == nrm_sq(0.5 * (v1 * ps1 + rev(ps1) * v1)));
        CHECK(wdg(ps1, v1) == nrm_sq(0.5 * (ps1 * v1 + v1 * rev(ps1))));


        CHECK((mvec3d{s2, v2, b2, pscalar3d{0.0}} << mvec3d{s1}) ==
              mvec3d{value_t(s1 * s2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

        CHECK((mvec3d{v2} << mvec3d{s1}) ==
              mvec3d{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

        CHECK((wdg(v2, v3) >> v1) == wdg(v1 << v2, v3) + wdg(gr_inv(v2), v1 << v3));
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
        // 1st with vectors lifted to full multivectors
        CHECK(V1 * M2 == (M2 >> V1) + wdg(V1, M2));
        CHECK(M2 * V1 == (V1 << M2) + wdg(M2, V1));
        CHECK((M2 >> V1) == -(V1 << gr_inv(M2))); // replacing the right contraction
        CHECK(V1 * M2 == -(V1 << gr_inv(M2)) + wdg(V1, M2)); // yep, replacing works

        // 2nd with vectors and multivectors directly
        CHECK(v1 * M2 == (M2 >> v1) + wdg(v1, M2));
        CHECK(M2 * v1 == (v1 << M2) + wdg(M2, v1));
        CHECK((M2 >> v1) == -(v1 << gr_inv(M2)));
        CHECK(wdg(v1, M2) == 0.5 * (v1 * M2 + gr_inv(M2) * v1));
        CHECK(wdg(M2, v1) == 0.5 * (M2 * v1 + v1 * gr_inv(M2)));

        // 3rd with vectors and bivectors directly
        CHECK(v1 * b1 == (b1 >> v1) + wdg(v1, b1));
        CHECK(b1 * v1 == (v1 << b1) + wdg(b1, v1));
        CHECK((b1 >> v1) == -(v1 << gr_inv(b1)));


        // express contraction and wedge product through the geometric product
        CHECK(wdg(V1, M2) == 0.5 * (V1 * M2 + gr_inv(M2) * V1));
        CHECK(wdg(M2, V1) == 0.5 * (M2 * V1 + V1 * gr_inv(M2)));

        // product dualities
        CHECK(dual((v1 << v2)) == wdg(v1, dual(v2)));
        CHECK(dual(wdg(v1, v2)) == (v1 << dual(v2)));
        CHECK(dual((M1 << M2)) == wdg(M1, dual(M2)));
        CHECK(dual(wdg(M1, M2)) == (M1 << dual(M2)));

        // cross-check direct implementation of rwdg by comparing with wdg
        CHECK(rwdg(M1, M2) == cmpl(wdg(cmpl(M1), cmpl(M2))));

        // commutator and contraction are equal for a specific case with bivectors
        CHECK(cmt(v1, b1) == (b1 >> v1));
        CHECK(cmt(b1, v1) == (v1 << b1));
        CHECK(cmt(v1, b1) == -cmt(b1, v1));
    }

    TEST_CASE("MVec3d: simple applications, complements, contraction, expansions")
    {
        fmt::println("MVec3d: simple applications, complements, contraction, expansions");

        auto s1 = scalar3d{2.0};
        auto v1 = vec3d{1.0, -3.0, 2.0};
        auto b1 = bivec3d{2.0, -4.0, 1.0};
        auto ps1 = pscalar3d{-2.0};
        auto M1 = mvec3d{s1, v1, b1, ps1};


        // simple projections
        auto v = vec3d{4.0, 3.5, 5.0};
        // auto B = e31_3d;
        auto B = e31_3d + 0.1 * e12_3d;
        // auto B = e31_3d + 0.01 * e12_3d + 0.1 * e23_3d;


        auto v_in_B = project_onto(v, B);
        auto v_perp_B = reject_from(v, B);

        // fmt::println("   v         = {:.3f}", v);
        // fmt::println("   B         = {:.3f}", B);
        // fmt::println("");
        // fmt::println("   v_in_B    = {:.3f}", v_in_B);
        // fmt::println("   v_perp_B  = {:.3f}", v_perp_B);
        // fmt::println("   v << B    = {:.3f}", v << B);
        // fmt::println("   B >> v    = {:.3f}", B >> v);
        // fmt::println("");

        // x^b = 0 for every point x in B
        CHECK(nrm_sq(wdg(v_in_B, B)) < eps);

        // v_perp_B should be proportional to the normal vector n of B
        // n = dual(B) or n = cmpl(B)
        // thus, wdg(v_perp_B, n) == 0 is required
        CHECK(nrm_sq(wdg(v_perp_B, dual(B))) < eps);
        CHECK(nrm_sq(wdg(v_perp_B, cmpl(B))) < eps);

        // v_in_B and v_perp_B should be perpendicular to each other
        CHECK(nrm_sq(dot(v_in_B, v_perp_B)) < eps);

        // v should be the sum of v_in_B and v_perp_B
        CHECK(v == v_in_B + v_perp_B);

        // check complements
        // fmt::println("");
        // auto s = scalar3d(1);
        // auto e1 = e1_3d;
        // auto e2 = e2_3d;
        // auto e3 = e3_3d;
        // auto e23 = e23_3d;
        // auto e31 = e31_3d;
        // auto e12 = e12_3d;
        // auto ps = I_3d;
        // fmt::println("   u:  s, e1, e2, e3, e23, e31, e12, ps \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}, \n"
        //              "{:.3f}, {:.3f}, {:.3f}, {:.3f}",
        //              s, e1, e2, e3, e23, e31, e12, ps);
        // fmt::println("");
        // fmt::println("   cmpl(u):  cmpl(s), cmpl(e1), cmpl(e2), cmpl(e3) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              cmpl(s), cmpl(e1), cmpl(e2), cmpl(e3));
        // fmt::println("");
        // fmt::println("   cmpl(u):  cmpl(e23), cmpl(e31), cmpl(e12), cmpl(ps) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              cmpl(e23), cmpl(e31), cmpl(e12), cmpl(ps));
        // fmt::println("");
        // fmt::println("   rev(u)*I_3d:  rev(s)*I_3d, rev(e1)*I_3d, rev(e2)*I_3d, "
        //              "rev(e3)*I_3d \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              rev(s) * I_3d, rev(e1) * I_3d, rev(e2) * I_3d, rev(e3) *
        //              I_3d);
        // fmt::println("");
        // fmt::println("   rev(u)*I_3d:  rev(e23)*I_3d, rev(e31)*I_3d, "
        //              "rev(e12)*I_3d, rev(ps)*I_3d \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              rev(e23) * I_3d, rev(e31) * I_3d, rev(e12) * I_3d, rev(ps) *
        //              I_3d);
        // fmt::println("");
        // fmt::println("   I_3d*rev(u):  I_3d*rev(s), I_3d*rev(e1), I_3d*rev(e2), "
        //              "I_3d*rev(e3) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              I_3d * rev(s), I_3d * rev(e1), I_3d * rev(e2), I_3d *
        //              rev(e3));
        // fmt::println("");
        // fmt::println("   I_3d*rev(u):  I_3d*rev(e23), I_3d*rev(e31), "
        //              "I_3d*rev(e13), I_3d*rev(ps) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              I_3d * rev(e23), I_3d * rev(e31), I_3d * rev(e12), I_3d *
        //              rev(ps));
        // fmt::println("");

        CHECK(cmpl(scalar3d(1.0)) == I_3d);
        CHECK(cmpl(scalar3d(1.0)) == rev(scalar3d(1.0)) * I_3d);
        CHECK(cmpl(e1_3d) == e23_3d);
        CHECK(cmpl(e1_3d) == rev(e1_3d) * I_3d);
        CHECK(cmpl(e2_3d) == e31_3d);
        CHECK(cmpl(e2_3d) == rev(e2_3d) * I_3d);
        CHECK(cmpl(e3_3d) == e12_3d);
        CHECK(cmpl(e3_3d) == rev(e3_3d) * I_3d);
        CHECK(cmpl(e23_3d) == e1_3d);
        CHECK(cmpl(e23_3d) == rev(e23_3d) * I_3d);
        CHECK(cmpl(e31_3d) == e2_3d);
        CHECK(cmpl(e31_3d) == rev(e31_3d) * I_3d);
        CHECK(cmpl(e12_3d) == e3_3d);
        CHECK(cmpl(e12_3d) == rev(e12_3d) * I_3d);
        CHECK(cmpl(I_3d) == scalar3d(1.0));
        CHECK(cmpl(I_3d) == rev(I_3d) * I_3d);
        //
        CHECK(cmpl(cmpl(M1)) == M1);
        CHECK(cmpl(cmpl(mvec3d_e(s1, b1))) == mvec3d_e(s1, b1));
        CHECK(cmpl(cmpl(mvec3d_u(v1, ps1))) == mvec3d_u(v1, ps1));
        //
        // complements are defined for basis elements only
        // => magnitude has to be covered separately for non-normalized elements

        // left complements = complements in spaces of odd dimension
        CHECK(wdg(cmpl(s1), s1) / nrm_sq(s1) == I_3d);
        CHECK(wdg(cmpl(v1), v1) / nrm_sq(v1) == I_3d);
        CHECK(wdg(cmpl(b1), b1) / nrm_sq(b1) == I_3d);
        CHECK(wdg(cmpl(ps1), ps1) / nrm_sq(ps1) == I_3d);

        // right complements = complements in spaces of odd dimension
        CHECK(wdg(s1, cmpl(s1)) / nrm_sq(s1) == I_3d);
        CHECK(wdg(v1, cmpl(v1)) / nrm_sq(v1) == I_3d);
        CHECK(wdg(b1, cmpl(b1)) / nrm_sq(b1) == I_3d);
        CHECK(wdg(ps1, cmpl(ps1)) / nrm_sq(ps1) == I_3d);

        // correspondence of complements with geometric products:
        // duals correspond to complements in ega3d (non-degenerate metric = identity)
        CHECK((I_3d * rev(s1)) == dual(s1));
        CHECK((I_3d * rev(v1)) == dual(v1));
        CHECK((I_3d * rev(b1)) == dual(b1));
        CHECK((I_3d * rev(ps1)) == dual(ps1));

        CHECK((rev(s1) * I_3d) == dual(s1));
        CHECK((rev(v1) * I_3d) == dual(v1));
        CHECK((rev(b1) * I_3d) == dual(b1));
        CHECK((rev(ps1) * I_3d) == dual(ps1));

        // check contractions: <<, >> and rwdg( u, compl(v) )
        // fmt::println("   v         = {:.3f}", v);
        // fmt::println("   cmpl(v)   = {:.3f}", cmpl(v));
        // fmt::println("");
        // fmt::println("   v << B = {:.3f}", v << b);
        // fmt::println("   B >> v = {:.3f}", b >> v);
        // fmt::println("");
        // fmt::println("   rwdg(b, cmpl(v))  = {:.3f}", rwdg(b, cmpl(v)));
        // fmt::println("   rwdg(cmpl(v), b)  = {:.3f}", rwdg(cmpl(v), b));
        // fmt::println("");
        CHECK((B >> v) == rwdg(B, cmpl(v)));
        CHECK((v << B) == rwdg(cmpl(v), B));

        // check expansions: v ^ cmpl(b)  and  cmpl(b) ^ v
        // (create new bivector that contains v and is perpendicular to b)
        //
        // also checks the duality correspondence:
        //      cmpl(v >> b) == v ^ cmpl(b)
        //      cmpl(v >> b) == cmpl(b) ^ v
        // fmt::println("   v      = {:.3f}", v);
        // fmt::println("   b      = {:.3f}", b);
        // fmt::println("   v << b = {:.3f}", v << b);
        // fmt::println("   b >> v = {:.3f}", b >> v);
        // fmt::println("");
        // fmt::println("   wdg(cmpl(b), v)         = {:.3f}", wdg(cmpl(b), v));
        // fmt::println("   n=cmpl(wdg(cmpl(b), v)) = {:.3f}", cmpl(wdg(cmpl(b), v)));
        // fmt::println("");
        // fmt::println("   wdg(v, cmpl(b))         = {:.3f}", wdg(v, cmpl(b)));
        // fmt::println("   n=cmpl(wdg(v, cmpl(b))) = {:.3f}", cmpl(wdg(v, cmpl(b))));
        // fmt::println("");
        // fmt::println("   v              = {:.3f}", v);
        // fmt::println("   B              = {:.3f}", B);
        // fmt::println("   l_expand(B, v) = {:.3f}", l_expand(B, v));
        // fmt::println("   r_expand(v, B) = {:.3f}", r_expand(v, B));
        // fmt::println("");

        // check l_expand an r_expand
        CHECK(l_expand(B, v) == wdg(cmpl(B), v));
        CHECK(r_expand(v, B) == wdg(v, cmpl(B)));
        // vector v is in resulting plane defined by the expansion
        CHECK(wdg(v, l_expand(B, v)) == 0.0);
        CHECK(wdg(v, r_expand(v, B)) == 0.0);
        // resulting planes (i.e. their bivectors) are perpendicular to B
        CHECK(dot(B, l_expand(B, v)) == 0.0);
        CHECK(dot(B, r_expand(v, B)) == 0.0);

        // duality of the contraction and the wedge product (based on dual)
        CHECK(dual(v << B) == wdg(v, dual(B)));
        CHECK(dual(B >> v) == wdg(dual(B), v));
        // duality of the contraction and the wedge product (based on complement)
        CHECK(cmpl(v << B) == wdg(v, cmpl(B)));
        CHECK(cmpl(B >> v) == wdg(cmpl(B), v));
    }

    TEST_CASE("EGA3D: congruence tests")
    {
        fmt::println("EGA3D: congruence tests");

        // Test scalars
        scalar3d s1{5.0};
        scalar3d s2{-3.0}; // different sign
        scalar3d s3{2.5};  // same sign as s1
        scalar3d s_zero{0.0};

        CHECK(is_congruent(s1, s2) == true);         // different signs are congruent
        CHECK(is_congruent(s1, s3) == true);         // same signs are congruent
        CHECK(is_congruent(s1, s_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(s_zero, s_zero) == true); // zero vs zero

        // Test vectors
        vec3d v1{1.0, 0.0, 0.0};  // x-axis
        vec3d v2{0.0, 1.0, 0.0};  // y-axis (perpendicular)
        vec3d v3{2.0, 0.0, 0.0};  // parallel to v1
        vec3d v4{-1.5, 0.0, 0.0}; // antiparallel to v1
        vec3d v_zero{0.0, 0.0, 0.0};

        CHECK(is_congruent(v1, v2) == false);        // perpendicular vectors
        CHECK(is_congruent(v1, v3) == true);         // parallel vectors
        CHECK(is_congruent(v1, v4) == true);         // antiparallel vectors
        CHECK(is_congruent(v1, v_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(v_zero, v_zero) == true); // zero vs zero

        // Test bivectors
        bivec3d b1{1.0, 0.0, 0.0};  // e23 plane
        bivec3d b2{0.0, 1.0, 0.0};  // e31 plane (different)
        bivec3d b3{2.0, 0.0, 0.0};  // e23 plane (parallel to b1)
        bivec3d b4{-1.5, 0.0, 0.0}; // e23 plane (antiparallel to b1)
        bivec3d b_zero{0.0, 0.0, 0.0};

        CHECK(is_congruent(b1, b2) == false);        // different planes
        CHECK(is_congruent(b1, b3) == true);         // same plane
        CHECK(is_congruent(b1, b4) == true);         // same plane, opposite orientation
        CHECK(is_congruent(b1, b_zero) == false);    // zero vs non-zero
        CHECK(is_congruent(b_zero, b_zero) == true); // zero vs zero

        // Test pseudoscalars
        pscalar3d p1{1.0};
        pscalar3d p2{2.0};  // different magnitude
        pscalar3d p3{-1.0}; // opposite sign
        pscalar3d p_zero{0.0};

        CHECK(is_congruent(p1, p2) == true);      // all non-zero pseudoscalars congruent
        CHECK(is_congruent(p1, p3) == true);      // different signs still congruent
        CHECK(is_congruent(p1, p_zero) == false); // zero vs non-zero
        CHECK(is_congruent(p_zero, p_zero) == true); // zero vs zero

        // fmt::println("   All EGA3D congruence tests passed");
    }

    TEST_CASE("G<3,0,0>: MVec3d_E, MVec3d_U, and MVec3d formatting tests")
    {
        fmt::println("G<3,0,0>: MVec3d_E, MVec3d_U, and MVec3d formatting tests");

        // Test MVec3d_E formatting (MVec4_t with 4 components: c0, c1, c2, c3)
        MVec3d_E<double> mvec3d_e_val{1.1, 2.2, 3.3, 4.4};

        // Basic output
        std::string mvec3d_e_basic = fmt::format("{}", mvec3d_e_val);
        CHECK(mvec3d_e_basic == "MVec3d_E(1.1,2.2,3.3,4.4)");

        // Two decimal places
        std::string mvec3d_e_two_decimals = fmt::format("{:.2f}", mvec3d_e_val);
        CHECK(mvec3d_e_two_decimals == "MVec3d_E(1.10,2.20,3.30,4.40)");

        // Scientific notation
        std::string mvec3d_e_scientific = fmt::format("{:.2e}", mvec3d_e_val);
        CHECK(mvec3d_e_scientific == "MVec3d_E(1.10e+00,2.20e+00,3.30e+00,4.40e+00)");

        // Contextual usage
        std::string mvec3d_e_contextual = fmt::format("Even MV: {:.2f}", mvec3d_e_val);
        CHECK(mvec3d_e_contextual == "Even MV: MVec3d_E(1.10,2.20,3.30,4.40)");

        // Test MVec3d_U formatting (MVec4_t with 4 components: c0,c1,c2,c3)
        MVec3d_U<double> mvec3d_u_val{5.5, 6.6, 7.7, 8.8};

        // Basic output
        std::string mvec3d_u_basic = fmt::format("{}", mvec3d_u_val);
        CHECK(mvec3d_u_basic == "MVec3d_U(5.5,6.6,7.7,8.8)");

        // Two decimal places
        std::string mvec3d_u_two_decimals = fmt::format("{:.2f}", mvec3d_u_val);
        CHECK(mvec3d_u_two_decimals == "MVec3d_U(5.50,6.60,7.70,8.80)");

        // Scientific notation
        std::string mvec3d_u_scientific = fmt::format("{:.2e}", mvec3d_u_val);
        CHECK(mvec3d_u_scientific == "MVec3d_U(5.50e+00,6.60e+00,7.70e+00,8.80e+00)");

        // Contextual usage
        std::string mvec3d_u_contextual = fmt::format("Odd MV: {:.2f}", mvec3d_u_val);
        CHECK(mvec3d_u_contextual == "Odd MV: MVec3d_U(5.50,6.60,7.70,8.80)");

        // Test MVec3d formatting (MVec8_t with 8 components: c0 through c7)
        MVec3d<double> mvec3d_val{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

        // Basic output
        std::string mvec3d_basic = fmt::format("{}", mvec3d_val);
        CHECK(mvec3d_basic == "MVec3d(1,2,3,4,5,6,7,8)");

        // Two decimal places
        std::string mvec3d_two_decimals = fmt::format("{:.2f}", mvec3d_val);
        CHECK(mvec3d_two_decimals == "MVec3d(1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00)");

        // Scientific notation
        std::string mvec3d_scientific = fmt::format("{:.2e}", mvec3d_val);
        CHECK(mvec3d_scientific == "MVec3d(1.00e+00,2.00e+00,3.00e+00,4.00e+00,"
                                   "5.00e+00,6.00e+00,7.00e+00,8.00e+00)");

        // Contextual usage
        std::string mvec3d_contextual = fmt::format("Full MV: {:.1f}", mvec3d_val);
        CHECK(mvec3d_contextual == "Full MV: MVec3d(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0)");

        fmt::println("   MVec3d_E basic: {}", mvec3d_e_basic);
        fmt::println("   MVec3d_E 2-decimal: {}", mvec3d_e_two_decimals);
        fmt::println("   MVec3d_E scientific: {}", mvec3d_e_scientific);
        fmt::println("   MVec3d_E contextual: {}", mvec3d_e_contextual);

        fmt::println("   MVec3d_U basic: {}", mvec3d_u_basic);
        fmt::println("   MVec3d_U 2-decimal: {}", mvec3d_u_two_decimals);
        fmt::println("   MVec3d_U scientific: {}", mvec3d_u_scientific);
        fmt::println("   MVec3d_U contextual: {}", mvec3d_u_contextual);

        fmt::println("   MVec3d basic: {}", mvec3d_basic);
        fmt::println("   MVec3d 2-decimal: {}", mvec3d_two_decimals);
        fmt::println("   MVec3d scientific: {}", mvec3d_scientific);
        fmt::println("   MVec3d contextual: {}", mvec3d_contextual);
    }

    TEST_CASE("G<3,0,0>: transwedge and regressive transwedge products")
    {
        fmt::println("G<3,0,0>: transwedge and regressive transwedge products");

        auto s = scalar3d{4.0};
        auto v = vec3d{1.0, 2.0, 3.0};
        auto v2 = vec3d{-1.0, 3.0, 4.0};
        auto B = bivec3d{10.0, 20.0, 30.0};
        auto B2 = bivec3d{-20.0, -40.0, -60.0};
        auto ps = pscalar3d{-3.0};

        // check identities of transwedge product
        CHECK(twdg1(ps, v) == ps * v);
        CHECK(twdg1(ps, v) == (v << ps));
        CHECK(twdg1(ps, v) == (ps >> v));
        CHECK(twdg1(ps, v) == rwdg(ps, dual(v))); // contraction expressed as rwdg prod.

        CHECK(twdg1(ps, v) == twdg1(v, ps));

        CHECK(twdg1(v, ps) == v * ps);
        CHECK(twdg1(v, ps) == (ps >> v));
        CHECK(twdg1(v, ps) == (v << ps));
        CHECK(twdg1(v, ps) == rwdg(dual(v), ps)); // contraction expressed as rwdg prod.

        CHECK(twdg1(B, B2) == cmt(B, B2));

        CHECK(twdg1(B, v) == (v << B));
        CHECK(twdg1(v, B) == (B >> v));

        CHECK(twdg1(v, v2) == dot(v, v2));
        CHECK(twdg1(v, v2) == dot(v2, v));

        // check identities of regressive transwedge product
        CHECK(rtwdg1(B, B2) == pscalar3d(value_t(dot(B, B2))));
        CHECK(rtwdg1(B, B2) == pscalar3d(value_t(dot(B2, B))));

        CHECK(rtwdg1(B, v) == wdg(v, dual(B))); // expansion expressed as wdg prod.
        CHECK(rtwdg1(B, v) == -rtwdg1(v, B));   // vector and bivector anticommute
        CHECK(rtwdg1(v, B) == wdg(dual(B), v)); // expansion expressed as wdg prod.

        CHECK(rtwdg1(B, s) == wdg(s, dual(B))); // expansion expressed as wdg prod.
        CHECK(rtwdg1(s, B) == rtwdg1(B, s));    // bivector and scalar commute
        CHECK(rtwdg1(s, B) == wdg(dual(B), s)); // expansion expressed as wdg prod.

        CHECK(rtwdg1(v, v2) == cross(v2, v));  // identity to cross product
        CHECK(rtwdg1(v, v2) == -cross(v, v2)); // identity to cross product
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Dual Mathematical Properties Tests
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<3,0,0>: dual involution properties")
    {
        fmt::println("G<3,0,0>: dual involution properties");

        // Test dual involution: dual(dual(u)) = ±u
        // For EGA3D (odd-dimensional, Euclidean): dual(dual(u)) = +u

        scalar3d s{3.0};
        vec3d v{2.0, 5.0, 7.0};
        bivec3d B{1.0, 2.0, 3.0};
        pscalar3d ps{11.0};
        mvec3d_e M_e{scalar3d{1.0}, bivec3d{2.0, 3.0, 4.0}};
        mvec3d M{scalar3d{1.0}, vec3d{2.0, 3.0, 4.0}, bivec3d{5.0, 6.0, 7.0},
                 pscalar3d{8.0}};

        // dual involution for odd-dimensional: dual(dual(u)) = +u
        CHECK(dual(dual(s)) == s);
        CHECK(dual(dual(v)) == v);
        CHECK(dual(dual(B)) == B);
        CHECK(dual(dual(ps)) == ps);
        CHECK(dual(dual(M_e)) == M_e);
        CHECK(dual(dual(M)) == M);

        fmt::println("dual involution: dual(dual(u)) = +u for odd-dimensional");
    }

    TEST_CASE("G<3,0,0>: complement involution properties")
    {
        fmt::println("G<3,0,0>: complement involution properties");

        // For odd-dimensional algebras: cmpl(cmpl(u)) = +u
        scalar3d s{3.0};
        vec3d v{2.0, 5.0, 7.0};
        bivec3d B{1.0, 2.0, 3.0};
        pscalar3d ps{11.0};
        mvec3d_e M_e{scalar3d{1.0}, bivec3d{2.0, 3.0, 4.0}};
        mvec3d M{scalar3d{1.0}, vec3d{2.0, 3.0, 4.0}, bivec3d{5.0, 6.0, 7.0},
                 pscalar3d{8.0}};

        CHECK(cmpl(cmpl(s)) == s);
        CHECK(cmpl(cmpl(v)) == v);
        CHECK(cmpl(cmpl(B)) == B);
        CHECK(cmpl(cmpl(ps)) == ps);
        CHECK(cmpl(cmpl(M_e)) == M_e);
        CHECK(cmpl(cmpl(M)) == M);

        fmt::println("complement involution: cmpl(cmpl(u)) = +u for odd-dimensional");
    }

    TEST_CASE("G<3,0,0>: complement-dual relationship")
    {
        fmt::println("G<3,0,0>: complement-dual relationship");

        // For EGA3D with identity metric: dual = complement
        scalar3d s{3.0};
        vec3d v{2.0, 5.0, 7.0};
        bivec3d B{1.0, 2.0, 3.0};
        pscalar3d ps{11.0};
        mvec3d_e M_e{scalar3d{1.0}, bivec3d{2.0, 3.0, 4.0}};
        mvec3d M{scalar3d{1.0}, vec3d{2.0, 3.0, 4.0}, bivec3d{5.0, 6.0, 7.0},
                 pscalar3d{8.0}};

        // dual = complement (identity metric)
        CHECK(dual(s) == cmpl(s));
        CHECK(dual(v) == cmpl(v));
        CHECK(dual(B) == cmpl(B));
        CHECK(dual(ps) == cmpl(ps));
        CHECK(dual(M_e) == cmpl(M_e));
        CHECK(dual(M) == cmpl(M));

        fmt::println("dual = complement for Euclidean algebras with identity metric");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Extended Metric Tests - Verify correctness using stored metric matrices
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<3,0,0>: extended metric matrix validation")
    {
        fmt::println("G<3,0,0>: extended metric matrix validation");

        // Get the extended metric matrix from stored constants
        auto G = ega3d_metric_view();

        // Basis order: 1, e1, e2, e3, e23, e31, e12, e123
        // Expected values (all +1 for Euclidean with identity metric):
        // - Scalar (1):    G[0,0] = 1
        // - Vectors:       G[1,1] = G[2,2] = G[3,3] = 1
        // - Bivectors:     G[4,4] = G[5,5] = G[6,6] = 1
        // - Pseudoscalar:  G[7,7] = 1

        fmt::println("  Extended metric diagonal values:");
        fmt::println("    G[0,0] = {} (scalar 1)", G[0, 0]);
        fmt::println("    G[1,1] = {} (e1)", G[1, 1]);
        fmt::println("    G[2,2] = {} (e2)", G[2, 2]);
        fmt::println("    G[3,3] = {} (e3)", G[3, 3]);
        fmt::println("    G[4,4] = {} (e23)", G[4, 4]);
        fmt::println("    G[5,5] = {} (e31)", G[5, 5]);
        fmt::println("    G[6,6] = {} (e12)", G[6, 6]);
        fmt::println("    G[7,7] = {} (e123)", G[7, 7]);

        // Verify all diagonal elements are +1
        CHECK(G[0, 0] == 1); // scalar
        CHECK(G[1, 1] == 1); // e1
        CHECK(G[2, 2] == 1); // e2
        CHECK(G[3, 3] == 1); // e3
        CHECK(G[4, 4] == 1); // e23
        CHECK(G[5, 5] == 1); // e31
        CHECK(G[6, 6] == 1); // e12
        CHECK(G[7, 7] == 1); // e123

        // Verify off-diagonal elements are zero (orthogonal basis)
        for (size_t i = 0; i < 8; ++i) {
            for (size_t j = 0; j < 8; ++j) {
                if (i != j) {
                    CHECK(G[i, j] == 0);
                }
            }
        }

        fmt::println("All extended metric values correct for EGA3D");
        fmt::println("Conforming property: G(eᵢ ∧ eⱼ) = G(eᵢ) · G(eⱼ) = (+1)·(+1) = +1");
    }

    TEST_CASE("G<3,0,0>: extended metric recursive extraction via wedge products")
    {
        fmt::println("G<3,0,0>: extended metric recursive extraction via wedge products");

        auto G = ega3d_metric_view();

        // Level 0: Scalar (always 1)
        CHECK(G[0, 0] == 1);

        // Level 1: Vectors (extract from dot products)
        value_t g_e1 = value_t(dot(e1_3d, e1_3d));
        value_t g_e2 = value_t(dot(e2_3d, e2_3d));
        value_t g_e3 = value_t(dot(e3_3d, e3_3d));
        CHECK(abs(g_e1 - value_t(G[1, 1])) < eps);
        CHECK(abs(g_e2 - value_t(G[2, 2])) < eps);
        CHECK(abs(g_e3 - value_t(G[3, 3])) < eps);

        // Level 2: Bivectors (extract from wedge products + dot)
        // Basis order: e23, e31, e12
        auto e23_constructed = wdg(e2_3d, e3_3d);
        auto e31_constructed = wdg(e3_3d, e1_3d);
        auto e12_constructed = wdg(e1_3d, e2_3d);

        value_t g_e23 = value_t(dot(e23_constructed, e23_constructed));
        value_t g_e31 = value_t(dot(e31_constructed, e31_constructed));
        value_t g_e12 = value_t(dot(e12_constructed, e12_constructed));

        CHECK(abs(g_e23 - value_t(G[4, 4])) < eps);
        CHECK(abs(g_e31 - value_t(G[5, 5])) < eps);
        CHECK(abs(g_e12 - value_t(G[6, 6])) < eps);

        // Verify constructed bivectors match canonical ones
        CHECK(e23_constructed == e23_3d);
        CHECK(e31_constructed == e31_3d);
        CHECK(e12_constructed == e12_3d);

        // Level 3: Trivector/Pseudoscalar (extract from wedge products + dot)
        auto e123_constructed = wdg(wdg(e1_3d, e2_3d), e3_3d);
        value_t g_e123 = value_t(dot(e123_constructed, e123_constructed));
        CHECK(abs(g_e123 - value_t(G[7, 7])) < eps);
        CHECK(e123_constructed == I_3d);

        fmt::println("Recursive extraction: scalar → vectors → bivectors → trivector "
                     "all match");
    }

    TEST_CASE("G<3,0,0>: extended metric vs dot products")
    {
        fmt::println("G<3,0,0>: extended metric vs dot products");

        // The extended metric diagonal should match dot(basis, basis)
        auto G = ega3d_metric_view();

        // Vectors
        CHECK(abs(value_t(dot(e1_3d, e1_3d)) - value_t(G[1, 1])) < eps);
        CHECK(abs(value_t(dot(e2_3d, e2_3d)) - value_t(G[2, 2])) < eps);
        CHECK(abs(value_t(dot(e3_3d, e3_3d)) - value_t(G[3, 3])) < eps);

        // Bivectors
        CHECK(abs(value_t(dot(e23_3d, e23_3d)) - value_t(G[4, 4])) < eps);
        CHECK(abs(value_t(dot(e31_3d, e31_3d)) - value_t(G[5, 5])) < eps);
        CHECK(abs(value_t(dot(e12_3d, e12_3d)) - value_t(G[6, 6])) < eps);

        // Pseudoscalar
        CHECK(abs(value_t(dot(I_3d, I_3d)) - value_t(G[7, 7])) < eps);

        fmt::println("Dot products match extended metric diagonal");
    }

    TEST_CASE("G<3,0,0>: exponential function R = exp(bivec)")
    {
        fmt::println("G<3,0,0>: exponential function R = exp(bivec)");

        // rotation by 30° in the e12 plane (e1^e2)
        auto u = vec3d{1.0, 0.0, 0.0};
        auto v = vec3d(std::cos(pi / 6.0), std::sin(pi / 6.0), 0.0); // unit vec +30°
        auto r = nrm(u) * nrm(v);

        auto angle_uv = angle(u, v);
        auto angle_vu = angle(v, u);
        auto B = -0.5 * angle_uv * e12_3d; // B = -(phi/2)*e12_3d for rotation by phi

        auto uv = u * v;
        auto vu = v * u;

        auto R = get_rotor(e12_3d, deg2rad(30));
        auto R2 = exp(B);
        auto R3 = sqrt(vu);

        auto v2 = rotate(u, R);
        auto v3 = rotate(u, R2);
        auto v4 = rotate(u, R3);

        fmt::println("u                  = {}", u);
        fmt::println("v                  = {}", v);
        fmt::println("u*v                = {}", uv);
        fmt::println("v*u                = {}", vu);
        fmt::println("(u,v) = {}, (u,v)° = {}", angle_uv, rad2deg(angle_uv));
        fmt::println("(v,u) = {}, (v,u)° = {}", angle_vu, rad2deg(angle_vu));
        fmt::println("");
        fmt::println("B = -0.5*(u,v)*e12 = {}", B);
        fmt::println("R = {}, R*rev(R)   = {}", R, R * rev(R));
        fmt::println("R2 = exp(B)        = {}", R2);
        fmt::println("R3 = sqrt(v*u)     = {}", R3);
        fmt::println("");
        fmt::println("v2 (by rotor R)    = {}", v2);
        fmt::println("v3 (by rotor R2)   = {}", v3);
        fmt::println("v4 (by rotor R3)   = {}", v4);

        CHECK(r == 1.0);
        CHECK(v == v2); // +30° rotation confirmed
        CHECK(v == v3); // +30° rotation confirmed
        CHECK(v == v4); // +30° rotation confirmed

        // we have a rotor R
        CHECK(gr0(R * rev(R)) == 1.0);
        CHECK(nrm(gr2(R * rev(R))) <= eps);

        // we have a rotor R2
        CHECK(gr0(R2 * rev(R2)) == 1.0);
        CHECK(nrm(gr2(R2 * rev(R2))) <= eps);

        // we have a rotor R3
        CHECK(gr0(R3 * rev(R3)) == 1.0);
        CHECK(nrm(gr2(R3 * rev(R3))) <= eps);

        CHECK(v == rotate(u, sqrt(vu)));
        CHECK(u == rotate(v, sqrt(uv)));

        CHECK(exp(e12_3d * pi) == mvec3d_e{scalar3d{-1.0}});

        fmt::println("");
    }

    TEST_CASE("G<3,0,0>: exponential function R = exp(bivec) for 24 angles")
    {
        fmt::println("G<3,0,0>: exponential function R = exp(bivec) for 24 angles");

        // test exp() for 24 uniformly distributed angles in [-pi, pi)
        // rotation plane: e12_3d (e1^e2), reference vector: e1_3d
        int const n_angles = 24;
        for (int k = 0; k < n_angles; ++k) {

            value_t const phi = -pi + k * (2.0 * pi / n_angles); // step: pi/12 = 15°
            auto const B =
                -(phi / 2.0) * e12_3d; // B = -(phi/2)*e12_3d for rotation by phi

            auto const R = exp(B);
            auto const R_ref = get_rotor(e12_3d, phi);

            // exp(B) must be a valid rotor: R * rev(R) == 1 (scalar), 0 (bivector)
            CHECK(std::abs(to_val(gr0(R * rev(R))) - 1.0) <= eps);
            CHECK(nrm(gr2(R * rev(R))) <= eps);

            // exp(-(phi/2) * e12_3d) must match get_rotor(e12_3d, phi) component-wise
            CHECK(std::abs(to_val(gr0(R)) - to_val(gr0(R_ref))) <= eps);
            CHECK(nrm(gr2(R) - gr2(R_ref)) <= eps);

            // rotating e1 by R must give (cos(phi), sin(phi), 0)
            auto const v_rot = rotate(e1_3d, R);
            CHECK(std::abs(v_rot.x - std::cos(phi)) <= eps);
            CHECK(std::abs(v_rot.y - std::sin(phi)) <= eps);
            CHECK(std::abs(v_rot.z) <= eps);
        }

        fmt::println("");
    }

    TEST_CASE("G<3,0,0>: square root function S = sqrt(rotor) for 24 angles")
    {
        fmt::println("G<3,0,0>: square root function S = sqrt(rotor) for 24 angles");

        // test sqrt() for 24 uniformly distributed angles in [-pi, pi)
        // R = get_rotor(e12_3d, phi) = exp(-phi/2 * e12_3d) encodes rotation by phi.
        // sqrt() uses normalize(1 + M), which halves the half-angle:
        //   sqrt(R) = exp(-phi/4 * e12_3d) = get_rotor(e12_3d, phi/2)
        // => sqrt(R) encodes rotation by phi/2.
        int const n_angles = 24;
        for (int k = 0; k < n_angles; ++k) {

            value_t const phi = -pi + k * (2.0 * pi / n_angles); // step: pi/12 = 15°
            auto const R = get_rotor(e12_3d, phi);               // exp(-phi/2 * e12_3d)
            auto const S = sqrt(R);
            auto const S_ref = get_rotor(e12_3d, phi / 2.0); // exp(-phi/4 * e12_3d)

            // sqrt(R) must be a valid rotor: S * rev(S) == 1 (scalar), 0 (bivector)
            CHECK(std::abs(to_val(gr0(S * rev(S))) - 1.0) <= eps);
            CHECK(nrm(gr2(S * rev(S))) <= eps);

            // sqrt(exp(-phi/2 * e12_3d)) must equal exp(-phi/4 * e12_3d) component-wise
            CHECK(std::abs(to_val(gr0(S)) - to_val(gr0(S_ref))) <= eps);
            CHECK(nrm(gr2(S) - gr2(S_ref)) <= eps);

            // rotating e1 by S = exp(-phi/4 * e12_3d) must give (cos(phi/2), sin(phi/2),
            // 0)
            auto const v_rot = rotate(e1_3d, S);
            CHECK(std::abs(v_rot.x - std::cos(phi / 2.0)) <= eps);
            CHECK(std::abs(v_rot.y - std::sin(phi / 2.0)) <= eps);
            CHECK(std::abs(v_rot.z) <= eps);
        }

        fmt::println("");
    }

    TEST_CASE("G<3,0,0>: logarithm function B = log(rotor) for 24 angles")
    {
        fmt::println("G<3,0,0>: logarithm function B = log(rotor) for 24 angles");

        // log() is the inverse of exp(): log(exp(B)) == B (within the principal range
        // |B| < pi) and exp(log(R)) == R. A tilted rotation plane exercises the full
        // bivector (all three components), over 24 angles avoiding the +-pi endpoints.
        auto const B_hat = normalize(bivec3d{1.0, -2.0, 2.0}); // unit rotation plane
        int const n_angles = 24;
        for (int k = 0; k < n_angles; ++k) {

            value_t const phi = -pi + (k + 0.5) * (2.0 * pi / n_angles); // in (-pi, pi)
            auto const B = phi * B_hat;
            auto const R = exp(B);
            auto const B_back = log(R);

            // log recovers the generator bivector
            CHECK(nrm(B_back - B) <= eps);

            // exp(log(R)) reproduces R component-wise
            auto const R_back = exp(B_back);
            CHECK(std::abs(to_val(gr0(R_back)) - to_val(gr0(R))) <= eps);
            CHECK(nrm(gr2(R_back) - gr2(R)) <= eps);
        }

        fmt::println("");
    }

    TEST_CASE("MVec3d: geometric product replacements")
    {
        fmt::println("MVec3d: geometric product replacements");

        auto s = scalar3d(5.0);
        auto v = vec3d{1.0, 2.0, 1.0};
        auto B = bivec3d{-1.0, 2.0, 1.0};
        auto ps = pscalar3d(-5.0);

        auto s1 = scalar3d(-3.0);
        auto v1 = vec3d{-3.0, -1.5, 2.0};
        auto B1 = bivec3d{-3.0, -1.5, 2.0};
        auto ps1 = pscalar3d(2.0);

        auto mv0 = mvec3d{scalar3d{}, v, B, ps}; // contains no scalar part
        auto mv1 = mvec3d{s, v, B, ps};          // contains a scalar part

        fmt::println("   mv0  = {}", mv0);
        fmt::println("   mv1  = {}", mv1);

        // all following tests must only work for vector multivector combinations
        //
        // check: v * mv == (mv >> v) + wdg(v, mv)
        // check: mv * v == (v << mv) + wdg(mv, v)
        //
        REQUIRE_FALSE(ps1 * mv0 == (mv0 >> ps1) + wdg(ps1, mv0));
        REQUIRE_FALSE(ps1 * mv1 == (mv1 >> ps1) + wdg(ps1, mv1));
        //
        REQUIRE_FALSE(B1 * mv0 == (mv0 >> B1) + wdg(B1, mv0));
        REQUIRE_FALSE(B1 * mv1 == (mv1 >> B1) + wdg(B1, mv1));
        //
        CHECK(v1 * mv0 == (mv0 >> v1) + wdg(v1, mv0));
        CHECK(v1 * mv1 == (mv1 >> v1) + wdg(v1, mv1));
        //
        REQUIRE_FALSE(s1 * mv0 == (mv0 >> s1) + wdg(s1, mv0));
        REQUIRE_FALSE(s1 * mv1 == (mv1 >> s1) + wdg(s1, mv1));
        //
        // now with reversed arguments (only for vectors)
        CHECK(mv0 * v1 == (v1 << mv0) + wdg(mv0, v1));
        CHECK(mv1 * v1 == (v1 << mv1) + wdg(mv1, v1));
    }

    TEST_CASE("MVec3d: dualization - complement vs. pseudoscalar-multiplication")
    {
        fmt::println("MVec3d: dualization - complement vs. pseudoscalar-multiplication");
        fmt::println("");

        ////////////////////////////////////////////////////////////////////////////////
        // Two schools of thought for dualizing a blade / multivector A in ega3d:
        //
        //   (1) THIS LIBRARY (metric + complement):
        //       dual(A) = cmpl(G * A), where G is the extended metric. In ega the metric
        //       is the identity, so dual(A) == cmpl(A). The complement is defined via the
        //       outer (wedge) product:  A ^ cmpl(A) = I_3d. Equivalently the right dual
        //       satisfies  dual(A) = rev(A) * I_3d  (note the reversion of A).
        //
        //   (2) OTHER SCHOOL (multiplication by the pseudoscalar):
        //       to_dual(A)   = A * I_3d           (right-multiply by the pseudoscalar)
        //       from_dual(A) = A * inv(I_3d)      (divide by the pseudoscalar,
        //                                          i.e. multiply by its inverse)
        //       Here NO reversion is applied to A before multiplying.
        //
        // The only difference between the two forward operations is the reversion of A:
        //
        //       dual(A) = rev(A) * I_3d  vs.  to_dual(A) = A * I_3d
        //
        // For a pure grade-k blade  rev(A) = (-1)^(k(k-1)/2) * A, so the two results
        // differ exactly by the reversion sign of grade k:
        //
        //       grade k        0    1    2    3
        //       (-1)^(k(k-1)/2) +    +    -    -
        //
        // => scalar (k=0) and vector (k=1) results agree, while bivector (k=2) and
        //    pseudoscalar (k=3) results have opposite signs.
        //
        // FORWARD dualization (primal -> dual) for every basis element:
        //
        //   primal A   grade   ours: dual(A)=rev(A)*I_3d   theirs: A*I_3d   sign
        //   --------   -----   ------------------------    --------------   ----
        //      1         0              e123                    e123          +
        //      e1        1              e23                     e23           +
        //      e2        1              e31                     e31           +
        //      e3        1              e12                     e12           +
        //      e23       2              e1                     -e1            -
        //      e31       2              e2                     -e2            -
        //      e12       2              e3                     -e3            -
        //      e123      3              1                      -1             -
        ////////////////////////////////////////////////////////////////////////////////

        auto const s = scalar3d{3.0};
        auto const v = vec3d{1.0, 2.0, 3.0};
        auto const B = bivec3d{1.0, 2.0, 3.0};
        auto const ps = pscalar3d{4.0};

        // forward dualization: ours (complement) vs. theirs (A * I_3d)
        auto const s_dual_ours = dual(s);   // -> pscalar3d
        auto const s_dual_mul = s * I_3d;   // -> pscalar3d
        auto const v_dual_ours = dual(v);   // -> bivec3d
        auto const v_dual_mul = v * I_3d;   // -> bivec3d
        auto const B_dual_ours = dual(B);   // -> vec3d
        auto const B_dual_mul = B * I_3d;   // -> vec3d
        auto const ps_dual_ours = dual(ps); // -> scalar3d
        auto const ps_dual_mul = ps * I_3d; // -> scalar3d

        // The GA formatters use fmt::nested_formatter, so a width spec like "{:<17}"
        // would pad the nested number, not the whole "Scalar3d(...)" token. Stringify
        // first, then pad the string, so columns line up regardless of the type printed.
        auto const sf = [](auto const& x) { return fmt::format("{}", x); };

        // one row of the forward/backward comparison table (column widths fixed)
        auto const dual_row = [&sf](std::string_view grade, auto const& operand,
                                    auto const& ours, auto const& theirs,
                                    std::string_view relation) {
            fmt::println("   {:^5} | {:<16} | {:<16} | {:<19} | {}", grade, sf(operand),
                         sf(ours), sf(theirs), relation);
        };
        auto const rule = [] { fmt::println("   {:-<82}", ""); };

        fmt::println("   FORWARD dualization  (primal A  ->  dual):");
        rule();
        fmt::println("   {:^5} | {:<16} | {:<16} | {:<19} | {}", "grade", "primal A",
                     "ours: dual(A)", "theirs: A*I_3d", "relation");
        rule();
        dual_row("0", s, s_dual_ours, s_dual_mul, "equal");
        dual_row("1", v, v_dual_ours, v_dual_mul, "equal");
        dual_row("2", B, B_dual_ours, B_dual_mul, "opposite sign");
        dual_row("3", ps, ps_dual_ours, ps_dual_mul, "opposite sign");
        fmt::println("");

        // grade 0 and grade 1: both schools agree
        CHECK(s_dual_mul == s_dual_ours);
        CHECK(v_dual_mul == v_dual_ours);
        // grade 2 and grade 3: results differ by the reversion sign (opposite signs)
        CHECK(B_dual_mul == -B_dual_ours);
        CHECK(ps_dual_mul == -ps_dual_ours);

        // explicit forward results for every basis element (the FORWARD table above):
        // grade 0 (scalar) and grade 1 (vector) -> both schools identical
        CHECK(dual(one_3d) == I_3d);
        CHECK(one_3d * I_3d == I_3d);
        CHECK(dual(e1_3d) == e23_3d);
        CHECK(e1_3d * I_3d == e23_3d);
        CHECK(dual(e2_3d) == e31_3d);
        CHECK(e2_3d * I_3d == e31_3d);
        CHECK(dual(e3_3d) == e12_3d);
        CHECK(e3_3d * I_3d == e12_3d);
        // grade 2 (bivector) and grade 3 (pseudoscalar) -> theirs is negated
        CHECK(dual(e23_3d) == e1_3d);
        CHECK(e23_3d * I_3d == -e1_3d);
        CHECK(dual(e31_3d) == e2_3d);
        CHECK(e31_3d * I_3d == -e2_3d);
        CHECK(dual(e12_3d) == e3_3d);
        CHECK(e12_3d * I_3d == -e3_3d);
        CHECK(dual(I_3d) == one_3d);
        CHECK(I_3d * I_3d == -one_3d);

        ////////////////////////////////////////////////////////////////////////////////
        // same comparison for a full multivector. dual() maps grade k -> grade (3-k),
        // so the input-grade sign table above shows up as: the result components coming
        // from the input scalar (-> gr3) and input vector (-> gr2) agree, while those
        // coming from the input bivector (-> gr1) and input pseudoscalar (-> gr0) flip.
        ////////////////////////////////////////////////////////////////////////////////

        auto const M = mvec3d{3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
        auto const M_dual_ours = dual(M);    // complement based
        auto const M_dual_mul = M * I_3d_mv; // pseudoscalar multiplication based

        fmt::println("   full multivector  M = {}:", M);
        fmt::println("     {:<22} = {}", "ours : dual(M)", sf(M_dual_ours));
        fmt::println("     {:<22} = {}", "them : M * I_3d_mv", sf(M_dual_mul));
        fmt::println("     (gr0, gr1 flip sign; gr2, gr3 agree)");
        fmt::println("");

        // grade-wise: gr0 (from input ps) and gr1 (from input bivec) flip sign,
        //             gr2 (from input vec) and gr3 (from input scalar) are equal
        CHECK(gr0(M_dual_mul) == -gr0(M_dual_ours));
        CHECK(gr1(M_dual_mul) == -gr1(M_dual_ours));
        CHECK(gr2(M_dual_mul) == gr2(M_dual_ours));
        CHECK(gr3(M_dual_mul) == gr3(M_dual_ours));

        ////////////////////////////////////////////////////////////////////////////////
        // inverse dualization (from dual) and a structural consequence:
        //
        //   - our complement-based dual is an INVOLUTION in odd-dimensional ega3d:
        //         dual(dual(A)) == A
        //
        //   - right-multiplication by the pseudoscalar is NOT an involution, because
        //         I_3d * I_3d == -1   (so (A * I_3d) * I_3d == -A).
        //     Hence the other school must explicitly divide by I_3d to undualize:
        //         (A * I_3d) * inv(I_3d) == A,     with inv(I_3d) == -I_3d.
        //
        // BACKWARD dualization (dual -> primal) for every basis element:
        //   ours re-applies dual()  (involution);  theirs is A*inv(I_3d) == -(A*I_3d).
        //
        //   dual A    grade   ours: dual(A)   theirs: A*inv(I_3d)   sign
        //   ------    -----   -------------   -------------------   ----
        //      1        0          e123             -e123            -
        //      e1       1          e23              -e23             -
        //      e2       1          e31              -e31             -
        //      e3       1          e12              -e12             -
        //      e23      2          e1                e1              +
        //      e31      2          e2                e2              +
        //      e12      2          e3                e3              +
        //      e123     3          1                 1               +
        ////////////////////////////////////////////////////////////////////////////////

        fmt::println("   BACKWARD dualization  (dual A  ->  primal):");
        rule();
        fmt::println("   {:^5} | {:<16} | {:<16} | {:<19} | {}", "grade", "dual A",
                     "ours: dual(A)", "theirs: A*inv(I_3d)", "relation");
        rule();
        dual_row("0", s, dual(s), s * inv(I_3d), "opposite sign");
        dual_row("1", v, dual(v), v * inv(I_3d), "opposite sign");
        dual_row("2", B, dual(B), B * inv(I_3d), "equal");
        dual_row("3", ps, dual(ps), ps * inv(I_3d), "equal");
        fmt::println("");

        fmt::println("   pseudoscalar facts (why the two backward operations differ):");
        fmt::println("     {:<18} = {:<14} ({})", "I_3d * I_3d", sf(I_3d * I_3d),
                     "I_3d^2 = -1  =>  *I_3d is not an involution");
        fmt::println("     {:<18} = {:<14} ({})", "inv(I_3d)", sf(inv(I_3d)),
                     "== rev(I_3d) for the unit pseudoscalar");
        fmt::println("     {:<18} = {}", "rev(I_3d)", sf(rev(I_3d)));
        fmt::println("     {:<18} = {:<14} ({})", "I_3d * inv(I_3d)",
                     sf(I_3d * inv(I_3d)), "round-trip identity");
        fmt::println("");

        // check inv == rev
        CHECK(inv(I_3d) == rev(I_3d));

        // explicit backward results for every basis element (the BACKWARD table above):
        // grade 0 (scalar) and grade 1 (vector): theirs is negated relative to ours
        CHECK(dual(one_3d) == I_3d);
        CHECK(one_3d * inv(I_3d) == -I_3d);
        CHECK(dual(e1_3d) == e23_3d);
        CHECK(e1_3d * inv(I_3d) == -e23_3d);
        CHECK(dual(e2_3d) == e31_3d);
        CHECK(e2_3d * inv(I_3d) == -e31_3d);
        CHECK(dual(e3_3d) == e12_3d);
        CHECK(e3_3d * inv(I_3d) == -e12_3d);
        // grade 2 (bivector) and grade 3 (pseudoscalar): both schools identical
        CHECK(dual(e23_3d) == e1_3d);
        CHECK(e23_3d * inv(I_3d) == e1_3d);
        CHECK(dual(e31_3d) == e2_3d);
        CHECK(e31_3d * inv(I_3d) == e2_3d);
        CHECK(dual(e12_3d) == e3_3d);
        CHECK(e12_3d * inv(I_3d) == e3_3d);
        CHECK(dual(I_3d) == one_3d);
        CHECK(I_3d * inv(I_3d) == one_3d);

        // our dual undoes itself
        CHECK(dual(dual(s)) == s);
        CHECK(dual(dual(v)) == v);
        CHECK(dual(dual(B)) == B);
        CHECK(dual(dual(ps)) == ps);

        // multiplying by I_3d twice negates (I_3d^2 == -1)
        CHECK(I_3d * I_3d == scalar3d{-1.0});
        CHECK((v * I_3d) * I_3d == -v);
        CHECK((B * I_3d) * I_3d == -B);

        // the other school recovers the original by dividing by the pseudoscalar
        CHECK(inv(I_3d) == -I_3d);
        CHECK((s * I_3d) * inv(I_3d) == s);
        CHECK((v * I_3d) * inv(I_3d) == v);
        CHECK((B * I_3d) * inv(I_3d) == B);
        CHECK((ps * I_3d) * inv(I_3d) == ps);

        ////////////////////////////////////////////////////////////////////////////////
        // concrete round-trip examples, with the explicit intermediate (dual) values
        // spelled out for all four grades: s = 3, v = (1,2,3), B = (1,2,3), ps = 4.
        ////////////////////////////////////////////////////////////////////////////////

        auto const rt_row = [&sf](std::string_view input, std::string_view school,
                                  auto const& fwd, auto const& bwd) {
            fmt::println("   {:<16} | {:<6} | {:<16} | {}", input, school, sf(fwd),
                         sf(bwd));
        };
        fmt::println("   concrete round trips (forward, then backward recovers input):");
        fmt::println("   {:-<70}", "");
        fmt::println("   {:<16} | {:<6} | {:<16} | {}", "input", "school",
                     "forward (dual)", "backward (-> input)");
        fmt::println("   {:-<70}", "");
        rt_row(fmt::format("s={}", sf(s)), "ours", dual(s), dual(dual(s)));
        rt_row("", "them", s * I_3d, (s * I_3d) * inv(I_3d));
        rt_row(fmt::format("v={}", sf(v)), "ours", dual(v), dual(dual(v)));
        rt_row("", "them", v * I_3d, (v * I_3d) * inv(I_3d));
        rt_row(fmt::format("B={}", sf(B)), "ours", dual(B), dual(dual(B)));
        rt_row("", "them", B * I_3d, (B * I_3d) * inv(I_3d));
        rt_row(fmt::format("ps={}", sf(ps)), "ours", dual(ps), dual(dual(ps)));
        rt_row("", "them", ps * I_3d, (ps * I_3d) * inv(I_3d));
        fmt::println("");

        // scalar s=3: forward dual agrees, both schools round-trip to s
        CHECK(dual(s) == pscalar3d{3.0});               // ours, forward
        CHECK(s * I_3d == pscalar3d{3.0});              // theirs, forward (equal)
        CHECK(dual(dual(s)) == scalar3d{3.0});          // ours, backward -> s
        CHECK((s * I_3d) * inv(I_3d) == scalar3d{3.0}); // theirs, backward -> s

        // vector v=(1,2,3): forward dual agrees, both schools round-trip to v
        CHECK(dual(v) == bivec3d{1.0, 2.0, 3.0});              // ours, forward
        CHECK(v * I_3d == bivec3d{1.0, 2.0, 3.0});             // theirs, forward (equal)
        CHECK(dual(dual(v)) == vec3d{1.0, 2.0, 3.0});          // ours, backward -> v
        CHECK((v * I_3d) * inv(I_3d) == vec3d{1.0, 2.0, 3.0}); // theirs, backward -> v

        // bivector B=(1,2,3): forward duals have opposite sign, both round-trip to B
        CHECK(dual(B) == vec3d{1.0, 2.0, 3.0});         // ours, forward
        CHECK(B * I_3d == vec3d{-1.0, -2.0, -3.0});     // theirs, forward (negated)
        CHECK(dual(dual(B)) == bivec3d{1.0, 2.0, 3.0}); // ours, backward -> B
        CHECK((B * I_3d) * inv(I_3d) == bivec3d{1.0, 2.0, 3.0}); // theirs, backward -> B

        // pseudoscalar ps=4: forward duals have opposite sign, both round-trip to ps
        CHECK(dual(ps) == scalar3d{4.0});                 // ours, forward
        CHECK(ps * I_3d == scalar3d{-4.0});               // theirs, forward (negated)
        CHECK(dual(dual(ps)) == pscalar3d{4.0});          // ours, backward -> ps
        CHECK((ps * I_3d) * inv(I_3d) == pscalar3d{4.0}); // theirs, backward -> ps
    }


    TEST_CASE("ega3d: is_close and is_same_rotation")
    {
        fmt::println("ega3d: is_close and is_same_rotation");

        // is_close: equality within a RELATIVE tolerance, where operator== cannot help.
        // At coordinates of order 1e6 one ulp already exceeds the absolute eps that
        // operator== measures against, so values agreeing to every digit compare unequal.
        auto const big = vec3d{1.0e6, 2.0e6, 3.0e6};
        auto const big_1ulp = vec3d{std::nextafter(big.x, 1e9), big.y, big.z};
        CHECK(big != big_1ulp);
        CHECK(is_close(big, big_1ulp));

        // not merely permissive: a relative 1e-9 is still rejected at the default
        CHECK(!is_close(big, vec3d{big.x * (1.0 + 1e-9), big.y, big.z}));

        // and distinct from is_congruent, which allows an arbitrary scale factor
        CHECK(is_congruent(vec3d{1.0, 2.0, 3.0}, vec3d{2.0, 4.0, 6.0}));
        CHECK(!is_close(vec3d{1.0, 2.0, 3.0}, vec3d{2.0, 4.0, 6.0}));

        // is_same_rotation: rotors DOUBLE-COVER the rotations, so R and -R rotate
        // identically while comparing unequal component by component
        auto const R = exp(-e12_3d * 0.6);
        auto const R_neg = mvec3d_e{-R};
        auto const v = vec3d{1.0, 2.0, 3.0};

        CHECK(R != R_neg);                 // representations differ ...
        CHECK(is_same_rotation(R, R_neg)); // ... the rotation does not
        CHECK(rotate(v, R) == rotate(v, R_neg));

        // a 2*pi turn is the identity rotation but returns the NEGATED rotor
        auto const R_2pi = exp(-e12_3d * pi);
        CHECK(is_same_rotation(R_2pi, mvec3d_e{Scalar3d<double>(1.0)}));
        CHECK(rotate(v, R_2pi) == v);

        // it still separates rotations that really differ
        CHECK(!is_same_rotation(R, exp(-e12_3d * 0.7)));
        CHECK(!is_same_rotation(exp(-e12_3d * 0.6), exp(-e23_3d * 0.6)));

        // a scaled rotor is not the same rotation -- the sandwich scales by |R|^2
        CHECK(!is_same_rotation(R, mvec3d_e{2.0 * R}));

        fmt::println("");
    }


    TEST_CASE("MVec3d: is_close for the odd and the full multivector")
    {
        fmt::println("MVec3d: is_close for the odd and the full multivector");

        auto A = mvec3d{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        auto B = A;
        CHECK(is_close(A, B));
        B.c0 += 1.0e-6;
        CHECK(!is_close(A, B));

        auto U = mvec3d_u{1.0, 2.0, 3.0, 4.0};
        auto V = U;
        CHECK(is_close(U, V));
        V.c0 += 1.0e-6;
        CHECK(!is_close(U, V));

        // relative, not absolute: cf. operator==
        auto big = mvec3d{1.0e9, 2.0e9, 3.0e9, 4.0e9, 5.0e9, 6.0e9, 7.0e9, 8.0e9};
        auto near = big;
        near.c3 += 1.0e-4;
        CHECK(is_close(big, near));
        CHECK(!(big == near));
    }

    TEST_CASE("PScalar3d: exp of a pseudoscalar (duality rotation, NOT a rotor)")
    {
        fmt::println("PScalar3d: exp of a pseudoscalar (duality rotation, NOT a rotor)");

        auto const al = 0.63;
        auto E = exp(pscalar3d(al));

        // cos(alpha) + I sin(alpha), and nothing at grades 1 and 2
        CHECK(std::abs(value_t(gr0(E)) - std::cos(al)) < eps);
        CHECK(std::abs(value_t(gr3(E)) - std::sin(al)) < eps);
        CHECK(nrm(gr1(E)) < eps);
        CHECK(nrm(gr2(E)) < eps);

        // composes like a complex phase and has unit magnitude
        CHECK(is_close(exp(pscalar3d(al)) * exp(pscalar3d(-0.41)),
                       exp(pscalar3d(al - 0.41))));
        CHECK(std::abs(nrm(E) - 1.0) < eps);

        // I_3d is CENTRAL, so the sandwich returns its argument unchanged: this is a
        // duality rotation acting by multiplication, not a motion
        auto X = mvec3d{scalar3d(0.5), vec3d{1.0, 2.0, 3.0}, bivec3d{-2.0, 1.0, 0.5},
                        pscalar3d(1.3)};
        CHECK(is_close(E * X, X * E));
        CHECK(is_close(E * X * rev(E), X));

        // alpha = pi/2 is the Hodge map: e -> -b, b -> e
        auto ev = vec3d{0.7, -0.3, 0.5};
        auto bv = vec3d{0.2, 0.9, -0.4};
        auto f = mvec3d{ev} + mvec3d{cmpl(bv)};
        auto fd = exp(pscalar3d(0.5 * pi)) * f;
        CHECK(is_close(gr1(fd), -bv));
        CHECK(is_close(gr2(fd), cmpl(ev)));
    }

    TEST_CASE("MVec3d: exp, cos and sin of a general multivector")
    {
        fmt::println("MVec3d: exp, cos and sin of a general multivector");

        auto one = []() { return mvec3d{scalar3d(1.0)}; };
        auto series_exp = [&](mvec3d const& Z) {
            mvec3d acc = one(), term = one();
            for (int m = 1; m <= 80; ++m) {
                term = term * Z * (1.0 / m);
                acc = acc + term;
            }
            return acc;
        };
        auto series_trig = [&](mvec3d const& Z, bool sine) {
            mvec3d acc{}, p = one();
            double f = 1.0;
            for (int m = 0; m <= 80; ++m) {
                if (m > 0) {
                    p = p * Z;
                    f *= m;
                }
                bool const take = sine ? (m % 2 == 1) : (m % 2 == 0);
                if (take) {
                    int const k = sine ? (m - 1) / 2 : m / 2;
                    acc = acc + ((k % 2 ? -1.0 : 1.0) / f) * p;
                }
            }
            return acc;
        };

        auto gen = mvec3d{scalar3d(0.4), vec3d{1.1, -0.7, 0.9}, bivec3d{0.3, -1.2, 0.6},
                          pscalar3d(1.5)};
        // an element of the subalgebra A(n) = span{1, n, I n, I} (the plane-wave case)
        auto n = normalize(vec3d{1.0, 2.0, -0.5});
        auto An = mvec3d{scalar3d(0.4)} + mvec3d{-0.7 * n} + mvec3d{0.9 * cmpl(n)} +
                  mvec3d{pscalar3d(1.3)};
        // a null W: |v| = |B| and v ^ B = 0, so W^2 = 0 and the series terminate
        auto nul = mvec3d{e1_3d} + mvec3d{cmpl(e1_3d)};
        auto zero = mvec3d{};

        for (auto const& Z : {gen, An, nul, zero, 3.7 * gen}) {
            CHECK(is_close(exp(Z), series_exp(Z)));
            CHECK(is_close(cos(Z), series_trig(Z, false)));
            CHECK(is_close(sin(Z), series_trig(Z, true)));
        }

        // consistent with the single-grade overloads
        auto B = bivec3d{0.3, -1.2, 0.6};
        CHECK(is_close(exp(mvec3d{B}), mvec3d{exp(B)}));
        CHECK(is_close(exp(mvec3d{pscalar3d(1.5)}), exp(pscalar3d(1.5))));

        // identities
        CHECK(is_close(exp(zero), one()));
        CHECK(is_close(cos(gen) * cos(gen) + sin(gen) * sin(gen), one()));
        CHECK(is_close(exp(gen) * exp(-gen), one()));
        // Euler with the CENTRAL pseudoscalar
        CHECK(is_close(exp(mvec3d{I_3d} * gen), cos(gen) + mvec3d{I_3d} * sin(gen)));
        // addition theorems hold for COMMUTING arguments
        auto A2 = mvec3d{I_3d} * An;
        CHECK(is_close(An * A2, A2 * An));
        CHECK(is_close(exp(An + A2), exp(An) * exp(A2)));
        CHECK(is_close(cos(An + A2), cos(An) * cos(A2) - sin(An) * sin(A2)));
    }


    TEST_CASE("ega3d calculus: the four derivative types of nabla")
    {
        fmt::println("ega3d calculus: the four derivative types of nabla");

        auto r0 = vec3d{0.37, -0.62, 0.81};
        auto A = [](vec3d const& r) {
            return vec3d{r.y * r.z, r.x * r.x - r.z, std::sin(r.x) * r.y};
        };
        auto Bf = [&](vec3d const& r) { return cmpl(A(r)); };
        auto phi = [](vec3d const& r) {
            return scalar3d(std::sin(r.x) * r.y * r.y + r.z * r.z);
        };
        auto Tf = [](vec3d const& r) {
            return pscalar3d(std::cos(r.x) * r.y + r.z * r.z);
        };

        // classical references, computed independently
        auto curlA = [&](vec3d const& q) {
            vec3d c{0.0, 0.0, 0.0};
            double const h = 1.0e-5;
            for (int i = 0; i < 3; ++i) {
                auto e = i == 0 ? e1_3d : (i == 1 ? e2_3d : e3_3d);
                c = c + cross(e, (A(q + h * e) - A(q - h * e)) * (1.0 / (2 * h)));
            }
            return c;
        };
        auto divA = [&](vec3d const& q) {
            double d = 0.0;
            double const h = 1.0e-5;
            for (int i = 0; i < 3; ++i) {
                auto e = i == 0 ? e1_3d : (i == 1 ? e2_3d : e3_3d);
                d += value_t(dot(e, (A(q + h * e) - A(q - h * e)) * (1.0 / (2 * h))));
            }
            return d;
        };

        // scalar field: the outer derivative is the gradient, the inner one vanishes
        CHECK(nrm(nabla_wdg(phi, r0) - vec3d{std::cos(r0.x) * r0.y * r0.y,
                                             2 * std::sin(r0.x) * r0.y, 2 * r0.z}) <
              1.0e-6);
        CHECK(std::abs(value_t(nabla_dot(phi, r0))) < eps);

        // vector field: inner = div, outer = dual(curl)
        CHECK(std::abs(value_t(nabla_dot(A, r0)) - divA(r0)) < 1.0e-6);
        CHECK(nrm(nabla_wdg(A, r0) - dual(curlA(r0))) < 1.0e-6);

        // the relations here are pinned with dual(), the map that carries the metric and
        // so survives a port to sta4ds, where cmpl and dual differ. ega3d cannot tell
        // the two apart -- that coincidence is what this pins, so a change to either map
        // fails here rather than silently in another algebra
        CHECK(is_close(cmpl(curlA(r0)), dual(curlA(r0))));
        CHECK(is_close(cmpl(nabla_wdg(A, r0)), dual(nabla_wdg(A, r0))));

        // bivector field: inner = -curl of its dual, outer = div of its dual
        CHECK(nrm(nabla_dot(Bf, r0) - (-curlA(r0))) < 1.0e-6);
        CHECK(std::abs(value_t(nabla_wdg(Bf, r0)) - divA(r0)) < 1.0e-6);

        // trivector field: inner = dual(grad of its coeff), outer vanishes (no grade 4)
        auto gt = vec3d{-std::sin(r0.x) * r0.y, std::cos(r0.x), 2 * r0.z};
        CHECK(nrm(nabla_dot(Tf, r0) - dual(gt)) < 1.0e-6);
        CHECK(std::abs(value_t(nabla_wdg(Tf, r0))) < eps);

        // the Laplacian is a scalar operator: it preserves the grade of the field
        CHECK(nrm(gr0(laplacian(A, r0))) < 1.0e-6);
        CHECK(nrm(gr2(laplacian(A, r0))) < 1.0e-6);
    }

    TEST_CASE("ega3d calculus: nabla identities and rotation covariance")
    {
        fmt::println("ega3d calculus: nabla identities and rotation covariance");

        auto r0 = vec3d{0.37, -0.62, 0.81};
        auto A = [](vec3d const& r) {
            return vec3d{r.y * r.z, r.x * r.x - r.z, std::sin(r.x) * r.y};
        };
        auto Bf = [&](vec3d const& r) { return cmpl(A(r)); };

        // nabla ^ (nabla ^ M) = 0 and nabla . (nabla . M) = 0 -- these contain
        // curl(grad) = 0 and div(curl) = 0 as their classical cases
        auto outer = [&](vec3d const& q) { return nabla_wdg(A, q, 1.0e-3); };
        auto inner = [&](vec3d const& q) { return nabla_dot(Bf, q, 1.0e-3); };
        CHECK(std::abs(value_t(nabla_wdg(outer, r0, 1.0e-3))) < 1.0e-4);
        CHECK(std::abs(value_t(nabla_dot(inner, r0, 1.0e-3))) < 1.0e-4);

        // nabla . (nabla ^ X) + nabla ^ (nabla . X) = laplacian X: the vector derivative
        // is a square root of the Laplacian
        auto cA = [&](vec3d const& q) { return nabla_wdg(A, q, 1.0e-3); };
        auto dA = [&](vec3d const& q) { return nabla_dot(A, q, 1.0e-3); };
        CHECK(nrm(nabla_dot(cA, r0, 1.0e-3) + nabla_wdg(dA, r0, 1.0e-3) -
                  gr1(laplacian(A, r0, 1.0e-3))) < 1.0e-3);

        // nabla f' = (nabla f)': the operator is rotation covariant
        auto R = exp(wdg(e1_3d, e2_3d), 0.4) * exp(wdg(e2_3d, e3_3d), 0.25);
        auto fr = [&](vec3d const& q) { return rotate(A(rotate(q, rev(R))), R); };
        auto N = nabla(A, rotate(r0, rev(R)));
        CHECK(std::abs(value_t(gr0(nabla(fr, r0))) - value_t(gr0(N))) < 1.0e-6);
        CHECK(nrm(gr2(nabla(fr, r0)) - rotate(gr2(N), R)) < 1.0e-6);

        // book problems: nabla . r = 3, nabla ^ r = 0, laplacian(1/r) = 0, nabla x_1 = e1
        auto rf = [](vec3d const& q) { return q; };
        CHECK(std::abs(value_t(nabla_dot(rf, r0)) - 3.0) < 1.0e-8);
        CHECK(nrm(nabla_wdg(rf, r0)) < 1.0e-8);
        auto inv_r = [](vec3d const& q) { return scalar3d(1.0 / nrm(q)); };
        CHECK(std::abs(value_t(gr0(laplacian(inv_r, r0)))) < 1.0e-4);
        auto x1 = [](vec3d const& q) { return scalar3d(q.x); };
        CHECK(nrm(nabla_wdg(x1, r0) - e1_3d) < 1.0e-8);
    }

    TEST_CASE("ega3d calculus: stencil schemes deliver their advertised order")
    {
        fmt::println("ega3d calculus: stencil schemes deliver their advertised order");

        auto r0 = vec3d{0.37, -0.62, 0.81};
        auto A = [](vec3d const& r) {
            return vec3d{std::sin(r.x) * r.y, std::exp(0.3 * r.y) * r.z,
                         std::cos(r.z) * r.x};
        };
        auto div_exact = [&](vec3d const& r) {
            return std::cos(r.x) * r.y + 0.3 * std::exp(0.3 * r.y) * r.z -
                   std::sin(r.z) * r.x;
        };

        // the generator reports what it achieved -- at EVERY order, not only the low
        // ones. This used to stop at 6: order 8 and up came back reporting order 0,
        // which sent fd_step() to machine eps and made nabla wrong by ~45 % with no
        // error raised. Asking for more accuracy must never give less.
        for (int ord = 2; ord <= 12; ord += 2) {
            CAPTURE(ord);
            CHECK(central_scheme(1, ord).order == ord);
            CHECK(central_scheme(2, ord).order == ord);
        }
        // and the classical central weights come out
        auto d2 = central_scheme(2, 2);
        CHECK(std::abs(d2.weights[0] - 1.0) < eps);
        CHECK(std::abs(d2.weights[1] + 2.0) < eps);
        CHECK(std::abs(d2.weights[2] - 1.0) < eps);

        // halving h must reduce the error by 2^order -- this is what routing through
        // stencil_t buys over a hard-coded 3-point difference
        for (int ord : {2, 4, 6}) {
            auto sc = central_scheme(1, ord);
            double const h1 = 0.1, h2 = 0.05;
            double const e1 =
                std::abs(value_t(gr0(nabla(A, r0, sc, h1))) - div_exact(r0));
            double const e2 =
                std::abs(value_t(gr0(nabla(A, r0, sc, h2))) - div_exact(r0));
            double const observed = std::log2(e1 / e2);
            CHECK(observed > ord - 0.3);
            CHECK(observed < ord + 0.3);
        }

        // the default step tracks the scheme's order, so a higher order is more accurate
        // at its own default step too
        double const err2 = std::abs(value_t(gr0(nabla(A, r0))) - div_exact(r0));
        double const err6 =
            std::abs(value_t(gr0(nabla(A, r0, central_scheme(1, 6)))) - div_exact(r0));
        CHECK(err6 < err2);

        // ... and that must still hold past order 6, at the DEFAULT step: this is the
        // end-to-end form of the regression above, where the divergence came back as
        // 0.83 instead of 1.50 because the step had collapsed
        for (int ord : {8, 10}) {
            CAPTURE(ord);
            double const err = std::abs(
                value_t(gr0(nabla(A, r0, central_scheme(1, ord)))) - div_exact(r0));
            CHECK(err < err2);
            CHECK(err < 1.0e-9);
        }

        // a staggered scheme samples on the half-offsets and develops between them,
        // so it never evaluates the field at r0 itself -- and still gets the answer
        double const err_stag =
            std::abs(value_t(gr0(nabla(A, r0, staggered_scheme(1, 4)))) - div_exact(r0));
        CHECK(err_stag < 1.0e-9);
    }

    TEST_CASE("ega3d calculus: the single Maxwell equation on a plane wave")
    {
        fmt::println("ega3d calculus: the single Maxwell equation on a plane wave");

        // a vacuum plane wave as the electromagnetic multivector f = e + I b, with
        // e along e1, b along e2 and propagation along e3 (natural units, sqrt(eps mu) =
        // 1)
        double const w = 1.3, E0 = 0.8;
        auto wave = [&](vec3d const& q, double t) {
            double const ph = w * (q.z - t);
            return mvec3d(E0 * std::cos(ph) * e1_3d) +
                   mvec3d(cmpl(E0 * std::cos(ph) * e2_3d));
        };
        auto r0 = vec3d{0.37, -0.62, 0.81};
        double const t0 = 0.29;

        // D f = (nabla + sqrt(eps mu) d/dt) f = 0
        auto Df = nabla(wave, r0, t0, 1.0e-5) + d_dt(wave, r0, t0, 1.0e-5);
        CHECK(nrm(Df) < 1.0e-7);

        // the wave equation follows: box f = laplacian f - eps mu d2f/dt2 = 0
        auto lap = laplacian([&](vec3d const& q) { return wave(q, t0); }, r0, 1.0e-3);
        auto d2t = d_dt(wave, r0, t0, central_scheme(2, 2), 1.0e-3);
        CHECK(nrm(lap - d2t) < 1.0e-4);
    }

    TEST_CASE("ega3d: Griffiths -- the Gibbs dictionary in one algebra")
    {
        fmt::println("ega3d: Griffiths -- the Gibbs dictionary in one algebra");

        // The translation table between classical vector algebra (Griffiths,
        // Introduction to Electrodynamics 4th ed., ch. 1) and this library. Every row
        // is checked against a reference written out by hand, not against another
        // library call.
        vec3d const a{1.3, -0.7, 2.1}, b{-0.4, 1.9, 0.6}, c{0.8, 0.5, -1.7};
        auto cross_ref = [](vec3d const& u, vec3d const& v) {
            return vec3d{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z,
                         u.x * v.y - u.y * v.x};
        };

        // (1.4)/(1.12)/(1.14): the cross product is the DUAL of the wedge. Both
        // spellings exist; wdg is the primitive and carries the plane, cross the axis
        CHECK(nrm(cross(a, b) - cross_ref(a, b)) < eps);
        CHECK(nrm(dual(wdg(a, b)) - cross_ref(a, b)) < eps);
        CHECK(dual(wdg(e1_3d, e2_3d)) == e3_3d);
        // the magnitude |a x b| = AB sin(theta) is the AREA -- no dual involved
        CHECK(nrm(cross(a, b)) == doctest::Approx(value_t(nrm(wdg(a, b)))));

        // (1.15)/(1.16): the scalar triple product is the trivector, cyclic by the
        // antisymmetry of wdg rather than by a rule to remember
        CHECK(value_t(wdg(a, wdg(b, c))) ==
              doctest::Approx(value_t(dot(a, cross_ref(b, c)))));
        CHECK(value_t(wdg(a, wdg(b, c))) == doctest::Approx(value_t(wdg(b, wdg(c, a)))));

        // (1.17) BAC-CAB is not a mnemonic but a vector contracted into a blade
        CHECK(nrm(cross_ref(a, cross_ref(b, c)) - (a << wdg(b, c))) < 1.0e-14);
        CHECK(nrm(cross_ref(a, cross_ref(b, c)) + (wdg(b, c) >> a)) < 1.0e-14);
        // and it agrees with the book's right-hand side
        CHECK(nrm((a << wdg(b, c)) - (b * value_t(dot(a, c)) - c * value_t(dot(a, b)))) <
              1.0e-14);

        // (11.56): the double cross product that runs through the radiation chapters
        // is a REJECTION -- the part of X perpendicular to the line of sight
        vec3d const n = normalize(c);
        vec3d const X{0.3, 0.9, -0.2};
        CHECK(nrm(cross_ref(n, cross_ref(n, X)) + reject_from(X, n)) < 1.0e-14);

        // (3.104)/(5.89): the field of a point dipole. Griffiths sets both as
        // problems (3.36, 5.34) because the coordinate-free form is the non-obvious
        // one. The GA form is a single REFLECTION in the separation direction, and
        // this pins the identity that rewrites it, at vector grade:
        //
        //     3 (p.rhat) rhat - p  ==  1/2 (p + 3 rhat (x) p (x) rhat)
        //
        // together with the split reading that says where the "3" went: writing
        // p = p_parallel + p_perp, the bracket is just 2 p_parallel - p_perp -- twice
        // the axial part minus the equatorial one. Gibbs notation has no name for
        // p_perp, which is why it has to write 3(p.rhat)rhat and subtract all of p.
        vec3d const p{0.7, -1.3, 0.4};
        vec3d const gibbs = 3.0 * value_t(dot(p, n)) * n - p;
        CHECK(nrm(0.5 * (p + 3.0 * gr1(n * p * n)) - gibbs) < 1.0e-14);
        CHECK(nrm(2.0 * project_onto(p, n) - reject_from(p, n) - gibbs) < 1.0e-14);
        // the sandwich itself is the axis reflection, 2(p.rhat)rhat - p
        CHECK(nrm(gr1(n * p * n) - (2.0 * value_t(dot(p, n)) * n - p)) < 1.0e-14);
        CHECK(nrm(gr1(n * p * n) - reflect_on_vec(p, n)) < eps);

        fmt::println("Gibbs dictionary: cross = dual(wdg), triple product = trivector, "
                     "BAC-CAB = contraction, double cross = rejection");
    }

    TEST_CASE("ega3d calculus: Maxwell WITH sources, and the field invariants")
    {
        fmt::println("ega3d calculus: Maxwell WITH sources, and the field invariants");

        // The plane-wave case above solves the source-free equation. Here the sources
        // are switched on, and the FOUR GRADES of one product are Griffiths' four
        // Maxwell equations (7.40). Natural units: c = eps0 = mu0 = 1.
        //
        //   gr0(nabla f) = div E          = rho / eps0        (2.14) Gauss
        //   gr1(nabla f) = -c curl B      = -mu0 c J          (5.56) Ampere
        //   gr2(nabla f) = I curl E       = 0     (static)    (7.16) Faraday
        //   gr3(nabla f) = I c div B      = 0                 (5.50) no monopoles
        //
        // fields with hand-computable sources: E = (rho/3) r has div E = rho and no
        // curl; B = (1/2) J x r has curl B = J and no divergence
        double const rho = 0.9;
        vec3d const J{0.0, 0.0, 1.4};
        auto cross_ref = [](vec3d const& u, vec3d const& v) {
            return vec3d{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z,
                         u.x * v.y - u.y * v.x};
        };
        auto f = [&](vec3d const& q) {
            vec3d const E = (rho / 3.0) * q;
            vec3d const B = 0.5 * cross_ref(J, q);
            return mvec3d(E) + mvec3d(dual(B)); // f = e + I b
        };
        vec3d const r0{0.6, -0.3, 0.8};
        auto const N = nabla(f, r0, 1.0e-4);
        CHECK(value_t(gr0(N)) == doctest::Approx(rho)); // Gauss
        CHECK(nrm(gr1(N) + J) < 1.0e-8);                // Ampere
        CHECK(nrm(gr2(N)) < 1.0e-8);                    // Faraday (static)
        CHECK(std::abs(value_t(gr3(N))) < 1.0e-8);      // no magnetic monopoles

        // falsify: a field with the wrong 1/3 factor reports a different charge density
        auto fbad = [&](vec3d const& q) {
            return mvec3d(0.5 * rho * q) + mvec3d(dual(0.5 * cross_ref(J, q)));
        };
        CHECK(value_t(gr0(nabla(fbad, r0, 1.0e-4))) != doctest::Approx(rho));

        // ---- the two field invariants and the two energy-momentum scalars ----------
        // Griffiths Prob. 12.47 asks the reader to show that (E.B) and (E^2 - c^2 B^2)
        // are relativistically invariant; they are the two grade parts of f * f
        vec3d const E{0.3, -1.2, 0.8}, B{1.1, 0.4, -0.6};
        auto const F = mvec3d(E) + mvec3d(dual(B));
        auto const FF = F * F;
        CHECK(value_t(gr0(FF)) ==
              doctest::Approx(value_t(nrm_sq(E)) - value_t(nrm_sq(B))));
        CHECK(value_t(gr3(FF)) == doctest::Approx(2.0 * value_t(dot(E, B))));
        // and f * rev(f) carries the energy density (8.5) and the Poynting vector
        // (8.10) -- the two halves of Poynting's theorem, in one product
        auto const FR = F * rev(F);
        CHECK(value_t(gr0(FR)) ==
              doctest::Approx(value_t(nrm_sq(E)) + value_t(nrm_sq(B))));
        CHECK(nrm(gr1(FR) - 2.0 * cross_ref(E, B)) < 1.0e-13);

        fmt::println("nabla f = rho - J: four grades, four Maxwell equations; f*f and "
                     "f*rev(f) carry the invariants and the energy/Poynting pair");
    }

    TEST_CASE("ega3d calculus: the pointwise nabla cannot see a delta function")
    {
        fmt::println("ega3d calculus: the pointwise nabla cannot see a delta function");

        // Griffiths sect. 1.5.1 opens with the paradox that motivates the Dirac delta:
        // E = rhat/r^2 has div E = 0 EVERYWHERE it can be evaluated, yet its flux
        // through any enclosing surface is 4 pi. The resolution is (1.99),
        // div(rhat/r^2) = 4 pi delta^3(r).
        //
        // The library has no delta and needs none -- a distribution is not an element of
        // any geometric algebra, and both reference books use delta symbolically. But
        // the consequence for THIS layer is a live trap and is pinned here: a field's
        // singular sources are invisible to the pointwise operator, so a numerical Gauss
        // law must be evaluated as a FLUX, never as a volume integral of nabla_dot.
        auto E = [](vec3d const& r) { return r / std::pow(value_t(nrm(r)), 3.0); };

        // (a) the pointwise divergence, summed over a ball that CONTAINS the source
        value_t const R = 1.0;
        int const N = 20;
        value_t const hc = 2.0 * R / N;
        value_t vol_sum = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    vec3d const c{-R + (i + 0.5) * hc, -R + (j + 0.5) * hc,
                                  -R + (k + 0.5) * hc};
                    if (value_t(nrm(c)) > R) continue;
                    vol_sum += value_t(nabla_dot(E, c, 1.0e-4)) * hc * hc * hc;
                }
            }
        }
        CHECK(std::abs(vol_sum) < 1.0e-2); // ... reports essentially nothing

        // (b) the flux of the same field through the bounding sphere
        int const Nt = 100, Np = 200;
        value_t flux = 0.0;
        for (int a = 0; a < Nt; ++a) {
            for (int b = 0; b < Np; ++b) {
                value_t const th = (a + 0.5) * pi / Nt, ph = (b + 0.5) * 2.0 * pi / Np;
                vec3d const n{std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph),
                              std::cos(th)};
                flux += value_t(dot(E(R * n), n)) * R * R * std::sin(th) * (pi / Nt) *
                        (2.0 * pi / Np);
            }
        }
        CHECK(flux == doctest::Approx(4.0 * pi).epsilon(1.0e-4)); // ... the whole charge

        // (c) the entire source sits in the difference -- that difference IS the delta
        CHECK(flux - vol_sum == doctest::Approx(4.0 * pi).epsilon(1.0e-3));

        // (d) and sampling AT the singularity returns a finite-difference artefact of
        // order 1/h^2, not a delta: ~3e12 for h = 1e-4. Never read it as a value.
        CHECK(std::abs(value_t(nabla_dot(E, vec3d{0.0, 0.0, 0.0}, 1.0e-4))) > 1.0e10);

        fmt::println("delta trap: div = 0 pointwise, flux = 4 pi; a numerical Gauss law "
                     "must be a flux, never a volume integral of nabla_dot");
    }

    TEST_CASE("ega3d calculus: Jancewicz P28 product rules and P29")
    {
        fmt::println("ega3d calculus: Jancewicz P28 product rules and P29");

        // Jancewicz Problems 0.28 (i)-(vi) and 0.29 (pp. 66-67) -- the product rules for
        // nabla acting on products of a vector field with a vector or a bivector field.
        // A, B are vector fields; Ab, Bb are bivector fields ("hatted" in the book):
        //
        //   (i)   grad.(A.Bb) = (grad^A).Bb - A.(grad.Bb)
        //   (ii)  grad^(A.B)  = (A.grad)B + (B.grad)A - B.(grad^A) - A.(grad^B)
        //   (iii) grad.(A^B)  = B(grad.A) - A(grad.B) + (A.grad)B - (B.grad)A
        //   (iv)  grad^(Ab Bb) = (Ab^grad).Bb + (Bb^grad).Ab - Bb.(grad.Ab) -
        //   Ab.(grad.Bb) (v)   grad^(A.Bb) = Bb(grad.A) - A(grad^Bb) + (A.grad)Bb -
        //   (Bb^grad)A (vi)  grad.(A^Bb) = (A.grad)Bb + (Bb^grad)A - Bb (x) (grad^A) -
        //   A^(grad.Bb) (29)  grad(A Bb)  = (grad A)Bb - A(grad Bb) + 2(A.grad)Bb
        //
        // The book's "." is gr_{|r-s|}(X * Y), so the terms are built as graded parts of
        // a geometric product rather than through the library contractions -- those
        // differ from it by a side swap at grade 1x2 and by a reversion sign at grades
        // 2 and 3.
        //
        // (vi)'s third term, written "(x)" above, is the MINGLED product (a wedge with a
        // dot over it in the book, eq. 31) = cmt. Both operands are bivectors there, so a
        // plain wedge would be grade 4 -- identically zero in 3d -- and a dot would be a
        // scalar, while the term has to be a bivector to match the left side.
        //
        // (iv) holds on the GRADE-1 part only. Its left side differentiates a cliffor
        // carrying grade 0 and grade 2, so the result has a trivector part as well
        // (measured |gr3| = 0.98 here, not a rounding artefact) which no term on the
        // right accounts for. Assert the vector part, not the whole.

        using Mv = MVec3d<value_t>;

        auto A = [](vec3d const& r) {
            return vec3d{std::sin(r.x) * r.y, r.x * r.z * r.z,
                         std::exp(0.3 * r.z) * r.x + r.y};
        };
        auto B = [](vec3d const& r) {
            return vec3d{std::cos(r.x) + r.z, r.x * r.y, r.y * r.y - r.z};
        };
        auto Bb = [](vec3d const& r) {
            return bivec3d{r.x * r.y + r.z, std::cos(r.z) * r.x, r.y * r.z + r.x * r.x};
        };
        auto Ab = [](vec3d const& r) {
            return bivec3d{std::sin(r.x), r.y * r.z, r.x + 2.0 * r.z};
        };
        vec3d const r{0.3, -0.4, 0.55};

        auto pd = [](auto&& f, vec3d const& p, int i, value_t h) {
            return detail::ega3d_axis_deriv(f, p, i, detail::ega3d_default_scheme(1), h);
        };
        auto cmp = [](vec3d const& v, int i) {
            return (i == 0) ? v.x : ((i == 1) ? v.y : v.z);
        };
        auto e = [](int i) { return detail::ega3d_basis<value_t>(i); };

        // |lhs - rhs| for all seven identities at step size h
        auto residuals = [&](value_t h) {
            auto Ar = A(r);
            auto Br = B(r);
            auto Bbr = Bb(r);
            auto Abr = Ab(r);
            auto curlA = nabla_wdg(A, r, h);
            auto divA = value_t(nabla_dot(A, r, h));
            auto curlB = nabla_wdg(B, r, h);
            auto divB = value_t(nabla_dot(B, r, h));
            auto divBb = nabla_dot(Bb, r, h); // vector
            auto wdgBb = nabla_wdg(Bb, r, h); // pseudoscalar
            auto divAb = nabla_dot(Ab, r, h);

            vec3d t1{}, t2{}, u1{}, u2{};
            bivec3d s1{};
            Mv s2{};
            for (int i = 0; i < 3; ++i) {
                t1 = t1 + cmp(Ar, i) * gr1(pd(B, r, i, h));             // (A.grad)B
                t2 = t2 + cmp(Br, i) * gr1(pd(A, r, i, h));             // (B.grad)A
                s1 = s1 + cmp(Ar, i) * gr2(pd(Bb, r, i, h));            // (A.grad)Bb
                s2 = s2 + Mv(wdg(Bbr, e(i))) * Mv(gr1(pd(A, r, i, h))); // (Bb^grad)A
                u1 = u1 +
                     gr1(Mv(wdg(Abr, e(i))) * Mv(gr2(pd(Bb, r, i, h)))); // (Ab^grad).Bb
                u2 = u2 +
                     gr1(Mv(wdg(Bbr, e(i))) * Mv(gr2(pd(Ab, r, i, h)))); // (Bb^grad).Ab
            }
            bivec3d const s2b = gr2(s2);

            auto AdotBb = [&](vec3d const& p) { return gr1(Mv(A(p)) * Mv(Bb(p))); };

            std::vector<value_t> d;
            d.push_back(
                std::abs(value_t(nabla_dot(AdotBb, r, h)) -
                         (value_t(gr0(Mv(curlA) * Mv(Bbr))) - value_t(dot(Ar, divBb)))));
            d.push_back(nrm(
                nabla_wdg(
                    [&](vec3d const& p) { return scalar3d(value_t(dot(A(p), B(p)))); }, r,
                    h) -
                (t1 + t2 - gr1(Mv(Br) * Mv(curlA)) - gr1(Mv(Ar) * Mv(curlB)))));
            d.push_back(
                nrm(nabla_dot([&](vec3d const& p) { return wdg(A(p), B(p)); }, r, h) -
                    (divA * Br - divB * Ar + t1 - t2)));

            Mv lhs_iv{};
            for (int i = 0; i < 3; ++i)
                lhs_iv =
                    lhs_iv + wdg(Mv(e(i)),
                                 pd([&](vec3d const& p) { return Mv(Ab(p)) * Mv(Bb(p)); },
                                    r, i, h));
            d.push_back(nrm(gr1(lhs_iv) - (u1 + u2 - gr1(Mv(Bbr) * Mv(divAb)) -
                                           gr1(Mv(Abr) * Mv(divBb)))));

            d.push_back(nrm(nabla_wdg(AdotBb, r, h) -
                            (divA * Bbr - gr2(Mv(Ar) * Mv(wdgBb)) + s1 - s2b)));
            d.push_back(
                nrm(nabla_dot([&](vec3d const& p) { return wdg(A(p), Bb(p)); }, r, h) -
                    (s1 + s2b - cmt(Bbr, curlA) - wdg(Ar, divBb))));

            Mv lhs29{}, gA{}, gB{};
            for (int i = 0; i < 3; ++i) {
                lhs29 =
                    lhs29 +
                    Mv(e(i)) *
                        pd([&](vec3d const& p) { return Mv(A(p)) * Mv(Bb(p)); }, r, i, h);
                gA = gA + Mv(e(i)) * pd(A, r, i, h);
                gB = gB + Mv(e(i)) * pd(Bb, r, i, h);
            }
            d.push_back(nrm(lhs29 - (gA * Mv(Bbr) - Mv(Ar) * gB + 2.0 * Mv(s1))));
            return d;
        };

        auto const d1 = residuals(1.0e-3);
        auto const d2 = residuals(5.0e-4);

        for (std::size_t i = 0; i < d1.size(); ++i) {
            CHECK(d1[i] < 1.0e-5);

            // The identities are EXACT, so the residual is nothing but the 2nd-order
            // scheme's truncation error and halving h must QUARTER it. Without this the
            // threshold above is nearly vacuous: at h = 1e-3 the truncation error IS
            // ~1e-6, so scaling one term of (vi) by 1.00002 -- a relative error of 2e-5,
            // a real mistake -- leaves a residual of 1.96e-6 that sails through the
            // threshold, while this ratio collapses from 4.00 to 1.04 and catches it.
            //
            // Measured 4.00 for all seven, so the band is tight. Its detection floor is
            // around 1e-6 relative error (a 2e-6 perturbation of (vi) gives 3.74); below
            // that the perturbation is lost in the truncation error and only a smaller h
            // would separate them.
            CHECK(d1[i] / d2[i] > 3.8);
            CHECK(d1[i] / d2[i] < 4.2);
        }
    }

    TEST_CASE("ega3d calculus: Jancewicz P34-P40 -- nabla on the standard fields")
    {
        fmt::println("ega3d calculus: Jancewicz P34-P40 -- nabla on the standard fields");

        // Jancewicz Problems 0.34-0.40 (p. 67-68): the standard closed-form derivatives,
        // with a a constant vector, Bc a constant bivector and r the position field.
        // These are the reference values a reader checks a calculus layer against, and
        // they reach parts the product-rule case does not: the radial fields r^k and 1/r
        // (no polynomial, so the stencil's truncation error is live), the Laplacian, and
        // the FULL Clifford derivative nabla(f) rather than only its graded parts.
        //
        // The cases split by whether the stencil is exact for the field:
        //
        //  - a 2nd-order central difference is EXACT for polynomials of degree <= 2, so
        //    every linear/bilinear field below is reproduced to roundoff (~1e-15) and is
        //    asserted tightly. That is a property of the scheme, not a coincidence.
        //  - the radial fields are not polynomial, so their residual IS the truncation
        //    error and only its h^2 decay says the closed form is right. Same reasoning
        //    as the P28 case.

        using Mv = MVec3d<value_t>;
        vec3d const r{0.42, -0.31, 0.57};
        vec3d const a{0.7, 1.3, -0.5};
        bivec3d const Bc{0.6, -1.1, 0.35};
        value_t const h = 1.0e-3;
        auto rad = [](vec3d const& p) { return std::sqrt(value_t(dot(p, p))); };
        value_t const R = rad(r);

        // ---- exact for the scheme: polynomial fields, asserted to roundoff ----
        value_t const tight = 1.0e-10;

        // P34: grad(r.a) = a ; P40: grad x_j = e_j ; (a.grad)r = a is the same statement
        CHECK(nrm(nabla_wdg([&](vec3d const& p) { return scalar3d(value_t(dot(p, a))); },
                            r, h) -
                  a) < tight);
        CHECK(nrm(nabla_wdg([](vec3d const& p) { return scalar3d(p.y); }, r, h) -
                  detail::ega3d_basis<value_t>(1)) < tight);

        // P34: div(a^r) = -2a and curl(a^r) = 0 -- note div is NOT zero although a is
        // constant: the r factor carries the divergence
        auto awr = [&](vec3d const& p) { return wdg(a, p); };
        CHECK(nrm(nabla_dot(awr, r, h) - (-2.0) * a) < tight);
        CHECK(std::abs(value_t(nabla_wdg(awr, r, h))) < tight);

        // P34 with a constant BIVECTOR: div(Bc^r) = Bc and curl(r.Bc) = 2 Bc, div(r.Bc) =
        // 0. The factors 1 and 2 follow from the contraction convention; they are exact
        // rationals here, which is why they can be pinned at roundoff.
        CHECK(
            nrm(gr2(nabla_dot([&](vec3d const& p) { return wdg(Mv(Bc), Mv(p)); }, r, h)) -
                Bc) < tight);
        auto rdB = [&](vec3d const& p) { return gr1(Mv(p) * Mv(Bc)); };
        CHECK(nrm(nabla_wdg(rdB, r, h) - 2.0 * Bc) < tight);
        CHECK(std::abs(value_t(nabla_dot(rdB, r, h))) < tight);

        // P35: the full Clifford derivatives. grad r = 3 is the divergence of the
        // position field; grad(ar) and grad(ra) differ because a and r do not commute.
        CHECK(nrm(nabla([](vec3d const& p) { return p; }, r, h) - Mv(scalar3d(3.0))) <
              tight);
        CHECK(nrm(nabla([&](vec3d const& p) { return Mv(a) * Mv(p); }, r, h) - Mv(-a)) <
              tight);
        CHECK(nrm(nabla([&](vec3d const& p) { return Mv(p) * Mv(a); }, r, h) -
                  Mv(3.0 * a)) < tight);
        CHECK(nrm(nabla(awr, r, h) - Mv((-2.0) * a)) < tight);

        // P36: grad r^k = k r^(k-2) r -- exact for k = 2, truncation-limited otherwise
        CHECK(nrm(nabla_wdg([&](vec3d const& p) { return scalar3d(rad(p) * rad(p)); }, r,
                            h) -
                  2.0 * r) < tight);

        // ---- truncation-limited: radial fields, pinned by h^2 decay ----
        auto residuals = [&](value_t hh) {
            std::vector<value_t> d;
            // P34: grad((r.a)/r^3)
            d.push_back(nrm(
                nabla_wdg(
                    [&](vec3d const& p) {
                        return scalar3d(value_t(dot(p, a)) / std::pow(rad(p), 3));
                    },
                    r, hh) -
                (a / std::pow(R, 3) - 3.0 * value_t(dot(r, a)) / std::pow(R, 5) * r)));
            // P36: k = 3 and k = -1
            for (int k : {3, -1})
                d.push_back(nrm(nabla_wdg(
                                    [&](vec3d const& p) {
                                        return scalar3d(std::pow(rad(p), value_t(k)));
                                    },
                                    r, hh) -
                                value_t(k) * std::pow(R, value_t(k - 2)) * r));
            // P39: the Laplacian of the fundamental solution vanishes off the source
            d.push_back(std::abs(value_t(gr0(laplacian(
                [&](vec3d const& p) { return scalar3d(1.0 / rad(p)); }, r, hh)))));
            return d;
        };
        auto const d1 = residuals(1.0e-3);
        auto const d2 = residuals(5.0e-4);
        for (std::size_t i = 0; i < d1.size(); ++i) {
            CHECK(d1[i] < 1.0e-4);
            CHECK(d1[i] / d2[i] > 3.8);
            CHECK(d1[i] / d2[i] < 4.2);
        }
    }


    TEST_CASE("ega3d calculus: compact (Pade) schemes and grid derivatives")
    {
        fmt::println("ega3d calculus: compact (Pade) schemes and grid derivatives");

        // the compact schemes on i-1, i, i+1 are the textbook Pade forms:
        //   f'_{i-1} + 4 f'_i + f'_{i+1}      = 3 (f_{i+1} - f_{i-1}) / h
        //   f''_{i-1} + 10 f''_i + f''_{i+1}  = 12 (f_{i-1} - 2 f_i + f_{i+1}) / h^2
        // normalised so the lhs weights sum to 1
        auto c1 = compact_scheme(1);
        CHECK(c1.is_compact());
        CHECK(c1.order == 4);
        CHECK(std::abs(c1.lhs_weights[0] - 1.0 / 6.0) < eps);
        CHECK(std::abs(c1.lhs_weights[1] - 2.0 / 3.0) < eps);
        CHECK(std::abs(c1.lhs_weights[2] - 1.0 / 6.0) < eps);
        CHECK(std::abs(c1.weights[0] + 0.5) < eps);
        CHECK(std::abs(c1.weights[2] - 0.5) < eps);

        auto c2 = compact_scheme(2);
        CHECK(c2.is_compact());
        CHECK(c2.order == 4);
        CHECK(std::abs(c2.lhs_weights[0] - 1.0 / 12.0) < eps);
        CHECK(std::abs(c2.lhs_weights[1] - 10.0 / 12.0) < eps);
        CHECK(std::abs(c2.weights[1] + 2.0) < eps);

        // an explicit scheme is not compact, and the point operators reject compact ones
        CHECK(!central_scheme(1, 2).is_compact());
        CHECK_THROWS_AS(nabla([](vec3d const& q) { return scalar3d(q.x); },
                              vec3d{1.0, 0.0, 0.0}, compact_scheme(1)),
                        std::invalid_argument);

        // a test function with known derivatives
        auto fn = [](double s) { return std::sin(3 * s) * std::exp(0.5 * s); };
        auto d1 = [](double s) {
            return std::exp(0.5 * s) * (3 * std::cos(3 * s) + 0.5 * std::sin(3 * s));
        };
        auto max_err = [&](std::vector<double> const& x, fd_kind kind) {
            std::vector<double> f;
            for (auto s : x)
                f.push_back(fn(s));
            auto d = fd_derivative(x, f, 1, kind, 2);
            double m = 0.0;
            for (std::size_t i = 2; i + 2 < x.size(); ++i) {
                m = std::max(m, std::abs(d[i] - d1(x[i])));
            }
            return m;
        };
        auto uniform = [](int n) {
            std::vector<double> x;
            for (int i = 0; i < n; ++i)
                x.push_back(double(i) / (n - 1));
            return x;
        };
        auto stretched = [](int n) { // smooth geometric stretching
            std::vector<double> x;
            for (int i = 0; i < n; ++i) {
                double const s = double(i) / (n - 1);
                x.push_back(s * s * (3 - 2 * s));
            }
            return x;
        };

        // EQUIDISTANT: compact reaches 4th order from a three-point footprint where the
        // explicit stencil of the same width reaches only 2nd
        double const e41 = max_err(uniform(41), fd_kind::explicit_fd);
        double const c41 = max_err(uniform(41), fd_kind::compact_fd);
        double const c81 = max_err(uniform(81), fd_kind::compact_fd);
        CHECK(c41 < e41 / 100.0);
        double const observed = std::log2(c41 / c81);
        CHECK(observed > 3.7);
        CHECK(observed < 4.3);

        // NON-EQUIDISTANT: the grid coordinates are the input, so the weights differ
        // from node to node
        auto xs = stretched(41);
        CHECK(max_err(xs, fd_kind::compact_fd) <
              max_err(xs, fd_kind::explicit_fd) / 10.0);
        auto s5 = scheme_at(xs, 5, 1, fd_kind::compact_fd);
        auto s30 = scheme_at(xs, 30, 1, fd_kind::compact_fd);
        CHECK(std::abs(s5.lhs_weights[0] - s30.lhs_weights[0]) > 1.0e-6);
        // the first-derivative compact scheme keeps 4th order on a stretched grid,
        // the second-derivative one drops to 3rd -- ask, do not assume
        CHECK(s5.order == 4);
        CHECK(scheme_at(xs, 5, 2, fd_kind::compact_fd).order == 3);

        // the grid path is value-type generic: multivector samples work unchanged
        auto xu = uniform(41);
        std::vector<mvec3d> mf;
        for (auto s : xu) {
            mf.push_back(mvec3d(scalar3d(std::sin(3 * s))) +
                         mvec3d(std::cos(2 * s) * e1_3d));
        }
        auto md = fd_derivative(xu, mf, 1, fd_kind::compact_fd);
        double mm = 0.0;
        for (std::size_t i = 2; i + 2 < xu.size(); ++i) {
            mm = std::max(mm, std::abs(value_t(gr0(md[i])) - 3 * std::cos(3 * xu[i])));
            mm = std::max(mm, nrm(gr1(md[i]) - (-2 * std::sin(2 * xu[i])) * e1_3d));
        }
        CHECK(mm < 1.0e-5);
    }


    TEST_CASE("Vec3d: reciprocal_frame")
    {
        fmt::println("Vec3d: reciprocal_frame");

        auto a1 = vec3d{1.0, 0.0, 0.0};
        auto a2 = vec3d{1.0, 1.0, 0.0};
        auto a3 = vec3d{0.0, 1.0, 2.0}; // skewed
        auto r = reciprocal_frame(a1, a2, a3);
        vec3d const a[3] = {a1, a2, a3};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                CHECK(std::abs(value_t(dot(r[i], a[j])) - (i == j ? 1.0 : 0.0)) < eps);
            }
        }
        // orthonormal: the reciprocal frame is the frame itself
        auto ro = reciprocal_frame(e1_3d, e2_3d, e3_3d);
        CHECK(is_close(ro[0], e1_3d));
        CHECK(is_close(ro[1], e2_3d));
        CHECK(is_close(ro[2], e3_3d));
        // merely ORTHOGONAL: the shortcut a^i = inv(a_i) must agree
        auto b1 = 2.0 * e1_3d, b2 = 3.0 * e2_3d, b3 = 0.5 * e3_3d;
        auto rb = reciprocal_frame(b1, b2, b3);
        CHECK(is_close(rb[0], inv(b1)));
        CHECK(is_close(rb[1], inv(b2)));
        CHECK(is_close(rb[2], inv(b3)));
        // a degenerate frame is rejected explicitly, not left to produce infinities
        CHECK_THROWS_AS(reciprocal_frame(e1_3d, e2_3d, e1_3d + e2_3d),
                        std::invalid_argument);
    }

} // EGA 3D Tests
