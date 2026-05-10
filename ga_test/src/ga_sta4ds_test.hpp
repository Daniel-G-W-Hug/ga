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

    TEST_CASE("G<3,1,0>: sta4ds")
    {
        fmt::println("G<3,1,0>: sta4ds");
        // Space-Time Algebra with signature (+,+,+,-)
        const algebra<3, 1, 0> alg;
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
        pscalar4ds ps{3.0};

        // same tag comparisons should work
        CHECK(s1 == s2);
        CHECK(s1 != scalar4ds{4.0});
        CHECK(s1 < scalar4ds{4.0});
        CHECK(s1 <= scalar4ds{3.0});
        CHECK(scalar4ds{4.0} > s1);
        CHECK(scalar4ds{3.0} >= s1);

        // different floating point types with same tag should work
        Scalar_t<float, scalar4ds_tag> sf{3.0f};
        CHECK(s1 == sf); // double vs float with same tag
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Dual Mathematical Properties Tests
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: dual composition properties")
    {
        fmt::println("G<1,3,0>: dual composition properties");

        // For STA4ds (even-dimensional): l_dual(r_dual(u)) =
        // r_dual(l_dual(u)) = u This follows from
        // l_cmpl(r_cmpl(u)) = u

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

        // l_dual(r_dual(u)) = u
        CHECK(l_dual(r_dual(s)) == s);
        CHECK(l_dual(r_dual(v)) == v);
        CHECK(l_dual(r_dual(B)) == B);
        CHECK(l_dual(r_dual(T)) == T);
        CHECK(l_dual(r_dual(ps)) == ps);
        CHECK(l_dual(r_dual(M_e)) == M_e);
        CHECK(l_dual(r_dual(M)) == M);

        // r_dual(l_dual(u)) = u
        CHECK(r_dual(l_dual(s)) == s);
        CHECK(r_dual(l_dual(v)) == v);
        CHECK(r_dual(l_dual(B)) == B);
        CHECK(r_dual(l_dual(T)) == T);
        CHECK(r_dual(l_dual(ps)) == ps);
        CHECK(r_dual(l_dual(M_e)) == M_e);
        CHECK(r_dual(l_dual(M)) == M);

        fmt::println("l_dual(r_dual(u)) = r_dual(l_dual(u)) = u");
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

    TEST_CASE("G<3,1,0>: MVec4ds_E and MVec4ds_U formatting tests")
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

    TEST_CASE("G<3,1,0>: sta4ds")
    {
        fmt::println("");
        fmt::println("G<3,1,0>: sta4ds");
        // 4ds euklidean geometric algebra
        const algebra<3, 1, 0> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 1);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   sta4ds: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   sta4ds: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    TEST_CASE("G<3,1,0>: defining basic types and ctor checks")
    {
        fmt::println("G<3,1,0>: defining basic types and ctor checks");

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

    TEST_CASE("G<3,1,0>: Scalar4ds and PScalar4ds formatting tests")
    {
        fmt::println("G<3,1,0>: Scalar4ds and PScalar4ds formatting tests");

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

    TEST_CASE("G<4,0,0>: Vec4ds and TriVec4ds formatting tests")
    {
        fmt::println("G<4,0,0>: Vec4ds and TriVec4ds formatting tests");

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

    TEST_CASE("G<4,0,0>: BiVec4ds formatting tests")
    {
        fmt::println("G<4,0,0>: BiVec4ds formatting tests");

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

    TEST_CASE("G<4,0,0>: MVec4ds formatting tests")
    {
        fmt::println("G<4,0,0>: MVec4ds formatting tests");

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

} // STA 3D Tests
