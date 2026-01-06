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


/////////////////////////////////////////////////////////////////////////////////////////
// Euclideon geometric algebra of four-dimensional space: ega4d
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA 4D Tests")
{

    TEST_CASE("G<4,0,0>: ega4d")
    {
        fmt::println("");
        fmt::println("G<4,0,0>: ega4d");
        // 4d euklidean geometric algebra
        const algebra<4> alg;
        CHECK(alg.p() == 4);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   ega4d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   ega4d: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    TEST_CASE("G<4,0,0>: defining basic types and ctor checks")
    {
        fmt::println("G<4,0,0>: defining basic types and ctor checks");

        auto s1 = scalar4d{5.0};
        auto v1 = vec4d{1.0, 2.0, 3.0, 4.0};
        auto b1 = bivec4d{-1.0, 2.0, -3.0, 4.0, 5.0, -6.0};
        auto t1 = trivec4d{-1.0, -2.0, -3.0, -4.0};
        auto ps1 = pscalar4d{-5.0};

        auto mv1 = mvec4d{s1};
        auto mv2 = mvec4d{v1};
        auto mv3 = mvec4d{b1};
        auto mv4 = mvec4d{t1};
        auto mv5 = mvec4d{ps1};

        auto mv6a = mvec4d{s1, b1, ps1};
        auto mv6b = mvec4d_e{s1};
        auto mv6c = mvec4d_e{b1};
        auto mv6d = mvec4d_e{ps1};
        auto mv6e = mvec4d_e{s1, b1, ps1};
        auto mv6f = mvec4d{mv6e};

        auto mv7a = mvec4d{v1, t1};
        auto mv7b = mvec4d_u{v1};
        auto mv7c = mvec4d_u{t1};
        auto mv7d = mvec4d_u{v1, t1};
        auto mv7e = mvec4d{mv7d};

        auto mv8 = mvec4d{s1, v1, b1, t1, ps1};

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
        CHECK(gr1(mv6a) == vec4d{});
        CHECK(gr2(mv6a) == b1);
        CHECK(gr3(mv6a) == trivec4d{});
        CHECK(gr4(mv6a) == ps1);

        CHECK(gr0(mv6b) == s1);
        CHECK(gr2(mv6b) == bivec4d{});
        CHECK(gr4(mv6b) == pscalar4d{});

        CHECK(gr0(mv6c) == scalar4d{});
        CHECK(gr2(mv6c) == b1);
        CHECK(gr4(mv6c) == pscalar4d{});

        CHECK(gr0(mv6d) == scalar4d{});
        CHECK(gr2(mv6d) == bivec4d{});
        CHECK(gr4(mv6d) == ps1);

        CHECK(gr0(mv6e) == s1);
        CHECK(gr2(mv6e) == b1);
        CHECK(gr4(mv6e) == ps1);

        CHECK(gr0(mv6f) == s1);
        CHECK(gr1(mv6f) == vec4d{});
        CHECK(gr2(mv6f) == b1);
        CHECK(gr3(mv6f) == trivec4d{});
        CHECK(gr4(mv6f) == ps1);

        CHECK(gr0(mv7a) == scalar4d{});
        CHECK(gr1(mv7a) == v1);
        CHECK(gr2(mv7a) == bivec4d{});
        CHECK(gr3(mv7a) == t1);
        CHECK(gr4(mv7a) == pscalar4d{});

        CHECK(gr1(mv7b) == v1);
        CHECK(gr3(mv7b) == trivec4d{});

        CHECK(gr1(mv7c) == vec4d{});
        CHECK(gr3(mv7c) == t1);

        CHECK(gr1(mv7d) == v1);
        CHECK(gr3(mv7d) == t1);

        CHECK(gr0(mv7e) == scalar4d{});
        CHECK(gr1(mv7e) == v1);
        CHECK(gr2(mv7e) == bivec4d{});
        CHECK(gr3(mv7e) == t1);
        CHECK(gr4(mv7e) == pscalar4d{});

        CHECK(gr0(mv8) == s1);
        CHECK(gr1(mv8) == v1);
        CHECK(gr2(mv8) == b1);
        CHECK(gr3(mv8) == t1);
        CHECK(gr4(mv8) == ps1);
    }

    TEST_CASE("G<4,0,0>: Scalar4d and PScalar4d formatting tests")
    {
        fmt::println("G<4,0,0>: Scalar4d and PScalar4d formatting tests");

        // Test Scalar4d formatting
        Scalar4d<double> scalar_val{2.23607};

        // Basic output
        std::string scalar_basic = fmt::format("{}", scalar_val);
        CHECK(scalar_basic == "Scalar4d(2.23607)");

        // Two decimal places
        std::string scalar_two_decimals = fmt::format("{:.2f}", scalar_val);
        CHECK(scalar_two_decimals == "Scalar4d(2.24)");

        // Scientific notation
        std::string scalar_scientific = fmt::format("{:.2e}", scalar_val);
        CHECK(scalar_scientific == "Scalar4d(2.24e+00)");

        // Contextual usage
        std::string scalar_contextual = fmt::format("Norm: {:.3f}", scalar_val);
        CHECK(scalar_contextual == "Norm: Scalar4d(2.236)");

        // Test PScalar4d formatting
        PScalar4d<double> pscalar_val{2.44949};

        // Basic output
        std::string pscalar_basic = fmt::format("{}", pscalar_val);
        CHECK(pscalar_basic == "PScalar4d(2.44949)");

        // Two decimal places
        std::string pscalar_two_decimals = fmt::format("{:.2f}", pscalar_val);
        CHECK(pscalar_two_decimals == "PScalar4d(2.45)");

        // Scientific notation
        std::string pscalar_scientific = fmt::format("{:.2e}", pscalar_val);
        CHECK(pscalar_scientific == "PScalar4d(2.45e+00)");

        // Contextual usage
        std::string pscalar_contextual = fmt::format("4D Volume: {:.3f}", pscalar_val);
        CHECK(pscalar_contextual == "4D Volume: PScalar4d(2.449)");

        fmt::println("   Scalar4d basic: {}", scalar_basic);
        fmt::println("   Scalar4d 2-decimal: {}", scalar_two_decimals);
        fmt::println("   Scalar4d scientific: {}", scalar_scientific);
        fmt::println("   Scalar4d contextual: {}", scalar_contextual);

        fmt::println("   PScalar4d basic: {}", pscalar_basic);
        fmt::println("   PScalar4d 2-decimal: {}", pscalar_two_decimals);
        fmt::println("   PScalar4d scientific: {}", pscalar_scientific);
        fmt::println("   PScalar4d contextual: {}", pscalar_contextual);
    }

    TEST_CASE("G<4,0,0>: Vec4d and TriVec4d formatting tests")
    {
        fmt::println("G<4,0,0>: Vec4d and TriVec4d formatting tests");

        // Test Vec4d formatting
        Vec4d<double> vec4d_val{3.14159, 2.71828, 1.41421, 1.73205};

        // Basic output
        std::string basic = fmt::format("{}", vec4d_val);
        CHECK(basic == "Vec4d(3.14159,2.71828,1.41421,1.73205)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", vec4d_val);
        CHECK(two_decimals == "Vec4d(3.14,2.72,1.41,1.73)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", vec4d_val);
        CHECK(scientific == "Vec4d(3.14e+00,2.72e+00,1.41e+00,1.73e+00)");

        // Contextual usage
        std::string contextual = fmt::format("Position: {:.2f}", vec4d_val);
        CHECK(contextual == "Position: Vec4d(3.14,2.72,1.41,1.73)");

        // Test TriVec4d formatting
        TriVec4d<double> trivec4d_val{1.2, 2.3, 3.4, 4.5};

        // Basic output
        std::string trivec_basic = fmt::format("{}", trivec4d_val);
        CHECK(trivec_basic == "TriVec4d(1.2,2.3,3.4,4.5)");

        // Two decimal places
        std::string trivec_two_decimals = fmt::format("{:.2f}", trivec4d_val);
        CHECK(trivec_two_decimals == "TriVec4d(1.20,2.30,3.40,4.50)");

        // Scientific notation
        std::string trivec_scientific = fmt::format("{:.2e}", trivec4d_val);
        CHECK(trivec_scientific == "TriVec4d(1.20e+00,2.30e+00,3.40e+00,4.50e+00)");

        // Contextual usage
        std::string trivec_contextual = fmt::format("Volume: {:.2f}", trivec4d_val);
        CHECK(trivec_contextual == "Volume: TriVec4d(1.20,2.30,3.40,4.50)");

        fmt::println("   Vec4d basic: {}", basic);
        fmt::println("   Vec4d 2-decimal: {}", two_decimals);
        fmt::println("   Vec4d scientific: {}", scientific);
        fmt::println("   Vec4d contextual: {}", contextual);

        fmt::println("   TriVec4d basic: {}", trivec_basic);
        fmt::println("   TriVec4d 2-decimal: {}", trivec_two_decimals);
        fmt::println("   TriVec4d scientific: {}", trivec_scientific);
        fmt::println("   TriVec4d contextual: {}", trivec_contextual);
    }

    TEST_CASE("G<4,0,0>: BiVec4d formatting tests")
    {
        fmt::println("G<4,0,0>: BiVec4d formatting tests");

        // Test BiVec4d formatting (BVec6_t with 6 components: vx, vy, vz, mx, my, mz)
        BiVec4d<double> bivec4d_val{1.1, 2.2, 3.3, 4.4, 5.5, 6.6};

        // Basic output
        std::string basic = fmt::format("{}", bivec4d_val);
        CHECK(basic == "BiVec4d(1.1,2.2,3.3,4.4,5.5,6.6)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", bivec4d_val);
        CHECK(two_decimals == "BiVec4d(1.10,2.20,3.30,4.40,5.50,6.60)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", bivec4d_val);
        CHECK(scientific ==
              "BiVec4d(1.10e+00,2.20e+00,3.30e+00,4.40e+00,5.50e+00,6.60e+00)");

        // Contextual usage
        std::string contextual = fmt::format("Plane: {:.2f}", bivec4d_val);
        CHECK(contextual == "Plane: BiVec4d(1.10,2.20,3.30,4.40,5.50,6.60)");

        fmt::println("   BiVec4d basic: {}", basic);
        fmt::println("   BiVec4d 2-decimal: {}", two_decimals);
        fmt::println("   BiVec4d scientific: {}", scientific);
        fmt::println("   BiVec4d contextual: {}", contextual);
    }

    TEST_CASE("G<4,0,0>: MVec4d formatting tests")
    {
        fmt::println("G<4,0,0>: MVec4d formatting tests");

        // Test MVec4d formatting (MVec16_t with 16 components: c0 through c15)
        MVec4d<double> mvec4d_val{1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,
                                  9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

        // Basic output
        std::string basic = fmt::format("{}", mvec4d_val);
        CHECK(basic == "MVec4d(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)");

        // Two decimal places
        std::string two_decimals = fmt::format("{:.2f}", mvec4d_val);
        CHECK(two_decimals == "MVec4d(1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,"
                              "9.00,10.00,11.00,12.00,13.00,14.00,15.00,16.00)");

        // Scientific notation
        std::string scientific = fmt::format("{:.2e}", mvec4d_val);
        CHECK(scientific == "MVec4d(1.00e+00,2.00e+00,3.00e+00,4.00e+00,5.00e+00,"
                            "6.00e+00,7.00e+00,8.00e+00,9.00e+00,1.00e+01,1.10e+01,"
                            "1.20e+01,1.30e+01,1.40e+01,1.50e+01,1.60e+01)");

        // Contextual usage
        std::string contextual = fmt::format("Full multivector: {:.1f}", mvec4d_val);
        CHECK(contextual == "Full multivector: MVec4d(1.0,2.0,3.0,4.0,5.0,6.0,7.0,"
                            "8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0)");

        fmt::println("   MVec4d basic: {}", basic);
        fmt::println("   MVec4d 2-decimal: {}", two_decimals);
        fmt::println("   MVec4d scientific: {}", scientific);
        fmt::println("   MVec4d contextual: {}", contextual);
    }

} // EGA 4D Tests