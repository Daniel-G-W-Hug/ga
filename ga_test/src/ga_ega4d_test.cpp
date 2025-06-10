// Copyright 2024-2025, Daniel Hug. All rights reserved.

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

} // EGA 4D Tests