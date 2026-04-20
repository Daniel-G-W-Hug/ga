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
// Space-Time Algebra: STA4D (G<1,3,0>)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("STA 3D Tests")
{

    TEST_CASE("G<1,3,0>: sta4d")
    {
        fmt::println("G<1,3,0>: sta4d");
        // Space-Time Algebra with signature (+,-,-,-)
        const algebra<1, 3> alg;
        CHECK(alg.p() == 1);
        CHECK(alg.n() == 3);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   sta4d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   sta4d: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    TEST_CASE("ScalarSta: comparison")
    {
        fmt::println("ScalarSta: comparison");
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
        Scalar_t<float, scalarsta4d_tag> sf{3.0f};
        CHECK(s1 == sf); // double vs float with same tag
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Dual Mathematical Properties Tests
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("G<1,3,0>: dual composition properties")
    {
        fmt::println("G<1,3,0>: dual composition properties");

        // For STA4D (even-dimensional): l_dual(r_dual(u)) =
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

} // STA 3D Tests
