// author: Daniel Hug, 2024

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_pga.hpp"

using namespace hd::ga; // use ga types, constants, etc.
// using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)

TEST_SUITE("Projective Geometric Algebra (PGA)")
{

    TEST_CASE("algebra<2,0,1> - pga_2d")
    {
        fmt::println("");
        fmt::println("algebra<2, 0, 1> - pga_2d:");
        // 2d projective geometric algebra (modelled in a 3d algebra)
        const algebra<2, 0, 1> alg;
        CHECK(alg.p() == 2);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("   pga_2d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga_2d: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("algebra<3, 0, 1> pga_3d")
    {
        fmt::println("");
        fmt::println("algebra<3, 0, 1> - pga_3d:");
        // 3d projective geometric algebra (modelled in a 4d algebra)
        const algebra<3, 0, 1> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   pga_3d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga_3d: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("pga_2d")
    {
        fmt::println("");
        fmt::println("pga_2d:");
        CHECK(true);
    }
} // TEST_SUITE("Projective Geometric Algebra (PGA)")