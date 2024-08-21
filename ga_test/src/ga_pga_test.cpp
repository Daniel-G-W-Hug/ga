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

    TEST_CASE("pga_2d")
    {
        fmt::println("");
        fmt::println("pga_2d:");
        CHECK(true);
    }
} // TEST_SUITE("Projective Geometric Algebra (PGA)")