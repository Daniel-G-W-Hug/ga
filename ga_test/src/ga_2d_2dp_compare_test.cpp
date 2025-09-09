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
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// Euclidean geometric algebra of two-dimensional space: ega2d, pga2dp comparison
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA2D PGA2DP Comparison Tests")
{

    TEST_CASE("ega2d pga2dp lines")
    {

        fmt::println("");

        auto u = vec2d{2, 0.5};
        auto a = vec2d{2, 5};

        auto M = wdg(a, u);
        auto d = M * inv(u);

        auto L = u + M;

        fmt::println("u = {}, nrm(u) = {}", u, nrm(u));
        fmt::println("a = {}, nrm(a) = {}", a, nrm(a));
        fmt::println("M = wdg(a,u) = {}", M);
        fmt::println("d = M * inv(u) = {}, nrm(d) = {}", d, nrm(d));
        fmt::println("nrm(u) * nrm(d) = {}", nrm(u) * nrm(d));
        fmt::println("L = {}", L);


        fmt::println("");
        auto up = vec2dp{2, 0.5, 0.0};
        auto ap = vec2dp{2, 5, 1.0};

        auto Lp = wdg(ap, up);
        auto dp = unitize(support2dp(Lp));


        fmt::println("up = {}, bulk_nrm(up) = {}", up, bulk_nrm(up));
        fmt::println("ap = {}, bulk_nrm(up) = {}", ap, bulk_nrm(ap));
        fmt::println("dp = unitize(support2dp(Lp)) = {}, bulk_nrm(dp) = {}", dp,
                     bulk_nrm(dp));
        fmt::println("Lp = {}", Lp);
        fmt::println("att(Lp) = {}", att(Lp));   // direction vector
        fmt::println("bulk(Lp) = {}", bulk(Lp)); // resulting moment w.r.t. origion

        fmt::println("");
    }

} // EGA2D PGA2DP Comparison Tests
