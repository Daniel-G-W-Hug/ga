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
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// projective geometric algebra pga2dp application tests
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA 2DP Application Tests")
{

    TEST_CASE("pga2dp: steering angle application")
    {
        fmt::println("pga2dp: steering angle application");

        // wheel positions rear (pr) and front (pf)
        constexpr double track = 2.0;
        constexpr double wheelbase = 3.0;
        // rear wheel positions
        auto prl = point2dp{-0.5 * track, 0, 1};
        auto prr = point2dp{0.5 * track, 0, 1};
        // front wheel postitions
        auto pfl = point2dp{-0.5 * track, wheelbase, 1};
        auto pfr = point2dp{0.5 * track, wheelbase, 1};

        // rear axis - connecting left and right rear wheels
        auto axis_r = unitize(wdg(prl, prr));

        const int step_max = 10;
        const double max_angle_deg = 30;

        // right front wheel direction defines turning angle
        for (int n = -step_max; n <= step_max; ++n) {
            // for (int n = 0; n <= step_max; ++n) {

            // steering angle phi of right front wheel
            double phi = n * deg2rad(max_angle_deg) / step_max;

            // axis of right front wheel defined by postion and direction
            auto axis_fr = unitize(wdg(pfr, vec2dp(cos(phi), sin(phi), 0)));

            // turning center is intersection point of both axis
            auto turn_ctr = rwdg(axis_fr, axis_r);
            // unitize to get projective point for phi != 0
            if (std::abs(turn_ctr.z) > eps) {
                turn_ctr = unitize(turn_ctr);
            }

            // now create axis_fl from turning point and position of left front wheel
            auto axis_fl = unitize(wdg(turn_ctr, pfl));

            fmt::print("n = {:>3}", n);
            if (phi < 0.0) {
                fmt::print(", phi_l = {:>-5.1f}°", rad2deg(angle(axis_r, axis_fl)) - 180);
            }
            else {
                fmt::print(", phi_l = {:>-5.1f}°", rad2deg(angle(axis_r, axis_fl)));
            }
            fmt::print(", phi_r (deg) = {:>-5.1f}°", rad2deg(phi));
            // fmt::print(", att(axis_fr) = {:>-6.3f}", att(axis_fr));
            // fmt::print(", axis_fr = {:>-8.3f}", axis_fr);
            // fmt::print(", axis_fl = {:>-8.3f}", axis_fl);
            fmt::println(", turn_ctr = {:>-8.3f}", turn_ctr);
        }
        fmt::println("");
    }

} // PGA 2DP Application Tests