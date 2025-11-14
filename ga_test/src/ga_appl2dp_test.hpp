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

TEST_SUITE("PGA2DP: application tests")
{

    TEST_CASE("pga2dp: application tests - Ackermann steering angle application")
    {
        fmt::println("pga2dp: application tests - Ackermann steering angle application");

        // wheel positions rear (pr) and front (pf)
        constexpr double track = 2.0;
        constexpr double wheelbase = 3.0;
        // wheel positions (r: rear, f: front, l:left, r:right )
        auto p_rl = point2dp{-0.5 * track, 0, 1};
        auto p_rr = point2dp{0.5 * track, 0, 1};
        auto p_fl = point2dp{-0.5 * track, wheelbase, 1};
        auto p_fr = point2dp{0.5 * track, wheelbase, 1};

        // rear axis - connecting left and right rear wheels
        auto axis_r = unitize(wdg(p_rl, p_rr));

        const int step_max = 10;
        const double max_angle_deg = 25;

        // right front wheel direction defines turning angle
        for (int n = -step_max; n <= step_max; ++n) {

            // steering angle phi_fr of right front wheel
            double phi_fr = n * deg2rad(max_angle_deg) / step_max;

            // axis of right front wheel defined by postion and direction
            auto axis_fr = unitize(wdg(p_fr, vec2dp(cos(phi_fr), sin(phi_fr), 0)));

            // turning center is intersection point of both axis
            auto turn_ctr = rwdg(axis_r, axis_fr);
            // unitize to get projective point for phi_fr != 0 (i.e. not at infinity)
            if (abs(turn_ctr.z) > eps) {
                turn_ctr = unitize(turn_ctr);
            }

            // now create axis_fl from turning point and position of left front wheel
            auto axis_fl = unitize(wdg(turn_ctr, p_fl));

            fmt::print("n = {:>3}", n);
            if (phi_fr <= 0.0) {
                fmt::print(", phi_fl = {:>-5.1f}°",
                           rad2deg(angle(axis_r, axis_fl)) - 180);
            }
            else {
                fmt::print(", phi_fl = {:>-5.1f}°", rad2deg(angle(axis_r, axis_fl)));
            }
            fmt::print(", phi_fr = {:>-5.1f}°", rad2deg(phi_fr));
            fmt::println(", turn_ctr = {:>-8.3f}", turn_ctr);
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: application tests - Four wheel steering angle application")
    {
        fmt::println("pga2dp: application tests - Four wheel steering angle application");

        // wheel positions rear (pr) and front (pf)
        constexpr double track = 2.0;
        constexpr double wheelbase = 3.0;
        // wheel positions (r: rear, f: front, l:left, r:right )
        auto p_rl = point2dp{-0.5 * track, 0, 1};
        auto p_rr = point2dp{0.5 * track, 0, 1};
        auto p_fl = point2dp{-0.5 * track, wheelbase, 1};
        auto p_fr = point2dp{0.5 * track, wheelbase, 1};

        // rear axis - virtually set to 20% of wheelbase in front of rear wheels
        auto axis_r = unitize(wdg(point2dp{0, 0.2 * wheelbase, 1}, vec2dp{1, 0, 0}));

        const int step_max = 10;
        const double max_angle_deg = 25;

        // right front wheel direction defines turning angle
        for (int n = -step_max; n <= step_max; ++n) {

            // steering angle phi_fr of right front wheel
            double phi_fr = n * deg2rad(max_angle_deg) / step_max;

            // axis of right front wheel defined by postion and direction
            auto axis_fr = unitize(wdg(p_fr, vec2dp(cos(phi_fr), sin(phi_fr), 0)));

            // turning center is intersection point of both axis
            auto turn_ctr = rwdg(axis_r, axis_fr);
            // unitize to get projective point for phi_fr != 0 (i.e. not at infinity)
            if (abs(turn_ctr.z) > eps) {
                turn_ctr = unitize(turn_ctr);
            }

            // now create axis from turning point and position of remaining wheels
            auto axis_fl = unitize(wdg(turn_ctr, p_fl));
            auto axis_rl = unitize(wdg(turn_ctr, p_rl));
            auto axis_rr = unitize(wdg(turn_ctr, p_rr));

            fmt::print("n = {:>3}", n);
            if (phi_fr <= 0.0) {
                fmt::print(", phi_fl = {:>-5.1f}°",
                           rad2deg(angle(axis_r, axis_fl)) - 180);
            }
            else {
                fmt::print(", phi_fl = {:>-5.1f}°", rad2deg(angle(axis_r, axis_fl)));
            }
            fmt::print(", phi_fr = {:>-5.1f}°", rad2deg(phi_fr));
            fmt::println(", turn_ctr = {:>-8.3f}", turn_ctr);

            if (phi_fr <= 0.0) {
                fmt::print("         phi_rl = {:>-5.1f}°",
                           -rad2deg(angle(axis_r, axis_rl)) + 180);
                fmt::println(", phi_rr = {:>-5.1f}°",
                             -rad2deg(angle(axis_r, axis_rr)) + 180);
            }
            else {
                fmt::print("         phi_rl = {:>-5.1f}°",
                           -rad2deg(angle(axis_r, axis_rl)));
                fmt::println(", phi_rr = {:>-5.1f}°", -rad2deg(angle(axis_r, axis_rr)));
            }
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: application tests - force and moment")
    {
        fmt::println("pga2dp: application tests - force and moment");

        // a line modelling a force and its moment w.r.t. an arbitrary point:
        // ------------------------------------------------------------------

        // point P the force is acting at (position as projective point)
        auto P = vec2dp{2, 1, 1};

        // force vector (contains direction and magnitude)
        auto f = vec2dp{3, 1, 0};
        f = f / to_val(bulk_nrm(f)) * 2.0; // force of 2 N in direction of f

        // resulting force line L modelling force and moment (bivector = bound line)
        auto L = wdg(P, f);

        // reference point R for moment calculation
        auto R = vec2dp{4, 0.5, 1};

        // a free bivector modelling a free moment of 4 Nm, e.g. generated by a force pair
        auto B = bivec2dp(0, 0, 4);

        // operations to retrieve the force f from force line L:
        // -----------------------------------------------------

        auto F = att(L); // direction and magnitude of force

        // moment of force line at reference point R, just from line and reference point
        // distance of R to L calculated from orthogonal projection of R onto L
        // moment = force ^ lever arm
        auto d = R - ortho_proj2dp(R, L);
        auto M_LR = wdg(d, F); // lever arm ^ force

        // alternative approach: calculate moment from support point D and f
        // => L = P ^ f = |P|*sin(P,f)*|f| = |D|*|f| = D ^ f
        //                                   (with D as perpendicular distance to point)
        // However, this needs a known point D on the line (best option if D is known)
        // or D could be obtained from: L = D ^ f via the inverse of f = f^(-1)
        // => D = L * f^(-1) = gr1(L * inv(att(L))
        //    D is the support, if R is the origin, but R can be an arbitrary point
        //
        // The resulting moment on the point R is M = d ^ f = (D - R) ^ att(L)
        //
        auto M_LR_calc = wdg(R - gr1(L * inv(F)), F);

        fmt::println("point P = {}", P);
        fmt::println("f       = {}", f);
        fmt::println("force line L = {}", L);
        fmt::println("support(L) = {}", support(L));
        fmt::println("");
        fmt::println("reference point R  = {}", R);
        fmt::println("distance to line d = {}", d);
        fmt::println("");
        fmt::println("retrieved force F  = att(L) = {}", F);
        fmt::println("retrieved D = L*inv(att(L)) = {}", L * inv(att(L)));
        fmt::println("retrieved D = gr1(L*inv(att(L))) = {}", gr1(L * inv(att(L))));
        fmt::println("retrieved d = to_val(bulk_nrm(D)) = {}",
                     to_val(bulk_nrm(gr1(L * inv(att(L))))));
        fmt::println("");
        fmt::println("moment at R:  M_LR = (R - proj(R, L)) ^ F = {}", M_LR);
        fmt::println("altern. calc. M_LR_calc                   = {}", M_LR_calc);
        fmt::println("");
        double M_LR_calc_free = 0.0;
        if (att(B) == vec2dp{}) { // identified a free moment, which acts on any point
            M_LR_calc_free = B.z;
        }
        fmt::println("calc. of free moment at R from bivector B = {}", M_LR_calc_free);


        CHECK(gr0(f * inv(f)) == scalar2dp{1});
        CHECK(gr2(f * inv(f)) == bivec2dp{});

        CHECK(f == F);
        CHECK(support(L) == gr1(L * inv(att(L))));
        CHECK(pscalar2dp{0.0} == gr3(L * inv(att(L))));

        fmt::println("");
    }

} // PGA 2DP Application Tests