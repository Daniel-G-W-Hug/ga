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

TEST_SUITE("PGA3DP: application tests")
{

    TEST_CASE("pga3dp: reference and tumbling plane (for grinding application)")
    {

        fmt::println("pga3dp: reference and tumbling plane (for grinding application)");
        fmt::println("");

        // define a plane by 3 measuring points A,B,C in coordinate system (e1,e2,e3,e4)
        // with e4 as the projective dimension (e4 = 0 for directions, e4 = 1 for points)
        // (here idealized plane close to e1, e2 plane, derived from fictive measurement)
        auto A = vec3dp(1.0, 0.0, 0.0, 1);
        auto B = vec3dp(1.0, 1.0, 0.01, 1);
        auto C = vec3dp(0.0, 1.0, 0.0, 1);

        // compute the unitized plane and its normal:
        auto p_ref = unitize(wdg(wdg(A, B), C)); // plane equation
        auto p_n = left_weight_dual(p_ref);      // normal of plane

        fmt::println("A = {}", A);
        fmt::println("B = {}", B);
        fmt::println("C = {}", C);
        fmt::println("           p_ref  = {:>-8.5f}", p_ref);
        fmt::println("           p_n    =    {:>-8.5f}, |p_n| = {}", p_n,
                     to_val(bulk_nrm(p_n)));
        fmt::println("");

        // define a reference point T (=center of rotating grinding tool)
        // located close to our measured reference plane
        auto T = vec3dp(3.0, 1.0, 0.0, 1);

        // normal of the tool plane is vector z in e3 direction
        // z is potentially tumbling with radius r at tip in e1, e2 plane
        auto z = vec3dp(0.0, 0.0, 1.0, 0);

        for (auto phi = 0.0; phi < 2 * pi; phi += deg2rad(15)) {

            auto r_dim = 0.01; // radius of tumbling (1% of normal vector)
            auto r = r_dim * vec3dp(cos(phi), sin(phi), 0.0, 0);

            // Point R is at tip of tumbling vector, z_n is normalized vector towards R
            auto R = T + z + r;
            auto z_n = (R - T) / to_val(bulk_nrm(R - T)); // normalized vector from T to
                                                          // R which is the normal of
                                                          // the current plane p_tool

            fmt::println("phi = {:>3.0f}: z_n         =    {:>-8.5f}", rad2deg(phi), z_n);

            // T and z_n form a line from which we can calculate the plane
            // 1st step: line through T in direction of z_n
            auto l_h = unitize(wdg(z_n, T));
            // 2nd step: plane containing T, perpendicular to l
            auto p_tool = unitize(wdg(T, right_weight_dual(l_h)));
            fmt::println("           p_tool      = {:>-8.5f}", p_tool);

            // both planes intersect in a line l
            auto l = rwdg(p_tool, p_ref);
            // print direction vector of line l -> att(l)
            fmt::println("           att(l)      =    {:>-8.5f}", att(l));
            fmt::println("           angle(l,e1) =    {:>-5.2f}°",
                         rad2deg(angle(att(l), e1_3dp)));

            fmt::println("");
        }
    }

    TEST_CASE("pga3dp: intersection of non-intersecting lines")
    {

        fmt::println("pga3dp: intersection of non-intersecting lines");

        auto l1 = x_axis_3dp;
        auto l2 = wdg(vec3dp{0, 1, 0, 1}, vec3dp{0, 1, -1, 1});

        auto l3 = rcmt(l1, l2);
        CHECK(l3 == y_axis_3dp);

        fmt::println("l1 = {}", l1);
        fmt::println("l2 = {}", l2);
        fmt::println("l3 = {}", l3);

        fmt::println("");
    }
}
