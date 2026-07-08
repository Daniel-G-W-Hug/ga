// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <algorithm> // std::min, std::max, std::clamp
#include <array>     // std::array (Cai HTM 4x4 reference matrices, Phase F)
#include <cmath>     // std::sin, std::cos, std::atan
#include <limits>    // std::numeric_limits
#include <random>    // std::mt19937, distributions (grain-size randomness)
#include <utility>   // std::pair
#include <vector>    // std::vector

// include functions to be tested
#include "ga/ga_ega.hpp" // for cross product
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)
using namespace hd::ga::ega; // use specific operations of EGA (for cross product)


/////////////////////////////////////////////////////////////////////////////////////////
// projective geometric algebra pga3dp application tests
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: application tests")
{

    TEST_CASE("pga3dp: line perpendicular to non-intersecting lines (rcmt)")
    {

        fmt::println("");
        fmt::println("pga3dp: line perpendicular to non-intersecting lines (rcmt)");
        fmt::println("");

        auto l1 = x_axis_3dp;
        auto l2 = wdg(vec3dp{0, 1, 0, 1}, vec3dp{0, 1, -1, 1});

        auto l3 = rcmt(l2, rrev(l1));
        CHECK(l3 == y_axis_3dp);

        fmt::println("l1 = {}", l1);
        fmt::println("l2 = {}", l2);
        fmt::println("l3 = {}", l3);
        fmt::println("");

        // use the opportunity to test the sign function
        fmt::println("sign(val) yields +1.0 for val >= 0.0 and -1.0 for val < 0.0).");
        fmt::println("val==0 is always regarded as positive sign:");
        fmt::println("sign(l2.vx) = {} (postive zero)", sign(l2.vx));
        fmt::println("sign(l2.vz) = {} (negative one)", sign(l2.vz));
        fmt::println("sign(l3.vx) = {} (negative zero)", sign(l3.vx));
        fmt::println("sign(l3.vy) = {} (postive one)", sign(l3.vy));
        CHECK(sign(l2.vx) == 1.0);
        CHECK(sign(l2.vz) == -1.0);
        CHECK(sign(l3.vx) == 1.0);
        CHECK(sign(l3.vy) == 1.0);
        fmt::println("");
    }

    TEST_CASE("pga: bivector rate of change")
    {

        fmt::println("");
        fmt::println("pga: bivector rate of change - translation:");

        auto v2dp = vec2dp{1, 2, 0};    // translation vector 2d
        auto v3dp = vec3dp{1, 2, 3, 0}; // translation vector 3d

        auto omega_tra_2dp = att(bulk_dual(v2dp));   // rate of change is a vec (2D)
        auto omega_tra_3dp = att(r_bulk_dual(v3dp)); // rate of change is a bivec (3D)

        CHECK(omega_tra_2dp == -cmpl(wdg(e3_2dp, v2dp)));   // -cmpl(bivec)  = vec
        CHECK(omega_tra_3dp == -r_cmpl(wdg(e4_3dp, v3dp))); // -r_cmpl(bivec) = bivec

        fmt::println("v2dp          = {}", v2dp);
        fmt::println("v3dp          = {}", v3dp);
        fmt::println("omega_tra_2dp = {}", omega_tra_2dp);
        fmt::println("omega_tra_3dp = {}", omega_tra_3dp);
        fmt::println("");
        // fmt::println("att(bulk_dual(v2dp)) = {}", att(bulk_dual(v2dp)));
        // fmt::println("");

        fmt::println("pga: bivector rate of change - rotation:");

        auto omeg2dp = scalar2dp{3};       // omega as angular speed which is a scalar
        auto omeg3dp = vec3dp{1, 2, 3, 0}; // rotation rate omega as vec

        auto omega_rot_2dp = wdg(e3_2dp, omeg2dp); // rate of change is a vec (2D)
        auto omega_rot_3dp = wdg(e4_3dp, omeg3dp); // rate of change is a bivec (3D)

        CHECK(omega_rot_2dp == wdg(e3_2dp, weight_dual(bulk_dual(omeg2dp))));
        CHECK(omega_rot_3dp == wdg(e4_3dp, l_weight_dual(r_bulk_dual(omeg3dp))));

        fmt::println("omeg2dp       = {}", omeg2dp);
        fmt::println("omeg3dp       = {}", omeg3dp);
        fmt::println("omega_rot_2dp = {}", omega_rot_2dp);
        // fmt::println("plr_2dp       = {}", bulk_dual(omeg2dp));
        fmt::println("omega_rot_3dp = {}", omega_rot_3dp);
        // fmt::println("plr_3dp       = {}", r_bulk_dual(omeg3dp));
        fmt::println("");

        fmt::println("pga: bivector rate of change - checking translational speed.");

        auto X1_2dp = vec2dp{1, 0, 1};
        auto X2_2dp = vec2dp{2, 0, 1};
        auto X3_2dp = vec2dp{0, 1, 1};

        CHECK(rcmt(omega_tra_2dp, X1_2dp) == v2dp);
        CHECK(rcmt(omega_tra_2dp, X2_2dp) == v2dp);
        CHECK(rcmt(omega_tra_2dp, X3_2dp) == v2dp);

        auto X1_3dp = vec3dp{1, 0, 0, 1};
        auto X2_3dp = vec3dp{2, 0, 0, 1};
        auto X3_3dp = vec3dp{0, 1, 0, 1};

        CHECK(rcmt(omega_tra_3dp, X1_3dp) == v3dp);
        CHECK(rcmt(omega_tra_3dp, X2_3dp) == v3dp);
        CHECK(rcmt(omega_tra_3dp, X3_3dp) == v3dp);

        fmt::println("pga: bivector rate of change - checking rotational speed.");

        CHECK(bulk_nrm(rcmt(omega_rot_2dp, X1_2dp)) == omeg2dp * bulk_nrm(X1_2dp));
        CHECK(bulk_nrm(rcmt(omega_rot_2dp, X2_2dp)) == omeg2dp * bulk_nrm(X2_2dp));
        CHECK(bulk_nrm(rcmt(omega_rot_2dp, X3_2dp)) == omeg2dp * bulk_nrm(X3_2dp));

        auto L = wdg(O_3dp, omeg3dp); // rotational axis
        CHECK(bulk_nrm(rcmt(omega_rot_3dp, X1_3dp)) ==
              bulk_nrm(omeg3dp) * bulk_nrm(ortho_proj3dp(X1_3dp, L) - X1_3dp));
        CHECK(bulk_nrm(rcmt(omega_rot_3dp, X2_3dp)) ==
              bulk_nrm(omeg3dp) * bulk_nrm(ortho_proj3dp(X2_3dp, L) - X2_3dp));
        CHECK(bulk_nrm(rcmt(omega_rot_3dp, X3_3dp)) ==
              bulk_nrm(omeg3dp) * bulk_nrm(ortho_proj3dp(X3_3dp, L) - X3_3dp));

        fmt::println("");
    }
}
