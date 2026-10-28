#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// THE EQUIVARIANCE PATTERN for a quantity read off motors.
//
// A motor composes correctly for any number of turns; a NUMBER read off it -- an angle,
// a heading, a controller's error -- goes through the logarithm, a chart that is one-to-
// one only near the identity (ga_docu.pdf, "Reading an angle off a motor"). A test that
// runs at heading zero sees none of the chart's failure classes: the wrap and the double
// cover, the non-additivity of logarithms, the ill-conditioning near a half turn. So the
// WHOLE scene is turned -- every motor M of it becomes G (x) M -- by headings that cover
// them: 0, 90, 180 and 270 degrees, a full turn plus a little and nearly two, about a
// vertical axis through a point off the origin. A reading then has to be
//
//   - INVARIANT if it is expressed in a frame of the scene (a body-frame error), and
//   - TURNED WITH G, move(reading, G), if it is expressed in world axes.
//
// The checkers return the largest deviation over the headings, so a gate states its own
// tolerance and can print the number. A reading here is any callable taking the scene's
// motors (a std::vector) and returning a twist.
/////////////////////////////////////////////////////////////////////////////////////////

#include "ga/ga_pga.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace equivariance {

using namespace hd::ga;
using namespace hd::ga::pga;

// the headings [deg] at which the chart's failure classes show
inline std::vector<value_t> const& headings_deg()
{
    static std::vector<value_t> const h{0.0, 90.0, 180.0, 270.0, 365.0, 719.0};
    return h;
}

// a turn by `deg` about a vertical axis through a point off the origin (3D), about a
// pivot off the origin (2D)
inline mvec3dp_e turn3dp(value_t deg)
{
    return get_motor(
        unitize(wdg(vec3dp{0.7, -0.4, 0.0, 1.0}, vec3dp{0.7, -0.4, 1.0, 1.0})),
        deg2rad(deg));
}
inline mvec2dp_u turn2dp(value_t deg)
{
    return get_motor(vec2dp{0.7, -0.4, 1.0}, deg2rad(deg));
}

inline value_t max_diff(bivec3dp const& a, bivec3dp const& b)
{
    return std::max({std::abs(a.vx - b.vx), std::abs(a.vy - b.vy), std::abs(a.vz - b.vz),
                     std::abs(a.mx - b.mx), std::abs(a.my - b.my),
                     std::abs(a.mz - b.mz)});
}
inline value_t max_diff(vec2dp const& a, vec2dp const& b)
{
    return std::max({std::abs(a.x - b.x), std::abs(a.y - b.y), std::abs(a.z - b.z)});
}

// the largest change of a BODY-frame reading when the scene is turned
template <class F>
value_t body_invariance3dp(F&& reading, std::vector<mvec3dp_e> const& scene)
{
    auto const r0 = reading(scene);
    value_t worst = 0.0;
    for (value_t const h : headings_deg()) {
        auto const G = turn3dp(h);
        std::vector<mvec3dp_e> turned;
        for (auto const& M : scene)
            turned.push_back(rgpr(G, M));
        worst = std::max(worst, max_diff(reading(turned), r0));
    }
    return worst;
}

// the largest deviation of a WORLD-frame reading from the reading turned with the scene
template <class F>
value_t world_covariance3dp(F&& reading, std::vector<mvec3dp_e> const& scene)
{
    auto const r0 = reading(scene);
    value_t worst = 0.0;
    for (value_t const h : headings_deg()) {
        auto const G = turn3dp(h);
        std::vector<mvec3dp_e> turned;
        for (auto const& M : scene)
            turned.push_back(rgpr(G, M));
        worst = std::max(worst, max_diff(reading(turned), move3dp(r0, G)));
    }
    return worst;
}

template <class F>
value_t body_invariance2dp(F&& reading, std::vector<mvec2dp_u> const& scene)
{
    auto const r0 = reading(scene);
    value_t worst = 0.0;
    for (value_t const h : headings_deg()) {
        auto const G = turn2dp(h);
        std::vector<mvec2dp_u> turned;
        for (auto const& M : scene)
            turned.push_back(rgpr(G, M));
        worst = std::max(worst, max_diff(reading(turned), r0));
    }
    return worst;
}

template <class F>
value_t world_covariance2dp(F&& reading, std::vector<mvec2dp_u> const& scene)
{
    auto const r0 = reading(scene);
    value_t worst = 0.0;
    for (value_t const h : headings_deg()) {
        auto const G = turn2dp(h);
        std::vector<mvec2dp_u> turned;
        for (auto const& M : scene)
            turned.push_back(rgpr(G, M));
        worst = std::max(worst, max_diff(reading(turned), move2dp(r0, G)));
    }
    return worst;
}

} // namespace equivariance
