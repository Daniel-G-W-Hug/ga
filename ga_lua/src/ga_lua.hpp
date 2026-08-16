#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// GCC compatibility: Include <climits> before Lua headers to ensure LLONG_MAX is defined
#include <climits>

#define SOL_ALL_SAFETIES_ON 1
#include "sol/sol.hpp"

#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"
#include "ga/ga_sta.hpp"

////////////////////////////////////////////////////////////////////////////////
// register basic types, geometric operations and constants
// of user-defined types with lua
////////////////////////////////////////////////////////////////////////////////
// EGA type registration functions
void register_2d_types(sol::state& lua);
void register_3d_types(sol::state& lua);
// PGA type registration functions
void register_2dp_types(sol::state& lua);
void register_3dp_types(sol::state& lua);
// STA type registration (G(1,3,0): scalar4ds ... mvec4ds, dualnum4ds)
void register_4ds_types(sol::state& lua);
// PGA geometric convenience types (point / vector / line / plane subclasses)
void register_convenience_types(sol::state& lua);
// PGA physics pure-data structs (pose / kin_state / joint_state / loop_constraint)
// and their scoped enums (joint / constraint)
void register_physics_pods(sol::state& lua);
// Function and constant registration
void register_functions(sol::state& lua);
void register_constants(sol::state& lua);
// Lua-level forwarders (templates that cannot be sol::resolve-bound)
void register_forwarders(sol::state& lua);


////////////////////////////////////////////////////////////////////////////////
// 2d types
////////////////////////////////////////////////////////////////////////////////
void register_2d_types(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::ega;

    lua.new_usertype<scalar2d>(
        "scalar2d",
        sol::constructors<scalar2d(), scalar2d(value_t const&), scalar2d(value_t&&)>(),
        "copy", [](const scalar2d& s) { return scalar2d(s); },
        sol::meta_function::to_string,
        [](const scalar2d& s) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Scalar2d({})", suppress_negative_zero(double(s)));
        },
        sol::meta_function::unary_minus, sol::resolve<scalar2d(scalar2d)>(operator-),
        sol::meta_function::addition,
        sol::overload(sol::resolve<scalar2d(scalar2d, scalar2d)>(operator+),
                      sol::resolve<mvec2d_e(scalar2d, pscalar2d)>(operator+),
                      sol::resolve<mvec2d_e(scalar2d, mvec2d_e const&)>(operator+),
                      sol::resolve<mvec2d(scalar2d, vec2d const&)>(operator+),
                      sol::resolve<mvec2d(pscalar2d, vec2d const&)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<scalar2d(scalar2d, scalar2d)>(operator-),
                      sol::resolve<mvec2d_e(scalar2d, pscalar2d)>(operator-),
                      sol::resolve<mvec2d_e(scalar2d, mvec2d_e const&)>(operator-),
                      sol::resolve<mvec2d(scalar2d, vec2d const&)>(operator-),
                      sol::resolve<mvec2d(pscalar2d, vec2d const&)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<scalar2d(scalar2d, scalar2d)>(wdg),
                      sol::resolve<vec2d(scalar2d, vec2d const&)>(wdg),
                      sol::resolve<pscalar2d(scalar2d, pscalar2d)>(wdg),
                      sol::resolve<mvec2d_e(scalar2d, mvec2d_e const&)>(wdg),
                      sol::resolve<mvec2d(scalar2d, mvec2d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<scalar2d(scalar2d, scalar2d)>(operator*),
                      sol::resolve<vec2d(scalar2d, vec2d const&)>(operator*),
                      sol::resolve<pscalar2d(scalar2d, pscalar2d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<scalar2d(scalar2d, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::overload(sol::resolve<vec2d(scalar2d, vec2d const&)>(operator<<),
                      sol::resolve<pscalar2d(scalar2d, pscalar2d)>(operator<<),
                      sol::resolve<mvec2d_e(scalar2d, mvec2d_e const&)>(operator<<),
                      sol::resolve<mvec2d(scalar2d, mvec2d const&)>(operator<<)),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<scalar2d(scalar2d, pscalar2d)>(operator>>),
                      sol::resolve<scalar2d(scalar2d, vec2d const&)>(operator>>)));


    lua.new_usertype<vec2d>(
        "vec2d",
        sol::constructors<vec2d(), vec2d(value_t, value_t), vec2d(vec2d const&),
                          vec2d(vec2d&&)>(),
        "copy", [](const vec2d& v) { return vec2d(v); },
        // component access
        "x", &vec2d::x, "y", &vec2d::y, sol::meta_function::to_string,
        [](const vec2d& v) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Vec2d({},{})", suppress_negative_zero(v.x),
                               suppress_negative_zero(v.y));
        },
        sol::meta_function::unary_minus, sol::resolve<vec2d(vec2d const&)>(operator-),
        sol::meta_function::addition,
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(operator+),
                      sol::resolve<mvec2d(vec2d const&, scalar2d)>(operator+),
                      sol::resolve<mvec2d(vec2d const&, pscalar2d)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(operator-),
                      sol::resolve<mvec2d(vec2d const&, scalar2d)>(operator-),
                      sol::resolve<mvec2d(vec2d const&, pscalar2d)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<vec2d(vec2d const&, scalar2d)>(wdg),
                      sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(wdg),
                      sol::resolve<scalar2d(vec2d const&, pscalar2d)>(wdg),
                      sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(wdg),
                      sol::resolve<mvec2d(vec2d const&, mvec2d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec2d(vec2d const&, value_t)>(operator*),
                      sol::resolve<vec2d(value_t, vec2d const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(operator*),
                      sol::resolve<mvec2d(vec2d const&, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d_e(vec2d const&, vec2d const&)>(operator*),
                      sol::resolve<vec2d(vec2d const&, pscalar2d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec2d(vec2d const&, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::overload(sol::resolve<scalar2d(vec2d const&, vec2d const&)>(operator<<),
                      sol::resolve<vec2d(vec2d const&, pscalar2d)>(operator<<),
                      sol::resolve<scalar2d(vec2d const&, scalar2d)>(operator<<),
                      sol::resolve<mvec2d(vec2d const&, mvec2d const&)>(operator<<)),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<vec2d(vec2d const&, scalar2d)>(operator>>),
                      sol::resolve<scalar2d(vec2d const&, vec2d const&)>(operator>>),
                      sol::resolve<scalar2d(vec2d const&, pscalar2d)>(operator>>)));


    lua.new_usertype<pscalar2d>(
        "pscalar2d",
        sol::constructors<pscalar2d(), pscalar2d(value_t const&), pscalar2d(value_t&&)>(),
        "copy", [](const pscalar2d& ps) { return pscalar2d(ps); },
        sol::meta_function::to_string,
        [](const pscalar2d& ps) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("PScalar2d({})", suppress_negative_zero(double(ps)));
        },
        sol::meta_function::unary_minus, sol::resolve<pscalar2d(pscalar2d)>(operator-),
        sol::meta_function::addition,
        sol::overload(sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(operator+),
                      sol::resolve<mvec2d_e(pscalar2d, scalar2d)>(operator+),
                      sol::resolve<mvec2d_e(pscalar2d, mvec2d_e const&)>(operator+),
                      sol::resolve<mvec2d(pscalar2d, vec2d const&)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(operator-),
                      sol::resolve<mvec2d_e(pscalar2d, scalar2d)>(operator-),
                      sol::resolve<mvec2d_e(pscalar2d, mvec2d_e const&)>(operator-),
                      sol::resolve<mvec2d(pscalar2d, vec2d const&)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<pscalar2d(pscalar2d, scalar2d)>(wdg),
                      sol::resolve<scalar2d(pscalar2d, vec2d const&)>(wdg),
                      sol::resolve<scalar2d(pscalar2d, pscalar2d)>(wdg),
                      sol::resolve<pscalar2d(pscalar2d, mvec2d_e const&)>(wdg),
                      sol::resolve<pscalar2d(pscalar2d, mvec2d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<pscalar2d(pscalar2d, value_t)>(operator*),
                      sol::resolve<pscalar2d(value_t, pscalar2d)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<mvec2d(pscalar2d, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d_e(pscalar2d, mvec2d_e const&)>(operator*),
                      sol::resolve<vec2d(pscalar2d, vec2d const&)>(operator*),
                      sol::resolve<pscalar2d(pscalar2d, scalar2d)>(operator*),
                      sol::resolve<scalar2d(pscalar2d, pscalar2d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar2d(pscalar2d, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::overload(sol::resolve<scalar2d(pscalar2d, scalar2d)>(operator<<),
                      sol::resolve<scalar2d(pscalar2d, vec2d const&)>(operator<<)),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<vec2d(pscalar2d, vec2d const&)>(operator>>),
                      sol::resolve<pscalar2d(pscalar2d, scalar2d)>(operator>>)));


    lua.new_usertype<mvec2d_e>(
        "mvec2d_e",
        sol::constructors<mvec2d_e(), mvec2d_e(value_t, value_t),
                          mvec2d_e(mvec2d_e const&), mvec2d_e(mvec2d_e&&),
                          mvec2d_e(scalar2d), mvec2d_e(pscalar2d),
                          mvec2d_e(scalar2d, pscalar2d)>(),
        "copy", [](const mvec2d_e& obj) { return mvec2d_e(obj); },
        // component access
        "c0", &mvec2d_e::c0, "c1", &mvec2d_e::c1, sol::meta_function::to_string,
        [](const mvec2d_e& mv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("MVec2d_E({},{})", suppress_negative_zero(mv.c0),
                               suppress_negative_zero(mv.c1));
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec2d_e(mvec2d_e const&)>(operator-), sol::meta_function::addition,
        sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator+),
                      sol::resolve<mvec2d_e(mvec2d_e const&, scalar2d)>(operator+),
                      sol::resolve<mvec2d_e(mvec2d_e const&, pscalar2d)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator-),
                      sol::resolve<mvec2d_e(mvec2d_e const&, scalar2d)>(operator-),
                      sol::resolve<mvec2d_e(mvec2d_e const&, pscalar2d)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&, scalar2d)>(wdg),
                      sol::resolve<vec2d(mvec2d_e const&, vec2d const&)>(wdg),
                      sol::resolve<pscalar2d(mvec2d_e const&, pscalar2d)>(wdg),
                      sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(wdg),
                      sol::resolve<mvec2d(mvec2d_e const&, mvec2d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&, value_t)>(operator*),
                      sol::resolve<mvec2d_e(value_t, mvec2d_e const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<vec2d(mvec2d_e const&, vec2d const&)>(operator*),
                      sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator*),
                      sol::resolve<mvec2d_e(mvec2d_e const&, pscalar2d)>(operator*),
                      sol::resolve<mvec2d(mvec2d_e const&, mvec2d const&)>(operator*),
                      sol::resolve<vec2d(mvec2d_e const&, vec2d const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2d_e(mvec2d_e const&, value_t)>(operator/),
        sol::meta_function::bitwise_right_shift,
        sol::resolve<mvec2d_e(mvec2d_e const&, scalar2d)>(operator>>));


    lua.new_usertype<mvec2d>(
        "mvec2d",
        sol::constructors<mvec2d(), mvec2d(value_t, value_t, value_t, value_t),
                          mvec2d(mvec2d const&), mvec2d(mvec2d&&), mvec2d(scalar2d),
                          mvec2d(vec2d const&), mvec2d(pscalar2d),
                          mvec2d(scalar2d, pscalar2d), mvec2d(mvec2d_e const&),
                          mvec2d(scalar2d, vec2d const&, pscalar2d)>(),
        "copy", [](const mvec2d& obj) { return mvec2d(obj); },
        // component access
        "c0", &mvec2d::c0, "c1", &mvec2d::c1, "c2", &mvec2d::c2, "c3", &mvec2d::c3,
        sol::meta_function::to_string,
        [](const mvec2d& mv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("MVec2d({},{},{},{})", suppress_negative_zero(mv.c0),
                               suppress_negative_zero(mv.c1),
                               suppress_negative_zero(mv.c2),
                               suppress_negative_zero(mv.c3));
        },
        sol::meta_function::unary_minus, sol::resolve<mvec2d(mvec2d const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator-),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(wdg),
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(wdg)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<mvec2d(mvec2d const&, scalar2d)>(wdg),
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(wdg),
                      sol::resolve<pscalar2d(mvec2d const&, pscalar2d)>(wdg),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(wdg),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec2d(mvec2d const&, value_t)>(operator*),
                      sol::resolve<mvec2d(value_t, mvec2d const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, pscalar2d)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2d(mvec2d const&, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator<<),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<mvec2d(mvec2d const&, scalar2d)>(operator>>),
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(operator>>),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator>>)));
}


////////////////////////////////////////////////////////////////////////////////
// 3d types
////////////////////////////////////////////////////////////////////////////////
void register_3d_types(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::ega;

    lua.new_usertype<scalar3d>(
        "scalar3d",
        sol::constructors<scalar3d(), scalar3d(value_t const&), scalar3d(value_t&&)>(),
        "copy", [](const scalar3d& obj) { return scalar3d(obj); },
        sol::meta_function::to_string,
        [](const scalar3d& s) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Scalar3d({})", suppress_negative_zero(double(s)));
        },
        sol::meta_function::unary_minus, sol::resolve<scalar3d(scalar3d)>(operator-),
        sol::meta_function::addition,
        sol::overload(sol::resolve<scalar3d(scalar3d, scalar3d)>(operator+),
                      sol::resolve<mvec3d_e(scalar3d, bivec3d const&)>(operator+),
                      sol::resolve<mvec3d_e(scalar3d, mvec3d_e const&)>(operator+),
                      sol::resolve<mvec3d(scalar3d, vec3d const&)>(operator+),
                      sol::resolve<mvec3d(scalar3d, pscalar3d)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<scalar3d(scalar3d, scalar3d)>(operator-),
                      sol::resolve<mvec3d_e(scalar3d, bivec3d const&)>(operator-),
                      sol::resolve<mvec3d_e(scalar3d, mvec3d_e const&)>(operator-),
                      sol::resolve<mvec3d(scalar3d, vec3d const&)>(operator-),
                      sol::resolve<mvec3d(scalar3d, pscalar3d)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<scalar3d(scalar3d, scalar3d)>(wdg),
                      sol::resolve<vec3d(scalar3d, vec3d const&)>(wdg),
                      sol::resolve<bivec3d(scalar3d, bivec3d const&)>(wdg),
                      sol::resolve<pscalar3d(scalar3d, pscalar3d)>(wdg),
                      sol::resolve<mvec3d(scalar3d, mvec3d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<scalar3d(scalar3d, scalar3d)>(operator*),
                      sol::resolve<vec3d(scalar3d, vec3d const&)>(operator*),
                      sol::resolve<bivec3d(scalar3d, bivec3d const&)>(operator*),
                      sol::resolve<pscalar3d(scalar3d, pscalar3d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<scalar3d(scalar3d, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::overload(sol::resolve<scalar3d(scalar3d, scalar3d)>(operator<<),
                      sol::resolve<vec3d(scalar3d, vec3d const&)>(operator<<),
                      sol::resolve<bivec3d(scalar3d, bivec3d const&)>(operator<<),
                      sol::resolve<pscalar3d(scalar3d, pscalar3d)>(operator<<),
                      sol::resolve<mvec3d_e(scalar3d, mvec3d_e const&)>(operator<<),
                      sol::resolve<mvec3d_u(scalar3d, mvec3d_u const&)>(operator<<),
                      sol::resolve<mvec3d(scalar3d, mvec3d const&)>(operator<<)),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<scalar3d(scalar3d, pscalar3d)>(operator>>),
                      sol::resolve<scalar3d(scalar3d, bivec3d const&)>(operator>>),
                      sol::resolve<scalar3d(scalar3d, vec3d const&)>(operator>>)));


    lua.new_usertype<vec3d>(
        "vec3d",
        sol::constructors<vec3d(), vec3d(value_t, value_t, value_t), vec3d(vec3d const&),
                          vec3d(vec3d&&)>(),
        "copy", [](const vec3d& obj) { return vec3d(obj); },
        // component access
        "x", &vec3d::x, "y", &vec3d::y, "z", &vec3d::z, sol::meta_function::to_string,
        [](const vec3d& v) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Vec3d({},{},{})", suppress_negative_zero(v.x),
                               suppress_negative_zero(v.y), suppress_negative_zero(v.z));
        },
        sol::meta_function::unary_minus, sol::resolve<vec3d(vec3d const&)>(operator-),
        sol::meta_function::addition,
        sol::overload(sol::resolve<vec3d(vec3d const&, vec3d const&)>(operator+),
                      sol::resolve<mvec3d_u(vec3d const&, pscalar3d)>(operator+),
                      sol::resolve<mvec3d_u(vec3d const&, mvec3d_u const&)>(operator+),
                      sol::resolve<mvec3d(vec3d const&, scalar3d)>(operator+),
                      sol::resolve<mvec3d(vec3d const&, bivec3d const&)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<vec3d(vec3d const&, vec3d const&)>(operator-),
                      sol::resolve<mvec3d_u(vec3d const&, pscalar3d)>(operator-),
                      sol::resolve<mvec3d_u(vec3d const&, mvec3d_u const&)>(operator-),
                      sol::resolve<mvec3d(vec3d const&, scalar3d)>(operator-),
                      sol::resolve<mvec3d(vec3d const&, bivec3d const&)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<vec3d(vec3d const&, scalar3d)>(wdg),
                      sol::resolve<bivec3d(vec3d const&, vec3d const&)>(wdg),
                      sol::resolve<pscalar3d(vec3d const&, bivec3d const&)>(wdg),
                      sol::resolve<scalar3d(vec3d const&, pscalar3d)>(wdg),
                      sol::resolve<mvec3d(vec3d const&, mvec3d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec3d(vec3d const&, value_t)>(operator*),
                      sol::resolve<vec3d(value_t, vec3d const&)>(operator*),
                      sol::resolve<mvec3d_e(vec3d const&, vec3d const&)>(operator*),
                      sol::resolve<mvec3d_u(vec3d const&, bivec3d const&)>(operator*),
                      sol::resolve<bivec3d(vec3d const&, pscalar3d)>(operator*),
                      sol::resolve<mvec3d(vec3d const&, mvec3d const&)>(operator*),
                      sol::resolve<mvec3d_u(vec3d const&, mvec3d_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec3d(vec3d const&, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::overload(sol::resolve<scalar3d(vec3d const&, vec3d const&)>(operator<<),
                      sol::resolve<vec3d(vec3d const&, bivec3d const&)>(operator<<),
                      sol::resolve<bivec3d(vec3d const&, pscalar3d)>(operator<<),
                      sol::resolve<scalar3d(vec3d const&, scalar3d)>(operator<<),
                      sol::resolve<mvec3d(vec3d const&, mvec3d const&)>(operator<<)),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<vec3d(vec3d const&, scalar3d)>(operator>>),
                      sol::resolve<scalar3d(vec3d const&, vec3d const&)>(operator>>),
                      sol::resolve<scalar3d(vec3d const&, pscalar3d)>(operator>>),
                      sol::resolve<scalar3d(vec3d const&, bivec3d const&)>(operator>>)));


    lua.new_usertype<bivec3d>(
        "bivec3d",
        sol::constructors<bivec3d(), bivec3d(value_t, value_t, value_t),
                          bivec3d(bivec3d const&), bivec3d(bivec3d&&)>(),
        "copy", [](const bivec3d& obj) { return bivec3d(obj); },
        // component access
        "x", &bivec3d::x, "y", &bivec3d::y, "z", &bivec3d::z,
        sol::meta_function::to_string,
        [](const bivec3d& bv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("BiVec3d({},{},{})", suppress_negative_zero(bv.x),
                               suppress_negative_zero(bv.y),
                               suppress_negative_zero(bv.z));
        },
        sol::meta_function::unary_minus, sol::resolve<bivec3d(bivec3d const&)>(operator-),
        sol::meta_function::addition,
        sol::overload(sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(operator+),
                      sol::resolve<mvec3d_e(bivec3d const&, scalar3d)>(operator+),
                      sol::resolve<mvec3d_e(bivec3d const&, mvec3d_e const&)>(operator+),
                      sol::resolve<mvec3d(bivec3d const&, vec3d const&)>(operator+),
                      sol::resolve<mvec3d(bivec3d const&, pscalar3d)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(operator-),
                      sol::resolve<mvec3d_e(bivec3d const&, scalar3d)>(operator-),
                      sol::resolve<mvec3d_e(bivec3d const&, mvec3d_e const&)>(operator-),
                      sol::resolve<mvec3d(bivec3d const&, vec3d const&)>(operator-),
                      sol::resolve<mvec3d(bivec3d const&, pscalar3d)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<bivec3d(bivec3d const&, scalar3d)>(wdg),
                      sol::resolve<pscalar3d(bivec3d const&, vec3d const&)>(wdg),
                      sol::resolve<scalar3d(bivec3d const&, bivec3d const&)>(wdg),
                      sol::resolve<scalar3d(bivec3d const&, pscalar3d)>(wdg),
                      sol::resolve<mvec3d(bivec3d const&, mvec3d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<bivec3d(bivec3d const&, value_t)>(operator*),
                      sol::resolve<bivec3d(value_t, bivec3d const&)>(operator*),
                      sol::resolve<mvec3d_u(bivec3d const&, vec3d const&)>(operator*),
                      sol::resolve<mvec3d_e(bivec3d const&, bivec3d const&)>(operator*),
                      sol::resolve<vec3d(bivec3d const&, pscalar3d)>(operator*),
                      sol::resolve<mvec3d_e(bivec3d const&, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d(bivec3d const&, mvec3d const&)>(operator*),
                      sol::resolve<mvec3d_u(bivec3d const&, mvec3d_u const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<bivec3d(bivec3d const&, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::overload(sol::resolve<scalar3d(bivec3d const&, bivec3d const&)>(operator<<),
                      sol::resolve<vec3d(bivec3d const&, pscalar3d)>(operator<<),
                      sol::resolve<scalar3d(bivec3d const&, scalar3d)>(operator<<),
                      sol::resolve<scalar3d(bivec3d const&, vec3d const&)>(operator<<),
                      sol::resolve<mvec3d(bivec3d const&, mvec3d const&)>(operator<<)),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<bivec3d(bivec3d const&, scalar3d)>(operator>>),
                      sol::resolve<scalar3d(bivec3d const&, bivec3d const&)>(operator>>),
                      sol::resolve<vec3d(bivec3d const&, vec3d const&)>(operator>>),
                      sol::resolve<scalar3d(bivec3d const&, pscalar3d)>(operator>>)));

    lua.new_usertype<pscalar3d>(
        "pscalar3d",
        sol::constructors<pscalar3d(), pscalar3d(value_t const&), pscalar3d(value_t&&)>(),
        "copy", [](const pscalar3d& obj) { return pscalar3d(obj); },
        sol::meta_function::to_string,
        [](const pscalar3d& ps) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("PScalar3d({})", suppress_negative_zero(double(ps)));
        },
        sol::meta_function::unary_minus, sol::resolve<pscalar3d(pscalar3d)>(operator-),
        sol::meta_function::addition,
        sol::overload(sol::resolve<pscalar3d(pscalar3d, pscalar3d)>(operator+),
                      sol::resolve<mvec3d_u(pscalar3d, vec3d const&)>(operator+),
                      sol::resolve<mvec3d_u(pscalar3d, mvec3d_u const&)>(operator+),
                      sol::resolve<mvec3d(pscalar3d, scalar3d)>(operator+),
                      sol::resolve<mvec3d(pscalar3d, bivec3d const&)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<pscalar3d(pscalar3d, pscalar3d)>(operator-),
                      sol::resolve<mvec3d_u(pscalar3d, vec3d const&)>(operator-),
                      sol::resolve<mvec3d_u(pscalar3d, mvec3d_u const&)>(operator-),
                      sol::resolve<mvec3d(pscalar3d, scalar3d)>(operator-),
                      sol::resolve<mvec3d(pscalar3d, bivec3d const&)>(operator-)),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<pscalar3d(pscalar3d, scalar3d)>(wdg),
                      sol::resolve<scalar3d(pscalar3d, vec3d const&)>(wdg),
                      sol::resolve<scalar3d(pscalar3d, bivec3d const&)>(wdg),
                      sol::resolve<scalar3d(pscalar3d, pscalar3d)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<pscalar3d(pscalar3d, value_t)>(operator*),
                      sol::resolve<pscalar3d(value_t, pscalar3d)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<bivec3d(pscalar3d, vec3d const&)>(operator*),
                      sol::resolve<vec3d(pscalar3d, bivec3d const&)>(operator*),
                      sol::resolve<mvec3d_e(pscalar3d, mvec3d_u const&)>(operator*),
                      sol::resolve<mvec3d_u(pscalar3d, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d(pscalar3d, mvec3d const&)>(operator*),
                      sol::resolve<scalar3d(pscalar3d, pscalar3d)>(operator*),
                      sol::resolve<pscalar3d(pscalar3d, scalar3d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar3d(pscalar3d, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::overload(sol::resolve<scalar3d(pscalar3d, scalar3d)>(operator<<),
                      sol::resolve<scalar3d(pscalar3d, vec3d const&)>(operator<<),
                      sol::resolve<scalar3d(pscalar3d, bivec3d const&)>(operator<<)),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<pscalar3d(pscalar3d, scalar3d)>(operator>>),
                      sol::resolve<bivec3d(pscalar3d, vec3d const&)>(operator>>),
                      sol::resolve<vec3d(pscalar3d, bivec3d const&)>(operator>>)));


    lua.new_usertype<mvec3d_e>(
        "mvec3d_e",
        sol::constructors<mvec3d_e(), mvec3d_e(value_t, value_t, value_t, value_t),
                          mvec3d_e(mvec3d_e const&), mvec3d_e(mvec3d_e&&),
                          mvec3d_e(scalar3d), mvec3d_e(bivec3d const&),
                          mvec3d_e(scalar3d, bivec3d const&)>(),
        "copy", [](const mvec3d_e& obj) { return mvec3d_e(obj); },
        // component access
        "c0", &mvec3d_e::c0, "c1", &mvec3d_e::c1, "c2", &mvec3d_e::c2, "c3",
        &mvec3d_e::c3, sol::meta_function::to_string,
        [](const mvec3d_e& mv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("MVec3d_E({},{},{},{})", suppress_negative_zero(mv.c0),
                               suppress_negative_zero(mv.c1),
                               suppress_negative_zero(mv.c2),
                               suppress_negative_zero(mv.c3));
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec3d_e(mvec3d_e const&)>(operator-), sol::meta_function::addition,
        sol::overload(sol::resolve<mvec3d_e(mvec3d_e const&, mvec3d_e const&)>(operator+),
                      sol::resolve<mvec3d_e(mvec3d_e const&, scalar3d)>(operator+),
                      sol::resolve<mvec3d_e(mvec3d_e const&, bivec3d const&)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<mvec3d_e(mvec3d_e const&, mvec3d_e const&)>(operator-),
                      sol::resolve<mvec3d_e(mvec3d_e const&, scalar3d)>(operator-),
                      sol::resolve<mvec3d_e(mvec3d_e const&, bivec3d const&)>(operator-)),
        sol::meta_function::multiplication,
        sol::overload(
            sol::resolve<mvec3d_e(mvec3d_e const&, value_t)>(operator*),
            sol::resolve<mvec3d_e(value_t, mvec3d_e const&)>(operator*),
            sol::resolve<mvec3d_u(mvec3d_e const&, pscalar3d)>(operator*),
            sol::resolve<mvec3d_u(mvec3d_e const&, vec3d const&)>(operator*),
            sol::resolve<mvec3d_e(mvec3d_e const&, bivec3d const&)>(operator*),
            sol::resolve<mvec3d_e(mvec3d_e const&, mvec3d_e const&)>(operator*),
            sol::resolve<mvec3d(mvec3d_e const&, mvec3d const&)>(operator*),
            sol::resolve<mvec3d_u(mvec3d_e const&, mvec3d_u const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3d_e(mvec3d_e const&, value_t)>(operator/),
        sol::meta_function::bitwise_right_shift,
        sol::resolve<mvec3d_e(mvec3d_e const&, scalar3d)>(operator>>));


    lua.new_usertype<mvec3d_u>(
        "mvec3d_u",
        sol::constructors<mvec3d_u(), mvec3d_u(value_t, value_t, value_t, value_t),
                          mvec3d_u(mvec3d_u const&), mvec3d_u(mvec3d_u&&),
                          mvec3d_u(vec3d const&), mvec3d_u(pscalar3d),
                          mvec3d_u(vec3d const&, pscalar3d)>(),
        "copy", [](const mvec3d_u& obj) { return mvec3d_u(obj); },
#// component access
        "c0", &mvec3d_u::c0, "c1", &mvec3d_u::c1, "c2", &mvec3d_u::c2, "c3",
        &mvec3d_u::c3, sol::meta_function::to_string,
        [](const mvec3d_u& mv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("MVec3d_U({},{},{},{})", suppress_negative_zero(mv.c0),
                               suppress_negative_zero(mv.c1),
                               suppress_negative_zero(mv.c2),
                               suppress_negative_zero(mv.c3));
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec3d_u(mvec3d_u const&)>(operator-), sol::meta_function::addition,
        sol::overload(sol::resolve<mvec3d_u(mvec3d_u const&, mvec3d_u const&)>(operator+),
                      sol::resolve<mvec3d_u(mvec3d_u const&, pscalar3d)>(operator+),
                      sol::resolve<mvec3d_u(mvec3d_u const&, vec3d const&)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<mvec3d_u(mvec3d_u const&, mvec3d_u const&)>(operator-),
                      sol::resolve<mvec3d_u(mvec3d_u const&, pscalar3d)>(operator-),
                      sol::resolve<mvec3d_u(mvec3d_u const&, vec3d const&)>(operator-)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec3d_u(mvec3d_u const&, value_t)>(operator*),
                      sol::resolve<mvec3d_u(value_t, mvec3d_u const&)>(operator*),
                      sol::resolve<mvec3d_e(mvec3d_u const&, pscalar3d)>(operator*),
                      sol::resolve<mvec3d_e(mvec3d_u const&, mvec3d_u const&)>(operator*),
                      sol::resolve<mvec3d_u(mvec3d_u const&, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d_u(mvec3d_u const&, bivec3d const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3d_u(mvec3d_u const&, value_t)>(operator/),
        sol::meta_function::bitwise_right_shift,
        sol::resolve<mvec3d_u(mvec3d_u const&, scalar3d)>(operator>>));


    lua.new_usertype<mvec3d>(
        "mvec3d",
        sol::constructors<mvec3d(),
                          mvec3d(value_t, value_t, value_t, value_t, value_t, value_t,
                                 value_t, value_t),
                          mvec3d(mvec3d const&), mvec3d(mvec3d&&), mvec3d(scalar3d),
                          mvec3d(vec3d const&), mvec3d(bivec3d const&), mvec3d(pscalar3d),
                          mvec3d(scalar3d, bivec3d const&), mvec3d(mvec3d_e const&),
                          mvec3d(vec3d const&, pscalar3d), mvec3d(mvec3d_u const&),
                          mvec3d(scalar3d, vec3d const&, bivec3d const&, pscalar3d)>(),
        "copy", [](const mvec3d& obj) { return mvec3d(obj); },
        // component access
        "c0", &mvec3d::c0, "c1", &mvec3d::c1, "c2", &mvec3d::c2, "c3", &mvec3d::c3, "c4",
        &mvec3d::c4, "c5", &mvec3d::c5, "c6", &mvec3d::c6, "c7", &mvec3d::c7,
        sol::meta_function::to_string,
        [](const mvec3d& mv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "MVec3d({},{},{},{},{},{},{},{})", suppress_negative_zero(mv.c0),
                suppress_negative_zero(mv.c1), suppress_negative_zero(mv.c2),
                suppress_negative_zero(mv.c3), suppress_negative_zero(mv.c4),
                suppress_negative_zero(mv.c5), suppress_negative_zero(mv.c6),
                suppress_negative_zero(mv.c7));
        },
        sol::meta_function::unary_minus, sol::resolve<mvec3d(mvec3d const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator-),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<mvec3d(mvec3d const&, vec3d const&)>(wdg),
                      sol::resolve<mvec3d(mvec3d const&, bivec3d const&)>(wdg),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec3d(mvec3d const&, value_t)>(operator*),
                      sol::resolve<mvec3d(value_t, mvec3d const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d_u const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, vec3d const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, bivec3d const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, pscalar3d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3d(mvec3d const&, value_t)>(operator/),
        sol::meta_function::bitwise_left_shift,
        sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator<<),
        sol::meta_function::bitwise_right_shift,
        sol::overload(sol::resolve<mvec3d(mvec3d const&, scalar3d)>(operator>>),
                      sol::resolve<mvec3d(mvec3d const&, vec3d const&)>(operator>>),
                      sol::resolve<mvec3d(mvec3d const&, bivec3d const&)>(operator>>),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator>>)));
}


////////////////////////////////////////////////////////////////////////////////
// 2dp PGA types
////////////////////////////////////////////////////////////////////////////////
void register_2dp_types(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::pga;

    // Basic scalar type
    lua.new_usertype<scalar2dp>(
        "scalar2dp",
        sol::constructors<scalar2dp(), scalar2dp(value_t const&), scalar2dp(value_t&&)>(),
        "copy", [](const scalar2dp& obj) { return scalar2dp(obj); },
        sol::meta_function::to_string,
        [](const scalar2dp& s) { return fmt::format("Scalar2dp({})", double(s)); },
        sol::meta_function::unary_minus, sol::resolve<scalar2dp(scalar2dp)>(operator-),
        sol::meta_function::addition,
        sol::resolve<scalar2dp(scalar2dp, scalar2dp)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<scalar2dp(scalar2dp, scalar2dp)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<scalar2dp(scalar2dp, scalar2dp)>(operator*),
                      sol::resolve<scalar2dp(scalar2dp, value_t)>(operator*),
                      sol::resolve<scalar2dp(value_t, scalar2dp)>(operator*)),
        sol::meta_function::division,
        sol::resolve<scalar2dp(scalar2dp, value_t)>(operator/));

    // Vector type with x,y,z components
    lua.new_usertype<vec2dp>(
        "vec2dp",
        sol::constructors<vec2dp(), vec2dp(value_t, value_t, value_t),
                          vec2dp(vec2dp const&), vec2dp(vec2dp&&)>(),
        "copy", [](const vec2dp& obj) { return vec2dp(obj); },
        // component access
        "x", &vec2dp::x, "y", &vec2dp::y, "z", &vec2dp::z, sol::meta_function::to_string,
        [](const vec2dp& v) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Vec2dp({},{},{})", suppress_negative_zero(v.x),
                               suppress_negative_zero(v.y), suppress_negative_zero(v.z));
        },
        sol::meta_function::unary_minus, sol::resolve<vec2dp(vec2dp const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<vec2dp(vec2dp const&, vec2dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<vec2dp(vec2dp const&, vec2dp const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec2dp(vec2dp const&, value_t)>(operator*),
                      sol::resolve<vec2dp(value_t, vec2dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec2dp(vec2dp const&, value_t)>(operator/),
        sol::meta_function::power_of,
        sol::resolve<bivec2dp(vec2dp const&, vec2dp const&)>(wdg));

    // Bivector type with x,y,z components
    lua.new_usertype<bivec2dp>(
        "bivec2dp",
        sol::constructors<bivec2dp(), bivec2dp(value_t, value_t, value_t),
                          bivec2dp(bivec2dp const&), bivec2dp(bivec2dp&&)>(),
        "copy", [](const bivec2dp& obj) { return bivec2dp(obj); },
        // component access
        "x", &bivec2dp::x, "y", &bivec2dp::y, "z", &bivec2dp::z,
        sol::meta_function::to_string,
        [](const bivec2dp& bv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("BiVec2dp({},{},{})", suppress_negative_zero(bv.x),
                               suppress_negative_zero(bv.y),
                               suppress_negative_zero(bv.z));
        },
        sol::meta_function::unary_minus,
        sol::resolve<bivec2dp(bivec2dp const&)>(operator-), sol::meta_function::addition,
        sol::resolve<bivec2dp(bivec2dp const&, bivec2dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<bivec2dp(bivec2dp const&, bivec2dp const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<bivec2dp(bivec2dp const&, value_t)>(operator*),
                      sol::resolve<bivec2dp(value_t, bivec2dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<bivec2dp(bivec2dp const&, value_t)>(operator/));

    // Pseudoscalar type
    lua.new_usertype<pscalar2dp>(
        "pscalar2dp",
        sol::constructors<pscalar2dp(), pscalar2dp(value_t const&),
                          pscalar2dp(value_t&&)>(),
        "copy", [](const pscalar2dp& obj) { return pscalar2dp(obj); },
        sol::meta_function::to_string,
        [](const pscalar2dp& ps) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("PScalar2dp({})", suppress_negative_zero(double(ps)));
        },
        sol::meta_function::unary_minus, sol::resolve<pscalar2dp(pscalar2dp)>(operator-),
        sol::meta_function::addition,
        sol::resolve<pscalar2dp(pscalar2dp, pscalar2dp)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<pscalar2dp(pscalar2dp, pscalar2dp)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<pscalar2dp(pscalar2dp, value_t)>(operator*),
                      sol::resolve<pscalar2dp(value_t, pscalar2dp)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar2dp(pscalar2dp, value_t)>(operator/));

    // Dual numbers (for PGA-specific calculations)
    lua.new_usertype<dualnum2dp>(
        "dualnum2dp",
        sol::constructors<dualnum2dp(), dualnum2dp(value_t, value_t),
                          dualnum2dp(scalar2dp), dualnum2dp(pscalar2dp),
                          dualnum2dp(scalar2dp, pscalar2dp),
                          dualnum2dp(dualnum2dp const&), dualnum2dp(dualnum2dp&&)>(),
        "copy", [](const dualnum2dp& obj) { return dualnum2dp(obj); },
        // component access
        "c0", &dualnum2dp::c0, "c1", &dualnum2dp::c1, sol::meta_function::to_string,
        [](const dualnum2dp& dn) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("DualNum2dp({},{})", suppress_negative_zero(dn.c0),
                               suppress_negative_zero(dn.c1));
        },
        sol::meta_function::unary_minus,
        sol::resolve<dualnum2dp(dualnum2dp const&)>(operator-),
        sol::meta_function::addition,
        sol::overload(
            sol::resolve<dualnum2dp(dualnum2dp const&, dualnum2dp const&)>(operator+),
            sol::resolve<dualnum2dp(scalar2dp, dualnum2dp const&)>(operator+),
            sol::resolve<dualnum2dp(dualnum2dp const&, scalar2dp)>(operator+),
            sol::resolve<dualnum2dp(pscalar2dp, dualnum2dp const&)>(operator+),
            sol::resolve<dualnum2dp(dualnum2dp const&, pscalar2dp)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(
            sol::resolve<dualnum2dp(dualnum2dp const&, dualnum2dp const&)>(operator-),
            sol::resolve<dualnum2dp(scalar2dp, dualnum2dp const&)>(operator-),
            sol::resolve<dualnum2dp(dualnum2dp const&, scalar2dp)>(operator-),
            sol::resolve<dualnum2dp(pscalar2dp, dualnum2dp const&)>(operator-),
            sol::resolve<dualnum2dp(dualnum2dp const&, pscalar2dp)>(operator-)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<dualnum2dp(dualnum2dp const&, value_t)>(operator*),
                      sol::resolve<dualnum2dp(value_t, dualnum2dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<dualnum2dp(dualnum2dp const&, value_t)>(operator/));

    // PGA 2DP multivector types
    lua.new_usertype<mvec2dp_e>(
        "mvec2dp_e",
        sol::constructors<mvec2dp_e(), mvec2dp_e(value_t, value_t, value_t, value_t),
                          mvec2dp_e(mvec2dp_e const&), mvec2dp_e(mvec2dp_e&&),
                          mvec2dp_e(scalar2dp), mvec2dp_e(bivec2dp const&),
                          mvec2dp_e(scalar2dp, bivec2dp const&)>(),
        "copy", [](const mvec2dp_e& obj) { return mvec2dp_e(obj); },
        // component access
        "c0", &mvec2dp_e::c0, "c1", &mvec2dp_e::c1, "c2", &mvec2dp_e::c2, "c3",
        &mvec2dp_e::c3, sol::meta_function::to_string,
        [](const mvec2dp_e& mve) {
            return fmt::format("MVec2dp_E({},{},{},{})", mve.c0, mve.c1, mve.c2, mve.c3);
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec2dp_e(mvec2dp_e const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec2dp_e(mvec2dp_e const&, mvec2dp_e const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2dp_e(mvec2dp_e const&, mvec2dp_e const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(
            sol::resolve<mvec2dp_e(mvec2dp_e const&, value_t)>(operator*),
            sol::resolve<mvec2dp_e(value_t, mvec2dp_e const&)>(operator*),
            sol::resolve<mvec2dp_e(mvec2dp_e const&, mvec2dp_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2dp_e(mvec2dp_e const&, value_t)>(operator/));

    lua.new_usertype<mvec2dp_u>(
        "mvec2dp_u",
        sol::constructors<mvec2dp_u(), mvec2dp_u(value_t, value_t, value_t, value_t),
                          mvec2dp_u(mvec2dp_u const&), mvec2dp_u(mvec2dp_u&&),
                          mvec2dp_u(vec2dp const&), mvec2dp_u(pscalar2dp),
                          mvec2dp_u(vec2dp const&, pscalar2dp)>(),
        "copy", [](const mvec2dp_u& obj) { return mvec2dp_u(obj); },
        // component access
        "c0", &mvec2dp_u::c0, "c1", &mvec2dp_u::c1, "c2", &mvec2dp_u::c2, "c3",
        &mvec2dp_u::c3, sol::meta_function::to_string,
        [](const mvec2dp_u& mvu) {
            return fmt::format("MVec2dp_U({},{},{},{})", mvu.c0, mvu.c1, mvu.c2, mvu.c3);
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec2dp_u(mvec2dp_u const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec2dp_u(mvec2dp_u const&, mvec2dp_u const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2dp_u(mvec2dp_u const&, mvec2dp_u const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec2dp_u(mvec2dp_u const&, value_t)>(operator*),
                      sol::resolve<mvec2dp_u(value_t, mvec2dp_u const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2dp_u(mvec2dp_u const&, value_t)>(operator/));

    lua.new_usertype<mvec2dp>(
        "mvec2dp",
        sol::constructors<mvec2dp(),
                          mvec2dp(value_t, value_t, value_t, value_t, value_t, value_t,
                                  value_t, value_t),
                          mvec2dp(mvec2dp const&), mvec2dp(mvec2dp&&), mvec2dp(scalar2dp),
                          mvec2dp(vec2dp const&), mvec2dp(bivec2dp const&),
                          mvec2dp(pscalar2dp), mvec2dp(mvec2dp_e const&),
                          mvec2dp(mvec2dp_u const&)>(),
        "copy", [](const mvec2dp& obj) { return mvec2dp(obj); },
        // component access
        "c0", &mvec2dp::c0, "c1", &mvec2dp::c1, "c2", &mvec2dp::c2, "c3", &mvec2dp::c3,
        "c4", &mvec2dp::c4, "c5", &mvec2dp::c5, "c6", &mvec2dp::c6, "c7", &mvec2dp::c7,
        sol::meta_function::to_string,
        [](const mvec2dp& mv) {
            return fmt::format("MVec2dp({},{},{},{},{},{},{},{})", mv.c0, mv.c1, mv.c2,
                               mv.c3, mv.c4, mv.c5, mv.c6, mv.c7);
        },
        sol::meta_function::unary_minus, sol::resolve<mvec2dp(mvec2dp const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec2dp(mvec2dp const&, value_t)>(operator*),
                      sol::resolve<mvec2dp(value_t, mvec2dp const&)>(operator*),
                      sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2dp(mvec2dp const&, value_t)>(operator/));
}


////////////////////////////////////////////////////////////////////////////////
// 3dp PGA types
////////////////////////////////////////////////////////////////////////////////
void register_3dp_types(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::pga;

    // Basic scalar type
    lua.new_usertype<scalar3dp>(
        "scalar3dp",
        sol::constructors<scalar3dp(), scalar3dp(value_t const&), scalar3dp(value_t&&)>(),
        "copy", [](const scalar3dp& obj) { return scalar3dp(obj); },
        sol::meta_function::to_string,
        [](const scalar3dp& s) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Scalar3dp({})", suppress_negative_zero(double(s)));
        },
        sol::meta_function::unary_minus, sol::resolve<scalar3dp(scalar3dp)>(operator-),
        sol::meta_function::addition,
        sol::resolve<scalar3dp(scalar3dp, scalar3dp)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<scalar3dp(scalar3dp, scalar3dp)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<scalar3dp(scalar3dp, scalar3dp)>(operator*),
                      sol::resolve<scalar3dp(scalar3dp, value_t)>(operator*),
                      sol::resolve<scalar3dp(value_t, scalar3dp)>(operator*)),
        sol::meta_function::division,
        sol::resolve<scalar3dp(scalar3dp, value_t)>(operator/));

    // Vector type with x,y,z,w components
    lua.new_usertype<vec3dp>(
        "vec3dp",
        sol::constructors<vec3dp(), vec3dp(value_t, value_t, value_t, value_t),
                          vec3dp(vec3dp const&), vec3dp(vec3dp&&)>(),
        "copy", [](const vec3dp& obj) { return vec3dp(obj); },
        // component access
        "x", &vec3dp::x, "y", &vec3dp::y, "z", &vec3dp::z, "w", &vec3dp::w,
        sol::meta_function::to_string,
        [](const vec3dp& v) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Vec3dp({},{},{},{})", suppress_negative_zero(v.x),
                               suppress_negative_zero(v.y), suppress_negative_zero(v.z),
                               suppress_negative_zero(v.w));
        },
        sol::meta_function::unary_minus, sol::resolve<vec3dp(vec3dp const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<vec3dp(vec3dp const&, vec3dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<vec3dp(vec3dp const&, vec3dp const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec3dp(vec3dp const&, value_t)>(operator*),
                      sol::resolve<vec3dp(value_t, vec3dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec3dp(vec3dp const&, value_t)>(operator/),
        sol::meta_function::power_of,
        sol::resolve<bivec3dp(vec3dp const&, vec3dp const&)>(wdg));

    // Bivector type with vx,vy,vz,mx,my,mz components
    lua.new_usertype<bivec3dp>(
        "bivec3dp",
        sol::constructors<bivec3dp(),
                          bivec3dp(value_t, value_t, value_t, value_t, value_t, value_t),
                          bivec3dp(bivec3dp const&), bivec3dp(bivec3dp&&)>(),
        "copy", [](const bivec3dp& obj) { return bivec3dp(obj); },
        // component access
        "vx", &bivec3dp::vx, "vy", &bivec3dp::vy, "vz", &bivec3dp::vz, "mx",
        &bivec3dp::mx, "my", &bivec3dp::my, "mz", &bivec3dp::mz,
        sol::meta_function::to_string,
        [](const bivec3dp& bv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "BiVec3dp({},{},{},{},{},{})", suppress_negative_zero(bv.vx),
                suppress_negative_zero(bv.vy), suppress_negative_zero(bv.vz),
                suppress_negative_zero(bv.mx), suppress_negative_zero(bv.my),
                suppress_negative_zero(bv.mz));
        },
        sol::meta_function::unary_minus,
        sol::resolve<bivec3dp(bivec3dp const&)>(operator-), sol::meta_function::addition,
        sol::resolve<bivec3dp(bivec3dp const&, bivec3dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<bivec3dp(bivec3dp const&, bivec3dp const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<bivec3dp(bivec3dp const&, value_t)>(operator*),
                      sol::resolve<bivec3dp(value_t, bivec3dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<bivec3dp(bivec3dp const&, value_t)>(operator/));

    // Trivector type with x,y,z,w components
    lua.new_usertype<trivec3dp>(
        "trivec3dp",
        sol::constructors<trivec3dp(), trivec3dp(value_t, value_t, value_t, value_t),
                          trivec3dp(trivec3dp const&), trivec3dp(trivec3dp&&)>(),
        "copy", [](const trivec3dp& obj) { return trivec3dp(obj); },
        // component access
        "x", &trivec3dp::x, "y", &trivec3dp::y, "z", &trivec3dp::z, "w", &trivec3dp::w,
        sol::meta_function::to_string,
        [](const trivec3dp& tv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("TriVec3dp({},{},{},{})", suppress_negative_zero(tv.x),
                               suppress_negative_zero(tv.y), suppress_negative_zero(tv.z),
                               suppress_negative_zero(tv.w));
        },
        sol::meta_function::unary_minus,
        sol::resolve<trivec3dp(trivec3dp const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<trivec3dp(trivec3dp const&, trivec3dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<trivec3dp(trivec3dp const&, trivec3dp const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<trivec3dp(trivec3dp const&, value_t)>(operator*),
                      sol::resolve<trivec3dp(value_t, trivec3dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<trivec3dp(trivec3dp const&, value_t)>(operator/));

    // Pseudoscalar type
    lua.new_usertype<pscalar3dp>(
        "pscalar3dp",
        sol::constructors<pscalar3dp(), pscalar3dp(value_t const&),
                          pscalar3dp(value_t&&)>(),
        "copy", [](const pscalar3dp& obj) { return pscalar3dp(obj); },
        sol::meta_function::to_string,
        [](const pscalar3dp& ps) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("PScalar3dp({})", suppress_negative_zero(double(ps)));
        },
        sol::meta_function::unary_minus, sol::resolve<pscalar3dp(pscalar3dp)>(operator-),
        sol::meta_function::addition,
        sol::resolve<pscalar3dp(pscalar3dp, pscalar3dp)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<pscalar3dp(pscalar3dp, pscalar3dp)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<pscalar3dp(pscalar3dp, value_t)>(operator*),
                      sol::resolve<pscalar3dp(value_t, pscalar3dp)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar3dp(pscalar3dp, value_t)>(operator/));

    // Dual numbers (for PGA-specific calculations)
    lua.new_usertype<dualnum3dp>(
        "dualnum3dp",
        sol::constructors<dualnum3dp(), dualnum3dp(value_t, value_t),
                          dualnum3dp(scalar3dp), dualnum3dp(pscalar3dp),
                          dualnum3dp(scalar3dp, pscalar3dp),
                          dualnum3dp(dualnum3dp const&), dualnum3dp(dualnum3dp&&)>(),
        "copy", [](const dualnum3dp& obj) { return dualnum3dp(obj); },
        // component access
        "c0", &dualnum3dp::c0, "c1", &dualnum3dp::c1, sol::meta_function::to_string,
        [](const dualnum3dp& dn) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("DualNum3dp({},{})", suppress_negative_zero(dn.c0),
                               suppress_negative_zero(dn.c1));
        },
        sol::meta_function::unary_minus,
        sol::resolve<dualnum3dp(dualnum3dp const&)>(operator-),
        sol::meta_function::addition,
        sol::overload(
            sol::resolve<dualnum3dp(dualnum3dp const&, dualnum3dp const&)>(operator+),
            sol::resolve<dualnum3dp(scalar3dp, dualnum3dp const&)>(operator+),
            sol::resolve<dualnum3dp(dualnum3dp const&, scalar3dp)>(operator+),
            sol::resolve<dualnum3dp(pscalar3dp, dualnum3dp const&)>(operator+),
            sol::resolve<dualnum3dp(dualnum3dp const&, pscalar3dp)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(
            sol::resolve<dualnum3dp(dualnum3dp const&, dualnum3dp const&)>(operator-),
            sol::resolve<dualnum3dp(scalar3dp, dualnum3dp const&)>(operator-),
            sol::resolve<dualnum3dp(dualnum3dp const&, scalar3dp)>(operator-),
            sol::resolve<dualnum3dp(pscalar3dp, dualnum3dp const&)>(operator-),
            sol::resolve<dualnum3dp(dualnum3dp const&, pscalar3dp)>(operator-)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<dualnum3dp(dualnum3dp const&, value_t)>(operator*),
                      sol::resolve<dualnum3dp(value_t, dualnum3dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<dualnum3dp(dualnum3dp const&, value_t)>(operator/));

    // PGA 3DP multivector types
    lua.new_usertype<mvec3dp_e>(
        "mvec3dp_e",
        sol::constructors<mvec3dp_e(),
                          mvec3dp_e(value_t, value_t, value_t, value_t, value_t, value_t,
                                    value_t, value_t),
                          mvec3dp_e(mvec3dp_e const&), mvec3dp_e(mvec3dp_e&&),
                          mvec3dp_e(scalar3dp), mvec3dp_e(bivec3dp const&),
                          mvec3dp_e(pscalar3dp), mvec3dp_e(scalar3dp, bivec3dp const&),
                          mvec3dp_e(bivec3dp const&, pscalar3dp),
                          mvec3dp_e(scalar3dp, pscalar3dp),
                          mvec3dp_e(scalar3dp, bivec3dp const&, pscalar3dp)>(),
        "copy", [](const mvec3dp_e& obj) { return mvec3dp_e(obj); },
        // component access
        "c0", &mvec3dp_e::c0, "c1", &mvec3dp_e::c1, "c2", &mvec3dp_e::c2, "c3",
        &mvec3dp_e::c3, "c4", &mvec3dp_e::c4, "c5", &mvec3dp_e::c5, "c6", &mvec3dp_e::c6,
        "c7", &mvec3dp_e::c7, sol::meta_function::to_string,
        [](const mvec3dp_e& mve) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "MVec3dp_E({},{},{},{},{},{},{},{})", suppress_negative_zero(mve.c0),
                suppress_negative_zero(mve.c1), suppress_negative_zero(mve.c2),
                suppress_negative_zero(mve.c3), suppress_negative_zero(mve.c4),
                suppress_negative_zero(mve.c5), suppress_negative_zero(mve.c6),
                suppress_negative_zero(mve.c7));
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec3dp_e(mvec3dp_e const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec3dp_e(mvec3dp_e const&, mvec3dp_e const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec3dp_e(mvec3dp_e const&, mvec3dp_e const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(
            sol::resolve<mvec3dp_e(mvec3dp_e const&, value_t)>(operator*),
            sol::resolve<mvec3dp_e(value_t, mvec3dp_e const&)>(operator*),
            sol::resolve<mvec3dp_e(mvec3dp_e const&, mvec3dp_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3dp_e(mvec3dp_e const&, value_t)>(operator/));

    lua.new_usertype<mvec3dp_u>(
        "mvec3dp_u",
        sol::constructors<mvec3dp_u(),
                          mvec3dp_u(value_t, value_t, value_t, value_t, value_t, value_t,
                                    value_t, value_t),
                          mvec3dp_u(mvec3dp_u const&), mvec3dp_u(mvec3dp_u&&),
                          mvec3dp_u(vec3dp const&), mvec3dp_u(trivec3dp const&),
                          mvec3dp_u(vec3dp const&, trivec3dp const&)>(),
        "copy", [](const mvec3dp_u& obj) { return mvec3dp_u(obj); },
        // component access
        "c0", &mvec3dp_u::c0, "c1", &mvec3dp_u::c1, "c2", &mvec3dp_u::c2, "c3",
        &mvec3dp_u::c3, "c4", &mvec3dp_u::c4, "c5", &mvec3dp_u::c5, "c6", &mvec3dp_u::c6,
        "c7", &mvec3dp_u::c7, sol::meta_function::to_string,
        [](const mvec3dp_u& mvu) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "MVec3dp_U({},{},{},{},{},{},{},{})", suppress_negative_zero(mvu.c0),
                suppress_negative_zero(mvu.c1), suppress_negative_zero(mvu.c2),
                suppress_negative_zero(mvu.c3), suppress_negative_zero(mvu.c4),
                suppress_negative_zero(mvu.c5), suppress_negative_zero(mvu.c6),
                suppress_negative_zero(mvu.c7));
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec3dp_u(mvec3dp_u const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec3dp_u(mvec3dp_u const&, mvec3dp_u const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec3dp_u(mvec3dp_u const&, mvec3dp_u const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec3dp_u(mvec3dp_u const&, value_t)>(operator*),
                      sol::resolve<mvec3dp_u(value_t, mvec3dp_u const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3dp_u(mvec3dp_u const&, value_t)>(operator/));

    lua.new_usertype<mvec3dp>(
        "mvec3dp",
        sol::constructors<mvec3dp(),
                          mvec3dp(value_t, value_t, value_t, value_t, value_t, value_t,
                                  value_t, value_t, value_t, value_t, value_t, value_t,
                                  value_t, value_t, value_t, value_t),
                          mvec3dp(mvec3dp const&), mvec3dp(mvec3dp&&), mvec3dp(scalar3dp),
                          mvec3dp(vec3dp const&), mvec3dp(bivec3dp const&),
                          mvec3dp(trivec3dp const&), mvec3dp(pscalar3dp),
                          mvec3dp(mvec3dp_e const&), mvec3dp(mvec3dp_u const&)>(),
        "copy", [](const mvec3dp& obj) { return mvec3dp(obj); },
        // component access
        "c0", &mvec3dp::c0, "c1", &mvec3dp::c1, "c2", &mvec3dp::c2, "c3", &mvec3dp::c3,
        "c4", &mvec3dp::c4, "c5", &mvec3dp::c5, "c6", &mvec3dp::c6, "c7", &mvec3dp::c7,
        "c8", &mvec3dp::c8, "c9", &mvec3dp::c9, "c10", &mvec3dp::c10, "c11",
        &mvec3dp::c11, "c12", &mvec3dp::c12, "c13", &mvec3dp::c13, "c14", &mvec3dp::c14,
        "c15", &mvec3dp::c15, sol::meta_function::to_string,
        [](const mvec3dp& mv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "MVec3dp({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})",
                suppress_negative_zero(mv.c0), suppress_negative_zero(mv.c1),
                suppress_negative_zero(mv.c2), suppress_negative_zero(mv.c3),
                suppress_negative_zero(mv.c4), suppress_negative_zero(mv.c5),
                suppress_negative_zero(mv.c6), suppress_negative_zero(mv.c7),
                suppress_negative_zero(mv.c8), suppress_negative_zero(mv.c9),
                suppress_negative_zero(mv.c10), suppress_negative_zero(mv.c11),
                suppress_negative_zero(mv.c12), suppress_negative_zero(mv.c13),
                suppress_negative_zero(mv.c14), suppress_negative_zero(mv.c15));
        },
        sol::meta_function::unary_minus, sol::resolve<mvec3dp(mvec3dp const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec3dp(mvec3dp const&, value_t)>(operator*),
                      sol::resolve<mvec3dp(value_t, mvec3dp const&)>(operator*),
                      sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3dp(mvec3dp const&, value_t)>(operator/));
}


////////////////////////////////////////////////////////////////////////////////
// functions commonly used across 2d and 3d types
////////////////////////////////////////////////////////////////////////////////
void register_functions(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::ega;
    using namespace hd::ga::pga;
    using namespace hd::ga::sta;

    lua.set_function("nrm_sq",
                     sol::overload( // EGA 2D types
                         sol::resolve<value_t(scalar2d)>(nrm_sq),
                         sol::resolve<value_t(vec2d const&)>(nrm_sq),
                         sol::resolve<value_t(pscalar2d)>(nrm_sq),
                         sol::resolve<value_t(mvec2d_e const&)>(nrm_sq),
                         sol::resolve<value_t(mvec2d const&)>(nrm_sq),
                         // EGA 3D types
                         sol::resolve<value_t(scalar3d)>(nrm_sq),
                         sol::resolve<value_t(vec3d const&)>(nrm_sq),
                         sol::resolve<value_t(bivec3d const&)>(nrm_sq),
                         sol::resolve<value_t(pscalar3d)>(nrm_sq),
                         sol::resolve<value_t(mvec3d_e const&)>(nrm_sq),
                         sol::resolve<value_t(mvec3d_u const&)>(nrm_sq),
                         sol::resolve<value_t(mvec3d const&)>(nrm_sq),
                         // sta
                         sol::resolve<value_t(scalar4ds)>(nrm_sq),
                         sol::resolve<value_t(vec4ds const&)>(nrm_sq),
                         sol::resolve<value_t(bivec4ds const&)>(nrm_sq),
                         sol::resolve<value_t(trivec4ds const&)>(nrm_sq),
                         sol::resolve<value_t(pscalar4ds)>(nrm_sq),
                         sol::resolve<value_t(mvec4ds_e const&)>(nrm_sq),
                         sol::resolve<value_t(mvec4ds_u const&)>(nrm_sq),
                         sol::resolve<value_t(mvec4ds const&)>(nrm_sq)));

    lua.set_function("nrm", sol::overload( // EGA 2D types
                                sol::resolve<value_t(scalar2d)>(nrm),
                                sol::resolve<value_t(vec2d const&)>(nrm),
                                sol::resolve<value_t(pscalar2d)>(nrm),
                                sol::resolve<value_t(mvec2d_e const&)>(nrm),
                                sol::resolve<value_t(mvec2d const&)>(nrm),
                                // EGA 3D types
                                sol::resolve<value_t(scalar3d)>(nrm),
                                sol::resolve<value_t(vec3d const&)>(nrm),
                                sol::resolve<value_t(bivec3d const&)>(nrm),
                                sol::resolve<value_t(pscalar3d)>(nrm),
                                sol::resolve<value_t(mvec3d_e const&)>(nrm),
                                sol::resolve<value_t(mvec3d_u const&)>(nrm),
                                sol::resolve<value_t(mvec3d const&)>(nrm),
                                // sta
                                sol::resolve<scalar4ds(scalar4ds)>(nrm),
                                sol::resolve<scalar4ds(vec4ds const&)>(nrm),
                                sol::resolve<scalar4ds(bivec4ds const&)>(nrm),
                                sol::resolve<scalar4ds(trivec4ds const&)>(nrm),
                                sol::resolve<scalar4ds(pscalar4ds)>(nrm),
                                sol::resolve<scalar4ds(mvec4ds_e const&)>(nrm),
                                sol::resolve<scalar4ds(mvec4ds_u const&)>(nrm),
                                sol::resolve<scalar4ds(mvec4ds const&)>(nrm)));

    lua.set_function("bulk_nrm_sq",
                     sol::overload(
                         // PGA 2DP types: use bulk_nrm_sq()
                         sol::resolve<value_t(scalar2dp)>(bulk_nrm_sq),
                         sol::resolve<value_t(vec2dp const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(bivec2dp const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(pscalar2dp)>(bulk_nrm_sq),
                         sol::resolve<value_t(mvec2dp_e const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(mvec2dp_u const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(mvec2dp const&)>(bulk_nrm_sq),
                         // sol::resolve<value_t(dualnum2dp const&)>(bulk_nrm_sq),
                         // PGA 3DP types: use bulk_nrm_sq()
                         sol::resolve<value_t(scalar3dp)>(bulk_nrm_sq),
                         sol::resolve<value_t(vec3dp const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(bivec3dp const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(trivec3dp const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(pscalar3dp)>(bulk_nrm_sq),
                         sol::resolve<value_t(mvec3dp_e const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(mvec3dp_u const&)>(bulk_nrm_sq),
                         sol::resolve<value_t(mvec3dp const&)>(bulk_nrm_sq)
                         // , sol::resolve<value_t(dualnum3dp const&)>(bulk_nrm_sq)
                         ));

    lua.set_function("weight_nrm_sq",
                     sol::overload(
                         // PGA 2DP types: use weight_nrm_sq()
                         sol::resolve<value_t(scalar2dp)>(weight_nrm_sq),
                         sol::resolve<value_t(vec2dp const&)>(weight_nrm_sq),
                         sol::resolve<value_t(bivec2dp const&)>(weight_nrm_sq),
                         sol::resolve<value_t(pscalar2dp)>(weight_nrm_sq),
                         sol::resolve<value_t(mvec2dp_e const&)>(weight_nrm_sq),
                         sol::resolve<value_t(mvec2dp_u const&)>(weight_nrm_sq),
                         sol::resolve<value_t(mvec2dp const&)>(weight_nrm_sq),
                         // sol::resolve<value_t(dualnum2dp const&)>(weight_nrm_sq),
                         // PGA 3DP types: use weight_nrm_sq()
                         sol::resolve<value_t(scalar3dp)>(weight_nrm_sq),
                         sol::resolve<value_t(vec3dp const&)>(weight_nrm_sq),
                         sol::resolve<value_t(bivec3dp const&)>(weight_nrm_sq),
                         sol::resolve<value_t(trivec3dp const&)>(weight_nrm_sq),
                         sol::resolve<value_t(pscalar3dp)>(weight_nrm_sq),
                         sol::resolve<value_t(mvec3dp_e const&)>(weight_nrm_sq),
                         sol::resolve<value_t(mvec3dp_u const&)>(weight_nrm_sq),
                         sol::resolve<value_t(mvec3dp const&)>(weight_nrm_sq)
                         // , sol::resolve<value_t(dualnum3dp const&)>(weight_nrm_sq)
                         ));

    lua.set_function("geom_nrm_sq",
                     sol::overload(
                         // PGA 2DP types: use geom_nrm_sq()
                         sol::resolve<value_t(scalar2dp)>(geom_nrm_sq),
                         sol::resolve<value_t(vec2dp const&)>(geom_nrm_sq),
                         sol::resolve<value_t(bivec2dp const&)>(geom_nrm_sq),
                         sol::resolve<value_t(pscalar2dp)>(geom_nrm_sq),
                         sol::resolve<value_t(mvec2dp_e const&)>(geom_nrm_sq),
                         sol::resolve<value_t(mvec2dp_u const&)>(geom_nrm_sq),
                         sol::resolve<value_t(mvec2dp const&)>(geom_nrm_sq),
                         // sol::resolve<value_t(dualnum2dp const&)>(geom_nrm_sq),
                         // PGA 3DP types: use geom_nrm_sq()
                         sol::resolve<value_t(scalar3dp)>(geom_nrm_sq),
                         sol::resolve<value_t(vec3dp const&)>(geom_nrm_sq),
                         sol::resolve<value_t(bivec3dp const&)>(geom_nrm_sq),
                         sol::resolve<value_t(trivec3dp const&)>(geom_nrm_sq),
                         sol::resolve<value_t(pscalar3dp)>(geom_nrm_sq),
                         sol::resolve<value_t(mvec3dp_e const&)>(geom_nrm_sq),
                         sol::resolve<value_t(mvec3dp_u const&)>(geom_nrm_sq),
                         sol::resolve<value_t(mvec3dp const&)>(geom_nrm_sq)
                         // , sol::resolve<value_t(dualnum3dp const&)>(geom_nrm_sq)
                         ));


    lua.set_function("bulk_nrm",
                     sol::overload(
                         // PGA 2DP types: use bulk_nrm()
                         sol::resolve<scalar2dp(scalar2dp)>(bulk_nrm),
                         sol::resolve<scalar2dp(vec2dp const&)>(bulk_nrm),
                         sol::resolve<scalar2dp(bivec2dp const&)>(bulk_nrm),
                         sol::resolve<scalar2dp(pscalar2dp)>(bulk_nrm),
                         sol::resolve<scalar2dp(mvec2dp_e const&)>(bulk_nrm),
                         sol::resolve<scalar2dp(mvec2dp_u const&)>(bulk_nrm),
                         sol::resolve<scalar2dp(mvec2dp const&)>(bulk_nrm),
                         // sol::resolve<scalar2dp(dualnum2dp const&)>(bulk_nrm),
                         // PGA 3DP types: use bulk_nrm()
                         sol::resolve<scalar3dp(scalar3dp)>(bulk_nrm),
                         sol::resolve<scalar3dp(vec3dp const&)>(bulk_nrm),
                         sol::resolve<scalar3dp(bivec3dp const&)>(bulk_nrm),
                         sol::resolve<scalar3dp(trivec3dp const&)>(bulk_nrm),
                         sol::resolve<scalar3dp(pscalar3dp)>(bulk_nrm),
                         sol::resolve<scalar3dp(mvec3dp_e const&)>(bulk_nrm),
                         sol::resolve<scalar3dp(mvec3dp_u const&)>(bulk_nrm),
                         sol::resolve<scalar3dp(mvec3dp const&)>(bulk_nrm)
                         // , sol::resolve<scalar3dp(dualnum3dp const&)>(bulk_nrm)
                         ));

    lua.set_function("weight_nrm",
                     sol::overload(
                         // PGA 2DP types: use weight_nrm()
                         sol::resolve<pscalar2dp(scalar2dp)>(weight_nrm),
                         sol::resolve<pscalar2dp(vec2dp const&)>(weight_nrm),
                         sol::resolve<pscalar2dp(bivec2dp const&)>(weight_nrm),
                         sol::resolve<pscalar2dp(pscalar2dp)>(weight_nrm),
                         sol::resolve<pscalar2dp(mvec2dp_e const&)>(weight_nrm),
                         sol::resolve<pscalar2dp(mvec2dp_u const&)>(weight_nrm),
                         sol::resolve<pscalar2dp(mvec2dp const&)>(weight_nrm),
                         //  sol::resolve<pscalar2dp(dualnum2dp const&)>(weight_nrm),
                         // PGA 3DP types: use weight_nrm()
                         sol::resolve<pscalar3dp(scalar3dp)>(weight_nrm),
                         sol::resolve<pscalar3dp(vec3dp const&)>(weight_nrm),
                         sol::resolve<pscalar3dp(bivec3dp const&)>(weight_nrm),
                         sol::resolve<pscalar3dp(trivec3dp const&)>(weight_nrm),
                         sol::resolve<pscalar3dp(pscalar3dp)>(weight_nrm),
                         sol::resolve<pscalar3dp(mvec3dp_e const&)>(weight_nrm),
                         sol::resolve<pscalar3dp(mvec3dp_u const&)>(weight_nrm),
                         sol::resolve<pscalar3dp(mvec3dp const&)>(weight_nrm)
                         // , sol::resolve<pscalar3dp(dualnum3dp const&)>(weight_nrm)
                         ));

    lua.set_function("geom_nrm",
                     sol::overload(
                         // PGA 2DP types: use geom_nrm()
                         sol::resolve<dualnum2dp(scalar2dp)>(geom_nrm),
                         sol::resolve<dualnum2dp(vec2dp const&)>(geom_nrm),
                         sol::resolve<dualnum2dp(bivec2dp const&)>(geom_nrm),
                         sol::resolve<dualnum2dp(pscalar2dp)>(geom_nrm),
                         sol::resolve<dualnum2dp(mvec2dp_e const&)>(geom_nrm),
                         sol::resolve<dualnum2dp(mvec2dp_u const&)>(geom_nrm),
                         sol::resolve<dualnum2dp(mvec2dp const&)>(geom_nrm),
                         //  sol::resolve<dualnum2dp(dualnum2dp const&)>(geom_nrm),
                         // PGA 3DP types: use geom_nrm()
                         sol::resolve<dualnum3dp(scalar3dp)>(geom_nrm),
                         sol::resolve<dualnum3dp(vec3dp const&)>(geom_nrm),
                         sol::resolve<dualnum3dp(bivec3dp const&)>(geom_nrm),
                         sol::resolve<dualnum3dp(trivec3dp const&)>(geom_nrm),
                         sol::resolve<dualnum3dp(pscalar3dp)>(geom_nrm),
                         sol::resolve<dualnum3dp(mvec3dp_e const&)>(geom_nrm),
                         sol::resolve<dualnum3dp(mvec3dp_u const&)>(geom_nrm),
                         sol::resolve<dualnum3dp(mvec3dp const&)>(geom_nrm)
                         // , sol::resolve<dualnum3dp(dualnum3dp const&)>(geom_nrm)
                         ));

    // normalization in EGA
    lua.set_function("normalize",
                     sol::overload(
                         // EGA 2D types
                         sol::resolve<scalar2d(scalar2d)>(normalize),
                         sol::resolve<vec2d(vec2d const&)>(normalize),
                         sol::resolve<pscalar2d(pscalar2d)>(normalize),
                         // EGA 3D types
                         sol::resolve<scalar3d(scalar3d)>(normalize),
                         sol::resolve<vec3d(vec3d const&)>(normalize),
                         sol::resolve<bivec3d(bivec3d const&)>(normalize),
                         sol::resolve<pscalar3d(pscalar3d)>(normalize),
                         // sta
                         sol::resolve<scalar4ds(scalar4ds)>(normalize),
                         sol::resolve<vec4ds(vec4ds const&)>(normalize),
                         sol::resolve<bivec4ds(bivec4ds const&)>(normalize),
                         sol::resolve<trivec4ds(trivec4ds const&)>(normalize)));

    // bulk_normalization in PGA
    lua.set_function("bulk_normalize",
                     sol::overload(
                         // PGA 2DP types
                         sol::resolve<scalar2dp(scalar2dp)>(bulk_normalize),
                         sol::resolve<vec2dp(vec2dp const&)>(bulk_normalize),
                         sol::resolve<bivec2dp(bivec2dp const&)>(bulk_normalize),
                         //  sol::resolve<pscalar2dp(pscalar2dp)>(bulk_normalize),
                         sol::resolve<mvec2dp_e(mvec2dp_e const&)>(bulk_normalize),
                         sol::resolve<mvec2dp_u(mvec2dp_u const&)>(bulk_normalize),
                         sol::resolve<mvec2dp(mvec2dp const&)>(bulk_normalize),
                         // PGA 3DP types
                         sol::resolve<scalar3dp(scalar3dp)>(bulk_normalize),
                         sol::resolve<vec3dp(vec3dp const&)>(bulk_normalize),
                         sol::resolve<bivec3dp(bivec3dp const&)>(bulk_normalize),
                         sol::resolve<trivec3dp(trivec3dp const&)>(bulk_normalize),
                         //  sol::resolve<pscalar3dp(pscalar3dp)>(bulk_normalize),
                         sol::resolve<mvec3dp_e(mvec3dp_e const&)>(bulk_normalize),
                         sol::resolve<mvec3dp_u(mvec3dp_u const&)>(bulk_normalize),
                         sol::resolve<mvec3dp(mvec3dp const&)>(bulk_normalize)));

    // PGA unitize functions (analogue to EGA normalize, but for weight == 1.0)
    // = weight_normalize()
    lua.set_function("unitize",
                     sol::overload(
                         // PGA 2DP types
                         sol::resolve<vec2dp(vec2dp const&)>(unitize),
                         sol::resolve<bivec2dp(bivec2dp const&)>(unitize),
                         sol::resolve<pscalar2dp(pscalar2dp)>(unitize),
                         sol::resolve<mvec2dp_e(mvec2dp_e const&)>(unitize),
                         sol::resolve<mvec2dp_u(mvec2dp_u const&)>(unitize),
                         sol::resolve<mvec2dp(mvec2dp const&)>(unitize),
                         sol::resolve<point2dp(point2dp const&)>(unitize),
                         sol::resolve<line2d(line2d const&)>(unitize),
                         sol::resolve<dualnum2dp(dualnum2dp const&)>(unitize),
                         // PGA 3DP types
                         sol::resolve<vec3dp(vec3dp const&)>(unitize),
                         sol::resolve<bivec3dp(bivec3dp const&)>(unitize),
                         sol::resolve<trivec3dp(trivec3dp const&)>(unitize),
                         sol::resolve<pscalar3dp(pscalar3dp)>(unitize),
                         sol::resolve<mvec3dp_e(mvec3dp_e const&)>(unitize),
                         sol::resolve<mvec3dp_u(mvec3dp_u const&)>(unitize),
                         sol::resolve<mvec3dp(mvec3dp const&)>(unitize),
                         sol::resolve<point3dp(point3dp const&)>(unitize),
                         sol::resolve<line3d(line3d const&)>(unitize),
                         sol::resolve<plane3d(plane3d const&)>(unitize),
                         sol::resolve<dualnum3dp(dualnum3dp const&)>(unitize)));

    ////////////////////////////////////////////////////////////////////////////////
    // PGA-specific bulk and weight operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("bulk", sol::overload(
                                 // PGA 2DP bulk operations
                                 sol::resolve<scalar2dp(scalar2dp)>(bulk),
                                 sol::resolve<vec2dp(vec2dp const&)>(bulk),
                                 sol::resolve<bivec2dp(bivec2dp const&)>(bulk),
                                 sol::resolve<pscalar2dp(pscalar2dp)>(bulk),
                                 sol::resolve<mvec2dp_e(mvec2dp_e const&)>(bulk),
                                 sol::resolve<mvec2dp_u(mvec2dp_u const&)>(bulk),
                                 sol::resolve<mvec2dp(mvec2dp const&)>(bulk),
                                 // PGA 3DP bulk operations
                                 sol::resolve<scalar3dp(scalar3dp)>(bulk),
                                 sol::resolve<vec3dp(vec3dp const&)>(bulk),
                                 sol::resolve<bivec3dp(bivec3dp const&)>(bulk),
                                 sol::resolve<trivec3dp(trivec3dp const&)>(bulk),
                                 sol::resolve<pscalar3dp(pscalar3dp)>(bulk),
                                 sol::resolve<mvec3dp_e(mvec3dp_e const&)>(bulk),
                                 sol::resolve<mvec3dp_u(mvec3dp_u const&)>(bulk),
                                 sol::resolve<mvec3dp(mvec3dp const&)>(bulk)));

    lua.set_function("weight", sol::overload(
                                   // PGA 2DP weight operations
                                   sol::resolve<scalar2dp(scalar2dp)>(weight),
                                   sol::resolve<vec2dp(vec2dp const&)>(weight),
                                   sol::resolve<bivec2dp(bivec2dp const&)>(weight),
                                   sol::resolve<pscalar2dp(pscalar2dp)>(weight),
                                   sol::resolve<mvec2dp_e(mvec2dp_e const&)>(weight),
                                   sol::resolve<mvec2dp_u(mvec2dp_u const&)>(weight),
                                   sol::resolve<mvec2dp(mvec2dp const&)>(weight),
                                   // PGA 3DP weight operations
                                   sol::resolve<scalar3dp(scalar3dp)>(weight),
                                   sol::resolve<vec3dp(vec3dp const&)>(weight),
                                   sol::resolve<bivec3dp(bivec3dp const&)>(weight),
                                   sol::resolve<trivec3dp(trivec3dp const&)>(weight),
                                   sol::resolve<pscalar3dp(pscalar3dp)>(weight),
                                   sol::resolve<mvec3dp_e(mvec3dp_e const&)>(weight),
                                   sol::resolve<mvec3dp_u(mvec3dp_u const&)>(weight),
                                   sol::resolve<mvec3dp(mvec3dp const&)>(weight)));


    ////////////////////////////////////////////////////////////////////////////////
    // PGA-specific dual operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("bulk_dual", sol::overload(
                                      // PGA 2DP bulk_dual operations
                                      sol::resolve<pscalar2dp(scalar2dp)>(bulk_dual),
                                      sol::resolve<bivec2dp(vec2dp const&)>(bulk_dual),
                                      sol::resolve<vec2dp(bivec2dp const&)>(bulk_dual),
                                      sol::resolve<scalar2dp(pscalar2dp)>(bulk_dual)));

    lua.set_function("weight_dual",
                     sol::overload(
                         // PGA 2DP weight_dual operations
                         sol::resolve<pscalar2dp(scalar2dp)>(weight_dual),
                         sol::resolve<bivec2dp(vec2dp const&)>(weight_dual),
                         sol::resolve<vec2dp(bivec2dp const&)>(weight_dual),
                         sol::resolve<scalar2dp(pscalar2dp)>(weight_dual)));


    lua.set_function("r_bulk_dual",
                     sol::overload(
                         // PGA 3DP bulk_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(r_bulk_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(r_bulk_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(r_bulk_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(r_bulk_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(r_bulk_dual)));

    lua.set_function("r_weight_dual",
                     sol::overload(
                         // PGA 3DP weight_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(r_weight_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(r_weight_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(r_weight_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(r_weight_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(r_weight_dual)));

    lua.set_function("l_bulk_dual",
                     sol::overload(
                         // PGA 3DP bulk_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(l_bulk_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(l_bulk_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(l_bulk_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(l_bulk_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(l_bulk_dual)));

    lua.set_function("l_weight_dual",
                     sol::overload(
                         // PGA 3DP weight_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(l_weight_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(l_weight_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(l_weight_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(l_weight_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(l_weight_dual)));


    ////////////////////////////////////////////////////////////////////////////////
    // grade operations on multivectors
    ////////////////////////////////////////////////////////////////////////////////
    lua.set_function("gr0", sol::overload(
                                // EGA grade 0 operations
                                sol::resolve<scalar2d(mvec2d_e const&)>(gr0),
                                sol::resolve<scalar2d(mvec2d const&)>(gr0),
                                sol::resolve<scalar3d(mvec3d_e const&)>(gr0),
                                sol::resolve<scalar3d(mvec3d const&)>(gr0),
                                // PGA grade 0 operations
                                sol::resolve<scalar2dp(mvec2dp_e const&)>(gr0),
                                sol::resolve<scalar2dp(mvec2dp const&)>(gr0),
                                sol::resolve<scalar2dp(dualnum2dp const&)>(gr0),
                                sol::resolve<scalar3dp(mvec3dp_e const&)>(gr0),
                                sol::resolve<scalar3dp(mvec3dp const&)>(gr0),
                                sol::resolve<scalar3dp(dualnum3dp const&)>(gr0)));

    lua.set_function("gr1", sol::overload(
                                // EGA grade 1 operations
                                sol::resolve<vec2d(mvec2d const&)>(gr1),
                                sol::resolve<vec3d(mvec3d_u const&)>(gr1),
                                sol::resolve<vec3d(mvec3d const&)>(gr1),
                                // PGA grade 1 operations
                                sol::resolve<vec2dp(mvec2dp_u const&)>(gr1),
                                sol::resolve<vec2dp(mvec2dp const&)>(gr1),
                                sol::resolve<vec3dp(mvec3dp_u const&)>(gr1),
                                sol::resolve<vec3dp(mvec3dp const&)>(gr1)));

    lua.set_function("gr2", sol::overload(
                                // EGA grade 2 operations
                                sol::resolve<pscalar2d(mvec2d_e const&)>(gr2),
                                sol::resolve<pscalar2d(mvec2d const&)>(gr2),
                                sol::resolve<bivec3d(mvec3d_e const&)>(gr2),
                                sol::resolve<bivec3d(mvec3d const&)>(gr2),
                                // PGA grade 2 operations
                                sol::resolve<bivec2dp(mvec2dp_e const&)>(gr2),
                                sol::resolve<bivec2dp(mvec2dp const&)>(gr2),
                                sol::resolve<bivec3dp(mvec3dp_e const&)>(gr2),
                                sol::resolve<bivec3dp(mvec3dp const&)>(gr2)));

    lua.set_function("gr3", sol::overload(
                                // EGA grade 3 operations
                                sol::resolve<pscalar3d(mvec3d_u const&)>(gr3),
                                sol::resolve<pscalar3d(mvec3d const&)>(gr3),
                                // PGA 2DP grade 3 operations
                                sol::resolve<pscalar2dp(mvec2dp_u const&)>(gr3),
                                sol::resolve<pscalar2dp(mvec2dp const&)>(gr3),
                                sol::resolve<pscalar2dp(dualnum2dp const&)>(gr3),
                                // PGA 3DP grade 3 operations
                                sol::resolve<trivec3dp(mvec3dp_u const&)>(gr3),
                                sol::resolve<trivec3dp(mvec3dp const&)>(gr3)));

    // PGA-specific grade 4 operation (only for 3DP)
    lua.set_function("gr4", sol::overload(
                                // PGA grade 4 operations
                                sol::resolve<pscalar3dp(mvec3dp_e const&)>(gr4),
                                sol::resolve<pscalar3dp(mvec3dp const&)>(gr4),
                                sol::resolve<pscalar3dp(dualnum3dp const&)>(gr4)));


    ////////////////////////////////////////////////////////////////////////////////
    // grade operations on multivectors (gr_inv, rev, inv)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("gr_inv", sol::overload(
                                   // EGA grade inversion operations
                                   sol::resolve<scalar2d(scalar2d)>(gr_inv),
                                   sol::resolve<vec2d(vec2d const&)>(gr_inv),
                                   sol::resolve<pscalar2d(pscalar2d)>(gr_inv),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(gr_inv),
                                   sol::resolve<mvec2d(mvec2d const&)>(gr_inv),
                                   sol::resolve<scalar3d(scalar3d)>(gr_inv),
                                   sol::resolve<vec3d(vec3d const&)>(gr_inv),
                                   sol::resolve<bivec3d(bivec3d const&)>(gr_inv),
                                   sol::resolve<pscalar3d(pscalar3d)>(gr_inv),
                                   sol::resolve<mvec3d_e(mvec3d_e const&)>(gr_inv),
                                   sol::resolve<mvec3d_u(mvec3d_u const&)>(gr_inv),
                                   sol::resolve<mvec3d(mvec3d const&)>(gr_inv),
                                   // PGA grade inversion operations
                                   sol::resolve<scalar2dp(scalar2dp)>(gr_inv),
                                   sol::resolve<vec2dp(vec2dp const&)>(gr_inv),
                                   sol::resolve<bivec2dp(bivec2dp const&)>(gr_inv),
                                   sol::resolve<pscalar2dp(pscalar2dp)>(gr_inv),
                                   sol::resolve<mvec2dp_e(mvec2dp_e const&)>(gr_inv),
                                   sol::resolve<mvec2dp_u(mvec2dp_u const&)>(gr_inv),
                                   sol::resolve<mvec2dp(mvec2dp const&)>(gr_inv),
                                   sol::resolve<scalar3dp(scalar3dp)>(gr_inv),
                                   sol::resolve<vec3dp(vec3dp const&)>(gr_inv),
                                   sol::resolve<bivec3dp(bivec3dp const&)>(gr_inv),
                                   sol::resolve<trivec3dp(trivec3dp const&)>(gr_inv),
                                   sol::resolve<pscalar3dp(pscalar3dp)>(gr_inv),
                                   sol::resolve<mvec3dp_e(mvec3dp_e const&)>(gr_inv),
                                   sol::resolve<mvec3dp_u(mvec3dp_u const&)>(gr_inv),
                                   sol::resolve<mvec3dp(mvec3dp const&)>(gr_inv),
                                   // sta
                                   sol::resolve<scalar4ds(scalar4ds)>(gr_inv),
                                   sol::resolve<vec4ds(vec4ds const&)>(gr_inv),
                                   sol::resolve<bivec4ds(bivec4ds const&)>(gr_inv),
                                   sol::resolve<trivec4ds(trivec4ds const&)>(gr_inv),
                                   sol::resolve<pscalar4ds(pscalar4ds)>(gr_inv),
                                   sol::resolve<mvec4ds_e(mvec4ds_e const&)>(gr_inv),
                                   sol::resolve<mvec4ds_u(mvec4ds_u const&)>(gr_inv),
                                   sol::resolve<mvec4ds(mvec4ds const&)>(gr_inv)));

    lua.set_function("rev", sol::overload(
                                // EGA reverse operations
                                sol::resolve<scalar2d(scalar2d)>(rev),
                                sol::resolve<vec2d(vec2d const&)>(rev),
                                sol::resolve<pscalar2d(pscalar2d)>(rev),
                                sol::resolve<mvec2d_e(mvec2d_e const&)>(rev),
                                sol::resolve<mvec2d(mvec2d const&)>(rev),
                                sol::resolve<scalar3d(scalar3d)>(rev),
                                sol::resolve<vec3d(vec3d const&)>(rev),
                                sol::resolve<bivec3d(bivec3d const&)>(rev),
                                sol::resolve<pscalar3d(pscalar3d)>(rev),
                                sol::resolve<mvec3d_e(mvec3d_e const&)>(rev),
                                sol::resolve<mvec3d_u(mvec3d_u const&)>(rev),
                                sol::resolve<mvec3d(mvec3d const&)>(rev),
                                // PGA reverse operations
                                sol::resolve<scalar2dp(scalar2dp)>(rev),
                                sol::resolve<vec2dp(vec2dp const&)>(rev),
                                sol::resolve<bivec2dp(bivec2dp const&)>(rev),
                                sol::resolve<pscalar2dp(pscalar2dp)>(rev),
                                sol::resolve<mvec2dp_e(mvec2dp_e const&)>(rev),
                                sol::resolve<mvec2dp_u(mvec2dp_u const&)>(rev),
                                sol::resolve<mvec2dp(mvec2dp const&)>(rev),
                                sol::resolve<scalar3dp(scalar3dp)>(rev),
                                sol::resolve<vec3dp(vec3dp const&)>(rev),
                                sol::resolve<bivec3dp(bivec3dp const&)>(rev),
                                sol::resolve<trivec3dp(trivec3dp const&)>(rev),
                                sol::resolve<pscalar3dp(pscalar3dp)>(rev),
                                sol::resolve<mvec3dp_e(mvec3dp_e const&)>(rev),
                                sol::resolve<mvec3dp_u(mvec3dp_u const&)>(rev),
                                sol::resolve<mvec3dp(mvec3dp const&)>(rev),
                                // sta
                                sol::resolve<scalar4ds(scalar4ds)>(rev),
                                sol::resolve<vec4ds(vec4ds const&)>(rev),
                                sol::resolve<bivec4ds(bivec4ds const&)>(rev),
                                sol::resolve<trivec4ds(trivec4ds const&)>(rev),
                                sol::resolve<pscalar4ds(pscalar4ds)>(rev),
                                sol::resolve<mvec4ds_e(mvec4ds_e const&)>(rev),
                                sol::resolve<mvec4ds_u(mvec4ds_u const&)>(rev),
                                sol::resolve<mvec4ds(mvec4ds const&)>(rev)));

    lua.set_function("conj", sol::overload(
                                 // EGA conjugation operations
                                 sol::resolve<scalar2d(scalar2d)>(conj),
                                 sol::resolve<vec2d(vec2d const&)>(conj),
                                 sol::resolve<pscalar2d(pscalar2d)>(conj),
                                 sol::resolve<mvec2d_e(mvec2d_e const&)>(conj),
                                 sol::resolve<mvec2d(mvec2d const&)>(conj),
                                 sol::resolve<scalar3d(scalar3d)>(conj),
                                 sol::resolve<vec3d(vec3d const&)>(conj),
                                 sol::resolve<bivec3d(bivec3d const&)>(conj),
                                 sol::resolve<pscalar3d(pscalar3d)>(conj),
                                 sol::resolve<mvec3d_e(mvec3d_e const&)>(conj),
                                 sol::resolve<mvec3d_u(mvec3d_u const&)>(conj),
                                 sol::resolve<mvec3d(mvec3d const&)>(conj),
                                 // PGA conjugation operations
                                 sol::resolve<scalar2dp(scalar2dp)>(conj),
                                 sol::resolve<vec2dp(vec2dp const&)>(conj),
                                 sol::resolve<bivec2dp(bivec2dp const&)>(conj),
                                 sol::resolve<pscalar2dp(pscalar2dp)>(conj),
                                 sol::resolve<mvec2dp_e(mvec2dp_e const&)>(conj),
                                 sol::resolve<mvec2dp_u(mvec2dp_u const&)>(conj),
                                 sol::resolve<mvec2dp(mvec2dp const&)>(conj),
                                 sol::resolve<scalar3dp(scalar3dp)>(conj),
                                 sol::resolve<vec3dp(vec3dp const&)>(conj),
                                 sol::resolve<bivec3dp(bivec3dp const&)>(conj),
                                 sol::resolve<trivec3dp(trivec3dp const&)>(conj),
                                 sol::resolve<pscalar3dp(pscalar3dp)>(conj),
                                 sol::resolve<mvec3dp_e(mvec3dp_e const&)>(conj),
                                 sol::resolve<mvec3dp_u(mvec3dp_u const&)>(conj),
                                 sol::resolve<mvec3dp(mvec3dp const&)>(conj),
                                 // sta
                                 sol::resolve<scalar4ds(scalar4ds)>(conj),
                                 sol::resolve<vec4ds(vec4ds const&)>(conj),
                                 sol::resolve<bivec4ds(bivec4ds const&)>(conj),
                                 sol::resolve<trivec4ds(trivec4ds const&)>(conj),
                                 sol::resolve<pscalar4ds(pscalar4ds)>(conj),
                                 sol::resolve<mvec4ds_e(mvec4ds_e const&)>(conj),
                                 sol::resolve<mvec4ds_u(mvec4ds_u const&)>(conj),
                                 sol::resolve<mvec4ds(mvec4ds const&)>(conj)));


    ////////////////////////////////////////////////////////////////////////////////
    // scalar product (dot)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("dot",
                     sol::overload(
                         // EGA dot products
                         sol::resolve<scalar2d(scalar2d, scalar2d)>(dot),
                         sol::resolve<scalar2d(vec2d const&, vec2d const&)>(dot),
                         sol::resolve<scalar2d(pscalar2d, pscalar2d)>(dot),
                         sol::resolve<scalar2d(mvec2d const&, mvec2d const&)>(dot),
                         sol::resolve<scalar3d(scalar3d, scalar3d)>(dot),
                         sol::resolve<scalar3d(vec3d const&, vec3d const&)>(dot),
                         sol::resolve<scalar3d(bivec3d const&, bivec3d const&)>(dot),
                         sol::resolve<scalar3d(pscalar3d, pscalar3d)>(dot),
                         sol::resolve<scalar3d(mvec3d const&, mvec3d const&)>(dot),
                         // PGA dot products (basic operations)
                         sol::resolve<scalar2dp(scalar2dp, scalar2dp)>(dot),
                         sol::resolve<scalar2dp(vec2dp const&, vec2dp const&)>(dot),
                         sol::resolve<scalar2dp(bivec2dp const&, bivec2dp const&)>(dot),
                         sol::resolve<scalar2dp(pscalar2dp, pscalar2dp)>(dot),
                         sol::resolve<scalar3dp(scalar3dp, scalar3dp)>(dot),
                         sol::resolve<scalar3dp(vec3dp const&, vec3dp const&)>(dot),
                         sol::resolve<scalar3dp(bivec3dp const&, bivec3dp const&)>(dot),
                         sol::resolve<scalar3dp(trivec3dp const&, trivec3dp const&)>(dot),
                         sol::resolve<scalar3dp(pscalar3dp, pscalar3dp)>(dot),
                         // sta
                         sol::resolve<scalar4ds(mvec4ds const&, mvec4ds const&)>(dot),
                         sol::resolve<scalar4ds(mvec4ds_e const&, mvec4ds_e const&)>(dot),
                         sol::resolve<scalar4ds(mvec4ds_u const&, mvec4ds_u const&)>(dot),
                         sol::resolve<scalar4ds(pscalar4ds, pscalar4ds)>(dot),
                         sol::resolve<scalar4ds(trivec4ds const&, trivec4ds const&)>(dot),
                         sol::resolve<scalar4ds(bivec4ds const&, bivec4ds const&)>(dot),
                         sol::resolve<scalar4ds(vec4ds const&, vec4ds const&)>(dot),
                         sol::resolve<scalar4ds(scalar4ds, scalar4ds)>(dot)));


    ////////////////////////////////////////////////////////////////////////////////
    // wedge product (= outer product) (wdg)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "wdg", sol::overload(
                   // EGA 2D wedge products
                   sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(wdg),
                   sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(wdg),
                   sol::resolve<mvec2d(mvec2d_e const&, mvec2d const&)>(wdg),
                   sol::resolve<pscalar2d(mvec2d const&, pscalar2d)>(wdg),
                   sol::resolve<pscalar2d(pscalar2d, mvec2d const&)>(wdg),
                   sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(wdg),
                   sol::resolve<mvec2d(vec2d const&, mvec2d const&)>(wdg),
                   sol::resolve<mvec2d(mvec2d const&, scalar2d)>(wdg),
                   sol::resolve<mvec2d(scalar2d, mvec2d const&)>(wdg),
                   sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(wdg),
                   sol::resolve<pscalar2d(mvec2d_e const&, pscalar2d)>(wdg),
                   sol::resolve<pscalar2d(pscalar2d, mvec2d_e const&)>(wdg),
                   sol::resolve<vec2d(mvec2d_e const&, vec2d const&)>(wdg),
                   sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(wdg),
                   sol::resolve<mvec2d_e(mvec2d_e const&, scalar2d)>(wdg),
                   sol::resolve<mvec2d_e(scalar2d, mvec2d_e const&)>(wdg),
                   sol::resolve<scalar2d(pscalar2d, pscalar2d)>(wdg),
                   sol::resolve<scalar2d(pscalar2d, vec2d const&)>(wdg),
                   sol::resolve<scalar2d(vec2d const&, pscalar2d)>(wdg),
                   sol::resolve<pscalar2d(pscalar2d, scalar2d)>(wdg),
                   sol::resolve<pscalar2d(scalar2d, pscalar2d)>(wdg),
                   sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(wdg),
                   sol::resolve<vec2d(vec2d const&, scalar2d)>(wdg),
                   sol::resolve<vec2d(scalar2d, vec2d const&)>(wdg),
                   sol::resolve<scalar2d(scalar2d, scalar2d)>(wdg),
                   // EGA 3D wedge products
                   sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(wdg),
                   sol::resolve<mvec3d(mvec3d const&, bivec3d const&)>(wdg),
                   sol::resolve<mvec3d(bivec3d const&, mvec3d const&)>(wdg),
                   sol::resolve<mvec3d(mvec3d const&, vec3d const&)>(wdg),
                   sol::resolve<mvec3d(vec3d const&, mvec3d const&)>(wdg),
                   sol::resolve<mvec3d(mvec3d const&, scalar3d)>(wdg),
                   sol::resolve<mvec3d(scalar3d, mvec3d const&)>(wdg),
                   sol::resolve<scalar3d(pscalar3d, pscalar3d)>(wdg),
                   sol::resolve<scalar3d(pscalar3d, bivec3d const&)>(wdg),
                   sol::resolve<scalar3d(bivec3d const&, pscalar3d)>(wdg),
                   sol::resolve<scalar3d(pscalar3d, vec3d const&)>(wdg),
                   sol::resolve<scalar3d(vec3d const&, pscalar3d)>(wdg),
                   sol::resolve<pscalar3d(pscalar3d, scalar3d)>(wdg),
                   sol::resolve<pscalar3d(scalar3d, pscalar3d)>(wdg),
                   sol::resolve<scalar3d(bivec3d const&, bivec3d const&)>(wdg),
                   sol::resolve<pscalar3d(bivec3d const&, vec3d const&)>(wdg),
                   sol::resolve<pscalar3d(vec3d const&, bivec3d const&)>(wdg),
                   sol::resolve<bivec3d(bivec3d const&, scalar3d)>(wdg),
                   sol::resolve<bivec3d(scalar3d, bivec3d const&)>(wdg),
                   sol::resolve<bivec3d(vec3d const&, vec3d const&)>(wdg),
                   sol::resolve<vec3d(vec3d const&, scalar3d)>(wdg),
                   sol::resolve<vec3d(scalar3d, vec3d const&)>(wdg),
                   sol::resolve<scalar3d(scalar3d, scalar3d)>(wdg),
                   // PGA 2DP wedge products
                   sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(wdg),
                   sol::resolve<mvec2dp_e(mvec2dp_e const&, mvec2dp_e const&)>(wdg),
                   sol::resolve<bivec2dp(mvec2dp_u const&, mvec2dp_u const&)>(wdg),
                   sol::resolve<scalar2dp(pscalar2dp, pscalar2dp)>(wdg),
                   sol::resolve<scalar2dp(pscalar2dp, bivec2dp const&)>(wdg),
                   sol::resolve<scalar2dp(bivec2dp const&, pscalar2dp)>(wdg),
                   sol::resolve<scalar2dp(pscalar2dp, vec2dp const&)>(wdg),
                   sol::resolve<scalar2dp(vec2dp const&, pscalar2dp)>(wdg),
                   sol::resolve<pscalar2dp(pscalar2dp, scalar2dp)>(wdg),
                   sol::resolve<pscalar2dp(scalar2dp, pscalar2dp)>(wdg),
                   sol::resolve<scalar2dp(bivec2dp const&, bivec2dp const&)>(wdg),
                   sol::resolve<pscalar2dp(bivec2dp const&, vec2dp const&)>(wdg),
                   sol::resolve<pscalar2dp(vec2dp const&, bivec2dp const&)>(wdg),
                   sol::resolve<bivec2dp(bivec2dp const&, scalar2dp)>(wdg),
                   sol::resolve<bivec2dp(scalar2dp, bivec2dp const&)>(wdg),
                   sol::resolve<bivec2dp(vec2dp const&, vec2dp const&)>(wdg),
                   sol::resolve<vec2dp(vec2dp const&, scalar2dp)>(wdg),
                   sol::resolve<vec2dp(scalar2dp, vec2dp const&)>(wdg),
                   sol::resolve<scalar2dp(scalar2dp, scalar2dp)>(wdg),
                   // Geometric alias functions pga2dp
                   sol::resolve<line2d(point2d const&, point2d const&)>(wdg),
                   // PGA 3DP wedge products
                   sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(wdg),
                   sol::resolve<mvec3dp_e(mvec3dp_e const&, mvec3dp_e const&)>(wdg),
                   sol::resolve<mvec3dp_e(mvec3dp_u const&, mvec3dp_u const&)>(wdg),
                   sol::resolve<scalar3dp(pscalar3dp, pscalar3dp)>(wdg),
                   sol::resolve<scalar3dp(pscalar3dp, trivec3dp const&)>(wdg),
                   sol::resolve<scalar3dp(trivec3dp const&, pscalar3dp)>(wdg),
                   sol::resolve<scalar3dp(pscalar3dp, bivec3dp const&)>(wdg),
                   sol::resolve<scalar3dp(bivec3dp const&, pscalar3dp)>(wdg),
                   sol::resolve<scalar3dp(pscalar3dp, vec3dp const&)>(wdg),
                   sol::resolve<scalar3dp(vec3dp const&, pscalar3dp)>(wdg),
                   sol::resolve<pscalar3dp(pscalar3dp, scalar3dp)>(wdg),
                   sol::resolve<pscalar3dp(scalar3dp, pscalar3dp)>(wdg),
                   sol::resolve<scalar3dp(trivec3dp const&, trivec3dp const&)>(wdg),
                   sol::resolve<scalar3dp(trivec3dp const&, bivec3dp const&)>(wdg),
                   sol::resolve<scalar3dp(bivec3dp const&, trivec3dp const&)>(wdg),
                   sol::resolve<pscalar3dp(trivec3dp const&, vec3dp const&)>(wdg),
                   sol::resolve<pscalar3dp(vec3dp const&, trivec3dp const&)>(wdg),
                   sol::resolve<trivec3dp(trivec3dp const&, scalar3dp)>(wdg),
                   sol::resolve<trivec3dp(scalar3dp, trivec3dp const&)>(wdg),
                   sol::resolve<pscalar3dp(bivec3dp const&, bivec3dp const&)>(wdg),
                   sol::resolve<trivec3dp(bivec3dp const&, vec3dp const&)>(wdg),
                   sol::resolve<trivec3dp(vec3dp const&, bivec3dp const&)>(wdg),
                   sol::resolve<bivec3dp(bivec3dp const&, scalar3dp)>(wdg),
                   sol::resolve<bivec3dp(scalar3dp, bivec3dp const&)>(wdg),
                   sol::resolve<bivec3dp(vec3dp const&, vec3dp const&)>(wdg),
                   sol::resolve<vec3dp(vec3dp const&, scalar3dp)>(wdg),
                   sol::resolve<vec3dp(scalar3dp, vec3dp const&)>(wdg),
                   sol::resolve<scalar3dp(scalar3dp, scalar3dp)>(wdg),
                   // Geometric alias functions pga3dp
                   sol::resolve<plane3d(line3d const&, point3d const&)>(wdg),
                   sol::resolve<plane3d(point3d const&, line3d const&)>(wdg),
                   sol::resolve<line3d(point3d const&, point3d const&)>(wdg),
                   // sta
                   sol::resolve<mvec4ds(mvec4ds const&, mvec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds const&, mvec4ds_e const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds_e const&, mvec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds const&, mvec4ds_u const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds_u const&, mvec4ds const&)>(wdg),
                   sol::resolve<pscalar4ds(mvec4ds const&, pscalar4ds)>(wdg),
                   sol::resolve<pscalar4ds(pscalar4ds, mvec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds const&, trivec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(trivec4ds const&, mvec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds const&, bivec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(bivec4ds const&, mvec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds const&, vec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(vec4ds const&, mvec4ds const&)>(wdg),
                   sol::resolve<mvec4ds(mvec4ds const&, scalar4ds)>(wdg),
                   sol::resolve<mvec4ds(scalar4ds, mvec4ds const&)>(wdg),
                   sol::resolve<mvec4ds_e(mvec4ds_e const&, mvec4ds_e const&)>(wdg),
                   sol::resolve<mvec4ds_u(mvec4ds_e const&, mvec4ds_u const&)>(wdg),
                   sol::resolve<mvec4ds_u(mvec4ds_u const&, mvec4ds_e const&)>(wdg),
                   sol::resolve<pscalar4ds(mvec4ds_e const&, pscalar4ds)>(wdg),
                   sol::resolve<pscalar4ds(pscalar4ds, mvec4ds_e const&)>(wdg),
                   sol::resolve<trivec4ds(mvec4ds_e const&, trivec4ds const&)>(wdg),
                   sol::resolve<trivec4ds(trivec4ds const&, mvec4ds_e const&)>(wdg),
                   sol::resolve<mvec4ds_e(mvec4ds_e const&, bivec4ds const&)>(wdg),
                   sol::resolve<mvec4ds_e(bivec4ds const&, mvec4ds_e const&)>(wdg),
                   sol::resolve<mvec4ds_u(mvec4ds_e const&, vec4ds const&)>(wdg),
                   sol::resolve<mvec4ds_u(vec4ds const&, mvec4ds_e const&)>(wdg),
                   sol::resolve<mvec4ds_e(mvec4ds_e const&, scalar4ds)>(wdg),
                   sol::resolve<mvec4ds_e(scalar4ds, mvec4ds_e const&)>(wdg),
                   sol::resolve<mvec4ds_e(mvec4ds_u const&, mvec4ds_u const&)>(wdg),
                   sol::resolve<scalar4ds(mvec4ds_u const&, pscalar4ds)>(wdg),
                   sol::resolve<scalar4ds(pscalar4ds, mvec4ds_u const&)>(wdg),
                   sol::resolve<pscalar4ds(mvec4ds_u const&, trivec4ds const&)>(wdg),
                   sol::resolve<pscalar4ds(trivec4ds const&, mvec4ds_u const&)>(wdg),
                   sol::resolve<trivec4ds(mvec4ds_u const&, bivec4ds const&)>(wdg),
                   sol::resolve<trivec4ds(bivec4ds const&, mvec4ds_u const&)>(wdg),
                   sol::resolve<mvec4ds_e(mvec4ds_u const&, vec4ds const&)>(wdg),
                   sol::resolve<mvec4ds_e(vec4ds const&, mvec4ds_u const&)>(wdg),
                   sol::resolve<mvec4ds_u(mvec4ds_u const&, scalar4ds)>(wdg),
                   sol::resolve<mvec4ds_u(scalar4ds, mvec4ds_u const&)>(wdg),
                   sol::resolve<scalar4ds(pscalar4ds, pscalar4ds)>(wdg),
                   sol::resolve<scalar4ds(pscalar4ds, trivec4ds const&)>(wdg),
                   sol::resolve<scalar4ds(trivec4ds const&, pscalar4ds)>(wdg),
                   sol::resolve<scalar4ds(pscalar4ds, bivec4ds const&)>(wdg),
                   sol::resolve<scalar4ds(bivec4ds const&, pscalar4ds)>(wdg),
                   sol::resolve<scalar4ds(pscalar4ds, vec4ds const&)>(wdg),
                   sol::resolve<scalar4ds(vec4ds const&, pscalar4ds)>(wdg),
                   sol::resolve<pscalar4ds(pscalar4ds, scalar4ds)>(wdg),
                   sol::resolve<pscalar4ds(scalar4ds, pscalar4ds)>(wdg),
                   sol::resolve<scalar4ds(trivec4ds const&, trivec4ds const&)>(wdg),
                   sol::resolve<scalar4ds(trivec4ds const&, bivec4ds const&)>(wdg),
                   sol::resolve<scalar4ds(bivec4ds const&, trivec4ds const&)>(wdg),
                   sol::resolve<pscalar4ds(trivec4ds const&, vec4ds const&)>(wdg),
                   sol::resolve<pscalar4ds(vec4ds const&, trivec4ds const&)>(wdg),
                   sol::resolve<trivec4ds(trivec4ds const&, scalar4ds)>(wdg),
                   sol::resolve<trivec4ds(scalar4ds, trivec4ds const&)>(wdg),
                   sol::resolve<pscalar4ds(bivec4ds const&, bivec4ds const&)>(wdg),
                   sol::resolve<trivec4ds(bivec4ds const&, vec4ds const&)>(wdg),
                   sol::resolve<trivec4ds(vec4ds const&, bivec4ds const&)>(wdg),
                   sol::resolve<bivec4ds(bivec4ds const&, scalar4ds)>(wdg),
                   sol::resolve<bivec4ds(scalar4ds, bivec4ds const&)>(wdg),
                   sol::resolve<bivec4ds(vec4ds const&, vec4ds const&)>(wdg),
                   sol::resolve<vec4ds(vec4ds const&, scalar4ds)>(wdg),
                   sol::resolve<vec4ds(scalar4ds, vec4ds const&)>(wdg),
                   sol::resolve<scalar4ds(scalar4ds, scalar4ds)>(wdg)));

    // join as convenience alias for wdg()
    lua.set_function("join",
                     sol::overload(
                         // join as convenience alias for pga2dp
                         sol::resolve<pscalar2dp(bivec2dp const&, vec2dp const&)>(join),
                         sol::resolve<pscalar2dp(vec2dp const&, bivec2dp const&)>(join),
                         sol::resolve<bivec2dp(vec2dp const&, vec2dp const&)>(join),
                         sol::resolve<line2d(point2d const&, point2d const&)>(join),
                         // join as convenience alias for pga3dp
                         sol::resolve<pscalar3dp(trivec3dp const&, vec3dp const&)>(join),
                         sol::resolve<pscalar3dp(vec3dp const&, trivec3dp const&)>(join),
                         sol::resolve<trivec3dp(bivec3dp const&, vec3dp const&)>(join),
                         sol::resolve<trivec3dp(vec3dp const&, bivec3dp const&)>(join),
                         sol::resolve<plane3d(line3d const&, point3d const&)>(join),
                         sol::resolve<plane3d(point3d const&, line3d const&)>(join),
                         sol::resolve<line3d(point3d const&, point3d const&)>(join),
                         sol::resolve<bivec3dp(vec3dp const&, vec3dp const&)>(join)));

    ////////////////////////////////////////////////////////////////////////////////
    // regressive wedge product (= outer product for complements) (rwdg)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "rwdg", sol::overload(
                    // EGA 2D regressive wedge products
                    sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(rwdg),
                    sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(rwdg),
                    sol::resolve<vec2d(pscalar2d, vec2d const&)>(rwdg),
                    sol::resolve<vec2d(vec2d const&, pscalar2d)>(rwdg),
                    sol::resolve<scalar2d(pscalar2d, scalar2d)>(rwdg),
                    sol::resolve<scalar2d(scalar2d, pscalar2d)>(rwdg),
                    sol::resolve<scalar2d(vec2d const&, vec2d const&)>(rwdg),
                    sol::resolve<scalar2d(vec2d const&, scalar2d)>(rwdg),
                    sol::resolve<scalar2d(scalar2d, vec2d const&)>(rwdg),
                    sol::resolve<scalar2d(scalar2d, scalar2d)>(rwdg),
                    // EGA 3D regressive wedge products
                    sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(rwdg),
                    sol::resolve<pscalar3d(pscalar3d, pscalar3d)>(rwdg),
                    sol::resolve<bivec3d(pscalar3d, bivec3d const&)>(rwdg),
                    sol::resolve<bivec3d(bivec3d const&, pscalar3d)>(rwdg),
                    sol::resolve<vec3d(pscalar3d, vec3d const&)>(rwdg),
                    sol::resolve<vec3d(vec3d const&, pscalar3d)>(rwdg),
                    sol::resolve<scalar3d(pscalar3d, scalar3d)>(rwdg),
                    sol::resolve<scalar3d(scalar3d, pscalar3d)>(rwdg),
                    sol::resolve<vec3d(bivec3d const&, bivec3d const&)>(rwdg),
                    sol::resolve<scalar3d(bivec3d const&, vec3d const&)>(rwdg),
                    sol::resolve<scalar3d(vec3d const&, bivec3d const&)>(rwdg),
                    sol::resolve<scalar3d(bivec3d const&, scalar3d)>(rwdg),
                    sol::resolve<scalar3d(scalar3d, bivec3d const&)>(rwdg),
                    sol::resolve<scalar3d(vec3d const&, vec3d const&)>(rwdg),
                    sol::resolve<scalar3d(vec3d const&, scalar3d)>(rwdg),
                    sol::resolve<scalar3d(scalar3d, vec3d const&)>(rwdg),
                    sol::resolve<scalar3d(scalar3d, scalar3d)>(rwdg),
                    // PGA 2DP regressive wedge products
                    sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(rwdg),
                    sol::resolve<vec2dp(mvec2dp_e const&, mvec2dp_e const&)>(rwdg),
                    sol::resolve<mvec2dp_u(mvec2dp_u const&, mvec2dp_u const&)>(rwdg),
                    //
                    sol::resolve<pscalar2dp(pscalar2dp, pscalar2dp)>(rwdg),
                    sol::resolve<bivec2dp(pscalar2dp, bivec2dp const&)>(rwdg),
                    sol::resolve<bivec2dp(bivec2dp const&, pscalar2dp)>(rwdg),
                    sol::resolve<vec2dp(pscalar2dp, vec2dp const&)>(rwdg),
                    sol::resolve<vec2dp(vec2dp const&, pscalar2dp)>(rwdg),
                    sol::resolve<scalar2dp(pscalar2dp, scalar2dp)>(rwdg),
                    sol::resolve<scalar2dp(scalar2dp, pscalar2dp)>(rwdg),
                    //
                    sol::resolve<vec2dp(bivec2dp const&, bivec2dp const&)>(rwdg),
                    sol::resolve<scalar2dp(bivec2dp const&, vec2dp const&)>(rwdg),
                    sol::resolve<scalar2dp(vec2dp const&, bivec2dp const&)>(rwdg),
                    sol::resolve<scalar2dp(bivec2dp const&, scalar2dp)>(rwdg),
                    sol::resolve<scalar2dp(scalar2dp, bivec2dp const&)>(rwdg),
                    sol::resolve<scalar2dp(vec2dp const&, vec2dp const&)>(rwdg),
                    sol::resolve<scalar2dp(vec2dp const&, scalar2dp)>(rwdg),
                    sol::resolve<scalar2dp(scalar2dp, vec2dp const&)>(rwdg),
                    sol::resolve<scalar2dp(scalar2dp, scalar2dp)>(rwdg),
                    // PGA 3DP regressive wedge products
                    sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(rwdg),
                    sol::resolve<mvec3dp_e(mvec3dp_e const&, mvec3dp_e const&)>(rwdg),
                    sol::resolve<mvec3dp_e(mvec3dp_u const&, mvec3dp_u const&)>(rwdg),
                    //
                    sol::resolve<pscalar3dp(pscalar3dp, pscalar3dp)>(rwdg),
                    sol::resolve<trivec3dp(pscalar3dp, trivec3dp const&)>(rwdg),
                    sol::resolve<trivec3dp(trivec3dp const&, pscalar3dp)>(rwdg),
                    sol::resolve<bivec3dp(pscalar3dp, bivec3dp const&)>(rwdg),
                    sol::resolve<bivec3dp(bivec3dp const&, pscalar3dp)>(rwdg),
                    sol::resolve<vec3dp(pscalar3dp, vec3dp const&)>(rwdg),
                    sol::resolve<vec3dp(vec3dp const&, pscalar3dp)>(rwdg),
                    sol::resolve<scalar3dp(pscalar3dp, scalar3dp)>(rwdg),
                    sol::resolve<scalar3dp(scalar3dp, pscalar3dp)>(rwdg),
                    //
                    sol::resolve<bivec3dp(trivec3dp const&, trivec3dp const&)>(rwdg),
                    sol::resolve<vec3dp(trivec3dp const&, bivec3dp const&)>(rwdg),
                    sol::resolve<vec3dp(bivec3dp const&, trivec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(trivec3dp const&, vec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(vec3dp const&, trivec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(trivec3dp const&, scalar3dp)>(rwdg),
                    sol::resolve<scalar3dp(scalar3dp, trivec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(bivec3dp const&, bivec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(bivec3dp const&, vec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(vec3dp const&, bivec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(vec3dp const&, vec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(vec3dp const&, scalar3dp)>(rwdg),
                    sol::resolve<scalar3dp(scalar3dp, vec3dp const&)>(rwdg),
                    sol::resolve<scalar3dp(scalar3dp, scalar3dp)>(rwdg),
                    // sta
                    sol::resolve<mvec4ds(mvec4ds const&, mvec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds const&, mvec4ds_e const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds_e const&, mvec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds const&, mvec4ds_u const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds_u const&, mvec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds const&, pscalar4ds)>(rwdg),
                    sol::resolve<mvec4ds(pscalar4ds, mvec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds const&, trivec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(trivec4ds const&, mvec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds const&, bivec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(bivec4ds const&, mvec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(mvec4ds const&, vec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds(vec4ds const&, mvec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(mvec4ds const&, scalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, mvec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds_e(mvec4ds_e const&, mvec4ds_e const&)>(rwdg),
                    sol::resolve<mvec4ds_u(mvec4ds_e const&, mvec4ds_u const&)>(rwdg),
                    sol::resolve<mvec4ds_u(mvec4ds_u const&, mvec4ds_e const&)>(rwdg),
                    sol::resolve<mvec4ds_e(mvec4ds_e const&, pscalar4ds)>(rwdg),
                    sol::resolve<mvec4ds_e(pscalar4ds, mvec4ds_e const&)>(rwdg),
                    sol::resolve<mvec4ds_u(mvec4ds_e const&, trivec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds_u(trivec4ds const&, mvec4ds_e const&)>(rwdg),
                    sol::resolve<mvec4ds_e(mvec4ds_e const&, bivec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds_e(bivec4ds const&, mvec4ds_e const&)>(rwdg),
                    sol::resolve<vec4ds(mvec4ds_e const&, vec4ds const&)>(rwdg),
                    sol::resolve<vec4ds(vec4ds const&, mvec4ds_e const&)>(rwdg),
                    sol::resolve<scalar4ds(mvec4ds_e const&, scalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, mvec4ds_e const&)>(rwdg),
                    sol::resolve<mvec4ds_e(mvec4ds_u const&, mvec4ds_u const&)>(rwdg),
                    sol::resolve<mvec4ds_u(mvec4ds_u const&, pscalar4ds)>(rwdg),
                    sol::resolve<mvec4ds_u(pscalar4ds, mvec4ds_u const&)>(rwdg),
                    sol::resolve<mvec4ds_e(mvec4ds_u const&, trivec4ds const&)>(rwdg),
                    sol::resolve<mvec4ds_e(trivec4ds const&, mvec4ds_u const&)>(rwdg),
                    sol::resolve<vec4ds(mvec4ds_u const&, bivec4ds const&)>(rwdg),
                    sol::resolve<vec4ds(bivec4ds const&, mvec4ds_u const&)>(rwdg),
                    sol::resolve<scalar4ds(mvec4ds_u const&, vec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(vec4ds const&, mvec4ds_u const&)>(rwdg),
                    sol::resolve<scalar4ds(mvec4ds_u const&, scalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, mvec4ds_u const&)>(rwdg),
                    sol::resolve<pscalar4ds(pscalar4ds, pscalar4ds)>(rwdg),
                    sol::resolve<trivec4ds(pscalar4ds, trivec4ds const&)>(rwdg),
                    sol::resolve<trivec4ds(trivec4ds const&, pscalar4ds)>(rwdg),
                    sol::resolve<bivec4ds(pscalar4ds, bivec4ds const&)>(rwdg),
                    sol::resolve<bivec4ds(bivec4ds const&, pscalar4ds)>(rwdg),
                    sol::resolve<vec4ds(pscalar4ds, vec4ds const&)>(rwdg),
                    sol::resolve<vec4ds(vec4ds const&, pscalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(pscalar4ds, scalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, pscalar4ds)>(rwdg),
                    sol::resolve<bivec4ds(trivec4ds const&, trivec4ds const&)>(rwdg),
                    sol::resolve<vec4ds(trivec4ds const&, bivec4ds const&)>(rwdg),
                    sol::resolve<vec4ds(bivec4ds const&, trivec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(trivec4ds const&, vec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(vec4ds const&, trivec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(trivec4ds const&, scalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, trivec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(bivec4ds const&, bivec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(bivec4ds const&, vec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(vec4ds const&, bivec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(bivec4ds const&, scalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, bivec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(vec4ds const&, vec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(vec4ds const&, scalar4ds)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, vec4ds const&)>(rwdg),
                    sol::resolve<scalar4ds(scalar4ds, scalar4ds)>(rwdg)));

    // meet as convenience alias for rwdg()
    lua.set_function("meet",
                     sol::overload(
                         // PGA 2DP
                         sol::resolve<vec2dp(bivec2dp const&, bivec2dp const&)>(meet),
                         sol::resolve<vec2dp(line2d const&, line2d const&)>(meet),
                         // PGA 3DP
                         sol::resolve<bivec3dp(trivec3dp const&, trivec3dp const&)>(meet),
                         sol::resolve<line3d(plane3d const&, plane3d const&)>(meet),
                         sol::resolve<vec3dp(trivec3dp const&, bivec3dp const&)>(meet),
                         sol::resolve<vec3dp(bivec3dp const&, trivec3dp const&)>(meet),
                         sol::resolve<vec3dp(plane3d const&, line3d const&)>(meet),
                         sol::resolve<vec3dp(line3d const&, plane3d const&)>(meet)));

    ////////////////////////////////////////////////////////////////////////////////
    // other products (inner, fat_dot, cmt, cross)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("cmt",
                     sol::overload(
                         // 2D commutator operations
                         sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(cmt),
                         sol::resolve<vec2d(pscalar2d, vec2d const&)>(cmt),
                         sol::resolve<vec2d(vec2d const&, pscalar2d)>(cmt),
                         sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(cmt),
                         // 3D commutator operations
                         sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(cmt),
                         sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(cmt),
                         sol::resolve<vec3d(bivec3d const&, vec3d const&)>(cmt),
                         sol::resolve<vec3d(vec3d const&, bivec3d const&)>(cmt),
                         sol::resolve<bivec3d(vec3d const&, vec3d const&)>(cmt),
                         // PGA 2DP commutator operations
                         sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(cmt),
                         sol::resolve<bivec2dp(bivec2dp const&, bivec2dp const&)>(cmt),
                         sol::resolve<vec2dp(bivec2dp const&, vec2dp const&)>(cmt),
                         sol::resolve<vec2dp(vec2dp const&, bivec2dp const&)>(cmt),
                         sol::resolve<bivec2dp(vec2dp const&, vec2dp const&)>(cmt),
                         // PGA 3DP commutator operations
                         sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(cmt),
                         sol::resolve<bivec3dp(trivec3dp const&, trivec3dp const&)>(cmt),
                         sol::resolve<trivec3dp(trivec3dp const&, bivec3dp const&)>(cmt),
                         sol::resolve<trivec3dp(bivec3dp const&, trivec3dp const&)>(cmt),
                         sol::resolve<pscalar3dp(trivec3dp const&, vec3dp const&)>(cmt),
                         sol::resolve<pscalar3dp(vec3dp const&, trivec3dp const&)>(cmt),
                         sol::resolve<bivec3dp(bivec3dp const&, bivec3dp const&)>(cmt),
                         sol::resolve<vec3dp(bivec3dp const&, vec3dp const&)>(cmt),
                         sol::resolve<vec3dp(vec3dp const&, bivec3dp const&)>(cmt),
                         sol::resolve<bivec3dp(vec3dp const&, vec3dp const&)>(cmt),
                         // sta
                         sol::resolve<mvec4ds(mvec4ds const&, mvec4ds const&)>(cmt),
                         sol::resolve<scalar4ds(pscalar4ds, pscalar4ds)>(cmt),
                         sol::resolve<vec4ds(pscalar4ds, trivec4ds const&)>(cmt),
                         sol::resolve<vec4ds(trivec4ds const&, pscalar4ds)>(cmt),
                         sol::resolve<scalar4ds(pscalar4ds, bivec4ds const&)>(cmt),
                         sol::resolve<scalar4ds(bivec4ds const&, pscalar4ds)>(cmt),
                         sol::resolve<trivec4ds(pscalar4ds, vec4ds const&)>(cmt),
                         sol::resolve<trivec4ds(vec4ds const&, pscalar4ds)>(cmt),
                         sol::resolve<pscalar4ds(pscalar4ds, scalar4ds)>(cmt),
                         sol::resolve<pscalar4ds(scalar4ds, pscalar4ds)>(cmt),
                         sol::resolve<bivec4ds(trivec4ds const&, trivec4ds const&)>(cmt),
                         sol::resolve<trivec4ds(trivec4ds const&, bivec4ds const&)>(cmt),
                         sol::resolve<trivec4ds(bivec4ds const&, trivec4ds const&)>(cmt),
                         sol::resolve<pscalar4ds(trivec4ds const&, vec4ds const&)>(cmt),
                         sol::resolve<pscalar4ds(vec4ds const&, trivec4ds const&)>(cmt),
                         sol::resolve<scalar4ds(trivec4ds const&, scalar4ds)>(cmt),
                         sol::resolve<scalar4ds(scalar4ds, trivec4ds const&)>(cmt),
                         sol::resolve<bivec4ds(bivec4ds const&, bivec4ds const&)>(cmt),
                         sol::resolve<vec4ds(bivec4ds const&, vec4ds const&)>(cmt),
                         sol::resolve<vec4ds(vec4ds const&, bivec4ds const&)>(cmt),
                         sol::resolve<scalar4ds(bivec4ds const&, scalar4ds)>(cmt),
                         sol::resolve<scalar4ds(scalar4ds, bivec4ds const&)>(cmt),
                         sol::resolve<bivec4ds(vec4ds const&, vec4ds const&)>(cmt),
                         sol::resolve<scalar4ds(vec4ds const&, scalar4ds)>(cmt),
                         sol::resolve<scalar4ds(scalar4ds, vec4ds const&)>(cmt),
                         sol::resolve<scalar4ds(scalar4ds, scalar4ds)>(cmt)));

    lua.set_function("cross", sol::resolve<vec3d(vec3d const&, vec3d const&)>(cross));

    ////////////////////////////////////////////////////////////////////////////////
    // attitude and support for pga
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("att", sol::overload(
                                // PGA 2DP attitude
                                sol::resolve<scalar2dp(vec2dp const&)>(att),
                                sol::resolve<vec2dp(bivec2dp const&)>(att),
                                sol::resolve<bivec2dp(pscalar2dp)>(att),
                                // PGA 3DP attitude
                                sol::resolve<scalar3dp(vec3dp const&)>(att),
                                sol::resolve<vec3dp(bivec3dp const&)>(att),
                                sol::resolve<bivec3dp(trivec3dp const&)>(att),
                                sol::resolve<trivec3dp(pscalar3dp)>(att)));

    lua.set_function("sup", sol::overload(
                                // PGA 2DP support
                                sol::resolve<vec2dp(bivec2dp const&)>(sup),
                                // PGA 3DP support
                                sol::resolve<vec3dp(bivec3dp const&)>(sup),
                                sol::resolve<vec3dp(trivec3dp const&)>(sup)));

    ////////////////////////////////////////////////////////////////////////////////
    // angles and rotations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("angle",
                     sol::overload( // EGA angle functions
                         sol::resolve<value_t(vec2d const&, vec2d const&)>(angle),
                         sol::resolve<value_t(vec3d const&, vec3d const&)>(angle),
                         sol::resolve<value_t(bivec3d const&, bivec3d const&)>(angle),
                         sol::resolve<value_t(vec3d const&, bivec3d const&)>(angle),
                         sol::resolve<value_t(bivec3d const&, vec3d const&)>(angle),
                         // PGA angle functions
                         sol::resolve<value_t(vec2dp const&, vec2dp const&)>(angle),
                         sol::resolve<value_t(bivec2dp const&, bivec2dp const&)>(angle),
                         sol::resolve<value_t(vec3dp const&, vec3dp const&)>(angle),
                         sol::resolve<value_t(bivec3dp const&, bivec3dp const&)>(angle),
                         sol::resolve<value_t(trivec3dp const&, bivec3dp const&)>(angle),
                         sol::resolve<value_t(bivec3dp const&, trivec3dp const&)>(angle),
                         sol::resolve<value_t(trivec3dp const&, trivec3dp const&)>(angle),
                         // sta
                         sol::resolve<value_t(vec4ds const&, vec4ds const&)>(angle)));

    lua.set_function("angle_to_re", sol::resolve<value_t(mvec2d_e const&)>(angle_to_re));

    lua.set_function("exp",
                     sol::overload(sol::resolve<mvec2d_e(pscalar2d, value_t)>(exp),
                                   sol::resolve<mvec3d_e(bivec3d const&, value_t)>(exp),
                                   // sta
                                   sol::resolve<mvec4ds_e(bivec4ds const&)>(exp)));

    lua.set_function(
        "get_rotor",
        sol::overload(sol::resolve<mvec2d_e(pscalar2d, value_t)>(get_rotor),
                      sol::resolve<mvec3d_e(bivec3d const&, value_t)>(get_rotor),
                      // sta
                      sol::resolve<mvec4ds_e(bivec4ds const&, value_t)>(get_rotor)));

    lua.set_function(
        "rotate",
        sol::overload(sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(rotate),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(rotate),
                      sol::resolve<vec3d(vec3d const&, mvec3d_e const&)>(rotate),
                      sol::resolve<bivec3d(bivec3d const&, mvec3d_e const&)>(rotate),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d_e const&)>(rotate)));

    // PGA motor functions (analogue to EGA get_rotor)
    lua.set_function("get_motor",
                     sol::overload(
                         // PGA 2DP motor functions
                         sol::resolve<mvec2dp_u(vec2dp const&, value_t)>(get_motor),
                         sol::resolve<mvec2dp_u(vec2dp const&)>(get_motor),
                         // PGA 3DP motor functions
                         sol::resolve<mvec3dp_e(bivec3dp const&, value_t)>(get_motor),
                         sol::resolve<mvec3dp_e(vec3dp const&)>(get_motor)));

    lua.set_function(
        "get_motor_from_lines",
        sol::resolve<mvec2dp_u(bivec2dp const&, bivec2dp const&)>(get_motor_from_lines));

    lua.set_function("get_motor_from_planes",
                     sol::resolve<mvec3dp_e(trivec3dp const&, trivec3dp const&)>(
                         get_motor_from_planes));

    // PGA movement functions (analogue to EGA rotate)
    lua.set_function(
        "move2dp",
        sol::overload(
            sol::resolve<vec2dp(vec2dp const&, mvec2dp_u const&)>(move2dp),
            sol::resolve<bivec2dp(bivec2dp const&, mvec2dp_u const&)>(move2dp)));

    lua.set_function(
        "move3dp",
        sol::overload(
            sol::resolve<vec3dp(vec3dp const&, mvec3dp_e const&)>(move3dp),
            sol::resolve<bivec3dp(bivec3dp const&, mvec3dp_e const&)>(move3dp),
            sol::resolve<trivec3dp(trivec3dp const&, mvec3dp_e const&)>(move3dp)));

    ////////////////////////////////////////////////////////////////////////////////
    // dualization operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("dual", sol::overload(sol::resolve<scalar3d(pscalar3d)>(dual),
                                           sol::resolve<vec3d(bivec3d const&)>(dual),
                                           sol::resolve<bivec3d(vec3d const&)>(dual),
                                           sol::resolve<pscalar3d(scalar3d)>(dual),
                                           sol::resolve<mvec3d_u(mvec3d_e const&)>(dual),
                                           sol::resolve<mvec3d_e(mvec3d_u const&)>(dual),
                                           sol::resolve<mvec3d(mvec3d const&)>(dual)));

    lua.set_function("l_dual",
                     sol::overload(sol::resolve<scalar2d(pscalar2d)>(l_dual),
                                   sol::resolve<vec2d(vec2d const&)>(l_dual),
                                   sol::resolve<pscalar2d(scalar2d)>(l_dual),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(l_dual),
                                   sol::resolve<mvec2d(mvec2d const&)>(l_dual),
                                   // sta
                                   sol::resolve<pscalar4ds(scalar4ds)>(l_dual),
                                   sol::resolve<trivec4ds(vec4ds const&)>(l_dual),
                                   sol::resolve<bivec4ds(bivec4ds const&)>(l_dual),
                                   sol::resolve<vec4ds(trivec4ds const&)>(l_dual),
                                   sol::resolve<scalar4ds(pscalar4ds)>(l_dual),
                                   sol::resolve<mvec4ds_e(mvec4ds_e const&)>(l_dual),
                                   sol::resolve<mvec4ds_u(mvec4ds_u const&)>(l_dual),
                                   sol::resolve<mvec4ds(mvec4ds const&)>(l_dual)));

    lua.set_function("r_dual",
                     sol::overload(sol::resolve<scalar2d(pscalar2d)>(r_dual),
                                   sol::resolve<vec2d(vec2d const&)>(r_dual),
                                   sol::resolve<pscalar2d(scalar2d)>(r_dual),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(r_dual),
                                   sol::resolve<mvec2d(mvec2d const&)>(r_dual),
                                   // sta
                                   sol::resolve<pscalar4ds(scalar4ds)>(r_dual),
                                   sol::resolve<trivec4ds(vec4ds const&)>(r_dual),
                                   sol::resolve<bivec4ds(bivec4ds const&)>(r_dual),
                                   sol::resolve<vec4ds(trivec4ds const&)>(r_dual),
                                   sol::resolve<scalar4ds(pscalar4ds)>(r_dual),
                                   sol::resolve<mvec4ds_e(mvec4ds_e const&)>(r_dual),
                                   sol::resolve<mvec4ds_u(mvec4ds_u const&)>(r_dual),
                                   sol::resolve<mvec4ds(mvec4ds const&)>(r_dual)));

    ////////////////////////////////////////////////////////////////////////////////
    // complement operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("r_cmpl",
                     sol::overload( // EGA 2D r_cmpl functions
                         sol::resolve<pscalar2d(scalar2d)>(r_cmpl),
                         sol::resolve<vec2d(vec2d const&)>(r_cmpl),
                         sol::resolve<scalar2d(pscalar2d)>(r_cmpl),
                         sol::resolve<mvec2d_e(mvec2d_e const&)>(r_cmpl),
                         sol::resolve<mvec2d(mvec2d const&)>(r_cmpl),
                         // PGA 3DP r_cmpl functions (8 total)
                         sol::resolve<pscalar3dp(scalar3dp)>(r_cmpl),
                         sol::resolve<trivec3dp(vec3dp const&)>(r_cmpl),
                         sol::resolve<bivec3dp(bivec3dp const&)>(r_cmpl),
                         sol::resolve<vec3dp(trivec3dp const&)>(r_cmpl),
                         sol::resolve<scalar3dp(pscalar3dp)>(r_cmpl),
                         sol::resolve<mvec3dp_u(mvec3dp_u const&)>(r_cmpl),
                         sol::resolve<mvec3dp_e(mvec3dp_e const&)>(r_cmpl),
                         sol::resolve<mvec3dp(mvec3dp const&)>(r_cmpl),
                         // sta
                         sol::resolve<pscalar4ds(scalar4ds)>(r_cmpl),
                         sol::resolve<trivec4ds(vec4ds const&)>(r_cmpl),
                         sol::resolve<bivec4ds(bivec4ds const&)>(r_cmpl),
                         sol::resolve<vec4ds(trivec4ds const&)>(r_cmpl),
                         sol::resolve<scalar4ds(pscalar4ds)>(r_cmpl),
                         sol::resolve<mvec4ds_u(mvec4ds_u const&)>(r_cmpl),
                         sol::resolve<mvec4ds_e(mvec4ds_e const&)>(r_cmpl),
                         sol::resolve<mvec4ds(mvec4ds const&)>(r_cmpl)));

    lua.set_function("l_cmpl",
                     sol::overload( // EGA 2D l_cmpl functions
                         sol::resolve<pscalar2d(scalar2d)>(l_cmpl),
                         sol::resolve<vec2d(vec2d const&)>(l_cmpl),
                         sol::resolve<scalar2d(pscalar2d)>(l_cmpl),
                         sol::resolve<mvec2d_e(mvec2d_e const&)>(l_cmpl),
                         sol::resolve<mvec2d(mvec2d const&)>(l_cmpl),
                         // PGA 3DP l_cmpl functions (8 total)
                         sol::resolve<pscalar3dp(scalar3dp)>(l_cmpl),
                         sol::resolve<trivec3dp(vec3dp const&)>(l_cmpl),
                         sol::resolve<bivec3dp(bivec3dp const&)>(l_cmpl),
                         sol::resolve<vec3dp(trivec3dp const&)>(l_cmpl),
                         sol::resolve<scalar3dp(pscalar3dp)>(l_cmpl),
                         sol::resolve<mvec3dp_u(mvec3dp_u const&)>(l_cmpl),
                         sol::resolve<mvec3dp_e(mvec3dp_e const&)>(l_cmpl),
                         sol::resolve<mvec3dp(mvec3dp const&)>(l_cmpl),
                         // sta
                         sol::resolve<pscalar4ds(scalar4ds)>(l_cmpl),
                         sol::resolve<trivec4ds(vec4ds const&)>(l_cmpl),
                         sol::resolve<bivec4ds(bivec4ds const&)>(l_cmpl),
                         sol::resolve<vec4ds(trivec4ds const&)>(l_cmpl),
                         sol::resolve<scalar4ds(pscalar4ds)>(l_cmpl),
                         sol::resolve<mvec4ds_u(mvec4ds_u const&)>(l_cmpl),
                         sol::resolve<mvec4ds_e(mvec4ds_e const&)>(l_cmpl),
                         sol::resolve<mvec4ds(mvec4ds const&)>(l_cmpl)));

    lua.set_function("cmpl", sol::overload( // EGA 3D cmpl functions
                                 sol::resolve<pscalar3d(scalar3d)>(cmpl),
                                 sol::resolve<bivec3d(vec3d const&)>(cmpl),
                                 sol::resolve<vec3d(bivec3d const&)>(cmpl),
                                 sol::resolve<scalar3d(pscalar3d)>(cmpl),
                                 sol::resolve<mvec3d_u(mvec3d_e const&)>(cmpl),
                                 sol::resolve<mvec3d_e(mvec3d_u const&)>(cmpl),
                                 sol::resolve<mvec3d(mvec3d const&)>(cmpl),
                                 // PGA 2DP cmpl functions (7 total)
                                 sol::resolve<pscalar2dp(scalar2dp)>(cmpl),
                                 sol::resolve<bivec2dp(vec2dp const&)>(cmpl),
                                 sol::resolve<vec2dp(bivec2dp const&)>(cmpl),
                                 sol::resolve<scalar2dp(pscalar2dp)>(cmpl),
                                 sol::resolve<mvec2dp_u(mvec2dp_e const&)>(cmpl),
                                 sol::resolve<mvec2dp_e(mvec2dp_u const&)>(cmpl),
                                 sol::resolve<mvec2dp(mvec2dp const&)>(cmpl)));

    ////////////////////////////////////////////////////////////////////////////////
    // projections, rejections and reflections
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "project_onto",
        sol::overload( // EGA projection functions
            sol::resolve<vec2d(vec2d const&, vec2d const&)>(project_onto),
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(project_onto),
            sol::resolve<vec3d(vec3d const&, bivec3d const&)>(project_onto),
            // PGA 2DP projection functions
            sol::resolve<vec2dp(vec2dp const&, vec2dp const&)>(project_onto),
            sol::resolve<vec2dp(vec2dp const&, bivec2dp const&)>(project_onto),
            // PGA 3DP projection functions
            sol::resolve<vec3dp(vec3dp const&, vec3dp const&)>(project_onto),
            sol::resolve<vec3dp(vec3dp const&, bivec3dp const&)>(project_onto),
            sol::resolve<vec3dp(vec3dp const&, trivec3dp const&)>(project_onto),
            // sta
            sol::resolve<vec4ds(vec4ds const&, vec4ds const&)>(project_onto),
            sol::resolve<vec4ds(vec4ds const&, bivec4ds const&)>(project_onto)));

    lua.set_function(
        "reject_from",
        sol::overload( // EGA rejection functions
            sol::resolve<vec2d(vec2d const&, vec2d const&)>(reject_from),
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(reject_from),
            sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reject_from),
            // PGA 2DP rejection functions
            sol::resolve<vec2dp(vec2dp const&, vec2dp const&)>(reject_from),
            sol::resolve<vec2dp(vec2dp const&, bivec2dp const&)>(reject_from),
            // PGA 3DP rejection functions
            sol::resolve<vec3dp(vec3dp const&, vec3dp const&)>(reject_from),
            sol::resolve<vec3dp(vec3dp const&, bivec3dp const&)>(reject_from),
            sol::resolve<vec3dp(vec3dp const&, trivec3dp const&)>(reject_from),
            // sta
            sol::resolve<vec4ds(vec4ds const&, vec4ds const&)>(reject_from),
            sol::resolve<vec4ds(vec4ds const&, bivec4ds const&)>(reject_from)));

    lua.set_function(
        "reflect_on",
        sol::overload(
            // EGA 2D reflection functions
            sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on),
            // EGA 3D reflection functions
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on),
            sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reflect_on),
            sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(reflect_on),
            // PGA 2DP reflection functions
            sol::resolve<vec2dp(vec2dp const&, bivec2dp const&)>(reflect_on),
            sol::resolve<bivec2dp(bivec2dp const&, bivec2dp const&)>(reflect_on),
            // PGA 3DP reflection functions
            sol::resolve<vec3dp(vec3dp const&, trivec3dp const&)>(reflect_on),
            sol::resolve<bivec3dp(bivec3dp const&, trivec3dp const&)>(reflect_on),
            sol::resolve<trivec3dp(trivec3dp const&, trivec3dp const&)>(reflect_on),
            // sta
            sol::resolve<vec4ds(vec4ds const&, vec4ds const&)>(reflect_on),
            sol::resolve<vec4ds(vec4ds const&, trivec4ds const&)>(reflect_on),
            sol::resolve<vec4ds(vec4ds const&, bivec4ds const&)>(reflect_on),
            sol::resolve<bivec4ds(bivec4ds const&, bivec4ds const&)>(reflect_on)));

    lua.set_function(
        "reflect_on_vec",
        sol::overload(
            sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on_vec),
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on_vec),
            // sta
            sol::resolve<vec4ds(vec4ds const&, vec4ds const&)>(reflect_on_vec)));

    ////////////////////////////////////////////////////////////////////////////////
    // testing for congruence
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "is_congruent",
        sol::overload(
            // EGA 2D congruence (with explicit tolerance)
            sol::resolve<bool(scalar2d, scalar2d, value_t)>(is_congruent),
            sol::resolve<bool(vec2d const&, vec2d const&, value_t)>(is_congruent),
            sol::resolve<bool(pscalar2d, pscalar2d, value_t)>(is_congruent),
            // EGA 2D congruence (with default tolerance)
            [](scalar2d a, scalar2d b) { return is_congruent(a, b); },
            [](vec2d const& a, vec2d const& b) { return is_congruent(a, b); },
            [](pscalar2d a, pscalar2d b) { return is_congruent(a, b); },
            // EGA 3D congruence (with explicit tolerance)
            sol::resolve<bool(scalar3d, scalar3d, value_t)>(is_congruent),
            sol::resolve<bool(vec3d const&, vec3d const&, value_t)>(is_congruent),
            sol::resolve<bool(bivec3d const&, bivec3d const&, value_t)>(is_congruent),
            sol::resolve<bool(pscalar3d, pscalar3d, value_t)>(is_congruent),
            // EGA 3D congruence (with default tolerance)
            [](scalar3d a, scalar3d b) { return is_congruent(a, b); },
            [](vec3d const& a, vec3d const& b) { return is_congruent(a, b); },
            [](bivec3d const& a, bivec3d const& b) { return is_congruent(a, b); },
            [](pscalar3d a, pscalar3d b) { return is_congruent(a, b); },
            // PGA 2DP congruence (with explicit tolerance)
            sol::resolve<bool(scalar2dp, scalar2dp, value_t)>(is_congruent),
            sol::resolve<bool(vec2dp const&, vec2dp const&, value_t)>(is_congruent),
            sol::resolve<bool(bivec2dp const&, bivec2dp const&, value_t)>(is_congruent),
            sol::resolve<bool(pscalar2dp, pscalar2dp, value_t)>(is_congruent),
            // PGA 2DP congruence (with default tolerance)
            [](scalar2dp a, scalar2dp b) { return is_congruent(a, b); },
            [](vec2dp const& a, vec2dp const& b) { return is_congruent(a, b); },
            [](bivec2dp const& a, bivec2dp const& b) { return is_congruent(a, b); },
            [](pscalar2dp a, pscalar2dp b) { return is_congruent(a, b); },
            // PGA 3DP congruence (with explicit tolerance)
            sol::resolve<bool(scalar3dp, scalar3dp, value_t)>(is_congruent),
            sol::resolve<bool(vec3dp const&, vec3dp const&, value_t)>(is_congruent),
            sol::resolve<bool(bivec3dp const&, bivec3dp const&, value_t)>(is_congruent),
            sol::resolve<bool(trivec3dp const&, trivec3dp const&, value_t)>(is_congruent),
            sol::resolve<bool(pscalar3dp, pscalar3dp, value_t)>(is_congruent),
            // PGA 3DP congruence (with default tolerance)
            [](scalar3dp a, scalar3dp b) { return is_congruent(a, b); },
            [](vec3dp const& a, vec3dp const& b) { return is_congruent(a, b); },
            [](bivec3dp const& a, bivec3dp const& b) { return is_congruent(a, b); },
            [](trivec3dp const& a, trivec3dp const& b) { return is_congruent(a, b); },
            [](pscalar3dp a, pscalar3dp b) { return is_congruent(a, b); },
            // sta
            sol::resolve<bool(scalar4ds, scalar4ds, value_t)>(is_congruent),
            sol::resolve<bool(vec4ds const&, vec4ds const&, value_t)>(is_congruent),
            sol::resolve<bool(bivec4ds const&, bivec4ds const&, value_t)>(is_congruent),
            sol::resolve<bool(trivec4ds const&, trivec4ds const&, value_t)>(is_congruent),
            sol::resolve<bool(pscalar4ds, pscalar4ds, value_t)>(is_congruent)));

    ////////////////////////////////////////////////////////////////////////////////
    // common helper functions for scripting in lua
    ////////////////////////////////////////////////////////////////////////////////

    // convert scalars & pscalars into numeric values for further calculations
    lua.set_function("to_val", sol::overload( // EGA to_val operations
                                   sol::resolve<value_t(scalar2d)>(to_val),
                                   sol::resolve<value_t(scalar3d)>(to_val),
                                   sol::resolve<value_t(pscalar2d)>(to_val),
                                   sol::resolve<value_t(pscalar3d)>(to_val),
                                   // PGA to_val operations
                                   sol::resolve<value_t(scalar2dp)>(to_val),
                                   sol::resolve<value_t(scalar3dp)>(to_val),
                                   sol::resolve<value_t(pscalar2dp)>(to_val),
                                   sol::resolve<value_t(pscalar3dp)>(to_val)));

    lua.set_function("rad2deg", &rad2deg);
    lua.set_function("deg2rad", &deg2rad);

    ////////////////////////////////////////////////////////////////////////////////
    // Critical missing functions - Only verified existing functions
    ////////////////////////////////////////////////////////////////////////////////

    // inv() - Multiplicative inverse (only existing overloads verified)
    lua.set_function("inv", sol::overload(
                                // EGA 2D types
                                sol::resolve<scalar2d(scalar2d)>(inv),
                                sol::resolve<vec2d(vec2d const&)>(inv),
                                sol::resolve<pscalar2d(pscalar2d)>(inv),
                                sol::resolve<mvec2d_e(mvec2d_e const&)>(inv),
                                sol::resolve<mvec2d(mvec2d const&)>(inv),
                                // EGA 3D types
                                sol::resolve<scalar3d(scalar3d)>(inv),
                                sol::resolve<vec3d(vec3d const&)>(inv),
                                sol::resolve<bivec3d(bivec3d const&)>(inv),
                                sol::resolve<pscalar3d(pscalar3d)>(inv),
                                sol::resolve<mvec3d_e(mvec3d_e const&)>(inv),
                                sol::resolve<mvec3d_u(mvec3d_u const&)>(inv),
                                sol::resolve<mvec3d(mvec3d const&)>(inv),
                                // PGA 2DP types (NO inv(pscalar2dp))
                                sol::resolve<scalar2dp(scalar2dp)>(inv),
                                sol::resolve<vec2dp(vec2dp const&)>(inv),
                                sol::resolve<bivec2dp(bivec2dp const&)>(inv),
                                sol::resolve<mvec2dp_e(mvec2dp_e const&)>(inv),
                                sol::resolve<mvec2dp_u(mvec2dp_u const&)>(inv),
                                sol::resolve<mvec2dp(mvec2dp const&)>(inv),
                                // PGA 3DP types (NO inv(pscalar3dp))
                                sol::resolve<scalar3dp(scalar3dp)>(inv),
                                sol::resolve<vec3dp(vec3dp const&)>(inv),
                                sol::resolve<bivec3dp(bivec3dp const&)>(inv),
                                sol::resolve<trivec3dp(trivec3dp const&)>(inv),
                                sol::resolve<mvec3dp_e(mvec3dp_e const&)>(inv),
                                sol::resolve<mvec3dp_u(mvec3dp_u const&)>(inv),
                                sol::resolve<mvec3dp(mvec3dp const&)>(inv),
                                // sta
                                sol::resolve<scalar4ds(scalar4ds)>(inv),
                                sol::resolve<vec4ds(vec4ds const&)>(inv),
                                sol::resolve<bivec4ds(bivec4ds const&)>(inv),
                                sol::resolve<trivec4ds(trivec4ds const&)>(inv),
                                sol::resolve<pscalar4ds(pscalar4ds)>(inv),
                                sol::resolve<mvec4ds_e(mvec4ds_e const&)>(inv),
                                sol::resolve<mvec4ds_u(mvec4ds_u const&)>(inv),
                                sol::resolve<mvec4ds(mvec4ds const&)>(inv)));


    // rrev() - Regressive reversion (verified existing signatures - PGA only)
    lua.set_function("rrev", sol::overload(
                                 // PGA 2DP types (verified signatures)
                                 sol::resolve<scalar2dp(scalar2dp)>(rrev),
                                 sol::resolve<vec2dp(vec2dp const&)>(rrev),
                                 sol::resolve<bivec2dp(bivec2dp const&)>(rrev),
                                 sol::resolve<pscalar2dp(pscalar2dp)>(rrev),
                                 sol::resolve<mvec2dp_e(mvec2dp_e const&)>(rrev),
                                 sol::resolve<mvec2dp_u(mvec2dp_u const&)>(rrev),
                                 sol::resolve<mvec2dp(mvec2dp const&)>(rrev),
                                 // PGA 3DP types (verified signatures)
                                 sol::resolve<scalar3dp(scalar3dp)>(rrev),
                                 sol::resolve<vec3dp(vec3dp const&)>(rrev),
                                 sol::resolve<bivec3dp(bivec3dp const&)>(rrev),
                                 sol::resolve<trivec3dp(trivec3dp const&)>(rrev),
                                 sol::resolve<pscalar3dp(pscalar3dp)>(rrev),
                                 sol::resolve<mvec3dp_e(mvec3dp_e const&)>(rrev),
                                 sol::resolve<mvec3dp_u(mvec3dp_u const&)>(rrev),
                                 sol::resolve<mvec3dp(mvec3dp const&)>(rrev),
                                 // sta
                                 sol::resolve<scalar4ds(scalar4ds)>(rrev),
                                 sol::resolve<vec4ds(vec4ds const&)>(rrev),
                                 sol::resolve<bivec4ds(bivec4ds const&)>(rrev),
                                 sol::resolve<trivec4ds(trivec4ds const&)>(rrev),
                                 sol::resolve<pscalar4ds(pscalar4ds)>(rrev),
                                 sol::resolve<mvec4ds_e(mvec4ds_e const&)>(rrev),
                                 sol::resolve<mvec4ds_u(mvec4ds_u const&)>(rrev),
                                 sol::resolve<mvec4ds(mvec4ds const&)>(rrev)));

    lua.set_function(
        "rdot", sol::overload(
                    // PGA 2DP types (verified signatures)
                    sol::resolve<pscalar2dp(mvec2dp const&, mvec2dp const&)>(rdot),
                    sol::resolve<pscalar2dp(mvec2dp_e const&, mvec2dp_e const&)>(rdot),
                    sol::resolve<pscalar2dp(mvec2dp_u const&, mvec2dp_u const&)>(rdot),
                    sol::resolve<pscalar2dp(pscalar2dp, pscalar2dp)>(rdot),
                    sol::resolve<pscalar2dp(bivec2dp const&, bivec2dp const&)>(rdot),
                    sol::resolve<pscalar2dp(vec2dp const&, vec2dp const&)>(rdot),
                    // PGA 3DP types (verified signatures)
                    sol::resolve<pscalar3dp(mvec3dp const&, mvec3dp const&)>(rdot),
                    sol::resolve<pscalar3dp(mvec3dp_e const&, mvec3dp_e const&)>(rdot),
                    sol::resolve<pscalar3dp(mvec3dp_u const&, mvec3dp_u const&)>(rdot),
                    sol::resolve<pscalar3dp(pscalar3dp, pscalar3dp)>(rdot),
                    sol::resolve<pscalar3dp(trivec3dp const&, trivec3dp const&)>(rdot),
                    sol::resolve<pscalar3dp(bivec3dp const&, bivec3dp const&)>(rdot),
                    sol::resolve<pscalar3dp(vec3dp const&, vec3dp const&)>(rdot)));

    lua.set_function("rcmt",
                     sol::overload(
                         // PGA 2DP types (verified signatures)
                         sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(rcmt),
                         sol::resolve<vec2dp(bivec2dp const&, bivec2dp const&)>(rcmt),
                         sol::resolve<bivec2dp(bivec2dp const&, vec2dp const&)>(rcmt),
                         sol::resolve<bivec2dp(vec2dp const&, bivec2dp const&)>(rcmt),
                         sol::resolve<vec2dp(vec2dp const&, vec2dp const&)>(rcmt),
                         // PGA 3DP types (verified signatures)
                         sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(rcmt),
                         sol::resolve<bivec3dp(trivec3dp const&, trivec3dp const&)>(rcmt),
                         sol::resolve<trivec3dp(trivec3dp const&, bivec3dp const&)>(rcmt),
                         sol::resolve<trivec3dp(bivec3dp const&, trivec3dp const&)>(rcmt),
                         sol::resolve<bivec3dp(bivec3dp const&, bivec3dp const&)>(rcmt),
                         sol::resolve<vec3dp(bivec3dp const&, vec3dp const&)>(rcmt),
                         sol::resolve<vec3dp(vec3dp const&, bivec3dp const&)>(rcmt),
                         sol::resolve<bivec3dp(vec3dp const&, vec3dp const&)>(rcmt)));

    lua.set_function(
        "rgpr", sol::overload(
                    // PGA 2DP types (verified signatures - starting with basic ones)
                    sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp(mvec2dp const&, mvec2dp_u const&)>(rgpr),
                    sol::resolve<mvec2dp(mvec2dp_u const&, mvec2dp const&)>(rgpr),
                    //
                    sol::resolve<mvec2dp_u(mvec2dp_e const&, mvec2dp_e const&)>(rgpr),
                    sol::resolve<mvec2dp_u(mvec2dp_u const&, mvec2dp_u const&)>(rgpr),
                    sol::resolve<mvec2dp_e(mvec2dp_e const&, mvec2dp_u const&)>(rgpr),
                    sol::resolve<mvec2dp_e(mvec2dp_u const&, mvec2dp_e const&)>(rgpr),
                    sol::resolve<mvec2dp_u(mvec2dp_e const&, bivec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp_u(bivec2dp const&, mvec2dp_e const&)>(rgpr),
                    sol::resolve<mvec2dp_u(mvec2dp_e const&, bivec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp_e(mvec2dp_u const&, bivec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp_e(bivec2dp const&, mvec2dp_u const&)>(rgpr),
                    sol::resolve<mvec2dp_e(mvec2dp_e const&, vec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp_e(vec2dp const&, mvec2dp_e const&)>(rgpr),
                    sol::resolve<mvec2dp_u(mvec2dp_u const&, vec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp_u(vec2dp const&, mvec2dp_u const&)>(rgpr),
                    sol::resolve<vec2dp(mvec2dp_e const&, scalar2dp)>(rgpr),
                    sol::resolve<vec2dp(scalar2dp, mvec2dp_e const&)>(rgpr),
                    sol::resolve<mvec2dp_e(mvec2dp_u const&, scalar2dp)>(rgpr),
                    sol::resolve<mvec2dp_e(scalar2dp, mvec2dp_u const&)>(rgpr),
                    //
                    sol::resolve<pscalar2dp(pscalar2dp, pscalar2dp)>(rgpr),
                    sol::resolve<bivec2dp(pscalar2dp, bivec2dp const&)>(rgpr),
                    sol::resolve<bivec2dp(bivec2dp const&, pscalar2dp)>(rgpr),
                    sol::resolve<vec2dp(pscalar2dp, vec2dp const&)>(rgpr),
                    sol::resolve<vec2dp(vec2dp const&, pscalar2dp)>(rgpr),
                    sol::resolve<scalar2dp(pscalar2dp, scalar2dp)>(rgpr),
                    sol::resolve<scalar2dp(scalar2dp, pscalar2dp)>(rgpr),
                    //
                    sol::resolve<mvec2dp_u(bivec2dp const&, bivec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp_e(bivec2dp const&, vec2dp const&)>(rgpr),
                    sol::resolve<mvec2dp_e(vec2dp const&, bivec2dp const&)>(rgpr),
                    sol::resolve<vec2dp(bivec2dp const&, scalar2dp)>(rgpr),
                    sol::resolve<vec2dp(scalar2dp, bivec2dp const&)>(rgpr),
                    //
                    sol::resolve<mvec2dp_u(vec2dp const&, vec2dp const&)>(rgpr),
                    sol::resolve<bivec2dp(vec2dp const&, scalar2dp)>(rgpr),
                    sol::resolve<bivec2dp(scalar2dp, vec2dp const&)>(rgpr),
                    //
                    sol::resolve<pscalar2dp(scalar2dp, scalar2dp)>(rgpr),
                    // PGA 3DP types (verified signatures - starting with basic ones)
                    sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp(mvec3dp_e const&, mvec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp(mvec3dp const&, mvec3dp_e const&)>(rgpr),
                    //
                    sol::resolve<mvec3dp_e(mvec3dp_e const&, mvec3dp_e const&)>(rgpr),
                    sol::resolve<mvec3dp_u(mvec3dp_u const&, mvec3dp_e const&)>(rgpr),
                    sol::resolve<mvec3dp_u(mvec3dp_e const&, trivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_u(trivec3dp const&, mvec3dp_e const&)>(rgpr),
                    sol::resolve<mvec3dp_e(mvec3dp_u const&, trivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_e(trivec3dp const&, mvec3dp_u const&)>(rgpr),
                    sol::resolve<mvec3dp_e(mvec3dp_e const&, bivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_e(bivec3dp const&, mvec3dp_e const&)>(rgpr),
                    sol::resolve<mvec3dp_u(mvec3dp_u const&, bivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_u(bivec3dp const&, mvec3dp_u const&)>(rgpr),
                    sol::resolve<mvec3dp_u(mvec3dp_e const&, vec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_u(vec3dp const&, mvec3dp_e const&)>(rgpr),
                    sol::resolve<mvec3dp_e(mvec3dp_e const&, scalar3dp)>(rgpr),
                    sol::resolve<mvec3dp_e(scalar3dp, mvec3dp_e const&)>(rgpr),
                    sol::resolve<mvec3dp_u(mvec3dp_u const&, scalar3dp)>(rgpr),
                    sol::resolve<mvec3dp_u(scalar3dp, mvec3dp_u const&)>(rgpr),
                    //
                    sol::resolve<pscalar3dp(pscalar3dp, pscalar3dp)>(rgpr),
                    sol::resolve<trivec3dp(pscalar3dp, trivec3dp const&)>(rgpr),
                    sol::resolve<trivec3dp(trivec3dp const&, pscalar3dp)>(rgpr),
                    sol::resolve<bivec3dp(pscalar3dp, bivec3dp const&)>(rgpr),
                    sol::resolve<bivec3dp(bivec3dp const&, pscalar3dp)>(rgpr),
                    sol::resolve<vec3dp(pscalar3dp, vec3dp const&)>(rgpr),
                    sol::resolve<vec3dp(vec3dp const&, pscalar3dp)>(rgpr),
                    sol::resolve<scalar3dp(pscalar3dp, scalar3dp)>(rgpr),
                    sol::resolve<scalar3dp(scalar3dp, pscalar3dp)>(rgpr),
                    //
                    sol::resolve<mvec3dp_e(trivec3dp const&, trivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_u(trivec3dp const&, bivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_u(bivec3dp const&, trivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_e(trivec3dp const&, vec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_e(vec3dp const&, trivec3dp const&)>(rgpr),
                    sol::resolve<vec3dp(trivec3dp const&, scalar3dp)>(rgpr),
                    sol::resolve<vec3dp(scalar3dp, trivec3dp const&)>(rgpr),
                    //
                    sol::resolve<mvec3dp_e(bivec3dp const&, bivec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_u(bivec3dp const&, vec3dp const&)>(rgpr),
                    sol::resolve<mvec3dp_u(vec3dp const&, bivec3dp const&)>(rgpr),
                    sol::resolve<bivec3dp(bivec3dp const&, scalar3dp)>(rgpr),
                    sol::resolve<bivec3dp(scalar3dp, bivec3dp const&)>(rgpr),
                    //
                    sol::resolve<mvec3dp_e(vec3dp const&, vec3dp const&)>(rgpr),
                    sol::resolve<trivec3dp(vec3dp const&, scalar3dp)>(rgpr),
                    sol::resolve<trivec3dp(scalar3dp, vec3dp const&)>(rgpr),
                    //
                    sol::resolve<pscalar3dp(scalar3dp, scalar3dp)>(rgpr)));

    // ---- Phase 1: EGA function parity with ga_py (generated via
    // utilities/gen_lua_overloads.py) ----
    lua.set_function(
        "gs_orthogonal",
        sol::overload(
            sol::resolve<std::vector<vec2d>(vec2d const&, vec2d const&)>(gs_orthogonal),
            sol::resolve<std::vector<vec3d>(vec3d const&, vec3d const&)>(gs_orthogonal),
            sol::resolve<std::vector<vec3d>(vec3d const&, vec3d const&, vec3d const&)>(
                gs_orthogonal)));

    lua.set_function(
        "gs_orthonormal",
        sol::overload(
            sol::resolve<std::vector<vec2d>(vec2d const&, vec2d const&)>(gs_orthonormal),
            sol::resolve<std::vector<vec3d>(vec3d const&, vec3d const&)>(gs_orthonormal),
            sol::resolve<std::vector<vec3d>(vec3d const&, vec3d const&, vec3d const&)>(
                gs_orthonormal)));

    lua.set_function(
        "l_expand",
        sol::overload(sol::resolve<bivec3d(bivec3d const&, vec3d const&)>(l_expand)));

    lua.set_function(
        "r_expand",
        sol::overload(sol::resolve<bivec3d(vec3d const&, bivec3d const&)>(r_expand)));

    lua.set_function(
        "rotate_opt",
        sol::overload(
            sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(rotate_opt),
            sol::resolve<std::vector<vec2d>(std::vector<vec2d> const&, mvec2d_e const&)>(
                rotate_opt),
            sol::resolve<vec3d(vec3d const&, mvec3d_e const&)>(rotate_opt),
            sol::resolve<std::vector<vec3d>(std::vector<vec3d> const&, mvec3d_e const&)>(
                rotate_opt),
            sol::resolve<bivec3d(bivec3d const&, mvec3d_e const&)>(rotate_opt),
            sol::resolve<std::vector<bivec3d>(std::vector<bivec3d> const&,
                                              mvec3d_e const&)>(rotate_opt)));

    lua.set_function("log", sol::overload(sol::resolve<pscalar2d(mvec2d_e const&)>(log),
                                          sol::resolve<bivec3d(mvec3d_e const&)>(log),
                                          // sta
                                          sol::resolve<bivec4ds(mvec4ds_e const&)>(log)));

    lua.set_function("sqrt",
                     sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&)>(sqrt),
                                   sol::resolve<mvec3d_e(mvec3d_e const&)>(sqrt),
                                   // sta
                                   sol::resolve<mvec4ds_e(mvec4ds_e const&)>(sqrt)));

    // pga versor exponential/logarithm/sqrt w.r.t. the regressive geometric product
    lua.set_function("rexp",
                     sol::overload(sol::resolve<mvec2dp_u(vec2dp const&)>(rexp),
                                   sol::resolve<mvec3dp_e(bivec3dp const&)>(rexp)));
    lua.set_function("rlog",
                     sol::overload(sol::resolve<vec2dp(mvec2dp_u const&)>(rlog),
                                   sol::resolve<bivec3dp(mvec3dp_e const&)>(rlog)));
    lua.set_function("rsqrt",
                     sol::overload(sol::resolve<mvec2dp_u(mvec2dp_u const&)>(rsqrt),
                                   sol::resolve<mvec3dp_e(mvec3dp_e const&)>(rsqrt)));

    lua.set_function(
        "twdg1",
        sol::overload(sol::resolve<vec2d(pscalar2d, vec2d const&)>(twdg1),
                      sol::resolve<vec2d(vec2d const&, pscalar2d)>(twdg1),
                      sol::resolve<scalar2d(vec2d const&, vec2d const&)>(twdg1),
                      sol::resolve<bivec3d(pscalar3d, vec3d const&)>(twdg1),
                      sol::resolve<bivec3d(vec3d const&, pscalar3d)>(twdg1),
                      sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(twdg1),
                      sol::resolve<vec3d(bivec3d const&, vec3d const&)>(twdg1),
                      sol::resolve<vec3d(vec3d const&, bivec3d const&)>(twdg1),
                      sol::resolve<scalar3d(vec3d const&, vec3d const&)>(twdg1),
                      sol::resolve<bivec2dp(pscalar2dp, vec2dp const&)>(twdg1),
                      sol::resolve<bivec2dp(vec2dp const&, pscalar2dp)>(twdg1),
                      sol::resolve<bivec2dp(bivec2dp const&, bivec2dp const&)>(twdg1),
                      sol::resolve<vec2dp(bivec2dp const&, vec2dp const&)>(twdg1),
                      sol::resolve<vec2dp(vec2dp const&, bivec2dp const&)>(twdg1),
                      sol::resolve<scalar2dp(vec2dp const&, vec2dp const&)>(twdg1),
                      sol::resolve<trivec3dp(pscalar3dp, vec3dp const&)>(twdg1),
                      sol::resolve<trivec3dp(vec3dp const&, pscalar3dp)>(twdg1),
                      sol::resolve<bivec3dp(trivec3dp const&, vec3dp const&)>(twdg1),
                      sol::resolve<bivec3dp(vec3dp const&, trivec3dp const&)>(twdg1),
                      sol::resolve<trivec3dp(trivec3dp const&, bivec3dp const&)>(twdg1),
                      sol::resolve<trivec3dp(bivec3dp const&, trivec3dp const&)>(twdg1),
                      sol::resolve<bivec3dp(bivec3dp const&, bivec3dp const&)>(twdg1),
                      sol::resolve<vec3dp(bivec3dp const&, vec3dp const&)>(twdg1),
                      sol::resolve<vec3dp(vec3dp const&, bivec3dp const&)>(twdg1),
                      sol::resolve<scalar3dp(vec3dp const&, vec3dp const&)>(twdg1)));

    lua.set_function(
        "rtwdg1",
        sol::overload(
            sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(rtwdg1),
            sol::resolve<vec2d(vec2d const&, scalar2d)>(rtwdg1),
            sol::resolve<vec2d(scalar2d, vec2d const&)>(rtwdg1),
            sol::resolve<pscalar3d(bivec3d const&, bivec3d const&)>(rtwdg1),
            sol::resolve<bivec3d(bivec3d const&, vec3d const&)>(rtwdg1),
            sol::resolve<bivec3d(vec3d const&, bivec3d const&)>(rtwdg1),
            sol::resolve<vec3d(bivec3d const&, scalar3d)>(rtwdg1),
            sol::resolve<vec3d(scalar3d, bivec3d const&)>(rtwdg1),
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(rtwdg1),
            sol::resolve<pscalar2dp(bivec2dp const&, bivec2dp const&)>(rtwdg1),
            sol::resolve<bivec2dp(bivec2dp const&, vec2dp const&)>(rtwdg1),
            sol::resolve<bivec2dp(vec2dp const&, bivec2dp const&)>(rtwdg1),
            sol::resolve<vec2dp(bivec2dp const&, scalar2dp)>(rtwdg1),
            sol::resolve<vec2dp(scalar2dp, bivec2dp const&)>(rtwdg1),
            sol::resolve<vec2dp(vec2dp const&, vec2dp const&)>(rtwdg1),
            sol::resolve<pscalar3dp(trivec3dp const&, trivec3dp const&)>(rtwdg1),
            sol::resolve<trivec3dp(trivec3dp const&, bivec3dp const&)>(rtwdg1),
            sol::resolve<trivec3dp(bivec3dp const&, trivec3dp const&)>(rtwdg1),
            sol::resolve<bivec3dp(trivec3dp const&, vec3dp const&)>(rtwdg1),
            sol::resolve<bivec3dp(vec3dp const&, trivec3dp const&)>(rtwdg1),
            sol::resolve<vec3dp(trivec3dp const&, scalar3dp)>(rtwdg1),
            sol::resolve<vec3dp(scalar3dp, trivec3dp const&)>(rtwdg1),
            sol::resolve<bivec3dp(bivec3dp const&, bivec3dp const&)>(rtwdg1),
            sol::resolve<vec3dp(bivec3dp const&, vec3dp const&)>(rtwdg1),
            sol::resolve<vec3dp(vec3dp const&, bivec3dp const&)>(rtwdg1)));

    // ---- Phase 2a: PGA function parity with ga_py (generated via
    // utilities/gen_lua_overloads.py) ----

    // gr(): grade of a graded object (top-level hd::ga; also used by the dist
    // forwarders in register_forwarders())
    lua.set_function(
        "gr",
        sol::overload(
            // EGA
            sol::resolve<size_t(scalar2d)>(gr), sol::resolve<size_t(vec2d const&)>(gr),
            sol::resolve<size_t(pscalar2d)>(gr), sol::resolve<size_t(scalar3d)>(gr),
            sol::resolve<size_t(vec3d const&)>(gr),
            sol::resolve<size_t(bivec3d const&)>(gr), sol::resolve<size_t(pscalar3d)>(gr),
            // PGA
            sol::resolve<size_t(scalar2dp)>(gr), sol::resolve<size_t(vec2dp const&)>(gr),
            sol::resolve<size_t(bivec2dp const&)>(gr),
            sol::resolve<size_t(pscalar2dp)>(gr), sol::resolve<size_t(scalar3dp)>(gr),
            sol::resolve<size_t(vec3dp const&)>(gr),
            sol::resolve<size_t(bivec3dp const&)>(gr),
            sol::resolve<size_t(trivec3dp const&)>(gr),
            sol::resolve<size_t(pscalar3dp)>(gr)));

    lua.set_function("rinv",
                     sol::overload(sol::resolve<vec2dp(vec2dp const&)>(rinv),
                                   sol::resolve<bivec2dp(bivec2dp const&)>(rinv),
                                   sol::resolve<pscalar2dp(pscalar2dp)>(rinv),
                                   sol::resolve<mvec2dp_e(mvec2dp_e const&)>(rinv),
                                   sol::resolve<mvec2dp_u(mvec2dp_u const&)>(rinv),
                                   sol::resolve<mvec2dp(mvec2dp const&)>(rinv),
                                   sol::resolve<vec3dp(vec3dp const&)>(rinv),
                                   sol::resolve<bivec3dp(bivec3dp const&)>(rinv),
                                   sol::resolve<trivec3dp(trivec3dp const&)>(rinv),
                                   sol::resolve<pscalar3dp(pscalar3dp)>(rinv),
                                   sol::resolve<mvec3dp_e(mvec3dp_e const&)>(rinv),
                                   sol::resolve<mvec3dp_u(mvec3dp_u const&)>(rinv),
                                   sol::resolve<mvec3dp(mvec3dp const&)>(rinv)));

    lua.set_function(
        "move2dp_opt",
        sol::overload(
            sol::resolve<vec2dp(vec2dp const&, mvec2dp_u const&)>(move2dp_opt),
            sol::resolve<bivec2dp(bivec2dp const&, mvec2dp_u const&)>(move2dp_opt)));

    lua.set_function(
        "move3dp_opt",
        sol::overload(
            sol::resolve<vec3dp(vec3dp const&, mvec3dp_e const&)>(move3dp_opt),
            sol::resolve<bivec3dp(bivec3dp const&, mvec3dp_e const&)>(move3dp_opt),
            sol::resolve<trivec3dp(trivec3dp const&, mvec3dp_e const&)>(move3dp_opt)));

    lua.set_function(
        "invert_on",
        sol::overload(
            sol::resolve<vec2dp(vec2dp const&, vec2dp const&)>(invert_on),
            sol::resolve<bivec2dp(bivec2dp const&, vec2dp const&)>(invert_on),
            sol::resolve<vec3dp(vec3dp const&, vec3dp const&)>(invert_on),
            sol::resolve<bivec3dp(bivec3dp const&, vec3dp const&)>(invert_on),
            sol::resolve<trivec3dp(trivec3dp const&, vec3dp const&)>(invert_on)));

    // expand(): join of a point with a line/plane into the higher-grade flat
    // (operates on the convenience types registered in register_convenience_types)
    lua.set_function(
        "expand",
        sol::overload(sol::resolve<line2d(point2d const&, line2d const&)>(expand),
                      sol::resolve<line3d(point3d const&, plane3d const&)>(expand),
                      sol::resolve<plane3d(point3d const&, line3d const&)>(expand),
                      sol::resolve<plane3d(line3d const&, plane3d const&)>(expand)));

    // pose <-> motor conversions (operate on the pose PODs)
    lua.set_function("motor_from_pose2dp",
                     sol::resolve<mvec2dp_u(pose2dp const&)>(motor_from_pose2dp));
    lua.set_function("motor_from_pose3dp",
                     sol::resolve<mvec3dp_e(pose3dp const&)>(motor_from_pose3dp));
    lua.set_function("pose2dp_from_motor",
                     sol::resolve<pose2dp(mvec2dp_u const&)>(pose2dp_from_motor));
    lua.set_function("pose3dp_from_motor",
                     sol::resolve<pose3dp(mvec3dp_e const&)>(pose3dp_from_motor));

    // ---- Phase 3c: STA-only functions (generated via utilities/gen_lua_overloads.py)
    // ----
    lua.set_function(
        "get_boost",
        sol::overload(sol::resolve<mvec4ds_e(bivec4ds const&, value_t)>(get_boost)));

    lua.set_function("is_lightlike",
                     sol::overload(sol::resolve<bool(vec4ds const&)>(is_lightlike),
                                   sol::resolve<bool(bivec4ds const&)>(is_lightlike),
                                   sol::resolve<bool(trivec4ds const&)>(is_lightlike)));

    lua.set_function("is_spacelike",
                     sol::overload(sol::resolve<bool(vec4ds const&)>(is_spacelike),
                                   sol::resolve<bool(bivec4ds const&)>(is_spacelike),
                                   sol::resolve<bool(trivec4ds const&)>(is_spacelike)));

    lua.set_function("is_timelike",
                     sol::overload(sol::resolve<bool(vec4ds const&)>(is_timelike),
                                   sol::resolve<bool(bivec4ds const&)>(is_timelike),
                                   sol::resolve<bool(trivec4ds const&)>(is_timelike)));

    lua.set_function("l_undual",
                     sol::overload(sol::resolve<pscalar4ds(scalar4ds)>(l_undual),
                                   sol::resolve<trivec4ds(vec4ds const&)>(l_undual),
                                   sol::resolve<bivec4ds(bivec4ds const&)>(l_undual),
                                   sol::resolve<vec4ds(trivec4ds const&)>(l_undual),
                                   sol::resolve<scalar4ds(pscalar4ds)>(l_undual),
                                   sol::resolve<mvec4ds_e(mvec4ds_e const&)>(l_undual),
                                   sol::resolve<mvec4ds_u(mvec4ds_u const&)>(l_undual),
                                   sol::resolve<mvec4ds(mvec4ds const&)>(l_undual)));

    lua.set_function("r_undual",
                     sol::overload(sol::resolve<pscalar4ds(scalar4ds)>(r_undual),
                                   sol::resolve<trivec4ds(vec4ds const&)>(r_undual),
                                   sol::resolve<bivec4ds(bivec4ds const&)>(r_undual),
                                   sol::resolve<vec4ds(trivec4ds const&)>(r_undual),
                                   sol::resolve<scalar4ds(pscalar4ds)>(r_undual),
                                   sol::resolve<mvec4ds_e(mvec4ds_e const&)>(r_undual),
                                   sol::resolve<mvec4ds_u(mvec4ds_u const&)>(r_undual),
                                   sol::resolve<mvec4ds(mvec4ds const&)>(r_undual)));

    lua.set_function(
        "rapidity",
        sol::overload(sol::resolve<value_t(vec4ds const&, vec4ds const&)>(rapidity)));

    lua.set_function("rel_bivec_split",
                     sol::overload(sol::resolve<bivec4ds(bivec4ds const&, vec4ds const&)>(
                         rel_bivec_split)));

    lua.set_function("rel_vec_split",
                     sol::overload(sol::resolve<bivec4ds(bivec4ds const&, vec4ds const&)>(
                         rel_vec_split)));

    lua.set_function(
        "space_split",
        sol::overload(sol::resolve<bivec4ds(vec4ds const&, vec4ds const&)>(space_split)));

    lua.set_function(
        "time_split",
        sol::overload(sol::resolve<scalar4ds(vec4ds const&, vec4ds const&)>(time_split)));

    lua.set_function(
        "transform",
        sol::overload(
            sol::resolve<vec4ds(vec4ds const&, mvec4ds_e const&)>(transform),
            sol::resolve<bivec4ds(bivec4ds const&, mvec4ds_e const&)>(transform),
            sol::resolve<trivec4ds(trivec4ds const&, mvec4ds_e const&)>(transform),
            sol::resolve<mvec4ds(mvec4ds const&, mvec4ds_e const&)>(transform)));

    lua.set_function(
        "transform_opt",
        sol::overload(
            sol::resolve<vec4ds(vec4ds const&, mvec4ds_e const&)>(transform_opt),
            sol::resolve<trivec4ds(trivec4ds const&, mvec4ds_e const&)>(transform_opt),
            sol::resolve<bivec4ds(bivec4ds const&, mvec4ds_e const&)>(transform_opt),
            sol::resolve<std::vector<vec4ds>(std::vector<vec4ds> const&,
                                             mvec4ds_e const&)>(transform_opt),
            sol::resolve<std::vector<trivec4ds>(std::vector<trivec4ds> const&,
                                                mvec4ds_e const&)>(transform_opt),
            sol::resolve<std::vector<bivec4ds>(std::vector<bivec4ds> const&,
                                               mvec4ds_e const&)>(transform_opt)));

    // ---- Phase 4: top-level helpers (unit conversions, rgr, sign). The
    // solver/integrator free functions lu_solve/lstsq_solve/rk4_* are skipped
    // (matrix/vector-span args don't map well to a Lua REPL). ----
    lua.set_function("Hz2radps", sol::overload(sol::resolve<value_t(value_t)>(Hz2radps)));

    lua.set_function("radps2Hz", sol::overload(sol::resolve<value_t(value_t)>(radps2Hz)));

    lua.set_function("radps2rpm",
                     sol::overload(sol::resolve<value_t(value_t)>(radps2rpm)));

    lua.set_function("rpm2radps",
                     sol::overload(sol::resolve<value_t(value_t)>(rpm2radps)));

    lua.set_function("rgr", sol::overload(sol::resolve<size_t(scalar2dp)>(rgr),
                                          sol::resolve<size_t(vec2dp const&)>(rgr),
                                          sol::resolve<size_t(bivec2dp const&)>(rgr),
                                          sol::resolve<size_t(pscalar2dp)>(rgr),
                                          sol::resolve<size_t(scalar3dp)>(rgr),
                                          sol::resolve<size_t(vec3dp const&)>(rgr),
                                          sol::resolve<size_t(bivec3dp const&)>(rgr),
                                          sol::resolve<size_t(trivec3dp const&)>(rgr),
                                          sol::resolve<size_t(pscalar3dp)>(rgr),
                                          sol::resolve<size_t(scalar4ds)>(rgr),
                                          sol::resolve<size_t(vec4ds const&)>(rgr),
                                          sol::resolve<size_t(bivec4ds const&)>(rgr),
                                          sol::resolve<size_t(trivec4ds const&)>(rgr),
                                          sol::resolve<size_t(pscalar4ds)>(rgr)));

    lua.set_function("sign", sol::resolve<value_t(value_t)>(sign));

    // Note:
    // - Geometric product is only available as operator*, not as
    // gpr() function
    // - Left/right contractions are only available as operator<<
    // and operator>>
    // - Many operations have both function and operator forms
    // available
}

////////////////////////////////////////////////////////////////////////////////
// make defined constants available as global variables in lua
////////////////////////////////////////////////////////////////////////////////
void register_constants(sol::state& lua)
{
    // TODO: manipulate global lua metatable such that these entries
    //       cannot be modified by the user

    using namespace hd::ga;
    using namespace hd::ga::ega;

    // general constants
    lua["eps"] = eps; // this is from ga_value_t.hpp,
                      // while all else is from ga_usr_consts.hpp

    // 2d constants
    lua["e1_2d"] = e1_2d; // as 2d vector
    lua["e2_2d"] = e2_2d;
    lua["x_dir_2d"] = x_dir_2d;
    lua["y_dir_2d"] = y_dir_2d;
    lua["O_2d"] = O_2d;

    lua["e1_2d_mv"] = e1_2d_mv; // as 2d multivector
    lua["e2_2d_mv"] = e2_2d_mv;

    lua["I_2d"] = I_2d; // as pscalar2d
    lua["e12_2d"] = I_2d;
    lua["I_2d_mv"] = I_2d_mv;     // as 2d multivector
    lua["I_2d_mv_e"] = I_2d_mv_e; // as even-grade 2d multivector

    lua["one_2d"] = one_2d;           // scalar 1
    lua["one_2d_mv"] = one_2d_mv;     // 1 as 2d multivector
    lua["one_2d_mv_e"] = one_2d_mv_e; // 1 as even-grade 2d multivector

    // 3d constants
    lua["e1_3d"] = e1_3d; // as 3d vector
    lua["e2_3d"] = e2_3d;
    lua["e3_3d"] = e3_3d;
    lua["x_dir_3d"] = x_dir_3d;
    lua["y_dir_3d"] = y_dir_3d;
    lua["z_dir_3d"] = z_dir_3d;
    lua["O_3d"] = O_3d;

    lua["e1_3d_mv"] = e1_3d_mv; // as 3d multivector
    lua["e2_3d_mv"] = e2_3d_mv;
    lua["e3_3d_mv"] = e3_3d_mv;

    lua["e23_3d"] = e23_3d; // as 3d bivector
    lua["e31_3d"] = e31_3d;
    lua["e12_3d"] = e12_3d;
    lua["yz_3d"] = yz_3d;
    lua["zx_3d"] = zx_3d;
    lua["xy_3d"] = xy_3d;


    lua["e23_3d_mv"] = e23_3d_mv; // as 3d multivector
    lua["e31_3d_mv"] = e31_3d_mv;
    lua["e12_3d_mv"] = e12_3d_mv;
    lua["e23_3d_mv_e"] = e23_3d_mv_e; // as even-grade 3d multivector
    lua["e31_3d_mv_e"] = e31_3d_mv_e;
    lua["e12_3d_mv_e"] = e12_3d_mv_e;

    lua["I_3d"] = I_3d; // as pscalar3d
    lua["e123_3d"] = e123_3d;
    lua["I_3d_mv"] = I_3d_mv;     // as 3d multivector
    lua["I_3d_mv_u"] = I_3d_mv_u; // as odd-grade 3d multivector

    lua["one_3d"] = one_3d;           // scalar 1
    lua["one_3d_mv"] = one_3d_mv;     // 1 as 3d multivector
    lua["one_3d_mv_e"] = one_3d_mv_e; // 1 as even-grade 3d multivector

    ////////////////////////////////////////////////////////////////////////////////
    // PGA 2dp constants
    ////////////////////////////////////////////////////////////////////////////////
    using namespace hd::ga::pga;

    // 2dp basis vectors
    lua["e1_2dp"] = e1_2dp;
    lua["e2_2dp"] = e2_2dp;
    lua["e3_2dp"] = e3_2dp;
    lua["O_2dp"] = O_2dp;

    // 2dp direction vectors
    lua["x_dir_2dp"] = x_dir_2dp;
    lua["y_dir_2dp"] = y_dir_2dp;

    // 2dp multivector versions
    lua["e1_2dp_mv"] = e1_2dp_mv;
    lua["e2_2dp_mv"] = e2_2dp_mv;
    lua["e3_2dp_mv"] = e3_2dp_mv;

    // 2dp bivectors
    lua["e31_2dp"] = e31_2dp;
    lua["e32_2dp"] = e32_2dp;
    lua["e12_2dp"] = e12_2dp;

    // 2dp axes as bivectors
    lua["x_axis_2dp"] = x_axis_2dp;
    lua["y_axis_2dp"] = y_axis_2dp;
    lua["H_2dp"] = H_2dp;

    // 2dp multivector / even-grade forms of bivectors
    lua["e31_2dp_mv"] = e31_2dp_mv;
    lua["e31_2dp_mv_e"] = e31_2dp_mv_e;
    lua["e32_2dp_mv"] = e32_2dp_mv;
    lua["e32_2dp_mv_e"] = e32_2dp_mv_e;
    lua["e12_2dp_mv"] = e12_2dp_mv;
    lua["e12_2dp_mv_e"] = e12_2dp_mv_e;

    // 2dp pseudoscalar
    lua["I_2dp"] = I_2dp;
    lua["e321_2dp"] = e321_2dp;
    lua["I_2dp_mv"] = I_2dp_mv;
    lua["I_2dp_mv_u"] = I_2dp_mv_u;

    // 2dp scalar one
    lua["one_2dp"] = one_2dp;
    lua["one_2dp_mv"] = one_2dp_mv;
    lua["one_2dp_mv_e"] = one_2dp_mv_e;

    ////////////////////////////////////////////////////////////////////////////////
    // PGA 3dp constants
    ////////////////////////////////////////////////////////////////////////////////

    // 3dp basis vectors
    lua["e1_3dp"] = e1_3dp;
    lua["e2_3dp"] = e2_3dp;
    lua["e3_3dp"] = e3_3dp;
    lua["e4_3dp"] = e4_3dp;
    lua["O_3dp"] = O_3dp;

    // 3dp direction vectors
    lua["x_dir_3dp"] = x_dir_3dp;
    lua["y_dir_3dp"] = y_dir_3dp;
    lua["z_dir_3dp"] = z_dir_3dp;

    // 3dp multivector versions
    lua["e1_3dp_mv"] = e1_3dp_mv;
    lua["e2_3dp_mv"] = e2_3dp_mv;
    lua["e3_3dp_mv"] = e3_3dp_mv;
    lua["e4_3dp_mv"] = e4_3dp_mv;

    // 3dp bivectors
    lua["e41_3dp"] = e41_3dp;
    lua["e42_3dp"] = e42_3dp;
    lua["e43_3dp"] = e43_3dp;
    lua["e23_3dp"] = e23_3dp;
    lua["e31_3dp"] = e31_3dp;
    lua["e12_3dp"] = e12_3dp;

    // 3dp axes as bivectors
    lua["x_axis_3dp"] = x_axis_3dp;
    lua["y_axis_3dp"] = y_axis_3dp;
    lua["z_axis_3dp"] = z_axis_3dp;

    // 3dp trivectors
    lua["e423_3dp"] = e423_3dp;
    lua["e431_3dp"] = e431_3dp;
    lua["e412_3dp"] = e412_3dp;
    lua["e321_3dp"] = e321_3dp;

    // 3dp planes
    lua["yz_3dp"] = yz_3dp;
    lua["zx_3dp"] = zx_3dp;
    lua["xy_3dp"] = xy_3dp;
    lua["H_3dp"] = H_3dp;

    // 3dp pseudoscalar
    lua["I_3dp"] = I_3dp;
    lua["e1234_3dp"] = e1234_3dp;

    // 3dp multivector / even-grade forms of bivectors
    lua["e41_3dp_mv"] = e41_3dp_mv;
    lua["e41_3dp_mv_e"] = e41_3dp_mv_e;
    lua["e42_3dp_mv"] = e42_3dp_mv;
    lua["e42_3dp_mv_e"] = e42_3dp_mv_e;
    lua["e43_3dp_mv"] = e43_3dp_mv;
    lua["e43_3dp_mv_e"] = e43_3dp_mv_e;
    lua["e23_3dp_mv"] = e23_3dp_mv;
    lua["e23_3dp_mv_e"] = e23_3dp_mv_e;
    lua["e31_3dp_mv"] = e31_3dp_mv;
    lua["e31_3dp_mv_e"] = e31_3dp_mv_e;
    lua["e12_3dp_mv"] = e12_3dp_mv;
    lua["e12_3dp_mv_e"] = e12_3dp_mv_e;

    // 3dp multivector forms of trivectors
    lua["e423_3dp_mv"] = e423_3dp_mv;
    lua["e431_3dp_mv"] = e431_3dp_mv;
    lua["e412_3dp_mv"] = e412_3dp_mv;
    lua["e321_3dp_mv"] = e321_3dp_mv;

    // 3dp pseudoscalar multivector forms
    lua["I_3dp_mv"] = I_3dp_mv;
    lua["I_3dp_mv_e"] = I_3dp_mv_e;

    // 3dp scalar one
    lua["one_3dp"] = one_3dp;
    lua["one_3dp_mv"] = one_3dp_mv;
    lua["one_3dp_mv_e"] = one_3dp_mv_e;

    ////////////////////////////////////////////////////////////////////////////////
    // STA constants (G(1,3,0))
    ////////////////////////////////////////////////////////////////////////////////
    using namespace hd::ga::sta;

    // STA basis vectors (g1..g4 spacelike/timelike)
    lua["g1_4ds"] = g1_4ds;
    lua["g2_4ds"] = g2_4ds;
    lua["g3_4ds"] = g3_4ds;
    lua["g4_4ds"] = g4_4ds;

    // STA basis vectors as multivectors
    lua["g1_4ds_mv"] = g1_4ds_mv;
    lua["g2_4ds_mv"] = g2_4ds_mv;
    lua["g3_4ds_mv"] = g3_4ds_mv;
    lua["g4_4ds_mv"] = g4_4ds_mv;

    // STA direction vectors
    lua["x_dir_4ds"] = x_dir_4ds;
    lua["y_dir_4ds"] = y_dir_4ds;
    lua["z_dir_4ds"] = z_dir_4ds;
    lua["t_dir_4ds"] = t_dir_4ds;

    // STA bivectors
    lua["g14_4ds"] = g14_4ds;
    lua["g24_4ds"] = g24_4ds;
    lua["g34_4ds"] = g34_4ds;
    lua["g23_4ds"] = g23_4ds;
    lua["g31_4ds"] = g31_4ds;
    lua["g12_4ds"] = g12_4ds;

    // STA bivectors as multivectors
    lua["g14_4ds_mv"] = g14_4ds_mv;
    lua["g24_4ds_mv"] = g24_4ds_mv;
    lua["g34_4ds_mv"] = g34_4ds_mv;
    lua["g23_4ds_mv"] = g23_4ds_mv;
    lua["g31_4ds_mv"] = g31_4ds_mv;
    lua["g12_4ds_mv"] = g12_4ds_mv;

    // STA bivectors as even-grade multivectors
    lua["g14_4ds_mv_e"] = g14_4ds_mv_e;
    lua["g24_4ds_mv_e"] = g24_4ds_mv_e;
    lua["g34_4ds_mv_e"] = g34_4ds_mv_e;
    lua["g23_4ds_mv_e"] = g23_4ds_mv_e;
    lua["g31_4ds_mv_e"] = g31_4ds_mv_e;
    lua["g12_4ds_mv_e"] = g12_4ds_mv_e;

    // STA trivectors
    lua["g234_4ds"] = g234_4ds;
    lua["g314_4ds"] = g314_4ds;
    lua["g124_4ds"] = g124_4ds;
    lua["g123_4ds"] = g123_4ds;

    // STA trivectors as multivectors
    lua["g234_4ds_mv"] = g234_4ds_mv;
    lua["g314_4ds_mv"] = g314_4ds_mv;
    lua["g124_4ds_mv"] = g124_4ds_mv;
    lua["g123_4ds_mv"] = g123_4ds_mv;

    // STA pseudoscalar
    lua["I_4ds"] = I_4ds;
    lua["g1234_4ds"] = g1234_4ds;
    lua["I_4ds_mv"] = I_4ds_mv;
    lua["I_4ds_mv_e"] = I_4ds_mv_e;

    // STA scalar one
    lua["one_4ds"] = one_4ds;
    lua["one_4ds_mv"] = one_4ds_mv;
    lua["one_4ds_mv_e"] = one_4ds_mv_e;
}

////////////////////////////////////////////////////////////////////////////////
// Lua-level forwarders for the PGA contraction / expansion / projection /
// distance products. These are generic C++ `if constexpr` templates that cannot
// be bound by sol::resolve, so - exactly as the Python wrapper does in
// ga_py/python/ga_py/__init__.py - they are defined as thin Lua functions over the
// already-registered primitives (rwdg, wdg, the bulk/weight duals, unitize, the
// norms, att, gr). One definition works for every operand type because those
// primitives are themselves overloaded.
//
// pga2dp is odd-dimensional (single bulk_dual / weight_dual); pga3dp is
// even-dimensional (left/right l_bulk_dual / r_bulk_dual / l_weight_dual /
// r_weight_dual). The contraction/expansion side must match the operand side.
////////////////////////////////////////////////////////////////////////////////
void register_forwarders(sol::state& lua)
{
    lua.script(R"lua(
        -- PGA 2dp bulk/weight contractions:  rwdg(dual(a), b) / rwdg(a, dual(b))
        function l_bulk_contract2dp(a, b)   return rwdg(bulk_dual(a), b)   end
        function l_weight_contract2dp(a, b) return rwdg(weight_dual(a), b) end
        function r_bulk_contract2dp(a, b)   return rwdg(a, bulk_dual(b))   end
        function r_weight_contract2dp(a, b) return rwdg(a, weight_dual(b)) end

        -- PGA 2dp bulk/weight expansions:  wdg(dual(a), b) / wdg(a, dual(b))
        function l_bulk_expand2dp(a, b)   return wdg(bulk_dual(a), b)   end
        function l_weight_expand2dp(a, b) return wdg(weight_dual(a), b) end
        function r_bulk_expand2dp(a, b)   return wdg(a, bulk_dual(b))   end
        function r_weight_expand2dp(a, b) return wdg(a, weight_dual(b)) end

        -- PGA 3dp bulk/weight contractions (left/right complement variants)
        function l_bulk_contract3dp(a, b)   return rwdg(l_bulk_dual(a), b)   end
        function l_weight_contract3dp(a, b) return rwdg(l_weight_dual(a), b) end
        function r_bulk_contract3dp(a, b)   return rwdg(a, r_bulk_dual(b))   end
        function r_weight_contract3dp(a, b) return rwdg(a, r_weight_dual(b)) end

        -- PGA 3dp bulk/weight expansions
        function l_bulk_expand3dp(a, b)   return wdg(l_bulk_dual(a), b)   end
        function l_weight_expand3dp(a, b) return wdg(l_weight_dual(a), b) end
        function r_bulk_expand3dp(a, b)   return wdg(a, r_bulk_dual(b))   end
        function r_weight_expand3dp(a, b) return wdg(a, r_weight_dual(b)) end

        -- "is the weight squared norm essentially 1?" guard (matches ga_py's _EPS)
        local _EPS = 1e-9
        local function _unitize_if_needed(p)
            local n = weight_nrm_sq(p)
            if n > _EPS and n ~= 1.0 then p = unitize(p) end
            return p
        end

        -- projections of the lower-grade a onto the larger-grade b  (gr(a) < gr(b))
        function ortho_proj2dp(a, b)   return _unitize_if_needed(rwdg(b, r_weight_expand2dp(a, b))) end
        function central_proj2dp(a, b) return _unitize_if_needed(rwdg(b, r_bulk_expand2dp(a, b)))   end
        function ortho_antiproj2dp(a, b) return wdg(b, r_weight_contract2dp(a, b)) end
        function ortho_proj3dp(a, b)   return _unitize_if_needed(rwdg(b, r_weight_expand3dp(a, b))) end
        function central_proj3dp(a, b) return _unitize_if_needed(rwdg(b, r_bulk_expand3dp(a, b)))   end
        function ortho_antiproj3dp(a, b) return wdg(b, r_weight_contract3dp(a, b)) end

        -- Euclidean distance -> dualnum(homogeneous_magnitude, weight). The C++
        -- `if constexpr (gr(a)+gr(b) == n)` is dispatched here via gr().
        function dist2dp(a, b)
            local c1 = to_val(weight_nrm(wdg(a, att(b))))
            local c0
            if gr(a) + gr(b) == 3 then c0 = to_val(rwdg(a, b))
            else                       c0 = to_val(bulk_nrm(att(wdg(a, b)))) end
            return dualnum2dp.new(c0, c1)
        end
        function dist3dp(a, b)
            local c1 = to_val(weight_nrm(wdg(a, att(b))))
            local c0
            if gr(a) + gr(b) == 4 then c0 = to_val(rwdg(a, b))
            else                       c0 = to_val(bulk_nrm(att(wdg(a, b)))) end
            return dualnum3dp.new(c0, c1)
        end

        -- STA4D expansions (generic C++ templates, mirroring ga_py's forwarders):
        --   l_expand4ds(a,b) = wdg(l_dual(a), b) ;  r_expand4ds(a,b) = wdg(a, r_dual(b))
        function l_expand4ds(a, b) return wdg(l_dual(a), b) end
        function r_expand4ds(a, b) return wdg(a, r_dual(b)) end
    )lua");
}

////////////////////////////////////////////////////////////////////////////////
// PGA geometric convenience types: point / vector / line / plane. Each is a
// subclass of an already-registered value type (e.g. Point2dp : Vec2dp); it adds
// convenient coordinate constructors, and sol::base_classes lets every operation
// bound for the base accept the derived type (results come back as the base type,
// matching GA semantics). Must run AFTER the base value types are registered.
////////////////////////////////////////////////////////////////////////////////
void register_convenience_types(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::ega;
    using namespace hd::ga::pga;

    // 2d points / vectors (subclasses of the EGA Vec2d) and the 2dp line
    lua.new_usertype<point2d>(
        "point2d", sol::constructors<point2d(), point2d(value_t, value_t)>(),
        sol::base_classes, sol::bases<vec2d>(), sol::meta_function::to_string,
        [](point2d const& p) { return fmt::format("{}", p); });

    lua.new_usertype<vector2d>(
        "vector2d", sol::constructors<vector2d(), vector2d(value_t, value_t)>(),
        sol::base_classes, sol::bases<vec2d>(), sol::meta_function::to_string,
        [](vector2d const& v) { return fmt::format("{}", v); });

    lua.new_usertype<point2dp>(
        "point2dp", sol::constructors<point2dp(), point2dp(value_t, value_t, value_t)>(),
        sol::base_classes, sol::bases<vec2dp>(), sol::meta_function::to_string,
        [](point2dp const& p) { return fmt::format("{}", p); });

    lua.new_usertype<line2d>(
        "line2d", sol::constructors<line2d(), line2d(value_t, value_t, value_t)>(),
        sol::base_classes, sol::bases<bivec2dp>(), sol::meta_function::to_string,
        [](line2d const& l) { return fmt::format("{}", l); });

    // 3d points / vectors (subclasses of the EGA Vec3d), the 3dp line and plane
    lua.new_usertype<point3d>(
        "point3d", sol::constructors<point3d(), point3d(value_t, value_t, value_t)>(),
        sol::base_classes, sol::bases<vec3d>(), sol::meta_function::to_string,
        [](point3d const& p) { return fmt::format("{}", p); });

    lua.new_usertype<vector3d>(
        "vector3d", sol::constructors<vector3d(), vector3d(value_t, value_t, value_t)>(),
        sol::base_classes, sol::bases<vec3d>(), sol::meta_function::to_string,
        [](vector3d const& v) { return fmt::format("{}", v); });

    lua.new_usertype<point3dp>(
        "point3dp",
        sol::constructors<point3dp(), point3dp(value_t, value_t, value_t, value_t)>(),
        sol::base_classes, sol::bases<vec3dp>(), sol::meta_function::to_string,
        [](point3dp const& p) { return fmt::format("{}", p); });

    lua.new_usertype<line3d>(
        "line3d",
        sol::constructors<line3d(),
                          line3d(value_t, value_t, value_t, value_t, value_t, value_t)>(),
        sol::base_classes, sol::bases<bivec3dp>(), sol::meta_function::to_string,
        [](line3d const& l) { return fmt::format("{}", l); });

    lua.new_usertype<plane3d>(
        "plane3d",
        sol::constructors<plane3d(), plane3d(value_t, value_t, value_t, value_t)>(),
        sol::base_classes, sol::bases<trivec3dp>(), sol::meta_function::to_string,
        [](plane3d const& p) { return fmt::format("{}", p); });
}
////////////////////////////////////////////////////////////////////////////////
// PGA physics pure-data structs (mirrors the ga_py POD bindings): pose,
// kin_state, joint_state, loop_constraint - each exposed with its constructors,
// read/write fields and a to_string. The scoped enums used as fields (joint type,
// constraint type) are registered first. The stateful *_system classes and the
// inertia matrix type are NOT bound (matching ga_py's POD-only policy), and the
// dynamics free functions (make_*_body / get_*_inertia / compute_omega_dot) are
// intentionally skipped.
////////////////////////////////////////////////////////////////////////////////
void register_physics_pods(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::pga;

    // scoped enums used as POD fields
    lua.new_enum<joint2dp>("joint2dp", {{"free", joint2dp::free},
                                        {"revolute", joint2dp::revolute},
                                        {"prismatic", joint2dp::prismatic}});
    lua.new_enum<joint3dp>("joint3dp", {{"free", joint3dp::free},
                                        {"revolute", joint3dp::revolute},
                                        {"prismatic", joint3dp::prismatic}});
    lua.new_enum<constraint2dp>("constraint2dp",
                                {{"coincidence", constraint2dp::coincidence}});
    lua.new_enum<constraint3dp>("constraint3dp",
                                {{"coincidence", constraint3dp::coincidence}});

    // pose: rigid placement (origin + orientation)
    lua.new_usertype<pose2dp>(
        "pose2dp", sol::constructors<pose2dp(), pose2dp(vec2dp const&, value_t)>(),
        "origin", &pose2dp::origin, "phi", &pose2dp::phi, sol::meta_function::to_string,
        [](pose2dp const& p) { return fmt::format("{}", p); });

    lua.new_usertype<pose3dp>(
        "pose3dp", sol::constructors<pose3dp(), pose3dp(vec3dp const&, vec3dp const&)>(),
        "origin", &pose3dp::origin, "rot", &pose3dp::rot, sol::meta_function::to_string,
        [](pose3dp const& p) { return fmt::format("{}", p); });

    // kin_state: velocity / acceleration kinematic state
    lua.new_usertype<kin_state2dp>(
        "kin_state2dp",
        sol::constructors<kin_state2dp(),
                          kin_state2dp(vec2dp const&, vec2dp const&, value_t, value_t)>(),
        "vel", &kin_state2dp::vel, "acc", &kin_state2dp::acc, "omega",
        &kin_state2dp::omega, "alpha", &kin_state2dp::alpha,
        sol::meta_function::to_string,
        [](kin_state2dp const& k) { return fmt::format("{}", k); });

    lua.new_usertype<kin_state3dp>(
        "kin_state3dp",
        sol::constructors<kin_state3dp(), kin_state3dp(vec3dp const&, vec3dp const&,
                                                       vec3dp const&, vec3dp const&)>(),
        "vel", &kin_state3dp::vel, "acc", &kin_state3dp::acc, "omega",
        &kin_state3dp::omega, "alpha", &kin_state3dp::alpha,
        sol::meta_function::to_string,
        [](kin_state3dp const& k) { return fmt::format("{}", k); });

    // joint_state: per-joint configuration + dynamics parameters
    lua.new_usertype<joint_state2dp>(
        "joint_state2dp",
        sol::constructors<joint_state2dp(),
                          joint_state2dp(joint2dp const&, vec2dp const&, mvec2dp_u const&,
                                         value_t, value_t, value_t, value_t, value_t)>(),
        "type", &joint_state2dp::type, "screw_b", &joint_state2dp::screw_b, "rest",
        &joint_state2dp::rest, "phi", &joint_state2dp::phi, "omega",
        &joint_state2dp::omega, "stiffness", &joint_state2dp::stiffness, "damping",
        &joint_state2dp::damping, "q_rest", &joint_state2dp::q_rest,
        sol::meta_function::to_string,
        [](joint_state2dp const& j) { return fmt::format("{}", j); });

    lua.new_usertype<joint_state3dp>(
        "joint_state3dp",
        sol::constructors<joint_state3dp(),
                          joint_state3dp(joint3dp const&, bivec3dp const&,
                                         mvec3dp_e const&, value_t, value_t, value_t,
                                         value_t, value_t)>(),
        "type", &joint_state3dp::type, "screw_b", &joint_state3dp::screw_b, "rest",
        &joint_state3dp::rest, "phi", &joint_state3dp::phi, "omega",
        &joint_state3dp::omega, "stiffness", &joint_state3dp::stiffness, "damping",
        &joint_state3dp::damping, "q_rest", &joint_state3dp::q_rest,
        sol::meta_function::to_string,
        [](joint_state3dp const& j) { return fmt::format("{}", j); });

    // loop_constraint: kinematic loop-closure spec
    lua.new_usertype<loop_constraint2dp>(
        "loop_constraint2dp",
        sol::constructors<loop_constraint2dp(),
                          loop_constraint2dp(int, vec2dp const&, int, vec2dp const&,
                                             constraint2dp const&)>(),
        "frame_a", &loop_constraint2dp::frame_a, "anchor_a",
        &loop_constraint2dp::anchor_a, "frame_b", &loop_constraint2dp::frame_b,
        "anchor_b", &loop_constraint2dp::anchor_b, "type", &loop_constraint2dp::type,
        sol::meta_function::to_string,
        [](loop_constraint2dp const& c) { return fmt::format("{}", c); });

    lua.new_usertype<loop_constraint3dp>(
        "loop_constraint3dp",
        sol::constructors<loop_constraint3dp(),
                          loop_constraint3dp(int, vec3dp const&, int, vec3dp const&,
                                             constraint3dp const&)>(),
        "frame_a", &loop_constraint3dp::frame_a, "anchor_a",
        &loop_constraint3dp::anchor_a, "frame_b", &loop_constraint3dp::frame_b,
        "anchor_b", &loop_constraint3dp::anchor_b, "type", &loop_constraint3dp::type,
        sol::meta_function::to_string,
        [](loop_constraint3dp const& c) { return fmt::format("{}", c); });
}

void register_4ds_types(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::sta;

    // Basic scalar type
    lua.new_usertype<scalar4ds>(
        "scalar4ds",
        sol::constructors<scalar4ds(), scalar4ds(value_t const&), scalar4ds(value_t&&)>(),
        "copy", [](const scalar4ds& obj) { return scalar4ds(obj); },
        sol::meta_function::to_string,
        [](const scalar4ds& s) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Scalar4ds({})", suppress_negative_zero(double(s)));
        },
        sol::meta_function::unary_minus, sol::resolve<scalar4ds(scalar4ds)>(operator-),
        sol::meta_function::addition,
        sol::resolve<scalar4ds(scalar4ds, scalar4ds)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<scalar4ds(scalar4ds, scalar4ds)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<scalar4ds(scalar4ds, scalar4ds)>(operator*),
                      sol::resolve<scalar4ds(scalar4ds, value_t)>(operator*),
                      sol::resolve<scalar4ds(value_t, scalar4ds)>(operator*)),
        sol::meta_function::division,
        sol::resolve<scalar4ds(scalar4ds, value_t)>(operator/));

    // Vector type with x,y,z,w components
    lua.new_usertype<vec4ds>(
        "vec4ds",
        sol::constructors<vec4ds(), vec4ds(value_t, value_t, value_t, value_t),
                          vec4ds(vec4ds const&), vec4ds(vec4ds&&)>(),
        "copy", [](const vec4ds& obj) { return vec4ds(obj); },
        // component access
        "x", &vec4ds::x, "y", &vec4ds::y, "z", &vec4ds::z, "w", &vec4ds::w,
        sol::meta_function::to_string,
        [](const vec4ds& v) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("Vec4ds({},{},{},{})", suppress_negative_zero(v.x),
                               suppress_negative_zero(v.y), suppress_negative_zero(v.z),
                               suppress_negative_zero(v.w));
        },
        sol::meta_function::unary_minus, sol::resolve<vec4ds(vec4ds const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<vec4ds(vec4ds const&, vec4ds const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<vec4ds(vec4ds const&, vec4ds const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec4ds(vec4ds const&, value_t)>(operator*),
                      sol::resolve<vec4ds(value_t, vec4ds const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec4ds(vec4ds const&, value_t)>(operator/),
        sol::meta_function::power_of,
        sol::resolve<bivec4ds(vec4ds const&, vec4ds const&)>(wdg));

    // Bivector type with vx,vy,vz,mx,my,mz components
    lua.new_usertype<bivec4ds>(
        "bivec4ds",
        sol::constructors<bivec4ds(),
                          bivec4ds(value_t, value_t, value_t, value_t, value_t, value_t),
                          bivec4ds(bivec4ds const&), bivec4ds(bivec4ds&&)>(),
        "copy", [](const bivec4ds& obj) { return bivec4ds(obj); },
        // component access
        "vx", &bivec4ds::vx, "vy", &bivec4ds::vy, "vz", &bivec4ds::vz, "mx",
        &bivec4ds::mx, "my", &bivec4ds::my, "mz", &bivec4ds::mz,
        sol::meta_function::to_string,
        [](const bivec4ds& bv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "BiVec4ds({},{},{},{},{},{})", suppress_negative_zero(bv.vx),
                suppress_negative_zero(bv.vy), suppress_negative_zero(bv.vz),
                suppress_negative_zero(bv.mx), suppress_negative_zero(bv.my),
                suppress_negative_zero(bv.mz));
        },
        sol::meta_function::unary_minus,
        sol::resolve<bivec4ds(bivec4ds const&)>(operator-), sol::meta_function::addition,
        sol::resolve<bivec4ds(bivec4ds const&, bivec4ds const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<bivec4ds(bivec4ds const&, bivec4ds const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<bivec4ds(bivec4ds const&, value_t)>(operator*),
                      sol::resolve<bivec4ds(value_t, bivec4ds const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<bivec4ds(bivec4ds const&, value_t)>(operator/));

    // Trivector type with x,y,z,w components
    lua.new_usertype<trivec4ds>(
        "trivec4ds",
        sol::constructors<trivec4ds(), trivec4ds(value_t, value_t, value_t, value_t),
                          trivec4ds(trivec4ds const&), trivec4ds(trivec4ds&&)>(),
        "copy", [](const trivec4ds& obj) { return trivec4ds(obj); },
        // component access
        "x", &trivec4ds::x, "y", &trivec4ds::y, "z", &trivec4ds::z, "w", &trivec4ds::w,
        sol::meta_function::to_string,
        [](const trivec4ds& tv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("TriVec4ds({},{},{},{})", suppress_negative_zero(tv.x),
                               suppress_negative_zero(tv.y), suppress_negative_zero(tv.z),
                               suppress_negative_zero(tv.w));
        },
        sol::meta_function::unary_minus,
        sol::resolve<trivec4ds(trivec4ds const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<trivec4ds(trivec4ds const&, trivec4ds const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<trivec4ds(trivec4ds const&, trivec4ds const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<trivec4ds(trivec4ds const&, value_t)>(operator*),
                      sol::resolve<trivec4ds(value_t, trivec4ds const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<trivec4ds(trivec4ds const&, value_t)>(operator/));

    // Pseudoscalar type
    lua.new_usertype<pscalar4ds>(
        "pscalar4ds",
        sol::constructors<pscalar4ds(), pscalar4ds(value_t const&),
                          pscalar4ds(value_t&&)>(),
        "copy", [](const pscalar4ds& obj) { return pscalar4ds(obj); },
        sol::meta_function::to_string,
        [](const pscalar4ds& ps) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("PScalar4ds({})", suppress_negative_zero(double(ps)));
        },
        sol::meta_function::unary_minus, sol::resolve<pscalar4ds(pscalar4ds)>(operator-),
        sol::meta_function::addition,
        sol::resolve<pscalar4ds(pscalar4ds, pscalar4ds)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<pscalar4ds(pscalar4ds, pscalar4ds)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<pscalar4ds(pscalar4ds, value_t)>(operator*),
                      sol::resolve<pscalar4ds(value_t, pscalar4ds)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar4ds(pscalar4ds, value_t)>(operator/));

    // Dual numbers (for PGA-specific calculations)
    lua.new_usertype<dualnum4ds>(
        "dualnum4ds",
        sol::constructors<dualnum4ds(), dualnum4ds(value_t, value_t),
                          dualnum4ds(scalar4ds), dualnum4ds(pscalar4ds),
                          dualnum4ds(scalar4ds, pscalar4ds),
                          dualnum4ds(dualnum4ds const&), dualnum4ds(dualnum4ds&&)>(),
        "copy", [](const dualnum4ds& obj) { return dualnum4ds(obj); },
        // component access
        "c0", &dualnum4ds::c0, "c1", &dualnum4ds::c1, sol::meta_function::to_string,
        [](const dualnum4ds& dn) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format("DualNum4ds({},{})", suppress_negative_zero(dn.c0),
                               suppress_negative_zero(dn.c1));
        },
        sol::meta_function::unary_minus,
        sol::resolve<dualnum4ds(dualnum4ds const&)>(operator-),
        sol::meta_function::addition,
        sol::overload(
            sol::resolve<dualnum4ds(dualnum4ds const&, dualnum4ds const&)>(operator+),
            sol::resolve<dualnum4ds(scalar4ds, dualnum4ds const&)>(operator+),
            sol::resolve<dualnum4ds(dualnum4ds const&, scalar4ds)>(operator+),
            sol::resolve<dualnum4ds(pscalar4ds, dualnum4ds const&)>(operator+),
            sol::resolve<dualnum4ds(dualnum4ds const&, pscalar4ds)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(
            sol::resolve<dualnum4ds(dualnum4ds const&, dualnum4ds const&)>(operator-),
            sol::resolve<dualnum4ds(scalar4ds, dualnum4ds const&)>(operator-),
            sol::resolve<dualnum4ds(dualnum4ds const&, scalar4ds)>(operator-),
            sol::resolve<dualnum4ds(pscalar4ds, dualnum4ds const&)>(operator-),
            sol::resolve<dualnum4ds(dualnum4ds const&, pscalar4ds)>(operator-)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<dualnum4ds(dualnum4ds const&, value_t)>(operator*),
                      sol::resolve<dualnum4ds(value_t, dualnum4ds const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<dualnum4ds(dualnum4ds const&, value_t)>(operator/));

    // PGA 3DP multivector types
    lua.new_usertype<mvec4ds_e>(
        "mvec4ds_e",
        sol::constructors<mvec4ds_e(),
                          mvec4ds_e(value_t, value_t, value_t, value_t, value_t, value_t,
                                    value_t, value_t),
                          mvec4ds_e(mvec4ds_e const&), mvec4ds_e(mvec4ds_e&&),
                          mvec4ds_e(scalar4ds), mvec4ds_e(bivec4ds const&),
                          mvec4ds_e(pscalar4ds), mvec4ds_e(scalar4ds, bivec4ds const&),
                          mvec4ds_e(bivec4ds const&, pscalar4ds),
                          mvec4ds_e(scalar4ds, pscalar4ds),
                          mvec4ds_e(scalar4ds, bivec4ds const&, pscalar4ds)>(),
        "copy", [](const mvec4ds_e& obj) { return mvec4ds_e(obj); },
        // component access
        "c0", &mvec4ds_e::c0, "c1", &mvec4ds_e::c1, "c2", &mvec4ds_e::c2, "c3",
        &mvec4ds_e::c3, "c4", &mvec4ds_e::c4, "c5", &mvec4ds_e::c5, "c6", &mvec4ds_e::c6,
        "c7", &mvec4ds_e::c7, sol::meta_function::to_string,
        [](const mvec4ds_e& mve) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "MVec4ds_E({},{},{},{},{},{},{},{})", suppress_negative_zero(mve.c0),
                suppress_negative_zero(mve.c1), suppress_negative_zero(mve.c2),
                suppress_negative_zero(mve.c3), suppress_negative_zero(mve.c4),
                suppress_negative_zero(mve.c5), suppress_negative_zero(mve.c6),
                suppress_negative_zero(mve.c7));
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec4ds_e(mvec4ds_e const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec4ds_e(mvec4ds_e const&, mvec4ds_e const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec4ds_e(mvec4ds_e const&, mvec4ds_e const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(
            sol::resolve<mvec4ds_e(mvec4ds_e const&, value_t)>(operator*),
            sol::resolve<mvec4ds_e(value_t, mvec4ds_e const&)>(operator*),
            sol::resolve<mvec4ds_e(mvec4ds_e const&, mvec4ds_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec4ds_e(mvec4ds_e const&, value_t)>(operator/));

    lua.new_usertype<mvec4ds_u>(
        "mvec4ds_u",
        sol::constructors<mvec4ds_u(),
                          mvec4ds_u(value_t, value_t, value_t, value_t, value_t, value_t,
                                    value_t, value_t),
                          mvec4ds_u(mvec4ds_u const&), mvec4ds_u(mvec4ds_u&&),
                          mvec4ds_u(vec4ds const&), mvec4ds_u(trivec4ds const&),
                          mvec4ds_u(vec4ds const&, trivec4ds const&)>(),
        "copy", [](const mvec4ds_u& obj) { return mvec4ds_u(obj); },
        // component access
        "c0", &mvec4ds_u::c0, "c1", &mvec4ds_u::c1, "c2", &mvec4ds_u::c2, "c3",
        &mvec4ds_u::c3, "c4", &mvec4ds_u::c4, "c5", &mvec4ds_u::c5, "c6", &mvec4ds_u::c6,
        "c7", &mvec4ds_u::c7, sol::meta_function::to_string,
        [](const mvec4ds_u& mvu) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "MVec4ds_U({},{},{},{},{},{},{},{})", suppress_negative_zero(mvu.c0),
                suppress_negative_zero(mvu.c1), suppress_negative_zero(mvu.c2),
                suppress_negative_zero(mvu.c3), suppress_negative_zero(mvu.c4),
                suppress_negative_zero(mvu.c5), suppress_negative_zero(mvu.c6),
                suppress_negative_zero(mvu.c7));
        },
        sol::meta_function::unary_minus,
        sol::resolve<mvec4ds_u(mvec4ds_u const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec4ds_u(mvec4ds_u const&, mvec4ds_u const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec4ds_u(mvec4ds_u const&, mvec4ds_u const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec4ds_u(mvec4ds_u const&, value_t)>(operator*),
                      sol::resolve<mvec4ds_u(value_t, mvec4ds_u const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec4ds_u(mvec4ds_u const&, value_t)>(operator/));

    lua.new_usertype<mvec4ds>(
        "mvec4ds",
        sol::constructors<mvec4ds(),
                          mvec4ds(value_t, value_t, value_t, value_t, value_t, value_t,
                                  value_t, value_t, value_t, value_t, value_t, value_t,
                                  value_t, value_t, value_t, value_t),
                          mvec4ds(mvec4ds const&), mvec4ds(mvec4ds&&), mvec4ds(scalar4ds),
                          mvec4ds(vec4ds const&), mvec4ds(bivec4ds const&),
                          mvec4ds(trivec4ds const&), mvec4ds(pscalar4ds),
                          mvec4ds(mvec4ds_e const&), mvec4ds(mvec4ds_u const&)>(),
        "copy", [](const mvec4ds& obj) { return mvec4ds(obj); },
        // component access
        "c0", &mvec4ds::c0, "c1", &mvec4ds::c1, "c2", &mvec4ds::c2, "c3", &mvec4ds::c3,
        "c4", &mvec4ds::c4, "c5", &mvec4ds::c5, "c6", &mvec4ds::c6, "c7", &mvec4ds::c7,
        "c8", &mvec4ds::c8, "c9", &mvec4ds::c9, "c10", &mvec4ds::c10, "c11",
        &mvec4ds::c11, "c12", &mvec4ds::c12, "c13", &mvec4ds::c13, "c14", &mvec4ds::c14,
        "c15", &mvec4ds::c15, sol::meta_function::to_string,
        [](const mvec4ds& mv) {
            using hd::ga::detail::suppress_negative_zero;
            return fmt::format(
                "MVec4ds({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})",
                suppress_negative_zero(mv.c0), suppress_negative_zero(mv.c1),
                suppress_negative_zero(mv.c2), suppress_negative_zero(mv.c3),
                suppress_negative_zero(mv.c4), suppress_negative_zero(mv.c5),
                suppress_negative_zero(mv.c6), suppress_negative_zero(mv.c7),
                suppress_negative_zero(mv.c8), suppress_negative_zero(mv.c9),
                suppress_negative_zero(mv.c10), suppress_negative_zero(mv.c11),
                suppress_negative_zero(mv.c12), suppress_negative_zero(mv.c13),
                suppress_negative_zero(mv.c14), suppress_negative_zero(mv.c15));
        },
        sol::meta_function::unary_minus, sol::resolve<mvec4ds(mvec4ds const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec4ds(mvec4ds const&, mvec4ds const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec4ds(mvec4ds const&, mvec4ds const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec4ds(mvec4ds const&, value_t)>(operator*),
                      sol::resolve<mvec4ds(value_t, mvec4ds const&)>(operator*),
                      sol::resolve<mvec4ds(mvec4ds const&, mvec4ds const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec4ds(mvec4ds const&, value_t)>(operator/));
}
