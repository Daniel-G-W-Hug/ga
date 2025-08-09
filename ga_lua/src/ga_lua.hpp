#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#define SOL_ALL_SAFETIES_ON 1
#include "sol/sol.hpp"

#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

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
// Function and constant registration
void register_functions(sol::state& lua);
void register_constants(sol::state& lua);


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
        sol::meta_function::to_string,
        [](const scalar2d& s) { return fmt::format("Scalar2d({})", double(s)); },
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
        "x", &vec2d::x, "y", &vec2d::y, sol::meta_function::to_string,
        [](const vec2d& v) { return fmt::format("Vec2d({}, {})", v.x, v.y); },
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
        sol::meta_function::to_string,
        [](const pscalar2d& ps) { return fmt::format("PScalar2d({})", double(ps)); },
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
                      sol::resolve<scalar2d(pscalar2d, pscalar2d)>(wdg)),
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
        "c0", &mvec2d_e::c0, "c1", &mvec2d_e::c1, sol::meta_function::to_string,
        [](const mvec2d_e& mv) { return fmt::format("MVec2d_E({}, {})", mv.c0, mv.c1); },
        sol::meta_function::unary_minus,
        sol::resolve<mvec2d_e(mvec2d_e const&)>(operator-), sol::meta_function::addition,
        sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator+),
                      sol::resolve<mvec2d_e(mvec2d_e const&, scalar2d)>(operator+),
                      sol::resolve<mvec2d_e(mvec2d_e const&, pscalar2d)>(operator+)),
        sol::meta_function::subtraction,
        sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator-),
                      sol::resolve<mvec2d_e(mvec2d_e const&, scalar2d)>(operator-),
                      sol::resolve<mvec2d_e(mvec2d_e const&, pscalar2d)>(operator-)),
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
        "c0", &mvec2d::c0, "c1", &mvec2d::c1, "c2", &mvec2d::c2, "c3", &mvec2d::c3,
        sol::meta_function::to_string,
        [](const mvec2d& mv) {
            return fmt::format("MVec2d({}, {}, {}, {})", mv.c0, mv.c1, mv.c2, mv.c3);
        },
        sol::meta_function::unary_minus, sol::resolve<mvec2d(mvec2d const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator-),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(wdg),
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(wdg)),
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
        sol::meta_function::to_string,
        [](const scalar3d& s) { return fmt::format("Scalar3d({})", double(s)); },
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
        "x", &vec3d::x, "y", &vec3d::y, "z", &vec3d::z, sol::meta_function::to_string,
        [](const vec3d& v) { return fmt::format("Vec3d({}, {}, {})", v.x, v.y, v.z); },
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
        "x", &bivec3d::x, "y", &bivec3d::y, "z", &bivec3d::z,
        sol::meta_function::to_string,
        [](const bivec3d& bv) {
            return fmt::format("BiVec3d({}, {}, {})", bv.x, bv.y, bv.z);
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
        sol::meta_function::to_string,
        [](const pscalar3d& ps) { return fmt::format("PScalar3d({})", double(ps)); },
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
        "c0", &mvec3d_e::c0, "c1", &mvec3d_e::c1, "c2", &mvec3d_e::c2, "c3",
        &mvec3d_e::c3, sol::meta_function::to_string,
        [](const mvec3d_e& mv) {
            return fmt::format("MVec3d_E({}, {}, {}, {})", mv.c0, mv.c1, mv.c2, mv.c3);
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
        "c0", &mvec3d_u::c0, "c1", &mvec3d_u::c1, "c2", &mvec3d_u::c2, "c3",
        &mvec3d_u::c3, sol::meta_function::to_string,
        [](const mvec3d_u& mv) {
            return fmt::format("MVec3d_U({}, {}, {}, {})", mv.c0, mv.c1, mv.c2, mv.c3);
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
        "c0", &mvec3d::c0, "c1", &mvec3d::c1, "c2", &mvec3d::c2, "c3", &mvec3d::c3, "c4",
        &mvec3d::c4, "c5", &mvec3d::c5, "c6", &mvec3d::c6, "c7", &mvec3d::c7,
        sol::meta_function::to_string,
        [](const mvec3d& mv) {
            return fmt::format("MVec3d({}, {}, {}, {}, {}, {}, {}, {})", mv.c0, mv.c1,
                               mv.c2, mv.c3, mv.c4, mv.c5, mv.c6, mv.c7);
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
        "x", &vec2dp::x, "y", &vec2dp::y, "z", &vec2dp::z, sol::meta_function::to_string,
        [](const vec2dp& v) { return fmt::format("Vec2dp({}, {}, {})", v.x, v.y, v.z); },
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
        "x", &bivec2dp::x, "y", &bivec2dp::y, "z", &bivec2dp::z,
        sol::meta_function::to_string,
        [](const bivec2dp& bv) {
            return fmt::format("BiVec2dp({}, {}, {})", bv.x, bv.y, bv.z);
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
        sol::meta_function::to_string,
        [](const pscalar2dp& ps) { return fmt::format("PScalar2dp({})", double(ps)); },
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
                          dualnum2dp(dualnum2dp const&), dualnum2dp(dualnum2dp&&)>(),
        "c0", &dualnum2dp::c0, "c1", &dualnum2dp::c1, sol::meta_function::to_string,
        [](const dualnum2dp& dn) {
            return fmt::format("DualNum2dp({}, {})", dn.c0, dn.c1);
        },
        sol::meta_function::unary_minus,
        sol::resolve<dualnum2dp(dualnum2dp const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<dualnum2dp(dualnum2dp const&, dualnum2dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<dualnum2dp(dualnum2dp const&, dualnum2dp const&)>(operator-),
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
        "c0", &mvec2dp_e::c0, "c1", &mvec2dp_e::c1, "c2", &mvec2dp_e::c2, "c3",
        &mvec2dp_e::c3, sol::meta_function::to_string,
        [](const mvec2dp_e& mve) {
            return fmt::format("MVec2dp_E({}, {}, {}, {})", mve.c0, mve.c1, mve.c2,
                               mve.c3);
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
        "c0", &mvec2dp_u::c0, "c1", &mvec2dp_u::c1, "c2", &mvec2dp_u::c2, "c3",
        &mvec2dp_u::c3, sol::meta_function::to_string,
        [](const mvec2dp_u& mvu) {
            return fmt::format("MVec2dp_U({}, {}, {}, {})", mvu.c0, mvu.c1, mvu.c2,
                               mvu.c3);
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
        "c0", &mvec2dp::c0, "c1", &mvec2dp::c1, "c2", &mvec2dp::c2, "c3", &mvec2dp::c3,
        "c4", &mvec2dp::c4, "c5", &mvec2dp::c5, "c6", &mvec2dp::c6, "c7", &mvec2dp::c7,
        sol::meta_function::to_string,
        [](const mvec2dp& mv) {
            return fmt::format("MVec2dp({}, {}, {}, {}, {}, {}, {}, {})", mv.c0, mv.c1,
                               mv.c2, mv.c3, mv.c4, mv.c5, mv.c6, mv.c7);
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
        sol::meta_function::to_string,
        [](const scalar3dp& s) { return fmt::format("Scalar3dp({})", double(s)); },
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
        "x", &vec3dp::x, "y", &vec3dp::y, "z", &vec3dp::z, "w", &vec3dp::w,
        sol::meta_function::to_string,
        [](const vec3dp& v) {
            return fmt::format("Vec3dp({}, {}, {}, {})", v.x, v.y, v.z, v.w);
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
        "vx", &bivec3dp::vx, "vy", &bivec3dp::vy, "vz", &bivec3dp::vz, "mx",
        &bivec3dp::mx, "my", &bivec3dp::my, "mz", &bivec3dp::mz,
        sol::meta_function::to_string,
        [](const bivec3dp& bv) {
            return fmt::format("BiVec3dp({}, {}, {}, {}, {}, {})", bv.vx, bv.vy, bv.vz,
                               bv.mx, bv.my, bv.mz);
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
        "x", &trivec3dp::x, "y", &trivec3dp::y, "z", &trivec3dp::z, "w", &trivec3dp::w,
        sol::meta_function::to_string,
        [](const trivec3dp& tv) {
            return fmt::format("TriVec3dp({}, {}, {}, {})", tv.x, tv.y, tv.z, tv.w);
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
        sol::meta_function::to_string,
        [](const pscalar3dp& ps) { return fmt::format("PScalar3dp({})", double(ps)); },
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
                          dualnum3dp(dualnum3dp const&), dualnum3dp(dualnum3dp&&)>(),
        "c0", &dualnum3dp::c0, "c1", &dualnum3dp::c1, sol::meta_function::to_string,
        [](const dualnum3dp& dn) {
            return fmt::format("DualNum3dp({}, {})", dn.c0, dn.c1);
        },
        sol::meta_function::unary_minus,
        sol::resolve<dualnum3dp(dualnum3dp const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<dualnum3dp(dualnum3dp const&, dualnum3dp const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<dualnum3dp(dualnum3dp const&, dualnum3dp const&)>(operator-),
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
        "c0", &mvec3dp_e::c0, "c1", &mvec3dp_e::c1, "c2", &mvec3dp_e::c2, "c3",
        &mvec3dp_e::c3, "c4", &mvec3dp_e::c4, "c5", &mvec3dp_e::c5, "c6", &mvec3dp_e::c6,
        "c7", &mvec3dp_e::c7, sol::meta_function::to_string,
        [](const mvec3dp_e& mve) {
            return fmt::format("MVec3dp_E({}, {}, {}, {}, {}, {}, {}, {})", mve.c0,
                               mve.c1, mve.c2, mve.c3, mve.c4, mve.c5, mve.c6, mve.c7);
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
        "c0", &mvec3dp_u::c0, "c1", &mvec3dp_u::c1, "c2", &mvec3dp_u::c2, "c3",
        &mvec3dp_u::c3, "c4", &mvec3dp_u::c4, "c5", &mvec3dp_u::c5, "c6", &mvec3dp_u::c6,
        "c7", &mvec3dp_u::c7, sol::meta_function::to_string,
        [](const mvec3dp_u& mvu) {
            return fmt::format("MVec3dp_U({}, {}, {}, {}, {}, {}, {}, {})", mvu.c0,
                               mvu.c1, mvu.c2, mvu.c3, mvu.c4, mvu.c5, mvu.c6, mvu.c7);
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
        "c0", &mvec3dp::c0, "c1", &mvec3dp::c1, "c2", &mvec3dp::c2, "c3", &mvec3dp::c3,
        "c4", &mvec3dp::c4, "c5", &mvec3dp::c5, "c6", &mvec3dp::c6, "c7", &mvec3dp::c7,
        "c8", &mvec3dp::c8, "c9", &mvec3dp::c9, "c10", &mvec3dp::c10, "c11",
        &mvec3dp::c11, "c12", &mvec3dp::c12, "c13", &mvec3dp::c13, "c14", &mvec3dp::c14,
        "c15", &mvec3dp::c15, sol::meta_function::to_string,
        [](const mvec3dp& mv) {
            return fmt::format(
                "MVec3dp({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
                mv.c0, mv.c1, mv.c2, mv.c3, mv.c4, mv.c5, mv.c6, mv.c7, mv.c8, mv.c9,
                mv.c10, mv.c11, mv.c12, mv.c13, mv.c14, mv.c15);
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
                         // PGA 2DP types
                         sol::resolve<value_t(scalar2dp)>(nrm_sq),
                         sol::resolve<value_t(vec2dp const&)>(nrm_sq),
                         sol::resolve<value_t(bivec2dp const&)>(nrm_sq),
                         sol::resolve<value_t(pscalar2dp)>(nrm_sq),
                         sol::resolve<value_t(mvec2dp_e const&)>(nrm_sq),
                         sol::resolve<value_t(mvec2dp_u const&)>(nrm_sq),
                         sol::resolve<value_t(mvec2dp const&)>(nrm_sq),
                         sol::resolve<value_t(dualnum2dp const&)>(nrm_sq),
                         // PGA 3DP types
                         sol::resolve<value_t(scalar3dp)>(nrm_sq),
                         sol::resolve<value_t(vec3dp const&)>(nrm_sq),
                         sol::resolve<value_t(bivec3dp const&)>(nrm_sq),
                         sol::resolve<value_t(trivec3dp const&)>(nrm_sq),
                         sol::resolve<value_t(pscalar3dp)>(nrm_sq),
                         sol::resolve<value_t(mvec3dp_e const&)>(nrm_sq),
                         sol::resolve<value_t(mvec3dp_u const&)>(nrm_sq),
                         sol::resolve<value_t(mvec3dp const&)>(nrm_sq),
                         sol::resolve<value_t(dualnum3dp const&)>(nrm_sq)));

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
                                // PGA 2DP types
                                sol::resolve<value_t(scalar2dp)>(nrm),
                                sol::resolve<value_t(vec2dp const&)>(nrm),
                                sol::resolve<value_t(bivec2dp const&)>(nrm),
                                sol::resolve<value_t(pscalar2dp)>(nrm),
                                sol::resolve<value_t(mvec2dp_e const&)>(nrm),
                                sol::resolve<value_t(mvec2dp_u const&)>(nrm),
                                sol::resolve<value_t(mvec2dp const&)>(nrm),
                                sol::resolve<value_t(dualnum2dp const&)>(nrm),
                                // PGA 3DP types
                                sol::resolve<value_t(scalar3dp)>(nrm),
                                sol::resolve<value_t(vec3dp const&)>(nrm),
                                sol::resolve<value_t(bivec3dp const&)>(nrm),
                                sol::resolve<value_t(trivec3dp const&)>(nrm),
                                sol::resolve<value_t(pscalar3dp)>(nrm),
                                sol::resolve<value_t(mvec3dp_e const&)>(nrm),
                                sol::resolve<value_t(mvec3dp_u const&)>(nrm),
                                sol::resolve<value_t(mvec3dp const&)>(nrm),
                                sol::resolve<value_t(dualnum3dp const&)>(nrm)));

    lua.set_function("normalize",
                     sol::overload(sol::resolve<scalar2d(scalar2d)>(normalize),
                                   sol::resolve<vec2d(vec2d const&)>(normalize),
                                   sol::resolve<pscalar2d(pscalar2d)>(normalize),
                                   sol::resolve<scalar3d(scalar3d)>(normalize),
                                   sol::resolve<vec3d(vec3d const&)>(normalize),
                                   sol::resolve<bivec3d(bivec3d const&)>(normalize),
                                   sol::resolve<pscalar3d(pscalar3d)>(normalize)));

    // PGA unitize functions (analogue to EGA normalize)
    lua.set_function("unitize",
                     sol::overload(sol::resolve<vec2dp(vec2dp const&)>(unitize),
                                   sol::resolve<bivec2dp(bivec2dp const&)>(unitize),
                                   sol::resolve<dualnum2dp(dualnum2dp const&)>(unitize),
                                   sol::resolve<vec3dp(vec3dp const&)>(unitize),
                                   sol::resolve<bivec3dp(bivec3dp const&)>(unitize),
                                   sol::resolve<trivec3dp(trivec3dp const&)>(unitize),
                                   sol::resolve<dualnum3dp(dualnum3dp const&)>(unitize)));

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


    lua.set_function("right_bulk_dual",
                     sol::overload(
                         // PGA 3DP bulk_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(right_bulk_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(right_bulk_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(right_bulk_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(right_bulk_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(right_bulk_dual)));

    lua.set_function("right_weight_dual",
                     sol::overload(
                         // PGA 3DP weight_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(right_weight_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(right_weight_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(right_weight_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(right_weight_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(right_weight_dual)));

    lua.set_function("left_bulk_dual",
                     sol::overload(
                         // PGA 3DP bulk_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(left_bulk_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(left_bulk_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(left_bulk_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(left_bulk_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(left_bulk_dual)));

    lua.set_function("left_weight_dual",
                     sol::overload(
                         // PGA 3DP weight_dual operations
                         sol::resolve<pscalar3dp(scalar3dp)>(left_weight_dual),
                         sol::resolve<trivec3dp(vec3dp const&)>(left_weight_dual),
                         sol::resolve<bivec3dp(bivec3dp const&)>(left_weight_dual),
                         sol::resolve<vec3dp(trivec3dp const&)>(left_weight_dual),
                         sol::resolve<scalar3dp(pscalar3dp)>(left_weight_dual)));
    ////////////////////////////////////////////////////////////////////////////////
    // PGA-specific support operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("support2dp", sol::resolve<vec2dp(bivec2dp const&)>(support2dp));
    lua.set_function("support3dp", sol::resolve<vec3dp(bivec3dp const&)>(support3dp));


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
                                sol::resolve<scalar3dp(mvec3dp_e const&)>(gr0),
                                sol::resolve<scalar3dp(mvec3dp const&)>(gr0)));

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
                                // PGA 3DP grade 3 operations
                                sol::resolve<trivec3dp(mvec3dp_u const&)>(gr3),
                                sol::resolve<trivec3dp(mvec3dp const&)>(gr3)));

    // PGA-specific grade 4 operation (only for 3DP)
    lua.set_function("gr4", sol::overload(
                                // PGA grade 4 operations
                                sol::resolve<pscalar3dp(mvec3dp_e const&)>(gr4),
                                sol::resolve<pscalar3dp(mvec3dp const&)>(gr4)));


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
                                   sol::resolve<mvec3dp(mvec3dp const&)>(gr_inv)));

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
                                sol::resolve<mvec3dp(mvec3dp const&)>(rev)));

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
                                 sol::resolve<mvec3dp(mvec3dp const&)>(conj)));

    lua.set_function(
        "inv",
        sol::overload(
            // EGA inverse operations
            sol::resolve<scalar2d(scalar2d)>(inv), sol::resolve<vec2d(vec2d const&)>(inv),
            sol::resolve<pscalar2d(pscalar2d)>(inv),
            sol::resolve<mvec2d_e(mvec2d_e const&)>(inv),
            sol::resolve<mvec2d(mvec2d const&)>(inv),
            sol::resolve<scalar3d(scalar3d)>(inv), sol::resolve<vec3d(vec3d const&)>(inv),
            sol::resolve<bivec3d(bivec3d const&)>(inv),
            sol::resolve<pscalar3d(pscalar3d)>(inv),
            sol::resolve<mvec3d_e(mvec3d_e const&)>(inv),
            sol::resolve<mvec3d_u(mvec3d_u const&)>(inv),
            sol::resolve<mvec3d(mvec3d const&)>(inv),
            // PGA inverse operations (conservative - basic types only)
            sol::resolve<scalar2dp(scalar2dp)>(inv),
            sol::resolve<scalar3dp(scalar3dp)>(inv)));


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
                         sol::resolve<scalar3dp(pscalar3dp, pscalar3dp)>(dot)));


    ////////////////////////////////////////////////////////////////////////////////
    // wedge product (= outer product) (wdg)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "wdg",
        sol::overload(sol::resolve<scalar2d(scalar2d, scalar2d)>(wdg),
                      sol::resolve<vec2d(scalar2d, vec2d const&)>(wdg),
                      sol::resolve<vec2d(vec2d const&, scalar2d)>(wdg),
                      sol::resolve<pscalar2d(scalar2d, pscalar2d)>(wdg),
                      sol::resolve<pscalar2d(pscalar2d, scalar2d)>(wdg),
                      sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(wdg),
                      sol::resolve<scalar2d(vec2d const&, pscalar2d)>(wdg),
                      sol::resolve<scalar2d(pscalar2d, vec2d const&)>(wdg),
                      sol::resolve<scalar2d(pscalar2d, vec2d const&)>(wdg),
                      sol::resolve<mvec2d(scalar2d, mvec2d const&)>(wdg),
                      sol::resolve<mvec2d(vec2d const&, mvec2d const&)>(wdg),
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(wdg),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(wdg),
                      sol::resolve<scalar3d(scalar3d, scalar3d)>(wdg),
                      sol::resolve<vec3d(scalar3d, vec3d const&)>(wdg),
                      sol::resolve<vec3d(vec3d const&, scalar3d)>(wdg),
                      sol::resolve<bivec3d(scalar3d, bivec3d const&)>(wdg),
                      sol::resolve<bivec3d(bivec3d const&, scalar3d)>(wdg),
                      sol::resolve<pscalar3d(scalar3d, pscalar3d)>(wdg),
                      sol::resolve<pscalar3d(pscalar3d, scalar3d)>(wdg),
                      sol::resolve<bivec3d(vec3d const&, vec3d const&)>(wdg),
                      sol::resolve<pscalar3d(vec3d const&, bivec3d const&)>(wdg),
                      sol::resolve<pscalar3d(bivec3d const&, vec3d const&)>(wdg),
                      sol::resolve<scalar3d(bivec3d const&, bivec3d const&)>(wdg),
                      sol::resolve<scalar3d(vec3d const&, pscalar3d)>(wdg),
                      sol::resolve<scalar3d(pscalar3d, vec3d const&)>(wdg),
                      sol::resolve<scalar3d(bivec3d const&, pscalar3d)>(wdg),
                      sol::resolve<scalar3d(pscalar3d, bivec3d const&)>(wdg),
                      sol::resolve<scalar3d(pscalar3d, pscalar3d)>(wdg),
                      sol::resolve<mvec3d(scalar3d, mvec3d const&)>(wdg),
                      sol::resolve<mvec3d(vec3d const&, mvec3d const&)>(wdg),
                      sol::resolve<mvec3d(bivec3d const&, mvec3d const&)>(wdg),
                      sol::resolve<mvec3d(mvec3d const&, vec3d const&)>(wdg),
                      sol::resolve<mvec3d(mvec3d const&, bivec3d const&)>(wdg),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(wdg),
                      // PGA 2DP wedge products (comprehensive)
                      sol::resolve<scalar2dp(scalar2dp, scalar2dp)>(wdg),
                      sol::resolve<vec2dp(scalar2dp, vec2dp const&)>(wdg),
                      sol::resolve<vec2dp(vec2dp const&, scalar2dp)>(wdg),
                      sol::resolve<bivec2dp(scalar2dp, bivec2dp const&)>(wdg),
                      sol::resolve<bivec2dp(bivec2dp const&, scalar2dp)>(wdg),
                      sol::resolve<pscalar2dp(scalar2dp, pscalar2dp)>(wdg),
                      sol::resolve<pscalar2dp(pscalar2dp, scalar2dp)>(wdg),
                      sol::resolve<bivec2dp(vec2dp const&, vec2dp const&)>(wdg),
                      sol::resolve<pscalar2dp(bivec2dp const&, vec2dp const&)>(wdg),
                      sol::resolve<pscalar2dp(vec2dp const&, bivec2dp const&)>(wdg),
                      sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(wdg),
                      sol::resolve<scalar2dp(pscalar2dp, pscalar2dp)>(wdg),
                      sol::resolve<scalar2dp(pscalar2dp, bivec2dp const&)>(wdg),
                      sol::resolve<scalar2dp(bivec2dp const&, pscalar2dp)>(wdg),
                      sol::resolve<scalar2dp(pscalar2dp, vec2dp const&)>(wdg),
                      sol::resolve<scalar2dp(vec2dp const&, pscalar2dp)>(wdg),
                      sol::resolve<scalar2dp(bivec2dp const&, bivec2dp const&)>(wdg),
                      // PGA 3DP wedge products (comprehensive)
                      sol::resolve<scalar3dp(scalar3dp, scalar3dp)>(wdg),
                      sol::resolve<vec3dp(scalar3dp, vec3dp const&)>(wdg),
                      sol::resolve<vec3dp(vec3dp const&, scalar3dp)>(wdg),
                      sol::resolve<bivec3dp(scalar3dp, bivec3dp const&)>(wdg),
                      sol::resolve<bivec3dp(bivec3dp const&, scalar3dp)>(wdg),
                      sol::resolve<trivec3dp(scalar3dp, trivec3dp const&)>(wdg),
                      sol::resolve<trivec3dp(trivec3dp const&, scalar3dp)>(wdg),
                      sol::resolve<pscalar3dp(scalar3dp, pscalar3dp)>(wdg),
                      sol::resolve<pscalar3dp(pscalar3dp, scalar3dp)>(wdg),
                      sol::resolve<bivec3dp(vec3dp const&, vec3dp const&)>(wdg),
                      sol::resolve<trivec3dp(bivec3dp const&, vec3dp const&)>(wdg),
                      sol::resolve<trivec3dp(vec3dp const&, bivec3dp const&)>(wdg),
                      sol::resolve<pscalar3dp(bivec3dp const&, bivec3dp const&)>(wdg),
                      sol::resolve<pscalar3dp(trivec3dp const&, vec3dp const&)>(wdg),
                      sol::resolve<pscalar3dp(vec3dp const&, trivec3dp const&)>(wdg),
                      sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(wdg)));


    ////////////////////////////////////////////////////////////////////////////////
    // regressive wedge product (= outer product for complements) (rwdg)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "rwdg",
        sol::overload(sol::resolve<scalar2d(vec2d const&, vec2d const&)>(rwdg),
                      sol::resolve<vec2d(pscalar2d, vec2d const&)>(rwdg),
                      sol::resolve<vec2d(vec2d const&, pscalar2d)>(rwdg),
                      sol::resolve<scalar2d(scalar2d, pscalar2d)>(rwdg),
                      sol::resolve<scalar2d(pscalar2d, scalar2d)>(rwdg),
                      sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(rwdg),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(rwdg),
                      sol::resolve<scalar3d(vec3d const&, bivec3d const&)>(rwdg),
                      sol::resolve<scalar3d(bivec3d const&, vec3d const&)>(rwdg),
                      sol::resolve<vec3d(bivec3d const&, bivec3d const&)>(rwdg),
                      sol::resolve<vec3d(pscalar3d, vec3d const&)>(rwdg),
                      sol::resolve<vec3d(vec3d const&, pscalar3d)>(rwdg),
                      sol::resolve<bivec3d(pscalar3d, bivec3d const&)>(rwdg),
                      sol::resolve<bivec3d(bivec3d const&, pscalar3d)>(rwdg),
                      sol::resolve<scalar3d(scalar3d, pscalar3d)>(rwdg),
                      sol::resolve<scalar3d(pscalar3d, scalar3d)>(rwdg),
                      sol::resolve<pscalar3d(pscalar3d, pscalar3d)>(rwdg),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(rwdg),
                      // PGA 2DP regressive wedge products (complete - all 7 functions)
                      sol::resolve<mvec2dp(mvec2dp const&, mvec2dp const&)>(rwdg),
                      sol::resolve<bivec2dp(pscalar2dp, bivec2dp const&)>(rwdg),
                      sol::resolve<bivec2dp(bivec2dp const&, pscalar2dp)>(rwdg),
                      sol::resolve<vec2dp(bivec2dp const&, bivec2dp const&)>(rwdg),
                      sol::resolve<scalar2dp(bivec2dp const&, vec2dp const&)>(rwdg),
                      sol::resolve<scalar2dp(vec2dp const&, bivec2dp const&)>(rwdg),
                      sol::resolve<scalar2dp(vec2dp const&, vec2dp const&)>(rwdg),
                      // PGA 3DP regressive wedge products (complete - all 12 functions)
                      sol::resolve<mvec3dp(mvec3dp const&, mvec3dp const&)>(rwdg),
                      sol::resolve<trivec3dp(pscalar3dp, trivec3dp const&)>(rwdg),
                      sol::resolve<trivec3dp(trivec3dp const&, pscalar3dp)>(rwdg),
                      sol::resolve<bivec3dp(trivec3dp const&, trivec3dp const&)>(rwdg),
                      sol::resolve<vec3dp(trivec3dp const&, bivec3dp const&)>(rwdg),
                      sol::resolve<vec3dp(bivec3dp const&, trivec3dp const&)>(rwdg),
                      sol::resolve<scalar3dp(trivec3dp const&, vec3dp const&)>(rwdg),
                      sol::resolve<scalar3dp(vec3dp const&, trivec3dp const&)>(rwdg),
                      sol::resolve<scalar3dp(bivec3dp const&, bivec3dp const&)>(rwdg),
                      sol::resolve<scalar3dp(bivec3dp const&, vec3dp const&)>(rwdg),
                      sol::resolve<scalar3dp(vec3dp const&, bivec3dp const&)>(rwdg),
                      sol::resolve<scalar3dp(vec3dp const&, vec3dp const&)>(rwdg)));

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
                         sol::resolve<bivec3dp(bivec3dp const&, bivec3dp const&)>(cmt),
                         sol::resolve<vec3dp(bivec3dp const&, vec3dp const&)>(cmt),
                         sol::resolve<vec3dp(vec3dp const&, bivec3dp const&)>(cmt),
                         sol::resolve<bivec3dp(vec3dp const&, vec3dp const&)>(cmt)));

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

    lua.set_function("support2dp", sol::resolve<vec2dp(bivec2dp const&)>(support2dp));
    lua.set_function("support3dp",
                     sol::overload(sol::resolve<vec3dp(bivec3dp const&)>(support3dp),
                                   sol::resolve<vec3dp(trivec3dp const&)>(support3dp)));

    ////////////////////////////////////////////////////////////////////////////////
    // angles and rotations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "angle",
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
            sol::resolve<value_t(trivec3dp const&, trivec3dp const&)>(angle)));

    lua.set_function("angle_to_re", sol::resolve<value_t(mvec2d_e const&)>(angle_to_re));

    lua.set_function("exp",
                     sol::overload(sol::resolve<mvec2d_e(pscalar2d, value_t)>(exp),
                                   sol::resolve<mvec3d_e(bivec3d const&, value_t)>(exp)));

    lua.set_function(
        "get_rotor",
        sol::overload(sol::resolve<mvec2d_e(pscalar2d, value_t)>(get_rotor),
                      sol::resolve<mvec3d_e(bivec3d const&, value_t)>(get_rotor)));

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

    lua.set_function("left_dual",
                     sol::overload(sol::resolve<scalar2d(pscalar2d)>(left_dual),
                                   sol::resolve<vec2d(vec2d const&)>(left_dual),
                                   sol::resolve<pscalar2d(scalar2d)>(left_dual),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(left_dual),
                                   sol::resolve<mvec2d(mvec2d const&)>(left_dual)));

    lua.set_function("right_dual",
                     sol::overload(sol::resolve<scalar2d(pscalar2d)>(right_dual),
                                   sol::resolve<vec2d(vec2d const&)>(right_dual),
                                   sol::resolve<pscalar2d(scalar2d)>(right_dual),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(right_dual),
                                   sol::resolve<mvec2d(mvec2d const&)>(right_dual)));

    ////////////////////////////////////////////////////////////////////////////////
    // complement operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("rcmpl",
                     sol::overload( // EGA 2D rcmpl functions
                         sol::resolve<pscalar2d(scalar2d)>(rcmpl),
                         sol::resolve<vec2d(vec2d const&)>(rcmpl),
                         sol::resolve<scalar2d(pscalar2d)>(rcmpl),
                         sol::resolve<mvec2d_e(mvec2d_e const&)>(rcmpl),
                         sol::resolve<mvec2d(mvec2d const&)>(rcmpl),
                         // PGA 3DP rcmpl functions (8 total)
                         sol::resolve<pscalar3dp(scalar3dp)>(rcmpl),
                         sol::resolve<trivec3dp(vec3dp const&)>(rcmpl),
                         sol::resolve<bivec3dp(bivec3dp const&)>(rcmpl),
                         sol::resolve<vec3dp(trivec3dp const&)>(rcmpl),
                         sol::resolve<scalar3dp(pscalar3dp)>(rcmpl),
                         sol::resolve<mvec3dp_u(mvec3dp_u const&)>(rcmpl),
                         sol::resolve<mvec3dp_e(mvec3dp_e const&)>(rcmpl),
                         sol::resolve<mvec3dp(mvec3dp const&)>(rcmpl)));

    lua.set_function("lcmpl",
                     sol::overload( // EGA 2D lcmpl functions
                         sol::resolve<pscalar2d(scalar2d)>(lcmpl),
                         sol::resolve<vec2d(vec2d const&)>(lcmpl),
                         sol::resolve<scalar2d(pscalar2d)>(lcmpl),
                         sol::resolve<mvec2d_e(mvec2d_e const&)>(lcmpl),
                         sol::resolve<mvec2d(mvec2d const&)>(lcmpl),
                         // PGA 3DP lcmpl functions (8 total)
                         sol::resolve<pscalar3dp(scalar3dp)>(lcmpl),
                         sol::resolve<trivec3dp(vec3dp const&)>(lcmpl),
                         sol::resolve<bivec3dp(bivec3dp const&)>(lcmpl),
                         sol::resolve<vec3dp(trivec3dp const&)>(lcmpl),
                         sol::resolve<scalar3dp(pscalar3dp)>(lcmpl),
                         sol::resolve<mvec3dp_u(mvec3dp_u const&)>(lcmpl),
                         sol::resolve<mvec3dp_e(mvec3dp_e const&)>(lcmpl),
                         sol::resolve<mvec3dp(mvec3dp const&)>(lcmpl)));

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
            sol::resolve<vec3dp(vec3dp const&, trivec3dp const&)>(project_onto)));

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
            sol::resolve<vec3dp(vec3dp const&, trivec3dp const&)>(reject_from)));

    lua.set_function(
        "reflect_on",
        sol::overload( // EGA reflection functions
            sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on),
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on),
            sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reflect_on),
            sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(reflect_on),
            // PGA 2DP reflection functions
            sol::resolve<vec2dp(vec2dp const&, bivec2dp const&)>(reflect_on),
            sol::resolve<bivec2dp(bivec2dp const&, bivec2dp const&)>(reflect_on),
            // PGA 3DP reflection functions
            sol::resolve<vec3dp(vec3dp const&, trivec3dp const&)>(reflect_on),
            sol::resolve<bivec3dp(bivec3dp const&, trivec3dp const&)>(reflect_on),
            sol::resolve<trivec3dp(trivec3dp const&, trivec3dp const&)>(reflect_on)));

    lua.set_function(
        "reflect_on_vec",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on_vec),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on_vec)));

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
    lua["x_axis_2d"] = x_axis_2d;
    lua["y_axis_2d"] = y_axis_2d;
    lua["origin_2d"] = origin_2d;

    lua["e1m_2d"] = e1m_2d; // as 2d multivector
    lua["e2m_2d"] = e2m_2d;

    lua["I_2d"] = I_2d; // as pscalar2d
    lua["e12_2d"] = I_2d;
    lua["Im_2d"] = Im_2d;     // as 2d multivector
    lua["Im_2d_E"] = Im_2d_E; // as even grade 2d multivector

    // 3d constants
    lua["e1_3d"] = e1_3d; // as 3d vector
    lua["e2_3d"] = e2_3d;
    lua["e3_3d"] = e3_3d;
    lua["x_axis_3d"] = x_axis_3d;
    lua["y_axis_3d"] = y_axis_3d;
    lua["z_axis_3d"] = z_axis_3d;
    lua["origin_3d"] = origin_3d;

    lua["e1m_3d"] = e1m_3d; // as 3d multivector
    lua["e2m_3d"] = e2m_3d;
    lua["e3m_3d"] = e3m_3d;

    lua["e23_3d"] = e23_3d; // as 3d bivector
    lua["e31_3d"] = e31_3d;
    lua["e12_3d"] = e12_3d;
    lua["yz_plane_3d"] = yz_plane_3d;
    lua["zx_plane_3d"] = zx_plane_3d;
    lua["xy_plane_3d"] = xy_plane_3d;


    lua["e23m_3d"] = e23m_3d; // as 3d multivector
    lua["e31m_3d"] = e31m_3d;
    lua["e12m_3d"] = e12m_3d;
    lua["e23me_3d"] = e23me_3d; // as even grade 3d multivector
    lua["e31me_3d"] = e31me_3d;
    lua["e12me_3d"] = e12me_3d;

    lua["I_3d"] = I_3d; // as pscalar3d
    lua["e123_3d"] = e123_3d;
    lua["Im_3d"] = Im_3d;     // as 3d multivector
    lua["Im_3d_U"] = Im_3d_U; // as odd grade 3d multivector

    ////////////////////////////////////////////////////////////////////////////////
    // PGA 2dp constants
    ////////////////////////////////////////////////////////////////////////////////
    using namespace hd::ga::pga;

    // 2dp basis vectors
    lua["e1_2dp"] = e1_2dp;
    lua["e2_2dp"] = e2_2dp;
    lua["e3_2dp"] = e3_2dp;
    lua["origin_2dp"] = origin_2dp;

    // 2dp direction vectors
    lua["x_axis_direction_2dp"] = x_axis_direction_2dp;
    lua["y_axis_direction_2dp"] = y_axis_direction_2dp;

    // 2dp multivector versions
    lua["e1m_2dp"] = e1m_2dp;
    lua["e2m_2dp"] = e2m_2dp;
    lua["e3m_2dp"] = e3m_2dp;

    // 2dp bivectors
    lua["e23_2dp"] = e23_2dp;
    lua["e31_2dp"] = e31_2dp;
    lua["e12_2dp"] = e12_2dp;

    // 2dp axes as bivectors
    lua["x_axis_2dp"] = x_axis_2dp;
    lua["y_axis_2dp"] = y_axis_2dp;
    lua["horizon_2dp"] = horizon_2dp;

    ////////////////////////////////////////////////////////////////////////////////
    // PGA 3dp constants
    ////////////////////////////////////////////////////////////////////////////////

    // 3dp basis vectors
    lua["e1_3dp"] = e1_3dp;
    lua["e2_3dp"] = e2_3dp;
    lua["e3_3dp"] = e3_3dp;
    lua["e4_3dp"] = e4_3dp;
    lua["origin_3dp"] = origin_3dp;

    // 3dp direction vectors
    lua["x_axis_direction_3dp"] = x_axis_direction_3dp;
    lua["y_axis_direction_3dp"] = y_axis_direction_3dp;
    lua["z_axis_direction_3dp"] = z_axis_direction_3dp;

    // 3dp multivector versions
    lua["e1m_3dp"] = e1m_3dp;
    lua["e2m_3dp"] = e2m_3dp;
    lua["e3m_3dp"] = e3m_3dp;
    lua["e4m_3dp"] = e4m_3dp;

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
    lua["yz_plane_3dp"] = yz_plane_3dp;
    lua["zx_plane_3dp"] = zx_plane_3dp;
    lua["xy_plane_3dp"] = xy_plane_3dp;
    lua["horizon_3dp"] = horizon_3dp;

    // 3dp pseudoscalar
    lua["I_3dp"] = I_3dp;
    lua["e1234_3dp"] = e1234_3dp;
}