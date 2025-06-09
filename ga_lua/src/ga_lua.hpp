#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#define SOL_ALL_SAFETIES_ON 1
#include "sol/sol.hpp"

#include "ga/ga_ega.hpp"

////////////////////////////////////////////////////////////////////////////////
// register basic types, geometric operations and constants
// of user-defined types with lua
////////////////////////////////////////////////////////////////////////////////
void register_2d_types(sol::state& lua);
void register_3d_types(sol::state& lua);
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
// functions commonly used across 2d and 3d types
////////////////////////////////////////////////////////////////////////////////
void register_functions(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::ega;

    lua.set_function("nrm_sq",
                     sol::overload(sol::resolve<value_t(scalar2d)>(nrm_sq),
                                   sol::resolve<value_t(vec2d const&)>(nrm_sq),
                                   sol::resolve<value_t(pscalar2d)>(nrm_sq),
                                   sol::resolve<value_t(mvec2d_e const&)>(nrm_sq),
                                   sol::resolve<value_t(mvec2d const&)>(nrm_sq),
                                   sol::resolve<value_t(scalar3d)>(nrm_sq),
                                   sol::resolve<value_t(vec3d const&)>(nrm_sq),
                                   sol::resolve<value_t(bivec3d const&)>(nrm_sq),
                                   sol::resolve<value_t(pscalar3d)>(nrm_sq),
                                   sol::resolve<value_t(mvec3d_e const&)>(nrm_sq),
                                   sol::resolve<value_t(mvec3d_u const&)>(nrm_sq),
                                   sol::resolve<value_t(mvec3d const&)>(nrm_sq)));

    lua.set_function("nrm", sol::overload(sol::resolve<value_t(scalar2d)>(nrm),
                                          sol::resolve<value_t(vec2d const&)>(nrm),
                                          sol::resolve<value_t(pscalar2d)>(nrm),
                                          sol::resolve<value_t(mvec2d_e const&)>(nrm),
                                          sol::resolve<value_t(mvec2d const&)>(nrm),
                                          sol::resolve<value_t(scalar3d)>(nrm),
                                          sol::resolve<value_t(vec3d const&)>(nrm),
                                          sol::resolve<value_t(bivec3d const&)>(nrm),
                                          sol::resolve<value_t(pscalar3d)>(nrm),
                                          sol::resolve<value_t(mvec3d_e const&)>(nrm),
                                          sol::resolve<value_t(mvec3d_u const&)>(nrm),
                                          sol::resolve<value_t(mvec3d const&)>(nrm)));

    lua.set_function("normalize",
                     sol::overload(sol::resolve<scalar2d(scalar2d)>(normalize),
                                   sol::resolve<vec2d(vec2d const&)>(normalize),
                                   sol::resolve<pscalar2d(pscalar2d)>(normalize),
                                   sol::resolve<scalar3d(scalar3d)>(normalize),
                                   sol::resolve<vec3d(vec3d const&)>(normalize),
                                   sol::resolve<bivec3d(bivec3d const&)>(normalize),
                                   sol::resolve<pscalar3d(pscalar3d)>(normalize)));


    ////////////////////////////////////////////////////////////////////////////////
    // grade operations on multivectors
    ////////////////////////////////////////////////////////////////////////////////
    lua.set_function("gr0", sol::overload(sol::resolve<scalar2d(mvec2d_e const&)>(gr0),
                                          sol::resolve<scalar2d(mvec2d const&)>(gr0),
                                          sol::resolve<scalar3d(mvec3d_e const&)>(gr0),
                                          sol::resolve<scalar3d(mvec3d const&)>(gr0)));

    lua.set_function("gr1", sol::overload(sol::resolve<vec2d(mvec2d const&)>(gr1),
                                          sol::resolve<vec3d(mvec3d_u const&)>(gr1),
                                          sol::resolve<vec3d(mvec3d const&)>(gr1)));

    lua.set_function("gr2", sol::overload(sol::resolve<pscalar2d(mvec2d_e const&)>(gr2),
                                          sol::resolve<pscalar2d(mvec2d const&)>(gr2),
                                          sol::resolve<bivec3d(mvec3d_e const&)>(gr2),
                                          sol::resolve<bivec3d(mvec3d const&)>(gr2)));

    lua.set_function("gr3", sol::overload(sol::resolve<pscalar3d(mvec3d_u const&)>(gr3),
                                          sol::resolve<pscalar3d(mvec3d const&)>(gr3)));


    ////////////////////////////////////////////////////////////////////////////////
    // grade operations on multivectors (gr_inv, rev, inv)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("gr_inv",
                     sol::overload(sol::resolve<scalar2d(scalar2d)>(gr_inv),
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
                                   sol::resolve<mvec3d(mvec3d const&)>(gr_inv)));

    lua.set_function("rev", sol::overload(sol::resolve<scalar2d(scalar2d)>(rev),
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
                                          sol::resolve<mvec3d(mvec3d const&)>(rev)));

    lua.set_function("conj", sol::overload(sol::resolve<scalar2d(scalar2d)>(conj),
                                           sol::resolve<vec2d(vec2d const&)>(conj),
                                           sol::resolve<pscalar2d(pscalar2d)>(conj),
                                           sol::resolve<mvec2d_e(mvec2d_e const&)>(conj),
                                           sol::resolve<mvec2d(mvec2d const&)>(conj),
                                           sol::resolve<scalar3d(scalar3d)>(conj),
                                           sol::resolve<vec3d(vec3d const&)>(conj),
                                           sol::resolve<pscalar3d(pscalar3d)>(conj),
                                           sol::resolve<mvec3d_e(mvec3d_e const&)>(conj),
                                           sol::resolve<mvec3d_u(mvec3d_u const&)>(conj),
                                           sol::resolve<mvec3d(mvec3d const&)>(conj)));

    lua.set_function("inv", sol::overload(sol::resolve<scalar2d(scalar2d)>(inv),
                                          sol::resolve<vec2d(vec2d const&)>(inv),
                                          sol::resolve<pscalar2d(pscalar2d)>(inv),
                                          sol::resolve<mvec2d_e(mvec2d_e const&)>(inv),
                                          sol::resolve<mvec2d(mvec2d const&)>(inv),
                                          sol::resolve<scalar3d(scalar3d)>(inv),
                                          sol::resolve<vec3d(vec3d const&)>(inv),
                                          sol::resolve<bivec3d(bivec3d const&)>(inv),
                                          sol::resolve<pscalar3d(pscalar3d)>(inv),
                                          sol::resolve<mvec3d_e(mvec3d_e const&)>(inv),
                                          sol::resolve<mvec3d(mvec3d const&)>(inv)));


    ////////////////////////////////////////////////////////////////////////////////
    // scalar product (dot)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "dot", sol::overload(sol::resolve<scalar2d(scalar2d, scalar2d)>(dot),
                             sol::resolve<scalar2d(vec2d const&, vec2d const&)>(dot),
                             sol::resolve<scalar2d(pscalar2d, pscalar2d)>(dot),
                             sol::resolve<scalar2d(mvec2d const&, mvec2d const&)>(dot),
                             // now left and right contraction in 2d:
                             //  sol::resolve<vec2d(vec2d const&, pscalar2d)>(dot),
                             //  sol::resolve<vec2d(pscalar2d, vec2d const&)>(dot),
                             sol::resolve<scalar3d(scalar3d, scalar3d)>(dot),
                             sol::resolve<scalar3d(vec3d const&, vec3d const&)>(dot),
                             sol::resolve<scalar3d(bivec3d const&, bivec3d const&)>(dot),
                             sol::resolve<scalar3d(pscalar3d, pscalar3d)>(dot),
                             sol::resolve<scalar3d(mvec3d const&, mvec3d const&)>(dot)
                             // now left and right contraction in 3d:
                             // , sol::resolve<vec3d(vec3d const&, bivec3d const&)>(dot),
                             // sol::resolve<vec3d(bivec3d const&, vec3d const&)>(dot)
                             ));


    ////////////////////////////////////////////////////////////////////////////////
    // wedge product (= outer product) (wdg)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "wdg", sol::overload(sol::resolve<scalar2d(scalar2d, scalar2d)>(wdg),
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
                             sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(wdg)));


    ////////////////////////////////////////////////////////////////////////////////
    // regressive wedge product (= outer product for complements) (rwdg)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "rwdg", sol::overload(sol::resolve<scalar2d(vec2d const&, vec2d const&)>(rwdg),
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
                              sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(rwdg)));

    ////////////////////////////////////////////////////////////////////////////////
    // other products (inner, fat_dot, cmt, cross)
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("cmt", sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(cmt));

    lua.set_function("cross", sol::resolve<vec3d(vec3d const&, vec3d const&)>(cross));


    ////////////////////////////////////////////////////////////////////////////////
    // angles and rotations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "angle",
        sol::overload(sol::resolve<value_t(vec2d const&, vec2d const&)>(angle),
                      sol::resolve<value_t(vec3d const&, vec3d const&)>(angle),
                      sol::resolve<value_t(bivec3d const&, bivec3d const&)>(angle),
                      sol::resolve<value_t(vec3d const&, bivec3d const&)>(angle),
                      sol::resolve<value_t(bivec3d const&, vec3d const&)>(angle)));

    lua.set_function("angle_to_re", sol::resolve<value_t(mvec2d_e const&)>(angle_to_re));

    lua.set_function("exp",
                     sol::overload(sol::resolve<mvec2d_e(pscalar2d, value_t)>(exp),
                                   sol::resolve<mvec3d_e(bivec3d const&, value_t)>(exp)));

    lua.set_function(
        "rotor", sol::overload(sol::resolve<mvec2d_e(pscalar2d, value_t)>(rotor),
                               sol::resolve<mvec3d_e(bivec3d const&, value_t)>(rotor)));

    lua.set_function(
        "rotate",
        sol::overload(sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(rotate),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(rotate),
                      sol::resolve<vec3d(vec3d const&, mvec3d_e const&)>(rotate),
                      sol::resolve<bivec3d(bivec3d const&, mvec3d_e const&)>(rotate),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d_e const&)>(rotate)));

    ////////////////////////////////////////////////////////////////////////////////
    // dualization operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("dual", sol::overload(sol::resolve<scalar2d(pscalar2d)>(dual),
                                           sol::resolve<pscalar2d(scalar2d)>(dual),
                                           sol::resolve<vec2d(vec2d const&)>(dual),
                                           sol::resolve<mvec2d_e(mvec2d_e const&)>(dual),
                                           sol::resolve<mvec2d(mvec2d const&)>(dual),
                                           sol::resolve<scalar3d(pscalar3d)>(dual),
                                           sol::resolve<pscalar3d(scalar3d)>(dual),
                                           sol::resolve<bivec3d(vec3d const&)>(dual),
                                           sol::resolve<vec3d(bivec3d const&)>(dual),
                                           sol::resolve<mvec3d_u(mvec3d_e const&)>(dual),
                                           sol::resolve<mvec3d_e(mvec3d_u const&)>(dual),
                                           sol::resolve<mvec3d(mvec3d const&)>(dual)));

    ////////////////////////////////////////////////////////////////////////////////
    // complement operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("rcmpl",
                     sol::overload(sol::resolve<pscalar2d(scalar2d)>(rcmpl),
                                   sol::resolve<vec2d(vec2d const&)>(rcmpl),
                                   sol::resolve<scalar2d(pscalar2d)>(rcmpl),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(rcmpl),
                                   sol::resolve<mvec2d(mvec2d const&)>(rcmpl)));

    lua.set_function("lcmpl",
                     sol::overload(sol::resolve<pscalar2d(scalar2d)>(lcmpl),
                                   sol::resolve<vec2d(vec2d const&)>(lcmpl),
                                   sol::resolve<scalar2d(pscalar2d)>(lcmpl),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(lcmpl),
                                   sol::resolve<mvec2d(mvec2d const&)>(lcmpl)));

    lua.set_function("cmpl", sol::overload(sol::resolve<pscalar3d(scalar3d)>(cmpl),
                                           sol::resolve<bivec3d(vec3d const&)>(cmpl),
                                           sol::resolve<vec3d(bivec3d const&)>(cmpl),
                                           sol::resolve<scalar3d(pscalar3d)>(cmpl),
                                           sol::resolve<mvec3d_u(mvec3d_e const&)>(cmpl),
                                           sol::resolve<mvec3d_e(mvec3d_u const&)>(cmpl),
                                           sol::resolve<mvec3d(mvec3d const&)>(cmpl)));

    ////////////////////////////////////////////////////////////////////////////////
    // projections, rejections and reflections
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "project_onto",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(project_onto),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(project_onto),
                      sol::resolve<vec3d(vec3d const&, bivec3d const&)>(project_onto)));

    lua.set_function(
        "reject_from",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(reject_from),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(reject_from),
                      sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reject_from)));

    lua.set_function(
        "reflect_on",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on),
                      sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reflect_on),
                      sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(reflect_on)));

    lua.set_function(
        "reflect_on_vec",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on_vec),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on_vec)));

    ////////////////////////////////////////////////////////////////////////////////
    // common helper functions for scripting in lua
    ////////////////////////////////////////////////////////////////////////////////

    // convert scalars & pscalars into numeric values for further calculations
    lua.set_function("to_val", sol::overload(sol::resolve<value_t(scalar2d)>(to_val),
                                             sol::resolve<value_t(scalar3d)>(to_val),
                                             sol::resolve<value_t(pscalar2d)>(to_val),
                                             sol::resolve<value_t(pscalar3d)>(to_val)));

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
    lua["Im_3d_U"] = Im_3d_U; // as uneven grade 3d multivector
}