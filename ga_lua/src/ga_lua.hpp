#pragma once

// author: Daniel Hug, 2024

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

    lua.new_usertype<scalar>(
        "scalar",
        sol::constructors<scalar(), scalar(value_t const&), scalar(value_t&&)>());


    lua.new_usertype<pscalar2d>(
        "pscalar2d",
        sol::constructors<pscalar2d(), pscalar2d(value_t const&), pscalar2d(value_t&&)>(),
        sol::meta_function::unary_minus, sol::resolve<pscalar2d(pscalar2d)>(operator-),
        sol::meta_function::addition,
        sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<pscalar2d(pscalar2d, value_t)>(operator*),
                      sol::resolve<pscalar2d(value_t, pscalar2d)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<mvec2d(pscalar2d, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d_e(pscalar2d, mvec2d_e const&)>(operator*),
                      sol::resolve<vec2d(pscalar2d, vec2d const&)>(operator*),
                      sol::resolve<value_t(pscalar2d, pscalar2d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar2d(pscalar2d, value_t)>(operator/));


    lua.new_usertype<vec2d>(
        "vec2d", sol::constructors<vec2d(), vec2d(value_t, value_t), vec2d(vec2d)>(), "x",
        &vec2d::x, "y", &vec2d::y, sol::meta_function::unary_minus,
        sol::resolve<vec2d(vec2d const&)>(operator-), sol::meta_function::addition,
        sol::resolve<vec2d(vec2d const&, vec2d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<vec2d(vec2d const&, vec2d const&)>(operator-),
        sol::meta_function::power_of,
        sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(wdg),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec2d(vec2d const&, value_t)>(operator*),
                      sol::resolve<vec2d(value_t, vec2d const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(operator*),
                      sol::resolve<mvec2d(vec2d const&, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d_e(vec2d const&, vec2d const&)>(operator*),
                      sol::resolve<vec2d(vec2d const&, pscalar2d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec2d(vec2d const&, value_t)>(operator/));


    lua.new_usertype<mvec2d_e>(
        "mvec2d_e",
        sol::constructors<mvec2d_e(), mvec2d_e(value_t, value_t), mvec2d_e(scalar),
                          mvec2d_e(pscalar2d), mvec2d_e(scalar, pscalar2d),
                          mvec2d_e(mvec2d_e)>(),
        "c0", &mvec2d_e::c0, "c1", &mvec2d_e::c1, sol::meta_function::unary_minus,
        sol::resolve<mvec2d_e(mvec2d_e const&)>(operator-), sol::meta_function::addition,
        sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator-),
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
        sol::resolve<mvec2d_e(mvec2d_e const&, value_t)>(operator/));


    lua.new_usertype<mvec2d>(
        "mvec2d",
        sol::constructors<mvec2d(), mvec2d(value_t, value_t, value_t, value_t),
                          mvec2d(scalar), mvec2d(vec2d), mvec2d(pscalar2d),
                          mvec2d(scalar, pscalar2d), mvec2d(mvec2d_e), mvec2d(mvec2d)>(),
        "c0", &mvec2d::c0, "c1", &mvec2d::c1, "c2", &mvec2d::c2, "c3", &mvec2d::c3,
        sol::meta_function::unary_minus, sol::resolve<mvec2d(mvec2d const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec2d(mvec2d const&, value_t)>(operator*),
                      sol::resolve<mvec2d(value_t, mvec2d const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, pscalar2d)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2d(mvec2d const&, value_t)>(operator/));
}


////////////////////////////////////////////////////////////////////////////////
// 3d types
////////////////////////////////////////////////////////////////////////////////
void register_3d_types(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::ega;

    lua.new_usertype<pscalar3d>(
        "pscalar3d",
        sol::constructors<pscalar3d(), pscalar3d(value_t const&), pscalar3d(value_t&&)>(),
        sol::meta_function::unary_minus, sol::resolve<pscalar3d(pscalar3d)>(operator-),
        sol::meta_function::addition,
        sol::resolve<pscalar3d(pscalar3d, pscalar3d)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<pscalar3d(pscalar3d, pscalar3d)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<pscalar3d(pscalar3d, value_t)>(operator*),
                      sol::resolve<pscalar3d(value_t, pscalar3d)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<bivec3d(pscalar3d, vec3d const&)>(operator*),
                      sol::resolve<vec3d(pscalar3d, bivec3d const&)>(operator*),
                      sol::resolve<mvec3d_e(pscalar3d, mvec3d_u const&)>(operator*),
                      sol::resolve<mvec3d_u(pscalar3d, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d(pscalar3d, mvec3d const&)>(operator*),
                      sol::resolve<value_t(pscalar3d, pscalar3d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar3d(pscalar3d, value_t)>(operator/));


    lua.new_usertype<vec3d>(
        "vec3d",
        sol::constructors<vec3d(), vec3d(value_t, value_t, value_t), vec3d(vec3d)>(), "x",
        &vec3d::x, "y", &vec3d::x, "z", &vec3d::y, sol::meta_function::unary_minus,
        sol::resolve<vec3d(vec3d const&)>(operator-), sol::meta_function::addition,
        sol::resolve<vec3d(vec3d const&, vec3d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<vec3d(vec3d const&, vec3d const&)>(operator-),
        sol::meta_function::power_of,
        sol::overload(sol::resolve<bivec3d(vec3d const&, vec3d const&)>(wdg),
                      sol::resolve<pscalar3d(vec3d const&, bivec3d const&)>(wdg)),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec3d(vec3d const&, value_t)>(operator*),
                      sol::resolve<vec3d(value_t, vec3d const&)>(operator*),
                      sol::resolve<mvec3d_e(vec3d const&, vec3d const&)>(operator*),
                      sol::resolve<mvec3d_u(vec3d const&, bivec3d const&)>(operator*),
                      sol::resolve<bivec3d(vec3d const&, pscalar3d)>(operator*),
                      sol::resolve<mvec3d_u(vec3d const&, mvec3d_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec3d(vec3d const&, value_t)>(operator/));


    lua.new_usertype<bivec3d>(
        "bivec3d",
        sol::constructors<bivec3d(), bivec3d(value_t, value_t, value_t),
                          bivec3d(bivec3d)>(),
        "x", &bivec3d::x, "y", &bivec3d::x, "z", &bivec3d::y,
        sol::meta_function::unary_minus, sol::resolve<bivec3d(bivec3d const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(operator-),
        sol::meta_function::power_of,
        sol::resolve<pscalar3d(bivec3d const&, vec3d const&)>(wdg),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<bivec3d(bivec3d const&, value_t)>(operator*),
                      sol::resolve<bivec3d(value_t, bivec3d const&)>(operator*),
                      sol::resolve<mvec3d_u(bivec3d const&, vec3d const&)>(operator*),
                      sol::resolve<mvec3d_e(bivec3d const&, bivec3d const&)>(operator*),
                      sol::resolve<vec3d(bivec3d const&, pscalar3d)>(operator*),
                      sol::resolve<mvec3d_e(bivec3d const&, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d_u(bivec3d const&, mvec3d_u const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<bivec3d(bivec3d const&, value_t)>(operator/));


    lua.new_usertype<mvec3d_e>(
        "mvec3d_e",
        sol::constructors<mvec3d_e(), mvec3d_e(value_t, value_t, value_t, value_t),
                          mvec3d_e(scalar), mvec3d_e(bivec3d), mvec3d_e(scalar, bivec3d),
                          mvec3d_e(mvec3d_e)>(),
        "c0", &mvec3d_e::c0, "c1", &mvec3d_e::c1, "c2", &mvec3d_e::c2, "c3",
        &mvec3d_e::c3, sol::meta_function::unary_minus,
        sol::resolve<mvec3d_e(mvec3d_e const&)>(operator-), sol::meta_function::addition,
        sol::resolve<mvec3d_e(mvec3d_e const&, mvec3d_e const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec3d_e(mvec3d_e const&, mvec3d_e const&)>(operator-),
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
        sol::resolve<mvec3d_e(mvec3d_e const&, value_t)>(operator/));


    lua.new_usertype<mvec3d_u>(
        "mvec3d_u",
        sol::constructors<mvec3d_u(), mvec3d_u(value_t, value_t, value_t, value_t),
                          mvec3d_u(pscalar3d), mvec3d_u(vec3d),
                          mvec3d_u(vec3d, pscalar3d), mvec3d_u(mvec3d_u)>(),
        "c0", &mvec3d_u::c0, "c1", &mvec3d_u::c1, "c2", &mvec3d_u::c2, "c3",
        &mvec3d_u::c3, sol::meta_function::unary_minus,
        sol::resolve<mvec3d_u(mvec3d_u const&)>(operator-), sol::meta_function::addition,
        sol::resolve<mvec3d_u(mvec3d_u const&, mvec3d_u const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec3d_u(mvec3d_u const&, mvec3d_u const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec3d_u(mvec3d_u const&, value_t)>(operator*),
                      sol::resolve<mvec3d_u(value_t, mvec3d_u const&)>(operator*),
                      sol::resolve<mvec3d_e(mvec3d_u const&, pscalar3d)>(operator*),
                      sol::resolve<mvec3d_e(mvec3d_u const&, mvec3d_u const&)>(operator*),
                      sol::resolve<mvec3d_u(mvec3d_u const&, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d_u(mvec3d_u const&, bivec3d const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3d_u(mvec3d_u const&, value_t)>(operator/));


    lua.new_usertype<mvec3d>(
        "mvec3d",
        sol::constructors<mvec3d(),
                          mvec3d(value_t, value_t, value_t, value_t, value_t, value_t,
                                 value_t, value_t),
                          mvec3d(scalar), mvec3d(vec3d), mvec3d(bivec3d),
                          mvec3d(scalar, bivec3d), mvec3d(mvec3d_e), mvec3d(mvec3d_u),
                          mvec3d(vec3d, pscalar3d), mvec3d(pscalar3d), mvec3d(mvec3d)>(),
        "c0", &mvec3d::c0, "c1", &mvec3d::c1, "c2", &mvec3d::c2, "c3", &mvec3d::c3, "c4",
        &mvec3d::c4, "c5", &mvec3d::c5, "c6", &mvec3d::c6, "c7", &mvec3d::c7,
        sol::meta_function::unary_minus, sol::resolve<mvec3d(mvec3d const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec3d(mvec3d const&, value_t)>(operator*),
                      sol::resolve<mvec3d(value_t, mvec3d const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, mvec3d_e const&)>(operator*),
                      sol::resolve<mvec3d(mvec3d const&, pscalar3d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec3d(mvec3d const&, value_t)>(operator/));
}


////////////////////////////////////////////////////////////////////////////////
// functions commonly used across 2d and 3d types
////////////////////////////////////////////////////////////////////////////////
void register_functions(sol::state& lua)
{
    using namespace hd::ga;
    using namespace hd::ga::ega;

    lua.set_function(
        "dot", sol::overload(sol::resolve<value_t(vec2d const&, vec2d const&)>(dot),
                             sol::resolve<vec2d(pscalar2d, vec2d const&)>(dot),
                             sol::resolve<vec2d(vec2d const&, pscalar2d)>(dot),
                             sol::resolve<value_t(vec3d const&, vec3d const&)>(dot),
                             sol::resolve<value_t(bivec3d const&, bivec3d const&)>(dot),
                             sol::resolve<vec3d(bivec3d const&, vec3d const&)>(dot),
                             sol::resolve<vec3d(vec3d const&, bivec3d const&)>(dot)));

    lua.set_function("cmt", sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(cmt));

    lua.set_function("sq_nrm",
                     sol::overload(sol::resolve<value_t(vec2d const&)>(sq_nrm),
                                   sol::resolve<value_t(mvec2d_e const&)>(sq_nrm),
                                   sol::resolve<value_t(mvec2d const&)>(sq_nrm),
                                   sol::resolve<value_t(pscalar2d)>(sq_nrm),
                                   sol::resolve<value_t(vec3d const&)>(sq_nrm),
                                   sol::resolve<value_t(bivec3d const&)>(sq_nrm),
                                   sol::resolve<value_t(mvec3d_e const&)>(sq_nrm),
                                   sol::resolve<value_t(mvec3d const&)>(sq_nrm),
                                   sol::resolve<value_t(pscalar3d)>(sq_nrm)));

    lua.set_function("nrm", sol::overload(sol::resolve<value_t(vec2d const&)>(nrm),
                                          sol::resolve<value_t(mvec2d_e const&)>(nrm),
                                          sol::resolve<value_t(mvec2d const&)>(nrm),
                                          sol::resolve<value_t(pscalar2d)>(nrm),
                                          sol::resolve<value_t(vec3d const&)>(nrm),
                                          sol::resolve<value_t(bivec3d const&)>(nrm),
                                          sol::resolve<value_t(mvec3d_e const&)>(nrm),
                                          sol::resolve<value_t(mvec3d const&)>(nrm),
                                          sol::resolve<value_t(pscalar3d)>(nrm)));

    lua.set_function("unitized",
                     sol::overload(sol::resolve<vec2d(vec2d const&)>(unitized),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(unitized),
                                   sol::resolve<mvec2d(mvec2d const&)>(unitized),
                                   sol::resolve<vec3d(vec3d const&)>(unitized),
                                   sol::resolve<bivec3d(bivec3d const&)>(unitized),
                                   sol::resolve<mvec3d_e(mvec3d_e const&)>(unitized),
                                   sol::resolve<mvec3d(mvec3d const&)>(unitized)));

    lua.set_function("inv", sol::overload(sol::resolve<vec2d(vec2d const&)>(inv),
                                          sol::resolve<mvec2d_e(mvec2d_e const&)>(inv),
                                          sol::resolve<mvec2d(mvec2d const&)>(inv),
                                          sol::resolve<pscalar2d(pscalar2d)>(inv),
                                          sol::resolve<vec3d(vec3d const&)>(inv),
                                          sol::resolve<bivec3d(bivec3d const&)>(inv),
                                          sol::resolve<mvec3d(mvec3d const&)>(inv),
                                          sol::resolve<pscalar3d(pscalar3d)>(inv)));

    lua.set_function(
        "wdg", sol::overload(sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(wdg),
                             sol::resolve<bivec3d(vec3d const&, vec3d const&)>(wdg),
                             sol::resolve<pscalar3d(vec3d const&, bivec3d const&)>(wdg),
                             sol::resolve<pscalar3d(bivec3d const&, vec3d const&)>(wdg)));

    lua.set_function("cross", sol::resolve<vec3d(vec3d const&, vec3d const&)>(cross));

    lua.set_function("rev", sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&)>(rev),
                                          sol::resolve<mvec2d(mvec2d const&)>(rev),
                                          sol::resolve<pscalar2d(pscalar2d)>(rev),
                                          sol::resolve<bivec3d(bivec3d const&)>(rev),
                                          sol::resolve<mvec3d(mvec3d const&)>(rev),
                                          sol::resolve<mvec3d_e(mvec3d_e const&)>(rev),
                                          sol::resolve<mvec3d_u(mvec3d_u const&)>(rev),
                                          sol::resolve<pscalar3d(pscalar3d)>(rev)));

    lua.set_function("conj", sol::overload(sol::resolve<mvec2d(mvec2d const&)>(conj),
                                           sol::resolve<mvec3d(mvec3d const&)>(conj)));

    ////////////////////////////////////////////////////////////////////////////////
    // projections, rejections and reflections
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function(
        "project_onto",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(project_onto),
                      sol::resolve<vec2d(vec2d const&, pscalar2d)>(project_onto),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(project_onto),
                      sol::resolve<vec3d(vec3d const&, bivec3d const&)>(project_onto)));

    lua.set_function(
        "project_onto_unitized",
        sol::overload(
            sol::resolve<vec2d(vec2d const&, vec2d const&)>(project_onto_unitized),
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(project_onto_unitized),
            sol::resolve<vec3d(vec3d const&, bivec3d const&)>(project_onto_unitized)));

    lua.set_function(
        "reject_from",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(reject_from),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(reject_from),
                      sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reject_from)));

    lua.set_function(
        "reject_from_unitized",
        sol::overload(
            sol::resolve<vec2d(vec2d const&, vec2d const&)>(reject_from_unitized),
            sol::resolve<vec3d(vec3d const&, vec3d const&)>(reject_from_unitized),
            sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reject_from_unitized)));

    lua.set_function(
        "reflect_on_hyp",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on_hyp),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on_hyp)));

    lua.set_function(
        "reflect_on_vec",
        sol::overload(sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on_vec),
                      sol::resolve<vec3d(vec3d const&, vec3d const&)>(reflect_on_vec)));

    lua.set_function(
        "reflect_on",
        sol::overload(sol::resolve<vec3d(vec3d const&, bivec3d const&)>(reflect_on),
                      sol::resolve<bivec3d(bivec3d const&, bivec3d const&)>(reflect_on)));

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

    lua.set_function("dual2d",
                     sol::overload(sol::resolve<scalar(pscalar2d)>(dual2d),
                                   sol::resolve<pscalar2d(scalar)>(dual2d),
                                   sol::resolve<vec2d(vec2d const&)>(dual2d),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(dual2d),
                                   sol::resolve<mvec2d(mvec2d const&)>(dual2d)));

    lua.set_function("dual3d",
                     sol::overload(sol::resolve<scalar(pscalar3d)>(dual3d),
                                   sol::resolve<pscalar3d(scalar)>(dual3d),
                                   sol::resolve<bivec3d(vec3d const&)>(dual3d),
                                   sol::resolve<vec3d(bivec3d const&)>(dual3d),
                                   sol::resolve<mvec3d_u(mvec3d_e const&)>(dual3d),
                                   sol::resolve<mvec3d_e(mvec3d_u const&)>(dual3d),
                                   sol::resolve<mvec3d(mvec3d const&)>(dual3d)));

    ////////////////////////////////////////////////////////////////////////////////
    // grade operations on multivectors
    ////////////////////////////////////////////////////////////////////////////////
    lua.set_function("gr0", sol::overload(sol::resolve<scalar(mvec2d_e const&)>(gr0),
                                          sol::resolve<scalar(mvec2d const&)>(gr0),
                                          sol::resolve<scalar(mvec3d const&)>(gr0),
                                          sol::resolve<scalar(mvec3d_e const&)>(gr0)));

    lua.set_function("gr1", sol::overload(sol::resolve<vec2d(mvec2d const&)>(gr1),
                                          sol::resolve<vec3d(mvec3d const&)>(gr1),
                                          sol::resolve<vec3d(mvec3d_u const&)>(gr1)));

    lua.set_function("gr2", sol::overload(sol::resolve<pscalar2d(mvec2d_e const&)>(gr2),
                                          sol::resolve<pscalar2d(mvec2d const&)>(gr2),
                                          sol::resolve<bivec3d(mvec3d const&)>(gr2),
                                          sol::resolve<bivec3d(mvec3d_e const&)>(gr2)));

    lua.set_function("gr3", sol::overload(sol::resolve<pscalar3d(mvec3d const&)>(gr3),
                                          sol::resolve<pscalar3d(mvec3d_u const&)>(gr3)));

    ////////////////////////////////////////////////////////////////////////////////
    // common helper functions for scripting in lua
    ////////////////////////////////////////////////////////////////////////////////

    // convert scalars & pscalars into numeric values for further calculations
    lua.set_function("to_val", sol::overload(sol::resolve<value_t(scalar)>(to_val),
                                             sol::resolve<value_t(pscalar2d)>(to_val),
                                             sol::resolve<value_t(pscalar3d)>(to_val)));

    lua.set_function("rad_to_deg", &rad_to_deg);
    lua.set_function("deg_to_rad", &deg_to_rad);
}

////////////////////////////////////////////////////////////////////////////////
// make defined constants available as global variables in lua
////////////////////////////////////////////////////////////////////////////////
void register_constants(sol::state& lua)
{
    // TODO: manipulate global lua metatable such that these entries
    //       cannot be modified by the user
    //       (currently these values can be changed)

    using namespace hd::ga;
    using namespace hd::ga::ega;

    // general constants
    lua["eps"] = eps;

    // 2d constants
    lua["e1_2d"] = e1_2d; // as 2d vector
    lua["e2_2d"] = e2_2d;
    lua["e1m_2d"] = e1m_2d; // as 2d multivector
    lua["e2m_2d"] = e2m_2d;

    lua["I_2d"] = I_2d;       // as pscalar2d
    lua["Im_2d"] = Im_2d;     // as 2d multivector
    lua["Im_2d_E"] = Im_2d_E; // as even grade 2d multivector

    // 3d constants
    lua["e1_3d"] = e1_3d; // as 3d vector
    lua["e2_3d"] = e2_3d;
    lua["e3_3d"] = e3_3d;
    lua["e1m_3d"] = e1m_3d; // as 3d multivector
    lua["e2m_3d"] = e2m_3d;
    lua["e3m_3d"] = e3m_3d;

    lua["e23_3d"] = e23_3d; // as 3d bivector
    lua["e31_3d"] = e31_3d;
    lua["e12_3d"] = e12_3d;
    lua["e23m_3d"] = e23m_3d; // as 3d multivector
    lua["e31m_3d"] = e31m_3d;
    lua["e12m_3d"] = e12m_3d;
    lua["e23me_3d"] = e23me_3d; // as even grade 3d multivector
    lua["e31me_3d"] = e31me_3d;
    lua["e12me_3d"] = e12me_3d;

    lua["I_3d"] = I_3d;       // as pscalar3d
    lua["Im_3d"] = Im_3d;     // as 3d multivector
    lua["Im_3d_U"] = Im_3d_U; // as uneven grade 3d multivector
}