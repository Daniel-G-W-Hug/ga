#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_2dp.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec2dp_U<T> M = (c0 * e1 + c1 * e2 + c2 * e3) + c4 * e3^e2^e1
//
// with the term in brackets being the vector c0 * e1 + c1 * e2 + c2 * e3
// and the term c4 e1^e2^e3 being the trivector
//
// This is the definition of a multivector in the uneven subalgebra of G<2,0,1>,
// i.e. it models only multivectors with even grades 1 and 3
// This subalgebra is used to store intermediate results for rotations.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.

////////////////////////////////////////////////////////////////////////////////
// use MVec4_t including its ctors and add specific ctors for MVec4_t<T, Tag>
// by using partial template specialization for the Tag=mvec2dp_u_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec4_t<T, mvec2dp_u_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a vector part exclusively (other grades = 0)
    MVec4_t(Vec2dp<T> const& v) : MVec4_t(v.x, v.y, v.z, T(0.0)) {}

    // assign a scalar part exclusively (other grades = 0)
    MVec4_t(PScalar2dp<T> ps) : MVec4_t(T(0.0), T(0.0), T(0.0), T(ps)) {}

    // assign vector and pseudoscalar parts
    MVec4_t(Vec2dp<T> const& v, PScalar2dp<T> ps) : MVec4_t(v.x, v.y, v.z, T(ps)) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec2dp_u_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of the uneven multivector
//
// grade 1: gr1() - vector
// grade 2: gr3() - pseudoscalar

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> gr1(MVec2dp_U<T> const& M)
{
    return Vec2dp<T>(M.c0, M.c1, M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> gr3(MVec2dp_U<T> const& M)
{
    return PScalar2dp<T>(M.c3);
}

} // namespace hd::ga