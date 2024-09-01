#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_mvec4_t.hpp"

#include "type_t/ga_type_0d.hpp"
#include "type_t/ga_type_2d.hpp"

#include "ga_mvec2d_e.hpp"


namespace hd::ga {

template <typename T> using MVec2d = MVec4_t<T, mvec2d_tag>;

////////////////////////////////////////////////////////////////////////////////
// use MVec4_t including its ctors and add specific ctors for MVec4_t<T, Tag>
// by using partial template specialization for the Tag=mvec2d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec4_t<T, mvec2d_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec4_t(Scalar<T> s) : MVec4_t(T(s), T(0.0), T(0.0), T(0.0)) {}

    // assign a vector part exclusively (other grades = 0)
    MVec4_t(Vec2d<T> const& v) : MVec4_t(T(0.0), v.x, v.y, T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec4_t(PScalar2d<T> ps) : MVec4_t(T(0.0), T(0.0), T(0.0), T(ps)) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    MVec4_t(Scalar<T> s, PScalar2d<T> ps) : MVec4_t(T(s), T(0.0), T(0.0), T(ps)) {}

    // assign from a complex number, i.e. from the even subalgebra
    MVec4_t(MVec2d_E<T> v) : MVec4_t(v.c0, T(0.0), T(0.0), v.c1) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec2d_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T> inline constexpr Scalar<T> gr0(MVec2d<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr Vec2d<T> gr1(MVec2d<T> const& v)
{
    return Vec2d<T>(v.c1, v.c2);
}

template <typename T> inline constexpr PScalar2d<T> gr2(MVec2d<T> const& v)
{
    return PScalar2d<T>(v.c3);
}

} // namespace hd::ga