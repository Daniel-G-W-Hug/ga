#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_2d.hpp"

#include "ga_mvec2d_e.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// use MVec4_t including its ctors and add specific ctors for MVec4_t<T, Tag>
// by using partial template specialization for the Tag=mvec2d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec4_t<T, mvec2d_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec4_t(Scalar2d<T> s) : MVec4_t(T(s), T(0.0), T(0.0), T(0.0)) {}

    // assign a vector part exclusively (other grades = 0)
    MVec4_t(Vec2d<T> const& v) : MVec4_t(T(0.0), v.x, v.y, T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec4_t(PScalar2d<T> ps) : MVec4_t(T(0.0), T(0.0), T(0.0), T(ps)) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    MVec4_t(Scalar2d<T> s, PScalar2d<T> ps) : MVec4_t(T(s), T(0.0), T(0.0), T(ps)) {}

    // assign from a complex number, i.e. from the even subalgebra
    MVec4_t(MVec2d_E<T> const& M) : MVec4_t(M.c0, T(0.0), T(0.0), M.c1) {}

    // assign a full multivector
    MVec4_t(Scalar2d<T> s, Vec2d<T> const& v, PScalar2d<T> ps) :
        MVec4_t(T(s), v.x, v.y, T(ps))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec2d_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> gr0(MVec2d<T> const& M)
{
    return Scalar2d<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> gr1(MVec2d<T> const& M)
{
    return Vec2d<T>(M.c1, M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> gr2(MVec2d<T> const& M)
{
    return PScalar2d<T>(M.c3);
}

} // namespace hd::ga