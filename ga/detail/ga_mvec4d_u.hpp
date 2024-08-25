#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_mvec8_t.hpp"

#include "type_t/ga_type_0d.hpp"
#include "type_t/ga_type_4d.hpp"


namespace hd::ga {

template <typename T> using MVec4d_U = MVec8_t<T, mvec4d_u_tag>;

////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec4d_u_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec8_t(MVec8_t<U, mvec4d_u_tag> const& v) :
        MVec8_t(v.c0, v.c1, v.c2, v.c3, v.c4, v.c5, v.c6, v.c7)
    {
    }

    // assign the vector part exclusively (other grades = 0)
    MVec8_t(Vec4d<T> v) : MVec8_t(v.x, v.y, v.z, v.w, 0.0, 0.0, 0.0, 0.0) {}

    // assign the trivector part exclusively (other grades = 0)
    MVec8_t(TriVec4d<T> const& t) : MVec8_t(0.0, 0.0, 0.0, 0.0, t.x, t.y, t.z, t.w) {}

    // assign from a vector and a trivector directly
    MVec8_t(Vec4d<T> const& v, TriVec4d<T> const& t) :
        MVec8_t(v.x, v.y, v.z, v.w, t.x, t.y, t.z, t.w)
    {
    }

    // assign from the uneven subalgebra directly
    MVec8_t(MVec4d_U<T> const& v) :
        MVec8_t(v.c0, v.c1, v.c2, v.c3, v.c4, v.c5, v.c6, v.c7)
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec4d_u_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a 4d uneven multivector
//
// grade 1: gr1() - vector
// grade 3: gr3() - trivector

template <typename T> inline constexpr Vec4d<T> gr1(MVec4d_U<T> const& v)
{
    return Vec4d<T>(v.c0, v.c1, v.c2, v.c3);
}

template <typename T> inline constexpr TriVec4d<T> gr3(MVec4d_U<T> const& v)
{
    return TriVec4d<T>(v.c4, v.c5, v.c6, v.c7);
}

} // namespace hd::ga