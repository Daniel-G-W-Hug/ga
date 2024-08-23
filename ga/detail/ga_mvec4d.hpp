#pragma once

// author: Daniel Hug, 2024

#include "ga_type_tags.hpp"

#include "ga_mvec16_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using MVec4d = MVec16_t<T, mvec4d_tag>;

////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec4d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec16_t<T, mvec4d_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t(MVec16_t<U, mvec4d_tag> const& v) :
        MVec16_t(v.c0, v.c1, v.c2, v.c3, v.c4, v.c5, v.c6, v.c7, v.c8, v.c9, v.c10, v.c11,
                 v.c12, v.c13, v.c14, v.c15)
    {
    }

    ////////////////////////////////////////////////////////////////////////////
    // ctors for MVec4d<T>
    ////////////////////////////////////////////////////////////////////////////

    // assign a scalar part exclusively (other grades = 0)
    MVec16_t(Scalar<T> s) :
        MVec16_t(s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0)
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec16_t(Vec4d<T> const& v) :
        MVec16_t(0.0, v.x, v.y, v.z, v.w, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0)
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec16_t(BiVec4d<T> const& v) :
        MVec16_t(0.0, 0.0, 0.0, 0.0, 0.0, v.vx, v.vy, v.vz, v.mx, v.my, v.mz, 0.0, 0.0,
                 0.0, 0.0, 0.0)
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    MVec16_t(TriVec4d<T> const& t) :
        MVec16_t(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, t.x, t.y, t.z,
                 t.w, 0.0)
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec16_t(PScalar4d<T> ps) :
        MVec16_t(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, ps)
    {
    }

    // assign components of an even grade subvector
    MVec16_t(Scalar<T> s, BiVec4d<T> const& v, PScalar4d<T> ps) :
        MVec16_t(s, 0.0, 0.0, 0.0, 0.0, v.vx, v.vy, v.vz, v.mx, v.my, v.mz, 0.0, 0.0, 0.0,
                 0.0, ps)
    {
    }

    // assign components of an uneven grade subvector
    MVec16_t(Vec4d<T> const& v, TriVec4d<T> const& t) :
        MVec16_t(0.0, v.x, v.y, v.z, v.w, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, t.x, t.y, t.z,
                 t.w, 0.0)
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec16_t<T, mvec4d_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector
// grade 4: gr4() - quadvector (= pseudoscalar in 4d)

template <typename T> inline constexpr Scalar<T> gr0(MVec4d<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr Vec4d<T> gr1(MVec4d<T> const& v)
{
    return Vec4d<T>(v.c1, v.c2, v.c3, v.c4);
}

template <typename T> inline constexpr BiVec3d<T> gr2(MVec4d<T> const& v)
{
    return BiVec4d<T>(v.c5, v.c6, v.c7, v.c8, v.c9, v.c10);
}

template <typename T> inline constexpr TriVec4d<T> gr3(MVec4d<T> const& v)
{
    return TriVec4d<T>(v.c11, v.c12, v.c13, v.c14);
}

template <typename T> inline constexpr PScalar4d<T> gr4(MVec4d<T> const& v)
{
    return PScalar4d<T>(v.c15);
}


} // namespace hd::ga