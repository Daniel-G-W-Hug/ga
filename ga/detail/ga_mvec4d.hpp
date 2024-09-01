#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_mvec16_t.hpp"

#include "type_t/ga_type_0d.hpp"
#include "type_t/ga_type_4d.hpp"

#include "ga_mvec4d_e.hpp"
#include "ga_mvec4d_u.hpp"


namespace hd::ga {

template <typename T> using MVec4d = MVec16_t<T, mvec4d_tag>;

////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec4d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec16_t<T, mvec4d_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    ////////////////////////////////////////////////////////////////////////////
    // ctors for MVec4d<T>
    ////////////////////////////////////////////////////////////////////////////

    // assign a scalar part exclusively (other grades = 0)
    MVec16_t(Scalar<T> s) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec16_t(Vec4d<T> const& v) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec16_t(BiVec4d<T> const& v) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), v.vx, v.vy, v.vz, v.mx, v.my,
                 v.mz, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    MVec16_t(TriVec4d<T> const& t) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec16_t(PScalar4d<T> ps) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade subvector
    MVec16_t(Scalar<T> s, BiVec4d<T> const& v, PScalar4d<T> ps) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), v.vx, v.vy, v.vz, v.mx, v.my, v.mz,
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade subvector
    MVec16_t(MVec4d_E<T> v) :
        MVec16_t(v.c0, T(0.0), T(0.0), T(0.0), T(0.0), v.c1, v.c2, v.c3, v.c4, v.c5, v.c6,
                 T(0.0), T(0.0), T(0.0), T(0.0), v.c7)
    {
    }

    // assign components of an uneven grade subvector
    MVec16_t(Vec4d<T> const& v, TriVec4d<T> const& t) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign components of an uneven grade subvector
    MVec16_t(MVec4d_U<T> v) :
        MVec16_t(T(0.0), v.c0, v.c1, v.c2, v.c3, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), v.c4, v.c5, v.c6, v.c7, T(0.0))
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

template <typename T> inline constexpr BiVec4d<T> gr2(MVec4d<T> const& v)
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