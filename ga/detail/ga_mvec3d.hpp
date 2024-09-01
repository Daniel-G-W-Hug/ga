#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_mvec8_t.hpp"

#include "type_t/ga_type_0d.hpp"
#include "type_t/ga_type_3d.hpp"

#include "ga_mvec3d_e.hpp"
#include "ga_mvec3d_u.hpp"


namespace hd::ga {

template <typename T> using MVec3d = MVec8_t<T, mvec3d_tag>;

////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec3d_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec8_t(Scalar<T> s) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec8_t(Vec3d<T> const& v) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec8_t(BiVec3d<T> const& v) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), v.x, v.y, v.z, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar3d<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) or via dot(v1,v2) and cmt(v1,v2) directly
    // (other grades = 0)
    MVec8_t(Scalar<T> s, BiVec3d<T> const& v) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), v.x, v.y, v.z, T(0.0))
    {
    }

    // assign from a quaternion, i.e. from the even subalgebra
    MVec8_t(MVec3d_E<T> const& v) :
        MVec8_t(v.c0, T(0.0), T(0.0), T(0.0), v.c1, v.c2, v.c3, T(0.0))
    {
    }

    // assign a geometric product resulting from a product of a vector and a bivector
    MVec8_t(Vec3d<T> const& v, PScalar3d<T> ps) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign from the uneven subalgebra
    MVec8_t(MVec3d_U<T> const& v) :
        MVec8_t(T(0.0), v.c0, v.c1, v.c2, T(0.0), T(0.0), T(0.0), v.c3)
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec3d_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a 3d multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector (= pseudoscalar in 3d)

template <typename T> inline constexpr Scalar<T> gr0(MVec3d<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr Vec3d<T> gr1(MVec3d<T> const& v)
{
    return Vec3d<T>(v.c1, v.c2, v.c3);
}

template <typename T> inline constexpr BiVec3d<T> gr2(MVec3d<T> const& v)
{
    return BiVec3d<T>(v.c4, v.c5, v.c6);
}

template <typename T> inline constexpr PScalar3d<T> gr3(MVec3d<T> const& v)
{
    return PScalar3d<T>(v.c7);
}

} // namespace hd::ga