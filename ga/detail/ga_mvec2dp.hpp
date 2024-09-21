#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_2dp.hpp"

#include "ga_mvec2dp_e.hpp"
#include "ga_mvec2dp_u.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec2dp_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec8_t(Scalar2dp<T> s) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec8_t(Vec2dp<T> const& v) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec8_t(BiVec2dp<T> const& v) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), v.x, v.y, v.z, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar2dp<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) or via dot(v1,v2) and cmt(v1,v2) directly
    // (other grades = 0)
    MVec8_t(Scalar2dp<T> s, BiVec2dp<T> const& b) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), b.x, b.y, b.z, T(0.0))
    {
    }

    // assign from a quaternion, i.e. from the even subalgebra
    MVec8_t(MVec2dp_E<T> const& v) :
        MVec8_t(v.c0, T(0.0), T(0.0), T(0.0), v.c1, v.c2, v.c3, T(0.0))
    {
    }

    // assign a geometric product resulting from a product of a vector and a bivector
    MVec8_t(Vec2dp<T> const& v, PScalar2dp<T> ps) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign from the uneven subalgebra
    MVec8_t(MVec2dp_U<T> const& v) :
        MVec8_t(T(0.0), v.c0, v.c1, v.c2, T(0.0), T(0.0), T(0.0), v.c3)
    {
    }

    // assign a full multivector
    MVec8_t(Scalar2dp<T> s, Vec2dp<T> const& v, BiVec2dp<T> const& b, PScalar2dp<T> ps) :
        MVec8_t(T(s), v.x, v.y, v.z, b.x, b.y, b.z, T(ps))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec2dp_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a 3d multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector (= pseudoscalar in 3d)

template <typename T> inline constexpr Scalar2dp<T> gr0(MVec2dp<T> const& v)
{
    return Scalar2dp<T>(v.c0);
}

template <typename T> inline constexpr Vec2dp<T> gr1(MVec2dp<T> const& v)
{
    return Vec2dp<T>(v.c1, v.c2, v.c3);
}

template <typename T> inline constexpr BiVec2dp<T> gr2(MVec2dp<T> const& v)
{
    return BiVec2dp<T>(v.c4, v.c5, v.c6);
}

template <typename T> inline constexpr PScalar2dp<T> gr3(MVec2dp<T> const& v)
{
    return PScalar2dp<T>(v.c7);
}

} // namespace hd::ga