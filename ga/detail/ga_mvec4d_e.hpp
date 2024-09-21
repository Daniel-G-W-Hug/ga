#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_4d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec4d_e_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec4d_e_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec8_t(Scalar4d<T> s) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec8_t(BiVec4d<T> const& v) :
        MVec8_t(T(0.0), v.vx, v.vy, v.vz, v.mx, v.my, v.mz, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar4d<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2), cmt(v1,v2) and wdg(v1,v2)
    MVec8_t(Scalar4d<T> s, BiVec4d<T> const& v, PScalar4d<T> ps) :
        MVec8_t(T(s), v.vx, v.vy, v.vz, v.mx, v.my, v.mz, T(ps))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec4d_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of an even 4d multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector
// grade 3: gr4() - quadvector (= pseudoscalar in 4d)

template <typename T> inline constexpr Scalar4d<T> gr0(MVec4d_E<T> const& v)
{
    return Scalar4d<T>(v.c0);
}

template <typename T> inline constexpr BiVec4d<T> gr2(MVec4d_E<T> const& v)
{
    return BiVec4d<T>(v.c1, v.c2, v.c3, v.c4, v.c5, v.c6);
}

template <typename T> inline constexpr PScalar4d<T> gr3(MVec4d_E<T> const& v)
{
    return PScalar4d<T>(v.c7);
}

} // namespace hd::ga