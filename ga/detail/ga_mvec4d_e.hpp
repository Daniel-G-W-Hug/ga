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
    MVec8_t(BiVec4d<T> const& B) :
        MVec8_t(T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar4d<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2), cmt(v1,v2) and wdg(v1,v2)
    MVec8_t(Scalar4d<T> s, BiVec4d<T> const& B, PScalar4d<T> ps) :
        MVec8_t(T(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(ps))
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

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar4d<T> gr0(MVec4d_E<T> const& M)
{
    return Scalar4d<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec4d<T> gr2(MVec4d_E<T> const& M)
{
    return BiVec4d<T>(M.c1, M.c2, M.c3, M.c4, M.c5, M.c6);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar4d<T> gr4(MVec4d_E<T> const& M)
{
    return PScalar4d<T>(M.c7);
}

} // namespace hd::ga