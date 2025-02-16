#pragma once

// author: Daniel Hug, 2024 & 2025

#include "type_t/ga_type_3dp.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dp_e_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec3dp_e_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec8_t(Scalar3dp<T> s) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec8_t(BiVec3dp<T> const& B) :
        MVec8_t(T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar3dp<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2), cmt(v1,v2) and wdg(v1,v2)
    MVec8_t(Scalar3dp<T> s, BiVec3dp<T> const& B, PScalar3dp<T> ps) :
        MVec8_t(T(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(ps))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec3dp_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of an even 3dp multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector
// grade 3: gr4() - quadvector (= pseudoscalar in 3dp)

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> gr0(MVec3dp_E<T> const& M)
{
    return Scalar3dp<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> gr2(MVec3dp_E<T> const& M)
{
    return BiVec3dp<T>(M.c1, M.c2, M.c3, M.c4, M.c5, M.c6);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> gr4(MVec3dp_E<T> const& M)
{
    return PScalar3dp<T>(M.c7);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + bivector => even grade multivector (in 3dp imply a zero pseudoscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(Scalar3dp<T> s,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, B, Pscalar3dp<ctype>(0.0));
}

// bivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(BiVec3dp<T> const& B,
                                                               Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, B, Pscalar3dp<ctype>(0.0));
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + bivector => even grade multivector (in 3dp imply a zero pseudoscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(Scalar3dp<T> s,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, -B, Pscalar3dp<ctype>(0.0));
}

// bivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(BiVec3dp<T> const& B,
                                                               Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-s, B, Pscalar3dp<ctype>(0.0));
}

} // namespace hd::ga