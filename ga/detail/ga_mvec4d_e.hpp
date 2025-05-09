#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

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

    // assign the scalar and the bivector part only, assuming the pseudoscalar part as 0.0
    MVec8_t(Scalar4d<T> s, BiVec4d<T> const& B) :
        MVec8_t(T(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(0.0))
    {
    }

    // assign the bivector and the pseudoscalar part only, assuming the scalar part as 0.0
    MVec8_t(BiVec4d<T> const& B, PScalar4d<T> ps) :
        MVec8_t(T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(ps))
    {
    }

    // assign the scalar and the pseudoscalar part only, assuming the bivector part as 0.0
    MVec8_t(Scalar4d<T> s, PScalar4d<T> ps) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign all three parts, the scalar, the bivector and the pseudoscalar explicitly
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

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + bivector => even grade multivector (in 4d imply a zero pseudoscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4d_E<std::common_type_t<T, U>> operator+(Scalar4d<T> s,
                                                              BiVec4d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_E<ctype>(s, B, Pscalar4d<ctype>(0.0));
}

// bivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4d_E<std::common_type_t<T, U>> operator+(BiVec4d<T> const& B,
                                                              Scalar4d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_E<ctype>(s, B, Pscalar4d<ctype>(0.0));
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + bivector => even grade multivector (in 4d imply a zero pseudoscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4d_E<std::common_type_t<T, U>> operator-(Scalar4d<T> s,
                                                              BiVec4d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_E<ctype>(s, -B, Pscalar4d<ctype>(0.0));
}

// bivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4d_E<std::common_type_t<T, U>> operator-(BiVec4d<T> const& B,
                                                              Scalar4d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_E<ctype>(-s, B, Pscalar4d<ctype>(0.0));
}

} // namespace hd::ga