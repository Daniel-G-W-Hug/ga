#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type4d.hpp" // type aliases for scalars, vector, pseudoscalar, etc.

#include "ga_mvec4d_e.hpp" // even grade multivector
#include "ga_mvec4d_u.hpp" // uneven (or odd) grade multivector


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec4d_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec16_t<T, mvec4d_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    ////////////////////////////////////////////////////////////////////////////
    // ctors for MVec4d<T>
    ////////////////////////////////////////////////////////////////////////////

    // assign a scalar part exclusively (other grades = 0)
    MVec16_t(Scalar4d<T> s) :
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
    MVec16_t(BiVec4d<T> const& B) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my,
                 B.mz, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
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
    MVec16_t(Scalar4d<T> s, BiVec4d<T> const& B, PScalar4d<T> ps) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz,
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade subvector
    MVec16_t(MVec4d_E<T> const& M) :
        MVec16_t(M.c0, T(0.0), T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6,
                 T(0.0), T(0.0), T(0.0), T(0.0), M.c7)
    {
    }

    // assign components of an odd grade subvector
    MVec16_t(Vec4d<T> const& v, TriVec4d<T> const& t) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign components of an odd grade subvector
    MVec16_t(MVec4d_U<T> const& M) :
        MVec16_t(T(0.0), M.c0, M.c1, M.c2, M.c3, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), M.c4, M.c5, M.c6, M.c7, T(0.0))
    {
    }

    // assign a full multivector
    // assign components of an even grade subvector
    MVec16_t(Scalar4d<T> s, Vec4d<T> const& v, BiVec4d<T> const& B, TriVec4d<T> const& t,
             PScalar4d<T> ps) :
        MVec16_t(T(s), v.x, v.y, v.z, v.w, B.vx, B.vy, B.vz, B.mx, B.my, B.mz, t.x, t.y,
                 t.z, t.w, T(ps))
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

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar4d<T> gr0(MVec4d<T> const& M)
{
    return Scalar4d<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec4d<T> gr1(MVec4d<T> const& M)
{
    return Vec4d<T>(M.c1, M.c2, M.c3, M.c4);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec4d<T> gr2(MVec4d<T> const& M)
{
    return BiVec4d<T>(M.c5, M.c6, M.c7, M.c8, M.c9, M.c10);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec4d<T> gr3(MVec4d<T> const& M)
{
    return TriVec4d<T>(M.c11, M.c12, M.c13, M.c14);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar4d<T> gr4(MVec4d<T> const& M)
{
    return PScalar4d<T>(M.c15);
}

// return the grades of the basic types

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] Scalar4d<T>)
{
    return 0;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] Vec4d<T> const&)
{
    return 1;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] BiVec4d<T> const&)
{
    return 2;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] TriVec4d<T> const&)
{
    return 3;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] PScalar4d<T>)
{
    return 4;
}

} // namespace hd::ga