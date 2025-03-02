#pragma once

// author: Daniel Hug, 2024 & 2025

#include "type_t/ga_type_3dp.hpp"

#include "ga_mvec3dp_e.hpp"
#include "ga_mvec3dp_u.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dp_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec16_t<T, mvec3dp_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    ////////////////////////////////////////////////////////////////////////////
    // ctors for MVec3dp<T>
    ////////////////////////////////////////////////////////////////////////////

    // assign a scalar part exclusively (other grades = 0)
    MVec16_t(Scalar3dp<T> s) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec16_t(Vec3dp<T> const& v) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec16_t(BiVec3dp<T> const& B) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my,
                 B.mz, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    MVec16_t(TriVec3dp<T> const& t) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec16_t(PScalar3dp<T> ps) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade subvector
    MVec16_t(Scalar3dp<T> s, BiVec3dp<T> const& B, PScalar3dp<T> ps) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz,
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade subvector
    MVec16_t(MVec3dp_E<T> const& M) :
        MVec16_t(M.c0, T(0.0), T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6,
                 T(0.0), T(0.0), T(0.0), T(0.0), M.c7)
    {
    }

    // assign components of an uneven grade subvector
    MVec16_t(Vec3dp<T> const& v, TriVec3dp<T> const& t) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign components of an uneven grade subvector
    MVec16_t(MVec3dp_U<T> const& M) :
        MVec16_t(T(0.0), M.c0, M.c1, M.c2, M.c3, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), M.c4, M.c5, M.c6, M.c7, T(0.0))
    {
    }

    // assign a full multivector
    // assign components of an even grade subvector
    MVec16_t(Scalar3dp<T> s, Vec3dp<T> const& v, BiVec3dp<T> const& B,
             TriVec3dp<T> const& t, PScalar3dp<T> ps) :
        MVec16_t(T(s), v.x, v.y, v.z, v.w, B.vx, B.vy, B.vz, B.mx, B.my, B.mz, t.x, t.y,
                 t.z, t.w, T(ps))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec16_t<T, mvec3dp_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector
// grade 4: gr4() - quadvector (= pseudoscalar in 3dp)

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> gr0(MVec3dp<T> const& M)
{
    return Scalar3dp<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> gr1(MVec3dp<T> const& M)
{
    return Vec3dp<T>(M.c1, M.c2, M.c3, M.c4);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> gr2(MVec3dp<T> const& M)
{
    return BiVec3dp<T>(M.c5, M.c6, M.c7, M.c8, M.c9, M.c10);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> gr3(MVec3dp<T> const& M)
{
    return TriVec3dp<T>(M.c11, M.c12, M.c13, M.c14);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> gr4(MVec3dp<T> const& M)
{
    return PScalar3dp<T>(M.c15);
}

// return the grades of the basic types

template <typename T>
    requires(std::floating_point<T>)
inline constexpr size_t gr([[maybe_unused]] Scalar3dp<T>)
{
    return 0;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr size_t gr([[maybe_unused]] Vec3dp<T> const&)
{
    return 1;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr size_t gr([[maybe_unused]] BiVec3dp<T> const&)
{
    return 2;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr size_t gr([[maybe_unused]] TriVec3dp<T> const&)
{
    return 3;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr size_t gr([[maybe_unused]] PScalar3dp<T>)
{
    return 4;
}

} // namespace hd::ga