#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "type_t/ga_type3dp.hpp" // type aliases for scalars, vector, pseudoscalar, etc.

#include "ga_mvec3dp_e.hpp" // even grade multivector
#include "ga_mvec3dp_u.hpp" // uneven (or odd) grade multivector


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dp_tag
/////////////////////////////////////////////////////////////////////////////////////////

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

    // assign components of an even grade multivector
    MVec16_t(Scalar3dp<T> s, BiVec3dp<T> const& B, PScalar3dp<T> ps) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz,
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade multivector
    MVec16_t(MVec3dp_E<T> const& M) :
        MVec16_t(M.c0, T(0.0), T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6,
                 T(0.0), T(0.0), T(0.0), T(0.0), M.c7)
    {
    }

    // assign components of an odd grade multivector
    MVec16_t(Vec3dp<T> const& v, TriVec3dp<T> const& t) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign components of an odd grade multivector
    MVec16_t(MVec3dp_U<T> const& M) :
        MVec16_t(T(0.0), M.c0, M.c1, M.c2, M.c3, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), M.c4, M.c5, M.c6, M.c7, T(0.0))
    {
    }

    // assign a full multivector
    // assign components of an even grade multivector
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

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec3dp_e.hpp and mvec3dp_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(Scalar3dp<T> s,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, v, BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// vector + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, v, BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// scalar + trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(Scalar3dp<T> s,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), t,
                          PScalar3dp<ctype>(0.0));
}

// trivector + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), t,
                          PScalar3dp<ctype>(0.0));
}

// scalar + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(Scalar3dp<T> s,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), ps);
}

// pseudoscalar + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), ps);
}

// vector + bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), v, B,
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// bivector + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(BiVec3dp<T> const& B,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), v, B,
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// vector + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), v,
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), ps);
}

// pseudoscalar + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), v,
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), ps);
}

// bivector + trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(BiVec3dp<T> const& B,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0), B, t,
                          PScalar3dp<ctype>(0.0));
}

// trivector + bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0), B, t,
                          PScalar3dp<ctype>(0.0));
}

// pseudoscalar + trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), t, ps);
}

// trivector + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), t, ps);
}

// mulitivector + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(MVec3dp<T> const& M,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M) + s, gr1(M), gr2(M), gr3(M), gr4(M));
}

// mulitivector + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(MVec3dp<T> const& M,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M) + v, gr2(M), gr3(M), gr4(M));
}

// mulitivector + bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(MVec3dp<T> const& M,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M), gr2(M) + B, gr3(M), gr4(M));
}

// mulitivector + trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(MVec3dp<T> const& M,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) + t, gr4(M));
}

// mulitivector + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator+(MVec3dp<T> const& M,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M), gr2(M), gr3(M), gr4(M) + ps);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec3dp_e.hpp and mvec3dp_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(Scalar3dp<T> s,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, -v, BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// vector - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(-s, v, BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// scalar - trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(Scalar3dp<T> s,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), -t,
                          PScalar3dp<ctype>(0.0));
}

// trivector - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(-s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), t,
                          PScalar3dp<ctype>(0.0));
}

// scalar - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(Scalar3dp<T> s,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), -ps);
}

// pseudoscalar - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(-s, Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), ps);
}

// vector - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), v, -B,
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// bivector - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(BiVec3dp<T> const& B,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), -v, B,
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), PScalar3dp<ctype>(0.0));
}

// vector - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), v,
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), -ps);
}

// pseudoscalar - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), -v,
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          TriVec3dp<ctype>(0.0, 0.0, 0.0, 0.0), ps);
}

// bivector - trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(BiVec3dp<T> const& B,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0), B, -t,
                          PScalar3dp<ctype>(0.0));
}

// trivector - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0), -B, t,
                          PScalar3dp<ctype>(0.0));
}

// pseudoscalar - trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), -t, ps);
}

// trivector - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(Scalar3dp<ctype>(0.0), Vec3dp<ctype>(0.0, 0.0, 0.0, 0.0),
                          BiVec3dp<ctype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), t, -ps);
}

// mulitivector - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(MVec3dp<T> const& M,
                                                             Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M) - s, gr1(M), gr2(M), gr3(M), gr4(M));
}

// mulitivector - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(MVec3dp<T> const& M,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M) - v, gr2(M), gr3(M), gr4(M));
}

// mulitivector - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(MVec3dp<T> const& M,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M), gr2(M) - B, gr3(M), gr4(M));
}

// mulitivector - trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(MVec3dp<T> const& M,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) - t, gr4(M));
}

// mulitivector - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator-(MVec3dp<T> const& M,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(gr0(M), gr1(M), gr2(M), gr3(M), gr4(M) - ps);
}

} // namespace hd::ga