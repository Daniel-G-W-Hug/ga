#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type3d.hpp" // type aliases for scalars, vector, pseudoscalar, etc.

#include "ga_mvec3d_e.hpp" // even grade multivector
#include "ga_mvec3d_u.hpp" // uneven (or odd) grade multivector


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3d_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec3d_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec8_t(Scalar3d<T> s) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec8_t(Vec3d<T> const& v) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec8_t(BiVec3d<T> const& B) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), B.x, B.y, B.z, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar3d<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign a scalar and a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(Scalar3d<T> s, PScalar3d<T> ps) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) or via dot(v1,v2) and cmt(v1,v2) directly
    // (other grades = 0)
    MVec8_t(Scalar3d<T> s, BiVec3d<T> const& B) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), B.x, B.y, B.z, T(0.0))
    {
    }

    // assign from a quaternion, i.e. from the even subalgebra
    MVec8_t(MVec3d_E<T> const& M) :
        MVec8_t(M.c0, T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, T(0.0))
    {
    }

    // assign a geometric product resulting from a product of a vector and a bivector
    MVec8_t(Vec3d<T> const& v, PScalar3d<T> ps) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign from the odd subalgebra
    MVec8_t(MVec3d_U<T> const& M) :
        MVec8_t(T(0.0), M.c0, M.c1, M.c2, T(0.0), T(0.0), T(0.0), M.c3)
    {
    }

    // assign a full multivector
    MVec8_t(Scalar3d<T> s, Vec3d<T> const& v, BiVec3d<T> const& B, PScalar3d<T> ps) :
        MVec8_t(T(s), v.x, v.y, v.z, B.x, B.y, B.z, T(ps))
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

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3d<T> gr0(MVec3d<T> const& M)
{
    return Scalar3d<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3d<T> gr1(MVec3d<T> const& M)
{
    return Vec3d<T>(M.c1, M.c2, M.c3);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3d<T> gr2(MVec3d<T> const& M)
{
    return BiVec3d<T>(M.c4, M.c5, M.c6);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3d<T> gr3(MVec3d<T> const& M)
{
    return PScalar3d<T>(M.c7);
}

// return the grades of the basic types

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] Scalar3d<T>)
{
    return 0;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] Vec3d<T> const&)
{
    return 1;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] BiVec3d<T> const&)
{
    return 2;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] PScalar3d<T>)
{
    return 3;
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec3d_e.hpp and mvec3d_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar + vector => multivector
// (bivector and pseudoscalar are set to zero)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(Scalar3d<T> s, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(s, v, BiVec3d<ctype>(0.0, 0.0, 0.0), PScalar3d<ctype>(0.0));
}

// vector + scalar  => multivector
// (bivector and pseudoscalar are set to zero)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(Vec3d<T> const& v, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(s, v, BiVec3d<ctype>(0.0, 0.0, 0.0), PScalar3d<ctype>(0.0));
}

// scalar + pseudoscalar (=trivector) => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(s, Vec3d<ctype>(0.0, 0.0, 0.0), BiVec3d<ctype>(0.0, 0.0, 0.0),
                         ps);
}

// pseudoscalar (=trivector) + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(s, Vec3d<ctype>(0.0, 0.0, 0.0), BiVec3d<ctype>(0.0, 0.0, 0.0),
                         ps);
}

// vector + bivector => multivector
// (scalar and pseudoscalar are set to zero)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(Vec3d<T> const& v,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), v, B, PScalar3d<ctype>(0.0));
}

// bivector + vector  => multivector
// (scalar and pseudoscalar are set to zero)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(BiVec3d<T> const& B,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), v, B, PScalar3d<ctype>(0.0));
}

// pseudoscalar + bivector => multivector
// (scalar and vector are set to zero)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(PScalar3d<T> ps, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), Vec3d<ctype>(0.0, 0.0, 0.0), B, ps);
}

// bivector + pseudoscalar  => multivector
// (scalar and vector are set to zero)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(BiVec3d<T> const& B, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), Vec3d<ctype>(0.0, 0.0, 0.0), B, ps);
}

// mulitivector + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(MVec3d<T> const& M, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M) + s, gr1(M), gr2(M), gr3(M));
}

// mulitivector + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(MVec3d<T> const& M,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M), gr1(M) + v, gr2(M), gr3(M));
}

// mulitivector + bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(MVec3d<T> const& M,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M), gr1(M), gr2(M) + B, gr3(M));
}

// mulitivector + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator+(MVec3d<T> const& M, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) + ps);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec3d_e.hpp and mvec3d_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(Scalar3d<T> s, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(s, -v, BiVec3d<ctype>(0.0, 0.0, 0.0), PScalar3d<ctype>(0.0));
}

// vector - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(Vec3d<T> const& v, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(-s, v, BiVec3d<ctype>(0.0, 0.0, 0.0), PScalar3d<ctype>(0.0));
}

// scalar - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(s, Vec3d<ctype>(0.0, 0.0, 0.0), BiVec3d<ctype>(0.0, 0.0, 0.0),
                         -ps);
}

// pseudoscalar - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(-s, Vec3d<ctype>(0.0, 0.0, 0.0), BiVec3d<ctype>(0.0, 0.0, 0.0),
                         ps);
}

// vector - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(Vec3d<T> const& v,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), v, -B, PScalar3d<ctype>(0.0));
}

// bivector - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(BiVec3d<T> const& B,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), -v, B, PScalar3d<ctype>(0.0));
}

// pseudoscalar - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(PScalar3d<T> ps, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), Vec3d<ctype>(0.0, 0.0, 0.0), -B, ps);
}

// bivector - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(BiVec3d<T> const& B, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar3d<ctype>(0.0), Vec3d<ctype>(0.0, 0.0, 0.0), B, -ps);
}

// mulitivector - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(MVec3d<T> const& M, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M) - s, gr1(M), gr2(M), gr3(M));
}

// mulitivector - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(MVec3d<T> const& M,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M), gr1(M) - v, gr2(M), gr3(M));
}

// mulitivector - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(MVec3d<T> const& M,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M), gr1(M), gr2(M) - B, gr3(M));
}

// mulitivector - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator-(MVec3d<T> const& M, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) - ps);
}

} // namespace hd::ga