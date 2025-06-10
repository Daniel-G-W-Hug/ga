#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "type_t/ga_type3dp.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dp_u_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec3dp_u_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign the vector part exclusively (other grades = 0)
    MVec8_t(Vec3dp<T> const& v) :
        MVec8_t(v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign the trivector part exclusively (other grades = 0)
    MVec8_t(TriVec3dp<T> const& t) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), t.x, t.y, t.z, t.w)
    {
    }

    // assign from a vector and a trivector directly
    MVec8_t(Vec3dp<T> const& v, TriVec3dp<T> const& t) :
        MVec8_t(v.x, v.y, v.z, v.w, t.x, t.y, t.z, t.w)
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec3dp_u_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a 3dp odd multivector
//
// grade 1: gr1() - vector
// grade 3: gr3() - trivector

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> gr1(MVec3dp_U<T> const& M)
{
    return Vec3dp<T>(M.c0, M.c1, M.c2, M.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> gr3(MVec3dp_U<T> const& M)
{
    return TriVec3dp<T>(M.c4, M.c5, M.c6, M.c7);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine vectors and trivectors to odd grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector + trivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, t);
}

// trivector + vector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, t);
}


// odd grade multivector + trivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(MVec3dp_U<T> const& M,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4 + t.c0, M.c5 + t.c1, M.c6 + t.c2,
                            M.c7 + t.c3);
}

// trivector + odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                               MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4 + t.c0, M.c5 + t.c1, M.c6 + t.c2,
                            M.c7 + t.c3);
}

// odd grade multivector + vector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(MVec3dp_U<T> const& M,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0 + v.x, M.c1 + v.y, M.c2 + v.z, M.c3 + v.w, M.c4, M.c5,
                            M.c6, M.c7);
}

// vector + odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
                                                               MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0 + v.x, M.c1 + v.y, M.c2 + v.z, M.c3 + v.w, M.c4, M.c5,
                            M.c6, M.c7);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine vectors and trivectors to odd grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector - trivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, -t);
}

// trivector - vector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-v, t);
}

// odd grade multivector - trivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(MVec3dp_U<T> const& M,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4 - t.c0, M.c5 - t.c1, M.c6 - t.c2,
                            M.c7 - t.c3);
}

// trivector - odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                               MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4 + t.c0, -M.c5 + t.c1,
                            -M.c6 + t.c2, -M.c7 + t.c3);
}

// odd grade multivector - vector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(MVec3dp_U<T> const& M,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0 - v.x, M.c1 - v.y, M.c2 - v.z, M.c3 - v.w, M.c4, M.c5,
                            M.c6, M.c7);
}

// vector - odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                               MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-M.c0 + v.x, -M.c1 + v.y, -M.c2 + v.z, -M.c3 + v.w, -M.c4,
                            -M.c5, -M.c6, -M.c7);
}

} // namespace hd::ga