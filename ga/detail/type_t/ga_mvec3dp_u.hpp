#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type3dp.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dp_u_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec8_t<T, mvec3dp_u_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign the vector part exclusively (other grades = 0)
    constexpr MVec8_t(Vec3dp<T> const& v) :
        MVec8_t(v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign the trivector part exclusively (other grades = 0)
    constexpr MVec8_t(TriVec3dp<T> const& t) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), t.x, t.y, t.z, t.w)
    {
    }

    // assign from a vector and a trivector directly
    constexpr MVec8_t(Vec3dp<T> const& v, TriVec3dp<T> const& t) :
        MVec8_t(v.x, v.y, v.z, v.w, t.x, t.y, t.z, t.w)
    {
    }

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons
    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator+=(MVec8_t<U, mvec3dp_u_tag> const& v) noexcept
    {
        this->c0 += v.c0;
        this->c1 += v.c1;
        this->c2 += v.c2;
        this->c3 += v.c3;
        this->c4 += v.c4;
        this->c5 += v.c5;
        this->c6 += v.c6;
        this->c7 += v.c7;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator-=(MVec8_t<U, mvec3dp_u_tag> const& v) noexcept
    {
        this->c0 -= v.c0;
        this->c1 -= v.c1;
        this->c2 -= v.c2;
        this->c3 -= v.c3;
        this->c4 -= v.c4;
        this->c5 -= v.c5;
        this->c6 -= v.c6;
        this->c7 -= v.c7;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator*=(U s) noexcept
    {
        this->c0 *= s;
        this->c1 *= s;
        this->c2 *= s;
        this->c3 *= s;
        this->c4 *= s;
        this->c5 *= s;
        this->c6 *= s;
        this->c7 *= s;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 8 comp.");
        this->c0 /= s;
        this->c1 /= s;
        this->c2 /= s;
        this->c3 /= s;
        this->c4 /= s;
        this->c5 /= s;
        this->c6 /= s;
        this->c7 /= s;
        return *this;
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
    requires(numeric_type<T>)
constexpr Vec3dp<T> gr1(MVec3dp_U<T> const& M)
{
    return Vec3dp<T>(M.c0, M.c1, M.c2, M.c3);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> gr3(MVec3dp_U<T> const& M)
{
    return TriVec3dp<T>(M.c4, M.c5, M.c6, M.c7);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine vectors and trivectors to odd grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector + trivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, t);
}

// trivector + vector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, t);
}


// odd grade multivector + trivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(MVec3dp_U<T> const& M,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4 + t.x, M.c5 + t.y, M.c6 + t.z,
                            M.c7 + t.w);
}

// trivector + odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                        MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4 + t.x, M.c5 + t.y, M.c6 + t.z,
                            M.c7 + t.w);
}

// odd grade multivector + vector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(MVec3dp_U<T> const& M,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0 + v.x, M.c1 + v.y, M.c2 + v.z, M.c3 + v.w, M.c4, M.c5,
                            M.c6, M.c7);
}

// vector + odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
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
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, -t);
}

// trivector - vector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-v, t);
}

// odd grade multivector - trivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(MVec3dp_U<T> const& M,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4 - t.x, M.c5 - t.y, M.c6 - t.z,
                            M.c7 - t.w);
}

// trivector - odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                        MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4 + t.x, -M.c5 + t.y,
                            -M.c6 + t.z, -M.c7 + t.w);
}

// odd grade multivector - vector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(MVec3dp_U<T> const& M,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(M.c0 - v.x, M.c1 - v.y, M.c2 - v.z, M.c3 - v.w, M.c4, M.c5,
                            M.c6, M.c7);
}

// vector - odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                        MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-M.c0 + v.x, -M.c1 + v.y, -M.c2 + v.z, -M.c3 + v.w, -M.c4,
                            -M.c5, -M.c6, -M.c7);
}

} // namespace hd::ga