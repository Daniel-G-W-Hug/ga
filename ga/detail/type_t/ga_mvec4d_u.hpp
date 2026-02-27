#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type4d.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec4d_u_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec8_t<T, mvec4d_u_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign the vector part exclusively (other grades = 0)
    constexpr MVec8_t(Vec4d<T> const& v) :
        MVec8_t(v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign the trivector part exclusively (other grades = 0)
    constexpr MVec8_t(TriVec4d<T> const& t) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), t.x, t.y, t.z, t.w)
    {
    }

    // assign from a vector and a trivector directly
    constexpr MVec8_t(Vec4d<T> const& v, TriVec4d<T> const& t) :
        MVec8_t(v.x, v.y, v.z, v.w, t.x, t.y, t.z, t.w)
    {
    }

    // assign from the odd subalgebra directly
    constexpr MVec8_t(MVec4d_U<T> const& M) :
        MVec8_t(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7)
    {
    }

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons
    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator+=(MVec8_t<U, mvec4d_u_tag> const& v) noexcept
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
    MVec8_t& operator-=(MVec8_t<U, mvec4d_u_tag> const& v) noexcept
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
// define grade operations for partial specialization MVec8_t<T, mvec4d_u_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a 4d odd multivector
//
// grade 1: gr1() - vector
// grade 3: gr3() - trivector

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4d<T> gr1(MVec4d_U<T> const& M)
{
    return Vec4d<T>(M.c0, M.c1, M.c2, M.c3);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4d<T> gr3(MVec4d_U<T> const& M)
{
    return TriVec4d<T>(M.c4, M.c5, M.c6, M.c7);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine vectors and trivectors to odd grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector + trivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4d_U<std::common_type_t<T, U>> operator+(Vec4d<T> const& v,
                                                       TriVec4d<T> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_U<ctype>(v, t);
}

// trivector + vector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4d_U<std::common_type_t<T, U>> operator+(TriVec4d<T> const& t,
                                                       Vec4d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_U<ctype>(v, t);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine vectors and trivectors to odd grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector + trivector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4d_U<std::common_type_t<T, U>> operator-(Vec4d<T> const& v,
                                                       TriVec4d<T> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_U<ctype>(v, -t);
}

// trivector + vector => odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4d_U<std::common_type_t<T, U>> operator-(TriVec4d<T> const& t,
                                                       Vec4d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4d_U<ctype>(-v, t);
}

} // namespace hd::ga