#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::ostream

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> definition of a multivector with 16 components c0, ..., c15
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(std::floating_point<T>)
struct MVec16_t {

    // ctors

    // (all grades = 0)
    MVec16_t() = default;

    // assign all components directly
    constexpr MVec16_t(T s, T vx, T vy, T vz, T vw, T bu, T bv, T bw, T bx, T by, T bz,
                       T tx, T ty, T tz, T tw, T ps) :
        c0(s), c1(vx), c2(vy), c3(vz), c4(vw), c5(bu), c6(bv), c7(bw), c8(bx), c9(by),
        c10(bz), c11(tx), c12(ty), c13(tz), c14(tw), c15(ps)
    {
    }

    // cp & mv ctor
    MVec16_t(MVec16_t const&) = default;
    MVec16_t(MVec16_t&&) noexcept = default;

    // cp & mv assignment
    MVec16_t& operator=(MVec16_t const&) = default;
    MVec16_t& operator=(MVec16_t&&) noexcept = default;

    // dtor
    ~MVec16_t() = default;

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    constexpr MVec16_t(MVec16_t<U, Tag> const& v) :
        MVec16_t(v.c0, v.c1, v.c2, v.c3, v.c4, v.c5, v.c6, v.c7, v.c8, v.c9, v.c10, v.c11,
                 v.c12, v.c13, v.c14, v.c15)
    {
    }

    friend void swap(MVec16_t& lhs, MVec16_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.c0, rhs.c0);
        swap(lhs.c1, rhs.c1);
        swap(lhs.c2, rhs.c2);
        swap(lhs.c3, rhs.c3);
        swap(lhs.c4, rhs.c4);
        swap(lhs.c5, rhs.c5);
        swap(lhs.c6, rhs.c6);
        swap(lhs.c7, rhs.c7);
        swap(lhs.c8, rhs.c8);
        swap(lhs.c9, rhs.c9);
        swap(lhs.c10, rhs.c10);
        swap(lhs.c11, rhs.c11);
        swap(lhs.c12, rhs.c12);
        swap(lhs.c13, rhs.c13);
        swap(lhs.c14, rhs.c14);
        swap(lhs.c15, rhs.c15);
    }

    ////////////////////////////////////////////////////////////////////////////
    // component definition
    ////////////////////////////////////////////////////////////////////////////

    T c0{};  //     scalar
             //
    T c1{};  //     vector 4d, 1st component    - maps to basis vector     e1
    T c2{};  //     vector 4d, 2nd component    - maps to basis vector     e2
    T c3{};  //     vector 4d, 3rd component    - maps to basis vector     e3
    T c4{};  //     vector 4d, 4th component    - maps to basis vector     e4
             //
    T c5{};  //   bivector 4d, 1st component    - maps to basis bivector   e4^e1
    T c6{};  //   bivector 4d, 2nd component    - maps to basis bivector   e4^e2
    T c7{};  //   bivector 4d, 3rd component    - maps to basis bivector   e4^e3
    T c8{};  //   bivector 4d, 4th component    - maps to basis bivector   e2^e3
    T c9{};  //   bivector 4d, 5th component    - maps to basis bivector   e3^e1
    T c10{}; //   bivector 4d, 6th component    - maps to basis bivector   e1^e2
             //
    T c11{}; //  trivector 4d, 1st component    - maps to basis trivector  e4^e2^e3
    T c12{}; //  trivector 4d, 2nd component    - maps to basis trivector  e4^e3^e1
    T c13{}; //  trivector 4d, 3rd component    - maps to basis trivector  e4^e1^e2
    T c14{}; //  trivector 4d, 4th component    - maps to basis trivector  e3^e2^e1
             //
    T c15{}; // quadvector 4d = 4d pseudoscalar - maps to basis quadvector e1^e2^e3^e4

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec16_t<U, Tag> const& rhs) const
    {
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto abs_delta_c2 = std::abs(rhs.c2 - c2);
        auto abs_delta_c3 = std::abs(rhs.c3 - c3);
        auto abs_delta_c4 = std::abs(rhs.c4 - c4);
        auto abs_delta_c5 = std::abs(rhs.c5 - c5);
        auto abs_delta_c6 = std::abs(rhs.c6 - c6);
        auto abs_delta_c7 = std::abs(rhs.c7 - c7);
        auto abs_delta_c8 = std::abs(rhs.c8 - c8);
        auto abs_delta_c9 = std::abs(rhs.c9 - c9);
        auto abs_delta_c10 = std::abs(rhs.c10 - c10);
        auto abs_delta_c11 = std::abs(rhs.c11 - c11);
        auto abs_delta_c12 = std::abs(rhs.c12 - c12);
        auto abs_delta_c13 = std::abs(rhs.c13 - c13);
        auto abs_delta_c14 = std::abs(rhs.c14 - c14);
        auto abs_delta_c15 = std::abs(rhs.c15 - c15);
        auto constexpr delta_eps = detail::safe_epsilon<T, U>();
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps &&
            abs_delta_c4 < delta_eps && abs_delta_c5 < delta_eps &&
            abs_delta_c6 < delta_eps && abs_delta_c7 < delta_eps &&
            abs_delta_c8 < delta_eps && abs_delta_c9 < delta_eps &&
            abs_delta_c10 < delta_eps && abs_delta_c11 < delta_eps &&
            abs_delta_c12 < delta_eps && abs_delta_c13 < delta_eps &&
            abs_delta_c14 < delta_eps && abs_delta_c15 < delta_eps)
            return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator+=(MVec16_t<U, Tag> const& v) noexcept
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        c4 += v.c4;
        c5 += v.c5;
        c6 += v.c6;
        c7 += v.c7;
        c8 += v.c8;
        c9 += v.c9;
        c10 += v.c10;
        c11 += v.c11;
        c12 += v.c12;
        c13 += v.c13;
        c14 += v.c14;
        c15 += v.c15;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator-=(MVec16_t<U, Tag> const& v) noexcept
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        c4 -= v.c4;
        c5 -= v.c5;
        c6 -= v.c6;
        c7 -= v.c7;
        c8 -= v.c8;
        c9 -= v.c9;
        c10 -= v.c10;
        c11 -= v.c11;
        c12 -= v.c12;
        c13 -= v.c13;
        c14 -= v.c14;
        c15 -= v.c15;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator*=(U s) noexcept
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        c4 *= s;
        c5 *= s;
        c6 *= s;
        c7 *= s;
        c8 *= s;
        c9 *= s;
        c10 *= s;
        c11 *= s;
        c12 *= s;
        c13 *= s;
        c14 *= s;
        c15 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 16 comp.");
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        c4 /= s;
        c5 /= s;
        c6 /= s;
        c7 /= s;
        c8 /= s;
        c9 /= s;
        c10 /= s;
        c11 /= s;
        c12 /= s;
        c13 /= s;
        c14 /= s;
        c15 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr MVec16_t<T, Tag> operator-(MVec16_t<T, Tag> const& v)
{
    return MVec16_t<T, Tag>(-v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, -v.c7, -v.c8,
                            -v.c9, -v.c10, -v.c11, -v.c12, -v.c13, -v.c14, -v.c15);
}

// add multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator+(MVec16_t<T, Tag> const& v1, MVec16_t<U, Tag> const& v2)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2, v1.c3 + v2.c3, v1.c4 + v2.c4,
        v1.c5 + v2.c5, v1.c6 + v2.c6, v1.c7 + v2.c7, v1.c8 + v2.c8, v1.c9 + v2.c9,
        v1.c10 + v2.c10, v1.c11 + v2.c11, v1.c12 + v2.c12, v1.c13 + v2.c13,
        v1.c14 + v2.c14, v1.c15 + v2.c15);
}

// substract multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator-(MVec16_t<T, Tag> const& v1, MVec16_t<U, Tag> const& v2)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2, v1.c3 - v2.c3, v1.c4 - v2.c4,
        v1.c5 - v2.c5, v1.c6 - v2.c6, v1.c7 - v2.c7, v1.c8 - v2.c8, v1.c9 - v2.c9,
        v1.c10 - v2.c10, v1.c11 - v2.c11, v1.c12 - v2.c12, v1.c13 - v2.c13,
        v1.c14 - v2.c14, v1.c15 - v2.c15);
}

// multiply a multivector with a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator*(MVec16_t<T, Tag> const& v, U s)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s, v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s,
        v.c8 * s, v.c9 * s, v.c10 * s, v.c11 * s, v.c12 * s, v.c13 * s, v.c14 * s,
        v.c15 * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator*(T s, MVec16_t<U, Tag> const& v)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s, v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s,
        v.c8 * s, v.c9 * s, v.c10 * s, v.c11 * s, v.c12 * s, v.c13 * s, v.c14 * s,
        v.c15 * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline MVec16_t<std::common_type_t<T, U>, Tag> operator/(MVec16_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "multivector division");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec16_t<ctype, Tag>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv,
                                v.c4 * inv, v.c5 * inv, v.c6 * inv, v.c7 * inv,
                                v.c8 * inv, v.c9 * inv, v.c10 * inv, v.c11 * inv,
                                v.c12 * inv, v.c13 * inv, v.c14 * inv, v.c15 * inv);
}

// magnitude of the sum of a k-vector (in representational space)
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm_sq(MVec16_t<T, Tag> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3 + v.c4 * v.c4 +
           v.c5 * v.c5 + v.c6 * v.c6 + v.c7 * v.c7 + v.c8 * v.c8 + v.c9 * v.c9 +
           v.c10 * v.c10 + v.c11 * v.c11 + v.c12 * v.c12 + v.c13 * v.c13 + v.c14 * v.c14 +
           v.c15 * v.c15;
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(MVec16_t<T, Tag> const& v)
{
    return std::sqrt(nrm_sq(v));
}

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec16_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << "," << v.c4 << ","
       << v.c5 << "," << v.c6 << "," << v.c7 << "," << v.c8 << "," << v.c9 << "," << v.c10
       << "," << v.c11 << "," << v.c12 << "," << v.c13 << "," << v.c14 << "," << v.c15
       << ")";
    return os;
}

} // namespace hd::ga