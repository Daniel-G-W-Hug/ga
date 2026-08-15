#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <algorithm> // std::max
#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // numeric_type<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// MVec32_t<T, Tag> definition of a multivector with 32 components c0, ..., c31
// (the full multivector of a 5-dimensional algebra; which blade each component
// maps to is defined by the algebra the Tag belongs to -- the component
// comments below state the cga3dc layout)
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(numeric_type<T>)
struct MVec32_t {

    // ctors

    // (all grades = 0)
    MVec32_t() = default;

    // assign all components directly
    constexpr MVec32_t(T s, T vx, T vy, T vz, T vw, T vu, T bvx, T bvy, T bvz, T bmx,
                       T bmy, T bmz, T bpx, T bpy, T bpz, T bpw, T tvx, T tvy, T tvz,
                       T tmx, T tmy, T tmz, T tpx, T tpy, T tpz, T tpw, T qx, T qy, T qz,
                       T qw, T qu, T ps) :
        c0(s), c1(vx), c2(vy), c3(vz), c4(vw), c5(vu), c6(bvx), c7(bvy), c8(bvz), c9(bmx),
        c10(bmy), c11(bmz), c12(bpx), c13(bpy), c14(bpz), c15(bpw), c16(tvx), c17(tvy),
        c18(tvz), c19(tmx), c20(tmy), c21(tmz), c22(tpx), c23(tpy), c24(tpz), c25(tpw),
        c26(qx), c27(qy), c28(qz), c29(qw), c30(qu), c31(ps)
    {
    }

    // cp & mv ctor
    MVec32_t(MVec32_t const&) = default;
    MVec32_t(MVec32_t&&) noexcept = default;

    // cp & mv assignment
    MVec32_t& operator=(MVec32_t const&) = default;
    MVec32_t& operator=(MVec32_t&&) noexcept = default;

    // dtor
    ~MVec32_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr MVec32_t(MVec32_t<U, Tag> const& v) :
        MVec32_t(v.c0, v.c1, v.c2, v.c3, v.c4, v.c5, v.c6, v.c7, v.c8, v.c9, v.c10, v.c11,
                 v.c12, v.c13, v.c14, v.c15, v.c16, v.c17, v.c18, v.c19, v.c20, v.c21,
                 v.c22, v.c23, v.c24, v.c25, v.c26, v.c27, v.c28, v.c29, v.c30, v.c31)
    {
    }

    friend void swap(MVec32_t& lhs, MVec32_t& rhs) noexcept
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
        swap(lhs.c16, rhs.c16);
        swap(lhs.c17, rhs.c17);
        swap(lhs.c18, rhs.c18);
        swap(lhs.c19, rhs.c19);
        swap(lhs.c20, rhs.c20);
        swap(lhs.c21, rhs.c21);
        swap(lhs.c22, rhs.c22);
        swap(lhs.c23, rhs.c23);
        swap(lhs.c24, rhs.c24);
        swap(lhs.c25, rhs.c25);
        swap(lhs.c26, rhs.c26);
        swap(lhs.c27, rhs.c27);
        swap(lhs.c28, rhs.c28);
        swap(lhs.c29, rhs.c29);
        swap(lhs.c30, rhs.c30);
        swap(lhs.c31, rhs.c31);
    }

    ////////////////////////////////////////////////////////////////////////////
    // component definition (cga3dc layout)
    ////////////////////////////////////////////////////////////////////////////

    T c0{};  // scalar
    T c1{};  // vector,       maps to basis blade e1
    T c2{};  // vector,       maps to basis blade e2
    T c3{};  // vector,       maps to basis blade e3
    T c4{};  // vector,       maps to basis blade e4
    T c5{};  // vector,       maps to basis blade e5
    T c6{};  // bivector,     maps to basis blade e41
    T c7{};  // bivector,     maps to basis blade e42
    T c8{};  // bivector,     maps to basis blade e43
    T c9{};  // bivector,     maps to basis blade e23
    T c10{}; // bivector,     maps to basis blade e31
    T c11{}; // bivector,     maps to basis blade e12
    T c12{}; // bivector,     maps to basis blade e15
    T c13{}; // bivector,     maps to basis blade e25
    T c14{}; // bivector,     maps to basis blade e35
    T c15{}; // bivector,     maps to basis blade e45
    T c16{}; // trivector,    maps to basis blade e415
    T c17{}; // trivector,    maps to basis blade e425
    T c18{}; // trivector,    maps to basis blade e435
    T c19{}; // trivector,    maps to basis blade e235
    T c20{}; // trivector,    maps to basis blade e315
    T c21{}; // trivector,    maps to basis blade e125
    T c22{}; // trivector,    maps to basis blade e423
    T c23{}; // trivector,    maps to basis blade e431
    T c24{}; // trivector,    maps to basis blade e412
    T c25{}; // trivector,    maps to basis blade e321
    T c26{}; // quadvector,   maps to basis blade e4235
    T c27{}; // quadvector,   maps to basis blade e4315
    T c28{}; // quadvector,   maps to basis blade e4125
    T c29{}; // quadvector,   maps to basis blade e3215
    T c30{}; // quadvector,   maps to basis blade e1234
    T c31{}; // pseudoscalar, maps to basis blade e12345

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator+=(MVec32_t<U, Tag> const& v) noexcept
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
        c16 += v.c16;
        c17 += v.c17;
        c18 += v.c18;
        c19 += v.c19;
        c20 += v.c20;
        c21 += v.c21;
        c22 += v.c22;
        c23 += v.c23;
        c24 += v.c24;
        c25 += v.c25;
        c26 += v.c26;
        c27 += v.c27;
        c28 += v.c28;
        c29 += v.c29;
        c30 += v.c30;
        c31 += v.c31;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator-=(MVec32_t<U, Tag> const& v) noexcept
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
        c16 -= v.c16;
        c17 -= v.c17;
        c18 -= v.c18;
        c19 -= v.c19;
        c20 -= v.c20;
        c21 -= v.c21;
        c22 -= v.c22;
        c23 -= v.c23;
        c24 -= v.c24;
        c25 -= v.c25;
        c26 -= v.c26;
        c27 -= v.c27;
        c28 -= v.c28;
        c29 -= v.c29;
        c30 -= v.c30;
        c31 -= v.c31;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator*=(U s) noexcept
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
        c16 *= s;
        c17 *= s;
        c18 *= s;
        c19 *= s;
        c20 *= s;
        c21 *= s;
        c22 *= s;
        c23 *= s;
        c24 *= s;
        c25 *= s;
        c26 *= s;
        c27 *= s;
        c28 *= s;
        c29 *= s;
        c30 *= s;
        c31 *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 32 comp.");
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
        c16 /= s;
        c17 /= s;
        c18 /= s;
        c19 /= s;
        c20 /= s;
        c21 /= s;
        c22 /= s;
        c23 /= s;
        c24 /= s;
        c25 /= s;
        c26 /= s;
        c27 /= s;
        c28 /= s;
        c29 /= s;
        c30 /= s;
        c31 /= s;
        return (*this);
    }
};


////////////////////////////////////////////////////////////////////////////////
// MVec32_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(MVec32_t<T, Tag> const& lhs, MVec32_t<U, Tag> const& rhs)
{
    // componentwise comparison
    // comparison is not exact, but accepts epsilon deviations
    auto constexpr delta_eps = detail::safe_epsilon<T, U>();
    return (
        std::abs(lhs.c0 - rhs.c0) < delta_eps && std::abs(lhs.c1 - rhs.c1) < delta_eps &&
        std::abs(lhs.c2 - rhs.c2) < delta_eps && std::abs(lhs.c3 - rhs.c3) < delta_eps &&
        std::abs(lhs.c4 - rhs.c4) < delta_eps && std::abs(lhs.c5 - rhs.c5) < delta_eps &&
        std::abs(lhs.c6 - rhs.c6) < delta_eps && std::abs(lhs.c7 - rhs.c7) < delta_eps &&
        std::abs(lhs.c8 - rhs.c8) < delta_eps && std::abs(lhs.c9 - rhs.c9) < delta_eps &&
        std::abs(lhs.c10 - rhs.c10) < delta_eps &&
        std::abs(lhs.c11 - rhs.c11) < delta_eps &&
        std::abs(lhs.c12 - rhs.c12) < delta_eps &&
        std::abs(lhs.c13 - rhs.c13) < delta_eps &&
        std::abs(lhs.c14 - rhs.c14) < delta_eps &&
        std::abs(lhs.c15 - rhs.c15) < delta_eps &&
        std::abs(lhs.c16 - rhs.c16) < delta_eps &&
        std::abs(lhs.c17 - rhs.c17) < delta_eps &&
        std::abs(lhs.c18 - rhs.c18) < delta_eps &&
        std::abs(lhs.c19 - rhs.c19) < delta_eps &&
        std::abs(lhs.c20 - rhs.c20) < delta_eps &&
        std::abs(lhs.c21 - rhs.c21) < delta_eps &&
        std::abs(lhs.c22 - rhs.c22) < delta_eps &&
        std::abs(lhs.c23 - rhs.c23) < delta_eps &&
        std::abs(lhs.c24 - rhs.c24) < delta_eps &&
        std::abs(lhs.c25 - rhs.c25) < delta_eps &&
        std::abs(lhs.c26 - rhs.c26) < delta_eps &&
        std::abs(lhs.c27 - rhs.c27) < delta_eps &&
        std::abs(lhs.c28 - rhs.c28) < delta_eps &&
        std::abs(lhs.c29 - rhs.c29) < delta_eps &&
        std::abs(lhs.c30 - rhs.c30) < delta_eps &&
        std::abs(lhs.c31 - rhs.c31) < delta_eps);
}

// inequality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(MVec32_t<T, Tag> const& lhs, MVec32_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr MVec32_t<T, Tag> operator-(MVec32_t<T, Tag> const& v)
{
    return MVec32_t<T, Tag>(-v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, -v.c7, -v.c8,
                            -v.c9, -v.c10, -v.c11, -v.c12, -v.c13, -v.c14, -v.c15, -v.c16,
                            -v.c17, -v.c18, -v.c19, -v.c20, -v.c21, -v.c22, -v.c23,
                            -v.c24, -v.c25, -v.c26, -v.c27, -v.c28, -v.c29, -v.c30,
                            -v.c31);
}

// adding multivectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec32_t<std::common_type_t<T, U>, Tag> operator+(MVec32_t<T, Tag> const& v1,
                                                            MVec32_t<U, Tag> const& v2)
{
    return MVec32_t<std::common_type_t<T, U>, Tag>(
        v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2, v1.c3 + v2.c3, v1.c4 + v2.c4,
        v1.c5 + v2.c5, v1.c6 + v2.c6, v1.c7 + v2.c7, v1.c8 + v2.c8, v1.c9 + v2.c9,
        v1.c10 + v2.c10, v1.c11 + v2.c11, v1.c12 + v2.c12, v1.c13 + v2.c13,
        v1.c14 + v2.c14, v1.c15 + v2.c15, v1.c16 + v2.c16, v1.c17 + v2.c17,
        v1.c18 + v2.c18, v1.c19 + v2.c19, v1.c20 + v2.c20, v1.c21 + v2.c21,
        v1.c22 + v2.c22, v1.c23 + v2.c23, v1.c24 + v2.c24, v1.c25 + v2.c25,
        v1.c26 + v2.c26, v1.c27 + v2.c27, v1.c28 + v2.c28, v1.c29 + v2.c29,
        v1.c30 + v2.c30, v1.c31 + v2.c31);
}

// substracting multivectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec32_t<std::common_type_t<T, U>, Tag> operator-(MVec32_t<T, Tag> const& v1,
                                                            MVec32_t<U, Tag> const& v2)
{
    return MVec32_t<std::common_type_t<T, U>, Tag>(
        v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2, v1.c3 - v2.c3, v1.c4 - v2.c4,
        v1.c5 - v2.c5, v1.c6 - v2.c6, v1.c7 - v2.c7, v1.c8 - v2.c8, v1.c9 - v2.c9,
        v1.c10 - v2.c10, v1.c11 - v2.c11, v1.c12 - v2.c12, v1.c13 - v2.c13,
        v1.c14 - v2.c14, v1.c15 - v2.c15, v1.c16 - v2.c16, v1.c17 - v2.c17,
        v1.c18 - v2.c18, v1.c19 - v2.c19, v1.c20 - v2.c20, v1.c21 - v2.c21,
        v1.c22 - v2.c22, v1.c23 - v2.c23, v1.c24 - v2.c24, v1.c25 - v2.c25,
        v1.c26 - v2.c26, v1.c27 - v2.c27, v1.c28 - v2.c28, v1.c29 - v2.c29,
        v1.c30 - v2.c30, v1.c31 - v2.c31);
}

// multiply with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec32_t<std::common_type_t<T, U>, Tag> operator*(MVec32_t<T, Tag> const& v,
                                                            U s)
{
    return MVec32_t<std::common_type_t<T, U>, Tag>(
        v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s, v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s,
        v.c8 * s, v.c9 * s, v.c10 * s, v.c11 * s, v.c12 * s, v.c13 * s, v.c14 * s,
        v.c15 * s, v.c16 * s, v.c17 * s, v.c18 * s, v.c19 * s, v.c20 * s, v.c21 * s,
        v.c22 * s, v.c23 * s, v.c24 * s, v.c25 * s, v.c26 * s, v.c27 * s, v.c28 * s,
        v.c29 * s, v.c30 * s, v.c31 * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec32_t<std::common_type_t<T, U>, Tag> operator*(T s,
                                                            MVec32_t<U, Tag> const& v)
{
    return MVec32_t<std::common_type_t<T, U>, Tag>(
        v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s, v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s,
        v.c8 * s, v.c9 * s, v.c10 * s, v.c11 * s, v.c12 * s, v.c13 * s, v.c14 * s,
        v.c15 * s, v.c16 * s, v.c17 * s, v.c18 * s, v.c19 * s, v.c20 * s, v.c21 * s,
        v.c22 * s, v.c23 * s, v.c24 * s, v.c25 * s, v.c26 * s, v.c27 * s, v.c28 * s,
        v.c29 * s, v.c30 * s, v.c31 * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline MVec32_t<std::common_type_t<T, U>, Tag> operator/(MVec32_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "multivector division 32 comp.");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec32_t<ctype, Tag>(
        v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv, v.c4 * inv, v.c5 * inv,
        v.c6 * inv, v.c7 * inv, v.c8 * inv, v.c9 * inv, v.c10 * inv, v.c11 * inv,
        v.c12 * inv, v.c13 * inv, v.c14 * inv, v.c15 * inv, v.c16 * inv, v.c17 * inv,
        v.c18 * inv, v.c19 * inv, v.c20 * inv, v.c21 * inv, v.c22 * inv, v.c23 * inv,
        v.c24 * inv, v.c25 * inv, v.c26 * inv, v.c27 * inv, v.c28 * inv, v.c29 * inv,
        v.c30 * inv, v.c31 * inv);
}

////////////////////////////////////////////////////////////////////////////////
// MVec32_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, MVec32_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << "," << v.c4 << ","
       << v.c5 << "," << v.c6 << "," << v.c7 << "," << v.c8 << "," << v.c9 << "," << v.c10
       << "," << v.c11 << "," << v.c12 << "," << v.c13 << "," << v.c14 << "," << v.c15
       << "," << v.c16 << "," << v.c17 << "," << v.c18 << "," << v.c19 << "," << v.c20
       << "," << v.c21 << "," << v.c22 << "," << v.c23 << "," << v.c24 << "," << v.c25
       << "," << v.c26 << "," << v.c27 << "," << v.c28 << "," << v.c29 << "," << v.c30
       << "," << v.c31 << ")";
    return os;
}

} // namespace hd::ga
