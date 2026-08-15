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
// BVec10_t<T, Tag> definition of a 10-component type with fields
// vx, vy, vz, mx, my, mz, px, py, pz, pw
//
// used for the grade-2 AND grade-3 elements of 5-dimensional algebras (shared
// via tags); the field names follow the conformal object structure (dossier
// convention): v/m carry the flat-line part, p the position part, pw the
// pure-Euclidean component.
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires numeric_type<T>
struct BVec10_t {

    // ctors

    // (all grades = 0)
    BVec10_t() = default;

    // assign all components directly
    constexpr BVec10_t(T vx_in, T vy_in, T vz_in, T mx_in, T my_in, T mz_in, T px_in,
                       T py_in, T pz_in, T pw_in) :
        vx(vx_in), vy(vy_in), vz(vz_in), mx(mx_in), my(my_in), mz(mz_in), px(px_in),
        py(py_in), pz(pz_in), pw(pw_in)
    {
    }

    // cp & mv ctor
    BVec10_t(BVec10_t const&) = default;
    BVec10_t(BVec10_t&&) noexcept = default;

    // cp & mv assignment
    BVec10_t& operator=(BVec10_t const&) = default;
    BVec10_t& operator=(BVec10_t&&) noexcept = default;

    // dtor
    ~BVec10_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr BVec10_t(BVec10_t<U, Tag> const& B) :
        vx(B.vx), vy(B.vy), vz(B.vz), mx(B.mx), my(B.my), mz(B.mz), px(B.px), py(B.py),
        pz(B.pz), pw(B.pw)
    {
    }

    friend void swap(BVec10_t& lhs, BVec10_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.vx, rhs.vx);
        swap(lhs.vy, rhs.vy);
        swap(lhs.vz, rhs.vz);
        swap(lhs.mx, rhs.mx);
        swap(lhs.my, rhs.my);
        swap(lhs.mz, rhs.mz);
        swap(lhs.px, rhs.px);
        swap(lhs.py, rhs.py);
        swap(lhs.pz, rhs.pz);
        swap(lhs.pw, rhs.pw);
    }

    T vx{};
    T vy{};
    T vz{};
    T mx{};
    T my{};
    T mz{};
    T px{};
    T py{};
    T pz{};
    T pw{};

    template <typename U>
        requires(numeric_type<U>)
    BVec10_t& operator+=(BVec10_t<U, Tag> const& B) noexcept
    {
        vx += B.vx;
        vy += B.vy;
        vz += B.vz;
        mx += B.mx;
        my += B.my;
        mz += B.mz;
        px += B.px;
        py += B.py;
        pz += B.pz;
        pw += B.pw;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    BVec10_t& operator-=(BVec10_t<U, Tag> const& B) noexcept
    {
        vx -= B.vx;
        vy -= B.vy;
        vz -= B.vz;
        mx -= B.mx;
        my -= B.my;
        mz -= B.mz;
        px -= B.px;
        py -= B.py;
        pz -= B.pz;
        pw -= B.pw;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    BVec10_t& operator*=(U s) noexcept
    {
        vx *= s;
        vy *= s;
        vz *= s;
        mx *= s;
        my *= s;
        mz *= s;
        px *= s;
        py *= s;
        pz *= s;
        pw *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    BVec10_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "division 10 comp.");
        vx /= s;
        vy /= s;
        vz /= s;
        mx /= s;
        my /= s;
        mz /= s;
        px /= s;
        py /= s;
        pz /= s;
        pw /= s;
        return (*this);
    }
};


////////////////////////////////////////////////////////////////////////////////
// BVec10_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(BVec10_t<T, Tag> const& lhs, BVec10_t<U, Tag> const& rhs)
{
    // componentwise comparison
    // comparison is not exact, but accepts epsilon deviations
    auto constexpr delta_eps = detail::safe_epsilon<T, U>();
    return (
        std::abs(lhs.vx - rhs.vx) < delta_eps && std::abs(lhs.vy - rhs.vy) < delta_eps &&
        std::abs(lhs.vz - rhs.vz) < delta_eps && std::abs(lhs.mx - rhs.mx) < delta_eps &&
        std::abs(lhs.my - rhs.my) < delta_eps && std::abs(lhs.mz - rhs.mz) < delta_eps &&
        std::abs(lhs.px - rhs.px) < delta_eps && std::abs(lhs.py - rhs.py) < delta_eps &&
        std::abs(lhs.pz - rhs.pz) < delta_eps && std::abs(lhs.pw - rhs.pw) < delta_eps);
}

// inequality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(BVec10_t<T, Tag> const& lhs, BVec10_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr BVec10_t<T, Tag> operator-(BVec10_t<T, Tag> const& B)
{
    return BVec10_t<T, Tag>(-B.vx, -B.vy, -B.vz, -B.mx, -B.my, -B.mz, -B.px, -B.py, -B.pz,
                            -B.pw);
}

// addition
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec10_t<std::common_type_t<T, U>, Tag> operator+(BVec10_t<T, Tag> const& B1,
                                                            BVec10_t<U, Tag> const& B2)
{
    return BVec10_t<std::common_type_t<T, U>, Tag>(
        B1.vx + B2.vx, B1.vy + B2.vy, B1.vz + B2.vz, B1.mx + B2.mx, B1.my + B2.my,
        B1.mz + B2.mz, B1.px + B2.px, B1.py + B2.py, B1.pz + B2.pz, B1.pw + B2.pw);
}

// subtraction
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec10_t<std::common_type_t<T, U>, Tag> operator-(BVec10_t<T, Tag> const& B1,
                                                            BVec10_t<U, Tag> const& B2)
{
    return BVec10_t<std::common_type_t<T, U>, Tag>(
        B1.vx - B2.vx, B1.vy - B2.vy, B1.vz - B2.vz, B1.mx - B2.mx, B1.my - B2.my,
        B1.mz - B2.mz, B1.px - B2.px, B1.py - B2.py, B1.pz - B2.pz, B1.pw - B2.pw);
}

// multiply with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec10_t<std::common_type_t<T, U>, Tag> operator*(BVec10_t<T, Tag> const& B,
                                                            U s)
{
    return BVec10_t<std::common_type_t<T, U>, Tag>(B.vx * s, B.vy * s, B.vz * s, B.mx * s,
                                                   B.my * s, B.mz * s, B.px * s, B.py * s,
                                                   B.pz * s, B.pw * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec10_t<std::common_type_t<T, U>, Tag> operator*(T s,
                                                            BVec10_t<U, Tag> const& B)
{
    return BVec10_t<std::common_type_t<T, U>, Tag>(B.vx * s, B.vy * s, B.vz * s, B.mx * s,
                                                   B.my * s, B.mz * s, B.px * s, B.py * s,
                                                   B.pz * s, B.pw * s);
}

// divide by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline BVec10_t<std::common_type_t<T, U>, Tag> operator/(BVec10_t<T, Tag> const& B, U s)
{
    detail::check_division_by_zero<T, U>(s, "division 10 comp.");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return BVec10_t<ctype, Tag>(B.vx * inv, B.vy * inv, B.vz * inv, B.mx * inv,
                                B.my * inv, B.mz * inv, B.px * inv, B.py * inv,
                                B.pz * inv, B.pw * inv);
}

////////////////////////////////////////////////////////////////////////////////
// BVec10_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, BVec10_t<T, Tag> const& B)
{
    os << "(" << B.vx << "," << B.vy << "," << B.vz << "," << B.mx << "," << B.my << ","
       << B.mz << "," << B.px << "," << B.py << "," << B.pz << "," << B.pw << ")";
    return os;
}

} // namespace hd::ga
