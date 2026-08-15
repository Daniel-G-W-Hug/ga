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
// Vec5_t<T, Tag> definition
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(numeric_type<T>)
struct Vec5_t {

    // assumes a vector basis {e1, e2, e3, e4, e5}
    // using components {x, y, z, w, u}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3 + w * e4 + u * e5
    //
    // quadvectors share this template with POSITIONAL field reuse (dossier
    // convention): the five basis quadvectors are stored in the same fields
    // in the label order x, y, z, w, u given by the omitted basis vector,
    // e.g. for cga3dc:
    //
    // Q.x <=> quadvector component on e4235  (omits e1)
    // Q.y <=> quadvector component on e4315  (omits e2)
    // Q.z <=> quadvector component on e4125  (omits e3)
    // Q.w <=> quadvector component on e3215  (omits e4)
    // Q.u <=> quadvector component on e1234  (omits e5)

    // ctors
    Vec5_t() = default;

    constexpr Vec5_t(T x_in, T y_in, T z_in, T w_in, T u_in) :
        x(x_in), y(y_in), z(z_in), w(w_in), u(u_in)
    {
    }

    // cp & mv ctor
    Vec5_t(Vec5_t const&) = default;
    Vec5_t(Vec5_t&&) noexcept = default;

    // cp & mv assignment
    Vec5_t& operator=(Vec5_t const&) = default;
    Vec5_t& operator=(Vec5_t&&) noexcept = default;

    // dtor
    ~Vec5_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr Vec5_t(Vec5_t<U, Tag> const& v) : x(v.x), y(v.y), z(v.z), w(v.w), u(v.u)
    {
    }

    friend void swap(Vec5_t& lhs, Vec5_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.x, rhs.x);
        swap(lhs.y, rhs.y);
        swap(lhs.z, rhs.z);
        swap(lhs.w, rhs.w);
        swap(lhs.u, rhs.u);
    }

    T x{};
    T y{};
    T z{};
    T w{};
    T u{};

    template <typename U>
        requires(numeric_type<U>)
    Vec5_t& operator+=(Vec5_t<U, Tag> const& v) noexcept
    {
        x += v.x;
        y += v.y;
        z += v.z;
        w += v.w;
        u += v.u;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec5_t& operator-=(Vec5_t<U, Tag> const& v) noexcept
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        w -= v.w;
        u -= v.u;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec5_t& operator*=(U s) noexcept
    {
        x *= s;
        y *= s;
        z *= s;
        w *= s;
        u *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec5_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "vector division 5 comp.");
        x /= s;
        y /= s;
        z /= s;
        w /= s;
        u /= s;
        return (*this);
    }
};


////////////////////////////////////////////////////////////////////////////////
// Vec5_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(Vec5_t<T, Tag> const& lhs, Vec5_t<U, Tag> const& rhs)
{
    // componentwise comparison
    // equality implies same magnitude and direction
    // comparison is not exact, but accepts epsilon deviations
    auto abs_delta_x = std::abs(lhs.x - rhs.x);
    auto abs_delta_y = std::abs(lhs.y - rhs.y);
    auto abs_delta_z = std::abs(lhs.z - rhs.z);
    auto abs_delta_w = std::abs(lhs.w - rhs.w);
    auto abs_delta_u = std::abs(lhs.u - rhs.u);
    auto constexpr delta_eps = detail::safe_epsilon<T, U>();
    return (abs_delta_x < delta_eps && abs_delta_y < delta_eps &&
            abs_delta_z < delta_eps && abs_delta_w < delta_eps &&
            abs_delta_u < delta_eps);
}

// inequality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(Vec5_t<T, Tag> const& lhs, Vec5_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr Vec5_t<T, Tag> operator-(Vec5_t<T, Tag> const& v)
{
    return Vec5_t<T, Tag>(-v.x, -v.y, -v.z, -v.w, -v.u);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec5_t<std::common_type_t<T, U>, Tag> operator+(Vec5_t<T, Tag> const& v1,
                                                          Vec5_t<U, Tag> const& v2)
{
    return Vec5_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z,
                                                 v1.w + v2.w, v1.u + v2.u);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec5_t<std::common_type_t<T, U>, Tag> operator-(Vec5_t<T, Tag> const& v1,
                                                          Vec5_t<U, Tag> const& v2)
{
    return Vec5_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z,
                                                 v1.w - v2.w, v1.u - v2.u);
}

// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec5_t<std::common_type_t<T, U>, Tag> operator*(Vec5_t<T, Tag> const& v, U s)
{
    return Vec5_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s, v.w * s,
                                                 v.u * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec5_t<std::common_type_t<T, U>, Tag> operator*(T s, Vec5_t<U, Tag> const& v)
{
    return Vec5_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s, v.w * s,
                                                 v.u * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline Vec5_t<std::common_type_t<T, U>, Tag> operator/(Vec5_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "vector division 5 comp.");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec5_t<ctype, Tag>(v.x * inv, v.y * inv, v.z * inv, v.w * inv, v.u * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec5_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, Vec5_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << "," << v.w << "," << v.u << ")";
    return os;
}

} // namespace hd::ga
