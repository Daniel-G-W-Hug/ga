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
// Vec4_t<T, Tag> definition
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(numeric_type<T>)
struct Vec4_t {

    // assumes a right-handed orthonormal vector basis {e1, e2, e3, e4}
    // using components {x, y, z, w}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3 + w * e4
    //
    // and for each trivector tv:
    // t = wyz * e4^e2^e3 + wzx * e4^e3^e1 + wxy * e4^e1^e2 + zyx * e3^e2^e1
    //   = t_x * e4^e2^e3 + t_y * e4^e3^e1 + t_z * e4^e1^e2 + t_w * e3^e2^e1
    // (same names like in vector (x, y, z, w),
    //  but semantic as above (wyz, wzx, wxy, zyx)

    // this is a mapping of the components
    // such that vector components x, y, z and w correspond to the
    // trivector elements represented by wyz, wzx, wxy and zyx
    // i.e. they can by converted to each other by a duality transformation
    //
    // T.x <=> trivector wyz <=> Vec4_t<T>::x; // maps to basis trivector e4^e2^e3
    // T.y <=> trivector wzx <=> Vec4_t<T>::y; // maps to basis trivector e4^e3^e1
    // T.z <=> trivector wxy <=> Vec4_t<T>::z; // maps to basis trivector e4^e1^e2
    // T.w <=> trivector zyx <=> Vec4_t<T>::w; // maps to basis trivector e3^e2^e1

    // ctors
    Vec4_t() = default;

    constexpr Vec4_t(T x_in, T y_in, T z_in, T w_in) : x(x_in), y(y_in), z(z_in), w(w_in)
    {
    }

    // cp & mv ctor
    Vec4_t(Vec4_t const&) = default;
    Vec4_t(Vec4_t&&) noexcept = default;

    // cp & mv assignment
    Vec4_t& operator=(Vec4_t const&) = default;
    Vec4_t& operator=(Vec4_t&&) noexcept = default;

    // dtor
    ~Vec4_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr Vec4_t(Vec4_t<U, Tag> const& v) : x(v.x), y(v.y), z(v.z), w(v.w)
    {
    }

    friend void swap(Vec4_t& lhs, Vec4_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.x, rhs.x);
        swap(lhs.y, rhs.y);
        swap(lhs.z, rhs.z);
        swap(lhs.w, rhs.w);
    }

    T x{};
    T y{};
    T z{};
    T w{};

    template <typename U>
        requires(numeric_type<U>)
    Vec4_t& operator+=(Vec4_t<U, Tag> const& v) noexcept
    {
        x += v.x;
        y += v.y;
        z += v.z;
        w += v.w;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec4_t& operator-=(Vec4_t<U, Tag> const& v) noexcept
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        w -= v.w;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec4_t& operator*=(U s) noexcept
    {
        x *= s;
        y *= s;
        z *= s;
        w *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec4_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "vector division 4 comp.");
        x /= s;
        y /= s;
        z /= s;
        w /= s;
        return (*this);
    }
};


////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(Vec4_t<T, Tag> const& lhs, Vec4_t<U, Tag> const& rhs)
{
    // componentwise comparison
    // equality implies same magnitude and direction
    // comparison is not exact, but accepts epsilon deviations
    auto abs_delta_x = std::abs(lhs.x - rhs.x);
    auto abs_delta_y = std::abs(lhs.y - rhs.y);
    auto abs_delta_z = std::abs(lhs.z - rhs.z);
    auto abs_delta_w = std::abs(lhs.w - rhs.w);
    auto constexpr delta_eps = detail::safe_epsilon<T, U>();
    return (abs_delta_x < delta_eps && abs_delta_y < delta_eps &&
            abs_delta_z < delta_eps && abs_delta_w < delta_eps);
}

// inequality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(Vec4_t<T, Tag> const& lhs, Vec4_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr Vec4_t<T, Tag> operator-(Vec4_t<T, Tag> const& v)
{
    return Vec4_t<T, Tag>(-v.x, -v.y, -v.z, -v.w);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator+(Vec4_t<T, Tag> const& v1,
                                                          Vec4_t<U, Tag> const& v2)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z,
                                                 v1.w + v2.w);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator-(Vec4_t<T, Tag> const& v1,
                                                          Vec4_t<U, Tag> const& v2)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z,
                                                 v1.w - v2.w);
}

// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator*(Vec4_t<T, Tag> const& v, U s)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s, v.w * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator*(T s, Vec4_t<U, Tag> const& v)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s, v.w * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline Vec4_t<std::common_type_t<T, U>, Tag> operator/(Vec4_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "vector division 4 comp.");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec4_t<ctype, Tag>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, Vec4_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << "," << v.w << ")";
    return os;
}

} // namespace hd::ga