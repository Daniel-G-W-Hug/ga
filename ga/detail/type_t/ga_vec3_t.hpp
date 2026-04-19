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
// Vec3_t<T, Tag> definition
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(numeric_type<T>)
struct Vec3_t {

    // assumes a right-handed orthonormal vector basis {e1, e2, e3}
    // using components {x, y, z}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3
    //
    // and for each bivector bv:
    // b = yz * e2^e3 + zx * e3^e1 + xy * e1^e2
    //   =  x * e2^e3 +  y * e3^e1 +  z * e1^e2
    // (same names like in vector (x, y, z), but semantic as one line above (yz, zx, xy)

    // this is a mapping of the components
    // such that vector components x, y and z correspond to the
    // normals of the corresponding plane elements represented by
    // bivector components yz, zx and xy
    // i.e. they can by converted to each other by a duality transformation
    //
    // T.x <=> bivector yz <=> Vec3_t<T>::x; // maps to basis bivector e2^e3
    // T.y <=> bivector zx <=> Vec3_t<T>::y; // maps to basis bivector e3^e1
    // T.z <=> bivector xy <=> Vec3_t<T>::z; // maps to basis bivector e1^e2

    // ctors
    Vec3_t() = default;

    constexpr Vec3_t(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {}

    // cp & mv ctor
    Vec3_t(Vec3_t const&) = default;
    Vec3_t(Vec3_t&&) noexcept = default;

    // cp & mv assignment
    Vec3_t& operator=(Vec3_t const&) = default;
    Vec3_t& operator=(Vec3_t&&) noexcept = default;

    // dtor
    ~Vec3_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr Vec3_t(Vec3_t<U, Tag> const& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    friend void swap(Vec3_t& lhs, Vec3_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.x, rhs.x);
        swap(lhs.y, rhs.y);
        swap(lhs.z, rhs.z);
    }

    T x{};
    T y{};
    T z{};

    template <typename U>
        requires(numeric_type<U>)
    Vec3_t& operator+=(Vec3_t<U, Tag> const& v) noexcept
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec3_t& operator-=(Vec3_t<U, Tag> const& v) noexcept
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec3_t& operator*=(U s) noexcept
    {
        x *= s;
        y *= s;
        z *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec3_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "vector division 3 comp.");
        x /= s;
        y /= s;
        z /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(Vec3_t<T, Tag> const& lhs, Vec3_t<U, Tag> const& rhs)
{
    // componentwise comparison
    // equality implies same magnitude and direction
    // comparison is not exact, but accepts epsilon deviations
    auto abs_delta_x = std::abs(lhs.x - rhs.x);
    auto abs_delta_y = std::abs(lhs.y - rhs.y);
    auto abs_delta_z = std::abs(lhs.z - rhs.z);
    auto constexpr delta_eps = detail::safe_epsilon<T, U>();
    return (abs_delta_x < delta_eps && abs_delta_y < delta_eps &&
            abs_delta_z < delta_eps);
}

// inequality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(Vec3_t<T, Tag> const& lhs, Vec3_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr Vec3_t<T, Tag> operator-(Vec3_t<T, Tag> const& v)
{
    return Vec3_t<T, Tag>(-v.x, -v.y, -v.z);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator+(Vec3_t<T, Tag> const& v1,
                                                          Vec3_t<U, Tag> const& v2)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator-(Vec3_t<T, Tag> const& v1,
                                                          Vec3_t<U, Tag> const& v2)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator*(Vec3_t<T, Tag> const& v, U s)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator*(T s, Vec3_t<U, Tag> const& v)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline Vec3_t<std::common_type_t<T, U>, Tag> operator/(Vec3_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "vector division 3 comp.");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec3_t<ctype, Tag>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, Vec3_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

} // namespace hd::ga