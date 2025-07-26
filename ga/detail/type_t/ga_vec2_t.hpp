#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <algorithm> // std::max
#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T, Tag> definition
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(std::floating_point<T>)
struct Vec2_t {

    // assumes a right-handed orthonormal vector basis {e1, e2} in EGA
    // using components {x, y}, such that for each vector v:
    // v = x * e1 + y * e2

    // ctors
    Vec2_t() = default;

    constexpr Vec2_t(T x_in, T y_in) : x(x_in), y(y_in) {}

    // cp & mv ctor
    Vec2_t(Vec2_t const&) = default;
    Vec2_t(Vec2_t&&) noexcept = default;

    // cp & mv assignment
    Vec2_t& operator=(Vec2_t const&) = default;
    Vec2_t& operator=(Vec2_t&&) noexcept = default;

    // dtor
    ~Vec2_t() = default;

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    constexpr Vec2_t(Vec2_t<U, Tag> const& v) : x(v.x), y(v.y)
    {
    }

    friend void swap(Vec2_t& lhs, Vec2_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.x, rhs.x);
        swap(lhs.y, rhs.y);
    }

    T x{};
    T y{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(Vec2_t<U, Tag> const& rhs) const
    {
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_x = std::abs(x - rhs.x);
        auto abs_delta_y = std::abs(y - rhs.y);
        auto constexpr delta_eps = detail::safe_epsilon<T, U>();
        return (abs_delta_x < delta_eps && abs_delta_y < delta_eps);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec2_t& operator+=(Vec2_t<U, Tag> const& v) noexcept
    {
        x += v.x;
        y += v.y;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec2_t& operator-=(Vec2_t<U, Tag> const& v) noexcept
    {
        x -= v.x;
        y -= v.y;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec2_t& operator*=(U s) noexcept
    {
        x *= s;
        y *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec2_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_normalization<U>(std::abs(s), "vector division 2 comp.");
        x /= s;
        y /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr Vec2_t<T, Tag> operator-(Vec2_t<T, Tag> const& v)
{
    return Vec2_t<T, Tag>(-v.x, -v.y);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator+(Vec2_t<T, Tag> const& v1,
                                                          Vec2_t<U, Tag> const& v2)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator-(Vec2_t<T, Tag> const& v1,
                                                          Vec2_t<U, Tag> const& v2)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator*(Vec2_t<T, Tag> const& v, U s)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator*(T s, Vec2_t<U, Tag> const& v)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec2_t<std::common_type_t<T, U>, Tag> operator/(Vec2_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "vector division 2 comp.");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec2_t<std::common_type_t<T, U>, Tag>(v.x * inv, v.y * inv);
}

// magnitude of the k-vector (in representational space)
//
// magnitude is always defined in the representational space, i.e. without
// covering the target metric of the representation
//
// -> in ega magnitude (magn) and norm (nrm) are identical
//
// -> in pga magnitude is defined as norm for the representational space
//    assuming the corresponding identity matrix as metric for that space:
//
//    G<3,0,0> as representational space for modelling 2d Euclidean space using G<2,0,1>
//             representational space:        e1^2=+1, e2^2=+1, e3^2=+1
//             modelled 2d Euclidean space:   e1^2=+1, e2^2=+1, e3^2= 0
//
//    G<4,0,0> as representational space for modelling 3d Euclidean space using G<3,0,1>
//             representational space:        e1^2=+1, e2^2=+1, e3^2=+1, e4^2=+1
//             modelled 3d Euclidean space:   e1^2=+1, e2^2=+1, e3^2=+1, e4^2= 0
//
//    norm is always defined for the modelled space incl. the target metric,
//    by using the exomorphisim matrix G as defined by Lengyel in the book
//    "Projective geometric algebra illuminated"
//
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T nrm_sq(Vec2_t<T, Tag> const& v)
{
    // implements the scalar product as defined by the geometric product *
    // |v|^2 = gr0( v*v ) = gr0( dot(v,v) + wdg(v,v) ) = dot(v,v)
    return v.x * v.x + v.y * v.y;
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T nrm(Vec2_t<T, Tag> const& v)
{
    return std::sqrt(nrm_sq(v));
}

// return a vector v normalized to nrm(v) == 1.0
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline Vec2_t<T, Tag> normalize(Vec2_t<T, Tag> const& v)
{
    T m = nrm(v);
    detail::check_normalization<T>(m, "vector");
    T inv = T(1.0) / m; // for multiplication with inverse of norm
    return Vec2_t<T, Tag>(v.x * inv, v.y * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Vec2_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << ")";
    return os;
}

} // namespace hd::ga