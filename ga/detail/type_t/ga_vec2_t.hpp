#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <array>    // std::array (componentwise comparison)
#include <concepts> // numeric_type<T>
#include <ostream>  // std::ostream

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T, Tag> definition
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(numeric_type<T>)
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
        requires(numeric_type<U>)
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

    template <typename U>
        requires(numeric_type<U>)
    Vec2_t& operator+=(Vec2_t<U, Tag> const& v) noexcept
    {
        x += v.x;
        y += v.y;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec2_t& operator-=(Vec2_t<U, Tag> const& v) noexcept
    {
        x -= v.x;
        y -= v.y;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec2_t& operator*=(U s) noexcept
    {
        x *= s;
        y *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    Vec2_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "vector division 2 comp.");
        x /= s;
        y /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(Vec2_t<T, Tag> const& lhs, Vec2_t<U, Tag> const& rhs)
{
    // componentwise comparison against the library's one comparison rule
    // (detail::coeffs_equal, in ga_error_handling.hpp -- see the note there)
    return detail::coeffs_equal(std::array{lhs.x, lhs.y}, std::array{rhs.x, rhs.y});
}

// inequality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(Vec2_t<T, Tag> const& lhs, Vec2_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr Vec2_t<T, Tag> operator-(Vec2_t<T, Tag> const& v)
{
    return Vec2_t<T, Tag>(-v.x, -v.y);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator+(Vec2_t<T, Tag> const& v1,
                                                          Vec2_t<U, Tag> const& v2)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator-(Vec2_t<T, Tag> const& v1,
                                                          Vec2_t<U, Tag> const& v2)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator*(Vec2_t<T, Tag> const& v, U s)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2_t<std::common_type_t<T, U>, Tag> operator*(T s, Vec2_t<U, Tag> const& v)
{
    return Vec2_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline Vec2_t<std::common_type_t<T, U>, Tag> operator/(Vec2_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "vector division 2 comp.");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec2_t<std::common_type_t<T, U>, Tag>(v.x * inv, v.y * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, Vec2_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << ")";
    return os;
}


////////////////////////////////////////////////////////////////////////////////
// coeff_sq: the sum of the squared coefficients -- a METRIC-FREE gauge of "how big is
// this object", used to judge degeneracy relative to the object itself (see
// detail::check_invertible in ga_error_handling.hpp). It is NOT a norm: no metric enters,
// which is the point -- nrm_sq vanishes on a null blade in a Lorentzian or conformal
// metric and on an ideal one in a degenerate metric, exactly where a gauge is needed.
////////////////////////////////////////////////////////////////////////////////

namespace detail {

template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T coeff_sq(Vec2_t<T, Tag> const& v)
{
    return v.x * v.x + v.y * v.y;
}

} // namespace detail

} // namespace hd::ga