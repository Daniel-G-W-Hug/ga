#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::ostream

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T, Tag> definition of a multivector with 4 components c0, ..., c3
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t, typename Tag = default_tag>
    requires(std::floating_point<T>)
struct MVec4_t {

    // ctors

    // (all grades = 0)
    MVec4_t() = default;

    // assign all components
    constexpr MVec4_t(T s, T x, T y, T ps) : c0(s), c1(x), c2(y), c3(ps) {}

    // cp & mv ctor
    MVec4_t(MVec4_t const&) = default;
    MVec4_t(MVec4_t&&) noexcept = default;

    // cp & mv assignment
    MVec4_t& operator=(MVec4_t const&) = default;
    MVec4_t& operator=(MVec4_t&&) noexcept = default;

    // dtor
    ~MVec4_t() = default;

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    constexpr MVec4_t(MVec4_t<U, Tag> const& v) : MVec4_t(v.c0, v.c1, v.c2, v.c3)
    {
    }

    friend void swap(MVec4_t& lhs, MVec4_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.c0, rhs.c0);
        swap(lhs.c1, rhs.c1);
        swap(lhs.c2, rhs.c2);
        swap(lhs.c3, rhs.c3);
    }

    ////////////////////////////////////////////////////////////////////////////
    // component definition
    ////////////////////////////////////////////////////////////////////////////

    // 2d multivector:   // 3d even multivector:   // 3d odd mulitvector:
    T c0{}; // 1 scalar  // 1 scalar               // e1
    T c1{}; // e1        // e2^e3 = e23 = yz       // e2
    T c2{}; // e2        // e3^e1 = e31 = zy       // e3
    T c3{}; // e12       // e1^e2 = e12 = xy       // e123 (pseudoscalar)

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec4_t<U, Tag> const& rhs) const
    {
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto abs_delta_c2 = std::abs(rhs.c2 - c2);
        auto abs_delta_c3 = std::abs(rhs.c3 - c3);
        auto constexpr delta_eps = detail::safe_epsilon<T, U>();
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps)
            return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec4_t& operator+=(MVec4_t<U, Tag> const& v) noexcept
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec4_t& operator-=(MVec4_t<U, Tag> const& v) noexcept
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec4_t& operator*=(U s) noexcept
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec4_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 4 comp.");
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr MVec4_t<T, Tag> operator-(MVec4_t<T, Tag> const& v)
{
    return MVec4_t<T, Tag>(-v.c0, -v.c1, -v.c2, -v.c3);
}

// add multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4_t<std::common_type_t<T, U>, Tag>
operator+(MVec4_t<T, Tag> const& v1, MVec4_t<U, Tag> const& v2)
{
    return MVec4_t<std::common_type_t<T, U>, Tag>(v1.c0 + v2.c0, v1.c1 + v2.c1,
                                                  v1.c2 + v2.c2, v1.c3 + v2.c3);
}

// substract multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4_t<std::common_type_t<T, U>, Tag>
operator-(MVec4_t<T, Tag> const& v1, MVec4_t<U, Tag> const& v2)
{
    return MVec4_t<std::common_type_t<T, U>, Tag>(v1.c0 - v2.c0, v1.c1 - v2.c1,
                                                  v1.c2 - v2.c2, v1.c3 - v2.c3);
}

// multiply a multivector with a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4_t<std::common_type_t<T, U>, Tag>
operator*(MVec4_t<T, Tag> const& v, U s)
{
    return MVec4_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec4_t<std::common_type_t<T, U>, Tag>
operator*(T s, MVec4_t<U, Tag> const& v)
{
    return MVec4_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline MVec4_t<std::common_type_t<T, U>, Tag> operator/(MVec4_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "multivector division");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplication with inverse value
    return MVec4_t<ctype, Tag>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv);
}

// magnitude of a k-vector (in representational space)
// |Z|^2 = Z rev(Z) = c0^2 + c1^2 + c2^2 + c3^2
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm_sq(MVec4_t<T, Tag> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3;
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(MVec4_t<T, Tag> const& v)
{
    return std::sqrt(nrm_sq(v));
}

////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec4_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << ")";
    return os;
}

} // namespace hd::ga