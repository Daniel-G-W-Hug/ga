#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <algorithm> // std::max
#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_type_tags.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> definition of a multivector with 2 components c0, c1
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(std::floating_point<T>)
struct MVec2_t {

    // ctors

    // (all grades = 0)
    MVec2_t() = default;

    // assign all components
    constexpr MVec2_t(T s, T ps) : c0(s), c1(ps) {}

    // cp & mv ctor
    MVec2_t(MVec2_t const&) = default;
    MVec2_t(MVec2_t&&) noexcept = default;

    // cp & mv assignment
    MVec2_t& operator=(MVec2_t const&) = default;
    MVec2_t& operator=(MVec2_t&&) noexcept = default;

    // dtor
    ~MVec2_t() = default;

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    constexpr MVec2_t(MVec2_t<U, Tag> const& v) : MVec2_t(v.c0, v.c1)
    {
    }

    friend void swap(MVec2_t& lhs, MVec2_t& rhs) noexcept
    {
        using std::swap;
        swap(static_cast<T&>(lhs), static_cast<T&>(rhs));
    }

    ////////////////////////////////////////////////////////////////////////////
    // component definition
    ////////////////////////////////////////////////////////////////////////////

    T c0{}; // scalar component
    T c1{}; // bivector component (2d pseudoscalar)

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec2_t<U, Tag> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps) return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2_t& operator+=(MVec2_t<U, Tag> const& v)
    {
        c0 += v.c0;
        c1 += v.c1;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2_t& operator-=(MVec2_t<U, Tag> const& v)
    {
        c0 -= v.c0;
        c1 -= v.c1;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2_t& operator*=(U s)
    {
        c0 *= s;
        c1 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2_t& operator/=(U s)
    {
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (s < std::numeric_limits<U>::epsilon()) {
            throw std::runtime_error("scalar too small for division " +
                                     std::to_string(s) + "\n");
        }
#endif
        c0 /= s;
        c1 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus for multivectors from the even subalgebra
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr MVec2_t<T, Tag> operator-(MVec2_t<T, Tag> const& v)
{
    return MVec2_t<T, Tag>(-v.c0, -v.c1);
}

// add multivectors from the even subalgebra
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2_t<std::common_type_t<T, U>, Tag>
operator+(MVec2_t<T, Tag> const& v1, MVec2_t<U, Tag> const& v2)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v1.c0 + v2.c0, v1.c1 + v2.c1);
}

// substract multivectors multivectors from the even subalgebra
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2_t<std::common_type_t<T, U>, Tag>
operator-(MVec2_t<T, Tag> const& v1, MVec2_t<U, Tag> const& v2)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v1.c0 - v2.c0, v1.c1 - v2.c1);
}


// multiply a multivector multivectors from the even subalgebra with a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2_t<std::common_type_t<T, U>, Tag>
operator*(MVec2_t<T, Tag> const& v, U s)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2_t<std::common_type_t<T, U>, Tag>
operator*(T s, MVec2_t<U, Tag> const& v)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s);
}

// devide an even multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2_t<std::common_type_t<T, U>, Tag>
operator/(MVec2_t<T, Tag> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(s) < std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                      std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero " +
                                 std::to_string(s) + "\n");
    }
#endif
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec2_t<ctype, Tag>(v.c0 * inv, v.c1 * inv);
}

// magnitude of a k-vector
// return squared magnitude (e.g. of complex number)
// |Z|^2 = Z rev(Z) = c0^2 + c1^2
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm_sq(MVec2_t<T, Tag> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1;
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(MVec2_t<T, Tag> const& v)
{
    return std::sqrt(nrm_sq(v));
}

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec2_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << ")";
    return os;
}

} // namespace hd::ga