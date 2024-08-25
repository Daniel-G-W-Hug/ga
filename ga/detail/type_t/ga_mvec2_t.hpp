#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::max
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_type_tags.hpp"

namespace hd::ga {

template <typename T, typename Tag = default_tag>
    requires(std::floating_point<T>)
struct MVec2_t {

    // ctors

    // (all grades = 0)
    MVec2_t() = default;

    // assign all components
    MVec2_t(T s, T ps) : c0(s), c1(ps) {}

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
    if (std::abs(s) < ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                                   std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec2_t<ctype, Tag>(v.c0 * inv, v.c1 * inv);
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