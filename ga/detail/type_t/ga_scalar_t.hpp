#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// strong type inspired by Jonathan Boccara's blog:
// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
//
// and by the blog of foonathan
// https://www.foonathan.net/2016/10/strong-typedefs/
/////////////////////////////////////////////////////////////////////////////////////////


#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // std::floating_point
#include <iostream>  // std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

// provide common type for Scalar2d<T>, Scalar3d<T>, ..., PScalar2d<T>, PScalar3d<T>, ...

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"

namespace hd::ga {

template <typename T, typename Tag = default_tag>
    requires(std::floating_point<T>)
class Scalar_t {

  public:

    Scalar_t() = default; // default initialization
    explicit Scalar_t(T const& val) : value(val) {}
    explicit Scalar_t(T&& val) noexcept(std::is_nothrow_move_constructible<T>::value) :
        value(std::move(val))
    {
    }

    // These operators bring implicit conversion to the underlying type T,
    // if not made explicit!
    // HINT: This is unfortunate, because implicit conversion to floating point type
    // allows for comparison of scalar and pscalar types, or for comparison of 2d and 3d
    // scalars, despite the comparison operator requiring same types
    // (hint: this must be suppressed with marking the operator explicit for type safety!)
    explicit operator T&() noexcept { return value; }
    explicit operator T const&() const noexcept { return value; }

    friend void swap(Scalar_t& lhs, Scalar_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.value, rhs.value);
    }

    // equality (only compare, if same tag, i.e. don't mix up scalars and pseudoscalars,
    //           or 2d and 3d scalars, for example)
    template <typename U, typename U_Tag>
        requires(std::floating_point<U> && std::is_same_v<Tag, U_Tag>)
    bool operator==(Scalar_t<U, U_Tag> rhs) const
    {
        // equality implies same magnitude
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta = std::abs(T(value) - U(rhs));
        auto constexpr delta_eps = detail::safe_epsilon<T, U>();
        return (abs_delta < delta_eps);
    }

    // inequality
    template <typename U, typename U_Tag>
        requires(std::floating_point<U> && std::is_same_v<Tag, U_Tag>)
    bool operator!=(Scalar_t<U, U_Tag> rhs) const
    {
        return !(*this == rhs);
    }

    // less than
    template <typename U, typename U_Tag>
        requires(std::floating_point<U> && std::is_same_v<Tag, U_Tag>)
    bool operator<(Scalar_t<U, U_Tag> rhs) const
    {
        return (value < U(rhs));
    }

    // less than or equal
    template <typename U, typename U_Tag>
        requires(std::floating_point<U> && std::is_same_v<Tag, U_Tag>)
    bool operator<=(Scalar_t<U, U_Tag> rhs) const
    {
        return (*this < rhs) || (*this == rhs);
    }

    // greater than
    template <typename U, typename U_Tag>
        requires(std::floating_point<U> && std::is_same_v<Tag, U_Tag>)
    bool operator>(Scalar_t<U, U_Tag> rhs) const
    {
        return rhs < *this;
    }

    // greater than or equal
    template <typename U, typename U_Tag>
        requires(std::floating_point<U> && std::is_same_v<Tag, U_Tag>)
    bool operator>=(Scalar_t<U, U_Tag> rhs) const
    {
        return (*this > rhs) || (*this == rhs);
    }

    // Comparison operators with raw numeric types (floating point and integral)
    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    bool operator==(U rhs) const
    {
        // equality implies same magnitude
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta = std::abs(T(value) - U(rhs));
        auto constexpr delta_eps = detail::safe_epsilon<T, U>();
        return (abs_delta < delta_eps);
    }

    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    bool operator!=(U rhs) const
    {
        return !(*this == rhs);
    }

    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    bool operator<(U rhs) const
    {
        return value < rhs;
    }

    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    bool operator<=(U rhs) const
    {
        return (*this < rhs) || (*this == rhs);
    }

    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    bool operator>(U rhs) const
    {
        return value > rhs;
    }

    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    bool operator>=(U rhs) const
    {
        return (*this > rhs) || (*this == rhs);
    }

    // Arithmetic operators with raw numeric types
    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    Scalar_t<std::common_type_t<T, U>, Tag> operator+(U rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        return Scalar_t<ctype, Tag>(ctype(value) + ctype(rhs));
    }

    template <typename U>
        requires(std::floating_point<U> || std::integral<U>)
    Scalar_t<std::common_type_t<T, U>, Tag> operator-(U rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        return Scalar_t<ctype, Tag>(ctype(value) - ctype(rhs));
    }


    ////////////////////////////////////////////////////////////////////////////////
    // other operators for raw floating point types with Scalar_t
    ////////////////////////////////////////////////////////////////////////////////

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator+=(Scalar_t<U, Tag> s) noexcept
    {
        value += s.value;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator-=(Scalar_t<U, Tag> s) noexcept
    {
        value -= s.value;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator*=(U s) noexcept
    {
        value *= s.value;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s.value, "scalar division");
        value /= s.value;
        return (*this);
    }

  private:

    T value{};
};

////////////////////////////////////////////////////////////////////////////////
// Reverse operators for raw floating point types with Scalar_t
////////////////////////////////////////////////////////////////////////////////

// Comparison operators: raw type with Scalar_t
template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
bool operator==(T lhs, Scalar_t<U, Tag> rhs)
{
    return rhs == lhs; // Use the member operator
}

template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
bool operator!=(T lhs, Scalar_t<U, Tag> rhs)
{
    return rhs != lhs; // Use the member operator
}

template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
bool operator<(T lhs, Scalar_t<U, Tag> rhs)
{
    return lhs < U(rhs);
}

template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
bool operator<=(T lhs, Scalar_t<U, Tag> rhs)
{
    return (lhs < rhs) || (lhs == rhs);
}

template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
bool operator>(T lhs, Scalar_t<U, Tag> rhs)
{
    return lhs > U(rhs);
}

template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
bool operator>=(T lhs, Scalar_t<U, Tag> rhs)
{
    return (lhs > rhs) || (lhs == rhs);
}

// Arithmetic operators: raw type with Scalar_t
template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator+(T lhs, Scalar_t<U, Tag> rhs)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(lhs) + ctype(rhs));
}

template <typename T, typename U, typename Tag>
    requires((std::floating_point<T> || std::integral<T>) && std::floating_point<U>)
constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator-(T lhs, Scalar_t<U, Tag> rhs)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(lhs) - ctype(rhs));
}


////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr Scalar_t<T, Tag> operator-(Scalar_t<T, Tag> s)
{
    return Scalar_t<T, Tag>(-T(s));
}

// adding pseudoscalars
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator+(Scalar_t<T, Tag> s1,
                                                            Scalar_t<U, Tag> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(s1) + ctype(s2));
}

// substracting pseudoscalars
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator-(Scalar_t<T, Tag> s1,
                                                            Scalar_t<U, Tag> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(s1) - ctype(s2));
}

// multiply a pseudoscalar with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator*(Scalar_t<T, Tag> s1, U s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(s1) * s2);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator*(T s1, Scalar_t<U, Tag> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(s1 * ctype(s2));
}

// devide a pseudoscalar by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Scalar_t<std::common_type_t<T, U>, Tag> operator/(Scalar_t<T, Tag> s1, U s2)
{
    detail::check_division_by_zero<T, U>(s2, "scalar division");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s2; // for multiplicaton with inverse value
    return Scalar_t<ctype, Tag>(ctype(s1) * inv);
}

// magnitude squared
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T nrm_sq(Scalar_t<T, Tag> s)
{
    return T(s) * T(s);
}

// magnitude
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T nrm(Scalar_t<T, Tag> s)
{
    return sqrt(nrm_sq(s));
}

// return a scalar normalized to nrm(s) == 1.0
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline Scalar_t<T, Tag> normalize(Scalar_t<T, Tag> s)
{
    T m = nrm(s);
    detail::check_normalization<T>(m, "scalar or pseudoscalar");
    T inv = T(1.0) / m; // for multiplication with inverse of norm
    return Scalar_t<T, Tag>(T(s) * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Scalar_t<T, Tag> s)
{
    os << "(" << T(s) << ")";
    return os;
}

// return the scalar value as value_t (for use in scripting)
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T to_val(Scalar_t<T, Tag> s)
{
    return T(s);
}

////////////////////////////////////////////////////////////////////////////////
// abs and sqrt for Scalar_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
// Unified mathematical functions for consistent interface throughout the GA library
//
// IMPORTANT: The GA library uses unqualified calls to mathematical functions (abs, sqrt,
// etc.) instead of std::abs, std::sqrt, etc. This allows these unified functions to be
// found via ADL (Argument Dependent Lookup) when called with Scalar_t types, while still
// working correctly with raw numeric types.
//
// This design provides:
// - Consistent interface for both raw types (double, float) and Scalar_t types
// - Seamless integration throughout the library without explicit conversions
// - Natural mathematical notation in user code
//
// Usage: Use unqualified function calls like abs(x) and sqrt(x) instead of std::abs(x),
// sqrt(x)

// Unified abs function that works for both raw numeric types and Scalar_t types
template <typename T>
    requires(std::floating_point<T> || std::integral<T>)
constexpr T abs(T value)
{
    return std::abs(value);
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T abs(Scalar_t<T, Tag> s)
{
    return std::abs(T(s));
}

// Unified sqrt function that works for both raw numeric types and Scalar_t types
template <typename T>
    requires(std::floating_point<T> || std::integral<T>)
constexpr auto sqrt(T value) -> decltype(std::sqrt(value))
{
    return std::sqrt(value);
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T sqrt(Scalar_t<T, Tag> s)
{
    return std::sqrt(T(s));
}

} // namespace hd::ga