#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

// strong type inspired by Jonathan Boccara's blog:
// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
//
// and by the blog of foonathan
// https://www.foonathan.net/2016/10/strong-typedefs/


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

    // these operators bring implicit conversion to the underlying type T
    // (hint: can be suppressed with marking the operator explicit)
    operator T&() noexcept { return value; }
    operator T const&() const noexcept { return value; }

    friend void swap(Scalar_t& lhs, Scalar_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.value, rhs.value);
    }

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
        value *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "scalar division");
        value /= s;
        return (*this);
    }

  private:

    T value{};
};

////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Scalar_t<T, Tag> operator-(Scalar_t<T, Tag> v)
{
    return Scalar_t<T, Tag>(-T(v));
}

// adding pseudoscalars
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator+(Scalar_t<T, Tag> v1,
                                                                   Scalar_t<U, Tag> v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(v1) + ctype(v2));
}

// substracting pseudoscalars
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator-(Scalar_t<T, Tag> v1,
                                                                   Scalar_t<U, Tag> v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(v1) - ctype(v2));
}

// multiply a pseudoscalar with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator*(Scalar_t<T, Tag> v,
                                                                   U s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(v) * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator*(T s,
                                                                   Scalar_t<U, Tag> v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(s * ctype(v));
}

// devide a pseudoscalar by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Scalar_t<std::common_type_t<T, U>, Tag> operator/(Scalar_t<T, Tag> v, U s)
{
    detail::check_division_by_zero<T, U>(s, "scalar division");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Scalar_t<ctype, Tag>(ctype(v) * inv);
}

// magnitude squared
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm_sq(Scalar_t<T, Tag> s)
{
    return T(s) * T(s);
}

// magnitude
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(Scalar_t<T, Tag> s)
{
    return std::sqrt(nrm_sq(s));
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
std::ostream& operator<<(std::ostream& os, Scalar_t<T, Tag> v)
{
    os << "(" << T(v) << ")";
    return os;
}

// return the scalar value as value_t (for use in scripting)
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T to_val(Scalar_t<T, Tag> s)
{
    return T(s);
}

} // namespace hd::ga