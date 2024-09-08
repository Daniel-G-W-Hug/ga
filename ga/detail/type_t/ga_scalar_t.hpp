#pragma once

// author: Daniel Hug, 2024

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
        swap(static_cast<T&>(lhs), static_cast<T&>(rhs));
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator+=(Scalar_t<U, Tag> s)
    {
        value += s.value;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator-=(Scalar_t<U, Tag> s)
    {
        value -= s.value;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator*=(U s)
    {
        value *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Scalar_t& operator/=(U s)
    {
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (s < std::numeric_limits<U>::epsilon()) {
            throw std::runtime_error("scalar too small for divison" + std::to_string(s) +
                                     "\n");
        }
#endif
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
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator/(Scalar_t<T, Tag> v,
                                                                   U s)
{
    using ctype = std::common_type_t<T, U>;
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(s) < std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                      std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
#endif
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Scalar_t<ctype, Tag>(ctype(v) * inv);
}

// magnitude squared
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T magn_sq(Scalar_t<T, Tag> s)
{
    return T(s) * T(s);
}

// magnitude
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T magn(Scalar_t<T, Tag> s)
{
    return std::sqrt(magn_sq(s));
}

// return a scalar normalized to magn(s) == 1.0
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Scalar_t<T, Tag> normalize(Scalar_t<T, Tag> s)
{
    T m = magn(s);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (m < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("magnitude too small for normalization" +
                                 std::to_string(m) + "\n");
    }
#endif
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

////////////////////////////////////////////////////////////////////////////////
// Scalar types - basic operations (identical for all algebras)
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm_sq(Scalar_t<T, Tag> s)
{
    return T(s) * T(s);
}

// return absolute magnitude
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(Scalar_t<T, Tag> s)
{
    return std::abs(T(s));
}

// return the scalar value as value_t (for use in scripting)
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T to_val(Scalar_t<T, Tag> s)
{
    return T(s);
}

} // namespace hd::ga