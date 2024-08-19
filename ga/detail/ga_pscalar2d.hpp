#pragma once

// author: Daniel Hug, 2024

#include <concepts> // std::floating_point
#include <iostream>

#include "ga_value_t.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// PScalar2d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> operator-(PScalar2d<T> v)
{
    return PScalar2d<T>(-value_t(v));
}

// adding pseudoscalars
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator+(PScalar2d<T> v1,
                                                               PScalar2d<U> v2)
{
    return PScalar2d<std::common_type_t<T, U>>(value_t(v1) + value_t(v2));
}

// substracting pseudoscalars
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator-(PScalar2d<T> v1,
                                                               PScalar2d<U> v2)
{
    return PScalar2d<std::common_type_t<T, U>>(value_t(v1) - value_t(v2));
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> v, U s)
{
    return PScalar2d<std::common_type_t<T, U>>(value_t(v) * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator*(T s, PScalar2d<U> v)
{
    return PScalar2d<std::common_type_t<T, U>>(s * value_t(v));
}

// devide a pseudoscalar by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator/(PScalar2d<T> v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return PScalar2d<ctype>(value_t(v) * inv);
}

////////////////////////////////////////////////////////////////////////////////
// PScalar2d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, PScalar2d<T> v)
{
    os << "(" << T(v) << ")";
    return os;
}

} // namespace hd::ga