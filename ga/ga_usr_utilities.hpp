#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <numbers> // math constants like pi

#include "detail/type_t/ga_scalar_t.hpp"
#include "ga_value_t.hpp"

namespace hd::ga {

using std::numbers::pi; // make pi available for users

constexpr value_t deg2rad(value_t value) { return value / 180.0 * pi; }
constexpr value_t rad2deg(value_t value) { return value * 180.0 / pi; }

constexpr value_t sign(value_t value) { return (value >= 0.0) ? 1.0 : -1.0; }

// sign function overload for Scalar_t types
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T sign(Scalar_t<T, Tag> s)
{
    return sign(T(s));
}

constexpr value_t kronecker(size_t i, size_t j) { return (i == j) ? 1.0 : 0.0; }

// Templates for is_even and is_odd work with any integer type
template <typename T> bool is_even(T arg)
{
    static_assert(std::is_integral_v<T>, "is_even: argument must be an integral type");
    return arg % 2 == 0;
}

template <typename T> bool is_odd(T arg)
{
    static_assert(std::is_integral_v<T>, "is_odd: argument must be an integral type");
    return arg % 2 != 0;
}

} // namespace hd::ga