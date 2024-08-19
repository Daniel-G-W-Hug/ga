#pragma once

// author: Daniel Hug, 2024

#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_value_t.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> definition (used for implementation of algebra<2,0,0>)
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct Vec2d {

    // assumes a right-handed orthonormal vector basis {e1, e2}
    // using components {x, y}, such that for each vector v:
    // v = x * e1 + y * e2

    // ctors
    Vec2d() = default;

    Vec2d(T x_in, T y_in) : x(x_in), y(y_in) {}

    template <typename U>
        requires(std::floating_point<U>)
    Vec2d(Vec2d<U> const& v) : x(v.x), y(v.y)
    {
    }

    T x{};
    T y{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(Vec2d<U> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_x = std::abs(rhs.x - x);
        auto abs_delta_y = std::abs(rhs.y - y);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_x < delta_eps && abs_delta_y < delta_eps) return true;
        return false;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, Vec2d<U> const& v);
};

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> operator-(Vec2d<T> const& v)
{
    return Vec2d<T>(-v.x, -v.y);
}

// adding vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator+(Vec2d<T> const& v1,
                                                           Vec2d<U> const& v2)
{
    return Vec2d<std::common_type_t<T, U>>(v1.x + v2.x, v1.y + v2.y);
}

// substracting vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator-(Vec2d<T> const& v1,
                                                           Vec2d<U> const& v2)
{
    return Vec2d<std::common_type_t<T, U>>(v1.x - v2.x, v1.y - v2.y);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v, U s)
{
    return Vec2d<std::common_type_t<T, U>>(v.x * s, v.y * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(T s, Vec2d<U> const& v)
{
    return Vec2d<std::common_type_t<T, U>>(v.x * s, v.y * s);
}

// devide a vector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator/(Vec2d<T> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec2d<std::common_type_t<T, U>>(v.x * inv, v.y * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Vec2d<T> const& v)
{
    os << "(" << v.x << "," << v.y << ")";
    return os;
}

} // namespace hd::ga