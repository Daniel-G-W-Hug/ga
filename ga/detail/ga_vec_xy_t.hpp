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
// Vec_xy_t<T, Tag> definition
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag>
    requires(std::floating_point<T>)
struct Vec_xy_t {

    // assumes a right-handed orthonormal vector basis {e1, e2} in EGA
    // using components {x, y}, such that for each vector v:
    // v = x * e1 + y * e2

    // ctors
    Vec_xy_t() = default;

    Vec_xy_t(T x_in, T y_in) : x(x_in), y(y_in) {}

    template <typename U>
        requires(std::floating_point<U>)
    Vec_xy_t(Vec_xy_t<U, Tag> const& v) : x(v.x), y(v.y)
    {
    }

    T x{};
    T y{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(Vec_xy_t<U, Tag> const& rhs) const
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
};

////////////////////////////////////////////////////////////////////////////////
// Vec_xy_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Vec_xy_t<T, Tag> operator-(Vec_xy_t<T, Tag> const& v)
{
    return Vec_xy_t<T, Tag>(-v.x, -v.y);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec_xy_t<std::common_type_t<T, U>, Tag>
operator+(Vec_xy_t<T, Tag> const& v1, Vec_xy_t<U, Tag> const& v2)
{
    return Vec_xy_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec_xy_t<std::common_type_t<T, U>, Tag>
operator-(Vec_xy_t<T, Tag> const& v1, Vec_xy_t<U, Tag> const& v2)
{
    return Vec_xy_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec_xy_t<std::common_type_t<T, U>, Tag>
operator*(Vec_xy_t<T, Tag> const& v, U s)
{
    return Vec_xy_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec_xy_t<std::common_type_t<T, U>, Tag>
operator*(T s, Vec_xy_t<U, Tag> const& v)
{
    return Vec_xy_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec_xy_t<std::common_type_t<T, U>, Tag>
operator/(Vec_xy_t<T, Tag> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                                   std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec_xy_t<std::common_type_t<T, U>, Tag>(v.x * inv, v.y * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec_xy_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Vec_xy_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << ")";
    return os;
}

} // namespace hd::ga