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
// Vec3d<T> definition (used for implementation of algebra<3,0,0>)
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct Vec3d {

    // assumes a right-handed orthonormal vector basis {e1, e2, e3}
    // using components {x, y, z}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3

    // ctors
    Vec3d() = default;

    Vec3d(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {}

    template <typename U>
        requires(std::floating_point<U>)
    Vec3d(Vec3d<U> const& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    T x{};
    T y{};
    T z{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(Vec3d<U> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_x = std::abs(rhs.x - x);
        auto abs_delta_y = std::abs(rhs.y - y);
        auto abs_delta_z = std::abs(rhs.z - z);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_x < delta_eps && abs_delta_y < delta_eps && abs_delta_z < delta_eps)
            return true;
        return false;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, Vec3d<U> const& v);
};

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> operator-(Vec3d<T> const& v)
{
    return Vec3d<T>(-v.x, -v.y, -v.z);
}

// adding vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator+(Vec3d<T> const& v1,
                                                           Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

// substracting vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator-(Vec3d<T> const& v1,
                                                           Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v, U s)
{
    return Vec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(T s, Vec3d<U> const& v)
{
    return Vec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

// devide a vector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator/(Vec3d<T> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec3d<std::common_type_t<T, U>>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Vec3d<T> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

} // namespace hd::ga