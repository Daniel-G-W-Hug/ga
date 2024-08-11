#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::clamp
#include <cmath>     // std::abs, std::sin, std::cos
#include <concepts>  // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "ga_cfg_value_t.hpp"


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
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return Vec3d<std::common_type_t<T, U>>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return the dot product of two vectors (= a scalar)
// coordinate free definition: dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    // this implementation is only valid in an orthonormal basis
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// return squared magnitude of vector
template <typename T> inline T sq_nrm(Vec3d<T> const& v) { return dot(v, v); }

// return magnitude of vector
template <typename T> inline T nrm(Vec3d<T> const& v) { return std::sqrt(dot(v, v)); }

// return a vector unitized to nrm(v) == 1.0
template <typename T> inline Vec3d<T> unitized(Vec3d<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec3d<T> inv(Vec3d<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // inverse of squared norm for a vector
    return Vec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    ctype nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


// unsuccessful try to extend angle range to -pi ... pi,
// because orientation is not defined uniquely:
//
// // return the angle between of two vectors
// // range of angle: -pi <= angle <= pi
// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U> const& v2)
// {
//     using ctype = std::common_type_t<T, U>;
//     using std::numbers::pi;

//     ctype nrm_prod = nrm(v1) * nrm(v2);
//     if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
//         throw std::runtime_error(
//             "vector norm product too small for calculation of angle" +
//             std::to_string(nrm_prod) + "\n");
//     }

// auto cos_angle = std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
// auto sin_angle = std::clamp(ctype(nrm(wdg(v1, v2))) / nrm_prod, ctype(-1.0),
// ctype(1.0));
//     // wdg() does contain magnitude, but no unique value of orientation
//     // so we chose one arbitrarily => but would deliver only pos. angles!

//     fmt::println("   c = {: .4f}, s = {: .4f}, wdg = {: .4f}, nrm_wdg = {: .4f}",
//                  cos_angle, sin_angle, wdg(v1, v2), nrm(wdg(v1, v2)));

//     if (cos_angle >= 0.0) {
//         // quadrant I or IV
//         return std::asin(sin_angle);
//     }
//     else if (cos_angle < 0.0 && sin_angle >= 0.0) {
//         // quadrant II
//         return pi - std::asin(sin_angle);
//     }
//     else {
//         // cos_angle < 0.0 && sin_angle < 0.0)
//         // quadrant III
//         return -pi - std::asin(sin_angle);
//     }
// }

// cross-product between two vectors (only defined in 3d)
//  => returns a 3d vector
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> cross(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}
//     double cross product identity:
//     a x (b x c) = -dot(a, wdg(b,c)) = -a.(b^c)

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


////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T> struct fmt::formatter<hd::ga::Vec3d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec3d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};

// Usage:
//
// std::vector<Vec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// Vec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));
