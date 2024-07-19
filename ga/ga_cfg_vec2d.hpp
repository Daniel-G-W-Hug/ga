#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::clamp
#include <cmath>     // std::abs, std::sin, std::cos
#include <concepts>  // std::floating_point<T>
#include <iostream>
#include <limits>
#include <numbers> // math constants like pi
#include <stdexcept>
#include <string>

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "ga_cfg_value_t.hpp"


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
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return Vec2d<std::common_type_t<T, U>>(v.x * inv, v.y * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return dot-product of two vectors
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    // this implementation is only valid in an orthonormal basis
    return v1.x * v2.x + v1.y * v2.y;
}

// return squared magnitude of vector
template <typename T> inline T sq_nrm(Vec2d<T> const& v) { return dot(v, v); }

// return magnitude of vector
template <typename T> inline T nrm(Vec2d<T> const& v) { return std::sqrt(dot(v, v)); }

// return a vector unitized to nrm(v) == 1.0
template <typename T> inline Vec2d<T> unitized(Vec2d<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec2d<T> inv(Vec2d<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // inverse of squared norm for a vector
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// wedge product (returns a bivector, which is the pseudoscalar in 2d)
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    return PScalar2d<std::common_type_t<T, U>>(v1.x * v2.y - v1.y * v2.x);
}

// return the angle between of two vectors
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    using std::numbers::pi;

    ctype nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }

    // std::clamp must be used to take care of numerical inaccuracies
    auto cos_angle = std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
    auto sin_angle = std::clamp(ctype(wdg(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
    // wdg() in 2d contains magnitude and orientation, but works this easy only in 2d,
    // because it is already a scalar value
    // (for 3d to be as effective, the 3d vectors would need to be transformed
    //  to a plane, the angle measured w.r.t. to the pseudoscalar of the plane)

    // fmt::println("   c = {: .4f}, s = {: .4f}, wdg = {: .4f}, nrm_wdg = {: .4f}",
    //              cos_angle, sin_angle, wdg(v1, v2), nrm(wdg(v1, v2)));

    if (cos_angle >= 0.0) {
        // quadrant I or IV
        return std::asin(sin_angle);
    }
    else if (cos_angle < 0.0 && sin_angle >= 0.0) {
        // quadrant II
        return pi - std::asin(sin_angle);
    }
    else {
        // cos_angle < 0.0 && sin_angle < 0.0)
        // quadrant III
        return -pi - std::asin(sin_angle);
    }
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


////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T> struct fmt::formatter<hd::ga::Vec2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{})", nested(v.x), nested(v.y));
    }
};

// Usage:
//
// std::vector<Vec2d<double>> vp1{{1.0, 1.0}, {1.5, 2.0}};
// Vec2d p{1.0, 2.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));
