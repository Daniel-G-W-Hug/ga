#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::clamp
#include <cmath>     // std::abs, std::sin, std::cos
#include <concepts>  // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "ga_value_t.hpp"
#include "ga_vec3d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> definition (used for implementation of algebra<3,0,0>)
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct BiVec3d {

    // assumes a right-handed orthonormal vector basis {e1, e2, e3}
    // using components {x, y, z}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3
    //
    // and for each bivector bv:
    // bv = yz * e2^e3 + zx * e3^e1 + xy * e1^e2
    //    =  x * e2^e3 +  y * e3^e1 +  z * e1^e2
    // (same names like in vector (x, y, z), but semantic as one live above (yz, zx, xy)

    // this is a mapping of the components
    // such that vector components x, y and z correspond to the
    // normals of the corresponding plane elements represented by
    // bivector components yz, zx and xy
    // i.e. they can by converted to each other by a duality transformation
    //
    // T.x <=> bivector yz <=> Vec3d<T>::x; // maps to basis bivector e2^e3
    // T.y <=> bivector zx <=> Vec3d<T>::y; // maps to basis bivector e3^e1
    // T.z <=> bivector xy <=> Vec3d<T>::z; // maps to basis bivector e1^e2

    // => everything otherwise is identical to Vec3d<T> w/o modification.

    BiVec3d() = default;

    BiVec3d(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {}

    template <typename U>
        requires(std::floating_point<U>)
    BiVec3d(BiVec3d<U> const& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    T x{};
    T y{};
    T z{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(BiVec3d<U> const& rhs) const
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
    friend std::ostream& operator<<(std::ostream& os, BiVec3d<U> const& v);
};

////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> operator-(BiVec3d<T> const& v)
{
    return BiVec3d<T>(-v.x, -v.y, -v.z);
}

// adding bivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator+(BiVec3d<T> const& v1,
                                                             BiVec3d<U> const& v2)
{
    return BiVec3d<std::common_type_t<T, U>>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

// substracting bivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator-(BiVec3d<T> const& v1,
                                                             BiVec3d<U> const& v2)
{
    return BiVec3d<std::common_type_t<T, U>>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}


// multiply a bivector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& v, U s)
{
    return BiVec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(T s, BiVec3d<U> const& v)
{
    return BiVec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

// devide a bivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator/(BiVec3d<T> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return BiVec3d<ctype>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, BiVec3d<T> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

} // namespace hd::ga