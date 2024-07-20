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
#include "ga_cfg_vec3d.hpp"


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
    // (same names (x, y, z), but same semantic as one live above (yz, zx, xy)

    // this is a mapping of the components as aliases
    // such that vector components x, y and z correspond to the
    // normals of the corresponding plane elements represented by
    // bivector components yz, zx and xy
    // i.e. they can by converted to each other by a duality transformation
    //
    // T.x <=> bivector yz <=> Vec3d<T>::x; // maps to basis bivector e2^e3
    // T.y <=> bivector zx <=> Vec3d<T>::y; // maps to basis bivector e3^e1
    // T.z <=> bivector xy <=> Vec3d<T>::z; // maps to basis bivector e1^e2

    // duality operations (as defined in Macdonald, "Linear and geometric algebra"):
    //
    // e2^e3 * rev(I_3d) = e2^e3 * e3^e2^e1 = e_23321 = e_1           = e1
    // e3^e1 * rev(I_3d) = e3^e1 * e3^e2^e1 = e_31321 = e_33112 = e_2 = e2
    // e1^e2 * rev(I_3d) = e1^e2 * e3^e2^e1 = e_12321 = e_11223 = e_3 = e3

    // duality operations (as defined in Doran, Lasenby, "GA for physicists"):
    // (same subspace as above, but other sign)
    //
    // I_3d * e2^e3 = e1^e2^e3 * e2^e3 = e_12323 = -e_1            = -e1
    // I_3d * e3^e1 = e1^e2^e3 * e3^e1 = e_12331 = -e_11233 = -e_2 = -e2
    // I_3d * e1^e2 = e1^e2^e3 * e1^e2 = e_12312 = -e_11223 = -e_3 = -e3

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
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return BiVec3d<std::common_type_t<T, U>>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return dot product of two bivectors A and B (= a scalar)
// dot(A,B) = gr0(A * B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> dot(BiVec3d<T> const& A, BiVec3d<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    return -A.x * B.x - A.y * B.y - A.z * B.z;
}

// return squared magnitude of bivector
template <typename T> inline constexpr T sq_nrm(BiVec3d<T> const& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

// return magnitude of bivector
template <typename T> inline constexpr T nrm(BiVec3d<T> const& v)
{
    return std::sqrt(sq_nrm(v));
}

// return a bivector unitized to nrm(v) == 1.0
template <typename T> inline constexpr BiVec3d<T> unitized(BiVec3d<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return BiVec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the bivector
template <typename T> inline constexpr BiVec3d<T> inv(BiVec3d<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = -1.0 / sq_n; // negative inverse of squared norm for a bivector
    return BiVec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return conjugate complex of a bivector
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline constexpr BiVec3d<T> rev(BiVec3d<T> const& v)
{
    // all bivector parts switch sign
    return BiVec3d<T>(-v.x, -v.y, -v.z);
}

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(BiVec3d<T> const& v1, BiVec3d<U> const& v2)
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


////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::BiVec3d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::BiVec3d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};

// Usage:
//
// std::vector<BiVec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// BiVec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));
