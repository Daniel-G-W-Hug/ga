#pragma once

// author: Daniel Hug, 2024 & 2025

#include <algorithm> // std::max
#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_type_tags.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T, Tag> definition
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(std::floating_point<T>)
struct Vec4_t {

    // assumes a right-handed orthonormal vector basis {e1, e2, e3, e4}
    // using components {x, y, z, w}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3 + w * e4
    //
    // and for each trivector tv:
    // tv = wyz * e4^e2^e3 + wzx * e4^e3^e1 + wxy * e4^e1^e2 + zyx * e3^e2^e1
    //    =  x  * e4^e2^e3 +  y  * e4^e3^e1 +  z  * e4^e1^e2 +  w  * e^4
    // (same names like in vector (x, y, z, w),
    //  but semantic as one live above (wyz, wzx, wxy, zyx)

    // this is a mapping of the components
    // such that vector components x, y, z and w correspond to the
    // trivector elements represented by wyz, wzx, wxy and zyx
    // i.e. they can by converted to each other by a duality transformation
    //
    // T.x <=> trivector wyz <=> Vec4_t<T>::x; // maps to basis trivector e4^e2^e3
    // T.y <=> trivector wzx <=> Vec4_t<T>::y; // maps to basis trivector e4^e3^e1
    // T.z <=> trivector wxy <=> Vec4_t<T>::z; // maps to basis trivector e4^e1^e2
    // T.w <=> trivector zyx <=> Vec4_t<T>::w; // maps to basis trivector e3^e2^e1

    // ctors
    Vec4_t() = default;

    Vec4_t(T x_in, T y_in, T z_in, T w_in) : x(x_in), y(y_in), z(z_in), w(w_in) {}

    // cp & mv ctor
    Vec4_t(Vec4_t const&) = default;
    Vec4_t(Vec4_t&&) noexcept = default;

    // cp & mv assignment
    Vec4_t& operator=(Vec4_t const&) = default;
    Vec4_t& operator=(Vec4_t&&) noexcept = default;

    // dtor
    ~Vec4_t() = default;

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    Vec4_t(Vec4_t<U, Tag> const& v) : x(v.x), y(v.y), z(v.z), w(v.w)
    {
    }

    friend void swap(Vec4_t& lhs, Vec4_t& rhs) noexcept
    {
        using std::swap;
        swap(static_cast<T&>(lhs), static_cast<T&>(rhs));
    }

    T x{};
    T y{};
    T z{};
    T w{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(Vec4_t<U, Tag> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_x = std::abs(rhs.x - x);
        auto abs_delta_y = std::abs(rhs.y - y);
        auto abs_delta_z = std::abs(rhs.z - z);
        auto abs_delta_w = std::abs(rhs.w - w);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_x < delta_eps && abs_delta_y < delta_eps &&
            abs_delta_z < delta_eps && abs_delta_w < delta_eps)
            return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec4_t& operator+=(Vec4_t<U, Tag> const& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        w += v.w;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec4_t& operator-=(Vec4_t<U, Tag> const& v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        w -= v.w;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec4_t& operator*=(U s)
    {
        x *= s;
        y *= s;
        z *= s;
        w *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec4_t& operator/=(U s)
    {
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (s < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("vector norm too small for normalization " +
                                     std::to_string(s) + "\n");
        }
#endif
        x /= s;
        y /= s;
        z /= s;
        w /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Vec4_t<T, Tag> operator-(Vec4_t<T, Tag> const& v)
{
    return Vec4_t<T, Tag>(-v.x, -v.y, -v.z, -v.w);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator+(Vec4_t<T, Tag> const& v1,
                                                                 Vec4_t<U, Tag> const& v2)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z,
                                                 v1.w + v2.w);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator-(Vec4_t<T, Tag> const& v1,
                                                                 Vec4_t<U, Tag> const& v2)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z,
                                                 v1.w - v2.w);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator*(Vec4_t<T, Tag> const& v,
                                                                 U s)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s, v.w * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator*(T s,
                                                                 Vec4_t<U, Tag> const& v)
{
    return Vec4_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s, v.w * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec4_t<std::common_type_t<T, U>, Tag> operator/(Vec4_t<T, Tag> const& v,
                                                                 U s)
{
    using ctype = std::common_type_t<T, U>;
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(s) < std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                      std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero " +
                                 std::to_string(s) + "\n");
    }
#endif
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec4_t<ctype, Tag>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

// magnitude of the k-vector (in representational space)
//
// magnitude is always defined in the representational space, i.e. without
// covering the target metric of the representation
//
// -> in ega magnitude (magn) and norm (nrm) are identical
//
// -> in pga magnitude is defined as norm for the representational space
//    assuming the corresponding identity matrix as metric for that space:
//
//    G<3,0,0> as representational space for modelling 2d Euclidean space using G<2,0,1>
//             representational space:        e1^2=+1, e2^2=+1, e3^2=+1
//             modelled 2d Euclidean space:   e1^2=+1, e2^2=+1, e3^2= 0
//
//    G<4,0,0> as representational space for modelling 3d Euclidean space using G<3,0,1>
//             representational space:        e1^2=+1, e2^2=+1, e3^2=+1, e4^2=+1
//             modelled 3d Euclidean space:   e1^2=+1, e2^2=+1, e3^2=+1, e4^2= 0
//
//    norm is always defined for the modelled space incl. the target metric,
//    by using the exomorphisim matrix G as defined by Lengyel in the book
//    "Projective geometric algebra illuminated"
//
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm_sq(Vec4_t<T, Tag> const& v)
{
    // implements the scalar product as defined by the geometric product *
    // for the vector and trivector case in 4d and for the modeling space in 3dp:
    // |v|^2 = gr0( v*rev(v))
    return v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(Vec4_t<T, Tag> const& v)
{
    return std::sqrt(nrm_sq(v));
}

// return a vector v normalized to nrm(v) == 1.0
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Vec4_t<T, Tag> normalize(Vec4_t<T, Tag> const& v)
{
    T m = nrm(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (m < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization " +
                                 std::to_string(m) + "\n");
    }
#endif
    T inv = T(1.0) / m; // for multiplication with inverse of norm
    return Vec4_t<T, Tag>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Vec4_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << "," << v.w << ")";
    return os;
}

} // namespace hd::ga