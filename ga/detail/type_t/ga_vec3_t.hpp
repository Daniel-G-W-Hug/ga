#pragma once

// author: Daniel Hug, 2024

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
// Vec3_t<T, Tag> definition
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(std::floating_point<T>)
struct Vec3_t {

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
    // T.x <=> bivector yz <=> Vec3_t<T>::x; // maps to basis bivector e2^e3
    // T.y <=> bivector zx <=> Vec3_t<T>::y; // maps to basis bivector e3^e1
    // T.z <=> bivector xy <=> Vec3_t<T>::z; // maps to basis bivector e1^e2

    // ctors
    Vec3_t() = default;

    Vec3_t(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {}

    // cp & mv ctor
    Vec3_t(Vec3_t const&) = default;
    Vec3_t(Vec3_t&&) noexcept = default;

    // cp & mv assignment
    Vec3_t& operator=(Vec3_t const&) = default;
    Vec3_t& operator=(Vec3_t&&) noexcept = default;

    // dtor
    ~Vec3_t() = default;

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    Vec3_t(Vec3_t<U, Tag> const& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    friend void swap(Vec3_t& lhs, Vec3_t& rhs) noexcept
    {
        using std::swap;
        swap(static_cast<T&>(lhs), static_cast<T&>(rhs));
    }

    T x{};
    T y{};
    T z{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(Vec3_t<U, Tag> const& rhs) const
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
        requires(std::floating_point<U>)
    Vec3_t& operator+=(Vec3_t<U, Tag> const& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec3_t& operator-=(Vec3_t<U, Tag> const& v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec3_t& operator*=(U s)
    {
        x *= s;
        y *= s;
        z *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    Vec3_t& operator/=(U s)
    {
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (s < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("vector norm too small for normalization" +
                                     std::to_string(s) + "\n");
        }
#endif
        x /= s;
        y /= s;
        z /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Vec3_t<T, Tag> operator-(Vec3_t<T, Tag> const& v)
{
    return Vec3_t<T, Tag>(-v.x, -v.y, -v.z);
}

// adding vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator+(Vec3_t<T, Tag> const& v1,
                                                                 Vec3_t<U, Tag> const& v2)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

// substracting vectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator-(Vec3_t<T, Tag> const& v1,
                                                                 Vec3_t<U, Tag> const& v2)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator*(Vec3_t<T, Tag> const& v,
                                                                 U s)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator*(T s,
                                                                 Vec3_t<U, Tag> const& v)
{
    return Vec3_t<std::common_type_t<T, U>, Tag>(v.x * s, v.y * s, v.z * s);
}

// devide a vector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3_t<std::common_type_t<T, U>, Tag> operator/(Vec3_t<T, Tag> const& v,
                                                                 U s)
{
    using ctype = std::common_type_t<T, U>;
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(s) < std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                      std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
#endif
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Vec3_t<ctype, Tag>(v.x * inv, v.y * inv, v.z * inv);
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
inline constexpr T nrm_sq(Vec3_t<T, Tag> const& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(Vec3_t<T, Tag> const& v)
{
    return std::sqrt(nrm_sq(v));
}

// return a vector v normalized to nrm(v) == 1.0
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Vec3_t<T, Tag> normalize(Vec3_t<T, Tag> const& v)
{
    T m = nrm(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (m < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(m) + "\n");
    }
#endif
    T inv = T(1.0) / m; // for multiplication with inverse of norm
    return Vec3_t<T, Tag>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Vec3_t<T, Tag> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

} // namespace hd::ga