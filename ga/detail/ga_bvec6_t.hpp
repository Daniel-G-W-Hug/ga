#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::max
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_type_3d.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> definition of a bivector with 6 components vx, vy, vz, lx, ly, lz
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag>
    requires(std::floating_point<T>)
struct BVec6_t {

    // ctors

    // (all grades = 0)
    BVec6_t() = default;

    // assign all components directly
    BVec6_t(T vx, T vy, T vz, T mx, T my, T mz) :
        vx(vx), vy(vy), vz(vz), mx(mx), my(my), mz(mz)
    {
    }

    // assign the vector parts separately
    BVec6_t(Vec3d<T> const& v, Vec3d<T> const& m) :
        vx(v.x), vy(v.y), vz(v.z), mx(m.x), my(m.y), mz(m.z)
    {
    }

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    BVec6_t(BVec6_t<U, Tag> const& v) :
        vx(v.vx), vy(v.vy), vz(v.vz), mx(v.mx), my(v.my), mz(v.mz)
    {
    }


    T vx{}; // as BiVemyd<T> maps to basis bivector e4^e1 - as Line3dp<T> to vx
    T vy{}; // as BiVemyd<T> maps to basis bivector e4^e2 - as Line3dp<T> to vy
    T vz{}; // as BiVemyd<T> maps to basis bivector e4^e3 - as Line3dp<T> to vz
    T mx{}; // as BiVemyd<T> maps to basis bivector e2^e3 - as Line3dp<T> to mx
    T my{}; // as BiVemyd<T> maps to basis bivector e3^e1 - as Line3dp<T> to my
    T mz{}; // as BiVemyd<T> maps to basis bivector e1^e2 - as Line3dp<T> to mz

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(BVec6_t<U, Tag> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_vx = std::abs(rhs.vx - vx);
        auto abs_delta_vy = std::abs(rhs.vy - vy);
        auto abs_delta_vz = std::abs(rhs.vz - vz);
        auto abs_delta_mx = std::abs(rhs.mx - mx);
        auto abs_delta_my = std::abs(rhs.my - my);
        auto abs_delta_mz = std::abs(rhs.mz - mz);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_vx < delta_eps && abs_delta_vy < delta_eps &&
            abs_delta_vz < delta_eps && abs_delta_mx < delta_eps &&
            abs_delta_my < delta_eps && abs_delta_mz < delta_eps)
            return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    BVec6_t& operator+=(BVec6_t<U, Tag> const& v)
    {
        vx += v.vx;
        vy += v.vy;
        vz += v.vz;
        mx += v.mx;
        my += v.my;
        mz += v.mz;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    BVec6_t& operator-=(BVec6_t<U, Tag> const& v)
    {
        vx -= v.vx;
        vy -= v.vy;
        vz -= v.vz;
        mx -= v.mx;
        my -= v.my;
        mz -= v.mz;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    BVec6_t& operator*=(U s)
    {
        vx *= s;
        vy *= s;
        vz *= s;
        mx *= s;
        my *= s;
        mz *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    BVec6_t& operator/=(U s)
    {
        vx /= s;
        vy /= s;
        vz /= s;
        mx /= s;
        my /= s;
        mz /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr BVec6_t<T, Tag> operator-(BVec6_t<T, Tag> const& v)
{
    return BVec6_t<T, Tag>(-v.vx, -v.vy, -v.vz, -v.mx, -v.my, -v.mz);
}

// add multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BVec6_t<std::common_type_t<T, U>, Tag>
operator+(BVec6_t<T, Tag> const& v1, BVec6_t<U, Tag> const& v2)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(v1.vx + v2.vx, v1.vy + v2.vy,
                                                  v1.vz + v2.vz, v1.mx + v2.mx,
                                                  v1.my + v2.my, v1.mz + v2.mz);
}

// substract multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BVec6_t<std::common_type_t<T, U>, Tag>
operator-(BVec6_t<T, Tag> const& v1, BVec6_t<U, Tag> const& v2)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(v1.vx - v2.vx, v1.vy - v2.vy,
                                                  v1.vz - v2.vz, v1.mx - v2.mx,
                                                  v1.my - v2.my, v1.mz - v2.mz);
}

// multiply a multivector with a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BVec6_t<std::common_type_t<T, U>, Tag>
operator*(BVec6_t<T, Tag> const& v, U s)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(v.vx * s, v.vy * s, v.vz * s, v.mx * s,
                                                  v.my * s, v.mz * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BVec6_t<std::common_type_t<T, U>, Tag>
operator*(T s, BVec6_t<U, Tag> const& v)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(v.vx * s, v.vy * s, v.vz * s, v.mx * s,
                                                  v.my * s, v.mz * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BVec6_t<std::common_type_t<T, U>, Tag>
operator/(BVec6_t<T, Tag> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                                   std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return BVec6_t<ctype, Tag>(v.vx * inv, v.vy * inv, v.vz * inv, v.mx * inv, v.my * inv,
                               v.mz * inv);
}

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, BVec6_t<T, Tag> const& v)
{
    os << "(" << v.vx << "," << v.vy << "," << v.vz << "," << v.mx << "," << v.my << ","
       << v.mz << ")";
    return os;
}

} // namespace hd::ga