#pragma once

// author: Daniel Hug, 2024

#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_type_3d.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec6_t<T, Tag> definition of a multivector with 6 components c0, ..., c5
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag>
    requires(std::floating_point<T>)
struct MVec6_t {

    // ctors

    // (all grades = 0)
    MVec6_t() = default;

    // assign all components directly
    MVec6_t(T vx, T vy, T vz, T mx, T my, T mz) :
        c0(vx), c1(vy), c2(vz), c3(mx), c4(my), c5(mz)
    {
    }

    // assign the vector parts separately
    MVec6_t(Vec3d<T> const& v, Vec3d<T> const& m) :
        c0(v.x), c1(v.y), c2(v.z), c3(m.x), c4(m.y), c5(m.z)
    {
    }

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec6_t(MVec6_t<U, Tag> const& v) :
        c0(v.c0), c1(v.c1), c2(v.c2), c3(v.c3), c4(v.c4), c5(v.c5)
    {
    }


    T c0{}; // as BiVec4d<T> maps to basis bivector e4^e1 - as Line3dp<T> to vx
    T c1{}; // as BiVec4d<T> maps to basis bivector e4^e2 - as Line3dp<T> to vy
    T c2{}; // as BiVec4d<T> maps to basis bivector e4^e3 - as Line3dp<T> to vz
    T c3{}; // as BiVec4d<T> maps to basis bivector e2^e3 - as Line3dp<T> to mx
    T c4{}; // as BiVec4d<T> maps to basis bivector e3^e1 - as Line3dp<T> to my
    T c5{}; // as BiVec4d<T> maps to basis bivector e1^e2 - as Line3dp<T> to mz

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec6_t<U, Tag> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto abs_delta_c2 = std::abs(rhs.c2 - c2);
        auto abs_delta_c3 = std::abs(rhs.c3 - c3);
        auto abs_delta_c4 = std::abs(rhs.c4 - c4);
        auto abs_delta_c5 = std::abs(rhs.c5 - c5);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps &&
            abs_delta_c4 < delta_eps && abs_delta_c5 < delta_eps)
            return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec6_t& operator+=(MVec6_t<U, Tag> const& v)
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        c4 += v.c4;
        c5 += v.c5;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec6_t& operator-=(MVec6_t<U, Tag> const& v)
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        c4 -= v.c4;
        c5 -= v.c5;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec6_t& operator*=(U s)
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        c4 *= s;
        c5 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec6_t& operator/=(U s)
    {
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        c4 /= s;
        c5 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec6_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr MVec6_t<T, Tag> operator-(MVec6_t<T, Tag> const& v)
{
    return MVec6_t<T, Tag>(-v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5);
}

// add multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec6_t<std::common_type_t<T, U>, Tag>
operator+(MVec6_t<T, Tag> const& v1, MVec6_t<U, Tag> const& v2)
{
    return MVec6_t<std::common_type_t<T, U>, Tag>(v1.c0 + v2.c0, v1.c1 + v2.c1,
                                                  v1.c2 + v2.c2, v1.c3 + v2.c3,
                                                  v1.c4 + v2.c4, v1.c5 + v2.c5);
}

// substract multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec6_t<std::common_type_t<T, U>, Tag>
operator-(MVec6_t<T, Tag> const& v1, MVec6_t<U, Tag> const& v2)
{
    return MVec6_t<std::common_type_t<T, U>, Tag>(v1.c0 - v2.c0, v1.c1 - v2.c1,
                                                  v1.c2 - v2.c2, v1.c3 - v2.c3,
                                                  v1.c4 - v2.c4, v1.c5 - v2.c5);
}

// multiply a multivector with a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec6_t<std::common_type_t<T, U>, Tag>
operator*(MVec6_t<T, Tag> const& v, U s)
{
    return MVec6_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                                  v.c4 * s, v.c5 * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec6_t<std::common_type_t<T, U>, Tag>
operator*(T s, MVec6_t<U, Tag> const& v)
{
    return MVec6_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                                  v.c4 * s, v.c5 * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec6_t<std::common_type_t<T, U>, Tag>
operator/(MVec6_t<T, Tag> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                                   std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec6_t<ctype, Tag>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv, v.c4 * inv,
                               v.c5 * inv);
}

////////////////////////////////////////////////////////////////////////////////
// MVec6_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec6_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << "," << v.c4 << ","
       << v.c5 << ")";
    return os;
}

} // namespace hd::ga