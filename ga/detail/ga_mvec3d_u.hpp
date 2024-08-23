#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::max
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_value_t.hpp"

#include "ga_type_0d.hpp"
#include "ga_type_3d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec3d_U<T> M = (c0 * e1 + c1 * e2 + c2 * e3) + c4 e1^e2^e3
//
// with the term in brackets being the vector c0 * e1 + c1 * e2 + c2 * e3
// and the term c4 e1^e2^e3 being the trivector
//
// This is the definition of a multivector in the uneven subalgebra of G<3,0,0>,
// i.e. it models only multivectors with even grades 1 and 3
// This subalgebra is used to store intermediate results for rotations.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.


template <typename T = value_t>
    requires(std::floating_point<T>)
struct MVec3d_U {

    // ctors

    // (all grades = 0)
    MVec3d_U() = default;

    // assign all components
    MVec3d_U(T x, T y, T z, T ps) : c0(x), c1(y), c2(z), c3(ps) {}

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_U(MVec3d_U<U> const& v) : c0(v.c0), c1(v.c1), c2(v.c2), c3(v.c3)
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec3d_U(Vec3d<T> v) : c0(v.x), c1(v.y), c2(v.z) {}

    // assign a scalar part exclusively (other grades = 0)
    MVec3d_U(PScalar3d<T> ps) : c3(ps) {}

    // assign vector and pseudoscalar parts
    MVec3d_U(Vec3d<T> v, PScalar3d<T> ps) : c0(v.x), c1(v.y), c2(v.z), c3(ps) {}

    ////////////////////////////////////////////////////////////////////////////
    // component definition
    ////////////////////////////////////////////////////////////////////////////

    T c0{}; // vector x component (maps to e1)
    T c1{}; // vector y component (maps to e2)
    T c2{}; // vector z component (maps to e3)
    T c3{}; // pseudoscalar component

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec3d_U<U> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto abs_delta_c2 = std::abs(rhs.c2 - c2);
        auto abs_delta_c3 = std::abs(rhs.c3 - c3);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps)
            return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_U& operator+=(MVec3d_U<U> const& v)
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_U& operator-=(MVec3d_U<U> const& v)
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_U& operator*=(U s)
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_U& operator/=(U s)
    {
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec3d_U<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> operator-(MVec3d_U<T> const& v)
{
    return MVec3d_U<T>(-v.c0, -v.c1, -v.c2, -v.c3);
}

// add multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator+(MVec3d_U<T> const& v1,
                                                              MVec3d_U<U> const& v2)
{
    return MVec3d_U<std::common_type_t<T, U>>(v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2,
                                              v1.c3 + v2.c3);
}

// substract multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator-(MVec3d_U<T> const& v1,
                                                              MVec3d_U<U> const& v2)
{
    return MVec3d_U<std::common_type_t<T, U>>(v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2,
                                              v1.c3 - v2.c3);
}


// multiply a multivector with a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& v, U s)
{
    return MVec3d_U<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(T s, MVec3d_U<U> const& v)
{
    return MVec3d_U<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

// devide a an uneven multivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator/(MVec3d_U<T> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec3d_U<std::common_type_t<T, U>>(v.c0 * inv, v.c1 * inv, v.c2 * inv,
                                              v.c3 * inv);
}

// returning various grades of the uneven multivector
//
// grade 1: gr1() - vector
// grade 2: gr3() - pseudoscalar

template <typename T> inline constexpr Vec3d<T> gr1(MVec3d_U<T> const& v)
{
    return Vec3d<T>(v.c0, v.c1, v.c2);
}

template <typename T> inline constexpr PScalar3d<T> gr3(MVec3d_U<T> const& v)
{
    return PScalar3d<T>(v.c3);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d_U<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec3d_U<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << ")";
    return os;
}

} // namespace hd::ga