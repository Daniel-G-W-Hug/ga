#pragma once

// author: Daniel Hug, 2024

#include <algorithm>
#include <cmath>    // abs, sqrt, acos
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "ga_value_t.hpp"

#include "ga_mvec2d_e.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec2d<T> definition of a multivector
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct MVec2d {

    // ctors

    // (all grades = 0)
    MVec2d() = default;

    // assign all components
    MVec2d(T s, T x, T y, T ps) : c0(s), c1(x), c2(y), c3(ps) {}

    // assign a scalar part exclusively (other grades = 0)
    MVec2d(Scalar<T> s) : c0(s) {}

    // assign a vector part exclusively (other grades = 0)
    MVec2d(Vec2d<T> const& v) : c1(v.x), c2(v.y) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec2d(PScalar2d<T> ps) : c3(ps) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    MVec2d(Scalar<T> s, PScalar2d<T> ps) : c0(s), c3(ps) {}

    // assign from a complex number, i.e. from the even subalgebra
    MVec2d(MVec2d_E<T> v) : c0(v.c0), c3(v.c1) {}

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec2d(MVec2d<U> const& v) : c0(v.c0), c1(v.c1), c2(v.c2), c3(v.c3)
    {
    }


    T c0{}; // scalar
    T c1{}; // vector 2d, 1st component
    T c2{}; // vector 2d, 2nd component
    T c3{}; // 2d pseudoscalar

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec2d<U> const& rhs) const
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
    MVec2d& operator+=(MVec2d<U> const& v)
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2d& operator-=(MVec2d<U> const& v)
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2d& operator*=(U s)
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2d& operator/=(U s)
    {
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> operator-(MVec2d<T> const& v)
{
    return MVec2d<T>(-v.c0, -v.c1, -v.c2, -v.c3);
}

// add multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator+(MVec2d<T> const& v1,
                                                            MVec2d<U> const& v2)
{
    return MVec2d<std::common_type_t<T, U>>(v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2,
                                            v1.c3 + v2.c3);
}

// substract multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator-(MVec2d<T> const& v1,
                                                            MVec2d<U> const& v2)
{
    return MVec2d<std::common_type_t<T, U>>(v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2,
                                            v1.c3 - v2.c3);
}


// multiply a multivector with a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& v, U s)
{
    return MVec2d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(T s, MVec2d<U> const& v)
{
    return MVec2d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

// devide a multivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator/(MVec2d<T> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec2d<ctype>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv);
}

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T> inline constexpr Scalar<T> gr0(MVec2d<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr Vec2d<T> gr1(MVec2d<T> const& v)
{
    return Vec2d<T>(v.c1, v.c2);
}

template <typename T> inline constexpr PScalar2d<T> gr2(MVec2d<T> const& v)
{
    return PScalar2d<T>(v.c3);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec2d<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << ")";
    return os;
}

} // namespace hd::ga