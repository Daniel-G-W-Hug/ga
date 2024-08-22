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

#include "ga_type_0d.hpp"
#include "ga_type_3d.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
//
// with the term in brackets being the bivector c1 e2^e3 + c2 e3^e1 + c3 e1^e2
// that models a plane in 3D defined by the bivector coordinates (c1, c2, c3)
//
// This is the definition of a multivector in the even subalgebra of G<3,0,0>,
// i.e. it models only multivectors with even grades 0 and 2
// which are quaternions.
// This subalgebra is closed under addition and multiplication.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.


template <typename T = value_t>
    requires(std::floating_point<T>)
struct MVec3d_E {

    // ctors

    // (all grades = 0)
    MVec3d_E() = default;

    // assign all components
    MVec3d_E(T s, T yz, T zx, T xy) : c0(s), c1(yz), c2(zx), c3(xy) {}

    // assign a scalar part exclusively (other grades = 0)
    MVec3d_E(Scalar<T> s) : c0(s) {}

    // assign a bivector part exclusively (other grades = 0)
    MVec3d_E(BiVec3d<T> b) : c1(b.x), c2(b.y), c3(b.z) {}

    // assign scalar and bivector parts
    MVec3d_E(Scalar<T> s, BiVec3d<T> b) : c0(s), c1(b.x), c2(b.y), c3(b.z) {}

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_E(MVec3d_E<U> const& v) : c0(v.c0), c1(v.c1), c2(v.c2), c3(v.c3)
    {
    }


    T c0{}; // scalar component
    T c1{}; // bivector x component (maps to e2^e3 = yz)
    T c2{}; // bivector y component (maps to e3^e1 = zy)
    T c3{}; // bivector z component (maps to e1^e2 = xy)

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec3d_E<U> const& rhs) const
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
    MVec3d_E& operator+=(MVec3d_E<U> const& v)
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_E& operator-=(MVec3d_E<U> const& v)
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_E& operator*=(U s)
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec3d_E& operator/=(U s)
    {
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        return (*this);
    }
};


////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> operator-(MVec3d_E<T> const& v)
{
    return MVec3d_E<T>(-v.c0, -v.c1, -v.c2, -v.c3);
}

// add multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator+(MVec3d_E<T> const& v1,
                                                              MVec3d_E<U> const& v2)
{
    return MVec3d_E<std::common_type_t<T, U>>(v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2,
                                              v1.c3 + v2.c3);
}

// substract multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator-(MVec3d_E<T> const& v1,
                                                              MVec3d_E<U> const& v2)
{
    return MVec3d_E<std::common_type_t<T, U>>(v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2,
                                              v1.c3 - v2.c3);
}


// multiply a multivector with a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& v, U s)
{
    return MVec3d_E<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(T s, MVec3d_E<U> const& v)
{
    return MVec3d_E<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s);
}

// devide an even multivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator/(MVec3d_E<T> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec3d_E<ctype>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv);
}

// returning various grades of the even multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector

template <typename T> inline constexpr Scalar<T> gr0(MVec3d_E<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr BiVec3d<T> gr2(MVec3d_E<T> const& v)
{
    return BiVec3d<T>(v.c1, v.c2, v.c3);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec3d_E<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << ")";
    return os;
}

} // namespace hd::ga