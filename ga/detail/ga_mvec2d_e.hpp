#pragma once

// author: Daniel Hug, 2024

#include <algorithm>
#include <cmath>    // abs, sqrt, acos
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <numbers> // math constants like pi
#include <stdexcept>
#include <string>

#include "ga_value_t.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec2d_E<T> M = c0 + c1 I (with I being the bivector of the plane e1^e2)
//
// This is the definition of a multivector in the even subalgebra of G<2,0,0>,
// i.e. it models only multivectors with even grades 0 and 2 in the plane e1^e2
// (complex numbers).
// This subalgebra is closed under addition and multiplication.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.


template <typename T = value_t>
    requires(std::floating_point<T>)
struct MVec2d_E {

    // ctors

    // (all grades = 0)
    MVec2d_E() = default;

    // assign all components
    MVec2d_E(T s, T ps) : c0(s), c1(ps) {}

    // assign a scalar part exclusively (other grades = 0)
    MVec2d_E(Scalar<T> s) : c0(s) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec2d_E(PScalar2d<T> ps) : c1(ps) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    MVec2d_E(Scalar<T> s, PScalar2d<T> ps) : c0(s), c1(ps) {}

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec2d_E(MVec2d_E<U> const& v) : c0(v.c0), c1(v.c1)
    {
    }


    T c0{}; // scalar component
    T c1{}; // bivector component (2d pseudoscalar)

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec2d_E<U> const& rhs) const
    {
        using ctype = std::common_type_t<T, U>;
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps) return true;
        return false;
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2d_E<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus for multivectors from the even subalgebra
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> operator-(MVec2d_E<T> const& v)
{
    return MVec2d_E<T>(-v.c0, -v.c1);
}

// add multivectors from the even subalgebra
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator+(MVec2d_E<T> const& v1,
                                                              MVec2d_E<U> const& v2)
{
    return MVec2d_E<std::common_type_t<T, U>>(v1.c0 + v2.c0, v1.c1 + v2.c1);
}

// substract multivectors multivectors from the even subalgebra
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator-(MVec2d_E<T> const& v1,
                                                              MVec2d_E<U> const& v2)
{
    return MVec2d_E<std::common_type_t<T, U>>(v1.c0 - v2.c0, v1.c1 - v2.c1);
}


// multiply a multivector multivectors from the even subalgebra with a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& v, U s)
{
    return MVec2d_E<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(T s, MVec2d_E<U> const& v)
{
    return MVec2d_E<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s);
}

// devide an even multivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator/(MVec2d_E<T> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < eps) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec2d_E<ctype>(v.c0 * inv, v.c1 * inv);
}

// returning various grades of a multivector from the even subalgebra
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T> inline constexpr Scalar<T> gr0(MVec2d_E<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr PScalar2d<T> gr2(MVec2d_E<T> const& v)
{
    return PScalar2d<T>(v.c1);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2d_E<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec2d_E<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << ")";
    return os;
}

} // namespace hd::ga