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

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "ga_cfg_value_t.hpp"
#include "ga_cfg_vec2d.hpp"


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
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto constexpr delta_eps =
            std::common_type_t<T, U>(5.0) *
            std::max<std::common_type_t<T, U>>(std::numeric_limits<T>::epsilon(),
                                               std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps) return true;
        return false;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, MVec2d_E<U> const& v);
};

// for printing via iostream
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec2d_E<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << ")";
    return os;
}

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

// devide a multivector multivectors from the even subalgebra by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator/(MVec2d_E<T> const& v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return MVec2d_E<std::common_type_t<T, U>>(v.c0 * inv, v.c1 * inv);
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
// MVec2d_E<T> basic operations for complex numbers
//                                  (= multivectors from the even subalgebra)
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of complex number
// |Z|^2 = Z rev(Z) = c0^2 + c1^2
template <typename T> inline T sq_nrm(MVec2d_E<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1;
}

// return magnitude of complex number
template <typename T> inline T nrm(MVec2d_E<T> const& v) { return std::sqrt(sq_nrm(v)); }

// return conjugate complex of a complex number,
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline MVec2d_E<T> rev(MVec2d_E<T> const& v)
{
    return MVec2d_E<T>(v.c0, -v.c1);
}

// return a complex number unitized to nrm(v) == 1.0
template <typename T> inline MVec2d_E<T> unitized(MVec2d_E<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return MVec2d_E<T>(v.c0 * inv, v.c1 * inv);
}

// return the multiplicative inverse of the complex number (inv(z) = 1/sq_nrm(z)*rev(z))
// with rev(z) being the complex conjugate
template <typename T> inline MVec2d_E<T> inv(MVec2d_E<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // inverse of squared norm for a vector
    return inv * rev(v);
}

// return the angle of the complex number w.r.t. the real axis
// range of angle: -pi <= angle <= pi
template <typename T>
    requires(std::floating_point<T>)
inline T angle_to_re(MVec2d_E<T> const& v)
{
    using std::numbers::pi;
    if (v.c0 > 0.0) {
        // quadrant I & IV
        return std::atan(v.c1 / v.c0);
    }
    if (v.c0 < 0.0 && v.c1 >= 0.0) {
        // quadrant II
        return std::atan(v.c1 / v.c0) + pi;
    }
    if (v.c0 < 0.0 && v.c1 < 0.0) {
        // quadrant III
        return std::atan(v.c1 / v.c0) - pi;
    }
    if (v.c0 == 0.0) {
        // on y-axis
        if (v.c1 > 0.0) return pi / 2.0;
        if (v.c1 < 0.0) return -pi / 2.0;
    }
    return 0.0; // zero as input => define 0 as corresponding angle
}

} // namespace hd::ga

////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec2d_E<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec2d_E<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{})", nested(v.c0), nested(v.c1));
    }
};
