#pragma once

// author: Daniel Hug, 2024

#include <algorithm>
#include <cmath>    // abs, sqrt, acos
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "ga_cfg_value_t.hpp"

#include "ga_cfg_vec2d.hpp"

#include "ga_cfg_mvec2d_e.hpp"


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
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto abs_delta_c2 = std::abs(rhs.c2 - c2);
        auto abs_delta_c3 = std::abs(rhs.c3 - c3);
        auto constexpr delta_eps =
            std::common_type_t<T, U>(5.0) *
            std::max<std::common_type_t<T, U>>(std::numeric_limits<T>::epsilon(),
                                               std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps)
            return true;
        return false;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, MVec2d<U> const& v);
};

// for printing via iostream
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec2d<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << ")";
    return os;
}

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
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return MVec2d<std::common_type_t<T, U>>(v.c0 * inv, v.c1 * inv, v.c2 * inv,
                                            v.c3 * inv);
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
// MVec2d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude
// |M|^2 = M rev(M) = (M.c0)^2 + (M.c1)^2 + (M.c2)^2 + (M.c3)^3
template <typename T> inline T sq_nrm(MVec2d<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3;
}

// return magnitude
template <typename T> inline T nrm(MVec2d<T> const& v) { return std::sqrt(sq_nrm(v)); }

// return the reverse
template <typename T> inline MVec2d<T> rev(MVec2d<T> const& v)
{
    return MVec2d<T>(v.c0, v.c1, v.c2, -v.c3);
}

// return the Clifford conjugate
template <typename T> inline MVec2d<T> conj(MVec2d<T> const& v)
{
    return MVec2d<T>(v.c0, -v.c1, -v.c2, -v.c3);
}

// return a multivector unitized to nrm(v) == 1.0
template <typename T> inline MVec2d<T> unitized(MVec2d<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return MVec2d<T>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv);
}

// return the multiplicative inverse of the multivector
// inv(M) = 1/( M*conj(M) ) * conj(M)  with M*conj(M) being a scalar value
template <typename T> inline MVec2d<T> inv(MVec2d<T> const& v)
{
    // from manual calculation of M*conj(M) in 2d:
    T m_conjm = v.c0 * v.c0 + v.c3 * v.c3 - sq_nrm(Vec2d<T>(v.c1, v.c2));
    //
    // alternative, but with slightly more computational effort:
    // T m_conjm = gr0(v * conj(v));
    //
    if (std::abs(m_conjm) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("multivector norm too small for inversion " +
                                 std::to_string(m_conjm) + "\n");
        // example: MVec2D(1,1,1,1) is not invertible
    }
    T inv = 1.0 / m_conjm; // inverse of squared norm for a vector
    return inv * conj(v);
}

} // namespace hd::ga


// ////////////////////////////////////////////////////////////////////////////////
// // printing support via fmt library
// ////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{})", nested(v.c0), nested(v.c1),
                              nested(v.c2), nested(v.c3));
    }
};
