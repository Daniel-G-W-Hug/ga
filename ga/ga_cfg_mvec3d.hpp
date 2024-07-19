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

#include "ga_cfg_vec3d.hpp"

#include "ga_cfg_bivec3d.hpp"

#include "ga_cfg_mvec3d_e.hpp"
#include "ga_cfg_mvec3d_u.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> definition of a multivector
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct MVec3d {

    // ctors

    // (all grades = 0)
    MVec3d() = default;

    // assign all components directly
    MVec3d(T s, T x, T y, T z, T yz, T zx, T xy, T ps) :
        c0(s), c1(x), c2(y), c3(z), c4(yz), c5(zx), c6(xy), c7(ps)
    {
    }

    // assign a scalar part exclusively (other grades = 0)
    MVec3d(Scalar<T> s) : c0(s) {}

    // assign a vector part exclusively (other grades = 0)
    MVec3d(Vec3d<T> const& v) : c1(v.x), c2(v.y), c3(v.z) {}

    // assign a bivector part exclusively (other grades = 0)
    MVec3d(BiVec3d<T> const& v) : c4(v.x), c5(v.y), c6(v.z) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) or via dot(v1,v2) and cmt(v1,v2) directly
    // (other grades = 0)
    MVec3d(Scalar<T> s, BiVec3d<T> const& v) : c0(s), c4(v.x), c5(v.y), c6(v.z) {}

    // assign from a quaternion, i.e. from the even subalgebra
    MVec3d(MVec3d_E<T> const& v) : c0(v.c0), c4(v.c1), c5(v.c2), c6(v.c3) {}

    // assign from the uneven subalgebra
    MVec3d(MVec3d_U<T> const& v) : c1(v.c0), c2(v.c1), c3(v.c2), c7(v.c3) {}

    // assign a geometric product resulting from a product of a vector and a bivector
    MVec3d(Vec3d<T> const& v, PScalar3d<T> ps) : c1(v.x), c2(v.y), c3(v.z), c7(ps) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec3d(PScalar3d<T> ps) : c7(ps) {}

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec3d(MVec3d<U> const& v) :
        c0(v.c0), c1(v.c1), c2(v.c2), c3(v.c3), c4(v.c4), c5(v.c5), c6(v.c6), c7(v.c7)
    {
    }

    T c0{}; // scalar
    T c1{}; // vector 3d, 1st component   (x)  - maps to basis bivector  e1
    T c2{}; // vector 3d, 2nd component   (y)  - maps to basis bivector  e2
    T c3{}; // vector 3d, 3rd component   (z)  - maps to basis bivector  e3
    T c4{}; // bivector 3d, 1st component (yz) - maps to basis bivector  e2^e3
    T c5{}; // bivector 3d, 2nd component (zx) - maps to basis bivector  e3^e1
    T c6{}; // bivector 3d, 3rd component (xy) - maps to basis bivector  e1^e2
    T c7{}; // trivector 3d = 3d pseudoscalar  - maps to basis trivector e1^e2^e3

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec3d<U> const& rhs) const
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
        auto abs_delta_c6 = std::abs(rhs.c6 - c6);
        auto abs_delta_c7 = std::abs(rhs.c7 - c7);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps &&
            abs_delta_c4 < delta_eps && abs_delta_c5 < delta_eps &&
            abs_delta_c6 < delta_eps && abs_delta_c7 < delta_eps)
            return true;
        return false;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, MVec3d<U> const& v);
};

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> operator-(MVec3d<T> const& v)
{
    return MVec3d<T>(-v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, -v.c7);
}

// add multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator+(MVec3d<T> const& v1,
                                                            MVec3d<U> const& v2)
{
    return MVec3d<std::common_type_t<T, U>>(v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2,
                                            v1.c3 + v2.c3, v1.c4 + v2.c4, v1.c5 + v2.c5,
                                            v1.c6 + v2.c6, v1.c7 + v2.c7);
}

// substract multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator-(MVec3d<T> const& v1,
                                                            MVec3d<U> const& v2)
{
    return MVec3d<std::common_type_t<T, U>>(v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2,
                                            v1.c3 - v2.c3, v1.c4 - v2.c4, v1.c5 - v2.c5,
                                            v1.c6 - v2.c6, v1.c7 - v2.c7);
}

// multiply a multivector with a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& v, U s)
{
    return MVec3d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                            v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(T s, MVec3d<U> const& v)
{
    return MVec3d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                            v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s);
}

// devide a multivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator/(MVec3d<T> const& v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return MVec3d<std::common_type_t<T, U>>(v.c0 * inv, v.c1 * inv, v.c2 * inv,
                                            v.c3 * inv, v.c4 * inv, v.c5 * inv,
                                            v.c6 * inv, v.c7 * inv);
}

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector (= pseudoscalar in 3d)

template <typename T> inline constexpr Scalar<T> gr0(MVec3d<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr Vec3d<T> gr1(MVec3d<T> const& v)
{
    return Vec3d<T>(v.c1, v.c2, v.c3);
}
template <typename T> inline constexpr BiVec3d<T> gr2(MVec3d<T> const& v)
{
    return BiVec3d<T>(v.c4, v.c5, v.c6);
}

template <typename T> inline constexpr PScalar3d<T> gr3(MVec3d<T> const& v)
{
    return PScalar3d<T>(v.c7);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude
// |M|^2 = M rev(M) = (M.c0)^2 + (M.c1)^2 + (M.c2)^2 + (M.c3)^3
//                  + (M.c4)^2 + (M.c5)^2 + (M.c6)^2 + (M.c7)^3
template <typename T> inline T sq_nrm(MVec3d<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3 + v.c4 * v.c4 +
           v.c5 * v.c5 + v.c6 * v.c6 + v.c7 * v.c7;
}

// return magnitude of complex number
template <typename T> inline T nrm(MVec3d<T> const& v) { return std::sqrt(sq_nrm(v)); }


// return the reverse
template <typename T> inline constexpr MVec3d<T> rev(MVec3d<T> const& v)
{
    // only bivector and trivector parts switch signs
    return MVec3d<T>(v.c0, v.c1, v.c2, v.c3, -v.c4, -v.c5, -v.c6, -v.c7);
}

// return the Clifford conjugate
template <typename T> inline constexpr MVec3d<T> conj(MVec3d<T> const& v)
{
    // only vector and bivector parts switch signs
    return MVec3d<T>(v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, v.c7);
}


// return a multivector unitized to nrm(v) == 1.0
template <typename T> inline MVec3d<T> unitized(MVec3d<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return MVec3d<T>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv, v.c4 * inv,
                     v.c5 * inv, v.c6 * inv, v.c7 * inv);
}

// return the multiplicative inverse of the multivector
// inv(M) = 1/( M*conj(M) ) * conj(M)  with M*conj(M) being a scalar value
template <typename T> inline MVec3d<T> inv(MVec3d<T> const& v)
{
    // alternative, but with slightly more computational effort:
    T m_conjm = gr0(v * conj(v));
    if (std::abs(m_conjm) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("multivector norm too small for inversion " +
                                 std::to_string(m_conjm) + "\n");
    }
    T inv = 1.0 / m_conjm; // inverse of squared norm for a vector
    return inv * conj(v);
}


// ATTENTION: there is a left and a right inverse (see paper of Hitzer, Sangwine)

// for printing via iostream
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec3d<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << "," << v.c4 << ","
       << v.c5 << "," << v.c6 << "," << v.c7 << ")";
    return os;
}

} // namespace hd::ga


// ////////////////////////////////////////////////////////////////////////////////
// // printing support via fmt library
// ////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec3d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec3d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{},{},{},{},{})", nested(v.c0),
                              nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                              nested(v.c5), nested(v.c6), nested(v.c7));
    }
};
