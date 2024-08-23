#pragma once

// author: Daniel Hug, 2024

#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_type_0d.hpp"
#include "ga_type_4d.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> definition of a multivector with 16 components c0, ..., c15
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag>
    requires(std::floating_point<T>)
struct MVec16_t {

    // ctors

    // (all grades = 0)
    MVec16_t() = default;

    // assign all components directly
    MVec16_t(T s, T vx, T vy, T vz, T vw, T bu, T bv, T bw, T bx, T by, T bz, T tx, T ty,
             T tz, T tw, T ps) :
        c0(s), c1(vx), c2(vy), c3(vz), c4(vw), c5(bu), c6(bv), c7(bw), c8(bx), c9(by),
        c10(bz), c11(tx), c12(ty), c13(tz), c14(tw), c15(ps)
    {
    }

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t(MVec16_t<U, Tag> const& v) :
        c0(v.c0), c1(v.c1), c2(v.c2), c3(v.c3), c4(v.c4), c5(v.c5), c6(v.c6), c7(v.c7),
        c8(v.c8), c9(v.c9), c10(v.c10), c11(v.c11), c12(v.c12), c13(v.c13), c14(v.c14),
        c15(v.c15)
    {
    }

    ////////////////////////////////////////////////////////////////////////////
    // ctors for MVec4d<T>
    ////////////////////////////////////////////////////////////////////////////

    // assign a scalar part exclusively (other grades = 0)
    MVec16_t(Scalar<T> s) : c0(s) {}

    // assign a vector part exclusively (other grades = 0)
    MVec16_t(Vec4d<T> const& v) : c1(v.x), c2(v.y), c3(v.z), c4(v.w) {}

    // assign a bivector part exclusively (other grades = 0)
    MVec16_t(BiVec4d<T> const& v) :
        c5(v.c0), c6(v.c1), c7(v.c2), c8(v.c3), c9(v.c4), c10(v.c5)
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    MVec16_t(TriVec4d<T> const& t) : c11(t.x), c12(t.y), c13(t.z), c14(t.w) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec16_t(PScalar4d<T> ps) : c15(ps) {}

    // assign components of an even grade subvector
    MVec16_t(Scalar<T> s, BiVec4d<T> const& v, PScalar4d<T> ps) :
        c0(s), c5(v.c0), c6(v.c1), c7(v.c2), c8(v.c3), c9(v.c4), c10(v.c5), c15(ps)
    {
    }

    // assign components of an uneven grade subvector
    MVec16_t(Vec4d<T> const& v, TriVec4d<T> const& t) :
        c1(v.x), c2(v.y), c3(v.z), c4(v.w), c11(t.x), c12(t.y), c13(t.z), c14(t.x)
    {
    }

    ////////////////////////////////////////////////////////////////////////////
    // component definition
    ////////////////////////////////////////////////////////////////////////////

    T c0{};  //     scalar
             //
    T c1{};  //     vector 4d, 1st component    - maps to basis vector     e1
    T c2{};  //     vector 4d, 2nd component    - maps to basis vector     e2
    T c3{};  //     vector 4d, 3rd component    - maps to basis vector     e3
    T c4{};  //     vector 4d, 4th component    - maps to basis vector     e4
             //
    T c5{};  //   bivector 4d, 1st component    - maps to basis bivector   e4^e1
    T c6{};  //   bivector 4d, 2nd component    - maps to basis bivector   e4^e2
    T c7{};  //   bivector 4d, 3rd component    - maps to basis bivector   e4^e3
    T c8{};  //   bivector 4d, 4th component    - maps to basis bivector   e2^e3
    T c9{};  //   bivector 4d, 5th component    - maps to basis bivector   e3^e1
    T c10{}; //   bivector 4d, 6th component    - maps to basis bivector   e1^e2
             //
    T c11{}; //  trivector 4d, 1st component    - maps to basis trivector  e4^e2^e3
    T c12{}; //  trivector 4d, 2nd component    - maps to basis trivector  e4^e3^e1
    T c13{}; //  trivector 4d, 3rd component    - maps to basis trivector  e4^e1^e2
    T c14{}; //  trivector 4d, 4th component    - maps to basis trivector  e3^e2^e1
             //
    T c15{}; // quadvector 4d = 4d pseudoscalar - maps to basis quadvector e1^e2^e3^e4

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MVec16_t<U, Tag> const& rhs) const
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
        auto abs_delta_c8 = std::abs(rhs.c8 - c8);
        auto abs_delta_c9 = std::abs(rhs.c9 - c9);
        auto abs_delta_c10 = std::abs(rhs.c10 - c10);
        auto abs_delta_c11 = std::abs(rhs.c11 - c11);
        auto abs_delta_c12 = std::abs(rhs.c12 - c12);
        auto abs_delta_c13 = std::abs(rhs.c13 - c13);
        auto abs_delta_c14 = std::abs(rhs.c14 - c14);
        auto abs_delta_c15 = std::abs(rhs.c15 - c15);
        auto constexpr delta_eps =
            ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                         std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps &&
            abs_delta_c4 < delta_eps && abs_delta_c5 < delta_eps &&
            abs_delta_c6 < delta_eps && abs_delta_c7 < delta_eps &&
            abs_delta_c8 < delta_eps && abs_delta_c9 < delta_eps &&
            abs_delta_c10 < delta_eps && abs_delta_c11 < delta_eps &&
            abs_delta_c12 < delta_eps && abs_delta_c13 < delta_eps &&
            abs_delta_c14 < delta_eps && abs_delta_c15 < delta_eps)
            return true;
        return false;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator+=(MVec16_t<U, Tag> const& v)
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        c4 += v.c4;
        c5 += v.c5;
        c6 += v.c6;
        c7 += v.c7;
        c8 += v.c8;
        c9 += v.c9;
        c10 += v.c10;
        c11 += v.c11;
        c12 += v.c12;
        c13 += v.c13;
        c14 += v.c14;
        c15 += v.c15;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator-=(MVec16_t<U, Tag> const& v)
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        c4 -= v.c4;
        c5 -= v.c5;
        c6 -= v.c6;
        c7 -= v.c7;
        c8 -= v.c8;
        c9 -= v.c9;
        c10 -= v.c10;
        c11 -= v.c11;
        c12 -= v.c12;
        c13 -= v.c13;
        c14 -= v.c14;
        c15 -= v.c15;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator*=(U s)
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        c4 *= s;
        c5 *= s;
        c6 *= s;
        c7 *= s;
        c8 *= s;
        c9 *= s;
        c10 *= s;
        c11 *= s;
        c12 *= s;
        c13 *= s;
        c14 *= s;
        c15 *= s;
        return (*this);
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec16_t& operator/=(U s)
    {
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        c4 /= s;
        c5 /= s;
        c6 /= s;
        c7 /= s;
        c8 /= s;
        c9 /= s;
        c10 /= s;
        c11 /= s;
        c12 /= s;
        c13 /= s;
        c14 /= s;
        c15 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr MVec16_t<T, Tag> operator-(MVec16_t<T, Tag> const& v)
{
    return MVec16_t<T, Tag>(-v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, -v.c7, -v.c8,
                            -v.c9, -v.c10, -v.c11, -v.c12, -v.c13, -v.c14, -v.c15);
}

// add multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator+(MVec16_t<T, Tag> const& v1, MVec16_t<U, Tag> const& v2)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2, v1.c3 + v2.c3, v1.c4 + v2.c4,
        v1.c5 + v2.c5, v1.c6 + v2.c6, v1.c7 + v2.c7, v1.c8 + v2.c8, v1.c9 + v2.c9,
        v1.c10 + v2.c10, v1.c11 + v2.c11, v1.c12 + v2.c12, v1.c13 + v2.c13,
        v1.c14 + v2.c14, v1.c15 + v2.c15);
}

// substract multivectors
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator-(MVec16_t<T, Tag> const& v1, MVec16_t<U, Tag> const& v2)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2, v1.c3 - v2.c3, v1.c4 - v2.c4,
        v1.c5 - v2.c5, v1.c6 - v2.c6, v1.c7 - v2.c7, v1.c8 - v2.c8, v1.c9 - v2.c9,
        v1.c10 - v2.c10, v1.c11 - v2.c11, v1.c12 - v2.c12, v1.c13 - v2.c13,
        v1.c14 - v2.c14, v1.c15 - v2.c15);
}

// multiply a multivector with a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator*(MVec16_t<T, Tag> const& v, U s)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s, v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s,
        v.c8 * s, v.c9 * s, v.c10 * s, v.c11 * s, v.c12 * s, v.c13 * s, v.c14 * s,
        v.c15 * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator*(T s, MVec16_t<U, Tag> const& v)
{
    return MVec16_t<std::common_type_t<T, U>, Tag>(
        v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s, v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s,
        v.c8 * s, v.c9 * s, v.c10 * s, v.c11 * s, v.c12 * s, v.c13 * s, v.c14 * s,
        v.c15 * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec16_t<std::common_type_t<T, U>, Tag>
operator/(MVec16_t<T, Tag> const& v, U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < ctype(5.0) * std::max<ctype>(std::numeric_limits<T>::epsilon(),
                                                   std::numeric_limits<U>::epsilon())) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec16_t<ctype, Tag>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv,
                                v.c4 * inv, v.c5 * inv, v.c6 * inv, v.c7 * inv,
                                v.c8 * inv, v.c9 * inv, v.c10 * inv, v.c11 * inv,
                                v.c12 * inv, v.c13 * inv, v.c14 * inv, v.c15 * inv);
}

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector
// grade 4: gr4() - quadvector (= pseudoscalar in 4d)

template <typename T, typename Tag>
inline constexpr Scalar<T> gr0(MVec16_t<T, Tag> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T, typename Tag>
inline constexpr Vec4d<T> gr1(MVec16_t<T, Tag> const& v)
{
    return Vec4d<T>(v.c1, v.c2, v.c3, v.c4);
}

template <typename T, typename Tag>
inline constexpr BiVec3d<T> gr2(MVec16_t<T, Tag> const& v)
{
    return BiVec4d<T>(v.c5, v.c6, v.c7, v.c8, v.c9, v.c10);
}

template <typename T, typename Tag>
inline constexpr TriVec4d<T> gr3(MVec16_t<T, Tag> const& v)
{
    return TriVec4d<T>(v.c11, v.c12, v.c13, v.c14);
}

template <typename T, typename Tag>
inline constexpr PScalar4d<T> gr4(MVec16_t<T, Tag> const& v)
{
    return PScalar4d<T>(v.c15);
}

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MVec16_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << "," << v.c4 << ","
       << v.c5 << "," << v.c6 << "," << v.c7 << "," << v.c8 << "," << v.c9 << "," << v.c10
       << "," << v.c11 << "," << v.c12 << "," << v.c13 << "," << v.c14 << "," << v.c15
       << ")";
    return os;
}

} // namespace hd::ga