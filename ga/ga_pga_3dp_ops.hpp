#pragma once

// author: Daniel Hug, 2024 & 2025

#include <algorithm> // std::max
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <numbers>   // math constants like pi
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "detail/ga_mvec3dp.hpp" // inclusion of multivector imports all component types
#include "detail/ga_pga_3dp_objects.hpp" // Point3dp, Vector3d, Point3d, Line3dp


namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> gr_inv(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> gr_inv(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> gr_inv(BiVec3dp<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> gr_inv(TriVec3dp<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> gr_inv(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> gr_inv(MVec3dp_E<T> const& E)
{
    // grade 0, 2, and 4: no sign change
    return E;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> gr_inv(MVec3dp_U<T> const& U)
{
    // grade 1 and 3: sign reversal
    return -U;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> gr_inv(MVec3dp<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    // grade 1 and 3: sign reversal
    return MVec3dp<T>(M.c0, -M.c1, -M.c2, -M.c3, -M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                      M.c10, -M.c11, -M.c12, -M.c13, -M.c14, M.c15);
}

////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> rev(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> rev(Vec3dp<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> rev(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> rev(TriVec3dp<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> rev(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> rev(MVec3dp_E<T> const& E)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec3dp_E<T>(E.c0, -E.c1, -E.c2, -E.c3, -E.c4, -E.c5, -E.c6, E.c7);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> rev(MVec3dp_U<T> const& U)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec3dp_U<T>(U.c0, U.c1, U.c2, U.c3, -U.c4, -U.c5, -U.c6, -U.c7);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> rev(MVec3dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    // grade 4: no sign change
    return MVec3dp<T>(M.c0, M.c1, M.c2, M.c3, M.c4, -M.c5, -M.c6, -M.c7, -M.c8, -M.c9,
                      -M.c10, -M.c11, -M.c12, -M.c13, -M.c14, M.c15);
}

////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> conj(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> conj(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> conj(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> conj(TriVec3dp<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> conj(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> conj(MVec3dp_E<T> const& E)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec3dp_E<T>(E.c0, -E.c1, -E.c2, -E.c3, -E.c4, -E.c5, -E.c6, E.c7);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> conj(MVec3dp_U<T> const& U)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3dp_U<T>(-U.c0, -U.c1, -U.c2, -U.c3, U.c4, U.c5, U.c6, U.c7);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> conj(MVec3dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec3dp<T>(M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, -M.c7, -M.c8, -M.c9,
                      -M.c10, M.c11, M.c12, M.c13, M.c14, M.c15);
}


////////////////////////////////////////////////////////////////////////////////
// dot products for 3dp (= defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = +1, e4^2 = 0
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot(MVec3dp<T> const& A,
                                                         MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 -
                            A.c8 * B.c8 - A.c9 * B.c9 - A.c10 * B.c10 - A.c14 * B.c14 +
                            A.c15 * B.c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot(PScalar3dp<T> ps1,
                                                         PScalar3dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(ps1) * ctype(ps2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(-t1.z * t2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot(BiVec3dp<T> const& B1,
                                                         BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(-B1.mx * B2.mx - B1.my * B2.my - B1.mz * B2.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot(Vec3dp<T> const& v1,
                                                         Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive dot products for 3dp (= defined for equal grades exclusively)
//
// rdot(v1,v2) = lcmpl( dot(rcmpl(v1),rcmpl(v2)) )
//
// returns a pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> rdot(PScalar3dp<T> ps1,
                                                           PScalar3dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps1) * ctype(ps2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> rdot(TriVec3dp<T> const& t1,
                                                           TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(t1.x * t2.x + t1.y * t2.y + t1.z * t2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> rdot(BiVec3dp<T> const& B1,
                                                           BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-B1.vx * B2.vx - B1.vy * B2.vy - B1.vz * B2.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> rdot(Vec3dp<T> const& v1,
                                                           Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-v1.w * v2.w);
}

// returns 0.0 due to degenerate metric with e3^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> rdot(Scalar3dp<T> s1,
                                                           Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B1,
                                                         BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz -
                             B1.mx * B2.vx - B1.my * B2.vy - B1.mz * B2.vz);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline TriVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(
        B.vy * v.z - B.vz * v.y + B.mx * v.w, -B.vx * v.z + B.vz * v.x + B.my * v.w,
        B.vx * v.y - B.vy * v.x + B.mz * v.w, -B.mx * v.x - B.my * v.y - B.mz * v.z);
}

// wedge product between a vector a and a bivector B
// => returns a trivector (i.e. the pseudoscalar for 3dp)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline TriVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(
        -v.y * B.vz + v.z * B.vy + v.w * B.mx, v.x * B.vz - v.z * B.vx + v.w * B.my,
        -v.x * B.vy + v.y * B.vx + v.w * B.mz, -v.x * B.mx - v.y * B.my - v.z * B.mz);
}

// wedge product between two vectors
// => returns a bivector
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    return BiVec3dp<std::common_type_t<T, U>>(
        v1.w * v2.x - v1.x * v2.w, v1.w * v2.y - v1.y * v2.w, v1.w * v2.z - v1.z * v2.w,
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}


////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product for complements) and meet operations
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product, just depending on the outer product (wdg)
// as well as the complement and thus the pseudoscalar of the space
// (this definition does NOT connect directly to the geometric product,
// but to the outer product exclusively)
//
//        rwdg(ul, ur) = lcmpl(wdg(rcmpl(ul),rcmpl(ur))) = lcmpl(rcmpl(ul) ^ rcmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t1,
                                             TriVec3dp<U> const& t2)
{
    return BiVec3dp<std::common_type_t<T, U>>(
        t1.z * t2.y - t1.y * t2.z, t1.x * t2.z - t1.z * t2.x, t1.y * t2.x - t1.x * t2.y,
        t1.x * t2.w - t1.w * t2.x, t1.y * t2.w - t1.w * t2.y, t1.z * t2.w - t1.w * t2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t, BiVec3dp<U> const& B)
{
    return Vec3dp<std::common_type_t<T, U>>(
        -t.y * B.mz + t.z * B.my + t.w * B.vx, t.x * B.mz - t.z * B.mx + t.w * B.vy,
        -t.x * B.my + t.y * B.mx + t.w * B.vz, -t.x * B.vx - t.y * B.vy - t.z * B.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B, TriVec3dp<U> const& t)
{
    return Vec3dp<std::common_type_t<T, U>>(
        B.vx * t.w + B.my * t.z - B.mz * t.y, B.vy * t.w - B.mx * t.z + B.mz * t.x,
        B.vz * t.w + B.mx * t.y - B.my * t.x, -B.vx * t.x - B.vy * t.y - B.vz * t.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Scalar3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B1,
                                                BiVec3dp<U> const& B2)
{
    return Scalar3dp<std::common_type_t<T, U>>(-B1.vx * B2.mx - B1.vy * B2.my -
                                               B1.vz * B2.mz - B1.mx * B2.vx -
                                               B1.my * B2.vy - B1.mz * B2.vz);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// L. Dorst: The contraction A << B of an a-blade A onto a b-blade B is a sub-blade
// of B of grade b-a which is perpendicular to A, and linear in both arguments
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator<<(TriVec3dp<T> const& t,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(t.x, t.y, t.z, t.w) * ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B,
                                                               PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(B.mx, B.my, B.mz, B.vx, B.vy, B.vz) * ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                                PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(v.x, v.y, v.z, v.w) * ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(
        B.vy * t.z - B.vz * t.y + B.mx * t.w, -B.vx * t.z + B.vz * t.x + B.my * t.w,
        B.vx * t.y - B.vy * t.x + B.mz * t.w, -B.mx * t.x - B.my * t.y - B.mz * t.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(v.y * t.z - v.z * t.y, -v.x * t.z + v.z * t.x,
                           v.x * t.y - v.y * t.x, -v.x * t.w + v.w * t.x,
                           -v.y * t.w + v.w * t.y, -v.z * t.w + v.w * t.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(
        -v.y * B.mz + v.z * B.my + v.w * B.vx, v.x * B.mz - v.z * B.mx + v.w * B.vy,
        -v.x * B.my + v.y * B.mx + v.w * B.vz, -v.x * B.vx - v.y * B.vy - v.z * B.vz);
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The resulting object lies in A and is perpendicular to B
////////////////////////////////////////////////////////////////////////////////


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                             TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * Vec3dp<ctype>(t.x, t.y, t.z, t.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * BiVec3dp<ctype>(B.mx, B.my, B.mz, B.vx, B.vy, B.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                                Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * TriVec3dp<ctype>(v.x, v.y, v.z, v.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(
        -t.y * B.vz + t.z * B.vy + t.w * B.mx, t.x * B.vz - t.z * B.vx + t.w * B.my,
        -t.x * B.vy + t.y * B.vx + t.w * B.mz, -t.x * B.mx - t.y * B.my - t.z * B.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-t.y * v.z + t.z * v.y, t.x * v.z - t.z * v.x,
                           -t.x * v.y + t.y * v.x, t.x * v.w - t.w * v.x,
                           t.y * v.w - t.w * v.y, t.z * v.w - t.w * v.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator>>(BiVec3dp<T> const& B,
                                                             Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(
        -B.vx * v.w - B.my * v.z + B.mz * v.y, -B.vy * v.w + B.mx * v.z - B.mz * v.x,
        -B.vz * v.w - B.mx * v.y + B.my * v.x, B.vx * v.x + B.vy * v.y + B.vz * v.z);
}

// ////////////////////////////////////////////////////////////////////////////////
// // alternative multivector products (in use instead of contractions)
// ////////////////////////////////////////////////////////////////////////////////

// // return commutator product cmt(A,B) of two bivectors A and B (= a bivector)
// // cmt(A,B) = 0.5*(AB-BA) = gr2(A * B)
// // the commutator product is antisymmetric, i.e. it is zero when a bivector is
// // multiplied by itself, i.e. in that case only the dot product remains
// // as the symmetric part
// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline BiVec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& A, BiVec3dp<U> const&
// B)
// {
//     // this implementation is only valid in an orthonormal basis
//     using ctype = std::common_type_t<T, U>;
//     return BiVec3dp<ctype>(A.z * B.y - A.y * B.z, A.x * B.z - A.z * B.x, ctype(0.0));
// }


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 3dp multivector
// Expensive! - Don't use if you don't have to! (48x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// a * b = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector*vector = scalar + bivector)
// A * b = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector*vector = vector + trivector)
// a * B = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)  (vector*bivector = vector + trivector)
//
// multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                             MVec3dp<U> const& B)
{
    // geometric product of two fully populated 3dp multivectors
    // => due to the degenerate algebra 16 terms are not present in G<2,0,1> compared to
    // G<3,0,0>
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c6 * B.c6;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c6 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c6 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
               A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0 -
               A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 +
               A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c0;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
               A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of a multivector with an even grade multivector
// multivector * even grade multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                             MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c6 * B.c3;
    ctype c1 = A.c1 * B.c0 - A.c2 * B.c3;
    ctype c2 = A.c1 * B.c3 + A.c2 * B.c0;
    ctype c3 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0 - A.c7 * B.c3;
    ctype c4 = A.c0 * B.c1 + A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype c5 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1;
    ctype c6 = A.c0 * B.c3 + A.c6 * B.c0;
    ctype c7 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c7 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of an even grade multivector with a multivector
// even grade multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                             MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c3 * B.c6;
    ctype c1 = A.c0 * B.c1 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c7;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4;
    ctype c6 = A.c0 * B.c6 + A.c3 * B.c0;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// multivector * uneven multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                             MVec3dp_U<U> const& B)
{
    // geometric product of two fully populated 3dp multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c1 * B.c0 + A.c2 * B.c1;
    ctype c1 = A.c0 * B.c0 + A.c6 * B.c1;
    ctype c2 = A.c0 * B.c1 - A.c6 * B.c0;
    ctype c3 = A.c0 * B.c2 - A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c3;
    ctype c4 = A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1 + A.c7 * B.c0;
    ctype c5 = -A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0 + A.c7 * B.c1;
    ctype c6 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype c7 = A.c0 * B.c3 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// uneven multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                             MVec3dp<U> const& B)
{
    // geometric product of two fully populated 3dp multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c1 + A.c1 * B.c2;
    ctype c1 = A.c0 * B.c0 - A.c1 * B.c6;
    ctype c2 = A.c0 * B.c6 + A.c1 * B.c0;
    ctype c3 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c0 - A.c3 * B.c6;
    ctype c4 = A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1;
    ctype c5 = -A.c0 * B.c3 + A.c1 * B.c7 + A.c2 * B.c1 + A.c3 * B.c2;
    ctype c6 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product A * ps of a multivector A with the trivector ps from the right
// multivector * trivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                             PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), -A.c6, A.c1, A.c2,
                          ctype(0.0), A.c0) *
           ctype(ps);
}

// geometric product ps * B of a trivector A (=3dp pseudoscalar) multiplied from the left
// to the multivector B
// trivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                             MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), -B.c6, B.c1,
                                      B.c2, ctype(0.0), B.c0);
}

// geometric product A * B for two multivectors from the even subalgebra (3d case)
// even grade multivector * even grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                               MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
        Scalar3dp<ctype>(A.c0 * B.c0 - A.c3 * B.c3),
        BiVec3dp<ctype>(A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2,
                        A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1,
                        A.c0 * B.c3 + A.c3 * B.c0));
}

// geometric product A * B for two multivectors from the uneven subalgebra (3d case)
// uneven grade multivector * uneven grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                               MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
        Scalar3dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1),
        BiVec3dp<ctype>(A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0,
                        -A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1,
                        A.c0 * B.c1 - A.c1 * B.c0));
}

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a multivector B of the uneven subalgebra
// even grade multivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                               MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(A.c0 * B.c0 + A.c3 * B.c1, A.c0 * B.c1 - A.c3 * B.c0,
                      A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3dp<ctype>(A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2));
}

// geometric product A * B of a multivector A from the uneven subalgebra (3d case)
// with a multivector B of the even subalgebra
// uneven grade multivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                               MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(A.c0 * B.c0 - A.c1 * B.c3, A.c0 * B.c3 + A.c1 * B.c0,
                      -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3dp<ctype>(A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0));
}

// geometric product A * ps of an even multivector A multiplied from the right
// by the trivector ps (=3dp pseudoscalar)
// even grade multivector * trivector => uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<U> const& A,
                                                               PScalar3dp<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(ctype(0.0), ctype(0.0), -A.c3),
                            PScalar3dp<ctype>(A.c0)) *
           ctype(ps);
}

// geometric product ps * B of a trivector A (=3dp pseudoscalar) multiplied from the left
// to the even grade multivector B
// trivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                               MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec3dp_U<ctype>(Vec3dp<ctype>(ctype(0.0), ctype(0.0), -B.c3),
                                        PScalar3dp<ctype>(B.c0));
}

// geometric product A * B of an uneven grade multivector A multiplied from the right
// by the trivector B (=3dp pseudoscalar)
// uneven grade multivector * trivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_U<U> const& A,
                                                               PScalar3dp<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(ctype(0.0)),
                            BiVec3dp<ctype>(A.c0, A.c1, ctype(0.0))) *
           ctype(ps);
}

// geometric product A * B of a trivector A (=3dp pseudoscalar) multiplied from the left
// to the uneven grade multivector B
// trivector * uneven grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                               MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec3dp_E<ctype>(Scalar3dp<ctype>(ctype(0.0)),
                                        BiVec3dp<ctype>(B.c0, B.c1, ctype(0.0)));
}

// geometric product M * B of a multivector M from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& M,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(-M.c3 * B.z),
                            BiVec3dp<ctype>(M.c0 * B.x - M.c2 * B.z + M.c3 * B.y,
                                            M.c0 * B.y + M.c1 * B.z - M.c3 * B.x,
                                            M.c0 * B.z));
}

// geometric product A * M of a bivector A with an even grade multivector M
// from the right
// bivector * even grade multivector => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                               MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(-B.z * M.c3),
                            BiVec3dp<ctype>(B.x * M.c0 - B.y * M.c3 + B.z * M.c2,
                                            B.x * M.c3 + B.y * M.c0 - B.z * M.c1,
                                            B.z * M.c0));
}

// geometric product M * B of an uneven grade multivector M with a bivector B
// uneven grade multivector * bivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& M,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(-M.c1 * B.z, M.c0 * B.z, -M.c0 * B.y + M.c1 * B.x - M.c3 * B.z),
        PScalar3dp<ctype>(M.c0 * B.x + M.c1 * B.y + M.c2 * B.z));
}

// geometric product B * M of a bivector A with an uneven grade multivector B
// bivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                               MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(B.z * M.c1, -B.z * M.c0, -B.x * M.c1 + B.y * M.c0 - B.z * M.c3),
        PScalar3dp<ctype>(B.x * M.c0 + B.y * M.c1 + B.z * M.c2));
}

// geometric product A * v of an even grade multivector A with a vector v
// from the right
// even grade multivector * vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(A.c0 * v.x + A.c3 * v.y,
                                          A.c0 * v.y - A.c3 * v.x,
                                          A.c0 * v.z - A.c1 * v.y + A.c2 * v.x),
                            PScalar3dp<ctype>(A.c1 * v.x + A.c2 * v.y + A.c3 * v.z));
}

// geometric product v * B of a vector v and an even grade multivector B from the left
// vector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                               MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(v.x * B.c0 - v.y * B.c3,
                                          v.x * B.c3 + v.y * B.c0,
                                          -v.x * B.c2 + v.y * B.c1 + v.z * B.c0),
                            PScalar3dp<ctype>(v.x * B.c1 + v.y * B.c2 + v.z * B.c3));
}

// geometric product of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>>
operator*([[maybe_unused]] PScalar3dp<T>, [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0); // trivectors in pga3dp contain e3^2 = 0
}

// geometric product ps * B of a trivector ps (=3dp pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                            BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * Vec3dp<ctype>(ctype(0.0), ctype(0.0), B.z);
}

// geometric product B * ps of bivector A multiplied by a trivector ps from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                            PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(ctype(0.0), ctype(0.0), -B.z) * ctype(ps);
}

// geometric product ps * v of a trivector ps (=3dp pseudoscalar) multiplied from the
// left to the vector v trivector * vector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                              Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * BiVec3dp<ctype>(v.x, v.y, ctype(0.0));
}

// geometric product v * ps of a vector a multiplied with a trivector ps from the right
// vector * trivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                              PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(v.x, v.y, ctype(0.0)) * ctype(ps);
}

// geometric product A * B between two bivectors
// A * B = gr0(A * B) + gr2(A * B)
//
// the full geometric bivector product only exists in >= 4d spaces:
// A * B = A * B + cmt(A,B) + wdg(A,B) = gr0(A * B) + gr2(A * B) + gr4(A * B)
// In 3D we don't have gr4(A * B) and thus only the terms up to grade 3 remain.
// The bivector product AxB = cmt(A,B) = 0.5*(ab-ba)is called the commutator product.
//
// A * B = dot(A,B) + cmt(A,B) + wdg(A,B)  (in 4d and higher dimensional spaces)
// A * B = dot(A,B) + cmt(A,B)             (in 3d)
//
// => bivector * bivector = scalar + bivector = even grade multivector (in 3d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B1,
                                                               BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(dot(B1, B2)), cmt(B1, B2));
}

// geometric product B * v for a bivector B and a vector v
// old: B * v = dot(B,v) + wdg(B,v) = gr1(B * v) + gr3(B * v)
// new: B * v = (B >> v) + wdg(B,v) = gr1(B * v) + gr3(B * v)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => bivector * vector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    // return MVec3dp_U<ctype>(dot(B, v), wdg(B, v));
    return MVec3dp_U<ctype>((B >> v), wdg(B, v));
}

// geometric product v * B for a vector v and a bivector B
// old: v * B = dot(v,B) + wdg(v,B) = gr1(v * B) + gr3(v * B)
// new: v * B = (v << B) + wdg(v,B) = gr1(v * B) + gr3(v * B)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => vector * bivector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // return MVec3dp_U<ctype>(dot(v, B), wdg(v, B));
    return MVec3dp_U<ctype>((v << B), wdg(v, B));
}

// geometric product a * b of two vectors
// a * b = dot(a,b) + wdg(a,b) = gr0(a * b) + gr2(a * b)
//
// HINT: if a full 3d multivector is required as result, it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// vector * vector => even grade multivector (= scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Vec3dp<T> const& a,
                                                               Vec3dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(dot(a, b), wdg(a, b));
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s1,
                                                               Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors
// w.r.t. the geometric product:
// for k-blades: A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ... (from reversion)
////////////////////////////////////////////////////////////////////////////////

// return the multiplicative inverse of the vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> inv(Vec3dp<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    // using |v|^2 = gr0(rev(v)*v) = dot(v,v)
    //
    T sq_v = dot(v, v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_v < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector dot product too small for inversion " +
                                 std::to_string(sq_v) + "\n");
    }
#endif
    T inv = T(1.0) / sq_v; // inverse of squared norm for a vector
    return Vec3dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> inv(BiVec3dp<T> const& B)
{
    // B^(-1) = rev(B)/|B|^2
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    // using |B|^2 = gr0(rev(B)*B) = gr0(-B*B) = -gr0(B*B) = -dot(B,B)
    // => B^(-1) = (-B)/(-dot(B,B)) = B/dot(B,B)
    T sq_n = -dot(B, B);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = -T(1.0) / sq_n; // negative inverse of squared norm for a bivector
    return BiVec3dp<T>(B.x * inv, B.y * inv, B.z * inv);
}


////////////////////////////////////////////////////////////////////////////////
// angle operations
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    ctype nrm_prod = bulk_nrm(v1) * bulk_nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between of a vector and a bivector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v) * nrm(B);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(
        std::clamp(ctype(nrm(dot(v, B))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between of a bivector and a vector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(BiVec3dp<T> const& B, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(B) * nrm(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(
        std::clamp(ctype(nrm(dot(B, v))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(BiVec3dp<T> const& B1,
                                                BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(B1) * nrm(B2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "bivector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(B1, B2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


////////////////////////////////////////////////////////////////////////////////
// 3dp rotation operations
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
// exponential function with bivector as argument for setup of quaternions
// as geometric multivector with a scalar and a bivector part
// MVec3dp_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
//
// quaternion: q = a + b I with I being the bivector in brackets above
//             representing a plane in the algebra G^3
//
//             a rotation in 3D is represented by the plane and the
//             size of the rotation, the later is given by the angle
//             theta, which is the magnitude of the bivector
//
// Inputs:
//         - an arbitray bivector representing the oriented plane of rotation
//           (does not need to be normalized)
//         - a rotation angle
// Output:
//         - a rotor representing the rotation
//
// HINT:     For a rotation around an axis n (n = normalize(Vec3dp<T>))
//           use the bivector B = n*I_3d  => B = Vec3dp<T> * PScalar3dp<T>
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> exp(BiVec3dp<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(std::cos(theta)),
                            normalize(I) * std::sin(theta));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//       1.) an arbitray bivector representing the oriented plane of rotation
//           (does not need to be normalized, defines what is a posive rotation angle)
//       2.) a rotation angle in that plane
// Output:
//           a rotor representing the requested rotation,
//           for applying the sandwich product as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// for a rotation about an axis n (n = normalized vector) choose the ansatz n*B = I_3d
// and multiply both sides with n from the left (remember n*n = |n|^2 = 1)
//
// => choose: B = n*I_3d
//
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> rotor(BiVec3dp<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(std::cos(half_angle)),
                            normalize(I) * std::sin(half_angle));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> rotate(Vec3dp<T> const& v,
                                                         MVec3dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3dp_E<ctype> rr = rev(rotor);
    // MVec3dp_U<ctype> tmp = rotor * v;
    // MVec3dp_U<ctype> res = tmp * rr;

    // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the pseudoscalar part
    // which will be zero anyway
    return Vec3dp<ctype>(gr1<ctype>(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> rotate(BiVec3dp<T> const& v,
                                                           MVec3dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3dp_E<ctype> rr = rev(rotor);
    // MVec3dp_E<ctype> tmp = rotor * v;
    // MVec3dp_E<ctype> res = tmp * rr;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway
    return BiVec3dp<ctype>(gr2<ctype>(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> rotate(MVec3dp<T> const& v,
                                                          MVec3dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(rotor * v * rev(rotor));
}

////////////////////////////////////////////////////////////////////////////////
// 3dp complement operations
// (the concept of complement is defined w.r.t. the outer product)
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^2 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors which are in the k-blade u with
// the basis vectors which are NOT contained in the k-blade u
// and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  l_cmpl(u) ^ u  = I_3dp = e1^e2^e3
// right complement: u ^ r_cmpl(u)  = I_3dp = e1^e2^e3
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs

// complement operation with e1^e2^e3 as the pseudoscalar

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> cmpl(Scalar3dp<T> s)
{
    // u ^ cmpl(u) = e1^e2^e3
    // u = 1:
    //     1 ^ cmpl(u) = e1^e2^e3 => cmpl(u) = s e1^e2^e3
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline BiVec3dp<T> cmpl(Vec3dp<T> const& v)
{
    // u ^ compl(u) = e1^e2^e3
    // u = v.x e1 + v.y e2 + v.z e3:
    //     u ^ cmpl(u) = e1^e2^e3 =>
    //     u = e1 => cmpl(u) = e23
    //     u = e2 => cmpl(u) = e31
    //     u = e3 => cmpl(u) = e12
    return BiVec3dp<T>(v.x, v.y, v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline Vec3dp<T> cmpl(BiVec3dp<T> const& B)
{
    // u ^ compl(u) = e1^e2^e3
    // u = b.x e23 + b.y e31 + b.z e12:
    //     u ^ cmpl(u) = e1^e2^e3 =>
    //     u = e23 => cmpl(u) = e1
    //     u = e31 => cmpl(u) = e2
    //     u = e12 => cmpl(u) = e3
    return Vec3dp<T>(B.x, B.y, B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> cmpl(PScalar3dp<T> ps)
{
    // u ^ compl(u) = e1^e2^e3
    // u = e1^e2^e3:
    //     e1^e2^e3 ^ cmpl(u) = e1^e2^e3 => cmpl(u) = ps 1
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> cmpl(MVec3dp_E<T> const& M)
{
    // use the component complements directly
    return MVec3dp_U<T>(cmpl(gr2(M)), cmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> cmpl(MVec3dp_U<T> const& M)
{
    // use the component complements directly
    return MVec3dp_E<T>(cmpl(gr3(M)), cmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> cmpl(MVec3dp<T> const& M)
{
    // use the component complements directly
    return MVec3dp<T>(cmpl(gr3(M)), cmpl(gr2(M)), cmpl(gr1(M)), cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> project_onto(Vec3dp<T> const& v1,
                                                               Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * Vec3dp<ctype>(inv(v2));
}

// projection of v1 onto v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>>
project_onto_normalized(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * v2; // v2 is already its own reverse
}

// projection of a vector v1 onto a bivector v2
// v_parallel = dot(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> project_onto(Vec3dp<T> const& v,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // // Vec3dp<ctype> a = dot(v, B);
    // Vec3dp<ctype> a = (v << B);
    // BiVec3dp<ctype> Bi = inv(B);
    // // use the formular equivalent to the geometric product to save computational cost
    // // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // // v_parallel = gr1(a * Bi) = dot(a,Bi)
    // // return Vec3dp<ctype>(dot(a, Bi));
    // return Vec3dp<ctype>((a << Bi));

    return Vec3dp<ctype>((v << inv(B)) << B);
    // return Vec3dp<ctype>(gr1((v << inv(B)) * B));
}

// projection of a vector v1 onto a normalized bivector v2
// u_parallel = gr1(dot(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>>
project_onto_normalized(Vec3dp<T> const& v1, BiVec3dp<U> const& v2)
{
    // requires v2 to be normalized

    using ctype = std::common_type_t<T, U>;
    Vec3dp<ctype> a = dot(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3dp<ctype> Bi = -v2;
    // use the formular equivalent to the geometric product to save computational cost
    // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(a * Bi) = dot(a,Bi)
    return Vec3dp<ctype>(dot(a, Bi));
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> reject_from(Vec3dp<T> const& v1,
                                                              Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    BiVec3dp<ctype> B = wdg(v1, v2);
    Vec3dp<ctype> v2_inv = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    // return Vec3dp<ctype>(dot(B, v2_inv));
    return Vec3dp<ctype>(B >> v2_inv);
}

// rejection of vector v1 from a normalized vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>>
reject_from_normalized(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    // requires v2 to be normalized
    using ctype = std::common_type_t<T, U>;
    BiVec3dp<ctype> B = wdg(v1, v2);
    Vec3dp<ctype> v2_inv = v2; // v2 is its own inverse, if normalized
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    // return Vec3dp<ctype>(dot(B, v2_inv));
    return Vec3dp<ctype>(B >> v2_inv);
}

// rejection of vector v1 from a bivector v2
// u_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> reject_from(Vec3dp<T> const& v1,
                                                              BiVec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3dp<ctype> A = wdg(v1, v2);
    BiVec3dp<ctype> B = inv(v2);
    // trivector * bivector = vector
    return A * B;
}

// rejection of vector v1 from a normalized bivector v2
// u_perp = wdg(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>>
reject_from_normalized(Vec3dp<T> const& v1, BiVec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3dp<ctype> a = wdg(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3dp<ctype> B = -v2;
    // trivector * bivector = vector (derived from full geometric product to save
    // costs)
    return Vec3dp<ctype>(-a * B.x, -a * B.y, -a * B.z);
}

// reflect a vector u on a hyperplane B orthogonal to vector b
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n (a line in 2d space)
// orthogonal to vector b: the hyperplane is dual to b (i.e. a one dimensional subspace)
//
// HINT: choose b * B = I_3dp  =>  B = b * I_3dp  (for normalized b)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> reflect_on_hyp(Vec3dp<T> const& u,
                                                                 Vec3dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(gr1(-b * u * inv(b)));
}

// reflect a vector u in an arbitrary bivector, i.e. a plane
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> reflect_on(Vec3dp<T> const& u,
                                                             BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(gr1(-B * u * inv(B)));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling planes)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> reflect_on(BiVec3dp<T> const& UB,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(gr2(B * UB * inv(B)));
}

// reflect a vector u another vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> reflect_on_vec(Vec3dp<T> const& u,
                                                                 Vec3dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(gr1(b * u * inv(b)));
}

} // namespace hd::ga::pga