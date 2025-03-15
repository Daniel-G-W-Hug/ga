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
// regressive reversion operation: reverse of complement + backtrafo
// rrev(A_r) = lcmpl((-1)^(r*(r-1)/2) rcmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=4, r = 0, 1, 2, 3, ...: + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> rrev(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> rrev(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> rrev(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> rrev(TriVec3dp<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> rrev(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> rrev(MVec3dp_E<T> const& E)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec3dp_E<T>(E.c0, -E.c1, -E.c2, -E.c3, -E.c4, -E.c5, -E.c6, E.c7);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> rrev(MVec3dp_U<T> const& U)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3dp_U<T>(-U.c0, -U.c1, -U.c2, -U.c3, U.c4, U.c5, U.c6, U.c7);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> rrev(MVec3dp<T> const& M)
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

// returns 0.0 due to degenerate metric with e4^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> dot(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(-t1.w * t2.w);
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

// returns 0.0 due to degenerate metric with e4^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> rdot([[maybe_unused]] Scalar3dp<T>,
                                                           [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////

// wedge product extended to fully populated multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& A,
                                                       MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c4 + A.c4 * B.c1 + A.c5 * B.c0;
    ctype c6 = A.c0 * B.c6 - A.c2 * B.c4 + A.c4 * B.c2 + A.c6 * B.c0;
    ctype c7 = A.c0 * B.c7 - A.c3 * B.c4 + A.c4 * B.c3 + A.c7 * B.c0;
    ctype c8 = A.c0 * B.c8 + A.c2 * B.c3 - A.c3 * B.c2 + A.c8 * B.c0;
    ctype c9 = A.c0 * B.c9 - A.c1 * B.c3 + A.c3 * B.c1 + A.c9 * B.c0;
    ctype c10 = A.c0 * B.c10 + A.c1 * B.c2 - A.c2 * B.c1 + A.c10 * B.c0;
    ctype c11 = A.c0 * B.c11 - A.c2 * B.c7 + A.c3 * B.c6 + A.c4 * B.c8 + A.c6 * B.c3 -
                A.c7 * B.c2 + A.c8 * B.c4 + A.c11 * B.c0;
    ctype c12 = A.c0 * B.c12 + A.c1 * B.c7 - A.c3 * B.c5 + A.c4 * B.c9 - A.c5 * B.c3 +
                A.c7 * B.c1 + A.c9 * B.c4 + A.c12 * B.c0;
    ctype c13 = A.c0 * B.c13 - A.c1 * B.c6 + A.c2 * B.c5 + A.c4 * B.c10 + A.c5 * B.c2 -
                A.c6 * B.c1 + A.c10 * B.c4 + A.c13 * B.c0;
    ctype c14 = A.c0 * B.c14 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 - A.c8 * B.c1 -
                A.c9 * B.c2 - A.c10 * B.c3 + A.c14 * B.c0;
    ctype c15 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
                A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 - A.c8 * B.c5 - A.c9 * B.c6 -
                A.c10 * B.c7 - A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 +
                A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3dp<std::common_type_t<T, U>> wdg(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * ctype(s));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) * ctype(ps));
}

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
inline constexpr TriVec3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t,
                                                         Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return t * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s,
                                                         TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * t;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B1,
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

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline TriVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(
        -v.y * B.vz + v.z * B.vy + v.w * B.mx, v.x * B.vz - v.z * B.vx + v.w * B.my,
        -v.x * B.vy + v.y * B.vx + v.w * B.mz, -v.x * B.mx - v.y * B.my - v.z * B.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
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

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg for Point2d
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Plane3d<std::common_type_t<T, U>> wdg(Line3d<T> const& l, Point3d<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return Plane3d<ctype>(l.vy * p.z - l.vz * p.y + l.mx, -l.vx * p.z + l.vz * p.x + l.my,
                          l.vx * p.y - l.vy * p.x + l.mz,
                          -l.mx * p.x - l.my * p.y - l.mz * p.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Plane3d<std::common_type_t<T, U>> wdg(Point3d<T> const& p, Line3d<U> const& l)
{
    using ctype = std::common_type_t<T, U>;
    return Plane3d<ctype>(-p.y * l.vz + p.z * l.vy + l.mx, p.x * l.vz - p.z * l.vx + l.my,
                          -p.x * l.vy + p.y * l.vx + l.mz,
                          -p.x * l.mx - p.y * l.my - p.z * l.mz);
}

// wedge product between two points (aka vectors with implicit p.z == 1)
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Line3d<std::common_type_t<T, U>> wdg(Point3d<T> const& p,
                                                      Point3d<U> const& q)
{
    return Line3d<std::common_type_t<T, U>>(q.x - p.x, q.y - p.y, q.z - p.z,
                                            p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z,
                                            p.x * q.y - p.y * q.x);
}

// expand to a new line with goes through point p and is perpendicular to the plane
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Line3d<std::common_type_t<T, U>> expand(Point3d<T> const& point,
                                                         Plane3d<U> const& plane)
{
    return weight_expansion(point3dp{point}, plane);
}

// expand to a new plane with contains the point p and is orthogonal to the line
// => returns a plane (aka a trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Plane3d<std::common_type_t<T, U>> expand(Point3d<T> const& point,
                                                          Line3d<U> const& line)
{
    return weight_expansion(point3dp{point}, line);
}

// expand to a new plane with contains the line and is perpendicular to the plane
// => returns a plane (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Plane3d<std::common_type_t<T, U>> expand(Line3d<T> const& line,
                                                          Plane3d<U> const& plane)
{
    return weight_expansion(line, plane);
}

////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg -> join
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> join(TriVec3dp<T> const& t,
                                                           Vec3dp<U> const& v)
{
    return wdg(t, v);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v,
                                                           TriVec3dp<U> const& t)
{
    return wdg(v, t);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> join(BiVec3dp<T> const& B,
                                                          Vec3dp<U> const& v)
{
    return wdg(B, v);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v,
                                                          BiVec3dp<U> const& B)
{
    return wdg(v, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Plane3d<std::common_type_t<T, U>> join(Line3d<T> const& l,
                                                        Point3d<U> const& p)
{
    return wdg(l, p);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Plane3d<std::common_type_t<T, U>> join(Point3d<T> const& p,
                                                        Line3d<U> const& l)
{
    return wdg(p, l);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v1,
                                                         Vec3dp<U> const& v2)
{
    return wdg(v1, v2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Line3d<std::common_type_t<T, U>> join(Point3d<T> const& p,
                                                       Point3d<U> const& q)
{
    return wdg(p, q);
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
inline constexpr TriVec3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps,
                                                          TriVec3dp<U> const& t)
{
    return TriVec3dp<std::common_type_t<T, U>>(ps * t.x, ps * t.y, ps * t.z, ps * t.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                          PScalar3dp<U> ps)
{
    return Vec3dp<std::common_type_t<T, U>>(t.x * ps, t.y * ps, t.z * ps, t.w * ps);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    return BiVec3dp<std::common_type_t<T, U>>(
        t1.z * t2.y - t1.y * t2.z, t1.x * t2.z - t1.z * t2.x, t1.y * t2.x - t1.x * t2.y,
        t1.x * t2.w - t1.w * t2.x, t1.y * t2.w - t1.w * t2.y, t1.z * t2.w - t1.w * t2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                       BiVec3dp<U> const& B)
{
    return Vec3dp<std::common_type_t<T, U>>(
        -t.y * B.mz + t.z * B.my + t.w * B.vx, t.x * B.mz - t.z * B.mx + t.w * B.vy,
        -t.x * B.my + t.y * B.mx + t.w * B.vz, -t.x * B.vx - t.y * B.vy - t.z * B.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B,
                                                       TriVec3dp<U> const& t)
{
    return Vec3dp<std::common_type_t<T, U>>(
        B.vx * t.w + B.my * t.z - B.mz * t.y, B.vy * t.w - B.mx * t.z + B.mz * t.x,
        B.vz * t.w + B.mx * t.y - B.my * t.x, -B.vx * t.x - B.vy * t.y - B.vz * t.z);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                          Vec3dp<U> const& v)
{
    return Scalar3dp<std::common_type_t<T, U>>(-t.x * v.x - t.y * v.y - t.z * v.z -
                                               t.w * v.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(Vec3dp<T> const& v,
                                                          TriVec3dp<U> const& t)
{
    return Scalar3dp<std::common_type_t<T, U>>(v.x * t.x + v.y * t.y + v.z * t.z +
                                               v.w * t.w);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B1,
                                                          BiVec3dp<U> const& B2)
{
    return Scalar3dp<std::common_type_t<T, U>>(-B1.vx * B2.mx - B1.vy * B2.my -
                                               B1.vz * B2.mz - B1.mx * B2.vx -
                                               B1.my * B2.vy - B1.mz * B2.vz);
}

// required to be present for dist3dp (to complile, even if not used)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>>
rwdg([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] Vec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>>
rwdg([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] BiVec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3dp<std::common_type_t<T, U>>
rwdg([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] Vec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions rwdg -> meet
////////////////////////////////////////////////////////////////////////////////

// intersection of two planes -> return a line
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> meet(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    return rwdg(t1, t2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Line3d<std::common_type_t<T, U>> meet(Plane3d<T> const& p1,
                                                       Plane3d<U> const& p2)
{
    return rwdg(static_cast<TriVec3dp<T> const&>(p1),
                static_cast<TriVec3dp<T> const&>(p2));
}

// intersection of plane and line -> return a point
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> meet(TriVec3dp<T> const& t,
                                                       BiVec2dp<U> const& B)
{
    return rwdg(t, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> meet(BiVec3dp<T> const& B,
                                                       TriVec3dp<U> const& t)
{
    return rwdg(B, t);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> meet(Plane3d<T> const& p,
                                                       Line3d<U> const& l)
{
    return rwdg(static_cast<TriVec3dp<T> const&>(p), static_cast<BiVec3dp<T> const&>(l));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> meet(Line3d<T> const& l,
                                                       Plane3d<U> const& p)
{
    return rwdg(static_cast<BiVec3dp<T> const&>(l), static_cast<TriVec3dp<T> const&>(p));
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

// return commutator product cmt(A,B) of two bivectors A and B (= a bivector)
// cmt(A,B) = 0.5*(AB-BA) = gr2(A * B)
// the commutator product is antisymmetric, i.e. it is zero when a bivector is
// multiplied by itself, i.e. in that case only the dot product remains
// as the symmetric part
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& B1,
                                              BiVec3dp<U> const& B2)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz + B1.mz * B2.vy,
                           B1.vx * B2.mz - B1.vz * B2.mx + B1.mx * B2.vz - B1.mz * B2.vx,
                           -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx,
                           -B1.my * B2.mz + B1.mz * B2.my, B1.mx * B2.mz - B1.mz * B2.mx,
                           -B1.mx * B2.my + B1.my * B2.mx);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////


// geometric product A*B for fully populated 3dp multivector
// REALLY VERY Expensive! - Don't use if you don't have to!
//
// Use equivalent formulae instead for not fully populated multivectors:
// a * b = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector*vector = scalar + bivector)
// A * b = (A>>b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector*vector = vector + trivector)
// a * B = (a<<B) + wdg(a,B) = gr1(aB) + gr3(aB)  (vector*bivector = vector + trivector)
//
// multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                             MVec3dp<U> const& B)
{
    // geometric product of two fully populated 3dp multivectors
    // => due to the degenerate algebra some terms are not present in G<3,0,1> compared to
    // G<4,0,0>
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c8 * B.c8 -
               A.c9 * B.c9 - A.c10 * B.c10 - A.c14 * B.c14;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c10 + A.c3 * B.c9 + A.c8 * B.c14 -
               A.c9 * B.c3 + A.c10 * B.c2 + A.c14 * B.c8;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c10 + A.c2 * B.c0 - A.c3 * B.c8 + A.c8 * B.c3 +
               A.c9 * B.c14 - A.c10 * B.c1 + A.c14 * B.c9;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c0 - A.c8 * B.c2 +
               A.c9 * B.c1 + A.c10 * B.c14 + A.c14 * B.c10;
    ctype c4 = A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 + A.c4 * B.c0 +
               A.c5 * B.c1 + A.c6 * B.c2 + A.c7 * B.c3 - A.c8 * B.c11 - A.c9 * B.c12 -
               A.c10 * B.c13 - A.c11 * B.c8 - A.c12 * B.c9 - A.c13 * B.c10 +
               A.c14 * B.c15 - A.c15 * B.c14;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c4 + A.c2 * B.c13 - A.c3 * B.c12 + A.c4 * B.c1 +
               A.c5 * B.c0 - A.c6 * B.c10 + A.c7 * B.c9 + A.c8 * B.c15 - A.c9 * B.c7 +
               A.c10 * B.c6 + A.c11 * B.c14 - A.c12 * B.c3 + A.c13 * B.c2 -
               A.c14 * B.c11 + A.c15 * B.c8;
    ctype c6 = A.c0 * B.c6 - A.c1 * B.c13 - A.c2 * B.c4 + A.c3 * B.c11 + A.c4 * B.c2 +
               A.c5 * B.c10 + A.c6 * B.c0 - A.c7 * B.c8 + A.c8 * B.c7 + A.c9 * B.c15 -
               A.c10 * B.c5 + A.c11 * B.c3 + A.c12 * B.c14 - A.c13 * B.c1 -
               A.c14 * B.c12 + A.c15 * B.c9;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c12 - A.c2 * B.c11 - A.c3 * B.c4 + A.c4 * B.c3 -
               A.c5 * B.c9 + A.c6 * B.c8 + A.c7 * B.c0 - A.c8 * B.c6 + A.c9 * B.c5 +
               A.c10 * B.c15 - A.c11 * B.c2 + A.c12 * B.c1 + A.c13 * B.c14 -
               A.c14 * B.c13 + A.c15 * B.c10;
    ctype c8 = A.c0 * B.c8 - A.c1 * B.c14 + A.c2 * B.c3 - A.c3 * B.c2 + A.c8 * B.c0 -
               A.c9 * B.c10 + A.c10 * B.c9 - A.c14 * B.c1;
    ctype c9 = A.c0 * B.c9 - A.c1 * B.c3 - A.c2 * B.c14 + A.c3 * B.c1 + A.c8 * B.c10 +
               A.c9 * B.c0 - A.c10 * B.c8 - A.c14 * B.c2;
    ctype c10 = A.c0 * B.c10 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c14 - A.c8 * B.c9 +
                A.c9 * B.c8 + A.c10 * B.c0 - A.c14 * B.c3;
    ctype c11 = A.c0 * B.c11 + A.c1 * B.c15 - A.c2 * B.c7 + A.c3 * B.c6 + A.c4 * B.c8 -
                A.c5 * B.c14 + A.c6 * B.c3 - A.c7 * B.c2 + A.c8 * B.c4 - A.c9 * B.c13 +
                A.c10 * B.c12 + A.c11 * B.c0 - A.c12 * B.c10 + A.c13 * B.c9 +
                A.c14 * B.c5 - A.c15 * B.c1;
    ctype c12 = A.c0 * B.c12 + A.c1 * B.c7 + A.c2 * B.c15 - A.c3 * B.c5 + A.c4 * B.c9 -
                A.c5 * B.c3 - A.c6 * B.c14 + A.c7 * B.c1 + A.c8 * B.c13 + A.c9 * B.c4 -
                A.c10 * B.c11 + A.c11 * B.c10 + A.c12 * B.c0 - A.c13 * B.c8 +
                A.c14 * B.c6 - A.c15 * B.c2;
    ctype c13 = A.c0 * B.c13 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c15 + A.c4 * B.c10 +
                A.c5 * B.c2 - A.c6 * B.c1 - A.c7 * B.c14 - A.c8 * B.c12 + A.c9 * B.c11 +
                A.c10 * B.c4 - A.c11 * B.c9 + A.c12 * B.c8 + A.c13 * B.c0 + A.c14 * B.c7 -
                A.c15 * B.c3;
    ctype c14 = A.c0 * B.c14 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 - A.c8 * B.c1 -
                A.c9 * B.c2 - A.c10 * B.c3 + A.c14 * B.c0;
    ctype c15 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
                A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 - A.c8 * B.c5 - A.c9 * B.c6 -
                A.c10 * B.c7 - A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 +
                A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
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
        Scalar3dp<ctype>(A.c0 * B.c0 - A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6),
        BiVec3dp<ctype>(A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 +
                            A.c4 * B.c7 - A.c5 * B.c3 + A.c6 * B.c2 + A.c7 * B.c4,
                        A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 +
                            A.c4 * B.c3 + A.c5 * B.c7 - A.c6 * B.c1 + A.c7 * B.c5,
                        A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 -
                            A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6,
                        A.c0 * B.c4 + A.c4 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5,
                        A.c0 * B.c5 + A.c4 * B.c6 + A.c5 * B.c0 - A.c6 * B.c4,
                        A.c0 * B.c6 - A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c0),
        PScalar3dp<ctype>(A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 -
                          A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0));
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
        Scalar3dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c7 * B.c7),
        BiVec3dp<ctype>(-A.c0 * B.c3 + A.c1 * B.c6 - A.c2 * B.c5 + A.c3 * B.c0 +
                            A.c4 * B.c7 - A.c5 * B.c2 + A.c6 * B.c1 - A.c7 * B.c4,
                        -A.c0 * B.c6 - A.c1 * B.c3 + A.c2 * B.c4 + A.c3 * B.c1 +
                            A.c4 * B.c2 + A.c5 * B.c7 - A.c6 * B.c0 - A.c7 * B.c5,
                        A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c3 + A.c3 * B.c2 -
                            A.c4 * B.c1 + A.c5 * B.c0 + A.c6 * B.c7 - A.c7 * B.c6,
                        -A.c0 * B.c7 + A.c1 * B.c2 - A.c2 * B.c1 - A.c7 * B.c0,
                        -A.c0 * B.c2 - A.c1 * B.c7 + A.c2 * B.c0 - A.c7 * B.c1,
                        A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c7 - A.c7 * B.c2),
        PScalar3dp<ctype>(A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 -
                          A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3));
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
        Vec3dp<ctype>(A.c0 * B.c0 + A.c4 * B.c7 - A.c5 * B.c2 + A.c6 * B.c1,
                      A.c0 * B.c1 + A.c4 * B.c2 + A.c5 * B.c7 - A.c6 * B.c0,
                      A.c0 * B.c2 - A.c4 * B.c1 + A.c5 * B.c0 + A.c6 * B.c7,
                      A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 -
                          A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7),
        TriVec3dp<ctype>(A.c0 * B.c4 - A.c1 * B.c7 + A.c2 * B.c2 - A.c3 * B.c1 +
                             A.c4 * B.c3 - A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c0,
                         A.c0 * B.c5 - A.c1 * B.c2 - A.c2 * B.c7 + A.c3 * B.c0 +
                             A.c4 * B.c6 + A.c5 * B.c3 - A.c6 * B.c4 - A.c7 * B.c1,
                         A.c0 * B.c6 + A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c7 -
                             A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c3 - A.c7 * B.c2,
                         A.c0 * B.c7 - A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2));
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
        Vec3dp<ctype>(A.c0 * B.c0 - A.c1 * B.c6 + A.c2 * B.c5 + A.c7 * B.c4,
                      A.c0 * B.c6 + A.c1 * B.c0 - A.c2 * B.c4 + A.c7 * B.c5,
                      -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c0 + A.c7 * B.c6,
                      -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0 -
                          A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 + A.c7 * B.c7),
        TriVec3dp<ctype>(A.c0 * B.c7 - A.c1 * B.c3 + A.c2 * B.c2 + A.c3 * B.c4 +
                             A.c4 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1,
                         A.c0 * B.c3 + A.c1 * B.c7 - A.c2 * B.c1 + A.c3 * B.c5 +
                             A.c4 * B.c6 + A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2,
                         -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c7 + A.c3 * B.c6 -
                             A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3,
                         -A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 + A.c7 * B.c0));
}

// geometric product M * B of a multivector M from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& M,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(M.c4 * t.w, M.c5 * t.w, M.c6 * t.w,
                      -M.c4 * t.x - M.c5 * t.y - M.c6 * t.z - M.c7 * t.w),
        TriVec3dp<ctype>(M.c0 * t.x - M.c1 * t.w - M.c5 * t.z + M.c6 * t.y,
                         M.c0 * t.y - M.c2 * t.w + M.c4 * t.z - M.c6 * t.x,
                         M.c0 * t.z - M.c3 * t.w - M.c4 * t.y + M.c5 * t.x, M.c0 * t.w));
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(-M.c4 * B.mx - M.c5 * B.my - M.c6 * B.mz),
                            BiVec3dp<ctype>(M.c0 * B.vx - M.c2 * B.mz + M.c3 * B.my -
                                                M.c5 * B.vz + M.c6 * B.vy + M.c7 * B.mx,
                                            M.c0 * B.vy + M.c1 * B.mz - M.c3 * B.mx +
                                                M.c4 * B.vz - M.c6 * B.vx + M.c7 * B.my,
                                            M.c0 * B.vz - M.c1 * B.my + M.c2 * B.mx -
                                                M.c4 * B.vy + M.c5 * B.vx + M.c7 * B.mz,
                                            M.c0 * B.mx - M.c5 * B.mz + M.c6 * B.my,
                                            M.c0 * B.my + M.c4 * B.mz - M.c6 * B.mx,
                                            M.c0 * B.mz - M.c4 * B.my + M.c5 * B.mx),
                            PScalar3dp<ctype>(-M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz -
                                              M.c4 * B.vx - M.c5 * B.vy - M.c6 * B.vz));
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
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(A.c0 * v.x - A.c5 * v.z + A.c6 * v.y,
                      A.c0 * v.y + A.c4 * v.z - A.c6 * v.x,
                      A.c0 * v.z - A.c4 * v.y + A.c5 * v.x,
                      A.c0 * v.w + A.c1 * v.x + A.c2 * v.y + A.c3 * v.z),
        TriVec3dp<ctype>(A.c2 * v.z - A.c3 * v.y + A.c4 * v.w - A.c7 * v.x,
                         -A.c1 * v.z + A.c3 * v.x + A.c5 * v.w - A.c7 * v.y,
                         A.c1 * v.y - A.c2 * v.x + A.c6 * v.w - A.c7 * v.z,
                         -A.c4 * v.x - A.c5 * v.y - A.c6 * v.z));
}

// geometric product ps * s of a trivector ps and a scalar s
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                                Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) *
                             ctype(s)); // scalar multiplication with a trivector
}

// geometric product s * ps of a scalar s and a trivector ps
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                                PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) *
                             ctype(ps)); // scalar multiplication with a trivector
}

// geometric product of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t1,
                                                               TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(-t1.w * t2.w),
                            BiVec3dp<ctype>(t1.x * t2.w - t1.w * t2.x,
                                            t1.y * t2.w - t1.w * t2.y,
                                            t1.z * t2.w - t1.w * t2.z, 0.0, 0.0, 0.0));
}

// geometric product ps * B of a trivector ps (=3dp pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                               BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(t.w * B.mx, t.w * B.my, t.w * B.mz,
                                          -t.x * B.mx - t.y * B.my - t.z * B.mz),
                            TriVec3dp<ctype>(-t.y * B.mz + t.z * B.my + t.w * B.vx,
                                             t.x * B.mz - t.z * B.mx + t.w * B.vy,
                                             -t.x * B.my + t.y * B.mx + t.w * B.vz, 0.0));
}

// geometric product B * ps of bivector A multiplied by a trivector ps from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(B.mx * t.w, B.my * t.w, B.mz * t.w,
                                          -B.mx * t.x - B.my * t.y - B.mz * t.z),
                            TriVec3dp<ctype>(-B.vx * t.w - B.my * t.z + B.mz * t.y,
                                             -B.vy * t.w + B.mx * t.z - B.mz * t.x,
                                             -B.vz * t.w - B.mx * t.y + B.my * t.x, 0.0));
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
        BiVec3dp<ctype>(-t.y * v.z + t.z * v.y, t.x * v.z - t.z * v.x,
                        -t.x * v.y + t.y * v.x, -t.w * v.x, -t.w * v.y, -t.w * v.z),
        PScalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
        BiVec3dp<ctype>(v.y * t.z - v.z * t.y, -v.x * t.z + v.z * t.x,
                        v.x * t.y - v.y * t.x, -v.x * t.w, -v.y * t.w, -v.z * t.w),
        PScalar3dp<ctype>(v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w));
}

// (geometric) product t * s of trivector t multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                               Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return t * ctype(s);
}

// (geometric) product s * t of a scalar s multiplied from the left to the trivector t
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                               TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * t;
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
    return MVec3dp_E<ctype>(dot(B1, B2), cmt(B1, B2), wdg(B1, B2));
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                               Vec3dp<U> const& v)
{
    // ATTENTION:
    // formula B*v = (B>>v) + wdg(B,v) does NOT hold for a degenerate metric!
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(-B.my * v.z + B.mz * v.y, B.mx * v.z - B.mz * v.x,
                      -B.mx * v.y + B.my * v.x, B.vx * v.x + B.vy * v.y + B.vz * v.z),
        TriVec3dp<ctype>(
            B.vy * v.z - B.vz * v.y + B.mx * v.w, -B.vx * v.z + B.vz * v.x + B.my * v.w,
            B.vx * v.y - B.vy * v.x + B.mz * v.w, -B.mx * v.x - B.my * v.y - B.mz * v.z));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                               BiVec3dp<U> const& B)
{
    // ATTENTION:
    // formula B*v = (v<<B) + wdg(v,B) does NOT hold for a degenerate metric!
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(-v.y * B.mz + v.z * B.my, v.x * B.mz - v.z * B.mx,
                      -v.x * B.my + v.y * B.mx, -v.x * B.vx - v.y * B.vy - v.z * B.vz),
        TriVec3dp<ctype>(-v.y * B.vz + v.z * B.vy + v.w * B.mx,
                         v.x * B.vz - v.z * B.vx + v.w * B.my,
                         -v.x * B.vy + v.y * B.vx + v.w * B.mz,
                         -v.x * B.mx - v.y * B.my - v.z * B.mz));
}

// (geometric) product B * s of bivector B multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                              Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// (geometric) product s * B of a scalar s multiplied from the left to the bivector B
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                              BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
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

// (geometric) product v * s of vector v multiplied with scalar s from the right
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                            Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// (geometric) product s * v of scalar s multiplied from the left to vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                            Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
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

// TODO: for degenerate metric: check whether this "partial inverse" is a
//       meaningful concept at all.
//       (A*A^(-1) = 1) is a requirement, which is only met if the e4 components
//       are not present
//       => only some objects have inverses?


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
    return Vec3dp<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
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
    return BiVec3dp<T>(0.0, 0.0, 0.0, B.mx * inv, B.my * inv, B.mz * inv);
}

// return the multiplicative inverse of the bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> inv(TriVec3dp<T> const& t)
{
    // t^(-1) = rev(t)/|t|^2
    // using rev(t) = (-1)^[k(k-1)/2] B for a k-blade: 3-blade => rev(t) = -t
    // using |t|^2 = gr0(rev(t)*t) = gr0(-t*t) = -gr0(t*t) = -dot(t,t)
    // => B^(-1) = (-t)/(-dot(t,t)) = t/dot(t,t)
    T sq_n = -dot(t, t);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("trivector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = -T(1.0) / sq_n; // negative inverse of squared norm for a bivector
    return TriVec3dp<T>(0.0, 0.0, 0.0, t.w * inv);
}


////////////////////////////////////////////////////////////////////////////////
// 3dp euclidean distance
////////////////////////////////////////////////////////////////////////////////

// returns the euclidean distance between objects as homogeneous magnitude
template <typename arg1, typename arg2>
inline constexpr DualNum3dp<value_t> dist3dp(arg1&& a, arg2&& b)
{
    if (gr(a) + gr(b) == 4) {
        return DualNum3dp<value_t>(rwdg(a, b), weight_nrm(wdg(a, att(b))));
    }
    else {
        return DualNum3dp<value_t>(bulk_nrm(att(wdg(a, b))), weight_nrm(wdg(a, att(b))));
    }
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
inline constexpr MVec3dp_E<std::common_type_t<T, U>> exp(BiVec3dp<T> const& B, U theta)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(std::cos(theta)),
                            normalize(B) * std::sin(theta));
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

// TODO: add equivalent formulas from symbolic manipulations to replace the sandwich
// product by a matrix * vector solution with minimal computational effort
//
// see ga_prdxpr

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_E<std::common_type_t<T, U>> rotor(BiVec3dp<T> const& B, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(std::cos(half_angle)),
                            normalize(B) * std::sin(half_angle));
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
inline constexpr Vec3dp<std::common_type_t<T, U>> rotate_opt(Vec3dp<T> const& v,
                                                             MVec3dp_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = R.c0 * v.x - R.c5 * v.z + R.c6 * v.y;
    ctype k2 = R.c0 * v.y + R.c4 * v.z - R.c6 * v.x;
    ctype k3 = R.c0 * v.z - R.c4 * v.y + R.c5 * v.x;
    ctype k4 = R.c0 * v.w + R.c1 * v.x + R.c2 * v.y + R.c3 * v.z;
    ctype k5 = R.c2 * v.z - R.c3 * v.y + R.c4 * v.w - R.c7 * v.x;
    ctype k6 = -R.c1 * v.z + R.c3 * v.x + R.c5 * v.w - R.c7 * v.y;
    ctype k7 = R.c1 * v.y - R.c2 * v.x + R.c6 * v.w - R.c7 * v.z;
    ctype k8 = -R.c4 * v.x - R.c5 * v.y - R.c6 * v.z;
    return Vec3dp<ctype>(k1 * R.c0 + k2 * R.c6 - k3 * R.c5 - k8 * R.c4,
                         -k1 * R.c6 + k2 * R.c0 + k3 * R.c4 - k8 * R.c5,
                         k1 * R.c5 - k2 * R.c4 + k3 * R.c0 - k8 * R.c6,
                         k1 * R.c1 + k2 * R.c2 + k3 * R.c3 + k4 * R.c0 + k5 * R.c4 +
                             k6 * R.c5 + k7 * R.c6 + k8 * R.c7);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> rotate(BiVec3dp<T> const& B,
                                                           MVec3dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3dp_E<ctype> rr = rev(rotor);
    // MVec3dp_E<ctype> tmp = rotor * B;
    // MVec3dp_E<ctype> res = tmp * rr;

    // scalar and pseudoscalar part of res are 0 due to symmetric product
    // rotor * B * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway
    return BiVec3dp<ctype>(gr2<ctype>(rotor * B * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3dp<std::common_type_t<T, U>> rotate_opt(BiVec3dp<T> const& B,
                                                               MVec3dp_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = -R.c4 * B.mx - R.c5 * B.my - R.c6 * B.mz;
    ctype k2 =
        R.c0 * B.vx - R.c2 * B.mz + R.c3 * B.my - R.c5 * B.vz + R.c6 * B.vy + R.c7 * B.mx;
    ctype k3 =
        R.c0 * B.vy + R.c1 * B.mz - R.c3 * B.mx + R.c4 * B.vz - R.c6 * B.vx + R.c7 * B.my;
    ctype k4 =
        R.c0 * B.vz - R.c1 * B.my + R.c2 * B.mx - R.c4 * B.vy + R.c5 * B.vx + R.c7 * B.mz;
    ctype k5 = R.c0 * B.mx - R.c5 * B.mz + R.c6 * B.my;
    ctype k6 = R.c0 * B.my + R.c4 * B.mz - R.c6 * B.mx;
    ctype k7 = R.c0 * B.mz - R.c4 * B.my + R.c5 * B.mx;
    ctype k8 = -R.c1 * B.mx - R.c2 * B.my - R.c3 * B.mz - R.c4 * B.vx - R.c5 * B.vy -
               R.c6 * B.vz;

    return BiVec3dp<ctype>(k1 * R.c1 + k2 * R.c0 + k3 * R.c6 - k4 * R.c5 + k5 * R.c7 +
                               k6 * R.c3 - k7 * R.c2 - k8 * R.c4,
                           k1 * R.c2 - k2 * R.c6 + k3 * R.c0 + k4 * R.c4 - k5 * R.c3 +
                               k6 * R.c7 + k7 * R.c1 - k8 * R.c5,
                           k1 * R.c3 + k2 * R.c5 - k3 * R.c4 + k4 * R.c0 + k5 * R.c2 -
                               k6 * R.c1 + k7 * R.c7 - k8 * R.c6,
                           -k1 * R.c4 + k5 * R.c0 + k6 * R.c6 - k7 * R.c5,
                           k1 * R.c5 - k5 * R.c6 + k6 * R.c0 + k7 * R.c4,
                           k1 * R.c6 + k5 * R.c5 - k6 * R.c4 + k7 * R.c0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> rotate(TriVec3dp<T> const& t,
                                                            MVec3dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3dp_E<ctype> rr = rev(rotor);
    // MVec3dp_E<ctype> tmp = rotor * t;
    // MVec3dp_E<ctype> res = tmp * rr;

    // vector part of res is 0 due to symmetric product  rotor * t * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the vector part
    // which will be zero anyway
    return TriVec3dp<ctype>(gr3<ctype>(rotor * t * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr TriVec3dp<std::common_type_t<T, U>> rotate_opt(TriVec3dp<T> const& t,
                                                                MVec3dp_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = R.c4 * t.w;
    ctype k2 = R.c5 * t.w;
    ctype k3 = R.c6 * t.w;
    ctype k4 = -R.c4 * t.x - R.c5 * t.y - R.c6 * t.z - R.c7 * t.w;
    ctype k5 = R.c0 * t.x - R.c1 * t.w - R.c5 * t.z + R.c6 * t.y;
    ctype k6 = R.c0 * t.y - R.c2 * t.w + R.c4 * t.z - R.c6 * t.x;
    ctype k7 = R.c0 * t.z - R.c3 * t.w - R.c4 * t.y + R.c5 * t.x;
    ctype k8 = R.c0 * t.w;

    return TriVec3dp<ctype>(k1 * R.c0 + k2 * R.c6 - k3 * R.c5 - k8 * R.c4,
                            -k1 * R.c6 + k2 * R.c0 + k3 * R.c4 - k8 * R.c5,
                            k1 * R.c5 - k2 * R.c4 + k3 * R.c0 - k8 * R.c6,
                            k1 * R.c1 + k2 * R.c2 + k3 * R.c3 + k4 * R.c0 + k5 * R.c4 +
                                k6 * R.c5 + k7 * R.c6 + k8 * R.c7);
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
// left complement:  lcmpl(u) ^ u  = I_3dp = e1^e2^e3^e4
// right complement: u ^ rcmpl(u)  = I_3dp = e1^e2^e3^e4
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs
//
// complement operation with e1^e2^e3^e4 as the pseudoscalar

// right complements

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> rcmpl(Scalar3dp<T> s)
{
    // u ^ rcmpl(u) = e1^e2^e3^e4
    // u = 1:
    //     1 ^ rcmpl(u) = e1^e2^e3^e4 => rcmpl(u) = s e1^e2^e3^e4
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline TriVec3dp<T> rcmpl(Vec3dp<T> const& v)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = v.x e1 + v.y e2 + v.z e3 + v.w e4:
    //     u ^ rcmpl(u) = e1^e2^e3^e4 =>
    //     u = e1 => rcmpl(u) = e423
    //     u = e2 => rcmpl(u) = e431
    //     u = e3 => rcmpl(u) = e412
    //     u = e4 => rcmpl(u) = e321
    return TriVec3dp<T>(v.x, v.y, v.z, v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline BiVec3dp<T> rcmpl(BiVec3dp<T> const& B)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     u ^ rcmpl(u) = e1^e2^e3^e4 =>
    //     u = e41 => rcmpl(u) = -e23
    //     u = e42 => rcmpl(u) = -e31
    //     u = e43 => rcmpl(u) = -e12
    //     u = e23 => rcmpl(u) = -e41
    //     u = e31 => rcmpl(u) = -e42
    //     u = e12 => rcmpl(u) = -e43
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
inline Vec3dp<T> rcmpl(TriVec3dp<T> const& t)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     u ^ rcmpl(u) = e1^e2^e3^e4 =>
    //     u = e423 => rcmpl(u) = -e1
    //     u = e431 => rcmpl(u) = -e2
    //     u = e412 => rcmpl(u) = -e3
    //     u = e321 => rcmpl(u) = -e4
    return Vec3dp<T>(-t.x, -t.y, -t.z, -t.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> rcmpl(PScalar3dp<T> ps)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     e1^e2^e3^e4 ^ rcmpl(u) = e1^e2^e3^e4 => rcmpl(u) = ps 1
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> rcmpl(MVec3dp_U<T> const& M)
{
    // use the component complements directly
    return MVec3dp_U<T>(rcmpl(gr3(M)), rcmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> rcmpl(MVec3dp_E<T> const& M)
{
    // use the component complements directly
    return MVec3dp_E<T>(rcmpl(gr4(M)), rcmpl(gr2(M)), rcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> rcmpl(MVec3dp<T> const& M)
{
    // use the component complements directly
    return MVec3dp<T>(rcmpl(gr4(M)), rcmpl(gr3(M)), rcmpl(gr2(M)), rcmpl(gr1(M)),
                      rcmpl(gr0(M)));
}


// left complements

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> lcmpl(Scalar3dp<T> s)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = 1:
    //     lcmpl(u) ^ 1 = e1^e2^e3^e4 => lcmpl(u) = s e1^e2^e3^e4
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline TriVec3dp<T> lcmpl(Vec3dp<T> const& v)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = v.x e1 + v.y e2 + v.z e3 + v.w e4:
    //     lcmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e1 => lcmpl(u) = -e423
    //     u = e2 => lcmpl(u) = -e431
    //     u = e3 => lcmpl(u) = -e412
    //     u = e4 => lcmpl(u) = -e321
    return TriVec3dp<T>(-v.x, -v.y, -v.z, -v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline BiVec3dp<T> lcmpl(BiVec3dp<T> const& B)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     lcmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e41 => lcmpl(u) = -e23
    //     u = e42 => lcmpl(u) = -e31
    //     u = e43 => lcmpl(u) = -e12
    //     u = e23 => lcmpl(u) = -e41
    //     u = e31 => lcmpl(u) = -e42
    //     u = e12 => lcmpl(u) = -e43
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
inline Vec3dp<T> lcmpl(TriVec3dp<T> const& t)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     lcmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e423 => lcmpl(u) = e1
    //     u = e431 => lcmpl(u) = e2
    //     u = e412 => lcmpl(u) = e3
    //     u = e321 => lcmpl(u) = e4
    return Vec3dp<T>(t.x, t.y, t.z, t.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> lcmpl(PScalar3dp<T> ps)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     lcmpl(u) ^ e1^e2^e3^e4 = e1^e2^e3^e4 => lcmpl(u) = ps 1
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> lcmpl(MVec3dp_U<T> const& M)
{
    // use the component complements directly
    return MVec3dp_U<T>(lcmpl(gr3(M)), lcmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> lcmpl(MVec3dp_E<T> const& M)
{
    // use the component complements directly
    return MVec3dp_E<T>(lcmpl(gr4(M)), lcmpl(gr2(M)), lcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> lcmpl(MVec3dp<T> const& M)
{
    // use the component complements directly
    return MVec3dp<T>(lcmpl(gr4(M)), lcmpl(gr3(M)), lcmpl(gr2(M)), lcmpl(gr1(M)),
                      lcmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
// returns component of v1 parallel to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> project_onto(Vec3dp<T> const& v1,
                                                               Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * Vec3dp<ctype>(inv(v2));
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2)) = v1 - project_onto(v1, v2)
// returns component of v1 perpendicular to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3dp<std::common_type_t<T, U>> reject_from(Vec3dp<T> const& v1,
                                                              Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec3dp<ctype>(gr1(wdg(v1, v2) * inv(v2)));

    return Vec3dp<ctype>(v1 - project_onto(v1, v2));
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