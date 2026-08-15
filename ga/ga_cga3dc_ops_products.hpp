#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_cga3dc_ops_basics.hpp"

#include "detail/ga_error_handling.hpp" // inv()/rinv() null checks


namespace hd::ga::cga {

////////////////////////////////////////////////////////////////////////////////
// provides cga3dc product operations:
//
// - dot()                   -> dot product
// - rdot()                  -> regressive dot product
//
// - wdg()                   -> wedge product
// - rwdg()                  -> regressive wedge product
//
// - operator<<()            -> left contraction
// - operator>>()            -> right contraction
//
// - cmt()                   -> commutator product (= asymmetric part of gpr)
// - rcmt()                  -> regressive commutator product
//                              (= asymmetric part of rgpr)
//
// - operator*()             -> geometric product (= gpr)
// - rgpr()                  -> regressive geometric product
//
// Expansions are the wdg-based duals of the corresponding contractions.
//
// - l_expand3dc()           -> left expansion
// - r_expand3dc()           -> right expansion
//
// - inv()                   -> inversion operation (w.r.t. geometric product;
//                              Hitzer-Sangwine n = 5 closed form)
// - rinv()                  -> inversion operation (w.r.t. regressive geometric
//                              product; full type coverage, see the inverses
//                              section: the metric is non-degenerate, so neither
//                              of the pga exclusions applies)
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// dot products (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here the cga3dc vector metric is the non-orthogonal null-pair metric
//
//     e1.e1 = e2.e2 = e3.e3 = +1,  e4.e4 = e5.e5 = 0,  e4.e5 = e5.e4 = -1
//
// (det g = -1, non-degenerate). Via the exomorphism G(a ^ b) = G(a) ^ G(b) the
// extended metric G is defined for all grades; for cga3dc it is a SIGNED
// PERMUTATION that pairs each projective (w-bearing) blade with a round
// (u-bearing) partner blade, so the dot rules are off-diagonal but stay
// single-term (e.g. dot(e4,e5) = -1).
////////////////////////////////////////////////////////////////////////////////

// cga3dc dot :: dot(mv,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(MVec3dc<T> const& A,
                                                  MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(
        A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c5 -
        A.c5 * B.c4 + A.c6 * B.c12 + A.c7 * B.c13 + A.c8 * B.c14 + A.c9 * B.c9 +
        A.c10 * B.c10 + A.c11 * B.c11 + A.c12 * B.c6 + A.c13 * B.c7 + A.c14 * B.c8 -
        A.c15 * B.c15 - A.c16 * B.c16 - A.c17 * B.c17 - A.c18 * B.c18 - A.c19 * B.c22 -
        A.c20 * B.c23 - A.c21 * B.c24 - A.c22 * B.c19 - A.c23 * B.c20 - A.c24 * B.c21 +
        A.c25 * B.c25 - A.c26 * B.c26 - A.c27 * B.c27 - A.c28 * B.c28 + A.c29 * B.c30 +
        A.c30 * B.c29 - A.c31 * B.c31);
}

// cga3dc dot :: dot(mv_e,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(MVec3dc_E<T> const& A,
                                                  MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(A.c0 * B.c0 + A.c1 * B.c7 + A.c2 * B.c8 + A.c3 * B.c9 +
                            A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c1 +
                            A.c8 * B.c2 + A.c9 * B.c3 - A.c10 * B.c10 - A.c11 * B.c11 -
                            A.c12 * B.c12 - A.c13 * B.c13 + A.c14 * B.c15 +
                            A.c15 * B.c14);
}

// cga3dc dot :: dot(mv_u,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(MVec3dc_U<T> const& A,
                                                  MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 -
                            A.c4 * B.c3 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7 -
                            A.c8 * B.c11 - A.c9 * B.c12 - A.c10 * B.c13 - A.c11 * B.c8 -
                            A.c12 * B.c9 - A.c13 * B.c10 + A.c14 * B.c14 - A.c15 * B.c15);
}

// cga3dc dot :: dot(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(PScalar3dc<T> ps1, PScalar3dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga3dc dot :: dot(quadvec,quadvec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(QuadVec3dc<T> const& Q1,
                                                  QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-Q1.x * Q2.x - Q1.y * Q2.y - Q1.z * Q2.z + Q1.w * Q2.u +
                            Q1.u * Q2.w);
}

// cga3dc dot :: dot(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(TriVec3dc<T> const& t1,
                                                  TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(
        -t1.vx * t2.vx - t1.vy * t2.vy - t1.vz * t2.vz - t1.mx * t2.px - t1.my * t2.py -
        t1.mz * t2.pz - t1.px * t2.mx - t1.py * t2.my - t1.pz * t2.mz + t1.pw * t2.pw);
}

// cga3dc dot :: dot(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(BiVec3dc<T> const& B1,
                                                  BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(
        B1.vx * B2.px + B1.vy * B2.py + B1.vz * B2.pz + B1.mx * B2.mx + B1.my * B2.my +
        B1.mz * B2.mz + B1.px * B2.vx + B1.py * B2.vy + B1.pz * B2.vz - B1.pw * B2.pw);
}

// cga3dc dot :: dot(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(Vec3dc<T> const& v1,
                                                  Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z - v1.w * v2.u -
                            v1.u * v2.w);
}

// cga3dc dot :: dot(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> dot(Scalar3dc<T> s1, Scalar3dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// regressive dot products
// (= inner product of the complements, defined for equal grades exclusively)
////////////////////////////////////////////////////////////////////////////////

// cga3dc rdot :: rdot(mv,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(MVec3dc<T> const& A,
                                                    MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(
        -A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c4 * B.c5 +
        A.c5 * B.c4 - A.c6 * B.c12 - A.c7 * B.c13 - A.c8 * B.c14 - A.c9 * B.c9 -
        A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c6 - A.c13 * B.c7 - A.c14 * B.c8 +
        A.c15 * B.c15 + A.c16 * B.c16 + A.c17 * B.c17 + A.c18 * B.c18 + A.c19 * B.c22 +
        A.c20 * B.c23 + A.c21 * B.c24 + A.c22 * B.c19 + A.c23 * B.c20 + A.c24 * B.c21 -
        A.c25 * B.c25 + A.c26 * B.c26 + A.c27 * B.c27 + A.c28 * B.c28 - A.c29 * B.c30 -
        A.c30 * B.c29 + A.c31 * B.c31);
}

// cga3dc rdot :: rdot(mv_e,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(MVec3dc_E<T> const& A,
                                                    MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(-A.c0 * B.c0 - A.c1 * B.c7 - A.c2 * B.c8 - A.c3 * B.c9 -
                             A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c1 -
                             A.c8 * B.c2 - A.c9 * B.c3 + A.c10 * B.c10 + A.c11 * B.c11 +
                             A.c12 * B.c12 + A.c13 * B.c13 - A.c14 * B.c15 -
                             A.c15 * B.c14);
}

// cga3dc rdot :: rdot(mv_u,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(MVec3dc_U<T> const& A,
                                                    MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(-A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c4 +
                             A.c4 * B.c3 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7 +
                             A.c8 * B.c11 + A.c9 * B.c12 + A.c10 * B.c13 + A.c11 * B.c8 +
                             A.c12 * B.c9 + A.c13 * B.c10 - A.c14 * B.c14 +
                             A.c15 * B.c15);
}

// cga3dc rdot :: rdot(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(PScalar3dc<T> ps1, PScalar3dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(ps1) * ctype(ps2));
}

// cga3dc rdot :: rdot(quadvec,quadvec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(QuadVec3dc<T> const& Q1,
                                                    QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(Q1.x * Q2.x + Q1.y * Q2.y + Q1.z * Q2.z - Q1.w * Q2.u -
                             Q1.u * Q2.w);
}

// cga3dc rdot :: rdot(trivec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(TriVec3dc<T> const& t1,
                                                    TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(
        t1.vx * t2.vx + t1.vy * t2.vy + t1.vz * t2.vz + t1.mx * t2.px + t1.my * t2.py +
        t1.mz * t2.pz + t1.px * t2.mx + t1.py * t2.my + t1.pz * t2.mz - t1.pw * t2.pw);
}

// cga3dc rdot :: rdot(bivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(BiVec3dc<T> const& B1,
                                                    BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(
        -B1.vx * B2.px - B1.vy * B2.py - B1.vz * B2.pz - B1.mx * B2.mx - B1.my * B2.my -
        B1.mz * B2.mz - B1.px * B2.vx - B1.py * B2.vy - B1.pz * B2.vz + B1.pw * B2.pw);
}

// cga3dc rdot :: rdot(vec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(Vec3dc<T> const& v1,
                                                    Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(-v1.x * v2.x - v1.y * v2.y - v1.z * v2.z + v1.w * v2.u +
                             v1.u * v2.w);
}

// cga3dc rdot :: rdot(s,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rdot(Scalar3dc<T> s1, Scalar3dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(-ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product)
////////////////////////////////////////////////////////////////////////////////

// cga3dc wdg :: wdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& A, MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c4 + A.c4 * B.c1 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 - A.c2 * B.c4 + A.c4 * B.c2 + A.c7 * B.c0;
    ctype const c8 = A.c0 * B.c8 - A.c3 * B.c4 + A.c4 * B.c3 + A.c8 * B.c0;
    ctype const c9 = A.c0 * B.c9 + A.c2 * B.c3 - A.c3 * B.c2 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c10 - A.c1 * B.c3 + A.c3 * B.c1 + A.c10 * B.c0;
    ctype const c11 = A.c0 * B.c11 + A.c1 * B.c2 - A.c2 * B.c1 + A.c11 * B.c0;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c5 - A.c5 * B.c1 + A.c12 * B.c0;
    ctype const c13 = A.c0 * B.c13 + A.c2 * B.c5 - A.c5 * B.c2 + A.c13 * B.c0;
    ctype const c14 = A.c0 * B.c14 + A.c3 * B.c5 - A.c5 * B.c3 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c15 + A.c4 * B.c5 - A.c5 * B.c4 + A.c15 * B.c0;
    ctype const c16 = A.c0 * B.c16 - A.c1 * B.c15 + A.c4 * B.c12 + A.c5 * B.c6 +
                      A.c6 * B.c5 + A.c12 * B.c4 - A.c15 * B.c1 + A.c16 * B.c0;
    ctype const c17 = A.c0 * B.c17 - A.c2 * B.c15 + A.c4 * B.c13 + A.c5 * B.c7 +
                      A.c7 * B.c5 + A.c13 * B.c4 - A.c15 * B.c2 + A.c17 * B.c0;
    ctype const c18 = A.c0 * B.c18 - A.c3 * B.c15 + A.c4 * B.c14 + A.c5 * B.c8 +
                      A.c8 * B.c5 + A.c14 * B.c4 - A.c15 * B.c3 + A.c18 * B.c0;
    ctype const c19 = A.c0 * B.c19 + A.c2 * B.c14 - A.c3 * B.c13 + A.c5 * B.c9 +
                      A.c9 * B.c5 - A.c13 * B.c3 + A.c14 * B.c2 + A.c19 * B.c0;
    ctype const c20 = A.c0 * B.c20 - A.c1 * B.c14 + A.c3 * B.c12 + A.c5 * B.c10 +
                      A.c10 * B.c5 + A.c12 * B.c3 - A.c14 * B.c1 + A.c20 * B.c0;
    ctype const c21 = A.c0 * B.c21 + A.c1 * B.c13 - A.c2 * B.c12 + A.c5 * B.c11 +
                      A.c11 * B.c5 - A.c12 * B.c2 + A.c13 * B.c1 + A.c21 * B.c0;
    ctype const c22 = A.c0 * B.c22 - A.c2 * B.c8 + A.c3 * B.c7 + A.c4 * B.c9 +
                      A.c7 * B.c3 - A.c8 * B.c2 + A.c9 * B.c4 + A.c22 * B.c0;
    ctype const c23 = A.c0 * B.c23 + A.c1 * B.c8 - A.c3 * B.c6 + A.c4 * B.c10 -
                      A.c6 * B.c3 + A.c8 * B.c1 + A.c10 * B.c4 + A.c23 * B.c0;
    ctype const c24 = A.c0 * B.c24 - A.c1 * B.c7 + A.c2 * B.c6 + A.c4 * B.c11 +
                      A.c6 * B.c2 - A.c7 * B.c1 + A.c11 * B.c4 + A.c24 * B.c0;
    ctype const c25 = A.c0 * B.c25 - A.c1 * B.c9 - A.c2 * B.c10 - A.c3 * B.c11 -
                      A.c9 * B.c1 - A.c10 * B.c2 - A.c11 * B.c3 + A.c25 * B.c0;
    ctype const c26 = A.c0 * B.c26 - A.c2 * B.c18 + A.c3 * B.c17 + A.c4 * B.c19 -
                      A.c5 * B.c22 + A.c7 * B.c14 - A.c8 * B.c13 + A.c9 * B.c15 -
                      A.c13 * B.c8 + A.c14 * B.c7 + A.c15 * B.c9 - A.c17 * B.c3 +
                      A.c18 * B.c2 - A.c19 * B.c4 + A.c22 * B.c5 + A.c26 * B.c0;
    ctype const c27 = A.c0 * B.c27 + A.c1 * B.c18 - A.c3 * B.c16 + A.c4 * B.c20 -
                      A.c5 * B.c23 - A.c6 * B.c14 + A.c8 * B.c12 + A.c10 * B.c15 +
                      A.c12 * B.c8 - A.c14 * B.c6 + A.c15 * B.c10 + A.c16 * B.c3 -
                      A.c18 * B.c1 - A.c20 * B.c4 + A.c23 * B.c5 + A.c27 * B.c0;
    ctype const c28 = A.c0 * B.c28 - A.c1 * B.c17 + A.c2 * B.c16 + A.c4 * B.c21 -
                      A.c5 * B.c24 + A.c6 * B.c13 - A.c7 * B.c12 + A.c11 * B.c15 -
                      A.c12 * B.c7 + A.c13 * B.c6 + A.c15 * B.c11 - A.c16 * B.c2 +
                      A.c17 * B.c1 - A.c21 * B.c4 + A.c24 * B.c5 + A.c28 * B.c0;
    ctype const c29 = A.c0 * B.c29 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                      A.c5 * B.c25 - A.c9 * B.c12 - A.c10 * B.c13 - A.c11 * B.c14 -
                      A.c12 * B.c9 - A.c13 * B.c10 - A.c14 * B.c11 + A.c19 * B.c1 +
                      A.c20 * B.c2 + A.c21 * B.c3 + A.c25 * B.c5 + A.c29 * B.c0;
    ctype const c30 = A.c0 * B.c30 + A.c1 * B.c22 + A.c2 * B.c23 + A.c3 * B.c24 +
                      A.c4 * B.c25 - A.c6 * B.c9 - A.c7 * B.c10 - A.c8 * B.c11 -
                      A.c9 * B.c6 - A.c10 * B.c7 - A.c11 * B.c8 - A.c22 * B.c1 -
                      A.c23 * B.c2 - A.c24 * B.c3 - A.c25 * B.c4 + A.c30 * B.c0;
    ctype const c31 = A.c0 * B.c31 + A.c1 * B.c26 + A.c2 * B.c27 + A.c3 * B.c28 +
                      A.c4 * B.c29 + A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                      A.c8 * B.c21 - A.c9 * B.c16 - A.c10 * B.c17 - A.c11 * B.c18 -
                      A.c12 * B.c22 - A.c13 * B.c23 - A.c14 * B.c24 - A.c15 * B.c25 -
                      A.c16 * B.c9 - A.c17 * B.c10 - A.c18 * B.c11 - A.c19 * B.c6 -
                      A.c20 * B.c7 - A.c21 * B.c8 - A.c22 * B.c12 - A.c23 * B.c13 -
                      A.c24 * B.c14 - A.c25 * B.c15 + A.c26 * B.c1 + A.c27 * B.c2 +
                      A.c28 * B.c3 + A.c29 * B.c4 + A.c30 * B.c5 + A.c31 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& A,
                                                MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = A.c4 * B.c0;
    ctype const c5 = A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c1 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c2 + A.c7 * B.c0;
    ctype const c8 = A.c0 * B.c3 + A.c8 * B.c0;
    ctype const c9 = A.c0 * B.c4 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c5 + A.c10 * B.c0;
    ctype const c11 = A.c0 * B.c6 + A.c11 * B.c0;
    ctype const c12 = A.c0 * B.c7 + A.c12 * B.c0;
    ctype const c13 = A.c0 * B.c8 + A.c13 * B.c0;
    ctype const c14 = A.c0 * B.c9 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c10 + A.c15 * B.c0;
    ctype const c16 = -A.c1 * B.c10 + A.c4 * B.c7 + A.c5 * B.c1 + A.c16 * B.c0;
    ctype const c17 = -A.c2 * B.c10 + A.c4 * B.c8 + A.c5 * B.c2 + A.c17 * B.c0;
    ctype const c18 = -A.c3 * B.c10 + A.c4 * B.c9 + A.c5 * B.c3 + A.c18 * B.c0;
    ctype const c19 = A.c2 * B.c9 - A.c3 * B.c8 + A.c5 * B.c4 + A.c19 * B.c0;
    ctype const c20 = -A.c1 * B.c9 + A.c3 * B.c7 + A.c5 * B.c5 + A.c20 * B.c0;
    ctype const c21 = A.c1 * B.c8 - A.c2 * B.c7 + A.c5 * B.c6 + A.c21 * B.c0;
    ctype const c22 = -A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c4 + A.c22 * B.c0;
    ctype const c23 = A.c1 * B.c3 - A.c3 * B.c1 + A.c4 * B.c5 + A.c23 * B.c0;
    ctype const c24 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c4 * B.c6 + A.c24 * B.c0;
    ctype const c25 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 + A.c25 * B.c0;
    ctype const c26 = A.c0 * B.c11 + A.c7 * B.c9 - A.c8 * B.c8 + A.c9 * B.c10 -
                      A.c13 * B.c3 + A.c14 * B.c2 + A.c15 * B.c4 + A.c26 * B.c0;
    ctype const c27 = A.c0 * B.c12 - A.c6 * B.c9 + A.c8 * B.c7 + A.c10 * B.c10 +
                      A.c12 * B.c3 - A.c14 * B.c1 + A.c15 * B.c5 + A.c27 * B.c0;
    ctype const c28 = A.c0 * B.c13 + A.c6 * B.c8 - A.c7 * B.c7 + A.c11 * B.c10 -
                      A.c12 * B.c2 + A.c13 * B.c1 + A.c15 * B.c6 + A.c28 * B.c0;
    ctype const c29 = A.c0 * B.c14 - A.c9 * B.c7 - A.c10 * B.c8 - A.c11 * B.c9 -
                      A.c12 * B.c4 - A.c13 * B.c5 - A.c14 * B.c6 + A.c29 * B.c0;
    ctype const c30 = A.c0 * B.c15 - A.c6 * B.c4 - A.c7 * B.c5 - A.c8 * B.c6 -
                      A.c9 * B.c1 - A.c10 * B.c2 - A.c11 * B.c3 + A.c30 * B.c0;
    ctype const c31 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 +
                      A.c5 * B.c15 - A.c16 * B.c4 - A.c17 * B.c5 - A.c18 * B.c6 -
                      A.c19 * B.c1 - A.c20 * B.c2 - A.c21 * B.c3 - A.c22 * B.c7 -
                      A.c23 * B.c8 - A.c24 * B.c9 - A.c25 * B.c10 + A.c31 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& A,
                                                MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4;
    ctype const c5 = A.c0 * B.c5;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c2 * B.c0;
    ctype const c8 = A.c0 * B.c8 + A.c3 * B.c0;
    ctype const c9 = A.c0 * B.c9 + A.c4 * B.c0;
    ctype const c10 = A.c0 * B.c10 + A.c5 * B.c0;
    ctype const c11 = A.c0 * B.c11 + A.c6 * B.c0;
    ctype const c12 = A.c0 * B.c12 + A.c7 * B.c0;
    ctype const c13 = A.c0 * B.c13 + A.c8 * B.c0;
    ctype const c14 = A.c0 * B.c14 + A.c9 * B.c0;
    ctype const c15 = A.c0 * B.c15 + A.c10 * B.c0;
    ctype const c16 = A.c0 * B.c16 + A.c1 * B.c5 + A.c7 * B.c4 - A.c10 * B.c1;
    ctype const c17 = A.c0 * B.c17 + A.c2 * B.c5 + A.c8 * B.c4 - A.c10 * B.c2;
    ctype const c18 = A.c0 * B.c18 + A.c3 * B.c5 + A.c9 * B.c4 - A.c10 * B.c3;
    ctype const c19 = A.c0 * B.c19 + A.c4 * B.c5 - A.c8 * B.c3 + A.c9 * B.c2;
    ctype const c20 = A.c0 * B.c20 + A.c5 * B.c5 + A.c7 * B.c3 - A.c9 * B.c1;
    ctype const c21 = A.c0 * B.c21 + A.c6 * B.c5 - A.c7 * B.c2 + A.c8 * B.c1;
    ctype const c22 = A.c0 * B.c22 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c4;
    ctype const c23 = A.c0 * B.c23 - A.c1 * B.c3 + A.c3 * B.c1 + A.c5 * B.c4;
    ctype const c24 = A.c0 * B.c24 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c4;
    ctype const c25 = A.c0 * B.c25 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype const c26 = A.c0 * B.c26 + A.c2 * B.c14 - A.c3 * B.c13 + A.c4 * B.c15 -
                      A.c8 * B.c8 + A.c9 * B.c7 + A.c10 * B.c9 + A.c11 * B.c0;
    ctype const c27 = A.c0 * B.c27 - A.c1 * B.c14 + A.c3 * B.c12 + A.c5 * B.c15 +
                      A.c7 * B.c8 - A.c9 * B.c6 + A.c10 * B.c10 + A.c12 * B.c0;
    ctype const c28 = A.c0 * B.c28 + A.c1 * B.c13 - A.c2 * B.c12 + A.c6 * B.c15 -
                      A.c7 * B.c7 + A.c8 * B.c6 + A.c10 * B.c11 + A.c13 * B.c0;
    ctype const c29 = A.c0 * B.c29 - A.c4 * B.c12 - A.c5 * B.c13 - A.c6 * B.c14 -
                      A.c7 * B.c9 - A.c8 * B.c10 - A.c9 * B.c11 + A.c14 * B.c0;
    ctype const c30 = A.c0 * B.c30 - A.c1 * B.c9 - A.c2 * B.c10 - A.c3 * B.c11 -
                      A.c4 * B.c6 - A.c5 * B.c7 - A.c6 * B.c8 + A.c15 * B.c0;
    ctype const c31 = A.c0 * B.c31 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                      A.c4 * B.c16 - A.c5 * B.c17 - A.c6 * B.c18 - A.c7 * B.c22 -
                      A.c8 * B.c23 - A.c9 * B.c24 - A.c10 * B.c25 + A.c11 * B.c1 +
                      A.c12 * B.c2 + A.c13 * B.c3 + A.c14 * B.c4 + A.c15 * B.c5;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& A,
                                                MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c0 * B.c1;
    ctype const c3 = A.c0 * B.c2;
    ctype const c4 = A.c0 * B.c3;
    ctype const c5 = A.c0 * B.c4;
    ctype const c6 = -A.c1 * B.c3 + A.c4 * B.c0;
    ctype const c7 = -A.c2 * B.c3 + A.c4 * B.c1;
    ctype const c8 = -A.c3 * B.c3 + A.c4 * B.c2;
    ctype const c9 = A.c2 * B.c2 - A.c3 * B.c1;
    ctype const c10 = -A.c1 * B.c2 + A.c3 * B.c0;
    ctype const c11 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype const c12 = A.c1 * B.c4 - A.c5 * B.c0;
    ctype const c13 = A.c2 * B.c4 - A.c5 * B.c1;
    ctype const c14 = A.c3 * B.c4 - A.c5 * B.c2;
    ctype const c15 = A.c4 * B.c4 - A.c5 * B.c3;
    ctype const c16 = A.c0 * B.c5 + A.c6 * B.c4 + A.c12 * B.c3 - A.c15 * B.c0;
    ctype const c17 = A.c0 * B.c6 + A.c7 * B.c4 + A.c13 * B.c3 - A.c15 * B.c1;
    ctype const c18 = A.c0 * B.c7 + A.c8 * B.c4 + A.c14 * B.c3 - A.c15 * B.c2;
    ctype const c19 = A.c0 * B.c8 + A.c9 * B.c4 - A.c13 * B.c2 + A.c14 * B.c1;
    ctype const c20 = A.c0 * B.c9 + A.c10 * B.c4 + A.c12 * B.c2 - A.c14 * B.c0;
    ctype const c21 = A.c0 * B.c10 + A.c11 * B.c4 - A.c12 * B.c1 + A.c13 * B.c0;
    ctype const c22 = A.c0 * B.c11 + A.c7 * B.c2 - A.c8 * B.c1 + A.c9 * B.c3;
    ctype const c23 = A.c0 * B.c12 - A.c6 * B.c2 + A.c8 * B.c0 + A.c10 * B.c3;
    ctype const c24 = A.c0 * B.c13 + A.c6 * B.c1 - A.c7 * B.c0 + A.c11 * B.c3;
    ctype const c25 = A.c0 * B.c14 - A.c9 * B.c0 - A.c10 * B.c1 - A.c11 * B.c2;
    ctype const c26 = -A.c2 * B.c7 + A.c3 * B.c6 + A.c4 * B.c8 - A.c5 * B.c11 -
                      A.c17 * B.c2 + A.c18 * B.c1 - A.c19 * B.c3 + A.c22 * B.c4;
    ctype const c27 = A.c1 * B.c7 - A.c3 * B.c5 + A.c4 * B.c9 - A.c5 * B.c12 +
                      A.c16 * B.c2 - A.c18 * B.c0 - A.c20 * B.c3 + A.c23 * B.c4;
    ctype const c28 = -A.c1 * B.c6 + A.c2 * B.c5 + A.c4 * B.c10 - A.c5 * B.c13 -
                      A.c16 * B.c1 + A.c17 * B.c0 - A.c21 * B.c3 + A.c24 * B.c4;
    ctype const c29 = -A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 - A.c5 * B.c14 +
                      A.c19 * B.c0 + A.c20 * B.c1 + A.c21 * B.c2 + A.c25 * B.c4;
    ctype const c30 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
                      A.c22 * B.c0 - A.c23 * B.c1 - A.c24 * B.c2 - A.c25 * B.c3;
    ctype const c31 = A.c0 * B.c15 - A.c6 * B.c8 - A.c7 * B.c9 - A.c8 * B.c10 -
                      A.c9 * B.c5 - A.c10 * B.c6 - A.c11 * B.c7 - A.c12 * B.c11 -
                      A.c13 * B.c12 - A.c14 * B.c13 - A.c15 * B.c14 + A.c26 * B.c0 +
                      A.c27 * B.c1 + A.c28 * B.c2 + A.c29 * B.c3 + A.c30 * B.c4;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& A,
                                                MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c1 * B.c0;
    ctype const c3 = A.c2 * B.c0;
    ctype const c4 = A.c3 * B.c0;
    ctype const c5 = A.c4 * B.c0;
    ctype const c6 = -A.c0 * B.c4 + A.c3 * B.c1;
    ctype const c7 = -A.c1 * B.c4 + A.c3 * B.c2;
    ctype const c8 = -A.c2 * B.c4 + A.c3 * B.c3;
    ctype const c9 = A.c1 * B.c3 - A.c2 * B.c2;
    ctype const c10 = -A.c0 * B.c3 + A.c2 * B.c1;
    ctype const c11 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c12 = A.c0 * B.c5 - A.c4 * B.c1;
    ctype const c13 = A.c1 * B.c5 - A.c4 * B.c2;
    ctype const c14 = A.c2 * B.c5 - A.c4 * B.c3;
    ctype const c15 = A.c3 * B.c5 - A.c4 * B.c4;
    ctype const c16 = -A.c0 * B.c15 + A.c3 * B.c12 + A.c4 * B.c6 + A.c5 * B.c0;
    ctype const c17 = -A.c1 * B.c15 + A.c3 * B.c13 + A.c4 * B.c7 + A.c6 * B.c0;
    ctype const c18 = -A.c2 * B.c15 + A.c3 * B.c14 + A.c4 * B.c8 + A.c7 * B.c0;
    ctype const c19 = A.c1 * B.c14 - A.c2 * B.c13 + A.c4 * B.c9 + A.c8 * B.c0;
    ctype const c20 = -A.c0 * B.c14 + A.c2 * B.c12 + A.c4 * B.c10 + A.c9 * B.c0;
    ctype const c21 = A.c0 * B.c13 - A.c1 * B.c12 + A.c4 * B.c11 + A.c10 * B.c0;
    ctype const c22 = -A.c1 * B.c8 + A.c2 * B.c7 + A.c3 * B.c9 + A.c11 * B.c0;
    ctype const c23 = A.c0 * B.c8 - A.c2 * B.c6 + A.c3 * B.c10 + A.c12 * B.c0;
    ctype const c24 = -A.c0 * B.c7 + A.c1 * B.c6 + A.c3 * B.c11 + A.c13 * B.c0;
    ctype const c25 = -A.c0 * B.c9 - A.c1 * B.c10 - A.c2 * B.c11 + A.c14 * B.c0;
    ctype const c26 = -A.c1 * B.c18 + A.c2 * B.c17 + A.c3 * B.c19 - A.c4 * B.c22 -
                      A.c6 * B.c3 + A.c7 * B.c2 - A.c8 * B.c4 + A.c11 * B.c5;
    ctype const c27 = A.c0 * B.c18 - A.c2 * B.c16 + A.c3 * B.c20 - A.c4 * B.c23 +
                      A.c5 * B.c3 - A.c7 * B.c1 - A.c9 * B.c4 + A.c12 * B.c5;
    ctype const c28 = -A.c0 * B.c17 + A.c1 * B.c16 + A.c3 * B.c21 - A.c4 * B.c24 -
                      A.c5 * B.c2 + A.c6 * B.c1 - A.c10 * B.c4 + A.c13 * B.c5;
    ctype const c29 = -A.c0 * B.c19 - A.c1 * B.c20 - A.c2 * B.c21 - A.c4 * B.c25 +
                      A.c8 * B.c1 + A.c9 * B.c2 + A.c10 * B.c3 + A.c14 * B.c5;
    ctype const c30 = A.c0 * B.c22 + A.c1 * B.c23 + A.c2 * B.c24 + A.c3 * B.c25 -
                      A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4;
    ctype const c31 = A.c0 * B.c26 + A.c1 * B.c27 + A.c2 * B.c28 + A.c3 * B.c29 +
                      A.c4 * B.c30 - A.c5 * B.c9 - A.c6 * B.c10 - A.c7 * B.c11 -
                      A.c8 * B.c6 - A.c9 * B.c7 - A.c10 * B.c8 - A.c11 * B.c12 -
                      A.c12 * B.c13 - A.c13 * B.c14 - A.c14 * B.c15 + A.c15 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& M, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = M.c0 * ctype(ps);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(PScalar3dc<T> ps, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = ctype(ps) * M.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,quadvec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& M,
                                                QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = M.c0 * Q.x;
    ctype const c27 = M.c0 * Q.y;
    ctype const c28 = M.c0 * Q.z;
    ctype const c29 = M.c0 * Q.w;
    ctype const c30 = M.c0 * Q.u;
    ctype const c31 = M.c1 * Q.x + M.c2 * Q.y + M.c3 * Q.z + M.c4 * Q.w + M.c5 * Q.u;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(quadvec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(QuadVec3dc<T> const& Q,
                                                MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = Q.x * M.c0;
    ctype const c27 = Q.y * M.c0;
    ctype const c28 = Q.z * M.c0;
    ctype const c29 = Q.w * M.c0;
    ctype const c30 = Q.u * M.c0;
    ctype const c31 = Q.x * M.c1 + Q.y * M.c2 + Q.z * M.c3 + Q.w * M.c4 + Q.u * M.c5;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& M,
                                                TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = M.c0 * t.vx;
    ctype const c17 = M.c0 * t.vy;
    ctype const c18 = M.c0 * t.vz;
    ctype const c19 = M.c0 * t.mx;
    ctype const c20 = M.c0 * t.my;
    ctype const c21 = M.c0 * t.mz;
    ctype const c22 = M.c0 * t.px;
    ctype const c23 = M.c0 * t.py;
    ctype const c24 = M.c0 * t.pz;
    ctype const c25 = M.c0 * t.pw;
    ctype const c26 = -M.c2 * t.vz + M.c3 * t.vy + M.c4 * t.mx - M.c5 * t.px;
    ctype const c27 = M.c1 * t.vz - M.c3 * t.vx + M.c4 * t.my - M.c5 * t.py;
    ctype const c28 = -M.c1 * t.vy + M.c2 * t.vx + M.c4 * t.mz - M.c5 * t.pz;
    ctype const c29 = -M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz - M.c5 * t.pw;
    ctype const c30 = M.c1 * t.px + M.c2 * t.py + M.c3 * t.pz + M.c4 * t.pw;
    ctype const c31 = -M.c6 * t.mx - M.c7 * t.my - M.c8 * t.mz - M.c9 * t.vx -
                      M.c10 * t.vy - M.c11 * t.vz - M.c12 * t.px - M.c13 * t.py -
                      M.c14 * t.pz - M.c15 * t.pw;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(TriVec3dc<T> const& t,
                                                MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = t.vx * M.c0;
    ctype const c17 = t.vy * M.c0;
    ctype const c18 = t.vz * M.c0;
    ctype const c19 = t.mx * M.c0;
    ctype const c20 = t.my * M.c0;
    ctype const c21 = t.mz * M.c0;
    ctype const c22 = t.px * M.c0;
    ctype const c23 = t.py * M.c0;
    ctype const c24 = t.pz * M.c0;
    ctype const c25 = t.pw * M.c0;
    ctype const c26 = -t.vy * M.c3 + t.vz * M.c2 - t.mx * M.c4 + t.px * M.c5;
    ctype const c27 = t.vx * M.c3 - t.vz * M.c1 - t.my * M.c4 + t.py * M.c5;
    ctype const c28 = -t.vx * M.c2 + t.vy * M.c1 - t.mz * M.c4 + t.pz * M.c5;
    ctype const c29 = t.mx * M.c1 + t.my * M.c2 + t.mz * M.c3 + t.pw * M.c5;
    ctype const c30 = -t.px * M.c1 - t.py * M.c2 - t.pz * M.c3 - t.pw * M.c4;
    ctype const c31 = -t.vx * M.c9 - t.vy * M.c10 - t.vz * M.c11 - t.mx * M.c6 -
                      t.my * M.c7 - t.mz * M.c8 - t.px * M.c12 - t.py * M.c13 -
                      t.pz * M.c14 - t.pw * M.c15;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& M, BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = M.c0 * B.vx;
    ctype const c7 = M.c0 * B.vy;
    ctype const c8 = M.c0 * B.vz;
    ctype const c9 = M.c0 * B.mx;
    ctype const c10 = M.c0 * B.my;
    ctype const c11 = M.c0 * B.mz;
    ctype const c12 = M.c0 * B.px;
    ctype const c13 = M.c0 * B.py;
    ctype const c14 = M.c0 * B.pz;
    ctype const c15 = M.c0 * B.pw;
    ctype const c16 = -M.c1 * B.pw + M.c4 * B.px + M.c5 * B.vx;
    ctype const c17 = -M.c2 * B.pw + M.c4 * B.py + M.c5 * B.vy;
    ctype const c18 = -M.c3 * B.pw + M.c4 * B.pz + M.c5 * B.vz;
    ctype const c19 = M.c2 * B.pz - M.c3 * B.py + M.c5 * B.mx;
    ctype const c20 = -M.c1 * B.pz + M.c3 * B.px + M.c5 * B.my;
    ctype const c21 = M.c1 * B.py - M.c2 * B.px + M.c5 * B.mz;
    ctype const c22 = -M.c2 * B.vz + M.c3 * B.vy + M.c4 * B.mx;
    ctype const c23 = M.c1 * B.vz - M.c3 * B.vx + M.c4 * B.my;
    ctype const c24 = -M.c1 * B.vy + M.c2 * B.vx + M.c4 * B.mz;
    ctype const c25 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz;
    ctype const c26 = M.c7 * B.pz - M.c8 * B.py + M.c9 * B.pw - M.c13 * B.vz +
                      M.c14 * B.vy + M.c15 * B.mx;
    ctype const c27 = -M.c6 * B.pz + M.c8 * B.px + M.c10 * B.pw + M.c12 * B.vz -
                      M.c14 * B.vx + M.c15 * B.my;
    ctype const c28 = M.c6 * B.py - M.c7 * B.px + M.c11 * B.pw - M.c12 * B.vy +
                      M.c13 * B.vx + M.c15 * B.mz;
    ctype const c29 = -M.c9 * B.px - M.c10 * B.py - M.c11 * B.pz - M.c12 * B.mx -
                      M.c13 * B.my - M.c14 * B.mz;
    ctype const c30 = -M.c6 * B.mx - M.c7 * B.my - M.c8 * B.mz - M.c9 * B.vx -
                      M.c10 * B.vy - M.c11 * B.vz;
    ctype const c31 = -M.c16 * B.mx - M.c17 * B.my - M.c18 * B.mz - M.c19 * B.vx -
                      M.c20 * B.vy - M.c21 * B.vz - M.c22 * B.px - M.c23 * B.py -
                      M.c24 * B.pz - M.c25 * B.pw;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(BiVec3dc<T> const& B, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = B.vx * M.c0;
    ctype const c7 = B.vy * M.c0;
    ctype const c8 = B.vz * M.c0;
    ctype const c9 = B.mx * M.c0;
    ctype const c10 = B.my * M.c0;
    ctype const c11 = B.mz * M.c0;
    ctype const c12 = B.px * M.c0;
    ctype const c13 = B.py * M.c0;
    ctype const c14 = B.pz * M.c0;
    ctype const c15 = B.pw * M.c0;
    ctype const c16 = B.vx * M.c5 + B.px * M.c4 - B.pw * M.c1;
    ctype const c17 = B.vy * M.c5 + B.py * M.c4 - B.pw * M.c2;
    ctype const c18 = B.vz * M.c5 + B.pz * M.c4 - B.pw * M.c3;
    ctype const c19 = B.mx * M.c5 - B.py * M.c3 + B.pz * M.c2;
    ctype const c20 = B.my * M.c5 + B.px * M.c3 - B.pz * M.c1;
    ctype const c21 = B.mz * M.c5 - B.px * M.c2 + B.py * M.c1;
    ctype const c22 = B.vy * M.c3 - B.vz * M.c2 + B.mx * M.c4;
    ctype const c23 = -B.vx * M.c3 + B.vz * M.c1 + B.my * M.c4;
    ctype const c24 = B.vx * M.c2 - B.vy * M.c1 + B.mz * M.c4;
    ctype const c25 = -B.mx * M.c1 - B.my * M.c2 - B.mz * M.c3;
    ctype const c26 = B.vy * M.c14 - B.vz * M.c13 + B.mx * M.c15 - B.py * M.c8 +
                      B.pz * M.c7 + B.pw * M.c9;
    ctype const c27 = -B.vx * M.c14 + B.vz * M.c12 + B.my * M.c15 + B.px * M.c8 -
                      B.pz * M.c6 + B.pw * M.c10;
    ctype const c28 = B.vx * M.c13 - B.vy * M.c12 + B.mz * M.c15 - B.px * M.c7 +
                      B.py * M.c6 + B.pw * M.c11;
    ctype const c29 = -B.mx * M.c12 - B.my * M.c13 - B.mz * M.c14 - B.px * M.c9 -
                      B.py * M.c10 - B.pz * M.c11;
    ctype const c30 = -B.vx * M.c9 - B.vy * M.c10 - B.vz * M.c11 - B.mx * M.c6 -
                      B.my * M.c7 - B.mz * M.c8;
    ctype const c31 = -B.vx * M.c19 - B.vy * M.c20 - B.vz * M.c21 - B.mx * M.c16 -
                      B.my * M.c17 - B.mz * M.c18 - B.px * M.c22 - B.py * M.c23 -
                      B.pz * M.c24 - B.pw * M.c25;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& M, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = M.c0 * v.w;
    ctype const c5 = M.c0 * v.u;
    ctype const c6 = -M.c1 * v.w + M.c4 * v.x;
    ctype const c7 = -M.c2 * v.w + M.c4 * v.y;
    ctype const c8 = -M.c3 * v.w + M.c4 * v.z;
    ctype const c9 = M.c2 * v.z - M.c3 * v.y;
    ctype const c10 = -M.c1 * v.z + M.c3 * v.x;
    ctype const c11 = M.c1 * v.y - M.c2 * v.x;
    ctype const c12 = M.c1 * v.u - M.c5 * v.x;
    ctype const c13 = M.c2 * v.u - M.c5 * v.y;
    ctype const c14 = M.c3 * v.u - M.c5 * v.z;
    ctype const c15 = M.c4 * v.u - M.c5 * v.w;
    ctype const c16 = M.c6 * v.u + M.c12 * v.w - M.c15 * v.x;
    ctype const c17 = M.c7 * v.u + M.c13 * v.w - M.c15 * v.y;
    ctype const c18 = M.c8 * v.u + M.c14 * v.w - M.c15 * v.z;
    ctype const c19 = M.c9 * v.u - M.c13 * v.z + M.c14 * v.y;
    ctype const c20 = M.c10 * v.u + M.c12 * v.z - M.c14 * v.x;
    ctype const c21 = M.c11 * v.u - M.c12 * v.y + M.c13 * v.x;
    ctype const c22 = M.c7 * v.z - M.c8 * v.y + M.c9 * v.w;
    ctype const c23 = -M.c6 * v.z + M.c8 * v.x + M.c10 * v.w;
    ctype const c24 = M.c6 * v.y - M.c7 * v.x + M.c11 * v.w;
    ctype const c25 = -M.c9 * v.x - M.c10 * v.y - M.c11 * v.z;
    ctype const c26 = -M.c17 * v.z + M.c18 * v.y - M.c19 * v.w + M.c22 * v.u;
    ctype const c27 = M.c16 * v.z - M.c18 * v.x - M.c20 * v.w + M.c23 * v.u;
    ctype const c28 = -M.c16 * v.y + M.c17 * v.x - M.c21 * v.w + M.c24 * v.u;
    ctype const c29 = M.c19 * v.x + M.c20 * v.y + M.c21 * v.z + M.c25 * v.u;
    ctype const c30 = -M.c22 * v.x - M.c23 * v.y - M.c24 * v.z - M.c25 * v.w;
    ctype const c31 = M.c26 * v.x + M.c27 * v.y + M.c28 * v.z + M.c29 * v.w + M.c30 * v.u;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = v.w * M.c0;
    ctype const c5 = v.u * M.c0;
    ctype const c6 = -v.x * M.c4 + v.w * M.c1;
    ctype const c7 = -v.y * M.c4 + v.w * M.c2;
    ctype const c8 = -v.z * M.c4 + v.w * M.c3;
    ctype const c9 = v.y * M.c3 - v.z * M.c2;
    ctype const c10 = -v.x * M.c3 + v.z * M.c1;
    ctype const c11 = v.x * M.c2 - v.y * M.c1;
    ctype const c12 = v.x * M.c5 - v.u * M.c1;
    ctype const c13 = v.y * M.c5 - v.u * M.c2;
    ctype const c14 = v.z * M.c5 - v.u * M.c3;
    ctype const c15 = v.w * M.c5 - v.u * M.c4;
    ctype const c16 = -v.x * M.c15 + v.w * M.c12 + v.u * M.c6;
    ctype const c17 = -v.y * M.c15 + v.w * M.c13 + v.u * M.c7;
    ctype const c18 = -v.z * M.c15 + v.w * M.c14 + v.u * M.c8;
    ctype const c19 = v.y * M.c14 - v.z * M.c13 + v.u * M.c9;
    ctype const c20 = -v.x * M.c14 + v.z * M.c12 + v.u * M.c10;
    ctype const c21 = v.x * M.c13 - v.y * M.c12 + v.u * M.c11;
    ctype const c22 = -v.y * M.c8 + v.z * M.c7 + v.w * M.c9;
    ctype const c23 = v.x * M.c8 - v.z * M.c6 + v.w * M.c10;
    ctype const c24 = -v.x * M.c7 + v.y * M.c6 + v.w * M.c11;
    ctype const c25 = -v.x * M.c9 - v.y * M.c10 - v.z * M.c11;
    ctype const c26 = -v.y * M.c18 + v.z * M.c17 + v.w * M.c19 - v.u * M.c22;
    ctype const c27 = v.x * M.c18 - v.z * M.c16 + v.w * M.c20 - v.u * M.c23;
    ctype const c28 = -v.x * M.c17 + v.y * M.c16 + v.w * M.c21 - v.u * M.c24;
    ctype const c29 = -v.x * M.c19 - v.y * M.c20 - v.z * M.c21 - v.u * M.c25;
    ctype const c30 = v.x * M.c22 + v.y * M.c23 + v.z * M.c24 + v.w * M.c25;
    ctype const c31 = v.x * M.c26 + v.y * M.c27 + v.z * M.c28 + v.w * M.c29 + v.u * M.c30;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(MVec3dc<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    ctype const c16 = M.c16 * ctype(s);
    ctype const c17 = M.c17 * ctype(s);
    ctype const c18 = M.c18 * ctype(s);
    ctype const c19 = M.c19 * ctype(s);
    ctype const c20 = M.c20 * ctype(s);
    ctype const c21 = M.c21 * ctype(s);
    ctype const c22 = M.c22 * ctype(s);
    ctype const c23 = M.c23 * ctype(s);
    ctype const c24 = M.c24 * ctype(s);
    ctype const c25 = M.c25 * ctype(s);
    ctype const c26 = M.c26 * ctype(s);
    ctype const c27 = M.c27 * ctype(s);
    ctype const c28 = M.c28 * ctype(s);
    ctype const c29 = M.c29 * ctype(s);
    ctype const c30 = M.c30 * ctype(s);
    ctype const c31 = M.c31 * ctype(s);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    ctype const c16 = ctype(s) * M.c16;
    ctype const c17 = ctype(s) * M.c17;
    ctype const c18 = ctype(s) * M.c18;
    ctype const c19 = ctype(s) * M.c19;
    ctype const c20 = ctype(s) * M.c20;
    ctype const c21 = ctype(s) * M.c21;
    ctype const c22 = ctype(s) * M.c22;
    ctype const c23 = ctype(s) * M.c23;
    ctype const c24 = ctype(s) * M.c24;
    ctype const c25 = ctype(s) * M.c25;
    ctype const c26 = ctype(s) * M.c26;
    ctype const c27 = ctype(s) * M.c27;
    ctype const c28 = ctype(s) * M.c28;
    ctype const c29 = ctype(s) * M.c29;
    ctype const c30 = ctype(s) * M.c30;
    ctype const c31 = ctype(s) * M.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc wdg :: wdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& A,
                                                  MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c7 * B.c0;
    ctype const c8 = A.c0 * B.c8 + A.c8 * B.c0;
    ctype const c9 = A.c0 * B.c9 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c10 + A.c10 * B.c0;
    ctype const c11 = A.c0 * B.c11 + A.c2 * B.c9 - A.c3 * B.c8 + A.c4 * B.c10 -
                      A.c8 * B.c3 + A.c9 * B.c2 + A.c10 * B.c4 + A.c11 * B.c0;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c9 + A.c3 * B.c7 + A.c5 * B.c10 +
                      A.c7 * B.c3 - A.c9 * B.c1 + A.c10 * B.c5 + A.c12 * B.c0;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c8 - A.c2 * B.c7 + A.c6 * B.c10 -
                      A.c7 * B.c2 + A.c8 * B.c1 + A.c10 * B.c6 + A.c13 * B.c0;
    ctype const c14 = A.c0 * B.c14 - A.c4 * B.c7 - A.c5 * B.c8 - A.c6 * B.c9 -
                      A.c7 * B.c4 - A.c8 * B.c5 - A.c9 * B.c6 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 -
                      A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 + A.c15 * B.c0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& A,
                                                  MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c4 + A.c7 * B.c3 - A.c10 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c2 * B.c4 + A.c8 * B.c3 - A.c10 * B.c1;
    ctype const c7 = A.c0 * B.c7 + A.c3 * B.c4 + A.c9 * B.c3 - A.c10 * B.c2;
    ctype const c8 = A.c0 * B.c8 + A.c4 * B.c4 - A.c8 * B.c2 + A.c9 * B.c1;
    ctype const c9 = A.c0 * B.c9 + A.c5 * B.c4 + A.c7 * B.c2 - A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c10 + A.c6 * B.c4 - A.c7 * B.c1 + A.c8 * B.c0;
    ctype const c11 = A.c0 * B.c11 + A.c2 * B.c2 - A.c3 * B.c1 + A.c4 * B.c3;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c2 + A.c3 * B.c0 + A.c5 * B.c3;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c1 - A.c2 * B.c0 + A.c6 * B.c3;
    ctype const c14 = A.c0 * B.c14 - A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                      A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 - A.c7 * B.c11 -
                      A.c8 * B.c12 - A.c9 * B.c13 - A.c10 * B.c14 + A.c11 * B.c0 +
                      A.c12 * B.c1 + A.c13 * B.c2 + A.c14 * B.c3 + A.c15 * B.c4;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& A,
                                                  MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = A.c4 * B.c0;
    ctype const c5 = -A.c0 * B.c10 + A.c3 * B.c7 + A.c4 * B.c1 + A.c5 * B.c0;
    ctype const c6 = -A.c1 * B.c10 + A.c3 * B.c8 + A.c4 * B.c2 + A.c6 * B.c0;
    ctype const c7 = -A.c2 * B.c10 + A.c3 * B.c9 + A.c4 * B.c3 + A.c7 * B.c0;
    ctype const c8 = A.c1 * B.c9 - A.c2 * B.c8 + A.c4 * B.c4 + A.c8 * B.c0;
    ctype const c9 = -A.c0 * B.c9 + A.c2 * B.c7 + A.c4 * B.c5 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c8 - A.c1 * B.c7 + A.c4 * B.c6 + A.c10 * B.c0;
    ctype const c11 = -A.c1 * B.c3 + A.c2 * B.c2 + A.c3 * B.c4 + A.c11 * B.c0;
    ctype const c12 = A.c0 * B.c3 - A.c2 * B.c1 + A.c3 * B.c5 + A.c12 * B.c0;
    ctype const c13 = -A.c0 * B.c2 + A.c1 * B.c1 + A.c3 * B.c6 + A.c13 * B.c0;
    ctype const c14 = -A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                      A.c4 * B.c15 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 -
                      A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 - A.c11 * B.c7 -
                      A.c12 * B.c8 - A.c13 * B.c9 - A.c14 * B.c10 + A.c15 * B.c0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(mv_e,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& M, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = M.c0 * ctype(ps);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(ps,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(PScalar3dc<T> ps, MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = ctype(ps) * M.c0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(mv_e,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& M,
                                                  QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = M.c0 * Q.x;
    ctype const c12 = M.c0 * Q.y;
    ctype const c13 = M.c0 * Q.z;
    ctype const c14 = M.c0 * Q.w;
    ctype const c15 = M.c0 * Q.u;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(quadvec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(QuadVec3dc<T> const& Q,
                                                  MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = Q.x * M.c0;
    ctype const c12 = Q.y * M.c0;
    ctype const c13 = Q.z * M.c0;
    ctype const c14 = Q.w * M.c0;
    ctype const c15 = Q.u * M.c0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& M,
                                                  TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = M.c0 * t.vx;
    ctype const c6 = M.c0 * t.vy;
    ctype const c7 = M.c0 * t.vz;
    ctype const c8 = M.c0 * t.mx;
    ctype const c9 = M.c0 * t.my;
    ctype const c10 = M.c0 * t.mz;
    ctype const c11 = M.c0 * t.px;
    ctype const c12 = M.c0 * t.py;
    ctype const c13 = M.c0 * t.pz;
    ctype const c14 = M.c0 * t.pw;
    ctype const c15 = -M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz - M.c4 * t.vx -
                      M.c5 * t.vy - M.c6 * t.vz - M.c7 * t.px - M.c8 * t.py -
                      M.c9 * t.pz - M.c10 * t.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(TriVec3dc<T> const& t,
                                                  MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = t.vx * M.c0;
    ctype const c6 = t.vy * M.c0;
    ctype const c7 = t.vz * M.c0;
    ctype const c8 = t.mx * M.c0;
    ctype const c9 = t.my * M.c0;
    ctype const c10 = t.mz * M.c0;
    ctype const c11 = t.px * M.c0;
    ctype const c12 = t.py * M.c0;
    ctype const c13 = t.pz * M.c0;
    ctype const c14 = t.pw * M.c0;
    ctype const c15 = -t.vx * M.c4 - t.vy * M.c5 - t.vz * M.c6 - t.mx * M.c1 -
                      t.my * M.c2 - t.mz * M.c3 - t.px * M.c7 - t.py * M.c8 -
                      t.pz * M.c9 - t.pw * M.c10;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& M,
                                                  BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * B.vx;
    ctype const c2 = M.c0 * B.vy;
    ctype const c3 = M.c0 * B.vz;
    ctype const c4 = M.c0 * B.mx;
    ctype const c5 = M.c0 * B.my;
    ctype const c6 = M.c0 * B.mz;
    ctype const c7 = M.c0 * B.px;
    ctype const c8 = M.c0 * B.py;
    ctype const c9 = M.c0 * B.pz;
    ctype const c10 = M.c0 * B.pw;
    ctype const c11 = M.c2 * B.pz - M.c3 * B.py + M.c4 * B.pw - M.c8 * B.vz +
                      M.c9 * B.vy + M.c10 * B.mx;
    ctype const c12 = -M.c1 * B.pz + M.c3 * B.px + M.c5 * B.pw + M.c7 * B.vz -
                      M.c9 * B.vx + M.c10 * B.my;
    ctype const c13 = M.c1 * B.py - M.c2 * B.px + M.c6 * B.pw - M.c7 * B.vy +
                      M.c8 * B.vx + M.c10 * B.mz;
    ctype const c14 = -M.c4 * B.px - M.c5 * B.py - M.c6 * B.pz - M.c7 * B.mx -
                      M.c8 * B.my - M.c9 * B.mz;
    ctype const c15 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz - M.c4 * B.vx -
                      M.c5 * B.vy - M.c6 * B.vz;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(BiVec3dc<T> const& B,
                                                  MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = B.vx * M.c0;
    ctype const c2 = B.vy * M.c0;
    ctype const c3 = B.vz * M.c0;
    ctype const c4 = B.mx * M.c0;
    ctype const c5 = B.my * M.c0;
    ctype const c6 = B.mz * M.c0;
    ctype const c7 = B.px * M.c0;
    ctype const c8 = B.py * M.c0;
    ctype const c9 = B.pz * M.c0;
    ctype const c10 = B.pw * M.c0;
    ctype const c11 = B.vy * M.c9 - B.vz * M.c8 + B.mx * M.c10 - B.py * M.c3 +
                      B.pz * M.c2 + B.pw * M.c4;
    ctype const c12 = -B.vx * M.c9 + B.vz * M.c7 + B.my * M.c10 + B.px * M.c3 -
                      B.pz * M.c1 + B.pw * M.c5;
    ctype const c13 = B.vx * M.c8 - B.vy * M.c7 + B.mz * M.c10 - B.px * M.c2 +
                      B.py * M.c1 + B.pw * M.c6;
    ctype const c14 = -B.mx * M.c7 - B.my * M.c8 - B.mz * M.c9 - B.px * M.c4 -
                      B.py * M.c5 - B.pz * M.c6;
    ctype const c15 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6 - B.mx * M.c1 -
                      B.my * M.c2 - B.mz * M.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& M,
                                                  Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    ctype const c3 = M.c0 * v.w;
    ctype const c4 = M.c0 * v.u;
    ctype const c5 = M.c1 * v.u + M.c7 * v.w - M.c10 * v.x;
    ctype const c6 = M.c2 * v.u + M.c8 * v.w - M.c10 * v.y;
    ctype const c7 = M.c3 * v.u + M.c9 * v.w - M.c10 * v.z;
    ctype const c8 = M.c4 * v.u - M.c8 * v.z + M.c9 * v.y;
    ctype const c9 = M.c5 * v.u + M.c7 * v.z - M.c9 * v.x;
    ctype const c10 = M.c6 * v.u - M.c7 * v.y + M.c8 * v.x;
    ctype const c11 = M.c2 * v.z - M.c3 * v.y + M.c4 * v.w;
    ctype const c12 = -M.c1 * v.z + M.c3 * v.x + M.c5 * v.w;
    ctype const c13 = M.c1 * v.y - M.c2 * v.x + M.c6 * v.w;
    ctype const c14 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z;
    ctype const c15 = M.c11 * v.x + M.c12 * v.y + M.c13 * v.z + M.c14 * v.w + M.c15 * v.u;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v,
                                                  MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    ctype const c3 = v.w * M.c0;
    ctype const c4 = v.u * M.c0;
    ctype const c5 = -v.x * M.c10 + v.w * M.c7 + v.u * M.c1;
    ctype const c6 = -v.y * M.c10 + v.w * M.c8 + v.u * M.c2;
    ctype const c7 = -v.z * M.c10 + v.w * M.c9 + v.u * M.c3;
    ctype const c8 = v.y * M.c9 - v.z * M.c8 + v.u * M.c4;
    ctype const c9 = -v.x * M.c9 + v.z * M.c7 + v.u * M.c5;
    ctype const c10 = v.x * M.c8 - v.y * M.c7 + v.u * M.c6;
    ctype const c11 = -v.y * M.c3 + v.z * M.c2 + v.w * M.c4;
    ctype const c12 = v.x * M.c3 - v.z * M.c1 + v.w * M.c5;
    ctype const c13 = -v.x * M.c2 + v.y * M.c1 + v.w * M.c6;
    ctype const c14 = -v.x * M.c4 - v.y * M.c5 - v.z * M.c6;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14 + v.u * M.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(MVec3dc_E<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& A,
                                                  MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c2 = -A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c3 = -A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c4 = A.c1 * B.c2 - A.c2 * B.c1;
    ctype const c5 = -A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c6 = A.c0 * B.c1 - A.c1 * B.c0;
    ctype const c7 = A.c0 * B.c4 - A.c4 * B.c0;
    ctype const c8 = A.c1 * B.c4 - A.c4 * B.c1;
    ctype const c9 = A.c2 * B.c4 - A.c4 * B.c2;
    ctype const c10 = A.c3 * B.c4 - A.c4 * B.c3;
    ctype const c11 = -A.c1 * B.c7 + A.c2 * B.c6 + A.c3 * B.c8 - A.c4 * B.c11 -
                      A.c6 * B.c2 + A.c7 * B.c1 - A.c8 * B.c3 + A.c11 * B.c4;
    ctype const c12 = A.c0 * B.c7 - A.c2 * B.c5 + A.c3 * B.c9 - A.c4 * B.c12 +
                      A.c5 * B.c2 - A.c7 * B.c0 - A.c9 * B.c3 + A.c12 * B.c4;
    ctype const c13 = -A.c0 * B.c6 + A.c1 * B.c5 + A.c3 * B.c10 - A.c4 * B.c13 -
                      A.c5 * B.c1 + A.c6 * B.c0 - A.c10 * B.c3 + A.c13 * B.c4;
    ctype const c14 = -A.c0 * B.c8 - A.c1 * B.c9 - A.c2 * B.c10 - A.c4 * B.c14 +
                      A.c8 * B.c0 + A.c9 * B.c1 + A.c10 * B.c2 + A.c14 * B.c4;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 -
                      A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c2 - A.c14 * B.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(mv_u,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] MVec3dc_U<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(ps,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] MVec3dc_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(mv_u,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& M,
                                                  QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = M.c0 * Q.x + M.c1 * Q.y + M.c2 * Q.z + M.c3 * Q.w + M.c4 * Q.u;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(quadvec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(QuadVec3dc<T> const& Q,
                                                  MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = Q.x * M.c0 + Q.y * M.c1 + Q.z * M.c2 + Q.w * M.c3 + Q.u * M.c4;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& M,
                                                  TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = -M.c1 * t.vz + M.c2 * t.vy + M.c3 * t.mx - M.c4 * t.px;
    ctype const c12 = M.c0 * t.vz - M.c2 * t.vx + M.c3 * t.my - M.c4 * t.py;
    ctype const c13 = -M.c0 * t.vy + M.c1 * t.vx + M.c3 * t.mz - M.c4 * t.pz;
    ctype const c14 = -M.c0 * t.mx - M.c1 * t.my - M.c2 * t.mz - M.c4 * t.pw;
    ctype const c15 = M.c0 * t.px + M.c1 * t.py + M.c2 * t.pz + M.c3 * t.pw;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(TriVec3dc<T> const& t,
                                                  MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = -t.vy * M.c2 + t.vz * M.c1 - t.mx * M.c3 + t.px * M.c4;
    ctype const c12 = t.vx * M.c2 - t.vz * M.c0 - t.my * M.c3 + t.py * M.c4;
    ctype const c13 = -t.vx * M.c1 + t.vy * M.c0 - t.mz * M.c3 + t.pz * M.c4;
    ctype const c14 = t.mx * M.c0 + t.my * M.c1 + t.mz * M.c2 + t.pw * M.c4;
    ctype const c15 = -t.px * M.c0 - t.py * M.c1 - t.pz * M.c2 - t.pw * M.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& M,
                                                  BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = -M.c0 * B.pw + M.c3 * B.px + M.c4 * B.vx;
    ctype const c6 = -M.c1 * B.pw + M.c3 * B.py + M.c4 * B.vy;
    ctype const c7 = -M.c2 * B.pw + M.c3 * B.pz + M.c4 * B.vz;
    ctype const c8 = M.c1 * B.pz - M.c2 * B.py + M.c4 * B.mx;
    ctype const c9 = -M.c0 * B.pz + M.c2 * B.px + M.c4 * B.my;
    ctype const c10 = M.c0 * B.py - M.c1 * B.px + M.c4 * B.mz;
    ctype const c11 = -M.c1 * B.vz + M.c2 * B.vy + M.c3 * B.mx;
    ctype const c12 = M.c0 * B.vz - M.c2 * B.vx + M.c3 * B.my;
    ctype const c13 = -M.c0 * B.vy + M.c1 * B.vx + M.c3 * B.mz;
    ctype const c14 = -M.c0 * B.mx - M.c1 * B.my - M.c2 * B.mz;
    ctype const c15 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                      M.c9 * B.vy - M.c10 * B.vz - M.c11 * B.px - M.c12 * B.py -
                      M.c13 * B.pz - M.c14 * B.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(BiVec3dc<T> const& B,
                                                  MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = B.vx * M.c4 + B.px * M.c3 - B.pw * M.c0;
    ctype const c6 = B.vy * M.c4 + B.py * M.c3 - B.pw * M.c1;
    ctype const c7 = B.vz * M.c4 + B.pz * M.c3 - B.pw * M.c2;
    ctype const c8 = B.mx * M.c4 - B.py * M.c2 + B.pz * M.c1;
    ctype const c9 = B.my * M.c4 + B.px * M.c2 - B.pz * M.c0;
    ctype const c10 = B.mz * M.c4 - B.px * M.c1 + B.py * M.c0;
    ctype const c11 = B.vy * M.c2 - B.vz * M.c1 + B.mx * M.c3;
    ctype const c12 = -B.vx * M.c2 + B.vz * M.c0 + B.my * M.c3;
    ctype const c13 = B.vx * M.c1 - B.vy * M.c0 + B.mz * M.c3;
    ctype const c14 = -B.mx * M.c0 - B.my * M.c1 - B.mz * M.c2;
    ctype const c15 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                      B.my * M.c6 - B.mz * M.c7 - B.px * M.c11 - B.py * M.c12 -
                      B.pz * M.c13 - B.pw * M.c14;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& M,
                                                  Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -M.c0 * v.w + M.c3 * v.x;
    ctype const c2 = -M.c1 * v.w + M.c3 * v.y;
    ctype const c3 = -M.c2 * v.w + M.c3 * v.z;
    ctype const c4 = M.c1 * v.z - M.c2 * v.y;
    ctype const c5 = -M.c0 * v.z + M.c2 * v.x;
    ctype const c6 = M.c0 * v.y - M.c1 * v.x;
    ctype const c7 = M.c0 * v.u - M.c4 * v.x;
    ctype const c8 = M.c1 * v.u - M.c4 * v.y;
    ctype const c9 = M.c2 * v.u - M.c4 * v.z;
    ctype const c10 = M.c3 * v.u - M.c4 * v.w;
    ctype const c11 = -M.c6 * v.z + M.c7 * v.y - M.c8 * v.w + M.c11 * v.u;
    ctype const c12 = M.c5 * v.z - M.c7 * v.x - M.c9 * v.w + M.c12 * v.u;
    ctype const c13 = -M.c5 * v.y + M.c6 * v.x - M.c10 * v.w + M.c13 * v.u;
    ctype const c14 = M.c8 * v.x + M.c9 * v.y + M.c10 * v.z + M.c14 * v.u;
    ctype const c15 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z - M.c14 * v.w;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v,
                                                  MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -v.x * M.c3 + v.w * M.c0;
    ctype const c2 = -v.y * M.c3 + v.w * M.c1;
    ctype const c3 = -v.z * M.c3 + v.w * M.c2;
    ctype const c4 = v.y * M.c2 - v.z * M.c1;
    ctype const c5 = -v.x * M.c2 + v.z * M.c0;
    ctype const c6 = v.x * M.c1 - v.y * M.c0;
    ctype const c7 = v.x * M.c4 - v.u * M.c0;
    ctype const c8 = v.y * M.c4 - v.u * M.c1;
    ctype const c9 = v.z * M.c4 - v.u * M.c2;
    ctype const c10 = v.w * M.c4 - v.u * M.c3;
    ctype const c11 = -v.y * M.c7 + v.z * M.c6 + v.w * M.c8 - v.u * M.c11;
    ctype const c12 = v.x * M.c7 - v.z * M.c5 + v.w * M.c9 - v.u * M.c12;
    ctype const c13 = -v.x * M.c6 + v.y * M.c5 + v.w * M.c10 - v.u * M.c13;
    ctype const c14 = -v.x * M.c8 - v.y * M.c9 - v.z * M.c10 - v.u * M.c14;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc wdg :: wdg(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(MVec3dc_U<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc wdg :: wdg(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(ps,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(quadvec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] QuadVec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> wdg(PScalar3dc<T> ps, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(ps) * ctype(s));
}

// cga3dc wdg :: wdg(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(s) * ctype(ps));
}

// cga3dc wdg :: wdg(quadvec,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] QuadVec3dc<T> const&,
                                                  [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(quadvec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] QuadVec3dc<T> const&,
                                                  [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(trivec,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dc<T> const&,
                                                  [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(quadvec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] QuadVec3dc<T> const&,
                                                  [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(bivec,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3dc<T> const&,
                                                  [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(quadvec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> wdg(QuadVec3dc<T> const& Q,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(Q.x * v.x + Q.y * v.y + Q.z * v.z + Q.w * v.w + Q.u * v.u);
}

// cga3dc wdg :: wdg(vec,quadvec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(v.x * Q.x + v.y * Q.y + v.z * Q.z + v.w * Q.w + v.u * Q.u);
}

// cga3dc wdg :: wdg(quadvec,s) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> wdg(QuadVec3dc<T> const& Q, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * ctype(s);
    ctype const c1 = Q.y * ctype(s);
    ctype const c2 = Q.z * ctype(s);
    ctype const c3 = Q.w * ctype(s);
    ctype const c4 = Q.u * ctype(s);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc wdg :: wdg(s,quadvec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * Q.x;
    ctype const c1 = ctype(s) * Q.y;
    ctype const c2 = ctype(s) * Q.z;
    ctype const c3 = ctype(s) * Q.w;
    ctype const c4 = ctype(s) * Q.u;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc wdg :: wdg(trivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dc<T> const&,
                                                  [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc wdg :: wdg(trivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> wdg(TriVec3dc<T> const& t,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(-t.vx * B.mx - t.vy * B.my - t.vz * B.mz - t.mx * B.vx -
                             t.my * B.vy - t.mz * B.vz - t.px * B.px - t.py * B.py -
                             t.pz * B.pz - t.pw * B.pw);
}

// cga3dc wdg :: wdg(bivec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> wdg(BiVec3dc<T> const& B,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(-B.vx * t.mx - B.vy * t.my - B.vz * t.mz - B.mx * t.vx -
                             B.my * t.vy - B.mz * t.vz - B.px * t.px - B.py * t.py -
                             B.pz * t.pz - B.pw * t.pw);
}

// cga3dc wdg :: wdg(trivec,vec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> wdg(TriVec3dc<T> const& t,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vy * v.z + t.vz * v.y - t.mx * v.w + t.px * v.u;
    ctype const c1 = t.vx * v.z - t.vz * v.x - t.my * v.w + t.py * v.u;
    ctype const c2 = -t.vx * v.y + t.vy * v.x - t.mz * v.w + t.pz * v.u;
    ctype const c3 = t.mx * v.x + t.my * v.y + t.mz * v.z + t.pw * v.u;
    ctype const c4 = -t.px * v.x - t.py * v.y - t.pz * v.z - t.pw * v.w;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc wdg :: wdg(vec,trivec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * t.vz + v.z * t.vy + v.w * t.mx - v.u * t.px;
    ctype const c1 = v.x * t.vz - v.z * t.vx + v.w * t.my - v.u * t.py;
    ctype const c2 = -v.x * t.vy + v.y * t.vx + v.w * t.mz - v.u * t.pz;
    ctype const c3 = -v.x * t.mx - v.y * t.my - v.z * t.mz - v.u * t.pw;
    ctype const c4 = v.x * t.px + v.y * t.py + v.z * t.pz + v.w * t.pw;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc wdg :: wdg(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> wdg(TriVec3dc<T> const& t, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * ctype(s);
    ctype const c1 = t.vy * ctype(s);
    ctype const c2 = t.vz * ctype(s);
    ctype const c3 = t.mx * ctype(s);
    ctype const c4 = t.my * ctype(s);
    ctype const c5 = t.mz * ctype(s);
    ctype const c6 = t.px * ctype(s);
    ctype const c7 = t.py * ctype(s);
    ctype const c8 = t.pz * ctype(s);
    ctype const c9 = t.pw * ctype(s);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc wdg :: wdg(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.vx;
    ctype const c1 = ctype(s) * t.vy;
    ctype const c2 = ctype(s) * t.vz;
    ctype const c3 = ctype(s) * t.mx;
    ctype const c4 = ctype(s) * t.my;
    ctype const c5 = ctype(s) * t.mz;
    ctype const c6 = ctype(s) * t.px;
    ctype const c7 = ctype(s) * t.py;
    ctype const c8 = ctype(s) * t.pz;
    ctype const c9 = ctype(s) * t.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc wdg :: wdg(bivec,bivec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> wdg(BiVec3dc<T> const& B1,
                                                   BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vy * B2.pz - B1.vz * B2.py + B1.mx * B2.pw - B1.py * B2.vz +
                     B1.pz * B2.vy + B1.pw * B2.mx;
    ctype const c1 = -B1.vx * B2.pz + B1.vz * B2.px + B1.my * B2.pw + B1.px * B2.vz -
                     B1.pz * B2.vx + B1.pw * B2.my;
    ctype const c2 = B1.vx * B2.py - B1.vy * B2.px + B1.mz * B2.pw - B1.px * B2.vy +
                     B1.py * B2.vx + B1.pw * B2.mz;
    ctype const c3 = -B1.mx * B2.px - B1.my * B2.py - B1.mz * B2.pz - B1.px * B2.mx -
                     B1.py * B2.my - B1.pz * B2.mz;
    ctype const c4 = -B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz - B1.mx * B2.vx -
                     B1.my * B2.vy - B1.mz * B2.vz;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc wdg :: wdg(bivec,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> wdg(BiVec3dc<T> const& B,
                                                  Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.u + B.px * v.w - B.pw * v.x;
    ctype const c1 = B.vy * v.u + B.py * v.w - B.pw * v.y;
    ctype const c2 = B.vz * v.u + B.pz * v.w - B.pw * v.z;
    ctype const c3 = B.mx * v.u - B.py * v.z + B.pz * v.y;
    ctype const c4 = B.my * v.u + B.px * v.z - B.pz * v.x;
    ctype const c5 = B.mz * v.u - B.px * v.y + B.py * v.x;
    ctype const c6 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c7 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c8 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c9 = -B.mx * v.x - B.my * v.y - B.mz * v.z;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc wdg :: wdg(vec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v,
                                                  BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * B.pw + v.w * B.px + v.u * B.vx;
    ctype const c1 = -v.y * B.pw + v.w * B.py + v.u * B.vy;
    ctype const c2 = -v.z * B.pw + v.w * B.pz + v.u * B.vz;
    ctype const c3 = v.y * B.pz - v.z * B.py + v.u * B.mx;
    ctype const c4 = -v.x * B.pz + v.z * B.px + v.u * B.my;
    ctype const c5 = v.x * B.py - v.y * B.px + v.u * B.mz;
    ctype const c6 = -v.y * B.vz + v.z * B.vy + v.w * B.mx;
    ctype const c7 = v.x * B.vz - v.z * B.vx + v.w * B.my;
    ctype const c8 = -v.x * B.vy + v.y * B.vx + v.w * B.mz;
    ctype const c9 = -v.x * B.mx - v.y * B.my - v.z * B.mz;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc wdg :: wdg(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> wdg(BiVec3dc<T> const& B, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    ctype const c6 = B.px * ctype(s);
    ctype const c7 = B.py * ctype(s);
    ctype const c8 = B.pz * ctype(s);
    ctype const c9 = B.pw * ctype(s);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc wdg :: wdg(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    ctype const c6 = ctype(s) * B.px;
    ctype const c7 = ctype(s) * B.py;
    ctype const c8 = ctype(s) * B.pz;
    ctype const c9 = ctype(s) * B.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc wdg :: wdg(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v1, Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c1 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c2 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c3 = v1.y * v2.z - v1.z * v2.y;
    ctype const c4 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c5 = v1.x * v2.y - v1.y * v2.x;
    ctype const c6 = v1.x * v2.u - v1.u * v2.x;
    ctype const c7 = v1.y * v2.u - v1.u * v2.y;
    ctype const c8 = v1.z * v2.u - v1.u * v2.z;
    ctype const c9 = v1.w * v2.u - v1.u * v2.w;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc wdg :: wdg(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> wdg(Vec3dc<T> const& v, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    ctype const c4 = v.u * ctype(s);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc wdg :: wdg(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> wdg(Scalar3dc<T> s, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    ctype const c4 = ctype(s) * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc wdg :: wdg(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> wdg(Scalar3dc<T> s1, Scalar3dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product for complements)
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product, just depending on the outer product (wdg)
// as well as the complement and thus the pseudoscalar of the space
// (this definition does NOT connect directly to the geometric product,
// but to the outer product exclusively)
//
//     rwdg(ul, ur) = cmpl(wdg(cmpl(ul),cmpl(ur))) = cmpl(cmpl(ul) ^ cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// cga3dc rwdg :: rwdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& A, MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c31 + A.c1 * B.c26 + A.c2 * B.c27 + A.c3 * B.c28 +
                     A.c4 * B.c29 + A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                     A.c8 * B.c21 - A.c9 * B.c16 - A.c10 * B.c17 - A.c11 * B.c18 -
                     A.c12 * B.c22 - A.c13 * B.c23 - A.c14 * B.c24 - A.c15 * B.c25 -
                     A.c16 * B.c9 - A.c17 * B.c10 - A.c18 * B.c11 - A.c19 * B.c6 -
                     A.c20 * B.c7 - A.c21 * B.c8 - A.c22 * B.c12 - A.c23 * B.c13 -
                     A.c24 * B.c14 - A.c25 * B.c15 + A.c26 * B.c1 + A.c27 * B.c2 +
                     A.c28 * B.c3 + A.c29 * B.c4 + A.c30 * B.c5 + A.c31 * B.c0;
    ctype const c1 = A.c1 * B.c31 + A.c6 * B.c29 + A.c10 * B.c28 - A.c11 * B.c27 -
                     A.c12 * B.c30 + A.c16 * B.c25 + A.c20 * B.c24 - A.c21 * B.c23 -
                     A.c23 * B.c21 + A.c24 * B.c20 + A.c25 * B.c16 + A.c27 * B.c11 -
                     A.c28 * B.c10 - A.c29 * B.c6 + A.c30 * B.c12 + A.c31 * B.c1;
    ctype const c2 = A.c2 * B.c31 + A.c7 * B.c29 - A.c9 * B.c28 + A.c11 * B.c26 -
                     A.c13 * B.c30 + A.c17 * B.c25 - A.c19 * B.c24 + A.c21 * B.c22 +
                     A.c22 * B.c21 - A.c24 * B.c19 + A.c25 * B.c17 - A.c26 * B.c11 +
                     A.c28 * B.c9 - A.c29 * B.c7 + A.c30 * B.c13 + A.c31 * B.c2;
    ctype const c3 = A.c3 * B.c31 + A.c8 * B.c29 + A.c9 * B.c27 - A.c10 * B.c26 -
                     A.c14 * B.c30 + A.c18 * B.c25 + A.c19 * B.c23 - A.c20 * B.c22 -
                     A.c22 * B.c20 + A.c23 * B.c19 + A.c25 * B.c18 + A.c26 * B.c10 -
                     A.c27 * B.c9 - A.c29 * B.c8 + A.c30 * B.c14 + A.c31 * B.c3;
    ctype const c4 = A.c4 * B.c31 - A.c6 * B.c26 - A.c7 * B.c27 - A.c8 * B.c28 -
                     A.c15 * B.c30 - A.c16 * B.c22 - A.c17 * B.c23 - A.c18 * B.c24 -
                     A.c22 * B.c16 - A.c23 * B.c17 - A.c24 * B.c18 + A.c26 * B.c6 +
                     A.c27 * B.c7 + A.c28 * B.c8 + A.c30 * B.c15 + A.c31 * B.c4;
    ctype const c5 = A.c5 * B.c31 + A.c12 * B.c26 + A.c13 * B.c27 + A.c14 * B.c28 +
                     A.c15 * B.c29 - A.c16 * B.c19 - A.c17 * B.c20 - A.c18 * B.c21 -
                     A.c19 * B.c16 - A.c20 * B.c17 - A.c21 * B.c18 - A.c26 * B.c12 -
                     A.c27 * B.c13 - A.c28 * B.c14 - A.c29 * B.c15 + A.c31 * B.c5;
    ctype const c6 = A.c6 * B.c31 + A.c16 * B.c30 - A.c23 * B.c28 + A.c24 * B.c27 +
                     A.c27 * B.c24 - A.c28 * B.c23 + A.c30 * B.c16 + A.c31 * B.c6;
    ctype const c7 = A.c7 * B.c31 + A.c17 * B.c30 + A.c22 * B.c28 - A.c24 * B.c26 -
                     A.c26 * B.c24 + A.c28 * B.c22 + A.c30 * B.c17 + A.c31 * B.c7;
    ctype const c8 = A.c8 * B.c31 + A.c18 * B.c30 - A.c22 * B.c27 + A.c23 * B.c26 +
                     A.c26 * B.c23 - A.c27 * B.c22 + A.c30 * B.c18 + A.c31 * B.c8;
    ctype const c9 = A.c9 * B.c31 + A.c19 * B.c30 + A.c22 * B.c29 - A.c25 * B.c26 -
                     A.c26 * B.c25 + A.c29 * B.c22 + A.c30 * B.c19 + A.c31 * B.c9;
    ctype const c10 = A.c10 * B.c31 + A.c20 * B.c30 + A.c23 * B.c29 - A.c25 * B.c27 -
                      A.c27 * B.c25 + A.c29 * B.c23 + A.c30 * B.c20 + A.c31 * B.c10;
    ctype const c11 = A.c11 * B.c31 + A.c21 * B.c30 + A.c24 * B.c29 - A.c25 * B.c28 -
                      A.c28 * B.c25 + A.c29 * B.c24 + A.c30 * B.c21 + A.c31 * B.c11;
    ctype const c12 = A.c12 * B.c31 + A.c16 * B.c29 + A.c20 * B.c28 - A.c21 * B.c27 -
                      A.c27 * B.c21 + A.c28 * B.c20 + A.c29 * B.c16 + A.c31 * B.c12;
    ctype const c13 = A.c13 * B.c31 + A.c17 * B.c29 - A.c19 * B.c28 + A.c21 * B.c26 +
                      A.c26 * B.c21 - A.c28 * B.c19 + A.c29 * B.c17 + A.c31 * B.c13;
    ctype const c14 = A.c14 * B.c31 + A.c18 * B.c29 + A.c19 * B.c27 - A.c20 * B.c26 -
                      A.c26 * B.c20 + A.c27 * B.c19 + A.c29 * B.c18 + A.c31 * B.c14;
    ctype const c15 = A.c15 * B.c31 - A.c16 * B.c26 - A.c17 * B.c27 - A.c18 * B.c28 -
                      A.c26 * B.c16 - A.c27 * B.c17 - A.c28 * B.c18 + A.c31 * B.c15;
    ctype const c16 = A.c16 * B.c31 - A.c27 * B.c28 + A.c28 * B.c27 + A.c31 * B.c16;
    ctype const c17 = A.c17 * B.c31 + A.c26 * B.c28 - A.c28 * B.c26 + A.c31 * B.c17;
    ctype const c18 = A.c18 * B.c31 - A.c26 * B.c27 + A.c27 * B.c26 + A.c31 * B.c18;
    ctype const c19 = A.c19 * B.c31 + A.c26 * B.c29 - A.c29 * B.c26 + A.c31 * B.c19;
    ctype const c20 = A.c20 * B.c31 + A.c27 * B.c29 - A.c29 * B.c27 + A.c31 * B.c20;
    ctype const c21 = A.c21 * B.c31 + A.c28 * B.c29 - A.c29 * B.c28 + A.c31 * B.c21;
    ctype const c22 = A.c22 * B.c31 - A.c26 * B.c30 + A.c30 * B.c26 + A.c31 * B.c22;
    ctype const c23 = A.c23 * B.c31 - A.c27 * B.c30 + A.c30 * B.c27 + A.c31 * B.c23;
    ctype const c24 = A.c24 * B.c31 - A.c28 * B.c30 + A.c30 * B.c28 + A.c31 * B.c24;
    ctype const c25 = A.c25 * B.c31 - A.c29 * B.c30 + A.c30 * B.c29 + A.c31 * B.c25;
    ctype const c26 = A.c26 * B.c31 + A.c31 * B.c26;
    ctype const c27 = A.c27 * B.c31 + A.c31 * B.c27;
    ctype const c28 = A.c28 * B.c31 + A.c31 * B.c28;
    ctype const c29 = A.c29 * B.c31 + A.c31 * B.c29;
    ctype const c30 = A.c30 * B.c31 + A.c31 * B.c30;
    ctype const c31 = A.c31 * B.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& A,
                                                 MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 +
                     A.c5 * B.c15 - A.c16 * B.c4 - A.c17 * B.c5 - A.c18 * B.c6 -
                     A.c19 * B.c1 - A.c20 * B.c2 - A.c21 * B.c3 - A.c22 * B.c7 -
                     A.c23 * B.c8 - A.c24 * B.c9 - A.c25 * B.c10 + A.c31 * B.c0;
    ctype const c1 = A.c6 * B.c14 + A.c10 * B.c13 - A.c11 * B.c12 - A.c12 * B.c15 +
                     A.c27 * B.c6 - A.c28 * B.c5 - A.c29 * B.c1 + A.c30 * B.c7;
    ctype const c2 = A.c7 * B.c14 - A.c9 * B.c13 + A.c11 * B.c11 - A.c13 * B.c15 -
                     A.c26 * B.c6 + A.c28 * B.c4 - A.c29 * B.c2 + A.c30 * B.c8;
    ctype const c3 = A.c8 * B.c14 + A.c9 * B.c12 - A.c10 * B.c11 - A.c14 * B.c15 +
                     A.c26 * B.c5 - A.c27 * B.c4 - A.c29 * B.c3 + A.c30 * B.c9;
    ctype const c4 = -A.c6 * B.c11 - A.c7 * B.c12 - A.c8 * B.c13 - A.c15 * B.c15 +
                     A.c26 * B.c1 + A.c27 * B.c2 + A.c28 * B.c3 + A.c30 * B.c10;
    ctype const c5 = A.c12 * B.c11 + A.c13 * B.c12 + A.c14 * B.c13 + A.c15 * B.c14 -
                     A.c26 * B.c7 - A.c27 * B.c8 - A.c28 * B.c9 - A.c29 * B.c10;
    ctype const c6 = A.c16 * B.c15 - A.c23 * B.c13 + A.c24 * B.c12 + A.c31 * B.c1;
    ctype const c7 = A.c17 * B.c15 + A.c22 * B.c13 - A.c24 * B.c11 + A.c31 * B.c2;
    ctype const c8 = A.c18 * B.c15 - A.c22 * B.c12 + A.c23 * B.c11 + A.c31 * B.c3;
    ctype const c9 = A.c19 * B.c15 + A.c22 * B.c14 - A.c25 * B.c11 + A.c31 * B.c4;
    ctype const c10 = A.c20 * B.c15 + A.c23 * B.c14 - A.c25 * B.c12 + A.c31 * B.c5;
    ctype const c11 = A.c21 * B.c15 + A.c24 * B.c14 - A.c25 * B.c13 + A.c31 * B.c6;
    ctype const c12 = A.c16 * B.c14 + A.c20 * B.c13 - A.c21 * B.c12 + A.c31 * B.c7;
    ctype const c13 = A.c17 * B.c14 - A.c19 * B.c13 + A.c21 * B.c11 + A.c31 * B.c8;
    ctype const c14 = A.c18 * B.c14 + A.c19 * B.c12 - A.c20 * B.c11 + A.c31 * B.c9;
    ctype const c15 = -A.c16 * B.c11 - A.c17 * B.c12 - A.c18 * B.c13 + A.c31 * B.c10;
    ctype const c16 = -A.c27 * B.c13 + A.c28 * B.c12;
    ctype const c17 = A.c26 * B.c13 - A.c28 * B.c11;
    ctype const c18 = -A.c26 * B.c12 + A.c27 * B.c11;
    ctype const c19 = A.c26 * B.c14 - A.c29 * B.c11;
    ctype const c20 = A.c27 * B.c14 - A.c29 * B.c12;
    ctype const c21 = A.c28 * B.c14 - A.c29 * B.c13;
    ctype const c22 = -A.c26 * B.c15 + A.c30 * B.c11;
    ctype const c23 = -A.c27 * B.c15 + A.c30 * B.c12;
    ctype const c24 = -A.c28 * B.c15 + A.c30 * B.c13;
    ctype const c25 = -A.c29 * B.c15 + A.c30 * B.c14;
    ctype const c26 = A.c31 * B.c11;
    ctype const c27 = A.c31 * B.c12;
    ctype const c28 = A.c31 * B.c13;
    ctype const c29 = A.c31 * B.c14;
    ctype const c30 = A.c31 * B.c15;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& A,
                                                 MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c31 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                     A.c4 * B.c16 - A.c5 * B.c17 - A.c6 * B.c18 - A.c7 * B.c22 -
                     A.c8 * B.c23 - A.c9 * B.c24 - A.c10 * B.c25 + A.c11 * B.c1 +
                     A.c12 * B.c2 + A.c13 * B.c3 + A.c14 * B.c4 + A.c15 * B.c5;
    ctype const c1 = A.c1 * B.c29 + A.c5 * B.c28 - A.c6 * B.c27 - A.c7 * B.c30 +
                     A.c12 * B.c11 - A.c13 * B.c10 - A.c14 * B.c6 + A.c15 * B.c12;
    ctype const c2 = A.c2 * B.c29 - A.c4 * B.c28 + A.c6 * B.c26 - A.c8 * B.c30 -
                     A.c11 * B.c11 + A.c13 * B.c9 - A.c14 * B.c7 + A.c15 * B.c13;
    ctype const c3 = A.c3 * B.c29 + A.c4 * B.c27 - A.c5 * B.c26 - A.c9 * B.c30 +
                     A.c11 * B.c10 - A.c12 * B.c9 - A.c14 * B.c8 + A.c15 * B.c14;
    ctype const c4 = -A.c1 * B.c26 - A.c2 * B.c27 - A.c3 * B.c28 - A.c10 * B.c30 +
                     A.c11 * B.c6 + A.c12 * B.c7 + A.c13 * B.c8 + A.c15 * B.c15;
    ctype const c5 = A.c7 * B.c26 + A.c8 * B.c27 + A.c9 * B.c28 + A.c10 * B.c29 -
                     A.c11 * B.c12 - A.c12 * B.c13 - A.c13 * B.c14 - A.c14 * B.c15;
    ctype const c6 = A.c1 * B.c31 + A.c12 * B.c24 - A.c13 * B.c23 + A.c15 * B.c16;
    ctype const c7 = A.c2 * B.c31 - A.c11 * B.c24 + A.c13 * B.c22 + A.c15 * B.c17;
    ctype const c8 = A.c3 * B.c31 + A.c11 * B.c23 - A.c12 * B.c22 + A.c15 * B.c18;
    ctype const c9 = A.c4 * B.c31 - A.c11 * B.c25 + A.c14 * B.c22 + A.c15 * B.c19;
    ctype const c10 = A.c5 * B.c31 - A.c12 * B.c25 + A.c14 * B.c23 + A.c15 * B.c20;
    ctype const c11 = A.c6 * B.c31 - A.c13 * B.c25 + A.c14 * B.c24 + A.c15 * B.c21;
    ctype const c12 = A.c7 * B.c31 - A.c12 * B.c21 + A.c13 * B.c20 + A.c14 * B.c16;
    ctype const c13 = A.c8 * B.c31 + A.c11 * B.c21 - A.c13 * B.c19 + A.c14 * B.c17;
    ctype const c14 = A.c9 * B.c31 - A.c11 * B.c20 + A.c12 * B.c19 + A.c14 * B.c18;
    ctype const c15 = A.c10 * B.c31 - A.c11 * B.c16 - A.c12 * B.c17 - A.c13 * B.c18;
    ctype const c16 = -A.c12 * B.c28 + A.c13 * B.c27;
    ctype const c17 = A.c11 * B.c28 - A.c13 * B.c26;
    ctype const c18 = -A.c11 * B.c27 + A.c12 * B.c26;
    ctype const c19 = A.c11 * B.c29 - A.c14 * B.c26;
    ctype const c20 = A.c12 * B.c29 - A.c14 * B.c27;
    ctype const c21 = A.c13 * B.c29 - A.c14 * B.c28;
    ctype const c22 = -A.c11 * B.c30 + A.c15 * B.c26;
    ctype const c23 = -A.c12 * B.c30 + A.c15 * B.c27;
    ctype const c24 = -A.c13 * B.c30 + A.c15 * B.c28;
    ctype const c25 = -A.c14 * B.c30 + A.c15 * B.c29;
    ctype const c26 = A.c11 * B.c31;
    ctype const c27 = A.c12 * B.c31;
    ctype const c28 = A.c13 * B.c31;
    ctype const c29 = A.c14 * B.c31;
    ctype const c30 = A.c15 * B.c31;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& A,
                                                 MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c6 * B.c8 - A.c7 * B.c9 - A.c8 * B.c10 -
                     A.c9 * B.c5 - A.c10 * B.c6 - A.c11 * B.c7 - A.c12 * B.c11 -
                     A.c13 * B.c12 - A.c14 * B.c13 - A.c15 * B.c14 + A.c26 * B.c0 +
                     A.c27 * B.c1 + A.c28 * B.c2 + A.c29 * B.c3 + A.c30 * B.c4;
    ctype const c1 = A.c1 * B.c15 + A.c16 * B.c14 + A.c20 * B.c13 - A.c21 * B.c12 -
                     A.c23 * B.c10 + A.c24 * B.c9 + A.c25 * B.c5 + A.c31 * B.c0;
    ctype const c2 = A.c2 * B.c15 + A.c17 * B.c14 - A.c19 * B.c13 + A.c21 * B.c11 +
                     A.c22 * B.c10 - A.c24 * B.c8 + A.c25 * B.c6 + A.c31 * B.c1;
    ctype const c3 = A.c3 * B.c15 + A.c18 * B.c14 + A.c19 * B.c12 - A.c20 * B.c11 -
                     A.c22 * B.c9 + A.c23 * B.c8 + A.c25 * B.c7 + A.c31 * B.c2;
    ctype const c4 = A.c4 * B.c15 - A.c16 * B.c11 - A.c17 * B.c12 - A.c18 * B.c13 -
                     A.c22 * B.c5 - A.c23 * B.c6 - A.c24 * B.c7 + A.c31 * B.c3;
    ctype const c5 = A.c5 * B.c15 - A.c16 * B.c8 - A.c17 * B.c9 - A.c18 * B.c10 -
                     A.c19 * B.c5 - A.c20 * B.c6 - A.c21 * B.c7 + A.c31 * B.c4;
    ctype const c6 = A.c6 * B.c15 + A.c27 * B.c13 - A.c28 * B.c12 + A.c30 * B.c5;
    ctype const c7 = A.c7 * B.c15 - A.c26 * B.c13 + A.c28 * B.c11 + A.c30 * B.c6;
    ctype const c8 = A.c8 * B.c15 + A.c26 * B.c12 - A.c27 * B.c11 + A.c30 * B.c7;
    ctype const c9 = A.c9 * B.c15 - A.c26 * B.c14 + A.c29 * B.c11 + A.c30 * B.c8;
    ctype const c10 = A.c10 * B.c15 - A.c27 * B.c14 + A.c29 * B.c12 + A.c30 * B.c9;
    ctype const c11 = A.c11 * B.c15 - A.c28 * B.c14 + A.c29 * B.c13 + A.c30 * B.c10;
    ctype const c12 = A.c12 * B.c15 - A.c27 * B.c10 + A.c28 * B.c9 + A.c29 * B.c5;
    ctype const c13 = A.c13 * B.c15 + A.c26 * B.c10 - A.c28 * B.c8 + A.c29 * B.c6;
    ctype const c14 = A.c14 * B.c15 - A.c26 * B.c9 + A.c27 * B.c8 + A.c29 * B.c7;
    ctype const c15 = A.c15 * B.c15 - A.c26 * B.c5 - A.c27 * B.c6 - A.c28 * B.c7;
    ctype const c16 = A.c16 * B.c15 + A.c31 * B.c5;
    ctype const c17 = A.c17 * B.c15 + A.c31 * B.c6;
    ctype const c18 = A.c18 * B.c15 + A.c31 * B.c7;
    ctype const c19 = A.c19 * B.c15 + A.c31 * B.c8;
    ctype const c20 = A.c20 * B.c15 + A.c31 * B.c9;
    ctype const c21 = A.c21 * B.c15 + A.c31 * B.c10;
    ctype const c22 = A.c22 * B.c15 + A.c31 * B.c11;
    ctype const c23 = A.c23 * B.c15 + A.c31 * B.c12;
    ctype const c24 = A.c24 * B.c15 + A.c31 * B.c13;
    ctype const c25 = A.c25 * B.c15 + A.c31 * B.c14;
    ctype const c26 = A.c26 * B.c15;
    ctype const c27 = A.c27 * B.c15;
    ctype const c28 = A.c28 * B.c15;
    ctype const c29 = A.c29 * B.c15;
    ctype const c30 = A.c30 * B.c15;
    ctype const c31 = A.c31 * B.c15;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& A,
                                                 MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c26 + A.c1 * B.c27 + A.c2 * B.c28 + A.c3 * B.c29 +
                     A.c4 * B.c30 - A.c5 * B.c9 - A.c6 * B.c10 - A.c7 * B.c11 -
                     A.c8 * B.c6 - A.c9 * B.c7 - A.c10 * B.c8 - A.c11 * B.c12 -
                     A.c12 * B.c13 - A.c13 * B.c14 - A.c14 * B.c15 + A.c15 * B.c0;
    ctype const c1 = A.c0 * B.c31 + A.c5 * B.c25 + A.c9 * B.c24 - A.c10 * B.c23 -
                     A.c12 * B.c21 + A.c13 * B.c20 + A.c14 * B.c16 + A.c15 * B.c1;
    ctype const c2 = A.c1 * B.c31 + A.c6 * B.c25 - A.c8 * B.c24 + A.c10 * B.c22 +
                     A.c11 * B.c21 - A.c13 * B.c19 + A.c14 * B.c17 + A.c15 * B.c2;
    ctype const c3 = A.c2 * B.c31 + A.c7 * B.c25 + A.c8 * B.c23 - A.c9 * B.c22 -
                     A.c11 * B.c20 + A.c12 * B.c19 + A.c14 * B.c18 + A.c15 * B.c3;
    ctype const c4 = A.c3 * B.c31 - A.c5 * B.c22 - A.c6 * B.c23 - A.c7 * B.c24 -
                     A.c11 * B.c16 - A.c12 * B.c17 - A.c13 * B.c18 + A.c15 * B.c4;
    ctype const c5 = A.c4 * B.c31 - A.c5 * B.c19 - A.c6 * B.c20 - A.c7 * B.c21 -
                     A.c8 * B.c16 - A.c9 * B.c17 - A.c10 * B.c18 + A.c15 * B.c5;
    ctype const c6 = A.c5 * B.c30 - A.c12 * B.c28 + A.c13 * B.c27 + A.c15 * B.c6;
    ctype const c7 = A.c6 * B.c30 + A.c11 * B.c28 - A.c13 * B.c26 + A.c15 * B.c7;
    ctype const c8 = A.c7 * B.c30 - A.c11 * B.c27 + A.c12 * B.c26 + A.c15 * B.c8;
    ctype const c9 = A.c8 * B.c30 + A.c11 * B.c29 - A.c14 * B.c26 + A.c15 * B.c9;
    ctype const c10 = A.c9 * B.c30 + A.c12 * B.c29 - A.c14 * B.c27 + A.c15 * B.c10;
    ctype const c11 = A.c10 * B.c30 + A.c13 * B.c29 - A.c14 * B.c28 + A.c15 * B.c11;
    ctype const c12 = A.c5 * B.c29 + A.c9 * B.c28 - A.c10 * B.c27 + A.c15 * B.c12;
    ctype const c13 = A.c6 * B.c29 - A.c8 * B.c28 + A.c10 * B.c26 + A.c15 * B.c13;
    ctype const c14 = A.c7 * B.c29 + A.c8 * B.c27 - A.c9 * B.c26 + A.c15 * B.c14;
    ctype const c15 = -A.c5 * B.c26 - A.c6 * B.c27 - A.c7 * B.c28 + A.c15 * B.c15;
    ctype const c16 = A.c5 * B.c31 + A.c15 * B.c16;
    ctype const c17 = A.c6 * B.c31 + A.c15 * B.c17;
    ctype const c18 = A.c7 * B.c31 + A.c15 * B.c18;
    ctype const c19 = A.c8 * B.c31 + A.c15 * B.c19;
    ctype const c20 = A.c9 * B.c31 + A.c15 * B.c20;
    ctype const c21 = A.c10 * B.c31 + A.c15 * B.c21;
    ctype const c22 = A.c11 * B.c31 + A.c15 * B.c22;
    ctype const c23 = A.c12 * B.c31 + A.c15 * B.c23;
    ctype const c24 = A.c13 * B.c31 + A.c15 * B.c24;
    ctype const c25 = A.c14 * B.c31 + A.c15 * B.c25;
    ctype const c26 = A.c15 * B.c26;
    ctype const c27 = A.c15 * B.c27;
    ctype const c28 = A.c15 * B.c28;
    ctype const c29 = A.c15 * B.c29;
    ctype const c30 = A.c15 * B.c30;
    ctype const c31 = A.c15 * B.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& M, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    ctype const c4 = M.c4 * ctype(ps);
    ctype const c5 = M.c5 * ctype(ps);
    ctype const c6 = M.c6 * ctype(ps);
    ctype const c7 = M.c7 * ctype(ps);
    ctype const c8 = M.c8 * ctype(ps);
    ctype const c9 = M.c9 * ctype(ps);
    ctype const c10 = M.c10 * ctype(ps);
    ctype const c11 = M.c11 * ctype(ps);
    ctype const c12 = M.c12 * ctype(ps);
    ctype const c13 = M.c13 * ctype(ps);
    ctype const c14 = M.c14 * ctype(ps);
    ctype const c15 = M.c15 * ctype(ps);
    ctype const c16 = M.c16 * ctype(ps);
    ctype const c17 = M.c17 * ctype(ps);
    ctype const c18 = M.c18 * ctype(ps);
    ctype const c19 = M.c19 * ctype(ps);
    ctype const c20 = M.c20 * ctype(ps);
    ctype const c21 = M.c21 * ctype(ps);
    ctype const c22 = M.c22 * ctype(ps);
    ctype const c23 = M.c23 * ctype(ps);
    ctype const c24 = M.c24 * ctype(ps);
    ctype const c25 = M.c25 * ctype(ps);
    ctype const c26 = M.c26 * ctype(ps);
    ctype const c27 = M.c27 * ctype(ps);
    ctype const c28 = M.c28 * ctype(ps);
    ctype const c29 = M.c29 * ctype(ps);
    ctype const c30 = M.c30 * ctype(ps);
    ctype const c31 = M.c31 * ctype(ps);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    ctype const c4 = ctype(ps) * M.c4;
    ctype const c5 = ctype(ps) * M.c5;
    ctype const c6 = ctype(ps) * M.c6;
    ctype const c7 = ctype(ps) * M.c7;
    ctype const c8 = ctype(ps) * M.c8;
    ctype const c9 = ctype(ps) * M.c9;
    ctype const c10 = ctype(ps) * M.c10;
    ctype const c11 = ctype(ps) * M.c11;
    ctype const c12 = ctype(ps) * M.c12;
    ctype const c13 = ctype(ps) * M.c13;
    ctype const c14 = ctype(ps) * M.c14;
    ctype const c15 = ctype(ps) * M.c15;
    ctype const c16 = ctype(ps) * M.c16;
    ctype const c17 = ctype(ps) * M.c17;
    ctype const c18 = ctype(ps) * M.c18;
    ctype const c19 = ctype(ps) * M.c19;
    ctype const c20 = ctype(ps) * M.c20;
    ctype const c21 = ctype(ps) * M.c21;
    ctype const c22 = ctype(ps) * M.c22;
    ctype const c23 = ctype(ps) * M.c23;
    ctype const c24 = ctype(ps) * M.c24;
    ctype const c25 = ctype(ps) * M.c25;
    ctype const c26 = ctype(ps) * M.c26;
    ctype const c27 = ctype(ps) * M.c27;
    ctype const c28 = ctype(ps) * M.c28;
    ctype const c29 = ctype(ps) * M.c29;
    ctype const c30 = ctype(ps) * M.c30;
    ctype const c31 = ctype(ps) * M.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,quadvec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& M,
                                                 QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * Q.x + M.c2 * Q.y + M.c3 * Q.z + M.c4 * Q.w + M.c5 * Q.u;
    ctype const c1 = M.c6 * Q.w + M.c10 * Q.z - M.c11 * Q.y - M.c12 * Q.u;
    ctype const c2 = M.c7 * Q.w - M.c9 * Q.z + M.c11 * Q.x - M.c13 * Q.u;
    ctype const c3 = M.c8 * Q.w + M.c9 * Q.y - M.c10 * Q.x - M.c14 * Q.u;
    ctype const c4 = -M.c6 * Q.x - M.c7 * Q.y - M.c8 * Q.z - M.c15 * Q.u;
    ctype const c5 = M.c12 * Q.x + M.c13 * Q.y + M.c14 * Q.z + M.c15 * Q.w;
    ctype const c6 = M.c16 * Q.u - M.c23 * Q.z + M.c24 * Q.y;
    ctype const c7 = M.c17 * Q.u + M.c22 * Q.z - M.c24 * Q.x;
    ctype const c8 = M.c18 * Q.u - M.c22 * Q.y + M.c23 * Q.x;
    ctype const c9 = M.c19 * Q.u + M.c22 * Q.w - M.c25 * Q.x;
    ctype const c10 = M.c20 * Q.u + M.c23 * Q.w - M.c25 * Q.y;
    ctype const c11 = M.c21 * Q.u + M.c24 * Q.w - M.c25 * Q.z;
    ctype const c12 = M.c16 * Q.w + M.c20 * Q.z - M.c21 * Q.y;
    ctype const c13 = M.c17 * Q.w - M.c19 * Q.z + M.c21 * Q.x;
    ctype const c14 = M.c18 * Q.w + M.c19 * Q.y - M.c20 * Q.x;
    ctype const c15 = -M.c16 * Q.x - M.c17 * Q.y - M.c18 * Q.z;
    ctype const c16 = -M.c27 * Q.z + M.c28 * Q.y;
    ctype const c17 = M.c26 * Q.z - M.c28 * Q.x;
    ctype const c18 = -M.c26 * Q.y + M.c27 * Q.x;
    ctype const c19 = M.c26 * Q.w - M.c29 * Q.x;
    ctype const c20 = M.c27 * Q.w - M.c29 * Q.y;
    ctype const c21 = M.c28 * Q.w - M.c29 * Q.z;
    ctype const c22 = -M.c26 * Q.u + M.c30 * Q.x;
    ctype const c23 = -M.c27 * Q.u + M.c30 * Q.y;
    ctype const c24 = -M.c28 * Q.u + M.c30 * Q.z;
    ctype const c25 = -M.c29 * Q.u + M.c30 * Q.w;
    ctype const c26 = M.c31 * Q.x;
    ctype const c27 = M.c31 * Q.y;
    ctype const c28 = M.c31 * Q.z;
    ctype const c29 = M.c31 * Q.w;
    ctype const c30 = M.c31 * Q.u;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(quadvec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q,
                                                 MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * M.c1 + Q.y * M.c2 + Q.z * M.c3 + Q.w * M.c4 + Q.u * M.c5;
    ctype const c1 = Q.y * M.c11 - Q.z * M.c10 - Q.w * M.c6 + Q.u * M.c12;
    ctype const c2 = -Q.x * M.c11 + Q.z * M.c9 - Q.w * M.c7 + Q.u * M.c13;
    ctype const c3 = Q.x * M.c10 - Q.y * M.c9 - Q.w * M.c8 + Q.u * M.c14;
    ctype const c4 = Q.x * M.c6 + Q.y * M.c7 + Q.z * M.c8 + Q.u * M.c15;
    ctype const c5 = -Q.x * M.c12 - Q.y * M.c13 - Q.z * M.c14 - Q.w * M.c15;
    ctype const c6 = Q.y * M.c24 - Q.z * M.c23 + Q.u * M.c16;
    ctype const c7 = -Q.x * M.c24 + Q.z * M.c22 + Q.u * M.c17;
    ctype const c8 = Q.x * M.c23 - Q.y * M.c22 + Q.u * M.c18;
    ctype const c9 = -Q.x * M.c25 + Q.w * M.c22 + Q.u * M.c19;
    ctype const c10 = -Q.y * M.c25 + Q.w * M.c23 + Q.u * M.c20;
    ctype const c11 = -Q.z * M.c25 + Q.w * M.c24 + Q.u * M.c21;
    ctype const c12 = -Q.y * M.c21 + Q.z * M.c20 + Q.w * M.c16;
    ctype const c13 = Q.x * M.c21 - Q.z * M.c19 + Q.w * M.c17;
    ctype const c14 = -Q.x * M.c20 + Q.y * M.c19 + Q.w * M.c18;
    ctype const c15 = -Q.x * M.c16 - Q.y * M.c17 - Q.z * M.c18;
    ctype const c16 = -Q.y * M.c28 + Q.z * M.c27;
    ctype const c17 = Q.x * M.c28 - Q.z * M.c26;
    ctype const c18 = -Q.x * M.c27 + Q.y * M.c26;
    ctype const c19 = Q.x * M.c29 - Q.w * M.c26;
    ctype const c20 = Q.y * M.c29 - Q.w * M.c27;
    ctype const c21 = Q.z * M.c29 - Q.w * M.c28;
    ctype const c22 = -Q.x * M.c30 + Q.u * M.c26;
    ctype const c23 = -Q.y * M.c30 + Q.u * M.c27;
    ctype const c24 = -Q.z * M.c30 + Q.u * M.c28;
    ctype const c25 = -Q.w * M.c30 + Q.u * M.c29;
    ctype const c26 = Q.x * M.c31;
    ctype const c27 = Q.y * M.c31;
    ctype const c28 = Q.z * M.c31;
    ctype const c29 = Q.w * M.c31;
    ctype const c30 = Q.u * M.c31;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& M,
                                                 TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c6 * t.mx - M.c7 * t.my - M.c8 * t.mz - M.c9 * t.vx -
                     M.c10 * t.vy - M.c11 * t.vz - M.c12 * t.px - M.c13 * t.py -
                     M.c14 * t.pz - M.c15 * t.pw;
    ctype const c1 = M.c16 * t.pw + M.c20 * t.pz - M.c21 * t.py - M.c23 * t.mz +
                     M.c24 * t.my + M.c25 * t.vx;
    ctype const c2 = M.c17 * t.pw - M.c19 * t.pz + M.c21 * t.px + M.c22 * t.mz -
                     M.c24 * t.mx + M.c25 * t.vy;
    ctype const c3 = M.c18 * t.pw + M.c19 * t.py - M.c20 * t.px - M.c22 * t.my +
                     M.c23 * t.mx + M.c25 * t.vz;
    ctype const c4 = -M.c16 * t.px - M.c17 * t.py - M.c18 * t.pz - M.c22 * t.vx -
                     M.c23 * t.vy - M.c24 * t.vz;
    ctype const c5 = -M.c16 * t.mx - M.c17 * t.my - M.c18 * t.mz - M.c19 * t.vx -
                     M.c20 * t.vy - M.c21 * t.vz;
    ctype const c6 = M.c27 * t.pz - M.c28 * t.py + M.c30 * t.vx;
    ctype const c7 = -M.c26 * t.pz + M.c28 * t.px + M.c30 * t.vy;
    ctype const c8 = M.c26 * t.py - M.c27 * t.px + M.c30 * t.vz;
    ctype const c9 = -M.c26 * t.pw + M.c29 * t.px + M.c30 * t.mx;
    ctype const c10 = -M.c27 * t.pw + M.c29 * t.py + M.c30 * t.my;
    ctype const c11 = -M.c28 * t.pw + M.c29 * t.pz + M.c30 * t.mz;
    ctype const c12 = -M.c27 * t.mz + M.c28 * t.my + M.c29 * t.vx;
    ctype const c13 = M.c26 * t.mz - M.c28 * t.mx + M.c29 * t.vy;
    ctype const c14 = -M.c26 * t.my + M.c27 * t.mx + M.c29 * t.vz;
    ctype const c15 = -M.c26 * t.vx - M.c27 * t.vy - M.c28 * t.vz;
    ctype const c16 = M.c31 * t.vx;
    ctype const c17 = M.c31 * t.vy;
    ctype const c18 = M.c31 * t.vz;
    ctype const c19 = M.c31 * t.mx;
    ctype const c20 = M.c31 * t.my;
    ctype const c21 = M.c31 * t.mz;
    ctype const c22 = M.c31 * t.px;
    ctype const c23 = M.c31 * t.py;
    ctype const c24 = M.c31 * t.pz;
    ctype const c25 = M.c31 * t.pw;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(TriVec3dc<T> const& t,
                                                 MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c9 - t.vy * M.c10 - t.vz * M.c11 - t.mx * M.c6 -
                     t.my * M.c7 - t.mz * M.c8 - t.px * M.c12 - t.py * M.c13 -
                     t.pz * M.c14 - t.pw * M.c15;
    ctype const c1 = t.vx * M.c25 + t.my * M.c24 - t.mz * M.c23 - t.py * M.c21 +
                     t.pz * M.c20 + t.pw * M.c16;
    ctype const c2 = t.vy * M.c25 - t.mx * M.c24 + t.mz * M.c22 + t.px * M.c21 -
                     t.pz * M.c19 + t.pw * M.c17;
    ctype const c3 = t.vz * M.c25 + t.mx * M.c23 - t.my * M.c22 - t.px * M.c20 +
                     t.py * M.c19 + t.pw * M.c18;
    ctype const c4 = -t.vx * M.c22 - t.vy * M.c23 - t.vz * M.c24 - t.px * M.c16 -
                     t.py * M.c17 - t.pz * M.c18;
    ctype const c5 = -t.vx * M.c19 - t.vy * M.c20 - t.vz * M.c21 - t.mx * M.c16 -
                     t.my * M.c17 - t.mz * M.c18;
    ctype const c6 = t.vx * M.c30 - t.py * M.c28 + t.pz * M.c27;
    ctype const c7 = t.vy * M.c30 + t.px * M.c28 - t.pz * M.c26;
    ctype const c8 = t.vz * M.c30 - t.px * M.c27 + t.py * M.c26;
    ctype const c9 = t.mx * M.c30 + t.px * M.c29 - t.pw * M.c26;
    ctype const c10 = t.my * M.c30 + t.py * M.c29 - t.pw * M.c27;
    ctype const c11 = t.mz * M.c30 + t.pz * M.c29 - t.pw * M.c28;
    ctype const c12 = t.vx * M.c29 + t.my * M.c28 - t.mz * M.c27;
    ctype const c13 = t.vy * M.c29 - t.mx * M.c28 + t.mz * M.c26;
    ctype const c14 = t.vz * M.c29 + t.mx * M.c27 - t.my * M.c26;
    ctype const c15 = -t.vx * M.c26 - t.vy * M.c27 - t.vz * M.c28;
    ctype const c16 = t.vx * M.c31;
    ctype const c17 = t.vy * M.c31;
    ctype const c18 = t.vz * M.c31;
    ctype const c19 = t.mx * M.c31;
    ctype const c20 = t.my * M.c31;
    ctype const c21 = t.mz * M.c31;
    ctype const c22 = t.px * M.c31;
    ctype const c23 = t.py * M.c31;
    ctype const c24 = t.pz * M.c31;
    ctype const c25 = t.pw * M.c31;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& M,
                                                 BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c16 * B.mx - M.c17 * B.my - M.c18 * B.mz - M.c19 * B.vx -
                     M.c20 * B.vy - M.c21 * B.vz - M.c22 * B.px - M.c23 * B.py -
                     M.c24 * B.pz - M.c25 * B.pw;
    ctype const c1 = M.c27 * B.mz - M.c28 * B.my - M.c29 * B.vx + M.c30 * B.px;
    ctype const c2 = -M.c26 * B.mz + M.c28 * B.mx - M.c29 * B.vy + M.c30 * B.py;
    ctype const c3 = M.c26 * B.my - M.c27 * B.mx - M.c29 * B.vz + M.c30 * B.pz;
    ctype const c4 = M.c26 * B.vx + M.c27 * B.vy + M.c28 * B.vz + M.c30 * B.pw;
    ctype const c5 = -M.c26 * B.px - M.c27 * B.py - M.c28 * B.pz - M.c29 * B.pw;
    ctype const c6 = M.c31 * B.vx;
    ctype const c7 = M.c31 * B.vy;
    ctype const c8 = M.c31 * B.vz;
    ctype const c9 = M.c31 * B.mx;
    ctype const c10 = M.c31 * B.my;
    ctype const c11 = M.c31 * B.mz;
    ctype const c12 = M.c31 * B.px;
    ctype const c13 = M.c31 * B.py;
    ctype const c14 = M.c31 * B.pz;
    ctype const c15 = M.c31 * B.pw;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(BiVec3dc<T> const& B,
                                                 MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c19 - B.vy * M.c20 - B.vz * M.c21 - B.mx * M.c16 -
                     B.my * M.c17 - B.mz * M.c18 - B.px * M.c22 - B.py * M.c23 -
                     B.pz * M.c24 - B.pw * M.c25;
    ctype const c1 = B.vx * M.c29 + B.my * M.c28 - B.mz * M.c27 - B.px * M.c30;
    ctype const c2 = B.vy * M.c29 - B.mx * M.c28 + B.mz * M.c26 - B.py * M.c30;
    ctype const c3 = B.vz * M.c29 + B.mx * M.c27 - B.my * M.c26 - B.pz * M.c30;
    ctype const c4 = -B.vx * M.c26 - B.vy * M.c27 - B.vz * M.c28 - B.pw * M.c30;
    ctype const c5 = B.px * M.c26 + B.py * M.c27 + B.pz * M.c28 + B.pw * M.c29;
    ctype const c6 = B.vx * M.c31;
    ctype const c7 = B.vy * M.c31;
    ctype const c8 = B.vz * M.c31;
    ctype const c9 = B.mx * M.c31;
    ctype const c10 = B.my * M.c31;
    ctype const c11 = B.mz * M.c31;
    ctype const c12 = B.px * M.c31;
    ctype const c13 = B.py * M.c31;
    ctype const c14 = B.pz * M.c31;
    ctype const c15 = B.pw * M.c31;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& M, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c26 * v.x + M.c27 * v.y + M.c28 * v.z + M.c29 * v.w + M.c30 * v.u;
    ctype const c1 = M.c31 * v.x;
    ctype const c2 = M.c31 * v.y;
    ctype const c3 = M.c31 * v.z;
    ctype const c4 = M.c31 * v.w;
    ctype const c5 = M.c31 * v.u;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(Vec3dc<T> const& v, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c26 + v.y * M.c27 + v.z * M.c28 + v.w * M.c29 + v.u * M.c30;
    ctype const c1 = v.x * M.c31;
    ctype const c2 = v.y * M.c31;
    ctype const c3 = v.z * M.c31;
    ctype const c4 = v.w * M.c31;
    ctype const c5 = v.u * M.c31;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(MVec3dc<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c31 * ctype(s);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rwdg(Scalar3dc<T> s, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c31;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rwdg :: rwdg(mv_e,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& A,
                                                   MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c14 + A.c5 * B.c13 - A.c6 * B.c12 - A.c7 * B.c15 +
                     A.c12 * B.c6 - A.c13 * B.c5 - A.c14 * B.c1 + A.c15 * B.c7;
    ctype const c1 = A.c2 * B.c14 - A.c4 * B.c13 + A.c6 * B.c11 - A.c8 * B.c15 -
                     A.c11 * B.c6 + A.c13 * B.c4 - A.c14 * B.c2 + A.c15 * B.c8;
    ctype const c2 = A.c3 * B.c14 + A.c4 * B.c12 - A.c5 * B.c11 - A.c9 * B.c15 +
                     A.c11 * B.c5 - A.c12 * B.c4 - A.c14 * B.c3 + A.c15 * B.c9;
    ctype const c3 = -A.c1 * B.c11 - A.c2 * B.c12 - A.c3 * B.c13 - A.c10 * B.c15 +
                     A.c11 * B.c1 + A.c12 * B.c2 + A.c13 * B.c3 + A.c15 * B.c10;
    ctype const c4 = A.c7 * B.c11 + A.c8 * B.c12 + A.c9 * B.c13 + A.c10 * B.c14 -
                     A.c11 * B.c7 - A.c12 * B.c8 - A.c13 * B.c9 - A.c14 * B.c10;
    ctype const c5 = -A.c12 * B.c13 + A.c13 * B.c12;
    ctype const c6 = A.c11 * B.c13 - A.c13 * B.c11;
    ctype const c7 = -A.c11 * B.c12 + A.c12 * B.c11;
    ctype const c8 = A.c11 * B.c14 - A.c14 * B.c11;
    ctype const c9 = A.c12 * B.c14 - A.c14 * B.c12;
    ctype const c10 = A.c13 * B.c14 - A.c14 * B.c13;
    ctype const c11 = -A.c11 * B.c15 + A.c15 * B.c11;
    ctype const c12 = -A.c12 * B.c15 + A.c15 * B.c12;
    ctype const c13 = -A.c13 * B.c15 + A.c15 * B.c13;
    ctype const c14 = -A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(mv_e,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& A,
                                                   MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                     A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 - A.c7 * B.c11 -
                     A.c8 * B.c12 - A.c9 * B.c13 - A.c10 * B.c14 + A.c11 * B.c0 +
                     A.c12 * B.c1 + A.c13 * B.c2 + A.c14 * B.c3 + A.c15 * B.c4;
    ctype const c1 = A.c1 * B.c15 + A.c12 * B.c13 - A.c13 * B.c12 + A.c15 * B.c5;
    ctype const c2 = A.c2 * B.c15 - A.c11 * B.c13 + A.c13 * B.c11 + A.c15 * B.c6;
    ctype const c3 = A.c3 * B.c15 + A.c11 * B.c12 - A.c12 * B.c11 + A.c15 * B.c7;
    ctype const c4 = A.c4 * B.c15 - A.c11 * B.c14 + A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c5 = A.c5 * B.c15 - A.c12 * B.c14 + A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c6 = A.c6 * B.c15 - A.c13 * B.c14 + A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c7 = A.c7 * B.c15 - A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5;
    ctype const c8 = A.c8 * B.c15 + A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c6;
    ctype const c9 = A.c9 * B.c15 - A.c11 * B.c9 + A.c12 * B.c8 + A.c14 * B.c7;
    ctype const c10 = A.c10 * B.c15 - A.c11 * B.c5 - A.c12 * B.c6 - A.c13 * B.c7;
    ctype const c11 = A.c11 * B.c15;
    ctype const c12 = A.c12 * B.c15;
    ctype const c13 = A.c13 * B.c15;
    ctype const c14 = A.c14 * B.c15;
    ctype const c15 = A.c15 * B.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(mv_u,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& A,
                                                   MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                     A.c4 * B.c15 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 -
                     A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 - A.c11 * B.c7 -
                     A.c12 * B.c8 - A.c13 * B.c9 - A.c14 * B.c10 + A.c15 * B.c0;
    ctype const c1 = A.c5 * B.c15 - A.c12 * B.c13 + A.c13 * B.c12 + A.c15 * B.c1;
    ctype const c2 = A.c6 * B.c15 + A.c11 * B.c13 - A.c13 * B.c11 + A.c15 * B.c2;
    ctype const c3 = A.c7 * B.c15 - A.c11 * B.c12 + A.c12 * B.c11 + A.c15 * B.c3;
    ctype const c4 = A.c8 * B.c15 + A.c11 * B.c14 - A.c14 * B.c11 + A.c15 * B.c4;
    ctype const c5 = A.c9 * B.c15 + A.c12 * B.c14 - A.c14 * B.c12 + A.c15 * B.c5;
    ctype const c6 = A.c10 * B.c15 + A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c6;
    ctype const c7 = A.c5 * B.c14 + A.c9 * B.c13 - A.c10 * B.c12 + A.c15 * B.c7;
    ctype const c8 = A.c6 * B.c14 - A.c8 * B.c13 + A.c10 * B.c11 + A.c15 * B.c8;
    ctype const c9 = A.c7 * B.c14 + A.c8 * B.c12 - A.c9 * B.c11 + A.c15 * B.c9;
    ctype const c10 = -A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 + A.c15 * B.c10;
    ctype const c11 = A.c15 * B.c11;
    ctype const c12 = A.c15 * B.c12;
    ctype const c13 = A.c15 * B.c13;
    ctype const c14 = A.c15 * B.c14;
    ctype const c15 = A.c15 * B.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& M,
                                                   PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    ctype const c4 = M.c4 * ctype(ps);
    ctype const c5 = M.c5 * ctype(ps);
    ctype const c6 = M.c6 * ctype(ps);
    ctype const c7 = M.c7 * ctype(ps);
    ctype const c8 = M.c8 * ctype(ps);
    ctype const c9 = M.c9 * ctype(ps);
    ctype const c10 = M.c10 * ctype(ps);
    ctype const c11 = M.c11 * ctype(ps);
    ctype const c12 = M.c12 * ctype(ps);
    ctype const c13 = M.c13 * ctype(ps);
    ctype const c14 = M.c14 * ctype(ps);
    ctype const c15 = M.c15 * ctype(ps);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    ctype const c4 = ctype(ps) * M.c4;
    ctype const c5 = ctype(ps) * M.c5;
    ctype const c6 = ctype(ps) * M.c6;
    ctype const c7 = ctype(ps) * M.c7;
    ctype const c8 = ctype(ps) * M.c8;
    ctype const c9 = ctype(ps) * M.c9;
    ctype const c10 = ctype(ps) * M.c10;
    ctype const c11 = ctype(ps) * M.c11;
    ctype const c12 = ctype(ps) * M.c12;
    ctype const c13 = ctype(ps) * M.c13;
    ctype const c14 = ctype(ps) * M.c14;
    ctype const c15 = ctype(ps) * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(mv_e,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& M,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * Q.w + M.c5 * Q.z - M.c6 * Q.y - M.c7 * Q.u;
    ctype const c1 = M.c2 * Q.w - M.c4 * Q.z + M.c6 * Q.x - M.c8 * Q.u;
    ctype const c2 = M.c3 * Q.w + M.c4 * Q.y - M.c5 * Q.x - M.c9 * Q.u;
    ctype const c3 = -M.c1 * Q.x - M.c2 * Q.y - M.c3 * Q.z - M.c10 * Q.u;
    ctype const c4 = M.c7 * Q.x + M.c8 * Q.y + M.c9 * Q.z + M.c10 * Q.w;
    ctype const c5 = -M.c12 * Q.z + M.c13 * Q.y;
    ctype const c6 = M.c11 * Q.z - M.c13 * Q.x;
    ctype const c7 = -M.c11 * Q.y + M.c12 * Q.x;
    ctype const c8 = M.c11 * Q.w - M.c14 * Q.x;
    ctype const c9 = M.c12 * Q.w - M.c14 * Q.y;
    ctype const c10 = M.c13 * Q.w - M.c14 * Q.z;
    ctype const c11 = -M.c11 * Q.u + M.c15 * Q.x;
    ctype const c12 = -M.c12 * Q.u + M.c15 * Q.y;
    ctype const c13 = -M.c13 * Q.u + M.c15 * Q.z;
    ctype const c14 = -M.c14 * Q.u + M.c15 * Q.w;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(quadvec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * M.c6 - Q.z * M.c5 - Q.w * M.c1 + Q.u * M.c7;
    ctype const c1 = -Q.x * M.c6 + Q.z * M.c4 - Q.w * M.c2 + Q.u * M.c8;
    ctype const c2 = Q.x * M.c5 - Q.y * M.c4 - Q.w * M.c3 + Q.u * M.c9;
    ctype const c3 = Q.x * M.c1 + Q.y * M.c2 + Q.z * M.c3 + Q.u * M.c10;
    ctype const c4 = -Q.x * M.c7 - Q.y * M.c8 - Q.z * M.c9 - Q.w * M.c10;
    ctype const c5 = -Q.y * M.c13 + Q.z * M.c12;
    ctype const c6 = Q.x * M.c13 - Q.z * M.c11;
    ctype const c7 = -Q.x * M.c12 + Q.y * M.c11;
    ctype const c8 = Q.x * M.c14 - Q.w * M.c11;
    ctype const c9 = Q.y * M.c14 - Q.w * M.c12;
    ctype const c10 = Q.z * M.c14 - Q.w * M.c13;
    ctype const c11 = -Q.x * M.c15 + Q.u * M.c11;
    ctype const c12 = -Q.y * M.c15 + Q.u * M.c12;
    ctype const c13 = -Q.z * M.c15 + Q.u * M.c13;
    ctype const c14 = -Q.w * M.c15 + Q.u * M.c14;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(mv_e,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& M,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz - M.c4 * t.vx -
                     M.c5 * t.vy - M.c6 * t.vz - M.c7 * t.px - M.c8 * t.py - M.c9 * t.pz -
                     M.c10 * t.pw;
    ctype const c1 = M.c12 * t.pz - M.c13 * t.py + M.c15 * t.vx;
    ctype const c2 = -M.c11 * t.pz + M.c13 * t.px + M.c15 * t.vy;
    ctype const c3 = M.c11 * t.py - M.c12 * t.px + M.c15 * t.vz;
    ctype const c4 = -M.c11 * t.pw + M.c14 * t.px + M.c15 * t.mx;
    ctype const c5 = -M.c12 * t.pw + M.c14 * t.py + M.c15 * t.my;
    ctype const c6 = -M.c13 * t.pw + M.c14 * t.pz + M.c15 * t.mz;
    ctype const c7 = -M.c12 * t.mz + M.c13 * t.my + M.c14 * t.vx;
    ctype const c8 = M.c11 * t.mz - M.c13 * t.mx + M.c14 * t.vy;
    ctype const c9 = -M.c11 * t.my + M.c12 * t.mx + M.c14 * t.vz;
    ctype const c10 = -M.c11 * t.vx - M.c12 * t.vy - M.c13 * t.vz;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(trivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(TriVec3dc<T> const& t,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c4 - t.vy * M.c5 - t.vz * M.c6 - t.mx * M.c1 -
                     t.my * M.c2 - t.mz * M.c3 - t.px * M.c7 - t.py * M.c8 - t.pz * M.c9 -
                     t.pw * M.c10;
    ctype const c1 = t.vx * M.c15 - t.py * M.c13 + t.pz * M.c12;
    ctype const c2 = t.vy * M.c15 + t.px * M.c13 - t.pz * M.c11;
    ctype const c3 = t.vz * M.c15 - t.px * M.c12 + t.py * M.c11;
    ctype const c4 = t.mx * M.c15 + t.px * M.c14 - t.pw * M.c11;
    ctype const c5 = t.my * M.c15 + t.py * M.c14 - t.pw * M.c12;
    ctype const c6 = t.mz * M.c15 + t.pz * M.c14 - t.pw * M.c13;
    ctype const c7 = t.vx * M.c14 + t.my * M.c13 - t.mz * M.c12;
    ctype const c8 = t.vy * M.c14 - t.mx * M.c13 + t.mz * M.c11;
    ctype const c9 = t.vz * M.c14 + t.mx * M.c12 - t.my * M.c11;
    ctype const c10 = -t.vx * M.c11 - t.vy * M.c12 - t.vz * M.c13;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(mv_e,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& M,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c12 * B.mz - M.c13 * B.my - M.c14 * B.vx + M.c15 * B.px;
    ctype const c1 = -M.c11 * B.mz + M.c13 * B.mx - M.c14 * B.vy + M.c15 * B.py;
    ctype const c2 = M.c11 * B.my - M.c12 * B.mx - M.c14 * B.vz + M.c15 * B.pz;
    ctype const c3 = M.c11 * B.vx + M.c12 * B.vy + M.c13 * B.vz + M.c15 * B.pw;
    ctype const c4 = -M.c11 * B.px - M.c12 * B.py - M.c13 * B.pz - M.c14 * B.pw;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(bivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(BiVec3dc<T> const& B,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c14 + B.my * M.c13 - B.mz * M.c12 - B.px * M.c15;
    ctype const c1 = B.vy * M.c14 - B.mx * M.c13 + B.mz * M.c11 - B.py * M.c15;
    ctype const c2 = B.vz * M.c14 + B.mx * M.c12 - B.my * M.c11 - B.pz * M.c15;
    ctype const c3 = -B.vx * M.c11 - B.vy * M.c12 - B.vz * M.c13 - B.pw * M.c15;
    ctype const c4 = B.px * M.c11 + B.py * M.c12 + B.pz * M.c13 + B.pw * M.c14;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(mv_e,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_E<T> const& M,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c11 * v.x + M.c12 * v.y + M.c13 * v.z + M.c14 * v.w + M.c15 * v.u;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(vec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(Vec3dc<T> const& v,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14 + v.u * M.c15;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(mv_e,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] MVec3dc_E<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(s,mv_e) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] MVec3dc_E<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(mv_u,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& A,
                                                   MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 + A.c5 * B.c14 + A.c9 * B.c13 - A.c10 * B.c12 -
                     A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 + A.c15 * B.c0;
    ctype const c1 = A.c1 * B.c15 + A.c6 * B.c14 - A.c8 * B.c13 + A.c10 * B.c11 +
                     A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c6 + A.c15 * B.c1;
    ctype const c2 = A.c2 * B.c15 + A.c7 * B.c14 + A.c8 * B.c12 - A.c9 * B.c11 -
                     A.c11 * B.c9 + A.c12 * B.c8 + A.c14 * B.c7 + A.c15 * B.c2;
    ctype const c3 = A.c3 * B.c15 - A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 -
                     A.c11 * B.c5 - A.c12 * B.c6 - A.c13 * B.c7 + A.c15 * B.c3;
    ctype const c4 = A.c4 * B.c15 - A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 -
                     A.c8 * B.c5 - A.c9 * B.c6 - A.c10 * B.c7 + A.c15 * B.c4;
    ctype const c5 = A.c5 * B.c15 + A.c15 * B.c5;
    ctype const c6 = A.c6 * B.c15 + A.c15 * B.c6;
    ctype const c7 = A.c7 * B.c15 + A.c15 * B.c7;
    ctype const c8 = A.c8 * B.c15 + A.c15 * B.c8;
    ctype const c9 = A.c9 * B.c15 + A.c15 * B.c9;
    ctype const c10 = A.c10 * B.c15 + A.c15 * B.c10;
    ctype const c11 = A.c11 * B.c15 + A.c15 * B.c11;
    ctype const c12 = A.c12 * B.c15 + A.c15 * B.c12;
    ctype const c13 = A.c13 * B.c15 + A.c15 * B.c13;
    ctype const c14 = A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = A.c15 * B.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& M,
                                                   PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    ctype const c4 = M.c4 * ctype(ps);
    ctype const c5 = M.c5 * ctype(ps);
    ctype const c6 = M.c6 * ctype(ps);
    ctype const c7 = M.c7 * ctype(ps);
    ctype const c8 = M.c8 * ctype(ps);
    ctype const c9 = M.c9 * ctype(ps);
    ctype const c10 = M.c10 * ctype(ps);
    ctype const c11 = M.c11 * ctype(ps);
    ctype const c12 = M.c12 * ctype(ps);
    ctype const c13 = M.c13 * ctype(ps);
    ctype const c14 = M.c14 * ctype(ps);
    ctype const c15 = M.c15 * ctype(ps);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    ctype const c4 = ctype(ps) * M.c4;
    ctype const c5 = ctype(ps) * M.c5;
    ctype const c6 = ctype(ps) * M.c6;
    ctype const c7 = ctype(ps) * M.c7;
    ctype const c8 = ctype(ps) * M.c8;
    ctype const c9 = ctype(ps) * M.c9;
    ctype const c10 = ctype(ps) * M.c10;
    ctype const c11 = ctype(ps) * M.c11;
    ctype const c12 = ctype(ps) * M.c12;
    ctype const c13 = ctype(ps) * M.c13;
    ctype const c14 = ctype(ps) * M.c14;
    ctype const c15 = ctype(ps) * M.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(mv_u,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& M,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * Q.x + M.c1 * Q.y + M.c2 * Q.z + M.c3 * Q.w + M.c4 * Q.u;
    ctype const c1 = M.c5 * Q.u - M.c12 * Q.z + M.c13 * Q.y;
    ctype const c2 = M.c6 * Q.u + M.c11 * Q.z - M.c13 * Q.x;
    ctype const c3 = M.c7 * Q.u - M.c11 * Q.y + M.c12 * Q.x;
    ctype const c4 = M.c8 * Q.u + M.c11 * Q.w - M.c14 * Q.x;
    ctype const c5 = M.c9 * Q.u + M.c12 * Q.w - M.c14 * Q.y;
    ctype const c6 = M.c10 * Q.u + M.c13 * Q.w - M.c14 * Q.z;
    ctype const c7 = M.c5 * Q.w + M.c9 * Q.z - M.c10 * Q.y;
    ctype const c8 = M.c6 * Q.w - M.c8 * Q.z + M.c10 * Q.x;
    ctype const c9 = M.c7 * Q.w + M.c8 * Q.y - M.c9 * Q.x;
    ctype const c10 = -M.c5 * Q.x - M.c6 * Q.y - M.c7 * Q.z;
    ctype const c11 = M.c15 * Q.x;
    ctype const c12 = M.c15 * Q.y;
    ctype const c13 = M.c15 * Q.z;
    ctype const c14 = M.c15 * Q.w;
    ctype const c15 = M.c15 * Q.u;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(quadvec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * M.c0 + Q.y * M.c1 + Q.z * M.c2 + Q.w * M.c3 + Q.u * M.c4;
    ctype const c1 = Q.y * M.c13 - Q.z * M.c12 + Q.u * M.c5;
    ctype const c2 = -Q.x * M.c13 + Q.z * M.c11 + Q.u * M.c6;
    ctype const c3 = Q.x * M.c12 - Q.y * M.c11 + Q.u * M.c7;
    ctype const c4 = -Q.x * M.c14 + Q.w * M.c11 + Q.u * M.c8;
    ctype const c5 = -Q.y * M.c14 + Q.w * M.c12 + Q.u * M.c9;
    ctype const c6 = -Q.z * M.c14 + Q.w * M.c13 + Q.u * M.c10;
    ctype const c7 = -Q.y * M.c10 + Q.z * M.c9 + Q.w * M.c5;
    ctype const c8 = Q.x * M.c10 - Q.z * M.c8 + Q.w * M.c6;
    ctype const c9 = -Q.x * M.c9 + Q.y * M.c8 + Q.w * M.c7;
    ctype const c10 = -Q.x * M.c5 - Q.y * M.c6 - Q.z * M.c7;
    ctype const c11 = Q.x * M.c15;
    ctype const c12 = Q.y * M.c15;
    ctype const c13 = Q.z * M.c15;
    ctype const c14 = Q.w * M.c15;
    ctype const c15 = Q.u * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(mv_u,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& M,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * t.pw + M.c9 * t.pz - M.c10 * t.py - M.c12 * t.mz +
                     M.c13 * t.my + M.c14 * t.vx;
    ctype const c1 = M.c6 * t.pw - M.c8 * t.pz + M.c10 * t.px + M.c11 * t.mz -
                     M.c13 * t.mx + M.c14 * t.vy;
    ctype const c2 = M.c7 * t.pw + M.c8 * t.py - M.c9 * t.px - M.c11 * t.my +
                     M.c12 * t.mx + M.c14 * t.vz;
    ctype const c3 = -M.c5 * t.px - M.c6 * t.py - M.c7 * t.pz - M.c11 * t.vx -
                     M.c12 * t.vy - M.c13 * t.vz;
    ctype const c4 = -M.c5 * t.mx - M.c6 * t.my - M.c7 * t.mz - M.c8 * t.vx -
                     M.c9 * t.vy - M.c10 * t.vz;
    ctype const c5 = M.c15 * t.vx;
    ctype const c6 = M.c15 * t.vy;
    ctype const c7 = M.c15 * t.vz;
    ctype const c8 = M.c15 * t.mx;
    ctype const c9 = M.c15 * t.my;
    ctype const c10 = M.c15 * t.mz;
    ctype const c11 = M.c15 * t.px;
    ctype const c12 = M.c15 * t.py;
    ctype const c13 = M.c15 * t.pz;
    ctype const c14 = M.c15 * t.pw;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(trivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(TriVec3dc<T> const& t,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * M.c14 + t.my * M.c13 - t.mz * M.c12 - t.py * M.c10 +
                     t.pz * M.c9 + t.pw * M.c5;
    ctype const c1 = t.vy * M.c14 - t.mx * M.c13 + t.mz * M.c11 + t.px * M.c10 -
                     t.pz * M.c8 + t.pw * M.c6;
    ctype const c2 = t.vz * M.c14 + t.mx * M.c12 - t.my * M.c11 - t.px * M.c9 +
                     t.py * M.c8 + t.pw * M.c7;
    ctype const c3 = -t.vx * M.c11 - t.vy * M.c12 - t.vz * M.c13 - t.px * M.c5 -
                     t.py * M.c6 - t.pz * M.c7;
    ctype const c4 = -t.vx * M.c8 - t.vy * M.c9 - t.vz * M.c10 - t.mx * M.c5 -
                     t.my * M.c6 - t.mz * M.c7;
    ctype const c5 = t.vx * M.c15;
    ctype const c6 = t.vy * M.c15;
    ctype const c7 = t.vz * M.c15;
    ctype const c8 = t.mx * M.c15;
    ctype const c9 = t.my * M.c15;
    ctype const c10 = t.mz * M.c15;
    ctype const c11 = t.px * M.c15;
    ctype const c12 = t.py * M.c15;
    ctype const c13 = t.pz * M.c15;
    ctype const c14 = t.pw * M.c15;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(mv_u,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& M,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                     M.c9 * B.vy - M.c10 * B.vz - M.c11 * B.px - M.c12 * B.py -
                     M.c13 * B.pz - M.c14 * B.pw;
    ctype const c1 = M.c15 * B.vx;
    ctype const c2 = M.c15 * B.vy;
    ctype const c3 = M.c15 * B.vz;
    ctype const c4 = M.c15 * B.mx;
    ctype const c5 = M.c15 * B.my;
    ctype const c6 = M.c15 * B.mz;
    ctype const c7 = M.c15 * B.px;
    ctype const c8 = M.c15 * B.py;
    ctype const c9 = M.c15 * B.pz;
    ctype const c10 = M.c15 * B.pw;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(bivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(BiVec3dc<T> const& B,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                     B.my * M.c6 - B.mz * M.c7 - B.px * M.c11 - B.py * M.c12 -
                     B.pz * M.c13 - B.pw * M.c14;
    ctype const c1 = B.vx * M.c15;
    ctype const c2 = B.vy * M.c15;
    ctype const c3 = B.vz * M.c15;
    ctype const c4 = B.mx * M.c15;
    ctype const c5 = B.my * M.c15;
    ctype const c6 = B.mz * M.c15;
    ctype const c7 = B.px * M.c15;
    ctype const c8 = B.py * M.c15;
    ctype const c9 = B.pz * M.c15;
    ctype const c10 = B.pw * M.c15;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(mv_u,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& M,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c15 * v.x;
    ctype const c1 = M.c15 * v.y;
    ctype const c2 = M.c15 * v.z;
    ctype const c3 = M.c15 * v.w;
    ctype const c4 = M.c15 * v.u;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(vec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rwdg(Vec3dc<T> const& v,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c15;
    ctype const c1 = v.y * M.c15;
    ctype const c2 = v.z * M.c15;
    ctype const c3 = v.w * M.c15;
    ctype const c4 = v.u * M.c15;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rwdg :: rwdg(mv_u,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(MVec3dc_U<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c15 * ctype(s);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(s,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rwdg(Scalar3dc<T> s, MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c15;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rwdg :: rwdg(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps1, PScalar3dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(ps1) * ctype(ps2));
}

// cga3dc rwdg :: rwdg(ps,quadvec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps,
                                                    QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * Q.x;
    ctype const c1 = ctype(ps) * Q.y;
    ctype const c2 = ctype(ps) * Q.z;
    ctype const c3 = ctype(ps) * Q.w;
    ctype const c4 = ctype(ps) * Q.u;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rwdg :: rwdg(quadvec,ps) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q,
                                                    PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * ctype(ps);
    ctype const c1 = Q.y * ctype(ps);
    ctype const c2 = Q.z * ctype(ps);
    ctype const c3 = Q.w * ctype(ps);
    ctype const c4 = Q.u * ctype(ps);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rwdg :: rwdg(ps,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * t.vx;
    ctype const c1 = ctype(ps) * t.vy;
    ctype const c2 = ctype(ps) * t.vz;
    ctype const c3 = ctype(ps) * t.mx;
    ctype const c4 = ctype(ps) * t.my;
    ctype const c5 = ctype(ps) * t.mz;
    ctype const c6 = ctype(ps) * t.px;
    ctype const c7 = ctype(ps) * t.py;
    ctype const c8 = ctype(ps) * t.pz;
    ctype const c9 = ctype(ps) * t.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rwdg :: rwdg(trivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rwdg(TriVec3dc<T> const& t,
                                                   PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * ctype(ps);
    ctype const c1 = t.vy * ctype(ps);
    ctype const c2 = t.vz * ctype(ps);
    ctype const c3 = t.mx * ctype(ps);
    ctype const c4 = t.my * ctype(ps);
    ctype const c5 = t.mz * ctype(ps);
    ctype const c6 = t.px * ctype(ps);
    ctype const c7 = t.py * ctype(ps);
    ctype const c8 = t.pz * ctype(ps);
    ctype const c9 = t.pw * ctype(ps);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rwdg :: rwdg(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps, BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.vx;
    ctype const c1 = ctype(ps) * B.vy;
    ctype const c2 = ctype(ps) * B.vz;
    ctype const c3 = ctype(ps) * B.mx;
    ctype const c4 = ctype(ps) * B.my;
    ctype const c5 = ctype(ps) * B.mz;
    ctype const c6 = ctype(ps) * B.px;
    ctype const c7 = ctype(ps) * B.py;
    ctype const c8 = ctype(ps) * B.pz;
    ctype const c9 = ctype(ps) * B.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rwdg :: rwdg(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rwdg(BiVec3dc<T> const& B, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(ps);
    ctype const c1 = B.vy * ctype(ps);
    ctype const c2 = B.vz * ctype(ps);
    ctype const c3 = B.mx * ctype(ps);
    ctype const c4 = B.my * ctype(ps);
    ctype const c5 = B.mz * ctype(ps);
    ctype const c6 = B.px * ctype(ps);
    ctype const c7 = B.py * ctype(ps);
    ctype const c8 = B.pz * ctype(ps);
    ctype const c9 = B.pw * ctype(ps);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rwdg :: rwdg(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = ctype(ps) * v.w;
    ctype const c4 = ctype(ps) * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rwdg :: rwdg(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rwdg(Vec3dc<T> const& v, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = v.w * ctype(ps);
    ctype const c4 = v.u * ctype(ps);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rwdg :: rwdg(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg(PScalar3dc<T> ps, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(ps) * ctype(s));
}

// cga3dc rwdg :: rwdg(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg(Scalar3dc<T> s, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(s) * ctype(ps));
}

// cga3dc rwdg :: rwdg(quadvec,quadvec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q1,
                                                   QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q1.y * Q2.z + Q1.z * Q2.y;
    ctype const c1 = Q1.x * Q2.z - Q1.z * Q2.x;
    ctype const c2 = -Q1.x * Q2.y + Q1.y * Q2.x;
    ctype const c3 = Q1.x * Q2.w - Q1.w * Q2.x;
    ctype const c4 = Q1.y * Q2.w - Q1.w * Q2.y;
    ctype const c5 = Q1.z * Q2.w - Q1.w * Q2.z;
    ctype const c6 = -Q1.x * Q2.u + Q1.u * Q2.x;
    ctype const c7 = -Q1.y * Q2.u + Q1.u * Q2.y;
    ctype const c8 = -Q1.z * Q2.u + Q1.u * Q2.z;
    ctype const c9 = -Q1.w * Q2.u + Q1.u * Q2.w;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rwdg :: rwdg(quadvec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q,
                                                  TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * t.pz - Q.z * t.py + Q.u * t.vx;
    ctype const c1 = -Q.x * t.pz + Q.z * t.px + Q.u * t.vy;
    ctype const c2 = Q.x * t.py - Q.y * t.px + Q.u * t.vz;
    ctype const c3 = -Q.x * t.pw + Q.w * t.px + Q.u * t.mx;
    ctype const c4 = -Q.y * t.pw + Q.w * t.py + Q.u * t.my;
    ctype const c5 = -Q.z * t.pw + Q.w * t.pz + Q.u * t.mz;
    ctype const c6 = -Q.y * t.mz + Q.z * t.my + Q.w * t.vx;
    ctype const c7 = Q.x * t.mz - Q.z * t.mx + Q.w * t.vy;
    ctype const c8 = -Q.x * t.my + Q.y * t.mx + Q.w * t.vz;
    ctype const c9 = -Q.x * t.vx - Q.y * t.vy - Q.z * t.vz;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rwdg :: rwdg(trivec,quadvec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rwdg(TriVec3dc<T> const& t,
                                                  QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * Q.u - t.py * Q.z + t.pz * Q.y;
    ctype const c1 = t.vy * Q.u + t.px * Q.z - t.pz * Q.x;
    ctype const c2 = t.vz * Q.u - t.px * Q.y + t.py * Q.x;
    ctype const c3 = t.mx * Q.u + t.px * Q.w - t.pw * Q.x;
    ctype const c4 = t.my * Q.u + t.py * Q.w - t.pw * Q.y;
    ctype const c5 = t.mz * Q.u + t.pz * Q.w - t.pw * Q.z;
    ctype const c6 = t.vx * Q.w + t.my * Q.z - t.mz * Q.y;
    ctype const c7 = t.vy * Q.w - t.mx * Q.z + t.mz * Q.x;
    ctype const c8 = t.vz * Q.w + t.mx * Q.y - t.my * Q.x;
    ctype const c9 = -t.vx * Q.x - t.vy * Q.y - t.vz * Q.z;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rwdg :: rwdg(quadvec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q,
                                                BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * B.mz - Q.z * B.my - Q.w * B.vx + Q.u * B.px;
    ctype const c1 = -Q.x * B.mz + Q.z * B.mx - Q.w * B.vy + Q.u * B.py;
    ctype const c2 = Q.x * B.my - Q.y * B.mx - Q.w * B.vz + Q.u * B.pz;
    ctype const c3 = Q.x * B.vx + Q.y * B.vy + Q.z * B.vz + Q.u * B.pw;
    ctype const c4 = -Q.x * B.px - Q.y * B.py - Q.z * B.pz - Q.w * B.pw;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rwdg :: rwdg(bivec,quadvec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rwdg(BiVec3dc<T> const& B,
                                                QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * Q.w + B.my * Q.z - B.mz * Q.y - B.px * Q.u;
    ctype const c1 = B.vy * Q.w - B.mx * Q.z + B.mz * Q.x - B.py * Q.u;
    ctype const c2 = B.vz * Q.w + B.mx * Q.y - B.my * Q.x - B.pz * Q.u;
    ctype const c3 = -B.vx * Q.x - B.vy * Q.y - B.vz * Q.z - B.pw * Q.u;
    ctype const c4 = B.px * Q.x + B.py * Q.y + B.pz * Q.z + B.pw * Q.w;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rwdg :: rwdg(quadvec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg(QuadVec3dc<T> const& Q,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(Q.x * v.x + Q.y * v.y + Q.z * v.z + Q.w * v.w + Q.u * v.u);
}

// cga3dc rwdg :: rwdg(vec,quadvec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg(Vec3dc<T> const& v,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(v.x * Q.x + v.y * Q.y + v.z * Q.z + v.w * Q.w + v.u * Q.u);
}

// cga3dc rwdg :: rwdg(quadvec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] QuadVec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(s,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(trivec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rwdg(TriVec3dc<T> const& t1,
                                                TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.vx * t2.pw + t1.my * t2.pz - t1.mz * t2.py - t1.py * t2.mz +
                     t1.pz * t2.my + t1.pw * t2.vx;
    ctype const c1 = t1.vy * t2.pw - t1.mx * t2.pz + t1.mz * t2.px + t1.px * t2.mz -
                     t1.pz * t2.mx + t1.pw * t2.vy;
    ctype const c2 = t1.vz * t2.pw + t1.mx * t2.py - t1.my * t2.px - t1.px * t2.my +
                     t1.py * t2.mx + t1.pw * t2.vz;
    ctype const c3 = -t1.vx * t2.px - t1.vy * t2.py - t1.vz * t2.pz - t1.px * t2.vx -
                     t1.py * t2.vy - t1.pz * t2.vz;
    ctype const c4 = -t1.vx * t2.mx - t1.vy * t2.my - t1.vz * t2.mz - t1.mx * t2.vx -
                     t1.my * t2.vy - t1.mz * t2.vz;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rwdg :: rwdg(trivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg(TriVec3dc<T> const& t,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-t.vx * B.mx - t.vy * B.my - t.vz * B.mz - t.mx * B.vx -
                            t.my * B.vy - t.mz * B.vz - t.px * B.px - t.py * B.py -
                            t.pz * B.pz - t.pw * B.pw);
}

// cga3dc rwdg :: rwdg(bivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg(BiVec3dc<T> const& B,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-B.vx * t.mx - B.vy * t.my - B.vz * t.mz - B.mx * t.vx -
                            B.my * t.vy - B.mz * t.vz - B.px * t.px - B.py * t.py -
                            B.pz * t.pz - B.pw * t.pw);
}

// cga3dc rwdg :: rwdg(trivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] TriVec3dc<T> const&,
                                                   [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(vec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dc<T> const&,
                                                   [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] TriVec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(bivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec3dc<T> const&,
                                                   [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(bivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec3dc<T> const&,
                                                   [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(vec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dc<T> const&,
                                                   [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(vec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dc<T> const&,
                                                   [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rwdg :: rwdg(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The result lies in B and is perpendicular to A. It is the left interior
// (inner) product: the grade-lowering part <A B>_(grade(B)-grade(A)) of the
// geometric product. On equal grades it reduces to the metric inner product:
//
//     (A << B) == dot(A, B)
//
// Equivalently, as defined in "PGA Illuminated" (E. Lengyel):
//
//     operator<<(a,b) = l_contract(a,b) = rwdg( dual(a), b )
//
// (holds exactly at every grade; the cga3dc metric is non-degenerate, so there
// is a single contraction pair, no bulk/weight split).
////////////////////////////////////////////////////////////////////////////////

// cga3dc l_contract :: l_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& A,
                                                       MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c12 + A.c7 * B.c13 + A.c8 * B.c14 +
                     A.c9 * B.c9 + A.c10 * B.c10 + A.c11 * B.c11 + A.c12 * B.c6 +
                     A.c13 * B.c7 + A.c14 * B.c8 - A.c15 * B.c15 - A.c16 * B.c16 -
                     A.c17 * B.c17 - A.c18 * B.c18 - A.c19 * B.c22 - A.c20 * B.c23 -
                     A.c21 * B.c24 - A.c22 * B.c19 - A.c23 * B.c20 - A.c24 * B.c21 +
                     A.c25 * B.c25 - A.c26 * B.c26 - A.c27 * B.c27 - A.c28 * B.c28 +
                     A.c29 * B.c30 + A.c30 * B.c29 - A.c31 * B.c31;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c11 - A.c3 * B.c10 - A.c4 * B.c12 +
                     A.c5 * B.c6 + A.c7 * B.c21 - A.c8 * B.c20 - A.c9 * B.c25 -
                     A.c13 * B.c24 + A.c14 * B.c23 + A.c15 * B.c16 + A.c17 * B.c28 -
                     A.c18 * B.c27 - A.c19 * B.c30 + A.c22 * B.c29 - A.c26 * B.c31;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c11 + A.c3 * B.c9 - A.c4 * B.c13 +
                     A.c5 * B.c7 - A.c6 * B.c21 + A.c8 * B.c19 - A.c10 * B.c25 +
                     A.c12 * B.c24 - A.c14 * B.c22 + A.c15 * B.c17 - A.c16 * B.c28 +
                     A.c18 * B.c26 - A.c20 * B.c30 + A.c23 * B.c29 - A.c27 * B.c31;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c10 - A.c2 * B.c9 - A.c4 * B.c14 +
                     A.c5 * B.c8 + A.c6 * B.c20 - A.c7 * B.c19 - A.c11 * B.c25 -
                     A.c12 * B.c23 + A.c13 * B.c22 + A.c15 * B.c18 + A.c16 * B.c27 -
                     A.c17 * B.c26 - A.c21 * B.c30 + A.c24 * B.c29 - A.c28 * B.c31;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c6 + A.c2 * B.c7 + A.c3 * B.c8 -
                     A.c4 * B.c15 + A.c6 * B.c16 + A.c7 * B.c17 + A.c8 * B.c18 +
                     A.c9 * B.c22 + A.c10 * B.c23 + A.c11 * B.c24 - A.c22 * B.c26 -
                     A.c23 * B.c27 - A.c24 * B.c28 + A.c25 * B.c30 + A.c30 * B.c31;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c12 - A.c2 * B.c13 - A.c3 * B.c14 +
                     A.c5 * B.c15 + A.c9 * B.c19 + A.c10 * B.c20 + A.c11 * B.c21 +
                     A.c12 * B.c16 + A.c13 * B.c17 + A.c14 * B.c18 + A.c19 * B.c26 +
                     A.c20 * B.c27 + A.c21 * B.c28 - A.c25 * B.c29 + A.c29 * B.c31;
    ctype const c6 = A.c0 * B.c6 + A.c2 * B.c24 - A.c3 * B.c23 - A.c4 * B.c16 +
                     A.c7 * B.c28 - A.c8 * B.c27 - A.c9 * B.c30 + A.c22 * B.c31;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c24 + A.c3 * B.c22 - A.c4 * B.c17 -
                     A.c6 * B.c28 + A.c8 * B.c26 - A.c10 * B.c30 + A.c23 * B.c31;
    ctype const c8 = A.c0 * B.c8 + A.c1 * B.c23 - A.c2 * B.c22 - A.c4 * B.c18 +
                     A.c6 * B.c27 - A.c7 * B.c26 - A.c11 * B.c30 + A.c24 * B.c31;
    ctype const c9 = A.c0 * B.c9 - A.c1 * B.c25 - A.c4 * B.c19 - A.c5 * B.c22 -
                     A.c6 * B.c29 - A.c12 * B.c30 - A.c15 * B.c26 + A.c16 * B.c31;
    ctype const c10 = A.c0 * B.c10 - A.c2 * B.c25 - A.c4 * B.c20 - A.c5 * B.c23 -
                      A.c7 * B.c29 - A.c13 * B.c30 - A.c15 * B.c27 + A.c17 * B.c31;
    ctype const c11 = A.c0 * B.c11 - A.c3 * B.c25 - A.c4 * B.c21 - A.c5 * B.c24 -
                      A.c8 * B.c29 - A.c14 * B.c30 - A.c15 * B.c28 + A.c18 * B.c31;
    ctype const c12 = A.c0 * B.c12 - A.c2 * B.c21 + A.c3 * B.c20 - A.c5 * B.c16 -
                      A.c9 * B.c29 - A.c13 * B.c28 + A.c14 * B.c27 + A.c19 * B.c31;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c21 - A.c3 * B.c19 - A.c5 * B.c17 -
                      A.c10 * B.c29 + A.c12 * B.c28 - A.c14 * B.c26 + A.c20 * B.c31;
    ctype const c14 = A.c0 * B.c14 - A.c1 * B.c20 + A.c2 * B.c19 - A.c5 * B.c18 -
                      A.c11 * B.c29 - A.c12 * B.c27 + A.c13 * B.c26 + A.c21 * B.c31;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c16 - A.c2 * B.c17 - A.c3 * B.c18 +
                      A.c9 * B.c26 + A.c10 * B.c27 + A.c11 * B.c28 - A.c25 * B.c31;
    ctype const c16 = A.c0 * B.c16 - A.c2 * B.c28 + A.c3 * B.c27 - A.c9 * B.c31;
    ctype const c17 = A.c0 * B.c17 + A.c1 * B.c28 - A.c3 * B.c26 - A.c10 * B.c31;
    ctype const c18 = A.c0 * B.c18 - A.c1 * B.c27 + A.c2 * B.c26 - A.c11 * B.c31;
    ctype const c19 = A.c0 * B.c19 + A.c1 * B.c29 + A.c5 * B.c26 - A.c12 * B.c31;
    ctype const c20 = A.c0 * B.c20 + A.c2 * B.c29 + A.c5 * B.c27 - A.c13 * B.c31;
    ctype const c21 = A.c0 * B.c21 + A.c3 * B.c29 + A.c5 * B.c28 - A.c14 * B.c31;
    ctype const c22 = A.c0 * B.c22 - A.c1 * B.c30 - A.c4 * B.c26 - A.c6 * B.c31;
    ctype const c23 = A.c0 * B.c23 - A.c2 * B.c30 - A.c4 * B.c27 - A.c7 * B.c31;
    ctype const c24 = A.c0 * B.c24 - A.c3 * B.c30 - A.c4 * B.c28 - A.c8 * B.c31;
    ctype const c25 = A.c0 * B.c25 - A.c4 * B.c29 + A.c5 * B.c30 + A.c15 * B.c31;
    ctype const c26 = A.c0 * B.c26 + A.c1 * B.c31;
    ctype const c27 = A.c0 * B.c27 + A.c2 * B.c31;
    ctype const c28 = A.c0 * B.c28 + A.c3 * B.c31;
    ctype const c29 = A.c0 * B.c29 - A.c5 * B.c31;
    ctype const c30 = A.c0 * B.c30 - A.c4 * B.c31;
    ctype const c31 = A.c0 * B.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& A,
                                                       MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c6 * B.c7 + A.c7 * B.c8 + A.c8 * B.c9 + A.c9 * B.c4 +
                     A.c10 * B.c5 + A.c11 * B.c6 + A.c12 * B.c1 + A.c13 * B.c2 +
                     A.c14 * B.c3 - A.c15 * B.c10 - A.c26 * B.c11 - A.c27 * B.c12 -
                     A.c28 * B.c13 + A.c29 * B.c15 + A.c30 * B.c14;
    ctype const c1 = A.c2 * B.c6 - A.c3 * B.c5 - A.c4 * B.c7 + A.c5 * B.c1 +
                     A.c17 * B.c13 - A.c18 * B.c12 - A.c19 * B.c15 + A.c22 * B.c14;
    ctype const c2 = -A.c1 * B.c6 + A.c3 * B.c4 - A.c4 * B.c8 + A.c5 * B.c2 -
                     A.c16 * B.c13 + A.c18 * B.c11 - A.c20 * B.c15 + A.c23 * B.c14;
    ctype const c3 = A.c1 * B.c5 - A.c2 * B.c4 - A.c4 * B.c9 + A.c5 * B.c3 +
                     A.c16 * B.c12 - A.c17 * B.c11 - A.c21 * B.c15 + A.c24 * B.c14;
    ctype const c4 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c10 -
                     A.c22 * B.c11 - A.c23 * B.c12 - A.c24 * B.c13 + A.c25 * B.c15;
    ctype const c5 = -A.c1 * B.c7 - A.c2 * B.c8 - A.c3 * B.c9 + A.c5 * B.c10 +
                     A.c19 * B.c11 + A.c20 * B.c12 + A.c21 * B.c13 - A.c25 * B.c14;
    ctype const c6 = A.c0 * B.c1 + A.c7 * B.c13 - A.c8 * B.c12 - A.c9 * B.c15;
    ctype const c7 = A.c0 * B.c2 - A.c6 * B.c13 + A.c8 * B.c11 - A.c10 * B.c15;
    ctype const c8 = A.c0 * B.c3 + A.c6 * B.c12 - A.c7 * B.c11 - A.c11 * B.c15;
    ctype const c9 = A.c0 * B.c4 - A.c6 * B.c14 - A.c12 * B.c15 - A.c15 * B.c11;
    ctype const c10 = A.c0 * B.c5 - A.c7 * B.c14 - A.c13 * B.c15 - A.c15 * B.c12;
    ctype const c11 = A.c0 * B.c6 - A.c8 * B.c14 - A.c14 * B.c15 - A.c15 * B.c13;
    ctype const c12 = A.c0 * B.c7 - A.c9 * B.c14 - A.c13 * B.c13 + A.c14 * B.c12;
    ctype const c13 = A.c0 * B.c8 - A.c10 * B.c14 + A.c12 * B.c13 - A.c14 * B.c11;
    ctype const c14 = A.c0 * B.c9 - A.c11 * B.c14 - A.c12 * B.c12 + A.c13 * B.c11;
    ctype const c15 = A.c0 * B.c10 + A.c9 * B.c11 + A.c10 * B.c12 + A.c11 * B.c13;
    ctype const c16 = -A.c2 * B.c13 + A.c3 * B.c12;
    ctype const c17 = A.c1 * B.c13 - A.c3 * B.c11;
    ctype const c18 = -A.c1 * B.c12 + A.c2 * B.c11;
    ctype const c19 = A.c1 * B.c14 + A.c5 * B.c11;
    ctype const c20 = A.c2 * B.c14 + A.c5 * B.c12;
    ctype const c21 = A.c3 * B.c14 + A.c5 * B.c13;
    ctype const c22 = -A.c1 * B.c15 - A.c4 * B.c11;
    ctype const c23 = -A.c2 * B.c15 - A.c4 * B.c12;
    ctype const c24 = -A.c3 * B.c15 - A.c4 * B.c13;
    ctype const c25 = -A.c4 * B.c14 + A.c5 * B.c15;
    ctype const c26 = A.c0 * B.c11;
    ctype const c27 = A.c0 * B.c12;
    ctype const c28 = A.c0 * B.c13;
    ctype const c29 = A.c0 * B.c14;
    ctype const c30 = A.c0 * B.c15;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& A,
                                                       MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                     A.c4 * B.c9 + A.c5 * B.c10 + A.c6 * B.c11 + A.c7 * B.c6 +
                     A.c8 * B.c7 + A.c9 * B.c8 - A.c10 * B.c15 - A.c11 * B.c26 -
                     A.c12 * B.c27 - A.c13 * B.c28 + A.c14 * B.c30 + A.c15 * B.c29;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c21 - A.c3 * B.c20 - A.c4 * B.c25 -
                     A.c8 * B.c24 + A.c9 * B.c23 + A.c10 * B.c16 - A.c11 * B.c31;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c21 + A.c3 * B.c19 - A.c5 * B.c25 +
                     A.c7 * B.c24 - A.c9 * B.c22 + A.c10 * B.c17 - A.c12 * B.c31;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c20 - A.c2 * B.c19 - A.c6 * B.c25 -
                     A.c7 * B.c23 + A.c8 * B.c22 + A.c10 * B.c18 - A.c13 * B.c31;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c16 + A.c2 * B.c17 + A.c3 * B.c18 +
                     A.c4 * B.c22 + A.c5 * B.c23 + A.c6 * B.c24 + A.c15 * B.c31;
    ctype const c5 = A.c0 * B.c5 + A.c4 * B.c19 + A.c5 * B.c20 + A.c6 * B.c21 +
                     A.c7 * B.c16 + A.c8 * B.c17 + A.c9 * B.c18 + A.c14 * B.c31;
    ctype const c6 = A.c0 * B.c6 + A.c2 * B.c28 - A.c3 * B.c27 - A.c4 * B.c30;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c28 + A.c3 * B.c26 - A.c5 * B.c30;
    ctype const c8 = A.c0 * B.c8 + A.c1 * B.c27 - A.c2 * B.c26 - A.c6 * B.c30;
    ctype const c9 = A.c0 * B.c9 - A.c1 * B.c29 - A.c7 * B.c30 - A.c10 * B.c26;
    ctype const c10 = A.c0 * B.c10 - A.c2 * B.c29 - A.c8 * B.c30 - A.c10 * B.c27;
    ctype const c11 = A.c0 * B.c11 - A.c3 * B.c29 - A.c9 * B.c30 - A.c10 * B.c28;
    ctype const c12 = A.c0 * B.c12 - A.c4 * B.c29 - A.c8 * B.c28 + A.c9 * B.c27;
    ctype const c13 = A.c0 * B.c13 - A.c5 * B.c29 + A.c7 * B.c28 - A.c9 * B.c26;
    ctype const c14 = A.c0 * B.c14 - A.c6 * B.c29 - A.c7 * B.c27 + A.c8 * B.c26;
    ctype const c15 = A.c0 * B.c15 + A.c4 * B.c26 + A.c5 * B.c27 + A.c6 * B.c28;
    ctype const c16 = A.c0 * B.c16 - A.c4 * B.c31;
    ctype const c17 = A.c0 * B.c17 - A.c5 * B.c31;
    ctype const c18 = A.c0 * B.c18 - A.c6 * B.c31;
    ctype const c19 = A.c0 * B.c19 - A.c7 * B.c31;
    ctype const c20 = A.c0 * B.c20 - A.c8 * B.c31;
    ctype const c21 = A.c0 * B.c21 - A.c9 * B.c31;
    ctype const c22 = A.c0 * B.c22 - A.c1 * B.c31;
    ctype const c23 = A.c0 * B.c23 - A.c2 * B.c31;
    ctype const c24 = A.c0 * B.c24 - A.c3 * B.c31;
    ctype const c25 = A.c0 * B.c25 + A.c10 * B.c31;
    ctype const c26 = A.c0 * B.c26;
    ctype const c27 = A.c0 * B.c27;
    ctype const c28 = A.c0 * B.c28;
    ctype const c29 = A.c0 * B.c29;
    ctype const c30 = A.c0 * B.c30;
    ctype const c31 = A.c0 * B.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& A,
                                                       MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c4 * B.c4 - A.c5 * B.c3 -
                     A.c16 * B.c5 - A.c17 * B.c6 - A.c18 * B.c7 - A.c19 * B.c11 -
                     A.c20 * B.c12 - A.c21 * B.c13 - A.c22 * B.c8 - A.c23 * B.c9 -
                     A.c24 * B.c10 + A.c25 * B.c14 - A.c31 * B.c15;
    ctype const c1 = A.c0 * B.c0 + A.c7 * B.c10 - A.c8 * B.c9 - A.c9 * B.c14 -
                     A.c13 * B.c13 + A.c14 * B.c12 + A.c15 * B.c5 - A.c26 * B.c15;
    ctype const c2 = A.c0 * B.c1 - A.c6 * B.c10 + A.c8 * B.c8 - A.c10 * B.c14 +
                     A.c12 * B.c13 - A.c14 * B.c11 + A.c15 * B.c6 - A.c27 * B.c15;
    ctype const c3 = A.c0 * B.c2 + A.c6 * B.c9 - A.c7 * B.c8 - A.c11 * B.c14 -
                     A.c12 * B.c12 + A.c13 * B.c11 + A.c15 * B.c7 - A.c28 * B.c15;
    ctype const c4 = A.c0 * B.c3 + A.c6 * B.c5 + A.c7 * B.c6 + A.c8 * B.c7 +
                     A.c9 * B.c11 + A.c10 * B.c12 + A.c11 * B.c13 + A.c30 * B.c15;
    ctype const c5 = A.c0 * B.c4 + A.c9 * B.c8 + A.c10 * B.c9 + A.c11 * B.c10 +
                     A.c12 * B.c5 + A.c13 * B.c6 + A.c14 * B.c7 + A.c29 * B.c15;
    ctype const c6 = A.c2 * B.c13 - A.c3 * B.c12 - A.c4 * B.c5 + A.c22 * B.c15;
    ctype const c7 = -A.c1 * B.c13 + A.c3 * B.c11 - A.c4 * B.c6 + A.c23 * B.c15;
    ctype const c8 = A.c1 * B.c12 - A.c2 * B.c11 - A.c4 * B.c7 + A.c24 * B.c15;
    ctype const c9 = -A.c1 * B.c14 - A.c4 * B.c8 - A.c5 * B.c11 + A.c16 * B.c15;
    ctype const c10 = -A.c2 * B.c14 - A.c4 * B.c9 - A.c5 * B.c12 + A.c17 * B.c15;
    ctype const c11 = -A.c3 * B.c14 - A.c4 * B.c10 - A.c5 * B.c13 + A.c18 * B.c15;
    ctype const c12 = -A.c2 * B.c10 + A.c3 * B.c9 - A.c5 * B.c5 + A.c19 * B.c15;
    ctype const c13 = A.c1 * B.c10 - A.c3 * B.c8 - A.c5 * B.c6 + A.c20 * B.c15;
    ctype const c14 = -A.c1 * B.c9 + A.c2 * B.c8 - A.c5 * B.c7 + A.c21 * B.c15;
    ctype const c15 = -A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 - A.c25 * B.c15;
    ctype const c16 = A.c0 * B.c5 - A.c9 * B.c15;
    ctype const c17 = A.c0 * B.c6 - A.c10 * B.c15;
    ctype const c18 = A.c0 * B.c7 - A.c11 * B.c15;
    ctype const c19 = A.c0 * B.c8 - A.c12 * B.c15;
    ctype const c20 = A.c0 * B.c9 - A.c13 * B.c15;
    ctype const c21 = A.c0 * B.c10 - A.c14 * B.c15;
    ctype const c22 = A.c0 * B.c11 - A.c6 * B.c15;
    ctype const c23 = A.c0 * B.c12 - A.c7 * B.c15;
    ctype const c24 = A.c0 * B.c13 - A.c8 * B.c15;
    ctype const c25 = A.c0 * B.c14 + A.c15 * B.c15;
    ctype const c26 = A.c1 * B.c15;
    ctype const c27 = A.c2 * B.c15;
    ctype const c28 = A.c3 * B.c15;
    ctype const c29 = -A.c5 * B.c15;
    ctype const c30 = -A.c4 * B.c15;
    ctype const c31 = A.c0 * B.c15;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& A,
                                                       MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c3 * B.c5 - A.c4 * B.c4 -
                     A.c5 * B.c16 - A.c6 * B.c17 - A.c7 * B.c18 - A.c8 * B.c22 -
                     A.c9 * B.c23 - A.c10 * B.c24 - A.c11 * B.c19 - A.c12 * B.c20 -
                     A.c13 * B.c21 + A.c14 * B.c25 - A.c15 * B.c31;
    ctype const c1 = A.c1 * B.c11 - A.c2 * B.c10 - A.c3 * B.c12 + A.c4 * B.c6 +
                     A.c6 * B.c28 - A.c7 * B.c27 - A.c8 * B.c30 + A.c11 * B.c29;
    ctype const c2 = -A.c0 * B.c11 + A.c2 * B.c9 - A.c3 * B.c13 + A.c4 * B.c7 -
                     A.c5 * B.c28 + A.c7 * B.c26 - A.c9 * B.c30 + A.c12 * B.c29;
    ctype const c3 = A.c0 * B.c10 - A.c1 * B.c9 - A.c3 * B.c14 + A.c4 * B.c8 +
                     A.c5 * B.c27 - A.c6 * B.c26 - A.c10 * B.c30 + A.c13 * B.c29;
    ctype const c4 = A.c0 * B.c6 + A.c1 * B.c7 + A.c2 * B.c8 - A.c3 * B.c15 -
                     A.c11 * B.c26 - A.c12 * B.c27 - A.c13 * B.c28 + A.c14 * B.c30;
    ctype const c5 = -A.c0 * B.c12 - A.c1 * B.c13 - A.c2 * B.c14 + A.c4 * B.c15 +
                     A.c8 * B.c26 + A.c9 * B.c27 + A.c10 * B.c28 - A.c14 * B.c29;
    ctype const c6 = A.c1 * B.c24 - A.c2 * B.c23 - A.c3 * B.c16 + A.c11 * B.c31;
    ctype const c7 = -A.c0 * B.c24 + A.c2 * B.c22 - A.c3 * B.c17 + A.c12 * B.c31;
    ctype const c8 = A.c0 * B.c23 - A.c1 * B.c22 - A.c3 * B.c18 + A.c13 * B.c31;
    ctype const c9 = -A.c0 * B.c25 - A.c3 * B.c19 - A.c4 * B.c22 + A.c5 * B.c31;
    ctype const c10 = -A.c1 * B.c25 - A.c3 * B.c20 - A.c4 * B.c23 + A.c6 * B.c31;
    ctype const c11 = -A.c2 * B.c25 - A.c3 * B.c21 - A.c4 * B.c24 + A.c7 * B.c31;
    ctype const c12 = -A.c1 * B.c21 + A.c2 * B.c20 - A.c4 * B.c16 + A.c8 * B.c31;
    ctype const c13 = A.c0 * B.c21 - A.c2 * B.c19 - A.c4 * B.c17 + A.c9 * B.c31;
    ctype const c14 = -A.c0 * B.c20 + A.c1 * B.c19 - A.c4 * B.c18 + A.c10 * B.c31;
    ctype const c15 = -A.c0 * B.c16 - A.c1 * B.c17 - A.c2 * B.c18 - A.c14 * B.c31;
    ctype const c16 = -A.c1 * B.c28 + A.c2 * B.c27;
    ctype const c17 = A.c0 * B.c28 - A.c2 * B.c26;
    ctype const c18 = -A.c0 * B.c27 + A.c1 * B.c26;
    ctype const c19 = A.c0 * B.c29 + A.c4 * B.c26;
    ctype const c20 = A.c1 * B.c29 + A.c4 * B.c27;
    ctype const c21 = A.c2 * B.c29 + A.c4 * B.c28;
    ctype const c22 = -A.c0 * B.c30 - A.c3 * B.c26;
    ctype const c23 = -A.c1 * B.c30 - A.c3 * B.c27;
    ctype const c24 = -A.c2 * B.c30 - A.c3 * B.c28;
    ctype const c25 = -A.c3 * B.c29 + A.c4 * B.c30;
    ctype const c26 = A.c0 * B.c31;
    ctype const c27 = A.c1 * B.c31;
    ctype const c28 = A.c2 * B.c31;
    ctype const c29 = -A.c4 * B.c31;
    ctype const c30 = -A.c3 * B.c31;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& M,
                                                       PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c31 * ctype(ps);
    ctype const c1 = -M.c26 * ctype(ps);
    ctype const c2 = -M.c27 * ctype(ps);
    ctype const c3 = -M.c28 * ctype(ps);
    ctype const c4 = M.c30 * ctype(ps);
    ctype const c5 = M.c29 * ctype(ps);
    ctype const c6 = M.c22 * ctype(ps);
    ctype const c7 = M.c23 * ctype(ps);
    ctype const c8 = M.c24 * ctype(ps);
    ctype const c9 = M.c16 * ctype(ps);
    ctype const c10 = M.c17 * ctype(ps);
    ctype const c11 = M.c18 * ctype(ps);
    ctype const c12 = M.c19 * ctype(ps);
    ctype const c13 = M.c20 * ctype(ps);
    ctype const c14 = M.c21 * ctype(ps);
    ctype const c15 = -M.c25 * ctype(ps);
    ctype const c16 = -M.c9 * ctype(ps);
    ctype const c17 = -M.c10 * ctype(ps);
    ctype const c18 = -M.c11 * ctype(ps);
    ctype const c19 = -M.c12 * ctype(ps);
    ctype const c20 = -M.c13 * ctype(ps);
    ctype const c21 = -M.c14 * ctype(ps);
    ctype const c22 = -M.c6 * ctype(ps);
    ctype const c23 = -M.c7 * ctype(ps);
    ctype const c24 = -M.c8 * ctype(ps);
    ctype const c25 = M.c15 * ctype(ps);
    ctype const c26 = M.c1 * ctype(ps);
    ctype const c27 = M.c2 * ctype(ps);
    ctype const c28 = M.c3 * ctype(ps);
    ctype const c29 = -M.c5 * ctype(ps);
    ctype const c30 = -M.c4 * ctype(ps);
    ctype const c31 = M.c0 * ctype(ps);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(PScalar3dc<T> ps,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c31;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,quadvec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& M,
                                                       QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c26 * Q.x - M.c27 * Q.y - M.c28 * Q.z + M.c29 * Q.u + M.c30 * Q.w;
    ctype const c1 = M.c17 * Q.z - M.c18 * Q.y - M.c19 * Q.u + M.c22 * Q.w;
    ctype const c2 = -M.c16 * Q.z + M.c18 * Q.x - M.c20 * Q.u + M.c23 * Q.w;
    ctype const c3 = M.c16 * Q.y - M.c17 * Q.x - M.c21 * Q.u + M.c24 * Q.w;
    ctype const c4 = -M.c22 * Q.x - M.c23 * Q.y - M.c24 * Q.z + M.c25 * Q.u;
    ctype const c5 = M.c19 * Q.x + M.c20 * Q.y + M.c21 * Q.z - M.c25 * Q.w;
    ctype const c6 = M.c7 * Q.z - M.c8 * Q.y - M.c9 * Q.u;
    ctype const c7 = -M.c6 * Q.z + M.c8 * Q.x - M.c10 * Q.u;
    ctype const c8 = M.c6 * Q.y - M.c7 * Q.x - M.c11 * Q.u;
    ctype const c9 = -M.c6 * Q.w - M.c12 * Q.u - M.c15 * Q.x;
    ctype const c10 = -M.c7 * Q.w - M.c13 * Q.u - M.c15 * Q.y;
    ctype const c11 = -M.c8 * Q.w - M.c14 * Q.u - M.c15 * Q.z;
    ctype const c12 = -M.c9 * Q.w - M.c13 * Q.z + M.c14 * Q.y;
    ctype const c13 = -M.c10 * Q.w + M.c12 * Q.z - M.c14 * Q.x;
    ctype const c14 = -M.c11 * Q.w - M.c12 * Q.y + M.c13 * Q.x;
    ctype const c15 = M.c9 * Q.x + M.c10 * Q.y + M.c11 * Q.z;
    ctype const c16 = -M.c2 * Q.z + M.c3 * Q.y;
    ctype const c17 = M.c1 * Q.z - M.c3 * Q.x;
    ctype const c18 = -M.c1 * Q.y + M.c2 * Q.x;
    ctype const c19 = M.c1 * Q.w + M.c5 * Q.x;
    ctype const c20 = M.c2 * Q.w + M.c5 * Q.y;
    ctype const c21 = M.c3 * Q.w + M.c5 * Q.z;
    ctype const c22 = -M.c1 * Q.u - M.c4 * Q.x;
    ctype const c23 = -M.c2 * Q.u - M.c4 * Q.y;
    ctype const c24 = -M.c3 * Q.u - M.c4 * Q.z;
    ctype const c25 = -M.c4 * Q.w + M.c5 * Q.u;
    ctype const c26 = M.c0 * Q.x;
    ctype const c27 = M.c0 * Q.y;
    ctype const c28 = M.c0 * Q.z;
    ctype const c29 = M.c0 * Q.w;
    ctype const c30 = M.c0 * Q.u;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(quadvec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(QuadVec3dc<T> const& Q,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c26 - Q.y * M.c27 - Q.z * M.c28 + Q.w * M.c30 + Q.u * M.c29;
    ctype const c1 = -Q.x * M.c31;
    ctype const c2 = -Q.y * M.c31;
    ctype const c3 = -Q.z * M.c31;
    ctype const c4 = Q.u * M.c31;
    ctype const c5 = Q.w * M.c31;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& M,
                                                       TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c16 * t.vx - M.c17 * t.vy - M.c18 * t.vz - M.c19 * t.px -
                     M.c20 * t.py - M.c21 * t.pz - M.c22 * t.mx - M.c23 * t.my -
                     M.c24 * t.mz + M.c25 * t.pw;
    ctype const c1 = M.c7 * t.mz - M.c8 * t.my - M.c9 * t.pw - M.c13 * t.pz +
                     M.c14 * t.py + M.c15 * t.vx;
    ctype const c2 = -M.c6 * t.mz + M.c8 * t.mx - M.c10 * t.pw + M.c12 * t.pz -
                     M.c14 * t.px + M.c15 * t.vy;
    ctype const c3 = M.c6 * t.my - M.c7 * t.mx - M.c11 * t.pw - M.c12 * t.py +
                     M.c13 * t.px + M.c15 * t.vz;
    ctype const c4 = M.c6 * t.vx + M.c7 * t.vy + M.c8 * t.vz + M.c9 * t.px +
                     M.c10 * t.py + M.c11 * t.pz;
    ctype const c5 = M.c9 * t.mx + M.c10 * t.my + M.c11 * t.mz + M.c12 * t.vx +
                     M.c13 * t.vy + M.c14 * t.vz;
    ctype const c6 = M.c2 * t.pz - M.c3 * t.py - M.c4 * t.vx;
    ctype const c7 = -M.c1 * t.pz + M.c3 * t.px - M.c4 * t.vy;
    ctype const c8 = M.c1 * t.py - M.c2 * t.px - M.c4 * t.vz;
    ctype const c9 = -M.c1 * t.pw - M.c4 * t.mx - M.c5 * t.px;
    ctype const c10 = -M.c2 * t.pw - M.c4 * t.my - M.c5 * t.py;
    ctype const c11 = -M.c3 * t.pw - M.c4 * t.mz - M.c5 * t.pz;
    ctype const c12 = -M.c2 * t.mz + M.c3 * t.my - M.c5 * t.vx;
    ctype const c13 = M.c1 * t.mz - M.c3 * t.mx - M.c5 * t.vy;
    ctype const c14 = -M.c1 * t.my + M.c2 * t.mx - M.c5 * t.vz;
    ctype const c15 = -M.c1 * t.vx - M.c2 * t.vy - M.c3 * t.vz;
    ctype const c16 = M.c0 * t.vx;
    ctype const c17 = M.c0 * t.vy;
    ctype const c18 = M.c0 * t.vz;
    ctype const c19 = M.c0 * t.mx;
    ctype const c20 = M.c0 * t.my;
    ctype const c21 = M.c0 * t.mz;
    ctype const c22 = M.c0 * t.px;
    ctype const c23 = M.c0 * t.py;
    ctype const c24 = M.c0 * t.pz;
    ctype const c25 = M.c0 * t.pw;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(TriVec3dc<T> const& t,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c16 - t.vy * M.c17 - t.vz * M.c18 - t.mx * M.c22 -
                     t.my * M.c23 - t.mz * M.c24 - t.px * M.c19 - t.py * M.c20 -
                     t.pz * M.c21 + t.pw * M.c25;
    ctype const c1 = t.vy * M.c28 - t.vz * M.c27 - t.mx * M.c30 + t.px * M.c29;
    ctype const c2 = -t.vx * M.c28 + t.vz * M.c26 - t.my * M.c30 + t.py * M.c29;
    ctype const c3 = t.vx * M.c27 - t.vy * M.c26 - t.mz * M.c30 + t.pz * M.c29;
    ctype const c4 = -t.px * M.c26 - t.py * M.c27 - t.pz * M.c28 + t.pw * M.c30;
    ctype const c5 = t.mx * M.c26 + t.my * M.c27 + t.mz * M.c28 - t.pw * M.c29;
    ctype const c6 = t.px * M.c31;
    ctype const c7 = t.py * M.c31;
    ctype const c8 = t.pz * M.c31;
    ctype const c9 = t.vx * M.c31;
    ctype const c10 = t.vy * M.c31;
    ctype const c11 = t.vz * M.c31;
    ctype const c12 = t.mx * M.c31;
    ctype const c13 = t.my * M.c31;
    ctype const c14 = t.mz * M.c31;
    ctype const c15 = -t.pw * M.c31;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& M,
                                                       BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c6 * B.px + M.c7 * B.py + M.c8 * B.pz + M.c9 * B.mx +
                     M.c10 * B.my + M.c11 * B.mz + M.c12 * B.vx + M.c13 * B.vy +
                     M.c14 * B.vz - M.c15 * B.pw;
    ctype const c1 = M.c2 * B.mz - M.c3 * B.my - M.c4 * B.px + M.c5 * B.vx;
    ctype const c2 = -M.c1 * B.mz + M.c3 * B.mx - M.c4 * B.py + M.c5 * B.vy;
    ctype const c3 = M.c1 * B.my - M.c2 * B.mx - M.c4 * B.pz + M.c5 * B.vz;
    ctype const c4 = M.c1 * B.vx + M.c2 * B.vy + M.c3 * B.vz - M.c4 * B.pw;
    ctype const c5 = -M.c1 * B.px - M.c2 * B.py - M.c3 * B.pz + M.c5 * B.pw;
    ctype const c6 = M.c0 * B.vx;
    ctype const c7 = M.c0 * B.vy;
    ctype const c8 = M.c0 * B.vz;
    ctype const c9 = M.c0 * B.mx;
    ctype const c10 = M.c0 * B.my;
    ctype const c11 = M.c0 * B.mz;
    ctype const c12 = M.c0 * B.px;
    ctype const c13 = M.c0 * B.py;
    ctype const c14 = M.c0 * B.pz;
    ctype const c15 = M.c0 * B.pw;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(BiVec3dc<T> const& B,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c12 + B.vy * M.c13 + B.vz * M.c14 + B.mx * M.c9 +
                     B.my * M.c10 + B.mz * M.c11 + B.px * M.c6 + B.py * M.c7 +
                     B.pz * M.c8 - B.pw * M.c15;
    ctype const c1 = B.vy * M.c21 - B.vz * M.c20 - B.mx * M.c25 - B.py * M.c24 +
                     B.pz * M.c23 + B.pw * M.c16;
    ctype const c2 = -B.vx * M.c21 + B.vz * M.c19 - B.my * M.c25 + B.px * M.c24 -
                     B.pz * M.c22 + B.pw * M.c17;
    ctype const c3 = B.vx * M.c20 - B.vy * M.c19 - B.mz * M.c25 - B.px * M.c23 +
                     B.py * M.c22 + B.pw * M.c18;
    ctype const c4 = B.vx * M.c16 + B.vy * M.c17 + B.vz * M.c18 + B.mx * M.c22 +
                     B.my * M.c23 + B.mz * M.c24;
    ctype const c5 = B.mx * M.c19 + B.my * M.c20 + B.mz * M.c21 + B.px * M.c16 +
                     B.py * M.c17 + B.pz * M.c18;
    ctype const c6 = B.vy * M.c28 - B.vz * M.c27 - B.mx * M.c30;
    ctype const c7 = -B.vx * M.c28 + B.vz * M.c26 - B.my * M.c30;
    ctype const c8 = B.vx * M.c27 - B.vy * M.c26 - B.mz * M.c30;
    ctype const c9 = -B.vx * M.c29 - B.px * M.c30 - B.pw * M.c26;
    ctype const c10 = -B.vy * M.c29 - B.py * M.c30 - B.pw * M.c27;
    ctype const c11 = -B.vz * M.c29 - B.pz * M.c30 - B.pw * M.c28;
    ctype const c12 = -B.mx * M.c29 - B.py * M.c28 + B.pz * M.c27;
    ctype const c13 = -B.my * M.c29 + B.px * M.c28 - B.pz * M.c26;
    ctype const c14 = -B.mz * M.c29 - B.px * M.c27 + B.py * M.c26;
    ctype const c15 = B.mx * M.c26 + B.my * M.c27 + B.mz * M.c28;
    ctype const c16 = -B.mx * M.c31;
    ctype const c17 = -B.my * M.c31;
    ctype const c18 = -B.mz * M.c31;
    ctype const c19 = -B.px * M.c31;
    ctype const c20 = -B.py * M.c31;
    ctype const c21 = -B.pz * M.c31;
    ctype const c22 = -B.vx * M.c31;
    ctype const c23 = -B.vy * M.c31;
    ctype const c24 = -B.vz * M.c31;
    ctype const c25 = B.pw * M.c31;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& M,
                                                       Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z - M.c4 * v.u - M.c5 * v.w;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = M.c0 * v.w;
    ctype const c5 = M.c0 * v.u;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3 - v.w * M.c5 - v.u * M.c4;
    ctype const c1 = v.y * M.c11 - v.z * M.c10 - v.w * M.c12 + v.u * M.c6;
    ctype const c2 = -v.x * M.c11 + v.z * M.c9 - v.w * M.c13 + v.u * M.c7;
    ctype const c3 = v.x * M.c10 - v.y * M.c9 - v.w * M.c14 + v.u * M.c8;
    ctype const c4 = v.x * M.c6 + v.y * M.c7 + v.z * M.c8 - v.w * M.c15;
    ctype const c5 = -v.x * M.c12 - v.y * M.c13 - v.z * M.c14 + v.u * M.c15;
    ctype const c6 = v.y * M.c24 - v.z * M.c23 - v.w * M.c16;
    ctype const c7 = -v.x * M.c24 + v.z * M.c22 - v.w * M.c17;
    ctype const c8 = v.x * M.c23 - v.y * M.c22 - v.w * M.c18;
    ctype const c9 = -v.x * M.c25 - v.w * M.c19 - v.u * M.c22;
    ctype const c10 = -v.y * M.c25 - v.w * M.c20 - v.u * M.c23;
    ctype const c11 = -v.z * M.c25 - v.w * M.c21 - v.u * M.c24;
    ctype const c12 = -v.y * M.c21 + v.z * M.c20 - v.u * M.c16;
    ctype const c13 = v.x * M.c21 - v.z * M.c19 - v.u * M.c17;
    ctype const c14 = -v.x * M.c20 + v.y * M.c19 - v.u * M.c18;
    ctype const c15 = -v.x * M.c16 - v.y * M.c17 - v.z * M.c18;
    ctype const c16 = -v.y * M.c28 + v.z * M.c27;
    ctype const c17 = v.x * M.c28 - v.z * M.c26;
    ctype const c18 = -v.x * M.c27 + v.y * M.c26;
    ctype const c19 = v.x * M.c29 + v.u * M.c26;
    ctype const c20 = v.y * M.c29 + v.u * M.c27;
    ctype const c21 = v.z * M.c29 + v.u * M.c28;
    ctype const c22 = -v.x * M.c30 - v.w * M.c26;
    ctype const c23 = -v.y * M.c30 - v.w * M.c27;
    ctype const c24 = -v.z * M.c30 - v.w * M.c28;
    ctype const c25 = -v.w * M.c29 + v.u * M.c30;
    ctype const c26 = v.x * M.c31;
    ctype const c27 = v.y * M.c31;
    ctype const c28 = v.z * M.c31;
    ctype const c29 = -v.u * M.c31;
    ctype const c30 = -v.w * M.c31;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(MVec3dc<T> const& M,
                                                       Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    ctype const c16 = ctype(s) * M.c16;
    ctype const c17 = ctype(s) * M.c17;
    ctype const c18 = ctype(s) * M.c18;
    ctype const c19 = ctype(s) * M.c19;
    ctype const c20 = ctype(s) * M.c20;
    ctype const c21 = ctype(s) * M.c21;
    ctype const c22 = ctype(s) * M.c22;
    ctype const c23 = ctype(s) * M.c23;
    ctype const c24 = ctype(s) * M.c24;
    ctype const c25 = ctype(s) * M.c25;
    ctype const c26 = ctype(s) * M.c26;
    ctype const c27 = ctype(s) * M.c27;
    ctype const c28 = ctype(s) * M.c28;
    ctype const c29 = ctype(s) * M.c29;
    ctype const c30 = ctype(s) * M.c30;
    ctype const c31 = ctype(s) * M.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc l_contract :: l_contract(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& A,
                                                         MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c7 + A.c2 * B.c8 + A.c3 * B.c9 + A.c4 * B.c4 +
                     A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c1 + A.c8 * B.c2 + A.c9 * B.c3 -
                     A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c12 - A.c13 * B.c13 +
                     A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c13 - A.c3 * B.c12 - A.c4 * B.c15;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c13 + A.c3 * B.c11 - A.c5 * B.c15;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c12 - A.c2 * B.c11 - A.c6 * B.c15;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c14 - A.c7 * B.c15 - A.c10 * B.c11;
    ctype const c5 = A.c0 * B.c5 - A.c2 * B.c14 - A.c8 * B.c15 - A.c10 * B.c12;
    ctype const c6 = A.c0 * B.c6 - A.c3 * B.c14 - A.c9 * B.c15 - A.c10 * B.c13;
    ctype const c7 = A.c0 * B.c7 - A.c4 * B.c14 - A.c8 * B.c13 + A.c9 * B.c12;
    ctype const c8 = A.c0 * B.c8 - A.c5 * B.c14 + A.c7 * B.c13 - A.c9 * B.c11;
    ctype const c9 = A.c0 * B.c9 - A.c6 * B.c14 - A.c7 * B.c12 + A.c8 * B.c11;
    ctype const c10 = A.c0 * B.c10 + A.c4 * B.c11 + A.c5 * B.c12 + A.c6 * B.c13;
    ctype const c11 = A.c0 * B.c11;
    ctype const c12 = A.c0 * B.c12;
    ctype const c13 = A.c0 * B.c13;
    ctype const c14 = A.c0 * B.c14;
    ctype const c15 = A.c0 * B.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& A,
                                                         MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c2 * B.c10 - A.c3 * B.c9 - A.c4 * B.c14 -
                     A.c8 * B.c13 + A.c9 * B.c12 + A.c10 * B.c5 - A.c11 * B.c15;
    ctype const c1 = A.c0 * B.c1 - A.c1 * B.c10 + A.c3 * B.c8 - A.c5 * B.c14 +
                     A.c7 * B.c13 - A.c9 * B.c11 + A.c10 * B.c6 - A.c12 * B.c15;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c9 - A.c2 * B.c8 - A.c6 * B.c14 -
                     A.c7 * B.c12 + A.c8 * B.c11 + A.c10 * B.c7 - A.c13 * B.c15;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 +
                     A.c4 * B.c11 + A.c5 * B.c12 + A.c6 * B.c13 + A.c15 * B.c15;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c8 + A.c5 * B.c9 + A.c6 * B.c10 +
                     A.c7 * B.c5 + A.c8 * B.c6 + A.c9 * B.c7 + A.c14 * B.c15;
    ctype const c5 = A.c0 * B.c5 - A.c4 * B.c15;
    ctype const c6 = A.c0 * B.c6 - A.c5 * B.c15;
    ctype const c7 = A.c0 * B.c7 - A.c6 * B.c15;
    ctype const c8 = A.c0 * B.c8 - A.c7 * B.c15;
    ctype const c9 = A.c0 * B.c9 - A.c8 * B.c15;
    ctype const c10 = A.c0 * B.c10 - A.c9 * B.c15;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c15;
    ctype const c12 = A.c0 * B.c12 - A.c2 * B.c15;
    ctype const c13 = A.c0 * B.c13 - A.c3 * B.c15;
    ctype const c14 = A.c0 * B.c14 + A.c10 * B.c15;
    ctype const c15 = A.c0 * B.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& A,
                                                         MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c7 + A.c4 * B.c1 +
                     A.c6 * B.c13 - A.c7 * B.c12 - A.c8 * B.c15 + A.c11 * B.c14;
    ctype const c1 = -A.c0 * B.c6 + A.c2 * B.c4 - A.c3 * B.c8 + A.c4 * B.c2 -
                     A.c5 * B.c13 + A.c7 * B.c11 - A.c9 * B.c15 + A.c12 * B.c14;
    ctype const c2 = A.c0 * B.c5 - A.c1 * B.c4 - A.c3 * B.c9 + A.c4 * B.c3 +
                     A.c5 * B.c12 - A.c6 * B.c11 - A.c10 * B.c15 + A.c13 * B.c14;
    ctype const c3 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c3 * B.c10 -
                     A.c11 * B.c11 - A.c12 * B.c12 - A.c13 * B.c13 + A.c14 * B.c15;
    ctype const c4 = -A.c0 * B.c7 - A.c1 * B.c8 - A.c2 * B.c9 + A.c4 * B.c10 +
                     A.c8 * B.c11 + A.c9 * B.c12 + A.c10 * B.c13 - A.c14 * B.c14;
    ctype const c5 = -A.c1 * B.c13 + A.c2 * B.c12;
    ctype const c6 = A.c0 * B.c13 - A.c2 * B.c11;
    ctype const c7 = -A.c0 * B.c12 + A.c1 * B.c11;
    ctype const c8 = A.c0 * B.c14 + A.c4 * B.c11;
    ctype const c9 = A.c1 * B.c14 + A.c4 * B.c12;
    ctype const c10 = A.c2 * B.c14 + A.c4 * B.c13;
    ctype const c11 = -A.c0 * B.c15 - A.c3 * B.c11;
    ctype const c12 = -A.c1 * B.c15 - A.c3 * B.c12;
    ctype const c13 = -A.c2 * B.c15 - A.c3 * B.c13;
    ctype const c14 = -A.c3 * B.c14 + A.c4 * B.c15;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(mv_e,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& M,
                                                         PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * ctype(ps);
    ctype const c1 = -M.c12 * ctype(ps);
    ctype const c2 = -M.c13 * ctype(ps);
    ctype const c3 = M.c15 * ctype(ps);
    ctype const c4 = M.c14 * ctype(ps);
    ctype const c5 = -M.c4 * ctype(ps);
    ctype const c6 = -M.c5 * ctype(ps);
    ctype const c7 = -M.c6 * ctype(ps);
    ctype const c8 = -M.c7 * ctype(ps);
    ctype const c9 = -M.c8 * ctype(ps);
    ctype const c10 = -M.c9 * ctype(ps);
    ctype const c11 = -M.c1 * ctype(ps);
    ctype const c12 = -M.c2 * ctype(ps);
    ctype const c13 = -M.c3 * ctype(ps);
    ctype const c14 = M.c10 * ctype(ps);
    ctype const c15 = M.c0 * ctype(ps);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(ps,mv_e) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dc<T>, [[maybe_unused]] MVec3dc_E<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(mv_e,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& M,
                                                         QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * Q.x - M.c12 * Q.y - M.c13 * Q.z + M.c14 * Q.u + M.c15 * Q.w;
    ctype const c1 = M.c2 * Q.z - M.c3 * Q.y - M.c4 * Q.u;
    ctype const c2 = -M.c1 * Q.z + M.c3 * Q.x - M.c5 * Q.u;
    ctype const c3 = M.c1 * Q.y - M.c2 * Q.x - M.c6 * Q.u;
    ctype const c4 = -M.c1 * Q.w - M.c7 * Q.u - M.c10 * Q.x;
    ctype const c5 = -M.c2 * Q.w - M.c8 * Q.u - M.c10 * Q.y;
    ctype const c6 = -M.c3 * Q.w - M.c9 * Q.u - M.c10 * Q.z;
    ctype const c7 = -M.c4 * Q.w - M.c8 * Q.z + M.c9 * Q.y;
    ctype const c8 = -M.c5 * Q.w + M.c7 * Q.z - M.c9 * Q.x;
    ctype const c9 = -M.c6 * Q.w - M.c7 * Q.y + M.c8 * Q.x;
    ctype const c10 = M.c4 * Q.x + M.c5 * Q.y + M.c6 * Q.z;
    ctype const c11 = M.c0 * Q.x;
    ctype const c12 = M.c0 * Q.y;
    ctype const c13 = M.c0 * Q.z;
    ctype const c14 = M.c0 * Q.w;
    ctype const c15 = M.c0 * Q.u;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(quadvec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(QuadVec3dc<T> const& Q,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c11 - Q.y * M.c12 - Q.z * M.c13 + Q.w * M.c15 + Q.u * M.c14;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& M,
                                                         TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * t.mz - M.c3 * t.my - M.c4 * t.pw - M.c8 * t.pz + M.c9 * t.py +
                     M.c10 * t.vx;
    ctype const c1 = -M.c1 * t.mz + M.c3 * t.mx - M.c5 * t.pw + M.c7 * t.pz -
                     M.c9 * t.px + M.c10 * t.vy;
    ctype const c2 = M.c1 * t.my - M.c2 * t.mx - M.c6 * t.pw - M.c7 * t.py + M.c8 * t.px +
                     M.c10 * t.vz;
    ctype const c3 =
        M.c1 * t.vx + M.c2 * t.vy + M.c3 * t.vz + M.c4 * t.px + M.c5 * t.py + M.c6 * t.pz;
    ctype const c4 =
        M.c4 * t.mx + M.c5 * t.my + M.c6 * t.mz + M.c7 * t.vx + M.c8 * t.vy + M.c9 * t.vz;
    ctype const c5 = M.c0 * t.vx;
    ctype const c6 = M.c0 * t.vy;
    ctype const c7 = M.c0 * t.vz;
    ctype const c8 = M.c0 * t.mx;
    ctype const c9 = M.c0 * t.my;
    ctype const c10 = M.c0 * t.mz;
    ctype const c11 = M.c0 * t.px;
    ctype const c12 = M.c0 * t.py;
    ctype const c13 = M.c0 * t.pz;
    ctype const c14 = M.c0 * t.pw;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(TriVec3dc<T> const& t,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vy * M.c13 - t.vz * M.c12 - t.mx * M.c15 + t.px * M.c14;
    ctype const c1 = -t.vx * M.c13 + t.vz * M.c11 - t.my * M.c15 + t.py * M.c14;
    ctype const c2 = t.vx * M.c12 - t.vy * M.c11 - t.mz * M.c15 + t.pz * M.c14;
    ctype const c3 = -t.px * M.c11 - t.py * M.c12 - t.pz * M.c13 + t.pw * M.c15;
    ctype const c4 = t.mx * M.c11 + t.my * M.c12 + t.mz * M.c13 - t.pw * M.c14;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& M,
                                                         BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.px + M.c2 * B.py + M.c3 * B.pz + M.c4 * B.mx + M.c5 * B.my +
                     M.c6 * B.mz + M.c7 * B.vx + M.c8 * B.vy + M.c9 * B.vz - M.c10 * B.pw;
    ctype const c1 = M.c0 * B.vx;
    ctype const c2 = M.c0 * B.vy;
    ctype const c3 = M.c0 * B.vz;
    ctype const c4 = M.c0 * B.mx;
    ctype const c5 = M.c0 * B.my;
    ctype const c6 = M.c0 * B.mz;
    ctype const c7 = M.c0 * B.px;
    ctype const c8 = M.c0 * B.py;
    ctype const c9 = M.c0 * B.pz;
    ctype const c10 = M.c0 * B.pw;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(BiVec3dc<T> const& B,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c7 + B.vy * M.c8 + B.vz * M.c9 + B.mx * M.c4 + B.my * M.c5 +
                     B.mz * M.c6 + B.px * M.c1 + B.py * M.c2 + B.pz * M.c3 - B.pw * M.c10;
    ctype const c1 = B.vy * M.c13 - B.vz * M.c12 - B.mx * M.c15;
    ctype const c2 = -B.vx * M.c13 + B.vz * M.c11 - B.my * M.c15;
    ctype const c3 = B.vx * M.c12 - B.vy * M.c11 - B.mz * M.c15;
    ctype const c4 = -B.vx * M.c14 - B.px * M.c15 - B.pw * M.c11;
    ctype const c5 = -B.vy * M.c14 - B.py * M.c15 - B.pw * M.c12;
    ctype const c6 = -B.vz * M.c14 - B.pz * M.c15 - B.pw * M.c13;
    ctype const c7 = -B.mx * M.c14 - B.py * M.c13 + B.pz * M.c12;
    ctype const c8 = -B.my * M.c14 + B.px * M.c13 - B.pz * M.c11;
    ctype const c9 = -B.mz * M.c14 - B.px * M.c12 + B.py * M.c11;
    ctype const c10 = B.mx * M.c11 + B.my * M.c12 + B.mz * M.c13;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& M,
                                                         Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    ctype const c3 = M.c0 * v.w;
    ctype const c4 = M.c0 * v.u;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * M.c6 - v.z * M.c5 - v.w * M.c7 + v.u * M.c1;
    ctype const c1 = -v.x * M.c6 + v.z * M.c4 - v.w * M.c8 + v.u * M.c2;
    ctype const c2 = v.x * M.c5 - v.y * M.c4 - v.w * M.c9 + v.u * M.c3;
    ctype const c3 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3 - v.w * M.c10;
    ctype const c4 = -v.x * M.c7 - v.y * M.c8 - v.z * M.c9 + v.u * M.c10;
    ctype const c5 = -v.y * M.c13 + v.z * M.c12;
    ctype const c6 = v.x * M.c13 - v.z * M.c11;
    ctype const c7 = -v.x * M.c12 + v.y * M.c11;
    ctype const c8 = v.x * M.c14 + v.u * M.c11;
    ctype const c9 = v.y * M.c14 + v.u * M.c12;
    ctype const c10 = v.z * M.c14 + v.u * M.c13;
    ctype const c11 = -v.x * M.c15 - v.w * M.c11;
    ctype const c12 = -v.y * M.c15 - v.w * M.c12;
    ctype const c13 = -v.z * M.c15 - v.w * M.c13;
    ctype const c14 = -v.w * M.c14 + v.u * M.c15;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_E<T> const& M,
                                                         Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& A,
                                                         MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 - A.c4 * B.c3 -
                     A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7 - A.c8 * B.c11 -
                     A.c9 * B.c12 - A.c10 * B.c13 - A.c11 * B.c8 - A.c12 * B.c9 -
                     A.c13 * B.c10 + A.c14 * B.c14 - A.c15 * B.c15;
    ctype const c1 = A.c1 * B.c13 - A.c2 * B.c12 - A.c3 * B.c5 + A.c11 * B.c15;
    ctype const c2 = -A.c0 * B.c13 + A.c2 * B.c11 - A.c3 * B.c6 + A.c12 * B.c15;
    ctype const c3 = A.c0 * B.c12 - A.c1 * B.c11 - A.c3 * B.c7 + A.c13 * B.c15;
    ctype const c4 = -A.c0 * B.c14 - A.c3 * B.c8 - A.c4 * B.c11 + A.c5 * B.c15;
    ctype const c5 = -A.c1 * B.c14 - A.c3 * B.c9 - A.c4 * B.c12 + A.c6 * B.c15;
    ctype const c6 = -A.c2 * B.c14 - A.c3 * B.c10 - A.c4 * B.c13 + A.c7 * B.c15;
    ctype const c7 = -A.c1 * B.c10 + A.c2 * B.c9 - A.c4 * B.c5 + A.c8 * B.c15;
    ctype const c8 = A.c0 * B.c10 - A.c2 * B.c8 - A.c4 * B.c6 + A.c9 * B.c15;
    ctype const c9 = -A.c0 * B.c9 + A.c1 * B.c8 - A.c4 * B.c7 + A.c10 * B.c15;
    ctype const c10 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c7 - A.c14 * B.c15;
    ctype const c11 = A.c0 * B.c15;
    ctype const c12 = A.c1 * B.c15;
    ctype const c13 = A.c2 * B.c15;
    ctype const c14 = -A.c4 * B.c15;
    ctype const c15 = -A.c3 * B.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_u,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& M,
                                                         PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c15 * ctype(ps);
    ctype const c1 = M.c11 * ctype(ps);
    ctype const c2 = M.c12 * ctype(ps);
    ctype const c3 = M.c13 * ctype(ps);
    ctype const c4 = M.c5 * ctype(ps);
    ctype const c5 = M.c6 * ctype(ps);
    ctype const c6 = M.c7 * ctype(ps);
    ctype const c7 = M.c8 * ctype(ps);
    ctype const c8 = M.c9 * ctype(ps);
    ctype const c9 = M.c10 * ctype(ps);
    ctype const c10 = -M.c14 * ctype(ps);
    ctype const c11 = M.c0 * ctype(ps);
    ctype const c12 = M.c1 * ctype(ps);
    ctype const c13 = M.c2 * ctype(ps);
    ctype const c14 = -M.c4 * ctype(ps);
    ctype const c15 = -M.c3 * ctype(ps);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(ps,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(PScalar3dc<T> ps,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c15;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_u,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& M,
                                                         QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c6 * Q.z - M.c7 * Q.y - M.c8 * Q.u + M.c11 * Q.w;
    ctype const c1 = -M.c5 * Q.z + M.c7 * Q.x - M.c9 * Q.u + M.c12 * Q.w;
    ctype const c2 = M.c5 * Q.y - M.c6 * Q.x - M.c10 * Q.u + M.c13 * Q.w;
    ctype const c3 = -M.c11 * Q.x - M.c12 * Q.y - M.c13 * Q.z + M.c14 * Q.u;
    ctype const c4 = M.c8 * Q.x + M.c9 * Q.y + M.c10 * Q.z - M.c14 * Q.w;
    ctype const c5 = -M.c1 * Q.z + M.c2 * Q.y;
    ctype const c6 = M.c0 * Q.z - M.c2 * Q.x;
    ctype const c7 = -M.c0 * Q.y + M.c1 * Q.x;
    ctype const c8 = M.c0 * Q.w + M.c4 * Q.x;
    ctype const c9 = M.c1 * Q.w + M.c4 * Q.y;
    ctype const c10 = M.c2 * Q.w + M.c4 * Q.z;
    ctype const c11 = -M.c0 * Q.u - M.c3 * Q.x;
    ctype const c12 = -M.c1 * Q.u - M.c3 * Q.y;
    ctype const c13 = -M.c2 * Q.u - M.c3 * Q.z;
    ctype const c14 = -M.c3 * Q.w + M.c4 * Q.u;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(quadvec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(QuadVec3dc<T> const& Q,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c15;
    ctype const c1 = -Q.y * M.c15;
    ctype const c2 = -Q.z * M.c15;
    ctype const c3 = Q.u * M.c15;
    ctype const c4 = Q.w * M.c15;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& M,
                                                         TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * t.vx - M.c6 * t.vy - M.c7 * t.vz - M.c8 * t.px -
                     M.c9 * t.py - M.c10 * t.pz - M.c11 * t.mx - M.c12 * t.my -
                     M.c13 * t.mz + M.c14 * t.pw;
    ctype const c1 = M.c1 * t.pz - M.c2 * t.py - M.c3 * t.vx;
    ctype const c2 = -M.c0 * t.pz + M.c2 * t.px - M.c3 * t.vy;
    ctype const c3 = M.c0 * t.py - M.c1 * t.px - M.c3 * t.vz;
    ctype const c4 = -M.c0 * t.pw - M.c3 * t.mx - M.c4 * t.px;
    ctype const c5 = -M.c1 * t.pw - M.c3 * t.my - M.c4 * t.py;
    ctype const c6 = -M.c2 * t.pw - M.c3 * t.mz - M.c4 * t.pz;
    ctype const c7 = -M.c1 * t.mz + M.c2 * t.my - M.c4 * t.vx;
    ctype const c8 = M.c0 * t.mz - M.c2 * t.mx - M.c4 * t.vy;
    ctype const c9 = -M.c0 * t.my + M.c1 * t.mx - M.c4 * t.vz;
    ctype const c10 = -M.c0 * t.vx - M.c1 * t.vy - M.c2 * t.vz;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(TriVec3dc<T> const& t,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c5 - t.vy * M.c6 - t.vz * M.c7 - t.mx * M.c11 -
                     t.my * M.c12 - t.mz * M.c13 - t.px * M.c8 - t.py * M.c9 -
                     t.pz * M.c10 + t.pw * M.c14;
    ctype const c1 = t.px * M.c15;
    ctype const c2 = t.py * M.c15;
    ctype const c3 = t.pz * M.c15;
    ctype const c4 = t.vx * M.c15;
    ctype const c5 = t.vy * M.c15;
    ctype const c6 = t.vz * M.c15;
    ctype const c7 = t.mx * M.c15;
    ctype const c8 = t.my * M.c15;
    ctype const c9 = t.mz * M.c15;
    ctype const c10 = -t.pw * M.c15;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& M,
                                                         BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.mz - M.c2 * B.my - M.c3 * B.px + M.c4 * B.vx;
    ctype const c1 = -M.c0 * B.mz + M.c2 * B.mx - M.c3 * B.py + M.c4 * B.vy;
    ctype const c2 = M.c0 * B.my - M.c1 * B.mx - M.c3 * B.pz + M.c4 * B.vz;
    ctype const c3 = M.c0 * B.vx + M.c1 * B.vy + M.c2 * B.vz - M.c3 * B.pw;
    ctype const c4 = -M.c0 * B.px - M.c1 * B.py - M.c2 * B.pz + M.c4 * B.pw;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(BiVec3dc<T> const& B,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * M.c10 - B.vz * M.c9 - B.mx * M.c14 - B.py * M.c13 +
                     B.pz * M.c12 + B.pw * M.c5;
    ctype const c1 = -B.vx * M.c10 + B.vz * M.c8 - B.my * M.c14 + B.px * M.c13 -
                     B.pz * M.c11 + B.pw * M.c6;
    ctype const c2 = B.vx * M.c9 - B.vy * M.c8 - B.mz * M.c14 - B.px * M.c12 +
                     B.py * M.c11 + B.pw * M.c7;
    ctype const c3 = B.vx * M.c5 + B.vy * M.c6 + B.vz * M.c7 + B.mx * M.c11 +
                     B.my * M.c12 + B.mz * M.c13;
    ctype const c4 = B.mx * M.c8 + B.my * M.c9 + B.mz * M.c10 + B.px * M.c5 +
                     B.py * M.c6 + B.pz * M.c7;
    ctype const c5 = -B.mx * M.c15;
    ctype const c6 = -B.my * M.c15;
    ctype const c7 = -B.mz * M.c15;
    ctype const c8 = -B.px * M.c15;
    ctype const c9 = -B.py * M.c15;
    ctype const c10 = -B.pz * M.c15;
    ctype const c11 = -B.vx * M.c15;
    ctype const c12 = -B.vy * M.c15;
    ctype const c13 = -B.vz * M.c15;
    ctype const c14 = B.pw * M.c15;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(MVec3dc_U<T> const& M,
                                                         Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y + M.c2 * v.z - M.c3 * v.u - M.c4 * v.w;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 + v.z * M.c2 - v.w * M.c4 - v.u * M.c3;
    ctype const c1 = v.y * M.c13 - v.z * M.c12 - v.w * M.c5;
    ctype const c2 = -v.x * M.c13 + v.z * M.c11 - v.w * M.c6;
    ctype const c3 = v.x * M.c12 - v.y * M.c11 - v.w * M.c7;
    ctype const c4 = -v.x * M.c14 - v.w * M.c8 - v.u * M.c11;
    ctype const c5 = -v.y * M.c14 - v.w * M.c9 - v.u * M.c12;
    ctype const c6 = -v.z * M.c14 - v.w * M.c10 - v.u * M.c13;
    ctype const c7 = -v.y * M.c10 + v.z * M.c9 - v.u * M.c5;
    ctype const c8 = v.x * M.c10 - v.z * M.c8 - v.u * M.c6;
    ctype const c9 = -v.x * M.c9 + v.y * M.c8 - v.u * M.c7;
    ctype const c10 = -v.x * M.c5 - v.y * M.c6 - v.z * M.c7;
    ctype const c11 = v.x * M.c15;
    ctype const c12 = v.y * M.c15;
    ctype const c13 = v.z * M.c15;
    ctype const c14 = -v.u * M.c15;
    ctype const c15 = -v.w * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc l_contract :: l_contract(mv_u,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] MVec3dc_U<T> const&, [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc l_contract :: l_contract(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator<<(PScalar3dc<T> ps1,
                                                         PScalar3dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga3dc l_contract :: l_contract(ps,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dc<T>, [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(quadvec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator<<(QuadVec3dc<T> const& Q,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * ctype(ps);
    ctype const c1 = -Q.y * ctype(ps);
    ctype const c2 = -Q.z * ctype(ps);
    ctype const c3 = Q.u * ctype(ps);
    ctype const c4 = Q.w * ctype(ps);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc l_contract :: l_contract(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dc<T>, [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(trivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator<<(TriVec3dc<T> const& t,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.px * ctype(ps);
    ctype const c1 = t.py * ctype(ps);
    ctype const c2 = t.pz * ctype(ps);
    ctype const c3 = t.vx * ctype(ps);
    ctype const c4 = t.vy * ctype(ps);
    ctype const c5 = t.vz * ctype(ps);
    ctype const c6 = t.mx * ctype(ps);
    ctype const c7 = t.my * ctype(ps);
    ctype const c8 = t.mz * ctype(ps);
    ctype const c9 = -t.pw * ctype(ps);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc l_contract :: l_contract(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dc<T>, [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(bivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator<<(BiVec3dc<T> const& B,
                                                         PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.mx * ctype(ps);
    ctype const c1 = -B.my * ctype(ps);
    ctype const c2 = -B.mz * ctype(ps);
    ctype const c3 = -B.px * ctype(ps);
    ctype const c4 = -B.py * ctype(ps);
    ctype const c5 = -B.pz * ctype(ps);
    ctype const c6 = -B.vx * ctype(ps);
    ctype const c7 = -B.vy * ctype(ps);
    ctype const c8 = -B.vz * ctype(ps);
    ctype const c9 = B.pw * ctype(ps);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc l_contract :: l_contract(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dc<T>, [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(vec,ps) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v,
                                                          PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = -v.u * ctype(ps);
    ctype const c4 = -v.w * ctype(ps);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc l_contract :: l_contract(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3dc<T>,
                                                         [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s,
                                                          PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(s) * ctype(ps));
}

// cga3dc l_contract :: l_contract(quadvec,quadvec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator<<(QuadVec3dc<T> const& Q1,
                                                         QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-Q1.x * Q2.x - Q1.y * Q2.y - Q1.z * Q2.z + Q1.w * Q2.u +
                            Q1.u * Q2.w);
}

// cga3dc l_contract :: l_contract(quadvec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] QuadVec3dc<T> const&, [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(trivec,quadvec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator<<(TriVec3dc<T> const& t,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vy * Q.z - t.vz * Q.y - t.mx * Q.u + t.px * Q.w;
    ctype const c1 = -t.vx * Q.z + t.vz * Q.x - t.my * Q.u + t.py * Q.w;
    ctype const c2 = t.vx * Q.y - t.vy * Q.x - t.mz * Q.u + t.pz * Q.w;
    ctype const c3 = -t.px * Q.x - t.py * Q.y - t.pz * Q.z + t.pw * Q.u;
    ctype const c4 = t.mx * Q.x + t.my * Q.y + t.mz * Q.z - t.pw * Q.w;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc l_contract :: l_contract(quadvec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] QuadVec3dc<T> const&, [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(bivec,quadvec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator<<(BiVec3dc<T> const& B,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * Q.z - B.vz * Q.y - B.mx * Q.u;
    ctype const c1 = -B.vx * Q.z + B.vz * Q.x - B.my * Q.u;
    ctype const c2 = B.vx * Q.y - B.vy * Q.x - B.mz * Q.u;
    ctype const c3 = -B.vx * Q.w - B.px * Q.u - B.pw * Q.x;
    ctype const c4 = -B.vy * Q.w - B.py * Q.u - B.pw * Q.y;
    ctype const c5 = -B.vz * Q.w - B.pz * Q.u - B.pw * Q.z;
    ctype const c6 = -B.mx * Q.w - B.py * Q.z + B.pz * Q.y;
    ctype const c7 = -B.my * Q.w + B.px * Q.z - B.pz * Q.x;
    ctype const c8 = -B.mz * Q.w - B.px * Q.y + B.py * Q.x;
    ctype const c9 = B.mx * Q.x + B.my * Q.y + B.mz * Q.z;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc l_contract :: l_contract(quadvec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] QuadVec3dc<T> const&, [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(vec,quadvec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v,
                                                         QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * Q.z + v.z * Q.y;
    ctype const c1 = v.x * Q.z - v.z * Q.x;
    ctype const c2 = -v.x * Q.y + v.y * Q.x;
    ctype const c3 = v.x * Q.w + v.u * Q.x;
    ctype const c4 = v.y * Q.w + v.u * Q.y;
    ctype const c5 = v.z * Q.w + v.u * Q.z;
    ctype const c6 = -v.x * Q.u - v.w * Q.x;
    ctype const c7 = -v.y * Q.u - v.w * Q.y;
    ctype const c8 = -v.z * Q.u - v.w * Q.z;
    ctype const c9 = -v.w * Q.w + v.u * Q.u;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc l_contract :: l_contract(quadvec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] QuadVec3dc<T> const&, [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(s,quadvec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s,
                                                          QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * Q.x;
    ctype const c1 = ctype(s) * Q.y;
    ctype const c2 = ctype(s) * Q.z;
    ctype const c3 = ctype(s) * Q.w;
    ctype const c4 = ctype(s) * Q.u;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc l_contract :: l_contract(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator<<(TriVec3dc<T> const& t1,
                                                         TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(
        -t1.vx * t2.vx - t1.vy * t2.vy - t1.vz * t2.vz - t1.mx * t2.px - t1.my * t2.py -
        t1.mz * t2.pz - t1.px * t2.mx - t1.py * t2.my - t1.pz * t2.mz + t1.pw * t2.pw);
}

// cga3dc l_contract :: l_contract(trivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dc<T> const&, [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(bivec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator<<(BiVec3dc<T> const& B,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vy * t.mz - B.vz * t.my - B.mx * t.pw - B.py * t.pz + B.pz * t.py + B.pw * t.vx;
    ctype const c1 = -B.vx * t.mz + B.vz * t.mx - B.my * t.pw + B.px * t.pz -
                     B.pz * t.px + B.pw * t.vy;
    ctype const c2 =
        B.vx * t.my - B.vy * t.mx - B.mz * t.pw - B.px * t.py + B.py * t.px + B.pw * t.vz;
    ctype const c3 =
        B.vx * t.vx + B.vy * t.vy + B.vz * t.vz + B.mx * t.px + B.my * t.py + B.mz * t.pz;
    ctype const c4 =
        B.mx * t.mx + B.my * t.my + B.mz * t.mz + B.px * t.vx + B.py * t.vy + B.pz * t.vz;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc l_contract :: l_contract(trivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dc<T> const&, [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(vec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * t.pz - v.z * t.py - v.w * t.vx;
    ctype const c1 = -v.x * t.pz + v.z * t.px - v.w * t.vy;
    ctype const c2 = v.x * t.py - v.y * t.px - v.w * t.vz;
    ctype const c3 = -v.x * t.pw - v.w * t.mx - v.u * t.px;
    ctype const c4 = -v.y * t.pw - v.w * t.my - v.u * t.py;
    ctype const c5 = -v.z * t.pw - v.w * t.mz - v.u * t.pz;
    ctype const c6 = -v.y * t.mz + v.z * t.my - v.u * t.vx;
    ctype const c7 = v.x * t.mz - v.z * t.mx - v.u * t.vy;
    ctype const c8 = -v.x * t.my + v.y * t.mx - v.u * t.vz;
    ctype const c9 = -v.x * t.vx - v.y * t.vy - v.z * t.vz;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc l_contract :: l_contract(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dc<T> const&, [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s,
                                                         TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.vx;
    ctype const c1 = ctype(s) * t.vy;
    ctype const c2 = ctype(s) * t.vz;
    ctype const c3 = ctype(s) * t.mx;
    ctype const c4 = ctype(s) * t.my;
    ctype const c5 = ctype(s) * t.mz;
    ctype const c6 = ctype(s) * t.px;
    ctype const c7 = ctype(s) * t.py;
    ctype const c8 = ctype(s) * t.pz;
    ctype const c9 = ctype(s) * t.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc l_contract :: l_contract(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator<<(BiVec3dc<T> const& B1,
                                                         BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(
        B1.vx * B2.px + B1.vy * B2.py + B1.vz * B2.pz + B1.mx * B2.mx + B1.my * B2.my +
        B1.mz * B2.mz + B1.px * B2.vx + B1.py * B2.vy + B1.pz * B2.vz - B1.pw * B2.pw);
}

// cga3dc l_contract :: l_contract(bivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3dc<T> const&, [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.mz - v.z * B.my - v.w * B.px + v.u * B.vx;
    ctype const c1 = -v.x * B.mz + v.z * B.mx - v.w * B.py + v.u * B.vy;
    ctype const c2 = v.x * B.my - v.y * B.mx - v.w * B.pz + v.u * B.vz;
    ctype const c3 = v.x * B.vx + v.y * B.vy + v.z * B.vz - v.w * B.pw;
    ctype const c4 = -v.x * B.px - v.y * B.py - v.z * B.pz + v.u * B.pw;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc l_contract :: l_contract(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3dc<T> const&, [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    ctype const c6 = ctype(s) * B.px;
    ctype const c7 = ctype(s) * B.py;
    ctype const c8 = ctype(s) * B.pz;
    ctype const c9 = ctype(s) * B.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc l_contract :: l_contract(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator<<(Vec3dc<T> const& v1,
                                                         Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z - v1.w * v2.u -
                            v1.u * v2.w);
}

// cga3dc l_contract :: l_contract(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec3dc<T> const&, [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc l_contract :: l_contract(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    ctype const c4 = ctype(s) * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc l_contract :: l_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator<<(Scalar3dc<T> s1, Scalar3dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The result lies in A and is perpendicular to B. It is the right interior
// (inner) product: the grade-lowering part <A B>_(grade(A)-grade(B)) of the
// geometric product. On equal grades it reduces to the metric inner product:
//
//     (A >> B) == dot(A, B)
//
// Equivalently, as defined in "PGA Illuminated" (E. Lengyel):
//
//     operator>>(a,b) = r_contract(a,b) = rwdg( a, dual(b) )
//
// (holds exactly at every grade).
////////////////////////////////////////////////////////////////////////////////

// cga3dc r_contract :: r_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& A,
                                                       MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c12 + A.c7 * B.c13 + A.c8 * B.c14 +
                     A.c9 * B.c9 + A.c10 * B.c10 + A.c11 * B.c11 + A.c12 * B.c6 +
                     A.c13 * B.c7 + A.c14 * B.c8 - A.c15 * B.c15 - A.c16 * B.c16 -
                     A.c17 * B.c17 - A.c18 * B.c18 - A.c19 * B.c22 - A.c20 * B.c23 -
                     A.c21 * B.c24 - A.c22 * B.c19 - A.c23 * B.c20 - A.c24 * B.c21 +
                     A.c25 * B.c25 - A.c26 * B.c26 - A.c27 * B.c27 - A.c28 * B.c28 +
                     A.c29 * B.c30 + A.c30 * B.c29 - A.c31 * B.c31;
    ctype const c1 = A.c1 * B.c0 - A.c6 * B.c5 + A.c10 * B.c3 - A.c11 * B.c2 +
                     A.c12 * B.c4 + A.c16 * B.c15 - A.c20 * B.c8 + A.c21 * B.c7 +
                     A.c23 * B.c14 - A.c24 * B.c13 - A.c25 * B.c9 + A.c27 * B.c18 -
                     A.c28 * B.c17 - A.c29 * B.c22 + A.c30 * B.c19 - A.c31 * B.c26;
    ctype const c2 = A.c2 * B.c0 - A.c7 * B.c5 - A.c9 * B.c3 + A.c11 * B.c1 +
                     A.c13 * B.c4 + A.c17 * B.c15 + A.c19 * B.c8 - A.c21 * B.c6 -
                     A.c22 * B.c14 + A.c24 * B.c12 - A.c25 * B.c10 - A.c26 * B.c18 +
                     A.c28 * B.c16 - A.c29 * B.c23 + A.c30 * B.c20 - A.c31 * B.c27;
    ctype const c3 = A.c3 * B.c0 - A.c8 * B.c5 + A.c9 * B.c2 - A.c10 * B.c1 +
                     A.c14 * B.c4 + A.c18 * B.c15 - A.c19 * B.c7 + A.c20 * B.c6 +
                     A.c22 * B.c13 - A.c23 * B.c12 - A.c25 * B.c11 + A.c26 * B.c17 -
                     A.c27 * B.c16 - A.c29 * B.c24 + A.c30 * B.c21 - A.c31 * B.c28;
    ctype const c4 = A.c4 * B.c0 - A.c6 * B.c1 - A.c7 * B.c2 - A.c8 * B.c3 +
                     A.c15 * B.c4 + A.c16 * B.c6 + A.c17 * B.c7 + A.c18 * B.c8 +
                     A.c22 * B.c9 + A.c23 * B.c10 + A.c24 * B.c11 + A.c26 * B.c22 +
                     A.c27 * B.c23 + A.c28 * B.c24 - A.c30 * B.c25 + A.c31 * B.c30;
    ctype const c5 = A.c5 * B.c0 + A.c12 * B.c1 + A.c13 * B.c2 + A.c14 * B.c3 -
                     A.c15 * B.c5 + A.c16 * B.c12 + A.c17 * B.c13 + A.c18 * B.c14 +
                     A.c19 * B.c9 + A.c20 * B.c10 + A.c21 * B.c11 - A.c26 * B.c19 -
                     A.c27 * B.c20 - A.c28 * B.c21 + A.c29 * B.c25 + A.c31 * B.c29;
    ctype const c6 = A.c6 * B.c0 - A.c16 * B.c4 - A.c23 * B.c3 + A.c24 * B.c2 -
                     A.c27 * B.c8 + A.c28 * B.c7 - A.c30 * B.c9 + A.c31 * B.c22;
    ctype const c7 = A.c7 * B.c0 - A.c17 * B.c4 + A.c22 * B.c3 - A.c24 * B.c1 +
                     A.c26 * B.c8 - A.c28 * B.c6 - A.c30 * B.c10 + A.c31 * B.c23;
    ctype const c8 = A.c8 * B.c0 - A.c18 * B.c4 - A.c22 * B.c2 + A.c23 * B.c1 -
                     A.c26 * B.c7 + A.c27 * B.c6 - A.c30 * B.c11 + A.c31 * B.c24;
    ctype const c9 = A.c9 * B.c0 - A.c19 * B.c4 - A.c22 * B.c5 - A.c25 * B.c1 -
                     A.c26 * B.c15 - A.c29 * B.c6 - A.c30 * B.c12 + A.c31 * B.c16;
    ctype const c10 = A.c10 * B.c0 - A.c20 * B.c4 - A.c23 * B.c5 - A.c25 * B.c2 -
                      A.c27 * B.c15 - A.c29 * B.c7 - A.c30 * B.c13 + A.c31 * B.c17;
    ctype const c11 = A.c11 * B.c0 - A.c21 * B.c4 - A.c24 * B.c5 - A.c25 * B.c3 -
                      A.c28 * B.c15 - A.c29 * B.c8 - A.c30 * B.c14 + A.c31 * B.c18;
    ctype const c12 = A.c12 * B.c0 - A.c16 * B.c5 + A.c20 * B.c3 - A.c21 * B.c2 +
                      A.c27 * B.c14 - A.c28 * B.c13 - A.c29 * B.c9 + A.c31 * B.c19;
    ctype const c13 = A.c13 * B.c0 - A.c17 * B.c5 - A.c19 * B.c3 + A.c21 * B.c1 -
                      A.c26 * B.c14 + A.c28 * B.c12 - A.c29 * B.c10 + A.c31 * B.c20;
    ctype const c14 = A.c14 * B.c0 - A.c18 * B.c5 + A.c19 * B.c2 - A.c20 * B.c1 +
                      A.c26 * B.c13 - A.c27 * B.c12 - A.c29 * B.c11 + A.c31 * B.c21;
    ctype const c15 = A.c15 * B.c0 - A.c16 * B.c1 - A.c17 * B.c2 - A.c18 * B.c3 +
                      A.c26 * B.c9 + A.c27 * B.c10 + A.c28 * B.c11 - A.c31 * B.c25;
    ctype const c16 = A.c16 * B.c0 - A.c27 * B.c3 + A.c28 * B.c2 - A.c31 * B.c9;
    ctype const c17 = A.c17 * B.c0 + A.c26 * B.c3 - A.c28 * B.c1 - A.c31 * B.c10;
    ctype const c18 = A.c18 * B.c0 - A.c26 * B.c2 + A.c27 * B.c1 - A.c31 * B.c11;
    ctype const c19 = A.c19 * B.c0 - A.c26 * B.c5 - A.c29 * B.c1 - A.c31 * B.c12;
    ctype const c20 = A.c20 * B.c0 - A.c27 * B.c5 - A.c29 * B.c2 - A.c31 * B.c13;
    ctype const c21 = A.c21 * B.c0 - A.c28 * B.c5 - A.c29 * B.c3 - A.c31 * B.c14;
    ctype const c22 = A.c22 * B.c0 + A.c26 * B.c4 + A.c30 * B.c1 - A.c31 * B.c6;
    ctype const c23 = A.c23 * B.c0 + A.c27 * B.c4 + A.c30 * B.c2 - A.c31 * B.c7;
    ctype const c24 = A.c24 * B.c0 + A.c28 * B.c4 + A.c30 * B.c3 - A.c31 * B.c8;
    ctype const c25 = A.c25 * B.c0 + A.c29 * B.c4 - A.c30 * B.c5 + A.c31 * B.c15;
    ctype const c26 = A.c26 * B.c0 + A.c31 * B.c1;
    ctype const c27 = A.c27 * B.c0 + A.c31 * B.c2;
    ctype const c28 = A.c28 * B.c0 + A.c31 * B.c3;
    ctype const c29 = A.c29 * B.c0 - A.c31 * B.c5;
    ctype const c30 = A.c30 * B.c0 - A.c31 * B.c4;
    ctype const c31 = A.c31 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& A,
                                                       MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c6 * B.c7 + A.c7 * B.c8 + A.c8 * B.c9 + A.c9 * B.c4 +
                     A.c10 * B.c5 + A.c11 * B.c6 + A.c12 * B.c1 + A.c13 * B.c2 +
                     A.c14 * B.c3 - A.c15 * B.c10 - A.c26 * B.c11 - A.c27 * B.c12 -
                     A.c28 * B.c13 + A.c29 * B.c15 + A.c30 * B.c14;
    ctype const c1 = A.c1 * B.c0 + A.c16 * B.c10 - A.c20 * B.c3 + A.c21 * B.c2 +
                     A.c23 * B.c9 - A.c24 * B.c8 - A.c25 * B.c4 - A.c31 * B.c11;
    ctype const c2 = A.c2 * B.c0 + A.c17 * B.c10 + A.c19 * B.c3 - A.c21 * B.c1 -
                     A.c22 * B.c9 + A.c24 * B.c7 - A.c25 * B.c5 - A.c31 * B.c12;
    ctype const c3 = A.c3 * B.c0 + A.c18 * B.c10 - A.c19 * B.c2 + A.c20 * B.c1 +
                     A.c22 * B.c8 - A.c23 * B.c7 - A.c25 * B.c6 - A.c31 * B.c13;
    ctype const c4 = A.c4 * B.c0 + A.c16 * B.c1 + A.c17 * B.c2 + A.c18 * B.c3 +
                     A.c22 * B.c4 + A.c23 * B.c5 + A.c24 * B.c6 + A.c31 * B.c15;
    ctype const c5 = A.c5 * B.c0 + A.c16 * B.c7 + A.c17 * B.c8 + A.c18 * B.c9 +
                     A.c19 * B.c4 + A.c20 * B.c5 + A.c21 * B.c6 + A.c31 * B.c14;
    ctype const c6 = A.c6 * B.c0 - A.c27 * B.c3 + A.c28 * B.c2 - A.c30 * B.c4;
    ctype const c7 = A.c7 * B.c0 + A.c26 * B.c3 - A.c28 * B.c1 - A.c30 * B.c5;
    ctype const c8 = A.c8 * B.c0 - A.c26 * B.c2 + A.c27 * B.c1 - A.c30 * B.c6;
    ctype const c9 = A.c9 * B.c0 - A.c26 * B.c10 - A.c29 * B.c1 - A.c30 * B.c7;
    ctype const c10 = A.c10 * B.c0 - A.c27 * B.c10 - A.c29 * B.c2 - A.c30 * B.c8;
    ctype const c11 = A.c11 * B.c0 - A.c28 * B.c10 - A.c29 * B.c3 - A.c30 * B.c9;
    ctype const c12 = A.c12 * B.c0 + A.c27 * B.c9 - A.c28 * B.c8 - A.c29 * B.c4;
    ctype const c13 = A.c13 * B.c0 - A.c26 * B.c9 + A.c28 * B.c7 - A.c29 * B.c5;
    ctype const c14 = A.c14 * B.c0 + A.c26 * B.c8 - A.c27 * B.c7 - A.c29 * B.c6;
    ctype const c15 = A.c15 * B.c0 + A.c26 * B.c4 + A.c27 * B.c5 + A.c28 * B.c6;
    ctype const c16 = A.c16 * B.c0 - A.c31 * B.c4;
    ctype const c17 = A.c17 * B.c0 - A.c31 * B.c5;
    ctype const c18 = A.c18 * B.c0 - A.c31 * B.c6;
    ctype const c19 = A.c19 * B.c0 - A.c31 * B.c7;
    ctype const c20 = A.c20 * B.c0 - A.c31 * B.c8;
    ctype const c21 = A.c21 * B.c0 - A.c31 * B.c9;
    ctype const c22 = A.c22 * B.c0 - A.c31 * B.c1;
    ctype const c23 = A.c23 * B.c0 - A.c31 * B.c2;
    ctype const c24 = A.c24 * B.c0 - A.c31 * B.c3;
    ctype const c25 = A.c25 * B.c0 + A.c31 * B.c10;
    ctype const c26 = A.c26 * B.c0;
    ctype const c27 = A.c27 * B.c0;
    ctype const c28 = A.c28 * B.c0;
    ctype const c29 = A.c29 * B.c0;
    ctype const c30 = A.c30 * B.c0;
    ctype const c31 = A.c31 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& A,
                                                       MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                     A.c4 * B.c9 + A.c5 * B.c10 + A.c6 * B.c11 + A.c7 * B.c6 +
                     A.c8 * B.c7 + A.c9 * B.c8 - A.c10 * B.c15 - A.c11 * B.c26 -
                     A.c12 * B.c27 - A.c13 * B.c28 + A.c14 * B.c30 + A.c15 * B.c29;
    ctype const c1 = -A.c1 * B.c5 + A.c5 * B.c3 - A.c6 * B.c2 + A.c7 * B.c4 +
                     A.c12 * B.c18 - A.c13 * B.c17 - A.c14 * B.c22 + A.c15 * B.c19;
    ctype const c2 = -A.c2 * B.c5 - A.c4 * B.c3 + A.c6 * B.c1 + A.c8 * B.c4 -
                     A.c11 * B.c18 + A.c13 * B.c16 - A.c14 * B.c23 + A.c15 * B.c20;
    ctype const c3 = -A.c3 * B.c5 + A.c4 * B.c2 - A.c5 * B.c1 + A.c9 * B.c4 +
                     A.c11 * B.c17 - A.c12 * B.c16 - A.c14 * B.c24 + A.c15 * B.c21;
    ctype const c4 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c10 * B.c4 +
                     A.c11 * B.c22 + A.c12 * B.c23 + A.c13 * B.c24 - A.c15 * B.c25;
    ctype const c5 = A.c7 * B.c1 + A.c8 * B.c2 + A.c9 * B.c3 - A.c10 * B.c5 -
                     A.c11 * B.c19 - A.c12 * B.c20 - A.c13 * B.c21 + A.c14 * B.c25;
    ctype const c6 = A.c1 * B.c0 - A.c12 * B.c8 + A.c13 * B.c7 - A.c15 * B.c9;
    ctype const c7 = A.c2 * B.c0 + A.c11 * B.c8 - A.c13 * B.c6 - A.c15 * B.c10;
    ctype const c8 = A.c3 * B.c0 - A.c11 * B.c7 + A.c12 * B.c6 - A.c15 * B.c11;
    ctype const c9 = A.c4 * B.c0 - A.c11 * B.c15 - A.c14 * B.c6 - A.c15 * B.c12;
    ctype const c10 = A.c5 * B.c0 - A.c12 * B.c15 - A.c14 * B.c7 - A.c15 * B.c13;
    ctype const c11 = A.c6 * B.c0 - A.c13 * B.c15 - A.c14 * B.c8 - A.c15 * B.c14;
    ctype const c12 = A.c7 * B.c0 + A.c12 * B.c14 - A.c13 * B.c13 - A.c14 * B.c9;
    ctype const c13 = A.c8 * B.c0 - A.c11 * B.c14 + A.c13 * B.c12 - A.c14 * B.c10;
    ctype const c14 = A.c9 * B.c0 + A.c11 * B.c13 - A.c12 * B.c12 - A.c14 * B.c11;
    ctype const c15 = A.c10 * B.c0 + A.c11 * B.c9 + A.c12 * B.c10 + A.c13 * B.c11;
    ctype const c16 = -A.c12 * B.c3 + A.c13 * B.c2;
    ctype const c17 = A.c11 * B.c3 - A.c13 * B.c1;
    ctype const c18 = -A.c11 * B.c2 + A.c12 * B.c1;
    ctype const c19 = -A.c11 * B.c5 - A.c14 * B.c1;
    ctype const c20 = -A.c12 * B.c5 - A.c14 * B.c2;
    ctype const c21 = -A.c13 * B.c5 - A.c14 * B.c3;
    ctype const c22 = A.c11 * B.c4 + A.c15 * B.c1;
    ctype const c23 = A.c12 * B.c4 + A.c15 * B.c2;
    ctype const c24 = A.c13 * B.c4 + A.c15 * B.c3;
    ctype const c25 = A.c14 * B.c4 - A.c15 * B.c5;
    ctype const c26 = A.c11 * B.c0;
    ctype const c27 = A.c12 * B.c0;
    ctype const c28 = A.c13 * B.c0;
    ctype const c29 = A.c14 * B.c0;
    ctype const c30 = A.c15 * B.c0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& A,
                                                       MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c4 * B.c4 - A.c5 * B.c3 -
                     A.c16 * B.c5 - A.c17 * B.c6 - A.c18 * B.c7 - A.c19 * B.c11 -
                     A.c20 * B.c12 - A.c21 * B.c13 - A.c22 * B.c8 - A.c23 * B.c9 -
                     A.c24 * B.c10 + A.c25 * B.c14 - A.c31 * B.c15;
    ctype const c1 = -A.c6 * B.c4 + A.c10 * B.c2 - A.c11 * B.c1 + A.c12 * B.c3 +
                     A.c27 * B.c7 - A.c28 * B.c6 - A.c29 * B.c11 + A.c30 * B.c8;
    ctype const c2 = -A.c7 * B.c4 - A.c9 * B.c2 + A.c11 * B.c0 + A.c13 * B.c3 -
                     A.c26 * B.c7 + A.c28 * B.c5 - A.c29 * B.c12 + A.c30 * B.c9;
    ctype const c3 = -A.c8 * B.c4 + A.c9 * B.c1 - A.c10 * B.c0 + A.c14 * B.c3 +
                     A.c26 * B.c6 - A.c27 * B.c5 - A.c29 * B.c13 + A.c30 * B.c10;
    ctype const c4 = -A.c6 * B.c0 - A.c7 * B.c1 - A.c8 * B.c2 + A.c15 * B.c3 +
                     A.c26 * B.c11 + A.c27 * B.c12 + A.c28 * B.c13 - A.c30 * B.c14;
    ctype const c5 = A.c12 * B.c0 + A.c13 * B.c1 + A.c14 * B.c2 - A.c15 * B.c4 -
                     A.c26 * B.c8 - A.c27 * B.c9 - A.c28 * B.c10 + A.c29 * B.c14;
    ctype const c6 = -A.c16 * B.c3 - A.c23 * B.c2 + A.c24 * B.c1 + A.c31 * B.c11;
    ctype const c7 = -A.c17 * B.c3 + A.c22 * B.c2 - A.c24 * B.c0 + A.c31 * B.c12;
    ctype const c8 = -A.c18 * B.c3 - A.c22 * B.c1 + A.c23 * B.c0 + A.c31 * B.c13;
    ctype const c9 = -A.c19 * B.c3 - A.c22 * B.c4 - A.c25 * B.c0 + A.c31 * B.c5;
    ctype const c10 = -A.c20 * B.c3 - A.c23 * B.c4 - A.c25 * B.c1 + A.c31 * B.c6;
    ctype const c11 = -A.c21 * B.c3 - A.c24 * B.c4 - A.c25 * B.c2 + A.c31 * B.c7;
    ctype const c12 = -A.c16 * B.c4 + A.c20 * B.c2 - A.c21 * B.c1 + A.c31 * B.c8;
    ctype const c13 = -A.c17 * B.c4 - A.c19 * B.c2 + A.c21 * B.c0 + A.c31 * B.c9;
    ctype const c14 = -A.c18 * B.c4 + A.c19 * B.c1 - A.c20 * B.c0 + A.c31 * B.c10;
    ctype const c15 = -A.c16 * B.c0 - A.c17 * B.c1 - A.c18 * B.c2 - A.c31 * B.c14;
    ctype const c16 = -A.c27 * B.c2 + A.c28 * B.c1;
    ctype const c17 = A.c26 * B.c2 - A.c28 * B.c0;
    ctype const c18 = -A.c26 * B.c1 + A.c27 * B.c0;
    ctype const c19 = -A.c26 * B.c4 - A.c29 * B.c0;
    ctype const c20 = -A.c27 * B.c4 - A.c29 * B.c1;
    ctype const c21 = -A.c28 * B.c4 - A.c29 * B.c2;
    ctype const c22 = A.c26 * B.c3 + A.c30 * B.c0;
    ctype const c23 = A.c27 * B.c3 + A.c30 * B.c1;
    ctype const c24 = A.c28 * B.c3 + A.c30 * B.c2;
    ctype const c25 = A.c29 * B.c3 - A.c30 * B.c4;
    ctype const c26 = A.c31 * B.c0;
    ctype const c27 = A.c31 * B.c1;
    ctype const c28 = A.c31 * B.c2;
    ctype const c29 = -A.c31 * B.c4;
    ctype const c30 = -A.c31 * B.c3;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& A,
                                                       MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c3 * B.c5 - A.c4 * B.c4 -
                     A.c5 * B.c16 - A.c6 * B.c17 - A.c7 * B.c18 - A.c8 * B.c22 -
                     A.c9 * B.c23 - A.c10 * B.c24 - A.c11 * B.c19 - A.c12 * B.c20 -
                     A.c13 * B.c21 + A.c14 * B.c25 - A.c15 * B.c31;
    ctype const c1 = A.c0 * B.c0 + A.c5 * B.c15 - A.c9 * B.c8 + A.c10 * B.c7 +
                     A.c12 * B.c14 - A.c13 * B.c13 - A.c14 * B.c9 - A.c15 * B.c26;
    ctype const c2 = A.c1 * B.c0 + A.c6 * B.c15 + A.c8 * B.c8 - A.c10 * B.c6 -
                     A.c11 * B.c14 + A.c13 * B.c12 - A.c14 * B.c10 - A.c15 * B.c27;
    ctype const c3 = A.c2 * B.c0 + A.c7 * B.c15 - A.c8 * B.c7 + A.c9 * B.c6 +
                     A.c11 * B.c13 - A.c12 * B.c12 - A.c14 * B.c11 - A.c15 * B.c28;
    ctype const c4 = A.c3 * B.c0 + A.c5 * B.c6 + A.c6 * B.c7 + A.c7 * B.c8 +
                     A.c11 * B.c9 + A.c12 * B.c10 + A.c13 * B.c11 + A.c15 * B.c30;
    ctype const c5 = A.c4 * B.c0 + A.c5 * B.c12 + A.c6 * B.c13 + A.c7 * B.c14 +
                     A.c8 * B.c9 + A.c9 * B.c10 + A.c10 * B.c11 + A.c15 * B.c29;
    ctype const c6 = -A.c5 * B.c4 - A.c12 * B.c3 + A.c13 * B.c2 + A.c15 * B.c22;
    ctype const c7 = -A.c6 * B.c4 + A.c11 * B.c3 - A.c13 * B.c1 + A.c15 * B.c23;
    ctype const c8 = -A.c7 * B.c4 - A.c11 * B.c2 + A.c12 * B.c1 + A.c15 * B.c24;
    ctype const c9 = -A.c8 * B.c4 - A.c11 * B.c5 - A.c14 * B.c1 + A.c15 * B.c16;
    ctype const c10 = -A.c9 * B.c4 - A.c12 * B.c5 - A.c14 * B.c2 + A.c15 * B.c17;
    ctype const c11 = -A.c10 * B.c4 - A.c13 * B.c5 - A.c14 * B.c3 + A.c15 * B.c18;
    ctype const c12 = -A.c5 * B.c5 + A.c9 * B.c3 - A.c10 * B.c2 + A.c15 * B.c19;
    ctype const c13 = -A.c6 * B.c5 - A.c8 * B.c3 + A.c10 * B.c1 + A.c15 * B.c20;
    ctype const c14 = -A.c7 * B.c5 + A.c8 * B.c2 - A.c9 * B.c1 + A.c15 * B.c21;
    ctype const c15 = -A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 - A.c15 * B.c25;
    ctype const c16 = A.c5 * B.c0 - A.c15 * B.c9;
    ctype const c17 = A.c6 * B.c0 - A.c15 * B.c10;
    ctype const c18 = A.c7 * B.c0 - A.c15 * B.c11;
    ctype const c19 = A.c8 * B.c0 - A.c15 * B.c12;
    ctype const c20 = A.c9 * B.c0 - A.c15 * B.c13;
    ctype const c21 = A.c10 * B.c0 - A.c15 * B.c14;
    ctype const c22 = A.c11 * B.c0 - A.c15 * B.c6;
    ctype const c23 = A.c12 * B.c0 - A.c15 * B.c7;
    ctype const c24 = A.c13 * B.c0 - A.c15 * B.c8;
    ctype const c25 = A.c14 * B.c0 + A.c15 * B.c15;
    ctype const c26 = A.c15 * B.c1;
    ctype const c27 = A.c15 * B.c2;
    ctype const c28 = A.c15 * B.c3;
    ctype const c29 = -A.c15 * B.c5;
    ctype const c30 = -A.c15 * B.c4;
    ctype const c31 = A.c15 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& M,
                                                       PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c31 * ctype(ps);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c31;
    ctype const c1 = -ctype(ps) * M.c26;
    ctype const c2 = -ctype(ps) * M.c27;
    ctype const c3 = -ctype(ps) * M.c28;
    ctype const c4 = ctype(ps) * M.c30;
    ctype const c5 = ctype(ps) * M.c29;
    ctype const c6 = ctype(ps) * M.c22;
    ctype const c7 = ctype(ps) * M.c23;
    ctype const c8 = ctype(ps) * M.c24;
    ctype const c9 = ctype(ps) * M.c16;
    ctype const c10 = ctype(ps) * M.c17;
    ctype const c11 = ctype(ps) * M.c18;
    ctype const c12 = ctype(ps) * M.c19;
    ctype const c13 = ctype(ps) * M.c20;
    ctype const c14 = ctype(ps) * M.c21;
    ctype const c15 = -ctype(ps) * M.c25;
    ctype const c16 = -ctype(ps) * M.c9;
    ctype const c17 = -ctype(ps) * M.c10;
    ctype const c18 = -ctype(ps) * M.c11;
    ctype const c19 = -ctype(ps) * M.c12;
    ctype const c20 = -ctype(ps) * M.c13;
    ctype const c21 = -ctype(ps) * M.c14;
    ctype const c22 = -ctype(ps) * M.c6;
    ctype const c23 = -ctype(ps) * M.c7;
    ctype const c24 = -ctype(ps) * M.c8;
    ctype const c25 = ctype(ps) * M.c15;
    ctype const c26 = ctype(ps) * M.c1;
    ctype const c27 = ctype(ps) * M.c2;
    ctype const c28 = ctype(ps) * M.c3;
    ctype const c29 = -ctype(ps) * M.c5;
    ctype const c30 = -ctype(ps) * M.c4;
    ctype const c31 = ctype(ps) * M.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,quadvec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& M,
                                                       QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c26 * Q.x - M.c27 * Q.y - M.c28 * Q.z + M.c29 * Q.u + M.c30 * Q.w;
    ctype const c1 = -M.c31 * Q.x;
    ctype const c2 = -M.c31 * Q.y;
    ctype const c3 = -M.c31 * Q.z;
    ctype const c4 = M.c31 * Q.u;
    ctype const c5 = M.c31 * Q.w;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(quadvec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c26 - Q.y * M.c27 - Q.z * M.c28 + Q.w * M.c30 + Q.u * M.c29;
    ctype const c1 = Q.y * M.c18 - Q.z * M.c17 - Q.w * M.c22 + Q.u * M.c19;
    ctype const c2 = -Q.x * M.c18 + Q.z * M.c16 - Q.w * M.c23 + Q.u * M.c20;
    ctype const c3 = Q.x * M.c17 - Q.y * M.c16 - Q.w * M.c24 + Q.u * M.c21;
    ctype const c4 = Q.x * M.c22 + Q.y * M.c23 + Q.z * M.c24 - Q.u * M.c25;
    ctype const c5 = -Q.x * M.c19 - Q.y * M.c20 - Q.z * M.c21 + Q.w * M.c25;
    ctype const c6 = -Q.y * M.c8 + Q.z * M.c7 - Q.u * M.c9;
    ctype const c7 = Q.x * M.c8 - Q.z * M.c6 - Q.u * M.c10;
    ctype const c8 = -Q.x * M.c7 + Q.y * M.c6 - Q.u * M.c11;
    ctype const c9 = -Q.x * M.c15 - Q.w * M.c6 - Q.u * M.c12;
    ctype const c10 = -Q.y * M.c15 - Q.w * M.c7 - Q.u * M.c13;
    ctype const c11 = -Q.z * M.c15 - Q.w * M.c8 - Q.u * M.c14;
    ctype const c12 = Q.y * M.c14 - Q.z * M.c13 - Q.w * M.c9;
    ctype const c13 = -Q.x * M.c14 + Q.z * M.c12 - Q.w * M.c10;
    ctype const c14 = Q.x * M.c13 - Q.y * M.c12 - Q.w * M.c11;
    ctype const c15 = Q.x * M.c9 + Q.y * M.c10 + Q.z * M.c11;
    ctype const c16 = -Q.y * M.c3 + Q.z * M.c2;
    ctype const c17 = Q.x * M.c3 - Q.z * M.c1;
    ctype const c18 = -Q.x * M.c2 + Q.y * M.c1;
    ctype const c19 = -Q.x * M.c5 - Q.w * M.c1;
    ctype const c20 = -Q.y * M.c5 - Q.w * M.c2;
    ctype const c21 = -Q.z * M.c5 - Q.w * M.c3;
    ctype const c22 = Q.x * M.c4 + Q.u * M.c1;
    ctype const c23 = Q.y * M.c4 + Q.u * M.c2;
    ctype const c24 = Q.z * M.c4 + Q.u * M.c3;
    ctype const c25 = Q.w * M.c4 - Q.u * M.c5;
    ctype const c26 = Q.x * M.c0;
    ctype const c27 = Q.y * M.c0;
    ctype const c28 = Q.z * M.c0;
    ctype const c29 = Q.w * M.c0;
    ctype const c30 = Q.u * M.c0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& M,
                                                       TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c16 * t.vx - M.c17 * t.vy - M.c18 * t.vz - M.c19 * t.px -
                     M.c20 * t.py - M.c21 * t.pz - M.c22 * t.mx - M.c23 * t.my -
                     M.c24 * t.mz + M.c25 * t.pw;
    ctype const c1 = M.c27 * t.vz - M.c28 * t.vy - M.c29 * t.px + M.c30 * t.mx;
    ctype const c2 = -M.c26 * t.vz + M.c28 * t.vx - M.c29 * t.py + M.c30 * t.my;
    ctype const c3 = M.c26 * t.vy - M.c27 * t.vx - M.c29 * t.pz + M.c30 * t.mz;
    ctype const c4 = M.c26 * t.px + M.c27 * t.py + M.c28 * t.pz - M.c30 * t.pw;
    ctype const c5 = -M.c26 * t.mx - M.c27 * t.my - M.c28 * t.mz + M.c29 * t.pw;
    ctype const c6 = M.c31 * t.px;
    ctype const c7 = M.c31 * t.py;
    ctype const c8 = M.c31 * t.pz;
    ctype const c9 = M.c31 * t.vx;
    ctype const c10 = M.c31 * t.vy;
    ctype const c11 = M.c31 * t.vz;
    ctype const c12 = M.c31 * t.mx;
    ctype const c13 = M.c31 * t.my;
    ctype const c14 = M.c31 * t.mz;
    ctype const c15 = -M.c31 * t.pw;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(TriVec3dc<T> const& t,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c16 - t.vy * M.c17 - t.vz * M.c18 - t.mx * M.c22 -
                     t.my * M.c23 - t.mz * M.c24 - t.px * M.c19 - t.py * M.c20 -
                     t.pz * M.c21 + t.pw * M.c25;
    ctype const c1 = t.vx * M.c15 - t.my * M.c8 + t.mz * M.c7 + t.py * M.c14 -
                     t.pz * M.c13 - t.pw * M.c9;
    ctype const c2 = t.vy * M.c15 + t.mx * M.c8 - t.mz * M.c6 - t.px * M.c14 +
                     t.pz * M.c12 - t.pw * M.c10;
    ctype const c3 = t.vz * M.c15 - t.mx * M.c7 + t.my * M.c6 + t.px * M.c13 -
                     t.py * M.c12 - t.pw * M.c11;
    ctype const c4 = t.vx * M.c6 + t.vy * M.c7 + t.vz * M.c8 + t.px * M.c9 +
                     t.py * M.c10 + t.pz * M.c11;
    ctype const c5 = t.vx * M.c12 + t.vy * M.c13 + t.vz * M.c14 + t.mx * M.c9 +
                     t.my * M.c10 + t.mz * M.c11;
    ctype const c6 = -t.vx * M.c4 - t.py * M.c3 + t.pz * M.c2;
    ctype const c7 = -t.vy * M.c4 + t.px * M.c3 - t.pz * M.c1;
    ctype const c8 = -t.vz * M.c4 - t.px * M.c2 + t.py * M.c1;
    ctype const c9 = -t.mx * M.c4 - t.px * M.c5 - t.pw * M.c1;
    ctype const c10 = -t.my * M.c4 - t.py * M.c5 - t.pw * M.c2;
    ctype const c11 = -t.mz * M.c4 - t.pz * M.c5 - t.pw * M.c3;
    ctype const c12 = -t.vx * M.c5 + t.my * M.c3 - t.mz * M.c2;
    ctype const c13 = -t.vy * M.c5 - t.mx * M.c3 + t.mz * M.c1;
    ctype const c14 = -t.vz * M.c5 + t.mx * M.c2 - t.my * M.c1;
    ctype const c15 = -t.vx * M.c1 - t.vy * M.c2 - t.vz * M.c3;
    ctype const c16 = t.vx * M.c0;
    ctype const c17 = t.vy * M.c0;
    ctype const c18 = t.vz * M.c0;
    ctype const c19 = t.mx * M.c0;
    ctype const c20 = t.my * M.c0;
    ctype const c21 = t.mz * M.c0;
    ctype const c22 = t.px * M.c0;
    ctype const c23 = t.py * M.c0;
    ctype const c24 = t.pz * M.c0;
    ctype const c25 = t.pw * M.c0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& M,
                                                       BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c6 * B.px + M.c7 * B.py + M.c8 * B.pz + M.c9 * B.mx +
                     M.c10 * B.my + M.c11 * B.mz + M.c12 * B.vx + M.c13 * B.vy +
                     M.c14 * B.vz - M.c15 * B.pw;
    ctype const c1 = M.c16 * B.pw - M.c20 * B.vz + M.c21 * B.vy + M.c23 * B.pz -
                     M.c24 * B.py - M.c25 * B.mx;
    ctype const c2 = M.c17 * B.pw + M.c19 * B.vz - M.c21 * B.vx - M.c22 * B.pz +
                     M.c24 * B.px - M.c25 * B.my;
    ctype const c3 = M.c18 * B.pw - M.c19 * B.vy + M.c20 * B.vx + M.c22 * B.py -
                     M.c23 * B.px - M.c25 * B.mz;
    ctype const c4 = M.c16 * B.vx + M.c17 * B.vy + M.c18 * B.vz + M.c22 * B.mx +
                     M.c23 * B.my + M.c24 * B.mz;
    ctype const c5 = M.c16 * B.px + M.c17 * B.py + M.c18 * B.pz + M.c19 * B.mx +
                     M.c20 * B.my + M.c21 * B.mz;
    ctype const c6 = -M.c27 * B.vz + M.c28 * B.vy - M.c30 * B.mx;
    ctype const c7 = M.c26 * B.vz - M.c28 * B.vx - M.c30 * B.my;
    ctype const c8 = -M.c26 * B.vy + M.c27 * B.vx - M.c30 * B.mz;
    ctype const c9 = -M.c26 * B.pw - M.c29 * B.vx - M.c30 * B.px;
    ctype const c10 = -M.c27 * B.pw - M.c29 * B.vy - M.c30 * B.py;
    ctype const c11 = -M.c28 * B.pw - M.c29 * B.vz - M.c30 * B.pz;
    ctype const c12 = M.c27 * B.pz - M.c28 * B.py - M.c29 * B.mx;
    ctype const c13 = -M.c26 * B.pz + M.c28 * B.px - M.c29 * B.my;
    ctype const c14 = M.c26 * B.py - M.c27 * B.px - M.c29 * B.mz;
    ctype const c15 = M.c26 * B.mx + M.c27 * B.my + M.c28 * B.mz;
    ctype const c16 = -M.c31 * B.mx;
    ctype const c17 = -M.c31 * B.my;
    ctype const c18 = -M.c31 * B.mz;
    ctype const c19 = -M.c31 * B.px;
    ctype const c20 = -M.c31 * B.py;
    ctype const c21 = -M.c31 * B.pz;
    ctype const c22 = -M.c31 * B.vx;
    ctype const c23 = -M.c31 * B.vy;
    ctype const c24 = -M.c31 * B.vz;
    ctype const c25 = M.c31 * B.pw;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(BiVec3dc<T> const& B,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c12 + B.vy * M.c13 + B.vz * M.c14 + B.mx * M.c9 +
                     B.my * M.c10 + B.mz * M.c11 + B.px * M.c6 + B.py * M.c7 +
                     B.pz * M.c8 - B.pw * M.c15;
    ctype const c1 = -B.vx * M.c5 + B.my * M.c3 - B.mz * M.c2 + B.px * M.c4;
    ctype const c2 = -B.vy * M.c5 - B.mx * M.c3 + B.mz * M.c1 + B.py * M.c4;
    ctype const c3 = -B.vz * M.c5 + B.mx * M.c2 - B.my * M.c1 + B.pz * M.c4;
    ctype const c4 = -B.vx * M.c1 - B.vy * M.c2 - B.vz * M.c3 + B.pw * M.c4;
    ctype const c5 = B.px * M.c1 + B.py * M.c2 + B.pz * M.c3 - B.pw * M.c5;
    ctype const c6 = B.vx * M.c0;
    ctype const c7 = B.vy * M.c0;
    ctype const c8 = B.vz * M.c0;
    ctype const c9 = B.mx * M.c0;
    ctype const c10 = B.my * M.c0;
    ctype const c11 = B.mz * M.c0;
    ctype const c12 = B.px * M.c0;
    ctype const c13 = B.py * M.c0;
    ctype const c14 = B.pz * M.c0;
    ctype const c15 = B.pw * M.c0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& M,
                                                       Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z - M.c4 * v.u - M.c5 * v.w;
    ctype const c1 = -M.c6 * v.u + M.c10 * v.z - M.c11 * v.y + M.c12 * v.w;
    ctype const c2 = -M.c7 * v.u - M.c9 * v.z + M.c11 * v.x + M.c13 * v.w;
    ctype const c3 = -M.c8 * v.u + M.c9 * v.y - M.c10 * v.x + M.c14 * v.w;
    ctype const c4 = -M.c6 * v.x - M.c7 * v.y - M.c8 * v.z + M.c15 * v.w;
    ctype const c5 = M.c12 * v.x + M.c13 * v.y + M.c14 * v.z - M.c15 * v.u;
    ctype const c6 = -M.c16 * v.w - M.c23 * v.z + M.c24 * v.y;
    ctype const c7 = -M.c17 * v.w + M.c22 * v.z - M.c24 * v.x;
    ctype const c8 = -M.c18 * v.w - M.c22 * v.y + M.c23 * v.x;
    ctype const c9 = -M.c19 * v.w - M.c22 * v.u - M.c25 * v.x;
    ctype const c10 = -M.c20 * v.w - M.c23 * v.u - M.c25 * v.y;
    ctype const c11 = -M.c21 * v.w - M.c24 * v.u - M.c25 * v.z;
    ctype const c12 = -M.c16 * v.u + M.c20 * v.z - M.c21 * v.y;
    ctype const c13 = -M.c17 * v.u - M.c19 * v.z + M.c21 * v.x;
    ctype const c14 = -M.c18 * v.u + M.c19 * v.y - M.c20 * v.x;
    ctype const c15 = -M.c16 * v.x - M.c17 * v.y - M.c18 * v.z;
    ctype const c16 = -M.c27 * v.z + M.c28 * v.y;
    ctype const c17 = M.c26 * v.z - M.c28 * v.x;
    ctype const c18 = -M.c26 * v.y + M.c27 * v.x;
    ctype const c19 = -M.c26 * v.u - M.c29 * v.x;
    ctype const c20 = -M.c27 * v.u - M.c29 * v.y;
    ctype const c21 = -M.c28 * v.u - M.c29 * v.z;
    ctype const c22 = M.c26 * v.w + M.c30 * v.x;
    ctype const c23 = M.c27 * v.w + M.c30 * v.y;
    ctype const c24 = M.c28 * v.w + M.c30 * v.z;
    ctype const c25 = M.c29 * v.w - M.c30 * v.u;
    ctype const c26 = M.c31 * v.x;
    ctype const c27 = M.c31 * v.y;
    ctype const c28 = M.c31 * v.z;
    ctype const c29 = -M.c31 * v.u;
    ctype const c30 = -M.c31 * v.w;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(Vec3dc<T> const& v,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3 - v.w * M.c5 - v.u * M.c4;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = v.w * M.c0;
    ctype const c5 = v.u * M.c0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(MVec3dc<T> const& M,
                                                       Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    ctype const c16 = M.c16 * ctype(s);
    ctype const c17 = M.c17 * ctype(s);
    ctype const c18 = M.c18 * ctype(s);
    ctype const c19 = M.c19 * ctype(s);
    ctype const c20 = M.c20 * ctype(s);
    ctype const c21 = M.c21 * ctype(s);
    ctype const c22 = M.c22 * ctype(s);
    ctype const c23 = M.c23 * ctype(s);
    ctype const c24 = M.c24 * ctype(s);
    ctype const c25 = M.c25 * ctype(s);
    ctype const c26 = M.c26 * ctype(s);
    ctype const c27 = M.c27 * ctype(s);
    ctype const c28 = M.c28 * ctype(s);
    ctype const c29 = M.c29 * ctype(s);
    ctype const c30 = M.c30 * ctype(s);
    ctype const c31 = M.c31 * ctype(s);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator>>(Scalar3dc<T> s,
                                                       MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    ctype const c16 = 0.0;
    ctype const c17 = 0.0;
    ctype const c18 = 0.0;
    ctype const c19 = 0.0;
    ctype const c20 = 0.0;
    ctype const c21 = 0.0;
    ctype const c22 = 0.0;
    ctype const c23 = 0.0;
    ctype const c24 = 0.0;
    ctype const c25 = 0.0;
    ctype const c26 = 0.0;
    ctype const c27 = 0.0;
    ctype const c28 = 0.0;
    ctype const c29 = 0.0;
    ctype const c30 = 0.0;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc r_contract :: r_contract(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& A,
                                                         MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c7 + A.c2 * B.c8 + A.c3 * B.c9 + A.c4 * B.c4 +
                     A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c1 + A.c8 * B.c2 + A.c9 * B.c3 -
                     A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c12 - A.c13 * B.c13 +
                     A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c1 = A.c1 * B.c0 - A.c12 * B.c3 + A.c13 * B.c2 - A.c15 * B.c4;
    ctype const c2 = A.c2 * B.c0 + A.c11 * B.c3 - A.c13 * B.c1 - A.c15 * B.c5;
    ctype const c3 = A.c3 * B.c0 - A.c11 * B.c2 + A.c12 * B.c1 - A.c15 * B.c6;
    ctype const c4 = A.c4 * B.c0 - A.c11 * B.c10 - A.c14 * B.c1 - A.c15 * B.c7;
    ctype const c5 = A.c5 * B.c0 - A.c12 * B.c10 - A.c14 * B.c2 - A.c15 * B.c8;
    ctype const c6 = A.c6 * B.c0 - A.c13 * B.c10 - A.c14 * B.c3 - A.c15 * B.c9;
    ctype const c7 = A.c7 * B.c0 + A.c12 * B.c9 - A.c13 * B.c8 - A.c14 * B.c4;
    ctype const c8 = A.c8 * B.c0 - A.c11 * B.c9 + A.c13 * B.c7 - A.c14 * B.c5;
    ctype const c9 = A.c9 * B.c0 + A.c11 * B.c8 - A.c12 * B.c7 - A.c14 * B.c6;
    ctype const c10 = A.c10 * B.c0 + A.c11 * B.c4 + A.c12 * B.c5 + A.c13 * B.c6;
    ctype const c11 = A.c11 * B.c0;
    ctype const c12 = A.c12 * B.c0;
    ctype const c13 = A.c13 * B.c0;
    ctype const c14 = A.c14 * B.c0;
    ctype const c15 = A.c15 * B.c0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& A,
                                                         MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * B.c4 + A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c3 +
                     A.c12 * B.c7 - A.c13 * B.c6 - A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c1 = -A.c2 * B.c4 - A.c4 * B.c2 + A.c6 * B.c0 + A.c8 * B.c3 -
                     A.c11 * B.c7 + A.c13 * B.c5 - A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c2 = -A.c3 * B.c4 + A.c4 * B.c1 - A.c5 * B.c0 + A.c9 * B.c3 +
                     A.c11 * B.c6 - A.c12 * B.c5 - A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c3 = -A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2 + A.c10 * B.c3 +
                     A.c11 * B.c11 + A.c12 * B.c12 + A.c13 * B.c13 - A.c15 * B.c14;
    ctype const c4 = A.c7 * B.c0 + A.c8 * B.c1 + A.c9 * B.c2 - A.c10 * B.c4 -
                     A.c11 * B.c8 - A.c12 * B.c9 - A.c13 * B.c10 + A.c14 * B.c14;
    ctype const c5 = -A.c12 * B.c2 + A.c13 * B.c1;
    ctype const c6 = A.c11 * B.c2 - A.c13 * B.c0;
    ctype const c7 = -A.c11 * B.c1 + A.c12 * B.c0;
    ctype const c8 = -A.c11 * B.c4 - A.c14 * B.c0;
    ctype const c9 = -A.c12 * B.c4 - A.c14 * B.c1;
    ctype const c10 = -A.c13 * B.c4 - A.c14 * B.c2;
    ctype const c11 = A.c11 * B.c3 + A.c15 * B.c0;
    ctype const c12 = A.c12 * B.c3 + A.c15 * B.c1;
    ctype const c13 = A.c13 * B.c3 + A.c15 * B.c2;
    ctype const c14 = A.c14 * B.c3 - A.c15 * B.c4;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& A,
                                                         MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c5 * B.c10 - A.c9 * B.c3 + A.c10 * B.c2 +
                     A.c12 * B.c9 - A.c13 * B.c8 - A.c14 * B.c4 - A.c15 * B.c11;
    ctype const c1 = A.c1 * B.c0 + A.c6 * B.c10 + A.c8 * B.c3 - A.c10 * B.c1 -
                     A.c11 * B.c9 + A.c13 * B.c7 - A.c14 * B.c5 - A.c15 * B.c12;
    ctype const c2 = A.c2 * B.c0 + A.c7 * B.c10 - A.c8 * B.c2 + A.c9 * B.c1 +
                     A.c11 * B.c8 - A.c12 * B.c7 - A.c14 * B.c6 - A.c15 * B.c13;
    ctype const c3 = A.c3 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2 + A.c7 * B.c3 +
                     A.c11 * B.c4 + A.c12 * B.c5 + A.c13 * B.c6 + A.c15 * B.c15;
    ctype const c4 = A.c4 * B.c0 + A.c5 * B.c7 + A.c6 * B.c8 + A.c7 * B.c9 + A.c8 * B.c4 +
                     A.c9 * B.c5 + A.c10 * B.c6 + A.c15 * B.c14;
    ctype const c5 = A.c5 * B.c0 - A.c15 * B.c4;
    ctype const c6 = A.c6 * B.c0 - A.c15 * B.c5;
    ctype const c7 = A.c7 * B.c0 - A.c15 * B.c6;
    ctype const c8 = A.c8 * B.c0 - A.c15 * B.c7;
    ctype const c9 = A.c9 * B.c0 - A.c15 * B.c8;
    ctype const c10 = A.c10 * B.c0 - A.c15 * B.c9;
    ctype const c11 = A.c11 * B.c0 - A.c15 * B.c1;
    ctype const c12 = A.c12 * B.c0 - A.c15 * B.c2;
    ctype const c13 = A.c13 * B.c0 - A.c15 * B.c3;
    ctype const c14 = A.c14 * B.c0 + A.c15 * B.c10;
    ctype const c15 = A.c15 * B.c0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(mv_e,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] MVec3dc_E<T> const&, [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(ps,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c11;
    ctype const c1 = -ctype(ps) * M.c12;
    ctype const c2 = -ctype(ps) * M.c13;
    ctype const c3 = ctype(ps) * M.c15;
    ctype const c4 = ctype(ps) * M.c14;
    ctype const c5 = -ctype(ps) * M.c4;
    ctype const c6 = -ctype(ps) * M.c5;
    ctype const c7 = -ctype(ps) * M.c6;
    ctype const c8 = -ctype(ps) * M.c7;
    ctype const c9 = -ctype(ps) * M.c8;
    ctype const c10 = -ctype(ps) * M.c9;
    ctype const c11 = -ctype(ps) * M.c1;
    ctype const c12 = -ctype(ps) * M.c2;
    ctype const c13 = -ctype(ps) * M.c3;
    ctype const c14 = ctype(ps) * M.c10;
    ctype const c15 = ctype(ps) * M.c0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(mv_e,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& M,
                                                         QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * Q.x - M.c12 * Q.y - M.c13 * Q.z + M.c14 * Q.u + M.c15 * Q.w;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(quadvec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c11 - Q.y * M.c12 - Q.z * M.c13 + Q.w * M.c15 + Q.u * M.c14;
    ctype const c1 = -Q.y * M.c3 + Q.z * M.c2 - Q.u * M.c4;
    ctype const c2 = Q.x * M.c3 - Q.z * M.c1 - Q.u * M.c5;
    ctype const c3 = -Q.x * M.c2 + Q.y * M.c1 - Q.u * M.c6;
    ctype const c4 = -Q.x * M.c10 - Q.w * M.c1 - Q.u * M.c7;
    ctype const c5 = -Q.y * M.c10 - Q.w * M.c2 - Q.u * M.c8;
    ctype const c6 = -Q.z * M.c10 - Q.w * M.c3 - Q.u * M.c9;
    ctype const c7 = Q.y * M.c9 - Q.z * M.c8 - Q.w * M.c4;
    ctype const c8 = -Q.x * M.c9 + Q.z * M.c7 - Q.w * M.c5;
    ctype const c9 = Q.x * M.c8 - Q.y * M.c7 - Q.w * M.c6;
    ctype const c10 = Q.x * M.c4 + Q.y * M.c5 + Q.z * M.c6;
    ctype const c11 = Q.x * M.c0;
    ctype const c12 = Q.y * M.c0;
    ctype const c13 = Q.z * M.c0;
    ctype const c14 = Q.w * M.c0;
    ctype const c15 = Q.u * M.c0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& M,
                                                         TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c12 * t.vz - M.c13 * t.vy - M.c14 * t.px + M.c15 * t.mx;
    ctype const c1 = -M.c11 * t.vz + M.c13 * t.vx - M.c14 * t.py + M.c15 * t.my;
    ctype const c2 = M.c11 * t.vy - M.c12 * t.vx - M.c14 * t.pz + M.c15 * t.mz;
    ctype const c3 = M.c11 * t.px + M.c12 * t.py + M.c13 * t.pz - M.c15 * t.pw;
    ctype const c4 = -M.c11 * t.mx - M.c12 * t.my - M.c13 * t.mz + M.c14 * t.pw;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(TriVec3dc<T> const& t,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * M.c10 - t.my * M.c3 + t.mz * M.c2 + t.py * M.c9 -
                     t.pz * M.c8 - t.pw * M.c4;
    ctype const c1 = t.vy * M.c10 + t.mx * M.c3 - t.mz * M.c1 - t.px * M.c9 +
                     t.pz * M.c7 - t.pw * M.c5;
    ctype const c2 = t.vz * M.c10 - t.mx * M.c2 + t.my * M.c1 + t.px * M.c8 -
                     t.py * M.c7 - t.pw * M.c6;
    ctype const c3 =
        t.vx * M.c1 + t.vy * M.c2 + t.vz * M.c3 + t.px * M.c4 + t.py * M.c5 + t.pz * M.c6;
    ctype const c4 =
        t.vx * M.c7 + t.vy * M.c8 + t.vz * M.c9 + t.mx * M.c4 + t.my * M.c5 + t.mz * M.c6;
    ctype const c5 = t.vx * M.c0;
    ctype const c6 = t.vy * M.c0;
    ctype const c7 = t.vz * M.c0;
    ctype const c8 = t.mx * M.c0;
    ctype const c9 = t.my * M.c0;
    ctype const c10 = t.mz * M.c0;
    ctype const c11 = t.px * M.c0;
    ctype const c12 = t.py * M.c0;
    ctype const c13 = t.pz * M.c0;
    ctype const c14 = t.pw * M.c0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& M,
                                                         BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.px + M.c2 * B.py + M.c3 * B.pz + M.c4 * B.mx + M.c5 * B.my +
                     M.c6 * B.mz + M.c7 * B.vx + M.c8 * B.vy + M.c9 * B.vz - M.c10 * B.pw;
    ctype const c1 = -M.c12 * B.vz + M.c13 * B.vy - M.c15 * B.mx;
    ctype const c2 = M.c11 * B.vz - M.c13 * B.vx - M.c15 * B.my;
    ctype const c3 = -M.c11 * B.vy + M.c12 * B.vx - M.c15 * B.mz;
    ctype const c4 = -M.c11 * B.pw - M.c14 * B.vx - M.c15 * B.px;
    ctype const c5 = -M.c12 * B.pw - M.c14 * B.vy - M.c15 * B.py;
    ctype const c6 = -M.c13 * B.pw - M.c14 * B.vz - M.c15 * B.pz;
    ctype const c7 = M.c12 * B.pz - M.c13 * B.py - M.c14 * B.mx;
    ctype const c8 = -M.c11 * B.pz + M.c13 * B.px - M.c14 * B.my;
    ctype const c9 = M.c11 * B.py - M.c12 * B.px - M.c14 * B.mz;
    ctype const c10 = M.c11 * B.mx + M.c12 * B.my + M.c13 * B.mz;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(BiVec3dc<T> const& B,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c7 + B.vy * M.c8 + B.vz * M.c9 + B.mx * M.c4 + B.my * M.c5 +
                     B.mz * M.c6 + B.px * M.c1 + B.py * M.c2 + B.pz * M.c3 - B.pw * M.c10;
    ctype const c1 = B.vx * M.c0;
    ctype const c2 = B.vy * M.c0;
    ctype const c3 = B.vz * M.c0;
    ctype const c4 = B.mx * M.c0;
    ctype const c5 = B.my * M.c0;
    ctype const c6 = B.mz * M.c0;
    ctype const c7 = B.px * M.c0;
    ctype const c8 = B.py * M.c0;
    ctype const c9 = B.pz * M.c0;
    ctype const c10 = B.pw * M.c0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& M,
                                                         Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * v.u + M.c5 * v.z - M.c6 * v.y + M.c7 * v.w;
    ctype const c1 = -M.c2 * v.u - M.c4 * v.z + M.c6 * v.x + M.c8 * v.w;
    ctype const c2 = -M.c3 * v.u + M.c4 * v.y - M.c5 * v.x + M.c9 * v.w;
    ctype const c3 = -M.c1 * v.x - M.c2 * v.y - M.c3 * v.z + M.c10 * v.w;
    ctype const c4 = M.c7 * v.x + M.c8 * v.y + M.c9 * v.z - M.c10 * v.u;
    ctype const c5 = -M.c12 * v.z + M.c13 * v.y;
    ctype const c6 = M.c11 * v.z - M.c13 * v.x;
    ctype const c7 = -M.c11 * v.y + M.c12 * v.x;
    ctype const c8 = -M.c11 * v.u - M.c14 * v.x;
    ctype const c9 = -M.c12 * v.u - M.c14 * v.y;
    ctype const c10 = -M.c13 * v.u - M.c14 * v.z;
    ctype const c11 = M.c11 * v.w + M.c15 * v.x;
    ctype const c12 = M.c12 * v.w + M.c15 * v.y;
    ctype const c13 = M.c13 * v.w + M.c15 * v.z;
    ctype const c14 = M.c14 * v.w - M.c15 * v.u;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(Vec3dc<T> const& v,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    ctype const c3 = v.w * M.c0;
    ctype const c4 = v.u * M.c0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_E<T> const& M,
                                                         Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(Scalar3dc<T> s,
                                                         MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& A,
                                                         MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 - A.c4 * B.c3 -
                     A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7 - A.c8 * B.c11 -
                     A.c9 * B.c12 - A.c10 * B.c13 - A.c11 * B.c8 - A.c12 * B.c9 -
                     A.c13 * B.c10 + A.c14 * B.c14 - A.c15 * B.c15;
    ctype const c1 = -A.c5 * B.c3 - A.c12 * B.c2 + A.c13 * B.c1 + A.c15 * B.c11;
    ctype const c2 = -A.c6 * B.c3 + A.c11 * B.c2 - A.c13 * B.c0 + A.c15 * B.c12;
    ctype const c3 = -A.c7 * B.c3 - A.c11 * B.c1 + A.c12 * B.c0 + A.c15 * B.c13;
    ctype const c4 = -A.c8 * B.c3 - A.c11 * B.c4 - A.c14 * B.c0 + A.c15 * B.c5;
    ctype const c5 = -A.c9 * B.c3 - A.c12 * B.c4 - A.c14 * B.c1 + A.c15 * B.c6;
    ctype const c6 = -A.c10 * B.c3 - A.c13 * B.c4 - A.c14 * B.c2 + A.c15 * B.c7;
    ctype const c7 = -A.c5 * B.c4 + A.c9 * B.c2 - A.c10 * B.c1 + A.c15 * B.c8;
    ctype const c8 = -A.c6 * B.c4 - A.c8 * B.c2 + A.c10 * B.c0 + A.c15 * B.c9;
    ctype const c9 = -A.c7 * B.c4 + A.c8 * B.c1 - A.c9 * B.c0 + A.c15 * B.c10;
    ctype const c10 = -A.c5 * B.c0 - A.c6 * B.c1 - A.c7 * B.c2 - A.c15 * B.c14;
    ctype const c11 = A.c15 * B.c0;
    ctype const c12 = A.c15 * B.c1;
    ctype const c13 = A.c15 * B.c2;
    ctype const c14 = -A.c15 * B.c4;
    ctype const c15 = -A.c15 * B.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_u,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& M,
                                                         PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c15 * ctype(ps);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(ps,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c15;
    ctype const c1 = ctype(ps) * M.c11;
    ctype const c2 = ctype(ps) * M.c12;
    ctype const c3 = ctype(ps) * M.c13;
    ctype const c4 = ctype(ps) * M.c5;
    ctype const c5 = ctype(ps) * M.c6;
    ctype const c6 = ctype(ps) * M.c7;
    ctype const c7 = ctype(ps) * M.c8;
    ctype const c8 = ctype(ps) * M.c9;
    ctype const c9 = ctype(ps) * M.c10;
    ctype const c10 = -ctype(ps) * M.c14;
    ctype const c11 = ctype(ps) * M.c0;
    ctype const c12 = ctype(ps) * M.c1;
    ctype const c13 = ctype(ps) * M.c2;
    ctype const c14 = -ctype(ps) * M.c4;
    ctype const c15 = -ctype(ps) * M.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_u,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& M,
                                                         QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c15 * Q.x;
    ctype const c1 = -M.c15 * Q.y;
    ctype const c2 = -M.c15 * Q.z;
    ctype const c3 = M.c15 * Q.u;
    ctype const c4 = M.c15 * Q.w;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(quadvec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * M.c7 - Q.z * M.c6 - Q.w * M.c11 + Q.u * M.c8;
    ctype const c1 = -Q.x * M.c7 + Q.z * M.c5 - Q.w * M.c12 + Q.u * M.c9;
    ctype const c2 = Q.x * M.c6 - Q.y * M.c5 - Q.w * M.c13 + Q.u * M.c10;
    ctype const c3 = Q.x * M.c11 + Q.y * M.c12 + Q.z * M.c13 - Q.u * M.c14;
    ctype const c4 = -Q.x * M.c8 - Q.y * M.c9 - Q.z * M.c10 + Q.w * M.c14;
    ctype const c5 = -Q.y * M.c2 + Q.z * M.c1;
    ctype const c6 = Q.x * M.c2 - Q.z * M.c0;
    ctype const c7 = -Q.x * M.c1 + Q.y * M.c0;
    ctype const c8 = -Q.x * M.c4 - Q.w * M.c0;
    ctype const c9 = -Q.y * M.c4 - Q.w * M.c1;
    ctype const c10 = -Q.z * M.c4 - Q.w * M.c2;
    ctype const c11 = Q.x * M.c3 + Q.u * M.c0;
    ctype const c12 = Q.y * M.c3 + Q.u * M.c1;
    ctype const c13 = Q.z * M.c3 + Q.u * M.c2;
    ctype const c14 = Q.w * M.c3 - Q.u * M.c4;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& M,
                                                         TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * t.vx - M.c6 * t.vy - M.c7 * t.vz - M.c8 * t.px -
                     M.c9 * t.py - M.c10 * t.pz - M.c11 * t.mx - M.c12 * t.my -
                     M.c13 * t.mz + M.c14 * t.pw;
    ctype const c1 = M.c15 * t.px;
    ctype const c2 = M.c15 * t.py;
    ctype const c3 = M.c15 * t.pz;
    ctype const c4 = M.c15 * t.vx;
    ctype const c5 = M.c15 * t.vy;
    ctype const c6 = M.c15 * t.vz;
    ctype const c7 = M.c15 * t.mx;
    ctype const c8 = M.c15 * t.my;
    ctype const c9 = M.c15 * t.mz;
    ctype const c10 = -M.c15 * t.pw;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(TriVec3dc<T> const& t,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c5 - t.vy * M.c6 - t.vz * M.c7 - t.mx * M.c11 -
                     t.my * M.c12 - t.mz * M.c13 - t.px * M.c8 - t.py * M.c9 -
                     t.pz * M.c10 + t.pw * M.c14;
    ctype const c1 = -t.vx * M.c3 - t.py * M.c2 + t.pz * M.c1;
    ctype const c2 = -t.vy * M.c3 + t.px * M.c2 - t.pz * M.c0;
    ctype const c3 = -t.vz * M.c3 - t.px * M.c1 + t.py * M.c0;
    ctype const c4 = -t.mx * M.c3 - t.px * M.c4 - t.pw * M.c0;
    ctype const c5 = -t.my * M.c3 - t.py * M.c4 - t.pw * M.c1;
    ctype const c6 = -t.mz * M.c3 - t.pz * M.c4 - t.pw * M.c2;
    ctype const c7 = -t.vx * M.c4 + t.my * M.c2 - t.mz * M.c1;
    ctype const c8 = -t.vy * M.c4 - t.mx * M.c2 + t.mz * M.c0;
    ctype const c9 = -t.vz * M.c4 + t.mx * M.c1 - t.my * M.c0;
    ctype const c10 = -t.vx * M.c0 - t.vy * M.c1 - t.vz * M.c2;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& M,
                                                         BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * B.pw - M.c9 * B.vz + M.c10 * B.vy + M.c12 * B.pz -
                     M.c13 * B.py - M.c14 * B.mx;
    ctype const c1 = M.c6 * B.pw + M.c8 * B.vz - M.c10 * B.vx - M.c11 * B.pz +
                     M.c13 * B.px - M.c14 * B.my;
    ctype const c2 = M.c7 * B.pw - M.c8 * B.vy + M.c9 * B.vx + M.c11 * B.py -
                     M.c12 * B.px - M.c14 * B.mz;
    ctype const c3 = M.c5 * B.vx + M.c6 * B.vy + M.c7 * B.vz + M.c11 * B.mx +
                     M.c12 * B.my + M.c13 * B.mz;
    ctype const c4 = M.c5 * B.px + M.c6 * B.py + M.c7 * B.pz + M.c8 * B.mx + M.c9 * B.my +
                     M.c10 * B.mz;
    ctype const c5 = -M.c15 * B.mx;
    ctype const c6 = -M.c15 * B.my;
    ctype const c7 = -M.c15 * B.mz;
    ctype const c8 = -M.c15 * B.px;
    ctype const c9 = -M.c15 * B.py;
    ctype const c10 = -M.c15 * B.pz;
    ctype const c11 = -M.c15 * B.vx;
    ctype const c12 = -M.c15 * B.vy;
    ctype const c13 = -M.c15 * B.vz;
    ctype const c14 = M.c15 * B.pw;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(BiVec3dc<T> const& B,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c4 + B.my * M.c2 - B.mz * M.c1 + B.px * M.c3;
    ctype const c1 = -B.vy * M.c4 - B.mx * M.c2 + B.mz * M.c0 + B.py * M.c3;
    ctype const c2 = -B.vz * M.c4 + B.mx * M.c1 - B.my * M.c0 + B.pz * M.c3;
    ctype const c3 = -B.vx * M.c0 - B.vy * M.c1 - B.vz * M.c2 + B.pw * M.c3;
    ctype const c4 = B.px * M.c0 + B.py * M.c1 + B.pz * M.c2 - B.pw * M.c4;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& M,
                                                         Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y + M.c2 * v.z - M.c3 * v.u - M.c4 * v.w;
    ctype const c1 = -M.c5 * v.w - M.c12 * v.z + M.c13 * v.y;
    ctype const c2 = -M.c6 * v.w + M.c11 * v.z - M.c13 * v.x;
    ctype const c3 = -M.c7 * v.w - M.c11 * v.y + M.c12 * v.x;
    ctype const c4 = -M.c8 * v.w - M.c11 * v.u - M.c14 * v.x;
    ctype const c5 = -M.c9 * v.w - M.c12 * v.u - M.c14 * v.y;
    ctype const c6 = -M.c10 * v.w - M.c13 * v.u - M.c14 * v.z;
    ctype const c7 = -M.c5 * v.u + M.c9 * v.z - M.c10 * v.y;
    ctype const c8 = -M.c6 * v.u - M.c8 * v.z + M.c10 * v.x;
    ctype const c9 = -M.c7 * v.u + M.c8 * v.y - M.c9 * v.x;
    ctype const c10 = -M.c5 * v.x - M.c6 * v.y - M.c7 * v.z;
    ctype const c11 = M.c15 * v.x;
    ctype const c12 = M.c15 * v.y;
    ctype const c13 = M.c15 * v.z;
    ctype const c14 = -M.c15 * v.u;
    ctype const c15 = -M.c15 * v.w;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator>>(Vec3dc<T> const& v,
                                                         MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 + v.z * M.c2 - v.w * M.c4 - v.u * M.c3;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc r_contract :: r_contract(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator>>(MVec3dc_U<T> const& M,
                                                         Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc r_contract :: r_contract(s,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dc<T>, [[maybe_unused]] MVec3dc_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps1,
                                                         PScalar3dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga3dc r_contract :: r_contract(ps,quadvec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * Q.x;
    ctype const c1 = -ctype(ps) * Q.y;
    ctype const c2 = -ctype(ps) * Q.z;
    ctype const c3 = ctype(ps) * Q.u;
    ctype const c4 = ctype(ps) * Q.w;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc r_contract :: r_contract(quadvec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] QuadVec3dc<T> const&, [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(ps,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * t.px;
    ctype const c1 = ctype(ps) * t.py;
    ctype const c2 = ctype(ps) * t.pz;
    ctype const c3 = ctype(ps) * t.vx;
    ctype const c4 = ctype(ps) * t.vy;
    ctype const c5 = ctype(ps) * t.vz;
    ctype const c6 = ctype(ps) * t.mx;
    ctype const c7 = ctype(ps) * t.my;
    ctype const c8 = ctype(ps) * t.mz;
    ctype const c9 = -ctype(ps) * t.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc r_contract :: r_contract(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] TriVec3dc<T> const&, [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(ps,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                         BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.mx;
    ctype const c1 = -ctype(ps) * B.my;
    ctype const c2 = -ctype(ps) * B.mz;
    ctype const c3 = -ctype(ps) * B.px;
    ctype const c4 = -ctype(ps) * B.py;
    ctype const c5 = -ctype(ps) * B.pz;
    ctype const c6 = -ctype(ps) * B.vx;
    ctype const c7 = -ctype(ps) * B.vy;
    ctype const c8 = -ctype(ps) * B.vz;
    ctype const c9 = ctype(ps) * B.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc r_contract :: r_contract(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3dc<T> const&, [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(ps,vec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                          Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = -ctype(ps) * v.u;
    ctype const c4 = -ctype(ps) * v.w;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc r_contract :: r_contract(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dc<T> const&, [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> operator>>(PScalar3dc<T> ps,
                                                          Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(ps) * ctype(s));
}

// cga3dc r_contract :: r_contract(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar3dc<T>,
                                                         [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(quadvec,quadvec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q1,
                                                         QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-Q1.x * Q2.x - Q1.y * Q2.y - Q1.z * Q2.z + Q1.w * Q2.u +
                            Q1.u * Q2.w);
}

// cga3dc r_contract :: r_contract(quadvec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * t.vz - Q.z * t.vy - Q.w * t.px + Q.u * t.mx;
    ctype const c1 = -Q.x * t.vz + Q.z * t.vx - Q.w * t.py + Q.u * t.my;
    ctype const c2 = Q.x * t.vy - Q.y * t.vx - Q.w * t.pz + Q.u * t.mz;
    ctype const c3 = Q.x * t.px + Q.y * t.py + Q.z * t.pz - Q.u * t.pw;
    ctype const c4 = -Q.x * t.mx - Q.y * t.my - Q.z * t.mz + Q.w * t.pw;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc r_contract :: r_contract(trivec,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] TriVec3dc<T> const&, [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(quadvec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.y * B.vz + Q.z * B.vy - Q.u * B.mx;
    ctype const c1 = Q.x * B.vz - Q.z * B.vx - Q.u * B.my;
    ctype const c2 = -Q.x * B.vy + Q.y * B.vx - Q.u * B.mz;
    ctype const c3 = -Q.x * B.pw - Q.w * B.vx - Q.u * B.px;
    ctype const c4 = -Q.y * B.pw - Q.w * B.vy - Q.u * B.py;
    ctype const c5 = -Q.z * B.pw - Q.w * B.vz - Q.u * B.pz;
    ctype const c6 = Q.y * B.pz - Q.z * B.py - Q.w * B.mx;
    ctype const c7 = -Q.x * B.pz + Q.z * B.px - Q.w * B.my;
    ctype const c8 = Q.x * B.py - Q.y * B.px - Q.w * B.mz;
    ctype const c9 = Q.x * B.mx + Q.y * B.my + Q.z * B.mz;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc r_contract :: r_contract(bivec,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3dc<T> const&, [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(quadvec,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q,
                                                         Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.y * v.z + Q.z * v.y;
    ctype const c1 = Q.x * v.z - Q.z * v.x;
    ctype const c2 = -Q.x * v.y + Q.y * v.x;
    ctype const c3 = -Q.x * v.u - Q.w * v.x;
    ctype const c4 = -Q.y * v.u - Q.w * v.y;
    ctype const c5 = -Q.z * v.u - Q.w * v.z;
    ctype const c6 = Q.x * v.w + Q.u * v.x;
    ctype const c7 = Q.y * v.w + Q.u * v.y;
    ctype const c8 = Q.z * v.w + Q.u * v.z;
    ctype const c9 = Q.w * v.w - Q.u * v.u;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc r_contract :: r_contract(vec,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dc<T> const&, [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(quadvec,s) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator>>(QuadVec3dc<T> const& Q,
                                                          Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * ctype(s);
    ctype const c1 = Q.y * ctype(s);
    ctype const c2 = Q.z * ctype(s);
    ctype const c3 = Q.w * ctype(s);
    ctype const c4 = Q.u * ctype(s);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc r_contract :: r_contract(s,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dc<T>, [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator>>(TriVec3dc<T> const& t1,
                                                         TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(
        -t1.vx * t2.vx - t1.vy * t2.vy - t1.vz * t2.vz - t1.mx * t2.px - t1.my * t2.py -
        t1.mz * t2.pz - t1.px * t2.mx - t1.py * t2.my - t1.pz * t2.mz + t1.pw * t2.pw);
}

// cga3dc r_contract :: r_contract(trivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator>>(TriVec3dc<T> const& t,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        t.vx * B.pw - t.my * B.vz + t.mz * B.vy + t.py * B.pz - t.pz * B.py - t.pw * B.mx;
    ctype const c1 =
        t.vy * B.pw + t.mx * B.vz - t.mz * B.vx - t.px * B.pz + t.pz * B.px - t.pw * B.my;
    ctype const c2 =
        t.vz * B.pw - t.mx * B.vy + t.my * B.vx + t.px * B.py - t.py * B.px - t.pw * B.mz;
    ctype const c3 =
        t.vx * B.vx + t.vy * B.vy + t.vz * B.vz + t.px * B.mx + t.py * B.my + t.pz * B.mz;
    ctype const c4 =
        t.vx * B.px + t.vy * B.py + t.vz * B.pz + t.mx * B.mx + t.my * B.my + t.mz * B.mz;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc r_contract :: r_contract(bivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3dc<T> const&, [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(trivec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator>>(TriVec3dc<T> const& t,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * v.w - t.py * v.z + t.pz * v.y;
    ctype const c1 = -t.vy * v.w + t.px * v.z - t.pz * v.x;
    ctype const c2 = -t.vz * v.w - t.px * v.y + t.py * v.x;
    ctype const c3 = -t.mx * v.w - t.px * v.u - t.pw * v.x;
    ctype const c4 = -t.my * v.w - t.py * v.u - t.pw * v.y;
    ctype const c5 = -t.mz * v.w - t.pz * v.u - t.pw * v.z;
    ctype const c6 = -t.vx * v.u + t.my * v.z - t.mz * v.y;
    ctype const c7 = -t.vy * v.u - t.mx * v.z + t.mz * v.x;
    ctype const c8 = -t.vz * v.u + t.mx * v.y - t.my * v.x;
    ctype const c9 = -t.vx * v.x - t.vy * v.y - t.vz * v.z;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc r_contract :: r_contract(vec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dc<T> const&, [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator>>(TriVec3dc<T> const& t,
                                                         Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * ctype(s);
    ctype const c1 = t.vy * ctype(s);
    ctype const c2 = t.vz * ctype(s);
    ctype const c3 = t.mx * ctype(s);
    ctype const c4 = t.my * ctype(s);
    ctype const c5 = t.mz * ctype(s);
    ctype const c6 = t.px * ctype(s);
    ctype const c7 = t.py * ctype(s);
    ctype const c8 = t.pz * ctype(s);
    ctype const c9 = t.pw * ctype(s);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc r_contract :: r_contract(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dc<T>, [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator>>(BiVec3dc<T> const& B1,
                                                         BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(
        B1.vx * B2.px + B1.vy * B2.py + B1.vz * B2.pz + B1.mx * B2.mx + B1.my * B2.my +
        B1.mz * B2.mz + B1.px * B2.vx + B1.py * B2.vy + B1.pz * B2.vz - B1.pw * B2.pw);
}

// cga3dc r_contract :: r_contract(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator>>(BiVec3dc<T> const& B,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * v.u + B.my * v.z - B.mz * v.y + B.px * v.w;
    ctype const c1 = -B.vy * v.u - B.mx * v.z + B.mz * v.x + B.py * v.w;
    ctype const c2 = -B.vz * v.u + B.mx * v.y - B.my * v.x + B.pz * v.w;
    ctype const c3 = -B.vx * v.x - B.vy * v.y - B.vz * v.z + B.pw * v.w;
    ctype const c4 = B.px * v.x + B.py * v.y + B.pz * v.z - B.pw * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc r_contract :: r_contract(vec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dc<T> const&, [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator>>(BiVec3dc<T> const& B,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    ctype const c6 = B.px * ctype(s);
    ctype const c7 = B.py * ctype(s);
    ctype const c8 = B.pz * ctype(s);
    ctype const c9 = B.pw * ctype(s);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc r_contract :: r_contract(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dc<T>, [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator>>(Vec3dc<T> const& v1,
                                                         Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z - v1.w * v2.u -
                            v1.u * v2.w);
}

// cga3dc r_contract :: r_contract(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator>>(Vec3dc<T> const& v, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    ctype const c4 = v.u * ctype(s);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc r_contract :: r_contract(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dc<T>, [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc r_contract :: r_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator>>(Scalar3dc<T> s1, Scalar3dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

// cga3dc cmt :: cmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> cmt(MVec3dc<T> const& A, MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c2 * B.c11 + A.c3 * B.c10 + A.c4 * B.c12 - A.c5 * B.c6 +
                     A.c6 * B.c5 - A.c10 * B.c3 + A.c11 * B.c2 - A.c12 * B.c4 +
                     A.c17 * B.c28 - A.c18 * B.c27 - A.c19 * B.c30 + A.c22 * B.c29 +
                     A.c27 * B.c18 - A.c28 * B.c17 - A.c29 * B.c22 + A.c30 * B.c19;
    ctype const c2 = A.c1 * B.c11 - A.c3 * B.c9 + A.c4 * B.c13 - A.c5 * B.c7 +
                     A.c7 * B.c5 + A.c9 * B.c3 - A.c11 * B.c1 - A.c13 * B.c4 -
                     A.c16 * B.c28 + A.c18 * B.c26 - A.c20 * B.c30 + A.c23 * B.c29 -
                     A.c26 * B.c18 + A.c28 * B.c16 - A.c29 * B.c23 + A.c30 * B.c20;
    ctype const c3 = -A.c1 * B.c10 + A.c2 * B.c9 + A.c4 * B.c14 - A.c5 * B.c8 +
                     A.c8 * B.c5 - A.c9 * B.c2 + A.c10 * B.c1 - A.c14 * B.c4 +
                     A.c16 * B.c27 - A.c17 * B.c26 - A.c21 * B.c30 + A.c24 * B.c29 +
                     A.c26 * B.c17 - A.c27 * B.c16 - A.c29 * B.c24 + A.c30 * B.c21;
    ctype const c4 = -A.c1 * B.c6 - A.c2 * B.c7 - A.c3 * B.c8 + A.c4 * B.c15 +
                     A.c6 * B.c1 + A.c7 * B.c2 + A.c8 * B.c3 - A.c15 * B.c4 -
                     A.c22 * B.c26 - A.c23 * B.c27 - A.c24 * B.c28 + A.c25 * B.c30 +
                     A.c26 * B.c22 + A.c27 * B.c23 + A.c28 * B.c24 - A.c30 * B.c25;
    ctype const c5 = A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 - A.c5 * B.c15 -
                     A.c12 * B.c1 - A.c13 * B.c2 - A.c14 * B.c3 + A.c15 * B.c5 +
                     A.c19 * B.c26 + A.c20 * B.c27 + A.c21 * B.c28 - A.c25 * B.c29 -
                     A.c26 * B.c19 - A.c27 * B.c20 - A.c28 * B.c21 + A.c29 * B.c25;
    ctype const c6 = -A.c1 * B.c4 + A.c4 * B.c1 + A.c6 * B.c15 - A.c7 * B.c11 +
                     A.c8 * B.c10 - A.c10 * B.c8 + A.c11 * B.c7 - A.c15 * B.c6 +
                     A.c17 * B.c24 - A.c18 * B.c23 + A.c22 * B.c25 + A.c23 * B.c18 -
                     A.c24 * B.c17 - A.c25 * B.c22 - A.c26 * B.c30 + A.c30 * B.c26;
    ctype const c7 = -A.c2 * B.c4 + A.c4 * B.c2 + A.c6 * B.c11 + A.c7 * B.c15 -
                     A.c8 * B.c9 + A.c9 * B.c8 - A.c11 * B.c6 - A.c15 * B.c7 -
                     A.c16 * B.c24 + A.c18 * B.c22 - A.c22 * B.c18 + A.c23 * B.c25 +
                     A.c24 * B.c16 - A.c25 * B.c23 - A.c27 * B.c30 + A.c30 * B.c27;
    ctype const c8 = -A.c3 * B.c4 + A.c4 * B.c3 - A.c6 * B.c10 + A.c7 * B.c9 +
                     A.c8 * B.c15 - A.c9 * B.c7 + A.c10 * B.c6 - A.c15 * B.c8 +
                     A.c16 * B.c23 - A.c17 * B.c22 + A.c22 * B.c17 - A.c23 * B.c16 +
                     A.c24 * B.c25 - A.c25 * B.c24 - A.c28 * B.c30 + A.c30 * B.c28;
    ctype const c9 = A.c2 * B.c3 - A.c3 * B.c2 - A.c7 * B.c14 + A.c8 * B.c13 -
                     A.c10 * B.c11 + A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c7 +
                     A.c17 * B.c18 - A.c18 * B.c17 + A.c20 * B.c24 - A.c21 * B.c23 +
                     A.c23 * B.c21 - A.c24 * B.c20 - A.c27 * B.c28 + A.c28 * B.c27;
    ctype const c10 = -A.c1 * B.c3 + A.c3 * B.c1 + A.c6 * B.c14 - A.c8 * B.c12 +
                      A.c9 * B.c11 - A.c11 * B.c9 + A.c12 * B.c8 - A.c14 * B.c6 -
                      A.c16 * B.c18 + A.c18 * B.c16 - A.c19 * B.c24 + A.c21 * B.c22 -
                      A.c22 * B.c21 + A.c24 * B.c19 + A.c26 * B.c28 - A.c28 * B.c26;
    ctype const c11 = A.c1 * B.c2 - A.c2 * B.c1 - A.c6 * B.c13 + A.c7 * B.c12 -
                      A.c9 * B.c10 + A.c10 * B.c9 - A.c12 * B.c7 + A.c13 * B.c6 +
                      A.c16 * B.c17 - A.c17 * B.c16 + A.c19 * B.c23 - A.c20 * B.c22 +
                      A.c22 * B.c20 - A.c23 * B.c19 - A.c26 * B.c27 + A.c27 * B.c26;
    ctype const c12 = A.c1 * B.c5 - A.c5 * B.c1 - A.c10 * B.c14 + A.c11 * B.c13 -
                      A.c12 * B.c15 - A.c13 * B.c11 + A.c14 * B.c10 + A.c15 * B.c12 +
                      A.c17 * B.c21 - A.c18 * B.c20 - A.c19 * B.c25 + A.c20 * B.c18 -
                      A.c21 * B.c17 + A.c25 * B.c19 + A.c26 * B.c29 - A.c29 * B.c26;
    ctype const c13 = A.c2 * B.c5 - A.c5 * B.c2 + A.c9 * B.c14 - A.c11 * B.c12 +
                      A.c12 * B.c11 - A.c13 * B.c15 - A.c14 * B.c9 + A.c15 * B.c13 -
                      A.c16 * B.c21 + A.c18 * B.c19 - A.c19 * B.c18 - A.c20 * B.c25 +
                      A.c21 * B.c16 + A.c25 * B.c20 + A.c27 * B.c29 - A.c29 * B.c27;
    ctype const c14 = A.c3 * B.c5 - A.c5 * B.c3 - A.c9 * B.c13 + A.c10 * B.c12 -
                      A.c12 * B.c10 + A.c13 * B.c9 - A.c14 * B.c15 + A.c15 * B.c14 +
                      A.c16 * B.c20 - A.c17 * B.c19 + A.c19 * B.c17 - A.c20 * B.c16 -
                      A.c21 * B.c25 + A.c25 * B.c21 + A.c28 * B.c29 - A.c29 * B.c28;
    ctype const c15 = A.c4 * B.c5 - A.c5 * B.c4 + A.c6 * B.c12 + A.c7 * B.c13 +
                      A.c8 * B.c14 - A.c12 * B.c6 - A.c13 * B.c7 - A.c14 * B.c8 +
                      A.c19 * B.c22 + A.c20 * B.c23 + A.c21 * B.c24 - A.c22 * B.c19 -
                      A.c23 * B.c20 - A.c24 * B.c21 + A.c29 * B.c30 - A.c30 * B.c29;
    ctype const c16 = A.c2 * B.c28 - A.c3 * B.c27 - A.c7 * B.c21 + A.c8 * B.c20 -
                      A.c10 * B.c18 + A.c11 * B.c17 - A.c13 * B.c24 + A.c14 * B.c23 -
                      A.c17 * B.c11 + A.c18 * B.c10 - A.c20 * B.c8 + A.c21 * B.c7 -
                      A.c23 * B.c14 + A.c24 * B.c13 + A.c27 * B.c3 - A.c28 * B.c2;
    ctype const c17 = -A.c1 * B.c28 + A.c3 * B.c26 + A.c6 * B.c21 - A.c8 * B.c19 +
                      A.c9 * B.c18 - A.c11 * B.c16 + A.c12 * B.c24 - A.c14 * B.c22 +
                      A.c16 * B.c11 - A.c18 * B.c9 + A.c19 * B.c8 - A.c21 * B.c6 +
                      A.c22 * B.c14 - A.c24 * B.c12 - A.c26 * B.c3 + A.c28 * B.c1;
    ctype const c18 = A.c1 * B.c27 - A.c2 * B.c26 - A.c6 * B.c20 + A.c7 * B.c19 -
                      A.c9 * B.c17 + A.c10 * B.c16 - A.c12 * B.c23 + A.c13 * B.c22 -
                      A.c16 * B.c10 + A.c17 * B.c9 - A.c19 * B.c7 + A.c20 * B.c6 -
                      A.c22 * B.c13 + A.c23 * B.c12 + A.c26 * B.c2 - A.c27 * B.c1;
    ctype const c19 = -A.c1 * B.c29 - A.c5 * B.c26 - A.c10 * B.c21 + A.c11 * B.c20 +
                      A.c12 * B.c25 - A.c13 * B.c18 + A.c14 * B.c17 + A.c15 * B.c19 -
                      A.c17 * B.c14 + A.c18 * B.c13 - A.c19 * B.c15 - A.c20 * B.c11 +
                      A.c21 * B.c10 - A.c25 * B.c12 + A.c26 * B.c5 + A.c29 * B.c1;
    ctype const c20 = -A.c2 * B.c29 - A.c5 * B.c27 + A.c9 * B.c21 - A.c11 * B.c19 +
                      A.c12 * B.c18 + A.c13 * B.c25 - A.c14 * B.c16 + A.c15 * B.c20 +
                      A.c16 * B.c14 - A.c18 * B.c12 + A.c19 * B.c11 - A.c20 * B.c15 -
                      A.c21 * B.c9 - A.c25 * B.c13 + A.c27 * B.c5 + A.c29 * B.c2;
    ctype const c21 = -A.c3 * B.c29 - A.c5 * B.c28 - A.c9 * B.c20 + A.c10 * B.c19 -
                      A.c12 * B.c17 + A.c13 * B.c16 + A.c14 * B.c25 + A.c15 * B.c21 -
                      A.c16 * B.c13 + A.c17 * B.c12 - A.c19 * B.c10 + A.c20 * B.c9 -
                      A.c21 * B.c15 - A.c25 * B.c14 + A.c28 * B.c5 + A.c29 * B.c3;
    ctype const c22 = A.c1 * B.c30 + A.c4 * B.c26 - A.c6 * B.c25 - A.c7 * B.c18 +
                      A.c8 * B.c17 - A.c10 * B.c24 + A.c11 * B.c23 - A.c15 * B.c22 -
                      A.c17 * B.c8 + A.c18 * B.c7 + A.c22 * B.c15 - A.c23 * B.c11 +
                      A.c24 * B.c10 + A.c25 * B.c6 - A.c26 * B.c4 - A.c30 * B.c1;
    ctype const c23 = A.c2 * B.c30 + A.c4 * B.c27 + A.c6 * B.c18 - A.c7 * B.c25 -
                      A.c8 * B.c16 + A.c9 * B.c24 - A.c11 * B.c22 - A.c15 * B.c23 +
                      A.c16 * B.c8 - A.c18 * B.c6 + A.c22 * B.c11 + A.c23 * B.c15 -
                      A.c24 * B.c9 + A.c25 * B.c7 - A.c27 * B.c4 - A.c30 * B.c2;
    ctype const c24 = A.c3 * B.c30 + A.c4 * B.c28 - A.c6 * B.c17 + A.c7 * B.c16 -
                      A.c8 * B.c25 - A.c9 * B.c23 + A.c10 * B.c22 - A.c15 * B.c24 -
                      A.c16 * B.c7 + A.c17 * B.c6 - A.c22 * B.c10 + A.c23 * B.c9 +
                      A.c24 * B.c15 + A.c25 * B.c8 - A.c28 * B.c4 - A.c30 * B.c3;
    ctype const c25 = A.c4 * B.c29 - A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                      A.c8 * B.c21 + A.c12 * B.c22 + A.c13 * B.c23 + A.c14 * B.c24 +
                      A.c19 * B.c6 + A.c20 * B.c7 + A.c21 * B.c8 - A.c22 * B.c12 -
                      A.c23 * B.c13 - A.c24 * B.c14 - A.c29 * B.c4 + A.c30 * B.c5;
    ctype const c26 = -A.c2 * B.c18 + A.c3 * B.c17 + A.c4 * B.c19 - A.c5 * B.c22 -
                      A.c6 * B.c29 - A.c10 * B.c28 + A.c11 * B.c27 + A.c12 * B.c30 -
                      A.c17 * B.c3 + A.c18 * B.c2 - A.c19 * B.c4 + A.c22 * B.c5 -
                      A.c27 * B.c11 + A.c28 * B.c10 + A.c29 * B.c6 - A.c30 * B.c12;
    ctype const c27 = A.c1 * B.c18 - A.c3 * B.c16 + A.c4 * B.c20 - A.c5 * B.c23 -
                      A.c7 * B.c29 + A.c9 * B.c28 - A.c11 * B.c26 + A.c13 * B.c30 +
                      A.c16 * B.c3 - A.c18 * B.c1 - A.c20 * B.c4 + A.c23 * B.c5 +
                      A.c26 * B.c11 - A.c28 * B.c9 + A.c29 * B.c7 - A.c30 * B.c13;
    ctype const c28 = -A.c1 * B.c17 + A.c2 * B.c16 + A.c4 * B.c21 - A.c5 * B.c24 -
                      A.c8 * B.c29 - A.c9 * B.c27 + A.c10 * B.c26 + A.c14 * B.c30 -
                      A.c16 * B.c2 + A.c17 * B.c1 - A.c21 * B.c4 + A.c24 * B.c5 -
                      A.c26 * B.c10 + A.c27 * B.c9 + A.c29 * B.c8 - A.c30 * B.c14;
    ctype const c29 = -A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 - A.c5 * B.c25 +
                      A.c12 * B.c26 + A.c13 * B.c27 + A.c14 * B.c28 + A.c15 * B.c29 +
                      A.c19 * B.c1 + A.c20 * B.c2 + A.c21 * B.c3 + A.c25 * B.c5 -
                      A.c26 * B.c12 - A.c27 * B.c13 - A.c28 * B.c14 - A.c29 * B.c15;
    ctype const c30 = A.c1 * B.c22 + A.c2 * B.c23 + A.c3 * B.c24 + A.c4 * B.c25 -
                      A.c6 * B.c26 - A.c7 * B.c27 - A.c8 * B.c28 - A.c15 * B.c30 -
                      A.c22 * B.c1 - A.c23 * B.c2 - A.c24 * B.c3 - A.c25 * B.c4 +
                      A.c26 * B.c6 + A.c27 * B.c7 + A.c28 * B.c8 + A.c30 * B.c15;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc cmt :: cmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(ps,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(quadvec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] QuadVec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] TriVec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec3dc<T> const&,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dc<T>,
                                                  [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dc<T>,
                                                  [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(quadvec,quadvec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> cmt(QuadVec3dc<T> const& Q1,
                                                 QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q1.x * Q2.u + Q1.u * Q2.x;
    ctype const c1 = -Q1.y * Q2.u + Q1.u * Q2.y;
    ctype const c2 = -Q1.z * Q2.u + Q1.u * Q2.z;
    ctype const c3 = -Q1.y * Q2.z + Q1.z * Q2.y;
    ctype const c4 = Q1.x * Q2.z - Q1.z * Q2.x;
    ctype const c5 = -Q1.x * Q2.y + Q1.y * Q2.x;
    ctype const c6 = Q1.x * Q2.w - Q1.w * Q2.x;
    ctype const c7 = Q1.y * Q2.w - Q1.w * Q2.y;
    ctype const c8 = Q1.z * Q2.w - Q1.w * Q2.z;
    ctype const c9 = Q1.w * Q2.u - Q1.u * Q2.w;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(quadvec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> cmt(QuadVec3dc<T> const& Q,
                                               TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * t.vz - Q.z * t.vy - Q.w * t.px + Q.u * t.mx;
    ctype const c1 = -Q.x * t.vz + Q.z * t.vx - Q.w * t.py + Q.u * t.my;
    ctype const c2 = Q.x * t.vy - Q.y * t.vx - Q.w * t.pz + Q.u * t.mz;
    ctype const c3 = Q.x * t.px + Q.y * t.py + Q.z * t.pz - Q.u * t.pw;
    ctype const c4 = -Q.x * t.mx - Q.y * t.my - Q.z * t.mz + Q.w * t.pw;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(trivec,quadvec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> cmt(TriVec3dc<T> const& t,
                                               QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vy * Q.z - t.vz * Q.y - t.mx * Q.u + t.px * Q.w;
    ctype const c1 = -t.vx * Q.z + t.vz * Q.x - t.my * Q.u + t.py * Q.w;
    ctype const c2 = t.vx * Q.y - t.vy * Q.x - t.mz * Q.u + t.pz * Q.w;
    ctype const c3 = -t.px * Q.x - t.py * Q.y - t.pz * Q.z + t.pw * Q.u;
    ctype const c4 = t.mx * Q.x + t.my * Q.y + t.mz * Q.z - t.pw * Q.w;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(quadvec,bivec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> cmt(QuadVec3dc<T> const& Q,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.y * B.mz + Q.z * B.my + Q.w * B.vx - Q.u * B.px;
    ctype const c1 = Q.x * B.mz - Q.z * B.mx + Q.w * B.vy - Q.u * B.py;
    ctype const c2 = -Q.x * B.my + Q.y * B.mx + Q.w * B.vz - Q.u * B.pz;
    ctype const c3 = -Q.x * B.px - Q.y * B.py - Q.z * B.pz - Q.w * B.pw;
    ctype const c4 = Q.x * B.vx + Q.y * B.vy + Q.z * B.vz + Q.u * B.pw;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(bivec,quadvec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> cmt(BiVec3dc<T> const& B,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * Q.w - B.my * Q.z + B.mz * Q.y + B.px * Q.u;
    ctype const c1 = -B.vy * Q.w + B.mx * Q.z - B.mz * Q.x + B.py * Q.u;
    ctype const c2 = -B.vz * Q.w - B.mx * Q.y + B.my * Q.x + B.pz * Q.u;
    ctype const c3 = B.px * Q.x + B.py * Q.y + B.pz * Q.z + B.pw * Q.w;
    ctype const c4 = -B.vx * Q.x - B.vy * Q.y - B.vz * Q.z - B.pw * Q.u;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(quadvec,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> cmt(QuadVec3dc<T> const& Q,
                                                  Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * v.z - Q.z * v.y;
    ctype const c1 = -Q.x * v.z + Q.z * v.x;
    ctype const c2 = Q.x * v.y - Q.y * v.x;
    ctype const c3 = Q.x * v.u + Q.w * v.x;
    ctype const c4 = Q.y * v.u + Q.w * v.y;
    ctype const c5 = Q.z * v.u + Q.w * v.z;
    ctype const c6 = -Q.x * v.w - Q.u * v.x;
    ctype const c7 = -Q.y * v.w - Q.u * v.y;
    ctype const c8 = -Q.z * v.w - Q.u * v.z;
    ctype const c9 = -Q.w * v.w + Q.u * v.u;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(vec,quadvec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> cmt(Vec3dc<T> const& v,
                                                  QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * Q.z - v.z * Q.y;
    ctype const c1 = -v.x * Q.z + v.z * Q.x;
    ctype const c2 = v.x * Q.y - v.y * Q.x;
    ctype const c3 = -v.x * Q.w - v.u * Q.x;
    ctype const c4 = -v.y * Q.w - v.u * Q.y;
    ctype const c5 = -v.z * Q.w - v.u * Q.z;
    ctype const c6 = v.x * Q.u + v.w * Q.x;
    ctype const c7 = v.y * Q.u + v.w * Q.y;
    ctype const c8 = v.z * Q.u + v.w * Q.z;
    ctype const c9 = v.w * Q.w - v.u * Q.u;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(quadvec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] QuadVec3dc<T> const&,
                                                  [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(s,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dc<T>,
                                                  [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(trivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> cmt(TriVec3dc<T> const& t1,
                                                 TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.vy * t2.pz - t1.vz * t2.py + t1.px * t2.pw + t1.py * t2.vz -
                     t1.pz * t2.vy - t1.pw * t2.px;
    ctype const c1 = -t1.vx * t2.pz + t1.vz * t2.px - t1.px * t2.vz + t1.py * t2.pw +
                     t1.pz * t2.vx - t1.pw * t2.py;
    ctype const c2 = t1.vx * t2.py - t1.vy * t2.px + t1.px * t2.vy - t1.py * t2.vx +
                     t1.pz * t2.pw - t1.pw * t2.pz;
    ctype const c3 = t1.vy * t2.vz - t1.vz * t2.vy + t1.my * t2.pz - t1.mz * t2.py +
                     t1.py * t2.mz - t1.pz * t2.my;
    ctype const c4 = -t1.vx * t2.vz + t1.vz * t2.vx - t1.mx * t2.pz + t1.mz * t2.px -
                     t1.px * t2.mz + t1.pz * t2.mx;
    ctype const c5 = t1.vx * t2.vy - t1.vy * t2.vx + t1.mx * t2.py - t1.my * t2.px +
                     t1.px * t2.my - t1.py * t2.mx;
    ctype const c6 = t1.vy * t2.mz - t1.vz * t2.my - t1.mx * t2.pw + t1.my * t2.vz -
                     t1.mz * t2.vy + t1.pw * t2.mx;
    ctype const c7 = -t1.vx * t2.mz + t1.vz * t2.mx - t1.mx * t2.vz - t1.my * t2.pw +
                     t1.mz * t2.vx + t1.pw * t2.my;
    ctype const c8 = t1.vx * t2.my - t1.vy * t2.mx + t1.mx * t2.vy - t1.my * t2.vx -
                     t1.mz * t2.pw + t1.pw * t2.mz;
    ctype const c9 = t1.mx * t2.px + t1.my * t2.py + t1.mz * t2.pz - t1.px * t2.mx -
                     t1.py * t2.my - t1.pz * t2.mz;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(trivec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> cmt(TriVec3dc<T> const& t,
                                                  BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vy * B.mz + t.vz * B.my - t.my * B.vz + t.mz * B.vy -
                     t.py * B.pz + t.pz * B.py;
    ctype const c1 =
        t.vx * B.mz - t.vz * B.mx + t.mx * B.vz - t.mz * B.vx + t.px * B.pz - t.pz * B.px;
    ctype const c2 = -t.vx * B.my + t.vy * B.mx - t.mx * B.vy + t.my * B.vx -
                     t.px * B.py + t.py * B.px;
    ctype const c3 = -t.vy * B.pz + t.vz * B.py - t.mx * B.pw - t.my * B.mz +
                     t.mz * B.my - t.pw * B.px;
    ctype const c4 =
        t.vx * B.pz - t.vz * B.px + t.mx * B.mz - t.my * B.pw - t.mz * B.mx - t.pw * B.py;
    ctype const c5 = -t.vx * B.py + t.vy * B.px - t.mx * B.my + t.my * B.mx -
                     t.mz * B.pw - t.pw * B.pz;
    ctype const c6 = -t.vy * B.vz + t.vz * B.vy + t.px * B.pw - t.py * B.mz +
                     t.pz * B.my + t.pw * B.vx;
    ctype const c7 =
        t.vx * B.vz - t.vz * B.vx + t.px * B.mz + t.py * B.pw - t.pz * B.mx + t.pw * B.vy;
    ctype const c8 = -t.vx * B.vy + t.vy * B.vx - t.px * B.my + t.py * B.mx +
                     t.pz * B.pw + t.pw * B.vz;
    ctype const c9 =
        t.mx * B.vx + t.my * B.vy + t.mz * B.vz - t.px * B.px - t.py * B.py - t.pz * B.pz;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(bivec,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> cmt(BiVec3dc<T> const& B,
                                                  TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * t.mz + B.vz * t.my - B.my * t.vz + B.mz * t.vy -
                     B.py * t.pz + B.pz * t.py;
    ctype const c1 =
        B.vx * t.mz - B.vz * t.mx + B.mx * t.vz - B.mz * t.vx + B.px * t.pz - B.pz * t.px;
    ctype const c2 = -B.vx * t.my + B.vy * t.mx - B.mx * t.vy + B.my * t.vx -
                     B.px * t.py + B.py * t.px;
    ctype const c3 = -B.my * t.mz + B.mz * t.my + B.px * t.pw - B.py * t.vz +
                     B.pz * t.vy + B.pw * t.mx;
    ctype const c4 =
        B.mx * t.mz - B.mz * t.mx + B.px * t.vz + B.py * t.pw - B.pz * t.vx + B.pw * t.my;
    ctype const c5 = -B.mx * t.my + B.my * t.mx - B.px * t.vy + B.py * t.vx +
                     B.pz * t.pw + B.pw * t.mz;
    ctype const c6 = -B.vx * t.pw - B.vy * t.vz + B.vz * t.vy - B.my * t.pz +
                     B.mz * t.py - B.pw * t.px;
    ctype const c7 =
        B.vx * t.vz - B.vy * t.pw - B.vz * t.vx + B.mx * t.pz - B.mz * t.px - B.pw * t.py;
    ctype const c8 = -B.vx * t.vy + B.vy * t.vx - B.vz * t.pw - B.mx * t.py +
                     B.my * t.px - B.pw * t.pz;
    ctype const c9 = -B.vx * t.mx - B.vy * t.my - B.vz * t.mz + B.px * t.px +
                     B.py * t.py + B.pz * t.pz;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(trivec,vec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> cmt(TriVec3dc<T> const& t,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vy * v.z + t.vz * v.y - t.mx * v.w + t.px * v.u;
    ctype const c1 = t.vx * v.z - t.vz * v.x - t.my * v.w + t.py * v.u;
    ctype const c2 = -t.vx * v.y + t.vy * v.x - t.mz * v.w + t.pz * v.u;
    ctype const c3 = t.mx * v.x + t.my * v.y + t.mz * v.z + t.pw * v.u;
    ctype const c4 = -t.px * v.x - t.py * v.y - t.pz * v.z - t.pw * v.w;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(vec,trivec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> cmt(Vec3dc<T> const& v,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * t.vz + v.z * t.vy + v.w * t.mx - v.u * t.px;
    ctype const c1 = v.x * t.vz - v.z * t.vx + v.w * t.my - v.u * t.py;
    ctype const c2 = -v.x * t.vy + v.y * t.vx + v.w * t.mz - v.u * t.pz;
    ctype const c3 = -v.x * t.mx - v.y * t.my - v.z * t.mz - v.u * t.pw;
    ctype const c4 = v.x * t.px + v.y * t.py + v.z * t.pz + v.w * t.pw;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] TriVec3dc<T> const&,
                                                  [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dc<T>,
                                                  [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(bivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> cmt(BiVec3dc<T> const& B1,
                                                 BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vx * B2.pw - B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz +
                     B1.mz * B2.vy - B1.pw * B2.vx;
    ctype const c1 = B1.vx * B2.mz + B1.vy * B2.pw - B1.vz * B2.mx + B1.mx * B2.vz -
                     B1.mz * B2.vx - B1.pw * B2.vy;
    ctype const c2 = -B1.vx * B2.my + B1.vy * B2.mx + B1.vz * B2.pw - B1.mx * B2.vy +
                     B1.my * B2.vx - B1.pw * B2.vz;
    ctype const c3 = -B1.vy * B2.pz + B1.vz * B2.py - B1.my * B2.mz + B1.mz * B2.my -
                     B1.py * B2.vz + B1.pz * B2.vy;
    ctype const c4 = B1.vx * B2.pz - B1.vz * B2.px + B1.mx * B2.mz - B1.mz * B2.mx +
                     B1.px * B2.vz - B1.pz * B2.vx;
    ctype const c5 = -B1.vx * B2.py + B1.vy * B2.px - B1.mx * B2.my + B1.my * B2.mx -
                     B1.px * B2.vy + B1.py * B2.vx;
    ctype const c6 = -B1.my * B2.pz + B1.mz * B2.py - B1.px * B2.pw - B1.py * B2.mz +
                     B1.pz * B2.my + B1.pw * B2.px;
    ctype const c7 = B1.mx * B2.pz - B1.mz * B2.px + B1.px * B2.mz - B1.py * B2.pw -
                     B1.pz * B2.mx + B1.pw * B2.py;
    ctype const c8 = -B1.mx * B2.py + B1.my * B2.px - B1.px * B2.my + B1.py * B2.mx -
                     B1.pz * B2.pw + B1.pw * B2.pz;
    ctype const c9 = B1.vx * B2.px + B1.vy * B2.py + B1.vz * B2.pz - B1.px * B2.vx -
                     B1.py * B2.vy - B1.pz * B2.vz;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> cmt(BiVec3dc<T> const& B, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.u - B.my * v.z + B.mz * v.y - B.px * v.w;
    ctype const c1 = B.vy * v.u + B.mx * v.z - B.mz * v.x - B.py * v.w;
    ctype const c2 = B.vz * v.u - B.mx * v.y + B.my * v.x - B.pz * v.w;
    ctype const c3 = B.vx * v.x + B.vy * v.y + B.vz * v.z - B.pw * v.w;
    ctype const c4 = -B.px * v.x - B.py * v.y - B.pz * v.z + B.pw * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> cmt(Vec3dc<T> const& v, BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my + v.w * B.px - v.u * B.vx;
    ctype const c1 = v.x * B.mz - v.z * B.mx + v.w * B.py - v.u * B.vy;
    ctype const c2 = -v.x * B.my + v.y * B.mx + v.w * B.pz - v.u * B.vz;
    ctype const c3 = -v.x * B.vx - v.y * B.vy - v.z * B.vz + v.w * B.pw;
    ctype const c4 = v.x * B.px + v.y * B.py + v.z * B.pz - v.u * B.pw;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc cmt :: cmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec3dc<T> const&,
                                                  [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dc<T>,
                                                  [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> cmt(Vec3dc<T> const& v1, Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c1 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c2 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c3 = v1.y * v2.z - v1.z * v2.y;
    ctype const c4 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c5 = v1.x * v2.y - v1.y * v2.x;
    ctype const c6 = v1.x * v2.u - v1.u * v2.x;
    ctype const c7 = v1.y * v2.u - v1.u * v2.y;
    ctype const c8 = v1.z * v2.u - v1.u * v2.z;
    ctype const c9 = v1.w * v2.u - v1.u * v2.w;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc cmt :: cmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec3dc<T> const&,
                                                  [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dc<T>,
                                                  [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc cmt :: cmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dc<T>,
                                                  [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////

// cga3dc rcmt :: rcmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rcmt(MVec3dc<T> const& A, MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c2 * B.c18 - A.c3 * B.c17 - A.c4 * B.c19 + A.c5 * B.c22 +
                     A.c6 * B.c29 + A.c10 * B.c28 - A.c11 * B.c27 - A.c12 * B.c30 +
                     A.c17 * B.c3 - A.c18 * B.c2 + A.c19 * B.c4 - A.c22 * B.c5 +
                     A.c27 * B.c11 - A.c28 * B.c10 - A.c29 * B.c6 + A.c30 * B.c12;
    ctype const c2 = -A.c1 * B.c18 + A.c3 * B.c16 - A.c4 * B.c20 + A.c5 * B.c23 +
                     A.c7 * B.c29 - A.c9 * B.c28 + A.c11 * B.c26 - A.c13 * B.c30 -
                     A.c16 * B.c3 + A.c18 * B.c1 + A.c20 * B.c4 - A.c23 * B.c5 -
                     A.c26 * B.c11 + A.c28 * B.c9 - A.c29 * B.c7 + A.c30 * B.c13;
    ctype const c3 = A.c1 * B.c17 - A.c2 * B.c16 - A.c4 * B.c21 + A.c5 * B.c24 +
                     A.c8 * B.c29 + A.c9 * B.c27 - A.c10 * B.c26 - A.c14 * B.c30 +
                     A.c16 * B.c2 - A.c17 * B.c1 + A.c21 * B.c4 - A.c24 * B.c5 +
                     A.c26 * B.c10 - A.c27 * B.c9 - A.c29 * B.c8 + A.c30 * B.c14;
    ctype const c4 = A.c1 * B.c22 + A.c2 * B.c23 + A.c3 * B.c24 + A.c4 * B.c25 -
                     A.c6 * B.c26 - A.c7 * B.c27 - A.c8 * B.c28 - A.c15 * B.c30 -
                     A.c22 * B.c1 - A.c23 * B.c2 - A.c24 * B.c3 - A.c25 * B.c4 +
                     A.c26 * B.c6 + A.c27 * B.c7 + A.c28 * B.c8 + A.c30 * B.c15;
    ctype const c5 = -A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 - A.c5 * B.c25 +
                     A.c12 * B.c26 + A.c13 * B.c27 + A.c14 * B.c28 + A.c15 * B.c29 +
                     A.c19 * B.c1 + A.c20 * B.c2 + A.c21 * B.c3 + A.c25 * B.c5 -
                     A.c26 * B.c12 - A.c27 * B.c13 - A.c28 * B.c14 - A.c29 * B.c15;
    ctype const c6 = -A.c1 * B.c30 - A.c4 * B.c26 + A.c6 * B.c25 + A.c7 * B.c18 -
                     A.c8 * B.c17 + A.c10 * B.c24 - A.c11 * B.c23 + A.c15 * B.c22 +
                     A.c17 * B.c8 - A.c18 * B.c7 - A.c22 * B.c15 + A.c23 * B.c11 -
                     A.c24 * B.c10 - A.c25 * B.c6 + A.c26 * B.c4 + A.c30 * B.c1;
    ctype const c7 = -A.c2 * B.c30 - A.c4 * B.c27 - A.c6 * B.c18 + A.c7 * B.c25 +
                     A.c8 * B.c16 - A.c9 * B.c24 + A.c11 * B.c22 + A.c15 * B.c23 -
                     A.c16 * B.c8 + A.c18 * B.c6 - A.c22 * B.c11 - A.c23 * B.c15 +
                     A.c24 * B.c9 - A.c25 * B.c7 + A.c27 * B.c4 + A.c30 * B.c2;
    ctype const c8 = -A.c3 * B.c30 - A.c4 * B.c28 + A.c6 * B.c17 - A.c7 * B.c16 +
                     A.c8 * B.c25 + A.c9 * B.c23 - A.c10 * B.c22 + A.c15 * B.c24 +
                     A.c16 * B.c7 - A.c17 * B.c6 + A.c22 * B.c10 - A.c23 * B.c9 -
                     A.c24 * B.c15 - A.c25 * B.c8 + A.c28 * B.c4 + A.c30 * B.c3;
    ctype const c9 = -A.c2 * B.c28 + A.c3 * B.c27 + A.c7 * B.c21 - A.c8 * B.c20 +
                     A.c10 * B.c18 - A.c11 * B.c17 + A.c13 * B.c24 - A.c14 * B.c23 +
                     A.c17 * B.c11 - A.c18 * B.c10 + A.c20 * B.c8 - A.c21 * B.c7 +
                     A.c23 * B.c14 - A.c24 * B.c13 - A.c27 * B.c3 + A.c28 * B.c2;
    ctype const c10 = A.c1 * B.c28 - A.c3 * B.c26 - A.c6 * B.c21 + A.c8 * B.c19 -
                      A.c9 * B.c18 + A.c11 * B.c16 - A.c12 * B.c24 + A.c14 * B.c22 -
                      A.c16 * B.c11 + A.c18 * B.c9 - A.c19 * B.c8 + A.c21 * B.c6 -
                      A.c22 * B.c14 + A.c24 * B.c12 + A.c26 * B.c3 - A.c28 * B.c1;
    ctype const c11 = -A.c1 * B.c27 + A.c2 * B.c26 + A.c6 * B.c20 - A.c7 * B.c19 +
                      A.c9 * B.c17 - A.c10 * B.c16 + A.c12 * B.c23 - A.c13 * B.c22 +
                      A.c16 * B.c10 - A.c17 * B.c9 + A.c19 * B.c7 - A.c20 * B.c6 +
                      A.c22 * B.c13 - A.c23 * B.c12 - A.c26 * B.c2 + A.c27 * B.c1;
    ctype const c12 = A.c1 * B.c29 + A.c5 * B.c26 + A.c10 * B.c21 - A.c11 * B.c20 -
                      A.c12 * B.c25 + A.c13 * B.c18 - A.c14 * B.c17 - A.c15 * B.c19 +
                      A.c17 * B.c14 - A.c18 * B.c13 + A.c19 * B.c15 + A.c20 * B.c11 -
                      A.c21 * B.c10 + A.c25 * B.c12 - A.c26 * B.c5 - A.c29 * B.c1;
    ctype const c13 = A.c2 * B.c29 + A.c5 * B.c27 - A.c9 * B.c21 + A.c11 * B.c19 -
                      A.c12 * B.c18 - A.c13 * B.c25 + A.c14 * B.c16 - A.c15 * B.c20 -
                      A.c16 * B.c14 + A.c18 * B.c12 - A.c19 * B.c11 + A.c20 * B.c15 +
                      A.c21 * B.c9 + A.c25 * B.c13 - A.c27 * B.c5 - A.c29 * B.c2;
    ctype const c14 = A.c3 * B.c29 + A.c5 * B.c28 + A.c9 * B.c20 - A.c10 * B.c19 +
                      A.c12 * B.c17 - A.c13 * B.c16 - A.c14 * B.c25 - A.c15 * B.c21 +
                      A.c16 * B.c13 - A.c17 * B.c12 + A.c19 * B.c10 - A.c20 * B.c9 +
                      A.c21 * B.c15 + A.c25 * B.c14 - A.c28 * B.c5 - A.c29 * B.c3;
    ctype const c15 = A.c4 * B.c29 - A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                      A.c8 * B.c21 + A.c12 * B.c22 + A.c13 * B.c23 + A.c14 * B.c24 +
                      A.c19 * B.c6 + A.c20 * B.c7 + A.c21 * B.c8 - A.c22 * B.c12 -
                      A.c23 * B.c13 - A.c24 * B.c14 - A.c29 * B.c4 + A.c30 * B.c5;
    ctype const c16 = A.c2 * B.c3 - A.c3 * B.c2 - A.c7 * B.c14 + A.c8 * B.c13 -
                      A.c10 * B.c11 + A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c7 +
                      A.c17 * B.c18 - A.c18 * B.c17 + A.c20 * B.c24 - A.c21 * B.c23 +
                      A.c23 * B.c21 - A.c24 * B.c20 - A.c27 * B.c28 + A.c28 * B.c27;
    ctype const c17 = -A.c1 * B.c3 + A.c3 * B.c1 + A.c6 * B.c14 - A.c8 * B.c12 +
                      A.c9 * B.c11 - A.c11 * B.c9 + A.c12 * B.c8 - A.c14 * B.c6 -
                      A.c16 * B.c18 + A.c18 * B.c16 - A.c19 * B.c24 + A.c21 * B.c22 -
                      A.c22 * B.c21 + A.c24 * B.c19 + A.c26 * B.c28 - A.c28 * B.c26;
    ctype const c18 = A.c1 * B.c2 - A.c2 * B.c1 - A.c6 * B.c13 + A.c7 * B.c12 -
                      A.c9 * B.c10 + A.c10 * B.c9 - A.c12 * B.c7 + A.c13 * B.c6 +
                      A.c16 * B.c17 - A.c17 * B.c16 + A.c19 * B.c23 - A.c20 * B.c22 +
                      A.c22 * B.c20 - A.c23 * B.c19 - A.c26 * B.c27 + A.c27 * B.c26;
    ctype const c19 = A.c1 * B.c5 - A.c5 * B.c1 - A.c10 * B.c14 + A.c11 * B.c13 -
                      A.c12 * B.c15 - A.c13 * B.c11 + A.c14 * B.c10 + A.c15 * B.c12 +
                      A.c17 * B.c21 - A.c18 * B.c20 - A.c19 * B.c25 + A.c20 * B.c18 -
                      A.c21 * B.c17 + A.c25 * B.c19 + A.c26 * B.c29 - A.c29 * B.c26;
    ctype const c20 = A.c2 * B.c5 - A.c5 * B.c2 + A.c9 * B.c14 - A.c11 * B.c12 +
                      A.c12 * B.c11 - A.c13 * B.c15 - A.c14 * B.c9 + A.c15 * B.c13 -
                      A.c16 * B.c21 + A.c18 * B.c19 - A.c19 * B.c18 - A.c20 * B.c25 +
                      A.c21 * B.c16 + A.c25 * B.c20 + A.c27 * B.c29 - A.c29 * B.c27;
    ctype const c21 = A.c3 * B.c5 - A.c5 * B.c3 - A.c9 * B.c13 + A.c10 * B.c12 -
                      A.c12 * B.c10 + A.c13 * B.c9 - A.c14 * B.c15 + A.c15 * B.c14 +
                      A.c16 * B.c20 - A.c17 * B.c19 + A.c19 * B.c17 - A.c20 * B.c16 -
                      A.c21 * B.c25 + A.c25 * B.c21 + A.c28 * B.c29 - A.c29 * B.c28;
    ctype const c22 = -A.c1 * B.c4 + A.c4 * B.c1 + A.c6 * B.c15 - A.c7 * B.c11 +
                      A.c8 * B.c10 - A.c10 * B.c8 + A.c11 * B.c7 - A.c15 * B.c6 +
                      A.c17 * B.c24 - A.c18 * B.c23 + A.c22 * B.c25 + A.c23 * B.c18 -
                      A.c24 * B.c17 - A.c25 * B.c22 - A.c26 * B.c30 + A.c30 * B.c26;
    ctype const c23 = -A.c2 * B.c4 + A.c4 * B.c2 + A.c6 * B.c11 + A.c7 * B.c15 -
                      A.c8 * B.c9 + A.c9 * B.c8 - A.c11 * B.c6 - A.c15 * B.c7 -
                      A.c16 * B.c24 + A.c18 * B.c22 - A.c22 * B.c18 + A.c23 * B.c25 +
                      A.c24 * B.c16 - A.c25 * B.c23 - A.c27 * B.c30 + A.c30 * B.c27;
    ctype const c24 = -A.c3 * B.c4 + A.c4 * B.c3 - A.c6 * B.c10 + A.c7 * B.c9 +
                      A.c8 * B.c15 - A.c9 * B.c7 + A.c10 * B.c6 - A.c15 * B.c8 +
                      A.c16 * B.c23 - A.c17 * B.c22 + A.c22 * B.c17 - A.c23 * B.c16 +
                      A.c24 * B.c25 - A.c25 * B.c24 - A.c28 * B.c30 + A.c30 * B.c28;
    ctype const c25 = -A.c4 * B.c5 + A.c5 * B.c4 - A.c6 * B.c12 - A.c7 * B.c13 -
                      A.c8 * B.c14 + A.c12 * B.c6 + A.c13 * B.c7 + A.c14 * B.c8 -
                      A.c19 * B.c22 - A.c20 * B.c23 - A.c21 * B.c24 + A.c22 * B.c19 +
                      A.c23 * B.c20 + A.c24 * B.c21 - A.c29 * B.c30 + A.c30 * B.c29;
    ctype const c26 = -A.c2 * B.c11 + A.c3 * B.c10 + A.c4 * B.c12 - A.c5 * B.c6 +
                      A.c6 * B.c5 - A.c10 * B.c3 + A.c11 * B.c2 - A.c12 * B.c4 +
                      A.c17 * B.c28 - A.c18 * B.c27 - A.c19 * B.c30 + A.c22 * B.c29 +
                      A.c27 * B.c18 - A.c28 * B.c17 - A.c29 * B.c22 + A.c30 * B.c19;
    ctype const c27 = A.c1 * B.c11 - A.c3 * B.c9 + A.c4 * B.c13 - A.c5 * B.c7 +
                      A.c7 * B.c5 + A.c9 * B.c3 - A.c11 * B.c1 - A.c13 * B.c4 -
                      A.c16 * B.c28 + A.c18 * B.c26 - A.c20 * B.c30 + A.c23 * B.c29 -
                      A.c26 * B.c18 + A.c28 * B.c16 - A.c29 * B.c23 + A.c30 * B.c20;
    ctype const c28 = -A.c1 * B.c10 + A.c2 * B.c9 + A.c4 * B.c14 - A.c5 * B.c8 +
                      A.c8 * B.c5 - A.c9 * B.c2 + A.c10 * B.c1 - A.c14 * B.c4 +
                      A.c16 * B.c27 - A.c17 * B.c26 - A.c21 * B.c30 + A.c24 * B.c29 +
                      A.c26 * B.c17 - A.c27 * B.c16 - A.c29 * B.c24 + A.c30 * B.c21;
    ctype const c29 = -A.c1 * B.c12 - A.c2 * B.c13 - A.c3 * B.c14 + A.c5 * B.c15 +
                      A.c12 * B.c1 + A.c13 * B.c2 + A.c14 * B.c3 - A.c15 * B.c5 -
                      A.c19 * B.c26 - A.c20 * B.c27 - A.c21 * B.c28 + A.c25 * B.c29 +
                      A.c26 * B.c19 + A.c27 * B.c20 + A.c28 * B.c21 - A.c29 * B.c25;
    ctype const c30 = A.c1 * B.c6 + A.c2 * B.c7 + A.c3 * B.c8 - A.c4 * B.c15 -
                      A.c6 * B.c1 - A.c7 * B.c2 - A.c8 * B.c3 + A.c15 * B.c4 +
                      A.c22 * B.c26 + A.c23 * B.c27 + A.c24 * B.c28 - A.c25 * B.c30 -
                      A.c26 * B.c22 - A.c27 * B.c23 - A.c28 * B.c24 + A.c30 * B.c25;
    ctype const c31 = 0.0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rcmt :: rcmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dc<T>,
                                                   [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(ps,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dc<T>,
                                                   [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(quadvec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] QuadVec3dc<T> const&,
                                                   [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dc<T>,
                                                   [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] TriVec3dc<T> const&,
                                                   [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dc<T>,
                                                   [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec3dc<T> const&,
                                                   [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dc<T>,
                                                   [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Vec3dc<T> const&,
                                                   [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dc<T>,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] PScalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(quadvec,quadvec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rcmt(QuadVec3dc<T> const& Q1,
                                                   QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q1.y * Q2.z + Q1.z * Q2.y;
    ctype const c1 = Q1.x * Q2.z - Q1.z * Q2.x;
    ctype const c2 = -Q1.x * Q2.y + Q1.y * Q2.x;
    ctype const c3 = Q1.x * Q2.w - Q1.w * Q2.x;
    ctype const c4 = Q1.y * Q2.w - Q1.w * Q2.y;
    ctype const c5 = Q1.z * Q2.w - Q1.w * Q2.z;
    ctype const c6 = -Q1.x * Q2.u + Q1.u * Q2.x;
    ctype const c7 = -Q1.y * Q2.u + Q1.u * Q2.y;
    ctype const c8 = -Q1.z * Q2.u + Q1.u * Q2.z;
    ctype const c9 = -Q1.w * Q2.u + Q1.u * Q2.w;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(quadvec,trivec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rcmt(QuadVec3dc<T> const& Q,
                                                    TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * t.vz - Q.z * t.vy - Q.w * t.px + Q.u * t.mx;
    ctype const c1 = -Q.x * t.vz + Q.z * t.vx - Q.w * t.py + Q.u * t.my;
    ctype const c2 = Q.x * t.vy - Q.y * t.vx - Q.w * t.pz + Q.u * t.mz;
    ctype const c3 = Q.x * t.mx + Q.y * t.my + Q.z * t.mz - Q.w * t.pw;
    ctype const c4 = -Q.x * t.px - Q.y * t.py - Q.z * t.pz + Q.u * t.pw;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(trivec,quadvec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rcmt(TriVec3dc<T> const& t,
                                                    QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vy * Q.z - t.vz * Q.y - t.mx * Q.u + t.px * Q.w;
    ctype const c1 = -t.vx * Q.z + t.vz * Q.x - t.my * Q.u + t.py * Q.w;
    ctype const c2 = t.vx * Q.y - t.vy * Q.x - t.mz * Q.u + t.pz * Q.w;
    ctype const c3 = -t.mx * Q.x - t.my * Q.y - t.mz * Q.z + t.pw * Q.w;
    ctype const c4 = t.px * Q.x + t.py * Q.y + t.pz * Q.z - t.pw * Q.u;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(quadvec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rcmt(QuadVec3dc<T> const& Q,
                                                BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * B.mz - Q.z * B.my - Q.w * B.vx + Q.u * B.px;
    ctype const c1 = -Q.x * B.mz + Q.z * B.mx - Q.w * B.vy + Q.u * B.py;
    ctype const c2 = Q.x * B.my - Q.y * B.mx - Q.w * B.vz + Q.u * B.pz;
    ctype const c3 = Q.x * B.vx + Q.y * B.vy + Q.z * B.vz + Q.u * B.pw;
    ctype const c4 = -Q.x * B.px - Q.y * B.py - Q.z * B.pz - Q.w * B.pw;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(bivec,quadvec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rcmt(BiVec3dc<T> const& B,
                                                QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * Q.w + B.my * Q.z - B.mz * Q.y - B.px * Q.u;
    ctype const c1 = B.vy * Q.w - B.mx * Q.z + B.mz * Q.x - B.py * Q.u;
    ctype const c2 = B.vz * Q.w + B.mx * Q.y - B.my * Q.x - B.pz * Q.u;
    ctype const c3 = -B.vx * Q.x - B.vy * Q.y - B.vz * Q.z - B.pw * Q.u;
    ctype const c4 = B.px * Q.x + B.py * Q.y + B.pz * Q.z + B.pw * Q.w;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(quadvec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rcmt(QuadVec3dc<T> const& Q,
                                                  Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * v.w + Q.u * v.x;
    ctype const c1 = Q.y * v.w + Q.u * v.y;
    ctype const c2 = Q.z * v.w + Q.u * v.z;
    ctype const c3 = -Q.y * v.z + Q.z * v.y;
    ctype const c4 = Q.x * v.z - Q.z * v.x;
    ctype const c5 = -Q.x * v.y + Q.y * v.x;
    ctype const c6 = -Q.x * v.u - Q.w * v.x;
    ctype const c7 = -Q.y * v.u - Q.w * v.y;
    ctype const c8 = -Q.z * v.u - Q.w * v.z;
    ctype const c9 = -Q.w * v.w + Q.u * v.u;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(vec,quadvec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rcmt(Vec3dc<T> const& v,
                                                  QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * Q.u - v.w * Q.x;
    ctype const c1 = -v.y * Q.u - v.w * Q.y;
    ctype const c2 = -v.z * Q.u - v.w * Q.z;
    ctype const c3 = -v.y * Q.z + v.z * Q.y;
    ctype const c4 = v.x * Q.z - v.z * Q.x;
    ctype const c5 = -v.x * Q.y + v.y * Q.x;
    ctype const c6 = v.x * Q.w + v.u * Q.x;
    ctype const c7 = v.y * Q.w + v.u * Q.y;
    ctype const c8 = v.z * Q.w + v.u * Q.z;
    ctype const c9 = v.w * Q.w - v.u * Q.u;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(quadvec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] QuadVec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(s,quadvec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] QuadVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(trivec,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rcmt(TriVec3dc<T> const& t1,
                                                   TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.vy * t2.vz - t1.vz * t2.vy + t1.my * t2.pz - t1.mz * t2.py +
                     t1.py * t2.mz - t1.pz * t2.my;
    ctype const c1 = -t1.vx * t2.vz + t1.vz * t2.vx - t1.mx * t2.pz + t1.mz * t2.px -
                     t1.px * t2.mz + t1.pz * t2.mx;
    ctype const c2 = t1.vx * t2.vy - t1.vy * t2.vx + t1.mx * t2.py - t1.my * t2.px +
                     t1.px * t2.my - t1.py * t2.mx;
    ctype const c3 = t1.vy * t2.mz - t1.vz * t2.my - t1.mx * t2.pw + t1.my * t2.vz -
                     t1.mz * t2.vy + t1.pw * t2.mx;
    ctype const c4 = -t1.vx * t2.mz + t1.vz * t2.mx - t1.mx * t2.vz - t1.my * t2.pw +
                     t1.mz * t2.vx + t1.pw * t2.my;
    ctype const c5 = t1.vx * t2.my - t1.vy * t2.mx + t1.mx * t2.vy - t1.my * t2.vx -
                     t1.mz * t2.pw + t1.pw * t2.mz;
    ctype const c6 = t1.vy * t2.pz - t1.vz * t2.py + t1.px * t2.pw + t1.py * t2.vz -
                     t1.pz * t2.vy - t1.pw * t2.px;
    ctype const c7 = -t1.vx * t2.pz + t1.vz * t2.px - t1.px * t2.vz + t1.py * t2.pw +
                     t1.pz * t2.vx - t1.pw * t2.py;
    ctype const c8 = t1.vx * t2.py - t1.vy * t2.px + t1.px * t2.vy - t1.py * t2.vx +
                     t1.pz * t2.pw - t1.pw * t2.pz;
    ctype const c9 = -t1.mx * t2.px - t1.my * t2.py - t1.mz * t2.pz + t1.px * t2.mx +
                     t1.py * t2.my + t1.pz * t2.mz;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(trivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rcmt(TriVec3dc<T> const& t,
                                                  BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        t.vy * B.vz - t.vz * B.vy - t.px * B.pw + t.py * B.mz - t.pz * B.my - t.pw * B.vx;
    ctype const c1 = -t.vx * B.vz + t.vz * B.vx - t.px * B.mz - t.py * B.pw +
                     t.pz * B.mx - t.pw * B.vy;
    ctype const c2 =
        t.vx * B.vy - t.vy * B.vx + t.px * B.my - t.py * B.mx - t.pz * B.pw - t.pw * B.vz;
    ctype const c3 =
        t.vy * B.mz - t.vz * B.my + t.my * B.vz - t.mz * B.vy + t.py * B.pz - t.pz * B.py;
    ctype const c4 = -t.vx * B.mz + t.vz * B.mx - t.mx * B.vz + t.mz * B.vx -
                     t.px * B.pz + t.pz * B.px;
    ctype const c5 =
        t.vx * B.my - t.vy * B.mx + t.mx * B.vy - t.my * B.vx + t.px * B.py - t.py * B.px;
    ctype const c6 =
        t.vy * B.pz - t.vz * B.py + t.mx * B.pw + t.my * B.mz - t.mz * B.my + t.pw * B.px;
    ctype const c7 = -t.vx * B.pz + t.vz * B.px - t.mx * B.mz + t.my * B.pw +
                     t.mz * B.mx + t.pw * B.py;
    ctype const c8 =
        t.vx * B.py - t.vy * B.px + t.mx * B.my - t.my * B.mx + t.mz * B.pw + t.pw * B.pz;
    ctype const c9 =
        t.mx * B.vx + t.my * B.vy + t.mz * B.vz - t.px * B.px - t.py * B.py - t.pz * B.pz;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(bivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rcmt(BiVec3dc<T> const& B,
                                                  TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * t.pw + B.vy * t.vz - B.vz * t.vy + B.my * t.pz - B.mz * t.py + B.pw * t.px;
    ctype const c1 = -B.vx * t.vz + B.vy * t.pw + B.vz * t.vx - B.mx * t.pz +
                     B.mz * t.px + B.pw * t.py;
    ctype const c2 =
        B.vx * t.vy - B.vy * t.vx + B.vz * t.pw + B.mx * t.py - B.my * t.px + B.pw * t.pz;
    ctype const c3 =
        B.vy * t.mz - B.vz * t.my + B.my * t.vz - B.mz * t.vy + B.py * t.pz - B.pz * t.py;
    ctype const c4 = -B.vx * t.mz + B.vz * t.mx - B.mx * t.vz + B.mz * t.vx -
                     B.px * t.pz + B.pz * t.px;
    ctype const c5 =
        B.vx * t.my - B.vy * t.mx + B.mx * t.vy - B.my * t.vx + B.px * t.py - B.py * t.px;
    ctype const c6 =
        B.my * t.mz - B.mz * t.my - B.px * t.pw + B.py * t.vz - B.pz * t.vy - B.pw * t.mx;
    ctype const c7 = -B.mx * t.mz + B.mz * t.mx - B.px * t.vz - B.py * t.pw +
                     B.pz * t.vx - B.pw * t.my;
    ctype const c8 =
        B.mx * t.my - B.my * t.mx + B.px * t.vy - B.py * t.vx - B.pz * t.pw - B.pw * t.mz;
    ctype const c9 = -B.vx * t.mx - B.vy * t.my - B.vz * t.mz + B.px * t.px +
                     B.py * t.py + B.pz * t.pz;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(trivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rcmt(TriVec3dc<T> const& t, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vy * v.z - t.vz * v.y + t.mx * v.w - t.px * v.u;
    ctype const c1 = -t.vx * v.z + t.vz * v.x + t.my * v.w - t.py * v.u;
    ctype const c2 = t.vx * v.y - t.vy * v.x + t.mz * v.w - t.pz * v.u;
    ctype const c3 = -t.px * v.x - t.py * v.y - t.pz * v.z - t.pw * v.w;
    ctype const c4 = t.mx * v.x + t.my * v.y + t.mz * v.z + t.pw * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(vec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rcmt(Vec3dc<T> const& v, TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * t.vz - v.z * t.vy - v.w * t.mx + v.u * t.px;
    ctype const c1 = -v.x * t.vz + v.z * t.vx - v.w * t.my + v.u * t.py;
    ctype const c2 = v.x * t.vy - v.y * t.vx - v.w * t.mz + v.u * t.pz;
    ctype const c3 = v.x * t.px + v.y * t.py + v.z * t.pz + v.w * t.pw;
    ctype const c4 = -v.x * t.mx - v.y * t.my - v.z * t.mz - v.u * t.pw;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] TriVec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] TriVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(bivec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rcmt(BiVec3dc<T> const& B1,
                                                   BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vy * B2.pz + B1.vz * B2.py - B1.my * B2.mz + B1.mz * B2.my -
                     B1.py * B2.vz + B1.pz * B2.vy;
    ctype const c1 = B1.vx * B2.pz - B1.vz * B2.px + B1.mx * B2.mz - B1.mz * B2.mx +
                     B1.px * B2.vz - B1.pz * B2.vx;
    ctype const c2 = -B1.vx * B2.py + B1.vy * B2.px - B1.mx * B2.my + B1.my * B2.mx -
                     B1.px * B2.vy + B1.py * B2.vx;
    ctype const c3 = -B1.my * B2.pz + B1.mz * B2.py - B1.px * B2.pw - B1.py * B2.mz +
                     B1.pz * B2.my + B1.pw * B2.px;
    ctype const c4 = B1.mx * B2.pz - B1.mz * B2.px + B1.px * B2.mz - B1.py * B2.pw -
                     B1.pz * B2.mx + B1.pw * B2.py;
    ctype const c5 = -B1.mx * B2.py + B1.my * B2.px - B1.px * B2.my + B1.py * B2.mx -
                     B1.pz * B2.pw + B1.pw * B2.pz;
    ctype const c6 = B1.vx * B2.pw - B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz +
                     B1.mz * B2.vy - B1.pw * B2.vx;
    ctype const c7 = B1.vx * B2.mz + B1.vy * B2.pw - B1.vz * B2.mx + B1.mx * B2.vz -
                     B1.mz * B2.vx - B1.pw * B2.vy;
    ctype const c8 = -B1.vx * B2.my + B1.vy * B2.mx + B1.vz * B2.pw - B1.mx * B2.vy +
                     B1.my * B2.vx - B1.pw * B2.vz;
    ctype const c9 = -B1.vx * B2.px - B1.vy * B2.py - B1.vz * B2.pz + B1.px * B2.vx +
                     B1.py * B2.vy + B1.pz * B2.vz;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(bivec,vec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rcmt(BiVec3dc<T> const& B,
                                                    Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.u - B.my * v.z + B.mz * v.y - B.px * v.w;
    ctype const c1 = B.vy * v.u + B.mx * v.z - B.mz * v.x - B.py * v.w;
    ctype const c2 = B.vz * v.u - B.mx * v.y + B.my * v.x - B.pz * v.w;
    ctype const c3 = B.px * v.x + B.py * v.y + B.pz * v.z - B.pw * v.u;
    ctype const c4 = -B.vx * v.x - B.vy * v.y - B.vz * v.z + B.pw * v.w;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(vec,bivec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rcmt(Vec3dc<T> const& v,
                                                    BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my + v.w * B.px - v.u * B.vx;
    ctype const c1 = v.x * B.mz - v.z * B.mx + v.w * B.py - v.u * B.vy;
    ctype const c2 = -v.x * B.my + v.y * B.mx + v.w * B.pz - v.u * B.vz;
    ctype const c3 = -v.x * B.px - v.y * B.py - v.z * B.pz + v.u * B.pw;
    ctype const c4 = v.x * B.vx + v.y * B.vy + v.z * B.vz - v.w * B.pw;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rcmt :: rcmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] BiVec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(vec,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rcmt(Vec3dc<T> const& v1,
                                                   Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.y * v2.z - v1.z * v2.y;
    ctype const c1 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c2 = v1.x * v2.y - v1.y * v2.x;
    ctype const c3 = v1.x * v2.u - v1.u * v2.x;
    ctype const c4 = v1.y * v2.u - v1.u * v2.y;
    ctype const c5 = v1.z * v2.u - v1.u * v2.z;
    ctype const c6 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c7 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c8 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c9 = -v1.w * v2.u + v1.u * v2.w;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rcmt :: rcmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Vec3dc<T> const&,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] Vec3dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

// cga3dc rcmt :: rcmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dc<T>,
                                                   [[maybe_unused]] Scalar3dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(0.0);
}

////////////////////////////////////////////////////////////////////////////////
// geometric products
//
// The null-pair metric makes single basis products multi-term:
//
//     e4 * e5 = -1 + e45        e45 * e45 = 1
//
// (a diagonal metric never mixes a scalar into a basis-blade product). The
// expressions below are generated via the diagonal-basis detour with an
// always-on direct null-basis cross-check (see ga_prdxpr).
////////////////////////////////////////////////////////////////////////////////

// cga3dc gpr :: gpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& A,
                                                      MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c5 -
                     A.c5 * B.c4 - A.c6 * B.c12 - A.c7 * B.c13 - A.c8 * B.c14 -
                     A.c9 * B.c9 - A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c6 -
                     A.c13 * B.c7 - A.c14 * B.c8 + A.c15 * B.c15 + A.c16 * B.c16 +
                     A.c17 * B.c17 + A.c18 * B.c18 + A.c19 * B.c22 + A.c20 * B.c23 +
                     A.c21 * B.c24 + A.c22 * B.c19 + A.c23 * B.c20 + A.c24 * B.c21 -
                     A.c25 * B.c25 - A.c26 * B.c26 - A.c27 * B.c27 - A.c28 * B.c28 +
                     A.c29 * B.c30 + A.c30 * B.c29 - A.c31 * B.c31;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c11 + A.c3 * B.c10 +
                     A.c4 * B.c12 - A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c21 +
                     A.c8 * B.c20 + A.c9 * B.c25 - A.c10 * B.c3 + A.c11 * B.c2 -
                     A.c12 * B.c4 + A.c13 * B.c24 - A.c14 * B.c23 - A.c15 * B.c16 -
                     A.c16 * B.c15 + A.c17 * B.c28 - A.c18 * B.c27 - A.c19 * B.c30 +
                     A.c20 * B.c8 - A.c21 * B.c7 + A.c22 * B.c29 - A.c23 * B.c14 +
                     A.c24 * B.c13 + A.c25 * B.c9 - A.c26 * B.c31 + A.c27 * B.c18 -
                     A.c28 * B.c17 - A.c29 * B.c22 + A.c30 * B.c19 - A.c31 * B.c26;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c11 + A.c2 * B.c0 - A.c3 * B.c9 +
                     A.c4 * B.c13 - A.c5 * B.c7 + A.c6 * B.c21 + A.c7 * B.c5 -
                     A.c8 * B.c19 + A.c9 * B.c3 + A.c10 * B.c25 - A.c11 * B.c1 -
                     A.c12 * B.c24 - A.c13 * B.c4 + A.c14 * B.c22 - A.c15 * B.c17 -
                     A.c16 * B.c28 - A.c17 * B.c15 + A.c18 * B.c26 - A.c19 * B.c8 -
                     A.c20 * B.c30 + A.c21 * B.c6 + A.c22 * B.c14 + A.c23 * B.c29 -
                     A.c24 * B.c12 + A.c25 * B.c10 - A.c26 * B.c18 - A.c27 * B.c31 +
                     A.c28 * B.c16 - A.c29 * B.c23 + A.c30 * B.c20 - A.c31 * B.c27;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c10 + A.c2 * B.c9 + A.c3 * B.c0 +
                     A.c4 * B.c14 - A.c5 * B.c8 - A.c6 * B.c20 + A.c7 * B.c19 +
                     A.c8 * B.c5 - A.c9 * B.c2 + A.c10 * B.c1 + A.c11 * B.c25 +
                     A.c12 * B.c23 - A.c13 * B.c22 - A.c14 * B.c4 - A.c15 * B.c18 +
                     A.c16 * B.c27 - A.c17 * B.c26 - A.c18 * B.c15 + A.c19 * B.c7 -
                     A.c20 * B.c6 - A.c21 * B.c30 - A.c22 * B.c13 + A.c23 * B.c12 +
                     A.c24 * B.c29 + A.c25 * B.c11 + A.c26 * B.c17 - A.c27 * B.c16 -
                     A.c28 * B.c31 - A.c29 * B.c24 + A.c30 * B.c21 - A.c31 * B.c28;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c6 - A.c2 * B.c7 - A.c3 * B.c8 + A.c4 * B.c0 +
                     A.c4 * B.c15 + A.c6 * B.c1 - A.c6 * B.c16 + A.c7 * B.c2 -
                     A.c7 * B.c17 + A.c8 * B.c3 - A.c8 * B.c18 - A.c9 * B.c22 -
                     A.c10 * B.c23 - A.c11 * B.c24 - A.c15 * B.c4 - A.c16 * B.c6 -
                     A.c17 * B.c7 - A.c18 * B.c8 - A.c22 * B.c9 - A.c22 * B.c26 -
                     A.c23 * B.c10 - A.c23 * B.c27 - A.c24 * B.c11 - A.c24 * B.c28 +
                     A.c25 * B.c30 + A.c26 * B.c22 + A.c27 * B.c23 + A.c28 * B.c24 -
                     A.c30 * B.c25 + A.c30 * B.c31 + A.c31 * B.c30;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                     A.c5 * B.c0 - A.c5 * B.c15 - A.c9 * B.c19 - A.c10 * B.c20 -
                     A.c11 * B.c21 - A.c12 * B.c1 - A.c12 * B.c16 - A.c13 * B.c2 -
                     A.c13 * B.c17 - A.c14 * B.c3 - A.c14 * B.c18 + A.c15 * B.c5 -
                     A.c16 * B.c12 - A.c17 * B.c13 - A.c18 * B.c14 - A.c19 * B.c9 +
                     A.c19 * B.c26 - A.c20 * B.c10 + A.c20 * B.c27 - A.c21 * B.c11 +
                     A.c21 * B.c28 - A.c25 * B.c29 - A.c26 * B.c19 - A.c27 * B.c20 -
                     A.c28 * B.c21 + A.c29 * B.c25 + A.c29 * B.c31 + A.c31 * B.c29;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c4 + A.c2 * B.c24 - A.c3 * B.c23 +
                     A.c4 * B.c1 - A.c4 * B.c16 + A.c6 * B.c0 + A.c6 * B.c15 -
                     A.c7 * B.c11 - A.c7 * B.c28 + A.c8 * B.c10 + A.c8 * B.c27 +
                     A.c9 * B.c30 - A.c10 * B.c8 + A.c11 * B.c7 - A.c15 * B.c6 -
                     A.c16 * B.c4 + A.c17 * B.c24 - A.c18 * B.c23 + A.c22 * B.c25 -
                     A.c22 * B.c31 - A.c23 * B.c3 + A.c23 * B.c18 + A.c24 * B.c2 -
                     A.c24 * B.c17 - A.c25 * B.c22 - A.c26 * B.c30 + A.c27 * B.c8 -
                     A.c28 * B.c7 + A.c30 * B.c9 + A.c30 * B.c26 - A.c31 * B.c22;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c24 - A.c2 * B.c4 + A.c3 * B.c22 +
                     A.c4 * B.c2 - A.c4 * B.c17 + A.c6 * B.c11 + A.c6 * B.c28 +
                     A.c7 * B.c0 + A.c7 * B.c15 - A.c8 * B.c9 - A.c8 * B.c26 +
                     A.c9 * B.c8 + A.c10 * B.c30 - A.c11 * B.c6 - A.c15 * B.c7 -
                     A.c16 * B.c24 - A.c17 * B.c4 + A.c18 * B.c22 + A.c22 * B.c3 -
                     A.c22 * B.c18 + A.c23 * B.c25 - A.c23 * B.c31 - A.c24 * B.c1 +
                     A.c24 * B.c16 - A.c25 * B.c23 - A.c26 * B.c8 - A.c27 * B.c30 +
                     A.c28 * B.c6 + A.c30 * B.c10 + A.c30 * B.c27 - A.c31 * B.c23;
    ctype const c8 = A.c0 * B.c8 + A.c1 * B.c23 - A.c2 * B.c22 - A.c3 * B.c4 +
                     A.c4 * B.c3 - A.c4 * B.c18 - A.c6 * B.c10 - A.c6 * B.c27 +
                     A.c7 * B.c9 + A.c7 * B.c26 + A.c8 * B.c0 + A.c8 * B.c15 -
                     A.c9 * B.c7 + A.c10 * B.c6 + A.c11 * B.c30 - A.c15 * B.c8 +
                     A.c16 * B.c23 - A.c17 * B.c22 - A.c18 * B.c4 - A.c22 * B.c2 +
                     A.c22 * B.c17 + A.c23 * B.c1 - A.c23 * B.c16 + A.c24 * B.c25 -
                     A.c24 * B.c31 - A.c25 * B.c24 + A.c26 * B.c7 - A.c27 * B.c6 -
                     A.c28 * B.c30 + A.c30 * B.c11 + A.c30 * B.c28 - A.c31 * B.c24;
    ctype const c9 = A.c0 * B.c9 - A.c1 * B.c25 + A.c2 * B.c3 - A.c3 * B.c2 -
                     A.c4 * B.c19 - A.c5 * B.c22 + A.c6 * B.c29 - A.c7 * B.c14 +
                     A.c8 * B.c13 + A.c9 * B.c0 - A.c10 * B.c11 + A.c11 * B.c10 +
                     A.c12 * B.c30 - A.c13 * B.c8 + A.c14 * B.c7 + A.c15 * B.c26 -
                     A.c16 * B.c31 + A.c17 * B.c18 - A.c18 * B.c17 - A.c19 * B.c4 +
                     A.c20 * B.c24 - A.c21 * B.c23 - A.c22 * B.c5 + A.c23 * B.c21 -
                     A.c24 * B.c20 - A.c25 * B.c1 + A.c26 * B.c15 - A.c27 * B.c28 +
                     A.c28 * B.c27 + A.c29 * B.c6 + A.c30 * B.c12 - A.c31 * B.c16;
    ctype const c10 = A.c0 * B.c10 - A.c1 * B.c3 - A.c2 * B.c25 + A.c3 * B.c1 -
                      A.c4 * B.c20 - A.c5 * B.c23 + A.c6 * B.c14 + A.c7 * B.c29 -
                      A.c8 * B.c12 + A.c9 * B.c11 + A.c10 * B.c0 - A.c11 * B.c9 +
                      A.c12 * B.c8 + A.c13 * B.c30 - A.c14 * B.c6 + A.c15 * B.c27 -
                      A.c16 * B.c18 - A.c17 * B.c31 + A.c18 * B.c16 - A.c19 * B.c24 -
                      A.c20 * B.c4 + A.c21 * B.c22 - A.c22 * B.c21 - A.c23 * B.c5 +
                      A.c24 * B.c19 - A.c25 * B.c2 + A.c26 * B.c28 + A.c27 * B.c15 -
                      A.c28 * B.c26 + A.c29 * B.c7 + A.c30 * B.c13 - A.c31 * B.c17;
    ctype const c11 = A.c0 * B.c11 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c25 -
                      A.c4 * B.c21 - A.c5 * B.c24 - A.c6 * B.c13 + A.c7 * B.c12 +
                      A.c8 * B.c29 - A.c9 * B.c10 + A.c10 * B.c9 + A.c11 * B.c0 -
                      A.c12 * B.c7 + A.c13 * B.c6 + A.c14 * B.c30 + A.c15 * B.c28 +
                      A.c16 * B.c17 - A.c17 * B.c16 - A.c18 * B.c31 + A.c19 * B.c23 -
                      A.c20 * B.c22 - A.c21 * B.c4 + A.c22 * B.c20 - A.c23 * B.c19 -
                      A.c24 * B.c5 - A.c25 * B.c3 - A.c26 * B.c27 + A.c27 * B.c26 +
                      A.c28 * B.c15 + A.c29 * B.c8 + A.c30 * B.c14 - A.c31 * B.c18;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c5 - A.c2 * B.c21 + A.c3 * B.c20 -
                      A.c5 * B.c1 - A.c5 * B.c16 + A.c9 * B.c29 - A.c10 * B.c14 +
                      A.c11 * B.c13 + A.c12 * B.c0 - A.c12 * B.c15 - A.c13 * B.c11 +
                      A.c13 * B.c28 + A.c14 * B.c10 - A.c14 * B.c27 + A.c15 * B.c12 -
                      A.c16 * B.c5 + A.c17 * B.c21 - A.c18 * B.c20 - A.c19 * B.c25 -
                      A.c19 * B.c31 + A.c20 * B.c3 + A.c20 * B.c18 - A.c21 * B.c2 -
                      A.c21 * B.c17 + A.c25 * B.c19 + A.c26 * B.c29 - A.c27 * B.c14 +
                      A.c28 * B.c13 + A.c29 * B.c9 - A.c29 * B.c26 - A.c31 * B.c19;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c21 + A.c2 * B.c5 - A.c3 * B.c19 -
                      A.c5 * B.c2 - A.c5 * B.c17 + A.c9 * B.c14 + A.c10 * B.c29 -
                      A.c11 * B.c12 + A.c12 * B.c11 - A.c12 * B.c28 + A.c13 * B.c0 -
                      A.c13 * B.c15 - A.c14 * B.c9 + A.c14 * B.c26 + A.c15 * B.c13 -
                      A.c16 * B.c21 - A.c17 * B.c5 + A.c18 * B.c19 - A.c19 * B.c3 -
                      A.c19 * B.c18 - A.c20 * B.c25 - A.c20 * B.c31 + A.c21 * B.c1 +
                      A.c21 * B.c16 + A.c25 * B.c20 + A.c26 * B.c14 + A.c27 * B.c29 -
                      A.c28 * B.c12 + A.c29 * B.c10 - A.c29 * B.c27 - A.c31 * B.c20;
    ctype const c14 = A.c0 * B.c14 - A.c1 * B.c20 + A.c2 * B.c19 + A.c3 * B.c5 -
                      A.c5 * B.c3 - A.c5 * B.c18 - A.c9 * B.c13 + A.c10 * B.c12 +
                      A.c11 * B.c29 - A.c12 * B.c10 + A.c12 * B.c27 + A.c13 * B.c9 -
                      A.c13 * B.c26 + A.c14 * B.c0 - A.c14 * B.c15 + A.c15 * B.c14 +
                      A.c16 * B.c20 - A.c17 * B.c19 - A.c18 * B.c5 + A.c19 * B.c2 +
                      A.c19 * B.c17 - A.c20 * B.c1 - A.c20 * B.c16 - A.c21 * B.c25 -
                      A.c21 * B.c31 + A.c25 * B.c21 - A.c26 * B.c13 + A.c27 * B.c12 +
                      A.c28 * B.c29 + A.c29 * B.c11 - A.c29 * B.c28 - A.c31 * B.c21;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c16 - A.c2 * B.c17 - A.c3 * B.c18 +
                      A.c4 * B.c5 - A.c5 * B.c4 + A.c6 * B.c12 + A.c7 * B.c13 +
                      A.c8 * B.c14 - A.c9 * B.c26 - A.c10 * B.c27 - A.c11 * B.c28 -
                      A.c12 * B.c6 - A.c13 * B.c7 - A.c14 * B.c8 + A.c15 * B.c0 -
                      A.c16 * B.c1 - A.c17 * B.c2 - A.c18 * B.c3 + A.c19 * B.c22 +
                      A.c20 * B.c23 + A.c21 * B.c24 - A.c22 * B.c19 - A.c23 * B.c20 -
                      A.c24 * B.c21 + A.c25 * B.c31 - A.c26 * B.c9 - A.c27 * B.c10 -
                      A.c28 * B.c11 + A.c29 * B.c30 - A.c30 * B.c29 + A.c31 * B.c25;
    ctype const c16 = A.c0 * B.c16 - A.c1 * B.c15 + A.c2 * B.c28 - A.c3 * B.c27 +
                      A.c4 * B.c12 + A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c21 +
                      A.c8 * B.c20 + A.c9 * B.c31 - A.c10 * B.c18 + A.c11 * B.c17 +
                      A.c12 * B.c4 - A.c13 * B.c24 + A.c14 * B.c23 - A.c15 * B.c1 +
                      A.c16 * B.c0 - A.c17 * B.c11 + A.c18 * B.c10 + A.c19 * B.c30 -
                      A.c20 * B.c8 + A.c21 * B.c7 + A.c22 * B.c29 - A.c23 * B.c14 +
                      A.c24 * B.c13 - A.c25 * B.c26 - A.c26 * B.c25 + A.c27 * B.c3 -
                      A.c28 * B.c2 + A.c29 * B.c22 + A.c30 * B.c19 + A.c31 * B.c9;
    ctype const c17 = A.c0 * B.c17 - A.c1 * B.c28 - A.c2 * B.c15 + A.c3 * B.c26 +
                      A.c4 * B.c13 + A.c5 * B.c7 + A.c6 * B.c21 + A.c7 * B.c5 -
                      A.c8 * B.c19 + A.c9 * B.c18 + A.c10 * B.c31 - A.c11 * B.c16 +
                      A.c12 * B.c24 + A.c13 * B.c4 - A.c14 * B.c22 - A.c15 * B.c2 +
                      A.c16 * B.c11 + A.c17 * B.c0 - A.c18 * B.c9 + A.c19 * B.c8 +
                      A.c20 * B.c30 - A.c21 * B.c6 + A.c22 * B.c14 + A.c23 * B.c29 -
                      A.c24 * B.c12 - A.c25 * B.c27 - A.c26 * B.c3 - A.c27 * B.c25 +
                      A.c28 * B.c1 + A.c29 * B.c23 + A.c30 * B.c20 + A.c31 * B.c10;
    ctype const c18 = A.c0 * B.c18 + A.c1 * B.c27 - A.c2 * B.c26 - A.c3 * B.c15 +
                      A.c4 * B.c14 + A.c5 * B.c8 - A.c6 * B.c20 + A.c7 * B.c19 +
                      A.c8 * B.c5 - A.c9 * B.c17 + A.c10 * B.c16 + A.c11 * B.c31 -
                      A.c12 * B.c23 + A.c13 * B.c22 + A.c14 * B.c4 - A.c15 * B.c3 -
                      A.c16 * B.c10 + A.c17 * B.c9 + A.c18 * B.c0 - A.c19 * B.c7 +
                      A.c20 * B.c6 + A.c21 * B.c30 - A.c22 * B.c13 + A.c23 * B.c12 +
                      A.c24 * B.c29 - A.c25 * B.c28 + A.c26 * B.c2 - A.c27 * B.c1 -
                      A.c28 * B.c25 + A.c29 * B.c24 + A.c30 * B.c21 + A.c31 * B.c11;
    ctype const c19 = A.c0 * B.c19 - A.c1 * B.c29 + A.c2 * B.c14 - A.c3 * B.c13 +
                      A.c5 * B.c9 - A.c5 * B.c26 + A.c9 * B.c5 - A.c10 * B.c21 +
                      A.c11 * B.c20 + A.c12 * B.c25 + A.c12 * B.c31 - A.c13 * B.c3 -
                      A.c13 * B.c18 + A.c14 * B.c2 + A.c14 * B.c17 + A.c15 * B.c19 +
                      A.c16 * B.c29 - A.c17 * B.c14 + A.c18 * B.c13 + A.c19 * B.c0 -
                      A.c19 * B.c15 - A.c20 * B.c11 + A.c20 * B.c28 + A.c21 * B.c10 -
                      A.c21 * B.c27 - A.c25 * B.c12 + A.c26 * B.c5 - A.c27 * B.c21 +
                      A.c28 * B.c20 + A.c29 * B.c1 + A.c29 * B.c16 + A.c31 * B.c12;
    ctype const c20 = A.c0 * B.c20 - A.c1 * B.c14 - A.c2 * B.c29 + A.c3 * B.c12 +
                      A.c5 * B.c10 - A.c5 * B.c27 + A.c9 * B.c21 + A.c10 * B.c5 -
                      A.c11 * B.c19 + A.c12 * B.c3 + A.c12 * B.c18 + A.c13 * B.c25 +
                      A.c13 * B.c31 - A.c14 * B.c1 - A.c14 * B.c16 + A.c15 * B.c20 +
                      A.c16 * B.c14 + A.c17 * B.c29 - A.c18 * B.c12 + A.c19 * B.c11 -
                      A.c19 * B.c28 + A.c20 * B.c0 - A.c20 * B.c15 - A.c21 * B.c9 +
                      A.c21 * B.c26 - A.c25 * B.c13 + A.c26 * B.c21 + A.c27 * B.c5 -
                      A.c28 * B.c19 + A.c29 * B.c2 + A.c29 * B.c17 + A.c31 * B.c13;
    ctype const c21 = A.c0 * B.c21 + A.c1 * B.c13 - A.c2 * B.c12 - A.c3 * B.c29 +
                      A.c5 * B.c11 - A.c5 * B.c28 - A.c9 * B.c20 + A.c10 * B.c19 +
                      A.c11 * B.c5 - A.c12 * B.c2 - A.c12 * B.c17 + A.c13 * B.c1 +
                      A.c13 * B.c16 + A.c14 * B.c25 + A.c14 * B.c31 + A.c15 * B.c21 -
                      A.c16 * B.c13 + A.c17 * B.c12 + A.c18 * B.c29 - A.c19 * B.c10 +
                      A.c19 * B.c27 + A.c20 * B.c9 - A.c20 * B.c26 + A.c21 * B.c0 -
                      A.c21 * B.c15 - A.c25 * B.c14 - A.c26 * B.c20 + A.c27 * B.c19 +
                      A.c28 * B.c5 + A.c29 * B.c3 + A.c29 * B.c18 + A.c31 * B.c14;
    ctype const c22 = A.c0 * B.c22 + A.c1 * B.c30 - A.c2 * B.c8 + A.c3 * B.c7 +
                      A.c4 * B.c9 + A.c4 * B.c26 - A.c6 * B.c25 + A.c6 * B.c31 +
                      A.c7 * B.c3 - A.c7 * B.c18 - A.c8 * B.c2 + A.c8 * B.c17 +
                      A.c9 * B.c4 - A.c10 * B.c24 + A.c11 * B.c23 - A.c15 * B.c22 +
                      A.c16 * B.c30 - A.c17 * B.c8 + A.c18 * B.c7 + A.c22 * B.c0 +
                      A.c22 * B.c15 - A.c23 * B.c11 - A.c23 * B.c28 + A.c24 * B.c10 +
                      A.c24 * B.c27 + A.c25 * B.c6 - A.c26 * B.c4 + A.c27 * B.c24 -
                      A.c28 * B.c23 - A.c30 * B.c1 + A.c30 * B.c16 + A.c31 * B.c6;
    ctype const c23 = A.c0 * B.c23 + A.c1 * B.c8 + A.c2 * B.c30 - A.c3 * B.c6 +
                      A.c4 * B.c10 + A.c4 * B.c27 - A.c6 * B.c3 + A.c6 * B.c18 -
                      A.c7 * B.c25 + A.c7 * B.c31 + A.c8 * B.c1 - A.c8 * B.c16 +
                      A.c9 * B.c24 + A.c10 * B.c4 - A.c11 * B.c22 - A.c15 * B.c23 +
                      A.c16 * B.c8 + A.c17 * B.c30 - A.c18 * B.c6 + A.c22 * B.c11 +
                      A.c22 * B.c28 + A.c23 * B.c0 + A.c23 * B.c15 - A.c24 * B.c9 -
                      A.c24 * B.c26 + A.c25 * B.c7 - A.c26 * B.c24 - A.c27 * B.c4 +
                      A.c28 * B.c22 - A.c30 * B.c2 + A.c30 * B.c17 + A.c31 * B.c7;
    ctype const c24 = A.c0 * B.c24 - A.c1 * B.c7 + A.c2 * B.c6 + A.c3 * B.c30 +
                      A.c4 * B.c11 + A.c4 * B.c28 + A.c6 * B.c2 - A.c6 * B.c17 -
                      A.c7 * B.c1 + A.c7 * B.c16 - A.c8 * B.c25 + A.c8 * B.c31 -
                      A.c9 * B.c23 + A.c10 * B.c22 + A.c11 * B.c4 - A.c15 * B.c24 -
                      A.c16 * B.c7 + A.c17 * B.c6 + A.c18 * B.c30 - A.c22 * B.c10 -
                      A.c22 * B.c27 + A.c23 * B.c9 + A.c23 * B.c26 + A.c24 * B.c0 +
                      A.c24 * B.c15 + A.c25 * B.c8 + A.c26 * B.c23 - A.c27 * B.c22 -
                      A.c28 * B.c4 - A.c30 * B.c3 + A.c30 * B.c18 + A.c31 * B.c8;
    ctype const c25 = A.c0 * B.c25 - A.c1 * B.c9 - A.c2 * B.c10 - A.c3 * B.c11 +
                      A.c4 * B.c29 - A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                      A.c8 * B.c21 - A.c9 * B.c1 - A.c10 * B.c2 - A.c11 * B.c3 +
                      A.c12 * B.c22 + A.c13 * B.c23 + A.c14 * B.c24 - A.c15 * B.c31 +
                      A.c16 * B.c26 + A.c17 * B.c27 + A.c18 * B.c28 + A.c19 * B.c6 +
                      A.c20 * B.c7 + A.c21 * B.c8 - A.c22 * B.c12 - A.c23 * B.c13 -
                      A.c24 * B.c14 + A.c25 * B.c0 + A.c26 * B.c16 + A.c27 * B.c17 +
                      A.c28 * B.c18 - A.c29 * B.c4 + A.c30 * B.c5 - A.c31 * B.c15;
    ctype const c26 = A.c0 * B.c26 + A.c1 * B.c31 - A.c2 * B.c18 + A.c3 * B.c17 +
                      A.c4 * B.c19 - A.c5 * B.c22 - A.c6 * B.c29 + A.c7 * B.c14 -
                      A.c8 * B.c13 + A.c9 * B.c15 - A.c10 * B.c28 + A.c11 * B.c27 +
                      A.c12 * B.c30 - A.c13 * B.c8 + A.c14 * B.c7 + A.c15 * B.c9 +
                      A.c16 * B.c25 - A.c17 * B.c3 + A.c18 * B.c2 - A.c19 * B.c4 +
                      A.c20 * B.c24 - A.c21 * B.c23 + A.c22 * B.c5 - A.c23 * B.c21 +
                      A.c24 * B.c20 + A.c25 * B.c16 + A.c26 * B.c0 - A.c27 * B.c11 +
                      A.c28 * B.c10 + A.c29 * B.c6 - A.c30 * B.c12 + A.c31 * B.c1;
    ctype const c27 = A.c0 * B.c27 + A.c1 * B.c18 + A.c2 * B.c31 - A.c3 * B.c16 +
                      A.c4 * B.c20 - A.c5 * B.c23 - A.c6 * B.c14 - A.c7 * B.c29 +
                      A.c8 * B.c12 + A.c9 * B.c28 + A.c10 * B.c15 - A.c11 * B.c26 +
                      A.c12 * B.c8 + A.c13 * B.c30 - A.c14 * B.c6 + A.c15 * B.c10 +
                      A.c16 * B.c3 + A.c17 * B.c25 - A.c18 * B.c1 - A.c19 * B.c24 -
                      A.c20 * B.c4 + A.c21 * B.c22 + A.c22 * B.c21 + A.c23 * B.c5 -
                      A.c24 * B.c19 + A.c25 * B.c17 + A.c26 * B.c11 + A.c27 * B.c0 -
                      A.c28 * B.c9 + A.c29 * B.c7 - A.c30 * B.c13 + A.c31 * B.c2;
    ctype const c28 = A.c0 * B.c28 - A.c1 * B.c17 + A.c2 * B.c16 + A.c3 * B.c31 +
                      A.c4 * B.c21 - A.c5 * B.c24 + A.c6 * B.c13 - A.c7 * B.c12 -
                      A.c8 * B.c29 - A.c9 * B.c27 + A.c10 * B.c26 + A.c11 * B.c15 -
                      A.c12 * B.c7 + A.c13 * B.c6 + A.c14 * B.c30 + A.c15 * B.c11 -
                      A.c16 * B.c2 + A.c17 * B.c1 + A.c18 * B.c25 + A.c19 * B.c23 -
                      A.c20 * B.c22 - A.c21 * B.c4 - A.c22 * B.c20 + A.c23 * B.c19 +
                      A.c24 * B.c5 + A.c25 * B.c18 - A.c26 * B.c10 + A.c27 * B.c9 +
                      A.c28 * B.c0 + A.c29 * B.c8 - A.c30 * B.c14 + A.c31 * B.c3;
    ctype const c29 = A.c0 * B.c29 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                      A.c5 * B.c25 - A.c5 * B.c31 - A.c9 * B.c12 - A.c10 * B.c13 -
                      A.c11 * B.c14 - A.c12 * B.c9 + A.c12 * B.c26 - A.c13 * B.c10 +
                      A.c13 * B.c27 - A.c14 * B.c11 + A.c14 * B.c28 + A.c15 * B.c29 +
                      A.c16 * B.c19 + A.c17 * B.c20 + A.c18 * B.c21 + A.c19 * B.c1 +
                      A.c19 * B.c16 + A.c20 * B.c2 + A.c20 * B.c17 + A.c21 * B.c3 +
                      A.c21 * B.c18 + A.c25 * B.c5 - A.c26 * B.c12 - A.c27 * B.c13 -
                      A.c28 * B.c14 + A.c29 * B.c0 - A.c29 * B.c15 - A.c31 * B.c5;
    ctype const c30 = A.c0 * B.c30 + A.c1 * B.c22 + A.c2 * B.c23 + A.c3 * B.c24 +
                      A.c4 * B.c25 - A.c4 * B.c31 - A.c6 * B.c9 - A.c6 * B.c26 -
                      A.c7 * B.c10 - A.c7 * B.c27 - A.c8 * B.c11 - A.c8 * B.c28 -
                      A.c9 * B.c6 - A.c10 * B.c7 - A.c11 * B.c8 - A.c15 * B.c30 +
                      A.c16 * B.c22 + A.c17 * B.c23 + A.c18 * B.c24 - A.c22 * B.c1 +
                      A.c22 * B.c16 - A.c23 * B.c2 + A.c23 * B.c17 - A.c24 * B.c3 +
                      A.c24 * B.c18 - A.c25 * B.c4 + A.c26 * B.c6 + A.c27 * B.c7 +
                      A.c28 * B.c8 + A.c30 * B.c0 + A.c30 * B.c15 - A.c31 * B.c4;
    ctype const c31 = A.c0 * B.c31 + A.c1 * B.c26 + A.c2 * B.c27 + A.c3 * B.c28 +
                      A.c4 * B.c29 + A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                      A.c8 * B.c21 - A.c9 * B.c16 - A.c10 * B.c17 - A.c11 * B.c18 -
                      A.c12 * B.c22 - A.c13 * B.c23 - A.c14 * B.c24 - A.c15 * B.c25 -
                      A.c16 * B.c9 - A.c17 * B.c10 - A.c18 * B.c11 - A.c19 * B.c6 -
                      A.c20 * B.c7 - A.c21 * B.c8 - A.c22 * B.c12 - A.c23 * B.c13 -
                      A.c24 * B.c14 - A.c25 * B.c15 + A.c26 * B.c1 + A.c27 * B.c2 +
                      A.c28 * B.c3 + A.c29 * B.c4 + A.c30 * B.c5 + A.c31 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& A,
                                                      MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c6 * B.c7 - A.c7 * B.c8 - A.c8 * B.c9 - A.c9 * B.c4 -
                     A.c10 * B.c5 - A.c11 * B.c6 - A.c12 * B.c1 - A.c13 * B.c2 -
                     A.c14 * B.c3 + A.c15 * B.c10 - A.c26 * B.c11 - A.c27 * B.c12 -
                     A.c28 * B.c13 + A.c29 * B.c15 + A.c30 * B.c14;
    ctype const c1 = A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c7 - A.c5 * B.c1 -
                     A.c16 * B.c10 + A.c17 * B.c13 - A.c18 * B.c12 - A.c19 * B.c15 +
                     A.c20 * B.c3 - A.c21 * B.c2 + A.c22 * B.c14 - A.c23 * B.c9 +
                     A.c24 * B.c8 + A.c25 * B.c4 - A.c31 * B.c11;
    ctype const c2 = A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c4 * B.c8 - A.c5 * B.c2 -
                     A.c16 * B.c13 - A.c17 * B.c10 + A.c18 * B.c11 - A.c19 * B.c3 -
                     A.c20 * B.c15 + A.c21 * B.c1 + A.c22 * B.c9 + A.c23 * B.c14 -
                     A.c24 * B.c7 + A.c25 * B.c5 - A.c31 * B.c12;
    ctype const c3 = -A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 + A.c4 * B.c9 -
                     A.c5 * B.c3 + A.c16 * B.c12 - A.c17 * B.c11 - A.c18 * B.c10 +
                     A.c19 * B.c2 - A.c20 * B.c1 - A.c21 * B.c15 - A.c22 * B.c8 +
                     A.c23 * B.c7 + A.c24 * B.c14 + A.c25 * B.c6 - A.c31 * B.c13;
    ctype const c4 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c4 * B.c0 +
                     A.c4 * B.c10 - A.c16 * B.c1 - A.c17 * B.c2 - A.c18 * B.c3 -
                     A.c22 * B.c4 - A.c22 * B.c11 - A.c23 * B.c5 - A.c23 * B.c12 -
                     A.c24 * B.c6 - A.c24 * B.c13 + A.c25 * B.c15 + A.c31 * B.c15;
    ctype const c5 = A.c1 * B.c7 + A.c2 * B.c8 + A.c3 * B.c9 + A.c5 * B.c0 -
                     A.c5 * B.c10 - A.c16 * B.c7 - A.c17 * B.c8 - A.c18 * B.c9 -
                     A.c19 * B.c4 + A.c19 * B.c11 - A.c20 * B.c5 + A.c20 * B.c12 -
                     A.c21 * B.c6 + A.c21 * B.c13 - A.c25 * B.c14 + A.c31 * B.c14;
    ctype const c6 = A.c0 * B.c1 + A.c6 * B.c0 + A.c6 * B.c10 - A.c7 * B.c6 -
                     A.c7 * B.c13 + A.c8 * B.c5 + A.c8 * B.c12 + A.c9 * B.c15 -
                     A.c10 * B.c3 + A.c11 * B.c2 - A.c15 * B.c1 - A.c26 * B.c15 +
                     A.c27 * B.c3 - A.c28 * B.c2 + A.c30 * B.c4 + A.c30 * B.c11;
    ctype const c7 = A.c0 * B.c2 + A.c6 * B.c6 + A.c6 * B.c13 + A.c7 * B.c0 +
                     A.c7 * B.c10 - A.c8 * B.c4 - A.c8 * B.c11 + A.c9 * B.c3 +
                     A.c10 * B.c15 - A.c11 * B.c1 - A.c15 * B.c2 - A.c26 * B.c3 -
                     A.c27 * B.c15 + A.c28 * B.c1 + A.c30 * B.c5 + A.c30 * B.c12;
    ctype const c8 = A.c0 * B.c3 - A.c6 * B.c5 - A.c6 * B.c12 + A.c7 * B.c4 +
                     A.c7 * B.c11 + A.c8 * B.c0 + A.c8 * B.c10 - A.c9 * B.c2 +
                     A.c10 * B.c1 + A.c11 * B.c15 - A.c15 * B.c3 + A.c26 * B.c2 -
                     A.c27 * B.c1 - A.c28 * B.c15 + A.c30 * B.c6 + A.c30 * B.c13;
    ctype const c9 = A.c0 * B.c4 + A.c6 * B.c14 - A.c7 * B.c9 + A.c8 * B.c8 +
                     A.c9 * B.c0 - A.c10 * B.c6 + A.c11 * B.c5 + A.c12 * B.c15 -
                     A.c13 * B.c3 + A.c14 * B.c2 + A.c15 * B.c11 + A.c26 * B.c10 -
                     A.c27 * B.c13 + A.c28 * B.c12 + A.c29 * B.c1 + A.c30 * B.c7;
    ctype const c10 = A.c0 * B.c5 + A.c6 * B.c9 + A.c7 * B.c14 - A.c8 * B.c7 +
                      A.c9 * B.c6 + A.c10 * B.c0 - A.c11 * B.c4 + A.c12 * B.c3 +
                      A.c13 * B.c15 - A.c14 * B.c1 + A.c15 * B.c12 + A.c26 * B.c13 +
                      A.c27 * B.c10 - A.c28 * B.c11 + A.c29 * B.c2 + A.c30 * B.c8;
    ctype const c11 = A.c0 * B.c6 - A.c6 * B.c8 + A.c7 * B.c7 + A.c8 * B.c14 -
                      A.c9 * B.c5 + A.c10 * B.c4 + A.c11 * B.c0 - A.c12 * B.c2 +
                      A.c13 * B.c1 + A.c14 * B.c15 + A.c15 * B.c13 - A.c26 * B.c12 +
                      A.c27 * B.c11 + A.c28 * B.c10 + A.c29 * B.c3 + A.c30 * B.c9;
    ctype const c12 = A.c0 * B.c7 + A.c9 * B.c14 - A.c10 * B.c9 + A.c11 * B.c8 +
                      A.c12 * B.c0 - A.c12 * B.c10 - A.c13 * B.c6 + A.c13 * B.c13 +
                      A.c14 * B.c5 - A.c14 * B.c12 + A.c15 * B.c7 + A.c26 * B.c14 -
                      A.c27 * B.c9 + A.c28 * B.c8 + A.c29 * B.c4 - A.c29 * B.c11;
    ctype const c13 = A.c0 * B.c8 + A.c9 * B.c9 + A.c10 * B.c14 - A.c11 * B.c7 +
                      A.c12 * B.c6 - A.c12 * B.c13 + A.c13 * B.c0 - A.c13 * B.c10 -
                      A.c14 * B.c4 + A.c14 * B.c11 + A.c15 * B.c8 + A.c26 * B.c9 +
                      A.c27 * B.c14 - A.c28 * B.c7 + A.c29 * B.c5 - A.c29 * B.c12;
    ctype const c14 = A.c0 * B.c9 - A.c9 * B.c8 + A.c10 * B.c7 + A.c11 * B.c14 -
                      A.c12 * B.c5 + A.c12 * B.c12 + A.c13 * B.c4 - A.c13 * B.c11 +
                      A.c14 * B.c0 - A.c14 * B.c10 + A.c15 * B.c9 - A.c26 * B.c8 +
                      A.c27 * B.c7 + A.c28 * B.c14 + A.c29 * B.c6 - A.c29 * B.c13;
    ctype const c15 = A.c0 * B.c10 + A.c6 * B.c7 + A.c7 * B.c8 + A.c8 * B.c9 -
                      A.c9 * B.c11 - A.c10 * B.c12 - A.c11 * B.c13 - A.c12 * B.c1 -
                      A.c13 * B.c2 - A.c14 * B.c3 + A.c15 * B.c0 - A.c26 * B.c4 -
                      A.c27 * B.c5 - A.c28 * B.c6 + A.c29 * B.c15 - A.c30 * B.c14;
    ctype const c16 = -A.c1 * B.c10 + A.c2 * B.c13 - A.c3 * B.c12 + A.c4 * B.c7 +
                      A.c5 * B.c1 + A.c16 * B.c0 - A.c17 * B.c6 + A.c18 * B.c5 +
                      A.c19 * B.c15 - A.c20 * B.c3 + A.c21 * B.c2 + A.c22 * B.c14 -
                      A.c23 * B.c9 + A.c24 * B.c8 - A.c25 * B.c11 + A.c31 * B.c4;
    ctype const c17 = -A.c1 * B.c13 - A.c2 * B.c10 + A.c3 * B.c11 + A.c4 * B.c8 +
                      A.c5 * B.c2 + A.c16 * B.c6 + A.c17 * B.c0 - A.c18 * B.c4 +
                      A.c19 * B.c3 + A.c20 * B.c15 - A.c21 * B.c1 + A.c22 * B.c9 +
                      A.c23 * B.c14 - A.c24 * B.c7 - A.c25 * B.c12 + A.c31 * B.c5;
    ctype const c18 = A.c1 * B.c12 - A.c2 * B.c11 - A.c3 * B.c10 + A.c4 * B.c9 +
                      A.c5 * B.c3 - A.c16 * B.c5 + A.c17 * B.c4 + A.c18 * B.c0 -
                      A.c19 * B.c2 + A.c20 * B.c1 + A.c21 * B.c15 - A.c22 * B.c8 +
                      A.c23 * B.c7 + A.c24 * B.c14 - A.c25 * B.c13 + A.c31 * B.c6;
    ctype const c19 = -A.c1 * B.c14 + A.c2 * B.c9 - A.c3 * B.c8 + A.c5 * B.c4 -
                      A.c5 * B.c11 + A.c16 * B.c14 - A.c17 * B.c9 + A.c18 * B.c8 +
                      A.c19 * B.c0 - A.c19 * B.c10 - A.c20 * B.c6 + A.c20 * B.c13 +
                      A.c21 * B.c5 - A.c21 * B.c12 - A.c25 * B.c7 + A.c31 * B.c7;
    ctype const c20 = -A.c1 * B.c9 - A.c2 * B.c14 + A.c3 * B.c7 + A.c5 * B.c5 -
                      A.c5 * B.c12 + A.c16 * B.c9 + A.c17 * B.c14 - A.c18 * B.c7 +
                      A.c19 * B.c6 - A.c19 * B.c13 + A.c20 * B.c0 - A.c20 * B.c10 -
                      A.c21 * B.c4 + A.c21 * B.c11 - A.c25 * B.c8 + A.c31 * B.c8;
    ctype const c21 = A.c1 * B.c8 - A.c2 * B.c7 - A.c3 * B.c14 + A.c5 * B.c6 -
                      A.c5 * B.c13 - A.c16 * B.c8 + A.c17 * B.c7 + A.c18 * B.c14 -
                      A.c19 * B.c5 + A.c19 * B.c12 + A.c20 * B.c4 - A.c20 * B.c11 +
                      A.c21 * B.c0 - A.c21 * B.c10 - A.c25 * B.c9 + A.c31 * B.c9;
    ctype const c22 = A.c1 * B.c15 - A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c4 +
                      A.c4 * B.c11 + A.c16 * B.c15 - A.c17 * B.c3 + A.c18 * B.c2 +
                      A.c22 * B.c0 + A.c22 * B.c10 - A.c23 * B.c6 - A.c23 * B.c13 +
                      A.c24 * B.c5 + A.c24 * B.c12 + A.c25 * B.c1 + A.c31 * B.c1;
    ctype const c23 = A.c1 * B.c3 + A.c2 * B.c15 - A.c3 * B.c1 + A.c4 * B.c5 +
                      A.c4 * B.c12 + A.c16 * B.c3 + A.c17 * B.c15 - A.c18 * B.c1 +
                      A.c22 * B.c6 + A.c22 * B.c13 + A.c23 * B.c0 + A.c23 * B.c10 -
                      A.c24 * B.c4 - A.c24 * B.c11 + A.c25 * B.c2 + A.c31 * B.c2;
    ctype const c24 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c15 + A.c4 * B.c6 +
                      A.c4 * B.c13 - A.c16 * B.c2 + A.c17 * B.c1 + A.c18 * B.c15 -
                      A.c22 * B.c5 - A.c22 * B.c12 + A.c23 * B.c4 + A.c23 * B.c11 +
                      A.c24 * B.c0 + A.c24 * B.c10 + A.c25 * B.c3 + A.c31 * B.c3;
    ctype const c25 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 + A.c4 * B.c14 -
                      A.c5 * B.c15 + A.c16 * B.c11 + A.c17 * B.c12 + A.c18 * B.c13 +
                      A.c19 * B.c1 + A.c20 * B.c2 + A.c21 * B.c3 - A.c22 * B.c7 -
                      A.c23 * B.c8 - A.c24 * B.c9 + A.c25 * B.c0 - A.c31 * B.c10;
    ctype const c26 = A.c0 * B.c11 - A.c6 * B.c14 + A.c7 * B.c9 - A.c8 * B.c8 +
                      A.c9 * B.c10 - A.c10 * B.c13 + A.c11 * B.c12 + A.c12 * B.c15 -
                      A.c13 * B.c3 + A.c14 * B.c2 + A.c15 * B.c4 + A.c26 * B.c0 -
                      A.c27 * B.c6 + A.c28 * B.c5 + A.c29 * B.c1 - A.c30 * B.c7;
    ctype const c27 = A.c0 * B.c12 - A.c6 * B.c9 - A.c7 * B.c14 + A.c8 * B.c7 +
                      A.c9 * B.c13 + A.c10 * B.c10 - A.c11 * B.c11 + A.c12 * B.c3 +
                      A.c13 * B.c15 - A.c14 * B.c1 + A.c15 * B.c5 + A.c26 * B.c6 +
                      A.c27 * B.c0 - A.c28 * B.c4 + A.c29 * B.c2 - A.c30 * B.c8;
    ctype const c28 = A.c0 * B.c13 + A.c6 * B.c8 - A.c7 * B.c7 - A.c8 * B.c14 -
                      A.c9 * B.c12 + A.c10 * B.c11 + A.c11 * B.c10 - A.c12 * B.c2 +
                      A.c13 * B.c1 + A.c14 * B.c15 + A.c15 * B.c6 - A.c26 * B.c5 +
                      A.c27 * B.c4 + A.c28 * B.c0 + A.c29 * B.c3 - A.c30 * B.c9;
    ctype const c29 = A.c0 * B.c14 - A.c9 * B.c7 - A.c10 * B.c8 - A.c11 * B.c9 -
                      A.c12 * B.c4 + A.c12 * B.c11 - A.c13 * B.c5 + A.c13 * B.c12 -
                      A.c14 * B.c6 + A.c14 * B.c13 + A.c15 * B.c14 - A.c26 * B.c7 -
                      A.c27 * B.c8 - A.c28 * B.c9 + A.c29 * B.c0 - A.c29 * B.c10;
    ctype const c30 = A.c0 * B.c15 - A.c6 * B.c4 - A.c6 * B.c11 - A.c7 * B.c5 -
                      A.c7 * B.c12 - A.c8 * B.c6 - A.c8 * B.c13 - A.c9 * B.c1 -
                      A.c10 * B.c2 - A.c11 * B.c3 - A.c15 * B.c15 + A.c26 * B.c1 +
                      A.c27 * B.c2 + A.c28 * B.c3 + A.c30 * B.c0 + A.c30 * B.c10;
    ctype const c31 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 +
                      A.c5 * B.c15 - A.c16 * B.c4 - A.c17 * B.c5 - A.c18 * B.c6 -
                      A.c19 * B.c1 - A.c20 * B.c2 - A.c21 * B.c3 - A.c22 * B.c7 -
                      A.c23 * B.c8 - A.c24 * B.c9 - A.c25 * B.c10 + A.c31 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& A,
                                                      MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c12 - A.c2 * B.c13 - A.c3 * B.c14 -
                     A.c4 * B.c9 - A.c5 * B.c10 - A.c6 * B.c11 - A.c7 * B.c6 -
                     A.c8 * B.c7 - A.c9 * B.c8 + A.c10 * B.c15 - A.c11 * B.c26 -
                     A.c12 * B.c27 - A.c13 * B.c28 + A.c14 * B.c30 + A.c15 * B.c29;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c5 - A.c2 * B.c21 + A.c3 * B.c20 +
                     A.c4 * B.c25 - A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c4 +
                     A.c8 * B.c24 - A.c9 * B.c23 - A.c10 * B.c16 - A.c11 * B.c31 +
                     A.c12 * B.c18 - A.c13 * B.c17 - A.c14 * B.c22 + A.c15 * B.c19;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c21 + A.c2 * B.c5 - A.c3 * B.c19 +
                     A.c4 * B.c3 + A.c5 * B.c25 - A.c6 * B.c1 - A.c7 * B.c24 -
                     A.c8 * B.c4 + A.c9 * B.c22 - A.c10 * B.c17 - A.c11 * B.c18 -
                     A.c12 * B.c31 + A.c13 * B.c16 - A.c14 * B.c23 + A.c15 * B.c20;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c20 + A.c2 * B.c19 + A.c3 * B.c5 -
                     A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c25 + A.c7 * B.c23 -
                     A.c8 * B.c22 - A.c9 * B.c4 - A.c10 * B.c18 + A.c11 * B.c17 -
                     A.c12 * B.c16 - A.c13 * B.c31 - A.c14 * B.c24 + A.c15 * B.c21;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c1 - A.c1 * B.c16 + A.c2 * B.c2 -
                     A.c2 * B.c17 + A.c3 * B.c3 - A.c3 * B.c18 - A.c4 * B.c22 -
                     A.c5 * B.c23 - A.c6 * B.c24 - A.c10 * B.c4 + A.c11 * B.c22 +
                     A.c12 * B.c23 + A.c13 * B.c24 - A.c15 * B.c25 + A.c15 * B.c31;
    ctype const c5 = A.c0 * B.c5 - A.c4 * B.c19 - A.c5 * B.c20 - A.c6 * B.c21 -
                     A.c7 * B.c1 - A.c7 * B.c16 - A.c8 * B.c2 - A.c8 * B.c17 -
                     A.c9 * B.c3 - A.c9 * B.c18 + A.c10 * B.c5 - A.c11 * B.c19 -
                     A.c12 * B.c20 - A.c13 * B.c21 + A.c14 * B.c25 + A.c14 * B.c31;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c0 + A.c1 * B.c15 - A.c2 * B.c11 -
                     A.c2 * B.c28 + A.c3 * B.c10 + A.c3 * B.c27 + A.c4 * B.c30 -
                     A.c5 * B.c8 + A.c6 * B.c7 - A.c10 * B.c6 - A.c11 * B.c30 +
                     A.c12 * B.c8 - A.c13 * B.c7 + A.c15 * B.c9 + A.c15 * B.c26;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c11 + A.c1 * B.c28 + A.c2 * B.c0 +
                     A.c2 * B.c15 - A.c3 * B.c9 - A.c3 * B.c26 + A.c4 * B.c8 +
                     A.c5 * B.c30 - A.c6 * B.c6 - A.c10 * B.c7 - A.c11 * B.c8 -
                     A.c12 * B.c30 + A.c13 * B.c6 + A.c15 * B.c10 + A.c15 * B.c27;
    ctype const c8 = A.c0 * B.c8 - A.c1 * B.c10 - A.c1 * B.c27 + A.c2 * B.c9 +
                     A.c2 * B.c26 + A.c3 * B.c0 + A.c3 * B.c15 - A.c4 * B.c7 +
                     A.c5 * B.c6 + A.c6 * B.c30 - A.c10 * B.c8 + A.c11 * B.c7 -
                     A.c12 * B.c6 - A.c13 * B.c30 + A.c15 * B.c11 + A.c15 * B.c28;
    ctype const c9 = A.c0 * B.c9 + A.c1 * B.c29 - A.c2 * B.c14 + A.c3 * B.c13 +
                     A.c4 * B.c0 - A.c5 * B.c11 + A.c6 * B.c10 + A.c7 * B.c30 -
                     A.c8 * B.c8 + A.c9 * B.c7 + A.c10 * B.c26 + A.c11 * B.c15 -
                     A.c12 * B.c28 + A.c13 * B.c27 + A.c14 * B.c6 + A.c15 * B.c12;
    ctype const c10 = A.c0 * B.c10 + A.c1 * B.c14 + A.c2 * B.c29 - A.c3 * B.c12 +
                      A.c4 * B.c11 + A.c5 * B.c0 - A.c6 * B.c9 + A.c7 * B.c8 +
                      A.c8 * B.c30 - A.c9 * B.c6 + A.c10 * B.c27 + A.c11 * B.c28 +
                      A.c12 * B.c15 - A.c13 * B.c26 + A.c14 * B.c7 + A.c15 * B.c13;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c13 + A.c2 * B.c12 + A.c3 * B.c29 -
                      A.c4 * B.c10 + A.c5 * B.c9 + A.c6 * B.c0 - A.c7 * B.c7 +
                      A.c8 * B.c6 + A.c9 * B.c30 + A.c10 * B.c28 - A.c11 * B.c27 +
                      A.c12 * B.c26 + A.c13 * B.c15 + A.c14 * B.c8 + A.c15 * B.c14;
    ctype const c12 = A.c0 * B.c12 + A.c4 * B.c29 - A.c5 * B.c14 + A.c6 * B.c13 +
                      A.c7 * B.c0 - A.c7 * B.c15 - A.c8 * B.c11 + A.c8 * B.c28 +
                      A.c9 * B.c10 - A.c9 * B.c27 + A.c10 * B.c12 + A.c11 * B.c29 -
                      A.c12 * B.c14 + A.c13 * B.c13 + A.c14 * B.c9 - A.c14 * B.c26;
    ctype const c13 = A.c0 * B.c13 + A.c4 * B.c14 + A.c5 * B.c29 - A.c6 * B.c12 +
                      A.c7 * B.c11 - A.c7 * B.c28 + A.c8 * B.c0 - A.c8 * B.c15 -
                      A.c9 * B.c9 + A.c9 * B.c26 + A.c10 * B.c13 + A.c11 * B.c14 +
                      A.c12 * B.c29 - A.c13 * B.c12 + A.c14 * B.c10 - A.c14 * B.c27;
    ctype const c14 = A.c0 * B.c14 - A.c4 * B.c13 + A.c5 * B.c12 + A.c6 * B.c29 -
                      A.c7 * B.c10 + A.c7 * B.c27 + A.c8 * B.c9 - A.c8 * B.c26 +
                      A.c9 * B.c0 - A.c9 * B.c15 + A.c10 * B.c14 - A.c11 * B.c13 +
                      A.c12 * B.c12 + A.c13 * B.c29 + A.c14 * B.c11 - A.c14 * B.c28;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 -
                      A.c4 * B.c26 - A.c5 * B.c27 - A.c6 * B.c28 - A.c7 * B.c6 -
                      A.c8 * B.c7 - A.c9 * B.c8 + A.c10 * B.c0 - A.c11 * B.c9 -
                      A.c12 * B.c10 - A.c13 * B.c11 + A.c14 * B.c30 - A.c15 * B.c29;
    ctype const c16 = A.c0 * B.c16 + A.c1 * B.c5 - A.c2 * B.c21 + A.c3 * B.c20 +
                      A.c4 * B.c31 - A.c5 * B.c18 + A.c6 * B.c17 + A.c7 * B.c4 -
                      A.c8 * B.c24 + A.c9 * B.c23 - A.c10 * B.c1 - A.c11 * B.c25 +
                      A.c12 * B.c3 - A.c13 * B.c2 + A.c14 * B.c22 + A.c15 * B.c19;
    ctype const c17 = A.c0 * B.c17 + A.c1 * B.c21 + A.c2 * B.c5 - A.c3 * B.c19 +
                      A.c4 * B.c18 + A.c5 * B.c31 - A.c6 * B.c16 + A.c7 * B.c24 +
                      A.c8 * B.c4 - A.c9 * B.c22 - A.c10 * B.c2 - A.c11 * B.c3 -
                      A.c12 * B.c25 + A.c13 * B.c1 + A.c14 * B.c23 + A.c15 * B.c20;
    ctype const c18 = A.c0 * B.c18 - A.c1 * B.c20 + A.c2 * B.c19 + A.c3 * B.c5 -
                      A.c4 * B.c17 + A.c5 * B.c16 + A.c6 * B.c31 - A.c7 * B.c23 +
                      A.c8 * B.c22 + A.c9 * B.c4 - A.c10 * B.c3 + A.c11 * B.c2 -
                      A.c12 * B.c1 - A.c13 * B.c25 + A.c14 * B.c24 + A.c15 * B.c21;
    ctype const c19 = A.c0 * B.c19 + A.c4 * B.c5 - A.c5 * B.c21 + A.c6 * B.c20 +
                      A.c7 * B.c25 + A.c7 * B.c31 - A.c8 * B.c3 - A.c8 * B.c18 +
                      A.c9 * B.c2 + A.c9 * B.c17 + A.c10 * B.c19 + A.c11 * B.c5 -
                      A.c12 * B.c21 + A.c13 * B.c20 + A.c14 * B.c1 + A.c14 * B.c16;
    ctype const c20 = A.c0 * B.c20 + A.c4 * B.c21 + A.c5 * B.c5 - A.c6 * B.c19 +
                      A.c7 * B.c3 + A.c7 * B.c18 + A.c8 * B.c25 + A.c8 * B.c31 -
                      A.c9 * B.c1 - A.c9 * B.c16 + A.c10 * B.c20 + A.c11 * B.c21 +
                      A.c12 * B.c5 - A.c13 * B.c19 + A.c14 * B.c2 + A.c14 * B.c17;
    ctype const c21 = A.c0 * B.c21 - A.c4 * B.c20 + A.c5 * B.c19 + A.c6 * B.c5 -
                      A.c7 * B.c2 - A.c7 * B.c17 + A.c8 * B.c1 + A.c8 * B.c16 +
                      A.c9 * B.c25 + A.c9 * B.c31 + A.c10 * B.c21 - A.c11 * B.c20 +
                      A.c12 * B.c19 + A.c13 * B.c5 + A.c14 * B.c3 + A.c14 * B.c18;
    ctype const c22 = A.c0 * B.c22 - A.c1 * B.c25 + A.c1 * B.c31 + A.c2 * B.c3 -
                      A.c2 * B.c18 - A.c3 * B.c2 + A.c3 * B.c17 + A.c4 * B.c4 -
                      A.c5 * B.c24 + A.c6 * B.c23 - A.c10 * B.c22 - A.c11 * B.c4 +
                      A.c12 * B.c24 - A.c13 * B.c23 - A.c15 * B.c1 + A.c15 * B.c16;
    ctype const c23 = A.c0 * B.c23 - A.c1 * B.c3 + A.c1 * B.c18 - A.c2 * B.c25 +
                      A.c2 * B.c31 + A.c3 * B.c1 - A.c3 * B.c16 + A.c4 * B.c24 +
                      A.c5 * B.c4 - A.c6 * B.c22 - A.c10 * B.c23 - A.c11 * B.c24 -
                      A.c12 * B.c4 + A.c13 * B.c22 - A.c15 * B.c2 + A.c15 * B.c17;
    ctype const c24 = A.c0 * B.c24 + A.c1 * B.c2 - A.c1 * B.c17 - A.c2 * B.c1 +
                      A.c2 * B.c16 - A.c3 * B.c25 + A.c3 * B.c31 - A.c4 * B.c23 +
                      A.c5 * B.c22 + A.c6 * B.c4 - A.c10 * B.c24 + A.c11 * B.c23 -
                      A.c12 * B.c22 - A.c13 * B.c4 - A.c15 * B.c3 + A.c15 * B.c18;
    ctype const c25 = A.c0 * B.c25 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                      A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c22 +
                      A.c8 * B.c23 + A.c9 * B.c24 - A.c10 * B.c31 + A.c11 * B.c16 +
                      A.c12 * B.c17 + A.c13 * B.c18 - A.c14 * B.c4 + A.c15 * B.c5;
    ctype const c26 = A.c0 * B.c26 - A.c1 * B.c29 + A.c2 * B.c14 - A.c3 * B.c13 +
                      A.c4 * B.c15 - A.c5 * B.c28 + A.c6 * B.c27 + A.c7 * B.c30 -
                      A.c8 * B.c8 + A.c9 * B.c7 + A.c10 * B.c9 + A.c11 * B.c0 -
                      A.c12 * B.c11 + A.c13 * B.c10 + A.c14 * B.c6 - A.c15 * B.c12;
    ctype const c27 = A.c0 * B.c27 - A.c1 * B.c14 - A.c2 * B.c29 + A.c3 * B.c12 +
                      A.c4 * B.c28 + A.c5 * B.c15 - A.c6 * B.c26 + A.c7 * B.c8 +
                      A.c8 * B.c30 - A.c9 * B.c6 + A.c10 * B.c10 + A.c11 * B.c11 +
                      A.c12 * B.c0 - A.c13 * B.c9 + A.c14 * B.c7 - A.c15 * B.c13;
    ctype const c28 = A.c0 * B.c28 + A.c1 * B.c13 - A.c2 * B.c12 - A.c3 * B.c29 -
                      A.c4 * B.c27 + A.c5 * B.c26 + A.c6 * B.c15 - A.c7 * B.c7 +
                      A.c8 * B.c6 + A.c9 * B.c30 + A.c10 * B.c11 - A.c11 * B.c10 +
                      A.c12 * B.c9 + A.c13 * B.c0 + A.c14 * B.c8 - A.c15 * B.c14;
    ctype const c29 = A.c0 * B.c29 - A.c4 * B.c12 - A.c5 * B.c13 - A.c6 * B.c14 -
                      A.c7 * B.c9 + A.c7 * B.c26 - A.c8 * B.c10 + A.c8 * B.c27 -
                      A.c9 * B.c11 + A.c9 * B.c28 + A.c10 * B.c29 - A.c11 * B.c12 -
                      A.c12 * B.c13 - A.c13 * B.c14 + A.c14 * B.c0 - A.c14 * B.c15;
    ctype const c30 = A.c0 * B.c30 - A.c1 * B.c9 - A.c1 * B.c26 - A.c2 * B.c10 -
                      A.c2 * B.c27 - A.c3 * B.c11 - A.c3 * B.c28 - A.c4 * B.c6 -
                      A.c5 * B.c7 - A.c6 * B.c8 - A.c10 * B.c30 + A.c11 * B.c6 +
                      A.c12 * B.c7 + A.c13 * B.c8 + A.c15 * B.c0 + A.c15 * B.c15;
    ctype const c31 = A.c0 * B.c31 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                      A.c4 * B.c16 - A.c5 * B.c17 - A.c6 * B.c18 - A.c7 * B.c22 -
                      A.c8 * B.c23 - A.c9 * B.c24 - A.c10 * B.c25 + A.c11 * B.c1 +
                      A.c12 * B.c2 + A.c13 * B.c3 + A.c14 * B.c4 + A.c15 * B.c5;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& A,
                                                      MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c4 * B.c4 - A.c5 * B.c3 +
                     A.c16 * B.c5 + A.c17 * B.c6 + A.c18 * B.c7 + A.c19 * B.c11 +
                     A.c20 * B.c12 + A.c21 * B.c13 + A.c22 * B.c8 + A.c23 * B.c9 +
                     A.c24 * B.c10 - A.c25 * B.c14 - A.c31 * B.c15;
    ctype const c1 = A.c0 * B.c0 + A.c6 * B.c4 - A.c7 * B.c10 + A.c8 * B.c9 +
                     A.c9 * B.c14 - A.c10 * B.c2 + A.c11 * B.c1 - A.c12 * B.c3 +
                     A.c13 * B.c13 - A.c14 * B.c12 - A.c15 * B.c5 - A.c26 * B.c15 +
                     A.c27 * B.c7 - A.c28 * B.c6 - A.c29 * B.c11 + A.c30 * B.c8;
    ctype const c2 = A.c0 * B.c1 + A.c6 * B.c10 + A.c7 * B.c4 - A.c8 * B.c8 +
                     A.c9 * B.c2 + A.c10 * B.c14 - A.c11 * B.c0 - A.c12 * B.c13 -
                     A.c13 * B.c3 + A.c14 * B.c11 - A.c15 * B.c6 - A.c26 * B.c7 -
                     A.c27 * B.c15 + A.c28 * B.c5 - A.c29 * B.c12 + A.c30 * B.c9;
    ctype const c3 = A.c0 * B.c2 - A.c6 * B.c9 + A.c7 * B.c8 + A.c8 * B.c4 - A.c9 * B.c1 +
                     A.c10 * B.c0 + A.c11 * B.c14 + A.c12 * B.c12 - A.c13 * B.c11 -
                     A.c14 * B.c3 - A.c15 * B.c7 + A.c26 * B.c6 - A.c27 * B.c5 -
                     A.c28 * B.c15 - A.c29 * B.c13 + A.c30 * B.c10;
    ctype const c4 = A.c0 * B.c3 + A.c6 * B.c0 - A.c6 * B.c5 + A.c7 * B.c1 - A.c7 * B.c6 +
                     A.c8 * B.c2 - A.c8 * B.c7 - A.c9 * B.c11 - A.c10 * B.c12 -
                     A.c11 * B.c13 - A.c15 * B.c3 + A.c26 * B.c11 + A.c27 * B.c12 +
                     A.c28 * B.c13 - A.c30 * B.c14 + A.c30 * B.c15;
    ctype const c5 = A.c0 * B.c4 - A.c9 * B.c8 - A.c10 * B.c9 - A.c11 * B.c10 -
                     A.c12 * B.c0 - A.c12 * B.c5 - A.c13 * B.c1 - A.c13 * B.c6 -
                     A.c14 * B.c2 - A.c14 * B.c7 + A.c15 * B.c4 - A.c26 * B.c8 -
                     A.c27 * B.c9 - A.c28 * B.c10 + A.c29 * B.c14 + A.c29 * B.c15;
    ctype const c6 = -A.c1 * B.c3 + A.c2 * B.c13 - A.c3 * B.c12 + A.c4 * B.c0 -
                     A.c4 * B.c5 - A.c16 * B.c3 + A.c17 * B.c13 - A.c18 * B.c12 +
                     A.c22 * B.c14 - A.c22 * B.c15 - A.c23 * B.c2 + A.c23 * B.c7 +
                     A.c24 * B.c1 - A.c24 * B.c6 - A.c25 * B.c11 - A.c31 * B.c11;
    ctype const c7 = -A.c1 * B.c13 - A.c2 * B.c3 + A.c3 * B.c11 + A.c4 * B.c1 -
                     A.c4 * B.c6 - A.c16 * B.c13 - A.c17 * B.c3 + A.c18 * B.c11 +
                     A.c22 * B.c2 - A.c22 * B.c7 + A.c23 * B.c14 - A.c23 * B.c15 -
                     A.c24 * B.c0 + A.c24 * B.c5 - A.c25 * B.c12 - A.c31 * B.c12;
    ctype const c8 = A.c1 * B.c12 - A.c2 * B.c11 - A.c3 * B.c3 + A.c4 * B.c2 -
                     A.c4 * B.c7 + A.c16 * B.c12 - A.c17 * B.c11 - A.c18 * B.c3 -
                     A.c22 * B.c1 + A.c22 * B.c6 + A.c23 * B.c0 - A.c23 * B.c5 +
                     A.c24 * B.c14 - A.c24 * B.c15 - A.c25 * B.c13 - A.c31 * B.c13;
    ctype const c9 = -A.c1 * B.c14 + A.c2 * B.c2 - A.c3 * B.c1 - A.c4 * B.c8 -
                     A.c5 * B.c11 - A.c16 * B.c15 + A.c17 * B.c7 - A.c18 * B.c6 -
                     A.c19 * B.c3 + A.c20 * B.c13 - A.c21 * B.c12 - A.c22 * B.c4 +
                     A.c23 * B.c10 - A.c24 * B.c9 - A.c25 * B.c0 - A.c31 * B.c5;
    ctype const c10 = -A.c1 * B.c2 - A.c2 * B.c14 + A.c3 * B.c0 - A.c4 * B.c9 -
                      A.c5 * B.c12 - A.c16 * B.c7 - A.c17 * B.c15 + A.c18 * B.c5 -
                      A.c19 * B.c13 - A.c20 * B.c3 + A.c21 * B.c11 - A.c22 * B.c10 -
                      A.c23 * B.c4 + A.c24 * B.c8 - A.c25 * B.c1 - A.c31 * B.c6;
    ctype const c11 = A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c14 - A.c4 * B.c10 -
                      A.c5 * B.c13 + A.c16 * B.c6 - A.c17 * B.c5 - A.c18 * B.c15 +
                      A.c19 * B.c12 - A.c20 * B.c11 - A.c21 * B.c3 + A.c22 * B.c9 -
                      A.c23 * B.c8 - A.c24 * B.c4 - A.c25 * B.c2 - A.c31 * B.c7;
    ctype const c12 = A.c1 * B.c4 - A.c2 * B.c10 + A.c3 * B.c9 - A.c5 * B.c0 -
                      A.c5 * B.c5 - A.c16 * B.c4 + A.c17 * B.c10 - A.c18 * B.c9 -
                      A.c19 * B.c14 - A.c19 * B.c15 + A.c20 * B.c2 + A.c20 * B.c7 -
                      A.c21 * B.c1 - A.c21 * B.c6 + A.c25 * B.c8 - A.c31 * B.c8;
    ctype const c13 = A.c1 * B.c10 + A.c2 * B.c4 - A.c3 * B.c8 - A.c5 * B.c1 -
                      A.c5 * B.c6 - A.c16 * B.c10 - A.c17 * B.c4 + A.c18 * B.c8 -
                      A.c19 * B.c2 - A.c19 * B.c7 - A.c20 * B.c14 - A.c20 * B.c15 +
                      A.c21 * B.c0 + A.c21 * B.c5 + A.c25 * B.c9 - A.c31 * B.c9;
    ctype const c14 = -A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c4 - A.c5 * B.c2 -
                      A.c5 * B.c7 + A.c16 * B.c9 - A.c17 * B.c8 - A.c18 * B.c4 +
                      A.c19 * B.c1 + A.c19 * B.c6 - A.c20 * B.c0 - A.c20 * B.c5 -
                      A.c21 * B.c14 - A.c21 * B.c15 + A.c25 * B.c10 - A.c31 * B.c10;
    ctype const c15 = -A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 + A.c4 * B.c4 -
                      A.c5 * B.c3 - A.c16 * B.c0 - A.c17 * B.c1 - A.c18 * B.c2 +
                      A.c19 * B.c11 + A.c20 * B.c12 + A.c21 * B.c13 - A.c22 * B.c8 -
                      A.c23 * B.c9 - A.c24 * B.c10 + A.c25 * B.c15 + A.c31 * B.c14;
    ctype const c16 = A.c0 * B.c5 + A.c6 * B.c4 - A.c7 * B.c10 + A.c8 * B.c9 +
                      A.c9 * B.c15 - A.c10 * B.c7 + A.c11 * B.c6 + A.c12 * B.c3 -
                      A.c13 * B.c13 + A.c14 * B.c12 - A.c15 * B.c0 - A.c26 * B.c14 +
                      A.c27 * B.c2 - A.c28 * B.c1 + A.c29 * B.c11 + A.c30 * B.c8;
    ctype const c17 = A.c0 * B.c6 + A.c6 * B.c10 + A.c7 * B.c4 - A.c8 * B.c8 +
                      A.c9 * B.c7 + A.c10 * B.c15 - A.c11 * B.c5 + A.c12 * B.c13 +
                      A.c13 * B.c3 - A.c14 * B.c11 - A.c15 * B.c1 - A.c26 * B.c2 -
                      A.c27 * B.c14 + A.c28 * B.c0 + A.c29 * B.c12 + A.c30 * B.c9;
    ctype const c18 = A.c0 * B.c7 - A.c6 * B.c9 + A.c7 * B.c8 + A.c8 * B.c4 -
                      A.c9 * B.c6 + A.c10 * B.c5 + A.c11 * B.c15 - A.c12 * B.c12 +
                      A.c13 * B.c11 + A.c14 * B.c3 - A.c15 * B.c2 + A.c26 * B.c1 -
                      A.c27 * B.c0 - A.c28 * B.c14 + A.c29 * B.c13 + A.c30 * B.c10;
    ctype const c19 = A.c0 * B.c8 + A.c9 * B.c4 - A.c10 * B.c10 + A.c11 * B.c9 +
                      A.c12 * B.c14 + A.c12 * B.c15 - A.c13 * B.c2 - A.c13 * B.c7 +
                      A.c14 * B.c1 + A.c14 * B.c6 + A.c15 * B.c8 + A.c26 * B.c4 -
                      A.c27 * B.c10 + A.c28 * B.c9 + A.c29 * B.c0 + A.c29 * B.c5;
    ctype const c20 = A.c0 * B.c9 + A.c9 * B.c10 + A.c10 * B.c4 - A.c11 * B.c8 +
                      A.c12 * B.c2 + A.c12 * B.c7 + A.c13 * B.c14 + A.c13 * B.c15 -
                      A.c14 * B.c0 - A.c14 * B.c5 + A.c15 * B.c9 + A.c26 * B.c10 +
                      A.c27 * B.c4 - A.c28 * B.c8 + A.c29 * B.c1 + A.c29 * B.c6;
    ctype const c21 = A.c0 * B.c10 - A.c9 * B.c9 + A.c10 * B.c8 + A.c11 * B.c4 -
                      A.c12 * B.c1 - A.c12 * B.c6 + A.c13 * B.c0 + A.c13 * B.c5 +
                      A.c14 * B.c14 + A.c14 * B.c15 + A.c15 * B.c10 - A.c26 * B.c9 +
                      A.c27 * B.c8 + A.c28 * B.c4 + A.c29 * B.c2 + A.c29 * B.c7;
    ctype const c22 = A.c0 * B.c11 - A.c6 * B.c14 + A.c6 * B.c15 + A.c7 * B.c2 -
                      A.c7 * B.c7 - A.c8 * B.c1 + A.c8 * B.c6 + A.c9 * B.c3 -
                      A.c10 * B.c13 + A.c11 * B.c12 - A.c15 * B.c11 - A.c26 * B.c3 +
                      A.c27 * B.c13 - A.c28 * B.c12 - A.c30 * B.c0 + A.c30 * B.c5;
    ctype const c23 = A.c0 * B.c12 - A.c6 * B.c2 + A.c6 * B.c7 - A.c7 * B.c14 +
                      A.c7 * B.c15 + A.c8 * B.c0 - A.c8 * B.c5 + A.c9 * B.c13 +
                      A.c10 * B.c3 - A.c11 * B.c11 - A.c15 * B.c12 - A.c26 * B.c13 -
                      A.c27 * B.c3 + A.c28 * B.c11 - A.c30 * B.c1 + A.c30 * B.c6;
    ctype const c24 = A.c0 * B.c13 + A.c6 * B.c1 - A.c6 * B.c6 - A.c7 * B.c0 +
                      A.c7 * B.c5 - A.c8 * B.c14 + A.c8 * B.c15 - A.c9 * B.c12 +
                      A.c10 * B.c11 + A.c11 * B.c3 - A.c15 * B.c13 + A.c26 * B.c12 -
                      A.c27 * B.c11 - A.c28 * B.c3 - A.c30 * B.c2 + A.c30 * B.c7;
    ctype const c25 = A.c0 * B.c14 - A.c6 * B.c8 - A.c7 * B.c9 - A.c8 * B.c10 -
                      A.c9 * B.c0 - A.c10 * B.c1 - A.c11 * B.c2 + A.c12 * B.c11 +
                      A.c13 * B.c12 + A.c14 * B.c13 - A.c15 * B.c15 + A.c26 * B.c5 +
                      A.c27 * B.c6 + A.c28 * B.c7 - A.c29 * B.c3 + A.c30 * B.c4;
    ctype const c26 = A.c1 * B.c15 - A.c2 * B.c7 + A.c3 * B.c6 + A.c4 * B.c8 -
                      A.c5 * B.c11 + A.c16 * B.c14 - A.c17 * B.c2 + A.c18 * B.c1 -
                      A.c19 * B.c3 + A.c20 * B.c13 - A.c21 * B.c12 + A.c22 * B.c4 -
                      A.c23 * B.c10 + A.c24 * B.c9 + A.c25 * B.c5 + A.c31 * B.c0;
    ctype const c27 = A.c1 * B.c7 + A.c2 * B.c15 - A.c3 * B.c5 + A.c4 * B.c9 -
                      A.c5 * B.c12 + A.c16 * B.c2 + A.c17 * B.c14 - A.c18 * B.c0 -
                      A.c19 * B.c13 - A.c20 * B.c3 + A.c21 * B.c11 + A.c22 * B.c10 +
                      A.c23 * B.c4 - A.c24 * B.c8 + A.c25 * B.c6 + A.c31 * B.c1;
    ctype const c28 = -A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c15 + A.c4 * B.c10 -
                      A.c5 * B.c13 - A.c16 * B.c1 + A.c17 * B.c0 + A.c18 * B.c14 +
                      A.c19 * B.c12 - A.c20 * B.c11 - A.c21 * B.c3 - A.c22 * B.c9 +
                      A.c23 * B.c8 + A.c24 * B.c4 + A.c25 * B.c7 + A.c31 * B.c2;
    ctype const c29 = -A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 - A.c5 * B.c14 -
                      A.c5 * B.c15 + A.c16 * B.c8 + A.c17 * B.c9 + A.c18 * B.c10 +
                      A.c19 * B.c0 + A.c19 * B.c5 + A.c20 * B.c1 + A.c20 * B.c6 +
                      A.c21 * B.c2 + A.c21 * B.c7 + A.c25 * B.c4 - A.c31 * B.c4;
    ctype const c30 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
                      A.c4 * B.c15 + A.c16 * B.c11 + A.c17 * B.c12 + A.c18 * B.c13 -
                      A.c22 * B.c0 + A.c22 * B.c5 - A.c23 * B.c1 + A.c23 * B.c6 -
                      A.c24 * B.c2 + A.c24 * B.c7 - A.c25 * B.c3 - A.c31 * B.c3;
    ctype const c31 = A.c0 * B.c15 - A.c6 * B.c8 - A.c7 * B.c9 - A.c8 * B.c10 -
                      A.c9 * B.c5 - A.c10 * B.c6 - A.c11 * B.c7 - A.c12 * B.c11 -
                      A.c13 * B.c12 - A.c14 * B.c13 - A.c15 * B.c14 + A.c26 * B.c0 +
                      A.c27 * B.c1 + A.c28 * B.c2 + A.c29 * B.c3 + A.c30 * B.c4;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& A,
                                                      MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c3 * B.c5 - A.c4 * B.c4 +
                     A.c5 * B.c16 + A.c6 * B.c17 + A.c7 * B.c18 + A.c8 * B.c22 +
                     A.c9 * B.c23 + A.c10 * B.c24 + A.c11 * B.c19 + A.c12 * B.c20 +
                     A.c13 * B.c21 - A.c14 * B.c25 - A.c15 * B.c31;
    ctype const c1 = A.c0 * B.c0 - A.c1 * B.c11 + A.c2 * B.c10 + A.c3 * B.c12 -
                     A.c4 * B.c6 - A.c5 * B.c15 + A.c6 * B.c28 - A.c7 * B.c27 -
                     A.c8 * B.c30 + A.c9 * B.c8 - A.c10 * B.c7 + A.c11 * B.c29 -
                     A.c12 * B.c14 + A.c13 * B.c13 + A.c14 * B.c9 - A.c15 * B.c26;
    ctype const c2 = A.c0 * B.c11 + A.c1 * B.c0 - A.c2 * B.c9 + A.c3 * B.c13 -
                     A.c4 * B.c7 - A.c5 * B.c28 - A.c6 * B.c15 + A.c7 * B.c26 -
                     A.c8 * B.c8 - A.c9 * B.c30 + A.c10 * B.c6 + A.c11 * B.c14 +
                     A.c12 * B.c29 - A.c13 * B.c12 + A.c14 * B.c10 - A.c15 * B.c27;
    ctype const c3 = -A.c0 * B.c10 + A.c1 * B.c9 + A.c2 * B.c0 + A.c3 * B.c14 -
                     A.c4 * B.c8 + A.c5 * B.c27 - A.c6 * B.c26 - A.c7 * B.c15 +
                     A.c8 * B.c7 - A.c9 * B.c6 - A.c10 * B.c30 - A.c11 * B.c13 +
                     A.c12 * B.c12 + A.c13 * B.c29 + A.c14 * B.c11 - A.c15 * B.c28;
    ctype const c4 = -A.c0 * B.c6 - A.c1 * B.c7 - A.c2 * B.c8 + A.c3 * B.c0 +
                     A.c3 * B.c15 - A.c5 * B.c6 - A.c6 * B.c7 - A.c7 * B.c8 -
                     A.c11 * B.c9 - A.c11 * B.c26 - A.c12 * B.c10 - A.c12 * B.c27 -
                     A.c13 * B.c11 - A.c13 * B.c28 + A.c14 * B.c30 + A.c15 * B.c30;
    ctype const c5 = A.c0 * B.c12 + A.c1 * B.c13 + A.c2 * B.c14 + A.c4 * B.c0 -
                     A.c4 * B.c15 - A.c5 * B.c12 - A.c6 * B.c13 - A.c7 * B.c14 -
                     A.c8 * B.c9 + A.c8 * B.c26 - A.c9 * B.c10 + A.c9 * B.c27 -
                     A.c10 * B.c11 + A.c10 * B.c28 - A.c14 * B.c29 + A.c15 * B.c29;
    ctype const c6 = -A.c0 * B.c4 + A.c1 * B.c24 - A.c2 * B.c23 + A.c3 * B.c1 -
                     A.c3 * B.c16 - A.c5 * B.c4 + A.c6 * B.c24 - A.c7 * B.c23 +
                     A.c11 * B.c25 - A.c11 * B.c31 - A.c12 * B.c3 + A.c12 * B.c18 +
                     A.c13 * B.c2 - A.c13 * B.c17 - A.c14 * B.c22 - A.c15 * B.c22;
    ctype const c7 = -A.c0 * B.c24 - A.c1 * B.c4 + A.c2 * B.c22 + A.c3 * B.c2 -
                     A.c3 * B.c17 - A.c5 * B.c24 - A.c6 * B.c4 + A.c7 * B.c22 +
                     A.c11 * B.c3 - A.c11 * B.c18 + A.c12 * B.c25 - A.c12 * B.c31 -
                     A.c13 * B.c1 + A.c13 * B.c16 - A.c14 * B.c23 - A.c15 * B.c23;
    ctype const c8 = A.c0 * B.c23 - A.c1 * B.c22 - A.c2 * B.c4 + A.c3 * B.c3 -
                     A.c3 * B.c18 + A.c5 * B.c23 - A.c6 * B.c22 - A.c7 * B.c4 -
                     A.c11 * B.c2 + A.c11 * B.c17 + A.c12 * B.c1 - A.c12 * B.c16 +
                     A.c13 * B.c25 - A.c13 * B.c31 - A.c14 * B.c24 - A.c15 * B.c24;
    ctype const c9 = -A.c0 * B.c25 + A.c1 * B.c3 - A.c2 * B.c2 - A.c3 * B.c19 -
                     A.c4 * B.c22 - A.c5 * B.c31 + A.c6 * B.c18 - A.c7 * B.c17 -
                     A.c8 * B.c4 + A.c9 * B.c24 - A.c10 * B.c23 - A.c11 * B.c5 +
                     A.c12 * B.c21 - A.c13 * B.c20 - A.c14 * B.c1 - A.c15 * B.c16;
    ctype const c10 = -A.c0 * B.c3 - A.c1 * B.c25 + A.c2 * B.c1 - A.c3 * B.c20 -
                      A.c4 * B.c23 - A.c5 * B.c18 - A.c6 * B.c31 + A.c7 * B.c16 -
                      A.c8 * B.c24 - A.c9 * B.c4 + A.c10 * B.c22 - A.c11 * B.c21 -
                      A.c12 * B.c5 + A.c13 * B.c19 - A.c14 * B.c2 - A.c15 * B.c17;
    ctype const c11 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c25 - A.c3 * B.c21 -
                      A.c4 * B.c24 + A.c5 * B.c17 - A.c6 * B.c16 - A.c7 * B.c31 +
                      A.c8 * B.c23 - A.c9 * B.c22 - A.c10 * B.c4 + A.c11 * B.c20 -
                      A.c12 * B.c19 - A.c13 * B.c5 - A.c14 * B.c3 - A.c15 * B.c18;
    ctype const c12 = A.c0 * B.c5 - A.c1 * B.c21 + A.c2 * B.c20 - A.c4 * B.c1 -
                      A.c4 * B.c16 - A.c5 * B.c5 + A.c6 * B.c21 - A.c7 * B.c20 -
                      A.c8 * B.c25 - A.c8 * B.c31 + A.c9 * B.c3 + A.c9 * B.c18 -
                      A.c10 * B.c2 - A.c10 * B.c17 + A.c14 * B.c19 - A.c15 * B.c19;
    ctype const c13 = A.c0 * B.c21 + A.c1 * B.c5 - A.c2 * B.c19 - A.c4 * B.c2 -
                      A.c4 * B.c17 - A.c5 * B.c21 - A.c6 * B.c5 + A.c7 * B.c19 -
                      A.c8 * B.c3 - A.c8 * B.c18 - A.c9 * B.c25 - A.c9 * B.c31 +
                      A.c10 * B.c1 + A.c10 * B.c16 + A.c14 * B.c20 - A.c15 * B.c20;
    ctype const c14 = -A.c0 * B.c20 + A.c1 * B.c19 + A.c2 * B.c5 - A.c4 * B.c3 -
                      A.c4 * B.c18 + A.c5 * B.c20 - A.c6 * B.c19 - A.c7 * B.c5 +
                      A.c8 * B.c2 + A.c8 * B.c17 - A.c9 * B.c1 - A.c9 * B.c16 -
                      A.c10 * B.c25 - A.c10 * B.c31 + A.c14 * B.c21 - A.c15 * B.c21;
    ctype const c15 = -A.c0 * B.c16 - A.c1 * B.c17 - A.c2 * B.c18 + A.c3 * B.c5 -
                      A.c4 * B.c4 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 +
                      A.c8 * B.c22 + A.c9 * B.c23 + A.c10 * B.c24 - A.c11 * B.c19 -
                      A.c12 * B.c20 - A.c13 * B.c21 + A.c14 * B.c31 + A.c15 * B.c25;
    ctype const c16 = -A.c0 * B.c15 + A.c1 * B.c28 - A.c2 * B.c27 + A.c3 * B.c12 +
                      A.c4 * B.c6 + A.c5 * B.c0 - A.c6 * B.c11 + A.c7 * B.c10 +
                      A.c8 * B.c30 - A.c9 * B.c8 + A.c10 * B.c7 + A.c11 * B.c29 -
                      A.c12 * B.c14 + A.c13 * B.c13 - A.c14 * B.c26 + A.c15 * B.c9;
    ctype const c17 = -A.c0 * B.c28 - A.c1 * B.c15 + A.c2 * B.c26 + A.c3 * B.c13 +
                      A.c4 * B.c7 + A.c5 * B.c11 + A.c6 * B.c0 - A.c7 * B.c9 +
                      A.c8 * B.c8 + A.c9 * B.c30 - A.c10 * B.c6 + A.c11 * B.c14 +
                      A.c12 * B.c29 - A.c13 * B.c12 - A.c14 * B.c27 + A.c15 * B.c10;
    ctype const c18 = A.c0 * B.c27 - A.c1 * B.c26 - A.c2 * B.c15 + A.c3 * B.c14 +
                      A.c4 * B.c8 - A.c5 * B.c10 + A.c6 * B.c9 + A.c7 * B.c0 -
                      A.c8 * B.c7 + A.c9 * B.c6 + A.c10 * B.c30 - A.c11 * B.c13 +
                      A.c12 * B.c12 + A.c13 * B.c29 - A.c14 * B.c28 + A.c15 * B.c11;
    ctype const c19 = -A.c0 * B.c29 + A.c1 * B.c14 - A.c2 * B.c13 + A.c4 * B.c9 -
                      A.c4 * B.c26 + A.c5 * B.c29 - A.c6 * B.c14 + A.c7 * B.c13 +
                      A.c8 * B.c0 - A.c8 * B.c15 - A.c9 * B.c11 + A.c9 * B.c28 +
                      A.c10 * B.c10 - A.c10 * B.c27 - A.c14 * B.c12 + A.c15 * B.c12;
    ctype const c20 = -A.c0 * B.c14 - A.c1 * B.c29 + A.c2 * B.c12 + A.c4 * B.c10 -
                      A.c4 * B.c27 + A.c5 * B.c14 + A.c6 * B.c29 - A.c7 * B.c12 +
                      A.c8 * B.c11 - A.c8 * B.c28 + A.c9 * B.c0 - A.c9 * B.c15 -
                      A.c10 * B.c9 + A.c10 * B.c26 - A.c14 * B.c13 + A.c15 * B.c13;
    ctype const c21 = A.c0 * B.c13 - A.c1 * B.c12 - A.c2 * B.c29 + A.c4 * B.c11 -
                      A.c4 * B.c28 - A.c5 * B.c13 + A.c6 * B.c12 + A.c7 * B.c29 -
                      A.c8 * B.c10 + A.c8 * B.c27 + A.c9 * B.c9 - A.c9 * B.c26 +
                      A.c10 * B.c0 - A.c10 * B.c15 - A.c14 * B.c14 + A.c15 * B.c14;
    ctype const c22 = A.c0 * B.c30 - A.c1 * B.c8 + A.c2 * B.c7 + A.c3 * B.c9 +
                      A.c3 * B.c26 + A.c5 * B.c30 - A.c6 * B.c8 + A.c7 * B.c7 +
                      A.c11 * B.c0 + A.c11 * B.c15 - A.c12 * B.c11 - A.c12 * B.c28 +
                      A.c13 * B.c10 + A.c13 * B.c27 + A.c14 * B.c6 + A.c15 * B.c6;
    ctype const c23 = A.c0 * B.c8 + A.c1 * B.c30 - A.c2 * B.c6 + A.c3 * B.c10 +
                      A.c3 * B.c27 + A.c5 * B.c8 + A.c6 * B.c30 - A.c7 * B.c6 +
                      A.c11 * B.c11 + A.c11 * B.c28 + A.c12 * B.c0 + A.c12 * B.c15 -
                      A.c13 * B.c9 - A.c13 * B.c26 + A.c14 * B.c7 + A.c15 * B.c7;
    ctype const c24 = -A.c0 * B.c7 + A.c1 * B.c6 + A.c2 * B.c30 + A.c3 * B.c11 +
                      A.c3 * B.c28 - A.c5 * B.c7 + A.c6 * B.c6 + A.c7 * B.c30 -
                      A.c11 * B.c10 - A.c11 * B.c27 + A.c12 * B.c9 + A.c12 * B.c26 +
                      A.c13 * B.c0 + A.c13 * B.c15 + A.c14 * B.c8 + A.c15 * B.c8;
    ctype const c25 = -A.c0 * B.c9 - A.c1 * B.c10 - A.c2 * B.c11 + A.c3 * B.c29 -
                      A.c4 * B.c30 + A.c5 * B.c26 + A.c6 * B.c27 + A.c7 * B.c28 +
                      A.c8 * B.c6 + A.c9 * B.c7 + A.c10 * B.c8 - A.c11 * B.c12 -
                      A.c12 * B.c13 - A.c13 * B.c14 + A.c14 * B.c0 - A.c15 * B.c15;
    ctype const c26 = A.c0 * B.c31 - A.c1 * B.c18 + A.c2 * B.c17 + A.c3 * B.c19 -
                      A.c4 * B.c22 + A.c5 * B.c25 - A.c6 * B.c3 + A.c7 * B.c2 -
                      A.c8 * B.c4 + A.c9 * B.c24 - A.c10 * B.c23 + A.c11 * B.c5 -
                      A.c12 * B.c21 + A.c13 * B.c20 + A.c14 * B.c16 + A.c15 * B.c1;
    ctype const c27 = A.c0 * B.c18 + A.c1 * B.c31 - A.c2 * B.c16 + A.c3 * B.c20 -
                      A.c4 * B.c23 + A.c5 * B.c3 + A.c6 * B.c25 - A.c7 * B.c1 -
                      A.c8 * B.c24 - A.c9 * B.c4 + A.c10 * B.c22 + A.c11 * B.c21 +
                      A.c12 * B.c5 - A.c13 * B.c19 + A.c14 * B.c17 + A.c15 * B.c2;
    ctype const c28 = -A.c0 * B.c17 + A.c1 * B.c16 + A.c2 * B.c31 + A.c3 * B.c21 -
                      A.c4 * B.c24 - A.c5 * B.c2 + A.c6 * B.c1 + A.c7 * B.c25 +
                      A.c8 * B.c23 - A.c9 * B.c22 - A.c10 * B.c4 - A.c11 * B.c20 +
                      A.c12 * B.c19 + A.c13 * B.c5 + A.c14 * B.c18 + A.c15 * B.c3;
    ctype const c29 = -A.c0 * B.c19 - A.c1 * B.c20 - A.c2 * B.c21 - A.c4 * B.c25 -
                      A.c4 * B.c31 + A.c5 * B.c19 + A.c6 * B.c20 + A.c7 * B.c21 +
                      A.c8 * B.c1 + A.c8 * B.c16 + A.c9 * B.c2 + A.c9 * B.c17 +
                      A.c10 * B.c3 + A.c10 * B.c18 + A.c14 * B.c5 - A.c15 * B.c5;
    ctype const c30 = A.c0 * B.c22 + A.c1 * B.c23 + A.c2 * B.c24 + A.c3 * B.c25 -
                      A.c3 * B.c31 + A.c5 * B.c22 + A.c6 * B.c23 + A.c7 * B.c24 -
                      A.c11 * B.c1 + A.c11 * B.c16 - A.c12 * B.c2 + A.c12 * B.c17 -
                      A.c13 * B.c3 + A.c13 * B.c18 - A.c14 * B.c4 - A.c15 * B.c4;
    ctype const c31 = A.c0 * B.c26 + A.c1 * B.c27 + A.c2 * B.c28 + A.c3 * B.c29 +
                      A.c4 * B.c30 - A.c5 * B.c9 - A.c6 * B.c10 - A.c7 * B.c11 -
                      A.c8 * B.c6 - A.c9 * B.c7 - A.c10 * B.c8 - A.c11 * B.c12 -
                      A.c12 * B.c13 - A.c13 * B.c14 - A.c14 * B.c15 + A.c15 * B.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& M,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c31 * ctype(ps);
    ctype const c1 = -M.c26 * ctype(ps);
    ctype const c2 = -M.c27 * ctype(ps);
    ctype const c3 = -M.c28 * ctype(ps);
    ctype const c4 = M.c30 * ctype(ps);
    ctype const c5 = M.c29 * ctype(ps);
    ctype const c6 = -M.c22 * ctype(ps);
    ctype const c7 = -M.c23 * ctype(ps);
    ctype const c8 = -M.c24 * ctype(ps);
    ctype const c9 = -M.c16 * ctype(ps);
    ctype const c10 = -M.c17 * ctype(ps);
    ctype const c11 = -M.c18 * ctype(ps);
    ctype const c12 = -M.c19 * ctype(ps);
    ctype const c13 = -M.c20 * ctype(ps);
    ctype const c14 = -M.c21 * ctype(ps);
    ctype const c15 = M.c25 * ctype(ps);
    ctype const c16 = M.c9 * ctype(ps);
    ctype const c17 = M.c10 * ctype(ps);
    ctype const c18 = M.c11 * ctype(ps);
    ctype const c19 = M.c12 * ctype(ps);
    ctype const c20 = M.c13 * ctype(ps);
    ctype const c21 = M.c14 * ctype(ps);
    ctype const c22 = M.c6 * ctype(ps);
    ctype const c23 = M.c7 * ctype(ps);
    ctype const c24 = M.c8 * ctype(ps);
    ctype const c25 = -M.c15 * ctype(ps);
    ctype const c26 = M.c1 * ctype(ps);
    ctype const c27 = M.c2 * ctype(ps);
    ctype const c28 = M.c3 * ctype(ps);
    ctype const c29 = -M.c5 * ctype(ps);
    ctype const c30 = -M.c4 * ctype(ps);
    ctype const c31 = M.c0 * ctype(ps);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c31;
    ctype const c1 = -ctype(ps) * M.c26;
    ctype const c2 = -ctype(ps) * M.c27;
    ctype const c3 = -ctype(ps) * M.c28;
    ctype const c4 = ctype(ps) * M.c30;
    ctype const c5 = ctype(ps) * M.c29;
    ctype const c6 = -ctype(ps) * M.c22;
    ctype const c7 = -ctype(ps) * M.c23;
    ctype const c8 = -ctype(ps) * M.c24;
    ctype const c9 = -ctype(ps) * M.c16;
    ctype const c10 = -ctype(ps) * M.c17;
    ctype const c11 = -ctype(ps) * M.c18;
    ctype const c12 = -ctype(ps) * M.c19;
    ctype const c13 = -ctype(ps) * M.c20;
    ctype const c14 = -ctype(ps) * M.c21;
    ctype const c15 = ctype(ps) * M.c25;
    ctype const c16 = ctype(ps) * M.c9;
    ctype const c17 = ctype(ps) * M.c10;
    ctype const c18 = ctype(ps) * M.c11;
    ctype const c19 = ctype(ps) * M.c12;
    ctype const c20 = ctype(ps) * M.c13;
    ctype const c21 = ctype(ps) * M.c14;
    ctype const c22 = ctype(ps) * M.c6;
    ctype const c23 = ctype(ps) * M.c7;
    ctype const c24 = ctype(ps) * M.c8;
    ctype const c25 = -ctype(ps) * M.c15;
    ctype const c26 = ctype(ps) * M.c1;
    ctype const c27 = ctype(ps) * M.c2;
    ctype const c28 = ctype(ps) * M.c3;
    ctype const c29 = -ctype(ps) * M.c5;
    ctype const c30 = -ctype(ps) * M.c4;
    ctype const c31 = ctype(ps) * M.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,quadvec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& M,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c26 * Q.x - M.c27 * Q.y - M.c28 * Q.z + M.c29 * Q.u + M.c30 * Q.w;
    ctype const c1 = M.c17 * Q.z - M.c18 * Q.y - M.c19 * Q.u + M.c22 * Q.w - M.c31 * Q.x;
    ctype const c2 = -M.c16 * Q.z + M.c18 * Q.x - M.c20 * Q.u + M.c23 * Q.w - M.c31 * Q.y;
    ctype const c3 = M.c16 * Q.y - M.c17 * Q.x - M.c21 * Q.u + M.c24 * Q.w - M.c31 * Q.z;
    ctype const c4 = -M.c22 * Q.x - M.c23 * Q.y - M.c24 * Q.z + M.c25 * Q.u + M.c31 * Q.u;
    ctype const c5 = M.c19 * Q.x + M.c20 * Q.y + M.c21 * Q.z - M.c25 * Q.w + M.c31 * Q.w;
    ctype const c6 = -M.c7 * Q.z + M.c8 * Q.y + M.c9 * Q.u - M.c26 * Q.u + M.c30 * Q.x;
    ctype const c7 = M.c6 * Q.z - M.c8 * Q.x + M.c10 * Q.u - M.c27 * Q.u + M.c30 * Q.y;
    ctype const c8 = -M.c6 * Q.y + M.c7 * Q.x + M.c11 * Q.u - M.c28 * Q.u + M.c30 * Q.z;
    ctype const c9 = M.c6 * Q.w + M.c12 * Q.u + M.c15 * Q.x - M.c27 * Q.z + M.c28 * Q.y;
    ctype const c10 = M.c7 * Q.w + M.c13 * Q.u + M.c15 * Q.y + M.c26 * Q.z - M.c28 * Q.x;
    ctype const c11 = M.c8 * Q.w + M.c14 * Q.u + M.c15 * Q.z - M.c26 * Q.y + M.c27 * Q.x;
    ctype const c12 = M.c9 * Q.w + M.c13 * Q.z - M.c14 * Q.y + M.c26 * Q.w - M.c29 * Q.x;
    ctype const c13 = M.c10 * Q.w - M.c12 * Q.z + M.c14 * Q.x + M.c27 * Q.w - M.c29 * Q.y;
    ctype const c14 = M.c11 * Q.w + M.c12 * Q.y - M.c13 * Q.x + M.c28 * Q.w - M.c29 * Q.z;
    ctype const c15 = -M.c9 * Q.x - M.c10 * Q.y - M.c11 * Q.z + M.c29 * Q.u - M.c30 * Q.w;
    ctype const c16 = M.c2 * Q.z - M.c3 * Q.y + M.c19 * Q.u + M.c22 * Q.w - M.c25 * Q.x;
    ctype const c17 = -M.c1 * Q.z + M.c3 * Q.x + M.c20 * Q.u + M.c23 * Q.w - M.c25 * Q.y;
    ctype const c18 = M.c1 * Q.y - M.c2 * Q.x + M.c21 * Q.u + M.c24 * Q.w - M.c25 * Q.z;
    ctype const c19 = -M.c1 * Q.w - M.c5 * Q.x + M.c16 * Q.w + M.c20 * Q.z - M.c21 * Q.y;
    ctype const c20 = -M.c2 * Q.w - M.c5 * Q.y + M.c17 * Q.w - M.c19 * Q.z + M.c21 * Q.x;
    ctype const c21 = -M.c3 * Q.w - M.c5 * Q.z + M.c18 * Q.w + M.c19 * Q.y - M.c20 * Q.x;
    ctype const c22 = M.c1 * Q.u + M.c4 * Q.x + M.c16 * Q.u - M.c23 * Q.z + M.c24 * Q.y;
    ctype const c23 = M.c2 * Q.u + M.c4 * Q.y + M.c17 * Q.u + M.c22 * Q.z - M.c24 * Q.x;
    ctype const c24 = M.c3 * Q.u + M.c4 * Q.z + M.c18 * Q.u - M.c22 * Q.y + M.c23 * Q.x;
    ctype const c25 = M.c4 * Q.w - M.c5 * Q.u + M.c16 * Q.x + M.c17 * Q.y + M.c18 * Q.z;
    ctype const c26 = M.c0 * Q.x - M.c6 * Q.w - M.c10 * Q.z + M.c11 * Q.y + M.c12 * Q.u;
    ctype const c27 = M.c0 * Q.y - M.c7 * Q.w + M.c9 * Q.z - M.c11 * Q.x + M.c13 * Q.u;
    ctype const c28 = M.c0 * Q.z - M.c8 * Q.w - M.c9 * Q.y + M.c10 * Q.x + M.c14 * Q.u;
    ctype const c29 = M.c0 * Q.w + M.c12 * Q.x + M.c13 * Q.y + M.c14 * Q.z + M.c15 * Q.w;
    ctype const c30 = M.c0 * Q.u - M.c6 * Q.x - M.c7 * Q.y - M.c8 * Q.z - M.c15 * Q.u;
    ctype const c31 = M.c1 * Q.x + M.c2 * Q.y + M.c3 * Q.z + M.c4 * Q.w + M.c5 * Q.u;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(quadvec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c26 - Q.y * M.c27 - Q.z * M.c28 + Q.w * M.c30 + Q.u * M.c29;
    ctype const c1 = -Q.x * M.c31 + Q.y * M.c18 - Q.z * M.c17 - Q.w * M.c22 + Q.u * M.c19;
    ctype const c2 = -Q.x * M.c18 - Q.y * M.c31 + Q.z * M.c16 - Q.w * M.c23 + Q.u * M.c20;
    ctype const c3 = Q.x * M.c17 - Q.y * M.c16 - Q.z * M.c31 - Q.w * M.c24 + Q.u * M.c21;
    ctype const c4 = Q.x * M.c22 + Q.y * M.c23 + Q.z * M.c24 - Q.u * M.c25 + Q.u * M.c31;
    ctype const c5 = -Q.x * M.c19 - Q.y * M.c20 - Q.z * M.c21 + Q.w * M.c25 + Q.w * M.c31;
    ctype const c6 = -Q.x * M.c30 + Q.y * M.c8 - Q.z * M.c7 + Q.u * M.c9 + Q.u * M.c26;
    ctype const c7 = -Q.x * M.c8 - Q.y * M.c30 + Q.z * M.c6 + Q.u * M.c10 + Q.u * M.c27;
    ctype const c8 = Q.x * M.c7 - Q.y * M.c6 - Q.z * M.c30 + Q.u * M.c11 + Q.u * M.c28;
    ctype const c9 = Q.x * M.c15 - Q.y * M.c28 + Q.z * M.c27 + Q.w * M.c6 + Q.u * M.c12;
    ctype const c10 = Q.x * M.c28 + Q.y * M.c15 - Q.z * M.c26 + Q.w * M.c7 + Q.u * M.c13;
    ctype const c11 = -Q.x * M.c27 + Q.y * M.c26 + Q.z * M.c15 + Q.w * M.c8 + Q.u * M.c14;
    ctype const c12 = Q.x * M.c29 - Q.y * M.c14 + Q.z * M.c13 + Q.w * M.c9 - Q.w * M.c26;
    ctype const c13 = Q.x * M.c14 + Q.y * M.c29 - Q.z * M.c12 + Q.w * M.c10 - Q.w * M.c27;
    ctype const c14 =
        -Q.x * M.c13 + Q.y * M.c12 + Q.z * M.c29 + Q.w * M.c11 - Q.w * M.c28;
    ctype const c15 = -Q.x * M.c9 - Q.y * M.c10 - Q.z * M.c11 + Q.w * M.c30 - Q.u * M.c29;
    ctype const c16 = -Q.x * M.c25 + Q.y * M.c3 - Q.z * M.c2 + Q.w * M.c22 + Q.u * M.c19;
    ctype const c17 = -Q.x * M.c3 - Q.y * M.c25 + Q.z * M.c1 + Q.w * M.c23 + Q.u * M.c20;
    ctype const c18 = Q.x * M.c2 - Q.y * M.c1 - Q.z * M.c25 + Q.w * M.c24 + Q.u * M.c21;
    ctype const c19 = Q.x * M.c5 - Q.y * M.c21 + Q.z * M.c20 + Q.w * M.c1 + Q.w * M.c16;
    ctype const c20 = Q.x * M.c21 + Q.y * M.c5 - Q.z * M.c19 + Q.w * M.c2 + Q.w * M.c17;
    ctype const c21 = -Q.x * M.c20 + Q.y * M.c19 + Q.z * M.c5 + Q.w * M.c3 + Q.w * M.c18;
    ctype const c22 = -Q.x * M.c4 + Q.y * M.c24 - Q.z * M.c23 - Q.u * M.c1 + Q.u * M.c16;
    ctype const c23 = -Q.x * M.c24 - Q.y * M.c4 + Q.z * M.c22 - Q.u * M.c2 + Q.u * M.c17;
    ctype const c24 = Q.x * M.c23 - Q.y * M.c22 - Q.z * M.c4 - Q.u * M.c3 + Q.u * M.c18;
    ctype const c25 = Q.x * M.c16 + Q.y * M.c17 + Q.z * M.c18 - Q.w * M.c4 + Q.u * M.c5;
    ctype const c26 = Q.x * M.c0 - Q.y * M.c11 + Q.z * M.c10 + Q.w * M.c6 - Q.u * M.c12;
    ctype const c27 = Q.x * M.c11 + Q.y * M.c0 - Q.z * M.c9 + Q.w * M.c7 - Q.u * M.c13;
    ctype const c28 = -Q.x * M.c10 + Q.y * M.c9 + Q.z * M.c0 + Q.w * M.c8 - Q.u * M.c14;
    ctype const c29 = -Q.x * M.c12 - Q.y * M.c13 - Q.z * M.c14 + Q.w * M.c0 - Q.w * M.c15;
    ctype const c30 = Q.x * M.c6 + Q.y * M.c7 + Q.z * M.c8 + Q.u * M.c0 + Q.u * M.c15;
    ctype const c31 = Q.x * M.c1 + Q.y * M.c2 + Q.z * M.c3 + Q.w * M.c4 + Q.u * M.c5;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& M,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c16 * t.vx + M.c17 * t.vy + M.c18 * t.vz + M.c19 * t.px +
                     M.c20 * t.py + M.c21 * t.pz + M.c22 * t.mx + M.c23 * t.my +
                     M.c24 * t.mz - M.c25 * t.pw;
    ctype const c1 = -M.c7 * t.mz + M.c8 * t.my + M.c9 * t.pw + M.c13 * t.pz -
                     M.c14 * t.py - M.c15 * t.vx + M.c27 * t.vz - M.c28 * t.vy -
                     M.c29 * t.px + M.c30 * t.mx;
    ctype const c2 = M.c6 * t.mz - M.c8 * t.mx + M.c10 * t.pw - M.c12 * t.pz +
                     M.c14 * t.px - M.c15 * t.vy - M.c26 * t.vz + M.c28 * t.vx -
                     M.c29 * t.py + M.c30 * t.my;
    ctype const c3 = -M.c6 * t.my + M.c7 * t.mx + M.c11 * t.pw + M.c12 * t.py -
                     M.c13 * t.px - M.c15 * t.vz + M.c26 * t.vy - M.c27 * t.vx -
                     M.c29 * t.pz + M.c30 * t.mz;
    ctype const c4 = -M.c6 * t.vx - M.c7 * t.vy - M.c8 * t.vz - M.c9 * t.px -
                     M.c10 * t.py - M.c11 * t.pz + M.c26 * t.px + M.c27 * t.py +
                     M.c28 * t.pz - M.c30 * t.pw;
    ctype const c5 = -M.c9 * t.mx - M.c10 * t.my - M.c11 * t.mz - M.c12 * t.vx -
                     M.c13 * t.vy - M.c14 * t.vz - M.c26 * t.mx - M.c27 * t.my -
                     M.c28 * t.mz + M.c29 * t.pw;
    ctype const c6 = M.c2 * t.pz - M.c3 * t.py - M.c4 * t.vx + M.c17 * t.pz -
                     M.c18 * t.py + M.c22 * t.pw + M.c23 * t.vz - M.c24 * t.vy -
                     M.c25 * t.px - M.c31 * t.px;
    ctype const c7 = -M.c1 * t.pz + M.c3 * t.px - M.c4 * t.vy - M.c16 * t.pz +
                     M.c18 * t.px - M.c22 * t.vz + M.c23 * t.pw + M.c24 * t.vx -
                     M.c25 * t.py - M.c31 * t.py;
    ctype const c8 = M.c1 * t.py - M.c2 * t.px - M.c4 * t.vz + M.c16 * t.py -
                     M.c17 * t.px + M.c22 * t.vy - M.c23 * t.vx + M.c24 * t.pw -
                     M.c25 * t.pz - M.c31 * t.pz;
    ctype const c9 = -M.c1 * t.pw - M.c4 * t.mx - M.c5 * t.px + M.c17 * t.vz -
                     M.c18 * t.vy + M.c20 * t.pz - M.c21 * t.py + M.c23 * t.mz -
                     M.c24 * t.my - M.c31 * t.vx;
    ctype const c10 = -M.c2 * t.pw - M.c4 * t.my - M.c5 * t.py - M.c16 * t.vz +
                      M.c18 * t.vx - M.c19 * t.pz + M.c21 * t.px - M.c22 * t.mz +
                      M.c24 * t.mx - M.c31 * t.vy;
    ctype const c11 = -M.c3 * t.pw - M.c4 * t.mz - M.c5 * t.pz + M.c16 * t.vy -
                      M.c17 * t.vx + M.c19 * t.py - M.c20 * t.px + M.c22 * t.my -
                      M.c23 * t.mx - M.c31 * t.vz;
    ctype const c12 = -M.c2 * t.mz + M.c3 * t.my - M.c5 * t.vx + M.c17 * t.mz -
                      M.c18 * t.my - M.c19 * t.pw + M.c20 * t.vz - M.c21 * t.vy +
                      M.c25 * t.mx - M.c31 * t.mx;
    ctype const c13 = M.c1 * t.mz - M.c3 * t.mx - M.c5 * t.vy - M.c16 * t.mz +
                      M.c18 * t.mx - M.c19 * t.vz - M.c20 * t.pw + M.c21 * t.vx +
                      M.c25 * t.my - M.c31 * t.my;
    ctype const c14 = -M.c1 * t.my + M.c2 * t.mx - M.c5 * t.vz + M.c16 * t.my -
                      M.c17 * t.mx + M.c19 * t.vy - M.c20 * t.vx - M.c21 * t.pw +
                      M.c25 * t.mz - M.c31 * t.mz;
    ctype const c15 = -M.c1 * t.vx - M.c2 * t.vy - M.c3 * t.vz + M.c19 * t.px +
                      M.c20 * t.py + M.c21 * t.pz - M.c22 * t.mx - M.c23 * t.my -
                      M.c24 * t.mz + M.c31 * t.pw;
    ctype const c16 = M.c0 * t.vx - M.c7 * t.mz + M.c8 * t.my - M.c10 * t.vz +
                      M.c11 * t.vy - M.c13 * t.pz + M.c14 * t.py - M.c26 * t.pw +
                      M.c29 * t.px + M.c30 * t.mx;
    ctype const c17 = M.c0 * t.vy + M.c6 * t.mz - M.c8 * t.mx + M.c9 * t.vz -
                      M.c11 * t.vx + M.c12 * t.pz - M.c14 * t.px - M.c27 * t.pw +
                      M.c29 * t.py + M.c30 * t.my;
    ctype const c18 = M.c0 * t.vz - M.c6 * t.my + M.c7 * t.mx - M.c9 * t.vy +
                      M.c10 * t.vx - M.c12 * t.py + M.c13 * t.px - M.c28 * t.pw +
                      M.c29 * t.pz + M.c30 * t.mz;
    ctype const c19 = M.c0 * t.mx - M.c10 * t.mz + M.c11 * t.my + M.c12 * t.pw -
                      M.c13 * t.vz + M.c14 * t.vy + M.c15 * t.mx - M.c27 * t.mz +
                      M.c28 * t.my + M.c29 * t.vx;
    ctype const c20 = M.c0 * t.my + M.c9 * t.mz - M.c11 * t.mx + M.c12 * t.vz +
                      M.c13 * t.pw - M.c14 * t.vx + M.c15 * t.my + M.c26 * t.mz -
                      M.c28 * t.mx + M.c29 * t.vy;
    ctype const c21 = M.c0 * t.mz - M.c9 * t.my + M.c10 * t.mx - M.c12 * t.vy +
                      M.c13 * t.vx + M.c14 * t.pw + M.c15 * t.mz - M.c26 * t.my +
                      M.c27 * t.mx + M.c29 * t.vz;
    ctype const c22 = M.c0 * t.px - M.c6 * t.pw - M.c7 * t.vz + M.c8 * t.vy -
                      M.c10 * t.pz + M.c11 * t.py - M.c15 * t.px + M.c27 * t.pz -
                      M.c28 * t.py + M.c30 * t.vx;
    ctype const c23 = M.c0 * t.py + M.c6 * t.vz - M.c7 * t.pw - M.c8 * t.vx +
                      M.c9 * t.pz - M.c11 * t.px - M.c15 * t.py - M.c26 * t.pz +
                      M.c28 * t.px + M.c30 * t.vy;
    ctype const c24 = M.c0 * t.pz - M.c6 * t.vy + M.c7 * t.vx - M.c8 * t.pw -
                      M.c9 * t.py + M.c10 * t.px - M.c15 * t.pz + M.c26 * t.py -
                      M.c27 * t.px + M.c30 * t.vz;
    ctype const c25 = M.c0 * t.pw - M.c6 * t.mx - M.c7 * t.my - M.c8 * t.mz +
                      M.c12 * t.px + M.c13 * t.py + M.c14 * t.pz + M.c26 * t.vx +
                      M.c27 * t.vy + M.c28 * t.vz;
    ctype const c26 = -M.c2 * t.vz + M.c3 * t.vy + M.c4 * t.mx - M.c5 * t.px +
                      M.c16 * t.pw + M.c20 * t.pz - M.c21 * t.py - M.c23 * t.mz +
                      M.c24 * t.my + M.c25 * t.vx;
    ctype const c27 = M.c1 * t.vz - M.c3 * t.vx + M.c4 * t.my - M.c5 * t.py +
                      M.c17 * t.pw - M.c19 * t.pz + M.c21 * t.px + M.c22 * t.mz -
                      M.c24 * t.mx + M.c25 * t.vy;
    ctype const c28 = -M.c1 * t.vy + M.c2 * t.vx + M.c4 * t.mz - M.c5 * t.pz +
                      M.c18 * t.pw + M.c19 * t.py - M.c20 * t.px - M.c22 * t.my +
                      M.c23 * t.mx + M.c25 * t.vz;
    ctype const c29 = -M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz - M.c5 * t.pw +
                      M.c16 * t.mx + M.c17 * t.my + M.c18 * t.mz + M.c19 * t.vx +
                      M.c20 * t.vy + M.c21 * t.vz;
    ctype const c30 = M.c1 * t.px + M.c2 * t.py + M.c3 * t.pz + M.c4 * t.pw +
                      M.c16 * t.px + M.c17 * t.py + M.c18 * t.pz + M.c22 * t.vx +
                      M.c23 * t.vy + M.c24 * t.vz;
    ctype const c31 = -M.c6 * t.mx - M.c7 * t.my - M.c8 * t.mz - M.c9 * t.vx -
                      M.c10 * t.vy - M.c11 * t.vz - M.c12 * t.px - M.c13 * t.py -
                      M.c14 * t.pz - M.c15 * t.pw;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * M.c16 + t.vy * M.c17 + t.vz * M.c18 + t.mx * M.c22 +
                     t.my * M.c23 + t.mz * M.c24 + t.px * M.c19 + t.py * M.c20 +
                     t.pz * M.c21 - t.pw * M.c25;
    ctype const c1 = -t.vx * M.c15 + t.vy * M.c28 - t.vz * M.c27 - t.mx * M.c30 +
                     t.my * M.c8 - t.mz * M.c7 + t.px * M.c29 - t.py * M.c14 +
                     t.pz * M.c13 + t.pw * M.c9;
    ctype const c2 = -t.vx * M.c28 - t.vy * M.c15 + t.vz * M.c26 - t.mx * M.c8 -
                     t.my * M.c30 + t.mz * M.c6 + t.px * M.c14 + t.py * M.c29 -
                     t.pz * M.c12 + t.pw * M.c10;
    ctype const c3 = t.vx * M.c27 - t.vy * M.c26 - t.vz * M.c15 + t.mx * M.c7 -
                     t.my * M.c6 - t.mz * M.c30 - t.px * M.c13 + t.py * M.c12 +
                     t.pz * M.c29 + t.pw * M.c11;
    ctype const c4 = -t.vx * M.c6 - t.vy * M.c7 - t.vz * M.c8 - t.px * M.c9 -
                     t.px * M.c26 - t.py * M.c10 - t.py * M.c27 - t.pz * M.c11 -
                     t.pz * M.c28 + t.pw * M.c30;
    ctype const c5 = -t.vx * M.c12 - t.vy * M.c13 - t.vz * M.c14 - t.mx * M.c9 +
                     t.mx * M.c26 - t.my * M.c10 + t.my * M.c27 - t.mz * M.c11 +
                     t.mz * M.c28 - t.pw * M.c29;
    ctype const c6 = -t.vx * M.c4 + t.vy * M.c24 - t.vz * M.c23 + t.px * M.c25 -
                     t.px * M.c31 - t.py * M.c3 + t.py * M.c18 + t.pz * M.c2 -
                     t.pz * M.c17 - t.pw * M.c22;
    ctype const c7 = -t.vx * M.c24 - t.vy * M.c4 + t.vz * M.c22 + t.px * M.c3 -
                     t.px * M.c18 + t.py * M.c25 - t.py * M.c31 - t.pz * M.c1 +
                     t.pz * M.c16 - t.pw * M.c23;
    ctype const c8 = t.vx * M.c23 - t.vy * M.c22 - t.vz * M.c4 - t.px * M.c2 +
                     t.px * M.c17 + t.py * M.c1 - t.py * M.c16 + t.pz * M.c25 -
                     t.pz * M.c31 - t.pw * M.c24;
    ctype const c9 = -t.vx * M.c31 + t.vy * M.c18 - t.vz * M.c17 - t.mx * M.c4 +
                     t.my * M.c24 - t.mz * M.c23 - t.px * M.c5 + t.py * M.c21 -
                     t.pz * M.c20 - t.pw * M.c1;
    ctype const c10 = -t.vx * M.c18 - t.vy * M.c31 + t.vz * M.c16 - t.mx * M.c24 -
                      t.my * M.c4 + t.mz * M.c22 - t.px * M.c21 - t.py * M.c5 +
                      t.pz * M.c19 - t.pw * M.c2;
    ctype const c11 = t.vx * M.c17 - t.vy * M.c16 - t.vz * M.c31 + t.mx * M.c23 -
                      t.my * M.c22 - t.mz * M.c4 + t.px * M.c20 - t.py * M.c19 -
                      t.pz * M.c5 - t.pw * M.c3;
    ctype const c12 = -t.vx * M.c5 + t.vy * M.c21 - t.vz * M.c20 - t.mx * M.c25 -
                      t.mx * M.c31 + t.my * M.c3 + t.my * M.c18 - t.mz * M.c2 -
                      t.mz * M.c17 + t.pw * M.c19;
    ctype const c13 = -t.vx * M.c21 - t.vy * M.c5 + t.vz * M.c19 - t.mx * M.c3 -
                      t.mx * M.c18 - t.my * M.c25 - t.my * M.c31 + t.mz * M.c1 +
                      t.mz * M.c16 + t.pw * M.c20;
    ctype const c14 = t.vx * M.c20 - t.vy * M.c19 - t.vz * M.c5 + t.mx * M.c2 +
                      t.mx * M.c17 - t.my * M.c1 - t.my * M.c16 - t.mz * M.c25 -
                      t.mz * M.c31 + t.pw * M.c21;
    ctype const c15 = -t.vx * M.c1 - t.vy * M.c2 - t.vz * M.c3 + t.mx * M.c22 +
                      t.my * M.c23 + t.mz * M.c24 - t.px * M.c19 - t.py * M.c20 -
                      t.pz * M.c21 + t.pw * M.c31;
    ctype const c16 = t.vx * M.c0 - t.vy * M.c11 + t.vz * M.c10 + t.mx * M.c30 -
                      t.my * M.c8 + t.mz * M.c7 + t.px * M.c29 - t.py * M.c14 +
                      t.pz * M.c13 - t.pw * M.c26;
    ctype const c17 = t.vx * M.c11 + t.vy * M.c0 - t.vz * M.c9 + t.mx * M.c8 +
                      t.my * M.c30 - t.mz * M.c6 + t.px * M.c14 + t.py * M.c29 -
                      t.pz * M.c12 - t.pw * M.c27;
    ctype const c18 = -t.vx * M.c10 + t.vy * M.c9 + t.vz * M.c0 - t.mx * M.c7 +
                      t.my * M.c6 + t.mz * M.c30 - t.px * M.c13 + t.py * M.c12 +
                      t.pz * M.c29 - t.pw * M.c28;
    ctype const c19 = t.vx * M.c29 - t.vy * M.c14 + t.vz * M.c13 + t.mx * M.c0 -
                      t.mx * M.c15 - t.my * M.c11 + t.my * M.c28 + t.mz * M.c10 -
                      t.mz * M.c27 - t.pw * M.c12;
    ctype const c20 = t.vx * M.c14 + t.vy * M.c29 - t.vz * M.c12 + t.mx * M.c11 -
                      t.mx * M.c28 + t.my * M.c0 - t.my * M.c15 - t.mz * M.c9 +
                      t.mz * M.c26 - t.pw * M.c13;
    ctype const c21 = -t.vx * M.c13 + t.vy * M.c12 + t.vz * M.c29 - t.mx * M.c10 +
                      t.mx * M.c27 + t.my * M.c9 - t.my * M.c26 + t.mz * M.c0 -
                      t.mz * M.c15 - t.pw * M.c14;
    ctype const c22 = t.vx * M.c30 - t.vy * M.c8 + t.vz * M.c7 + t.px * M.c0 +
                      t.px * M.c15 - t.py * M.c11 - t.py * M.c28 + t.pz * M.c10 +
                      t.pz * M.c27 + t.pw * M.c6;
    ctype const c23 = t.vx * M.c8 + t.vy * M.c30 - t.vz * M.c6 + t.px * M.c11 +
                      t.px * M.c28 + t.py * M.c0 + t.py * M.c15 - t.pz * M.c9 -
                      t.pz * M.c26 + t.pw * M.c7;
    ctype const c24 = -t.vx * M.c7 + t.vy * M.c6 + t.vz * M.c30 - t.px * M.c10 -
                      t.px * M.c27 + t.py * M.c9 + t.py * M.c26 + t.pz * M.c0 +
                      t.pz * M.c15 + t.pw * M.c8;
    ctype const c25 = t.vx * M.c26 + t.vy * M.c27 + t.vz * M.c28 + t.mx * M.c6 +
                      t.my * M.c7 + t.mz * M.c8 - t.px * M.c12 - t.py * M.c13 -
                      t.pz * M.c14 + t.pw * M.c0;
    ctype const c26 = t.vx * M.c25 - t.vy * M.c3 + t.vz * M.c2 - t.mx * M.c4 +
                      t.my * M.c24 - t.mz * M.c23 + t.px * M.c5 - t.py * M.c21 +
                      t.pz * M.c20 + t.pw * M.c16;
    ctype const c27 = t.vx * M.c3 + t.vy * M.c25 - t.vz * M.c1 - t.mx * M.c24 -
                      t.my * M.c4 + t.mz * M.c22 + t.px * M.c21 + t.py * M.c5 -
                      t.pz * M.c19 + t.pw * M.c17;
    ctype const c28 = -t.vx * M.c2 + t.vy * M.c1 + t.vz * M.c25 + t.mx * M.c23 -
                      t.my * M.c22 - t.mz * M.c4 - t.px * M.c20 + t.py * M.c19 +
                      t.pz * M.c5 + t.pw * M.c18;
    ctype const c29 = t.vx * M.c19 + t.vy * M.c20 + t.vz * M.c21 + t.mx * M.c1 +
                      t.mx * M.c16 + t.my * M.c2 + t.my * M.c17 + t.mz * M.c3 +
                      t.mz * M.c18 + t.pw * M.c5;
    ctype const c30 = t.vx * M.c22 + t.vy * M.c23 + t.vz * M.c24 - t.px * M.c1 +
                      t.px * M.c16 - t.py * M.c2 + t.py * M.c17 - t.pz * M.c3 +
                      t.pz * M.c18 - t.pw * M.c4;
    ctype const c31 = -t.vx * M.c9 - t.vy * M.c10 - t.vz * M.c11 - t.mx * M.c6 -
                      t.my * M.c7 - t.mz * M.c8 - t.px * M.c12 - t.py * M.c13 -
                      t.pz * M.c14 - t.pw * M.c15;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& M,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c6 * B.px - M.c7 * B.py - M.c8 * B.pz - M.c9 * B.mx -
                     M.c10 * B.my - M.c11 * B.mz - M.c12 * B.vx - M.c13 * B.vy -
                     M.c14 * B.vz + M.c15 * B.pw;
    ctype const c1 = -M.c2 * B.mz + M.c3 * B.my + M.c4 * B.px - M.c5 * B.vx -
                     M.c16 * B.pw + M.c20 * B.vz - M.c21 * B.vy - M.c23 * B.pz +
                     M.c24 * B.py + M.c25 * B.mx;
    ctype const c2 = M.c1 * B.mz - M.c3 * B.mx + M.c4 * B.py - M.c5 * B.vy -
                     M.c17 * B.pw - M.c19 * B.vz + M.c21 * B.vx + M.c22 * B.pz -
                     M.c24 * B.px + M.c25 * B.my;
    ctype const c3 = -M.c1 * B.my + M.c2 * B.mx + M.c4 * B.pz - M.c5 * B.vz -
                     M.c18 * B.pw + M.c19 * B.vy - M.c20 * B.vx - M.c22 * B.py +
                     M.c23 * B.px + M.c25 * B.mz;
    ctype const c4 = -M.c1 * B.vx - M.c2 * B.vy - M.c3 * B.vz + M.c4 * B.pw -
                     M.c16 * B.vx - M.c17 * B.vy - M.c18 * B.vz - M.c22 * B.mx -
                     M.c23 * B.my - M.c24 * B.mz;
    ctype const c5 = M.c1 * B.px + M.c2 * B.py + M.c3 * B.pz - M.c5 * B.pw -
                     M.c16 * B.px - M.c17 * B.py - M.c18 * B.pz - M.c19 * B.mx -
                     M.c20 * B.my - M.c21 * B.mz;
    ctype const c6 = M.c0 * B.vx + M.c6 * B.pw - M.c7 * B.mz + M.c8 * B.my -
                     M.c10 * B.vz + M.c11 * B.vy - M.c15 * B.vx + M.c27 * B.vz -
                     M.c28 * B.vy + M.c30 * B.mx;
    ctype const c7 = M.c0 * B.vy + M.c6 * B.mz + M.c7 * B.pw - M.c8 * B.mx + M.c9 * B.vz -
                     M.c11 * B.vx - M.c15 * B.vy - M.c26 * B.vz + M.c28 * B.vx +
                     M.c30 * B.my;
    ctype const c8 = M.c0 * B.vz - M.c6 * B.my + M.c7 * B.mx + M.c8 * B.pw - M.c9 * B.vy +
                     M.c10 * B.vx - M.c15 * B.vz + M.c26 * B.vy - M.c27 * B.vx +
                     M.c30 * B.mz;
    ctype const c9 = M.c0 * B.mx - M.c7 * B.pz + M.c8 * B.py - M.c10 * B.mz +
                     M.c11 * B.my - M.c13 * B.vz + M.c14 * B.vy + M.c26 * B.pw +
                     M.c29 * B.vx + M.c30 * B.px;
    ctype const c10 = M.c0 * B.my + M.c6 * B.pz - M.c8 * B.px + M.c9 * B.mz -
                      M.c11 * B.mx + M.c12 * B.vz - M.c14 * B.vx + M.c27 * B.pw +
                      M.c29 * B.vy + M.c30 * B.py;
    ctype const c11 = M.c0 * B.mz - M.c6 * B.py + M.c7 * B.px - M.c9 * B.my +
                      M.c10 * B.mx - M.c12 * B.vy + M.c13 * B.vx + M.c28 * B.pw +
                      M.c29 * B.vz + M.c30 * B.pz;
    ctype const c12 = M.c0 * B.px - M.c10 * B.pz + M.c11 * B.py - M.c12 * B.pw -
                      M.c13 * B.mz + M.c14 * B.my + M.c15 * B.px - M.c27 * B.pz +
                      M.c28 * B.py + M.c29 * B.mx;
    ctype const c13 = M.c0 * B.py + M.c9 * B.pz - M.c11 * B.px + M.c12 * B.mz -
                      M.c13 * B.pw - M.c14 * B.mx + M.c15 * B.py + M.c26 * B.pz -
                      M.c28 * B.px + M.c29 * B.my;
    ctype const c14 = M.c0 * B.pz - M.c9 * B.py + M.c10 * B.px - M.c12 * B.my +
                      M.c13 * B.mx - M.c14 * B.pw + M.c15 * B.pz - M.c26 * B.py +
                      M.c27 * B.px + M.c29 * B.mz;
    ctype const c15 = M.c0 * B.pw + M.c6 * B.px + M.c7 * B.py + M.c8 * B.pz -
                      M.c12 * B.vx - M.c13 * B.vy - M.c14 * B.vz - M.c26 * B.mx -
                      M.c27 * B.my - M.c28 * B.mz;
    ctype const c16 = -M.c1 * B.pw + M.c4 * B.px + M.c5 * B.vx - M.c17 * B.mz +
                      M.c18 * B.my - M.c20 * B.vz + M.c21 * B.vy - M.c23 * B.pz +
                      M.c24 * B.py + M.c31 * B.mx;
    ctype const c17 = -M.c2 * B.pw + M.c4 * B.py + M.c5 * B.vy + M.c16 * B.mz -
                      M.c18 * B.mx + M.c19 * B.vz - M.c21 * B.vx + M.c22 * B.pz -
                      M.c24 * B.px + M.c31 * B.my;
    ctype const c18 = -M.c3 * B.pw + M.c4 * B.pz + M.c5 * B.vz - M.c16 * B.my +
                      M.c17 * B.mx - M.c19 * B.vy + M.c20 * B.vx - M.c22 * B.py +
                      M.c23 * B.px + M.c31 * B.mz;
    ctype const c19 = M.c2 * B.pz - M.c3 * B.py + M.c5 * B.mx - M.c17 * B.pz +
                      M.c18 * B.py - M.c19 * B.pw - M.c20 * B.mz + M.c21 * B.my -
                      M.c25 * B.px + M.c31 * B.px;
    ctype const c20 = -M.c1 * B.pz + M.c3 * B.px + M.c5 * B.my + M.c16 * B.pz -
                      M.c18 * B.px + M.c19 * B.mz - M.c20 * B.pw - M.c21 * B.mx -
                      M.c25 * B.py + M.c31 * B.py;
    ctype const c21 = M.c1 * B.py - M.c2 * B.px + M.c5 * B.mz - M.c16 * B.py +
                      M.c17 * B.px - M.c19 * B.my + M.c20 * B.mx - M.c21 * B.pw -
                      M.c25 * B.pz + M.c31 * B.pz;
    ctype const c22 = -M.c2 * B.vz + M.c3 * B.vy + M.c4 * B.mx - M.c17 * B.vz +
                      M.c18 * B.vy + M.c22 * B.pw - M.c23 * B.mz + M.c24 * B.my +
                      M.c25 * B.vx + M.c31 * B.vx;
    ctype const c23 = M.c1 * B.vz - M.c3 * B.vx + M.c4 * B.my + M.c16 * B.vz -
                      M.c18 * B.vx + M.c22 * B.mz + M.c23 * B.pw - M.c24 * B.mx +
                      M.c25 * B.vy + M.c31 * B.vy;
    ctype const c24 = -M.c1 * B.vy + M.c2 * B.vx + M.c4 * B.mz - M.c16 * B.vy +
                      M.c17 * B.vx - M.c22 * B.my + M.c23 * B.mx + M.c24 * B.pw +
                      M.c25 * B.vz + M.c31 * B.vz;
    ctype const c25 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz + M.c19 * B.vx +
                      M.c20 * B.vy + M.c21 * B.vz - M.c22 * B.px - M.c23 * B.py -
                      M.c24 * B.pz - M.c31 * B.pw;
    ctype const c26 = M.c7 * B.pz - M.c8 * B.py + M.c9 * B.pw - M.c13 * B.vz +
                      M.c14 * B.vy + M.c15 * B.mx - M.c27 * B.mz + M.c28 * B.my +
                      M.c29 * B.vx - M.c30 * B.px;
    ctype const c27 = -M.c6 * B.pz + M.c8 * B.px + M.c10 * B.pw + M.c12 * B.vz -
                      M.c14 * B.vx + M.c15 * B.my + M.c26 * B.mz - M.c28 * B.mx +
                      M.c29 * B.vy - M.c30 * B.py;
    ctype const c28 = M.c6 * B.py - M.c7 * B.px + M.c11 * B.pw - M.c12 * B.vy +
                      M.c13 * B.vx + M.c15 * B.mz - M.c26 * B.my + M.c27 * B.mx +
                      M.c29 * B.vz - M.c30 * B.pz;
    ctype const c29 = -M.c9 * B.px - M.c10 * B.py - M.c11 * B.pz - M.c12 * B.mx -
                      M.c13 * B.my - M.c14 * B.mz - M.c26 * B.px - M.c27 * B.py -
                      M.c28 * B.pz - M.c29 * B.pw;
    ctype const c30 = -M.c6 * B.mx - M.c7 * B.my - M.c8 * B.mz - M.c9 * B.vx -
                      M.c10 * B.vy - M.c11 * B.vz + M.c26 * B.vx + M.c27 * B.vy +
                      M.c28 * B.vz + M.c30 * B.pw;
    ctype const c31 = -M.c16 * B.mx - M.c17 * B.my - M.c18 * B.mz - M.c19 * B.vx -
                      M.c20 * B.vy - M.c21 * B.vz - M.c22 * B.px - M.c23 * B.py -
                      M.c24 * B.pz - M.c25 * B.pw;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c12 - B.vy * M.c13 - B.vz * M.c14 - B.mx * M.c9 -
                     B.my * M.c10 - B.mz * M.c11 - B.px * M.c6 - B.py * M.c7 -
                     B.pz * M.c8 + B.pw * M.c15;
    ctype const c1 = B.vx * M.c5 - B.vy * M.c21 + B.vz * M.c20 + B.mx * M.c25 -
                     B.my * M.c3 + B.mz * M.c2 - B.px * M.c4 + B.py * M.c24 -
                     B.pz * M.c23 - B.pw * M.c16;
    ctype const c2 = B.vx * M.c21 + B.vy * M.c5 - B.vz * M.c19 + B.mx * M.c3 +
                     B.my * M.c25 - B.mz * M.c1 - B.px * M.c24 - B.py * M.c4 +
                     B.pz * M.c22 - B.pw * M.c17;
    ctype const c3 = -B.vx * M.c20 + B.vy * M.c19 + B.vz * M.c5 - B.mx * M.c2 +
                     B.my * M.c1 + B.mz * M.c25 + B.px * M.c23 - B.py * M.c22 -
                     B.pz * M.c4 - B.pw * M.c18;
    ctype const c4 = B.vx * M.c1 - B.vx * M.c16 + B.vy * M.c2 - B.vy * M.c17 +
                     B.vz * M.c3 - B.vz * M.c18 - B.mx * M.c22 - B.my * M.c23 -
                     B.mz * M.c24 - B.pw * M.c4;
    ctype const c5 = -B.mx * M.c19 - B.my * M.c20 - B.mz * M.c21 - B.px * M.c1 -
                     B.px * M.c16 - B.py * M.c2 - B.py * M.c17 - B.pz * M.c3 -
                     B.pz * M.c18 + B.pw * M.c5;
    ctype const c6 = B.vx * M.c0 + B.vx * M.c15 - B.vy * M.c11 - B.vy * M.c28 +
                     B.vz * M.c10 + B.vz * M.c27 + B.mx * M.c30 - B.my * M.c8 +
                     B.mz * M.c7 - B.pw * M.c6;
    ctype const c7 = B.vx * M.c11 + B.vx * M.c28 + B.vy * M.c0 + B.vy * M.c15 -
                     B.vz * M.c9 - B.vz * M.c26 + B.mx * M.c8 + B.my * M.c30 -
                     B.mz * M.c6 - B.pw * M.c7;
    ctype const c8 = -B.vx * M.c10 - B.vx * M.c27 + B.vy * M.c9 + B.vy * M.c26 +
                     B.vz * M.c0 + B.vz * M.c15 - B.mx * M.c7 + B.my * M.c6 +
                     B.mz * M.c30 - B.pw * M.c8;
    ctype const c9 = B.vx * M.c29 - B.vy * M.c14 + B.vz * M.c13 + B.mx * M.c0 -
                     B.my * M.c11 + B.mz * M.c10 + B.px * M.c30 - B.py * M.c8 +
                     B.pz * M.c7 + B.pw * M.c26;
    ctype const c10 = B.vx * M.c14 + B.vy * M.c29 - B.vz * M.c12 + B.mx * M.c11 +
                      B.my * M.c0 - B.mz * M.c9 + B.px * M.c8 + B.py * M.c30 -
                      B.pz * M.c6 + B.pw * M.c27;
    ctype const c11 = -B.vx * M.c13 + B.vy * M.c12 + B.vz * M.c29 - B.mx * M.c10 +
                      B.my * M.c9 + B.mz * M.c0 - B.px * M.c7 + B.py * M.c6 +
                      B.pz * M.c30 + B.pw * M.c28;
    ctype const c12 = B.mx * M.c29 - B.my * M.c14 + B.mz * M.c13 + B.px * M.c0 -
                      B.px * M.c15 - B.py * M.c11 + B.py * M.c28 + B.pz * M.c10 -
                      B.pz * M.c27 + B.pw * M.c12;
    ctype const c13 = B.mx * M.c14 + B.my * M.c29 - B.mz * M.c12 + B.px * M.c11 -
                      B.px * M.c28 + B.py * M.c0 - B.py * M.c15 - B.pz * M.c9 +
                      B.pz * M.c26 + B.pw * M.c13;
    ctype const c14 = -B.mx * M.c13 + B.my * M.c12 + B.mz * M.c29 - B.px * M.c10 +
                      B.px * M.c27 + B.py * M.c9 - B.py * M.c26 + B.pz * M.c0 -
                      B.pz * M.c15 + B.pw * M.c14;
    ctype const c15 = B.vx * M.c12 + B.vy * M.c13 + B.vz * M.c14 - B.mx * M.c26 -
                      B.my * M.c27 - B.mz * M.c28 - B.px * M.c6 - B.py * M.c7 -
                      B.pz * M.c8 + B.pw * M.c0;
    ctype const c16 = B.vx * M.c5 - B.vy * M.c21 + B.vz * M.c20 + B.mx * M.c31 -
                      B.my * M.c18 + B.mz * M.c17 + B.px * M.c4 - B.py * M.c24 +
                      B.pz * M.c23 - B.pw * M.c1;
    ctype const c17 = B.vx * M.c21 + B.vy * M.c5 - B.vz * M.c19 + B.mx * M.c18 +
                      B.my * M.c31 - B.mz * M.c16 + B.px * M.c24 + B.py * M.c4 -
                      B.pz * M.c22 - B.pw * M.c2;
    ctype const c18 = -B.vx * M.c20 + B.vy * M.c19 + B.vz * M.c5 - B.mx * M.c17 +
                      B.my * M.c16 + B.mz * M.c31 - B.px * M.c23 + B.py * M.c22 +
                      B.pz * M.c4 - B.pw * M.c3;
    ctype const c19 = B.mx * M.c5 - B.my * M.c21 + B.mz * M.c20 + B.px * M.c25 +
                      B.px * M.c31 - B.py * M.c3 - B.py * M.c18 + B.pz * M.c2 +
                      B.pz * M.c17 + B.pw * M.c19;
    ctype const c20 = B.mx * M.c21 + B.my * M.c5 - B.mz * M.c19 + B.px * M.c3 +
                      B.px * M.c18 + B.py * M.c25 + B.py * M.c31 - B.pz * M.c1 -
                      B.pz * M.c16 + B.pw * M.c20;
    ctype const c21 = -B.mx * M.c20 + B.my * M.c19 + B.mz * M.c5 - B.px * M.c2 -
                      B.px * M.c17 + B.py * M.c1 + B.py * M.c16 + B.pz * M.c25 +
                      B.pz * M.c31 + B.pw * M.c21;
    ctype const c22 = -B.vx * M.c25 + B.vx * M.c31 + B.vy * M.c3 - B.vy * M.c18 -
                      B.vz * M.c2 + B.vz * M.c17 + B.mx * M.c4 - B.my * M.c24 +
                      B.mz * M.c23 - B.pw * M.c22;
    ctype const c23 = -B.vx * M.c3 + B.vx * M.c18 - B.vy * M.c25 + B.vy * M.c31 +
                      B.vz * M.c1 - B.vz * M.c16 + B.mx * M.c24 + B.my * M.c4 -
                      B.mz * M.c22 - B.pw * M.c23;
    ctype const c24 = B.vx * M.c2 - B.vx * M.c17 - B.vy * M.c1 + B.vy * M.c16 -
                      B.vz * M.c25 + B.vz * M.c31 - B.mx * M.c23 + B.my * M.c22 +
                      B.mz * M.c4 - B.pw * M.c24;
    ctype const c25 = -B.vx * M.c19 - B.vy * M.c20 - B.vz * M.c21 - B.mx * M.c1 -
                      B.my * M.c2 - B.mz * M.c3 + B.px * M.c22 + B.py * M.c23 +
                      B.pz * M.c24 - B.pw * M.c31;
    ctype const c26 = -B.vx * M.c29 + B.vy * M.c14 - B.vz * M.c13 + B.mx * M.c15 -
                      B.my * M.c28 + B.mz * M.c27 + B.px * M.c30 - B.py * M.c8 +
                      B.pz * M.c7 + B.pw * M.c9;
    ctype const c27 = -B.vx * M.c14 - B.vy * M.c29 + B.vz * M.c12 + B.mx * M.c28 +
                      B.my * M.c15 - B.mz * M.c26 + B.px * M.c8 + B.py * M.c30 -
                      B.pz * M.c6 + B.pw * M.c10;
    ctype const c28 = B.vx * M.c13 - B.vy * M.c12 - B.vz * M.c29 - B.mx * M.c27 +
                      B.my * M.c26 + B.mz * M.c15 - B.px * M.c7 + B.py * M.c6 +
                      B.pz * M.c30 + B.pw * M.c11;
    ctype const c29 = -B.mx * M.c12 - B.my * M.c13 - B.mz * M.c14 - B.px * M.c9 +
                      B.px * M.c26 - B.py * M.c10 + B.py * M.c27 - B.pz * M.c11 +
                      B.pz * M.c28 + B.pw * M.c29;
    ctype const c30 = -B.vx * M.c9 - B.vx * M.c26 - B.vy * M.c10 - B.vy * M.c27 -
                      B.vz * M.c11 - B.vz * M.c28 - B.mx * M.c6 - B.my * M.c7 -
                      B.mz * M.c8 - B.pw * M.c30;
    ctype const c31 = -B.vx * M.c19 - B.vy * M.c20 - B.vz * M.c21 - B.mx * M.c16 -
                      B.my * M.c17 - B.mz * M.c18 - B.px * M.c22 - B.py * M.c23 -
                      B.pz * M.c24 - B.pw * M.c25;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& M,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z - M.c4 * v.u - M.c5 * v.w;
    ctype const c1 = M.c0 * v.x + M.c6 * v.u - M.c10 * v.z + M.c11 * v.y - M.c12 * v.w;
    ctype const c2 = M.c0 * v.y + M.c7 * v.u + M.c9 * v.z - M.c11 * v.x - M.c13 * v.w;
    ctype const c3 = M.c0 * v.z + M.c8 * v.u - M.c9 * v.y + M.c10 * v.x - M.c14 * v.w;
    ctype const c4 = M.c0 * v.w + M.c6 * v.x + M.c7 * v.y + M.c8 * v.z - M.c15 * v.w;
    ctype const c5 = M.c0 * v.u - M.c12 * v.x - M.c13 * v.y - M.c14 * v.z + M.c15 * v.u;
    ctype const c6 = -M.c1 * v.w + M.c4 * v.x - M.c16 * v.w - M.c23 * v.z + M.c24 * v.y;
    ctype const c7 = -M.c2 * v.w + M.c4 * v.y - M.c17 * v.w + M.c22 * v.z - M.c24 * v.x;
    ctype const c8 = -M.c3 * v.w + M.c4 * v.z - M.c18 * v.w - M.c22 * v.y + M.c23 * v.x;
    ctype const c9 = M.c2 * v.z - M.c3 * v.y - M.c19 * v.w - M.c22 * v.u - M.c25 * v.x;
    ctype const c10 = -M.c1 * v.z + M.c3 * v.x - M.c20 * v.w - M.c23 * v.u - M.c25 * v.y;
    ctype const c11 = M.c1 * v.y - M.c2 * v.x - M.c21 * v.w - M.c24 * v.u - M.c25 * v.z;
    ctype const c12 = M.c1 * v.u - M.c5 * v.x - M.c16 * v.u + M.c20 * v.z - M.c21 * v.y;
    ctype const c13 = M.c2 * v.u - M.c5 * v.y - M.c17 * v.u - M.c19 * v.z + M.c21 * v.x;
    ctype const c14 = M.c3 * v.u - M.c5 * v.z - M.c18 * v.u + M.c19 * v.y - M.c20 * v.x;
    ctype const c15 = M.c4 * v.u - M.c5 * v.w - M.c16 * v.x - M.c17 * v.y - M.c18 * v.z;
    ctype const c16 = M.c6 * v.u + M.c12 * v.w - M.c15 * v.x + M.c27 * v.z - M.c28 * v.y;
    ctype const c17 = M.c7 * v.u + M.c13 * v.w - M.c15 * v.y - M.c26 * v.z + M.c28 * v.x;
    ctype const c18 = M.c8 * v.u + M.c14 * v.w - M.c15 * v.z + M.c26 * v.y - M.c27 * v.x;
    ctype const c19 = M.c9 * v.u - M.c13 * v.z + M.c14 * v.y + M.c26 * v.u + M.c29 * v.x;
    ctype const c20 = M.c10 * v.u + M.c12 * v.z - M.c14 * v.x + M.c27 * v.u + M.c29 * v.y;
    ctype const c21 = M.c11 * v.u - M.c12 * v.y + M.c13 * v.x + M.c28 * v.u + M.c29 * v.z;
    ctype const c22 = M.c7 * v.z - M.c8 * v.y + M.c9 * v.w - M.c26 * v.w - M.c30 * v.x;
    ctype const c23 = -M.c6 * v.z + M.c8 * v.x + M.c10 * v.w - M.c27 * v.w - M.c30 * v.y;
    ctype const c24 = M.c6 * v.y - M.c7 * v.x + M.c11 * v.w - M.c28 * v.w - M.c30 * v.z;
    ctype const c25 = -M.c9 * v.x - M.c10 * v.y - M.c11 * v.z - M.c29 * v.w + M.c30 * v.u;
    ctype const c26 =
        -M.c17 * v.z + M.c18 * v.y - M.c19 * v.w + M.c22 * v.u + M.c31 * v.x;
    ctype const c27 = M.c16 * v.z - M.c18 * v.x - M.c20 * v.w + M.c23 * v.u + M.c31 * v.y;
    ctype const c28 =
        -M.c16 * v.y + M.c17 * v.x - M.c21 * v.w + M.c24 * v.u + M.c31 * v.z;
    ctype const c29 = M.c19 * v.x + M.c20 * v.y + M.c21 * v.z + M.c25 * v.u - M.c31 * v.u;
    ctype const c30 =
        -M.c22 * v.x - M.c23 * v.y - M.c24 * v.z - M.c25 * v.w - M.c31 * v.w;
    ctype const c31 = M.c26 * v.x + M.c27 * v.y + M.c28 * v.z + M.c29 * v.w + M.c30 * v.u;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3 - v.w * M.c5 - v.u * M.c4;
    ctype const c1 = v.x * M.c0 - v.y * M.c11 + v.z * M.c10 + v.w * M.c12 - v.u * M.c6;
    ctype const c2 = v.x * M.c11 + v.y * M.c0 - v.z * M.c9 + v.w * M.c13 - v.u * M.c7;
    ctype const c3 = -v.x * M.c10 + v.y * M.c9 + v.z * M.c0 + v.w * M.c14 - v.u * M.c8;
    ctype const c4 = -v.x * M.c6 - v.y * M.c7 - v.z * M.c8 + v.w * M.c0 + v.w * M.c15;
    ctype const c5 = v.x * M.c12 + v.y * M.c13 + v.z * M.c14 + v.u * M.c0 - v.u * M.c15;
    ctype const c6 = -v.x * M.c4 + v.y * M.c24 - v.z * M.c23 + v.w * M.c1 - v.w * M.c16;
    ctype const c7 = -v.x * M.c24 - v.y * M.c4 + v.z * M.c22 + v.w * M.c2 - v.w * M.c17;
    ctype const c8 = v.x * M.c23 - v.y * M.c22 - v.z * M.c4 + v.w * M.c3 - v.w * M.c18;
    ctype const c9 = -v.x * M.c25 + v.y * M.c3 - v.z * M.c2 - v.w * M.c19 - v.u * M.c22;
    ctype const c10 = -v.x * M.c3 - v.y * M.c25 + v.z * M.c1 - v.w * M.c20 - v.u * M.c23;
    ctype const c11 = v.x * M.c2 - v.y * M.c1 - v.z * M.c25 - v.w * M.c21 - v.u * M.c24;
    ctype const c12 = v.x * M.c5 - v.y * M.c21 + v.z * M.c20 - v.u * M.c1 - v.u * M.c16;
    ctype const c13 = v.x * M.c21 + v.y * M.c5 - v.z * M.c19 - v.u * M.c2 - v.u * M.c17;
    ctype const c14 = -v.x * M.c20 + v.y * M.c19 + v.z * M.c5 - v.u * M.c3 - v.u * M.c18;
    ctype const c15 = -v.x * M.c16 - v.y * M.c17 - v.z * M.c18 + v.w * M.c5 - v.u * M.c4;
    ctype const c16 = -v.x * M.c15 + v.y * M.c28 - v.z * M.c27 + v.w * M.c12 + v.u * M.c6;
    ctype const c17 = -v.x * M.c28 - v.y * M.c15 + v.z * M.c26 + v.w * M.c13 + v.u * M.c7;
    ctype const c18 = v.x * M.c27 - v.y * M.c26 - v.z * M.c15 + v.w * M.c14 + v.u * M.c8;
    ctype const c19 = -v.x * M.c29 + v.y * M.c14 - v.z * M.c13 + v.u * M.c9 - v.u * M.c26;
    ctype const c20 =
        -v.x * M.c14 - v.y * M.c29 + v.z * M.c12 + v.u * M.c10 - v.u * M.c27;
    ctype const c21 = v.x * M.c13 - v.y * M.c12 - v.z * M.c29 + v.u * M.c11 - v.u * M.c28;
    ctype const c22 = v.x * M.c30 - v.y * M.c8 + v.z * M.c7 + v.w * M.c9 + v.w * M.c26;
    ctype const c23 = v.x * M.c8 + v.y * M.c30 - v.z * M.c6 + v.w * M.c10 + v.w * M.c27;
    ctype const c24 = -v.x * M.c7 + v.y * M.c6 + v.z * M.c30 + v.w * M.c11 + v.w * M.c28;
    ctype const c25 = -v.x * M.c9 - v.y * M.c10 - v.z * M.c11 + v.w * M.c29 - v.u * M.c30;
    ctype const c26 = v.x * M.c31 - v.y * M.c18 + v.z * M.c17 + v.w * M.c19 - v.u * M.c22;
    ctype const c27 = v.x * M.c18 + v.y * M.c31 - v.z * M.c16 + v.w * M.c20 - v.u * M.c23;
    ctype const c28 =
        -v.x * M.c17 + v.y * M.c16 + v.z * M.c31 + v.w * M.c21 - v.u * M.c24;
    ctype const c29 =
        -v.x * M.c19 - v.y * M.c20 - v.z * M.c21 - v.u * M.c25 - v.u * M.c31;
    ctype const c30 = v.x * M.c22 + v.y * M.c23 + v.z * M.c24 + v.w * M.c25 - v.w * M.c31;
    ctype const c31 = v.x * M.c26 + v.y * M.c27 + v.z * M.c28 + v.w * M.c29 + v.u * M.c30;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(MVec3dc<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    ctype const c16 = M.c16 * ctype(s);
    ctype const c17 = M.c17 * ctype(s);
    ctype const c18 = M.c18 * ctype(s);
    ctype const c19 = M.c19 * ctype(s);
    ctype const c20 = M.c20 * ctype(s);
    ctype const c21 = M.c21 * ctype(s);
    ctype const c22 = M.c22 * ctype(s);
    ctype const c23 = M.c23 * ctype(s);
    ctype const c24 = M.c24 * ctype(s);
    ctype const c25 = M.c25 * ctype(s);
    ctype const c26 = M.c26 * ctype(s);
    ctype const c27 = M.c27 * ctype(s);
    ctype const c28 = M.c28 * ctype(s);
    ctype const c29 = M.c29 * ctype(s);
    ctype const c30 = M.c30 * ctype(s);
    ctype const c31 = M.c31 * ctype(s);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator*(Scalar3dc<T> s, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    ctype const c16 = ctype(s) * M.c16;
    ctype const c17 = ctype(s) * M.c17;
    ctype const c18 = ctype(s) * M.c18;
    ctype const c19 = ctype(s) * M.c19;
    ctype const c20 = ctype(s) * M.c20;
    ctype const c21 = ctype(s) * M.c21;
    ctype const c22 = ctype(s) * M.c22;
    ctype const c23 = ctype(s) * M.c23;
    ctype const c24 = ctype(s) * M.c24;
    ctype const c25 = ctype(s) * M.c25;
    ctype const c26 = ctype(s) * M.c26;
    ctype const c27 = ctype(s) * M.c27;
    ctype const c28 = ctype(s) * M.c28;
    ctype const c29 = ctype(s) * M.c29;
    ctype const c30 = ctype(s) * M.c30;
    ctype const c31 = ctype(s) * M.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc gpr :: gpr(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& A,
                                                        MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c7 - A.c2 * B.c8 - A.c3 * B.c9 - A.c4 * B.c4 -
                     A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c1 - A.c8 * B.c2 - A.c9 * B.c3 +
                     A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c12 - A.c13 * B.c13 +
                     A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 + A.c1 * B.c10 - A.c2 * B.c6 -
                     A.c2 * B.c13 + A.c3 * B.c5 + A.c3 * B.c12 + A.c4 * B.c15 -
                     A.c5 * B.c3 + A.c6 * B.c2 - A.c10 * B.c1 - A.c11 * B.c15 +
                     A.c12 * B.c3 - A.c13 * B.c2 + A.c15 * B.c4 + A.c15 * B.c11;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c1 * B.c13 + A.c2 * B.c0 +
                     A.c2 * B.c10 - A.c3 * B.c4 - A.c3 * B.c11 + A.c4 * B.c3 +
                     A.c5 * B.c15 - A.c6 * B.c1 - A.c10 * B.c2 - A.c11 * B.c3 -
                     A.c12 * B.c15 + A.c13 * B.c1 + A.c15 * B.c5 + A.c15 * B.c12;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c5 - A.c1 * B.c12 + A.c2 * B.c4 +
                     A.c2 * B.c11 + A.c3 * B.c0 + A.c3 * B.c10 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c15 - A.c10 * B.c3 + A.c11 * B.c2 -
                     A.c12 * B.c1 - A.c13 * B.c15 + A.c15 * B.c6 + A.c15 * B.c13;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c14 - A.c2 * B.c9 + A.c3 * B.c8 +
                     A.c4 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c15 -
                     A.c8 * B.c3 + A.c9 * B.c2 + A.c10 * B.c11 + A.c11 * B.c10 -
                     A.c12 * B.c13 + A.c13 * B.c12 + A.c14 * B.c1 + A.c15 * B.c7;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c9 + A.c2 * B.c14 - A.c3 * B.c7 +
                     A.c4 * B.c6 + A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c3 +
                     A.c8 * B.c15 - A.c9 * B.c1 + A.c10 * B.c12 + A.c11 * B.c13 +
                     A.c12 * B.c10 - A.c13 * B.c11 + A.c14 * B.c2 + A.c15 * B.c8;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c8 + A.c2 * B.c7 + A.c3 * B.c14 -
                     A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c0 - A.c7 * B.c2 + A.c8 * B.c1 +
                     A.c9 * B.c15 + A.c10 * B.c13 - A.c11 * B.c12 + A.c12 * B.c11 +
                     A.c13 * B.c10 + A.c14 * B.c3 + A.c15 * B.c9;
    ctype const c7 = A.c0 * B.c7 + A.c4 * B.c14 - A.c5 * B.c9 + A.c6 * B.c8 +
                     A.c7 * B.c0 - A.c7 * B.c10 - A.c8 * B.c6 + A.c8 * B.c13 +
                     A.c9 * B.c5 - A.c9 * B.c12 + A.c10 * B.c7 + A.c11 * B.c14 -
                     A.c12 * B.c9 + A.c13 * B.c8 + A.c14 * B.c4 - A.c14 * B.c11;
    ctype const c8 = A.c0 * B.c8 + A.c4 * B.c9 + A.c5 * B.c14 - A.c6 * B.c7 +
                     A.c7 * B.c6 - A.c7 * B.c13 + A.c8 * B.c0 - A.c8 * B.c10 -
                     A.c9 * B.c4 + A.c9 * B.c11 + A.c10 * B.c8 + A.c11 * B.c9 +
                     A.c12 * B.c14 - A.c13 * B.c7 + A.c14 * B.c5 - A.c14 * B.c12;
    ctype const c9 = A.c0 * B.c9 - A.c4 * B.c8 + A.c5 * B.c7 + A.c6 * B.c14 -
                     A.c7 * B.c5 + A.c7 * B.c12 + A.c8 * B.c4 - A.c8 * B.c11 +
                     A.c9 * B.c0 - A.c9 * B.c10 + A.c10 * B.c9 - A.c11 * B.c8 +
                     A.c12 * B.c7 + A.c13 * B.c14 + A.c14 * B.c6 - A.c14 * B.c13;
    ctype const c10 = A.c0 * B.c10 + A.c1 * B.c7 + A.c2 * B.c8 + A.c3 * B.c9 -
                      A.c4 * B.c11 - A.c5 * B.c12 - A.c6 * B.c13 - A.c7 * B.c1 -
                      A.c8 * B.c2 - A.c9 * B.c3 + A.c10 * B.c0 - A.c11 * B.c4 -
                      A.c12 * B.c5 - A.c13 * B.c6 + A.c14 * B.c15 - A.c15 * B.c14;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c14 + A.c2 * B.c9 - A.c3 * B.c8 +
                      A.c4 * B.c10 - A.c5 * B.c13 + A.c6 * B.c12 + A.c7 * B.c15 -
                      A.c8 * B.c3 + A.c9 * B.c2 + A.c10 * B.c4 + A.c11 * B.c0 -
                      A.c12 * B.c6 + A.c13 * B.c5 + A.c14 * B.c1 - A.c15 * B.c7;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c9 - A.c2 * B.c14 + A.c3 * B.c7 +
                      A.c4 * B.c13 + A.c5 * B.c10 - A.c6 * B.c11 + A.c7 * B.c3 +
                      A.c8 * B.c15 - A.c9 * B.c1 + A.c10 * B.c5 + A.c11 * B.c6 +
                      A.c12 * B.c0 - A.c13 * B.c4 + A.c14 * B.c2 - A.c15 * B.c8;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c8 - A.c2 * B.c7 - A.c3 * B.c14 -
                      A.c4 * B.c12 + A.c5 * B.c11 + A.c6 * B.c10 - A.c7 * B.c2 +
                      A.c8 * B.c1 + A.c9 * B.c15 + A.c10 * B.c6 - A.c11 * B.c5 +
                      A.c12 * B.c4 + A.c13 * B.c0 + A.c14 * B.c3 - A.c15 * B.c9;
    ctype const c14 = A.c0 * B.c14 - A.c4 * B.c7 - A.c5 * B.c8 - A.c6 * B.c9 -
                      A.c7 * B.c4 + A.c7 * B.c11 - A.c8 * B.c5 + A.c8 * B.c12 -
                      A.c9 * B.c6 + A.c9 * B.c13 + A.c10 * B.c14 - A.c11 * B.c7 -
                      A.c12 * B.c8 - A.c13 * B.c9 + A.c14 * B.c0 - A.c14 * B.c10;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c4 - A.c1 * B.c11 - A.c2 * B.c5 -
                      A.c2 * B.c12 - A.c3 * B.c6 - A.c3 * B.c13 - A.c4 * B.c1 -
                      A.c5 * B.c2 - A.c6 * B.c3 - A.c10 * B.c15 + A.c11 * B.c1 +
                      A.c12 * B.c2 + A.c13 * B.c3 + A.c15 * B.c0 + A.c15 * B.c10;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& A,
                                                        MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c4 - A.c2 * B.c10 + A.c3 * B.c9 +
                     A.c4 * B.c14 - A.c5 * B.c2 + A.c6 * B.c1 - A.c7 * B.c3 +
                     A.c8 * B.c13 - A.c9 * B.c12 - A.c10 * B.c5 - A.c11 * B.c15 +
                     A.c12 * B.c7 - A.c13 * B.c6 - A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c10 + A.c2 * B.c4 - A.c3 * B.c8 +
                     A.c4 * B.c2 + A.c5 * B.c14 - A.c6 * B.c0 - A.c7 * B.c13 -
                     A.c8 * B.c3 + A.c9 * B.c11 - A.c10 * B.c6 - A.c11 * B.c7 -
                     A.c12 * B.c15 + A.c13 * B.c5 - A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c4 - A.c4 * B.c1 +
                     A.c5 * B.c0 + A.c6 * B.c14 + A.c7 * B.c12 - A.c8 * B.c11 -
                     A.c9 * B.c3 - A.c10 * B.c7 + A.c11 * B.c6 - A.c12 * B.c5 -
                     A.c13 * B.c15 - A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c0 - A.c1 * B.c5 + A.c2 * B.c1 - A.c2 * B.c6 +
                     A.c3 * B.c2 - A.c3 * B.c7 - A.c4 * B.c11 - A.c5 * B.c12 -
                     A.c6 * B.c13 - A.c10 * B.c3 + A.c11 * B.c11 + A.c12 * B.c12 +
                     A.c13 * B.c13 - A.c15 * B.c14 + A.c15 * B.c15;
    ctype const c4 = A.c0 * B.c4 - A.c4 * B.c8 - A.c5 * B.c9 - A.c6 * B.c10 -
                     A.c7 * B.c0 - A.c7 * B.c5 - A.c8 * B.c1 - A.c8 * B.c6 - A.c9 * B.c2 -
                     A.c9 * B.c7 + A.c10 * B.c4 - A.c11 * B.c8 - A.c12 * B.c9 -
                     A.c13 * B.c10 + A.c14 * B.c14 + A.c14 * B.c15;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c4 - A.c2 * B.c10 + A.c3 * B.c9 +
                     A.c4 * B.c15 - A.c5 * B.c7 + A.c6 * B.c6 + A.c7 * B.c3 -
                     A.c8 * B.c13 + A.c9 * B.c12 - A.c10 * B.c0 - A.c11 * B.c14 +
                     A.c12 * B.c2 - A.c13 * B.c1 + A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c10 + A.c2 * B.c4 - A.c3 * B.c8 +
                     A.c4 * B.c7 + A.c5 * B.c15 - A.c6 * B.c5 + A.c7 * B.c13 +
                     A.c8 * B.c3 - A.c9 * B.c11 - A.c10 * B.c1 - A.c11 * B.c2 -
                     A.c12 * B.c14 + A.c13 * B.c0 + A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c4 - A.c4 * B.c6 +
                     A.c5 * B.c5 + A.c6 * B.c15 - A.c7 * B.c12 + A.c8 * B.c11 +
                     A.c9 * B.c3 - A.c10 * B.c2 + A.c11 * B.c1 - A.c12 * B.c0 -
                     A.c13 * B.c14 + A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c8 = A.c0 * B.c8 + A.c4 * B.c4 - A.c5 * B.c10 + A.c6 * B.c9 +
                     A.c7 * B.c14 + A.c7 * B.c15 - A.c8 * B.c2 - A.c8 * B.c7 +
                     A.c9 * B.c1 + A.c9 * B.c6 + A.c10 * B.c8 + A.c11 * B.c4 -
                     A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c0 + A.c14 * B.c5;
    ctype const c9 = A.c0 * B.c9 + A.c4 * B.c10 + A.c5 * B.c4 - A.c6 * B.c8 +
                     A.c7 * B.c2 + A.c7 * B.c7 + A.c8 * B.c14 + A.c8 * B.c15 -
                     A.c9 * B.c0 - A.c9 * B.c5 + A.c10 * B.c9 + A.c11 * B.c10 +
                     A.c12 * B.c4 - A.c13 * B.c8 + A.c14 * B.c1 + A.c14 * B.c6;
    ctype const c10 = A.c0 * B.c10 - A.c4 * B.c9 + A.c5 * B.c8 + A.c6 * B.c4 -
                      A.c7 * B.c1 - A.c7 * B.c6 + A.c8 * B.c0 + A.c8 * B.c5 +
                      A.c9 * B.c14 + A.c9 * B.c15 + A.c10 * B.c10 - A.c11 * B.c9 +
                      A.c12 * B.c8 + A.c13 * B.c4 + A.c14 * B.c2 + A.c14 * B.c7;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c14 + A.c1 * B.c15 + A.c2 * B.c2 -
                      A.c2 * B.c7 - A.c3 * B.c1 + A.c3 * B.c6 + A.c4 * B.c3 -
                      A.c5 * B.c13 + A.c6 * B.c12 - A.c10 * B.c11 - A.c11 * B.c3 +
                      A.c12 * B.c13 - A.c13 * B.c12 - A.c15 * B.c0 + A.c15 * B.c5;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c2 + A.c1 * B.c7 - A.c2 * B.c14 +
                      A.c2 * B.c15 + A.c3 * B.c0 - A.c3 * B.c5 + A.c4 * B.c13 +
                      A.c5 * B.c3 - A.c6 * B.c11 - A.c10 * B.c12 - A.c11 * B.c13 -
                      A.c12 * B.c3 + A.c13 * B.c11 - A.c15 * B.c1 + A.c15 * B.c6;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c1 - A.c1 * B.c6 - A.c2 * B.c0 +
                      A.c2 * B.c5 - A.c3 * B.c14 + A.c3 * B.c15 - A.c4 * B.c12 +
                      A.c5 * B.c11 + A.c6 * B.c3 - A.c10 * B.c13 + A.c11 * B.c12 -
                      A.c12 * B.c11 - A.c13 * B.c3 - A.c15 * B.c2 + A.c15 * B.c7;
    ctype const c14 = A.c0 * B.c14 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                      A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2 + A.c7 * B.c11 +
                      A.c8 * B.c12 + A.c9 * B.c13 - A.c10 * B.c15 + A.c11 * B.c5 +
                      A.c12 * B.c6 + A.c13 * B.c7 - A.c14 * B.c3 + A.c15 * B.c4;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                      A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 - A.c7 * B.c11 -
                      A.c8 * B.c12 - A.c9 * B.c13 - A.c10 * B.c14 + A.c11 * B.c0 +
                      A.c12 * B.c1 + A.c13 * B.c2 + A.c14 * B.c3 + A.c15 * B.c4;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& A,
                                                        MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c7 - A.c4 * B.c1 -
                     A.c5 * B.c10 + A.c6 * B.c13 - A.c7 * B.c12 - A.c8 * B.c15 +
                     A.c9 * B.c3 - A.c10 * B.c2 + A.c11 * B.c14 - A.c12 * B.c9 +
                     A.c13 * B.c8 + A.c14 * B.c4 - A.c15 * B.c11;
    ctype const c1 = A.c0 * B.c6 + A.c1 * B.c0 - A.c2 * B.c4 + A.c3 * B.c8 - A.c4 * B.c2 -
                     A.c5 * B.c13 - A.c6 * B.c10 + A.c7 * B.c11 - A.c8 * B.c3 -
                     A.c9 * B.c15 + A.c10 * B.c1 + A.c11 * B.c9 + A.c12 * B.c14 -
                     A.c13 * B.c7 + A.c14 * B.c5 - A.c15 * B.c12;
    ctype const c2 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c0 + A.c3 * B.c9 -
                     A.c4 * B.c3 + A.c5 * B.c12 - A.c6 * B.c11 - A.c7 * B.c10 +
                     A.c8 * B.c2 - A.c9 * B.c1 - A.c10 * B.c15 - A.c11 * B.c8 +
                     A.c12 * B.c7 + A.c13 * B.c14 + A.c14 * B.c6 - A.c15 * B.c13;
    ctype const c3 = -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0 +
                     A.c3 * B.c10 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 -
                     A.c11 * B.c4 - A.c11 * B.c11 - A.c12 * B.c5 - A.c12 * B.c12 -
                     A.c13 * B.c6 - A.c13 * B.c13 + A.c14 * B.c15 + A.c15 * B.c15;
    ctype const c4 = A.c0 * B.c7 + A.c1 * B.c8 + A.c2 * B.c9 + A.c4 * B.c0 -
                     A.c4 * B.c10 - A.c5 * B.c7 - A.c6 * B.c8 - A.c7 * B.c9 -
                     A.c8 * B.c4 + A.c8 * B.c11 - A.c9 * B.c5 + A.c9 * B.c12 -
                     A.c10 * B.c6 + A.c10 * B.c13 - A.c14 * B.c14 + A.c15 * B.c14;
    ctype const c5 = -A.c0 * B.c10 + A.c1 * B.c13 - A.c2 * B.c12 + A.c3 * B.c7 +
                     A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c6 + A.c7 * B.c5 +
                     A.c8 * B.c15 - A.c9 * B.c3 + A.c10 * B.c2 + A.c11 * B.c14 -
                     A.c12 * B.c9 + A.c13 * B.c8 - A.c14 * B.c11 + A.c15 * B.c4;
    ctype const c6 = -A.c0 * B.c13 - A.c1 * B.c10 + A.c2 * B.c11 + A.c3 * B.c8 +
                     A.c4 * B.c2 + A.c5 * B.c6 + A.c6 * B.c0 - A.c7 * B.c4 + A.c8 * B.c3 +
                     A.c9 * B.c15 - A.c10 * B.c1 + A.c11 * B.c9 + A.c12 * B.c14 -
                     A.c13 * B.c7 - A.c14 * B.c12 + A.c15 * B.c5;
    ctype const c7 = A.c0 * B.c12 - A.c1 * B.c11 - A.c2 * B.c10 + A.c3 * B.c9 +
                     A.c4 * B.c3 - A.c5 * B.c5 + A.c6 * B.c4 + A.c7 * B.c0 - A.c8 * B.c2 +
                     A.c9 * B.c1 + A.c10 * B.c15 - A.c11 * B.c8 + A.c12 * B.c7 +
                     A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c6;
    ctype const c8 = -A.c0 * B.c14 + A.c1 * B.c9 - A.c2 * B.c8 + A.c4 * B.c4 -
                     A.c4 * B.c11 + A.c5 * B.c14 - A.c6 * B.c9 + A.c7 * B.c8 +
                     A.c8 * B.c0 - A.c8 * B.c10 - A.c9 * B.c6 + A.c9 * B.c13 +
                     A.c10 * B.c5 - A.c10 * B.c12 - A.c14 * B.c7 + A.c15 * B.c7;
    ctype const c9 = -A.c0 * B.c9 - A.c1 * B.c14 + A.c2 * B.c7 + A.c4 * B.c5 -
                     A.c4 * B.c12 + A.c5 * B.c9 + A.c6 * B.c14 - A.c7 * B.c7 +
                     A.c8 * B.c6 - A.c8 * B.c13 + A.c9 * B.c0 - A.c9 * B.c10 -
                     A.c10 * B.c4 + A.c10 * B.c11 - A.c14 * B.c8 + A.c15 * B.c8;
    ctype const c10 = A.c0 * B.c8 - A.c1 * B.c7 - A.c2 * B.c14 + A.c4 * B.c6 -
                      A.c4 * B.c13 - A.c5 * B.c8 + A.c6 * B.c7 + A.c7 * B.c14 -
                      A.c8 * B.c5 + A.c8 * B.c12 + A.c9 * B.c4 - A.c9 * B.c11 +
                      A.c10 * B.c0 - A.c10 * B.c10 - A.c14 * B.c9 + A.c15 * B.c9;
    ctype const c11 = A.c0 * B.c15 - A.c1 * B.c3 + A.c2 * B.c2 + A.c3 * B.c4 +
                      A.c3 * B.c11 + A.c5 * B.c15 - A.c6 * B.c3 + A.c7 * B.c2 +
                      A.c11 * B.c0 + A.c11 * B.c10 - A.c12 * B.c6 - A.c12 * B.c13 +
                      A.c13 * B.c5 + A.c13 * B.c12 + A.c14 * B.c1 + A.c15 * B.c1;
    ctype const c12 = A.c0 * B.c3 + A.c1 * B.c15 - A.c2 * B.c1 + A.c3 * B.c5 +
                      A.c3 * B.c12 + A.c5 * B.c3 + A.c6 * B.c15 - A.c7 * B.c1 +
                      A.c11 * B.c6 + A.c11 * B.c13 + A.c12 * B.c0 + A.c12 * B.c10 -
                      A.c13 * B.c4 - A.c13 * B.c11 + A.c14 * B.c2 + A.c15 * B.c2;
    ctype const c13 = -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c15 + A.c3 * B.c6 +
                      A.c3 * B.c13 - A.c5 * B.c2 + A.c6 * B.c1 + A.c7 * B.c15 -
                      A.c11 * B.c5 - A.c11 * B.c12 + A.c12 * B.c4 + A.c12 * B.c11 +
                      A.c13 * B.c0 + A.c13 * B.c10 + A.c14 * B.c3 + A.c15 * B.c3;
    ctype const c14 = -A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 + A.c3 * B.c14 -
                      A.c4 * B.c15 + A.c5 * B.c11 + A.c6 * B.c12 + A.c7 * B.c13 +
                      A.c8 * B.c1 + A.c9 * B.c2 + A.c10 * B.c3 - A.c11 * B.c7 -
                      A.c12 * B.c8 - A.c13 * B.c9 + A.c14 * B.c0 - A.c15 * B.c10;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                      A.c4 * B.c15 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 -
                      A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 - A.c11 * B.c7 -
                      A.c12 * B.c8 - A.c13 * B.c9 - A.c14 * B.c10 + A.c15 * B.c0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(mv_e,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& M,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * ctype(ps);
    ctype const c1 = -M.c12 * ctype(ps);
    ctype const c2 = -M.c13 * ctype(ps);
    ctype const c3 = M.c15 * ctype(ps);
    ctype const c4 = M.c14 * ctype(ps);
    ctype const c5 = M.c4 * ctype(ps);
    ctype const c6 = M.c5 * ctype(ps);
    ctype const c7 = M.c6 * ctype(ps);
    ctype const c8 = M.c7 * ctype(ps);
    ctype const c9 = M.c8 * ctype(ps);
    ctype const c10 = M.c9 * ctype(ps);
    ctype const c11 = M.c1 * ctype(ps);
    ctype const c12 = M.c2 * ctype(ps);
    ctype const c13 = M.c3 * ctype(ps);
    ctype const c14 = -M.c10 * ctype(ps);
    ctype const c15 = M.c0 * ctype(ps);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(ps,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c11;
    ctype const c1 = -ctype(ps) * M.c12;
    ctype const c2 = -ctype(ps) * M.c13;
    ctype const c3 = ctype(ps) * M.c15;
    ctype const c4 = ctype(ps) * M.c14;
    ctype const c5 = ctype(ps) * M.c4;
    ctype const c6 = ctype(ps) * M.c5;
    ctype const c7 = ctype(ps) * M.c6;
    ctype const c8 = ctype(ps) * M.c7;
    ctype const c9 = ctype(ps) * M.c8;
    ctype const c10 = ctype(ps) * M.c9;
    ctype const c11 = ctype(ps) * M.c1;
    ctype const c12 = ctype(ps) * M.c2;
    ctype const c13 = ctype(ps) * M.c3;
    ctype const c14 = -ctype(ps) * M.c10;
    ctype const c15 = ctype(ps) * M.c0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(mv_e,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& M,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * Q.x - M.c12 * Q.y - M.c13 * Q.z + M.c14 * Q.u + M.c15 * Q.w;
    ctype const c1 = -M.c2 * Q.z + M.c3 * Q.y + M.c4 * Q.u - M.c11 * Q.u + M.c15 * Q.x;
    ctype const c2 = M.c1 * Q.z - M.c3 * Q.x + M.c5 * Q.u - M.c12 * Q.u + M.c15 * Q.y;
    ctype const c3 = -M.c1 * Q.y + M.c2 * Q.x + M.c6 * Q.u - M.c13 * Q.u + M.c15 * Q.z;
    ctype const c4 = M.c1 * Q.w + M.c7 * Q.u + M.c10 * Q.x - M.c12 * Q.z + M.c13 * Q.y;
    ctype const c5 = M.c2 * Q.w + M.c8 * Q.u + M.c10 * Q.y + M.c11 * Q.z - M.c13 * Q.x;
    ctype const c6 = M.c3 * Q.w + M.c9 * Q.u + M.c10 * Q.z - M.c11 * Q.y + M.c12 * Q.x;
    ctype const c7 = M.c4 * Q.w + M.c8 * Q.z - M.c9 * Q.y + M.c11 * Q.w - M.c14 * Q.x;
    ctype const c8 = M.c5 * Q.w - M.c7 * Q.z + M.c9 * Q.x + M.c12 * Q.w - M.c14 * Q.y;
    ctype const c9 = M.c6 * Q.w + M.c7 * Q.y - M.c8 * Q.x + M.c13 * Q.w - M.c14 * Q.z;
    ctype const c10 = -M.c4 * Q.x - M.c5 * Q.y - M.c6 * Q.z + M.c14 * Q.u - M.c15 * Q.w;
    ctype const c11 = M.c0 * Q.x - M.c1 * Q.w - M.c5 * Q.z + M.c6 * Q.y + M.c7 * Q.u;
    ctype const c12 = M.c0 * Q.y - M.c2 * Q.w + M.c4 * Q.z - M.c6 * Q.x + M.c8 * Q.u;
    ctype const c13 = M.c0 * Q.z - M.c3 * Q.w - M.c4 * Q.y + M.c5 * Q.x + M.c9 * Q.u;
    ctype const c14 = M.c0 * Q.w + M.c7 * Q.x + M.c8 * Q.y + M.c9 * Q.z + M.c10 * Q.w;
    ctype const c15 = M.c0 * Q.u - M.c1 * Q.x - M.c2 * Q.y - M.c3 * Q.z - M.c10 * Q.u;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(quadvec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c11 - Q.y * M.c12 - Q.z * M.c13 + Q.w * M.c15 + Q.u * M.c14;
    ctype const c1 = -Q.x * M.c15 + Q.y * M.c3 - Q.z * M.c2 + Q.u * M.c4 + Q.u * M.c11;
    ctype const c2 = -Q.x * M.c3 - Q.y * M.c15 + Q.z * M.c1 + Q.u * M.c5 + Q.u * M.c12;
    ctype const c3 = Q.x * M.c2 - Q.y * M.c1 - Q.z * M.c15 + Q.u * M.c6 + Q.u * M.c13;
    ctype const c4 = Q.x * M.c10 - Q.y * M.c13 + Q.z * M.c12 + Q.w * M.c1 + Q.u * M.c7;
    ctype const c5 = Q.x * M.c13 + Q.y * M.c10 - Q.z * M.c11 + Q.w * M.c2 + Q.u * M.c8;
    ctype const c6 = -Q.x * M.c12 + Q.y * M.c11 + Q.z * M.c10 + Q.w * M.c3 + Q.u * M.c9;
    ctype const c7 = Q.x * M.c14 - Q.y * M.c9 + Q.z * M.c8 + Q.w * M.c4 - Q.w * M.c11;
    ctype const c8 = Q.x * M.c9 + Q.y * M.c14 - Q.z * M.c7 + Q.w * M.c5 - Q.w * M.c12;
    ctype const c9 = -Q.x * M.c8 + Q.y * M.c7 + Q.z * M.c14 + Q.w * M.c6 - Q.w * M.c13;
    ctype const c10 = -Q.x * M.c4 - Q.y * M.c5 - Q.z * M.c6 + Q.w * M.c15 - Q.u * M.c14;
    ctype const c11 = Q.x * M.c0 - Q.y * M.c6 + Q.z * M.c5 + Q.w * M.c1 - Q.u * M.c7;
    ctype const c12 = Q.x * M.c6 + Q.y * M.c0 - Q.z * M.c4 + Q.w * M.c2 - Q.u * M.c8;
    ctype const c13 = -Q.x * M.c5 + Q.y * M.c4 + Q.z * M.c0 + Q.w * M.c3 - Q.u * M.c9;
    ctype const c14 = -Q.x * M.c7 - Q.y * M.c8 - Q.z * M.c9 + Q.w * M.c0 - Q.w * M.c10;
    ctype const c15 = Q.x * M.c1 + Q.y * M.c2 + Q.z * M.c3 + Q.u * M.c0 + Q.u * M.c10;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& M,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c2 * t.mz + M.c3 * t.my + M.c4 * t.pw + M.c8 * t.pz -
                     M.c9 * t.py - M.c10 * t.vx + M.c12 * t.vz - M.c13 * t.vy -
                     M.c14 * t.px + M.c15 * t.mx;
    ctype const c1 = M.c1 * t.mz - M.c3 * t.mx + M.c5 * t.pw - M.c7 * t.pz + M.c9 * t.px -
                     M.c10 * t.vy - M.c11 * t.vz + M.c13 * t.vx - M.c14 * t.py +
                     M.c15 * t.my;
    ctype const c2 = -M.c1 * t.my + M.c2 * t.mx + M.c6 * t.pw + M.c7 * t.py -
                     M.c8 * t.px - M.c10 * t.vz + M.c11 * t.vy - M.c12 * t.vx -
                     M.c14 * t.pz + M.c15 * t.mz;
    ctype const c3 = -M.c1 * t.vx - M.c2 * t.vy - M.c3 * t.vz - M.c4 * t.px -
                     M.c5 * t.py - M.c6 * t.pz + M.c11 * t.px + M.c12 * t.py +
                     M.c13 * t.pz - M.c15 * t.pw;
    ctype const c4 = -M.c4 * t.mx - M.c5 * t.my - M.c6 * t.mz - M.c7 * t.vx -
                     M.c8 * t.vy - M.c9 * t.vz - M.c11 * t.mx - M.c12 * t.my -
                     M.c13 * t.mz + M.c14 * t.pw;
    ctype const c5 = M.c0 * t.vx - M.c2 * t.mz + M.c3 * t.my - M.c5 * t.vz + M.c6 * t.vy -
                     M.c8 * t.pz + M.c9 * t.py - M.c11 * t.pw + M.c14 * t.px +
                     M.c15 * t.mx;
    ctype const c6 = M.c0 * t.vy + M.c1 * t.mz - M.c3 * t.mx + M.c4 * t.vz - M.c6 * t.vx +
                     M.c7 * t.pz - M.c9 * t.px - M.c12 * t.pw + M.c14 * t.py +
                     M.c15 * t.my;
    ctype const c7 = M.c0 * t.vz - M.c1 * t.my + M.c2 * t.mx - M.c4 * t.vy + M.c5 * t.vx -
                     M.c7 * t.py + M.c8 * t.px - M.c13 * t.pw + M.c14 * t.pz +
                     M.c15 * t.mz;
    ctype const c8 = M.c0 * t.mx - M.c5 * t.mz + M.c6 * t.my + M.c7 * t.pw - M.c8 * t.vz +
                     M.c9 * t.vy + M.c10 * t.mx - M.c12 * t.mz + M.c13 * t.my +
                     M.c14 * t.vx;
    ctype const c9 = M.c0 * t.my + M.c4 * t.mz - M.c6 * t.mx + M.c7 * t.vz + M.c8 * t.pw -
                     M.c9 * t.vx + M.c10 * t.my + M.c11 * t.mz - M.c13 * t.mx +
                     M.c14 * t.vy;
    ctype const c10 = M.c0 * t.mz - M.c4 * t.my + M.c5 * t.mx - M.c7 * t.vy +
                      M.c8 * t.vx + M.c9 * t.pw + M.c10 * t.mz - M.c11 * t.my +
                      M.c12 * t.mx + M.c14 * t.vz;
    ctype const c11 = M.c0 * t.px - M.c1 * t.pw - M.c2 * t.vz + M.c3 * t.vy -
                      M.c5 * t.pz + M.c6 * t.py - M.c10 * t.px + M.c12 * t.pz -
                      M.c13 * t.py + M.c15 * t.vx;
    ctype const c12 = M.c0 * t.py + M.c1 * t.vz - M.c2 * t.pw - M.c3 * t.vx +
                      M.c4 * t.pz - M.c6 * t.px - M.c10 * t.py - M.c11 * t.pz +
                      M.c13 * t.px + M.c15 * t.vy;
    ctype const c13 = M.c0 * t.pz - M.c1 * t.vy + M.c2 * t.vx - M.c3 * t.pw -
                      M.c4 * t.py + M.c5 * t.px - M.c10 * t.pz + M.c11 * t.py -
                      M.c12 * t.px + M.c15 * t.vz;
    ctype const c14 = M.c0 * t.pw - M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz +
                      M.c7 * t.px + M.c8 * t.py + M.c9 * t.pz + M.c11 * t.vx +
                      M.c12 * t.vy + M.c13 * t.vz;
    ctype const c15 = -M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz - M.c4 * t.vx -
                      M.c5 * t.vy - M.c6 * t.vz - M.c7 * t.px - M.c8 * t.py -
                      M.c9 * t.pz - M.c10 * t.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c10 + t.vy * M.c13 - t.vz * M.c12 - t.mx * M.c15 +
                     t.my * M.c3 - t.mz * M.c2 + t.px * M.c14 - t.py * M.c9 +
                     t.pz * M.c8 + t.pw * M.c4;
    ctype const c1 = -t.vx * M.c13 - t.vy * M.c10 + t.vz * M.c11 - t.mx * M.c3 -
                     t.my * M.c15 + t.mz * M.c1 + t.px * M.c9 + t.py * M.c14 -
                     t.pz * M.c7 + t.pw * M.c5;
    ctype const c2 = t.vx * M.c12 - t.vy * M.c11 - t.vz * M.c10 + t.mx * M.c2 -
                     t.my * M.c1 - t.mz * M.c15 - t.px * M.c8 + t.py * M.c7 +
                     t.pz * M.c14 + t.pw * M.c6;
    ctype const c3 = -t.vx * M.c1 - t.vy * M.c2 - t.vz * M.c3 - t.px * M.c4 -
                     t.px * M.c11 - t.py * M.c5 - t.py * M.c12 - t.pz * M.c6 -
                     t.pz * M.c13 + t.pw * M.c15;
    ctype const c4 = -t.vx * M.c7 - t.vy * M.c8 - t.vz * M.c9 - t.mx * M.c4 +
                     t.mx * M.c11 - t.my * M.c5 + t.my * M.c12 - t.mz * M.c6 +
                     t.mz * M.c13 - t.pw * M.c14;
    ctype const c5 = t.vx * M.c0 - t.vy * M.c6 + t.vz * M.c5 + t.mx * M.c15 -
                     t.my * M.c3 + t.mz * M.c2 + t.px * M.c14 - t.py * M.c9 +
                     t.pz * M.c8 - t.pw * M.c11;
    ctype const c6 = t.vx * M.c6 + t.vy * M.c0 - t.vz * M.c4 + t.mx * M.c3 +
                     t.my * M.c15 - t.mz * M.c1 + t.px * M.c9 + t.py * M.c14 -
                     t.pz * M.c7 - t.pw * M.c12;
    ctype const c7 = -t.vx * M.c5 + t.vy * M.c4 + t.vz * M.c0 - t.mx * M.c2 +
                     t.my * M.c1 + t.mz * M.c15 - t.px * M.c8 + t.py * M.c7 +
                     t.pz * M.c14 - t.pw * M.c13;
    ctype const c8 = t.vx * M.c14 - t.vy * M.c9 + t.vz * M.c8 + t.mx * M.c0 -
                     t.mx * M.c10 - t.my * M.c6 + t.my * M.c13 + t.mz * M.c5 -
                     t.mz * M.c12 - t.pw * M.c7;
    ctype const c9 = t.vx * M.c9 + t.vy * M.c14 - t.vz * M.c7 + t.mx * M.c6 -
                     t.mx * M.c13 + t.my * M.c0 - t.my * M.c10 - t.mz * M.c4 +
                     t.mz * M.c11 - t.pw * M.c8;
    ctype const c10 = -t.vx * M.c8 + t.vy * M.c7 + t.vz * M.c14 - t.mx * M.c5 +
                      t.mx * M.c12 + t.my * M.c4 - t.my * M.c11 + t.mz * M.c0 -
                      t.mz * M.c10 - t.pw * M.c9;
    ctype const c11 = t.vx * M.c15 - t.vy * M.c3 + t.vz * M.c2 + t.px * M.c0 +
                      t.px * M.c10 - t.py * M.c6 - t.py * M.c13 + t.pz * M.c5 +
                      t.pz * M.c12 + t.pw * M.c1;
    ctype const c12 = t.vx * M.c3 + t.vy * M.c15 - t.vz * M.c1 + t.px * M.c6 +
                      t.px * M.c13 + t.py * M.c0 + t.py * M.c10 - t.pz * M.c4 -
                      t.pz * M.c11 + t.pw * M.c2;
    ctype const c13 = -t.vx * M.c2 + t.vy * M.c1 + t.vz * M.c15 - t.px * M.c5 -
                      t.px * M.c12 + t.py * M.c4 + t.py * M.c11 + t.pz * M.c0 +
                      t.pz * M.c10 + t.pw * M.c3;
    ctype const c14 = t.vx * M.c11 + t.vy * M.c12 + t.vz * M.c13 + t.mx * M.c1 +
                      t.my * M.c2 + t.mz * M.c3 - t.px * M.c7 - t.py * M.c8 -
                      t.pz * M.c9 + t.pw * M.c0;
    ctype const c15 = -t.vx * M.c4 - t.vy * M.c5 - t.vz * M.c6 - t.mx * M.c1 -
                      t.my * M.c2 - t.mz * M.c3 - t.px * M.c7 - t.py * M.c8 -
                      t.pz * M.c9 - t.pw * M.c10;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& M,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.px - M.c2 * B.py - M.c3 * B.pz - M.c4 * B.mx -
                     M.c5 * B.my - M.c6 * B.mz - M.c7 * B.vx - M.c8 * B.vy - M.c9 * B.vz +
                     M.c10 * B.pw;
    ctype const c1 = M.c0 * B.vx + M.c1 * B.pw - M.c2 * B.mz + M.c3 * B.my - M.c5 * B.vz +
                     M.c6 * B.vy - M.c10 * B.vx + M.c12 * B.vz - M.c13 * B.vy +
                     M.c15 * B.mx;
    ctype const c2 = M.c0 * B.vy + M.c1 * B.mz + M.c2 * B.pw - M.c3 * B.mx + M.c4 * B.vz -
                     M.c6 * B.vx - M.c10 * B.vy - M.c11 * B.vz + M.c13 * B.vx +
                     M.c15 * B.my;
    ctype const c3 = M.c0 * B.vz - M.c1 * B.my + M.c2 * B.mx + M.c3 * B.pw - M.c4 * B.vy +
                     M.c5 * B.vx - M.c10 * B.vz + M.c11 * B.vy - M.c12 * B.vx +
                     M.c15 * B.mz;
    ctype const c4 = M.c0 * B.mx - M.c2 * B.pz + M.c3 * B.py - M.c5 * B.mz + M.c6 * B.my -
                     M.c8 * B.vz + M.c9 * B.vy + M.c11 * B.pw + M.c14 * B.vx +
                     M.c15 * B.px;
    ctype const c5 = M.c0 * B.my + M.c1 * B.pz - M.c3 * B.px + M.c4 * B.mz - M.c6 * B.mx +
                     M.c7 * B.vz - M.c9 * B.vx + M.c12 * B.pw + M.c14 * B.vy +
                     M.c15 * B.py;
    ctype const c6 = M.c0 * B.mz - M.c1 * B.py + M.c2 * B.px - M.c4 * B.my + M.c5 * B.mx -
                     M.c7 * B.vy + M.c8 * B.vx + M.c13 * B.pw + M.c14 * B.vz +
                     M.c15 * B.pz;
    ctype const c7 = M.c0 * B.px - M.c5 * B.pz + M.c6 * B.py - M.c7 * B.pw - M.c8 * B.mz +
                     M.c9 * B.my + M.c10 * B.px - M.c12 * B.pz + M.c13 * B.py +
                     M.c14 * B.mx;
    ctype const c8 = M.c0 * B.py + M.c4 * B.pz - M.c6 * B.px + M.c7 * B.mz - M.c8 * B.pw -
                     M.c9 * B.mx + M.c10 * B.py + M.c11 * B.pz - M.c13 * B.px +
                     M.c14 * B.my;
    ctype const c9 = M.c0 * B.pz - M.c4 * B.py + M.c5 * B.px - M.c7 * B.my + M.c8 * B.mx -
                     M.c9 * B.pw + M.c10 * B.pz - M.c11 * B.py + M.c12 * B.px +
                     M.c14 * B.mz;
    ctype const c10 = M.c0 * B.pw + M.c1 * B.px + M.c2 * B.py + M.c3 * B.pz -
                      M.c7 * B.vx - M.c8 * B.vy - M.c9 * B.vz - M.c11 * B.mx -
                      M.c12 * B.my - M.c13 * B.mz;
    ctype const c11 = M.c2 * B.pz - M.c3 * B.py + M.c4 * B.pw - M.c8 * B.vz +
                      M.c9 * B.vy + M.c10 * B.mx - M.c12 * B.mz + M.c13 * B.my +
                      M.c14 * B.vx - M.c15 * B.px;
    ctype const c12 = -M.c1 * B.pz + M.c3 * B.px + M.c5 * B.pw + M.c7 * B.vz -
                      M.c9 * B.vx + M.c10 * B.my + M.c11 * B.mz - M.c13 * B.mx +
                      M.c14 * B.vy - M.c15 * B.py;
    ctype const c13 = M.c1 * B.py - M.c2 * B.px + M.c6 * B.pw - M.c7 * B.vy +
                      M.c8 * B.vx + M.c10 * B.mz - M.c11 * B.my + M.c12 * B.mx +
                      M.c14 * B.vz - M.c15 * B.pz;
    ctype const c14 = -M.c4 * B.px - M.c5 * B.py - M.c6 * B.pz - M.c7 * B.mx -
                      M.c8 * B.my - M.c9 * B.mz - M.c11 * B.px - M.c12 * B.py -
                      M.c13 * B.pz - M.c14 * B.pw;
    ctype const c15 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz - M.c4 * B.vx -
                      M.c5 * B.vy - M.c6 * B.vz + M.c11 * B.vx + M.c12 * B.vy +
                      M.c13 * B.vz + M.c15 * B.pw;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c7 - B.vy * M.c8 - B.vz * M.c9 - B.mx * M.c4 -
                     B.my * M.c5 - B.mz * M.c6 - B.px * M.c1 - B.py * M.c2 - B.pz * M.c3 +
                     B.pw * M.c10;
    ctype const c1 = B.vx * M.c0 + B.vx * M.c10 - B.vy * M.c6 - B.vy * M.c13 +
                     B.vz * M.c5 + B.vz * M.c12 + B.mx * M.c15 - B.my * M.c3 +
                     B.mz * M.c2 - B.pw * M.c1;
    ctype const c2 = B.vx * M.c6 + B.vx * M.c13 + B.vy * M.c0 + B.vy * M.c10 -
                     B.vz * M.c4 - B.vz * M.c11 + B.mx * M.c3 + B.my * M.c15 -
                     B.mz * M.c1 - B.pw * M.c2;
    ctype const c3 = -B.vx * M.c5 - B.vx * M.c12 + B.vy * M.c4 + B.vy * M.c11 +
                     B.vz * M.c0 + B.vz * M.c10 - B.mx * M.c2 + B.my * M.c1 +
                     B.mz * M.c15 - B.pw * M.c3;
    ctype const c4 = B.vx * M.c14 - B.vy * M.c9 + B.vz * M.c8 + B.mx * M.c0 -
                     B.my * M.c6 + B.mz * M.c5 + B.px * M.c15 - B.py * M.c3 +
                     B.pz * M.c2 + B.pw * M.c11;
    ctype const c5 = B.vx * M.c9 + B.vy * M.c14 - B.vz * M.c7 + B.mx * M.c6 +
                     B.my * M.c0 - B.mz * M.c4 + B.px * M.c3 + B.py * M.c15 -
                     B.pz * M.c1 + B.pw * M.c12;
    ctype const c6 = -B.vx * M.c8 + B.vy * M.c7 + B.vz * M.c14 - B.mx * M.c5 +
                     B.my * M.c4 + B.mz * M.c0 - B.px * M.c2 + B.py * M.c1 +
                     B.pz * M.c15 + B.pw * M.c13;
    ctype const c7 = B.mx * M.c14 - B.my * M.c9 + B.mz * M.c8 + B.px * M.c0 -
                     B.px * M.c10 - B.py * M.c6 + B.py * M.c13 + B.pz * M.c5 -
                     B.pz * M.c12 + B.pw * M.c7;
    ctype const c8 = B.mx * M.c9 + B.my * M.c14 - B.mz * M.c7 + B.px * M.c6 -
                     B.px * M.c13 + B.py * M.c0 - B.py * M.c10 - B.pz * M.c4 +
                     B.pz * M.c11 + B.pw * M.c8;
    ctype const c9 = -B.mx * M.c8 + B.my * M.c7 + B.mz * M.c14 - B.px * M.c5 +
                     B.px * M.c12 + B.py * M.c4 - B.py * M.c11 + B.pz * M.c0 -
                     B.pz * M.c10 + B.pw * M.c9;
    ctype const c10 = B.vx * M.c7 + B.vy * M.c8 + B.vz * M.c9 - B.mx * M.c11 -
                      B.my * M.c12 - B.mz * M.c13 - B.px * M.c1 - B.py * M.c2 -
                      B.pz * M.c3 + B.pw * M.c0;
    ctype const c11 = -B.vx * M.c14 + B.vy * M.c9 - B.vz * M.c8 + B.mx * M.c10 -
                      B.my * M.c13 + B.mz * M.c12 + B.px * M.c15 - B.py * M.c3 +
                      B.pz * M.c2 + B.pw * M.c4;
    ctype const c12 = -B.vx * M.c9 - B.vy * M.c14 + B.vz * M.c7 + B.mx * M.c13 +
                      B.my * M.c10 - B.mz * M.c11 + B.px * M.c3 + B.py * M.c15 -
                      B.pz * M.c1 + B.pw * M.c5;
    ctype const c13 = B.vx * M.c8 - B.vy * M.c7 - B.vz * M.c14 - B.mx * M.c12 +
                      B.my * M.c11 + B.mz * M.c10 - B.px * M.c2 + B.py * M.c1 +
                      B.pz * M.c15 + B.pw * M.c6;
    ctype const c14 = -B.mx * M.c7 - B.my * M.c8 - B.mz * M.c9 - B.px * M.c4 +
                      B.px * M.c11 - B.py * M.c5 + B.py * M.c12 - B.pz * M.c6 +
                      B.pz * M.c13 + B.pw * M.c14;
    ctype const c15 = -B.vx * M.c4 - B.vx * M.c11 - B.vy * M.c5 - B.vy * M.c12 -
                      B.vz * M.c6 - B.vz * M.c13 - B.mx * M.c1 - B.my * M.c2 -
                      B.mz * M.c3 - B.pw * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& M,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.u - M.c5 * v.z + M.c6 * v.y - M.c7 * v.w;
    ctype const c1 = M.c0 * v.y + M.c2 * v.u + M.c4 * v.z - M.c6 * v.x - M.c8 * v.w;
    ctype const c2 = M.c0 * v.z + M.c3 * v.u - M.c4 * v.y + M.c5 * v.x - M.c9 * v.w;
    ctype const c3 = M.c0 * v.w + M.c1 * v.x + M.c2 * v.y + M.c3 * v.z - M.c10 * v.w;
    ctype const c4 = M.c0 * v.u - M.c7 * v.x - M.c8 * v.y - M.c9 * v.z + M.c10 * v.u;
    ctype const c5 = M.c1 * v.u + M.c7 * v.w - M.c10 * v.x + M.c12 * v.z - M.c13 * v.y;
    ctype const c6 = M.c2 * v.u + M.c8 * v.w - M.c10 * v.y - M.c11 * v.z + M.c13 * v.x;
    ctype const c7 = M.c3 * v.u + M.c9 * v.w - M.c10 * v.z + M.c11 * v.y - M.c12 * v.x;
    ctype const c8 = M.c4 * v.u - M.c8 * v.z + M.c9 * v.y + M.c11 * v.u + M.c14 * v.x;
    ctype const c9 = M.c5 * v.u + M.c7 * v.z - M.c9 * v.x + M.c12 * v.u + M.c14 * v.y;
    ctype const c10 = M.c6 * v.u - M.c7 * v.y + M.c8 * v.x + M.c13 * v.u + M.c14 * v.z;
    ctype const c11 = M.c2 * v.z - M.c3 * v.y + M.c4 * v.w - M.c11 * v.w - M.c15 * v.x;
    ctype const c12 = -M.c1 * v.z + M.c3 * v.x + M.c5 * v.w - M.c12 * v.w - M.c15 * v.y;
    ctype const c13 = M.c1 * v.y - M.c2 * v.x + M.c6 * v.w - M.c13 * v.w - M.c15 * v.z;
    ctype const c14 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z - M.c14 * v.w + M.c15 * v.u;
    ctype const c15 = M.c11 * v.x + M.c12 * v.y + M.c13 * v.z + M.c14 * v.w + M.c15 * v.u;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 - v.y * M.c6 + v.z * M.c5 + v.w * M.c7 - v.u * M.c1;
    ctype const c1 = v.x * M.c6 + v.y * M.c0 - v.z * M.c4 + v.w * M.c8 - v.u * M.c2;
    ctype const c2 = -v.x * M.c5 + v.y * M.c4 + v.z * M.c0 + v.w * M.c9 - v.u * M.c3;
    ctype const c3 = -v.x * M.c1 - v.y * M.c2 - v.z * M.c3 + v.w * M.c0 + v.w * M.c10;
    ctype const c4 = v.x * M.c7 + v.y * M.c8 + v.z * M.c9 + v.u * M.c0 - v.u * M.c10;
    ctype const c5 = -v.x * M.c10 + v.y * M.c13 - v.z * M.c12 + v.w * M.c7 + v.u * M.c1;
    ctype const c6 = -v.x * M.c13 - v.y * M.c10 + v.z * M.c11 + v.w * M.c8 + v.u * M.c2;
    ctype const c7 = v.x * M.c12 - v.y * M.c11 - v.z * M.c10 + v.w * M.c9 + v.u * M.c3;
    ctype const c8 = -v.x * M.c14 + v.y * M.c9 - v.z * M.c8 + v.u * M.c4 - v.u * M.c11;
    ctype const c9 = -v.x * M.c9 - v.y * M.c14 + v.z * M.c7 + v.u * M.c5 - v.u * M.c12;
    ctype const c10 = v.x * M.c8 - v.y * M.c7 - v.z * M.c14 + v.u * M.c6 - v.u * M.c13;
    ctype const c11 = v.x * M.c15 - v.y * M.c3 + v.z * M.c2 + v.w * M.c4 + v.w * M.c11;
    ctype const c12 = v.x * M.c3 + v.y * M.c15 - v.z * M.c1 + v.w * M.c5 + v.w * M.c12;
    ctype const c13 = -v.x * M.c2 + v.y * M.c1 + v.z * M.c15 + v.w * M.c6 + v.w * M.c13;
    ctype const c14 = -v.x * M.c4 - v.y * M.c5 - v.z * M.c6 + v.w * M.c14 - v.u * M.c15;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14 + v.u * M.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_E<T> const& M,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(Scalar3dc<T> s,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& A,
                                                        MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 - A.c4 * B.c3 +
                     A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7 + A.c8 * B.c11 +
                     A.c9 * B.c12 + A.c10 * B.c13 + A.c11 * B.c8 + A.c12 * B.c9 +
                     A.c13 * B.c10 - A.c14 * B.c14 - A.c15 * B.c15;
    ctype const c1 = -A.c0 * B.c3 + A.c1 * B.c13 - A.c2 * B.c12 + A.c3 * B.c0 -
                     A.c3 * B.c5 - A.c5 * B.c3 + A.c6 * B.c13 - A.c7 * B.c12 +
                     A.c11 * B.c14 - A.c11 * B.c15 - A.c12 * B.c2 + A.c12 * B.c7 +
                     A.c13 * B.c1 - A.c13 * B.c6 - A.c14 * B.c11 - A.c15 * B.c11;
    ctype const c2 = -A.c0 * B.c13 - A.c1 * B.c3 + A.c2 * B.c11 + A.c3 * B.c1 -
                     A.c3 * B.c6 - A.c5 * B.c13 - A.c6 * B.c3 + A.c7 * B.c11 +
                     A.c11 * B.c2 - A.c11 * B.c7 + A.c12 * B.c14 - A.c12 * B.c15 -
                     A.c13 * B.c0 + A.c13 * B.c5 - A.c14 * B.c12 - A.c15 * B.c12;
    ctype const c3 = A.c0 * B.c12 - A.c1 * B.c11 - A.c2 * B.c3 + A.c3 * B.c2 -
                     A.c3 * B.c7 + A.c5 * B.c12 - A.c6 * B.c11 - A.c7 * B.c3 -
                     A.c11 * B.c1 + A.c11 * B.c6 + A.c12 * B.c0 - A.c12 * B.c5 +
                     A.c13 * B.c14 - A.c13 * B.c15 - A.c14 * B.c13 - A.c15 * B.c13;
    ctype const c4 = -A.c0 * B.c14 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c8 -
                     A.c4 * B.c11 - A.c5 * B.c15 + A.c6 * B.c7 - A.c7 * B.c6 -
                     A.c8 * B.c3 + A.c9 * B.c13 - A.c10 * B.c12 - A.c11 * B.c4 +
                     A.c12 * B.c10 - A.c13 * B.c9 - A.c14 * B.c0 - A.c15 * B.c5;
    ctype const c5 = -A.c0 * B.c2 - A.c1 * B.c14 + A.c2 * B.c0 - A.c3 * B.c9 -
                     A.c4 * B.c12 - A.c5 * B.c7 - A.c6 * B.c15 + A.c7 * B.c5 -
                     A.c8 * B.c13 - A.c9 * B.c3 + A.c10 * B.c11 - A.c11 * B.c10 -
                     A.c12 * B.c4 + A.c13 * B.c8 - A.c14 * B.c1 - A.c15 * B.c6;
    ctype const c6 = A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c14 - A.c3 * B.c10 -
                     A.c4 * B.c13 + A.c5 * B.c6 - A.c6 * B.c5 - A.c7 * B.c15 +
                     A.c8 * B.c12 - A.c9 * B.c11 - A.c10 * B.c3 + A.c11 * B.c9 -
                     A.c12 * B.c8 - A.c13 * B.c4 - A.c14 * B.c2 - A.c15 * B.c7;
    ctype const c7 = A.c0 * B.c4 - A.c1 * B.c10 + A.c2 * B.c9 - A.c4 * B.c0 -
                     A.c4 * B.c5 - A.c5 * B.c4 + A.c6 * B.c10 - A.c7 * B.c9 -
                     A.c8 * B.c14 - A.c8 * B.c15 + A.c9 * B.c2 + A.c9 * B.c7 -
                     A.c10 * B.c1 - A.c10 * B.c6 + A.c14 * B.c8 - A.c15 * B.c8;
    ctype const c8 = A.c0 * B.c10 + A.c1 * B.c4 - A.c2 * B.c8 - A.c4 * B.c1 -
                     A.c4 * B.c6 - A.c5 * B.c10 - A.c6 * B.c4 + A.c7 * B.c8 -
                     A.c8 * B.c2 - A.c8 * B.c7 - A.c9 * B.c14 - A.c9 * B.c15 +
                     A.c10 * B.c0 + A.c10 * B.c5 + A.c14 * B.c9 - A.c15 * B.c9;
    ctype const c9 = -A.c0 * B.c9 + A.c1 * B.c8 + A.c2 * B.c4 - A.c4 * B.c2 -
                     A.c4 * B.c7 + A.c5 * B.c9 - A.c6 * B.c8 - A.c7 * B.c4 + A.c8 * B.c1 +
                     A.c8 * B.c6 - A.c9 * B.c0 - A.c9 * B.c5 - A.c10 * B.c14 -
                     A.c10 * B.c15 + A.c14 * B.c10 - A.c15 * B.c10;
    ctype const c10 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c7 + A.c3 * B.c4 -
                      A.c4 * B.c3 - A.c5 * B.c0 - A.c6 * B.c1 - A.c7 * B.c2 +
                      A.c8 * B.c11 + A.c9 * B.c12 + A.c10 * B.c13 - A.c11 * B.c8 -
                      A.c12 * B.c9 - A.c13 * B.c10 + A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c11 = A.c0 * B.c15 - A.c1 * B.c7 + A.c2 * B.c6 + A.c3 * B.c8 -
                      A.c4 * B.c11 + A.c5 * B.c14 - A.c6 * B.c2 + A.c7 * B.c1 -
                      A.c8 * B.c3 + A.c9 * B.c13 - A.c10 * B.c12 + A.c11 * B.c4 -
                      A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 + A.c15 * B.c0;
    ctype const c12 = A.c0 * B.c7 + A.c1 * B.c15 - A.c2 * B.c5 + A.c3 * B.c9 -
                      A.c4 * B.c12 + A.c5 * B.c2 + A.c6 * B.c14 - A.c7 * B.c0 -
                      A.c8 * B.c13 - A.c9 * B.c3 + A.c10 * B.c11 + A.c11 * B.c10 +
                      A.c12 * B.c4 - A.c13 * B.c8 + A.c14 * B.c6 + A.c15 * B.c1;
    ctype const c13 = -A.c0 * B.c6 + A.c1 * B.c5 + A.c2 * B.c15 + A.c3 * B.c10 -
                      A.c4 * B.c13 - A.c5 * B.c1 + A.c6 * B.c0 + A.c7 * B.c14 +
                      A.c8 * B.c12 - A.c9 * B.c11 - A.c10 * B.c3 - A.c11 * B.c9 +
                      A.c12 * B.c8 + A.c13 * B.c4 + A.c14 * B.c7 + A.c15 * B.c2;
    ctype const c14 = -A.c0 * B.c8 - A.c1 * B.c9 - A.c2 * B.c10 - A.c4 * B.c14 -
                      A.c4 * B.c15 + A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c10 +
                      A.c8 * B.c0 + A.c8 * B.c5 + A.c9 * B.c1 + A.c9 * B.c6 +
                      A.c10 * B.c2 + A.c10 * B.c7 + A.c14 * B.c4 - A.c15 * B.c4;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 -
                      A.c3 * B.c15 + A.c5 * B.c11 + A.c6 * B.c12 + A.c7 * B.c13 -
                      A.c11 * B.c0 + A.c11 * B.c5 - A.c12 * B.c1 + A.c12 * B.c6 -
                      A.c13 * B.c2 + A.c13 * B.c7 - A.c14 * B.c3 - A.c15 * B.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_u,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& M,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c15 * ctype(ps);
    ctype const c1 = -M.c11 * ctype(ps);
    ctype const c2 = -M.c12 * ctype(ps);
    ctype const c3 = -M.c13 * ctype(ps);
    ctype const c4 = -M.c5 * ctype(ps);
    ctype const c5 = -M.c6 * ctype(ps);
    ctype const c6 = -M.c7 * ctype(ps);
    ctype const c7 = -M.c8 * ctype(ps);
    ctype const c8 = -M.c9 * ctype(ps);
    ctype const c9 = -M.c10 * ctype(ps);
    ctype const c10 = M.c14 * ctype(ps);
    ctype const c11 = M.c0 * ctype(ps);
    ctype const c12 = M.c1 * ctype(ps);
    ctype const c13 = M.c2 * ctype(ps);
    ctype const c14 = -M.c4 * ctype(ps);
    ctype const c15 = -M.c3 * ctype(ps);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(ps,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c15;
    ctype const c1 = -ctype(ps) * M.c11;
    ctype const c2 = -ctype(ps) * M.c12;
    ctype const c3 = -ctype(ps) * M.c13;
    ctype const c4 = -ctype(ps) * M.c5;
    ctype const c5 = -ctype(ps) * M.c6;
    ctype const c6 = -ctype(ps) * M.c7;
    ctype const c7 = -ctype(ps) * M.c8;
    ctype const c8 = -ctype(ps) * M.c9;
    ctype const c9 = -ctype(ps) * M.c10;
    ctype const c10 = ctype(ps) * M.c14;
    ctype const c11 = ctype(ps) * M.c0;
    ctype const c12 = ctype(ps) * M.c1;
    ctype const c13 = ctype(ps) * M.c2;
    ctype const c14 = -ctype(ps) * M.c4;
    ctype const c15 = -ctype(ps) * M.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_u,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& M,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c6 * Q.z - M.c7 * Q.y - M.c8 * Q.u + M.c11 * Q.w - M.c15 * Q.x;
    ctype const c1 = -M.c5 * Q.z + M.c7 * Q.x - M.c9 * Q.u + M.c12 * Q.w - M.c15 * Q.y;
    ctype const c2 = M.c5 * Q.y - M.c6 * Q.x - M.c10 * Q.u + M.c13 * Q.w - M.c15 * Q.z;
    ctype const c3 = -M.c11 * Q.x - M.c12 * Q.y - M.c13 * Q.z + M.c14 * Q.u + M.c15 * Q.u;
    ctype const c4 = M.c8 * Q.x + M.c9 * Q.y + M.c10 * Q.z - M.c14 * Q.w + M.c15 * Q.w;
    ctype const c5 = M.c1 * Q.z - M.c2 * Q.y + M.c8 * Q.u + M.c11 * Q.w - M.c14 * Q.x;
    ctype const c6 = -M.c0 * Q.z + M.c2 * Q.x + M.c9 * Q.u + M.c12 * Q.w - M.c14 * Q.y;
    ctype const c7 = M.c0 * Q.y - M.c1 * Q.x + M.c10 * Q.u + M.c13 * Q.w - M.c14 * Q.z;
    ctype const c8 = -M.c0 * Q.w - M.c4 * Q.x + M.c5 * Q.w + M.c9 * Q.z - M.c10 * Q.y;
    ctype const c9 = -M.c1 * Q.w - M.c4 * Q.y + M.c6 * Q.w - M.c8 * Q.z + M.c10 * Q.x;
    ctype const c10 = -M.c2 * Q.w - M.c4 * Q.z + M.c7 * Q.w + M.c8 * Q.y - M.c9 * Q.x;
    ctype const c11 = M.c0 * Q.u + M.c3 * Q.x + M.c5 * Q.u - M.c12 * Q.z + M.c13 * Q.y;
    ctype const c12 = M.c1 * Q.u + M.c3 * Q.y + M.c6 * Q.u + M.c11 * Q.z - M.c13 * Q.x;
    ctype const c13 = M.c2 * Q.u + M.c3 * Q.z + M.c7 * Q.u - M.c11 * Q.y + M.c12 * Q.x;
    ctype const c14 = M.c3 * Q.w - M.c4 * Q.u + M.c5 * Q.x + M.c6 * Q.y + M.c7 * Q.z;
    ctype const c15 = M.c0 * Q.x + M.c1 * Q.y + M.c2 * Q.z + M.c3 * Q.w + M.c4 * Q.u;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(quadvec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * M.c15 + Q.y * M.c7 - Q.z * M.c6 - Q.w * M.c11 + Q.u * M.c8;
    ctype const c1 = -Q.x * M.c7 - Q.y * M.c15 + Q.z * M.c5 - Q.w * M.c12 + Q.u * M.c9;
    ctype const c2 = Q.x * M.c6 - Q.y * M.c5 - Q.z * M.c15 - Q.w * M.c13 + Q.u * M.c10;
    ctype const c3 = Q.x * M.c11 + Q.y * M.c12 + Q.z * M.c13 - Q.u * M.c14 + Q.u * M.c15;
    ctype const c4 = -Q.x * M.c8 - Q.y * M.c9 - Q.z * M.c10 + Q.w * M.c14 + Q.w * M.c15;
    ctype const c5 = -Q.x * M.c14 + Q.y * M.c2 - Q.z * M.c1 + Q.w * M.c11 + Q.u * M.c8;
    ctype const c6 = -Q.x * M.c2 - Q.y * M.c14 + Q.z * M.c0 + Q.w * M.c12 + Q.u * M.c9;
    ctype const c7 = Q.x * M.c1 - Q.y * M.c0 - Q.z * M.c14 + Q.w * M.c13 + Q.u * M.c10;
    ctype const c8 = Q.x * M.c4 - Q.y * M.c10 + Q.z * M.c9 + Q.w * M.c0 + Q.w * M.c5;
    ctype const c9 = Q.x * M.c10 + Q.y * M.c4 - Q.z * M.c8 + Q.w * M.c1 + Q.w * M.c6;
    ctype const c10 = -Q.x * M.c9 + Q.y * M.c8 + Q.z * M.c4 + Q.w * M.c2 + Q.w * M.c7;
    ctype const c11 = -Q.x * M.c3 + Q.y * M.c13 - Q.z * M.c12 - Q.u * M.c0 + Q.u * M.c5;
    ctype const c12 = -Q.x * M.c13 - Q.y * M.c3 + Q.z * M.c11 - Q.u * M.c1 + Q.u * M.c6;
    ctype const c13 = Q.x * M.c12 - Q.y * M.c11 - Q.z * M.c3 - Q.u * M.c2 + Q.u * M.c7;
    ctype const c14 = Q.x * M.c5 + Q.y * M.c6 + Q.z * M.c7 - Q.w * M.c3 + Q.u * M.c4;
    ctype const c15 = Q.x * M.c0 + Q.y * M.c1 + Q.z * M.c2 + Q.w * M.c3 + Q.u * M.c4;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& M,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * t.vx + M.c6 * t.vy + M.c7 * t.vz + M.c8 * t.px + M.c9 * t.py +
                     M.c10 * t.pz + M.c11 * t.mx + M.c12 * t.my + M.c13 * t.mz -
                     M.c14 * t.pw;
    ctype const c1 = M.c1 * t.pz - M.c2 * t.py - M.c3 * t.vx + M.c6 * t.pz - M.c7 * t.py +
                     M.c11 * t.pw + M.c12 * t.vz - M.c13 * t.vy - M.c14 * t.px -
                     M.c15 * t.px;
    ctype const c2 = -M.c0 * t.pz + M.c2 * t.px - M.c3 * t.vy - M.c5 * t.pz +
                     M.c7 * t.px - M.c11 * t.vz + M.c12 * t.pw + M.c13 * t.vx -
                     M.c14 * t.py - M.c15 * t.py;
    ctype const c3 = M.c0 * t.py - M.c1 * t.px - M.c3 * t.vz + M.c5 * t.py - M.c6 * t.px +
                     M.c11 * t.vy - M.c12 * t.vx + M.c13 * t.pw - M.c14 * t.pz -
                     M.c15 * t.pz;
    ctype const c4 = -M.c0 * t.pw - M.c3 * t.mx - M.c4 * t.px + M.c6 * t.vz -
                     M.c7 * t.vy + M.c9 * t.pz - M.c10 * t.py + M.c12 * t.mz -
                     M.c13 * t.my - M.c15 * t.vx;
    ctype const c5 = -M.c1 * t.pw - M.c3 * t.my - M.c4 * t.py - M.c5 * t.vz +
                     M.c7 * t.vx - M.c8 * t.pz + M.c10 * t.px - M.c11 * t.mz +
                     M.c13 * t.mx - M.c15 * t.vy;
    ctype const c6 = -M.c2 * t.pw - M.c3 * t.mz - M.c4 * t.pz + M.c5 * t.vy -
                     M.c6 * t.vx + M.c8 * t.py - M.c9 * t.px + M.c11 * t.my -
                     M.c12 * t.mx - M.c15 * t.vz;
    ctype const c7 = -M.c1 * t.mz + M.c2 * t.my - M.c4 * t.vx + M.c6 * t.mz -
                     M.c7 * t.my - M.c8 * t.pw + M.c9 * t.vz - M.c10 * t.vy +
                     M.c14 * t.mx - M.c15 * t.mx;
    ctype const c8 = M.c0 * t.mz - M.c2 * t.mx - M.c4 * t.vy - M.c5 * t.mz + M.c7 * t.mx -
                     M.c8 * t.vz - M.c9 * t.pw + M.c10 * t.vx + M.c14 * t.my -
                     M.c15 * t.my;
    ctype const c9 = -M.c0 * t.my + M.c1 * t.mx - M.c4 * t.vz + M.c5 * t.my -
                     M.c6 * t.mx + M.c8 * t.vy - M.c9 * t.vx - M.c10 * t.pw +
                     M.c14 * t.mz - M.c15 * t.mz;
    ctype const c10 = -M.c0 * t.vx - M.c1 * t.vy - M.c2 * t.vz + M.c8 * t.px +
                      M.c9 * t.py + M.c10 * t.pz - M.c11 * t.mx - M.c12 * t.my -
                      M.c13 * t.mz + M.c15 * t.pw;
    ctype const c11 = -M.c1 * t.vz + M.c2 * t.vy + M.c3 * t.mx - M.c4 * t.px +
                      M.c5 * t.pw + M.c9 * t.pz - M.c10 * t.py - M.c12 * t.mz +
                      M.c13 * t.my + M.c14 * t.vx;
    ctype const c12 = M.c0 * t.vz - M.c2 * t.vx + M.c3 * t.my - M.c4 * t.py +
                      M.c6 * t.pw - M.c8 * t.pz + M.c10 * t.px + M.c11 * t.mz -
                      M.c13 * t.mx + M.c14 * t.vy;
    ctype const c13 = -M.c0 * t.vy + M.c1 * t.vx + M.c3 * t.mz - M.c4 * t.pz +
                      M.c7 * t.pw + M.c8 * t.py - M.c9 * t.px - M.c11 * t.my +
                      M.c12 * t.mx + M.c14 * t.vz;
    ctype const c14 = -M.c0 * t.mx - M.c1 * t.my - M.c2 * t.mz - M.c4 * t.pw +
                      M.c5 * t.mx + M.c6 * t.my + M.c7 * t.mz + M.c8 * t.vx +
                      M.c9 * t.vy + M.c10 * t.vz;
    ctype const c15 = M.c0 * t.px + M.c1 * t.py + M.c2 * t.pz + M.c3 * t.pw +
                      M.c5 * t.px + M.c6 * t.py + M.c7 * t.pz + M.c11 * t.vx +
                      M.c12 * t.vy + M.c13 * t.vz;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * M.c5 + t.vy * M.c6 + t.vz * M.c7 + t.mx * M.c11 +
                     t.my * M.c12 + t.mz * M.c13 + t.px * M.c8 + t.py * M.c9 +
                     t.pz * M.c10 - t.pw * M.c14;
    ctype const c1 = -t.vx * M.c3 + t.vy * M.c13 - t.vz * M.c12 + t.px * M.c14 -
                     t.px * M.c15 - t.py * M.c2 + t.py * M.c7 + t.pz * M.c1 -
                     t.pz * M.c6 - t.pw * M.c11;
    ctype const c2 = -t.vx * M.c13 - t.vy * M.c3 + t.vz * M.c11 + t.px * M.c2 -
                     t.px * M.c7 + t.py * M.c14 - t.py * M.c15 - t.pz * M.c0 +
                     t.pz * M.c5 - t.pw * M.c12;
    ctype const c3 = t.vx * M.c12 - t.vy * M.c11 - t.vz * M.c3 - t.px * M.c1 +
                     t.px * M.c6 + t.py * M.c0 - t.py * M.c5 + t.pz * M.c14 -
                     t.pz * M.c15 - t.pw * M.c13;
    ctype const c4 = -t.vx * M.c15 + t.vy * M.c7 - t.vz * M.c6 - t.mx * M.c3 +
                     t.my * M.c13 - t.mz * M.c12 - t.px * M.c4 + t.py * M.c10 -
                     t.pz * M.c9 - t.pw * M.c0;
    ctype const c5 = -t.vx * M.c7 - t.vy * M.c15 + t.vz * M.c5 - t.mx * M.c13 -
                     t.my * M.c3 + t.mz * M.c11 - t.px * M.c10 - t.py * M.c4 +
                     t.pz * M.c8 - t.pw * M.c1;
    ctype const c6 = t.vx * M.c6 - t.vy * M.c5 - t.vz * M.c15 + t.mx * M.c12 -
                     t.my * M.c11 - t.mz * M.c3 + t.px * M.c9 - t.py * M.c8 -
                     t.pz * M.c4 - t.pw * M.c2;
    ctype const c7 = -t.vx * M.c4 + t.vy * M.c10 - t.vz * M.c9 - t.mx * M.c14 -
                     t.mx * M.c15 + t.my * M.c2 + t.my * M.c7 - t.mz * M.c1 -
                     t.mz * M.c6 + t.pw * M.c8;
    ctype const c8 = -t.vx * M.c10 - t.vy * M.c4 + t.vz * M.c8 - t.mx * M.c2 -
                     t.mx * M.c7 - t.my * M.c14 - t.my * M.c15 + t.mz * M.c0 +
                     t.mz * M.c5 + t.pw * M.c9;
    ctype const c9 = t.vx * M.c9 - t.vy * M.c8 - t.vz * M.c4 + t.mx * M.c1 + t.mx * M.c6 -
                     t.my * M.c0 - t.my * M.c5 - t.mz * M.c14 - t.mz * M.c15 +
                     t.pw * M.c10;
    ctype const c10 = -t.vx * M.c0 - t.vy * M.c1 - t.vz * M.c2 + t.mx * M.c11 +
                      t.my * M.c12 + t.mz * M.c13 - t.px * M.c8 - t.py * M.c9 -
                      t.pz * M.c10 + t.pw * M.c15;
    ctype const c11 = t.vx * M.c14 - t.vy * M.c2 + t.vz * M.c1 - t.mx * M.c3 +
                      t.my * M.c13 - t.mz * M.c12 + t.px * M.c4 - t.py * M.c10 +
                      t.pz * M.c9 + t.pw * M.c5;
    ctype const c12 = t.vx * M.c2 + t.vy * M.c14 - t.vz * M.c0 - t.mx * M.c13 -
                      t.my * M.c3 + t.mz * M.c11 + t.px * M.c10 + t.py * M.c4 -
                      t.pz * M.c8 + t.pw * M.c6;
    ctype const c13 = -t.vx * M.c1 + t.vy * M.c0 + t.vz * M.c14 + t.mx * M.c12 -
                      t.my * M.c11 - t.mz * M.c3 - t.px * M.c9 + t.py * M.c8 +
                      t.pz * M.c4 + t.pw * M.c7;
    ctype const c14 = t.vx * M.c8 + t.vy * M.c9 + t.vz * M.c10 + t.mx * M.c0 +
                      t.mx * M.c5 + t.my * M.c1 + t.my * M.c6 + t.mz * M.c2 +
                      t.mz * M.c7 + t.pw * M.c4;
    ctype const c15 = t.vx * M.c11 + t.vy * M.c12 + t.vz * M.c13 - t.px * M.c0 +
                      t.px * M.c5 - t.py * M.c1 + t.py * M.c6 - t.pz * M.c2 +
                      t.pz * M.c7 - t.pw * M.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& M,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.mz + M.c2 * B.my + M.c3 * B.px - M.c4 * B.vx -
                     M.c5 * B.pw + M.c9 * B.vz - M.c10 * B.vy - M.c12 * B.pz +
                     M.c13 * B.py + M.c14 * B.mx;
    ctype const c1 = M.c0 * B.mz - M.c2 * B.mx + M.c3 * B.py - M.c4 * B.vy - M.c6 * B.pw -
                     M.c8 * B.vz + M.c10 * B.vx + M.c11 * B.pz - M.c13 * B.px +
                     M.c14 * B.my;
    ctype const c2 = -M.c0 * B.my + M.c1 * B.mx + M.c3 * B.pz - M.c4 * B.vz -
                     M.c7 * B.pw + M.c8 * B.vy - M.c9 * B.vx - M.c11 * B.py +
                     M.c12 * B.px + M.c14 * B.mz;
    ctype const c3 = -M.c0 * B.vx - M.c1 * B.vy - M.c2 * B.vz + M.c3 * B.pw -
                     M.c5 * B.vx - M.c6 * B.vy - M.c7 * B.vz - M.c11 * B.mx -
                     M.c12 * B.my - M.c13 * B.mz;
    ctype const c4 = M.c0 * B.px + M.c1 * B.py + M.c2 * B.pz - M.c4 * B.pw - M.c5 * B.px -
                     M.c6 * B.py - M.c7 * B.pz - M.c8 * B.mx - M.c9 * B.my - M.c10 * B.mz;
    ctype const c5 = -M.c0 * B.pw + M.c3 * B.px + M.c4 * B.vx - M.c6 * B.mz +
                     M.c7 * B.my - M.c9 * B.vz + M.c10 * B.vy - M.c12 * B.pz +
                     M.c13 * B.py + M.c15 * B.mx;
    ctype const c6 = -M.c1 * B.pw + M.c3 * B.py + M.c4 * B.vy + M.c5 * B.mz -
                     M.c7 * B.mx + M.c8 * B.vz - M.c10 * B.vx + M.c11 * B.pz -
                     M.c13 * B.px + M.c15 * B.my;
    ctype const c7 = -M.c2 * B.pw + M.c3 * B.pz + M.c4 * B.vz - M.c5 * B.my +
                     M.c6 * B.mx - M.c8 * B.vy + M.c9 * B.vx - M.c11 * B.py +
                     M.c12 * B.px + M.c15 * B.mz;
    ctype const c8 = M.c1 * B.pz - M.c2 * B.py + M.c4 * B.mx - M.c6 * B.pz + M.c7 * B.py -
                     M.c8 * B.pw - M.c9 * B.mz + M.c10 * B.my - M.c14 * B.px +
                     M.c15 * B.px;
    ctype const c9 = -M.c0 * B.pz + M.c2 * B.px + M.c4 * B.my + M.c5 * B.pz -
                     M.c7 * B.px + M.c8 * B.mz - M.c9 * B.pw - M.c10 * B.mx -
                     M.c14 * B.py + M.c15 * B.py;
    ctype const c10 = M.c0 * B.py - M.c1 * B.px + M.c4 * B.mz - M.c5 * B.py +
                      M.c6 * B.px - M.c8 * B.my + M.c9 * B.mx - M.c10 * B.pw -
                      M.c14 * B.pz + M.c15 * B.pz;
    ctype const c11 = -M.c1 * B.vz + M.c2 * B.vy + M.c3 * B.mx - M.c6 * B.vz +
                      M.c7 * B.vy + M.c11 * B.pw - M.c12 * B.mz + M.c13 * B.my +
                      M.c14 * B.vx + M.c15 * B.vx;
    ctype const c12 = M.c0 * B.vz - M.c2 * B.vx + M.c3 * B.my + M.c5 * B.vz -
                      M.c7 * B.vx + M.c11 * B.mz + M.c12 * B.pw - M.c13 * B.mx +
                      M.c14 * B.vy + M.c15 * B.vy;
    ctype const c13 = -M.c0 * B.vy + M.c1 * B.vx + M.c3 * B.mz - M.c5 * B.vy +
                      M.c6 * B.vx - M.c11 * B.my + M.c12 * B.mx + M.c13 * B.pw +
                      M.c14 * B.vz + M.c15 * B.vz;
    ctype const c14 = -M.c0 * B.mx - M.c1 * B.my - M.c2 * B.mz + M.c8 * B.vx +
                      M.c9 * B.vy + M.c10 * B.vz - M.c11 * B.px - M.c12 * B.py -
                      M.c13 * B.pz - M.c15 * B.pw;
    ctype const c15 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                      M.c9 * B.vy - M.c10 * B.vz - M.c11 * B.px - M.c12 * B.py -
                      M.c13 * B.pz - M.c14 * B.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c4 - B.vy * M.c10 + B.vz * M.c9 + B.mx * M.c14 -
                     B.my * M.c2 + B.mz * M.c1 - B.px * M.c3 + B.py * M.c13 -
                     B.pz * M.c12 - B.pw * M.c5;
    ctype const c1 = B.vx * M.c10 + B.vy * M.c4 - B.vz * M.c8 + B.mx * M.c2 +
                     B.my * M.c14 - B.mz * M.c0 - B.px * M.c13 - B.py * M.c3 +
                     B.pz * M.c11 - B.pw * M.c6;
    ctype const c2 = -B.vx * M.c9 + B.vy * M.c8 + B.vz * M.c4 - B.mx * M.c1 +
                     B.my * M.c0 + B.mz * M.c14 + B.px * M.c12 - B.py * M.c11 -
                     B.pz * M.c3 - B.pw * M.c7;
    ctype const c3 = B.vx * M.c0 - B.vx * M.c5 + B.vy * M.c1 - B.vy * M.c6 + B.vz * M.c2 -
                     B.vz * M.c7 - B.mx * M.c11 - B.my * M.c12 - B.mz * M.c13 -
                     B.pw * M.c3;
    ctype const c4 = -B.mx * M.c8 - B.my * M.c9 - B.mz * M.c10 - B.px * M.c0 -
                     B.px * M.c5 - B.py * M.c1 - B.py * M.c6 - B.pz * M.c2 - B.pz * M.c7 +
                     B.pw * M.c4;
    ctype const c5 = B.vx * M.c4 - B.vy * M.c10 + B.vz * M.c9 + B.mx * M.c15 -
                     B.my * M.c7 + B.mz * M.c6 + B.px * M.c3 - B.py * M.c13 +
                     B.pz * M.c12 - B.pw * M.c0;
    ctype const c6 = B.vx * M.c10 + B.vy * M.c4 - B.vz * M.c8 + B.mx * M.c7 +
                     B.my * M.c15 - B.mz * M.c5 + B.px * M.c13 + B.py * M.c3 -
                     B.pz * M.c11 - B.pw * M.c1;
    ctype const c7 = -B.vx * M.c9 + B.vy * M.c8 + B.vz * M.c4 - B.mx * M.c6 +
                     B.my * M.c5 + B.mz * M.c15 - B.px * M.c12 + B.py * M.c11 +
                     B.pz * M.c3 - B.pw * M.c2;
    ctype const c8 = B.mx * M.c4 - B.my * M.c10 + B.mz * M.c9 + B.px * M.c14 +
                     B.px * M.c15 - B.py * M.c2 - B.py * M.c7 + B.pz * M.c1 +
                     B.pz * M.c6 + B.pw * M.c8;
    ctype const c9 = B.mx * M.c10 + B.my * M.c4 - B.mz * M.c8 + B.px * M.c2 +
                     B.px * M.c7 + B.py * M.c14 + B.py * M.c15 - B.pz * M.c0 -
                     B.pz * M.c5 + B.pw * M.c9;
    ctype const c10 = -B.mx * M.c9 + B.my * M.c8 + B.mz * M.c4 - B.px * M.c1 -
                      B.px * M.c6 + B.py * M.c0 + B.py * M.c5 + B.pz * M.c14 +
                      B.pz * M.c15 + B.pw * M.c10;
    ctype const c11 = -B.vx * M.c14 + B.vx * M.c15 + B.vy * M.c2 - B.vy * M.c7 -
                      B.vz * M.c1 + B.vz * M.c6 + B.mx * M.c3 - B.my * M.c13 +
                      B.mz * M.c12 - B.pw * M.c11;
    ctype const c12 = -B.vx * M.c2 + B.vx * M.c7 - B.vy * M.c14 + B.vy * M.c15 +
                      B.vz * M.c0 - B.vz * M.c5 + B.mx * M.c13 + B.my * M.c3 -
                      B.mz * M.c11 - B.pw * M.c12;
    ctype const c13 = B.vx * M.c1 - B.vx * M.c6 - B.vy * M.c0 + B.vy * M.c5 -
                      B.vz * M.c14 + B.vz * M.c15 - B.mx * M.c12 + B.my * M.c11 +
                      B.mz * M.c3 - B.pw * M.c13;
    ctype const c14 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c0 -
                      B.my * M.c1 - B.mz * M.c2 + B.px * M.c11 + B.py * M.c12 +
                      B.pz * M.c13 - B.pw * M.c15;
    ctype const c15 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                      B.my * M.c6 - B.mz * M.c7 - B.px * M.c11 - B.py * M.c12 -
                      B.pz * M.c13 - B.pw * M.c14;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& M,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y + M.c2 * v.z - M.c3 * v.u - M.c4 * v.w;
    ctype const c1 = -M.c0 * v.w + M.c3 * v.x - M.c5 * v.w - M.c12 * v.z + M.c13 * v.y;
    ctype const c2 = -M.c1 * v.w + M.c3 * v.y - M.c6 * v.w + M.c11 * v.z - M.c13 * v.x;
    ctype const c3 = -M.c2 * v.w + M.c3 * v.z - M.c7 * v.w - M.c11 * v.y + M.c12 * v.x;
    ctype const c4 = M.c1 * v.z - M.c2 * v.y - M.c8 * v.w - M.c11 * v.u - M.c14 * v.x;
    ctype const c5 = -M.c0 * v.z + M.c2 * v.x - M.c9 * v.w - M.c12 * v.u - M.c14 * v.y;
    ctype const c6 = M.c0 * v.y - M.c1 * v.x - M.c10 * v.w - M.c13 * v.u - M.c14 * v.z;
    ctype const c7 = M.c0 * v.u - M.c4 * v.x - M.c5 * v.u + M.c9 * v.z - M.c10 * v.y;
    ctype const c8 = M.c1 * v.u - M.c4 * v.y - M.c6 * v.u - M.c8 * v.z + M.c10 * v.x;
    ctype const c9 = M.c2 * v.u - M.c4 * v.z - M.c7 * v.u + M.c8 * v.y - M.c9 * v.x;
    ctype const c10 = M.c3 * v.u - M.c4 * v.w - M.c5 * v.x - M.c6 * v.y - M.c7 * v.z;
    ctype const c11 = -M.c6 * v.z + M.c7 * v.y - M.c8 * v.w + M.c11 * v.u + M.c15 * v.x;
    ctype const c12 = M.c5 * v.z - M.c7 * v.x - M.c9 * v.w + M.c12 * v.u + M.c15 * v.y;
    ctype const c13 = -M.c5 * v.y + M.c6 * v.x - M.c10 * v.w + M.c13 * v.u + M.c15 * v.z;
    ctype const c14 = M.c8 * v.x + M.c9 * v.y + M.c10 * v.z + M.c14 * v.u - M.c15 * v.u;
    ctype const c15 =
        -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z - M.c14 * v.w - M.c15 * v.w;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 + v.z * M.c2 - v.w * M.c4 - v.u * M.c3;
    ctype const c1 = -v.x * M.c3 + v.y * M.c13 - v.z * M.c12 + v.w * M.c0 - v.w * M.c5;
    ctype const c2 = -v.x * M.c13 - v.y * M.c3 + v.z * M.c11 + v.w * M.c1 - v.w * M.c6;
    ctype const c3 = v.x * M.c12 - v.y * M.c11 - v.z * M.c3 + v.w * M.c2 - v.w * M.c7;
    ctype const c4 = -v.x * M.c14 + v.y * M.c2 - v.z * M.c1 - v.w * M.c8 - v.u * M.c11;
    ctype const c5 = -v.x * M.c2 - v.y * M.c14 + v.z * M.c0 - v.w * M.c9 - v.u * M.c12;
    ctype const c6 = v.x * M.c1 - v.y * M.c0 - v.z * M.c14 - v.w * M.c10 - v.u * M.c13;
    ctype const c7 = v.x * M.c4 - v.y * M.c10 + v.z * M.c9 - v.u * M.c0 - v.u * M.c5;
    ctype const c8 = v.x * M.c10 + v.y * M.c4 - v.z * M.c8 - v.u * M.c1 - v.u * M.c6;
    ctype const c9 = -v.x * M.c9 + v.y * M.c8 + v.z * M.c4 - v.u * M.c2 - v.u * M.c7;
    ctype const c10 = -v.x * M.c5 - v.y * M.c6 - v.z * M.c7 + v.w * M.c4 - v.u * M.c3;
    ctype const c11 = v.x * M.c15 - v.y * M.c7 + v.z * M.c6 + v.w * M.c8 - v.u * M.c11;
    ctype const c12 = v.x * M.c7 + v.y * M.c15 - v.z * M.c5 + v.w * M.c9 - v.u * M.c12;
    ctype const c13 = -v.x * M.c6 + v.y * M.c5 + v.z * M.c15 + v.w * M.c10 - v.u * M.c13;
    ctype const c14 = -v.x * M.c8 - v.y * M.c9 - v.z * M.c10 - v.u * M.c14 - v.u * M.c15;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14 - v.w * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(MVec3dc_U<T> const& M,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    ctype const c4 = M.c4 * ctype(s);
    ctype const c5 = M.c5 * ctype(s);
    ctype const c6 = M.c6 * ctype(s);
    ctype const c7 = M.c7 * ctype(s);
    ctype const c8 = M.c8 * ctype(s);
    ctype const c9 = M.c9 * ctype(s);
    ctype const c10 = M.c10 * ctype(s);
    ctype const c11 = M.c11 * ctype(s);
    ctype const c12 = M.c12 * ctype(s);
    ctype const c13 = M.c13 * ctype(s);
    ctype const c14 = M.c14 * ctype(s);
    ctype const c15 = M.c15 * ctype(s);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(Scalar3dc<T> s,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    ctype const c4 = ctype(s) * M.c4;
    ctype const c5 = ctype(s) * M.c5;
    ctype const c6 = ctype(s) * M.c6;
    ctype const c7 = ctype(s) * M.c7;
    ctype const c8 = ctype(s) * M.c8;
    ctype const c9 = ctype(s) * M.c9;
    ctype const c10 = ctype(s) * M.c10;
    ctype const c11 = ctype(s) * M.c11;
    ctype const c12 = ctype(s) * M.c12;
    ctype const c13 = ctype(s) * M.c13;
    ctype const c14 = ctype(s) * M.c14;
    ctype const c15 = ctype(s) * M.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps1,
                                                        PScalar3dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga3dc gpr :: gpr(ps,quadvec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps,
                                                     QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * Q.x;
    ctype const c1 = -ctype(ps) * Q.y;
    ctype const c2 = -ctype(ps) * Q.z;
    ctype const c3 = ctype(ps) * Q.u;
    ctype const c4 = ctype(ps) * Q.w;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(quadvec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                     PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q.x * ctype(ps);
    ctype const c1 = -Q.y * ctype(ps);
    ctype const c2 = -Q.z * ctype(ps);
    ctype const c3 = Q.u * ctype(ps);
    ctype const c4 = Q.w * ctype(ps);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(ps,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps,
                                                       TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * t.px;
    ctype const c1 = -ctype(ps) * t.py;
    ctype const c2 = -ctype(ps) * t.pz;
    ctype const c3 = -ctype(ps) * t.vx;
    ctype const c4 = -ctype(ps) * t.vy;
    ctype const c5 = -ctype(ps) * t.vz;
    ctype const c6 = -ctype(ps) * t.mx;
    ctype const c7 = -ctype(ps) * t.my;
    ctype const c8 = -ctype(ps) * t.mz;
    ctype const c9 = ctype(ps) * t.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(trivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                       PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.px * ctype(ps);
    ctype const c1 = -t.py * ctype(ps);
    ctype const c2 = -t.pz * ctype(ps);
    ctype const c3 = -t.vx * ctype(ps);
    ctype const c4 = -t.vy * ctype(ps);
    ctype const c5 = -t.vz * ctype(ps);
    ctype const c6 = -t.mx * ctype(ps);
    ctype const c7 = -t.my * ctype(ps);
    ctype const c8 = -t.mz * ctype(ps);
    ctype const c9 = t.pw * ctype(ps);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(ps,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.mx;
    ctype const c1 = ctype(ps) * B.my;
    ctype const c2 = ctype(ps) * B.mz;
    ctype const c3 = ctype(ps) * B.px;
    ctype const c4 = ctype(ps) * B.py;
    ctype const c5 = ctype(ps) * B.pz;
    ctype const c6 = ctype(ps) * B.vx;
    ctype const c7 = ctype(ps) * B.vy;
    ctype const c8 = ctype(ps) * B.vz;
    ctype const c9 = -ctype(ps) * B.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(bivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.mx * ctype(ps);
    ctype const c1 = B.my * ctype(ps);
    ctype const c2 = B.mz * ctype(ps);
    ctype const c3 = B.px * ctype(ps);
    ctype const c4 = B.py * ctype(ps);
    ctype const c5 = B.pz * ctype(ps);
    ctype const c6 = B.vx * ctype(ps);
    ctype const c7 = B.vy * ctype(ps);
    ctype const c8 = B.vz * ctype(ps);
    ctype const c9 = -B.pw * ctype(ps);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(ps,vec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps,
                                                         Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = -ctype(ps) * v.u;
    ctype const c4 = -ctype(ps) * v.w;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(vec,ps) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v,
                                                         PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = -v.u * ctype(ps);
    ctype const c4 = -v.w * ctype(ps);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> operator*(PScalar3dc<T> ps, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(ps) * ctype(s));
}

// cga3dc gpr :: gpr(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> operator*(Scalar3dc<T> s, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(s) * ctype(ps));
}

// cga3dc gpr :: gpr(quadvec,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q1,
                                                        QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -Q1.x * Q2.x - Q1.y * Q2.y - Q1.z * Q2.z + Q1.w * Q2.u + Q1.u * Q2.w;
    ctype const c1 = -Q1.x * Q2.u + Q1.u * Q2.x;
    ctype const c2 = -Q1.y * Q2.u + Q1.u * Q2.y;
    ctype const c3 = -Q1.z * Q2.u + Q1.u * Q2.z;
    ctype const c4 = -Q1.y * Q2.z + Q1.z * Q2.y;
    ctype const c5 = Q1.x * Q2.z - Q1.z * Q2.x;
    ctype const c6 = -Q1.x * Q2.y + Q1.y * Q2.x;
    ctype const c7 = Q1.x * Q2.w - Q1.w * Q2.x;
    ctype const c8 = Q1.y * Q2.w - Q1.w * Q2.y;
    ctype const c9 = Q1.z * Q2.w - Q1.w * Q2.z;
    ctype const c10 = Q1.w * Q2.u - Q1.u * Q2.w;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(quadvec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * t.vz - Q.z * t.vy - Q.w * t.px + Q.u * t.mx;
    ctype const c1 = -Q.x * t.vz + Q.z * t.vx - Q.w * t.py + Q.u * t.my;
    ctype const c2 = Q.x * t.vy - Q.y * t.vx - Q.w * t.pz + Q.u * t.mz;
    ctype const c3 = Q.x * t.px + Q.y * t.py + Q.z * t.pz - Q.u * t.pw;
    ctype const c4 = -Q.x * t.mx - Q.y * t.my - Q.z * t.mz + Q.w * t.pw;
    ctype const c5 = -Q.x * t.pw + Q.w * t.px + Q.u * t.mx;
    ctype const c6 = -Q.y * t.pw + Q.w * t.py + Q.u * t.my;
    ctype const c7 = -Q.z * t.pw + Q.w * t.pz + Q.u * t.mz;
    ctype const c8 = -Q.y * t.mz + Q.z * t.my + Q.w * t.vx;
    ctype const c9 = Q.x * t.mz - Q.z * t.mx + Q.w * t.vy;
    ctype const c10 = -Q.x * t.my + Q.y * t.mx + Q.w * t.vz;
    ctype const c11 = Q.y * t.pz - Q.z * t.py + Q.u * t.vx;
    ctype const c12 = -Q.x * t.pz + Q.z * t.px + Q.u * t.vy;
    ctype const c13 = Q.x * t.py - Q.y * t.px + Q.u * t.vz;
    ctype const c14 = Q.x * t.vx + Q.y * t.vy + Q.z * t.vz;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(trivec,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vy * Q.z - t.vz * Q.y - t.mx * Q.u + t.px * Q.w;
    ctype const c1 = -t.vx * Q.z + t.vz * Q.x - t.my * Q.u + t.py * Q.w;
    ctype const c2 = t.vx * Q.y - t.vy * Q.x - t.mz * Q.u + t.pz * Q.w;
    ctype const c3 = -t.px * Q.x - t.py * Q.y - t.pz * Q.z + t.pw * Q.u;
    ctype const c4 = t.mx * Q.x + t.my * Q.y + t.mz * Q.z - t.pw * Q.w;
    ctype const c5 = t.mx * Q.u + t.px * Q.w - t.pw * Q.x;
    ctype const c6 = t.my * Q.u + t.py * Q.w - t.pw * Q.y;
    ctype const c7 = t.mz * Q.u + t.pz * Q.w - t.pw * Q.z;
    ctype const c8 = t.vx * Q.w + t.my * Q.z - t.mz * Q.y;
    ctype const c9 = t.vy * Q.w - t.mx * Q.z + t.mz * Q.x;
    ctype const c10 = t.vz * Q.w + t.mx * Q.y - t.my * Q.x;
    ctype const c11 = t.vx * Q.u - t.py * Q.z + t.pz * Q.y;
    ctype const c12 = t.vy * Q.u + t.px * Q.z - t.pz * Q.x;
    ctype const c13 = t.vz * Q.u - t.px * Q.y + t.py * Q.x;
    ctype const c14 = t.vx * Q.x + t.vy * Q.y + t.vz * Q.z;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(quadvec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = Q.y * B.vz - Q.z * B.vy + Q.u * B.mx;
    ctype const c2 = -Q.x * B.vz + Q.z * B.vx + Q.u * B.my;
    ctype const c3 = Q.x * B.vy - Q.y * B.vx + Q.u * B.mz;
    ctype const c4 = Q.x * B.pw + Q.w * B.vx + Q.u * B.px;
    ctype const c5 = Q.y * B.pw + Q.w * B.vy + Q.u * B.py;
    ctype const c6 = Q.z * B.pw + Q.w * B.vz + Q.u * B.pz;
    ctype const c7 = -Q.y * B.pz + Q.z * B.py + Q.w * B.mx;
    ctype const c8 = Q.x * B.pz - Q.z * B.px + Q.w * B.my;
    ctype const c9 = -Q.x * B.py + Q.y * B.px + Q.w * B.mz;
    ctype const c10 = -Q.x * B.mx - Q.y * B.my - Q.z * B.mz;
    ctype const c11 = -Q.y * B.mz + Q.z * B.my + Q.w * B.vx - Q.u * B.px;
    ctype const c12 = Q.x * B.mz - Q.z * B.mx + Q.w * B.vy - Q.u * B.py;
    ctype const c13 = -Q.x * B.my + Q.y * B.mx + Q.w * B.vz - Q.u * B.pz;
    ctype const c14 = -Q.x * B.px - Q.y * B.py - Q.z * B.pz - Q.w * B.pw;
    ctype const c15 = Q.x * B.vx + Q.y * B.vy + Q.z * B.vz + Q.u * B.pw;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(bivec,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -B.vy * Q.z + B.vz * Q.y + B.mx * Q.u;
    ctype const c2 = B.vx * Q.z - B.vz * Q.x + B.my * Q.u;
    ctype const c3 = -B.vx * Q.y + B.vy * Q.x + B.mz * Q.u;
    ctype const c4 = B.vx * Q.w + B.px * Q.u + B.pw * Q.x;
    ctype const c5 = B.vy * Q.w + B.py * Q.u + B.pw * Q.y;
    ctype const c6 = B.vz * Q.w + B.pz * Q.u + B.pw * Q.z;
    ctype const c7 = B.mx * Q.w + B.py * Q.z - B.pz * Q.y;
    ctype const c8 = B.my * Q.w - B.px * Q.z + B.pz * Q.x;
    ctype const c9 = B.mz * Q.w + B.px * Q.y - B.py * Q.x;
    ctype const c10 = -B.mx * Q.x - B.my * Q.y - B.mz * Q.z;
    ctype const c11 = -B.vx * Q.w - B.my * Q.z + B.mz * Q.y + B.px * Q.u;
    ctype const c12 = -B.vy * Q.w + B.mx * Q.z - B.mz * Q.x + B.py * Q.u;
    ctype const c13 = -B.vz * Q.w - B.mx * Q.y + B.my * Q.x + B.pz * Q.u;
    ctype const c14 = B.px * Q.x + B.py * Q.y + B.pz * Q.z + B.pw * Q.w;
    ctype const c15 = -B.vx * Q.x - B.vy * Q.y - B.vz * Q.z - B.pw * Q.u;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(quadvec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = Q.y * v.z - Q.z * v.y;
    ctype const c6 = -Q.x * v.z + Q.z * v.x;
    ctype const c7 = Q.x * v.y - Q.y * v.x;
    ctype const c8 = Q.x * v.u + Q.w * v.x;
    ctype const c9 = Q.y * v.u + Q.w * v.y;
    ctype const c10 = Q.z * v.u + Q.w * v.z;
    ctype const c11 = -Q.x * v.w - Q.u * v.x;
    ctype const c12 = -Q.y * v.w - Q.u * v.y;
    ctype const c13 = -Q.z * v.w - Q.u * v.z;
    ctype const c14 = -Q.w * v.w + Q.u * v.u;
    ctype const c15 = Q.x * v.x + Q.y * v.y + Q.z * v.z + Q.w * v.w + Q.u * v.u;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(vec,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = v.y * Q.z - v.z * Q.y;
    ctype const c6 = -v.x * Q.z + v.z * Q.x;
    ctype const c7 = v.x * Q.y - v.y * Q.x;
    ctype const c8 = -v.x * Q.w - v.u * Q.x;
    ctype const c9 = -v.y * Q.w - v.u * Q.y;
    ctype const c10 = -v.z * Q.w - v.u * Q.z;
    ctype const c11 = v.x * Q.u + v.w * Q.x;
    ctype const c12 = v.y * Q.u + v.w * Q.y;
    ctype const c13 = v.z * Q.u + v.w * Q.z;
    ctype const c14 = v.w * Q.w - v.u * Q.u;
    ctype const c15 = v.x * Q.x + v.y * Q.y + v.z * Q.z + v.w * Q.w + v.u * Q.u;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(quadvec,s) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator*(QuadVec3dc<T> const& Q,
                                                         Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * ctype(s);
    ctype const c1 = Q.y * ctype(s);
    ctype const c2 = Q.z * ctype(s);
    ctype const c3 = Q.w * ctype(s);
    ctype const c4 = Q.u * ctype(s);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(s,quadvec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> operator*(Scalar3dc<T> s,
                                                         QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * Q.x;
    ctype const c1 = ctype(s) * Q.y;
    ctype const c2 = ctype(s) * Q.z;
    ctype const c3 = ctype(s) * Q.w;
    ctype const c4 = ctype(s) * Q.u;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(trivec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t1,
                                                        TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.vx * t2.vx + t1.vy * t2.vy + t1.vz * t2.vz + t1.mx * t2.px +
                     t1.my * t2.py + t1.mz * t2.pz + t1.px * t2.mx + t1.py * t2.my +
                     t1.pz * t2.mz - t1.pw * t2.pw;
    ctype const c1 = t1.vy * t2.pz - t1.vz * t2.py + t1.px * t2.pw + t1.py * t2.vz -
                     t1.pz * t2.vy - t1.pw * t2.px;
    ctype const c2 = -t1.vx * t2.pz + t1.vz * t2.px - t1.px * t2.vz + t1.py * t2.pw +
                     t1.pz * t2.vx - t1.pw * t2.py;
    ctype const c3 = t1.vx * t2.py - t1.vy * t2.px + t1.px * t2.vy - t1.py * t2.vx +
                     t1.pz * t2.pw - t1.pw * t2.pz;
    ctype const c4 = t1.vy * t2.vz - t1.vz * t2.vy + t1.my * t2.pz - t1.mz * t2.py +
                     t1.py * t2.mz - t1.pz * t2.my;
    ctype const c5 = -t1.vx * t2.vz + t1.vz * t2.vx - t1.mx * t2.pz + t1.mz * t2.px -
                     t1.px * t2.mz + t1.pz * t2.mx;
    ctype const c6 = t1.vx * t2.vy - t1.vy * t2.vx + t1.mx * t2.py - t1.my * t2.px +
                     t1.px * t2.my - t1.py * t2.mx;
    ctype const c7 = t1.vy * t2.mz - t1.vz * t2.my - t1.mx * t2.pw + t1.my * t2.vz -
                     t1.mz * t2.vy + t1.pw * t2.mx;
    ctype const c8 = -t1.vx * t2.mz + t1.vz * t2.mx - t1.mx * t2.vz - t1.my * t2.pw +
                     t1.mz * t2.vx + t1.pw * t2.my;
    ctype const c9 = t1.vx * t2.my - t1.vy * t2.mx + t1.mx * t2.vy - t1.my * t2.vx -
                     t1.mz * t2.pw + t1.pw * t2.mz;
    ctype const c10 = t1.mx * t2.px + t1.my * t2.py + t1.mz * t2.pz - t1.px * t2.mx -
                      t1.py * t2.my - t1.pz * t2.mz;
    ctype const c11 = t1.vx * t2.pw + t1.my * t2.pz - t1.mz * t2.py - t1.py * t2.mz +
                      t1.pz * t2.my + t1.pw * t2.vx;
    ctype const c12 = t1.vy * t2.pw - t1.mx * t2.pz + t1.mz * t2.px + t1.px * t2.mz -
                      t1.pz * t2.mx + t1.pw * t2.vy;
    ctype const c13 = t1.vz * t2.pw + t1.mx * t2.py - t1.my * t2.px - t1.px * t2.my +
                      t1.py * t2.mx + t1.pw * t2.vz;
    ctype const c14 = t1.vx * t2.mx + t1.vy * t2.my + t1.vz * t2.mz + t1.mx * t2.vx +
                      t1.my * t2.vy + t1.mz * t2.vz;
    ctype const c15 = t1.vx * t2.px + t1.vy * t2.py + t1.vz * t2.pz + t1.px * t2.vx +
                      t1.py * t2.vy + t1.pz * t2.vz;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(trivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * B.pw + t.my * B.vz - t.mz * B.vy - t.py * B.pz +
                     t.pz * B.py + t.pw * B.mx;
    ctype const c1 = -t.vy * B.pw - t.mx * B.vz + t.mz * B.vx + t.px * B.pz -
                     t.pz * B.px + t.pw * B.my;
    ctype const c2 = -t.vz * B.pw + t.mx * B.vy - t.my * B.vx - t.px * B.py +
                     t.py * B.px + t.pw * B.mz;
    ctype const c3 = -t.vx * B.vx - t.vy * B.vy - t.vz * B.vz - t.px * B.mx -
                     t.py * B.my - t.pz * B.mz;
    ctype const c4 = -t.vx * B.px - t.vy * B.py - t.vz * B.pz - t.mx * B.mx -
                     t.my * B.my - t.mz * B.mz;
    ctype const c5 = -t.vy * B.mz + t.vz * B.my - t.my * B.vz + t.mz * B.vy -
                     t.py * B.pz + t.pz * B.py;
    ctype const c6 =
        t.vx * B.mz - t.vz * B.mx + t.mx * B.vz - t.mz * B.vx + t.px * B.pz - t.pz * B.px;
    ctype const c7 = -t.vx * B.my + t.vy * B.mx - t.mx * B.vy + t.my * B.vx -
                     t.px * B.py + t.py * B.px;
    ctype const c8 = -t.vy * B.pz + t.vz * B.py - t.mx * B.pw - t.my * B.mz +
                     t.mz * B.my - t.pw * B.px;
    ctype const c9 =
        t.vx * B.pz - t.vz * B.px + t.mx * B.mz - t.my * B.pw - t.mz * B.mx - t.pw * B.py;
    ctype const c10 = -t.vx * B.py + t.vy * B.px - t.mx * B.my + t.my * B.mx -
                      t.mz * B.pw - t.pw * B.pz;
    ctype const c11 = -t.vy * B.vz + t.vz * B.vy + t.px * B.pw - t.py * B.mz +
                      t.pz * B.my + t.pw * B.vx;
    ctype const c12 =
        t.vx * B.vz - t.vz * B.vx + t.px * B.mz + t.py * B.pw - t.pz * B.mx + t.pw * B.vy;
    ctype const c13 = -t.vx * B.vy + t.vy * B.vx - t.px * B.my + t.py * B.mx +
                      t.pz * B.pw + t.pw * B.vz;
    ctype const c14 =
        t.mx * B.vx + t.my * B.vy + t.mz * B.vz - t.px * B.px - t.py * B.py - t.pz * B.pz;
    ctype const c15 = -t.vx * B.mx - t.vy * B.my - t.vz * B.mz - t.mx * B.vx -
                      t.my * B.vy - t.mz * B.vz - t.px * B.px - t.py * B.py -
                      t.pz * B.pz - t.pw * B.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(bivec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * t.mz + B.vz * t.my + B.mx * t.pw + B.py * t.pz -
                     B.pz * t.py - B.pw * t.vx;
    ctype const c1 =
        B.vx * t.mz - B.vz * t.mx + B.my * t.pw - B.px * t.pz + B.pz * t.px - B.pw * t.vy;
    ctype const c2 = -B.vx * t.my + B.vy * t.mx + B.mz * t.pw + B.px * t.py -
                     B.py * t.px - B.pw * t.vz;
    ctype const c3 = -B.vx * t.vx - B.vy * t.vy - B.vz * t.vz - B.mx * t.px -
                     B.my * t.py - B.mz * t.pz;
    ctype const c4 = -B.mx * t.mx - B.my * t.my - B.mz * t.mz - B.px * t.vx -
                     B.py * t.vy - B.pz * t.vz;
    ctype const c5 = -B.vy * t.mz + B.vz * t.my - B.my * t.vz + B.mz * t.vy -
                     B.py * t.pz + B.pz * t.py;
    ctype const c6 =
        B.vx * t.mz - B.vz * t.mx + B.mx * t.vz - B.mz * t.vx + B.px * t.pz - B.pz * t.px;
    ctype const c7 = -B.vx * t.my + B.vy * t.mx - B.mx * t.vy + B.my * t.vx -
                     B.px * t.py + B.py * t.px;
    ctype const c8 = -B.my * t.mz + B.mz * t.my + B.px * t.pw - B.py * t.vz +
                     B.pz * t.vy + B.pw * t.mx;
    ctype const c9 =
        B.mx * t.mz - B.mz * t.mx + B.px * t.vz + B.py * t.pw - B.pz * t.vx + B.pw * t.my;
    ctype const c10 = -B.mx * t.my + B.my * t.mx - B.px * t.vy + B.py * t.vx +
                      B.pz * t.pw + B.pw * t.mz;
    ctype const c11 = -B.vx * t.pw - B.vy * t.vz + B.vz * t.vy - B.my * t.pz +
                      B.mz * t.py - B.pw * t.px;
    ctype const c12 =
        B.vx * t.vz - B.vy * t.pw - B.vz * t.vx + B.mx * t.pz - B.mz * t.px - B.pw * t.py;
    ctype const c13 = -B.vx * t.vy + B.vy * t.vx - B.vz * t.pw - B.mx * t.py +
                      B.my * t.px - B.pw * t.pz;
    ctype const c14 = -B.vx * t.mx - B.vy * t.my - B.vz * t.mz + B.px * t.px +
                      B.py * t.py + B.pz * t.pz;
    ctype const c15 = -B.vx * t.mx - B.vy * t.my - B.vz * t.mz - B.mx * t.vx -
                      B.my * t.vy - B.mz * t.vz - B.px * t.px - B.py * t.py -
                      B.pz * t.pz - B.pw * t.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(trivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -t.vx * v.w - t.py * v.z + t.pz * v.y;
    ctype const c2 = -t.vy * v.w + t.px * v.z - t.pz * v.x;
    ctype const c3 = -t.vz * v.w - t.px * v.y + t.py * v.x;
    ctype const c4 = -t.mx * v.w - t.px * v.u - t.pw * v.x;
    ctype const c5 = -t.my * v.w - t.py * v.u - t.pw * v.y;
    ctype const c6 = -t.mz * v.w - t.pz * v.u - t.pw * v.z;
    ctype const c7 = -t.vx * v.u + t.my * v.z - t.mz * v.y;
    ctype const c8 = -t.vy * v.u - t.mx * v.z + t.mz * v.x;
    ctype const c9 = -t.vz * v.u + t.mx * v.y - t.my * v.x;
    ctype const c10 = -t.vx * v.x - t.vy * v.y - t.vz * v.z;
    ctype const c11 = -t.vy * v.z + t.vz * v.y - t.mx * v.w + t.px * v.u;
    ctype const c12 = t.vx * v.z - t.vz * v.x - t.my * v.w + t.py * v.u;
    ctype const c13 = -t.vx * v.y + t.vy * v.x - t.mz * v.w + t.pz * v.u;
    ctype const c14 = t.mx * v.x + t.my * v.y + t.mz * v.z + t.pw * v.u;
    ctype const c15 = -t.px * v.x - t.py * v.y - t.pz * v.z - t.pw * v.w;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(vec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.y * t.pz - v.z * t.py - v.w * t.vx;
    ctype const c2 = -v.x * t.pz + v.z * t.px - v.w * t.vy;
    ctype const c3 = v.x * t.py - v.y * t.px - v.w * t.vz;
    ctype const c4 = -v.x * t.pw - v.w * t.mx - v.u * t.px;
    ctype const c5 = -v.y * t.pw - v.w * t.my - v.u * t.py;
    ctype const c6 = -v.z * t.pw - v.w * t.mz - v.u * t.pz;
    ctype const c7 = -v.y * t.mz + v.z * t.my - v.u * t.vx;
    ctype const c8 = v.x * t.mz - v.z * t.mx - v.u * t.vy;
    ctype const c9 = -v.x * t.my + v.y * t.mx - v.u * t.vz;
    ctype const c10 = -v.x * t.vx - v.y * t.vy - v.z * t.vz;
    ctype const c11 = -v.y * t.vz + v.z * t.vy + v.w * t.mx - v.u * t.px;
    ctype const c12 = v.x * t.vz - v.z * t.vx + v.w * t.my - v.u * t.py;
    ctype const c13 = -v.x * t.vy + v.y * t.vx + v.w * t.mz - v.u * t.pz;
    ctype const c14 = -v.x * t.mx - v.y * t.my - v.z * t.mz - v.u * t.pw;
    ctype const c15 = v.x * t.px + v.y * t.py + v.z * t.pz + v.w * t.pw;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator*(TriVec3dc<T> const& t,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * ctype(s);
    ctype const c1 = t.vy * ctype(s);
    ctype const c2 = t.vz * ctype(s);
    ctype const c3 = t.mx * ctype(s);
    ctype const c4 = t.my * ctype(s);
    ctype const c5 = t.mz * ctype(s);
    ctype const c6 = t.px * ctype(s);
    ctype const c7 = t.py * ctype(s);
    ctype const c8 = t.pz * ctype(s);
    ctype const c9 = t.pw * ctype(s);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> operator*(Scalar3dc<T> s,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.vx;
    ctype const c1 = ctype(s) * t.vy;
    ctype const c2 = ctype(s) * t.vz;
    ctype const c3 = ctype(s) * t.mx;
    ctype const c4 = ctype(s) * t.my;
    ctype const c5 = ctype(s) * t.mz;
    ctype const c6 = ctype(s) * t.px;
    ctype const c7 = ctype(s) * t.py;
    ctype const c8 = ctype(s) * t.pz;
    ctype const c9 = ctype(s) * t.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(bivec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B1,
                                                        BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vx * B2.px - B1.vy * B2.py - B1.vz * B2.pz - B1.mx * B2.mx -
                     B1.my * B2.my - B1.mz * B2.mz - B1.px * B2.vx - B1.py * B2.vy -
                     B1.pz * B2.vz + B1.pw * B2.pw;
    ctype const c1 = B1.vx * B2.pw - B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz +
                     B1.mz * B2.vy - B1.pw * B2.vx;
    ctype const c2 = B1.vx * B2.mz + B1.vy * B2.pw - B1.vz * B2.mx + B1.mx * B2.vz -
                     B1.mz * B2.vx - B1.pw * B2.vy;
    ctype const c3 = -B1.vx * B2.my + B1.vy * B2.mx + B1.vz * B2.pw - B1.mx * B2.vy +
                     B1.my * B2.vx - B1.pw * B2.vz;
    ctype const c4 = -B1.vy * B2.pz + B1.vz * B2.py - B1.my * B2.mz + B1.mz * B2.my -
                     B1.py * B2.vz + B1.pz * B2.vy;
    ctype const c5 = B1.vx * B2.pz - B1.vz * B2.px + B1.mx * B2.mz - B1.mz * B2.mx +
                     B1.px * B2.vz - B1.pz * B2.vx;
    ctype const c6 = -B1.vx * B2.py + B1.vy * B2.px - B1.mx * B2.my + B1.my * B2.mx -
                     B1.px * B2.vy + B1.py * B2.vx;
    ctype const c7 = -B1.my * B2.pz + B1.mz * B2.py - B1.px * B2.pw - B1.py * B2.mz +
                     B1.pz * B2.my + B1.pw * B2.px;
    ctype const c8 = B1.mx * B2.pz - B1.mz * B2.px + B1.px * B2.mz - B1.py * B2.pw -
                     B1.pz * B2.mx + B1.pw * B2.py;
    ctype const c9 = -B1.mx * B2.py + B1.my * B2.px - B1.px * B2.my + B1.py * B2.mx -
                     B1.pz * B2.pw + B1.pw * B2.pz;
    ctype const c10 = B1.vx * B2.px + B1.vy * B2.py + B1.vz * B2.pz - B1.px * B2.vx -
                      B1.py * B2.vy - B1.pz * B2.vz;
    ctype const c11 = B1.vy * B2.pz - B1.vz * B2.py + B1.mx * B2.pw - B1.py * B2.vz +
                      B1.pz * B2.vy + B1.pw * B2.mx;
    ctype const c12 = -B1.vx * B2.pz + B1.vz * B2.px + B1.my * B2.pw + B1.px * B2.vz -
                      B1.pz * B2.vx + B1.pw * B2.my;
    ctype const c13 = B1.vx * B2.py - B1.vy * B2.px + B1.mz * B2.pw - B1.px * B2.vy +
                      B1.py * B2.vx + B1.pw * B2.mz;
    ctype const c14 = -B1.mx * B2.px - B1.my * B2.py - B1.mz * B2.pz - B1.px * B2.mx -
                      B1.py * B2.my - B1.pz * B2.mz;
    ctype const c15 = -B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz - B1.mx * B2.vx -
                      B1.my * B2.vy - B1.mz * B2.vz;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(bivec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.u - B.my * v.z + B.mz * v.y - B.px * v.w;
    ctype const c1 = B.vy * v.u + B.mx * v.z - B.mz * v.x - B.py * v.w;
    ctype const c2 = B.vz * v.u - B.mx * v.y + B.my * v.x - B.pz * v.w;
    ctype const c3 = B.vx * v.x + B.vy * v.y + B.vz * v.z - B.pw * v.w;
    ctype const c4 = -B.px * v.x - B.py * v.y - B.pz * v.z + B.pw * v.u;
    ctype const c5 = B.vx * v.u + B.px * v.w - B.pw * v.x;
    ctype const c6 = B.vy * v.u + B.py * v.w - B.pw * v.y;
    ctype const c7 = B.vz * v.u + B.pz * v.w - B.pw * v.z;
    ctype const c8 = B.mx * v.u - B.py * v.z + B.pz * v.y;
    ctype const c9 = B.my * v.u + B.px * v.z - B.pz * v.x;
    ctype const c10 = B.mz * v.u - B.px * v.y + B.py * v.x;
    ctype const c11 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c12 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c13 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c14 = -B.mx * v.x - B.my * v.y - B.mz * v.z;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(vec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my + v.w * B.px - v.u * B.vx;
    ctype const c1 = v.x * B.mz - v.z * B.mx + v.w * B.py - v.u * B.vy;
    ctype const c2 = -v.x * B.my + v.y * B.mx + v.w * B.pz - v.u * B.vz;
    ctype const c3 = -v.x * B.vx - v.y * B.vy - v.z * B.vz + v.w * B.pw;
    ctype const c4 = v.x * B.px + v.y * B.py + v.z * B.pz - v.u * B.pw;
    ctype const c5 = -v.x * B.pw + v.w * B.px + v.u * B.vx;
    ctype const c6 = -v.y * B.pw + v.w * B.py + v.u * B.vy;
    ctype const c7 = -v.z * B.pw + v.w * B.pz + v.u * B.vz;
    ctype const c8 = v.y * B.pz - v.z * B.py + v.u * B.mx;
    ctype const c9 = -v.x * B.pz + v.z * B.px + v.u * B.my;
    ctype const c10 = v.x * B.py - v.y * B.px + v.u * B.mz;
    ctype const c11 = -v.y * B.vz + v.z * B.vy + v.w * B.mx;
    ctype const c12 = v.x * B.vz - v.z * B.vx + v.w * B.my;
    ctype const c13 = -v.x * B.vy + v.y * B.vx + v.w * B.mz;
    ctype const c14 = -v.x * B.mx - v.y * B.my - v.z * B.mz;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc gpr :: gpr(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator*(BiVec3dc<T> const& B,
                                                       Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    ctype const c6 = B.px * ctype(s);
    ctype const c7 = B.py * ctype(s);
    ctype const c8 = B.pz * ctype(s);
    ctype const c9 = B.pw * ctype(s);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> operator*(Scalar3dc<T> s,
                                                       BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    ctype const c6 = ctype(s) * B.px;
    ctype const c7 = ctype(s) * B.py;
    ctype const c8 = ctype(s) * B.pz;
    ctype const c9 = ctype(s) * B.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc gpr :: gpr(vec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v1,
                                                        Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z - v1.w * v2.u - v1.u * v2.w;
    ctype const c1 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c2 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c3 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c4 = v1.y * v2.z - v1.z * v2.y;
    ctype const c5 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c6 = v1.x * v2.y - v1.y * v2.x;
    ctype const c7 = v1.x * v2.u - v1.u * v2.x;
    ctype const c8 = v1.y * v2.u - v1.u * v2.y;
    ctype const c9 = v1.z * v2.u - v1.u * v2.z;
    ctype const c10 = v1.w * v2.u - v1.u * v2.w;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc gpr :: gpr(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator*(Vec3dc<T> const& v, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    ctype const c4 = v.u * ctype(s);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> operator*(Scalar3dc<T> s, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    ctype const c4 = ctype(s) * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc gpr :: gpr(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> operator*(Scalar3dc<T> s1, Scalar3dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
//
// rgpr(ul, ur) = cmpl(gpr(cmpl(ul),cmpl(ur)))
//
// The identity element of rgpr is the pseudoscalar: rgpr(X, I_3dc) = X.
////////////////////////////////////////////////////////////////////////////////

// cga3dc rgpr :: rgpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& A, MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c31 + A.c1 * B.c26 + A.c2 * B.c27 + A.c3 * B.c28 +
                     A.c4 * B.c29 + A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                     A.c8 * B.c21 - A.c9 * B.c16 - A.c10 * B.c17 - A.c11 * B.c18 -
                     A.c12 * B.c22 - A.c13 * B.c23 - A.c14 * B.c24 - A.c15 * B.c25 -
                     A.c16 * B.c9 - A.c17 * B.c10 - A.c18 * B.c11 - A.c19 * B.c6 -
                     A.c20 * B.c7 - A.c21 * B.c8 - A.c22 * B.c12 - A.c23 * B.c13 -
                     A.c24 * B.c14 - A.c25 * B.c15 + A.c26 * B.c1 + A.c27 * B.c2 +
                     A.c28 * B.c3 + A.c29 * B.c4 + A.c30 * B.c5 + A.c31 * B.c0;
    ctype const c1 = A.c0 * B.c26 + A.c1 * B.c31 + A.c2 * B.c18 - A.c3 * B.c17 -
                     A.c4 * B.c19 + A.c5 * B.c22 + A.c6 * B.c29 + A.c7 * B.c14 -
                     A.c8 * B.c13 + A.c9 * B.c15 + A.c10 * B.c28 - A.c11 * B.c27 -
                     A.c12 * B.c30 - A.c13 * B.c8 + A.c14 * B.c7 + A.c15 * B.c9 +
                     A.c16 * B.c25 + A.c17 * B.c3 - A.c18 * B.c2 + A.c19 * B.c4 +
                     A.c20 * B.c24 - A.c21 * B.c23 - A.c22 * B.c5 - A.c23 * B.c21 +
                     A.c24 * B.c20 + A.c25 * B.c16 + A.c26 * B.c0 + A.c27 * B.c11 -
                     A.c28 * B.c10 - A.c29 * B.c6 + A.c30 * B.c12 + A.c31 * B.c1;
    ctype const c2 = A.c0 * B.c27 - A.c1 * B.c18 + A.c2 * B.c31 + A.c3 * B.c16 -
                     A.c4 * B.c20 + A.c5 * B.c23 - A.c6 * B.c14 + A.c7 * B.c29 +
                     A.c8 * B.c12 - A.c9 * B.c28 + A.c10 * B.c15 + A.c11 * B.c26 +
                     A.c12 * B.c8 - A.c13 * B.c30 - A.c14 * B.c6 + A.c15 * B.c10 -
                     A.c16 * B.c3 + A.c17 * B.c25 + A.c18 * B.c1 - A.c19 * B.c24 +
                     A.c20 * B.c4 + A.c21 * B.c22 + A.c22 * B.c21 - A.c23 * B.c5 -
                     A.c24 * B.c19 + A.c25 * B.c17 - A.c26 * B.c11 + A.c27 * B.c0 +
                     A.c28 * B.c9 - A.c29 * B.c7 + A.c30 * B.c13 + A.c31 * B.c2;
    ctype const c3 = A.c0 * B.c28 + A.c1 * B.c17 - A.c2 * B.c16 + A.c3 * B.c31 -
                     A.c4 * B.c21 + A.c5 * B.c24 + A.c6 * B.c13 - A.c7 * B.c12 +
                     A.c8 * B.c29 + A.c9 * B.c27 - A.c10 * B.c26 + A.c11 * B.c15 -
                     A.c12 * B.c7 + A.c13 * B.c6 - A.c14 * B.c30 + A.c15 * B.c11 +
                     A.c16 * B.c2 - A.c17 * B.c1 + A.c18 * B.c25 + A.c19 * B.c23 -
                     A.c20 * B.c22 + A.c21 * B.c4 - A.c22 * B.c20 + A.c23 * B.c19 -
                     A.c24 * B.c5 + A.c25 * B.c18 + A.c26 * B.c10 - A.c27 * B.c9 +
                     A.c28 * B.c0 - A.c29 * B.c8 + A.c30 * B.c14 + A.c31 * B.c3;
    ctype const c4 = -A.c0 * B.c30 + A.c1 * B.c22 + A.c2 * B.c23 + A.c3 * B.c24 +
                     A.c4 * B.c25 + A.c4 * B.c31 + A.c6 * B.c9 - A.c6 * B.c26 +
                     A.c7 * B.c10 - A.c7 * B.c27 + A.c8 * B.c11 - A.c8 * B.c28 +
                     A.c9 * B.c6 + A.c10 * B.c7 + A.c11 * B.c8 - A.c15 * B.c30 -
                     A.c16 * B.c22 - A.c17 * B.c23 - A.c18 * B.c24 - A.c22 * B.c1 -
                     A.c22 * B.c16 - A.c23 * B.c2 - A.c23 * B.c17 - A.c24 * B.c3 -
                     A.c24 * B.c18 - A.c25 * B.c4 + A.c26 * B.c6 + A.c27 * B.c7 +
                     A.c28 * B.c8 - A.c30 * B.c0 + A.c30 * B.c15 + A.c31 * B.c4;
    ctype const c5 = -A.c0 * B.c29 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                     A.c5 * B.c25 + A.c5 * B.c31 + A.c9 * B.c12 + A.c10 * B.c13 +
                     A.c11 * B.c14 + A.c12 * B.c9 + A.c12 * B.c26 + A.c13 * B.c10 +
                     A.c13 * B.c27 + A.c14 * B.c11 + A.c14 * B.c28 + A.c15 * B.c29 -
                     A.c16 * B.c19 - A.c17 * B.c20 - A.c18 * B.c21 + A.c19 * B.c1 -
                     A.c19 * B.c16 + A.c20 * B.c2 - A.c20 * B.c17 + A.c21 * B.c3 -
                     A.c21 * B.c18 + A.c25 * B.c5 - A.c26 * B.c12 - A.c27 * B.c13 -
                     A.c28 * B.c14 - A.c29 * B.c0 - A.c29 * B.c15 + A.c31 * B.c5;
    ctype const c6 = A.c0 * B.c22 - A.c1 * B.c30 - A.c2 * B.c8 + A.c3 * B.c7 +
                     A.c4 * B.c9 - A.c4 * B.c26 + A.c6 * B.c25 + A.c6 * B.c31 +
                     A.c7 * B.c3 + A.c7 * B.c18 - A.c8 * B.c2 - A.c8 * B.c17 +
                     A.c9 * B.c4 + A.c10 * B.c24 - A.c11 * B.c23 + A.c15 * B.c22 +
                     A.c16 * B.c30 + A.c17 * B.c8 - A.c18 * B.c7 + A.c22 * B.c0 -
                     A.c22 * B.c15 + A.c23 * B.c11 - A.c23 * B.c28 - A.c24 * B.c10 +
                     A.c24 * B.c27 - A.c25 * B.c6 + A.c26 * B.c4 + A.c27 * B.c24 -
                     A.c28 * B.c23 + A.c30 * B.c1 + A.c30 * B.c16 + A.c31 * B.c6;
    ctype const c7 = A.c0 * B.c23 + A.c1 * B.c8 - A.c2 * B.c30 - A.c3 * B.c6 +
                     A.c4 * B.c10 - A.c4 * B.c27 - A.c6 * B.c3 - A.c6 * B.c18 +
                     A.c7 * B.c25 + A.c7 * B.c31 + A.c8 * B.c1 + A.c8 * B.c16 -
                     A.c9 * B.c24 + A.c10 * B.c4 + A.c11 * B.c22 + A.c15 * B.c23 -
                     A.c16 * B.c8 + A.c17 * B.c30 + A.c18 * B.c6 - A.c22 * B.c11 +
                     A.c22 * B.c28 + A.c23 * B.c0 - A.c23 * B.c15 + A.c24 * B.c9 -
                     A.c24 * B.c26 - A.c25 * B.c7 - A.c26 * B.c24 + A.c27 * B.c4 +
                     A.c28 * B.c22 + A.c30 * B.c2 + A.c30 * B.c17 + A.c31 * B.c7;
    ctype const c8 = A.c0 * B.c24 - A.c1 * B.c7 + A.c2 * B.c6 - A.c3 * B.c30 +
                     A.c4 * B.c11 - A.c4 * B.c28 + A.c6 * B.c2 + A.c6 * B.c17 -
                     A.c7 * B.c1 - A.c7 * B.c16 + A.c8 * B.c25 + A.c8 * B.c31 +
                     A.c9 * B.c23 - A.c10 * B.c22 + A.c11 * B.c4 + A.c15 * B.c24 +
                     A.c16 * B.c7 - A.c17 * B.c6 + A.c18 * B.c30 + A.c22 * B.c10 -
                     A.c22 * B.c27 - A.c23 * B.c9 + A.c23 * B.c26 + A.c24 * B.c0 -
                     A.c24 * B.c15 - A.c25 * B.c8 + A.c26 * B.c23 - A.c27 * B.c22 +
                     A.c28 * B.c4 + A.c30 * B.c3 + A.c30 * B.c18 + A.c31 * B.c8;
    ctype const c9 = A.c0 * B.c16 - A.c1 * B.c15 - A.c2 * B.c28 + A.c3 * B.c27 +
                     A.c4 * B.c12 + A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c21 -
                     A.c8 * B.c20 + A.c9 * B.c31 + A.c10 * B.c18 - A.c11 * B.c17 +
                     A.c12 * B.c4 + A.c13 * B.c24 - A.c14 * B.c23 - A.c15 * B.c1 +
                     A.c16 * B.c0 + A.c17 * B.c11 - A.c18 * B.c10 + A.c19 * B.c30 +
                     A.c20 * B.c8 - A.c21 * B.c7 + A.c22 * B.c29 + A.c23 * B.c14 -
                     A.c24 * B.c13 - A.c25 * B.c26 - A.c26 * B.c25 - A.c27 * B.c3 +
                     A.c28 * B.c2 + A.c29 * B.c22 + A.c30 * B.c19 + A.c31 * B.c9;
    ctype const c10 = A.c0 * B.c17 + A.c1 * B.c28 - A.c2 * B.c15 - A.c3 * B.c26 +
                      A.c4 * B.c13 + A.c5 * B.c7 - A.c6 * B.c21 + A.c7 * B.c5 +
                      A.c8 * B.c19 - A.c9 * B.c18 + A.c10 * B.c31 + A.c11 * B.c16 -
                      A.c12 * B.c24 + A.c13 * B.c4 + A.c14 * B.c22 - A.c15 * B.c2 -
                      A.c16 * B.c11 + A.c17 * B.c0 + A.c18 * B.c9 - A.c19 * B.c8 +
                      A.c20 * B.c30 + A.c21 * B.c6 - A.c22 * B.c14 + A.c23 * B.c29 +
                      A.c24 * B.c12 - A.c25 * B.c27 + A.c26 * B.c3 - A.c27 * B.c25 -
                      A.c28 * B.c1 + A.c29 * B.c23 + A.c30 * B.c20 + A.c31 * B.c10;
    ctype const c11 = A.c0 * B.c18 - A.c1 * B.c27 + A.c2 * B.c26 - A.c3 * B.c15 +
                      A.c4 * B.c14 + A.c5 * B.c8 + A.c6 * B.c20 - A.c7 * B.c19 +
                      A.c8 * B.c5 + A.c9 * B.c17 - A.c10 * B.c16 + A.c11 * B.c31 +
                      A.c12 * B.c23 - A.c13 * B.c22 + A.c14 * B.c4 - A.c15 * B.c3 +
                      A.c16 * B.c10 - A.c17 * B.c9 + A.c18 * B.c0 + A.c19 * B.c7 -
                      A.c20 * B.c6 + A.c21 * B.c30 + A.c22 * B.c13 - A.c23 * B.c12 +
                      A.c24 * B.c29 - A.c25 * B.c28 - A.c26 * B.c2 + A.c27 * B.c1 -
                      A.c28 * B.c25 + A.c29 * B.c24 + A.c30 * B.c21 + A.c31 * B.c11;
    ctype const c12 = A.c0 * B.c19 + A.c1 * B.c29 + A.c2 * B.c14 - A.c3 * B.c13 +
                      A.c5 * B.c9 + A.c5 * B.c26 + A.c9 * B.c5 + A.c10 * B.c21 -
                      A.c11 * B.c20 - A.c12 * B.c25 + A.c12 * B.c31 - A.c13 * B.c3 +
                      A.c13 * B.c18 + A.c14 * B.c2 - A.c14 * B.c17 - A.c15 * B.c19 +
                      A.c16 * B.c29 + A.c17 * B.c14 - A.c18 * B.c13 + A.c19 * B.c0 +
                      A.c19 * B.c15 + A.c20 * B.c11 + A.c20 * B.c28 - A.c21 * B.c10 -
                      A.c21 * B.c27 + A.c25 * B.c12 - A.c26 * B.c5 - A.c27 * B.c21 +
                      A.c28 * B.c20 - A.c29 * B.c1 + A.c29 * B.c16 + A.c31 * B.c12;
    ctype const c13 = A.c0 * B.c20 - A.c1 * B.c14 + A.c2 * B.c29 + A.c3 * B.c12 +
                      A.c5 * B.c10 + A.c5 * B.c27 - A.c9 * B.c21 + A.c10 * B.c5 +
                      A.c11 * B.c19 + A.c12 * B.c3 - A.c12 * B.c18 - A.c13 * B.c25 +
                      A.c13 * B.c31 - A.c14 * B.c1 + A.c14 * B.c16 - A.c15 * B.c20 -
                      A.c16 * B.c14 + A.c17 * B.c29 + A.c18 * B.c12 - A.c19 * B.c11 -
                      A.c19 * B.c28 + A.c20 * B.c0 + A.c20 * B.c15 + A.c21 * B.c9 +
                      A.c21 * B.c26 + A.c25 * B.c13 + A.c26 * B.c21 - A.c27 * B.c5 -
                      A.c28 * B.c19 - A.c29 * B.c2 + A.c29 * B.c17 + A.c31 * B.c13;
    ctype const c14 = A.c0 * B.c21 + A.c1 * B.c13 - A.c2 * B.c12 + A.c3 * B.c29 +
                      A.c5 * B.c11 + A.c5 * B.c28 + A.c9 * B.c20 - A.c10 * B.c19 +
                      A.c11 * B.c5 - A.c12 * B.c2 + A.c12 * B.c17 + A.c13 * B.c1 -
                      A.c13 * B.c16 - A.c14 * B.c25 + A.c14 * B.c31 - A.c15 * B.c21 +
                      A.c16 * B.c13 - A.c17 * B.c12 + A.c18 * B.c29 + A.c19 * B.c10 +
                      A.c19 * B.c27 - A.c20 * B.c9 - A.c20 * B.c26 + A.c21 * B.c0 +
                      A.c21 * B.c15 + A.c25 * B.c14 - A.c26 * B.c20 + A.c27 * B.c19 -
                      A.c28 * B.c5 - A.c29 * B.c3 + A.c29 * B.c18 + A.c31 * B.c14;
    ctype const c15 = -A.c0 * B.c25 + A.c1 * B.c9 + A.c2 * B.c10 + A.c3 * B.c11 +
                      A.c4 * B.c29 - A.c5 * B.c30 - A.c6 * B.c19 - A.c7 * B.c20 -
                      A.c8 * B.c21 + A.c9 * B.c1 + A.c10 * B.c2 + A.c11 * B.c3 +
                      A.c12 * B.c22 + A.c13 * B.c23 + A.c14 * B.c24 + A.c15 * B.c31 -
                      A.c16 * B.c26 - A.c17 * B.c27 - A.c18 * B.c28 + A.c19 * B.c6 +
                      A.c20 * B.c7 + A.c21 * B.c8 - A.c22 * B.c12 - A.c23 * B.c13 -
                      A.c24 * B.c14 - A.c25 * B.c0 - A.c26 * B.c16 - A.c27 * B.c17 -
                      A.c28 * B.c18 - A.c29 * B.c4 + A.c30 * B.c5 + A.c31 * B.c15;
    ctype const c16 = -A.c0 * B.c9 + A.c1 * B.c25 + A.c2 * B.c3 - A.c3 * B.c2 +
                      A.c4 * B.c19 + A.c5 * B.c22 - A.c6 * B.c29 - A.c7 * B.c14 +
                      A.c8 * B.c13 - A.c9 * B.c0 - A.c10 * B.c11 + A.c11 * B.c10 -
                      A.c12 * B.c30 - A.c13 * B.c8 + A.c14 * B.c7 - A.c15 * B.c26 +
                      A.c16 * B.c31 + A.c17 * B.c18 - A.c18 * B.c17 + A.c19 * B.c4 +
                      A.c20 * B.c24 - A.c21 * B.c23 + A.c22 * B.c5 + A.c23 * B.c21 -
                      A.c24 * B.c20 + A.c25 * B.c1 - A.c26 * B.c15 - A.c27 * B.c28 +
                      A.c28 * B.c27 - A.c29 * B.c6 - A.c30 * B.c12 + A.c31 * B.c16;
    ctype const c17 = -A.c0 * B.c10 - A.c1 * B.c3 + A.c2 * B.c25 + A.c3 * B.c1 +
                      A.c4 * B.c20 + A.c5 * B.c23 + A.c6 * B.c14 - A.c7 * B.c29 -
                      A.c8 * B.c12 + A.c9 * B.c11 - A.c10 * B.c0 - A.c11 * B.c9 +
                      A.c12 * B.c8 - A.c13 * B.c30 - A.c14 * B.c6 - A.c15 * B.c27 -
                      A.c16 * B.c18 + A.c17 * B.c31 + A.c18 * B.c16 - A.c19 * B.c24 +
                      A.c20 * B.c4 + A.c21 * B.c22 - A.c22 * B.c21 + A.c23 * B.c5 +
                      A.c24 * B.c19 + A.c25 * B.c2 + A.c26 * B.c28 - A.c27 * B.c15 -
                      A.c28 * B.c26 - A.c29 * B.c7 - A.c30 * B.c13 + A.c31 * B.c17;
    ctype const c18 = -A.c0 * B.c11 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c25 +
                      A.c4 * B.c21 + A.c5 * B.c24 - A.c6 * B.c13 + A.c7 * B.c12 -
                      A.c8 * B.c29 - A.c9 * B.c10 + A.c10 * B.c9 - A.c11 * B.c0 -
                      A.c12 * B.c7 + A.c13 * B.c6 - A.c14 * B.c30 - A.c15 * B.c28 +
                      A.c16 * B.c17 - A.c17 * B.c16 + A.c18 * B.c31 + A.c19 * B.c23 -
                      A.c20 * B.c22 + A.c21 * B.c4 + A.c22 * B.c20 - A.c23 * B.c19 +
                      A.c24 * B.c5 + A.c25 * B.c3 - A.c26 * B.c27 + A.c27 * B.c26 -
                      A.c28 * B.c15 - A.c29 * B.c8 - A.c30 * B.c14 + A.c31 * B.c18;
    ctype const c19 = -A.c0 * B.c12 + A.c1 * B.c5 + A.c2 * B.c21 - A.c3 * B.c20 -
                      A.c5 * B.c1 + A.c5 * B.c16 - A.c9 * B.c29 - A.c10 * B.c14 +
                      A.c11 * B.c13 - A.c12 * B.c0 - A.c12 * B.c15 - A.c13 * B.c11 -
                      A.c13 * B.c28 + A.c14 * B.c10 + A.c14 * B.c27 + A.c15 * B.c12 +
                      A.c16 * B.c5 + A.c17 * B.c21 - A.c18 * B.c20 - A.c19 * B.c25 +
                      A.c19 * B.c31 - A.c20 * B.c3 + A.c20 * B.c18 + A.c21 * B.c2 -
                      A.c21 * B.c17 + A.c25 * B.c19 + A.c26 * B.c29 + A.c27 * B.c14 -
                      A.c28 * B.c13 - A.c29 * B.c9 - A.c29 * B.c26 + A.c31 * B.c19;
    ctype const c20 = -A.c0 * B.c13 - A.c1 * B.c21 + A.c2 * B.c5 + A.c3 * B.c19 -
                      A.c5 * B.c2 + A.c5 * B.c17 + A.c9 * B.c14 - A.c10 * B.c29 -
                      A.c11 * B.c12 + A.c12 * B.c11 + A.c12 * B.c28 - A.c13 * B.c0 -
                      A.c13 * B.c15 - A.c14 * B.c9 - A.c14 * B.c26 + A.c15 * B.c13 -
                      A.c16 * B.c21 + A.c17 * B.c5 + A.c18 * B.c19 + A.c19 * B.c3 -
                      A.c19 * B.c18 - A.c20 * B.c25 + A.c20 * B.c31 - A.c21 * B.c1 +
                      A.c21 * B.c16 + A.c25 * B.c20 - A.c26 * B.c14 + A.c27 * B.c29 +
                      A.c28 * B.c12 - A.c29 * B.c10 - A.c29 * B.c27 + A.c31 * B.c20;
    ctype const c21 = -A.c0 * B.c14 + A.c1 * B.c20 - A.c2 * B.c19 + A.c3 * B.c5 -
                      A.c5 * B.c3 + A.c5 * B.c18 - A.c9 * B.c13 + A.c10 * B.c12 -
                      A.c11 * B.c29 - A.c12 * B.c10 - A.c12 * B.c27 + A.c13 * B.c9 +
                      A.c13 * B.c26 - A.c14 * B.c0 - A.c14 * B.c15 + A.c15 * B.c14 +
                      A.c16 * B.c20 - A.c17 * B.c19 + A.c18 * B.c5 - A.c19 * B.c2 +
                      A.c19 * B.c17 + A.c20 * B.c1 - A.c20 * B.c16 - A.c21 * B.c25 +
                      A.c21 * B.c31 + A.c25 * B.c21 + A.c26 * B.c13 - A.c27 * B.c12 +
                      A.c28 * B.c29 - A.c29 * B.c11 - A.c29 * B.c28 + A.c31 * B.c21;
    ctype const c22 = -A.c0 * B.c6 - A.c1 * B.c4 - A.c2 * B.c24 + A.c3 * B.c23 +
                      A.c4 * B.c1 + A.c4 * B.c16 - A.c6 * B.c0 + A.c6 * B.c15 -
                      A.c7 * B.c11 + A.c7 * B.c28 + A.c8 * B.c10 - A.c8 * B.c27 -
                      A.c9 * B.c30 - A.c10 * B.c8 + A.c11 * B.c7 - A.c15 * B.c6 +
                      A.c16 * B.c4 + A.c17 * B.c24 - A.c18 * B.c23 + A.c22 * B.c25 +
                      A.c22 * B.c31 + A.c23 * B.c3 + A.c23 * B.c18 - A.c24 * B.c2 -
                      A.c24 * B.c17 - A.c25 * B.c22 - A.c26 * B.c30 - A.c27 * B.c8 +
                      A.c28 * B.c7 - A.c30 * B.c9 + A.c30 * B.c26 + A.c31 * B.c22;
    ctype const c23 = -A.c0 * B.c7 + A.c1 * B.c24 - A.c2 * B.c4 - A.c3 * B.c22 +
                      A.c4 * B.c2 + A.c4 * B.c17 + A.c6 * B.c11 - A.c6 * B.c28 -
                      A.c7 * B.c0 + A.c7 * B.c15 - A.c8 * B.c9 + A.c8 * B.c26 +
                      A.c9 * B.c8 - A.c10 * B.c30 - A.c11 * B.c6 - A.c15 * B.c7 -
                      A.c16 * B.c24 + A.c17 * B.c4 + A.c18 * B.c22 - A.c22 * B.c3 -
                      A.c22 * B.c18 + A.c23 * B.c25 + A.c23 * B.c31 + A.c24 * B.c1 +
                      A.c24 * B.c16 - A.c25 * B.c23 + A.c26 * B.c8 - A.c27 * B.c30 -
                      A.c28 * B.c6 - A.c30 * B.c10 + A.c30 * B.c27 + A.c31 * B.c23;
    ctype const c24 = -A.c0 * B.c8 - A.c1 * B.c23 + A.c2 * B.c22 - A.c3 * B.c4 +
                      A.c4 * B.c3 + A.c4 * B.c18 - A.c6 * B.c10 + A.c6 * B.c27 +
                      A.c7 * B.c9 - A.c7 * B.c26 - A.c8 * B.c0 + A.c8 * B.c15 -
                      A.c9 * B.c7 + A.c10 * B.c6 - A.c11 * B.c30 - A.c15 * B.c8 +
                      A.c16 * B.c23 - A.c17 * B.c22 + A.c18 * B.c4 + A.c22 * B.c2 +
                      A.c22 * B.c17 - A.c23 * B.c1 - A.c23 * B.c16 + A.c24 * B.c25 +
                      A.c24 * B.c31 - A.c25 * B.c24 - A.c26 * B.c7 + A.c27 * B.c6 -
                      A.c28 * B.c30 - A.c30 * B.c11 + A.c30 * B.c28 + A.c31 * B.c24;
    ctype const c25 = A.c0 * B.c15 - A.c1 * B.c16 - A.c2 * B.c17 - A.c3 * B.c18 -
                      A.c4 * B.c5 + A.c5 * B.c4 - A.c6 * B.c12 - A.c7 * B.c13 -
                      A.c8 * B.c14 - A.c9 * B.c26 - A.c10 * B.c27 - A.c11 * B.c28 +
                      A.c12 * B.c6 + A.c13 * B.c7 + A.c14 * B.c8 + A.c15 * B.c0 -
                      A.c16 * B.c1 - A.c17 * B.c2 - A.c18 * B.c3 - A.c19 * B.c22 -
                      A.c20 * B.c23 - A.c21 * B.c24 + A.c22 * B.c19 + A.c23 * B.c20 +
                      A.c24 * B.c21 + A.c25 * B.c31 - A.c26 * B.c9 - A.c27 * B.c10 -
                      A.c28 * B.c11 - A.c29 * B.c30 + A.c30 * B.c29 + A.c31 * B.c25;
    ctype const c26 = -A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c11 + A.c3 * B.c10 +
                      A.c4 * B.c12 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c21 -
                      A.c8 * B.c20 - A.c9 * B.c25 - A.c10 * B.c3 + A.c11 * B.c2 -
                      A.c12 * B.c4 - A.c13 * B.c24 + A.c14 * B.c23 + A.c15 * B.c16 +
                      A.c16 * B.c15 + A.c17 * B.c28 - A.c18 * B.c27 - A.c19 * B.c30 -
                      A.c20 * B.c8 + A.c21 * B.c7 + A.c22 * B.c29 + A.c23 * B.c14 -
                      A.c24 * B.c13 - A.c25 * B.c9 + A.c26 * B.c31 + A.c27 * B.c18 -
                      A.c28 * B.c17 - A.c29 * B.c22 + A.c30 * B.c19 + A.c31 * B.c26;
    ctype const c27 = -A.c0 * B.c2 + A.c1 * B.c11 - A.c2 * B.c0 - A.c3 * B.c9 +
                      A.c4 * B.c13 - A.c5 * B.c7 - A.c6 * B.c21 + A.c7 * B.c5 +
                      A.c8 * B.c19 + A.c9 * B.c3 - A.c10 * B.c25 - A.c11 * B.c1 +
                      A.c12 * B.c24 - A.c13 * B.c4 - A.c14 * B.c22 + A.c15 * B.c17 -
                      A.c16 * B.c28 + A.c17 * B.c15 + A.c18 * B.c26 + A.c19 * B.c8 -
                      A.c20 * B.c30 - A.c21 * B.c6 - A.c22 * B.c14 + A.c23 * B.c29 +
                      A.c24 * B.c12 - A.c25 * B.c10 - A.c26 * B.c18 + A.c27 * B.c31 +
                      A.c28 * B.c16 - A.c29 * B.c23 + A.c30 * B.c20 + A.c31 * B.c27;
    ctype const c28 = -A.c0 * B.c3 - A.c1 * B.c10 + A.c2 * B.c9 - A.c3 * B.c0 +
                      A.c4 * B.c14 - A.c5 * B.c8 + A.c6 * B.c20 - A.c7 * B.c19 +
                      A.c8 * B.c5 - A.c9 * B.c2 + A.c10 * B.c1 - A.c11 * B.c25 -
                      A.c12 * B.c23 + A.c13 * B.c22 - A.c14 * B.c4 + A.c15 * B.c18 +
                      A.c16 * B.c27 - A.c17 * B.c26 + A.c18 * B.c15 - A.c19 * B.c7 +
                      A.c20 * B.c6 - A.c21 * B.c30 + A.c22 * B.c13 - A.c23 * B.c12 +
                      A.c24 * B.c29 - A.c25 * B.c11 + A.c26 * B.c17 - A.c27 * B.c16 +
                      A.c28 * B.c31 - A.c29 * B.c24 + A.c30 * B.c21 + A.c31 * B.c28;
    ctype const c29 = A.c0 * B.c5 - A.c1 * B.c12 - A.c2 * B.c13 - A.c3 * B.c14 +
                      A.c5 * B.c0 + A.c5 * B.c15 - A.c9 * B.c19 - A.c10 * B.c20 -
                      A.c11 * B.c21 + A.c12 * B.c1 - A.c12 * B.c16 + A.c13 * B.c2 -
                      A.c13 * B.c17 + A.c14 * B.c3 - A.c14 * B.c18 - A.c15 * B.c5 -
                      A.c16 * B.c12 - A.c17 * B.c13 - A.c18 * B.c14 - A.c19 * B.c9 -
                      A.c19 * B.c26 - A.c20 * B.c10 - A.c20 * B.c27 - A.c21 * B.c11 -
                      A.c21 * B.c28 + A.c25 * B.c29 + A.c26 * B.c19 + A.c27 * B.c20 +
                      A.c28 * B.c21 - A.c29 * B.c25 + A.c29 * B.c31 + A.c31 * B.c29;
    ctype const c30 = A.c0 * B.c4 + A.c1 * B.c6 + A.c2 * B.c7 + A.c3 * B.c8 +
                      A.c4 * B.c0 - A.c4 * B.c15 - A.c6 * B.c1 - A.c6 * B.c16 -
                      A.c7 * B.c2 - A.c7 * B.c17 - A.c8 * B.c3 - A.c8 * B.c18 -
                      A.c9 * B.c22 - A.c10 * B.c23 - A.c11 * B.c24 + A.c15 * B.c4 -
                      A.c16 * B.c6 - A.c17 * B.c7 - A.c18 * B.c8 - A.c22 * B.c9 +
                      A.c22 * B.c26 - A.c23 * B.c10 + A.c23 * B.c27 - A.c24 * B.c11 +
                      A.c24 * B.c28 - A.c25 * B.c30 - A.c26 * B.c22 - A.c27 * B.c23 -
                      A.c28 * B.c24 + A.c30 * B.c25 + A.c30 * B.c31 + A.c31 * B.c30;
    ctype const c31 = -A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 +
                      A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c12 + A.c7 * B.c13 +
                      A.c8 * B.c14 + A.c9 * B.c9 + A.c10 * B.c10 + A.c11 * B.c11 +
                      A.c12 * B.c6 + A.c13 * B.c7 + A.c14 * B.c8 - A.c15 * B.c15 -
                      A.c16 * B.c16 - A.c17 * B.c17 - A.c18 * B.c18 - A.c19 * B.c22 -
                      A.c20 * B.c23 - A.c21 * B.c24 - A.c22 * B.c19 - A.c23 * B.c20 -
                      A.c24 * B.c21 + A.c25 * B.c25 + A.c26 * B.c26 + A.c27 * B.c27 +
                      A.c28 * B.c28 - A.c29 * B.c30 - A.c30 * B.c29 + A.c31 * B.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& A,
                                                 MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 +
                     A.c5 * B.c15 - A.c16 * B.c4 - A.c17 * B.c5 - A.c18 * B.c6 -
                     A.c19 * B.c1 - A.c20 * B.c2 - A.c21 * B.c3 - A.c22 * B.c7 -
                     A.c23 * B.c8 - A.c24 * B.c9 - A.c25 * B.c10 + A.c31 * B.c0;
    ctype const c1 = A.c0 * B.c11 + A.c6 * B.c14 + A.c7 * B.c9 - A.c8 * B.c8 +
                     A.c9 * B.c10 + A.c10 * B.c13 - A.c11 * B.c12 - A.c12 * B.c15 -
                     A.c13 * B.c3 + A.c14 * B.c2 + A.c15 * B.c4 + A.c26 * B.c0 +
                     A.c27 * B.c6 - A.c28 * B.c5 - A.c29 * B.c1 + A.c30 * B.c7;
    ctype const c2 = A.c0 * B.c12 - A.c6 * B.c9 + A.c7 * B.c14 + A.c8 * B.c7 -
                     A.c9 * B.c13 + A.c10 * B.c10 + A.c11 * B.c11 + A.c12 * B.c3 -
                     A.c13 * B.c15 - A.c14 * B.c1 + A.c15 * B.c5 - A.c26 * B.c6 +
                     A.c27 * B.c0 + A.c28 * B.c4 - A.c29 * B.c2 + A.c30 * B.c8;
    ctype const c3 = A.c0 * B.c13 + A.c6 * B.c8 - A.c7 * B.c7 + A.c8 * B.c14 +
                     A.c9 * B.c12 - A.c10 * B.c11 + A.c11 * B.c10 - A.c12 * B.c2 +
                     A.c13 * B.c1 - A.c14 * B.c15 + A.c15 * B.c6 + A.c26 * B.c5 -
                     A.c27 * B.c4 + A.c28 * B.c0 - A.c29 * B.c3 + A.c30 * B.c9;
    ctype const c4 = -A.c0 * B.c15 + A.c6 * B.c4 - A.c6 * B.c11 + A.c7 * B.c5 -
                     A.c7 * B.c12 + A.c8 * B.c6 - A.c8 * B.c13 + A.c9 * B.c1 +
                     A.c10 * B.c2 + A.c11 * B.c3 - A.c15 * B.c15 + A.c26 * B.c1 +
                     A.c27 * B.c2 + A.c28 * B.c3 - A.c30 * B.c0 + A.c30 * B.c10;
    ctype const c5 = -A.c0 * B.c14 + A.c9 * B.c7 + A.c10 * B.c8 + A.c11 * B.c9 +
                     A.c12 * B.c4 + A.c12 * B.c11 + A.c13 * B.c5 + A.c13 * B.c12 +
                     A.c14 * B.c6 + A.c14 * B.c13 + A.c15 * B.c14 - A.c26 * B.c7 -
                     A.c27 * B.c8 - A.c28 * B.c9 - A.c29 * B.c0 - A.c29 * B.c10;
    ctype const c6 = -A.c1 * B.c15 - A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c4 -
                     A.c4 * B.c11 + A.c16 * B.c15 + A.c17 * B.c3 - A.c18 * B.c2 +
                     A.c22 * B.c0 - A.c22 * B.c10 + A.c23 * B.c6 - A.c23 * B.c13 -
                     A.c24 * B.c5 + A.c24 * B.c12 - A.c25 * B.c1 + A.c31 * B.c1;
    ctype const c7 = A.c1 * B.c3 - A.c2 * B.c15 - A.c3 * B.c1 + A.c4 * B.c5 -
                     A.c4 * B.c12 - A.c16 * B.c3 + A.c17 * B.c15 + A.c18 * B.c1 -
                     A.c22 * B.c6 + A.c22 * B.c13 + A.c23 * B.c0 - A.c23 * B.c10 +
                     A.c24 * B.c4 - A.c24 * B.c11 - A.c25 * B.c2 + A.c31 * B.c2;
    ctype const c8 = -A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c15 + A.c4 * B.c6 -
                     A.c4 * B.c13 + A.c16 * B.c2 - A.c17 * B.c1 + A.c18 * B.c15 +
                     A.c22 * B.c5 - A.c22 * B.c12 - A.c23 * B.c4 + A.c23 * B.c11 +
                     A.c24 * B.c0 - A.c24 * B.c10 - A.c25 * B.c3 + A.c31 * B.c3;
    ctype const c9 = -A.c1 * B.c10 - A.c2 * B.c13 + A.c3 * B.c12 + A.c4 * B.c7 +
                     A.c5 * B.c1 + A.c16 * B.c0 + A.c17 * B.c6 - A.c18 * B.c5 +
                     A.c19 * B.c15 + A.c20 * B.c3 - A.c21 * B.c2 + A.c22 * B.c14 +
                     A.c23 * B.c9 - A.c24 * B.c8 - A.c25 * B.c11 + A.c31 * B.c4;
    ctype const c10 = A.c1 * B.c13 - A.c2 * B.c10 - A.c3 * B.c11 + A.c4 * B.c8 +
                      A.c5 * B.c2 - A.c16 * B.c6 + A.c17 * B.c0 + A.c18 * B.c4 -
                      A.c19 * B.c3 + A.c20 * B.c15 + A.c21 * B.c1 - A.c22 * B.c9 +
                      A.c23 * B.c14 + A.c24 * B.c7 - A.c25 * B.c12 + A.c31 * B.c5;
    ctype const c11 = -A.c1 * B.c12 + A.c2 * B.c11 - A.c3 * B.c10 + A.c4 * B.c9 +
                      A.c5 * B.c3 + A.c16 * B.c5 - A.c17 * B.c4 + A.c18 * B.c0 +
                      A.c19 * B.c2 - A.c20 * B.c1 + A.c21 * B.c15 + A.c22 * B.c8 -
                      A.c23 * B.c7 + A.c24 * B.c14 - A.c25 * B.c13 + A.c31 * B.c6;
    ctype const c12 = A.c1 * B.c14 + A.c2 * B.c9 - A.c3 * B.c8 + A.c5 * B.c4 +
                      A.c5 * B.c11 + A.c16 * B.c14 + A.c17 * B.c9 - A.c18 * B.c8 +
                      A.c19 * B.c0 + A.c19 * B.c10 + A.c20 * B.c6 + A.c20 * B.c13 -
                      A.c21 * B.c5 - A.c21 * B.c12 + A.c25 * B.c7 + A.c31 * B.c7;
    ctype const c13 = -A.c1 * B.c9 + A.c2 * B.c14 + A.c3 * B.c7 + A.c5 * B.c5 +
                      A.c5 * B.c12 - A.c16 * B.c9 + A.c17 * B.c14 + A.c18 * B.c7 -
                      A.c19 * B.c6 - A.c19 * B.c13 + A.c20 * B.c0 + A.c20 * B.c10 +
                      A.c21 * B.c4 + A.c21 * B.c11 + A.c25 * B.c8 + A.c31 * B.c8;
    ctype const c14 = A.c1 * B.c8 - A.c2 * B.c7 + A.c3 * B.c14 + A.c5 * B.c6 +
                      A.c5 * B.c13 + A.c16 * B.c8 - A.c17 * B.c7 + A.c18 * B.c14 +
                      A.c19 * B.c5 + A.c19 * B.c12 - A.c20 * B.c4 - A.c20 * B.c11 +
                      A.c21 * B.c0 + A.c21 * B.c10 + A.c25 * B.c9 + A.c31 * B.c9;
    ctype const c15 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c14 -
                      A.c5 * B.c15 - A.c16 * B.c11 - A.c17 * B.c12 - A.c18 * B.c13 +
                      A.c19 * B.c1 + A.c20 * B.c2 + A.c21 * B.c3 - A.c22 * B.c7 -
                      A.c23 * B.c8 - A.c24 * B.c9 - A.c25 * B.c0 + A.c31 * B.c10;
    ctype const c16 = -A.c0 * B.c4 - A.c6 * B.c14 - A.c7 * B.c9 + A.c8 * B.c8 -
                      A.c9 * B.c0 - A.c10 * B.c6 + A.c11 * B.c5 - A.c12 * B.c15 -
                      A.c13 * B.c3 + A.c14 * B.c2 - A.c15 * B.c11 - A.c26 * B.c10 -
                      A.c27 * B.c13 + A.c28 * B.c12 - A.c29 * B.c1 - A.c30 * B.c7;
    ctype const c17 = -A.c0 * B.c5 + A.c6 * B.c9 - A.c7 * B.c14 - A.c8 * B.c7 +
                      A.c9 * B.c6 - A.c10 * B.c0 - A.c11 * B.c4 + A.c12 * B.c3 -
                      A.c13 * B.c15 - A.c14 * B.c1 - A.c15 * B.c12 + A.c26 * B.c13 -
                      A.c27 * B.c10 - A.c28 * B.c11 - A.c29 * B.c2 - A.c30 * B.c8;
    ctype const c18 = -A.c0 * B.c6 - A.c6 * B.c8 + A.c7 * B.c7 - A.c8 * B.c14 -
                      A.c9 * B.c5 + A.c10 * B.c4 - A.c11 * B.c0 - A.c12 * B.c2 +
                      A.c13 * B.c1 - A.c14 * B.c15 - A.c15 * B.c13 - A.c26 * B.c12 +
                      A.c27 * B.c11 - A.c28 * B.c10 - A.c29 * B.c3 - A.c30 * B.c9;
    ctype const c19 = -A.c0 * B.c7 - A.c9 * B.c14 - A.c10 * B.c9 + A.c11 * B.c8 -
                      A.c12 * B.c0 - A.c12 * B.c10 - A.c13 * B.c6 - A.c13 * B.c13 +
                      A.c14 * B.c5 + A.c14 * B.c12 + A.c15 * B.c7 + A.c26 * B.c14 +
                      A.c27 * B.c9 - A.c28 * B.c8 - A.c29 * B.c4 - A.c29 * B.c11;
    ctype const c20 = -A.c0 * B.c8 + A.c9 * B.c9 - A.c10 * B.c14 - A.c11 * B.c7 +
                      A.c12 * B.c6 + A.c12 * B.c13 - A.c13 * B.c0 - A.c13 * B.c10 -
                      A.c14 * B.c4 - A.c14 * B.c11 + A.c15 * B.c8 - A.c26 * B.c9 +
                      A.c27 * B.c14 + A.c28 * B.c7 - A.c29 * B.c5 - A.c29 * B.c12;
    ctype const c21 = -A.c0 * B.c9 - A.c9 * B.c8 + A.c10 * B.c7 - A.c11 * B.c14 -
                      A.c12 * B.c5 - A.c12 * B.c12 + A.c13 * B.c4 + A.c13 * B.c11 -
                      A.c14 * B.c0 - A.c14 * B.c10 + A.c15 * B.c9 + A.c26 * B.c8 -
                      A.c27 * B.c7 + A.c28 * B.c14 - A.c29 * B.c6 - A.c29 * B.c13;
    ctype const c22 = -A.c0 * B.c1 - A.c6 * B.c0 + A.c6 * B.c10 - A.c7 * B.c6 +
                      A.c7 * B.c13 + A.c8 * B.c5 - A.c8 * B.c12 - A.c9 * B.c15 -
                      A.c10 * B.c3 + A.c11 * B.c2 - A.c15 * B.c1 - A.c26 * B.c15 -
                      A.c27 * B.c3 + A.c28 * B.c2 - A.c30 * B.c4 + A.c30 * B.c11;
    ctype const c23 = -A.c0 * B.c2 + A.c6 * B.c6 - A.c6 * B.c13 - A.c7 * B.c0 +
                      A.c7 * B.c10 - A.c8 * B.c4 + A.c8 * B.c11 + A.c9 * B.c3 -
                      A.c10 * B.c15 - A.c11 * B.c1 - A.c15 * B.c2 + A.c26 * B.c3 -
                      A.c27 * B.c15 - A.c28 * B.c1 - A.c30 * B.c5 + A.c30 * B.c12;
    ctype const c24 = -A.c0 * B.c3 - A.c6 * B.c5 + A.c6 * B.c12 + A.c7 * B.c4 -
                      A.c7 * B.c11 - A.c8 * B.c0 + A.c8 * B.c10 - A.c9 * B.c2 +
                      A.c10 * B.c1 - A.c11 * B.c15 - A.c15 * B.c3 - A.c26 * B.c2 +
                      A.c27 * B.c1 - A.c28 * B.c15 - A.c30 * B.c6 + A.c30 * B.c13;
    ctype const c25 = A.c0 * B.c10 - A.c6 * B.c7 - A.c7 * B.c8 - A.c8 * B.c9 -
                      A.c9 * B.c11 - A.c10 * B.c12 - A.c11 * B.c13 + A.c12 * B.c1 +
                      A.c13 * B.c2 + A.c14 * B.c3 + A.c15 * B.c0 - A.c26 * B.c4 -
                      A.c27 * B.c5 - A.c28 * B.c6 - A.c29 * B.c15 + A.c30 * B.c14;
    ctype const c26 = -A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c7 -
                      A.c5 * B.c1 + A.c16 * B.c10 + A.c17 * B.c13 - A.c18 * B.c12 -
                      A.c19 * B.c15 - A.c20 * B.c3 + A.c21 * B.c2 + A.c22 * B.c14 +
                      A.c23 * B.c9 - A.c24 * B.c8 - A.c25 * B.c4 + A.c31 * B.c11;
    ctype const c27 = A.c1 * B.c6 - A.c2 * B.c0 - A.c3 * B.c4 + A.c4 * B.c8 -
                      A.c5 * B.c2 - A.c16 * B.c13 + A.c17 * B.c10 + A.c18 * B.c11 +
                      A.c19 * B.c3 - A.c20 * B.c15 - A.c21 * B.c1 - A.c22 * B.c9 +
                      A.c23 * B.c14 + A.c24 * B.c7 - A.c25 * B.c5 + A.c31 * B.c12;
    ctype const c28 = -A.c1 * B.c5 + A.c2 * B.c4 - A.c3 * B.c0 + A.c4 * B.c9 -
                      A.c5 * B.c3 + A.c16 * B.c12 - A.c17 * B.c11 + A.c18 * B.c10 -
                      A.c19 * B.c2 + A.c20 * B.c1 - A.c21 * B.c15 + A.c22 * B.c8 -
                      A.c23 * B.c7 + A.c24 * B.c14 - A.c25 * B.c6 + A.c31 * B.c13;
    ctype const c29 = -A.c1 * B.c7 - A.c2 * B.c8 - A.c3 * B.c9 + A.c5 * B.c0 +
                      A.c5 * B.c10 - A.c16 * B.c7 - A.c17 * B.c8 - A.c18 * B.c9 -
                      A.c19 * B.c4 - A.c19 * B.c11 - A.c20 * B.c5 - A.c20 * B.c12 -
                      A.c21 * B.c6 - A.c21 * B.c13 + A.c25 * B.c14 + A.c31 * B.c14;
    ctype const c30 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c4 * B.c0 -
                      A.c4 * B.c10 - A.c16 * B.c1 - A.c17 * B.c2 - A.c18 * B.c3 -
                      A.c22 * B.c4 + A.c22 * B.c11 - A.c23 * B.c5 + A.c23 * B.c12 -
                      A.c24 * B.c6 + A.c24 * B.c13 - A.c25 * B.c15 + A.c31 * B.c15;
    ctype const c31 = -A.c0 * B.c0 + A.c6 * B.c7 + A.c7 * B.c8 + A.c8 * B.c9 +
                      A.c9 * B.c4 + A.c10 * B.c5 + A.c11 * B.c6 + A.c12 * B.c1 +
                      A.c13 * B.c2 + A.c14 * B.c3 - A.c15 * B.c10 + A.c26 * B.c11 +
                      A.c27 * B.c12 + A.c28 * B.c13 - A.c29 * B.c15 - A.c30 * B.c14;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& A,
                                                 MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c31 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 -
                     A.c4 * B.c16 - A.c5 * B.c17 - A.c6 * B.c18 - A.c7 * B.c22 -
                     A.c8 * B.c23 - A.c9 * B.c24 - A.c10 * B.c25 + A.c11 * B.c1 +
                     A.c12 * B.c2 + A.c13 * B.c3 + A.c14 * B.c4 + A.c15 * B.c5;
    ctype const c1 = A.c0 * B.c26 + A.c1 * B.c29 + A.c2 * B.c14 - A.c3 * B.c13 +
                     A.c4 * B.c15 + A.c5 * B.c28 - A.c6 * B.c27 - A.c7 * B.c30 -
                     A.c8 * B.c8 + A.c9 * B.c7 + A.c10 * B.c9 + A.c11 * B.c0 +
                     A.c12 * B.c11 - A.c13 * B.c10 - A.c14 * B.c6 + A.c15 * B.c12;
    ctype const c2 = A.c0 * B.c27 - A.c1 * B.c14 + A.c2 * B.c29 + A.c3 * B.c12 -
                     A.c4 * B.c28 + A.c5 * B.c15 + A.c6 * B.c26 + A.c7 * B.c8 -
                     A.c8 * B.c30 - A.c9 * B.c6 + A.c10 * B.c10 - A.c11 * B.c11 +
                     A.c12 * B.c0 + A.c13 * B.c9 - A.c14 * B.c7 + A.c15 * B.c13;
    ctype const c3 = A.c0 * B.c28 + A.c1 * B.c13 - A.c2 * B.c12 + A.c3 * B.c29 +
                     A.c4 * B.c27 - A.c5 * B.c26 + A.c6 * B.c15 - A.c7 * B.c7 +
                     A.c8 * B.c6 - A.c9 * B.c30 + A.c10 * B.c11 + A.c11 * B.c10 -
                     A.c12 * B.c9 + A.c13 * B.c0 - A.c14 * B.c8 + A.c15 * B.c14;
    ctype const c4 = -A.c0 * B.c30 + A.c1 * B.c9 - A.c1 * B.c26 + A.c2 * B.c10 -
                     A.c2 * B.c27 + A.c3 * B.c11 - A.c3 * B.c28 + A.c4 * B.c6 +
                     A.c5 * B.c7 + A.c6 * B.c8 - A.c10 * B.c30 + A.c11 * B.c6 +
                     A.c12 * B.c7 + A.c13 * B.c8 - A.c15 * B.c0 + A.c15 * B.c15;
    ctype const c5 = -A.c0 * B.c29 + A.c4 * B.c12 + A.c5 * B.c13 + A.c6 * B.c14 +
                     A.c7 * B.c9 + A.c7 * B.c26 + A.c8 * B.c10 + A.c8 * B.c27 +
                     A.c9 * B.c11 + A.c9 * B.c28 + A.c10 * B.c29 - A.c11 * B.c12 -
                     A.c12 * B.c13 - A.c13 * B.c14 - A.c14 * B.c0 - A.c14 * B.c15;
    ctype const c6 = A.c0 * B.c22 + A.c1 * B.c25 + A.c1 * B.c31 + A.c2 * B.c3 +
                     A.c2 * B.c18 - A.c3 * B.c2 - A.c3 * B.c17 + A.c4 * B.c4 +
                     A.c5 * B.c24 - A.c6 * B.c23 + A.c10 * B.c22 + A.c11 * B.c4 +
                     A.c12 * B.c24 - A.c13 * B.c23 + A.c15 * B.c1 + A.c15 * B.c16;
    ctype const c7 = A.c0 * B.c23 - A.c1 * B.c3 - A.c1 * B.c18 + A.c2 * B.c25 +
                     A.c2 * B.c31 + A.c3 * B.c1 + A.c3 * B.c16 - A.c4 * B.c24 +
                     A.c5 * B.c4 + A.c6 * B.c22 + A.c10 * B.c23 - A.c11 * B.c24 +
                     A.c12 * B.c4 + A.c13 * B.c22 + A.c15 * B.c2 + A.c15 * B.c17;
    ctype const c8 = A.c0 * B.c24 + A.c1 * B.c2 + A.c1 * B.c17 - A.c2 * B.c1 -
                     A.c2 * B.c16 + A.c3 * B.c25 + A.c3 * B.c31 + A.c4 * B.c23 -
                     A.c5 * B.c22 + A.c6 * B.c4 + A.c10 * B.c24 + A.c11 * B.c23 -
                     A.c12 * B.c22 + A.c13 * B.c4 + A.c15 * B.c3 + A.c15 * B.c18;
    ctype const c9 = A.c0 * B.c16 + A.c1 * B.c5 + A.c2 * B.c21 - A.c3 * B.c20 +
                     A.c4 * B.c31 + A.c5 * B.c18 - A.c6 * B.c17 + A.c7 * B.c4 +
                     A.c8 * B.c24 - A.c9 * B.c23 - A.c10 * B.c1 - A.c11 * B.c25 -
                     A.c12 * B.c3 + A.c13 * B.c2 + A.c14 * B.c22 + A.c15 * B.c19;
    ctype const c10 = A.c0 * B.c17 - A.c1 * B.c21 + A.c2 * B.c5 + A.c3 * B.c19 -
                      A.c4 * B.c18 + A.c5 * B.c31 + A.c6 * B.c16 - A.c7 * B.c24 +
                      A.c8 * B.c4 + A.c9 * B.c22 - A.c10 * B.c2 + A.c11 * B.c3 -
                      A.c12 * B.c25 - A.c13 * B.c1 + A.c14 * B.c23 + A.c15 * B.c20;
    ctype const c11 = A.c0 * B.c18 + A.c1 * B.c20 - A.c2 * B.c19 + A.c3 * B.c5 +
                      A.c4 * B.c17 - A.c5 * B.c16 + A.c6 * B.c31 + A.c7 * B.c23 -
                      A.c8 * B.c22 + A.c9 * B.c4 - A.c10 * B.c3 - A.c11 * B.c2 +
                      A.c12 * B.c1 - A.c13 * B.c25 + A.c14 * B.c24 + A.c15 * B.c21;
    ctype const c12 = A.c0 * B.c19 + A.c4 * B.c5 + A.c5 * B.c21 - A.c6 * B.c20 -
                      A.c7 * B.c25 + A.c7 * B.c31 - A.c8 * B.c3 + A.c8 * B.c18 +
                      A.c9 * B.c2 - A.c9 * B.c17 - A.c10 * B.c19 - A.c11 * B.c5 -
                      A.c12 * B.c21 + A.c13 * B.c20 - A.c14 * B.c1 + A.c14 * B.c16;
    ctype const c13 = A.c0 * B.c20 - A.c4 * B.c21 + A.c5 * B.c5 + A.c6 * B.c19 +
                      A.c7 * B.c3 - A.c7 * B.c18 - A.c8 * B.c25 + A.c8 * B.c31 -
                      A.c9 * B.c1 + A.c9 * B.c16 - A.c10 * B.c20 + A.c11 * B.c21 -
                      A.c12 * B.c5 - A.c13 * B.c19 - A.c14 * B.c2 + A.c14 * B.c17;
    ctype const c14 = A.c0 * B.c21 + A.c4 * B.c20 - A.c5 * B.c19 + A.c6 * B.c5 -
                      A.c7 * B.c2 + A.c7 * B.c17 + A.c8 * B.c1 - A.c8 * B.c16 -
                      A.c9 * B.c25 + A.c9 * B.c31 - A.c10 * B.c21 - A.c11 * B.c20 +
                      A.c12 * B.c19 - A.c13 * B.c5 - A.c14 * B.c3 + A.c14 * B.c18;
    ctype const c15 = -A.c0 * B.c25 - A.c1 * B.c19 - A.c2 * B.c20 - A.c3 * B.c21 +
                      A.c4 * B.c1 + A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c22 +
                      A.c8 * B.c23 + A.c9 * B.c24 + A.c10 * B.c31 - A.c11 * B.c16 -
                      A.c12 * B.c17 - A.c13 * B.c18 - A.c14 * B.c4 + A.c15 * B.c5;
    ctype const c16 = -A.c0 * B.c9 - A.c1 * B.c29 - A.c2 * B.c14 + A.c3 * B.c13 -
                      A.c4 * B.c0 - A.c5 * B.c11 + A.c6 * B.c10 - A.c7 * B.c30 -
                      A.c8 * B.c8 + A.c9 * B.c7 - A.c10 * B.c26 - A.c11 * B.c15 -
                      A.c12 * B.c28 + A.c13 * B.c27 - A.c14 * B.c6 - A.c15 * B.c12;
    ctype const c17 = -A.c0 * B.c10 + A.c1 * B.c14 - A.c2 * B.c29 - A.c3 * B.c12 +
                      A.c4 * B.c11 - A.c5 * B.c0 - A.c6 * B.c9 + A.c7 * B.c8 -
                      A.c8 * B.c30 - A.c9 * B.c6 - A.c10 * B.c27 + A.c11 * B.c28 -
                      A.c12 * B.c15 - A.c13 * B.c26 - A.c14 * B.c7 - A.c15 * B.c13;
    ctype const c18 = -A.c0 * B.c11 - A.c1 * B.c13 + A.c2 * B.c12 - A.c3 * B.c29 -
                      A.c4 * B.c10 + A.c5 * B.c9 - A.c6 * B.c0 - A.c7 * B.c7 +
                      A.c8 * B.c6 - A.c9 * B.c30 - A.c10 * B.c28 - A.c11 * B.c27 +
                      A.c12 * B.c26 - A.c13 * B.c15 - A.c14 * B.c8 - A.c15 * B.c14;
    ctype const c19 = -A.c0 * B.c12 - A.c4 * B.c29 - A.c5 * B.c14 + A.c6 * B.c13 -
                      A.c7 * B.c0 - A.c7 * B.c15 - A.c8 * B.c11 - A.c8 * B.c28 +
                      A.c9 * B.c10 + A.c9 * B.c27 + A.c10 * B.c12 + A.c11 * B.c29 +
                      A.c12 * B.c14 - A.c13 * B.c13 - A.c14 * B.c9 - A.c14 * B.c26;
    ctype const c20 = -A.c0 * B.c13 + A.c4 * B.c14 - A.c5 * B.c29 - A.c6 * B.c12 +
                      A.c7 * B.c11 + A.c7 * B.c28 - A.c8 * B.c0 - A.c8 * B.c15 -
                      A.c9 * B.c9 - A.c9 * B.c26 + A.c10 * B.c13 - A.c11 * B.c14 +
                      A.c12 * B.c29 + A.c13 * B.c12 - A.c14 * B.c10 - A.c14 * B.c27;
    ctype const c21 = -A.c0 * B.c14 - A.c4 * B.c13 + A.c5 * B.c12 - A.c6 * B.c29 -
                      A.c7 * B.c10 - A.c7 * B.c27 + A.c8 * B.c9 + A.c8 * B.c26 -
                      A.c9 * B.c0 - A.c9 * B.c15 + A.c10 * B.c14 + A.c11 * B.c13 -
                      A.c12 * B.c12 + A.c13 * B.c29 - A.c14 * B.c11 - A.c14 * B.c28;
    ctype const c22 = -A.c0 * B.c6 - A.c1 * B.c0 + A.c1 * B.c15 - A.c2 * B.c11 +
                      A.c2 * B.c28 + A.c3 * B.c10 - A.c3 * B.c27 - A.c4 * B.c30 -
                      A.c5 * B.c8 + A.c6 * B.c7 - A.c10 * B.c6 - A.c11 * B.c30 -
                      A.c12 * B.c8 + A.c13 * B.c7 - A.c15 * B.c9 + A.c15 * B.c26;
    ctype const c23 = -A.c0 * B.c7 + A.c1 * B.c11 - A.c1 * B.c28 - A.c2 * B.c0 +
                      A.c2 * B.c15 - A.c3 * B.c9 + A.c3 * B.c26 + A.c4 * B.c8 -
                      A.c5 * B.c30 - A.c6 * B.c6 - A.c10 * B.c7 + A.c11 * B.c8 -
                      A.c12 * B.c30 - A.c13 * B.c6 - A.c15 * B.c10 + A.c15 * B.c27;
    ctype const c24 = -A.c0 * B.c8 - A.c1 * B.c10 + A.c1 * B.c27 + A.c2 * B.c9 -
                      A.c2 * B.c26 - A.c3 * B.c0 + A.c3 * B.c15 - A.c4 * B.c7 +
                      A.c5 * B.c6 - A.c6 * B.c30 - A.c10 * B.c8 - A.c11 * B.c7 +
                      A.c12 * B.c6 - A.c13 * B.c30 - A.c15 * B.c11 + A.c15 * B.c28;
    ctype const c25 = A.c0 * B.c15 - A.c1 * B.c12 - A.c2 * B.c13 - A.c3 * B.c14 -
                      A.c4 * B.c26 - A.c5 * B.c27 - A.c6 * B.c28 + A.c7 * B.c6 +
                      A.c8 * B.c7 + A.c9 * B.c8 + A.c10 * B.c0 - A.c11 * B.c9 -
                      A.c12 * B.c10 - A.c13 * B.c11 - A.c14 * B.c30 + A.c15 * B.c29;
    ctype const c26 = -A.c0 * B.c1 + A.c1 * B.c5 + A.c2 * B.c21 - A.c3 * B.c20 -
                      A.c4 * B.c25 - A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c4 -
                      A.c8 * B.c24 + A.c9 * B.c23 + A.c10 * B.c16 + A.c11 * B.c31 +
                      A.c12 * B.c18 - A.c13 * B.c17 - A.c14 * B.c22 + A.c15 * B.c19;
    ctype const c27 = -A.c0 * B.c2 - A.c1 * B.c21 + A.c2 * B.c5 + A.c3 * B.c19 +
                      A.c4 * B.c3 - A.c5 * B.c25 - A.c6 * B.c1 + A.c7 * B.c24 -
                      A.c8 * B.c4 - A.c9 * B.c22 + A.c10 * B.c17 - A.c11 * B.c18 +
                      A.c12 * B.c31 + A.c13 * B.c16 - A.c14 * B.c23 + A.c15 * B.c20;
    ctype const c28 = -A.c0 * B.c3 + A.c1 * B.c20 - A.c2 * B.c19 + A.c3 * B.c5 -
                      A.c4 * B.c2 + A.c5 * B.c1 - A.c6 * B.c25 - A.c7 * B.c23 +
                      A.c8 * B.c22 - A.c9 * B.c4 + A.c10 * B.c18 + A.c11 * B.c17 -
                      A.c12 * B.c16 + A.c13 * B.c31 - A.c14 * B.c24 + A.c15 * B.c21;
    ctype const c29 = A.c0 * B.c5 - A.c4 * B.c19 - A.c5 * B.c20 - A.c6 * B.c21 +
                      A.c7 * B.c1 - A.c7 * B.c16 + A.c8 * B.c2 - A.c8 * B.c17 +
                      A.c9 * B.c3 - A.c9 * B.c18 - A.c10 * B.c5 + A.c11 * B.c19 +
                      A.c12 * B.c20 + A.c13 * B.c21 - A.c14 * B.c25 + A.c14 * B.c31;
    ctype const c30 = A.c0 * B.c4 - A.c1 * B.c1 - A.c1 * B.c16 - A.c2 * B.c2 -
                      A.c2 * B.c17 - A.c3 * B.c3 - A.c3 * B.c18 - A.c4 * B.c22 -
                      A.c5 * B.c23 - A.c6 * B.c24 + A.c10 * B.c4 - A.c11 * B.c22 -
                      A.c12 * B.c23 - A.c13 * B.c24 + A.c15 * B.c25 + A.c15 * B.c31;
    ctype const c31 = -A.c0 * B.c0 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                      A.c4 * B.c9 + A.c5 * B.c10 + A.c6 * B.c11 + A.c7 * B.c6 +
                      A.c8 * B.c7 + A.c9 * B.c8 - A.c10 * B.c15 + A.c11 * B.c26 +
                      A.c12 * B.c27 + A.c13 * B.c28 - A.c14 * B.c30 - A.c15 * B.c29;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& A,
                                                 MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c6 * B.c8 - A.c7 * B.c9 - A.c8 * B.c10 -
                     A.c9 * B.c5 - A.c10 * B.c6 - A.c11 * B.c7 - A.c12 * B.c11 -
                     A.c13 * B.c12 - A.c14 * B.c13 - A.c15 * B.c14 + A.c26 * B.c0 +
                     A.c27 * B.c1 + A.c28 * B.c2 + A.c29 * B.c3 + A.c30 * B.c4;
    ctype const c1 = A.c1 * B.c15 + A.c2 * B.c7 - A.c3 * B.c6 - A.c4 * B.c8 +
                     A.c5 * B.c11 + A.c16 * B.c14 + A.c17 * B.c2 - A.c18 * B.c1 +
                     A.c19 * B.c3 + A.c20 * B.c13 - A.c21 * B.c12 - A.c22 * B.c4 -
                     A.c23 * B.c10 + A.c24 * B.c9 + A.c25 * B.c5 + A.c31 * B.c0;
    ctype const c2 = -A.c1 * B.c7 + A.c2 * B.c15 + A.c3 * B.c5 - A.c4 * B.c9 +
                     A.c5 * B.c12 - A.c16 * B.c2 + A.c17 * B.c14 + A.c18 * B.c0 -
                     A.c19 * B.c13 + A.c20 * B.c3 + A.c21 * B.c11 + A.c22 * B.c10 -
                     A.c23 * B.c4 - A.c24 * B.c8 + A.c25 * B.c6 + A.c31 * B.c1;
    ctype const c3 = A.c1 * B.c6 - A.c2 * B.c5 + A.c3 * B.c15 - A.c4 * B.c10 +
                     A.c5 * B.c13 + A.c16 * B.c1 - A.c17 * B.c0 + A.c18 * B.c14 +
                     A.c19 * B.c12 - A.c20 * B.c11 + A.c21 * B.c3 - A.c22 * B.c9 +
                     A.c23 * B.c8 - A.c24 * B.c4 + A.c25 * B.c7 + A.c31 * B.c2;
    ctype const c4 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 +
                     A.c4 * B.c15 - A.c16 * B.c11 - A.c17 * B.c12 - A.c18 * B.c13 -
                     A.c22 * B.c0 - A.c22 * B.c5 - A.c23 * B.c1 - A.c23 * B.c6 -
                     A.c24 * B.c2 - A.c24 * B.c7 - A.c25 * B.c3 + A.c31 * B.c3;
    ctype const c5 = -A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 - A.c5 * B.c14 +
                     A.c5 * B.c15 - A.c16 * B.c8 - A.c17 * B.c9 - A.c18 * B.c10 +
                     A.c19 * B.c0 - A.c19 * B.c5 + A.c20 * B.c1 - A.c20 * B.c6 +
                     A.c21 * B.c2 - A.c21 * B.c7 + A.c25 * B.c4 + A.c31 * B.c4;
    ctype const c6 = A.c0 * B.c11 + A.c6 * B.c14 + A.c6 * B.c15 + A.c7 * B.c2 +
                     A.c7 * B.c7 - A.c8 * B.c1 - A.c8 * B.c6 + A.c9 * B.c3 +
                     A.c10 * B.c13 - A.c11 * B.c12 + A.c15 * B.c11 + A.c26 * B.c3 +
                     A.c27 * B.c13 - A.c28 * B.c12 + A.c30 * B.c0 + A.c30 * B.c5;
    ctype const c7 = A.c0 * B.c12 - A.c6 * B.c2 - A.c6 * B.c7 + A.c7 * B.c14 +
                     A.c7 * B.c15 + A.c8 * B.c0 + A.c8 * B.c5 - A.c9 * B.c13 +
                     A.c10 * B.c3 + A.c11 * B.c11 + A.c15 * B.c12 - A.c26 * B.c13 +
                     A.c27 * B.c3 + A.c28 * B.c11 + A.c30 * B.c1 + A.c30 * B.c6;
    ctype const c8 = A.c0 * B.c13 + A.c6 * B.c1 + A.c6 * B.c6 - A.c7 * B.c0 -
                     A.c7 * B.c5 + A.c8 * B.c14 + A.c8 * B.c15 + A.c9 * B.c12 -
                     A.c10 * B.c11 + A.c11 * B.c3 + A.c15 * B.c13 + A.c26 * B.c12 -
                     A.c27 * B.c11 + A.c28 * B.c3 + A.c30 * B.c2 + A.c30 * B.c7;
    ctype const c9 = A.c0 * B.c5 + A.c6 * B.c4 + A.c7 * B.c10 - A.c8 * B.c9 +
                     A.c9 * B.c15 + A.c10 * B.c7 - A.c11 * B.c6 + A.c12 * B.c3 +
                     A.c13 * B.c13 - A.c14 * B.c12 - A.c15 * B.c0 - A.c26 * B.c14 -
                     A.c27 * B.c2 + A.c28 * B.c1 + A.c29 * B.c11 + A.c30 * B.c8;
    ctype const c10 = A.c0 * B.c6 - A.c6 * B.c10 + A.c7 * B.c4 + A.c8 * B.c8 -
                      A.c9 * B.c7 + A.c10 * B.c15 + A.c11 * B.c5 - A.c12 * B.c13 +
                      A.c13 * B.c3 + A.c14 * B.c11 - A.c15 * B.c1 + A.c26 * B.c2 -
                      A.c27 * B.c14 - A.c28 * B.c0 + A.c29 * B.c12 + A.c30 * B.c9;
    ctype const c11 = A.c0 * B.c7 + A.c6 * B.c9 - A.c7 * B.c8 + A.c8 * B.c4 +
                      A.c9 * B.c6 - A.c10 * B.c5 + A.c11 * B.c15 + A.c12 * B.c12 -
                      A.c13 * B.c11 + A.c14 * B.c3 - A.c15 * B.c2 - A.c26 * B.c1 +
                      A.c27 * B.c0 - A.c28 * B.c14 + A.c29 * B.c13 + A.c30 * B.c10;
    ctype const c12 = A.c0 * B.c8 + A.c9 * B.c4 + A.c10 * B.c10 - A.c11 * B.c9 -
                      A.c12 * B.c14 + A.c12 * B.c15 - A.c13 * B.c2 + A.c13 * B.c7 +
                      A.c14 * B.c1 - A.c14 * B.c6 - A.c15 * B.c8 - A.c26 * B.c4 -
                      A.c27 * B.c10 + A.c28 * B.c9 - A.c29 * B.c0 + A.c29 * B.c5;
    ctype const c13 = A.c0 * B.c9 - A.c9 * B.c10 + A.c10 * B.c4 + A.c11 * B.c8 +
                      A.c12 * B.c2 - A.c12 * B.c7 - A.c13 * B.c14 + A.c13 * B.c15 -
                      A.c14 * B.c0 + A.c14 * B.c5 - A.c15 * B.c9 + A.c26 * B.c10 -
                      A.c27 * B.c4 - A.c28 * B.c8 - A.c29 * B.c1 + A.c29 * B.c6;
    ctype const c14 = A.c0 * B.c10 + A.c9 * B.c9 - A.c10 * B.c8 + A.c11 * B.c4 -
                      A.c12 * B.c1 + A.c12 * B.c6 + A.c13 * B.c0 - A.c13 * B.c5 -
                      A.c14 * B.c14 + A.c14 * B.c15 - A.c15 * B.c10 - A.c26 * B.c9 +
                      A.c27 * B.c8 - A.c28 * B.c4 - A.c29 * B.c2 + A.c29 * B.c7;
    ctype const c15 = -A.c0 * B.c14 - A.c6 * B.c8 - A.c7 * B.c9 - A.c8 * B.c10 +
                      A.c9 * B.c0 + A.c10 * B.c1 + A.c11 * B.c2 + A.c12 * B.c11 +
                      A.c13 * B.c12 + A.c14 * B.c13 + A.c15 * B.c15 - A.c26 * B.c5 -
                      A.c27 * B.c6 - A.c28 * B.c7 - A.c29 * B.c3 + A.c30 * B.c4;
    ctype const c16 = A.c1 * B.c14 + A.c2 * B.c2 - A.c3 * B.c1 + A.c4 * B.c8 +
                      A.c5 * B.c11 + A.c16 * B.c15 + A.c17 * B.c7 - A.c18 * B.c6 +
                      A.c19 * B.c3 + A.c20 * B.c13 - A.c21 * B.c12 + A.c22 * B.c4 +
                      A.c23 * B.c10 - A.c24 * B.c9 + A.c25 * B.c0 + A.c31 * B.c5;
    ctype const c17 = -A.c1 * B.c2 + A.c2 * B.c14 + A.c3 * B.c0 + A.c4 * B.c9 +
                      A.c5 * B.c12 - A.c16 * B.c7 + A.c17 * B.c15 + A.c18 * B.c5 -
                      A.c19 * B.c13 + A.c20 * B.c3 + A.c21 * B.c11 - A.c22 * B.c10 +
                      A.c23 * B.c4 + A.c24 * B.c8 + A.c25 * B.c1 + A.c31 * B.c6;
    ctype const c18 = A.c1 * B.c1 - A.c2 * B.c0 + A.c3 * B.c14 + A.c4 * B.c10 +
                      A.c5 * B.c13 + A.c16 * B.c6 - A.c17 * B.c5 + A.c18 * B.c15 +
                      A.c19 * B.c12 - A.c20 * B.c11 + A.c21 * B.c3 + A.c22 * B.c9 -
                      A.c23 * B.c8 + A.c24 * B.c4 + A.c25 * B.c2 + A.c31 * B.c7;
    ctype const c19 = A.c1 * B.c4 + A.c2 * B.c10 - A.c3 * B.c9 - A.c5 * B.c0 +
                      A.c5 * B.c5 + A.c16 * B.c4 + A.c17 * B.c10 - A.c18 * B.c9 -
                      A.c19 * B.c14 + A.c19 * B.c15 - A.c20 * B.c2 + A.c20 * B.c7 +
                      A.c21 * B.c1 - A.c21 * B.c6 + A.c25 * B.c8 + A.c31 * B.c8;
    ctype const c20 = -A.c1 * B.c10 + A.c2 * B.c4 + A.c3 * B.c8 - A.c5 * B.c1 +
                      A.c5 * B.c6 - A.c16 * B.c10 + A.c17 * B.c4 + A.c18 * B.c8 +
                      A.c19 * B.c2 - A.c19 * B.c7 - A.c20 * B.c14 + A.c20 * B.c15 -
                      A.c21 * B.c0 + A.c21 * B.c5 + A.c25 * B.c9 + A.c31 * B.c9;
    ctype const c21 = A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c4 - A.c5 * B.c2 +
                      A.c5 * B.c7 + A.c16 * B.c9 - A.c17 * B.c8 + A.c18 * B.c4 -
                      A.c19 * B.c1 + A.c19 * B.c6 + A.c20 * B.c0 - A.c20 * B.c5 -
                      A.c21 * B.c14 + A.c21 * B.c15 + A.c25 * B.c10 + A.c31 * B.c10;
    ctype const c22 = -A.c1 * B.c3 - A.c2 * B.c13 + A.c3 * B.c12 + A.c4 * B.c0 +
                      A.c4 * B.c5 + A.c16 * B.c3 + A.c17 * B.c13 - A.c18 * B.c12 +
                      A.c22 * B.c14 + A.c22 * B.c15 + A.c23 * B.c2 + A.c23 * B.c7 -
                      A.c24 * B.c1 - A.c24 * B.c6 - A.c25 * B.c11 + A.c31 * B.c11;
    ctype const c23 = A.c1 * B.c13 - A.c2 * B.c3 - A.c3 * B.c11 + A.c4 * B.c1 +
                      A.c4 * B.c6 - A.c16 * B.c13 + A.c17 * B.c3 + A.c18 * B.c11 -
                      A.c22 * B.c2 - A.c22 * B.c7 + A.c23 * B.c14 + A.c23 * B.c15 +
                      A.c24 * B.c0 + A.c24 * B.c5 - A.c25 * B.c12 + A.c31 * B.c12;
    ctype const c24 = -A.c1 * B.c12 + A.c2 * B.c11 - A.c3 * B.c3 + A.c4 * B.c2 +
                      A.c4 * B.c7 + A.c16 * B.c12 - A.c17 * B.c11 + A.c18 * B.c3 +
                      A.c22 * B.c1 + A.c22 * B.c6 - A.c23 * B.c0 - A.c23 * B.c5 +
                      A.c24 * B.c14 + A.c24 * B.c15 - A.c25 * B.c13 + A.c31 * B.c13;
    ctype const c25 = -A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 - A.c4 * B.c4 +
                      A.c5 * B.c3 - A.c16 * B.c0 - A.c17 * B.c1 - A.c18 * B.c2 -
                      A.c19 * B.c11 - A.c20 * B.c12 - A.c21 * B.c13 + A.c22 * B.c8 +
                      A.c23 * B.c9 + A.c24 * B.c10 + A.c25 * B.c15 + A.c31 * B.c14;
    ctype const c26 = -A.c0 * B.c0 + A.c6 * B.c4 + A.c7 * B.c10 - A.c8 * B.c9 -
                      A.c9 * B.c14 - A.c10 * B.c2 + A.c11 * B.c1 - A.c12 * B.c3 -
                      A.c13 * B.c13 + A.c14 * B.c12 + A.c15 * B.c5 + A.c26 * B.c15 +
                      A.c27 * B.c7 - A.c28 * B.c6 - A.c29 * B.c11 + A.c30 * B.c8;
    ctype const c27 = -A.c0 * B.c1 - A.c6 * B.c10 + A.c7 * B.c4 + A.c8 * B.c8 +
                      A.c9 * B.c2 - A.c10 * B.c14 - A.c11 * B.c0 + A.c12 * B.c13 -
                      A.c13 * B.c3 - A.c14 * B.c11 + A.c15 * B.c6 - A.c26 * B.c7 +
                      A.c27 * B.c15 + A.c28 * B.c5 - A.c29 * B.c12 + A.c30 * B.c9;
    ctype const c28 = -A.c0 * B.c2 + A.c6 * B.c9 - A.c7 * B.c8 + A.c8 * B.c4 -
                      A.c9 * B.c1 + A.c10 * B.c0 - A.c11 * B.c14 - A.c12 * B.c12 +
                      A.c13 * B.c11 - A.c14 * B.c3 + A.c15 * B.c7 + A.c26 * B.c6 -
                      A.c27 * B.c5 + A.c28 * B.c15 - A.c29 * B.c13 + A.c30 * B.c10;
    ctype const c29 = A.c0 * B.c4 - A.c9 * B.c8 - A.c10 * B.c9 - A.c11 * B.c10 +
                      A.c12 * B.c0 - A.c12 * B.c5 + A.c13 * B.c1 - A.c13 * B.c6 +
                      A.c14 * B.c2 - A.c14 * B.c7 - A.c15 * B.c4 + A.c26 * B.c8 +
                      A.c27 * B.c9 + A.c28 * B.c10 - A.c29 * B.c14 + A.c29 * B.c15;
    ctype const c30 = A.c0 * B.c3 - A.c6 * B.c0 - A.c6 * B.c5 - A.c7 * B.c1 -
                      A.c7 * B.c6 - A.c8 * B.c2 - A.c8 * B.c7 - A.c9 * B.c11 -
                      A.c10 * B.c12 - A.c11 * B.c13 + A.c15 * B.c3 - A.c26 * B.c11 -
                      A.c27 * B.c12 - A.c28 * B.c13 + A.c30 * B.c14 + A.c30 * B.c15;
    ctype const c31 = -A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2 + A.c4 * B.c4 +
                      A.c5 * B.c3 - A.c16 * B.c5 - A.c17 * B.c6 - A.c18 * B.c7 -
                      A.c19 * B.c11 - A.c20 * B.c12 - A.c21 * B.c13 - A.c22 * B.c8 -
                      A.c23 * B.c9 - A.c24 * B.c10 + A.c25 * B.c14 + A.c31 * B.c15;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& A,
                                                 MVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c26 + A.c1 * B.c27 + A.c2 * B.c28 + A.c3 * B.c29 +
                     A.c4 * B.c30 - A.c5 * B.c9 - A.c6 * B.c10 - A.c7 * B.c11 -
                     A.c8 * B.c6 - A.c9 * B.c7 - A.c10 * B.c8 - A.c11 * B.c12 -
                     A.c12 * B.c13 - A.c13 * B.c14 - A.c14 * B.c15 + A.c15 * B.c0;
    ctype const c1 = A.c0 * B.c31 + A.c1 * B.c18 - A.c2 * B.c17 - A.c3 * B.c19 +
                     A.c4 * B.c22 + A.c5 * B.c25 + A.c6 * B.c3 - A.c7 * B.c2 +
                     A.c8 * B.c4 + A.c9 * B.c24 - A.c10 * B.c23 - A.c11 * B.c5 -
                     A.c12 * B.c21 + A.c13 * B.c20 + A.c14 * B.c16 + A.c15 * B.c1;
    ctype const c2 = -A.c0 * B.c18 + A.c1 * B.c31 + A.c2 * B.c16 - A.c3 * B.c20 +
                     A.c4 * B.c23 - A.c5 * B.c3 + A.c6 * B.c25 + A.c7 * B.c1 -
                     A.c8 * B.c24 + A.c9 * B.c4 + A.c10 * B.c22 + A.c11 * B.c21 -
                     A.c12 * B.c5 - A.c13 * B.c19 + A.c14 * B.c17 + A.c15 * B.c2;
    ctype const c3 = A.c0 * B.c17 - A.c1 * B.c16 + A.c2 * B.c31 - A.c3 * B.c21 +
                     A.c4 * B.c24 + A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c25 +
                     A.c8 * B.c23 - A.c9 * B.c22 + A.c10 * B.c4 - A.c11 * B.c20 +
                     A.c12 * B.c19 - A.c13 * B.c5 + A.c14 * B.c18 + A.c15 * B.c3;
    ctype const c4 = A.c0 * B.c22 + A.c1 * B.c23 + A.c2 * B.c24 + A.c3 * B.c25 +
                     A.c3 * B.c31 - A.c5 * B.c22 - A.c6 * B.c23 - A.c7 * B.c24 -
                     A.c11 * B.c1 - A.c11 * B.c16 - A.c12 * B.c2 - A.c12 * B.c17 -
                     A.c13 * B.c3 - A.c13 * B.c18 - A.c14 * B.c4 + A.c15 * B.c4;
    ctype const c5 = -A.c0 * B.c19 - A.c1 * B.c20 - A.c2 * B.c21 - A.c4 * B.c25 +
                     A.c4 * B.c31 - A.c5 * B.c19 - A.c6 * B.c20 - A.c7 * B.c21 +
                     A.c8 * B.c1 - A.c8 * B.c16 + A.c9 * B.c2 - A.c9 * B.c17 +
                     A.c10 * B.c3 - A.c10 * B.c18 + A.c14 * B.c5 + A.c15 * B.c5;
    ctype const c6 = -A.c0 * B.c30 - A.c1 * B.c8 + A.c2 * B.c7 + A.c3 * B.c9 -
                     A.c3 * B.c26 + A.c5 * B.c30 + A.c6 * B.c8 - A.c7 * B.c7 +
                     A.c11 * B.c0 - A.c11 * B.c15 + A.c12 * B.c11 - A.c12 * B.c28 -
                     A.c13 * B.c10 + A.c13 * B.c27 - A.c14 * B.c6 + A.c15 * B.c6;
    ctype const c7 = A.c0 * B.c8 - A.c1 * B.c30 - A.c2 * B.c6 + A.c3 * B.c10 -
                     A.c3 * B.c27 - A.c5 * B.c8 + A.c6 * B.c30 + A.c7 * B.c6 -
                     A.c11 * B.c11 + A.c11 * B.c28 + A.c12 * B.c0 - A.c12 * B.c15 +
                     A.c13 * B.c9 - A.c13 * B.c26 - A.c14 * B.c7 + A.c15 * B.c7;
    ctype const c8 = -A.c0 * B.c7 + A.c1 * B.c6 - A.c2 * B.c30 + A.c3 * B.c11 -
                     A.c3 * B.c28 + A.c5 * B.c7 - A.c6 * B.c6 + A.c7 * B.c30 +
                     A.c11 * B.c10 - A.c11 * B.c27 - A.c12 * B.c9 + A.c12 * B.c26 +
                     A.c13 * B.c0 - A.c13 * B.c15 - A.c14 * B.c8 + A.c15 * B.c8;
    ctype const c9 = -A.c0 * B.c15 - A.c1 * B.c28 + A.c2 * B.c27 + A.c3 * B.c12 +
                     A.c4 * B.c6 + A.c5 * B.c0 + A.c6 * B.c11 - A.c7 * B.c10 +
                     A.c8 * B.c30 + A.c9 * B.c8 - A.c10 * B.c7 + A.c11 * B.c29 +
                     A.c12 * B.c14 - A.c13 * B.c13 - A.c14 * B.c26 + A.c15 * B.c9;
    ctype const c10 = A.c0 * B.c28 - A.c1 * B.c15 - A.c2 * B.c26 + A.c3 * B.c13 +
                      A.c4 * B.c7 - A.c5 * B.c11 + A.c6 * B.c0 + A.c7 * B.c9 -
                      A.c8 * B.c8 + A.c9 * B.c30 + A.c10 * B.c6 - A.c11 * B.c14 +
                      A.c12 * B.c29 + A.c13 * B.c12 - A.c14 * B.c27 + A.c15 * B.c10;
    ctype const c11 = -A.c0 * B.c27 + A.c1 * B.c26 - A.c2 * B.c15 + A.c3 * B.c14 +
                      A.c4 * B.c8 + A.c5 * B.c10 - A.c6 * B.c9 + A.c7 * B.c0 +
                      A.c8 * B.c7 - A.c9 * B.c6 + A.c10 * B.c30 + A.c11 * B.c13 -
                      A.c12 * B.c12 + A.c13 * B.c29 - A.c14 * B.c28 + A.c15 * B.c11;
    ctype const c12 = A.c0 * B.c29 + A.c1 * B.c14 - A.c2 * B.c13 + A.c4 * B.c9 +
                      A.c4 * B.c26 + A.c5 * B.c29 + A.c6 * B.c14 - A.c7 * B.c13 +
                      A.c8 * B.c0 + A.c8 * B.c15 + A.c9 * B.c11 + A.c9 * B.c28 -
                      A.c10 * B.c10 - A.c10 * B.c27 + A.c14 * B.c12 + A.c15 * B.c12;
    ctype const c13 = -A.c0 * B.c14 + A.c1 * B.c29 + A.c2 * B.c12 + A.c4 * B.c10 +
                      A.c4 * B.c27 - A.c5 * B.c14 + A.c6 * B.c29 + A.c7 * B.c12 -
                      A.c8 * B.c11 - A.c8 * B.c28 + A.c9 * B.c0 + A.c9 * B.c15 +
                      A.c10 * B.c9 + A.c10 * B.c26 + A.c14 * B.c13 + A.c15 * B.c13;
    ctype const c14 = A.c0 * B.c13 - A.c1 * B.c12 + A.c2 * B.c29 + A.c4 * B.c11 +
                      A.c4 * B.c28 + A.c5 * B.c13 - A.c6 * B.c12 + A.c7 * B.c29 +
                      A.c8 * B.c10 + A.c8 * B.c27 - A.c9 * B.c9 - A.c9 * B.c26 +
                      A.c10 * B.c0 + A.c10 * B.c15 + A.c14 * B.c14 + A.c15 * B.c14;
    ctype const c15 = A.c0 * B.c9 + A.c1 * B.c10 + A.c2 * B.c11 + A.c3 * B.c29 -
                      A.c4 * B.c30 - A.c5 * B.c26 - A.c6 * B.c27 - A.c7 * B.c28 +
                      A.c8 * B.c6 + A.c9 * B.c7 + A.c10 * B.c8 - A.c11 * B.c12 -
                      A.c12 * B.c13 - A.c13 * B.c14 - A.c14 * B.c0 + A.c15 * B.c15;
    ctype const c16 = A.c0 * B.c25 + A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c19 +
                      A.c4 * B.c22 + A.c5 * B.c31 + A.c6 * B.c18 - A.c7 * B.c17 +
                      A.c8 * B.c4 + A.c9 * B.c24 - A.c10 * B.c23 + A.c11 * B.c5 +
                      A.c12 * B.c21 - A.c13 * B.c20 + A.c14 * B.c1 + A.c15 * B.c16;
    ctype const c17 = -A.c0 * B.c3 + A.c1 * B.c25 + A.c2 * B.c1 + A.c3 * B.c20 +
                      A.c4 * B.c23 - A.c5 * B.c18 + A.c6 * B.c31 + A.c7 * B.c16 -
                      A.c8 * B.c24 + A.c9 * B.c4 + A.c10 * B.c22 - A.c11 * B.c21 +
                      A.c12 * B.c5 + A.c13 * B.c19 + A.c14 * B.c2 + A.c15 * B.c17;
    ctype const c18 = A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c25 + A.c3 * B.c21 +
                      A.c4 * B.c24 + A.c5 * B.c17 - A.c6 * B.c16 + A.c7 * B.c31 +
                      A.c8 * B.c23 - A.c9 * B.c22 + A.c10 * B.c4 + A.c11 * B.c20 -
                      A.c12 * B.c19 + A.c13 * B.c5 + A.c14 * B.c3 + A.c15 * B.c18;
    ctype const c19 = A.c0 * B.c5 + A.c1 * B.c21 - A.c2 * B.c20 - A.c4 * B.c1 +
                      A.c4 * B.c16 + A.c5 * B.c5 + A.c6 * B.c21 - A.c7 * B.c20 -
                      A.c8 * B.c25 + A.c8 * B.c31 - A.c9 * B.c3 + A.c9 * B.c18 +
                      A.c10 * B.c2 - A.c10 * B.c17 + A.c14 * B.c19 + A.c15 * B.c19;
    ctype const c20 = -A.c0 * B.c21 + A.c1 * B.c5 + A.c2 * B.c19 - A.c4 * B.c2 +
                      A.c4 * B.c17 - A.c5 * B.c21 + A.c6 * B.c5 + A.c7 * B.c19 +
                      A.c8 * B.c3 - A.c8 * B.c18 - A.c9 * B.c25 + A.c9 * B.c31 -
                      A.c10 * B.c1 + A.c10 * B.c16 + A.c14 * B.c20 + A.c15 * B.c20;
    ctype const c21 = A.c0 * B.c20 - A.c1 * B.c19 + A.c2 * B.c5 - A.c4 * B.c3 +
                      A.c4 * B.c18 + A.c5 * B.c20 - A.c6 * B.c19 + A.c7 * B.c5 -
                      A.c8 * B.c2 + A.c8 * B.c17 + A.c9 * B.c1 - A.c9 * B.c16 -
                      A.c10 * B.c25 + A.c10 * B.c31 + A.c14 * B.c21 + A.c15 * B.c21;
    ctype const c22 = -A.c0 * B.c4 - A.c1 * B.c24 + A.c2 * B.c23 + A.c3 * B.c1 +
                      A.c3 * B.c16 + A.c5 * B.c4 + A.c6 * B.c24 - A.c7 * B.c23 +
                      A.c11 * B.c25 + A.c11 * B.c31 + A.c12 * B.c3 + A.c12 * B.c18 -
                      A.c13 * B.c2 - A.c13 * B.c17 - A.c14 * B.c22 + A.c15 * B.c22;
    ctype const c23 = A.c0 * B.c24 - A.c1 * B.c4 - A.c2 * B.c22 + A.c3 * B.c2 +
                      A.c3 * B.c17 - A.c5 * B.c24 + A.c6 * B.c4 + A.c7 * B.c22 -
                      A.c11 * B.c3 - A.c11 * B.c18 + A.c12 * B.c25 + A.c12 * B.c31 +
                      A.c13 * B.c1 + A.c13 * B.c16 - A.c14 * B.c23 + A.c15 * B.c23;
    ctype const c24 = -A.c0 * B.c23 + A.c1 * B.c22 - A.c2 * B.c4 + A.c3 * B.c3 +
                      A.c3 * B.c18 + A.c5 * B.c23 - A.c6 * B.c22 + A.c7 * B.c4 +
                      A.c11 * B.c2 + A.c11 * B.c17 - A.c12 * B.c1 - A.c12 * B.c16 +
                      A.c13 * B.c25 + A.c13 * B.c31 - A.c14 * B.c24 + A.c15 * B.c24;
    ctype const c25 = -A.c0 * B.c16 - A.c1 * B.c17 - A.c2 * B.c18 - A.c3 * B.c5 +
                      A.c4 * B.c4 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 -
                      A.c8 * B.c22 - A.c9 * B.c23 - A.c10 * B.c24 + A.c11 * B.c19 +
                      A.c12 * B.c20 + A.c13 * B.c21 + A.c14 * B.c31 + A.c15 * B.c25;
    ctype const c26 = -A.c0 * B.c0 - A.c1 * B.c11 + A.c2 * B.c10 + A.c3 * B.c12 -
                      A.c4 * B.c6 + A.c5 * B.c15 + A.c6 * B.c28 - A.c7 * B.c27 -
                      A.c8 * B.c30 - A.c9 * B.c8 + A.c10 * B.c7 + A.c11 * B.c29 +
                      A.c12 * B.c14 - A.c13 * B.c13 - A.c14 * B.c9 + A.c15 * B.c26;
    ctype const c27 = A.c0 * B.c11 - A.c1 * B.c0 - A.c2 * B.c9 + A.c3 * B.c13 -
                      A.c4 * B.c7 - A.c5 * B.c28 + A.c6 * B.c15 + A.c7 * B.c26 +
                      A.c8 * B.c8 - A.c9 * B.c30 - A.c10 * B.c6 - A.c11 * B.c14 +
                      A.c12 * B.c29 + A.c13 * B.c12 - A.c14 * B.c10 + A.c15 * B.c27;
    ctype const c28 = -A.c0 * B.c10 + A.c1 * B.c9 - A.c2 * B.c0 + A.c3 * B.c14 -
                      A.c4 * B.c8 + A.c5 * B.c27 - A.c6 * B.c26 + A.c7 * B.c15 -
                      A.c8 * B.c7 + A.c9 * B.c6 - A.c10 * B.c30 + A.c11 * B.c13 -
                      A.c12 * B.c12 + A.c13 * B.c29 - A.c14 * B.c11 + A.c15 * B.c28;
    ctype const c29 = -A.c0 * B.c12 - A.c1 * B.c13 - A.c2 * B.c14 + A.c4 * B.c0 +
                      A.c4 * B.c15 - A.c5 * B.c12 - A.c6 * B.c13 - A.c7 * B.c14 -
                      A.c8 * B.c9 - A.c8 * B.c26 - A.c9 * B.c10 - A.c9 * B.c27 -
                      A.c10 * B.c11 - A.c10 * B.c28 + A.c14 * B.c29 + A.c15 * B.c29;
    ctype const c30 = A.c0 * B.c6 + A.c1 * B.c7 + A.c2 * B.c8 + A.c3 * B.c0 -
                      A.c3 * B.c15 - A.c5 * B.c6 - A.c6 * B.c7 - A.c7 * B.c8 -
                      A.c11 * B.c9 + A.c11 * B.c26 - A.c12 * B.c10 + A.c12 * B.c27 -
                      A.c13 * B.c11 + A.c13 * B.c28 - A.c14 * B.c30 + A.c15 * B.c30;
    ctype const c31 = -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c5 +
                      A.c4 * B.c4 - A.c5 * B.c16 - A.c6 * B.c17 - A.c7 * B.c18 -
                      A.c8 * B.c22 - A.c9 * B.c23 - A.c10 * B.c24 - A.c11 * B.c19 -
                      A.c12 * B.c20 - A.c13 * B.c21 + A.c14 * B.c25 + A.c15 * B.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& M, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    ctype const c4 = M.c4 * ctype(ps);
    ctype const c5 = M.c5 * ctype(ps);
    ctype const c6 = M.c6 * ctype(ps);
    ctype const c7 = M.c7 * ctype(ps);
    ctype const c8 = M.c8 * ctype(ps);
    ctype const c9 = M.c9 * ctype(ps);
    ctype const c10 = M.c10 * ctype(ps);
    ctype const c11 = M.c11 * ctype(ps);
    ctype const c12 = M.c12 * ctype(ps);
    ctype const c13 = M.c13 * ctype(ps);
    ctype const c14 = M.c14 * ctype(ps);
    ctype const c15 = M.c15 * ctype(ps);
    ctype const c16 = M.c16 * ctype(ps);
    ctype const c17 = M.c17 * ctype(ps);
    ctype const c18 = M.c18 * ctype(ps);
    ctype const c19 = M.c19 * ctype(ps);
    ctype const c20 = M.c20 * ctype(ps);
    ctype const c21 = M.c21 * ctype(ps);
    ctype const c22 = M.c22 * ctype(ps);
    ctype const c23 = M.c23 * ctype(ps);
    ctype const c24 = M.c24 * ctype(ps);
    ctype const c25 = M.c25 * ctype(ps);
    ctype const c26 = M.c26 * ctype(ps);
    ctype const c27 = M.c27 * ctype(ps);
    ctype const c28 = M.c28 * ctype(ps);
    ctype const c29 = M.c29 * ctype(ps);
    ctype const c30 = M.c30 * ctype(ps);
    ctype const c31 = M.c31 * ctype(ps);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    ctype const c4 = ctype(ps) * M.c4;
    ctype const c5 = ctype(ps) * M.c5;
    ctype const c6 = ctype(ps) * M.c6;
    ctype const c7 = ctype(ps) * M.c7;
    ctype const c8 = ctype(ps) * M.c8;
    ctype const c9 = ctype(ps) * M.c9;
    ctype const c10 = ctype(ps) * M.c10;
    ctype const c11 = ctype(ps) * M.c11;
    ctype const c12 = ctype(ps) * M.c12;
    ctype const c13 = ctype(ps) * M.c13;
    ctype const c14 = ctype(ps) * M.c14;
    ctype const c15 = ctype(ps) * M.c15;
    ctype const c16 = ctype(ps) * M.c16;
    ctype const c17 = ctype(ps) * M.c17;
    ctype const c18 = ctype(ps) * M.c18;
    ctype const c19 = ctype(ps) * M.c19;
    ctype const c20 = ctype(ps) * M.c20;
    ctype const c21 = ctype(ps) * M.c21;
    ctype const c22 = ctype(ps) * M.c22;
    ctype const c23 = ctype(ps) * M.c23;
    ctype const c24 = ctype(ps) * M.c24;
    ctype const c25 = ctype(ps) * M.c25;
    ctype const c26 = ctype(ps) * M.c26;
    ctype const c27 = ctype(ps) * M.c27;
    ctype const c28 = ctype(ps) * M.c28;
    ctype const c29 = ctype(ps) * M.c29;
    ctype const c30 = ctype(ps) * M.c30;
    ctype const c31 = ctype(ps) * M.c31;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,quadvec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& M,
                                                 QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * Q.x + M.c2 * Q.y + M.c3 * Q.z + M.c4 * Q.w + M.c5 * Q.u;
    ctype const c1 = M.c0 * Q.x + M.c6 * Q.w + M.c10 * Q.z - M.c11 * Q.y - M.c12 * Q.u;
    ctype const c2 = M.c0 * Q.y + M.c7 * Q.w - M.c9 * Q.z + M.c11 * Q.x - M.c13 * Q.u;
    ctype const c3 = M.c0 * Q.z + M.c8 * Q.w + M.c9 * Q.y - M.c10 * Q.x - M.c14 * Q.u;
    ctype const c4 = -M.c0 * Q.u - M.c6 * Q.x - M.c7 * Q.y - M.c8 * Q.z - M.c15 * Q.u;
    ctype const c5 = -M.c0 * Q.w + M.c12 * Q.x + M.c13 * Q.y + M.c14 * Q.z + M.c15 * Q.w;
    ctype const c6 = -M.c1 * Q.u - M.c4 * Q.x + M.c16 * Q.u - M.c23 * Q.z + M.c24 * Q.y;
    ctype const c7 = -M.c2 * Q.u - M.c4 * Q.y + M.c17 * Q.u + M.c22 * Q.z - M.c24 * Q.x;
    ctype const c8 = -M.c3 * Q.u - M.c4 * Q.z + M.c18 * Q.u - M.c22 * Q.y + M.c23 * Q.x;
    ctype const c9 = -M.c2 * Q.z + M.c3 * Q.y + M.c19 * Q.u + M.c22 * Q.w - M.c25 * Q.x;
    ctype const c10 = M.c1 * Q.z - M.c3 * Q.x + M.c20 * Q.u + M.c23 * Q.w - M.c25 * Q.y;
    ctype const c11 = -M.c1 * Q.y + M.c2 * Q.x + M.c21 * Q.u + M.c24 * Q.w - M.c25 * Q.z;
    ctype const c12 = M.c1 * Q.w + M.c5 * Q.x + M.c16 * Q.w + M.c20 * Q.z - M.c21 * Q.y;
    ctype const c13 = M.c2 * Q.w + M.c5 * Q.y + M.c17 * Q.w - M.c19 * Q.z + M.c21 * Q.x;
    ctype const c14 = M.c3 * Q.w + M.c5 * Q.z + M.c18 * Q.w + M.c19 * Q.y - M.c20 * Q.x;
    ctype const c15 = M.c4 * Q.w - M.c5 * Q.u - M.c16 * Q.x - M.c17 * Q.y - M.c18 * Q.z;
    ctype const c16 = -M.c6 * Q.w - M.c12 * Q.u - M.c15 * Q.x - M.c27 * Q.z + M.c28 * Q.y;
    ctype const c17 = -M.c7 * Q.w - M.c13 * Q.u - M.c15 * Q.y + M.c26 * Q.z - M.c28 * Q.x;
    ctype const c18 = -M.c8 * Q.w - M.c14 * Q.u - M.c15 * Q.z - M.c26 * Q.y + M.c27 * Q.x;
    ctype const c19 = -M.c9 * Q.w - M.c13 * Q.z + M.c14 * Q.y + M.c26 * Q.w - M.c29 * Q.x;
    ctype const c20 =
        -M.c10 * Q.w + M.c12 * Q.z - M.c14 * Q.x + M.c27 * Q.w - M.c29 * Q.y;
    ctype const c21 =
        -M.c11 * Q.w - M.c12 * Q.y + M.c13 * Q.x + M.c28 * Q.w - M.c29 * Q.z;
    ctype const c22 = M.c7 * Q.z - M.c8 * Q.y - M.c9 * Q.u - M.c26 * Q.u + M.c30 * Q.x;
    ctype const c23 = -M.c6 * Q.z + M.c8 * Q.x - M.c10 * Q.u - M.c27 * Q.u + M.c30 * Q.y;
    ctype const c24 = M.c6 * Q.y - M.c7 * Q.x - M.c11 * Q.u - M.c28 * Q.u + M.c30 * Q.z;
    ctype const c25 = -M.c9 * Q.x - M.c10 * Q.y - M.c11 * Q.z - M.c29 * Q.u + M.c30 * Q.w;
    ctype const c26 = M.c17 * Q.z - M.c18 * Q.y - M.c19 * Q.u + M.c22 * Q.w + M.c31 * Q.x;
    ctype const c27 =
        -M.c16 * Q.z + M.c18 * Q.x - M.c20 * Q.u + M.c23 * Q.w + M.c31 * Q.y;
    ctype const c28 = M.c16 * Q.y - M.c17 * Q.x - M.c21 * Q.u + M.c24 * Q.w + M.c31 * Q.z;
    ctype const c29 =
        -M.c19 * Q.x - M.c20 * Q.y - M.c21 * Q.z + M.c25 * Q.w + M.c31 * Q.w;
    ctype const c30 = M.c22 * Q.x + M.c23 * Q.y + M.c24 * Q.z - M.c25 * Q.u + M.c31 * Q.u;
    ctype const c31 = M.c26 * Q.x + M.c27 * Q.y + M.c28 * Q.z - M.c29 * Q.u - M.c30 * Q.w;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(quadvec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q,
                                                 MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * M.c1 + Q.y * M.c2 + Q.z * M.c3 + Q.w * M.c4 + Q.u * M.c5;
    ctype const c1 = Q.x * M.c0 + Q.y * M.c11 - Q.z * M.c10 - Q.w * M.c6 + Q.u * M.c12;
    ctype const c2 = -Q.x * M.c11 + Q.y * M.c0 + Q.z * M.c9 - Q.w * M.c7 + Q.u * M.c13;
    ctype const c3 = Q.x * M.c10 - Q.y * M.c9 + Q.z * M.c0 - Q.w * M.c8 + Q.u * M.c14;
    ctype const c4 = Q.x * M.c6 + Q.y * M.c7 + Q.z * M.c8 - Q.u * M.c0 + Q.u * M.c15;
    ctype const c5 = -Q.x * M.c12 - Q.y * M.c13 - Q.z * M.c14 - Q.w * M.c0 - Q.w * M.c15;
    ctype const c6 = Q.x * M.c4 + Q.y * M.c24 - Q.z * M.c23 + Q.u * M.c1 + Q.u * M.c16;
    ctype const c7 = -Q.x * M.c24 + Q.y * M.c4 + Q.z * M.c22 + Q.u * M.c2 + Q.u * M.c17;
    ctype const c8 = Q.x * M.c23 - Q.y * M.c22 + Q.z * M.c4 + Q.u * M.c3 + Q.u * M.c18;
    ctype const c9 = -Q.x * M.c25 - Q.y * M.c3 + Q.z * M.c2 + Q.w * M.c22 + Q.u * M.c19;
    ctype const c10 = Q.x * M.c3 - Q.y * M.c25 - Q.z * M.c1 + Q.w * M.c23 + Q.u * M.c20;
    ctype const c11 = -Q.x * M.c2 + Q.y * M.c1 - Q.z * M.c25 + Q.w * M.c24 + Q.u * M.c21;
    ctype const c12 = -Q.x * M.c5 - Q.y * M.c21 + Q.z * M.c20 - Q.w * M.c1 + Q.w * M.c16;
    ctype const c13 = Q.x * M.c21 - Q.y * M.c5 - Q.z * M.c19 - Q.w * M.c2 + Q.w * M.c17;
    ctype const c14 = -Q.x * M.c20 + Q.y * M.c19 - Q.z * M.c5 - Q.w * M.c3 + Q.w * M.c18;
    ctype const c15 = -Q.x * M.c16 - Q.y * M.c17 - Q.z * M.c18 - Q.w * M.c4 + Q.u * M.c5;
    ctype const c16 = -Q.x * M.c15 - Q.y * M.c28 + Q.z * M.c27 - Q.w * M.c6 - Q.u * M.c12;
    ctype const c17 = Q.x * M.c28 - Q.y * M.c15 - Q.z * M.c26 - Q.w * M.c7 - Q.u * M.c13;
    ctype const c18 = -Q.x * M.c27 + Q.y * M.c26 - Q.z * M.c15 - Q.w * M.c8 - Q.u * M.c14;
    ctype const c19 = Q.x * M.c29 + Q.y * M.c14 - Q.z * M.c13 - Q.w * M.c9 - Q.w * M.c26;
    ctype const c20 =
        -Q.x * M.c14 + Q.y * M.c29 + Q.z * M.c12 - Q.w * M.c10 - Q.w * M.c27;
    ctype const c21 = Q.x * M.c13 - Q.y * M.c12 + Q.z * M.c29 - Q.w * M.c11 - Q.w * M.c28;
    ctype const c22 = -Q.x * M.c30 - Q.y * M.c8 + Q.z * M.c7 - Q.u * M.c9 + Q.u * M.c26;
    ctype const c23 = Q.x * M.c8 - Q.y * M.c30 - Q.z * M.c6 - Q.u * M.c10 + Q.u * M.c27;
    ctype const c24 = -Q.x * M.c7 + Q.y * M.c6 - Q.z * M.c30 - Q.u * M.c11 + Q.u * M.c28;
    ctype const c25 = -Q.x * M.c9 - Q.y * M.c10 - Q.z * M.c11 - Q.w * M.c30 + Q.u * M.c29;
    ctype const c26 = Q.x * M.c31 + Q.y * M.c18 - Q.z * M.c17 - Q.w * M.c22 + Q.u * M.c19;
    ctype const c27 =
        -Q.x * M.c18 + Q.y * M.c31 + Q.z * M.c16 - Q.w * M.c23 + Q.u * M.c20;
    ctype const c28 = Q.x * M.c17 - Q.y * M.c16 + Q.z * M.c31 - Q.w * M.c24 + Q.u * M.c21;
    ctype const c29 = Q.x * M.c19 + Q.y * M.c20 + Q.z * M.c21 - Q.w * M.c25 + Q.w * M.c31;
    ctype const c30 =
        -Q.x * M.c22 - Q.y * M.c23 - Q.z * M.c24 + Q.u * M.c25 + Q.u * M.c31;
    ctype const c31 = Q.x * M.c26 + Q.y * M.c27 + Q.z * M.c28 - Q.w * M.c30 - Q.u * M.c29;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& M,
                                                 TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c6 * t.mx - M.c7 * t.my - M.c8 * t.mz - M.c9 * t.vx -
                     M.c10 * t.vy - M.c11 * t.vz - M.c12 * t.px - M.c13 * t.py -
                     M.c14 * t.pz - M.c15 * t.pw;
    ctype const c1 = M.c2 * t.vz - M.c3 * t.vy - M.c4 * t.mx + M.c5 * t.px +
                     M.c16 * t.pw + M.c20 * t.pz - M.c21 * t.py - M.c23 * t.mz +
                     M.c24 * t.my + M.c25 * t.vx;
    ctype const c2 = -M.c1 * t.vz + M.c3 * t.vx - M.c4 * t.my + M.c5 * t.py +
                     M.c17 * t.pw - M.c19 * t.pz + M.c21 * t.px + M.c22 * t.mz -
                     M.c24 * t.mx + M.c25 * t.vy;
    ctype const c3 = M.c1 * t.vy - M.c2 * t.vx - M.c4 * t.mz + M.c5 * t.pz +
                     M.c18 * t.pw + M.c19 * t.py - M.c20 * t.px - M.c22 * t.my +
                     M.c23 * t.mx + M.c25 * t.vz;
    ctype const c4 = M.c1 * t.px + M.c2 * t.py + M.c3 * t.pz + M.c4 * t.pw -
                     M.c16 * t.px - M.c17 * t.py - M.c18 * t.pz - M.c22 * t.vx -
                     M.c23 * t.vy - M.c24 * t.vz;
    ctype const c5 = -M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz - M.c5 * t.pw -
                     M.c16 * t.mx - M.c17 * t.my - M.c18 * t.mz - M.c19 * t.vx -
                     M.c20 * t.vy - M.c21 * t.vz;
    ctype const c6 = M.c0 * t.px + M.c6 * t.pw + M.c7 * t.vz - M.c8 * t.vy +
                     M.c10 * t.pz - M.c11 * t.py + M.c15 * t.px + M.c27 * t.pz -
                     M.c28 * t.py + M.c30 * t.vx;
    ctype const c7 = M.c0 * t.py - M.c6 * t.vz + M.c7 * t.pw + M.c8 * t.vx - M.c9 * t.pz +
                     M.c11 * t.px + M.c15 * t.py - M.c26 * t.pz + M.c28 * t.px +
                     M.c30 * t.vy;
    ctype const c8 = M.c0 * t.pz + M.c6 * t.vy - M.c7 * t.vx + M.c8 * t.pw + M.c9 * t.py -
                     M.c10 * t.px + M.c15 * t.pz + M.c26 * t.py - M.c27 * t.px +
                     M.c30 * t.vz;
    ctype const c9 = M.c0 * t.vx + M.c7 * t.mz - M.c8 * t.my + M.c10 * t.vz -
                     M.c11 * t.vy + M.c13 * t.pz - M.c14 * t.py - M.c26 * t.pw +
                     M.c29 * t.px + M.c30 * t.mx;
    ctype const c10 = M.c0 * t.vy - M.c6 * t.mz + M.c8 * t.mx - M.c9 * t.vz +
                      M.c11 * t.vx - M.c12 * t.pz + M.c14 * t.px - M.c27 * t.pw +
                      M.c29 * t.py + M.c30 * t.my;
    ctype const c11 = M.c0 * t.vz + M.c6 * t.my - M.c7 * t.mx + M.c9 * t.vy -
                      M.c10 * t.vx + M.c12 * t.py - M.c13 * t.px - M.c28 * t.pw +
                      M.c29 * t.pz + M.c30 * t.mz;
    ctype const c12 = M.c0 * t.mx + M.c10 * t.mz - M.c11 * t.my - M.c12 * t.pw +
                      M.c13 * t.vz - M.c14 * t.vy - M.c15 * t.mx - M.c27 * t.mz +
                      M.c28 * t.my + M.c29 * t.vx;
    ctype const c13 = M.c0 * t.my - M.c9 * t.mz + M.c11 * t.mx - M.c12 * t.vz -
                      M.c13 * t.pw + M.c14 * t.vx - M.c15 * t.my + M.c26 * t.mz -
                      M.c28 * t.mx + M.c29 * t.vy;
    ctype const c14 = M.c0 * t.mz + M.c9 * t.my - M.c10 * t.mx + M.c12 * t.vy -
                      M.c13 * t.vx - M.c14 * t.pw - M.c15 * t.mz - M.c26 * t.my +
                      M.c27 * t.mx + M.c29 * t.vz;
    ctype const c15 = -M.c0 * t.pw - M.c6 * t.mx - M.c7 * t.my - M.c8 * t.mz +
                      M.c12 * t.px + M.c13 * t.py + M.c14 * t.pz - M.c26 * t.vx -
                      M.c27 * t.vy - M.c28 * t.vz;
    ctype const c16 = M.c1 * t.pw + M.c4 * t.mx + M.c5 * t.px + M.c17 * t.vz -
                      M.c18 * t.vy + M.c20 * t.pz - M.c21 * t.py + M.c23 * t.mz -
                      M.c24 * t.my + M.c31 * t.vx;
    ctype const c17 = M.c2 * t.pw + M.c4 * t.my + M.c5 * t.py - M.c16 * t.vz +
                      M.c18 * t.vx - M.c19 * t.pz + M.c21 * t.px - M.c22 * t.mz +
                      M.c24 * t.mx + M.c31 * t.vy;
    ctype const c18 = M.c3 * t.pw + M.c4 * t.mz + M.c5 * t.pz + M.c16 * t.vy -
                      M.c17 * t.vx + M.c19 * t.py - M.c20 * t.px + M.c22 * t.my -
                      M.c23 * t.mx + M.c31 * t.vz;
    ctype const c19 = M.c2 * t.mz - M.c3 * t.my + M.c5 * t.vx + M.c17 * t.mz -
                      M.c18 * t.my - M.c19 * t.pw + M.c20 * t.vz - M.c21 * t.vy +
                      M.c25 * t.mx + M.c31 * t.mx;
    ctype const c20 = -M.c1 * t.mz + M.c3 * t.mx + M.c5 * t.vy - M.c16 * t.mz +
                      M.c18 * t.mx - M.c19 * t.vz - M.c20 * t.pw + M.c21 * t.vx +
                      M.c25 * t.my + M.c31 * t.my;
    ctype const c21 = M.c1 * t.my - M.c2 * t.mx + M.c5 * t.vz + M.c16 * t.my -
                      M.c17 * t.mx + M.c19 * t.vy - M.c20 * t.vx - M.c21 * t.pw +
                      M.c25 * t.mz + M.c31 * t.mz;
    ctype const c22 = -M.c2 * t.pz + M.c3 * t.py + M.c4 * t.vx + M.c17 * t.pz -
                      M.c18 * t.py + M.c22 * t.pw + M.c23 * t.vz - M.c24 * t.vy -
                      M.c25 * t.px + M.c31 * t.px;
    ctype const c23 = M.c1 * t.pz - M.c3 * t.px + M.c4 * t.vy - M.c16 * t.pz +
                      M.c18 * t.px - M.c22 * t.vz + M.c23 * t.pw + M.c24 * t.vx -
                      M.c25 * t.py + M.c31 * t.py;
    ctype const c24 = -M.c1 * t.py + M.c2 * t.px + M.c4 * t.vz + M.c16 * t.py -
                      M.c17 * t.px + M.c22 * t.vy - M.c23 * t.vx + M.c24 * t.pw -
                      M.c25 * t.pz + M.c31 * t.pz;
    ctype const c25 = -M.c1 * t.vx - M.c2 * t.vy - M.c3 * t.vz - M.c19 * t.px -
                      M.c20 * t.py - M.c21 * t.pz + M.c22 * t.mx + M.c23 * t.my +
                      M.c24 * t.mz + M.c31 * t.pw;
    ctype const c26 = M.c7 * t.mz - M.c8 * t.my - M.c9 * t.pw - M.c13 * t.pz +
                      M.c14 * t.py + M.c15 * t.vx + M.c27 * t.vz - M.c28 * t.vy -
                      M.c29 * t.px + M.c30 * t.mx;
    ctype const c27 = -M.c6 * t.mz + M.c8 * t.mx - M.c10 * t.pw + M.c12 * t.pz -
                      M.c14 * t.px + M.c15 * t.vy - M.c26 * t.vz + M.c28 * t.vx -
                      M.c29 * t.py + M.c30 * t.my;
    ctype const c28 = M.c6 * t.my - M.c7 * t.mx - M.c11 * t.pw - M.c12 * t.py +
                      M.c13 * t.px + M.c15 * t.vz + M.c26 * t.vy - M.c27 * t.vx -
                      M.c29 * t.pz + M.c30 * t.mz;
    ctype const c29 = -M.c9 * t.mx - M.c10 * t.my - M.c11 * t.mz - M.c12 * t.vx -
                      M.c13 * t.vy - M.c14 * t.vz + M.c26 * t.mx + M.c27 * t.my +
                      M.c28 * t.mz - M.c29 * t.pw;
    ctype const c30 = -M.c6 * t.vx - M.c7 * t.vy - M.c8 * t.vz - M.c9 * t.px -
                      M.c10 * t.py - M.c11 * t.pz - M.c26 * t.px - M.c27 * t.py -
                      M.c28 * t.pz + M.c30 * t.pw;
    ctype const c31 = -M.c16 * t.vx - M.c17 * t.vy - M.c18 * t.vz - M.c19 * t.px -
                      M.c20 * t.py - M.c21 * t.pz - M.c22 * t.mx - M.c23 * t.my -
                      M.c24 * t.mz + M.c25 * t.pw;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t,
                                                 MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c9 - t.vy * M.c10 - t.vz * M.c11 - t.mx * M.c6 -
                     t.my * M.c7 - t.mz * M.c8 - t.px * M.c12 - t.py * M.c13 -
                     t.pz * M.c14 - t.pw * M.c15;
    ctype const c1 = t.vx * M.c25 + t.vy * M.c3 - t.vz * M.c2 + t.mx * M.c4 +
                     t.my * M.c24 - t.mz * M.c23 - t.px * M.c5 - t.py * M.c21 +
                     t.pz * M.c20 + t.pw * M.c16;
    ctype const c2 = -t.vx * M.c3 + t.vy * M.c25 + t.vz * M.c1 - t.mx * M.c24 +
                     t.my * M.c4 + t.mz * M.c22 + t.px * M.c21 - t.py * M.c5 -
                     t.pz * M.c19 + t.pw * M.c17;
    ctype const c3 = t.vx * M.c2 - t.vy * M.c1 + t.vz * M.c25 + t.mx * M.c23 -
                     t.my * M.c22 + t.mz * M.c4 - t.px * M.c20 + t.py * M.c19 -
                     t.pz * M.c5 + t.pw * M.c18;
    ctype const c4 = -t.vx * M.c22 - t.vy * M.c23 - t.vz * M.c24 - t.px * M.c1 -
                     t.px * M.c16 - t.py * M.c2 - t.py * M.c17 - t.pz * M.c3 -
                     t.pz * M.c18 - t.pw * M.c4;
    ctype const c5 = -t.vx * M.c19 - t.vy * M.c20 - t.vz * M.c21 + t.mx * M.c1 -
                     t.mx * M.c16 + t.my * M.c2 - t.my * M.c17 + t.mz * M.c3 -
                     t.mz * M.c18 + t.pw * M.c5;
    ctype const c6 = t.vx * M.c30 + t.vy * M.c8 - t.vz * M.c7 + t.px * M.c0 -
                     t.px * M.c15 + t.py * M.c11 - t.py * M.c28 - t.pz * M.c10 +
                     t.pz * M.c27 - t.pw * M.c6;
    ctype const c7 = -t.vx * M.c8 + t.vy * M.c30 + t.vz * M.c6 - t.px * M.c11 +
                     t.px * M.c28 + t.py * M.c0 - t.py * M.c15 + t.pz * M.c9 -
                     t.pz * M.c26 - t.pw * M.c7;
    ctype const c8 = t.vx * M.c7 - t.vy * M.c6 + t.vz * M.c30 + t.px * M.c10 -
                     t.px * M.c27 - t.py * M.c9 + t.py * M.c26 + t.pz * M.c0 -
                     t.pz * M.c15 - t.pw * M.c8;
    ctype const c9 = t.vx * M.c0 + t.vy * M.c11 - t.vz * M.c10 + t.mx * M.c30 +
                     t.my * M.c8 - t.mz * M.c7 + t.px * M.c29 + t.py * M.c14 -
                     t.pz * M.c13 - t.pw * M.c26;
    ctype const c10 = -t.vx * M.c11 + t.vy * M.c0 + t.vz * M.c9 - t.mx * M.c8 +
                      t.my * M.c30 + t.mz * M.c6 - t.px * M.c14 + t.py * M.c29 +
                      t.pz * M.c12 - t.pw * M.c27;
    ctype const c11 = t.vx * M.c10 - t.vy * M.c9 + t.vz * M.c0 + t.mx * M.c7 -
                      t.my * M.c6 + t.mz * M.c30 + t.px * M.c13 - t.py * M.c12 +
                      t.pz * M.c29 - t.pw * M.c28;
    ctype const c12 = t.vx * M.c29 + t.vy * M.c14 - t.vz * M.c13 + t.mx * M.c0 +
                      t.mx * M.c15 + t.my * M.c11 + t.my * M.c28 - t.mz * M.c10 -
                      t.mz * M.c27 + t.pw * M.c12;
    ctype const c13 = -t.vx * M.c14 + t.vy * M.c29 + t.vz * M.c12 - t.mx * M.c11 -
                      t.mx * M.c28 + t.my * M.c0 + t.my * M.c15 + t.mz * M.c9 +
                      t.mz * M.c26 + t.pw * M.c13;
    ctype const c14 = t.vx * M.c13 - t.vy * M.c12 + t.vz * M.c29 + t.mx * M.c10 +
                      t.mx * M.c27 - t.my * M.c9 - t.my * M.c26 + t.mz * M.c0 +
                      t.mz * M.c15 + t.pw * M.c14;
    ctype const c15 = -t.vx * M.c26 - t.vy * M.c27 - t.vz * M.c28 + t.mx * M.c6 +
                      t.my * M.c7 + t.mz * M.c8 - t.px * M.c12 - t.py * M.c13 -
                      t.pz * M.c14 - t.pw * M.c0;
    ctype const c16 = t.vx * M.c31 + t.vy * M.c18 - t.vz * M.c17 + t.mx * M.c4 +
                      t.my * M.c24 - t.mz * M.c23 + t.px * M.c5 + t.py * M.c21 -
                      t.pz * M.c20 + t.pw * M.c1;
    ctype const c17 = -t.vx * M.c18 + t.vy * M.c31 + t.vz * M.c16 - t.mx * M.c24 +
                      t.my * M.c4 + t.mz * M.c22 - t.px * M.c21 + t.py * M.c5 +
                      t.pz * M.c19 + t.pw * M.c2;
    ctype const c18 = t.vx * M.c17 - t.vy * M.c16 + t.vz * M.c31 + t.mx * M.c23 -
                      t.my * M.c22 + t.mz * M.c4 + t.px * M.c20 - t.py * M.c19 +
                      t.pz * M.c5 + t.pw * M.c3;
    ctype const c19 = t.vx * M.c5 + t.vy * M.c21 - t.vz * M.c20 - t.mx * M.c25 +
                      t.mx * M.c31 - t.my * M.c3 + t.my * M.c18 + t.mz * M.c2 -
                      t.mz * M.c17 + t.pw * M.c19;
    ctype const c20 = -t.vx * M.c21 + t.vy * M.c5 + t.vz * M.c19 + t.mx * M.c3 -
                      t.mx * M.c18 - t.my * M.c25 + t.my * M.c31 - t.mz * M.c1 +
                      t.mz * M.c16 + t.pw * M.c20;
    ctype const c21 = t.vx * M.c20 - t.vy * M.c19 + t.vz * M.c5 - t.mx * M.c2 +
                      t.mx * M.c17 + t.my * M.c1 - t.my * M.c16 - t.mz * M.c25 +
                      t.mz * M.c31 + t.pw * M.c21;
    ctype const c22 = t.vx * M.c4 + t.vy * M.c24 - t.vz * M.c23 + t.px * M.c25 +
                      t.px * M.c31 + t.py * M.c3 + t.py * M.c18 - t.pz * M.c2 -
                      t.pz * M.c17 - t.pw * M.c22;
    ctype const c23 = -t.vx * M.c24 + t.vy * M.c4 + t.vz * M.c22 - t.px * M.c3 -
                      t.px * M.c18 + t.py * M.c25 + t.py * M.c31 + t.pz * M.c1 +
                      t.pz * M.c16 - t.pw * M.c23;
    ctype const c24 = t.vx * M.c23 - t.vy * M.c22 + t.vz * M.c4 + t.px * M.c2 +
                      t.px * M.c17 - t.py * M.c1 - t.py * M.c16 + t.pz * M.c25 +
                      t.pz * M.c31 - t.pw * M.c24;
    ctype const c25 = -t.vx * M.c1 - t.vy * M.c2 - t.vz * M.c3 - t.mx * M.c22 -
                      t.my * M.c23 - t.mz * M.c24 + t.px * M.c19 + t.py * M.c20 +
                      t.pz * M.c21 + t.pw * M.c31;
    ctype const c26 = t.vx * M.c15 + t.vy * M.c28 - t.vz * M.c27 - t.mx * M.c30 -
                      t.my * M.c8 + t.mz * M.c7 + t.px * M.c29 + t.py * M.c14 -
                      t.pz * M.c13 - t.pw * M.c9;
    ctype const c27 = -t.vx * M.c28 + t.vy * M.c15 + t.vz * M.c26 + t.mx * M.c8 -
                      t.my * M.c30 - t.mz * M.c6 - t.px * M.c14 + t.py * M.c29 +
                      t.pz * M.c12 - t.pw * M.c10;
    ctype const c28 = t.vx * M.c27 - t.vy * M.c26 + t.vz * M.c15 - t.mx * M.c7 +
                      t.my * M.c6 - t.mz * M.c30 + t.px * M.c13 - t.py * M.c12 +
                      t.pz * M.c29 - t.pw * M.c11;
    ctype const c29 = -t.vx * M.c12 - t.vy * M.c13 - t.vz * M.c14 - t.mx * M.c9 -
                      t.mx * M.c26 - t.my * M.c10 - t.my * M.c27 - t.mz * M.c11 -
                      t.mz * M.c28 + t.pw * M.c29;
    ctype const c30 = -t.vx * M.c6 - t.vy * M.c7 - t.vz * M.c8 - t.px * M.c9 +
                      t.px * M.c26 - t.py * M.c10 + t.py * M.c27 - t.pz * M.c11 +
                      t.pz * M.c28 - t.pw * M.c30;
    ctype const c31 = -t.vx * M.c16 - t.vy * M.c17 - t.vz * M.c18 - t.mx * M.c22 -
                      t.my * M.c23 - t.mz * M.c24 - t.px * M.c19 - t.py * M.c20 -
                      t.pz * M.c21 + t.pw * M.c25;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& M,
                                                 BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c16 * B.mx - M.c17 * B.my - M.c18 * B.mz - M.c19 * B.vx -
                     M.c20 * B.vy - M.c21 * B.vz - M.c22 * B.px - M.c23 * B.py -
                     M.c24 * B.pz - M.c25 * B.pw;
    ctype const c1 = M.c7 * B.pz - M.c8 * B.py + M.c9 * B.pw - M.c13 * B.vz +
                     M.c14 * B.vy + M.c15 * B.mx + M.c27 * B.mz - M.c28 * B.my -
                     M.c29 * B.vx + M.c30 * B.px;
    ctype const c2 = -M.c6 * B.pz + M.c8 * B.px + M.c10 * B.pw + M.c12 * B.vz -
                     M.c14 * B.vx + M.c15 * B.my - M.c26 * B.mz + M.c28 * B.mx -
                     M.c29 * B.vy + M.c30 * B.py;
    ctype const c3 = M.c6 * B.py - M.c7 * B.px + M.c11 * B.pw - M.c12 * B.vy +
                     M.c13 * B.vx + M.c15 * B.mz + M.c26 * B.my - M.c27 * B.mx -
                     M.c29 * B.vz + M.c30 * B.pz;
    ctype const c4 = M.c6 * B.mx + M.c7 * B.my + M.c8 * B.mz + M.c9 * B.vx +
                     M.c10 * B.vy + M.c11 * B.vz + M.c26 * B.vx + M.c27 * B.vy +
                     M.c28 * B.vz + M.c30 * B.pw;
    ctype const c5 = M.c9 * B.px + M.c10 * B.py + M.c11 * B.pz + M.c12 * B.mx +
                     M.c13 * B.my + M.c14 * B.mz - M.c26 * B.px - M.c27 * B.py -
                     M.c28 * B.pz - M.c29 * B.pw;
    ctype const c6 = -M.c2 * B.vz + M.c3 * B.vy + M.c4 * B.mx + M.c17 * B.vz -
                     M.c18 * B.vy - M.c22 * B.pw + M.c23 * B.mz - M.c24 * B.my -
                     M.c25 * B.vx + M.c31 * B.vx;
    ctype const c7 = M.c1 * B.vz - M.c3 * B.vx + M.c4 * B.my - M.c16 * B.vz +
                     M.c18 * B.vx - M.c22 * B.mz - M.c23 * B.pw + M.c24 * B.mx -
                     M.c25 * B.vy + M.c31 * B.vy;
    ctype const c8 = -M.c1 * B.vy + M.c2 * B.vx + M.c4 * B.mz + M.c16 * B.vy -
                     M.c17 * B.vx + M.c22 * B.my - M.c23 * B.mx - M.c24 * B.pw -
                     M.c25 * B.vz + M.c31 * B.vz;
    ctype const c9 = -M.c1 * B.pw + M.c4 * B.px + M.c5 * B.vx + M.c17 * B.mz -
                     M.c18 * B.my + M.c20 * B.vz - M.c21 * B.vy + M.c23 * B.pz -
                     M.c24 * B.py + M.c31 * B.mx;
    ctype const c10 = -M.c2 * B.pw + M.c4 * B.py + M.c5 * B.vy - M.c16 * B.mz +
                      M.c18 * B.mx - M.c19 * B.vz + M.c21 * B.vx - M.c22 * B.pz +
                      M.c24 * B.px + M.c31 * B.my;
    ctype const c11 = -M.c3 * B.pw + M.c4 * B.pz + M.c5 * B.vz + M.c16 * B.my -
                      M.c17 * B.mx + M.c19 * B.vy - M.c20 * B.vx + M.c22 * B.py -
                      M.c23 * B.px + M.c31 * B.mz;
    ctype const c12 = M.c2 * B.pz - M.c3 * B.py + M.c5 * B.mx + M.c17 * B.pz -
                      M.c18 * B.py + M.c19 * B.pw + M.c20 * B.mz - M.c21 * B.my +
                      M.c25 * B.px + M.c31 * B.px;
    ctype const c13 = -M.c1 * B.pz + M.c3 * B.px + M.c5 * B.my - M.c16 * B.pz +
                      M.c18 * B.px - M.c19 * B.mz + M.c20 * B.pw + M.c21 * B.mx +
                      M.c25 * B.py + M.c31 * B.py;
    ctype const c14 = M.c1 * B.py - M.c2 * B.px + M.c5 * B.mz + M.c16 * B.py -
                      M.c17 * B.px + M.c19 * B.my - M.c20 * B.mx + M.c21 * B.pw +
                      M.c25 * B.pz + M.c31 * B.pz;
    ctype const c15 = M.c1 * B.mx + M.c2 * B.my + M.c3 * B.mz + M.c19 * B.vx +
                      M.c20 * B.vy + M.c21 * B.vz - M.c22 * B.px - M.c23 * B.py -
                      M.c24 * B.pz + M.c31 * B.pw;
    ctype const c16 = -M.c0 * B.mx - M.c7 * B.pz + M.c8 * B.py - M.c10 * B.mz +
                      M.c11 * B.my - M.c13 * B.vz + M.c14 * B.vy - M.c26 * B.pw -
                      M.c29 * B.vx - M.c30 * B.px;
    ctype const c17 = -M.c0 * B.my + M.c6 * B.pz - M.c8 * B.px + M.c9 * B.mz -
                      M.c11 * B.mx + M.c12 * B.vz - M.c14 * B.vx - M.c27 * B.pw -
                      M.c29 * B.vy - M.c30 * B.py;
    ctype const c18 = -M.c0 * B.mz - M.c6 * B.py + M.c7 * B.px - M.c9 * B.my +
                      M.c10 * B.mx - M.c12 * B.vy + M.c13 * B.vx - M.c28 * B.pw -
                      M.c29 * B.vz - M.c30 * B.pz;
    ctype const c19 = -M.c0 * B.px - M.c10 * B.pz + M.c11 * B.py - M.c12 * B.pw -
                      M.c13 * B.mz + M.c14 * B.my + M.c15 * B.px + M.c27 * B.pz -
                      M.c28 * B.py - M.c29 * B.mx;
    ctype const c20 = -M.c0 * B.py + M.c9 * B.pz - M.c11 * B.px + M.c12 * B.mz -
                      M.c13 * B.pw - M.c14 * B.mx + M.c15 * B.py - M.c26 * B.pz +
                      M.c28 * B.px - M.c29 * B.my;
    ctype const c21 = -M.c0 * B.pz - M.c9 * B.py + M.c10 * B.px - M.c12 * B.my +
                      M.c13 * B.mx - M.c14 * B.pw + M.c15 * B.pz + M.c26 * B.py -
                      M.c27 * B.px - M.c29 * B.mz;
    ctype const c22 = -M.c0 * B.vx + M.c6 * B.pw - M.c7 * B.mz + M.c8 * B.my -
                      M.c10 * B.vz + M.c11 * B.vy - M.c15 * B.vx - M.c27 * B.vz +
                      M.c28 * B.vy - M.c30 * B.mx;
    ctype const c23 = -M.c0 * B.vy + M.c6 * B.mz + M.c7 * B.pw - M.c8 * B.mx +
                      M.c9 * B.vz - M.c11 * B.vx - M.c15 * B.vy + M.c26 * B.vz -
                      M.c28 * B.vx - M.c30 * B.my;
    ctype const c24 = -M.c0 * B.vz - M.c6 * B.my + M.c7 * B.mx + M.c8 * B.pw -
                      M.c9 * B.vy + M.c10 * B.vx - M.c15 * B.vz - M.c26 * B.vy +
                      M.c27 * B.vx - M.c30 * B.mz;
    ctype const c25 = M.c0 * B.pw - M.c6 * B.px - M.c7 * B.py - M.c8 * B.pz +
                      M.c12 * B.vx + M.c13 * B.vy + M.c14 * B.vz - M.c26 * B.mx -
                      M.c27 * B.my - M.c28 * B.mz;
    ctype const c26 = -M.c2 * B.mz + M.c3 * B.my + M.c4 * B.px - M.c5 * B.vx +
                      M.c16 * B.pw - M.c20 * B.vz + M.c21 * B.vy + M.c23 * B.pz -
                      M.c24 * B.py - M.c25 * B.mx;
    ctype const c27 = M.c1 * B.mz - M.c3 * B.mx + M.c4 * B.py - M.c5 * B.vy +
                      M.c17 * B.pw + M.c19 * B.vz - M.c21 * B.vx - M.c22 * B.pz +
                      M.c24 * B.px - M.c25 * B.my;
    ctype const c28 = -M.c1 * B.my + M.c2 * B.mx + M.c4 * B.pz - M.c5 * B.vz +
                      M.c18 * B.pw - M.c19 * B.vy + M.c20 * B.vx + M.c22 * B.py -
                      M.c23 * B.px - M.c25 * B.mz;
    ctype const c29 = -M.c1 * B.px - M.c2 * B.py - M.c3 * B.pz + M.c5 * B.pw -
                      M.c16 * B.px - M.c17 * B.py - M.c18 * B.pz - M.c19 * B.mx -
                      M.c20 * B.my - M.c21 * B.mz;
    ctype const c30 = M.c1 * B.vx + M.c2 * B.vy + M.c3 * B.vz - M.c4 * B.pw -
                      M.c16 * B.vx - M.c17 * B.vy - M.c18 * B.vz - M.c22 * B.mx -
                      M.c23 * B.my - M.c24 * B.mz;
    ctype const c31 = M.c6 * B.px + M.c7 * B.py + M.c8 * B.pz + M.c9 * B.mx +
                      M.c10 * B.my + M.c11 * B.mz + M.c12 * B.vx + M.c13 * B.vy +
                      M.c14 * B.vz - M.c15 * B.pw;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B,
                                                 MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c19 - B.vy * M.c20 - B.vz * M.c21 - B.mx * M.c16 -
                     B.my * M.c17 - B.mz * M.c18 - B.px * M.c22 - B.py * M.c23 -
                     B.pz * M.c24 - B.pw * M.c25;
    ctype const c1 = B.vx * M.c29 + B.vy * M.c14 - B.vz * M.c13 + B.mx * M.c15 +
                     B.my * M.c28 - B.mz * M.c27 - B.px * M.c30 - B.py * M.c8 +
                     B.pz * M.c7 + B.pw * M.c9;
    ctype const c2 = -B.vx * M.c14 + B.vy * M.c29 + B.vz * M.c12 - B.mx * M.c28 +
                     B.my * M.c15 + B.mz * M.c26 + B.px * M.c8 - B.py * M.c30 -
                     B.pz * M.c6 + B.pw * M.c10;
    ctype const c3 = B.vx * M.c13 - B.vy * M.c12 + B.vz * M.c29 + B.mx * M.c27 -
                     B.my * M.c26 + B.mz * M.c15 - B.px * M.c7 + B.py * M.c6 -
                     B.pz * M.c30 + B.pw * M.c11;
    ctype const c4 = B.vx * M.c9 - B.vx * M.c26 + B.vy * M.c10 - B.vy * M.c27 +
                     B.vz * M.c11 - B.vz * M.c28 + B.mx * M.c6 + B.my * M.c7 +
                     B.mz * M.c8 - B.pw * M.c30;
    ctype const c5 = B.mx * M.c12 + B.my * M.c13 + B.mz * M.c14 + B.px * M.c9 +
                     B.px * M.c26 + B.py * M.c10 + B.py * M.c27 + B.pz * M.c11 +
                     B.pz * M.c28 + B.pw * M.c29;
    ctype const c6 = B.vx * M.c25 + B.vx * M.c31 + B.vy * M.c3 + B.vy * M.c18 -
                     B.vz * M.c2 - B.vz * M.c17 + B.mx * M.c4 + B.my * M.c24 -
                     B.mz * M.c23 + B.pw * M.c22;
    ctype const c7 = -B.vx * M.c3 - B.vx * M.c18 + B.vy * M.c25 + B.vy * M.c31 +
                     B.vz * M.c1 + B.vz * M.c16 - B.mx * M.c24 + B.my * M.c4 +
                     B.mz * M.c22 + B.pw * M.c23;
    ctype const c8 = B.vx * M.c2 + B.vx * M.c17 - B.vy * M.c1 - B.vy * M.c16 +
                     B.vz * M.c25 + B.vz * M.c31 + B.mx * M.c23 - B.my * M.c22 +
                     B.mz * M.c4 + B.pw * M.c24;
    ctype const c9 = B.vx * M.c5 + B.vy * M.c21 - B.vz * M.c20 + B.mx * M.c31 +
                     B.my * M.c18 - B.mz * M.c17 + B.px * M.c4 + B.py * M.c24 -
                     B.pz * M.c23 - B.pw * M.c1;
    ctype const c10 = -B.vx * M.c21 + B.vy * M.c5 + B.vz * M.c19 - B.mx * M.c18 +
                      B.my * M.c31 + B.mz * M.c16 - B.px * M.c24 + B.py * M.c4 +
                      B.pz * M.c22 - B.pw * M.c2;
    ctype const c11 = B.vx * M.c20 - B.vy * M.c19 + B.vz * M.c5 + B.mx * M.c17 -
                      B.my * M.c16 + B.mz * M.c31 + B.px * M.c23 - B.py * M.c22 +
                      B.pz * M.c4 - B.pw * M.c3;
    ctype const c12 = B.mx * M.c5 + B.my * M.c21 - B.mz * M.c20 - B.px * M.c25 +
                      B.px * M.c31 - B.py * M.c3 + B.py * M.c18 + B.pz * M.c2 -
                      B.pz * M.c17 - B.pw * M.c19;
    ctype const c13 = -B.mx * M.c21 + B.my * M.c5 + B.mz * M.c19 + B.px * M.c3 -
                      B.px * M.c18 - B.py * M.c25 + B.py * M.c31 - B.pz * M.c1 +
                      B.pz * M.c16 - B.pw * M.c20;
    ctype const c14 = B.mx * M.c20 - B.my * M.c19 + B.mz * M.c5 - B.px * M.c2 +
                      B.px * M.c17 + B.py * M.c1 - B.py * M.c16 - B.pz * M.c25 +
                      B.pz * M.c31 - B.pw * M.c21;
    ctype const c15 = -B.vx * M.c19 - B.vy * M.c20 - B.vz * M.c21 + B.mx * M.c1 +
                      B.my * M.c2 + B.mz * M.c3 + B.px * M.c22 + B.py * M.c23 +
                      B.pz * M.c24 + B.pw * M.c31;
    ctype const c16 = -B.vx * M.c29 - B.vy * M.c14 + B.vz * M.c13 - B.mx * M.c0 -
                      B.my * M.c11 + B.mz * M.c10 - B.px * M.c30 - B.py * M.c8 +
                      B.pz * M.c7 - B.pw * M.c26;
    ctype const c17 = B.vx * M.c14 - B.vy * M.c29 - B.vz * M.c12 + B.mx * M.c11 -
                      B.my * M.c0 - B.mz * M.c9 + B.px * M.c8 - B.py * M.c30 -
                      B.pz * M.c6 - B.pw * M.c27;
    ctype const c18 = -B.vx * M.c13 + B.vy * M.c12 - B.vz * M.c29 - B.mx * M.c10 +
                      B.my * M.c9 - B.mz * M.c0 - B.px * M.c7 + B.py * M.c6 -
                      B.pz * M.c30 - B.pw * M.c28;
    ctype const c19 = -B.mx * M.c29 - B.my * M.c14 + B.mz * M.c13 - B.px * M.c0 -
                      B.px * M.c15 - B.py * M.c11 - B.py * M.c28 + B.pz * M.c10 +
                      B.pz * M.c27 + B.pw * M.c12;
    ctype const c20 = B.mx * M.c14 - B.my * M.c29 - B.mz * M.c12 + B.px * M.c11 +
                      B.px * M.c28 - B.py * M.c0 - B.py * M.c15 - B.pz * M.c9 -
                      B.pz * M.c26 + B.pw * M.c13;
    ctype const c21 = -B.mx * M.c13 + B.my * M.c12 - B.mz * M.c29 - B.px * M.c10 -
                      B.px * M.c27 + B.py * M.c9 + B.py * M.c26 - B.pz * M.c0 -
                      B.pz * M.c15 + B.pw * M.c14;
    ctype const c22 = -B.vx * M.c0 + B.vx * M.c15 - B.vy * M.c11 + B.vy * M.c28 +
                      B.vz * M.c10 - B.vz * M.c27 - B.mx * M.c30 - B.my * M.c8 +
                      B.mz * M.c7 - B.pw * M.c6;
    ctype const c23 = B.vx * M.c11 - B.vx * M.c28 - B.vy * M.c0 + B.vy * M.c15 -
                      B.vz * M.c9 + B.vz * M.c26 + B.mx * M.c8 - B.my * M.c30 -
                      B.mz * M.c6 - B.pw * M.c7;
    ctype const c24 = -B.vx * M.c10 + B.vx * M.c27 + B.vy * M.c9 - B.vy * M.c26 -
                      B.vz * M.c0 + B.vz * M.c15 - B.mx * M.c7 + B.my * M.c6 -
                      B.mz * M.c30 - B.pw * M.c8;
    ctype const c25 = -B.vx * M.c12 - B.vy * M.c13 - B.vz * M.c14 - B.mx * M.c26 -
                      B.my * M.c27 - B.mz * M.c28 + B.px * M.c6 + B.py * M.c7 +
                      B.pz * M.c8 + B.pw * M.c0;
    ctype const c26 = B.vx * M.c5 + B.vy * M.c21 - B.vz * M.c20 - B.mx * M.c25 -
                      B.my * M.c3 + B.mz * M.c2 - B.px * M.c4 - B.py * M.c24 +
                      B.pz * M.c23 + B.pw * M.c16;
    ctype const c27 = -B.vx * M.c21 + B.vy * M.c5 + B.vz * M.c19 + B.mx * M.c3 -
                      B.my * M.c25 - B.mz * M.c1 + B.px * M.c24 - B.py * M.c4 -
                      B.pz * M.c22 + B.pw * M.c17;
    ctype const c28 = B.vx * M.c20 - B.vy * M.c19 + B.vz * M.c5 - B.mx * M.c2 +
                      B.my * M.c1 - B.mz * M.c25 - B.px * M.c23 + B.py * M.c22 -
                      B.pz * M.c4 + B.pw * M.c18;
    ctype const c29 = -B.mx * M.c19 - B.my * M.c20 - B.mz * M.c21 + B.px * M.c1 -
                      B.px * M.c16 + B.py * M.c2 - B.py * M.c17 + B.pz * M.c3 -
                      B.pz * M.c18 - B.pw * M.c5;
    ctype const c30 = -B.vx * M.c1 - B.vx * M.c16 - B.vy * M.c2 - B.vy * M.c17 -
                      B.vz * M.c3 - B.vz * M.c18 - B.mx * M.c22 - B.my * M.c23 -
                      B.mz * M.c24 + B.pw * M.c4;
    ctype const c31 = B.vx * M.c12 + B.vy * M.c13 + B.vz * M.c14 + B.mx * M.c9 +
                      B.my * M.c10 + B.mz * M.c11 + B.px * M.c6 + B.py * M.c7 +
                      B.pz * M.c8 - B.pw * M.c15;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& M, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c26 * v.x + M.c27 * v.y + M.c28 * v.z + M.c29 * v.w + M.c30 * v.u;
    ctype const c1 = M.c17 * v.z - M.c18 * v.y + M.c19 * v.w - M.c22 * v.u + M.c31 * v.x;
    ctype const c2 = -M.c16 * v.z + M.c18 * v.x + M.c20 * v.w - M.c23 * v.u + M.c31 * v.y;
    ctype const c3 = M.c16 * v.y - M.c17 * v.x + M.c21 * v.w - M.c24 * v.u + M.c31 * v.z;
    ctype const c4 = -M.c22 * v.x - M.c23 * v.y - M.c24 * v.z - M.c25 * v.w + M.c31 * v.w;
    ctype const c5 = M.c19 * v.x + M.c20 * v.y + M.c21 * v.z + M.c25 * v.u + M.c31 * v.u;
    ctype const c6 = M.c7 * v.z - M.c8 * v.y + M.c9 * v.w + M.c26 * v.w + M.c30 * v.x;
    ctype const c7 = -M.c6 * v.z + M.c8 * v.x + M.c10 * v.w + M.c27 * v.w + M.c30 * v.y;
    ctype const c8 = M.c6 * v.y - M.c7 * v.x + M.c11 * v.w + M.c28 * v.w + M.c30 * v.z;
    ctype const c9 = M.c6 * v.u + M.c12 * v.w - M.c15 * v.x - M.c27 * v.z + M.c28 * v.y;
    ctype const c10 = M.c7 * v.u + M.c13 * v.w - M.c15 * v.y + M.c26 * v.z - M.c28 * v.x;
    ctype const c11 = M.c8 * v.u + M.c14 * v.w - M.c15 * v.z - M.c26 * v.y + M.c27 * v.x;
    ctype const c12 = M.c9 * v.u - M.c13 * v.z + M.c14 * v.y - M.c26 * v.u - M.c29 * v.x;
    ctype const c13 = M.c10 * v.u + M.c12 * v.z - M.c14 * v.x - M.c27 * v.u - M.c29 * v.y;
    ctype const c14 = M.c11 * v.u - M.c12 * v.y + M.c13 * v.x - M.c28 * v.u - M.c29 * v.z;
    ctype const c15 = M.c9 * v.x + M.c10 * v.y + M.c11 * v.z - M.c29 * v.w + M.c30 * v.u;
    ctype const c16 = M.c2 * v.z - M.c3 * v.y + M.c19 * v.w + M.c22 * v.u + M.c25 * v.x;
    ctype const c17 = -M.c1 * v.z + M.c3 * v.x + M.c20 * v.w + M.c23 * v.u + M.c25 * v.y;
    ctype const c18 = M.c1 * v.y - M.c2 * v.x + M.c21 * v.w + M.c24 * v.u + M.c25 * v.z;
    ctype const c19 = M.c1 * v.u - M.c5 * v.x + M.c16 * v.u - M.c20 * v.z + M.c21 * v.y;
    ctype const c20 = M.c2 * v.u - M.c5 * v.y + M.c17 * v.u + M.c19 * v.z - M.c21 * v.x;
    ctype const c21 = M.c3 * v.u - M.c5 * v.z + M.c18 * v.u - M.c19 * v.y + M.c20 * v.x;
    ctype const c22 = -M.c1 * v.w + M.c4 * v.x + M.c16 * v.w + M.c23 * v.z - M.c24 * v.y;
    ctype const c23 = -M.c2 * v.w + M.c4 * v.y + M.c17 * v.w - M.c22 * v.z + M.c24 * v.x;
    ctype const c24 = -M.c3 * v.w + M.c4 * v.z + M.c18 * v.w + M.c22 * v.y - M.c23 * v.x;
    ctype const c25 = -M.c4 * v.u + M.c5 * v.w - M.c16 * v.x - M.c17 * v.y - M.c18 * v.z;
    ctype const c26 = -M.c0 * v.x + M.c6 * v.u - M.c10 * v.z + M.c11 * v.y - M.c12 * v.w;
    ctype const c27 = -M.c0 * v.y + M.c7 * v.u + M.c9 * v.z - M.c11 * v.x - M.c13 * v.w;
    ctype const c28 = -M.c0 * v.z + M.c8 * v.u - M.c9 * v.y + M.c10 * v.x - M.c14 * v.w;
    ctype const c29 = M.c0 * v.u + M.c12 * v.x + M.c13 * v.y + M.c14 * v.z - M.c15 * v.u;
    ctype const c30 = M.c0 * v.w - M.c6 * v.x - M.c7 * v.y - M.c8 * v.z + M.c15 * v.w;
    ctype const c31 = -M.c1 * v.x - M.c2 * v.y - M.c3 * v.z + M.c4 * v.u + M.c5 * v.w;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c26 + v.y * M.c27 + v.z * M.c28 + v.w * M.c29 + v.u * M.c30;
    ctype const c1 = v.x * M.c31 + v.y * M.c18 - v.z * M.c17 - v.w * M.c19 + v.u * M.c22;
    ctype const c2 = -v.x * M.c18 + v.y * M.c31 + v.z * M.c16 - v.w * M.c20 + v.u * M.c23;
    ctype const c3 = v.x * M.c17 - v.y * M.c16 + v.z * M.c31 - v.w * M.c21 + v.u * M.c24;
    ctype const c4 = v.x * M.c22 + v.y * M.c23 + v.z * M.c24 + v.w * M.c25 + v.w * M.c31;
    ctype const c5 = -v.x * M.c19 - v.y * M.c20 - v.z * M.c21 - v.u * M.c25 + v.u * M.c31;
    ctype const c6 = -v.x * M.c30 - v.y * M.c8 + v.z * M.c7 + v.w * M.c9 - v.w * M.c26;
    ctype const c7 = v.x * M.c8 - v.y * M.c30 - v.z * M.c6 + v.w * M.c10 - v.w * M.c27;
    ctype const c8 = -v.x * M.c7 + v.y * M.c6 - v.z * M.c30 + v.w * M.c11 - v.w * M.c28;
    ctype const c9 = -v.x * M.c15 - v.y * M.c28 + v.z * M.c27 + v.w * M.c12 + v.u * M.c6;
    ctype const c10 = v.x * M.c28 - v.y * M.c15 - v.z * M.c26 + v.w * M.c13 + v.u * M.c7;
    ctype const c11 = -v.x * M.c27 + v.y * M.c26 - v.z * M.c15 + v.w * M.c14 + v.u * M.c8;
    ctype const c12 = v.x * M.c29 + v.y * M.c14 - v.z * M.c13 + v.u * M.c9 + v.u * M.c26;
    ctype const c13 =
        -v.x * M.c14 + v.y * M.c29 + v.z * M.c12 + v.u * M.c10 + v.u * M.c27;
    ctype const c14 = v.x * M.c13 - v.y * M.c12 + v.z * M.c29 + v.u * M.c11 + v.u * M.c28;
    ctype const c15 = v.x * M.c9 + v.y * M.c10 + v.z * M.c11 + v.w * M.c29 - v.u * M.c30;
    ctype const c16 = v.x * M.c25 + v.y * M.c3 - v.z * M.c2 + v.w * M.c19 + v.u * M.c22;
    ctype const c17 = -v.x * M.c3 + v.y * M.c25 + v.z * M.c1 + v.w * M.c20 + v.u * M.c23;
    ctype const c18 = v.x * M.c2 - v.y * M.c1 + v.z * M.c25 + v.w * M.c21 + v.u * M.c24;
    ctype const c19 = v.x * M.c5 + v.y * M.c21 - v.z * M.c20 - v.u * M.c1 + v.u * M.c16;
    ctype const c20 = -v.x * M.c21 + v.y * M.c5 + v.z * M.c19 - v.u * M.c2 + v.u * M.c17;
    ctype const c21 = v.x * M.c20 - v.y * M.c19 + v.z * M.c5 - v.u * M.c3 + v.u * M.c18;
    ctype const c22 = -v.x * M.c4 - v.y * M.c24 + v.z * M.c23 + v.w * M.c1 + v.w * M.c16;
    ctype const c23 = v.x * M.c24 - v.y * M.c4 - v.z * M.c22 + v.w * M.c2 + v.w * M.c17;
    ctype const c24 = -v.x * M.c23 + v.y * M.c22 - v.z * M.c4 + v.w * M.c3 + v.w * M.c18;
    ctype const c25 = -v.x * M.c16 - v.y * M.c17 - v.z * M.c18 - v.w * M.c5 + v.u * M.c4;
    ctype const c26 = -v.x * M.c0 - v.y * M.c11 + v.z * M.c10 + v.w * M.c12 - v.u * M.c6;
    ctype const c27 = v.x * M.c11 - v.y * M.c0 - v.z * M.c9 + v.w * M.c13 - v.u * M.c7;
    ctype const c28 = -v.x * M.c10 + v.y * M.c9 - v.z * M.c0 + v.w * M.c14 - v.u * M.c8;
    ctype const c29 = -v.x * M.c12 - v.y * M.c13 - v.z * M.c14 + v.u * M.c0 + v.u * M.c15;
    ctype const c30 = v.x * M.c6 + v.y * M.c7 + v.z * M.c8 + v.w * M.c0 - v.w * M.c15;
    ctype const c31 = -v.x * M.c1 - v.y * M.c2 - v.z * M.c3 + v.w * M.c5 + v.u * M.c4;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(MVec3dc<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c31 * ctype(s);
    ctype const c1 = M.c26 * ctype(s);
    ctype const c2 = M.c27 * ctype(s);
    ctype const c3 = M.c28 * ctype(s);
    ctype const c4 = -M.c30 * ctype(s);
    ctype const c5 = -M.c29 * ctype(s);
    ctype const c6 = M.c22 * ctype(s);
    ctype const c7 = M.c23 * ctype(s);
    ctype const c8 = M.c24 * ctype(s);
    ctype const c9 = M.c16 * ctype(s);
    ctype const c10 = M.c17 * ctype(s);
    ctype const c11 = M.c18 * ctype(s);
    ctype const c12 = M.c19 * ctype(s);
    ctype const c13 = M.c20 * ctype(s);
    ctype const c14 = M.c21 * ctype(s);
    ctype const c15 = -M.c25 * ctype(s);
    ctype const c16 = -M.c9 * ctype(s);
    ctype const c17 = -M.c10 * ctype(s);
    ctype const c18 = -M.c11 * ctype(s);
    ctype const c19 = -M.c12 * ctype(s);
    ctype const c20 = -M.c13 * ctype(s);
    ctype const c21 = -M.c14 * ctype(s);
    ctype const c22 = -M.c6 * ctype(s);
    ctype const c23 = -M.c7 * ctype(s);
    ctype const c24 = -M.c8 * ctype(s);
    ctype const c25 = M.c15 * ctype(s);
    ctype const c26 = -M.c1 * ctype(s);
    ctype const c27 = -M.c2 * ctype(s);
    ctype const c28 = -M.c3 * ctype(s);
    ctype const c29 = M.c5 * ctype(s);
    ctype const c30 = M.c4 * ctype(s);
    ctype const c31 = -M.c0 * ctype(s);
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c31;
    ctype const c1 = ctype(s) * M.c26;
    ctype const c2 = ctype(s) * M.c27;
    ctype const c3 = ctype(s) * M.c28;
    ctype const c4 = -ctype(s) * M.c30;
    ctype const c5 = -ctype(s) * M.c29;
    ctype const c6 = ctype(s) * M.c22;
    ctype const c7 = ctype(s) * M.c23;
    ctype const c8 = ctype(s) * M.c24;
    ctype const c9 = ctype(s) * M.c16;
    ctype const c10 = ctype(s) * M.c17;
    ctype const c11 = ctype(s) * M.c18;
    ctype const c12 = ctype(s) * M.c19;
    ctype const c13 = ctype(s) * M.c20;
    ctype const c14 = ctype(s) * M.c21;
    ctype const c15 = -ctype(s) * M.c25;
    ctype const c16 = -ctype(s) * M.c9;
    ctype const c17 = -ctype(s) * M.c10;
    ctype const c18 = -ctype(s) * M.c11;
    ctype const c19 = -ctype(s) * M.c12;
    ctype const c20 = -ctype(s) * M.c13;
    ctype const c21 = -ctype(s) * M.c14;
    ctype const c22 = -ctype(s) * M.c6;
    ctype const c23 = -ctype(s) * M.c7;
    ctype const c24 = -ctype(s) * M.c8;
    ctype const c25 = ctype(s) * M.c15;
    ctype const c26 = -ctype(s) * M.c1;
    ctype const c27 = -ctype(s) * M.c2;
    ctype const c28 = -ctype(s) * M.c3;
    ctype const c29 = ctype(s) * M.c5;
    ctype const c30 = ctype(s) * M.c4;
    ctype const c31 = -ctype(s) * M.c0;
    return MVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27,
                          c28, c29, c30, c31);
}

// cga3dc rgpr :: rgpr(mv_e,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& A,
                                                   MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c11 + A.c1 * B.c14 + A.c2 * B.c9 - A.c3 * B.c8 +
                     A.c4 * B.c10 + A.c5 * B.c13 - A.c6 * B.c12 - A.c7 * B.c15 -
                     A.c8 * B.c3 + A.c9 * B.c2 + A.c10 * B.c4 + A.c11 * B.c0 +
                     A.c12 * B.c6 - A.c13 * B.c5 - A.c14 * B.c1 + A.c15 * B.c7;
    ctype const c1 = A.c0 * B.c12 - A.c1 * B.c9 + A.c2 * B.c14 + A.c3 * B.c7 -
                     A.c4 * B.c13 + A.c5 * B.c10 + A.c6 * B.c11 + A.c7 * B.c3 -
                     A.c8 * B.c15 - A.c9 * B.c1 + A.c10 * B.c5 - A.c11 * B.c6 +
                     A.c12 * B.c0 + A.c13 * B.c4 - A.c14 * B.c2 + A.c15 * B.c8;
    ctype const c2 = A.c0 * B.c13 + A.c1 * B.c8 - A.c2 * B.c7 + A.c3 * B.c14 +
                     A.c4 * B.c12 - A.c5 * B.c11 + A.c6 * B.c10 - A.c7 * B.c2 +
                     A.c8 * B.c1 - A.c9 * B.c15 + A.c10 * B.c6 + A.c11 * B.c5 -
                     A.c12 * B.c4 + A.c13 * B.c0 - A.c14 * B.c3 + A.c15 * B.c9;
    ctype const c3 = -A.c0 * B.c15 + A.c1 * B.c4 - A.c1 * B.c11 + A.c2 * B.c5 -
                     A.c2 * B.c12 + A.c3 * B.c6 - A.c3 * B.c13 + A.c4 * B.c1 +
                     A.c5 * B.c2 + A.c6 * B.c3 - A.c10 * B.c15 + A.c11 * B.c1 +
                     A.c12 * B.c2 + A.c13 * B.c3 - A.c15 * B.c0 + A.c15 * B.c10;
    ctype const c4 = -A.c0 * B.c14 + A.c4 * B.c7 + A.c5 * B.c8 + A.c6 * B.c9 +
                     A.c7 * B.c4 + A.c7 * B.c11 + A.c8 * B.c5 + A.c8 * B.c12 +
                     A.c9 * B.c6 + A.c9 * B.c13 + A.c10 * B.c14 - A.c11 * B.c7 -
                     A.c12 * B.c8 - A.c13 * B.c9 - A.c14 * B.c0 - A.c14 * B.c10;
    ctype const c5 = -A.c0 * B.c4 - A.c1 * B.c14 - A.c2 * B.c9 + A.c3 * B.c8 -
                     A.c4 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c15 -
                     A.c8 * B.c3 + A.c9 * B.c2 - A.c10 * B.c11 - A.c11 * B.c10 -
                     A.c12 * B.c13 + A.c13 * B.c12 - A.c14 * B.c1 - A.c15 * B.c7;
    ctype const c6 = -A.c0 * B.c5 + A.c1 * B.c9 - A.c2 * B.c14 - A.c3 * B.c7 +
                     A.c4 * B.c6 - A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c3 -
                     A.c8 * B.c15 - A.c9 * B.c1 - A.c10 * B.c12 + A.c11 * B.c13 -
                     A.c12 * B.c10 - A.c13 * B.c11 - A.c14 * B.c2 - A.c15 * B.c8;
    ctype const c7 = -A.c0 * B.c6 - A.c1 * B.c8 + A.c2 * B.c7 - A.c3 * B.c14 -
                     A.c4 * B.c5 + A.c5 * B.c4 - A.c6 * B.c0 - A.c7 * B.c2 + A.c8 * B.c1 -
                     A.c9 * B.c15 - A.c10 * B.c13 - A.c11 * B.c12 + A.c12 * B.c11 -
                     A.c13 * B.c10 - A.c14 * B.c3 - A.c15 * B.c9;
    ctype const c8 = -A.c0 * B.c7 - A.c4 * B.c14 - A.c5 * B.c9 + A.c6 * B.c8 -
                     A.c7 * B.c0 - A.c7 * B.c10 - A.c8 * B.c6 - A.c8 * B.c13 +
                     A.c9 * B.c5 + A.c9 * B.c12 + A.c10 * B.c7 + A.c11 * B.c14 +
                     A.c12 * B.c9 - A.c13 * B.c8 - A.c14 * B.c4 - A.c14 * B.c11;
    ctype const c9 = -A.c0 * B.c8 + A.c4 * B.c9 - A.c5 * B.c14 - A.c6 * B.c7 +
                     A.c7 * B.c6 + A.c7 * B.c13 - A.c8 * B.c0 - A.c8 * B.c10 -
                     A.c9 * B.c4 - A.c9 * B.c11 + A.c10 * B.c8 - A.c11 * B.c9 +
                     A.c12 * B.c14 + A.c13 * B.c7 - A.c14 * B.c5 - A.c14 * B.c12;
    ctype const c10 = -A.c0 * B.c9 - A.c4 * B.c8 + A.c5 * B.c7 - A.c6 * B.c14 -
                      A.c7 * B.c5 - A.c7 * B.c12 + A.c8 * B.c4 + A.c8 * B.c11 -
                      A.c9 * B.c0 - A.c9 * B.c10 + A.c10 * B.c9 + A.c11 * B.c8 -
                      A.c12 * B.c7 + A.c13 * B.c14 - A.c14 * B.c6 - A.c14 * B.c13;
    ctype const c11 = -A.c0 * B.c1 - A.c1 * B.c0 + A.c1 * B.c10 - A.c2 * B.c6 +
                      A.c2 * B.c13 + A.c3 * B.c5 - A.c3 * B.c12 - A.c4 * B.c15 -
                      A.c5 * B.c3 + A.c6 * B.c2 - A.c10 * B.c1 - A.c11 * B.c15 -
                      A.c12 * B.c3 + A.c13 * B.c2 - A.c15 * B.c4 + A.c15 * B.c11;
    ctype const c12 = -A.c0 * B.c2 + A.c1 * B.c6 - A.c1 * B.c13 - A.c2 * B.c0 +
                      A.c2 * B.c10 - A.c3 * B.c4 + A.c3 * B.c11 + A.c4 * B.c3 -
                      A.c5 * B.c15 - A.c6 * B.c1 - A.c10 * B.c2 + A.c11 * B.c3 -
                      A.c12 * B.c15 - A.c13 * B.c1 - A.c15 * B.c5 + A.c15 * B.c12;
    ctype const c13 = -A.c0 * B.c3 - A.c1 * B.c5 + A.c1 * B.c12 + A.c2 * B.c4 -
                      A.c2 * B.c11 - A.c3 * B.c0 + A.c3 * B.c10 - A.c4 * B.c2 +
                      A.c5 * B.c1 - A.c6 * B.c15 - A.c10 * B.c3 - A.c11 * B.c2 +
                      A.c12 * B.c1 - A.c13 * B.c15 - A.c15 * B.c6 + A.c15 * B.c13;
    ctype const c14 = A.c0 * B.c10 - A.c1 * B.c7 - A.c2 * B.c8 - A.c3 * B.c9 -
                      A.c4 * B.c11 - A.c5 * B.c12 - A.c6 * B.c13 + A.c7 * B.c1 +
                      A.c8 * B.c2 + A.c9 * B.c3 + A.c10 * B.c0 - A.c11 * B.c4 -
                      A.c12 * B.c5 - A.c13 * B.c6 - A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = -A.c0 * B.c0 + A.c1 * B.c7 + A.c2 * B.c8 + A.c3 * B.c9 +
                      A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c1 +
                      A.c8 * B.c2 + A.c9 * B.c3 - A.c10 * B.c10 + A.c11 * B.c11 +
                      A.c12 * B.c12 + A.c13 * B.c13 - A.c14 * B.c15 - A.c15 * B.c14;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_e,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& A,
                                                   MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                     A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 - A.c7 * B.c11 -
                     A.c8 * B.c12 - A.c9 * B.c13 - A.c10 * B.c14 + A.c11 * B.c0 +
                     A.c12 * B.c1 + A.c13 * B.c2 + A.c14 * B.c3 + A.c15 * B.c4;
    ctype const c1 = A.c0 * B.c11 + A.c1 * B.c14 + A.c1 * B.c15 + A.c2 * B.c2 +
                     A.c2 * B.c7 - A.c3 * B.c1 - A.c3 * B.c6 + A.c4 * B.c3 +
                     A.c5 * B.c13 - A.c6 * B.c12 + A.c10 * B.c11 + A.c11 * B.c3 +
                     A.c12 * B.c13 - A.c13 * B.c12 + A.c15 * B.c0 + A.c15 * B.c5;
    ctype const c2 = A.c0 * B.c12 - A.c1 * B.c2 - A.c1 * B.c7 + A.c2 * B.c14 +
                     A.c2 * B.c15 + A.c3 * B.c0 + A.c3 * B.c5 - A.c4 * B.c13 +
                     A.c5 * B.c3 + A.c6 * B.c11 + A.c10 * B.c12 - A.c11 * B.c13 +
                     A.c12 * B.c3 + A.c13 * B.c11 + A.c15 * B.c1 + A.c15 * B.c6;
    ctype const c3 = A.c0 * B.c13 + A.c1 * B.c1 + A.c1 * B.c6 - A.c2 * B.c0 -
                     A.c2 * B.c5 + A.c3 * B.c14 + A.c3 * B.c15 + A.c4 * B.c12 -
                     A.c5 * B.c11 + A.c6 * B.c3 + A.c10 * B.c13 + A.c11 * B.c12 -
                     A.c12 * B.c11 + A.c13 * B.c3 + A.c15 * B.c2 + A.c15 * B.c7;
    ctype const c4 = A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c10 - A.c3 * B.c9 +
                     A.c4 * B.c15 + A.c5 * B.c7 - A.c6 * B.c6 + A.c7 * B.c3 +
                     A.c8 * B.c13 - A.c9 * B.c12 - A.c10 * B.c0 - A.c11 * B.c14 -
                     A.c12 * B.c2 + A.c13 * B.c1 + A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c5 = A.c0 * B.c6 - A.c1 * B.c10 + A.c2 * B.c4 + A.c3 * B.c8 -
                     A.c4 * B.c7 + A.c5 * B.c15 + A.c6 * B.c5 - A.c7 * B.c13 +
                     A.c8 * B.c3 + A.c9 * B.c11 - A.c10 * B.c1 + A.c11 * B.c2 -
                     A.c12 * B.c14 - A.c13 * B.c0 + A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c6 = A.c0 * B.c7 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c4 + A.c4 * B.c6 -
                     A.c5 * B.c5 + A.c6 * B.c15 + A.c7 * B.c12 - A.c8 * B.c11 +
                     A.c9 * B.c3 - A.c10 * B.c2 - A.c11 * B.c1 + A.c12 * B.c0 -
                     A.c13 * B.c14 + A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c7 = A.c0 * B.c8 + A.c4 * B.c4 + A.c5 * B.c10 - A.c6 * B.c9 -
                     A.c7 * B.c14 + A.c7 * B.c15 - A.c8 * B.c2 + A.c8 * B.c7 +
                     A.c9 * B.c1 - A.c9 * B.c6 - A.c10 * B.c8 - A.c11 * B.c4 -
                     A.c12 * B.c10 + A.c13 * B.c9 - A.c14 * B.c0 + A.c14 * B.c5;
    ctype const c8 = A.c0 * B.c9 - A.c4 * B.c10 + A.c5 * B.c4 + A.c6 * B.c8 +
                     A.c7 * B.c2 - A.c7 * B.c7 - A.c8 * B.c14 + A.c8 * B.c15 -
                     A.c9 * B.c0 + A.c9 * B.c5 - A.c10 * B.c9 + A.c11 * B.c10 -
                     A.c12 * B.c4 - A.c13 * B.c8 - A.c14 * B.c1 + A.c14 * B.c6;
    ctype const c9 = A.c0 * B.c10 + A.c4 * B.c9 - A.c5 * B.c8 + A.c6 * B.c4 -
                     A.c7 * B.c1 + A.c7 * B.c6 + A.c8 * B.c0 - A.c8 * B.c5 -
                     A.c9 * B.c14 + A.c9 * B.c15 - A.c10 * B.c10 - A.c11 * B.c9 +
                     A.c12 * B.c8 - A.c13 * B.c4 - A.c14 * B.c2 + A.c14 * B.c7;
    ctype const c10 = -A.c0 * B.c14 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 +
                      A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2 + A.c7 * B.c11 +
                      A.c8 * B.c12 + A.c9 * B.c13 + A.c10 * B.c15 - A.c11 * B.c5 -
                      A.c12 * B.c6 - A.c13 * B.c7 - A.c14 * B.c3 + A.c15 * B.c4;
    ctype const c11 = -A.c0 * B.c0 + A.c1 * B.c4 + A.c2 * B.c10 - A.c3 * B.c9 -
                      A.c4 * B.c14 - A.c5 * B.c2 + A.c6 * B.c1 - A.c7 * B.c3 -
                      A.c8 * B.c13 + A.c9 * B.c12 + A.c10 * B.c5 + A.c11 * B.c15 +
                      A.c12 * B.c7 - A.c13 * B.c6 - A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c12 = -A.c0 * B.c1 - A.c1 * B.c10 + A.c2 * B.c4 + A.c3 * B.c8 +
                      A.c4 * B.c2 - A.c5 * B.c14 - A.c6 * B.c0 + A.c7 * B.c13 -
                      A.c8 * B.c3 - A.c9 * B.c11 + A.c10 * B.c6 - A.c11 * B.c7 +
                      A.c12 * B.c15 + A.c13 * B.c5 - A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c13 = -A.c0 * B.c2 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c4 -
                      A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c14 - A.c7 * B.c12 +
                      A.c8 * B.c11 - A.c9 * B.c3 + A.c10 * B.c7 + A.c11 * B.c6 -
                      A.c12 * B.c5 + A.c13 * B.c15 - A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c14 = A.c0 * B.c4 - A.c4 * B.c8 - A.c5 * B.c9 - A.c6 * B.c10 +
                      A.c7 * B.c0 - A.c7 * B.c5 + A.c8 * B.c1 - A.c8 * B.c6 +
                      A.c9 * B.c2 - A.c9 * B.c7 - A.c10 * B.c4 + A.c11 * B.c8 +
                      A.c12 * B.c9 + A.c13 * B.c10 - A.c14 * B.c14 + A.c14 * B.c15;
    ctype const c15 = A.c0 * B.c3 - A.c1 * B.c0 - A.c1 * B.c5 - A.c2 * B.c1 -
                      A.c2 * B.c6 - A.c3 * B.c2 - A.c3 * B.c7 - A.c4 * B.c11 -
                      A.c5 * B.c12 - A.c6 * B.c13 + A.c10 * B.c3 - A.c11 * B.c11 -
                      A.c12 * B.c12 - A.c13 * B.c13 + A.c15 * B.c14 + A.c15 * B.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(mv_u,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& A,
                                                   MVec3dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                     A.c4 * B.c15 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 -
                     A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 - A.c11 * B.c7 -
                     A.c12 * B.c8 - A.c13 * B.c9 - A.c14 * B.c10 + A.c15 * B.c0;
    ctype const c1 = -A.c0 * B.c15 - A.c1 * B.c3 + A.c2 * B.c2 + A.c3 * B.c4 -
                     A.c3 * B.c11 + A.c5 * B.c15 + A.c6 * B.c3 - A.c7 * B.c2 +
                     A.c11 * B.c0 - A.c11 * B.c10 + A.c12 * B.c6 - A.c12 * B.c13 -
                     A.c13 * B.c5 + A.c13 * B.c12 - A.c14 * B.c1 + A.c15 * B.c1;
    ctype const c2 = A.c0 * B.c3 - A.c1 * B.c15 - A.c2 * B.c1 + A.c3 * B.c5 -
                     A.c3 * B.c12 - A.c5 * B.c3 + A.c6 * B.c15 + A.c7 * B.c1 -
                     A.c11 * B.c6 + A.c11 * B.c13 + A.c12 * B.c0 - A.c12 * B.c10 +
                     A.c13 * B.c4 - A.c13 * B.c11 - A.c14 * B.c2 + A.c15 * B.c2;
    ctype const c3 = -A.c0 * B.c2 + A.c1 * B.c1 - A.c2 * B.c15 + A.c3 * B.c6 -
                     A.c3 * B.c13 + A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c15 +
                     A.c11 * B.c5 - A.c11 * B.c12 - A.c12 * B.c4 + A.c12 * B.c11 +
                     A.c13 * B.c0 - A.c13 * B.c10 - A.c14 * B.c3 + A.c15 * B.c3;
    ctype const c4 = -A.c0 * B.c10 - A.c1 * B.c13 + A.c2 * B.c12 + A.c3 * B.c7 +
                     A.c4 * B.c1 + A.c5 * B.c0 + A.c6 * B.c6 - A.c7 * B.c5 +
                     A.c8 * B.c15 + A.c9 * B.c3 - A.c10 * B.c2 + A.c11 * B.c14 +
                     A.c12 * B.c9 - A.c13 * B.c8 - A.c14 * B.c11 + A.c15 * B.c4;
    ctype const c5 = A.c0 * B.c13 - A.c1 * B.c10 - A.c2 * B.c11 + A.c3 * B.c8 +
                     A.c4 * B.c2 - A.c5 * B.c6 + A.c6 * B.c0 + A.c7 * B.c4 - A.c8 * B.c3 +
                     A.c9 * B.c15 + A.c10 * B.c1 - A.c11 * B.c9 + A.c12 * B.c14 +
                     A.c13 * B.c7 - A.c14 * B.c12 + A.c15 * B.c5;
    ctype const c6 = -A.c0 * B.c12 + A.c1 * B.c11 - A.c2 * B.c10 + A.c3 * B.c9 +
                     A.c4 * B.c3 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c0 + A.c8 * B.c2 -
                     A.c9 * B.c1 + A.c10 * B.c15 + A.c11 * B.c8 - A.c12 * B.c7 +
                     A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c6;
    ctype const c7 = A.c0 * B.c14 + A.c1 * B.c9 - A.c2 * B.c8 + A.c4 * B.c4 +
                     A.c4 * B.c11 + A.c5 * B.c14 + A.c6 * B.c9 - A.c7 * B.c8 +
                     A.c8 * B.c0 + A.c8 * B.c10 + A.c9 * B.c6 + A.c9 * B.c13 -
                     A.c10 * B.c5 - A.c10 * B.c12 + A.c14 * B.c7 + A.c15 * B.c7;
    ctype const c8 = -A.c0 * B.c9 + A.c1 * B.c14 + A.c2 * B.c7 + A.c4 * B.c5 +
                     A.c4 * B.c12 - A.c5 * B.c9 + A.c6 * B.c14 + A.c7 * B.c7 -
                     A.c8 * B.c6 - A.c8 * B.c13 + A.c9 * B.c0 + A.c9 * B.c10 +
                     A.c10 * B.c4 + A.c10 * B.c11 + A.c14 * B.c8 + A.c15 * B.c8;
    ctype const c9 = A.c0 * B.c8 - A.c1 * B.c7 + A.c2 * B.c14 + A.c4 * B.c6 +
                     A.c4 * B.c13 + A.c5 * B.c8 - A.c6 * B.c7 + A.c7 * B.c14 +
                     A.c8 * B.c5 + A.c8 * B.c12 - A.c9 * B.c4 - A.c9 * B.c11 +
                     A.c10 * B.c0 + A.c10 * B.c10 + A.c14 * B.c9 + A.c15 * B.c9;
    ctype const c10 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c14 -
                      A.c4 * B.c15 - A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 +
                      A.c8 * B.c1 + A.c9 * B.c2 + A.c10 * B.c3 - A.c11 * B.c7 -
                      A.c12 * B.c8 - A.c13 * B.c9 - A.c14 * B.c0 + A.c15 * B.c10;
    ctype const c11 = -A.c0 * B.c0 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c7 -
                      A.c4 * B.c1 + A.c5 * B.c10 + A.c6 * B.c13 - A.c7 * B.c12 -
                      A.c8 * B.c15 - A.c9 * B.c3 + A.c10 * B.c2 + A.c11 * B.c14 +
                      A.c12 * B.c9 - A.c13 * B.c8 - A.c14 * B.c4 + A.c15 * B.c11;
    ctype const c12 = A.c0 * B.c6 - A.c1 * B.c0 - A.c2 * B.c4 + A.c3 * B.c8 -
                      A.c4 * B.c2 - A.c5 * B.c13 + A.c6 * B.c10 + A.c7 * B.c11 +
                      A.c8 * B.c3 - A.c9 * B.c15 - A.c10 * B.c1 - A.c11 * B.c9 +
                      A.c12 * B.c14 + A.c13 * B.c7 - A.c14 * B.c5 + A.c15 * B.c12;
    ctype const c13 = -A.c0 * B.c5 + A.c1 * B.c4 - A.c2 * B.c0 + A.c3 * B.c9 -
                      A.c4 * B.c3 + A.c5 * B.c12 - A.c6 * B.c11 + A.c7 * B.c10 -
                      A.c8 * B.c2 + A.c9 * B.c1 - A.c10 * B.c15 + A.c11 * B.c8 -
                      A.c12 * B.c7 + A.c13 * B.c14 - A.c14 * B.c6 + A.c15 * B.c13;
    ctype const c14 = -A.c0 * B.c7 - A.c1 * B.c8 - A.c2 * B.c9 + A.c4 * B.c0 +
                      A.c4 * B.c10 - A.c5 * B.c7 - A.c6 * B.c8 - A.c7 * B.c9 -
                      A.c8 * B.c4 - A.c8 * B.c11 - A.c9 * B.c5 - A.c9 * B.c12 -
                      A.c10 * B.c6 - A.c10 * B.c13 + A.c14 * B.c14 + A.c15 * B.c14;
    ctype const c15 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0 -
                      A.c3 * B.c10 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 -
                      A.c11 * B.c4 + A.c11 * B.c11 - A.c12 * B.c5 + A.c12 * B.c12 -
                      A.c13 * B.c6 + A.c13 * B.c13 - A.c14 * B.c15 + A.c15 * B.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& M,
                                                   PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    ctype const c4 = M.c4 * ctype(ps);
    ctype const c5 = M.c5 * ctype(ps);
    ctype const c6 = M.c6 * ctype(ps);
    ctype const c7 = M.c7 * ctype(ps);
    ctype const c8 = M.c8 * ctype(ps);
    ctype const c9 = M.c9 * ctype(ps);
    ctype const c10 = M.c10 * ctype(ps);
    ctype const c11 = M.c11 * ctype(ps);
    ctype const c12 = M.c12 * ctype(ps);
    ctype const c13 = M.c13 * ctype(ps);
    ctype const c14 = M.c14 * ctype(ps);
    ctype const c15 = M.c15 * ctype(ps);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    ctype const c4 = ctype(ps) * M.c4;
    ctype const c5 = ctype(ps) * M.c5;
    ctype const c6 = ctype(ps) * M.c6;
    ctype const c7 = ctype(ps) * M.c7;
    ctype const c8 = ctype(ps) * M.c8;
    ctype const c9 = ctype(ps) * M.c9;
    ctype const c10 = ctype(ps) * M.c10;
    ctype const c11 = ctype(ps) * M.c11;
    ctype const c12 = ctype(ps) * M.c12;
    ctype const c13 = ctype(ps) * M.c13;
    ctype const c14 = ctype(ps) * M.c14;
    ctype const c15 = ctype(ps) * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(mv_e,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& M,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * Q.x + M.c1 * Q.w + M.c5 * Q.z - M.c6 * Q.y - M.c7 * Q.u;
    ctype const c1 = M.c0 * Q.y + M.c2 * Q.w - M.c4 * Q.z + M.c6 * Q.x - M.c8 * Q.u;
    ctype const c2 = M.c0 * Q.z + M.c3 * Q.w + M.c4 * Q.y - M.c5 * Q.x - M.c9 * Q.u;
    ctype const c3 = -M.c0 * Q.u - M.c1 * Q.x - M.c2 * Q.y - M.c3 * Q.z - M.c10 * Q.u;
    ctype const c4 = -M.c0 * Q.w + M.c7 * Q.x + M.c8 * Q.y + M.c9 * Q.z + M.c10 * Q.w;
    ctype const c5 = -M.c1 * Q.w - M.c7 * Q.u - M.c10 * Q.x - M.c12 * Q.z + M.c13 * Q.y;
    ctype const c6 = -M.c2 * Q.w - M.c8 * Q.u - M.c10 * Q.y + M.c11 * Q.z - M.c13 * Q.x;
    ctype const c7 = -M.c3 * Q.w - M.c9 * Q.u - M.c10 * Q.z - M.c11 * Q.y + M.c12 * Q.x;
    ctype const c8 = -M.c4 * Q.w - M.c8 * Q.z + M.c9 * Q.y + M.c11 * Q.w - M.c14 * Q.x;
    ctype const c9 = -M.c5 * Q.w + M.c7 * Q.z - M.c9 * Q.x + M.c12 * Q.w - M.c14 * Q.y;
    ctype const c10 = -M.c6 * Q.w - M.c7 * Q.y + M.c8 * Q.x + M.c13 * Q.w - M.c14 * Q.z;
    ctype const c11 = M.c2 * Q.z - M.c3 * Q.y - M.c4 * Q.u - M.c11 * Q.u + M.c15 * Q.x;
    ctype const c12 = -M.c1 * Q.z + M.c3 * Q.x - M.c5 * Q.u - M.c12 * Q.u + M.c15 * Q.y;
    ctype const c13 = M.c1 * Q.y - M.c2 * Q.x - M.c6 * Q.u - M.c13 * Q.u + M.c15 * Q.z;
    ctype const c14 = -M.c4 * Q.x - M.c5 * Q.y - M.c6 * Q.z - M.c14 * Q.u + M.c15 * Q.w;
    ctype const c15 = M.c11 * Q.x + M.c12 * Q.y + M.c13 * Q.z - M.c14 * Q.u - M.c15 * Q.w;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(quadvec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * M.c0 + Q.y * M.c6 - Q.z * M.c5 - Q.w * M.c1 + Q.u * M.c7;
    ctype const c1 = -Q.x * M.c6 + Q.y * M.c0 + Q.z * M.c4 - Q.w * M.c2 + Q.u * M.c8;
    ctype const c2 = Q.x * M.c5 - Q.y * M.c4 + Q.z * M.c0 - Q.w * M.c3 + Q.u * M.c9;
    ctype const c3 = Q.x * M.c1 + Q.y * M.c2 + Q.z * M.c3 - Q.u * M.c0 + Q.u * M.c10;
    ctype const c4 = -Q.x * M.c7 - Q.y * M.c8 - Q.z * M.c9 - Q.w * M.c0 - Q.w * M.c10;
    ctype const c5 = -Q.x * M.c10 - Q.y * M.c13 + Q.z * M.c12 - Q.w * M.c1 - Q.u * M.c7;
    ctype const c6 = Q.x * M.c13 - Q.y * M.c10 - Q.z * M.c11 - Q.w * M.c2 - Q.u * M.c8;
    ctype const c7 = -Q.x * M.c12 + Q.y * M.c11 - Q.z * M.c10 - Q.w * M.c3 - Q.u * M.c9;
    ctype const c8 = Q.x * M.c14 + Q.y * M.c9 - Q.z * M.c8 - Q.w * M.c4 - Q.w * M.c11;
    ctype const c9 = -Q.x * M.c9 + Q.y * M.c14 + Q.z * M.c7 - Q.w * M.c5 - Q.w * M.c12;
    ctype const c10 = Q.x * M.c8 - Q.y * M.c7 + Q.z * M.c14 - Q.w * M.c6 - Q.w * M.c13;
    ctype const c11 = -Q.x * M.c15 - Q.y * M.c3 + Q.z * M.c2 - Q.u * M.c4 + Q.u * M.c11;
    ctype const c12 = Q.x * M.c3 - Q.y * M.c15 - Q.z * M.c1 - Q.u * M.c5 + Q.u * M.c12;
    ctype const c13 = -Q.x * M.c2 + Q.y * M.c1 - Q.z * M.c15 - Q.u * M.c6 + Q.u * M.c13;
    ctype const c14 = -Q.x * M.c4 - Q.y * M.c5 - Q.z * M.c6 - Q.w * M.c15 + Q.u * M.c14;
    ctype const c15 = Q.x * M.c11 + Q.y * M.c12 + Q.z * M.c13 - Q.w * M.c15 - Q.u * M.c14;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_e,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& M,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz - M.c4 * t.vx -
                     M.c5 * t.vy - M.c6 * t.vz - M.c7 * t.px - M.c8 * t.py - M.c9 * t.pz -
                     M.c10 * t.pw;
    ctype const c1 = M.c0 * t.px + M.c1 * t.pw + M.c2 * t.vz - M.c3 * t.vy + M.c5 * t.pz -
                     M.c6 * t.py + M.c10 * t.px + M.c12 * t.pz - M.c13 * t.py +
                     M.c15 * t.vx;
    ctype const c2 = M.c0 * t.py - M.c1 * t.vz + M.c2 * t.pw + M.c3 * t.vx - M.c4 * t.pz +
                     M.c6 * t.px + M.c10 * t.py - M.c11 * t.pz + M.c13 * t.px +
                     M.c15 * t.vy;
    ctype const c3 = M.c0 * t.pz + M.c1 * t.vy - M.c2 * t.vx + M.c3 * t.pw + M.c4 * t.py -
                     M.c5 * t.px + M.c10 * t.pz + M.c11 * t.py - M.c12 * t.px +
                     M.c15 * t.vz;
    ctype const c4 = M.c0 * t.vx + M.c2 * t.mz - M.c3 * t.my + M.c5 * t.vz - M.c6 * t.vy +
                     M.c8 * t.pz - M.c9 * t.py - M.c11 * t.pw + M.c14 * t.px +
                     M.c15 * t.mx;
    ctype const c5 = M.c0 * t.vy - M.c1 * t.mz + M.c3 * t.mx - M.c4 * t.vz + M.c6 * t.vx -
                     M.c7 * t.pz + M.c9 * t.px - M.c12 * t.pw + M.c14 * t.py +
                     M.c15 * t.my;
    ctype const c6 = M.c0 * t.vz + M.c1 * t.my - M.c2 * t.mx + M.c4 * t.vy - M.c5 * t.vx +
                     M.c7 * t.py - M.c8 * t.px - M.c13 * t.pw + M.c14 * t.pz +
                     M.c15 * t.mz;
    ctype const c7 = M.c0 * t.mx + M.c5 * t.mz - M.c6 * t.my - M.c7 * t.pw + M.c8 * t.vz -
                     M.c9 * t.vy - M.c10 * t.mx - M.c12 * t.mz + M.c13 * t.my +
                     M.c14 * t.vx;
    ctype const c8 = M.c0 * t.my - M.c4 * t.mz + M.c6 * t.mx - M.c7 * t.vz - M.c8 * t.pw +
                     M.c9 * t.vx - M.c10 * t.my + M.c11 * t.mz - M.c13 * t.mx +
                     M.c14 * t.vy;
    ctype const c9 = M.c0 * t.mz + M.c4 * t.my - M.c5 * t.mx + M.c7 * t.vy - M.c8 * t.vx -
                     M.c9 * t.pw - M.c10 * t.mz - M.c11 * t.my + M.c12 * t.mx +
                     M.c14 * t.vz;
    ctype const c10 = -M.c0 * t.pw - M.c1 * t.mx - M.c2 * t.my - M.c3 * t.mz +
                      M.c7 * t.px + M.c8 * t.py + M.c9 * t.pz - M.c11 * t.vx -
                      M.c12 * t.vy - M.c13 * t.vz;
    ctype const c11 = M.c2 * t.mz - M.c3 * t.my - M.c4 * t.pw - M.c8 * t.pz +
                      M.c9 * t.py + M.c10 * t.vx + M.c12 * t.vz - M.c13 * t.vy -
                      M.c14 * t.px + M.c15 * t.mx;
    ctype const c12 = -M.c1 * t.mz + M.c3 * t.mx - M.c5 * t.pw + M.c7 * t.pz -
                      M.c9 * t.px + M.c10 * t.vy - M.c11 * t.vz + M.c13 * t.vx -
                      M.c14 * t.py + M.c15 * t.my;
    ctype const c13 = M.c1 * t.my - M.c2 * t.mx - M.c6 * t.pw - M.c7 * t.py +
                      M.c8 * t.px + M.c10 * t.vz + M.c11 * t.vy - M.c12 * t.vx -
                      M.c14 * t.pz + M.c15 * t.mz;
    ctype const c14 = -M.c4 * t.mx - M.c5 * t.my - M.c6 * t.mz - M.c7 * t.vx -
                      M.c8 * t.vy - M.c9 * t.vz + M.c11 * t.mx + M.c12 * t.my +
                      M.c13 * t.mz - M.c14 * t.pw;
    ctype const c15 = -M.c1 * t.vx - M.c2 * t.vy - M.c3 * t.vz - M.c4 * t.px -
                      M.c5 * t.py - M.c6 * t.pz - M.c11 * t.px - M.c12 * t.py -
                      M.c13 * t.pz + M.c15 * t.pw;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(trivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * M.c4 - t.vy * M.c5 - t.vz * M.c6 - t.mx * M.c1 -
                     t.my * M.c2 - t.mz * M.c3 - t.px * M.c7 - t.py * M.c8 - t.pz * M.c9 -
                     t.pw * M.c10;
    ctype const c1 = t.vx * M.c15 + t.vy * M.c3 - t.vz * M.c2 + t.px * M.c0 -
                     t.px * M.c10 + t.py * M.c6 - t.py * M.c13 - t.pz * M.c5 +
                     t.pz * M.c12 - t.pw * M.c1;
    ctype const c2 = -t.vx * M.c3 + t.vy * M.c15 + t.vz * M.c1 - t.px * M.c6 +
                     t.px * M.c13 + t.py * M.c0 - t.py * M.c10 + t.pz * M.c4 -
                     t.pz * M.c11 - t.pw * M.c2;
    ctype const c3 = t.vx * M.c2 - t.vy * M.c1 + t.vz * M.c15 + t.px * M.c5 -
                     t.px * M.c12 - t.py * M.c4 + t.py * M.c11 + t.pz * M.c0 -
                     t.pz * M.c10 - t.pw * M.c3;
    ctype const c4 = t.vx * M.c0 + t.vy * M.c6 - t.vz * M.c5 + t.mx * M.c15 +
                     t.my * M.c3 - t.mz * M.c2 + t.px * M.c14 + t.py * M.c9 -
                     t.pz * M.c8 - t.pw * M.c11;
    ctype const c5 = -t.vx * M.c6 + t.vy * M.c0 + t.vz * M.c4 - t.mx * M.c3 +
                     t.my * M.c15 + t.mz * M.c1 - t.px * M.c9 + t.py * M.c14 +
                     t.pz * M.c7 - t.pw * M.c12;
    ctype const c6 = t.vx * M.c5 - t.vy * M.c4 + t.vz * M.c0 + t.mx * M.c2 - t.my * M.c1 +
                     t.mz * M.c15 + t.px * M.c8 - t.py * M.c7 + t.pz * M.c14 -
                     t.pw * M.c13;
    ctype const c7 = t.vx * M.c14 + t.vy * M.c9 - t.vz * M.c8 + t.mx * M.c0 +
                     t.mx * M.c10 + t.my * M.c6 + t.my * M.c13 - t.mz * M.c5 -
                     t.mz * M.c12 + t.pw * M.c7;
    ctype const c8 = -t.vx * M.c9 + t.vy * M.c14 + t.vz * M.c7 - t.mx * M.c6 -
                     t.mx * M.c13 + t.my * M.c0 + t.my * M.c10 + t.mz * M.c4 +
                     t.mz * M.c11 + t.pw * M.c8;
    ctype const c9 = t.vx * M.c8 - t.vy * M.c7 + t.vz * M.c14 + t.mx * M.c5 +
                     t.mx * M.c12 - t.my * M.c4 - t.my * M.c11 + t.mz * M.c0 +
                     t.mz * M.c10 + t.pw * M.c9;
    ctype const c10 = -t.vx * M.c11 - t.vy * M.c12 - t.vz * M.c13 + t.mx * M.c1 +
                      t.my * M.c2 + t.mz * M.c3 - t.px * M.c7 - t.py * M.c8 -
                      t.pz * M.c9 - t.pw * M.c0;
    ctype const c11 = t.vx * M.c10 + t.vy * M.c13 - t.vz * M.c12 - t.mx * M.c15 -
                      t.my * M.c3 + t.mz * M.c2 + t.px * M.c14 + t.py * M.c9 -
                      t.pz * M.c8 - t.pw * M.c4;
    ctype const c12 = -t.vx * M.c13 + t.vy * M.c10 + t.vz * M.c11 + t.mx * M.c3 -
                      t.my * M.c15 - t.mz * M.c1 - t.px * M.c9 + t.py * M.c14 +
                      t.pz * M.c7 - t.pw * M.c5;
    ctype const c13 = t.vx * M.c12 - t.vy * M.c11 + t.vz * M.c10 - t.mx * M.c2 +
                      t.my * M.c1 - t.mz * M.c15 + t.px * M.c8 - t.py * M.c7 +
                      t.pz * M.c14 - t.pw * M.c6;
    ctype const c14 = -t.vx * M.c7 - t.vy * M.c8 - t.vz * M.c9 - t.mx * M.c4 -
                      t.mx * M.c11 - t.my * M.c5 - t.my * M.c12 - t.mz * M.c6 -
                      t.mz * M.c13 + t.pw * M.c14;
    ctype const c15 = -t.vx * M.c1 - t.vy * M.c2 - t.vz * M.c3 - t.px * M.c4 +
                      t.px * M.c11 - t.py * M.c5 + t.py * M.c12 - t.pz * M.c6 +
                      t.pz * M.c13 - t.pw * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(mv_e,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& M,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * B.pz - M.c3 * B.py + M.c4 * B.pw - M.c8 * B.vz + M.c9 * B.vy +
                     M.c10 * B.mx + M.c12 * B.mz - M.c13 * B.my - M.c14 * B.vx +
                     M.c15 * B.px;
    ctype const c1 = -M.c1 * B.pz + M.c3 * B.px + M.c5 * B.pw + M.c7 * B.vz -
                     M.c9 * B.vx + M.c10 * B.my - M.c11 * B.mz + M.c13 * B.mx -
                     M.c14 * B.vy + M.c15 * B.py;
    ctype const c2 = M.c1 * B.py - M.c2 * B.px + M.c6 * B.pw - M.c7 * B.vy + M.c8 * B.vx +
                     M.c10 * B.mz + M.c11 * B.my - M.c12 * B.mx - M.c14 * B.vz +
                     M.c15 * B.pz;
    ctype const c3 = M.c1 * B.mx + M.c2 * B.my + M.c3 * B.mz + M.c4 * B.vx + M.c5 * B.vy +
                     M.c6 * B.vz + M.c11 * B.vx + M.c12 * B.vy + M.c13 * B.vz +
                     M.c15 * B.pw;
    ctype const c4 = M.c4 * B.px + M.c5 * B.py + M.c6 * B.pz + M.c7 * B.mx + M.c8 * B.my +
                     M.c9 * B.mz - M.c11 * B.px - M.c12 * B.py - M.c13 * B.pz -
                     M.c14 * B.pw;
    ctype const c5 = -M.c0 * B.mx - M.c2 * B.pz + M.c3 * B.py - M.c5 * B.mz +
                     M.c6 * B.my - M.c8 * B.vz + M.c9 * B.vy - M.c11 * B.pw -
                     M.c14 * B.vx - M.c15 * B.px;
    ctype const c6 = -M.c0 * B.my + M.c1 * B.pz - M.c3 * B.px + M.c4 * B.mz -
                     M.c6 * B.mx + M.c7 * B.vz - M.c9 * B.vx - M.c12 * B.pw -
                     M.c14 * B.vy - M.c15 * B.py;
    ctype const c7 = -M.c0 * B.mz - M.c1 * B.py + M.c2 * B.px - M.c4 * B.my +
                     M.c5 * B.mx - M.c7 * B.vy + M.c8 * B.vx - M.c13 * B.pw -
                     M.c14 * B.vz - M.c15 * B.pz;
    ctype const c8 = -M.c0 * B.px - M.c5 * B.pz + M.c6 * B.py - M.c7 * B.pw -
                     M.c8 * B.mz + M.c9 * B.my + M.c10 * B.px + M.c12 * B.pz -
                     M.c13 * B.py - M.c14 * B.mx;
    ctype const c9 = -M.c0 * B.py + M.c4 * B.pz - M.c6 * B.px + M.c7 * B.mz -
                     M.c8 * B.pw - M.c9 * B.mx + M.c10 * B.py - M.c11 * B.pz +
                     M.c13 * B.px - M.c14 * B.my;
    ctype const c10 = -M.c0 * B.pz - M.c4 * B.py + M.c5 * B.px - M.c7 * B.my +
                      M.c8 * B.mx - M.c9 * B.pw + M.c10 * B.pz + M.c11 * B.py -
                      M.c12 * B.px - M.c14 * B.mz;
    ctype const c11 = -M.c0 * B.vx + M.c1 * B.pw - M.c2 * B.mz + M.c3 * B.my -
                      M.c5 * B.vz + M.c6 * B.vy - M.c10 * B.vx - M.c12 * B.vz +
                      M.c13 * B.vy - M.c15 * B.mx;
    ctype const c12 = -M.c0 * B.vy + M.c1 * B.mz + M.c2 * B.pw - M.c3 * B.mx +
                      M.c4 * B.vz - M.c6 * B.vx - M.c10 * B.vy + M.c11 * B.vz -
                      M.c13 * B.vx - M.c15 * B.my;
    ctype const c13 = -M.c0 * B.vz - M.c1 * B.my + M.c2 * B.mx + M.c3 * B.pw -
                      M.c4 * B.vy + M.c5 * B.vx - M.c10 * B.vz - M.c11 * B.vy +
                      M.c12 * B.vx - M.c15 * B.mz;
    ctype const c14 = M.c0 * B.pw - M.c1 * B.px - M.c2 * B.py - M.c3 * B.pz +
                      M.c7 * B.vx + M.c8 * B.vy + M.c9 * B.vz - M.c11 * B.mx -
                      M.c12 * B.my - M.c13 * B.mz;
    ctype const c15 = M.c1 * B.px + M.c2 * B.py + M.c3 * B.pz + M.c4 * B.mx +
                      M.c5 * B.my + M.c6 * B.mz + M.c7 * B.vx + M.c8 * B.vy +
                      M.c9 * B.vz - M.c10 * B.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(bivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c14 + B.vy * M.c9 - B.vz * M.c8 + B.mx * M.c10 +
                     B.my * M.c13 - B.mz * M.c12 - B.px * M.c15 - B.py * M.c3 +
                     B.pz * M.c2 + B.pw * M.c4;
    ctype const c1 = -B.vx * M.c9 + B.vy * M.c14 + B.vz * M.c7 - B.mx * M.c13 +
                     B.my * M.c10 + B.mz * M.c11 + B.px * M.c3 - B.py * M.c15 -
                     B.pz * M.c1 + B.pw * M.c5;
    ctype const c2 = B.vx * M.c8 - B.vy * M.c7 + B.vz * M.c14 + B.mx * M.c12 -
                     B.my * M.c11 + B.mz * M.c10 - B.px * M.c2 + B.py * M.c1 -
                     B.pz * M.c15 + B.pw * M.c6;
    ctype const c3 = B.vx * M.c4 - B.vx * M.c11 + B.vy * M.c5 - B.vy * M.c12 +
                     B.vz * M.c6 - B.vz * M.c13 + B.mx * M.c1 + B.my * M.c2 +
                     B.mz * M.c3 - B.pw * M.c15;
    ctype const c4 = B.mx * M.c7 + B.my * M.c8 + B.mz * M.c9 + B.px * M.c4 +
                     B.px * M.c11 + B.py * M.c5 + B.py * M.c12 + B.pz * M.c6 +
                     B.pz * M.c13 + B.pw * M.c14;
    ctype const c5 = -B.vx * M.c14 - B.vy * M.c9 + B.vz * M.c8 - B.mx * M.c0 -
                     B.my * M.c6 + B.mz * M.c5 - B.px * M.c15 - B.py * M.c3 +
                     B.pz * M.c2 - B.pw * M.c11;
    ctype const c6 = B.vx * M.c9 - B.vy * M.c14 - B.vz * M.c7 + B.mx * M.c6 -
                     B.my * M.c0 - B.mz * M.c4 + B.px * M.c3 - B.py * M.c15 -
                     B.pz * M.c1 - B.pw * M.c12;
    ctype const c7 = -B.vx * M.c8 + B.vy * M.c7 - B.vz * M.c14 - B.mx * M.c5 +
                     B.my * M.c4 - B.mz * M.c0 - B.px * M.c2 + B.py * M.c1 -
                     B.pz * M.c15 - B.pw * M.c13;
    ctype const c8 = -B.mx * M.c14 - B.my * M.c9 + B.mz * M.c8 - B.px * M.c0 -
                     B.px * M.c10 - B.py * M.c6 - B.py * M.c13 + B.pz * M.c5 +
                     B.pz * M.c12 + B.pw * M.c7;
    ctype const c9 = B.mx * M.c9 - B.my * M.c14 - B.mz * M.c7 + B.px * M.c6 +
                     B.px * M.c13 - B.py * M.c0 - B.py * M.c10 - B.pz * M.c4 -
                     B.pz * M.c11 + B.pw * M.c8;
    ctype const c10 = -B.mx * M.c8 + B.my * M.c7 - B.mz * M.c14 - B.px * M.c5 -
                      B.px * M.c12 + B.py * M.c4 + B.py * M.c11 - B.pz * M.c0 -
                      B.pz * M.c10 + B.pw * M.c9;
    ctype const c11 = -B.vx * M.c0 + B.vx * M.c10 - B.vy * M.c6 + B.vy * M.c13 +
                      B.vz * M.c5 - B.vz * M.c12 - B.mx * M.c15 - B.my * M.c3 +
                      B.mz * M.c2 - B.pw * M.c1;
    ctype const c12 = B.vx * M.c6 - B.vx * M.c13 - B.vy * M.c0 + B.vy * M.c10 -
                      B.vz * M.c4 + B.vz * M.c11 + B.mx * M.c3 - B.my * M.c15 -
                      B.mz * M.c1 - B.pw * M.c2;
    ctype const c13 = -B.vx * M.c5 + B.vx * M.c12 + B.vy * M.c4 - B.vy * M.c11 -
                      B.vz * M.c0 + B.vz * M.c10 - B.mx * M.c2 + B.my * M.c1 -
                      B.mz * M.c15 - B.pw * M.c3;
    ctype const c14 = -B.vx * M.c7 - B.vy * M.c8 - B.vz * M.c9 - B.mx * M.c11 -
                      B.my * M.c12 - B.mz * M.c13 + B.px * M.c1 + B.py * M.c2 +
                      B.pz * M.c3 + B.pw * M.c0;
    ctype const c15 = B.vx * M.c7 + B.vy * M.c8 + B.vz * M.c9 + B.mx * M.c4 +
                      B.my * M.c5 + B.mz * M.c6 + B.px * M.c1 + B.py * M.c2 +
                      B.pz * M.c3 - B.pw * M.c10;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_e,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& M,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c11 * v.x + M.c12 * v.y + M.c13 * v.z + M.c14 * v.w + M.c15 * v.u;
    ctype const c1 = M.c2 * v.z - M.c3 * v.y + M.c4 * v.w + M.c11 * v.w + M.c15 * v.x;
    ctype const c2 = -M.c1 * v.z + M.c3 * v.x + M.c5 * v.w + M.c12 * v.w + M.c15 * v.y;
    ctype const c3 = M.c1 * v.y - M.c2 * v.x + M.c6 * v.w + M.c13 * v.w + M.c15 * v.z;
    ctype const c4 = M.c1 * v.u + M.c7 * v.w - M.c10 * v.x - M.c12 * v.z + M.c13 * v.y;
    ctype const c5 = M.c2 * v.u + M.c8 * v.w - M.c10 * v.y + M.c11 * v.z - M.c13 * v.x;
    ctype const c6 = M.c3 * v.u + M.c9 * v.w - M.c10 * v.z - M.c11 * v.y + M.c12 * v.x;
    ctype const c7 = M.c4 * v.u - M.c8 * v.z + M.c9 * v.y - M.c11 * v.u - M.c14 * v.x;
    ctype const c8 = M.c5 * v.u + M.c7 * v.z - M.c9 * v.x - M.c12 * v.u - M.c14 * v.y;
    ctype const c9 = M.c6 * v.u - M.c7 * v.y + M.c8 * v.x - M.c13 * v.u - M.c14 * v.z;
    ctype const c10 = M.c4 * v.x + M.c5 * v.y + M.c6 * v.z - M.c14 * v.w + M.c15 * v.u;
    ctype const c11 = -M.c0 * v.x + M.c1 * v.u - M.c5 * v.z + M.c6 * v.y - M.c7 * v.w;
    ctype const c12 = -M.c0 * v.y + M.c2 * v.u + M.c4 * v.z - M.c6 * v.x - M.c8 * v.w;
    ctype const c13 = -M.c0 * v.z + M.c3 * v.u - M.c4 * v.y + M.c5 * v.x - M.c9 * v.w;
    ctype const c14 = M.c0 * v.u + M.c7 * v.x + M.c8 * v.y + M.c9 * v.z - M.c10 * v.u;
    ctype const c15 = M.c0 * v.w - M.c1 * v.x - M.c2 * v.y - M.c3 * v.z + M.c10 * v.w;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(vec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v,
                                                   MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14 + v.u * M.c15;
    ctype const c1 = -v.x * M.c15 - v.y * M.c3 + v.z * M.c2 + v.w * M.c4 - v.w * M.c11;
    ctype const c2 = v.x * M.c3 - v.y * M.c15 - v.z * M.c1 + v.w * M.c5 - v.w * M.c12;
    ctype const c3 = -v.x * M.c2 + v.y * M.c1 - v.z * M.c15 + v.w * M.c6 - v.w * M.c13;
    ctype const c4 = -v.x * M.c10 - v.y * M.c13 + v.z * M.c12 + v.w * M.c7 + v.u * M.c1;
    ctype const c5 = v.x * M.c13 - v.y * M.c10 - v.z * M.c11 + v.w * M.c8 + v.u * M.c2;
    ctype const c6 = -v.x * M.c12 + v.y * M.c11 - v.z * M.c10 + v.w * M.c9 + v.u * M.c3;
    ctype const c7 = v.x * M.c14 + v.y * M.c9 - v.z * M.c8 + v.u * M.c4 + v.u * M.c11;
    ctype const c8 = -v.x * M.c9 + v.y * M.c14 + v.z * M.c7 + v.u * M.c5 + v.u * M.c12;
    ctype const c9 = v.x * M.c8 - v.y * M.c7 + v.z * M.c14 + v.u * M.c6 + v.u * M.c13;
    ctype const c10 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6 + v.w * M.c14 - v.u * M.c15;
    ctype const c11 = -v.x * M.c0 - v.y * M.c6 + v.z * M.c5 + v.w * M.c7 - v.u * M.c1;
    ctype const c12 = v.x * M.c6 - v.y * M.c0 - v.z * M.c4 + v.w * M.c8 - v.u * M.c2;
    ctype const c13 = -v.x * M.c5 + v.y * M.c4 - v.z * M.c0 + v.w * M.c9 - v.u * M.c3;
    ctype const c14 = -v.x * M.c7 - v.y * M.c8 - v.z * M.c9 + v.u * M.c0 + v.u * M.c10;
    ctype const c15 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3 + v.w * M.c0 - v.w * M.c10;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(mv_e,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_E<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c11 * ctype(s);
    ctype const c1 = M.c12 * ctype(s);
    ctype const c2 = M.c13 * ctype(s);
    ctype const c3 = -M.c15 * ctype(s);
    ctype const c4 = -M.c14 * ctype(s);
    ctype const c5 = -M.c4 * ctype(s);
    ctype const c6 = -M.c5 * ctype(s);
    ctype const c7 = -M.c6 * ctype(s);
    ctype const c8 = -M.c7 * ctype(s);
    ctype const c9 = -M.c8 * ctype(s);
    ctype const c10 = -M.c9 * ctype(s);
    ctype const c11 = -M.c1 * ctype(s);
    ctype const c12 = -M.c2 * ctype(s);
    ctype const c13 = -M.c3 * ctype(s);
    ctype const c14 = M.c10 * ctype(s);
    ctype const c15 = -M.c0 * ctype(s);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(s,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c11;
    ctype const c1 = ctype(s) * M.c12;
    ctype const c2 = ctype(s) * M.c13;
    ctype const c3 = -ctype(s) * M.c15;
    ctype const c4 = -ctype(s) * M.c14;
    ctype const c5 = -ctype(s) * M.c4;
    ctype const c6 = -ctype(s) * M.c5;
    ctype const c7 = -ctype(s) * M.c6;
    ctype const c8 = -ctype(s) * M.c7;
    ctype const c9 = -ctype(s) * M.c8;
    ctype const c10 = -ctype(s) * M.c9;
    ctype const c11 = -ctype(s) * M.c1;
    ctype const c12 = -ctype(s) * M.c2;
    ctype const c13 = -ctype(s) * M.c3;
    ctype const c14 = ctype(s) * M.c10;
    ctype const c15 = -ctype(s) * M.c0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_u,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& A,
                                                   MVec3dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 + A.c1 * B.c7 - A.c2 * B.c6 - A.c3 * B.c8 +
                     A.c4 * B.c11 + A.c5 * B.c14 + A.c6 * B.c2 - A.c7 * B.c1 +
                     A.c8 * B.c3 + A.c9 * B.c13 - A.c10 * B.c12 - A.c11 * B.c4 -
                     A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 + A.c15 * B.c0;
    ctype const c1 = -A.c0 * B.c7 + A.c1 * B.c15 + A.c2 * B.c5 - A.c3 * B.c9 +
                     A.c4 * B.c12 - A.c5 * B.c2 + A.c6 * B.c14 + A.c7 * B.c0 -
                     A.c8 * B.c13 + A.c9 * B.c3 + A.c10 * B.c11 + A.c11 * B.c10 -
                     A.c12 * B.c4 - A.c13 * B.c8 + A.c14 * B.c6 + A.c15 * B.c1;
    ctype const c2 = A.c0 * B.c6 - A.c1 * B.c5 + A.c2 * B.c15 - A.c3 * B.c10 +
                     A.c4 * B.c13 + A.c5 * B.c1 - A.c6 * B.c0 + A.c7 * B.c14 +
                     A.c8 * B.c12 - A.c9 * B.c11 + A.c10 * B.c3 - A.c11 * B.c9 +
                     A.c12 * B.c8 - A.c13 * B.c4 + A.c14 * B.c7 + A.c15 * B.c2;
    ctype const c3 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 +
                     A.c3 * B.c15 - A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 -
                     A.c11 * B.c0 - A.c11 * B.c5 - A.c12 * B.c1 - A.c12 * B.c6 -
                     A.c13 * B.c2 - A.c13 * B.c7 - A.c14 * B.c3 + A.c15 * B.c3;
    ctype const c4 = -A.c0 * B.c8 - A.c1 * B.c9 - A.c2 * B.c10 - A.c4 * B.c14 +
                     A.c4 * B.c15 - A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 +
                     A.c8 * B.c0 - A.c8 * B.c5 + A.c9 * B.c1 - A.c9 * B.c6 +
                     A.c10 * B.c2 - A.c10 * B.c7 + A.c14 * B.c4 + A.c15 * B.c4;
    ctype const c5 = A.c0 * B.c14 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c8 +
                     A.c4 * B.c11 + A.c5 * B.c15 + A.c6 * B.c7 - A.c7 * B.c6 +
                     A.c8 * B.c3 + A.c9 * B.c13 - A.c10 * B.c12 + A.c11 * B.c4 +
                     A.c12 * B.c10 - A.c13 * B.c9 + A.c14 * B.c0 + A.c15 * B.c5;
    ctype const c6 = -A.c0 * B.c2 + A.c1 * B.c14 + A.c2 * B.c0 + A.c3 * B.c9 +
                     A.c4 * B.c12 - A.c5 * B.c7 + A.c6 * B.c15 + A.c7 * B.c5 -
                     A.c8 * B.c13 + A.c9 * B.c3 + A.c10 * B.c11 - A.c11 * B.c10 +
                     A.c12 * B.c4 + A.c13 * B.c8 + A.c14 * B.c1 + A.c15 * B.c6;
    ctype const c7 = A.c0 * B.c1 - A.c1 * B.c0 + A.c2 * B.c14 + A.c3 * B.c10 +
                     A.c4 * B.c13 + A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c15 +
                     A.c8 * B.c12 - A.c9 * B.c11 + A.c10 * B.c3 + A.c11 * B.c9 -
                     A.c12 * B.c8 + A.c13 * B.c4 + A.c14 * B.c2 + A.c15 * B.c7;
    ctype const c8 = A.c0 * B.c4 + A.c1 * B.c10 - A.c2 * B.c9 - A.c4 * B.c0 +
                     A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c10 - A.c7 * B.c9 -
                     A.c8 * B.c14 + A.c8 * B.c15 - A.c9 * B.c2 + A.c9 * B.c7 +
                     A.c10 * B.c1 - A.c10 * B.c6 + A.c14 * B.c8 + A.c15 * B.c8;
    ctype const c9 = -A.c0 * B.c10 + A.c1 * B.c4 + A.c2 * B.c8 - A.c4 * B.c1 +
                     A.c4 * B.c6 - A.c5 * B.c10 + A.c6 * B.c4 + A.c7 * B.c8 +
                     A.c8 * B.c2 - A.c8 * B.c7 - A.c9 * B.c14 + A.c9 * B.c15 -
                     A.c10 * B.c0 + A.c10 * B.c5 + A.c14 * B.c9 + A.c15 * B.c9;
    ctype const c10 = A.c0 * B.c9 - A.c1 * B.c8 + A.c2 * B.c4 - A.c4 * B.c2 +
                      A.c4 * B.c7 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c4 -
                      A.c8 * B.c1 + A.c8 * B.c6 + A.c9 * B.c0 - A.c9 * B.c5 -
                      A.c10 * B.c14 + A.c10 * B.c15 + A.c14 * B.c10 + A.c15 * B.c10;
    ctype const c11 = -A.c0 * B.c3 - A.c1 * B.c13 + A.c2 * B.c12 + A.c3 * B.c0 +
                      A.c3 * B.c5 + A.c5 * B.c3 + A.c6 * B.c13 - A.c7 * B.c12 +
                      A.c11 * B.c14 + A.c11 * B.c15 + A.c12 * B.c2 + A.c12 * B.c7 -
                      A.c13 * B.c1 - A.c13 * B.c6 - A.c14 * B.c11 + A.c15 * B.c11;
    ctype const c12 = A.c0 * B.c13 - A.c1 * B.c3 - A.c2 * B.c11 + A.c3 * B.c1 +
                      A.c3 * B.c6 - A.c5 * B.c13 + A.c6 * B.c3 + A.c7 * B.c11 -
                      A.c11 * B.c2 - A.c11 * B.c7 + A.c12 * B.c14 + A.c12 * B.c15 +
                      A.c13 * B.c0 + A.c13 * B.c5 - A.c14 * B.c12 + A.c15 * B.c12;
    ctype const c13 = -A.c0 * B.c12 + A.c1 * B.c11 - A.c2 * B.c3 + A.c3 * B.c2 +
                      A.c3 * B.c7 + A.c5 * B.c12 - A.c6 * B.c11 + A.c7 * B.c3 +
                      A.c11 * B.c1 + A.c11 * B.c6 - A.c12 * B.c0 - A.c12 * B.c5 +
                      A.c13 * B.c14 + A.c13 * B.c15 - A.c14 * B.c13 + A.c15 * B.c13;
    ctype const c14 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c7 - A.c3 * B.c4 +
                      A.c4 * B.c3 - A.c5 * B.c0 - A.c6 * B.c1 - A.c7 * B.c2 -
                      A.c8 * B.c11 - A.c9 * B.c12 - A.c10 * B.c13 + A.c11 * B.c8 +
                      A.c12 * B.c9 + A.c13 * B.c10 + A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = -A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c4 +
                      A.c4 * B.c3 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7 -
                      A.c8 * B.c11 - A.c9 * B.c12 - A.c10 * B.c13 - A.c11 * B.c8 -
                      A.c12 * B.c9 - A.c13 * B.c10 + A.c14 * B.c14 + A.c15 * B.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& M,
                                                   PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    ctype const c4 = M.c4 * ctype(ps);
    ctype const c5 = M.c5 * ctype(ps);
    ctype const c6 = M.c6 * ctype(ps);
    ctype const c7 = M.c7 * ctype(ps);
    ctype const c8 = M.c8 * ctype(ps);
    ctype const c9 = M.c9 * ctype(ps);
    ctype const c10 = M.c10 * ctype(ps);
    ctype const c11 = M.c11 * ctype(ps);
    ctype const c12 = M.c12 * ctype(ps);
    ctype const c13 = M.c13 * ctype(ps);
    ctype const c14 = M.c14 * ctype(ps);
    ctype const c15 = M.c15 * ctype(ps);
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    ctype const c4 = ctype(ps) * M.c4;
    ctype const c5 = ctype(ps) * M.c5;
    ctype const c6 = ctype(ps) * M.c6;
    ctype const c7 = ctype(ps) * M.c7;
    ctype const c8 = ctype(ps) * M.c8;
    ctype const c9 = ctype(ps) * M.c9;
    ctype const c10 = ctype(ps) * M.c10;
    ctype const c11 = ctype(ps) * M.c11;
    ctype const c12 = ctype(ps) * M.c12;
    ctype const c13 = ctype(ps) * M.c13;
    ctype const c14 = ctype(ps) * M.c14;
    ctype const c15 = ctype(ps) * M.c15;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_u,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& M,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * Q.x + M.c1 * Q.y + M.c2 * Q.z + M.c3 * Q.w + M.c4 * Q.u;
    ctype const c1 = -M.c0 * Q.u - M.c3 * Q.x + M.c5 * Q.u - M.c12 * Q.z + M.c13 * Q.y;
    ctype const c2 = -M.c1 * Q.u - M.c3 * Q.y + M.c6 * Q.u + M.c11 * Q.z - M.c13 * Q.x;
    ctype const c3 = -M.c2 * Q.u - M.c3 * Q.z + M.c7 * Q.u - M.c11 * Q.y + M.c12 * Q.x;
    ctype const c4 = -M.c1 * Q.z + M.c2 * Q.y + M.c8 * Q.u + M.c11 * Q.w - M.c14 * Q.x;
    ctype const c5 = M.c0 * Q.z - M.c2 * Q.x + M.c9 * Q.u + M.c12 * Q.w - M.c14 * Q.y;
    ctype const c6 = -M.c0 * Q.y + M.c1 * Q.x + M.c10 * Q.u + M.c13 * Q.w - M.c14 * Q.z;
    ctype const c7 = M.c0 * Q.w + M.c4 * Q.x + M.c5 * Q.w + M.c9 * Q.z - M.c10 * Q.y;
    ctype const c8 = M.c1 * Q.w + M.c4 * Q.y + M.c6 * Q.w - M.c8 * Q.z + M.c10 * Q.x;
    ctype const c9 = M.c2 * Q.w + M.c4 * Q.z + M.c7 * Q.w + M.c8 * Q.y - M.c9 * Q.x;
    ctype const c10 = M.c3 * Q.w - M.c4 * Q.u - M.c5 * Q.x - M.c6 * Q.y - M.c7 * Q.z;
    ctype const c11 = M.c6 * Q.z - M.c7 * Q.y - M.c8 * Q.u + M.c11 * Q.w + M.c15 * Q.x;
    ctype const c12 = -M.c5 * Q.z + M.c7 * Q.x - M.c9 * Q.u + M.c12 * Q.w + M.c15 * Q.y;
    ctype const c13 = M.c5 * Q.y - M.c6 * Q.x - M.c10 * Q.u + M.c13 * Q.w + M.c15 * Q.z;
    ctype const c14 = -M.c8 * Q.x - M.c9 * Q.y - M.c10 * Q.z + M.c14 * Q.w + M.c15 * Q.w;
    ctype const c15 = M.c11 * Q.x + M.c12 * Q.y + M.c13 * Q.z - M.c14 * Q.u + M.c15 * Q.u;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(quadvec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * M.c0 + Q.y * M.c1 + Q.z * M.c2 + Q.w * M.c3 + Q.u * M.c4;
    ctype const c1 = Q.x * M.c3 + Q.y * M.c13 - Q.z * M.c12 + Q.u * M.c0 + Q.u * M.c5;
    ctype const c2 = -Q.x * M.c13 + Q.y * M.c3 + Q.z * M.c11 + Q.u * M.c1 + Q.u * M.c6;
    ctype const c3 = Q.x * M.c12 - Q.y * M.c11 + Q.z * M.c3 + Q.u * M.c2 + Q.u * M.c7;
    ctype const c4 = -Q.x * M.c14 - Q.y * M.c2 + Q.z * M.c1 + Q.w * M.c11 + Q.u * M.c8;
    ctype const c5 = Q.x * M.c2 - Q.y * M.c14 - Q.z * M.c0 + Q.w * M.c12 + Q.u * M.c9;
    ctype const c6 = -Q.x * M.c1 + Q.y * M.c0 - Q.z * M.c14 + Q.w * M.c13 + Q.u * M.c10;
    ctype const c7 = -Q.x * M.c4 - Q.y * M.c10 + Q.z * M.c9 - Q.w * M.c0 + Q.w * M.c5;
    ctype const c8 = Q.x * M.c10 - Q.y * M.c4 - Q.z * M.c8 - Q.w * M.c1 + Q.w * M.c6;
    ctype const c9 = -Q.x * M.c9 + Q.y * M.c8 - Q.z * M.c4 - Q.w * M.c2 + Q.w * M.c7;
    ctype const c10 = -Q.x * M.c5 - Q.y * M.c6 - Q.z * M.c7 - Q.w * M.c3 + Q.u * M.c4;
    ctype const c11 = Q.x * M.c15 + Q.y * M.c7 - Q.z * M.c6 - Q.w * M.c11 + Q.u * M.c8;
    ctype const c12 = -Q.x * M.c7 + Q.y * M.c15 + Q.z * M.c5 - Q.w * M.c12 + Q.u * M.c9;
    ctype const c13 = Q.x * M.c6 - Q.y * M.c5 + Q.z * M.c15 - Q.w * M.c13 + Q.u * M.c10;
    ctype const c14 = Q.x * M.c8 + Q.y * M.c9 + Q.z * M.c10 - Q.w * M.c14 + Q.w * M.c15;
    ctype const c15 =
        -Q.x * M.c11 - Q.y * M.c12 - Q.z * M.c13 + Q.u * M.c14 + Q.u * M.c15;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(mv_u,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& M,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * t.vz - M.c2 * t.vy - M.c3 * t.mx + M.c4 * t.px + M.c5 * t.pw +
                     M.c9 * t.pz - M.c10 * t.py - M.c12 * t.mz + M.c13 * t.my +
                     M.c14 * t.vx;
    ctype const c1 = -M.c0 * t.vz + M.c2 * t.vx - M.c3 * t.my + M.c4 * t.py +
                     M.c6 * t.pw - M.c8 * t.pz + M.c10 * t.px + M.c11 * t.mz -
                     M.c13 * t.mx + M.c14 * t.vy;
    ctype const c2 = M.c0 * t.vy - M.c1 * t.vx - M.c3 * t.mz + M.c4 * t.pz + M.c7 * t.pw +
                     M.c8 * t.py - M.c9 * t.px - M.c11 * t.my + M.c12 * t.mx +
                     M.c14 * t.vz;
    ctype const c3 = M.c0 * t.px + M.c1 * t.py + M.c2 * t.pz + M.c3 * t.pw - M.c5 * t.px -
                     M.c6 * t.py - M.c7 * t.pz - M.c11 * t.vx - M.c12 * t.vy -
                     M.c13 * t.vz;
    ctype const c4 = -M.c0 * t.mx - M.c1 * t.my - M.c2 * t.mz - M.c4 * t.pw -
                     M.c5 * t.mx - M.c6 * t.my - M.c7 * t.mz - M.c8 * t.vx - M.c9 * t.vy -
                     M.c10 * t.vz;
    ctype const c5 = M.c0 * t.pw + M.c3 * t.mx + M.c4 * t.px + M.c6 * t.vz - M.c7 * t.vy +
                     M.c9 * t.pz - M.c10 * t.py + M.c12 * t.mz - M.c13 * t.my +
                     M.c15 * t.vx;
    ctype const c6 = M.c1 * t.pw + M.c3 * t.my + M.c4 * t.py - M.c5 * t.vz + M.c7 * t.vx -
                     M.c8 * t.pz + M.c10 * t.px - M.c11 * t.mz + M.c13 * t.mx +
                     M.c15 * t.vy;
    ctype const c7 = M.c2 * t.pw + M.c3 * t.mz + M.c4 * t.pz + M.c5 * t.vy - M.c6 * t.vx +
                     M.c8 * t.py - M.c9 * t.px + M.c11 * t.my - M.c12 * t.mx +
                     M.c15 * t.vz;
    ctype const c8 = M.c1 * t.mz - M.c2 * t.my + M.c4 * t.vx + M.c6 * t.mz - M.c7 * t.my -
                     M.c8 * t.pw + M.c9 * t.vz - M.c10 * t.vy + M.c14 * t.mx +
                     M.c15 * t.mx;
    ctype const c9 = -M.c0 * t.mz + M.c2 * t.mx + M.c4 * t.vy - M.c5 * t.mz +
                     M.c7 * t.mx - M.c8 * t.vz - M.c9 * t.pw + M.c10 * t.vx +
                     M.c14 * t.my + M.c15 * t.my;
    ctype const c10 = M.c0 * t.my - M.c1 * t.mx + M.c4 * t.vz + M.c5 * t.my -
                      M.c6 * t.mx + M.c8 * t.vy - M.c9 * t.vx - M.c10 * t.pw +
                      M.c14 * t.mz + M.c15 * t.mz;
    ctype const c11 = -M.c1 * t.pz + M.c2 * t.py + M.c3 * t.vx + M.c6 * t.pz -
                      M.c7 * t.py + M.c11 * t.pw + M.c12 * t.vz - M.c13 * t.vy -
                      M.c14 * t.px + M.c15 * t.px;
    ctype const c12 = M.c0 * t.pz - M.c2 * t.px + M.c3 * t.vy - M.c5 * t.pz +
                      M.c7 * t.px - M.c11 * t.vz + M.c12 * t.pw + M.c13 * t.vx -
                      M.c14 * t.py + M.c15 * t.py;
    ctype const c13 = -M.c0 * t.py + M.c1 * t.px + M.c3 * t.vz + M.c5 * t.py -
                      M.c6 * t.px + M.c11 * t.vy - M.c12 * t.vx + M.c13 * t.pw -
                      M.c14 * t.pz + M.c15 * t.pz;
    ctype const c14 = -M.c0 * t.vx - M.c1 * t.vy - M.c2 * t.vz - M.c8 * t.px -
                      M.c9 * t.py - M.c10 * t.pz + M.c11 * t.mx + M.c12 * t.my +
                      M.c13 * t.mz + M.c15 * t.pw;
    ctype const c15 = -M.c5 * t.vx - M.c6 * t.vy - M.c7 * t.vz - M.c8 * t.px -
                      M.c9 * t.py - M.c10 * t.pz - M.c11 * t.mx - M.c12 * t.my -
                      M.c13 * t.mz + M.c14 * t.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(trivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * M.c14 + t.vy * M.c2 - t.vz * M.c1 + t.mx * M.c3 +
                     t.my * M.c13 - t.mz * M.c12 - t.px * M.c4 - t.py * M.c10 +
                     t.pz * M.c9 + t.pw * M.c5;
    ctype const c1 = -t.vx * M.c2 + t.vy * M.c14 + t.vz * M.c0 - t.mx * M.c13 +
                     t.my * M.c3 + t.mz * M.c11 + t.px * M.c10 - t.py * M.c4 -
                     t.pz * M.c8 + t.pw * M.c6;
    ctype const c2 = t.vx * M.c1 - t.vy * M.c0 + t.vz * M.c14 + t.mx * M.c12 -
                     t.my * M.c11 + t.mz * M.c3 - t.px * M.c9 + t.py * M.c8 -
                     t.pz * M.c4 + t.pw * M.c7;
    ctype const c3 = -t.vx * M.c11 - t.vy * M.c12 - t.vz * M.c13 - t.px * M.c0 -
                     t.px * M.c5 - t.py * M.c1 - t.py * M.c6 - t.pz * M.c2 - t.pz * M.c7 -
                     t.pw * M.c3;
    ctype const c4 = -t.vx * M.c8 - t.vy * M.c9 - t.vz * M.c10 + t.mx * M.c0 -
                     t.mx * M.c5 + t.my * M.c1 - t.my * M.c6 + t.mz * M.c2 - t.mz * M.c7 +
                     t.pw * M.c4;
    ctype const c5 = t.vx * M.c15 + t.vy * M.c7 - t.vz * M.c6 + t.mx * M.c3 +
                     t.my * M.c13 - t.mz * M.c12 + t.px * M.c4 + t.py * M.c10 -
                     t.pz * M.c9 + t.pw * M.c0;
    ctype const c6 = -t.vx * M.c7 + t.vy * M.c15 + t.vz * M.c5 - t.mx * M.c13 +
                     t.my * M.c3 + t.mz * M.c11 - t.px * M.c10 + t.py * M.c4 +
                     t.pz * M.c8 + t.pw * M.c1;
    ctype const c7 = t.vx * M.c6 - t.vy * M.c5 + t.vz * M.c15 + t.mx * M.c12 -
                     t.my * M.c11 + t.mz * M.c3 + t.px * M.c9 - t.py * M.c8 +
                     t.pz * M.c4 + t.pw * M.c2;
    ctype const c8 = t.vx * M.c4 + t.vy * M.c10 - t.vz * M.c9 - t.mx * M.c14 +
                     t.mx * M.c15 - t.my * M.c2 + t.my * M.c7 + t.mz * M.c1 -
                     t.mz * M.c6 + t.pw * M.c8;
    ctype const c9 = -t.vx * M.c10 + t.vy * M.c4 + t.vz * M.c8 + t.mx * M.c2 -
                     t.mx * M.c7 - t.my * M.c14 + t.my * M.c15 - t.mz * M.c0 +
                     t.mz * M.c5 + t.pw * M.c9;
    ctype const c10 = t.vx * M.c9 - t.vy * M.c8 + t.vz * M.c4 - t.mx * M.c1 +
                      t.mx * M.c6 + t.my * M.c0 - t.my * M.c5 - t.mz * M.c14 +
                      t.mz * M.c15 + t.pw * M.c10;
    ctype const c11 = t.vx * M.c3 + t.vy * M.c13 - t.vz * M.c12 + t.px * M.c14 +
                      t.px * M.c15 + t.py * M.c2 + t.py * M.c7 - t.pz * M.c1 -
                      t.pz * M.c6 - t.pw * M.c11;
    ctype const c12 = -t.vx * M.c13 + t.vy * M.c3 + t.vz * M.c11 - t.px * M.c2 -
                      t.px * M.c7 + t.py * M.c14 + t.py * M.c15 + t.pz * M.c0 +
                      t.pz * M.c5 - t.pw * M.c12;
    ctype const c13 = t.vx * M.c12 - t.vy * M.c11 + t.vz * M.c3 + t.px * M.c1 +
                      t.px * M.c6 - t.py * M.c0 - t.py * M.c5 + t.pz * M.c14 +
                      t.pz * M.c15 - t.pw * M.c13;
    ctype const c14 = -t.vx * M.c0 - t.vy * M.c1 - t.vz * M.c2 - t.mx * M.c11 -
                      t.my * M.c12 - t.mz * M.c13 + t.px * M.c8 + t.py * M.c9 +
                      t.pz * M.c10 + t.pw * M.c15;
    ctype const c15 = -t.vx * M.c5 - t.vy * M.c6 - t.vz * M.c7 - t.mx * M.c11 -
                      t.my * M.c12 - t.mz * M.c13 - t.px * M.c8 - t.py * M.c9 -
                      t.pz * M.c10 + t.pw * M.c14;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_u,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& M,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                     M.c9 * B.vy - M.c10 * B.vz - M.c11 * B.px - M.c12 * B.py -
                     M.c13 * B.pz - M.c14 * B.pw;
    ctype const c1 = -M.c1 * B.vz + M.c2 * B.vy + M.c3 * B.mx + M.c6 * B.vz -
                     M.c7 * B.vy - M.c11 * B.pw + M.c12 * B.mz - M.c13 * B.my -
                     M.c14 * B.vx + M.c15 * B.vx;
    ctype const c2 = M.c0 * B.vz - M.c2 * B.vx + M.c3 * B.my - M.c5 * B.vz + M.c7 * B.vx -
                     M.c11 * B.mz - M.c12 * B.pw + M.c13 * B.mx - M.c14 * B.vy +
                     M.c15 * B.vy;
    ctype const c3 = -M.c0 * B.vy + M.c1 * B.vx + M.c3 * B.mz + M.c5 * B.vy -
                     M.c6 * B.vx + M.c11 * B.my - M.c12 * B.mx - M.c13 * B.pw -
                     M.c14 * B.vz + M.c15 * B.vz;
    ctype const c4 = -M.c0 * B.pw + M.c3 * B.px + M.c4 * B.vx + M.c6 * B.mz -
                     M.c7 * B.my + M.c9 * B.vz - M.c10 * B.vy + M.c12 * B.pz -
                     M.c13 * B.py + M.c15 * B.mx;
    ctype const c5 = -M.c1 * B.pw + M.c3 * B.py + M.c4 * B.vy - M.c5 * B.mz +
                     M.c7 * B.mx - M.c8 * B.vz + M.c10 * B.vx - M.c11 * B.pz +
                     M.c13 * B.px + M.c15 * B.my;
    ctype const c6 = -M.c2 * B.pw + M.c3 * B.pz + M.c4 * B.vz + M.c5 * B.my -
                     M.c6 * B.mx + M.c8 * B.vy - M.c9 * B.vx + M.c11 * B.py -
                     M.c12 * B.px + M.c15 * B.mz;
    ctype const c7 = M.c1 * B.pz - M.c2 * B.py + M.c4 * B.mx + M.c6 * B.pz - M.c7 * B.py +
                     M.c8 * B.pw + M.c9 * B.mz - M.c10 * B.my + M.c14 * B.px +
                     M.c15 * B.px;
    ctype const c8 = -M.c0 * B.pz + M.c2 * B.px + M.c4 * B.my - M.c5 * B.pz +
                     M.c7 * B.px - M.c8 * B.mz + M.c9 * B.pw + M.c10 * B.mx +
                     M.c14 * B.py + M.c15 * B.py;
    ctype const c9 = M.c0 * B.py - M.c1 * B.px + M.c4 * B.mz + M.c5 * B.py - M.c6 * B.px +
                     M.c8 * B.my - M.c9 * B.mx + M.c10 * B.pw + M.c14 * B.pz +
                     M.c15 * B.pz;
    ctype const c10 = M.c0 * B.mx + M.c1 * B.my + M.c2 * B.mz + M.c8 * B.vx +
                      M.c9 * B.vy + M.c10 * B.vz - M.c11 * B.px - M.c12 * B.py -
                      M.c13 * B.pz + M.c15 * B.pw;
    ctype const c11 = -M.c1 * B.mz + M.c2 * B.my + M.c3 * B.px - M.c4 * B.vx +
                      M.c5 * B.pw - M.c9 * B.vz + M.c10 * B.vy + M.c12 * B.pz -
                      M.c13 * B.py - M.c14 * B.mx;
    ctype const c12 = M.c0 * B.mz - M.c2 * B.mx + M.c3 * B.py - M.c4 * B.vy +
                      M.c6 * B.pw + M.c8 * B.vz - M.c10 * B.vx - M.c11 * B.pz +
                      M.c13 * B.px - M.c14 * B.my;
    ctype const c13 = -M.c0 * B.my + M.c1 * B.mx + M.c3 * B.pz - M.c4 * B.vz +
                      M.c7 * B.pw - M.c8 * B.vy + M.c9 * B.vx + M.c11 * B.py -
                      M.c12 * B.px - M.c14 * B.mz;
    ctype const c14 = -M.c0 * B.px - M.c1 * B.py - M.c2 * B.pz + M.c4 * B.pw -
                      M.c5 * B.px - M.c6 * B.py - M.c7 * B.pz - M.c8 * B.mx -
                      M.c9 * B.my - M.c10 * B.mz;
    ctype const c15 = M.c0 * B.vx + M.c1 * B.vy + M.c2 * B.vz - M.c3 * B.pw -
                      M.c5 * B.vx - M.c6 * B.vy - M.c7 * B.vz - M.c11 * B.mx -
                      M.c12 * B.my - M.c13 * B.mz;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(bivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                     B.my * M.c6 - B.mz * M.c7 - B.px * M.c11 - B.py * M.c12 -
                     B.pz * M.c13 - B.pw * M.c14;
    ctype const c1 = B.vx * M.c14 + B.vx * M.c15 + B.vy * M.c2 + B.vy * M.c7 -
                     B.vz * M.c1 - B.vz * M.c6 + B.mx * M.c3 + B.my * M.c13 -
                     B.mz * M.c12 + B.pw * M.c11;
    ctype const c2 = -B.vx * M.c2 - B.vx * M.c7 + B.vy * M.c14 + B.vy * M.c15 +
                     B.vz * M.c0 + B.vz * M.c5 - B.mx * M.c13 + B.my * M.c3 +
                     B.mz * M.c11 + B.pw * M.c12;
    ctype const c3 = B.vx * M.c1 + B.vx * M.c6 - B.vy * M.c0 - B.vy * M.c5 +
                     B.vz * M.c14 + B.vz * M.c15 + B.mx * M.c12 - B.my * M.c11 +
                     B.mz * M.c3 + B.pw * M.c13;
    ctype const c4 = B.vx * M.c4 + B.vy * M.c10 - B.vz * M.c9 + B.mx * M.c15 +
                     B.my * M.c7 - B.mz * M.c6 + B.px * M.c3 + B.py * M.c13 -
                     B.pz * M.c12 - B.pw * M.c0;
    ctype const c5 = -B.vx * M.c10 + B.vy * M.c4 + B.vz * M.c8 - B.mx * M.c7 +
                     B.my * M.c15 + B.mz * M.c5 - B.px * M.c13 + B.py * M.c3 +
                     B.pz * M.c11 - B.pw * M.c1;
    ctype const c6 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c4 + B.mx * M.c6 - B.my * M.c5 +
                     B.mz * M.c15 + B.px * M.c12 - B.py * M.c11 + B.pz * M.c3 -
                     B.pw * M.c2;
    ctype const c7 = B.mx * M.c4 + B.my * M.c10 - B.mz * M.c9 - B.px * M.c14 +
                     B.px * M.c15 - B.py * M.c2 + B.py * M.c7 + B.pz * M.c1 -
                     B.pz * M.c6 - B.pw * M.c8;
    ctype const c8 = -B.mx * M.c10 + B.my * M.c4 + B.mz * M.c8 + B.px * M.c2 -
                     B.px * M.c7 - B.py * M.c14 + B.py * M.c15 - B.pz * M.c0 +
                     B.pz * M.c5 - B.pw * M.c9;
    ctype const c9 = B.mx * M.c9 - B.my * M.c8 + B.mz * M.c4 - B.px * M.c1 + B.px * M.c6 +
                     B.py * M.c0 - B.py * M.c5 - B.pz * M.c14 + B.pz * M.c15 -
                     B.pw * M.c10;
    ctype const c10 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 + B.mx * M.c0 +
                      B.my * M.c1 + B.mz * M.c2 + B.px * M.c11 + B.py * M.c12 +
                      B.pz * M.c13 + B.pw * M.c15;
    ctype const c11 = B.vx * M.c4 + B.vy * M.c10 - B.vz * M.c9 - B.mx * M.c14 -
                      B.my * M.c2 + B.mz * M.c1 - B.px * M.c3 - B.py * M.c13 +
                      B.pz * M.c12 + B.pw * M.c5;
    ctype const c12 = -B.vx * M.c10 + B.vy * M.c4 + B.vz * M.c8 + B.mx * M.c2 -
                      B.my * M.c14 - B.mz * M.c0 + B.px * M.c13 - B.py * M.c3 -
                      B.pz * M.c11 + B.pw * M.c6;
    ctype const c13 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c4 - B.mx * M.c1 +
                      B.my * M.c0 - B.mz * M.c14 - B.px * M.c12 + B.py * M.c11 -
                      B.pz * M.c3 + B.pw * M.c7;
    ctype const c14 = -B.mx * M.c8 - B.my * M.c9 - B.mz * M.c10 + B.px * M.c0 -
                      B.px * M.c5 + B.py * M.c1 - B.py * M.c6 + B.pz * M.c2 -
                      B.pz * M.c7 - B.pw * M.c4;
    ctype const c15 = -B.vx * M.c0 - B.vx * M.c5 - B.vy * M.c1 - B.vy * M.c6 -
                      B.vz * M.c2 - B.vz * M.c7 - B.mx * M.c11 - B.my * M.c12 -
                      B.mz * M.c13 + B.pw * M.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(mv_u,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& M,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c6 * v.z - M.c7 * v.y + M.c8 * v.w - M.c11 * v.u + M.c15 * v.x;
    ctype const c1 = -M.c5 * v.z + M.c7 * v.x + M.c9 * v.w - M.c12 * v.u + M.c15 * v.y;
    ctype const c2 = M.c5 * v.y - M.c6 * v.x + M.c10 * v.w - M.c13 * v.u + M.c15 * v.z;
    ctype const c3 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z - M.c14 * v.w + M.c15 * v.w;
    ctype const c4 = M.c8 * v.x + M.c9 * v.y + M.c10 * v.z + M.c14 * v.u + M.c15 * v.u;
    ctype const c5 = M.c1 * v.z - M.c2 * v.y + M.c8 * v.w + M.c11 * v.u + M.c14 * v.x;
    ctype const c6 = -M.c0 * v.z + M.c2 * v.x + M.c9 * v.w + M.c12 * v.u + M.c14 * v.y;
    ctype const c7 = M.c0 * v.y - M.c1 * v.x + M.c10 * v.w + M.c13 * v.u + M.c14 * v.z;
    ctype const c8 = M.c0 * v.u - M.c4 * v.x + M.c5 * v.u - M.c9 * v.z + M.c10 * v.y;
    ctype const c9 = M.c1 * v.u - M.c4 * v.y + M.c6 * v.u + M.c8 * v.z - M.c10 * v.x;
    ctype const c10 = M.c2 * v.u - M.c4 * v.z + M.c7 * v.u - M.c8 * v.y + M.c9 * v.x;
    ctype const c11 = -M.c0 * v.w + M.c3 * v.x + M.c5 * v.w + M.c12 * v.z - M.c13 * v.y;
    ctype const c12 = -M.c1 * v.w + M.c3 * v.y + M.c6 * v.w - M.c11 * v.z + M.c13 * v.x;
    ctype const c13 = -M.c2 * v.w + M.c3 * v.z + M.c7 * v.w + M.c11 * v.y - M.c12 * v.x;
    ctype const c14 = -M.c3 * v.u + M.c4 * v.w - M.c5 * v.x - M.c6 * v.y - M.c7 * v.z;
    ctype const c15 = -M.c0 * v.x - M.c1 * v.y - M.c2 * v.z + M.c3 * v.u + M.c4 * v.w;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(vec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v,
                                                   MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c15 + v.y * M.c7 - v.z * M.c6 - v.w * M.c8 + v.u * M.c11;
    ctype const c1 = -v.x * M.c7 + v.y * M.c15 + v.z * M.c5 - v.w * M.c9 + v.u * M.c12;
    ctype const c2 = v.x * M.c6 - v.y * M.c5 + v.z * M.c15 - v.w * M.c10 + v.u * M.c13;
    ctype const c3 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14 + v.w * M.c15;
    ctype const c4 = -v.x * M.c8 - v.y * M.c9 - v.z * M.c10 - v.u * M.c14 + v.u * M.c15;
    ctype const c5 = v.x * M.c14 + v.y * M.c2 - v.z * M.c1 + v.w * M.c8 + v.u * M.c11;
    ctype const c6 = -v.x * M.c2 + v.y * M.c14 + v.z * M.c0 + v.w * M.c9 + v.u * M.c12;
    ctype const c7 = v.x * M.c1 - v.y * M.c0 + v.z * M.c14 + v.w * M.c10 + v.u * M.c13;
    ctype const c8 = v.x * M.c4 + v.y * M.c10 - v.z * M.c9 - v.u * M.c0 + v.u * M.c5;
    ctype const c9 = -v.x * M.c10 + v.y * M.c4 + v.z * M.c8 - v.u * M.c1 + v.u * M.c6;
    ctype const c10 = v.x * M.c9 - v.y * M.c8 + v.z * M.c4 - v.u * M.c2 + v.u * M.c7;
    ctype const c11 = -v.x * M.c3 - v.y * M.c13 + v.z * M.c12 + v.w * M.c0 + v.w * M.c5;
    ctype const c12 = v.x * M.c13 - v.y * M.c3 - v.z * M.c11 + v.w * M.c1 + v.w * M.c6;
    ctype const c13 = -v.x * M.c12 + v.y * M.c11 - v.z * M.c3 + v.w * M.c2 + v.w * M.c7;
    ctype const c14 = -v.x * M.c5 - v.y * M.c6 - v.z * M.c7 - v.w * M.c4 + v.u * M.c3;
    ctype const c15 = -v.x * M.c0 - v.y * M.c1 - v.z * M.c2 + v.w * M.c4 + v.u * M.c3;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(mv_u,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(MVec3dc_U<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c15 * ctype(s);
    ctype const c1 = M.c11 * ctype(s);
    ctype const c2 = M.c12 * ctype(s);
    ctype const c3 = M.c13 * ctype(s);
    ctype const c4 = M.c5 * ctype(s);
    ctype const c5 = M.c6 * ctype(s);
    ctype const c6 = M.c7 * ctype(s);
    ctype const c7 = M.c8 * ctype(s);
    ctype const c8 = M.c9 * ctype(s);
    ctype const c9 = M.c10 * ctype(s);
    ctype const c10 = -M.c14 * ctype(s);
    ctype const c11 = -M.c0 * ctype(s);
    ctype const c12 = -M.c1 * ctype(s);
    ctype const c13 = -M.c2 * ctype(s);
    ctype const c14 = M.c4 * ctype(s);
    ctype const c15 = M.c3 * ctype(s);
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(s,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c15;
    ctype const c1 = ctype(s) * M.c11;
    ctype const c2 = ctype(s) * M.c12;
    ctype const c3 = ctype(s) * M.c13;
    ctype const c4 = ctype(s) * M.c5;
    ctype const c5 = ctype(s) * M.c6;
    ctype const c6 = ctype(s) * M.c7;
    ctype const c7 = ctype(s) * M.c8;
    ctype const c8 = ctype(s) * M.c9;
    ctype const c9 = ctype(s) * M.c10;
    ctype const c10 = -ctype(s) * M.c14;
    ctype const c11 = -ctype(s) * M.c0;
    ctype const c12 = -ctype(s) * M.c1;
    ctype const c13 = -ctype(s) * M.c2;
    ctype const c14 = ctype(s) * M.c4;
    ctype const c15 = ctype(s) * M.c3;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps1, PScalar3dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(ctype(ps1) * ctype(ps2));
}

// cga3dc rgpr :: rgpr(ps,quadvec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps,
                                                    QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * Q.x;
    ctype const c1 = ctype(ps) * Q.y;
    ctype const c2 = ctype(ps) * Q.z;
    ctype const c3 = ctype(ps) * Q.w;
    ctype const c4 = ctype(ps) * Q.u;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(quadvec,ps) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q,
                                                    PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * ctype(ps);
    ctype const c1 = Q.y * ctype(ps);
    ctype const c2 = Q.z * ctype(ps);
    ctype const c3 = Q.w * ctype(ps);
    ctype const c4 = Q.u * ctype(ps);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(ps,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * t.vx;
    ctype const c1 = ctype(ps) * t.vy;
    ctype const c2 = ctype(ps) * t.vz;
    ctype const c3 = ctype(ps) * t.mx;
    ctype const c4 = ctype(ps) * t.my;
    ctype const c5 = ctype(ps) * t.mz;
    ctype const c6 = ctype(ps) * t.px;
    ctype const c7 = ctype(ps) * t.py;
    ctype const c8 = ctype(ps) * t.pz;
    ctype const c9 = ctype(ps) * t.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(trivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t,
                                                   PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vx * ctype(ps);
    ctype const c1 = t.vy * ctype(ps);
    ctype const c2 = t.vz * ctype(ps);
    ctype const c3 = t.mx * ctype(ps);
    ctype const c4 = t.my * ctype(ps);
    ctype const c5 = t.mz * ctype(ps);
    ctype const c6 = t.px * ctype(ps);
    ctype const c7 = t.py * ctype(ps);
    ctype const c8 = t.pz * ctype(ps);
    ctype const c9 = t.pw * ctype(ps);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps, BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.vx;
    ctype const c1 = ctype(ps) * B.vy;
    ctype const c2 = ctype(ps) * B.vz;
    ctype const c3 = ctype(ps) * B.mx;
    ctype const c4 = ctype(ps) * B.my;
    ctype const c5 = ctype(ps) * B.mz;
    ctype const c6 = ctype(ps) * B.px;
    ctype const c7 = ctype(ps) * B.py;
    ctype const c8 = ctype(ps) * B.pz;
    ctype const c9 = ctype(ps) * B.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(ps);
    ctype const c1 = B.vy * ctype(ps);
    ctype const c2 = B.vz * ctype(ps);
    ctype const c3 = B.mx * ctype(ps);
    ctype const c4 = B.my * ctype(ps);
    ctype const c5 = B.mz * ctype(ps);
    ctype const c6 = B.px * ctype(ps);
    ctype const c7 = B.py * ctype(ps);
    ctype const c8 = B.pz * ctype(ps);
    ctype const c9 = B.pw * ctype(ps);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = ctype(ps) * v.w;
    ctype const c4 = ctype(ps) * v.u;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = v.w * ctype(ps);
    ctype const c4 = v.u * ctype(ps);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rgpr(PScalar3dc<T> ps, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(ps) * ctype(s));
}

// cga3dc rgpr :: rgpr(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dc<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dc<ctype>(ctype(s) * ctype(ps));
}

// cga3dc rgpr :: rgpr(quadvec,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q1,
                                                   QuadVec3dc<U> const& Q2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = -Q1.y * Q2.z + Q1.z * Q2.y;
    ctype const c6 = Q1.x * Q2.z - Q1.z * Q2.x;
    ctype const c7 = -Q1.x * Q2.y + Q1.y * Q2.x;
    ctype const c8 = Q1.x * Q2.w - Q1.w * Q2.x;
    ctype const c9 = Q1.y * Q2.w - Q1.w * Q2.y;
    ctype const c10 = Q1.z * Q2.w - Q1.w * Q2.z;
    ctype const c11 = -Q1.x * Q2.u + Q1.u * Q2.x;
    ctype const c12 = -Q1.y * Q2.u + Q1.u * Q2.y;
    ctype const c13 = -Q1.z * Q2.u + Q1.u * Q2.z;
    ctype const c14 = -Q1.w * Q2.u + Q1.u * Q2.w;
    ctype const c15 = Q1.x * Q2.x + Q1.y * Q2.y + Q1.z * Q2.z - Q1.w * Q2.u - Q1.u * Q2.w;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(quadvec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = Q.y * t.pz - Q.z * t.py + Q.u * t.vx;
    ctype const c2 = -Q.x * t.pz + Q.z * t.px + Q.u * t.vy;
    ctype const c3 = Q.x * t.py - Q.y * t.px + Q.u * t.vz;
    ctype const c4 = -Q.x * t.pw + Q.w * t.px + Q.u * t.mx;
    ctype const c5 = -Q.y * t.pw + Q.w * t.py + Q.u * t.my;
    ctype const c6 = -Q.z * t.pw + Q.w * t.pz + Q.u * t.mz;
    ctype const c7 = -Q.y * t.mz + Q.z * t.my + Q.w * t.vx;
    ctype const c8 = Q.x * t.mz - Q.z * t.mx + Q.w * t.vy;
    ctype const c9 = -Q.x * t.my + Q.y * t.mx + Q.w * t.vz;
    ctype const c10 = -Q.x * t.vx - Q.y * t.vy - Q.z * t.vz;
    ctype const c11 = Q.y * t.vz - Q.z * t.vy - Q.w * t.px + Q.u * t.mx;
    ctype const c12 = -Q.x * t.vz + Q.z * t.vx - Q.w * t.py + Q.u * t.my;
    ctype const c13 = Q.x * t.vy - Q.y * t.vx - Q.w * t.pz + Q.u * t.mz;
    ctype const c14 = Q.x * t.mx + Q.y * t.my + Q.z * t.mz - Q.w * t.pw;
    ctype const c15 = -Q.x * t.px - Q.y * t.py - Q.z * t.pz + Q.u * t.pw;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(trivec,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = t.vx * Q.u - t.py * Q.z + t.pz * Q.y;
    ctype const c2 = t.vy * Q.u + t.px * Q.z - t.pz * Q.x;
    ctype const c3 = t.vz * Q.u - t.px * Q.y + t.py * Q.x;
    ctype const c4 = t.mx * Q.u + t.px * Q.w - t.pw * Q.x;
    ctype const c5 = t.my * Q.u + t.py * Q.w - t.pw * Q.y;
    ctype const c6 = t.mz * Q.u + t.pz * Q.w - t.pw * Q.z;
    ctype const c7 = t.vx * Q.w + t.my * Q.z - t.mz * Q.y;
    ctype const c8 = t.vy * Q.w - t.mx * Q.z + t.mz * Q.x;
    ctype const c9 = t.vz * Q.w + t.mx * Q.y - t.my * Q.x;
    ctype const c10 = -t.vx * Q.x - t.vy * Q.y - t.vz * Q.z;
    ctype const c11 = t.vy * Q.z - t.vz * Q.y - t.mx * Q.u + t.px * Q.w;
    ctype const c12 = -t.vx * Q.z + t.vz * Q.x - t.my * Q.u + t.py * Q.w;
    ctype const c13 = t.vx * Q.y - t.vy * Q.x - t.mz * Q.u + t.pz * Q.w;
    ctype const c14 = -t.mx * Q.x - t.my * Q.y - t.mz * Q.z + t.pw * Q.w;
    ctype const c15 = t.px * Q.x + t.py * Q.y + t.pz * Q.z - t.pw * Q.u;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(quadvec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.y * B.mz - Q.z * B.my - Q.w * B.vx + Q.u * B.px;
    ctype const c1 = -Q.x * B.mz + Q.z * B.mx - Q.w * B.vy + Q.u * B.py;
    ctype const c2 = Q.x * B.my - Q.y * B.mx - Q.w * B.vz + Q.u * B.pz;
    ctype const c3 = Q.x * B.vx + Q.y * B.vy + Q.z * B.vz + Q.u * B.pw;
    ctype const c4 = -Q.x * B.px - Q.y * B.py - Q.z * B.pz - Q.w * B.pw;
    ctype const c5 = -Q.x * B.pw - Q.w * B.vx - Q.u * B.px;
    ctype const c6 = -Q.y * B.pw - Q.w * B.vy - Q.u * B.py;
    ctype const c7 = -Q.z * B.pw - Q.w * B.vz - Q.u * B.pz;
    ctype const c8 = Q.y * B.pz - Q.z * B.py - Q.w * B.mx;
    ctype const c9 = -Q.x * B.pz + Q.z * B.px - Q.w * B.my;
    ctype const c10 = Q.x * B.py - Q.y * B.px - Q.w * B.mz;
    ctype const c11 = -Q.y * B.vz + Q.z * B.vy - Q.u * B.mx;
    ctype const c12 = Q.x * B.vz - Q.z * B.vx - Q.u * B.my;
    ctype const c13 = -Q.x * B.vy + Q.y * B.vx - Q.u * B.mz;
    ctype const c14 = -Q.x * B.mx - Q.y * B.my - Q.z * B.mz;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(bivec,quadvec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * Q.w + B.my * Q.z - B.mz * Q.y - B.px * Q.u;
    ctype const c1 = B.vy * Q.w - B.mx * Q.z + B.mz * Q.x - B.py * Q.u;
    ctype const c2 = B.vz * Q.w + B.mx * Q.y - B.my * Q.x - B.pz * Q.u;
    ctype const c3 = -B.vx * Q.x - B.vy * Q.y - B.vz * Q.z - B.pw * Q.u;
    ctype const c4 = B.px * Q.x + B.py * Q.y + B.pz * Q.z + B.pw * Q.w;
    ctype const c5 = -B.vx * Q.w - B.px * Q.u - B.pw * Q.x;
    ctype const c6 = -B.vy * Q.w - B.py * Q.u - B.pw * Q.y;
    ctype const c7 = -B.vz * Q.w - B.pz * Q.u - B.pw * Q.z;
    ctype const c8 = -B.mx * Q.w - B.py * Q.z + B.pz * Q.y;
    ctype const c9 = -B.my * Q.w + B.px * Q.z - B.pz * Q.x;
    ctype const c10 = -B.mz * Q.w - B.px * Q.y + B.py * Q.x;
    ctype const c11 = B.vy * Q.z - B.vz * Q.y - B.mx * Q.u;
    ctype const c12 = -B.vx * Q.z + B.vz * Q.x - B.my * Q.u;
    ctype const c13 = B.vx * Q.y - B.vy * Q.x - B.mz * Q.u;
    ctype const c14 = -B.mx * Q.x - B.my * Q.y - B.mz * Q.z;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(quadvec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * v.x + Q.y * v.y + Q.z * v.z + Q.w * v.w + Q.u * v.u;
    ctype const c1 = Q.x * v.w + Q.u * v.x;
    ctype const c2 = Q.y * v.w + Q.u * v.y;
    ctype const c3 = Q.z * v.w + Q.u * v.z;
    ctype const c4 = -Q.y * v.z + Q.z * v.y;
    ctype const c5 = Q.x * v.z - Q.z * v.x;
    ctype const c6 = -Q.x * v.y + Q.y * v.x;
    ctype const c7 = -Q.x * v.u - Q.w * v.x;
    ctype const c8 = -Q.y * v.u - Q.w * v.y;
    ctype const c9 = -Q.z * v.u - Q.w * v.z;
    ctype const c10 = -Q.w * v.w + Q.u * v.u;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(vec,quadvec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v,
                                                   QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * Q.x + v.y * Q.y + v.z * Q.z + v.w * Q.w + v.u * Q.u;
    ctype const c1 = -v.x * Q.u - v.w * Q.x;
    ctype const c2 = -v.y * Q.u - v.w * Q.y;
    ctype const c3 = -v.z * Q.u - v.w * Q.z;
    ctype const c4 = -v.y * Q.z + v.z * Q.y;
    ctype const c5 = v.x * Q.z - v.z * Q.x;
    ctype const c6 = -v.x * Q.y + v.y * Q.x;
    ctype const c7 = v.x * Q.w + v.u * Q.x;
    ctype const c8 = v.y * Q.w + v.u * Q.y;
    ctype const c9 = v.z * Q.w + v.u * Q.z;
    ctype const c10 = v.w * Q.w - v.u * Q.u;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(quadvec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rgpr(QuadVec3dc<T> const& Q, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = Q.x * ctype(s);
    ctype const c1 = Q.y * ctype(s);
    ctype const c2 = Q.z * ctype(s);
    ctype const c3 = -Q.u * ctype(s);
    ctype const c4 = -Q.w * ctype(s);
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(s,quadvec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * Q.x;
    ctype const c1 = ctype(s) * Q.y;
    ctype const c2 = ctype(s) * Q.z;
    ctype const c3 = -ctype(s) * Q.u;
    ctype const c4 = -ctype(s) * Q.w;
    return Vec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(trivec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t1,
                                                   TriVec3dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.vx * t2.pw + t1.my * t2.pz - t1.mz * t2.py - t1.py * t2.mz +
                     t1.pz * t2.my + t1.pw * t2.vx;
    ctype const c1 = t1.vy * t2.pw - t1.mx * t2.pz + t1.mz * t2.px + t1.px * t2.mz -
                     t1.pz * t2.mx + t1.pw * t2.vy;
    ctype const c2 = t1.vz * t2.pw + t1.mx * t2.py - t1.my * t2.px - t1.px * t2.my +
                     t1.py * t2.mx + t1.pw * t2.vz;
    ctype const c3 = -t1.vx * t2.px - t1.vy * t2.py - t1.vz * t2.pz - t1.px * t2.vx -
                     t1.py * t2.vy - t1.pz * t2.vz;
    ctype const c4 = -t1.vx * t2.mx - t1.vy * t2.my - t1.vz * t2.mz - t1.mx * t2.vx -
                     t1.my * t2.vy - t1.mz * t2.vz;
    ctype const c5 = t1.vy * t2.vz - t1.vz * t2.vy + t1.my * t2.pz - t1.mz * t2.py +
                     t1.py * t2.mz - t1.pz * t2.my;
    ctype const c6 = -t1.vx * t2.vz + t1.vz * t2.vx - t1.mx * t2.pz + t1.mz * t2.px -
                     t1.px * t2.mz + t1.pz * t2.mx;
    ctype const c7 = t1.vx * t2.vy - t1.vy * t2.vx + t1.mx * t2.py - t1.my * t2.px +
                     t1.px * t2.my - t1.py * t2.mx;
    ctype const c8 = t1.vy * t2.mz - t1.vz * t2.my - t1.mx * t2.pw + t1.my * t2.vz -
                     t1.mz * t2.vy + t1.pw * t2.mx;
    ctype const c9 = -t1.vx * t2.mz + t1.vz * t2.mx - t1.mx * t2.vz - t1.my * t2.pw +
                     t1.mz * t2.vx + t1.pw * t2.my;
    ctype const c10 = t1.vx * t2.my - t1.vy * t2.mx + t1.mx * t2.vy - t1.my * t2.vx -
                      t1.mz * t2.pw + t1.pw * t2.mz;
    ctype const c11 = t1.vy * t2.pz - t1.vz * t2.py + t1.px * t2.pw + t1.py * t2.vz -
                      t1.pz * t2.vy - t1.pw * t2.px;
    ctype const c12 = -t1.vx * t2.pz + t1.vz * t2.px - t1.px * t2.vz + t1.py * t2.pw +
                      t1.pz * t2.vx - t1.pw * t2.py;
    ctype const c13 = t1.vx * t2.py - t1.vy * t2.px + t1.px * t2.vy - t1.py * t2.vx +
                      t1.pz * t2.pw - t1.pw * t2.pz;
    ctype const c14 = -t1.mx * t2.px - t1.my * t2.py - t1.mz * t2.pz + t1.px * t2.mx +
                      t1.py * t2.my + t1.pz * t2.mz;
    ctype const c15 = -t1.vx * t2.vx - t1.vy * t2.vy - t1.vz * t2.vz - t1.mx * t2.px -
                      t1.my * t2.py - t1.mz * t2.pz - t1.px * t2.mx - t1.py * t2.my -
                      t1.pz * t2.mz + t1.pw * t2.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(trivec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.vx * B.mx - t.vy * B.my - t.vz * B.mz - t.mx * B.vx -
                     t.my * B.vy - t.mz * B.vz - t.px * B.px - t.py * B.py - t.pz * B.pz -
                     t.pw * B.pw;
    ctype const c1 =
        t.vy * B.vz - t.vz * B.vy - t.px * B.pw + t.py * B.mz - t.pz * B.my - t.pw * B.vx;
    ctype const c2 = -t.vx * B.vz + t.vz * B.vx - t.px * B.mz - t.py * B.pw +
                     t.pz * B.mx - t.pw * B.vy;
    ctype const c3 =
        t.vx * B.vy - t.vy * B.vx + t.px * B.my - t.py * B.mx - t.pz * B.pw - t.pw * B.vz;
    ctype const c4 =
        t.vy * B.mz - t.vz * B.my + t.my * B.vz - t.mz * B.vy + t.py * B.pz - t.pz * B.py;
    ctype const c5 = -t.vx * B.mz + t.vz * B.mx - t.mx * B.vz + t.mz * B.vx -
                     t.px * B.pz + t.pz * B.px;
    ctype const c6 =
        t.vx * B.my - t.vy * B.mx + t.mx * B.vy - t.my * B.vx + t.px * B.py - t.py * B.px;
    ctype const c7 =
        t.vy * B.pz - t.vz * B.py + t.mx * B.pw + t.my * B.mz - t.mz * B.my + t.pw * B.px;
    ctype const c8 = -t.vx * B.pz + t.vz * B.px - t.mx * B.mz + t.my * B.pw +
                     t.mz * B.mx + t.pw * B.py;
    ctype const c9 =
        t.vx * B.py - t.vy * B.px + t.mx * B.my - t.my * B.mx + t.mz * B.pw + t.pw * B.pz;
    ctype const c10 =
        t.mx * B.vx + t.my * B.vy + t.mz * B.vz - t.px * B.px - t.py * B.py - t.pz * B.pz;
    ctype const c11 =
        t.vx * B.pw - t.my * B.vz + t.mz * B.vy + t.py * B.pz - t.pz * B.py - t.pw * B.mx;
    ctype const c12 =
        t.vy * B.pw + t.mx * B.vz - t.mz * B.vx - t.px * B.pz + t.pz * B.px - t.pw * B.my;
    ctype const c13 =
        t.vz * B.pw - t.mx * B.vy + t.my * B.vx + t.px * B.py - t.py * B.px - t.pw * B.mz;
    ctype const c14 = -t.vx * B.px - t.vy * B.py - t.vz * B.pz - t.mx * B.mx -
                      t.my * B.my - t.mz * B.mz;
    ctype const c15 = -t.vx * B.vx - t.vy * B.vy - t.vz * B.vz - t.px * B.mx -
                      t.py * B.my - t.pz * B.mz;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(bivec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * t.mx - B.vy * t.my - B.vz * t.mz - B.mx * t.vx -
                     B.my * t.vy - B.mz * t.vz - B.px * t.px - B.py * t.py - B.pz * t.pz -
                     B.pw * t.pw;
    ctype const c1 =
        B.vx * t.pw + B.vy * t.vz - B.vz * t.vy + B.my * t.pz - B.mz * t.py + B.pw * t.px;
    ctype const c2 = -B.vx * t.vz + B.vy * t.pw + B.vz * t.vx - B.mx * t.pz +
                     B.mz * t.px + B.pw * t.py;
    ctype const c3 =
        B.vx * t.vy - B.vy * t.vx + B.vz * t.pw + B.mx * t.py - B.my * t.px + B.pw * t.pz;
    ctype const c4 =
        B.vy * t.mz - B.vz * t.my + B.my * t.vz - B.mz * t.vy + B.py * t.pz - B.pz * t.py;
    ctype const c5 = -B.vx * t.mz + B.vz * t.mx - B.mx * t.vz + B.mz * t.vx -
                     B.px * t.pz + B.pz * t.px;
    ctype const c6 =
        B.vx * t.my - B.vy * t.mx + B.mx * t.vy - B.my * t.vx + B.px * t.py - B.py * t.px;
    ctype const c7 =
        B.my * t.mz - B.mz * t.my - B.px * t.pw + B.py * t.vz - B.pz * t.vy - B.pw * t.mx;
    ctype const c8 = -B.mx * t.mz + B.mz * t.mx - B.px * t.vz - B.py * t.pw +
                     B.pz * t.vx - B.pw * t.my;
    ctype const c9 =
        B.mx * t.my - B.my * t.mx + B.px * t.vy - B.py * t.vx - B.pz * t.pw - B.pw * t.mz;
    ctype const c10 = -B.vx * t.mx - B.vy * t.my - B.vz * t.mz + B.px * t.px +
                      B.py * t.py + B.pz * t.pz;
    ctype const c11 =
        B.vy * t.mz - B.vz * t.my - B.mx * t.pw - B.py * t.pz + B.pz * t.py + B.pw * t.vx;
    ctype const c12 = -B.vx * t.mz + B.vz * t.mx - B.my * t.pw + B.px * t.pz -
                      B.pz * t.px + B.pw * t.vy;
    ctype const c13 =
        B.vx * t.my - B.vy * t.mx - B.mz * t.pw - B.px * t.py + B.py * t.px + B.pw * t.vz;
    ctype const c14 = -B.mx * t.mx - B.my * t.my - B.mz * t.mz - B.px * t.vx -
                      B.py * t.vy - B.pz * t.vz;
    ctype const c15 = -B.vx * t.vx - B.vy * t.vy - B.vz * t.vz - B.mx * t.px -
                      B.my * t.py - B.mz * t.pz;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(trivec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.vy * v.z - t.vz * v.y + t.mx * v.w - t.px * v.u;
    ctype const c1 = -t.vx * v.z + t.vz * v.x + t.my * v.w - t.py * v.u;
    ctype const c2 = t.vx * v.y - t.vy * v.x + t.mz * v.w - t.pz * v.u;
    ctype const c3 = -t.px * v.x - t.py * v.y - t.pz * v.z - t.pw * v.w;
    ctype const c4 = t.mx * v.x + t.my * v.y + t.mz * v.z + t.pw * v.u;
    ctype const c5 = t.mx * v.w + t.px * v.u + t.pw * v.x;
    ctype const c6 = t.my * v.w + t.py * v.u + t.pw * v.y;
    ctype const c7 = t.mz * v.w + t.pz * v.u + t.pw * v.z;
    ctype const c8 = t.vx * v.u - t.my * v.z + t.mz * v.y;
    ctype const c9 = t.vy * v.u + t.mx * v.z - t.mz * v.x;
    ctype const c10 = t.vz * v.u - t.mx * v.y + t.my * v.x;
    ctype const c11 = t.vx * v.w + t.py * v.z - t.pz * v.y;
    ctype const c12 = t.vy * v.w - t.px * v.z + t.pz * v.x;
    ctype const c13 = t.vz * v.w + t.px * v.y - t.py * v.x;
    ctype const c14 = -t.vx * v.x - t.vy * v.y - t.vz * v.z;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(vec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v,
                                                   TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * t.vz - v.z * t.vy - v.w * t.mx + v.u * t.px;
    ctype const c1 = -v.x * t.vz + v.z * t.vx - v.w * t.my + v.u * t.py;
    ctype const c2 = v.x * t.vy - v.y * t.vx - v.w * t.mz + v.u * t.pz;
    ctype const c3 = v.x * t.px + v.y * t.py + v.z * t.pz + v.w * t.pw;
    ctype const c4 = -v.x * t.mx - v.y * t.my - v.z * t.mz - v.u * t.pw;
    ctype const c5 = v.x * t.pw + v.w * t.mx + v.u * t.px;
    ctype const c6 = v.y * t.pw + v.w * t.my + v.u * t.py;
    ctype const c7 = v.z * t.pw + v.w * t.mz + v.u * t.pz;
    ctype const c8 = v.y * t.mz - v.z * t.my + v.u * t.vx;
    ctype const c9 = -v.x * t.mz + v.z * t.mx + v.u * t.vy;
    ctype const c10 = v.x * t.my - v.y * t.mx + v.u * t.vz;
    ctype const c11 = -v.y * t.pz + v.z * t.py + v.w * t.vx;
    ctype const c12 = v.x * t.pz - v.z * t.px + v.w * t.vy;
    ctype const c13 = -v.x * t.py + v.y * t.px + v.w * t.vz;
    ctype const c14 = -v.x * t.vx - v.y * t.vy - v.z * t.vz;
    ctype const c15 = 0.0;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(trivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rgpr(TriVec3dc<T> const& t, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.px * ctype(s);
    ctype const c1 = t.py * ctype(s);
    ctype const c2 = t.pz * ctype(s);
    ctype const c3 = t.vx * ctype(s);
    ctype const c4 = t.vy * ctype(s);
    ctype const c5 = t.vz * ctype(s);
    ctype const c6 = t.mx * ctype(s);
    ctype const c7 = t.my * ctype(s);
    ctype const c8 = t.mz * ctype(s);
    ctype const c9 = -t.pw * ctype(s);
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(s,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.px;
    ctype const c1 = ctype(s) * t.py;
    ctype const c2 = ctype(s) * t.pz;
    ctype const c3 = ctype(s) * t.vx;
    ctype const c4 = ctype(s) * t.vy;
    ctype const c5 = ctype(s) * t.vz;
    ctype const c6 = ctype(s) * t.mx;
    ctype const c7 = ctype(s) * t.my;
    ctype const c8 = ctype(s) * t.mz;
    ctype const c9 = -ctype(s) * t.pw;
    return BiVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(bivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B1,
                                                   BiVec3dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vy * B2.pz - B1.vz * B2.py + B1.mx * B2.pw - B1.py * B2.vz +
                     B1.pz * B2.vy + B1.pw * B2.mx;
    ctype const c1 = -B1.vx * B2.pz + B1.vz * B2.px + B1.my * B2.pw + B1.px * B2.vz -
                     B1.pz * B2.vx + B1.pw * B2.my;
    ctype const c2 = B1.vx * B2.py - B1.vy * B2.px + B1.mz * B2.pw - B1.px * B2.vy +
                     B1.py * B2.vx + B1.pw * B2.mz;
    ctype const c3 = B1.vx * B2.mx + B1.vy * B2.my + B1.vz * B2.mz + B1.mx * B2.vx +
                     B1.my * B2.vy + B1.mz * B2.vz;
    ctype const c4 = B1.mx * B2.px + B1.my * B2.py + B1.mz * B2.pz + B1.px * B2.mx +
                     B1.py * B2.my + B1.pz * B2.mz;
    ctype const c5 = -B1.vy * B2.pz + B1.vz * B2.py - B1.my * B2.mz + B1.mz * B2.my -
                     B1.py * B2.vz + B1.pz * B2.vy;
    ctype const c6 = B1.vx * B2.pz - B1.vz * B2.px + B1.mx * B2.mz - B1.mz * B2.mx +
                     B1.px * B2.vz - B1.pz * B2.vx;
    ctype const c7 = -B1.vx * B2.py + B1.vy * B2.px - B1.mx * B2.my + B1.my * B2.mx -
                     B1.px * B2.vy + B1.py * B2.vx;
    ctype const c8 = -B1.my * B2.pz + B1.mz * B2.py - B1.px * B2.pw - B1.py * B2.mz +
                     B1.pz * B2.my + B1.pw * B2.px;
    ctype const c9 = B1.mx * B2.pz - B1.mz * B2.px + B1.px * B2.mz - B1.py * B2.pw -
                     B1.pz * B2.mx + B1.pw * B2.py;
    ctype const c10 = -B1.mx * B2.py + B1.my * B2.px - B1.px * B2.my + B1.py * B2.mx -
                      B1.pz * B2.pw + B1.pw * B2.pz;
    ctype const c11 = B1.vx * B2.pw - B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz +
                      B1.mz * B2.vy - B1.pw * B2.vx;
    ctype const c12 = B1.vx * B2.mz + B1.vy * B2.pw - B1.vz * B2.mx + B1.mx * B2.vz -
                      B1.mz * B2.vx - B1.pw * B2.vy;
    ctype const c13 = -B1.vx * B2.my + B1.vy * B2.mx + B1.vz * B2.pw - B1.mx * B2.vy +
                      B1.my * B2.vx - B1.pw * B2.vz;
    ctype const c14 = -B1.vx * B2.px - B1.vy * B2.py - B1.vz * B2.pz + B1.px * B2.vx +
                      B1.py * B2.vy + B1.pz * B2.vz;
    ctype const c15 = B1.vx * B2.px + B1.vy * B2.py + B1.vz * B2.pz + B1.mx * B2.mx +
                      B1.my * B2.my + B1.mz * B2.mz + B1.px * B2.vx + B1.py * B2.vy +
                      B1.pz * B2.vz - B1.pw * B2.pw;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(bivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B,
                                                   Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c2 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c3 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c4 = B.vx * v.u + B.px * v.w - B.pw * v.x;
    ctype const c5 = B.vy * v.u + B.py * v.w - B.pw * v.y;
    ctype const c6 = B.vz * v.u + B.pz * v.w - B.pw * v.z;
    ctype const c7 = B.mx * v.u - B.py * v.z + B.pz * v.y;
    ctype const c8 = B.my * v.u + B.px * v.z - B.pz * v.x;
    ctype const c9 = B.mz * v.u - B.px * v.y + B.py * v.x;
    ctype const c10 = B.mx * v.x + B.my * v.y + B.mz * v.z;
    ctype const c11 = B.vx * v.u - B.my * v.z + B.mz * v.y - B.px * v.w;
    ctype const c12 = B.vy * v.u + B.mx * v.z - B.mz * v.x - B.py * v.w;
    ctype const c13 = B.vz * v.u - B.mx * v.y + B.my * v.x - B.pz * v.w;
    ctype const c14 = B.px * v.x + B.py * v.y + B.pz * v.z - B.pw * v.u;
    ctype const c15 = -B.vx * v.x - B.vy * v.y - B.vz * v.z + B.pw * v.w;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(vec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v,
                                                   BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -v.y * B.vz + v.z * B.vy + v.w * B.mx;
    ctype const c2 = v.x * B.vz - v.z * B.vx + v.w * B.my;
    ctype const c3 = -v.x * B.vy + v.y * B.vx + v.w * B.mz;
    ctype const c4 = -v.x * B.pw + v.w * B.px + v.u * B.vx;
    ctype const c5 = -v.y * B.pw + v.w * B.py + v.u * B.vy;
    ctype const c6 = -v.z * B.pw + v.w * B.pz + v.u * B.vz;
    ctype const c7 = v.y * B.pz - v.z * B.py + v.u * B.mx;
    ctype const c8 = -v.x * B.pz + v.z * B.px + v.u * B.my;
    ctype const c9 = v.x * B.py - v.y * B.px + v.u * B.mz;
    ctype const c10 = v.x * B.mx + v.y * B.my + v.z * B.mz;
    ctype const c11 = -v.y * B.mz + v.z * B.my + v.w * B.px - v.u * B.vx;
    ctype const c12 = v.x * B.mz - v.z * B.mx + v.w * B.py - v.u * B.vy;
    ctype const c13 = -v.x * B.my + v.y * B.mx + v.w * B.pz - v.u * B.vz;
    ctype const c14 = -v.x * B.px - v.y * B.py - v.z * B.pz + v.u * B.pw;
    ctype const c15 = v.x * B.vx + v.y * B.vy + v.z * B.vz - v.w * B.pw;
    return MVec3dc_E<ctype>(Scalar3dc<ctype>(c0),
                            BiVec3dc<ctype>(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10),
                            QuadVec3dc<ctype>(c11, c12, c13, c14, c15));
}

// cga3dc rgpr :: rgpr(bivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rgpr(BiVec3dc<T> const& B, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.mx * ctype(s);
    ctype const c1 = -B.my * ctype(s);
    ctype const c2 = -B.mz * ctype(s);
    ctype const c3 = -B.px * ctype(s);
    ctype const c4 = -B.py * ctype(s);
    ctype const c5 = -B.pz * ctype(s);
    ctype const c6 = -B.vx * ctype(s);
    ctype const c7 = -B.vy * ctype(s);
    ctype const c8 = -B.vz * ctype(s);
    ctype const c9 = B.pw * ctype(s);
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(s,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * B.mx;
    ctype const c1 = -ctype(s) * B.my;
    ctype const c2 = -ctype(s) * B.mz;
    ctype const c3 = -ctype(s) * B.px;
    ctype const c4 = -ctype(s) * B.py;
    ctype const c5 = -ctype(s) * B.pz;
    ctype const c6 = -ctype(s) * B.vx;
    ctype const c7 = -ctype(s) * B.vy;
    ctype const c8 = -ctype(s) * B.vz;
    ctype const c9 = ctype(s) * B.pw;
    return TriVec3dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
}

// cga3dc rgpr :: rgpr(vec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v1,
                                                   Vec3dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = v1.y * v2.z - v1.z * v2.y;
    ctype const c6 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c7 = v1.x * v2.y - v1.y * v2.x;
    ctype const c8 = v1.x * v2.u - v1.u * v2.x;
    ctype const c9 = v1.y * v2.u - v1.u * v2.y;
    ctype const c10 = v1.z * v2.u - v1.u * v2.z;
    ctype const c11 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c12 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c13 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c14 = -v1.w * v2.u + v1.u * v2.w;
    ctype const c15 =
        -v1.x * v2.x - v1.y * v2.y - v1.z * v2.z + v1.w * v2.u + v1.u * v2.w;
    return MVec3dc_U<ctype>(Vec3dc<ctype>(c0, c1, c2, c3, c4),
                            TriVec3dc<ctype>(c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
                            PScalar3dc<ctype>(c15));
}

// cga3dc rgpr :: rgpr(vec,s) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rgpr(Vec3dc<T> const& v, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * ctype(s);
    ctype const c1 = -v.y * ctype(s);
    ctype const c2 = -v.z * ctype(s);
    ctype const c3 = v.u * ctype(s);
    ctype const c4 = v.w * ctype(s);
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(s,vec) -> quadvec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * v.x;
    ctype const c1 = -ctype(s) * v.y;
    ctype const c2 = -ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.u;
    ctype const c4 = ctype(s) * v.w;
    return QuadVec3dc<ctype>(c0, c1, c2, c3, c4);
}

// cga3dc rgpr :: rgpr(s,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dc<std::common_type_t<T, U>> rgpr(Scalar3dc<T> s1, Scalar3dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dc<ctype>(-ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors w.r.t. the
// geometric product: u * inv(u) = inv(u) * u = 1
//
// The cga3dc metric is non-degenerate (det G = -1), so every element with
// non-zero geometric square is invertible. NULL elements are not -- and in a
// conformal algebra these are geometrically meaningful, not exotic: every
// embedded point q satisfies q * q = 0, so inv() of a round point (and of the
// null basis vectors e4, e5 themselves) throws. Invert only elements known to
// be non-null (e.g. versors).
//
// for k-blades:  A^(-1) = rev(A) / <A rev(A)>_0
//   - scalar / vector: rev = +A, so A^(-1) = A / gr0(A A)
//   - higher grades and the mixed-grade multivectors: closed-form inverse of
//     Hitzer & Sangwine, "Multivector and multivector matrix inverses in real
//     Clifford algebras" (2016), section 8 (n = 5): with the three-step chain
//
//         w = x conj(x),   y = w rev(w) = x conj(x) gr_inv(x) rev(x),
//         z = y m14(y),
//
//     where the grade-sign map m14() negates grades 1 and 4 only, z is a real
//     scalar and  x^(-1) = conj(x) gr_inv(x) rev(x) m14(y) / z.  The formula
//     is an algebraic identity in the algebra's own geometric product, so it
//     holds for the non-orthogonal null-basis metric unchanged. (The paper's
//     scalarness proof leaves one isotropic special case to numerics, Remark
//     8.1; the u * inv(u) = 1 test gate covers it here.)
//
// The squared norm is SIGNED here (e.g. I_3dc * I_3dc = -1), so the null check
// uses std::abs(denominator): it rejects only (near-)null elements, never
// merely negative squares.
////////////////////////////////////////////////////////////////////////////////
// HINT: inv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Scalar3dc<T> inv(Scalar3dc<T> s)
{
    T sq_n = T(s) * T(s);
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "scalar");
    return Scalar3dc<T>(T(s) / sq_n);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec3dc<T> inv(Vec3dc<T> const& v)
{
    // v^(-1) = rev(v)/<v v>_0 = v/dot(v,v); throws for null vectors -- which
    // includes every embedded point and the null basis vectors e4, e5
    T sq_n = T(dot(v, v));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "vector");
    T inv = T(1.0) / sq_n;
    return Vec3dc<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv, v.u * inv);
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016 (n = 5 chain, see the section comment)
template <typename T>
    requires(numeric_type<T>)
inline BiVec3dc<T> inv(BiVec3dc<T> const& B)
{
    auto y = B * conj(B) * gr_inv(B) * rev(B);         // even: grades 0, 2, 4
    auto ymap = MVec3dc_E<T>(gr0(y), gr2(y), -gr4(y)); // m14: negate grades 1, 4
    T sq_n = T(gr0(y * ymap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "bivector");
    return gr2(conj(B) * gr_inv(B) * rev(B) * ymap) / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016 (n = 5 chain, see the section comment)
template <typename T>
    requires(numeric_type<T>)
inline TriVec3dc<T> inv(TriVec3dc<T> const& t)
{
    auto y = t * conj(t) * gr_inv(t) * rev(t);         // even: grades 0, 2, 4
    auto ymap = MVec3dc_E<T>(gr0(y), gr2(y), -gr4(y)); // m14: negate grades 1, 4
    T sq_n = T(gr0(y * ymap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "trivector");
    return gr3(conj(t) * gr_inv(t) * rev(t) * ymap) / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016 (n = 5 chain, see the section comment)
template <typename T>
    requires(numeric_type<T>)
inline QuadVec3dc<T> inv(QuadVec3dc<T> const& Q)
{
    auto y = Q * conj(Q) * gr_inv(Q) * rev(Q);         // even: grades 0, 2, 4
    auto ymap = MVec3dc_E<T>(gr0(y), gr2(y), -gr4(y)); // m14: negate grades 1, 4
    T sq_n = T(gr0(y * ymap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "quadvector");
    return gr4(conj(Q) * gr_inv(Q) * rev(Q) * ymap) / sq_n;
}

// the cga3dc metric is non-degenerate, so (unlike the pga pseudoscalars) the
// pseudoscalar is invertible: I_3dc * I_3dc = -1, hence inv(I_3dc) = -I_3dc
template <typename T>
    requires(numeric_type<T>)
inline PScalar3dc<T> inv(PScalar3dc<T> ps)
{
    T sq_n = T(ps * ps); // = -ps^2 (signed geometric square)
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "pseudoscalar");
    return PScalar3dc<T>(T(ps) / sq_n);
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016 (n = 5 chain, see the section comment)
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_E<T> inv(MVec3dc_E<T> const& E)
{
    auto y = E * conj(E) * gr_inv(E) * rev(E);         // even: grades 0, 2, 4
    auto ymap = MVec3dc_E<T>(gr0(y), gr2(y), -gr4(y)); // m14: negate grades 1, 4
    T sq_n = T(gr0(y * ymap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "even-grade multivector");
    return conj(E) * gr_inv(E) * rev(E) * ymap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016 (n = 5 chain, see the section comment)
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> inv(MVec3dc_U<T> const& U)
{
    auto y = U * conj(U) * gr_inv(U) * rev(U);         // even: grades 0, 2, 4
    auto ymap = MVec3dc_E<T>(gr0(y), gr2(y), -gr4(y)); // m14: negate grades 1, 4
    T sq_n = T(gr0(y * ymap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "odd-grade multivector");
    return conj(U) * gr_inv(U) * rev(U) * ymap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016 (n = 5 chain, see the section comment);
// a left and a right inverse are the same (see the paper)
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc<T> inv(MVec3dc<T> const& M)
{
    auto y = M * conj(M) * gr_inv(M) * rev(M);
    auto ymap = MVec3dc<T>(gr0(y), -gr1(y), gr2(y), gr3(y), -gr4(y), gr5(y)); // m14 map
    T sq_n = T(gr0(y * ymap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "multivector");
    return conj(M) * gr_inv(M) * rev(M) * ymap / sq_n;
}


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses w.r.t. the regressive geometric product:
//
//     u (v) rinv(u) = rinv(u) (v) u = I_3dc   (the pseudoscalar, identity of rgpr)
//
// defined as the geometric inverse carried through the complement map
// (odd-dimensional algebra, so the complement is singular and involutive):
//
//     rinv(u) = cmpl(inv(cmpl(u)))
//
// Unlike pga (where the scalar has no rinv, mirroring the pseudoscalar having
// no inv), the non-degenerate cga3dc metric admits BOTH inverses for every
// non-null element, so rinv covers all types -- e.g. rinv(1) = -1 (a scalar:
// rgpr(1,-1) = I_3dc), the complement-dual of inv(I_3dc) = -I_3dc. Null
// elements throw, as for inv().
////////////////////////////////////////////////////////////////////////////////
// HINT: rinv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Scalar3dc<T> rinv(Scalar3dc<T> s)
{
    return cmpl(inv(cmpl(s)));
}

template <typename T>
    requires(numeric_type<T>)
inline Vec3dc<T> rinv(Vec3dc<T> const& v)
{
    return cmpl(inv(cmpl(v)));
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec3dc<T> rinv(BiVec3dc<T> const& B)
{
    return cmpl(inv(cmpl(B)));
}

template <typename T>
    requires(numeric_type<T>)
inline TriVec3dc<T> rinv(TriVec3dc<T> const& t)
{
    return cmpl(inv(cmpl(t)));
}

template <typename T>
    requires(numeric_type<T>)
inline QuadVec3dc<T> rinv(QuadVec3dc<T> const& Q)
{
    return cmpl(inv(cmpl(Q)));
}

template <typename T>
    requires(numeric_type<T>)
inline PScalar3dc<T> rinv(PScalar3dc<T> ps)
{
    return cmpl(inv(cmpl(ps)));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_E<T> rinv(MVec3dc_E<T> const& E)
{
    return cmpl(inv(cmpl(E)));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> rinv(MVec3dc_U<T> const& U)
{
    return cmpl(inv(cmpl(U)));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dc<T> rinv(MVec3dc<T> const& M)
{
    return cmpl(inv(cmpl(M)));
}


////////////////////////////////////////////////////////////////////////////////
// left and right expansions
//
// Expansions are the wdg-based duals of the contractions (which use rwdg):
//
//     l_expand3dc(a, b) = wdg( dual(a), b )
//     r_expand3dc(a, b) = wdg( a, dual(b) )
//
// (the odd-dimensional cga3dc has a single dual, so the left/right split is
// carried by the operand position alone). Implemented as generic wrappers over
// the existing dual and wedge operations, so every grade-pair combination is
// covered automatically. The requires-clause admits exactly the operand types
// for which the composition is well-formed (i.e. the cga3dc blade and
// multivector types).
////////////////////////////////////////////////////////////////////////////////

template <typename A, typename B>
    requires requires(A const& a, B const& b) { wdg(dual(a), b); }
constexpr auto l_expand3dc(A const& a, B const& b)
{
    return wdg(dual(a), b);
}

template <typename A, typename B>
    requires requires(A const& a, B const& b) { wdg(a, dual(b)); }
constexpr auto r_expand3dc(A const& a, B const& b)
{
    return wdg(a, dual(b));
}

} // namespace hd::ga::cga
