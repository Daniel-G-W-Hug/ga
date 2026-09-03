#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_cga2dc_ops_basics.hpp"

#include "detail/ga_error_handling.hpp" // inv()/rinv() null checks


namespace hd::ga::cga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides cga2dc product operations:
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
// - twdg1()                 -> transwedge product (k=1)
// - rtwdg1()                -> regressive transwedge product (k=1)
//
// - operator*()             -> geometric product (= gpr)
// - rgpr()                  -> regressive geometric product
//
// Expansions are the wdg-based duals of the corresponding contractions.
//
// - l_expand2dc()           -> left expansion
// - r_expand2dc()           -> right expansion
//
// - inv()                   -> inversion operation (w.r.t. geometric product)
// - rinv()                  -> inversion operation (w.r.t. regressive geometric
//                              product; full type coverage, see the inverses
//                              section: the metric is non-degenerate, so neither
//                              of the pga exclusions applies)
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot products (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here the cga2dc vector metric is the non-orthogonal null-pair metric
//
//     e1.e1 = e2.e2 = +1,  e3.e3 = e4.e4 = 0,  e3.e4 = e4.e3 = -1
//
// (det g = -1, non-degenerate). Via the exomorphism G(a ^ b) = G(a) ^ G(b) the
// extended metric G is defined for all grades; for cga2dc it is a SIGNED
// PERMUTATION that pairs each projective (w-bearing) blade with a round
// (u-bearing) partner blade, so the dot rules are off-diagonal but stay
// single-term (e.g. dot(e3,e4) = -1).
////////////////////////////////////////////////////////////////////////////////

// cga2dc dot :: dot(mv,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(MVec2dc<T> const& A,
                                                  MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 -
                            A.c4 * B.c3 + A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c7 +
                            A.c8 * B.c5 + A.c9 * B.c6 - A.c10 * B.c10 - A.c11 * B.c11 -
                            A.c12 * B.c12 + A.c13 * B.c14 + A.c14 * B.c13 -
                            A.c15 * B.c15);
}

// cga2dc dot :: dot(mv_e,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(MVec2dc_E<T> const& A,
                                                  MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(A.c0 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c3 +
                            A.c4 * B.c1 + A.c5 * B.c2 - A.c6 * B.c6 - A.c7 * B.c7);
}

// cga2dc dot :: dot(mv_u,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(MVec2dc_U<T> const& A,
                                                  MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(A.c0 * B.c0 + A.c1 * B.c1 - A.c2 * B.c3 - A.c3 * B.c2 -
                            A.c4 * B.c4 - A.c5 * B.c5 + A.c6 * B.c7 + A.c7 * B.c6);
}

// cga2dc dot :: dot(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(PScalar2dc<T> ps1, PScalar2dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga2dc dot :: dot(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(TriVec2dc<T> const& t1,
                                                  TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-t1.x * t2.x - t1.y * t2.y + t1.z * t2.w + t1.w * t2.z);
}

// cga2dc dot :: dot(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(BiVec2dc<T> const& B1,
                                                  BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(B1.vx * B2.mx + B1.vy * B2.my + B1.vz * B2.vz +
                            B1.mx * B2.vx + B1.my * B2.vy - B1.mz * B2.mz);
}

// cga2dc dot :: dot(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(Vec2dc<T> const& v1,
                                                  Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(v1.x * v2.x + v1.y * v2.y - v1.z * v2.w - v1.w * v2.z);
}

// cga2dc dot :: dot(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> dot(Scalar2dc<T> s1, Scalar2dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s1) * ctype(s2));
}
////////////////////////////////////////////////////////////////////////////////
// regressive dot products
// (= inner product of the complements, defined for equal grades exclusively)
////////////////////////////////////////////////////////////////////////////////

// cga2dc rdot :: rdot(mv,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(MVec2dc<T> const& A,
                                                    MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c4 +
                             A.c4 * B.c3 - A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c7 -
                             A.c8 * B.c5 - A.c9 * B.c6 + A.c10 * B.c10 + A.c11 * B.c11 +
                             A.c12 * B.c12 - A.c13 * B.c14 - A.c14 * B.c13 +
                             A.c15 * B.c15);
}

// cga2dc rdot :: rdot(mv_e,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(MVec2dc_E<T> const& A,
                                                    MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-A.c0 * B.c0 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c3 -
                             A.c4 * B.c1 - A.c5 * B.c2 + A.c6 * B.c6 + A.c7 * B.c7);
}

// cga2dc rdot :: rdot(mv_u,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(MVec2dc_U<T> const& A,
                                                    MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-A.c0 * B.c0 - A.c1 * B.c1 + A.c2 * B.c3 + A.c3 * B.c2 +
                             A.c4 * B.c4 + A.c5 * B.c5 - A.c6 * B.c7 - A.c7 * B.c6);
}

// cga2dc rdot :: rdot(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(PScalar2dc<T> ps1, PScalar2dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps1) * ctype(ps2));
}

// cga2dc rdot :: rdot(trivec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(TriVec2dc<T> const& t1,
                                                    TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(t1.x * t2.x + t1.y * t2.y - t1.z * t2.w - t1.w * t2.z);
}

// cga2dc rdot :: rdot(bivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(BiVec2dc<T> const& B1,
                                                    BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.vz -
                             B1.mx * B2.vx - B1.my * B2.vy + B1.mz * B2.mz);
}

// cga2dc rdot :: rdot(vec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(Vec2dc<T> const& v1,
                                                    Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-v1.x * v2.x - v1.y * v2.y + v1.z * v2.w + v1.w * v2.z);
}

// cga2dc rdot :: rdot(s,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rdot(Scalar2dc<T> s1, Scalar2dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-ctype(s1) * ctype(s2));
}
////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product)
////////////////////////////////////////////////////////////////////////////////

// cga2dc wdg :: wdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& A, MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c3 * B.c1 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 - A.c2 * B.c3 + A.c3 * B.c2 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c2 - A.c2 * B.c1 + A.c7 * B.c0;
    ctype const c8 = A.c0 * B.c8 + A.c1 * B.c4 - A.c4 * B.c1 + A.c8 * B.c0;
    ctype const c9 = A.c0 * B.c9 + A.c2 * B.c4 - A.c4 * B.c2 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c10 + A.c3 * B.c4 - A.c4 * B.c3 + A.c10 * B.c0;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c10 + A.c3 * B.c8 + A.c4 * B.c5 +
                      A.c5 * B.c4 + A.c8 * B.c3 - A.c10 * B.c1 + A.c11 * B.c0;
    ctype const c12 = A.c0 * B.c12 - A.c2 * B.c10 + A.c3 * B.c9 + A.c4 * B.c6 +
                      A.c6 * B.c4 + A.c9 * B.c3 - A.c10 * B.c2 + A.c12 * B.c0;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c9 - A.c2 * B.c8 + A.c4 * B.c7 +
                      A.c7 * B.c4 - A.c8 * B.c2 + A.c9 * B.c1 + A.c13 * B.c0;
    ctype const c14 = A.c0 * B.c14 + A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c7 -
                      A.c5 * B.c2 + A.c6 * B.c1 - A.c7 * B.c3 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c12 + A.c2 * B.c11 + A.c3 * B.c13 +
                      A.c4 * B.c14 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c10 -
                      A.c8 * B.c6 + A.c9 * B.c5 + A.c10 * B.c7 - A.c11 * B.c2 +
                      A.c12 * B.c1 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& A,
                                                MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c1 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c2 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c3 + A.c7 * B.c0;
    ctype const c8 = A.c0 * B.c4 + A.c8 * B.c0;
    ctype const c9 = A.c0 * B.c5 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c6 + A.c10 * B.c0;
    ctype const c11 = -A.c1 * B.c6 + A.c3 * B.c4 + A.c4 * B.c1 + A.c11 * B.c0;
    ctype const c12 = -A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c2 + A.c12 * B.c0;
    ctype const c13 = A.c1 * B.c5 - A.c2 * B.c4 + A.c4 * B.c3 + A.c13 * B.c0;
    ctype const c14 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c3 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c7 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c6 -
                      A.c8 * B.c2 + A.c9 * B.c1 + A.c10 * B.c3 + A.c15 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& A,
                                                MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c2 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c3 * B.c0;
    ctype const c8 = A.c0 * B.c8 + A.c4 * B.c0;
    ctype const c9 = A.c0 * B.c9 + A.c5 * B.c0;
    ctype const c10 = A.c0 * B.c10 + A.c6 * B.c0;
    ctype const c11 = A.c0 * B.c11 + A.c1 * B.c4 + A.c4 * B.c3 - A.c6 * B.c1;
    ctype const c12 = A.c0 * B.c12 + A.c2 * B.c4 + A.c5 * B.c3 - A.c6 * B.c2;
    ctype const c13 = A.c0 * B.c13 + A.c3 * B.c4 - A.c4 * B.c2 + A.c5 * B.c1;
    ctype const c14 = A.c0 * B.c14 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c3;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c10 -
                      A.c4 * B.c6 + A.c5 * B.c5 + A.c6 * B.c7 + A.c7 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& A,
                                                MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c0 * B.c1;
    ctype const c3 = A.c0 * B.c2;
    ctype const c4 = A.c0 * B.c3;
    ctype const c5 = -A.c1 * B.c2 + A.c3 * B.c0;
    ctype const c6 = -A.c2 * B.c2 + A.c3 * B.c1;
    ctype const c7 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype const c8 = A.c1 * B.c3 - A.c4 * B.c0;
    ctype const c9 = A.c2 * B.c3 - A.c4 * B.c1;
    ctype const c10 = A.c3 * B.c3 - A.c4 * B.c2;
    ctype const c11 = A.c0 * B.c4 + A.c5 * B.c3 + A.c8 * B.c2 - A.c10 * B.c0;
    ctype const c12 = A.c0 * B.c5 + A.c6 * B.c3 + A.c9 * B.c2 - A.c10 * B.c1;
    ctype const c13 = A.c0 * B.c6 + A.c7 * B.c3 - A.c8 * B.c1 + A.c9 * B.c0;
    ctype const c14 = A.c0 * B.c7 - A.c5 * B.c1 + A.c6 * B.c0 - A.c7 * B.c2;
    ctype const c15 = -A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c6 + A.c4 * B.c7 -
                      A.c11 * B.c1 + A.c12 * B.c0 - A.c13 * B.c2 - A.c14 * B.c3;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc_U<T> const& A,
                                                MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c1 * B.c0;
    ctype const c3 = A.c2 * B.c0;
    ctype const c4 = A.c3 * B.c0;
    ctype const c5 = -A.c0 * B.c3 + A.c2 * B.c1;
    ctype const c6 = -A.c1 * B.c3 + A.c2 * B.c2;
    ctype const c7 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c8 = A.c0 * B.c4 - A.c3 * B.c1;
    ctype const c9 = A.c1 * B.c4 - A.c3 * B.c2;
    ctype const c10 = A.c2 * B.c4 - A.c3 * B.c3;
    ctype const c11 = -A.c0 * B.c10 + A.c2 * B.c8 + A.c3 * B.c5 + A.c4 * B.c0;
    ctype const c12 = -A.c1 * B.c10 + A.c2 * B.c9 + A.c3 * B.c6 + A.c5 * B.c0;
    ctype const c13 = A.c0 * B.c9 - A.c1 * B.c8 + A.c3 * B.c7 + A.c6 * B.c0;
    ctype const c14 = A.c0 * B.c6 - A.c1 * B.c5 - A.c2 * B.c7 + A.c7 * B.c0;
    ctype const c15 = -A.c0 * B.c12 + A.c1 * B.c11 + A.c2 * B.c13 + A.c3 * B.c14 -
                      A.c4 * B.c2 + A.c5 * B.c1 - A.c6 * B.c3 - A.c7 * B.c4;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& M, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(M.c0 * ctype(ps));
}

// cga2dc wdg :: wdg(ps,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(PScalar2dc<T> ps, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps) * M.c0);
}

// cga2dc wdg :: wdg(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& M,
                                                TriVec2dc<U> const& t)
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
    ctype const c11 = M.c0 * t.x;
    ctype const c12 = M.c0 * t.y;
    ctype const c13 = M.c0 * t.z;
    ctype const c14 = M.c0 * t.w;
    ctype const c15 = -M.c1 * t.y + M.c2 * t.x + M.c3 * t.z + M.c4 * t.w;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(TriVec2dc<T> const& t,
                                                MVec2dc<U> const& M)
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
    ctype const c11 = t.x * M.c0;
    ctype const c12 = t.y * M.c0;
    ctype const c13 = t.z * M.c0;
    ctype const c14 = t.w * M.c0;
    ctype const c15 = -t.x * M.c2 + t.y * M.c1 - t.z * M.c3 - t.w * M.c4;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& M, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = M.c0 * B.vx;
    ctype const c6 = M.c0 * B.vy;
    ctype const c7 = M.c0 * B.vz;
    ctype const c8 = M.c0 * B.mx;
    ctype const c9 = M.c0 * B.my;
    ctype const c10 = M.c0 * B.mz;
    ctype const c11 = -M.c1 * B.mz + M.c3 * B.mx + M.c4 * B.vx;
    ctype const c12 = -M.c2 * B.mz + M.c3 * B.my + M.c4 * B.vy;
    ctype const c13 = M.c1 * B.my - M.c2 * B.mx + M.c4 * B.vz;
    ctype const c14 = M.c1 * B.vy - M.c2 * B.vx - M.c3 * B.vz;
    ctype const c15 = M.c5 * B.my - M.c6 * B.mx + M.c7 * B.mz - M.c8 * B.vy +
                      M.c9 * B.vx + M.c10 * B.vz;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(BiVec2dc<T> const& B, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = B.vx * M.c0;
    ctype const c6 = B.vy * M.c0;
    ctype const c7 = B.vz * M.c0;
    ctype const c8 = B.mx * M.c0;
    ctype const c9 = B.my * M.c0;
    ctype const c10 = B.mz * M.c0;
    ctype const c11 = B.vx * M.c4 + B.mx * M.c3 - B.mz * M.c1;
    ctype const c12 = B.vy * M.c4 + B.my * M.c3 - B.mz * M.c2;
    ctype const c13 = B.vz * M.c4 - B.mx * M.c2 + B.my * M.c1;
    ctype const c14 = -B.vx * M.c2 + B.vy * M.c1 - B.vz * M.c3;
    ctype const c15 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c10 - B.mx * M.c6 +
                      B.my * M.c5 + B.mz * M.c7;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& M, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = M.c0 * v.w;
    ctype const c5 = -M.c1 * v.z + M.c3 * v.x;
    ctype const c6 = -M.c2 * v.z + M.c3 * v.y;
    ctype const c7 = M.c1 * v.y - M.c2 * v.x;
    ctype const c8 = M.c1 * v.w - M.c4 * v.x;
    ctype const c9 = M.c2 * v.w - M.c4 * v.y;
    ctype const c10 = M.c3 * v.w - M.c4 * v.z;
    ctype const c11 = M.c5 * v.w + M.c8 * v.z - M.c10 * v.x;
    ctype const c12 = M.c6 * v.w + M.c9 * v.z - M.c10 * v.y;
    ctype const c13 = M.c7 * v.w - M.c8 * v.y + M.c9 * v.x;
    ctype const c14 = -M.c5 * v.y + M.c6 * v.x - M.c7 * v.z;
    ctype const c15 = -M.c11 * v.y + M.c12 * v.x - M.c13 * v.z - M.c14 * v.w;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(Vec2dc<T> const& v, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = v.w * M.c0;
    ctype const c5 = -v.x * M.c3 + v.z * M.c1;
    ctype const c6 = -v.y * M.c3 + v.z * M.c2;
    ctype const c7 = v.x * M.c2 - v.y * M.c1;
    ctype const c8 = v.x * M.c4 - v.w * M.c1;
    ctype const c9 = v.y * M.c4 - v.w * M.c2;
    ctype const c10 = v.z * M.c4 - v.w * M.c3;
    ctype const c11 = -v.x * M.c10 + v.z * M.c8 + v.w * M.c5;
    ctype const c12 = -v.y * M.c10 + v.z * M.c9 + v.w * M.c6;
    ctype const c13 = v.x * M.c9 - v.y * M.c8 + v.w * M.c7;
    ctype const c14 = v.x * M.c6 - v.y * M.c5 - v.z * M.c7;
    ctype const c15 = -v.x * M.c12 + v.y * M.c11 + v.z * M.c13 + v.w * M.c14;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(MVec2dc<T> const& M, Scalar2dc<U> s)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> wdg(Scalar2dc<T> s, MVec2dc<U> const& M)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc wdg :: wdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& A,
                                                  MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& A,
                                                  MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c3 + A.c4 * B.c2 - A.c6 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c2 * B.c3 + A.c5 * B.c2 - A.c6 * B.c1;
    ctype const c6 = A.c0 * B.c6 + A.c3 * B.c3 - A.c4 * B.c1 + A.c5 * B.c0;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c2;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc wdg :: wdg(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> wdg(MVec2dc_U<T> const& A,
                                                  MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = -A.c0 * B.c6 + A.c2 * B.c4 + A.c3 * B.c1 + A.c4 * B.c0;
    ctype const c5 = -A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c2 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c5 - A.c1 * B.c4 + A.c3 * B.c3 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c3 + A.c7 * B.c0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc wdg :: wdg(mv_e,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& M,
                                                   PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(M.c0 * ctype(ps));
}

// cga2dc wdg :: wdg(ps,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(PScalar2dc<T> ps,
                                                   MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps) * M.c0);
}

// cga2dc wdg :: wdg(mv_e,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& M,
                                                  TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * t.x;
    ctype const c1 = M.c0 * t.y;
    ctype const c2 = M.c0 * t.z;
    ctype const c3 = M.c0 * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(trivec,mv_e) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(TriVec2dc<T> const& t,
                                                  MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * M.c0;
    ctype const c1 = t.y * M.c0;
    ctype const c2 = t.z * M.c0;
    ctype const c3 = t.w * M.c0;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& M,
                                                  BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * B.vx;
    ctype const c2 = M.c0 * B.vy;
    ctype const c3 = M.c0 * B.vz;
    ctype const c4 = M.c0 * B.mx;
    ctype const c5 = M.c0 * B.my;
    ctype const c6 = M.c0 * B.mz;
    ctype const c7 =
        M.c1 * B.my - M.c2 * B.mx + M.c3 * B.mz - M.c4 * B.vy + M.c5 * B.vx + M.c6 * B.vz;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(BiVec2dc<T> const& B,
                                                  MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = B.vx * M.c0;
    ctype const c2 = B.vy * M.c0;
    ctype const c3 = B.vz * M.c0;
    ctype const c4 = B.mx * M.c0;
    ctype const c5 = B.my * M.c0;
    ctype const c6 = B.mz * M.c0;
    ctype const c7 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c6 - B.mx * M.c2 + B.my * M.c1 + B.mz * M.c3;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& M,
                                                  Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    ctype const c3 = M.c0 * v.w;
    ctype const c4 = M.c1 * v.w + M.c4 * v.z - M.c6 * v.x;
    ctype const c5 = M.c2 * v.w + M.c5 * v.z - M.c6 * v.y;
    ctype const c6 = M.c3 * v.w - M.c4 * v.y + M.c5 * v.x;
    ctype const c7 = -M.c1 * v.y + M.c2 * v.x - M.c3 * v.z;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc wdg :: wdg(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> wdg(Vec2dc<T> const& v,
                                                  MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    ctype const c3 = v.w * M.c0;
    ctype const c4 = -v.x * M.c6 + v.z * M.c4 + v.w * M.c1;
    ctype const c5 = -v.y * M.c6 + v.z * M.c5 + v.w * M.c2;
    ctype const c6 = v.x * M.c5 - v.y * M.c4 + v.w * M.c3;
    ctype const c7 = v.x * M.c2 - v.y * M.c1 - v.z * M.c3;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc wdg :: wdg(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(MVec2dc_E<T> const& M, Scalar2dc<U> s)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(Scalar2dc<T> s, MVec2dc_E<U> const& M)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(MVec2dc_U<T> const& A,
                                                  MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c2 = -A.c1 * B.c2 + A.c2 * B.c1;
    ctype const c3 = A.c0 * B.c1 - A.c1 * B.c0;
    ctype const c4 = A.c0 * B.c3 - A.c3 * B.c0;
    ctype const c5 = A.c1 * B.c3 - A.c3 * B.c1;
    ctype const c6 = A.c2 * B.c3 - A.c3 * B.c2;
    ctype const c7 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c6 + A.c3 * B.c7 -
                     A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2 - A.c7 * B.c3;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(mv_u,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] MVec2dc_U<T> const&,
                                                  [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(ps,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dc<T>,
                                                  [[maybe_unused]] MVec2dc_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(mv_u,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(MVec2dc_U<T> const& M,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-M.c0 * t.y + M.c1 * t.x + M.c2 * t.z + M.c3 * t.w);
}

// cga2dc wdg :: wdg(trivec,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(TriVec2dc<T> const& t,
                                                   MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-t.x * M.c1 + t.y * M.c0 - t.z * M.c2 - t.w * M.c3);
}

// cga2dc wdg :: wdg(mv_u,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(MVec2dc_U<T> const& M,
                                                  BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c0 * B.mz + M.c2 * B.mx + M.c3 * B.vx;
    ctype const c1 = -M.c1 * B.mz + M.c2 * B.my + M.c3 * B.vy;
    ctype const c2 = M.c0 * B.my - M.c1 * B.mx + M.c3 * B.vz;
    ctype const c3 = M.c0 * B.vy - M.c1 * B.vx - M.c2 * B.vz;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(bivec,mv_u) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(BiVec2dc<T> const& B,
                                                  MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c3 + B.mx * M.c2 - B.mz * M.c0;
    ctype const c1 = B.vy * M.c3 + B.my * M.c2 - B.mz * M.c1;
    ctype const c2 = B.vz * M.c3 - B.mx * M.c1 + B.my * M.c0;
    ctype const c3 = -B.vx * M.c1 + B.vy * M.c0 - B.vz * M.c2;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(MVec2dc_U<T> const& M,
                                                  Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -M.c0 * v.z + M.c2 * v.x;
    ctype const c2 = -M.c1 * v.z + M.c2 * v.y;
    ctype const c3 = M.c0 * v.y - M.c1 * v.x;
    ctype const c4 = M.c0 * v.w - M.c3 * v.x;
    ctype const c5 = M.c1 * v.w - M.c3 * v.y;
    ctype const c6 = M.c2 * v.w - M.c3 * v.z;
    ctype const c7 = -M.c4 * v.y + M.c5 * v.x - M.c6 * v.z - M.c7 * v.w;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> wdg(Vec2dc<T> const& v,
                                                  MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -v.x * M.c2 + v.z * M.c0;
    ctype const c2 = -v.y * M.c2 + v.z * M.c1;
    ctype const c3 = v.x * M.c1 - v.y * M.c0;
    ctype const c4 = v.x * M.c3 - v.w * M.c0;
    ctype const c5 = v.y * M.c3 - v.w * M.c1;
    ctype const c6 = v.z * M.c3 - v.w * M.c2;
    ctype const c7 = -v.x * M.c5 + v.y * M.c4 + v.z * M.c6 + v.w * M.c7;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc wdg :: wdg(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> wdg(MVec2dc_U<T> const& M, Scalar2dc<U> s)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc wdg :: wdg(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> wdg(Scalar2dc<T> s, MVec2dc_U<U> const& M)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc wdg :: wdg(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dc<T>,
                                                  [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dc<T>,
                                                  [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec2dc<T> const&,
                                                  [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dc<T>,
                                                  [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec2dc<T> const&,
                                                  [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dc<T>,
                                                  [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec2dc<T> const&,
                                                  [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(PScalar2dc<T> ps, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps) * ctype(s));
}

// cga2dc wdg :: wdg(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(Scalar2dc<T> s, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(s) * ctype(ps));
}

// cga2dc wdg :: wdg(trivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec2dc<T> const&,
                                                  [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(trivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec2dc<T> const&,
                                                  [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(bivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec2dc<T> const&,
                                                  [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc wdg :: wdg(trivec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(TriVec2dc<T> const& t,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-t.x * v.y + t.y * v.x - t.z * v.z - t.w * v.w);
}

// cga2dc wdg :: wdg(vec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(Vec2dc<T> const& v,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-v.x * t.y + v.y * t.x + v.z * t.z + v.w * t.w);
}

// cga2dc wdg :: wdg(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(TriVec2dc<T> const& t, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(Scalar2dc<T> s, TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(bivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> wdg(BiVec2dc<T> const& B1,
                                                   BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(B1.vx * B2.my - B1.vy * B2.mx + B1.vz * B2.mz -
                             B1.mx * B2.vy + B1.my * B2.vx + B1.mz * B2.vz);
}

// cga2dc wdg :: wdg(bivec,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(BiVec2dc<T> const& B,
                                                  Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.w + B.mx * v.z - B.mz * v.x;
    ctype const c1 = B.vy * v.w + B.my * v.z - B.mz * v.y;
    ctype const c2 = B.vz * v.w - B.mx * v.y + B.my * v.x;
    ctype const c3 = -B.vx * v.y + B.vy * v.x - B.vz * v.z;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(vec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> wdg(Vec2dc<T> const& v,
                                                  BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * B.mz + v.z * B.mx + v.w * B.vx;
    ctype const c1 = -v.y * B.mz + v.z * B.my + v.w * B.vy;
    ctype const c2 = v.x * B.my - v.y * B.mx + v.w * B.vz;
    ctype const c3 = v.x * B.vy - v.y * B.vx - v.z * B.vz;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> wdg(BiVec2dc<T> const& B, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc wdg :: wdg(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> wdg(Scalar2dc<T> s, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc wdg :: wdg(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> wdg(Vec2dc<T> const& v1, Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c1 = -v1.y * v2.z + v1.z * v2.y;
    ctype const c2 = v1.x * v2.y - v1.y * v2.x;
    ctype const c3 = v1.x * v2.w - v1.w * v2.x;
    ctype const c4 = v1.y * v2.w - v1.w * v2.y;
    ctype const c5 = v1.z * v2.w - v1.w * v2.z;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc wdg :: wdg(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> wdg(Vec2dc<T> const& v, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> wdg(Scalar2dc<T> s, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc wdg :: wdg(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> wdg(Scalar2dc<T> s1, Scalar2dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s1) * ctype(s2));
}
////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product for complements)
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product, just depending on the outer product (wdg)
// as well as the complement and thus the pseudoscalar of the space
// (this definition does NOT connect directly to the geometric product,
// but to the outer product exclusively)
//
//     rwdg(ul, ur) = l_cmpl(wdg(r_cmpl(ul),r_cmpl(ur))) = l_cmpl(r_cmpl(ul) ^ r_cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// cga2dc rwdg :: rwdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& A, MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c1 * B.c12 + A.c2 * B.c11 + A.c3 * B.c13 +
                     A.c4 * B.c14 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c10 -
                     A.c8 * B.c6 + A.c9 * B.c5 + A.c10 * B.c7 - A.c11 * B.c2 +
                     A.c12 * B.c1 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    ctype const c1 = A.c1 * B.c15 + A.c5 * B.c13 - A.c7 * B.c11 - A.c8 * B.c14 -
                     A.c11 * B.c7 + A.c13 * B.c5 - A.c14 * B.c8 + A.c15 * B.c1;
    ctype const c2 = A.c2 * B.c15 + A.c6 * B.c13 - A.c7 * B.c12 - A.c9 * B.c14 -
                     A.c12 * B.c7 + A.c13 * B.c6 - A.c14 * B.c9 + A.c15 * B.c2;
    ctype const c3 = A.c3 * B.c15 + A.c5 * B.c12 - A.c6 * B.c11 - A.c10 * B.c14 -
                     A.c11 * B.c6 + A.c12 * B.c5 - A.c14 * B.c10 + A.c15 * B.c3;
    ctype const c4 = A.c4 * B.c15 - A.c8 * B.c12 + A.c9 * B.c11 + A.c10 * B.c13 +
                     A.c11 * B.c9 - A.c12 * B.c8 + A.c13 * B.c10 + A.c15 * B.c4;
    ctype const c5 = A.c5 * B.c15 + A.c11 * B.c14 - A.c14 * B.c11 + A.c15 * B.c5;
    ctype const c6 = A.c6 * B.c15 + A.c12 * B.c14 - A.c14 * B.c12 + A.c15 * B.c6;
    ctype const c7 = A.c7 * B.c15 + A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c7;
    ctype const c8 = A.c8 * B.c15 + A.c11 * B.c13 - A.c13 * B.c11 + A.c15 * B.c8;
    ctype const c9 = A.c9 * B.c15 + A.c12 * B.c13 - A.c13 * B.c12 + A.c15 * B.c9;
    ctype const c10 = A.c10 * B.c15 + A.c11 * B.c12 - A.c12 * B.c11 + A.c15 * B.c10;
    ctype const c11 = A.c11 * B.c15 + A.c15 * B.c11;
    ctype const c12 = A.c12 * B.c15 + A.c15 * B.c12;
    ctype const c13 = A.c13 * B.c15 + A.c15 * B.c13;
    ctype const c14 = A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = A.c15 * B.c15;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& A,
                                                 MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c6 - A.c8 * B.c2 +
                     A.c9 * B.c1 + A.c10 * B.c3 + A.c15 * B.c0;
    ctype const c1 = A.c1 * B.c7 - A.c11 * B.c3 + A.c13 * B.c1 - A.c14 * B.c4;
    ctype const c2 = A.c2 * B.c7 - A.c12 * B.c3 + A.c13 * B.c2 - A.c14 * B.c5;
    ctype const c3 = A.c3 * B.c7 - A.c11 * B.c2 + A.c12 * B.c1 - A.c14 * B.c6;
    ctype const c4 = A.c4 * B.c7 + A.c11 * B.c5 - A.c12 * B.c4 + A.c13 * B.c6;
    ctype const c5 = A.c5 * B.c7 + A.c15 * B.c1;
    ctype const c6 = A.c6 * B.c7 + A.c15 * B.c2;
    ctype const c7 = A.c7 * B.c7 + A.c15 * B.c3;
    ctype const c8 = A.c8 * B.c7 + A.c15 * B.c4;
    ctype const c9 = A.c9 * B.c7 + A.c15 * B.c5;
    ctype const c10 = A.c10 * B.c7 + A.c15 * B.c6;
    ctype const c11 = A.c11 * B.c7;
    ctype const c12 = A.c12 * B.c7;
    ctype const c13 = A.c13 * B.c7;
    ctype const c14 = A.c14 * B.c7;
    ctype const c15 = A.c15 * B.c7;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& A,
                                                 MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c10 -
                     A.c4 * B.c6 + A.c5 * B.c5 + A.c6 * B.c7 + A.c7 * B.c0;
    ctype const c1 = A.c1 * B.c13 - A.c3 * B.c11 - A.c4 * B.c14 + A.c7 * B.c1;
    ctype const c2 = A.c2 * B.c13 - A.c3 * B.c12 - A.c5 * B.c14 + A.c7 * B.c2;
    ctype const c3 = A.c1 * B.c12 - A.c2 * B.c11 - A.c6 * B.c14 + A.c7 * B.c3;
    ctype const c4 = -A.c4 * B.c12 + A.c5 * B.c11 + A.c6 * B.c13 + A.c7 * B.c4;
    ctype const c5 = A.c1 * B.c15 + A.c7 * B.c5;
    ctype const c6 = A.c2 * B.c15 + A.c7 * B.c6;
    ctype const c7 = A.c3 * B.c15 + A.c7 * B.c7;
    ctype const c8 = A.c4 * B.c15 + A.c7 * B.c8;
    ctype const c9 = A.c5 * B.c15 + A.c7 * B.c9;
    ctype const c10 = A.c6 * B.c15 + A.c7 * B.c10;
    ctype const c11 = A.c7 * B.c11;
    ctype const c12 = A.c7 * B.c12;
    ctype const c13 = A.c7 * B.c13;
    ctype const c14 = A.c7 * B.c14;
    ctype const c15 = A.c7 * B.c15;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& A,
                                                 MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c6 + A.c4 * B.c7 -
                     A.c11 * B.c1 + A.c12 * B.c0 - A.c13 * B.c2 - A.c14 * B.c3;
    ctype const c1 = A.c5 * B.c6 - A.c7 * B.c4 - A.c8 * B.c7 + A.c15 * B.c0;
    ctype const c2 = A.c6 * B.c6 - A.c7 * B.c5 - A.c9 * B.c7 + A.c15 * B.c1;
    ctype const c3 = A.c5 * B.c5 - A.c6 * B.c4 - A.c10 * B.c7 + A.c15 * B.c2;
    ctype const c4 = -A.c8 * B.c5 + A.c9 * B.c4 + A.c10 * B.c6 + A.c15 * B.c3;
    ctype const c5 = A.c11 * B.c7 - A.c14 * B.c4;
    ctype const c6 = A.c12 * B.c7 - A.c14 * B.c5;
    ctype const c7 = A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c8 = A.c11 * B.c6 - A.c13 * B.c4;
    ctype const c9 = A.c12 * B.c6 - A.c13 * B.c5;
    ctype const c10 = A.c11 * B.c5 - A.c12 * B.c4;
    ctype const c11 = A.c15 * B.c4;
    ctype const c12 = A.c15 * B.c5;
    ctype const c13 = A.c15 * B.c6;
    ctype const c14 = A.c15 * B.c7;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc_U<T> const& A,
                                                 MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c12 + A.c1 * B.c11 + A.c2 * B.c13 + A.c3 * B.c14 -
                     A.c4 * B.c2 + A.c5 * B.c1 - A.c6 * B.c3 - A.c7 * B.c4;
    ctype const c1 = A.c0 * B.c15 - A.c4 * B.c7 + A.c6 * B.c5 - A.c7 * B.c8;
    ctype const c2 = A.c1 * B.c15 - A.c5 * B.c7 + A.c6 * B.c6 - A.c7 * B.c9;
    ctype const c3 = A.c2 * B.c15 - A.c4 * B.c6 + A.c5 * B.c5 - A.c7 * B.c10;
    ctype const c4 = A.c3 * B.c15 + A.c4 * B.c9 - A.c5 * B.c8 + A.c6 * B.c10;
    ctype const c5 = A.c4 * B.c14 - A.c7 * B.c11;
    ctype const c6 = A.c5 * B.c14 - A.c7 * B.c12;
    ctype const c7 = A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c8 = A.c4 * B.c13 - A.c6 * B.c11;
    ctype const c9 = A.c5 * B.c13 - A.c6 * B.c12;
    ctype const c10 = A.c4 * B.c12 - A.c5 * B.c11;
    ctype const c11 = A.c4 * B.c15;
    ctype const c12 = A.c5 * B.c15;
    ctype const c13 = A.c6 * B.c15;
    ctype const c14 = A.c7 * B.c15;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& M, PScalar2dc<U> ps)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps, MVec2dc<U> const& M)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& M,
                                                 TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * t.y + M.c2 * t.x + M.c3 * t.z + M.c4 * t.w;
    ctype const c1 = M.c5 * t.z - M.c7 * t.x - M.c8 * t.w;
    ctype const c2 = M.c6 * t.z - M.c7 * t.y - M.c9 * t.w;
    ctype const c3 = M.c5 * t.y - M.c6 * t.x - M.c10 * t.w;
    ctype const c4 = -M.c8 * t.y + M.c9 * t.x + M.c10 * t.z;
    ctype const c5 = M.c11 * t.w - M.c14 * t.x;
    ctype const c6 = M.c12 * t.w - M.c14 * t.y;
    ctype const c7 = M.c13 * t.w - M.c14 * t.z;
    ctype const c8 = M.c11 * t.z - M.c13 * t.x;
    ctype const c9 = M.c12 * t.z - M.c13 * t.y;
    ctype const c10 = M.c11 * t.y - M.c12 * t.x;
    ctype const c11 = M.c15 * t.x;
    ctype const c12 = M.c15 * t.y;
    ctype const c13 = M.c15 * t.z;
    ctype const c14 = M.c15 * t.w;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(TriVec2dc<T> const& t,
                                                 MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c2 + t.y * M.c1 - t.z * M.c3 - t.w * M.c4;
    ctype const c1 = -t.x * M.c7 + t.z * M.c5 - t.w * M.c8;
    ctype const c2 = -t.y * M.c7 + t.z * M.c6 - t.w * M.c9;
    ctype const c3 = -t.x * M.c6 + t.y * M.c5 - t.w * M.c10;
    ctype const c4 = t.x * M.c9 - t.y * M.c8 + t.z * M.c10;
    ctype const c5 = t.x * M.c14 - t.w * M.c11;
    ctype const c6 = t.y * M.c14 - t.w * M.c12;
    ctype const c7 = t.z * M.c14 - t.w * M.c13;
    ctype const c8 = t.x * M.c13 - t.z * M.c11;
    ctype const c9 = t.y * M.c13 - t.z * M.c12;
    ctype const c10 = t.x * M.c12 - t.y * M.c11;
    ctype const c11 = t.x * M.c15;
    ctype const c12 = t.y * M.c15;
    ctype const c13 = t.z * M.c15;
    ctype const c14 = t.w * M.c15;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& M,
                                                 BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * B.my - M.c6 * B.mx + M.c7 * B.mz - M.c8 * B.vy + M.c9 * B.vx +
                     M.c10 * B.vz;
    ctype const c1 = -M.c11 * B.vz + M.c13 * B.vx - M.c14 * B.mx;
    ctype const c2 = -M.c12 * B.vz + M.c13 * B.vy - M.c14 * B.my;
    ctype const c3 = -M.c11 * B.vy + M.c12 * B.vx - M.c14 * B.mz;
    ctype const c4 = M.c11 * B.my - M.c12 * B.mx + M.c13 * B.mz;
    ctype const c5 = M.c15 * B.vx;
    ctype const c6 = M.c15 * B.vy;
    ctype const c7 = M.c15 * B.vz;
    ctype const c8 = M.c15 * B.mx;
    ctype const c9 = M.c15 * B.my;
    ctype const c10 = M.c15 * B.mz;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(BiVec2dc<T> const& B,
                                                 MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c10 - B.mx * M.c6 +
                     B.my * M.c5 + B.mz * M.c7;
    ctype const c1 = B.vx * M.c13 - B.vz * M.c11 - B.mx * M.c14;
    ctype const c2 = B.vy * M.c13 - B.vz * M.c12 - B.my * M.c14;
    ctype const c3 = B.vx * M.c12 - B.vy * M.c11 - B.mz * M.c14;
    ctype const c4 = -B.mx * M.c12 + B.my * M.c11 + B.mz * M.c13;
    ctype const c5 = B.vx * M.c15;
    ctype const c6 = B.vy * M.c15;
    ctype const c7 = B.vz * M.c15;
    ctype const c8 = B.mx * M.c15;
    ctype const c9 = B.my * M.c15;
    ctype const c10 = B.mz * M.c15;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& M, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * v.y + M.c12 * v.x - M.c13 * v.z - M.c14 * v.w;
    ctype const c1 = M.c15 * v.x;
    ctype const c2 = M.c15 * v.y;
    ctype const c3 = M.c15 * v.z;
    ctype const c4 = M.c15 * v.w;
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rwdg(Vec2dc<T> const& v, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * M.c12 + v.y * M.c11 + v.z * M.c13 + v.w * M.c14;
    ctype const c1 = v.x * M.c15;
    ctype const c2 = v.y * M.c15;
    ctype const c3 = v.z * M.c15;
    ctype const c4 = v.w * M.c15;
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rwdg :: rwdg(mv,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(MVec2dc<T> const& M, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(M.c15 * ctype(s));
}

// cga2dc rwdg :: rwdg(s,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(Scalar2dc<T> s, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s) * M.c15);
}

// cga2dc rwdg :: rwdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& A,
                                                   MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c1 * B.c7 + A.c7 * B.c1;
    ctype const c2 = A.c2 * B.c7 + A.c7 * B.c2;
    ctype const c3 = A.c3 * B.c7 + A.c7 * B.c3;
    ctype const c4 = A.c4 * B.c7 + A.c7 * B.c4;
    ctype const c5 = A.c5 * B.c7 + A.c7 * B.c5;
    ctype const c6 = A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = A.c7 * B.c7;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& A,
                                                   MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c6 - A.c3 * B.c4 - A.c4 * B.c7 + A.c7 * B.c0;
    ctype const c1 = A.c2 * B.c6 - A.c3 * B.c5 - A.c5 * B.c7 + A.c7 * B.c1;
    ctype const c2 = A.c1 * B.c5 - A.c2 * B.c4 - A.c6 * B.c7 + A.c7 * B.c2;
    ctype const c3 = -A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c6 + A.c7 * B.c3;
    ctype const c4 = A.c7 * B.c4;
    ctype const c5 = A.c7 * B.c5;
    ctype const c6 = A.c7 * B.c6;
    ctype const c7 = A.c7 * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rwdg :: rwdg(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rwdg(MVec2dc_U<T> const& A,
                                                   MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c4 * B.c3 + A.c6 * B.c1 - A.c7 * B.c4;
    ctype const c1 = A.c1 * B.c7 - A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c5;
    ctype const c2 = A.c2 * B.c7 - A.c4 * B.c2 + A.c5 * B.c1 - A.c7 * B.c6;
    ctype const c3 = A.c3 * B.c7 + A.c4 * B.c5 - A.c5 * B.c4 + A.c6 * B.c6;
    ctype const c4 = A.c4 * B.c7;
    ctype const c5 = A.c5 * B.c7;
    ctype const c6 = A.c6 * B.c7;
    ctype const c7 = A.c7 * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rwdg :: rwdg(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& M,
                                                   PScalar2dc<U> ps)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps,
                                                   MVec2dc_E<U> const& M)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& M,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * t.z - M.c3 * t.x - M.c4 * t.w;
    ctype const c1 = M.c2 * t.z - M.c3 * t.y - M.c5 * t.w;
    ctype const c2 = M.c1 * t.y - M.c2 * t.x - M.c6 * t.w;
    ctype const c3 = -M.c4 * t.y + M.c5 * t.x + M.c6 * t.z;
    ctype const c4 = M.c7 * t.x;
    ctype const c5 = M.c7 * t.y;
    ctype const c6 = M.c7 * t.z;
    ctype const c7 = M.c7 * t.w;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rwdg :: rwdg(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rwdg(TriVec2dc<T> const& t,
                                                   MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c3 + t.z * M.c1 - t.w * M.c4;
    ctype const c1 = -t.y * M.c3 + t.z * M.c2 - t.w * M.c5;
    ctype const c2 = -t.x * M.c2 + t.y * M.c1 - t.w * M.c6;
    ctype const c3 = t.x * M.c5 - t.y * M.c4 + t.z * M.c6;
    ctype const c4 = t.x * M.c7;
    ctype const c5 = t.y * M.c7;
    ctype const c6 = t.z * M.c7;
    ctype const c7 = t.w * M.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rwdg :: rwdg(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& M,
                                                   BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        M.c1 * B.my - M.c2 * B.mx + M.c3 * B.mz - M.c4 * B.vy + M.c5 * B.vx + M.c6 * B.vz;
    ctype const c1 = M.c7 * B.vx;
    ctype const c2 = M.c7 * B.vy;
    ctype const c3 = M.c7 * B.vz;
    ctype const c4 = M.c7 * B.mx;
    ctype const c5 = M.c7 * B.my;
    ctype const c6 = M.c7 * B.mz;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(BiVec2dc<T> const& B,
                                                   MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c6 - B.mx * M.c2 + B.my * M.c1 + B.mz * M.c3;
    ctype const c1 = B.vx * M.c7;
    ctype const c2 = B.vy * M.c7;
    ctype const c3 = B.vz * M.c7;
    ctype const c4 = B.mx * M.c7;
    ctype const c5 = B.my * M.c7;
    ctype const c6 = B.mz * M.c7;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(mv_e,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& M, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c7 * v.x;
    ctype const c1 = M.c7 * v.y;
    ctype const c2 = M.c7 * v.z;
    ctype const c3 = M.c7 * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(vec,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(Vec2dc<T> const& v, MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c7;
    ctype const c1 = v.y * M.c7;
    ctype const c2 = v.z * M.c7;
    ctype const c3 = v.w * M.c7;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(mv_e,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(MVec2dc_E<T> const& M, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(M.c7 * ctype(s));
}

// cga2dc rwdg :: rwdg(s,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(Scalar2dc<T> s, MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s) * M.c7);
}

// cga2dc rwdg :: rwdg(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(MVec2dc_U<T> const& A,
                                                   MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c6 + A.c3 * B.c7 -
                     A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2 - A.c7 * B.c3;
    ctype const c1 = A.c4 * B.c7 - A.c7 * B.c4;
    ctype const c2 = A.c5 * B.c7 - A.c7 * B.c5;
    ctype const c3 = A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c4 = A.c4 * B.c6 - A.c6 * B.c4;
    ctype const c5 = A.c5 * B.c6 - A.c6 * B.c5;
    ctype const c6 = A.c4 * B.c5 - A.c5 * B.c4;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rwdg(MVec2dc_U<T> const& M,
                                                   PScalar2dc<U> ps)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rwdg :: rwdg(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps,
                                                   MVec2dc_U<U> const& M)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rwdg :: rwdg(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(MVec2dc_U<T> const& M,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c0 * t.y + M.c1 * t.x + M.c2 * t.z + M.c3 * t.w;
    ctype const c1 = M.c4 * t.w - M.c7 * t.x;
    ctype const c2 = M.c5 * t.w - M.c7 * t.y;
    ctype const c3 = M.c6 * t.w - M.c7 * t.z;
    ctype const c4 = M.c4 * t.z - M.c6 * t.x;
    ctype const c5 = M.c5 * t.z - M.c6 * t.y;
    ctype const c6 = M.c4 * t.y - M.c5 * t.x;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rwdg(TriVec2dc<T> const& t,
                                                   MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c1 + t.y * M.c0 - t.z * M.c2 - t.w * M.c3;
    ctype const c1 = t.x * M.c7 - t.w * M.c4;
    ctype const c2 = t.y * M.c7 - t.w * M.c5;
    ctype const c3 = t.z * M.c7 - t.w * M.c6;
    ctype const c4 = t.x * M.c6 - t.z * M.c4;
    ctype const c5 = t.y * M.c6 - t.z * M.c5;
    ctype const c6 = t.x * M.c5 - t.y * M.c4;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rwdg :: rwdg(mv_u,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(MVec2dc_U<T> const& M,
                                                BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * B.vz + M.c6 * B.vx - M.c7 * B.mx;
    ctype const c1 = -M.c5 * B.vz + M.c6 * B.vy - M.c7 * B.my;
    ctype const c2 = -M.c4 * B.vy + M.c5 * B.vx - M.c7 * B.mz;
    ctype const c3 = M.c4 * B.my - M.c5 * B.mx + M.c6 * B.mz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(bivec,mv_u) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(BiVec2dc<T> const& B,
                                                MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c6 - B.vz * M.c4 - B.mx * M.c7;
    ctype const c1 = B.vy * M.c6 - B.vz * M.c5 - B.my * M.c7;
    ctype const c2 = B.vx * M.c5 - B.vy * M.c4 - B.mz * M.c7;
    ctype const c3 = -B.mx * M.c5 + B.my * M.c4 + B.mz * M.c6;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(mv_u,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(MVec2dc_U<T> const& M,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-M.c4 * v.y + M.c5 * v.x - M.c6 * v.z - M.c7 * v.w);
}

// cga2dc rwdg :: rwdg(vec,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(Vec2dc<T> const& v,
                                                   MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-v.x * M.c5 + v.y * M.c4 + v.z * M.c6 + v.w * M.c7);
}

// cga2dc rwdg :: rwdg(mv_u,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] MVec2dc_U<T> const&,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(s,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] MVec2dc_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps1, PScalar2dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps1) * ctype(ps2));
}

// cga2dc rwdg :: rwdg(ps,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * t.x;
    ctype const c1 = ctype(ps) * t.y;
    ctype const c2 = ctype(ps) * t.z;
    ctype const c3 = ctype(ps) * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(trivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rwdg(TriVec2dc<T> const& t,
                                                   PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(ps);
    ctype const c1 = t.y * ctype(ps);
    ctype const c2 = t.z * ctype(ps);
    ctype const c3 = t.w * ctype(ps);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.vx;
    ctype const c1 = ctype(ps) * B.vy;
    ctype const c2 = ctype(ps) * B.vz;
    ctype const c3 = ctype(ps) * B.mx;
    ctype const c4 = ctype(ps) * B.my;
    ctype const c5 = ctype(ps) * B.mz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rwdg :: rwdg(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rwdg(BiVec2dc<T> const& B, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(ps);
    ctype const c1 = B.vy * ctype(ps);
    ctype const c2 = B.vz * ctype(ps);
    ctype const c3 = B.mx * ctype(ps);
    ctype const c4 = B.my * ctype(ps);
    ctype const c5 = B.mz * ctype(ps);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rwdg :: rwdg(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = ctype(ps) * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(Vec2dc<T> const& v, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = v.w * ctype(ps);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(PScalar2dc<T> ps, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(ps) * ctype(s));
}

// cga2dc rwdg :: rwdg(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(Scalar2dc<T> s, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s) * ctype(ps));
}

// cga2dc rwdg :: rwdg(trivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rwdg(TriVec2dc<T> const& t1,
                                                  TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.x * t2.w - t1.w * t2.x;
    ctype const c1 = t1.y * t2.w - t1.w * t2.y;
    ctype const c2 = t1.z * t2.w - t1.w * t2.z;
    ctype const c3 = t1.x * t2.z - t1.z * t2.x;
    ctype const c4 = t1.y * t2.z - t1.z * t2.y;
    ctype const c5 = t1.x * t2.y - t1.y * t2.x;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rwdg :: rwdg(trivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(TriVec2dc<T> const& t,
                                                BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * B.vz + t.z * B.vx - t.w * B.mx;
    ctype const c1 = -t.y * B.vz + t.z * B.vy - t.w * B.my;
    ctype const c2 = -t.x * B.vy + t.y * B.vx - t.w * B.mz;
    ctype const c3 = t.x * B.my - t.y * B.mx + t.z * B.mz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(bivec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rwdg(BiVec2dc<T> const& B,
                                                TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * t.z - B.vz * t.x - B.mx * t.w;
    ctype const c1 = B.vy * t.z - B.vz * t.y - B.my * t.w;
    ctype const c2 = B.vx * t.y - B.vy * t.x - B.mz * t.w;
    ctype const c3 = -B.mx * t.y + B.my * t.x + B.mz * t.z;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rwdg :: rwdg(trivec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(TriVec2dc<T> const& t,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-t.x * v.y + t.y * v.x - t.z * v.z - t.w * v.w);
}

// cga2dc rwdg :: rwdg(vec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(Vec2dc<T> const& v,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-v.x * t.y + v.y * t.x + v.z * t.z + v.w * t.w);
}

// cga2dc rwdg :: rwdg(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] TriVec2dc<T> const&,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg(BiVec2dc<T> const& B1,
                                                   BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(B1.vx * B2.my - B1.vy * B2.mx + B1.vz * B2.mz -
                            B1.mx * B2.vy + B1.my * B2.vx + B1.mz * B2.vz);
}

// cga2dc rwdg :: rwdg(bivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec2dc<T> const&,
                                                   [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(vec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2dc<T> const&,
                                                   [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec2dc<T> const&,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(vec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2dc<T> const&,
                                                   [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2dc<T> const&,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rwdg :: rwdg(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}
////////////////////////////////////////////////////////////////////////////////
// transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// twdg1(ps,vec) = trivector
//     -> identical to the geometric product gpr(ps,vec)
//     -> identical to the commutator product cmt(ps,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> twdg1(PScalar2dc<T> ps, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.y;
    ctype const c1 = ctype(ps) * v.x;
    ctype const c2 = ctype(ps) * v.w;
    ctype const c3 = ctype(ps) * v.z;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// twdg1(vec,ps) = trivector
//     -> identical to the geometric product gpr(vec,ps)
//     -> identical to the commutator product cmt(vec,ps)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> twdg1(Vec2dc<T> const& v, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * ctype(ps);
    ctype const c1 = -v.x * ctype(ps);
    ctype const c2 = -v.w * ctype(ps);
    ctype const c3 = -v.z * ctype(ps);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// twdg1(trivec,vec) = bivector
//     -> identical to the right contraction trivec >> vec
//     -> identical to the left contraction vec << trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> twdg1(TriVec2dc<T> const& t,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * v.z - t.w * v.y;
    ctype const c1 = -t.y * v.z + t.w * v.x;
    ctype const c2 = -t.z * v.z + t.w * v.w;
    ctype const c3 = -t.x * v.w - t.z * v.y;
    ctype const c4 = -t.y * v.w + t.z * v.x;
    ctype const c5 = -t.x * v.x - t.y * v.y;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// twdg1(vec,trivec) = bivector
//     -> identical to the left contraction vec << trivec
//     -> identical to the right contraction trivec >> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> twdg1(Vec2dc<T> const& v,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * t.w - v.z * t.x;
    ctype const c1 = v.x * t.w - v.z * t.y;
    ctype const c2 = -v.z * t.z + v.w * t.w;
    ctype const c3 = -v.y * t.z - v.w * t.x;
    ctype const c4 = v.x * t.z - v.w * t.y;
    ctype const c5 = -v.x * t.x - v.y * t.y;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// twdg1(trivec,bivec) = trivector
//     -> identical to the commutator product cmt(trivec,bivec)
//     -> identical to the negated commutator product cmt(bivec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> twdg1(TriVec2dc<T> const& t,
                                                    BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.vz + t.z * B.vy - t.w * B.my;
    ctype const c1 = t.x * B.vz - t.z * B.vx + t.w * B.mx;
    ctype const c2 = -t.x * B.my + t.y * B.mx - t.z * B.mz;
    ctype const c3 = t.x * B.vy - t.y * B.vx + t.w * B.mz;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// twdg1(bivec,trivec) = trivector
//     -> identical to the commutator product cmt(bivec,trivec)
//     -> identical to the negated commutator product cmt(trivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> twdg1(BiVec2dc<T> const& B,
                                                    TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * t.z + B.vz * t.y + B.my * t.w;
    ctype const c1 = B.vx * t.z - B.vz * t.x - B.mx * t.w;
    ctype const c2 = -B.mx * t.y + B.my * t.x + B.mz * t.z;
    ctype const c3 = B.vx * t.y - B.vy * t.x - B.mz * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// twdg1(bivec,bivec) = bivector
//     -> identical to the commutator product cmt(bivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> twdg1(BiVec2dc<T> const& B1,
                                                   BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vx * B2.mz - B1.vy * B2.vz + B1.vz * B2.vy - B1.mz * B2.vx;
    ctype const c1 = B1.vx * B2.vz + B1.vy * B2.mz - B1.vz * B2.vx - B1.mz * B2.vy;
    ctype const c2 = -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx;
    ctype const c3 = B1.vz * B2.my - B1.mx * B2.mz - B1.my * B2.vz + B1.mz * B2.mx;
    ctype const c4 = -B1.vz * B2.mx + B1.mx * B2.vz - B1.my * B2.mz + B1.mz * B2.my;
    ctype const c5 = B1.vx * B2.mx + B1.vy * B2.my - B1.mx * B2.vx - B1.my * B2.vy;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// twdg1(bivec,vec) = vector
//     -> identical to the commutator product cmt(bivec,vec)
//     -> identical to the negated right contraction bivec >> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> twdg1(BiVec2dc<T> const& B, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.w + B.vz * v.y - B.mx * v.z;
    ctype const c1 = B.vy * v.w - B.vz * v.x - B.my * v.z;
    ctype const c2 = B.vx * v.x + B.vy * v.y - B.mz * v.z;
    ctype const c3 = -B.mx * v.x - B.my * v.y + B.mz * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// twdg1(vec,bivec) = vector
//     -> identical to the commutator product cmt(vec,bivec)
//     -> identical to the negated left contraction vec << bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> twdg1(Vec2dc<T> const& v, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.vz + v.z * B.mx - v.w * B.vx;
    ctype const c1 = v.x * B.vz + v.z * B.my - v.w * B.vy;
    ctype const c2 = -v.x * B.vx - v.y * B.vy + v.z * B.mz;
    ctype const c3 = v.x * B.mx + v.y * B.my - v.w * B.mz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// twdg1(vec,vec) = scalar
//     -> identical to the dot product dot(vec,vec)
//     -> identical to the left contraction vec << vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> twdg1(Vec2dc<T> const& v1,
                                                    Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(v1.x * v2.x + v1.y * v2.y - v1.z * v2.w - v1.w * v2.z);
}

////////////////////////////////////////////////////////////////////////////////
// regressive transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// rtwdg1(trivec,trivec) = pseudoscalar
//     -> identical to the regressive dot product rdot(trivec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rtwdg1(TriVec2dc<T> const& t1,
                                                      TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(t1.x * t2.x + t1.y * t2.y - t1.z * t2.w - t1.w * t2.z);
}

// rtwdg1(trivec,bivec) = trivector
//     -> identical to the left expansion l_expand2dc(trivec,bivec)
//     -> identical to the regressive commutator product rcmt(trivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rtwdg1(TriVec2dc<T> const& t,
                                                     BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.mz + t.z * B.vx + t.w * B.mx;
    ctype const c1 = t.x * B.mz + t.z * B.vy + t.w * B.my;
    ctype const c2 = t.x * B.mx + t.y * B.my + t.z * B.vz;
    ctype const c3 = t.x * B.vx + t.y * B.vy - t.w * B.vz;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// rtwdg1(bivec,trivec) = trivector
//     -> identical to the right expansion r_expand2dc(bivec,trivec)
//     -> identical to the regressive commutator product rcmt(bivec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rtwdg1(BiVec2dc<T> const& B,
                                                     TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * t.z - B.mx * t.w + B.mz * t.y;
    ctype const c1 = -B.vy * t.z - B.my * t.w - B.mz * t.x;
    ctype const c2 = -B.vz * t.z - B.mx * t.x - B.my * t.y;
    ctype const c3 = -B.vx * t.x - B.vy * t.y + B.vz * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// rtwdg1(trivec,vec) = bivector
//     -> identical to the negated left expansion l_expand2dc(trivec,vec)
//     -> identical to the negated right expansion r_expand2dc(vec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rtwdg1(TriVec2dc<T> const& t,
                                                    Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * v.z - t.w * v.x;
    ctype const c1 = -t.x * v.z - t.w * v.y;
    ctype const c2 = -t.x * v.x - t.y * v.y;
    ctype const c3 = -t.y * v.w + t.z * v.x;
    ctype const c4 = t.x * v.w + t.z * v.y;
    ctype const c5 = t.z * v.z - t.w * v.w;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// rtwdg1(vec,trivec) = bivector
//     -> identical to the negated right expansion r_expand2dc(vec,trivec)
//     -> identical to the negated left expansion l_expand2dc(trivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rtwdg1(Vec2dc<T> const& v,
                                                    TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * t.w + v.z * t.y;
    ctype const c1 = -v.y * t.w - v.z * t.x;
    ctype const c2 = -v.x * t.x - v.y * t.y;
    ctype const c3 = v.x * t.z - v.w * t.y;
    ctype const c4 = v.y * t.z + v.w * t.x;
    ctype const c5 = v.z * t.z - v.w * t.w;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// rtwdg1(trivec,s) = vector
//     -> identical to the left expansion l_expand2dc(trivec,s)
//     -> identical to the regressive geometric product rgpr(trivec,s)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rtwdg1(TriVec2dc<T> const& t, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * ctype(s);
    ctype const c1 = -t.x * ctype(s);
    ctype const c2 = t.w * ctype(s);
    ctype const c3 = t.z * ctype(s);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// rtwdg1(s,trivec) = vector
//     -> identical to the right expansion r_expand2dc(s,trivec)
//     -> identical to the regressive geometric product rgpr(s,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rtwdg1(Scalar2dc<T> s, TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * t.y;
    ctype const c1 = ctype(s) * t.x;
    ctype const c2 = -ctype(s) * t.w;
    ctype const c3 = -ctype(s) * t.z;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// rtwdg1(bivec,bivec) = bivector
//     -> identical to the regressive commutator product rcmt(bivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rtwdg1(BiVec2dc<T> const& B1,
                                                    BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vx * B2.vz - B1.vy * B2.mz + B1.vz * B2.vx + B1.mz * B2.vy;
    ctype const c1 = B1.vx * B2.mz - B1.vy * B2.vz + B1.vz * B2.vy - B1.mz * B2.vx;
    ctype const c2 = B1.vx * B2.mx + B1.vy * B2.my - B1.mx * B2.vx - B1.my * B2.vy;
    ctype const c3 = -B1.vz * B2.mx + B1.mx * B2.vz - B1.my * B2.mz + B1.mz * B2.my;
    ctype const c4 = -B1.vz * B2.my + B1.mx * B2.mz + B1.my * B2.vz - B1.mz * B2.mx;
    ctype const c5 = B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// rtwdg1(bivec,vec) = vector
//     -> identical to the regressive commutator product rcmt(bivec,vec)
//     -> identical to the negated regressive commutator product rcmt(vec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rtwdg1(BiVec2dc<T> const& B,
                                                  Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * v.w - B.my * v.z + B.mz * v.y;
    ctype const c1 = B.vx * v.w + B.mx * v.z - B.mz * v.x;
    ctype const c2 = B.vx * v.y - B.vy * v.x + B.vz * v.z;
    ctype const c3 = -B.vz * v.w + B.mx * v.y - B.my * v.x;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// rtwdg1(vec,bivec) = vector
//     -> identical to the regressive commutator product rcmt(vec,bivec)
//     -> identical to the negated regressive commutator product rcmt(bivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rtwdg1(Vec2dc<T> const& v,
                                                  BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my + v.w * B.vy;
    ctype const c1 = v.x * B.mz - v.z * B.mx - v.w * B.vx;
    ctype const c2 = v.x * B.vy - v.y * B.vx - v.z * B.vz;
    ctype const c3 = v.x * B.my - v.y * B.mx + v.w * B.vz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
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
//     operator<<(a,b) = l_contract(a,b) = rwdg( l_dual(a), b )
//
// (holds exactly at every grade; the cga2dc metric is non-degenerate, so there
// is a single contraction pair, no bulk/weight split).
////////////////////////////////////////////////////////////////////////////////

// cga2dc l_contract :: l_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& A,
                                                       MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 - A.c4 * B.c3 +
                     A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c7 + A.c8 * B.c5 + A.c9 * B.c6 -
                     A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c12 + A.c13 * B.c14 +
                     A.c14 * B.c13 - A.c15 * B.c15;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c7 - A.c3 * B.c8 + A.c4 * B.c5 +
                     A.c6 * B.c13 + A.c9 * B.c14 + A.c10 * B.c11 + A.c12 * B.c15;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c7 - A.c3 * B.c9 + A.c4 * B.c6 -
                     A.c5 * B.c13 - A.c8 * B.c14 + A.c10 * B.c12 - A.c11 * B.c15;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c5 + A.c2 * B.c6 - A.c3 * B.c10 +
                     A.c5 * B.c11 + A.c6 * B.c12 - A.c7 * B.c14 + A.c14 * B.c15;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c8 - A.c2 * B.c9 + A.c4 * B.c10 +
                     A.c7 * B.c13 + A.c8 * B.c11 + A.c9 * B.c12 + A.c13 * B.c15;
    ctype const c5 = A.c0 * B.c5 - A.c2 * B.c14 - A.c3 * B.c11 + A.c6 * B.c15;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c14 - A.c3 * B.c12 - A.c5 * B.c15;
    ctype const c7 = A.c0 * B.c7 - A.c3 * B.c13 + A.c4 * B.c14 - A.c10 * B.c15;
    ctype const c8 = A.c0 * B.c8 - A.c2 * B.c13 - A.c4 * B.c11 - A.c9 * B.c15;
    ctype const c9 = A.c0 * B.c9 + A.c1 * B.c13 - A.c4 * B.c12 + A.c8 * B.c15;
    ctype const c10 = A.c0 * B.c10 - A.c1 * B.c11 - A.c2 * B.c12 + A.c7 * B.c15;
    ctype const c11 = A.c0 * B.c11 - A.c2 * B.c15;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c15;
    ctype const c13 = A.c0 * B.c13 + A.c4 * B.c15;
    ctype const c14 = A.c0 * B.c14 + A.c3 * B.c15;
    ctype const c15 = A.c0 * B.c15;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& A,
                                                       MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c5 * B.c4 + A.c6 * B.c5 + A.c7 * B.c3 + A.c8 * B.c1 +
                     A.c9 * B.c2 - A.c10 * B.c6 - A.c15 * B.c7;
    ctype const c1 = A.c2 * B.c3 - A.c3 * B.c4 + A.c4 * B.c1 + A.c12 * B.c7;
    ctype const c2 = -A.c1 * B.c3 - A.c3 * B.c5 + A.c4 * B.c2 - A.c11 * B.c7;
    ctype const c3 = A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c6 + A.c14 * B.c7;
    ctype const c4 = -A.c1 * B.c4 - A.c2 * B.c5 + A.c4 * B.c6 + A.c13 * B.c7;
    ctype const c5 = A.c0 * B.c1 + A.c6 * B.c7;
    ctype const c6 = A.c0 * B.c2 - A.c5 * B.c7;
    ctype const c7 = A.c0 * B.c3 - A.c10 * B.c7;
    ctype const c8 = A.c0 * B.c4 - A.c9 * B.c7;
    ctype const c9 = A.c0 * B.c5 + A.c8 * B.c7;
    ctype const c10 = A.c0 * B.c6 + A.c7 * B.c7;
    ctype const c11 = -A.c2 * B.c7;
    ctype const c12 = A.c1 * B.c7;
    ctype const c13 = A.c4 * B.c7;
    ctype const c14 = A.c3 * B.c7;
    ctype const c15 = A.c0 * B.c7;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& A,
                                                       MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c7 + A.c4 * B.c5 +
                     A.c5 * B.c6 - A.c6 * B.c10 - A.c7 * B.c15;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c13 + A.c5 * B.c14 + A.c6 * B.c11;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c13 - A.c4 * B.c14 + A.c6 * B.c12;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c11 + A.c2 * B.c12 - A.c3 * B.c14;
    ctype const c4 = A.c0 * B.c4 + A.c3 * B.c13 + A.c4 * B.c11 + A.c5 * B.c12;
    ctype const c5 = A.c0 * B.c5 + A.c2 * B.c15;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c15;
    ctype const c7 = A.c0 * B.c7 - A.c6 * B.c15;
    ctype const c8 = A.c0 * B.c8 - A.c5 * B.c15;
    ctype const c9 = A.c0 * B.c9 + A.c4 * B.c15;
    ctype const c10 = A.c0 * B.c10 + A.c3 * B.c15;
    ctype const c11 = A.c0 * B.c11;
    ctype const c12 = A.c0 * B.c12;
    ctype const c13 = A.c0 * B.c13;
    ctype const c14 = A.c0 * B.c14;
    ctype const c15 = A.c0 * B.c15;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& A,
                                                       MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 - A.c3 * B.c3 - A.c4 * B.c2 -
                     A.c11 * B.c4 - A.c12 * B.c5 + A.c13 * B.c7 + A.c14 * B.c6;
    ctype const c1 = A.c0 * B.c0 + A.c6 * B.c6 + A.c9 * B.c7 + A.c10 * B.c4;
    ctype const c2 = A.c0 * B.c1 - A.c5 * B.c6 - A.c8 * B.c7 + A.c10 * B.c5;
    ctype const c3 = A.c0 * B.c2 + A.c5 * B.c4 + A.c6 * B.c5 - A.c7 * B.c7;
    ctype const c4 = A.c0 * B.c3 + A.c7 * B.c6 + A.c8 * B.c4 + A.c9 * B.c5;
    ctype const c5 = -A.c2 * B.c7 - A.c3 * B.c4;
    ctype const c6 = A.c1 * B.c7 - A.c3 * B.c5;
    ctype const c7 = -A.c3 * B.c6 + A.c4 * B.c7;
    ctype const c8 = -A.c2 * B.c6 - A.c4 * B.c4;
    ctype const c9 = A.c1 * B.c6 - A.c4 * B.c5;
    ctype const c10 = -A.c1 * B.c4 - A.c2 * B.c5;
    ctype const c11 = A.c0 * B.c4;
    ctype const c12 = A.c0 * B.c5;
    ctype const c13 = A.c0 * B.c6;
    ctype const c14 = A.c0 * B.c7;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc_U<T> const& A,
                                                       MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c4 - A.c3 * B.c3 -
                     A.c4 * B.c11 - A.c5 * B.c12 + A.c6 * B.c14 + A.c7 * B.c13;
    ctype const c1 = A.c1 * B.c7 - A.c2 * B.c8 + A.c3 * B.c5 + A.c5 * B.c15;
    ctype const c2 = -A.c0 * B.c7 - A.c2 * B.c9 + A.c3 * B.c6 - A.c4 * B.c15;
    ctype const c3 = A.c0 * B.c5 + A.c1 * B.c6 - A.c2 * B.c10 + A.c7 * B.c15;
    ctype const c4 = -A.c0 * B.c8 - A.c1 * B.c9 + A.c3 * B.c10 + A.c6 * B.c15;
    ctype const c5 = -A.c1 * B.c14 - A.c2 * B.c11;
    ctype const c6 = A.c0 * B.c14 - A.c2 * B.c12;
    ctype const c7 = -A.c2 * B.c13 + A.c3 * B.c14;
    ctype const c8 = -A.c1 * B.c13 - A.c3 * B.c11;
    ctype const c9 = A.c0 * B.c13 - A.c3 * B.c12;
    ctype const c10 = -A.c0 * B.c11 - A.c1 * B.c12;
    ctype const c11 = -A.c1 * B.c15;
    ctype const c12 = A.c0 * B.c15;
    ctype const c13 = A.c3 * B.c15;
    ctype const c14 = A.c2 * B.c15;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& M,
                                                       PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c15 * ctype(ps);
    ctype const c1 = M.c12 * ctype(ps);
    ctype const c2 = -M.c11 * ctype(ps);
    ctype const c3 = M.c14 * ctype(ps);
    ctype const c4 = M.c13 * ctype(ps);
    ctype const c5 = M.c6 * ctype(ps);
    ctype const c6 = -M.c5 * ctype(ps);
    ctype const c7 = -M.c10 * ctype(ps);
    ctype const c8 = -M.c9 * ctype(ps);
    ctype const c9 = M.c8 * ctype(ps);
    ctype const c10 = M.c7 * ctype(ps);
    ctype const c11 = -M.c2 * ctype(ps);
    ctype const c12 = M.c1 * ctype(ps);
    ctype const c13 = M.c4 * ctype(ps);
    ctype const c14 = M.c3 * ctype(ps);
    ctype const c15 = M.c0 * ctype(ps);
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(PScalar2dc<T> ps,
                                                       MVec2dc<U> const& M)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& M,
                                                       TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * t.x - M.c12 * t.y + M.c13 * t.w + M.c14 * t.z;
    ctype const c1 = M.c6 * t.z + M.c9 * t.w + M.c10 * t.x;
    ctype const c2 = -M.c5 * t.z - M.c8 * t.w + M.c10 * t.y;
    ctype const c3 = M.c5 * t.x + M.c6 * t.y - M.c7 * t.w;
    ctype const c4 = M.c7 * t.z + M.c8 * t.x + M.c9 * t.y;
    ctype const c5 = -M.c2 * t.w - M.c3 * t.x;
    ctype const c6 = M.c1 * t.w - M.c3 * t.y;
    ctype const c7 = -M.c3 * t.z + M.c4 * t.w;
    ctype const c8 = -M.c2 * t.z - M.c4 * t.x;
    ctype const c9 = M.c1 * t.z - M.c4 * t.y;
    ctype const c10 = -M.c1 * t.x - M.c2 * t.y;
    ctype const c11 = M.c0 * t.x;
    ctype const c12 = M.c0 * t.y;
    ctype const c13 = M.c0 * t.z;
    ctype const c14 = M.c0 * t.w;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(TriVec2dc<T> const& t,
                                                       MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c11 - t.y * M.c12 + t.z * M.c14 + t.w * M.c13;
    ctype const c1 = t.y * M.c15;
    ctype const c2 = -t.x * M.c15;
    ctype const c3 = t.w * M.c15;
    ctype const c4 = t.z * M.c15;
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& M,
                                                       BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * B.mx + M.c6 * B.my + M.c7 * B.vz + M.c8 * B.vx + M.c9 * B.vy -
                     M.c10 * B.mz;
    ctype const c1 = M.c2 * B.vz - M.c3 * B.mx + M.c4 * B.vx;
    ctype const c2 = -M.c1 * B.vz - M.c3 * B.my + M.c4 * B.vy;
    ctype const c3 = M.c1 * B.vx + M.c2 * B.vy - M.c3 * B.mz;
    ctype const c4 = -M.c1 * B.mx - M.c2 * B.my + M.c4 * B.mz;
    ctype const c5 = M.c0 * B.vx;
    ctype const c6 = M.c0 * B.vy;
    ctype const c7 = M.c0 * B.vz;
    ctype const c8 = M.c0 * B.mx;
    ctype const c9 = M.c0 * B.my;
    ctype const c10 = M.c0 * B.mz;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(BiVec2dc<T> const& B,
                                                       MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c8 + B.vy * M.c9 + B.vz * M.c7 + B.mx * M.c5 + B.my * M.c6 -
                     B.mz * M.c10;
    ctype const c1 = B.vy * M.c13 + B.my * M.c14 + B.mz * M.c11;
    ctype const c2 = -B.vx * M.c13 - B.mx * M.c14 + B.mz * M.c12;
    ctype const c3 = B.vx * M.c11 + B.vy * M.c12 - B.vz * M.c14;
    ctype const c4 = B.vz * M.c13 + B.mx * M.c11 + B.my * M.c12;
    ctype const c5 = B.vy * M.c15;
    ctype const c6 = -B.vx * M.c15;
    ctype const c7 = -B.mz * M.c15;
    ctype const c8 = -B.my * M.c15;
    ctype const c9 = B.mx * M.c15;
    ctype const c10 = B.vz * M.c15;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& M,
                                                       Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y - M.c3 * v.w - M.c4 * v.z;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = M.c0 * v.w;
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(Vec2dc<T> const& v,
                                                       MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 - v.z * M.c4 - v.w * M.c3;
    ctype const c1 = v.y * M.c7 - v.z * M.c8 + v.w * M.c5;
    ctype const c2 = -v.x * M.c7 - v.z * M.c9 + v.w * M.c6;
    ctype const c3 = v.x * M.c5 + v.y * M.c6 - v.z * M.c10;
    ctype const c4 = -v.x * M.c8 - v.y * M.c9 + v.w * M.c10;
    ctype const c5 = -v.y * M.c14 - v.z * M.c11;
    ctype const c6 = v.x * M.c14 - v.z * M.c12;
    ctype const c7 = -v.z * M.c13 + v.w * M.c14;
    ctype const c8 = -v.y * M.c13 - v.w * M.c11;
    ctype const c9 = v.x * M.c13 - v.w * M.c12;
    ctype const c10 = -v.x * M.c11 - v.y * M.c12;
    ctype const c11 = -v.y * M.c15;
    ctype const c12 = v.x * M.c15;
    ctype const c13 = v.w * M.c15;
    ctype const c14 = v.z * M.c15;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(MVec2dc<T> const& M,
                                                       Scalar2dc<U> s)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s,
                                                       MVec2dc<U> const& M)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc l_contract :: l_contract(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& A,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c3 + A.c4 * B.c1 +
                     A.c5 * B.c2 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c7;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c7;
    ctype const c3 = A.c0 * B.c3 - A.c6 * B.c7;
    ctype const c4 = A.c0 * B.c4 - A.c5 * B.c7;
    ctype const c5 = A.c0 * B.c5 + A.c4 * B.c7;
    ctype const c6 = A.c0 * B.c6 + A.c3 * B.c7;
    ctype const c7 = A.c0 * B.c7;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& A,
                                                         MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c2 * B.c6 + A.c5 * B.c7 + A.c6 * B.c4;
    ctype const c1 = A.c0 * B.c1 - A.c1 * B.c6 - A.c4 * B.c7 + A.c6 * B.c5;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c4 + A.c2 * B.c5 - A.c3 * B.c7;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c6 + A.c4 * B.c4 + A.c5 * B.c5;
    ctype const c4 = A.c0 * B.c4;
    ctype const c5 = A.c0 * B.c5;
    ctype const c6 = A.c0 * B.c6;
    ctype const c7 = A.c0 * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(MVec2dc_U<T> const& A,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c3 - A.c2 * B.c4 + A.c3 * B.c1 + A.c5 * B.c7;
    ctype const c1 = -A.c0 * B.c3 - A.c2 * B.c5 + A.c3 * B.c2 - A.c4 * B.c7;
    ctype const c2 = A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c6 + A.c7 * B.c7;
    ctype const c3 = -A.c0 * B.c4 - A.c1 * B.c5 + A.c3 * B.c6 + A.c6 * B.c7;
    ctype const c4 = -A.c1 * B.c7;
    ctype const c5 = A.c0 * B.c7;
    ctype const c6 = A.c3 * B.c7;
    ctype const c7 = A.c2 * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& A,
                                                         PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c7 * ctype(ps);
    ctype const c1 = A.c2 * ctype(ps);
    ctype const c2 = -A.c1 * ctype(ps);
    ctype const c3 = -A.c6 * ctype(ps);
    ctype const c4 = -A.c5 * ctype(ps);
    ctype const c5 = A.c4 * ctype(ps);
    ctype const c6 = A.c3 * ctype(ps);
    ctype const c7 = A.c0 * ctype(ps);
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(PScalar2dc<T> ps,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c7;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& M,
                                                         TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * t.z + M.c5 * t.w + M.c6 * t.x;
    ctype const c1 = -M.c1 * t.z - M.c4 * t.w + M.c6 * t.y;
    ctype const c2 = M.c1 * t.x + M.c2 * t.y - M.c3 * t.w;
    ctype const c3 = M.c3 * t.z + M.c4 * t.x + M.c5 * t.y;
    ctype const c4 = M.c0 * t.x;
    ctype const c5 = M.c0 * t.y;
    ctype const c6 = M.c0 * t.z;
    ctype const c7 = M.c0 * t.w;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(TriVec2dc<T> const& t,
                                                         MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * M.c7;
    ctype const c1 = -t.x * M.c7;
    ctype const c2 = t.w * M.c7;
    ctype const c3 = t.z * M.c7;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& M,
                                                         BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        M.c1 * B.mx + M.c2 * B.my + M.c3 * B.vz + M.c4 * B.vx + M.c5 * B.vy - M.c6 * B.mz;
    ctype const c1 = M.c0 * B.vx;
    ctype const c2 = M.c0 * B.vy;
    ctype const c3 = M.c0 * B.vz;
    ctype const c4 = M.c0 * B.mx;
    ctype const c5 = M.c0 * B.my;
    ctype const c6 = M.c0 * B.mz;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(BiVec2dc<T> const& B,
                                                         MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c4 + B.vy * M.c5 + B.vz * M.c3 + B.mx * M.c1 + B.my * M.c2 - B.mz * M.c6;
    ctype const c1 = B.vy * M.c7;
    ctype const c2 = -B.vx * M.c7;
    ctype const c3 = -B.mz * M.c7;
    ctype const c4 = -B.my * M.c7;
    ctype const c5 = B.mx * M.c7;
    ctype const c6 = B.vz * M.c7;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& A,
                                                         Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * v.x;
    ctype const c1 = A.c0 * v.y;
    ctype const c2 = A.c0 * v.z;
    ctype const c3 = A.c0 * v.w;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(Vec2dc<T> const& v,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.c3 - v.z * B.c4 + v.w * B.c1;
    ctype const c1 = -v.x * B.c3 - v.z * B.c5 + v.w * B.c2;
    ctype const c2 = v.x * B.c1 + v.y * B.c2 - v.z * B.c6;
    ctype const c3 = -v.x * B.c4 - v.y * B.c5 + v.w * B.c6;
    ctype const c4 = -v.y * B.c7;
    ctype const c5 = v.x * B.c7;
    ctype const c6 = v.w * B.c7;
    ctype const c7 = v.z * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(MVec2dc_E<T> const& M,
                                                         Scalar2dc<U> s)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s,
                                                         MVec2dc_E<U> const& M)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(MVec2dc_U<T> const& A,
                                                         MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 - A.c2 * B.c3 - A.c3 * B.c2 - A.c4 * B.c4 -
                     A.c5 * B.c5 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c1 = -A.c1 * B.c7 - A.c2 * B.c4;
    ctype const c2 = A.c0 * B.c7 - A.c2 * B.c5;
    ctype const c3 = -A.c2 * B.c6 + A.c3 * B.c7;
    ctype const c4 = -A.c1 * B.c6 - A.c3 * B.c4;
    ctype const c5 = A.c0 * B.c6 - A.c3 * B.c5;
    ctype const c6 = -A.c0 * B.c4 - A.c1 * B.c5;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(MVec2dc_U<T> const& A,
                                                         PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c5 * ctype(ps);
    ctype const c1 = -A.c4 * ctype(ps);
    ctype const c2 = A.c7 * ctype(ps);
    ctype const c3 = A.c6 * ctype(ps);
    ctype const c4 = -A.c1 * ctype(ps);
    ctype const c5 = A.c0 * ctype(ps);
    ctype const c6 = A.c3 * ctype(ps);
    ctype const c7 = A.c2 * ctype(ps);
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(ps,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dc<T>, [[maybe_unused]] MVec2dc_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(MVec2dc_U<T> const& M,
                                                         TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * t.x - M.c5 * t.y + M.c6 * t.w + M.c7 * t.z;
    ctype const c1 = -M.c1 * t.w - M.c2 * t.x;
    ctype const c2 = M.c0 * t.w - M.c2 * t.y;
    ctype const c3 = -M.c2 * t.z + M.c3 * t.w;
    ctype const c4 = -M.c1 * t.z - M.c3 * t.x;
    ctype const c5 = M.c0 * t.z - M.c3 * t.y;
    ctype const c6 = -M.c0 * t.x - M.c1 * t.y;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(TriVec2dc<T> const& t,
                                                         MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c4 - t.y * M.c5 + t.z * M.c7 + t.w * M.c6;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(MVec2dc_U<T> const& M,
                                                         BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.vz - M.c2 * B.mx + M.c3 * B.vx;
    ctype const c1 = -M.c0 * B.vz - M.c2 * B.my + M.c3 * B.vy;
    ctype const c2 = M.c0 * B.vx + M.c1 * B.vy - M.c2 * B.mz;
    ctype const c3 = -M.c0 * B.mx - M.c1 * B.my + M.c3 * B.mz;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(BiVec2dc<T> const& B,
                                                         MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * M.c6 + B.my * M.c7 + B.mz * M.c4;
    ctype const c1 = -B.vx * M.c6 - B.mx * M.c7 + B.mz * M.c5;
    ctype const c2 = B.vx * M.c4 + B.vy * M.c5 - B.vz * M.c7;
    ctype const c3 = B.vz * M.c6 + B.mx * M.c4 + B.my * M.c5;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(MVec2dc_U<T> const& M,
                                                         Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y - M.c2 * v.w - M.c3 * v.z;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator<<(Vec2dc<T> const& v,
                                                         MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 - v.z * M.c3 - v.w * M.c2;
    ctype const c1 = -v.y * M.c7 - v.z * M.c4;
    ctype const c2 = v.x * M.c7 - v.z * M.c5;
    ctype const c3 = -v.z * M.c6 + v.w * M.c7;
    ctype const c4 = -v.y * M.c6 - v.w * M.c4;
    ctype const c5 = v.x * M.c6 - v.w * M.c5;
    ctype const c6 = -v.x * M.c4 - v.y * M.c5;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc l_contract :: l_contract(mv_u,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] MVec2dc_U<T> const&, [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s,
                                                         MVec2dc_U<U> const& M)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc l_contract :: l_contract(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator<<(PScalar2dc<T> ps1,
                                                         PScalar2dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga2dc l_contract :: l_contract(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dc<T>, [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(trivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator<<(TriVec2dc<T> const& t,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * ctype(ps);
    ctype const c1 = -t.x * ctype(ps);
    ctype const c2 = t.w * ctype(ps);
    ctype const c3 = t.z * ctype(ps);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc l_contract :: l_contract(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dc<T>, [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator<<(BiVec2dc<T> const& B,
                                                        PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * ctype(ps);
    ctype const c1 = -B.vx * ctype(ps);
    ctype const c2 = -B.mz * ctype(ps);
    ctype const c3 = -B.my * ctype(ps);
    ctype const c4 = B.mx * ctype(ps);
    ctype const c5 = B.vz * ctype(ps);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc l_contract :: l_contract(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dc<T>, [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(vec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator<<(Vec2dc<T> const& v,
                                                         PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(ps);
    ctype const c1 = v.x * ctype(ps);
    ctype const c2 = v.w * ctype(ps);
    ctype const c3 = v.z * ctype(ps);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc l_contract :: l_contract(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2dc<T>,
                                                         [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s,
                                                          PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(s) * ctype(ps));
}

// cga2dc l_contract :: l_contract(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator<<(TriVec2dc<T> const& t1,
                                                         TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-t1.x * t2.x - t1.y * t2.y + t1.z * t2.w + t1.w * t2.z);
}

// cga2dc l_contract :: l_contract(trivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec2dc<T> const&, [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(bivec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator<<(BiVec2dc<T> const& B,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * t.z + B.my * t.w + B.mz * t.x;
    ctype const c1 = -B.vx * t.z - B.mx * t.w + B.mz * t.y;
    ctype const c2 = B.vx * t.x + B.vy * t.y - B.vz * t.w;
    ctype const c3 = B.vz * t.z + B.mx * t.x + B.my * t.y;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc l_contract :: l_contract(trivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec2dc<T> const&, [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(vec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator<<(Vec2dc<T> const& v,
                                                        TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * t.w - v.z * t.x;
    ctype const c1 = v.x * t.w - v.z * t.y;
    ctype const c2 = -v.z * t.z + v.w * t.w;
    ctype const c3 = -v.y * t.z - v.w * t.x;
    ctype const c4 = v.x * t.z - v.w * t.y;
    ctype const c5 = -v.x * t.x - v.y * t.y;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc l_contract :: l_contract(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec2dc<T> const&, [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s,
                                                         TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc l_contract :: l_contract(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator<<(BiVec2dc<T> const& B1,
                                                         BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(B1.vx * B2.mx + B1.vy * B2.my + B1.vz * B2.vz +
                            B1.mx * B2.vx + B1.my * B2.vy - B1.mz * B2.mz);
}

// cga2dc l_contract :: l_contract(bivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec2dc<T> const&, [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator<<(Vec2dc<T> const& v,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.vz - v.z * B.mx + v.w * B.vx;
    ctype const c1 = -v.x * B.vz - v.z * B.my + v.w * B.vy;
    ctype const c2 = v.x * B.vx + v.y * B.vy - v.z * B.mz;
    ctype const c3 = -v.x * B.mx - v.y * B.my + v.w * B.mz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc l_contract :: l_contract(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec2dc<T> const&, [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s,
                                                        BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc l_contract :: l_contract(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator<<(Vec2dc<T> const& v1,
                                                         Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(v1.x * v2.x + v1.y * v2.y - v1.z * v2.w - v1.w * v2.z);
}

// cga2dc l_contract :: l_contract(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec2dc<T> const&, [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc l_contract :: l_contract(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc l_contract :: l_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator<<(Scalar2dc<T> s1, Scalar2dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s1) * ctype(s2));
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
//     operator>>(a,b) = r_contract(a,b) = rwdg( a, r_dual(b) )
//
// (holds exactly at every grade).
////////////////////////////////////////////////////////////////////////////////

// cga2dc r_contract :: r_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& A,
                                                       MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 - A.c4 * B.c3 +
                     A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c7 + A.c8 * B.c5 + A.c9 * B.c6 -
                     A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c12 + A.c13 * B.c14 +
                     A.c14 * B.c13 - A.c15 * B.c15;
    ctype const c1 = A.c1 * B.c0 - A.c5 * B.c4 - A.c7 * B.c2 + A.c8 * B.c3 +
                     A.c11 * B.c10 + A.c13 * B.c6 + A.c14 * B.c9 - A.c15 * B.c12;
    ctype const c2 = A.c2 * B.c0 - A.c6 * B.c4 + A.c7 * B.c1 + A.c9 * B.c3 +
                     A.c12 * B.c10 - A.c13 * B.c5 - A.c14 * B.c8 + A.c15 * B.c11;
    ctype const c3 = A.c3 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2 + A.c10 * B.c3 +
                     A.c11 * B.c5 + A.c12 * B.c6 - A.c14 * B.c7 - A.c15 * B.c14;
    ctype const c4 = A.c4 * B.c0 + A.c8 * B.c1 + A.c9 * B.c2 - A.c10 * B.c4 +
                     A.c11 * B.c8 + A.c12 * B.c9 + A.c13 * B.c7 - A.c15 * B.c13;
    ctype const c5 = A.c5 * B.c0 - A.c11 * B.c3 - A.c14 * B.c2 + A.c15 * B.c6;
    ctype const c6 = A.c6 * B.c0 - A.c12 * B.c3 + A.c14 * B.c1 - A.c15 * B.c5;
    ctype const c7 = A.c7 * B.c0 - A.c13 * B.c3 + A.c14 * B.c4 - A.c15 * B.c10;
    ctype const c8 = A.c8 * B.c0 - A.c11 * B.c4 - A.c13 * B.c2 - A.c15 * B.c9;
    ctype const c9 = A.c9 * B.c0 - A.c12 * B.c4 + A.c13 * B.c1 + A.c15 * B.c8;
    ctype const c10 = A.c10 * B.c0 - A.c11 * B.c1 - A.c12 * B.c2 + A.c15 * B.c7;
    ctype const c11 = A.c11 * B.c0 + A.c15 * B.c2;
    ctype const c12 = A.c12 * B.c0 - A.c15 * B.c1;
    ctype const c13 = A.c13 * B.c0 - A.c15 * B.c4;
    ctype const c14 = A.c14 * B.c0 - A.c15 * B.c3;
    ctype const c15 = A.c15 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& A,
                                                       MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c5 * B.c4 + A.c6 * B.c5 + A.c7 * B.c3 + A.c8 * B.c1 +
                     A.c9 * B.c2 - A.c10 * B.c6 - A.c15 * B.c7;
    ctype const c1 = A.c1 * B.c0 + A.c11 * B.c6 + A.c13 * B.c2 + A.c14 * B.c5;
    ctype const c2 = A.c2 * B.c0 + A.c12 * B.c6 - A.c13 * B.c1 - A.c14 * B.c4;
    ctype const c3 = A.c3 * B.c0 + A.c11 * B.c1 + A.c12 * B.c2 - A.c14 * B.c3;
    ctype const c4 = A.c4 * B.c0 + A.c11 * B.c4 + A.c12 * B.c5 + A.c13 * B.c3;
    ctype const c5 = A.c5 * B.c0 + A.c15 * B.c2;
    ctype const c6 = A.c6 * B.c0 - A.c15 * B.c1;
    ctype const c7 = A.c7 * B.c0 - A.c15 * B.c6;
    ctype const c8 = A.c8 * B.c0 - A.c15 * B.c5;
    ctype const c9 = A.c9 * B.c0 + A.c15 * B.c4;
    ctype const c10 = A.c10 * B.c0 + A.c15 * B.c3;
    ctype const c11 = A.c11 * B.c0;
    ctype const c12 = A.c12 * B.c0;
    ctype const c13 = A.c13 * B.c0;
    ctype const c14 = A.c14 * B.c0;
    ctype const c15 = A.c15 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& A,
                                                       MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c7 + A.c4 * B.c5 +
                     A.c5 * B.c6 - A.c6 * B.c10 - A.c7 * B.c15;
    ctype const c1 = -A.c1 * B.c4 - A.c3 * B.c2 + A.c4 * B.c3 - A.c7 * B.c12;
    ctype const c2 = -A.c2 * B.c4 + A.c3 * B.c1 + A.c5 * B.c3 + A.c7 * B.c11;
    ctype const c3 = -A.c1 * B.c1 - A.c2 * B.c2 + A.c6 * B.c3 - A.c7 * B.c14;
    ctype const c4 = A.c4 * B.c1 + A.c5 * B.c2 - A.c6 * B.c4 - A.c7 * B.c13;
    ctype const c5 = A.c1 * B.c0 + A.c7 * B.c6;
    ctype const c6 = A.c2 * B.c0 - A.c7 * B.c5;
    ctype const c7 = A.c3 * B.c0 - A.c7 * B.c10;
    ctype const c8 = A.c4 * B.c0 - A.c7 * B.c9;
    ctype const c9 = A.c5 * B.c0 + A.c7 * B.c8;
    ctype const c10 = A.c6 * B.c0 + A.c7 * B.c7;
    ctype const c11 = A.c7 * B.c2;
    ctype const c12 = -A.c7 * B.c1;
    ctype const c13 = -A.c7 * B.c4;
    ctype const c14 = -A.c7 * B.c3;
    ctype const c15 = A.c7 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& A,
                                                       MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 - A.c3 * B.c3 - A.c4 * B.c2 -
                     A.c11 * B.c4 - A.c12 * B.c5 + A.c13 * B.c7 + A.c14 * B.c6;
    ctype const c1 = -A.c5 * B.c3 - A.c7 * B.c1 + A.c8 * B.c2 - A.c15 * B.c5;
    ctype const c2 = -A.c6 * B.c3 + A.c7 * B.c0 + A.c9 * B.c2 + A.c15 * B.c4;
    ctype const c3 = -A.c5 * B.c0 - A.c6 * B.c1 + A.c10 * B.c2 - A.c15 * B.c7;
    ctype const c4 = A.c8 * B.c0 + A.c9 * B.c1 - A.c10 * B.c3 - A.c15 * B.c6;
    ctype const c5 = -A.c11 * B.c2 - A.c14 * B.c1;
    ctype const c6 = -A.c12 * B.c2 + A.c14 * B.c0;
    ctype const c7 = -A.c13 * B.c2 + A.c14 * B.c3;
    ctype const c8 = -A.c11 * B.c3 - A.c13 * B.c1;
    ctype const c9 = -A.c12 * B.c3 + A.c13 * B.c0;
    ctype const c10 = -A.c11 * B.c0 - A.c12 * B.c1;
    ctype const c11 = A.c15 * B.c1;
    ctype const c12 = -A.c15 * B.c0;
    ctype const c13 = -A.c15 * B.c3;
    ctype const c14 = -A.c15 * B.c2;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc_U<T> const& A,
                                                       MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c4 - A.c3 * B.c3 -
                     A.c4 * B.c11 - A.c5 * B.c12 + A.c6 * B.c14 + A.c7 * B.c13;
    ctype const c1 = A.c0 * B.c0 + A.c4 * B.c10 + A.c6 * B.c6 + A.c7 * B.c9;
    ctype const c2 = A.c1 * B.c0 + A.c5 * B.c10 - A.c6 * B.c5 - A.c7 * B.c8;
    ctype const c3 = A.c2 * B.c0 + A.c4 * B.c5 + A.c5 * B.c6 - A.c7 * B.c7;
    ctype const c4 = A.c3 * B.c0 + A.c4 * B.c8 + A.c5 * B.c9 + A.c6 * B.c7;
    ctype const c5 = -A.c4 * B.c3 - A.c7 * B.c2;
    ctype const c6 = -A.c5 * B.c3 + A.c7 * B.c1;
    ctype const c7 = -A.c6 * B.c3 + A.c7 * B.c4;
    ctype const c8 = -A.c4 * B.c4 - A.c6 * B.c2;
    ctype const c9 = -A.c5 * B.c4 + A.c6 * B.c1;
    ctype const c10 = -A.c4 * B.c1 - A.c5 * B.c2;
    ctype const c11 = A.c4 * B.c0;
    ctype const c12 = A.c5 * B.c0;
    ctype const c13 = A.c6 * B.c0;
    ctype const c14 = A.c7 * B.c0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& M,
                                                       PScalar2dc<U> ps)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps,
                                                       MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c15;
    ctype const c1 = -ctype(ps) * M.c12;
    ctype const c2 = ctype(ps) * M.c11;
    ctype const c3 = -ctype(ps) * M.c14;
    ctype const c4 = -ctype(ps) * M.c13;
    ctype const c5 = ctype(ps) * M.c6;
    ctype const c6 = -ctype(ps) * M.c5;
    ctype const c7 = -ctype(ps) * M.c10;
    ctype const c8 = -ctype(ps) * M.c9;
    ctype const c9 = ctype(ps) * M.c8;
    ctype const c10 = ctype(ps) * M.c7;
    ctype const c11 = ctype(ps) * M.c2;
    ctype const c12 = -ctype(ps) * M.c1;
    ctype const c13 = -ctype(ps) * M.c4;
    ctype const c14 = -ctype(ps) * M.c3;
    ctype const c15 = ctype(ps) * M.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& M,
                                                       TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * t.x - M.c12 * t.y + M.c13 * t.w + M.c14 * t.z;
    ctype const c1 = -M.c15 * t.y;
    ctype const c2 = M.c15 * t.x;
    ctype const c3 = -M.c15 * t.w;
    ctype const c4 = -M.c15 * t.z;
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(TriVec2dc<T> const& t,
                                                       MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c11 - t.y * M.c12 + t.z * M.c14 + t.w * M.c13;
    ctype const c1 = t.x * M.c10 + t.z * M.c6 + t.w * M.c9;
    ctype const c2 = t.y * M.c10 - t.z * M.c5 - t.w * M.c8;
    ctype const c3 = t.x * M.c5 + t.y * M.c6 - t.w * M.c7;
    ctype const c4 = t.x * M.c8 + t.y * M.c9 + t.z * M.c7;
    ctype const c5 = -t.x * M.c3 - t.w * M.c2;
    ctype const c6 = -t.y * M.c3 + t.w * M.c1;
    ctype const c7 = -t.z * M.c3 + t.w * M.c4;
    ctype const c8 = -t.x * M.c4 - t.z * M.c2;
    ctype const c9 = -t.y * M.c4 + t.z * M.c1;
    ctype const c10 = -t.x * M.c1 - t.y * M.c2;
    ctype const c11 = t.x * M.c0;
    ctype const c12 = t.y * M.c0;
    ctype const c13 = t.z * M.c0;
    ctype const c14 = t.w * M.c0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& M,
                                                       BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * B.mx + M.c6 * B.my + M.c7 * B.vz + M.c8 * B.vx + M.c9 * B.vy -
                     M.c10 * B.mz;
    ctype const c1 = M.c11 * B.mz + M.c13 * B.vy + M.c14 * B.my;
    ctype const c2 = M.c12 * B.mz - M.c13 * B.vx - M.c14 * B.mx;
    ctype const c3 = M.c11 * B.vx + M.c12 * B.vy - M.c14 * B.vz;
    ctype const c4 = M.c11 * B.mx + M.c12 * B.my + M.c13 * B.vz;
    ctype const c5 = M.c15 * B.vy;
    ctype const c6 = -M.c15 * B.vx;
    ctype const c7 = -M.c15 * B.mz;
    ctype const c8 = -M.c15 * B.my;
    ctype const c9 = M.c15 * B.mx;
    ctype const c10 = M.c15 * B.vz;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(BiVec2dc<T> const& B,
                                                       MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c8 + B.vy * M.c9 + B.vz * M.c7 + B.mx * M.c5 + B.my * M.c6 -
                     B.mz * M.c10;
    ctype const c1 = -B.vx * M.c4 - B.vz * M.c2 + B.mx * M.c3;
    ctype const c2 = -B.vy * M.c4 + B.vz * M.c1 + B.my * M.c3;
    ctype const c3 = -B.vx * M.c1 - B.vy * M.c2 + B.mz * M.c3;
    ctype const c4 = B.mx * M.c1 + B.my * M.c2 - B.mz * M.c4;
    ctype const c5 = B.vx * M.c0;
    ctype const c6 = B.vy * M.c0;
    ctype const c7 = B.vz * M.c0;
    ctype const c8 = B.mx * M.c0;
    ctype const c9 = B.my * M.c0;
    ctype const c10 = B.mz * M.c0;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = 0.0;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& M,
                                                       Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y - M.c3 * v.w - M.c4 * v.z;
    ctype const c1 = -M.c5 * v.w - M.c7 * v.y + M.c8 * v.z;
    ctype const c2 = -M.c6 * v.w + M.c7 * v.x + M.c9 * v.z;
    ctype const c3 = -M.c5 * v.x - M.c6 * v.y + M.c10 * v.z;
    ctype const c4 = M.c8 * v.x + M.c9 * v.y - M.c10 * v.w;
    ctype const c5 = -M.c11 * v.z - M.c14 * v.y;
    ctype const c6 = -M.c12 * v.z + M.c14 * v.x;
    ctype const c7 = -M.c13 * v.z + M.c14 * v.w;
    ctype const c8 = -M.c11 * v.w - M.c13 * v.y;
    ctype const c9 = -M.c12 * v.w + M.c13 * v.x;
    ctype const c10 = -M.c11 * v.x - M.c12 * v.y;
    ctype const c11 = M.c15 * v.y;
    ctype const c12 = -M.c15 * v.x;
    ctype const c13 = -M.c15 * v.w;
    ctype const c14 = -M.c15 * v.z;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(Vec2dc<T> const& v,
                                                       MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 - v.z * M.c4 - v.w * M.c3;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = v.w * M.c0;
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(MVec2dc<T> const& M,
                                                       Scalar2dc<U> s)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator>>(Scalar2dc<T> s,
                                                       MVec2dc<U> const& M)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc r_contract :: r_contract(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& A,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c3 + A.c4 * B.c1 +
                     A.c5 * B.c2 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype const c1 = A.c1 * B.c0 + A.c7 * B.c2;
    ctype const c2 = A.c2 * B.c0 - A.c7 * B.c1;
    ctype const c3 = A.c3 * B.c0 - A.c7 * B.c6;
    ctype const c4 = A.c4 * B.c0 - A.c7 * B.c5;
    ctype const c5 = A.c5 * B.c0 + A.c7 * B.c4;
    ctype const c6 = A.c6 * B.c0 + A.c7 * B.c3;
    ctype const c7 = A.c7 * B.c0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& A,
                                                         MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * B.c3 - A.c3 * B.c1 + A.c4 * B.c2 - A.c7 * B.c5;
    ctype const c1 = -A.c2 * B.c3 + A.c3 * B.c0 + A.c5 * B.c2 + A.c7 * B.c4;
    ctype const c2 = -A.c1 * B.c0 - A.c2 * B.c1 + A.c6 * B.c2 - A.c7 * B.c7;
    ctype const c3 = A.c4 * B.c0 + A.c5 * B.c1 - A.c6 * B.c3 - A.c7 * B.c6;
    ctype const c4 = A.c7 * B.c1;
    ctype const c5 = -A.c7 * B.c0;
    ctype const c6 = -A.c7 * B.c3;
    ctype const c7 = -A.c7 * B.c2;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(MVec2dc_U<T> const& A,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c4 * B.c6 + A.c6 * B.c2 + A.c7 * B.c5;
    ctype const c1 = A.c1 * B.c0 + A.c5 * B.c6 - A.c6 * B.c1 - A.c7 * B.c4;
    ctype const c2 = A.c2 * B.c0 + A.c4 * B.c1 + A.c5 * B.c2 - A.c7 * B.c3;
    ctype const c3 = A.c3 * B.c0 + A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c3;
    ctype const c4 = A.c4 * B.c0;
    ctype const c5 = A.c5 * B.c0;
    ctype const c6 = A.c6 * B.c0;
    ctype const c7 = A.c7 * B.c0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& A,
                                                         PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c7 * ctype(ps);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c7;
    ctype const c1 = ctype(ps) * B.c2;
    ctype const c2 = -ctype(ps) * B.c1;
    ctype const c3 = -ctype(ps) * B.c6;
    ctype const c4 = -ctype(ps) * B.c5;
    ctype const c5 = ctype(ps) * B.c4;
    ctype const c6 = ctype(ps) * B.c3;
    ctype const c7 = ctype(ps) * B.c0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& M,
                                                         TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c7 * t.y;
    ctype const c1 = M.c7 * t.x;
    ctype const c2 = -M.c7 * t.w;
    ctype const c3 = -M.c7 * t.z;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(TriVec2dc<T> const& t,
                                                         MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * M.c6 + t.z * M.c2 + t.w * M.c5;
    ctype const c1 = t.y * M.c6 - t.z * M.c1 - t.w * M.c4;
    ctype const c2 = t.x * M.c1 + t.y * M.c2 - t.w * M.c3;
    ctype const c3 = t.x * M.c4 + t.y * M.c5 + t.z * M.c3;
    ctype const c4 = t.x * M.c0;
    ctype const c5 = t.y * M.c0;
    ctype const c6 = t.z * M.c0;
    ctype const c7 = t.w * M.c0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& M,
                                                         BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        M.c1 * B.mx + M.c2 * B.my + M.c3 * B.vz + M.c4 * B.vx + M.c5 * B.vy - M.c6 * B.mz;
    ctype const c1 = M.c7 * B.vy;
    ctype const c2 = -M.c7 * B.vx;
    ctype const c3 = -M.c7 * B.mz;
    ctype const c4 = -M.c7 * B.my;
    ctype const c5 = M.c7 * B.mx;
    ctype const c6 = M.c7 * B.vz;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(BiVec2dc<T> const& B,
                                                         MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c4 + B.vy * M.c5 + B.vz * M.c3 + B.mx * M.c1 + B.my * M.c2 - B.mz * M.c6;
    ctype const c1 = B.vx * M.c0;
    ctype const c2 = B.vy * M.c0;
    ctype const c3 = B.vz * M.c0;
    ctype const c4 = B.mx * M.c0;
    ctype const c5 = B.my * M.c0;
    ctype const c6 = B.mz * M.c0;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& A,
                                                         Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * v.w - A.c3 * v.y + A.c4 * v.z;
    ctype const c1 = -A.c2 * v.w + A.c3 * v.x + A.c5 * v.z;
    ctype const c2 = -A.c1 * v.x - A.c2 * v.y + A.c6 * v.z;
    ctype const c3 = A.c4 * v.x + A.c5 * v.y - A.c6 * v.w;
    ctype const c4 = A.c7 * v.y;
    ctype const c5 = -A.c7 * v.x;
    ctype const c6 = -A.c7 * v.w;
    ctype const c7 = -A.c7 * v.z;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(Vec2dc<T> const& v,
                                                         MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * B.c0;
    ctype const c1 = v.y * B.c0;
    ctype const c2 = v.z * B.c0;
    ctype const c3 = v.w * B.c0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(MVec2dc_E<T> const& M,
                                                         Scalar2dc<U> s)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(Scalar2dc<T> s,
                                                         MVec2dc_E<U> const& M)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(MVec2dc_U<T> const& A,
                                                         MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 - A.c2 * B.c3 - A.c3 * B.c2 - A.c4 * B.c4 -
                     A.c5 * B.c5 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c1 = -A.c4 * B.c2 - A.c7 * B.c1;
    ctype const c2 = -A.c5 * B.c2 + A.c7 * B.c0;
    ctype const c3 = -A.c6 * B.c2 + A.c7 * B.c3;
    ctype const c4 = -A.c4 * B.c3 - A.c6 * B.c1;
    ctype const c5 = -A.c5 * B.c3 + A.c6 * B.c0;
    ctype const c6 = -A.c4 * B.c0 - A.c5 * B.c1;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(mv_u,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] MVec2dc_U<T> const&, [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps,
                                                         MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c5;
    ctype const c1 = ctype(ps) * B.c4;
    ctype const c2 = -ctype(ps) * B.c7;
    ctype const c3 = -ctype(ps) * B.c6;
    ctype const c4 = ctype(ps) * B.c1;
    ctype const c5 = -ctype(ps) * B.c0;
    ctype const c6 = -ctype(ps) * B.c3;
    ctype const c7 = -ctype(ps) * B.c2;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(MVec2dc_U<T> const& M,
                                                         TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * t.x - M.c5 * t.y + M.c6 * t.w + M.c7 * t.z;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(TriVec2dc<T> const& t,
                                                         MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c4 - t.y * M.c5 + t.z * M.c7 + t.w * M.c6;
    ctype const c1 = -t.x * M.c2 - t.w * M.c1;
    ctype const c2 = -t.y * M.c2 + t.w * M.c0;
    ctype const c3 = -t.z * M.c2 + t.w * M.c3;
    ctype const c4 = -t.x * M.c3 - t.z * M.c1;
    ctype const c5 = -t.y * M.c3 + t.z * M.c0;
    ctype const c6 = -t.x * M.c0 - t.y * M.c1;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(MVec2dc_U<T> const& M,
                                                         BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c4 * B.mz + M.c6 * B.vy + M.c7 * B.my;
    ctype const c1 = M.c5 * B.mz - M.c6 * B.vx - M.c7 * B.mx;
    ctype const c2 = M.c4 * B.vx + M.c5 * B.vy - M.c7 * B.vz;
    ctype const c3 = M.c4 * B.mx + M.c5 * B.my + M.c6 * B.vz;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(BiVec2dc<T> const& B,
                                                         MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c3 - B.vz * M.c1 + B.mx * M.c2;
    ctype const c1 = -B.vy * M.c3 + B.vz * M.c0 + B.my * M.c2;
    ctype const c2 = -B.vx * M.c0 - B.vy * M.c1 + B.mz * M.c2;
    ctype const c3 = B.mx * M.c0 + B.my * M.c1 - B.mz * M.c3;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(MVec2dc_U<T> const& M,
                                                         Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y - M.c2 * v.w - M.c3 * v.z;
    ctype const c1 = -M.c4 * v.z - M.c7 * v.y;
    ctype const c2 = -M.c5 * v.z + M.c7 * v.x;
    ctype const c3 = -M.c6 * v.z + M.c7 * v.w;
    ctype const c4 = -M.c4 * v.w - M.c6 * v.y;
    ctype const c5 = -M.c5 * v.w + M.c6 * v.x;
    ctype const c6 = -M.c4 * v.x - M.c5 * v.y;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator>>(Vec2dc<T> const& v,
                                                         MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 - v.z * M.c3 - v.w * M.c2;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc r_contract :: r_contract(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator>>(MVec2dc_U<T> const& M,
                                                         Scalar2dc<U> s)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc r_contract :: r_contract(s,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dc<T>, [[maybe_unused]] MVec2dc_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps1,
                                                         PScalar2dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga2dc r_contract :: r_contract(ps,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * t.y;
    ctype const c1 = ctype(ps) * t.x;
    ctype const c2 = -ctype(ps) * t.w;
    ctype const c3 = -ctype(ps) * t.z;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc r_contract :: r_contract(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] TriVec2dc<T> const&, [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps,
                                                        BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.vy;
    ctype const c1 = -ctype(ps) * B.vx;
    ctype const c2 = -ctype(ps) * B.mz;
    ctype const c3 = -ctype(ps) * B.my;
    ctype const c4 = ctype(ps) * B.mx;
    ctype const c5 = ctype(ps) * B.vz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc r_contract :: r_contract(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec2dc<T> const&, [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(ps,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps,
                                                         Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.y;
    ctype const c1 = -ctype(ps) * v.x;
    ctype const c2 = -ctype(ps) * v.w;
    ctype const c3 = -ctype(ps) * v.z;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc r_contract :: r_contract(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dc<T> const&, [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> operator>>(PScalar2dc<T> ps,
                                                          Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps) * ctype(s));
}

// cga2dc r_contract :: r_contract(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar2dc<T>,
                                                         [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator>>(TriVec2dc<T> const& t1,
                                                         TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-t1.x * t2.x - t1.y * t2.y + t1.z * t2.w + t1.w * t2.z);
}

// cga2dc r_contract :: r_contract(trivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator>>(TriVec2dc<T> const& t,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * B.mz + t.z * B.vy + t.w * B.my;
    ctype const c1 = t.y * B.mz - t.z * B.vx - t.w * B.mx;
    ctype const c2 = t.x * B.vx + t.y * B.vy - t.w * B.vz;
    ctype const c3 = t.x * B.mx + t.y * B.my + t.z * B.vz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc r_contract :: r_contract(bivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec2dc<T> const&, [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(trivec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator>>(TriVec2dc<T> const& t,
                                                        Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * v.z - t.w * v.y;
    ctype const c1 = -t.y * v.z + t.w * v.x;
    ctype const c2 = -t.z * v.z + t.w * v.w;
    ctype const c3 = -t.x * v.w - t.z * v.y;
    ctype const c4 = -t.y * v.w + t.z * v.x;
    ctype const c5 = -t.x * v.x - t.y * v.y;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc r_contract :: r_contract(vec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dc<T> const&, [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator>>(TriVec2dc<T> const& t,
                                                         Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc r_contract :: r_contract(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dc<T>, [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator>>(BiVec2dc<T> const& B1,
                                                         BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(B1.vx * B2.mx + B1.vy * B2.my + B1.vz * B2.vz +
                            B1.mx * B2.vx + B1.my * B2.vy - B1.mz * B2.mz);
}

// cga2dc r_contract :: r_contract(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator>>(BiVec2dc<T> const& B,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * v.w - B.vz * v.y + B.mx * v.z;
    ctype const c1 = -B.vy * v.w + B.vz * v.x + B.my * v.z;
    ctype const c2 = -B.vx * v.x - B.vy * v.y + B.mz * v.z;
    ctype const c3 = B.mx * v.x + B.my * v.y - B.mz * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc r_contract :: r_contract(vec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dc<T> const&, [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator>>(BiVec2dc<T> const& B,
                                                        Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc r_contract :: r_contract(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dc<T>, [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator>>(Vec2dc<T> const& v1,
                                                         Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(v1.x * v2.x + v1.y * v2.y - v1.z * v2.w - v1.w * v2.z);
}

// cga2dc r_contract :: r_contract(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator>>(Vec2dc<T> const& v, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc r_contract :: r_contract(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dc<T>, [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc r_contract :: r_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator>>(Scalar2dc<T> s1, Scalar2dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s1) * ctype(s2));
}

// SKIP cga2dc l_expand :: l_expand(ps,trivec) -> trivec -- non-primitive product
// 'l_expand' SKIP cga2dc l_expand :: l_expand(ps,bivec) -> bivec -- non-primitive product
// 'l_expand' SKIP cga2dc l_expand :: l_expand(ps,vec) -> vec -- non-primitive product
// 'l_expand' SKIP cga2dc l_expand :: l_expand(trivec,bivec) -> trivec -- non-primitive
// product 'l_expand' SKIP cga2dc l_expand :: l_expand(trivec,vec) -> bivec --
// non-primitive product 'l_expand' SKIP cga2dc l_expand :: l_expand(trivec,s) -> vec --
// non-primitive product 'l_expand' SKIP cga2dc l_expand :: l_expand(bivec,vec) -> trivec
// -- non-primitive product 'l_expand' SKIP cga2dc l_expand :: l_expand(bivec,s) -> bivec
// -- non-primitive product 'l_expand' SKIP cga2dc l_expand :: l_expand(vec,s) -> trivec
// -- non-primitive product 'l_expand' SKIP cga2dc r_expand :: r_expand(trivec,ps) ->
// trivec -- non-primitive product 'r_expand' SKIP cga2dc r_expand ::
// r_expand(bivec,trivec) -> trivec -- non-primitive product 'r_expand' SKIP cga2dc
// r_expand :: r_expand(bivec,ps) -> bivec -- non-primitive product 'r_expand' SKIP cga2dc
// r_expand :: r_expand(vec,bivec) -> trivec -- non-primitive product 'r_expand' SKIP
// cga2dc r_expand :: r_expand(vec,trivec) -> bivec -- non-primitive product 'r_expand'
// SKIP cga2dc r_expand :: r_expand(vec,ps) -> vec -- non-primitive product 'r_expand'
// SKIP cga2dc r_expand :: r_expand(s,vec) -> trivec -- non-primitive product 'r_expand'
// SKIP cga2dc r_expand :: r_expand(s,bivec) -> bivec -- non-primitive product 'r_expand'
// SKIP cga2dc r_expand :: r_expand(s,trivec) -> vec -- non-primitive product 'r_expand'
////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

// cga2dc cmt :: cmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> cmt(MVec2dc<T> const& A, MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c2 * B.c7 + A.c3 * B.c8 - A.c4 * B.c5 + A.c5 * B.c4 +
                     A.c7 * B.c2 - A.c8 * B.c3 + A.c12 * B.c15 - A.c15 * B.c12;
    ctype const c2 = A.c1 * B.c7 + A.c3 * B.c9 - A.c4 * B.c6 + A.c6 * B.c4 - A.c7 * B.c1 -
                     A.c9 * B.c3 - A.c11 * B.c15 + A.c15 * B.c11;
    ctype const c3 = -A.c1 * B.c5 - A.c2 * B.c6 + A.c3 * B.c10 + A.c5 * B.c1 +
                     A.c6 * B.c2 - A.c10 * B.c3 + A.c14 * B.c15 - A.c15 * B.c14;
    ctype const c4 = A.c1 * B.c8 + A.c2 * B.c9 - A.c4 * B.c10 - A.c8 * B.c1 -
                     A.c9 * B.c2 + A.c10 * B.c4 + A.c13 * B.c15 - A.c15 * B.c13;
    ctype const c5 = -A.c1 * B.c3 + A.c3 * B.c1 + A.c5 * B.c10 - A.c6 * B.c7 +
                     A.c7 * B.c6 - A.c10 * B.c5 - A.c12 * B.c14 + A.c14 * B.c12;
    ctype const c6 = -A.c2 * B.c3 + A.c3 * B.c2 + A.c5 * B.c7 + A.c6 * B.c10 -
                     A.c7 * B.c5 - A.c10 * B.c6 + A.c11 * B.c14 - A.c14 * B.c11;
    ctype const c7 = A.c1 * B.c2 - A.c2 * B.c1 - A.c5 * B.c9 + A.c6 * B.c8 - A.c8 * B.c6 +
                     A.c9 * B.c5 + A.c11 * B.c12 - A.c12 * B.c11;
    ctype const c8 = A.c1 * B.c4 - A.c4 * B.c1 + A.c7 * B.c9 - A.c8 * B.c10 -
                     A.c9 * B.c7 + A.c10 * B.c8 + A.c12 * B.c13 - A.c13 * B.c12;
    ctype const c9 = A.c2 * B.c4 - A.c4 * B.c2 - A.c7 * B.c8 + A.c8 * B.c7 -
                     A.c9 * B.c10 + A.c10 * B.c9 - A.c11 * B.c13 + A.c13 * B.c11;
    ctype const c10 = A.c3 * B.c4 - A.c4 * B.c3 + A.c5 * B.c8 + A.c6 * B.c9 -
                      A.c8 * B.c5 - A.c9 * B.c6 - A.c13 * B.c14 + A.c14 * B.c13;
    ctype const c11 = A.c2 * B.c15 - A.c6 * B.c13 + A.c7 * B.c12 + A.c9 * B.c14 -
                      A.c12 * B.c7 + A.c13 * B.c6 - A.c14 * B.c9 - A.c15 * B.c2;
    ctype const c12 = -A.c1 * B.c15 + A.c5 * B.c13 - A.c7 * B.c11 - A.c8 * B.c14 +
                      A.c11 * B.c7 - A.c13 * B.c5 + A.c14 * B.c8 + A.c15 * B.c1;
    ctype const c13 = -A.c4 * B.c15 - A.c8 * B.c12 + A.c9 * B.c11 + A.c10 * B.c13 -
                      A.c11 * B.c9 + A.c12 * B.c8 - A.c13 * B.c10 + A.c15 * B.c4;
    ctype const c14 = -A.c3 * B.c15 + A.c5 * B.c12 - A.c6 * B.c11 - A.c10 * B.c14 +
                      A.c11 * B.c6 - A.c12 * B.c5 + A.c14 * B.c10 + A.c15 * B.c3;
    ctype const c15 = -A.c1 * B.c12 + A.c2 * B.c11 + A.c3 * B.c13 + A.c4 * B.c14 -
                      A.c11 * B.c2 + A.c12 * B.c1 - A.c13 * B.c3 - A.c14 * B.c4;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc cmt :: cmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2dc<T>,
                                                  [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(ps,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> cmt(PScalar2dc<T> ps, TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * t.y;
    ctype const c1 = ctype(ps) * t.x;
    ctype const c2 = -ctype(ps) * t.w;
    ctype const c3 = -ctype(ps) * t.z;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(trivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> cmt(TriVec2dc<T> const& t, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * ctype(ps);
    ctype const c1 = -t.x * ctype(ps);
    ctype const c2 = t.w * ctype(ps);
    ctype const c3 = t.z * ctype(ps);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2dc<T>,
                                                  [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec2dc<T> const&,
                                                  [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(ps,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> cmt(PScalar2dc<T> ps, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.y;
    ctype const c1 = ctype(ps) * v.x;
    ctype const c2 = ctype(ps) * v.w;
    ctype const c3 = ctype(ps) * v.z;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(vec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> cmt(Vec2dc<T> const& v, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * ctype(ps);
    ctype const c1 = -v.x * ctype(ps);
    ctype const c2 = -v.w * ctype(ps);
    ctype const c3 = -v.z * ctype(ps);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(ps,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2dc<T>,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(s,ps) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(trivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> cmt(TriVec2dc<T> const& t1,
                                                 TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t1.y * t2.w + t1.w * t2.y;
    ctype const c1 = t1.x * t2.w - t1.w * t2.x;
    ctype const c2 = t1.x * t2.y - t1.y * t2.x;
    ctype const c3 = t1.y * t2.z - t1.z * t2.y;
    ctype const c4 = -t1.x * t2.z + t1.z * t2.x;
    ctype const c5 = -t1.z * t2.w + t1.w * t2.z;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc cmt :: cmt(trivec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> cmt(TriVec2dc<T> const& t,
                                                  BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.vz + t.z * B.vy - t.w * B.my;
    ctype const c1 = t.x * B.vz - t.z * B.vx + t.w * B.mx;
    ctype const c2 = -t.x * B.my + t.y * B.mx - t.z * B.mz;
    ctype const c3 = t.x * B.vy - t.y * B.vx + t.w * B.mz;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(bivec,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> cmt(BiVec2dc<T> const& B,
                                                  TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * t.z + B.vz * t.y + B.my * t.w;
    ctype const c1 = B.vx * t.z - B.vz * t.x - B.mx * t.w;
    ctype const c2 = -B.mx * t.y + B.my * t.x + B.mz * t.z;
    ctype const c3 = B.vx * t.y - B.vy * t.x - B.mz * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(trivec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> cmt(TriVec2dc<T> const& t,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-t.x * v.y + t.y * v.x - t.z * v.z - t.w * v.w);
}

// cga2dc cmt :: cmt(vec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> cmt(Vec2dc<T> const& v,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-v.x * t.y + v.y * t.x + v.z * t.z + v.w * t.w);
}

// cga2dc cmt :: cmt(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] TriVec2dc<T> const&,
                                                  [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dc<T>,
                                                  [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(bivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> cmt(BiVec2dc<T> const& B1,
                                                 BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vx * B2.mz - B1.vy * B2.vz + B1.vz * B2.vy - B1.mz * B2.vx;
    ctype const c1 = B1.vx * B2.vz + B1.vy * B2.mz - B1.vz * B2.vx - B1.mz * B2.vy;
    ctype const c2 = -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx;
    ctype const c3 = B1.vz * B2.my - B1.mx * B2.mz - B1.my * B2.vz + B1.mz * B2.mx;
    ctype const c4 = -B1.vz * B2.mx + B1.mx * B2.vz - B1.my * B2.mz + B1.mz * B2.my;
    ctype const c5 = B1.vx * B2.mx + B1.vy * B2.my - B1.mx * B2.vx - B1.my * B2.vy;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc cmt :: cmt(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> cmt(BiVec2dc<T> const& B, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.w + B.vz * v.y - B.mx * v.z;
    ctype const c1 = B.vy * v.w - B.vz * v.x - B.my * v.z;
    ctype const c2 = B.vx * v.x + B.vy * v.y - B.mz * v.z;
    ctype const c3 = -B.mx * v.x - B.my * v.y + B.mz * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> cmt(Vec2dc<T> const& v, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.vz + v.z * B.mx - v.w * B.vx;
    ctype const c1 = v.x * B.vz + v.z * B.my - v.w * B.vy;
    ctype const c2 = -v.x * B.vx - v.y * B.vy + v.z * B.mz;
    ctype const c3 = v.x * B.mx + v.y * B.my - v.w * B.mz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc cmt :: cmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec2dc<T> const&,
                                                  [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dc<T>,
                                                  [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> cmt(Vec2dc<T> const& v1, Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c1 = -v1.y * v2.z + v1.z * v2.y;
    ctype const c2 = v1.x * v2.y - v1.y * v2.x;
    ctype const c3 = v1.x * v2.w - v1.w * v2.x;
    ctype const c4 = v1.y * v2.w - v1.w * v2.y;
    ctype const c5 = v1.z * v2.w - v1.w * v2.z;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc cmt :: cmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec2dc<T> const&,
                                                  [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dc<T>,
                                                  [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc cmt :: cmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dc<T>,
                                                  [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}
////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////

// cga2dc rcmt :: rcmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rcmt(MVec2dc<T> const& A, MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * B.c12 + A.c2 * B.c11 + A.c3 * B.c13 + A.c4 * B.c14 -
                     A.c11 * B.c2 + A.c12 * B.c1 - A.c13 * B.c3 - A.c14 * B.c4;
    ctype const c1 = -A.c0 * B.c12 - A.c2 * B.c10 + A.c3 * B.c9 + A.c4 * B.c6 -
                     A.c6 * B.c4 - A.c9 * B.c3 + A.c10 * B.c2 + A.c12 * B.c0;
    ctype const c2 = A.c0 * B.c11 + A.c1 * B.c10 - A.c3 * B.c8 - A.c4 * B.c5 +
                     A.c5 * B.c4 + A.c8 * B.c3 - A.c10 * B.c1 - A.c11 * B.c0;
    ctype const c3 = -A.c0 * B.c14 + A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c7 +
                     A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c3 + A.c14 * B.c0;
    ctype const c4 = -A.c0 * B.c13 + A.c1 * B.c9 - A.c2 * B.c8 + A.c4 * B.c7 -
                     A.c7 * B.c4 + A.c8 * B.c2 - A.c9 * B.c1 + A.c13 * B.c0;
    ctype const c5 = -A.c2 * B.c3 + A.c3 * B.c2 - A.c5 * B.c7 - A.c6 * B.c10 +
                     A.c7 * B.c5 + A.c10 * B.c6 + A.c11 * B.c14 - A.c14 * B.c11;
    ctype const c6 = A.c1 * B.c3 - A.c3 * B.c1 + A.c5 * B.c10 - A.c6 * B.c7 +
                     A.c7 * B.c6 - A.c10 * B.c5 + A.c12 * B.c14 - A.c14 * B.c12;
    ctype const c7 = -A.c3 * B.c4 + A.c4 * B.c3 + A.c5 * B.c8 + A.c6 * B.c9 -
                     A.c8 * B.c5 - A.c9 * B.c6 + A.c13 * B.c14 - A.c14 * B.c13;
    ctype const c8 = -A.c2 * B.c4 + A.c4 * B.c2 - A.c7 * B.c8 + A.c8 * B.c7 -
                     A.c9 * B.c10 + A.c10 * B.c9 + A.c11 * B.c13 - A.c13 * B.c11;
    ctype const c9 = A.c1 * B.c4 - A.c4 * B.c1 - A.c7 * B.c9 + A.c8 * B.c10 +
                     A.c9 * B.c7 - A.c10 * B.c8 + A.c12 * B.c13 - A.c13 * B.c12;
    ctype const c10 = A.c1 * B.c2 - A.c2 * B.c1 + A.c5 * B.c9 - A.c6 * B.c8 +
                      A.c8 * B.c6 - A.c9 * B.c5 + A.c11 * B.c12 - A.c12 * B.c11;
    ctype const c11 = -A.c0 * B.c2 + A.c2 * B.c0 - A.c5 * B.c13 - A.c8 * B.c14 +
                      A.c10 * B.c12 - A.c12 * B.c10 + A.c13 * B.c5 + A.c14 * B.c8;
    ctype const c12 = A.c0 * B.c1 - A.c1 * B.c0 - A.c6 * B.c13 - A.c9 * B.c14 -
                      A.c10 * B.c11 + A.c11 * B.c10 + A.c13 * B.c6 + A.c14 * B.c9;
    ctype const c13 = A.c0 * B.c4 - A.c4 * B.c0 - A.c7 * B.c13 - A.c8 * B.c11 -
                      A.c9 * B.c12 + A.c11 * B.c8 + A.c12 * B.c9 + A.c13 * B.c7;
    ctype const c14 = A.c0 * B.c3 - A.c3 * B.c0 - A.c5 * B.c11 - A.c6 * B.c12 +
                      A.c7 * B.c14 + A.c11 * B.c5 + A.c12 * B.c6 - A.c14 * B.c7;
    ctype const c15 = 0.0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rcmt :: rcmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dc<T>,
                                                   [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dc<T>,
                                                   [[maybe_unused]] TriVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] TriVec2dc<T> const&,
                                                   [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dc<T>,
                                                   [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec2dc<T> const&,
                                                   [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dc<T>,
                                                   [[maybe_unused]] Vec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Vec2dc<T> const&,
                                                   [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dc<T>,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] PScalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(trivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rcmt(TriVec2dc<T> const& t1,
                                                  TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.x * t2.w - t1.w * t2.x;
    ctype const c1 = t1.y * t2.w - t1.w * t2.y;
    ctype const c2 = t1.z * t2.w - t1.w * t2.z;
    ctype const c3 = t1.x * t2.z - t1.z * t2.x;
    ctype const c4 = t1.y * t2.z - t1.z * t2.y;
    ctype const c5 = t1.x * t2.y - t1.y * t2.x;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rcmt :: rcmt(trivec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rcmt(TriVec2dc<T> const& t,
                                                   BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.mz + t.z * B.vx + t.w * B.mx;
    ctype const c1 = t.x * B.mz + t.z * B.vy + t.w * B.my;
    ctype const c2 = t.x * B.mx + t.y * B.my + t.z * B.vz;
    ctype const c3 = t.x * B.vx + t.y * B.vy - t.w * B.vz;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(bivec,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rcmt(BiVec2dc<T> const& B,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * t.z - B.mx * t.w + B.mz * t.y;
    ctype const c1 = -B.vy * t.z - B.my * t.w - B.mz * t.x;
    ctype const c2 = -B.vz * t.z - B.mx * t.x - B.my * t.y;
    ctype const c3 = -B.vx * t.x - B.vy * t.y + B.vz * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(trivec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt(TriVec2dc<T> const& t,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-t.x * v.y + t.y * v.x - t.z * v.z - t.w * v.w);
}

// cga2dc rcmt :: rcmt(vec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt(Vec2dc<T> const& v,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-v.x * t.y + v.y * t.x + v.z * t.z + v.w * t.w);
}

// cga2dc rcmt :: rcmt(trivec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rcmt(TriVec2dc<T> const& t, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * ctype(s);
    ctype const c1 = -t.x * ctype(s);
    ctype const c2 = t.w * ctype(s);
    ctype const c3 = t.z * ctype(s);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(s,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rcmt(Scalar2dc<T> s, TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * t.y;
    ctype const c1 = ctype(s) * t.x;
    ctype const c2 = -ctype(s) * t.w;
    ctype const c3 = -ctype(s) * t.z;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(bivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rcmt(BiVec2dc<T> const& B1,
                                                  BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vx * B2.vz - B1.vy * B2.mz + B1.vz * B2.vx + B1.mz * B2.vy;
    ctype const c1 = B1.vx * B2.mz - B1.vy * B2.vz + B1.vz * B2.vy - B1.mz * B2.vx;
    ctype const c2 = B1.vx * B2.mx + B1.vy * B2.my - B1.mx * B2.vx - B1.my * B2.vy;
    ctype const c3 = -B1.vz * B2.mx + B1.mx * B2.vz - B1.my * B2.mz + B1.mz * B2.my;
    ctype const c4 = -B1.vz * B2.my + B1.mx * B2.mz + B1.my * B2.vz - B1.mz * B2.mx;
    ctype const c5 = B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rcmt :: rcmt(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rcmt(BiVec2dc<T> const& B, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * v.w - B.my * v.z + B.mz * v.y;
    ctype const c1 = B.vx * v.w + B.mx * v.z - B.mz * v.x;
    ctype const c2 = B.vx * v.y - B.vy * v.x + B.vz * v.z;
    ctype const c3 = -B.vz * v.w + B.mx * v.y - B.my * v.x;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rcmt(Vec2dc<T> const& v, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my + v.w * B.vy;
    ctype const c1 = v.x * B.mz - v.z * B.mx - v.w * B.vx;
    ctype const c2 = v.x * B.vy - v.y * B.vx - v.z * B.vz;
    ctype const c3 = v.x * B.my - v.y * B.mx + v.w * B.vz;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec2dc<T> const&,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] BiVec2dc<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}

// cga2dc rcmt :: rcmt(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rcmt(Vec2dc<T> const& v1,
                                                  Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.y * v2.z + v1.z * v2.y;
    ctype const c1 = v1.x * v2.z - v1.z * v2.x;
    ctype const c2 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c3 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c4 = v1.x * v2.w - v1.w * v2.x;
    ctype const c5 = v1.x * v2.y - v1.y * v2.x;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rcmt :: rcmt(vec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rcmt(Vec2dc<T> const& v, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * ctype(s);
    ctype const c1 = -v.x * ctype(s);
    ctype const c2 = -v.w * ctype(s);
    ctype const c3 = -v.z * ctype(s);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(s,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rcmt(Scalar2dc<T> s, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * v.y;
    ctype const c1 = ctype(s) * v.x;
    ctype const c2 = ctype(s) * v.w;
    ctype const c3 = ctype(s) * v.z;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rcmt :: rcmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar2dc<T>,
                                                   [[maybe_unused]] Scalar2dc<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(0.0);
}
////////////////////////////////////////////////////////////////////////////////
// geometric products
//
// The null-pair metric makes single basis products multi-term:
//
//     e3 * e4 = -1 + e34        e34 * e34 = 1
//
// (a diagonal metric never mixes a scalar into a basis-blade product). The
// expressions below are generated via the diagonal-basis detour with an
// always-on direct null-basis cross-check (see ga_prdxpr).
////////////////////////////////////////////////////////////////////////////////

// cga2dc gpr :: gpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& A,
                                                      MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 - A.c4 * B.c3 -
                     A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c7 - A.c8 * B.c5 - A.c9 * B.c6 +
                     A.c10 * B.c10 + A.c11 * B.c11 + A.c12 * B.c12 - A.c13 * B.c14 -
                     A.c14 * B.c13 - A.c15 * B.c15;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c7 + A.c3 * B.c8 - A.c4 * B.c5 +
                     A.c5 * B.c4 - A.c6 * B.c13 + A.c7 * B.c2 - A.c8 * B.c3 -
                     A.c9 * B.c14 - A.c10 * B.c11 - A.c11 * B.c10 + A.c12 * B.c15 -
                     A.c13 * B.c6 - A.c14 * B.c9 - A.c15 * B.c12;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c7 + A.c2 * B.c0 + A.c3 * B.c9 - A.c4 * B.c6 +
                     A.c5 * B.c13 + A.c6 * B.c4 - A.c7 * B.c1 + A.c8 * B.c14 -
                     A.c9 * B.c3 - A.c10 * B.c12 - A.c11 * B.c15 - A.c12 * B.c10 +
                     A.c13 * B.c5 + A.c14 * B.c8 + A.c15 * B.c11;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c5 - A.c2 * B.c6 + A.c3 * B.c0 +
                     A.c3 * B.c10 + A.c5 * B.c1 - A.c5 * B.c11 + A.c6 * B.c2 -
                     A.c6 * B.c12 + A.c7 * B.c14 - A.c10 * B.c3 - A.c11 * B.c5 -
                     A.c12 * B.c6 + A.c14 * B.c7 + A.c14 * B.c15 - A.c15 * B.c14;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c8 + A.c2 * B.c9 + A.c4 * B.c0 -
                     A.c4 * B.c10 - A.c7 * B.c13 - A.c8 * B.c1 - A.c8 * B.c11 -
                     A.c9 * B.c2 - A.c9 * B.c12 + A.c10 * B.c4 - A.c11 * B.c8 -
                     A.c12 * B.c9 - A.c13 * B.c7 + A.c13 * B.c15 - A.c15 * B.c13;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c3 - A.c2 * B.c14 + A.c3 * B.c1 -
                     A.c3 * B.c11 + A.c5 * B.c0 + A.c5 * B.c10 - A.c6 * B.c7 -
                     A.c6 * B.c15 + A.c7 * B.c6 - A.c10 * B.c5 - A.c11 * B.c3 -
                     A.c12 * B.c14 - A.c14 * B.c2 + A.c14 * B.c12 - A.c15 * B.c6;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c14 - A.c2 * B.c3 + A.c3 * B.c2 -
                     A.c3 * B.c12 + A.c5 * B.c7 + A.c5 * B.c15 + A.c6 * B.c0 +
                     A.c6 * B.c10 - A.c7 * B.c5 - A.c10 * B.c6 + A.c11 * B.c14 -
                     A.c12 * B.c3 + A.c14 * B.c1 - A.c14 * B.c11 + A.c15 * B.c5;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c13 +
                     A.c4 * B.c14 - A.c5 * B.c9 + A.c6 * B.c8 + A.c7 * B.c0 -
                     A.c8 * B.c6 + A.c9 * B.c5 + A.c10 * B.c15 + A.c11 * B.c12 -
                     A.c12 * B.c11 - A.c13 * B.c3 + A.c14 * B.c4 + A.c15 * B.c10;
    ctype const c8 = A.c0 * B.c8 + A.c1 * B.c4 - A.c2 * B.c13 - A.c4 * B.c1 -
                     A.c4 * B.c11 + A.c7 * B.c9 + A.c8 * B.c0 - A.c8 * B.c10 -
                     A.c9 * B.c7 + A.c9 * B.c15 + A.c10 * B.c8 - A.c11 * B.c4 +
                     A.c12 * B.c13 - A.c13 * B.c2 - A.c13 * B.c12 + A.c15 * B.c9;
    ctype const c9 = A.c0 * B.c9 + A.c1 * B.c13 + A.c2 * B.c4 - A.c4 * B.c2 -
                     A.c4 * B.c12 - A.c7 * B.c8 + A.c8 * B.c7 - A.c8 * B.c15 +
                     A.c9 * B.c0 - A.c9 * B.c10 + A.c10 * B.c9 - A.c11 * B.c13 -
                     A.c12 * B.c4 + A.c13 * B.c1 + A.c13 * B.c11 - A.c15 * B.c8;
    ctype const c10 = A.c0 * B.c10 - A.c1 * B.c11 - A.c2 * B.c12 + A.c3 * B.c4 -
                      A.c4 * B.c3 + A.c5 * B.c8 + A.c6 * B.c9 - A.c7 * B.c15 -
                      A.c8 * B.c5 - A.c9 * B.c6 + A.c10 * B.c0 - A.c11 * B.c1 -
                      A.c12 * B.c2 - A.c13 * B.c14 + A.c14 * B.c13 - A.c15 * B.c7;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c10 + A.c2 * B.c15 + A.c3 * B.c8 +
                      A.c4 * B.c5 + A.c5 * B.c4 - A.c6 * B.c13 + A.c7 * B.c12 +
                      A.c8 * B.c3 + A.c9 * B.c14 - A.c10 * B.c1 + A.c11 * B.c0 -
                      A.c12 * B.c7 + A.c13 * B.c6 - A.c14 * B.c9 - A.c15 * B.c2;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c15 - A.c2 * B.c10 + A.c3 * B.c9 +
                      A.c4 * B.c6 + A.c5 * B.c13 + A.c6 * B.c4 - A.c7 * B.c11 -
                      A.c8 * B.c14 + A.c9 * B.c3 - A.c10 * B.c2 + A.c11 * B.c7 +
                      A.c12 * B.c0 - A.c13 * B.c5 + A.c14 * B.c8 + A.c15 * B.c1;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c9 - A.c2 * B.c8 + A.c4 * B.c7 -
                      A.c4 * B.c15 + A.c7 * B.c4 - A.c8 * B.c2 - A.c8 * B.c12 +
                      A.c9 * B.c1 + A.c9 * B.c11 + A.c10 * B.c13 - A.c11 * B.c9 +
                      A.c12 * B.c8 + A.c13 * B.c0 - A.c13 * B.c10 + A.c15 * B.c4;
    ctype const c14 = A.c0 * B.c14 + A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c7 -
                      A.c3 * B.c15 - A.c5 * B.c2 + A.c5 * B.c12 + A.c6 * B.c1 -
                      A.c6 * B.c11 - A.c7 * B.c3 - A.c10 * B.c14 + A.c11 * B.c6 -
                      A.c12 * B.c5 + A.c14 * B.c0 + A.c14 * B.c10 + A.c15 * B.c3;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c12 + A.c2 * B.c11 + A.c3 * B.c13 +
                      A.c4 * B.c14 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c10 -
                      A.c8 * B.c6 + A.c9 * B.c5 + A.c10 * B.c7 - A.c11 * B.c2 +
                      A.c12 * B.c1 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& A,
                                                      MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c3 - A.c8 * B.c1 -
                     A.c9 * B.c2 + A.c10 * B.c6 - A.c15 * B.c7;
    ctype const c1 = A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c4 - A.c4 * B.c1 -
                     A.c11 * B.c6 + A.c12 * B.c7 - A.c13 * B.c2 - A.c14 * B.c5;
    ctype const c2 = A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c5 - A.c4 * B.c2 -
                     A.c11 * B.c7 - A.c12 * B.c6 + A.c13 * B.c1 + A.c14 * B.c4;
    ctype const c3 = -A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c0 + A.c3 * B.c6 -
                     A.c11 * B.c1 - A.c12 * B.c2 + A.c14 * B.c3 + A.c14 * B.c7;
    ctype const c4 = A.c1 * B.c4 + A.c2 * B.c5 + A.c4 * B.c0 - A.c4 * B.c6 -
                     A.c11 * B.c4 - A.c12 * B.c5 - A.c13 * B.c3 + A.c13 * B.c7;
    ctype const c5 = A.c0 * B.c1 + A.c5 * B.c0 + A.c5 * B.c6 - A.c6 * B.c3 - A.c6 * B.c7 +
                     A.c7 * B.c2 - A.c10 * B.c1 - A.c15 * B.c2;
    ctype const c6 = A.c0 * B.c2 + A.c5 * B.c3 + A.c5 * B.c7 + A.c6 * B.c0 + A.c6 * B.c6 -
                     A.c7 * B.c1 - A.c10 * B.c2 + A.c15 * B.c1;
    ctype const c7 = A.c0 * B.c3 - A.c5 * B.c5 + A.c6 * B.c4 + A.c7 * B.c0 - A.c8 * B.c2 +
                     A.c9 * B.c1 + A.c10 * B.c7 + A.c15 * B.c6;
    ctype const c8 = A.c0 * B.c4 + A.c7 * B.c5 + A.c8 * B.c0 - A.c8 * B.c6 - A.c9 * B.c3 +
                     A.c9 * B.c7 + A.c10 * B.c4 + A.c15 * B.c5;
    ctype const c9 = A.c0 * B.c5 - A.c7 * B.c4 + A.c8 * B.c3 - A.c8 * B.c7 + A.c9 * B.c0 -
                     A.c9 * B.c6 + A.c10 * B.c5 - A.c15 * B.c4;
    ctype const c10 = A.c0 * B.c6 + A.c5 * B.c4 + A.c6 * B.c5 - A.c7 * B.c7 -
                      A.c8 * B.c1 - A.c9 * B.c2 + A.c10 * B.c0 - A.c15 * B.c3;
    ctype const c11 = -A.c1 * B.c6 + A.c2 * B.c7 + A.c3 * B.c4 + A.c4 * B.c1 +
                      A.c11 * B.c0 - A.c12 * B.c3 + A.c13 * B.c2 - A.c14 * B.c5;
    ctype const c12 = -A.c1 * B.c7 - A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c2 +
                      A.c11 * B.c3 + A.c12 * B.c0 - A.c13 * B.c1 + A.c14 * B.c4;
    ctype const c13 = A.c1 * B.c5 - A.c2 * B.c4 + A.c4 * B.c3 - A.c4 * B.c7 -
                      A.c11 * B.c5 + A.c12 * B.c4 + A.c13 * B.c0 - A.c13 * B.c6;
    ctype const c14 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c3 - A.c3 * B.c7 +
                      A.c11 * B.c2 - A.c12 * B.c1 + A.c14 * B.c0 + A.c14 * B.c6;
    ctype const c15 = A.c0 * B.c7 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c6 -
                      A.c8 * B.c2 + A.c9 * B.c1 + A.c10 * B.c3 + A.c15 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& A,
                                                      MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c7 - A.c4 * B.c5 -
                     A.c5 * B.c6 + A.c6 * B.c10 - A.c7 * B.c15;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c4 - A.c2 * B.c13 + A.c3 * B.c2 -
                     A.c4 * B.c3 - A.c5 * B.c14 - A.c6 * B.c11 - A.c7 * B.c12;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c13 + A.c2 * B.c4 - A.c3 * B.c1 +
                     A.c4 * B.c14 - A.c5 * B.c3 - A.c6 * B.c12 + A.c7 * B.c11;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c1 - A.c1 * B.c11 + A.c2 * B.c2 -
                     A.c2 * B.c12 + A.c3 * B.c14 - A.c6 * B.c3 - A.c7 * B.c14;
    ctype const c4 = A.c0 * B.c4 - A.c3 * B.c13 - A.c4 * B.c1 - A.c4 * B.c11 -
                     A.c5 * B.c2 - A.c5 * B.c12 + A.c6 * B.c4 - A.c7 * B.c13;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c0 + A.c1 * B.c10 - A.c2 * B.c7 -
                     A.c2 * B.c15 + A.c3 * B.c6 - A.c6 * B.c5 - A.c7 * B.c6;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c7 + A.c1 * B.c15 + A.c2 * B.c0 +
                     A.c2 * B.c10 - A.c3 * B.c5 - A.c6 * B.c6 + A.c7 * B.c5;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c0 - A.c4 * B.c6 +
                     A.c5 * B.c5 + A.c6 * B.c15 + A.c7 * B.c10;
    ctype const c8 = A.c0 * B.c8 + A.c3 * B.c9 + A.c4 * B.c0 - A.c4 * B.c10 -
                     A.c5 * B.c7 + A.c5 * B.c15 + A.c6 * B.c8 + A.c7 * B.c9;
    ctype const c9 = A.c0 * B.c9 - A.c3 * B.c8 + A.c4 * B.c7 - A.c4 * B.c15 +
                     A.c5 * B.c0 - A.c5 * B.c10 + A.c6 * B.c9 - A.c7 * B.c8;
    ctype const c10 = A.c0 * B.c10 + A.c1 * B.c8 + A.c2 * B.c9 - A.c3 * B.c15 -
                      A.c4 * B.c5 - A.c5 * B.c6 + A.c6 * B.c0 - A.c7 * B.c7;
    ctype const c11 = A.c0 * B.c11 + A.c1 * B.c4 - A.c2 * B.c13 + A.c3 * B.c12 +
                      A.c4 * B.c3 + A.c5 * B.c14 - A.c6 * B.c1 - A.c7 * B.c2;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c13 + A.c2 * B.c4 - A.c3 * B.c11 -
                      A.c4 * B.c14 + A.c5 * B.c3 - A.c6 * B.c2 + A.c7 * B.c1;
    ctype const c13 = A.c0 * B.c13 + A.c3 * B.c4 - A.c4 * B.c2 - A.c4 * B.c12 +
                      A.c5 * B.c1 + A.c5 * B.c11 + A.c6 * B.c13 + A.c7 * B.c4;
    ctype const c14 = A.c0 * B.c14 - A.c1 * B.c2 + A.c1 * B.c12 + A.c2 * B.c1 -
                      A.c2 * B.c11 - A.c3 * B.c3 - A.c6 * B.c14 + A.c7 * B.c3;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c10 -
                      A.c4 * B.c6 + A.c5 * B.c5 + A.c6 * B.c7 + A.c7 * B.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& A,
                                                      MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 - A.c3 * B.c3 - A.c4 * B.c2 +
                     A.c11 * B.c4 + A.c12 * B.c5 - A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c1 = A.c0 * B.c0 + A.c5 * B.c3 - A.c6 * B.c6 + A.c7 * B.c1 - A.c8 * B.c2 -
                     A.c9 * B.c7 - A.c10 * B.c4 - A.c15 * B.c5;
    ctype const c2 = A.c0 * B.c1 + A.c5 * B.c6 + A.c6 * B.c3 - A.c7 * B.c0 + A.c8 * B.c7 -
                     A.c9 * B.c2 - A.c10 * B.c5 + A.c15 * B.c4;
    ctype const c3 = A.c0 * B.c2 + A.c5 * B.c0 - A.c5 * B.c4 + A.c6 * B.c1 - A.c6 * B.c5 +
                     A.c7 * B.c7 - A.c10 * B.c2 - A.c15 * B.c7;
    ctype const c4 = A.c0 * B.c3 - A.c7 * B.c6 - A.c8 * B.c0 - A.c8 * B.c4 - A.c9 * B.c1 -
                     A.c9 * B.c5 + A.c10 * B.c3 - A.c15 * B.c6;
    ctype const c5 = -A.c1 * B.c2 - A.c2 * B.c7 + A.c3 * B.c0 - A.c3 * B.c4 -
                     A.c11 * B.c2 - A.c12 * B.c7 - A.c14 * B.c1 + A.c14 * B.c5;
    ctype const c6 = A.c1 * B.c7 - A.c2 * B.c2 + A.c3 * B.c1 - A.c3 * B.c5 +
                     A.c11 * B.c7 - A.c12 * B.c2 + A.c14 * B.c0 - A.c14 * B.c4;
    ctype const c7 = A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c6 + A.c4 * B.c7 +
                     A.c11 * B.c5 - A.c12 * B.c4 - A.c13 * B.c2 + A.c14 * B.c3;
    ctype const c8 = A.c1 * B.c3 - A.c2 * B.c6 - A.c4 * B.c0 - A.c4 * B.c4 -
                     A.c11 * B.c3 + A.c12 * B.c6 - A.c13 * B.c1 - A.c13 * B.c5;
    ctype const c9 = A.c1 * B.c6 + A.c2 * B.c3 - A.c4 * B.c1 - A.c4 * B.c5 -
                     A.c11 * B.c6 - A.c12 * B.c3 + A.c13 * B.c0 + A.c13 * B.c4;
    ctype const c10 = -A.c1 * B.c4 - A.c2 * B.c5 + A.c3 * B.c3 - A.c4 * B.c2 -
                      A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c7 + A.c14 * B.c6;
    ctype const c11 = A.c0 * B.c4 + A.c5 * B.c3 - A.c6 * B.c6 + A.c7 * B.c5 +
                      A.c8 * B.c2 + A.c9 * B.c7 - A.c10 * B.c0 - A.c15 * B.c1;
    ctype const c12 = A.c0 * B.c5 + A.c5 * B.c6 + A.c6 * B.c3 - A.c7 * B.c4 -
                      A.c8 * B.c7 + A.c9 * B.c2 - A.c10 * B.c1 + A.c15 * B.c0;
    ctype const c13 = A.c0 * B.c6 + A.c7 * B.c3 - A.c8 * B.c1 - A.c8 * B.c5 +
                      A.c9 * B.c0 + A.c9 * B.c4 + A.c10 * B.c6 + A.c15 * B.c3;
    ctype const c14 = A.c0 * B.c7 - A.c5 * B.c1 + A.c5 * B.c5 + A.c6 * B.c0 -
                      A.c6 * B.c4 - A.c7 * B.c2 - A.c10 * B.c7 + A.c15 * B.c2;
    ctype const c15 = -A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c6 + A.c4 * B.c7 -
                      A.c11 * B.c1 + A.c12 * B.c0 - A.c13 * B.c2 - A.c14 * B.c3;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& A,
                                                      MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c4 - A.c3 * B.c3 +
                     A.c4 * B.c11 + A.c5 * B.c12 - A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c1 = A.c0 * B.c0 - A.c1 * B.c7 + A.c2 * B.c8 - A.c3 * B.c5 -
                     A.c4 * B.c10 + A.c5 * B.c15 - A.c6 * B.c6 - A.c7 * B.c9;
    ctype const c2 = A.c0 * B.c7 + A.c1 * B.c0 + A.c2 * B.c9 - A.c3 * B.c6 -
                     A.c4 * B.c15 - A.c5 * B.c10 + A.c6 * B.c5 + A.c7 * B.c8;
    ctype const c3 = -A.c0 * B.c5 - A.c1 * B.c6 + A.c2 * B.c0 + A.c2 * B.c10 -
                     A.c4 * B.c5 - A.c5 * B.c6 + A.c7 * B.c7 + A.c7 * B.c15;
    ctype const c4 = A.c0 * B.c8 + A.c1 * B.c9 + A.c3 * B.c0 - A.c3 * B.c10 -
                     A.c4 * B.c8 - A.c5 * B.c9 - A.c6 * B.c7 + A.c6 * B.c15;
    ctype const c5 = -A.c0 * B.c3 - A.c1 * B.c14 + A.c2 * B.c1 - A.c2 * B.c11 -
                     A.c4 * B.c3 - A.c5 * B.c14 - A.c7 * B.c2 + A.c7 * B.c12;
    ctype const c6 = A.c0 * B.c14 - A.c1 * B.c3 + A.c2 * B.c2 - A.c2 * B.c12 +
                     A.c4 * B.c14 - A.c5 * B.c3 + A.c7 * B.c1 - A.c7 * B.c11;
    ctype const c7 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c13 + A.c3 * B.c14 +
                     A.c4 * B.c12 - A.c5 * B.c11 - A.c6 * B.c3 + A.c7 * B.c4;
    ctype const c8 = A.c0 * B.c4 - A.c1 * B.c13 - A.c3 * B.c1 - A.c3 * B.c11 -
                     A.c4 * B.c4 + A.c5 * B.c13 - A.c6 * B.c2 - A.c6 * B.c12;
    ctype const c9 = A.c0 * B.c13 + A.c1 * B.c4 - A.c3 * B.c2 - A.c3 * B.c12 -
                     A.c4 * B.c13 - A.c5 * B.c4 + A.c6 * B.c1 + A.c6 * B.c11;
    ctype const c10 = -A.c0 * B.c11 - A.c1 * B.c12 + A.c2 * B.c4 - A.c3 * B.c3 -
                      A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c14 + A.c7 * B.c13;
    ctype const c11 = -A.c0 * B.c10 + A.c1 * B.c15 + A.c2 * B.c8 + A.c3 * B.c5 +
                      A.c4 * B.c0 - A.c5 * B.c7 + A.c6 * B.c6 - A.c7 * B.c9;
    ctype const c12 = -A.c0 * B.c15 - A.c1 * B.c10 + A.c2 * B.c9 + A.c3 * B.c6 +
                      A.c4 * B.c7 + A.c5 * B.c0 - A.c6 * B.c5 + A.c7 * B.c8;
    ctype const c13 = A.c0 * B.c9 - A.c1 * B.c8 + A.c3 * B.c7 - A.c3 * B.c15 -
                      A.c4 * B.c9 + A.c5 * B.c8 + A.c6 * B.c0 - A.c6 * B.c10;
    ctype const c14 = A.c0 * B.c6 - A.c1 * B.c5 - A.c2 * B.c7 - A.c2 * B.c15 +
                      A.c4 * B.c6 - A.c5 * B.c5 + A.c7 * B.c0 + A.c7 * B.c10;
    ctype const c15 = -A.c0 * B.c12 + A.c1 * B.c11 + A.c2 * B.c13 + A.c3 * B.c14 -
                      A.c4 * B.c2 + A.c5 * B.c1 - A.c6 * B.c3 - A.c7 * B.c4;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& M,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c15 * ctype(ps);
    ctype const c1 = M.c12 * ctype(ps);
    ctype const c2 = -M.c11 * ctype(ps);
    ctype const c3 = M.c14 * ctype(ps);
    ctype const c4 = M.c13 * ctype(ps);
    ctype const c5 = -M.c6 * ctype(ps);
    ctype const c6 = M.c5 * ctype(ps);
    ctype const c7 = M.c10 * ctype(ps);
    ctype const c8 = M.c9 * ctype(ps);
    ctype const c9 = -M.c8 * ctype(ps);
    ctype const c10 = -M.c7 * ctype(ps);
    ctype const c11 = M.c2 * ctype(ps);
    ctype const c12 = -M.c1 * ctype(ps);
    ctype const c13 = -M.c4 * ctype(ps);
    ctype const c14 = -M.c3 * ctype(ps);
    ctype const c15 = M.c0 * ctype(ps);
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c15;
    ctype const c1 = -ctype(ps) * M.c12;
    ctype const c2 = ctype(ps) * M.c11;
    ctype const c3 = -ctype(ps) * M.c14;
    ctype const c4 = -ctype(ps) * M.c13;
    ctype const c5 = -ctype(ps) * M.c6;
    ctype const c6 = ctype(ps) * M.c5;
    ctype const c7 = ctype(ps) * M.c10;
    ctype const c8 = ctype(ps) * M.c9;
    ctype const c9 = -ctype(ps) * M.c8;
    ctype const c10 = -ctype(ps) * M.c7;
    ctype const c11 = -ctype(ps) * M.c2;
    ctype const c12 = ctype(ps) * M.c1;
    ctype const c13 = ctype(ps) * M.c4;
    ctype const c14 = ctype(ps) * M.c3;
    ctype const c15 = ctype(ps) * M.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& M,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c11 * t.x + M.c12 * t.y - M.c13 * t.w - M.c14 * t.z;
    ctype const c1 = -M.c6 * t.z - M.c9 * t.w - M.c10 * t.x - M.c15 * t.y;
    ctype const c2 = M.c5 * t.z + M.c8 * t.w - M.c10 * t.y + M.c15 * t.x;
    ctype const c3 = -M.c5 * t.x - M.c6 * t.y + M.c7 * t.w - M.c15 * t.w;
    ctype const c4 = -M.c7 * t.z - M.c8 * t.x - M.c9 * t.y - M.c15 * t.z;
    ctype const c5 = -M.c2 * t.w - M.c3 * t.x - M.c12 * t.w + M.c14 * t.y;
    ctype const c6 = M.c1 * t.w - M.c3 * t.y + M.c11 * t.w - M.c14 * t.x;
    ctype const c7 = -M.c3 * t.z + M.c4 * t.w + M.c11 * t.y - M.c12 * t.x;
    ctype const c8 = -M.c2 * t.z - M.c4 * t.x + M.c12 * t.z - M.c13 * t.y;
    ctype const c9 = M.c1 * t.z - M.c4 * t.y - M.c11 * t.z + M.c13 * t.x;
    ctype const c10 = -M.c1 * t.x - M.c2 * t.y - M.c13 * t.w + M.c14 * t.z;
    ctype const c11 = M.c0 * t.x - M.c6 * t.z + M.c7 * t.y + M.c9 * t.w;
    ctype const c12 = M.c0 * t.y + M.c5 * t.z - M.c7 * t.x - M.c8 * t.w;
    ctype const c13 = M.c0 * t.z - M.c8 * t.y + M.c9 * t.x + M.c10 * t.z;
    ctype const c14 = M.c0 * t.w + M.c5 * t.y - M.c6 * t.x - M.c10 * t.w;
    ctype const c15 = -M.c1 * t.y + M.c2 * t.x + M.c3 * t.z + M.c4 * t.w;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * M.c11 + t.y * M.c12 - t.z * M.c14 - t.w * M.c13;
    ctype const c1 = -t.x * M.c10 + t.y * M.c15 - t.z * M.c6 - t.w * M.c9;
    ctype const c2 = -t.x * M.c15 - t.y * M.c10 + t.z * M.c5 + t.w * M.c8;
    ctype const c3 = -t.x * M.c5 - t.y * M.c6 + t.w * M.c7 + t.w * M.c15;
    ctype const c4 = -t.x * M.c8 - t.y * M.c9 - t.z * M.c7 + t.z * M.c15;
    ctype const c5 = -t.x * M.c3 - t.y * M.c14 - t.w * M.c2 + t.w * M.c12;
    ctype const c6 = t.x * M.c14 - t.y * M.c3 + t.w * M.c1 - t.w * M.c11;
    ctype const c7 = t.x * M.c12 - t.y * M.c11 - t.z * M.c3 + t.w * M.c4;
    ctype const c8 = -t.x * M.c4 + t.y * M.c13 - t.z * M.c2 - t.z * M.c12;
    ctype const c9 = -t.x * M.c13 - t.y * M.c4 + t.z * M.c1 + t.z * M.c11;
    ctype const c10 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c14 + t.w * M.c13;
    ctype const c11 = t.x * M.c0 - t.y * M.c7 + t.z * M.c6 - t.w * M.c9;
    ctype const c12 = t.x * M.c7 + t.y * M.c0 - t.z * M.c5 + t.w * M.c8;
    ctype const c13 = -t.x * M.c9 + t.y * M.c8 + t.z * M.c0 - t.z * M.c10;
    ctype const c14 = t.x * M.c6 - t.y * M.c5 + t.w * M.c0 + t.w * M.c10;
    ctype const c15 = -t.x * M.c2 + t.y * M.c1 - t.z * M.c3 - t.w * M.c4;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& M,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.vz - M.c8 * B.vx -
                     M.c9 * B.vy + M.c10 * B.mz;
    ctype const c1 = -M.c2 * B.vz + M.c3 * B.mx - M.c4 * B.vx - M.c11 * B.mz -
                     M.c13 * B.vy - M.c14 * B.my;
    ctype const c2 = M.c1 * B.vz + M.c3 * B.my - M.c4 * B.vy - M.c12 * B.mz +
                     M.c13 * B.vx + M.c14 * B.mx;
    ctype const c3 = -M.c1 * B.vx - M.c2 * B.vy + M.c3 * B.mz - M.c11 * B.vx -
                     M.c12 * B.vy + M.c14 * B.vz;
    ctype const c4 = M.c1 * B.mx + M.c2 * B.my - M.c4 * B.mz - M.c11 * B.mx -
                     M.c12 * B.my - M.c13 * B.vz;
    ctype const c5 = M.c0 * B.vx + M.c5 * B.mz - M.c6 * B.vz + M.c7 * B.vy -
                     M.c10 * B.vx - M.c15 * B.vy;
    ctype const c6 = M.c0 * B.vy + M.c5 * B.vz + M.c6 * B.mz - M.c7 * B.vx -
                     M.c10 * B.vy + M.c15 * B.vx;
    ctype const c7 = M.c0 * B.vz - M.c5 * B.my + M.c6 * B.mx - M.c8 * B.vy + M.c9 * B.vx +
                     M.c15 * B.mz;
    ctype const c8 = M.c0 * B.mx + M.c7 * B.my - M.c8 * B.mz - M.c9 * B.vz +
                     M.c10 * B.mx + M.c15 * B.my;
    ctype const c9 = M.c0 * B.my - M.c7 * B.mx + M.c8 * B.vz - M.c9 * B.mz +
                     M.c10 * B.my - M.c15 * B.mx;
    ctype const c10 = M.c0 * B.mz + M.c5 * B.mx + M.c6 * B.my - M.c8 * B.vx -
                      M.c9 * B.vy - M.c15 * B.vz;
    ctype const c11 = -M.c1 * B.mz + M.c3 * B.mx + M.c4 * B.vx - M.c12 * B.vz +
                      M.c13 * B.vy - M.c14 * B.my;
    ctype const c12 = -M.c2 * B.mz + M.c3 * B.my + M.c4 * B.vy + M.c11 * B.vz -
                      M.c13 * B.vx + M.c14 * B.mx;
    ctype const c13 = M.c1 * B.my - M.c2 * B.mx + M.c4 * B.vz - M.c11 * B.my +
                      M.c12 * B.mx - M.c13 * B.mz;
    ctype const c14 = M.c1 * B.vy - M.c2 * B.vx - M.c3 * B.vz + M.c11 * B.vy -
                      M.c12 * B.vx + M.c14 * B.mz;
    ctype const c15 = M.c5 * B.my - M.c6 * B.mx + M.c7 * B.mz - M.c8 * B.vy +
                      M.c9 * B.vx + M.c10 * B.vz;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c7 - B.mx * M.c5 -
                     B.my * M.c6 + B.mz * M.c10;
    ctype const c1 = B.vx * M.c4 - B.vy * M.c13 + B.vz * M.c2 - B.mx * M.c3 -
                     B.my * M.c14 - B.mz * M.c11;
    ctype const c2 = B.vx * M.c13 + B.vy * M.c4 - B.vz * M.c1 + B.mx * M.c14 -
                     B.my * M.c3 - B.mz * M.c12;
    ctype const c3 = B.vx * M.c1 - B.vx * M.c11 + B.vy * M.c2 - B.vy * M.c12 +
                     B.vz * M.c14 - B.mz * M.c3;
    ctype const c4 = -B.vz * M.c13 - B.mx * M.c1 - B.mx * M.c11 - B.my * M.c2 -
                     B.my * M.c12 + B.mz * M.c4;
    ctype const c5 = B.vx * M.c0 + B.vx * M.c10 - B.vy * M.c7 - B.vy * M.c15 +
                     B.vz * M.c6 - B.mz * M.c5;
    ctype const c6 = B.vx * M.c7 + B.vx * M.c15 + B.vy * M.c0 + B.vy * M.c10 -
                     B.vz * M.c5 - B.mz * M.c6;
    ctype const c7 = -B.vx * M.c9 + B.vy * M.c8 + B.vz * M.c0 - B.mx * M.c6 +
                     B.my * M.c5 + B.mz * M.c15;
    ctype const c8 = B.vz * M.c9 + B.mx * M.c0 - B.mx * M.c10 - B.my * M.c7 +
                     B.my * M.c15 + B.mz * M.c8;
    ctype const c9 = -B.vz * M.c8 + B.mx * M.c7 - B.mx * M.c15 + B.my * M.c0 -
                     B.my * M.c10 + B.mz * M.c9;
    ctype const c10 = B.vx * M.c8 + B.vy * M.c9 - B.vz * M.c15 - B.mx * M.c5 -
                      B.my * M.c6 + B.mz * M.c0;
    ctype const c11 = B.vx * M.c4 - B.vy * M.c13 + B.vz * M.c12 + B.mx * M.c3 +
                      B.my * M.c14 - B.mz * M.c1;
    ctype const c12 = B.vx * M.c13 + B.vy * M.c4 - B.vz * M.c11 - B.mx * M.c14 +
                      B.my * M.c3 - B.mz * M.c2;
    ctype const c13 = B.vz * M.c4 - B.mx * M.c2 - B.mx * M.c12 + B.my * M.c1 +
                      B.my * M.c11 + B.mz * M.c13;
    ctype const c14 = -B.vx * M.c2 + B.vx * M.c12 + B.vy * M.c1 - B.vy * M.c11 -
                      B.vz * M.c3 - B.mz * M.c14;
    ctype const c15 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c10 - B.mx * M.c6 +
                      B.my * M.c5 + B.mz * M.c7;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& M,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y - M.c3 * v.w - M.c4 * v.z;
    ctype const c1 = M.c0 * v.x + M.c5 * v.w + M.c7 * v.y - M.c8 * v.z;
    ctype const c2 = M.c0 * v.y + M.c6 * v.w - M.c7 * v.x - M.c9 * v.z;
    ctype const c3 = M.c0 * v.z + M.c5 * v.x + M.c6 * v.y - M.c10 * v.z;
    ctype const c4 = M.c0 * v.w - M.c8 * v.x - M.c9 * v.y + M.c10 * v.w;
    ctype const c5 = -M.c1 * v.z + M.c3 * v.x - M.c11 * v.z - M.c14 * v.y;
    ctype const c6 = -M.c2 * v.z + M.c3 * v.y - M.c12 * v.z + M.c14 * v.x;
    ctype const c7 = M.c1 * v.y - M.c2 * v.x - M.c13 * v.z + M.c14 * v.w;
    ctype const c8 = M.c1 * v.w - M.c4 * v.x - M.c11 * v.w - M.c13 * v.y;
    ctype const c9 = M.c2 * v.w - M.c4 * v.y - M.c12 * v.w + M.c13 * v.x;
    ctype const c10 = M.c3 * v.w - M.c4 * v.z - M.c11 * v.x - M.c12 * v.y;
    ctype const c11 = M.c5 * v.w + M.c8 * v.z - M.c10 * v.x - M.c15 * v.y;
    ctype const c12 = M.c6 * v.w + M.c9 * v.z - M.c10 * v.y + M.c15 * v.x;
    ctype const c13 = M.c7 * v.w - M.c8 * v.y + M.c9 * v.x + M.c15 * v.w;
    ctype const c14 = -M.c5 * v.y + M.c6 * v.x - M.c7 * v.z + M.c15 * v.z;
    ctype const c15 = -M.c11 * v.y + M.c12 * v.x - M.c13 * v.z - M.c14 * v.w;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 - v.z * M.c4 - v.w * M.c3;
    ctype const c1 = v.x * M.c0 - v.y * M.c7 + v.z * M.c8 - v.w * M.c5;
    ctype const c2 = v.x * M.c7 + v.y * M.c0 + v.z * M.c9 - v.w * M.c6;
    ctype const c3 = -v.x * M.c5 - v.y * M.c6 + v.z * M.c0 + v.z * M.c10;
    ctype const c4 = v.x * M.c8 + v.y * M.c9 + v.w * M.c0 - v.w * M.c10;
    ctype const c5 = -v.x * M.c3 - v.y * M.c14 + v.z * M.c1 - v.z * M.c11;
    ctype const c6 = v.x * M.c14 - v.y * M.c3 + v.z * M.c2 - v.z * M.c12;
    ctype const c7 = v.x * M.c2 - v.y * M.c1 - v.z * M.c13 + v.w * M.c14;
    ctype const c8 = v.x * M.c4 - v.y * M.c13 - v.w * M.c1 - v.w * M.c11;
    ctype const c9 = v.x * M.c13 + v.y * M.c4 - v.w * M.c2 - v.w * M.c12;
    ctype const c10 = -v.x * M.c11 - v.y * M.c12 + v.z * M.c4 - v.w * M.c3;
    ctype const c11 = -v.x * M.c10 + v.y * M.c15 + v.z * M.c8 + v.w * M.c5;
    ctype const c12 = -v.x * M.c15 - v.y * M.c10 + v.z * M.c9 + v.w * M.c6;
    ctype const c13 = v.x * M.c9 - v.y * M.c8 + v.w * M.c7 - v.w * M.c15;
    ctype const c14 = v.x * M.c6 - v.y * M.c5 - v.z * M.c7 - v.z * M.c15;
    ctype const c15 = -v.x * M.c12 + v.y * M.c11 + v.z * M.c13 + v.w * M.c14;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(MVec2dc<T> const& M, Scalar2dc<U> s)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator*(Scalar2dc<T> s, MVec2dc<U> const& M)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc gpr :: gpr(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& A,
                                                        MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c3 - A.c4 * B.c1 -
                     A.c5 * B.c2 + A.c6 * B.c6 - A.c7 * B.c7;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 + A.c1 * B.c6 - A.c2 * B.c3 - A.c2 * B.c7 +
                     A.c3 * B.c2 - A.c6 * B.c1 - A.c7 * B.c2;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c1 * B.c7 + A.c2 * B.c0 + A.c2 * B.c6 -
                     A.c3 * B.c1 - A.c6 * B.c2 + A.c7 * B.c1;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c4 = A.c0 * B.c4 + A.c3 * B.c5 + A.c4 * B.c0 - A.c4 * B.c6 - A.c5 * B.c3 +
                     A.c5 * B.c7 + A.c6 * B.c4 + A.c7 * B.c5;
    ctype const c5 = A.c0 * B.c5 - A.c3 * B.c4 + A.c4 * B.c3 - A.c4 * B.c7 + A.c5 * B.c0 -
                     A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c4;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c4 + A.c2 * B.c5 - A.c3 * B.c7 - A.c4 * B.c1 -
                     A.c5 * B.c2 + A.c6 * B.c0 - A.c7 * B.c3;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& A,
                                                        MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c3 - A.c2 * B.c6 + A.c3 * B.c1 - A.c4 * B.c2 -
                     A.c5 * B.c7 - A.c6 * B.c4 - A.c7 * B.c5;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c6 + A.c2 * B.c3 - A.c3 * B.c0 + A.c4 * B.c7 -
                     A.c5 * B.c2 - A.c6 * B.c5 + A.c7 * B.c4;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c0 - A.c1 * B.c4 + A.c2 * B.c1 - A.c2 * B.c5 +
                     A.c3 * B.c7 - A.c6 * B.c2 - A.c7 * B.c7;
    ctype const c3 = A.c0 * B.c3 - A.c3 * B.c6 - A.c4 * B.c0 - A.c4 * B.c4 - A.c5 * B.c1 -
                     A.c5 * B.c5 + A.c6 * B.c3 - A.c7 * B.c6;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c3 - A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c2 +
                     A.c5 * B.c7 - A.c6 * B.c0 - A.c7 * B.c1;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c3 - A.c3 * B.c4 - A.c4 * B.c7 +
                     A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c3 * B.c3 - A.c4 * B.c1 - A.c4 * B.c5 + A.c5 * B.c0 +
                     A.c5 * B.c4 + A.c6 * B.c6 + A.c7 * B.c3;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c1 + A.c1 * B.c5 + A.c2 * B.c0 - A.c2 * B.c4 -
                     A.c3 * B.c2 - A.c6 * B.c7 + A.c7 * B.c2;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& A,
                                                        MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c3 + A.c2 * B.c4 - A.c3 * B.c1 - A.c4 * B.c6 +
                     A.c5 * B.c7 - A.c6 * B.c2 - A.c7 * B.c5;
    ctype const c1 = A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c5 - A.c3 * B.c2 - A.c4 * B.c7 -
                     A.c5 * B.c6 + A.c6 * B.c1 + A.c7 * B.c4;
    ctype const c2 = -A.c0 * B.c1 - A.c1 * B.c2 + A.c2 * B.c0 + A.c2 * B.c6 -
                     A.c4 * B.c1 - A.c5 * B.c2 + A.c7 * B.c3 + A.c7 * B.c7;
    ctype const c3 = A.c0 * B.c4 + A.c1 * B.c5 + A.c3 * B.c0 - A.c3 * B.c6 - A.c4 * B.c4 -
                     A.c5 * B.c5 - A.c6 * B.c3 + A.c6 * B.c7;
    ctype const c4 = -A.c0 * B.c6 + A.c1 * B.c7 + A.c2 * B.c4 + A.c3 * B.c1 +
                     A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c5;
    ctype const c5 = -A.c0 * B.c7 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c2 +
                     A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1 + A.c7 * B.c4;
    ctype const c6 = A.c0 * B.c5 - A.c1 * B.c4 + A.c3 * B.c3 - A.c3 * B.c7 - A.c4 * B.c5 +
                     A.c5 * B.c4 + A.c6 * B.c0 - A.c6 * B.c6;
    ctype const c7 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c3 - A.c2 * B.c7 + A.c4 * B.c2 -
                     A.c5 * B.c1 + A.c7 * B.c0 + A.c7 * B.c6;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& A,
                                                        PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c7 * ctype(ps);
    ctype const c1 = -A.c2 * ctype(ps);
    ctype const c2 = A.c1 * ctype(ps);
    ctype const c3 = A.c6 * ctype(ps);
    ctype const c4 = A.c5 * ctype(ps);
    ctype const c5 = -A.c4 * ctype(ps);
    ctype const c6 = -A.c3 * ctype(ps);
    ctype const c7 = A.c0 * ctype(ps);
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps,
                                                        MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c7;
    ctype const c1 = -ctype(ps) * B.c2;
    ctype const c2 = ctype(ps) * B.c1;
    ctype const c3 = ctype(ps) * B.c6;
    ctype const c4 = ctype(ps) * B.c5;
    ctype const c5 = -ctype(ps) * B.c4;
    ctype const c6 = -ctype(ps) * B.c3;
    ctype const c7 = ctype(ps) * B.c0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& M,
                                                        TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c2 * t.z - M.c5 * t.w - M.c6 * t.x - M.c7 * t.y;
    ctype const c1 = M.c1 * t.z + M.c4 * t.w - M.c6 * t.y + M.c7 * t.x;
    ctype const c2 = -M.c1 * t.x - M.c2 * t.y + M.c3 * t.w - M.c7 * t.w;
    ctype const c3 = -M.c3 * t.z - M.c4 * t.x - M.c5 * t.y - M.c7 * t.z;
    ctype const c4 = M.c0 * t.x - M.c2 * t.z + M.c3 * t.y + M.c5 * t.w;
    ctype const c5 = M.c0 * t.y + M.c1 * t.z - M.c3 * t.x - M.c4 * t.w;
    ctype const c6 = M.c0 * t.z - M.c4 * t.y + M.c5 * t.x + M.c6 * t.z;
    ctype const c7 = M.c0 * t.w + M.c1 * t.y - M.c2 * t.x - M.c6 * t.w;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t,
                                                        MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c6 + t.y * M.c7 - t.z * M.c2 - t.w * M.c5;
    ctype const c1 = -t.x * M.c7 - t.y * M.c6 + t.z * M.c1 + t.w * M.c4;
    ctype const c2 = -t.x * M.c1 - t.y * M.c2 + t.w * M.c3 + t.w * M.c7;
    ctype const c3 = -t.x * M.c4 - t.y * M.c5 - t.z * M.c3 + t.z * M.c7;
    ctype const c4 = t.x * M.c0 - t.y * M.c3 + t.z * M.c2 - t.w * M.c5;
    ctype const c5 = t.x * M.c3 + t.y * M.c0 - t.z * M.c1 + t.w * M.c4;
    ctype const c6 = -t.x * M.c5 + t.y * M.c4 + t.z * M.c0 - t.z * M.c6;
    ctype const c7 = t.x * M.c2 - t.y * M.c1 + t.w * M.c0 + t.w * M.c6;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& M,
                                                        BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.vz - M.c4 * B.vx -
                     M.c5 * B.vy + M.c6 * B.mz;
    ctype const c1 =
        M.c0 * B.vx + M.c1 * B.mz - M.c2 * B.vz + M.c3 * B.vy - M.c6 * B.vx - M.c7 * B.vy;
    ctype const c2 =
        M.c0 * B.vy + M.c1 * B.vz + M.c2 * B.mz - M.c3 * B.vx - M.c6 * B.vy + M.c7 * B.vx;
    ctype const c3 =
        M.c0 * B.vz - M.c1 * B.my + M.c2 * B.mx - M.c4 * B.vy + M.c5 * B.vx + M.c7 * B.mz;
    ctype const c4 =
        M.c0 * B.mx + M.c3 * B.my - M.c4 * B.mz - M.c5 * B.vz + M.c6 * B.mx + M.c7 * B.my;
    ctype const c5 =
        M.c0 * B.my - M.c3 * B.mx + M.c4 * B.vz - M.c5 * B.mz + M.c6 * B.my - M.c7 * B.mx;
    ctype const c6 =
        M.c0 * B.mz + M.c1 * B.mx + M.c2 * B.my - M.c4 * B.vx - M.c5 * B.vy - M.c7 * B.vz;
    ctype const c7 =
        M.c1 * B.my - M.c2 * B.mx + M.c3 * B.mz - M.c4 * B.vy + M.c5 * B.vx + M.c6 * B.vz;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B,
                                                        MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c3 - B.mx * M.c1 -
                     B.my * M.c2 + B.mz * M.c6;
    ctype const c1 =
        B.vx * M.c0 + B.vx * M.c6 - B.vy * M.c3 - B.vy * M.c7 + B.vz * M.c2 - B.mz * M.c1;
    ctype const c2 =
        B.vx * M.c3 + B.vx * M.c7 + B.vy * M.c0 + B.vy * M.c6 - B.vz * M.c1 - B.mz * M.c2;
    ctype const c3 = -B.vx * M.c5 + B.vy * M.c4 + B.vz * M.c0 - B.mx * M.c2 +
                     B.my * M.c1 + B.mz * M.c7;
    ctype const c4 =
        B.vz * M.c5 + B.mx * M.c0 - B.mx * M.c6 - B.my * M.c3 + B.my * M.c7 + B.mz * M.c4;
    ctype const c5 = -B.vz * M.c4 + B.mx * M.c3 - B.mx * M.c7 + B.my * M.c0 -
                     B.my * M.c6 + B.mz * M.c5;
    ctype const c6 =
        B.vx * M.c4 + B.vy * M.c5 - B.vz * M.c7 - B.mx * M.c1 - B.my * M.c2 + B.mz * M.c0;
    ctype const c7 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c6 - B.mx * M.c2 + B.my * M.c1 + B.mz * M.c3;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& A,
                                                        Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * v.x + A.c1 * v.w + A.c3 * v.y - A.c4 * v.z;
    ctype const c1 = A.c0 * v.y + A.c2 * v.w - A.c3 * v.x - A.c5 * v.z;
    ctype const c2 = A.c0 * v.z + A.c1 * v.x + A.c2 * v.y - A.c6 * v.z;
    ctype const c3 = A.c0 * v.w - A.c4 * v.x - A.c5 * v.y + A.c6 * v.w;
    ctype const c4 = A.c1 * v.w + A.c4 * v.z - A.c6 * v.x - A.c7 * v.y;
    ctype const c5 = A.c2 * v.w + A.c5 * v.z - A.c6 * v.y + A.c7 * v.x;
    ctype const c6 = A.c3 * v.w - A.c4 * v.y + A.c5 * v.x + A.c7 * v.w;
    ctype const c7 = -A.c1 * v.y + A.c2 * v.x - A.c3 * v.z + A.c7 * v.z;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v,
                                                        MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * B.c0 - v.y * B.c3 + v.z * B.c4 - v.w * B.c1;
    ctype const c1 = v.x * B.c3 + v.y * B.c0 + v.z * B.c5 - v.w * B.c2;
    ctype const c2 = -v.x * B.c1 - v.y * B.c2 + v.z * B.c0 + v.z * B.c6;
    ctype const c3 = v.x * B.c4 + v.y * B.c5 + v.w * B.c0 - v.w * B.c6;
    ctype const c4 = -v.x * B.c6 + v.y * B.c7 + v.z * B.c4 + v.w * B.c1;
    ctype const c5 = -v.x * B.c7 - v.y * B.c6 + v.z * B.c5 + v.w * B.c2;
    ctype const c6 = v.x * B.c5 - v.y * B.c4 + v.w * B.c3 - v.w * B.c7;
    ctype const c7 = v.x * B.c2 - v.y * B.c1 - v.z * B.c3 - v.z * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(MVec2dc_E<T> const& M,
                                                        Scalar2dc<U> s)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(Scalar2dc<T> s,
                                                        MVec2dc_E<U> const& M)
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
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& A,
                                                        MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 - A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c4 +
                     A.c5 * B.c5 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c1 = -A.c0 * B.c2 - A.c1 * B.c7 + A.c2 * B.c0 - A.c2 * B.c4 -
                     A.c4 * B.c2 - A.c5 * B.c7 - A.c7 * B.c1 + A.c7 * B.c5;
    ctype const c2 = A.c0 * B.c7 - A.c1 * B.c2 + A.c2 * B.c1 - A.c2 * B.c5 + A.c4 * B.c7 -
                     A.c5 * B.c2 + A.c7 * B.c0 - A.c7 * B.c4;
    ctype const c3 = A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c7 + A.c4 * B.c5 -
                     A.c5 * B.c4 - A.c6 * B.c2 + A.c7 * B.c3;
    ctype const c4 = A.c0 * B.c3 - A.c1 * B.c6 - A.c3 * B.c0 - A.c3 * B.c4 - A.c4 * B.c3 +
                     A.c5 * B.c6 - A.c6 * B.c1 - A.c6 * B.c5;
    ctype const c5 = A.c0 * B.c6 + A.c1 * B.c3 - A.c3 * B.c1 - A.c3 * B.c5 - A.c4 * B.c6 -
                     A.c5 * B.c3 + A.c6 * B.c0 + A.c6 * B.c4;
    ctype const c6 = -A.c0 * B.c4 - A.c1 * B.c5 + A.c2 * B.c3 - A.c3 * B.c2 -
                     A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c6 + A.c3 * B.c7 -
                     A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2 - A.c7 * B.c3;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& A,
                                                        PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c5 * ctype(ps);
    ctype const c1 = -A.c4 * ctype(ps);
    ctype const c2 = A.c7 * ctype(ps);
    ctype const c3 = A.c6 * ctype(ps);
    ctype const c4 = A.c1 * ctype(ps);
    ctype const c5 = -A.c0 * ctype(ps);
    ctype const c6 = -A.c3 * ctype(ps);
    ctype const c7 = -A.c2 * ctype(ps);
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps,
                                                        MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c5;
    ctype const c1 = ctype(ps) * B.c4;
    ctype const c2 = -ctype(ps) * B.c7;
    ctype const c3 = -ctype(ps) * B.c6;
    ctype const c4 = -ctype(ps) * B.c1;
    ctype const c5 = ctype(ps) * B.c0;
    ctype const c6 = ctype(ps) * B.c3;
    ctype const c7 = ctype(ps) * B.c2;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& M,
                                                        TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c4 * t.x + M.c5 * t.y - M.c6 * t.w - M.c7 * t.z;
    ctype const c1 = -M.c1 * t.w - M.c2 * t.x - M.c5 * t.w + M.c7 * t.y;
    ctype const c2 = M.c0 * t.w - M.c2 * t.y + M.c4 * t.w - M.c7 * t.x;
    ctype const c3 = -M.c2 * t.z + M.c3 * t.w + M.c4 * t.y - M.c5 * t.x;
    ctype const c4 = -M.c1 * t.z - M.c3 * t.x + M.c5 * t.z - M.c6 * t.y;
    ctype const c5 = M.c0 * t.z - M.c3 * t.y - M.c4 * t.z + M.c6 * t.x;
    ctype const c6 = -M.c0 * t.x - M.c1 * t.y - M.c6 * t.w + M.c7 * t.z;
    ctype const c7 = -M.c0 * t.y + M.c1 * t.x + M.c2 * t.z + M.c3 * t.w;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t,
                                                        MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * M.c4 + t.y * M.c5 - t.z * M.c7 - t.w * M.c6;
    ctype const c1 = -t.x * M.c2 - t.y * M.c7 - t.w * M.c1 + t.w * M.c5;
    ctype const c2 = t.x * M.c7 - t.y * M.c2 + t.w * M.c0 - t.w * M.c4;
    ctype const c3 = t.x * M.c5 - t.y * M.c4 - t.z * M.c2 + t.w * M.c3;
    ctype const c4 = -t.x * M.c3 + t.y * M.c6 - t.z * M.c1 - t.z * M.c5;
    ctype const c5 = -t.x * M.c6 - t.y * M.c3 + t.z * M.c0 + t.z * M.c4;
    ctype const c6 = -t.x * M.c0 - t.y * M.c1 - t.z * M.c7 + t.w * M.c6;
    ctype const c7 = -t.x * M.c1 + t.y * M.c0 - t.z * M.c2 - t.w * M.c3;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& M,
                                                        BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.vz + M.c2 * B.mx - M.c3 * B.vx - M.c4 * B.mz -
                     M.c6 * B.vy - M.c7 * B.my;
    ctype const c1 =
        M.c0 * B.vz + M.c2 * B.my - M.c3 * B.vy - M.c5 * B.mz + M.c6 * B.vx + M.c7 * B.mx;
    ctype const c2 = -M.c0 * B.vx - M.c1 * B.vy + M.c2 * B.mz - M.c4 * B.vx -
                     M.c5 * B.vy + M.c7 * B.vz;
    ctype const c3 =
        M.c0 * B.mx + M.c1 * B.my - M.c3 * B.mz - M.c4 * B.mx - M.c5 * B.my - M.c6 * B.vz;
    ctype const c4 = -M.c0 * B.mz + M.c2 * B.mx + M.c3 * B.vx - M.c5 * B.vz +
                     M.c6 * B.vy - M.c7 * B.my;
    ctype const c5 = -M.c1 * B.mz + M.c2 * B.my + M.c3 * B.vy + M.c4 * B.vz -
                     M.c6 * B.vx + M.c7 * B.mx;
    ctype const c6 =
        M.c0 * B.my - M.c1 * B.mx + M.c3 * B.vz - M.c4 * B.my + M.c5 * B.mx - M.c6 * B.mz;
    ctype const c7 =
        M.c0 * B.vy - M.c1 * B.vx - M.c2 * B.vz + M.c4 * B.vy - M.c5 * B.vx + M.c7 * B.mz;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B,
                                                        MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c3 - B.vy * M.c6 + B.vz * M.c1 - B.mx * M.c2 - B.my * M.c7 - B.mz * M.c4;
    ctype const c1 =
        B.vx * M.c6 + B.vy * M.c3 - B.vz * M.c0 + B.mx * M.c7 - B.my * M.c2 - B.mz * M.c5;
    ctype const c2 =
        B.vx * M.c0 - B.vx * M.c4 + B.vy * M.c1 - B.vy * M.c5 + B.vz * M.c7 - B.mz * M.c2;
    ctype const c3 = -B.vz * M.c6 - B.mx * M.c0 - B.mx * M.c4 - B.my * M.c1 -
                     B.my * M.c5 + B.mz * M.c3;
    ctype const c4 =
        B.vx * M.c3 - B.vy * M.c6 + B.vz * M.c5 + B.mx * M.c2 + B.my * M.c7 - B.mz * M.c0;
    ctype const c5 =
        B.vx * M.c6 + B.vy * M.c3 - B.vz * M.c4 - B.mx * M.c7 + B.my * M.c2 - B.mz * M.c1;
    ctype const c6 =
        B.vz * M.c3 - B.mx * M.c1 - B.mx * M.c5 + B.my * M.c0 + B.my * M.c4 + B.mz * M.c6;
    ctype const c7 = -B.vx * M.c1 + B.vx * M.c5 + B.vy * M.c0 - B.vy * M.c4 -
                     B.vz * M.c2 - B.mz * M.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& M,
                                                        Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y - M.c2 * v.w - M.c3 * v.z;
    ctype const c1 = -M.c0 * v.z + M.c2 * v.x - M.c4 * v.z - M.c7 * v.y;
    ctype const c2 = -M.c1 * v.z + M.c2 * v.y - M.c5 * v.z + M.c7 * v.x;
    ctype const c3 = M.c0 * v.y - M.c1 * v.x - M.c6 * v.z + M.c7 * v.w;
    ctype const c4 = M.c0 * v.w - M.c3 * v.x - M.c4 * v.w - M.c6 * v.y;
    ctype const c5 = M.c1 * v.w - M.c3 * v.y - M.c5 * v.w + M.c6 * v.x;
    ctype const c6 = M.c2 * v.w - M.c3 * v.z - M.c4 * v.x - M.c5 * v.y;
    ctype const c7 = -M.c4 * v.y + M.c5 * v.x - M.c6 * v.z - M.c7 * v.w;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v,
                                                        MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 - v.z * M.c3 - v.w * M.c2;
    ctype const c1 = -v.x * M.c2 - v.y * M.c7 + v.z * M.c0 - v.z * M.c4;
    ctype const c2 = v.x * M.c7 - v.y * M.c2 + v.z * M.c1 - v.z * M.c5;
    ctype const c3 = v.x * M.c1 - v.y * M.c0 - v.z * M.c6 + v.w * M.c7;
    ctype const c4 = v.x * M.c3 - v.y * M.c6 - v.w * M.c0 - v.w * M.c4;
    ctype const c5 = v.x * M.c6 + v.y * M.c3 - v.w * M.c1 - v.w * M.c5;
    ctype const c6 = -v.x * M.c4 - v.y * M.c5 + v.z * M.c3 - v.w * M.c2;
    ctype const c7 = -v.x * M.c5 + v.y * M.c4 + v.z * M.c6 + v.w * M.c7;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(MVec2dc_U<T> const& M,
                                                        Scalar2dc<U> s)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(Scalar2dc<T> s,
                                                        MVec2dc_U<U> const& M)
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
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps1,
                                                        PScalar2dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(-ctype(ps1) * ctype(ps2));
}

// cga2dc gpr :: gpr(ps,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps,
                                                     TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * t.y;
    ctype const c1 = ctype(ps) * t.x;
    ctype const c2 = -ctype(ps) * t.w;
    ctype const c3 = -ctype(ps) * t.z;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(trivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t,
                                                     PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * ctype(ps);
    ctype const c1 = -t.x * ctype(ps);
    ctype const c2 = t.w * ctype(ps);
    ctype const c3 = t.z * ctype(ps);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps,
                                                       BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.vy;
    ctype const c1 = ctype(ps) * B.vx;
    ctype const c2 = ctype(ps) * B.mz;
    ctype const c3 = ctype(ps) * B.my;
    ctype const c4 = -ctype(ps) * B.mx;
    ctype const c5 = -ctype(ps) * B.vz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc gpr :: gpr(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B,
                                                       PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * ctype(ps);
    ctype const c1 = B.vx * ctype(ps);
    ctype const c2 = B.mz * ctype(ps);
    ctype const c3 = B.my * ctype(ps);
    ctype const c4 = -B.mx * ctype(ps);
    ctype const c5 = -B.vz * ctype(ps);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc gpr :: gpr(ps,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps,
                                                        Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.y;
    ctype const c1 = ctype(ps) * v.x;
    ctype const c2 = ctype(ps) * v.w;
    ctype const c3 = ctype(ps) * v.z;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(vec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v,
                                                        PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * ctype(ps);
    ctype const c1 = -v.x * ctype(ps);
    ctype const c2 = -v.w * ctype(ps);
    ctype const c3 = -v.z * ctype(ps);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> operator*(PScalar2dc<T> ps, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps) * ctype(s));
}

// cga2dc gpr :: gpr(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> operator*(Scalar2dc<T> s, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(s) * ctype(ps));
}

// cga2dc gpr :: gpr(trivec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t1,
                                                        TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.x * t2.x + t1.y * t2.y - t1.z * t2.w - t1.w * t2.z;
    ctype const c1 = -t1.y * t2.w + t1.w * t2.y;
    ctype const c2 = t1.x * t2.w - t1.w * t2.x;
    ctype const c3 = t1.x * t2.y - t1.y * t2.x;
    ctype const c4 = t1.y * t2.z - t1.z * t2.y;
    ctype const c5 = -t1.x * t2.z + t1.z * t2.x;
    ctype const c6 = -t1.z * t2.w + t1.w * t2.z;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(trivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t,
                                                        BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * B.mz - t.z * B.vy - t.w * B.my;
    ctype const c1 = -t.y * B.mz + t.z * B.vx + t.w * B.mx;
    ctype const c2 = -t.x * B.vx - t.y * B.vy + t.w * B.vz;
    ctype const c3 = -t.x * B.mx - t.y * B.my - t.z * B.vz;
    ctype const c4 = -t.y * B.vz + t.z * B.vy - t.w * B.my;
    ctype const c5 = t.x * B.vz - t.z * B.vx + t.w * B.mx;
    ctype const c6 = -t.x * B.my + t.y * B.mx - t.z * B.mz;
    ctype const c7 = t.x * B.vy - t.y * B.vx + t.w * B.mz;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(bivec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B,
                                                        TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * t.z - B.my * t.w - B.mz * t.x;
    ctype const c1 = B.vx * t.z + B.mx * t.w - B.mz * t.y;
    ctype const c2 = -B.vx * t.x - B.vy * t.y + B.vz * t.w;
    ctype const c3 = -B.vz * t.z - B.mx * t.x - B.my * t.y;
    ctype const c4 = -B.vy * t.z + B.vz * t.y + B.my * t.w;
    ctype const c5 = B.vx * t.z - B.vz * t.x - B.mx * t.w;
    ctype const c6 = -B.mx * t.y + B.my * t.x + B.mz * t.z;
    ctype const c7 = B.vx * t.y - B.vy * t.x - B.mz * t.w;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(trivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t,
                                                        Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -t.x * v.z - t.w * v.y;
    ctype const c2 = -t.y * v.z + t.w * v.x;
    ctype const c3 = -t.z * v.z + t.w * v.w;
    ctype const c4 = -t.x * v.w - t.z * v.y;
    ctype const c5 = -t.y * v.w + t.z * v.x;
    ctype const c6 = -t.x * v.x - t.y * v.y;
    ctype const c7 = -t.x * v.y + t.y * v.x - t.z * v.z - t.w * v.w;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(vec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v,
                                                        TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -v.y * t.w - v.z * t.x;
    ctype const c2 = v.x * t.w - v.z * t.y;
    ctype const c3 = -v.z * t.z + v.w * t.w;
    ctype const c4 = -v.y * t.z - v.w * t.x;
    ctype const c5 = v.x * t.z - v.w * t.y;
    ctype const c6 = -v.x * t.x - v.y * t.y;
    ctype const c7 = -v.x * t.y + v.y * t.x + v.z * t.z + v.w * t.w;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator*(TriVec2dc<T> const& t,
                                                        Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> operator*(Scalar2dc<T> s,
                                                        TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(bivec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B1,
                                                        BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.vz - B1.mx * B2.vx -
                     B1.my * B2.vy + B1.mz * B2.mz;
    ctype const c1 = B1.vx * B2.mz - B1.vy * B2.vz + B1.vz * B2.vy - B1.mz * B2.vx;
    ctype const c2 = B1.vx * B2.vz + B1.vy * B2.mz - B1.vz * B2.vx - B1.mz * B2.vy;
    ctype const c3 = -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx;
    ctype const c4 = B1.vz * B2.my - B1.mx * B2.mz - B1.my * B2.vz + B1.mz * B2.mx;
    ctype const c5 = -B1.vz * B2.mx + B1.mx * B2.vz - B1.my * B2.mz + B1.mz * B2.my;
    ctype const c6 = B1.vx * B2.mx + B1.vy * B2.my - B1.mx * B2.vx - B1.my * B2.vy;
    ctype const c7 = B1.vx * B2.my - B1.vy * B2.mx + B1.vz * B2.mz - B1.mx * B2.vy +
                     B1.my * B2.vx + B1.mz * B2.vz;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(bivec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B,
                                                        Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.w + B.vz * v.y - B.mx * v.z;
    ctype const c1 = B.vy * v.w - B.vz * v.x - B.my * v.z;
    ctype const c2 = B.vx * v.x + B.vy * v.y - B.mz * v.z;
    ctype const c3 = -B.mx * v.x - B.my * v.y + B.mz * v.w;
    ctype const c4 = B.vx * v.w + B.mx * v.z - B.mz * v.x;
    ctype const c5 = B.vy * v.w + B.my * v.z - B.mz * v.y;
    ctype const c6 = B.vz * v.w - B.mx * v.y + B.my * v.x;
    ctype const c7 = -B.vx * v.y + B.vy * v.x - B.vz * v.z;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(vec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v,
                                                        BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.vz + v.z * B.mx - v.w * B.vx;
    ctype const c1 = v.x * B.vz + v.z * B.my - v.w * B.vy;
    ctype const c2 = -v.x * B.vx - v.y * B.vy + v.z * B.mz;
    ctype const c3 = v.x * B.mx + v.y * B.my - v.w * B.mz;
    ctype const c4 = -v.x * B.mz + v.z * B.mx + v.w * B.vx;
    ctype const c5 = -v.y * B.mz + v.z * B.my + v.w * B.vy;
    ctype const c6 = v.x * B.my - v.y * B.mx + v.w * B.vz;
    ctype const c7 = v.x * B.vy - v.y * B.vx - v.z * B.vz;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc gpr :: gpr(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator*(BiVec2dc<T> const& B,
                                                       Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc gpr :: gpr(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> operator*(Scalar2dc<T> s,
                                                       BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc gpr :: gpr(vec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v1,
                                                        Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.x * v2.x + v1.y * v2.y - v1.z * v2.w - v1.w * v2.z;
    ctype const c1 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c2 = -v1.y * v2.z + v1.z * v2.y;
    ctype const c3 = v1.x * v2.y - v1.y * v2.x;
    ctype const c4 = v1.x * v2.w - v1.w * v2.x;
    ctype const c5 = v1.y * v2.w - v1.w * v2.y;
    ctype const c6 = v1.z * v2.w - v1.w * v2.z;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc gpr :: gpr(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator*(Vec2dc<T> const& v, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> operator*(Scalar2dc<T> s, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc gpr :: gpr(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> operator*(Scalar2dc<T> s1, Scalar2dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s1) * ctype(s2));
}
////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
//
// rgpr(ul, ur) = l_cmpl(gpr(r_cmpl(ul),r_cmpl(ur)))
//
// The identity element of rgpr is the pseudoscalar: rgpr(X, I_2dc) = X.
////////////////////////////////////////////////////////////////////////////////

// cga2dc rgpr :: rgpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& A, MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c1 * B.c12 + A.c2 * B.c11 + A.c3 * B.c13 +
                     A.c4 * B.c14 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c10 -
                     A.c8 * B.c6 + A.c9 * B.c5 + A.c10 * B.c7 - A.c11 * B.c2 +
                     A.c12 * B.c1 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    ctype const c1 = -A.c0 * B.c12 + A.c1 * B.c15 - A.c2 * B.c10 + A.c3 * B.c9 +
                     A.c4 * B.c6 + A.c5 * B.c13 - A.c6 * B.c4 - A.c7 * B.c11 -
                     A.c8 * B.c14 - A.c9 * B.c3 + A.c10 * B.c2 - A.c11 * B.c7 +
                     A.c12 * B.c0 + A.c13 * B.c5 - A.c14 * B.c8 + A.c15 * B.c1;
    ctype const c2 = A.c0 * B.c11 + A.c1 * B.c10 + A.c2 * B.c15 - A.c3 * B.c8 -
                     A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c13 - A.c7 * B.c12 +
                     A.c8 * B.c3 - A.c9 * B.c14 - A.c10 * B.c1 - A.c11 * B.c0 -
                     A.c12 * B.c7 + A.c13 * B.c6 - A.c14 * B.c9 + A.c15 * B.c2;
    ctype const c3 = -A.c0 * B.c14 + A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c7 +
                     A.c3 * B.c15 + A.c5 * B.c2 + A.c5 * B.c12 - A.c6 * B.c1 -
                     A.c6 * B.c11 + A.c7 * B.c3 - A.c10 * B.c14 - A.c11 * B.c6 +
                     A.c12 * B.c5 + A.c14 * B.c0 - A.c14 * B.c10 + A.c15 * B.c3;
    ctype const c4 = -A.c0 * B.c13 + A.c1 * B.c9 - A.c2 * B.c8 + A.c4 * B.c7 +
                     A.c4 * B.c15 - A.c7 * B.c4 + A.c8 * B.c2 - A.c8 * B.c12 -
                     A.c9 * B.c1 + A.c9 * B.c11 + A.c10 * B.c13 + A.c11 * B.c9 -
                     A.c12 * B.c8 + A.c13 * B.c0 + A.c13 * B.c10 + A.c15 * B.c4;
    ctype const c5 = A.c0 * B.c6 - A.c1 * B.c14 - A.c2 * B.c3 + A.c3 * B.c2 +
                     A.c3 * B.c12 - A.c5 * B.c7 + A.c5 * B.c15 + A.c6 * B.c0 -
                     A.c6 * B.c10 + A.c7 * B.c5 + A.c10 * B.c6 + A.c11 * B.c14 +
                     A.c12 * B.c3 - A.c14 * B.c1 - A.c14 * B.c11 + A.c15 * B.c5;
    ctype const c6 = -A.c0 * B.c5 + A.c1 * B.c3 - A.c2 * B.c14 - A.c3 * B.c1 -
                     A.c3 * B.c11 - A.c5 * B.c0 + A.c5 * B.c10 - A.c6 * B.c7 +
                     A.c6 * B.c15 + A.c7 * B.c6 - A.c10 * B.c5 - A.c11 * B.c3 +
                     A.c12 * B.c14 - A.c14 * B.c2 - A.c14 * B.c12 + A.c15 * B.c6;
    ctype const c7 = -A.c0 * B.c10 - A.c1 * B.c11 - A.c2 * B.c12 - A.c3 * B.c4 +
                     A.c4 * B.c3 + A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c15 -
                     A.c8 * B.c5 - A.c9 * B.c6 - A.c10 * B.c0 - A.c11 * B.c1 -
                     A.c12 * B.c2 + A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c7;
    ctype const c8 = -A.c0 * B.c9 + A.c1 * B.c13 - A.c2 * B.c4 + A.c4 * B.c2 -
                     A.c4 * B.c12 - A.c7 * B.c8 + A.c8 * B.c7 + A.c8 * B.c15 -
                     A.c9 * B.c0 - A.c9 * B.c10 + A.c10 * B.c9 + A.c11 * B.c13 -
                     A.c12 * B.c4 + A.c13 * B.c1 - A.c13 * B.c11 + A.c15 * B.c8;
    ctype const c9 = A.c0 * B.c8 + A.c1 * B.c4 + A.c2 * B.c13 - A.c4 * B.c1 +
                     A.c4 * B.c11 - A.c7 * B.c9 + A.c8 * B.c0 + A.c8 * B.c10 +
                     A.c9 * B.c7 + A.c9 * B.c15 - A.c10 * B.c8 + A.c11 * B.c4 +
                     A.c12 * B.c13 + A.c13 * B.c2 - A.c13 * B.c12 + A.c15 * B.c9;
    ctype const c10 = A.c0 * B.c7 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c13 -
                      A.c4 * B.c14 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c0 +
                      A.c8 * B.c6 - A.c9 * B.c5 + A.c10 * B.c15 + A.c11 * B.c12 -
                      A.c12 * B.c11 + A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c10;
    ctype const c11 = -A.c0 * B.c2 - A.c1 * B.c7 + A.c2 * B.c0 - A.c3 * B.c9 +
                      A.c4 * B.c6 - A.c5 * B.c13 + A.c6 * B.c4 - A.c7 * B.c1 -
                      A.c8 * B.c14 - A.c9 * B.c3 + A.c10 * B.c12 + A.c11 * B.c15 -
                      A.c12 * B.c10 + A.c13 * B.c5 + A.c14 * B.c8 + A.c15 * B.c11;
    ctype const c12 = A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c7 + A.c3 * B.c8 -
                      A.c4 * B.c5 - A.c5 * B.c4 - A.c6 * B.c13 - A.c7 * B.c2 +
                      A.c8 * B.c3 - A.c9 * B.c14 - A.c10 * B.c11 + A.c11 * B.c10 +
                      A.c12 * B.c15 + A.c13 * B.c6 + A.c14 * B.c9 + A.c15 * B.c12;
    ctype const c13 = A.c0 * B.c4 + A.c1 * B.c8 + A.c2 * B.c9 - A.c4 * B.c0 -
                      A.c4 * B.c10 - A.c7 * B.c13 + A.c8 * B.c1 - A.c8 * B.c11 +
                      A.c9 * B.c2 - A.c9 * B.c12 - A.c10 * B.c4 + A.c11 * B.c8 +
                      A.c12 * B.c9 + A.c13 * B.c7 + A.c13 * B.c15 + A.c15 * B.c13;
    ctype const c14 = A.c0 * B.c3 - A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c0 +
                      A.c3 * B.c10 - A.c5 * B.c1 - A.c5 * B.c11 - A.c6 * B.c2 -
                      A.c6 * B.c12 + A.c7 * B.c14 + A.c10 * B.c3 + A.c11 * B.c5 +
                      A.c12 * B.c6 - A.c14 * B.c7 + A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = -A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c4 -
                      A.c4 * B.c3 + A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c7 +
                      A.c8 * B.c5 + A.c9 * B.c6 - A.c10 * B.c10 + A.c11 * B.c11 +
                      A.c12 * B.c12 - A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c15;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& A,
                                                 MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c6 - A.c8 * B.c2 +
                     A.c9 * B.c1 + A.c10 * B.c3 + A.c15 * B.c0;
    ctype const c1 = A.c1 * B.c7 - A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c2 -
                     A.c11 * B.c3 + A.c12 * B.c0 + A.c13 * B.c1 - A.c14 * B.c4;
    ctype const c2 = A.c1 * B.c6 + A.c2 * B.c7 - A.c3 * B.c4 - A.c4 * B.c1 -
                     A.c11 * B.c0 - A.c12 * B.c3 + A.c13 * B.c2 - A.c14 * B.c5;
    ctype const c3 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c3 + A.c3 * B.c7 -
                     A.c11 * B.c2 + A.c12 * B.c1 + A.c14 * B.c0 - A.c14 * B.c6;
    ctype const c4 = A.c1 * B.c5 - A.c2 * B.c4 + A.c4 * B.c3 + A.c4 * B.c7 +
                     A.c11 * B.c5 - A.c12 * B.c4 + A.c13 * B.c0 + A.c13 * B.c6;
    ctype const c5 = A.c0 * B.c2 - A.c5 * B.c3 + A.c5 * B.c7 + A.c6 * B.c0 - A.c6 * B.c6 +
                     A.c7 * B.c1 + A.c10 * B.c2 + A.c15 * B.c1;
    ctype const c6 = -A.c0 * B.c1 - A.c5 * B.c0 + A.c5 * B.c6 - A.c6 * B.c3 +
                     A.c6 * B.c7 + A.c7 * B.c2 - A.c10 * B.c1 + A.c15 * B.c2;
    ctype const c7 = -A.c0 * B.c6 + A.c5 * B.c4 + A.c6 * B.c5 + A.c7 * B.c7 -
                     A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c0 + A.c15 * B.c3;
    ctype const c8 = -A.c0 * B.c5 - A.c7 * B.c4 + A.c8 * B.c3 + A.c8 * B.c7 -
                     A.c9 * B.c0 - A.c9 * B.c6 + A.c10 * B.c5 + A.c15 * B.c4;
    ctype const c9 = A.c0 * B.c4 - A.c7 * B.c5 + A.c8 * B.c0 + A.c8 * B.c6 + A.c9 * B.c3 +
                     A.c9 * B.c7 - A.c10 * B.c4 + A.c15 * B.c5;
    ctype const c10 = A.c0 * B.c3 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c0 +
                      A.c8 * B.c2 - A.c9 * B.c1 + A.c10 * B.c7 + A.c15 * B.c6;
    ctype const c11 = -A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c5 + A.c4 * B.c2 +
                      A.c11 * B.c7 - A.c12 * B.c6 + A.c13 * B.c1 + A.c14 * B.c4;
    ctype const c12 = -A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c4 - A.c4 * B.c1 +
                      A.c11 * B.c6 + A.c12 * B.c7 + A.c13 * B.c2 + A.c14 * B.c5;
    ctype const c13 = A.c1 * B.c4 + A.c2 * B.c5 - A.c4 * B.c0 - A.c4 * B.c6 +
                      A.c11 * B.c4 + A.c12 * B.c5 + A.c13 * B.c3 + A.c13 * B.c7;
    ctype const c14 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c0 + A.c3 * B.c6 +
                      A.c11 * B.c1 + A.c12 * B.c2 - A.c14 * B.c3 + A.c14 * B.c7;
    ctype const c15 = -A.c0 * B.c0 + A.c5 * B.c4 + A.c6 * B.c5 + A.c7 * B.c3 +
                      A.c8 * B.c1 + A.c9 * B.c2 - A.c10 * B.c6 + A.c15 * B.c7;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& A,
                                                 MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c10 -
                     A.c4 * B.c6 + A.c5 * B.c5 + A.c6 * B.c7 + A.c7 * B.c0;
    ctype const c1 = -A.c0 * B.c12 + A.c1 * B.c13 - A.c2 * B.c4 - A.c3 * B.c11 -
                     A.c4 * B.c14 - A.c5 * B.c3 + A.c6 * B.c2 + A.c7 * B.c1;
    ctype const c2 = A.c0 * B.c11 + A.c1 * B.c4 + A.c2 * B.c13 - A.c3 * B.c12 +
                     A.c4 * B.c3 - A.c5 * B.c14 - A.c6 * B.c1 + A.c7 * B.c2;
    ctype const c3 = -A.c0 * B.c14 + A.c1 * B.c2 + A.c1 * B.c12 - A.c2 * B.c1 -
                     A.c2 * B.c11 + A.c3 * B.c3 - A.c6 * B.c14 + A.c7 * B.c3;
    ctype const c4 = -A.c0 * B.c13 - A.c3 * B.c4 + A.c4 * B.c2 - A.c4 * B.c12 -
                     A.c5 * B.c1 + A.c5 * B.c11 + A.c6 * B.c13 + A.c7 * B.c4;
    ctype const c5 = A.c0 * B.c6 - A.c1 * B.c7 + A.c1 * B.c15 + A.c2 * B.c0 -
                     A.c2 * B.c10 + A.c3 * B.c5 + A.c6 * B.c6 + A.c7 * B.c5;
    ctype const c6 = -A.c0 * B.c5 - A.c1 * B.c0 + A.c1 * B.c10 - A.c2 * B.c7 +
                     A.c2 * B.c15 + A.c3 * B.c6 - A.c6 * B.c5 + A.c7 * B.c6;
    ctype const c7 = -A.c0 * B.c10 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c15 -
                     A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c0 + A.c7 * B.c7;
    ctype const c8 = -A.c0 * B.c9 - A.c3 * B.c8 + A.c4 * B.c7 + A.c4 * B.c15 -
                     A.c5 * B.c0 - A.c5 * B.c10 + A.c6 * B.c9 + A.c7 * B.c8;
    ctype const c9 = A.c0 * B.c8 - A.c3 * B.c9 + A.c4 * B.c0 + A.c4 * B.c10 +
                     A.c5 * B.c7 + A.c5 * B.c15 - A.c6 * B.c8 + A.c7 * B.c9;
    ctype const c10 = A.c0 * B.c7 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c0 +
                      A.c4 * B.c6 - A.c5 * B.c5 + A.c6 * B.c15 + A.c7 * B.c10;
    ctype const c11 = -A.c0 * B.c2 - A.c1 * B.c13 + A.c2 * B.c4 - A.c3 * B.c1 -
                      A.c4 * B.c14 - A.c5 * B.c3 + A.c6 * B.c12 + A.c7 * B.c11;
    ctype const c12 = A.c0 * B.c1 - A.c1 * B.c4 - A.c2 * B.c13 - A.c3 * B.c2 +
                      A.c4 * B.c3 - A.c5 * B.c14 - A.c6 * B.c11 + A.c7 * B.c12;
    ctype const c13 = A.c0 * B.c4 - A.c3 * B.c13 + A.c4 * B.c1 - A.c4 * B.c11 +
                      A.c5 * B.c2 - A.c5 * B.c12 - A.c6 * B.c4 + A.c7 * B.c13;
    ctype const c14 = A.c0 * B.c3 - A.c1 * B.c1 - A.c1 * B.c11 - A.c2 * B.c2 -
                      A.c2 * B.c12 + A.c3 * B.c14 + A.c6 * B.c3 + A.c7 * B.c14;
    ctype const c15 = -A.c0 * B.c0 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c7 +
                      A.c4 * B.c5 + A.c5 * B.c6 - A.c6 * B.c10 + A.c7 * B.c15;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& A,
                                                 MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c6 + A.c4 * B.c7 -
                     A.c11 * B.c1 + A.c12 * B.c0 - A.c13 * B.c2 - A.c14 * B.c3;
    ctype const c1 = -A.c0 * B.c5 + A.c5 * B.c6 - A.c6 * B.c3 - A.c7 * B.c4 -
                     A.c8 * B.c7 - A.c9 * B.c2 + A.c10 * B.c1 + A.c15 * B.c0;
    ctype const c2 = A.c0 * B.c4 + A.c5 * B.c3 + A.c6 * B.c6 - A.c7 * B.c5 + A.c8 * B.c2 -
                     A.c9 * B.c7 - A.c10 * B.c0 + A.c15 * B.c1;
    ctype const c3 = -A.c0 * B.c7 + A.c5 * B.c1 + A.c5 * B.c5 - A.c6 * B.c0 -
                     A.c6 * B.c4 + A.c7 * B.c2 - A.c10 * B.c7 + A.c15 * B.c2;
    ctype const c4 = -A.c0 * B.c6 - A.c7 * B.c3 + A.c8 * B.c1 - A.c8 * B.c5 -
                     A.c9 * B.c0 + A.c9 * B.c4 + A.c10 * B.c6 + A.c15 * B.c3;
    ctype const c5 = -A.c1 * B.c7 - A.c2 * B.c2 + A.c3 * B.c1 + A.c3 * B.c5 +
                     A.c11 * B.c7 + A.c12 * B.c2 - A.c14 * B.c0 - A.c14 * B.c4;
    ctype const c6 = A.c1 * B.c2 - A.c2 * B.c7 - A.c3 * B.c0 - A.c3 * B.c4 -
                     A.c11 * B.c2 + A.c12 * B.c7 - A.c14 * B.c1 - A.c14 * B.c5;
    ctype const c7 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c3 + A.c4 * B.c2 -
                     A.c11 * B.c0 - A.c12 * B.c1 + A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c8 = A.c1 * B.c6 - A.c2 * B.c3 + A.c4 * B.c1 - A.c4 * B.c5 +
                     A.c11 * B.c6 - A.c12 * B.c3 + A.c13 * B.c0 - A.c13 * B.c4;
    ctype const c9 = A.c1 * B.c3 + A.c2 * B.c6 - A.c4 * B.c0 + A.c4 * B.c4 +
                     A.c11 * B.c3 + A.c12 * B.c6 + A.c13 * B.c1 - A.c13 * B.c5;
    ctype const c10 = A.c1 * B.c1 - A.c2 * B.c0 + A.c3 * B.c6 - A.c4 * B.c7 +
                      A.c11 * B.c5 - A.c12 * B.c4 + A.c13 * B.c2 - A.c14 * B.c3;
    ctype const c11 = -A.c0 * B.c1 - A.c5 * B.c6 + A.c6 * B.c3 - A.c7 * B.c0 -
                      A.c8 * B.c7 - A.c9 * B.c2 + A.c10 * B.c5 + A.c15 * B.c4;
    ctype const c12 = A.c0 * B.c0 - A.c5 * B.c3 - A.c6 * B.c6 - A.c7 * B.c1 +
                      A.c8 * B.c2 - A.c9 * B.c7 - A.c10 * B.c4 + A.c15 * B.c5;
    ctype const c13 = A.c0 * B.c3 - A.c7 * B.c6 + A.c8 * B.c0 - A.c8 * B.c4 +
                      A.c9 * B.c1 - A.c9 * B.c5 - A.c10 * B.c3 + A.c15 * B.c6;
    ctype const c14 = A.c0 * B.c2 - A.c5 * B.c0 - A.c5 * B.c4 - A.c6 * B.c1 -
                      A.c6 * B.c5 + A.c7 * B.c7 + A.c10 * B.c2 + A.c15 * B.c7;
    ctype const c15 = A.c1 * B.c0 + A.c2 * B.c1 - A.c3 * B.c3 - A.c4 * B.c2 +
                      A.c11 * B.c4 + A.c12 * B.c5 - A.c13 * B.c7 - A.c14 * B.c6;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& A,
                                                 MVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c12 + A.c1 * B.c11 + A.c2 * B.c13 + A.c3 * B.c14 -
                     A.c4 * B.c2 + A.c5 * B.c1 - A.c6 * B.c3 - A.c7 * B.c4;
    ctype const c1 = A.c0 * B.c15 - A.c1 * B.c10 + A.c2 * B.c9 + A.c3 * B.c6 -
                     A.c4 * B.c7 + A.c5 * B.c0 + A.c6 * B.c5 - A.c7 * B.c8;
    ctype const c2 = A.c0 * B.c10 + A.c1 * B.c15 - A.c2 * B.c8 - A.c3 * B.c5 -
                     A.c4 * B.c0 - A.c5 * B.c7 + A.c6 * B.c6 - A.c7 * B.c9;
    ctype const c3 = A.c0 * B.c6 - A.c1 * B.c5 - A.c2 * B.c7 + A.c2 * B.c15 -
                     A.c4 * B.c6 + A.c5 * B.c5 + A.c7 * B.c0 - A.c7 * B.c10;
    ctype const c4 = A.c0 * B.c9 - A.c1 * B.c8 + A.c3 * B.c7 + A.c3 * B.c15 +
                     A.c4 * B.c9 - A.c5 * B.c8 + A.c6 * B.c0 + A.c6 * B.c10;
    ctype const c5 = -A.c0 * B.c14 - A.c1 * B.c3 + A.c2 * B.c2 + A.c2 * B.c12 +
                     A.c4 * B.c14 + A.c5 * B.c3 - A.c7 * B.c1 - A.c7 * B.c11;
    ctype const c6 = A.c0 * B.c3 - A.c1 * B.c14 - A.c2 * B.c1 - A.c2 * B.c11 -
                     A.c4 * B.c3 + A.c5 * B.c14 - A.c7 * B.c2 - A.c7 * B.c12;
    ctype const c7 = -A.c0 * B.c11 - A.c1 * B.c12 - A.c2 * B.c4 + A.c3 * B.c3 -
                     A.c4 * B.c1 - A.c5 * B.c2 + A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c8 = A.c0 * B.c13 - A.c1 * B.c4 + A.c3 * B.c2 - A.c3 * B.c12 +
                     A.c4 * B.c13 - A.c5 * B.c4 + A.c6 * B.c1 - A.c6 * B.c11;
    ctype const c9 = A.c0 * B.c4 + A.c1 * B.c13 - A.c3 * B.c1 + A.c3 * B.c11 +
                     A.c4 * B.c4 + A.c5 * B.c13 + A.c6 * B.c2 - A.c6 * B.c12;
    ctype const c10 = A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c13 - A.c3 * B.c14 +
                      A.c4 * B.c12 - A.c5 * B.c11 + A.c6 * B.c3 - A.c7 * B.c4;
    ctype const c11 = -A.c0 * B.c7 + A.c1 * B.c0 - A.c2 * B.c9 + A.c3 * B.c6 +
                      A.c4 * B.c15 - A.c5 * B.c10 + A.c6 * B.c5 + A.c7 * B.c8;
    ctype const c12 = -A.c0 * B.c0 - A.c1 * B.c7 + A.c2 * B.c8 - A.c3 * B.c5 +
                      A.c4 * B.c10 + A.c5 * B.c15 + A.c6 * B.c6 + A.c7 * B.c9;
    ctype const c13 = A.c0 * B.c8 + A.c1 * B.c9 - A.c3 * B.c0 - A.c3 * B.c10 +
                      A.c4 * B.c8 + A.c5 * B.c9 + A.c6 * B.c7 + A.c6 * B.c15;
    ctype const c14 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c0 + A.c2 * B.c10 +
                      A.c4 * B.c5 + A.c5 * B.c6 - A.c7 * B.c7 + A.c7 * B.c15;
    ctype const c15 = A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c4 - A.c3 * B.c3 +
                      A.c4 * B.c11 + A.c5 * B.c12 - A.c6 * B.c14 - A.c7 * B.c13;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& M, PScalar2dc<U> ps)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps, MVec2dc<U> const& M)
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
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& M,
                                                 TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * t.y + M.c2 * t.x + M.c3 * t.z + M.c4 * t.w;
    ctype const c1 = -M.c0 * t.y + M.c5 * t.z - M.c7 * t.x - M.c8 * t.w;
    ctype const c2 = M.c0 * t.x + M.c6 * t.z - M.c7 * t.y - M.c9 * t.w;
    ctype const c3 = -M.c0 * t.w + M.c5 * t.y - M.c6 * t.x - M.c10 * t.w;
    ctype const c4 = -M.c0 * t.z - M.c8 * t.y + M.c9 * t.x + M.c10 * t.z;
    ctype const c5 = -M.c1 * t.w + M.c3 * t.y + M.c11 * t.w - M.c14 * t.x;
    ctype const c6 = -M.c2 * t.w - M.c3 * t.x + M.c12 * t.w - M.c14 * t.y;
    ctype const c7 = -M.c1 * t.x - M.c2 * t.y + M.c13 * t.w - M.c14 * t.z;
    ctype const c8 = M.c1 * t.z - M.c4 * t.y + M.c11 * t.z - M.c13 * t.x;
    ctype const c9 = M.c2 * t.z + M.c4 * t.x + M.c12 * t.z - M.c13 * t.y;
    ctype const c10 = M.c3 * t.z - M.c4 * t.w + M.c11 * t.y - M.c12 * t.x;
    ctype const c11 = -M.c5 * t.z - M.c8 * t.w + M.c10 * t.y + M.c15 * t.x;
    ctype const c12 = -M.c6 * t.z - M.c9 * t.w - M.c10 * t.x + M.c15 * t.y;
    ctype const c13 = -M.c7 * t.z - M.c8 * t.x - M.c9 * t.y + M.c15 * t.z;
    ctype const c14 = -M.c5 * t.x - M.c6 * t.y + M.c7 * t.w + M.c15 * t.w;
    ctype const c15 = M.c11 * t.x + M.c12 * t.y - M.c13 * t.w - M.c14 * t.z;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t,
                                                 MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c2 + t.y * M.c1 - t.z * M.c3 - t.w * M.c4;
    ctype const c1 = -t.x * M.c7 + t.y * M.c0 + t.z * M.c5 - t.w * M.c8;
    ctype const c2 = -t.x * M.c0 - t.y * M.c7 + t.z * M.c6 - t.w * M.c9;
    ctype const c3 = -t.x * M.c6 + t.y * M.c5 + t.w * M.c0 - t.w * M.c10;
    ctype const c4 = t.x * M.c9 - t.y * M.c8 + t.z * M.c0 + t.z * M.c10;
    ctype const c5 = t.x * M.c14 + t.y * M.c3 - t.w * M.c1 - t.w * M.c11;
    ctype const c6 = -t.x * M.c3 + t.y * M.c14 - t.w * M.c2 - t.w * M.c12;
    ctype const c7 = -t.x * M.c1 - t.y * M.c2 + t.z * M.c14 - t.w * M.c13;
    ctype const c8 = t.x * M.c13 - t.y * M.c4 + t.z * M.c1 - t.z * M.c11;
    ctype const c9 = t.x * M.c4 + t.y * M.c13 + t.z * M.c2 - t.z * M.c12;
    ctype const c10 = t.x * M.c12 - t.y * M.c11 + t.z * M.c3 - t.w * M.c4;
    ctype const c11 = t.x * M.c15 - t.y * M.c10 + t.z * M.c5 + t.w * M.c8;
    ctype const c12 = t.x * M.c10 + t.y * M.c15 + t.z * M.c6 + t.w * M.c9;
    ctype const c13 = t.x * M.c8 + t.y * M.c9 + t.z * M.c7 + t.z * M.c15;
    ctype const c14 = t.x * M.c5 + t.y * M.c6 - t.w * M.c7 + t.w * M.c15;
    ctype const c15 = t.x * M.c11 + t.y * M.c12 - t.z * M.c14 - t.w * M.c13;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& M,
                                                 BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * B.my - M.c6 * B.mx + M.c7 * B.mz - M.c8 * B.vy + M.c9 * B.vx +
                     M.c10 * B.vz;
    ctype const c1 = -M.c2 * B.mz + M.c3 * B.my + M.c4 * B.vy - M.c11 * B.vz +
                     M.c13 * B.vx - M.c14 * B.mx;
    ctype const c2 = M.c1 * B.mz - M.c3 * B.mx - M.c4 * B.vx - M.c12 * B.vz +
                     M.c13 * B.vy - M.c14 * B.my;
    ctype const c3 = M.c1 * B.vy - M.c2 * B.vx - M.c3 * B.vz - M.c11 * B.vy +
                     M.c12 * B.vx - M.c14 * B.mz;
    ctype const c4 = M.c1 * B.my - M.c2 * B.mx + M.c4 * B.vz + M.c11 * B.my -
                     M.c12 * B.mx + M.c13 * B.mz;
    ctype const c5 = M.c0 * B.vy - M.c5 * B.vz - M.c6 * B.mz + M.c7 * B.vx +
                     M.c10 * B.vy + M.c15 * B.vx;
    ctype const c6 = -M.c0 * B.vx + M.c5 * B.mz - M.c6 * B.vz + M.c7 * B.vy -
                     M.c10 * B.vx + M.c15 * B.vy;
    ctype const c7 = -M.c0 * B.mz + M.c5 * B.mx + M.c6 * B.my - M.c8 * B.vx -
                     M.c9 * B.vy + M.c15 * B.vz;
    ctype const c8 = -M.c0 * B.my - M.c7 * B.mx + M.c8 * B.vz - M.c9 * B.mz +
                     M.c10 * B.my + M.c15 * B.mx;
    ctype const c9 = M.c0 * B.mx - M.c7 * B.my + M.c8 * B.mz + M.c9 * B.vz -
                     M.c10 * B.mx + M.c15 * B.my;
    ctype const c10 = M.c0 * B.vz + M.c5 * B.my - M.c6 * B.mx + M.c8 * B.vy -
                      M.c9 * B.vx + M.c15 * B.mz;
    ctype const c11 = -M.c1 * B.vz - M.c3 * B.my + M.c4 * B.vy - M.c12 * B.mz +
                      M.c13 * B.vx + M.c14 * B.mx;
    ctype const c12 = -M.c2 * B.vz + M.c3 * B.mx - M.c4 * B.vx + M.c11 * B.mz +
                      M.c13 * B.vy + M.c14 * B.my;
    ctype const c13 = M.c1 * B.mx + M.c2 * B.my - M.c4 * B.mz + M.c11 * B.mx +
                      M.c12 * B.my + M.c13 * B.vz;
    ctype const c14 = -M.c1 * B.vx - M.c2 * B.vy + M.c3 * B.mz + M.c11 * B.vx +
                      M.c12 * B.vy - M.c14 * B.vz;
    ctype const c15 = M.c5 * B.mx + M.c6 * B.my + M.c7 * B.vz + M.c8 * B.vx +
                      M.c9 * B.vy - M.c10 * B.mz;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B,
                                                 MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c10 - B.mx * M.c6 +
                     B.my * M.c5 + B.mz * M.c7;
    ctype const c1 = B.vx * M.c13 - B.vy * M.c4 - B.vz * M.c11 - B.mx * M.c14 -
                     B.my * M.c3 + B.mz * M.c2;
    ctype const c2 = B.vx * M.c4 + B.vy * M.c13 - B.vz * M.c12 + B.mx * M.c3 -
                     B.my * M.c14 - B.mz * M.c1;
    ctype const c3 = B.vx * M.c2 + B.vx * M.c12 - B.vy * M.c1 - B.vy * M.c11 +
                     B.vz * M.c3 - B.mz * M.c14;
    ctype const c4 = -B.vz * M.c4 + B.mx * M.c2 - B.mx * M.c12 - B.my * M.c1 +
                     B.my * M.c11 + B.mz * M.c13;
    ctype const c5 = -B.vx * M.c7 + B.vx * M.c15 + B.vy * M.c0 - B.vy * M.c10 +
                     B.vz * M.c5 + B.mz * M.c6;
    ctype const c6 = -B.vx * M.c0 + B.vx * M.c10 - B.vy * M.c7 + B.vy * M.c15 +
                     B.vz * M.c6 - B.mz * M.c5;
    ctype const c7 = B.vx * M.c8 + B.vy * M.c9 + B.vz * M.c15 - B.mx * M.c5 -
                     B.my * M.c6 - B.mz * M.c0;
    ctype const c8 = -B.vz * M.c8 + B.mx * M.c7 + B.mx * M.c15 - B.my * M.c0 -
                     B.my * M.c10 + B.mz * M.c9;
    ctype const c9 = -B.vz * M.c9 + B.mx * M.c0 + B.mx * M.c10 + B.my * M.c7 +
                     B.my * M.c15 - B.mz * M.c8;
    ctype const c10 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c0 + B.mx * M.c6 -
                      B.my * M.c5 + B.mz * M.c15;
    ctype const c11 = -B.vx * M.c13 + B.vy * M.c4 - B.vz * M.c1 - B.mx * M.c14 -
                      B.my * M.c3 + B.mz * M.c12;
    ctype const c12 = -B.vx * M.c4 - B.vy * M.c13 - B.vz * M.c2 + B.mx * M.c3 -
                      B.my * M.c14 - B.mz * M.c11;
    ctype const c13 = -B.vz * M.c13 + B.mx * M.c1 - B.mx * M.c11 + B.my * M.c2 -
                      B.my * M.c12 - B.mz * M.c4;
    ctype const c14 = -B.vx * M.c1 - B.vx * M.c11 - B.vy * M.c2 - B.vy * M.c12 +
                      B.vz * M.c14 + B.mz * M.c3;
    ctype const c15 = B.vx * M.c8 + B.vy * M.c9 + B.vz * M.c7 + B.mx * M.c5 +
                      B.my * M.c6 - B.mz * M.c10;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& M, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * v.y + M.c12 * v.x - M.c13 * v.z - M.c14 * v.w;
    ctype const c1 = -M.c6 * v.w - M.c9 * v.z + M.c10 * v.y + M.c15 * v.x;
    ctype const c2 = M.c5 * v.w + M.c8 * v.z - M.c10 * v.x + M.c15 * v.y;
    ctype const c3 = M.c5 * v.y - M.c6 * v.x + M.c7 * v.z + M.c15 * v.z;
    ctype const c4 = -M.c7 * v.w + M.c8 * v.y - M.c9 * v.x + M.c15 * v.w;
    ctype const c5 = -M.c2 * v.z + M.c3 * v.y + M.c12 * v.z - M.c14 * v.x;
    ctype const c6 = M.c1 * v.z - M.c3 * v.x - M.c11 * v.z - M.c14 * v.y;
    ctype const c7 = -M.c3 * v.w + M.c4 * v.z - M.c11 * v.x - M.c12 * v.y;
    ctype const c8 = -M.c2 * v.w + M.c4 * v.y - M.c12 * v.w + M.c13 * v.x;
    ctype const c9 = M.c1 * v.w - M.c4 * v.x + M.c11 * v.w + M.c13 * v.y;
    ctype const c10 = M.c1 * v.y - M.c2 * v.x + M.c13 * v.z - M.c14 * v.w;
    ctype const c11 = -M.c0 * v.y + M.c6 * v.w - M.c7 * v.x - M.c9 * v.z;
    ctype const c12 = M.c0 * v.x - M.c5 * v.w - M.c7 * v.y + M.c8 * v.z;
    ctype const c13 = M.c0 * v.w + M.c8 * v.x + M.c9 * v.y - M.c10 * v.w;
    ctype const c14 = M.c0 * v.z - M.c5 * v.x - M.c6 * v.y + M.c10 * v.z;
    ctype const c15 = M.c1 * v.x + M.c2 * v.y - M.c3 * v.w - M.c4 * v.z;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * M.c12 + v.y * M.c11 + v.z * M.c13 + v.w * M.c14;
    ctype const c1 = v.x * M.c15 - v.y * M.c10 + v.z * M.c9 + v.w * M.c6;
    ctype const c2 = v.x * M.c10 + v.y * M.c15 - v.z * M.c8 - v.w * M.c5;
    ctype const c3 = v.x * M.c6 - v.y * M.c5 - v.z * M.c7 + v.z * M.c15;
    ctype const c4 = v.x * M.c9 - v.y * M.c8 + v.w * M.c7 + v.w * M.c15;
    ctype const c5 = -v.x * M.c14 - v.y * M.c3 + v.z * M.c2 + v.z * M.c12;
    ctype const c6 = v.x * M.c3 - v.y * M.c14 - v.z * M.c1 - v.z * M.c11;
    ctype const c7 = -v.x * M.c11 - v.y * M.c12 - v.z * M.c4 + v.w * M.c3;
    ctype const c8 = v.x * M.c13 - v.y * M.c4 + v.w * M.c2 - v.w * M.c12;
    ctype const c9 = v.x * M.c4 + v.y * M.c13 - v.w * M.c1 + v.w * M.c11;
    ctype const c10 = v.x * M.c2 - v.y * M.c1 + v.z * M.c13 - v.w * M.c14;
    ctype const c11 = -v.x * M.c7 + v.y * M.c0 - v.z * M.c9 + v.w * M.c6;
    ctype const c12 = -v.x * M.c0 - v.y * M.c7 + v.z * M.c8 - v.w * M.c5;
    ctype const c13 = v.x * M.c8 + v.y * M.c9 - v.w * M.c0 - v.w * M.c10;
    ctype const c14 = -v.x * M.c5 - v.y * M.c6 - v.z * M.c0 + v.z * M.c10;
    ctype const c15 = v.x * M.c1 + v.y * M.c2 - v.z * M.c4 - v.w * M.c3;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(MVec2dc<T> const& M, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c15 * ctype(s);
    ctype const c1 = M.c12 * ctype(s);
    ctype const c2 = -M.c11 * ctype(s);
    ctype const c3 = M.c14 * ctype(s);
    ctype const c4 = M.c13 * ctype(s);
    ctype const c5 = M.c6 * ctype(s);
    ctype const c6 = -M.c5 * ctype(s);
    ctype const c7 = -M.c10 * ctype(s);
    ctype const c8 = -M.c9 * ctype(s);
    ctype const c9 = M.c8 * ctype(s);
    ctype const c10 = M.c7 * ctype(s);
    ctype const c11 = M.c2 * ctype(s);
    ctype const c12 = -M.c1 * ctype(s);
    ctype const c13 = -M.c4 * ctype(s);
    ctype const c14 = -M.c3 * ctype(s);
    ctype const c15 = -M.c0 * ctype(s);
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c15;
    ctype const c1 = -ctype(s) * M.c12;
    ctype const c2 = ctype(s) * M.c11;
    ctype const c3 = -ctype(s) * M.c14;
    ctype const c4 = -ctype(s) * M.c13;
    ctype const c5 = ctype(s) * M.c6;
    ctype const c6 = -ctype(s) * M.c5;
    ctype const c7 = -ctype(s) * M.c10;
    ctype const c8 = -ctype(s) * M.c9;
    ctype const c9 = ctype(s) * M.c8;
    ctype const c10 = ctype(s) * M.c7;
    ctype const c11 = -ctype(s) * M.c2;
    ctype const c12 = ctype(s) * M.c1;
    ctype const c13 = ctype(s) * M.c4;
    ctype const c14 = ctype(s) * M.c3;
    ctype const c15 = -ctype(s) * M.c0;
    return MVec2dc<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// cga2dc rgpr :: rgpr(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& A,
                                                   MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c0 * B.c2 - A.c1 * B.c3 + A.c1 * B.c7 + A.c2 * B.c0 - A.c2 * B.c6 +
                     A.c3 * B.c1 + A.c6 * B.c2 + A.c7 * B.c1;
    ctype const c2 = -A.c0 * B.c1 - A.c1 * B.c0 + A.c1 * B.c6 - A.c2 * B.c3 +
                     A.c2 * B.c7 + A.c3 * B.c2 - A.c6 * B.c1 + A.c7 * B.c2;
    ctype const c3 = -A.c0 * B.c6 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c7 -
                     A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c0 + A.c7 * B.c3;
    ctype const c4 = -A.c0 * B.c5 - A.c3 * B.c4 + A.c4 * B.c3 + A.c4 * B.c7 -
                     A.c5 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c4;
    ctype const c5 = A.c0 * B.c4 - A.c3 * B.c5 + A.c4 * B.c0 + A.c4 * B.c6 + A.c5 * B.c3 +
                     A.c5 * B.c7 - A.c6 * B.c4 + A.c7 * B.c5;
    ctype const c6 = A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c0 + A.c4 * B.c2 -
                     A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = -A.c0 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c3 +
                     A.c4 * B.c1 + A.c5 * B.c2 - A.c6 * B.c6 + A.c7 * B.c7;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& A,
                                                   MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c5 + A.c1 * B.c6 - A.c2 * B.c3 - A.c3 * B.c4 -
                     A.c4 * B.c7 - A.c5 * B.c2 + A.c6 * B.c1 + A.c7 * B.c0;
    ctype const c1 = A.c0 * B.c4 + A.c1 * B.c3 + A.c2 * B.c6 - A.c3 * B.c5 + A.c4 * B.c2 -
                     A.c5 * B.c7 - A.c6 * B.c0 + A.c7 * B.c1;
    ctype const c2 = -A.c0 * B.c7 + A.c1 * B.c1 + A.c1 * B.c5 - A.c2 * B.c0 -
                     A.c2 * B.c4 + A.c3 * B.c2 - A.c6 * B.c7 + A.c7 * B.c2;
    ctype const c3 = -A.c0 * B.c6 - A.c3 * B.c3 + A.c4 * B.c1 - A.c4 * B.c5 -
                     A.c5 * B.c0 + A.c5 * B.c4 + A.c6 * B.c6 + A.c7 * B.c3;
    ctype const c4 = -A.c0 * B.c1 - A.c1 * B.c6 + A.c2 * B.c3 - A.c3 * B.c0 -
                     A.c4 * B.c7 - A.c5 * B.c2 + A.c6 * B.c5 + A.c7 * B.c4;
    ctype const c5 = A.c0 * B.c0 - A.c1 * B.c3 - A.c2 * B.c6 - A.c3 * B.c1 + A.c4 * B.c2 -
                     A.c5 * B.c7 - A.c6 * B.c4 + A.c7 * B.c5;
    ctype const c6 = A.c0 * B.c3 - A.c3 * B.c6 + A.c4 * B.c0 - A.c4 * B.c4 + A.c5 * B.c1 -
                     A.c5 * B.c5 - A.c6 * B.c3 + A.c7 * B.c6;
    ctype const c7 = A.c0 * B.c2 - A.c1 * B.c0 - A.c1 * B.c4 - A.c2 * B.c1 - A.c2 * B.c5 +
                     A.c3 * B.c7 + A.c6 * B.c2 + A.c7 * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& A,
                                                   MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c2 - A.c4 * B.c3 +
                     A.c5 * B.c0 + A.c6 * B.c1 - A.c7 * B.c4;
    ctype const c1 = A.c0 * B.c6 + A.c1 * B.c7 - A.c2 * B.c4 - A.c3 * B.c1 - A.c4 * B.c0 -
                     A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c5;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c3 + A.c2 * B.c7 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c7 * B.c0 - A.c7 * B.c6;
    ctype const c3 = A.c0 * B.c5 - A.c1 * B.c4 + A.c3 * B.c3 + A.c3 * B.c7 + A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c0 + A.c6 * B.c6;
    ctype const c4 = -A.c0 * B.c3 + A.c1 * B.c0 - A.c2 * B.c5 + A.c3 * B.c2 +
                     A.c4 * B.c7 - A.c5 * B.c6 + A.c6 * B.c1 + A.c7 * B.c4;
    ctype const c5 = -A.c0 * B.c0 - A.c1 * B.c3 + A.c2 * B.c4 - A.c3 * B.c1 +
                     A.c4 * B.c6 + A.c5 * B.c7 + A.c6 * B.c2 + A.c7 * B.c5;
    ctype const c6 = A.c0 * B.c4 + A.c1 * B.c5 - A.c3 * B.c0 - A.c3 * B.c6 + A.c4 * B.c4 +
                     A.c5 * B.c5 + A.c6 * B.c3 + A.c6 * B.c7;
    ctype const c7 = -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c0 + A.c2 * B.c6 +
                     A.c4 * B.c1 + A.c5 * B.c2 - A.c7 * B.c3 + A.c7 * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& A,
                                                   PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * ctype(ps);
    ctype const c1 = A.c1 * ctype(ps);
    ctype const c2 = A.c2 * ctype(ps);
    ctype const c3 = A.c3 * ctype(ps);
    ctype const c4 = A.c4 * ctype(ps);
    ctype const c5 = A.c5 * ctype(ps);
    ctype const c6 = A.c6 * ctype(ps);
    ctype const c7 = A.c7 * ctype(ps);
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps,
                                                   MVec2dc_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.c0;
    ctype const c1 = ctype(ps) * B.c1;
    ctype const c2 = ctype(ps) * B.c2;
    ctype const c3 = ctype(ps) * B.c3;
    ctype const c4 = ctype(ps) * B.c4;
    ctype const c5 = ctype(ps) * B.c5;
    ctype const c6 = ctype(ps) * B.c6;
    ctype const c7 = ctype(ps) * B.c7;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& M,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c0 * t.y + M.c1 * t.z - M.c3 * t.x - M.c4 * t.w;
    ctype const c1 = M.c0 * t.x + M.c2 * t.z - M.c3 * t.y - M.c5 * t.w;
    ctype const c2 = -M.c0 * t.w + M.c1 * t.y - M.c2 * t.x - M.c6 * t.w;
    ctype const c3 = -M.c0 * t.z - M.c4 * t.y + M.c5 * t.x + M.c6 * t.z;
    ctype const c4 = -M.c1 * t.z - M.c4 * t.w + M.c6 * t.y + M.c7 * t.x;
    ctype const c5 = -M.c2 * t.z - M.c5 * t.w - M.c6 * t.x + M.c7 * t.y;
    ctype const c6 = -M.c3 * t.z - M.c4 * t.x - M.c5 * t.y + M.c7 * t.z;
    ctype const c7 = -M.c1 * t.x - M.c2 * t.y + M.c3 * t.w + M.c7 * t.w;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t,
                                                   MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c3 + t.y * M.c0 + t.z * M.c1 - t.w * M.c4;
    ctype const c1 = -t.x * M.c0 - t.y * M.c3 + t.z * M.c2 - t.w * M.c5;
    ctype const c2 = -t.x * M.c2 + t.y * M.c1 + t.w * M.c0 - t.w * M.c6;
    ctype const c3 = t.x * M.c5 - t.y * M.c4 + t.z * M.c0 + t.z * M.c6;
    ctype const c4 = t.x * M.c7 - t.y * M.c6 + t.z * M.c1 + t.w * M.c4;
    ctype const c5 = t.x * M.c6 + t.y * M.c7 + t.z * M.c2 + t.w * M.c5;
    ctype const c6 = t.x * M.c4 + t.y * M.c5 + t.z * M.c3 + t.z * M.c7;
    ctype const c7 = t.x * M.c1 + t.y * M.c2 - t.w * M.c3 + t.w * M.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& M,
                                                   BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        M.c1 * B.my - M.c2 * B.mx + M.c3 * B.mz - M.c4 * B.vy + M.c5 * B.vx + M.c6 * B.vz;
    ctype const c1 =
        M.c0 * B.vy - M.c1 * B.vz - M.c2 * B.mz + M.c3 * B.vx + M.c6 * B.vy + M.c7 * B.vx;
    ctype const c2 = -M.c0 * B.vx + M.c1 * B.mz - M.c2 * B.vz + M.c3 * B.vy -
                     M.c6 * B.vx + M.c7 * B.vy;
    ctype const c3 = -M.c0 * B.mz + M.c1 * B.mx + M.c2 * B.my - M.c4 * B.vx -
                     M.c5 * B.vy + M.c7 * B.vz;
    ctype const c4 = -M.c0 * B.my - M.c3 * B.mx + M.c4 * B.vz - M.c5 * B.mz +
                     M.c6 * B.my + M.c7 * B.mx;
    ctype const c5 =
        M.c0 * B.mx - M.c3 * B.my + M.c4 * B.mz + M.c5 * B.vz - M.c6 * B.mx + M.c7 * B.my;
    ctype const c6 =
        M.c0 * B.vz + M.c1 * B.my - M.c2 * B.mx + M.c4 * B.vy - M.c5 * B.vx + M.c7 * B.mz;
    ctype const c7 =
        M.c1 * B.mx + M.c2 * B.my + M.c3 * B.vz + M.c4 * B.vx + M.c5 * B.vy - M.c6 * B.mz;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B,
                                                   MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c6 - B.mx * M.c2 + B.my * M.c1 + B.mz * M.c3;
    ctype const c1 = -B.vx * M.c3 + B.vx * M.c7 + B.vy * M.c0 - B.vy * M.c6 +
                     B.vz * M.c1 + B.mz * M.c2;
    ctype const c2 = -B.vx * M.c0 + B.vx * M.c6 - B.vy * M.c3 + B.vy * M.c7 +
                     B.vz * M.c2 - B.mz * M.c1;
    ctype const c3 =
        B.vx * M.c4 + B.vy * M.c5 + B.vz * M.c7 - B.mx * M.c1 - B.my * M.c2 - B.mz * M.c0;
    ctype const c4 = -B.vz * M.c4 + B.mx * M.c3 + B.mx * M.c7 - B.my * M.c0 -
                     B.my * M.c6 + B.mz * M.c5;
    ctype const c5 = -B.vz * M.c5 + B.mx * M.c0 + B.mx * M.c6 + B.my * M.c3 +
                     B.my * M.c7 - B.mz * M.c4;
    ctype const c6 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c0 + B.mx * M.c2 - B.my * M.c1 + B.mz * M.c7;
    ctype const c7 =
        B.vx * M.c4 + B.vy * M.c5 + B.vz * M.c3 + B.mx * M.c1 + B.my * M.c2 - B.mz * M.c6;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& M,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c2 * v.w - M.c5 * v.z + M.c6 * v.y + M.c7 * v.x;
    ctype const c1 = M.c1 * v.w + M.c4 * v.z - M.c6 * v.x + M.c7 * v.y;
    ctype const c2 = M.c1 * v.y - M.c2 * v.x + M.c3 * v.z + M.c7 * v.z;
    ctype const c3 = -M.c3 * v.w + M.c4 * v.y - M.c5 * v.x + M.c7 * v.w;
    ctype const c4 = -M.c0 * v.y + M.c2 * v.w - M.c3 * v.x - M.c5 * v.z;
    ctype const c5 = M.c0 * v.x - M.c1 * v.w - M.c3 * v.y + M.c4 * v.z;
    ctype const c6 = M.c0 * v.w + M.c4 * v.x + M.c5 * v.y - M.c6 * v.w;
    ctype const c7 = M.c0 * v.z - M.c1 * v.x - M.c2 * v.y + M.c6 * v.z;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v,
                                                   MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c7 - v.y * M.c6 + v.z * M.c5 + v.w * M.c2;
    ctype const c1 = v.x * M.c6 + v.y * M.c7 - v.z * M.c4 - v.w * M.c1;
    ctype const c2 = v.x * M.c2 - v.y * M.c1 - v.z * M.c3 + v.z * M.c7;
    ctype const c3 = v.x * M.c5 - v.y * M.c4 + v.w * M.c3 + v.w * M.c7;
    ctype const c4 = -v.x * M.c3 + v.y * M.c0 - v.z * M.c5 + v.w * M.c2;
    ctype const c5 = -v.x * M.c0 - v.y * M.c3 + v.z * M.c4 - v.w * M.c1;
    ctype const c6 = v.x * M.c4 + v.y * M.c5 - v.w * M.c0 - v.w * M.c6;
    ctype const c7 = -v.x * M.c1 - v.y * M.c2 - v.z * M.c0 + v.z * M.c6;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(MVec2dc_E<T> const& M, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c7 * ctype(s);
    ctype const c1 = M.c2 * ctype(s);
    ctype const c2 = -M.c1 * ctype(s);
    ctype const c3 = -M.c6 * ctype(s);
    ctype const c4 = -M.c5 * ctype(s);
    ctype const c5 = M.c4 * ctype(s);
    ctype const c6 = M.c3 * ctype(s);
    ctype const c7 = -M.c0 * ctype(s);
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s, MVec2dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c7;
    ctype const c1 = ctype(s) * M.c2;
    ctype const c2 = -ctype(s) * M.c1;
    ctype const c3 = -ctype(s) * M.c6;
    ctype const c4 = -ctype(s) * M.c5;
    ctype const c5 = ctype(s) * M.c4;
    ctype const c6 = ctype(s) * M.c3;
    ctype const c7 = -ctype(s) * M.c0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& A,
                                                   MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c6 + A.c3 * B.c7 -
                     A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2 - A.c7 * B.c3;
    ctype const c1 = -A.c0 * B.c7 - A.c1 * B.c2 + A.c2 * B.c1 + A.c2 * B.c5 +
                     A.c4 * B.c7 + A.c5 * B.c2 - A.c7 * B.c0 - A.c7 * B.c4;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c7 - A.c2 * B.c0 - A.c2 * B.c4 - A.c4 * B.c2 +
                     A.c5 * B.c7 - A.c7 * B.c1 - A.c7 * B.c5;
    ctype const c3 = -A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c3 + A.c3 * B.c2 -
                     A.c4 * B.c0 - A.c5 * B.c1 + A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c4 = A.c0 * B.c6 - A.c1 * B.c3 + A.c3 * B.c1 - A.c3 * B.c5 + A.c4 * B.c6 -
                     A.c5 * B.c3 + A.c6 * B.c0 - A.c6 * B.c4;
    ctype const c5 = A.c0 * B.c3 + A.c1 * B.c6 - A.c3 * B.c0 + A.c3 * B.c4 + A.c4 * B.c3 +
                     A.c5 * B.c6 + A.c6 * B.c1 - A.c6 * B.c5;
    ctype const c6 = A.c0 * B.c1 - A.c1 * B.c0 + A.c2 * B.c6 - A.c3 * B.c7 + A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c2 - A.c7 * B.c3;
    ctype const c7 = A.c0 * B.c0 + A.c1 * B.c1 - A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c4 +
                     A.c5 * B.c5 - A.c6 * B.c7 - A.c7 * B.c6;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& A,
                                                   PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * ctype(ps);
    ctype const c1 = A.c1 * ctype(ps);
    ctype const c2 = A.c2 * ctype(ps);
    ctype const c3 = A.c3 * ctype(ps);
    ctype const c4 = A.c4 * ctype(ps);
    ctype const c5 = A.c5 * ctype(ps);
    ctype const c6 = A.c6 * ctype(ps);
    ctype const c7 = A.c7 * ctype(ps);
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps,
                                                   MVec2dc_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.c0;
    ctype const c1 = ctype(ps) * B.c1;
    ctype const c2 = ctype(ps) * B.c2;
    ctype const c3 = ctype(ps) * B.c3;
    ctype const c4 = ctype(ps) * B.c4;
    ctype const c5 = ctype(ps) * B.c5;
    ctype const c6 = ctype(ps) * B.c6;
    ctype const c7 = ctype(ps) * B.c7;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& M,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c0 * t.y + M.c1 * t.x + M.c2 * t.z + M.c3 * t.w;
    ctype const c1 = -M.c0 * t.w + M.c2 * t.y + M.c4 * t.w - M.c7 * t.x;
    ctype const c2 = -M.c1 * t.w - M.c2 * t.x + M.c5 * t.w - M.c7 * t.y;
    ctype const c3 = -M.c0 * t.x - M.c1 * t.y + M.c6 * t.w - M.c7 * t.z;
    ctype const c4 = M.c0 * t.z - M.c3 * t.y + M.c4 * t.z - M.c6 * t.x;
    ctype const c5 = M.c1 * t.z + M.c3 * t.x + M.c5 * t.z - M.c6 * t.y;
    ctype const c6 = M.c2 * t.z - M.c3 * t.w + M.c4 * t.y - M.c5 * t.x;
    ctype const c7 = M.c4 * t.x + M.c5 * t.y - M.c6 * t.w - M.c7 * t.z;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t,
                                                   MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c1 + t.y * M.c0 - t.z * M.c2 - t.w * M.c3;
    ctype const c1 = t.x * M.c7 + t.y * M.c2 - t.w * M.c0 - t.w * M.c4;
    ctype const c2 = -t.x * M.c2 + t.y * M.c7 - t.w * M.c1 - t.w * M.c5;
    ctype const c3 = -t.x * M.c0 - t.y * M.c1 + t.z * M.c7 - t.w * M.c6;
    ctype const c4 = t.x * M.c6 - t.y * M.c3 + t.z * M.c0 - t.z * M.c4;
    ctype const c5 = t.x * M.c3 + t.y * M.c6 + t.z * M.c1 - t.z * M.c5;
    ctype const c6 = t.x * M.c5 - t.y * M.c4 + t.z * M.c2 - t.w * M.c3;
    ctype const c7 = t.x * M.c4 + t.y * M.c5 - t.z * M.c7 - t.w * M.c6;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& M,
                                                   BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.mz + M.c2 * B.my + M.c3 * B.vy - M.c4 * B.vz +
                     M.c6 * B.vx - M.c7 * B.mx;
    ctype const c1 =
        M.c0 * B.mz - M.c2 * B.mx - M.c3 * B.vx - M.c5 * B.vz + M.c6 * B.vy - M.c7 * B.my;
    ctype const c2 =
        M.c0 * B.vy - M.c1 * B.vx - M.c2 * B.vz - M.c4 * B.vy + M.c5 * B.vx - M.c7 * B.mz;
    ctype const c3 =
        M.c0 * B.my - M.c1 * B.mx + M.c3 * B.vz + M.c4 * B.my - M.c5 * B.mx + M.c6 * B.mz;
    ctype const c4 = -M.c0 * B.vz - M.c2 * B.my + M.c3 * B.vy - M.c5 * B.mz +
                     M.c6 * B.vx + M.c7 * B.mx;
    ctype const c5 = -M.c1 * B.vz + M.c2 * B.mx - M.c3 * B.vx + M.c4 * B.mz +
                     M.c6 * B.vy + M.c7 * B.my;
    ctype const c6 =
        M.c0 * B.mx + M.c1 * B.my - M.c3 * B.mz + M.c4 * B.mx + M.c5 * B.my + M.c6 * B.vz;
    ctype const c7 = -M.c0 * B.vx - M.c1 * B.vy + M.c2 * B.mz + M.c4 * B.vx +
                     M.c5 * B.vy - M.c7 * B.vz;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B,
                                                   MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c6 - B.vy * M.c3 - B.vz * M.c4 - B.mx * M.c7 - B.my * M.c2 + B.mz * M.c1;
    ctype const c1 =
        B.vx * M.c3 + B.vy * M.c6 - B.vz * M.c5 + B.mx * M.c2 - B.my * M.c7 - B.mz * M.c0;
    ctype const c2 =
        B.vx * M.c1 + B.vx * M.c5 - B.vy * M.c0 - B.vy * M.c4 + B.vz * M.c2 - B.mz * M.c7;
    ctype const c3 = -B.vz * M.c3 + B.mx * M.c1 - B.mx * M.c5 - B.my * M.c0 +
                     B.my * M.c4 + B.mz * M.c6;
    ctype const c4 = -B.vx * M.c6 + B.vy * M.c3 - B.vz * M.c0 - B.mx * M.c7 -
                     B.my * M.c2 + B.mz * M.c5;
    ctype const c5 = -B.vx * M.c3 - B.vy * M.c6 - B.vz * M.c1 + B.mx * M.c2 -
                     B.my * M.c7 - B.mz * M.c4;
    ctype const c6 = -B.vz * M.c6 + B.mx * M.c0 - B.mx * M.c4 + B.my * M.c1 -
                     B.my * M.c5 - B.mz * M.c3;
    ctype const c7 = -B.vx * M.c0 - B.vx * M.c4 - B.vy * M.c1 - B.vy * M.c5 +
                     B.vz * M.c7 + B.mz * M.c2;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& M,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * v.y + M.c5 * v.x - M.c6 * v.z - M.c7 * v.w;
    ctype const c1 = -M.c1 * v.z + M.c2 * v.y + M.c5 * v.z - M.c7 * v.x;
    ctype const c2 = M.c0 * v.z - M.c2 * v.x - M.c4 * v.z - M.c7 * v.y;
    ctype const c3 = -M.c2 * v.w + M.c3 * v.z - M.c4 * v.x - M.c5 * v.y;
    ctype const c4 = -M.c1 * v.w + M.c3 * v.y - M.c5 * v.w + M.c6 * v.x;
    ctype const c5 = M.c0 * v.w - M.c3 * v.x + M.c4 * v.w + M.c6 * v.y;
    ctype const c6 = M.c0 * v.y - M.c1 * v.x + M.c6 * v.z - M.c7 * v.w;
    ctype const c7 = M.c0 * v.x + M.c1 * v.y - M.c2 * v.w - M.c3 * v.z;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v,
                                                   MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * M.c5 + v.y * M.c4 + v.z * M.c6 + v.w * M.c7;
    ctype const c1 = -v.x * M.c7 - v.y * M.c2 + v.z * M.c1 + v.z * M.c5;
    ctype const c2 = v.x * M.c2 - v.y * M.c7 - v.z * M.c0 - v.z * M.c4;
    ctype const c3 = -v.x * M.c4 - v.y * M.c5 - v.z * M.c3 + v.w * M.c2;
    ctype const c4 = v.x * M.c6 - v.y * M.c3 + v.w * M.c1 - v.w * M.c5;
    ctype const c5 = v.x * M.c3 + v.y * M.c6 - v.w * M.c0 + v.w * M.c4;
    ctype const c6 = v.x * M.c1 - v.y * M.c0 + v.z * M.c6 - v.w * M.c7;
    ctype const c7 = v.x * M.c0 + v.y * M.c1 - v.z * M.c3 - v.w * M.c2;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(MVec2dc_U<T> const& M, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * ctype(s);
    ctype const c1 = -M.c4 * ctype(s);
    ctype const c2 = M.c7 * ctype(s);
    ctype const c3 = M.c6 * ctype(s);
    ctype const c4 = M.c1 * ctype(s);
    ctype const c5 = -M.c0 * ctype(s);
    ctype const c6 = -M.c3 * ctype(s);
    ctype const c7 = -M.c2 * ctype(s);
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s, MVec2dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * M.c5;
    ctype const c1 = ctype(s) * M.c4;
    ctype const c2 = -ctype(s) * M.c7;
    ctype const c3 = -ctype(s) * M.c6;
    ctype const c4 = -ctype(s) * M.c1;
    ctype const c5 = ctype(s) * M.c0;
    ctype const c6 = ctype(s) * M.c3;
    ctype const c7 = ctype(s) * M.c2;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps1, PScalar2dc<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(ctype(ps1) * ctype(ps2));
}

// cga2dc rgpr :: rgpr(ps,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * t.x;
    ctype const c1 = ctype(ps) * t.y;
    ctype const c2 = ctype(ps) * t.z;
    ctype const c3 = ctype(ps) * t.w;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(trivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t,
                                                   PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(ps);
    ctype const c1 = t.y * ctype(ps);
    ctype const c2 = t.z * ctype(ps);
    ctype const c3 = t.w * ctype(ps);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.vx;
    ctype const c1 = ctype(ps) * B.vy;
    ctype const c2 = ctype(ps) * B.vz;
    ctype const c3 = ctype(ps) * B.mx;
    ctype const c4 = ctype(ps) * B.my;
    ctype const c5 = ctype(ps) * B.mz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rgpr :: rgpr(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(ps);
    ctype const c1 = B.vy * ctype(ps);
    ctype const c2 = B.vz * ctype(ps);
    ctype const c3 = B.mx * ctype(ps);
    ctype const c4 = B.my * ctype(ps);
    ctype const c5 = B.mz * ctype(ps);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rgpr :: rgpr(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = ctype(ps) * v.w;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = v.w * ctype(ps);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rgpr(PScalar2dc<T> ps, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(ps) * ctype(s));
}

// cga2dc rgpr :: rgpr(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dc<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s, PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dc<ctype>(ctype(s) * ctype(ps));
}

// cga2dc rgpr :: rgpr(trivec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t1,
                                                   TriVec2dc<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = t1.x * t2.w - t1.w * t2.x;
    ctype const c2 = t1.y * t2.w - t1.w * t2.y;
    ctype const c3 = t1.z * t2.w - t1.w * t2.z;
    ctype const c4 = t1.x * t2.z - t1.z * t2.x;
    ctype const c5 = t1.y * t2.z - t1.z * t2.y;
    ctype const c6 = t1.x * t2.y - t1.y * t2.x;
    ctype const c7 = t1.x * t2.x + t1.y * t2.y - t1.z * t2.w - t1.w * t2.z;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(trivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t,
                                                   BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * B.vz + t.z * B.vx - t.w * B.mx;
    ctype const c1 = -t.y * B.vz + t.z * B.vy - t.w * B.my;
    ctype const c2 = -t.x * B.vy + t.y * B.vx - t.w * B.mz;
    ctype const c3 = t.x * B.my - t.y * B.mx + t.z * B.mz;
    ctype const c4 = -t.y * B.mz + t.z * B.vx + t.w * B.mx;
    ctype const c5 = t.x * B.mz + t.z * B.vy + t.w * B.my;
    ctype const c6 = t.x * B.mx + t.y * B.my + t.z * B.vz;
    ctype const c7 = t.x * B.vx + t.y * B.vy - t.w * B.vz;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(bivec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * t.z - B.vz * t.x - B.mx * t.w;
    ctype const c1 = B.vy * t.z - B.vz * t.y - B.my * t.w;
    ctype const c2 = B.vx * t.y - B.vy * t.x - B.mz * t.w;
    ctype const c3 = -B.mx * t.y + B.my * t.x + B.mz * t.z;
    ctype const c4 = -B.vx * t.z - B.mx * t.w + B.mz * t.y;
    ctype const c5 = -B.vy * t.z - B.my * t.w - B.mz * t.x;
    ctype const c6 = -B.vz * t.z - B.mx * t.x - B.my * t.y;
    ctype const c7 = -B.vx * t.x - B.vy * t.y + B.vz * t.w;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(trivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * v.y + t.y * v.x - t.z * v.z - t.w * v.w;
    ctype const c1 = t.y * v.z - t.w * v.x;
    ctype const c2 = -t.x * v.z - t.w * v.y;
    ctype const c3 = -t.x * v.x - t.y * v.y;
    ctype const c4 = -t.y * v.w + t.z * v.x;
    ctype const c5 = t.x * v.w + t.z * v.y;
    ctype const c6 = t.z * v.z - t.w * v.w;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(vec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v,
                                                   TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * t.y + v.y * t.x + v.z * t.z + v.w * t.w;
    ctype const c1 = -v.x * t.w + v.z * t.y;
    ctype const c2 = -v.y * t.w - v.z * t.x;
    ctype const c3 = -v.x * t.x - v.y * t.y;
    ctype const c4 = v.x * t.z - v.w * t.y;
    ctype const c5 = v.y * t.z + v.w * t.x;
    ctype const c6 = v.z * t.z - v.w * t.w;
    ctype const c7 = 0.0;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(trivec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rgpr(TriVec2dc<T> const& t, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * ctype(s);
    ctype const c1 = -t.x * ctype(s);
    ctype const c2 = t.w * ctype(s);
    ctype const c3 = t.z * ctype(s);
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(s,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s, TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * t.y;
    ctype const c1 = ctype(s) * t.x;
    ctype const c2 = -ctype(s) * t.w;
    ctype const c3 = -ctype(s) * t.z;
    return Vec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(bivec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B1,
                                                   BiVec2dc<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vx * B2.my - B1.vy * B2.mx + B1.vz * B2.mz - B1.mx * B2.vy +
                     B1.my * B2.vx + B1.mz * B2.vz;
    ctype const c1 = -B1.vx * B2.vz - B1.vy * B2.mz + B1.vz * B2.vx + B1.mz * B2.vy;
    ctype const c2 = B1.vx * B2.mz - B1.vy * B2.vz + B1.vz * B2.vy - B1.mz * B2.vx;
    ctype const c3 = B1.vx * B2.mx + B1.vy * B2.my - B1.mx * B2.vx - B1.my * B2.vy;
    ctype const c4 = -B1.vz * B2.mx + B1.mx * B2.vz - B1.my * B2.mz + B1.mz * B2.my;
    ctype const c5 = -B1.vz * B2.my + B1.mx * B2.mz + B1.my * B2.vz - B1.mz * B2.mx;
    ctype const c6 = B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx;
    ctype const c7 = B1.vx * B2.mx + B1.vy * B2.my + B1.vz * B2.vz + B1.mx * B2.vx +
                     B1.my * B2.vy - B1.mz * B2.mz;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(bivec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B,
                                                   Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * v.w - B.my * v.z + B.mz * v.y;
    ctype const c1 = B.vx * v.w + B.mx * v.z - B.mz * v.x;
    ctype const c2 = B.vx * v.y - B.vy * v.x + B.vz * v.z;
    ctype const c3 = -B.vz * v.w + B.mx * v.y - B.my * v.x;
    ctype const c4 = B.vy * v.w - B.vz * v.x - B.my * v.z;
    ctype const c5 = -B.vx * v.w - B.vz * v.y + B.mx * v.z;
    ctype const c6 = B.mx * v.x + B.my * v.y - B.mz * v.w;
    ctype const c7 = -B.vx * v.x - B.vy * v.y + B.mz * v.z;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(vec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_U<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v,
                                                   BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my + v.w * B.vy;
    ctype const c1 = v.x * B.mz - v.z * B.mx - v.w * B.vx;
    ctype const c2 = v.x * B.vy - v.y * B.vx - v.z * B.vz;
    ctype const c3 = v.x * B.my - v.y * B.mx + v.w * B.vz;
    ctype const c4 = -v.x * B.vz - v.z * B.my + v.w * B.vy;
    ctype const c5 = -v.y * B.vz + v.z * B.mx - v.w * B.vx;
    ctype const c6 = v.x * B.mx + v.y * B.my - v.w * B.mz;
    ctype const c7 = -v.x * B.vx - v.y * B.vy + v.z * B.mz;
    return MVec2dc_U<ctype>(Vec2dc<ctype>(c0, c1, c2, c3),
                            TriVec2dc<ctype>(c4, c5, c6, c7));
}

// cga2dc rgpr :: rgpr(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rgpr(BiVec2dc<T> const& B, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * ctype(s);
    ctype const c1 = -B.vx * ctype(s);
    ctype const c2 = -B.mz * ctype(s);
    ctype const c3 = -B.my * ctype(s);
    ctype const c4 = B.mx * ctype(s);
    ctype const c5 = B.vz * ctype(s);
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rgpr :: rgpr(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s, BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vy;
    ctype const c1 = -ctype(s) * B.vx;
    ctype const c2 = -ctype(s) * B.mz;
    ctype const c3 = -ctype(s) * B.my;
    ctype const c4 = ctype(s) * B.mx;
    ctype const c5 = ctype(s) * B.vz;
    return BiVec2dc<ctype>(c0, c1, c2, c3, c4, c5);
}

// cga2dc rgpr :: rgpr(vec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc_E<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v1,
                                                   Vec2dc<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -v1.y * v2.z + v1.z * v2.y;
    ctype const c2 = v1.x * v2.z - v1.z * v2.x;
    ctype const c3 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c4 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c5 = v1.x * v2.w - v1.w * v2.x;
    ctype const c6 = v1.x * v2.y - v1.y * v2.x;
    ctype const c7 = v1.x * v2.x + v1.y * v2.y - v1.z * v2.w - v1.w * v2.z;
    return MVec2dc_E<ctype>(Scalar2dc<ctype>(c0), BiVec2dc<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar2dc<ctype>(c7));
}

// cga2dc rgpr :: rgpr(vec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rgpr(Vec2dc<T> const& v, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * ctype(s);
    ctype const c1 = -v.x * ctype(s);
    ctype const c2 = -v.w * ctype(s);
    ctype const c3 = -v.z * ctype(s);
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(s,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(s) * v.y;
    ctype const c1 = ctype(s) * v.x;
    ctype const c2 = ctype(s) * v.w;
    ctype const c3 = ctype(s) * v.z;
    return TriVec2dc<ctype>(c0, c1, c2, c3);
}

// cga2dc rgpr :: rgpr(s,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dc<std::common_type_t<T, U>> rgpr(Scalar2dc<T> s1, Scalar2dc<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dc<ctype>(-ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors w.r.t. the
// geometric product: u * inv(u) = inv(u) * u = 1
//
// The cga2dc metric is non-degenerate (det G = -1), so every element with
// non-zero geometric square is invertible. NULL elements are not -- and in a
// conformal algebra these are geometrically meaningful, not exotic: every
// embedded point q satisfies q * q = 0, so inv() of a round point (and of the
// null basis vectors e3, e4 themselves) throws. Invert only elements known to
// be non-null (e.g. versors).
//
// for k-blades:  A^(-1) = rev(A) / <A rev(A)>_0
//   - scalar / vector: rev = +A, so A^(-1) = A / gr0(A A)
//   - bivector / trivector and the mixed-grade multivectors: closed-form
//     inverse of Hitzer & Sangwine, "Multivector and multivector matrix
//     inverses in real Clifford algebras" (2016) -- A^(-1) = conj(A) * map /
//     <A conj(A) map>_0 with the dim-4 grade-sign maps. The formula is an
//     algebraic identity in the algebra's own geometric product, so it holds
//     for the non-orthogonal null-basis metric unchanged.
//
// The squared norm is SIGNED here (e.g. I_2dc * I_2dc = -1), so the null check
// uses std::abs(denominator): it rejects only (near-)null elements, never
// merely negative squares.
////////////////////////////////////////////////////////////////////////////////
// HINT: inv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Scalar2dc<T> inv(Scalar2dc<T> s)
{
    T sq_n = T(s) * T(s);
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "scalar");
    return Scalar2dc<T>(T(s) / sq_n);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec2dc<T> inv(Vec2dc<T> const& v)
{
    // v^(-1) = rev(v)/<v v>_0 = v/dot(v,v); throws for null vectors -- which
    // includes every embedded point and the null basis vectors e3, e4
    T sq_n = T(dot(v, v));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "vector");
    T inv = T(1.0) / sq_n;
    return Vec2dc<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline BiVec2dc<T> inv(BiVec2dc<T> const& B)
{
    auto bc = B * conj(B);
    auto bcmap = gr0(bc) + gr2(bc) - gr4(bc);
    T sq_n = T(gr0(bc * bcmap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "bivector");
    return gr2(conj(B) * bcmap) / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline TriVec2dc<T> inv(TriVec2dc<T> const& t)
{
    auto tc = t * conj(t);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "trivector");
    return gr3(conj(t) * tcmap) / sq_n;
}

// the cga2dc metric is non-degenerate, so (unlike the pga pseudoscalars) the
// pseudoscalar is invertible: I_2dc * I_2dc = -1, hence inv(I_2dc) = -I_2dc
template <typename T>
    requires(numeric_type<T>)
inline PScalar2dc<T> inv(PScalar2dc<T> ps)
{
    T sq_n = T(ps * ps); // = -ps^2 (signed geometric square)
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "pseudoscalar");
    return PScalar2dc<T>(T(ps) / sq_n);
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_E<T> inv(MVec2dc_E<T> const& E)
{
    auto tc = E * conj(E);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "even-grade multivector");
    return conj(E) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_U<T> inv(MVec2dc_U<T> const& U)
{
    auto tc = U * conj(U);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "odd-grade multivector");
    return conj(U) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Clifford
// algebras", Hitzer, Sangwine, 2016; a left and a right inverse are the same
// (see the paper)
template <typename T>
    requires(numeric_type<T>)
inline MVec2dc<T> inv(MVec2dc<T> const& M)
{
    auto tc = M * conj(M);
    auto tcmap = gr0(tc) + gr1(tc) + gr2(tc) - gr3(tc) - gr4(tc);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(std::abs(sq_n), "multivector");
    return conj(M) * tcmap / sq_n;
}


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses w.r.t. the regressive geometric product:
//
//     u (v) rinv(u) = rinv(u) (v) u = I_2dc   (the pseudoscalar, identity of rgpr)
//
// defined as the geometric inverse carried through the complement map
// (even-dimensional algebra, so l_/r_ complements pair up):
//
//     rinv(u) = l_cmpl(inv(r_cmpl(u))) = r_cmpl(inv(l_cmpl(u)))
//
// Unlike pga (where the scalar has no rinv, mirroring the pseudoscalar having
// no inv), the non-degenerate cga2dc metric admits BOTH inverses for every
// non-null element, so rinv covers all types -- e.g. rinv(1) = -1 (a scalar:
// rgpr(1,-1) = I_2dc), the complement-dual of inv(I_2dc) = -I_2dc. Null
// elements throw, as for inv().
////////////////////////////////////////////////////////////////////////////////
// HINT: rinv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Scalar2dc<T> rinv(Scalar2dc<T> s)
{
    return l_cmpl(inv(r_cmpl(s)));
}

template <typename T>
    requires(numeric_type<T>)
inline Vec2dc<T> rinv(Vec2dc<T> const& v)
{
    return l_cmpl(inv(r_cmpl(v)));
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec2dc<T> rinv(BiVec2dc<T> const& B)
{
    return l_cmpl(inv(r_cmpl(B)));
}

template <typename T>
    requires(numeric_type<T>)
inline TriVec2dc<T> rinv(TriVec2dc<T> const& t)
{
    return l_cmpl(inv(r_cmpl(t)));
}

template <typename T>
    requires(numeric_type<T>)
inline PScalar2dc<T> rinv(PScalar2dc<T> ps)
{
    return l_cmpl(inv(r_cmpl(ps)));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_E<T> rinv(MVec2dc_E<T> const& E)
{
    return l_cmpl(inv(r_cmpl(E)));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_U<T> rinv(MVec2dc_U<T> const& U)
{
    return l_cmpl(inv(r_cmpl(U)));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2dc<T> rinv(MVec2dc<T> const& M)
{
    return l_cmpl(inv(r_cmpl(M)));
}


////////////////////////////////////////////////////////////////////////////////
// left and right expansions
//
// Expansions are the wdg-based duals of the contractions (which use rwdg):
//
//     l_expand2dc(a, b) = wdg( l_dual(a), b )
//     r_expand2dc(a, b) = wdg( a, r_dual(b) )
//
// Implemented as generic wrappers over the existing dual and wedge operations,
// so every grade-pair combination is covered automatically. The requires-clause
// admits exactly the operand types for which the composition is well-formed
// (i.e. the cga2dc blade and multivector types).
////////////////////////////////////////////////////////////////////////////////

template <typename A, typename B>
    requires requires(A const& a, B const& b) { wdg(l_dual(a), b); }
constexpr auto l_expand2dc(A const& a, B const& b)
{
    return wdg(l_dual(a), b);
}

template <typename A, typename B>
    requires requires(A const& a, B const& b) { wdg(a, r_dual(b)); }
constexpr auto r_expand2dc(A const& a, B const& b)
{
    return wdg(a, r_dual(b));
}

} // namespace hd::ga::cga
