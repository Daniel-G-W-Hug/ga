#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_sta4ds_ops_basics.hpp"


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga4ds product operations:
//
// - dot()                   -> dot product
// - rdot()                  -> regressive dot product
// - wdg()                   -> wedge product (join as convenience interface)
// - rwdg()                  -> regressive wedge product (meet as convenience interface)
// - operator<<()            -> left contraction
// - operator>>()            -> right contraction
// - cmt()                   -> commutator product (= asymmetric part of gpr)
// - rcmt()                  -> regressive commutator product (=asymmetric part of rgpr)
// - operator*()             -> geometric product (= gpr)
// - rgpr()                  -> regressive geometric product
// - inv()                   -> inversion operation (w.r.t. geometric product)
//
// Expansions are the regressive versions of the corresponding contractions.
//
// - l_expand4ds()           -> left expansion
// - r_expand4ds()           -> right expansion
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot products for 4ds (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume g1^2 = -1, g2^2 = -1, g3^2 = -1, g4^2 = +1
//
// Via the exomorphism G(a ^ b) = G(a) ^ G(b) the extended metric G is defined as a 2^n by
// 2^n matrix via the outer product, which in turn allows for definition of an extended
// dot product for all grades of a full multi-vector from scalar to pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

// sta4ds dot :: dot(mv,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(MVec4ds<T> const& A,
                                                  MVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 +
                            A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7 -
                            A.c8 * B.c8 - A.c9 * B.c9 - A.c10 * B.c10 - A.c11 * B.c11 -
                            A.c12 * B.c12 - A.c13 * B.c13 + A.c14 * B.c14 -
                            A.c15 * B.c15);
}

// sta4ds dot :: dot(mv_e,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(MVec4ds_E<T> const& A,
                                                  MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 -
                            A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7);
}

// sta4ds dot :: dot(mv_u,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(MVec4ds_U<T> const& A,
                                                  MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(-A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c3 -
                            A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 + A.c7 * B.c7);
}

// sta4ds dot :: dot(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(PScalar4ds<T> ps1, PScalar4ds<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(-ctype(ps1) * ctype(ps2));
}

// sta4ds dot :: dot(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(TriVec4ds<T> const& t1,
                                                  TriVec4ds<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(-t1.x * t2.x - t1.y * t2.y - t1.z * t2.z + t1.w * t2.w);
}

// sta4ds dot :: dot(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(BiVec4ds<T> const& B1,
                                                  BiVec4ds<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(B1.vx * B2.vx + B1.vy * B2.vy + B1.vz * B2.vz -
                            B1.mx * B2.mx - B1.my * B2.my - B1.mz * B2.mz);
}

// sta4ds dot :: dot(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(Vec4ds<T> const& v1,
                                                  Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(-v1.x * v2.x - v1.y * v2.y - v1.z * v2.z + v1.w * v2.w);
}

// sta4ds dot :: dot(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> dot(Scalar4ds<T> s1, Scalar4ds<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive dot products for 4ds (= defined for equal grades exclusively)
//
// rdot(v1,v2) = l_cmpl( dot(r_cmpl(v1),r_cmpl(v2)) )
//
// returns a pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product for complements) and meet operations
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product, just depending on the outer product (wdg)
// as well as the complement and thus the pseudoscalar of the space
// (this definition does NOT connect directly to the geometric product,
// but to the outer product exclusively)
//
//        rwdg(ul, ur) = l_cmpl(wdg(r_cmpl(ul),r_cmpl(ur))) = l_cmpl(r_cmpl(ul) ^
//        r_cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// sta4ds wdg :: wdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& A, MVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c4 - A.c4 * B.c1 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c2 * B.c4 - A.c4 * B.c2 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c3 * B.c4 - A.c4 * B.c3 + A.c7 * B.c0;
    ctype const c8 = A.c0 * B.c8 + A.c2 * B.c3 - A.c3 * B.c2 + A.c8 * B.c0;
    ctype const c9 = A.c0 * B.c9 - A.c1 * B.c3 + A.c3 * B.c1 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c10 + A.c1 * B.c2 - A.c2 * B.c1 + A.c10 * B.c0;
    ctype const c11 = A.c0 * B.c11 + A.c2 * B.c7 - A.c3 * B.c6 + A.c4 * B.c8 -
                      A.c6 * B.c3 + A.c7 * B.c2 + A.c8 * B.c4 + A.c11 * B.c0;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c7 + A.c3 * B.c5 + A.c4 * B.c9 +
                      A.c5 * B.c3 - A.c7 * B.c1 + A.c9 * B.c4 + A.c12 * B.c0;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c6 - A.c2 * B.c5 + A.c4 * B.c10 -
                      A.c5 * B.c2 + A.c6 * B.c1 + A.c10 * B.c4 + A.c13 * B.c0;
    ctype const c14 = A.c0 * B.c14 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c10 +
                      A.c8 * B.c1 + A.c9 * B.c2 + A.c10 * B.c3 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 -
                      A.c4 * B.c14 + A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c10 +
                      A.c8 * B.c5 + A.c9 * B.c6 + A.c10 * B.c7 - A.c11 * B.c1 -
                      A.c12 * B.c2 - A.c13 * B.c3 + A.c14 * B.c4 + A.c15 * B.c0;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& A,
                                                MVec4ds_E<U> const& B)
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
    ctype const c11 = A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c4 + A.c11 * B.c0;
    ctype const c12 = -A.c1 * B.c3 + A.c3 * B.c1 + A.c4 * B.c5 + A.c12 * B.c0;
    ctype const c13 = A.c1 * B.c2 - A.c2 * B.c1 + A.c4 * B.c6 + A.c13 * B.c0;
    ctype const c14 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c7 + A.c5 * B.c4 + A.c6 * B.c5 + A.c7 * B.c6 +
                      A.c8 * B.c1 + A.c9 * B.c2 + A.c10 * B.c3 + A.c15 * B.c0;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& A,
                                                MVec4ds<U> const& B)
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
    ctype const c11 = A.c0 * B.c11 - A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c4;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c3 - A.c3 * B.c1 + A.c5 * B.c4;
    ctype const c13 = A.c0 * B.c13 - A.c1 * B.c2 + A.c2 * B.c1 + A.c6 * B.c4;
    ctype const c14 = A.c0 * B.c14 + A.c4 * B.c1 + A.c5 * B.c2 + A.c6 * B.c3;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c10 +
                      A.c4 * B.c5 + A.c5 * B.c6 + A.c6 * B.c7 + A.c7 * B.c0;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& A,
                                                MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c0 * B.c1;
    ctype const c3 = A.c0 * B.c2;
    ctype const c4 = A.c0 * B.c3;
    ctype const c5 = A.c1 * B.c3 - A.c4 * B.c0;
    ctype const c6 = A.c2 * B.c3 - A.c4 * B.c1;
    ctype const c7 = A.c3 * B.c3 - A.c4 * B.c2;
    ctype const c8 = A.c2 * B.c2 - A.c3 * B.c1;
    ctype const c9 = -A.c1 * B.c2 + A.c3 * B.c0;
    ctype const c10 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype const c11 = A.c0 * B.c4 - A.c6 * B.c2 + A.c7 * B.c1 + A.c8 * B.c3;
    ctype const c12 = A.c0 * B.c5 + A.c5 * B.c2 - A.c7 * B.c0 + A.c9 * B.c3;
    ctype const c13 = A.c0 * B.c6 - A.c5 * B.c1 + A.c6 * B.c0 + A.c10 * B.c3;
    ctype const c14 = A.c0 * B.c7 + A.c8 * B.c0 + A.c9 * B.c1 + A.c10 * B.c2;
    ctype const c15 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 - A.c4 * B.c7 -
                      A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c2 + A.c14 * B.c3;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds_U<T> const& A,
                                                MVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c1 * B.c0;
    ctype const c3 = A.c2 * B.c0;
    ctype const c4 = A.c3 * B.c0;
    ctype const c5 = A.c0 * B.c4 - A.c3 * B.c1;
    ctype const c6 = A.c1 * B.c4 - A.c3 * B.c2;
    ctype const c7 = A.c2 * B.c4 - A.c3 * B.c3;
    ctype const c8 = A.c1 * B.c3 - A.c2 * B.c2;
    ctype const c9 = -A.c0 * B.c3 + A.c2 * B.c1;
    ctype const c10 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c11 = A.c1 * B.c7 - A.c2 * B.c6 + A.c3 * B.c8 + A.c4 * B.c0;
    ctype const c12 = -A.c0 * B.c7 + A.c2 * B.c5 + A.c3 * B.c9 + A.c5 * B.c0;
    ctype const c13 = A.c0 * B.c6 - A.c1 * B.c5 + A.c3 * B.c10 + A.c6 * B.c0;
    ctype const c14 = A.c0 * B.c8 + A.c1 * B.c9 + A.c2 * B.c10 + A.c7 * B.c0;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 - A.c3 * B.c14 -
                      A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c4;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& M, PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(M.c0 * ctype(ps));
}

// sta4ds wdg :: wdg(ps,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(PScalar4ds<T> ps, MVec4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(ctype(ps) * M.c0);
}

// sta4ds wdg :: wdg(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& M,
                                                TriVec4ds<U> const& t)
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
    ctype const c15 = M.c1 * t.x + M.c2 * t.y + M.c3 * t.z - M.c4 * t.w;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(TriVec4ds<T> const& t,
                                                MVec4ds<U> const& M)
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
    ctype const c15 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3 + t.w * M.c4;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& M, BiVec4ds<U> const& B)
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
    ctype const c11 = M.c2 * B.vz - M.c3 * B.vy + M.c4 * B.mx;
    ctype const c12 = -M.c1 * B.vz + M.c3 * B.vx + M.c4 * B.my;
    ctype const c13 = M.c1 * B.vy - M.c2 * B.vx + M.c4 * B.mz;
    ctype const c14 = M.c1 * B.mx + M.c2 * B.my + M.c3 * B.mz;
    ctype const c15 = M.c5 * B.mx + M.c6 * B.my + M.c7 * B.mz + M.c8 * B.vx +
                      M.c9 * B.vy + M.c10 * B.vz;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(BiVec4ds<T> const& B, MVec4ds<U> const& M)
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
    ctype const c11 = -B.vy * M.c3 + B.vz * M.c2 + B.mx * M.c4;
    ctype const c12 = B.vx * M.c3 - B.vz * M.c1 + B.my * M.c4;
    ctype const c13 = -B.vx * M.c2 + B.vy * M.c1 + B.mz * M.c4;
    ctype const c14 = B.mx * M.c1 + B.my * M.c2 + B.mz * M.c3;
    ctype const c15 = B.vx * M.c8 + B.vy * M.c9 + B.vz * M.c10 + B.mx * M.c5 +
                      B.my * M.c6 + B.mz * M.c7;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& M, Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = M.c0 * v.w;
    ctype const c5 = M.c1 * v.w - M.c4 * v.x;
    ctype const c6 = M.c2 * v.w - M.c4 * v.y;
    ctype const c7 = M.c3 * v.w - M.c4 * v.z;
    ctype const c8 = M.c2 * v.z - M.c3 * v.y;
    ctype const c9 = -M.c1 * v.z + M.c3 * v.x;
    ctype const c10 = M.c1 * v.y - M.c2 * v.x;
    ctype const c11 = -M.c6 * v.z + M.c7 * v.y + M.c8 * v.w;
    ctype const c12 = M.c5 * v.z - M.c7 * v.x + M.c9 * v.w;
    ctype const c13 = -M.c5 * v.y + M.c6 * v.x + M.c10 * v.w;
    ctype const c14 = M.c8 * v.x + M.c9 * v.y + M.c10 * v.z;
    ctype const c15 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z + M.c14 * v.w;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(Vec4ds<T> const& v, MVec4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = v.w * M.c0;
    ctype const c5 = v.x * M.c4 - v.w * M.c1;
    ctype const c6 = v.y * M.c4 - v.w * M.c2;
    ctype const c7 = v.z * M.c4 - v.w * M.c3;
    ctype const c8 = v.y * M.c3 - v.z * M.c2;
    ctype const c9 = -v.x * M.c3 + v.z * M.c1;
    ctype const c10 = v.x * M.c2 - v.y * M.c1;
    ctype const c11 = v.y * M.c7 - v.z * M.c6 + v.w * M.c8;
    ctype const c12 = -v.x * M.c7 + v.z * M.c5 + v.w * M.c9;
    ctype const c13 = v.x * M.c6 - v.y * M.c5 + v.w * M.c10;
    ctype const c14 = v.x * M.c8 + v.y * M.c9 + v.z * M.c10;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 - v.w * M.c14;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(MVec4ds<T> const& M, Scalar4ds<U> s)
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
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> wdg(Scalar4ds<T> s, MVec4ds<U> const& M)
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
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds wdg :: wdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& A,
                                                  MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
                     A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& A,
                                                  MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4 - A.c2 * B.c2 + A.c3 * B.c1 + A.c4 * B.c3;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c2 - A.c3 * B.c0 + A.c5 * B.c3;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c1 + A.c2 * B.c0 + A.c6 * B.c3;
    ctype const c7 = A.c0 * B.c7 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds wdg :: wdg(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> wdg(MVec4ds_U<T> const& A,
                                                  MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c4 + A.c4 * B.c0;
    ctype const c5 = -A.c0 * B.c3 + A.c2 * B.c1 + A.c3 * B.c5 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c2 - A.c1 * B.c1 + A.c3 * B.c6 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c7 * B.c0;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds wdg :: wdg(mv_e,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& M,
                                                   PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(M.c0 * ctype(ps));
}

// sta4ds wdg :: wdg(ps,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(PScalar4ds<T> ps,
                                                   MVec4ds_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(ctype(ps) * M.c0);
}

// sta4ds wdg :: wdg(mv_e,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& M,
                                                  TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * t.x;
    ctype const c1 = M.c0 * t.y;
    ctype const c2 = M.c0 * t.z;
    ctype const c3 = M.c0 * t.w;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(trivec,mv_e) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(TriVec4ds<T> const& t,
                                                  MVec4ds_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * M.c0;
    ctype const c1 = t.y * M.c0;
    ctype const c2 = t.z * M.c0;
    ctype const c3 = t.w * M.c0;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& M,
                                                  BiVec4ds<U> const& B)
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
        M.c1 * B.mx + M.c2 * B.my + M.c3 * B.mz + M.c4 * B.vx + M.c5 * B.vy + M.c6 * B.vz;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(BiVec4ds<T> const& B,
                                                  MVec4ds_E<U> const& M)
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
        B.vx * M.c4 + B.vy * M.c5 + B.vz * M.c6 + B.mx * M.c1 + B.my * M.c2 + B.mz * M.c3;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& M,
                                                  Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    ctype const c3 = M.c0 * v.w;
    ctype const c4 = -M.c2 * v.z + M.c3 * v.y + M.c4 * v.w;
    ctype const c5 = M.c1 * v.z - M.c3 * v.x + M.c5 * v.w;
    ctype const c6 = -M.c1 * v.y + M.c2 * v.x + M.c6 * v.w;
    ctype const c7 = M.c4 * v.x + M.c5 * v.y + M.c6 * v.z;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds wdg :: wdg(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> wdg(Vec4ds<T> const& v,
                                                  MVec4ds_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    ctype const c3 = v.w * M.c0;
    ctype const c4 = v.y * M.c3 - v.z * M.c2 + v.w * M.c4;
    ctype const c5 = -v.x * M.c3 + v.z * M.c1 + v.w * M.c5;
    ctype const c6 = v.x * M.c2 - v.y * M.c1 + v.w * M.c6;
    ctype const c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds wdg :: wdg(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(MVec4ds_E<T> const& M, Scalar4ds<U> s)
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
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(Scalar4ds<T> s, MVec4ds_E<U> const& M)
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
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(MVec4ds_U<T> const& A,
                                                  MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c3 - A.c3 * B.c0;
    ctype const c2 = A.c1 * B.c3 - A.c3 * B.c1;
    ctype const c3 = A.c2 * B.c3 - A.c3 * B.c2;
    ctype const c4 = A.c1 * B.c2 - A.c2 * B.c1;
    ctype const c5 = -A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c6 = A.c0 * B.c1 - A.c1 * B.c0;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 - A.c3 * B.c7 - A.c4 * B.c0 -
                     A.c5 * B.c1 - A.c6 * B.c2 + A.c7 * B.c3;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(mv_u,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] MVec4ds_U<T> const&,
                                                  [[maybe_unused]] PScalar4ds<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(ps,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar4ds<T>,
                                                  [[maybe_unused]] MVec4ds_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(mv_u,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(MVec4ds_U<T> const& M,
                                                   TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(M.c0 * t.x + M.c1 * t.y + M.c2 * t.z - M.c3 * t.w);
}

// sta4ds wdg :: wdg(trivec,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(TriVec4ds<T> const& t,
                                                   MVec4ds_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(-t.x * M.c0 - t.y * M.c1 - t.z * M.c2 + t.w * M.c3);
}

// sta4ds wdg :: wdg(mv_u,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(MVec4ds_U<T> const& M,
                                                  BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.vz - M.c2 * B.vy + M.c3 * B.mx;
    ctype const c1 = -M.c0 * B.vz + M.c2 * B.vx + M.c3 * B.my;
    ctype const c2 = M.c0 * B.vy - M.c1 * B.vx + M.c3 * B.mz;
    ctype const c3 = M.c0 * B.mx + M.c1 * B.my + M.c2 * B.mz;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(bivec,mv_u) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(BiVec4ds<T> const& B,
                                                  MVec4ds_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * M.c2 + B.vz * M.c1 + B.mx * M.c3;
    ctype const c1 = B.vx * M.c2 - B.vz * M.c0 + B.my * M.c3;
    ctype const c2 = -B.vx * M.c1 + B.vy * M.c0 + B.mz * M.c3;
    ctype const c3 = B.mx * M.c0 + B.my * M.c1 + B.mz * M.c2;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(MVec4ds_U<T> const& M,
                                                  Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.w - M.c3 * v.x;
    ctype const c2 = M.c1 * v.w - M.c3 * v.y;
    ctype const c3 = M.c2 * v.w - M.c3 * v.z;
    ctype const c4 = M.c1 * v.z - M.c2 * v.y;
    ctype const c5 = -M.c0 * v.z + M.c2 * v.x;
    ctype const c6 = M.c0 * v.y - M.c1 * v.x;
    ctype const c7 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z + M.c7 * v.w;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> wdg(Vec4ds<T> const& v,
                                                  MVec4ds_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c3 - v.w * M.c0;
    ctype const c2 = v.y * M.c3 - v.w * M.c1;
    ctype const c3 = v.z * M.c3 - v.w * M.c2;
    ctype const c4 = v.y * M.c2 - v.z * M.c1;
    ctype const c5 = -v.x * M.c2 + v.z * M.c0;
    ctype const c6 = v.x * M.c1 - v.y * M.c0;
    ctype const c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6 - v.w * M.c7;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds wdg :: wdg(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> wdg(MVec4ds_U<T> const& M, Scalar4ds<U> s)
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
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds wdg :: wdg(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> wdg(Scalar4ds<T> s, MVec4ds_U<U> const& M)
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
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds wdg :: wdg(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar4ds<T>,
                                                  [[maybe_unused]] PScalar4ds<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar4ds<T>,
                                                  [[maybe_unused]] TriVec4ds<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec4ds<T> const&,
                                                  [[maybe_unused]] PScalar4ds<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar4ds<T>,
                                                  [[maybe_unused]] BiVec4ds<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec4ds<T> const&,
                                                  [[maybe_unused]] PScalar4ds<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar4ds<T>,
                                                  [[maybe_unused]] Vec4ds<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec4ds<T> const&,
                                                  [[maybe_unused]] PScalar4ds<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(PScalar4ds<T> ps, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(ctype(ps) * ctype(s));
}

// sta4ds wdg :: wdg(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(Scalar4ds<T> s, PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(ctype(s) * ctype(ps));
}

// sta4ds wdg :: wdg(trivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec4ds<T> const&,
                                                  [[maybe_unused]] TriVec4ds<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(trivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec4ds<T> const&,
                                                  [[maybe_unused]] BiVec4ds<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(bivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec4ds<T> const&,
                                                  [[maybe_unused]] TriVec4ds<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(0.0);
}

// sta4ds wdg :: wdg(trivec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(TriVec4ds<T> const& t,
                                                   Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z + t.w * v.w);
}

// sta4ds wdg :: wdg(vec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(Vec4ds<T> const& v,
                                                   TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(v.x * t.x + v.y * t.y + v.z * t.z - v.w * t.w);
}

// sta4ds wdg :: wdg(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(TriVec4ds<T> const& t, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(Scalar4ds<T> s, TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(bivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> wdg(BiVec4ds<T> const& B1,
                                                   BiVec4ds<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(B1.vx * B2.mx + B1.vy * B2.my + B1.vz * B2.mz +
                             B1.mx * B2.vx + B1.my * B2.vy + B1.mz * B2.vz);
}

// sta4ds wdg :: wdg(bivec,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(BiVec4ds<T> const& B,
                                                  Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vy * v.z + B.vz * v.y + B.mx * v.w;
    ctype const c1 = B.vx * v.z - B.vz * v.x + B.my * v.w;
    ctype const c2 = -B.vx * v.y + B.vy * v.x + B.mz * v.w;
    ctype const c3 = B.mx * v.x + B.my * v.y + B.mz * v.z;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(vec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> wdg(Vec4ds<T> const& v,
                                                  BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.vz - v.z * B.vy + v.w * B.mx;
    ctype const c1 = -v.x * B.vz + v.z * B.vx + v.w * B.my;
    ctype const c2 = v.x * B.vy - v.y * B.vx + v.w * B.mz;
    ctype const c3 = v.x * B.mx + v.y * B.my + v.z * B.mz;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> wdg(BiVec4ds<T> const& B, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec4ds<ctype>(c0, c1, c2, c3, c4, c5);
}

// sta4ds wdg :: wdg(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> wdg(Scalar4ds<T> s, BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec4ds<ctype>(c0, c1, c2, c3, c4, c5);
}

// sta4ds wdg :: wdg(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> wdg(Vec4ds<T> const& v1, Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.x * v2.w - v1.w * v2.x;
    ctype const c1 = v1.y * v2.w - v1.w * v2.y;
    ctype const c2 = v1.z * v2.w - v1.w * v2.z;
    ctype const c3 = v1.y * v2.z - v1.z * v2.y;
    ctype const c4 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c5 = v1.x * v2.y - v1.y * v2.x;
    return BiVec4ds<ctype>(c0, c1, c2, c3, c4, c5);
}

// sta4ds wdg :: wdg(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> wdg(Vec4ds<T> const& v, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> wdg(Scalar4ds<T> s, Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds wdg :: wdg(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> wdg(Scalar4ds<T> s1, Scalar4ds<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// Implements the left contraction as per "PGA Illuminated", E. Lengyel:
//
// operator<<(a,b) = l_contract(a,b) = rwdg( l_dual(a), b )
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The resulting object lies in A and is perpendicular to B
//
// Implements the right bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator>>(a,b) = r_contract(a,b) = rwdg( a, r_dual(b) )
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// sta4ds gpr :: gpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& A,
                                                      MVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c4 * B.c4 +
                     A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7 - A.c8 * B.c8 - A.c9 * B.c9 -
                     A.c10 * B.c10 - A.c11 * B.c11 - A.c12 * B.c12 - A.c13 * B.c13 +
                     A.c14 * B.c14 - A.c15 * B.c15;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 + A.c2 * B.c10 - A.c3 * B.c9 -
                     A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c13 - A.c7 * B.c12 -
                     A.c8 * B.c14 + A.c9 * B.c3 - A.c10 * B.c2 + A.c11 * B.c15 -
                     A.c12 * B.c7 + A.c13 * B.c6 - A.c14 * B.c8 - A.c15 * B.c11;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c10 + A.c2 * B.c0 + A.c3 * B.c8 -
                     A.c4 * B.c6 - A.c5 * B.c13 + A.c6 * B.c4 + A.c7 * B.c11 -
                     A.c8 * B.c3 - A.c9 * B.c14 + A.c10 * B.c1 + A.c11 * B.c7 +
                     A.c12 * B.c15 - A.c13 * B.c5 - A.c14 * B.c9 - A.c15 * B.c12;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c0 - A.c4 * B.c7 +
                     A.c5 * B.c12 - A.c6 * B.c11 + A.c7 * B.c4 + A.c8 * B.c2 -
                     A.c9 * B.c1 - A.c10 * B.c14 - A.c11 * B.c6 + A.c12 * B.c5 +
                     A.c13 * B.c15 - A.c14 * B.c10 - A.c15 * B.c13;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 + A.c4 * B.c0 +
                     A.c5 * B.c1 + A.c6 * B.c2 + A.c7 * B.c3 - A.c8 * B.c11 -
                     A.c9 * B.c12 - A.c10 * B.c13 - A.c11 * B.c8 - A.c12 * B.c9 -
                     A.c13 * B.c10 + A.c14 * B.c15 - A.c15 * B.c14;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c13 - A.c3 * B.c12 -
                     A.c4 * B.c1 + A.c5 * B.c0 + A.c6 * B.c10 - A.c7 * B.c9 -
                     A.c8 * B.c15 + A.c9 * B.c7 - A.c10 * B.c6 + A.c11 * B.c14 -
                     A.c12 * B.c3 + A.c13 * B.c2 - A.c14 * B.c11 - A.c15 * B.c8;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c13 + A.c2 * B.c4 + A.c3 * B.c11 -
                     A.c4 * B.c2 - A.c5 * B.c10 + A.c6 * B.c0 + A.c7 * B.c8 -
                     A.c8 * B.c7 - A.c9 * B.c15 + A.c10 * B.c5 + A.c11 * B.c3 +
                     A.c12 * B.c14 - A.c13 * B.c1 - A.c14 * B.c12 - A.c15 * B.c9;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c12 - A.c2 * B.c11 + A.c3 * B.c4 -
                     A.c4 * B.c3 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c0 + A.c8 * B.c6 -
                     A.c9 * B.c5 - A.c10 * B.c15 - A.c11 * B.c2 + A.c12 * B.c1 +
                     A.c13 * B.c14 - A.c14 * B.c13 - A.c15 * B.c10;
    ctype const c8 = A.c0 * B.c8 - A.c1 * B.c14 + A.c2 * B.c3 - A.c3 * B.c2 +
                     A.c4 * B.c11 + A.c5 * B.c15 - A.c6 * B.c7 + A.c7 * B.c6 +
                     A.c8 * B.c0 + A.c9 * B.c10 - A.c10 * B.c9 + A.c11 * B.c4 +
                     A.c12 * B.c13 - A.c13 * B.c12 - A.c14 * B.c1 + A.c15 * B.c5;
    ctype const c9 = A.c0 * B.c9 - A.c1 * B.c3 - A.c2 * B.c14 + A.c3 * B.c1 +
                     A.c4 * B.c12 + A.c5 * B.c7 + A.c6 * B.c15 - A.c7 * B.c5 -
                     A.c8 * B.c10 + A.c9 * B.c0 + A.c10 * B.c8 - A.c11 * B.c13 +
                     A.c12 * B.c4 + A.c13 * B.c11 - A.c14 * B.c2 + A.c15 * B.c6;
    ctype const c10 = A.c0 * B.c10 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c14 +
                      A.c4 * B.c13 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c15 +
                      A.c8 * B.c9 - A.c9 * B.c8 + A.c10 * B.c0 + A.c11 * B.c12 -
                      A.c12 * B.c11 + A.c13 * B.c4 - A.c14 * B.c3 + A.c15 * B.c7;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c15 + A.c2 * B.c7 - A.c3 * B.c6 +
                      A.c4 * B.c8 + A.c5 * B.c14 - A.c6 * B.c3 + A.c7 * B.c2 +
                      A.c8 * B.c4 + A.c9 * B.c13 - A.c10 * B.c12 + A.c11 * B.c0 +
                      A.c12 * B.c10 - A.c13 * B.c9 - A.c14 * B.c5 + A.c15 * B.c1;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c7 - A.c2 * B.c15 + A.c3 * B.c5 +
                      A.c4 * B.c9 + A.c5 * B.c3 + A.c6 * B.c14 - A.c7 * B.c1 -
                      A.c8 * B.c13 + A.c9 * B.c4 + A.c10 * B.c11 - A.c11 * B.c10 +
                      A.c12 * B.c0 + A.c13 * B.c8 - A.c14 * B.c6 + A.c15 * B.c2;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c15 +
                      A.c4 * B.c10 - A.c5 * B.c2 + A.c6 * B.c1 + A.c7 * B.c14 +
                      A.c8 * B.c12 - A.c9 * B.c11 + A.c10 * B.c4 + A.c11 * B.c9 -
                      A.c12 * B.c8 + A.c13 * B.c0 - A.c14 * B.c7 + A.c15 * B.c3;
    ctype const c14 = A.c0 * B.c14 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c10 -
                      A.c4 * B.c15 + A.c5 * B.c11 + A.c6 * B.c12 + A.c7 * B.c13 +
                      A.c8 * B.c1 + A.c9 * B.c2 + A.c10 * B.c3 - A.c11 * B.c5 -
                      A.c12 * B.c6 - A.c13 * B.c7 + A.c14 * B.c0 + A.c15 * B.c4;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 -
                      A.c4 * B.c14 + A.c5 * B.c8 + A.c6 * B.c9 + A.c7 * B.c10 +
                      A.c8 * B.c5 + A.c9 * B.c6 + A.c10 * B.c7 - A.c11 * B.c1 -
                      A.c12 * B.c2 - A.c13 * B.c3 + A.c14 * B.c4 + A.c15 * B.c0;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& A,
                                                      MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2 + A.c7 * B.c3 - A.c8 * B.c4 -
                     A.c9 * B.c5 - A.c10 * B.c6 - A.c15 * B.c7;
    ctype const c1 = A.c1 * B.c0 + A.c2 * B.c6 - A.c3 * B.c5 - A.c4 * B.c1 +
                     A.c11 * B.c7 - A.c12 * B.c3 + A.c13 * B.c2 - A.c14 * B.c4;
    ctype const c2 = -A.c1 * B.c6 + A.c2 * B.c0 + A.c3 * B.c4 - A.c4 * B.c2 +
                     A.c11 * B.c3 + A.c12 * B.c7 - A.c13 * B.c1 - A.c14 * B.c5;
    ctype const c3 = A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c3 -
                     A.c11 * B.c2 + A.c12 * B.c1 + A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c4 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c4 * B.c0 -
                     A.c11 * B.c4 - A.c12 * B.c5 - A.c13 * B.c6 + A.c14 * B.c7;
    ctype const c5 = A.c0 * B.c1 + A.c5 * B.c0 + A.c6 * B.c6 - A.c7 * B.c5 - A.c8 * B.c7 +
                     A.c9 * B.c3 - A.c10 * B.c2 - A.c15 * B.c4;
    ctype const c6 = A.c0 * B.c2 - A.c5 * B.c6 + A.c6 * B.c0 + A.c7 * B.c4 - A.c8 * B.c3 -
                     A.c9 * B.c7 + A.c10 * B.c1 - A.c15 * B.c5;
    ctype const c7 = A.c0 * B.c3 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c0 + A.c8 * B.c2 -
                     A.c9 * B.c1 - A.c10 * B.c7 - A.c15 * B.c6;
    ctype const c8 = A.c0 * B.c4 + A.c5 * B.c7 - A.c6 * B.c3 + A.c7 * B.c2 + A.c8 * B.c0 +
                     A.c9 * B.c6 - A.c10 * B.c5 + A.c15 * B.c1;
    ctype const c9 = A.c0 * B.c5 + A.c5 * B.c3 + A.c6 * B.c7 - A.c7 * B.c1 - A.c8 * B.c6 +
                     A.c9 * B.c0 + A.c10 * B.c4 + A.c15 * B.c2;
    ctype const c10 = A.c0 * B.c6 - A.c5 * B.c2 + A.c6 * B.c1 + A.c7 * B.c7 +
                      A.c8 * B.c5 - A.c9 * B.c4 + A.c10 * B.c0 + A.c15 * B.c3;
    ctype const c11 = -A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c4 +
                      A.c11 * B.c0 + A.c12 * B.c6 - A.c13 * B.c5 - A.c14 * B.c1;
    ctype const c12 = -A.c1 * B.c3 - A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c5 -
                      A.c11 * B.c6 + A.c12 * B.c0 + A.c13 * B.c4 - A.c14 * B.c2;
    ctype const c13 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c7 + A.c4 * B.c6 +
                      A.c11 * B.c5 - A.c12 * B.c4 + A.c13 * B.c0 - A.c14 * B.c3;
    ctype const c14 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 - A.c4 * B.c7 -
                      A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c7 + A.c5 * B.c4 + A.c6 * B.c5 + A.c7 * B.c6 +
                      A.c8 * B.c1 + A.c9 * B.c2 + A.c10 * B.c3 + A.c15 * B.c0;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& A,
                                                      MVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 - A.c4 * B.c8 -
                     A.c5 * B.c9 - A.c6 * B.c10 - A.c7 * B.c15;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c4 + A.c2 * B.c13 - A.c3 * B.c12 -
                     A.c4 * B.c14 + A.c5 * B.c3 - A.c6 * B.c2 - A.c7 * B.c11;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c13 + A.c2 * B.c4 + A.c3 * B.c11 -
                     A.c4 * B.c3 - A.c5 * B.c14 + A.c6 * B.c1 - A.c7 * B.c12;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c12 - A.c2 * B.c11 + A.c3 * B.c4 +
                     A.c4 * B.c2 - A.c5 * B.c1 - A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 -
                     A.c4 * B.c11 - A.c5 * B.c12 - A.c6 * B.c13 - A.c7 * B.c14;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c0 + A.c2 * B.c10 - A.c3 * B.c9 -
                     A.c4 * B.c15 + A.c5 * B.c7 - A.c6 * B.c6 - A.c7 * B.c8;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c10 + A.c2 * B.c0 + A.c3 * B.c8 -
                     A.c4 * B.c7 - A.c5 * B.c15 + A.c6 * B.c5 - A.c7 * B.c9;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c0 + A.c4 * B.c6 -
                     A.c5 * B.c5 - A.c6 * B.c15 - A.c7 * B.c10;
    ctype const c8 = A.c0 * B.c8 + A.c1 * B.c15 - A.c2 * B.c7 + A.c3 * B.c6 +
                     A.c4 * B.c0 + A.c5 * B.c10 - A.c6 * B.c9 + A.c7 * B.c5;
    ctype const c9 = A.c0 * B.c9 + A.c1 * B.c7 + A.c2 * B.c15 - A.c3 * B.c5 -
                     A.c4 * B.c10 + A.c5 * B.c0 + A.c6 * B.c8 + A.c7 * B.c6;
    ctype const c10 = A.c0 * B.c10 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c15 +
                      A.c4 * B.c9 - A.c5 * B.c8 + A.c6 * B.c0 + A.c7 * B.c7;
    ctype const c11 = A.c0 * B.c11 + A.c1 * B.c14 - A.c2 * B.c3 + A.c3 * B.c2 +
                      A.c4 * B.c4 + A.c5 * B.c13 - A.c6 * B.c12 + A.c7 * B.c1;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c3 + A.c2 * B.c14 - A.c3 * B.c1 -
                      A.c4 * B.c13 + A.c5 * B.c4 + A.c6 * B.c11 + A.c7 * B.c2;
    ctype const c13 = A.c0 * B.c13 - A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c14 +
                      A.c4 * B.c12 - A.c5 * B.c11 + A.c6 * B.c4 + A.c7 * B.c3;
    ctype const c14 = A.c0 * B.c14 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 +
                      A.c4 * B.c1 + A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c4;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c8 + A.c2 * B.c9 + A.c3 * B.c10 +
                      A.c4 * B.c5 + A.c5 * B.c6 + A.c6 * B.c7 + A.c7 * B.c0;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& A,
                                                      MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2 + A.c4 * B.c3 -
                     A.c11 * B.c4 - A.c12 * B.c5 - A.c13 * B.c6 + A.c14 * B.c7;
    ctype const c1 = A.c0 * B.c0 + A.c5 * B.c3 + A.c6 * B.c6 - A.c7 * B.c5 - A.c8 * B.c7 +
                     A.c9 * B.c2 - A.c10 * B.c1 - A.c15 * B.c4;
    ctype const c2 = A.c0 * B.c1 - A.c5 * B.c6 + A.c6 * B.c3 + A.c7 * B.c4 - A.c8 * B.c2 -
                     A.c9 * B.c7 + A.c10 * B.c0 - A.c15 * B.c5;
    ctype const c3 = A.c0 * B.c2 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c3 + A.c8 * B.c1 -
                     A.c9 * B.c0 - A.c10 * B.c7 - A.c15 * B.c6;
    ctype const c4 = A.c0 * B.c3 + A.c5 * B.c0 + A.c6 * B.c1 + A.c7 * B.c2 - A.c8 * B.c4 -
                     A.c9 * B.c5 - A.c10 * B.c6 - A.c15 * B.c7;
    ctype const c5 = A.c1 * B.c3 + A.c2 * B.c6 - A.c3 * B.c5 - A.c4 * B.c0 +
                     A.c11 * B.c7 - A.c12 * B.c2 + A.c13 * B.c1 - A.c14 * B.c4;
    ctype const c6 = -A.c1 * B.c6 + A.c2 * B.c3 + A.c3 * B.c4 - A.c4 * B.c1 +
                     A.c11 * B.c2 + A.c12 * B.c7 - A.c13 * B.c0 - A.c14 * B.c5;
    ctype const c7 = A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c3 - A.c4 * B.c2 -
                     A.c11 * B.c1 + A.c12 * B.c0 + A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c8 = -A.c1 * B.c7 + A.c2 * B.c2 - A.c3 * B.c1 + A.c4 * B.c4 +
                     A.c11 * B.c3 + A.c12 * B.c6 - A.c13 * B.c5 - A.c14 * B.c0;
    ctype const c9 = -A.c1 * B.c2 - A.c2 * B.c7 + A.c3 * B.c0 + A.c4 * B.c5 -
                     A.c11 * B.c6 + A.c12 * B.c3 + A.c13 * B.c4 - A.c14 * B.c1;
    ctype const c10 = A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c7 + A.c4 * B.c6 +
                      A.c11 * B.c5 - A.c12 * B.c4 + A.c13 * B.c3 - A.c14 * B.c2;
    ctype const c11 = A.c0 * B.c4 + A.c5 * B.c7 - A.c6 * B.c2 + A.c7 * B.c1 +
                      A.c8 * B.c3 + A.c9 * B.c6 - A.c10 * B.c5 + A.c15 * B.c0;
    ctype const c12 = A.c0 * B.c5 + A.c5 * B.c2 + A.c6 * B.c7 - A.c7 * B.c0 -
                      A.c8 * B.c6 + A.c9 * B.c3 + A.c10 * B.c4 + A.c15 * B.c1;
    ctype const c13 = A.c0 * B.c6 - A.c5 * B.c1 + A.c6 * B.c0 + A.c7 * B.c7 +
                      A.c8 * B.c5 - A.c9 * B.c4 + A.c10 * B.c3 + A.c15 * B.c2;
    ctype const c14 = A.c0 * B.c7 + A.c5 * B.c4 + A.c6 * B.c5 + A.c7 * B.c6 +
                      A.c8 * B.c0 + A.c9 * B.c1 + A.c10 * B.c2 + A.c15 * B.c3;
    ctype const c15 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 - A.c4 * B.c7 -
                      A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c2 + A.c14 * B.c3;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& A,
                                                      MVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c4 -
                     A.c4 * B.c11 - A.c5 * B.c12 - A.c6 * B.c13 + A.c7 * B.c14;
    ctype const c1 = A.c0 * B.c0 + A.c1 * B.c10 - A.c2 * B.c9 - A.c3 * B.c5 +
                     A.c4 * B.c15 - A.c5 * B.c7 + A.c6 * B.c6 - A.c7 * B.c8;
    ctype const c2 = -A.c0 * B.c10 + A.c1 * B.c0 + A.c2 * B.c8 - A.c3 * B.c6 +
                     A.c4 * B.c7 + A.c5 * B.c15 - A.c6 * B.c5 - A.c7 * B.c9;
    ctype const c3 = A.c0 * B.c9 - A.c1 * B.c8 + A.c2 * B.c0 - A.c3 * B.c7 - A.c4 * B.c6 +
                     A.c5 * B.c5 + A.c6 * B.c15 - A.c7 * B.c10;
    ctype const c4 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c7 + A.c3 * B.c0 -
                     A.c4 * B.c8 - A.c5 * B.c9 - A.c6 * B.c10 + A.c7 * B.c15;
    ctype const c5 = A.c0 * B.c4 + A.c1 * B.c13 - A.c2 * B.c12 - A.c3 * B.c1 +
                     A.c4 * B.c14 - A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c11;
    ctype const c6 = -A.c0 * B.c13 + A.c1 * B.c4 + A.c2 * B.c11 - A.c3 * B.c2 +
                     A.c4 * B.c3 + A.c5 * B.c14 - A.c6 * B.c1 - A.c7 * B.c12;
    ctype const c7 = A.c0 * B.c12 - A.c1 * B.c11 + A.c2 * B.c4 - A.c3 * B.c3 -
                     A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c8 = -A.c0 * B.c14 + A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c11 +
                     A.c4 * B.c4 + A.c5 * B.c13 - A.c6 * B.c12 - A.c7 * B.c1;
    ctype const c9 = -A.c0 * B.c3 - A.c1 * B.c14 + A.c2 * B.c1 + A.c3 * B.c12 -
                     A.c4 * B.c13 + A.c5 * B.c4 + A.c6 * B.c11 - A.c7 * B.c2;
    ctype const c10 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c14 + A.c3 * B.c13 +
                      A.c4 * B.c12 - A.c5 * B.c11 + A.c6 * B.c4 - A.c7 * B.c3;
    ctype const c11 = -A.c0 * B.c15 + A.c1 * B.c7 - A.c2 * B.c6 + A.c3 * B.c8 +
                      A.c4 * B.c0 + A.c5 * B.c10 - A.c6 * B.c9 - A.c7 * B.c5;
    ctype const c12 = -A.c0 * B.c7 - A.c1 * B.c15 + A.c2 * B.c5 + A.c3 * B.c9 -
                      A.c4 * B.c10 + A.c5 * B.c0 + A.c6 * B.c8 - A.c7 * B.c6;
    ctype const c13 = A.c0 * B.c6 - A.c1 * B.c5 - A.c2 * B.c15 + A.c3 * B.c10 +
                      A.c4 * B.c9 - A.c5 * B.c8 + A.c6 * B.c0 - A.c7 * B.c7;
    ctype const c14 = A.c0 * B.c8 + A.c1 * B.c9 + A.c2 * B.c10 - A.c3 * B.c15 -
                      A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 + A.c7 * B.c0;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 - A.c3 * B.c14 -
                      A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c4;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& M,
                                                      PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c15 * ctype(ps);
    ctype const c1 = M.c11 * ctype(ps);
    ctype const c2 = M.c12 * ctype(ps);
    ctype const c3 = M.c13 * ctype(ps);
    ctype const c4 = M.c14 * ctype(ps);
    ctype const c5 = -M.c8 * ctype(ps);
    ctype const c6 = -M.c9 * ctype(ps);
    ctype const c7 = -M.c10 * ctype(ps);
    ctype const c8 = M.c5 * ctype(ps);
    ctype const c9 = M.c6 * ctype(ps);
    ctype const c10 = M.c7 * ctype(ps);
    ctype const c11 = -M.c1 * ctype(ps);
    ctype const c12 = -M.c2 * ctype(ps);
    ctype const c13 = -M.c3 * ctype(ps);
    ctype const c14 = -M.c4 * ctype(ps);
    ctype const c15 = M.c0 * ctype(ps);
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps,
                                                      MVec4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c15;
    ctype const c1 = -ctype(ps) * M.c11;
    ctype const c2 = -ctype(ps) * M.c12;
    ctype const c3 = -ctype(ps) * M.c13;
    ctype const c4 = -ctype(ps) * M.c14;
    ctype const c5 = -ctype(ps) * M.c8;
    ctype const c6 = -ctype(ps) * M.c9;
    ctype const c7 = -ctype(ps) * M.c10;
    ctype const c8 = ctype(ps) * M.c5;
    ctype const c9 = ctype(ps) * M.c6;
    ctype const c10 = ctype(ps) * M.c7;
    ctype const c11 = ctype(ps) * M.c1;
    ctype const c12 = ctype(ps) * M.c2;
    ctype const c13 = ctype(ps) * M.c3;
    ctype const c14 = ctype(ps) * M.c4;
    ctype const c15 = ctype(ps) * M.c0;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& M,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * t.x - M.c12 * t.y - M.c13 * t.z + M.c14 * t.w;
    ctype const c1 = M.c6 * t.z - M.c7 * t.y - M.c8 * t.w - M.c15 * t.x;
    ctype const c2 = -M.c5 * t.z + M.c7 * t.x - M.c9 * t.w - M.c15 * t.y;
    ctype const c3 = M.c5 * t.y - M.c6 * t.x - M.c10 * t.w - M.c15 * t.z;
    ctype const c4 = -M.c8 * t.x - M.c9 * t.y - M.c10 * t.z - M.c15 * t.w;
    ctype const c5 = M.c2 * t.z - M.c3 * t.y + M.c11 * t.w - M.c14 * t.x;
    ctype const c6 = -M.c1 * t.z + M.c3 * t.x + M.c12 * t.w - M.c14 * t.y;
    ctype const c7 = M.c1 * t.y - M.c2 * t.x + M.c13 * t.w - M.c14 * t.z;
    ctype const c8 = -M.c1 * t.w + M.c4 * t.x + M.c12 * t.z - M.c13 * t.y;
    ctype const c9 = -M.c2 * t.w + M.c4 * t.y - M.c11 * t.z + M.c13 * t.x;
    ctype const c10 = -M.c3 * t.w + M.c4 * t.z + M.c11 * t.y - M.c12 * t.x;
    ctype const c11 = M.c0 * t.x + M.c5 * t.w + M.c9 * t.z - M.c10 * t.y;
    ctype const c12 = M.c0 * t.y + M.c6 * t.w - M.c8 * t.z + M.c10 * t.x;
    ctype const c13 = M.c0 * t.z + M.c7 * t.w + M.c8 * t.y - M.c9 * t.x;
    ctype const c14 = M.c0 * t.w + M.c5 * t.x + M.c6 * t.y + M.c7 * t.z;
    ctype const c15 = M.c1 * t.x + M.c2 * t.y + M.c3 * t.z - M.c4 * t.w;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t,
                                                      MVec4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c11 - t.y * M.c12 - t.z * M.c13 + t.w * M.c14;
    ctype const c1 = t.x * M.c15 - t.y * M.c7 + t.z * M.c6 - t.w * M.c8;
    ctype const c2 = t.x * M.c7 + t.y * M.c15 - t.z * M.c5 - t.w * M.c9;
    ctype const c3 = -t.x * M.c6 + t.y * M.c5 + t.z * M.c15 - t.w * M.c10;
    ctype const c4 = -t.x * M.c8 - t.y * M.c9 - t.z * M.c10 + t.w * M.c15;
    ctype const c5 = t.x * M.c14 - t.y * M.c3 + t.z * M.c2 - t.w * M.c11;
    ctype const c6 = t.x * M.c3 + t.y * M.c14 - t.z * M.c1 - t.w * M.c12;
    ctype const c7 = -t.x * M.c2 + t.y * M.c1 + t.z * M.c14 - t.w * M.c13;
    ctype const c8 = t.x * M.c4 + t.y * M.c13 - t.z * M.c12 - t.w * M.c1;
    ctype const c9 = -t.x * M.c13 + t.y * M.c4 + t.z * M.c11 - t.w * M.c2;
    ctype const c10 = t.x * M.c12 - t.y * M.c11 + t.z * M.c4 - t.w * M.c3;
    ctype const c11 = t.x * M.c0 + t.y * M.c10 - t.z * M.c9 - t.w * M.c5;
    ctype const c12 = -t.x * M.c10 + t.y * M.c0 + t.z * M.c8 - t.w * M.c6;
    ctype const c13 = t.x * M.c9 - t.y * M.c8 + t.z * M.c0 - t.w * M.c7;
    ctype const c14 = -t.x * M.c5 - t.y * M.c6 - t.z * M.c7 + t.w * M.c0;
    ctype const c15 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3 + t.w * M.c4;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& M,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c5 * B.vx + M.c6 * B.vy + M.c7 * B.vz - M.c8 * B.mx - M.c9 * B.my -
                     M.c10 * B.mz;
    ctype const c1 = M.c2 * B.mz - M.c3 * B.my - M.c4 * B.vx - M.c12 * B.vz +
                     M.c13 * B.vy - M.c14 * B.mx;
    ctype const c2 = -M.c1 * B.mz + M.c3 * B.mx - M.c4 * B.vy + M.c11 * B.vz -
                     M.c13 * B.vx - M.c14 * B.my;
    ctype const c3 = M.c1 * B.my - M.c2 * B.mx - M.c4 * B.vz - M.c11 * B.vy +
                     M.c12 * B.vx - M.c14 * B.mz;
    ctype const c4 = -M.c1 * B.vx - M.c2 * B.vy - M.c3 * B.vz - M.c11 * B.mx -
                     M.c12 * B.my - M.c13 * B.mz;
    ctype const c5 = M.c0 * B.vx + M.c6 * B.mz - M.c7 * B.my + M.c9 * B.vz -
                     M.c10 * B.vy - M.c15 * B.mx;
    ctype const c6 = M.c0 * B.vy - M.c5 * B.mz + M.c7 * B.mx - M.c8 * B.vz +
                     M.c10 * B.vx - M.c15 * B.my;
    ctype const c7 = M.c0 * B.vz + M.c5 * B.my - M.c6 * B.mx + M.c8 * B.vy - M.c9 * B.vx -
                     M.c15 * B.mz;
    ctype const c8 = M.c0 * B.mx - M.c6 * B.vz + M.c7 * B.vy + M.c9 * B.mz -
                     M.c10 * B.my + M.c15 * B.vx;
    ctype const c9 = M.c0 * B.my + M.c5 * B.vz - M.c7 * B.vx - M.c8 * B.mz +
                     M.c10 * B.mx + M.c15 * B.vy;
    ctype const c10 = M.c0 * B.mz - M.c5 * B.vy + M.c6 * B.vx + M.c8 * B.my -
                      M.c9 * B.mx + M.c15 * B.vz;
    ctype const c11 = M.c2 * B.vz - M.c3 * B.vy + M.c4 * B.mx + M.c12 * B.mz -
                      M.c13 * B.my - M.c14 * B.vx;
    ctype const c12 = -M.c1 * B.vz + M.c3 * B.vx + M.c4 * B.my - M.c11 * B.mz +
                      M.c13 * B.mx - M.c14 * B.vy;
    ctype const c13 = M.c1 * B.vy - M.c2 * B.vx + M.c4 * B.mz + M.c11 * B.my -
                      M.c12 * B.mx - M.c14 * B.vz;
    ctype const c14 = M.c1 * B.mx + M.c2 * B.my + M.c3 * B.mz - M.c11 * B.vx -
                      M.c12 * B.vy - M.c13 * B.vz;
    ctype const c15 = M.c5 * B.mx + M.c6 * B.my + M.c7 * B.mz + M.c8 * B.vx +
                      M.c9 * B.vy + M.c10 * B.vz;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B,
                                                      MVec4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c5 + B.vy * M.c6 + B.vz * M.c7 - B.mx * M.c8 - B.my * M.c9 -
                     B.mz * M.c10;
    ctype const c1 = B.vx * M.c4 + B.vy * M.c13 - B.vz * M.c12 - B.mx * M.c14 +
                     B.my * M.c3 - B.mz * M.c2;
    ctype const c2 = -B.vx * M.c13 + B.vy * M.c4 + B.vz * M.c11 - B.mx * M.c3 -
                     B.my * M.c14 + B.mz * M.c1;
    ctype const c3 = B.vx * M.c12 - B.vy * M.c11 + B.vz * M.c4 + B.mx * M.c2 -
                     B.my * M.c1 - B.mz * M.c14;
    ctype const c4 = B.vx * M.c1 + B.vy * M.c2 + B.vz * M.c3 - B.mx * M.c11 -
                     B.my * M.c12 - B.mz * M.c13;
    ctype const c5 = B.vx * M.c0 + B.vy * M.c10 - B.vz * M.c9 - B.mx * M.c15 +
                     B.my * M.c7 - B.mz * M.c6;
    ctype const c6 = -B.vx * M.c10 + B.vy * M.c0 + B.vz * M.c8 - B.mx * M.c7 -
                     B.my * M.c15 + B.mz * M.c5;
    ctype const c7 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c0 + B.mx * M.c6 - B.my * M.c5 -
                     B.mz * M.c15;
    ctype const c8 = B.vx * M.c15 - B.vy * M.c7 + B.vz * M.c6 + B.mx * M.c0 +
                     B.my * M.c10 - B.mz * M.c9;
    ctype const c9 = B.vx * M.c7 + B.vy * M.c15 - B.vz * M.c5 - B.mx * M.c10 +
                     B.my * M.c0 + B.mz * M.c8;
    ctype const c10 = -B.vx * M.c6 + B.vy * M.c5 + B.vz * M.c15 + B.mx * M.c9 -
                      B.my * M.c8 + B.mz * M.c0;
    ctype const c11 = B.vx * M.c14 - B.vy * M.c3 + B.vz * M.c2 + B.mx * M.c4 +
                      B.my * M.c13 - B.mz * M.c12;
    ctype const c12 = B.vx * M.c3 + B.vy * M.c14 - B.vz * M.c1 - B.mx * M.c13 +
                      B.my * M.c4 + B.mz * M.c11;
    ctype const c13 = -B.vx * M.c2 + B.vy * M.c1 + B.vz * M.c14 + B.mx * M.c12 -
                      B.my * M.c11 + B.mz * M.c4;
    ctype const c14 = B.vx * M.c11 + B.vy * M.c12 + B.vz * M.c13 + B.mx * M.c1 +
                      B.my * M.c2 + B.mz * M.c3;
    ctype const c15 = B.vx * M.c8 + B.vy * M.c9 + B.vz * M.c10 + B.mx * M.c5 +
                      B.my * M.c6 + B.mz * M.c7;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& M,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * v.x - M.c2 * v.y - M.c3 * v.z + M.c4 * v.w;
    ctype const c1 = M.c0 * v.x + M.c5 * v.w + M.c9 * v.z - M.c10 * v.y;
    ctype const c2 = M.c0 * v.y + M.c6 * v.w - M.c8 * v.z + M.c10 * v.x;
    ctype const c3 = M.c0 * v.z + M.c7 * v.w + M.c8 * v.y - M.c9 * v.x;
    ctype const c4 = M.c0 * v.w + M.c5 * v.x + M.c6 * v.y + M.c7 * v.z;
    ctype const c5 = M.c1 * v.w - M.c4 * v.x - M.c12 * v.z + M.c13 * v.y;
    ctype const c6 = M.c2 * v.w - M.c4 * v.y + M.c11 * v.z - M.c13 * v.x;
    ctype const c7 = M.c3 * v.w - M.c4 * v.z - M.c11 * v.y + M.c12 * v.x;
    ctype const c8 = M.c2 * v.z - M.c3 * v.y + M.c11 * v.w - M.c14 * v.x;
    ctype const c9 = -M.c1 * v.z + M.c3 * v.x + M.c12 * v.w - M.c14 * v.y;
    ctype const c10 = M.c1 * v.y - M.c2 * v.x + M.c13 * v.w - M.c14 * v.z;
    ctype const c11 = -M.c6 * v.z + M.c7 * v.y + M.c8 * v.w + M.c15 * v.x;
    ctype const c12 = M.c5 * v.z - M.c7 * v.x + M.c9 * v.w + M.c15 * v.y;
    ctype const c13 = -M.c5 * v.y + M.c6 * v.x + M.c10 * v.w + M.c15 * v.z;
    ctype const c14 = M.c8 * v.x + M.c9 * v.y + M.c10 * v.z + M.c15 * v.w;
    ctype const c15 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z + M.c14 * v.w;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v,
                                                      MVec4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * M.c1 - v.y * M.c2 - v.z * M.c3 + v.w * M.c4;
    ctype const c1 = v.x * M.c0 + v.y * M.c10 - v.z * M.c9 - v.w * M.c5;
    ctype const c2 = -v.x * M.c10 + v.y * M.c0 + v.z * M.c8 - v.w * M.c6;
    ctype const c3 = v.x * M.c9 - v.y * M.c8 + v.z * M.c0 - v.w * M.c7;
    ctype const c4 = -v.x * M.c5 - v.y * M.c6 - v.z * M.c7 + v.w * M.c0;
    ctype const c5 = v.x * M.c4 + v.y * M.c13 - v.z * M.c12 - v.w * M.c1;
    ctype const c6 = -v.x * M.c13 + v.y * M.c4 + v.z * M.c11 - v.w * M.c2;
    ctype const c7 = v.x * M.c12 - v.y * M.c11 + v.z * M.c4 - v.w * M.c3;
    ctype const c8 = -v.x * M.c14 + v.y * M.c3 - v.z * M.c2 + v.w * M.c11;
    ctype const c9 = -v.x * M.c3 - v.y * M.c14 + v.z * M.c1 + v.w * M.c12;
    ctype const c10 = v.x * M.c2 - v.y * M.c1 - v.z * M.c14 + v.w * M.c13;
    ctype const c11 = -v.x * M.c15 + v.y * M.c7 - v.z * M.c6 + v.w * M.c8;
    ctype const c12 = -v.x * M.c7 - v.y * M.c15 + v.z * M.c5 + v.w * M.c9;
    ctype const c13 = v.x * M.c6 - v.y * M.c5 - v.z * M.c15 + v.w * M.c10;
    ctype const c14 = v.x * M.c8 + v.y * M.c9 + v.z * M.c10 - v.w * M.c15;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 - v.w * M.c14;
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(MVec4ds<T> const& M, Scalar4ds<U> s)
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
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator*(Scalar4ds<T> s, MVec4ds<U> const& M)
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
    return MVec4ds<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// sta4ds gpr :: gpr(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& A,
                                                        MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
                     A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 + A.c2 * B.c6 - A.c3 * B.c5 - A.c4 * B.c7 +
                     A.c5 * B.c3 - A.c6 * B.c2 - A.c7 * B.c4;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c6 + A.c2 * B.c0 + A.c3 * B.c4 - A.c4 * B.c3 -
                     A.c5 * B.c7 + A.c6 * B.c1 - A.c7 * B.c5;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c0 + A.c4 * B.c2 -
                     A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c7 - A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c0 +
                     A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c1;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c3 + A.c2 * B.c7 - A.c3 * B.c1 - A.c4 * B.c6 +
                     A.c5 * B.c0 + A.c6 * B.c4 + A.c7 * B.c2;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c7 + A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
                     A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& A,
                                                        MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c3 + A.c2 * B.c6 - A.c3 * B.c5 - A.c4 * B.c7 +
                     A.c5 * B.c2 - A.c6 * B.c1 - A.c7 * B.c4;
    ctype const c1 = A.c0 * B.c1 - A.c1 * B.c6 + A.c2 * B.c3 + A.c3 * B.c4 - A.c4 * B.c2 -
                     A.c5 * B.c7 + A.c6 * B.c0 - A.c7 * B.c5;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c3 + A.c4 * B.c1 -
                     A.c5 * B.c0 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c4 * B.c4 -
                     A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c7 - A.c2 * B.c2 + A.c3 * B.c1 + A.c4 * B.c3 +
                     A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c2 + A.c2 * B.c7 - A.c3 * B.c0 - A.c4 * B.c6 +
                     A.c5 * B.c3 + A.c6 * B.c4 + A.c7 * B.c1;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c1 + A.c2 * B.c0 + A.c3 * B.c7 + A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c3 + A.c7 * B.c2;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c0 +
                     A.c5 * B.c1 + A.c6 * B.c2 + A.c7 * B.c3;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& A,
                                                        MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c1 + A.c4 * B.c7 -
                     A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c4;
    ctype const c1 = -A.c0 * B.c6 + A.c1 * B.c0 + A.c2 * B.c4 - A.c3 * B.c2 +
                     A.c4 * B.c3 + A.c5 * B.c7 - A.c6 * B.c1 - A.c7 * B.c5;
    ctype const c2 = A.c0 * B.c5 - A.c1 * B.c4 + A.c2 * B.c0 - A.c3 * B.c3 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c3 = -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0 -
                     A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 + A.c7 * B.c7;
    ctype const c4 = -A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c4 +
                     A.c4 * B.c0 + A.c5 * B.c6 - A.c6 * B.c5 - A.c7 * B.c1;
    ctype const c5 = -A.c0 * B.c3 - A.c1 * B.c7 + A.c2 * B.c1 + A.c3 * B.c5 -
                     A.c4 * B.c6 + A.c5 * B.c0 + A.c6 * B.c4 - A.c7 * B.c2;
    ctype const c6 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c7 + A.c3 * B.c6 + A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c0 - A.c7 * B.c3;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 - A.c3 * B.c7 - A.c4 * B.c1 -
                     A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& A,
                                                        PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c7 * ctype(ps);
    ctype const c1 = -A.c4 * ctype(ps);
    ctype const c2 = -A.c5 * ctype(ps);
    ctype const c3 = -A.c6 * ctype(ps);
    ctype const c4 = A.c1 * ctype(ps);
    ctype const c5 = A.c2 * ctype(ps);
    ctype const c6 = A.c3 * ctype(ps);
    ctype const c7 = A.c0 * ctype(ps);
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps,
                                                        MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c7;
    ctype const c1 = -ctype(ps) * B.c4;
    ctype const c2 = -ctype(ps) * B.c5;
    ctype const c3 = -ctype(ps) * B.c6;
    ctype const c4 = ctype(ps) * B.c1;
    ctype const c5 = ctype(ps) * B.c2;
    ctype const c6 = ctype(ps) * B.c3;
    ctype const c7 = ctype(ps) * B.c0;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& M,
                                                        TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * t.z - M.c3 * t.y - M.c4 * t.w - M.c7 * t.x;
    ctype const c1 = -M.c1 * t.z + M.c3 * t.x - M.c5 * t.w - M.c7 * t.y;
    ctype const c2 = M.c1 * t.y - M.c2 * t.x - M.c6 * t.w - M.c7 * t.z;
    ctype const c3 = -M.c4 * t.x - M.c5 * t.y - M.c6 * t.z - M.c7 * t.w;
    ctype const c4 = M.c0 * t.x + M.c1 * t.w + M.c5 * t.z - M.c6 * t.y;
    ctype const c5 = M.c0 * t.y + M.c2 * t.w - M.c4 * t.z + M.c6 * t.x;
    ctype const c6 = M.c0 * t.z + M.c3 * t.w + M.c4 * t.y - M.c5 * t.x;
    ctype const c7 = M.c0 * t.w + M.c1 * t.x + M.c2 * t.y + M.c3 * t.z;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t,
                                                        MVec4ds_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * M.c7 - t.y * M.c3 + t.z * M.c2 - t.w * M.c4;
    ctype const c1 = t.x * M.c3 + t.y * M.c7 - t.z * M.c1 - t.w * M.c5;
    ctype const c2 = -t.x * M.c2 + t.y * M.c1 + t.z * M.c7 - t.w * M.c6;
    ctype const c3 = -t.x * M.c4 - t.y * M.c5 - t.z * M.c6 + t.w * M.c7;
    ctype const c4 = t.x * M.c0 + t.y * M.c6 - t.z * M.c5 - t.w * M.c1;
    ctype const c5 = -t.x * M.c6 + t.y * M.c0 + t.z * M.c4 - t.w * M.c2;
    ctype const c6 = t.x * M.c5 - t.y * M.c4 + t.z * M.c0 - t.w * M.c3;
    ctype const c7 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3 + t.w * M.c0;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& M,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        M.c1 * B.vx + M.c2 * B.vy + M.c3 * B.vz - M.c4 * B.mx - M.c5 * B.my - M.c6 * B.mz;
    ctype const c1 =
        M.c0 * B.vx + M.c2 * B.mz - M.c3 * B.my + M.c5 * B.vz - M.c6 * B.vy - M.c7 * B.mx;
    ctype const c2 =
        M.c0 * B.vy - M.c1 * B.mz + M.c3 * B.mx - M.c4 * B.vz + M.c6 * B.vx - M.c7 * B.my;
    ctype const c3 =
        M.c0 * B.vz + M.c1 * B.my - M.c2 * B.mx + M.c4 * B.vy - M.c5 * B.vx - M.c7 * B.mz;
    ctype const c4 =
        M.c0 * B.mx - M.c2 * B.vz + M.c3 * B.vy + M.c5 * B.mz - M.c6 * B.my + M.c7 * B.vx;
    ctype const c5 =
        M.c0 * B.my + M.c1 * B.vz - M.c3 * B.vx - M.c4 * B.mz + M.c6 * B.mx + M.c7 * B.vy;
    ctype const c6 =
        M.c0 * B.mz - M.c1 * B.vy + M.c2 * B.vx + M.c4 * B.my - M.c5 * B.mx + M.c7 * B.vz;
    ctype const c7 =
        M.c1 * B.mx + M.c2 * B.my + M.c3 * B.mz + M.c4 * B.vx + M.c5 * B.vy + M.c6 * B.vz;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B,
                                                        MVec4ds_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c1 + B.vy * M.c2 + B.vz * M.c3 - B.mx * M.c4 - B.my * M.c5 - B.mz * M.c6;
    ctype const c1 =
        B.vx * M.c0 + B.vy * M.c6 - B.vz * M.c5 - B.mx * M.c7 + B.my * M.c3 - B.mz * M.c2;
    ctype const c2 = -B.vx * M.c6 + B.vy * M.c0 + B.vz * M.c4 - B.mx * M.c3 -
                     B.my * M.c7 + B.mz * M.c1;
    ctype const c3 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c0 + B.mx * M.c2 - B.my * M.c1 - B.mz * M.c7;
    ctype const c4 =
        B.vx * M.c7 - B.vy * M.c3 + B.vz * M.c2 + B.mx * M.c0 + B.my * M.c6 - B.mz * M.c5;
    ctype const c5 =
        B.vx * M.c3 + B.vy * M.c7 - B.vz * M.c1 - B.mx * M.c6 + B.my * M.c0 + B.mz * M.c4;
    ctype const c6 = -B.vx * M.c2 + B.vy * M.c1 + B.vz * M.c7 + B.mx * M.c5 -
                     B.my * M.c4 + B.mz * M.c0;
    ctype const c7 =
        B.vx * M.c4 + B.vy * M.c5 + B.vz * M.c6 + B.mx * M.c1 + B.my * M.c2 + B.mz * M.c3;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& A,
                                                        Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * v.x + A.c1 * v.w + A.c5 * v.z - A.c6 * v.y;
    ctype const c1 = A.c0 * v.y + A.c2 * v.w - A.c4 * v.z + A.c6 * v.x;
    ctype const c2 = A.c0 * v.z + A.c3 * v.w + A.c4 * v.y - A.c5 * v.x;
    ctype const c3 = A.c0 * v.w + A.c1 * v.x + A.c2 * v.y + A.c3 * v.z;
    ctype const c4 = -A.c2 * v.z + A.c3 * v.y + A.c4 * v.w + A.c7 * v.x;
    ctype const c5 = A.c1 * v.z - A.c3 * v.x + A.c5 * v.w + A.c7 * v.y;
    ctype const c6 = -A.c1 * v.y + A.c2 * v.x + A.c6 * v.w + A.c7 * v.z;
    ctype const c7 = A.c4 * v.x + A.c5 * v.y + A.c6 * v.z + A.c7 * v.w;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v,
                                                        MVec4ds_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * B.c0 + v.y * B.c6 - v.z * B.c5 - v.w * B.c1;
    ctype const c1 = -v.x * B.c6 + v.y * B.c0 + v.z * B.c4 - v.w * B.c2;
    ctype const c2 = v.x * B.c5 - v.y * B.c4 + v.z * B.c0 - v.w * B.c3;
    ctype const c3 = -v.x * B.c1 - v.y * B.c2 - v.z * B.c3 + v.w * B.c0;
    ctype const c4 = -v.x * B.c7 + v.y * B.c3 - v.z * B.c2 + v.w * B.c4;
    ctype const c5 = -v.x * B.c3 - v.y * B.c7 + v.z * B.c1 + v.w * B.c5;
    ctype const c6 = v.x * B.c2 - v.y * B.c1 - v.z * B.c7 + v.w * B.c6;
    ctype const c7 = v.x * B.c4 + v.y * B.c5 + v.z * B.c6 - v.w * B.c7;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(MVec4ds_E<T> const& M,
                                                        Scalar4ds<U> s)
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
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(Scalar4ds<T> s,
                                                        MVec4ds_E<U> const& M)
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
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& A,
                                                        MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c3 -
                     A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 + A.c7 * B.c7;
    ctype const c1 = A.c0 * B.c3 + A.c1 * B.c6 - A.c2 * B.c5 - A.c3 * B.c0 + A.c4 * B.c7 -
                     A.c5 * B.c2 + A.c6 * B.c1 - A.c7 * B.c4;
    ctype const c2 = -A.c0 * B.c6 + A.c1 * B.c3 + A.c2 * B.c4 - A.c3 * B.c1 +
                     A.c4 * B.c2 + A.c5 * B.c7 - A.c6 * B.c0 - A.c7 * B.c5;
    ctype const c3 = A.c0 * B.c5 - A.c1 * B.c4 + A.c2 * B.c3 - A.c3 * B.c2 - A.c4 * B.c1 +
                     A.c5 * B.c0 + A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c4 = -A.c0 * B.c7 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c4 +
                     A.c4 * B.c3 + A.c5 * B.c6 - A.c6 * B.c5 - A.c7 * B.c0;
    ctype const c5 = -A.c0 * B.c2 - A.c1 * B.c7 + A.c2 * B.c0 + A.c3 * B.c5 -
                     A.c4 * B.c6 + A.c5 * B.c3 + A.c6 * B.c4 - A.c7 * B.c1;
    ctype const c6 = A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c7 + A.c3 * B.c6 + A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c3 - A.c7 * B.c2;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 - A.c3 * B.c7 - A.c4 * B.c0 -
                     A.c5 * B.c1 - A.c6 * B.c2 + A.c7 * B.c3;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& A,
                                                        PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c4 * ctype(ps);
    ctype const c1 = A.c5 * ctype(ps);
    ctype const c2 = A.c6 * ctype(ps);
    ctype const c3 = A.c7 * ctype(ps);
    ctype const c4 = -A.c0 * ctype(ps);
    ctype const c5 = -A.c1 * ctype(ps);
    ctype const c6 = -A.c2 * ctype(ps);
    ctype const c7 = -A.c3 * ctype(ps);
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps,
                                                        MVec4ds_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c4;
    ctype const c1 = -ctype(ps) * B.c5;
    ctype const c2 = -ctype(ps) * B.c6;
    ctype const c3 = -ctype(ps) * B.c7;
    ctype const c4 = ctype(ps) * B.c0;
    ctype const c5 = ctype(ps) * B.c1;
    ctype const c6 = ctype(ps) * B.c2;
    ctype const c7 = ctype(ps) * B.c3;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& M,
                                                        TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * t.x - M.c5 * t.y - M.c6 * t.z + M.c7 * t.w;
    ctype const c1 = M.c1 * t.z - M.c2 * t.y + M.c4 * t.w - M.c7 * t.x;
    ctype const c2 = -M.c0 * t.z + M.c2 * t.x + M.c5 * t.w - M.c7 * t.y;
    ctype const c3 = M.c0 * t.y - M.c1 * t.x + M.c6 * t.w - M.c7 * t.z;
    ctype const c4 = -M.c0 * t.w + M.c3 * t.x + M.c5 * t.z - M.c6 * t.y;
    ctype const c5 = -M.c1 * t.w + M.c3 * t.y - M.c4 * t.z + M.c6 * t.x;
    ctype const c6 = -M.c2 * t.w + M.c3 * t.z + M.c4 * t.y - M.c5 * t.x;
    ctype const c7 = M.c0 * t.x + M.c1 * t.y + M.c2 * t.z - M.c3 * t.w;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t,
                                                        MVec4ds_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c4 - t.y * M.c5 - t.z * M.c6 + t.w * M.c7;
    ctype const c1 = t.x * M.c7 - t.y * M.c2 + t.z * M.c1 - t.w * M.c4;
    ctype const c2 = t.x * M.c2 + t.y * M.c7 - t.z * M.c0 - t.w * M.c5;
    ctype const c3 = -t.x * M.c1 + t.y * M.c0 + t.z * M.c7 - t.w * M.c6;
    ctype const c4 = t.x * M.c3 + t.y * M.c6 - t.z * M.c5 - t.w * M.c0;
    ctype const c5 = -t.x * M.c6 + t.y * M.c3 + t.z * M.c4 - t.w * M.c1;
    ctype const c6 = t.x * M.c5 - t.y * M.c4 + t.z * M.c3 - t.w * M.c2;
    ctype const c7 = -t.x * M.c0 - t.y * M.c1 - t.z * M.c2 + t.w * M.c3;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& M,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        M.c1 * B.mz - M.c2 * B.my - M.c3 * B.vx - M.c5 * B.vz + M.c6 * B.vy - M.c7 * B.mx;
    ctype const c1 = -M.c0 * B.mz + M.c2 * B.mx - M.c3 * B.vy + M.c4 * B.vz -
                     M.c6 * B.vx - M.c7 * B.my;
    ctype const c2 =
        M.c0 * B.my - M.c1 * B.mx - M.c3 * B.vz - M.c4 * B.vy + M.c5 * B.vx - M.c7 * B.mz;
    ctype const c3 = -M.c0 * B.vx - M.c1 * B.vy - M.c2 * B.vz - M.c4 * B.mx -
                     M.c5 * B.my - M.c6 * B.mz;
    ctype const c4 =
        M.c1 * B.vz - M.c2 * B.vy + M.c3 * B.mx + M.c5 * B.mz - M.c6 * B.my - M.c7 * B.vx;
    ctype const c5 = -M.c0 * B.vz + M.c2 * B.vx + M.c3 * B.my - M.c4 * B.mz +
                     M.c6 * B.mx - M.c7 * B.vy;
    ctype const c6 =
        M.c0 * B.vy - M.c1 * B.vx + M.c3 * B.mz + M.c4 * B.my - M.c5 * B.mx - M.c7 * B.vz;
    ctype const c7 =
        M.c0 * B.mx + M.c1 * B.my + M.c2 * B.mz - M.c4 * B.vx - M.c5 * B.vy - M.c6 * B.vz;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B,
                                                        MVec4ds_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c3 + B.vy * M.c6 - B.vz * M.c5 - B.mx * M.c7 + B.my * M.c2 - B.mz * M.c1;
    ctype const c1 = -B.vx * M.c6 + B.vy * M.c3 + B.vz * M.c4 - B.mx * M.c2 -
                     B.my * M.c7 + B.mz * M.c0;
    ctype const c2 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c3 + B.mx * M.c1 - B.my * M.c0 - B.mz * M.c7;
    ctype const c3 =
        B.vx * M.c0 + B.vy * M.c1 + B.vz * M.c2 - B.mx * M.c4 - B.my * M.c5 - B.mz * M.c6;
    ctype const c4 =
        B.vx * M.c7 - B.vy * M.c2 + B.vz * M.c1 + B.mx * M.c3 + B.my * M.c6 - B.mz * M.c5;
    ctype const c5 =
        B.vx * M.c2 + B.vy * M.c7 - B.vz * M.c0 - B.mx * M.c6 + B.my * M.c3 + B.mz * M.c4;
    ctype const c6 = -B.vx * M.c1 + B.vy * M.c0 + B.vz * M.c7 + B.mx * M.c5 -
                     B.my * M.c4 + B.mz * M.c3;
    ctype const c7 =
        B.vx * M.c4 + B.vy * M.c5 + B.vz * M.c6 + B.mx * M.c0 + B.my * M.c1 + B.mz * M.c2;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& M,
                                                        Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c0 * v.x - M.c1 * v.y - M.c2 * v.z + M.c3 * v.w;
    ctype const c1 = M.c0 * v.w - M.c3 * v.x - M.c5 * v.z + M.c6 * v.y;
    ctype const c2 = M.c1 * v.w - M.c3 * v.y + M.c4 * v.z - M.c6 * v.x;
    ctype const c3 = M.c2 * v.w - M.c3 * v.z - M.c4 * v.y + M.c5 * v.x;
    ctype const c4 = M.c1 * v.z - M.c2 * v.y + M.c4 * v.w - M.c7 * v.x;
    ctype const c5 = -M.c0 * v.z + M.c2 * v.x + M.c5 * v.w - M.c7 * v.y;
    ctype const c6 = M.c0 * v.y - M.c1 * v.x + M.c6 * v.w - M.c7 * v.z;
    ctype const c7 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z + M.c7 * v.w;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v,
                                                        MVec4ds_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * M.c0 - v.y * M.c1 - v.z * M.c2 + v.w * M.c3;
    ctype const c1 = v.x * M.c3 + v.y * M.c6 - v.z * M.c5 - v.w * M.c0;
    ctype const c2 = -v.x * M.c6 + v.y * M.c3 + v.z * M.c4 - v.w * M.c1;
    ctype const c3 = v.x * M.c5 - v.y * M.c4 + v.z * M.c3 - v.w * M.c2;
    ctype const c4 = -v.x * M.c7 + v.y * M.c2 - v.z * M.c1 + v.w * M.c4;
    ctype const c5 = -v.x * M.c2 - v.y * M.c7 + v.z * M.c0 + v.w * M.c5;
    ctype const c6 = v.x * M.c1 - v.y * M.c0 - v.z * M.c7 + v.w * M.c6;
    ctype const c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6 - v.w * M.c7;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(MVec4ds_U<T> const& M,
                                                        Scalar4ds<U> s)
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
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(Scalar4ds<T> s,
                                                        MVec4ds_U<U> const& M)
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
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps1,
                                                        PScalar4ds<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(-ctype(ps1) * ctype(ps2));
}

// sta4ds gpr :: gpr(ps,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps,
                                                     TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * t.x;
    ctype const c1 = -ctype(ps) * t.y;
    ctype const c2 = -ctype(ps) * t.z;
    ctype const c3 = -ctype(ps) * t.w;
    return Vec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(trivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t,
                                                     PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(ps);
    ctype const c1 = t.y * ctype(ps);
    ctype const c2 = t.z * ctype(ps);
    ctype const c3 = t.w * ctype(ps);
    return Vec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps,
                                                       BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.mx;
    ctype const c1 = -ctype(ps) * B.my;
    ctype const c2 = -ctype(ps) * B.mz;
    ctype const c3 = ctype(ps) * B.vx;
    ctype const c4 = ctype(ps) * B.vy;
    ctype const c5 = ctype(ps) * B.vz;
    return BiVec4ds<ctype>(c0, c1, c2, c3, c4, c5);
}

// sta4ds gpr :: gpr(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B,
                                                       PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.mx * ctype(ps);
    ctype const c1 = -B.my * ctype(ps);
    ctype const c2 = -B.mz * ctype(ps);
    ctype const c3 = B.vx * ctype(ps);
    ctype const c4 = B.vy * ctype(ps);
    ctype const c5 = B.vz * ctype(ps);
    return BiVec4ds<ctype>(c0, c1, c2, c3, c4, c5);
}

// sta4ds gpr :: gpr(ps,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps,
                                                        Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = ctype(ps) * v.w;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(vec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v,
                                                        PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * ctype(ps);
    ctype const c1 = -v.y * ctype(ps);
    ctype const c2 = -v.z * ctype(ps);
    ctype const c3 = -v.w * ctype(ps);
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> operator*(PScalar4ds<T> ps, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(ctype(ps) * ctype(s));
}

// sta4ds gpr :: gpr(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar4ds<std::common_type_t<T, U>> operator*(Scalar4ds<T> s, PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar4ds<ctype>(ctype(s) * ctype(ps));
}

// sta4ds gpr :: gpr(trivec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t1,
                                                        TriVec4ds<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t1.x * t2.x - t1.y * t2.y - t1.z * t2.z + t1.w * t2.w;
    ctype const c1 = t1.x * t2.w - t1.w * t2.x;
    ctype const c2 = t1.y * t2.w - t1.w * t2.y;
    ctype const c3 = t1.z * t2.w - t1.w * t2.z;
    ctype const c4 = t1.y * t2.z - t1.z * t2.y;
    ctype const c5 = -t1.x * t2.z + t1.z * t2.x;
    ctype const c6 = t1.x * t2.y - t1.y * t2.x;
    ctype const c7 = 0.0;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(trivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.vz + t.z * B.vy - t.w * B.mx;
    ctype const c1 = t.x * B.vz - t.z * B.vx - t.w * B.my;
    ctype const c2 = -t.x * B.vy + t.y * B.vx - t.w * B.mz;
    ctype const c3 = -t.x * B.mx - t.y * B.my - t.z * B.mz;
    ctype const c4 = t.y * B.mz - t.z * B.my - t.w * B.vx;
    ctype const c5 = -t.x * B.mz + t.z * B.mx - t.w * B.vy;
    ctype const c6 = t.x * B.my - t.y * B.mx - t.w * B.vz;
    ctype const c7 = -t.x * B.vx - t.y * B.vy - t.z * B.vz;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(bivec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B,
                                                        TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * t.z - B.vz * t.y - B.mx * t.w;
    ctype const c1 = -B.vx * t.z + B.vz * t.x - B.my * t.w;
    ctype const c2 = B.vx * t.y - B.vy * t.x - B.mz * t.w;
    ctype const c3 = -B.mx * t.x - B.my * t.y - B.mz * t.z;
    ctype const c4 = B.vx * t.w + B.my * t.z - B.mz * t.y;
    ctype const c5 = B.vy * t.w - B.mx * t.z + B.mz * t.x;
    ctype const c6 = B.vz * t.w + B.mx * t.y - B.my * t.x;
    ctype const c7 = B.vx * t.x + B.vy * t.y + B.vz * t.z;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(trivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t,
                                                        Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -t.y * v.z + t.z * v.y;
    ctype const c2 = t.x * v.z - t.z * v.x;
    ctype const c3 = -t.x * v.y + t.y * v.x;
    ctype const c4 = t.x * v.w - t.w * v.x;
    ctype const c5 = t.y * v.w - t.w * v.y;
    ctype const c6 = t.z * v.w - t.w * v.z;
    ctype const c7 = -t.x * v.x - t.y * v.y - t.z * v.z + t.w * v.w;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(vec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v,
                                                        TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.y * t.z - v.z * t.y;
    ctype const c2 = -v.x * t.z + v.z * t.x;
    ctype const c3 = v.x * t.y - v.y * t.x;
    ctype const c4 = -v.x * t.w + v.w * t.x;
    ctype const c5 = -v.y * t.w + v.w * t.y;
    ctype const c6 = -v.z * t.w + v.w * t.z;
    ctype const c7 = v.x * t.x + v.y * t.y + v.z * t.z - v.w * t.w;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> operator*(TriVec4ds<T> const& t,
                                                        Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> operator*(Scalar4ds<T> s,
                                                        TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(bivec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B1,
                                                        BiVec4ds<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vx * B2.vx + B1.vy * B2.vy + B1.vz * B2.vz - B1.mx * B2.mx -
                     B1.my * B2.my - B1.mz * B2.mz;
    ctype const c1 = B1.vy * B2.mz - B1.vz * B2.my + B1.my * B2.vz - B1.mz * B2.vy;
    ctype const c2 = -B1.vx * B2.mz + B1.vz * B2.mx - B1.mx * B2.vz + B1.mz * B2.vx;
    ctype const c3 = B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx;
    ctype const c4 = -B1.vy * B2.vz + B1.vz * B2.vy + B1.my * B2.mz - B1.mz * B2.my;
    ctype const c5 = B1.vx * B2.vz - B1.vz * B2.vx - B1.mx * B2.mz + B1.mz * B2.mx;
    ctype const c6 = -B1.vx * B2.vy + B1.vy * B2.vx + B1.mx * B2.my - B1.my * B2.mx;
    ctype const c7 = B1.vx * B2.mx + B1.vy * B2.my + B1.vz * B2.mz + B1.mx * B2.vx +
                     B1.my * B2.vy + B1.mz * B2.vz;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(bivec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B,
                                                        Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * v.w + B.my * v.z - B.mz * v.y;
    ctype const c1 = B.vy * v.w - B.mx * v.z + B.mz * v.x;
    ctype const c2 = B.vz * v.w + B.mx * v.y - B.my * v.x;
    ctype const c3 = B.vx * v.x + B.vy * v.y + B.vz * v.z;
    ctype const c4 = -B.vy * v.z + B.vz * v.y + B.mx * v.w;
    ctype const c5 = B.vx * v.z - B.vz * v.x + B.my * v.w;
    ctype const c6 = -B.vx * v.y + B.vy * v.x + B.mz * v.w;
    ctype const c7 = B.mx * v.x + B.my * v.y + B.mz * v.z;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(vec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_U<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.mz - v.z * B.my - v.w * B.vx;
    ctype const c1 = -v.x * B.mz + v.z * B.mx - v.w * B.vy;
    ctype const c2 = v.x * B.my - v.y * B.mx - v.w * B.vz;
    ctype const c3 = -v.x * B.vx - v.y * B.vy - v.z * B.vz;
    ctype const c4 = v.y * B.vz - v.z * B.vy + v.w * B.mx;
    ctype const c5 = -v.x * B.vz + v.z * B.vx + v.w * B.my;
    ctype const c6 = v.x * B.vy - v.y * B.vx + v.w * B.mz;
    ctype const c7 = v.x * B.mx + v.y * B.my + v.z * B.mz;
    return MVec4ds_U<ctype>(Vec4ds<ctype>(c0, c1, c2, c3),
                            TriVec4ds<ctype>(c4, c5, c6, c7));
}

// sta4ds gpr :: gpr(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> operator*(BiVec4ds<T> const& B,
                                                       Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec4ds<ctype>(c0, c1, c2, c3, c4, c5);
}

// sta4ds gpr :: gpr(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> operator*(Scalar4ds<T> s,
                                                       BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec4ds<ctype>(c0, c1, c2, c3, c4, c5);
}

// sta4ds gpr :: gpr(vec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v1,
                                                        Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.x - v1.y * v2.y - v1.z * v2.z + v1.w * v2.w;
    ctype const c1 = v1.x * v2.w - v1.w * v2.x;
    ctype const c2 = v1.y * v2.w - v1.w * v2.y;
    ctype const c3 = v1.z * v2.w - v1.w * v2.z;
    ctype const c4 = v1.y * v2.z - v1.z * v2.y;
    ctype const c5 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c6 = v1.x * v2.y - v1.y * v2.x;
    ctype const c7 = 0.0;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(c0), BiVec4ds<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar4ds<ctype>(c7));
}

// sta4ds gpr :: gpr(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> operator*(Vec4ds<T> const& v, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> operator*(Scalar4ds<T> s, Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec4ds<ctype>(c0, c1, c2, c3);
}

// sta4ds gpr :: gpr(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> operator*(Scalar4ds<T> s1, Scalar4ds<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar4ds<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors
// w.r.t. the geometric product:
// for k-blades: A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ... (from reversion)
////////////////////////////////////////////////////////////////////////////////
// HINT: inv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Scalar4ds<T> inv(Scalar4ds<T> s)
{
    T sq_n = bulk_nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar4ds<T>(rev(s) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec4ds<T> inv(Vec4ds<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v) = v/bulk_sq_nrm(v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    T sq_n = bulk_nrm_sq(v);
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return Vec4ds<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> inv(BiVec4ds<T> const& B)
{
    auto bc = B * conj(B);
    auto bcmap = gr0(bc) + gr2(bc) - gr4(bc);
    // fmt::println("B={}", B);
    // fmt::println("bc={}", bc);
    // fmt::println("bcmap={}", bcmap);
    // fmt::println("bc*bcmap={}", bc * bcmap);
    T sq_n = T(gr0(bc * bcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "bivector");
    return gr2(conj(B) * bcmap) / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline TriVec4ds<T> inv(TriVec4ds<T> const& t)
{
    auto tc = t * conj(t);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("t={}", t);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "trivector");
    return gr3(conj(t) * tcmap) / sq_n;
}

// TODO: invert the pseudoscalar

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> inv(MVec4ds_E<T> const& E)
{
    auto tc = E * conj(E);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("E={}", E);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "even-grade multivector");
    return conj(E) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_U<T> inv(MVec4ds_U<T> const& U)
{
    auto tc = U * conj(U);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("E={}", E);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "odd-grade multivector");
    return conj(U) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016 left and a right inverse are the same (see paper
// of Hitzer, Sangwine)
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds<T> inv(MVec4ds<T> const& M)
{
    auto tc = M * conj(M);
    auto tcmap = gr0(tc) + gr1(tc) + gr2(tc) - gr3(tc) - gr4(tc);
    // fmt::println("M={}", M);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "multivector");
    return conj(M) * tcmap / sq_n;
}

} // namespace hd::ga::sta