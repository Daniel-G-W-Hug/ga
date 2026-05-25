#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_pga3dp_ops_basics.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga3dp product operations:
//
// - dot()                   -> dot product
// - rdot()                  -> regressive dot product
// - wdg(), join()           -> wedge product (join as convenience interface)
// - rwdg(), meet()          -> regressive wedge product (meet as convenience interface)
// - twdg1()                 -> transwedge product (k=1)
// - rtwdg1()                -> regressive transwedge product (k=1)
// - operator<<()            -> left bulk contraction
// - operator>>()            -> right bulk contraction
// - cmt()                   -> commutator product (= asymmetric part of gpr)
// - rcmt()                  -> regressive commutator product (=asymmetric part of rgpr)
// - operator*()             -> geometric product (= gpr)
// - rgpr()                  -> regressive geometric product
// - inv()                   -> inversion operation (w.r.t. geometric product)
//
// - l_bulk_contract3dp()    -> left bulk contraction
// - l_weight_contract3dp()  -> left weight contraction
// - r_bulk_contract3dp()    -> right bulk contraction
// - r_weight_contract3dp()  -> right weight contraction
//
// Bulk expansions are the regressive versions of the corresponding weight contractions.
// Weight expansions are the regressive versions of the corresponding bulk contractions.
//
// - l_bulk_expand3dp()      -> left bulk expansion
// - l_weight_expand3dp()    -> left weight expansion
// - r_bulk_expand3dp()      -> right bulk expansion
// - r_weight_expand3dp()    -> right weight expansion
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot products for 3dp (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = +1, e4^2 = 0
//
// Via the exomorphism G(a ^ b) = G(a) ^ G(b) the extended metric G is defined as a 2^n by
// 2^n matrix via the outer product, which in turn allows for definition of an extended
// dot product for all grades of a full multi-vector from scalar to pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

// pga3dp dot :: dot(mv,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(MVec3dp<T> const& A,
                                                  MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 +
                            A.c8 * B.c8 + A.c9 * B.c9 + A.c10 * B.c10 + A.c14 * B.c14);
}

// pga3dp dot :: dot(mv_e,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(MVec3dp_E<T> const& A,
                                                  MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(A.c0 * B.c0 + A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6);
}

// pga3dp dot :: dot(mv_u,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(MVec3dp_U<T> const& A,
                                                  MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c7 * B.c7);
}

// pga3dp dot :: dot(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp dot :: dot(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(t1.w * t2.w);
}

// pga3dp dot :: dot(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(BiVec3dp<T> const& B1,
                                                  BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(B1.mx * B2.mx + B1.my * B2.my + B1.mz * B2.mz);
}

// pga3dp dot :: dot(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(Vec3dp<T> const& v1,
                                                  Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// pga3dp dot :: dot(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive dot products for 3dp (= defined for equal grades exclusively)
//
// rdot(v1,v2) = l_cmpl( dot(r_cmpl(v1),r_cmpl(v2)) )
//
// returns a pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

// pga3dp rdot :: rdot(mv,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(MVec3dp<T> const& A,
                                                    MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7 +
                             A.c11 * B.c11 + A.c12 * B.c12 + A.c13 * B.c13 +
                             A.c15 * B.c15);
}

// pga3dp rdot :: rdot(mv_e,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(MVec3dp_E<T> const& A,
                                                    MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c7 * B.c7);
}

// pga3dp rdot :: rdot(mv_u,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(MVec3dp_U<T> const& A,
                                                    MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(A.c3 * B.c3 + A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6);
}

// pga3dp rdot :: rdot(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(PScalar3dp<T> ps1, PScalar3dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps1) * ctype(ps2));
}

// pga3dp rdot :: rdot(trivec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(TriVec3dp<T> const& t1,
                                                    TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(t1.x * t2.x + t1.y * t2.y + t1.z * t2.z);
}

// pga3dp rdot :: rdot(bivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(BiVec3dp<T> const& B1,
                                                    BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(B1.vx * B2.vx + B1.vy * B2.vy + B1.vz * B2.vz);
}

// pga3dp rdot :: rdot(vec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(Vec3dp<T> const& v1,
                                                    Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(v1.w * v2.w);
}

// pga3dp rdot :: rdot(s,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot([[maybe_unused]] Scalar3dp<T>,
                                                    [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////

// pga3dp wdg :: wdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c4 + A.c4 * B.c1 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 - A.c2 * B.c4 + A.c4 * B.c2 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 - A.c3 * B.c4 + A.c4 * B.c3 + A.c7 * B.c0;
    ctype const c8 = A.c0 * B.c8 + A.c2 * B.c3 - A.c3 * B.c2 + A.c8 * B.c0;
    ctype const c9 = A.c0 * B.c9 - A.c1 * B.c3 + A.c3 * B.c1 + A.c9 * B.c0;
    ctype const c10 = A.c0 * B.c10 + A.c1 * B.c2 - A.c2 * B.c1 + A.c10 * B.c0;
    ctype const c11 = A.c0 * B.c11 - A.c2 * B.c7 + A.c3 * B.c6 + A.c4 * B.c8 +
                      A.c6 * B.c3 - A.c7 * B.c2 + A.c8 * B.c4 + A.c11 * B.c0;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c7 - A.c3 * B.c5 + A.c4 * B.c9 -
                      A.c5 * B.c3 + A.c7 * B.c1 + A.c9 * B.c4 + A.c12 * B.c0;
    ctype const c13 = A.c0 * B.c13 - A.c1 * B.c6 + A.c2 * B.c5 + A.c4 * B.c10 +
                      A.c5 * B.c2 - A.c6 * B.c1 + A.c10 * B.c4 + A.c13 * B.c0;
    ctype const c14 = A.c0 * B.c14 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                      A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 +
                      A.c4 * B.c14 - A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 -
                      A.c8 * B.c5 - A.c9 * B.c6 - A.c10 * B.c7 - A.c11 * B.c1 -
                      A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& A,
                                                MVec3dp_E<U> const& B)
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
    ctype const c11 = -A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c4 + A.c11 * B.c0;
    ctype const c12 = A.c1 * B.c3 - A.c3 * B.c1 + A.c4 * B.c5 + A.c12 * B.c0;
    ctype const c13 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c4 * B.c6 + A.c13 * B.c0;
    ctype const c14 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c7 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 -
                      A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 + A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& A,
                                                MVec3dp<U> const& B)
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
    ctype const c11 = A.c0 * B.c11 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c4;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c3 + A.c3 * B.c1 + A.c5 * B.c4;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c4;
    ctype const c14 = A.c0 * B.c14 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                      A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 + A.c7 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& A,
                                                MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c0 * B.c1;
    ctype const c3 = A.c0 * B.c2;
    ctype const c4 = A.c0 * B.c3;
    ctype const c5 = -A.c1 * B.c3 + A.c4 * B.c0;
    ctype const c6 = -A.c2 * B.c3 + A.c4 * B.c1;
    ctype const c7 = -A.c3 * B.c3 + A.c4 * B.c2;
    ctype const c8 = A.c2 * B.c2 - A.c3 * B.c1;
    ctype const c9 = -A.c1 * B.c2 + A.c3 * B.c0;
    ctype const c10 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype const c11 = A.c0 * B.c4 + A.c6 * B.c2 - A.c7 * B.c1 + A.c8 * B.c3;
    ctype const c12 = A.c0 * B.c5 - A.c5 * B.c2 + A.c7 * B.c0 + A.c9 * B.c3;
    ctype const c13 = A.c0 * B.c6 + A.c5 * B.c1 - A.c6 * B.c0 + A.c10 * B.c3;
    ctype const c14 = A.c0 * B.c7 - A.c8 * B.c0 - A.c9 * B.c1 - A.c10 * B.c2;
    ctype const c15 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c7 -
                      A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c2 - A.c14 * B.c3;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp_U<T> const& A,
                                                MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c1 * B.c0;
    ctype const c3 = A.c2 * B.c0;
    ctype const c4 = A.c3 * B.c0;
    ctype const c5 = -A.c0 * B.c4 + A.c3 * B.c1;
    ctype const c6 = -A.c1 * B.c4 + A.c3 * B.c2;
    ctype const c7 = -A.c2 * B.c4 + A.c3 * B.c3;
    ctype const c8 = A.c1 * B.c3 - A.c2 * B.c2;
    ctype const c9 = -A.c0 * B.c3 + A.c2 * B.c1;
    ctype const c10 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c11 = -A.c1 * B.c7 + A.c2 * B.c6 + A.c3 * B.c8 + A.c4 * B.c0;
    ctype const c12 = A.c0 * B.c7 - A.c2 * B.c5 + A.c3 * B.c9 + A.c5 * B.c0;
    ctype const c13 = -A.c0 * B.c6 + A.c1 * B.c5 + A.c3 * B.c10 + A.c6 * B.c0;
    ctype const c14 = -A.c0 * B.c8 - A.c1 * B.c9 - A.c2 * B.c10 + A.c7 * B.c0;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 -
                      A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 - A.c7 * B.c4;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& M, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(M.c0 * ctype(ps));
}

// pga3dp wdg :: wdg(ps,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(PScalar3dp<T> ps, MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * M.c0);
}

// pga3dp wdg :: wdg(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& M,
                                                TriVec3dp<U> const& t)
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
    ctype const c15 = M.c1 * t.x + M.c2 * t.y + M.c3 * t.z + M.c4 * t.w;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t,
                                                MVec3dp<U> const& M)
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
    ctype const c15 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3 - t.w * M.c4;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& M, BiVec3dp<U> const& B)
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
    ctype const c11 = -M.c2 * B.vz + M.c3 * B.vy + M.c4 * B.mx;
    ctype const c12 = M.c1 * B.vz - M.c3 * B.vx + M.c4 * B.my;
    ctype const c13 = -M.c1 * B.vy + M.c2 * B.vx + M.c4 * B.mz;
    ctype const c14 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz;
    ctype const c15 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                      M.c9 * B.vy - M.c10 * B.vz;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B, MVec3dp<U> const& M)
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
    ctype const c11 = B.vy * M.c3 - B.vz * M.c2 + B.mx * M.c4;
    ctype const c12 = -B.vx * M.c3 + B.vz * M.c1 + B.my * M.c4;
    ctype const c13 = B.vx * M.c2 - B.vy * M.c1 + B.mz * M.c4;
    ctype const c14 = -B.mx * M.c1 - B.my * M.c2 - B.mz * M.c3;
    ctype const c15 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                      B.my * M.c6 - B.mz * M.c7;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& M, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = M.c0 * v.w;
    ctype const c5 = -M.c1 * v.w + M.c4 * v.x;
    ctype const c6 = -M.c2 * v.w + M.c4 * v.y;
    ctype const c7 = -M.c3 * v.w + M.c4 * v.z;
    ctype const c8 = M.c2 * v.z - M.c3 * v.y;
    ctype const c9 = -M.c1 * v.z + M.c3 * v.x;
    ctype const c10 = M.c1 * v.y - M.c2 * v.x;
    ctype const c11 = M.c6 * v.z - M.c7 * v.y + M.c8 * v.w;
    ctype const c12 = -M.c5 * v.z + M.c7 * v.x + M.c9 * v.w;
    ctype const c13 = M.c5 * v.y - M.c6 * v.x + M.c10 * v.w;
    ctype const c14 = -M.c8 * v.x - M.c9 * v.y - M.c10 * v.z;
    ctype const c15 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z - M.c14 * v.w;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v, MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = v.w * M.c0;
    ctype const c5 = -v.x * M.c4 + v.w * M.c1;
    ctype const c6 = -v.y * M.c4 + v.w * M.c2;
    ctype const c7 = -v.z * M.c4 + v.w * M.c3;
    ctype const c8 = v.y * M.c3 - v.z * M.c2;
    ctype const c9 = -v.x * M.c3 + v.z * M.c1;
    ctype const c10 = v.x * M.c2 - v.y * M.c1;
    ctype const c11 = -v.y * M.c7 + v.z * M.c6 + v.w * M.c8;
    ctype const c12 = v.x * M.c7 - v.z * M.c5 + v.w * M.c9;
    ctype const c13 = -v.x * M.c6 + v.y * M.c5 + v.w * M.c10;
    ctype const c14 = -v.x * M.c8 - v.y * M.c9 - v.z * M.c10;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& M, Scalar3dp<U> s)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, MVec3dp<U> const& M)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp wdg :: wdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& A,
                                                  MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
                     A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& A,
                                                  MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4 + A.c2 * B.c2 - A.c3 * B.c1 + A.c4 * B.c3;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c2 + A.c3 * B.c0 + A.c5 * B.c3;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c1 - A.c2 * B.c0 + A.c6 * B.c3;
    ctype const c7 = A.c0 * B.c7 - A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp wdg :: wdg(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> wdg(MVec3dp_U<T> const& A,
                                                  MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = -A.c1 * B.c3 + A.c2 * B.c2 + A.c3 * B.c4 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c3 - A.c2 * B.c1 + A.c3 * B.c5 + A.c5 * B.c0;
    ctype const c6 = -A.c0 * B.c2 + A.c1 * B.c1 + A.c3 * B.c6 + A.c6 * B.c0;
    ctype const c7 = -A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 + A.c7 * B.c0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp wdg :: wdg(mv_e,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& M,
                                                   PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(M.c0 * ctype(ps));
}

// pga3dp wdg :: wdg(ps,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(PScalar3dp<T> ps,
                                                   MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * M.c0);
}

// pga3dp wdg :: wdg(mv_e,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& M,
                                                  TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * t.x;
    ctype const c1 = M.c0 * t.y;
    ctype const c2 = M.c0 * t.z;
    ctype const c3 = M.c0 * t.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(trivec,mv_e) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t,
                                                  MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * M.c0;
    ctype const c1 = t.y * M.c0;
    ctype const c2 = t.z * M.c0;
    ctype const c3 = t.w * M.c0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& M,
                                                  BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * B.vx;
    ctype const c2 = M.c0 * B.vy;
    ctype const c3 = M.c0 * B.vz;
    ctype const c4 = M.c0 * B.mx;
    ctype const c5 = M.c0 * B.my;
    ctype const c6 = M.c0 * B.mz;
    ctype const c7 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz - M.c4 * B.vx -
                     M.c5 * B.vy - M.c6 * B.vz;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B,
                                                  MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = B.vx * M.c0;
    ctype const c2 = B.vy * M.c0;
    ctype const c3 = B.vz * M.c0;
    ctype const c4 = B.mx * M.c0;
    ctype const c5 = B.my * M.c0;
    ctype const c6 = B.mz * M.c0;
    ctype const c7 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6 - B.mx * M.c1 -
                     B.my * M.c2 - B.mz * M.c3;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& M,
                                                  Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    ctype const c3 = M.c0 * v.w;
    ctype const c4 = M.c2 * v.z - M.c3 * v.y + M.c4 * v.w;
    ctype const c5 = -M.c1 * v.z + M.c3 * v.x + M.c5 * v.w;
    ctype const c6 = M.c1 * v.y - M.c2 * v.x + M.c6 * v.w;
    ctype const c7 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp wdg :: wdg(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v,
                                                  MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    ctype const c3 = v.w * M.c0;
    ctype const c4 = -v.y * M.c3 + v.z * M.c2 + v.w * M.c4;
    ctype const c5 = v.x * M.c3 - v.z * M.c1 + v.w * M.c5;
    ctype const c6 = -v.x * M.c2 + v.y * M.c1 + v.w * M.c6;
    ctype const c7 = -v.x * M.c4 - v.y * M.c5 - v.z * M.c6;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp wdg :: wdg(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(MVec3dp_E<T> const& M, Scalar3dp<U> s)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, MVec3dp_E<U> const& M)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(MVec3dp_U<T> const& A,
                                                  MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c2 = -A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c3 = -A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c4 = A.c1 * B.c2 - A.c2 * B.c1;
    ctype const c5 = -A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c6 = A.c0 * B.c1 - A.c1 * B.c0;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 - A.c4 * B.c0 -
                     A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(mv_u,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] MVec3dp_U<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(ps,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] MVec3dp_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(mv_u,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(MVec3dp_U<T> const& M,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(M.c0 * t.x + M.c1 * t.y + M.c2 * t.z + M.c3 * t.w);
}

// pga3dp wdg :: wdg(trivec,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t,
                                                   MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-t.x * M.c0 - t.y * M.c1 - t.z * M.c2 - t.w * M.c3);
}

// pga3dp wdg :: wdg(mv_u,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(MVec3dp_U<T> const& M,
                                                  BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.vz + M.c2 * B.vy + M.c3 * B.mx;
    ctype const c1 = M.c0 * B.vz - M.c2 * B.vx + M.c3 * B.my;
    ctype const c2 = -M.c0 * B.vy + M.c1 * B.vx + M.c3 * B.mz;
    ctype const c3 = -M.c0 * B.mx - M.c1 * B.my - M.c2 * B.mz;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(bivec,mv_u) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B,
                                                  MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * M.c2 - B.vz * M.c1 + B.mx * M.c3;
    ctype const c1 = -B.vx * M.c2 + B.vz * M.c0 + B.my * M.c3;
    ctype const c2 = B.vx * M.c1 - B.vy * M.c0 + B.mz * M.c3;
    ctype const c3 = -B.mx * M.c0 - B.my * M.c1 - B.mz * M.c2;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(MVec3dp_U<T> const& M,
                                                  Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -M.c0 * v.w + M.c3 * v.x;
    ctype const c2 = -M.c1 * v.w + M.c3 * v.y;
    ctype const c3 = -M.c2 * v.w + M.c3 * v.z;
    ctype const c4 = M.c1 * v.z - M.c2 * v.y;
    ctype const c5 = -M.c0 * v.z + M.c2 * v.x;
    ctype const c6 = M.c0 * v.y - M.c1 * v.x;
    ctype const c7 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z - M.c7 * v.w;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v,
                                                  MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -v.x * M.c3 + v.w * M.c0;
    ctype const c2 = -v.y * M.c3 + v.w * M.c1;
    ctype const c3 = -v.z * M.c3 + v.w * M.c2;
    ctype const c4 = v.y * M.c2 - v.z * M.c1;
    ctype const c5 = -v.x * M.c2 + v.z * M.c0;
    ctype const c6 = v.x * M.c1 - v.y * M.c0;
    ctype const c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6 + v.w * M.c7;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp wdg :: wdg(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> wdg(MVec3dp_U<T> const& M, Scalar3dp<U> s)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp wdg :: wdg(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, MVec3dp_U<U> const& M)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp wdg :: wdg(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dp<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3dp<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec3dp<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * ctype(s));
}

// pga3dp wdg :: wdg(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) * ctype(ps));
}

// pga3dp wdg :: wdg(trivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dp<T> const&,
                                                  [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(trivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dp<T> const&,
                                                  [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(bivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3dp<T> const&,
                                                  [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp wdg :: wdg(trivec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w);
}

// pga3dp wdg :: wdg(vec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w);
}

// pga3dp wdg :: wdg(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(bivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B1,
                                                   BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz -
                             B1.mx * B2.vx - B1.my * B2.vy - B1.mz * B2.vz);
}

// pga3dp wdg :: wdg(bivec,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B,
                                                  Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c1 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c2 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c3 = -B.mx * v.x - B.my * v.y - B.mz * v.z;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(vec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v,
                                                  BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.vz + v.z * B.vy + v.w * B.mx;
    ctype const c1 = v.x * B.vz - v.z * B.vx + v.w * B.my;
    ctype const c2 = -v.x * B.vy + v.y * B.vx + v.w * B.mz;
    ctype const c3 = -v.x * B.mx - v.y * B.my - v.z * B.mz;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp wdg :: wdg(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp wdg :: wdg(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c1 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c2 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c3 = v1.y * v2.z - v1.z * v2.y;
    ctype const c4 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c5 = v1.x * v2.y - v1.y * v2.x;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp wdg :: wdg(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp wdg :: wdg(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg for Point3d
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Plane3d<std::common_type_t<T, U>> wdg(Line3d<T> const& l, Point3d<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = l.vy * p.z - l.vz * p.y + l.mx;
    ctype const c1 = -l.vx * p.z + l.vz * p.x + l.my;
    ctype const c2 = l.vx * p.y - l.vy * p.x + l.mz;
    ctype const c3 = -l.mx * p.x - l.my * p.y - l.mz * p.z;
    return Plane3d<ctype>(c0, c1, c2, c3);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Plane3d<std::common_type_t<T, U>> wdg(Point3d<T> const& p, Line3d<U> const& l)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -p.y * l.vz + p.z * l.vy + l.mx;
    ctype const c1 = p.x * l.vz - p.z * l.vx + l.my;
    ctype const c2 = -p.x * l.vy + p.y * l.vx + l.mz;
    ctype const c3 = -p.x * l.mx - p.y * l.my - p.z * l.mz;
    return Plane3d<ctype>(c0, c1, c2, c3);
}

// wedge product between two points (aka vectors with implicit p.z == 1)
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Line3d<std::common_type_t<T, U>> wdg(Point3d<T> const& p, Point3d<U> const& q)
{
    return Line3d<std::common_type_t<T, U>>(q.x - p.x, q.y - p.y, q.z - p.z,
                                            p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z,
                                            p.x * q.y - p.y * q.x);
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg -> join
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> join(TriVec3dp<T> const& t,
                                                    Vec3dp<U> const& v)
{
    return wdg(t, v);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v,
                                                    TriVec3dp<U> const& t)
{
    return wdg(v, t);
}


template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> join(BiVec3dp<T> const& B,
                                                   Vec3dp<U> const& v)
{
    return wdg(B, v);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v,
                                                   BiVec3dp<U> const& B)
{
    return wdg(v, B);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Plane3d<std::common_type_t<T, U>> join(Line3d<T> const& l, Point3d<U> const& p)
{
    return wdg(l, p);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Plane3d<std::common_type_t<T, U>> join(Point3d<T> const& p, Line3d<U> const& l)
{
    return wdg(p, l);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Line3d<std::common_type_t<T, U>> join(Point3d<T> const& p, Point3d<U> const& q)
{
    return wdg(p, q);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v1,
                                                  Vec3dp<U> const& v2)
{
    return wdg(v1, v2);
}


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

// pga3dp rwdg :: rwdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 +
                     A.c4 * B.c14 - A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 -
                     A.c8 * B.c5 - A.c9 * B.c6 - A.c10 * B.c7 - A.c11 * B.c1 -
                     A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    ctype const c1 = A.c1 * B.c15 + A.c5 * B.c14 + A.c9 * B.c13 - A.c10 * B.c12 -
                     A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 + A.c15 * B.c1;
    ctype const c2 = A.c2 * B.c15 + A.c6 * B.c14 - A.c8 * B.c13 + A.c10 * B.c11 +
                     A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c6 + A.c15 * B.c2;
    ctype const c3 = A.c3 * B.c15 + A.c7 * B.c14 + A.c8 * B.c12 - A.c9 * B.c11 -
                     A.c11 * B.c9 + A.c12 * B.c8 + A.c14 * B.c7 + A.c15 * B.c3;
    ctype const c4 = A.c4 * B.c15 - A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 -
                     A.c11 * B.c5 - A.c12 * B.c6 - A.c13 * B.c7 + A.c15 * B.c4;
    ctype const c5 = A.c5 * B.c15 - A.c12 * B.c13 + A.c13 * B.c12 + A.c15 * B.c5;
    ctype const c6 = A.c6 * B.c15 + A.c11 * B.c13 - A.c13 * B.c11 + A.c15 * B.c6;
    ctype const c7 = A.c7 * B.c15 - A.c11 * B.c12 + A.c12 * B.c11 + A.c15 * B.c7;
    ctype const c8 = A.c8 * B.c15 + A.c11 * B.c14 - A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c9 = A.c9 * B.c15 + A.c12 * B.c14 - A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c10 = A.c10 * B.c15 + A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c11 = A.c11 * B.c15 + A.c15 * B.c11;
    ctype const c12 = A.c12 * B.c15 + A.c15 * B.c12;
    ctype const c13 = A.c13 * B.c15 + A.c15 * B.c13;
    ctype const c14 = A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = A.c15 * B.c15;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& A,
                                                 MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 - A.c8 * B.c1 -
                     A.c9 * B.c2 - A.c10 * B.c3 + A.c15 * B.c0;
    ctype const c1 = A.c1 * B.c7 - A.c12 * B.c6 + A.c13 * B.c5 + A.c14 * B.c1;
    ctype const c2 = A.c2 * B.c7 + A.c11 * B.c6 - A.c13 * B.c4 + A.c14 * B.c2;
    ctype const c3 = A.c3 * B.c7 - A.c11 * B.c5 + A.c12 * B.c4 + A.c14 * B.c3;
    ctype const c4 = A.c4 * B.c7 - A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3;
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& A,
                                                 MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                     A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 + A.c7 * B.c0;
    ctype const c1 = A.c1 * B.c14 + A.c5 * B.c13 - A.c6 * B.c12 + A.c7 * B.c1;
    ctype const c2 = A.c2 * B.c14 - A.c4 * B.c13 + A.c6 * B.c11 + A.c7 * B.c2;
    ctype const c3 = A.c3 * B.c14 + A.c4 * B.c12 - A.c5 * B.c11 + A.c7 * B.c3;
    ctype const c4 = -A.c1 * B.c11 - A.c2 * B.c12 - A.c3 * B.c13 + A.c7 * B.c4;
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& A,
                                                 MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c7 -
                     A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c2 - A.c14 * B.c3;
    ctype const c1 = A.c5 * B.c7 + A.c9 * B.c6 - A.c10 * B.c5 + A.c15 * B.c0;
    ctype const c2 = A.c6 * B.c7 - A.c8 * B.c6 + A.c10 * B.c4 + A.c15 * B.c1;
    ctype const c3 = A.c7 * B.c7 + A.c8 * B.c5 - A.c9 * B.c4 + A.c15 * B.c2;
    ctype const c4 = -A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 + A.c15 * B.c3;
    ctype const c5 = -A.c12 * B.c6 + A.c13 * B.c5;
    ctype const c6 = A.c11 * B.c6 - A.c13 * B.c4;
    ctype const c7 = -A.c11 * B.c5 + A.c12 * B.c4;
    ctype const c8 = A.c11 * B.c7 - A.c14 * B.c4;
    ctype const c9 = A.c12 * B.c7 - A.c14 * B.c5;
    ctype const c10 = A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c11 = A.c15 * B.c4;
    ctype const c12 = A.c15 * B.c5;
    ctype const c13 = A.c15 * B.c6;
    ctype const c14 = A.c15 * B.c7;
    ctype const c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp_U<T> const& A,
                                                 MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 -
                     A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 - A.c7 * B.c4;
    ctype const c1 = A.c0 * B.c15 - A.c5 * B.c10 + A.c6 * B.c9 + A.c7 * B.c5;
    ctype const c2 = A.c1 * B.c15 + A.c4 * B.c10 - A.c6 * B.c8 + A.c7 * B.c6;
    ctype const c3 = A.c2 * B.c15 - A.c4 * B.c9 + A.c5 * B.c8 + A.c7 * B.c7;
    ctype const c4 = A.c3 * B.c15 - A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7;
    ctype const c5 = -A.c5 * B.c13 + A.c6 * B.c12;
    ctype const c6 = A.c4 * B.c13 - A.c6 * B.c11;
    ctype const c7 = -A.c4 * B.c12 + A.c5 * B.c11;
    ctype const c8 = A.c4 * B.c14 - A.c7 * B.c11;
    ctype const c9 = A.c5 * B.c14 - A.c7 * B.c12;
    ctype const c10 = A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c11 = A.c4 * B.c15;
    ctype const c12 = A.c5 * B.c15;
    ctype const c13 = A.c6 * B.c15;
    ctype const c14 = A.c7 * B.c15;
    ctype const c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& M, PScalar3dp<U> ps)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps, MVec3dp<U> const& M)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& M,
                                                 TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * t.x + M.c2 * t.y + M.c3 * t.z + M.c4 * t.w;
    ctype const c1 = M.c5 * t.w + M.c9 * t.z - M.c10 * t.y;
    ctype const c2 = M.c6 * t.w - M.c8 * t.z + M.c10 * t.x;
    ctype const c3 = M.c7 * t.w + M.c8 * t.y - M.c9 * t.x;
    ctype const c4 = -M.c5 * t.x - M.c6 * t.y - M.c7 * t.z;
    ctype const c5 = -M.c12 * t.z + M.c13 * t.y;
    ctype const c6 = M.c11 * t.z - M.c13 * t.x;
    ctype const c7 = -M.c11 * t.y + M.c12 * t.x;
    ctype const c8 = M.c11 * t.w - M.c14 * t.x;
    ctype const c9 = M.c12 * t.w - M.c14 * t.y;
    ctype const c10 = M.c13 * t.w - M.c14 * t.z;
    ctype const c11 = M.c15 * t.x;
    ctype const c12 = M.c15 * t.y;
    ctype const c13 = M.c15 * t.z;
    ctype const c14 = M.c15 * t.w;
    ctype const c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                 MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3 - t.w * M.c4;
    ctype const c1 = -t.y * M.c10 + t.z * M.c9 + t.w * M.c5;
    ctype const c2 = t.x * M.c10 - t.z * M.c8 + t.w * M.c6;
    ctype const c3 = -t.x * M.c9 + t.y * M.c8 + t.w * M.c7;
    ctype const c4 = -t.x * M.c5 - t.y * M.c6 - t.z * M.c7;
    ctype const c5 = -t.y * M.c13 + t.z * M.c12;
    ctype const c6 = t.x * M.c13 - t.z * M.c11;
    ctype const c7 = -t.x * M.c12 + t.y * M.c11;
    ctype const c8 = t.x * M.c14 - t.w * M.c11;
    ctype const c9 = t.y * M.c14 - t.w * M.c12;
    ctype const c10 = t.z * M.c14 - t.w * M.c13;
    ctype const c11 = t.x * M.c15;
    ctype const c12 = t.y * M.c15;
    ctype const c13 = t.z * M.c15;
    ctype const c14 = t.w * M.c15;
    ctype const c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& M,
                                                 BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                     M.c9 * B.vy - M.c10 * B.vz;
    ctype const c1 = -M.c12 * B.mz + M.c13 * B.my + M.c14 * B.vx;
    ctype const c2 = M.c11 * B.mz - M.c13 * B.mx + M.c14 * B.vy;
    ctype const c3 = -M.c11 * B.my + M.c12 * B.mx + M.c14 * B.vz;
    ctype const c4 = -M.c11 * B.vx - M.c12 * B.vy - M.c13 * B.vz;
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B,
                                                 MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                     B.my * M.c6 - B.mz * M.c7;
    ctype const c1 = B.vx * M.c14 + B.my * M.c13 - B.mz * M.c12;
    ctype const c2 = B.vy * M.c14 - B.mx * M.c13 + B.mz * M.c11;
    ctype const c3 = B.vz * M.c14 + B.mx * M.c12 - B.my * M.c11;
    ctype const c4 = -B.vx * M.c11 - B.vy * M.c12 - B.vz * M.c13;
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& M, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z - M.c14 * v.w;
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(Vec3dp<T> const& v, MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14;
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rwdg :: rwdg(mv,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& M, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(M.c15 * ctype(s));
}

// pga3dp rwdg :: rwdg(s,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(Scalar3dp<T> s, MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s) * M.c15);
}

// pga3dp rwdg :: rwdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& A,
                                                   MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
                     A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c1 * B.c7 + A.c7 * B.c1;
    ctype const c2 = A.c2 * B.c7 + A.c7 * B.c2;
    ctype const c3 = A.c3 * B.c7 + A.c7 * B.c3;
    ctype const c4 = A.c4 * B.c7 + A.c7 * B.c4;
    ctype const c5 = A.c5 * B.c7 + A.c7 * B.c5;
    ctype const c6 = A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = A.c7 * B.c7;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& A,
                                                   MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c7 + A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c0;
    ctype const c1 = A.c2 * B.c7 - A.c4 * B.c6 + A.c6 * B.c4 + A.c7 * B.c1;
    ctype const c2 = A.c3 * B.c7 + A.c4 * B.c5 - A.c5 * B.c4 + A.c7 * B.c2;
    ctype const c3 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 + A.c7 * B.c3;
    ctype const c4 = A.c7 * B.c4;
    ctype const c5 = A.c7 * B.c5;
    ctype const c6 = A.c7 * B.c6;
    ctype const c7 = A.c7 * B.c7;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rwdg :: rwdg(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rwdg(MVec3dp_U<T> const& A,
                                                   MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype const c1 = A.c1 * B.c7 + A.c4 * B.c6 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype const c2 = A.c2 * B.c7 - A.c4 * B.c5 + A.c5 * B.c4 + A.c7 * B.c3;
    ctype const c3 = A.c3 * B.c7 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype const c4 = A.c4 * B.c7;
    ctype const c5 = A.c5 * B.c7;
    ctype const c6 = A.c6 * B.c7;
    ctype const c7 = A.c7 * B.c7;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rwdg :: rwdg(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& M,
                                                   PScalar3dp<U> ps)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps,
                                                   MVec3dp_E<U> const& M)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& M,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * t.w + M.c5 * t.z - M.c6 * t.y;
    ctype const c1 = M.c2 * t.w - M.c4 * t.z + M.c6 * t.x;
    ctype const c2 = M.c3 * t.w + M.c4 * t.y - M.c5 * t.x;
    ctype const c3 = -M.c1 * t.x - M.c2 * t.y - M.c3 * t.z;
    ctype const c4 = M.c7 * t.x;
    ctype const c5 = M.c7 * t.y;
    ctype const c6 = M.c7 * t.z;
    ctype const c7 = M.c7 * t.w;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rwdg :: rwdg(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                   MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * M.c6 + t.z * M.c5 + t.w * M.c1;
    ctype const c1 = t.x * M.c6 - t.z * M.c4 + t.w * M.c2;
    ctype const c2 = -t.x * M.c5 + t.y * M.c4 + t.w * M.c3;
    ctype const c3 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3;
    ctype const c4 = t.x * M.c7;
    ctype const c5 = t.y * M.c7;
    ctype const c6 = t.z * M.c7;
    ctype const c7 = t.w * M.c7;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rwdg :: rwdg(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& M,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz - M.c4 * B.vx -
                     M.c5 * B.vy - M.c6 * B.vz;
    ctype const c1 = M.c7 * B.vx;
    ctype const c2 = M.c7 * B.vy;
    ctype const c3 = M.c7 * B.vz;
    ctype const c4 = M.c7 * B.mx;
    ctype const c5 = M.c7 * B.my;
    ctype const c6 = M.c7 * B.mz;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B,
                                                   MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6 - B.mx * M.c1 -
                     B.my * M.c2 - B.mz * M.c3;
    ctype const c1 = B.vx * M.c7;
    ctype const c2 = B.vy * M.c7;
    ctype const c3 = B.vz * M.c7;
    ctype const c4 = B.mx * M.c7;
    ctype const c5 = B.my * M.c7;
    ctype const c6 = B.mz * M.c7;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(mv_e,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& M, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c7 * v.x;
    ctype const c1 = M.c7 * v.y;
    ctype const c2 = M.c7 * v.z;
    ctype const c3 = M.c7 * v.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(vec,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(Vec3dp<T> const& v, MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c7;
    ctype const c1 = v.y * M.c7;
    ctype const c2 = v.z * M.c7;
    ctype const c3 = v.w * M.c7;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(mv_e,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(MVec3dp_E<T> const& M, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(M.c7 * ctype(s));
}

// pga3dp rwdg :: rwdg(s,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(Scalar3dp<T> s, MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s) * M.c7);
}

// pga3dp rwdg :: rwdg(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(MVec3dp_U<T> const& A,
                                                   MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 - A.c4 * B.c0 -
                     A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3;
    ctype const c1 = -A.c5 * B.c6 + A.c6 * B.c5;
    ctype const c2 = A.c4 * B.c6 - A.c6 * B.c4;
    ctype const c3 = -A.c4 * B.c5 + A.c5 * B.c4;
    ctype const c4 = A.c4 * B.c7 - A.c7 * B.c4;
    ctype const c5 = A.c5 * B.c7 - A.c7 * B.c5;
    ctype const c6 = A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rwdg(MVec3dp_U<T> const& M,
                                                   PScalar3dp<U> ps)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rwdg :: rwdg(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps,
                                                   MVec3dp_U<U> const& M)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rwdg :: rwdg(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(MVec3dp_U<T> const& M,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * t.x + M.c1 * t.y + M.c2 * t.z + M.c3 * t.w;
    ctype const c1 = -M.c5 * t.z + M.c6 * t.y;
    ctype const c2 = M.c4 * t.z - M.c6 * t.x;
    ctype const c3 = -M.c4 * t.y + M.c5 * t.x;
    ctype const c4 = M.c4 * t.w - M.c7 * t.x;
    ctype const c5 = M.c5 * t.w - M.c7 * t.y;
    ctype const c6 = M.c6 * t.w - M.c7 * t.z;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                   MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c0 - t.y * M.c1 - t.z * M.c2 - t.w * M.c3;
    ctype const c1 = -t.y * M.c6 + t.z * M.c5;
    ctype const c2 = t.x * M.c6 - t.z * M.c4;
    ctype const c3 = -t.x * M.c5 + t.y * M.c4;
    ctype const c4 = t.x * M.c7 - t.w * M.c4;
    ctype const c5 = t.y * M.c7 - t.w * M.c5;
    ctype const c6 = t.z * M.c7 - t.w * M.c6;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rwdg :: rwdg(mv_u,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(MVec3dp_U<T> const& M,
                                                BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * B.mz + M.c6 * B.my + M.c7 * B.vx;
    ctype const c1 = M.c4 * B.mz - M.c6 * B.mx + M.c7 * B.vy;
    ctype const c2 = -M.c4 * B.my + M.c5 * B.mx + M.c7 * B.vz;
    ctype const c3 = -M.c4 * B.vx - M.c5 * B.vy - M.c6 * B.vz;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(bivec,mv_u) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B,
                                                MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * M.c7 + B.my * M.c6 - B.mz * M.c5;
    ctype const c1 = B.vy * M.c7 - B.mx * M.c6 + B.mz * M.c4;
    ctype const c2 = B.vz * M.c7 + B.mx * M.c5 - B.my * M.c4;
    ctype const c3 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(mv_u,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(MVec3dp_U<T> const& M,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(-M.c4 * v.x - M.c5 * v.y - M.c6 * v.z - M.c7 * v.w);
}

// pga3dp rwdg :: rwdg(vec,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(Vec3dp<T> const& v,
                                                   MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v.x * M.c4 + v.y * M.c5 + v.z * M.c6 + v.w * M.c7);
}

// pga3dp rwdg :: rwdg(mv_u,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] MVec3dp_U<T> const&,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(s,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] MVec3dp_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps1, PScalar3dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps1) * ctype(ps2));
}

// pga3dp rwdg :: rwdg(ps,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * t.x;
    ctype const c1 = ctype(ps) * t.y;
    ctype const c2 = ctype(ps) * t.z;
    ctype const c3 = ctype(ps) * t.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(trivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                   PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(ps);
    ctype const c1 = t.y * ctype(ps);
    ctype const c2 = t.z * ctype(ps);
    ctype const c3 = t.w * ctype(ps);
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.vx;
    ctype const c1 = ctype(ps) * B.vy;
    ctype const c2 = ctype(ps) * B.vz;
    ctype const c3 = ctype(ps) * B.mx;
    ctype const c4 = ctype(ps) * B.my;
    ctype const c5 = ctype(ps) * B.mz;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rwdg :: rwdg(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(ps);
    ctype const c1 = B.vy * ctype(ps);
    ctype const c2 = B.vz * ctype(ps);
    ctype const c3 = B.mx * ctype(ps);
    ctype const c4 = B.my * ctype(ps);
    ctype const c5 = B.mz * ctype(ps);
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rwdg :: rwdg(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = ctype(ps) * v.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(Vec3dp<T> const& v, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = v.w * ctype(ps);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(ps) * ctype(s));
}

// pga3dp rwdg :: rwdg(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s) * ctype(ps));
}

// pga3dp rwdg :: rwdg(trivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t1.y * t2.z + t1.z * t2.y;
    ctype const c1 = t1.x * t2.z - t1.z * t2.x;
    ctype const c2 = -t1.x * t2.y + t1.y * t2.x;
    ctype const c3 = t1.x * t2.w - t1.w * t2.x;
    ctype const c4 = t1.y * t2.w - t1.w * t2.y;
    ctype const c5 = t1.z * t2.w - t1.w * t2.z;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rwdg :: rwdg(trivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.mz + t.z * B.my + t.w * B.vx;
    ctype const c1 = t.x * B.mz - t.z * B.mx + t.w * B.vy;
    ctype const c2 = -t.x * B.my + t.y * B.mx + t.w * B.vz;
    ctype const c3 = -t.x * B.vx - t.y * B.vy - t.z * B.vz;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(bivec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B,
                                                TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * t.w + B.my * t.z - B.mz * t.y;
    ctype const c1 = B.vy * t.w - B.mx * t.z + B.mz * t.x;
    ctype const c2 = B.vz * t.w + B.mx * t.y - B.my * t.x;
    ctype const c3 = -B.vx * t.x - B.vy * t.y - B.vz * t.z;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rwdg :: rwdg(trivec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w);
}

// pga3dp rwdg :: rwdg(vec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w);
}

// pga3dp rwdg :: rwdg(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] TriVec3dp<T> const&,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B1,
                                                   BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(-B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz -
                            B1.mx * B2.vx - B1.my * B2.vy - B1.mz * B2.vz);
}

// pga3dp rwdg :: rwdg(bivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec3dp<T> const&,
                                                   [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(vec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dp<T> const&,
                                                   [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec3dp<T> const&,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(vec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dp<T> const&,
                                                   [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dp<T> const&,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rwdg :: rwdg(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions rwdg -> meet
////////////////////////////////////////////////////////////////////////////////

// intersection of two planes -> return a line
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> meet(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    return rwdg(t1, t2);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Line3d<std::common_type_t<T, U>> meet(Plane3d<T> const& p1,
                                                Plane3d<U> const& p2)
{
    return rwdg(static_cast<TriVec3dp<T> const&>(p1),
                static_cast<TriVec3dp<T> const&>(p2));
}

// intersection of plane and line -> return a point
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(TriVec3dp<T> const& t,
                                                BiVec3dp<U> const& B)
{
    return rwdg(t, B);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(BiVec3dp<T> const& B,
                                                TriVec3dp<U> const& t)
{
    return rwdg(B, t);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(Plane3d<T> const& p, Line3d<U> const& l)
{
    return rwdg(static_cast<TriVec3dp<T> const&>(p), static_cast<BiVec3dp<T> const&>(l));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(Line3d<T> const& l, Plane3d<U> const& p)
{
    return rwdg(static_cast<BiVec3dp<T> const&>(l), static_cast<TriVec3dp<T> const&>(p));
}


////////////////////////////////////////////////////////////////////////////////
// transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// twdg1(ps,vec) = trivec -> identical to geometric product gpr(ps,vec)
//                        -> and identical to the left contraction vec << ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> twdg1(PScalar3dp<T> ps, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.x;
    ctype const c1 = -ctype(ps) * v.y;
    ctype const c2 = -ctype(ps) * v.z;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// twdg1(vec,ps) = trivec -> identical to geometric product gpr(vec,ps)
//                        -> and identical to the right contraction ps >> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> twdg1(Vec3dp<T> const& v, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// twdg1(trivec,vec) = bivector -> identical to left contraction operator<<(vec,trivec)
//                              -> identical to right contraction operator>>(trivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> twdg1(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * v.z + t.z * v.y;
    ctype const c1 = t.x * v.z - t.z * v.x;
    ctype const c2 = -t.x * v.y + t.y * v.x;
    ctype const c3 = -t.w * v.x;
    ctype const c4 = -t.w * v.y;
    ctype const c5 = -t.w * v.z;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// twdg1(vec,trivec) = bivector -> identical to right contraction operator>>(trivec,vec)
//                              -> identical to left contraction operator<<(vec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> twdg1(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * t.z - v.z * t.y;
    ctype const c1 = -v.x * t.z + v.z * t.x;
    ctype const c2 = v.x * t.y - v.y * t.x;
    ctype const c3 = -v.x * t.w;
    ctype const c4 = -v.y * t.w;
    ctype const c5 = -v.z * t.w;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// twdg1(trivec,bivec) = trivector -> identical to commutator product cmt(bivec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> twdg1(TriVec3dp<T> const& t,
                                                    BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.mz + t.z * B.my + t.w * B.vx;
    ctype const c1 = t.x * B.mz - t.z * B.mx + t.w * B.vy;
    ctype const c2 = -t.x * B.my + t.y * B.mx + t.w * B.vz;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// twdg1(bivec,trivec) = trivector -> identical to commutator product cmt(trivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> twdg1(BiVec3dp<T> const& B,
                                                    TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * t.w - B.my * t.z + B.mz * t.y;
    ctype const c1 = -B.vy * t.w + B.mx * t.z - B.mz * t.x;
    ctype const c2 = -B.vz * t.w - B.mx * t.y + B.my * t.x;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}


// twdg1(bivec,bivec) = bivector -> identical to commutator product cmt(bivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> twdg1(BiVec3dp<T> const& B1,
                                                   BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz + B1.mz * B2.vy;
    ctype const c1 = B1.vx * B2.mz - B1.vz * B2.mx + B1.mx * B2.vz - B1.mz * B2.vx;
    ctype const c2 = -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx;
    ctype const c3 = -B1.my * B2.mz + B1.mz * B2.my;
    ctype const c4 = B1.mx * B2.mz - B1.mz * B2.mx;
    ctype const c5 = -B1.mx * B2.my + B1.my * B2.mx;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// twdg1(bivec,vec) = vector -> identical to left contraction operator<<(vec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> twdg1(BiVec3dp<T> const& B, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.my * v.z + B.mz * v.y;
    ctype const c1 = B.mx * v.z - B.mz * v.x;
    ctype const c2 = -B.mx * v.y + B.my * v.x;
    ctype const c3 = B.vx * v.x + B.vy * v.y + B.vz * v.z;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// twdg1(vec,bivec) = vector -> identical to right contraction operator>>(bivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> twdg1(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my;
    ctype const c1 = v.x * B.mz - v.z * B.mx;
    ctype const c2 = -v.x * B.my + v.y * B.mx;
    ctype const c3 = -v.x * B.vx - v.y * B.vy - v.z * B.vz;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// twdg1(vec,vec) = scalar -> identical to dot product dot(vec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> twdg1(Vec3dp<T> const& v1,
                                                    Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}


////////////////////////////////////////////////////////////////////////////////
// regressive transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// rtwdg1(trivec,trivec) = pseudoscalar -> identical to regressive dot product
// rdot(trivec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rtwdg1(TriVec3dp<T> const& t1,
                                                      TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(t1.x * t2.x + t1.y * t2.y + t1.z * t2.z);
}

// rtwdg1(trivec,bivec) = trivector -> identical to regressive commutator product
// rcmt(trivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rtwdg1(TriVec3dp<T> const& t,
                                                     BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * B.vz - t.z * B.vy;
    ctype const c1 = -t.x * B.vz + t.z * B.vx;
    ctype const c2 = t.x * B.vy - t.y * B.vx;
    ctype const c3 = t.x * B.mx + t.y * B.my + t.z * B.mz;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// rtwdg1(bivec,trivec) = trivector -> identical to regressive commutator product
// rcmt(bivec,trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rtwdg1(BiVec3dp<T> const& B,
                                                     TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * t.z - B.vz * t.y;
    ctype const c1 = -B.vx * t.z + B.vz * t.x;
    ctype const c2 = B.vx * t.y - B.vy * t.x;
    ctype const c3 = -B.mx * t.x - B.my * t.y - B.mz * t.z;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// rtwdg1(trivec,vec) = bivec -> identical to r_weight_expand3dp(vec,trivec)
//                            ->            = vec ^ r_weight_dual(trivec)
//                            ->            = l_weight_dual(trivec) ^ vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rtwdg1(TriVec3dp<T> const& t,
                                                    Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * v.w;
    ctype const c1 = -t.y * v.w;
    ctype const c2 = -t.z * v.w;
    ctype const c3 = t.y * v.z - t.z * v.y;
    ctype const c4 = -t.x * v.z + t.z * v.x;
    ctype const c5 = t.x * v.y - t.y * v.x;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// rtwdg1(vec,trivec) = bivec -> identical to l_weight_expand3dp(trivec,vec)
//                            ->            = l_weight_dual(trivec) ^ vec
//                            ->            = vec ^ r_weight_dual(trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rtwdg1(Vec3dp<T> const& v,
                                                    TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.w * t.x;
    ctype const c1 = -v.w * t.y;
    ctype const c2 = -v.w * t.z;
    ctype const c3 = -v.y * t.z + v.z * t.y;
    ctype const c4 = v.x * t.z - v.z * t.x;
    ctype const c5 = -v.x * t.y + v.y * t.x;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// rtwdg1(trivec,s) = vec -> identical to r_weight_expand3dp(s,trivec)
//                        ->        = s ^ r_weight_dual(trivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rtwdg1(TriVec3dp<T> const& t, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * ctype(s);
    ctype const c1 = -t.y * ctype(s);
    ctype const c2 = -t.z * ctype(s);
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// rtwdg1(s,trivec) -> identical to l_weight_expand3dp(trivec,s)
//                  ->  = l_weight_dual(trivec) ^ s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rtwdg1(Scalar3dp<T> s, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// rtwdg1(bivec,bivec) = bivec -> identical to regressive commutator product
// rcmt(bivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rtwdg1(BiVec3dp<T> const& B1,
                                                    BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vy * B2.vz - B1.vz * B2.vy;
    ctype const c1 = -B1.vx * B2.vz + B1.vz * B2.vx;
    ctype const c2 = B1.vx * B2.vy - B1.vy * B2.vx;
    ctype const c3 = B1.vy * B2.mz - B1.vz * B2.my + B1.my * B2.vz - B1.mz * B2.vy;
    ctype const c4 = -B1.vx * B2.mz + B1.vz * B2.mx - B1.mx * B2.vz + B1.mz * B2.vx;
    ctype const c5 = B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// rtwdg1(bivec,vec) = vec ->            = rcmt(bivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rtwdg1(BiVec3dp<T> const& B,
                                                  Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c1 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c2 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// rtwdg1(vec,bivec) = vec ->            = rcmt(vec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rtwdg1(Vec3dp<T> const& v,
                                                  BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.vz - v.z * B.vy - v.w * B.mx;
    ctype const c1 = -v.x * B.vz + v.z * B.vx - v.w * B.my;
    ctype const c2 = v.x * B.vy - v.y * B.vx - v.w * B.mz;
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// Implements the left bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator<<(a,b) = l_bulk_contract(a,b) = rwdg( l_bulk_dual(a), b )
//
////////////////////////////////////////////////////////////////////////////////

// pga3dp l_bulk_contract :: l_bulk_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator<<(MVec3dp<T> const& A,
                                                       MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c8 * B.c8 +
                     A.c9 * B.c9 + A.c10 * B.c10 + A.c14 * B.c14;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c10 - A.c3 * B.c9 - A.c8 * B.c14;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c10 + A.c3 * B.c8 - A.c9 * B.c14;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c9 - A.c2 * B.c8 - A.c10 * B.c14;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 +
                     A.c8 * B.c11 + A.c9 * B.c12 + A.c10 * B.c13 + A.c14 * B.c15;
    ctype const c5 = A.c0 * B.c5 + A.c2 * B.c13 - A.c3 * B.c12 - A.c8 * B.c15;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c13 + A.c3 * B.c11 - A.c9 * B.c15;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c12 - A.c2 * B.c11 - A.c10 * B.c15;
    ctype const c8 = A.c0 * B.c8 - A.c1 * B.c14;
    ctype const c9 = A.c0 * B.c9 - A.c2 * B.c14;
    ctype const c10 = A.c0 * B.c10 - A.c3 * B.c14;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c15;
    ctype const c12 = A.c0 * B.c12 - A.c2 * B.c15;
    ctype const c13 = A.c0 * B.c13 - A.c3 * B.c15;
    ctype const c14 = A.c0 * B.c14;
    ctype const c15 = A.c0 * B.c15;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp l_bulk_contract :: l_bulk_contract(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dp<T>, [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(trivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(TriVec3dp<T> const& t,
                                                      PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = t.w * ctype(ps);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp l_bulk_contract :: l_bulk_contract(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dp<T>, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.mx * ctype(ps);
    ctype const c1 = -B.my * ctype(ps);
    ctype const c2 = -B.mz * ctype(ps);
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp l_bulk_contract :: l_bulk_contract(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dp<T>, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(vec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                         PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * ctype(ps);
    ctype const c1 = -v.y * ctype(ps);
    ctype const c2 = -v.z * ctype(ps);
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp l_bulk_contract :: l_bulk_contract(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s,
                                                          PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) * ctype(ps));
}

// pga3dp l_bulk_contract :: l_bulk_contract(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(t1.w * t2.w);
}

// pga3dp l_bulk_contract :: l_bulk_contract(trivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(bivec,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B,
                                                      TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.mx * t.w;
    ctype const c1 = -B.my * t.w;
    ctype const c2 = -B.mz * t.w;
    ctype const c3 = B.mx * t.x + B.my * t.y + B.mz * t.z;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp l_bulk_contract :: l_bulk_contract(trivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(vec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * t.z - v.z * t.y;
    ctype const c1 = -v.x * t.z + v.z * t.x;
    ctype const c2 = v.x * t.y - v.y * t.x;
    ctype const c3 = -v.x * t.w;
    ctype const c4 = -v.y * t.w;
    ctype const c5 = -v.z * t.w;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp l_bulk_contract :: l_bulk_contract(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s,
                                                         TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp l_bulk_contract :: l_bulk_contract(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B1,
                                                         BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(B1.mx * B2.mx + B1.my * B2.my + B1.mz * B2.mz);
}

// pga3dp l_bulk_contract :: l_bulk_contract(bivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// (v << B) == cmt(B,v)
// pga3dp l_bulk_contract :: l_bulk_contract(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                      BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.mz - v.z * B.my;
    ctype const c1 = -v.x * B.mz + v.z * B.mx;
    ctype const c2 = v.x * B.my - v.y * B.mx;
    ctype const c3 = v.x * B.vx + v.y * B.vy + v.z * B.vz;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp l_bulk_contract :: l_bulk_contract(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp l_bulk_contract :: l_bulk_contract(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v1,
                                                         Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// pga3dp l_bulk_contract :: l_bulk_contract(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp l_bulk_contract :: l_bulk_contract(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp l_bulk_contract :: l_bulk_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The resulting object lies in A and is perpendicular to B
//
// Implements the right bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator>>(a,b) = r_bulk_contract(a,b) = rwdg( a, r_bulk_dual(b) )
//
////////////////////////////////////////////////////////////////////////////////

// pga3dp r_bulk_contract :: r_bulk_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator>>(MVec3dp<T> const& A,
                                                       MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c8 * B.c8 +
                     A.c9 * B.c9 + A.c10 * B.c10 + A.c14 * B.c14;
    ctype const c1 = A.c1 * B.c0 + A.c9 * B.c3 - A.c10 * B.c2 - A.c14 * B.c8;
    ctype const c2 = A.c2 * B.c0 - A.c8 * B.c3 + A.c10 * B.c1 - A.c14 * B.c9;
    ctype const c3 = A.c3 * B.c0 + A.c8 * B.c2 - A.c9 * B.c1 - A.c14 * B.c10;
    ctype const c4 = A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 +
                     A.c11 * B.c8 + A.c12 * B.c9 + A.c13 * B.c10 - A.c15 * B.c14;
    ctype const c5 = A.c5 * B.c0 - A.c12 * B.c3 + A.c13 * B.c2 - A.c15 * B.c8;
    ctype const c6 = A.c6 * B.c0 + A.c11 * B.c3 - A.c13 * B.c1 - A.c15 * B.c9;
    ctype const c7 = A.c7 * B.c0 - A.c11 * B.c2 + A.c12 * B.c1 - A.c15 * B.c10;
    ctype const c8 = A.c8 * B.c0 - A.c14 * B.c1;
    ctype const c9 = A.c9 * B.c0 - A.c14 * B.c2;
    ctype const c10 = A.c10 * B.c0 - A.c14 * B.c3;
    ctype const c11 = A.c11 * B.c0 + A.c15 * B.c1;
    ctype const c12 = A.c12 * B.c0 + A.c15 * B.c2;
    ctype const c13 = A.c13 * B.c0 + A.c15 * B.c3;
    ctype const c14 = A.c14 * B.c0;
    ctype const c15 = A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp r_bulk_contract :: r_bulk_contract(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(ps,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                      TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = -ctype(ps) * t.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp r_bulk_contract :: r_bulk_contract(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.mx;
    ctype const c1 = -ctype(ps) * B.my;
    ctype const c2 = -ctype(ps) * B.mz;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp r_bulk_contract :: r_bulk_contract(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(ps,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                         Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp r_bulk_contract :: r_bulk_contract(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                          Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * ctype(s));
}

// pga3dp r_bulk_contract :: r_bulk_contract(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar3dp<T>,
                                                         [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(trivec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(t1.w * t2.w);
}

// pga3dp r_bulk_contract :: r_bulk_contract(trivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                      BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.w * B.mx;
    ctype const c1 = -t.w * B.my;
    ctype const c2 = -t.w * B.mz;
    ctype const c3 = t.x * B.mx + t.y * B.my + t.z * B.mz;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp r_bulk_contract :: r_bulk_contract(bivec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(trivec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * v.z + t.z * v.y;
    ctype const c1 = t.x * v.z - t.z * v.x;
    ctype const c2 = -t.x * v.y + t.y * v.x;
    ctype const c3 = -t.w * v.x;
    ctype const c4 = -t.w * v.y;
    ctype const c5 = -t.w * v.z;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp r_bulk_contract :: r_bulk_contract(vec,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                         Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp r_bulk_contract :: r_bulk_contract(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dp<T>, [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(BiVec3dp<T> const& B1,
                                                         BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(B1.mx * B2.mx + B1.my * B2.my + B1.mz * B2.mz);
}

// (B >> v) == cmt(v,B)
// pga3dp r_bulk_contract :: r_bulk_contract(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(BiVec3dp<T> const& B,
                                                      Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.my * v.z - B.mz * v.y;
    ctype const c1 = -B.mx * v.z + B.mz * v.x;
    ctype const c2 = B.mx * v.y - B.my * v.x;
    ctype const c3 = -B.vx * v.x - B.vy * v.y - B.vz * v.z;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp r_bulk_contract :: r_bulk_contract(vec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(BiVec3dp<T> const& B,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp r_bulk_contract :: r_bulk_contract(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dp<T>, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(Vec3dp<T> const& v1,
                                                         Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// pga3dp r_bulk_contract :: r_bulk_contract(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp r_bulk_contract :: r_bulk_contract(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dp<T>, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp r_bulk_contract :: r_bulk_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

// pga3dp cmt :: cmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> cmt(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c2 * B.c10 + A.c3 * B.c9 - A.c9 * B.c3 + A.c10 * B.c2;
    ctype const c2 = A.c1 * B.c10 - A.c3 * B.c8 + A.c8 * B.c3 - A.c10 * B.c1;
    ctype const c3 = -A.c1 * B.c9 + A.c2 * B.c8 - A.c8 * B.c2 + A.c9 * B.c1;
    ctype const c4 = -A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 + A.c5 * B.c1 +
                     A.c6 * B.c2 + A.c7 * B.c3 + A.c14 * B.c15 - A.c15 * B.c14;
    ctype const c5 = -A.c1 * B.c4 + A.c4 * B.c1 - A.c6 * B.c10 + A.c7 * B.c9 -
                     A.c9 * B.c7 + A.c10 * B.c6 + A.c11 * B.c14 - A.c14 * B.c11;
    ctype const c6 = -A.c2 * B.c4 + A.c4 * B.c2 + A.c5 * B.c10 - A.c7 * B.c8 +
                     A.c8 * B.c7 - A.c10 * B.c5 + A.c12 * B.c14 - A.c14 * B.c12;
    ctype const c7 = -A.c3 * B.c4 + A.c4 * B.c3 - A.c5 * B.c9 + A.c6 * B.c8 -
                     A.c8 * B.c6 + A.c9 * B.c5 + A.c13 * B.c14 - A.c14 * B.c13;
    ctype const c8 = A.c2 * B.c3 - A.c3 * B.c2 - A.c9 * B.c10 + A.c10 * B.c9;
    ctype const c9 = -A.c1 * B.c3 + A.c3 * B.c1 + A.c8 * B.c10 - A.c10 * B.c8;
    ctype const c10 = A.c1 * B.c2 - A.c2 * B.c1 - A.c8 * B.c9 + A.c9 * B.c8;
    ctype const c11 = A.c1 * B.c15 - A.c5 * B.c14 - A.c9 * B.c13 + A.c10 * B.c12 -
                      A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 - A.c15 * B.c1;
    ctype const c12 = A.c2 * B.c15 - A.c6 * B.c14 + A.c8 * B.c13 - A.c10 * B.c11 +
                      A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c6 - A.c15 * B.c2;
    ctype const c13 = A.c3 * B.c15 - A.c7 * B.c14 - A.c8 * B.c12 + A.c9 * B.c11 -
                      A.c11 * B.c9 + A.c12 * B.c8 + A.c14 * B.c7 - A.c15 * B.c3;
    ctype const c14 = 0.0;
    ctype const c15 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
                      A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp cmt :: cmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(ps,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> cmt(PScalar3dp<T> ps, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = -ctype(ps) * t.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp cmt :: cmt(trivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> cmt(TriVec3dp<T> const& t, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = t.w * ctype(ps);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp cmt :: cmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec3dp<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(ps,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> cmt(PScalar3dp<T> ps, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.x;
    ctype const c1 = -ctype(ps) * v.y;
    ctype const c2 = -ctype(ps) * v.z;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp cmt :: cmt(vec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> cmt(Vec3dp<T> const& v, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp cmt :: cmt(ps,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3dp<T>,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(s,ps) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(trivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> cmt(TriVec3dp<T> const& t1,
                                                 TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t1.x * t2.w - t1.w * t2.x;
    ctype const c1 = t1.y * t2.w - t1.w * t2.y;
    ctype const c2 = t1.z * t2.w - t1.w * t2.z;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp cmt :: cmt(trivec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> cmt(TriVec3dp<T> const& t,
                                                  BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.mz + t.z * B.my + t.w * B.vx;
    ctype const c1 = t.x * B.mz - t.z * B.mx + t.w * B.vy;
    ctype const c2 = -t.x * B.my + t.y * B.mx + t.w * B.vz;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp cmt :: cmt(bivec,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& B,
                                                  TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * t.w - B.my * t.z + B.mz * t.y;
    ctype const c1 = -B.vy * t.w + B.mx * t.z - B.mz * t.x;
    ctype const c2 = -B.vz * t.w - B.mx * t.y + B.my * t.x;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp cmt :: cmt(trivec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> cmt(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w);
}

// pga3dp cmt :: cmt(vec,trivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> cmt(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w);
}

// pga3dp cmt :: cmt(trivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] TriVec3dp<T> const&,
                                                  [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(s,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dp<T>,
                                                  [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(bivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& B1,
                                                 BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz + B1.mz * B2.vy;
    ctype const c1 = B1.vx * B2.mz - B1.vz * B2.mx + B1.mx * B2.vz - B1.mz * B2.vx;
    ctype const c2 = -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx;
    ctype const c3 = -B1.my * B2.mz + B1.mz * B2.my;
    ctype const c4 = B1.mx * B2.mz - B1.mz * B2.mx;
    ctype const c5 = -B1.mx * B2.my + B1.my * B2.mx;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// cmt(B,v) == -cmt(v,B) == (v << B)
// pga3dp cmt :: cmt(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& B, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.my * v.z + B.mz * v.y;
    ctype const c1 = B.mx * v.z - B.mz * v.x;
    ctype const c2 = -B.mx * v.y + B.my * v.x;
    ctype const c3 = B.vx * v.x + B.vy * v.y + B.vz * v.z;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// cmt(v,B) = -cmt(B,v) == (B >> v)
// pga3dp cmt :: cmt(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> cmt(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my;
    ctype const c1 = v.x * B.mz - v.z * B.mx;
    ctype const c2 = -v.x * B.my + v.y * B.mx;
    ctype const c3 = -v.x * B.vx - v.y * B.vy - v.z * B.vz;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp cmt :: cmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec3dp<T> const&,
                                                  [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dp<T>,
                                                  [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> cmt(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c1 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c2 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c3 = v1.y * v2.z - v1.z * v2.y;
    ctype const c4 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c5 = v1.x * v2.y - v1.y * v2.x;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp cmt :: cmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec3dp<T> const&,
                                                  [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dp<T>,
                                                  [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp cmt :: cmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3dp<T>,
                                                  [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////

// pga3dp rcmt :: rcmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rcmt(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
                     A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4;
    ctype const c1 = A.c0 * B.c11 + A.c2 * B.c7 - A.c3 * B.c6 - A.c4 * B.c8 +
                     A.c6 * B.c3 - A.c7 * B.c2 + A.c8 * B.c4 - A.c11 * B.c0;
    ctype const c2 = A.c0 * B.c12 - A.c1 * B.c7 + A.c3 * B.c5 - A.c4 * B.c9 -
                     A.c5 * B.c3 + A.c7 * B.c1 + A.c9 * B.c4 - A.c12 * B.c0;
    ctype const c3 = A.c0 * B.c13 + A.c1 * B.c6 - A.c2 * B.c5 - A.c4 * B.c10 +
                     A.c5 * B.c2 - A.c6 * B.c1 + A.c10 * B.c4 - A.c13 * B.c0;
    ctype const c4 = 0.0;
    ctype const c5 = A.c6 * B.c7 - A.c7 * B.c6 - A.c12 * B.c13 + A.c13 * B.c12;
    ctype const c6 = -A.c5 * B.c7 + A.c7 * B.c5 + A.c11 * B.c13 - A.c13 * B.c11;
    ctype const c7 = A.c5 * B.c6 - A.c6 * B.c5 - A.c11 * B.c12 + A.c12 * B.c11;
    ctype const c8 = -A.c1 * B.c4 + A.c4 * B.c1 + A.c6 * B.c10 - A.c7 * B.c9 +
                     A.c9 * B.c7 - A.c10 * B.c6 + A.c11 * B.c14 - A.c14 * B.c11;
    ctype const c9 = -A.c2 * B.c4 + A.c4 * B.c2 - A.c5 * B.c10 + A.c7 * B.c8 -
                     A.c8 * B.c7 + A.c10 * B.c5 + A.c12 * B.c14 - A.c14 * B.c12;
    ctype const c10 = -A.c3 * B.c4 + A.c4 * B.c3 + A.c5 * B.c9 - A.c6 * B.c8 +
                      A.c8 * B.c6 - A.c9 * B.c5 + A.c13 * B.c14 - A.c14 * B.c13;
    ctype const c11 = A.c6 * B.c13 - A.c7 * B.c12 + A.c12 * B.c7 - A.c13 * B.c6;
    ctype const c12 = -A.c5 * B.c13 + A.c7 * B.c11 - A.c11 * B.c7 + A.c13 * B.c5;
    ctype const c13 = A.c5 * B.c12 - A.c6 * B.c11 + A.c11 * B.c6 - A.c12 * B.c5;
    ctype const c14 = A.c0 * B.c4 - A.c4 * B.c0 - A.c8 * B.c11 - A.c9 * B.c12 -
                      A.c10 * B.c13 + A.c11 * B.c8 + A.c12 * B.c9 + A.c13 * B.c10;
    ctype const c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rcmt :: rcmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dp<T>,
                                                   [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(ps,trivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dp<T>,
                                                   [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(trivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] TriVec3dp<T> const&,
                                                   [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dp<T>,
                                                   [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec3dp<T> const&,
                                                   [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dp<T>,
                                                   [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Vec3dp<T> const&,
                                                   [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar3dp<T>,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(trivec,trivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rcmt(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t1.y * t2.z + t1.z * t2.y;
    ctype const c1 = t1.x * t2.z - t1.z * t2.x;
    ctype const c2 = -t1.x * t2.y + t1.y * t2.x;
    ctype const c3 = t1.x * t2.w - t1.w * t2.x;
    ctype const c4 = t1.y * t2.w - t1.w * t2.y;
    ctype const c5 = t1.z * t2.w - t1.w * t2.z;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rcmt :: rcmt(trivec,bivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rcmt(TriVec3dp<T> const& t,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.y * B.vz - t.z * B.vy;
    ctype const c1 = -t.x * B.vz + t.z * B.vx;
    ctype const c2 = t.x * B.vy - t.y * B.vx;
    ctype const c3 = t.x * B.mx + t.y * B.my + t.z * B.mz;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(bivec,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rcmt(BiVec3dp<T> const& B,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * t.z - B.vz * t.y;
    ctype const c1 = -B.vx * t.z + B.vz * t.x;
    ctype const c2 = B.vx * t.y - B.vy * t.x;
    ctype const c3 = -B.mx * t.x - B.my * t.y - B.mz * t.z;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(trivec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w);
}

// pga3dp rcmt :: rcmt(vec,trivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w);
}

// pga3dp rcmt :: rcmt(trivec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rcmt(TriVec3dp<T> const& t, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * ctype(s);
    ctype const c1 = -t.y * ctype(s);
    ctype const c2 = -t.z * ctype(s);
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(s,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rcmt(Scalar3dp<T> s, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(bivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rcmt(BiVec3dp<T> const& B1,
                                                  BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.vy * B2.vz - B1.vz * B2.vy;
    ctype const c1 = -B1.vx * B2.vz + B1.vz * B2.vx;
    ctype const c2 = B1.vx * B2.vy - B1.vy * B2.vx;
    ctype const c3 = B1.vy * B2.mz - B1.vz * B2.my + B1.my * B2.vz - B1.mz * B2.vy;
    ctype const c4 = -B1.vx * B2.mz + B1.vz * B2.mx - B1.mx * B2.vz + B1.mz * B2.vx;
    ctype const c5 = B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rcmt :: rcmt(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rcmt(BiVec3dp<T> const& B, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c1 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c2 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rcmt(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.vz - v.z * B.vy - v.w * B.mx;
    ctype const c1 = -v.x * B.vz + v.z * B.vx - v.w * B.my;
    ctype const c2 = v.x * B.vy - v.y * B.vx - v.w * B.mz;
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec3dp<T> const&,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp rcmt :: rcmt(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rcmt(Vec3dp<T> const& v1,
                                                  Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c4 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c5 = -v1.z * v2.w + v1.w * v2.z;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rcmt :: rcmt(vec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rcmt(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = -v.w * ctype(s);
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(s,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rcmt(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = ctype(s) * v.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rcmt :: rcmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 3dp multivector
// REALLY VERY Expensive! - Don't use if you don't have to!

// pga3dp gpr :: gpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                      MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c8 * B.c8 -
                     A.c9 * B.c9 - A.c10 * B.c10 - A.c14 * B.c14;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c10 + A.c3 * B.c9 +
                     A.c8 * B.c14 - A.c9 * B.c3 + A.c10 * B.c2 + A.c14 * B.c8;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c10 + A.c2 * B.c0 - A.c3 * B.c8 +
                     A.c8 * B.c3 + A.c9 * B.c14 - A.c10 * B.c1 + A.c14 * B.c9;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c0 - A.c8 * B.c2 +
                     A.c9 * B.c1 + A.c10 * B.c14 + A.c14 * B.c10;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 + A.c4 * B.c0 +
                     A.c5 * B.c1 + A.c6 * B.c2 + A.c7 * B.c3 - A.c8 * B.c11 -
                     A.c9 * B.c12 - A.c10 * B.c13 - A.c11 * B.c8 - A.c12 * B.c9 -
                     A.c13 * B.c10 + A.c14 * B.c15 - A.c15 * B.c14;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c4 + A.c2 * B.c13 - A.c3 * B.c12 +
                     A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c10 + A.c7 * B.c9 +
                     A.c8 * B.c15 - A.c9 * B.c7 + A.c10 * B.c6 + A.c11 * B.c14 -
                     A.c12 * B.c3 + A.c13 * B.c2 - A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c13 - A.c2 * B.c4 + A.c3 * B.c11 +
                     A.c4 * B.c2 + A.c5 * B.c10 + A.c6 * B.c0 - A.c7 * B.c8 +
                     A.c8 * B.c7 + A.c9 * B.c15 - A.c10 * B.c5 + A.c11 * B.c3 +
                     A.c12 * B.c14 - A.c13 * B.c1 - A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c12 - A.c2 * B.c11 - A.c3 * B.c4 +
                     A.c4 * B.c3 - A.c5 * B.c9 + A.c6 * B.c8 + A.c7 * B.c0 - A.c8 * B.c6 +
                     A.c9 * B.c5 + A.c10 * B.c15 - A.c11 * B.c2 + A.c12 * B.c1 +
                     A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c8 = A.c0 * B.c8 - A.c1 * B.c14 + A.c2 * B.c3 - A.c3 * B.c2 +
                     A.c8 * B.c0 - A.c9 * B.c10 + A.c10 * B.c9 - A.c14 * B.c1;
    ctype const c9 = A.c0 * B.c9 - A.c1 * B.c3 - A.c2 * B.c14 + A.c3 * B.c1 +
                     A.c8 * B.c10 + A.c9 * B.c0 - A.c10 * B.c8 - A.c14 * B.c2;
    ctype const c10 = A.c0 * B.c10 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c14 -
                      A.c8 * B.c9 + A.c9 * B.c8 + A.c10 * B.c0 - A.c14 * B.c3;
    ctype const c11 = A.c0 * B.c11 + A.c1 * B.c15 - A.c2 * B.c7 + A.c3 * B.c6 +
                      A.c4 * B.c8 - A.c5 * B.c14 + A.c6 * B.c3 - A.c7 * B.c2 +
                      A.c8 * B.c4 - A.c9 * B.c13 + A.c10 * B.c12 + A.c11 * B.c0 -
                      A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 - A.c15 * B.c1;
    ctype const c12 = A.c0 * B.c12 + A.c1 * B.c7 + A.c2 * B.c15 - A.c3 * B.c5 +
                      A.c4 * B.c9 - A.c5 * B.c3 - A.c6 * B.c14 + A.c7 * B.c1 +
                      A.c8 * B.c13 + A.c9 * B.c4 - A.c10 * B.c11 + A.c11 * B.c10 +
                      A.c12 * B.c0 - A.c13 * B.c8 + A.c14 * B.c6 - A.c15 * B.c2;
    ctype const c13 = A.c0 * B.c13 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c15 +
                      A.c4 * B.c10 + A.c5 * B.c2 - A.c6 * B.c1 - A.c7 * B.c14 -
                      A.c8 * B.c12 + A.c9 * B.c11 + A.c10 * B.c4 - A.c11 * B.c9 +
                      A.c12 * B.c8 + A.c13 * B.c0 + A.c14 * B.c7 - A.c15 * B.c3;
    ctype const c14 = A.c0 * B.c14 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                      A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 +
                      A.c4 * B.c14 - A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 -
                      A.c8 * B.c5 - A.c9 * B.c6 - A.c10 * B.c7 - A.c11 * B.c1 -
                      A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                      MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c8 * B.c4 - A.c9 * B.c5 - A.c10 * B.c6;
    ctype const c1 = A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 + A.c14 * B.c4;
    ctype const c2 = A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c14 * B.c5;
    ctype const c3 = -A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 + A.c14 * B.c6;
    ctype const c4 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c4 * B.c0 -
                     A.c11 * B.c4 - A.c12 * B.c5 - A.c13 * B.c6 + A.c14 * B.c7;
    ctype const c5 = A.c0 * B.c1 + A.c5 * B.c0 - A.c6 * B.c6 + A.c7 * B.c5 + A.c8 * B.c7 -
                     A.c9 * B.c3 + A.c10 * B.c2 + A.c15 * B.c4;
    ctype const c6 = A.c0 * B.c2 + A.c5 * B.c6 + A.c6 * B.c0 - A.c7 * B.c4 + A.c8 * B.c3 +
                     A.c9 * B.c7 - A.c10 * B.c1 + A.c15 * B.c5;
    ctype const c7 = A.c0 * B.c3 - A.c5 * B.c5 + A.c6 * B.c4 + A.c7 * B.c0 - A.c8 * B.c2 +
                     A.c9 * B.c1 + A.c10 * B.c7 + A.c15 * B.c6;
    ctype const c8 = A.c0 * B.c4 + A.c8 * B.c0 - A.c9 * B.c6 + A.c10 * B.c5;
    ctype const c9 = A.c0 * B.c5 + A.c8 * B.c6 + A.c9 * B.c0 - A.c10 * B.c4;
    ctype const c10 = A.c0 * B.c6 - A.c8 * B.c5 + A.c9 * B.c4 + A.c10 * B.c0;
    ctype const c11 = A.c1 * B.c7 - A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c4 +
                      A.c11 * B.c0 - A.c12 * B.c6 + A.c13 * B.c5 + A.c14 * B.c1;
    ctype const c12 = A.c1 * B.c3 + A.c2 * B.c7 - A.c3 * B.c1 + A.c4 * B.c5 +
                      A.c11 * B.c6 + A.c12 * B.c0 - A.c13 * B.c4 + A.c14 * B.c2;
    ctype const c13 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c7 + A.c4 * B.c6 -
                      A.c11 * B.c5 + A.c12 * B.c4 + A.c13 * B.c0 + A.c14 * B.c3;
    ctype const c14 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 + A.c14 * B.c0;
    ctype const c15 = A.c0 * B.c7 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 -
                      A.c8 * B.c1 - A.c9 * B.c2 - A.c10 * B.c3 + A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                      MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c4 * B.c8 - A.c5 * B.c9 - A.c6 * B.c10;
    ctype const c1 = A.c0 * B.c1 + A.c4 * B.c14 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype const c2 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c14 - A.c6 * B.c1;
    ctype const c3 = A.c0 * B.c3 - A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c14;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 -
                     A.c4 * B.c11 - A.c5 * B.c12 - A.c6 * B.c13 - A.c7 * B.c14;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c0 - A.c2 * B.c10 + A.c3 * B.c9 +
                     A.c4 * B.c15 - A.c5 * B.c7 + A.c6 * B.c6 + A.c7 * B.c8;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c10 + A.c2 * B.c0 - A.c3 * B.c8 +
                     A.c4 * B.c7 + A.c5 * B.c15 - A.c6 * B.c5 + A.c7 * B.c9;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c0 - A.c4 * B.c6 +
                     A.c5 * B.c5 + A.c6 * B.c15 + A.c7 * B.c10;
    ctype const c8 = A.c0 * B.c8 + A.c4 * B.c0 - A.c5 * B.c10 + A.c6 * B.c9;
    ctype const c9 = A.c0 * B.c9 + A.c4 * B.c10 + A.c5 * B.c0 - A.c6 * B.c8;
    ctype const c10 = A.c0 * B.c10 - A.c4 * B.c9 + A.c5 * B.c8 + A.c6 * B.c0;
    ctype const c11 = A.c0 * B.c11 - A.c1 * B.c14 + A.c2 * B.c3 - A.c3 * B.c2 +
                      A.c4 * B.c4 - A.c5 * B.c13 + A.c6 * B.c12 - A.c7 * B.c1;
    ctype const c12 = A.c0 * B.c12 - A.c1 * B.c3 - A.c2 * B.c14 + A.c3 * B.c1 +
                      A.c4 * B.c13 + A.c5 * B.c4 - A.c6 * B.c11 - A.c7 * B.c2;
    ctype const c13 = A.c0 * B.c13 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c14 -
                      A.c4 * B.c12 + A.c5 * B.c11 + A.c6 * B.c4 - A.c7 * B.c3;
    ctype const c14 = A.c0 * B.c14 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype const c15 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                      A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 + A.c7 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                      MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c14 * B.c7;
    ctype const c1 = A.c0 * B.c0 + A.c8 * B.c7 - A.c9 * B.c2 + A.c10 * B.c1;
    ctype const c2 = A.c0 * B.c1 + A.c8 * B.c2 + A.c9 * B.c7 - A.c10 * B.c0;
    ctype const c3 = A.c0 * B.c2 - A.c8 * B.c1 + A.c9 * B.c0 + A.c10 * B.c7;
    ctype const c4 = A.c0 * B.c3 + A.c5 * B.c0 + A.c6 * B.c1 + A.c7 * B.c2 - A.c8 * B.c4 -
                     A.c9 * B.c5 - A.c10 * B.c6 - A.c15 * B.c7;
    ctype const c5 = -A.c1 * B.c3 + A.c2 * B.c6 - A.c3 * B.c5 + A.c4 * B.c0 +
                     A.c11 * B.c7 - A.c12 * B.c2 + A.c13 * B.c1 - A.c14 * B.c4;
    ctype const c6 = -A.c1 * B.c6 - A.c2 * B.c3 + A.c3 * B.c4 + A.c4 * B.c1 +
                     A.c11 * B.c2 + A.c12 * B.c7 - A.c13 * B.c0 - A.c14 * B.c5;
    ctype const c7 = A.c1 * B.c5 - A.c2 * B.c4 - A.c3 * B.c3 + A.c4 * B.c2 -
                     A.c11 * B.c1 + A.c12 * B.c0 + A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c8 = -A.c1 * B.c7 + A.c2 * B.c2 - A.c3 * B.c1 - A.c14 * B.c0;
    ctype const c9 = -A.c1 * B.c2 - A.c2 * B.c7 + A.c3 * B.c0 - A.c14 * B.c1;
    ctype const c10 = A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c7 - A.c14 * B.c2;
    ctype const c11 = A.c0 * B.c4 - A.c5 * B.c7 + A.c6 * B.c2 - A.c7 * B.c1 +
                      A.c8 * B.c3 - A.c9 * B.c6 + A.c10 * B.c5 - A.c15 * B.c0;
    ctype const c12 = A.c0 * B.c5 - A.c5 * B.c2 - A.c6 * B.c7 + A.c7 * B.c0 +
                      A.c8 * B.c6 + A.c9 * B.c3 - A.c10 * B.c4 - A.c15 * B.c1;
    ctype const c13 = A.c0 * B.c6 + A.c5 * B.c1 - A.c6 * B.c0 - A.c7 * B.c7 -
                      A.c8 * B.c5 + A.c9 * B.c4 + A.c10 * B.c3 - A.c15 * B.c2;
    ctype const c14 = A.c0 * B.c7 - A.c8 * B.c0 - A.c9 * B.c1 - A.c10 * B.c2;
    ctype const c15 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c7 -
                      A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c2 - A.c14 * B.c3;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                      MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c7 * B.c14;
    ctype const c1 = A.c0 * B.c0 - A.c1 * B.c10 + A.c2 * B.c9 + A.c7 * B.c8;
    ctype const c2 = A.c0 * B.c10 + A.c1 * B.c0 - A.c2 * B.c8 + A.c7 * B.c9;
    ctype const c3 = -A.c0 * B.c9 + A.c1 * B.c8 + A.c2 * B.c0 + A.c7 * B.c10;
    ctype const c4 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c7 + A.c3 * B.c0 -
                     A.c4 * B.c8 - A.c5 * B.c9 - A.c6 * B.c10 + A.c7 * B.c15;
    ctype const c5 = -A.c0 * B.c4 + A.c1 * B.c13 - A.c2 * B.c12 + A.c3 * B.c1 +
                     A.c4 * B.c14 - A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c11;
    ctype const c6 = -A.c0 * B.c13 - A.c1 * B.c4 + A.c2 * B.c11 + A.c3 * B.c2 +
                     A.c4 * B.c3 + A.c5 * B.c14 - A.c6 * B.c1 - A.c7 * B.c12;
    ctype const c7 = A.c0 * B.c12 - A.c1 * B.c11 - A.c2 * B.c4 + A.c3 * B.c3 -
                     A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c8 = -A.c0 * B.c14 + A.c1 * B.c3 - A.c2 * B.c2 - A.c7 * B.c1;
    ctype const c9 = -A.c0 * B.c3 - A.c1 * B.c14 + A.c2 * B.c1 - A.c7 * B.c2;
    ctype const c10 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c14 - A.c7 * B.c3;
    ctype const c11 = A.c0 * B.c15 - A.c1 * B.c7 + A.c2 * B.c6 + A.c3 * B.c8 +
                      A.c4 * B.c0 - A.c5 * B.c10 + A.c6 * B.c9 + A.c7 * B.c5;
    ctype const c12 = A.c0 * B.c7 + A.c1 * B.c15 - A.c2 * B.c5 + A.c3 * B.c9 +
                      A.c4 * B.c10 + A.c5 * B.c0 - A.c6 * B.c8 + A.c7 * B.c6;
    ctype const c13 = -A.c0 * B.c6 + A.c1 * B.c5 + A.c2 * B.c15 + A.c3 * B.c10 -
                      A.c4 * B.c9 + A.c5 * B.c8 + A.c6 * B.c0 + A.c7 * B.c7;
    ctype const c14 = -A.c0 * B.c8 - A.c1 * B.c9 - A.c2 * B.c10 + A.c7 * B.c0;
    ctype const c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 -
                      A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 - A.c7 * B.c4;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& M,
                                                      PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = M.c14 * ctype(ps);
    ctype const c5 = M.c8 * ctype(ps);
    ctype const c6 = M.c9 * ctype(ps);
    ctype const c7 = M.c10 * ctype(ps);
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = M.c1 * ctype(ps);
    ctype const c12 = M.c2 * ctype(ps);
    ctype const c13 = M.c3 * ctype(ps);
    ctype const c14 = 0.0;
    ctype const c15 = M.c0 * ctype(ps);
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                      MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = -ctype(ps) * M.c14;
    ctype const c5 = ctype(ps) * M.c8;
    ctype const c6 = ctype(ps) * M.c9;
    ctype const c7 = ctype(ps) * M.c10;
    ctype const c8 = 0.0;
    ctype const c9 = 0.0;
    ctype const c10 = 0.0;
    ctype const c11 = -ctype(ps) * M.c1;
    ctype const c12 = -ctype(ps) * M.c2;
    ctype const c13 = -ctype(ps) * M.c3;
    ctype const c14 = 0.0;
    ctype const c15 = ctype(ps) * M.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& M,
                                                      TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c14 * t.w;
    ctype const c1 = M.c8 * t.w;
    ctype const c2 = M.c9 * t.w;
    ctype const c3 = M.c10 * t.w;
    ctype const c4 = -M.c8 * t.x - M.c9 * t.y - M.c10 * t.z - M.c15 * t.w;
    ctype const c5 = M.c2 * t.z - M.c3 * t.y + M.c11 * t.w - M.c14 * t.x;
    ctype const c6 = -M.c1 * t.z + M.c3 * t.x + M.c12 * t.w - M.c14 * t.y;
    ctype const c7 = M.c1 * t.y - M.c2 * t.x + M.c13 * t.w - M.c14 * t.z;
    ctype const c8 = -M.c1 * t.w;
    ctype const c9 = -M.c2 * t.w;
    ctype const c10 = -M.c3 * t.w;
    ctype const c11 = M.c0 * t.x - M.c5 * t.w - M.c9 * t.z + M.c10 * t.y;
    ctype const c12 = M.c0 * t.y - M.c6 * t.w + M.c8 * t.z - M.c10 * t.x;
    ctype const c13 = M.c0 * t.z - M.c7 * t.w - M.c8 * t.y + M.c9 * t.x;
    ctype const c14 = M.c0 * t.w;
    ctype const c15 = M.c1 * t.x + M.c2 * t.y + M.c3 * t.z + M.c4 * t.w;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                      MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.w * M.c14;
    ctype const c1 = t.w * M.c8;
    ctype const c2 = t.w * M.c9;
    ctype const c3 = t.w * M.c10;
    ctype const c4 = -t.x * M.c8 - t.y * M.c9 - t.z * M.c10 + t.w * M.c15;
    ctype const c5 = t.x * M.c14 - t.y * M.c3 + t.z * M.c2 - t.w * M.c11;
    ctype const c6 = t.x * M.c3 + t.y * M.c14 - t.z * M.c1 - t.w * M.c12;
    ctype const c7 = -t.x * M.c2 + t.y * M.c1 + t.z * M.c14 - t.w * M.c13;
    ctype const c8 = -t.w * M.c1;
    ctype const c9 = -t.w * M.c2;
    ctype const c10 = -t.w * M.c3;
    ctype const c11 = t.x * M.c0 - t.y * M.c10 + t.z * M.c9 + t.w * M.c5;
    ctype const c12 = t.x * M.c10 + t.y * M.c0 - t.z * M.c8 + t.w * M.c6;
    ctype const c13 = -t.x * M.c9 + t.y * M.c8 + t.z * M.c0 + t.w * M.c7;
    ctype const c14 = t.w * M.c0;
    ctype const c15 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3 - t.w * M.c4;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& M,
                                                      BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c8 * B.mx - M.c9 * B.my - M.c10 * B.mz;
    ctype const c1 = -M.c2 * B.mz + M.c3 * B.my + M.c14 * B.mx;
    ctype const c2 = M.c1 * B.mz - M.c3 * B.mx + M.c14 * B.my;
    ctype const c3 = -M.c1 * B.my + M.c2 * B.mx + M.c14 * B.mz;
    ctype const c4 = -M.c1 * B.vx - M.c2 * B.vy - M.c3 * B.vz - M.c11 * B.mx -
                     M.c12 * B.my - M.c13 * B.mz;
    ctype const c5 = M.c0 * B.vx - M.c6 * B.mz + M.c7 * B.my - M.c9 * B.vz +
                     M.c10 * B.vy + M.c15 * B.mx;
    ctype const c6 = M.c0 * B.vy + M.c5 * B.mz - M.c7 * B.mx + M.c8 * B.vz -
                     M.c10 * B.vx + M.c15 * B.my;
    ctype const c7 = M.c0 * B.vz - M.c5 * B.my + M.c6 * B.mx - M.c8 * B.vy + M.c9 * B.vx +
                     M.c15 * B.mz;
    ctype const c8 = M.c0 * B.mx - M.c9 * B.mz + M.c10 * B.my;
    ctype const c9 = M.c0 * B.my + M.c8 * B.mz - M.c10 * B.mx;
    ctype const c10 = M.c0 * B.mz - M.c8 * B.my + M.c9 * B.mx;
    ctype const c11 = -M.c2 * B.vz + M.c3 * B.vy + M.c4 * B.mx - M.c12 * B.mz +
                      M.c13 * B.my + M.c14 * B.vx;
    ctype const c12 = M.c1 * B.vz - M.c3 * B.vx + M.c4 * B.my + M.c11 * B.mz -
                      M.c13 * B.mx + M.c14 * B.vy;
    ctype const c13 = -M.c1 * B.vy + M.c2 * B.vx + M.c4 * B.mz - M.c11 * B.my +
                      M.c12 * B.mx + M.c14 * B.vz;
    ctype const c14 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz;
    ctype const c15 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                      M.c9 * B.vy - M.c10 * B.vz;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                      MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.mx * M.c8 - B.my * M.c9 - B.mz * M.c10;
    ctype const c1 = B.mx * M.c14 - B.my * M.c3 + B.mz * M.c2;
    ctype const c2 = B.mx * M.c3 + B.my * M.c14 - B.mz * M.c1;
    ctype const c3 = -B.mx * M.c2 + B.my * M.c1 + B.mz * M.c14;
    ctype const c4 = B.vx * M.c1 + B.vy * M.c2 + B.vz * M.c3 - B.mx * M.c11 -
                     B.my * M.c12 - B.mz * M.c13;
    ctype const c5 = B.vx * M.c0 - B.vy * M.c10 + B.vz * M.c9 + B.mx * M.c15 -
                     B.my * M.c7 + B.mz * M.c6;
    ctype const c6 = B.vx * M.c10 + B.vy * M.c0 - B.vz * M.c8 + B.mx * M.c7 +
                     B.my * M.c15 - B.mz * M.c5;
    ctype const c7 = -B.vx * M.c9 + B.vy * M.c8 + B.vz * M.c0 - B.mx * M.c6 +
                     B.my * M.c5 + B.mz * M.c15;
    ctype const c8 = B.mx * M.c0 - B.my * M.c10 + B.mz * M.c9;
    ctype const c9 = B.mx * M.c10 + B.my * M.c0 - B.mz * M.c8;
    ctype const c10 = -B.mx * M.c9 + B.my * M.c8 + B.mz * M.c0;
    ctype const c11 = -B.vx * M.c14 + B.vy * M.c3 - B.vz * M.c2 + B.mx * M.c4 -
                      B.my * M.c13 + B.mz * M.c12;
    ctype const c12 = -B.vx * M.c3 - B.vy * M.c14 + B.vz * M.c1 + B.mx * M.c13 +
                      B.my * M.c4 - B.mz * M.c11;
    ctype const c13 = B.vx * M.c2 - B.vy * M.c1 - B.vz * M.c14 - B.mx * M.c12 +
                      B.my * M.c11 + B.mz * M.c4;
    ctype const c14 = -B.mx * M.c1 - B.my * M.c2 - B.mz * M.c3;
    ctype const c15 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                      B.my * M.c6 - B.mz * M.c7;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& M,
                                                      Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    ctype const c1 = M.c0 * v.x - M.c9 * v.z + M.c10 * v.y;
    ctype const c2 = M.c0 * v.y + M.c8 * v.z - M.c10 * v.x;
    ctype const c3 = M.c0 * v.z - M.c8 * v.y + M.c9 * v.x;
    ctype const c4 = M.c0 * v.w + M.c5 * v.x + M.c6 * v.y + M.c7 * v.z;
    ctype const c5 = -M.c1 * v.w + M.c4 * v.x - M.c12 * v.z + M.c13 * v.y;
    ctype const c6 = -M.c2 * v.w + M.c4 * v.y + M.c11 * v.z - M.c13 * v.x;
    ctype const c7 = -M.c3 * v.w + M.c4 * v.z - M.c11 * v.y + M.c12 * v.x;
    ctype const c8 = M.c2 * v.z - M.c3 * v.y - M.c14 * v.x;
    ctype const c9 = -M.c1 * v.z + M.c3 * v.x - M.c14 * v.y;
    ctype const c10 = M.c1 * v.y - M.c2 * v.x - M.c14 * v.z;
    ctype const c11 = M.c6 * v.z - M.c7 * v.y + M.c8 * v.w - M.c15 * v.x;
    ctype const c12 = -M.c5 * v.z + M.c7 * v.x + M.c9 * v.w - M.c15 * v.y;
    ctype const c13 = M.c5 * v.y - M.c6 * v.x + M.c10 * v.w - M.c15 * v.z;
    ctype const c14 = -M.c8 * v.x - M.c9 * v.y - M.c10 * v.z;
    ctype const c15 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z - M.c14 * v.w;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                      MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    ctype const c1 = v.x * M.c0 - v.y * M.c10 + v.z * M.c9;
    ctype const c2 = v.x * M.c10 + v.y * M.c0 - v.z * M.c8;
    ctype const c3 = -v.x * M.c9 + v.y * M.c8 + v.z * M.c0;
    ctype const c4 = -v.x * M.c5 - v.y * M.c6 - v.z * M.c7 + v.w * M.c0;
    ctype const c5 = -v.x * M.c4 + v.y * M.c13 - v.z * M.c12 + v.w * M.c1;
    ctype const c6 = -v.x * M.c13 - v.y * M.c4 + v.z * M.c11 + v.w * M.c2;
    ctype const c7 = v.x * M.c12 - v.y * M.c11 - v.z * M.c4 + v.w * M.c3;
    ctype const c8 = -v.x * M.c14 + v.y * M.c3 - v.z * M.c2;
    ctype const c9 = -v.x * M.c3 - v.y * M.c14 + v.z * M.c1;
    ctype const c10 = v.x * M.c2 - v.y * M.c1 - v.z * M.c14;
    ctype const c11 = v.x * M.c15 - v.y * M.c7 + v.z * M.c6 + v.w * M.c8;
    ctype const c12 = v.x * M.c7 + v.y * M.c15 - v.z * M.c5 + v.w * M.c9;
    ctype const c13 = -v.x * M.c6 + v.y * M.c5 + v.z * M.c15 + v.w * M.c10;
    ctype const c14 = -v.x * M.c8 - v.y * M.c9 - v.z * M.c10;
    ctype const c15 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& M, Scalar3dp<U> s)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s, MVec3dp<U> const& M)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp gpr :: gpr(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                        MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c7 -
                     A.c5 * B.c3 + A.c6 * B.c2 + A.c7 * B.c4;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c4 * B.c3 +
                     A.c5 * B.c7 - A.c6 * B.c1 + A.c7 * B.c5;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
                     A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c4 = A.c0 * B.c4 + A.c4 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5;
    ctype const c5 = A.c0 * B.c5 + A.c4 * B.c6 + A.c5 * B.c0 - A.c6 * B.c4;
    ctype const c6 = A.c0 * B.c6 - A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
                     A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                        MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c4 * B.c7 - A.c5 * B.c2 + A.c6 * B.c1;
    ctype const c1 = A.c0 * B.c1 + A.c4 * B.c2 + A.c5 * B.c7 - A.c6 * B.c0;
    ctype const c2 = A.c0 * B.c2 - A.c4 * B.c1 + A.c5 * B.c0 + A.c6 * B.c7;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c4 * B.c4 -
                     A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c7 + A.c2 * B.c2 - A.c3 * B.c1 + A.c4 * B.c3 -
                     A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c0;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c2 - A.c2 * B.c7 + A.c3 * B.c0 + A.c4 * B.c6 +
                     A.c5 * B.c3 - A.c6 * B.c4 - A.c7 * B.c1;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c7 - A.c4 * B.c5 +
                     A.c5 * B.c4 + A.c6 * B.c3 - A.c7 * B.c2;
    ctype const c7 = A.c0 * B.c7 - A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                        MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c6 + A.c2 * B.c5 + A.c7 * B.c4;
    ctype const c1 = A.c0 * B.c6 + A.c1 * B.c0 - A.c2 * B.c4 + A.c7 * B.c5;
    ctype const c2 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c0 + A.c7 * B.c6;
    ctype const c3 = -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0 -
                     A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 + A.c7 * B.c7;
    ctype const c4 = A.c0 * B.c7 - A.c1 * B.c3 + A.c2 * B.c2 + A.c3 * B.c4 + A.c4 * B.c0 -
                     A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype const c5 = A.c0 * B.c3 + A.c1 * B.c7 - A.c2 * B.c1 + A.c3 * B.c5 + A.c4 * B.c6 +
                     A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype const c6 = -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c7 + A.c3 * B.c6 -
                     A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3;
    ctype const c7 = -A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 + A.c7 * B.c0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c4 * ctype(ps);
    ctype const c2 = A.c5 * ctype(ps);
    ctype const c3 = A.c6 * ctype(ps);
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = A.c0 * ctype(ps);
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                        MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = ctype(ps) * B.c4;
    ctype const c2 = ctype(ps) * B.c5;
    ctype const c3 = ctype(ps) * B.c6;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = ctype(ps) * B.c0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& M,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c4 * t.w;
    ctype const c1 = M.c5 * t.w;
    ctype const c2 = M.c6 * t.w;
    ctype const c3 = -M.c4 * t.x - M.c5 * t.y - M.c6 * t.z - M.c7 * t.w;
    ctype const c4 = M.c0 * t.x - M.c1 * t.w - M.c5 * t.z + M.c6 * t.y;
    ctype const c5 = M.c0 * t.y - M.c2 * t.w + M.c4 * t.z - M.c6 * t.x;
    ctype const c6 = M.c0 * t.z - M.c3 * t.w - M.c4 * t.y + M.c5 * t.x;
    ctype const c7 = M.c0 * t.w;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.w * M.c4;
    ctype const c1 = t.w * M.c5;
    ctype const c2 = t.w * M.c6;
    ctype const c3 = -t.x * M.c4 - t.y * M.c5 - t.z * M.c6 + t.w * M.c7;
    ctype const c4 = t.x * M.c0 - t.y * M.c6 + t.z * M.c5 + t.w * M.c1;
    ctype const c5 = t.x * M.c6 + t.y * M.c0 - t.z * M.c4 + t.w * M.c2;
    ctype const c6 = -t.x * M.c5 + t.y * M.c4 + t.z * M.c0 + t.w * M.c3;
    ctype const c7 = t.w * M.c0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& M,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * B.mx - M.c5 * B.my - M.c6 * B.mz;
    ctype const c1 =
        M.c0 * B.vx - M.c2 * B.mz + M.c3 * B.my - M.c5 * B.vz + M.c6 * B.vy + M.c7 * B.mx;
    ctype const c2 =
        M.c0 * B.vy + M.c1 * B.mz - M.c3 * B.mx + M.c4 * B.vz - M.c6 * B.vx + M.c7 * B.my;
    ctype const c3 =
        M.c0 * B.vz - M.c1 * B.my + M.c2 * B.mx - M.c4 * B.vy + M.c5 * B.vx + M.c7 * B.mz;
    ctype const c4 = M.c0 * B.mx - M.c5 * B.mz + M.c6 * B.my;
    ctype const c5 = M.c0 * B.my + M.c4 * B.mz - M.c6 * B.mx;
    ctype const c6 = M.c0 * B.mz - M.c4 * B.my + M.c5 * B.mx;
    ctype const c7 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz - M.c4 * B.vx -
                     M.c5 * B.vy - M.c6 * B.vz;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.mx * M.c4 - B.my * M.c5 - B.mz * M.c6;
    ctype const c1 =
        B.vx * M.c0 - B.vy * M.c6 + B.vz * M.c5 + B.mx * M.c7 - B.my * M.c3 + B.mz * M.c2;
    ctype const c2 =
        B.vx * M.c6 + B.vy * M.c0 - B.vz * M.c4 + B.mx * M.c3 + B.my * M.c7 - B.mz * M.c1;
    ctype const c3 = -B.vx * M.c5 + B.vy * M.c4 + B.vz * M.c0 - B.mx * M.c2 +
                     B.my * M.c1 + B.mz * M.c7;
    ctype const c4 = B.mx * M.c0 - B.my * M.c6 + B.mz * M.c5;
    ctype const c5 = B.mx * M.c6 + B.my * M.c0 - B.mz * M.c4;
    ctype const c6 = -B.mx * M.c5 + B.my * M.c4 + B.mz * M.c0;
    ctype const c7 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6 - B.mx * M.c1 -
                     B.my * M.c2 - B.mz * M.c3;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * v.x - A.c5 * v.z + A.c6 * v.y;
    ctype const c1 = A.c0 * v.y + A.c4 * v.z - A.c6 * v.x;
    ctype const c2 = A.c0 * v.z - A.c4 * v.y + A.c5 * v.x;
    ctype const c3 = A.c0 * v.w + A.c1 * v.x + A.c2 * v.y + A.c3 * v.z;
    ctype const c4 = A.c2 * v.z - A.c3 * v.y + A.c4 * v.w - A.c7 * v.x;
    ctype const c5 = -A.c1 * v.z + A.c3 * v.x + A.c5 * v.w - A.c7 * v.y;
    ctype const c6 = A.c1 * v.y - A.c2 * v.x + A.c6 * v.w - A.c7 * v.z;
    ctype const c7 = -A.c4 * v.x - A.c5 * v.y - A.c6 * v.z;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * B.c0 - v.y * B.c6 + v.z * B.c5;
    ctype const c1 = v.x * B.c6 + v.y * B.c0 - v.z * B.c4;
    ctype const c2 = -v.x * B.c5 + v.y * B.c4 + v.z * B.c0;
    ctype const c3 = -v.x * B.c1 - v.y * B.c2 - v.z * B.c3 + v.w * B.c0;
    ctype const c4 = v.x * B.c7 - v.y * B.c3 + v.z * B.c2 + v.w * B.c4;
    ctype const c5 = v.x * B.c3 + v.y * B.c7 - v.z * B.c1 + v.w * B.c5;
    ctype const c6 = -v.x * B.c2 + v.y * B.c1 + v.z * B.c7 + v.w * B.c6;
    ctype const c7 = -v.x * B.c4 - v.y * B.c5 - v.z * B.c6;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& M,
                                                        Scalar3dp<U> s)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                        MVec3dp_E<U> const& M)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                        MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c7 * B.c7;
    ctype const c1 = -A.c0 * B.c3 + A.c1 * B.c6 - A.c2 * B.c5 + A.c3 * B.c0 +
                     A.c4 * B.c7 - A.c5 * B.c2 + A.c6 * B.c1 - A.c7 * B.c4;
    ctype const c2 = -A.c0 * B.c6 - A.c1 * B.c3 + A.c2 * B.c4 + A.c3 * B.c1 +
                     A.c4 * B.c2 + A.c5 * B.c7 - A.c6 * B.c0 - A.c7 * B.c5;
    ctype const c3 = A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c3 + A.c3 * B.c2 - A.c4 * B.c1 +
                     A.c5 * B.c0 + A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c4 = -A.c0 * B.c7 + A.c1 * B.c2 - A.c2 * B.c1 - A.c7 * B.c0;
    ctype const c5 = -A.c0 * B.c2 - A.c1 * B.c7 + A.c2 * B.c0 - A.c7 * B.c1;
    ctype const c6 = A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c7 - A.c7 * B.c2;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 - A.c4 * B.c0 -
                     A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = A.c7 * ctype(ps);
    ctype const c4 = A.c0 * ctype(ps);
    ctype const c5 = A.c1 * ctype(ps);
    ctype const c6 = A.c2 * ctype(ps);
    ctype const c7 = 0.0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                        MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = -ctype(ps) * B.c7;
    ctype const c4 = -ctype(ps) * B.c0;
    ctype const c5 = -ctype(ps) * B.c1;
    ctype const c6 = -ctype(ps) * B.c2;
    ctype const c7 = 0.0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& M,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c7 * t.w;
    ctype const c1 = M.c1 * t.z - M.c2 * t.y + M.c4 * t.w - M.c7 * t.x;
    ctype const c2 = -M.c0 * t.z + M.c2 * t.x + M.c5 * t.w - M.c7 * t.y;
    ctype const c3 = M.c0 * t.y - M.c1 * t.x + M.c6 * t.w - M.c7 * t.z;
    ctype const c4 = -M.c0 * t.w;
    ctype const c5 = -M.c1 * t.w;
    ctype const c6 = -M.c2 * t.w;
    ctype const c7 = M.c0 * t.x + M.c1 * t.y + M.c2 * t.z + M.c3 * t.w;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.w * M.c7;
    ctype const c1 = t.x * M.c7 - t.y * M.c2 + t.z * M.c1 - t.w * M.c4;
    ctype const c2 = t.x * M.c2 + t.y * M.c7 - t.z * M.c0 - t.w * M.c5;
    ctype const c3 = -t.x * M.c1 + t.y * M.c0 + t.z * M.c7 - t.w * M.c6;
    ctype const c4 = -t.w * M.c0;
    ctype const c5 = -t.w * M.c1;
    ctype const c6 = -t.w * M.c2;
    ctype const c7 = -t.x * M.c0 - t.y * M.c1 - t.z * M.c2 - t.w * M.c3;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& M,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.mz + M.c2 * B.my + M.c7 * B.mx;
    ctype const c1 = M.c0 * B.mz - M.c2 * B.mx + M.c7 * B.my;
    ctype const c2 = -M.c0 * B.my + M.c1 * B.mx + M.c7 * B.mz;
    ctype const c3 = -M.c0 * B.vx - M.c1 * B.vy - M.c2 * B.vz - M.c4 * B.mx -
                     M.c5 * B.my - M.c6 * B.mz;
    ctype const c4 = -M.c1 * B.vz + M.c2 * B.vy + M.c3 * B.mx - M.c5 * B.mz +
                     M.c6 * B.my + M.c7 * B.vx;
    ctype const c5 =
        M.c0 * B.vz - M.c2 * B.vx + M.c3 * B.my + M.c4 * B.mz - M.c6 * B.mx + M.c7 * B.vy;
    ctype const c6 = -M.c0 * B.vy + M.c1 * B.vx + M.c3 * B.mz - M.c4 * B.my +
                     M.c5 * B.mx + M.c7 * B.vz;
    ctype const c7 = -M.c0 * B.mx - M.c1 * B.my - M.c2 * B.mz;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.mx * M.c7 - B.my * M.c2 + B.mz * M.c1;
    ctype const c1 = B.mx * M.c2 + B.my * M.c7 - B.mz * M.c0;
    ctype const c2 = -B.mx * M.c1 + B.my * M.c0 + B.mz * M.c7;
    ctype const c3 =
        B.vx * M.c0 + B.vy * M.c1 + B.vz * M.c2 - B.mx * M.c4 - B.my * M.c5 - B.mz * M.c6;
    ctype const c4 = -B.vx * M.c7 + B.vy * M.c2 - B.vz * M.c1 + B.mx * M.c3 -
                     B.my * M.c6 + B.mz * M.c5;
    ctype const c5 = -B.vx * M.c2 - B.vy * M.c7 + B.vz * M.c0 + B.mx * M.c6 +
                     B.my * M.c3 - B.mz * M.c4;
    ctype const c6 =
        B.vx * M.c1 - B.vy * M.c0 - B.vz * M.c7 - B.mx * M.c5 + B.my * M.c4 + B.mz * M.c3;
    ctype const c7 = -B.mx * M.c0 - B.my * M.c1 - B.mz * M.c2;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& M,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y + M.c2 * v.z;
    ctype const c1 = -M.c0 * v.w + M.c3 * v.x - M.c5 * v.z + M.c6 * v.y;
    ctype const c2 = -M.c1 * v.w + M.c3 * v.y + M.c4 * v.z - M.c6 * v.x;
    ctype const c3 = -M.c2 * v.w + M.c3 * v.z - M.c4 * v.y + M.c5 * v.x;
    ctype const c4 = M.c1 * v.z - M.c2 * v.y - M.c7 * v.x;
    ctype const c5 = -M.c0 * v.z + M.c2 * v.x - M.c7 * v.y;
    ctype const c6 = M.c0 * v.y - M.c1 * v.x - M.c7 * v.z;
    ctype const c7 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z - M.c7 * v.w;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 + v.z * M.c2;
    ctype const c1 = -v.x * M.c3 + v.y * M.c6 - v.z * M.c5 + v.w * M.c0;
    ctype const c2 = -v.x * M.c6 - v.y * M.c3 + v.z * M.c4 + v.w * M.c1;
    ctype const c3 = v.x * M.c5 - v.y * M.c4 - v.z * M.c3 + v.w * M.c2;
    ctype const c4 = -v.x * M.c7 + v.y * M.c2 - v.z * M.c1;
    ctype const c5 = -v.x * M.c2 - v.y * M.c7 + v.z * M.c0;
    ctype const c6 = v.x * M.c1 - v.y * M.c0 - v.z * M.c7;
    ctype const c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6 + v.w * M.c7;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& M,
                                                        Scalar3dp<U> s)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                        MVec3dp_U<U> const& M)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator*([[maybe_unused]] PScalar3dp<T>,
                                                        [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// pga3dp gpr :: gpr(ps,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                     TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = -ctype(ps) * t.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp gpr :: gpr(trivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                     PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = t.w * ctype(ps);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp gpr :: gpr(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                       BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.mx;
    ctype const c1 = ctype(ps) * B.my;
    ctype const c2 = ctype(ps) * B.mz;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp gpr :: gpr(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                       PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.mx * ctype(ps);
    ctype const c1 = B.my * ctype(ps);
    ctype const c2 = B.mz * ctype(ps);
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp gpr :: gpr(ps,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.x;
    ctype const c1 = -ctype(ps) * v.y;
    ctype const c2 = -ctype(ps) * v.z;
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp gpr :: gpr(vec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = 0.0;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp gpr :: gpr(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * ctype(s));
}

// pga3dp gpr :: gpr(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) * ctype(ps));
}

// t1 * t2 = -dot(t1, t2) + cmt(t1, t2)
// pga3dp gpr :: gpr(trivec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t1,
                                                        TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t1.w * t2.w;
    ctype const c1 = t1.x * t2.w - t1.w * t2.x;
    ctype const c2 = t1.y * t2.w - t1.w * t2.y;
    ctype const c3 = t1.z * t2.w - t1.w * t2.z;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// t * B = -(t >> B) + cmt(t, B)
// pga3dp gpr :: gpr(trivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.w * B.mx;
    ctype const c1 = t.w * B.my;
    ctype const c2 = t.w * B.mz;
    ctype const c3 = -t.x * B.mx - t.y * B.my - t.z * B.mz;
    ctype const c4 = -t.y * B.mz + t.z * B.my + t.w * B.vx;
    ctype const c5 = t.x * B.mz - t.z * B.mx + t.w * B.vy;
    ctype const c6 = -t.x * B.my + t.y * B.mx + t.w * B.vz;
    ctype const c7 = 0.0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// B * t = -(B << t) + cmt(B, t)
// pga3dp gpr :: gpr(bivec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.mx * t.w;
    ctype const c1 = B.my * t.w;
    ctype const c2 = B.mz * t.w;
    ctype const c3 = -B.mx * t.x - B.my * t.y - B.mz * t.z;
    ctype const c4 = -B.vx * t.w - B.my * t.z + B.mz * t.y;
    ctype const c5 = -B.vy * t.w + B.mx * t.z - B.mz * t.x;
    ctype const c6 = -B.vz * t.w - B.mx * t.y + B.my * t.x;
    ctype const c7 = 0.0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// t * v = (t >> v) + wdg(t, v)
// pga3dp gpr :: gpr(trivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -t.y * v.z + t.z * v.y;
    ctype const c2 = t.x * v.z - t.z * v.x;
    ctype const c3 = -t.x * v.y + t.y * v.x;
    ctype const c4 = -t.w * v.x;
    ctype const c5 = -t.w * v.y;
    ctype const c6 = -t.w * v.z;
    ctype const c7 = -t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// v * t = (v << t) + wdg(v,t)
// pga3dp gpr :: gpr(vec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.y * t.z - v.z * t.y;
    ctype const c2 = -v.x * t.z + v.z * t.x;
    ctype const c3 = v.x * t.y - v.y * t.x;
    ctype const c4 = -v.x * t.w;
    ctype const c5 = -v.y * t.w;
    ctype const c6 = -v.z * t.w;
    ctype const c7 = v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(trivec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(s);
    ctype const c1 = t.y * ctype(s);
    ctype const c2 = t.z * ctype(s);
    ctype const c3 = t.w * ctype(s);
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp gpr :: gpr(s,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = ctype(s) * t.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// A * B = -dot(A,B) + cmt(A,B) + wdg(A,B)  (in 4d and higher dimensional spaces)
// The bivector product AxB = cmt(A,B) = 0.5*(A*B-B*A) is the commutator product.
//
// pga3dp gpr :: gpr(bivec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B1,
                                                        BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.mx * B2.mx - B1.my * B2.my - B1.mz * B2.mz;
    ctype const c1 = -B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz + B1.mz * B2.vy;
    ctype const c2 = B1.vx * B2.mz - B1.vz * B2.mx + B1.mx * B2.vz - B1.mz * B2.vx;
    ctype const c3 = -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx;
    ctype const c4 = -B1.my * B2.mz + B1.mz * B2.my;
    ctype const c5 = B1.mx * B2.mz - B1.mz * B2.mx;
    ctype const c6 = -B1.mx * B2.my + B1.my * B2.mx;
    ctype const c7 = -B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz - B1.mx * B2.vx -
                     B1.my * B2.vy - B1.mz * B2.vz;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// B * v = (v << B) + wdg(B,v)
// pga3dp gpr :: gpr(bivec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.my * v.z + B.mz * v.y;
    ctype const c1 = B.mx * v.z - B.mz * v.x;
    ctype const c2 = -B.mx * v.y + B.my * v.x;
    ctype const c3 = B.vx * v.x + B.vy * v.y + B.vz * v.z;
    ctype const c4 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c5 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c6 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c7 = -B.mx * v.x - B.my * v.y - B.mz * v.z;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// v * B = (B >> v) + wdg(v,B)
// pga3dp gpr :: gpr(vec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.mz + v.z * B.my;
    ctype const c1 = v.x * B.mz - v.z * B.mx;
    ctype const c2 = -v.x * B.my + v.y * B.mx;
    ctype const c3 = -v.x * B.vx - v.y * B.vy - v.z * B.vz;
    ctype const c4 = -v.y * B.vz + v.z * B.vy + v.w * B.mx;
    ctype const c5 = v.x * B.vz - v.z * B.vx + v.w * B.my;
    ctype const c6 = -v.x * B.vy + v.y * B.vx + v.w * B.mz;
    ctype const c7 = -v.x * B.mx - v.y * B.my - v.z * B.mz;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp gpr :: gpr(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                       Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(s);
    ctype const c1 = B.vy * ctype(s);
    ctype const c2 = B.vz * ctype(s);
    ctype const c3 = B.mx * ctype(s);
    ctype const c4 = B.my * ctype(s);
    ctype const c5 = B.mz * ctype(s);
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp gpr :: gpr(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                       BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.vx;
    ctype const c1 = ctype(s) * B.vy;
    ctype const c2 = ctype(s) * B.vz;
    ctype const c3 = ctype(s) * B.mx;
    ctype const c4 = ctype(s) * B.my;
    ctype const c5 = ctype(s) * B.mz;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp gpr :: gpr(vec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v1,
                                                        Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    ctype const c1 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c2 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c3 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c4 = v1.y * v2.z - v1.z * v2.y;
    ctype const c5 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c6 = v1.x * v2.y - v1.y * v2.x;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp gpr :: gpr(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    ctype const c3 = v.w * ctype(s);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp gpr :: gpr(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    ctype const c3 = ctype(s) * v.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp gpr :: gpr(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
////////////////////////////////////////////////////////////////////////////////

// pga3dp rgpr :: rgpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 +
                     A.c4 * B.c14 - A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 -
                     A.c8 * B.c5 - A.c9 * B.c6 - A.c10 * B.c7 - A.c11 * B.c1 -
                     A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 + A.c15 * B.c0;
    ctype const c1 = A.c0 * B.c11 + A.c1 * B.c15 + A.c2 * B.c7 - A.c3 * B.c6 -
                     A.c4 * B.c8 + A.c5 * B.c14 + A.c6 * B.c3 - A.c7 * B.c2 +
                     A.c8 * B.c4 + A.c9 * B.c13 - A.c10 * B.c12 - A.c11 * B.c0 -
                     A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 + A.c15 * B.c1;
    ctype const c2 = A.c0 * B.c12 - A.c1 * B.c7 + A.c2 * B.c15 + A.c3 * B.c5 -
                     A.c4 * B.c9 - A.c5 * B.c3 + A.c6 * B.c14 + A.c7 * B.c1 -
                     A.c8 * B.c13 + A.c9 * B.c4 + A.c10 * B.c11 + A.c11 * B.c10 -
                     A.c12 * B.c0 - A.c13 * B.c8 + A.c14 * B.c6 + A.c15 * B.c2;
    ctype const c3 = A.c0 * B.c13 + A.c1 * B.c6 - A.c2 * B.c5 + A.c3 * B.c15 -
                     A.c4 * B.c10 + A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c14 +
                     A.c8 * B.c12 - A.c9 * B.c11 + A.c10 * B.c4 - A.c11 * B.c9 +
                     A.c12 * B.c8 - A.c13 * B.c0 + A.c14 * B.c7 + A.c15 * B.c3;
    ctype const c4 = A.c4 * B.c15 - A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 -
                     A.c11 * B.c5 - A.c12 * B.c6 - A.c13 * B.c7 + A.c15 * B.c4;
    ctype const c5 = -A.c4 * B.c11 + A.c5 * B.c15 + A.c6 * B.c7 - A.c7 * B.c6 -
                     A.c11 * B.c4 - A.c12 * B.c13 + A.c13 * B.c12 + A.c15 * B.c5;
    ctype const c6 = -A.c4 * B.c12 - A.c5 * B.c7 + A.c6 * B.c15 + A.c7 * B.c5 +
                     A.c11 * B.c13 - A.c12 * B.c4 - A.c13 * B.c11 + A.c15 * B.c6;
    ctype const c7 = -A.c4 * B.c13 + A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c15 -
                     A.c11 * B.c12 + A.c12 * B.c11 - A.c13 * B.c4 + A.c15 * B.c7;
    ctype const c8 = A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c13 + A.c3 * B.c12 +
                     A.c4 * B.c1 + A.c5 * B.c0 + A.c6 * B.c10 - A.c7 * B.c9 +
                     A.c8 * B.c15 + A.c9 * B.c7 - A.c10 * B.c6 + A.c11 * B.c14 +
                     A.c12 * B.c3 - A.c13 * B.c2 - A.c14 * B.c11 + A.c15 * B.c8;
    ctype const c9 = A.c0 * B.c6 + A.c1 * B.c13 - A.c2 * B.c4 - A.c3 * B.c11 +
                     A.c4 * B.c2 - A.c5 * B.c10 + A.c6 * B.c0 + A.c7 * B.c8 -
                     A.c8 * B.c7 + A.c9 * B.c15 + A.c10 * B.c5 - A.c11 * B.c3 +
                     A.c12 * B.c14 + A.c13 * B.c1 - A.c14 * B.c12 + A.c15 * B.c9;
    ctype const c10 = A.c0 * B.c7 - A.c1 * B.c12 + A.c2 * B.c11 - A.c3 * B.c4 +
                      A.c4 * B.c3 + A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c0 +
                      A.c8 * B.c6 - A.c9 * B.c5 + A.c10 * B.c15 + A.c11 * B.c2 -
                      A.c12 * B.c1 + A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c10;
    ctype const c11 = A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c13 - A.c7 * B.c12 +
                      A.c11 * B.c15 + A.c12 * B.c7 - A.c13 * B.c6 + A.c15 * B.c11;
    ctype const c12 = A.c4 * B.c6 - A.c5 * B.c13 + A.c6 * B.c4 + A.c7 * B.c11 -
                      A.c11 * B.c7 + A.c12 * B.c15 + A.c13 * B.c5 + A.c15 * B.c12;
    ctype const c13 = A.c4 * B.c7 + A.c5 * B.c12 - A.c6 * B.c11 + A.c7 * B.c4 +
                      A.c11 * B.c6 - A.c12 * B.c5 + A.c13 * B.c15 + A.c15 * B.c13;
    ctype const c14 = A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 -
                      A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 -
                      A.c8 * B.c11 - A.c9 * B.c12 - A.c10 * B.c13 + A.c11 * B.c8 +
                      A.c12 * B.c9 + A.c13 * B.c10 + A.c14 * B.c15 + A.c15 * B.c14;
    ctype const c15 = -A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7 +
                      A.c11 * B.c11 + A.c12 * B.c12 + A.c13 * B.c13 + A.c15 * B.c15;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& A,
                                                 MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 - A.c8 * B.c1 -
                     A.c9 * B.c2 - A.c10 * B.c3 + A.c15 * B.c0;
    ctype const c1 = A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 - A.c4 * B.c4 -
                     A.c11 * B.c0 - A.c12 * B.c6 + A.c13 * B.c5 + A.c14 * B.c1;
    ctype const c2 = -A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 - A.c4 * B.c5 +
                     A.c11 * B.c6 - A.c12 * B.c0 - A.c13 * B.c4 + A.c14 * B.c2;
    ctype const c3 = A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 - A.c4 * B.c6 -
                     A.c11 * B.c5 + A.c12 * B.c4 - A.c13 * B.c0 + A.c14 * B.c3;
    ctype const c4 = A.c4 * B.c7 - A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3;
    ctype const c5 = A.c5 * B.c7 + A.c6 * B.c3 - A.c7 * B.c2 + A.c15 * B.c1;
    ctype const c6 = -A.c5 * B.c3 + A.c6 * B.c7 + A.c7 * B.c1 + A.c15 * B.c2;
    ctype const c7 = A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c7 + A.c15 * B.c3;
    ctype const c8 = A.c0 * B.c1 + A.c5 * B.c0 + A.c6 * B.c6 - A.c7 * B.c5 + A.c8 * B.c7 +
                     A.c9 * B.c3 - A.c10 * B.c2 + A.c15 * B.c4;
    ctype const c9 = A.c0 * B.c2 - A.c5 * B.c6 + A.c6 * B.c0 + A.c7 * B.c4 - A.c8 * B.c3 +
                     A.c9 * B.c7 + A.c10 * B.c1 + A.c15 * B.c5;
    ctype const c10 = A.c0 * B.c3 + A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c0 +
                      A.c8 * B.c2 - A.c9 * B.c1 + A.c10 * B.c7 + A.c15 * B.c6;
    ctype const c11 = A.c4 * B.c1 + A.c11 * B.c7 + A.c12 * B.c3 - A.c13 * B.c2;
    ctype const c12 = A.c4 * B.c2 - A.c11 * B.c3 + A.c12 * B.c7 + A.c13 * B.c1;
    ctype const c13 = A.c4 * B.c3 + A.c11 * B.c2 - A.c12 * B.c1 + A.c13 * B.c7;
    ctype const c14 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 - A.c4 * B.c0 +
                      A.c11 * B.c4 + A.c12 * B.c5 + A.c13 * B.c6 + A.c14 * B.c7;
    ctype const c15 = -A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 + A.c15 * B.c7;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& A,
                                                 MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 -
                     A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7 + A.c7 * B.c0;
    ctype const c1 = A.c0 * B.c11 + A.c1 * B.c14 + A.c2 * B.c3 - A.c3 * B.c2 +
                     A.c4 * B.c4 + A.c5 * B.c13 - A.c6 * B.c12 + A.c7 * B.c1;
    ctype const c2 = A.c0 * B.c12 - A.c1 * B.c3 + A.c2 * B.c14 + A.c3 * B.c1 -
                     A.c4 * B.c13 + A.c5 * B.c4 + A.c6 * B.c11 + A.c7 * B.c2;
    ctype const c3 = A.c0 * B.c13 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c14 +
                     A.c4 * B.c12 - A.c5 * B.c11 + A.c6 * B.c4 + A.c7 * B.c3;
    ctype const c4 = -A.c1 * B.c11 - A.c2 * B.c12 - A.c3 * B.c13 + A.c7 * B.c4;
    ctype const c5 = A.c1 * B.c15 + A.c2 * B.c7 - A.c3 * B.c6 + A.c7 * B.c5;
    ctype const c6 = -A.c1 * B.c7 + A.c2 * B.c15 + A.c3 * B.c5 + A.c7 * B.c6;
    ctype const c7 = A.c1 * B.c6 - A.c2 * B.c5 + A.c3 * B.c15 + A.c7 * B.c7;
    ctype const c8 = A.c0 * B.c5 + A.c1 * B.c0 + A.c2 * B.c10 - A.c3 * B.c9 +
                     A.c4 * B.c15 + A.c5 * B.c7 - A.c6 * B.c6 + A.c7 * B.c8;
    ctype const c9 = A.c0 * B.c6 - A.c1 * B.c10 + A.c2 * B.c0 + A.c3 * B.c8 -
                     A.c4 * B.c7 + A.c5 * B.c15 + A.c6 * B.c5 + A.c7 * B.c9;
    ctype const c10 = A.c0 * B.c7 + A.c1 * B.c9 - A.c2 * B.c8 + A.c3 * B.c0 +
                      A.c4 * B.c6 - A.c5 * B.c5 + A.c6 * B.c15 + A.c7 * B.c10;
    ctype const c11 = A.c1 * B.c4 + A.c2 * B.c13 - A.c3 * B.c12 + A.c7 * B.c11;
    ctype const c12 = -A.c1 * B.c13 + A.c2 * B.c4 + A.c3 * B.c11 + A.c7 * B.c12;
    ctype const c13 = A.c1 * B.c12 - A.c2 * B.c11 + A.c3 * B.c4 + A.c7 * B.c13;
    ctype const c14 = A.c0 * B.c4 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 -
                      A.c4 * B.c11 - A.c5 * B.c12 - A.c6 * B.c13 + A.c7 * B.c14;
    ctype const c15 = -A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 + A.c7 * B.c15;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& A,
                                                 MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c7 -
                     A.c11 * B.c0 - A.c12 * B.c1 - A.c13 * B.c2 - A.c14 * B.c3;
    ctype const c1 = A.c0 * B.c4 + A.c5 * B.c7 + A.c6 * B.c2 - A.c7 * B.c1 + A.c8 * B.c3 +
                     A.c9 * B.c6 - A.c10 * B.c5 + A.c15 * B.c0;
    ctype const c2 = A.c0 * B.c5 - A.c5 * B.c2 + A.c6 * B.c7 + A.c7 * B.c0 - A.c8 * B.c6 +
                     A.c9 * B.c3 + A.c10 * B.c4 + A.c15 * B.c1;
    ctype const c3 = A.c0 * B.c6 + A.c5 * B.c1 - A.c6 * B.c0 + A.c7 * B.c7 + A.c8 * B.c5 -
                     A.c9 * B.c4 + A.c10 * B.c3 + A.c15 * B.c2;
    ctype const c4 = -A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 + A.c15 * B.c3;
    ctype const c5 = -A.c4 * B.c4 - A.c11 * B.c3 - A.c12 * B.c6 + A.c13 * B.c5;
    ctype const c6 = -A.c4 * B.c5 + A.c11 * B.c6 - A.c12 * B.c3 - A.c13 * B.c4;
    ctype const c7 = -A.c4 * B.c6 - A.c11 * B.c5 + A.c12 * B.c4 - A.c13 * B.c3;
    ctype const c8 = -A.c1 * B.c3 - A.c2 * B.c6 + A.c3 * B.c5 + A.c4 * B.c0 +
                     A.c11 * B.c7 + A.c12 * B.c2 - A.c13 * B.c1 - A.c14 * B.c4;
    ctype const c9 = A.c1 * B.c6 - A.c2 * B.c3 - A.c3 * B.c4 + A.c4 * B.c1 -
                     A.c11 * B.c2 + A.c12 * B.c7 + A.c13 * B.c0 - A.c14 * B.c5;
    ctype const c10 = -A.c1 * B.c5 + A.c2 * B.c4 - A.c3 * B.c3 + A.c4 * B.c2 +
                      A.c11 * B.c1 - A.c12 * B.c0 + A.c13 * B.c7 - A.c14 * B.c6;
    ctype const c11 = A.c5 * B.c3 + A.c6 * B.c6 - A.c7 * B.c5 + A.c15 * B.c4;
    ctype const c12 = -A.c5 * B.c6 + A.c6 * B.c3 + A.c7 * B.c4 + A.c15 * B.c5;
    ctype const c13 = A.c5 * B.c5 - A.c6 * B.c4 + A.c7 * B.c3 + A.c15 * B.c6;
    ctype const c14 = A.c0 * B.c3 - A.c5 * B.c0 - A.c6 * B.c1 - A.c7 * B.c2 -
                      A.c8 * B.c4 - A.c9 * B.c5 - A.c10 * B.c6 + A.c15 * B.c7;
    ctype const c15 = -A.c4 * B.c3 + A.c11 * B.c4 + A.c12 * B.c5 + A.c13 * B.c6;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& A,
                                                 MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 -
                     A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 - A.c7 * B.c4;
    ctype const c1 = A.c0 * B.c15 + A.c1 * B.c7 - A.c2 * B.c6 - A.c3 * B.c8 -
                     A.c4 * B.c0 - A.c5 * B.c10 + A.c6 * B.c9 + A.c7 * B.c5;
    ctype const c2 = -A.c0 * B.c7 + A.c1 * B.c15 + A.c2 * B.c5 - A.c3 * B.c9 +
                     A.c4 * B.c10 - A.c5 * B.c0 - A.c6 * B.c8 + A.c7 * B.c6;
    ctype const c3 = A.c0 * B.c6 - A.c1 * B.c5 + A.c2 * B.c15 - A.c3 * B.c10 -
                     A.c4 * B.c9 + A.c5 * B.c8 - A.c6 * B.c0 + A.c7 * B.c7;
    ctype const c4 = A.c3 * B.c15 - A.c4 * B.c5 - A.c5 * B.c6 - A.c6 * B.c7;
    ctype const c5 = -A.c3 * B.c11 - A.c4 * B.c4 - A.c5 * B.c13 + A.c6 * B.c12;
    ctype const c6 = -A.c3 * B.c12 + A.c4 * B.c13 - A.c5 * B.c4 - A.c6 * B.c11;
    ctype const c7 = -A.c3 * B.c13 - A.c4 * B.c12 + A.c5 * B.c11 - A.c6 * B.c4;
    ctype const c8 = -A.c0 * B.c4 - A.c1 * B.c13 + A.c2 * B.c12 + A.c3 * B.c1 +
                     A.c4 * B.c14 + A.c5 * B.c3 - A.c6 * B.c2 - A.c7 * B.c11;
    ctype const c9 = A.c0 * B.c13 - A.c1 * B.c4 - A.c2 * B.c11 + A.c3 * B.c2 -
                     A.c4 * B.c3 + A.c5 * B.c14 + A.c6 * B.c1 - A.c7 * B.c12;
    ctype const c10 = -A.c0 * B.c12 + A.c1 * B.c11 - A.c2 * B.c4 + A.c3 * B.c3 +
                      A.c4 * B.c2 - A.c5 * B.c1 + A.c6 * B.c14 - A.c7 * B.c13;
    ctype const c11 = A.c3 * B.c5 + A.c4 * B.c15 + A.c5 * B.c7 - A.c6 * B.c6;
    ctype const c12 = A.c3 * B.c6 - A.c4 * B.c7 + A.c5 * B.c15 + A.c6 * B.c5;
    ctype const c13 = A.c3 * B.c7 + A.c4 * B.c6 - A.c5 * B.c5 + A.c6 * B.c15;
    ctype const c14 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c7 - A.c3 * B.c0 +
                      A.c4 * B.c8 + A.c5 * B.c9 + A.c6 * B.c10 + A.c7 * B.c15;
    ctype const c15 = -A.c3 * B.c4 + A.c4 * B.c11 + A.c5 * B.c12 + A.c6 * B.c13;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& M, PScalar3dp<U> ps)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps, MVec3dp<U> const& M)
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
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv,trivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& M,
                                                 TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * t.x + M.c2 * t.y + M.c3 * t.z + M.c4 * t.w;
    ctype const c1 = M.c0 * t.x + M.c5 * t.w + M.c9 * t.z - M.c10 * t.y;
    ctype const c2 = M.c0 * t.y + M.c6 * t.w - M.c8 * t.z + M.c10 * t.x;
    ctype const c3 = M.c0 * t.z + M.c7 * t.w + M.c8 * t.y - M.c9 * t.x;
    ctype const c4 = -M.c5 * t.x - M.c6 * t.y - M.c7 * t.z;
    ctype const c5 = -M.c4 * t.x - M.c12 * t.z + M.c13 * t.y;
    ctype const c6 = -M.c4 * t.y + M.c11 * t.z - M.c13 * t.x;
    ctype const c7 = -M.c4 * t.z - M.c11 * t.y + M.c12 * t.x;
    ctype const c8 = -M.c2 * t.z + M.c3 * t.y + M.c11 * t.w - M.c14 * t.x;
    ctype const c9 = M.c1 * t.z - M.c3 * t.x + M.c12 * t.w - M.c14 * t.y;
    ctype const c10 = -M.c1 * t.y + M.c2 * t.x + M.c13 * t.w - M.c14 * t.z;
    ctype const c11 = M.c6 * t.z - M.c7 * t.y + M.c15 * t.x;
    ctype const c12 = -M.c5 * t.z + M.c7 * t.x + M.c15 * t.y;
    ctype const c13 = M.c5 * t.y - M.c6 * t.x + M.c15 * t.z;
    ctype const c14 = -M.c8 * t.x - M.c9 * t.y - M.c10 * t.z + M.c15 * t.w;
    ctype const c15 = M.c11 * t.x + M.c12 * t.y + M.c13 * t.z;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(trivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                 MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3 - t.w * M.c4;
    ctype const c1 = -t.x * M.c0 - t.y * M.c10 + t.z * M.c9 + t.w * M.c5;
    ctype const c2 = t.x * M.c10 - t.y * M.c0 - t.z * M.c8 + t.w * M.c6;
    ctype const c3 = -t.x * M.c9 + t.y * M.c8 - t.z * M.c0 + t.w * M.c7;
    ctype const c4 = -t.x * M.c5 - t.y * M.c6 - t.z * M.c7;
    ctype const c5 = -t.x * M.c4 - t.y * M.c13 + t.z * M.c12;
    ctype const c6 = t.x * M.c13 - t.y * M.c4 - t.z * M.c11;
    ctype const c7 = -t.x * M.c12 + t.y * M.c11 - t.z * M.c4;
    ctype const c8 = t.x * M.c14 + t.y * M.c3 - t.z * M.c2 - t.w * M.c11;
    ctype const c9 = -t.x * M.c3 + t.y * M.c14 + t.z * M.c1 - t.w * M.c12;
    ctype const c10 = t.x * M.c2 - t.y * M.c1 + t.z * M.c14 - t.w * M.c13;
    ctype const c11 = t.x * M.c15 + t.y * M.c7 - t.z * M.c6;
    ctype const c12 = -t.x * M.c7 + t.y * M.c15 + t.z * M.c5;
    ctype const c13 = t.x * M.c6 - t.y * M.c5 + t.z * M.c15;
    ctype const c14 = t.x * M.c8 + t.y * M.c9 + t.z * M.c10 + t.w * M.c15;
    ctype const c15 = t.x * M.c11 + t.y * M.c12 + t.z * M.c13;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& M,
                                                 BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c5 * B.mx - M.c6 * B.my - M.c7 * B.mz - M.c8 * B.vx -
                     M.c9 * B.vy - M.c10 * B.vz;
    ctype const c1 = M.c2 * B.vz - M.c3 * B.vy - M.c4 * B.mx - M.c12 * B.mz +
                     M.c13 * B.my + M.c14 * B.vx;
    ctype const c2 = -M.c1 * B.vz + M.c3 * B.vx - M.c4 * B.my + M.c11 * B.mz -
                     M.c13 * B.mx + M.c14 * B.vy;
    ctype const c3 = M.c1 * B.vy - M.c2 * B.vx - M.c4 * B.mz - M.c11 * B.my +
                     M.c12 * B.mx + M.c14 * B.vz;
    ctype const c4 = -M.c11 * B.vx - M.c12 * B.vy - M.c13 * B.vz;
    ctype const c5 = M.c6 * B.vz - M.c7 * B.vy + M.c15 * B.vx;
    ctype const c6 = -M.c5 * B.vz + M.c7 * B.vx + M.c15 * B.vy;
    ctype const c7 = M.c5 * B.vy - M.c6 * B.vx + M.c15 * B.vz;
    ctype const c8 = M.c0 * B.vx + M.c6 * B.mz - M.c7 * B.my + M.c9 * B.vz -
                     M.c10 * B.vy + M.c15 * B.mx;
    ctype const c9 = M.c0 * B.vy - M.c5 * B.mz + M.c7 * B.mx - M.c8 * B.vz +
                     M.c10 * B.vx + M.c15 * B.my;
    ctype const c10 = M.c0 * B.vz + M.c5 * B.my - M.c6 * B.mx + M.c8 * B.vy -
                      M.c9 * B.vx + M.c15 * B.mz;
    ctype const c11 = M.c4 * B.vx + M.c12 * B.vz - M.c13 * B.vy;
    ctype const c12 = M.c4 * B.vy - M.c11 * B.vz + M.c13 * B.vx;
    ctype const c13 = M.c4 * B.vz + M.c11 * B.vy - M.c12 * B.vx;
    ctype const c14 = -M.c1 * B.vx - M.c2 * B.vy - M.c3 * B.vz + M.c11 * B.mx +
                      M.c12 * B.my + M.c13 * B.mz;
    ctype const c15 = -M.c5 * B.vx - M.c6 * B.vy - M.c7 * B.vz;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B,
                                                 MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c8 - B.vy * M.c9 - B.vz * M.c10 - B.mx * M.c5 -
                     B.my * M.c6 - B.mz * M.c7;
    ctype const c1 = B.vx * M.c14 + B.vy * M.c3 - B.vz * M.c2 + B.mx * M.c4 +
                     B.my * M.c13 - B.mz * M.c12;
    ctype const c2 = -B.vx * M.c3 + B.vy * M.c14 + B.vz * M.c1 - B.mx * M.c13 +
                     B.my * M.c4 + B.mz * M.c11;
    ctype const c3 = B.vx * M.c2 - B.vy * M.c1 + B.vz * M.c14 + B.mx * M.c12 -
                     B.my * M.c11 + B.mz * M.c4;
    ctype const c4 = -B.vx * M.c11 - B.vy * M.c12 - B.vz * M.c13;
    ctype const c5 = B.vx * M.c15 + B.vy * M.c7 - B.vz * M.c6;
    ctype const c6 = -B.vx * M.c7 + B.vy * M.c15 + B.vz * M.c5;
    ctype const c7 = B.vx * M.c6 - B.vy * M.c5 + B.vz * M.c15;
    ctype const c8 = B.vx * M.c0 + B.vy * M.c10 - B.vz * M.c9 + B.mx * M.c15 +
                     B.my * M.c7 - B.mz * M.c6;
    ctype const c9 = -B.vx * M.c10 + B.vy * M.c0 + B.vz * M.c8 - B.mx * M.c7 +
                     B.my * M.c15 + B.mz * M.c5;
    ctype const c10 = B.vx * M.c9 - B.vy * M.c8 + B.vz * M.c0 + B.mx * M.c6 -
                      B.my * M.c5 + B.mz * M.c15;
    ctype const c11 = B.vx * M.c4 + B.vy * M.c13 - B.vz * M.c12;
    ctype const c12 = -B.vx * M.c13 + B.vy * M.c4 + B.vz * M.c11;
    ctype const c13 = B.vx * M.c12 - B.vy * M.c11 + B.vz * M.c4;
    ctype const c14 = -B.vx * M.c1 - B.vy * M.c2 - B.vz * M.c3 - B.mx * M.c11 -
                      B.my * M.c12 - B.mz * M.c13;
    ctype const c15 = -B.vx * M.c5 - B.vy * M.c6 - B.vz * M.c7;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& M, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c11 * v.x - M.c12 * v.y - M.c13 * v.z - M.c14 * v.w;
    ctype const c1 = M.c6 * v.z - M.c7 * v.y + M.c8 * v.w + M.c15 * v.x;
    ctype const c2 = -M.c5 * v.z + M.c7 * v.x + M.c9 * v.w + M.c15 * v.y;
    ctype const c3 = M.c5 * v.y - M.c6 * v.x + M.c10 * v.w + M.c15 * v.z;
    ctype const c4 = M.c15 * v.w;
    ctype const c5 = -M.c11 * v.w;
    ctype const c6 = -M.c12 * v.w;
    ctype const c7 = -M.c13 * v.w;
    ctype const c8 = -M.c1 * v.w + M.c4 * v.x + M.c12 * v.z - M.c13 * v.y;
    ctype const c9 = -M.c2 * v.w + M.c4 * v.y - M.c11 * v.z + M.c13 * v.x;
    ctype const c10 = -M.c3 * v.w + M.c4 * v.z + M.c11 * v.y - M.c12 * v.x;
    ctype const c11 = M.c5 * v.w;
    ctype const c12 = M.c6 * v.w;
    ctype const c13 = M.c7 * v.w;
    ctype const c14 = M.c0 * v.w - M.c5 * v.x - M.c6 * v.y - M.c7 * v.z;
    ctype const c15 = -M.c4 * v.w;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v, MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c11 + v.y * M.c12 + v.z * M.c13 + v.w * M.c14;
    ctype const c1 = v.x * M.c15 + v.y * M.c7 - v.z * M.c6 - v.w * M.c8;
    ctype const c2 = -v.x * M.c7 + v.y * M.c15 + v.z * M.c5 - v.w * M.c9;
    ctype const c3 = v.x * M.c6 - v.y * M.c5 + v.z * M.c15 - v.w * M.c10;
    ctype const c4 = v.w * M.c15;
    ctype const c5 = -v.w * M.c11;
    ctype const c6 = -v.w * M.c12;
    ctype const c7 = -v.w * M.c13;
    ctype const c8 = -v.x * M.c4 - v.y * M.c13 + v.z * M.c12 + v.w * M.c1;
    ctype const c9 = v.x * M.c13 - v.y * M.c4 - v.z * M.c11 + v.w * M.c2;
    ctype const c10 = -v.x * M.c12 + v.y * M.c11 - v.z * M.c4 + v.w * M.c3;
    ctype const c11 = v.w * M.c5;
    ctype const c12 = v.w * M.c6;
    ctype const c13 = v.w * M.c7;
    ctype const c14 = -v.x * M.c5 - v.y * M.c6 - v.z * M.c7 - v.w * M.c0;
    ctype const c15 = -v.w * M.c4;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& M, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c15 * ctype(s);
    ctype const c1 = -M.c11 * ctype(s);
    ctype const c2 = -M.c12 * ctype(s);
    ctype const c3 = -M.c13 * ctype(s);
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = M.c5 * ctype(s);
    ctype const c9 = M.c6 * ctype(s);
    ctype const c10 = M.c7 * ctype(s);
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = -M.c4 * ctype(s);
    ctype const c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(Scalar3dp<T> s, MVec3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c15;
    ctype const c1 = ctype(s) * M.c11;
    ctype const c2 = ctype(s) * M.c12;
    ctype const c3 = ctype(s) * M.c13;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    ctype const c8 = ctype(s) * M.c5;
    ctype const c9 = ctype(s) * M.c6;
    ctype const c10 = ctype(s) * M.c7;
    ctype const c11 = 0.0;
    ctype const c12 = 0.0;
    ctype const c13 = 0.0;
    ctype const c14 = ctype(s) * M.c4;
    ctype const c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

// pga3dp rgpr :: rgpr(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& A,
                                                   MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
                     A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c7 * B.c1;
    ctype const c2 = -A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c7 * B.c2;
    ctype const c3 = A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 + A.c7 * B.c3;
    ctype const c4 = A.c0 * B.c1 + A.c1 * B.c0 + A.c2 * B.c6 - A.c3 * B.c5 + A.c4 * B.c7 +
                     A.c5 * B.c3 - A.c6 * B.c2 + A.c7 * B.c4;
    ctype const c5 = A.c0 * B.c2 - A.c1 * B.c6 + A.c2 * B.c0 + A.c3 * B.c4 - A.c4 * B.c3 +
                     A.c5 * B.c7 + A.c6 * B.c1 + A.c7 * B.c5;
    ctype const c6 = A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c0 + A.c4 * B.c2 -
                     A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c7 * B.c7;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& A,
                                                   MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c2 - A.c3 * B.c1 + A.c4 * B.c3 +
                     A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c0;
    ctype const c1 = A.c0 * B.c5 - A.c1 * B.c2 + A.c2 * B.c7 + A.c3 * B.c0 - A.c4 * B.c6 +
                     A.c5 * B.c3 + A.c6 * B.c4 + A.c7 * B.c1;
    ctype const c2 = A.c0 * B.c6 + A.c1 * B.c1 - A.c2 * B.c0 + A.c3 * B.c7 + A.c4 * B.c5 -
                     A.c5 * B.c4 + A.c6 * B.c3 + A.c7 * B.c2;
    ctype const c3 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 + A.c7 * B.c3;
    ctype const c4 = A.c1 * B.c3 + A.c2 * B.c6 - A.c3 * B.c5 + A.c7 * B.c4;
    ctype const c5 = -A.c1 * B.c6 + A.c2 * B.c3 + A.c3 * B.c4 + A.c7 * B.c5;
    ctype const c6 = A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c3 + A.c7 * B.c6;
    ctype const c7 = A.c0 * B.c3 - A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2 - A.c4 * B.c4 -
                     A.c5 * B.c5 - A.c6 * B.c6 + A.c7 * B.c7;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& A,
                                                   MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 - A.c3 * B.c4 - A.c4 * B.c0 -
                     A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype const c1 = -A.c0 * B.c3 + A.c1 * B.c7 + A.c2 * B.c1 - A.c3 * B.c5 +
                     A.c4 * B.c6 - A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c7 - A.c3 * B.c6 - A.c4 * B.c5 +
                     A.c5 * B.c4 - A.c6 * B.c0 + A.c7 * B.c3;
    ctype const c3 = A.c3 * B.c7 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype const c4 = A.c3 * B.c1 + A.c4 * B.c7 + A.c5 * B.c3 - A.c6 * B.c2;
    ctype const c5 = A.c3 * B.c2 - A.c4 * B.c3 + A.c5 * B.c7 + A.c6 * B.c1;
    ctype const c6 = A.c3 * B.c3 + A.c4 * B.c2 - A.c5 * B.c1 + A.c6 * B.c7;
    ctype const c7 = -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 - A.c3 * B.c0 +
                     A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& A,
                                                   PScalar3dp<U> ps)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps,
                                                   MVec3dp_E<U> const& B)
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
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(mv_e,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& M,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * t.x + M.c1 * t.w + M.c5 * t.z - M.c6 * t.y;
    ctype const c1 = M.c0 * t.y + M.c2 * t.w - M.c4 * t.z + M.c6 * t.x;
    ctype const c2 = M.c0 * t.z + M.c3 * t.w + M.c4 * t.y - M.c5 * t.x;
    ctype const c3 = -M.c1 * t.x - M.c2 * t.y - M.c3 * t.z;
    ctype const c4 = M.c2 * t.z - M.c3 * t.y + M.c7 * t.x;
    ctype const c5 = -M.c1 * t.z + M.c3 * t.x + M.c7 * t.y;
    ctype const c6 = M.c1 * t.y - M.c2 * t.x + M.c7 * t.z;
    ctype const c7 = -M.c4 * t.x - M.c5 * t.y - M.c6 * t.z + M.c7 * t.w;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(trivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                   MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c0 - t.y * M.c6 + t.z * M.c5 + t.w * M.c1;
    ctype const c1 = t.x * M.c6 - t.y * M.c0 - t.z * M.c4 + t.w * M.c2;
    ctype const c2 = -t.x * M.c5 + t.y * M.c4 - t.z * M.c0 + t.w * M.c3;
    ctype const c3 = -t.x * M.c1 - t.y * M.c2 - t.z * M.c3;
    ctype const c4 = t.x * M.c7 + t.y * M.c3 - t.z * M.c2;
    ctype const c5 = -t.x * M.c3 + t.y * M.c7 + t.z * M.c1;
    ctype const c6 = t.x * M.c2 - t.y * M.c1 + t.z * M.c7;
    ctype const c7 = t.x * M.c4 + t.y * M.c5 + t.z * M.c6 + t.w * M.c7;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(mv_e,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& M,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz - M.c4 * B.vx -
                     M.c5 * B.vy - M.c6 * B.vz;
    ctype const c1 = M.c2 * B.vz - M.c3 * B.vy + M.c7 * B.vx;
    ctype const c2 = -M.c1 * B.vz + M.c3 * B.vx + M.c7 * B.vy;
    ctype const c3 = M.c1 * B.vy - M.c2 * B.vx + M.c7 * B.vz;
    ctype const c4 =
        M.c0 * B.vx + M.c2 * B.mz - M.c3 * B.my + M.c5 * B.vz - M.c6 * B.vy + M.c7 * B.mx;
    ctype const c5 =
        M.c0 * B.vy - M.c1 * B.mz + M.c3 * B.mx - M.c4 * B.vz + M.c6 * B.vx + M.c7 * B.my;
    ctype const c6 =
        M.c0 * B.vz + M.c1 * B.my - M.c2 * B.mx + M.c4 * B.vy - M.c5 * B.vx + M.c7 * B.mz;
    ctype const c7 = -M.c1 * B.vx - M.c2 * B.vy - M.c3 * B.vz;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(bivec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B,
                                                   MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6 - B.mx * M.c1 -
                     B.my * M.c2 - B.mz * M.c3;
    ctype const c1 = B.vx * M.c7 + B.vy * M.c3 - B.vz * M.c2;
    ctype const c2 = -B.vx * M.c3 + B.vy * M.c7 + B.vz * M.c1;
    ctype const c3 = B.vx * M.c2 - B.vy * M.c1 + B.vz * M.c7;
    ctype const c4 =
        B.vx * M.c0 + B.vy * M.c6 - B.vz * M.c5 + B.mx * M.c7 + B.my * M.c3 - B.mz * M.c2;
    ctype const c5 = -B.vx * M.c6 + B.vy * M.c0 + B.vz * M.c4 - B.mx * M.c3 +
                     B.my * M.c7 + B.mz * M.c1;
    ctype const c6 =
        B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c0 + B.mx * M.c2 - B.my * M.c1 + B.mz * M.c7;
    ctype const c7 = -B.vx * M.c1 - B.vy * M.c2 - B.vz * M.c3;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& M,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * v.z - M.c3 * v.y + M.c4 * v.w + M.c7 * v.x;
    ctype const c1 = -M.c1 * v.z + M.c3 * v.x + M.c5 * v.w + M.c7 * v.y;
    ctype const c2 = M.c1 * v.y - M.c2 * v.x + M.c6 * v.w + M.c7 * v.z;
    ctype const c3 = M.c7 * v.w;
    ctype const c4 = M.c1 * v.w;
    ctype const c5 = M.c2 * v.w;
    ctype const c6 = M.c3 * v.w;
    ctype const c7 = M.c0 * v.w - M.c1 * v.x - M.c2 * v.y - M.c3 * v.z;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v,
                                                   MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c7 + v.y * M.c3 - v.z * M.c2 - v.w * M.c4;
    ctype const c1 = -v.x * M.c3 + v.y * M.c7 + v.z * M.c1 - v.w * M.c5;
    ctype const c2 = v.x * M.c2 - v.y * M.c1 + v.z * M.c7 - v.w * M.c6;
    ctype const c3 = v.w * M.c7;
    ctype const c4 = v.w * M.c1;
    ctype const c5 = v.w * M.c2;
    ctype const c6 = v.w * M.c3;
    ctype const c7 = -v.x * M.c1 - v.y * M.c2 - v.z * M.c3 - v.w * M.c0;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& M, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c7 * ctype(s);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = M.c1 * ctype(s);
    ctype const c5 = M.c2 * ctype(s);
    ctype const c6 = M.c3 * ctype(s);
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(Scalar3dp<T> s, MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c7;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = ctype(s) * M.c1;
    ctype const c5 = ctype(s) * M.c2;
    ctype const c6 = ctype(s) * M.c3;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(mv_u,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& A,
                                                   MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 - A.c4 * B.c0 -
                     A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3;
    ctype const c1 = -A.c3 * B.c4 - A.c4 * B.c3 - A.c5 * B.c6 + A.c6 * B.c5;
    ctype const c2 = -A.c3 * B.c5 + A.c4 * B.c6 - A.c5 * B.c3 - A.c6 * B.c4;
    ctype const c3 = -A.c3 * B.c6 - A.c4 * B.c5 + A.c5 * B.c4 - A.c6 * B.c3;
    ctype const c4 = -A.c0 * B.c3 - A.c1 * B.c6 + A.c2 * B.c5 + A.c3 * B.c0 +
                     A.c4 * B.c7 + A.c5 * B.c2 - A.c6 * B.c1 - A.c7 * B.c4;
    ctype const c5 = A.c0 * B.c6 - A.c1 * B.c3 - A.c2 * B.c4 + A.c3 * B.c1 - A.c4 * B.c2 +
                     A.c5 * B.c7 + A.c6 * B.c0 - A.c7 * B.c5;
    ctype const c6 = -A.c0 * B.c5 + A.c1 * B.c4 - A.c2 * B.c3 + A.c3 * B.c2 +
                     A.c4 * B.c1 - A.c5 * B.c0 + A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c7 = -A.c3 * B.c3 + A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& A,
                                                   PScalar3dp<U> ps)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps,
                                                   MVec3dp_U<U> const& B)
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
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(mv_u,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& M,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * t.x + M.c1 * t.y + M.c2 * t.z + M.c3 * t.w;
    ctype const c1 = -M.c3 * t.x - M.c5 * t.z + M.c6 * t.y;
    ctype const c2 = -M.c3 * t.y + M.c4 * t.z - M.c6 * t.x;
    ctype const c3 = -M.c3 * t.z - M.c4 * t.y + M.c5 * t.x;
    ctype const c4 = -M.c1 * t.z + M.c2 * t.y + M.c4 * t.w - M.c7 * t.x;
    ctype const c5 = M.c0 * t.z - M.c2 * t.x + M.c5 * t.w - M.c7 * t.y;
    ctype const c6 = -M.c0 * t.y + M.c1 * t.x + M.c6 * t.w - M.c7 * t.z;
    ctype const c7 = M.c4 * t.x + M.c5 * t.y + M.c6 * t.z;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(trivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                   MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * M.c0 - t.y * M.c1 - t.z * M.c2 - t.w * M.c3;
    ctype const c1 = -t.x * M.c3 - t.y * M.c6 + t.z * M.c5;
    ctype const c2 = t.x * M.c6 - t.y * M.c3 - t.z * M.c4;
    ctype const c3 = -t.x * M.c5 + t.y * M.c4 - t.z * M.c3;
    ctype const c4 = t.x * M.c7 + t.y * M.c2 - t.z * M.c1 - t.w * M.c4;
    ctype const c5 = -t.x * M.c2 + t.y * M.c7 + t.z * M.c0 - t.w * M.c5;
    ctype const c6 = t.x * M.c1 - t.y * M.c0 + t.z * M.c7 - t.w * M.c6;
    ctype const c7 = t.x * M.c4 + t.y * M.c5 + t.z * M.c6;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(mv_u,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& M,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        M.c1 * B.vz - M.c2 * B.vy - M.c3 * B.mx - M.c5 * B.mz + M.c6 * B.my + M.c7 * B.vx;
    ctype const c1 = -M.c0 * B.vz + M.c2 * B.vx - M.c3 * B.my + M.c4 * B.mz -
                     M.c6 * B.mx + M.c7 * B.vy;
    ctype const c2 =
        M.c0 * B.vy - M.c1 * B.vx - M.c3 * B.mz - M.c4 * B.my + M.c5 * B.mx + M.c7 * B.vz;
    ctype const c3 = -M.c4 * B.vx - M.c5 * B.vy - M.c6 * B.vz;
    ctype const c4 = M.c3 * B.vx + M.c5 * B.vz - M.c6 * B.vy;
    ctype const c5 = M.c3 * B.vy - M.c4 * B.vz + M.c6 * B.vx;
    ctype const c6 = M.c3 * B.vz + M.c4 * B.vy - M.c5 * B.vx;
    ctype const c7 = -M.c0 * B.vx - M.c1 * B.vy - M.c2 * B.vz + M.c4 * B.mx +
                     M.c5 * B.my + M.c6 * B.mz;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(bivec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B,
                                                   MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 =
        B.vx * M.c7 + B.vy * M.c2 - B.vz * M.c1 + B.mx * M.c3 + B.my * M.c6 - B.mz * M.c5;
    ctype const c1 = -B.vx * M.c2 + B.vy * M.c7 + B.vz * M.c0 - B.mx * M.c6 +
                     B.my * M.c3 + B.mz * M.c4;
    ctype const c2 =
        B.vx * M.c1 - B.vy * M.c0 + B.vz * M.c7 + B.mx * M.c5 - B.my * M.c4 + B.mz * M.c3;
    ctype const c3 = -B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6;
    ctype const c4 = B.vx * M.c3 + B.vy * M.c6 - B.vz * M.c5;
    ctype const c5 = -B.vx * M.c6 + B.vy * M.c3 + B.vz * M.c4;
    ctype const c6 = B.vx * M.c5 - B.vy * M.c4 + B.vz * M.c3;
    ctype const c7 = -B.vx * M.c0 - B.vy * M.c1 - B.vz * M.c2 - B.mx * M.c4 -
                     B.my * M.c5 - B.mz * M.c6;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(mv_u,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& M,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * v.x - M.c5 * v.y - M.c6 * v.z - M.c7 * v.w;
    ctype const c1 = -M.c4 * v.w;
    ctype const c2 = -M.c5 * v.w;
    ctype const c3 = -M.c6 * v.w;
    ctype const c4 = -M.c0 * v.w + M.c3 * v.x + M.c5 * v.z - M.c6 * v.y;
    ctype const c5 = -M.c1 * v.w + M.c3 * v.y - M.c4 * v.z + M.c6 * v.x;
    ctype const c6 = -M.c2 * v.w + M.c3 * v.z + M.c4 * v.y - M.c5 * v.x;
    ctype const c7 = -M.c3 * v.w;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(vec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v,
                                                   MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6 + v.w * M.c7;
    ctype const c1 = -v.w * M.c4;
    ctype const c2 = -v.w * M.c5;
    ctype const c3 = -v.w * M.c6;
    ctype const c4 = -v.x * M.c3 - v.y * M.c6 + v.z * M.c5 + v.w * M.c0;
    ctype const c5 = v.x * M.c6 - v.y * M.c3 - v.z * M.c4 + v.w * M.c1;
    ctype const c6 = -v.x * M.c5 + v.y * M.c4 - v.z * M.c3 + v.w * M.c2;
    ctype const c7 = -v.w * M.c3;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& M, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * ctype(s);
    ctype const c1 = -M.c5 * ctype(s);
    ctype const c2 = -M.c6 * ctype(s);
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = -M.c3 * ctype(s);
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(Scalar3dp<T> s, MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c4;
    ctype const c1 = ctype(s) * M.c5;
    ctype const c2 = ctype(s) * M.c6;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = ctype(s) * M.c3;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps1, PScalar3dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps1) * ctype(ps2));
}

// pga3dp rgpr :: rgpr(ps,trivec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * t.x;
    ctype const c1 = ctype(ps) * t.y;
    ctype const c2 = ctype(ps) * t.z;
    ctype const c3 = ctype(ps) * t.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(trivec,ps) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                   PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = t.x * ctype(ps);
    ctype const c1 = t.y * ctype(ps);
    ctype const c2 = t.z * ctype(ps);
    ctype const c3 = t.w * ctype(ps);
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.vx;
    ctype const c1 = ctype(ps) * B.vy;
    ctype const c2 = ctype(ps) * B.vz;
    ctype const c3 = ctype(ps) * B.mx;
    ctype const c4 = ctype(ps) * B.my;
    ctype const c5 = ctype(ps) * B.mz;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rgpr :: rgpr(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * ctype(ps);
    ctype const c1 = B.vy * ctype(ps);
    ctype const c2 = B.vz * ctype(ps);
    ctype const c3 = B.mx * ctype(ps);
    ctype const c4 = B.my * ctype(ps);
    ctype const c5 = B.mz * ctype(ps);
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rgpr :: rgpr(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    ctype const c3 = ctype(ps) * v.w;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    ctype const c3 = v.w * ctype(ps);
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rgpr(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(ps) * ctype(s));
}

// pga3dp rgpr :: rgpr(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rgpr(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s) * ctype(ps));
}

// pga3dp rgpr :: rgpr(trivec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t1,
                                                   TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -t1.y * t2.z + t1.z * t2.y;
    ctype const c2 = t1.x * t2.z - t1.z * t2.x;
    ctype const c3 = -t1.x * t2.y + t1.y * t2.x;
    ctype const c4 = t1.x * t2.w - t1.w * t2.x;
    ctype const c5 = t1.y * t2.w - t1.w * t2.y;
    ctype const c6 = t1.z * t2.w - t1.w * t2.z;
    ctype const c7 = t1.x * t2.x + t1.y * t2.y + t1.z * t2.z;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(trivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.y * B.mz + t.z * B.my + t.w * B.vx;
    ctype const c1 = t.x * B.mz - t.z * B.mx + t.w * B.vy;
    ctype const c2 = -t.x * B.my + t.y * B.mx + t.w * B.vz;
    ctype const c3 = -t.x * B.vx - t.y * B.vy - t.z * B.vz;
    ctype const c4 = t.y * B.vz - t.z * B.vy;
    ctype const c5 = -t.x * B.vz + t.z * B.vx;
    ctype const c6 = t.x * B.vy - t.y * B.vx;
    ctype const c7 = t.x * B.mx + t.y * B.my + t.z * B.mz;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(bivec,trivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vx * t.w + B.my * t.z - B.mz * t.y;
    ctype const c1 = B.vy * t.w - B.mx * t.z + B.mz * t.x;
    ctype const c2 = B.vz * t.w + B.mx * t.y - B.my * t.x;
    ctype const c3 = -B.vx * t.x - B.vy * t.y - B.vz * t.z;
    ctype const c4 = B.vy * t.z - B.vz * t.y;
    ctype const c5 = -B.vx * t.z + B.vz * t.x;
    ctype const c6 = B.vx * t.y - B.vy * t.x;
    ctype const c7 = -B.mx * t.x - B.my * t.y - B.mz * t.z;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(trivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w;
    ctype const c1 = -t.x * v.w;
    ctype const c2 = -t.y * v.w;
    ctype const c3 = -t.z * v.w;
    ctype const c4 = t.y * v.z - t.z * v.y;
    ctype const c5 = -t.x * v.z + t.z * v.x;
    ctype const c6 = t.x * v.y - t.y * v.x;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(vec,trivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w;
    ctype const c1 = -v.w * t.x;
    ctype const c2 = -v.w * t.y;
    ctype const c3 = -v.w * t.z;
    ctype const c4 = -v.y * t.z + v.z * t.y;
    ctype const c5 = v.x * t.z - v.z * t.x;
    ctype const c6 = -v.x * t.y + v.y * t.x;
    ctype const c7 = 0.0;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(trivec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -t.x * ctype(s);
    ctype const c1 = -t.y * ctype(s);
    ctype const c2 = -t.z * ctype(s);
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(s,trivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rgpr(Scalar3dp<T> s, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * t.x;
    ctype const c1 = ctype(s) * t.y;
    ctype const c2 = ctype(s) * t.z;
    ctype const c3 = 0.0;
    return Vec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(bivec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B1,
                                                   BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz - B1.mx * B2.vx -
                     B1.my * B2.vy - B1.mz * B2.vz;
    ctype const c1 = B1.vy * B2.vz - B1.vz * B2.vy;
    ctype const c2 = -B1.vx * B2.vz + B1.vz * B2.vx;
    ctype const c3 = B1.vx * B2.vy - B1.vy * B2.vx;
    ctype const c4 = B1.vy * B2.mz - B1.vz * B2.my + B1.my * B2.vz - B1.mz * B2.vy;
    ctype const c5 = -B1.vx * B2.mz + B1.vz * B2.mx - B1.mx * B2.vz + B1.mz * B2.vx;
    ctype const c6 = B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx;
    ctype const c7 = -B1.vx * B2.vx - B1.vy * B2.vy - B1.vz * B2.vz;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(bivec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.vy * v.z - B.vz * v.y + B.mx * v.w;
    ctype const c1 = -B.vx * v.z + B.vz * v.x + B.my * v.w;
    ctype const c2 = B.vx * v.y - B.vy * v.x + B.mz * v.w;
    ctype const c3 = 0.0;
    ctype const c4 = B.vx * v.w;
    ctype const c5 = B.vy * v.w;
    ctype const c6 = B.vz * v.w;
    ctype const c7 = -B.vx * v.x - B.vy * v.y - B.vz * v.z;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(vec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.vz - v.z * B.vy - v.w * B.mx;
    ctype const c1 = -v.x * B.vz + v.z * B.vx - v.w * B.my;
    ctype const c2 = v.x * B.vy - v.y * B.vx - v.w * B.mz;
    ctype const c3 = 0.0;
    ctype const c4 = v.w * B.vx;
    ctype const c5 = v.w * B.vy;
    ctype const c6 = v.w * B.vz;
    ctype const c7 = -v.x * B.vx - v.y * B.vy - v.z * B.vz;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(c0, c1, c2, c3),
                            TriVec3dp<ctype>(c4, c5, c6, c7));
}

// pga3dp rgpr :: rgpr(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = B.vx * ctype(s);
    ctype const c4 = B.vy * ctype(s);
    ctype const c5 = B.vz * ctype(s);
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rgpr :: rgpr(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rgpr(Scalar3dp<T> s, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = ctype(s) * B.vx;
    ctype const c4 = ctype(s) * B.vy;
    ctype const c5 = ctype(s) * B.vz;
    return BiVec3dp<ctype>(c0, c1, c2, c3, c4, c5);
}

// pga3dp rgpr :: rgpr(vec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v1,
                                                   Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = -v1.x * v2.w + v1.w * v2.x;
    ctype const c5 = -v1.y * v2.w + v1.w * v2.y;
    ctype const c6 = -v1.z * v2.w + v1.w * v2.z;
    ctype const c7 = -v1.w * v2.w;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(c0), BiVec3dp<ctype>(c1, c2, c3, c4, c5, c6),
                            PScalar3dp<ctype>(c7));
}

// pga3dp rgpr :: rgpr(vec,s) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = -v.w * ctype(s);
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(s,vec) -> trivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rgpr(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = ctype(s) * v.w;
    return TriVec3dp<ctype>(c0, c1, c2, c3);
}

// pga3dp rgpr :: rgpr(s,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rgpr([[maybe_unused]] Scalar3dp<T>,
                                                    [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(0.0);
}


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
inline Scalar3dp<T> inv(Scalar3dp<T> s)
{
    T sq_n = bulk_nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar3dp<T>(rev(s) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec3dp<T> inv(Vec3dp<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v) = v/bulk_sq_nrm(v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    T sq_n = bulk_nrm_sq(v);
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return Vec3dp<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline BiVec3dp<T> inv(BiVec3dp<T> const& B)
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
inline TriVec3dp<T> inv(TriVec3dp<T> const& t)
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

// due to the degenerate metric the pseudoscalar does not have an inverse

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline MVec3dp_E<T> inv(MVec3dp_E<T> const& E)
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
inline MVec3dp_U<T> inv(MVec3dp_U<T> const& U)
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
inline MVec3dp<T> inv(MVec3dp<T> const& M)
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


/////////////////////////////////////////////////////////////////////////////////////////
// Projective contractions for 3dp:
//
// REQUIRES: the dualized operand must be of smaller grade, or the result is zero
//
// l_bulk_contract3dp(a,b) = rwdg(l_bulk_dual(a), b)
// l_weight_contract3dp(a,b) = rwdg(l_weight_dual(a), b)
//
// r_bulk_contract3dp(a,b) = rwdg(a, r_bulk_dual(b))
// r_weight_contract3dp(a,b) = rwdg(a, r_weight_dual(b))
//
// The contraction subracts the grades of the operands.
//
// When the metric is the identity, these two contractions are identical and are
// called the contraction (they are identical with the right contraction).
// When the metric is degenerate they produce different results.
//
// In general a contraction throws away parts that are perpendicular to each other.
// The result of r_bulk_contract(B,v) lies in B and is perpendicular to v.
/////////////////////////////////////////////////////////////////////////////////////////

// REQUIRES: the dualized operand must be of smaller grade, or the result is zero

template <typename arg1, typename arg2>
decltype(auto) l_bulk_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(l_bulk_dual(a), b);
    return rwdg(l_bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) l_weight_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(l_weight_dual(a), b);
    return rwdg(l_weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) r_bulk_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, r_bulk_dual(b));
    return rwdg(std::forward<arg1>(a), r_bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) r_weight_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, r_weight_dual(b));
    return rwdg(std::forward<arg1>(a), r_weight_dual(std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// Projective expansions for 3dp:
//
// REQUIRES: the dualized operand must be of larger grade, or the result is zero
//
// l_bulk_expand3dp(a,b)   = wdg(l_bulk_dual(a),b)   (dual to l_weight_contract)
// l_weight_expand3dp(a,b) = wdg(l_weight_dual(a),b) (dual to l_bulk_contract)
//
// r_bulk_expand3dp(a,b)   = wdg(a,r_bulk_dual(b))   (dual to r_weight_contract)
// r_weight_expand3dp(a,b) = wdg(a,r_weight_dual(b)) (dual to r_bulk_contract)
//
// The expansion subtracts the antigrades of the objects.
//
// When the metric is the identity, these two expansions are identical and are
// called the expansion (they are identical with the left complement of the left
// contraction). When the metric is degenerate they produce different results.
//
// The expansion takes the parts of one object that are parallel to the other
// object and combines them with the space that is perpendicular to that other
// dualized object.
////////////////////////////////////////////////////////////////////////////////

// REQUIRES: the dualized operand must be of larger grade, or the result is zero

template <typename arg1, typename arg2>
decltype(auto) l_bulk_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(l_bulk_dual(a), b);
    return wdg(l_bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) l_weight_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(l_weight_dual(a), b);
    return wdg(l_weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) r_bulk_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(a, r_bulk_dual(b));
    return wdg(std::forward<arg1>(a), r_bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) r_weight_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(a, r_weight_dual(b));
    return wdg(std::forward<arg1>(a), r_weight_dual(std::forward<arg2>(b)));
}

} // namespace hd::ga::pga