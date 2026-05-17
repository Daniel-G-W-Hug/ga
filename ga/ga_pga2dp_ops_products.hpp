#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_pga2dp_ops_basics.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga2dp product operations:
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
// - rcmt()                  -> regressive commutator product (= asymmetric part of rgpr)
// - operator*()             -> geometric product (=gpr())
// - rgpr()                  -> regressive geometric product
// - inv()                   -> inversion operation (w.r.t. geometric product)
//
// - l_bulk_contract2dp()    -> left bulk contraction
// - l_weight_contract2dp()  -> left weight contraction
// - r_bulk_contract2dp()    -> right bulk contraction
// - r_weight_contract2dp()  -> right weight contraction
//
// Bulk expansions are the regressive versions of the corresponding weight contractions.
// Weight expansions are the regressive versions of the corresponding bulk contractions.
//
// - l_bulk_expand2dp()      -> left bulk expansion
// - l_weight_expand2dp()    -> left weight expansion
// - r_bulk_expand2dp()      -> right bulk expansion
// - r_weight_expand2dp()    -> right weight expansion
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot products for 2dp (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = 0
//
// Via the exomorphism G(a ^ b) = G(a) ^ G(b) the extended metric G is defined as a 2^n by
// 2^n matrix via the outer product, which in turn allows for definition of an extended
// dot product for all grades of a full multi-vector from scalar to pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

// pga2dp dot :: dot(mv,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp<T> const& A,
                                                  MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c6 * B.c6);
}

// pga2dp dot :: dot(mv_e,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp_E<T> const& A,
                                                  MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(A.c0 * B.c0 + A.c3 * B.c3);
}

// pga2dp dot :: dot(mv_u,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp_U<T> const& A,
                                                  MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1);
}

// pga2dp dot :: dot(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp dot :: dot(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(BiVec2dp<T> const& B1,
                                                  BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(B1.z * B2.z);
}

// pga2dp dot :: dot(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(Vec2dp<T> const& v1,
                                                  Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// pga2dp dot :: dot(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive dot products for 2dp (= defined for equal grades exclusively)
//
// rdot(v1,v2) = cmpl( dot(cmpl(v1),cmpl(v2)) )
//
// returns a pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

// pga2dp rdot :: rdot(mv,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(MVec2dp<T> const& A,
                                                    MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(A.c3 * B.c3 + A.c4 * B.c4 + A.c5 * B.c5 + A.c7 * B.c7);
}

// pga2dp rdot :: rdot(mv_e,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(MVec2dp_E<T> const& A,
                                                    MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(A.c1 * B.c1 + A.c2 * B.c2);
}

// pga2dp rdot :: rdot(mv_u,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(MVec2dp_U<T> const& A,
                                                    MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(A.c2 * B.c2 + A.c3 * B.c3);
}

// pga2dp rdot :: rdot(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(PScalar2dp<T> ps1, PScalar2dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps1) * ctype(ps2));
}

// pga2dp rdot :: rdot(bivec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(BiVec2dp<T> const& B1,
                                                    BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(B1.x * B2.x + B1.y * B2.y);
}

// pga2dp rdot :: rdot(vec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(Vec2dp<T> const& v1,
                                                    Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(v1.z * v2.z);
}

// pga2dp rdot :: rdot(s,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot([[maybe_unused]] Scalar2dp<T>,
                                                    [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////

// pga2dp wdg :: wdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c3 + A.c3 * B.c1 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 - A.c2 * B.c3 + A.c3 * B.c2 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 - A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 - A.c6 * B.c3 + A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& A,
                                                MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c1 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c2 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c3 + A.c6 * B.c0;
    ctype const c7 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c3 + A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp_E<T> const& A,
                                                MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c2 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c3 * B.c0;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& A,
                                                MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c0 * B.c1;
    ctype const c3 = A.c0 * B.c2;
    ctype const c4 = -A.c1 * B.c2 + A.c3 * B.c0;
    ctype const c5 = -A.c2 * B.c2 + A.c3 * B.c1;
    ctype const c6 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype const c7 = A.c0 * B.c3 - A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp_U<T> const& A,
                                                MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c1 * B.c0;
    ctype const c3 = A.c2 * B.c0;
    ctype const c4 = -A.c0 * B.c3 + A.c2 * B.c1;
    ctype const c5 = -A.c1 * B.c3 + A.c2 * B.c2;
    ctype const c6 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c7 = A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c6 + A.c3 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& M, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(M.c0 * ctype(ps));
}

// pga2dp wdg :: wdg(ps,mv) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(PScalar2dp<T> ps, MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * M.c0);
}

// pga2dp wdg :: wdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& M, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = M.c0 * B.x;
    ctype const c5 = M.c0 * B.y;
    ctype const c6 = M.c0 * B.z;
    ctype const c7 = M.c1 * B.y - M.c2 * B.x - M.c3 * B.z;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B, MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = B.x * M.c0;
    ctype const c5 = B.y * M.c0;
    ctype const c6 = B.z * M.c0;
    ctype const c7 = -B.x * M.c2 + B.y * M.c1 - B.z * M.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& M, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = -M.c1 * v.z + M.c3 * v.x;
    ctype const c5 = -M.c2 * v.z + M.c3 * v.y;
    ctype const c6 = M.c1 * v.y - M.c2 * v.x;
    ctype const c7 = -M.c4 * v.y + M.c5 * v.x - M.c6 * v.z;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v, MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = -v.x * M.c3 + v.z * M.c1;
    ctype const c5 = -v.y * M.c3 + v.z * M.c2;
    ctype const c6 = v.x * M.c2 - v.y * M.c1;
    ctype const c7 = v.x * M.c5 - v.y * M.c4 - v.z * M.c6;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& M, Scalar2dp<U> s)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, MVec2dp<U> const& M)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp wdg :: wdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> wdg(MVec2dp_E<T> const& A,
                                                  MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp wdg :: wdg(mv_e,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> wdg(MVec2dp_E<T> const& A,
                                                  MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c2;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp wdg :: wdg(mv_u,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> wdg(MVec2dp_U<T> const& A,
                                                  MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c3 + A.c3 * B.c0;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp wdg :: wdg(mv_e,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(MVec2dp_E<T> const& M,
                                                   PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(M.c0 * ctype(ps));
}

// pga2dp wdg :: wdg(ps,mv_e) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(PScalar2dp<T> ps,
                                                   MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * M.c0);
}

// pga2dp wdg :: wdg(mv_e,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(MVec2dp_E<T> const& M,
                                                 BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * B.x;
    ctype const c1 = M.c0 * B.y;
    ctype const c2 = M.c0 * B.z;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(bivec,mv_e) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B,
                                                 MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c0;
    ctype const c1 = B.y * M.c0;
    ctype const c2 = B.z * M.c0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(mv_e,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> wdg(MVec2dp_E<T> const& M,
                                                  Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    ctype const c3 = -M.c1 * v.y + M.c2 * v.x - M.c3 * v.z;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp wdg :: wdg(vec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v,
                                                  MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    ctype const c3 = v.x * M.c2 - v.y * M.c1 - v.z * M.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp wdg :: wdg(mv_e,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> wdg(MVec2dp_E<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp wdg :: wdg(s,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp wdg :: wdg(mv_u,mv_u) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(MVec2dp_U<T> const& A,
                                                 MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c1 = -A.c1 * B.c2 + A.c2 * B.c1;
    ctype const c2 = A.c0 * B.c1 - A.c1 * B.c0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(mv_u,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] MVec2dp_U<T> const&,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(ps,mv_u) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] MVec2dp_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(mv_u,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(MVec2dp_U<T> const& M,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(M.c0 * B.y - M.c1 * B.x - M.c2 * B.z);
}

// pga2dp wdg :: wdg(bivec,mv_u) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B,
                                                   MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-B.x * M.c1 + B.y * M.c0 - B.z * M.c2);
}

// pga2dp wdg :: wdg(mv_u,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(MVec2dp_U<T> const& M,
                                                 Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c0 * v.z + M.c2 * v.x;
    ctype const c1 = -M.c1 * v.z + M.c2 * v.y;
    ctype const c2 = M.c0 * v.y - M.c1 * v.x;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(vec,mv_u) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v,
                                                 MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.x * M.c2 + v.z * M.c0;
    ctype const c1 = -v.y * M.c2 + v.z * M.c1;
    ctype const c2 = v.x * M.c1 - v.y * M.c0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(mv_u,s) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> wdg(MVec2dp_U<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp wdg :: wdg(s,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp wdg :: wdg(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec2dp<T> const&,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec2dp<T> const&,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// pga2dp wdg :: wdg(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// pga2dp wdg :: wdg(bivec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec2dp<T> const&,
                                                  [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp wdg :: wdg(bivec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-B.x * v.y + B.y * v.x - B.z * v.z);
}

// pga2dp wdg :: wdg(vec,bivec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(v.x * B.y - v.y * B.x - v.z * B.z);
}

// pga2dp wdg :: wdg(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(s);
    ctype const c1 = B.y * ctype(s);
    ctype const c2 = B.z * ctype(s);
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.x;
    ctype const c1 = ctype(s) * B.y;
    ctype const c2 = ctype(s) * B.z;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
// pga2dp wdg :: wdg(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c1 = -v1.y * v2.z + v1.z * v2.y;
    ctype const c2 = v1.x * v2.y - v1.y * v2.x;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp wdg :: wdg(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg for Point2d
////////////////////////////////////////////////////////////////////////////////

// wedge product between two points (aka vectors with implicit v.z == 1)
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Line2d<std::common_type_t<T, U>> wdg(Point2d<T> const& p, Point2d<U> const& q)
{
    return Line2d<std::common_type_t<T, U>>(q.x - p.x, q.y - p.y, p.x * q.y - p.y * q.x);
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg -> join
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> join(BiVec2dp<T> const& B,
                                                    Vec2dp<U> const& v)
{
    return wdg(B, v);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& v,
                                                    BiVec2dp<U> const& B)
{
    return wdg(v, B);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& v1,
                                                  Vec2dp<U> const& v2)
{
    return wdg(v1, v2);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Line2d<std::common_type_t<T, U>> join(Point2d<T> const& p, Point2d<U> const& q)
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
//        rwdg(ul, ur) = cmpl(wdg(cmpl(ul),cmpl(ur))) = cmpl(cmpl(ul) ^ cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// pga2dp rwdg :: rwdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 - A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 - A.c6 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c1 * B.c7 - A.c4 * B.c6 + A.c6 * B.c4 + A.c7 * B.c1;
    ctype const c2 = A.c2 * B.c7 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c2;
    ctype const c3 = A.c3 * B.c7 - A.c4 * B.c5 + A.c5 * B.c4 + A.c7 * B.c3;
    ctype const c4 = A.c4 * B.c7 + A.c7 * B.c4;
    ctype const c5 = A.c5 * B.c7 + A.c7 * B.c5;
    ctype const c6 = A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = A.c7 * B.c7;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& A,
                                                 MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c3 + A.c7 * B.c0;
    ctype const c1 = -A.c4 * B.c3 + A.c6 * B.c1;
    ctype const c2 = -A.c5 * B.c3 + A.c6 * B.c2;
    ctype const c3 = -A.c4 * B.c2 + A.c5 * B.c1;
    ctype const c4 = A.c7 * B.c1;
    ctype const c5 = A.c7 * B.c2;
    ctype const c6 = A.c7 * B.c3;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp_E<T> const& A,
                                                 MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c3;
    ctype const c1 = -A.c1 * B.c6 + A.c3 * B.c4;
    ctype const c2 = -A.c2 * B.c6 + A.c3 * B.c5;
    ctype const c3 = -A.c1 * B.c5 + A.c2 * B.c4;
    ctype const c4 = A.c1 * B.c7;
    ctype const c5 = A.c2 * B.c7;
    ctype const c6 = A.c3 * B.c7;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& A,
                                                 MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 - A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2;
    ctype const c1 = A.c1 * B.c3 + A.c7 * B.c0;
    ctype const c2 = A.c2 * B.c3 + A.c7 * B.c1;
    ctype const c3 = A.c3 * B.c3 + A.c7 * B.c2;
    ctype const c4 = A.c4 * B.c3;
    ctype const c5 = A.c5 * B.c3;
    ctype const c6 = A.c6 * B.c3;
    ctype const c7 = A.c7 * B.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp_U<T> const& A,
                                                 MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c6 + A.c3 * B.c0;
    ctype const c1 = A.c0 * B.c7 + A.c3 * B.c1;
    ctype const c2 = A.c1 * B.c7 + A.c3 * B.c2;
    ctype const c3 = A.c2 * B.c7 + A.c3 * B.c3;
    ctype const c4 = A.c3 * B.c4;
    ctype const c5 = A.c3 * B.c5;
    ctype const c6 = A.c3 * B.c6;
    ctype const c7 = A.c3 * B.c7;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& M, PScalar2dp<U> ps)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps, MVec2dp<U> const& M)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& M,
                                                 BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.y - M.c2 * B.x - M.c3 * B.z;
    ctype const c1 = -M.c4 * B.z + M.c6 * B.x;
    ctype const c2 = -M.c5 * B.z + M.c6 * B.y;
    ctype const c3 = -M.c4 * B.y + M.c5 * B.x;
    ctype const c4 = M.c7 * B.x;
    ctype const c5 = M.c7 * B.y;
    ctype const c6 = M.c7 * B.z;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B,
                                                 MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c2 + B.y * M.c1 - B.z * M.c3;
    ctype const c1 = -B.x * M.c6 + B.z * M.c4;
    ctype const c2 = -B.y * M.c6 + B.z * M.c5;
    ctype const c3 = -B.x * M.c5 + B.y * M.c4;
    ctype const c4 = B.x * M.c7;
    ctype const c5 = B.y * M.c7;
    ctype const c6 = B.z * M.c7;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& M, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * v.y + M.c5 * v.x - M.c6 * v.z;
    ctype const c1 = M.c7 * v.x;
    ctype const c2 = M.c7 * v.y;
    ctype const c3 = M.c7 * v.z;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& v, MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c5 - v.y * M.c4 - v.z * M.c6;
    ctype const c1 = v.x * M.c7;
    ctype const c2 = v.y * M.c7;
    ctype const c3 = v.z * M.c7;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rwdg :: rwdg(mv,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(M.c7 * ctype(s));
}

// pga2dp rwdg :: rwdg(s,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(Scalar2dp<T> s, MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s) * M.c7);
}

// pga2dp rwdg :: rwdg(mv_e,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(MVec2dp_E<T> const& A,
                                                MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c1 = -A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c2 = -A.c1 * B.c2 + A.c2 * B.c1;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(mv_e,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rwdg(MVec2dp_E<T> const& A,
                                                   MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c2;
    ctype const c1 = A.c1 * B.c3;
    ctype const c2 = A.c2 * B.c3;
    ctype const c3 = A.c3 * B.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rwdg :: rwdg(mv_u,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rwdg(MVec2dp_U<T> const& A,
                                                   MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c3 + A.c3 * B.c0;
    ctype const c1 = A.c3 * B.c1;
    ctype const c2 = A.c3 * B.c2;
    ctype const c3 = A.c3 * B.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rwdg :: rwdg(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rwdg(MVec2dp_E<T> const& M,
                                                   PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rwdg :: rwdg(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps,
                                                   MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rwdg :: rwdg(mv_e,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(MVec2dp_E<T> const& M,
                                                BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.z + M.c3 * B.x;
    ctype const c1 = -M.c2 * B.z + M.c3 * B.y;
    ctype const c2 = -M.c1 * B.y + M.c2 * B.x;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(bivec,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B,
                                                MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c3 + B.z * M.c1;
    ctype const c1 = -B.y * M.c3 + B.z * M.c2;
    ctype const c2 = -B.x * M.c2 + B.y * M.c1;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(mv_e,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(MVec2dp_E<T> const& M,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(-M.c1 * v.y + M.c2 * v.x - M.c3 * v.z);
}

// pga2dp rwdg :: rwdg(vec,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& v,
                                                   MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(v.x * M.c2 - v.y * M.c1 - v.z * M.c3);
}

// pga2dp rwdg :: rwdg(mv_e,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] MVec2dp_E<T> const&,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rwdg :: rwdg(s,mv_e) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] MVec2dp_E<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rwdg :: rwdg(mv_u,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rwdg(MVec2dp_U<T> const& A,
                                                   MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c1 = A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c2 = A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c3 = A.c3 * B.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rwdg :: rwdg(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rwdg(MVec2dp_U<T> const& M,
                                                   PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rwdg :: rwdg(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps,
                                                   MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rwdg :: rwdg(mv_u,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rwdg(MVec2dp_U<T> const& M,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * B.y - M.c1 * B.x - M.c2 * B.z;
    ctype const c1 = M.c3 * B.x;
    ctype const c2 = M.c3 * B.y;
    ctype const c3 = M.c3 * B.z;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rwdg :: rwdg(bivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B,
                                                   MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c1 + B.y * M.c0 - B.z * M.c2;
    ctype const c1 = B.x * M.c3;
    ctype const c2 = B.y * M.c3;
    ctype const c3 = B.z * M.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rwdg :: rwdg(mv_u,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(MVec2dp_U<T> const& M, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c3 * v.x;
    ctype const c1 = M.c3 * v.y;
    ctype const c2 = M.c3 * v.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(vec,mv_u) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& v, MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c3;
    ctype const c1 = v.y * M.c3;
    ctype const c2 = v.z * M.c3;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(mv_u,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(MVec2dp_U<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(M.c3 * ctype(s));
}

// pga2dp rwdg :: rwdg(s,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(Scalar2dp<T> s, MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s) * M.c3);
}

// pga2dp rwdg :: rwdg(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps1, PScalar2dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps1) * ctype(ps2));
}

// pga2dp rwdg :: rwdg(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.x;
    ctype const c1 = ctype(ps) * B.y;
    ctype const c2 = ctype(ps) * B.z;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(ps);
    ctype const c1 = B.y * ctype(ps);
    ctype const c2 = B.z * ctype(ps);
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& v, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(ps) * ctype(s));
}

// pga2dp rwdg :: rwdg(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s) * ctype(ps));
}

// pga2dp rwdg :: rwdg(bivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.x * B2.z + B1.z * B2.x;
    ctype const c1 = -B1.y * B2.z + B1.z * B2.y;
    ctype const c2 = -B1.x * B2.y + B1.y * B2.x;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rwdg :: rwdg(bivec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(-B.x * v.y + B.y * v.x - B.z * v.z);
}

// pga2dp rwdg :: rwdg(vec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& v,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(v.x * B.y - v.y * B.x - v.z * B.z);
}

// pga2dp rwdg :: rwdg(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec2dp<T> const&,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rwdg :: rwdg(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rwdg :: rwdg(vec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2dp<T> const&,
                                                   [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rwdg :: rwdg(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2dp<T> const&,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rwdg :: rwdg(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rwdg :: rwdg(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions rwdg -> meet
////////////////////////////////////////////////////////////////////////////////

// intersection of two lines -> return a point
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> meet(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{
    return rwdg(B1, B2);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> meet(Line2d<T> const& l1, Line2d<U> const& l2)
{
    return rwdg(static_cast<BiVec2dp<T> const&>(l1), static_cast<BiVec2dp<U> const&>(l2));
}


////////////////////////////////////////////////////////////////////////////////
// transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// twdg1(ps,vec) = bivec -> identical to geometric product gpr(ps,vec)
//                       -> and identical to the left contraction vec << ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> twdg1(PScalar2dp<T> ps, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.y;
    ctype const c1 = ctype(ps) * v.x;
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// twdg1(vec,ps) = bivec -> identical to geometric product gpr(vec,ps)
//                       -> and identical to the right contraction ps >> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> twdg1(Vec2dp<T> const& v, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(ps);
    ctype const c1 = v.x * ctype(ps);
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// twdg1(bivec,bivec) = bivector -> identical to commutator product cmt(bivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> twdg1(BiVec2dp<T> const& B1,
                                                   BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.y * B2.z + B1.z * B2.y;
    ctype const c1 = B1.x * B2.z - B1.z * B2.x;
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// twdg1(bivec,vec) = vector -> identical to left contraction operator<<(vec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> twdg1(BiVec2dp<T> const& B, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.z * v.y;
    ctype const c1 = -B.z * v.x;
    ctype const c2 = B.x * v.x + B.y * v.y;
    return Vec2dp<ctype>(c0, c1, c2);
}

// twdg1(vec,bivec) = vector -> identical to right contraction operator>>(bivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> twdg1(Vec2dp<T> const& v, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.z;
    ctype const c1 = v.x * B.z;
    ctype const c2 = -v.x * B.x - v.y * B.y;
    return Vec2dp<ctype>(c0, c1, c2);
}

// twdg1(vec,vec) = scalar -> identical to dot product dot(vec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> twdg1(Vec2dp<T> const& v1,
                                                    Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(v1.x * v2.x + v1.y * v2.y);
}


////////////////////////////////////////////////////////////////////////////////
// regressive transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// rtwdg1(vec,vec) = pseudoscalar -> identical to regressive dot product rdot(vec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rtwdg1(BiVec2dp<T> const& B1,
                                                      BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(B1.x * B2.x + B1.y * B2.y);
}

// rtwdg1(bivec,vec) = bivec -> identical to r_weight_expand2dp(vec,bivec)
//                           ->            = vec ^ weight_dual(bivec)
//                           ->            = rcmt(bivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rtwdg1(BiVec2dp<T> const& B,
                                                    Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * v.z;
    ctype const c1 = -B.x * v.z;
    ctype const c2 = -B.x * v.x - B.y * v.y;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// rtwdg1(vec,bivec) = bivec -> identical to l_weight_expand2dp(bivec,vec)
//                           ->            = weight_dual(bivec) ^ vec
//                           ->            = rcmt(vec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rtwdg1(Vec2dp<T> const& v,
                                                    BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.z * B.y;
    ctype const c1 = v.z * B.x;
    ctype const c2 = v.x * B.x + v.y * B.y;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// rtwdg1(bivec,s) = vec -> identical to r_weight_expand2dp(s,bivec)
//                           ->        = s ^ weight_dual(bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rtwdg1(BiVec2dp<T> const& B, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * ctype(s);
    ctype const c1 = -B.x * ctype(s);
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}

// rtwdg1(s,bivec) -> identical to l_weight_expand2dp(bivec,s)
//                           ->  = weight_dual(bivec) ^ s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rtwdg1(Scalar2dp<T> s, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.y;
    ctype const c1 = -ctype(s) * B.x;
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}

// rtwdg1(vec1,vec2) -> identical to regressive commutator product rcmt(vec1,vec2)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rtwdg1(Vec2dp<T> const& v1,
                                                  Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.y * v2.z - v1.z * v2.y;
    ctype const c1 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// Implements the left bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator<<(a,b) = l_bulk_contract(a,b) = rwdg( bulk_dual(a), b )
//
////////////////////////////////////////////////////////////////////////////////

// pga2dp l_bulk_contract :: l_bulk_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator<<(MVec2dp<T> const& A,
                                                       MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c6 * B.c6;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c6;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c6;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c4 + A.c2 * B.c5 - A.c6 * B.c7;
    ctype const c4 = A.c0 * B.c4 - A.c2 * B.c7;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c7;
    ctype const c6 = A.c0 * B.c6;
    ctype const c7 = A.c0 * B.c7;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp l_bulk_contract :: l_bulk_contract(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2dp<T>,
                                                         [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp l_bulk_contract :: l_bulk_contract(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp l_bulk_contract :: l_bulk_contract(bivec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator<<(BiVec2dp<T> const& B,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = -B.z * ctype(ps);
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp l_bulk_contract :: l_bulk_contract(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp l_bulk_contract :: l_bulk_contract(vec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v,
                                                        PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(ps);
    ctype const c1 = v.x * ctype(ps);
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp l_bulk_contract :: l_bulk_contract(ps,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2dp<T>,
                                                         [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp l_bulk_contract :: l_bulk_contract(s,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s,
                                                          PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// left contraction (B1 << B2) == dot(B1,B2) == (B1 >> B2))
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(BiVec2dp<T> const& B1,
                                                         BiVec2dp<U> const& B2)
{
    // using ctype = std::common_type_t<T, U>;
    // return Scalar2dp<ctype>(B1.z * B2.z);
    return dot(B1, B2);
}

// pga2dp l_bulk_contract :: l_bulk_contract(bivec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// (v << B) == cmt(B,v)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.z;
    ctype const c1 = -v.x * B.z;
    ctype const c2 = v.x * B.x + v.y * B.y;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp l_bulk_contract :: l_bulk_contract(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp l_bulk_contract :: l_bulk_contract(s,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.x;
    ctype const c1 = ctype(s) * B.y;
    ctype const c2 = ctype(s) * B.z;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// (v1 << v2) == dot(v1,v2) == (v1 >> v2))
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v1,
                                                         Vec2dp<U> const& v2)
{
    // using ctype = std::common_type_t<T, U>;
    // return Scalar2dp<ctype>(v1.x * v2.x + v1.y * v2.y);
    return dot(v1, v2);
}

// pga2dp l_bulk_contract :: l_bulk_contract(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// identical with scalar multiplication of the vector
// pga2dp l_bulk_contract :: l_bulk_contract(s,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp l_bulk_contract :: l_bulk_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The resulting object lies in A and is perpendicular to B
//
// Implements the right bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator>>(a,b) = r_bulk_contract(a,b) = rwdg( a, bulk_dual(b) )
//
////////////////////////////////////////////////////////////////////////////////

// pga2dp r_bulk_contract :: r_bulk_contract(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator>>(MVec2dp<T> const& A,
                                                       MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c6 * B.c6;
    ctype const c1 = A.c1 * B.c0 - A.c6 * B.c2;
    ctype const c2 = A.c2 * B.c0 + A.c6 * B.c1;
    ctype const c3 = A.c3 * B.c0 - A.c4 * B.c1 - A.c5 * B.c2 - A.c7 * B.c6;
    ctype const c4 = A.c4 * B.c0 - A.c7 * B.c2;
    ctype const c5 = A.c5 * B.c0 + A.c7 * B.c1;
    ctype const c6 = A.c6 * B.c0;
    ctype const c7 = A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp r_bulk_contract :: r_bulk_contract(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] PScalar2dp<T>,
                                                         [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp r_bulk_contract :: r_bulk_contract(ps,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = -ctype(ps) * B.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp r_bulk_contract :: r_bulk_contract(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp r_bulk_contract :: r_bulk_contract(ps,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                        Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.y;
    ctype const c1 = ctype(ps) * v.x;
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp r_bulk_contract :: r_bulk_contract(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp r_bulk_contract :: r_bulk_contract(ps,s) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                          Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// pga2dp r_bulk_contract :: r_bulk_contract(s,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar2dp<T>,
                                                         [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// (B1 >> B2) == product dot(B1,B2) == (B1 >> B2))
// pga2dp r_bulk_contract :: r_bulk_contract(bivec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B1,
                                                         BiVec2dp<U> const& B2)
{
    // using ctype = std::common_type_t<T, U>;
    // return Scalar2dp<ctype>(B1.z * B2.z);
    return dot(B1, B2);
}

// (B >> v) == cmt(v,B))
// pga2dp r_bulk_contract :: r_bulk_contract(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B,
                                                      Vec2dp<U> const& v)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.z * v.y;
    ctype const c1 = B.z * v.x;
    ctype const c2 = -B.x * v.x - B.y * v.y;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp r_bulk_contract :: r_bulk_contract(vec,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// identical with scalar multiplication of the bivector
// pga2dp r_bulk_contract :: r_bulk_contract(bivec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B,
                                                        Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(s);
    ctype const c1 = B.y * ctype(s);
    ctype const c2 = B.z * ctype(s);
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp r_bulk_contract :: r_bulk_contract(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dp<T>, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// (v1 >> v2) == dot(v1,v2) == (v1 << v2))
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(Vec2dp<T> const& v1,
                                                         Vec2dp<U> const& v2)
{
    // using ctype = std::common_type_t<T, U>;
    // return Scalar2dp<ctype>(v1.x * v2.x + v1.y * v2.y);
    return dot(v1, v2);
}

// identical with scalar multiplication of the vector
// pga2dp r_bulk_contract :: r_bulk_contract(vec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator>>(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    return Vec2dp<ctype>(c0, c1, c2);
    ;
}

// pga2dp r_bulk_contract :: r_bulk_contract(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dp<T>, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp r_bulk_contract :: r_bulk_contract(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

// pga2dp cmt :: cmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> cmt(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c2 * B.c6 + A.c6 * B.c2;
    ctype const c2 = A.c1 * B.c6 - A.c6 * B.c1;
    ctype const c3 = -A.c1 * B.c4 - A.c2 * B.c5 + A.c4 * B.c1 + A.c5 * B.c2;
    ctype const c4 = -A.c1 * B.c3 + A.c3 * B.c1 - A.c5 * B.c6 + A.c6 * B.c5;
    ctype const c5 = -A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c6 - A.c6 * B.c4;
    ctype const c6 = A.c1 * B.c2 - A.c2 * B.c1;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp cmt :: cmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec2dp<T> const&,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec2dp<T> const&,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(ps,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2dp<T>,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(s,ps) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(bivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> cmt(BiVec2dp<T> const& B1,
                                                 BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.y * B2.z + B1.z * B2.y;
    ctype const c1 = B1.x * B2.z - B1.z * B2.x;
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// cmt(B,v) == -cmt(v,B) == (v << B)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> cmt(BiVec2dp<T> const& B, Vec2dp<U> const& v)
{
    // using ctype = std::common_type_t<T, U>;
    // ctype const c0 = B.z * v.y;
    // ctype const c1 = -B.z * v.x;
    // ctype const c2 = B.x * v.x + B.y * v.y;
    // return Vec2dp<ctype>(c0, c1, c2);
    return (v << B);
}

// cmt(v,B) = -cmt(B,v) == (B >> v)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> cmt(Vec2dp<T> const& v, BiVec2dp<U> const& B)
{
    // using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -v.y * B.z;
    // ctype const c1 = v.x * B.z;
    // ctype const c2 = -v.x * B.x - v.y * B.y;
    // return Vec2dp<ctype>(c0, c1, c2);
    return (B >> v);
}

// pga2dp cmt :: cmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec2dp<T> const&,
                                                  [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dp<T>,
                                                  [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// cmt(v1,v2) = wdg(v1,v2)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> cmt(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    // using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -v1.x * v2.z + v1.z * v2.x;
    // ctype const c1 = -v1.y * v2.z + v1.z * v2.y;
    // ctype const c2 = v1.x * v2.y - v1.y * v2.x;
    // return BiVec2dp<ctype>(c0, c1, c2);
    return wdg(v1, v2);
}

// pga2dp cmt :: cmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec2dp<T> const&,
                                                  [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dp<T>,
                                                  [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp cmt :: cmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2dp<T>,
                                                  [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////

// pga2dp rcmt :: rcmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rcmt(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c2 * B.c3 - A.c3 * B.c2 - A.c4 * B.c6 + A.c6 * B.c4;
    ctype const c2 = -A.c1 * B.c3 + A.c3 * B.c1 - A.c5 * B.c6 + A.c6 * B.c5;
    ctype const c3 = -A.c4 * B.c5 + A.c5 * B.c4;
    ctype const c4 = -A.c3 * B.c5 + A.c5 * B.c3;
    ctype const c5 = A.c3 * B.c4 - A.c4 * B.c3;
    ctype const c6 = A.c1 * B.c4 + A.c2 * B.c5 - A.c4 * B.c1 - A.c5 * B.c2;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rcmt :: rcmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dp<T>,
                                                   [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dp<T>,
                                                   [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec2dp<T> const&,
                                                   [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dp<T>,
                                                   [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Vec2dp<T> const&,
                                                   [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(ps,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] PScalar2dp<T>,
                                                    [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(s,ps) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar2dp<T>,
                                                    [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
}

// rcmt(B1,B2) == rwdg(B2,B1)
// pga2dp rcmt :: rcmt(bivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rcmt(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{
    // using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -B1.x * B2.z + B1.z * B2.x;
    // ctype const c1 = -B1.y * B2.z + B1.z * B2.y;
    // ctype const c2 = -B1.x * B2.y + B1.y * B2.x;
    // return Vec2dp<ctype>(c0, c1, c2);
    return rwdg(B1, B2);
}

// pga2dp rcmt :: rcmt(bivec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rcmt(BiVec2dp<T> const& B,
                                                  Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * v.z;
    ctype const c1 = -B.x * v.z;
    ctype const c2 = -B.x * v.x - B.y * v.y;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rcmt :: rcmt(vec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rcmt(Vec2dp<T> const& v,
                                                  BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.z * B.y;
    ctype const c1 = v.z * B.x;
    ctype const c2 = v.x * B.x + v.y * B.y;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rcmt :: rcmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] BiVec2dp<T> const&,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(vec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rcmt(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.y * v2.z - v1.z * v2.y;
    ctype const c1 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rcmt :: rcmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Vec2dp<T> const&,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// pga2dp rcmt :: rcmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rcmt([[maybe_unused]] Scalar2dp<T>,
                                                   [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 2dp multivector
// Expensive! - Don't use if you don't have to! (48x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// a * b = dot(a,b) + wdg(a,b) (vector*vector = scalar + bivector)
// A * b = (b << A) + wdg(A,b) (bivector*vector = vector + trivector)
// a * B = (B >> a) + wdg(a,B) (vector*bivector = vector + trivector)
//
// multivector * multivector => multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                      MVec2dp<U> const& B)
{
    // geometric product of two fully populated 2dp multivectors
    // => due to the degenerate algebra 16 terms are not present in G<2,0,1> compared to
    // G<3,0,0>
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c6 * B.c6;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c6 * B.c2;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c6 * B.c1;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c4 - A.c2 * B.c5 + A.c3 * B.c0 + A.c4 * B.c1 +
                     A.c5 * B.c2 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c4 = A.c0 * B.c4 - A.c1 * B.c3 - A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c0 -
                     A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c2;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c7 - A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c6 +
                     A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c1;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 - A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 - A.c6 * B.c3 + A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of a multivector with an even-grade multivector
// multivector * even-grade multivector => multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                      MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c6 * B.c3;
    ctype const c1 = A.c1 * B.c0 - A.c2 * B.c3;
    ctype const c2 = A.c1 * B.c3 + A.c2 * B.c0;
    ctype const c3 = -A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c0 + A.c7 * B.c3;
    ctype const c4 = A.c0 * B.c1 + A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype const c5 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1;
    ctype const c6 = A.c0 * B.c3 + A.c6 * B.c0;
    ctype const c7 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c3 + A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of an even-grade multivector with a multivector
// even-grade multivector * multivector => multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                      MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c3 * B.c6;
    ctype const c1 = A.c0 * B.c1 + A.c3 * B.c2;
    ctype const c2 = A.c0 * B.c2 - A.c3 * B.c1;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c7;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4;
    ctype const c6 = A.c0 * B.c6 + A.c3 * B.c0;
    ctype const c7 = A.c0 * B.c7 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// multivector * odd multivector => multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                      MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1;
    ctype const c1 = A.c0 * B.c0 + A.c6 * B.c1;
    ctype const c2 = A.c0 * B.c1 - A.c6 * B.c0;
    ctype const c3 = A.c0 * B.c2 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c3;
    ctype const c4 = -A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0 - A.c7 * B.c1;
    ctype const c5 = A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1 + A.c7 * B.c0;
    ctype const c6 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype const c7 = A.c0 * B.c3 - A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// odd multivector * multivector => multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                      MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2;
    ctype const c1 = A.c0 * B.c0 - A.c1 * B.c6;
    ctype const c2 = A.c0 * B.c6 + A.c1 * B.c0;
    ctype const c3 = -A.c0 * B.c4 - A.c1 * B.c5 + A.c2 * B.c0 + A.c3 * B.c6;
    ctype const c4 = -A.c0 * B.c3 - A.c1 * B.c7 + A.c2 * B.c1 - A.c3 * B.c2;
    ctype const c5 = A.c0 * B.c7 - A.c1 * B.c3 + A.c2 * B.c2 + A.c3 * B.c1;
    ctype const c6 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c7 = A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c6 + A.c3 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& M,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = M.c6 * ctype(ps);
    ctype const c4 = -M.c2 * ctype(ps);
    ctype const c5 = M.c1 * ctype(ps);
    ctype const c6 = 0.0;
    ctype const c7 = M.c0 * ctype(ps);
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                      MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = ctype(ps) * M.c6;
    ctype const c4 = -ctype(ps) * M.c2;
    ctype const c5 = ctype(ps) * M.c1;
    ctype const c6 = 0.0;
    ctype const c7 = ctype(ps) * M.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& M,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c6 * B.z;
    ctype const c1 = -M.c2 * B.z;
    ctype const c2 = M.c1 * B.z;
    ctype const c3 = -M.c1 * B.x - M.c2 * B.y + M.c7 * B.z;
    ctype const c4 = M.c0 * B.x - M.c5 * B.z + M.c6 * B.y;
    ctype const c5 = M.c0 * B.y + M.c4 * B.z - M.c6 * B.x;
    ctype const c6 = M.c0 * B.z;
    ctype const c7 = M.c1 * B.y - M.c2 * B.x - M.c3 * B.z;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                      MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.z * M.c6;
    ctype const c1 = B.z * M.c2;
    ctype const c2 = -B.z * M.c1;
    ctype const c3 = B.x * M.c1 + B.y * M.c2 + B.z * M.c7;
    ctype const c4 = B.x * M.c0 - B.y * M.c6 + B.z * M.c5;
    ctype const c5 = B.x * M.c6 + B.y * M.c0 - B.z * M.c4;
    ctype const c6 = B.z * M.c0;
    ctype const c7 = -B.x * M.c2 + B.y * M.c1 - B.z * M.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& M,
                                                      Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y;
    ctype const c1 = M.c0 * v.x + M.c6 * v.y;
    ctype const c2 = M.c0 * v.y - M.c6 * v.x;
    ctype const c3 = M.c0 * v.z + M.c4 * v.x + M.c5 * v.y;
    ctype const c4 = -M.c1 * v.z + M.c3 * v.x - M.c7 * v.y;
    ctype const c5 = -M.c2 * v.z + M.c3 * v.y + M.c7 * v.x;
    ctype const c6 = M.c1 * v.y - M.c2 * v.x;
    ctype const c7 = -M.c4 * v.y + M.c5 * v.x - M.c6 * v.z;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                      MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2;
    ctype const c1 = v.x * M.c0 - v.y * M.c6;
    ctype const c2 = v.x * M.c6 + v.y * M.c0;
    ctype const c3 = -v.x * M.c4 - v.y * M.c5 + v.z * M.c0;
    ctype const c4 = -v.x * M.c3 - v.y * M.c7 + v.z * M.c1;
    ctype const c5 = v.x * M.c7 - v.y * M.c3 + v.z * M.c2;
    ctype const c6 = v.x * M.c2 - v.y * M.c1;
    ctype const c7 = v.x * M.c5 - v.y * M.c4 - v.z * M.c6;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& M, Scalar2dp<U> s)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s, MVec2dp<U> const& M)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product A * B for two multivectors from the even subalgebra (3d case)
// even-grade multivector * even-grade multivector = even-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// geometric product A * B for two multivectors from the odd subalgebra (3d case)
// odd-grade multivector * odd-grade multivector = even-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                        MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1;
    ctype const c1 = -A.c0 * B.c2 - A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype const c2 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0;
    ctype const c3 = A.c0 * B.c1 - A.c1 * B.c0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a multivector B of the odd subalgebra
// even-grade multivector * odd-grade multivector => odd-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype const c1 = A.c0 * B.c1 - A.c3 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c3;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c2;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product A * B of a multivector A from the odd subalgebra (3d case)
// with a multivector B of the even subalgebra
// odd-grade multivector * even-grade multivector => odd-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c3;
    ctype const c1 = A.c0 * B.c3 + A.c1 * B.c0;
    ctype const c2 = -A.c0 * B.c1 - A.c1 * B.c2 + A.c2 * B.c0 + A.c3 * B.c3;
    ctype const c3 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c3 + A.c3 * B.c0;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product A * ps of an even multivector A multiplied from the right
// by the trivector ps (=2dp pseudoscalar)
// even-grade multivector * trivector => odd multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = A.c3 * ctype(ps);
    ctype const c3 = A.c0 * ctype(ps);
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product ps * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the even-grade multivector B
// trivector * even-grade multivector => odd-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = ctype(ps) * B.c3;
    ctype const c3 = ctype(ps) * B.c0;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product A * B of an odd-grade multivector A multiplied from the right
// by the trivector B (=2dp pseudoscalar)
// odd-grade multivector * trivector => even-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                        PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c1 * ctype(ps);
    ctype const c2 = A.c0 * ctype(ps);
    ctype const c3 = 0.0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// geometric product A * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the odd-grade multivector B
// trivector * odd-grade multivector => even-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                        MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -ctype(ps) * B.c1;
    ctype const c2 = ctype(ps) * B.c0;
    ctype const c3 = 0.0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// geometric product M * B of a multivector M from the even subalgebra (3d case)
// with a bivector B from the right
// even-grade multivector * bivector => even-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& M,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c3 * B.z;
    ctype const c1 = M.c0 * B.x - M.c2 * B.z + M.c3 * B.y;
    ctype const c2 = M.c0 * B.y + M.c1 * B.z - M.c3 * B.x;
    ctype const c3 = M.c0 * B.z;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// geometric product A * M of a bivector A with an even-grade multivector M
// from the right
// bivector * even-grade multivector => returns an even-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                        MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.z * M.c3;
    ctype const c1 = B.x * M.c0 - B.y * M.c3 + B.z * M.c2;
    ctype const c2 = B.x * M.c3 + B.y * M.c0 - B.z * M.c1;
    ctype const c3 = B.z * M.c0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// geometric product M * B of an odd-grade multivector M with a bivector B
// odd-grade multivector * bivector => odd-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& M,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.z;
    ctype const c1 = M.c0 * B.z;
    ctype const c2 = -M.c0 * B.x - M.c1 * B.y + M.c3 * B.z;
    ctype const c3 = M.c0 * B.y - M.c1 * B.x - M.c2 * B.z;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product B * M of a bivector A with an odd-grade multivector B
// bivector * odd-grade multivector => odd-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                        MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.z * M.c1;
    ctype const c1 = -B.z * M.c0;
    ctype const c2 = B.x * M.c0 + B.y * M.c1 + B.z * M.c3;
    ctype const c3 = -B.x * M.c1 + B.y * M.c0 - B.z * M.c2;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product A * v of an even-grade multivector A with a vector v
// from the right
// even-grade multivector * vector => odd-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * v.x + A.c3 * v.y;
    ctype const c1 = A.c0 * v.y - A.c3 * v.x;
    ctype const c2 = A.c0 * v.z + A.c1 * v.x + A.c2 * v.y;
    ctype const c3 = -A.c1 * v.y + A.c2 * v.x - A.c3 * v.z;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product v * B of a vector v and an even-grade multivector B from the left
// vector * even-grade multivector => odd-grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * B.c0 - v.y * B.c3;
    ctype const c1 = v.x * B.c3 + v.y * B.c0;
    ctype const c2 = -v.x * B.c1 - v.y * B.c2 + v.z * B.c0;
    ctype const c3 = v.x * B.c2 - v.y * B.c1 - v.z * B.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& M,
                                                        Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y;
    ctype const c1 = -M.c0 * v.z + M.c2 * v.x - M.c3 * v.y;
    ctype const c2 = -M.c1 * v.z + M.c2 * v.y + M.c3 * v.x;
    ctype const c3 = M.c0 * v.y - M.c1 * v.x;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                        MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1;
    ctype const c1 = -v.x * M.c2 - v.y * M.c3 + v.z * M.c0;
    ctype const c2 = v.x * M.c3 - v.y * M.c2 + v.z * M.c1;
    ctype const c3 = v.x * M.c1 - v.y * M.c0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& M,
                                                        Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(Scalar2dp<T> s,
                                                        MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& M,
                                                        Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Scalar2dp<T> s,
                                                        MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// geometric product of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator*([[maybe_unused]] PScalar2dp<T>,
                                                        [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// geometric product ps * B of a trivector ps (=2dp pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                     BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = ctype(ps) * B.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// geometric product B * ps of bivector B multiplied by a trivector ps from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                     PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = B.z * ctype(ps);
    return Vec2dp<ctype>(c0, c1, c2);
}

// geometric product ps * v of a trivector ps (=2dp pseudoscalar) multiplied from the
// left to the vector v
// trivector * vector => bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                       Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.y;
    ctype const c1 = ctype(ps) * v.x;
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// geometric product v * ps of a vector v multiplied with a trivector ps from the right
// vector * trivector => bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                       PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(ps);
    ctype const c1 = v.x * ctype(ps);
    ctype const c2 = 0.0;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// geometric product ps * s of a trivector ps and a scalar s
// bivector * scalar => bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// geometric product s * ps of a scalar s and a trivector ps
// scalar * bivector => bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// geometric product A * B between two bivectors:
//
// the full geometric bivector product only exists in >= 4d spaces:
// A * B = gr0(A * B) + gr2(A * B) + gr4(A * B) = -dot(A,B) + cmt(A,B) + wdg(A,B)
// In 3D we don't have gr4(A * B) and thus only the terms up to grade 3 remain.
// The bivector product AxB = cmt(A,B) = 0.5*(A*B-B*A) is called the commutator product.
//
// A * B = -dot(A,B) + cmt(A,B)             (in 3d)
//
// => bivector*bivector = scalar + bivector = even-grade multivector (in 3d)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B1,
                                                        BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(-dot(B1, B2), cmt(B1, B2));
}

// geometric product B * v for a bivector B and a vector v
// B * v = (v << B) + wdg(B,v) = gr1(B * v) + gr3(B * v)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => bivector * vector = vector + trivector (= odd multivector)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                        Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>((v << B), wdg(B, v));
}

// geometric product v * B for a vector v and a bivector B
// v * B = (B >> v) + wdg(v,B) = gr1(v * B) + gr3(v * B)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => vector * bivector = vector + trivector (= odd multivector)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>((B >> v), wdg(v, B));
}

// (geometric) product B * s of bivector B multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                       Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(s);
    ctype const c1 = B.y * ctype(s);
    ctype const c2 = B.z * ctype(s);
    return BiVec2dp<ctype>(c0, c1, c2);
}

// (geometric) product s * B of a scalar s multiplied from the left to the bivector B
// scalar * bivector => bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s,
                                                       BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.x;
    ctype const c1 = ctype(s) * B.y;
    ctype const c2 = ctype(s) * B.z;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// geometric product a * b of two vectors
// a * b = dot(a,b) + wdg(a,b) = gr0(a * b) + gr2(a * b)
//
// HINT: if a full 3d multivector is required as result, it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// vector * vector => even-grade multivector (= scalar + bivector)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v1,
                                                        Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>((v2 >> v1), wdg(v1, v2));
}

// (geometric) product v * s of vector v multiplied with scalar s from the right
// vector * scalar => vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    return Vec2dp<ctype>(c0, c1, c2);
}

// (geometric) product s * v of scalar s multiplied from the left to vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
////////////////////////////////////////////////////////////////////////////////

// pga2dp rgpr :: rgpr(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c1 * B.c5 - A.c2 * B.c4 - A.c3 * B.c6 - A.c4 * B.c2 +
                     A.c5 * B.c1 - A.c6 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c0 * B.c5 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 - A.c4 * B.c6 +
                     A.c5 * B.c0 + A.c6 * B.c4 + A.c7 * B.c1;
    ctype const c2 = -A.c0 * B.c4 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 -
                     A.c4 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c2;
    ctype const c3 = A.c3 * B.c7 - A.c4 * B.c5 + A.c5 * B.c4 + A.c7 * B.c3;
    ctype const c4 = -A.c3 * B.c5 + A.c4 * B.c7 + A.c5 * B.c3 + A.c7 * B.c4;
    ctype const c5 = A.c3 * B.c4 - A.c4 * B.c3 + A.c5 * B.c7 + A.c7 * B.c5;
    ctype const c6 = A.c0 * B.c3 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c0 - A.c4 * B.c1 -
                     A.c5 * B.c2 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = -A.c3 * B.c3 + A.c4 * B.c4 + A.c5 * B.c5 + A.c7 * B.c7;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv,mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& A,
                                                 MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c0 * B.c2 - A.c4 * B.c3 + A.c5 * B.c0 + A.c6 * B.c1;
    ctype const c2 = -A.c0 * B.c1 - A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype const c3 = -A.c4 * B.c2 + A.c5 * B.c1;
    ctype const c4 = -A.c3 * B.c2 + A.c7 * B.c1;
    ctype const c5 = A.c3 * B.c1 + A.c7 * B.c2;
    ctype const c6 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c0 + A.c7 * B.c3;
    ctype const c7 = A.c4 * B.c1 + A.c5 * B.c2;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
                                                 MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c5 - A.c1 * B.c6 + A.c2 * B.c0 + A.c3 * B.c4;
    ctype const c2 = -A.c0 * B.c4 - A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype const c3 = -A.c1 * B.c5 + A.c2 * B.c4;
    ctype const c4 = A.c1 * B.c7 + A.c2 * B.c3;
    ctype const c5 = -A.c1 * B.c3 + A.c2 * B.c7;
    ctype const c6 = A.c0 * B.c3 - A.c1 * B.c1 - A.c2 * B.c2 + A.c3 * B.c7;
    ctype const c7 = A.c1 * B.c4 + A.c2 * B.c5;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& A,
                                                 MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 - A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c2;
    ctype const c1 = A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1 + A.c7 * B.c0;
    ctype const c2 = -A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0 + A.c7 * B.c1;
    ctype const c3 = A.c3 * B.c3 + A.c7 * B.c2;
    ctype const c4 = A.c4 * B.c3 + A.c5 * B.c2;
    ctype const c5 = -A.c4 * B.c2 + A.c5 * B.c3;
    ctype const c6 = A.c0 * B.c2 - A.c4 * B.c0 - A.c5 * B.c1 + A.c6 * B.c3;
    ctype const c7 = -A.c3 * B.c2 + A.c7 * B.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
                                                 MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c6 + A.c3 * B.c0;
    ctype const c1 = A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1;
    ctype const c2 = -A.c0 * B.c3 + A.c1 * B.c7 + A.c2 * B.c1 + A.c3 * B.c2;
    ctype const c3 = A.c2 * B.c7 + A.c3 * B.c3;
    ctype const c4 = -A.c2 * B.c5 + A.c3 * B.c4;
    ctype const c5 = A.c2 * B.c4 + A.c3 * B.c5;
    ctype const c6 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c0 + A.c3 * B.c6;
    ctype const c7 = -A.c2 * B.c3 + A.c3 * B.c7;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& M, PScalar2dp<U> ps)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(PScalar2dp<T> ps, MVec2dp<U> const& M)
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
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& M,
                                                 BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.y - M.c2 * B.x - M.c3 * B.z;
    ctype const c1 = M.c0 * B.y - M.c4 * B.z + M.c6 * B.x;
    ctype const c2 = -M.c0 * B.x - M.c5 * B.z + M.c6 * B.y;
    ctype const c3 = -M.c4 * B.y + M.c5 * B.x;
    ctype const c4 = -M.c3 * B.y + M.c7 * B.x;
    ctype const c5 = M.c3 * B.x + M.c7 * B.y;
    ctype const c6 = M.c1 * B.x + M.c2 * B.y + M.c7 * B.z;
    ctype const c7 = M.c4 * B.x + M.c5 * B.y;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B,
                                                 MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c2 + B.y * M.c1 - B.z * M.c3;
    ctype const c1 = -B.x * M.c6 + B.y * M.c0 + B.z * M.c4;
    ctype const c2 = -B.x * M.c0 - B.y * M.c6 + B.z * M.c5;
    ctype const c3 = -B.x * M.c5 + B.y * M.c4;
    ctype const c4 = B.x * M.c7 + B.y * M.c3;
    ctype const c5 = -B.x * M.c3 + B.y * M.c7;
    ctype const c6 = -B.x * M.c1 - B.y * M.c2 + B.z * M.c7;
    ctype const c7 = B.x * M.c4 + B.y * M.c5;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& M, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * v.y + M.c5 * v.x - M.c6 * v.z;
    ctype const c1 = M.c2 * v.z - M.c3 * v.y + M.c7 * v.x;
    ctype const c2 = -M.c1 * v.z + M.c3 * v.x + M.c7 * v.y;
    ctype const c3 = M.c7 * v.z;
    ctype const c4 = M.c5 * v.z;
    ctype const c5 = -M.c4 * v.z;
    ctype const c6 = M.c0 * v.z - M.c4 * v.x - M.c5 * v.y;
    ctype const c7 = -M.c3 * v.z;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v, MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c5 - v.y * M.c4 - v.z * M.c6;
    ctype const c1 = v.x * M.c7 + v.y * M.c3 - v.z * M.c2;
    ctype const c2 = -v.x * M.c3 + v.y * M.c7 + v.z * M.c1;
    ctype const c3 = v.z * M.c7;
    ctype const c4 = -v.z * M.c5;
    ctype const c5 = v.z * M.c4;
    ctype const c6 = v.x * M.c4 + v.y * M.c5 + v.z * M.c0;
    ctype const c7 = -v.z * M.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv,s) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c7 * ctype(s);
    ctype const c1 = M.c5 * ctype(s);
    ctype const c2 = -M.c4 * ctype(s);
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = M.c3 * ctype(s);
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(s,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(Scalar2dp<T> s, MVec2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c7;
    ctype const c1 = ctype(s) * M.c5;
    ctype const c2 = -ctype(s) * M.c4;
    ctype const c3 = 0.0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = ctype(s) * M.c3;
    ctype const c7 = 0.0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// pga2dp rgpr :: rgpr(mv_e,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
                                                   MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c2 - A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1;
    ctype const c1 = -A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c2 = -A.c1 * B.c2 + A.c2 * B.c1;
    ctype const c3 = A.c1 * B.c1 + A.c2 * B.c2;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(mv_e,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
                                                   MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c2;
    ctype const c1 = A.c1 * B.c3 + A.c2 * B.c2;
    ctype const c2 = -A.c1 * B.c2 + A.c2 * B.c3;
    ctype const c3 = A.c0 * B.c2 - A.c1 * B.c0 - A.c2 * B.c1 + A.c3 * B.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(mv_u,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
                                                   MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c3 + A.c3 * B.c0;
    ctype const c1 = -A.c2 * B.c2 + A.c3 * B.c1;
    ctype const c2 = A.c2 * B.c1 + A.c3 * B.c2;
    ctype const c3 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c0 + A.c3 * B.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
                                                   PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * ctype(ps);
    ctype const c1 = A.c1 * ctype(ps);
    ctype const c2 = A.c2 * ctype(ps);
    ctype const c3 = A.c3 * ctype(ps);
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(PScalar2dp<T> ps,
                                                   MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.c0;
    ctype const c1 = ctype(ps) * B.c1;
    ctype const c2 = ctype(ps) * B.c2;
    ctype const c3 = ctype(ps) * B.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(mv_e,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& M,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * B.y - M.c1 * B.z + M.c3 * B.x;
    ctype const c1 = -M.c0 * B.x - M.c2 * B.z + M.c3 * B.y;
    ctype const c2 = -M.c1 * B.y + M.c2 * B.x;
    ctype const c3 = M.c1 * B.x + M.c2 * B.y;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(bivec,mv_e) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B,
                                                   MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c3 + B.y * M.c0 + B.z * M.c1;
    ctype const c1 = -B.x * M.c0 - B.y * M.c3 + B.z * M.c2;
    ctype const c2 = -B.x * M.c2 + B.y * M.c1;
    ctype const c3 = B.x * M.c1 + B.y * M.c2;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(mv_e,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& M,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * v.y + M.c2 * v.x - M.c3 * v.z;
    ctype const c1 = M.c2 * v.z;
    ctype const c2 = -M.c1 * v.z;
    ctype const c3 = M.c0 * v.z - M.c1 * v.x - M.c2 * v.y;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(vec,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v,
                                                   MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c2 - v.y * M.c1 - v.z * M.c3;
    ctype const c1 = -v.z * M.c2;
    ctype const c2 = v.z * M.c1;
    ctype const c3 = v.x * M.c1 + v.y * M.c2 + v.z * M.c0;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(mv_e,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * ctype(s);
    ctype const c1 = -M.c1 * ctype(s);
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(s,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rgpr(Scalar2dp<T> s, MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c2;
    ctype const c1 = -ctype(s) * M.c1;
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(mv_u,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
                                                   MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    ctype const c1 = -A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1;
    ctype const c2 = A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c3 = -A.c2 * B.c2 + A.c3 * B.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
                                                   PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * ctype(ps);
    ctype const c1 = A.c1 * ctype(ps);
    ctype const c2 = A.c2 * ctype(ps);
    ctype const c3 = A.c3 * ctype(ps);
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(PScalar2dp<T> ps,
                                                   MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.c0;
    ctype const c1 = ctype(ps) * B.c1;
    ctype const c2 = ctype(ps) * B.c2;
    ctype const c3 = ctype(ps) * B.c3;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(mv_u,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& M,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * B.y - M.c1 * B.x - M.c2 * B.z;
    ctype const c1 = -M.c2 * B.y + M.c3 * B.x;
    ctype const c2 = M.c2 * B.x + M.c3 * B.y;
    ctype const c3 = M.c0 * B.x + M.c1 * B.y + M.c3 * B.z;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(bivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B,
                                                   MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c1 + B.y * M.c0 - B.z * M.c2;
    ctype const c1 = B.x * M.c3 + B.y * M.c2;
    ctype const c2 = -B.x * M.c2 + B.y * M.c3;
    ctype const c3 = -B.x * M.c0 - B.y * M.c1 + B.z * M.c3;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(mv_u,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& M,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.z - M.c2 * v.y + M.c3 * v.x;
    ctype const c1 = -M.c0 * v.z + M.c2 * v.x + M.c3 * v.y;
    ctype const c2 = M.c3 * v.z;
    ctype const c3 = -M.c2 * v.z;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(vec,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v,
                                                   MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c3 + v.y * M.c2 - v.z * M.c1;
    ctype const c1 = -v.x * M.c2 + v.y * M.c3 + v.z * M.c0;
    ctype const c2 = v.z * M.c3;
    ctype const c3 = -v.z * M.c2;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
}

// pga2dp rgpr :: rgpr(mv_u,s) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c3 * ctype(s);
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = M.c2 * ctype(s);
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(s,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(Scalar2dp<T> s, MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c3;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = ctype(s) * M.c2;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
}

// pga2dp rgpr :: rgpr(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rgpr(PScalar2dp<T> ps1, PScalar2dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps1) * ctype(ps2));
}

// pga2dp rgpr :: rgpr(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rgpr(PScalar2dp<T> ps, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.x;
    ctype const c1 = ctype(ps) * B.y;
    ctype const c2 = ctype(ps) * B.z;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(ps);
    ctype const c1 = B.y * ctype(ps);
    ctype const c2 = B.z * ctype(ps);
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rgpr(PScalar2dp<T> ps, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rgpr(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(ps) * ctype(s));
}

// pga2dp rgpr :: rgpr(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rgpr(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s) * ctype(ps));
}

// rgpr(B1,B2) = rwdg(B1, B2) + rdot(B1, B2)
// pga2dp rgpr :: rgpr(bivec,bivec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B1,
                                                   BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -B1.x * B2.z + B1.z * B2.x;
    // ctype const c1 = -B1.y * B2.z + B1.z * B2.y;
    // ctype const c2 = -B1.x * B2.y + B1.y * B2.x;
    // ctype const c3 = B1.x * B2.x + B1.y * B2.y;
    // return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
    return MVec2dp_U<ctype>(rwdg(B1, B2), rdot(B1, B2));
}

// rgpr(B,v) = rwdg(B, v) + rcmt(B, v)
// pga2dp rgpr :: rgpr(bivec,vec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -B.x * v.y + B.y * v.x - B.z * v.z;
    // ctype const c1 = B.y * v.z;
    // ctype const c2 = -B.x * v.z;
    // ctype const c3 = -B.x * v.x - B.y * v.y;
    // return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
    return MVec2dp_E<ctype>(rwdg(B, v), rcmt(B, v));
}

// rgpr(v,B) = rwdg(v, B) + rcmt(v, B)
// pga2dp rgpr :: rgpr(vec,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = v.x * B.y - v.y * B.x - v.z * B.z;
    // ctype const c1 = -v.z * B.y;
    // ctype const c2 = v.z * B.x;
    // ctype const c3 = v.x * B.x + v.y * B.y;
    // return MVec2dp_E<ctype>(Scalar2dp<ctype>(c0), BiVec2dp<ctype>(c1, c2, c3));
    return MVec2dp_E<ctype>(rwdg(v, B), rcmt(v, B));
}

// pga2dp rgpr :: rgpr(bivec,s) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * ctype(s);
    ctype const c1 = -B.x * ctype(s);
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(s,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rgpr(Scalar2dp<T> s, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.y;
    ctype const c1 = -ctype(s) * B.x;
    ctype const c2 = 0.0;
    return Vec2dp<ctype>(c0, c1, c2);
}

// rgpr(v1,v2) = rcmt(v1,v2) - rdot(v1,v2)
// pga2dp rgpr :: rgpr(vec,vec) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v1,
                                                   Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = v1.y * v2.z - v1.z * v2.y;
    // ctype const c1 = -v1.x * v2.z + v1.z * v2.x;
    // ctype const c2 = 0.0;
    // ctype const c3 = -v1.z * v2.z;
    // return MVec2dp_U<ctype>(Vec2dp<ctype>(c0, c1, c2), PScalar2dp<ctype>(c3));
    return MVec2dp_U<ctype>(rcmt(v1, v2), -rdot(v1, v2));
}

// pga2dp rgpr :: rgpr(vec,s) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = v.z * ctype(s);
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(s,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rgpr(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = ctype(s) * v.z;
    return BiVec2dp<ctype>(c0, c1, c2);
}

// pga2dp rgpr :: rgpr(s,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rgpr([[maybe_unused]] Scalar2dp<T>,
                                                    [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
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
inline Scalar2dp<T> inv(Scalar2dp<T> s)
{
    T sq_n = bulk_nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar2dp<T>(rev(s) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec2dp<T> inv(Vec2dp<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v) = v/bulk_sq_nrm(v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    T sq_n = bulk_nrm_sq(v);
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return Vec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec2dp<T> inv(BiVec2dp<T> const& B)
{
    // B^(-1) = rev(B)/|B|^2
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    T sq_n = bulk_nrm_sq(B);
    hd::ga::detail::check_normalization<T>(sq_n, "bivector");
    T inv = -T(1.0) / sq_n; // negative inverse of squared norm for a bivector
    return BiVec2dp<T>(B.x * inv, B.y * inv, B.z * inv);
}

// due to the degenerate metric the pseudoscalar does not have an inverse

template <typename T>
    requires(numeric_type<T>)
inline MVec2dp_E<T> inv(MVec2dp_E<T> const& E)
{
    T sq_n = bulk_nrm_sq(E);
    hd::ga::detail::check_normalization<T>(sq_n, "even-grade multivector");
    T inv = T(1.0) / sq_n;
    return MVec2dp_E<T>(rev(E) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2dp_U<T> inv(MVec2dp_U<T> const& U)
{
    T sq_n = bulk_nrm_sq(U);
    hd::ga::detail::check_normalization<T>(sq_n, "odd-grade multivector");
    T inv = T(1.0) / sq_n;
    return MVec2dp_U<T>(rev(U) * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016 left and a right inverse are the same (see paper
// of Hitzer, Sangwine)
template <typename T>
    requires(numeric_type<T>)
inline MVec2dp<T> inv(MVec2dp<T> const& M)
{
    T m_conjm = T(gr0(M * conj(M) * gr_inv(M) * rev(M)));
    hd::ga::detail::check_normalization<T>(std::abs(m_conjm), "multivector");
    T inv = T(1.0) / m_conjm;
    return MVec2dp<T>(conj(M) * gr_inv(M) * rev(M) * inv);
}


/////////////////////////////////////////////////////////////////////////////////////////
// Projective contractions for 2dp:
//
// REQUIRES: the dualized operand must be of smaller grade, or the result is zero
//
// l_bulk_contract2dp(a,b) = rwdg(bulk_dual(a), b)
// l_weight_contract2dp(a,b) = rwdg(weight_dual(a), b)
//
// r_bulk_contract2dp(a,b) = rwdg(a, bulk_dual(b))
// r_weight_contract2dp(a,b) = rwdg(a, weight_dual(b))
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
decltype(auto) l_bulk_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(bulk_dual(a), b);
    return rwdg(bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) l_weight_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(weight_dual(a), b);
    return rwdg(weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) r_bulk_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, bulk_dual(b));
    return rwdg(std::forward<arg1>(a), bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) r_weight_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, weight_dual(b));
    return rwdg(std::forward<arg1>(a), weight_dual(std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// Projective expansions for 2dp:
//
// REQUIRES: the dualized operand must be of larger grade, or the result is zero
//
// l_bulk_expand2dp(a,b) = wdg(bulk_dual(a), b)
// (dual to l_weight_contract)
// l_weight_expand2dp(a,b) = wdg(weight_dual(a), b)
// (dual to l_bulk_contract)
//
// r_bulk_expand2dp(a,b) = wdg(a, bulk_dual(b))
// (dual to r_weight_contract)
// r_weight_expand2dp(a,b) = wdg(a, weight_dual(b))
// (dual to r_bulk_contract)
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
decltype(auto) l_bulk_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(bulk_dual(a), b);
    return wdg(bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) l_weight_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(weight_dual(a), b);
    return wdg(weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) r_bulk_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(a, bulk_dual(b));
    return wdg(std::forward<arg1>(a), bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) r_weight_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(a, weight_dual(b));
    return wdg(std::forward<arg1>(a), weight_dual(std::forward<arg2>(b)));
}

} // namespace hd::ga::pga