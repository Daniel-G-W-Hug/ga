#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_ega2d_ops_basics.hpp" // ega2d ops basics


namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega2d product operations:
//
// - dot()              -> dot product
// - wdg()              -> wedge product
// - rwdg()             -> regressive wedge product
// - twdg1()            -> transwedge product (k=1)
// - rtwdg1()           -> regressive transwedge product (k=1)
// - operator<<()       -> left contraction
// - operator>>()       -> right contraction
// - cmt()              -> commutator product
// - operator*()        -> geometric product
// - inv()              -> inversion operation
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot product (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1
//
// Via the exomorphism G(a ^ b) = G(a) ^ G(b) the extended metric G is defined as a 2^n by
// 2^n matrix via the outer product, which in turn allows for definition of an extended
// dot product for all grades of a full multi-vector from scalar to pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

// ega2d dot :: dot(mv,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3);
}

// ega2d dot :: dot(mv_e,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(MVec2d_E<T> const& A,
                                                 MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(A.c0 * B.c0 + A.c1 * B.c1);
}

// ega2d dot :: dot(ps,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(PScalar2d<T> ps1, PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// return dot-product of two vectors in G<2,0,0>
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
// ega2d dot :: dot(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// ega2d dot :: dot(s,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// wedge product (= outer product)
////////////////////////////////////////////////////////////////////////////////

// wedge product extended to fully populated multivectors
// ega2d wdg :: mv ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d wdg :: mv ^ mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& A, MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d wdg :: mv_e ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d wdg :: mv ^ ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(M.c0 * ctype(ps));
}

// ega2d wdg :: ps ^ mv -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * M.c0);
}

// ega2d wdg :: mv ^ vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c1 * v.y - M.c2 * v.x;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d wdg :: vec ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.x * M.c2 - v.y * M.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d wdg :: mv ^ s -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d wdg :: s ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d wdg :: mv_e ^ mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& A,
                                                 MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d wdg :: mv_e ^ ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(M.c0 * ctype(ps));
}

// ega2d wdg :: ps ^ mv_e -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * M.c0);
}

// ega2d wdg :: mv_e ^ vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    return Vec2d<ctype>(c0, c1);
}

// ega2d wdg :: vec ^ mv_e -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    return Vec2d<ctype>(c0, c1);
}

// ega2d wdg :: mv_e ^ s -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d wdg :: s ^ mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(Scalar2d<T> s, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d wdg :: ps ^ ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2d<T>,
                                                 [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d wdg :: ps ^ vec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2d<T>,
                                                 [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d wdg :: vec ^ ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec2d<T> const&,
                                                 [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d wdg :: ps ^ s -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// ega2d wdg :: s ^ ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// ega2d wdg :: vec ^ vec -> ps
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(v1.x * v2.y - v1.y * v2.x);
}

// ega2d wdg :: vec ^ s -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    return Vec2d<ctype>(c0, c1);
}

// ega2d wdg :: s ^ vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    return Vec2d<ctype>(c0, c1);
}

// ega2d wdg :: s ^ s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product between complements)
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product (as is the case for the dual), but just depending
// on the outer product (wdg) as well as the complement of the multivector
// (i.e. this definition is does NOT connect directly to the geometric product, but to the
// outer product, exclusively)
//
// rwdg(ul, ur) = l_cmpl(wdg(r_cmpl(ul),r_cmpl(ur))) = l_cmpl(r_cmpl(ul)^r_cmpl(ur))
//              = r_cmpl(wdg(l_cmpl(ul),l_cmpl(ur))) = r_cmpl(l_cmpl(ul)^l_cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// regressive wedge product extended to fully populated multivectors
// ega2d rwdg :: rwdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    ctype const c1 = A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c2 = A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c3 = A.c3 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d rwdg :: rwdg(mv, mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(MVec2d<T> const& A, MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c3 * B.c0;
    ctype const c1 = A.c1 * B.c1;
    ctype const c2 = A.c2 * B.c1;
    ctype const c3 = A.c3 * B.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d rwdg :: rwdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(MVec2d_E<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 + A.c1 * B.c0;
    ctype const c1 = A.c1 * B.c1;
    ctype const c2 = A.c1 * B.c2;
    ctype const c3 = A.c1 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d rwdg :: rwdg(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(MVec2d<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d rwdg :: rwdg(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d rwdg :: rwdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(MVec2d<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.y - M.c2 * v.x;
    ctype const c1 = M.c3 * v.x;
    ctype const c2 = M.c3 * v.y;
    ctype const c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d rwdg :: rwdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c2 - v.y * M.c1;
    ctype const c1 = v.x * M.c3;
    ctype const c2 = v.y * M.c3;
    ctype const c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d rwdg :: rwdg(mv,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c3 * ctype(s));
}

// ega2d rwdg :: rwdg(s,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s) * M.c3);
}

// ega2d rwdg :: rwdg(mv_e,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> rwdg(MVec2d_E<T> const& A,
                                                  MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c1 = A.c1 * B.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d rwdg :: rwdg(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> rwdg(MVec2d_E<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d rwdg :: rwdg(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d rwdg :: rwdg(mv_e,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> rwdg(MVec2d_E<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x;
    ctype const c1 = M.c1 * v.y;
    return Vec2d<ctype>(c0, c1);
}

// ega2d rwdg :: rwdg(vec,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1;
    ctype const c1 = v.y * M.c1;
    return Vec2d<ctype>(c0, c1);
}

// ega2d rwdg :: rwdg(mv_e,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(MVec2d_E<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c1 * ctype(s));
}

// ega2d rwdg :: rwdg(s,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Scalar2d<T> s, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s) * M.c1);
}

// ega2d rwdg :: rwdg(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps1, PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// ega2d rwdg :: rwdg(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    return Vec2d<ctype>(c0, c1);
}

// ega2d rwdg :: rwdg(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    return Vec2d<ctype>(c0, c1);
}

// ega2d rwdg :: rwdg(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps) * ctype(s));
}

// ega2d rwdg :: rwdg(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s) * ctype(ps));
}

// ega2d rwdg :: rwdg(vec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.y - v1.y * v2.x);
}

// ega2d rwdg :: rwdg(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2d<T> const&,
                                                  [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d rwdg :: rwdg(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2d<T>,
                                                  [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d rwdg :: rwdg(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2d<T>,
                                                  [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// twdg1(ps,vec) = vec -> identical to geometric product gpr(ps,vec)
//                     -> and identical to the left contraction vec << ps
// (turn by -90°, i.e. against orientation of e12)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> twdg1(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.y;
    ctype const c1 = -ctype(ps) * v.x;
    return Vec2d<ctype>(c0, c1);
}

// twdg1(vec,ps) = vec -> identical to geometric product gpr(vec,ps)
//                     -> and identical to the right contraction ps >> vec
// (turn by +90°, i.e. into orientation of e12)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> twdg1(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(ps);
    ctype const c1 = v.x * ctype(ps);
    return Vec2d<ctype>(c0, c1);
}

// twdg1(vec,vec) = scalar -> identical to dot product dot(vec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> twdg1(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}


////////////////////////////////////////////////////////////////////////////////
// regressive transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// rtwdg1(vec,vec) = pseudoscalar -> identical to regressive dot product rdot(vec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> rtwdg1(Vec2d<T> const& v1,
                                                     Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// rtwdg1(vec,s) = vec -> returns a vector perpendicular to vec scaled by s
// identical to right expansions: -l_expand(vec,s) = r_expand(s,vec)
// => sign change makes sense since expansions are derived from the wedge product
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> rtwdg1(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(s);
    ctype const c1 = v.x * ctype(s);
    return Vec2d<ctype>(c0, c1);
}

// rtwdg1(s,vec) -> returns a vector perpendicular to vec scaled by s
// identical to expansions: -r_expand(s,vec) = l_expand(vec,s)
// => sign change makes sense since expansions are derived from the wedge product
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> rtwdg1(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.y;
    ctype const c1 = -ctype(s) * v.x;
    return Vec2d<ctype>(c0, c1);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B
//
// The resulting object is A lies in B and is perpendicular to A
//
// L. Dorst: The contraction A << B of an a-blade A onto a b-blade B is a sub-blade
// of B of grade b-a which is perpendicular to A, and linear in both arguments
////////////////////////////////////////////////////////////////////////////////

// ega2d l_contract :: mv << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c3;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c3;
    ctype const c3 = A.c0 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d l_contract :: mv << mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& A,
                                                      MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype const c1 = A.c2 * B.c1;
    ctype const c2 = -A.c1 * B.c1;
    ctype const c3 = A.c0 * B.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d l_contract :: mv_e << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c3;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d l_contract :: mv << ps -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c3 * ctype(ps);
    ctype const c1 = M.c2 * ctype(ps);
    ctype const c2 = -M.c1 * ctype(ps);
    ctype const c3 = M.c0 * ctype(ps);
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d l_contract :: ps << mv -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps,
                                                        MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps) * M.c3);
}

// ega2d l_contract :: mv << vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& M,
                                                      Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d l_contract :: vec << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v,
                                                      MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2;
    ctype const c1 = v.y * M.c3;
    ctype const c2 = -v.x * M.c3;
    ctype const c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d l_contract :: mv << s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c0 * ctype(s));
}

// ega2d l_contract :: s << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d l_contract :: mv_e << mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& A,
                                                        MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1;
    ctype const c1 = A.c0 * B.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d l_contract :: mv_e << ps -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                        PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * ctype(ps);
    ctype const c1 = M.c0 * ctype(ps);
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d l_contract :: ps << mv_e -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps) * M.c1);
}

// ega2d l_contract :: mv_e << vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    return Vec2d<ctype>(c0, c1);
}

// ega2d l_contract :: vec << mv_e -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v,
                                                     MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * M.c1;
    ctype const c1 = -v.x * M.c1;
    return Vec2d<ctype>(c0, c1);
}

// ega2d l_contract :: mv_e << s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                        Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c0 * ctype(s));
}

// ega2d l_contract :: s << mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(Scalar2d<T> s,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d l_contract :: ps << ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps1,
                                                        PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// ega2d l_contract :: ps << vec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2d<T>,
                                                        [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (v << ps) = v * ps = cmt(ps, v)
// ega2d l_contract :: vec << ps -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * ctype(ps);
    ctype const c1 = -v.x * ctype(ps);
    return Vec2d<ctype>(c0, c1);
}

// ega2d l_contract :: ps << s -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2d<T>,
                                                        [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d l_contract :: s << ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// left contraction (v1 << v2) = dot(v1,v2) = (v1 >> v2))
// ega2d l_contract :: vec << vec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v1,
                                                        Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// ega2d l_contract :: vec << s -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<([[maybe_unused]] Vec2d<T> const&,
                                                        [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d l_contract :: s << vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    return Vec2d<ctype>(c0, c1);
}

// ega2d l_contract :: s << s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B
//
// The resulting object lies in A and is perpendicular to B
////////////////////////////////////////////////////////////////////////////////

// ega2d r_contract :: mv >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype const c1 = A.c1 * B.c0 - A.c3 * B.c2;
    ctype const c2 = A.c2 * B.c0 + A.c3 * B.c1;
    ctype const c3 = A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d r_contract :: mv >> mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& A,
                                                      MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d r_contract :: mv_e >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c3;
    ctype const c1 = -A.c1 * B.c2;
    ctype const c2 = A.c1 * B.c1;
    ctype const c3 = A.c1 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d r_contract :: mv >> ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& M,
                                                        PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c3 * ctype(ps));
}

// ega2d r_contract :: ps >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c3;
    ctype const c1 = -ctype(ps) * M.c2;
    ctype const c2 = ctype(ps) * M.c1;
    ctype const c3 = ctype(ps) * M.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d r_contract :: mv >> vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& M,
                                                      Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y;
    ctype const c1 = -M.c3 * v.y;
    ctype const c2 = M.c3 * v.x;
    ctype const c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d r_contract :: vec >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v,
                                                      MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d r_contract :: mv >> s -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d r_contract :: s >> mv -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s) * M.c0);
}

// ega2d r_contract :: mv_e >> mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& A,
                                                        MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1;
    ctype const c1 = A.c1 * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d r_contract :: mv_e >> ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                        PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c1 * ctype(ps));
}

// ega2d r_contract :: ps >> mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c1;
    ctype const c1 = ctype(ps) * M.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d r_contract :: mv_e >> vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * v.y;
    ctype const c1 = M.c1 * v.x;
    return Vec2d<ctype>(c0, c1);
}

// ega2d r_contract :: vec >> mv_e -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v,
                                                     MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    return Vec2d<ctype>(c0, c1);
}

// ega2d r_contract :: mv_e >> s -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                        Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d r_contract :: s >> mv_e -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s) * M.c0);
}

// ega2d r_contract :: ps >> ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps1,
                                                        PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// right contraction (ps >> v) = ps * v = cmt(v,ps)
// ega2d r_contract :: ps >> vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * v.y;
    ctype const c1 = ctype(ps) * v.x;
    return Vec2d<ctype>(c0, c1);
}

// ega2d r_contract :: vec >> ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Vec2d<T> const&,
                                                        [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d r_contract :: ps >> s -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// ega2d r_contract :: s >> ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar2d<T>,
                                                        [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// right contraction (v1 >> v2) = dot(v1,v2) = (v1 << v2)
// ega2d r_contract :: vec >> vec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v1,
                                                        Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// ega2d r_contract :: vec >> s -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    return Vec2d<ctype>(c0, c1);
}

// ega2d r_contract :: s >> vec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar2d<T>,
                                                        [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d r_contract :: s >> s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

// ega2d cmt :: cmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> cmt(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c2 = A.c1 * B.c3 - A.c3 * B.c1;
    ctype const c3 = A.c1 * B.c2 - A.c2 * B.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d cmt :: cmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2d<T>,
                                                 [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// cmt(ps,v) == -cmt(v,ps) == (v << ps)
// ega2d cmt :: cmt(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> cmt(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.y;
    ctype const c1 = -ctype(ps) * v.x;
    return Vec2d<ctype>(c0, c1);
}

// cmt(v,ps) == -cmt(ps,v) == (ps >> v)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> cmt(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(ps);
    ctype const c1 = v.x * ctype(ps);
    return Vec2d<ctype>(c0, c1);
}

// ega2d cmt :: cmt(ps,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar2d<T>,
                                                  [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(0.0);
}

// ega2d cmt :: cmt(s,ps) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2d<T>,
                                                  [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(0.0);
}

// ega2d cmt :: cmt(vec,vec) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> cmt(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(v1.x * v2.y - v1.y * v2.x);
}

// ega2d cmt :: cmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec2d<T> const&,
                                                 [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d cmt :: cmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2d<T>,
                                                 [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// ega2d cmt :: cmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar2d<T>,
                                                 [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// ega2d gpr :: mv * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                     MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: mv * mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                     MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c3 * B.c1;
    ctype const c1 = A.c1 * B.c0 - A.c2 * B.c1;
    ctype const c2 = A.c1 * B.c1 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: mv_e * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                     MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c3;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c2;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: mv * ps -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c3 * ctype(ps);
    ctype const c1 = -A.c2 * ctype(ps);
    ctype const c2 = A.c1 * ctype(ps);
    ctype const c3 = A.c0 * ctype(ps);
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: ps * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c3;
    ctype const c1 = ctype(ps) * B.c2;
    ctype const c2 = -ctype(ps) * B.c1;
    ctype const c3 = ctype(ps) * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: mv * vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * v.x + A.c2 * v.y;
    ctype const c1 = A.c0 * v.x + A.c3 * v.y;
    ctype const c2 = A.c0 * v.y - A.c3 * v.x;
    ctype const c3 = A.c1 * v.y - A.c2 * v.x;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: vec * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v,
                                                     MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * B.c1 + v.y * B.c2;
    ctype const c1 = v.x * B.c0 - v.y * B.c3;
    ctype const c2 = v.x * B.c3 + v.y * B.c0;
    ctype const c3 = v.x * B.c2 - v.y * B.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: mv * s -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * ctype(s);
    ctype const c1 = A.c1 * ctype(s);
    ctype const c2 = A.c2 * ctype(s);
    ctype const c3 = A.c3 * ctype(s);
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: s * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.c0;
    ctype const c1 = ctype(s) * B.c1;
    ctype const c2 = ctype(s) * B.c2;
    ctype const c3 = ctype(s) * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// ega2d gpr :: mv_e * mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                       MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c1;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d gpr :: mv_e * ps -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                       PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -A.c1 * ctype(ps);
    ctype const c1 = A.c0 * ctype(ps);
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d gpr :: ps * mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(PScalar2d<T> ps,
                                                       MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.c1;
    ctype const c1 = ctype(ps) * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d gpr :: mv_e * vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                    Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * v.x + A.c1 * v.y;
    ctype const c1 = A.c0 * v.y - A.c1 * v.x;
    return Vec2d<ctype>(c0, c1);
}

// ega2d gpr :: vec * mv_e -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v,
                                                    MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * B.c0 - v.y * B.c1;
    ctype const c1 = v.x * B.c1 + v.y * B.c0;
    return Vec2d<ctype>(c0, c1);
}

// ega2d gpr :: mv_e * s -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                       Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * ctype(s);
    ctype const c1 = A.c1 * ctype(s);
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d gpr :: s * mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Scalar2d<T> s,
                                                       MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.c0;
    ctype const c1 = ctype(s) * B.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// ega2d gpr :: ps * ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps1, PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(-ctype(ps1) * ctype(ps2));
}

// ega2d gpr :: ps * vec -> vec
// this multiplication turns the vector by -90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.y;
    ctype const c1 = -ctype(ps) * v.x;
    return Vec2d<ctype>(c0, c1);
}

// ega2d gpr :: vec * ps -> vec
// this multiplication turns the vector by +90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * ctype(ps);
    ctype const c1 = v.x * ctype(ps);
    return Vec2d<ctype>(c0, c1);
}

// ega2d gpr :: ps * s -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// ega2d gpr :: s * ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// geometric product ab for two vectors (returns a multivector of the even subalgebra)
// v1*v2 = dot(v1,v2) + wdg(v1,v2)
// => vec * vec = s + bivec = mv_e
//
// HINT: if a full 2d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types

// ega2d gpr :: vec * vec -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Vec2d<T> const& v1,
                                                       Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = v1.x * v2.x + v1.y * v2.y;
    // ctype const c1 = v1.x * v2.y - v1.y * v2.x;
    // return MVec2d_E<ctype>(c0, c1);
    return MVec2d_E<ctype>(dot(v1, v2), wdg(v1, v2));
}

// ega2d gpr :: vec * s -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    return Vec2d<ctype>(c0, c1);
}

// ega2d gpr :: s * vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    return Vec2d<ctype>(c0, c1);
}

// ega2d gpr :: s * s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors
// w.r.t. the geometric product:
// for k-blades: A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
////////////////////////////////////////////////////////////////////////////////
// HINT: inv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Scalar2d<T> inv(Scalar2d<T> s)
{
    T sq_n = nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar2d<T>(rev(s) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec2d<T> inv(Vec2d<T> const& v)
{
    T sq_n = nrm_sq(v);
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n;
    return Vec2d<T>(rev(v) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline PScalar2d<T> inv(PScalar2d<T> ps)
{
    T sq_n = nrm_sq(ps);
    hd::ga::detail::check_normalization<T>(sq_n, "pseudoscalar");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector

    return PScalar2d<T>(rev(ps) * inv);
}

// return the multiplicative inverse of the complex number (inv(z) = 1/nrm_sq(z)*rev(z))
// with rev(z) being the complex conjugate
template <typename T>
    requires(numeric_type<T>)
inline MVec2d_E<T> inv(MVec2d_E<T> const& E)
{
    T sq_n = nrm_sq(E);
    hd::ga::detail::check_normalization<T>(sq_n, "even-grade multivector");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return MVec2d_E<T>(rev(E) * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
// left and a right inverse are the same (see paper of Hitzer, Sangwine)
template <typename T>
    requires(numeric_type<T>)
inline MVec2d<T> inv(MVec2d<T> const& M)
{
    T m_conjm = T(gr0(M * conj(M)));
    hd::ga::detail::check_normalization<T>(std::abs(m_conjm), "multivector");
    T inv = T(1.0) / m_conjm; // inverse of squared norm for a vector
    return MVec2d<T>(conj(M) * inv);
}

} // namespace hd::ga::ega