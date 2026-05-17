#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_ega3d_ops_basics.hpp" // ega3d ops basics


namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega3d product operations:
//
// - dot()              -> dot product
// - wdg()              -> wedge product
// - rwdg()             -> regressive wedge product
// - twdg1()            -> transwedge product (k=1)
// - rtwdg1()           -> regressive transwedge product (k=1)
// - operator<<()       -> left contraction
// - operator>>()       -> right contraction
// - cross()            -> cross product
// - cmt()              -> commutator product
// - operator*()        -> geometric product
// - inv()              -> inversion operation
//
// - l_expand()         -> left expansion
// - r_expand()         -> right expansion
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot product (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = +1
//
// Via the exomorphism G(a ^ b) = G(a) ^ G(b) the extended metric G is defined as a 2^n by
// 2^n matrix via the outer product, which in turn allows for definition of an extended
// dot product for all grades of a full multi-vector from scalar to pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> dot(MVec3d<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 +
                           A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> dot(MVec3d_E<T> const& A,
                                                 MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> dot(MVec3d_U<T> const& A,
                                                 MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> dot(PScalar3d<T> ps1, PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> dot(BiVec3d<T> const& B1,
                                                 BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B1.x * B2.x + B1.y * B2.y + B1.z * B2.z);
}

// coordinate free definition: dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> dot(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> dot(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// wedge product (= outer product)
////////////////////////////////////////////////////////////////////////////////

// wedge product extended to fully populated multivectors

// ega3d wdg :: mv ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c4 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c3 * B.c1 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
                     A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv ^ mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& A, MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    ctype const c4 = A.c0 * B.c1 + A.c4 * B.c0;
    ctype const c5 = A.c0 * B.c2 + A.c5 * B.c0;
    ctype const c6 = A.c0 * B.c3 + A.c6 * B.c0;
    ctype const c7 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv_e ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d_E<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c0;
    ctype const c5 = A.c0 * B.c5 + A.c2 * B.c0;
    ctype const c6 = A.c0 * B.c6 + A.c3 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv ^ mv_u -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& A, MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c0 * B.c1;
    ctype const c3 = A.c0 * B.c2;
    ctype const c4 = A.c2 * B.c2 - A.c3 * B.c1;
    ctype const c5 = -A.c1 * B.c2 + A.c3 * B.c0;
    ctype const c6 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype const c7 = A.c0 * B.c3 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv_u ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d_U<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = A.c0 * B.c0;
    ctype const c2 = A.c1 * B.c0;
    ctype const c3 = A.c2 * B.c0;
    ctype const c4 = A.c1 * B.c3 - A.c2 * B.c2;
    ctype const c5 = -A.c0 * B.c3 + A.c2 * B.c1;
    ctype const c6 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv ^ ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& M, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(M.c0 * ctype(ps));
}

// ega3d wdg :: ps ^ mv -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(PScalar3d<T> ps, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * M.c0);
}

// ega3d wdg :: mv ^ bivec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& M, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = M.c0 * B.x;
    ctype const c5 = M.c0 * B.y;
    ctype const c6 = M.c0 * B.z;
    ctype const c7 = M.c1 * B.x + M.c2 * B.y + M.c3 * B.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: bivec ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = 0.0;
    ctype const c2 = 0.0;
    ctype const c3 = 0.0;
    ctype const c4 = B.x * M.c0;
    ctype const c5 = B.y * M.c0;
    ctype const c6 = B.z * M.c0;
    ctype const c7 = B.x * M.c1 + B.y * M.c2 + B.z * M.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv ^ vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& M, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = M.c2 * v.z - M.c3 * v.y;
    ctype const c5 = -M.c1 * v.z + M.c3 * v.x;
    ctype const c6 = M.c1 * v.y - M.c2 * v.x;
    ctype const c7 = M.c4 * v.x + M.c5 * v.y + M.c6 * v.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: vec ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = v.y * M.c3 - v.z * M.c2;
    ctype const c5 = -v.x * M.c3 + v.z * M.c1;
    ctype const c6 = v.x * M.c2 - v.y * M.c1;
    ctype const c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv ^ s -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& M, Scalar3d<U> s)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: s ^ mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, MVec3d<U> const& M)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d wdg :: mv_e ^ mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> wdg(MVec3d_E<T> const& A,
                                                 MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c0;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d wdg :: mv_e ^ mv_u -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> wdg(MVec3d_E<T> const& A,
                                                 MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d wdg :: mv_u ^ mv_e -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> wdg(MVec3d_U<T> const& A,
                                                 MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d wdg :: mv_e ^ ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(MVec3d_E<T> const& M, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(M.c0 * ctype(ps));
}

// ega3d wdg :: ps ^ mv_e -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(PScalar3d<T> ps, MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * M.c0);
}

// ega3d wdg :: mv_e ^ bivec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(MVec3d_E<T> const& M, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * B.x;
    ctype const c1 = M.c0 * B.y;
    ctype const c2 = M.c0 * B.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: bivec ^ mv_e -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B, MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c0;
    ctype const c1 = B.y * M.c0;
    ctype const c2 = B.z * M.c0;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: mv_e ^ vec -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> wdg(MVec3d_E<T> const& M, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    ctype const c3 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d wdg :: vec ^ mv_e -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> wdg(Vec3d<T> const& v, MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    ctype const c3 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d wdg :: mv_e ^ s -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> wdg(MVec3d_E<T> const& M, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d wdg :: s ^ mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> wdg(Scalar3d<T> s, MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d wdg :: mv_u ^ mv_u -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(MVec3d_U<T> const& A,
                                                MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c2 - A.c2 * B.c1;
    ctype const c1 = -A.c0 * B.c2 + A.c2 * B.c0;
    ctype const c2 = A.c0 * B.c1 - A.c1 * B.c0;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: mv_u ^ ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] MVec3d_U<T> const&,
                                                 [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: ps ^ mv_u -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3d<T>,
                                                 [[maybe_unused]] MVec3d_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: mv_u ^ bivec -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(MVec3d_U<T> const& M,
                                                  BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(M.c0 * B.x + M.c1 * B.y + M.c2 * B.z);
}

// ega3d wdg :: bivec ^ mv_u -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B,
                                                  MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(B.x * M.c0 + B.y * M.c1 + B.z * M.c2);
}

// ega3d wdg :: mv_u ^ vec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(MVec3d_U<T> const& M, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.z - M.c2 * v.y;
    ctype const c1 = -M.c0 * v.z + M.c2 * v.x;
    ctype const c2 = M.c0 * v.y - M.c1 * v.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: vec ^ mv_u -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v, MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * M.c2 - v.z * M.c1;
    ctype const c1 = -v.x * M.c2 + v.z * M.c0;
    ctype const c2 = v.x * M.c1 - v.y * M.c0;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: mv_u ^ s -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> wdg(MVec3d_U<T> const& M, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d wdg :: s ^ mv_u -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> wdg(Scalar3d<T> s, MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d wdg :: ps ^ ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3d<T>,
                                                 [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: ps ^ bivec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3d<T>,
                                                 [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: bivec ^ ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3d<T> const&,
                                                 [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: ps ^ vec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3d<T>,
                                                 [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: vec ^ ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec3d<T> const&,
                                                 [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: ps ^ s -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * ctype(s));
}

// ega3d wdg :: s ^ ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(s) * ctype(ps));
}

// ega3d wdg :: bivec ^ bivec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3d<T> const&,
                                                 [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d wdg :: bivec ^ vec -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(B.x * v.x + B.y * v.y + B.z * v.z);
}

// ega3d wdg :: vec ^ bivec -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(v.x * B.x + v.y * B.y + v.z * B.z);
}

// ega3d wdg :: bivec ^ s -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(s);
    ctype const c1 = B.y * ctype(s);
    ctype const c2 = B.z * ctype(s);
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: s ^ bivec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.x;
    ctype const c1 = ctype(s) * B.y;
    ctype const c2 = ctype(s) * B.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: vec ^ vec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.y * v2.z - v1.z * v2.y;
    ctype const c1 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c2 = v1.x * v2.y - v1.y * v2.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: vec ^ s -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: s ^ vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d wdg :: s ^ s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product for complements)
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product, just depending on the outer product (wdg)
// as well as the complement (cmpl).
//
// This definition is NOT connect directly to the geometric product, but to the outer
// product, exclusively.
//
//        rwdg(ul, ur) = cmpl(wdg(cmpl(ul),cmpl(ur))) = cmpl(cmpl(ul) ^ cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// ega3d rwdg :: rwdg(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
                     A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c1 * B.c7 + A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c1;
    ctype const c2 = A.c2 * B.c7 - A.c4 * B.c6 + A.c6 * B.c4 + A.c7 * B.c2;
    ctype const c3 = A.c3 * B.c7 + A.c4 * B.c5 - A.c5 * B.c4 + A.c7 * B.c3;
    ctype const c4 = A.c4 * B.c7 + A.c7 * B.c4;
    ctype const c5 = A.c5 * B.c7 + A.c7 * B.c5;
    ctype const c6 = A.c6 * B.c7 + A.c7 * B.c6;
    ctype const c7 = A.c7 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv, mv_e) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& A, MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c7 * B.c0;
    ctype const c1 = A.c5 * B.c3 - A.c6 * B.c2;
    ctype const c2 = -A.c4 * B.c3 + A.c6 * B.c1;
    ctype const c3 = A.c4 * B.c2 - A.c5 * B.c1;
    ctype const c4 = A.c7 * B.c1;
    ctype const c5 = A.c7 * B.c2;
    ctype const c6 = A.c7 * B.c3;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv_e,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d_E<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c7 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype const c1 = A.c2 * B.c6 - A.c3 * B.c5;
    ctype const c2 = -A.c1 * B.c6 + A.c3 * B.c4;
    ctype const c3 = A.c1 * B.c5 - A.c2 * B.c4;
    ctype const c4 = A.c1 * B.c7;
    ctype const c5 = A.c2 * B.c7;
    ctype const c6 = A.c3 * B.c7;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv,mv_u) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& A, MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2;
    ctype const c1 = A.c1 * B.c3 + A.c7 * B.c0;
    ctype const c2 = A.c2 * B.c3 + A.c7 * B.c1;
    ctype const c3 = A.c3 * B.c3 + A.c7 * B.c2;
    ctype const c4 = A.c4 * B.c3;
    ctype const c5 = A.c5 * B.c3;
    ctype const c6 = A.c6 * B.c3;
    ctype const c7 = A.c7 * B.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv_u,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d_U<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c0;
    ctype const c1 = A.c0 * B.c7 + A.c3 * B.c1;
    ctype const c2 = A.c1 * B.c7 + A.c3 * B.c2;
    ctype const c3 = A.c2 * B.c7 + A.c3 * B.c3;
    ctype const c4 = A.c3 * B.c4;
    ctype const c5 = A.c3 * B.c5;
    ctype const c6 = A.c3 * B.c6;
    ctype const c7 = A.c3 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv,ps) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& M, PScalar3d<U> ps)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(ps,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, MVec3d<U> const& M)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv,bivec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& M, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.x + M.c2 * B.y + M.c3 * B.z;
    ctype const c1 = M.c5 * B.z - M.c6 * B.y;
    ctype const c2 = -M.c4 * B.z + M.c6 * B.x;
    ctype const c3 = M.c4 * B.y - M.c5 * B.x;
    ctype const c4 = M.c7 * B.x;
    ctype const c5 = M.c7 * B.y;
    ctype const c6 = M.c7 * B.z;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(bivec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c1 + B.y * M.c2 + B.z * M.c3;
    ctype const c1 = B.y * M.c6 - B.z * M.c5;
    ctype const c2 = -B.x * M.c6 + B.z * M.c4;
    ctype const c3 = B.x * M.c5 - B.y * M.c4;
    ctype const c4 = B.x * M.c7;
    ctype const c5 = B.y * M.c7;
    ctype const c6 = B.z * M.c7;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv,vec) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& M, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c4 * v.x + M.c5 * v.y + M.c6 * v.z;
    ctype const c1 = M.c7 * v.x;
    ctype const c2 = M.c7 * v.y;
    ctype const c3 = M.c7 * v.z;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(vec,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rwdg(Vec3d<T> const& v, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6;
    ctype const c1 = v.x * M.c7;
    ctype const c2 = v.y * M.c7;
    ctype const c3 = v.z * M.c7;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d rwdg :: rwdg(mv,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& M, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c7 * ctype(s));
}

// ega3d rwdg :: rwdg(s,mv) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Scalar3d<T> s, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s) * M.c7);
}

// ega3d rwdg :: rwdg(mv_e,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(MVec3d_E<T> const& A, MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c2 * B.c3 - A.c3 * B.c2;
    ctype const c1 = -A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c2 = A.c1 * B.c2 - A.c2 * B.c1;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(mv_e,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> rwdg(MVec3d_E<T> const& A,
                                                  MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2;
    ctype const c1 = A.c1 * B.c3;
    ctype const c2 = A.c2 * B.c3;
    ctype const c3 = A.c3 * B.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d rwdg :: rwdg(mv_u,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> rwdg(MVec3d_U<T> const& A,
                                                  MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0;
    ctype const c1 = A.c3 * B.c1;
    ctype const c2 = A.c3 * B.c2;
    ctype const c3 = A.c3 * B.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d rwdg :: rwdg(mv_e,ps) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> rwdg(MVec3d_E<T> const& M, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d rwdg :: rwdg(ps,mv_e) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d rwdg :: rwdg(mv_e,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(MVec3d_E<T> const& M, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * B.z - M.c3 * B.y;
    ctype const c1 = -M.c1 * B.z + M.c3 * B.x;
    ctype const c2 = M.c1 * B.y - M.c2 * B.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(bivec,mv_e) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B, MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * M.c3 - B.z * M.c2;
    ctype const c1 = -B.x * M.c3 + B.z * M.c1;
    ctype const c2 = B.x * M.c2 - B.y * M.c1;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(mv_e,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(MVec3d_E<T> const& M, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c1 * v.x + M.c2 * v.y + M.c3 * v.z);
}

// ega3d rwdg :: rwdg(vec,mv_e) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Vec3d<T> const& v, MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v.x * M.c1 + v.y * M.c2 + v.z * M.c3);
}

// ega3d rwdg :: rwdg(mv_e,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] MVec3d_E<T> const&,
                                                  [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d rwdg :: rwdg(s,mv_e) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3d<T>,
                                                  [[maybe_unused]] MVec3d_E<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d rwdg :: rwdg(mv_u,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> rwdg(MVec3d_U<T> const& A,
                                                  MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype const c1 = A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c2 = A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c3 = A.c3 * B.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d rwdg :: rwdg(mv_u,ps) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> rwdg(MVec3d_U<T> const& M, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(ps);
    ctype const c1 = M.c1 * ctype(ps);
    ctype const c2 = M.c2 * ctype(ps);
    ctype const c3 = M.c3 * ctype(ps);
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d rwdg :: rwdg(ps,mv_u) -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c0;
    ctype const c1 = ctype(ps) * M.c1;
    ctype const c2 = ctype(ps) * M.c2;
    ctype const c3 = ctype(ps) * M.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d rwdg :: rwdg(mv_u,bivec) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> rwdg(MVec3d_U<T> const& M,
                                                  BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * B.x + M.c1 * B.y + M.c2 * B.z;
    ctype const c1 = M.c3 * B.x;
    ctype const c2 = M.c3 * B.y;
    ctype const c3 = M.c3 * B.z;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d rwdg :: rwdg(bivec,mv_u) -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B,
                                                  MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c0 + B.y * M.c1 + B.z * M.c2;
    ctype const c1 = B.x * M.c3;
    ctype const c2 = B.y * M.c3;
    ctype const c3 = B.z * M.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d rwdg :: rwdg(mv_u,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(MVec3d_U<T> const& M, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c3 * v.x;
    ctype const c1 = M.c3 * v.y;
    ctype const c2 = M.c3 * v.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(vec,mv_u) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(Vec3d<T> const& v, MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c3;
    ctype const c1 = v.y * M.c3;
    ctype const c2 = v.z * M.c3;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(mv_u,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(MVec3d_U<T> const& M, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c3 * ctype(s));
}

// ega3d rwdg :: rwdg(s,mv_u) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Scalar3d<T> s, MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s) * M.c3);
}

// ega3d rwdg :: rwdg(ps,ps) -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps1, PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

// ega3d rwdg :: rwdg(ps,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.x;
    ctype const c1 = ctype(ps) * B.y;
    ctype const c2 = ctype(ps) * B.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(bivec,ps) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(ps);
    ctype const c1 = B.y * ctype(ps);
    ctype const c2 = B.z * ctype(ps);
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(ps,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(vec,ps) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(Vec3d<T> const& v, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(ps,s) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps) * ctype(s));
}

// ega3d rwdg :: rwdg(s,ps) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s) * ctype(ps));
}

// ega3d rwdg :: rwdg(bivec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B1, BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B1.y * B2.z - B1.z * B2.y;
    ctype const c1 = -B1.x * B2.z + B1.z * B2.x;
    ctype const c2 = B1.x * B2.y - B1.y * B2.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d rwdg :: rwdg(bivec,vec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B.x * v.x + B.y * v.y + B.z * v.z);
}

// ega3d rwdg :: rwdg(vec,bivec) -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v.x * B.x + v.y * B.y + v.z * B.z);
}

// ega3d rwdg :: rwdg(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec3d<T> const&,
                                                  [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d rwdg :: rwdg(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3d<T>,
                                                  [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d rwdg :: rwdg(vec,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3d<T> const&,
                                                  [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d rwdg :: rwdg(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3d<T> const&,
                                                  [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d rwdg :: rwdg(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3d<T>,
                                                  [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d rwdg :: rwdg(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3d<T>,
                                                  [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// twdg1(ps,vec) = vec -> identical to geometric product gpr(ps,vec)
//                     -> and identical to the left contraction vec << ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> twdg1(PScalar3d<T> ps, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// twdg1(vec,ps) = vec -> identical to geometric product gpr(vec,ps)
//                     -> and identical to the right contraction ps >> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> twdg1(Vec3d<T> const& v, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    return BiVec3d<ctype>(c0, c1, c2);
}

// twdg1(bivec,bivec) = bivector -> identical to commutator product cmt(bivec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> twdg1(BiVec3d<T> const& B1,
                                                  BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.y * B2.z + B1.z * B2.y;
    ctype const c1 = B1.x * B2.z - B1.z * B2.x;
    ctype const c2 = -B1.x * B2.y + B1.y * B2.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

// twdg1(bivec,vec) = vector -> identical to left contraction operator<<(vec,bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> twdg1(BiVec3d<T> const& B, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.y * v.z + B.z * v.y;
    ctype const c1 = B.x * v.z - B.z * v.x;
    ctype const c2 = -B.x * v.y + B.y * v.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// twdg1(vec,bivec) = vector -> identical to right contraction operator>>(bivec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> twdg1(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.z + v.z * B.y;
    ctype const c1 = v.x * B.z - v.z * B.x;
    ctype const c2 = -v.x * B.y + v.y * B.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// twdg1(vec,vec) = scalar -> identical to dot product dot(vec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> twdg1(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}


////////////////////////////////////////////////////////////////////////////////
// regressive transwedge product (k=1)
////////////////////////////////////////////////////////////////////////////////

// rtwdg1(vec,vec) = pseudoscalar -> identical to regressive dot product rdot(vec,vec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> rtwdg1(BiVec3d<T> const& B1,
                                                     BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(B1.x * B2.x + B1.y * B2.y + B1.z * B2.z);
}

// rtwdg1(bivec,vec) = bivec -> identical to right expansion r_expand(vec,bivec)
//                           ->                             = vec ^ dual(bivec)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> rtwdg1(BiVec3d<T> const& B, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.y * v.z + B.z * v.y;
    ctype const c1 = B.x * v.z - B.z * v.x;
    ctype const c2 = -B.x * v.y + B.y * v.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

// rtwdg1(vec,bivec) = bivec -> identical to left expansion l_expand(bivec,vec)
//                           ->                             = dual(bivec) ^ vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> rtwdg1(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.z + v.z * B.y;
    ctype const c1 = v.x * B.z - v.z * B.x;
    ctype const c2 = -v.x * B.y + v.y * B.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

// rtwdg1(bivec,s) = vec -> returns a vector perpendicular to bivec scaled by s
// identical to right expansion: -l_expand(bivec,s) = r_expand(s,bivec)
// => sign change makes sense since expansions are derived from the wedge product
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rtwdg1(BiVec3d<T> const& B, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(s);
    ctype const c1 = B.y * ctype(s);
    ctype const c2 = B.z * ctype(s);
    return Vec3d<ctype>(c0, c1, c2);
}

// rtwdg1(s,bivec) -> returns a vector perpendicular to bivec scaled by s
// identical to expansions: -r_expand(s,bivec) = l_expand(bivec,s)
// => sign change makes sense since expansions are derived from the wedge product
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rtwdg1(Scalar3d<T> s, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.x;
    ctype const c1 = ctype(s) * B.y;
    ctype const c2 = ctype(s) * B.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// rtwdg1(vec1,vec2) = vector -> identical to cross product cross(vec2,vec1)
//                                                       = -cross(vec1,vec2)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rtwdg1(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v1.y * v2.z + v1.z * v2.y;
    ctype const c1 = v1.x * v2.z - v1.z * v2.x;
    ctype const c2 = -v1.x * v2.y + v1.y * v2.x;
    return Vec3d<ctype>(c0, c1, c2);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// The contraction A << B of an a-blade A onto a b-blade B is a sub-blade of B
// of grade gr(b)-gr(a) which is perpendicular to A, and linear in both arguments
//
// i.e. the resulting object lies in B and is perpendicular to A
////////////////////////////////////////////////////////////////////////////////

// ega3d l_contract :: mv << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c4 * B.c4 +
                     A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c6 - A.c3 * B.c5 + A.c4 * B.c7;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c6 + A.c3 * B.c4 + A.c5 * B.c7;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 + A.c6 * B.c7;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c7;
    ctype const c5 = A.c0 * B.c5 + A.c2 * B.c7;
    ctype const c6 = A.c0 * B.c6 + A.c3 * B.c7;
    ctype const c7 = A.c0 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv << mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& A,
                                                      MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c4 * B.c1 + A.c5 * B.c2 + A.c6 * B.c3;
    ctype const c1 = A.c2 * B.c3 - A.c3 * B.c2;
    ctype const c2 = -A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c3 = A.c1 * B.c2 - A.c2 * B.c1;
    ctype const c4 = A.c0 * B.c1;
    ctype const c5 = A.c0 * B.c2;
    ctype const c6 = A.c0 * B.c3;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv_e << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d_E<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c7;
    ctype const c2 = A.c0 * B.c2 + A.c2 * B.c7;
    ctype const c3 = A.c0 * B.c3 + A.c3 * B.c7;
    ctype const c4 = A.c0 * B.c4;
    ctype const c5 = A.c0 * B.c5;
    ctype const c6 = A.c0 * B.c6;
    ctype const c7 = A.c0 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv << mv_u -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& A,
                                                      MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 + A.c7 * B.c3;
    ctype const c1 = A.c0 * B.c0 + A.c4 * B.c3;
    ctype const c2 = A.c0 * B.c1 + A.c5 * B.c3;
    ctype const c3 = A.c0 * B.c2 + A.c6 * B.c3;
    ctype const c4 = A.c1 * B.c3;
    ctype const c5 = A.c2 * B.c3;
    ctype const c6 = A.c3 * B.c3;
    ctype const c7 = A.c0 * B.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv_u << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d_U<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c7;
    ctype const c1 = A.c1 * B.c6 - A.c2 * B.c5;
    ctype const c2 = -A.c0 * B.c6 + A.c2 * B.c4;
    ctype const c3 = A.c0 * B.c5 - A.c1 * B.c4;
    ctype const c4 = A.c0 * B.c7;
    ctype const c5 = A.c1 * B.c7;
    ctype const c6 = A.c2 * B.c7;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv << ps -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& M, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c7 * ctype(ps);
    ctype const c1 = M.c4 * ctype(ps);
    ctype const c2 = M.c5 * ctype(ps);
    ctype const c3 = M.c6 * ctype(ps);
    ctype const c4 = M.c1 * ctype(ps);
    ctype const c5 = M.c2 * ctype(ps);
    ctype const c6 = M.c3 * ctype(ps);
    ctype const c7 = M.c0 * ctype(ps);
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: ps << mv -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(PScalar3d<T> ps,
                                                        MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps) * M.c7);
}

// ega3d l_contract :: mv << bivec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& M,
                                                      BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c4 * B.x + M.c5 * B.y + M.c6 * B.z;
    ctype const c1 = M.c2 * B.z - M.c3 * B.y;
    ctype const c2 = -M.c1 * B.z + M.c3 * B.x;
    ctype const c3 = M.c1 * B.y - M.c2 * B.x;
    ctype const c4 = M.c0 * B.x;
    ctype const c5 = M.c0 * B.y;
    ctype const c6 = M.c0 * B.z;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: bivec << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B,
                                                      MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c4 + B.y * M.c5 + B.z * M.c6;
    ctype const c1 = B.x * M.c7;
    ctype const c2 = B.y * M.c7;
    ctype const c3 = B.z * M.c7;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv << vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& M,
                                                      Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    ctype const c1 = M.c0 * v.x;
    ctype const c2 = M.c0 * v.y;
    ctype const c3 = M.c0 * v.z;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: vec << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v,
                                                      MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    ctype const c1 = v.y * M.c6 - v.z * M.c5;
    ctype const c2 = -v.x * M.c6 + v.z * M.c4;
    ctype const c3 = v.x * M.c5 - v.y * M.c4;
    ctype const c4 = v.x * M.c7;
    ctype const c5 = v.y * M.c7;
    ctype const c6 = v.z * M.c7;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv << s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& M, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c0 * ctype(s));
}

// ega3d l_contract :: s << mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s, MVec3d<U> const& M)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d l_contract :: mv_e << mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator<<(MVec3d_E<T> const& A,
                                                        MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c1;
    ctype const c2 = A.c0 * B.c2;
    ctype const c3 = A.c0 * B.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d l_contract :: mv_e << mv_u -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator<<(MVec3d_E<T> const& A,
                                                        MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c3;
    ctype const c1 = A.c0 * B.c1 + A.c2 * B.c3;
    ctype const c2 = A.c0 * B.c2 + A.c3 * B.c3;
    ctype const c3 = A.c0 * B.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d l_contract :: mv_u << mv_e -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(MVec3d_U<T> const& A,
                                                     MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c3 - A.c2 * B.c2;
    ctype const c1 = -A.c0 * B.c3 + A.c2 * B.c1;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c1;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: mv_e << ps -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator<<(MVec3d_E<T> const& M,
                                                        PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * ctype(ps);
    ctype const c1 = M.c2 * ctype(ps);
    ctype const c2 = M.c3 * ctype(ps);
    ctype const c3 = M.c0 * ctype(ps);
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d l_contract :: ps << mv_e -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3d<T>, [[maybe_unused]] MVec3d_E<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: mv_e << bivec -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator<<(MVec3d_E<T> const& M,
                                                        BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.x + M.c2 * B.y + M.c3 * B.z;
    ctype const c1 = M.c0 * B.x;
    ctype const c2 = M.c0 * B.y;
    ctype const c3 = M.c0 * B.z;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d l_contract :: bivec << mv_e -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B,
                                                        MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B.x * M.c1 + B.y * M.c2 + B.z * M.c3);
}

// ega3d l_contract :: mv_e << vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(MVec3d_E<T> const& M,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x;
    ctype const c1 = M.c0 * v.y;
    ctype const c2 = M.c0 * v.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: vec << mv_e -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v,
                                                     MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * M.c3 - v.z * M.c2;
    ctype const c1 = -v.x * M.c3 + v.z * M.c1;
    ctype const c2 = v.x * M.c2 - v.y * M.c1;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: mv_e << s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(MVec3d_E<T> const& M,
                                                        Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c0 * ctype(s));
}

// ega3d l_contract :: s << mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                        MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d l_contract :: mv_u << mv_u -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator<<(MVec3d_U<T> const& A,
                                                        MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c3;
    ctype const c2 = A.c1 * B.c3;
    ctype const c3 = A.c2 * B.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d l_contract :: mv_u << ps -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator<<(MVec3d_U<T> const& M,
                                                        PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c3 * ctype(ps);
    ctype const c1 = M.c0 * ctype(ps);
    ctype const c2 = M.c1 * ctype(ps);
    ctype const c3 = M.c2 * ctype(ps);
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d l_contract :: ps << mv_u -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(PScalar3d<T> ps,
                                                        MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps) * M.c3);
}

// ega3d l_contract :: mv_u << bivec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(MVec3d_U<T> const& M,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * B.z - M.c2 * B.y;
    ctype const c1 = -M.c0 * B.z + M.c2 * B.x;
    ctype const c2 = M.c0 * B.y - M.c1 * B.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: bivec << mv_u -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B,
                                                     MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c3;
    ctype const c1 = B.y * M.c3;
    ctype const c2 = B.z * M.c3;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: mv_u << vec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(MVec3d_U<T> const& M,
                                                        Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c0 * v.x + M.c1 * v.y + M.c2 * v.z);
}

// ega3d l_contract :: vec << mv_u -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v,
                                                        MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 + v.z * M.c2;
    ctype const c1 = v.x * M.c3;
    ctype const c2 = v.y * M.c3;
    ctype const c3 = v.z * M.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d l_contract :: mv_u << s -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] MVec3d_U<T> const&, [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: s << mv_u -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                        MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d l_contract :: ps << ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(PScalar3d<T> ps1,
                                                        PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

// ega3d l_contract :: ps << bivec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3d<T>, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: bivec << ps -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(ps);
    ctype const c1 = B.y * ctype(ps);
    ctype const c2 = B.z * ctype(ps);
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: ps << vec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3d<T>,
                                                        [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: vec << ps -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: ps << s -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3d<T>,
                                                        [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: s << ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(s) * ctype(ps));
}

// ega3d l_contract :: bivec << bivec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B1,
                                                        BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B1.x * B2.x + B1.y * B2.y + B1.z * B2.z);
}

// ega3d l_contract :: bivec << vec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: vec << bivec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.z - v.z * B.y;
    ctype const c1 = -v.x * B.z + v.z * B.x;
    ctype const c2 = v.x * B.y - v.y * B.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: bivec << s -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: s << bivec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.x;
    ctype const c1 = ctype(s) * B.y;
    ctype const c2 = ctype(s) * B.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: vec << vec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v1,
                                                        Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// ega3d l_contract :: vec << s -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<([[maybe_unused]] Vec3d<T> const&,
                                                        [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d l_contract :: s << vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d l_contract :: s << s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The resulting object lies in A and is perpendicular to B
////////////////////////////////////////////////////////////////////////////////

// ega3d r_contract :: mv >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c4 * B.c4 +
                     A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7;
    ctype const c1 = A.c1 * B.c0 + A.c5 * B.c3 - A.c6 * B.c2 + A.c7 * B.c4;
    ctype const c2 = A.c2 * B.c0 - A.c4 * B.c3 + A.c6 * B.c1 + A.c7 * B.c5;
    ctype const c3 = A.c3 * B.c0 + A.c4 * B.c2 - A.c5 * B.c1 + A.c7 * B.c6;
    ctype const c4 = A.c4 * B.c0 + A.c7 * B.c1;
    ctype const c5 = A.c5 * B.c0 + A.c7 * B.c2;
    ctype const c6 = A.c6 * B.c0 + A.c7 * B.c3;
    ctype const c7 = A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv >> mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& A,
                                                      MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c4 * B.c1 + A.c5 * B.c2 + A.c6 * B.c3;
    ctype const c1 = A.c1 * B.c0 + A.c7 * B.c1;
    ctype const c2 = A.c2 * B.c0 + A.c7 * B.c2;
    ctype const c3 = A.c3 * B.c0 + A.c7 * B.c3;
    ctype const c4 = A.c4 * B.c0;
    ctype const c5 = A.c5 * B.c0;
    ctype const c6 = A.c6 * B.c0;
    ctype const c7 = A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv_e >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d_E<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6;
    ctype const c1 = A.c2 * B.c3 - A.c3 * B.c2;
    ctype const c2 = -A.c1 * B.c3 + A.c3 * B.c1;
    ctype const c3 = A.c1 * B.c2 - A.c2 * B.c1;
    ctype const c4 = A.c1 * B.c0;
    ctype const c5 = A.c2 * B.c0;
    ctype const c6 = A.c3 * B.c0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv >> mv_u -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& A,
                                                      MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 + A.c7 * B.c3;
    ctype const c1 = A.c5 * B.c2 - A.c6 * B.c1;
    ctype const c2 = -A.c4 * B.c2 + A.c6 * B.c0;
    ctype const c3 = A.c4 * B.c1 - A.c5 * B.c0;
    ctype const c4 = A.c7 * B.c0;
    ctype const c5 = A.c7 * B.c1;
    ctype const c6 = A.c7 * B.c2;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv_u >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c7;
    ctype const c1 = A.c0 * B.c0 + A.c3 * B.c4;
    ctype const c2 = A.c1 * B.c0 + A.c3 * B.c5;
    ctype const c3 = A.c2 * B.c0 + A.c3 * B.c6;
    ctype const c4 = A.c3 * B.c1;
    ctype const c5 = A.c3 * B.c2;
    ctype const c6 = A.c3 * B.c3;
    ctype const c7 = A.c3 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv >> ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& M,
                                                        PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c7 * ctype(ps));
}

// ega3d r_contract :: ps >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c7;
    ctype const c1 = ctype(ps) * M.c4;
    ctype const c2 = ctype(ps) * M.c5;
    ctype const c3 = ctype(ps) * M.c6;
    ctype const c4 = ctype(ps) * M.c1;
    ctype const c5 = ctype(ps) * M.c2;
    ctype const c6 = ctype(ps) * M.c3;
    ctype const c7 = ctype(ps) * M.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv >> bivec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& M,
                                                      BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c4 * B.x + M.c5 * B.y + M.c6 * B.z;
    ctype const c1 = M.c7 * B.x;
    ctype const c2 = M.c7 * B.y;
    ctype const c3 = M.c7 * B.z;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: bivec >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B,
                                                      MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c4 + B.y * M.c5 + B.z * M.c6;
    ctype const c1 = B.y * M.c3 - B.z * M.c2;
    ctype const c2 = -B.x * M.c3 + B.z * M.c1;
    ctype const c3 = B.x * M.c2 - B.y * M.c1;
    ctype const c4 = B.x * M.c0;
    ctype const c5 = B.y * M.c0;
    ctype const c6 = B.z * M.c0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv >> vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& M,
                                                      Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    ctype const c1 = M.c5 * v.z - M.c6 * v.y;
    ctype const c2 = -M.c4 * v.z + M.c6 * v.x;
    ctype const c3 = M.c4 * v.y - M.c5 * v.x;
    ctype const c4 = M.c7 * v.x;
    ctype const c5 = M.c7 * v.y;
    ctype const c6 = M.c7 * v.z;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: vec >> mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(Vec3d<T> const& v,
                                                      MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    ctype const c1 = v.x * M.c0;
    ctype const c2 = v.y * M.c0;
    ctype const c3 = v.z * M.c0;
    ctype const c4 = 0.0;
    ctype const c5 = 0.0;
    ctype const c6 = 0.0;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: mv >> s -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& M, Scalar3d<U> s)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d r_contract :: s >> mv -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(Scalar3d<T> s, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s) * M.c0);
}

// ega3d r_contract :: mv_e >> mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator>>(MVec3d_E<T> const& A,
                                                        MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype const c1 = A.c1 * B.c0;
    ctype const c2 = A.c2 * B.c0;
    ctype const c3 = A.c3 * B.c0;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d r_contract :: mv_e >> mv_u -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(MVec3d_E<T> const& A,
                                                     MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c2 * B.c2 - A.c3 * B.c1;
    ctype const c1 = -A.c1 * B.c2 + A.c3 * B.c0;
    ctype const c2 = A.c1 * B.c1 - A.c2 * B.c0;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: mv_u >> mv_e -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& A,
                                                        MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype const c1 = A.c1 * B.c0 + A.c3 * B.c2;
    ctype const c2 = A.c2 * B.c0 + A.c3 * B.c3;
    ctype const c3 = A.c3 * B.c0;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d r_contract :: mv_e >> ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] MVec3d_E<T> const&, [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: ps >> mv_e -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps,
                                                        MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c1;
    ctype const c1 = ctype(ps) * M.c2;
    ctype const c2 = ctype(ps) * M.c3;
    ctype const c3 = ctype(ps) * M.c0;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d r_contract :: mv_e >> bivec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(MVec3d_E<T> const& M,
                                                        BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c1 * B.x + M.c2 * B.y + M.c3 * B.z);
}

// ega3d r_contract :: bivec >> mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B,
                                                        MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * M.c1 + B.y * M.c2 + B.z * M.c3;
    ctype const c1 = B.x * M.c0;
    ctype const c2 = B.y * M.c0;
    ctype const c3 = B.z * M.c0;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d r_contract :: mv_e >> vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(MVec3d_E<T> const& M,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c2 * v.z - M.c3 * v.y;
    ctype const c1 = -M.c1 * v.z + M.c3 * v.x;
    ctype const c2 = M.c1 * v.y - M.c2 * v.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: vec >> mv_e -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(Vec3d<T> const& v,
                                                     MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0;
    ctype const c1 = v.y * M.c0;
    ctype const c2 = v.z * M.c0;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: mv_e >> s -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator>>(MVec3d_E<T> const& M,
                                                        Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d r_contract :: s >> mv_e -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(Scalar3d<T> s,
                                                        MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s) * M.c0);
}

// ega3d r_contract :: mv_u >> mv_u -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& A,
                                                        MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype const c1 = A.c3 * B.c0;
    ctype const c2 = A.c3 * B.c1;
    ctype const c3 = A.c3 * B.c2;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d r_contract :: mv_u >> ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& M,
                                                        PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(M.c3 * ctype(ps));
}

// ega3d r_contract :: ps >> mv_u -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps,
                                                        MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * M.c3;
    ctype const c1 = ctype(ps) * M.c0;
    ctype const c2 = ctype(ps) * M.c1;
    ctype const c3 = ctype(ps) * M.c2;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d r_contract :: mv_u >> bivec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& M,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c3 * B.x;
    ctype const c1 = M.c3 * B.y;
    ctype const c2 = M.c3 * B.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: bivec >> mv_u -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B,
                                                     MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * M.c2 - B.z * M.c1;
    ctype const c1 = -B.x * M.c2 + B.z * M.c0;
    ctype const c2 = B.x * M.c1 - B.y * M.c0;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: mv_u >> vec -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& M,
                                                        Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y + M.c2 * v.z;
    ctype const c1 = M.c3 * v.x;
    ctype const c2 = M.c3 * v.y;
    ctype const c3 = M.c3 * v.z;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d r_contract :: vec >> mv_u -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(Vec3d<T> const& v,
                                                        MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v.x * M.c0 + v.y * M.c1 + v.z * M.c2);
}

// ega3d r_contract :: mv_u >> s -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& M,
                                                        Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d r_contract :: s >> mv_u -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3d<T>, [[maybe_unused]] MVec3d_U<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: ps >> ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps1,
                                                        PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

// ega3d r_contract :: ps >> bivec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * B.x;
    ctype const c1 = ctype(ps) * B.y;
    ctype const c2 = ctype(ps) * B.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: bivec >> ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: ps >> vec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: vec >> ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Vec3d<T> const&,
                                                        [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: ps >> s -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * ctype(s));
}

// ega3d r_contract :: s >> ps -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar3d<T>,
                                                        [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: bivec >> bivec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B1,
                                                        BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B1.x * B2.x + B1.y * B2.y + B1.z * B2.z);
}

// ega3d r_contract :: bivec >> vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * v.z - B.z * v.y;
    ctype const c1 = -B.x * v.z + B.z * v.x;
    ctype const c2 = B.x * v.y - B.y * v.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: vec >> bivec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3d<T> const&, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: bivec >> s -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(s);
    ctype const c1 = B.y * ctype(s);
    ctype const c2 = B.z * ctype(s);
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: s >> bivec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3d<T>, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: vec >> vec -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(Vec3d<T> const& v1,
                                                        Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// ega3d r_contract :: vec >> s -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator>>(Vec3d<T> const& v, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d r_contract :: s >> vec -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar3d<T>,
                                                        [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d r_contract :: s >> s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator>>(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// classical cross product (only defined in 3d)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> cross(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

// ega3d cmt :: cmt(mv,mv) -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> cmt(MVec3d<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = 0.0;
    ctype const c1 = -A.c2 * B.c6 + A.c3 * B.c5 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype const c2 = A.c1 * B.c6 - A.c3 * B.c4 + A.c4 * B.c3 - A.c6 * B.c1;
    ctype const c3 = -A.c1 * B.c5 + A.c2 * B.c4 - A.c4 * B.c2 + A.c5 * B.c1;
    ctype const c4 = A.c2 * B.c3 - A.c3 * B.c2 - A.c5 * B.c6 + A.c6 * B.c5;
    ctype const c5 = -A.c1 * B.c3 + A.c3 * B.c1 + A.c4 * B.c6 - A.c6 * B.c4;
    ctype const c6 = A.c1 * B.c2 - A.c2 * B.c1 - A.c4 * B.c5 + A.c5 * B.c4;
    ctype const c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d cmt :: cmt(ps,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3d<T>,
                                                 [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(ps,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3d<T>,
                                                 [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(bivec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec3d<T> const&,
                                                 [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(ps,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3d<T>,
                                                 [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(vec,ps) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec3d<T> const&,
                                                 [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(ps,s) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] PScalar3d<T>,
                                                  [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(s,ps) -> 0 ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3d<T>,
                                                  [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(bivec,bivec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> cmt(BiVec3d<T> const& B1,
                                                BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B1.y * B2.z + B1.z * B2.y;
    ctype const c1 = B1.x * B2.z - B1.z * B2.x;
    ctype const c2 = -B1.x * B2.y + B1.y * B2.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

// cmt(B,v) == -cmt(v,B) == (v << B)
// ega3d cmt :: cmt(bivec,vec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> cmt(BiVec3d<T> const& B, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.y * v.z + B.z * v.y;
    ctype const c1 = B.x * v.z - B.z * v.x;
    ctype const c2 = -B.x * v.y + B.y * v.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// cmt(v,B) == -cmt(B,v) == (B >> v)
// ega3d cmt :: cmt(vec,bivec) -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> cmt(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -v.y * B.z + v.z * B.y;
    ctype const c1 = v.x * B.z - v.z * B.x;
    ctype const c2 = -v.x * B.y + v.y * B.x;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d cmt :: cmt(bivec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] BiVec3d<T> const&,
                                                 [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(s,bivec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3d<T>,
                                                 [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(vec,vec) -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> cmt(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v1.y * v2.z - v1.z * v2.y;
    ctype const c1 = -v1.x * v2.z + v1.z * v2.x;
    ctype const c2 = v1.x * v2.y - v1.y * v2.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d cmt :: cmt(vec,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] Vec3d<T> const&,
                                                 [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(s,vec) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3d<T>,
                                                 [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// ega3d cmt :: cmt(s,s) -> 0
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> cmt([[maybe_unused]] Scalar3d<T>,
                                                 [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// ega3d gpr :: mv * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                     MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
                     A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 - A.c4 * B.c7 -
                     A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c4;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c4 * B.c3 -
                     A.c5 * B.c7 - A.c6 * B.c1 - A.c7 * B.c5;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
                     A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0 -
                     A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype const c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 +
                     A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype const c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 - A.c4 * B.c5 +
                     A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
                     A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv * mv_e -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                     MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype const c1 = A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2 - A.c7 * B.c1;
    ctype const c2 = A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1 - A.c7 * B.c2;
    ctype const c3 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0 - A.c7 * B.c3;
    ctype const c4 = A.c0 * B.c1 + A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype const c5 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1;
    ctype const c6 = A.c0 * B.c3 - A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c0;
    ctype const c7 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv_e * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                     MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6;
    ctype const c1 = A.c0 * B.c1 - A.c1 * B.c7 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c3 - A.c2 * B.c7 - A.c3 * B.c1;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c7;
    ctype const c4 = A.c0 * B.c4 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype const c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4;
    ctype const c6 = A.c0 * B.c6 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0;
    ctype const c7 = A.c0 * B.c7 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv * mv_u -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                     MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c7 * B.c3;
    ctype const c1 = A.c0 * B.c0 - A.c4 * B.c3 - A.c5 * B.c2 + A.c6 * B.c1;
    ctype const c2 = A.c0 * B.c1 + A.c4 * B.c2 - A.c5 * B.c3 - A.c6 * B.c0;
    ctype const c3 = A.c0 * B.c2 - A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c3;
    ctype const c4 = A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1 + A.c7 * B.c0;
    ctype const c5 = -A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0 + A.c7 * B.c1;
    ctype const c6 = A.c1 * B.c1 - A.c2 * B.c0 + A.c3 * B.c3 + A.c7 * B.c2;
    ctype const c7 = A.c0 * B.c3 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv_u * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                     MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c3 * B.c7;
    ctype const c1 = A.c0 * B.c0 - A.c1 * B.c6 + A.c2 * B.c5 - A.c3 * B.c4;
    ctype const c2 = A.c0 * B.c6 + A.c1 * B.c0 - A.c2 * B.c4 - A.c3 * B.c5;
    ctype const c3 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c0 - A.c3 * B.c6;
    ctype const c4 = A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1;
    ctype const c5 = -A.c0 * B.c3 + A.c1 * B.c7 + A.c2 * B.c1 + A.c3 * B.c2;
    ctype const c6 = A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c7 + A.c3 * B.c3;
    ctype const c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv * ps -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A, PScalar3d<U> ps)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: ps * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps, MVec3d<U> const& B)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv * bivec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& M,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c4 * B.x - M.c5 * B.y - M.c6 * B.z;
    ctype const c1 = -M.c2 * B.z + M.c3 * B.y - M.c7 * B.x;
    ctype const c2 = M.c1 * B.z - M.c3 * B.x - M.c7 * B.y;
    ctype const c3 = -M.c1 * B.y + M.c2 * B.x - M.c7 * B.z;
    ctype const c4 = M.c0 * B.x - M.c5 * B.z + M.c6 * B.y;
    ctype const c5 = M.c0 * B.y + M.c4 * B.z - M.c6 * B.x;
    ctype const c6 = M.c0 * B.z - M.c4 * B.y + M.c5 * B.x;
    ctype const c7 = M.c1 * B.x + M.c2 * B.y + M.c3 * B.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: bivec * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                     MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c4 - B.y * M.c5 - B.z * M.c6;
    ctype const c1 = -B.x * M.c7 - B.y * M.c3 + B.z * M.c2;
    ctype const c2 = B.x * M.c3 - B.y * M.c7 - B.z * M.c1;
    ctype const c3 = -B.x * M.c2 + B.y * M.c1 - B.z * M.c7;
    ctype const c4 = B.x * M.c0 - B.y * M.c6 + B.z * M.c5;
    ctype const c5 = B.x * M.c6 + B.y * M.c0 - B.z * M.c4;
    ctype const c6 = -B.x * M.c5 + B.y * M.c4 + B.z * M.c0;
    ctype const c7 = B.x * M.c1 + B.y * M.c2 + B.z * M.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv * vec -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& M,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    ctype const c1 = M.c0 * v.x - M.c5 * v.z + M.c6 * v.y;
    ctype const c2 = M.c0 * v.y + M.c4 * v.z - M.c6 * v.x;
    ctype const c3 = M.c0 * v.z - M.c4 * v.y + M.c5 * v.x;
    ctype const c4 = M.c2 * v.z - M.c3 * v.y + M.c7 * v.x;
    ctype const c5 = -M.c1 * v.z + M.c3 * v.x + M.c7 * v.y;
    ctype const c6 = M.c1 * v.y - M.c2 * v.x + M.c7 * v.z;
    ctype const c7 = M.c4 * v.x + M.c5 * v.y + M.c6 * v.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: vec * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                     MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    ctype const c1 = v.x * M.c0 - v.y * M.c6 + v.z * M.c5;
    ctype const c2 = v.x * M.c6 + v.y * M.c0 - v.z * M.c4;
    ctype const c3 = -v.x * M.c5 + v.y * M.c4 + v.z * M.c0;
    ctype const c4 = v.x * M.c7 + v.y * M.c3 - v.z * M.c2;
    ctype const c5 = -v.x * M.c3 + v.y * M.c7 + v.z * M.c1;
    ctype const c6 = v.x * M.c2 - v.y * M.c1 + v.z * M.c7;
    ctype const c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv * s -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& M, Scalar3d<U> s)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: s * mv -> mv
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s, MVec3d<U> const& M)
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
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// ega3d gpr :: mv_e * mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                       MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype const c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype const c3 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: mv_e * mv_u -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                       MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1;
    ctype const c1 = A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c3 - A.c3 * B.c0;
    ctype const c2 = A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3;
    ctype const c3 = A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: mv_u * mv_e -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                       MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 - A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1;
    ctype const c1 = A.c0 * B.c3 + A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2;
    ctype const c2 = -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3;
    ctype const c3 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: mv_e * ps -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& M,
                                                       PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * ctype(ps);
    ctype const c1 = -M.c2 * ctype(ps);
    ctype const c2 = -M.c3 * ctype(ps);
    ctype const c3 = M.c0 * ctype(ps);
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: ps * mv_e -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(PScalar3d<T> ps,
                                                       MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c1;
    ctype const c1 = -ctype(ps) * M.c2;
    ctype const c2 = -ctype(ps) * M.c3;
    ctype const c3 = ctype(ps) * M.c0;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: mv_e * bivec -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& M,
                                                       BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.x - M.c2 * B.y - M.c3 * B.z;
    ctype const c1 = M.c0 * B.x - M.c2 * B.z + M.c3 * B.y;
    ctype const c2 = M.c0 * B.y + M.c1 * B.z - M.c3 * B.x;
    ctype const c3 = M.c0 * B.z - M.c1 * B.y + M.c2 * B.x;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: bivec * mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                       MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c1 - B.y * M.c2 - B.z * M.c3;
    ctype const c1 = B.x * M.c0 - B.y * M.c3 + B.z * M.c2;
    ctype const c2 = B.x * M.c3 + B.y * M.c0 - B.z * M.c1;
    ctype const c3 = -B.x * M.c2 + B.y * M.c1 + B.z * M.c0;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: mv_e * vec -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& M,
                                                       Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x - M.c2 * v.z + M.c3 * v.y;
    ctype const c1 = M.c0 * v.y + M.c1 * v.z - M.c3 * v.x;
    ctype const c2 = M.c0 * v.z - M.c1 * v.y + M.c2 * v.x;
    ctype const c3 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: vec * mv_e -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                       MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 - v.y * M.c3 + v.z * M.c2;
    ctype const c1 = v.x * M.c3 + v.y * M.c0 - v.z * M.c1;
    ctype const c2 = -v.x * M.c2 + v.y * M.c1 + v.z * M.c0;
    ctype const c3 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: mv_e * s -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& M,
                                                       Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: s * mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                       MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: mv_u * mv_u -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                       MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3;
    ctype const c1 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    ctype const c2 = -A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1;
    ctype const c3 = A.c0 * B.c1 - A.c1 * B.c0 + A.c2 * B.c3 + A.c3 * B.c2;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: mv_u * ps -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& M,
                                                       PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c3 * ctype(ps);
    ctype const c1 = M.c0 * ctype(ps);
    ctype const c2 = M.c1 * ctype(ps);
    ctype const c3 = M.c2 * ctype(ps);
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: ps * mv_u -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(PScalar3d<T> ps,
                                                       MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * M.c3;
    ctype const c1 = ctype(ps) * M.c0;
    ctype const c2 = ctype(ps) * M.c1;
    ctype const c3 = ctype(ps) * M.c2;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: mv_u * bivec -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& M,
                                                       BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -M.c1 * B.z + M.c2 * B.y - M.c3 * B.x;
    ctype const c1 = M.c0 * B.z - M.c2 * B.x - M.c3 * B.y;
    ctype const c2 = -M.c0 * B.y + M.c1 * B.x - M.c3 * B.z;
    ctype const c3 = M.c0 * B.x + M.c1 * B.y + M.c2 * B.z;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: bivec * mv_u -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                       MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * M.c3 - B.y * M.c2 + B.z * M.c1;
    ctype const c1 = B.x * M.c2 - B.y * M.c3 - B.z * M.c0;
    ctype const c2 = -B.x * M.c1 + B.y * M.c0 - B.z * M.c3;
    ctype const c3 = B.x * M.c0 + B.y * M.c1 + B.z * M.c2;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: mv_u * vec -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& M,
                                                       Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * v.x + M.c1 * v.y + M.c2 * v.z;
    ctype const c1 = M.c1 * v.z - M.c2 * v.y + M.c3 * v.x;
    ctype const c2 = -M.c0 * v.z + M.c2 * v.x + M.c3 * v.y;
    ctype const c3 = M.c0 * v.y - M.c1 * v.x + M.c3 * v.z;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: vec * mv_u -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                       MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * M.c0 + v.y * M.c1 + v.z * M.c2;
    ctype const c1 = v.x * M.c3 + v.y * M.c2 - v.z * M.c1;
    ctype const c2 = -v.x * M.c2 + v.y * M.c3 + v.z * M.c0;
    ctype const c3 = v.x * M.c1 - v.y * M.c0 + v.z * M.c3;
    return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
}

// ega3d gpr :: mv_u * s -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& M,
                                                       Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = M.c0 * ctype(s);
    ctype const c1 = M.c1 * ctype(s);
    ctype const c2 = M.c2 * ctype(s);
    ctype const c3 = M.c3 * ctype(s);
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: s * mv_u -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                       MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * M.c0;
    ctype const c1 = ctype(s) * M.c1;
    ctype const c2 = ctype(s) * M.c2;
    ctype const c3 = ctype(s) * M.c3;
    return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
}

// ega3d gpr :: ps * ps -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps1, PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(-ctype(ps1) * ctype(ps2));
}

// ega3d gpr :: ps * bivec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -ctype(ps) * B.x;
    ctype const c1 = -ctype(ps) * B.y;
    ctype const c2 = -ctype(ps) * B.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d gpr :: bivec * ps -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = -B.x * ctype(ps);
    ctype const c1 = -B.y * ctype(ps);
    ctype const c2 = -B.z * ctype(ps);
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d gpr :: ps * vec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(ps) * v.x;
    ctype const c1 = ctype(ps) * v.y;
    ctype const c2 = ctype(ps) * v.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d gpr :: vec * ps -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(ps);
    ctype const c1 = v.y * ctype(ps);
    ctype const c2 = v.z * ctype(ps);
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d gpr :: ps * s -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * ctype(s));
}

// ega3d gpr :: s * ps -> ps
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr PScalar3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(s) * ctype(ps));
}

// geometric product A * B between two bivectors:
//
// the full geometric bivector product only exists in >= 4d spaces:
// A * B = gr0(A * B) + gr2(A * B) + gr4(A * B) = -dot(A,B) + cmt(A,B) + wdg(A,B)
// In 3D we don't have gr4(A * B) and thus only the terms up to grade 3 remain.
// The bivector product AxB = cmt(A,B) = 0.5*(A*B-B*A) is called the commutator
// product.
//
// A * B = -dot(A,B) + cmt(A,B)             (in 3d)
//
// ega3d gpr :: bivec * bivec -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B1,
                                                       BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -B1.x * B2.x - B1.y * B2.y - B1.z * B2.z;
    // ctype const c1 = -B1.y * B2.z + B1.z * B2.y;
    // ctype const c2 = B1.x * B2.z - B1.z * B2.x;
    // ctype const c3 = -B1.x * B2.y + B1.y * B2.x;
    // return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
    return MVec3d_E<ctype>(-dot(B1, B2), cmt(B1, B2));
}

// geometric product B * v for a bivector B and a vector v
// B * v = (v << B) + wdg(B,v) = gr1(B * v) + gr3(B * v)
//
// HINT: if a full 3d multivector is required as result it must be converted
//       explicitly, since C++ does not allow overloading on different return types
//
// ega3d gpr :: bivec * vec -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                       Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -B.y * v.z + B.z * v.y;
    // ctype const c1 = B.x * v.z - B.z * v.x;
    // ctype const c2 = -B.x * v.y + B.y * v.x;
    // ctype const c3 = B.x * v.x + B.y * v.y + B.z * v.z;
    // return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
    return MVec3d_U<ctype>((v << B), wdg(B, v));
}

// geometric product v * B for a vector v and a bivector B
// v * B = (B >> v) + wdg(v,B) = gr1(v * B) + gr3(v * B)
//
// HINT: if a full 3d multivector is required as result it must be converted
//       explicitly, since C++ does not allow overloading on different return types
//
// ega3d gpr :: vec * bivec -> mv_u
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                       BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = -v.y * B.z + v.z * B.y;
    // ctype const c1 = v.x * B.z - v.z * B.x;
    // ctype const c2 = -v.x * B.y + v.y * B.x;
    // ctype const c3 = v.x * B.x + v.y * B.y + v.z * B.z;
    // return MVec3d_U<ctype>(Vec3d<ctype>(c0, c1, c2), PScalar3d<ctype>(c3));
    return MVec3d_U<ctype>((B >> v), wdg(v, B));
}

// ega3d gpr :: bivec * s -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.x * ctype(s);
    ctype const c1 = B.y * ctype(s);
    ctype const c2 = B.z * ctype(s);
    return BiVec3d<ctype>(c0, c1, c2);
}

// ega3d gpr :: s * bivec -> bivec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * B.x;
    ctype const c1 = ctype(s) * B.y;
    ctype const c2 = ctype(s) * B.z;
    return BiVec3d<ctype>(c0, c1, c2);
}

// geometric product a * b of two vectors
// a * b = (b >> a) + wdg(a,b) = gr0(a * b) + gr2(a * b)
//
// HINT: if a full 3d multivector is required as result, it must be converted
//       explicitly, since C++ does not allow overloading on different return types
//
// ega3d gpr :: vec * vec -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Vec3d<T> const& v1,
                                                       Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    // ctype const c0 = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    // ctype const c1 = v1.y * v2.z - v1.z * v2.y;
    // ctype const c2 = -v1.x * v2.z + v1.z * v2.x;
    // ctype const c3 = v1.x * v2.y - v1.y * v2.x;
    // return MVec3d_E<ctype>(Scalar3d<ctype>(c0), BiVec3d<ctype>(c1, c2, c3));
    return MVec3d_E<ctype>((v2 >> v1), wdg(v1, v2));
}

// ega3d gpr :: vec * s -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.x * ctype(s);
    ctype const c1 = v.y * ctype(s);
    ctype const c2 = v.z * ctype(s);
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d gpr :: s * vec -> vec
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = ctype(s) * v.x;
    ctype const c1 = ctype(s) * v.y;
    ctype const c2 = ctype(s) * v.z;
    return Vec3d<ctype>(c0, c1, c2);
}

// ega3d gpr :: s * s -> s
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
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
inline Scalar3d<T> inv(Scalar3d<T> s)
{
    T sq_n = nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar3d<T>(rev(s) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec3d<T> inv(Vec3d<T> const& v)
{
    T sq_n = nrm_sq(v);
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n;
    return Vec3d<T>(rev(v) * inv);
}

// return the multiplicative inverse of the bivector
template <typename T>
    requires(numeric_type<T>)
inline BiVec3d<T> inv(BiVec3d<T> const& B)
{
    T sq_n = nrm_sq(B);
    hd::ga::detail::check_normalization<T>(sq_n, "bivector");
    T inv = T(1.0) / sq_n;
    return BiVec3d<T>(rev(B) * inv); // minus sign due to reversion
}

template <typename T>
    requires(numeric_type<T>)
inline PScalar3d<T> inv(PScalar3d<T> ps)
{
    T sq_n = nrm_sq(ps);
    hd::ga::detail::check_normalization<T>(sq_n, "pseudoscalar");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector

    return PScalar3d<T>(rev(ps) * inv); // minus sign due to reversion
}

// return the multiplicative inverse of the quaternion (inv(z) = 1/nrm_sq(z)*rev(z))
template <typename T>
    requires(numeric_type<T>)
inline MVec3d_E<T> inv(MVec3d_E<T> const& E)
{
    T sq_n = nrm_sq(E);
    hd::ga::detail::check_normalization<T>(sq_n, "even-grade-multivector");
    T inv = T(1.0) / sq_n;
    return MVec3d_E<T>(rev(E) * inv);
}

// return the multiplicative inverse of the odd-grade multivector
template <typename T>
    requires(numeric_type<T>)
inline MVec3d_U<T> inv(MVec3d_U<T> const& U)
{
    T sq_n = nrm_sq(U);
    hd::ga::detail::check_normalization<T>(sq_n, "odd-grade multivector");
    T inv = T(1.0) / sq_n;
    return MVec3d_U<T>(rev(U) * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016 left and a right inverse are the same (see paper
// of Hitzer, Sangwine)
template <typename T>
    requires(numeric_type<T>)
inline MVec3d<T> inv(MVec3d<T> const& M)
{
    T m_conjm = T(gr0(M * conj(M) * gr_inv(M) * rev(M)));
    hd::ga::detail::check_normalization<T>(std::abs(m_conjm), "multivector");
    T inv = T(1.0) / m_conjm;
    return MVec3d<T>(conj(M) * gr_inv(M) * rev(M) * inv);
}

////////////////////////////////////////////////////////////////////////////////
// left and right expansion
//
// given bivector B and vector v:
// return a bivector B' that represents a plane through the origin containing v
// and is perpendicular to B
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> l_expand(BiVec3d<T> const& B,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = B.y * v.z - B.z * v.y;
    ctype const c1 = -B.x * v.z + B.z * v.x;
    ctype const c2 = B.x * v.y - B.y * v.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> r_expand(Vec3d<T> const& v,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype const c0 = v.y * B.z - v.z * B.y;
    ctype const c1 = -v.x * B.z + v.z * B.x;
    ctype const c2 = v.x * B.y - v.y * B.x;
    return BiVec3d<ctype>(c0, c1, c2);
}

} // namespace hd::ga::ega