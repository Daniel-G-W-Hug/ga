#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_pga3dp_ops_basics.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga3dp product operations:
//
// - dot()                   -> dot product
// - rdot()                  -> regressive dot product
// - wdg(), join()           -> wedge product (join as convenience interface)
// - rwdg(), meet()          -> regressive wedge product (meet as convenience interface)
// - operator<<()            -> left contraction (= left bulk contraction)
// - operator>>()            -> right contraction (= right bulk contraction)
// - cmt()                   -> commutator product (= asymmetric part of gpr)
// - rcmt()                  -> regressive commutator product (=asymmetric part of rgpr)
// - operator*()             -> geometric product (= gpr)
// - rgpr()                  -> regressive geometric product
// - inv()                   -> inversion operation (w.r.t. geometric product)
//
// - left_bulk_contract3dp()    -> left bulk contraction
// - left_weight_contract3dp()  -> left weight contraction
// - right_bulk_contract3dp()   -> right bulk contraction
// - right_weight_contract3dp() -> right weight contraction
//
// - left_bulk_expand3dp()      -> left bulk expansion
// - left_weight_expand3dp()    -> left weight expansion
// - right_bulk_expand3dp()     -> right bulk expansion
// - right_weight_expand3dp()   -> right weight expansion
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot products for 3dp (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = +1, e4^2 = 0
//
////////////////////////////////////////////////////////////////////////////////

// returns 0.0 due to degenerate metric with e4^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(t1.w * t2.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(BiVec3dp<T> const& B1,
                                                  BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(B1.mx * B2.mx + B1.my * B2.my + B1.mz * B2.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(Vec3dp<T> const& v1,
                                                  Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}

// scalar product dot(a,b) = nrm_sq(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(MVec3dp<T> const& A,
                                                  MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // return Scalar3dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 +
    //                         A.c8 * B.c8 + A.c9 * B.c9 + A.c10 * B.c10 + A.c14 * B.c14);
    return Scalar3dp<ctype>(dot(gr0(A), gr0(B)) + dot(gr1(A), gr1(B)) +
                            dot(gr2(A), gr2(B)) + dot(gr3(A), gr3(B)) +
                            dot(gr4(A), gr4(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(MVec3dp_E<T> const& A,
                                                  MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(dot(gr0(A), gr0(B)) + dot(gr2(A), gr2(B)) +
                            dot(gr4(A), gr4(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> dot(MVec3dp_U<T> const& A,
                                                  MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(dot(gr1(A), gr1(B)) + dot(gr3(A), gr3(B)));
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
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(MVec3dp<T> const& A,
                                                    MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7 +
                             A.c11 * B.c11 + A.c12 * B.c12 + A.c13 * B.c13 +
                             A.c15 * B.c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(PScalar3dp<T> ps1, PScalar3dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps1) * ctype(ps2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(TriVec3dp<T> const& t1,
                                                    TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(t1.x * t2.x + t1.y * t2.y + t1.z * t2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(BiVec3dp<T> const& B1,
                                                    BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(B1.vx * B2.vx + B1.vy * B2.vy + B1.vz * B2.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(Vec3dp<T> const& v1,
                                                    Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(v1.w * v2.w);
}

// returns 0.0 due to degenerate metric with e4^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot([[maybe_unused]] Scalar3dp<T>,
                                                    [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(MVec3dp_E<T> const& A,
                                                    MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(rdot(gr0(A), gr0(B)) + rdot(gr2(A), gr2(B)) +
                             rdot(gr4(A), gr4(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rdot(MVec3dp_U<T> const& A,
                                                    MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(rdot(gr1(A), gr1(B)) + rdot(gr3(A), gr3(B)));
}


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////

// wedge product extended to fully populated multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> wdg(MVec3dp<T> const& A, MVec3dp<U> const& B)
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
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dp<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3dp<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3dp<T>,
                                                  [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec3dp<T> const&,
                                                  [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * ctype(s));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) * ctype(ps));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dp<T> const&,
                                                  [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] TriVec3dp<T> const&,
                                                  [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec3dp<T> const&,
                                                  [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(v.x * t.x + v.y * t.y + v.z * t.z + v.w * t.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(TriVec3dp<T> const& t, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return t * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * t;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B1,
                                                   BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(-B1.vx * B2.mx - B1.vy * B2.my - B1.vz * B2.mz -
                             B1.mx * B2.vx - B1.my * B2.vy - B1.mz * B2.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B,
                                                  Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(
        B.vy * v.z - B.vz * v.y + B.mx * v.w, -B.vx * v.z + B.vz * v.x + B.my * v.w,
        B.vx * v.y - B.vy * v.x + B.mz * v.w, -B.mx * v.x - B.my * v.y - B.mz * v.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v,
                                                  BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(
        -v.y * B.vz + v.z * B.vy + v.w * B.mx, v.x * B.vz - v.z * B.vx + v.w * B.my,
        -v.x * B.vy + v.y * B.vx + v.w * B.mz, -v.x * B.mx - v.y * B.my - v.z * B.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> wdg(BiVec3dp<T> const& B, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, BiVec3dp<U> const& B)
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
constexpr BiVec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    return BiVec3dp<std::common_type_t<T, U>>(
        v1.w * v2.x - v1.x * v2.w, v1.w * v2.y - v1.y * v2.w, v1.w * v2.z - v1.z * v2.w,
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> wdg(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> wdg(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg for Point3d
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Plane3d<std::common_type_t<T, U>> wdg(Line3d<T> const& l, Point3d<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return Plane3d<ctype>(l.vy * p.z - l.vz * p.y + l.mx, -l.vx * p.z + l.vz * p.x + l.my,
                          l.vx * p.y - l.vy * p.x + l.mz,
                          -l.mx * p.x - l.my * p.y - l.mz * p.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Plane3d<std::common_type_t<T, U>> wdg(Point3d<T> const& p, Line3d<U> const& l)
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
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> join(TriVec3dp<T> const& t,
                                                    Vec3dp<U> const& v)
{
    return wdg(t, v);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v,
                                                    TriVec3dp<U> const& t)
{
    return wdg(v, t);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> join(BiVec3dp<T> const& B,
                                                   Vec3dp<U> const& v)
{
    return wdg(B, v);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> join(Vec3dp<T> const& v,
                                                   BiVec3dp<U> const& B)
{
    return wdg(v, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Plane3d<std::common_type_t<T, U>> join(Line3d<T> const& l, Point3d<U> const& p)
{
    return wdg(l, p);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Plane3d<std::common_type_t<T, U>> join(Point3d<T> const& p, Line3d<U> const& l)
{
    return wdg(p, l);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Line3d<std::common_type_t<T, U>> join(Point3d<T> const& p, Point3d<U> const& q)
{
    return wdg(p, q);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
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
//        rwdg(ul, ur) = lcmpl(wdg(rcmpl(ul),rcmpl(ur))) = lcmpl(rcmpl(ul) ^ rcmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rwdg(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
               A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 - A.c8 * B.c5 - A.c9 * B.c6 -
               A.c10 * B.c7 - A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 +
               A.c15 * B.c0;
    ctype c1 = A.c1 * B.c15 + A.c5 * B.c14 + A.c9 * B.c13 - A.c10 * B.c12 -
               A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 + A.c15 * B.c1;
    ctype c2 = A.c2 * B.c15 + A.c6 * B.c14 - A.c8 * B.c13 + A.c10 * B.c11 +
               A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c6 + A.c15 * B.c2;
    ctype c3 = A.c3 * B.c15 + A.c7 * B.c14 + A.c8 * B.c12 - A.c9 * B.c11 - A.c11 * B.c9 +
               A.c12 * B.c8 + A.c14 * B.c7 + A.c15 * B.c3;
    ctype c4 = A.c4 * B.c15 - A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 - A.c11 * B.c5 -
               A.c12 * B.c6 - A.c13 * B.c7 + A.c15 * B.c4;
    ctype c5 = A.c5 * B.c15 - A.c12 * B.c13 + A.c13 * B.c12 + A.c15 * B.c5;
    ctype c6 = A.c6 * B.c15 + A.c11 * B.c13 - A.c13 * B.c11 + A.c15 * B.c6;
    ctype c7 = A.c7 * B.c15 - A.c11 * B.c12 + A.c12 * B.c11 + A.c15 * B.c7;
    ctype c8 = A.c8 * B.c15 + A.c11 * B.c14 - A.c14 * B.c11 + A.c15 * B.c8;
    ctype c9 = A.c9 * B.c15 + A.c12 * B.c14 - A.c14 * B.c12 + A.c15 * B.c9;
    ctype c10 = A.c10 * B.c15 + A.c13 * B.c14 - A.c14 * B.c13 + A.c15 * B.c10;
    ctype c11 = A.c11 * B.c15 + A.c15 * B.c11;
    ctype c12 = A.c12 * B.c15 + A.c15 * B.c12;
    ctype c13 = A.c13 * B.c15 + A.c15 * B.c13;
    ctype c14 = A.c14 * B.c15 + A.c15 * B.c14;
    ctype c15 = A.c15 * B.c15;
    return MVec3dp<std::common_type_t<T, U>>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,
                                             c11, c12, c13, c14, c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps1, PScalar3dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp(ctype(ps1) * ctype(ps2)); // convert to ctype before product!
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rwdg(PScalar3dp<T> ps,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * t;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                   PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return t * ctype(ps);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    return BiVec3dp<std::common_type_t<T, U>>(
        t1.z * t2.y - t1.y * t2.z, t1.x * t2.z - t1.z * t2.x, t1.y * t2.x - t1.x * t2.y,
        t1.x * t2.w - t1.w * t2.x, t1.y * t2.w - t1.w * t2.y, t1.z * t2.w - t1.w * t2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                BiVec3dp<U> const& B)
{
    return Vec3dp<std::common_type_t<T, U>>(
        -t.y * B.mz + t.z * B.my + t.w * B.vx, t.x * B.mz - t.z * B.mx + t.w * B.vy,
        -t.x * B.my + t.y * B.mx + t.w * B.vz, -t.x * B.vx - t.y * B.vy - t.z * B.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B,
                                                TriVec3dp<U> const& t)
{
    return Vec3dp<std::common_type_t<T, U>>(
        B.vx * t.w + B.my * t.z - B.mz * t.y, B.vy * t.w - B.mx * t.z + B.mz * t.x,
        B.vz * t.w + B.mx * t.y - B.my * t.x, -B.vx * t.x - B.vy * t.y - B.vz * t.z);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    return Scalar3dp<std::common_type_t<T, U>>(-t.x * v.x - t.y * v.y - t.z * v.z -
                                               t.w * v.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(Vec3dp<T> const& v,
                                                   TriVec3dp<U> const& t)
{
    return Scalar3dp<std::common_type_t<T, U>>(v.x * t.x + v.y * t.y + v.z * t.z +
                                               v.w * t.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] TriVec3dp<T> const&,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    return Scalar3dp<std::common_type_t<T, U>>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] TriVec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>{0.0};
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg(BiVec3dp<T> const& B1,
                                                   BiVec3dp<U> const& B2)
{
    return Scalar3dp<std::common_type_t<T, U>>(-B1.vx * B2.mx - B1.vy * B2.my -
                                               B1.vz * B2.mz - B1.mx * B2.vx -
                                               B1.my * B2.vy - B1.mz * B2.vz);
}

// required to be present for dist3dp (to complile, even if not used)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] BiVec3dp<T> const&,
                                                   [[maybe_unused]] Vec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dp<T> const&,
                                                   [[maybe_unused]] BiVec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dp<T> const&,
                                                   [[maybe_unused]] Vec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec3dp<T> const&,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] Vec3dp<U> const&)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar3dp<T>,
                                                   [[maybe_unused]] Scalar3dp<U>)
{
    return Scalar3dp<std::common_type_t<T, U>>(0.0);
}

////////////////////////////////////////////////////////////////////////////////
// convenience functions rwdg -> meet
////////////////////////////////////////////////////////////////////////////////

// intersection of two planes -> return a line
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> meet(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    return rwdg(t1, t2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Line3d<std::common_type_t<T, U>> meet(Plane3d<T> const& p1,
                                                Plane3d<U> const& p2)
{
    return rwdg(static_cast<TriVec3dp<T> const&>(p1),
                static_cast<TriVec3dp<T> const&>(p2));
}

// intersection of plane and line -> return a point
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(TriVec3dp<T> const& t,
                                                BiVec3dp<U> const& B)
{
    return rwdg(t, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(BiVec3dp<T> const& B,
                                                TriVec3dp<U> const& t)
{
    return rwdg(B, t);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(Plane3d<T> const& p, Line3d<U> const& l)
{
    return rwdg(static_cast<TriVec3dp<T> const&>(p), static_cast<BiVec3dp<T> const&>(l));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> meet(Line3d<T> const& l, Plane3d<U> const& p)
{
    return rwdg(static_cast<BiVec3dp<T> const&>(l), static_cast<TriVec3dp<T> const&>(p));
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// Implements the left bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator<<(a,b) = left_bulk_contract(a,b) = rwdg( left_bulk_dual(a), b )
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator<<(MVec3dp<T> const& A,
                                                       MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c8 * B.c8 +
               A.c9 * B.c9 + A.c10 * B.c10 + A.c14 * B.c14;
    ctype c1 = A.c0 * B.c1 + A.c2 * B.c10 - A.c3 * B.c9 - A.c8 * B.c14;
    ctype c2 = A.c0 * B.c2 - A.c1 * B.c10 + A.c3 * B.c8 - A.c9 * B.c14;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c9 - A.c2 * B.c8 - A.c10 * B.c14;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c7 + A.c8 * B.c11 +
               A.c9 * B.c12 + A.c10 * B.c13 + A.c14 * B.c15;
    ctype c5 = A.c0 * B.c5 + A.c2 * B.c13 - A.c3 * B.c12 - A.c8 * B.c15;
    ctype c6 = A.c0 * B.c6 - A.c1 * B.c13 + A.c3 * B.c11 - A.c9 * B.c15;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c12 - A.c2 * B.c11 - A.c10 * B.c15;
    ctype c8 = A.c0 * B.c8 - A.c1 * B.c14;
    ctype c9 = A.c0 * B.c9 - A.c2 * B.c14;
    ctype c10 = A.c0 * B.c10 - A.c3 * B.c14;
    ctype c11 = A.c0 * B.c11 - A.c1 * B.c15;
    ctype c12 = A.c0 * B.c12 - A.c2 * B.c15;
    ctype c13 = A.c0 * B.c13 - A.c3 * B.c15;
    ctype c14 = A.c0 * B.c14;
    ctype c15 = A.c0 * B.c15;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] TriVec3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(TriVec3dp<T> const& t,
                                                      PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(0.0, 0.0, 0.0, t.w * ctype(ps));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dp<T>, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-B.mx, -B.my, -B.mz, 0.0, 0.0, 0.0) * ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3dp<T>, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                         PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(-v.x, -v.y, -v.z, 0.0) * ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s,
                                                          PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) * ctype(ps));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    return dot(t1, t2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B,
                                                      TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(-B.mx * t.w, -B.my * t.w, -B.mz * t.w,
                         B.mx * t.x + B.my * t.y + B.mz * t.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(v.y * t.z - v.z * t.y, -v.x * t.z + v.z * t.x,
                           v.x * t.y - v.y * t.x, -v.x * t.w, -v.y * t.w, -v.z * t.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s,
                                                         TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * t;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<(BiVec3dp<T> const& B1,
                                                         BiVec3dp<U> const& B2)
{
    return dot(B1, B2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

// identical to cmt(B,v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v,
                                                      BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(v.y * B.mz - v.z * B.my, -v.x * B.mz + v.z * B.mx,
                         v.x * B.my - v.y * B.mx, v.x * B.vx + v.y * B.vy + v.z * B.vz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator<<(Vec3dp<T> const& v1,
                                                         Vec3dp<U> const& v2)
{
    return dot(v1, v2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] Scalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator<<(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
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
// operator>>(a,b) = right_bulk_contract(a,b) = rwdg( a, right_bulk_dual(b) )
//
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator>>(MVec3dp<T> const& A,
                                                       MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c8 * B.c8 +
               A.c9 * B.c9 + A.c10 * B.c10 + A.c14 * B.c14;
    ctype c1 = A.c1 * B.c0 + A.c9 * B.c3 - A.c10 * B.c2 - A.c14 * B.c8;
    ctype c2 = A.c2 * B.c0 - A.c8 * B.c3 + A.c10 * B.c1 - A.c14 * B.c9;
    ctype c3 = A.c3 * B.c0 + A.c8 * B.c2 - A.c9 * B.c1 - A.c14 * B.c10;
    ctype c4 = A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 + A.c11 * B.c8 +
               A.c12 * B.c9 + A.c13 * B.c10 - A.c15 * B.c14;
    ctype c5 = A.c5 * B.c0 - A.c12 * B.c3 + A.c13 * B.c2 - A.c15 * B.c8;
    ctype c6 = A.c6 * B.c0 + A.c11 * B.c3 - A.c13 * B.c1 - A.c15 * B.c9;
    ctype c7 = A.c7 * B.c0 - A.c11 * B.c2 + A.c12 * B.c1 - A.c15 * B.c10;
    ctype c8 = A.c8 * B.c0 - A.c14 * B.c1;
    ctype c9 = A.c9 * B.c0 - A.c14 * B.c2;
    ctype c10 = A.c10 * B.c0 - A.c14 * B.c3;
    ctype c11 = A.c11 * B.c0 + A.c15 * B.c1;
    ctype c12 = A.c12 * B.c0 + A.c15 * B.c2;
    ctype c13 = A.c13 * B.c0 + A.c15 * B.c3;
    ctype c14 = A.c14 * B.c0;
    ctype c15 = A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] PScalar3dp<T>,
                                                         [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                      TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(0.0, 0.0, 0.0, -ctype(ps) * t.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] TriVec3dp<T> const&, [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * BiVec3dp<ctype>(-B.mx, -B.my, -B.mz, 0.0, 0.0, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                         Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * TriVec3dp<ctype>(v.x, v.y, v.z, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator>>(PScalar3dp<T> ps,
                                                          Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) * ctype(s));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar3dp<T>,
                                                         [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    return dot(t1, t2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                      BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(-t.w * B.mx, -t.w * B.my, -t.w * B.mz,
                         t.x * B.mx + t.y * B.my + t.z * B.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3dp<T> const&, [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-t.y * v.z + t.z * v.y, t.x * v.z - t.z * v.x,
                           -t.x * v.y + t.y * v.x, -t.w * v.x, -t.w * v.y, -t.w * v.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator>>(TriVec3dp<T> const& t,
                                                         Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return t * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dp<T>, [[maybe_unused]] TriVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(BiVec3dp<T> const& B1,
                                                         BiVec3dp<U> const& B2)
{
    return dot(B1, B2);
}

// identical to cmt(v,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(BiVec3dp<T> const& B,
                                                      Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(B.my * v.z - B.mz * v.y, -B.mx * v.z + B.mz * v.x,
                         B.mx * v.y - B.my * v.x, -B.vx * v.x - B.vy * v.y - B.vz * v.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3dp<T> const&, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator>>(BiVec3dp<T> const& B,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dp<T>, [[maybe_unused]] BiVec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(Vec3dp<T> const& v1,
                                                         Vec3dp<U> const& v2)
{
    return dot(v1, v2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator>>(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3dp<T>, [[maybe_unused]] Vec3dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator>>(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> cmt(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(0.0, -A.c2 * B.c10 + A.c3 * B.c9 - A.c9 * B.c3 + A.c10 * B.c2,
                          A.c1 * B.c10 - A.c3 * B.c8 + A.c8 * B.c3 - A.c10 * B.c1,
                          -A.c1 * B.c9 + A.c2 * B.c8 - A.c8 * B.c2 + A.c9 * B.c1,
                          -A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 + A.c5 * B.c1 +
                              A.c6 * B.c2 + A.c7 * B.c3 + A.c14 * B.c15 - A.c15 * B.c14,
                          -A.c1 * B.c4 + A.c4 * B.c1 - A.c6 * B.c10 + A.c7 * B.c9 -
                              A.c9 * B.c7 + A.c10 * B.c6 + A.c11 * B.c14 - A.c14 * B.c11,
                          -A.c2 * B.c4 + A.c4 * B.c2 + A.c5 * B.c10 - A.c7 * B.c8 +
                              A.c8 * B.c7 - A.c10 * B.c5 + A.c12 * B.c14 - A.c14 * B.c12,
                          -A.c3 * B.c4 + A.c4 * B.c3 - A.c5 * B.c9 + A.c6 * B.c8 -
                              A.c8 * B.c6 + A.c9 * B.c5 + A.c13 * B.c14 - A.c14 * B.c13,
                          A.c2 * B.c3 - A.c3 * B.c2 - A.c9 * B.c10 + A.c10 * B.c9,
                          -A.c1 * B.c3 + A.c3 * B.c1 + A.c8 * B.c10 - A.c10 * B.c8,
                          A.c1 * B.c2 - A.c2 * B.c1 - A.c8 * B.c9 + A.c9 * B.c8,
                          A.c1 * B.c15 - A.c5 * B.c14 - A.c9 * B.c13 + A.c10 * B.c12 -
                              A.c12 * B.c10 + A.c13 * B.c9 + A.c14 * B.c5 - A.c15 * B.c1,
                          A.c2 * B.c15 - A.c6 * B.c14 + A.c8 * B.c13 - A.c10 * B.c11 +
                              A.c11 * B.c10 - A.c13 * B.c8 + A.c14 * B.c6 - A.c15 * B.c2,
                          A.c3 * B.c15 - A.c7 * B.c14 - A.c8 * B.c12 + A.c9 * B.c11 -
                              A.c11 * B.c9 + A.c12 * B.c8 + A.c14 * B.c7 - A.c15 * B.c3,
                          0.0,
                          A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
                              A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> cmt(TriVec3dp<T> const& t1,
                                                 TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(t1.x * t2.w - t1.w * t2.x, t1.y * t2.w - t1.w * t2.y,
                           t1.z * t2.w - t1.w * t2.z, 0.0, 0.0, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> cmt(TriVec3dp<T> const& t,
                                                  BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(-t.y * B.mz + t.z * B.my + t.w * B.vx,
                            t.x * B.mz - t.z * B.mx + t.w * B.vy,
                            -t.x * B.my + t.y * B.mx + t.w * B.vz, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& B,
                                                  TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(-B.vx * t.w - B.my * t.z + B.mz * t.y,
                            -B.vy * t.w + B.mx * t.z - B.mz * t.x,
                            -B.vz * t.w - B.mx * t.y + B.my * t.x, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& B1,
                                                 BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-B1.vy * B2.mz + B1.vz * B2.my - B1.my * B2.vz + B1.mz * B2.vy,
                           B1.vx * B2.mz - B1.vz * B2.mx + B1.mx * B2.vz - B1.mz * B2.vx,
                           -B1.vx * B2.my + B1.vy * B2.mx - B1.mx * B2.vy + B1.my * B2.vx,
                           -B1.my * B2.mz + B1.mz * B2.my, B1.mx * B2.mz - B1.mz * B2.mx,
                           -B1.mx * B2.my + B1.my * B2.mx);
}

// cmt(B,v) = -cmt(v,B)
// identical to (v << B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> cmt(BiVec3dp<T> const& B, Vec3dp<U> const& v)
{
    return (v << B);
}

// cmt(v,B) = -cmt(B,v)
// identical to (B >> v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> cmt(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    return (B >> v);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> cmt(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-v1.x * v2.w + v1.w * v2.x, -v1.y * v2.w + v1.w * v2.y,
                           -v1.z * v2.w + v1.w * v2.z, v1.y * v2.z - v1.z * v2.y,
                           -v1.x * v2.z + v1.z * v2.x, v1.x * v2.y - v1.y * v2.x);
}


////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rcmt(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 - A.c11 * B.c1 -
               A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4;
    ctype c1 = A.c0 * B.c11 + A.c2 * B.c7 - A.c3 * B.c6 - A.c4 * B.c8 + A.c6 * B.c3 -
               A.c7 * B.c2 + A.c8 * B.c4 - A.c11 * B.c0;
    ctype c2 = A.c0 * B.c12 - A.c1 * B.c7 + A.c3 * B.c5 - A.c4 * B.c9 - A.c5 * B.c3 +
               A.c7 * B.c1 + A.c9 * B.c4 - A.c12 * B.c0;
    ctype c3 = A.c0 * B.c13 + A.c1 * B.c6 - A.c2 * B.c5 - A.c4 * B.c10 + A.c5 * B.c2 -
               A.c6 * B.c1 + A.c10 * B.c4 - A.c13 * B.c0;
    ctype c4 = 0.0;
    ctype c5 = A.c6 * B.c7 - A.c7 * B.c6 - A.c12 * B.c13 + A.c13 * B.c12;
    ctype c6 = -A.c5 * B.c7 + A.c7 * B.c5 + A.c11 * B.c13 - A.c13 * B.c11;
    ctype c7 = A.c5 * B.c6 - A.c6 * B.c5 - A.c11 * B.c12 + A.c12 * B.c11;
    ctype c8 = -A.c1 * B.c4 + A.c4 * B.c1 + A.c6 * B.c10 - A.c7 * B.c9 + A.c9 * B.c7 -
               A.c10 * B.c6 + A.c11 * B.c14 - A.c14 * B.c11;
    ctype c9 = -A.c2 * B.c4 + A.c4 * B.c2 - A.c5 * B.c10 + A.c7 * B.c8 - A.c8 * B.c7 +
               A.c10 * B.c5 + A.c12 * B.c14 - A.c14 * B.c12;
    ctype c10 = -A.c3 * B.c4 + A.c4 * B.c3 + A.c5 * B.c9 - A.c6 * B.c8 + A.c8 * B.c6 -
                A.c9 * B.c5 + A.c13 * B.c14 - A.c14 * B.c13;
    ctype c11 = A.c6 * B.c13 - A.c7 * B.c12 + A.c12 * B.c7 - A.c13 * B.c6;
    ctype c12 = -A.c5 * B.c13 + A.c7 * B.c11 - A.c11 * B.c7 + A.c13 * B.c5;
    ctype c13 = A.c5 * B.c12 - A.c6 * B.c11 + A.c11 * B.c6 - A.c12 * B.c5;
    ctype c14 = A.c0 * B.c4 - A.c4 * B.c0 - A.c8 * B.c11 - A.c9 * B.c12 - A.c10 * B.c13 +
                A.c11 * B.c8 + A.c12 * B.c9 + A.c13 * B.c10;
    ctype c15 = 0.0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rcmt(TriVec3dp<T> const& t1,
                                                  TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-t1.y * t2.z + t1.z * t2.y, t1.x * t2.z - t1.z * t2.x,
                           -t1.x * t2.y + t1.y * t2.x, t1.x * t2.w - t1.w * t2.x,
                           t1.y * t2.w - t1.w * t2.y, t1.z * t2.w - t1.w * t2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rcmt(TriVec3dp<T> const& t,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(t.y * B.vz - t.z * B.vy, -t.x * B.vz + t.z * B.vx,
                            t.x * B.vy - t.y * B.vx,
                            t.x * B.mx + t.y * B.my + t.z * B.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> rcmt(BiVec3dp<T> const& B,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(B.vy * t.z - B.vz * t.y, -B.vx * t.z + B.vz * t.x,
                            B.vx * t.y - B.vy * t.x,
                            -B.mx * t.x - B.my * t.y - B.mz * t.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rcmt(BiVec3dp<T> const& B1,
                                                  BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(B1.vy * B2.vz - B1.vz * B2.vy, -B1.vx * B2.vz + B1.vz * B2.vx,
                           B1.vx * B2.vy - B1.vy * B2.vx,
                           B1.vy * B2.mz - B1.vz * B2.my + B1.my * B2.vz - B1.mz * B2.vy,
                           -B1.vx * B2.mz + B1.vz * B2.mx - B1.mx * B2.vz + B1.mz * B2.vx,
                           B1.vx * B2.my - B1.vy * B2.mx + B1.mx * B2.vy - B1.my * B2.vx);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rcmt(BiVec3dp<T> const& B, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(B.vy * v.z - B.vz * v.y + B.mx * v.w,
                         -B.vx * v.z + B.vz * v.x + B.my * v.w,
                         B.vx * v.y - B.vy * v.x + B.mz * v.w, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> rcmt(Vec3dp<T> const& v, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(v.y * B.vz - v.z * B.vy - v.w * B.mx,
                         -v.x * B.vz + v.z * B.vx - v.w * B.my,
                         v.x * B.vy - v.y * B.vx - v.w * B.mz, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> rcmt(Vec3dp<T> const& v1,
                                                  Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(0.0, 0.0, 0.0, -v1.x * v2.w + v1.w * v2.x,
                           -v1.y * v2.w + v1.w * v2.y, -v1.z * v2.w + v1.w * v2.z);
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
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
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

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                      MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c8 * B.c4 - A.c9 * B.c5 - A.c10 * B.c6;
    ctype c1 = A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 + A.c14 * B.c4;
    ctype c2 = A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c14 * B.c5;
    ctype c3 = -A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 + A.c14 * B.c6;
    ctype c4 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c4 * B.c0 - A.c11 * B.c4 -
               A.c12 * B.c5 - A.c13 * B.c6 + A.c14 * B.c7;
    ctype c5 = A.c0 * B.c1 + A.c5 * B.c0 - A.c6 * B.c6 + A.c7 * B.c5 + A.c8 * B.c7 -
               A.c9 * B.c3 + A.c10 * B.c2 + A.c15 * B.c4;
    ctype c6 = A.c0 * B.c2 + A.c5 * B.c6 + A.c6 * B.c0 - A.c7 * B.c4 + A.c8 * B.c3 +
               A.c9 * B.c7 - A.c10 * B.c1 + A.c15 * B.c5;
    ctype c7 = A.c0 * B.c3 - A.c5 * B.c5 + A.c6 * B.c4 + A.c7 * B.c0 - A.c8 * B.c2 +
               A.c9 * B.c1 + A.c10 * B.c7 + A.c15 * B.c6;
    ctype c8 = A.c0 * B.c4 + A.c8 * B.c0 - A.c9 * B.c6 + A.c10 * B.c5;
    ctype c9 = A.c0 * B.c5 + A.c8 * B.c6 + A.c9 * B.c0 - A.c10 * B.c4;
    ctype c10 = A.c0 * B.c6 - A.c8 * B.c5 + A.c9 * B.c4 + A.c10 * B.c0;
    ctype c11 = A.c1 * B.c7 - A.c2 * B.c3 + A.c3 * B.c2 + A.c4 * B.c4 + A.c11 * B.c0 -
                A.c12 * B.c6 + A.c13 * B.c5 + A.c14 * B.c1;
    ctype c12 = A.c1 * B.c3 + A.c2 * B.c7 - A.c3 * B.c1 + A.c4 * B.c5 + A.c11 * B.c6 +
                A.c12 * B.c0 - A.c13 * B.c4 + A.c14 * B.c2;
    ctype c13 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c7 + A.c4 * B.c6 - A.c11 * B.c5 +
                A.c12 * B.c4 + A.c13 * B.c0 + A.c14 * B.c3;
    ctype c14 = -A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 + A.c14 * B.c0;
    ctype c15 = A.c0 * B.c7 - A.c5 * B.c4 - A.c6 * B.c5 - A.c7 * B.c6 - A.c8 * B.c1 -
                A.c9 * B.c2 - A.c10 * B.c3 + A.c15 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
                                                      MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c4 * B.c8 - A.c5 * B.c9 - A.c6 * B.c10;
    ctype c1 = A.c0 * B.c1 + A.c4 * B.c14 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c14 - A.c6 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c14;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c11 -
               A.c5 * B.c12 - A.c6 * B.c13 - A.c7 * B.c14;
    ctype c5 = A.c0 * B.c5 + A.c1 * B.c0 - A.c2 * B.c10 + A.c3 * B.c9 + A.c4 * B.c15 -
               A.c5 * B.c7 + A.c6 * B.c6 + A.c7 * B.c8;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c10 + A.c2 * B.c0 - A.c3 * B.c8 + A.c4 * B.c7 +
               A.c5 * B.c15 - A.c6 * B.c5 + A.c7 * B.c9;
    ctype c7 = A.c0 * B.c7 - A.c1 * B.c9 + A.c2 * B.c8 + A.c3 * B.c0 - A.c4 * B.c6 +
               A.c5 * B.c5 + A.c6 * B.c15 + A.c7 * B.c10;
    ctype c8 = A.c0 * B.c8 + A.c4 * B.c0 - A.c5 * B.c10 + A.c6 * B.c9;
    ctype c9 = A.c0 * B.c9 + A.c4 * B.c10 + A.c5 * B.c0 - A.c6 * B.c8;
    ctype c10 = A.c0 * B.c10 - A.c4 * B.c9 + A.c5 * B.c8 + A.c6 * B.c0;
    ctype c11 = A.c0 * B.c11 - A.c1 * B.c14 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c4 -
                A.c5 * B.c13 + A.c6 * B.c12 - A.c7 * B.c1;
    ctype c12 = A.c0 * B.c12 - A.c1 * B.c3 - A.c2 * B.c14 + A.c3 * B.c1 + A.c4 * B.c13 +
                A.c5 * B.c4 - A.c6 * B.c11 - A.c7 * B.c2;
    ctype c13 = A.c0 * B.c13 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c14 - A.c4 * B.c12 +
                A.c5 * B.c11 + A.c6 * B.c4 - A.c7 * B.c3;
    ctype c14 = A.c0 * B.c14 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype c15 = A.c0 * B.c15 - A.c1 * B.c8 - A.c2 * B.c9 - A.c3 * B.c10 - A.c4 * B.c5 -
                A.c5 * B.c6 - A.c6 * B.c7 + A.c7 * B.c0;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                      MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c14 * B.c7;
    ctype c1 = A.c0 * B.c0 + A.c8 * B.c7 - A.c9 * B.c2 + A.c10 * B.c1;
    ctype c2 = A.c0 * B.c1 + A.c8 * B.c2 + A.c9 * B.c7 - A.c10 * B.c0;
    ctype c3 = A.c0 * B.c2 - A.c8 * B.c1 + A.c9 * B.c0 + A.c10 * B.c7;
    ctype c4 = A.c0 * B.c3 + A.c5 * B.c0 + A.c6 * B.c1 + A.c7 * B.c2 - A.c8 * B.c4 -
               A.c9 * B.c5 - A.c10 * B.c6 - A.c15 * B.c7;
    ctype c5 = -A.c1 * B.c3 + A.c2 * B.c6 - A.c3 * B.c5 + A.c4 * B.c0 + A.c11 * B.c7 -
               A.c12 * B.c2 + A.c13 * B.c1 - A.c14 * B.c4;
    ctype c6 = -A.c1 * B.c6 - A.c2 * B.c3 + A.c3 * B.c4 + A.c4 * B.c1 + A.c11 * B.c2 +
               A.c12 * B.c7 - A.c13 * B.c0 - A.c14 * B.c5;
    ctype c7 = A.c1 * B.c5 - A.c2 * B.c4 - A.c3 * B.c3 + A.c4 * B.c2 - A.c11 * B.c1 +
               A.c12 * B.c0 + A.c13 * B.c7 - A.c14 * B.c6;
    ctype c8 = -A.c1 * B.c7 + A.c2 * B.c2 - A.c3 * B.c1 - A.c14 * B.c0;
    ctype c9 = -A.c1 * B.c2 - A.c2 * B.c7 + A.c3 * B.c0 - A.c14 * B.c1;
    ctype c10 = A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c7 - A.c14 * B.c2;
    ctype c11 = A.c0 * B.c4 - A.c5 * B.c7 + A.c6 * B.c2 - A.c7 * B.c1 + A.c8 * B.c3 -
                A.c9 * B.c6 + A.c10 * B.c5 - A.c15 * B.c0;
    ctype c12 = A.c0 * B.c5 - A.c5 * B.c2 - A.c6 * B.c7 + A.c7 * B.c0 + A.c8 * B.c6 +
                A.c9 * B.c3 - A.c10 * B.c4 - A.c15 * B.c1;
    ctype c13 = A.c0 * B.c6 + A.c5 * B.c1 - A.c6 * B.c0 - A.c7 * B.c7 - A.c8 * B.c5 +
                A.c9 * B.c4 + A.c10 * B.c3 - A.c15 * B.c2;
    ctype c14 = A.c0 * B.c7 - A.c8 * B.c0 - A.c9 * B.c1 - A.c10 * B.c2;
    ctype c15 = A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c7 - A.c11 * B.c0 -
                A.c12 * B.c1 - A.c13 * B.c2 - A.c14 * B.c3;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
                                                      MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c7 * B.c14;
    ctype c1 = A.c0 * B.c0 - A.c1 * B.c10 + A.c2 * B.c9 + A.c7 * B.c8;
    ctype c2 = A.c0 * B.c10 + A.c1 * B.c0 - A.c2 * B.c8 + A.c7 * B.c9;
    ctype c3 = -A.c0 * B.c9 + A.c1 * B.c8 + A.c2 * B.c0 + A.c7 * B.c10;
    ctype c4 = -A.c0 * B.c5 - A.c1 * B.c6 - A.c2 * B.c7 + A.c3 * B.c0 - A.c4 * B.c8 -
               A.c5 * B.c9 - A.c6 * B.c10 + A.c7 * B.c15;
    ctype c5 = -A.c0 * B.c4 + A.c1 * B.c13 - A.c2 * B.c12 + A.c3 * B.c1 + A.c4 * B.c14 -
               A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c11;
    ctype c6 = -A.c0 * B.c13 - A.c1 * B.c4 + A.c2 * B.c11 + A.c3 * B.c2 + A.c4 * B.c3 +
               A.c5 * B.c14 - A.c6 * B.c1 - A.c7 * B.c12;
    ctype c7 = A.c0 * B.c12 - A.c1 * B.c11 - A.c2 * B.c4 + A.c3 * B.c3 - A.c4 * B.c2 +
               A.c5 * B.c1 + A.c6 * B.c14 - A.c7 * B.c13;
    ctype c8 = -A.c0 * B.c14 + A.c1 * B.c3 - A.c2 * B.c2 - A.c7 * B.c1;
    ctype c9 = -A.c0 * B.c3 - A.c1 * B.c14 + A.c2 * B.c1 - A.c7 * B.c2;
    ctype c10 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c14 - A.c7 * B.c3;
    ctype c11 = A.c0 * B.c15 - A.c1 * B.c7 + A.c2 * B.c6 + A.c3 * B.c8 + A.c4 * B.c0 -
                A.c5 * B.c10 + A.c6 * B.c9 + A.c7 * B.c5;
    ctype c12 = A.c0 * B.c7 + A.c1 * B.c15 - A.c2 * B.c5 + A.c3 * B.c9 + A.c4 * B.c10 +
                A.c5 * B.c0 - A.c6 * B.c8 + A.c7 * B.c6;
    ctype c13 = -A.c0 * B.c6 + A.c1 * B.c5 + A.c2 * B.c15 + A.c3 * B.c10 - A.c4 * B.c9 +
                A.c5 * B.c8 + A.c6 * B.c0 + A.c7 * B.c7;
    ctype c14 = -A.c0 * B.c8 - A.c1 * B.c9 - A.c2 * B.c10 + A.c7 * B.c0;
    ctype c15 = A.c0 * B.c11 + A.c1 * B.c12 + A.c2 * B.c13 + A.c3 * B.c14 - A.c4 * B.c1 -
                A.c5 * B.c2 - A.c6 * B.c3 - A.c7 * B.c4;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(MVec3dp<T> const& A,
                                                      PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), A.c14, A.c8,
                          A.c9, A.c10, ctype(0.0), ctype(0.0), ctype(0.0), A.c1, A.c2,
                          A.c3, ctype(0.0), A.c0) *
           ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                      MVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                                      -B.c14, B.c8, B.c9, B.c10, ctype(0.0), ctype(0.0),
                                      ctype(0.0), -B.c1, -B.c2, -B.c3, ctype(0.0), B.c0);
}


// geometric product A * B for two multivectors from the even subalgebra (3d case)
// even grade multivector * even grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
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

// geometric product A * B for two multivectors from the odd subalgebra (3d case)
// odd grade multivector * odd grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
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
// with a multivector B of the odd subalgebra
// even grade multivector * odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
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

// geometric product A * B of a multivector A from the odd subalgebra (3d case)
// with a multivector B of the even subalgebra
// odd grade multivector * even grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& A,
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

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<U> const& A,
                                                        PScalar3dp<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
               BiVec3dp<ctype>(A.c4, A.c5, A.c6, ctype(0.0), ctype(0.0), ctype(0.0)),
               PScalar3dp<ctype>(A.c0)) *
           ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                        MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec3dp_E<ctype>(BiVec3dp<ctype>(B.c4, B.c5, B.c6, ctype(0.0),
                                                        ctype(0.0), ctype(0.0)),
                                        PScalar3dp<ctype>(B.c0));
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<U> const& A,
                                                        PScalar3dp<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), A.c7),
                            TriVec3dp<ctype>(A.c0, A.c1, A.c2, ctype(0.0))) *
           ctype(ps);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                        MVec3dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) *
           MVec3dp_U<ctype>(Vec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), -B.c7),
                            TriVec3dp<ctype>(-B.c0, -B.c1, -B.c2, ctype(0.0)));
}


// geometric product M * B of a multivector M from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& M,
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

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(t.w * M.c4, t.w * M.c5, t.w * M.c6,
                      -t.x * M.c4 - t.y * M.c5 - t.z * M.c6 + t.w * M.c7),
        TriVec3dp<ctype>(t.x * M.c0 - t.y * M.c6 + t.z * M.c5 + t.w * M.c1,
                         t.x * M.c6 + t.y * M.c0 - t.z * M.c4 + t.w * M.c2,
                         -t.x * M.c5 + t.y * M.c4 + t.z * M.c0 + t.w * M.c3, t.w * M.c0));
}

// geometric product M * B of a multivector M from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& M,
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

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(-B.mx * M.c4 - B.my * M.c5 - B.mz * M.c6),
                            BiVec3dp<ctype>(B.vx * M.c0 - B.vy * M.c6 + B.vz * M.c5 +
                                                B.mx * M.c7 - B.my * M.c3 + B.mz * M.c2,
                                            B.vx * M.c6 + B.vy * M.c0 - B.vz * M.c4 +
                                                B.mx * M.c3 + B.my * M.c7 - B.mz * M.c1,
                                            -B.vx * M.c5 + B.vy * M.c4 + B.vz * M.c0 -
                                                B.mx * M.c2 + B.my * M.c1 + B.mz * M.c7,
                                            B.mx * M.c0 - B.my * M.c6 + B.mz * M.c5,
                                            B.mx * M.c6 + B.my * M.c0 - B.mz * M.c4,
                                            -B.mx * M.c5 + B.my * M.c4 + B.mz * M.c0),
                            PScalar3dp<ctype>(-B.vx * M.c4 - B.vy * M.c5 - B.vz * M.c6 -
                                              B.mx * M.c1 - B.my * M.c2 - B.mz * M.c3));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_U<T> const& M,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(-M.c1 * B.mz + M.c2 * B.my + M.c7 * B.mx,
                                          M.c0 * B.mz - M.c2 * B.mx + M.c7 * B.my,
                                          -M.c0 * B.my + M.c1 * B.mx + M.c7 * B.mz,
                                          -M.c0 * B.vx - M.c1 * B.vy - M.c2 * B.vz -
                                              M.c4 * B.mx - M.c5 * B.my - M.c6 * B.mz),
                            TriVec3dp<ctype>(-M.c1 * B.vz + M.c2 * B.vy + M.c3 * B.mx -
                                                 M.c5 * B.mz + M.c6 * B.my + M.c7 * B.vx,
                                             M.c0 * B.vz - M.c2 * B.vx + M.c3 * B.my +
                                                 M.c4 * B.mz - M.c6 * B.mx + M.c7 * B.vy,
                                             -M.c0 * B.vy + M.c1 * B.vx + M.c3 * B.mz -
                                                 M.c4 * B.my + M.c5 * B.mx + M.c7 * B.vz,
                                             -M.c0 * B.mx - M.c1 * B.my - M.c2 * B.mz));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        MVec3dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(B.mx * M.c7 - B.my * M.c2 + B.mz * M.c1,
                                          B.mx * M.c2 + B.my * M.c7 - B.mz * M.c0,
                                          -B.mx * M.c1 + B.my * M.c0 + B.mz * M.c7,
                                          B.vx * M.c0 + B.vy * M.c1 + B.vz * M.c2 -
                                              B.mx * M.c4 - B.my * M.c5 - B.mz * M.c6),
                            TriVec3dp<ctype>(-B.vx * M.c7 + B.vy * M.c2 - B.vz * M.c1 +
                                                 B.mx * M.c3 - B.my * M.c6 + B.mz * M.c5,
                                             -B.vx * M.c2 - B.vy * M.c7 + B.vz * M.c0 +
                                                 B.mx * M.c6 + B.my * M.c3 - B.mz * M.c4,
                                             B.vx * M.c1 - B.vy * M.c0 - B.vz * M.c7 -
                                                 B.mx * M.c5 + B.my * M.c4 + B.mz * M.c3,
                                             -B.mx * M.c0 - B.my * M.c1 - B.mz * M.c2));
}

// geometric product A * v of an even grade multivector A with a vector v
// from the right
// even grade multivector * vector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(MVec3dp_E<T> const& A,
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

// geometric product v * B of an even grade multivector B with a vector v
// from the left
// vector * even grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Vec3dp<U> const& v,
                                                        MVec3dp_E<T> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(v.x * B.c0 - v.y * B.c6 + v.z * B.c5,
                      v.x * B.c6 + v.y * B.c0 - v.z * B.c4,
                      -v.x * B.c5 + v.y * B.c4 + v.z * B.c0,
                      -v.x * B.c1 - v.y * B.c2 - v.z * B.c3 + v.w * B.c0),
        TriVec3dp<ctype>(v.x * B.c7 - v.y * B.c3 + v.z * B.c2 + v.w * B.c4,
                         v.x * B.c3 + v.y * B.c7 - v.z * B.c1 + v.w * B.c5,
                         -v.x * B.c2 + v.y * B.c1 + v.z * B.c7 + v.w * B.c6,
                         -v.x * B.c4 - v.y * B.c5 - v.z * B.c6));
}

// geometric product of two quadvectors (=3dp pseudoscalars)
// quadvector * quadvector => scalar == 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator*([[maybe_unused]] PScalar3dp<T>,
                                                        [[maybe_unused]] PScalar3dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(0.0); // pseudoscalar in pga3dp contains e4^2 = 0
}

// geometric product ps * t of a quadvector ps (=3dp pseudoscalar) multiplied from the
// left to the trivector t quadvector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                     TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), -t.w);
}

// geometric product t * ps of trivector t multiplied by a quadvector ps from the right
// trivector * quadvector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                     PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), t.w) * ctype(ps);
}

// geometric product ps * B of a quadvector ps (=3dp pseudoscalar) multiplied from the
// left to the bivector B quadvector * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                       BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) *
           BiVec3dp<ctype>(B.mx, B.my, B.mz, ctype(0.0), ctype(0.0), ctype(0.0));
}

// geometric product B * ps of bivector B multiplied by a quadvector ps from the right
// bivector * quadvector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                       PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(B.mx, B.my, B.mz, ctype(0.0), ctype(0.0), ctype(0.0)) *
           ctype(ps);
}

// geometric product ps * v of a quadvector ps (=3dp pseudoscalar) multiplied from the
// left to the vector v
// quadvector * vector => trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * TriVec3dp<ctype>(-v.x, -v.y, -v.z, ctype(0.0));
}

// geometric product v * ps of a vector v multiplied with a quadvector ps from the right
// vector * quadvector => trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(v.x, v.y, v.z, ctype(0.0)) * ctype(ps);
}

// geometric product ps * s of a quadvector ps and a scalar s
// quadvector * scalar => quadvector (=3dp pseudoscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator*(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(ps) *
                             ctype(s)); // scalar multiplication with a trivector
}

// geometric product s * ps of a scalar s and a quadvector ps
// scalar * quadvector => quadvector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3dp<ctype>(ctype(s) *
                             ctype(ps)); // scalar multiplication with a quadvector
}

// t1 * t2 = -dot(t1, t2) + cmt(t1, t2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t1,
                                                        TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-dot(t1, t2), cmt(t1, t2));
}

// t * B = -(t >> B) + cmt(t, B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-(t >> B), cmt(t, B));
}

// B * t = -(B << t) + cmt(B, t)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-(B << t), cmt(B, t));
}

// t * v = (t >> v) + wdg(t, v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>((t >> v), wdg(t, v));
}

// v * t = (v << t) + wdg(v,t)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>((v << t), wdg(v, t));
}

// (geometric) product t * s of trivector t multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(TriVec3dp<T> const& t,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return t * ctype(s);
}

// (geometric) product s * t of a scalar s multiplied from the left to the trivector t
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * t;
}

// geometric product A * B between two bivectors:
//
// the full geometric bivector product only exists in >= 4d spaces:
// A * B = gr0(A * B) + gr2(A * B) + gr4(A * B) = -dot(A,B) + cmt(A,B) + wdg(A,B)
// In 3D we don't have gr4(A * B) and thus only the terms up to grade 3 remain.
// The bivector product AxB = cmt(A,B) = 0.5*(A*B-B*A) is called the commutator product.
//
// A * B = -dot(A,B) + cmt(A,B) + wdg(A,B)  (in 4d and higher dimensional spaces)
//
// => bivector*bivector = scalar + bivector + quadvector = even grade multivector (in 4d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B1,
                                                        BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-dot(B1, B2), cmt(B1, B2), wdg(B1, B2));
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                        Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>((v << B), wdg(B, v));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>((B >> v), wdg(v, B));
}

// (geometric) product B * s of bivector B multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(BiVec3dp<T> const& B,
                                                       Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// (geometric) product s * B of a scalar s multiplied from the left to the bivector B
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s,
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
constexpr MVec3dp_E<std::common_type_t<T, U>> operator*(Vec3dp<T> const& a,
                                                        Vec3dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(dot(a, b), wdg(a, b));
}

// (geometric) product v * s of vector v multiplied with scalar s from the right
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(Vec3dp<T> const& v, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// (geometric) product s * v of scalar s multiplied from the left to vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s, Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar3dp<std::common_type_t<T, U>> operator*(Scalar3dp<T> s1, Scalar3dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> rgpr(MVec3dp<T> const& A, MVec3dp<U> const& B)
{
    // regressive geometric product of two fully populated 3dp multivectors
    // => due to the degenerate algebra some terms are not present in G<3,0,1> compared to
    // G<4,0,0>
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c15 + A.c1 * B.c11 + A.c2 * B.c12 + A.c3 * B.c13 + A.c4 * B.c14 -
               A.c5 * B.c8 - A.c6 * B.c9 - A.c7 * B.c10 - A.c8 * B.c5 - A.c9 * B.c6 -
               A.c10 * B.c7 - A.c11 * B.c1 - A.c12 * B.c2 - A.c13 * B.c3 - A.c14 * B.c4 +
               A.c15 * B.c0;
    ctype c1 = A.c0 * B.c11 + A.c1 * B.c15 + A.c2 * B.c7 - A.c3 * B.c6 - A.c4 * B.c8 +
               A.c5 * B.c14 + A.c6 * B.c3 - A.c7 * B.c2 + A.c8 * B.c4 + A.c9 * B.c13 -
               A.c10 * B.c12 - A.c11 * B.c0 - A.c12 * B.c10 + A.c13 * B.c9 +
               A.c14 * B.c5 + A.c15 * B.c1;
    ctype c2 = A.c0 * B.c12 - A.c1 * B.c7 + A.c2 * B.c15 + A.c3 * B.c5 - A.c4 * B.c9 -
               A.c5 * B.c3 + A.c6 * B.c14 + A.c7 * B.c1 - A.c8 * B.c13 + A.c9 * B.c4 +
               A.c10 * B.c11 + A.c11 * B.c10 - A.c12 * B.c0 - A.c13 * B.c8 +
               A.c14 * B.c6 + A.c15 * B.c2;
    ctype c3 = A.c0 * B.c13 + A.c1 * B.c6 - A.c2 * B.c5 + A.c3 * B.c15 - A.c4 * B.c10 +
               A.c5 * B.c2 - A.c6 * B.c1 + A.c7 * B.c14 + A.c8 * B.c12 - A.c9 * B.c11 +
               A.c10 * B.c4 - A.c11 * B.c9 + A.c12 * B.c8 - A.c13 * B.c0 + A.c14 * B.c7 +
               A.c15 * B.c3;
    ctype c4 = A.c4 * B.c15 - A.c5 * B.c11 - A.c6 * B.c12 - A.c7 * B.c13 - A.c11 * B.c5 -
               A.c12 * B.c6 - A.c13 * B.c7 + A.c15 * B.c4;
    ctype c5 = -A.c4 * B.c11 + A.c5 * B.c15 + A.c6 * B.c7 - A.c7 * B.c6 - A.c11 * B.c4 -
               A.c12 * B.c13 + A.c13 * B.c12 + A.c15 * B.c5;
    ctype c6 = -A.c4 * B.c12 - A.c5 * B.c7 + A.c6 * B.c15 + A.c7 * B.c5 + A.c11 * B.c13 -
               A.c12 * B.c4 - A.c13 * B.c11 + A.c15 * B.c6;
    ctype c7 = -A.c4 * B.c13 + A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c15 - A.c11 * B.c12 +
               A.c12 * B.c11 - A.c13 * B.c4 + A.c15 * B.c7;
    ctype c8 = A.c0 * B.c5 - A.c1 * B.c4 - A.c2 * B.c13 + A.c3 * B.c12 + A.c4 * B.c1 +
               A.c5 * B.c0 + A.c6 * B.c10 - A.c7 * B.c9 + A.c8 * B.c15 + A.c9 * B.c7 -
               A.c10 * B.c6 + A.c11 * B.c14 + A.c12 * B.c3 - A.c13 * B.c2 -
               A.c14 * B.c11 + A.c15 * B.c8;
    ctype c9 = A.c0 * B.c6 + A.c1 * B.c13 - A.c2 * B.c4 - A.c3 * B.c11 + A.c4 * B.c2 -
               A.c5 * B.c10 + A.c6 * B.c0 + A.c7 * B.c8 - A.c8 * B.c7 + A.c9 * B.c15 +
               A.c10 * B.c5 - A.c11 * B.c3 + A.c12 * B.c14 + A.c13 * B.c1 -
               A.c14 * B.c12 + A.c15 * B.c9;
    ctype c10 = A.c0 * B.c7 - A.c1 * B.c12 + A.c2 * B.c11 - A.c3 * B.c4 + A.c4 * B.c3 +
                A.c5 * B.c9 - A.c6 * B.c8 + A.c7 * B.c0 + A.c8 * B.c6 - A.c9 * B.c5 +
                A.c10 * B.c15 + A.c11 * B.c2 - A.c12 * B.c1 + A.c13 * B.c14 -
                A.c14 * B.c13 + A.c15 * B.c10;
    ctype c11 = A.c4 * B.c5 + A.c5 * B.c4 + A.c6 * B.c13 - A.c7 * B.c12 + A.c11 * B.c15 +
                A.c12 * B.c7 - A.c13 * B.c6 + A.c15 * B.c11;
    ctype c12 = A.c4 * B.c6 - A.c5 * B.c13 + A.c6 * B.c4 + A.c7 * B.c11 - A.c11 * B.c7 +
                A.c12 * B.c15 + A.c13 * B.c5 + A.c15 * B.c12;
    ctype c13 = A.c4 * B.c7 + A.c5 * B.c12 - A.c6 * B.c11 + A.c7 * B.c4 + A.c11 * B.c6 -
                A.c12 * B.c5 + A.c13 * B.c15 + A.c15 * B.c13;
    ctype c14 = A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 - A.c3 * B.c7 - A.c4 * B.c0 -
                A.c5 * B.c1 - A.c6 * B.c2 - A.c7 * B.c3 - A.c8 * B.c11 - A.c9 * B.c12 -
                A.c10 * B.c13 + A.c11 * B.c8 + A.c12 * B.c9 + A.c13 * B.c10 +
                A.c14 * B.c15 + A.c15 * B.c14;
    ctype c15 = -A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7 + A.c11 * B.c11 +
                A.c12 * B.c12 + A.c13 * B.c13 + A.c15 * B.c15;
    return MVec3dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,
                          c15);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& A,
                                                   MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
        Scalar3dp<ctype>(A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 -
                         A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0),
        BiVec3dp<ctype>(A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c7 * B.c1,
                        -A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c7 * B.c2,
                        A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 + A.c7 * B.c3,
                        A.c0 * B.c1 + A.c1 * B.c0 + A.c2 * B.c6 - A.c3 * B.c5 +
                            A.c4 * B.c7 + A.c5 * B.c3 - A.c6 * B.c2 + A.c7 * B.c4,
                        A.c0 * B.c2 - A.c1 * B.c6 + A.c2 * B.c0 + A.c3 * B.c4 -
                            A.c4 * B.c3 + A.c5 * B.c7 + A.c6 * B.c1 + A.c7 * B.c5,
                        A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c0 +
                            A.c4 * B.c2 - A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6),
        PScalar3dp<ctype>(-A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c7 * B.c7));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& A,
                                                   MVec3dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 - A.c3 * B.c4 -
                          A.c4 * B.c0 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1,
                      -A.c0 * B.c3 + A.c1 * B.c7 + A.c2 * B.c1 - A.c3 * B.c5 +
                          A.c4 * B.c6 - A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2,
                      A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c7 - A.c3 * B.c6 -
                          A.c4 * B.c5 + A.c5 * B.c4 - A.c6 * B.c0 + A.c7 * B.c3,
                      A.c3 * B.c7 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3),
        TriVec3dp<ctype>(A.c3 * B.c1 + A.c4 * B.c7 + A.c5 * B.c3 - A.c6 * B.c2,
                         A.c3 * B.c2 - A.c4 * B.c3 + A.c5 * B.c7 + A.c6 * B.c1,
                         A.c3 * B.c3 + A.c4 * B.c2 - A.c5 * B.c1 + A.c6 * B.c7,
                         -A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 - A.c3 * B.c0 +
                             A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& M,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(M.c0 * t.x + M.c1 * t.w + M.c5 * t.z - M.c6 * t.y,
                      M.c0 * t.y + M.c2 * t.w - M.c4 * t.z + M.c6 * t.x,
                      M.c0 * t.z + M.c3 * t.w + M.c4 * t.y - M.c5 * t.x,
                      -M.c1 * t.x - M.c2 * t.y - M.c3 * t.z),
        TriVec3dp<ctype>(M.c2 * t.z - M.c3 * t.y + M.c7 * t.x,
                         -M.c1 * t.z + M.c3 * t.x + M.c7 * t.y,
                         M.c1 * t.y - M.c2 * t.x + M.c7 * t.z,
                         -M.c4 * t.x - M.c5 * t.y - M.c6 * t.z + M.c7 * t.w));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& M,
                                                   TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
        Scalar3dp<ctype>(M.c0 * t.x + M.c1 * t.y + M.c2 * t.z + M.c3 * t.w),
        BiVec3dp<ctype>(-M.c3 * t.x - M.c5 * t.z + M.c6 * t.y,
                        -M.c3 * t.y + M.c4 * t.z - M.c6 * t.x,
                        -M.c3 * t.z - M.c4 * t.y + M.c5 * t.x,
                        -M.c1 * t.z + M.c2 * t.y + M.c4 * t.w - M.c7 * t.x,
                        M.c0 * t.z - M.c2 * t.x + M.c5 * t.w - M.c7 * t.y,
                        -M.c0 * t.y + M.c1 * t.x + M.c6 * t.w - M.c7 * t.z),
        PScalar3dp<ctype>(M.c4 * t.x + M.c5 * t.y + M.c6 * t.z));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& M,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(Scalar3dp<ctype>(-M.c1 * B.mx - M.c2 * B.my - M.c3 * B.mz -
                                             M.c4 * B.vx - M.c5 * B.vy - M.c6 * B.vz),
                            BiVec3dp<ctype>(M.c2 * B.vz - M.c3 * B.vy + M.c7 * B.vx,
                                            -M.c1 * B.vz + M.c3 * B.vx + M.c7 * B.vy,
                                            M.c1 * B.vy - M.c2 * B.vx + M.c7 * B.vz,
                                            M.c0 * B.vx + M.c2 * B.mz - M.c3 * B.my +
                                                M.c5 * B.vz - M.c6 * B.vy + M.c7 * B.mx,
                                            M.c0 * B.vy - M.c1 * B.mz + M.c3 * B.mx -
                                                M.c4 * B.vz + M.c6 * B.vx + M.c7 * B.my,
                                            M.c0 * B.vz + M.c1 * B.my - M.c2 * B.mx +
                                                M.c4 * B.vy - M.c5 * B.vx + M.c7 * B.mz),
                            PScalar3dp<ctype>(-M.c1 * B.vx - M.c2 * B.vy - M.c3 * B.vz));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_U<T> const& M,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(Vec3dp<ctype>(M.c1 * B.vz - M.c2 * B.vy - M.c3 * B.mx -
                                              M.c5 * B.mz + M.c6 * B.my + M.c7 * B.vx,
                                          -M.c0 * B.vz + M.c2 * B.vx - M.c3 * B.my +
                                              M.c4 * B.mz - M.c6 * B.mx + M.c7 * B.vy,
                                          M.c0 * B.vy - M.c1 * B.vx - M.c3 * B.mz -
                                              M.c4 * B.my + M.c5 * B.mx + M.c7 * B.vz,
                                          -M.c4 * B.vx - M.c5 * B.vy - M.c6 * B.vz),
                            TriVec3dp<ctype>(M.c3 * B.vx + M.c5 * B.vz - M.c6 * B.vy,
                                             M.c3 * B.vy - M.c4 * B.vz + M.c6 * B.vx,
                                             M.c3 * B.vz + M.c4 * B.vy - M.c5 * B.vx,
                                             -M.c0 * B.vx - M.c1 * B.vy - M.c2 * B.vz +
                                                 M.c4 * B.mx + M.c5 * B.my +
                                                 M.c6 * B.mz));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(MVec3dp_E<T> const& M,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(M.c2 * v.z - M.c3 * v.y + M.c4 * v.w + M.c7 * v.x,
                      -M.c1 * v.z + M.c3 * v.x + M.c5 * v.w + M.c7 * v.y,
                      M.c1 * v.y - M.c2 * v.x + M.c6 * v.w + M.c7 * v.z, M.c7 * v.w),
        TriVec3dp<ctype>(M.c1 * v.w, M.c2 * v.w, M.c3 * v.w,
                         M.c0 * v.w - M.c1 * v.x - M.c2 * v.y - M.c3 * v.z));
}

// rgpr(t1,t2 = rcmt(t1, t2) + rdot(t1, t2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t1,
                                                   TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(rcmt(t1, t2), rdot(t1, t2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_U<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                   BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(
        Vec3dp<ctype>(
            -t.y * B.mz + t.z * B.my + t.w * B.vx, t.x * B.mz - t.z * B.mx + t.w * B.vy,
            -t.x * B.my + t.y * B.mx + t.w * B.vz, -t.x * B.vx - t.y * B.vy - t.z * B.vz),
        TriVec3dp<ctype>(t.y * B.vz - t.z * B.vy, -t.x * B.vz + t.z * B.vx,
                         t.x * B.vy - t.y * B.vx, t.x * B.mx + t.y * B.my + t.z * B.mz));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(TriVec3dp<T> const& t,
                                                   Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(
        Scalar3dp<ctype>(-t.x * v.x - t.y * v.y - t.z * v.z - t.w * v.w),
        BiVec3dp<ctype>(-t.x * v.w, -t.y * v.w, -t.z * v.w, t.y * v.z - t.z * v.y,
                        -t.x * v.z + t.z * v.x, t.x * v.y - t.y * v.x));
}

// rgpr(B1,B2) = rwdg(B1, B2) + rcmt(B1, B2) - rdot(B1, B2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(BiVec3dp<T> const& B1,
                                                   BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(rwdg(B1, B2), rcmt(B1, B2), -rdot(B1, B2));
}

// rgpr(v1,v2) = rcmt(v1, v2) - rdot(v1, v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> rgpr(Vec3dp<T> const& v1,
                                                   Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(rcmt(v1, v2), -rdot(v1, v2));
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
    requires(std::floating_point<T>)
inline Scalar3dp<T> inv(Scalar3dp<T> s)
{
    T sq_n = T(bulk_nrm_sq(s));
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar3dp<T>(rev(s) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline Vec3dp<T> inv(Vec3dp<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v) = v/bulk_sq_nrm(v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    T sq_n = T(bulk_nrm_sq(v));
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return Vec3dp<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
template <typename T>
    requires(std::floating_point<T>)
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

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
template <typename T>
    requires(std::floating_point<T>)
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

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
template <typename T>
    requires(std::floating_point<T>)
inline MVec3dp_E<T> inv(MVec3dp_E<T> const& E)
{
    auto tc = E * conj(E);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("E={}", E);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "even grade multivector");
    return conj(E) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
template <typename T>
    requires(std::floating_point<T>)
inline MVec3dp_U<T> inv(MVec3dp_U<T> const& U)
{
    auto tc = U * conj(U);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("E={}", E);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "odd grade multivector");
    return conj(U) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
// left and a right inverse are the same (see paper of Hitzer, Sangwine)
template <typename T>
    requires(std::floating_point<T>)
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
// left_bulk_contract3dp(a,b) = rwdg(left_bulk_dual(a), b)
// left_weight_contract3dp(a,b) = rwdg(left_weight_dual(a), b)
//
// right_bulk_contract3dp(a,b) = rwdg(a, right_bulk_dual(b))
// right_weight_contract3dp(a,b) = rwdg(a, right_weight_dual(b))
//
// The contraction subracts the grades of the operands.
//
// When the metric is the identity, these two contractions are identical and are
// called the contraction (they are identical with the right contraction).
// When the metric is degenerate they produce different results.
//
// In general a contraction throws away parts that are perpendicular to each other.
// The result of right_bulk_contract(B,v) lies in B and is perpendicular to v.
/////////////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2>
decltype(auto) left_bulk_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(left_bulk_dual(a), b);
    return rwdg(left_bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) left_weight_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(left_weight_dual(a), b);
    return rwdg(left_weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) right_bulk_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, right_bulk_dual(b));
    return rwdg(std::forward<arg1>(a), right_bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) right_weight_contract3dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, right_weight_dual(b));
    return rwdg(std::forward<arg1>(a), right_weight_dual(std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// Projective expansions for 3dp:
//
// left_bulk_expand3dp(a,b) = wdg(left_bulk_dual(a), b)       (dual to
// left_weight_contract) left_weight_expand3dp(a,b) = wdg(left_weight_dual(a), b)   (dual
// to left_bulk_contract)
//
// right_bulk_expand3dp(a,b) = wdg(a, right_bulk_dual(b))       (dual to
// right_weight_contract) right_weight_expand3dp(a,b) = wdg(a, right_weight_dual(b)) (dual
// to right_bulk_contract)
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

template <typename arg1, typename arg2>
decltype(auto) left_bulk_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(left_bulk_dual(a), b);
    return wdg(left_bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) left_weight_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(left_weight_dual(a), b);
    return wdg(left_weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) right_bulk_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(a, right_bulk_dual(b));
    return wdg(std::forward<arg1>(a), right_bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) right_weight_expand3dp(arg1&& a, arg2&& b)
{
    // return wdg(a, right_weight_dual(b));
    return wdg(std::forward<arg1>(a), right_weight_dual(std::forward<arg2>(b)));
}

} // namespace hd::ga::pga