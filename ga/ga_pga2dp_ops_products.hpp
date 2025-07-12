#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_pga2dp_ops_basics.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga2dp product operations:
//
// - dot()                   -> dot product
// - rdot()                  -> regressive dot product
// - wdg(), join()           -> wedge product (join as convenience interface)
// - rwdg(), meet()          -> regressive wedge product (meet as convenience interface)
// - operator<<()            -> left contraction (= left bulk contraction)
// - operator>>()            -> right contraction (= right bulk contraction)
// - cmt()                   -> commutator product (= asymmetric part of gpr)
// - rcmt()                  -> regressive commutator product (= asymmetric part of rgpr)
// - operator*()             -> geometric product (=gpr())
// - rgpr()                  -> regressive geometric product
// - inv()                   -> inversion operation (w.r.t. geometric product)
//
// - left_bulk_contract2dp()    -> left bulk contraction
// - left_weight_contract2dp()  -> left weight contraction
// - right_bulk_contract2dp()   -> right bulk contraction
// - right_weight_contract2dp() -> right weight contraction
//
// - left_bulk_expand2dp()      -> left bulk expansion
// - left_weight_expand2dp()    -> left weight expansion
// - right_bulk_expand2dp()     -> right bulk expansion
// - right_weight_expand2dp()   -> right weight expansion
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot products for 2dp (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = 0
//
////////////////////////////////////////////////////////////////////////////////

// returns 0.0 due to degenerate metric with e3^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// return dot product of two bivectors A and B (= a scalar)
// dot(A,B) = gr0(A * B)
//
// arguments: BiVec2dp<T> or Line2d<T> : public BiVec2dp<T>
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(BiVec2dp<T> const& B1,
                                                  BiVec2dp<U> const& B2)
{
    // definition: dot(A, B) = gr0(A*B)
    // -> only the symmetric (i.e. scalar) part remains
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = 0
    // and dot(e23, e23) = 0, dot(e31,e31) = 0, dot(e12,e12) = -1
    // and all other dot(exy,ezw)=0
    // (every index combination containing an index twice becomes zero)
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(B1.z * B2.z);
}

// return dot-product of two vectors in the modelled space G<2,0,1>
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
// dot(v,v) = |v|^2 = gr0(v*rev(v))
//
// arguments: Vec2dp<T> or Point2dp<T> : public Vec2dp<T>
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(Vec2dp<T> const& v1,
                                                  Vec2dp<U> const& v2)
{
    // definition: dot(v1, v2) = (v1)^T g_12 v2 with the metric g_12
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = 0
    // as diagonal elements of g_12, thus the z-component does not matter
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(v1.x * v2.x + v1.y * v2.y);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp<T> const& A,
                                                  MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // return Scalar2dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c6 * B.c6);
    return Scalar2dp<ctype>(dot(gr0(A), gr0(B)) + dot(gr1(A), gr1(B)) +
                            dot(gr2(A), gr2(B)) + dot(gr3(A), gr3(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp_E<T> const& A,
                                                  MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(dot(gr0(A), gr0(B)) + dot(gr2(A), gr2(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp_U<T> const& A,
                                                  MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(dot(gr1(A), gr1(B)) + dot(gr3(A), gr3(B)));
}


////////////////////////////////////////////////////////////////////////////////
// regressive dot products for 2dp (= defined for equal grades exclusively)
//
// rdot(v1,v2) = cmpl( dot(cmpl(v1),cmpl(v2)) )
//
// returns a pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(PScalar2dp<T> ps1, PScalar2dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps1) * ctype(ps2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(BiVec2dp<T> const& B1,
                                                    BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(B1.x * B2.x + B1.y * B2.y);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(Vec2dp<T> const& v1,
                                                    Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(v1.z * v2.z);
}

// returns 0.0 due to degenerate metric with e3^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot([[maybe_unused]] Scalar2dp<T>,
                                                    [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(MVec2dp<T> const& A,
                                                    MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(rdot(gr0(A), gr0(B)) + rdot(gr1(A), gr1(B)) +
                             rdot(gr2(A), gr2(B)) + rdot(gr3(A), gr3(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(MVec2dp_E<T> const& A,
                                                    MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(rdot(gr0(A), gr0(B)) + rdot(gr2(A), gr2(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> rdot(MVec2dp_U<T> const& A,
                                                    MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(rdot(gr1(A), gr1(B)) + rdot(gr3(A), gr3(B)));
}

////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////

// wedge product extended to fully populated multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype c4 = A.c0 * B.c4 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c3 * B.c1 + A.c5 * B.c0;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c0;
    ctype c7 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
               A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// wedge product between two pseudoscalars
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a pseudoscalar and a bivector
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a bivector and a pseudoscalar
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec2dp<T> const&,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a pseudoscalar and a vector
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                  [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a vector and a pseudoscalar
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec2dp<T> const&,
                                                  [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a trivector ps and a scalar s
// wdg(ps,s) = ps*s*e123
// => returns a trivector (scalar s multiple of the trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// wedge product between a scalar s and a trivector ps
// wdg(s,ps) = s*ps*e123
// => returns a trivector (scalar s multiple of the trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// wedge product between two bivectors
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] BiVec2dp<T> const&,
                                                  [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between bivector B and vector v
// wdg(B,v) = gr3(B * v)
// => returns a trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-B.x * v.x - B.y * v.y - B.z * v.z);
}

// wedge product between a vector a and a bivector B
// => returns a trivector (i.e. the pseudoscalar for 2dp)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-v.x * B.x - v.y * B.y - v.z * B.z);
}

// wedge product with one scalar (returns a scaled bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// wedge product with one scalar (returns a scaled bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, BiVec2dp<U> const& B)
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
constexpr BiVec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    return BiVec2dp<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// wedge product between two scalars (returns a scalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
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
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Line2d<std::common_type_t<T, U>> wdg(Point2d<T> const& p, Point2d<U> const& q)
{
    return Line2d<std::common_type_t<T, U>>(p.y - q.y, q.x - p.x, p.x * q.y - p.y * q.x);
}


////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg -> join
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> join(BiVec2dp<T> const& B,
                                                    Vec2dp<U> const& v)
{
    return wdg(B, v);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& v,
                                                    BiVec2dp<U> const& B)
{
    return wdg(v, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& v1,
                                                  Vec2dp<U> const& v2)
{
    return wdg(v1, v2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
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

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
               A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    ctype c1 = A.c1 * B.c7 - A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype c2 = A.c2 * B.c7 + A.c4 * B.c6 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype c3 = A.c3 * B.c7 - A.c4 * B.c5 + A.c5 * B.c4 + A.c7 * B.c3;
    ctype c4 = A.c4 * B.c7 + A.c7 * B.c4;
    ctype c5 = A.c5 * B.c7 + A.c7 * B.c5;
    ctype c6 = A.c6 * B.c7 + A.c7 * B.c6;
    ctype c7 = A.c7 * B.c7;
    return MVec2dp<std::common_type_t<T, U>>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * B;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(ps);
}

// regressive wedge product of two bivectors
// rwdg(a,b) = !wgd(!a,!b) with ! as the complement operation
// => returns a vector
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{
    return Vec2dp<std::common_type_t<T, U>>(-B1.y * B2.z + B1.z * B2.y,
                                            B1.x * B2.z - B1.z * B2.x,
                                            -B1.x * B2.y + B1.y * B2.x);
}

// regressive wedge product between a vector and a bivector
// rwdg(A,b) = !wgd(!A,!b) with ! as the complement operation
// => returns a scalar
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(-B.x * v.x - B.y * v.y - B.z * v.z);
}

// regressive wedge product between a vector and a bivector
// rwdg(a,B) = !wgd(!a,!B) with ! as the complement operation
// => returns a scalar
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& v,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(-v.x * B.x - v.y * B.y - v.z * B.z);
}

// required to be present for dist2dp (to complile, even if not used)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2dp<T> const&,
                                                   [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

////////////////////////////////////////////////////////////////////////////////
// convenience functions rwdg -> meet
////////////////////////////////////////////////////////////////////////////////

// intersection of two lines -> return a point
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> meet(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{
    return rwdg(B1, B2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> meet(Line2d<T> const& l1, Line2d<U> const& l2)
{
    return rwdg(static_cast<BiVec2dp<T> const&>(l1), static_cast<BiVec2dp<T> const&>(l2));
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// Implements the left bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator<<(a,b) = left_bulk_contract(a,b) = rwdg( bulk_dual(a), b )
//
////////////////////////////////////////////////////////////////////////////////

// left contraction (A << B) of two fully populated multivectors
// returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator<<(MVec2dp<T> const& A,
                                                       MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c6 * B.c6;
    ctype c1 = A.c0 * B.c1 + A.c2 * B.c6;
    ctype c2 = A.c0 * B.c2 - A.c1 * B.c6;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 - A.c6 * B.c7;
    ctype c4 = A.c0 * B.c4 - A.c1 * B.c7;
    ctype c5 = A.c0 * B.c5 - A.c2 * B.c7;
    ctype c6 = A.c0 * B.c6;
    ctype c7 = A.c0 * B.c7;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// left contraction (ps1 << ps2) - pseudoscalar ps1 contracted onto pseudoscalar ps2
// returns a scalar (the product of the pseudoscalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2dp<T>,
                                                         [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction - pseudoscalar contracted onto bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction of a bivector B onto a pseudoscalar ps (=trivector)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator<<(BiVec2dp<T> const& B,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(0.0, 0.0, -B.z * ctype(ps));
}

// left contraction - pseudoscalar contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction of a vector v onto a pseudoscalar ps (=trivector)
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v,
                                                        PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(-v.x, -v.y, 0.0) * ctype(ps);
}

// left contraction - pseudoscalar contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2dp<T>,
                                                         [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction (s << ps) - scalar s contracted onto pseudoscalar ps
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s,
                                                          PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// left contraction (B1 << B2) - bivector B1 contracted onto B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(BiVec2dp<T> const& B1,
                                                         BiVec2dp<U> const& B2)
{
    return dot(B1, B2);
}

// left contraction - bivector contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// return the left contraction of a vector v onto a bivector B
// => returns a vector
// identical to cmt(B,v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v.y * B.z, -v.x * B.z, v.x * B.y - v.y * B.x);
}

// left contraction - bivector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction (s << B) - scalar s contracted onto bivector B
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// left contraction (v1 << v2) - vector v1 contracted onto v2
// (=identical to scalar product dot(v1,v2) of two vectors; identical to (v1 >> v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v1,
                                                         Vec2dp<U> const& v2)
{
    return dot(v1, v2);
}

// left contraction - vector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction (s << v) - scalar s contracted onto vector v
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// left contraction (s1 << s2) - scalar s1 contracted onto scalar s2
// returns a scalar (the product of the scalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
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
// operator>>(a,b) = right_bulk_contract(a,b) = rwdg( a, bulk_dual(b) )
//
////////////////////////////////////////////////////////////////////////////////

// right contraction (A >> B) of two fully populated multivectors
// returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator>>(MVec2dp<T> const& A,
                                                       MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c6 * B.c6;
    ctype c1 = A.c1 * B.c0 - A.c6 * B.c2;
    ctype c2 = A.c2 * B.c0 + A.c6 * B.c1;
    ctype c3 = A.c3 * B.c0 + A.c4 * B.c2 - A.c5 * B.c1 - A.c7 * B.c6;
    ctype c4 = A.c4 * B.c0 - A.c7 * B.c1;
    ctype c5 = A.c5 * B.c0 - A.c7 * B.c2;
    ctype c6 = A.c6 * B.c0;
    ctype c7 = A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// right contraction - pseudoscalar contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] PScalar2dp<T>,
                                                         [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction of pseudoscalar ps (=trivector) by bivector B "ps contracted by B"
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(ctype(0.0), ctype(0.0), -ctype(ps) * B.z);
}

// right contraction - bivector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction of pseudoscalar ps (=trivector) by vector v (" ps contracted by v")
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                        Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * BiVec2dp<ctype>(v.x, v.y, ctype(0.0));
}

// right contraction -  vector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (ps >> s) - pseudoscalar ps contracted by scalar s
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                          Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// right contraction - scalar contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar2dp<T>,
                                                         [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (B1 >> B2) - bivector B1 contracted by bivector B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B1,
                                                         BiVec2dp<U> const& B2)
{
    return dot(B1, B2);
}

// right contraction of bivector B by vector v ("B contracted by v")
// => returns a vector
// (identically to cmt(v,B))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B,
                                                      Vec2dp<U> const& v)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(-B.z * v.y, B.z * v.x, B.x * v.y - B.y * v.x);
}

// right contraction - vector contracted by a bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (B >> s) - bivector B contracted by scalar s
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B,
                                                        Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// right contraction - scalar contracted by a bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dp<T>, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (v1 >> v2) - vector v1 contracted by vector v2
// (=identical to scalar product dot(v1,v2) of two vectors; identical to (v1 << v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(Vec2dp<T> const& v1,
                                                         Vec2dp<U> const& v2)
{
    return dot(v1, v2);
}

// right contraction (v >> s) - vector v contracted by scalar s
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator>>(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// right contraction - scalar contracted by a vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dp<T>, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (s1 >> s2) - scalar s1 contracted by scalar s2
// returns a scalar (the product of the scalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> cmt(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(0.0, -A.c2 * B.c6 + A.c6 * B.c2, A.c1 * B.c6 - A.c6 * B.c1,
                          -A.c1 * B.c5 + A.c2 * B.c4 - A.c4 * B.c2 + A.c5 * B.c1,
                          A.c2 * B.c3 - A.c3 * B.c2 - A.c5 * B.c6 + A.c6 * B.c5,
                          -A.c1 * B.c3 + A.c3 * B.c1 + A.c4 * B.c6 - A.c6 * B.c4,
                          A.c1 * B.c2 - A.c2 * B.c1, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> cmt(BiVec2dp<T> const& B1,
                                                 BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(-B1.y * B2.z + B1.z * B2.y, B1.x * B2.z - B1.z * B2.x, 0.0);
}

// cmt(B,v) = -cmt(v,B)
// identical to (v << B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> cmt(BiVec2dp<T> const& B, Vec2dp<U> const& v)
{
    return (v << B);
}

// cmt(v,B) = -cmt(B,v)
// identcial to (B >> v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> cmt(Vec2dp<T> const& v, BiVec2dp<U> const& B)
{
    return (B >> v);
}

// cmt(v1,v2) = wdg(v1,v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> cmt(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    return wdg(v1, v2);
}


////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rcmt(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(0.0, A.c2 * B.c3 - A.c3 * B.c2 - A.c5 * B.c6 + A.c6 * B.c5,
                          -A.c1 * B.c3 + A.c3 * B.c1 + A.c4 * B.c6 - A.c6 * B.c4,
                          -A.c4 * B.c5 + A.c5 * B.c4, -A.c3 * B.c5 + A.c5 * B.c3,
                          A.c3 * B.c4 - A.c4 * B.c3,
                          A.c1 * B.c5 - A.c2 * B.c4 + A.c4 * B.c2 - A.c5 * B.c1, 0.0);
}

// rcmt(B1,B2) = rwdg(B1,B2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rcmt(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{
    return rwdg(B1, B2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rcmt(BiVec2dp<T> const& B,
                                                  Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(B.y * v.z, -B.x * v.z, B.x * v.y - B.y * v.x);
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> rcmt(Vec2dp<T> const& v,
                                                  BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(-v.z * B.y, v.z * B.x, v.x * B.y - v.y * B.x);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> rcmt(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v1.y * v2.z - v1.z * v2.y, -v1.x * v2.z + v1.z * v2.x, 0.0);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 2dp multivector
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
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                      MVec2dp<U> const& B)
{
    // geometric product of two fully populated 2dp multivectors
    // => due to the degenerate algebra 16 terms are not present in G<2,0,1> compared to
    // G<3,0,0>
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c6 * B.c6;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c6 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c6 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
               A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype c4 = A.c0 * B.c4 - A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0 -
               A.c5 * B.c6 + A.c6 * B.c5 - A.c7 * B.c1;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c3 - A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 +
               A.c5 * B.c0 - A.c6 * B.c4 - A.c7 * B.c2;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c6 * B.c0;
    ctype c7 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
               A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of a multivector with an even grade multivector
// multivector * even grade multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                      MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c6 * B.c3;
    ctype c1 = A.c1 * B.c0 - A.c2 * B.c3;
    ctype c2 = A.c1 * B.c3 + A.c2 * B.c0;
    ctype c3 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0 + A.c7 * B.c3;
    ctype c4 = A.c0 * B.c1 + A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype c5 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1;
    ctype c6 = A.c0 * B.c3 + A.c6 * B.c0;
    ctype c7 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 + A.c7 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of an even grade multivector with a multivector
// even grade multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                      MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c3 * B.c6;
    ctype c1 = A.c0 * B.c1 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c7;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4;
    ctype c6 = A.c0 * B.c6 + A.c3 * B.c0;
    ctype c7 = A.c0 * B.c7 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// multivector * odd multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                      MVec2dp_U<U> const& B)
{
    // geometric product of two fully populated 2dp multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c1 * B.c0 + A.c2 * B.c1;
    ctype c1 = A.c0 * B.c0 + A.c6 * B.c1;
    ctype c2 = A.c0 * B.c1 - A.c6 * B.c0;
    ctype c3 = A.c0 * B.c2 - A.c4 * B.c1 + A.c5 * B.c0 + A.c6 * B.c3;
    ctype c4 = -A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1 - A.c7 * B.c0;
    ctype c5 = -A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0 - A.c7 * B.c1;
    ctype c6 = A.c1 * B.c1 - A.c2 * B.c0;
    ctype c7 = A.c0 * B.c3 - A.c4 * B.c0 - A.c5 * B.c1 - A.c6 * B.c2;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// odd multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                      MVec2dp<U> const& B)
{
    // geometric product of two fully populated 2dp multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c1 + A.c1 * B.c2;
    ctype c1 = A.c0 * B.c0 - A.c1 * B.c6;
    ctype c2 = A.c0 * B.c6 + A.c1 * B.c0;
    ctype c3 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c0 + A.c3 * B.c6;
    ctype c4 = -A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 - A.c3 * B.c1;
    ctype c5 = -A.c0 * B.c3 - A.c1 * B.c7 + A.c2 * B.c1 - A.c3 * B.c2;
    ctype c6 = A.c0 * B.c2 - A.c1 * B.c1;
    ctype c7 = -A.c0 * B.c4 - A.c1 * B.c5 - A.c2 * B.c6 + A.c3 * B.c0;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product A * ps of a multivector A with the trivector ps from the right
// multivector * trivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), A.c6, -A.c1, -A.c2,
                          ctype(0.0), A.c0) *
           ctype(ps);
}

// geometric product ps * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the multivector B
// trivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                      MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec2dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), B.c6, -B.c1,
                                      -B.c2, ctype(0.0), B.c0);
}

// geometric product A * B for two multivectors from the even subalgebra (3d case)
// even grade multivector * even grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(A.c0 * B.c0 - A.c3 * B.c3),
        BiVec2dp<ctype>(A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2,
                        A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1,
                        A.c0 * B.c3 + A.c3 * B.c0));
}

// geometric product A * B for two multivectors from the odd subalgebra (3d case)
// odd grade multivector * odd grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                        MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1),
        BiVec2dp<ctype>(-A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 - A.c3 * B.c0,
                        -A.c0 * B.c2 - A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1,
                        A.c0 * B.c1 - A.c1 * B.c0));
}

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a multivector B of the odd subalgebra
// even grade multivector * odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(A.c0 * B.c0 + A.c3 * B.c1, A.c0 * B.c1 - A.c3 * B.c0,
                      A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 + A.c3 * B.c3),
        PScalar2dp<ctype>(A.c0 * B.c3 - A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2));
}

// geometric product A * B of a multivector A from the odd subalgebra (3d case)
// with a multivector B of the even subalgebra
// odd grade multivector * even grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(A.c0 * B.c0 - A.c1 * B.c3, A.c0 * B.c3 + A.c1 * B.c0,
                      -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c0 + A.c3 * B.c3),
        PScalar2dp<ctype>(-A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0));
}

// geometric product A * ps of an even multivector A multiplied from the right
// by the trivector ps (=2dp pseudoscalar)
// even grade multivector * trivector => odd multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<U> const& A,
                                                        PScalar2dp<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(ctype(0.0), ctype(0.0), A.c3),
                            PScalar2dp<ctype>(A.c0)) *
           ctype(ps);
}

// geometric product ps * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the even grade multivector B
// trivector * even grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec2dp_U<ctype>(Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.c3),
                                        PScalar2dp<ctype>(B.c0));
}

// geometric product A * B of an odd grade multivector A multiplied from the right
// by the trivector B (=2dp pseudoscalar)
// odd grade multivector * trivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<U> const& A,
                                                        PScalar2dp<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(ctype(0.0)),
                            BiVec2dp<ctype>(-A.c0, -A.c1, ctype(0.0))) *
           ctype(ps);
}

// geometric product A * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the odd grade multivector B
// trivector * odd grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                        MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * MVec2dp_E<ctype>(Scalar2dp<ctype>(ctype(0.0)),
                                         BiVec2dp<ctype>(B.c0, B.c1, ctype(0.0)));
}

// geometric product M * B of a multivector M from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& M,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(-M.c3 * B.z),
                            BiVec2dp<ctype>(M.c0 * B.x - M.c2 * B.z + M.c3 * B.y,
                                            M.c0 * B.y + M.c1 * B.z - M.c3 * B.x,
                                            M.c0 * B.z));
}

// geometric product A * M of a bivector A with an even grade multivector M
// from the right
// bivector * even grade multivector => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                        MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(-B.z * M.c3),
                            BiVec2dp<ctype>(B.x * M.c0 - B.y * M.c3 + B.z * M.c2,
                                            B.x * M.c3 + B.y * M.c0 - B.z * M.c1,
                                            B.z * M.c0));
}

// geometric product M * B of an odd grade multivector M with a bivector B
// odd grade multivector * bivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& M,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(-M.c1 * B.z, M.c0 * B.z, -M.c0 * B.y + M.c1 * B.x + M.c3 * B.z),
        PScalar2dp<ctype>(-M.c0 * B.x - M.c1 * B.y - M.c2 * B.z));
}

// geometric product B * M of a bivector A with an odd grade multivector B
// bivector * odd grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                        MVec2dp_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(B.z * M.c1, -B.z * M.c0, -B.x * M.c1 + B.y * M.c0 + B.z * M.c3),
        PScalar2dp<ctype>(-B.x * M.c0 - B.y * M.c1 - B.z * M.c2));
}

// geometric product A * v of an even grade multivector A with a vector v
// from the right
// even grade multivector * vector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(A.c0 * v.x + A.c3 * v.y,
                                          A.c0 * v.y - A.c3 * v.x,
                                          A.c0 * v.z - A.c1 * v.y + A.c2 * v.x),
                            PScalar2dp<ctype>(-A.c1 * v.x - A.c2 * v.y - A.c3 * v.z));
}

// geometric product v * B of a vector v and an even grade multivector B from the left
// vector * even grade multivector => odd grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(v.x * B.c0 - v.y * B.c3,
                                          v.x * B.c3 + v.y * B.c0,
                                          -v.x * B.c2 + v.y * B.c1 + v.z * B.c0),
                            PScalar2dp<ctype>(-v.x * B.c1 - v.y * B.c2 - v.z * B.c3));
}

// geometric product of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator*([[maybe_unused]] PScalar2dp<T>,
                                                        [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0); // trivectors in pga2dp contain e3^2 = 0
}

// geometric product ps * B of a trivector ps (=2dp pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                     BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.z);
}

// geometric product B * ps of bivector A multiplied by a trivector ps from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                     PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.z) * ctype(ps);
}

// geometric product ps * v of a trivector ps (=2dp pseudoscalar) multiplied from the
// left to the vector v trivector * vector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                       Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * BiVec2dp<ctype>(v.x, v.y, ctype(0.0));
}

// geometric product v * ps of a vector a multiplied with a trivector ps from the right
// vector * trivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                       PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(-v.x, -v.y, ctype(0.0)) * ctype(ps);
}

// geometric product ps * s of a trivector ps and a scalar s
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) *
                             ctype(s)); // scalar multiplication with a trivector
}

// geometric product s * ps of a scalar s and a trivector ps
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) *
                             ctype(ps)); // scalar multiplication with a trivector
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
// => bivector*bivector = scalar + bivector = even grade multivector (in 3d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
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
    requires(std::floating_point<T> && std::floating_point<U>)
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
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>((B >> v), wdg(v, B));
}

// (geometric) product B * s of bivector B multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                       Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// (geometric) product s * B of a scalar s multiplied from the left to the bivector B
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s,
                                                       BiVec2dp<U> const& B)
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
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(Vec2dp<T> const& a,
                                                        Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(dot(a, b), wdg(a, b));
}

// (geometric) product v * s of vector v multiplied with scalar s from the right
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// (geometric) product s * v of scalar s multiplied from the left to vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& A, MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6 - A.c4 * B.c1 -
               A.c5 * B.c2 - A.c6 * B.c3 + A.c7 * B.c0;
    ctype c1 = -A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 - A.c4 * B.c0 -
               A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype c2 = -A.c0 * B.c5 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 -
               A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype c3 = A.c3 * B.c7 - A.c4 * B.c5 + A.c5 * B.c4 + A.c7 * B.c3;
    ctype c4 = -A.c3 * B.c5 + A.c4 * B.c7 + A.c5 * B.c3 + A.c7 * B.c4;
    ctype c5 = A.c3 * B.c4 - A.c4 * B.c3 + A.c5 * B.c7 + A.c7 * B.c5;
    ctype c6 = A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 + A.c3 * B.c0 + A.c4 * B.c2 -
               A.c5 * B.c1 + A.c6 * B.c7 + A.c7 * B.c6;
    ctype c7 = -A.c3 * B.c3 + A.c4 * B.c4 + A.c5 * B.c5 + A.c7 * B.c7;
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
                                                   MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(-A.c0 * B.c1 - A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2,
                      -A.c0 * B.c2 + A.c1 * B.c3 - A.c2 * B.c0 - A.c3 * B.c1,
                      -A.c1 * B.c2 + A.c2 * B.c1),
        PScalar2dp<ctype>(A.c1 * B.c1 + A.c2 * B.c2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
                                                   MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0,
                      -A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1,
                      A.c2 * B.c3 + A.c3 * B.c2),
        PScalar2dp<ctype>(-A.c2 * B.c2 + A.c3 * B.c3));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
                                                   MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(A.c0 * B.c3 - A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2),
        BiVec2dp<ctype>(A.c1 * B.c3 + A.c2 * B.c2, -A.c1 * B.c2 + A.c2 * B.c3,
                        A.c0 * B.c2 + A.c1 * B.c1 - A.c2 * B.c0 + A.c3 * B.c3));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
                                                   MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(-A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 + A.c3 * B.c0),
        BiVec2dp<ctype>(-A.c2 * B.c2 + A.c3 * B.c1, A.c2 * B.c1 + A.c3 * B.c2,
                        A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 + A.c3 * B.c3));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& M,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(-M.c0 * B.x - M.c2 * B.z + M.c3 * B.y,
                                          -M.c0 * B.y + M.c1 * B.z - M.c3 * B.x,
                                          -M.c1 * B.y + M.c2 * B.x),
                            PScalar2dp<ctype>(M.c1 * B.x + M.c2 * B.y));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& M,
                                                   BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(-M.c0 * B.x - M.c1 * B.y - M.c2 * B.z),
                            BiVec2dp<ctype>(-M.c2 * B.y + M.c3 * B.x,
                                            M.c2 * B.x + M.c3 * B.y,
                                            M.c0 * B.y - M.c1 * B.x + M.c3 * B.z));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& M,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(-M.c1 * v.x - M.c2 * v.y - M.c3 * v.z),
        BiVec2dp<ctype>(M.c2 * v.z, -M.c1 * v.z, M.c0 * v.z + M.c1 * v.y - M.c2 * v.x));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& M,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(M.c1 * v.z - M.c2 * v.y + M.c3 * v.x,
                                          -M.c0 * v.z + M.c2 * v.x + M.c3 * v.y,
                                          M.c3 * v.z),
                            PScalar2dp<ctype>(-M.c2 * v.z));
}

// rgpr(B1,B2) = rwdg(B1, B2) + rdot(B1, B2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B1,
                                                   BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(rwdg(B1, B2), rdot(B1, B2));
}

// rgpr(B,v) = rwdg(B, v) + rcmt(B, v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B,
                                                   Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(rwdg(B, v), rcmt(B, v));
}

// rgpr(v1,v2) = rcmt(v1,v2) - rdot(v1,v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v1,
                                                   Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(rcmt(v1, v2), -rdot(v1, v2));
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
inline Scalar2dp<T> inv(Scalar2dp<T> s)
{
    T sq_n = bulk_nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar2dp<T>(rev(s) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
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
    requires(std::floating_point<T>)
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
    requires(std::floating_point<T>)
inline MVec2dp_E<T> inv(MVec2dp_E<T> const& E)
{
    T sq_n = bulk_nrm_sq(E);
    hd::ga::detail::check_normalization<T>(sq_n, "even grade multivector");
    T inv = T(1.0) / sq_n;
    return MVec2dp_E<T>(rev(E) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline MVec2dp_U<T> inv(MVec2dp_U<T> const& U)
{
    T sq_n = bulk_nrm_sq(U);
    hd::ga::detail::check_normalization<T>(sq_n, "odd grade multivector");
    T inv = T(1.0) / sq_n;
    return MVec2dp_U<T>(rev(U) * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
// left and a right inverse are the same (see paper of Hitzer, Sangwine)
template <typename T>
    requires(std::floating_point<T>)
inline MVec2dp<T> inv(MVec2dp<T> const& M)
{
    T m_conjm = gr0(M * conj(M) * gr_inv(M) * rev(M));
    hd::ga::detail::check_normalization<T>(std::abs(m_conjm), "multivector");
    T inv = T(1.0) / m_conjm;
    return MVec2dp<T>(conj(M) * gr_inv(M) * rev(M) * inv);
}


/////////////////////////////////////////////////////////////////////////////////////////
// Projective contractions for 2dp:
//
// left_bulk_contract2dp(a,b) = rwdg(bulk_dual(a), b)
// left_weight_contract2dp(a,b) = rwdg(weight_dual(a), b)
//
// right_bulk_contract2dp(a,b) = rwdg(a, bulk_dual(b))
// right_weight_contract2dp(a,b) = rwdg(a, weight_dual(b))
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
decltype(auto) left_bulk_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(bulk_dual(a), b);
    return rwdg(bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) left_weight_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(weight_dual(a), b);
    return rwdg(weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) right_bulk_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, bulk_dual(b));
    return rwdg(std::forward<arg1>(a), bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) right_weight_contract2dp(arg1&& a, arg2&& b)
{
    // return rwdg(a, weight_dual(b));
    return rwdg(std::forward<arg1>(a), weight_dual(std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// Projective expansions for 2dp:
//
// left_bulk_expand2dp(a,b) = wdg(bulk_dual(a), b)       (dual to left_weight_contract)
// left_weight_expand2dp(a,b) = wdg(weight_dual(a), b)   (dual to left_bulk_contract)
//
// right_bulk_expand2dp(a,b) = wdg(a, bulk_dual(b))       (dual to right_weight_contract)
// right_weight_expand2dp(a,b) = wdg(a, weight_dual(b))   (dual to right_bulk_contract)
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
decltype(auto) left_bulk_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(bulk_dual(a), b);
    return wdg(bulk_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) left_weight_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(weight_dual(a), b);
    return wdg(weight_dual(std::forward<arg1>(a)), std::forward<arg2>(b));
}

template <typename arg1, typename arg2>
decltype(auto) right_bulk_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(a, bulk_dual(b));
    return wdg(std::forward<arg1>(a), bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) right_weight_expand2dp(arg1&& a, arg2&& b)
{
    // return wdg(a, weight_dual(b));
    return wdg(std::forward<arg1>(a), weight_dual(std::forward<arg2>(b)));
}

} // namespace hd::ga::pga