#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "detail/ga_foundation.hpp"     // standard library headers and GA infrastructure
#include "detail/ga_mvec2dp.hpp"        // 2dp multivector types
#include "detail/type_t/ga_type2dp.hpp" // Point2dp, Vector2d, Point2d, Line2d

namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> gr_inv(Scalar2dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> gr_inv(Vec2dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> gr_inv(BiVec2dp<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> gr_inv(PScalar2dp<T> ps)
{
    // grade 3: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> gr_inv(MVec2dp_E<T> const& M)
{
    // grade 0 and 2: no sign change
    return M;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> gr_inv(MVec2dp_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> gr_inv(MVec2dp<T> const& M)
{
    // grade 0 and 2: no sign change
    // grade 1 and 3: sign reversal
    return MVec2dp<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)));
}

////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> rev(Scalar2dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> rev(Vec2dp<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> rev(BiVec2dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> rev(PScalar2dp<T> ps)
{
    // grade 3: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> rev(MVec2dp_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec2dp_E<T>(rev(gr0(M)), rev(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> rev(MVec2dp_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec2dp_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> rev(MVec2dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    return MVec2dp<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)));
}

////////////////////////////////////////////////////////////////////////////////
// regressive reversion operation: reverse applied to the complement + backtrafo
// rrev(A_r) = cmpl((-1)^(r*(r-1)/2) cmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=3, r = 0, 1, 2, 3, ...: - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> rrev(Scalar2dp<T> s)
{
    // grade 0: sign reversal
    return -s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> rrev(Vec2dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> rrev(BiVec2dp<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> rrev(PScalar2dp<T> ps)
{
    // grade 3: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> rrev(MVec2dp_E<T> const& M)
{
    // grade 0: sign reversal
    // grade 2: no sign change
    return MVec2dp_E<T>(rrev(gr0(M)), rrev(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> rrev(MVec2dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dp_U<T>(rrev(gr1(M)), rrev(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> rrev(MVec2dp<T> const& M)
{
    // grade 0: sign reversal
    // grade 1: sign reversal
    // grade 2: no sign change
    // grade 3: no sign change
    return MVec2dp<T>(rrev(gr0(M)), rrev(gr1(M)), rrev(gr2(M)), rrev(gr3(M)));
}

////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> conj(Scalar2dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> conj(Vec2dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> conj(BiVec2dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> conj(PScalar2dp<T> ps)
{
    // grade 3: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> conj(MVec2dp_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec2dp_E<T>(conj(gr0(M)), conj(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> conj(MVec2dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dp_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> conj(MVec2dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    return MVec2dp<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)));
}


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
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot([[maybe_unused]] PScalar2dp<T>,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot(BiVec2dp<T> const& B1,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot(Vec2dp<T> const& v1,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot(Scalar2dp<T> s1, Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp<T> const& A,
                                                         MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // return Scalar2dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c6 * B.c6);
    return Scalar2dp<ctype>(dot(gr0(A), gr0(B)) + dot(gr1(A), gr1(B)) +
                            dot(gr2(A), gr2(B)) + dot(gr3(A), gr3(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp_E<T> const& A,
                                                         MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(dot(gr0(A), gr0(B)) + dot(gr2(A), gr2(B)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp_U<T> const& A,
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
inline constexpr PScalar2dp<std::common_type_t<T, U>> rdot(PScalar2dp<T> ps1,
                                                           PScalar2dp<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps1) * ctype(ps2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> rdot(BiVec2dp<T> const& B1,
                                                           BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(B1.x * B2.x + B1.y * B2.y);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> rdot(Vec2dp<T> const& v1,
                                                           Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(v1.z * v2.z);
}

// returns 0.0 due to degenerate metric with e3^2 = 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> rdot([[maybe_unused]] Scalar2dp<T>,
                                                           [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////

// wedge product extended to fully populated multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> wdg(MVec2dp<T> const& A,
                                                       MVec2dp<U> const& B)
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
inline constexpr Scalar2dp<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2dp<T>,
                                                         [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a pseudoscalar and a bivector
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
wdg([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a bivector and a pseudoscalar
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
wdg([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a pseudoscalar and a vector
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
wdg([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a vector and a pseudoscalar
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
wdg([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between a trivector ps and a scalar s
// wdg(ps,s) = ps*s*e123
// => returns a trivector (scalar s multiple of the trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> wdg(PScalar2dp<T> ps,
                                                          Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// wedge product between a scalar s and a trivector ps
// wdg(s,ps) = s*ps*e123
// => returns a trivector (scalar s multiple of the trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s,
                                                          PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// wedge product between two bivectors
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
wdg([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// wedge product between bivector B and vector v
// wdg(B,v) = gr3(B * v)
// => returns a trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B,
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
inline constexpr PScalar2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v,
                                                          BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-v.x * B.x - v.y * B.y - v.z * B.z);
}

// wedge product with one scalar (returns a scaled bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B,
                                                        Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// wedge product with one scalar (returns a scaled bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s,
                                                        BiVec2dp<U> const& B)
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
inline constexpr BiVec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v1,
                                                        Vec2dp<U> const& v2)
{
    return BiVec2dp<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// wedge product between two scalars (returns a scalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s1, Scalar2dp<U> s2)
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
inline constexpr Line2d<std::common_type_t<T, U>> wdg(Point2d<T> const& p,
                                                      Point2d<U> const& q)
{
    return Line2d<std::common_type_t<T, U>>(p.y - q.y, q.x - p.x, p.x * q.y - p.y * q.x);
}

// expand to a new line with goes through point p and is perpendicular to line l
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Line2d<std::common_type_t<T, U>> expand(Point2d<T> const& p,
                                                         Line2d<U> const& l)
{
    return rweight_expansion(Point2dp<std::common_type_t<T, U>>(p), l);
}

////////////////////////////////////////////////////////////////////////////////
// convenience functions wdg -> join
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> join(BiVec2dp<T> const& B,
                                                           Vec2dp<U> const& v)
{
    return wdg(B, v);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& v,
                                                           BiVec2dp<U> const& B)
{
    return wdg(v, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& v1,
                                                         Vec2dp<U> const& v2)
{
    return wdg(v1, v2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Line2d<std::common_type_t<T, U>> join(Point2d<T> const& p,
                                                       Point2d<U> const& q)
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
inline constexpr MVec2dp<std::common_type_t<T, U>> rwdg(MVec2dp<T> const& A,
                                                        MVec2dp<U> const& B)
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
inline constexpr BiVec2dp<std::common_type_t<T, U>> rwdg(PScalar2dp<T> ps,
                                                         BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * B;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B,
                                                         PScalar2dp<U> ps)
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
inline constexpr Vec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B1,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& B,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& v,
                                                          BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(-v.x * B.x - v.y * B.y - v.z * B.z);
}

// required to be present for dist2dp (to complile, even if not used)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
rwdg([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] Vec2dp<U> const&)
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
inline constexpr Vec2dp<std::common_type_t<T, U>> meet(BiVec2dp<T> const& B1,
                                                       BiVec2dp<U> const& B2)
{
    return rwdg(B1, B2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> meet(Line2d<T> const& l1,
                                                       Line2d<U> const& l2)
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
// operator<<(a,b) = lbulk_contract(a,b) = rwdg( bulk_dual(a), b )
//
////////////////////////////////////////////////////////////////////////////////

// left contraction (A << B) of two fully populated multivectors
// returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> operator<<(MVec2dp<T> const& A,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction - pseudoscalar contracted onto bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction of a bivector B onto a pseudoscalar ps (=trivector)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator<<(BiVec2dp<T> const& B,
                                                             PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(0.0, 0.0, -B.z) * ctype(ps);
}

// left contraction - pseudoscalar contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction of a vector v onto a pseudoscalar ps (=trivector)
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v,
                                                               PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(-v.x, -v.y, 0.0) * ctype(ps);
}

// left contraction - pseudoscalar contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction (s << ps) - scalar s contracted onto pseudoscalar ps
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s,
                                                                 PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// left contraction (B1 << B2) - bivector B1 contracted onto B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(BiVec2dp<T> const& B1,
                                                                BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(B1.z * B2.z);
}

// left contraction - bivector contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
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
inline constexpr Vec2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v,
                                                             BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v.y * B.z, -v.x * B.z, v.x * B.y - v.y * B.x);
}

// left contraction - bivector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction (s << B) - scalar s contracted onto bivector B
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// left contraction (v1 << v2) - vector v1 contracted onto v2
// (=identical to scalar product dot(v1,v2) of two vectors; identical to (v1 >> v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v1,
                                                                Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// left contraction - vector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] Scalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// left contraction (s << v) - scalar s contracted onto vector v
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s,
                                                             Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// left contraction (s1 << s2) - scalar s1 contracted onto scalar s2
// returns a scalar (the product of the scalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator<<(Scalar2dp<T> s1,
                                                                Scalar2dp<U> s2)
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
// operator>>(a,b) = rbulk_contract(a,b) = rwdg( a, bulk_dual(b) )
//
////////////////////////////////////////////////////////////////////////////////

// right contraction (A >> B) of two fully populated multivectors
// returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> operator>>(MVec2dp<T> const& A,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction of pseudoscalar ps (=trivector) by bivector B "ps contracted by B"
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                             BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2dp<ctype>(ctype(0.0), ctype(0.0), -B.z);
}

// right contraction - bivector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction of pseudoscalar ps (=trivector) by vector v (" ps contracted by v")
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                               Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * BiVec2dp<ctype>(v.x, v.y, ctype(0.0));
}

// right contraction -  vector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (ps >> s) - pseudoscalar ps contracted by scalar s
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<T> ps,
                                                                 Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// right contraction - scalar contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dp<T>, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (B1 >> B2) - bivector B1 contracted by bivector B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B1,
                                                                BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(B1.z * B2.z);
}

// right contraction of bivector B by vector v ("B contracted by v")
// => returns a vector
// (identically to cmt(v,B))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2dp<T> const&, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (B >> s) - bivector B contracted by scalar s
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B,
                                                               Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// right contraction - scalar contracted by a bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dp<T>, [[maybe_unused]] BiVec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (v1 >> v2) - vector v1 contracted by vector v2
// (=identical to scalar product dot(v1,v2) of two vectors; identical to (v1 << v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(Vec2dp<T> const& v1,
                                                                Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// right contraction (v >> s) - vector v contracted by scalar s
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator>>(Vec2dp<T> const& v,
                                                             Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// right contraction - scalar contracted by a vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2dp<T>, [[maybe_unused]] Vec2dp<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0);
}

// right contraction (s1 >> s2) - scalar s1 contracted by scalar s2
// returns a scalar (the product of the scalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator>>(Scalar2dp<T> s1,
                                                                Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> cmt(MVec2dp<T> const& A,
                                                       MVec2dp<U> const& B)
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
inline constexpr BiVec2dp<std::common_type_t<T, U>> cmt(BiVec2dp<T> const& B1,
                                                        BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(-B1.y * B2.z + B1.z * B2.y, B1.x * B2.z - B1.z * B2.x, 0.0);
}

// cmt(B,v) = -cmt(v,B)
// identical to (v << B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> cmt(BiVec2dp<T> const& B,
                                                      Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(B.z * v.y, -B.z * v.x, -B.x * v.y + B.y * v.x);
}

// cmt(v,B) = -cmt(B,v)
// identcial to (B >> v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> cmt(Vec2dp<T> const& v,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(-v.y * B.z, v.x * B.z, -v.x * B.y + v.y * B.x);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> cmt(Vec2dp<T> const& v1,
                                                        Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(v1.y * v2.z - v1.z * v2.y, -v1.x * v2.z + v1.z * v2.x,
                           v1.x * v2.y - v1.y * v2.x);
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
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
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
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
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
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
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
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
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
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
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
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
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
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<U> const& A,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<U> const& A,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& M,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& M,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
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
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator*([[maybe_unused]] PScalar2dp<T>, [[maybe_unused]] PScalar2dp<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0); // trivectors in pga2dp contain e3^2 = 0
}

// geometric product ps * B of a trivector ps (=2dp pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                            BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.z);
}

// geometric product B * ps of bivector A multiplied by a trivector ps from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                            PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.z) * ctype(ps);
}

// geometric product ps * v of a trivector ps (=2dp pseudoscalar) multiplied from the
// left to the vector v trivector * vector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                              Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(ps) * BiVec2dp<ctype>(v.x, v.y, ctype(0.0));
}

// geometric product v * ps of a vector a multiplied with a trivector ps from the right
// vector * trivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                              PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(-v.x, -v.y, ctype(0.0)) * ctype(ps);
}

// geometric product ps * s of a trivector ps and a scalar s
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                                Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) *
                             ctype(s)); // scalar multiplication with a trivector
}

// geometric product s * ps of a scalar s and a trivector ps
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s,
                                                                PScalar2dp<U> ps)
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B1,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>((B >> v), wdg(v, B));
}

// (geometric) product B * s of bivector B multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                              Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// (geometric) product s * B of a scalar s multiplied from the left to the bivector B
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(Vec2dp<T> const& a,
                                                               Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(dot(a, b), wdg(a, b));
}

// (geometric) product v * s of vector v multiplied with scalar s from the right
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                            Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// (geometric) product s * v of scalar s multiplied from the left to vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s,
                                                            Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> s1,
                                                               Scalar2dp<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(s1) * ctype(s2));
}

////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> rgpr(MVec2dp<T> const& A,
                                                        MVec2dp<U> const& B)
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& A,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& A,
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
inline constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& M,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& M,
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(MVec2dp_E<T> const& M,
                                                          Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(-M.c1 * v.x - M.c2 * v.y - M.c3 * v.z),
        BiVec2dp<ctype>(M.c2 * v.z, -M.c1 * v.z, M.c0 * v.z + M.c1 * v.y - M.c2 * v.x));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(MVec2dp_U<T> const& M,
                                                          Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(M.c1 * v.z - M.c2 * v.y + M.c3 * v.x,
                                          -M.c0 * v.z + M.c2 * v.x + M.c3 * v.y,
                                          M.c3 * v.z),
                            PScalar2dp<ctype>(-M.c2 * v.z));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B1,
                                                          BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(-B1.y * B2.z + B1.z * B2.y,
                                          B1.x * B2.z - B1.z * B2.x,
                                          -B1.x * B2.y + B1.y * B2.x),
                            PScalar2dp<ctype>(B1.x * B2.x + B1.y * B2.y));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> rgpr(BiVec2dp<T> const& B,
                                                          Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(-B.x * v.x - B.y * v.y - B.z * v.z),
        BiVec2dp<ctype>(B.y * v.z, -B.x * v.z, B.x * v.y - B.y * v.x));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> rgpr(Vec2dp<T> const& v1,
                                                          Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(v1.y * v2.z - v1.z * v2.y, -v1.x * v2.z + v1.z * v2.x, ctype(0.0)),
        PScalar2dp<ctype>(-v1.z * v2.z));
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

////////////////////////////////////////////////////////////////////////////////
// 2dp euclidean distance
////////////////////////////////////////////////////////////////////////////////

// returns the euclidean distance between objects as homogeneous magnitude
template <typename arg1, typename arg2> DualNum2dp<value_t> dist2dp(arg1&& a, arg2&& b)
{
    if (gr(a) + gr(b) == 3) {
        return DualNum2dp<value_t>(rwdg(a, b), weight_nrm(wdg(a, att(b))));
    }
    else {
        return DualNum2dp<value_t>(bulk_nrm(att(wdg(a, b))), weight_nrm(wdg(a, att(b))));
    }
}

////////////////////////////////////////////////////////////////////////////////
// 2dp motor operations (translation and rotation)
//
// Every motor in pga2dp is an odd grade multivector MVec2dp_U.
//
// A proper isometry in 2dp has a fixed point p = p.x e1 + p.y e2 + p.z e3
// around which a rotation occurs with an angle phi.
//
// Every motor has the form: M = p sin(phi) + e321 cos(phi).
//
// (derived from the exponential function with respect to the regressive
// geometric product)
////////////////////////////////////////////////////////////////////////////////

// create a (unitized) motor from a fixed point and a turning angle
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> motor(Vec2dp<T> const& p, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = 0.5 * theta;
    return unitize(MVec2dp_U<ctype>(Vec2dp<ctype>(p * std::sin(half_angle)),
                                    PScalar2dp<ctype>(std::cos(half_angle))));
}

// create a translation motor from a direction vector (given as a Vec2dp)
// move in direction and by length of direction vector (lenght = its bulk_nrm)
// ATTENTION: the direction is assumed to be a direction vector, i.e. with z == 0
//            the z-component is ignored and only the x- and y-components are used
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> motor(Vec2dp<T> const& direction)
{
    return MVec2dp_U<T>(0.5 * Vec2dp<T>(-direction.y, direction.x, T(0.0)),
                        PScalar2dp<T>(1.0));
}

// create a translation motor from a direction vector (given as a Vector2d)
// move in direction and by length of direction vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> motor(Vector2d<T> const& direction)
{
    return MVec2dp_U<T>(0.5 * Vec2dp<T>(-direction.y, direction.x, T(0.0)),
                        PScalar2dp<T>(1.0));
}

// create a (unitized) motor directly from two (potentially intersecting) lines
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>>
motor_from_lines(BiVec2dp<T> const& B1, BiVec2dp<U> const& B2)
{
    // take lines (=bivectors) as input and return a motor R
    // 1st apply reflection across line B1, then across B2 to get a motor that rotates
    // (or translates) around the intersection point of lines B1 and B2
    //
    // for use of motor R either directly on object u (inefficient):
    //     auto v_moved = gr1( rgpr(rgpr(R, v), rrev(R)) );
    // or
    //     auto B_moved = gr2( rgpr(rgpr(R, B), rrev(R)) );
    // or
    //                                   // optimized for reduced effort
    //     auto v_moved = move2dp(v,R);  // moves v according to the motor R
    //     auto B_moved = move2dp(B,R);  // moves B according to the motor R
    //
    return unitize(rgpr(B2, B1)); // based on the regressive geometric product
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp(Vec2dp<T> const& v,
                                                          MVec2dp_U<U> const& R)
{
    // assumes: motor R is unitized

    // moves v (a vector representing a projective point) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(gr1(rgpr(rgpr(R, v), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp(BiVec2dp<T> const& B,
                                                            MVec2dp_U<U> const& R)
{
    // assumes: motor R is unitized

    // moves B (a bivector representing a line) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(rgpr(rgpr(R, B), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp_opt(Vec2dp<T> const& v,
                                                              MVec2dp_U<U> const& R)
{
    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    //
    // (could potentially be further optimized by exporting the matrix-representation
    // if many transformations should be calculated using the same rotor as
    // v' = matrix * v)

    using ctype = std::common_type_t<T, U>;
    auto k02 = R.c0 * R.c2;
    auto k03 = R.c0 * R.c3;
    auto k12 = R.c1 * R.c2;
    auto k13 = R.c1 * R.c3;
    auto k22 = R.c2 * R.c2;
    auto k23 = R.c2 * R.c3;
    auto k33 = R.c3 * R.c3;
    return Vec2dp<ctype>((-k22 + k33) * v.x - 2.0 * k23 * v.y + 2.0 * (k02 + k13) * v.z,
                         2.0 * k23 * v.x + (-k22 + k33) * v.y + 2.0 * (-k03 + k12) * v.z,
                         (k22 + k33) * v.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp_opt(BiVec2dp<T> const& B,
                                                                MVec2dp_U<U> const& R)
{
    // moves B (a bivector representing a line) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    //
    // (could potentially be further optimized by exporting the matrix-representation
    // if many transformations should be calculated using the same rotor as
    // B' = matrix * B)

    using ctype = std::common_type_t<T, U>;
    auto k02 = R.c0 * R.c2;
    auto k03 = R.c0 * R.c3;
    auto k12 = R.c1 * R.c2;
    auto k13 = R.c1 * R.c3;
    auto k22 = R.c2 * R.c2;
    auto k23 = R.c2 * R.c3;
    auto k33 = R.c3 * R.c3;
    return BiVec2dp<ctype>(
        (-k22 + k33) * B.x + 2.0 * (-k23) * B.y, 2.0 * k23 * B.x + (-k22 + k33) * B.y,
        2.0 * (k02 - k13) * B.x + 2.0 * (k12 + k03) * B.y + (k22 + k33) * B.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp_opt2(Vec2dp<T> const& v,
                                                               MVec2dp_U<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = R.c1 * v.z - R.c2 * v.y + R.c3 * v.x;
    ctype k2 = -R.c0 * v.z + R.c2 * v.x + R.c3 * v.y;
    ctype k3 = R.c3 * v.z;
    ctype k4 = -R.c2 * v.z;
    return Vec2dp<ctype>(k1 * R.c3 - k2 * R.c2 + k3 * R.c1 - k4 * R.c0,
                         k1 * R.c2 + k2 * R.c3 - k3 * R.c0 - k4 * R.c1,
                         k3 * R.c3 - k4 * R.c2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp_opt2(BiVec2dp<T> const& B,
                                                                 MVec2dp_U<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = -R.c0 * B.x - R.c1 * B.y - R.c2 * B.z;
    ctype k2 = -R.c2 * B.y + R.c3 * B.x;
    ctype k3 = R.c2 * B.x + R.c3 * B.y;
    ctype k4 = R.c0 * B.y - R.c1 * B.x + R.c3 * B.z;
    return BiVec2dp<ctype>(k2 * R.c3 - k3 * R.c2, k2 * R.c2 + k3 * R.c3,
                           -k1 * R.c2 - k2 * R.c1 + k3 * R.c0 + k4 * R.c3);
}

////////////////////////////////////////////////////////////////////////////////
// 2dp complement operations
// (the concept of complement is defined w.r.t. the outer product)
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^2 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors which are in the k-blade u with
// the basis vectors which are NOT contained in the k-blade u
// and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  lcmpl(u) ^ u  = I_2dp = e1^e2^e3
// right complement: u ^ rcmpl(u)  = I_2dp = e1^e2^e3
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs

// complement operation with e1^e2^e3 as the pseudoscalar

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> cmpl(Scalar2dp<T> s)
{
    // u ^ cmpl(u) = e3^e2^e1
    // u = 1:
    //     1 ^ cmpl(u) = e2^e2^e1 => cmpl(u) = s e3^e2^e1
    return PScalar2dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> cmpl(Vec2dp<T> const& v)
{
    // u ^ compl(u) = e3^e2^e1
    // u = v.x e1 + v.y e2 + v.z e3:
    //     u ^ cmpl(u) = e3^e2^e1 =>
    //     u = e1 => cmpl(u) = -e23
    //     u = e2 => cmpl(u) = -e31
    //     u = e3 => cmpl(u) = -e12
    return BiVec2dp<T>(-v.x, -v.y, -v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> cmpl(BiVec2dp<T> const& B)
{
    // u ^ compl(u) = e3^e2^e1
    // u = b.x e23 + b.y e31 + b.z e12:
    //     u ^ cmpl(u) = e3^e2^e1 =>
    //     u = e23 => cmpl(u) = -e1
    //     u = e31 => cmpl(u) = -e2
    //     u = e12 => cmpl(u) = -e3
    return Vec2dp<T>(-B.x, -B.y, -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> cmpl(PScalar2dp<T> ps)
{
    // u ^ compl(u) = e3^e2^e1
    // u = e3^e2^e1:
    //     e3^e2^e1 ^ cmpl(u) = e3^e2^e1 => cmpl(u) = ps 1
    return Scalar2dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> cmpl(MVec2dp_E<T> const& M)
{
    // use the component complements directly
    return MVec2dp_U<T>(cmpl(gr2(M)), cmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> cmpl(MVec2dp_U<T> const& M)
{
    // use the component complements directly
    return MVec2dp_E<T>(cmpl(gr3(M)), cmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> cmpl(MVec2dp<T> const& M)
{
    // use the component complements directly
    return MVec2dp<T>(cmpl(gr3(M)), cmpl(gr2(M)), cmpl(gr1(M)), cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// projections, rejections
////////////////////////////////////////////////////////////////////////////////

// TODO: check whether the vector-vector formulas make sense at all, since they model
//       the representational space an not the modelled subspace
//       => potentially only meaningful, when implemented based on type vector2dp

// projection of a vector v1 onto vector v2
// returns component of v1 parallel to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v1,
                                                               Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of vector v1 from a vector v2
// returns component of v1 perpendicular to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v1,
                                                              Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec2dp<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}


// projection of a vector v onto a bivector B (a line)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(rwdg(B, wdg(v, weight_dual(B)))); // ortho_proj
}

// rejection of vector v from a bivector B (a line)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v,
                                                              BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v - project_onto(v, B));
}


////////////////////////////////////////////////////////////////////////////////
// reflections
////////////////////////////////////////////////////////////////////////////////

// reflect a vector u in an arbitrary bivector, i.e. a line
// B must be unitized
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reflect_on(Vec2dp<T> const& v,
                                                             BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(-gr1(rgpr(rgpr(B, v), B)));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling lines)
// B must be unitized
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> reflect_on(BiVec2dp<T> const& UB,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(rgpr(rgpr(B, UB), B)));
}

} // namespace hd::ga::pga