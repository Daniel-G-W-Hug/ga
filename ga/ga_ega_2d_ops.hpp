#pragma once

// author: Daniel Hug, 2024 & 2025

#include <algorithm> // std::max, std::clamp
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <numbers>   // math constants like pi
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "detail/ga_mvec2d.hpp" // inclusion of multivector imports all component types

namespace hd::ga::ega {

////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> gr_inv(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> gr_inv(Vec2d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> gr_inv(PScalar2d<T> ps)
{
    // grade 2: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> gr_inv(MVec2d_E<T> const& M)
{
    // grade 0 and 2: no sign change
    return M;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> gr_inv(MVec2d<T> const& M)
{
    // grade 0 and 2: no sign change
    // grade 1: sign reversal
    return MVec2d<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)));
}

////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> rev(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> rev(Vec2d<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> rev(PScalar2d<T> ps)
{
    // grade 2: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> rev(MVec2d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign reversal
    return MVec2d_E<T>(rev(gr0(M)), rev(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> rev(MVec2d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign reversal
    return MVec2d<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)));
}

////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> conj(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> conj(Vec2d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> conj(PScalar2d<T> ps)
{
    // grade 2: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> conj(MVec2d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign reversal
    return MVec2d_E<T>(conj(gr0(M)), conj(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> conj(MVec2d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    return MVec2d<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)));
}


////////////////////////////////////////////////////////////////////////////////
// dot product (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> dot(MVec2d<T> const& A,
                                                        MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> dot(PScalar2d<T> ps1,
                                                        PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// return dot-product of two vectors in G<2,0,0>
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> dot(Vec2d<T> const& v1,
                                                        Vec2d<U> const& v2)
{
    // definition: dot(v1, v2) = (v1)^T g_12 v2 with the metric g_12
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1
    // as diagonal elements of g_12
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> dot(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// wedge product (= outer product)
////////////////////////////////////////////////////////////////////////////////

// wedge product between two fully populated multivectors A and B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between multivectors A and even grade multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& A,
                                                      MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c1 * B.c0;
    ctype c2 = A.c2 * B.c0;
    ctype c3 = A.c0 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between even grade multivector A and multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1;
    ctype c2 = A.c0 * B.c2;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between multivector M and pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M,
                                                         PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(M.c0 * ctype(ps));
}

// wedge product between pseudoscalar ps and multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps,
                                                         MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * M.c0);
}

// wedge product M ^ v between multivector M and vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M,
                                                      Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0.0;
    ctype c1 = M.c0 * v.x;
    ctype c2 = M.c0 * v.y;
    ctype c3 = M.c1 * v.y - M.c2 * v.x;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product v ^ M between vector v and multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v,
                                                      MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0.0;
    ctype c1 = v.x * M.c0;
    ctype c2 = v.y * M.c0;
    ctype c3 = v.x * M.c2 - v.y * M.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between multivector M and scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// wedge product between scalar s and multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// wedge product between two even grade multivectors A and B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& A,
                                                        MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// wedge product between even grade multivector M and pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M,
                                                         PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(M.c0 * ctype(ps));
}

// wedge product between pseudoscalar ps and even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps,
                                                         MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * M.c0);
}

// wedge product between even grade multivector M and vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M,
                                                        Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(M.c0 * v.x, M.c0 * v.y);
}

// wedge product between vector v and even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> wdg(Vec2d<T> const& v,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x * M.c0, v.y * M.c0);
}

// wedge product between even grade multivector M and scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M,
                                                        Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// wedge product between scalar s and even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> wdg(Scalar2d<T> s,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// wedge product of two pseudoscalars
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2d<T>,
                                                        [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// wedge product of a vector with a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec2d<T> const&,
                                                        [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// wedge product of a pseudoscalar with a vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2d<T>,
                                                        [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// wedge product with a pscalar (returns a scaled pscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// wedge product (returns a bivector)
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v1,
                                                         Vec2d<U> const& v2)
{
    return PScalar2d<std::common_type_t<T, U>>(v1.x * v2.y - v1.y * v2.x);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// wedge product between two scalars (returns a scalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s1, Scalar2d<U> s2)
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
// rwdg(ul, ur) = lcmpl(wdg(rcmpl(ul),rcmpl(ur))) = lcmpl(rcmpl(ul)^rcmpl(ur))
//              = rcmpl(wdg(lcmpl(ul),lcmpl(ur))) = rcmpl(lcmpl(ul)^lcmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// regressive wedge product extended to a fully populated multivectors
// => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> rwdg(MVec2d<T> const& A,
                                                       MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    ctype c1 = A.c1 * B.c3 + A.c3 * B.c1;
    ctype c2 = A.c2 * B.c3 + A.c3 * B.c2;
    ctype c3 = A.c3 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// regressive wedge product between to pseudoscalars
// => returns a scaled pseudoscalar (=bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps1,
                                                          PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// regressive wedge product between a pseudoscalar ps (=bivector) and a vector v
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2d<ctype>(v.x, v.y);
}

// regressive wedge product between a vector v and a pseudoscalar ps (=bivector)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x, v.y) * ctype(ps);
}

// regressive wedge product between a pseudoscalar ps (=bivector) and a scalar s
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps) * ctype(s));
}

// regressive wedge product between a scalar s and a pseudoscalar ps (=bivector)
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s) * ctype(ps));
}

// regressive wedge product between a vector v1 and a vector v2
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v1,
                                                         Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.y - v1.y * v2.x);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B
//
// The resulting object is A lies in B and is perpendicular to A
//
// L. Dorst: The contraction A << B of an a-blade A onto a b-blade B is a sub-blade
// of B of grade b-a which is perpendicular to A, and linear in both arguments
////////////////////////////////////////////////////////////////////////////////

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& A,
                                                             MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype c1 = A.c0 * B.c1 + A.c2 * B.c3;
    ctype c2 = A.c0 * B.c2 - A.c1 * B.c3;
    ctype c3 = A.c0 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& A,
                                                             MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype c1 = A.c2 * B.c1;
    ctype c2 = -A.c1 * B.c1;
    ctype c3 = A.c0 * B.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& A,
                                                             MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c3;
    ctype c1 = A.c0 * B.c1;
    ctype c2 = A.c0 * B.c2;
    ctype c3 = A.c0 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (M << ps) of mulivector A onto pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<U> const& M,
                                                             PScalar2d<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c3 * ps;
    ctype c1 = M.c2 * ps;
    ctype c2 = -M.c1 * ps;
    ctype c3 = M.c0 * ps;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (ps << M) of a pseudoscalar ps onto mulivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps,
                                                               MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ps * M.c3);
}

// left contraction (M << v) of multivector M onto vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& M,
                                                             Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * v.x + M.c2 * v.y;
    ctype c1 = M.c0 * v.x;
    ctype c2 = M.c0 * v.y;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (v << M) of vector v onto multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v,
                                                             MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = v.x * M.c1 + v.y * M.c2;
    ctype c1 = v.y * M.c3;
    ctype c2 = -v.x * M.c3;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (M << s) of multivector M onto scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& M,
                                                               Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c0 * s);
}

// left contraction (s << M) of scalar s onto multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s,
                                                             MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;

    return ctype(s) * M;
}

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& A,
                                                               MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1;
    ctype c1 = A.c0 * B.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// left contraction (M << ps) of even grade mulivector M with pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                               PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * ps;
    ctype c1 = M.c0 * ps;
    return MVec2d_E<ctype>(c0, c1);
}

// left contraction (ps << M) of pseudoscalar ps with an even grade mulivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps,
                                                               MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ps * M.c1);
}

// left contraction (M << v) - "M contracted onto v"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                            Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(M.c0 * v.x, M.c0 * v.y);
}

// left contraction (v << M) - "v contracted onto M"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v,
                                                               MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.y * M.c1, -v.x * M.c1);
}

// left contraction (M << s) of even grade mulivector M with scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                               Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c0 * s);
}

// left contraction (s << M) of a scalar s with an even grade mulivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(Scalar2d<T> s,
                                                               MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// left contraction - pseudoscalar ps1 contracted onto pseudoscalar ps2
// returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps1,
                                                               PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ps1 * ps2);
}

// left contraction -  pseudoscalar contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2d<T>, [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (v << ps) - vector v contracted onto pseudoscalar ps
// (identical with the geometric product v * B for this case)
// returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v,
                                                            PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.y, -v.x) * ctype(ps);
}

// left contraction - pseudoscalar contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar2d<T>, [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (s << ps) - scalar s contracted onto pseudoscalar ps
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s,
                                                                PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// left contraction (v1 << v2) - vector v1 contracted onto vector v2
// (=identical to scalar product dot(v1,v2) of two vectors and identical to (v1 >> v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v1,
                                                               Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// left contraction - vector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec2d<T> const&, [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (s << v) - scalar s contracted onto vector v
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s,
                                                            Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// left contraction (s1 << s2) - scalar s1 contracted onto scalar s2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s1,
                                                               Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B
//
// The resulting object lies in A and is perpendicular to B
////////////////////////////////////////////////////////////////////////////////

// right contraction (A >> B)  - "A contracted by B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<T> operator>>(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype c1 = A.c1 * B.c0 - A.c3 * B.c2;
    ctype c2 = A.c2 * B.c0 + A.c3 * B.c1;
    ctype c3 = A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (A >> B)  - "A contracted by B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<T> operator>>(MVec2d<T> const& A, MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype c1 = A.c1 * B.c0;
    ctype c2 = A.c2 * B.c0;
    ctype c3 = A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (A >> B)  - "A contracted by B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<T> operator>>(MVec2d_E<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c3;
    ctype c1 = -A.c1 * B.c2;
    ctype c2 = A.c1 * B.c1;
    ctype c3 = A.c1 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}
// right contraction (M >> ps) of a multivector M by a pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& M,
                                                               PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c3 * ps);
}

// right contraction (ps >> M) of a pseudoscalar ps by a multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps,
                                                             MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = ps * M.c3;
    ctype c1 = -ps * M.c2;
    ctype c2 = ps * M.c1;
    ctype c3 = ps * M.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (M >> v) of multivector A contracted by vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<T> operator>>(MVec2d<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * v.x + M.c2 * v.y;
    ctype c1 = -M.c3 * v.y;
    ctype c2 = M.c3 * v.x;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (v >> M) of vector v contracted by multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<T> operator>>(Vec2d<T> const& v, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = v.x * M.c1 + v.y * M.c2;
    ctype c1 = v.x * M.c0;
    ctype c2 = v.y * M.c0;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (M >> s) of multivector M by scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& M,
                                                             Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// right contraction (s >> M) of scalar s by multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s,
                                                               MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(s * M.c0);
}

// right contraction (A >> B) of an even grade multivector A by even grade multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& A,
                                                               MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1;
    ctype c1 = A.c1 * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// right contraction (M >> ps) of an even grade multivector M by a pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                               PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c1 * ps);
}

// right contraction (ps >> M) of pseudoscalar ps by even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps,
                                                               MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = ps * M.c1;
    ctype c1 = ps * M.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// left contraction (M >> v) - "M contracted by v"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                            Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-M.c1 * v.y, M.c1 * v.x);
}

// left contraction (v >> M) - "v contracted by M"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v,
                                                               MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x * M.c0, v.y * M.c0);
}

// right contraction (M >> s) of even grade multivector M by scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                               Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// right contraction (s >> M) of scalar s by even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s,
                                                               MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(s * M.c0);
}

// right contraction - pseudoscalar s1 contracted by pseudoscalar s2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps1,
                                                               PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// right contraction (ps >> v) - pseudoscalar ps contracted by vector v
// (identical with the geometric product ps * v for this case)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps,
                                                            Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2d<ctype>(-v.y, v.x);
}

// right contraction -  vector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec2d<T> const&, [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// right contraction (ps >> s) - pseudoscalar ps contracted by scalar s
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps,
                                                                Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// right contraction - scalar contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2d<T>, [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// right contraction (v1 >> v2) - vector v2 taken out of a vector v1
// (=identical to scalar product dot(v1,v2) of two vectors and identical to (v1 << v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v1,
                                                               Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

// right contraction (v >> s) - vector v contracted by scalar s
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator>>(Vec2d<U> const& v,
                                                            Scalar2d<T> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// right contraction - scalar contracted by a vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar2d<T>, [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (s1 >> s2) - scalar s1 contracted by scalar s2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s1,
                                                               Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A * B for fully populated 2d multivectors
// Expensive! - Don't use if you don't have to! (16x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors, e.g.:
// a*b = dot(a,b) + wdg(a,b) = gr0(a*b) + gr2(a*b) (vector * vector = scalar + bivector)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // geometric product A * B of two fully populated multivectors
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// geometric product A * B for a multivector A with an even grade multivector B
// multivector * even grade multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c0 * B.c0 - A.c3 * B.c1, A.c1 * B.c0 - A.c2 * B.c1,
                         A.c1 * B.c1 + A.c2 * B.c0, A.c0 * B.c1 + A.c3 * B.c0);
}

// geometric product A * B for an even grade multivector A with a full 2d multivector B
// from the right
// even grade multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c0 * B.c0 - A.c1 * B.c3, A.c0 * B.c1 + A.c1 * B.c2,
                         A.c0 * B.c2 - A.c1 * B.c1, A.c0 * B.c3 + A.c1 * B.c0);
}

// geometric product A * ps of a multivector A with the pseudoscalar ps (=bivector)
// multivector * pseudoscalar (=bivector) => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-A.c3, -A.c2, A.c1, A.c0) * ctype(ps);
}

// geometric product ps * B of pseudoscalar ps (=bivector) with the multivector B
// 2d pseudoscalar (=bivector) * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec2d<ctype>(-B.c3, B.c2, -B.c1, B.c0);
}

// geometric product A * v for a multivector A with a vector v
// multivector * vector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c1 * v.x + A.c2 * v.y, A.c0 * v.x + A.c3 * v.y,
                         A.c0 * v.y - A.c3 * v.x, A.c1 * v.y - A.c2 * v.x);
}

// geometric product v * B for a vector v with a full 2d multivector B
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(v.x * B.c1 + v.y * B.c2, v.x * B.c0 - v.y * B.c3,
                         v.x * B.c3 + v.y * B.c0, v.x * B.c2 - v.y * B.c1);
}

// geometric product A * s of a multivector A multiplied from the right by the scalar s
// multivector * scalar => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return A * ctype(s);
}

// geometric product s * B of a scalar s multiplied to the multivector B from left
// scalar * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// geometric product A * B for two multivectors from the even subalgebra (2d case)
// even grade multivector * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                              MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(A.c0 * B.c0 - A.c1 * B.c1, A.c0 * B.c1 + A.c1 * B.c0);
}

// geometric product A * ps of an even grade multivector A multiplied from the right
// with the 2d pseudoscalar (=bivector) ps
// even grade multivector * 2d pseudoscalar (=bivector) => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                              PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(-A.c1, A.c0) * ctype(ps);
}

// geometric product ps * B of a bivector A multiplied from the left
// to the multivector from the even subalgebra B (MVec2d_E)
// bivector * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(PScalar2d<T> ps,
                                                              MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec2d_E<ctype>(-B.c1, B.c0);
}

// geometric product A * v for a multivector from the even subalgebra A with a vector v
// even grade multivector * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                           Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(A.c0 * v.x + A.c1 * v.y, A.c0 * v.y - A.c1 * v.x);
}

// geometric product v * B of a vector v with an even grade multivector B
// vector * even grade multivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v,
                                                           MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x * B.c0 - v.y * B.c1, v.x * B.c1 + v.y * B.c0);
}

// geometric product A * s of an even grade multivector A multiplied with a scalar s
// even grade multivector * scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                              Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return A * ctype(s);
}

// geometric product s * B of a scalar with an even grade multivector B
// scalar * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Scalar2d<T> s,
                                                              MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// geometric product ps * ps of two pseudoscalars
// pseudoscalar * pseudoscalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps1,
                                                              PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(-ctype(ps1) * ctype(ps2));
}

// geometric product ps * v of a pseudoscalar s with the vector b
// pseudoscalar (=bivector) * vector => vector
// this multiplication turns the vector by -90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps,
                                                           Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2d<ctype>(v.y, -v.x);
}

// geometric product v * ps of the vector v multiplied with the pseudoscalar ps
// vector * pseudoscalar (=bivector) => vector
// this multiplication turns the vector by +90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v,
                                                           PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-v.y, v.x) * ctype(ps);
}

// geometric product ps * s of pseudscalar ps and scalar s
// pseudoscalar (=bivector) * scalar => pseudoscalar (=bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps,
                                                               Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// geometric product s * ps of scalar s and pseudoscalar ps
// scalar * pseudoscalar (=bivector) => pseudoscalar (=bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s,
                                                               PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// geometric product ab for two vectors (returns a multivector of the even subalgebra)
// v1*v2 = dot(v1,v2) + wdg(v1,v2) = gr0(v1*v2) + gr2(v1*v2)
// => vector * vector = scalar + bivector
//
// HINT: if a full 2d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Vec2d<T> const& v1,
                                                              Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(dot(v1, v2), wdg(v1, v2));
}

// geometric product v * s of vector v and scalar s
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<U> const& v,
                                                           Scalar2d<T> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// geometric product s * v of scalar s and vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s,
                                                           Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s1,
                                                              Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors
// w.r.t. the geometric product:
// for k-blades: A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> inv(Scalar2d<T> s)
{
    T sq_n = nrm_sq(s);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("scalar norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;

    return Scalar2d<T>(rev(s) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> inv(Vec2d<T> const& v)
{
    T sq_n = nrm_sq(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;
    return Vec2d<T>(rev(v) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> inv(PScalar2d<T> ps)
{
    T sq_n = nrm_sq(ps);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector

    return PScalar2d<T>(rev(ps) * inv);
}

// return the multiplicative inverse of the complex number (inv(z) = 1/nrm_sq(z)*rev(z))
// with rev(z) being the complex conjugate
template <typename T> inline constexpr MVec2d_E<T> inv(MVec2d_E<T> const& E)
{
    T sq_n = nrm_sq(E);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error(
            "complex norm of even grade multivector too small for inversion " +
            std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return MVec2d_E<T>(rev(E) * inv);
}

template <typename T> inline constexpr MVec2d<T> inv(MVec2d<T> const& M)
{
    // inv(M) = 1/( M*conj(M) ) * conj(M)  with M*conj(M) being a scalar value
    // from manual calculation of M*conj(M) in 2d:
    T m_conjm = M.c0 * M.c0 + M.c3 * M.c3 - nrm_sq(Vec2d<T>(M.c1, M.c2));
    //
    // alternative, but with slightly more computational effort:
    // T m_conjm = gr0(v * conj(v));

#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(m_conjm) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("multivector norm too small for inversion " +
                                 std::to_string(m_conjm) + "\n");
        // example: MVec2D(1,1,1,1) is not invertible
    }
#endif
    T inv = T(1.0) / m_conjm; // inverse of squared norm for a vector
    return MVec2d<T>(conj(M) * inv);
}


////////////////////////////////////////////////////////////////////////////////
// angle operations
////////////////////////////////////////////////////////////////////////////////

// return the angle between two vectors
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    using std::numbers::pi;

    ctype nrm_prod = nrm(v1) * nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif

    // std::clamp must be used to take care of numerical inaccuracies
    auto cos_angle = std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
    auto sin_angle = std::clamp(ctype(wdg(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
    // wdg() in 2d contains magnitude and orientation, but works this easy only in 2d,
    // because it is already a scalar value
    // (for 3d to be as effective, the 3d vectors would need to be transformed
    //  to a plane, the angle measured w.r.t. to the pseudoscalar of the plane)

    // fmt::println("   c = {: .4f}, s = {: .4f}, wdg = {: .4f}, nrm_wdg = {: .4f}",
    //              cos_angle, sin_angle, wdg(v1, v2), nrm(wdg(v1, v2)));

    if (cos_angle >= 0.0) {
        // quadrant I or IV
        return std::asin(sin_angle);
    }
    else if (cos_angle < 0.0 && sin_angle >= 0.0) {
        // quadrant II
        return pi - std::asin(sin_angle);
    }
    else {
        // cos_angle < 0.0 && sin_angle < 0.0)
        // quadrant III
        return -pi - std::asin(sin_angle);
    }
}

// MVec2d_E<T> is used to model complex numbers using multivectors from the
//             even subalgebra
// returns the angle of the complex number w.r.t. the real axis
// range of angle: -pi <= angle <= pi
template <typename T>
    requires(std::floating_point<T>)
inline constexpr T angle_to_re(MVec2d_E<T> const& v)
{
    using std::numbers::pi;
    if (v.c0 > 0.0) {
        // quadrant I & IV
        return std::atan(v.c1 / v.c0);
    }
    if (v.c0 < 0.0 && v.c1 >= 0.0) {
        // quadrant II
        return std::atan(v.c1 / v.c0) + pi;
    }
    if (v.c0 < 0.0 && v.c1 < 0.0) {
        // quadrant III
        return std::atan(v.c1 / v.c0) - pi;
    }
    if (v.c0 == 0.0) {
        // on y-axis
        if (v.c1 > 0.0) return pi / 2.0;
        if (v.c1 < 0.0) return -pi / 2.0;
    }
    return 0.0; // zero as input => define 0 as corresponding angle
}

////////////////////////////////////////////////////////////////////////////////
// 2d rotation operations
////////////////////////////////////////////////////////////////////////////////

// exponential function for setup of complex numbers and rotations
// as geometric multivectors with a scalar and a bivector part
//
// r = 1 is the vector length of the complex number in polar form
// theta is the bivector angle (i.e. a multiple of the bivector I_2d)
// such that uv = r exp(I_2d, theta) = a + I_2d b
// with r = |u| |v| = sqrt(a^2 + b^2) = 1
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> exp([[maybe_unused]] PScalar2d<T> I,
                                                        U theta)
{
    // PScalar2d<T> is just used here for a unique overload of exp()
    // and to make the function signature similar to the 3D case
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(Scalar2d<ctype>(std::cos(theta)),
                           PScalar2d<ctype>(std::sin(theta)));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//         - a 2d pseudoscalar representing the plane of 2d space
//         - a rotation angle in the plane of 2d space
// Output:
//         - a rotor representing the requested rotation,
//           when applying the sandwich product with the rotor as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// implemented here to make it formally the same with the 3d case (and potentially higher
// dimensional applications). In 2d the rotation can be directly expressed with less
// effort as
//
// exp(I_2d, -theta) * v = v * exp(I_2d, theta)  to express a 2d rotation of vector v by
//                                               the angle theta
//
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> rotor([[maybe_unused]] PScalar2d<T> I,
                                                          U theta)
{
    // PScalar2d<T> is just used here to be able to overload exp
    // and to make the function similar to the 3D case
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec2d_E<ctype>(Scalar2d<ctype>(std::cos(half_angle)),
                           PScalar2d<ctype>(std::sin(half_angle)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> rotate(Vec2d<T> const& v,
                                                        MVec2d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(rotor * v * rev(rotor));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> rotate(MVec2d<T> const& M,
                                                         MVec2d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(rotor * M * rev(rotor));
}


////////////////////////////////////////////////////////////////////////////////
// 2d complement operations (impact on basis vector exclusively!)
// (the concept of complement is defined w.r.t. the outer product)
// it depends on the pseudoscalar of the space
// here implemented for I_2d = e1^e2 = e12
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^2 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors (exclusively, no impact on magnitude),
// which are in the k-blade u with the basis vectors which are NOT contained in the
// k-blade u and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  lcmpl(u) ^ u  = I_2d = e1^e2  =>  lcmpl(u) = I_2d * rev(u)
// right complement: u ^ rcmpl(u)  = I_2d = e1^e2  =>  rcmpl(u) = rev(u) * I_2d
//
// (get formula on the rhs by multiplying with inv(u) from the right or left respectively)
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined lcmpl(u) = rcmpl(u) = cmpl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// complements have different signs

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> lcmpl(Scalar2d<T> s)
{
    // lcmpl(u) ^ u = e1^e2
    // u = s 1:
    //     lcmpl(u) ^ u = e1^e2 => lcmpl(u) = I_2d * rev(s) = s e1^e2
    return PScalar2d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> lcmpl(Vec2d<T> const& v)
{
    // lcmpl(u) ^ u = e1^e2
    // u = v.x e1 + v.y e2:
    //     lcmpl(u) ^ u = e1^e2 => cmpl(u) = I_2d * rev(v)
    //     lcmpl(u) = e12 * (v.x * e1 + v.y * e2) = -v.x * e2 + v.y * e1
    //                                            =  v.y * e1 - v.x * e2
    //     => lcmpl(u).x =  v.y e1
    //     => lcmpl(u).y = -v.x e2
    return Vec2d<T>(v.y, -v.x);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> lcmpl(PScalar2d<T> ps)
{
    // lcmpl(u) ^ u = e1^e2
    // u = ps e1^e2:
    //     lcmpl(u) ^ u = e1^e2 => lcmpl(u) = I_2d * rev(ps) = ps 1
    return Scalar2d<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> lcmpl(MVec2d_E<T> const& M)
{
    // use the component complements directly
    return MVec2d_E<T>(lcmpl(gr2(M)), lcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> lcmpl(MVec2d<T> const& M)
{
    // use the component complements directly
    return MVec2d<T>(lcmpl(gr2(M)), lcmpl(gr1(M)), lcmpl(gr0(M)));
}


template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> rcmpl(Scalar2d<T> s)
{
    // u ^ rcmpl(u) = e1^e2
    // u = s 1:
    //     u ^ rcmpl(u) = e1^e2 => rcmpl(u) = rev(s) * I_2d = s e1^e2
    return PScalar2d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> rcmpl(Vec2d<T> const& v)
{
    // u ^ rcmpl(u) = e1^e2
    // u = v.x e1 + v.y e2:
    //     u ^ rcmpl(u) = e1^e2 => rcmpl(u) = rev(v) * I_2d
    //     rcmpl(u) = (v.x * e1 + v.y * e2) * e12 =   v.x * e2 - v.y * e1
    //                                            = - v.y * e1 + v.x * e2
    //     => rcmpl(u).x = -v.y e1
    //     => rcmpl(u).y =  v.x e2

    return Vec2d<T>(-v.y, v.x);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> rcmpl(PScalar2d<T> ps)
{
    // u ^ rcmpl(u) = e1^e2
    // u = ps e1^e2:
    //     u ^ rcmpl(u) = e1^e2 => rcmpl(u) = rev(ps) * I_2d = ps 1
    return Scalar2d<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> rcmpl(MVec2d_E<T> const& M)
{
    // use the component complements directly
    return MVec2d_E<T>(rcmpl(gr2(M)), rcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> rcmpl(MVec2d<T> const& M)
{
    // use the component complements directly
    return MVec2d<T>(rcmpl(gr2(M)), rcmpl(gr1(M)), rcmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Lengyel, "PGA illuminated")
// is defined w.r.t. the outer product
////////////////////////////////////////////////////////////////////////////////
//
// if M represents the subspace B as subspace of R^2 then
// dual(M) represents the subspace orthorgonal to B
//
// dual(A) = cmpl(A) in spaces of uneven dimension
//         = rcmpl(A) in spaces of even dimension (right dual)
//
// this dual satisfies (right) dual(A) = rev(A) * I_n
//
// analogy: for the left cmpl/dual the (left) dual satisfies ldual(L) = I_n * rev(L)
//
// -> derived from the defining equation of the left and right complements


// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars, this can only be avoided, when the scalar
// type is uniquely defined for the corresponding algebra
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual(Scalar2d<T> s)
{
    return rcmpl(s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> dual(Vec2d<T> const& v)
{
    return rcmpl(v);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> dual(PScalar2d<T> ps)
{
    return rcmpl(ps);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> dual(MVec2d_E<T> const& M)
{
    return MVec2d_E<T>(rcmpl(gr2(M)), rcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> dual(MVec2d<T> const& M)
{
    return MVec2d<T>(rcmpl(gr2(M)), rcmpl(gr1(M)), rcmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of v1 onto v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto(Vec2d<T> const& v1,
                                                              Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of v1 from v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from(Vec2d<T> const& v1,
                                                             Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec3d<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}

// reflect a vector u on a hyperplane B orthogonal to vector b
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n
// (e.g. a line through the origin in 2d space)
// orthogonal to vector b: the hyperplane is dual to b
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reflect_on_hyp(Vec2d<T> const& u,
                                                                Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-b * u * inv(b));
}

// reflect a vector u another vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reflect_on_vec(Vec2d<T> const& u,
                                                                Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(b * u * inv(b));
}

////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt-Orthogonalization:
////////////////////////////////////////////////////////////////////////////////
// input:  two linear independent vectors u and v in 2d
// output: two orthogonal vectors with the first one being u and the second one a vector
// perpendicular to u in the orientation of v, both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::vector<Vec2d<std::common_type_t<T, U>>>
gs_orthogonal(Vec2d<T> const& u, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec2d<ctype>> basis;
    basis.push_back(u);
    basis.emplace_back(reject_from(v, u));
    return basis;
}

// input:  two linear independent vectors u and v in 2d
// output: two orthonormal vectors with the first one being normalize(u) and
// the second one a normalized vector perpendicular to u in the orientation of v,
// both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::vector<Vec2d<std::common_type_t<T, U>>>
gs_orthonormal(Vec2d<T> const& u, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec2d<ctype>> basis;
    Vec2d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from(v, u_unitized)));
    return basis;
}

} // namespace hd::ga::ega