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

#include "detail/ga_mvec2dp.hpp" // inclusion of multivector imports all component types
#include "detail/ga_pga_2dp_objects.hpp" // Point2dp, Vector2d, Point2d, Line2d


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
inline constexpr MVec2dp_E<T> gr_inv(MVec2dp_E<T> const& E)
{
    // grade 0 and 2: no sign change
    return E;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> gr_inv(MVec2dp_U<T> const& U)
{
    // grade 1 and 3: sign reversal
    return -U;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> gr_inv(MVec2dp<T> const& M)
{
    // grade 0 and 2: no sign change
    // grade 1 and 3: sign reversal
    return MVec2dp<T>(M.c0, -M.c1, -M.c2, -M.c3, M.c4, M.c5, M.c6, -M.c7);
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
inline constexpr MVec2dp_E<T> rev(MVec2dp_E<T> const& E)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec2dp_E<T>(E.c0, -E.c1, -E.c2, -E.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> rev(MVec2dp_U<T> const& U)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec2dp_U<T>(U.c0, U.c1, U.c2, -U.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> rev(MVec2dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    return MVec2dp<T>(M.c0, M.c1, M.c2, M.c3, -M.c4, -M.c5, -M.c6, -M.c7);
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
inline constexpr MVec2dp_E<T> rrev(MVec2dp_E<T> const& E)
{
    // grade 0: sign reversal
    // grade 2: no sign change
    return MVec2dp_E<T>(-E.c0, E.c1, E.c2, E.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> rrev(MVec2dp_U<T> const& U)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dp_U<T>(-U.c0, -U.c1, -U.c2, U.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> rrev(MVec2dp<T> const& M)
{
    // grade 0: sign reversal
    // grade 1: sign reversal
    // grade 2: no sign change
    // grade 3: no sign change
    return MVec2dp<T>(-M.c0, -M.c1, -M.c2, -M.c3, M.c4, M.c5, M.c6, M.c7);
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
inline constexpr MVec2dp_E<T> conj(MVec2dp_E<T> const& E)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec2dp_E<T>(E.c0, -E.c1, -E.c2, -E.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> conj(MVec2dp_U<T> const& U)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dp_U<T>(-U.c0, -U.c1, -U.c2, U.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> conj(MVec2dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    return MVec2dp<T>(M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, M.c7);
}


////////////////////////////////////////////////////////////////////////////////
// dot products for 2dp (= defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = 0
//
////////////////////////////////////////////////////////////////////////////////

// scalar product dot(a,b) (nrm_sq(a,b) = dot(a, rev(b)))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> dot(MVec2dp<T> const& A,
                                                         MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c6 * B.c6);
}

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
    return Scalar2dp<ctype>(-B1.z * B2.z);
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
    return PScalar2dp<ctype>(-v1.z * v2.z);
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
inline PScalar2dp<std::common_type_t<T, U>> wdg(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(ps) * ctype(s));
}

// wedge product between a scalar s and a trivector ps
// wdg(s,ps) = s*ps*e123
// => returns a trivector (scalar s multiple of the trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, PScalar2dp<U> ps)
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
inline PScalar2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& B, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-B.x * v.x - B.y * v.y - B.z * v.z);
}

// wedge product between a vector a and a bivector B
// => returns a trivector (i.e. the pseudoscalar for 2dp)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v, BiVec2dp<U> const& B)
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
inline BiVec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
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
    return weight_expansion(point2dp{p}, l);
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
// L. Dorst: The contraction A << B of an a-blade A onto a b-blade B is a sub-blade
// of B of grade b-a which is perpendicular to A, and linear in both arguments
////////////////////////////////////////////////////////////////////////////////

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
    return Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.z) * ctype(ps);
}

// left contraction -  pseudoscalar contracted onto vector
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
    return BiVec2dp<ctype>(-v.x, -v.y, ctype(0.0)) * ctype(ps);
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
    return Scalar2dp<ctype>(-B1.z * B2.z);
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
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator<<(Vec2dp<T> const& v,
                                                             BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(-v.y * B.z, v.x * B.z, -v.x * B.y + v.y * B.x);
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
////////////////////////////////////////////////////////////////////////////////

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
inline constexpr Vec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<U> ps,
                                                             BiVec2dp<T> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.z);
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
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<U> ps,
                                                               Vec2dp<T> const& v)
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
inline constexpr PScalar2dp<std::common_type_t<T, U>> operator>>(PScalar2dp<U> ps,
                                                                 Scalar2dp<T> s)
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
    return Scalar2dp<ctype>(-B1.z * B2.z);
}

// right contraction of bivector B by vector v ("B contracted by v")
// => returns a vector
// (originally implemented as dot-product)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<T> const& B,
                                                   Vec2dp<U> const& v)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(B.z * v.y, -B.z * v.x, -B.x * v.y + B.y * v.x);
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
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator>>(BiVec2dp<U> const& B,
                                                               Scalar2dp<T> s)
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
inline constexpr Vec2dp<std::common_type_t<T, U>> operator>>(Vec2dp<U> const& v,
                                                             Scalar2dp<T> s)
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
// alternative multivector products (in use instead of contractions)
////////////////////////////////////////////////////////////////////////////////

// return commutator product cmt(A,B) of two bivectors A and B (= a bivector)
// cmt(A,B) = 0.5*(AB-BA) = gr2(A * B)
// the commutator product is antisymmetric, i.e. it is zero when a bivector is
// multiplied by itself, i.e. in that case only the dot product remains
// as the symmetric part
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec2dp<std::common_type_t<T, U>> cmt(BiVec2dp<T> const& A, BiVec2dp<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(A.z * B.y - A.y * B.z, A.x * B.z - A.z * B.x, ctype(0.0));
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

// multivector * uneven multivector => multivector
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

// uneven multivector * multivector => multivector
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

// geometric product A * B for two multivectors from the uneven subalgebra (3d case)
// uneven grade multivector * uneven grade multivector = even grade multivector
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
// with a multivector B of the uneven subalgebra
// even grade multivector * uneven grade multivector => uneven grade multivector
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

// geometric product A * B of a multivector A from the uneven subalgebra (3d case)
// with a multivector B of the even subalgebra
// uneven grade multivector * even grade multivector => uneven grade multivector
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
// even grade multivector * trivector => uneven multivector
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
// trivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(PScalar2dp<T> ps,
                                                               MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec2dp_U<ctype>(Vec2dp<ctype>(ctype(0.0), ctype(0.0), B.c3),
                                        PScalar2dp<ctype>(B.c0));
}

// geometric product A * B of an uneven grade multivector A multiplied from the right
// by the trivector B (=2dp pseudoscalar)
// uneven grade multivector * trivector => even grade multivector
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
// to the uneven grade multivector B
// trivector * uneven grade multivector => even grade multivector
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

// geometric product M * B of an uneven grade multivector M with a bivector B
// uneven grade multivector * bivector => uneven grade multivector
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

// geometric product B * M of a bivector A with an uneven grade multivector B
// bivector * uneven grade multivector => uneven grade multivector
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
// even grade multivector * vector => uneven grade multivector
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
// vector * even grade multivector => uneven grade multivector
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

// geometric product A * B between two bivectors
// A * B = gr0(A * B) + gr2(A * B)
//
// the full geometric bivector product only exists in >= 4d spaces:
// A * B = A * B + cmt(A,B) + wdg(A,B) = gr0(A * B) + gr2(A * B) + gr4(A * B)
// In 3D we don't have gr4(A * B) and thus only the terms up to grade 3 remain.
// The bivector product AxB = cmt(A,B) = 0.5*(ab-ba)is called the commutator product.
//
// A * B = dot(A,B) + cmt(A,B) + wdg(A,B)  (in 4d and higher dimensional spaces)
// A * B = dot(A,B) + cmt(A,B)             (in 3d)
//
// => bivector * bivector = scalar + bivector = even grade multivector (in 3d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B1,
                                                               BiVec2dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(dot(B1, B2)), cmt(B1, B2));
}

// geometric product B * v for a bivector B and a vector v
// old: B * v = dot(B,v) + wdg(B,v) = gr1(B * v) + gr3(B * v)
// new: B * v = (B >> v) + wdg(B,v) = gr1(B * v) + gr3(B * v)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => bivector * vector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& B,
                                                               Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    // return MVec2dp_U<ctype>(dot(B, v), wdg(B, v));
    return MVec2dp_U<ctype>((B >> v), wdg(B, v));
}

// geometric product v * B for a vector v and a bivector B
// old: v * B = dot(v,B) + wdg(v,B) = gr1(v * B) + gr3(v * B)
// new: v * B = (v << B) + wdg(v,B) = gr1(v * B) + gr3(v * B)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => vector * bivector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& v,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // return MVec2dp_U<ctype>(dot(v, B), wdg(v, B));
    return MVec2dp_U<ctype>((v << B), wdg(v, B));
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

// return the multiplicative inverse of the vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> inv(Vec2dp<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    // using |v|^2 = gr0(rev(v)*v) = dot(v,v)
    //
    T sq_v = dot(v, v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_v < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector dot product too small for inversion " +
                                 std::to_string(sq_v) + "\n");
    }
#endif
    T inv = T(1.0) / sq_v; // inverse of squared norm for a vector
    return Vec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> inv(BiVec2dp<T> const& B)
{
    // B^(-1) = rev(B)/|B|^2
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    // using |B|^2 = gr0(rev(B)*B) = gr0(-B*B) = -gr0(B*B) = -dot(B,B)
    // => B^(-1) = (-B)/(-dot(B,B)) = B/dot(B,B)
    T sq_n = -dot(B, B);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = -T(1.0) / sq_n; // negative inverse of squared norm for a bivector
    return BiVec2dp<T>(B.x * inv, B.y * inv, B.z * inv);
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
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>>
motor2dp_from_ln(BiVec2dp<T> const& B1, BiVec2dp<U> const& B2)
{
    // take lines as input an return a motor
    // 1st apply reflection across line B1, then across B2 to get a motor that rotates
    // (or translates) around the intersection point of lines B1 and B2
    //
    // for use of motor m either directly on object u (inefficient):
    //     auto v_moved = gr1(rgpr(rgpr(R, v), rrev(R))));
    // or
    //     auto B_moved = gr2(rgpr(rgpr(R, B), rrev(R))));
    // or
    //                                   // optimized for reduced effort
    //     auto v_moved = move2dp(v,R);  // moves v according to the motor R
    //     auto B_moved = move2dp(B,R);  // moves B according to the motor R
    //
    return rgpr(B2, B1); // based on the regressive geometric product
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp_orig(Vec2dp<T> const& v,
                                                               MVec2dp_U<U> const& R)
{
    // moves v (a vector representing a projective point) according to the motor R
    return gr1(rgpr(rgpr(R, v), rrev(R))); // the inefficient original
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp_orig(BiVec2dp<T> const& B,
                                                                 MVec2dp_U<U> const& R)
{
    // moves B (a bivector representing a line) according to the motor R
    return gr2(rgpr(rgpr(R, B), rrev(R))); // the inefficient original
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp(Vec2dp<T> const& v,
                                                          MVec2dp_U<U> const& R)
{
    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    //
    // (could potentially be further optimized by exporting the matrix-representation if
    // many transformations should be calculated using the same rotor as v' = matrix * v)
    using ctype = std::common_type_t<T, U>;
    auto k02 = R.c0 * R.c2;
    auto k03 = R.c0 * R.c3;
    auto k12 = R.c1 * R.c2;
    auto k13 = R.c1 * R.c3;
    auto k22 = R.c2 * R.c2;
    auto k23 = R.c2 * R.c3;
    auto k33 = R.c3 * R.c3;
    return Vec2dp<ctype>((k22 - k33) * v.x + (2.0 * k23) * v.y + 2.0 * (-k02 - k13) * v.z,
                         -2.0 * k23 * v.x + (k22 - k33) * v.y + 2.0 * (k03 - k12) * v.z,
                         (-k22 - k33) * v.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp(BiVec2dp<T> const& B,
                                                            MVec2dp_U<U> const& R)
{
    // moves B (a bivector representing a line) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    //
    // (could potentially be further optimized by exporting the matrix-representation if
    // many transformations should be calculated using the same rotor as B' = matrix * B)
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

////////////////////////////////////////////////////////////////////////////////
// 2dp rotation operations
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> exp(BiVec2dp<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(std::cos(theta)),
                            normalize(I) * std::sin(theta));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//       1.) an arbitray bivector representing the oriented plane of rotation
//           (does not need to be normalized, defines what is a posive rotation angle)
//       2.) a rotation angle in that plane
// Output:
//           a rotor representing the requested rotation,
//           for applying the sandwich product as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// for a rotation about an axis n (n = normalized vector) choose the ansatz n*B = I_3d
// and multiply both sides with n from the left (remember n*n = |n|^2 = 1)
//
// => choose: B = n*I_3d
//
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> rotor(BiVec2dp<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(std::cos(half_angle)),
                            normalize(I) * std::sin(half_angle));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> rotate(Vec2dp<T> const& v,
                                                         MVec2dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec2dp_E<ctype> rr = rev(rotor);
    // MVec2dp_U<ctype> tmp = rotor * v;
    // MVec2dp_U<ctype> res = tmp * rr;

    // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the pseudoscalar part
    // which will be zero anyway
    return Vec2dp<ctype>(gr1<ctype>(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> rotate_opt(Vec2dp<T> const& v,
                                                             MVec2dp_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = R.c0 * v.x + R.c3 * v.y;
    ctype k2 = R.c0 * v.y - R.c3 * v.x;
    ctype k3 = R.c0 * v.z - R.c1 * v.y + R.c2 * v.x;
    ctype k4 = R.c1 * v.x + R.c2 * v.y + R.c3 * v.z;
    return Vec2dp<ctype>(k1 * R.c0 + k2 * R.c3, -k1 * R.c3 + k2 * R.c0,
                         k1 * R.c2 - k2 * R.c1 + k3 * R.c0 + k4 * R.c3);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> rotate(BiVec2dp<T> const& v,
                                                           MVec2dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec2dp_E<ctype> rr = rev(rotor);
    // MVec2dp_E<ctype> tmp = rotor * v;
    // MVec2dp_E<ctype> res = tmp * rr;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway
    return BiVec2dp<ctype>(gr2<ctype>(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> rotate_opt(BiVec2dp<T> const& B,
                                                               MVec2dp_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = -R.c3 * B.z;
    ctype k2 = R.c0 * B.x - R.c2 * B.z + R.c3 * B.y;
    ctype k3 = R.c0 * B.y + R.c1 * B.z - R.c3 * B.x;
    ctype k4 = R.c0 * B.z;
    return BiVec2dp<ctype>(-k1 * R.c1 + k2 * R.c0 + k3 * R.c3 - k4 * R.c2,
                           -k1 * R.c2 - k2 * R.c3 + k3 * R.c0 + k4 * R.c1,
                           -k1 * R.c3 + k4 * R.c0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> rotate(MVec2dp<T> const& v,
                                                          MVec2dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(rotor * v * rev(rotor));
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
inline BiVec2dp<T> cmpl(Vec2dp<T> const& v)
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
inline Vec2dp<T> cmpl(BiVec2dp<T> const& B)
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
// projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
// returns component of v1 parallel to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v1,
                                                               Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * Vec2dp<ctype>(inv(v2));
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2)) = v1 - project_onto(v1, v2)
// returns component of v1 perpendicular to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v1,
                                                              Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec2dp<ctype>(gr1(wdg(v1, v2) * inv(v2)));

    return Vec2dp<ctype>(v1 - project_onto(v1, v2));
}


// projection of a vector v1 onto a bivector v2
// v_parallel = dot(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // // Vec2dp<ctype> a = dot(v, B);
    // Vec2dp<ctype> a = (v << B);
    // BiVec2dp<ctype> Bi = inv(B);
    // // use the formular equivalent to the geometric product to save computational cost
    // // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // // v_parallel = gr1(a * Bi) = dot(a,Bi)
    // // return Vec2dp<ctype>(dot(a, Bi));
    // return Vec2dp<ctype>((a << Bi));

    return Vec2dp<ctype>((v << inv(B)) << B);
    // return Vec2dp<ctype>(gr1((v << inv(B)) * B));
}


// rejection of vector v1 from a bivector v2
// u_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v1,
                                                              BiVec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar2dp<ctype> A = wdg(v1, v2);
    BiVec2dp<ctype> B = inv(v2);
    // trivector * bivector = vector
    return A * B;
}

// reflect a vector u on a hyperplane B orthogonal to vector b
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n (a line in 2d space)
// orthogonal to vector b: the hyperplane is dual to b (i.e. a one dimensional subspace)
//
// HINT: choose b * B = I_2dp  =>  B = b * I_2dp  (for normalized b)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reflect_on_hyp(Vec2dp<T> const& u,
                                                                 Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(gr1(-b * u * inv(b)));
}

// reflect a vector u in an arbitrary bivector, i.e. a plane
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reflect_on(Vec2dp<T> const& u,
                                                             BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(gr1(-B * u * inv(B)));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling planes)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> reflect_on(BiVec2dp<T> const& UB,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(B * UB * inv(B)));
}

// reflect a vector u another vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reflect_on_vec(Vec2dp<T> const& u,
                                                                 Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(gr1(b * u * inv(b)));
}

} // namespace hd::ga::pga