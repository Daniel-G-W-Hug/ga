#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::max
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <numbers>   // math constants like pi
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "detail/ga_mvec2dp.hpp" // inclusion of multivector imports all component types
#include "detail/ga_pga_2dp_objects.hpp" // Point2dp, Vector2d, Point2d, Line2dp


namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// dot products for 2dp
////////////////////////////////////////////////////////////////////////////////

// return dot-product of two vectors in the modelled space G<2,0,1>
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
// dot(v,v) = |v|^2 = gr0(rev(v)*v)
//
// arguments: Vec2dp<T> or Point2dp<T> : public Vec2dp<T>
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    // definition: dot(v1, v2) = (v1)^T g_12 v2 with the metric g_12
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = 0
    // as diagonal elements of g_12, thus the z-component does not matter
    return v1.x * v2.x + v1.y * v2.y;
}

// return dot product of two bivectors A and B (= a scalar)
// dot(A,B) = gr0(A * B)
//
// arguments: BiVec2dp<T> or Line2dp<T> : public BiVec2dp<T>
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> dot(BiVec2dp<T> const& A, BiVec2dp<U> const& B)
{
    // definition: dot(A, B) = gr0(A*B)
    // -> only the symmetric (i.e. scalar) part remains
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = 0
    // and dot(e23, e23) = 0, dot(e31,e31) = 0, dot(e12,e12) = -1
    // and all other dot(exy,ezw)=0
    // (every index combination containing an index twice becomes zero)
    return -A.z * B.z;
}

////////////////////////////////////////////////////////////////////////////////
// Vec2dp<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared norm of vector
template <typename T> inline T bulk_nrm_sq(Vec2dp<T> const& v)
{
    // |v|^2 = gr0(rev(v)*v)
    return dot(v, v);
}

// return norm of vector
template <typename T> inline T bulk_nrm(Vec2dp<T> const& v)
{
    return std::sqrt(bulk_nrm_sq(v));
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec2dp<T> inv(Vec2dp<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    // using |v|^2 = gr0(rev(v)*v) = dot(v,v)
    //
    T sq_v = dot(v, v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_v < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector dot product too small for inversion" +
                                 std::to_string(sq_v) + "\n");
    }
#endif
    T inv = T(1.0) / sq_v; // inverse of squared norm for a vector
    return Vec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the angle between of two vectors
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    ctype nrm_prod = bulk_nrm(v1) * bulk_nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return a vector unitized to v.z == 1.0
template <typename T> inline Vec2dp<T> unitize(Vec2dp<T> const& v)
{
    T n = v.z;
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for unitization" +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec2dp<T>(v.x * inv, v.y * inv, T(1.0));
}

////////////////////////////////////////////////////////////////////////////////
// BiVec2dp<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of bivector
template <typename T> inline constexpr T bulk_nrm_sq(BiVec2dp<T> const& v)
{
    // |v|^2 = gr0(rev(v)*v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 2-blade => rev(v) = -v
    // using |v|^2 = gr0(rev(v)*v) = gr0(-v*v) = -gr0(v*v) = -dot(v,v)
    return -dot(v, v);
}

// return magnitude of bivector
template <typename T> inline constexpr T bulk_nrm(BiVec2dp<T> const& v)
{
    return std::sqrt(bulk_nrm_sq(v));
}

// return conjugate complex of a bivector
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline constexpr BiVec2dp<T> rev(BiVec2dp<T> const& v)
{
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 2-blade => rev(v) = -v
    // all bivector parts switch sign
    return BiVec2dp<T>(-v.x, -v.y, -v.z);
}

// return the multiplicative inverse of the bivector
template <typename T> inline constexpr BiVec2dp<T> inv(BiVec2dp<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 2-blade => rev(v) = -v
    // using |v|^2 = gr0(rev(v)*v) = gr0(-v*v) = -gr0(v*v) = -dot(v,v)
    // => v^(-1) = (-v)/(-dot(v,v)) = v/dot(v,v)
    T sq_n = -dot(v, v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = -T(1.0) / sq_n; // negative inverse of squared norm for a bivector
    return BiVec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(BiVec2dp<T> const& v1, BiVec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v1) * nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

////////////////////////////////////////////////////////////////////////////////
// Vec2dp<T> & BiVec2dp<T> mixed geometric operations
////////////////////////////////////////////////////////////////////////////////

// return the dot product of a bivector A and a vector b
// dot(A,b) = gr1(A * b)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec2dp<std::common_type_t<T, U>> dot(BiVec2dp<T> const& A, Vec2dp<U> const& b)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(A.z * b.y, -A.z * b.x, A.y * b.x - A.x * b.y);
}

// return the dot product of a vector a and a bivector B
// dot(a,B) = gr1(a * B)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec2dp<std::common_type_t<T, U>> dot(Vec2dp<T> const& a, BiVec2dp<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(-a.y * B.z, a.x * B.z, a.y * B.x - a.x * B.y);
}

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

// return the angle between of a vector and a bivector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec2dp<T> const& v1, BiVec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v1) * nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(
        std::clamp(ctype(nrm(dot(v1, v2))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between of a bivector and a vector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(BiVec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v1) * nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(
        std::clamp(ctype(nrm(dot(v1, v2))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

////////////////////////////////////////////////////////////////////////////////
// wedge products and join operations
////////////////////////////////////////////////////////////////////////////////

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

// wedge product between two points (aka vectors with implicit v.z=1)
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Line2dp<std::common_type_t<T, U>> wdg(Point2d<T> const& p, Point2d<U> const& q)
{
    return Line2dp<std::common_type_t<T, U>>(p.y - q.y, q.x - p.x, p.x * q.y - p.y * q.x);
}

// wedge product between a vector a and a bivector B
// => returns a trivector (i.e. the pseudoscalar for 2dp)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& a, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-a.x * B.x - a.y * B.y - a.z * B.z);
}

// wedge product between a bivector A and a vector b
// wdg(A,b) = gr3(A * b)
// => returns a trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> wdg(BiVec2dp<T> const& A, Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(-A.x * b.x - A.y * b.y - A.z * b.z);
}

// wedge product between a trivector ps and a scalar s
// wdg(ps,s) = s*ps*e321
// => returns a trivector (scalar s multiple of the trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> wdg(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

// wedge product between a scalar s and a trivector ps
// wdg(s,ps) = s*ps*e321
// => returns a trivector (scalar s multiple of the trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> wdg(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2dp<ctype>(ctype(s) * ctype(ps));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    return wdg(v1, v2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Line2dp<std::common_type_t<T, U>> join(Point2d<T> const& p, Point2d<U> const& q)
{
    return wdg(p, q);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> join(Vec2dp<T> const& a, BiVec2dp<U> const& B)
{
    return wdg(a, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2dp<std::common_type_t<T, U>> join(BiVec2dp<T> const& A, Vec2dp<U> const& b)
{
    return wdg(A, b);
}

////////////////////////////////////////////////////////////////////////////////
// regressive wedge products and meet operations
////////////////////////////////////////////////////////////////////////////////

// regressive wedge product of two bivectors
// rwdg(a,b) = !wgd(!a,!b) with ! as the complement operation
// => returns a vector
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& v1, BiVec2dp<U> const& v2)
{
    return Vec2dp<std::common_type_t<T, U>>(
        v1.z * v2.y - v1.y * v2.z, v1.x * v2.z - v1.z * v2.x, v1.y * v2.x - v1.x * v2.y);
}

// regressive wedge product between a vector and a bivector
// rwdg(a,B) = !wgd(!a,!B) with ! as the complement operation
// => returns a scalar
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Scalar2dp<std::common_type_t<T, U>> rwdg(Vec2dp<T> const& a, BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(-a.x * B.x - a.y * B.y - a.z * B.z);
}

// regressive wedge product between a vector and a bivector
// rwdg(A,b) = !wgd(!A,!b) with ! as the complement operation
// => returns a scalar
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Scalar2dp<std::common_type_t<T, U>> rwdg(BiVec2dp<T> const& A, Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(-A.x * b.x - A.y * b.y - A.z * b.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec2dp<std::common_type_t<T, U>> meet(BiVec2dp<T> const& v1, BiVec2dp<U> const& v2)
{
    return rwdg(v1, v2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Scalar2dp<std::common_type_t<T, U>> meet(Vec2dp<T> const& a, BiVec2dp<U> const& B)
{
    return rwdg(a, B);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Scalar2dp<std::common_type_t<T, U>> meet(BiVec2dp<T> const& A, Vec2dp<U> const& b)
{
    return rwdg(A, b);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2dp<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return the reverse
template <typename T> inline constexpr MVec2dp<T> rev(MVec2dp<T> const& v)
{
    // only bivector and trivector parts switch signs
    return MVec2dp<T>(v.c0, v.c1, v.c2, v.c3, -v.c4, -v.c5, -v.c6, -v.c7);
}

// return the Clifford conjugate
template <typename T> inline constexpr MVec2dp<T> conj(MVec2dp<T> const& v)
{
    // only vector and bivector parts switch signs
    return MVec2dp<T>(v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, v.c7);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2dp_E<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return conjugate complex of quaternion (MVec2dp_E<T>),
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline constexpr MVec2dp_E<T> rev(MVec2dp_E<T> const& v)
{
    // only the bivector part switches sign
    return MVec2dp_E<T>(v.c0, -v.c1, -v.c2, -v.c3);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2dp_U<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return the reverse
template <typename T> inline constexpr MVec2dp_U<T> rev(MVec2dp_U<T> const& v)
{
    // only the trivector part switches signs
    return MVec2dp_U<T>(v.c0, v.c1, v.c2, -v.c3);
}

////////////////////////////////////////////////////////////////////////////////
// PScalar2dp<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return reverse of a 2d bivector
template <typename T> inline PScalar2dp<T> rev(PScalar2dp<T> A)
{
    // the 2dp bivector switches sign on reversion
    return PScalar2dp<T>(-T(A));
}

////////////////////////////////////////////////////////////////////////////////
// MVec2dp<T> geometric products
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
    // => due to the degenerate algebra 16 terms are not present compared to G<3,0,0>
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
    ctype c7 = -(-A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
                 A.c5 * B.c2 + A.c6 * B.c3 - A.c7 * B.c0);
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
    ctype c7 = -(-A.c0 * B.c7 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3);
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
    ctype c7 = -(A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c7 * B.c0);
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
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
    ctype c6 = A.c0 * B.c2 - A.c1 * B.c1 - A.c2 * B.c7 - A.c3 * B.c3;
    ctype c7 = -(A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 - A.c3 * B.c0);
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
    ctype c6 = A.c1 * B.c1 - A.c2 * B.c0 - A.c3 * B.c3 - A.c7 * B.c2;
    ctype c7 = -(-A.c0 * B.c3 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2);
    return MVec2dp<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
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
        -PScalar2dp<ctype>(-A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2));
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
        -PScalar2dp<ctype>(A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c3 * B.c0));
}

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(-A.c3 * B.z),
                            BiVec2dp<ctype>(A.c0 * B.x - A.c2 * B.z + A.c3 * B.y,
                                            A.c0 * B.y + A.c1 * B.z - A.c3 * B.x,
                                            A.c0 * B.z));
}

// geometric product A * B of a bivector A with an even grade multivector B
// from the right
// bivector * even grade multivector => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& A,
                                                               MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(-A.z * B.c3),
                            BiVec2dp<ctype>(A.x * B.c0 - A.y * B.c3 + A.z * B.c2,
                                            A.x * B.c3 + A.y * B.c0 - A.z * B.c1,
                                            A.z * B.c0));
}

// geometric product A * b of an even grade multivector A with a vector b
// from the right
// even grade multivector * vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                               Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(A.c0 * b.x + A.c3 * b.y,
                                          A.c0 * b.y - A.c3 * b.x,
                                          A.c0 * b.z - A.c1 * b.y + A.c2 * b.x),
                            -PScalar2dp<ctype>(A.c1 * b.x + A.c2 * b.y + A.c3 * b.z));
}

// geometric product a * B of a vector a and an even grade multivector B from the left
// vector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& a,
                                                               MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(a.x * B.c0 - a.y * B.c3,
                                          a.x * B.c3 + a.y * B.c0,
                                          -a.x * B.c2 + a.y * B.c1 + a.z * B.c0),
                            -PScalar2dp<ctype>(a.x * B.c1 + a.y * B.c2 + a.z * B.c3));
}

// geometric product A * B of an uneven grade multivector A with a bivector B
// uneven grade multivector * bivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_U<T> const& A,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(-A.c1 * B.z, A.c0 * B.z, -A.c0 * B.y + A.c1 * B.x + A.c3 * B.z),
        -PScalar2dp<ctype>(A.c0 * B.x + A.c1 * B.y + A.c2 * B.z));
}

// geometric product A * B of a bivector A with an uneven grade multivector B
// bivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& A,
                                                               MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(
        Vec2dp<ctype>(A.z * B.c1, -A.z * B.c0, -A.x * B.c1 + A.y * B.c0 + A.z * B.c3),
        -PScalar2dp<ctype>(A.x * B.c0 + A.y * B.c1 + A.z * B.c2));
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
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& a,
                                                               BiVec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(dot(a, b)), cmt(a, b));
}

// geometric product A*b for a bivector and a vector
// A * b = dot(A,b) + wdg(A,b) = gr1(A * b) + gr3(A * b)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => bivector * vector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& A,
                                                               Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(dot(A, b), wdg(A, b));
}

// geometric product a * B for a vector a and a bivector B
// a * B = dot(a,B) + wdg(a,B) = gr1(a * B) + gr3(a * B)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => vector * bivector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(Vec2dp<T> const& a,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(dot(a, B), wdg(a, B));
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
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(dot(a, b)), wdg(a, b));
}

// geometric product A * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the multivector B
// trivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> A,
                                                             MVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * MVec2dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), -B.c6, B.c1,
                                      B.c2, ctype(0.0), -B.c0);
}

// geometric product A * B of a multivector A with the trivector B from the right
// multivector * trivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp<std::common_type_t<T, U>> operator*(MVec2dp<T> const& A,
                                                             PScalar2dp<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), -A.c6, A.c1, A.c2,
                          ctype(0.0), -A.c0) *
           -ctype(B);
}

// geometric product A * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the even grade multivector B
// trivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(PScalar2dp<T> A,
                                                               MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * MVec2dp_U<ctype>(Vec2dp<ctype>(ctype(0.0), ctype(0.0), -B.c3),
                                        PScalar2dp<ctype>(-B.c0));
}

// geometric product A * B of an even multivector A multiplied from the right
// by the trivector B (=2dp pseudoscalar)
// even grade multivector * trivector => uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<U> const& A,
                                                               PScalar2dp<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(Vec2dp<ctype>(ctype(0.0), ctype(0.0), -A.c3),
                            PScalar2dp<ctype>(-A.c0)) *
           -ctype(B);
}

// geometric product A * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the uneven grade multivector B
// trivector * uneven grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(PScalar2dp<T> A,
                                                               MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * MVec2dp_E<ctype>(Scalar2dp<ctype>(ctype(0.0)),
                                        BiVec2dp<ctype>(B.c0, B.c1, ctype(0.0)));
}

// geometric product A * B of an uneven grade multivector A multiplied from the right
// by the trivector B (=2dp pseudoscalar)
// uneven grade multivector * trivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<U> const& A,
                                                               PScalar2dp<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(Scalar2dp<ctype>(ctype(0.0)),
                            BiVec2dp<ctype>(A.c0, A.c1, ctype(0.0))) *
           -ctype(B);
}

// geometric product A * B of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> A,
                                                            BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * Vec2dp<ctype>(ctype(0.0), ctype(0.0), -B.z);
}

// geometric product A * B of bivector A multiplied by a trivector B from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> operator*(BiVec2dp<T> const& A,
                                                            PScalar2dp<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(ctype(0.0), ctype(0.0), -A.z) * -ctype(B);
}

// geometric product A * b of a trivector A (=2dp pseudoscalar) multiplied from the left
// to the vector b
// trivector * vector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator*(PScalar2dp<T> A,
                                                              Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * BiVec2dp<ctype>(b.x, b.y, ctype(0.0));
}

// geometric product a * B of a vector a multiplied with a trivector B from the right
// vector * trivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> operator*(Vec2dp<T> const& a,
                                                              PScalar2dp<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(a.x, a.y, ctype(0.0)) * -ctype(B);
}

// geometric product A * B of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>>
operator*([[maybe_unused]] PScalar2dp<T> A, [[maybe_unused]] PScalar2dp<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(0.0); // trivectors in 3d pga contain e3^2 = 0
}

// geometric product A * B of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2dp<std::common_type_t<T, U>> operator*(Scalar2dp<T> A,
                                                               Scalar2dp<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2dp<ctype>(ctype(A) * ctype(B));
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
// left complement:  l_complement(u) ^ u  = I_2dp = e3^e2^e1
// right complement: u ^ r_complement(u)  = I_2dp = e3^e2^e1
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> complement(Scalar2dp<T> s)
{
    // u ^ complement(u) = e3^e2^e1
    // u = 1:
    //     1 ^ complement(u) = e3^e2^e1 => complement(u) = e3^e2^e1
    return PScalar2dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline BiVec2dp<T> complement(Vec2dp<T> const& v)
{
    // u ^ compl(u) = e3^e2^e1
    // u = v.x e1 + v.y e2 + v.z e3:
    //     u ^ complement(u) = e3^e2^e1
    //     u = e1 => complement(u) = -e23
    //     u = e2 => complement(u) = -e31
    //     u = e3 => complement(u) = -e12
    return BiVec2dp<T>(-v.x, -v.y, -v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline Vec2dp<T> complement(BiVec2dp<T> const& b)
{
    // u ^ compl(u) = e3^e2^e1
    // u = b.x e23 + b.y e31 + b.z e12:
    //     u ^ complement(u) = e3^e2^e1
    //     u = e23 => complement(u) = -e1
    //     u = e31 => complement(u) = -e2
    //     u = e12 => complement(u) = -e3
    return Vec2dp<T>(-b.x, -b.y, -b.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> complement(PScalar2dp<T> ps)
{
    // u ^ compl(u) = e3^e2^e1
    // u = e3^e2^e1:
    //     e3^e2^e1 ^ complement(u) = e3^e2^e1 => complement(u) = 1
    return Scalar2dp<T>(T(ps));
}

////////////////////////////////////////////////////////////////////////////////
// 2dp rotation operations
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
// exponential function with bivector as argument for setup of quaternions
// as geometric multivector with a scalar and a bivector part
// MVec2dp_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
//
// quaternion: q = a + b I with I being the bivector in brackets above
//             representing a plane in the algebra G^3
//
//             a rotation in 3D is represented by the plane and the
//             size of the rotation, the later is given by the angle
//             theta, which is the magnitude of the bivector
//
// Inputs:
//         - an arbitray bivector representing the oriented plane of rotation
//           (does not need to be normalized)
//         - a rotation angle
// Output:
//         - a rotor representing the rotation
//
// HINT:     For a rotation around an axis n (n = normalize(Vec3d<T>))
//           use the bivector B = n*I_3d  => B = Vec3d<T> * PScalar3d<T>
//////////////////////////////////////////////////////////////////////////////////////////
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
inline constexpr MVec2dp<std::common_type_t<T, U>> rotate(MVec2dp<T> const& v,
                                                          MVec2dp_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(rotor * v * rev(rotor));
}

////////////////////////////////////////////////////////////////////////////////
// Vec2dp<T> and BiVec2dp<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v1,
                                                               Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return dot(v1, v2) * Vec2dp<ctype>(inv(v2));
}

// projection of v1 onto v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>>
project_onto_normalized(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    return dot(v1, v2) * v2; // v2 is already its own reverse
}

// projection of a vector v1 onto a bivector v2
// v_parallel = dot(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v1,
                                                               BiVec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    Vec2dp<ctype> a = dot(v1, v2);
    BiVec2dp<ctype> Bi = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(a * Bi) = dot(a,Bi)
    return Vec2dp<ctype>(dot(a, Bi));
}

// projection of a vector v1 onto a normalized bivector v2
// u_parallel = gr1(dot(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>>
project_onto_normalized(Vec2dp<T> const& v1, BiVec2dp<U> const& v2)
{
    // requires v2 to be normalized

    using ctype = std::common_type_t<T, U>;
    Vec2dp<ctype> a = dot(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec2dp<ctype> Bi = -v2;
    // use the formular equivalent to the geometric product to save computational cost
    // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(a * Bi) = dot(a,Bi)
    return Vec2dp<ctype>(dot(a, Bi));
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v1,
                                                              Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    BiVec2dp<ctype> B = wdg(v1, v2);
    Vec2dp<ctype> v2_inv = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    return Vec2dp<ctype>(dot(B, v2_inv));
}

// rejection of vector v1 from a normalized vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>>
reject_from_normalized(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    // requires v2 to be normalized
    using ctype = std::common_type_t<T, U>;
    BiVec2dp<ctype> B = wdg(v1, v2);
    Vec2dp<ctype> v2_inv = v2; // v2 is its own inverse, if normalized
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    return Vec2dp<ctype>(dot(B, v2_inv));
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

// rejection of vector v1 from a normalized bivector v2
// u_perp = wdg(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>>
reject_from_normalized(Vec2dp<T> const& v1, BiVec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar2dp<ctype> a = wdg(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec2dp<ctype> B = -v2;
    // trivector * bivector = vector (derived from full geometric product to save
    // costs)
    return Vec2dp<ctype>(-a * B.x, -a * B.y, -a * B.z);
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