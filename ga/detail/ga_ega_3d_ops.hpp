#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::max, std::clamp
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <numbers>   // math constants like pi
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "ga_mvec3d.hpp" // inclusion of multivector imports all component types

namespace hd::ga::ega {

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return the dot product of two vectors in G<3,0,0>
// coordinate free definition: dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    // definition: dot(v1, v2) = (v1)^T g_12 v2 using the metric g_12
    // definition: dot(A, B) = gr0(A*B) using the geometric product
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = +1
    // as diagonal elements of g_12
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// return squared magnitude of vector
template <typename T> inline T nrm_sq(Vec3d<T> const& v) { return magn_sq(v); }

// return magnitude of vector
template <typename T> inline T nrm(Vec3d<T> const& v) { return magn(v); }

// return the multiplicative inverse of the vector
template <typename T> inline Vec3d<T> inv(Vec3d<T> const& v)
{
    // v^(-1) = v/dot(v,v)
    T sq_v = dot(v, v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_v < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector dot produc too small for inversion" +
                                 std::to_string(sq_v) + "\n");
    }
#endif
    T inv = T(1.0) / sq_v; // inverse of squared norm for a vector
    return Vec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U> const& v2)
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


// unsuccessful try to extend angle range to -pi ... pi,
// because orientation is not defined uniquely:
//
// // return the angle between of two vectors
// // range of angle: -pi <= angle <= pi
// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U>
// const& v2)
// {
//     using ctype = std::common_type_t<T, U>;
//     using std::numbers::pi;

// ctype nrm_prod = nrm(v1) * nrm(v2);
// #if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
// if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
//     throw std::runtime_error("vector norm product too small for calculation "
//                              "of angle" +
//                              std::to_string(nrm_prod) + "\n");
// }
// #endif

// auto cos_angle = std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
// auto sin_angle = std::clamp(ctype(nrm(wdg(v1, v2))) / nrm_prod, ctype(-1.0),
// ctype(1.0));
//     // wdg() does contain magnitude, but no unique value of orientation
//     // so we chose one arbitrarily => but would deliver only pos. angles!

//     fmt::println("   c = {: .4f}, s = {: .4f}, wdg = {: .4f}, nrm_wdg = {: .4f}",
//                  cos_angle, sin_angle, wdg(v1, v2), nrm(wdg(v1, v2)));

//     if (cos_angle >= 0.0) {
//         // quadrant I or IV
//         return std::asin(sin_angle);
//     }
//     else if (cos_angle < 0.0 && sin_angle >= 0.0) {
//         // quadrant II
//         return pi - std::asin(sin_angle);
//     }
//     else {
//         // cos_angle < 0.0 && sin_angle < 0.0)
//         // quadrant III
//         return -pi - std::asin(sin_angle);
//     }
// }

// cross-product between two vectors (only defined in 3d)
//  => returns a 3d vector
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> cross(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}
//     double cross product identity:
//     a x (b x c) = -dot(a, wdg(b,c)) = -a.(b^c)


////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return dot product of two bivectors A and B (= a scalar)
// dot(A,B) = gr0(A * B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> dot(BiVec3d<T> const& A, BiVec3d<U> const& B)
{
    // definition: dot(A, B) = gr0(A*B)
    // -> only the symmetric (i.e. scalar) part remains
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = +1
    // and dot(e23, e23) = -1, dot(e31,e31) = -1, dot(e12,e12) = -1
    // and all other dot(exy,ezw)=0
    return -A.x * B.x - A.y * B.y - A.z * B.z;
}

// return squared magnitude of bivector
template <typename T> inline constexpr T nrm_sq(BiVec3d<T> const& v)
{
    // |v|^2 = gr0(rev(v)*v)
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

// return magnitude of bivector
template <typename T> inline constexpr T nrm(BiVec3d<T> const& v)
{
    return std::sqrt(nrm_sq(v));
}

// return conjugate complex of a bivector
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline constexpr BiVec3d<T> rev(BiVec3d<T> const& v)
{
    // definition: rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    // all bivector parts switch sign
    return BiVec3d<T>(-v.x, -v.y, -v.z);
}

// return the multiplicative inverse of the bivector
template <typename T> inline constexpr BiVec3d<T> inv(BiVec3d<T> const& v)
{
    T sq_n = nrm_sq(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = -T(1.0) / sq_n; // negative inverse of squared norm for a bivector
    return BiVec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(BiVec3d<T> const& v1, BiVec3d<U> const& v2)
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
// Vec3d<T> & BiVec3d<T> mixed geometric operations
////////////////////////////////////////////////////////////////////////////////

// return the dot product of a bivector A and a vector b
// dot(A,b) = gr1(A * b)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> dot(BiVec3d<T> const& A, Vec3d<U> const& b)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(A.z * b.y - A.y * b.z, A.x * b.z - A.z * b.x,
                        A.y * b.x - A.x * b.y);
}

// return the dot product of a vector a and a bivector B
// dot(a,B) = gr1(a * B)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> dot(Vec3d<T> const& a, BiVec3d<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(a.z * B.y - a.y * B.z, a.x * B.z - a.z * B.x,
                        a.y * B.x - a.x * B.y);
}

// return commutator product cmt(A,B) of two bivectors A and B (= a bivector)
// cmt(A,B) = 0.5*(AB-BA) = gr2(A * B)
// the commutator product is antisymmetric, i.e. it is zero when a bivector is
// multiplied by itself, i.e. in that case only the dot product remains
// as the symmetric part
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3d<std::common_type_t<T, U>> cmt(BiVec3d<T> const& A, BiVec3d<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(A.z * B.y - A.y * B.z, A.x * B.z - A.z * B.x,
                          A.y * B.x - A.x * B.y);
}

// return the angle between of a vector and a bivector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, BiVec3d<U> const& v2)
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
inline std::common_type_t<T, U> angle(BiVec3d<T> const& v1, Vec3d<U> const& v2)
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

// wedge product between two vectors (returns a bivector in 3d)
// coordinate-free definition: wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                          v1.x * v2.y - v1.y * v2.x);
}

// wedge product between a vector a and a bivector B
// wdg(a,B) = gr3(a * B)
// => returns a trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& a, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(a.x * B.x + a.y * B.y + a.z * B.z);
}

// wedge product between a bivector A and a vector b
// wdg(A,b) = gr3(A * b)
// => returns a trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& A, Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(A.x * b.x + A.y * b.y + A.z * b.z);
}

////////////////////////////////////////////////////////////////////////////////
// PScalar3d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return the reverse of a 3d trivector
template <typename T> inline PScalar3d<T> rev(PScalar3d<T> A)
{
    // the 3d trivector switches sign on reversion
    return PScalar3d<T>(-T(A));
}

// return inverse of the pseudoscalar (A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// k is the dimension of the space of the pseudoscalar formed by k orthogonal vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> inv(PScalar3d<T> ps)
{
    return -PScalar3d<T>(ps) / nrm_sq(ps);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude
// |M|^2 = M rev(M) = (M.c0)^2 + (M.c1)^2 + (M.c2)^2 + (M.c3)^3
//                  + (M.c4)^2 + (M.c5)^2 + (M.c6)^2 + (M.c7)^3
template <typename T> inline T nrm_sq(MVec3d<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3 + v.c4 * v.c4 +
           v.c5 * v.c5 + v.c6 * v.c6 + v.c7 * v.c7;
}

// return magnitude of complex number
template <typename T> inline T nrm(MVec3d<T> const& v) { return std::sqrt(nrm_sq(v)); }


// return the reverse
template <typename T> inline constexpr MVec3d<T> rev(MVec3d<T> const& v)
{
    // only bivector and trivector parts switch signs
    return MVec3d<T>(v.c0, v.c1, v.c2, v.c3, -v.c4, -v.c5, -v.c6, -v.c7);
}

// return the Clifford conjugate
template <typename T> inline constexpr MVec3d<T> conj(MVec3d<T> const& v)
{
    // only vector and bivector parts switch signs
    return MVec3d<T>(v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, v.c7);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of quaternion
// |Z|^2 = Z rev(Z) = c0^2 + c1^2 + c2^2 + c3^2
template <typename T> inline T nrm_sq(MVec3d_E<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3;
}

// return magnitude of quaternion
template <typename T> inline T nrm(MVec3d_E<T> const& v) { return std::sqrt(nrm_sq(v)); }

// return conjugate complex of quaternion (MVec3d_E<T>),
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline constexpr MVec3d_E<T> rev(MVec3d_E<T> const& v)
{
    // only the bivector part switches sign
    return MVec3d_E<T>(v.c0, -v.c1, -v.c2, -v.c3);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d_U<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return the reverse
template <typename T> inline constexpr MVec3d_U<T> rev(MVec3d_U<T> const& v)
{
    // only the trivector part switches signs
    return MVec3d_U<T>(v.c0, v.c1, v.c2, -v.c3);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 3d multivector
// Expensive! - Don't use if you don't have to! (64x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// a * b = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector*vector = scalar + bivector)
// A * b = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector*vector = vector + trivector)
// a * B = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)  (vector*bivector = vector + trivector)
//
// multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            MVec3d<U> const& B)
{
    // geometric product of two fully populated 3d multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
               A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 - A.c4 * B.c7 -
               A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c4;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c4 * B.c3 -
               A.c5 * B.c7 - A.c6 * B.c1 - A.c7 * B.c5;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
               A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0 -
               A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 +
               A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 - A.c4 * B.c5 +
               A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
               A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of an even grade multivector with a multivector
// even grade multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                            MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6;
    ctype c1 = A.c0 * B.c1 - A.c1 * B.c7 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c3 - A.c2 * B.c7 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c7;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4;
    ctype c6 = A.c0 * B.c6 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of a multivector with an even grade multivector
// multivector * even grade multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype c1 = A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2 - A.c7 * B.c1;
    ctype c2 = A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1 - A.c7 * B.c2;
    ctype c3 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0 - A.c7 * B.c3;
    ctype c4 = A.c0 * B.c1 + A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype c5 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1;
    ctype c6 = A.c0 * B.c3 - A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c0;
    ctype c7 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product A * B for two multivectors from the even subalgebra (3d case)
// even grade multivector * even grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(
        Scalar<ctype>(A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3),
        BiVec3d<ctype>(A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2,
                       A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1,
                       A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0));
}

// geometric product A * B for two multivectors from the uneven subalgebra (3d case)
// uneven grade multivector * uneven grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(
        Scalar<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3),
        BiVec3d<ctype>(A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0,
                       -A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1,
                       A.c0 * B.c1 - A.c1 * B.c0 + A.c2 * B.c3 + A.c3 * B.c2));
}

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a multivector B of the uneven subalgebra
// even grade multivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(
        Vec3d<ctype>(A.c0 * B.c0 - A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1,
                     A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c3 - A.c3 * B.c0,
                     A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3d<ctype>(A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2));
}

// geometric product A * B of a multivector A from the uneven subalgebra (3d case)
// with a multivector B of the even subalgebra
// uneven grade multivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(
        Vec3d<ctype>(A.c0 * B.c0 - A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1,
                     A.c0 * B.c3 + A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2,
                     -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3d<ctype>(A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0));
}

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar<ctype>(-A.c1 * B.x - A.c2 * B.y - A.c3 * B.z),
                           BiVec3d<ctype>(A.c0 * B.x - A.c2 * B.z + A.c3 * B.y,
                                          A.c0 * B.y + A.c1 * B.z - A.c3 * B.x,
                                          A.c0 * B.z - A.c1 * B.y + A.c2 * B.x));
}

// geometric product A * B of a bivector A with an even grade multivector B
// from the right
// bivector * even grade multivector => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar<ctype>(-A.x * B.c1 - A.y * B.c2 - A.z * B.c3),
                           BiVec3d<ctype>(A.x * B.c0 - A.y * B.c3 + A.z * B.c2,
                                          A.x * B.c3 + A.y * B.c0 - A.z * B.c1,
                                          -A.x * B.c2 + A.y * B.c1 + A.z * B.c0));
}

// geometric product A * b of an even grade multivector A with a vector b
// from the right
// even grade multivector * vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(A.c0 * b.x - A.c2 * b.z + A.c3 * b.y,
                                        A.c0 * b.y + A.c1 * b.z - A.c3 * b.x,
                                        A.c0 * b.z - A.c1 * b.y + A.c2 * b.x),
                           PScalar3d<ctype>(A.c1 * b.x + A.c2 * b.y + A.c3 * b.z));
}

// geometric product a * B of a vector a and an even grade multivector B from the left
// vector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(a.x * B.c0 - a.y * B.c3 + a.z * B.c2,
                                        a.x * B.c3 + a.y * B.c0 - a.z * B.c1,
                                        -a.x * B.c2 + a.y * B.c1 + a.z * B.c0),
                           PScalar3d<ctype>(a.x * B.c1 + a.y * B.c2 + a.z * B.c3));
}

// geometric product A * B of an uneven grade multivector A with a bivector B
// uneven grade multivector * bivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-A.c1 * B.z + A.c2 * B.y - A.c3 * B.x,
                                        A.c0 * B.z - A.c2 * B.x - A.c3 * B.y,
                                        -A.c0 * B.y + A.c1 * B.x - A.c3 * B.z),
                           PScalar3d<ctype>(A.c0 * B.x + A.c1 * B.y + A.c2 * B.z));
}

// geometric product A * B of a bivector A with an uneven grade multivector B
// bivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-A.x * B.c3 - A.y * B.c2 + A.z * B.c1,
                                        +A.x * B.c2 - A.y * B.c3 - A.z * B.c0,
                                        -A.x * B.c1 + A.y * B.c0 - A.z * B.c3),
                           PScalar3d<ctype>(A.x * B.c0 + A.y * B.c1 + A.z * B.c2));
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
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& a,
                                                              BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar<ctype>(dot(a, b)), cmt(a, b));
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
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(dot(A, b), wdg(A, b));
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
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(dot(a, B), wdg(a, B));
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
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar<ctype>(dot(a, b)), wdg(a, b));
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
// to the multivector B
// trivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                            MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec3d<ctype>(-B.c7, -B.c4, -B.c5, -B.c6, B.c1, B.c2, B.c3, B.c0);
}

// geometric product A * B of a multivector A with the trivector B from the right
// multivector * trivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(-A.c7, -A.c4, -A.c5, -A.c6, A.c1, A.c2, A.c3, A.c0) * ctype(B);
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
// to the even grade multivector B
// trivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) *
           MVec3d_U<ctype>(Vec3d<ctype>(-B.c1, -B.c2, -B.c3), PScalar3d<ctype>(B.c0));
}

// geometric product A * B of an even multivector A multiplied from the right
// by the trivector B (=3d pseudoscalar)
// even grade multivector * trivector => uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<U> const& A,
                                                              PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-A.c1, -A.c2, -A.c3), PScalar3d<ctype>(A.c0)) *
           ctype(B);
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
// to the uneven grade multivector B
// trivector * uneven grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) *
           MVec3d_E<ctype>(Scalar<ctype>(-B.c3), BiVec3d<ctype>(B.c0, B.c1, B.c2));
}

// geometric product A * B of an uneven grade multivector A multiplied from the right
// by the trivector B (=3d pseudoscalar)
// uneven grade multivector * trivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<U> const& A,
                                                              PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar<ctype>(-A.c3), BiVec3d<ctype>(A.c0, A.c1, A.c2)) *
           ctype(B);
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                           BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * Vec3d<ctype>(B.x, B.y, B.z);
}

// geometric product A * B of bivector A multiplied by a trivector B from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                           PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return -Vec3d<ctype>(A.x, A.y, A.z) * ctype(B);
}

// geometric product A * b of a trivector A (=3d pseudoscalar) multiplied from the left
// to the vector b
// trivector * vector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                             Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * BiVec3d<ctype>(b.x, b.y, b.z);
}

// geometric product a * B of a vector a multiplied with a trivector B from the right
// vector * trivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                             PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(a.x, a.y, a.z) * ctype(B);
}

// geometric product A * B of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> operator*(PScalar3d<T> A, PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * ctype(B); // trivectors in 3d square to -1
}

// return the multiplicative inverse of the multivector
// inv(M) = 1/( M*conj(M) ) * conj(M)  with M*conj(M) being a scalar value
template <typename T> inline MVec3d<T> inv(MVec3d<T> const& v)
{
    // alternative, but with slightly more computational effort:
    T m_conjm = gr0(v * conj(v));
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(m_conjm) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("multivector norm too small for inversion " +
                                 std::to_string(m_conjm) + "\n");
    }
#endif
    T inv = T(1.0) / m_conjm; // inverse of squared norm for a vector
    return inv * conj(v);
}
// ATTENTION: there is a left and a right inverse (see paper of Hitzer, Sangwine)

////////////////////////////////////////////////////////////////////////////////
// 3d rotation operations
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
// exponential function with bivector as argument for setup of quaternions
// as geometric multivector with a scalar and a bivector part
// MVec3d_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
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
inline constexpr MVec3d_E<std::common_type_t<T, U>> exp(BiVec3d<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar<ctype>(std::cos(theta)),
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
inline constexpr MVec3d_E<std::common_type_t<T, U>> rotor(BiVec3d<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec3d_E<ctype>(Scalar<ctype>(std::cos(half_angle)),
                           normalize(I) * std::sin(half_angle));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rotate(Vec3d<T> const& v,
                                                        MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> rr = rev(rotor);
    // MVec3d_U<ctype> tmp = rotor * v;
    // MVec3d_U<ctype> res = tmp * rr;

    // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the pseudoscalar part
    // which will be zero anyway
    // return Vec3d<ctype>(gr1<ctype>(rotor * v * rev(rotor)));

    // optimized version:
    MVec3d_E<ctype> rr = rev(rotor);
    MVec3d_U<ctype> tmp = rotor * v;
    // formular from operator*(MVec3d_U<T>, MVec3d_E<U>) - only vector part
    return Vec3d<ctype>(tmp.c0 * rr.c0 - tmp.c1 * rr.c3 + tmp.c2 * rr.c2 - tmp.c3 * rr.c1,
                        tmp.c0 * rr.c3 + tmp.c1 * rr.c0 - tmp.c2 * rr.c1 - tmp.c3 * rr.c2,
                        -tmp.c0 * rr.c2 + tmp.c1 * rr.c1 + tmp.c2 * rr.c0 -
                            tmp.c3 * rr.c3);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rotate(BiVec3d<T> const& v,
                                                          MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> rr = rev(rotor);
    // MVec3d_E<ctype> tmp = rotor * v;
    // MVec3d_E<ctype> res = tmp * rr;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway
    // return BiVec3d<ctype>(gr2<ctype>(rotor * v * rev(rotor)));

    // optimized version:
    MVec3d_E<ctype> rr = rev(rotor);
    MVec3d_E<ctype> tmp = rotor * v;
    // formular from operator*(MVec3d_E<T>, MVec3d_E<U>) - only bivector part
    return BiVec3d<ctype>(
        tmp.c0 * rr.c1 + tmp.c1 * rr.c0 - tmp.c2 * rr.c3 + tmp.c3 * rr.c2,
        tmp.c0 * rr.c2 + tmp.c1 * rr.c3 + tmp.c2 * rr.c0 - tmp.c3 * rr.c1,
        tmp.c0 * rr.c3 - tmp.c1 * rr.c2 + tmp.c2 * rr.c1 + tmp.c3 * rr.c0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> rotate(MVec3d<T> const& v,
                                                         MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(rotor * v * rev(rotor));
}

////////////////////////////////////////////////////////////////////////////////
// 3d duality operations
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B as subspace of R^3 then
// dual(M) represents the subspace orthorgonal to B
// => return the dual(M) of the multivector M

#if defined(_HD_GA_HESTENES_DORAN_LASENBY_DUALITY)
////////////////////////////////////////////////////////////////////////////////
// duality (as defined by Hestenes or by Doran, Lasenby in "GA for physicists"):
// (same subspace as for Macdonld's definition, but other sign)
// (will have influence on all formulae concerning duality)
////////////////////////////////////////////////////////////////////////////////
//
// dual(A) = I*A
//
// I_3d * 1 = e1^e2^e3 * 1 = e1^e2^e3
//
// I_3d * e1 = e1^e2^e3 * e1 = e_1231 =  e_1123 = e_23         = e_2^e_3
// I_3d * e2 = e1^e2^e3 * e2 = e_1232 = -e_1223 = -e_13 = e_31 = e_3^e_1
// I_3d * e3 = e1^e2^e3 * e3 = e_1233 =  e_12                  = e_1^e_2
//
// I_3d * e2^e3 = e1^e2^e3 * e2^e3 = e_12323 = -e_1            = -e1
// I_3d * e3^e1 = e1^e2^e3 * e3^e1 = e_12331 = -e_11233 = -e_2 = -e2
// I_3d * e1^e2 = e1^e2^e3 * e1^e2 = e_12312 = -e_11223 = -e_3 = -e3
//
// I_3d * e1^e2^e3 = e1^e2^e3 * e1^e2^e3 = e_123123 = -1

// using this duality definition, following identities hold for the cross-product:
//
//    a x b  = -dual(a ^ b) = -I_3d*(a ^ b) = -(a ^ b)*I_3d =  (a ^ b)*rev(I_3d)
//    a ^ b  =  dual(a x b) =  I_3d*(a x b) =  (a x b)*I_3d = -(a x b)*rev(I_3d)
//
//    the last 3 columns in both lines are valid independent of the duality definition!
//    the sign issue only occurs due to the duality definition itself

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar<T> dual3d(PScalar3d<T> ps)
{
    // dual(A) = I*A
    // e123 * (ps * e123) = -ps
    return Scalar<T>(-T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform to different pseudoscalars
// the 2d and 3d adders in the function name are required for disambiguation
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual3d(Scalar<T> s)
{
    // dual(A) = I*A
    // e123 * (s) = s * e123
    return PScalar3d<T>(T(s));
}

// this overload is provided to accept T directly as an alternative to Scalar<T>
// e.g. with T as a result of a dot product between two vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual3d(T s)
{
    // dual(A) = I*A
    // e123 * (s) = s * e123
    return PScalar3d<T>(s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> dual3d(Vec3d<T> const& v)
{
    // dual(A) = I*A
    // e123 * (v.x * e1  + v.y * e2  + v.z * e3)
    //      =  v.x * e23 + v.y * e31 + v.z * e12
    return BiVec3d<T>(v.x, v.y, v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> dual3d(BiVec3d<T> const& B)
{
    // dual(A) = I*A
    // e123 * (  b.x * e23 + b.y * e31 + b.z * e12)
    //      =  - b.x * e1  - b.y * e2  - b.z * e3
    return Vec3d<T>(-B.x, -B.y, -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> dual3d(MVec3d_E<T> const& M)
{
    // dual(A) = I*A
    // e123 * (s + b.x * e23 + b.y * e31 + b.z * e12)
    //      =    - b.x * e1  - b.y * e2  - b.z * e3 + s * e123
    return MVec3d_U<T>(-M.c1, -M.c2, -M.c3, M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> dual3d(MVec3d_U<T> const& M)
{
    // dual(A) = I*A
    // e123 * (      v.x * e1  + v.y * e2  + v.z * e3 + ps * e123)
    //      =  -ps + v.x * e23 + v.y * e31 + v.z * e12
    return MVec3d_E<T>(-M.c3, M.c0, M.c1, M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> dual3d(MVec3d<T> const& M)
{
    // dual(A) = I*A
    // e123 * (  s + v.x * e1  + v.y * e2  + v.z * e3
    //             + b.x * e23 + b.y * e31 + b.z * e12 + ps * e123)
    //      =  -ps - b.x * e1  - b.y * e2  - b.z * e3
    //             + v.x * e23 + v.y * e31 + v.z * e12 +  s * e123
    return MVec3d<T>(-M.c7, -M.c4, -M.c5, -M.c6, M.c1, M.c2, M.c3, M.c0);
}

#else
////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Macdonald, "Linear and geometric algebra"):
////////////////////////////////////////////////////////////////////////////////
//
// dual(A) = A/I = A*I^(-1) = A*rev(I)
//
// 1 * rev(I_3d) = 1 * e3^e2^e1 = e3^e2^e1 = e_321 = -e_123 = -e1^e2^e3
//
// e1 * rev(I_3d) = e1 * e3^e2^e1 = e_1321 =  e_1132 = -e_23 = -e2^e3
// e2 * rev(I_3d) = e2 * e3^e2^e1 = e_2321 = -e_2231 = -e_31 = -e3^e1
// e3 * rev(I_3d) = e3 * e3^e2^e1 = e_3321 = -e_3312 = -e_12 = -e1^e2
//
// e2^e3 * rev(I_3d) = e2^e3 * e3^e2^e1 = e_23321 = e_1           = e1
// e3^e1 * rev(I_3d) = e3^e1 * e3^e2^e1 = e_31321 = e_33112 = e_2 = e2
// e1^e2 * rev(I_3d) = e1^e2 * e3^e2^e1 = e_12321 = e_11223 = e_3 = e3
//
// e1^e2^e3 * rev(I_3d) = e1^e2^e3 * e3^e2^e1 = e_123321= 1

// using this duality definition, following duality properties hold
// (A. Macdonald, p. 110):
//
// a) dual(aA) = a dual(A)
// b) dual(A + B) = dual(A) + dual(B)
// c) dual(dual(A)) = (-1)^(n*(n-1)/2) A   (with n as dimension of the (sub)space)
// d) |dual(B)| = |B|
// e) if B is a j-blade then dual(B) is an (n-j)-blade
// f) if A is a j-vector then dual(A) is an (n-j)-vector
//    (remember: a j-vector is a sum of j-blades, which are outer products)

// using this duality definition, following identities hold for the cross-product:
//
//    a x b  =  dual(a ^ b) =  (a ^ b)*rev(I_3d) = -I_3d*(a ^ b) = -(a ^ b)*I_3d
//    a ^ b  = -dual(a x b) = -(a x b)*rev(I_3d) =  I_3d*(a x b) =  (a x b)*I_3d
//
//    the last 3 columns in both lines are valid independent of the duality definition!
//    the sign issue only occurs due to the duality definition itself

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar<T> dual3d(PScalar3d<T> ps)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    // (ps * e123) * e321 = ps
    return Scalar<T>(T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform to different pseudoscalars
// the 2d and 3d adders in the function name are required for disambiguation
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual3d(Scalar<T> s)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    // (s) * e321 = -s * e123
    return PScalar3d<T>(-T(s));
}

// this overload is provided to accept T directly as an alternative to Scalar<T>
// e.g. with T as a result of a dot product between two vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual3d(T s)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    // (s) * e321 = -s * e123
    return PScalar3d<T>(-s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> dual3d(Vec3d<T> const& v)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   ( v.x * e1  + v.y * e2  + v.z * e3) * e321
    // =  -v.x * e23 - v.y * e31 - v.z * e12
    return BiVec3d<T>(-v.x, -v.y, -v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> dual3d(BiVec3d<T> const& B)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (b.x * e23 + b.y * e31 + b.z * e12) * e321
    // =  b.x * e1  + b.y * e2  + b.z * e3
    return Vec3d<T>(B.x, B.y, B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> dual3d(MVec3d_E<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (s + b.x * e23 + b.y * e31 + b.z * e12 ) * e321
    // =      b.x * e1  + b.y * e2  + b.z * e3 - s * e123
    return MVec3d_U<T>(M.c1, M.c2, M.c3, -M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> dual3d(MVec3d_U<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (     v.x * e1  + v.y * e2  + v.z * e3 + ps * e123) * e321
    // =  ps - v.x * e23 - v.y * e31 - v.z * e12
    return MVec3d_E<T>(M.c3, -M.c0, -M.c1, -M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> dual3d(MVec3d<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //  (  s + v.x * e1  + v.y * e2  + v.z * e3
    //       + b.x * e23 + b.y * e31 + b.z * e12 + ps * e123) * e321
    // =  ps + b.x * e1  + b.y * e2  + b.z * e3
    //       - v.x * e23 - v.y * e31 - v.z * e12 -  s * e123
    return MVec3d<T>(M.c7, M.c4, M.c5, M.c6, -M.c1, -M.c2, -M.c3, -M.c0);
}
#endif


////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> and BiVec3d<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v1,
                                                              Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return dot(v1, v2) * Vec3d<ctype>(inv(v2));
}

// projection of v1 onto v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto_unitized(Vec3d<T> const& v1,
                                                                       Vec3d<U> const& v2)
{
    return dot(v1, v2) * v2; // v2 is already its own reverse
}

// projection of a vector v1 onto a bivector v2
// v_parallel = dot(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v1,
                                                              BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    Vec3d<ctype> a = dot(v1, v2);
    BiVec3d<ctype> Bi = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(a * Bi) = dot(a,Bi)
    return Vec3d<ctype>(dot(a, Bi));
}

// projection of a vector v1 onto a normalized bivector v2
// u_parallel = gr1(dot(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
project_onto_unitized(Vec3d<T> const& v1, BiVec3d<U> const& v2)
{
    // requires v2 to be normalized

    using ctype = std::common_type_t<T, U>;
    Vec3d<ctype> a = dot(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3d<ctype> Bi = -v2;
    // use the formular equivalent to the geometric product to save computational cost
    // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(a * Bi) = dot(a,Bi)
    return Vec3d<ctype>(dot(a, Bi));
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v1,
                                                             Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    BiVec3d<ctype> B = wdg(v1, v2);
    Vec3d<ctype> v2_inv = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    return Vec3d<ctype>(dot(B, v2_inv));
}

// rejection of vector v1 from a normalized vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from_unitized(Vec3d<T> const& v1,
                                                                      Vec3d<U> const& v2)
{
    // requires v2 to be normalized
    using ctype = std::common_type_t<T, U>;
    BiVec3d<ctype> B = wdg(v1, v2);
    Vec3d<ctype> v2_inv = v2; // v2 is its own inverse, if normalized
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    return Vec3d<ctype>(dot(B, v2_inv));
}

// rejection of vector v1 from a bivector v2
// u_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v1,
                                                             BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3d<ctype> A = wdg(v1, v2);
    BiVec3d<ctype> B = inv(v2);
    // trivector * bivector = vector
    return A * B;
}

// rejection of vector v1 from a normalized bivector v2
// u_perp = wdg(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
reject_from_unitized(Vec3d<T> const& v1, BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3d<ctype> a = wdg(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3d<ctype> B = -v2;
    // trivector * bivector = vector (derived from full geometric product to save
    // costs)
    return Vec3d<ctype>(-a * B.x, -a * B.y, -a * B.z);
}

// reflect a vector u on a hyperplane B orthogonal to vector b
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n (a line in 2d space)
// orthogonal to vector b: the hyperplane is dual to b (i.e. a one dimensional subspace)
//
// HINT: choose b * B = I_3d  =>  B = b * I_3d  (for normalized b)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reflect_on_hyp(Vec3d<T> const& u,
                                                                Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(-b * u * inv(b)));
}

// reflect a vector u in an arbitrary bivector, i.e. a plane
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reflect_on(Vec3d<T> const& u,
                                                            BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(-B * u * inv(B)));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling planes)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> reflect_on(BiVec3d<T> const& UB,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(gr2(B * UB * inv(B)));
}

// reflect a vector u another vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reflect_on_vec(Vec3d<T> const& u,
                                                                Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(b * u * inv(b)));
}

////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt-Orthogonalization part 1: 2d plane embedded in 3d space
////////////////////////////////////////////////////////////////////////////////
//
// input:  two linear independent vectors u and v in 3d defining a plane
// output: two orthogonal vectors with the first one being u and the second one a vector
// perpendicular to u in the orientation of v, both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
std::vector<Vec3d<std::common_type_t<T, U>>> gs_orthogonal(Vec3d<T> const& u,
                                                           Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec3d<ctype>> basis;
    basis.push_back(u);
    basis.emplace_back(reject_from(v, u));
    return basis;
}

// input:  two linear independent vectors u and v in 3d defining a plane
// output: two orthonormal vectors with the first one being normalize(u) and
// the second one a normalized vector perpendicular to u in the orientation of v,
// both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
std::vector<Vec3d<std::common_type_t<T, U>>> gs_orthonormal(Vec3d<T> const& u,
                                                            Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec3d<ctype>> basis;
    Vec3d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from_unitized(v, u_unitized)));
    return basis;
}

////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt-Orthogonalization part 2: 3d space
////////////////////////////////////////////////////////////////////////////////
//
// input:  three linear independent vectors u, v and w in 3d
// output: three orthogonal vectors with the first one being u and the second and third
// being  perpendicular to u and the plane spanned by u and v respectively.
// All three from an orthogonal system
template <typename U, typename V, typename W>
    requires(std::floating_point<U> && std::floating_point<V> && std::floating_point<W>)
std::vector<Vec3d<std::common_type_t<U, V, W>>>
gs_orthogonal(Vec3d<U> const& u, Vec3d<V> const& v, Vec3d<W> const& w)
{
    using ctype = std::common_type_t<U, V, W>;
    std::vector<Vec3d<ctype>> basis;
    basis.push_back(u);
    basis.emplace_back(reject_from(v, u));
    basis.emplace_back(reject_from(w, wdg(u, v)));
    return basis;
}

// input:  three linear independent vectors u, v and w in 3d
// output: three orthonormal vectors with the first one being normalized u and the second
// and third being normalized and perpendicular to u and the plane spanned by u and v
// respectively. All three from an orthogonal system
template <typename U, typename V, typename W>
    requires(std::floating_point<U> && std::floating_point<V> && std::floating_point<W>)
std::vector<Vec3d<std::common_type_t<U, V, W>>>
gs_orthonormal(Vec3d<U> const& u, Vec3d<V> const& v, Vec3d<W> const& w)
{
    using ctype = std::common_type_t<U, V, W>;
    std::vector<Vec3d<ctype>> basis;
    Vec3d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from_unitized(v, u_unitized)));
    basis.emplace_back(normalize(reject_from(w, wdg(u, v))));
    return basis;
}

} // namespace hd::ga::ega