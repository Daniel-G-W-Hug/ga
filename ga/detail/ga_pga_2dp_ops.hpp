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

#include "ga_value_t.hpp"

#include "ga_mvec2dp.hpp" // inclusion of multivector imports all component types


namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// Vec2dp<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return dot-product of two vectors
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    // this assumes an orhronormal basis with e1^2 = 1, e2^2 = 1, e3^2 = 0
    // thus the z-component does not matter
    return v1.x * v2.x + v1.y * v2.y;
}

// return squared magnitude of vector
template <typename T> inline T sq_nrm(Vec2dp<T> const& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

// return magnitude of vector
template <typename T> inline T nrm(Vec2dp<T> const& v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// return a vector normalized to nrm(v) == 1.0
template <typename T> inline Vec2dp<T> normalize(Vec2dp<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return a vector unitized to v.z == 1.0
template <typename T> inline Vec2dp<T> unitize(Vec2dp<T> const& v)
{
    T n = v.z;
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for unitization" +
                                 std::to_string(n) + "\n");
    }
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec2dp<T>(v.x * inv, v.y * inv, T(1.0));
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec2dp<T> inv(Vec2dp<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return Vec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// wedge product (returns a bivector)
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec2dp<std::common_type_t<T, U>> wdg(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    return BiVec2dp<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// return the angle between of two vectors
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    using std::numbers::pi;

    ctype nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }

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


////////////////////////////////////////////////////////////////////////////////
// BiVec3dp<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return dot product of two bivectors A and B (= a scalar)
// dot(A,B) = gr0(A * B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> dot(BiVec2dp<T> const& A, BiVec2dp<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    return -A.z * B.z;
}

// return squared magnitude of bivector
template <typename T> inline constexpr T sq_nrm(BiVec2dp<T> const& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

// return magnitude of bivector
template <typename T> inline constexpr T nrm(BiVec2dp<T> const& v)
{
    return std::sqrt(sq_nrm(v));
}

// return a bivector normalized to nrm(v) == 1.0
template <typename T> inline constexpr BiVec2dp<T> normalize(BiVec2dp<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return BiVec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the bivector
template <typename T> inline constexpr BiVec2dp<T> inv(BiVec2dp<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = -T(1.0) / sq_n; // negative inverse of squared norm for a bivector
    return BiVec2dp<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return conjugate complex of a bivector
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline constexpr BiVec2dp<T> rev(BiVec2dp<T> const& v)
{
    // all bivector parts switch sign
    return BiVec2dp<T>(-v.x, -v.y, -v.z);
}

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(BiVec2dp<T> const& v1, BiVec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
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
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
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
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(
        std::clamp(ctype(nrm(dot(v1, v2))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// wedge product between a vector a and a bivector B
// wdg(a,B) = gr3(a * B)
// => returns a trivector
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

////////////////////////////////////////////////////////////////////////////////
// MVec2dp<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude
// |M|^2 = M rev(M) = (M.c0)^2 + (M.c1)^2 + (M.c2)^2 + (M.c3)^3
//                  + (M.c4)^2 + (M.c5)^2 + (M.c6)^2 + (M.c7)^3
template <typename T> inline T sq_nrm(MVec2dp<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3 + v.c4 * v.c4 +
           v.c5 * v.c5 + v.c6 * v.c6 + v.c7 * v.c7;
}

// return magnitude of complex number
template <typename T> inline T nrm(MVec2dp<T> const& v) { return std::sqrt(sq_nrm(v)); }


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


// return a multivector normalized to nrm(v) == 1.0
template <typename T> inline MVec2dp<T> normalize(MVec2dp<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return MVec2dp<T>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv, v.c4 * inv,
                      v.c5 * inv, v.c6 * inv, v.c7 * inv);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of quaternion
// |Z|^2 = Z rev(Z) = c0^2 + c1^2 + c2^2 + c3^2
template <typename T> inline T sq_nrm(MVec2dp_E<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3;
}

// return magnitude of quaternion
template <typename T> inline T nrm(MVec2dp_E<T> const& v) { return std::sqrt(sq_nrm(v)); }

// return conjugate complex of quaternion (MVec2dp_E<T>),
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline constexpr MVec2dp_E<T> rev(MVec2dp_E<T> const& v)
{
    // only the bivector part switches sign
    return MVec2dp_E<T>(v.c0, -v.c1, -v.c2, -v.c3);
}

// return a complex normalized to nrm(v) == 1.0
template <typename T> inline MVec2dp_E<T> normalize(MVec2dp_E<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return MVec2dp_E<T>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv);
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
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c6 * B.c6;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c6 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c6 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
               A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0 -
               A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 +
               A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 - A.c4 * B.c5 +
               A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3;
    ctype c7 = -A.c0 * B.c7 - A.c1 * B.c4 - A.c2 * B.c5 - A.c6 * B.c3;
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
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c7;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4;
    ctype c6 = A.c0 * B.c6 + A.c3 * B.c0;
    ctype c7 = -A.c0 * B.c7 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3;
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
    ctype c3 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0 - A.c7 * B.c3;
    ctype c4 = A.c0 * B.c1 + A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype c5 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1;
    ctype c6 = A.c0 * B.c3 + A.c6 * B.c0;
    ctype c7 = -A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3 - A.c7 * B.c0;
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
        Scalar<ctype>(A.c0 * B.c0 - A.c3 * B.c3),
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
        Scalar<ctype>(A.c0 * B.c0 + A.c1 * B.c1),
        BiVec2dp<ctype>(A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0,
                        -A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1,
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
                      A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar2dp<ctype>(-A.c0 * B.c3 - A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2));
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
                      -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar2dp<ctype>(-A.c0 * B.c1 - A.c1 * B.c2 - A.c2 * B.c3 - A.c3 * B.c0));
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
    return MVec2dp_E<ctype>(Scalar<ctype>(-A.c3 * B.z),
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
    return MVec2dp_E<ctype>(Scalar<ctype>(-A.z * B.c3),
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
                            PScalar2dp<ctype>(-A.c1 * b.x - A.c2 * b.y - A.c3 * b.z));
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
                            PScalar2dp<ctype>(-a.x * B.c1 - a.y * B.c2 - a.z * B.c3));
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
        Vec2dp<ctype>(-A.c1 * B.z, A.c0 * B.z, -A.c0 * B.y + A.c1 * B.x - A.c3 * B.z),
        PScalar2dp<ctype>(-A.c0 * B.x - A.c1 * B.y - A.c2 * B.z));
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
        Vec2dp<ctype>(A.z * B.c1, -A.z * B.c0, -A.x * B.c1 + A.y * B.c0 - A.z * B.c3),
        PScalar2dp<ctype>(-A.x * B.c0 - A.y * B.c1 - A.z * B.c2));
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
    return MVec2dp_E<ctype>(Scalar<ctype>(dot(a, b)), cmt(a, b));
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
    return MVec2dp_E<ctype>(Scalar<ctype>(dot(a, b)), wdg(a, b));
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
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
    return -ctype(B) * MVec2dp<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), -A.c6, A.c1,
                                      A.c2, ctype(0.0), -A.c0);
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
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
// by the trivector B (=3d pseudoscalar)
// even grade multivector * trivector => uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> operator*(MVec2dp_E<U> const& A,
                                                               PScalar2dp<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(B) * MVec2dp_U<ctype>(Vec2dp<ctype>(ctype(0.0), ctype(0.0), -A.c3),
                                        PScalar2dp<ctype>(-A.c0));
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
// to the uneven grade multivector B
// trivector * uneven grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(PScalar2dp<T> A,
                                                               MVec2dp_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) *
           MVec2dp_E<ctype>(Scalar<ctype>(B.c3), BiVec2dp<ctype>(B.c0, B.c1, ctype(0.0)));
}

// geometric product A * B of an uneven grade multivector A multiplied from the right
// by the trivector B (=3d pseudoscalar)
// uneven grade multivector * trivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_U<U> const& A,
                                                               PScalar2dp<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(B) *
           MVec2dp_E<ctype>(Scalar<ctype>(A.c3), BiVec2dp<ctype>(A.c0, A.c1, ctype(0.0)));
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
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
    return -ctype(B) * Vec2dp<ctype>(ctype(0.0), ctype(0.0), -A.z);
}

// geometric product A * b of a trivector A (=3d pseudoscalar) multiplied from the left
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
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                             PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(B) * BiVec3d<ctype>(a.x, a.y, ctype(0.0));
}

// geometric product A * B of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> operator*([[maybe_unused]] PScalar3d<T> A,
                                                    [[maybe_unused]] PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(0.0); // trivectors in 3d pga contain e3^2 = 0
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
    return MVec2dp_E<ctype>(Scalar<ctype>(std::cos(theta)),
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
    return MVec2dp_E<ctype>(Scalar<ctype>(std::cos(half_angle)),
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

} // namespace hd::ga::pga