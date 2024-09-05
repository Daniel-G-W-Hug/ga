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

#include "ga_mvec2d.hpp" // inclusion of multivector imports all component types

namespace hd::ga::ega {

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return dot-product of two vectors in G<2,0,0>
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    // definition: dot(v1, v2) = (v1)^T g_12 v2 with the metric g_12
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1
    // as diagonal elements of g_12
    return v1.x * v2.x + v1.y * v2.y;
}

// return squared magnitude of vector
template <typename T> inline T nrm_sq(Vec2d<T> const& v) { return dot(v, v); }

// return magnitude of vector
template <typename T> inline T nrm(Vec2d<T> const& v) { return std::sqrt(dot(v, v)); }

// return the multiplicative inverse of the vector
template <typename T> inline Vec2d<T> inv(Vec2d<T> const& v)
{
    T sq_n = nrm_sq(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// wedge product (returns a bivector, which is the pseudoscalar in 2d)
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    return PScalar2d<std::common_type_t<T, U>>(v1.x * v2.y - v1.y * v2.x);
}

// return the angle between of two vectors
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    using std::numbers::pi;

    ctype nrm_prod = nrm(v1) * nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
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

////////////////////////////////////////////////////////////////////////////////
// MVec2d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude
// |M|^2 = M rev(M) = (M.c0)^2 + (M.c1)^2 + (M.c2)^2 + (M.c3)^3
template <typename T> inline T nrm_sq(MVec2d<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3;
}

// return magnitude
template <typename T> inline T nrm(MVec2d<T> const& v) { return std::sqrt(nrm_sq(v)); }

// return the reverse
template <typename T> inline MVec2d<T> rev(MVec2d<T> const& v)
{
    return MVec2d<T>(v.c0, v.c1, v.c2, -v.c3);
}

// return the Clifford conjugate
template <typename T> inline MVec2d<T> conj(MVec2d<T> const& v)
{
    return MVec2d<T>(v.c0, -v.c1, -v.c2, -v.c3);
}

// return the multiplicative inverse of the multivector
// inv(M) = 1/( M*conj(M) ) * conj(M)  with M*conj(M) being a scalar value
template <typename T> inline MVec2d<T> inv(MVec2d<T> const& v)
{
    // from manual calculation of M*conj(M) in 2d:
    T m_conjm = v.c0 * v.c0 + v.c3 * v.c3 - nrm_sq(Vec2d<T>(v.c1, v.c2));
    //
    // alternative, but with slightly more computational effort:
    // T m_conjm = gr0(v * conj(v));
    //

#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(m_conjm) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("multivector norm too small for inversion " +
                                 std::to_string(m_conjm) + "\n");
        // example: MVec2D(1,1,1,1) is not invertible
    }
#endif
    T inv = T(1.0) / m_conjm; // inverse of squared norm for a vector
    return inv * conj(v);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2d_E<T> geometric operations for complex numbers
//             (= multivectors from the even subalgebra)
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of complex number
// |Z|^2 = Z rev(Z) = c0^2 + c1^2
template <typename T> inline T nrm_sq(MVec2d_E<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1;
}

// return magnitude of complex number
template <typename T> inline T nrm(MVec2d_E<T> const& v) { return std::sqrt(nrm_sq(v)); }

// return conjugate complex of a complex number,
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline MVec2d_E<T> rev(MVec2d_E<T> const& v)
{
    return MVec2d_E<T>(v.c0, -v.c1);
}

// return the multiplicative inverse of the complex number (inv(z) = 1/nrm_sq(z)*rev(z))
// with rev(z) being the complex conjugate
template <typename T> inline MVec2d_E<T> inv(MVec2d_E<T> const& v)
{
    T sq_n = nrm_sq(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return inv * rev(v);
}

// return the angle of the complex number w.r.t. the real axis
// range of angle: -pi <= angle <= pi
template <typename T>
    requires(std::floating_point<T>)
inline T angle_to_re(MVec2d_E<T> const& v)
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
// PScalar2d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return reverse of a 2d bivector
template <typename T> inline PScalar2d<T> rev(PScalar2d<T> A)
{
    // the 2d bivector switches sign on reversion
    return PScalar2d<T>(-T(A));
}

// return inverse of the pseudoscalar (A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// k is the dimension of the space of the pseudoscalar formed by k orthogonal vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> inv(PScalar2d<T> ps)
{
    return -PScalar2d<T>(ps) / nrm_sq(ps);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2d<T> geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product ab for fully populated 2d multivector
// Expensive! - Don't use if you don't have to! (16x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors, e.g.:
// a*b = dot(a,b) + wdg(a,b) = gr0(a*b) + gr2(a*b) (vector vector = scalar + bivector)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // geometric product with a fully populated 2d multivector
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// geometric product ab for two vectors (returns a multivector of the even subalgebra)
// a*b = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
// => vector vector = scalar + bivector
//
// HINT: if a full 2d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                              Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(Scalar2d<ctype>(dot(a, b)), wdg(a, b));
}

// geometric product A*B of a 2d pseudoscalar (=bivector) A multiplied from the left
// to the multivector B
// 2d pseudoscalar (=bivector) * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec2d<ctype>(-B.c3, B.c2, -B.c1, B.c0);
}

// geometric product A*B of a bivector A multiplied from the left
// to the multivector from the even subalgebra B (MVec2d_E)
// bivector * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                              MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec2d_E<ctype>(-B.c1, B.c0);
}

// geometric product A*b of a bivector A multiplied from the left
// to the vector b
// bivector * vector => vector
// this multiplication turns the vector by -90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                           Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * Vec2d<ctype>(b.y, -b.x);
}

// geometric product A*B of a multivector A multiplied from the right by
// the pseudoscalar (=bivector) B
// multivector * bivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-A.c3, -A.c2, A.c1, A.c0) * ctype(B);
}

// geometric product A*B of an even grade multivector A multiplied from the right
// with the 2d pseudoscalar (=bivector) B
// even grade multivector * 2d pseudoscalar (=bivector) => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                              PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(-A.c1, A.c0) * ctype(B);
}

// geometric product a*B of the vector a multiplied from the right
// with the 2d pseudoscalar (=bivector) B
// vector * 2d pseudoscalar (=bivector) => vector
// this multiplication turns the vector by +90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                           PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-a.y, a.x) * ctype(B);
}

// geometric product a*B for a a vector a with a full 2d multivector B
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(a.x * B.c1 + a.y * B.c2, a.x * B.c0 - a.y * B.c3,
                         a.x * B.c3 + a.y * B.c0, a.x * B.c2 - a.y * B.c1);
}

// geometric product A*B for an even grade multivector A with a full 2d multivector B
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

// geometric product A*b for a multivector from the even subalgebra A
// with a vector b
// even grade multivector * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                           Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(A.c0 * b.x + A.c1 * b.y, A.c0 * b.y - A.c1 * b.x);
}

// geometric product A*B for a multivector A with an even grade multivector B
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

// geometric product A*b for a multivector A with a vector b
// multivector * vector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c1 * b.x + A.c2 * b.y, A.c0 * b.x + A.c3 * b.y,
                         -A.c3 * b.x + A.c0 * b.y, -A.c2 * b.x + A.c1 * b.y);
}

// geometric product a*B of a vector a with an even grade multivector B
// vector * even grade multivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                           MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(a.x * B.c0 - a.y * B.c1, a.x * B.c1 + a.y * B.c0);
}

// geometric product A*B for two multivectors from the even subalgebra (2d case)
// even grade multivector * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                              MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(A.c0 * B.c0 - A.c1 * B.c1, A.c0 * B.c1 + A.c1 * B.c0);
}

// geometric product A * B of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                              PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(-ctype(A) * ctype(B)); // bivectors in 2d square to -1
}

// geometric product A * B of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar2d<std::common_type_t<T, U>> operator*(Scalar2d<T> A,
                                                              Scalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(A) * ctype(B));
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
// 2d duality operations
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B as subspace of R^2 then
// dual(M) represents the subspace orthorgonal to B
// => return the dual(M) of the multivector M

#if defined(_HD_GA_HESTENES_DORAN_LASENBY_DUALITY)
////////////////////////////////////////////////////////////////////////////////
// duality as defined by Hestenes or by Doran, Lasenby in "GA for physicists":
// (same subspace as for Macdonld's definition below, but other resulting signs)
// (=> will have influence on formulae concerning duality directly)
////////////////////////////////////////////////////////////////////////////////
//
// dual(A) = I*A
//
// I_2d * 1 = e1^e2 * 1 = e1^e2
//
// I_2d * e1 = e1^e2 * e1 = e_121 = -e_112 = -e_2 = -e2
// I_2d * e2 = e1^e2 * e2 = e_122 =  e_122 =  e_1 =  e1
//
// I_2d * e1^e2 = e1^e2 * e1^e2 = e_1212 = -e_1122 = -1

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> dual2d(PScalar2d<T> ps)
{
    // dual(A) = I*A
    // e12 * (ps * e12) = -ps
    return Scalar2d<T>(-T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars
// the 2d and 3d adders in the function name are required for disambiguation
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual2d(Scalar2d<T> s)
{
    // dual(A) = I*A
    // e12 * (s) = s * e12
    return PScalar2d<T>(T(s));
}

// this overload is provided to accept T directly as an alternative to Scalar2d<T>
// e.g. with T as a result of a dot product between two vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual2d(T s)
{
    // dual(A) = I*A
    // e12 * (s) = s * e12
    return PScalar2d<T>(s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> dual2d(Vec2d<T> const& v)
{
    // dual(A) = I*A
    // e12 * (v.x * e1 + v.y * e2)
    //     =  v.y * e1 - v.x * e2
    return Vec2d<T>(v.y, -v.x);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> dual2d(MVec2d_E<T> const& M)
{
    // dual(A) = I*A
    // e12 * (  s + ps * e12)
    //     =  -ps +  s * e12
    return MVec2d_E<T>(-M.c1, M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> dual2d(MVec2d<T> const& M)
{
    // dual(A) = I*A
    // e12 * (  s + v.x * e1 + v.y * e2 + ps * e12)
    //     =  -ps + v.y * e1 - v.x * e2 + s * e12
    return MVec2d<T>(-M.c3, M.c2, -M.c1, M.c0);
}

#else
////////////////////////////////////////////////////////////////////////////////
// duality as defined in Macdonald, "Linear and geometric algebra", p. 109:
////////////////////////////////////////////////////////////////////////////////
//
// dual(A) = A/I = A*I^(-1) = A*rev(I)
//
// 1 * rev(I_2d) = 1 * e2^e1 = e2^e1 = e_21 = -e_12 = -e1^e2
//
// e1 * rev(I_2d) = e1 * e2^e1 = e_121 = -e_112 = -e_2 = -e2
// e2 * rev(I_2d) = e2 * e2^e1 = e_221 =  e_221 =  e_1 =  e1
//
// e1^e2 * rev(I_2d) = e1^e2 * e2^e1 = e_1221 = 1

// using this duality definition, following duality properties hold
// (A. Macdonald, "Linear and geometric algebra", p. 110):
//
// a) dual(aA) = a dual(A)
// b) dual(A + B) = dual(A) + dual(B)
// c) dual(dual(A)) = (-1)^(n*(n-1)/2) A   (with n as dimension of the (sub)space)
// d) |dual(B)| = |B|
// e) if B is a j-blade then dual(B) is an (n-j)-blade
// f) if A is a j-vector then dual(A) is an (n-j)-vector
//    (remember: a j-vector is a sum of j-blades, which are outer products)

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> dual2d(PScalar2d<T> ps)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (ps * e12) * e21
    // =  ps
    return Scalar2d<T>(T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars
// the 2d and 3d adders in the function name are required for disambiguation
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual2d(Scalar2d<T> s)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (s) * e21
    // =  -s * e12
    return PScalar2d<T>(-T(s));
}

// this overload is provided to accept T directly as an alternative to Scalar2d<T>
// e.g. with T as a result of a dot product between two vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual2d(T s)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (s) * e21
    // =  -s * e12
    return PScalar2d<T>(-s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> dual2d(Vec2d<T> const& v)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (v.x * e1 + v.y * e2) * e21
    // =  v.y * e1 - v.x * e2
    return Vec2d<T>(v.y, -v.x);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> dual2d(MVec2d_E<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (  s + ps * e12) * e21
    //   = ps -  s * e12
    return MVec2d_E<T>(M.c1, -M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> dual2d(MVec2d<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (  s + v.x * e1 + v.y * e2 + ps * e12) * e21
    //   = ps + v.y * e1 - v.x * e2 -  s * e12
    return MVec2d<T>(M.c3, M.c2, -M.c1, -M.c0);
}

#endif


////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> & PScalar2d<T> mixed operations using the geometric product
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> dot(PScalar2d<T> A, Vec2d<U> const& b)
{
    // the dot product is identical with the geometric product in this case (A^b = 0)
    // ATTENTION: the dot-product in NOT symmetric in G^n as it is in R^n
    return A * b;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> dot(Vec2d<T> const& a, PScalar2d<U> B)
{
    // the dot product is identical with the geometric product in this case (a^B = 0)
    // ATTENTION: the dot-product in NOT symmetric in G^n as it is in R^n
    return a * B;
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
    return dot(v1, v2) * inv(v2);
}

// projection of v1 onto v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>>
project_onto_normalized(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    // requires v2 to be normalized
    return dot(v1, v2) * v2;
}

// projection of v onto ps (returns the vector directly)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto(Vec2d<T> const& v,
                                                              PScalar2d<U> ps)
{
    return dot(v, ps) * inv(ps);
}

// rejection of v1 from v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from(Vec2d<T> const& v1,
                                                             Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    // version using geometric algebra wedge product manually computed
    // from "wdg(v1,v2)*inv(v2)"
    PScalar2d<ctype> w = wdg(v1, v2); // bivector with component e12
    ctype sq_n = nrm_sq(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    ctype w_sq_n_inv = w / sq_n;
    return Vec2d<ctype>(v2.y * w_sq_n_inv, -v2.x * w_sq_n_inv);
}

// rejection of v1 from v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>>
reject_from_normalized(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    // requires v2 to be normalized

    using ctype = std::common_type_t<T, U>;
    // version using geometric algebra wedge product manually computed
    // from "wdg(v1,v2)*inv(v2)" + v2 being already it's own inverse
    PScalar2d<ctype> w = wdg(v1, v2); // bivector with component e12
    return Vec2d<ctype>(v2.y * w, -v2.x * w);
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
std::vector<Vec2d<std::common_type_t<T, U>>> gs_orthogonal(Vec2d<T> const& u,
                                                           Vec2d<U> const& v)
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
std::vector<Vec2d<std::common_type_t<T, U>>> gs_orthonormal(Vec2d<T> const& u,
                                                            Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec2d<ctype>> basis;
    Vec2d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from_normalized(v, u_unitized)));
    return basis;
}

} // namespace hd::ga::ega