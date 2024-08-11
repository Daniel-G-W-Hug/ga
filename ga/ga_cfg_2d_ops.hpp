#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::clamp
#include <cmath>     // std::abs, std::sin, std::cos
#include <concepts>  // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "ga_cfg_value_t.hpp"

#include "ga_cfg_vec2d.hpp"

#include "ga_cfg_mvec2d.hpp"
#include "ga_cfg_mvec2d_e.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> & PScalar2d<T> mixed geometric operations
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

// projection of v1 onto v2 (v2 must already be unitized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto_unitized(Vec2d<T> const& v1,
                                                                       Vec2d<U> const& v2)
{
    // requires v2 to be unitized
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
    ctype sq_n = sq_nrm(v2);          //
    if (sq_n < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    ctype w_sq_n_inv = w / sq_n;
    return Vec2d<ctype>(v2.y * w_sq_n_inv, -v2.x * w_sq_n_inv);
}

// rejection of v1 from v2 (v2 must already be unitized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from_unitized(Vec2d<T> const& v1,
                                                                      Vec2d<U> const& v2)
{
    // requires v2 to be unitized

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
    return MVec2d_E<ctype>(Scalar<ctype>(dot(a, b)), wdg(a, b));
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

// geometric product A*B of two 2d pseudoscalars (=bivectors)
// bivector * bivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> operator*(PScalar2d<T> A, PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * ctype(B); // bivectors square to -1
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
    return MVec2d_E<ctype>(Scalar<ctype>(std::cos(theta)),
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
    return MVec2d_E<ctype>(Scalar<ctype>(std::cos(half_angle)),
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
inline constexpr Scalar<T> dual2d(PScalar2d<T> ps)
{
    // dual(A) = I*A
    // e12 * (ps * e12) = -ps
    return Scalar<T>(-T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars
// the 2d and 3d adders in the function name are required for disambiguation
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual2d(Scalar<T> s)
{
    // dual(A) = I*A
    // e12 * (s) = s * e12
    return PScalar2d<T>(T(s));
}

// this overload is provided to accept T directly as an alternative to Scalar<T>
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
inline constexpr Scalar<T> dual2d(PScalar2d<T> ps)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (ps * e12) * e21
    // =  ps
    return Scalar<T>(T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars
// the 2d and 3d adders in the function name are required for disambiguation
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual2d(Scalar<T> s)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (s) * e21
    // =  -s * e12
    return PScalar2d<T>(-T(s));
}

// this overload is provided to accept T directly as an alternative to Scalar<T>
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

} // namespace hd::ga