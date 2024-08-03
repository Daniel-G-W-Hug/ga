#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::clamp
#include <cmath>     // std::abs, std::sin, std::cos
#include <concepts>  // std::floating_point<T>
#include <iostream>
#include <limits>
#include <numbers> // math constants like pi
#include <stdexcept>
#include <string>

#include "ga_cfg_value_t.hpp"

#include "ga_cfg_vec3d.hpp"

#include "ga_cfg_bivec3d.hpp"

#include "ga_cfg_mvec3d.hpp"
#include "ga_cfg_mvec3d_e.hpp"
#include "ga_cfg_mvec3d_u.hpp"


namespace hd::ga {

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
inline std::common_type_t<T, U> angle(BiVec3d<T> const& v1, Vec3d<U> const& v2)
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

// projection of v1 onto v2 (v2 must already be unitized to nrm(v2) == 1)
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

// projection of a vector v1 onto a unitized bivector v2
// u_parallel = gr1(dot(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
project_onto_unitized(Vec3d<T> const& v1, BiVec3d<U> const& v2)
{
    // requires v2 to be unitized

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

// rejection of vector v1 from a unitized vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from_unitized(Vec3d<T> const& v1,
                                                                      Vec3d<U> const& v2)
{
    // requires v2 to be unitized
    using ctype = std::common_type_t<T, U>;
    BiVec3d<ctype> B = wdg(v1, v2);
    Vec3d<ctype> v2_inv = v2; // v2 is its own inverse, if unitized
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

// rejection of vector v1 from a unitized bivector v2
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
// MVec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 3d multivector
// Expensive! - Don't use if you don't have to! (64x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// a * b = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector vector = scalar + bivector)
// A * b = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector vector = vector + trivector)
// a * B = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)  (vector bivector = vector + trivector)
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
                           PScalar3d<ctype>(+A.x * B.c0 + A.y * B.c1 + A.z * B.c2));
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
    return -ctype(A) * ctype(B); // trivectors square to -1
}

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
//           (does not need to be unitized)
//         - a rotation angle
// Output:
//         - a rotor representing the rotation
//
// HINT:     For a rotation around an axis n (n = unitized(Vec3d<T>))
//           use the bivector B = n*I_3d  => B = Vec3d<T> * PScalar3d<T>
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> exp(BiVec3d<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar<ctype>(std::cos(theta)), unitized(I) * std::sin(theta));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//       1.) an arbitray bivector representing the oriented plane of rotation
//           (does not need to be unitized, defines what is a posive rotation angle)
//       2.) a rotation angle in that plane
// Output:
//           a rotor representing the requested rotation,
//           for applying the sandwich product as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// for a rotation about an axis n (n = unitized vector) choose the ansatz n*B = I_3d
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
                           unitized(I) * std::sin(half_angle));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rotate(Vec3d<T> const& v,
                                                        MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> reverse_rotor = rev(rotor);
    // MVec3d_U<ctype> tmp = rotor * v;
    // MVec3d_U<ctype> res = tmp * reverse rotor;

    // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the pseudoscalar part
    // which will be zero anyway
    return Vec3d<ctype>(gr1<ctype>(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rotate(BiVec3d<T> const& v,
                                                          MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> reverse_rotor = rev(rotor);
    // MVec3d_E<ctype> tmp = rotor * v;
    // MVec3d_E<ctype> res = tmp * reverse rotor;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway
    return BiVec3d<ctype>(gr2<ctype>(rotor * v * rev(rotor)));
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

// return the dual(M) of the multivector M
// if M represents the subspace B as subspace of R^2 then
// dual(M) represents the orthogonal subspace B^perp (perpendicular to B)
// => returns the orthogonal complement
//
// dual by left multiplication with Im_3d
// as defined in Doran/Lasenby "GA for physicists"
//
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar<T> dual3d(PScalar3d<T> ps)
{
    // e123 * ps * e123 = -ps
    return Scalar<T>(-T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform to different pseudoscalars
// the 2d and 3d adders in the function name are required for disambiguation
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual3d(Scalar<T> s)
{
    // e123 * s = s * e123
    return PScalar3d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> dual3d(Vec3d<T> const& v)
{
    // e123 * (v.x * e1  + v.y * e2  + v.z * e3)
    //      = (v.x * e23 + v.y * e31 + v.z * e12)
    return BiVec3d<T>(v.x, v.y, v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> dual3d(BiVec3d<T> const& B)
{
    // e123 * (  b.x * e23 + b.y * e31 + b.z * e12)
    //      = (- b.x * e1  - b.y * e2  - b.z * e3)
    return Vec3d<T>(-B.x, -B.y, -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> dual3d(MVec3d_E<T> const& M)
{
    // e123 * (s + b.x * e23 + b.y * e31 + b.z * e12)
    //      = (  - b.x * e1  - b.y * e2  - b.z * e3 + s * e123)
    return MVec3d_U<T>(-M.c1, -M.c2, -M.c3, M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> dual3d(MVec3d_U<T> const& M)
{
    // e123 * (      v.x * e1  + v.y * e2  + v.z * e3 + ps * e123)
    //      = (-ps + v.x * e23 + v.y * e31 + v.z * e12)
    return MVec3d_E<T>(-M.c3, M.c0, M.c1, M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> dual3d(MVec3d<T> const& M)
{
    // e123 * (  s + v.x * e1  + v.y * e2  + v.z * e3
    //             + b.x * e23 + b.y * e31 + b.z * e12 + ps * e123)
    //      = (-ps - b.x * e1  - b.y * e2  - b.z * e3
    //             + v.x * e23 + v.y * e31 + v.z * e12 +  s * e123)
    return MVec3d<T>(-M.c7, -M.c4, -M.c5, -M.c6, M.c1, M.c2, M.c3, M.c0);
}

} // namespace hd::ga
