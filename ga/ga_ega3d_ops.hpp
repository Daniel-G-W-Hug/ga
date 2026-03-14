#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_ega3d_ops_basics.hpp"   // ega3d ops basics
#include "ga_ega3d_ops_products.hpp" // ega3d ops products

#include <vector>


namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega3d functionality that is based on ega3d ops basics and products:
//
// - angle()                        -> angle operations
// - exp(bivec) -> rotor            -> exponential function (w.r.t. gpr)
// - sqrt(rotor) -> rotor           -> sqrt function (w.r.t. gpr) halves the rot. angle
// - get_rotor()                    -> provide a rotor
// - rotate(), rotate_opt()         -> rotate object with rotor (sandwich + optimized)
// - project_onto(), reject_from()  -> projection and rejection
// - reflect_on(), reflect_on_vec() -> reflections
// - gs_orthogonal()                -> Gram-Schmidt-orthogonalization
// - is_congruent()                 -> Same up to a scalar factor (is same subspace)
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// angle operations
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    ctype nrm_prod = nrm(v1) * nrm(v2);
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");

    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


// unsuccessful try to extend angle range to -pi ... pi,
// because orientation is not defined uniquely:
//
// // return the angle between of two vectors
// // range of angle: -pi <= angle <= pi
// template <typename T, typename U>
//     requires(numeric_type<T> && numeric_type<U>)
// inline constexpr std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U>
// const& v2)
// {
//     using ctype = std::common_type_t<T, U>;
//     using std::numbers::pi;

// ctype nrm_prod = nrm(v1) * nrm(v2);
// hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
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

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> angle(BiVec3d<T> const& v1, BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v1) * nrm(v2);
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


// return the angle between of a vector and a bivector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v) * nrm(B);
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
    // std::clamp must be used to take care of numerical inaccuracies
    // return std::acos(
    //     std::clamp(ctype(nrm(dot(v, B))) / nrm_prod, ctype(-1.0), ctype(1.0)));
    return std::acos(
        std::clamp(ctype(nrm((v << B))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between of a bivector and a vector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> angle(BiVec3d<T> const& B, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(B) * nrm(v);
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");

    // std::clamp must be used to take care of numerical inaccuracies
    // return std::acos(
    //     std::clamp(ctype(nrm(dot(B, v))) / nrm_prod, ctype(-1.0), ctype(1.0)));
    return std::acos(
        std::clamp(ctype(nrm((B >> v))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


////////////////////////////////////////////////////////////////////////////////
// 3d exponential function w.r.t. the geometric product
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3d_E<T> exp(BiVec3d<T> const& B)
{
    // exp(B) = cos(phi) + sin(phi) B/nrm(B) = cos(phi) + sin(phi) B_hat
    // with
    // B = phi * B_hat, where B_hat = B/B_nrm and B_hat^2 = -1
    // (B_hat is the bivector representing a plane, a 2D subspace of 3D space)
    //
    // ATTENTION: B = -phi * B_hat = phi * (-B_hat) are equivalent and indistinguishable

    T const B_nrm = nrm(B);
    if (B_nrm <= detail::safe_epsilon<T>()) {
        return MVec3d_E<T>(Scalar3d<T>(1.0),
                           BiVec3d<T>(0.0, 0.0, 0.0)); // return identity transformation
    }
    auto B_hat{B / B_nrm};

    auto phi = B_nrm;
    return MVec3d_E<T>(Scalar3d<T>(std::cos(phi)), B_hat * std::sin(phi));
}


////////////////////////////////////////////////////////////////////////////////
// sqrt(motor) w.r.t. geometric product
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3d_E<T> sqrt(MVec3d_E<T> const& R)
{
    auto nrmsq = nrm_sq(R);
    auto M{R};
    if (nrmsq != 1.0) {
        M = normalize(M); // motor must be normalized!
    }
    if (M.c0 == -1.0 && M.c1 == 0.0 && M.c2 == 0.0 && M.c3 == 0.0) { // rotation by 2*pi
        M.c0 = 1.0; // replace by identity transformation
    }

    // when M is a unit rotor, i.e. |M|=1, we have
    // M      = exp(-phi*e12_2d) = cos(phi) - e12_2d*sin(phi)
    // 1 + M  = (1+cos(phi)) - e12_2d * sin(phi)
    // Imagine this as adding two unit vectors. The sum points at direction -phi/2.
    // The normalization gives us the rotation with half the angle.
    // normalize(1 + M) = cos(phi/2) - e12_2d*sin(phi/2) = exp(-phi/2*e12_2d) = sqrt(M)

    return normalize(Scalar3d<T>(1.0) + M);
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
//           (does not need to be normalized)
//         - a rotation angle
// Output:
//         - a rotor representing the rotation
//
// HINT:     For a rotation around an axis n (n = normalize(Vec3d<T>))
//           use the bivector B = n*I_3d  => B = Vec3d<T> * PScalar3d<T>
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> exp(BiVec3d<T> const& B, U theta)
{
    using ctype = std::common_type_t<T, U>;
    if ((nrm(B) <= detail::safe_epsilon<T>()) ||
        (std::abs(theta) <= detail::safe_epsilon<T>())) {
        return MVec3d_E<T>(Scalar3d<T>(1.0),
                           BiVec3d<T>(0.0, 0.0, 0.0)); // return identity transformation
    }

    // B must be normalized or otherwise it will not only be scaled with angle theta
    return MVec3d_E<ctype>(Scalar3d<ctype>(std::cos(theta)),
                           normalize(B) * std::sin(theta));
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
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d_E<std::common_type_t<T, U>> get_rotor(BiVec3d<T> const& B, U theta)
{
    // attention: Bivector B which describes the plane of rotation is normalized
    //            to avoid unexpected surprises
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec3d_E<ctype>(Scalar3d<ctype>(std::cos(half_angle)),
                           normalize(B) * std::sin(half_angle));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rotate(Vec3d<T> const& v,
                                                 MVec3d_E<U> const& rotor)
{
    // rotate one vector
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> rr = rev(rotor);
    // MVec3d_U<ctype> tmp = rotor * v;
    // MVec3d_U<ctype> res = tmp * rr;

    // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the pseudoscalar part
    // which will be zero anyway
    // (typically done automatically with a good optimizer in -O3)
    return Vec3d<ctype>(gr1(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> rotate_opt(Vec3d<T> const& v,
                                                     MVec3d_E<U> const& R)
{
    // rotate one vector with rotor R (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (ega3d sandwich product)
    ctype const h0 = R.c0 * R.c0;
    ctype const h1 = R.c1 * R.c1;
    ctype const h2 = R.c2 * R.c2;
    ctype const h3 = R.c3 * R.c3;
    ctype const h4 = R.c0 * R.c3;
    ctype const h5 = R.c1 * R.c2;
    ctype const h6 = R.c0 * R.c2;
    ctype const h7 = R.c1 * R.c3;
    ctype const h8 = R.c0 * R.c1;
    ctype const h9 = R.c2 * R.c3;

    ctype const k11 = h0 + h1 - h2 - h3;
    ctype const k12 = 2.0 * (h4 + h5);
    ctype const k13 = 2.0 * (-h6 + h7);
    ctype const k21 = 2.0 * (-h4 + h5);
    ctype const k22 = h0 - h1 + h2 - h3;
    ctype const k23 = 2.0 * (h8 + h9);
    ctype const k31 = 2.0 * (h6 + h7);
    ctype const k32 = 2.0 * (-h8 + h9);
    ctype const k33 = h0 - h1 - h2 + h3;

    return Vec3d<ctype>(k11 * v.x + k12 * v.y + k13 * v.z,
                        k21 * v.x + k22 * v.y + k23 * v.z,
                        k31 * v.x + k32 * v.y + k33 * v.z);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::vector<Vec3d<std::common_type_t<T, U>>>
rotate_opt(std::vector<Vec3d<T>> const& vec, MVec3d_E<U> const& R)
{
    // rotate one vector with rotor R (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (ega3d sandwich product)
    ctype const h0 = R.c0 * R.c0;
    ctype const h1 = R.c1 * R.c1;
    ctype const h2 = R.c2 * R.c2;
    ctype const h3 = R.c3 * R.c3;
    ctype const h4 = R.c0 * R.c3;
    ctype const h5 = R.c1 * R.c2;
    ctype const h6 = R.c0 * R.c2;
    ctype const h7 = R.c1 * R.c3;
    ctype const h8 = R.c0 * R.c1;
    ctype const h9 = R.c2 * R.c3;

    ctype const k11 = h0 + h1 - h2 - h3;
    ctype const k12 = 2.0 * (h4 + h5);
    ctype const k13 = 2.0 * (-h6 + h7);
    ctype const k21 = 2.0 * (-h4 + h5);
    ctype const k22 = h0 - h1 + h2 - h3;
    ctype const k23 = 2.0 * (h8 + h9);
    ctype const k31 = 2.0 * (h6 + h7);
    ctype const k32 = 2.0 * (-h8 + h9);
    ctype const k33 = h0 - h1 - h2 + h3;

    std::vector<Vec3d<ctype>> result;
    result.reserve(vec.size());

    for (auto const& v : vec) {
        result.emplace_back(Vec3d<ctype>(k11 * v.x + k12 * v.y + k13 * v.z,
                                         k21 * v.x + k22 * v.y + k23 * v.z,
                                         k31 * v.x + k32 * v.y + k33 * v.z));
    }
    return result;
}


template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> rotate(BiVec3d<T> const& B,
                                                   MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway

    return BiVec3d<ctype>(gr2(rotor * B * rev(rotor)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> rotate_opt(BiVec3d<T> const& B,
                                                       MVec3d_E<U> const& R)
{
    // rotate one bivector B with rotor R (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (ega3d sandwich product)
    // coefficients kij are identical to coefficients used for vector transformation
    ctype const h0 = R.c0 * R.c0;
    ctype const h1 = R.c1 * R.c1;
    ctype const h2 = R.c2 * R.c2;
    ctype const h3 = R.c3 * R.c3;
    ctype const h4 = R.c0 * R.c3;
    ctype const h5 = R.c1 * R.c2;
    ctype const h6 = R.c0 * R.c2;
    ctype const h7 = R.c1 * R.c3;
    ctype const h8 = R.c0 * R.c1;
    ctype const h9 = R.c2 * R.c3;

    ctype const k11 = h0 + h1 - h2 - h3;
    ctype const k12 = 2.0 * (h4 + h5);
    ctype const k13 = 2.0 * (-h6 + h7);
    ctype const k21 = 2.0 * (-h4 + h5);
    ctype const k22 = h0 - h1 + h2 - h3;
    ctype const k23 = 2.0 * (h8 + h9);
    ctype const k31 = 2.0 * (h6 + h7);
    ctype const k32 = 2.0 * (-h8 + h9);
    ctype const k33 = h0 - h1 - h2 + h3;

    return BiVec3d<ctype>(k11 * B.x + k12 * B.y + k13 * B.z,
                          k21 * B.x + k22 * B.y + k23 * B.z,
                          k31 * B.x + k32 * B.y + k33 * B.z);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::vector<BiVec3d<std::common_type_t<T, U>>>
rotate_opt(std::vector<BiVec3d<T>> const& bvec, MVec3d_E<U> const& R)
{
    // rotate many bivectors bivec(B) with rotor R (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (ega3d sandwich product)
    // coefficients kij are identical to coefficients used for vector transformation
    ctype const h0 = R.c0 * R.c0;
    ctype const h1 = R.c1 * R.c1;
    ctype const h2 = R.c2 * R.c2;
    ctype const h3 = R.c3 * R.c3;
    ctype const h4 = R.c0 * R.c3;
    ctype const h5 = R.c1 * R.c2;
    ctype const h6 = R.c0 * R.c2;
    ctype const h7 = R.c1 * R.c3;
    ctype const h8 = R.c0 * R.c1;
    ctype const h9 = R.c2 * R.c3;

    ctype const k11 = h0 + h1 - h2 - h3;
    ctype const k12 = 2.0 * (h4 + h5);
    ctype const k13 = 2.0 * (-h6 + h7);
    ctype const k21 = 2.0 * (-h4 + h5);
    ctype const k22 = h0 - h1 + h2 - h3;
    ctype const k23 = 2.0 * (h8 + h9);
    ctype const k31 = 2.0 * (h6 + h7);
    ctype const k32 = 2.0 * (-h8 + h9);
    ctype const k33 = h0 - h1 - h2 + h3;

    std::vector<BiVec3d<ctype>> result;
    result.reserve(bvec.size());

    for (auto const& B : bvec) {
        result.emplace_back(BiVec3d<ctype>(k11 * B.x + k12 * B.y + k13 * B.z,
                                           k21 * B.x + k22 * B.y + k23 * B.z,
                                           k31 * B.x + k32 * B.y + k33 * B.z));
    }
    return result;
}


template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3d<std::common_type_t<T, U>> rotate(MVec3d<T> const& M,
                                                  MVec3d_E<U> const& rotor)
{
    // rotate one multivector M with rotor R
    // (only rotates the vector and bivector components of the 3d multivector)
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(rotor * M * rev(rotor));
}


////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> and BiVec3d<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
// v_parallel = dot(v1, v2)) * inv(v2)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v1,
                                                       Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v1,
                                                      Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec3d<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}

// projection of a vector v onto a bivector B
// v_parallel = gr1((B >> v) * inv(B))
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v,
                                                       BiVec3d<U> const& B)
{
    return gr1((B >> v) * inv(B));
}

// rejection of vector v1 from a bivector v2
// u_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v,
                                                      BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v - project_onto(v, B));

    // works, but is more effort compared to solution via projection and vector difference
    // return gr1(wdg(v, B) * inv(B));
}

////////////////////////////////////////////////////////////////////////////////
// reflections
////////////////////////////////////////////////////////////////////////////////
// Macdonald p. 129:
//
// Reflect a j-blade u in the k-dimensional subspace B is the blade
//
//   u_reflected = (-1)^[j*(k+1)]*B*u*inv(B)
//
////////////////////////////////////////////////////////////////////////////////

// reflect a vector v on a hyperplane B orthogonal to vector nB
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n (a plane in 3d space)
// orthogonal to vector nB: the hyperplane B is dual to nB (i.e. a two dimensional
// subspace)
//
// HINT: choose nB = dual(B) (for normalized nB)
//
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> reflect_on(Vec3d<T> const& v,
                                                     Vec3d<U> const& nB)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(-nB * v * inv(nB)));
}

// reflect a vector v in an arbitrary bivector B, i.e. a plane
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> reflect_on(Vec3d<T> const& v,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(-B * v * inv(B)));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling planes)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3d<std::common_type_t<T, U>> reflect_on(BiVec3d<T> const& UB,
                                                       BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(gr2(B * UB * inv(B)));
}

// reflect a vector v on another vector b
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3d<std::common_type_t<T, U>> reflect_on_vec(Vec3d<T> const& v,
                                                         Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(b * v * inv(b)));
}

////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt-Orthogonalization part 1: 2d plane embedded in 3d space
////////////////////////////////////////////////////////////////////////////////
//
// input:  two linear independent vectors u and v in 3d defining a plane
// output: two orthogonal vectors with the first one being u and the second one a vector
// perpendicular to u in the orientation of v, both forming an orthogonal system
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
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
    requires(numeric_type<T> && numeric_type<U>)
std::vector<Vec3d<std::common_type_t<T, U>>> gs_orthonormal(Vec3d<T> const& u,
                                                            Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec3d<ctype>> basis;
    Vec3d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from(v, u_unitized)));
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
    requires(numeric_type<U> && numeric_type<V> && numeric_type<W>)
constexpr std::vector<Vec3d<std::common_type_t<U, V, W>>>
gs_orthogonal(Vec3d<U> const& u, Vec3d<V> const& v, Vec3d<W> const& w)
{
    using ctype = std::common_type_t<U, V, W>;
    std::vector<Vec3d<ctype>> basis;
    basis.reserve(3);
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
    requires(numeric_type<U> && numeric_type<V> && numeric_type<W>)
constexpr std::vector<Vec3d<std::common_type_t<U, V, W>>>
gs_orthonormal(Vec3d<U> const& u, Vec3d<V> const& v, Vec3d<W> const& w)
{
    using ctype = std::common_type_t<U, V, W>;
    std::vector<Vec3d<ctype>> basis;
    basis.reserve(3);
    Vec3d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from(v, u_unitized)));
    basis.emplace_back(normalize(reject_from(w, wdg(u, v))));
    return basis;
}

////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar3d<T> a, Scalar3d<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero scalars are congruent (represent the same 0-dimensional subspace)
    return true;
}

// For vectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Vec3d<T> const& a, Vec3d<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.x) < tolerance) && (std::abs(a.y) < tolerance) &&
                     (std::abs(a.z) < tolerance);
    bool b_is_zero = (std::abs(b.x) < tolerance) && (std::abs(b.y) < tolerance) &&
                     (std::abs(b.z) < tolerance);

    if (a_is_zero && b_is_zero) {
        return true; // Both are effectively zero
    }
    if (a_is_zero || b_is_zero) {
        return false; // Only one is zero
    }

    // Find scale factor k where a = k*b, checking all components
    ctype k = 0.0;
    bool k_found = false;

    // Find first non-zero component pair to establish k
    if (std::abs(b.x) > tolerance) {
        k = a.x / b.x;
        k_found = true;
    }
    else if (std::abs(b.y) > tolerance) {
        k = a.y / b.y;
        k_found = true;
    }
    else if (std::abs(b.z) > tolerance) {
        k = a.z / b.z;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol =
        tolerance * std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z), std::abs(b.x),
                              std::abs(b.y), std::abs(b.z), value_t(1.0)});
    return (std::abs(a.x - k * b.x) < rel_tol) && (std::abs(a.y - k * b.y) < rel_tol) &&
           (std::abs(a.z - k * b.z) < rel_tol);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec3d<T> const& a, BiVec3d<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.x) < tolerance) && (std::abs(a.y) < tolerance) &&
                     (std::abs(a.z) < tolerance);
    bool b_is_zero = (std::abs(b.x) < tolerance) && (std::abs(b.y) < tolerance) &&
                     (std::abs(b.z) < tolerance);

    if (a_is_zero && b_is_zero) {
        return true; // Both are effectively zero
    }
    if (a_is_zero || b_is_zero) {
        return false; // Only one is zero
    }

    // Find scale factor k where a = k*b, checking all components
    ctype k = 0.0;
    bool k_found = false;

    // Find first non-zero component pair to establish k
    if (std::abs(b.x) > tolerance) {
        k = a.x / b.x;
        k_found = true;
    }
    else if (std::abs(b.y) > tolerance) {
        k = a.y / b.y;
        k_found = true;
    }
    else if (std::abs(b.z) > tolerance) {
        k = a.z / b.z;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol =
        tolerance * std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z), std::abs(b.x),
                              std::abs(b.y), std::abs(b.z), value_t(1.0)});
    return (std::abs(a.x - k * b.x) < rel_tol) && (std::abs(a.y - k * b.y) < rel_tol) &&
           (std::abs(a.z - k * b.z) < rel_tol);
}

// For pseudoscalars: all non-zero pseudoscalars in 3D represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar3d<T> a, PScalar3d<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars in 3D are congruent (represent the full 3D space)
    return true;
}

} // namespace hd::ga::ega