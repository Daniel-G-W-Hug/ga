#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_pga3dp_ops_basics.hpp"
#include "ga_pga3dp_ops_products.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides functionality that is based on pga3dp ops basics and products:
//
// - angle()                        -> angle operations
// - TODO: exp()                          -> exponential function
// - motor(), motor_from_planes()   -> provide a motor
// - move3dp(), move3dp_opt()       -> move object with motor
// - project_onto(), reject_from()  -> simple projection and rejection
// - expand()                       -> expansion: new line/plane through point/line
//                                                perpendicular to line/plane
// - ortho_proj3dp()                -> orthogonal projection onto object
// - central_proj3dp()              -> central projection towards origin onto object
// - ortho_antiproj3dp()            -> orthogonal antiprojection onto object
// - reflect_on()                   -> reflections
// - support3dp()                   -> point on line that is nearest to origin
// - att                            -> object attitude
// - dist3dp()                      -> Euclidean distance and homogeneous magnitude
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// angle operations 3dp
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors, i.e. directions to points at infinity
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    if ((weight_nrm_sq(v1) != 0.0) || (weight_nrm_sq(v2) != 0.0)) {
        // the angle between points not at infinity or points not at infinity and a
        // direction towards infinity is defined as zero
        return 0.0;
    }

    // angle is defined only between directions towards points a infinity

    ctype nrm_prod = bulk_nrm(v1) * bulk_nrm(v2);
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between two bivectors, i.e. lines
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::common_type_t<T, U> angle(BiVec3dp<T> const& B1, BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = right_weight_contract3dp(B1, B2);
    // hint: weight_nrm returns pscalar! ctype() required around each single result,
    // otherwise geometric product which evaluates to zero
    ctype nrm_prod = ctype(weight_nrm(B1) * ctype(weight_nrm(B2)));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

// return the angle between a trivector and a bivector, i.e. a plane and a line
// range of angle: 0 <= angle <= pi/2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::common_type_t<T, U> angle(TriVec3dp<T> const& t, BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(bulk_nrm(right_weight_contract3dp(t, B)));
    // hint: weight_nrm returns pscalar! ctype() required around each single result,
    // otherwise geometric product which evaluates to zero
    ctype nrm_prod = ctype(weight_nrm(t)) * ctype(weight_nrm(B));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

// return the angle between a trivector and a bivector, i.e. a plane and a line
// range of angle: 0 <= angle <= pi/2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::common_type_t<T, U> angle(BiVec3dp<T> const& B, TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(bulk_nrm(right_weight_contract3dp(t, B)));
    // hint: weight_nrm returns pscalar! ctype() required around each single result,
    // otherwise geometric product which evaluates to zero
    ctype nrm_prod = ctype(weight_nrm(t)) * ctype(weight_nrm(B));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

// return the angle between two trivectors, i.e. two planes
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::common_type_t<T, U> angle(TriVec3dp<T> const& t1, TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(right_weight_contract3dp(t1, t2));
    // hint: weight_nrm returns pscalar! ctype() required around each single result,
    // otherwise geometric product which evaluates to zero
    ctype nrm_prod = ctype(weight_nrm(t1)) * ctype(weight_nrm(t2));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}


////////////////////////////////////////////////////////////////////////////////
// 3dp motor operations (translation and rotation)
//
// Every motor in pga3dp is an even grade multivector MVec3dp_E (w/o scalar part).
//
// A proper isometry in 3dp has a fixed line l around which a rotation occurs
// with an angle phi.
//
// The rotation so modeled by two consecutive reflections across two planes which
// intersect in the line l around which the rotation occurs.
//
// In this case the motor has the form: M = l sin(phi) + e1234 cos(phi).
//
// (derived from the exponential function with respect to the regressive
// geometric product)
////////////////////////////////////////////////////////////////////////////////

// create a (unitized) motor from a fixed line of rotation and a turning angle
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> motor(BiVec3dp<T> const& l, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = 0.5 * theta;
    return unitize(MVec3dp_E<ctype>(BiVec3dp<ctype>(l * std::sin(half_angle)),
                                    PScalar3dp<ctype>(std::cos(half_angle))));
}

// create a translation motor from a direction vector (given as a Vec3dp)
// move in direction and by length of direction vector (length = its bulk_nrm)
// ATTENTION: the direction is assumed to be a direction vector, i.e. with w == 0
//            the w-component is ignored and only the x-, y- and z-components are used
template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> motor(Vec3dp<T> const& direction)
{
    return MVec3dp_E<T>(
        0.5 * BiVec3dp<T>(T(0.0), T(0.0), T(0.0), direction.x, direction.y, direction.z),
        PScalar3dp<T>(1.0));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> motor_from_planes(TriVec3dp<T> const& t1,
                                                                TriVec3dp<U> const& t2)
{
    // take planes as input and return a motor R
    // 1st apply reflection across plane t1, then across t2 to get a motor that rotates
    // (or translates) around the intersection line of planes t1 and t2
    //
    // for use of motor R either directly on object u (inefficient):
    //     auto v_moved = gr1(rgpr(rgpr(R, v), rrev(R))));
    // or
    //     auto B_moved = gr2(rgpr(rgpr(R, B), rrev(R))));
    // or
    //     auto t_moved = gr3(rgpr(rgpr(R, t), rrev(R))));
    // or
    //                                   // optimized for reduced effort
    //     auto v_moved = move3dp(v,R);  // moves v according to the motor R
    //     auto B_moved = move3dp(B,R);  // moves B according to the motor R
    //     auto t_moved = move3dp(t,R);  // moves t according to the motor R
    //
    return unitize(rgpr(t2, t1)); // based on the regressive geometric product
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> move3dp(Vec3dp<T> const& v,
                                                   MVec3dp_E<U> const& R)
{
    // assumes: motor R is unitized

    // moves v (a vector representing a projective point) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(gr1(rgpr(rgpr(R, v), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> move3dp(BiVec3dp<T> const& B,
                                                     MVec3dp_E<U> const& R)
{
    // assumes: motor R is unitized

    // moves B (a bivector representing a line) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(gr2(rgpr(rgpr(R, B), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> move3dp(TriVec3dp<T> const& t,
                                                      MVec3dp_E<U> const& R)
{
    // assumes: motor R is unitized

    // moves t (a trivector representing a plane) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(gr3(rgpr(rgpr(R, t), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> move3dp_opt(Vec3dp<T> const& v,
                                                       MVec3dp_E<U> const& R)
{
    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    //
    // (could potentially be further optimized by exporting the matrix-representation if
    // many transformations should be calculated using the same rotor as v' = matrix * v)
    using ctype = std::common_type_t<T, U>;
    auto k01 = R.c0 * R.c1;
    auto k02 = R.c0 * R.c2;
    auto k03 = R.c0 * R.c3;
    auto k11 = R.c1 * R.c1;
    auto k12 = R.c1 * R.c2;
    auto k13 = R.c1 * R.c3;
    auto k15 = R.c1 * R.c5;
    auto k16 = R.c1 * R.c6;
    auto k17 = R.c1 * R.c7;
    auto k22 = R.c2 * R.c2;
    auto k23 = R.c2 * R.c3;
    auto k24 = R.c2 * R.c4;
    auto k26 = R.c2 * R.c6;
    auto k27 = R.c2 * R.c7;
    auto k33 = R.c3 * R.c3;
    auto k34 = R.c3 * R.c4;
    auto k35 = R.c3 * R.c5;
    auto k37 = R.c3 * R.c7;
    auto k47 = R.c4 * R.c7;
    auto k57 = R.c5 * R.c7;
    auto k67 = R.c6 * R.c7;
    auto k77 = R.c7 * R.c7;
    return Vec3dp<ctype>(
        (k11 - k22 - k33 + k77) * v.x + (2.0 * (k12 - k37)) * v.y +
            2.0 * (k13 + k27) * v.z + 2.0 * (-k01 + k26 - k35 + k47) * v.w,
        (2.0 * (k12 + k37)) * v.x + (-k11 + k22 - k33 + k77) * v.y +
            2.0 * (-k17 + k23) * v.z + 2.0 * (-k02 - k16 + k34 + k57) * v.w,
        (2.0 * (k13 - k27)) * v.x + (2.0 * (k17 + k23)) * v.y +
            (-k11 - k22 + k33 + k77) * v.z + 2.0 * (-k03 + k15 - k24 + k67) * v.w,
        (k11 + k22 + k33 + k77) * v.w);
}

// TODO: implement the optimized versions (optional, since -O3 already has same effect)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> move3dp_opt(BiVec3dp<T> const& B,
                                                         MVec3dp_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(gr2(rgpr(rgpr(R, B), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> move3dp_opt(TriVec3dp<T> const& t,
                                                          MVec3dp_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(gr3(rgpr(rgpr(R, t), rrev(R))));
}


// ////////////////////////////////////////////////////////////////////////////////
// // 3dp rotation operations
// ////////////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////////////////////
// // exponential function with bivector as argument for setup of quaternions
// // as geometric multivector with a scalar and a bivector part
// // MVec3dp_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
// //
// // quaternion: q = a + b I with I being the bivector in brackets above
// //             representing a plane in the algebra G^3
// //
// //             a rotation in 3D is represented by the plane and the
// //             size of the rotation, the later is given by the angle
// //             theta, which is the magnitude of the bivector
// //
// // Inputs:
// //         - an arbitray bivector representing the oriented plane of rotation
// //           (does not need to be normalized)
// //         - a rotation angle
// // Output:
// //         - a rotor representing the rotation
// //
// // HINT:     For a rotation around an axis n (n = normalize(Vec3dp<T>))
// //           use the bivector B = n*I_3d  => B = Vec3dp<T> * PScalar3dp<T>
// //////////////////////////////////////////////////////////////////////////////////////////
// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr MVec3dp_E<std::common_type_t<T, U>> exp(BiVec3dp<T> const& B, U theta)
// {
//     using ctype = std::common_type_t<T, U>;
//     return MVec3dp_E<ctype>(Scalar3dp<ctype>(std::cos(theta)),
//                             normalize(B) * std::sin(theta));
// }

// //////////////////////////////////////////////////////////////////////////////////////////
// // Inputs:
// //       1.) an arbitray bivector representing the oriented plane of rotation
// //           (does not need to be normalized, defines what is a posive rotation angle)
// //       2.) a rotation angle in that plane
// // Output:
// //           a rotor representing the requested rotation,
// //           for applying the sandwich product as in rotate(v,rotor)
// //
// //////////////////////////////////////////////////////////////////////////////////////////
// //
// // for a rotation about an axis n (n = normalized vector) choose the ansatz n*B = I_3d
// // and multiply both sides with n from the left (remember n*n = |n|^2 = 1)
// //
// // => choose: B = n*I_3d
// //
// //////////////////////////////////////////////////////////////////////////////////////////

// // TODO: add equivalent formulas from symbolic manipulations to replace the sandwich
// // product by a matrix * vector solution with minimal computational effort
// //
// // see ga_prdxpr

// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr MVec3dp_E<std::common_type_t<T, U>> rotor(BiVec3dp<T> const& B, U
// theta)
// {
//     using ctype = std::common_type_t<T, U>;
//     ctype half_angle = -0.5 * theta;
//     return MVec3dp_E<ctype>(Scalar3dp<ctype>(std::cos(half_angle)),
//                             normalize(B) * std::sin(half_angle));
// }

// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr Vec3dp<std::common_type_t<T, U>> rotate(Vec3dp<T> const& v,
//                                                          MVec3dp_E<U> const& rotor)
// {
//     using ctype = std::common_type_t<T, U>;

//     // MVec3dp_E<ctype> rr = rev(rotor);
//     // MVec3dp_U<ctype> tmp = rotor * v;
//     // MVec3dp_U<ctype> res = tmp * rr;

//     // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
//     //
//     // optimization potential for sandwich product by replacing the second product
//     // with a specific operation that skips the calculation of the pseudoscalar part
//     // which will be zero anyway
//     return Vec3dp<ctype>(gr1<ctype>(rotor * v * rev(rotor)));
// }

// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr Vec3dp<std::common_type_t<T, U>> rotate_opt(Vec3dp<T> const& v,
//                                                              MVec3dp_E<U> const& R)
// {
//     using ctype = std::common_type_t<T, U>;

//     ctype k1 = R.c0 * v.x - R.c5 * v.z + R.c6 * v.y;
//     ctype k2 = R.c0 * v.y + R.c4 * v.z - R.c6 * v.x;
//     ctype k3 = R.c0 * v.z - R.c4 * v.y + R.c5 * v.x;
//     ctype k4 = R.c0 * v.w + R.c1 * v.x + R.c2 * v.y + R.c3 * v.z;
//     ctype k5 = R.c2 * v.z - R.c3 * v.y + R.c4 * v.w - R.c7 * v.x;
//     ctype k6 = -R.c1 * v.z + R.c3 * v.x + R.c5 * v.w - R.c7 * v.y;
//     ctype k7 = R.c1 * v.y - R.c2 * v.x + R.c6 * v.w - R.c7 * v.z;
//     ctype k8 = -R.c4 * v.x - R.c5 * v.y - R.c6 * v.z;
//     return Vec3dp<ctype>(k1 * R.c0 + k2 * R.c6 - k3 * R.c5 - k8 * R.c4,
//                          -k1 * R.c6 + k2 * R.c0 + k3 * R.c4 - k8 * R.c5,
//                          k1 * R.c5 - k2 * R.c4 + k3 * R.c0 - k8 * R.c6,
//                          k1 * R.c1 + k2 * R.c2 + k3 * R.c3 + k4 * R.c0 + k5 * R.c4 +
//                              k6 * R.c5 + k7 * R.c6 + k8 * R.c7);
// }

// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr BiVec3dp<std::common_type_t<T, U>> rotate(BiVec3dp<T> const& B,
//                                                            MVec3dp_E<U> const& rotor)
// {
//     using ctype = std::common_type_t<T, U>;

//     // MVec3dp_E<ctype> rr = rev(rotor);
//     // MVec3dp_E<ctype> tmp = rotor * B;
//     // MVec3dp_E<ctype> res = tmp * rr;

//     // scalar and pseudoscalar part of res are 0 due to symmetric product
//     // rotor * B * rev(rotor)
//     //
//     // optimization potential for sandwich product by replacing the second product
//     // with a specific operation that skips the calculation of the scalar part
//     // which will be zero anyway
//     return BiVec3dp<ctype>(gr2<ctype>(rotor * B * rev(rotor)));
// }

// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr BiVec3dp<std::common_type_t<T, U>> rotate_opt(BiVec3dp<T> const& B,
//                                                                MVec3dp_E<U> const& R)
// {
//     using ctype = std::common_type_t<T, U>;

//     ctype k1 = -R.c4 * B.mx - R.c5 * B.my - R.c6 * B.mz;
//     ctype k2 =
//         R.c0 * B.vx - R.c2 * B.mz + R.c3 * B.my - R.c5 * B.vz + R.c6 * B.vy + R.c7 *
//         B.mx;
//     ctype k3 =
//         R.c0 * B.vy + R.c1 * B.mz - R.c3 * B.mx + R.c4 * B.vz - R.c6 * B.vx + R.c7 *
//         B.my;
//     ctype k4 =
//         R.c0 * B.vz - R.c1 * B.my + R.c2 * B.mx - R.c4 * B.vy + R.c5 * B.vx + R.c7 *
//         B.mz;
//     ctype k5 = R.c0 * B.mx - R.c5 * B.mz + R.c6 * B.my;
//     ctype k6 = R.c0 * B.my + R.c4 * B.mz - R.c6 * B.mx;
//     ctype k7 = R.c0 * B.mz - R.c4 * B.my + R.c5 * B.mx;
//     ctype k8 = -R.c1 * B.mx - R.c2 * B.my - R.c3 * B.mz - R.c4 * B.vx - R.c5 * B.vy -
//                R.c6 * B.vz;

//     return BiVec3dp<ctype>(k1 * R.c1 + k2 * R.c0 + k3 * R.c6 - k4 * R.c5 + k5 * R.c7 +
//                                k6 * R.c3 - k7 * R.c2 - k8 * R.c4,
//                            k1 * R.c2 - k2 * R.c6 + k3 * R.c0 + k4 * R.c4 - k5 * R.c3 +
//                                k6 * R.c7 + k7 * R.c1 - k8 * R.c5,
//                            k1 * R.c3 + k2 * R.c5 - k3 * R.c4 + k4 * R.c0 + k5 * R.c2 -
//                                k6 * R.c1 + k7 * R.c7 - k8 * R.c6,
//                            -k1 * R.c4 + k5 * R.c0 + k6 * R.c6 - k7 * R.c5,
//                            k1 * R.c5 - k5 * R.c6 + k6 * R.c0 + k7 * R.c4,
//                            k1 * R.c6 + k5 * R.c5 - k6 * R.c4 + k7 * R.c0);
// }

// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr TriVec3dp<std::common_type_t<T, U>> rotate(TriVec3dp<T> const& t,
//                                                             MVec3dp_E<U> const& rotor)
// {
//     using ctype = std::common_type_t<T, U>;

//     // MVec3dp_E<ctype> rr = rev(rotor);
//     // MVec3dp_E<ctype> tmp = rotor * t;
//     // MVec3dp_E<ctype> res = tmp * rr;

//     // vector part of res is 0 due to symmetric product  rotor * t * rev(rotor)
//     //
//     // optimization potential for sandwich product by replacing the second product
//     // with a specific operation that skips the calculation of the vector part
//     // which will be zero anyway
//     return TriVec3dp<ctype>(gr3<ctype>(rotor * t * rev(rotor)));
// }

// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr TriVec3dp<std::common_type_t<T, U>> rotate_opt(TriVec3dp<T> const& t,
//                                                                 MVec3dp_E<U> const& R)
// {
//     using ctype = std::common_type_t<T, U>;

//     ctype k1 = R.c4 * t.w;
//     ctype k2 = R.c5 * t.w;
//     ctype k3 = R.c6 * t.w;
//     ctype k4 = -R.c4 * t.x - R.c5 * t.y - R.c6 * t.z - R.c7 * t.w;
//     ctype k5 = R.c0 * t.x - R.c1 * t.w - R.c5 * t.z + R.c6 * t.y;
//     ctype k6 = R.c0 * t.y - R.c2 * t.w + R.c4 * t.z - R.c6 * t.x;
//     ctype k7 = R.c0 * t.z - R.c3 * t.w - R.c4 * t.y + R.c5 * t.x;
//     ctype k8 = R.c0 * t.w;

//     return TriVec3dp<ctype>(k1 * R.c0 + k2 * R.c6 - k3 * R.c5 - k8 * R.c4,
//                             -k1 * R.c6 + k2 * R.c0 + k3 * R.c4 - k8 * R.c5,
//                             k1 * R.c5 - k2 * R.c4 + k3 * R.c0 - k8 * R.c6,
//                             k1 * R.c1 + k2 * R.c2 + k3 * R.c3 + k4 * R.c0 + k5 * R.c4 +
//                                 k6 * R.c5 + k7 * R.c6 + k8 * R.c7);
// }


////////////////////////////////////////////////////////////////////////////////
// projections, rejections
////////////////////////////////////////////////////////////////////////////////

// TODO: check whether the vector-vector formulas make sense at all, since they model
//       the representational space an not the modelled subspace

// projection of a vector v1 onto vector v2
// returns component of v1 parallel to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> project_onto(Vec3dp<T> const& v1,
                                                        Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of vector v1 from a vector v2
// returns component of v1 perpendicular to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> reject_from(Vec3dp<T> const& v1,
                                                       Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec3dp<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}


// projection of a vector v onto a bivector B (a line)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> project_onto(Vec3dp<T> const& v,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(rwdg(B, wdg(v, right_weight_dual(B)))); // ortho_proj3dp
}

// rejection of vector v from a bivector B (a line)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> reject_from(Vec3dp<T> const& v,
                                                       BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(v - project_onto(v, B));
}


// projection of a vector v onto a trivector t (a plane)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> project_onto(Vec3dp<T> const& v,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(rwdg(t, wdg(v, right_weight_dual(t)))); // ortho_proj3dp
}

// rejection of vector v from a bivector B (a line)
// u_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> reject_from(Vec3dp<T> const& v,
                                                       TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(v - project_onto(v, t));
}


// expand to a new line with goes through point p and is perpendicular to the plane
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Line3d<std::common_type_t<T, U>> expand(Point3d<T> const& point,
                                                  Plane3d<U> const& plane)
{
    return weight_expansion(point3dp{point}, plane);
}

// expand to a new plane with contains the point p and is orthogonal to the line
// => returns a plane (aka a trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Plane3d<std::common_type_t<T, U>> expand(Point3d<T> const& point,
                                                   Line3d<U> const& line)
{
    return weight_expansion(point3dp{point}, line);
}

// expand to a new plane which contains the line and is perpendicular to the plane
// => returns a plane (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Plane3d<std::common_type_t<T, U>> expand(Line3d<T> const& line,
                                                   Plane3d<U> const& plane)
{
    return weight_expansion(line, plane);
}


////////////////////////////////////////////////////////////////////////////////
// Projections for 3dp:
//
// ortho_proj3dp(a, b)     = rwdg(b, right_weight_expand3dp(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' contained in b)
// REQUIRES: gr(a) < gr(b)
//
//
// central_proj3dp(a, b)   = rwdg(b, right_bulk_expand3dp(a, b) )
// (a projected centrally towards origin onto b, effectively creating a new a'
// contained in b)
// REQUIRES: gr(a) < gr(b)
//
// ortho_antiproj3dp(a, b) = wdg(b, right_weight_contract3dp(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' containing b)
// REQUIRES: gr(a) > gr(b)
//
////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2> decltype(auto) ortho_proj3dp(arg1&& a, arg2&& b)
{

    // REQUIRES: gr(a) < gr(b), or does not compile!
    // project the smaller grade object onto to larger grade object
    return rwdg(std::forward<arg2>(b),
                right_weight_expand3dp(std::forward<arg1>(a), std::forward<arg2>(b)));
}

template <typename arg1, typename arg2> decltype(auto) central_proj3dp(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!
    // project the smaller grade object onto to larger grade object
    return rwdg(std::forward<arg2>(b),
                right_bulk_expand3dp(std::forward<arg1>(a), std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) ortho_antiproj3dp(arg1&& a, arg2&& b)
{
    return wdg(std::forward<arg2>(b),
               right_weight_contract3dp(std::forward<arg1>(a), std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// reflections of points, lines, and planes on planes
////////////////////////////////////////////////////////////////////////////////

// reflect a vector u in an arbitrary trivector, i.e. a plane t
// t must be unitized
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> reflect_on(Vec3dp<T> const& v,
                                                      TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(-gr1(rgpr(rgpr(t, v), t)));
}

// reflect a bivector B (a line) in an arbitrary trivector t
// t must be unitized
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> reflect_on(BiVec3dp<T> const& B,
                                                        TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-gr2(rgpr(rgpr(t, B), t)));
}

// reflect a trivector t1 (a plane) in an arbitrary trivector t2 (a unitized plane)
// t2 must be unitized
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> reflect_on(TriVec3dp<T> const& t1,
                                                         TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(-gr3(rgpr(rgpr(t2, t1), t2)));
}


////////////////////////////////////////////////////////////////////////////////
// support3dp: nearest point to the origin (by projecting origin onto line/plane)
////////////////////////////////////////////////////////////////////////////////

template <typename arg1> decltype(auto) support3dp(arg1&& a)
{
    // REQUIRES: a line (BiVec3dp) or a plane (TriVec3dp) as argument
    return ortho_proj3dp(origin_3dp, std::forward<arg1>(a));
}


////////////////////////////////////////////////////////////////////////////////
// attitude operations: att = rwdg( u, rcmpl(e4_3dp) ) = rwdg(u, horizon_3dp)
//
// (the attitude is the intersection of the object with the horizon)
// the result of att(object_with_grade_k) is an object with grade k-1
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> att(Vec3dp<T> const& v)
{
    return Scalar3dp<T>(v.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> att(Point3dp<T> const& p)
{
    return Scalar3dp<T>(p.w);
}

// return the attitude (i.e. the direction vector) of the line
template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> att(BiVec3dp<T> const& B)
{
    return Vec3dp<T>(B.vx, B.vy, B.vz, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> att(Line3d<T> const& l)
{
    return Vec3dp<T>(l.vx, l.vy, l.vz, T(0.0));
}

// return the attitude (i.e. the normal vector) of the plane
template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> att(TriVec3dp<T> const& t)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), t.x, t.y, t.z);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> att(Plane3d<T> const& l)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), l.x, l.y, l.z);
}


template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> att(PScalar3dp<T> ps)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), ps);
}


////////////////////////////////////////////////////////////////////////////////
// 3dp euclidean distance
////////////////////////////////////////////////////////////////////////////////

// returns the euclidean distance between objects as homogeneous magnitude
template <typename arg1, typename arg2>
constexpr DualNum3dp<value_t> dist3dp(arg1&& a, arg2&& b)
{
    if (gr(a) + gr(b) == 4) {
        return DualNum3dp<value_t>(rwdg(a, b), weight_nrm(wdg(a, att(b))));
    }
    else {
        return DualNum3dp<value_t>(bulk_nrm(att(wdg(a, b))), weight_nrm(wdg(a, att(b))));
    }
}

} // namespace hd::ga::pga