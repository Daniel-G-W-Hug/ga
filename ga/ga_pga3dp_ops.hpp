#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_pga3dp_ops_basics.hpp"
#include "ga_pga3dp_ops_products.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides functionality that is based on pga3dp ops basics and products:
//
// - angle()                              -> angle operations
// - TODO: exp()                          -> exponential function
// - get_motor(), get_motor_from_planes() -> provide a motor
// - move3dp(), move3dp_opt()             -> move object with motor
// - project_onto(), reject_from()        -> simple projection and rejection
// - expand()                             -> expansion: new line/plane through point/line
//                                                      perpendicular to line/plane
// - ortho_proj3dp()                      -> orthogonal projection onto object
// - central_proj3dp()                    -> central projection towards origin onto object
// - ortho_antiproj3dp()                  -> orthogonal antiprojection onto object
// - reflect_on()                         -> reflections
// - invert_on()                          -> inversions
// - support()                            -> point on line/plane that is nearest to origin
// - att()                                -> object attitude
// - dist3dp()                            -> Euclidean distance and homogeneous magnitude
// - is_congruent()                       -> Same up to a scalar factor (is same subspace)
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

    if ((ctype(weight_nrm_sq(v1)) != 0.0) || (ctype(weight_nrm_sq(v2) != 0.0))) {
        // the angle between points not at infinity or points not at infinity and a
        // direction towards infinity is defined as zero
        return 0.0;
    }

    // angle is defined only between directions towards points a infinity

    ctype nrm_prod = ctype(bulk_nrm(v1)) * ctype(bulk_nrm(v2));
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
    ctype contr = ctype(right_weight_contract3dp(B1, B2));
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
// with an angle phi. The line is represented by a bivector.
//
// The rotation so modelled by two consecutive reflections across two planes which
// intersect in the line l around which the rotation occurs.
//
// In this case the motor has the form: M = l sin(phi) + e1234 cos(phi).
//
// Translations can be covered by a rotation around a point at infinity.
//
// Combined rotations and translations must be created by concatenating the motion
// operations in the sequence gpr(motor_applied_last, motor_applied_first)
////////////////////////////////////////////////////////////////////////////////

// create a motor from a fixed line of rotation and a turning angle
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> get_motor(BiVec3dp<T> const& L, U theta)
{
    // line L must be unitized to avoid surprises
    auto l = unitize(L);

    using ctype = std::common_type_t<T, U>;
    ctype half_angle = 0.5 * theta;
    return MVec3dp_E<ctype>(BiVec3dp<ctype>(l * std::sin(half_angle)),
                            PScalar3dp<ctype>(std::cos(half_angle)));
}

// create a translation motor from a translation vector (given as a Vec3dp).
// Move in direction and by length of translation vector (length = its bulk_nrm)
// ATTENTION: The translation is assumed to be a vector, i.e. with w-component == 0.
//            The w-component is ignored and only the x-, y- and z-components are used.
//            Due to the application via the regressive sandwich product the vector needs
//            to be multiplied by 0.5 before application.
template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> get_motor(Vec3dp<T> const& translation)
{
    return MVec3dp_E<T>(0.5 * BiVec3dp<T>(T(0.0), T(0.0), T(0.0), translation.x,
                                          translation.y, translation.z),
                        PScalar3dp<T>(1.0));
}

// overload of translation motor:
// create a translation motor from a translation vector (given as a Vector3d)
// move in direction and by length of translation vector
template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> get_motor(Vector3d<T> const& translation)
{
    return MVec3dp_E<T>(0.5 * BiVec3dp<T>(T(0.0), T(0.0), T(0.0), translation.x,
                                          translation.y, translation.z),
                        PScalar3dp<T>(1.0));
}

// create a motor from a fixed line of rotation l, a turning angle theta, and a
// distance dist to move along the line
// => screw motion with axis l, turnging angle theta and movement along l by dist
// see e.g. E. Lengyel, "PGA Illuminated", p. 143 (equation 3.93)
template <typename T, typename U, typename V>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U, V>> get_motor(BiVec3dp<T> const& L, U theta,
                                                           V dist)
{
    // line L must be unitized to avoid surprises
    auto l = unitize(L);

    using ctype = std::common_type_t<T, U, V>;
    ctype half_angle = 0.5 * theta;
    ctype half_dist = 0.5 * dist;
    return MVec3dp_E<ctype>(
        Scalar3dp<ctype>(-half_dist * std::sin(half_angle)),
        BiVec3dp<ctype>(l * std::sin(half_angle) -
                        right_weight_dual(l) * half_dist * std::cos(half_angle)),
        PScalar3dp<ctype>(std::cos(half_angle)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>>
get_motor_from_planes(TriVec3dp<T> const& t1, TriVec3dp<U> const& t2)
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
                                                   MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves v (a vector representing a projective point) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(gr1(rgpr(rgpr(R, v), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> move3dp(BiVec3dp<T> const& B,
                                                     MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves B (a bivector representing a line) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(gr2(rgpr(rgpr(R, B), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> move3dp(TriVec3dp<T> const& t,
                                                      MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves t (a trivector representing a plane) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(gr3(rgpr(rgpr(R, t), rrev(R))));
}

// rotate a motor (required for robotics applications using kinematic chains)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> move3dp(MVec3dp_E<T> const& M_orig,
                                                      MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves M_orig (a motor) according to the motor R
    // e.g. kinematic chains in robotics application with coupled joints
    // effectively "rotating the rotation direction"
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(rgpr(rgpr(R, M_orig), rrev(R)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> move3dp_opt(Vec3dp<T> const& v,
                                                       MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;
    // coefficients calculated with ga_prdxpr (pga3dp regressive sandwich product)

    auto k11 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k12 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k13 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);

    auto k14 = 2.0 * (-R.c0 * R.c1 + R.c2 * R.c6 - R.c3 * R.c5 + R.c4 * R.c7);
    auto k21 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k22 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k23 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);
    auto k24 = 2.0 * (-R.c0 * R.c2 - R.c1 * R.c6 + R.c3 * R.c4 + R.c5 * R.c7);

    auto k31 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k32 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k33 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;
    auto k34 = 2.0 * (-R.c0 * R.c3 + R.c1 * R.c5 - R.c2 * R.c4 + R.c6 * R.c7);
    auto k44 = R.c1 * R.c1 + R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    return Vec3dp<ctype>(k11 * v.x + k12 * v.y + k13 * v.z + k14 * v.w,
                         k21 * v.x + k22 * v.y + k23 * v.z + k24 * v.w,
                         k31 * v.x + k32 * v.y + k33 * v.z + k34 * v.w, k44 * v.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::vector<Vec3dp<std::common_type_t<T, U>>>
move3dp(std::vector<Vec3dp<T>> const& vec, MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;
    // coefficients calculated with ga_prdxpr (pga3dp regressive sandwich product)

    auto k11 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k12 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k13 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);
    auto k14 = 2.0 * (-R.c0 * R.c1 + R.c2 * R.c6 - R.c3 * R.c5 + R.c4 * R.c7);

    auto k21 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k22 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k23 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);
    auto k24 = 2.0 * (-R.c0 * R.c2 - R.c1 * R.c6 + R.c3 * R.c4 + R.c5 * R.c7);

    auto k31 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k32 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k33 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;
    auto k34 = 2.0 * (-R.c0 * R.c3 + R.c1 * R.c5 - R.c2 * R.c4 + R.c6 * R.c7);

    auto k44 = R.c1 * R.c1 + R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    std::vector<Vec3dp<ctype>> result;
    result.reserve(vec.size());

    for (auto const& v : vec) {
        result.emplace_back(Vec3dp<ctype>(k11 * v.x + k12 * v.y + k13 * v.z + k14 * v.w,
                                          k21 * v.x + k22 * v.y + k23 * v.z + k24 * v.w,
                                          k31 * v.x + k32 * v.y + k33 * v.z + k34 * v.w,
                                          k44 * v.w));
    }
    return result;
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> move3dp_opt(BiVec3dp<T> const& B,
                                                         MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves B a bivector representing a line according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;
    // coefficients calculated with ga_prdxpr (pga3dp regressive sandwich product)

    auto k11 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k12 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k13 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);

    auto k21 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k22 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k23 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);

    auto k31 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k32 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k33 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    auto k41 = 2.0 * (R.c0 * R.c7 + R.c1 * R.c4 - R.c2 * R.c5 - R.c3 * R.c6);
    auto k42 = 2.0 * (-R.c0 * R.c3 + R.c1 * R.c5 + R.c2 * R.c4 - R.c6 * R.c7);
    auto k43 = 2.0 * (R.c0 * R.c2 + R.c1 * R.c6 + R.c3 * R.c4 + R.c5 * R.c7);
    auto k44 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k45 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k46 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);

    auto k51 = 2.0 * (R.c0 * R.c3 + R.c1 * R.c5 + R.c2 * R.c4 + R.c6 * R.c7);
    auto k52 = 2.0 * (R.c0 * R.c7 - R.c1 * R.c4 + R.c2 * R.c5 - R.c3 * R.c6);
    auto k53 = 2.0 * (-R.c0 * R.c1 + R.c2 * R.c6 + R.c3 * R.c5 - R.c4 * R.c7);
    auto k54 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k55 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k56 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);

    auto k61 = 2.0 * (-R.c0 * R.c2 + R.c1 * R.c6 + R.c3 * R.c4 - R.c5 * R.c7);
    auto k62 = 2.0 * (R.c0 * R.c1 + R.c2 * R.c6 + R.c3 * R.c5 + R.c4 * R.c7);
    auto k63 = 2.0 * (R.c0 * R.c7 - R.c1 * R.c4 - R.c2 * R.c5 + R.c3 * R.c6);
    auto k64 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k65 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k66 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    return BiVec3dp<ctype>(
        k11 * B.vx + k12 * B.vy + k13 * B.vz, k21 * B.vx + k22 * B.vy + k23 * B.vz,
        k31 * B.vx + k32 * B.vy + k33 * B.vz,
        k41 * B.vx + k42 * B.vy + k43 * B.vz + k44 * B.mx + k45 * B.my * k46 * B.mz,
        k51 * B.vx + k52 * B.vy + k53 * B.vz + k54 * B.mx + k55 * B.my * k56 * B.mz,
        k61 * B.vx + k62 * B.vy + k63 * B.vz + k64 * B.mx + k65 * B.my * k66 * B.mz);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::vector<BiVec3dp<std::common_type_t<T, U>>>
move3dp(std::vector<BiVec3dp<T>> const& bvec, MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves B a bivector representing a line according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;
    // coefficients calculated with ga_prdxpr (pga3dp regressive sandwich product)

    auto k11 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k12 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k13 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);

    auto k21 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k22 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k23 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);

    auto k31 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k32 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k33 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    auto k41 = 2.0 * (R.c0 * R.c7 + R.c1 * R.c4 - R.c2 * R.c5 - R.c3 * R.c6);
    auto k42 = 2.0 * (-R.c0 * R.c3 + R.c1 * R.c5 + R.c2 * R.c4 - R.c6 * R.c7);
    auto k43 = 2.0 * (R.c0 * R.c2 + R.c1 * R.c6 + R.c3 * R.c4 + R.c5 * R.c7);
    auto k44 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k45 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k46 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);

    auto k51 = 2.0 * (R.c0 * R.c3 + R.c1 * R.c5 + R.c2 * R.c4 + R.c6 * R.c7);
    auto k52 = 2.0 * (R.c0 * R.c7 - R.c1 * R.c4 + R.c2 * R.c5 - R.c3 * R.c6);
    auto k53 = 2.0 * (-R.c0 * R.c1 + R.c2 * R.c6 + R.c3 * R.c5 - R.c4 * R.c7);
    auto k54 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k55 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k56 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);

    auto k61 = 2.0 * (-R.c0 * R.c2 + R.c1 * R.c6 + R.c3 * R.c4 - R.c5 * R.c7);
    auto k62 = 2.0 * (R.c0 * R.c1 + R.c2 * R.c6 + R.c3 * R.c5 + R.c4 * R.c7);
    auto k63 = 2.0 * (R.c0 * R.c7 - R.c1 * R.c4 - R.c2 * R.c5 + R.c3 * R.c6);
    auto k64 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k65 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k66 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    std::vector<BiVec3dp<ctype>> result;
    result.reserve(bvec.size());

    for (auto const& B : bvec) {
        result.emplace_back(BiVec3dp<ctype>(
            k11 * B.vx + k12 * B.vy + k13 * B.vz, k21 * B.vx + k22 * B.vy + k23 * B.vz,
            k31 * B.vx + k32 * B.vy + k33 * B.vz,
            k41 * B.vx + k42 * B.vy + k43 * B.vz + k44 * B.mx + k45 * B.my * k46 * B.mz,
            k51 * B.vx + k52 * B.vy + k53 * B.vz + k54 * B.mx + k55 * B.my * k56 * B.mz,
            k61 * B.vx + k62 * B.vy + k63 * B.vz + k64 * B.mx + k65 * B.my * k66 * B.mz));
    }
    return result;
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> move3dp_opt(TriVec3dp<T> const& t,
                                                          MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves t (a trivector representing a plane) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;
    // coefficients calculated with ga_prdxpr (pga3dp regressive sandwich product)

    auto k11 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k12 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k13 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);

    auto k21 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k22 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k23 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);

    auto k31 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k32 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k33 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    auto k41 = 2.0 * (R.c0 * R.c1 + R.c2 * R.c6 - R.c3 * R.c5 - R.c4 * R.c7);
    auto k42 = 2.0 * (R.c0 * R.c2 - R.c1 * R.c6 + R.c3 * R.c4 - R.c5 * R.c7);
    auto k43 = 2.0 * (R.c0 * R.c3 + R.c1 * R.c5 - R.c2 * R.c4 - R.c6 * R.c7);
    auto k44 = R.c1 * R.c1 + R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    return TriVec3dp<ctype>(
        k11 * t.x + k12 * t.y + k13 * t.z, k21 * t.x + k22 * t.y + k23 * t.z,
        k31 * t.x + k32 * t.y + k33 * t.z, k41 * t.x + k42 * t.y + k43 * t.z + k44 * t.w);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::vector<TriVec3dp<std::common_type_t<T, U>>>
move3dp(std::vector<TriVec3dp<T>> const& tvec, MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;
    // coefficients calculated with ga_prdxpr (pga3dp regressive sandwich product)

    auto k11 = R.c1 * R.c1 - R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k12 = 2.0 * (R.c1 * R.c2 - R.c3 * R.c7);
    auto k13 = 2.0 * (R.c1 * R.c3 + R.c2 * R.c7);

    auto k21 = 2.0 * (R.c1 * R.c2 + R.c3 * R.c7);
    auto k22 = -R.c1 * R.c1 + R.c2 * R.c2 - R.c3 * R.c3 + R.c7 * R.c7;
    auto k23 = 2.0 * (-R.c1 * R.c7 + R.c2 * R.c3);

    auto k31 = 2.0 * (R.c1 * R.c3 - R.c2 * R.c7);
    auto k32 = 2.0 * (R.c1 * R.c7 + R.c2 * R.c3);
    auto k33 = -R.c1 * R.c1 - R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    auto k41 = 2.0 * (R.c0 * R.c1 + R.c2 * R.c6 - R.c3 * R.c5 - R.c4 * R.c7);
    auto k42 = 2.0 * (R.c0 * R.c2 - R.c1 * R.c6 + R.c3 * R.c4 - R.c5 * R.c7);
    auto k43 = 2.0 * (R.c0 * R.c3 + R.c1 * R.c5 - R.c2 * R.c4 - R.c6 * R.c7);
    auto k44 = R.c1 * R.c1 + R.c2 * R.c2 + R.c3 * R.c3 + R.c7 * R.c7;

    std::vector<TriVec3dp<ctype>> result;
    result.reserve(tvec.size());

    for (auto const& t : tvec) {
        result.emplace_back(TriVec3dp<ctype>(
            k11 * t.x + k12 * t.y + k13 * t.z, k21 * t.x + k22 * t.y + k23 * t.z,
            k31 * t.x + k32 * t.y + k33 * t.z,
            k41 * t.x + k42 * t.y + k43 * t.z + k44 * t.w));
    }
    return result;
}


template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp<std::common_type_t<T, U>> move3dp(MVec3dp<T> const& MV,
                                                    MVec3dp_E<U> const& M)
{
    // motor M must be unitized to avoid surprises
    auto R = unitize(M);

    // rotate one multivector M with motor R
    using ctype = std::common_type_t<T, U>;
    return MVec3dp<ctype>(rgpr(rgpr(R, MV), rrev(R)));
}


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
// Projections for 3dp: (HINT: unitize after projection, if not at infinity)
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
    auto p = rwdg(std::forward<arg2>(b),
                  right_weight_expand3dp(std::forward<arg1>(a), std::forward<arg2>(b)));

    // return a unitized object, if it is not located in the horizon
    if (weight_nrm_sq(p) != 0.0) {
        p = unitize(p);
    }

    return p;
}

template <typename arg1, typename arg2> decltype(auto) central_proj3dp(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!

    // project the smaller grade object onto to larger grade object
    auto p = rwdg(std::forward<arg2>(b),
                  right_bulk_expand3dp(std::forward<arg1>(a), std::forward<arg2>(b)));

    // return a unitized object, if it is not located in the horizon
    if (weight_nrm_sq(p) != 0.0) {
        p = unitize(p);
    }

    return p;
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
// inversion
// (=point reflection resuling from reflecting across three perpendicular planes)
// the resulting minus comes for the 3 consecutive reflections
////////////////////////////////////////////////////////////////////////////////

// (point-)reflect a point q in an arbitrary point p
// p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec3dp<std::common_type_t<T, U>> invert_on(Vec3dp<T> const& q,
                                                     Vec3dp<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3dp<ctype>(-gr1(rgpr(rgpr(p, q), rrev(p))));
}

// (point-)reflect a line l in an arbitrary point p
// p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr BiVec3dp<std::common_type_t<T, U>> invert_on(BiVec3dp<T> const& l,
                                                       Vec3dp<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3dp<ctype>(-gr2(rgpr(rgpr(p, l), rrev(p))));
}

// (point-)reflect a plane t in an arbitrary point p
// p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr TriVec3dp<std::common_type_t<T, U>> invert_on(TriVec3dp<T> const& t,
                                                        Vec3dp<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec3dp<ctype>(-gr3(rgpr(rgpr(p, t), rrev(p))));
}


////////////////////////////////////////////////////////////////////////////////
// support: nearest point to the origin (by projecting origin onto line/plane)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
Vec3dp<T> support(BiVec3dp<T> const& B)
{
    // REQUIRES: a line (BiVec3dp) as argument

    // project origin onto line
    return ortho_proj3dp(O_3dp, B);
}

template <typename T>
    requires(std::floating_point<T>)
Vec3dp<T> support(TriVec3dp<T> const& t)
{
    // REQUIRES: a plane (TriVec3dp) as argument

    // project origin onto line
    return ortho_proj3dp(O_3dp, t);
}


////////////////////////////////////////////////////////////////////////////////
// attitude operations: att = rwdg( u, rcmpl(e4_3dp) ) = rwdg(u, H_3dp)
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
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), T(ps));
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

////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent(Scalar3dp<T> a, Scalar3dp<U> b, value_t tolerance = eps)
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
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent(Vec3dp<T> const& a, Vec3dp<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check (bulk_nrm_sq can be zero for non-zero
    // PGA elements)
    bool a_is_zero = (std::abs(a.x) < tolerance) && (std::abs(a.y) < tolerance) &&
                     (std::abs(a.z) < tolerance) && (std::abs(a.w) < tolerance);
    bool b_is_zero = (std::abs(b.x) < tolerance) && (std::abs(b.y) < tolerance) &&
                     (std::abs(b.z) < tolerance) && (std::abs(b.w) < tolerance);

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
    else if (std::abs(b.w) > tolerance) {
        k = a.w / b.w;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol = tolerance * std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z),
                                            std::abs(a.w), std::abs(b.x), std::abs(b.y),
                                            std::abs(b.z), std::abs(b.w), value_t(1.0)});
    return (std::abs(a.x - k * b.x) < rel_tol) && (std::abs(a.y - k * b.y) < rel_tol) &&
           (std::abs(a.z - k * b.z) < rel_tol) && (std::abs(a.w - k * b.w) < rel_tol);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent(BiVec3dp<T> const& a, BiVec3dp<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.vx) < tolerance) && (std::abs(a.vy) < tolerance) &&
                     (std::abs(a.vz) < tolerance) && (std::abs(a.mx) < tolerance) &&
                     (std::abs(a.my) < tolerance) && (std::abs(a.mz) < tolerance);
    bool b_is_zero = (std::abs(b.vx) < tolerance) && (std::abs(b.vy) < tolerance) &&
                     (std::abs(b.vz) < tolerance) && (std::abs(b.mx) < tolerance) &&
                     (std::abs(b.my) < tolerance) && (std::abs(b.mz) < tolerance);

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
    if (std::abs(b.vx) > tolerance) {
        k = a.vx / b.vx;
        k_found = true;
    }
    else if (std::abs(b.vy) > tolerance) {
        k = a.vy / b.vy;
        k_found = true;
    }
    else if (std::abs(b.vz) > tolerance) {
        k = a.vz / b.vz;
        k_found = true;
    }
    else if (std::abs(b.mx) > tolerance) {
        k = a.mx / b.mx;
        k_found = true;
    }
    else if (std::abs(b.my) > tolerance) {
        k = a.my / b.my;
        k_found = true;
    }
    else if (std::abs(b.mz) > tolerance) {
        k = a.mz / b.mz;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol =
        tolerance *
        std::max({std::abs(a.vx), std::abs(a.vy), std::abs(a.vz), std::abs(a.mx),
                  std::abs(a.my), std::abs(a.mz), std::abs(b.vx), std::abs(b.vy),
                  std::abs(b.vz), std::abs(b.mx), std::abs(b.my), std::abs(b.mz),
                  value_t(1.0)});
    return (std::abs(a.vx - k * b.vx) < rel_tol) &&
           (std::abs(a.vy - k * b.vy) < rel_tol) &&
           (std::abs(a.vz - k * b.vz) < rel_tol) &&
           (std::abs(a.mx - k * b.mx) < rel_tol) &&
           (std::abs(a.my - k * b.my) < rel_tol) && (std::abs(a.mz - k * b.mz) < rel_tol);
}

// For trivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent(TriVec3dp<T> const& a, TriVec3dp<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.x) < tolerance) && (std::abs(a.y) < tolerance) &&
                     (std::abs(a.z) < tolerance) && (std::abs(a.w) < tolerance);
    bool b_is_zero = (std::abs(b.x) < tolerance) && (std::abs(b.y) < tolerance) &&
                     (std::abs(b.z) < tolerance) && (std::abs(b.w) < tolerance);

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
    else if (std::abs(b.w) > tolerance) {
        k = a.w / b.w;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol = tolerance * std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z),
                                            std::abs(a.w), std::abs(b.x), std::abs(b.y),
                                            std::abs(b.z), std::abs(b.w), value_t(1.0)});
    return (std::abs(a.x - k * b.x) < rel_tol) && (std::abs(a.y - k * b.y) < rel_tol) &&
           (std::abs(a.z - k * b.z) < rel_tol) && (std::abs(a.w - k * b.w) < rel_tol);
}

// For pseudoscalars: all non-zero pseudoscalars in 3DP represent the same subspace
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent(PScalar3dp<T> a, PScalar3dp<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars in 3DP are congruent (represent the full 4D space)
    return true;
}

} // namespace hd::ga::pga