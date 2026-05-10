#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_sta4ds_ops_basics.hpp"
#include "ga_sta4ds_ops_products.hpp"


/////////////////////////////////////////////////////////////////////////////////////////
// TODO: all functions here copy paste from pga3dp -> nees adaptation and testing!
/////////////////////////////////////////////////////////////////////////////////////////


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// provides functionality that is based on pga4ds ops basics and products:
//
// - angle()                              -> angle operations
// - exp()                                -> exponential (w.r.t. rgpr)
// TODO: - log()                                -> logarithm (w.r.t rgpr)
// - sqrt(M)                              -> sqrt of a motor (w.r.t. rgpr)
// - get_motor()                          -> provide a motor from (line, phi), or (delta),
//                                           or (line, phi, dist along line)
// - get_motor_from_planes()              -> provide a motor (from two plane reflections)
// - get_motor_from_lines()               -> provide a motor (from two lines moved
//                                                            into each other)
// - move4ds(), move4ds_opt()             -> move object with motor
// - project_onto(), reject_from()        -> simple projection and rejection
// - expand()                             -> expansion: new line/plane through point/line
//                                                      perpendicular to line/plane
// - ortho_proj4ds()                      -> orthogonal projection onto object
// - central_proj4ds()                    -> central projection towards origin onto object
// - ortho_antiproj4ds()                  -> orthogonal antiprojection onto object
// - reflect_on()                         -> reflections
// - invert_on()                          -> inversions
// - sup()                                -> point on line/plane that is nearest to origin
// - att()                                -> object attitude
// - dist4ds()                            -> Euclidean distance and homogeneous magnitude
// - is_congruent()                       -> Same up to a scalar factor (is same subspace)
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// angle operations 4ds
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors, i.e. directions to points at infinity
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> angle(Vec4ds<T> const& v1, Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    if ((weight_nrm_sq(v1) != 0.0) || (weight_nrm_sq(v2) != 0.0)) {
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
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::common_type_t<T, U> angle(BiVec4ds<T> const& B1, BiVec4ds<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(r_weight_contract4ds(B1, B2));
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
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::common_type_t<T, U> angle(TriVec4ds<T> const& t, BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(bulk_nrm(r_weight_contract4ds(t, B)));
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
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::common_type_t<T, U> angle(BiVec4ds<T> const& B, TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(bulk_nrm(r_weight_contract4ds(t, B)));
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
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::common_type_t<T, U> angle(TriVec4ds<T> const& t1, TriVec4ds<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(r_weight_contract4ds(t1, t2));
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
// exp() with respect to the regressive geometric product rgpr()
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> exp(BiVec4ds<T> const& B)
{
    T phi_sq = weight_nrm_sq(B); // rotation angle^2
    if (phi_sq == 0.0) {
        // pure translation

        // B = att(r_bulk_dual(delta)) to move in direction of vector delta
        // B does only contain a bulk part which encodes the translation
        return MVec4ds_E<T>(Scalar4ds<T>(0.0), B, PScalar4ds<T>(1.0));
    }

    // rotation with phi != 0
    T phi = std::sqrt(phi_sq); // rotation angle phi around unitized l
    T phi_inv = 1.0 / phi;
    T dot_lvlm = B.vx * B.mx + B.vy * B.my + B.vz * B.mz; // are Bv and Bm orthogonal?

    // retrieve line l (as in case of pure rotation)
    auto l = phi_inv * B;
    double dist{0.0}; // default: no translation, i.e. dist = 0

    if (std::abs(dot_lvlm) > eps) {
        // case where B contains a screw motion with dist * Bv encoded in Bm.
        // Thus we need to extract the unitized line l and the translation distance.

        // We get the unitized line by projecting the momentum vector lm onto
        // the direction vector lv and by removing the non-orthogonal part
        dist = dot_lvlm * phi_inv; // length of projection of lm in direction of lv
        l = l - dist * phi_inv * phi_inv * bivec4ds{0.0, 0.0, 0.0, B.vx, B.vy, B.vz};
    }

    return MVec4ds_E<T>(
        Scalar4ds<T>(-dist * std::sin(phi)),
        BiVec4ds<T>(l * std::sin(phi) - r_weight_dual(l) * dist * std::cos(phi)),
        PScalar4ds<T>(std::cos(phi)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> sqrt(MVec4ds_E<T> const& M)
{
    if (std::abs(M.c0) < eps) {
        // simple motor, if M.c0 == 0.0 (i.e. no scalar part)
        return unitize(M + PScalar4ds<T>(1.0));
    }

    // non-simple motor for other cases (s. Lengyel, "PGA illuminated", p. 151)
    return (M + PScalar4ds<T>(1.0)) / std::sqrt(2.0 + 2.0 * M.c7) -
           rgpr((M + PScalar4ds<T>(1.0)) / std::sqrt(2.0 + 2.0 * M.c7),
                Scalar4ds<T>(M.c0 / (2.0 + 2.0 * M.c7)));
}

////////////////////////////////////////////////////////////////////////////////
// 4ds motor operations (translation and rotation)
//
// Every motor in pga4ds is an even-grade multivector MVec4ds_E.
//
// A proper isometry in 4ds has a fixed line l around which a rotation occurs
// with an angle phi. The line is represented by a bivector.
//
// The rotation so modelled by two consecutive reflections across two planes which
// intersect in the line l around which the rotation occurs.
//
// In this case the motor has the form: M = l sin(phi) + e1234 cos(phi).
//
// Translations can be covered by a rotation around a line at infinity.
//
// Combined rotations and translations must be created by concatenating the motion
// operations in the sequence rgpr(motor_applied_last, motor_applied_first)
////////////////////////////////////////////////////////////////////////////////

// create a motor from a fixed line of rotation and a turning angle
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> get_motor(BiVec4ds<T> const& L, U theta)
{
    // line L must be unitized to avoid surprises
    auto L_weight_nrm_sq = weight_nrm_sq(L);
    if (L_weight_nrm_sq < eps) {
        throw std::invalid_argument(
            "get_motor: Cannot use ideal lines L with weight_nrm_sq(L) == 0.0");
    }
    auto l{L};
    if ((L_weight_nrm_sq > eps) && (L_weight_nrm_sq != 1.0)) {
        l = unitize(L);
    }

    using ctype = std::common_type_t<T, U>;
    ctype half_angle = 0.5 * theta;
    return MVec4ds_E<ctype>(BiVec4ds<ctype>(l * std::sin(half_angle)),
                            PScalar4ds<ctype>(std::cos(half_angle)));
}

// Create a translation motor from a translation
// delta = vec4ds(delta_x, delta_y, delta_z, 0)
// Move in direction and by length of translation vector (length = bulk_nrm(delta))
// ATTENTION: The translation is assumed to be a direction, i.e. with w-component == 0.
//            The w-component is ignored and only the x-, y- and z-components are used.
//            Due to the application via the regressive sandwich product the vector needs
//            to be multiplied by 0.5
//
// The function implements a motor calculated from exp_{\veedot}(B_{tra})
// with the bivector B_{tra} = 0.5 * bulk_dual(delta) \vee H_{4ds}
//                           = 0.5 * att(bulk_dual(delta))
template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> get_motor(Vec4ds<T> const& delta)
{
    return MVec4ds_E<T>(
        0.5 * BiVec4ds<T>(T(0.0), T(0.0), T(0.0), delta.x, delta.y, delta.z),
        PScalar4ds<T>(1.0));
}


// create a motor from a fixed line of rotation l, a turning angle theta, and a
// distance dist to move along the line
// => screw motion with axis l, turnging angle theta and movement along l by dist
// see e.g. E. Lengyel, "PGA Illuminated", p. 143 (equation 3.93)
template <typename T, typename U, typename V>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U, V>> get_motor(BiVec4ds<T> const& L, U theta,
                                                           V dist)
{
    // line L must be unitized to avoid surprises
    auto L_weight_nrm_sq = weight_nrm_sq(L);
    auto l{L};
    if ((L_weight_nrm_sq > eps) && (L_weight_nrm_sq != 1.0)) {
        l = unitize(L);
    }

    using ctype = std::common_type_t<T, U, V>;
    ctype half_angle = 0.5 * theta;
    ctype half_dist = 0.5 * dist;
    return MVec4ds_E<ctype>(
        Scalar4ds<ctype>(-half_dist * std::sin(half_angle)),
        BiVec4ds<ctype>(l * std::sin(half_angle) -
                        r_weight_dual(l) * half_dist * std::cos(half_angle)),
        PScalar4ds<ctype>(std::cos(half_angle)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>>
get_motor_from_planes(TriVec4ds<T> const& t1, TriVec4ds<U> const& t2)
{
    // take planes (=trivectors) as input and return a motor M
    // 1st apply reflection across plane t1, then across t2 to get a motor that rotates
    // around the intersection line of planes t1 and t2
    // or translates when t1 and t2 are parallel
    //
    // for use of motor M either directly on object u:
    //
    //     auto v_moved = gr1( rgpr(rgpr(M, v), rrev(M)) );
    // or
    //     auto B_moved = gr2( rgpr(rgpr(M, B), rrev(M)) );
    // or
    //     auto t_moved = gr3( rgpr(rgpr(M, t), rrev(M)) );
    // or
    //     auto v_moved = move4ds(v,M);  // moves v according to the motor M
    //     auto B_moved = move4ds(B,M);  // moves B according to the motor M
    //     auto t_moved = move4ds(t,M);  // moves t according to the motor M

    // planes t1 and t2 need to be unitized to avoid surprises
    auto w_nrm_sq = weight_nrm_sq(t1);
    auto tu1{t1};
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        tu1 = unitize(tu1);
    }
    w_nrm_sq = weight_nrm_sq(t2);
    auto tu2{t2};
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        tu2 = unitize(tu2);
    }

    // the resulting motor must be unitized to avoid surprises
    auto M{rgpr(tu2, tu1)};
    w_nrm_sq = weight_nrm_sq(M);
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        M = unitize(M);
    }
    return M; // based on the regressive geometric product
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> get_motor_from_lines(BiVec4ds<T> const& l1,
                                                                   BiVec4ds<U> const& l2)
{
    // Return a motor that moves line l1 to line l2.
    //
    // Line of rotation is l_rot = rcmt(l2, rrev(l1)).
    // The operator rotating by 2*phi and moving by 2*dist is rgpr(l2, rrev(l1))
    // Thus, sqrt(rgpr(l2, rrev(l1))) is the motor we look for rotating by phi and moving
    // by dist.
    // (see Lengyel, "PGA Illuminated", p. 152)

    // l1 and l2 need to be unitized to avoid surprises
    auto w_nrm_sq = weight_nrm_sq(l1);
    auto lu1{l1};
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        lu1 = unitize(lu1);
    }
    w_nrm_sq = weight_nrm_sq(l2);
    auto lu2{l2};
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        lu2 = unitize(lu2);
    }

    // the resulting motor must be unitized to avoid surprises
    auto M{sqrt(rgpr(l2, rrev(l1)))};
    w_nrm_sq = weight_nrm_sq(M);
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        M = unitize(M);
    }
    return M; // based on the regressive geometric product
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> move4ds(Vec4ds<T> const& v,
                                                   MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves v (a vector representing a projective point) according to the motor M
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(gr1(rgpr(rgpr(M, v), rrev(M))));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> move4ds(BiVec4ds<T> const& B,
                                                     MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves B (a bivector representing a line) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(gr2(rgpr(rgpr(M, B), rrev(M))));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> move4ds(TriVec4ds<T> const& t,
                                                      MVec4ds_E<U> const& M)
{
    // motor M must be unitized to avoid surprises

    // moves t (a trivector representing a plane) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return TriVec4ds<ctype>(gr3(rgpr(rgpr(M, t), rrev(M))));
}

// rotate a motor (required for robotics applications using kinematic chains)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> move4ds(MVec4ds_E<T> const& M_orig,
                                                      MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves M_orig (a motor) according to the motor M
    // e.g. kinematic chains in robotics application with coupled joints
    // effectively "rotating the rotation direction"
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(rgpr(rgpr(M, M_orig), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> move4ds(MVec4ds<T> const& MV,
                                                    MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // rotate one multivector M with motor M
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(rgpr(rgpr(M, MV), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> move4ds_opt(Vec4ds<T> const& v,
                                                       MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga4ds regressive sandwich product)
    ctype const h0 = M.c1 * M.c1;
    ctype const h1 = M.c2 * M.c2;
    ctype const h2 = M.c3 * M.c3;
    ctype const h3 = M.c7 * M.c7;
    ctype const h4 = M.c1 * M.c2;
    ctype const h5 = M.c3 * M.c7;
    ctype const h6 = M.c1 * M.c3;
    ctype const h7 = M.c2 * M.c7;
    ctype const h8 = M.c1 * M.c7;
    ctype const h9 = M.c2 * M.c3;

    ctype const k11 = h0 - h1 - h2 + h3;
    ctype const k12 = 2.0 * (h4 - h5);
    ctype const k13 = 2.0 * (h6 + h7);
    ctype const k14 = 2.0 * (-M.c0 * M.c1 + M.c2 * M.c6 - M.c3 * M.c5 + M.c4 * M.c7);
    ctype const k21 = 2.0 * (h4 + h5);
    ctype const k22 = -h0 + h1 - h2 + h3;
    ctype const k23 = 2.0 * (-h8 + h9);
    ctype const k24 = 2.0 * (-M.c0 * M.c2 - M.c1 * M.c6 + M.c3 * M.c4 + M.c5 * M.c7);
    ctype const k31 = 2.0 * (h6 - h7);
    ctype const k32 = 2.0 * (h8 + h9);
    ctype const k33 = -h0 - h1 + h2 + h3;
    ctype const k34 = 2.0 * (-M.c0 * M.c3 + M.c1 * M.c5 - M.c2 * M.c4 + M.c6 * M.c7);
    ctype const k44 = h0 + h1 + h2 + h3;

    return Vec4ds<ctype>(k11 * v.x + k12 * v.y + k13 * v.z + k14 * v.w,
                         k21 * v.x + k22 * v.y + k23 * v.z + k24 * v.w,
                         k31 * v.x + k32 * v.y + k33 * v.z + k34 * v.w, k44 * v.w);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::vector<Vec4ds<std::common_type_t<T, U>>>
move4ds(std::vector<Vec4ds<T>> const& vec, MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga4ds regressive sandwich product)
    ctype const h0 = M.c1 * M.c1;
    ctype const h1 = M.c2 * M.c2;
    ctype const h2 = M.c3 * M.c3;
    ctype const h3 = M.c7 * M.c7;
    ctype const h4 = M.c1 * M.c2;
    ctype const h5 = M.c3 * M.c7;
    ctype const h6 = M.c1 * M.c3;
    ctype const h7 = M.c2 * M.c7;
    ctype const h8 = M.c1 * M.c7;
    ctype const h9 = M.c2 * M.c3;

    ctype const k11 = h0 - h1 - h2 + h3;
    ctype const k12 = 2.0 * (h4 - h5);
    ctype const k13 = 2.0 * (h6 + h7);
    ctype const k14 = 2.0 * (-M.c0 * M.c1 + M.c2 * M.c6 - M.c3 * M.c5 + M.c4 * M.c7);
    ctype const k21 = 2.0 * (h4 + h5);
    ctype const k22 = -h0 + h1 - h2 + h3;
    ctype const k23 = 2.0 * (-h8 + h9);
    ctype const k24 = 2.0 * (-M.c0 * M.c2 - M.c1 * M.c6 + M.c3 * M.c4 + M.c5 * M.c7);
    ctype const k31 = 2.0 * (h6 - h7);
    ctype const k32 = 2.0 * (h8 + h9);
    ctype const k33 = -h0 - h1 + h2 + h3;
    ctype const k34 = 2.0 * (-M.c0 * M.c3 + M.c1 * M.c5 - M.c2 * M.c4 + M.c6 * M.c7);
    ctype const k44 = h0 + h1 + h2 + h3;

    std::vector<Vec4ds<ctype>> result;
    result.reserve(vec.size());

    for (auto const& v : vec) {
        result.emplace_back(Vec4ds<ctype>(k11 * v.x + k12 * v.y + k13 * v.z + k14 * v.w,
                                          k21 * v.x + k22 * v.y + k23 * v.z + k24 * v.w,
                                          k31 * v.x + k32 * v.y + k33 * v.z + k34 * v.w,
                                          k44 * v.w));
    }
    return result;
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> move4ds_opt(BiVec4ds<T> const& B,
                                                         MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves B a bivector representing a line according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga4ds regressive sandwich product)
    ctype const h0 = M.c1 * M.c1;
    ctype const h1 = M.c2 * M.c2;
    ctype const h2 = M.c3 * M.c3;
    ctype const h3 = M.c7 * M.c7;
    ctype const h4 = M.c1 * M.c2;
    ctype const h5 = M.c3 * M.c7;
    ctype const h6 = M.c1 * M.c3;
    ctype const h7 = M.c2 * M.c7;
    ctype const h8 = M.c1 * M.c7;
    ctype const h9 = M.c2 * M.c3;
    ctype const h10 = M.c0 * M.c7;
    ctype const h11 = M.c1 * M.c4;
    ctype const h12 = M.c2 * M.c5;
    ctype const h13 = M.c3 * M.c6;
    ctype const h14 = M.c0 * M.c3;
    ctype const h15 = M.c1 * M.c5;
    ctype const h16 = M.c2 * M.c4;
    ctype const h17 = M.c6 * M.c7;
    ctype const h18 = M.c0 * M.c2;
    ctype const h19 = M.c1 * M.c6;
    ctype const h20 = M.c3 * M.c4;
    ctype const h21 = M.c5 * M.c7;
    ctype const h22 = M.c0 * M.c1;
    ctype const h23 = M.c2 * M.c6;
    ctype const h24 = M.c3 * M.c5;
    ctype const h25 = M.c4 * M.c7;

    ctype const k11 = h0 - h1 - h2 + h3;
    ctype const k12 = 2.0 * (h4 - h5);
    ctype const k13 = 2.0 * (h6 + h7);
    ctype const k21 = 2.0 * (h4 + h5);
    ctype const k22 = -h0 + h1 - h2 + h3;
    ctype const k23 = 2.0 * (-h8 + h9);
    ctype const k31 = 2.0 * (h6 - h7);
    ctype const k32 = 2.0 * (h8 + h9);
    ctype const k33 = -h0 - h1 + h2 + h3;

    ctype const k41 = 2.0 * (h10 + h11 - h12 - h13);
    ctype const k42 = 2.0 * (-h14 + h15 + h16 - h17);
    ctype const k43 = 2.0 * (h18 + h19 + h20 + h21);
    ctype const k51 = 2.0 * (h14 + h15 + h16 + h17);
    ctype const k52 = 2.0 * (h10 - h11 + h12 - h13);
    ctype const k53 = 2.0 * (-h22 + h23 + h24 - h25);
    ctype const k61 = 2.0 * (-h18 + h19 + h20 - h21);
    ctype const k62 = 2.0 * (h22 + h23 + h24 + h25);
    ctype const k63 = 2.0 * (h10 - h11 - h12 + h13);

    return BiVec4ds<ctype>(
        k11 * B.vx + k12 * B.vy + k13 * B.vz, k21 * B.vx + k22 * B.vy + k23 * B.vz,
        k31 * B.vx + k32 * B.vy + k33 * B.vz,
        k41 * B.vx + k42 * B.vy + k43 * B.vz + k11 * B.mx + k12 * B.my + k13 * B.mz,
        k51 * B.vx + k52 * B.vy + k53 * B.vz + k21 * B.mx + k22 * B.my + k23 * B.mz,
        k61 * B.vx + k62 * B.vy + k63 * B.vz + k31 * B.mx + k32 * B.my + k33 * B.mz);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::vector<BiVec4ds<std::common_type_t<T, U>>>
move4ds(std::vector<BiVec4ds<T>> const& bvec, MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves B a bivector representing a line according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga4ds regressive sandwich product)
    ctype const h0 = M.c1 * M.c1;
    ctype const h1 = M.c2 * M.c2;
    ctype const h2 = M.c3 * M.c3;
    ctype const h3 = M.c7 * M.c7;
    ctype const h4 = M.c1 * M.c2;
    ctype const h5 = M.c3 * M.c7;
    ctype const h6 = M.c1 * M.c3;
    ctype const h7 = M.c2 * M.c7;
    ctype const h8 = M.c1 * M.c7;
    ctype const h9 = M.c2 * M.c3;
    ctype const h10 = M.c0 * M.c7;
    ctype const h11 = M.c1 * M.c4;
    ctype const h12 = M.c2 * M.c5;
    ctype const h13 = M.c3 * M.c6;
    ctype const h14 = M.c0 * M.c3;
    ctype const h15 = M.c1 * M.c5;
    ctype const h16 = M.c2 * M.c4;
    ctype const h17 = M.c6 * M.c7;
    ctype const h18 = M.c0 * M.c2;
    ctype const h19 = M.c1 * M.c6;
    ctype const h20 = M.c3 * M.c4;
    ctype const h21 = M.c5 * M.c7;
    ctype const h22 = M.c0 * M.c1;
    ctype const h23 = M.c2 * M.c6;
    ctype const h24 = M.c3 * M.c5;
    ctype const h25 = M.c4 * M.c7;

    ctype const k11 = h0 - h1 - h2 + h3;
    ctype const k12 = 2.0 * (h4 - h5);
    ctype const k13 = 2.0 * (h6 + h7);
    ctype const k21 = 2.0 * (h4 + h5);
    ctype const k22 = -h0 + h1 - h2 + h3;
    ctype const k23 = 2.0 * (-h8 + h9);
    ctype const k31 = 2.0 * (h6 - h7);
    ctype const k32 = 2.0 * (h8 + h9);
    ctype const k33 = -h0 - h1 + h2 + h3;

    ctype const k41 = 2.0 * (h10 + h11 - h12 - h13);
    ctype const k42 = 2.0 * (-h14 + h15 + h16 - h17);
    ctype const k43 = 2.0 * (h18 + h19 + h20 + h21);
    ctype const k51 = 2.0 * (h14 + h15 + h16 + h17);
    ctype const k52 = 2.0 * (h10 - h11 + h12 - h13);
    ctype const k53 = 2.0 * (-h22 + h23 + h24 - h25);
    ctype const k61 = 2.0 * (-h18 + h19 + h20 - h21);
    ctype const k62 = 2.0 * (h22 + h23 + h24 + h25);
    ctype const k63 = 2.0 * (h10 - h11 - h12 + h13);

    std::vector<BiVec4ds<ctype>> result;
    result.reserve(bvec.size());

    for (auto const& B : bvec) {
        result.emplace_back(BiVec4ds<ctype>(
            k11 * B.vx + k12 * B.vy + k13 * B.vz, k21 * B.vx + k22 * B.vy + k23 * B.vz,
            k31 * B.vx + k32 * B.vy + k33 * B.vz,
            k41 * B.vx + k42 * B.vy + k43 * B.vz + k11 * B.mx + k12 * B.my + k13 * B.mz,
            k51 * B.vx + k52 * B.vy + k53 * B.vz + k21 * B.mx + k22 * B.my + k23 * B.mz,
            k61 * B.vx + k62 * B.vy + k63 * B.vz + k31 * B.mx + k32 * B.my + k33 * B.mz));
    }
    return result;
}


template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> move4ds_opt(TriVec4ds<T> const& t,
                                                          MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves t (a trivector representing a plane) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga4ds regressive sandwich product)
    ctype const h0 = M.c1 * M.c1;
    ctype const h1 = M.c2 * M.c2;
    ctype const h2 = M.c3 * M.c3;
    ctype const h3 = M.c7 * M.c7;
    ctype const h4 = M.c1 * M.c2;
    ctype const h5 = M.c3 * M.c7;
    ctype const h6 = M.c1 * M.c3;
    ctype const h7 = M.c2 * M.c7;
    ctype const h8 = M.c1 * M.c7;
    ctype const h9 = M.c2 * M.c3;

    ctype const k11 = h0 - h1 - h2 + h3;
    ctype const k12 = 2.0 * (h4 - h5);
    ctype const k13 = 2.0 * (h6 + h7);
    ctype const k21 = 2.0 * (h4 + h5);
    ctype const k22 = -h0 + h1 - h2 + h3;
    ctype const k23 = 2.0 * (-h8 + h9);
    ctype const k31 = 2.0 * (h6 - h7);
    ctype const k32 = 2.0 * (h8 + h9);
    ctype const k33 = -h0 - h1 + h2 + h3;
    ctype const k41 = 2.0 * (M.c0 * M.c1 + M.c2 * M.c6 - M.c3 * M.c5 - M.c4 * M.c7);
    ctype const k42 = 2.0 * (M.c0 * M.c2 - M.c1 * M.c6 + M.c3 * M.c4 - M.c5 * M.c7);
    ctype const k43 = 2.0 * (M.c0 * M.c3 + M.c1 * M.c5 - M.c2 * M.c4 - M.c6 * M.c7);
    ctype const k44 = h0 + h1 + h2 + h3;

    return TriVec4ds<ctype>(
        k11 * t.x + k12 * t.y + k13 * t.z, k21 * t.x + k22 * t.y + k23 * t.z,
        k31 * t.x + k32 * t.y + k33 * t.z, k41 * t.x + k42 * t.y + k43 * t.z + k44 * t.w);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::vector<TriVec4ds<std::common_type_t<T, U>>>
move4ds(std::vector<TriVec4ds<T>> const& tvec, MVec4ds_E<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga4ds regressive sandwich product)
    ctype const h0 = M.c1 * M.c1;
    ctype const h1 = M.c2 * M.c2;
    ctype const h2 = M.c3 * M.c3;
    ctype const h3 = M.c7 * M.c7;
    ctype const h4 = M.c1 * M.c2;
    ctype const h5 = M.c3 * M.c7;
    ctype const h6 = M.c1 * M.c3;
    ctype const h7 = M.c2 * M.c7;
    ctype const h8 = M.c1 * M.c7;
    ctype const h9 = M.c2 * M.c3;

    ctype const k11 = h0 - h1 - h2 + h3;
    ctype const k12 = 2.0 * (h4 - h5);
    ctype const k13 = 2.0 * (h6 + h7);
    ctype const k21 = 2.0 * (h4 + h5);
    ctype const k22 = -h0 + h1 - h2 + h3;
    ctype const k23 = 2.0 * (-h8 + h9);
    ctype const k31 = 2.0 * (h6 - h7);
    ctype const k32 = 2.0 * (h8 + h9);
    ctype const k33 = -h0 - h1 + h2 + h3;
    ctype const k41 = 2.0 * (M.c0 * M.c1 + M.c2 * M.c6 - M.c3 * M.c5 - M.c4 * M.c7);
    ctype const k42 = 2.0 * (M.c0 * M.c2 - M.c1 * M.c6 + M.c3 * M.c4 - M.c5 * M.c7);
    ctype const k43 = 2.0 * (M.c0 * M.c3 + M.c1 * M.c5 - M.c2 * M.c4 - M.c6 * M.c7);
    ctype const k44 = h0 + h1 + h2 + h3;

    std::vector<TriVec4ds<ctype>> result;
    result.reserve(tvec.size());

    for (auto const& t : tvec) {
        result.emplace_back(TriVec4ds<ctype>(
            k11 * t.x + k12 * t.y + k13 * t.z, k21 * t.x + k22 * t.y + k23 * t.z,
            k31 * t.x + k32 * t.y + k33 * t.z,
            k41 * t.x + k42 * t.y + k43 * t.z + k44 * t.w));
    }
    return result;
}


////////////////////////////////////////////////////////////////////////////////
// projections, rejections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
// returns component of v1 parallel to v2
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> project_onto(Vec4ds<T> const& v1,
                                                        Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2); // works directly in representational space
}

// rejection of vector v1 from a vector v2
// returns component of v1 perpendicular to v2
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reject_from(Vec4ds<T> const& v1,
                                                       Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec4ds<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}


// orthogonal projection of a vector v onto a bivector B (a line)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> project_onto(Vec4ds<T> const& v,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(rwdg(B, wdg(v, r_weight_dual(B)))); // ortho_proj4ds
}

// rejection of vector v from a bivector B (a line)
// rejection = v - project_onto(v, B)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reject_from(Vec4ds<T> const& v,
                                                       BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(v - project_onto(v, B));
}


// orthogonal projection of a vector v onto a trivector t (a plane)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> project_onto(Vec4ds<T> const& v,
                                                        TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(rwdg(t, wdg(v, r_weight_dual(t)))); // ortho_proj4ds
}

// rejection of vector v from a trivector t (a plane)
// rejection = v - project_onto(v, t)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reject_from(Vec4ds<T> const& v,
                                                       TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(v - project_onto(v, t));
}

// orthogonal projection of a line l (a bivector) onto a plane t (a trivector)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> project_onto(BiVec4ds<T> const& B,
                                                          TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(rwdg(t, wdg(B, r_weight_dual(t)))); // ortho_proj4ds
}


////////////////////////////////////////////////////////////////////////////////
// Projections for 4ds: (HINT: unitize after projection, if not at infinity)
//
// ortho_proj4ds(a, b)     = rwdg(b, r_weight_expand4ds(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' contained in b)
// REQUIRES: gr(a) < gr(b)
//
//
// central_proj4ds(a, b)   = rwdg(b, r_bulk_expand4ds(a, b) )
// (a projected centrally towards origin onto b, effectively creating a new a'
// contained in b)
// REQUIRES: gr(a) < gr(b)
//
// ortho_antiproj4ds(a, b) = wdg(b, r_weight_contract4ds(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' containing b)
// REQUIRES: gr(a) > gr(b)
//
////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2> decltype(auto) ortho_proj4ds(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!

    // project the smaller grade object onto to larger grade object
    auto p = rwdg(std::forward<arg2>(b),
                  r_weight_expand4ds(std::forward<arg1>(a), std::forward<arg2>(b)));

    // return a unitized object, if it is not located in the horizon
    auto nrm_sq = weight_nrm_sq(p);
    if ((nrm_sq > eps) && (nrm_sq != 1.0)) {
        p = unitize(p);
    }

    return p;
}

template <typename arg1, typename arg2> decltype(auto) central_proj4ds(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!

    // project the smaller grade object onto to larger grade object
    auto p = rwdg(std::forward<arg2>(b),
                  r_bulk_expand4ds(std::forward<arg1>(a), std::forward<arg2>(b)));

    // return a unitized object, if it is not located in the horizon
    auto nrm_sq = weight_nrm_sq(p);
    if ((nrm_sq > eps) && (nrm_sq != 1.0)) {
        p = unitize(p);
    }

    return p;
}

template <typename arg1, typename arg2>
decltype(auto) ortho_antiproj4ds(arg1&& a, arg2&& b)
{
    return wdg(std::forward<arg2>(b),
               r_weight_contract4ds(std::forward<arg1>(a), std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// reflections of points, lines, and planes on planes
////////////////////////////////////////////////////////////////////////////////

// reflect a vector u in an arbitrary trivector, i.e. a plane t
// pre-condition: t must be unitized
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reflect_on(Vec4ds<T> const& v,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(gr1(rgpr(rgpr(t, v), rrev(t))));
}

// reflect a bivector B (a line) in an arbitrary trivector t
// pre-condition: t must be unitized
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> reflect_on(BiVec4ds<T> const& B,
                                                        TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(-gr2(rgpr(rgpr(t, B), rrev(t))));
}

// reflect a trivector t1 (a plane) in an arbitrary trivector t2 (a unitized plane)
// pre-condition: t2 must be unitized
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> reflect_on(TriVec4ds<T> const& t1,
                                                         TriVec4ds<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec4ds<ctype>(-gr3(rgpr(rgpr(t2, t1), rrev(t2))));
}

////////////////////////////////////////////////////////////////////////////////
// inversion
// (=point reflection resuling from reflecting across three perpendicular planes)
// the resulting minus comes for the 3 consecutive reflections
////////////////////////////////////////////////////////////////////////////////

// (point-)reflect a point q in an arbitrary point p
// pre-condition: p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> invert_on(Vec4ds<T> const& q,
                                                     Vec4ds<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(-gr1(rgpr(rgpr(p, q), rrev(p))));
}

// (point-)reflect a line l in an arbitrary point p
// pre-condition: p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> invert_on(BiVec4ds<T> const& l,
                                                       Vec4ds<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(-gr2(rgpr(rgpr(p, l), rrev(p))));
}

// (point-)reflect a plane t in an arbitrary point p
// pre-condition: p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> invert_on(TriVec4ds<T> const& t,
                                                        Vec4ds<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec4ds<ctype>(-gr3(rgpr(rgpr(p, t), rrev(p))));
}


////////////////////////////////////////////////////////////////////////////////
// 4ds euclidean distance
////////////////////////////////////////////////////////////////////////////////

// returns the euclidean distance between objects as homogeneous magnitude
template <typename arg1, typename arg2>
constexpr DualNum4ds<value_t> dist4ds(arg1&& a, arg2&& b)
{
    if constexpr (gr(std::decay_t<arg1>{}) + gr(std::decay_t<arg2>{}) == 4) {
        return DualNum4ds<value_t>(rwdg(a, b), weight_nrm(wdg(a, att(b))));
    }
    else {
        return DualNum4ds<value_t>(bulk_nrm(att(wdg(a, b))), weight_nrm(wdg(a, att(b))));
    }
}

////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar4ds<T> a, Scalar4ds<U> b, value_t tolerance = eps)
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
bool is_congruent(Vec4ds<T> const& a, Vec4ds<U> const& b, value_t tolerance = eps)
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
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec4ds<T> const& a, BiVec4ds<U> const& b, value_t tolerance = eps)
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
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(TriVec4ds<T> const& a, TriVec4ds<U> const& b, value_t tolerance = eps)
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

// For pseudoscalars: all non-zero pseudoscalars in 4ds represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar4ds<T> a, PScalar4ds<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars in 4ds are congruent (represent the full 4D space)
    return true;
}

} // namespace hd::ga::sta