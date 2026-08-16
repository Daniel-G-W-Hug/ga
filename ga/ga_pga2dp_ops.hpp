#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_pga2dp_ops_basics.hpp"
#include "ga_pga2dp_ops_products.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides functionality that is based on pga2dp ops basics and products:
//
// - angle()                             -> angle operations
// - rexp()                               -> exponential (w.r.t. rgpr)
// - rlog()                               -> logarithm (w.r.t. rgpr, inverse of rexp)
// - rsqrt(M)                             -> sqrt of a motor (w.r.t. rgpr)
// - get_motor()                         -> provide a motor from (point, phi), or (delta),
// - get_motor_from_lines()              -> provide a motor from (from two lines moved
//                                                                into each other)
// - move2dp(), move2dp_opt()            -> move object with motor
// - project_onto(), reject_from()       -> simple projection and rejection
// - expand()                            -> expansion: new line through point
//                                                     perpendicular to line
// - ortho_proj2dp()                     -> orthogonal projection onto object
// - central_proj2dp()                   -> central projection towards origin onto object
// - ortho_antiproj2dp()                 -> orthogonal antiprojection onto object
// - reflect_on()                        -> reflections
// - invert_on()                         -> inversions
// - sup()                               -> point on line that is nearest to origin
// - att()                               -> object attitude
// - dist2dp()                           -> Euclidean distance and homogeneous magnitude
//
// - is_congruent()                      -> Same up to a scalar factor (is same subspace)
// - is_close()                          -> Same value within a RELATIVE tolerance
// - is_same_motion()                    -> Do two motors describe the same rigid
//                                          motion? (motors double-cover them)
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// angle operations 2dp
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors, i.e. directions to points at infinity
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> angle(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    // angle is defined only between directions towards infinity
    // thus only the bulk part is considered
    // (allows angle measurement between points and directions)

    ctype nrm_prod = ctype(bulk_nrm(v1)) * ctype(bulk_nrm(v2));
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between two bivectors, i.e. between two lines
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::common_type_t<T, U> angle(BiVec2dp<T> const& B1, BiVec2dp<U> const& B2)
{

    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(r_weight_contract2dp(B1, B2));
    // hint: weight_nrm returns a pscalar2dp! ctype() required around each single result,
    // otherwise geometric product between pseudoscalars evaluated, which evaluates to
    // zero
    ctype nrm_prod = ctype(weight_nrm(B1)) * ctype(weight_nrm(B2));
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}


////////////////////////////////////////////////////////////////////////////////
// rexp(motor) and rsqrt(motor) with respect to the regressive geometric product rgpr()
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dp_U<T> rexp(Vec2dp<T> const& arg)
{
    auto w_sq = weight_nrm_sq(arg); // rotation angle^2
    if (w_sq == T(0.0)) {
        // pure translation

        // arg = att(bulk_dual(delta)) to move in direction of vector delta
        return MVec2dp_U<T>(arg, PScalar2dp<T>(1.0));
    }

    // rotation with angle != 0
    T phi = arg.z;                   // arg.z == rotation angle
    T scale = std::sin(phi) / arg.z; // unitize arg and multiply by sin(phi)
    return unitize(MVec2dp_U<T>(arg * scale, PScalar2dp<T>(std::cos(phi))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dp_U<T> rsqrt(MVec2dp_U<T> const& M)
{
    return unitize(M + PScalar2dp<T>(1.0));
}

////////////////////////////////////////////////////////////////////////////////
// rlog() with respect to rgpr() -- the inverse of rexp(): for a (unit) motor M it returns
// the vector generator arg with rexp(arg) == M.
//
// A 2dp motor is M = p*sin(phi) + e321*cos(phi) (rotation by phi about the point p), so
// for a unit motor the e3 (weight) component is sin(phi) and the pseudoscalar is
// cos(phi). rexp() builds the rotation branch as arg * sin(phi)/arg.z with arg.z = phi,
// so the inverse is arg = M_vec * (phi / sin(phi)) with arg.z = phi. A pure translation
// (phi ~ 0) has arg.z = 0 and M = (arg, 1), so arg is just the vector part.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dp<T> rlog(MVec2dp_U<T> const& M_in)
{
    auto const M = unitize(M_in);
    T const vz = M.c2;                // e3 (weight) component = sin(phi)
    T const ps = M.c3;                // pseudoscalar          = cos(phi)
    T const phi = std::atan2(vz, ps); // rotation angle
    T const s = std::sin(phi);
    if (std::abs(s) < eps) {
        // phi ~ 0: pure translation (or identity). rexp(arg) with arg.z = 0 gives (arg,
        // 1).
        return Vec2dp<T>(M.c0, M.c1, T(0.0));
    }
    T const scale = phi / s; // inverse of rexp()'s sin(phi)/phi
    return Vec2dp<T>(M.c0 * scale, M.c1 * scale, phi);
}


////////////////////////////////////////////////////////////////////////////////
// 2dp motor operations (translation and rotation)
//
// Every motor in pga2dp is an odd-grade multivector MVec2dp_U.
//
// A proper isometry in 2dp has a fixed point p = p.x e1 + p.y e2 + p.z e3
// around which a rotation occurs with an angle phi.
//
// Every rotating motor has the form: M = p sin(phi) + e321 cos(phi).
//
// Translations can be covered by a rotation around a point at infinity.
//
// Combined rotations and translations must be created by concatenating the motion
// operations in the sequence rgpr(motor_applied_last, motor_applied_first)
////////////////////////////////////////////////////////////////////////////////

// create a (unitized) motor from a fixed point and a turning angle
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> get_motor(Vec2dp<T> const& P, U theta)
{
    // point p must be unitized to avoid suprises
    auto P_weight_nrm_sq = weight_nrm_sq(P);
    if (P_weight_nrm_sq == 0.0) {
        throw std::invalid_argument(
            "get_motor: Cannot use ideal points P with P.z == 0.0");
    }
    auto Pn{P};
    if ((P_weight_nrm_sq > eps) && (P_weight_nrm_sq != 1.0)) {
        Pn = unitize(P);
    }

    using ctype = std::common_type_t<T, U>;
    ctype half_angle = 0.5 * theta;
    return unitize(MVec2dp_U<ctype>(Vec2dp<ctype>(Pn * std::sin(half_angle)),
                                    PScalar2dp<ctype>(std::cos(half_angle))));
}

// Create a translation motor from a translation delta = vec2dp(delta_x, delta_y, 0)
// Move in direction and by length of translation vector (length = bulk_nrm(delta))
// ATTENTION: The translation is assumed to be a direction, i.e. with z-component == 0.
//            The z-component is ignored and only the x- and y-components are used.
//            Due to the application via the regressive sandwich product the vector in
//            the motor needs to be multiplied by 0.5
//
// The function implements a motor calculated from exp_{\veedot}(P)
// with the vector P = 0.5 * bulk_dual(delta) \vee H_{2dp} = 0.5 * att(bulk_dual(delta))
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dp_U<T> get_motor(Vec2dp<T> const& delta)
{
    return MVec2dp_U<T>(0.5 * Vec2dp<T>(-delta.y, delta.x, T(0.0)), PScalar2dp<T>(1.0));
}

// overload of translation motor for a Vector2d:
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dp_U<T> get_motor(Vector2d<T> const& delta)
{
    return MVec2dp_U<T>(0.5 * Vec2dp<T>(-delta.y, delta.x, T(0.0)), PScalar2dp<T>(1.0));
}

// create a (unitized) motor directly from two consequtive reflections across two
// (potentially intersecting) lines
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> get_motor_from_lines(BiVec2dp<T> const& B1,
                                                                   BiVec2dp<U> const& B2)
{
    // take lines (=bivectors) as input and return a motor M
    // 1st apply reflection across line B1, then across B2 to get a motor that rotates
    // around the intersection point of lines B1 and B2
    // or translates if the lines are parallel
    //
    // for use of motor M either directly on object u:
    //
    //     auto v_moved = gr1( rgpr(rgpr(M, v), rrev(M)) );
    // or
    //     auto B_moved = gr2( rgpr(rgpr(M, B), rrev(M)) );
    // or
    //     auto v_moved = move2dp(v,M);  // moves v according to the motor M
    //     auto B_moved = move2dp(B,M);  // moves B according to the motor M

    // lines B1 and B2 need to be unitized to avoid surprises
    auto w_nrm_sq = weight_nrm_sq(B1);
    auto Bu1{B1};
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        Bu1 = unitize(Bu1);
    }
    w_nrm_sq = weight_nrm_sq(B2);
    auto Bu2{B2};
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        Bu2 = unitize(Bu2);
    }

    auto M = rgpr(Bu2, Bu1);
    w_nrm_sq = weight_nrm_sq(M);
    if ((w_nrm_sq > eps) && (w_nrm_sq != 1.0)) {
        M = unitize(M);
    }
    return M; // based on the regressive geometric product
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> move2dp(Vec2dp<T> const& v,
                                                   MVec2dp_U<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves v (a vector representing a projective point) according to the motor M
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(gr1(rgpr(rgpr(M, v), rrev(M))));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> move2dp(BiVec2dp<T> const& B,
                                                     MVec2dp_U<U> const& M)
{
    // motor M must be unitized to avoid surprises

    // moves B (a bivector representing a line) according to the motor M
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(rgpr(rgpr(M, B), rrev(M))));
}

// rotate a motor (required for robotics applications using kinematic chains)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_U<std::common_type_t<T, U>> move2dp(MVec2dp_U<T> const& M_orig,
                                                      MVec2dp_U<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // moves M_orig (a motor) according to the motor M
    // e.g. kinematic chains in robotics application with coupled joints
    // effectively "rotating the rotation direction"
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_U<ctype>(rgpr(rgpr(M, M_orig), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp<std::common_type_t<T, U>> move2dp(MVec2dp<T> const& MV,
                                                    MVec2dp_U<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // rotate one multivector MV with motor M
    // (only rotates the vector and bivector components of the 2dp multivector)
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(rgpr(rgpr(M, MV), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> move2dp_opt(Vec2dp<T> const& v,
                                                       MVec2dp_U<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // move one vector v with motor M (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga2dp regressive sandwich product)
    ctype const h0 = M.c2 * M.c2;
    ctype const h1 = M.c3 * M.c3;
    ctype const h2 = 2.0 * M.c2 * M.c3;

    ctype const k11 = -h0 + h1;
    ctype const k12 = -h2;
    ctype const k13 = 2.0 * (M.c0 * M.c2 + M.c1 * M.c3);

    ctype const k21 = h2;
    ctype const k22 = k11;
    ctype const k23 = 2.0 * (-M.c0 * M.c3 + M.c1 * M.c2);

    ctype const k33 = h0 + h1;

    return Vec2dp<ctype>(k11 * v.x + k12 * v.y + k13 * v.z,
                         k21 * v.x + k22 * v.y + k23 * v.z, k33 * v.z);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::vector<Vec2dp<std::common_type_t<T, U>>>
move2dp(std::vector<Vec2dp<T>> const& vec, MVec2dp_U<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // move many vectors vec(v) with motor M (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga2dp regressive sandwich product)
    ctype const h0 = M.c2 * M.c2;
    ctype const h1 = M.c3 * M.c3;
    ctype const h2 = 2.0 * M.c2 * M.c3;

    ctype const k11 = -h0 + h1;
    ctype const k12 = -h2;
    ctype const k13 = 2.0 * (M.c0 * M.c2 + M.c1 * M.c3);

    ctype const k21 = h2;
    ctype const k22 = k11;
    ctype const k23 = 2.0 * (-M.c0 * M.c3 + M.c1 * M.c2);

    ctype const k33 = h0 + h1;

    std::vector<Vec2dp<ctype>> result;
    result.reserve(vec.size());

    for (auto const& v : vec) {
        result.emplace_back(Vec2dp<ctype>(k11 * v.x + k12 * v.y + k13 * v.z,
                                          k21 * v.x + k22 * v.y + k23 * v.z, k33 * v.z));
    }
    return result;
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> move2dp_opt(BiVec2dp<T> const& B,
                                                         MVec2dp_U<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // move one bivector B with motor M (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga2dp regressive sandwich product)
    ctype const h0 = M.c2 * M.c2;
    ctype const h1 = M.c3 * M.c3;
    ctype const h2 = 2.0 * M.c2 * M.c3;

    ctype const k11 = -h0 + h1;
    ctype const k12 = -h2;
    ctype const k21 = h2;
    ctype const k22 = k11;

    ctype const k31 = 2.0 * (M.c0 * M.c3 + M.c1 * M.c2);
    ctype const k32 = 2.0 * (-M.c0 * M.c2 + M.c1 * M.c3);
    ctype const k33 = h0 + h1;

    return BiVec2dp<ctype>(k11 * B.x + k12 * B.y, k21 * B.x + k22 * B.y,
                           k31 * B.x + k32 * B.y + k33 * B.z);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr std::vector<BiVec2dp<std::common_type_t<T, U>>>
move2dp(std::vector<BiVec2dp<T>> const& bvec, MVec2dp_U<U> const& M)
{
    // pre: motor M must be unitized to avoid surprises

    // move many bivectors bvec(B) with motor M (optimized)
    using ctype = std::common_type_t<T, U>;

    // coefficients calculated with ga_prdxpr (pga2dp regressive sandwich product)
    ctype const h0 = M.c2 * M.c2;
    ctype const h1 = M.c3 * M.c3;
    ctype const h2 = 2.0 * M.c2 * M.c3;

    ctype const k11 = -h0 + h1;
    ctype const k12 = -h2;
    ctype const k21 = h2;
    ctype const k22 = k11;

    ctype const k31 = 2.0 * (M.c0 * M.c3 + M.c1 * M.c2);
    ctype const k32 = 2.0 * (-M.c0 * M.c2 + M.c1 * M.c3);
    ctype const k33 = h0 + h1;

    std::vector<BiVec2dp<ctype>> result;
    result.reserve(bvec.size());

    for (auto const& B : bvec) {
        result.emplace_back(BiVec2dp<ctype>(k11 * B.x + k12 * B.y, k21 * B.x + k22 * B.y,
                                            k31 * B.x + k32 * B.y + k33 * B.z));
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
constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v1,
                                                        Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of vector v1 from a vector v2
// returns component of v1 perpendicular to v2
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v1,
                                                       Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec2dp<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}


// projection of a vector v onto a bivector B (a line)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(rwdg(B, wdg(v, weight_dual(B)))); // ortho_proj
}

// rejection of vector v from a bivector B (a line)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v,
                                                       BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v - project_onto(v, B));
}

// expand to a new line with goes through point p and is perpendicular to line l
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Line2d<std::common_type_t<T, U>> expand(Point2d<T> const& p, Line2d<U> const& l)
{
    return r_weight_expand2dp(Point2dp<std::common_type_t<T, U>>(p), l);
}


////////////////////////////////////////////////////////////////////////////////
// Projections for 2dp: (HINT: unitize after projection, if not at infinity)
//
// ortho_proj2dp(a, b)     = rwdg(b, r_weight_expand2dp(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' contained in b)
// REQUIRES: gr(a) < gr(b)
//
// central_proj2dp(a, b)   = rwdg(b, r_bulk_expand2dp(a, b) )
// (a projected centrally towards origin onto b, effectively creating a new a'
// contained in b)
// REQUIRES: gr(a) < gr(b)
//
// ortho_antiproj2dp(a, b) = wdg(b, r_weight_contract2dp(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' containing b)
// REQUIRES: gr(a) > gr(b)
//
////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2> decltype(auto) ortho_proj2dp(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!

    // project the smaller grade object onto to larger grade object
    auto p = rwdg(std::forward<arg2>(b),
                  r_weight_expand2dp(std::forward<arg1>(a), std::forward<arg2>(b)));

    // return a unitized object, if it is not located in the horizon
    auto nrm_sq = weight_nrm_sq(p);
    if ((nrm_sq > eps) && (nrm_sq != 1.0)) {
        p = unitize(p);
    }

    return p;
}

template <typename arg1, typename arg2> decltype(auto) central_proj2dp(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!

    // project the smaller grade object onto to larger grade object
    auto p = rwdg(std::forward<arg2>(b),
                  r_bulk_expand2dp(std::forward<arg1>(a), std::forward<arg2>(b)));

    // return a unitized object, if it is not located in the horizon
    auto nrm_sq = weight_nrm_sq(p);
    if ((nrm_sq > eps) && (nrm_sq != 1.0)) {
        p = unitize(p);
    }

    return p;
}

template <typename arg1, typename arg2>
decltype(auto) ortho_antiproj2dp(arg1&& a, arg2&& b)
{
    return wdg(std::forward<arg2>(b),
               r_weight_contract2dp(std::forward<arg1>(a), std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// reflections
//
// see MacDonald, "Linear and geometric algebra", p. 129, theorem 7.10:
// The j-blade U reflected in the k-dimensional subspace B is
// M_B(U) = (-1)^[j*(k+1)] B ⟇ U ⟇ B^(-1) = (-1)^[j*(k+1)] B ⟇ U ⟇ rrev(B)
////////////////////////////////////////////////////////////////////////////////

// reflect a vector u in an arbitrary bivector, i.e. a line
// pre-condition: B must be unitized, or the object will be scaled as well!
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> reflect_on(Vec2dp<T> const& v,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(-gr1(rgpr(rgpr(B, v), rrev(B))));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling lines)
// pre-condition: B must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> reflect_on(BiVec2dp<T> const& UB,
                                                        BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(rgpr(rgpr(B, UB), rrev(B))));
}


////////////////////////////////////////////////////////////////////////////////
// inversion
// (=point reflection resuling from reflecting across two perpendicular lines)
// the resulting positive sign comes from the two consecutive reflections
////////////////////////////////////////////////////////////////////////////////

// (point-)reflect a point q in an arbitrary point p
// pre-condition: p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dp<std::common_type_t<T, U>> invert_on(Vec2dp<T> const& q,
                                                     Vec2dp<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(gr1(rgpr(rgpr(p, q), rrev(p))));
}

// (point-)reflect a line l in an arbitrary point p
// pre-condition: p must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dp<std::common_type_t<T, U>> invert_on(BiVec2dp<T> const& l,
                                                       Vec2dp<U> const& p)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(rgpr(rgpr(p, l), rrev(p))));
}


////////////////////////////////////////////////////////////////////////////////
// support: nearest point to the origin (by projecting origin onto line)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
Vec2dp<T> sup(BiVec2dp<T> const& B)
{
    // REQUIRES: a line (BiVec2dp) as argument

    // project origin onto line
    return ortho_proj2dp(O_2dp, B);
}


////////////////////////////////////////////////////////////////////////////////
// attitude operations: att = rwdg( u, cmpl(e3_2dp) ) = rwdg(u, H_2dp)
//
// (the attitude is the intersection of the object with the horizon)
// the result of att(object_with_grade_k) is an object with grade k-1
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dp<T> att(Vec2dp<T> const& v)
{
    return Scalar2dp<T>(v.z);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dp<T> att(Point2dp<T> const& p)
{
    return Scalar2dp<T>(p.z);
}

// return the attitude (i.e. the direction vector) of the line
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dp<T> att(BiVec2dp<T> const& B)
{
    return Vec2dp<T>(B.x, B.y, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dp<T> att(Line2d<T> const& l)
{
    return Vec2dp<T>(l.x, l.y, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dp<T> att(PScalar2dp<T> ps)
{
    return BiVec2dp<T>(T(0.0), T(0.0), T(ps));
}


////////////////////////////////////////////////////////////////////////////////
// 2dp euclidean distance
////////////////////////////////////////////////////////////////////////////////

// returns the euclidean distance between objects as homogeneous magnitude
template <typename arg1, typename arg2> DualNum2dp<value_t> dist2dp(arg1&& a, arg2&& b)
{
    if constexpr (gr(std::decay_t<arg1>{}) + gr(std::decay_t<arg2>{}) == 3) {
        return DualNum2dp<value_t>(rwdg(a, b), weight_nrm(wdg(a, att(b))));
    }
    else {
        return DualNum2dp<value_t>(bulk_nrm(att(wdg(a, b))), weight_nrm(wdg(a, att(b))));
    }
}

////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar2dp<T> a, Scalar2dp<U> b, value_t tolerance = eps)
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
bool is_congruent(Vec2dp<T> const& a, Vec2dp<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<3>({value_t(a.x), value_t(a.y), value_t(a.z)},
                                       {value_t(b.x), value_t(b.y), value_t(b.z)},
                                       tolerance);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec2dp<T> const& a, BiVec2dp<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<3>({value_t(a.x), value_t(a.y), value_t(a.z)},
                                       {value_t(b.x), value_t(b.y), value_t(b.z)},
                                       tolerance);
}

// For pseudoscalars: all non-zero pseudoscalars in 2DP represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar2dp<T> a, PScalar2dp<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars in 2DP are congruent (represent the full 3D space)
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// is_close(): equality within a RELATIVE tolerance
////////////////////////////////////////////////////////////////////////////////
//
// Same question as operator==, but with the threshold scaled by the operands (see
// detail::coeffs_close). Reach for it whenever the values carry a physical scale:
// operator== measures against an absolute eps, which cannot resolve anything once
// coordinates grow large, where a single ulp already exceeds it.
//
// Distinct from is_congruent(), which allows an arbitrary scale factor between the
// operands and so answers "same subspace", not "same value".

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Scalar2dp<T> a, Scalar2dp<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Vec2dp<T> const& a, Vec2dp<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<3>({value_t(a.x), value_t(a.y), value_t(a.z)},
                                   {value_t(b.x), value_t(b.y), value_t(b.z)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(BiVec2dp<T> const& a, BiVec2dp<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<3>({value_t(a.x), value_t(a.y), value_t(a.z)},
                                   {value_t(b.x), value_t(b.y), value_t(b.z)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(PScalar2dp<T> a, PScalar2dp<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec2dp_U<T> const& a, MVec2dp_U<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<4>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3)}, rel_tol);
}


////////////////////////////////////////////////////////////////////////////////
// is_same_motion(): do two motors describe the same rigid motion?
////////////////////////////////////////////////////////////////////////////////
//
// NOT the same question as M1 == M2. Motors DOUBLE-COVER the rigid motions: the sandwich
//
//     X -> M (x) X (x) rrev(M)
//
// is quadratic in M, so M and -M move every object identically while comparing unequal
// component by component. A rotation by 2*pi returns -1 rather than 1 for exactly this
// reason, and a motor obtained by a different route (a log/exp round trip, a chain of
// compositions, a frame tree walk) routinely comes back negated. Testing motors with
// operator== therefore asks about the REPRESENTATION; this asks about the MOTION.
//
// It is decided by what the motors DO -- where they send the origin and the two axis
// directions, which together pin a planar rigid motion uniquely -- rather than by
// comparing coefficients up to a sign. That keeps the answer right for motors that are
// not unit: a scaled motor leaves points where they were (the sandwich scales by |M|^2,
// which unitize divides out again) but stretches directions, so it is NOT the same
// motion, and the direction comparisons below say so.
//
// Tolerances are relative (is_close), since the origin's image carries whatever physical
// scale the translation has while the axis images are of order one.
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_same_motion(MVec2dp_U<T> const& M1, MVec2dp_U<U> const& M2,
                    value_t rel_tol = eps_congruent)
{
    if (!is_close(unitize(move2dp(O_2dp, M1)), unitize(move2dp(O_2dp, M2)), rel_tol)) {
        return false;
    }
    return is_close(move2dp(x_dir_2dp, M1), move2dp(x_dir_2dp, M2), rel_tol) &&
           is_close(move2dp(y_dir_2dp, M1), move2dp(y_dir_2dp, M2), rel_tol);
}

} // namespace hd::ga::pga