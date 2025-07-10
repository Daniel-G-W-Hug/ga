#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_pga2dp_ops_basics.hpp"
#include "ga_pga2dp_ops_products.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides functionality that is based on pga2dp ops basics and products:
//
// - angle()                        -> angle operations
// - TODO: exp()                          -> exponential function
// - motor(), motor_from_lines()    -> provide a motor
// - move2dp(), move2dp_opt()       -> move object with motor
// - project_onto(), reject_from()  -> simple projection and rejection
// - expand()                       -> expansion: new line through point perpend. to line
// - ortho_proj2dp()                -> orthogonal projection onto object
// - central_proj2dp()              -> central projection towards origin onto object
// - ortho_antiproj2dp()            -> orthogonal antiprojection onto object
// - reflect_on()                   -> reflections
// - support2dp()                   -> point on line that is nearest to origin
// - att                            -> object attitude
// - dist2dp()                      -> Euclidean distance and homogeneous magnitude
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// angle operations 2dp
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors, i.e. directions to points at infinity
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    if ((weight_nrm_sq(v1) != 0.0) || (weight_nrm_sq(v2) != 0.0)) {
        // the angle between points not at infinity or points not at infinity and a
        // direction towards infinity is defined as zero
        return 0.0;
    }

    // angle is defined only between directions towards points at infinity

    ctype nrm_prod = bulk_nrm(v1) * bulk_nrm(v2);
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between two bivectors, i.e. between two lines
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{

    using ctype = std::common_type_t<T, U>;
    ctype contr = right_weight_contract2dp(B1, B2);
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
// 2dp motor operations (translation and rotation)
//
// Every motor in pga2dp is an odd grade multivector MVec2dp_U.
//
// A proper isometry in 2dp has a fixed point p = p.x e1 + p.y e2 + p.z e3
// around which a rotation occurs with an angle phi.
//
// Every motor has the form: M = p sin(phi) + e321 cos(phi).
//
// (derived from the exponential function with respect to the regressive
// geometric product)
////////////////////////////////////////////////////////////////////////////////

// create a (unitized) motor from a fixed point and a turning angle
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>> motor(Vec2dp<T> const& p, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = 0.5 * theta;
    return unitize(MVec2dp_U<ctype>(Vec2dp<ctype>(p * std::sin(half_angle)),
                                    PScalar2dp<ctype>(std::cos(half_angle))));
}

// create a translation motor from a direction vector (given as a Vec2dp)
// move in direction and by length of direction vector (lenght = its bulk_nrm)
// ATTENTION: the direction is assumed to be a direction vector, i.e. with z == 0
//            the z-component is ignored and only the x- and y-components are used
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> motor(Vec2dp<T> const& direction)
{
    return MVec2dp_U<T>(0.5 * Vec2dp<T>(-direction.y, direction.x, T(0.0)),
                        PScalar2dp<T>(1.0));
}

// create a translation motor from a direction vector (given as a Vector2d)
// move in direction and by length of direction vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> motor(Vector2d<T> const& direction)
{
    return MVec2dp_U<T>(0.5 * Vec2dp<T>(-direction.y, direction.x, T(0.0)),
                        PScalar2dp<T>(1.0));
}

// create a (unitized) motor directly from two (potentially intersecting) lines
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_U<std::common_type_t<T, U>>
motor_from_lines(BiVec2dp<T> const& B1, BiVec2dp<U> const& B2)
{
    // take lines (=bivectors) as input and return a motor R
    // 1st apply reflection across line B1, then across B2 to get a motor that rotates
    // (or translates) around the intersection point of lines B1 and B2
    //
    // for use of motor R either directly on object u (inefficient):
    //     auto v_moved = gr1( rgpr(rgpr(R, v), rrev(R)) );
    // or
    //     auto B_moved = gr2( rgpr(rgpr(R, B), rrev(R)) );
    // or
    //                                   // optimized for reduced effort
    //     auto v_moved = move2dp(v,R);  // moves v according to the motor R
    //     auto B_moved = move2dp(B,R);  // moves B according to the motor R
    //
    return unitize(rgpr(B2, B1)); // based on the regressive geometric product
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp(Vec2dp<T> const& v,
                                                          MVec2dp_U<U> const& R)
{
    // assumes: motor R is unitized

    // moves v (a vector representing a projective point) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(gr1(rgpr(rgpr(R, v), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp(BiVec2dp<T> const& B,
                                                            MVec2dp_U<U> const& R)
{
    // assumes: motor R is unitized

    // moves B (a bivector representing a line) according to the motor R
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(rgpr(rgpr(R, B), rrev(R))));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp_opt(Vec2dp<T> const& v,
                                                              MVec2dp_U<U> const& R)
{
    // moves v (a vector representing a projective point) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    //
    // (could potentially be further optimized by exporting the matrix-representation
    // if many transformations should be calculated using the same rotor as
    // v' = matrix * v)

    using ctype = std::common_type_t<T, U>;
    auto k02 = R.c0 * R.c2;
    auto k03 = R.c0 * R.c3;
    auto k12 = R.c1 * R.c2;
    auto k13 = R.c1 * R.c3;
    auto k22 = R.c2 * R.c2;
    auto k23 = R.c2 * R.c3;
    auto k33 = R.c3 * R.c3;
    return Vec2dp<ctype>((-k22 + k33) * v.x - 2.0 * k23 * v.y + 2.0 * (k02 + k13) * v.z,
                         2.0 * k23 * v.x + (-k22 + k33) * v.y + 2.0 * (-k03 + k12) * v.z,
                         (k22 + k33) * v.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp_opt(BiVec2dp<T> const& B,
                                                                MVec2dp_U<U> const& R)
{
    // moves B (a bivector representing a line) according to the motor R
    // optimized by avoiding non-required calculations vs. original version
    //
    // (could potentially be further optimized by exporting the matrix-representation
    // if many transformations should be calculated using the same rotor as
    // B' = matrix * B)

    using ctype = std::common_type_t<T, U>;
    auto k02 = R.c0 * R.c2;
    auto k03 = R.c0 * R.c3;
    auto k12 = R.c1 * R.c2;
    auto k13 = R.c1 * R.c3;
    auto k22 = R.c2 * R.c2;
    auto k23 = R.c2 * R.c3;
    auto k33 = R.c3 * R.c3;
    return BiVec2dp<ctype>(
        (-k22 + k33) * B.x + 2.0 * (-k23) * B.y, 2.0 * k23 * B.x + (-k22 + k33) * B.y,
        2.0 * (k02 - k13) * B.x + 2.0 * (k12 + k03) * B.y + (k22 + k33) * B.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> move2dp_opt2(Vec2dp<T> const& v,
                                                               MVec2dp_U<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = R.c1 * v.z - R.c2 * v.y + R.c3 * v.x;
    ctype k2 = -R.c0 * v.z + R.c2 * v.x + R.c3 * v.y;
    ctype k3 = R.c3 * v.z;
    ctype k4 = -R.c2 * v.z;
    return Vec2dp<ctype>(k1 * R.c3 - k2 * R.c2 + k3 * R.c1 - k4 * R.c0,
                         k1 * R.c2 + k2 * R.c3 - k3 * R.c0 - k4 * R.c1,
                         k3 * R.c3 - k4 * R.c2);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> move2dp_opt2(BiVec2dp<T> const& B,
                                                                 MVec2dp_U<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = -R.c0 * B.x - R.c1 * B.y - R.c2 * B.z;
    ctype k2 = -R.c2 * B.y + R.c3 * B.x;
    ctype k3 = R.c2 * B.x + R.c3 * B.y;
    ctype k4 = R.c0 * B.y - R.c1 * B.x + R.c3 * B.z;
    return BiVec2dp<ctype>(k2 * R.c3 - k3 * R.c2, k2 * R.c2 + k3 * R.c3,
                           -k1 * R.c2 - k2 * R.c1 + k3 * R.c0 + k4 * R.c3);
}


////////////////////////////////////////////////////////////////////////////////
// projections, rejections
////////////////////////////////////////////////////////////////////////////////

// TODO: check whether the vector-vector formulas make sense at all, since they model
//       the representational space an not the modelled subspace
//       => potentially only meaningful, when implemented based on type vector2dp

// projection of a vector v1 onto vector v2
// returns component of v1 parallel to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v1,
                                                               Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of vector v1 from a vector v2
// returns component of v1 perpendicular to v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v1,
                                                              Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec2dp<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}


// projection of a vector v onto a bivector B (a line)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> project_onto(Vec2dp<T> const& v,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(rwdg(B, wdg(v, weight_dual(B)))); // ortho_proj
}

// rejection of vector v from a bivector B (a line)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reject_from(Vec2dp<T> const& v,
                                                              BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(v - project_onto(v, B));
}


// expand to a new line with goes through point p and is perpendicular to line l
// => returns a line (aka a bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Line2d<std::common_type_t<T, U>> expand(Point2d<T> const& p,
                                                         Line2d<U> const& l)
{
    return right_weight_expand2dp(Point2dp<std::common_type_t<T, U>>(p), l);
}


////////////////////////////////////////////////////////////////////////////////
// Projections for 2dp:
//
// ortho_proj2dp(a, b)     = rwdg(b, right_weight_expand2dp(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' contained in b)
// REQUIRES: gr(a) < gr(b)
//
// central_proj2dp(a, b)   = rwdg(b, right_bulk_expand2dp(a, b) )
// (a projected centrally towards origin onto b, effectively creating a new a'
// contained in b)
// REQUIRES: gr(a) < gr(b)
//
// ortho_antiproj2dp(a, b) = wdg(b, right_weight_contract2dp(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' containing b)
// REQUIRES: gr(a) > gr(b)
//
////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2> decltype(auto) ortho_proj2dp(arg1&& a, arg2&& b)
{

    // REQUIRES: gr(a) < gr(b), or does not compile!
    // project the smaller grade object onto to larger grade object
    return rwdg(std::forward<arg2>(b),
                right_weight_expand2dp(std::forward<arg1>(a), std::forward<arg2>(b)));
}

template <typename arg1, typename arg2> decltype(auto) central_proj2dp(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!
    // project the smaller grade object onto to larger grade object
    return rwdg(std::forward<arg2>(b),
                right_bulk_expand2dp(std::forward<arg1>(a), std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) ortho_antiproj2dp(arg1&& a, arg2&& b)
{
    return wdg(std::forward<arg2>(b),
               right_weight_contract2dp(std::forward<arg1>(a), std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// reflections
////////////////////////////////////////////////////////////////////////////////

// reflect a vector u in an arbitrary bivector, i.e. a line
// B must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2dp<std::common_type_t<T, U>> reflect_on(Vec2dp<T> const& v,
                                                             BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2dp<ctype>(-gr1(rgpr(rgpr(B, v), B)));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling lines)
// B must be unitized, or object will be scaled as well!
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec2dp<std::common_type_t<T, U>> reflect_on(BiVec2dp<T> const& UB,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec2dp<ctype>(gr2(rgpr(rgpr(B, UB), B)));
}


////////////////////////////////////////////////////////////////////////////////
// support2dp: nearest point to the origin (by projecting origin onto line)
////////////////////////////////////////////////////////////////////////////////

template <typename arg1> decltype(auto) support2dp(arg1&& a)
{
    // REQUIRES: a line (BiVec2dp) as argument
    return ortho_proj2dp(origin_2dp, std::forward<arg1>(a));
}


////////////////////////////////////////////////////////////////////////////////
// attitude operations: att = rwdg( u, cmpl(e3_2dp) ) = rwdg(u, horizon_2dp)
//
// (the attitude is the intersection of the object with the horizon)
// the result of att(object_with_grade_k) is an object with grade k-1
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> att(Vec2dp<T> const& v)
{
    return Scalar2dp<T>(v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> att(Point2dp<T> const& p)
{
    return Scalar2dp<T>(p.z);
}

// return the attitude (i.e. the direction vector) of the line
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> att(BiVec2dp<T> const& B)
{
    return Vec2dp<T>(B.y, -B.x, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> att(Line2d<T> const& l)
{
    return Vec2dp<T>(l.y, -l.x, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> att(PScalar2dp<T> ps)
{
    return BiVec2dp<T>(T(0.0), T(0.0), ps);
}


////////////////////////////////////////////////////////////////////////////////
// 2dp euclidean distance
////////////////////////////////////////////////////////////////////////////////

// returns the euclidean distance between objects as homogeneous magnitude
template <typename arg1, typename arg2> DualNum2dp<value_t> dist2dp(arg1&& a, arg2&& b)
{
    if (gr(a) + gr(b) == 3) {
        return DualNum2dp<value_t>(rwdg(a, b), weight_nrm(wdg(a, att(b))));
    }
    else {
        return DualNum2dp<value_t>(bulk_nrm(att(wdg(a, b))), weight_nrm(wdg(a, att(b))));
    }
}

} // namespace hd::ga::pga