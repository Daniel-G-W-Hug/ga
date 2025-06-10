#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "type_t/ga_type2dp.hpp"

#include "ga_error_handling.hpp"

namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// define pga2dp specific functions and special geometric object alias definitions:
//
// - bulk(), weight()           - return bulk and weight parts of objects
// - bulk_nrm_sq, bulk_nrm      - return bulk norm
// - weight_nrm_sq, weight_nrm  - return weight norm
// - geom_nrm_sq, geom_nrm      - return geometric norm
// - unitization                - return unitized object (weight_nrm scaled to 1.0)
// - bulk_dual, weight_dual     - return dual of respective part of object
// - att                        - return object attitude
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// bulk: u_bulk = G u (with G as the metric)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk(Scalar2dp<T> s)
{
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> bulk(Vec2dp<T> const& v)
{
    return Vec2dp<T>(v.x, v.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> bulk(BiVec2dp<T> const& B)
{
    return BiVec2dp<T>(T(0.0), T(0.0), B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> bulk([[maybe_unused]] PScalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> bulk(MVec2dp_E<T> const& M)
{
    return MVec2dp_E<T>(bulk(gr0(M)), bulk(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> bulk(MVec2dp_U<T> const& M)
{
    return MVec2dp_U<T>(bulk(gr1(M)), bulk(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> bulk(MVec2dp<T> const& M)
{
    return MVec2dp<T>(bulk(gr0(M)), bulk(gr1(M)), bulk(gr2(M)), bulk(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight: u_weight = lcmpl( G rcmpl(u) ) = rG u
//         (with G as the metric and rG as the anti-metric as given by Lengyel)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> weight([[maybe_unused]] Scalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> weight(Vec2dp<T> const& v)
{
    return Vec2dp<T>(T(0.0), T(0.0), v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> weight(BiVec2dp<T> const& B)
{
    return BiVec2dp<T>(B.x, B.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight(PScalar2dp<T> ps)
{
    return PScalar2dp<T>(ps);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> weight(MVec2dp_E<T> const& M)
{
    return MVec2dp_E<T>(weight(gr0(M)), weight(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> weight(MVec2dp_U<T> const& M)
{
    return MVec2dp_U<T>(weight(gr1(M)), weight(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> weight(MVec2dp<T> const& M)
{
    return MVec2dp<T>(weight(gr0(M)), weight(gr1(M)), weight(gr2(M)), weight(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk norm
////////////////////////////////////////////////////////////////////////////////

// return squared bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm_sq(Scalar2dp<T> s)
{
    // ||bulk(s)||^2 = dot(s,s) = s^(T) G s = gr0(s*rev(s)) = s*s
    // using rev(s) = (-1)^[k(k-1)/2] s for a k-blade: 0-blade => rev(s) = s
    return Scalar2dp<T>(T(s) * T(s));
}

// return bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm(Scalar2dp<T> s)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(s)));
}

// return squared bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm_sq(Vec2dp<T> const& v)
{
    // ||bulk(v)||^2 = dot(v,v) = v^(T) G v = gr0(v*rev(v)) = v.x * v.x + v.y * v.y
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    return Scalar2dp<T>(v.x * v.x + v.y * v.y);
}

// return bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm(Vec2dp<T> const& v)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(v)));
}

// return squared bulk magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm_sq(BiVec2dp<T> const& B)
{
    // ||bulk(B)||^2 = dot(B,B) = B^(T) G B = gr0(B*rev(B)) = B.z * B.z
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    return Scalar2dp<T>(B.z * B.z);
}

// return magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm(BiVec2dp<T> const& B)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(B)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm_sq([[maybe_unused]] PScalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm([[maybe_unused]] PScalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm_sq(MVec2dp_E<T> const& M)
{
    return Scalar2dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm(MVec2dp_E<T> const& M)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm_sq(MVec2dp_U<T> const& M)
{
    return Scalar2dp<T>(bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm(MVec2dp_U<T> const& M)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm_sq(MVec2dp<T> const& M)
{
    return Scalar2dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr2(M)) +
                        bulk_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_nrm(MVec2dp<T> const& M)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

////////////////////////////////////////////////////////////////////////////////
// weight norm
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> weight_nrm_sq([[maybe_unused]] Scalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> weight_nrm([[maybe_unused]] Scalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

// return squared weight norm of vector
// ||weight(v)||^2 = cmpl( dot(cmpl(v), cmpl(v)) ) = rdot(v,v) = v.z * v.z
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm_sq(Vec2dp<T> const& v)
{
    return PScalar2dp<T>(v.z * v.z);
}

// return weigth norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm(Vec2dp<T> const& v)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(v)));
}

// return squared weight norm of bivector
// ||weight(B)||^2 = cmpl( dot(cmpl(B),cmpl(B)) ) = rdot(B, B) = B.x * B.x + B.y * B.y
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm_sq(BiVec2dp<T> const& B)
{
    return PScalar2dp<T>(B.x * B.x + B.y * B.y);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm(BiVec2dp<T> const& B)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(B)));
}

// return squared weight norm of pseudoscalar
// ||weight(ps)||^2 = cmpl( dot(cmpl(ps),cmpl(ps)) ) = rdot(ps, ps) = ps * ps
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm_sq(PScalar2dp<T> ps)
{
    return PScalar2dp<T>(T(ps) * T(ps));
}

// return weight norm of pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm(PScalar2dp<T> ps)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(ps)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm_sq(MVec2dp_E<T> const& M)
{
    return PScalar2dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm(MVec2dp_E<T> const& M)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm_sq(MVec2dp_U<T> const& M)
{
    return PScalar2dp<T>(weight_nrm_sq(gr1(M)) + weight_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm(MVec2dp_U<T> const& M)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm_sq(MVec2dp<T> const& M)
{
    return PScalar2dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr1(M)) +
                         weight_nrm_sq(gr2(M)) + weight_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_nrm(MVec2dp<T> const& M)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(M)));
}


////////////////////////////////////////////////////////////////////////////////
// geometric norm: (perpendicular) distance to the origin (distance = c0/c1),
//                 i.e. the value c0 after unitization
//
// returns a dual number, for correct handling of objects at infinity (c1 == 0)
//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm_sq(Scalar2dp<T> s)
{
    return DualNum2dp<T>(bulk_nrm_sq(s), weight_nrm_sq(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm(Scalar2dp<T> s)
{
    return DualNum2dp<T>(bulk_nrm(s), weight_nrm(s));
}

// provide the distance of the point from the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm_sq(Vec2dp<T> const& v)
{
    return DualNum2dp<T>(bulk_nrm_sq(v), weight_nrm_sq(v));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm(Vec2dp<T> const& v)
{
    return DualNum2dp<T>(bulk_nrm(v), weight_nrm(v));
}

// provide the perpendicular distance of the line to the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm_sq(BiVec2dp<T> const& B)
{
    return DualNum2dp<T>(bulk_nrm_sq(B), weight_nrm_sq(B));
}

// return geometric norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm(BiVec2dp<T> const& B)
{
    return DualNum2dp<T>(bulk_nrm(B), weight_nrm(B));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm_sq(PScalar2dp<T> ps)
{
    return DualNum2dp<T>(bulk_nrm_sq(ps), weight_nrm_sq(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm(PScalar2dp<T> ps)
{
    return DualNum2dp<T>(bulk_nrm(ps), weight_nrm(ps));
}

// even grade multivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm_sq(MVec2dp_E<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm(MVec2dp_E<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm(M), weight_nrm(M));
}

// odd grade multivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm_sq(MVec2dp_U<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm(MVec2dp_U<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm(M), weight_nrm(M));
}

// multivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm_sq(MVec2dp<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum2dp<T> geom_nrm(MVec2dp<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm(M), weight_nrm(M));
}


////////////////////////////////////////////////////////////////////////////////
// unitization operations
////////////////////////////////////////////////////////////////////////////////

// return a DualNum2dp unitized to weight_nrm == 1.0
// => if the dual number represents a homogeneous norm,
//    the scalar part represents the geometric norm the after unitization
template <typename T>
    requires(std::floating_point<T>)
inline DualNum2dp<T> unitize(DualNum2dp<T> const& D)
{
    T n = D.c1; // the pseudoscalar part is the weight_nrm part
    hd::ga::detail::check_unitization<T>(std::abs(n), "dual number (2dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * D;
}

// return a vector unitized to v.z == 1.0  (implies weight_nrm(v) = 1.0)
template <typename T>
    requires(std::floating_point<T>)
inline Vec2dp<T> unitize(Vec2dp<T> const& v)
{
    T n = v.z;
    hd::ga::detail::check_unitization<T>(std::abs(n), "vector (2dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec2dp<T>(v.x * inv, v.y * inv, T(1.0));
}

// return a bivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline BiVec2dp<T> unitize(BiVec2dp<T> const& B)
{
    T n = weight_nrm(B);
    hd::ga::detail::check_unitization<T>(n, "bivector (2dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * B;
}

// return an even grade multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline MVec2dp_E<T> unitize(MVec2dp_E<T> const& M)
{
    T n = weight_nrm(M);
    hd::ga::detail::check_unitization<T>(n, "even grade multivector (2dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}

// return an odd grade multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline MVec2dp_U<T> unitize(MVec2dp_U<T> const& M)
{
    T n = weight_nrm(M);
    hd::ga::detail::check_unitization<T>(n, "odd grade multivector (2dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}

// return a multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline MVec2dp<T> unitize(MVec2dp<T> const& M)
{
    T n = weight_nrm(M);
    hd::ga::detail::check_unitization<T>(n, "multivector (2dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}


template <typename T>
    requires std::floating_point<T>
inline Point2dp<T> unitize(Point2dp<T> const& p)
{
    hd::ga::detail::check_unitization<T>(std::abs(p.z), "Point2dp");
    T inv = T(1.0) / p.z;
    return Point2dp<T>(p.x * inv, p.y * inv, T(1.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline Line2d<T> unitize(Line2d<T> const& l)
{
    // unitization for a 2d bivector means std::sqrt(x^2 + y^2) = 1
    // i.e. unitization of the direction vector of the line
    T wn = weight_nrm(l);
    hd::ga::detail::check_unitization<T>(wn, "Line2d");
    T inv = T(1.0) / wn;
    return Line2d<T>(l.x * inv, l.y * inv, l.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// bulk_dual (=complement operation applied to the bulk)
//
// the duals are defined uniquely w.r.t. the outer product
//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> bulk_dual(Scalar2dp<T> s)
{
    return PScalar2dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> bulk_dual(Vec2dp<T> const& v)
{
    return BiVec2dp<T>(-v.x, -v.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> bulk_dual(BiVec2dp<T> const& B)
{
    return Vec2dp<T>(T(0.0), T(0.0), -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> bulk_dual([[maybe_unused]] PScalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> bulk_dual(MVec2dp_E<T> const& M)
{
    return MVec2dp_U<T>(bulk_dual(gr2(M)), bulk_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> bulk_dual(MVec2dp_U<T> const& M)
{
    return MVec2dp_E<T>(bulk_dual(gr3(M)), bulk_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> bulk_dual(MVec2dp<T> const& M)
{
    return MVec2dp<T>(bulk_dual(gr3(M)), bulk_dual(gr2(M)), bulk_dual(gr1(M)),
                      bulk_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight_dual (=complement operation applied to the weight)
//
// the duals are defined uniquely w.r.t. the outer product
//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2dp<T> weight_dual([[maybe_unused]] Scalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> weight_dual(Vec2dp<T> const& v)
{
    return BiVec2dp<T>(T(0.0), T(0.0), -v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> weight_dual(BiVec2dp<T> const& B)
{
    return Vec2dp<T>(-B.x, -B.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> weight_dual(PScalar2dp<T> ps)
{
    return Scalar2dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> weight_dual(MVec2dp_E<T> const& M)
{
    return MVec2dp_U<T>(weight_dual(gr2(M)), weight_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> weight_dual(MVec2dp_U<T> const& M)
{
    return MVec2dp_E<T>(weight_dual(gr3(M)), weight_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> weight_dual(MVec2dp<T> const& M)
{
    return MVec2dp<T>(weight_dual(gr3(M)), weight_dual(gr2(M)), weight_dual(gr1(M)),
                      weight_dual(gr0(M)));
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

} // namespace hd::ga::pga