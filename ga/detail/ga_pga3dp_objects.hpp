#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "type_t/ga_type3dp.hpp"

#include "ga_error_handling.hpp"

namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// define pga3dp specific functions and special geometric object alias definitions:
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
inline constexpr Scalar3dp<T> bulk(Scalar3dp<T> s)
{
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> bulk(Vec3dp<T> const& v)
{
    return Vec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> bulk(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), B.mx, B.my, B.mz);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> bulk(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), t.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> bulk([[maybe_unused]] PScalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> bulk(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(bulk(gr0(M)), bulk(gr2(M)), bulk(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> bulk(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(bulk(gr1(M)), bulk(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> bulk(MVec3dp<T> const& M)
{
    return MVec3dp<T>(bulk(gr0(M)), bulk(gr1(M)), bulk(gr2(M)), bulk(gr3(M)),
                      bulk(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight: u_weight = lcmpl( G rcmpl(u) ) = rG u
//         (with G as the metric and rG as the anti-metric as given by Lengyel)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> weight([[maybe_unused]] Scalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> weight(Vec3dp<T> const& v)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> weight(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(B.vx, B.vy, B.vz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> weight(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(t.x, t.y, t.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(ps);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> weight(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(weight(gr0(M)), weight(gr2(M)), weight(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> weight(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(weight(gr1(M)), weight(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> weight(MVec3dp<T> const& M)
{
    return MVec3dp<T>(weight(gr0(M)), weight(gr1(M)), weight(gr2(M)), weight(gr3(M)),
                      weight(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk norm
////////////////////////////////////////////////////////////////////////////////

// return squared bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(Scalar3dp<T> s)
{
    // ||bulk(s)||^2 = dot(s,s) = s^(T) G s = gr0(s*rev(s)) = s*s
    // using rev(s) = (-1)^[k(k-1)/2] s for a k-blade: 0-blade => rev(s) = s
    return Scalar3dp<T>(T(s) * T(s));
}

// return bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(Scalar3dp<T> s)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(s)));
}

// return squared bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(Vec3dp<T> const& v)
{
    // ||bulk(v)||^2 = dot(v,v) = v^(T) G v = gr0(v*rev(v))
    //               = v.x * v.x + v.y * v.y + v.z * v.z
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    return Scalar3dp<T>(v.x * v.x + v.y * v.y + v.z * v.z);
}

// return bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(Vec3dp<T> const& v)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(v)));
}

// return squared bulk magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(BiVec3dp<T> const& B)
{
    // ||bulk(B)||^2 = dot(B,B) = B^(T) G B = gr0(B*rev(B))
    //               = B.mx * B.mx + B.my * B.my + B.mz * B.mz
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    return Scalar3dp<T>(B.mx * B.mx + B.my * B.my + B.mz * B.mz);
}

// return magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(BiVec3dp<T> const& B)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(B)));
}

// return squared bulk norm of trivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(TriVec3dp<T> const& t)
{
    // ||bulk(t)||^2 = dot(t,t) = t^(T) G t = gr0(t*rev(t)) = t.w * t.w
    // using rev(t) = (-1)^[k(k-1)/2] t for a k-blade: 3-blade => rev(t) = -t
    return Scalar3dp<T>(t.w * t.w);
}

// return bulk norm of trivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(TriVec3dp<T> const& t)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(t)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(MVec3dp_E<T> const& M)
{
    return Scalar3dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr2(M)) + bulk_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(MVec3dp_E<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(MVec3dp_U<T> const& M)
{
    return Scalar3dp<T>(bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(MVec3dp_U<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(MVec3dp<T> const& M)
{
    return Scalar3dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr2(M)) +
                        bulk_nrm_sq(gr3(M)) + bulk_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(MVec3dp<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight norm
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

// return squared weight norm of vector
// ||weight(v)||^2 = lcmpl( dot(rcmpl(v), rcmpl(v)) ) = rdot(v,v) = v.w * v.w
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(Vec3dp<T> const& v)
{
    return PScalar3dp<T>(v.w * v.w);
}

// return weigth norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(Vec3dp<T> const& v)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(v)));
}

// return squared weight norm of bivector
// ||weight(B)||^2 = lcmpl( dot(rcmpl(B),cmpl(B)) ) = rdot(B, B)
//                 = B.vx * B.vx + B.vy * B.vy + B.vz * B.vz
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(BiVec3dp<T> const& B)
{
    return PScalar3dp<T>(B.vx * B.vx + B.vy * B.vy + B.vz * B.vz);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(BiVec3dp<T> const& B)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(B)));
}

// return squared weight norm of trivector
// ||weight(t)||^2 = lcmpl( dot(rcmpl(t),cmpl(t)) ) = rdot(t, t)
//                 = t.x * t.x + t.y * t.y + t.z * t.z
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(TriVec3dp<T> const& t)
{
    return PScalar3dp<T>(t.x * t.x + t.y * t.y + t.z * t.z);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(TriVec3dp<T> const& t)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(t)));
}

// return squared weight norm of pseudoscalar
// ||weight(ps)||^2 = cmpl( dot(cmpl(ps),cmpl(ps)) ) = rdot(ps, ps) = ps * ps
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(T(ps) * T(ps));
}

// return weight norm of pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(ps)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(MVec3dp_E<T> const& M)
{
    return PScalar3dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr2(M)) +
                         weight_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(MVec3dp_E<T> const& M)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(MVec3dp_U<T> const& M)
{
    return PScalar3dp<T>(weight_nrm_sq(gr1(M)) + weight_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(MVec3dp_U<T> const& M)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(MVec3dp<T> const& M)
{
    return PScalar3dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr1(M)) +
                         weight_nrm_sq(gr2(M)) + weight_nrm_sq(gr3(M)) +
                         weight_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(MVec3dp<T> const& M)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(M)));
}


////////////////////////////////////////////////////////////////////////////////
// geometric norm: (perpendicular) distance to the origin (distance = c0/c1)
//                 i.e. the value c0 after unitization
//
// returns a dual number, for correct handling of objects at infinity (c1 == 0)
//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(Scalar3dp<T> s)
{
    return DualNum3dp<T>(bulk_nrm_sq(s), weight_nrm_sq(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(Scalar3dp<T> s)
{
    return DualNum3dp<T>(bulk_nrm(s), weight_nrm(s));
}

// provide the distance of the point from the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(Vec3dp<T> const& v)
{
    return DualNum3dp<T>(bulk_nrm_sq(v), weight_nrm_sq(v));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(Vec3dp<T> const& v)
{
    return DualNum3dp<T>(bulk_nrm(v), weight_nrm(v));
}

// provide the perpendicular distance of the line to the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(BiVec3dp<T> const& B)
{
    return DualNum3dp<T>(bulk_nrm_sq(B), weight_nrm_sq(B));
}

// return geometric norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(BiVec3dp<T> const& B)
{
    return DualNum3dp<T>(bulk_nrm(B), weight_nrm(B));
}

// provide the perpendicular distance of the plane to the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(TriVec3dp<T> const& t)
{
    return DualNum3dp<T>(bulk_nrm_sq(t), weight_nrm_sq(t));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(TriVec3dp<T> const& t)
{
    return DualNum3dp<T>(bulk_nrm(t), weight_nrm(t));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(PScalar3dp<T> ps)
{
    return DualNum3dp<T>(bulk_nrm_sq(ps), weight_nrm_sq(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(PScalar3dp<T> ps)
{
    return DualNum3dp<T>(bulk_nrm(ps), weight_nrm(ps));
}

// even grade multivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(MVec3dp_E<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(MVec3dp_E<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}

// odd grade multivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(MVec3dp_U<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(MVec3dp_U<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}

// multivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(MVec3dp<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(MVec3dp<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}


////////////////////////////////////////////////////////////////////////////////
// unitization operations
////////////////////////////////////////////////////////////////////////////////

// return a DualNum3dp unitized to weight_nrm == 1.0
// => if the dual number represents a homogeneous norm,
//    the scalar part represents the geometric norm the after unitization
template <typename T>
    requires(std::floating_point<T>)
inline DualNum3dp<T> unitize(DualNum3dp<T> const& D)
{
    T n = D.c1; // the pseudoscalar part is the weight_nrm part
    hd::ga::detail::check_unitization<T>(std::abs(n), "dual number (3dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * D;
}

// return a vector unitized to v.w == 1.0  (implies weight_nrm(v) = 1.0)
template <typename T>
    requires(std::floating_point<T>)
inline Vec3dp<T> unitize(Vec3dp<T> const& v)
{
    T n = v.w; // v.w == sign(v.w)*weight_nrm(v);
    hd::ga::detail::check_unitization<T>(std::abs(n), "vector (3dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec3dp<T>(v.x * inv, v.y * inv, v.z * inv, T(1.0));
}

// return a bivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline BiVec3dp<T> unitize(BiVec3dp<T> const& B)
{
    T n = weight_nrm(B);
    hd::ga::detail::check_unitization<T>(n, "bivector (3dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * B;
}

// return a trivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline TriVec3dp<T> unitize(TriVec3dp<T> const& t)
{
    T n = weight_nrm(t);
    hd::ga::detail::check_unitization<T>(n, "trivector (3dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * t;
}

// return an even grade multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline MVec3dp_E<T> unitize(MVec3dp_E<T> const& M)
{
    T n = weight_nrm(M);
    hd::ga::detail::check_unitization<T>(n, "even grade multivector (3dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}

// return an odd grade multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline MVec3dp_U<T> unitize(MVec3dp_U<T> const& M)
{
    T n = weight_nrm(M);
    hd::ga::detail::check_unitization<T>(n, "odd grade multivector (3dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}

// return a multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline MVec3dp<T> unitize(MVec3dp<T> const& M)
{
    T n = weight_nrm(M);
    hd::ga::detail::check_unitization<T>(n, "multivector (3dp)");
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}

template <typename T>
    requires std::floating_point<T>
inline Point3dp<T> unitize(Point3dp<T> const& p)
{
    hd::ga::detail::check_unitization<T>(std::abs(p.w), "Point3dp");
    T inv = T(1.0) / p.w;
    return Point3dp<T>(p.x * inv, p.y * inv, p.z * inv, 1.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline Line3d<T> unitize(Line3d<T> const& l)
{
    // unitization for a 3d bivector means a normalized direction vector
    // std::sqrt((l.vx)^2 + (l.vy)^2 + (l.vz)^2) = 1
    T wn = weight_nrm(l);
    hd::ga::detail::check_unitization<T>(wn, "Line3d");
    T inv = T(1.0) / wn;
    return Line3d<T>(l.vx * inv, l.vy * inv, l.vz * inv, l.mx * inv, l.my * inv,
                     l.mz * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline Plane3d<T> unitize(Plane3d<T> const& p)
{
    T wn = weight_nrm(p);
    hd::ga::detail::check_unitization<T>(wn, "Plane3d");
    T inv = T(1.0) / wn;
    return Plane3d<T>(p.x * inv, p.y * inv, p.z * inv, p.w * inv);
}


////////////////////////////////////////////////////////////////////////////////
// (right) bulk_dual (= right complement operation applied to the bulk)
//
// the duals are defined uniquely w.r.t. the outer product
//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> bulk_dual(Scalar3dp<T> s)
{
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> bulk_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> bulk_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> bulk_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), -t.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_dual([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> bulk_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(bulk_dual(gr4(M)), bulk_dual(gr2(M)), bulk_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> bulk_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(bulk_dual(gr3(M)), bulk_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> bulk_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(bulk_dual(gr4(M)), bulk_dual(gr3(M)), bulk_dual(gr2(M)),
                      bulk_dual(gr1(M)), bulk_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// (right) weight_dual (=right complement operation applied to the weight)
//
// the duals are defined uniquely w.r.t. the outer product
//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_dual([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> weight_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> weight_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> weight_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(-t.x, -t.y, -t.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> weight_dual(PScalar3dp<T> ps)
{
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_E<T> weight_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(weight_dual(gr4(M)), weight_dual(gr2(M)), weight_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp_U<T> weight_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(weight_dual(gr3(M)), weight_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3dp<T> weight_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(weight_dual(gr4(M)), weight_dual(gr3(M)), weight_dual(gr2(M)),
                      weight_dual(gr1(M)), weight_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// attitude operations: att = rwdg( u, rcmpl(e4_3dp) ) = rwdg(u, horizon_3dp)
//
// (the attitude is the intersection of the object with the horizon)
// the result of att(object_with_grade_k) is an object with grade k-1
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> att(Vec3dp<T> const& v)
{
    return Scalar3dp<T>(v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> att(Point3dp<T> const& p)
{
    return Scalar3dp<T>(p.w);
}

// return the attitude (i.e. the direction vector) of the line
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> att(BiVec3dp<T> const& B)
{
    return Vec3dp<T>(B.vx, B.vy, B.vz, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> att(Line3d<T> const& l)
{
    return Vec3dp<T>(l.vx, l.vy, l.vz, T(0.0));
}

// return the attitude (i.e. the normal vector) of the plane
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> att(TriVec3dp<T> const& t)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), t.x, t.y, t.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> att(Plane3d<T> const& l)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), l.x, l.y, l.z);
}


template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> att(PScalar3dp<T> ps)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), ps);
}

} // namespace hd::ga::pga