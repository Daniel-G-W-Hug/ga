#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "type_t/ga_type_2d.hpp"
#include "type_t/ga_type_2dp.hpp"
#include "type_t/ga_vec2_t.hpp"
#include "type_t/ga_vec3_t.hpp"


namespace hd::ga::pga {

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

// uneven grade multivector
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
inline constexpr DualNum2dp<T> unitize(DualNum2dp<T> const& D)
{
    T n = D.c1; // the pseudoscalar part is the weight_nrm part
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("DualNum2dp weight_norm too small for unitization " +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * D;
}

// return a vector unitized to v.z == 1.0  (implies weight_nrm(v) = 1.0)
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2dp<T> unitize(Vec2dp<T> const& v)
{
    T n = v.z;
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector weight_nrm too small for unitization " +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec2dp<T>(v.x * inv, v.y * inv, T(1.0));
}

// return a bivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> unitize(BiVec2dp<T> const& B)
{
    T n = weight_nrm(B);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector weight_norm too small for unitization " +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * B;
}

// return an even grade multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_E<T> unitize(MVec2dp_E<T> const& M)
{
    T n = weight_nrm(M);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error(
            "Even grade mulitivector weight_norm too small for unitization " +
            std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}

// return an uneven grade multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp_U<T> unitize(MVec2dp_U<T> const& M)
{
    T n = weight_nrm(M);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error(
            "Uneven grade mulitivector weight_norm too small for unitization " +
            std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
}

// return a multivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2dp<T> unitize(MVec2dp<T> const& M)
{
    T n = weight_nrm(M);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("Mulitivector weight_norm too small for unitization " +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * M;
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
// convenient type aliases
////////////////////////////////////////////////////////////////////////////////

// Vector2d: 2d vector of projective algebra storing only components x, y explicitly.
// The z component is assumed to be z = 0
template <typename T>
    requires(std::floating_point<T>)
struct Vector2d : public Vec2d<T> {

    using Vec2d<T>::Vec2d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2d<T>::x;
    using Vec2d<T>::y;

    Vector2d(Vec2d<T> const& v) : Vec2d<T>(v) {};
};

// Point2d: 2d point of projective algebra storing only components x, y explicitly.
// The z component is assumed to by z=1
template <typename T>
    requires(std::floating_point<T>)
struct Point2d : public Vec2d<T> {

    using Vec2d<T>::Vec2d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2d<T>::x;
    using Vec2d<T>::y;

    Point2d(Vec2d<T> const& v) : Vec2d<T>(v) {};
    Point2d(T x, T y) : Vec2d<T>(x, y) {};
};

// Point2dp: 2d point of projective algebra storing all three components x, y, z
//           explicitly
//
// a Point2dp is a Vec2dp, thus all operations defined for Vec2dp
// work directly for Point2dp - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Point2dp : public Vec2dp<T> {

    using Vec2dp<T>::Vec2dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2dp<T>::x;
    using Vec2dp<T>::y;
    using Vec2dp<T>::z;

    Point2dp(Vec2dp<T> const& v) : Vec2dp<T>(v) {};
    Point2dp(Point2d<T> const& p) : Vec2dp<T>(p.x, p.y, T(1.0)) {};
    Point2dp(Vec2d<T> const& v) : Vec2dp<T>(v.x, v.y, T(1.0)) {};

    Point2dp& unitize()
    {
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (std::abs(z) < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("z-component too small for unitization " +
                                     std::to_string(z) + "\n");
        }
#endif
        x /= z;
        y /= z;
        z = T(1.0);
        return (*this);
    }
};

template <typename T>
    requires std::floating_point<T>
inline constexpr Point2dp<T> unitize(Point2dp<T> const& p)
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(p.z) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("z-component too small for unitization " +
                                 std::to_string(p.z) + "\n");
    }
#endif
    T inv = T(1.0) / p.z;
    return Point2dp<T>(p.x * inv, p.y * inv, T(1.0));
}


// Line2d: 2d line of projective algebra storing all three components
//          explicitly as components x, y, z of a BiVec2dp
//
// a Line2d is a BiVec2dp, thus all operations defined for BiVec2dp
// work directly for Line2d - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Line2d : public BiVec2dp<T> {

    using BiVec2dp<T>::BiVec2dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using BiVec2dp<T>::x;
    using BiVec2dp<T>::y;
    using BiVec2dp<T>::z;

    Line2d() = default;
    Line2d(T x, T y, T z) : BiVec2dp<T>(x, y, z) {};
    Line2d(BiVec2dp<T> const& b) : BiVec2dp<T>(b) {};
    Line2d(Point2d<T> const& p, Point2d<T> const& q) :
        // Line2d = wdg(p,q), but wdg() cannot be used here to avoid circular dependency
        BiVec2dp<T>(p.y - q.y, q.x - p.x, p.x * q.y - p.y * q.x) {};
    Line2d(Point2d<T> const& p, Vec2d<T> const& v) :
        // Line2d constructed from a point and a direction vector
        BiVec2dp<T>(-v.y, v.x, p.x * v.y - p.y * v.x) {};

    Line2d& unitize()
    {
        // unitization for a 2d bivector means std::sqrt(x^2 + y^2) = 1
        // = weight_nrm of bivector
        T wn = weight_nrm(*this);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (wn < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("bivector weight norm too small for unitization " +
                                     std::to_string(wn) + "\n");
        }
#endif
        T inv = T(1.0) / wn;
        x *= inv;
        y *= inv;
        z *= inv;
        return (*this);
    }
};

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Line2d<T> unitize(Line2d<T> const& l)
{
    // unitization for a 2d bivector means std::sqrt(x^2 + y^2) = 1
    // i.e. unitization of the direction vector of the line
    T wn = weight_nrm(l);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (wn < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector weight norm too small for unitization " +
                                 std::to_string(wn) + "\n");
    }
#endif
    T inv = T(1.0) / wn;
    return Line2d<T>(l.x * inv, l.y * inv, l.z * inv);
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