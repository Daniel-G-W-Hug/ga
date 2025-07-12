#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "detail/ga_foundation.hpp"     // ga library headers and infrastructure
#include "detail/type_t/ga_mvec2dp.hpp" // 2dp multivector types
#include "detail/type_t/ga_type2dp.hpp" // Point2dp, Vector2d, Point2d, Line2d

#include "detail/ga_error_handling.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga2dp basic operations:
//
// - gr_inv()                     -> grade inversion
// - rev()                        -> reversion
// - rrev()                       -> regressive reversion
// - conj()                       -> conjugation
// - cmpl()                       -> complement
// - dual()                       -> dual
//
// - bulk(), weight()             -> return bulk and weight parts of objects
// - bulk_nrm_sq, bulk_nrm        -> return bulk norm
// - weight_nrm_sq, weight_nrm    -> return weight norm
// - geom_nrm_sq, geom_nrm        -> return geometric norm
// - unitization                  -> return unitized object (weight_nrm scaled to 1.0)
//
// - bulk_dual                    -> return bulk dual
// - weight_dual                  -> return weight dual
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> gr_inv(Scalar2dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> gr_inv(Vec2dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> gr_inv(BiVec2dp<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> gr_inv(PScalar2dp<T> ps)
{
    // grade 3: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> gr_inv(MVec2dp_E<T> const& M)
{
    // grade 0 and 2: no sign change
    return M;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> gr_inv(MVec2dp_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> gr_inv(MVec2dp<T> const& M)
{
    // grade 0 and 2: no sign change
    // grade 1 and 3: sign reversal
    return MVec2dp<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> rev(Scalar2dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> rev(Vec2dp<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> rev(BiVec2dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> rev(PScalar2dp<T> ps)
{
    // grade 3: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> rev(MVec2dp_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec2dp_E<T>(rev(gr0(M)), rev(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> rev(MVec2dp_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec2dp_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> rev(MVec2dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    return MVec2dp<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// regressive reversion operation: reverse applied to the complement + backtrafo
// rrev(A_r) = cmpl((-1)^(r*(r-1)/2) cmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=3, r = 0, 1, 2, 3, ...: - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> rrev(Scalar2dp<T> s)
{
    // grade 0: sign reversal
    return -s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> rrev(Vec2dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> rrev(BiVec2dp<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> rrev(PScalar2dp<T> ps)
{
    // grade 3: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> rrev(MVec2dp_E<T> const& M)
{
    // grade 0: sign reversal
    // grade 2: no sign change
    return MVec2dp_E<T>(rrev(gr0(M)), rrev(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> rrev(MVec2dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dp_U<T>(rrev(gr1(M)), rrev(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> rrev(MVec2dp<T> const& M)
{
    // grade 0: sign reversal
    // grade 1: sign reversal
    // grade 2: no sign change
    // grade 3: no sign change
    return MVec2dp<T>(rrev(gr0(M)), rrev(gr1(M)), rrev(gr2(M)), rrev(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> conj(Scalar2dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> conj(Vec2dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> conj(BiVec2dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> conj(PScalar2dp<T> ps)
{
    // grade 3: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> conj(MVec2dp_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec2dp_E<T>(conj(gr0(M)), conj(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> conj(MVec2dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dp_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> conj(MVec2dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    return MVec2dp<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// 2dp complement operations
// (the concept of complement is defined w.r.t. the outer product)
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^2 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors which are in the k-blade u with
// the basis vectors which are NOT contained in the k-blade u
// and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  lcmpl(u) ^ u  = I_2dp = e1^e2^e3
// right complement: u ^ rcmpl(u)  = I_2dp = e1^e2^e3
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs

// complement operation with e1^e2^e3 as the pseudoscalar

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> cmpl(Scalar2dp<T> s)
{
    // u ^ cmpl(u) = e3^e2^e1
    // u = 1:
    //     1 ^ cmpl(u) = e2^e2^e1 => cmpl(u) = s e3^e2^e1
    return PScalar2dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> cmpl(Vec2dp<T> const& v)
{
    // u ^ compl(u) = e3^e2^e1
    // u = v.x e1 + v.y e2 + v.z e3:
    //     u ^ cmpl(u) = e3^e2^e1 =>
    //     u = e1 => cmpl(u) = -e23
    //     u = e2 => cmpl(u) = -e31
    //     u = e3 => cmpl(u) = -e12
    return BiVec2dp<T>(-v.x, -v.y, -v.z);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> cmpl(BiVec2dp<T> const& B)
{
    // u ^ compl(u) = e3^e2^e1
    // u = b.x e23 + b.y e31 + b.z e12:
    //     u ^ cmpl(u) = e3^e2^e1 =>
    //     u = e23 => cmpl(u) = -e1
    //     u = e31 => cmpl(u) = -e2
    //     u = e12 => cmpl(u) = -e3
    return Vec2dp<T>(-B.x, -B.y, -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> cmpl(PScalar2dp<T> ps)
{
    // u ^ compl(u) = e3^e2^e1
    // u = e3^e2^e1:
    //     e3^e2^e1 ^ cmpl(u) = e3^e2^e1 => cmpl(u) = ps 1
    return Scalar2dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> cmpl(MVec2dp_E<T> const& M)
{
    // use the component complements directly
    return MVec2dp_U<T>(cmpl(gr2(M)), cmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> cmpl(MVec2dp_U<T> const& M)
{
    // use the component complements directly
    return MVec2dp_E<T>(cmpl(gr3(M)), cmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> cmpl(MVec2dp<T> const& M)
{
    // use the component complements directly
    return MVec2dp<T>(cmpl(gr3(M)), cmpl(gr2(M)), cmpl(gr1(M)), cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk: u_bulk = G u (with G as the metric)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk(Scalar2dp<T> s)
{
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> bulk(Vec2dp<T> const& v)
{
    return Vec2dp<T>(v.x, v.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> bulk(BiVec2dp<T> const& B)
{
    return BiVec2dp<T>(T(0.0), T(0.0), B.z);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> bulk([[maybe_unused]] PScalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> bulk(MVec2dp_E<T> const& M)
{
    return MVec2dp_E<T>(bulk(gr0(M)), bulk(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> bulk(MVec2dp_U<T> const& M)
{
    return MVec2dp_U<T>(bulk(gr1(M)), bulk(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> bulk(MVec2dp<T> const& M)
{
    return MVec2dp<T>(bulk(gr0(M)), bulk(gr1(M)), bulk(gr2(M)), bulk(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight: u_weight = lcmpl( G rcmpl(u) ) = rG u
//         (with G as the metric and rG as the anti-metric as given by Lengyel)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> weight([[maybe_unused]] Scalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> weight(Vec2dp<T> const& v)
{
    return Vec2dp<T>(T(0.0), T(0.0), v.z);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> weight(BiVec2dp<T> const& B)
{
    return BiVec2dp<T>(B.x, B.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight(PScalar2dp<T> ps)
{
    return PScalar2dp<T>(ps);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> weight(MVec2dp_E<T> const& M)
{
    return MVec2dp_E<T>(weight(gr0(M)), weight(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> weight(MVec2dp_U<T> const& M)
{
    return MVec2dp_U<T>(weight(gr1(M)), weight(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> weight(MVec2dp<T> const& M)
{
    return MVec2dp<T>(weight(gr0(M)), weight(gr1(M)), weight(gr2(M)), weight(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk norm
////////////////////////////////////////////////////////////////////////////////

// return squared bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm_sq(Scalar2dp<T> s)
{
    // ||bulk(s)||^2 = dot(s,s) = s^(T) G s = gr0(s*rev(s)) = s*s
    // using rev(s) = (-1)^[k(k-1)/2] s for a k-blade: 0-blade => rev(s) = s
    return Scalar2dp<T>(T(s) * T(s));
}

// return bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm(Scalar2dp<T> s)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(s)));
}

// return squared bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm_sq(Vec2dp<T> const& v)
{
    // ||bulk(v)||^2 = dot(v,v) = v^(T) G v = gr0(v*rev(v)) = v.x * v.x + v.y * v.y
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    return Scalar2dp<T>(v.x * v.x + v.y * v.y);
}

// return bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm(Vec2dp<T> const& v)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(v)));
}

// return squared bulk magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm_sq(BiVec2dp<T> const& B)
{
    // ||bulk(B)||^2 = dot(B,B) = B^(T) G B = gr0(B*rev(B)) = B.z * B.z
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    return Scalar2dp<T>(B.z * B.z);
}

// return magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm(BiVec2dp<T> const& B)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(B)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm_sq([[maybe_unused]] PScalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm([[maybe_unused]] PScalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm_sq(MVec2dp_E<T> const& M)
{
    return Scalar2dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm(MVec2dp_E<T> const& M)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm_sq(MVec2dp_U<T> const& M)
{
    return Scalar2dp<T>(bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm(MVec2dp_U<T> const& M)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm_sq(MVec2dp<T> const& M)
{
    return Scalar2dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr2(M)) +
                        bulk_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_nrm(MVec2dp<T> const& M)
{
    return Scalar2dp<T>(std::sqrt(bulk_nrm_sq(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight norm
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> weight_nrm_sq([[maybe_unused]] Scalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> weight_nrm([[maybe_unused]] Scalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

// return squared weight norm of vector
// ||weight(v)||^2 = cmpl( dot(cmpl(v), cmpl(v)) ) = rdot(v,v) = v.z * v.z
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm_sq(Vec2dp<T> const& v)
{
    return PScalar2dp<T>(v.z * v.z);
}

// return weigth norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm(Vec2dp<T> const& v)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(v)));
}

// return squared weight norm of bivector
// ||weight(B)||^2 = cmpl( dot(cmpl(B),cmpl(B)) ) = rdot(B, B) = B.x * B.x + B.y * B.y
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm_sq(BiVec2dp<T> const& B)
{
    return PScalar2dp<T>(B.x * B.x + B.y * B.y);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm(BiVec2dp<T> const& B)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(B)));
}

// return squared weight norm of pseudoscalar
// ||weight(ps)||^2 = cmpl( dot(cmpl(ps),cmpl(ps)) ) = rdot(ps, ps) = ps * ps
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm_sq(PScalar2dp<T> ps)
{
    return PScalar2dp<T>(T(ps) * T(ps));
}

// return weight norm of pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm(PScalar2dp<T> ps)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(ps)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm_sq(MVec2dp_E<T> const& M)
{
    return PScalar2dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm(MVec2dp_E<T> const& M)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm_sq(MVec2dp_U<T> const& M)
{
    return PScalar2dp<T>(weight_nrm_sq(gr1(M)) + weight_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm(MVec2dp_U<T> const& M)
{
    return PScalar2dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm_sq(MVec2dp<T> const& M)
{
    return PScalar2dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr1(M)) +
                         weight_nrm_sq(gr2(M)) + weight_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> weight_nrm(MVec2dp<T> const& M)
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
constexpr DualNum2dp<T> geom_nrm_sq(Scalar2dp<T> s)
{
    return DualNum2dp<T>(bulk_nrm_sq(s), weight_nrm_sq(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm(Scalar2dp<T> s)
{
    return DualNum2dp<T>(bulk_nrm(s), weight_nrm(s));
}

// provide the distance of the point from the origin
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm_sq(Vec2dp<T> const& v)
{
    return DualNum2dp<T>(bulk_nrm_sq(v), weight_nrm_sq(v));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm(Vec2dp<T> const& v)
{
    return DualNum2dp<T>(bulk_nrm(v), weight_nrm(v));
}

// provide the perpendicular distance of the line to the origin
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm_sq(BiVec2dp<T> const& B)
{
    return DualNum2dp<T>(bulk_nrm_sq(B), weight_nrm_sq(B));
}

// return geometric norm of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm(BiVec2dp<T> const& B)
{
    return DualNum2dp<T>(bulk_nrm(B), weight_nrm(B));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm_sq(PScalar2dp<T> ps)
{
    return DualNum2dp<T>(bulk_nrm_sq(ps), weight_nrm_sq(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm(PScalar2dp<T> ps)
{
    return DualNum2dp<T>(bulk_nrm(ps), weight_nrm(ps));
}

// even grade multivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm_sq(MVec2dp_E<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm(MVec2dp_E<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm(M), weight_nrm(M));
}

// odd grade multivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm_sq(MVec2dp_U<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm(MVec2dp_U<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm(M), weight_nrm(M));
}

// multivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm_sq(MVec2dp<T> const& M)
{
    return DualNum2dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum2dp<T> geom_nrm(MVec2dp<T> const& M)
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
// bulk_dual = cmpl(bulk(A))
//
// -> complement operation applied to the bulk
//
////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Lengyel, "PGA illuminated")
// is defined w.r.t. the outer product, because the complement is defined w.r.t.
// the outer product.
//
// The only difference between the dual and the complement is the fact that
// for the dual the argument is multiplied with the extended metric before
// taking the complement. However, in ega the metric is the identity matrix.
// Thus the dual and the complement are identical.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> bulk_dual(Scalar2dp<T> s)
{
    return PScalar2dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> bulk_dual(Vec2dp<T> const& v)
{
    return BiVec2dp<T>(-v.x, -v.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> bulk_dual(BiVec2dp<T> const& B)
{
    return Vec2dp<T>(T(0.0), T(0.0), -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> bulk_dual([[maybe_unused]] PScalar2dp<T>)
{
    return Scalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> bulk_dual(MVec2dp_E<T> const& M)
{
    return MVec2dp_U<T>(bulk_dual(gr2(M)), bulk_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> bulk_dual(MVec2dp_U<T> const& M)
{
    return MVec2dp_E<T>(bulk_dual(gr3(M)), bulk_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> bulk_dual(MVec2dp<T> const& M)
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
constexpr PScalar2dp<T> weight_dual([[maybe_unused]] Scalar2dp<T>)
{
    return PScalar2dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> weight_dual(Vec2dp<T> const& v)
{
    return BiVec2dp<T>(T(0.0), T(0.0), -v.z);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> weight_dual(BiVec2dp<T> const& B)
{
    return Vec2dp<T>(-B.x, -B.y, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> weight_dual(PScalar2dp<T> ps)
{
    return Scalar2dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_U<T> weight_dual(MVec2dp_E<T> const& M)
{
    return MVec2dp_U<T>(weight_dual(gr2(M)), weight_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp_E<T> weight_dual(MVec2dp_U<T> const& M)
{
    return MVec2dp_E<T>(weight_dual(gr3(M)), weight_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec2dp<T> weight_dual(MVec2dp<T> const& M)
{
    return MVec2dp<T>(weight_dual(gr3(M)), weight_dual(gr2(M)), weight_dual(gr1(M)),
                      weight_dual(gr0(M)));
}

} // namespace hd::ga::pga