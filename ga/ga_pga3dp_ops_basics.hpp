#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "detail/ga_foundation.hpp"     // ga library headers and infrastructure
#include "detail/type_t/ga_mvec3dp.hpp" // 3dp multivector types
#include "detail/type_t/ga_type3dp.hpp" // Point3dp, Vector3d, Point3d, Line3d, Plane3d

#include "detail/ga_error_handling.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga3dp basic operations:
//
// - gr_inv()                     -> grade inversion
// - rev()                        -> reversion
// - rrev()                       -> regressive reversion
// - conj()                       -> conjugation
// - lcmpl(), rcmpl()             -> left and right complement
//
// - bulk(), weight()             -> return bulk and weight parts of objects
// - bulk_nrm_sq, bulk_nrm        -> return bulk norm
// - weight_nrm_sq, weight_nrm    -> return weight norm
// - geom_nrm_sq, geom_nrm        -> return geometric norm
// - unitization                  -> return unitized object (weight_nrm scaled to 1.0)
//
// - right_bulk_dual              -> return right bulk dual
// - right_weight_dual            -> return right weight dual
// - left_bulk_dual               -> return right bulk dual
// - left_weight_dual             -> return right weight dual
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> gr_inv(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> gr_inv(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> gr_inv(BiVec3dp<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> gr_inv(TriVec3dp<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> gr_inv(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> gr_inv(MVec3dp_E<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    return M;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> gr_inv(MVec3dp_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> gr_inv(MVec3dp<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    // grade 1 and 3: sign reversal
    return MVec3dp<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)),
                      gr_inv(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> rev(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> rev(Vec3dp<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> rev(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> rev(TriVec3dp<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> rev(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> rev(MVec3dp_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec3dp_E<T>(rev(gr0(M)), rev(gr2(M)), rev(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> rev(MVec3dp_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec3dp_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> rev(MVec3dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    // grade 4: no sign change
    return MVec3dp<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)), rev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// regressive reversion operation: reverse of complement + backtrafo
// rrev(A_r) = lcmpl((-1)^(r*(r-1)/2) rcmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=4, r = 0, 1, 2, 3, ...: + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> rrev(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> rrev(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> rrev(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> rrev(TriVec3dp<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> rrev(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> rrev(MVec3dp_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec3dp_E<T>(rrev(gr0(M)), rrev(gr2(M)), rrev(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> rrev(MVec3dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3dp_U<T>(rrev(gr1(M)), rrev(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> rrev(MVec3dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec3dp<T>(rrev(gr0(M)), rrev(gr1(M)), rrev(gr2(M)), rrev(gr3(M)),
                      rrev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> conj(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> conj(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> conj(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> conj(TriVec3dp<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> conj(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> conj(MVec3dp_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec3dp_E<T>(conj(gr0(M)), conj(gr2(M)), conj(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> conj(MVec3dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3dp_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> conj(MVec3dp<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec3dp<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)),
                      conj(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// 3dp complement operations
// (the concept of complement is defined w.r.t. the outer product)
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^2 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors which are in the k-blade u with
// the basis vectors which are NOT contained in the k-blade u
// and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  lcmpl(u) ^ u  = I_3dp = e1^e2^e3^e4
// right complement: u ^ rcmpl(u)  = I_3dp = e1^e2^e3^e4
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs
//
// complement operation with e1^e2^e3^e4 as the pseudoscalar

// right complements

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> rcmpl(Scalar3dp<T> s)
{
    // u ^ rcmpl(u) = e1^e2^e3^e4
    // u = 1:
    //     1 ^ rcmpl(u) = e1^e2^e3^e4 => rcmpl(u) = s e1^e2^e3^e4
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> rcmpl(Vec3dp<T> const& v)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = v.x e1 + v.y e2 + v.z e3 + v.w e4:
    //     u ^ rcmpl(u) = e1^e2^e3^e4 =>
    //     u = e1 => rcmpl(u) = e423
    //     u = e2 => rcmpl(u) = e431
    //     u = e3 => rcmpl(u) = e412
    //     u = e4 => rcmpl(u) = e321
    return TriVec3dp<T>(v.x, v.y, v.z, v.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> rcmpl(BiVec3dp<T> const& B)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     u ^ rcmpl(u) = e1^e2^e3^e4 =>
    //     u = e41 => rcmpl(u) = -e23
    //     u = e42 => rcmpl(u) = -e31
    //     u = e43 => rcmpl(u) = -e12
    //     u = e23 => rcmpl(u) = -e41
    //     u = e31 => rcmpl(u) = -e42
    //     u = e12 => rcmpl(u) = -e43
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> rcmpl(TriVec3dp<T> const& t)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     u ^ rcmpl(u) = e1^e2^e3^e4 =>
    //     u = e423 => rcmpl(u) = -e1
    //     u = e431 => rcmpl(u) = -e2
    //     u = e412 => rcmpl(u) = -e3
    //     u = e321 => rcmpl(u) = -e4
    return Vec3dp<T>(-t.x, -t.y, -t.z, -t.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> rcmpl(PScalar3dp<T> ps)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     e1^e2^e3^e4 ^ rcmpl(u) = e1^e2^e3^e4 => rcmpl(u) = ps 1
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> rcmpl(MVec3dp_U<T> const& M)
{
    // use the component complements directly
    return MVec3dp_U<T>(rcmpl(gr3(M)), rcmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> rcmpl(MVec3dp_E<T> const& M)
{
    // use the component complements directly
    return MVec3dp_E<T>(rcmpl(gr4(M)), rcmpl(gr2(M)), rcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> rcmpl(MVec3dp<T> const& M)
{
    // use the component complements directly
    return MVec3dp<T>(rcmpl(gr4(M)), rcmpl(gr3(M)), rcmpl(gr2(M)), rcmpl(gr1(M)),
                      rcmpl(gr0(M)));
}


// left complements

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> lcmpl(Scalar3dp<T> s)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = 1:
    //     lcmpl(u) ^ 1 = e1^e2^e3^e4 => lcmpl(u) = s e1^e2^e3^e4
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> lcmpl(Vec3dp<T> const& v)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = v.x e1 + v.y e2 + v.z e3 + v.w e4:
    //     lcmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e1 => lcmpl(u) = -e423
    //     u = e2 => lcmpl(u) = -e431
    //     u = e3 => lcmpl(u) = -e412
    //     u = e4 => lcmpl(u) = -e321
    return TriVec3dp<T>(-v.x, -v.y, -v.z, -v.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> lcmpl(BiVec3dp<T> const& B)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     lcmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e41 => lcmpl(u) = -e23
    //     u = e42 => lcmpl(u) = -e31
    //     u = e43 => lcmpl(u) = -e12
    //     u = e23 => lcmpl(u) = -e41
    //     u = e31 => lcmpl(u) = -e42
    //     u = e12 => lcmpl(u) = -e43
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> lcmpl(TriVec3dp<T> const& t)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     lcmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e423 => lcmpl(u) = e1
    //     u = e431 => lcmpl(u) = e2
    //     u = e412 => lcmpl(u) = e3
    //     u = e321 => lcmpl(u) = e4
    return Vec3dp<T>(t.x, t.y, t.z, t.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> lcmpl(PScalar3dp<T> ps)
{
    // lcmpl(u) ^ u = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     lcmpl(u) ^ e1^e2^e3^e4 = e1^e2^e3^e4 => lcmpl(u) = ps 1
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> lcmpl(MVec3dp_U<T> const& M)
{
    // use the component complements directly
    return MVec3dp_U<T>(lcmpl(gr3(M)), lcmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> lcmpl(MVec3dp_E<T> const& M)
{
    // use the component complements directly
    return MVec3dp_E<T>(lcmpl(gr4(M)), lcmpl(gr2(M)), lcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> lcmpl(MVec3dp<T> const& M)
{
    // use the component complements directly
    return MVec3dp<T>(lcmpl(gr4(M)), lcmpl(gr3(M)), lcmpl(gr2(M)), lcmpl(gr1(M)),
                      lcmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk: u_bulk = G u (with G as the metric)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk(Scalar3dp<T> s)
{
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> bulk(Vec3dp<T> const& v)
{
    return Vec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> bulk(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), B.mx, B.my, B.mz);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> bulk(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), t.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> bulk([[maybe_unused]] PScalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> bulk(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(bulk(gr0(M)), bulk(gr2(M)), bulk(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> bulk(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(bulk(gr1(M)), bulk(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> bulk(MVec3dp<T> const& M)
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
constexpr Scalar3dp<T> weight([[maybe_unused]] Scalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> weight(Vec3dp<T> const& v)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> weight(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(B.vx, B.vy, B.vz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> weight(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(t.x, t.y, t.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(ps);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> weight(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(weight(gr0(M)), weight(gr2(M)), weight(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> weight(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(weight(gr1(M)), weight(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> weight(MVec3dp<T> const& M)
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
constexpr Scalar3dp<T> bulk_nrm_sq(Scalar3dp<T> s)
{
    // ||bulk(s)||^2 = dot(s,s) = s^(T) G s = gr0(s*rev(s)) = s*s
    // using rev(s) = (-1)^[k(k-1)/2] s for a k-blade: 0-blade => rev(s) = s
    return Scalar3dp<T>(T(s) * T(s));
}

// return bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm(Scalar3dp<T> s)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(s)));
}

// return squared bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm_sq(Vec3dp<T> const& v)
{
    // ||bulk(v)||^2 = dot(v,v) = v^(T) G v = gr0(v*rev(v))
    //               = v.x * v.x + v.y * v.y + v.z * v.z
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    return Scalar3dp<T>(v.x * v.x + v.y * v.y + v.z * v.z);
}

// return bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm(Vec3dp<T> const& v)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(v)));
}

// return squared bulk magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm_sq(BiVec3dp<T> const& B)
{
    // ||bulk(B)||^2 = dot(B,B) = B^(T) G B = gr0(B*rev(B))
    //               = B.mx * B.mx + B.my * B.my + B.mz * B.mz
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    return Scalar3dp<T>(B.mx * B.mx + B.my * B.my + B.mz * B.mz);
}

// return magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm(BiVec3dp<T> const& B)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(B)));
}

// return squared bulk norm of trivector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm_sq(TriVec3dp<T> const& t)
{
    // ||bulk(t)||^2 = dot(t,t) = t^(T) G t = gr0(t*rev(t)) = t.w * t.w
    // using rev(t) = (-1)^[k(k-1)/2] t for a k-blade: 3-blade => rev(t) = -t
    return Scalar3dp<T>(t.w * t.w);
}

// return bulk norm of trivector
template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm(TriVec3dp<T> const& t)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(t)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm_sq([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm_sq(MVec3dp_E<T> const& M)
{
    return Scalar3dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr2(M)) + bulk_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm(MVec3dp_E<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm_sq(MVec3dp_U<T> const& M)
{
    return Scalar3dp<T>(bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm(MVec3dp_U<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm_sq(MVec3dp<T> const& M)
{
    return Scalar3dp<T>(bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr2(M)) +
                        bulk_nrm_sq(gr3(M)) + bulk_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> bulk_nrm(MVec3dp<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight norm
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

// return squared weight norm of vector
// ||weight(v)||^2 = lcmpl( dot(rcmpl(v), rcmpl(v)) ) = rdot(v,v) = v.w * v.w
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq(Vec3dp<T> const& v)
{
    return PScalar3dp<T>(v.w * v.w);
}

// return weigth norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm(Vec3dp<T> const& v)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(v)));
}

// return squared weight norm of bivector
// ||weight(B)||^2 = lcmpl( dot(rcmpl(B),cmpl(B)) ) = rdot(B, B)
//                 = B.vx * B.vx + B.vy * B.vy + B.vz * B.vz
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq(BiVec3dp<T> const& B)
{
    return PScalar3dp<T>(B.vx * B.vx + B.vy * B.vy + B.vz * B.vz);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm(BiVec3dp<T> const& B)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(B)));
}

// return squared weight norm of trivector
// ||weight(t)||^2 = lcmpl( dot(rcmpl(t),cmpl(t)) ) = rdot(t, t)
//                 = t.x * t.x + t.y * t.y + t.z * t.z
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq(TriVec3dp<T> const& t)
{
    return PScalar3dp<T>(t.x * t.x + t.y * t.y + t.z * t.z);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm(TriVec3dp<T> const& t)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(t)));
}

// return squared weight norm of pseudoscalar
// ||weight(ps)||^2 = cmpl( dot(cmpl(ps),cmpl(ps)) ) = rdot(ps, ps) = ps * ps
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(T(ps) * T(ps));
}

// return weight norm of pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(ps)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq(MVec3dp_E<T> const& M)
{
    return PScalar3dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr2(M)) +
                         weight_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm(MVec3dp_E<T> const& M)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq(MVec3dp_U<T> const& M)
{
    return PScalar3dp<T>(weight_nrm_sq(gr1(M)) + weight_nrm_sq(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm(MVec3dp_U<T> const& M)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm_sq(MVec3dp<T> const& M)
{
    return PScalar3dp<T>(weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr1(M)) +
                         weight_nrm_sq(gr2(M)) + weight_nrm_sq(gr3(M)) +
                         weight_nrm_sq(gr4(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> weight_nrm(MVec3dp<T> const& M)
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
constexpr DualNum3dp<T> geom_nrm_sq(Scalar3dp<T> s)
{
    return DualNum3dp<T>(bulk_nrm_sq(s), weight_nrm_sq(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(Scalar3dp<T> s)
{
    return DualNum3dp<T>(bulk_nrm(s), weight_nrm(s));
}

// provide the distance of the point from the origin
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm_sq(Vec3dp<T> const& v)
{
    return DualNum3dp<T>(bulk_nrm_sq(v), weight_nrm_sq(v));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(Vec3dp<T> const& v)
{
    return DualNum3dp<T>(bulk_nrm(v), weight_nrm(v));
}

// provide the perpendicular distance of the line to the origin
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm_sq(BiVec3dp<T> const& B)
{
    return DualNum3dp<T>(bulk_nrm_sq(B), weight_nrm_sq(B));
}

// return geometric norm of bivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(BiVec3dp<T> const& B)
{
    return DualNum3dp<T>(bulk_nrm(B), weight_nrm(B));
}

// provide the perpendicular distance of the plane to the origin
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm_sq(TriVec3dp<T> const& t)
{
    return DualNum3dp<T>(bulk_nrm_sq(t), weight_nrm_sq(t));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(TriVec3dp<T> const& t)
{
    return DualNum3dp<T>(bulk_nrm(t), weight_nrm(t));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm_sq(PScalar3dp<T> ps)
{
    return DualNum3dp<T>(bulk_nrm_sq(ps), weight_nrm_sq(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(PScalar3dp<T> ps)
{
    return DualNum3dp<T>(bulk_nrm(ps), weight_nrm(ps));
}

// even grade multivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm_sq(MVec3dp_E<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(MVec3dp_E<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}

// odd grade multivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm_sq(MVec3dp_U<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(MVec3dp_U<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}

// multivector
template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm_sq(MVec3dp<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm_sq(M), weight_nrm_sq(M));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr DualNum3dp<T> geom_nrm(MVec3dp<T> const& M)
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
// right_bulk_dual = rcmpl(bulk(A))
//
// -> right complement operation applied to the bulk
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
// Due to the degenerate metric in pga some elements are zero (not the identity matrix!).
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> right_bulk_dual(Scalar3dp<T> s)
{
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> right_bulk_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> right_bulk_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> right_bulk_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), -t.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> right_bulk_dual([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> right_bulk_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(right_bulk_dual(gr4(M)), right_bulk_dual(gr2(M)),
                        right_bulk_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> right_bulk_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(right_bulk_dual(gr3(M)), right_bulk_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> right_bulk_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(right_bulk_dual(gr4(M)), right_bulk_dual(gr3(M)),
                      right_bulk_dual(gr2(M)), right_bulk_dual(gr1(M)),
                      right_bulk_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// right_weight_dual(A) = = rcmpl(weight(A))
//
// -> right complement operation applied to the weight
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
constexpr PScalar3dp<T> right_weight_dual([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> right_weight_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> right_weight_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> right_weight_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(-t.x, -t.y, -t.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> right_weight_dual(PScalar3dp<T> ps)
{
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> right_weight_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(right_weight_dual(gr4(M)), right_weight_dual(gr2(M)),
                        right_weight_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> right_weight_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(right_weight_dual(gr3(M)), right_weight_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> right_weight_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(right_weight_dual(gr4(M)), right_weight_dual(gr3(M)),
                      right_weight_dual(gr2(M)), right_weight_dual(gr1(M)),
                      right_weight_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// left_bulk_dual = lcmpl(bulk(A))
//
// -> left complement operation applied to the bulk
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
constexpr PScalar3dp<T> left_bulk_dual(Scalar3dp<T> s)
{
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> left_bulk_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(-v.x, -v.y, -v.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> left_bulk_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> left_bulk_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), t.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> left_bulk_dual([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> left_bulk_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(left_bulk_dual(gr4(M)), left_bulk_dual(gr2(M)),
                        left_bulk_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> left_bulk_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(left_bulk_dual(gr3(M)), left_bulk_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> left_bulk_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(left_bulk_dual(gr4(M)), left_bulk_dual(gr3(M)),
                      left_bulk_dual(gr2(M)), left_bulk_dual(gr1(M)),
                      left_bulk_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// left_weight_dual(A) = = lcmpl(weight(A))
//
// -> left complement operation applied to the weight
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
constexpr PScalar3dp<T> left_weight_dual([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr TriVec3dp<T> left_weight_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), -v.w);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> left_weight_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec3dp<T> left_weight_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(t.x, t.y, t.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> left_weight_dual(PScalar3dp<T> ps)
{
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_E<T> left_weight_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(left_weight_dual(gr4(M)), left_weight_dual(gr2(M)),
                        left_weight_dual(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp_U<T> left_weight_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(left_weight_dual(gr3(M)), left_weight_dual(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
constexpr MVec3dp<T> left_weight_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(left_weight_dual(gr4(M)), left_weight_dual(gr3(M)),
                      left_weight_dual(gr2(M)), left_weight_dual(gr1(M)),
                      left_weight_dual(gr0(M)));
}

} // namespace hd::ga::pga