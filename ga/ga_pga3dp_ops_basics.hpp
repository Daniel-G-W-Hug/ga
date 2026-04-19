#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_foundation.hpp"        // ga library headers and infrastructure
#include "detail/type_t/ga_dualnum3dp.hpp" // DualNum3dp
#include "detail/type_t/ga_mvec3dp.hpp"    // 3dp multivector types
#include "detail/type_t/ga_type3dp.hpp"    // Point3dp, Vector3d, Point3d, Line3d, Plane3d

#include "detail/ga_error_handling.hpp"


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga3dp basic operations:
//
// - gr_inv()                      -> grade inversion
// - rev()                         -> reversion
// - rrev()                        -> regressive reversion
// - conj()                        -> conjugation
//
// - l_cmpl(), r_cmpl()            -> left and right complement
//
// - bulk(), weight()              -> return bulk and weight parts of objects
// - bulk_nrm_sq(), bulk_nrm()     -> return bulk norm
// - weight_nrm_sq(), weight_nrm() -> return weight norm
// - geom_nrm_sq(), geom_nrm()     -> return geometric norm
//
// - bulk_normalize()              -> return normalized object (bulk_nrm scaled to 1.0)
// - unitize()                     -> return unitized object (weight_nrm scaled to 1.0)
//
// - r_bulk_dual()                 -> return right bulk dual
// - r_weight_dual()               -> return right weight dual
// - l_bulk_dual()                 -> return right bulk dual
// - l_weight_dual()               -> return right weight dual
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> gr_inv(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> gr_inv(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> gr_inv(BiVec3dp<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> gr_inv(TriVec3dp<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> gr_inv(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> gr_inv(MVec3dp_E<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    return M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> gr_inv(MVec3dp_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(numeric_type<T>)
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
    requires(numeric_type<T>)
constexpr Scalar3dp<T> rev(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> rev(Vec3dp<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> rev(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> rev(TriVec3dp<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> rev(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> rev(MVec3dp_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec3dp_E<T>(rev(gr0(M)), rev(gr2(M)), rev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> rev(MVec3dp_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec3dp_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
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
// rrev(A_r) = l_cmpl((-1)^(r*(r-1)/2) r_cmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=4, r = 0, 1, 2, 3, ...: + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> rrev(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> rrev(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> rrev(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> rrev(TriVec3dp<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> rrev(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> rrev(MVec3dp_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec3dp_E<T>(rrev(gr0(M)), rrev(gr2(M)), rrev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> rrev(MVec3dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3dp_U<T>(rrev(gr1(M)), rrev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
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
    requires(numeric_type<T>)
constexpr Scalar3dp<T> conj(Scalar3dp<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> conj(Vec3dp<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> conj(BiVec3dp<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> conj(TriVec3dp<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> conj(PScalar3dp<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> conj(MVec3dp_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec3dp_E<T>(conj(gr0(M)), conj(gr2(M)), conj(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> conj(MVec3dp_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3dp_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
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
// left complement:  l_cmpl(u) ^ u  = I_3dp = e1^e2^e3^e4
// right complement: u ^ r_cmpl(u)  = I_3dp = e1^e2^e3^e4
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs
//
// complement operation with e1^e2^e3^e4 as the pseudoscalar

////////////////////////////////////////////////////////////////////////////////
// right complements
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> r_cmpl(Scalar3dp<T> s)
{
    // u ^ r_cmpl(u) = e1^e2^e3^e4
    // u = 1:
    //     1 ^ r_cmpl(u) = e1^e2^e3^e4 => r_cmpl(u) = s e1^e2^e3^e4
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> r_cmpl(Vec3dp<T> const& v)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = v.x e1 + v.y e2 + v.z e3 + v.w e4:
    //     u ^ r_cmpl(u) = e1^e2^e3^e4 =>
    //     u = e1 => r_cmpl(u) = e423
    //     u = e2 => r_cmpl(u) = e431
    //     u = e3 => r_cmpl(u) = e412
    //     u = e4 => r_cmpl(u) = e321
    return TriVec3dp<T>(v.x, v.y, v.z, v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> r_cmpl(BiVec3dp<T> const& B)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     u ^ r_cmpl(u) = e1^e2^e3^e4 =>
    //     u = e41 => r_cmpl(u) = -e23
    //     u = e42 => r_cmpl(u) = -e31
    //     u = e43 => r_cmpl(u) = -e12
    //     u = e23 => r_cmpl(u) = -e41
    //     u = e31 => r_cmpl(u) = -e42
    //     u = e12 => r_cmpl(u) = -e43
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> r_cmpl(TriVec3dp<T> const& t)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     u ^ r_cmpl(u) = e1^e2^e3^e4 =>
    //     u = e423 => r_cmpl(u) = -e1
    //     u = e431 => r_cmpl(u) = -e2
    //     u = e412 => r_cmpl(u) = -e3
    //     u = e321 => r_cmpl(u) = -e4
    return Vec3dp<T>(-t.x, -t.y, -t.z, -t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> r_cmpl(PScalar3dp<T> ps)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     e1^e2^e3^e4 ^ r_cmpl(u) = e1^e2^e3^e4 => r_cmpl(u) = ps 1
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> r_cmpl(MVec3dp_U<T> const& M)
{
    // use the component complements directly
    return MVec3dp_U<T>(r_cmpl(gr3(M)), r_cmpl(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> r_cmpl(MVec3dp_E<T> const& M)
{
    // use the component complements directly
    return MVec3dp_E<T>(r_cmpl(gr4(M)), r_cmpl(gr2(M)), r_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> r_cmpl(MVec3dp<T> const& M)
{
    // use the component complements directly
    return MVec3dp<T>(r_cmpl(gr4(M)), r_cmpl(gr3(M)), r_cmpl(gr2(M)), r_cmpl(gr1(M)),
                      r_cmpl(gr0(M)));
}

////////////////////////////////////////////////////////////////////////////////
// left complements
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> l_cmpl(Scalar3dp<T> s)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = 1:
    //     l_cmpl(u) ^ 1 = e1^e2^e3^e4 => l_cmpl(u) = s e1^e2^e3^e4
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> l_cmpl(Vec3dp<T> const& v)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = v.x e1 + v.y e2 + v.z e3 + v.w e4:
    //     l_cmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e1 => l_cmpl(u) = -e423
    //     u = e2 => l_cmpl(u) = -e431
    //     u = e3 => l_cmpl(u) = -e412
    //     u = e4 => l_cmpl(u) = -e321
    return TriVec3dp<T>(-v.x, -v.y, -v.z, -v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> l_cmpl(BiVec3dp<T> const& B)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     l_cmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e41 => l_cmpl(u) = -e23
    //     u = e42 => l_cmpl(u) = -e31
    //     u = e43 => l_cmpl(u) = -e12
    //     u = e23 => l_cmpl(u) = -e41
    //     u = e31 => l_cmpl(u) = -e42
    //     u = e12 => l_cmpl(u) = -e43
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> l_cmpl(TriVec3dp<T> const& t)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     l_cmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e423 => l_cmpl(u) = e1
    //     u = e431 => l_cmpl(u) = e2
    //     u = e412 => l_cmpl(u) = e3
    //     u = e321 => l_cmpl(u) = e4
    return Vec3dp<T>(t.x, t.y, t.z, t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> l_cmpl(PScalar3dp<T> ps)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     l_cmpl(u) ^ e1^e2^e3^e4 = e1^e2^e3^e4 => l_cmpl(u) = ps 1
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> l_cmpl(MVec3dp_U<T> const& M)
{
    // use the component complements directly
    return MVec3dp_U<T>(l_cmpl(gr3(M)), l_cmpl(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> l_cmpl(MVec3dp_E<T> const& M)
{
    // use the component complements directly
    return MVec3dp_E<T>(l_cmpl(gr4(M)), l_cmpl(gr2(M)), l_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> l_cmpl(MVec3dp<T> const& M)
{
    // use the component complements directly
    return MVec3dp<T>(l_cmpl(gr4(M)), l_cmpl(gr3(M)), l_cmpl(gr2(M)), l_cmpl(gr1(M)),
                      l_cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk: u_bulk = G u (with G as the metric)
////////////////////////////////////////////////////////////////////////////////
// By the multiplication with the metric G (which is degenerate in this case)
// the bulk selects the components of the blade argument which dot NOT contain e4.
//
// => bulk():
//    - does not contain a factor of e4 (the basis component with e4*e4 = 0)
//    - contains positional information about u
//    - zero bulk means that u contains the origin (e4)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk(Scalar3dp<T> s)
{
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> bulk(Vec3dp<T> const& v)
{
    return Vec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> bulk(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), B.mx, B.my, B.mz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> bulk(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> bulk([[maybe_unused]] PScalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> bulk(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(bulk(gr0(M)), bulk(gr2(M)), bulk(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> bulk(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(bulk(gr1(M)), bulk(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> bulk(MVec3dp<T> const& M)
{
    return MVec3dp<T>(bulk(gr0(M)), bulk(gr1(M)), bulk(gr2(M)), bulk(gr3(M)),
                      bulk(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight: u_weight = l_cmpl( G r_cmpl(u) ) = rG u
//         (with G as the metric and rG as the anti-metric as given by Lengyel)
////////////////////////////////////////////////////////////////////////////////
// By the multiplication with the metric G (which is degenerate in this case)
// the weight selects the components of the blade argument which do contain e3.
//
// => weight():
//    - does contain a factor of e4 (the basis component with e4*e4 = 0)
//    - contains directional information about u
//    - zero weight means that u is contained in the horizon = r_cmpl(e4) = e321_3dp
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> weight([[maybe_unused]] Scalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> weight(Vec3dp<T> const& v)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> weight(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(B.vx, B.vy, B.vz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> weight(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(t.x, t.y, t.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(ps);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> weight(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(weight(gr0(M)), weight(gr2(M)), weight(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> weight(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(weight(gr1(M)), weight(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> weight(MVec3dp<T> const& M)
{
    return MVec3dp<T>(weight(gr0(M)), weight(gr1(M)), weight(gr2(M)), weight(gr3(M)),
                      weight(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk_nrm(u): return the bulk norm of u (a scalar value)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq(Scalar3dp<T> s)
{
    return T(s) * T(s);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm(Scalar3dp<T> s)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(s)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq(Vec3dp<T> const& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm(Vec3dp<T> const& v)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(v)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq(BiVec3dp<T> const& B)
{
    return B.mx * B.mx + B.my * B.my + B.mz * B.mz;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm(BiVec3dp<T> const& B)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(B)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq(TriVec3dp<T> const& t)
{
    return t.w * t.w;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm(TriVec3dp<T> const& t)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(t)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq([[maybe_unused]] PScalar3dp<T>)
{
    return T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq(MVec3dp_E<T> const& M)
{
    return bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr2(M)) + bulk_nrm_sq(gr4(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm(MVec3dp_E<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq(MVec3dp_U<T> const& M)
{
    return bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr3(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm(MVec3dp_U<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T bulk_nrm_sq(MVec3dp<T> const& M)
{
    return bulk_nrm_sq(gr0(M)) + bulk_nrm_sq(gr1(M)) + bulk_nrm_sq(gr2(M)) +
           bulk_nrm_sq(gr3(M)) + bulk_nrm_sq(gr4(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> bulk_nrm(MVec3dp<T> const& M)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(M)));
}


////////////////////////////////////////////////////////////////////////////////
// weight_nrm(u): return the weight norm of u (a pseudoscalar value)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq([[maybe_unused]] Scalar3dp<T>)
{
    return T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight_nrm([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq(Vec3dp<T> const& v)
{
    return v.w * v.w;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight_nrm(Vec3dp<T> const& v)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(v)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq(BiVec3dp<T> const& B)
{
    return B.vx * B.vx + B.vy * B.vy + B.vz * B.vz;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight_nrm(BiVec3dp<T> const& B)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(B)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq(TriVec3dp<T> const& t)
{
    return t.x * t.x + t.y * t.y + t.z * t.z;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight_nrm(TriVec3dp<T> const& t)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(t)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq(PScalar3dp<T> ps)
{
    return T(ps) * T(ps);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight_nrm(PScalar3dp<T> ps)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(ps)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq(MVec3dp_E<T> const& M)
{
    return weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr2(M)) + weight_nrm_sq(gr4(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight_nrm(MVec3dp_E<T> const& M)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq(MVec3dp_U<T> const& M)
{
    return weight_nrm_sq(gr1(M)) + weight_nrm_sq(gr3(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dp<T> weight_nrm(MVec3dp_U<T> const& M)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T weight_nrm_sq(MVec3dp<T> const& M)
{
    return weight_nrm_sq(gr0(M)) + weight_nrm_sq(gr1(M)) + weight_nrm_sq(gr2(M)) +
           weight_nrm_sq(gr3(M)) + weight_nrm_sq(gr4(M));
}

template <typename T>
    requires(numeric_type<T>)
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
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(Scalar3dp<T> s)
{
    return bulk_nrm_sq(s) + weight_nrm_sq(s);
}

template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(Scalar3dp<T> s)
{
    return DualNum3dp<T>(bulk_nrm(s), weight_nrm(s));
}

// provide the distance of the point from the origin
template <typename T>
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(Vec3dp<T> const& v)
{
    return bulk_nrm_sq(v) + weight_nrm_sq(v);
}

// return geometric norm of vector
template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(Vec3dp<T> const& v)
{
    return DualNum3dp<T>(bulk_nrm(v), weight_nrm(v));
}

// provide the perpendicular distance of the line to the origin
template <typename T>
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(BiVec3dp<T> const& B)
{
    return bulk_nrm_sq(B) + weight_nrm_sq(B);
}

// return geometric norm of bivector
template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(BiVec3dp<T> const& B)
{
    return DualNum3dp<T>(bulk_nrm(B), weight_nrm(B));
}

// provide the perpendicular distance of the plane to the origin
template <typename T>
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(TriVec3dp<T> const& t)
{
    return bulk_nrm_sq(t) + weight_nrm_sq(t);
}

// return geometric norm of vector
template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(TriVec3dp<T> const& t)
{
    return DualNum3dp<T>(bulk_nrm(t), weight_nrm(t));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(PScalar3dp<T> ps)
{
    return bulk_nrm_sq(ps) + weight_nrm_sq(ps);
}

template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(PScalar3dp<T> ps)
{
    return DualNum3dp<T>(bulk_nrm(ps), weight_nrm(ps));
}

// even-grade multivector
template <typename T>
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(MVec3dp_E<T> const& M)
{
    return bulk_nrm_sq(M) + weight_nrm_sq(M);
}

template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(MVec3dp_E<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}

// odd-grade multivector
template <typename T>
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(MVec3dp_U<T> const& M)
{
    return bulk_nrm_sq(M) + weight_nrm_sq(M);
}

template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(MVec3dp_U<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}

// multivector
template <typename T>
    requires(numeric_type<T>)
constexpr T geom_nrm_sq(MVec3dp<T> const& M)
{
    return bulk_nrm_sq(M) + weight_nrm_sq(M);
}

template <typename T>
    requires(numeric_type<T>)
constexpr DualNum3dp<T> geom_nrm(MVec3dp<T> const& M)
{
    return DualNum3dp<T>(bulk_nrm(M), weight_nrm(M));
}


////////////////////////////////////////////////////////////////////////////////
// bulk normalization operations:
// return an object normalized to bulk_nrm(u) == 1.0
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
inline Scalar3dp<T> bulk_normalize(Scalar3dp<T> s)
{
    T m = to_val(bulk_nrm(s));
    hd::ga::detail::check_normalization<T>(m, "scalar (3dp)");
    T scale = T(1.0) / m;
    return Scalar3dp<T>(scale * T(s));
}

template <typename T>
    requires(numeric_type<T>)
inline Vec3dp<T> bulk_normalize(Vec3dp<T> const& v)
{
    T n = to_val(bulk_nrm(v));
    hd::ga::detail::check_normalization<T>(std::abs(n), "vector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * v;
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec3dp<T> bulk_normalize(BiVec3dp<T> const& B)
{
    T n = to_val(bulk_nrm(B));
    hd::ga::detail::check_normalization<T>(n, "bivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * B;
}

template <typename T>
    requires(numeric_type<T>)
inline TriVec3dp<T> bulk_normalize(TriVec3dp<T> const& t)
{
    T n = to_val(bulk_nrm(t));
    hd::ga::detail::check_normalization<T>(n, "trivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * t;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dp_E<T> bulk_normalize(MVec3dp_E<T> const& M)
{
    T n = to_val(bulk_nrm(M));
    hd::ga::detail::check_normalization<T>(n, "even-grade multivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dp_U<T> bulk_normalize(MVec3dp_U<T> const& M)
{
    T n = to_val(bulk_nrm(M));
    hd::ga::detail::check_normalization<T>(n, "odd-grade multivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dp<T> bulk_normalize(MVec3dp<T> const& M)
{
    T n = to_val(bulk_nrm(M));
    hd::ga::detail::check_normalization<T>(n, "multivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}


////////////////////////////////////////////////////////////////////////////////
// unitization operations:
// return an object normalized to weight_nrm(u) == 1.0
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
inline Vec3dp<T> unitize(Vec3dp<T> const& v)
{
    T n = v.w; // v.w == sign(v.w) * weight_nrm(v)
    hd::ga::detail::check_unitization<T>(std::abs(n), "vector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return Vec3dp<T>(v.x * scale, v.y * scale, v.z * scale, T(1.0));
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec3dp<T> unitize(BiVec3dp<T> const& B)
{
    T n = T(weight_nrm(B));
    hd::ga::detail::check_unitization<T>(n, "bivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * B;
}

template <typename T>
    requires(numeric_type<T>)
inline TriVec3dp<T> unitize(TriVec3dp<T> const& t)
{
    T n = T(weight_nrm(t));
    hd::ga::detail::check_unitization<T>(n, "trivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * t;
}

template <typename T>
    requires(numeric_type<T>)
inline PScalar3dp<T> unitize(PScalar3dp<T> ps)
{
    T n = T(weight_nrm(ps));
    hd::ga::detail::check_unitization<T>(n, "pseudoscalar (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return PScalar3dp<T>(scale * T(ps));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dp_E<T> unitize(MVec3dp_E<T> const& M)
{
    T n = T(weight_nrm(M));
    hd::ga::detail::check_unitization<T>(n, "even-grade multivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dp_U<T> unitize(MVec3dp_U<T> const& M)
{
    T n = T(weight_nrm(M));
    hd::ga::detail::check_unitization<T>(n, "odd-grade multivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dp<T> unitize(MVec3dp<T> const& M)
{
    T n = T(weight_nrm(M));
    hd::ga::detail::check_unitization<T>(n, "multivector (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}

template <typename T>
    requires numeric_type<T>
inline Point3dp<T> unitize(Point3dp<T> const& p)
{
    hd::ga::detail::check_unitization<T>(std::abs(p.w), "Point3dp");
    T scale = T(1.0) / p.w;
    return Point3dp<T>(p.x * scale, p.y * scale, p.z * scale, T(1.0));
}

template <typename T>
    requires(numeric_type<T>)
inline Line3d<T> unitize(Line3d<T> const& l)
{
    // unitization for a 3dp bivector means a normalized direction vector
    // sqrt((l.vx)^2 + (l.vy)^2 + (l.vz)^2) = 1
    T wn = T(weight_nrm(l));
    hd::ga::detail::check_unitization<T>(wn, "Line3d");
    T scale = T(1.0) / wn;
    return Line3d<T>(l.vx * scale, l.vy * scale, l.vz * scale, l.mx * scale, l.my * scale,
                     l.mz * scale);
}

template <typename T>
    requires(numeric_type<T>)
inline Plane3d<T> unitize(Plane3d<T> const& p)
{
    T wn = T(weight_nrm(p));
    hd::ga::detail::check_unitization<T>(wn, "Plane3d");
    T scale = T(1.0) / wn;
    return Plane3d<T>(p.x * scale, p.y * scale, p.z * scale, p.w * scale);
}

template <typename T>
    requires(numeric_type<T>)
inline DualNum3dp<T> unitize(DualNum3dp<T> const& D)
{
    T n = D.c1; // the pseudoscalar part is the weight_nrm part
    hd::ga::detail::check_unitization<T>(std::abs(n), "dual number (3dp)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * D;
}


////////////////////////////////////////////////////////////////////////////////
// r_bulk_dual(A) = r_cmpl(bulk(A)) = r_cmpl( metric * A )
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
    requires(numeric_type<T>)
constexpr PScalar3dp<T> r_bulk_dual(Scalar3dp<T> s)
{
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> r_bulk_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> r_bulk_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> r_bulk_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), -t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> r_bulk_dual([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> r_bulk_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(r_bulk_dual(gr4(M)), r_bulk_dual(gr2(M)), r_bulk_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> r_bulk_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(r_bulk_dual(gr3(M)), r_bulk_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> r_bulk_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(r_bulk_dual(gr4(M)), r_bulk_dual(gr3(M)), r_bulk_dual(gr2(M)),
                      r_bulk_dual(gr1(M)), r_bulk_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// r_weight_dual(A) = r_cmpl(weight(A)) = r_cmpl( regressive_metric * A )
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
    requires(numeric_type<T>)
constexpr PScalar3dp<T> r_weight_dual([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> r_weight_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> r_weight_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> r_weight_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(-t.x, -t.y, -t.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> r_weight_dual(PScalar3dp<T> ps)
{
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> r_weight_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(r_weight_dual(gr4(M)), r_weight_dual(gr2(M)),
                        r_weight_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> r_weight_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(r_weight_dual(gr3(M)), r_weight_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> r_weight_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(r_weight_dual(gr4(M)), r_weight_dual(gr3(M)), r_weight_dual(gr2(M)),
                      r_weight_dual(gr1(M)), r_weight_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// l_bulk_dual(A) = l_cmpl(bulk(A)) = l_cmpl( metric * A )
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
    requires(numeric_type<T>)
constexpr PScalar3dp<T> l_bulk_dual(Scalar3dp<T> s)
{
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> l_bulk_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(-v.x, -v.y, -v.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> l_bulk_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> l_bulk_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> l_bulk_dual([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> l_bulk_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(l_bulk_dual(gr4(M)), l_bulk_dual(gr2(M)), l_bulk_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> l_bulk_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(l_bulk_dual(gr3(M)), l_bulk_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> l_bulk_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(l_bulk_dual(gr4(M)), l_bulk_dual(gr3(M)), l_bulk_dual(gr2(M)),
                      l_bulk_dual(gr1(M)), l_bulk_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// l_weight_dual(A) = = l_cmpl(weight(A)) = l_cmpl( regressive_metric * A )
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
    requires(numeric_type<T>)
constexpr PScalar3dp<T> l_weight_dual([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dp<T> l_weight_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), -v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dp<T> l_weight_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dp<T> l_weight_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(t.x, t.y, t.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dp<T> l_weight_dual(PScalar3dp<T> ps)
{
    return Scalar3dp<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_E<T> l_weight_dual(MVec3dp_E<T> const& M)
{
    return MVec3dp_E<T>(l_weight_dual(gr4(M)), l_weight_dual(gr2(M)),
                        l_weight_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp_U<T> l_weight_dual(MVec3dp_U<T> const& M)
{
    return MVec3dp_U<T>(l_weight_dual(gr3(M)), l_weight_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dp<T> l_weight_dual(MVec3dp<T> const& M)
{
    return MVec3dp<T>(l_weight_dual(gr4(M)), l_weight_dual(gr3(M)), l_weight_dual(gr2(M)),
                      l_weight_dual(gr1(M)), l_weight_dual(gr0(M)));
}

} // namespace hd::ga::pga