#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_foundation.hpp"        // ga library headers and infrastructure
#include "detail/type_t/ga_dualnum4ds.hpp" // DualNum4ds
#include "detail/type_t/ga_mvec4ds.hpp"    // 4ds multivector types
#include "detail/type_t/ga_type4ds.hpp"    // Point4ds, Vector3d, Point3d, Line3d, Plane3d

#include "detail/ga_error_handling.hpp"


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga4ds basic operations:
//
// - gr_inv()                      -> grade involution
// - rev()                         -> reversion
// - rrev()                        -> regressive reversion
// - conj()                        -> Clifford conjugation
//
// - r_cmpl()                      -> right complement (non-metric)
// - l_cmpl()                      -> left complement (non-metric)
//
// - nrm_sq()                      -> squared norm (based on metric)
// - nrm()                         -> magnitude = sqrt(|nrm_sq()|), always >= 0
//
// - is_timelike()                 -> true, if nrm_sq() > 0    (timelike part dominates)
// - is_spacelike()                -> true, if nrm_sq() < 0    (spacelike part dominates)
// - is_lightlike()                -> true, if nrm_sq() == 0   (lightlike)
//
// - normalize()                   -> return normalized object
//                                    nrm_sq scaled to +1.0, -1.0, 0.0 for
//                                    timelike, spacelike, lightlike objects
//
// - r_dual()                      -> right metric dual
// - l_dual()                      -> left metric dual
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for r = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> gr_inv(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> gr_inv(Vec4ds<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> gr_inv(BiVec4ds<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> gr_inv(TriVec4ds<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> gr_inv(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> gr_inv(MVec4ds_E<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    return M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> gr_inv(MVec4ds_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> gr_inv(MVec4ds<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    // grade 1 and 3: sign reversal
    return MVec4ds<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)),
                      gr_inv(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> rev(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> rev(Vec4ds<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> rev(BiVec4ds<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> rev(TriVec4ds<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> rev(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> rev(MVec4ds_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec4ds_E<T>(rev(gr0(M)), rev(gr2(M)), rev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> rev(MVec4ds_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec4ds_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> rev(MVec4ds<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    // grade 4: no sign change
    return MVec4ds<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)), rev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// regressive reversion operation: reverse of complement + backtrafo
// rrev(A_r) = l_cmpl((-1)^(r*(r-1)/2) r_cmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=4, r = 0, 1, 2, 3, ...: + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> rrev(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> rrev(Vec4ds<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> rrev(BiVec4ds<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> rrev(TriVec4ds<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> rrev(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> rrev(MVec4ds_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec4ds_E<T>(rrev(gr0(M)), rrev(gr2(M)), rrev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> rrev(MVec4ds_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec4ds_U<T>(rrev(gr1(M)), rrev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> rrev(MVec4ds<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec4ds<T>(rrev(gr0(M)), rrev(gr1(M)), rrev(gr2(M)), rrev(gr3(M)),
                      rrev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation (identical sign pattern to regressive reverse for pga4ds):
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> conj(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> conj(Vec4ds<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> conj(BiVec4ds<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> conj(TriVec4ds<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> conj(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> conj(MVec4ds_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec4ds_E<T>(conj(gr0(M)), conj(gr2(M)), conj(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> conj(MVec4ds_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec4ds_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> conj(MVec4ds<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec4ds<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)),
                      conj(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// 4ds complement operations
// (the concept of complement is defined w.r.t. the outer product)
////////////////////////////////////////////////////////////////////////////////

// If M represents the subspace B of the blade u then complement of M represents
// the subspace orthorgonal to B.
// The complement exchanges basis vectors which are in the k-blade u with
// the basis vectors which are NOT contained in the k-blade u,
// and are needed to fill the space completely up to the corresponding pseudoscalar.
//
// left complement:  l_cmpl(u) ^ u  = I_4ds = g1^g2^g3^g4
// right complement: u ^ r_cmpl(u)  = I_4ds = g1^g2^g3^g4
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_cmpl(u) = r_cmpl(u) = cmpl(u)

////////////////////////////////////////////////////////////////////////////////
// right complements
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> r_cmpl(Scalar4ds<T> s)
{
    // u ^ r_cmpl(u) = g1^g2^g3^g4
    // u = 1:
    //     1 ^ r_cmpl(u) = g1^g2^g3^g4 => r_cmpl(u) = s g1^g2^g3^g4
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> r_cmpl(Vec4ds<T> const& v)
{
    // u ^ r_compl(u) = g1^g2^g3^g4
    // u = v.x g1 + v.y g2 + v.z g3 + v.w g4:
    //     u ^ r_cmpl(u) = g1^g2^g3^g4 =>
    //     u = g1 => r_cmpl(u) =  g234
    //     u = g2 => r_cmpl(u) =  g314
    //     u = g3 => r_cmpl(u) =  g124
    //     u = g4 => r_cmpl(u) =  -g123
    return TriVec4ds<T>(v.x, v.y, v.z, -v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> r_cmpl(BiVec4ds<T> const& B)
{
    // u ^ r_compl(u) = g1^g2^g3^g4
    // u = B.vx g14 + B.vy g24 + B.vz g34 + B.mx g23 + B.my g31 + B.mz g12:
    //     u ^ r_cmpl(u) = g1^g2^g3^g4 =>
    //     u = g14 => r_cmpl(u) = g23
    //     u = g24 => r_cmpl(u) = g31
    //     u = g34 => r_cmpl(u) = g12
    //     u = g23 => r_cmpl(u) = g14
    //     u = g31 => r_cmpl(u) = g24
    //     u = g12 => r_cmpl(u) = g34
    return BiVec4ds<T>(B.mx, B.my, B.mz, B.vx, B.vy, B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> r_cmpl(TriVec4ds<T> const& t)
{
    // u ^ r_compl(u) = g1^g2^g3^g4
    // u = t.x g234 + t.y g314 + t.z g124 + t.w g123:
    //     u ^ r_cmpl(u) = g1^g2^g3^g4 =>
    //     u = g234 => r_cmpl(u) = -g1
    //     u = g314 => r_cmpl(u) = -g2
    //     u = g124 => r_cmpl(u) = -g3
    //     u = g123 => r_cmpl(u) =  g4
    return Vec4ds<T>(-t.x, -t.y, -t.z, t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> r_cmpl(PScalar4ds<T> ps)
{
    // u ^ r_compl(u) = g1^g2^g3^g4
    // u = g1^g2^g3^g4:
    //     g1^g2^g3^g4 ^ r_cmpl(u) = g1^g2^g3^g4 => r_cmpl(u) = ps 1
    return Scalar4ds<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> r_cmpl(MVec4ds_U<T> const& M)
{
    // use the component complements directly
    return MVec4ds_U<T>(r_cmpl(gr3(M)), r_cmpl(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> r_cmpl(MVec4ds_E<T> const& M)
{
    // use the component complements directly
    return MVec4ds_E<T>(r_cmpl(gr4(M)), r_cmpl(gr2(M)), r_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> r_cmpl(MVec4ds<T> const& M)
{
    // use the component complements directly
    return MVec4ds<T>(r_cmpl(gr4(M)), r_cmpl(gr3(M)), r_cmpl(gr2(M)), r_cmpl(gr1(M)),
                      r_cmpl(gr0(M)));
}

////////////////////////////////////////////////////////////////////////////////
// left complements
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> l_cmpl(Scalar4ds<T> s)
{
    // l_cmpl(u) ^ u = g1^g2^g3^g4
    // u = 1:
    //     l_cmpl(u) ^ 1 = g1^g2^g3^g4 => l_cmpl(u) = s g1^g2^g3^g4
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> l_cmpl(Vec4ds<T> const& v)
{
    // l_cmpl(u) ^ u = g1^g2^g3^g4
    // u = v.x g1 + v.y g2 + v.z g3 + v.w g4:
    //     l_cmpl(u) ^ u = g1^g2^g3^g4 =>
    //     u = g1 => l_cmpl(u) = -g234
    //     u = g2 => l_cmpl(u) = -g314
    //     u = g3 => l_cmpl(u) = -g124
    //     u = g4 => l_cmpl(u) =  g123
    return TriVec4ds<T>(-v.x, -v.y, -v.z, v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> l_cmpl(BiVec4ds<T> const& B)
{
    // l_cmpl(u) ^ u = g1^g2^g3^g4
    // u = B.vx g14 + B.vy g24 + B.vz g34 + B.mx g23 + B.my g31 + B.mz g12:
    //     l_cmpl(u) ^ u = g1^g2^g3^g4 =>
    //     u = g14 => l_cmpl(u) = g23
    //     u = g24 => l_cmpl(u) = g31
    //     u = g34 => l_cmpl(u) = g12
    //     u = g23 => l_cmpl(u) = g14
    //     u = g31 => l_cmpl(u) = g24
    //     u = g12 => l_cmpl(u) = g34
    return BiVec4ds<T>(B.mx, B.my, B.mz, B.vx, B.vy, B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> l_cmpl(TriVec4ds<T> const& t)
{
    // l_cmpl(u) ^ u = g1^g2^g3^g4
    // u = t.x g234 + t.y g314 + t.z g124 + t.w g123:
    //     l_cmpl(u) ^ u = g1^g2^g3^g4 =>
    //     u = g234 => l_cmpl(u) =  g1
    //     u = g314 => l_cmpl(u) =  g2
    //     u = g124 => l_cmpl(u) =  g3
    //     u = g123 => l_cmpl(u) = -g4
    return Vec4ds<T>(t.x, t.y, t.z, -t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> l_cmpl(PScalar4ds<T> ps)
{
    // l_cmpl(u) ^ u = g1^g2^g3^g4
    // u = g1^g2^g3^g4:
    //     l_cmpl(u) ^ g1^g2^g3^g4 = g1^g2^g3^g4 => l_cmpl(u) = ps 1
    return Scalar4ds<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> l_cmpl(MVec4ds_U<T> const& M)
{
    // use the component complements directly
    return MVec4ds_U<T>(l_cmpl(gr3(M)), l_cmpl(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> l_cmpl(MVec4ds_E<T> const& M)
{
    // use the component complements directly
    return MVec4ds_E<T>(l_cmpl(gr4(M)), l_cmpl(gr2(M)), l_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> l_cmpl(MVec4ds<T> const& M)
{
    // use the component complements directly
    return MVec4ds<T>(l_cmpl(gr4(M)), l_cmpl(gr3(M)), l_cmpl(gr2(M)), l_cmpl(gr1(M)),
                      l_cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// nrm_sq(u): return the squared norm of u
//            (a scalar value resulting from dot(u,u), i.e. reflects the metric)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(Scalar4ds<T> s)
{
    return T(s) * T(s);
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(Vec4ds<T> const& v)
{
    return -v.x * v.x - v.y * v.y - v.z * v.z + v.w * v.w;
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(BiVec4ds<T> const& B)
{
    return B.vx * B.vx + B.vy * B.vy + B.vz * B.vz - B.mx * B.mx - B.my * B.my -
           B.mz * B.mz;
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(TriVec4ds<T> const& t)
{
    return -t.x * t.x - t.y * t.y - t.z * t.z + t.w * t.w;
}


template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(PScalar4ds<T> ps)
{
    return -T(ps) * T(ps);
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec4ds_E<T> const& M)
{
    return nrm_sq(gr0(M)) + nrm_sq(gr2(M)) + nrm_sq(gr4(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec4ds_U<T> const& M)
{
    return nrm_sq(gr1(M)) + nrm_sq(gr3(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec4ds<T> const& M)
{
    return nrm_sq(gr0(M)) + nrm_sq(gr1(M)) + nrm_sq(gr2(M)) + nrm_sq(gr3(M)) +
           nrm_sq(gr4(M));
}


////////////////////////////////////////////////////////////////////////////////
// causal character of a k-vector (grades 1..n-1: Vec, BiVec, TriVec), based on
// the extended metric of G(1,3,0) in this library's (-,-,-,+) convention.
//
// Each basis blade carries a metric indicator G(e) = dot(e,e) in {-1,+1}, with
// +1 marking a timelike and -1 a spacelike basis blade:
//   timelike  (G = +1): g4 ; g14, g24, g34 ; g123
//   spacelike (G = -1): g1, g2, g3 ; g23, g31, g12 ; g234, g314, g124
// so within every grade 1..n-1 some basis blades are timelike, others spacelike.
//
// For a general k-vector the character is decided by which contribution
// dominates, which the metric quadratic form nrm_sq reports directly:
//   nrm_sq(u) = dot(u,u) = sum_(timelike) u_i^2  -  sum_(spacelike) u_i^2
// hence
//   is_timelike(u)  == true  if nrm_sq(u) > 0    (timelike part dominates)
//   is_spacelike(u) == true  if nrm_sq(u) < 0    (spacelike part dominates)
//   is_lightlike(u) == true  if nrm_sq(u) == 0   (balanced; on the light cone)
//
// The three predicates are mutually exclusive and exhaustive, and the rule is
// uniform across grades 1..n-1. The scalar (grade 0) and the pseudoscalar
// (grade n) are excluded by design; mixed-grade multivectors are omitted too.
////////////////////////////////////////////////////////////////////////////////

// is_timelike(u): nrm_sq(u) > 0
template <typename T>
    requires(numeric_type<T>)
constexpr bool is_timelike(Vec4ds<T> const& v)
{
    return nrm_sq(v) > T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_timelike(BiVec4ds<T> const& B)
{
    return nrm_sq(B) > T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_timelike(TriVec4ds<T> const& t)
{
    return nrm_sq(t) > T(0.0);
}

// is_spacelike(u): nrm_sq(u) < 0
template <typename T>
    requires(numeric_type<T>)
constexpr bool is_spacelike(Vec4ds<T> const& v)
{
    return nrm_sq(v) < T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_spacelike(BiVec4ds<T> const& B)
{
    return nrm_sq(B) < T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_spacelike(TriVec4ds<T> const& t)
{
    return nrm_sq(t) < T(0.0);
}

// is_lightlike(u): nrm_sq(u) == 0 (null)
template <typename T>
    requires(numeric_type<T>)
constexpr bool is_lightlike(Vec4ds<T> const& v)
{
    return nrm_sq(v) == T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_lightlike(BiVec4ds<T> const& B)
{
    return nrm_sq(B) == T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_lightlike(TriVec4ds<T> const& t)
{
    return nrm_sq(t) == T(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// nrm(u): Minkowski magnitude = sqrt(|nrm_sq(u)|), always >= 0.
//
// Real-valued for every element and exactly 0 for lightlike (null) blades.
// Dividing a non-null grade 1..n-1 blade by nrm() scales it to nrm_sq == +/-1,
// preserving its causal character (spacelike -> -1, timelike -> +1). The sign
// itself is dropped here (it is kept in nrm_sq); recover it via the is_*
// predicates.
//
// NOTE: for mixed-grade multivectors nrm_sq is a signed grade sum that need not
// be a meaningful "length"; nrm() is still defined mechanically there.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(Scalar4ds<T> s)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(s))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(Vec4ds<T> const& v)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(v))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(BiVec4ds<T> const& B)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(B))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(TriVec4ds<T> const& t)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(t))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(PScalar4ds<T> ps)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(ps))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(MVec4ds_E<T> const& M)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(M))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(MVec4ds_U<T> const& M)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(M))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(MVec4ds<T> const& M)
{
    return Scalar4ds<T>(std::sqrt(std::abs(nrm_sq(M))));
}


////////////////////////////////////////////////////////////////////////////////
// normalization operations: scale an object to unit magnitude (nrm == 1.0),
// i.e. nrm_sq == +1 for timelike and -1 for spacelike blades; lightlike blades
// (nrm_sq == 0) are returned unchanged.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
inline Scalar4ds<T> normalize(Scalar4ds<T> s)
{
    T m = to_val(nrm(s));
    hd::ga::detail::check_normalization<T>(m, "scalar (4ds)");
    T scale = T(1.0) / m;
    return Scalar4ds<T>(scale * T(s));
}

template <typename T>
    requires(numeric_type<T>)
inline Vec4ds<T> normalize(Vec4ds<T> const& v)
{
    // lightlike (null): nrm_sq == 0 already; cannot be scaled to +/-1
    if (is_lightlike(v)) return v;
    // spacelike -> -1, timelike -> +1; dividing by nrm = sqrt(|nrm_sq|)
    // yields the correctly signed unit blade automatically
    T n = to_val(nrm(v));
    hd::ga::detail::check_normalization<T>(n, "vector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * v;
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> normalize(BiVec4ds<T> const& B)
{
    // lightlike (null): nrm_sq == 0 already; cannot be scaled to +/-1
    if (is_lightlike(B)) return B;
    // spacelike -> -1, timelike -> +1; dividing by nrm = sqrt(|nrm_sq|)
    // yields the correctly signed unit blade automatically
    T n = to_val(nrm(B));
    hd::ga::detail::check_normalization<T>(n, "bivector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * B;
}

template <typename T>
    requires(numeric_type<T>)
inline TriVec4ds<T> normalize(TriVec4ds<T> const& t)
{
    // lightlike (null): nrm_sq == 0 already; cannot be scaled to +/-1
    if (is_lightlike(t)) return t;
    // spacelike -> -1, timelike -> +1; dividing by nrm = sqrt(|nrm_sq|)
    // yields the correctly signed unit blade automatically
    T n = to_val(nrm(t));
    hd::ga::detail::check_normalization<T>(n, "trivector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * t;
}

// NOTE: normalization of mixed-grade multivectors is disabled for now. In an
// indefinite (Minkowski) metric nrm_sq is a signed grade sum and does not give
// a well-defined blade norm, so the timelike/spacelike/lightlike case
// distinction does not apply.


////////////////////////////////////////////////////////////////////////////////
// metric right dual: r_dual(A) = r_cmpl( metric * A )
//
// -> right complement operation after multiplication with the metric
//
////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Lengyel, "PGA illuminated")
// is defined w.r.t. the outer product, because the complement is defined w.r.t.
// the outer product.
//
// The only difference between the dual and the complement is the fact that
// for the dual the argument is multiplied with the extended metric before
// taking the complement.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> r_dual(Scalar4ds<T> s)
{
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> r_dual(Vec4ds<T> const& v)
{
    return TriVec4ds<T>(-v.x, -v.y, -v.z, -v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> r_dual(BiVec4ds<T> const& B)
{
    return BiVec4ds<T>(-B.mx, -B.my, -B.mz, B.vx, B.vy, B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> r_dual(TriVec4ds<T> const& t)
{
    return Vec4ds<T>(t.x, t.y, t.z, t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> r_dual(PScalar4ds<T> ps)
{
    return Scalar4ds<T>(-T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> r_dual(MVec4ds_E<T> const& M)
{
    return MVec4ds_E<T>(r_dual(gr4(M)), r_dual(gr2(M)), r_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> r_dual(MVec4ds_U<T> const& M)
{
    return MVec4ds_U<T>(r_dual(gr3(M)), r_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> r_dual(MVec4ds<T> const& M)
{
    return MVec4ds<T>(r_dual(gr4(M)), r_dual(gr3(M)), r_dual(gr2(M)), r_dual(gr1(M)),
                      r_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// metric left dual: l_dual(A) = l_cmpl( metric * A )
//
// -> left complement operation after multiplication with the metric
//
////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Lengyel, "PGA illuminated")
// is defined w.r.t. the outer product, because the complement is defined w.r.t.
// the outer product.
//
// The only difference between the dual and the complement is the fact that
// for the dual the argument is multiplied with the extended metric before
// taking the complement.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> l_dual(Scalar4ds<T> s)
{
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> l_dual(Vec4ds<T> const& v)
{
    return TriVec4ds<T>(v.x, v.y, v.z, v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> l_dual(BiVec4ds<T> const& B)
{
    return BiVec4ds<T>(-B.mx, -B.my, -B.mz, B.vx, B.vy, B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> l_dual(TriVec4ds<T> const& t)
{
    return Vec4ds<T>(-t.x, -t.y, -t.z, -t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> l_dual(PScalar4ds<T> ps)
{
    return Scalar4ds<T>(-T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> l_dual(MVec4ds_E<T> const& M)
{
    return MVec4ds_E<T>(l_dual(gr4(M)), l_dual(gr2(M)), l_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> l_dual(MVec4ds_U<T> const& M)
{
    return MVec4ds_U<T>(l_dual(gr3(M)), l_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> l_dual(MVec4ds<T> const& M)
{
    return MVec4ds<T>(l_dual(gr4(M)), l_dual(gr3(M)), l_dual(gr2(M)), l_dual(gr1(M)),
                      l_dual(gr0(M)));
}

} // namespace hd::ga::sta