#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_foundation.hpp"     // ga library headers and infrastructure
#include "detail/type_t/ga_mvec4ds.hpp" // 4ds multivector types
#include "detail/type_t/ga_type4ds.hpp" // Point4ds, Vector3d, Point3d, Line3d, Plane3d

#include "detail/ga_error_handling.hpp"


namespace hd::ga::detail {

// Geometric square B^2 = gr0(X * X) of an sta4ds k-blade (grades 1..n-1), computed
// directly from the vector signature (g1,g2,g3 spacelike: -1; g4 timelike: +1).
//
// This is the blade's GEOMETRIC square -- the physically meaningful quantity for causal
// character: a boost (timelike) plane has B^2 > 0, a rotation (spacelike) plane B^2 < 0.
// It is intentionally DISTINCT from nrm_sq, which under the P-unify metric is the
// reverse-norm <X ~X>_0 = sigma(k) * B^2 (equal for vectors, sign-flipped for
// bi-/trivectors). Kept here so the causal-character predicates do not depend on nrm_sq
// and stay invariant under the metric's reverse-norm convention.
//
// MAINTENANCE: these are a hand-written transcription of gr0(X * X), the scalar part of
// the GEOMETRIC PRODUCT of the blade with itself. They are spelled out via components
// only because operator* is not yet available in this header. They MUST stay in sync
// with the geometric product (operator*): if the signature, the basis ordering, or the
// gpr generation changes, these signs change too. If causal-character / rotor tests
// (is_timelike/spacelike/lightlike, exp, get_rotor/get_boost) start failing
// unexpectedly, check these against `gr0(X * X)` first (e.g. assert
// sta4ds_geom_sq(X) == value_t(gr0(X * X)) per basis blade).
template <typename T>
    requires(numeric_type<T>)
constexpr T sta4ds_geom_sq(Vec4ds<T> const& v)
{
    return -v.x * v.x - v.y * v.y - v.z * v.z + v.w * v.w;
}

template <typename T>
    requires(numeric_type<T>)
constexpr T sta4ds_geom_sq(BiVec4ds<T> const& B)
{
    return B.vx * B.vx + B.vy * B.vy + B.vz * B.vz - B.mx * B.mx - B.my * B.my -
           B.mz * B.mz;
}

template <typename T>
    requires(numeric_type<T>)
constexpr T sta4ds_geom_sq(TriVec4ds<T> const& t)
{
    return -t.x * t.x - t.y * t.y - t.z * t.z + t.w * t.w;
}

} // namespace hd::ga::detail


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
// - nrm_sq()                      -> squared norm (the metric reverse-norm <X ~X>_0,
//                                    i.e. the extended metric P; see causal-character
//                                    note below for how this differs from B^2)
// - nrm()                         -> magnitude = sqrt(|nrm_sq()|), always >= 0
//
// - is_timelike()                 -> true, if B^2 = gr0(X*X) > 0   (timelike plane)
// - is_spacelike()                -> true, if B^2 = gr0(X*X) < 0   (spacelike plane)
// - is_lightlike()                -> true, if B^2 = gr0(X*X) == 0  (lightlike/null)
//                                    (causal character reads the GEOMETRIC square B^2,
//                                    not nrm_sq -- see detail::sta4ds_geom_sq)
//
// - normalize()                   -> return normalized object
//                                    nrm scales the blade to nrm_sq == +/-1 (the sign
//                                    follows the reverse-norm); null blades unchanged
//
// - r_dual()                      -> right metric dual
// - l_dual()                      -> left metric dual
// - r_undual()                    -> inverse of r_dual (metric on the dual side)
// - l_undual()                    -> inverse of l_dual (metric on the dual side)
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
    return -B.vx * B.vx - B.vy * B.vy - B.vz * B.vz + B.mx * B.mx + B.my * B.my +
           B.mz * B.mz;
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(TriVec4ds<T> const& t)
{
    return t.x * t.x + t.y * t.y + t.z * t.z - t.w * t.w;
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
// causal character of a k-vector (grades 1..n-1: Vec, BiVec, TriVec) in this
// library's (-,-,-,+) convention for G(1,3,0).
//
// Causal character is set by the GEOMETRIC square B^2 = gr0(X*X) (the blade times
// itself under the geometric product), NOT by nrm_sq. For each basis blade B^2 is
// +1 (timelike) or -1 (spacelike):
//   timelike  (B^2 = +1): g4 ; g14, g24, g34 ; g123
//   spacelike (B^2 = -1): g1, g2, g3 ; g23, g31, g12 ; g234, g314, g124
// so within every grade 1..n-1 some basis blades are timelike, others spacelike.
//
// For a general k-vector the character is decided by which contribution dominates:
//   B^2 = sum_(timelike) u_i^2  -  sum_(spacelike) u_i^2
// hence
//   is_timelike(u)  == true  if B^2 > 0    (timelike part dominates)
//   is_spacelike(u) == true  if B^2 < 0    (spacelike part dominates)
//   is_lightlike(u) == true  if B^2 == 0   (balanced; on the light cone)
//
// IMPORTANT: under the P-unify metric nrm_sq(u) = dot(u,u) is the reverse-norm
// <u ~u>_0 = sigma(k) * B^2 -- it AGREES with B^2 for vectors (sigma(1) = +1) but is
// SIGN-FLIPPED for bi-/trivectors (sigma(2) = sigma(3) = -1). Causal character therefore
// reads B^2 via detail::sta4ds_geom_sq(), never nrm_sq (whose sign no longer tracks the
// timelike/spacelike split at grades 2,3).
//
// The three predicates are mutually exclusive and exhaustive, and the rule is
// uniform across grades 1..n-1. The scalar (grade 0) and the pseudoscalar
// (grade n) are excluded by design; mixed-grade multivectors are omitted too.
////////////////////////////////////////////////////////////////////////////////

// is_timelike(u): B^2 = gr0(u*u) > 0
template <typename T>
    requires(numeric_type<T>)
constexpr bool is_timelike(Vec4ds<T> const& v)
{
    return detail::sta4ds_geom_sq(v) > T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_timelike(BiVec4ds<T> const& B)
{
    return detail::sta4ds_geom_sq(B) > T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_timelike(TriVec4ds<T> const& t)
{
    return detail::sta4ds_geom_sq(t) > T(0.0);
}

// is_spacelike(u): B^2 = gr0(u*u) < 0
template <typename T>
    requires(numeric_type<T>)
constexpr bool is_spacelike(Vec4ds<T> const& v)
{
    return detail::sta4ds_geom_sq(v) < T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_spacelike(BiVec4ds<T> const& B)
{
    return detail::sta4ds_geom_sq(B) < T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_spacelike(TriVec4ds<T> const& t)
{
    return detail::sta4ds_geom_sq(t) < T(0.0);
}

// is_lightlike(u): B^2 = gr0(u*u) == 0 (null)
template <typename T>
    requires(numeric_type<T>)
constexpr bool is_lightlike(Vec4ds<T> const& v)
{
    return detail::sta4ds_geom_sq(v) == T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_lightlike(BiVec4ds<T> const& B)
{
    return detail::sta4ds_geom_sq(B) == T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr bool is_lightlike(TriVec4ds<T> const& t)
{
    return detail::sta4ds_geom_sq(t) == T(0.0);
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
    return BiVec4ds<T>(B.mx, B.my, B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> r_dual(TriVec4ds<T> const& t)
{
    return Vec4ds<T>(-t.x, -t.y, -t.z, -t.w);
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
    return BiVec4ds<T>(B.mx, B.my, B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> l_dual(TriVec4ds<T> const& t)
{
    return Vec4ds<T>(t.x, t.y, t.z, t.w);
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


////////////////////////////////////////////////////////////////////////////////
// metric un-duals (inverse duals): the genuine inverse of the metric dual, with the
// metric applied on the DUAL side rather than the primal side (Lengyel: "antidual").
//
//   r_undual(D) = l_cmpl( r_dual( l_cmpl(D) ) )   inverts r_dual: r_undual(r_dual(A)) ==
//   A l_undual(D) = r_cmpl( l_dual( r_cmpl(D) ) )   inverts l_dual: l_undual(l_dual(A))
//   == A
//
// Even-dimensional STA has two duals, hence two un-duals -- one per handedness. They
// coincide on the even grades (where l_dual == r_dual) and differ by a sign on the odd
// grades, exactly as l_dual / r_dual do.
//
// NOTE the naive composition l_dual(r_dual(A)) is NOT the identity: it equals det · A
// (= -A for G(1,3,0)). The un-dual puts the metric on the dual side instead and recovers
// A exactly, at every grade and for any non-degenerate signature (it depends on the
// metric only through G² = I). Defined here by composing the (metric) dual with the
// (metric-free) complements, so it inherits their correctness rather than re-deriving
// sign tables.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> r_undual(Scalar4ds<T> s)
{
    return l_cmpl(r_dual(l_cmpl(s)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> r_undual(Vec4ds<T> const& v)
{
    return l_cmpl(r_dual(l_cmpl(v)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> r_undual(BiVec4ds<T> const& B)
{
    return l_cmpl(r_dual(l_cmpl(B)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> r_undual(TriVec4ds<T> const& t)
{
    return l_cmpl(r_dual(l_cmpl(t)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> r_undual(PScalar4ds<T> ps)
{
    return l_cmpl(r_dual(l_cmpl(ps)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> r_undual(MVec4ds_E<T> const& M)
{
    return l_cmpl(r_dual(l_cmpl(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> r_undual(MVec4ds_U<T> const& M)
{
    return l_cmpl(r_dual(l_cmpl(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> r_undual(MVec4ds<T> const& M)
{
    return l_cmpl(r_dual(l_cmpl(M)));
}


template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> l_undual(Scalar4ds<T> s)
{
    return r_cmpl(l_dual(r_cmpl(s)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> l_undual(Vec4ds<T> const& v)
{
    return r_cmpl(l_dual(r_cmpl(v)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> l_undual(BiVec4ds<T> const& B)
{
    return r_cmpl(l_dual(r_cmpl(B)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> l_undual(TriVec4ds<T> const& t)
{
    return r_cmpl(l_dual(r_cmpl(t)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> l_undual(PScalar4ds<T> ps)
{
    return r_cmpl(l_dual(r_cmpl(ps)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> l_undual(MVec4ds_E<T> const& M)
{
    return r_cmpl(l_dual(r_cmpl(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> l_undual(MVec4ds_U<T> const& M)
{
    return r_cmpl(l_dual(r_cmpl(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> l_undual(MVec4ds<T> const& M)
{
    return r_cmpl(l_dual(r_cmpl(M)));
}

} // namespace hd::ga::sta