#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_cga_types.hpp"  // cga3dc types incl. multivector wrappers
#include "detail/ga_foundation.hpp" // ga library headers and infrastructure

#include "detail/ga_error_handling.hpp" // unitize() zero-division checks


namespace hd::ga::cga {

////////////////////////////////////////////////////////////////////////////////
// provides cga3dc basic operations:
//
// - gr_inv()                      -> grade involution
// - rev()                         -> reversion
// - rrev()                        -> regressive reversion (== rev in 5d)
// - conj()                        -> Clifford conjugation
// - cconj()                       -> conformal conjugate (round part kept,
//                                    flat part negated)
//
// - cmpl()                        -> complement (non-metric; left == right in 5d)
//
// - dual()                        -> metric dual
// - antidual()                    -> metric antidual
//
// - round_bulk(), round_weight(), flat_bulk(), flat_weight()
//                                 -> the four-part bulk/weight split
// - round_bulk_nrm{,_sq}(), round_weight_nrm{,_sq}(),
//   flat_bulk_nrm{,_sq}(), flat_weight_nrm{,_sq}() -> the four part norms
// - center_nrm{,_sq}()            -> weighted distance origin <-> center
// - unitize()                     -> scale so the round weight norm equals one
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for r = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> gr_inv(Scalar3dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> gr_inv(Vec3dc<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> gr_inv(BiVec3dc<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> gr_inv(TriVec3dc<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> gr_inv(QuadVec3dc<T> const& Q)
{
    // grade 4: no sign change
    return Q;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> gr_inv(PScalar3dc<T> ps)
{
    // grade 5: sign reversal
    return -ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_E<T> gr_inv(MVec3dc_E<T> const& M)
{
    // grade 0, 2, 4: no sign change
    return M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_U<T> gr_inv(MVec3dc_U<T> const& M)
{
    // grade 1, 3, 5: sign reversal
    return -M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc<T> gr_inv(MVec3dc<T> const& M)
{
    // grade 0, 2, 4: no sign change
    // grade 1, 3, 5: sign reversal
    return MVec3dc<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)),
                      gr_inv(gr4(M)), gr_inv(gr5(M)));
}


////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> rev(Scalar3dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> rev(Vec3dc<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> rev(BiVec3dc<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> rev(TriVec3dc<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> rev(QuadVec3dc<T> const& Q)
{
    // grade 4: no sign change
    return Q;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> rev(PScalar3dc<T> ps)
{
    // grade 5: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_E<T> rev(MVec3dc_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec3dc_E<T>(rev(gr0(M)), rev(gr2(M)), rev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_U<T> rev(MVec3dc_U<T> const& M)
{
    // grade 1, 5: no sign change
    // grade 3: sign reversal
    return MVec3dc_U<T>(rev(gr1(M)), rev(gr3(M)), rev(gr5(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc<T> rev(MVec3dc<T> const& M)
{
    // grade 0, 1, 4, 5: no sign change
    // grade 2, 3: sign reversal
    return MVec3dc<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)), rev(gr4(M)),
                      rev(gr5(M)));
}


////////////////////////////////////////////////////////////////////////////////
// regressive reversion operation: reverse of complement + backtrafo
// rrev(A_r) = cmpl((-1)^(r*(r-1)/2) cmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=5, r = 0, 1, 2, 3, ...: + + - - + +
//
// For the odd-dimensional cga3dc this coincides with rev() at every grade
// (rev == rrev in 5d); both names are kept so that generic code reads the
// same across the algebras.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> rrev(Scalar3dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> rrev(Vec3dc<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> rrev(BiVec3dc<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> rrev(TriVec3dc<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> rrev(QuadVec3dc<T> const& Q)
{
    // grade 4: no sign change
    return Q;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> rrev(PScalar3dc<T> ps)
{
    // grade 5: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_E<T> rrev(MVec3dc_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec3dc_E<T>(rrev(gr0(M)), rrev(gr2(M)), rrev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_U<T> rrev(MVec3dc_U<T> const& M)
{
    // grade 1, 5: no sign change
    // grade 3: sign reversal
    return MVec3dc_U<T>(rrev(gr1(M)), rrev(gr3(M)), rrev(gr5(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc<T> rrev(MVec3dc<T> const& M)
{
    // grade 0, 1, 4, 5: no sign change
    // grade 2, 3: sign reversal
    return MVec3dc<T>(rrev(gr0(M)), rrev(gr1(M)), rrev(gr2(M)), rrev(gr3(M)),
                      rrev(gr4(M)), rrev(gr5(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> conj(Scalar3dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> conj(Vec3dc<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> conj(BiVec3dc<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> conj(TriVec3dc<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> conj(QuadVec3dc<T> const& Q)
{
    // grade 4: no sign change
    return Q;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> conj(PScalar3dc<T> ps)
{
    // grade 5: sign reversal
    return -ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_E<T> conj(MVec3dc_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec3dc_E<T>(conj(gr0(M)), conj(gr2(M)), conj(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_U<T> conj(MVec3dc_U<T> const& M)
{
    // grade 3: no sign change
    // grade 1, 5: sign reversal
    return MVec3dc_U<T>(conj(gr1(M)), conj(gr3(M)), conj(gr5(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc<T> conj(MVec3dc<T> const& M)
{
    // grade 0, 3, 4: no sign change
    // grade 1, 2, 5: sign reversal
    return MVec3dc<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)),
                      conj(gr4(M)), conj(gr5(M)));
}


////////////////////////////////////////////////////////////////////////////////
// conformal conjugate: keep the round part, negate the flat part (i.e. negate
// every component carrying a factor of the infinity dimension e5)
//
//     cconj(u) = round_bulk(u) + round_weight(u) - flat_bulk(u) - flat_weight(u)
//
// An involution. The conjugate of a round object has its center REFLECTED
// through the origin at the same distance and its radius adjusted such that
// for a unitized u
//
//     dot(u, cconj(u)) = squared distance origin <-> center
//
// which is the alternate form of the center norm. (Note: distinct from the
// Clifford conjugation conj() above, which is a grade involution.)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> cconj(Vec3dc<T> const& v)
{
    return Vec3dc<T>(v.x, v.y, v.z, v.w, -v.u);
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> cconj(BiVec3dc<T> const& B)
{
    return BiVec3dc<T>(B.vx, B.vy, B.vz, B.mx, B.my, B.mz, -B.px, -B.py, -B.pz, -B.pw);
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> cconj(TriVec3dc<T> const& t)
{
    return TriVec3dc<T>(-t.vx, -t.vy, -t.vz, -t.mx, -t.my, -t.mz, t.px, t.py, t.pz, t.pw);
}
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> cconj(QuadVec3dc<T> const& Q)
{
    return QuadVec3dc<T>(-Q.x, -Q.y, -Q.z, -Q.w, Q.u);
}


////////////////////////////////////////////////////////////////////////////////
// complement operations (non-metric)
//
// The complement maps a blade onto the blade that completes it to the
// pseudoscalar:
//
//     cmpl(u) ^ u = u ^ cmpl(u) = I_3dc
//
// For the odd-dimensional cga3dc the left and right complement coincide (as
// for ega3d, pga2dp), so there is a single cmpl() with cmpl(cmpl(u)) = u.
////////////////////////////////////////////////////////////////////////////////

// cga3dc cmpl :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> cmpl(Scalar3dc<T> s)
{
    return PScalar3dc<T>(T(s));
}

// cga3dc cmpl :: vec -> quadvec
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> cmpl(Vec3dc<T> const& v)
{
    return QuadVec3dc<T>(v.x, v.y, v.z, v.w, v.u);
}

// cga3dc cmpl :: bivec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> cmpl(BiVec3dc<T> const& B)
{
    return TriVec3dc<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz, -B.px, -B.py, -B.pz,
                        -B.pw);
}

// cga3dc cmpl :: trivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> cmpl(TriVec3dc<T> const& t)
{
    return BiVec3dc<T>(-t.mx, -t.my, -t.mz, -t.vx, -t.vy, -t.vz, -t.px, -t.py, -t.pz,
                       -t.pw);
}

// cga3dc cmpl :: quadvec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> cmpl(QuadVec3dc<T> const& Q)
{
    return Vec3dc<T>(Q.x, Q.y, Q.z, Q.w, Q.u);
}

// cga3dc cmpl :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> cmpl(PScalar3dc<T> ps)
{
    return Scalar3dc<T>(T(ps));
}

// cga3dc cmpl :: mv_e -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_U<T> cmpl(MVec3dc_E<T> const& M)
{
    return MVec3dc_U<T>(cmpl(gr4(M)), cmpl(gr2(M)), cmpl(gr0(M)));
}

// cga3dc cmpl :: mv_u -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_E<T> cmpl(MVec3dc_U<T> const& M)
{
    return MVec3dc_E<T>(cmpl(gr5(M)), cmpl(gr3(M)), cmpl(gr1(M)));
}

// cga3dc cmpl :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc<T> cmpl(MVec3dc<T> const& M)
{
    return MVec3dc<T>(cmpl(gr5(M)), cmpl(gr4(M)), cmpl(gr3(M)), cmpl(gr2(M)),
                      cmpl(gr1(M)), cmpl(gr0(M)));
}

////////////////////////////////////////////////////////////////////////////////
// dual (= complement of the metric image: dual(a) = cmpl(G a))
//
// The cga3dc metric is non-degenerate, so there is no bulk/weight dual split
// (as needed for the degenerate pga metrics): the extended metric G is a
// signed permutation that pairs every projective (w-bearing) blade with a
// round (u-bearing) partner blade, and the dual stays single-term.
////////////////////////////////////////////////////////////////////////////////

// cga3dc dual :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> dual(Scalar3dc<T> s)
{
    return PScalar3dc<T>(T(s));
}

// cga3dc dual :: vec -> quadvec
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> dual(Vec3dc<T> const& v)
{
    return QuadVec3dc<T>(v.x, v.y, v.z, -v.u, -v.w);
}

// cga3dc dual :: bivec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> dual(BiVec3dc<T> const& B)
{
    return TriVec3dc<T>(-B.mx, -B.my, -B.mz, -B.px, -B.py, -B.pz, -B.vx, -B.vy, -B.vz,
                        B.pw);
}

// cga3dc dual :: trivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> dual(TriVec3dc<T> const& t)
{
    return BiVec3dc<T>(t.px, t.py, t.pz, t.vx, t.vy, t.vz, t.mx, t.my, t.mz, -t.pw);
}

// cga3dc dual :: quadvec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> dual(QuadVec3dc<T> const& Q)
{
    return Vec3dc<T>(-Q.x, -Q.y, -Q.z, Q.u, Q.w);
}

// cga3dc dual :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> dual(PScalar3dc<T> ps)
{
    return Scalar3dc<T>(-T(ps));
}

// cga3dc dual :: mv_e -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_U<T> dual(MVec3dc_E<T> const& M)
{
    return MVec3dc_U<T>(dual(gr4(M)), dual(gr2(M)), dual(gr0(M)));
}

// cga3dc dual :: mv_u -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_E<T> dual(MVec3dc_U<T> const& M)
{
    return MVec3dc_E<T>(dual(gr5(M)), dual(gr3(M)), dual(gr1(M)));
}

// cga3dc dual :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc<T> dual(MVec3dc<T> const& M)
{
    return MVec3dc<T>(dual(gr5(M)), dual(gr4(M)), dual(gr3(M)), dual(gr2(M)),
                      dual(gr1(M)), dual(gr0(M)));
}

////////////////////////////////////////////////////////////////////////////////
// antidual (= complement of the anti-metric image: antidual(a) = cmpl(AG a))
//
// The antidual uses the metric anti-exomorphism AG (the extended metric built
// w.r.t. the antiwedge product). For cga3dc the anti-exomorphism equals -G
// (verified at rule-generation time and pinned in ga_cga_test), so the
// antidual is the negated dual:
//
//     antidual(a) = -dual(a)
////////////////////////////////////////////////////////////////////////////////

// cga3dc antidual :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> antidual(Scalar3dc<T> s)
{
    return PScalar3dc<T>(-T(s));
}

// cga3dc antidual :: vec -> quadvec
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> antidual(Vec3dc<T> const& v)
{
    return QuadVec3dc<T>(-v.x, -v.y, -v.z, v.u, v.w);
}

// cga3dc antidual :: bivec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> antidual(BiVec3dc<T> const& B)
{
    return TriVec3dc<T>(B.mx, B.my, B.mz, B.px, B.py, B.pz, B.vx, B.vy, B.vz, -B.pw);
}

// cga3dc antidual :: trivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> antidual(TriVec3dc<T> const& t)
{
    return BiVec3dc<T>(-t.px, -t.py, -t.pz, -t.vx, -t.vy, -t.vz, -t.mx, -t.my, -t.mz,
                       t.pw);
}

// cga3dc antidual :: quadvec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> antidual(QuadVec3dc<T> const& Q)
{
    return Vec3dc<T>(Q.x, Q.y, Q.z, -Q.u, -Q.w);
}

// cga3dc antidual :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> antidual(PScalar3dc<T> ps)
{
    return Scalar3dc<T>(T(ps));
}

// cga3dc antidual :: mv_e -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_U<T> antidual(MVec3dc_E<T> const& M)
{
    return MVec3dc_U<T>(antidual(gr4(M)), antidual(gr2(M)), antidual(gr0(M)));
}

// cga3dc antidual :: mv_u -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc_E<T> antidual(MVec3dc_U<T> const& M)
{
    return MVec3dc_E<T>(antidual(gr5(M)), antidual(gr3(M)), antidual(gr1(M)));
}

// cga3dc antidual :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec3dc<T> antidual(MVec3dc<T> const& M)
{
    return MVec3dc<T>(antidual(gr5(M)), antidual(gr4(M)), antidual(gr3(M)),
                      antidual(gr2(M)), antidual(gr1(M)), antidual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// bulk/weight split
//
// Every component belongs to exactly one of four parts, by whether the blade
// contains the infinity dimension e5 (flat part: yes, round part: no) and
// the origin e4 (weight: yes, bulk: no):
//
//     round bulk:   no e5, no e4      round weight: no e5, with e4
//     flat bulk:    with e5, no e4    flat weight:  with e5 and e4
//
// The round weight carries the homogeneous weight of a round object, so
// unitizing makes the round weight norm one.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> round_bulk(Vec3dc<T> const& v)
{
    return Vec3dc<T>(v.x, v.y, v.z, T(0.0), T(0.0)); // e1, e2, e3
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> round_bulk(BiVec3dc<T> const& B)
{
    return BiVec3dc<T>(T(0.0), T(0.0), T(0.0), B.mx, B.my, B.mz, T(0.0), T(0.0), T(0.0),
                       T(0.0)); // e23, e31, e12
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> round_bulk(TriVec3dc<T> const& t)
{
    return TriVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                        T(0.0), t.pw); // e321
}
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> round_bulk([[maybe_unused]] QuadVec3dc<T> const& Q)
{
    return QuadVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)); // (none at grade 4)
}
template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> round_weight(Vec3dc<T> const& v)
{
    return Vec3dc<T>(T(0.0), T(0.0), T(0.0), v.w, T(0.0)); // e4
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> round_weight(BiVec3dc<T> const& B)
{
    return BiVec3dc<T>(B.vx, B.vy, B.vz, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                       T(0.0)); // e41, e42, e43
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> round_weight(TriVec3dc<T> const& t)
{
    return TriVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), t.px, t.py, t.pz,
                        T(0.0)); // e423, e431, e412
}
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> round_weight(QuadVec3dc<T> const& Q)
{
    return QuadVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), Q.u); // e1234
}
template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> flat_bulk(Vec3dc<T> const& v)
{
    return Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), v.u); // e5
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> flat_bulk(BiVec3dc<T> const& B)
{
    return BiVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.px, B.py, B.pz,
                       T(0.0)); // e15, e25, e35
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> flat_bulk(TriVec3dc<T> const& t)
{
    return TriVec3dc<T>(T(0.0), T(0.0), T(0.0), t.mx, t.my, t.mz, T(0.0), T(0.0), T(0.0),
                        T(0.0)); // e235, e315, e125
}
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> flat_bulk(QuadVec3dc<T> const& Q)
{
    return QuadVec3dc<T>(T(0.0), T(0.0), T(0.0), Q.w, T(0.0)); // e3215
}
template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> flat_weight([[maybe_unused]] Vec3dc<T> const& v)
{
    return Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)); // (none at grade 1)
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> flat_weight(BiVec3dc<T> const& B)
{
    return BiVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                       T(0.0), B.pw); // e45
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> flat_weight(TriVec3dc<T> const& t)
{
    return TriVec3dc<T>(t.vx, t.vy, t.vz, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                        T(0.0)); // e415, e425, e435
}
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> flat_weight(QuadVec3dc<T> const& Q)
{
    return QuadVec3dc<T>(Q.x, Q.y, Q.z, T(0.0), T(0.0)); // e4235, e4315, e4125
}

////////////////////////////////////////////////////////////////////////////////
// the four part norms (Euclidean lengths of the four component groups) and the
// center norm
//
//     center_nrm(u) = sqrt( round_bulk_nrm^2 + flat_weight_nrm^2 )
//
// is the weighted distance between the origin and the object's center. (The
// radius norm is defined in the geometric layer -- it needs the product-layer
// antidot square rdot(u, u).)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm_sq(Vec3dc<T> const& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm_sq(BiVec3dc<T> const& B)
{
    return B.mx * B.mx + B.my * B.my + B.mz * B.mz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm_sq(TriVec3dc<T> const& t)
{
    return t.pw * t.pw;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm_sq([[maybe_unused]] QuadVec3dc<T> const& Q)
{
    return T(0.0);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm_sq(Vec3dc<T> const& v)
{
    return v.w * v.w;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm_sq(BiVec3dc<T> const& B)
{
    return B.vx * B.vx + B.vy * B.vy + B.vz * B.vz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm_sq(TriVec3dc<T> const& t)
{
    return t.px * t.px + t.py * t.py + t.pz * t.pz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm_sq(QuadVec3dc<T> const& Q)
{
    return Q.u * Q.u;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm_sq(Vec3dc<T> const& v)
{
    return v.u * v.u;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm_sq(BiVec3dc<T> const& B)
{
    return B.px * B.px + B.py * B.py + B.pz * B.pz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm_sq(TriVec3dc<T> const& t)
{
    return t.mx * t.mx + t.my * t.my + t.mz * t.mz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm_sq(QuadVec3dc<T> const& Q)
{
    return Q.w * Q.w;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm_sq([[maybe_unused]] Vec3dc<T> const& v)
{
    return T(0.0);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm_sq(BiVec3dc<T> const& B)
{
    return B.pw * B.pw;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm_sq(TriVec3dc<T> const& t)
{
    return t.vx * t.vx + t.vy * t.vy + t.vz * t.vz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm_sq(QuadVec3dc<T> const& Q)
{
    return Q.x * Q.x + Q.y * Q.y + Q.z * Q.z;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm(Vec3dc<T> const& v)
{
    return std::sqrt(round_bulk_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm(BiVec3dc<T> const& B)
{
    return std::sqrt(round_bulk_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm(TriVec3dc<T> const& t)
{
    return std::sqrt(round_bulk_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm(QuadVec3dc<T> const& Q)
{
    return std::sqrt(round_bulk_nrm_sq(Q));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm(Vec3dc<T> const& v)
{
    return std::sqrt(round_weight_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm(BiVec3dc<T> const& B)
{
    return std::sqrt(round_weight_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm(TriVec3dc<T> const& t)
{
    return std::sqrt(round_weight_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm(QuadVec3dc<T> const& Q)
{
    return std::sqrt(round_weight_nrm_sq(Q));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm(Vec3dc<T> const& v)
{
    return std::sqrt(flat_bulk_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm(BiVec3dc<T> const& B)
{
    return std::sqrt(flat_bulk_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm(TriVec3dc<T> const& t)
{
    return std::sqrt(flat_bulk_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm(QuadVec3dc<T> const& Q)
{
    return std::sqrt(flat_bulk_nrm_sq(Q));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm(Vec3dc<T> const& v)
{
    return std::sqrt(flat_weight_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm(BiVec3dc<T> const& B)
{
    return std::sqrt(flat_weight_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm(TriVec3dc<T> const& t)
{
    return std::sqrt(flat_weight_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm(QuadVec3dc<T> const& Q)
{
    return std::sqrt(flat_weight_nrm_sq(Q));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm_sq(Vec3dc<T> const& v)
{
    return round_bulk_nrm_sq(v) + flat_weight_nrm_sq(v);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm_sq(BiVec3dc<T> const& B)
{
    return round_bulk_nrm_sq(B) + flat_weight_nrm_sq(B);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm_sq(TriVec3dc<T> const& t)
{
    return round_bulk_nrm_sq(t) + flat_weight_nrm_sq(t);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm_sq(QuadVec3dc<T> const& Q)
{
    return round_bulk_nrm_sq(Q) + flat_weight_nrm_sq(Q);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm(Vec3dc<T> const& v)
{
    return std::sqrt(center_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm(BiVec3dc<T> const& B)
{
    return std::sqrt(center_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm(TriVec3dc<T> const& t)
{
    return std::sqrt(center_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm(QuadVec3dc<T> const& Q)
{
    return std::sqrt(center_nrm_sq(Q));
}

// scale a round object so its round weight norm becomes one (sign preserved;
// throws for flat objects)
// HINT: unitize() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Vec3dc<T> unitize(Vec3dc<T> const& v)
{
    T const wn = round_weight_nrm(v);
    hd::ga::detail::check_normalization<T>(wn, "round point (round weight)");
    T const inv = T(1.0) / wn;
    return Vec3dc<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv, v.u * inv);
}
template <typename T>
    requires(numeric_type<T>)
inline BiVec3dc<T> unitize(BiVec3dc<T> const& B)
{
    T const wn = round_weight_nrm(B);
    hd::ga::detail::check_normalization<T>(wn, "dipole (round weight)");
    T const inv = T(1.0) / wn;
    return BiVec3dc<T>(B.vx * inv, B.vy * inv, B.vz * inv, B.mx * inv, B.my * inv,
                       B.mz * inv, B.px * inv, B.py * inv, B.pz * inv, B.pw * inv);
}
template <typename T>
    requires(numeric_type<T>)
inline TriVec3dc<T> unitize(TriVec3dc<T> const& t)
{
    T const wn = round_weight_nrm(t);
    hd::ga::detail::check_normalization<T>(wn, "circle (round weight)");
    T const inv = T(1.0) / wn;
    return TriVec3dc<T>(t.vx * inv, t.vy * inv, t.vz * inv, t.mx * inv, t.my * inv,
                        t.mz * inv, t.px * inv, t.py * inv, t.pz * inv, t.pw * inv);
}
template <typename T>
    requires(numeric_type<T>)
inline QuadVec3dc<T> unitize(QuadVec3dc<T> const& Q)
{
    T const wn = round_weight_nrm(Q);
    hd::ga::detail::check_normalization<T>(wn, "sphere (round weight)");
    T const inv = T(1.0) / wn;
    return QuadVec3dc<T>(Q.x * inv, Q.y * inv, Q.z * inv, Q.w * inv, Q.u * inv);
}


} // namespace hd::ga::cga
