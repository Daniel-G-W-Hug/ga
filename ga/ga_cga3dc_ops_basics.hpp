#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_cga_types.hpp"  // cga3dc types incl. multivector wrappers
#include "detail/ga_foundation.hpp" // ga library headers and infrastructure


namespace hd::ga::cga {

////////////////////////////////////////////////////////////////////////////////
// provides cga3dc basic operations:
//
// - gr_inv()                      -> grade involution
// - rev()                         -> reversion
// - rrev()                        -> regressive reversion (== rev in 5d)
// - conj()                        -> Clifford conjugation
//
// - cmpl()                        -> complement (non-metric; left == right in 5d)
//
// - dual()                        -> metric dual
// - antidual()                    -> metric antidual
//
// bulk/weight, the norms and normalize are added together with the geometric
// operations layer (the round objects extend the norm story vs. the flat
// algebras)
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


} // namespace hd::ga::cga
