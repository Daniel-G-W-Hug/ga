#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_cga_types.hpp"  // cga2dc types incl. multivector wrappers
#include "detail/ga_foundation.hpp" // ga library headers and infrastructure


namespace hd::ga::cga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides cga2dc basic operations:
//
// - gr_inv()                      -> grade involution
// - rev()                         -> reversion
// - rrev()                        -> regressive reversion
// - conj()                        -> Clifford conjugation
//
// - r_cmpl()                      -> right complement (non-metric)
// - l_cmpl()                      -> left complement (non-metric)
//
// - r_dual()                      -> right metric dual
// - l_dual()                      -> left metric dual
// - r_antidual()                  -> right metric antidual
// - l_antidual()                  -> left metric antidual
//
// bulk/weight, the norms and normalize are added together with the geometric
// operations layer (the round objects extend the norm story vs. the flat
// algebras)
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for r = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> gr_inv(Scalar2dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> gr_inv(Vec2dc<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> gr_inv(BiVec2dc<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> gr_inv(TriVec2dc<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> gr_inv(PScalar2dc<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> gr_inv(MVec2dc_E<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    return M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> gr_inv(MVec2dc_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> gr_inv(MVec2dc<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    // grade 1 and 3: sign reversal
    return MVec2dc<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)),
                      gr_inv(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> rev(Scalar2dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> rev(Vec2dc<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> rev(BiVec2dc<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> rev(TriVec2dc<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> rev(PScalar2dc<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> rev(MVec2dc_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec2dc_E<T>(rev(gr0(M)), rev(gr2(M)), rev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> rev(MVec2dc_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec2dc_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> rev(MVec2dc<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    // grade 4: no sign change
    return MVec2dc<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)), rev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// regressive reversion operation: reverse of complement + backtrafo
// rrev(A_r) = l_cmpl((-1)^(r*(r-1)/2) r_cmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=4, r = 0, 1, 2, 3, ...: + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> rrev(Scalar2dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> rrev(Vec2dc<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> rrev(BiVec2dc<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> rrev(TriVec2dc<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> rrev(PScalar2dc<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> rrev(MVec2dc_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec2dc_E<T>(rrev(gr0(M)), rrev(gr2(M)), rrev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> rrev(MVec2dc_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dc_U<T>(rrev(gr1(M)), rrev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> rrev(MVec2dc<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec2dc<T>(rrev(gr0(M)), rrev(gr1(M)), rrev(gr2(M)), rrev(gr3(M)),
                      rrev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation (identical sign pattern to regressive reverse for cga2dc):
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> conj(Scalar2dc<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> conj(Vec2dc<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> conj(BiVec2dc<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> conj(TriVec2dc<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> conj(PScalar2dc<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> conj(MVec2dc_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec2dc_E<T>(conj(gr0(M)), conj(gr2(M)), conj(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> conj(MVec2dc_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec2dc_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> conj(MVec2dc<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec2dc<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)),
                      conj(gr4(M)));
}
////////////////////////////////////////////////////////////////////////////////
// complement operations (non-metric)
//
// The complement fills a blade up to the pseudoscalar of the space:
//
//     l_cmpl(u) ^ u = I_2dc        u ^ r_cmpl(u) = I_2dc
//
// For the even-dimensional cga2dc the left and right complement differ at the
// odd grades (as for ega2d, pga3dp, sta4ds).
////////////////////////////////////////////////////////////////////////////////

// cga2dc l_cmpl :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> l_cmpl(Scalar2dc<T> s)
{
    return PScalar2dc<T>(T(s));
}

// cga2dc l_cmpl :: vec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> l_cmpl(Vec2dc<T> const& v)
{
    return TriVec2dc<T>(-v.y, v.x, -v.z, -v.w);
}

// cga2dc l_cmpl :: bivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> l_cmpl(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(B.my, -B.mx, B.mz, -B.vy, B.vx, B.vz);
}

// cga2dc l_cmpl :: trivec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> l_cmpl(TriVec2dc<T> const& t)
{
    return Vec2dc<T>(-t.y, t.x, t.z, t.w);
}

// cga2dc l_cmpl :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> l_cmpl(PScalar2dc<T> ps)
{
    return Scalar2dc<T>(T(ps));
}

// cga2dc l_cmpl :: mv_e -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> l_cmpl(MVec2dc_E<T> const& M)
{
    return MVec2dc_E<T>(l_cmpl(gr4(M)), l_cmpl(gr2(M)), l_cmpl(gr0(M)));
}

// cga2dc l_cmpl :: mv_u -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> l_cmpl(MVec2dc_U<T> const& M)
{
    return MVec2dc_U<T>(l_cmpl(gr3(M)), l_cmpl(gr1(M)));
}

// cga2dc l_cmpl :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> l_cmpl(MVec2dc<T> const& M)
{
    return MVec2dc<T>(l_cmpl(gr4(M)), l_cmpl(gr3(M)), l_cmpl(gr2(M)), l_cmpl(gr1(M)),
                      l_cmpl(gr0(M)));
}

// cga2dc r_cmpl :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> r_cmpl(Scalar2dc<T> s)
{
    return PScalar2dc<T>(T(s));
}

// cga2dc r_cmpl :: vec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> r_cmpl(Vec2dc<T> const& v)
{
    return TriVec2dc<T>(v.y, -v.x, v.z, v.w);
}

// cga2dc r_cmpl :: bivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> r_cmpl(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(B.my, -B.mx, B.mz, -B.vy, B.vx, B.vz);
}

// cga2dc r_cmpl :: trivec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> r_cmpl(TriVec2dc<T> const& t)
{
    return Vec2dc<T>(t.y, -t.x, -t.z, -t.w);
}

// cga2dc r_cmpl :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> r_cmpl(PScalar2dc<T> ps)
{
    return Scalar2dc<T>(T(ps));
}

// cga2dc r_cmpl :: mv_e -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> r_cmpl(MVec2dc_E<T> const& M)
{
    return MVec2dc_E<T>(r_cmpl(gr4(M)), r_cmpl(gr2(M)), r_cmpl(gr0(M)));
}

// cga2dc r_cmpl :: mv_u -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> r_cmpl(MVec2dc_U<T> const& M)
{
    return MVec2dc_U<T>(r_cmpl(gr3(M)), r_cmpl(gr1(M)));
}

// cga2dc r_cmpl :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> r_cmpl(MVec2dc<T> const& M)
{
    return MVec2dc<T>(r_cmpl(gr4(M)), r_cmpl(gr3(M)), r_cmpl(gr2(M)), r_cmpl(gr1(M)),
                      r_cmpl(gr0(M)));
}
////////////////////////////////////////////////////////////////////////////////
// duals (= complement of the metric image: dual(a) = cmpl(G a))
//
// The cga2dc metric is non-degenerate, so there is no bulk/weight dual split
// (as needed for the degenerate pga metrics): the extended metric G is a
// signed permutation that pairs every projective (w-bearing) blade with a
// round (u-bearing) partner blade, and the duals stay single-term:
//
//     l_dual(a) = l_cmpl(G a)        r_dual(a) = r_cmpl(G a)
//
////////////////////////////////////////////////////////////////////////////////

// cga2dc l_dual :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> l_dual(Scalar2dc<T> s)
{
    return PScalar2dc<T>(T(s));
}

// cga2dc l_dual :: vec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> l_dual(Vec2dc<T> const& v)
{
    return TriVec2dc<T>(-v.y, v.x, v.w, v.z);
}

// cga2dc l_dual :: bivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> l_dual(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(B.vy, -B.vx, -B.mz, -B.my, B.mx, B.vz);
}

// cga2dc l_dual :: trivec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> l_dual(TriVec2dc<T> const& t)
{
    return Vec2dc<T>(t.y, -t.x, t.w, t.z);
}

// cga2dc l_dual :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> l_dual(PScalar2dc<T> ps)
{
    return Scalar2dc<T>(-T(ps));
}

// cga2dc l_dual :: mv_e -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> l_dual(MVec2dc_E<T> const& M)
{
    return MVec2dc_E<T>(l_dual(gr4(M)), l_dual(gr2(M)), l_dual(gr0(M)));
}

// cga2dc l_dual :: mv_u -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> l_dual(MVec2dc_U<T> const& M)
{
    return MVec2dc_U<T>(l_dual(gr3(M)), l_dual(gr1(M)));
}

// cga2dc l_dual :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> l_dual(MVec2dc<T> const& M)
{
    return MVec2dc<T>(l_dual(gr4(M)), l_dual(gr3(M)), l_dual(gr2(M)), l_dual(gr1(M)),
                      l_dual(gr0(M)));
}

// cga2dc r_dual :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> r_dual(Scalar2dc<T> s)
{
    return PScalar2dc<T>(T(s));
}

// cga2dc r_dual :: vec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> r_dual(Vec2dc<T> const& v)
{
    return TriVec2dc<T>(v.y, -v.x, -v.w, -v.z);
}

// cga2dc r_dual :: bivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> r_dual(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(B.vy, -B.vx, -B.mz, -B.my, B.mx, B.vz);
}

// cga2dc r_dual :: trivec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> r_dual(TriVec2dc<T> const& t)
{
    return Vec2dc<T>(-t.y, t.x, -t.w, -t.z);
}

// cga2dc r_dual :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> r_dual(PScalar2dc<T> ps)
{
    return Scalar2dc<T>(-T(ps));
}

// cga2dc r_dual :: mv_e -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> r_dual(MVec2dc_E<T> const& M)
{
    return MVec2dc_E<T>(r_dual(gr4(M)), r_dual(gr2(M)), r_dual(gr0(M)));
}

// cga2dc r_dual :: mv_u -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> r_dual(MVec2dc_U<T> const& M)
{
    return MVec2dc_U<T>(r_dual(gr3(M)), r_dual(gr1(M)));
}

// cga2dc r_dual :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> r_dual(MVec2dc<T> const& M)
{
    return MVec2dc<T>(r_dual(gr4(M)), r_dual(gr3(M)), r_dual(gr2(M)), r_dual(gr1(M)),
                      r_dual(gr0(M)));
}
////////////////////////////////////////////////////////////////////////////////
// antiduals (= complement of the anti-metric image: antidual(a) = cmpl(AG a))
//
// The antidual uses the metric anti-exomorphism AG (the extended metric built
// w.r.t. the antiwedge product). For cga2dc the anti-exomorphism equals -G
// (verified at rule-generation time), so the antidual is the negated dual:
//
//     l_antidual(a) = -l_dual(a)     r_antidual(a) = -r_dual(a)
//
////////////////////////////////////////////////////////////////////////////////

// cga2dc l_antidual :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> l_antidual(Scalar2dc<T> s)
{
    return PScalar2dc<T>(-T(s));
}

// cga2dc l_antidual :: vec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> l_antidual(Vec2dc<T> const& v)
{
    return TriVec2dc<T>(v.y, -v.x, -v.w, -v.z);
}

// cga2dc l_antidual :: bivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> l_antidual(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(-B.vy, B.vx, B.mz, B.my, -B.mx, -B.vz);
}

// cga2dc l_antidual :: trivec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> l_antidual(TriVec2dc<T> const& t)
{
    return Vec2dc<T>(-t.y, t.x, -t.w, -t.z);
}

// cga2dc l_antidual :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> l_antidual(PScalar2dc<T> ps)
{
    return Scalar2dc<T>(T(ps));
}

// cga2dc l_antidual :: mv_e -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> l_antidual(MVec2dc_E<T> const& M)
{
    return MVec2dc_E<T>(l_antidual(gr4(M)), l_antidual(gr2(M)), l_antidual(gr0(M)));
}

// cga2dc l_antidual :: mv_u -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> l_antidual(MVec2dc_U<T> const& M)
{
    return MVec2dc_U<T>(l_antidual(gr3(M)), l_antidual(gr1(M)));
}

// cga2dc l_antidual :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> l_antidual(MVec2dc<T> const& M)
{
    return MVec2dc<T>(l_antidual(gr4(M)), l_antidual(gr3(M)), l_antidual(gr2(M)),
                      l_antidual(gr1(M)), l_antidual(gr0(M)));
}

// cga2dc r_antidual :: s -> ps
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> r_antidual(Scalar2dc<T> s)
{
    return PScalar2dc<T>(-T(s));
}

// cga2dc r_antidual :: vec -> trivec
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> r_antidual(Vec2dc<T> const& v)
{
    return TriVec2dc<T>(-v.y, v.x, v.w, v.z);
}

// cga2dc r_antidual :: bivec -> bivec
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> r_antidual(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(-B.vy, B.vx, B.mz, B.my, -B.mx, -B.vz);
}

// cga2dc r_antidual :: trivec -> vec
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> r_antidual(TriVec2dc<T> const& t)
{
    return Vec2dc<T>(t.y, -t.x, t.w, t.z);
}

// cga2dc r_antidual :: ps -> s
template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> r_antidual(PScalar2dc<T> ps)
{
    return Scalar2dc<T>(T(ps));
}

// cga2dc r_antidual :: mv_e -> mv_e
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> r_antidual(MVec2dc_E<T> const& M)
{
    return MVec2dc_E<T>(r_antidual(gr4(M)), r_antidual(gr2(M)), r_antidual(gr0(M)));
}

// cga2dc r_antidual :: mv_u -> mv_u
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_U<T> r_antidual(MVec2dc_U<T> const& M)
{
    return MVec2dc_U<T>(r_antidual(gr3(M)), r_antidual(gr1(M)));
}

// cga2dc r_antidual :: mv -> mv
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc<T> r_antidual(MVec2dc<T> const& M)
{
    return MVec2dc<T>(r_antidual(gr4(M)), r_antidual(gr3(M)), r_antidual(gr2(M)),
                      r_antidual(gr1(M)), r_antidual(gr0(M)));
}
} // namespace hd::ga::cga
