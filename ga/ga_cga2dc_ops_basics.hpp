#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_cga_types.hpp"  // cga2dc types incl. multivector wrappers
#include "detail/ga_foundation.hpp" // ga library headers and infrastructure

#include "detail/ga_error_handling.hpp" // unitize() zero-division checks


namespace hd::ga::cga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides cga2dc basic operations:
//
// - gr_inv()                      -> grade involution
// - rev()                         -> reversion
// - rrev()                        -> regressive reversion
// - conj()                        -> Clifford conjugation
// - cconj()                       -> conformal conjugate (round part kept,
//                                    flat part negated)
//
// - r_cmpl()                      -> right complement (non-metric)
// - l_cmpl()                      -> left complement (non-metric)
//
// - r_dual()                      -> right metric dual
// - l_dual()                      -> left metric dual
// - r_antidual()                  -> right metric antidual
// - l_antidual()                  -> left metric antidual
//
// - round_bulk(), round_weight(), flat_bulk(), flat_weight()
//                                 -> the four-part bulk/weight split
// - round_bulk_nrm{,_sq}(), round_weight_nrm{,_sq}(),
//   flat_bulk_nrm{,_sq}(), flat_weight_nrm{,_sq}() -> the four part norms
// - center_nrm{,_sq}()            -> weighted distance origin <-> center
// - unitize()                     -> scale so the round weight norm equals one
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
// conformal conjugate: keep the round part, negate the flat part (i.e. negate
// every component carrying a factor of the infinity dimension e4)
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
constexpr Vec2dc<T> cconj(Vec2dc<T> const& v)
{
    return Vec2dc<T>(v.x, v.y, v.z, -v.w);
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> cconj(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(B.vx, B.vy, B.vz, -B.mx, -B.my, -B.mz);
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> cconj(TriVec2dc<T> const& t)
{
    return TriVec2dc<T>(-t.x, -t.y, -t.z, t.w);
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

////////////////////////////////////////////////////////////////////////////////
// bulk/weight split
//
// Every component belongs to exactly one of four parts, by whether the blade
// contains the infinity dimension e4 (flat part: yes, round part: no) and
// the origin e3 (weight: yes, bulk: no):
//
//     round bulk:   no e4, no e3      round weight: no e4, with e3
//     flat bulk:    with e4, no e3    flat weight:  with e4 and e3
//
// The round weight carries the homogeneous weight of a round object, so
// unitizing makes the round weight norm one.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> round_bulk(Vec2dc<T> const& v)
{
    return Vec2dc<T>(v.x, v.y, T(0.0), T(0.0)); // e1, e2
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> round_bulk(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(T(0.0), T(0.0), B.vz, T(0.0), T(0.0), T(0.0)); // e12
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> round_bulk([[maybe_unused]] TriVec2dc<T> const& t)
{
    return TriVec2dc<T>(T(0.0), T(0.0), T(0.0), T(0.0)); // (none at grade 3)
}
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> round_weight(Vec2dc<T> const& v)
{
    return Vec2dc<T>(T(0.0), T(0.0), v.z, T(0.0)); // e3
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> round_weight(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(B.vx, B.vy, T(0.0), T(0.0), T(0.0), T(0.0)); // e31, e32
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> round_weight(TriVec2dc<T> const& t)
{
    return TriVec2dc<T>(T(0.0), T(0.0), T(0.0), t.w); // e321
}
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> flat_bulk(Vec2dc<T> const& v)
{
    return Vec2dc<T>(T(0.0), T(0.0), T(0.0), v.w); // e4
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> flat_bulk(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(T(0.0), T(0.0), T(0.0), B.mx, B.my, T(0.0)); // e14, e24
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> flat_bulk(TriVec2dc<T> const& t)
{
    return TriVec2dc<T>(T(0.0), T(0.0), t.z, T(0.0)); // e124
}
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> flat_weight([[maybe_unused]] Vec2dc<T> const& v)
{
    return Vec2dc<T>(T(0.0), T(0.0), T(0.0), T(0.0)); // (no e34-type component)
}
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> flat_weight(BiVec2dc<T> const& B)
{
    return BiVec2dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.mz); // e34
}
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> flat_weight(TriVec2dc<T> const& t)
{
    return TriVec2dc<T>(t.x, t.y, T(0.0), T(0.0)); // e314, e324
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
constexpr T round_bulk_nrm_sq(Vec2dc<T> const& v)
{
    return v.x * v.x + v.y * v.y;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm_sq(BiVec2dc<T> const& B)
{
    return B.vz * B.vz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm_sq([[maybe_unused]] TriVec2dc<T> const& t)
{
    return T(0.0);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm_sq(Vec2dc<T> const& v)
{
    return v.z * v.z;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm_sq(BiVec2dc<T> const& B)
{
    return B.vx * B.vx + B.vy * B.vy;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm_sq(TriVec2dc<T> const& t)
{
    return t.w * t.w;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm_sq(Vec2dc<T> const& v)
{
    return v.w * v.w;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm_sq(BiVec2dc<T> const& B)
{
    return B.mx * B.mx + B.my * B.my;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm_sq(TriVec2dc<T> const& t)
{
    return t.z * t.z;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm_sq([[maybe_unused]] Vec2dc<T> const& v)
{
    return T(0.0);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm_sq(BiVec2dc<T> const& B)
{
    return B.mz * B.mz;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm_sq(TriVec2dc<T> const& t)
{
    return t.x * t.x + t.y * t.y;
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm(Vec2dc<T> const& v)
{
    return std::sqrt(round_bulk_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm(BiVec2dc<T> const& B)
{
    return std::sqrt(round_bulk_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_bulk_nrm(TriVec2dc<T> const& t)
{
    return std::sqrt(round_bulk_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm(Vec2dc<T> const& v)
{
    return std::sqrt(round_weight_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm(BiVec2dc<T> const& B)
{
    return std::sqrt(round_weight_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T round_weight_nrm(TriVec2dc<T> const& t)
{
    return std::sqrt(round_weight_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm(Vec2dc<T> const& v)
{
    return std::sqrt(flat_bulk_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm(BiVec2dc<T> const& B)
{
    return std::sqrt(flat_bulk_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_bulk_nrm(TriVec2dc<T> const& t)
{
    return std::sqrt(flat_bulk_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm(Vec2dc<T> const& v)
{
    return std::sqrt(flat_weight_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm(BiVec2dc<T> const& B)
{
    return std::sqrt(flat_weight_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T flat_weight_nrm(TriVec2dc<T> const& t)
{
    return std::sqrt(flat_weight_nrm_sq(t));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm_sq(Vec2dc<T> const& v)
{
    return round_bulk_nrm_sq(v) + flat_weight_nrm_sq(v);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm_sq(BiVec2dc<T> const& B)
{
    return round_bulk_nrm_sq(B) + flat_weight_nrm_sq(B);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm_sq(TriVec2dc<T> const& t)
{
    return round_bulk_nrm_sq(t) + flat_weight_nrm_sq(t);
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm(Vec2dc<T> const& v)
{
    return std::sqrt(center_nrm_sq(v));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm(BiVec2dc<T> const& B)
{
    return std::sqrt(center_nrm_sq(B));
}
template <typename T>
    requires(numeric_type<T>)
constexpr T center_nrm(TriVec2dc<T> const& t)
{
    return std::sqrt(center_nrm_sq(t));
}

// scale a round object so its round weight norm becomes one (sign preserved;
// throws for flat objects)
// HINT: unitize() cannot be constexpr due to the checks for division by zero
//       which might throw

// scale a round object so its round weight norm becomes one (sign preserved;
// throws for flat objects)
template <typename T>
    requires(numeric_type<T>)
inline Vec2dc<T> unitize(Vec2dc<T> const& v)
{
    T const wn = round_weight_nrm(v);
    hd::ga::detail::check_normalization<T>(wn, "round point (round weight)");
    T const inv = T(1.0) / wn;
    return Vec2dc<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}
template <typename T>
    requires(numeric_type<T>)
inline BiVec2dc<T> unitize(BiVec2dc<T> const& B)
{
    T const wn = round_weight_nrm(B);
    hd::ga::detail::check_normalization<T>(wn, "dipole (round weight)");
    T const inv = T(1.0) / wn;
    return BiVec2dc<T>(B.vx * inv, B.vy * inv, B.vz * inv, B.mx * inv, B.my * inv,
                       B.mz * inv);
}
template <typename T>
    requires(numeric_type<T>)
inline TriVec2dc<T> unitize(TriVec2dc<T> const& t)
{
    T const wn = round_weight_nrm(t);
    hd::ga::detail::check_normalization<T>(wn, "circle (round weight)");
    T const inv = T(1.0) / wn;
    return TriVec2dc<T>(t.x * inv, t.y * inv, t.z * inv, t.w * inv);
}


} // namespace hd::ga::cga
