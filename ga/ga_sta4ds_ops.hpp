#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_sta4ds_ops_basics.hpp"
#include "ga_sta4ds_ops_products.hpp"

#include <algorithm> // std::clamp, std::max
#include <array>     // std::array (transform_opt coefficient matrices)
#include <cmath>     // std::cos, std::sin, std::cosh, std::sinh, std::sqrt, std::abs
#include <stdexcept> // std::runtime_error


namespace hd::ga::detail {

// Closed-form rotor-sandwich coefficient matrices for sta::transform_opt(), kept in
// `detail` so the binding generator (which scans hd::ga / ::ega / ::pga / ::sta) skips
// them. Each is built once from the rotor coefficients R.c0..R.c7 and reused by both
// the scalar and the std::vector batch transform_opt overloads.
//
// The matrices are the "after transformation" coefficients emitted by ga_prdxpr, with
// every distinct product R.ci*R.cj collected into a named temporary once (squares reused
// on the diagonal, mixed products reused across symmetric off-diagonal entries).

// 4x4 matrix M(R) for the grade-1 (vector) AND grade-3 (trivector) sandwich -- these
// share identical coefficients in sta. Row-major: entry (row r, col c) is k[4*r + c],
// rows/cols ordered (g1,g2,g3,g4) for vectors / (g234,g314,g124,g123) for trivectors.
template <typename T> inline std::array<T, 16> sta_rotor_xf_mat_vec(MVec4ds_E<T> const& R)
{
    T const a0 = R.c0 * R.c0, a1 = R.c1 * R.c1, a2 = R.c2 * R.c2, a3 = R.c3 * R.c3;
    T const a4 = R.c4 * R.c4, a5 = R.c5 * R.c5, a6 = R.c6 * R.c6, a7 = R.c7 * R.c7;

    T const b01 = R.c0 * R.c1, b02 = R.c0 * R.c2, b03 = R.c0 * R.c3, b04 = R.c0 * R.c4;
    T const b05 = R.c0 * R.c5, b06 = R.c0 * R.c6;
    T const b12 = R.c1 * R.c2, b13 = R.c1 * R.c3, b15 = R.c1 * R.c5, b16 = R.c1 * R.c6;
    T const b17 = R.c1 * R.c7;
    T const b23 = R.c2 * R.c3, b24 = R.c2 * R.c4, b26 = R.c2 * R.c6, b27 = R.c2 * R.c7;
    T const b34 = R.c3 * R.c4, b35 = R.c3 * R.c5, b37 = R.c3 * R.c7;
    T const b45 = R.c4 * R.c5, b46 = R.c4 * R.c6, b47 = R.c4 * R.c7;
    T const b56 = R.c5 * R.c6, b57 = R.c5 * R.c7;
    T const b67 = R.c6 * R.c7;

    return {a0 + a1 - a2 - a3 + a4 - a5 - a6 + a7, T(2.0) * (-b06 + b12 + b37 + b45),
            T(2.0) * (b05 + b13 - b27 + b46),      T(2.0) * (b01 - b26 + b35 + b47),

            T(2.0) * (b06 + b12 - b37 + b45),      a0 - a1 + a2 - a3 - a4 + a5 - a6 + a7,
            T(2.0) * (-b04 + b17 + b23 + b56),     T(2.0) * (b02 + b16 - b34 + b57),

            T(2.0) * (-b05 + b13 + b27 + b46),     T(2.0) * (b04 - b17 + b23 + b56),
            a0 - a1 - a2 + a3 - a4 - a5 + a6 + a7, T(2.0) * (b03 - b15 + b24 + b67),

            T(2.0) * (b01 + b26 - b35 + b47),      T(2.0) * (b02 - b16 + b34 + b57),
            T(2.0) * (b03 + b15 - b24 + b67),      a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7};
}

// 6x6 matrix M(R) for the grade-2 (bivector) sandwich. Row-major: entry (row r, col c)
// is k[6*r + c]; rows/cols ordered (vx,vy,vz,mx,my,mz) = (g14,g24,g34,g23,g31,g12).
// All 28 distinct mixed products appear here (none is unused).
template <typename T>
inline std::array<T, 36> sta_rotor_xf_mat_bivec(MVec4ds_E<T> const& R)
{
    T const a0 = R.c0 * R.c0, a1 = R.c1 * R.c1, a2 = R.c2 * R.c2, a3 = R.c3 * R.c3;
    T const a4 = R.c4 * R.c4, a5 = R.c5 * R.c5, a6 = R.c6 * R.c6, a7 = R.c7 * R.c7;

    T const b01 = R.c0 * R.c1, b02 = R.c0 * R.c2, b03 = R.c0 * R.c3, b04 = R.c0 * R.c4;
    T const b05 = R.c0 * R.c5, b06 = R.c0 * R.c6, b07 = R.c0 * R.c7;
    T const b12 = R.c1 * R.c2, b13 = R.c1 * R.c3, b14 = R.c1 * R.c4, b15 = R.c1 * R.c5;
    T const b16 = R.c1 * R.c6, b17 = R.c1 * R.c7;
    T const b23 = R.c2 * R.c3, b24 = R.c2 * R.c4, b25 = R.c2 * R.c5, b26 = R.c2 * R.c6;
    T const b27 = R.c2 * R.c7;
    T const b34 = R.c3 * R.c4, b35 = R.c3 * R.c5, b36 = R.c3 * R.c6, b37 = R.c3 * R.c7;
    T const b45 = R.c4 * R.c5, b46 = R.c4 * R.c6, b47 = R.c4 * R.c7;
    T const b56 = R.c5 * R.c6, b57 = R.c5 * R.c7;
    T const b67 = R.c6 * R.c7;

    return {// g14 <- (vx,vy,vz,mx,my,mz)
            a0 - a1 + a2 + a3 + a4 - a5 - a6 - a7, T(2.0) * (-b06 - b12 - b37 + b45),
            T(2.0) * (b05 - b13 + b27 + b46), T(2.0) * (-b07 + b14 - b25 - b36),
            T(2.0) * (-b03 + b15 + b24 + b67), T(2.0) * (b02 + b16 + b34 - b57),
            // g24
            T(2.0) * (b06 - b12 + b37 + b45), a0 + a1 - a2 + a3 - a4 + a5 - a6 - a7,
            T(2.0) * (-b04 - b17 - b23 + b56), T(2.0) * (b03 + b15 + b24 - b67),
            T(2.0) * (-b07 - b14 + b25 - b36), T(2.0) * (-b01 + b26 + b35 + b47),
            // g34
            T(2.0) * (-b05 - b13 - b27 + b46), T(2.0) * (b04 + b17 - b23 + b56),
            a0 + a1 + a2 - a3 - a4 - a5 + a6 - a7, T(2.0) * (-b02 + b16 + b34 + b57),
            T(2.0) * (b01 + b26 + b35 - b47), T(2.0) * (-b07 - b14 - b25 + b36),
            // g23
            T(2.0) * (b07 - b14 + b25 + b36), T(2.0) * (b03 - b15 - b24 - b67),
            T(2.0) * (-b02 - b16 - b34 + b57), a0 - a1 + a2 + a3 + a4 - a5 - a6 - a7,
            T(2.0) * (-b06 - b12 - b37 + b45), T(2.0) * (b05 - b13 + b27 + b46),
            // g31
            T(2.0) * (-b03 - b15 - b24 + b67), T(2.0) * (b07 + b14 - b25 + b36),
            T(2.0) * (b01 - b26 - b35 - b47), T(2.0) * (b06 - b12 + b37 + b45),
            a0 + a1 - a2 + a3 - a4 + a5 - a6 - a7, T(2.0) * (-b04 - b17 - b23 + b56),
            // g12
            T(2.0) * (b02 - b16 - b34 - b57), T(2.0) * (-b01 - b26 - b35 + b47),
            T(2.0) * (b07 + b14 + b25 - b36), T(2.0) * (-b05 - b13 - b27 + b46),
            T(2.0) * (b04 + b17 - b23 + b56), a0 + a1 + a2 - a3 - a4 - a5 + a6 - a7};
}

/////////////////////////////////////////////////////////////////////////////////////////
// rotor exp / log / sqrt support helpers.
//
// The signature-agnostic renormalisation and invariant bivector decomposition follow
// De Keninck & Roelfs, "Normalization, Square Roots, and the Exponential and Logarithmic
// Maps in Geometric Algebras of Less than 6D" (2022, arXiv:2206.07496). STA is R(1,3)
// used here on the DUAL of the paper's R(3,1) PGA convention -- the geometric-product
// rotor sandwich is the same, and the pseudoscalar squares to I^2 = -1 (a complex-type
// Study number), so the paper's Listing 1 signs must be re-derived for this signature;
// the code below is written signature-agnostically from the paper's sec. 5-6 instead.
//
// These live in hd::ga::detail (like the transform_opt matrices above) so the binding
// generator -- which scans hd::ga / ::ega / ::pga / ::sta -- skips them. They use the
// sta:: geometric product and reversion, brought in with local using-declarations since
// ADL on the hd::ga-scoped operand types would not otherwise reach hd::ga::sta.
/////////////////////////////////////////////////////////////////////////////////////////

// exp of a SIMPLE bivector b whose geometric square is the scalar lambda = gr0(b*b):
//   lambda < 0 (spacelike / rotation):  cos(t)  + (sin(t)/t)  b,  t = sqrt(-lambda)
//   lambda > 0 (timelike  / boost):     cosh(t) + (sinh(t)/t) b,  t = sqrt( lambda)
//   lambda = 0 (lightlike / null):      1 + b   (exact: b*b == 0 kills all higher terms)
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> sta4ds_exp_simple(BiVec4ds<T> const& b, T lambda)
{
    if (lambda < -safe_epsilon<T>()) {
        T const t = std::sqrt(-lambda);
        return MVec4ds_E<T>(Scalar4ds<T>(std::cos(t)), (std::sin(t) / t) * b);
    }
    if (lambda > safe_epsilon<T>()) {
        T const t = std::sqrt(lambda);
        return MVec4ds_E<T>(Scalar4ds<T>(std::cosh(t)), (std::sinh(t) / t) * b);
    }
    return MVec4ds_E<T>(Scalar4ds<T>(1.0), b); // null / zero plane -> 1 + b (exact)
}

// I * B for a bivector B (geometric product by the unit pseudoscalar) -> bivector.
// Maps a rotation 2-plane to its orthogonal (dual) boost 2-plane and vice versa;
// (I b_hat)^2 = I^2 b_hat^2 flips the causal character (unit rotation <-> unit boost).
template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> sta4ds_biv_dual(BiVec4ds<T> const& B)
{
    using sta::operator*;
    return gr2(PScalar4ds<T>(1.0) * MVec4ds_E<T>(B));
}

// Invariant (Study-projector) decomposition of an STA bivector into two orthogonal,
// commuting, SIMPLE bivectors B = b_boost + b_rot (De Keninck & Roelfs sec. 6). Both
// the geometric square and the wedge square are proper Lorentz invariants, so the
// decomposition is observer-INDEPENDENT -- unlike the spacetime split rel_vec_split() /
// rel_bivec_split(), which is taken relative to a chosen time direction u.
//
// With I^2 = -1 the two squares straddle zero, lambda_boost >= 0 >= lambda_rot, so
// b_boost generates a Lorentz boost and b_rot a spatial rotation.
//
// A SIMPLE B (B ^ B == 0) already lies in a single plane; it is then returned in the
// slot matching its own causal character, with the other slot zero. This keeps the
// invariant "b_boost never squares negative, b_rot never squares positive" in all cases.
//
// Shared by exp() and by the public boost_part() / rot_part().
template <typename T> struct sta4ds_biv_parts {
    BiVec4ds<T> b_boost; // time-like part,  b_boost^2 = lambda_boost >= 0
    BiVec4ds<T> b_rot;   // space-like part, b_rot^2   = lambda_rot   <= 0
    T lambda_boost;
    T lambda_rot;
    T geom_sq;   // <B^2>_0, the geometric square of the input
    bool simple; // true if B ^ B == 0
};

template <typename T>
    requires(numeric_type<T>)
inline sta4ds_biv_parts<T> sta4ds_biv_decompose(BiVec4ds<T> const& B)
{
    using sta::operator*;
    BiVec4ds<T> const zero(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0));

    T const bb_s = sta4ds_geom_sq(B);            // <B^2>_0 (geometric square)
    T const bb_ps = T(gr4(MVec4ds_E<T>(B) * B)); // <B^2>_4 (wedge B^B)

    if (std::abs(bb_ps) <= safe_epsilon<T>()) {
        // SIMPLE: a single plane -- sort it by its own causal character
        if (bb_s < T(0.0)) return {zero, B, T(0.0), bb_s, bb_s, true};
        return {B, zero, bb_s, T(0.0), bb_s, true};
    }

    // NON-simple: Study projectors P_+/- = (1 +/- conj(B^2)/||B^2||_S)/2,
    // with conj(B^2) = bb_s - bb_ps I
    T const nS = std::sqrt(bb_s * bb_s + bb_ps * bb_ps); // ||B^2||_S
    BiVec4ds<T> const proj =
        (bb_s * B - bb_ps * sta4ds_biv_dual(B)) / nS; // conj(B^2) B / ||B^2||_S
    return {T(0.5) * (B + proj),
            T(0.5) * (B - proj),
            T(0.5) * (bb_s + nS),
            T(0.5) * (bb_s - nS),
            bb_s,
            false};
}

// Nearest rotor to an even element X (De Keninck & Roelfs eq. 24): the signature-agnostic
// renormalisation via the self-reverse Study number  X*rev(X) = s + t*I  (I^2 = -1),
//   R = (X rev(X))^{-1/2} X = alpha X + beta (I X),
// which stays a proper unit rotor even for NON-simple X (t != 0) -- the case a
// scalar-only norm |<X rev(X)>_0|^{-1/2} would silently drop. Reduces to X / sqrt(s) when
// t == 0.
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> sta4ds_nearest_rotor(MVec4ds_E<T> const& X)
{
    using sta::operator*;
    using sta::rev;
    MVec4ds_E<T> const XX = X * rev(X);
    T const s = T(gr0(XX));                // <X rev(X)>_0
    T const t = T(gr4(XX));                // <X rev(X)>_4  (coeff of g1234)
    T const nS = std::sqrt(s * s + t * t); // ||X rev(X)||_S
    if (nS <= safe_epsilon<T>()) {
        // singular (X ~ 0, e.g. X = 1 + R with R = -1): no unique nearest rotor
        return MVec4ds_E<T>(Scalar4ds<T>(1.0), BiVec4ds<T>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                            PScalar4ds<T>(0.0));
    }
    T const c = std::sqrt(T(0.5) * (s + nS));       // c+ branch (natural for I^2 = -1)
    T const denom = T(4.0) * c * c * c * c + t * t; // 4 c^4 - <XX~>_4^2 = 4 c^4 + t^2
    T const alpha = T(4.0) * c * c * c / denom;
    T const beta = -T(2.0) * c * t / denom;
    return alpha * X + beta * (PScalar4ds<T>(1.0) * X); // alpha X + beta (I X)
}

} // namespace hd::ga::detail


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// High-level STA (G(1,3,0)) operations built on the geometric product (gpr).
//
// Unlike pga3dp (which uses regressive motors for rigid motions), STA transforms via
// the geometric-product rotor sandwich  X' = R * X * rev(R),  with R = exp(bivector).
// The structure mirrors ega3d_ops.hpp (also gpr-based), extended for the Lorentzian
// signature: spatial bivectors  (g23, g31, g12; B^2 < 0) generate rotations,
//            timelike bivectors (g14, g24, g34; B^2 > 0) generate Lorentz boosts.
//
// Implemented:
//   - exp(BiVec)                      -> rotor exponential of a bivector (simple or not)
//   - exp(PScalar)                    -> cos + I sin: a duality-rotation factor for
//                                        bivector fields, NOT a rotor (see the block)
//   - log(rotor)                      -> bivector log of a rotor (inverse of exp)
//   - get_rotor(plane, angle)         -> rotor for a spatial rotation
//   - get_boost(plane, phi)           -> rotor for a Lorentz boost (rapidity phi)
//   - sqrt(rotor)                     -> rotor halving the rotation angle / rapidity
//   - angle() / rapidity()            -> separation of two vectors (spacelike / timelike)
//   - transform(X, R)                 -> apply a rotor via the sandwich R*X*rev(R)
//   - transform_opt(X, R)             -> closed-form transform, vec/bivec/trivec
//                                        (scalar + std::vector batch overloads)
//   - time_split() / space_split()    -> spacetime split of a vector (time + rel. space)
//   - rel_vec_split() / rel_bivec_split() -> spacetime split of a bivector (E / B parts)
//   - is_simple()                     -> does the bivector lie in a single plane?
//   - boost_part() / rot_part()       -> invariant (observer-independent) decomposition
//                                        of a bivector into its two orthogonal planes
//   - project_onto() / reject_from()  -> projection / rejection (onto vector or bivector)
//   - reflect_on() / reflect_on_vec() -> reflections (hyperplane, 2-plane, vector)
//
//   - is_congruent()                  -> same subspace up to a scalar factor
//   - is_close()                      -> same value within a RELATIVE tolerance
//   - is_same_transform()             -> do two rotors describe the same Lorentz
//                                        transformation? (rotors double-cover them)
//
// exp() / log() / sqrt() handle general (non-simple) rotors -- a Lorentz boost and a
// spatial rotation combined in dual planes -- via the invariant bivector decomposition
// and the signature-agnostic Study-number renormalisation (De Keninck & Roelfs 2022,
// arXiv:2206.07496); see the detail:: helpers at the top of this file.
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// rotor exponential w.r.t. the geometric product
//
// B*B is a self-reverse Study number  B*B = <B^2>_0 + <B^2>_4 = bb_s + bb_ps I:
//   bb_s  = gr0(B*B)  the geometric square (its sign = the plane's causal character)
//   bb_ps = gr4(B*B)  the wedge B^B  (zero iff B is SIMPLE, i.e. a single 2-plane)
//
// SIMPLE B (bb_ps == 0), with a = sqrt(|bb_s|):
//   B spacelike (bb_s < 0):  exp(B) = cos(a)  + (B/a) sin(a)    (rotation)
//   B timelike  (bb_s > 0):  exp(B) = cosh(a) + (B/a) sinh(a)   (Lorentz boost)
//   B lightlike (bb_s == 0): exp(B) = 1 + B   (exact -- B*B == 0 kills higher terms)
//
// NON-simple B (bb_ps != 0): the invariant decomposition splits B into two orthogonal
// commuting simple bivectors  B = b_boost + b_rot  (De Keninck & Roelfs sec. 6). With
// I^2 = -1 their squares straddle zero (lambda_boost >= 0 >= lambda_rot), so a general
// STA rotor is a Lorentz boost times a spatial rotation in dual planes:
//   exp(B) = exp(b_boost) * exp(b_rot)   (the factors commute).
//
// The causal branch reads the GEOMETRIC square gr0(B*B), NOT nrm_sq(B): under the
// P-unify metric nrm_sq(B) is the reverse-norm (= -B^2 for a bivector), so the causal
// character lives in the geometric square, while nrm(B) = sqrt(|bb_s|) gives |a|.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> exp(BiVec4ds<T> const& B)
{
    auto const p = detail::sta4ds_biv_decompose(B);

    if (p.simple) {
        // SIMPLE bivector -> closed-form generalised Euler formula
        return detail::sta4ds_exp_simple(B, p.geom_sq);
    }

    // NON-simple bivector: the two invariant planes commute, so the exponentials factor
    using sta::operator*;
    return detail::sta4ds_exp_simple(p.b_boost, p.lambda_boost) *
           detail::sta4ds_exp_simple(p.b_rot, p.lambda_rot);
}


////////////////////////////////////////////////////////////////////////////////
// exponential of a PSEUDOSCALAR argument
////////////////////////////////////////////////////////////////////////////////
//
// I_4ds^2 = -1, so the power series closes into
//
//     exp(alpha*I_4ds) = cos(alpha) + I_4ds sin(alpha)
//
// Scalar and pseudoscalar are both EVEN in four dimensions, so the result is an
// MVec4ds_E with a zero bivector part.
//
// Two warnings, both different from the 3d case:
//
//  - I_4ds is NOT central here. It commutes with even elements and ANTIcommutes
//    with odd ones, so exp(alpha*I) passes through bivectors but not through
//    vectors or trivectors.
//  - it is not a unit versor: rev(I_4ds) = I_4ds, so exp(a*I) rev(exp(a*I))
//    = exp(2a*I), not 1. Do not feed it to transform() expecting a motion.
//
// Its use is the DUALITY ROTATION of a bivector field: because it commutes with
// bivectors, F -> exp(alpha*I) F maps a bivector to a bivector, rotating the two
// parts of the spacetime split into each other. At alpha = pi/2 it is the Hodge
// map. ega3d's exp(PScalar3d) is the same transformation one algebra down.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> exp(PScalar4ds<T> ps)
{
    // I_4ds is a UNIT pseudoscalar, so the coefficient is the angle itself
    T const alpha = T(ps);
    return MVec4ds_E<T>(Scalar4ds<T>(std::cos(alpha)),
                        BiVec4ds<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)),
                        PScalar4ds<T>(std::sin(alpha)));
}


////////////////////////////////////////////////////////////////////////////////
// rotor logarithm w.r.t. the geometric product -- the inverse of exp()
//
// R's pseudoscalar part s4 = gr4(R) = sinh(f) sin(th) is the non-simplicity indicator:
// it is zero iff R is a SIMPLE rotor (a single rotation OR a single boost).
//
// SIMPLE rotor (|s4| ~ 0): gr2(R) shares B's plane and causal character, with magnitude
// nrm(gr2 R) = |sin(a)| (rotation) or sinh(a) (boost), s = gr0(R):
//   gr2(R) spacelike (rotation): R = cos(a)  + (sin(a)/a)  B  => a = atan2(|sin a|, cos
//   a) gr2(R) timelike  (boost):    R = cosh(a) + (sinh(a)/a) B  => a = asinh(sinh a)
//   gr2(R) lightlike / a ~ 0:    R = 1 + B (or identity)      => B = gr2(R) directly
// then B = (a / nrm(gr2 R)) * gr2(R).
//
// NON-simple rotor (|s4| > 0): R = exp(b_boost) exp(b_rot) is a boost (rapidity f) and a
// spatial rotation (angle th) in orthogonal (dual) planes. Decompose gr2(R) into its two
// invariant simple parts (De Keninck & Roelfs sec. 7); the rotation plane b_rot comes
// from gr2(R) directly, and the boost plane b_boost is recovered as the dual of the
// rotation plane, oriented so that exp reproduces R's g1234 part s4 -- this stays robust
// even at th == pi/2 (where the boost plane vanishes from gr2(R)).
//
// A rotation by ~pi (th -> 0 or pi) shares the usual rotation-logarithm axis ambiguity;
// there the non-simple split degenerates back to the simple branch.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> log(MVec4ds_E<T> const& R)
{
    BiVec4ds<T> const Bv = gr2(R);
    T const s0 = T(gr0(R));
    T const s4 = T(gr4(R)); // = sinh(f) sin(th); zero iff SIMPLE

    if (std::abs(s4) <= detail::safe_epsilon<T>()) {
        // SIMPLE rotor -> single-plane logarithm
        T const bn = T(nrm(Bv)); // |sin a| (rotation) / sinh a (boost) / 0
        if (bn <= detail::safe_epsilon<T>()) {
            // a ~ 0 -> identity (Bv ~ 0); or a lightlike (null) plane, exp(B) = 1 + B
            return Bv;
        }
        if (is_spacelike(Bv)) {
            T const a = std::atan2(bn, s0); // rotation by angle a in [0, pi]
            return (a / bn) * Bv;
        }
        // timelike plane -> Lorentz boost (rapidity a); asinh(sinh a) is well-conditioned
        T const a = std::asinh(bn);
        return (a / bn) * Bv;
    }

    // NON-simple rotor: split gr2(R) into its two invariant planes via the Study
    // projectors, exactly as exp() does for the generating bivector.
    T const bv_s = detail::sta4ds_geom_sq(Bv);
    T const bv_ps = T(gr4(MVec4ds_E<T>(Bv) * Bv));
    T const nSv = std::sqrt(bv_s * bv_s + bv_ps * bv_ps);
    BiVec4ds<T> const projv = (bv_s * Bv - bv_ps * detail::sta4ds_biv_dual(Bv)) / nSv;
    BiVec4ds<T> const sm = T(0.5) * (Bv - projv); // rotation-plane part, sm^2 <= 0
    T const P =
        std::sqrt(std::max(T(0.0), T(0.5) * (bv_s + nSv))); // |sp| = |cos th| sinh f
    T const Q =
        std::sqrt(std::max(T(0.0), T(0.5) * (nSv - bv_s))); // |sm| = cosh f  sin th

    // rotation: angle th in [0, pi], unit plane b_hat_rot = sm / Q  (Q > 0 when |s4| > 0)
    T const th = std::atan2(Q, s0);
    BiVec4ds<T> const bhat_rot = sm / Q;

    // boost: rapidity f from sinh f = sqrt(P^2 + s4^2); the boost plane is the dual of
    // the rotation plane, sign-fixed so exp reproduces R's g1234 coefficient s4.
    T const f = std::asinh(std::sqrt(P * P + s4 * s4));
    BiVec4ds<T> const bhat_boost_cand =
        detail::sta4ds_biv_dual(bhat_rot); // unit, ^2 = +1
    using sta::operator*;
    T const kappa = T(gr4(MVec4ds_E<T>(bhat_boost_cand) * bhat_rot)); // (b+ ^ b-) = +/-1
    T const sgn = (s4 * kappa >= T(0.0)) ? T(1.0) : T(-1.0);

    return f * (sgn * bhat_boost_cand) + th * bhat_rot;
}


////////////////////////////////////////////////////////////////////////////////
// rotor for a spatial rotation by angle theta in the oriented plane B
// (B a spatial bivector, B*B < 0; need not be normalized).
// Apply via the sandwich transform(X, R) = R * X * rev(R): rotates by +theta, i.e.
// it turns the plane's first basis direction toward its second (g1 -> g2 for g12),
// which is the same sense ega3d's get_rotor() produces (e1 -> e2 for e12).
//
// SIGN (this is where the metric bites): the half-angle enters POSITIVELY here,
// exp(+theta/2 * B_hat) -- the SAME sign get_boost() uses. ega3d needs an explicit
// -theta/2 for the equivalent sense; here that flip is already supplied by the metric,
// since the spatial basis vectors square to -1. Copying ega3d's literal -theta/2
// therefore REVERSES the rotation (it did, until 2026-08-16). Guarded by the test case
// "rotation and boost sense (convention gate vs. ega3d)".
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline MVec4ds_E<std::common_type_t<T, U>> get_rotor(BiVec4ds<T> const& B, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype const half = 0.5 * theta;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(std::cos(half)),
                            normalize(B) * std::sin(half));
}


////////////////////////////////////////////////////////////////////////////////
// rotor for a Lorentz boost by rapidity phi in the oriented time-space plane B
// (B a boost bivector, B*B > 0; need not be normalized).
// rapidity relates to velocity by  beta = tanh(phi),  gamma = cosh(phi).
// Apply via the sandwich transform(X, R) = R * X * rev(R): boosts toward the plane's
// first basis direction (g4 -> +g1 for g14), with gamma = cosh(phi).
// (The +phi/2 sign convention is validated against gamma = cosh(phi) in step 2.)
// Same POSITIVE half-argument as get_rotor() -- in STA both builders are
// exp(+p/2 * B_hat); see the SIGN note there for why ega3d's rotor needs -theta/2
// for the equivalent sense and this one does not.
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline MVec4ds_E<std::common_type_t<T, U>> get_boost(BiVec4ds<T> const& B, U phi)
{
    using ctype = std::common_type_t<T, U>;
    ctype const half = 0.5 * phi;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(std::cosh(half)),
                            normalize(B) * std::sinh(half));
}


////////////////////////////////////////////////////////////////////////////////
// sqrt(rotor) w.r.t. the geometric product -- halves the rotation angle / boost rapidity
//
// The principal square root sqrt(R) = Normalize(1 + R) (De Keninck & Roelfs eq. 11)
// bisects the versor arc from the identity to R:
//   spatial rotation  R = cos a  + sin a  B  (B^2 < 0)  ->  cos(a/2)  + sin(a/2)  B
//   Lorentz boost     R = cosh f + sinh f B  (B^2 > 0)  ->  cosh(f/2) + sinh(f/2) B
//
// Normalization is the signature-agnostic Study renormalisation (eq. 24, in
// detail::sta4ds_nearest_rotor): for a SIMPLE rotor 1 + R has a scalar norm and this is
// just (1 + R)/sqrt(<..>_0), but for a NON-simple rotor 1 + R has a non-scalar (Study)
// norm with a g1234 part -- dropping it (a scalar-only versor norm) would return a
// non-rotor, so the full Study normalisation is used. R is first renormalised to a unit
// rotor the same way, so non-unit inputs (e.g. integration drift) are handled.
//
// R = -1 (a 2*pi rotation) gives 1 + R = 0, which has no unique principal root ->
// nearest_rotor returns the identity rotor (the conventional choice).
//
// NOTE: the c+ Study branch is the principal root for <1 + R>_0 >= 0, which holds for any
// rotor near the identity (the renormalisation-drift use case). An extreme rotor
// combining a large boost with a rotation past pi/2 can push <1 + R>_0 < 0 onto the other
// (c-) branch; that regime shares the branch-cut ambiguity of any square root.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> sqrt(MVec4ds_E<T> const& R)
{
    MVec4ds_E<T> const Ru = detail::sta4ds_nearest_rotor(R);     // ensure a unit rotor
    return detail::sta4ds_nearest_rotor(Scalar4ds<T>(1.0) + Ru); // Normalize(1 + R)
}


////////////////////////////////////////////////////////////////////////////////
// separation between two vectors (Lorentzian -- split by causal character)
//
// In a Lorentzian metric there is no single "angle" between two vectors. The two
// physically meaningful cases are kept as separate, guarded functions:
//
//   angle(v1, v2)    -- ordinary EUCLIDEAN angle, both vectors SPACELIKE, [0, pi]
//   rapidity(u1, u2) -- HYPERBOLIC angle (relative rapidity), both TIMELIKE, [0, inf)
//
// The sign handling differs from ega3d because dot() carries the metric: for a
// spacelike vector dot(v,v) = -nrm(v)^2 < 0, while for a timelike vector
// dot(u,u) = +nrm(u)^2 > 0. Mixed (one spacelike, one timelike) and lightlike
// inputs have no well-defined separation here and are rejected.
////////////////////////////////////////////////////////////////////////////////

// Euclidean angle between two SPACELIKE vectors, range 0 <= angle <= pi.
//
// cos(angle) = -dot(v1, v2) / (nrm(v1) * nrm(v2))   (the leading minus undoes the
// negative spacelike metric, so two equal spacelike vectors give angle 0).
//
// PRE: both v1, v2 spacelike (nrm_sq < 0). The result is only geometrically a true
// angle when v1, v2 span a spacelike plane; if the plane they span is indefinite the
// ratio can exceed 1 and is clamped (the angle then loses its literal meaning).
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> angle(Vec4ds<T> const& v1, Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    if (!is_spacelike(v1) || !is_spacelike(v2)) {
        throw std::runtime_error(
            "GA Error: angle(v1, v2) is defined for spacelike vectors only "
            "(use rapidity() for timelike vectors)");
    }
    ctype const nrm_prod = ctype(nrm(v1)) * ctype(nrm(v2));
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
    return std::acos(std::clamp(-ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// Relative rapidity between two TIMELIKE vectors, range 0 <= rapidity < inf.
//
// For two unit future-pointing 4-velocities the Lorentz factor is gamma = dot(u1, u2)
// and the relative rapidity is acosh(gamma). Using |dot| / (nrm * nrm) makes this
// scale-independent and orientation-agnostic, with cosh >= 1 enforced by the clamp.
//
// PRE: both u1, u2 timelike (nrm_sq > 0).
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline std::common_type_t<T, U> rapidity(Vec4ds<T> const& u1, Vec4ds<U> const& u2)
{
    using ctype = std::common_type_t<T, U>;
    if (!is_timelike(u1) || !is_timelike(u2)) {
        throw std::runtime_error(
            "GA Error: rapidity(u1, u2) is defined for timelike vectors only "
            "(use angle() for spacelike vectors)");
    }
    ctype const nrm_prod = ctype(nrm(u1)) * ctype(nrm(u2));
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");
    return std::acosh(std::max(std::abs(ctype(dot(u1, u2))) / nrm_prod, ctype(1.0)));
}


////////////////////////////////////////////////////////////////////////////////
// apply a Lorentz transformation (rotation and/or boost) via the rotor sandwich
//
//     transform(X, R) = R * X * rev(R)
//
// with R a unit rotor (from get_rotor / get_boost / exp). The transformation is
// grade-preserving and preserves the spacetime interval nrm_sq(X).
//
// (direct geometric-product form. A closed-form transform_opt() built from the
//  ga_prdxpr sandwich coefficients follows below; for one-off transforms this
//  direct form is actually faster -- only the batch std::vector overload of
//  transform_opt() wins, by amortizing the rotor-only matrix over many vectors.)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> transform(Vec4ds<T> const& v,
                                                     MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(gr1(R * v * rev(R)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> transform(BiVec4ds<T> const& B,
                                                       MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(gr2(R * B * rev(R)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> transform(TriVec4ds<T> const& t,
                                                        MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    return TriVec4ds<ctype>(gr3(R * t * rev(R)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> transform(MVec4ds<T> const& M,
                                                      MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(R * M * rev(R));
}

// optimized closed-form Lorentz transformation:  X' = R * X * rev(R)
//
// The sandwich collapses to a single matrix M(R) acting on the blade's components,
// quadratic in the rotor coefficients R.c0..R.c7 (the "after transformation"
// coefficients emitted by ga_prdxpr). The matrix is built once by the helpers in
// hd::ga::detail (collecting every distinct R.ci*R.cj product), then applied:
//   - grade 1 (vector)    and grade 3 (trivector) share one 4x4 matrix
//   - grade 2 (bivector)  uses a 6x6 matrix
//
// validated against the direct transform() in the test suite. For a single one-off
// transform the direct transform() is faster (the matrix-build cost is not amortised);
// the closed form only pays off in the std::vector batch overloads below, where one
// matrix is reused across many blades. See ga_test/utilities/bench_sta4ds_transform.

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> transform_opt(Vec4ds<T> const& v,
                                                         MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    auto const k = detail::sta_rotor_xf_mat_vec<ctype>(R);
    return Vec4ds<ctype>(k[0] * v.x + k[1] * v.y + k[2] * v.z + k[3] * v.w,
                         k[4] * v.x + k[5] * v.y + k[6] * v.z + k[7] * v.w,
                         k[8] * v.x + k[9] * v.y + k[10] * v.z + k[11] * v.w,
                         k[12] * v.x + k[13] * v.y + k[14] * v.z + k[15] * v.w);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec4ds<std::common_type_t<T, U>> transform_opt(TriVec4ds<T> const& t,
                                                            MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    // trivectors transform with the same 4x4 matrix as vectors in sta
    auto const k = detail::sta_rotor_xf_mat_vec<ctype>(R);
    return TriVec4ds<ctype>(k[0] * t.x + k[1] * t.y + k[2] * t.z + k[3] * t.w,
                            k[4] * t.x + k[5] * t.y + k[6] * t.z + k[7] * t.w,
                            k[8] * t.x + k[9] * t.y + k[10] * t.z + k[11] * t.w,
                            k[12] * t.x + k[13] * t.y + k[14] * t.z + k[15] * t.w);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> transform_opt(BiVec4ds<T> const& B,
                                                           MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    auto const k = detail::sta_rotor_xf_mat_bivec<ctype>(R);
    return BiVec4ds<ctype>(k[0] * B.vx + k[1] * B.vy + k[2] * B.vz + k[3] * B.mx +
                               k[4] * B.my + k[5] * B.mz,
                           k[6] * B.vx + k[7] * B.vy + k[8] * B.vz + k[9] * B.mx +
                               k[10] * B.my + k[11] * B.mz,
                           k[12] * B.vx + k[13] * B.vy + k[14] * B.vz + k[15] * B.mx +
                               k[16] * B.my + k[17] * B.mz,
                           k[18] * B.vx + k[19] * B.vy + k[20] * B.vz + k[21] * B.mx +
                               k[22] * B.my + k[23] * B.mz,
                           k[24] * B.vx + k[25] * B.vy + k[26] * B.vz + k[27] * B.mx +
                               k[28] * B.my + k[29] * B.mz,
                           k[30] * B.vx + k[31] * B.vy + k[32] * B.vz + k[33] * B.mx +
                               k[34] * B.my + k[35] * B.mz);
}

// batch Lorentz transformation of many blades by the SAME rotor R.
//
// This is the variant where the closed form actually pays off: the matrix M(R) depends
// only on the rotor, so it is built once and then applied to every blade (one
// matrix-vector product each). For a single one-off transform prefer transform().

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
std::vector<Vec4ds<std::common_type_t<T, U>>>
transform_opt(std::vector<Vec4ds<T>> const& vecs, MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    auto const k = detail::sta_rotor_xf_mat_vec<ctype>(R);
    std::vector<Vec4ds<ctype>> res;
    res.reserve(vecs.size());
    for (auto const& v : vecs) {
        res.emplace_back(
            Vec4ds<ctype>(k[0] * v.x + k[1] * v.y + k[2] * v.z + k[3] * v.w,
                          k[4] * v.x + k[5] * v.y + k[6] * v.z + k[7] * v.w,
                          k[8] * v.x + k[9] * v.y + k[10] * v.z + k[11] * v.w,
                          k[12] * v.x + k[13] * v.y + k[14] * v.z + k[15] * v.w));
    }
    return res;
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
std::vector<TriVec4ds<std::common_type_t<T, U>>>
transform_opt(std::vector<TriVec4ds<T>> const& tris, MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    auto const k = detail::sta_rotor_xf_mat_vec<ctype>(R);
    std::vector<TriVec4ds<ctype>> res;
    res.reserve(tris.size());
    for (auto const& t : tris) {
        res.emplace_back(
            TriVec4ds<ctype>(k[0] * t.x + k[1] * t.y + k[2] * t.z + k[3] * t.w,
                             k[4] * t.x + k[5] * t.y + k[6] * t.z + k[7] * t.w,
                             k[8] * t.x + k[9] * t.y + k[10] * t.z + k[11] * t.w,
                             k[12] * t.x + k[13] * t.y + k[14] * t.z + k[15] * t.w));
    }
    return res;
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
std::vector<BiVec4ds<std::common_type_t<T, U>>>
transform_opt(std::vector<BiVec4ds<T>> const& bivecs, MVec4ds_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    auto const k = detail::sta_rotor_xf_mat_bivec<ctype>(R);
    std::vector<BiVec4ds<ctype>> res;
    res.reserve(bivecs.size());
    for (auto const& B : bivecs) {
        res.emplace_back(BiVec4ds<ctype>(k[0] * B.vx + k[1] * B.vy + k[2] * B.vz +
                                             k[3] * B.mx + k[4] * B.my + k[5] * B.mz,
                                         k[6] * B.vx + k[7] * B.vy + k[8] * B.vz +
                                             k[9] * B.mx + k[10] * B.my + k[11] * B.mz,
                                         k[12] * B.vx + k[13] * B.vy + k[14] * B.vz +
                                             k[15] * B.mx + k[16] * B.my + k[17] * B.mz,
                                         k[18] * B.vx + k[19] * B.vy + k[20] * B.vz +
                                             k[21] * B.mx + k[22] * B.my + k[23] * B.mz,
                                         k[24] * B.vx + k[25] * B.vy + k[26] * B.vz +
                                             k[27] * B.mx + k[28] * B.my + k[29] * B.mz,
                                         k[30] * B.vx + k[31] * B.vy + k[32] * B.vz +
                                             k[33] * B.mx + k[34] * B.my + k[35] * B.mz));
    }
    return res;
}


////////////////////////////////////////////////////////////////////////////////
// spacetime split of a vector x relative to a unit timelike observer u (u*u = +1):
//
//     x * u = (x . u) + (x ^ u)
//
//   time_split(x, u)  = x . u   (scalar)   -> time component of x measured by u
//   space_split(x, u) = x ^ u   (bivector) -> relative spatial vector, a "relative
//                                vector" (g_k4-type bivector) in the u-frame
//
// Recover x from its parts via  x = gr1( (time_split + space_split) * u )  (u*u=1).
// The standard observer is the time direction u = g4.
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Scalar4ds<std::common_type_t<T, U>> time_split(Vec4ds<T> const& x,
                                                         Vec4ds<U> const& u)
{
    return dot(x, u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> space_split(Vec4ds<T> const& x,
                                                         Vec4ds<U> const& u)
{
    return wdg(x, u);
}


////////////////////////////////////////////////////////////////////////////////
// spacetime split of a bivector F relative to a unit timelike observer u (u*u = +1).
//
// F decomposes into the part anti-commuting with u -- a "relative vector", spanned by
// the g_k4 boost bivectors (the electric-type part) -- and the part commuting with u
// -- a "relative bivector", spanned by the g_jk rotation bivectors (the magnetic-type
// part):
//
//   rel_vec_split(F, u)   = 0.5 * (F - gr2(u * F * u))   (anti-commuting, g_k4-type)
//   rel_bivec_split(F, u) = 0.5 * (F + gr2(u * F * u))   (commuting,      g_jk-type)
//
// Their sum is F. Under a boost the two parts mix (the electromagnetic field
// transformation). The standard observer is the time direction u = g4.
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> rel_vec_split(BiVec4ds<T> const& F,
                                                           Vec4ds<U> const& u)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(0.5 * (F - gr2(u * F * u)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> rel_bivec_split(BiVec4ds<T> const& F,
                                                             Vec4ds<U> const& u)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(0.5 * (F + gr2(u * F * u)));
}


////////////////////////////////////////////////////////////////////////////////
// invariant decomposition of a bivector (the OBSERVER-INDEPENDENT split)
////////////////////////////////////////////////////////////////////////////////
//
// Every bivector B splits uniquely into two orthogonal, commuting, SIMPLE bivectors
//
//     B = boost_part(B) + rot_part(B)
//
// whose squares straddle zero: boost_part(B)^2 >= 0 >= rot_part(B)^2. The time-like
// part generates a Lorentz boost, the space-like part a spatial rotation, and because
// they commute the rotor factors: exp(B) = exp(boost_part(B)) * exp(rot_part(B)).
//
// Both parts are built from the two Lorentz invariants of B -- its geometric square
// <B^2>_0 and its wedge square B^B -- so the decomposition depends on B ALONE. That is
// the difference to rel_vec_split() / rel_bivec_split() above, which split relative to a
// chosen time direction u and mix under a boost.
//
// A SIMPLE bivector (B ^ B == 0) already lies in one plane and is returned unchanged in
// the slot matching its causal character, the other part being zero.
////////////////////////////////////////////////////////////////////////////////

// true if the bivector lies in a single plane, i.e. B ^ B == 0 (always true in 3d, but
// NOT in 4d: a generic sta bivector is a sum of two volutors)
template <typename T>
    requires(numeric_type<T>)
inline bool is_simple(BiVec4ds<T> const& B)
{
    return std::abs(T(gr4(MVec4ds_E<T>(B) * B))) <= detail::safe_epsilon<T>();
}

// the time-like (boost-generating) part; squares to >= 0
template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> boost_part(BiVec4ds<T> const& B)
{
    return detail::sta4ds_biv_decompose(B).b_boost;
}

// the space-like (rotation-generating) part; squares to <= 0
template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> rot_part(BiVec4ds<T> const& B)
{
    return detail::sta4ds_biv_decompose(B).b_rot;
}


////////////////////////////////////////////////////////////////////////////////
// projections and rejections (geometric-product based, as in ega3d)
//
//   project_onto(a, b): component of a parallel to / lying in b
//   reject_from(a, b) = a - project_onto(a, b): component perpendicular to b
//
// PRE: the target b must be non-null (nrm_sq(b) != 0). Projection onto a null
// (lightlike) blade is undefined here (inv() divides by nrm_sq) -- there is no
// Euclidean analog. project_onto + reject_from == a, and project_onto(a,b) lies
// in b (wdg(project_onto(a,b), b) == 0).
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto a vector v2
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> project_onto(Vec4ds<T> const& v1,
                                                        Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of a vector v1 from a vector v2
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reject_from(Vec4ds<T> const& v1,
                                                       Vec4ds<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(v1 - project_onto(v1, v2));
}

// projection of a vector v onto a bivector B (a 2-plane)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> project_onto(Vec4ds<T> const& v,
                                                        BiVec4ds<U> const& B)
{
    return gr1((B >> v) * inv(B));
}

// rejection of a vector v from a bivector B (a 2-plane)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reject_from(Vec4ds<T> const& v,
                                                       BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(v - project_onto(v, B));
}


////////////////////////////////////////////////////////////////////////////////
// reflections (geometric-product based, as in ega3d)
////////////////////////////////////////////////////////////////////////////////
// Macdonald p. 129: reflecting a j-blade u in the k-dimensional subspace B gives
//
//   u_reflected = (-1)^[j*(k+1)] * B * u * inv(B)
//
// The reflecting subspace B must be non-null (nrm_sq(B) != 0), since inv(B)
// divides by nrm_sq(B); reflection in a lightlike subspace is undefined here.
// All these reflections preserve nrm_sq (they are sandwiches by an invertible
// versor) and are involutions (applying the same reflection twice is identity).
// In STA they also realise the discrete spacetime symmetries: reflecting on the
// spatial hyperplane (normal g4) is parity P, reflecting on the time vector g4
// flips the three spatial axes, etc.
////////////////////////////////////////////////////////////////////////////////

// reflect a vector v on a hyperplane given by its (non-null) normal vector nB
//
// hyperplane: an (n-1)-dimensional subspace; here the 3d subspace orthogonal to nB.
// The reflection flips the component of v along nB and keeps the in-plane part:
//   v_reflected = -nB * v * inv(nB)
//
// This is the complement of reflect_on_vec(v, b) below: the leading minus sign
// turns "reflect onto the line b" into "reflect across the hyperplane normal to nB".
//
// HINT: nB = dual(t) for a hyperplane trivector t (with normalized nB).
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reflect_on(Vec4ds<T> const& v,
                                                      Vec4ds<U> const& nB)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(gr1(-nB * v * inv(nB)));
}

// reflect a vector v on a hyperplane given directly as a (non-null) trivector t
//
// In 4d a hyperplane is 3-dimensional, i.e. a trivector. Macdonald: j=1, k=3 =>
// sign (-1)^[1*4] = +1, hence v_reflected = t * v * inv(t).
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reflect_on(Vec4ds<T> const& v,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(gr1(t * v * inv(t)));
}

// reflect a vector v in an arbitrary (non-null) bivector B, i.e. a 2-plane
// Macdonald: j=1, k=2 => sign (-1)^[1*3] = -1, hence v_reflected = -B * v * inv(B).
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reflect_on(Vec4ds<T> const& v,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(gr1(-B * v * inv(B)));
}

// reflect a bivector UB in an arbitrary (non-null) bivector B (both 2-planes)
// Macdonald: j=2, k=2 => sign (-1)^[2*3] = +1, hence UB_reflected = B * UB * inv(B).
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec4ds<std::common_type_t<T, U>> reflect_on(BiVec4ds<T> const& UB,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec4ds<ctype>(gr2(B * UB * inv(B)));
}

// reflect a vector v on another (non-null) vector b
// Macdonald: j=1, k=1 => sign (-1)^[1*2] = +1, hence v_reflected = b * v * inv(b).
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec4ds<std::common_type_t<T, U>> reflect_on_vec(Vec4ds<T> const& v,
                                                          Vec4ds<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec4ds<ctype>(gr1(b * v * inv(b)));
}


////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar4ds<T> a, Scalar4ds<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero scalars are congruent (represent the same 0-dimensional subspace)
    return true;
}

// For vectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Vec4ds<T> const& a, Vec4ds<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, tolerance);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec4ds<T> const& a, BiVec4ds<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<6>({value_t(a.vx), value_t(a.vy), value_t(a.vz),
                                        value_t(a.mx), value_t(a.my), value_t(a.mz)},
                                       {value_t(b.vx), value_t(b.vy), value_t(b.vz),
                                        value_t(b.mx), value_t(b.my), value_t(b.mz)},
                                       tolerance);
}

// For trivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(TriVec4ds<T> const& a, TriVec4ds<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, tolerance);
}

// For pseudoscalars: all non-zero pseudoscalars in 4ds represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar4ds<T> a, PScalar4ds<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars in 4ds are congruent (represent the full 4D space)
    return true;
}


////////////////////////////////////////////////////////////////////////////////
// is_close(): equality within a RELATIVE tolerance
////////////////////////////////////////////////////////////////////////////////
//
// Same question as operator==, but with the threshold scaled by the operands (see
// detail::coeffs_close). Reach for it whenever the values carry a physical scale:
// operator== measures against an absolute eps, which cannot resolve anything once
// coordinates grow large, where a single ulp already exceeds it.
//
// Distinct from is_congruent(), which allows an arbitrary scale factor between the
// operands and so answers "same subspace", not "same value".

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Scalar4ds<T> a, Scalar4ds<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Vec4ds<T> const& a, Vec4ds<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(BiVec4ds<T> const& a, BiVec4ds<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<6>({value_t(a.vx), value_t(a.vy), value_t(a.vz),
                                    value_t(a.mx), value_t(a.my), value_t(a.mz)},
                                   {value_t(b.vx), value_t(b.vy), value_t(b.vz),
                                    value_t(b.mx), value_t(b.my), value_t(b.mz)},
                                   rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(TriVec4ds<T> const& a, TriVec4ds<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(PScalar4ds<T> a, PScalar4ds<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec4ds_E<T> const& a, MVec4ds_E<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<8>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec4ds_U<T> const& a, MVec4ds_U<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<8>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec4ds<T> const& a, MVec4ds<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<16>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7), value_t(a.c8), value_t(a.c9),
         value_t(a.c10), value_t(a.c11), value_t(a.c12), value_t(a.c13), value_t(a.c14),
         value_t(a.c15)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7), value_t(b.c8), value_t(b.c9),
         value_t(b.c10), value_t(b.c11), value_t(b.c12), value_t(b.c13), value_t(b.c14),
         value_t(b.c15)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(DualNum4ds<T> const& a, DualNum4ds<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<2>({value_t(a.c0), value_t(a.c1)},
                                   {value_t(b.c0), value_t(b.c1)}, rel_tol);
}


////////////////////////////////////////////////////////////////////////////////
// is_same_transform(): do two rotors describe the same Lorentz transformation?
////////////////////////////////////////////////////////////////////////////////
//
// NOT the same question as R1 == R2. Rotors DOUBLE-COVER the Lorentz transformations: the
// sandwich
//
//     X -> R (x) X (x) rev(R)
//
// is quadratic in R, so R and -R act identically on every object while comparing unequal
// component by component. A rotation by 2*pi returns -1 rather than 1 for exactly this
// reason, and a rotor obtained by a different route (an exp/log round trip, a chain of
// compositions, a sqrt) routinely comes back negated. Testing rotors with operator==
// therefore asks about the REPRESENTATION; this asks about the transformation itself.
//
// It is decided by what the rotors DO -- where they send the basis directions, which
// together pin the transformation uniquely -- rather than by comparing coefficients up
// to a sign. That also keeps the answer right for rotors that are not unit: the sandwich
// scales by |R|^2, so a scaled rotor stretches what it acts on and is NOT the same
// Lorentz transformation; the comparisons below say so.
//
// Named transform rather than motion or rotation: an STA rotor covers boosts as well
// as spatial rotations, matching this algebra's transform() sandwich verb.
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_same_transform(MVec4ds_E<T> const& R1, MVec4ds_E<U> const& R2,
                       value_t rel_tol = eps_congruent)
{
    return is_close(transform(g1_4ds, R1), transform(g1_4ds, R2), rel_tol) &&
           is_close(transform(g2_4ds, R1), transform(g2_4ds, R2), rel_tol) &&
           is_close(transform(g3_4ds, R1), transform(g3_4ds, R2), rel_tol) &&
           is_close(transform(g4_4ds, R1), transform(g4_4ds, R2), rel_tol);
}

} // namespace hd::ga::sta
