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
//   - exp(BiVec)                      -> rotor exponential of a (simple) bivector
//   - log(rotor)                      -> bivector log of a (simple) rotor (exp inverse)
//   - get_rotor(plane, angle)         -> rotor for a spatial rotation
//   - get_boost(plane, phi)           -> rotor for a Lorentz boost (rapidity phi)
//   - sqrt(rotor)                     -> rotor halving the rotation angle / rapidity
//   - angle() / rapidity()            -> separation of two vectors (spacelike / timelike)
//   - transform(X, R)                 -> apply a rotor via the sandwich R*X*rev(R)
//   - transform_opt(X, R)             -> closed-form transform, vec/bivec/trivec
//                                        (scalar + std::vector batch overloads)
//   - time_split() / space_split()    -> spacetime split of a vector (time + rel. space)
//   - rel_vec_split() / rel_bivec_split() -> spacetime split of a bivector (E / B parts)
//   - project_onto() / reject_from()  -> projection / rejection (onto vector or bivector)
//   - reflect_on() / reflect_on_vec() -> reflections (hyperplane, 2-plane, vector)
//   - is_congruent()                  -> same subspace up to a scalar factor
//
// TODO (next steps): general (non-simple) rotor support in exp() / log() / sqrt() via the
//   invariant bivector decomposition (log() currently handles SIMPLE rotors, matching
//   exp()).
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// rotor exponential w.r.t. the geometric product
//
// For a SIMPLE bivector B (B^B == 0, so B*B is a pure scalar = the geometric square
// s = B^2 = gr0(B*B)). Its sign is the plane's causal character; the magnitude is
// a = nrm(B) = sqrt(|s|):
//   B spacelike (s < 0):  exp(B) = cos(a)  + (B/a) sin(a)    (rotation)
//   B timelike  (s > 0):  exp(B) = cosh(a) + (B/a) sinh(a)   (Lorentz boost)
//   B lightlike (s == 0): exp(B) = 1 + B
//
// The branch is taken from is_spacelike / is_timelike / is_lightlike (which read the
// geometric square gr0(B*B)), NOT from nrm_sq(B): under the P-unify metric nrm_sq(B)
// is the reverse-norm (= -B^2 for a bivector), so the causal character lives in the
// geometric square, while nrm(B) = sqrt(|nrm_sq(B)|) = sqrt(|s|) still gives |a|.
//
// TODO: general (non-simple) bivectors, where B*B also has a pseudoscalar part.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> exp(BiVec4ds<T> const& B)
{
    T const a =
        T(nrm(B)); // sqrt(|B^2|): rotation angle (spacelike) / rapidity (timelike)
    if (is_lightlike(B) || a <= detail::safe_epsilon<T>()) {
        // null/zero plane -> identity rotor (covers the lightlike case to first order)
        return MVec4ds_E<T>(Scalar4ds<T>(1.0), B);
    }
    if (is_spacelike(B)) {
        // spacelike plane -> circular (rotation)
        return MVec4ds_E<T>(Scalar4ds<T>(std::cos(a)), (std::sin(a) / a) * B);
    }
    // timelike plane -> hyperbolic (Lorentz boost)
    return MVec4ds_E<T>(Scalar4ds<T>(std::cosh(a)), (std::sinh(a) / a) * B);
}


////////////////////////////////////////////////////////////////////////////////
// rotor logarithm w.r.t. the geometric product -- the inverse of exp()
//
// For a SIMPLE rotor R = exp(B) with a simple bivector B (B^B == 0), recovers B. The
// rotor's bivector part gr2(R) shares B's plane and causal character, with magnitude
// nrm(gr2 R) = |sin(a)| (rotation) or sinh(a) (boost), where a = nrm(B). Branch on the
// causal character of gr2(R) (matching exp()), with s = gr0(R):
//   gr2(R) spacelike (rotation): R = cos(a)  + (sin(a)/a)  B  => a = atan2(|sin a|, cos
//   a) gr2(R) timelike  (boost):    R = cosh(a) + (sinh(a)/a) B  => a = acosh(cosh a)
//   gr2(R) lightlike / a ~ 0:    R = 1 + B (or identity)      => B = gr2(R) directly
// then B = (a / nrm(gr2 R)) * gr2(R). A rotation by ~pi leaves the plane orientation
// ambiguous (the usual rotation-logarithm axis ambiguity).
//
// TODO: general (non-simple) rotors (where gr2(R)^gr2(R) != 0, i.e. R also has a
// pseudoscalar part) need the invariant bivector decomposition; not handled here.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> log(MVec4ds_E<T> const& R)
{
    BiVec4ds<T> const Bv = gr2(R);
    T const bn =
        T(nrm(Bv)); // |sin(a)| (rotation) / sinh(a) (boost) / 0 (identity or null)
    if (bn <= detail::safe_epsilon<T>()) {
        // a ~ 0 -> identity (Bv ~ 0); or a lightlike (null) plane where exp(B) = 1 + B
        return Bv;
    }
    T const s = T(gr0(R)); // cos(a) (rotation) or cosh(a) (boost)
    if (is_spacelike(Bv)) {
        // spacelike plane -> circular (rotation by angle a)
        T const a = std::atan2(bn, s);
        return (a / bn) * Bv;
    }
    // timelike plane -> hyperbolic (Lorentz boost, rapidity a)
    T const a = std::acosh(s >= T(1.0) ? s : T(1.0));
    return (a / bn) * Bv;
}


////////////////////////////////////////////////////////////////////////////////
// rotor for a spatial rotation by angle theta in the oriented plane B
// (B a spatial bivector, B*B < 0; need not be normalized).
// Apply via the sandwich transform(X, R) = R * X * rev(R): rotates by +theta.
// Half-angle and sign convention match ega3d get_rotor().
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
inline MVec4ds_E<std::common_type_t<T, U>> get_rotor(BiVec4ds<T> const& B, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype const half = -0.5 * theta;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(std::cos(half)),
                            normalize(B) * std::sin(half));
}


////////////////////////////////////////////////////////////////////////////////
// rotor for a Lorentz boost by rapidity phi in the oriented time-space plane B
// (B a boost bivector, B*B > 0; need not be normalized).
// rapidity relates to velocity by  beta = tanh(phi),  gamma = cosh(phi).
// Apply via the sandwich transform(X, R) = R * X * rev(R).
// (The +phi/2 sign convention is validated against gamma = cosh(phi) in step 2.)
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
// For a SIMPLE unit rotor R the "versor average" sqrt(R) = normalize(1 + R) bisects the
// versor arc from the identity to R, then renormalises, yielding exp(0.5 * log R):
//   spatial rotation  R = cos a  + sin a  B  (B^2 < 0)  ->  cos(a/2)  + sin(a/2)  B
//   Lorentz boost     R = cosh f + sinh f B  (B^2 > 0)  ->  cosh(f/2) + sinh(f/2) B
//
// The renormalisation uses the VERSOR norm  |X| = sqrt(gr0(rev(X) * X))  (the scalar
// part of rev(X) * X), which equals 1 for a unit STA rotor. This is deliberately NOT
// the grade-wise nrm_sq(MVec4ds_E): for these rotors that grade sum is cos(2a) / cosh(2f)
// (signed Lorentzian bivector norm), not the versor norm, so it must not be used here.
//
// PRE: R is a simple rotor (a single rotation or boost plane) -- e.g. the output of
//      get_rotor / get_boost / exp of a simple bivector. General (non-simple) rotors,
//      whose bivector carries a g1234 part, are not handled (mirrors exp()'s scope).
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> sqrt(MVec4ds_E<T> const& R)
{
    MVec4ds_E<T> M{R};

    // ensure a unit rotor (the bisection assumes |R| = 1 in the versor norm)
    T const r_nsq = T(gr0(rev(M) * M));
    if (r_nsq > detail::safe_epsilon<T>() &&
        std::abs(r_nsq - T(1.0)) > detail::safe_epsilon<T>()) {
        M = M / std::sqrt(r_nsq);
    }

    // X = 1 + R, renormalised by the versor norm sqrt(gr0(rev(X) * X))
    MVec4ds_E<T> const X = Scalar4ds<T>(1.0) + M;
    T const n_sq = T(gr0(rev(X) * X));
    if (n_sq <= detail::safe_epsilon<T>()) {
        // degenerate: R = -1 (a 2*pi rotation) -> 1 + R = 0 has no unique direction;
        // return the identity rotor (the conventional principal square root)
        return MVec4ds_E<T>(Scalar4ds<T>(1.0), BiVec4ds<T>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                            PScalar4ds<T>(0.0));
    }
    return X / std::sqrt(n_sq);
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
bool is_congruent(Vec4ds<T> const& a, Vec4ds<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.x) < tolerance) && (std::abs(a.y) < tolerance) &&
                     (std::abs(a.z) < tolerance) && (std::abs(a.w) < tolerance);
    bool b_is_zero = (std::abs(b.x) < tolerance) && (std::abs(b.y) < tolerance) &&
                     (std::abs(b.z) < tolerance) && (std::abs(b.w) < tolerance);

    if (a_is_zero && b_is_zero) {
        return true; // Both are effectively zero
    }
    if (a_is_zero || b_is_zero) {
        return false; // Only one is zero
    }

    // Find scale factor k where a = k*b, checking all components
    ctype k = 0.0;
    bool k_found = false;

    // Find first non-zero component pair to establish k
    if (std::abs(b.x) > tolerance) {
        k = a.x / b.x;
        k_found = true;
    }
    else if (std::abs(b.y) > tolerance) {
        k = a.y / b.y;
        k_found = true;
    }
    else if (std::abs(b.z) > tolerance) {
        k = a.z / b.z;
        k_found = true;
    }
    else if (std::abs(b.w) > tolerance) {
        k = a.w / b.w;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol = tolerance * std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z),
                                            std::abs(a.w), std::abs(b.x), std::abs(b.y),
                                            std::abs(b.z), std::abs(b.w), value_t(1.0)});
    return (std::abs(a.x - k * b.x) < rel_tol) && (std::abs(a.y - k * b.y) < rel_tol) &&
           (std::abs(a.z - k * b.z) < rel_tol) && (std::abs(a.w - k * b.w) < rel_tol);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec4ds<T> const& a, BiVec4ds<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.vx) < tolerance) && (std::abs(a.vy) < tolerance) &&
                     (std::abs(a.vz) < tolerance) && (std::abs(a.mx) < tolerance) &&
                     (std::abs(a.my) < tolerance) && (std::abs(a.mz) < tolerance);
    bool b_is_zero = (std::abs(b.vx) < tolerance) && (std::abs(b.vy) < tolerance) &&
                     (std::abs(b.vz) < tolerance) && (std::abs(b.mx) < tolerance) &&
                     (std::abs(b.my) < tolerance) && (std::abs(b.mz) < tolerance);

    if (a_is_zero && b_is_zero) {
        return true; // Both are effectively zero
    }
    if (a_is_zero || b_is_zero) {
        return false; // Only one is zero
    }

    // Find scale factor k where a = k*b, checking all components
    ctype k = 0.0;
    bool k_found = false;

    // Find first non-zero component pair to establish k
    if (std::abs(b.vx) > tolerance) {
        k = a.vx / b.vx;
        k_found = true;
    }
    else if (std::abs(b.vy) > tolerance) {
        k = a.vy / b.vy;
        k_found = true;
    }
    else if (std::abs(b.vz) > tolerance) {
        k = a.vz / b.vz;
        k_found = true;
    }
    else if (std::abs(b.mx) > tolerance) {
        k = a.mx / b.mx;
        k_found = true;
    }
    else if (std::abs(b.my) > tolerance) {
        k = a.my / b.my;
        k_found = true;
    }
    else if (std::abs(b.mz) > tolerance) {
        k = a.mz / b.mz;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol =
        tolerance *
        std::max({std::abs(a.vx), std::abs(a.vy), std::abs(a.vz), std::abs(a.mx),
                  std::abs(a.my), std::abs(a.mz), std::abs(b.vx), std::abs(b.vy),
                  std::abs(b.vz), std::abs(b.mx), std::abs(b.my), std::abs(b.mz),
                  value_t(1.0)});
    return (std::abs(a.vx - k * b.vx) < rel_tol) &&
           (std::abs(a.vy - k * b.vy) < rel_tol) &&
           (std::abs(a.vz - k * b.vz) < rel_tol) &&
           (std::abs(a.mx - k * b.mx) < rel_tol) &&
           (std::abs(a.my - k * b.my) < rel_tol) && (std::abs(a.mz - k * b.mz) < rel_tol);
}

// For trivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(TriVec4ds<T> const& a, TriVec4ds<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.x) < tolerance) && (std::abs(a.y) < tolerance) &&
                     (std::abs(a.z) < tolerance) && (std::abs(a.w) < tolerance);
    bool b_is_zero = (std::abs(b.x) < tolerance) && (std::abs(b.y) < tolerance) &&
                     (std::abs(b.z) < tolerance) && (std::abs(b.w) < tolerance);

    if (a_is_zero && b_is_zero) {
        return true; // Both are effectively zero
    }
    if (a_is_zero || b_is_zero) {
        return false; // Only one is zero
    }

    // Find scale factor k where a = k*b, checking all components
    ctype k = 0.0;
    bool k_found = false;

    // Find first non-zero component pair to establish k
    if (std::abs(b.x) > tolerance) {
        k = a.x / b.x;
        k_found = true;
    }
    else if (std::abs(b.y) > tolerance) {
        k = a.y / b.y;
        k_found = true;
    }
    else if (std::abs(b.z) > tolerance) {
        k = a.z / b.z;
        k_found = true;
    }
    else if (std::abs(b.w) > tolerance) {
        k = a.w / b.w;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol = tolerance * std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z),
                                            std::abs(a.w), std::abs(b.x), std::abs(b.y),
                                            std::abs(b.z), std::abs(b.w), value_t(1.0)});
    return (std::abs(a.x - k * b.x) < rel_tol) && (std::abs(a.y - k * b.y) < rel_tol) &&
           (std::abs(a.z - k * b.z) < rel_tol) && (std::abs(a.w - k * b.w) < rel_tol);
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

} // namespace hd::ga::sta
