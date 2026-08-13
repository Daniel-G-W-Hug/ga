#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <algorithm> // std::max (detail::coeffs_congruent)
#include <array>     // detail::coeffs_congruent
#include <cmath>     // std::abs (detail::coeffs_congruent)
#include <cstddef>   // std::size_t
#include <limits>    // std::numeric_limits

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

// select the floating point type used for scalars, vector and bivector components, etc.
// here, i.e. choose one of the two following definitions:
// using value_t = float;
using value_t = double;

// eps is used internally for equality testing (5 * machine epsilon)
constexpr auto eps = 5.0 * std::numeric_limits<value_t>::epsilon();

// Relative tolerance for is_congruent / is_close. Congruence is almost always asked of
// two INDEPENDENTLY COMPUTED geometric quantities (results of wedge / dual / motor
// chains), whose difference has accumulated error far beyond a single rounding -- so
// these tests use a geometric relative tolerance, not machine eps. The factor is the
// error budget: about 4500 roundings' worth, i.e. some 3-4 significant digits given up
// to accumulation, while staying far below any real geometric signal.
//
// It MUST be expressed in machine epsilons rather than as a literal, because value_t is
// a switch (see above). The historical literal 1.0e-12 is right for double -- it is
// 4.5e3 machine epsilons -- but for float it lands 1e5 times BELOW machine epsilon, so a
// single ulp of rounding already exceeds it and is_congruent / is_close would answer
// "different" for every pair of independently computed values. Deriving it from
// numeric_limits keeps the intended error budget under either setting:
//
//     double: 4.5e3 * 2.22e-16 = 9.99e-13   (the historical value)
//     float : 4.5e3 * 1.19e-07 = 5.36e-04
//
constexpr auto eps_congruent = 4.5e3 * std::numeric_limits<value_t>::epsilon();

namespace detail {

// Congruence of two coefficient lists: true iff a == k*b for some non-zero scale k (the
// same direction, up to sign and magnitude), within the RELATIVE tolerance rel_tol. Zero
// lists are handled component-wise -- both zero -> congruent, exactly one zero -> not --
// using an absolute eps floor independent of rel_tol (robust for degenerate PGA metrics,
// where bulk_nrm_sq can vanish for a non-zero element).
//
// The scale factor k is pivoted on the LARGEST-magnitude entry of b: dividing by the
// dominant component keeps k well conditioned, so the rounding of a SMALL pivot is never
// amplified onto the dominant entries (the failure mode of a first-non-zero pivot, where
// the amplification factor max_component / pivot_component is unbounded).
template <std::size_t N>
inline bool coeffs_congruent(std::array<value_t, N> const& a,
                             std::array<value_t, N> const& b,
                             value_t rel_tol = eps_congruent)
{
    value_t a_max = 0.0;
    value_t b_max = 0.0;
    std::size_t i_piv = 0;
    for (std::size_t i = 0; i < N; ++i) {
        a_max = std::max(a_max, std::abs(a[i]));
        if (std::abs(b[i]) > b_max) {
            b_max = std::abs(b[i]);
            i_piv = i;
        }
    }

    bool const a_zero = a_max < eps;
    bool const b_zero = b_max < eps;
    if (a_zero && b_zero) return true;  // both effectively zero
    if (a_zero || b_zero) return false; // only one is zero

    value_t const k = a[i_piv] / b[i_piv];
    value_t const tol = rel_tol * std::max(std::max(a_max, b_max), value_t(1.0));
    for (std::size_t i = 0; i < N; ++i) {
        if (std::abs(a[i] - k * b[i]) >= tol) return false;
    }
    return true;
}

// Component-wise EQUALITY within a RELATIVE tolerance -- the k == 1 case of
// coeffs_congruent above (which asks the weaker question "equal up to some scale
// factor").
//
// Why this exists next to operator==: operator== compares against the ABSOLUTE eps, so
// it answers "are these the same" only while the operands are of order one. Two
// coordinates of order 1e6 differ by ~1e-10 in their last bit alone, so an absolute eps
// rejects values that agree to every digit a double can carry; two coordinates of order
// 1e-9 pass it no matter how different they are. Scaling the threshold by the operands
// removes both failure modes:
//
//     |a[i] - b[i]| < rel_tol * max(|a|_inf, |b|_inf, 1)
//
// The floor of 1 keeps the test absolute near zero, so a component that should vanish is
// still compared against rel_tol rather than against nothing.
//
// Use this whenever the operands carry a physical scale (positions in metres, a motor's
// translation part, anything derived from them); keep operator== for values known to be
// of order one, and is_congruent when only the direction matters, not the magnitude.
template <std::size_t N>
inline bool coeffs_close(std::array<value_t, N> const& a, std::array<value_t, N> const& b,
                         value_t rel_tol = eps_congruent)
{
    value_t scale = 1.0;
    for (std::size_t i = 0; i < N; ++i) {
        scale = std::max(scale, std::max(std::abs(a[i]), std::abs(b[i])));
    }

    value_t const tol = rel_tol * scale;
    for (std::size_t i = 0; i < N; ++i) {
        if (std::abs(a[i] - b[i]) >= tol) return false;
    }
    return true;
}

} // namespace detail

} // namespace hd::ga
