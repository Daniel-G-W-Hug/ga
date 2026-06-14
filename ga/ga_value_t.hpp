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

// Relative tolerance for is_congruent. Congruence is almost always asked of two
// INDEPENDENTLY COMPUTED geometric quantities (results of wedge / dual / motor chains),
// whose difference has accumulated error far beyond a single rounding -- so the
// congruence test uses a geometric relative tolerance, not machine eps. Still ~12 orders
// of magnitude below any real geometric signal.
constexpr auto eps_congruent = 1.0e-12;

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

} // namespace detail

} // namespace hd::ga
