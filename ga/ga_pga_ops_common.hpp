#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_foundation.hpp" // ga library headers and infrastructure

#include "ga_value_t.hpp" // value_t, eps


namespace hd::ga::detail {

// remove the TARGET's scale from a projection-like result that is quadratic in it
//
// The projective projection, reflection and antiprojection expressions are built from
// products and duals alone, and those can only ADD homogeneity degrees: the weight dual
// is linear in the target and the meet brings the target in a second time. No such
// expression can be degree 0 in the target, so the scale has to be divided out --
// otherwise scaling a plane by 3 scales the projection onto it by 9. EGA divides by the
// same quantity inside inv(); the degenerate metric has no inv() for these blades, so
// the division is written out here.
//
// The caller passes the target's weight_nrm_sq rather than the target itself: the types
// live in one namespace and the operations in another, so a lookup inside this helper
// would bind to whichever algebra happened to be parsed first.
//
// An IDEAL target has no scale to remove and is left untouched: the caller gets the
// degenerate answer rather than a division by zero. "Ideal" is decided exactly as
// unitize() decides it -- weight_nrm > safe_epsilon -- but written on the SQUARE, which
// is the same test without a square root:
//
//     weight_nrm > e   <=>   weight_nrm_sq > e * e
//
// so try_unitize(), unitize() and this guard never disagree about which objects carry
// weight. Note this is NOT wsq > e: that would be a threshold on the norm seven orders
// of magnitude larger, and it would skip the division for targets that are merely small
// rather than ideal -- silently returning a result still scaled by their weight.
//
// wsq == 1 is a fast path only. weight_nrm_sq returns a plain scalar, so the comparison
// is an exact one and an object that is unit only to within rounding still takes the
// division -- which changes nothing, because it divides by 1 +/- an ulp.
//
// _HD_GA_FAST_PROJECTION_ON_UNITIZED removes the step entirely, for a caller who
// guarantees unitized targets. Define it on the TARGET (project-wide), never per file:
// this is a header-only library, so two translation units compiled with different
// settings of it would define the same inline function differently -- an ODR violation
// that the linker resolves silently, in favour of whichever it sees first.
template <typename R, typename S>
inline R by_weight_sq([[maybe_unused]] R r, [[maybe_unused]] S const& wsq)
{
#if !defined(_HD_GA_FAST_PROJECTION_ON_UNITIZED)
    if (has_weight(wsq) && (wsq != S(1.0))) r = r / wsq;
#endif
    return r;
}

} // namespace hd::ga::detail
