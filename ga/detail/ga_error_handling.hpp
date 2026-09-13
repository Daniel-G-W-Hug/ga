#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <algorithm> // std::max
#include <array>     // std::array (componentwise comparison)
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string
#include <utility>   // std::index_sequence

/////////////////////////////////////////////////////////////////////////////////////////
// Centralized error handling utilities for GA library
// Provides consistent thresholds and error messages across all operations
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga::detail {

// Compile-time check for extended testing mode
// (enables conditional noexcept via pre-processor)
constexpr bool extended_testing_enabled() noexcept
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    return true;
#else
    return false;
#endif
}

// Consistent epsilon calculation for single type
template <typename T>
    requires(std::floating_point<T>)
constexpr T safe_epsilon() noexcept
{
    return T(5.0) * std::numeric_limits<T>::epsilon();
}

// Consistent epsilon calculation for mixed types
// hint: important to evaluate std::max<ctype>
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::common_type_t<T, U> safe_epsilon() noexcept
{
    using ctype = std::common_type_t<T, U>;
    constexpr ctype eps_t = ctype(std::numeric_limits<T>::epsilon());
    constexpr ctype eps_u = ctype(std::numeric_limits<U>::epsilon());
    return ctype(5.0) * std::max(eps_t, eps_u);
}

// Epsilon calculation for floating point with integral types
template <typename T, typename U>
    requires(std::floating_point<T> && std::integral<U>)
constexpr T safe_epsilon() noexcept
{
    return T(5.0) * std::numeric_limits<T>::epsilon();
}

template <typename T, typename U>
    requires(std::integral<T> && std::floating_point<U>)
constexpr U safe_epsilon() noexcept
{
    return U(5.0) * std::numeric_limits<U>::epsilon();
}

// Standardized error checking and throwing for division operations
template <typename T>
    requires(std::floating_point<T>)
inline void check_division_by_zero(T divisor, char const* operation_name = "division")
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(divisor) < safe_epsilon<T>()) {
        throw std::runtime_error(
            std::string("GA Error: ") + operation_name +
            " by value too small (near zero): " + std::to_string(divisor));
    }
#else
    // Suppress unused parameter warning when testing is disabled
    (void)divisor;
    (void)operation_name;
#endif
}

// Standardized error checking for mixed-type division operations
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline void check_division_by_zero(U divisor, char const* operation_name = "division")
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(divisor) < safe_epsilon<T, U>()) {
        throw std::runtime_error(
            std::string("GA Error: ") + operation_name +
            " by value too small (near zero): " + std::to_string(divisor));
    }
#else
    // Suppress unused parameter warning when testing is disabled
    (void)divisor;
    (void)operation_name;
#endif
}

// An absolute floor, for a DIMENSIONLESS divisor whose natural scale is one: a series
// coefficient, a dilation factor. Nothing is gained by judging such a quantity against
// the object it came from -- there is no object, and the number is already a ratio. For
// anything that scales with its operand use check_invertible / check_nonzero instead, and
// for the weight divisions check_unitization below.
template <typename T>
    requires(std::floating_point<T>)
inline void check_normalization(T magnitude, char const* object_type = "vector")
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (magnitude < safe_epsilon<T>()) {
        throw std::runtime_error(
            std::string("GA Error: ") + object_type +
            " magnitude too small for normalization: " + std::to_string(magnitude));
    }
#else
    // Suppress unused parameter warning when testing is disabled
    (void)magnitude;
    (void)object_type;
#endif
}

// The weight divisions: unitize() and the dehomogenizing accessors, which divide an
// object by its own weight. The floor is ABSOLUTE on purpose -- it flags a weight that is
// negligible on the scale the coordinates are expressed in, i.e. an object that is
// effectively ideal -- and no relative form separates "ideal" from "very far away", since
// the two differ exactly by that absolute comparison of weight against bulk.
template <typename T>
    requires(std::floating_point<T>)
inline void check_unitization(T weight_norm, char const* object_type = "multivector")
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (weight_norm < safe_epsilon<T>()) {
        throw std::runtime_error(
            std::string("GA Error: ") + object_type +
            " weight norm too small for unitization: " + std::to_string(weight_norm));
    }
#else
    // Suppress unused parameter warning when testing is disabled
    (void)weight_norm;
    (void)object_type;
#endif
}

// Standardized invertibility checks.
//
// WHY THESE ARE NOT ONE ABSOLUTE EPSILON. What makes a blade non-invertible is that it is
// DEGENERATE -- null in a Lorentzian metric, ideal in a degenerate one -- and that is a
// property of its direction, not of its size. An absolute bound on the divisor confuses
// the two: inv() of a vector of magnitude 1e-9 used to throw although its inverse (7e6)
// is perfectly representable, so the same geometry expressed in millimetres instead of
// metres failed, while a genuinely null blade at a large scale sailed through.
//
// So the divisor is judged against the SUM OF THE SQUARED COEFFICIENTS of the object it
// came from, which makes the test exactly scale-invariant -- both quantities are
// quadratic in the object. Where the metric is Euclidean the two are the SAME quantity,
// the ratio is 1 by construction, and the question collapses to "is the object zero":
// that is what check_nonzero below is for, and it is the honest spelling there.
template <typename T>
    requires(std::floating_point<T>)
inline void check_invertible(T divisor, T coeff_sq, char const* object_type = "blade")
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    // the negation also rejects NaN, which a plain <= would let through
    if (!(std::abs(divisor) > safe_epsilon<T>() * coeff_sq)) {
        throw std::runtime_error(std::string("GA Error: ") + object_type +
                                 " is degenerate (metric " + std::to_string(divisor) +
                                 " negligible against its own scale " +
                                 std::to_string(coeff_sq) + "), so it has no inverse");
    }
#else
    (void)divisor;
    (void)coeff_sq;
    (void)object_type;
#endif
}

// The Euclidean case of the above, and the case of a scalar: the divisor IS the object's
// own scale, so only a vanishing (or underflowed) one fails. A merely small object is
// invertible and must stay that way. The same spelling serves a precondition of the form
// "the zero element is not an object", where the quantity tested is likewise the object's
// own gauge and no relative form exists.
template <typename T>
    requires(std::floating_point<T>)
inline void check_nonzero(T divisor, char const* object_type = "blade")
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (!(std::abs(divisor) > T(0.0))) {
        throw std::runtime_error(std::string("GA Error: ") + object_type +
                                 " is zero -- no division possible");
    }
#else
    (void)divisor;
    (void)object_type;
#endif
}

// Standardized precondition check: the target of a projection or a reflection must be a
// BLADE, i.e. it must represent a subspace at all.
//
// Only 4d and up can violate this, and only at grade 2: a general bivector there is a sum
// of two blades and spans no plane (B ^ B != 0), so there is nothing to project onto and
// the expressions return an answer that is no projection -- not even idempotent. In 3d,
// and for a vector or a hyperplane target in any dimension, simplicity is automatic.
//
// The caller passes the ANSWER rather than the object: is_simple() lives in the algebra
// namespaces (hd::ga::pga, hd::ga::sta) while this helper lives in hd::ga::detail, and a
// call written here would bind to whichever algebra parsed first -- the same two-phase
// lookup trap that by_weight_sq() avoids the same way.
//
// Gated by its own switch so a consumer can keep the division guard while dropping this
// one (see ga/CMakeLists.txt): the check costs a wedge square, which the division guard
// does not.
inline void check_simple_target(bool target_is_simple, char const* operation_name)
{
#if defined(_HD_GA_EXTENDED_TEST_BLADE_TARGET)
    if (!target_is_simple) {
        throw std::runtime_error(
            std::string("GA Error: ") + operation_name +
            ": the target is not a blade (B ^ B != 0), so it represents no subspace");
    }
#else
    // Suppress unused parameter warnings when the check is disabled
    (void)target_is_simple;
    (void)operation_name;
#endif
}

// Componentwise equality of two value types, and the ONE place the comparison rule
// lives.
//
// Every operator== in the library delegates here, so the rule -- the tolerance, whether
// it is absolute or relative, and what it is measured against -- is stated once instead
// of once per type. A type definition therefore does not have to change when the rule
// does.
//
// The rule: two objects are equal when every component differs by less than
//
//     eps * max(|a|_inf, |b|_inf, 1)
//
// i.e. a RELATIVE tolerance of safe_epsilon measured against the larger of the two
// operands, with a floor of 1 so a component that should vanish is still compared
// against the tolerance itself rather than against nothing.
//
// In ulps: eps is ulp(1), so eps * scale lies between ulp(scale) and 2 * ulp(scale)
// depending on where scale sits inside its binade. The budget is therefore 5 to 10 ulps
// of the operand -- constant in relative terms, which is the property that matters, and
// the reason a measured threshold reads 5.5 at magnitude 1 and 9.8 at 1e3.
//
// The scale is taken over the WHOLE object, not per component. A point at (1e6, 0) must
// not require its second component to match to 1e-15 absolute: at that magnitude a
// computed zero carries the rounding of the coordinates it came from.
//
// This used to be an absolute tolerance, which is right near magnitude 1 and quietly
// wrong above it: the window stayed at 1.11e-15 while one ulp grew with the operands, so
// from about magnitude 1e3 the comparison was stricter than the floating-point grid and
// two values as close as doubles can be compared UNEQUAL. At an earth radius one ulp is
// ~9.3e-10, nearly a million times the old window. The relative form holds 5-10 ulps at
// every magnitude, and below 1 it reproduces the old behaviour exactly.
//
// is_close() is the same rule with a larger budget (eps_congruent, ~4500 ulps), for
// operands that were computed independently rather than by the same expression.
//
// The comparison is NOT exact and NOT transitive (a == b and b == c does not imply
// a == c), so it must not be used as an equivalence relation -- no sorting, no
// std::unique, no associative containers keyed on these types.
//
// Written as a fold rather than a loop on purpose: the loop short-circuits and costs a
// branch per component, while the fold compiles to the same branchless form the
// hand-written bodies had.
template <typename T, typename U, std::size_t N>
constexpr bool coeffs_equal(std::array<T, N> const& a, std::array<U, N> const& b)
{
    using ctype = std::common_type_t<T, U>;
    ctype scale = ctype(1);
    for (std::size_t i = 0; i < N; ++i) {
        scale = std::max({scale, ctype(std::abs(a[i])), ctype(std::abs(b[i]))});
    }
    ctype const tol = safe_epsilon<T, U>() * scale;
    // bitwise fold, not logical: && would short-circuit and put a branch between every
    // pair of components, where the hand-written bodies compared all of them branchlessly
    return [&]<std::size_t... I>(std::index_sequence<I...>) {
        return ((static_cast<unsigned>(std::abs(a[I] - b[I]) < tol) & ...)) != 0u;
    }(std::make_index_sequence<N>{});
}

// Does an object carry weight, i.e. is it NOT ideal (at infinity)?
//
// The question comes up in three places -- unitize()'s precondition (above), the
// conditional try_unitize(), and the projection normalization that divides a target's
// scale back out -- so it is answered once, here, and they cannot drift apart.
//
// The argument is the SQUARED weight norm (weight_nrm_sq), never the norm: unitize()
// divides by the norm and refuses below safe_epsilon, and the identical test on the
// square needs no square root,
//
//     weight_nrm > e   <=>   weight_nrm_sq > e * e
//
// A threshold placed on the square itself (wsq > e) would be seven orders of magnitude
// looser and would treat merely small objects as ideal.
template <typename T>
    requires(std::floating_point<T>)
constexpr bool has_weight(T weight_norm_sq) noexcept
{
    T const e = safe_epsilon<T>();
    return weight_norm_sq > e * e;
}


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