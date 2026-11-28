#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Tolerances for suites that also run with value_t = float (_HD_GA_VALUE_T_FLOAT, see
// ga/ga_value_t.hpp).
//
// Most of the suites' absolute tolerances (1.0e-12, 1.0e-14, ...) were written for
// double: a budget of so many roundings, stated in double's units. scaled_tol(t) keeps
// the budget and changes the unit -- it scales t by the ratio of the machine epsilons,
//
//     scaled_tol(t) = t * epsilon(value_t) / epsilon(double)
//
// so 1.0e-12 stays 1.0e-12 in a double build (the ratio is exactly 1, the comparison is
// the one it always was) and becomes 5.4e-4 in a float build.
//
// It is meant for ROUNDING budgets only. A bound on a TRUNCATION error -- a finite
// difference, whose float round-off at a step chosen for double swamps the truncation
// -- does not rescale this way; per_precision(d, f) states it once per precision, the
// float value measured. And a convergence order or a quantity pinned to double's
// resolution has no float counterpart at all; such checks are guarded with
// value_t_is_double.

#include <limits>      // std::numeric_limits
#include <type_traits> // std::is_same_v

#include "ga/ga_value_t.hpp" // value_t

inline constexpr bool value_t_is_double = std::is_same_v<hd::ga::value_t, double>;

constexpr hd::ga::value_t scaled_tol(double t)
{
    using hd::ga::value_t;
    return value_t(t * (double(std::numeric_limits<value_t>::epsilon()) /
                        std::numeric_limits<double>::epsilon()));
}

constexpr hd::ga::value_t per_precision(double for_double, double for_float)
{
    using hd::ga::value_t;
    return value_t_is_double ? value_t(for_double) : value_t(for_float);
}

static_assert(!value_t_is_double || double(scaled_tol(1.0e-12)) == 1.0e-12,
              "scaled_tol must be the identity in a double build");
