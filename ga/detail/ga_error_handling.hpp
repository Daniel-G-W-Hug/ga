#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <algorithm> // std::max
#include <concepts>  // std::floating_point
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

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
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::common_type_t<T, U> safe_epsilon() noexcept
{
    using ctype = std::common_type_t<T, U>;
    return ctype(5.0) *
           std::max(std::numeric_limits<T>::epsilon(), std::numeric_limits<U>::epsilon());
}

// Standardized error checking and throwing for division operations
template <typename T>
    requires(std::floating_point<T>)
inline void check_division_by_zero(T divisor, const char* operation_name = "division")
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
inline void check_division_by_zero(U divisor, const char* operation_name = "division")
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

// Standardized error checking for normalization operations
template <typename T>
    requires(std::floating_point<T>)
inline void check_normalization(T magnitude, const char* object_type = "vector")
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

// Standardized error checking for unitization operations
template <typename T>
    requires(std::floating_point<T>)
inline void check_unitization(T weight_norm, const char* object_type = "multivector")
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

} // namespace hd::ga::detail