#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// Foundation infrastructure for GA library
// Common dependencies used across all GA operations and types
/////////////////////////////////////////////////////////////////////////////////////////

// Standard library headers commonly used throughout GA operations
#include <algorithm> // std::max, sorting, algorithms
#include <cmath>     // mathematical functions
#include <concepts>  // C++20 concepts for type constraints
#include <iostream>  // I/O operations
#include <limits>    // std::numeric_limits
#include <numbers>   // std::numbers::pi, etc.
#include <stdexcept> // std::runtime_error
#include <string>    // std::string

// Core GA infrastructure
#include "../ga_algebra.hpp"     // Algebra definitions G(P,N,Z)
#include "../ga_value_t.hpp"     // Default scalar type (float/double)
#include "ga_error_handling.hpp" // Centralized error handling utilities

/////////////////////////////////////////////////////////////////////////////////////////
// This header provides the common foundation that all GA operations depend on.
// By including this single header, operations files get all the basic
// infrastructure they need without having to repeat these includes.
//
// Usage in operations files:
// #include "detail/ga_foundation.hpp"
//
// Now you have access to:
// - All common standard library facilities
// - value_t type and eps constant
// - algebra<P,N,Z> template
// - Error checking functions (check_division_by_zero, etc.)
/////////////////////////////////////////////////////////////////////////////////////////