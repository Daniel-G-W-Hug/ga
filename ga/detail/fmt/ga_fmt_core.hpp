#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "fmt/chrono.h" // chrono support
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples


#include "../type_t/ga_type_tags.hpp"
#include <type_traits>

/////////////////////////////////////////////////////////////////////////////////////////
// Core formatting infrastructure for GA library types
// Provides base formatter definitions and common utilities shared
// between EGA and PGA formatters
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga::detail {

// Helper function to suppress negative zero in formatting
// Converts -0.0 to +0.0 to prevent displaying "-0" in output
template <typename T> inline constexpr T suppress_negative_zero(T value)
{
    // Check if value is zero (both +0.0 and -0.0)
    // If so, return +0.0; otherwise return the original value
    return (value == T{0}) ? T{0} : value;
}

} // namespace hd::ga::detail
