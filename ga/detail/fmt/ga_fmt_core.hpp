#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "../type_t/ga_type_tags.hpp"
#include <type_traits>

////////////////////////////////////////////////////////////////////////////////
// Core formatting infrastructure for GA library types
// Provides base formatter definitions and common utilities
////////////////////////////////////////////////////////////////////////////////

// This file contains the core formatter infrastructure that is shared
// between EGA and PGA formatters, including basic template definitions
// and common formatting utilities.

// Usage example for array types:
//
// std::vector<BiVec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// BiVec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));