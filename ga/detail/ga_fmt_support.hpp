#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "fmt/ga_fmt_core.hpp" // Core infrastructure
#include "fmt/ga_fmt_ega.hpp"  // EGA type formatters
#include "fmt/ga_fmt_pga.hpp"  // PGA type formatters

/////////////////////////////////////////////////////////////////////////////////////////
// Unified formatting support for all GA library types (backward compatibility)
// This file maintains compatibility with existing includes while providing
// the option to use more granular includes for better compilation performance
/////////////////////////////////////////////////////////////////////////////////////////

// For new code, consider using the more specific includes:
// - #include "detail/fmt/ga_fmt_ega.hpp"  // EGA types only
// - #include "detail/fmt/ga_fmt_pga.hpp"  // PGA types only
// - #include "detail/fmt/ga_fmt_core.hpp" // Core infrastructure only
//
// Usage examples:
//
// bivec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
//
// std::vector<BiVec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));