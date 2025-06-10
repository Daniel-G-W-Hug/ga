#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "../type_t/ga_type_tags.hpp"
#include <type_traits>

/////////////////////////////////////////////////////////////////////////////////////////
// Core formatting infrastructure for GA library types
// Provides base formatter definitions and common utilities shared
// between EGA and PGA formatters
/////////////////////////////////////////////////////////////////////////////////////////
