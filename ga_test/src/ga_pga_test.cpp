// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// WORKAROUND: doctest falls back to __LINE__ for Clang >= 22, causing symbol
// collisions when multiple test headers share the same line number for TEST_CASE.
// Remove this define once doctest fixes the DOCTEST_CLANG >= 22 branch in:
//   /opt/homebrew/include/doctest/doctest.h, lines ~556-564
// Verify by checking: if (DOCTEST_CLANG >= DOCTEST_COMPILER(22,0,0)) is gone or
// uses __COUNTER__ again. Track at: https://github.com/doctest/doctest/issues
#define DOCTEST_COUNTER __COUNTER__
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)

// Include dimension-specific test files
#include "ga_pga2dp_test.hpp"
#include "ga_pga3dp_test.hpp"
