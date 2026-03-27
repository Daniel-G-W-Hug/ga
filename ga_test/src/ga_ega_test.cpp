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

#include <chrono>
#include <iostream>
#include <tuple>
#include <variant>
#include <vector>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_ega.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)

// Include dimension-specific test files
#include "ga_ega2d_test.hpp"
#include "ga_ega3d_test.hpp"
#include "ga_ega4d_test.hpp"


/////////////////////////////////////////////////////////////////////////////////////////
// Generic Euclidean geometric algebra (valid for ega independent of dimension)
/////////////////////////////////////////////////////////////////////////////////////////

// TEST_SUITE("EGA")
// {

// } // EGA tests
