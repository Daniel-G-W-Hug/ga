// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_sta.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::sta; // use specific operations of STA (Space-Time Algebra)

// Include dimension-specific test files (just prepared, not yet included)
#include "ga_sta4d_test.hpp"
