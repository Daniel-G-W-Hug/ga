// Copyright 2024-2025, Daniel Hug. All rights reserved.

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
#include "ga_ega2d_test.cpp"
#include "ga_ega3d_test.cpp"
#include "ga_ega4d_test.cpp"


/////////////////////////////////////////////////////////////////////////////////////////
// Generic Euclideon geometric algebra (valid for ega independent of dimension)
/////////////////////////////////////////////////////////////////////////////////////////

// TEST_SUITE("EGA")
// {

// }
