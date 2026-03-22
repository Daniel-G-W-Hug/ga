// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// Include PGA3DP application test files.
// IMPORTANT: doctest runs test suites ordered by filename (--order-by=file default).
// Any new header added here must be named so that it appears at the correct position
// in alphabetical order to preserve the intended test execution sequence:
//   ga_appl3dp_3d_cmp_test.hpp   ('3' < 'a') -> comparison tests first
//   ga_appl3dp_appl_test.hpp     ('a' < 'p') -> application tests second
//   ga_appl3dp_physics_test.hpp              -> physics tests last
#include "ga_appl3dp_3d_cmp_test.hpp"
#include "ga_appl3dp_appl_test.hpp"
#include "ga_appl3dp_physics_test.hpp"
