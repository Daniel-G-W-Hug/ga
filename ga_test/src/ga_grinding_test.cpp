// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// Wafer-grinding application test bundle (PGA3DP). Extracted from the shared
// ga_appl3dp_* test files so the grinding case is one self-contained target.
// doctest orders suites by filename (--order-by=file): _appl_ before _mechanics_.
#include "ga_grinding_appl_test.hpp"
#include "ga_grinding_mechanics_test.hpp"
