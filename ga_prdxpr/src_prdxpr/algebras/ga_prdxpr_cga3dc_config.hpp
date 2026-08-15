#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "core/ga_prdxpr_config_types.hpp"
#include "rules/ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga3dc algebra configuration
//
// Provides the AlgebraConfig for automatic rule generation (basis, the full
// non-orthogonal metric matrix of the null-vector basis, and the exact
// diagonalization used by the product computation). The generated complement,
// metric, wedge and geometric product tables were reviewed against the
// reference literature (see reference_output/). The generator-facing
// AlgebraData and ProductConfig builders follow with the main-generator
// registration (mirroring the cga2dc rollout).
////////////////////////////////////////////////////////////////////////////////

AlgebraConfig get_cga3dc_algebra_config();
