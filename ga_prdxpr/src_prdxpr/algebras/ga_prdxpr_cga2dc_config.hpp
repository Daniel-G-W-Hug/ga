#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "rules/ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga2dc algebra configuration
//
// Provides the AlgebraConfig for automatic rule generation: basis, the full
// (non-orthogonal) metric matrix of the null-vector basis, and the exact
// diagonalization used by the product computation (see the rule generator).
//
// ProductConfig builders (expression cases for the generator) follow in a
// later step, after the generated wedge/geometric product tables have been
// reviewed against the reference literature.
////////////////////////////////////////////////////////////////////////////////

AlgebraConfig get_cga2dc_algebra_config();
