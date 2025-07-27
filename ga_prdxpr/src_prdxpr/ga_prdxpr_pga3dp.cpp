// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_pga3dp.hpp"
#include "ga_prdxpr_pga3dp_config.hpp"
#include "ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - PGA3DP
//
// This file contains the automatic generation of product rules for PGA3DP.
// Rules are generated at program startup via static initialization.
//
// The algebra configuration is defined in ga_prdxpr_pga3dp_config.cpp
// The algebra interface is declared in ga_prdxpr_pga3dp.hpp
//
// This separation provides:
// - Clean separation between configuration (what) and generation (how)
// - Single initialization point for all PGA3DP product rules
// - Guaranteed consistency through algorithmic generation
////////////////////////////////////////////////////////////////////////////////

static auto pga3dp_generated_rules = generate_algebra_rules(get_pga3dp_algebra_config());

const prd_rules gpr_pga3dp_rules = pga3dp_generated_rules.geometric_product;
const prd_rules wdg_pga3dp_rules = pga3dp_generated_rules.wedge_product;
const prd_rules dot_pga3dp_rules = pga3dp_generated_rules.dot_product;