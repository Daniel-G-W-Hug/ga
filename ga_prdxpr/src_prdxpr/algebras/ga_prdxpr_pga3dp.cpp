// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_pga3dp.hpp"
#include "algebras/ga_prdxpr_pga3dp_config.hpp"
#include "rules/ga_prdxpr_rule_generator.hpp"

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

// Stage 1: Generate all algebra rules in one go
static auto pga3dp_generated_rules = generate_algebra_rules(get_pga3dp_algebra_config());

// Stage 2: Initialize product rules (no dependencies)
const prd_rules gpr_pga3dp_rules = pga3dp_generated_rules.geometric_product;
const prd_rules wdg_pga3dp_rules = pga3dp_generated_rules.wedge_product;
const prd_rules dot_pga3dp_rules = pga3dp_generated_rules.dot_product;

// Stage 3: Initialize complement rules (depend on wedge rules, but those are already
// generated)
const prd_rules l_cmpl_pga3dp_rules = pga3dp_generated_rules.l_cmpl;
const prd_rules r_cmpl_pga3dp_rules = pga3dp_generated_rules.r_cmpl;

// Stage 4: Initialize dual rules (generated from complement rules + extended metric)
// PGA3DP has l_bulk_dual, r_bulk_dual, l_weight_dual, r_weight_dual
const prd_rules l_bulk_dual_pga3dp_rules = pga3dp_generated_rules.l_bulk_dual;
const prd_rules r_bulk_dual_pga3dp_rules = pga3dp_generated_rules.r_bulk_dual;
const prd_rules l_weight_dual_pga3dp_rules = pga3dp_generated_rules.l_weight_dual;
const prd_rules r_weight_dual_pga3dp_rules = pga3dp_generated_rules.r_weight_dual;