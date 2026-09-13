// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_sta4ds.hpp"
#include "algebras/ga_prdxpr_sta4ds_config.hpp"

#include "rules/ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - sta4ds
//
// This file contains the automatic generation of product rules for sta4ds.
// Rules are generated at program startup via static initialization.
//
// The algebra configuration is defined in ga_prdxpr_sta4ds_config.cpp
// The algebra interface is declared in ga_prdxpr_sta4ds.hpp
//
// This separation provides:
// - Clean separation between configuration (what) and generation (how)
// - Single initialization point for all sta4ds product rules
// - Guaranteed consistency through algorithmic generation
////////////////////////////////////////////////////////////////////////////////

// Stage 1: Generate all algebra rules in one go
static auto sta4ds_generated_rules = generate_algebra_rules(get_sta4ds_algebra_config());

// Stage 2: Initialize product rules (no dependencies)
prd_rules const gpr_sta4ds_rules = sta4ds_generated_rules.geometric_product;
prd_rules const wdg_sta4ds_rules = sta4ds_generated_rules.wedge_product;
prd_rules const dot_sta4ds_rules = sta4ds_generated_rules.dot_product;

// Stage 3: Initialize complement rules (depend on wedge rules, but those are already
// generated)
prd_rules const l_cmpl_sta4ds_rules = sta4ds_generated_rules.l_cmpl;
prd_rules const r_cmpl_sta4ds_rules = sta4ds_generated_rules.r_cmpl;

// Stage 4: Initialize dual rules (generated from complement rules + extended metric)
prd_rules const l_dual_sta4ds_rules = sta4ds_generated_rules.l_dual;
prd_rules const r_dual_sta4ds_rules = sta4ds_generated_rules.r_dual;