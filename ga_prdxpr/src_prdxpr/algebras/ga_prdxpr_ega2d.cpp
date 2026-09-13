// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_ega2d.hpp"
#include "algebras/ga_prdxpr_ega2d_config.hpp"
#include "rules/ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - EGA2D
//
// This file contains the automatic generation of product rules for EGA2D.
// Rules are generated at program startup via static initialization.
//
// The algebra configuration is defined in ga_prdxpr_ega2d_config.cpp
// The algebra interface is declared in ga_prdxpr_ega2d.hpp
//
// This separation provides:
// - Clean separation between configuration (what) and generation (how)
// - Single initialization point for all EGA2D product rules
// - Guaranteed consistency through algorithmic generation
////////////////////////////////////////////////////////////////////////////////

// Stage 1: Generate all algebra rules in one go
static auto ega2d_generated_rules = generate_algebra_rules(get_ega2d_algebra_config());

// Stage 2: Initialize product rules (no dependencies)
prd_rules const gpr_ega2d_rules = ega2d_generated_rules.geometric_product;
prd_rules const wdg_ega2d_rules = ega2d_generated_rules.wedge_product;
prd_rules const dot_ega2d_rules = ega2d_generated_rules.dot_product;

// Stage 3: Initialize complement rules (depend on wedge rules, but those are already
// generated)
prd_rules const l_cmpl_ega2d_rules = ega2d_generated_rules.l_cmpl;
prd_rules const r_cmpl_ega2d_rules = ega2d_generated_rules.r_cmpl;

// Stage 4: Initialize dual rules (generated from complement rules + extended metric)
prd_rules const l_dual_ega2d_rules = ega2d_generated_rules.l_dual;
prd_rules const r_dual_ega2d_rules = ega2d_generated_rules.r_dual;