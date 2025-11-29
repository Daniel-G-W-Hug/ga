// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_sta4d.hpp"
#include "ga_prdxpr_sta4d_config.hpp"

#include "ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - sta4d
//
// This file contains the automatic generation of product rules for sta4d.
// Rules are generated at program startup via static initialization.
//
// The algebra configuration is defined in ga_prdxpr_sta4d_config.cpp
// The algebra interface is declared in ga_prdxpr_sta4d.hpp
//
// This separation provides:
// - Clean separation between configuration (what) and generation (how)
// - Single initialization point for all sta4d product rules
// - Guaranteed consistency through algorithmic generation
////////////////////////////////////////////////////////////////////////////////

// Stage 1: Generate all algebra rules in one go
static auto sta4d_generated_rules = generate_algebra_rules(get_sta4d_algebra_config());

// Stage 2: Initialize product rules (no dependencies)
const prd_rules gpr_sta4d_rules = sta4d_generated_rules.geometric_product;
const prd_rules wdg_sta4d_rules = sta4d_generated_rules.wedge_product;
const prd_rules dot_sta4d_rules = sta4d_generated_rules.dot_product;

// Stage 3: Initialize complement rules (depend on wedge rules, but those are already
// generated)
const prd_rules lcmpl_sta4d_rules = sta4d_generated_rules.left_complement;
const prd_rules rcmpl_sta4d_rules = sta4d_generated_rules.right_complement;

// Stage 4: Initialize dual rules (generated from complement rules + extended metric)
const prd_rules left_dual_sta4d_rules = sta4d_generated_rules.left_dual;
const prd_rules right_dual_sta4d_rules = sta4d_generated_rules.right_dual;