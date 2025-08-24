// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_ega3d.hpp"
#include "ga_prdxpr_ega3d_config.hpp"
#include "ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - EGA3D
//
// This file contains the automatic generation of product rules for EGA3D.
// Rules are generated at program startup via static initialization.
//
// The algebra configuration is defined in ga_prdxpr_ega3d_config.cpp
// The algebra interface is declared in ga_prdxpr_ega3d.hpp
//
// This separation provides:
// - Clean separation between configuration (what) and generation (how)
// - Single initialization point for all EGA3D product rules
// - Guaranteed consistency through algorithmic generation
////////////////////////////////////////////////////////////////////////////////

// Stage 1: Generate all algebra rules in one go
static auto ega3d_generated_rules = generate_algebra_rules(get_ega3d_algebra_config());

// Stage 2: Initialize product rules (no dependencies)
const prd_rules gpr_ega3d_rules = ega3d_generated_rules.geometric_product;
const prd_rules wdg_ega3d_rules = ega3d_generated_rules.wedge_product;
const prd_rules dot_ega3d_rules = ega3d_generated_rules.dot_product;

// Stage 3: Initialize complement rules (depend on wedge rules, but those are already
// generated)
const prd_rules cmpl_ega3d_rules = ega3d_generated_rules.complement;

// Stage 4: TODO - Initialize dual rules (will depend on complement rules + extended
// metric) For now, dual rules reference complement rules manually (metric is identity in
// EGA3D) Future: Generate dual rules from complement rules + extended metric information