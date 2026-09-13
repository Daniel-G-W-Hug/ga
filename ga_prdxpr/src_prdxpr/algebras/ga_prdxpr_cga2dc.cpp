// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_cga2dc.hpp"
#include "algebras/ga_prdxpr_cga2dc_config.hpp"

#include "rules/ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - cga2dc
//
// Rules are generated at program startup via static initialization (staged, no
// cross-references — see the static-init discipline in the project docs).
//
// The geometric product rules are MULTI-TERM (null-pair metric) and generated
// via the diagonal-basis detour with the direct-construction cross-check; any
// route disagreement throws at startup, so a rule inconsistency cannot go
// unnoticed.
////////////////////////////////////////////////////////////////////////////////

// Stage 1: Generate all algebra rules in one go
static auto cga2dc_generated_rules = generate_algebra_rules(get_cga2dc_algebra_config());

// Stage 2: Initialize product rules (no dependencies)
prd_rules_mt const gpr_cga2dc_rules_mt = cga2dc_generated_rules.geometric_product_mt;
prd_rules const wdg_cga2dc_rules = cga2dc_generated_rules.wedge_product;
prd_rules const dot_cga2dc_rules = cga2dc_generated_rules.dot_product;

// Stage 3: Initialize complement rules (depend on wedge rules, but those are
// already generated)
prd_rules const l_cmpl_cga2dc_rules = cga2dc_generated_rules.l_cmpl;
prd_rules const r_cmpl_cga2dc_rules = cga2dc_generated_rules.r_cmpl;

// Stage 4: Initialize dual rules (generated from complement rules + extended
// metric matrix; non-degenerate metric -> regular l/r duals)
prd_rules const l_dual_cga2dc_rules = cga2dc_generated_rules.l_dual;
prd_rules const r_dual_cga2dc_rules = cga2dc_generated_rules.r_dual;

// Stage 5: Initialize antidual rules (negated duals; 𝔾 = -G verified during
// generation)
prd_rules const l_antidual_cga2dc_rules = cga2dc_generated_rules.l_antidual;
prd_rules const r_antidual_cga2dc_rules = cga2dc_generated_rules.r_antidual;
