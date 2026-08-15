// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_cga3dc.hpp"
#include "algebras/ga_prdxpr_cga3dc_config.hpp"

#include "rules/ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - cga3dc
//
// Rules are generated at program startup via static initialization (staged, no
// cross-references — see the static-init discipline in the project docs).
//
// The geometric product rules are MULTI-TERM (null-pair metric) and generated
// via the diagonal-basis detour with the direct-construction cross-check; any
// route disagreement throws at startup, so a rule inconsistency cannot go
// unnoticed. (32x32 basis pairs, computed once per binary that links this TU.)
////////////////////////////////////////////////////////////////////////////////

// Stage 1: Generate all algebra rules in one go
static auto cga3dc_generated_rules = generate_algebra_rules(get_cga3dc_algebra_config());

// Stage 2: Initialize product rules (no dependencies)
const prd_rules_mt gpr_cga3dc_rules_mt = cga3dc_generated_rules.geometric_product_mt;
const prd_rules wdg_cga3dc_rules = cga3dc_generated_rules.wedge_product;
const prd_rules dot_cga3dc_rules = cga3dc_generated_rules.dot_product;

// Stage 3: Initialize complement rules (depend on wedge rules, but those are
// already generated; ODD dimension -> the generator fills the singular
// complement/dual/antidual fields, and left == right)
const prd_rules cmpl_cga3dc_rules = cga3dc_generated_rules.complement;

// Stage 4: Initialize dual rules (generated from complement rules + extended
// metric matrix; non-degenerate metric -> a single regular dual)
const prd_rules dual_cga3dc_rules = cga3dc_generated_rules.dual;

// Stage 5: Initialize antidual rules (negated dual; the anti-exomorphism
// identity is verified during generation)
const prd_rules antidual_cga3dc_rules = cga3dc_generated_rules.antidual;
