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

namespace configurable {

// cga3dc algebra data creation
AlgebraData create_cga3dc_algebra_data();

// cga3dc product configurations (initial set; further products follow the
// established per-product workflow, mirroring the cga2dc rollout)
ProductConfig get_cga3dc_gpr_config(); // geometric product (multi-term rules)
ProductConfig get_cga3dc_cmt_config(); // commutator product (multi-term rules)
ProductConfig get_cga3dc_wdg_config(); // wedge product
ProductConfig get_cga3dc_dot_config(); // inner product

ProductConfig get_cga3dc_l_contract_config(); // left contraction
ProductConfig get_cga3dc_r_contract_config(); // right contraction

// Expansions are the regressive versions of the corresponding contractions.
ProductConfig get_cga3dc_l_expand_config(); // left expansion
ProductConfig get_cga3dc_r_expand_config(); // right expansion

ProductConfig get_cga3dc_rgpr_config(); // regressive geometric product (multi-term)
ProductConfig get_cga3dc_rcmt_config(); // regressive commutator product (multi-term)
ProductConfig get_cga3dc_rwdg_config(); // regressive wedge product
ProductConfig get_cga3dc_rdot_config(); // regressive inner product
ProductConfig get_cga3dc_sandwich_rgpr_config(); // regressive sandwich product

} // namespace configurable
