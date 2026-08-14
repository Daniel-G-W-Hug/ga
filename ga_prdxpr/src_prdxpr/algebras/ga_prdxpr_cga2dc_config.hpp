#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "core/ga_prdxpr_config_types.hpp"
#include "rules/ga_prdxpr_rule_generator.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga2dc algebra configuration
//
// Provides the AlgebraConfig for automatic rule generation (basis, the full
// non-orthogonal metric matrix of the null-vector basis, and the exact
// diagonalization used by the product computation) plus the generator-facing
// AlgebraData and ProductConfig builders. The generated wedge/geometric
// product tables were reviewed against the reference literature before the
// expression cases were added (see reference_output/).
////////////////////////////////////////////////////////////////////////////////

AlgebraConfig get_cga2dc_algebra_config();

namespace configurable {

// cga2dc algebra data creation
AlgebraData create_cga2dc_algebra_data();

// cga2dc product configurations (initial set: the literature-reviewed
// products; further products follow the established per-product workflow)
ProductConfig get_cga2dc_gpr_config(); // geometric product (multi-term rules)
ProductConfig get_cga2dc_cmt_config(); // commutator product (multi-term rules)
ProductConfig get_cga2dc_wdg_config(); // wedge product
ProductConfig get_cga2dc_dot_config(); // inner product

ProductConfig get_cga2dc_l_contract_config(); // left contraction
ProductConfig get_cga2dc_r_contract_config(); // right contraction

// Expansions are the regressive versions of the corresponding contractions.
ProductConfig get_cga2dc_l_expand_config(); // left expansion
ProductConfig get_cga2dc_r_expand_config(); // right expansion

ProductConfig get_cga2dc_rgpr_config(); // regressive geometric product (multi-term)
ProductConfig get_cga2dc_rcmt_config(); // regressive commutator product (multi-term)
ProductConfig get_cga2dc_rwdg_config(); // regressive wedge product
ProductConfig get_cga2dc_rdot_config(); // regressive inner product

ProductConfig get_cga2dc_sandwich_rgpr_config(); // regressive sandwich product

} // namespace configurable
