// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#pragma once

#include "core/ga_prdxpr_config_types.hpp"
#include "algebras/ga_prdxpr_pga3dp.hpp"

// Algebra configuration function (used by rule generation)
AlgebraConfig get_pga3dp_algebra_config();

namespace configurable {

// PGA3DP algebra data creation
AlgebraData create_pga3dp_algebra_data();

// PGA3DP product configurations
ProductConfig get_pga3dp_gpr_config();     // geometric product
ProductConfig get_pga3dp_gpr_alt_config(); // geometric product (alternative def.)
ProductConfig get_pga3dp_twdg1_config();   // transwedge product for k=1
ProductConfig get_pga3dp_cmt_config();     // commutator product
ProductConfig get_pga3dp_wdg_config();     // wedge product
ProductConfig get_pga3dp_dot_config();     // inner product

ProductConfig get_pga3dp_l_bulk_contract_config();   // left bulk contraction
ProductConfig get_pga3dp_r_bulk_contract_config();   // right bulk contraction
ProductConfig get_pga3dp_l_weight_contract_config(); // left weight contraction
ProductConfig get_pga3dp_r_weight_contract_config(); // right weight contraction

// Bulk expansions are the regressive versions of the corresponding weight contractions.
// Weight expansions are the regressive versions of the corresponding bulk contractions.
//
ProductConfig get_pga3dp_l_bulk_expand_config();   // left bulk expansion
ProductConfig get_pga3dp_r_bulk_expand_config();   // right bulk expansion
ProductConfig get_pga3dp_l_weight_expand_config(); // left weight expansion
ProductConfig get_pga3dp_r_weight_expand_config(); // right weight expansion

ProductConfig get_pga3dp_rgpr_config();     // regressive geometric product
ProductConfig get_pga3dp_rgpr_alt_config(); // regressive geometric product (alt. def.)
ProductConfig get_pga3dp_rtwdg1_config();   // regressive transwedge product for k=1
ProductConfig get_pga3dp_rcmt_config();     // regressive commutator product
ProductConfig get_pga3dp_rwdg_config();     // regressive wedge product
ProductConfig get_pga3dp_rdot_config();     // regressive inner product

ProductConfig get_pga3dp_sandwich_rgpr_config(); // regressive sandwich product

} // namespace configurable