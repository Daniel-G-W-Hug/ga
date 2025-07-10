// Copyright 2024-2025, Daniel Hug. All rights reserved.

#pragma once

#include "ga_prdxpr_config_types.hpp"
#include "ga_prdxpr_pga3dp.hpp"

namespace configurable {

// PGA3DP algebra data creation
AlgebraData create_pga3dp_algebra_data();

// PGA3DP product configurations
ProductConfig get_pga3dp_gpr_config();                   // geometric product
ProductConfig get_pga3dp_cmt_config();                   // commutator product
ProductConfig get_pga3dp_wdg_config();                   // wedge product
ProductConfig get_pga3dp_dot_config();                   // inner product
ProductConfig get_pga3dp_rwdg_config();                  // regressive wedge product
ProductConfig get_pga3dp_rdot_config();                  // regressive inner product
ProductConfig get_pga3dp_rgpr_config();                  // regressive geometric product
ProductConfig get_pga3dp_right_bulk_contract_config();   // right bulk contraction
ProductConfig get_pga3dp_right_weight_contract_config(); // right weight contraction
ProductConfig get_pga3dp_left_bulk_contract_config();    // left bulk contraction
ProductConfig get_pga3dp_left_weight_contract_config();  // left weight contraction
ProductConfig get_pga3dp_right_bulk_expand_config();     // right bulk expansion
ProductConfig get_pga3dp_right_weight_expand_config();   // right weight expansion
ProductConfig get_pga3dp_left_bulk_expand_config();      // left bulk expansion
ProductConfig get_pga3dp_left_weight_expand_config();    // left weight expansion
ProductConfig get_pga3dp_sandwich_rgpr_config();         // regressive sandwich product

} // namespace configurable