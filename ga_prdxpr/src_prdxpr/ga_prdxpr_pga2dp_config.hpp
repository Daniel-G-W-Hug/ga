// Copyright 2024-2025, Daniel Hug. All rights reserved.

#pragma once

#include "ga_prdxpr_config_types.hpp"
#include "ga_prdxpr_pga2dp.hpp"

namespace configurable {

// PGA2DP Algebra configuration
AlgebraData create_pga2dp_algebra_data();

// PGA2DP Product configurations
ProductConfig get_pga2dp_gpr_config();              // geometric product
ProductConfig get_pga2dp_cmt_config();              // commutator product
ProductConfig get_pga2dp_wdg_config();              // wedge product
ProductConfig get_pga2dp_dot_config();              // inner product
ProductConfig get_pga2dp_rwdg_config();             // regressive wedge product
ProductConfig get_pga2dp_rdot_config();             // regressive inner product
ProductConfig get_pga2dp_rgpr_config();             // regressive geometric product
ProductConfig get_pga2dp_rbulk_contract_config();   // right bulk contraction
ProductConfig get_pga2dp_rweight_contract_config(); // right weight contraction
ProductConfig get_pga2dp_lbulk_contract_config();   // left bulk contraction
ProductConfig get_pga2dp_lweight_contract_config(); // left weight contraction
ProductConfig get_pga2dp_rbulk_expand_config();     // right bulk expansion
ProductConfig get_pga2dp_rweight_expand_config();   // right weight expansion
ProductConfig get_pga2dp_lbulk_expand_config();     // left bulk expansion
ProductConfig get_pga2dp_lweight_expand_config();   // left weight expansion
ProductConfig get_pga2dp_sandwich_rgpr_config();    // regressive sandwich product

} // namespace configurable