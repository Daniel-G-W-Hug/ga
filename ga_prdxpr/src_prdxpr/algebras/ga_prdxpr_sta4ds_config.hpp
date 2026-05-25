// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#pragma once

#include "algebras/ga_prdxpr_sta4ds.hpp"
#include "core/ga_prdxpr_config_types.hpp"

// Algebra configuration function (used by rule generation)
AlgebraConfig get_sta4ds_algebra_config();

namespace configurable {

// sta4ds algebra data creation
AlgebraData create_sta4ds_algebra_data();

// sta4ds product configurations
ProductConfig get_sta4ds_gpr_config(); // geometric product
ProductConfig get_sta4ds_cmt_config(); // commutator product
ProductConfig get_sta4ds_wdg_config(); // wedge product
ProductConfig get_sta4ds_dot_config(); // inner product

ProductConfig get_sta4ds_l_contract_config(); // left contraction
ProductConfig get_sta4ds_r_contract_config(); // right contraction

// Expansions are the regressive versions of the corresponding contractions.
ProductConfig get_sta4ds_l_expand_config(); // left expansion
ProductConfig get_sta4ds_r_expand_config(); // right expansion

ProductConfig get_sta4ds_rgpr_config(); // regressive geometric product
ProductConfig get_sta4ds_rcmt_config(); // regressive commutator product
ProductConfig get_sta4ds_rwdg_config(); // regressive wedge product
ProductConfig get_sta4ds_rdot_config(); // regressive inner product

ProductConfig get_sta4ds_sandwich_gpr_config(); // sandwich product

} // namespace configurable