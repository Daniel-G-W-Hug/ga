#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_prdxpr_config_types.hpp"
#include "ga_prdxpr_ega2d.hpp"

// Algebra configuration function (used by rule generation)
AlgebraConfig get_ega2d_algebra_config();

namespace configurable {

// Create EGA2D algebra data using existing mathematical definitions (unchanged)
AlgebraData create_ega2d_algebra_data();

ProductConfig get_ega2d_gpr_config();     // geometric product
ProductConfig get_ega2d_gpr_alt_config(); // geometric product (alternative def.)
ProductConfig get_ega2d_twdg1_config();   // transwedge product for k=1
ProductConfig get_ega2d_cmt_config();     // commutator product
ProductConfig get_ega2d_wdg_config();     // wedge product
ProductConfig get_ega2d_dot_config();     // inner product
ProductConfig get_ega2d_dot_alt_config(); // inner product (alternative definition)

ProductConfig get_ega2d_l_contract_config(); // left contraction
ProductConfig get_ega2d_r_contract_config(); // right contraction

ProductConfig get_ega2d_l_expand_config(); // left expansion
ProductConfig get_ega2d_r_expand_config(); // right expansion

ProductConfig get_ega2d_rgpr_config();     // regressive geometric product
ProductConfig get_ega2d_rgpr_alt_config(); // regressive geometric product (alt. def.)
ProductConfig get_ega2d_rtwdg1_config();   // regressive transwedge product for k=1
ProductConfig get_ega2d_rcmt_config();     // regressive commutator product
ProductConfig get_ega2d_rwdg_config();     // regressive wedge product
ProductConfig get_ega2d_rdot_config();     // regressive inner product

ProductConfig get_ega2d_sandwich_gpr_config(); // rotor geometric product (sandwich)

} // namespace configurable