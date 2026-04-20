// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#pragma once

#include "ga_prdxpr_config_types.hpp"
#include "ga_prdxpr_ega3d.hpp"

// Algebra configuration function (used by rule generation)
AlgebraConfig get_ega3d_algebra_config();

namespace configurable {

// EGA3D algebra configuration
AlgebraData create_ega3d_algebra_data();

// EGA3D product configurations
ProductConfig get_ega3d_gpr_config();     // geometric product
ProductConfig get_ega3d_gpr_alt_config(); // geometric product (alternative def.)
ProductConfig get_ega3d_twdg1_config();   // transwedge product for k=1
ProductConfig get_ega3d_cmt_config();     // commutator product
ProductConfig get_ega3d_wdg_config();     // wedge product
ProductConfig get_ega3d_dot_config();     // inner product

ProductConfig get_ega3d_l_contract_config(); // left contraction
ProductConfig get_ega3d_r_contract_config(); // right contraction

ProductConfig get_ega3d_l_expand_config(); // left expansion
ProductConfig get_ega3d_r_expand_config(); // right expansion

ProductConfig get_ega3d_rgpr_config();     // regressive geometric product
ProductConfig get_ega3d_rgpr_alt_config(); // regressive geometric product (alt. def.)
ProductConfig get_ega3d_rtwdg1_config();   // regressive transwedge product for k=1
ProductConfig get_ega3d_rcmt_config();     // regressive commutator product
ProductConfig get_ega3d_rwdg_config();     // regressive wedge
ProductConfig get_ega3d_rdot_config();     // regressive inner product
ProductConfig get_ega3d_rgpr_alt_config(); // regressive geometric product (alt. def.)

ProductConfig get_ega3d_sandwich_gpr_config(); // rotor geometric product (sandwich)

} // namespace configurable