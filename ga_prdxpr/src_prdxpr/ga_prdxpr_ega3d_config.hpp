// Copyright 2024-2025, Daniel Hug. All rights reserved.

#pragma once

#include "ga_prdxpr_config_types.hpp"
#include "ga_prdxpr_ega3d.hpp"

// Algebra configuration function (used by rule generation)
AlgebraConfig get_ega3d_algebra_config();

namespace configurable {

// EGA3D algebra configuration
AlgebraData create_ega3d_algebra_data();

// EGA3D product configurations
ProductConfig get_ega3d_gpr_config();          // geometric product
ProductConfig get_ega3d_cmt_config();          // commutator product
ProductConfig get_ega3d_wdg_config();          // wedge product
ProductConfig get_ega3d_dot_config();          // inner product
ProductConfig get_ega3d_rwdg_config();         // regressive wedge
ProductConfig get_ega3d_lcontract_config();    // left contraction
ProductConfig get_ega3d_rcontract_config();    // right contraction
ProductConfig get_ega3d_lexpand_config();      // left expansion
ProductConfig get_ega3d_rexpand_config();      // right expansion
ProductConfig get_ega3d_sandwich_gpr_config(); // rotor geometric product (sandwich)

} // namespace configurable