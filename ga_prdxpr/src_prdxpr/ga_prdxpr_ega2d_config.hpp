#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config_types.hpp"
#include "ga_prdxpr_ega2d.hpp"

// Algebra configuration function (used by rule generation)
AlgebraConfig get_ega2d_algebra_config();

namespace configurable {

// Create EGA2D algebra data using existing mathematical definitions (unchanged)
AlgebraData create_ega2d_algebra_data();

ProductConfig get_ega2d_gpr_config();          // geometric product
ProductConfig get_ega2d_cmt_config();          // commutator product
ProductConfig get_ega2d_wdg_config();          // wedge product
ProductConfig get_ega2d_dot_config();          // inner product
ProductConfig get_ega2d_dot_alt_config();      // inner product (alternative definition)
ProductConfig get_ega2d_rwdg_config();         // regressive wedge
ProductConfig get_ega2d_lcontract_config();    // left contraction
ProductConfig get_ega2d_rcontract_config();    // right contraction
ProductConfig get_ega2d_lexpand_config();      // left expansion
ProductConfig get_ega2d_rexpand_config();      // right expansion
ProductConfig get_ega2d_sandwich_gpr_config(); // rotor geometric product (sandwich)

} // namespace configurable