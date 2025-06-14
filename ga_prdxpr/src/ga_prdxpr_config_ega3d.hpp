#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config.hpp"

namespace configs {

// EGA3D Configuration and Products
extern const AlgebraConfig ega3d_config;
extern const std::vector<ProductDefinition> ega3d_products;

// EGA3D Product Cases
extern const std::vector<ProductCase> ega3d_geometric_cases;
extern const std::vector<ProductCase> ega3d_commutator_cases;
extern const std::vector<ProductCase> ega3d_wedge_cases;
extern const std::vector<ProductCase> ega3d_inner_cases;
extern const std::vector<ProductCase> ega3d_regressive_wedge_cases;
extern const std::vector<ProductCase> ega3d_left_contraction_cases;
extern const std::vector<ProductCase> ega3d_right_contraction_cases;
extern const std::vector<ProductCase> ega3d_left_expansion_cases;
extern const std::vector<ProductCase> ega3d_right_expansion_cases;

} // namespace configs