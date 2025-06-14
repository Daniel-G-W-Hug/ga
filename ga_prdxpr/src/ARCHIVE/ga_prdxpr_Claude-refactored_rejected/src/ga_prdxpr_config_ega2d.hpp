#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config.hpp"

namespace configs {

// EGA2D Configuration and Products
extern const AlgebraConfig ega2d_config;
extern const std::vector<ProductDefinition> ega2d_products;

// EGA2D Product Cases
extern const std::vector<ProductCase> ega2d_geometric_cases;
extern const std::vector<ProductCase> ega2d_commutator_cases;
extern const std::vector<ProductCase> ega2d_wedge_cases;
extern const std::vector<ProductCase> ega2d_inner_cases;
extern const std::vector<ProductCase> ega2d_regressive_wedge_cases;
extern const std::vector<ProductCase> ega2d_left_contraction_cases;
extern const std::vector<ProductCase> ega2d_right_contraction_cases;

} // namespace configs