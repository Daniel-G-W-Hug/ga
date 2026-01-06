// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#ifndef GA_PRDXPR_DUAL_CALC_HPP
#define GA_PRDXPR_DUAL_CALC_HPP

#include "ga_prdxpr_rule_generator.hpp"
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// Dual Rule Calculation via Extended Metric Matrices
////////////////////////////////////////////////////////////////////////////////

// Calculate dual rules using extended metric matrix G
// Mathematical definition: dual(basis[i]) = Σⱼ G[i,j] * complement(basis[j])
//
// For diagonal metrics (EGA, PGA, STA):
//   dual(basis[i]) = G[i,i] * complement(basis[i])
//
// For non-diagonal metrics (future CGA):
//   Full matrix-vector multiplication required
//
// Returns: Map of {basis_element → dual_expression}
prd_rules calculate_dual_rules(AlgebraConfig const& config,
                               std::vector<int> const& G_data,
                               prd_rules const& complement_rules);

// Calculate regressive dual rules using regressive metric matrix Ḡ
// Mathematical definition: regressive_dual(basis[i]) = Σⱼ Ḡ[i,j] * complement(basis[j])
//
// For PGA algebras with diagonal Ḡ:
//   regressive_dual(basis[i]) = Ḡ[i,i] * complement(basis[i])
//
// Used for weight-part operations in projective geometric algebra
//
// Returns: Map of {basis_element → regressive_dual_expression}
prd_rules calculate_regressive_dual_rules(AlgebraConfig const& config,
                                          std::vector<int> const& G_reg_data,
                                          prd_rules const& complement_rules);

#endif // GA_PRDXPR_DUAL_CALC_HPP
