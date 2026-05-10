#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "core/ga_prdxpr_common.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <set>

////////////////////////////////////////////////////////////////////////////////
// Automatic GA Rule Generation System
//
// HOW TO ADD A NEW GEOMETRIC ALGEBRA (e.g., STA4D)
//
// This system automatically generates mathematically correct product rules
// for any geometric algebra from its metric signature and basis vectors.
//
// STEP 1: Create Header File (e.g., ga_prdxpr_sta4ds.hpp)
//   - Define multivector basis: const mvec_coeff mvSTA3d_basis = {"1", "g0", "g1", "g2",
//   "g3", ...};
//   - Add extern declarations for product rules:
//     extern const prd_rules gpr_sta4ds_rules;  // Geometric product
//     extern const prd_rules wdg_sta4ds_rules;  // Wedge product
//     extern const prd_rules dot_sta4ds_rules;  // Dot product
//   - Define coefficient arrays and complement rules as needed
//
// STEP 2: Create Config File (e.g., ga_prdxpr_sta4ds_config.cpp)
//   - Implement AlgebraConfig function:
//     AlgebraConfig get_sta4ds_algebra_config() {
//         return {
//             {"g1", "g2", "g3", "g4"},        // basis_vectors
//             {+1, +1, +1, -1},                // metric_signature
//             mvsta3ds_basis,                    // multivector_basis (from header)
//             one_str(),                         // scalar_name
//             "g"                              // basis_prefix
//         };
//     }
//   - Generate rules with static initialization:
//     static auto sta4ds_generated_rules =
//     generate_algebra_rules(get_sta4ds_algebra_config()); const prd_rules
//     gpr_sta4ds_rules = sta4ds_generated_rules.geometric_product; const prd_rules
//     wdg_sta4ds_rules = sta4ds_generated_rules.wedge_product; const prd_rules
//     dot_sta4ds_rules = sta4ds_generated_rules.dot_product;
//
// STEP 3: Add to Build System (CMakeLists.txt)
//   - Add config file to SOURCES: src_prdxpr/ga_prdxpr_sta4ds_config.cpp
//   - Add header to HEADERS: src_prdxpr/ga_prdxpr_sta4ds.hpp
//   - Include in test program if needed
//
// STEP 4: Integration (ga_prdxpr_generator.cpp)
//   - Include the new header: #include "ga_prdxpr_sta4ds.hpp"
//   - Add switch case for STA4D in dimension handlers
//   - Add product generation calls as needed
//
// SUPPORTED ALGEBRAS:
//   - Euclidean: G(p,0,0) - metric signature: all +1
//   - Minkowski: G(p,q,0) - metric signature: p positive, q negative
//   - Projective: G(p,0,r) - metric signature: p positive, r zeros
//   - General: G(p,q,r) - any combination
//
// LIMITATION - ORTHOGONAL BASIS ONLY:
//   Current implementation assumes ORTHOGONAL basis vectors (diagonal metric)
//   ✓ Supports: EGA, PGA, STA (all have orthogonal bases)
//   ✗ NOT supported: CGA (Conformal GA) - has non-orthogonal null vectors
//
//   For CGA support, AlgebraConfig would need extension:
//     - Add: std::vector<std::vector<int>> metric_matrix; // Full metric matrix
//     - Current: metric_signature only provides diagonal: g_{ii}
//     - CGA needs: off-diagonal elements g_{ij} for i≠j
//     - Example CGA3D: e₊·e₋ = -1 (off-diagonal), e₊·e₊ = 0 (diagonal)
//
// INDEXING SYSTEMS:
//   - 1-based: e1, e2, e3 (traditional GA notation)
//   - 0-based: g0, g1, g2, g3 (physics/spacetime notation)
//   - System automatically detects and handles both
//
// The generation system handles all mathematical details:
//   - Basis element multiplication with correct signs
//   - Canonical ordering (e.g., e31 not e13)
//   - Zero results for degenerate directions
//   - Grade-ordered output for readability
////////////////////////////////////////////////////////////////////////////////

struct AlgebraConfig {
    std::vector<std::string> basis_vectors; // {"e1", "e2", "e3"}
    std::vector<int> metric_signature;      // {+1, +1, +1} or {+1, -1, -1, -1}
                                            // LIMITATION: Diagonal metric only!
                                            // For CGA, need full metric matrix
    mvec_coeff
        multivector_basis; // {"1", "e1", "e2", "e12"} - complete basis in canonical order
    std::string scalar_name = one_str(); // Use consistent "1" from common
    std::string basis_prefix = "e";

    // FUTURE CGA EXTENSION:
    // std::vector<std::vector<int>> metric_matrix; // Full n×n metric matrix
    //   - Would replace metric_signature for non-orthogonal algebras
    //   - metric_matrix[i][j] = e_i · e_j (dot product of basis vectors)
    //   - For orthogonal algebras: metric_matrix[i][j] = 0 for i≠j
    //   - For CGA: non-zero off-diagonal elements (e.g., e₊·e₋ = -1)
};

struct ProductRules {
    mvec_coeff basis;
    prd_rules geometric_product;
    prd_rules wedge_product;
    prd_rules dot_product;
    // Complement rules - generated based on algebra dimensionality
    prd_rules r_cmpl;     // For even algebras (EGA2D, PGA3DP)
    prd_rules l_cmpl;     // For even algebras (EGA2D, PGA3DP)
    prd_rules complement; // For odd algebras (EGA3D, PGA2DP)

    // Dual rules - generated from complement rules + extended metric
    prd_rules r_dual; // For even Euclidean/Minkowski algebras (EGA2D, STA4D)
    prd_rules l_dual; // For even Euclidean/Minkowski algebras (EGA2D, STA4D)
    prd_rules dual;   // For odd Euclidean algebras (EGA3D)

    // PGA-specific dual rules (for projective algebras with degenerate dimensions)
    // Odd-dimensional PGA (PGA2DP): single bulk_dual and weight_dual
    prd_rules bulk_dual;   // For PGA2DP (odd-dimensional)
    prd_rules weight_dual; // For PGA2DP (odd-dimensional)

    // Even-dimensional PGA (PGA3DP): left/right bulk_dual and weight_dual
    prd_rules l_bulk_dual;   // For PGA3DP (even-dimensional)
    prd_rules r_bulk_dual;   // For PGA3DP (even-dimensional)
    prd_rules l_weight_dual; // For PGA3DP (even-dimensional)
    prd_rules r_weight_dual; // For PGA3DP (even-dimensional)
};

////////////////////////////////////////////////////////////////////////////////
// Helper Functions (declarations)
////////////////////////////////////////////////////////////////////////////////

std::vector<int> parse_indices(std::string const& basis_element,
                               std::string const& prefix = "e");
std::string indices_to_basis(std::vector<int> const& indices,
                             std::string const& prefix = "e");
std::pair<std::string, int> multiply_basis_elements(std::string const& a,
                                                    std::string const& b,
                                                    AlgebraConfig const& config);

// Look up the metric signature value of the basis vector whose name is
// "<basis_prefix><digit>" by finding its position in config.basis_vectors and
// indexing config.metric_signature at that position. Slot-keyed, not
// digit-keyed: reordering config.basis_vectors automatically reorders the
// signature lookup. Returns 0 if the digit does not correspond to a known
// basis vector. Replaces the older "digit ± 1-based offset" lookup, which
// silently assumed slot order matches digit order.
int vector_metric_value(int digit, AlgebraConfig const& config);

////////////////////////////////////////////////////////////////////////////////
// Extended Metric Calculation (declarations)
////////////////////////////////////////////////////////////////////////////////

std::vector<int> calculate_extended_metric_matrix(AlgebraConfig const& config);
std::vector<int> calculate_extended_metric(AlgebraConfig const& config);
std::vector<int> calculate_extended_metric_recursive(AlgebraConfig const& config,
                                                     prd_rules const& wedge_rules,
                                                     prd_rules const& dot_rules);
std::vector<int> calculate_regressive_extended_metric(AlgebraConfig const& config);
std::vector<int> calculate_regressive_extended_metric_matrix(AlgebraConfig const& config);

////////////////////////////////////////////////////////////////////////////////
// Main Generation Functions (declarations)
////////////////////////////////////////////////////////////////////////////////

mvec_coeff generate_basis(AlgebraConfig const& config);
prd_rules generate_ordered_rules(
    AlgebraConfig const& config, std::string const& operator_str,
    std::function<std::pair<std::string, int>(std::string const&, std::string const&,
                                              AlgebraConfig const&)>
        multiply_func);
prd_rules generate_geometric_product_rules(AlgebraConfig const& config);
std::pair<std::string, int> multiply_basis_elements_wedge(std::string const& a,
                                                          std::string const& b,
                                                          AlgebraConfig const& config);
std::pair<std::string, int> multiply_basis_elements_dot(std::string const& a,
                                                        std::string const& b,
                                                        AlgebraConfig const& config);
prd_rules generate_wedge_product_rules(AlgebraConfig const& config);
prd_rules generate_dot_product_rules(AlgebraConfig const& config);
prd_rules generate_r_cmpl_rules(AlgebraConfig const& config,
                                prd_rules const& wedge_rules);
prd_rules generate_l_cmpl_rules(AlgebraConfig const& config,
                                prd_rules const& wedge_rules);
prd_rules generate_cmpl_rules(AlgebraConfig const& config, prd_rules const& wedge_rules);
prd_rules generate_l_dual_rules(AlgebraConfig const& config,
                                prd_rules const& l_cmpl_rules);
prd_rules generate_r_dual_rules(AlgebraConfig const& config,
                                prd_rules const& r_cmpl_rules);
prd_rules generate_dual_rules(AlgebraConfig const& config,
                              prd_rules const& complement_rules);

// PGA-specific dual generation (bulk and weight duals for projective algebras)
prd_rules generate_bulk_dual_rules(AlgebraConfig const& config,
                                   prd_rules const& complement_rules);
prd_rules generate_weight_dual_rules(AlgebraConfig const& config,
                                     prd_rules const& complement_rules);
prd_rules generate_l_bulk_dual_rules(AlgebraConfig const& config,
                                     prd_rules const& l_cmpl_rules);
prd_rules generate_r_bulk_dual_rules(AlgebraConfig const& config,
                                     prd_rules const& r_cmpl_rules);
prd_rules generate_l_weight_dual_rules(AlgebraConfig const& config,
                                       prd_rules const& l_cmpl_rules);
prd_rules generate_r_weight_dual_rules(AlgebraConfig const& config,
                                       prd_rules const& r_cmpl_rules);

ProductRules generate_algebra_rules(AlgebraConfig const& config);

////////////////////////////////////////////////////////////////////////////////
// Validation and Testing Functions (declarations)
////////////////////////////////////////////////////////////////////////////////

bool validate_rules(prd_rules const& generated, prd_rules const& reference);
void print_rule_comparison(prd_rules const& generated, prd_rules const& reference,
                           std::string const& product_name);