#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <set>

////////////////////////////////////////////////////////////////////////////////
// Automatic GA Rule Generation System
//
// HOW TO ADD A NEW GEOMETRIC ALGEBRA (e.g., STA3D)
//
// This system automatically generates mathematically correct product rules
// for any geometric algebra from its metric signature and basis vectors.
//
// STEP 1: Create Header File (e.g., ga_prdxpr_sta3d.hpp)
//   - Define multivector basis: const mvec_coeff mvSTA3d_basis = {"1", "g0", "g1", "g2",
//   "g3", ...};
//   - Add extern declarations for product rules:
//     extern const prd_rules gpr_sta3d_rules;  // Geometric product
//     extern const prd_rules wdg_sta3d_rules;  // Wedge product
//     extern const prd_rules dot_sta3d_rules;  // Dot product
//   - Define coefficient arrays and complement rules as needed
//
// STEP 2: Create Config File (e.g., ga_prdxpr_sta3d_config.cpp)
//   - Implement AlgebraConfig function:
//     AlgebraConfig get_sta3d_algebra_config() {
//         return {
//             {"g0", "g1", "g2", "g3"},        // basis_vectors
//             {+1, -1, -1, -1},                // metric_signature
//             mvSTA3d_basis,                    // multivector_basis (from header)
//             one_str(),                         // scalar_name
//             "g"                              // basis_prefix
//         };
//     }
//   - Generate rules with static initialization:
//     static auto sta3d_generated_rules =
//     generate_algebra_rules(get_sta3d_algebra_config()); const prd_rules gpr_sta3d_rules
//     = sta3d_generated_rules.geometric_product; const prd_rules wdg_sta3d_rules =
//     sta3d_generated_rules.wedge_product; const prd_rules dot_sta3d_rules =
//     sta3d_generated_rules.dot_product;
//
// STEP 3: Add to Build System (CMakeLists.txt)
//   - Add config file to SOURCES: src_prdxpr/ga_prdxpr_sta3d_config.cpp
//   - Add header to HEADERS: src_prdxpr/ga_prdxpr_sta3d.hpp
//   - Include in test program if needed
//
// STEP 4: Integration (ga_prdxpr_generator.cpp)
//   - Include the new header: #include "ga_prdxpr_sta3d.hpp"
//   - Add switch case for STA3D in dimension handlers
//   - Add product generation calls as needed
//
// SUPPORTED ALGEBRAS:
//   - Euclidean: G(p,0,0) - metric signature: all +1
//   - Minkowski: G(p,q,0) - metric signature: p positive, q negative
//   - Projective: G(p,0,r) - metric signature: p positive, r zeros
//   - General: G(p,q,r) - any combination
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
    mvec_coeff
        multivector_basis; // {"1", "e1", "e2", "e12"} - complete basis in canonical order
    std::string scalar_name = one_str(); // Use consistent "1" from common
    std::string basis_prefix = "e";
};

struct ProductRules {
    mvec_coeff basis;
    prd_rules geometric_product;
    prd_rules wedge_product;
    prd_rules dot_product;
    // Complement rules - generated based on algebra dimensionality
    prd_rules right_complement; // For even algebras (EGA2D, PGA3DP)
    prd_rules left_complement;  // For even algebras (EGA2D, PGA3DP)
    prd_rules complement;       // For odd algebras (EGA3D, PGA2DP)

    // TODO: Dual rules - will be generated from complement rules + extended metric
    // prd_rules right_dual;    // For even algebras (EGA2D, PGA3DP)
    // prd_rules left_dual;     // For even algebras (EGA2D, PGA3DP)
    // prd_rules dual;          // For odd algebras (EGA3D, PGA2DP)
};

////////////////////////////////////////////////////////////////////////////////
// Helper Functions (declarations)
////////////////////////////////////////////////////////////////////////////////

std::vector<int> parse_indices(const std::string& basis_element,
                               const std::string& prefix = "e");
std::string indices_to_basis(const std::vector<int>& indices,
                             const std::string& prefix = "e");
std::pair<std::string, int> multiply_basis_elements(const std::string& a,
                                                    const std::string& b,
                                                    const AlgebraConfig& config);

////////////////////////////////////////////////////////////////////////////////
// Extended Metric Calculation (declarations)
////////////////////////////////////////////////////////////////////////////////

std::vector<int> calculate_extended_metric(const AlgebraConfig& config);

////////////////////////////////////////////////////////////////////////////////
// Main Generation Functions (declarations)
////////////////////////////////////////////////////////////////////////////////

mvec_coeff generate_basis(const AlgebraConfig& config);
prd_rules generate_ordered_rules(
    const AlgebraConfig& config, const std::string& operator_str,
    std::function<std::pair<std::string, int>(const std::string&, const std::string&,
                                              const AlgebraConfig&)>
        multiply_func);
prd_rules generate_geometric_product_rules(const AlgebraConfig& config);
std::pair<std::string, int> multiply_basis_elements_wedge(const std::string& a,
                                                          const std::string& b,
                                                          const AlgebraConfig& config);
std::pair<std::string, int> multiply_basis_elements_dot(const std::string& a,
                                                        const std::string& b,
                                                        const AlgebraConfig& config);
prd_rules generate_wedge_product_rules(const AlgebraConfig& config);
prd_rules generate_dot_product_rules(const AlgebraConfig& config);
prd_rules generate_right_complement_rules(const AlgebraConfig& config,
                                          const prd_rules& wedge_rules);
prd_rules generate_left_complement_rules(const AlgebraConfig& config,
                                         const prd_rules& wedge_rules);
prd_rules generate_complement_rules(const AlgebraConfig& config,
                                    const prd_rules& wedge_rules);
ProductRules generate_algebra_rules(const AlgebraConfig& config);

////////////////////////////////////////////////////////////////////////////////
// Validation and Testing Functions (declarations)
////////////////////////////////////////////////////////////////////////////////

bool validate_rules(const prd_rules& generated, const prd_rules& reference);
void print_rule_comparison(const prd_rules& generated, const prd_rules& reference,
                           const std::string& product_name);