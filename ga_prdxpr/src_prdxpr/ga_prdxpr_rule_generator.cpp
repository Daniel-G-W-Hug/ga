// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_prdxpr_dual_calc.hpp"
#include "ga_prdxpr_metric_calc.hpp"
#include "ga_prdxpr_rule_generator.hpp"
#include <algorithm>
#include <functional>
#include <mdspan>
#include <numeric>
#include <set>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// Helper Functions Implementation
////////////////////////////////////////////////////////////////////////////////

// Parse basis element to extract indices (e.g., "e12" -> {1, 2})
std::vector<int> parse_indices(std::string const& basis_element,
                               std::string const& prefix)
{
    std::vector<int> indices;
    if (basis_element == one_str()) return indices; // scalar has no indices

    if (basis_element.substr(0, prefix.length()) == prefix) {
        std::string indices_str = basis_element.substr(prefix.length());
        for (char c : indices_str) {
            if (std::isdigit(c)) {
                indices.push_back(c - '0');
            }
        }
    }
    return indices;
}

// Create basis element from indices (simple concatenation)
std::string indices_to_basis(std::vector<int> const& indices, std::string const& prefix)
{
    if (indices.empty()) return one_str();

    std::string result = prefix;
    for (int idx : indices) {
        result += std::to_string(idx);
    }
    return result;
}

// Multiply two basis elements using geometric algebra rules
std::pair<std::string, int> multiply_basis_elements(std::string const& a,
                                                    std::string const& b,
                                                    AlgebraConfig const& config)
{
    // Handle scalar multiplication - use consistent string constants
    if (a == config.scalar_name) return {b, 1};
    if (b == config.scalar_name) return {a, 1};

    auto indices_a = parse_indices(a, config.basis_prefix);
    auto indices_b = parse_indices(b, config.basis_prefix);

    // Create combined sequence for geometric product calculation
    std::vector<int> sequence;
    sequence.insert(sequence.end(), indices_a.begin(), indices_a.end());
    sequence.insert(sequence.end(), indices_b.begin(), indices_b.end());

    int sign = 1;

    // Simplify the sequence by removing pairs of identical indices
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < sequence.size() && !changed; ++i) {
            for (size_t j = i + 1; j < sequence.size(); ++j) {
                if (sequence[i] == sequence[j]) {
                    // Found a pair - apply metric signature
                    int basis_index = sequence[i];
                    int metric_index = basis_index;

                    // Check if this is 1-based indexing vs 0-based indexing
                    bool is_one_based = true;
                    for (const auto& bv : config.basis_vectors) {
                        auto indices = parse_indices(bv, config.basis_prefix);
                        if (!indices.empty() && indices[0] == 0) {
                            is_one_based = false;
                            break;
                        }
                    }

                    if (is_one_based && basis_index > 0) {
                        metric_index = basis_index - 1; // Convert 1-based to 0-based
                    }

                    sign *= config.metric_signature[metric_index];

                    // Count swaps needed to bring indices together
                    int swaps = static_cast<int>(j - i - 1);
                    if (swaps % 2 == 1) {
                        sign *= -1;
                    }

                    // Remove both indices (remove larger index first to avoid shifting)
                    sequence.erase(sequence.begin() + j);
                    sequence.erase(sequence.begin() + i);
                    changed = true;
                    break;
                }
            }
        }
    }

    // If no indices remain, result is scalar
    if (sequence.empty()) {
        return {config.scalar_name, sign};
    }

    // Find the canonical form of this basis element in the user-provided basis
    for (const std::string& canonical : config.multivector_basis) {
        if (canonical == config.scalar_name) continue;

        auto canonical_indices = parse_indices(canonical, config.basis_prefix);

        // Check if same set of indices (regardless of order)
        if (sequence.size() == canonical_indices.size()) {
            std::vector<int> sorted_sequence = sequence;
            std::vector<int> sorted_canonical = canonical_indices;
            std::sort(sorted_sequence.begin(), sorted_sequence.end());
            std::sort(sorted_canonical.begin(), sorted_canonical.end());

            if (sorted_sequence == sorted_canonical) {
                // Count swaps needed to transform sequence to canonical_indices
                std::vector<int> temp_sequence = sequence;
                int swaps = 0;
                for (size_t i = 0; i < canonical_indices.size(); ++i) {
                    // Find position of canonical_indices[i] in temp_sequence
                    auto it = std::find(temp_sequence.begin() + i, temp_sequence.end(),
                                        canonical_indices[i]);
                    if (it != temp_sequence.begin() + i) {
                        // Swap to correct position
                        std::swap(temp_sequence[i], *it);
                        swaps++;
                    }
                }

                // Apply sign from swaps
                if (swaps % 2 == 1) {
                    sign *= -1;
                }

                return {canonical, sign};
            }
        }
    }

    // If no canonical form found, this shouldn't happen with correct input
    return {indices_to_basis(sequence, config.basis_prefix), sign};
}

////////////////////////////////////////////////////////////////////////////////
// Extended Metric Calculation Implementation
////////////////////////////////////////////////////////////////////////////////

// Calculate full extended metric matrix G (including off-diagonal elements)
// Returns a flattened row-major matrix stored in std::vector<int>
// The matrix is size n×n where n = basis.size()
//
// NEW IMPLEMENTATION: Uses systematic Gram matrix determinant calculation
// Mathematical definition: G[i,j] = ⟨basis[i], basis[j]⟩ = det(Gram)
//
// Properties:
// - Symmetric: G[i,j] = G[j,i]
// - Block-diagonal by grade (different grades → 0)
// - Diagonal for orthogonal bases (EGA, PGA, STA) → backward compatible
// - Non-diagonal for null bases (future CGA support)
//
// This implementation uses the new metric calculation module (ga_prdxpr_metric_calc.cpp)
// which computes ALL matrix elements (not just diagonal) via Gram determinants.
//
// For current orthogonal-basis algebras (EGA, PGA, STA):
// - Off-diagonal elements are automatically zero
// - Diagonal elements match previous implementation exactly
// - Results are 100% character-identical to old implementation
std::vector<int> calculate_extended_metric_matrix(AlgebraConfig const& config)
{
    // Use new full matrix computation from ga_prdxpr_metric_calc.cpp
    // This computes G[i,j] = ⟨basis[i], basis[j]⟩ for ALL i,j via Gram determinants
    return calculate_extended_metric_matrix_full(config);
}

std::vector<int> calculate_extended_metric(AlgebraConfig const& config)
{
    auto const& basis = config.multivector_basis;
    auto const& metric = config.metric_signature;
    std::vector<int> extended_metric(basis.size(), 0);

    // Scalar always has metric value 1
    extended_metric[0] = 1;

    // Check if this is 1-based indexing vs 0-based indexing
    bool is_one_based = true;
    for (const auto& bv : config.basis_vectors) {
        auto bv_indices = parse_indices(bv, config.basis_prefix);
        if (!bv_indices.empty() && bv_indices[0] == 0) {
            is_one_based = false;
            break;
        }
    }

    // Calculate determinant of metric for pseudoscalar
    int determinant = 1;
    for (int m : metric) {
        determinant *= m;
    }

    for (size_t i = 1; i < basis.size(); ++i) {
        const std::string& element = basis[i];
        auto indices = parse_indices(element, config.basis_prefix);

        if (indices.empty()) {
            extended_metric[i] = 0;
            continue;
        }

        // Check if this is the pseudoscalar (highest grade element)
        if (indices.size() == config.basis_vectors.size()) {
            // Pseudoscalar: use determinant of metric
            extended_metric[i] = determinant;
        }
        else {
            // For other elements: calculate using conforming property G(a ^ b) = G(a) ^
            // G(b)
            int metric_value = 1;
            int negative_count = 0; // Count negative metric values

            for (int idx : indices) {
                int metric_index = idx;
                if (is_one_based && idx > 0) {
                    metric_index = idx - 1; // Convert 1-based to 0-based
                }

                if (metric_index >= 0 && metric_index < static_cast<int>(metric.size())) {
                    int m = metric[metric_index];
                    metric_value *= m;
                    if (m < 0) negative_count++;
                }
                else {
                    metric_value = 0; // Invalid index
                    break;
                }
            }

            // Apply special rules for mixed signature algebras
            // For STA4D G(1,3,0): extended metric depends on presence of g0 (index 0)
            // For other algebras: use standard product rule
            if (config.metric_signature == std::vector<int>{+1, -1, -1, -1}) {
                // STA4D special case: check if g0 (index 0) is present
                bool has_g0 =
                    std::find(indices.begin(), indices.end(), 0) != indices.end();

                if (indices.size() == 2) {
                    // Bivectors: +1 if contains g0, -1 if only spatial indices
                    metric_value = has_g0 ? +1 : -1;
                }
                else if (indices.size() == 3) {
                    // Trivectors: -1 if contains g0, +1 if only spatial indices
                    metric_value = has_g0 ? -1 : +1;
                }
            }
            else {
                // For other algebras: apply conforming property G(a ^ b) = G(a) ^ G(b)
                // This ensures the extended metric follows the outer product structure
                if (indices.size() == 2 && negative_count == 1) {
                    // Mixed signature bivector (one positive, one negative): flip sign
                    metric_value *= -1;
                }
            }

            extended_metric[i] = metric_value;
        }
    }

    return extended_metric;
}

// Calculate extended metric using recursive extraction from product rules
// This approach avoids hardcoded special cases by extracting metric values
// directly from the algebra's own wedge and dot product operations.
//
// Algorithm:
// 1. Level 0 (scalar): always 1 by definition
// 2. Level 1 (vectors): extract from dot products g_i = dot(ei, ei)
// 3. Level k (k≥2): recursively construct from wedge products, then extract
//                    via dot product: g_basis = dot(basis_constructed, basis_constructed)
//
// This works for all algebras (EGA, STA, PGA) and naturally handles:
// - Mixed signatures (STA4D: one timelike among spacelike)
// - Null dimensions (PGA: degenerate directions)
// - Any metric signature without special-case code
//
// The conforming property G(a∧b) = (Ga)∧(Gb) ensures correctness.
std::vector<int> calculate_extended_metric_recursive(AlgebraConfig const& config,
                                                     prd_rules const& wedge_rules,
                                                     prd_rules const& dot_rules)
{
    auto const& basis = config.multivector_basis;
    std::vector<int> extended_metric(basis.size(), 0);

    // Level 0: Scalar always has metric value 1
    extended_metric[0] = 1;

    // Level 1: Vectors - extract from dot products
    // For each basis vector ei: g_i = dot(ei, ei)
    for (size_t i = 1; i < basis.size(); ++i) {
        std::string const& element = basis[i];
        auto indices = parse_indices(element, config.basis_prefix);

        if (indices.size() == 1) {
            // This is a vector - look up dot(element, element)
            std::string key = element + space_str() + mul_str() + space_str() + element;
            auto it = dot_rules.find(key);

            if (it != dot_rules.end()) {
                std::string const& result = it->second;
                // Parse the scalar coefficient (should be "+1", "-1", or "0")
                if (result == plus_str() + one_str() || result == one_str()) {
                    extended_metric[i] = +1;
                }
                else if (result == minus_str() + one_str()) {
                    extended_metric[i] = -1;
                }
                else if (result == zero_str()) {
                    extended_metric[i] = 0;
                }
                else {
                    // Unexpected result format
                    throw std::runtime_error(
                        fmt::format("Unexpected dot product result for vector {}: {}",
                                    element, result));
                }
            }
        }
    }

    // Level k (k≥2): Higher grades - recursively construct via wedge products
    // Strategy: For each basis element of grade k, find a decomposition into
    //           wedge product of two lower-grade elements, construct it,
    //           then extract metric via dot product
    //
    // For simplicity, decompose grade-k elements into: (grade-1) ∧ (grade-(k-1))
    // This works because wedge products are well-defined bottom-up.

    for (size_t i = 1; i < basis.size(); ++i) {
        std::string const& element = basis[i];
        auto indices = parse_indices(element, config.basis_prefix);

        if (indices.size() >= 2) {
            // Grade k ≥ 2: Try to find wedge decomposition
            // Strategy: Try all possible decompositions from the existing basis
            // For element with indices {i1, i2, ..., ik}, try wedging any lower-grade
            // basis elements that could produce it.
            //
            // Simpler approach: Search through all basis elements of grade (k-1)
            // and grade 1, find which wedge product produces current element.

            bool found_metric = false;
            for (auto const& basis_a : config.multivector_basis) {
                if (found_metric) break;

                auto indices_a = parse_indices(basis_a, config.basis_prefix);
                if (indices_a.empty()) continue; // Skip scalar

                for (auto const& basis_b : config.multivector_basis) {
                    if (found_metric) break;

                    auto indices_b = parse_indices(basis_b, config.basis_prefix);
                    if (indices_b.empty()) continue; // Skip scalar

                    // Check if grades add up correctly
                    if (indices_a.size() + indices_b.size() != indices.size()) continue;

                    // Look up wedge product: basis_a ∧ basis_b
                    std::string wedge_key =
                        basis_a + space_str() + wdg_str() + space_str() + basis_b;
                    auto wedge_it = wedge_rules.find(wedge_key);

                    if (wedge_it != wedge_rules.end()) {
                        std::string constructed = wedge_it->second;

                        // Remove sign prefix if present
                        std::string canonical_constructed = constructed;
                        if (constructed.length() > 0 &&
                            (constructed[0] == '+' || constructed[0] == '-')) {
                            canonical_constructed = constructed.substr(1);
                        }

                        // Check if this produces our target element
                        if (canonical_constructed == element) {
                            // Found the right decomposition!
                            // Now compute dot(element, element)
                            std::string dot_key =
                                element + space_str() + mul_str() + space_str() + element;
                            auto dot_it = dot_rules.find(dot_key);

                            if (dot_it != dot_rules.end()) {
                                std::string const& dot_result = dot_it->second;

                                // Parse the scalar result
                                if (dot_result == plus_str() + one_str() ||
                                    dot_result == one_str()) {
                                    extended_metric[i] = +1;
                                }
                                else if (dot_result == minus_str() + one_str()) {
                                    extended_metric[i] = -1;
                                }
                                else if (dot_result == zero_str()) {
                                    extended_metric[i] = 0;
                                }
                                else {
                                    throw std::runtime_error(fmt::format(
                                        "Unexpected dot product result for {}: {}",
                                        element, dot_result));
                                }

                                found_metric = true;
                                break;
                            }
                            else {
                                throw std::runtime_error(fmt::format(
                                    "Dot product rule not found for key: {}", dot_key));
                            }
                        }
                    }
                }
            }

            if (!found_metric) {
                throw std::runtime_error(fmt::format(
                    "Could not find wedge decomposition for basis element: {}", element));
            }
        }
    }

    return extended_metric;
}

// Calculate regressive extended metric (for PGA algebras)
// Regressive metric: Ḡ(u) is complementary to G(u)
// Where G(u) = 1, Ḡ(u) = 0; where G(u) = 0, Ḡ(u) = 1
//
// This is only meaningful for algebras with degenerate dimensions (PGA)
// Mathematical foundation: Ḡ(u) = cmpl(G(cmpl(u))) for PGA
// Simplified implementation: Ḡ[i] = 1 - G[i] (works for diagonal orthogonal metrics)
//
// LIMITATION: Only works for diagonal metric matrices
//   - PGA2DP, PGA3DP: diagonal with 0s and 1s ✓
//   - Future CGA: NON-diagonal, this simple inversion won't work!
//
// For CGA support, would need to implement:
//   - Full matrix complement transformation: Ḡ = cmpl · G · cmpl
//   - Where cmpl represents the complement operation as a matrix transformation
std::vector<int> calculate_regressive_extended_metric(AlgebraConfig const& config)
{
    auto const& basis = config.multivector_basis;

    // Get standard extended metric
    auto standard_metric = calculate_extended_metric(config);

    // Regressive metric is complementary: Ḡ = 1 - G (for values 0 and 1)
    std::vector<int> regressive_metric(basis.size(), 0);
    for (size_t i = 0; i < basis.size(); ++i) {
        regressive_metric[i] = (standard_metric[i] == 0) ? 1 : 0;
    }

    return regressive_metric;
}

// Calculate regressive extended metric matrix (for PGA algebras)
// Returns the full matrix where diagonal elements are inverted (0<->1)
// Off-diagonal elements remain 0 (both G and Ḡ are diagonal matrices)
//
// LIMITATION: Only inverts diagonal elements, assumes orthogonal basis
//   - Works for PGA2DP, PGA3DP (diagonal matrices with 0s and 1s) ✓
//   - Will NOT work for CGA with non-diagonal metric matrix!
//
// Mathematical relationship for diagonal PGA metrics:
//   - Extended metric G: diagonal with G[i,i] ∈ {0, 1}
//   - Regressive metric Ḡ: diagonal with Ḡ[i,i] = 1 - G[i,i]
//   - Both matrices are diagonal (all off-diagonal elements = 0)
//   - Relationship: G[i,i] + Ḡ[i,i] = 1 (complementary on diagonal)
//
// For future CGA support with non-diagonal metrics:
//   - Cannot use simple element-wise inversion
//   - Need matrix-based complement transformation
//   - Or compute via: Ḡ = cmpl_matrix · G · cmpl_matrix^T
//   - Where cmpl_matrix represents complement as linear transformation
std::vector<int> calculate_regressive_extended_metric_matrix(AlgebraConfig const& config)
{
    auto const& basis = config.multivector_basis;
    size_t const n = basis.size();

    // Get standard extended metric matrix
    auto standard_matrix = calculate_extended_metric_matrix(config);

    // Create mdspan view for standard matrix
    std::mdspan G{standard_matrix.data(), n, n};

    // Regressive metric matrix: invert only diagonal elements, off-diagonal stays 0
    std::vector<int> regressive_matrix(n * n, 0);
    std::mdspan G_reg{regressive_matrix.data(), n, n};

    // Invert only diagonal elements: Ḡ[i,i] = 1 - G[i,i]
    for (size_t i = 0; i < n; ++i) {
        G_reg[i, i] = (G[i, i] == 0) ? 1 : 0;
    }
    // Off-diagonal elements remain 0 (already initialized)
    // NOTE: For CGA, off-diagonal elements of both G and Ḡ would be non-zero!

    return regressive_matrix;
}

////////////////////////////////////////////////////////////////////////////////
// Main Generation Functions Implementation
////////////////////////////////////////////////////////////////////////////////

mvec_coeff generate_basis(AlgebraConfig const& config)
{
    return config.multivector_basis;
}

prd_rules generate_ordered_rules(
    AlgebraConfig const& config, std::string const& operator_str,
    std::function<std::pair<std::string, int>(std::string const&, std::string const&,
                                              AlgebraConfig const&)>
        multiply_func)
{
    // Use ordered map to maintain insertion order which follows our grade-ordered basis
    std::map<std::string, std::string> rules;
    auto basis = generate_basis(config);

    // Generate rules in grade order by iterating through ordered basis
    for (const auto& a : basis) {
        for (const auto& b : basis) {
            auto [result, sign] = multiply_func(a, b, config);

            // Create key using consistent operators from common header
            std::string key = a + space_str() + operator_str + space_str() + b;
            std::string value;

            if (sign == 1) {
                value = result;
            }
            else if (sign == -1) {
                value = minus_str() + result;
            }
            else {
                value = zero_str(); // For null results
            }

            rules[key] = value;
        }
    }

    return rules;
}

prd_rules generate_geometric_product_rules(AlgebraConfig const& config)
{
    return generate_ordered_rules(config, mul_str(), multiply_basis_elements);
}

std::pair<std::string, int> multiply_basis_elements_wedge(std::string const& a,
                                                          std::string const& b,
                                                          AlgebraConfig const& config)
{
    auto indices_a = parse_indices(a, config.basis_prefix);
    auto indices_b = parse_indices(b, config.basis_prefix);

    // Check for repeated indices (wedge product = 0 for repeated elements)
    std::set<int> set_a(indices_a.begin(), indices_a.end());
    std::set<int> set_b(indices_b.begin(), indices_b.end());

    std::vector<int> intersection;
    std::set_intersection(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                          std::back_inserter(intersection));

    if (!intersection.empty()) {
        return {zero_str(), 0};
    }
    else {
        // Use geometric product for non-intersecting indices
        return multiply_basis_elements(a, b, config);
    }
}

std::pair<std::string, int> multiply_basis_elements_dot(std::string const& a,
                                                        std::string const& b,
                                                        AlgebraConfig const& config)
{
    // Handle scalar cases directly
    if (a == config.scalar_name && b == config.scalar_name) {
        return {config.scalar_name, 1};
    }
    else if (a == config.scalar_name || b == config.scalar_name) {
        return {zero_str(), 0};
    }

    // Dot product is only non-zero when both elements are identical
    if (a != b) {
        return {zero_str(), 0};
    }

    // For identical elements: use the extended metric value
    auto extended_metric = calculate_extended_metric(config);

    // Find the index of element a in the basis
    auto it =
        std::find(config.multivector_basis.begin(), config.multivector_basis.end(), a);
    if (it == config.multivector_basis.end()) {
        return {zero_str(), 0}; // Element not found in basis
    }

    size_t index = std::distance(config.multivector_basis.begin(), it);
    if (index >= extended_metric.size()) {
        return {zero_str(), 0}; // Index out of range
    }

    int metric_value = extended_metric[index];
    if (metric_value == 0) {
        return {zero_str(), 0};
    }
    else {
        return {config.scalar_name, metric_value};
    }
}

prd_rules generate_wedge_product_rules(AlgebraConfig const& config)
{
    return generate_ordered_rules(config, wdg_str(), multiply_basis_elements_wedge);
}

prd_rules generate_dot_product_rules(AlgebraConfig const& config)
{
    return generate_ordered_rules(config, mul_str(), multiply_basis_elements_dot);
}

////////////////////////////////////////////////////////////////////////////////
// Complement Rule Generation Implementation
////////////////////////////////////////////////////////////////////////////////

// Generate complement rules from wedge product table
// Algorithm: For complement relationship u ^ cmpl(u) = I_n (right complement)
//            or cmpl(u) ^ u = I_n (left complement)
prd_rules generate_complement_from_wedge_table(AlgebraConfig const& config,
                                               prd_rules const& wedge_rules,
                                               bool is_left_complement)
{
    prd_rules complement_rules;

    // Get pseudoscalar (last element of multivector basis)
    if (config.multivector_basis.empty()) {
        throw std::runtime_error("Empty multivector basis in config");
    }

    std::string pseudoscalar = config.multivector_basis.back();

    // Create wedge product table from rules - this is the key insight!
    // Step 1: Create coefficient table structure
    prd_table wedge_table =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.multivector_basis,
                                                     config.multivector_basis, wdg_str()),
                           wedge_rules);

    size_t basis_size = config.multivector_basis.size();

    // For each basis element, find its complement using table lookup
    for (size_t i = 0; i < basis_size; ++i) {
        const std::string& basis_element = config.multivector_basis[i];
        bool found_complement = false;

        // Special case: scalar complement is always the pseudoscalar
        if (basis_element == config.scalar_name) {
            complement_rules[basis_element] = pseudoscalar;
            found_complement = true;
        }
        // Special case: pseudoscalar complement is always the scalar
        else if (basis_element == pseudoscalar) {
            complement_rules[basis_element] = config.scalar_name;
            found_complement = true;
        }
        else {
            if (is_left_complement) {
                // Left complement: lcmpl(u) ^ u = I_n
                // Search column i for pseudoscalar
                for (size_t row = 0; row < basis_size; ++row) {
                    if (row < wedge_table.size() && i < wedge_table[row].size()) {
                        const std::string& table_entry = wedge_table[row][i];
                        if (table_entry == pseudoscalar) {
                            complement_rules[basis_element] =
                                config.multivector_basis[row];
                            found_complement = true;
                            break;
                        }
                        else if (table_entry == "-" + pseudoscalar) {
                            complement_rules[basis_element] =
                                "-" + config.multivector_basis[row];
                            found_complement = true;
                            break;
                        }
                    }
                }
            }
            else {
                // Right complement: u ^ rcmpl(u) = I_n
                // Search row i for pseudoscalar
                if (i < wedge_table.size()) {
                    for (size_t col = 0; col < basis_size && col < wedge_table[i].size();
                         ++col) {
                        const std::string& table_entry = wedge_table[i][col];
                        if (table_entry == pseudoscalar) {
                            complement_rules[basis_element] =
                                config.multivector_basis[col];
                            found_complement = true;
                            break;
                        }
                        else if (table_entry == "-" + pseudoscalar) {
                            complement_rules[basis_element] =
                                "-" + config.multivector_basis[col];
                            found_complement = true;
                            break;
                        }
                    }
                }
            }
        }

        if (!found_complement) {
            throw std::runtime_error("Could not find complement for basis element: " +
                                     basis_element);
        }
    }

    return complement_rules;
}

prd_rules generate_right_complement_rules(AlgebraConfig const& config,
                                          prd_rules const& wedge_rules)
{
    return generate_complement_from_wedge_table(config, wedge_rules, false);
}

prd_rules generate_left_complement_rules(AlgebraConfig const& config,
                                         prd_rules const& wedge_rules)
{
    return generate_complement_from_wedge_table(config, wedge_rules, true);
}

prd_rules generate_complement_rules(AlgebraConfig const& config,
                                    prd_rules const& wedge_rules)
{
    // For odd algebras, left and right complements are the same
    return generate_complement_from_wedge_table(config, wedge_rules, false);
}

////////////////////////////////////////////////////////////////////////////////
// Dual Rule Generation Implementation
////////////////////////////////////////////////////////////////////////////////

// Generate left dual rules: left_dual(u) = left_complement(G × u)
// where G is the extended metric matrix and × is matrix-vector multiplication
prd_rules generate_left_dual_rules(AlgebraConfig const& config,
                                   prd_rules const& left_complement_rules)
{
    // Delegate to systematic matrix-vector multiplication approach
    // left_dual(basis[i]) = Σⱼ G[i,j] · left_complement(basis[j])
    auto G_data = calculate_extended_metric_matrix_full(config);
    return calculate_dual_rules(config, G_data, left_complement_rules);
}

// Generate right dual rules: right_dual(u) = right_complement(G × u)
prd_rules generate_right_dual_rules(AlgebraConfig const& config,
                                    prd_rules const& right_complement_rules)
{
    // Delegate to systematic matrix-vector multiplication approach
    // right_dual(basis[i]) = Σⱼ G[i,j] · right_complement(basis[j])
    auto G_data = calculate_extended_metric_matrix_full(config);
    return calculate_dual_rules(config, G_data, right_complement_rules);
}

// Generate dual rules for odd-dimensional algebras: dual(u) = complement(G × u)
prd_rules generate_dual_rules(AlgebraConfig const& config,
                              prd_rules const& complement_rules)
{
    // Step 1: Calculate extended metric matrix (full matrix, not just diagonal)
    auto G_data = calculate_extended_metric_matrix_full(config);

    // Step 2: Calculate dual rules via matrix-vector multiplication
    // dual(basis[i]) = Σⱼ G[i,j] * complement(basis[j])
    return calculate_dual_rules(config, G_data, complement_rules);
}

////////////////////////////////////////////////////////////////////////////////
// PGA-Specific Dual Generation (Bulk and Weight Duals)
////////////////////////////////////////////////////////////////////////////////

// Generate bulk_dual rules for odd-dimensional PGA (uses complement)
// bulk_dual(u) = complement(G × u) where G is the extended metric matrix
prd_rules generate_bulk_dual_rules(AlgebraConfig const& config,
                                   prd_rules const& complement_rules)
{
    // For odd-dimensional PGA, bulk_dual is the same as regular dual
    return generate_dual_rules(config, complement_rules);
}

// Generate weight_dual rules for odd-dimensional PGA
// weight_dual(u) = complement(Ḡ × u) where Ḡ = regressive metric
// Regressive metric: Ḡ × u = complement(G × complement(u))
// Therefore: weight_dual(u) = complement(complement(G × complement(u))) = G ×
// complement(u)
prd_rules generate_weight_dual_rules(AlgebraConfig const& config,
                                     prd_rules const& complement_rules)
{
    // For odd-dimensional PGA: weight_dual(u) = G × complement(u)
    mvec_coeff const& basis = config.multivector_basis;
    size_t const n = basis.size();
    auto G_data = calculate_extended_metric_matrix(config);
    std::mdspan G{G_data.data(), n, n};

    prd_rules weight_dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Step 1: Apply complement to input basis element
        auto it_cmpl = complement_rules.find(basis_element);
        if (it_cmpl == complement_rules.end()) {
            weight_dual_rules[basis_element] = zero_str();
            continue;
        }
        std::string cmpl_u = it_cmpl->second;

        if (cmpl_u == zero_str()) {
            weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 2: Find the index of complement(u) in the basis (strip sign if present)
        bool cmpl_has_minus = (cmpl_u.find(minus_str()) == 0);
        std::string cmpl_u_unsigned = cmpl_has_minus ? cmpl_u.substr(1) : cmpl_u;

        size_t cmpl_index = n; // invalid index
        for (size_t j = 0; j < n; ++j) {
            if (basis[j] == cmpl_u_unsigned) {
                cmpl_index = j;
                break;
            }
        }

        if (cmpl_index >= n) {
            // Complement not found in basis (shouldn't happen)
            weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 3: Multiply by extended metric G[complement_index, complement_index]
        // Result is G × complement(u)
        int metric_value = G[cmpl_index, cmpl_index];

        if (metric_value == 1) {
            // metric × complement = complement (preserve sign)
            weight_dual_rules[basis_element] = cmpl_u;
        }
        else if (metric_value == -1) {
            // metric × complement = -complement (flip sign)
            if (cmpl_has_minus) {
                weight_dual_rules[basis_element] = cmpl_u_unsigned; // remove minus
            }
            else {
                weight_dual_rules[basis_element] = minus_str() + cmpl_u; // add minus
            }
        }
        else if (metric_value == 0) {
            weight_dual_rules[basis_element] = zero_str();
        }
        else {
            // General case: metric_value × complement
            weight_dual_rules[basis_element] =
                std::to_string(metric_value) + mul_str() + cmpl_u;
        }
    }

    return weight_dual_rules;
}

// Generate left_bulk_dual rules for even-dimensional PGA
prd_rules generate_left_bulk_dual_rules(AlgebraConfig const& config,
                                        prd_rules const& left_complement_rules)
{
    return generate_left_dual_rules(config, left_complement_rules);
}

// Generate right_bulk_dual rules for even-dimensional PGA
prd_rules generate_right_bulk_dual_rules(AlgebraConfig const& config,
                                         prd_rules const& right_complement_rules)
{
    return generate_right_dual_rules(config, right_complement_rules);
}

// Generate left_weight_dual rules for even-dimensional PGA
// left_weight_dual(u) = left_complement(Ḡ × u) where Ḡ = regressive metric
// Regressive metric: Ḡ × u = left_complement(G × left_complement(u))
// Therefore: left_weight_dual(u) = left_complement(left_complement(G ×
// left_complement(u)))
//                                  = G × left_complement(u) (using left_complement
//                                  involution)
prd_rules generate_left_weight_dual_rules(AlgebraConfig const& config,
                                          prd_rules const& left_complement_rules)
{
    // For even-dimensional PGA: left_weight_dual(u) = G × left_complement(u)
    mvec_coeff const& basis = config.multivector_basis;
    size_t const n = basis.size();
    auto G_data = calculate_extended_metric_matrix(config);
    std::mdspan G{G_data.data(), n, n};

    prd_rules left_weight_dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Step 1: Apply left_complement to input basis element
        auto it_lcmpl = left_complement_rules.find(basis_element);
        if (it_lcmpl == left_complement_rules.end()) {
            left_weight_dual_rules[basis_element] = zero_str();
            continue;
        }
        std::string lcmpl_u = it_lcmpl->second;

        if (lcmpl_u == zero_str()) {
            left_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 2: Find the index of left_complement(u) in the basis (strip sign if
        // present)
        bool lcmpl_has_minus = (lcmpl_u.find(minus_str()) == 0);
        std::string lcmpl_u_unsigned = lcmpl_has_minus ? lcmpl_u.substr(1) : lcmpl_u;

        size_t lcmpl_index = n; // invalid index
        for (size_t j = 0; j < n; ++j) {
            if (basis[j] == lcmpl_u_unsigned) {
                lcmpl_index = j;
                break;
            }
        }

        if (lcmpl_index >= n) {
            left_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 3: Multiply by extended metric G[left_complement_index,
        // left_complement_index] Result is G × left_complement(u)
        int metric_value = G[lcmpl_index, lcmpl_index];

        if (metric_value == 1) {
            // metric × left_complement = left_complement (preserve sign)
            left_weight_dual_rules[basis_element] = lcmpl_u;
        }
        else if (metric_value == -1) {
            // metric × left_complement = -left_complement (flip sign)
            if (lcmpl_has_minus) {
                left_weight_dual_rules[basis_element] = lcmpl_u_unsigned; // remove minus
            }
            else {
                left_weight_dual_rules[basis_element] =
                    minus_str() + lcmpl_u; // add minus
            }
        }
        else if (metric_value == 0) {
            left_weight_dual_rules[basis_element] = zero_str();
        }
        else {
            // General case: metric_value × left_complement
            left_weight_dual_rules[basis_element] =
                std::to_string(metric_value) + mul_str() + lcmpl_u;
        }
    }

    return left_weight_dual_rules;
}

// Generate right_weight_dual rules for even-dimensional PGA
// right_weight_dual(u) = right_complement(Ḡ × u) where Ḡ = regressive metric
// Regressive metric: Ḡ × u = left_complement(G × right_complement(u))
// Therefore: right_weight_dual(u) = right_complement(left_complement(G ×
// right_complement(u)))
//                                  = G × right_complement(u) (using complement
//                                  involution)
prd_rules generate_right_weight_dual_rules(AlgebraConfig const& config,
                                           prd_rules const& right_complement_rules)
{
    // For even-dimensional PGA: right_weight_dual(u) = G × right_complement(u)
    mvec_coeff const& basis = config.multivector_basis;
    size_t const n = basis.size();
    auto G_data = calculate_extended_metric_matrix(config);
    std::mdspan G{G_data.data(), n, n};

    prd_rules right_weight_dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Step 1: Apply right_complement to input basis element
        auto it_rcmpl = right_complement_rules.find(basis_element);
        if (it_rcmpl == right_complement_rules.end()) {
            right_weight_dual_rules[basis_element] = zero_str();
            continue;
        }
        std::string rcmpl_u = it_rcmpl->second;

        if (rcmpl_u == zero_str()) {
            right_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 2: Find the index of right_complement(u) in the basis (strip sign if
        // present)
        bool rcmpl_has_minus = (rcmpl_u.find(minus_str()) == 0);
        std::string rcmpl_u_unsigned = rcmpl_has_minus ? rcmpl_u.substr(1) : rcmpl_u;

        size_t rcmpl_index = n; // invalid index
        for (size_t j = 0; j < n; ++j) {
            if (basis[j] == rcmpl_u_unsigned) {
                rcmpl_index = j;
                break;
            }
        }

        if (rcmpl_index >= n) {
            right_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 3: Multiply by extended metric G[right_complement_index,
        // right_complement_index] Result is G × right_complement(u)
        int metric_value = G[rcmpl_index, rcmpl_index];

        if (metric_value == 1) {
            // metric × right_complement = right_complement (preserve sign)
            right_weight_dual_rules[basis_element] = rcmpl_u;
        }
        else if (metric_value == -1) {
            // metric × right_complement = -right_complement (flip sign)
            if (rcmpl_has_minus) {
                right_weight_dual_rules[basis_element] = rcmpl_u_unsigned; // remove minus
            }
            else {
                right_weight_dual_rules[basis_element] =
                    minus_str() + rcmpl_u; // add minus
            }
        }
        else if (metric_value == 0) {
            right_weight_dual_rules[basis_element] = zero_str();
        }
        else {
            // General case: metric_value × right_complement
            right_weight_dual_rules[basis_element] =
                std::to_string(metric_value) + mul_str() + rcmpl_u;
        }
    }

    return right_weight_dual_rules;
}

ProductRules generate_algebra_rules(AlgebraConfig const& config)
{
    // Validate input
    if (config.basis_vectors.size() != config.metric_signature.size()) {
        throw std::runtime_error(
            "Number of basis vectors must match metric signature size");
    }

    ProductRules result;
    result.basis = generate_basis(config);
    result.geometric_product = generate_geometric_product_rules(config);
    result.wedge_product = generate_wedge_product_rules(config);
    result.dot_product = generate_dot_product_rules(config);

    // Generate complement rules based on algebra dimensionality
    // Even-dimensional algebras (EGA2D: 2D, STA4D: 4D, PGA3DP: 4D) have left and right
    // complements Odd-dimensional algebras (EGA3D: 3D, PGA2DP: 3D) have a single
    // complement
    size_t num_basis_vectors = config.basis_vectors.size();
    bool is_even_dimensional = (num_basis_vectors % 2 == 0);

    // Detect PGA algebras (have degenerate dimensions with metric signature 0)
    bool is_pga = false;
    for (int metric_val : config.metric_signature) {
        if (metric_val == 0) {
            is_pga = true;
            break;
        }
    }

    if (is_even_dimensional) {
        // Generate both left and right complements for even algebras
        result.left_complement =
            generate_left_complement_rules(config, result.wedge_product);
        result.right_complement =
            generate_right_complement_rules(config, result.wedge_product);

        // Generate left and right dual rules from complements (only for non-PGA)
        if (!is_pga) {
            result.left_dual = generate_left_dual_rules(config, result.left_complement);
            result.right_dual =
                generate_right_dual_rules(config, result.right_complement);
        }

        // For even-dimensional PGA, generate bulk and weight duals
        if (is_pga) {
            result.left_bulk_dual =
                generate_left_bulk_dual_rules(config, result.left_complement);
            result.right_bulk_dual =
                generate_right_bulk_dual_rules(config, result.right_complement);
            result.left_weight_dual =
                generate_left_weight_dual_rules(config, result.left_complement);
            result.right_weight_dual =
                generate_right_weight_dual_rules(config, result.right_complement);
        }
    }
    else {
        // Generate single complement for odd algebras
        result.complement = generate_complement_rules(config, result.wedge_product);

        // Generate dual rules from complement (only for non-PGA)
        if (!is_pga) {
            result.dual = generate_dual_rules(config, result.complement);
        }

        // For odd-dimensional PGA, generate bulk and weight duals
        if (is_pga) {
            result.bulk_dual = generate_bulk_dual_rules(config, result.complement);
            result.weight_dual = generate_weight_dual_rules(config, result.complement);
        }
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Validation and Testing Functions Implementation
////////////////////////////////////////////////////////////////////////////////

bool validate_rules(prd_rules const& generated, prd_rules const& reference)
{
    if (generated.size() != reference.size()) {
        return false;
    }

    for (const auto& [key, value] : reference) {
        auto it = generated.find(key);
        if (it == generated.end() || it->second != value) {
            return false;
        }
    }

    return true;
}

void print_rule_comparison(prd_rules const& generated, prd_rules const& reference,
                           std::string const& product_name)
{
    fmt::println("Validating {} rules:", product_name);

    bool all_match = true;
    for (const auto& [key, ref_value] : reference) {
        auto it = generated.find(key);
        if (it == generated.end()) {
            fmt::println("  MISSING: {} -> {}", key, ref_value);
            all_match = false;
        }
        else if (it->second != ref_value) {
            fmt::println("  MISMATCH: {} -> generated: '{}', reference: '{}'", key,
                         it->second, ref_value);
            all_match = false;
        }
    }

    if (all_match) {
        fmt::println("  ✓ All {} rules match reference", product_name);
    }
    else {
        fmt::println("  ✗ Some {} rules differ from reference", product_name);
    }
}