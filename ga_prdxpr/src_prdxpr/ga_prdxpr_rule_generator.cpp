// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_rule_generator.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <set>
#include <mdspan>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// Helper Functions Implementation
////////////////////////////////////////////////////////////////////////////////

// Parse basis element to extract indices (e.g., "e12" -> {1, 2})
std::vector<int> parse_indices(const std::string& basis_element,
                               const std::string& prefix)
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
std::string indices_to_basis(const std::vector<int>& indices, const std::string& prefix)
{
    if (indices.empty()) return one_str();

    std::string result = prefix;
    for (int idx : indices) {
        result += std::to_string(idx);
    }
    return result;
}

// Multiply two basis elements using geometric algebra rules
std::pair<std::string, int> multiply_basis_elements(const std::string& a,
                                                    const std::string& b,
                                                    const AlgebraConfig& config)
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

// Calculate full extended metric matrix G (not just diagonal)
// Returns a flattened row-major matrix stored in std::vector<int>
// The matrix is size n×n where n = basis.size()
//
// Mathematical definition:
// - G[0,0] = 1 (scalar), all other G[0,i] = G[i,0] = 0
// - For grade 1 vectors: G[i,i] = metric_signature[i-1] (classical metric)
// - For higher grades: G(e_i ∧ e_j) = G(e_i) ∧ G(e_j) (conforming property)
//
// The matrix is block diagonal by grade in canonical basis ordering
std::vector<int> calculate_extended_metric_matrix(const AlgebraConfig& config)
{
    const auto& basis = config.multivector_basis;
    const size_t n = basis.size();

    // Flattened row-major storage: matrix[i,j] = data[i*n + j]
    std::vector<int> matrix_data(n * n, 0);

    // Create mdspan view with square bracket access: G[i,j]
    std::mdspan G{matrix_data.data(), n, n};

    // Check if this is 1-based indexing vs 0-based indexing
    bool is_one_based = true;
    for (const auto& bv : config.basis_vectors) {
        auto bv_indices = parse_indices(bv, config.basis_prefix);
        if (!bv_indices.empty() && bv_indices[0] == 0) {
            is_one_based = false;
            break;
        }
    }

    // Grade 0 (scalar): G[0,0] = 1
    G[0, 0] = 1;

    // Calculate determinant of metric for pseudoscalar
    int determinant = 1;
    for (int m : config.metric_signature) {
        determinant *= m;
    }

    // For each basis element, calculate its extended metric value
    for (size_t i = 1; i < n; ++i) {
        const std::string& element = basis[i];
        auto indices = parse_indices(element, config.basis_prefix);

        if (indices.empty()) {
            // Shouldn't happen for valid basis elements
            G[i, i] = 0;
            continue;
        }

        // Check if this is the pseudoscalar (highest grade element)
        if (indices.size() == config.basis_vectors.size()) {
            // Pseudoscalar: use determinant of metric
            G[i, i] = determinant;
        }
        else if (indices.size() == 1) {
            // Grade 1 (vectors): use classical metric directly
            int idx = indices[0];
            int metric_index = is_one_based ? (idx - 1) : idx;

            if (metric_index >= 0 &&
                metric_index < static_cast<int>(config.metric_signature.size())) {
                G[i, i] = config.metric_signature[metric_index];
            }
            else {
                G[i, i] = 0; // Degenerate direction
            }
        }
        else {
            // Higher grades: calculate using conforming property G(a ∧ b) = G(a) ∧ G(b)
            // For the diagonal, this is the product of constituent vector metrics
            int metric_value = 1;

            for (int idx : indices) {
                int metric_index = is_one_based ? (idx - 1) : idx;

                if (metric_index >= 0 &&
                    metric_index < static_cast<int>(config.metric_signature.size())) {
                    metric_value *= config.metric_signature[metric_index];
                }
                else {
                    metric_value = 0; // Degenerate direction
                    break;
                }
            }

            // Apply special rules for mixed signature algebras
            // For STA3D G(1,3,0): extended metric depends on presence of g0 (index 0)
            if (config.metric_signature == std::vector<int>{+1, -1, -1, -1}) {
                // STA3D special case: check if g0 (index 0) is present
                bool has_g0 = std::find(indices.begin(), indices.end(), 0) != indices.end();

                if (indices.size() == 2) {
                    // Bivectors: +1 if contains g0, -1 if only spatial indices
                    metric_value = has_g0 ? +1 : -1;
                }
                else if (indices.size() == 3) {
                    // Trivectors: -1 if contains g0, +1 if only spatial indices
                    metric_value = has_g0 ? -1 : +1;
                }
            }

            G[i, i] = metric_value;
        }
    }

    // LIMITATION: Off-diagonal elements remain 0 (orthogonal basis assumption)
    // This is correct for:
    //   - EGA (Euclidean GA): orthonormal basis, G[i,j] = δ_ij for i≠j
    //   - PGA (Projective GA): orthogonal basis with degenerate directions
    //   - STA (Space-Time Algebra): orthogonal mixed signature basis
    //
    // NOT supported by current implementation:
    //   - CGA (Conformal GA): requires non-diagonal metric matrix
    //     Example CGA3D has null vectors e₊, e₋ with:
    //       e₊·e₊ = 0, e₋·e₋ = 0 (diagonal)
    //       e₊·e₋ = -1 (OFF-DIAGONAL, not supported!)
    //
    // To support CGA in the future:
    //   1. Extend AlgebraConfig to include full metric matrix (not just diagonal signature)
    //   2. Add off-diagonal calculation loop after diagonal calculation:
    //      for (size_t i = 0; i < n; ++i) {
    //        for (size_t j = i+1; j < n; ++j) {
    //          // Compute G[i,j] = basis[i] · basis[j] using geometric product
    //          // Then set G[j,i] = G[i,j] (symmetric matrix)
    //        }
    //      }
    //   3. Use multiply_basis_elements() to compute dot products between basis elements
    //   4. Update regressive metric calculation (can't use simple 1-G inversion)

    return matrix_data;
}

std::vector<int> calculate_extended_metric(const AlgebraConfig& config)
{
    const auto& basis = config.multivector_basis;
    const auto& metric = config.metric_signature;
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
            // For STA3D G(1,3,0): extended metric depends on presence of g0 (index 0)
            // For other algebras: use standard product rule
            if (config.metric_signature == std::vector<int>{+1, -1, -1, -1}) {
                // STA3D special case: check if g0 (index 0) is present
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
std::vector<int> calculate_regressive_extended_metric(const AlgebraConfig& config)
{
    const auto& basis = config.multivector_basis;

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
std::vector<int> calculate_regressive_extended_metric_matrix(const AlgebraConfig& config)
{
    const auto& basis = config.multivector_basis;
    const size_t n = basis.size();

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

mvec_coeff generate_basis(const AlgebraConfig& config)
{
    return config.multivector_basis;
}

prd_rules generate_ordered_rules(
    const AlgebraConfig& config, const std::string& operator_str,
    std::function<std::pair<std::string, int>(const std::string&, const std::string&,
                                              const AlgebraConfig&)>
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

prd_rules generate_geometric_product_rules(const AlgebraConfig& config)
{
    return generate_ordered_rules(config, mul_str(), multiply_basis_elements);
}

std::pair<std::string, int> multiply_basis_elements_wedge(const std::string& a,
                                                          const std::string& b,
                                                          const AlgebraConfig& config)
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

std::pair<std::string, int> multiply_basis_elements_dot(const std::string& a,
                                                        const std::string& b,
                                                        const AlgebraConfig& config)
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

prd_rules generate_wedge_product_rules(const AlgebraConfig& config)
{
    return generate_ordered_rules(config, wdg_str(), multiply_basis_elements_wedge);
}

prd_rules generate_dot_product_rules(const AlgebraConfig& config)
{
    return generate_ordered_rules(config, mul_str(), multiply_basis_elements_dot);
}

////////////////////////////////////////////////////////////////////////////////
// Complement Rule Generation Implementation
////////////////////////////////////////////////////////////////////////////////

// Generate complement rules from wedge product table
// Algorithm: For complement relationship u ^ cmpl(u) = I_n (right complement)
//            or cmpl(u) ^ u = I_n (left complement)
prd_rules generate_complement_from_wedge_table(const AlgebraConfig& config,
                                               const prd_rules& wedge_rules,
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

prd_rules generate_right_complement_rules(const AlgebraConfig& config,
                                          const prd_rules& wedge_rules)
{
    return generate_complement_from_wedge_table(config, wedge_rules, false);
}

prd_rules generate_left_complement_rules(const AlgebraConfig& config,
                                         const prd_rules& wedge_rules)
{
    return generate_complement_from_wedge_table(config, wedge_rules, true);
}

prd_rules generate_complement_rules(const AlgebraConfig& config,
                                    const prd_rules& wedge_rules)
{
    // For odd algebras, left and right complements are the same
    return generate_complement_from_wedge_table(config, wedge_rules, false);
}

////////////////////////////////////////////////////////////////////////////////
// Dual Rule Generation Implementation
////////////////////////////////////////////////////////////////////////////////

// Generate left dual rules: left_dual(u) = left_complement(G × u)
// where G is the extended metric matrix and × is matrix-vector multiplication
prd_rules generate_left_dual_rules(const AlgebraConfig& config,
                                    const prd_rules& left_complement_rules)
{
    const auto& basis = config.multivector_basis;
    const size_t n = basis.size();

    // Calculate extended metric matrix
    auto matrix_data = calculate_extended_metric_matrix(config);
    std::mdspan G{matrix_data.data(), n, n};

    prd_rules dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Perform matrix-vector multiplication: result = G × basis_element
        // The basis_element is treated as a vector with 1 at position i, 0 elsewhere
        // So G × basis_element gives column i of G

        // Find the non-zero result (should be G[j,i] where j is the row index)
        // Since G is block diagonal for orthogonal bases, only G[i,i] is non-zero
        int metric_value = G[i, i];

        std::string intermediate_result;
        if (metric_value == 1) {
            intermediate_result = basis_element;
        }
        else if (metric_value == -1) {
            intermediate_result = minus_str() + basis_element;
        }
        else if (metric_value == 0) {
            intermediate_result = zero_str();
        }
        else {
            // Shouldn't happen for standard GAs, but handle anyway
            intermediate_result = std::to_string(metric_value) + mul_str() + basis_element;
        }

        // Apply left complement to the intermediate result
        if (intermediate_result == zero_str()) {
            dual_rules[basis_element] = zero_str();
        }
        else {
            // Extract sign and basis element from intermediate result
            bool has_minus = (intermediate_result.find(minus_str()) == 0);
            std::string element_to_complement =
                has_minus ? intermediate_result.substr(1) : intermediate_result;

            // Look up complement
            auto it = left_complement_rules.find(element_to_complement);
            if (it != left_complement_rules.end()) {
                std::string complement_value = it->second;

                // Apply sign from intermediate result
                if (has_minus) {
                    if (complement_value.find(minus_str()) == 0) {
                        // Double negative becomes positive
                        complement_value = complement_value.substr(1);
                    }
                    else if (complement_value != zero_str()) {
                        // Add negative
                        complement_value = minus_str() + complement_value;
                    }
                }

                dual_rules[basis_element] = complement_value;
            }
            else {
                // Shouldn't happen
                dual_rules[basis_element] = zero_str();
            }
        }
    }

    return dual_rules;
}

// Generate right dual rules: right_dual(u) = right_complement(G × u)
prd_rules generate_right_dual_rules(const AlgebraConfig& config,
                                     const prd_rules& right_complement_rules)
{
    const auto& basis = config.multivector_basis;
    const size_t n = basis.size();

    // Calculate extended metric matrix
    auto matrix_data = calculate_extended_metric_matrix(config);
    std::mdspan G{matrix_data.data(), n, n};

    prd_rules dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Perform matrix-vector multiplication: result = G × basis_element
        int metric_value = G[i, i];

        std::string intermediate_result;
        if (metric_value == 1) {
            intermediate_result = basis_element;
        }
        else if (metric_value == -1) {
            intermediate_result = minus_str() + basis_element;
        }
        else if (metric_value == 0) {
            intermediate_result = zero_str();
        }
        else {
            intermediate_result = std::to_string(metric_value) + mul_str() + basis_element;
        }

        // Apply right complement to the intermediate result
        if (intermediate_result == zero_str()) {
            dual_rules[basis_element] = zero_str();
        }
        else {
            // Extract sign and basis element from intermediate result
            bool has_minus = (intermediate_result.find(minus_str()) == 0);
            std::string element_to_complement =
                has_minus ? intermediate_result.substr(1) : intermediate_result;

            // Look up complement
            auto it = right_complement_rules.find(element_to_complement);
            if (it != right_complement_rules.end()) {
                std::string complement_value = it->second;

                // Apply sign from intermediate result
                if (has_minus) {
                    if (complement_value.find(minus_str()) == 0) {
                        complement_value = complement_value.substr(1);
                    }
                    else if (complement_value != zero_str()) {
                        complement_value = minus_str() + complement_value;
                    }
                }

                dual_rules[basis_element] = complement_value;
            }
            else {
                dual_rules[basis_element] = zero_str();
            }
        }
    }

    return dual_rules;
}

// Generate dual rules for odd-dimensional algebras: dual(u) = complement(G × u)
prd_rules generate_dual_rules(const AlgebraConfig& config, const prd_rules& complement_rules)
{
    const auto& basis = config.multivector_basis;
    const size_t n = basis.size();

    // Calculate extended metric matrix
    auto matrix_data = calculate_extended_metric_matrix(config);
    std::mdspan G{matrix_data.data(), n, n};

    prd_rules dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Perform matrix-vector multiplication: result = G × basis_element
        int metric_value = G[i, i];

        std::string intermediate_result;
        if (metric_value == 1) {
            intermediate_result = basis_element;
        }
        else if (metric_value == -1) {
            intermediate_result = minus_str() + basis_element;
        }
        else if (metric_value == 0) {
            intermediate_result = zero_str();
        }
        else {
            intermediate_result = std::to_string(metric_value) + mul_str() + basis_element;
        }

        // Apply complement to the intermediate result
        if (intermediate_result == zero_str()) {
            dual_rules[basis_element] = zero_str();
        }
        else {
            // Extract sign and basis element from intermediate result
            bool has_minus = (intermediate_result.find(minus_str()) == 0);
            std::string element_to_complement =
                has_minus ? intermediate_result.substr(1) : intermediate_result;

            // Look up complement
            auto it = complement_rules.find(element_to_complement);
            if (it != complement_rules.end()) {
                std::string complement_value = it->second;

                // Apply sign from intermediate result
                if (has_minus) {
                    if (complement_value.find(minus_str()) == 0) {
                        complement_value = complement_value.substr(1);
                    }
                    else if (complement_value != zero_str()) {
                        complement_value = minus_str() + complement_value;
                    }
                }

                dual_rules[basis_element] = complement_value;
            }
            else {
                dual_rules[basis_element] = zero_str();
            }
        }
    }

    return dual_rules;
}

////////////////////////////////////////////////////////////////////////////////
// PGA-Specific Dual Generation (Bulk and Weight Duals)
////////////////////////////////////////////////////////////////////////////////

// Generate bulk_dual rules for odd-dimensional PGA (uses complement)
// bulk_dual(u) = complement(G × u) where G is the extended metric matrix
prd_rules generate_bulk_dual_rules(const AlgebraConfig& config,
                                   const prd_rules& complement_rules)
{
    // For odd-dimensional PGA, bulk_dual is the same as regular dual
    return generate_dual_rules(config, complement_rules);
}

// Generate weight_dual rules for odd-dimensional PGA
// weight_dual(u) = complement(Ḡ × u) where Ḡ = regressive metric
// Regressive metric: Ḡ × u = complement(G × complement(u))
// Therefore: weight_dual(u) = complement(complement(G × complement(u))) = G × complement(u)
prd_rules generate_weight_dual_rules(const AlgebraConfig& config,
                                     const prd_rules& complement_rules,
                                     const prd_rules& left_complement_rules,
                                     const prd_rules& right_complement_rules)
{
    // For odd-dimensional PGA: weight_dual(u) = G × complement(u)
    const mvec_coeff& basis = config.multivector_basis;
    const size_t n = basis.size();
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
prd_rules generate_left_bulk_dual_rules(const AlgebraConfig& config,
                                        const prd_rules& left_complement_rules)
{
    return generate_left_dual_rules(config, left_complement_rules);
}

// Generate right_bulk_dual rules for even-dimensional PGA
prd_rules generate_right_bulk_dual_rules(const AlgebraConfig& config,
                                         const prd_rules& right_complement_rules)
{
    return generate_right_dual_rules(config, right_complement_rules);
}

// Generate left_weight_dual rules for even-dimensional PGA
// left_weight_dual(u) = left_complement(Ḡ × u) where Ḡ = regressive metric
// Regressive metric: Ḡ × u = left_complement(G × left_complement(u))
// Therefore: left_weight_dual(u) = left_complement(left_complement(G × left_complement(u)))
//                                  = G × left_complement(u) (using left_complement involution)
prd_rules generate_left_weight_dual_rules(const AlgebraConfig& config,
                                          const prd_rules& left_complement_rules,
                                          const prd_rules& right_complement_rules)
{
    // For even-dimensional PGA: left_weight_dual(u) = G × left_complement(u)
    const mvec_coeff& basis = config.multivector_basis;
    const size_t n = basis.size();
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

        // Step 2: Find the index of left_complement(u) in the basis (strip sign if present)
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

        // Step 3: Multiply by extended metric G[left_complement_index, left_complement_index]
        // Result is G × left_complement(u)
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
                left_weight_dual_rules[basis_element] = minus_str() + lcmpl_u; // add minus
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
// Therefore: right_weight_dual(u) = right_complement(left_complement(G × right_complement(u)))
//                                  = G × right_complement(u) (using complement involution)
prd_rules generate_right_weight_dual_rules(const AlgebraConfig& config,
                                           const prd_rules& left_complement_rules,
                                           const prd_rules& right_complement_rules)
{
    // For even-dimensional PGA: right_weight_dual(u) = G × right_complement(u)
    const mvec_coeff& basis = config.multivector_basis;
    const size_t n = basis.size();
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

        // Step 2: Find the index of right_complement(u) in the basis (strip sign if present)
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

        // Step 3: Multiply by extended metric G[right_complement_index, right_complement_index]
        // Result is G × right_complement(u)
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
                right_weight_dual_rules[basis_element] = minus_str() + rcmpl_u; // add minus
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

ProductRules generate_algebra_rules(const AlgebraConfig& config)
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
    // Even-dimensional algebras (EGA2D: 2D, STA3D: 4D, PGA3DP: 4D) have left and right complements
    // Odd-dimensional algebras (EGA3D: 3D, PGA2DP: 3D) have a single complement
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
            result.right_dual = generate_right_dual_rules(config, result.right_complement);
        }

        // For even-dimensional PGA, generate bulk and weight duals
        if (is_pga) {
            result.left_bulk_dual =
                generate_left_bulk_dual_rules(config, result.left_complement);
            result.right_bulk_dual =
                generate_right_bulk_dual_rules(config, result.right_complement);
            result.left_weight_dual = generate_left_weight_dual_rules(
                config, result.left_complement, result.right_complement);
            result.right_weight_dual = generate_right_weight_dual_rules(
                config, result.left_complement, result.right_complement);
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
            result.weight_dual = generate_weight_dual_rules(
                config, result.complement, result.left_complement, result.right_complement);
        }
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Validation and Testing Functions Implementation
////////////////////////////////////////////////////////////////////////////////

bool validate_rules(const prd_rules& generated, const prd_rules& reference)
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

void print_rule_comparison(const prd_rules& generated, const prd_rules& reference,
                           const std::string& product_name)
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