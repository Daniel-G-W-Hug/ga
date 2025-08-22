// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_rule_generator.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <set>

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
                    int swaps = j - i - 1;
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
    // Even-dimensional algebras (EGA2D: 2D, PGA3DP: 4D) have left and right complements
    // Odd-dimensional algebras (EGA3D: 3D, PGA2DP: 3D) have a single complement
    size_t num_basis_vectors = config.basis_vectors.size();
    bool is_even_dimensional = (num_basis_vectors % 2 == 0);

    if (is_even_dimensional) {
        // Generate both left and right complements for even algebras
        result.left_complement =
            generate_left_complement_rules(config, result.wedge_product);
        result.right_complement =
            generate_right_complement_rules(config, result.wedge_product);
    }
    else {
        // Generate single complement for odd algebras
        result.complement = generate_complement_rules(config, result.wedge_product);
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