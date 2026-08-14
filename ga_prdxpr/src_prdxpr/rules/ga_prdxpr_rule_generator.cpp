// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "rules/ga_prdxpr_rule_generator.hpp"
#include "rules/ga_prdxpr_dual_calc.hpp"
#include "rules/ga_prdxpr_metric_calc.hpp"
#include <algorithm>
#include <cmath>
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

// See declaration in ga_prdxpr_rule_generator.hpp.
int vector_metric_value(int digit, AlgebraConfig const& config)
{
    std::string const name = config.basis_prefix + std::to_string(digit);
    auto it = std::find(config.basis_vectors.begin(), config.basis_vectors.end(), name);
    if (it == config.basis_vectors.end()) return 0;
    size_t const slot =
        static_cast<size_t>(std::distance(config.basis_vectors.begin(), it));
    if (slot >= config.metric_signature.size()) return 0;
    return config.metric_signature[slot];
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
                    // Found a pair - apply metric signature via slot lookup.
                    // vector_metric_value() finds the basis vector by name in
                    // config.basis_vectors and indexes config.metric_signature
                    // at that slot, so reordering the basis vector list works
                    // correctly regardless of digit/slot alignment.
                    sign *= vector_metric_value(sequence[i], config);

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

    // Non-orthogonal metric: the orthogonal shortcuts below (product of vector
    // metrics; signature product for the pseudoscalar) are wrong on null pairs —
    // take the diagonal of the full Gram-determinant matrix instead.
    if (config.has_metric_matrix()) {
        auto const full = calculate_extended_metric_matrix_full(config);
        size_t const n = basis.size();
        for (size_t i = 0; i < n; ++i) {
            extended_metric[i] = full[i * n + i];
        }
        return extended_metric;
    }

    // Scalar always has metric value 1
    extended_metric[0] = 1;

    // Calculate determinant of metric for pseudoscalar
    int determinant = 1;
    for (int m : metric) {
        determinant *= m;
    }

    for (size_t i = 1; i < basis.size(); ++i) {
        std::string const& element = basis[i];
        auto digits = parse_indices(element, config.basis_prefix);

        if (digits.empty()) {
            extended_metric[i] = 0;
            continue;
        }

        // Pseudoscalar (highest grade): use determinant of metric.
        if (digits.size() == config.basis_vectors.size()) {
            extended_metric[i] = determinant;
            continue;
        }

        // Translate basis-element digits into slot indices via basis_vectors,
        // then accumulate the product-of-vector-metrics. Any unknown digit
        // (no matching basis vector) collapses the result to zero.
        int metric_value = 1;
        bool valid = true;
        for (int d : digits) {
            std::string const name = config.basis_prefix + std::to_string(d);
            auto it =
                std::find(config.basis_vectors.begin(), config.basis_vectors.end(), name);
            if (it == config.basis_vectors.end()) {
                metric_value = 0;
                valid = false;
                break;
            }
            int const slot =
                static_cast<int>(std::distance(config.basis_vectors.begin(), it));
            metric_value *= metric[slot];
        }

        if (!valid) {
            extended_metric[i] = 0;
            continue;
        }

        // The extended metric is the unique metric exomorphism: it extends the
        // vector metric via G(a∧b) = G(a)∧G(b), giving G(e_S) = ∏_{i∈S} g_i (the
        // product-of-vector-metrics computed above). This holds for every algebra
        // — EGA (all g_i=+1), PGA (degenerate g_i∈{0,+1}), and STA (mixed ±1). The
        // blade square e_S·e_S = σ(k)·∏g_i is the *geometric square*, a distinct
        // object that belongs to the causal/rotor path, not the metric slot.

        extended_metric[i] = metric_value;
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
    if (config.has_metric_matrix()) {
        // the single-term pair-cancellation in multiply_basis_elements assumes an
        // orthogonal basis and would silently produce wrong rules here
        throw std::runtime_error("generate_geometric_product_rules: non-orthogonal "
                                 "metric — use generate_geometric_product_rules_mt");
    }
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

    // Non-orthogonal metric: the dot product of two basis blades is the extended
    // metric matrix entry — off-diagonal pairs can be non-zero (metric partners,
    // e.g. a null-vector pair). Computed via the Gram determinant.
    if (config.has_metric_matrix()) {
        int const metric_value = compute_blade_inner_product(a, b, config);
        if (metric_value == 0) {
            return {zero_str(), 0};
        }
        if (metric_value != 1 && metric_value != -1) {
            // generate_ordered_rules encodes results as (element, sign) with
            // sign in {+1,-1}; larger magnitudes would be silently dropped
            throw std::runtime_error("multiply_basis_elements_dot: extended metric "
                                     "value with |value| > 1 is not representable in "
                                     "prd_rules for '" +
                                     a + "' * '" + b + "'");
        }
        return {config.scalar_name, metric_value};
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
// Multi-Term Geometric Product (non-orthogonal metrics) Implementation
////////////////////////////////////////////////////////////////////////////////

namespace {

// map a basis-vector digit to its slot in config.basis_vectors (-1 if unknown)
int digit_to_slot(int digit, AlgebraConfig const& config)
{
    std::string const name = config.basis_prefix + std::to_string(digit);
    auto it = std::find(config.basis_vectors.begin(), config.basis_vectors.end(), name);
    if (it == config.basis_vectors.end()) return -1;
    return static_cast<int>(std::distance(config.basis_vectors.begin(), it));
}

// map a slot in config.basis_vectors back to its digit (-1 if unparsable)
int slot_to_digit(int slot, AlgebraConfig const& config)
{
    auto const digits = parse_indices(config.basis_vectors[static_cast<size_t>(slot)],
                                      config.basis_prefix);
    if (digits.size() != 1) return -1;
    return digits[0];
}

// find the canonical basis element carrying exactly the given digit set (any
// order, no repeats) and the permutation sign relating the given order to the
// canonical listed order. Returns (position in multivector_basis, sign);
// position -1 if no canonical element matches. Empty digits -> scalar.
std::pair<int, int> canonical_from_digits(std::vector<int> const& digits,
                                          AlgebraConfig const& config)
{
    if (digits.empty()) {
        return {0, 1}; // scalar is always at position 0
    }

    for (size_t pos = 0; pos < config.multivector_basis.size(); ++pos) {
        std::string const& canonical = config.multivector_basis[pos];
        if (canonical == config.scalar_name) continue;

        auto canonical_digits = parse_indices(canonical, config.basis_prefix);
        if (canonical_digits.size() != digits.size()) continue;

        std::vector<int> sorted_given = digits;
        std::vector<int> sorted_canonical = canonical_digits;
        std::sort(sorted_given.begin(), sorted_given.end());
        std::sort(sorted_canonical.begin(), sorted_canonical.end());
        if (sorted_given != sorted_canonical) continue;

        // count swaps needed to transform digits into canonical_digits
        std::vector<int> temp = digits;
        int swaps = 0;
        for (size_t i = 0; i < canonical_digits.size(); ++i) {
            auto it = std::find(temp.begin() + static_cast<long>(i), temp.end(),
                                canonical_digits[i]);
            if (it != temp.begin() + static_cast<long>(i)) {
                std::swap(temp[i], *it);
                ++swaps;
            }
        }
        int const sign = (swaps % 2 == 0) ? 1 : -1;
        return {static_cast<int>(pos), sign};
    }
    return {-1, 0};
}

// collect (basis position -> integer coefficient) into an ordered prd_terms
prd_terms terms_from_position_map(std::map<int, int> const& by_position,
                                  AlgebraConfig const& config)
{
    prd_terms terms;
    for (auto const& [pos, coeff] : by_position) {
        if (coeff == 0) continue;
        terms.push_back(
            prd_term{coeff, config.multivector_basis[static_cast<size_t>(pos)]});
    }
    return terms;
}

// wedge-multilinear expansion: expand the wedge of the given source vectors
// (rows of M, indexed by source slot, over target slots) into canonicalized
// (sorted, sign-folded) target index sequences with double coefficients.
// M[s][m] = coefficient of target vector m in source vector s.
std::map<std::vector<int>, double> wedge_expand(std::vector<int> const& source_slots,
                                                std::vector<std::vector<double>> const& M)
{
    std::vector<std::pair<double, std::vector<int>>> terms{{1.0, {}}};

    for (int s : source_slots) {
        auto const& row = M[static_cast<size_t>(s)];
        std::vector<std::pair<double, std::vector<int>>> next;
        for (auto const& [c, seq] : terms) {
            for (size_t m = 0; m < row.size(); ++m) {
                if (row[m] == 0.0) continue;
                if (std::find(seq.begin(), seq.end(), static_cast<int>(m)) != seq.end()) {
                    continue; // wedge: repeated vector annihilates the term
                }
                auto seq2 = seq;
                seq2.push_back(static_cast<int>(m));
                next.emplace_back(c * row[m], seq2);
            }
        }
        terms = std::move(next);
    }

    // canonicalize: sort each sequence ascending, folding the permutation sign
    std::map<std::vector<int>, double> acc;
    for (auto& [c, seq] : terms) {
        int sign = 1;
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            for (size_t j = 0; j + 1 < seq.size() - i; ++j) {
                if (seq[j] > seq[j + 1]) {
                    std::swap(seq[j], seq[j + 1]);
                    sign = -sign;
                }
            }
        }
        acc[seq] += sign * c;
    }
    std::erase_if(acc, [](auto const& kv) { return kv.second == 0.0; });
    return acc;
}

// geometric product of two ORTHOGONAL-basis blades given as strictly sorted
// slot sequences; sig holds the diagonal signature per slot. Returns
// (integer factor, sorted merged sequence); factor 0 means annihilation.
std::pair<int, std::vector<int>> diag_blade_gpr(std::vector<int> const& a,
                                                std::vector<int> const& b,
                                                std::vector<int> const& sig)
{
    std::vector<int> seq = a;
    seq.insert(seq.end(), b.begin(), b.end());
    int sign = 1;

    // cancel repeated slots (each appears at most twice) with the metric factor
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < seq.size() && !changed; ++i) {
            for (size_t j = i + 1; j < seq.size(); ++j) {
                if (seq[i] == seq[j]) {
                    sign *= sig[static_cast<size_t>(seq[i])];
                    if ((j - i - 1) % 2 == 1) sign = -sign;
                    seq.erase(seq.begin() + static_cast<long>(j));
                    seq.erase(seq.begin() + static_cast<long>(i));
                    changed = true;
                    break;
                }
            }
        }
        if (sign == 0) return {0, {}};
    }

    // sort the remaining distinct slots, folding the permutation sign
    for (size_t i = 0; i + 1 < seq.size(); ++i) {
        for (size_t j = 0; j + 1 < seq.size() - i; ++j) {
            if (seq[j] > seq[j + 1]) {
                std::swap(seq[j], seq[j + 1]);
                sign = -sign;
            }
        }
    }
    return {sign, seq};
}

// parse a basis element into basis-vector slots (throws on unknown digits)
std::vector<int> blade_slots(std::string const& blade, AlgebraConfig const& config)
{
    std::vector<int> slots;
    for (int d : parse_indices(blade, config.basis_prefix)) {
        int const s = digit_to_slot(d, config);
        if (s < 0) {
            throw std::runtime_error("unknown basis vector digit in blade '" + blade +
                                     "'");
        }
        slots.push_back(s);
    }
    return slots;
}

} // namespace

void validate_metric_matrix_config(AlgebraConfig const& config)
{
    if (!config.has_metric_matrix()) return;

    size_t const n = config.basis_vectors.size();
    auto const& G = config.metric_matrix;
    auto const& T = config.basis_change;
    auto const& Tinv = config.basis_change_inv;
    auto const& sig = config.diag_signature;

    auto require = [](bool cond, std::string const& msg) {
        if (!cond) throw std::runtime_error("metric matrix config invalid: " + msg);
    };

    require(G.size() == n, "metric_matrix must be n x n (row count)");
    for (auto const& row : G) {
        require(row.size() == n, "metric_matrix must be n x n (column count)");
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            require(G[i][j] == G[j][i], "metric_matrix must be symmetric");
        }
    }
    require(config.metric_signature.size() == n,
            "metric_signature must have one entry per basis vector");
    for (size_t i = 0; i < n; ++i) {
        require(config.metric_signature[i] == G[i][i],
                "metric_signature must equal the metric_matrix diagonal");
    }

    // Non-degeneracy, checked EXACTLY via the integer cofactor determinant
    // (compute_integer_determinant, n <= 5). The PGA (degenerate-metric)
    // classification is signature-based and deliberately skipped for matrix
    // configs — a degenerate full-matrix metric would silently take the wrong
    // dual path, so it is rejected here: degenerate algebras must be
    // configured via the orthogonal metric_signature path instead.
    {
        std::vector<int> G_flat(n * n, 0);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                G_flat[i * n + j] = G[i][j];
            }
        }
        std::mdspan<int const,
                    std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const
            G_view{G_flat.data(), n, n};
        require(compute_integer_determinant(G_view, static_cast<int>(n)) != 0,
                "metric_matrix is degenerate (det = 0) — degenerate algebras must use "
                "the orthogonal metric_signature path");
    }

    require(T.size() == n && Tinv.size() == n && sig.size() == n,
            "basis_change, basis_change_inv and diag_signature must be provided "
            "(n rows each) together with metric_matrix");
    for (auto const& row : T) {
        require(row.size() == n, "basis_change must be n x n");
    }
    for (auto const& row : Tinv) {
        require(row.size() == n, "basis_change_inv must be n x n");
    }

    // basis_change * basis_change_inv = I (dyadic rationals: exact in double)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < n; ++k) {
                sum += T[i][k] * Tinv[k][j];
            }
            require(sum == (i == j ? 1.0 : 0.0),
                    "basis_change * basis_change_inv must be the identity");
        }
    }

    // metric_matrix = basis_change * diag(diag_signature) * basis_change^T
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < n; ++k) {
                sum += T[i][k] * sig[k] * T[j][k];
            }
            require(sum == static_cast<double>(G[i][j]),
                    "metric_matrix must equal basis_change * diag(diag_signature) "
                    "* basis_change^T");
        }
    }
}

// A linear combination of strictly ascending slot SEQUENCES, each denoting a
// GEOMETRIC PRODUCT of basis vectors (NOT a wedge blade — for non-orthogonal
// metrics the two differ from grade 2 upward, e.g. e_w e_u = -1 + e_w∧e_u).
// Coefficients are dyadic rationals, exact in double.
using seq_combo = std::map<std::vector<int>, double>;

namespace {

// accumulate coeff * (normal-ordered form of the vector product sequence seq)
// into acc, using the non-orthogonal rewrite rules
//
//     e_i e_i = G[i,i]
//     e_i e_j = 2*G[i,j] - e_j e_i        (i != j)
//
// the worklist shrinks in (length, inversions), so the rewriting terminates
void normal_order_into(seq_combo& acc, double coeff, std::vector<int> seq,
                       AlgebraConfig const& config)
{
    auto const& G = config.metric_matrix;

    struct term_t {
        double c;
        std::vector<int> s;
    };
    std::vector<term_t> worklist{{coeff, std::move(seq)}};

    while (!worklist.empty()) {
        term_t t = std::move(worklist.back());
        worklist.pop_back();

        // find first adjacent violation of strict ascending order
        size_t p = 0;
        bool violation = false;
        for (; p + 1 < t.s.size(); ++p) {
            if (t.s[p] >= t.s[p + 1]) {
                violation = true;
                break;
            }
        }

        if (!violation) {
            acc[t.s] += t.c;
            continue;
        }

        int const i = t.s[p];
        int const j = t.s[p + 1];
        if (i == j) {
            // contract with the diagonal metric entry
            int const g = G[static_cast<size_t>(i)][static_cast<size_t>(i)];
            if (g != 0) {
                term_t t2{t.c * g, t.s};
                t2.s.erase(t2.s.begin() + static_cast<long>(p),
                           t2.s.begin() + static_cast<long>(p) + 2);
                worklist.push_back(std::move(t2));
            }
        }
        else {
            // swap term: -(e_j e_i)
            term_t t1{-t.c, t.s};
            std::swap(t1.s[p], t1.s[p + 1]);
            worklist.push_back(std::move(t1));

            // contraction term: 2*G[i,j] with both slots removed
            int const g = G[static_cast<size_t>(i)][static_cast<size_t>(j)];
            if (g != 0) {
                term_t t2{t.c * 2 * g, t.s};
                t2.s.erase(t2.s.begin() + static_cast<long>(p),
                           t2.s.begin() + static_cast<long>(p) + 2);
                worklist.push_back(std::move(t2));
            }
        }
    }
}

// express the wedge BLADE over the strictly ascending slot set as a combination
// of ordered product sequences, via the Riesz recursion
//
//     blade(S) = v ∧ blade(S\{v})
//     v ∧ A_k  = 1/2 * (v ⟑ A_k + (-1)^k A_k ⟑ v)        [A_k of grade k]
//
// unitriangular by grade: blade(S) = seq(S) + (strictly lower-grade sequences),
// which makes the inverse conversion a simple top-down elimination
seq_combo blade_to_seq_combo(std::vector<int> const& slots, AlgebraConfig const& config)
{
    seq_combo res;
    if (slots.size() <= 1) {
        res[slots] = 1.0;
        return res;
    }

    int const v = slots.front();
    std::vector<int> const rest(slots.begin() + 1, slots.end());
    auto const inner = blade_to_seq_combo(rest, config);

    // (-1)^k with k = grade of blade(rest) = slots.size() - 1
    double const sign = ((slots.size() - 1) % 2 == 0) ? 1.0 : -1.0;

    for (auto const& [s, c] : inner) {
        std::vector<int> vs{v};
        vs.insert(vs.end(), s.begin(), s.end());
        normal_order_into(res, 0.5 * c, std::move(vs), config);

        std::vector<int> sv{s};
        sv.push_back(v);
        normal_order_into(res, 0.5 * c * sign, std::move(sv), config);
    }
    std::erase_if(res, [](auto const& kv) { return kv.second == 0.0; });
    return res;
}

} // namespace

prd_terms multiply_basis_elements_general(std::string const& a, std::string const& b,
                                          AlgebraConfig const& config)
{
    // scalar shortcuts
    if (a == config.scalar_name && b == config.scalar_name) {
        return {prd_term{1, config.scalar_name}};
    }
    if (a == config.scalar_name) return {prd_term{1, b}};
    if (b == config.scalar_name) return {prd_term{1, a}};

    // sort each blade's slots ascending, folding the permutation sign (the
    // blade named in the basis may list its indices in any order)
    auto sorted_slots_with_sign = [&](std::string const& blade) {
        auto slots = blade_slots(blade, config);
        double sign = 1.0;
        for (size_t i = 0; i + 1 < slots.size(); ++i) {
            for (size_t j = 0; j + 1 < slots.size() - i; ++j) {
                if (slots[j] > slots[j + 1]) {
                    std::swap(slots[j], slots[j + 1]);
                    sign = -sign;
                }
            }
        }
        return std::pair{slots, sign};
    };

    auto const [slots_a, sign_a] = sorted_slots_with_sign(a);
    auto const [slots_b, sign_b] = sorted_slots_with_sign(b);

    // blade -> product-sequence combos, multiply pairwise, normal-order
    auto const combo_a = blade_to_seq_combo(slots_a, config);
    auto const combo_b = blade_to_seq_combo(slots_b, config);

    seq_combo prod;
    for (auto const& [sa, ca] : combo_a) {
        for (auto const& [sb, cb] : combo_b) {
            std::vector<int> seq{sa};
            seq.insert(seq.end(), sb.begin(), sb.end());
            normal_order_into(prod, sign_a * sign_b * ca * cb, std::move(seq), config);
        }
    }
    std::erase_if(prod, [](auto const& kv) { return kv.second == 0.0; });

    // convert the product-sequence combo back to blade coefficients by top-down
    // elimination (longest sequences first; blade(U) = seq(U) + lower terms)
    std::map<std::vector<int>, double> blades; // ascending slot set -> coefficient
    while (!prod.empty()) {
        auto longest = prod.begin();
        for (auto it = prod.begin(); it != prod.end(); ++it) {
            if (it->first.size() > longest->first.size()) longest = it;
        }
        auto const seq = longest->first;
        double const c = longest->second;

        auto const bc = blade_to_seq_combo(seq, config);

        // unitriangularity invariant: blade(U) = 1*seq(U) + lower-grade terms.
        // If it fails, the subtraction below would not eliminate seq(U) and the
        // loop would never terminate — fail loudly instead (guards against
        // implementation bugs in the blade/sequence conversion).
        auto const self = bc.find(seq);
        if (self == bc.end() || self->second != 1.0) {
            throw std::runtime_error("multiply_basis_elements_general: blade/sequence "
                                     "conversion lost unitriangularity for '" +
                                     a + "' * '" + b + "'");
        }

        blades[seq] += c;
        for (auto const& [s, sc] : bc) {
            prod[s] -= c * sc;
        }
        std::erase_if(prod, [](auto const& kv) { return kv.second == 0.0; });
    }

    // map onto canonical basis elements; coefficients must come out integral
    std::map<int, int> by_position;
    for (auto const& [seq, c] : blades) {
        if (c == 0.0) continue;
        double const rounded = std::round(c);
        if (c != rounded) {
            throw std::runtime_error("multiply_basis_elements_general: non-integer "
                                     "coefficient for a product term of '" +
                                     a + "' * '" + b + "'");
        }
        std::vector<int> digits;
        for (int s : seq) {
            digits.push_back(slot_to_digit(s, config));
        }
        auto const [pos, sign] = canonical_from_digits(digits, config);
        if (pos < 0) {
            throw std::runtime_error("multiply_basis_elements_general: no canonical "
                                     "basis element for a product term of '" +
                                     a + "' * '" + b + "'");
        }
        by_position[pos] += static_cast<int>(rounded) * sign;
    }
    return terms_from_position_map(by_position, config);
}

prd_terms multiply_basis_elements_via_diagonal(std::string const& a, std::string const& b,
                                               AlgebraConfig const& config)
{
    // scalar shortcuts
    if (a == config.scalar_name && b == config.scalar_name) {
        return {prd_term{1, config.scalar_name}};
    }
    if (a == config.scalar_name) return {prd_term{1, b}};
    if (b == config.scalar_name) return {prd_term{1, a}};

    // expand both blades in the internal orthogonal (diagonal) basis
    auto const ea = wedge_expand(blade_slots(a, config), config.basis_change);
    auto const eb = wedge_expand(blade_slots(b, config), config.basis_change);

    // multiply in the diagonal basis (single-term per blade pair)
    std::map<std::vector<int>, double> res_diag;
    for (auto const& [seq_a, ca] : ea) {
        for (auto const& [seq_b, cb] : eb) {
            auto const [factor, seq] =
                diag_blade_gpr(seq_a, seq_b, config.diag_signature);
            if (factor != 0) {
                res_diag[seq] += ca * cb * factor;
            }
        }
    }

    // transform the result back into the non-orthogonal basis
    std::map<std::vector<int>, double> res_null;
    for (auto const& [seq, c] : res_diag) {
        if (c == 0.0) continue;
        for (auto const& [nseq, nc] : wedge_expand(seq, config.basis_change_inv)) {
            res_null[nseq] += c * nc;
        }
    }

    // map onto canonical basis elements; coefficients must come out integral
    // (all arithmetic is dyadic-rational and therefore exact in double)
    std::map<int, int> by_position;
    for (auto const& [seq, c] : res_null) {
        if (c == 0.0) continue;
        double const rounded = std::round(c);
        if (c != rounded) {
            throw std::runtime_error(
                "multiply_basis_elements_via_diagonal: non-integer coefficient for "
                "a product term of '" +
                a + "' * '" + b + "'");
        }
        std::vector<int> digits;
        for (int s : seq) {
            digits.push_back(slot_to_digit(s, config));
        }
        auto const [pos, sign] = canonical_from_digits(digits, config);
        if (pos < 0) {
            throw std::runtime_error("multiply_basis_elements_via_diagonal: no "
                                     "canonical basis element for a product term of '" +
                                     a + "' * '" + b + "'");
        }
        by_position[pos] += static_cast<int>(rounded) * sign;
    }
    return terms_from_position_map(by_position, config);
}

prd_rules_mt generate_geometric_product_rules_mt(AlgebraConfig const& config)
{
    if (!config.has_metric_matrix()) {
        throw std::runtime_error("generate_geometric_product_rules_mt requires a "
                                 "metric_matrix (use generate_geometric_product_rules "
                                 "for orthogonal metrics)");
    }
    validate_metric_matrix_config(config);

    prd_rules_mt rules;
    for (auto const& a : config.multivector_basis) {
        for (auto const& b : config.multivector_basis) {
            // compute each pair via BOTH routes: the diagonal-basis detour is the
            // primary computation; the direct non-orthogonal rewrite must agree
            // (always-on transcription gate for the generator itself)
            auto const detour = multiply_basis_elements_via_diagonal(a, b, config);
            auto const direct = multiply_basis_elements_general(a, b, config);
            if (detour != direct) {
                throw std::runtime_error(fmt::format(
                    "generate_geometric_product_rules_mt: diagonal-detour "
                    "and direct-rewrite results disagree for '{} * {}': "
                    "'{}' vs '{}'",
                    a, b, prd_terms_to_string(detour), prd_terms_to_string(direct)));
            }
            rules[a + space_str() + mul_str() + space_str() + b] = detour;
        }
    }
    return rules;
}

////////////////////////////////////////////////////////////////////////////////
// Complement Rule Generation Implementation
////////////////////////////////////////////////////////////////////////////////

// Generate complement rules from wedge product table
// Algorithm: For complement relationship u ^ cmpl(u) = I_n (right complement)
//            or cmpl(u) ^ u = I_n (left complement)
prd_rules generate_cmpl_from_wedge_table(AlgebraConfig const& config,
                                         prd_rules const& wedge_rules, bool is_l_cmpl)
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
        bool found_cmpl = false;

        // Special case: scalar complement is always the pseudoscalar
        if (basis_element == config.scalar_name) {
            complement_rules[basis_element] = pseudoscalar;
            found_cmpl = true;
        }
        // Special case: pseudoscalar complement is always the scalar
        else if (basis_element == pseudoscalar) {
            complement_rules[basis_element] = config.scalar_name;
            found_cmpl = true;
        }
        else {
            if (is_l_cmpl) {
                // Left complement: l_cmpl(u) ^ u = I_n
                // Search column i for pseudoscalar
                for (size_t row = 0; row < basis_size; ++row) {
                    if (row < wedge_table.size() && i < wedge_table[row].size()) {
                        const std::string& table_entry = wedge_table[row][i];
                        if (table_entry == pseudoscalar) {
                            complement_rules[basis_element] =
                                config.multivector_basis[row];
                            found_cmpl = true;
                            break;
                        }
                        else if (table_entry == "-" + pseudoscalar) {
                            complement_rules[basis_element] =
                                "-" + config.multivector_basis[row];
                            found_cmpl = true;
                            break;
                        }
                    }
                }
            }
            else {
                // Right complement: u ^ r_cmpl(u) = I_n
                // Search row i for pseudoscalar
                if (i < wedge_table.size()) {
                    for (size_t col = 0; col < basis_size && col < wedge_table[i].size();
                         ++col) {
                        const std::string& table_entry = wedge_table[i][col];
                        if (table_entry == pseudoscalar) {
                            complement_rules[basis_element] =
                                config.multivector_basis[col];
                            found_cmpl = true;
                            break;
                        }
                        else if (table_entry == "-" + pseudoscalar) {
                            complement_rules[basis_element] =
                                "-" + config.multivector_basis[col];
                            found_cmpl = true;
                            break;
                        }
                    }
                }
            }
        }

        if (!found_cmpl) {
            throw std::runtime_error("Could not find complement for basis element: " +
                                     basis_element);
        }
    }

    return complement_rules;
}

prd_rules generate_r_cmpl_rules(AlgebraConfig const& config, prd_rules const& wedge_rules)
{
    return generate_cmpl_from_wedge_table(config, wedge_rules, false);
}

prd_rules generate_l_cmpl_rules(AlgebraConfig const& config, prd_rules const& wedge_rules)
{
    return generate_cmpl_from_wedge_table(config, wedge_rules, true);
}

prd_rules generate_cmpl_rules(AlgebraConfig const& config, prd_rules const& wedge_rules)
{
    // For odd algebras, left and right complements are the same
    return generate_cmpl_from_wedge_table(config, wedge_rules, false);
}

////////////////////////////////////////////////////////////////////////////////
// Dual Rule Generation Implementation
////////////////////////////////////////////////////////////////////////////////

// Generate left dual rules: l_dual(u) = l_cmpl(G × u)
// where G is the extended metric matrix and × is matrix-vector multiplication
prd_rules generate_l_dual_rules(AlgebraConfig const& config,
                                prd_rules const& l_cmpl_rules)
{
    // Delegate to systematic matrix-vector multiplication approach
    // l_dual(basis[i]) = Σⱼ G[i,j] · l_cmpl(basis[j])
    auto G_data = calculate_extended_metric_matrix_full(config);
    return calculate_dual_rules(config, G_data, l_cmpl_rules);
}

// Generate right dual rules: r_dual(u) = r_cmpl(G × u)
prd_rules generate_r_dual_rules(AlgebraConfig const& config,
                                prd_rules const& r_cmpl_rules)
{
    // Delegate to systematic matrix-vector multiplication approach
    // r_dual(basis[i]) = Σⱼ G[i,j] · r_cmpl(basis[j])
    auto G_data = calculate_extended_metric_matrix_full(config);
    return calculate_dual_rules(config, G_data, r_cmpl_rules);
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

// Generate l_bulk_dual rules for even-dimensional PGA
prd_rules generate_l_bulk_dual_rules(AlgebraConfig const& config,
                                     prd_rules const& l_cmpl_rules)
{
    return generate_l_dual_rules(config, l_cmpl_rules);
}

// Generate r_bulk_dual rules for even-dimensional PGA
prd_rules generate_r_bulk_dual_rules(AlgebraConfig const& config,
                                     prd_rules const& r_cmpl_rules)
{
    return generate_r_dual_rules(config, r_cmpl_rules);
}

// Generate l_weight_dual rules for even-dimensional PGA
// l_weight_dual(u) = l_cmpl(Ḡ × u) where Ḡ = regressive metric
// Regressive metric: Ḡ × u = l_cmpl(G × l_cmpl(u))
// Therefore: l_weight_dual(u) = l_cmpl(l_cmpl(G ×
// l_cmpl(u)))
//                                  = G × l_cmpl(u) (using l_cmpl
//                                  involution)
prd_rules generate_l_weight_dual_rules(AlgebraConfig const& config,
                                       prd_rules const& l_cmpl_rules)
{
    // For even-dimensional PGA: l_weight_dual(u) = G × l_cmpl(u)
    mvec_coeff const& basis = config.multivector_basis;
    size_t const n = basis.size();
    auto G_data = calculate_extended_metric_matrix(config);
    std::mdspan G{G_data.data(), n, n};

    prd_rules l_weight_dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Step 1: Apply l_cmpl to input basis element
        auto it_lcmpl = l_cmpl_rules.find(basis_element);
        if (it_lcmpl == l_cmpl_rules.end()) {
            l_weight_dual_rules[basis_element] = zero_str();
            continue;
        }
        std::string l_cmpl_u = it_lcmpl->second;

        if (l_cmpl_u == zero_str()) {
            l_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 2: Find the index of l_cmpl(u) in the basis (strip sign if
        // present)
        bool l_cmpl_has_minus = (l_cmpl_u.find(minus_str()) == 0);
        std::string l_cmpl_u_unsigned = l_cmpl_has_minus ? l_cmpl_u.substr(1) : l_cmpl_u;

        size_t l_cmpl_index = n; // invalid index
        for (size_t j = 0; j < n; ++j) {
            if (basis[j] == l_cmpl_u_unsigned) {
                l_cmpl_index = j;
                break;
            }
        }

        if (l_cmpl_index >= n) {
            l_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 3: Multiply by extended metric G[l_cmpl_index,
        // l_cmpl_index] Result is G × l_cmpl(u)
        int metric_value = G[l_cmpl_index, l_cmpl_index];

        if (metric_value == 1) {
            // metric × l_cmpl = l_cmpl (preserve sign)
            l_weight_dual_rules[basis_element] = l_cmpl_u;
        }
        else if (metric_value == -1) {
            // metric × l_cmpl = -l_cmpl (flip sign)
            if (l_cmpl_has_minus) {
                l_weight_dual_rules[basis_element] = l_cmpl_u_unsigned; // remove minus
            }
            else {
                l_weight_dual_rules[basis_element] = minus_str() + l_cmpl_u; // add minus
            }
        }
        else if (metric_value == 0) {
            l_weight_dual_rules[basis_element] = zero_str();
        }
        else {
            // General case: metric_value × l_cmpl
            l_weight_dual_rules[basis_element] =
                std::to_string(metric_value) + mul_str() + l_cmpl_u;
        }
    }

    return l_weight_dual_rules;
}

// Generate r_weight_dual rules for even-dimensional PGA
// r_weight_dual(u) = r_cmpl(Ḡ × u) where Ḡ = regressive metric
// Regressive metric: Ḡ × u = l_cmpl(G × r_cmpl(u))
// Therefore: r_weight_dual(u) = r_cmpl(l_cmpl(G ×
// r_cmpl(u)))
//                                  = G × r_cmpl(u) (using complement
//                                  involution)
prd_rules generate_r_weight_dual_rules(AlgebraConfig const& config,
                                       prd_rules const& r_cmpl_rules)
{
    // For even-dimensional PGA: r_weight_dual(u) = G × r_cmpl(u)
    mvec_coeff const& basis = config.multivector_basis;
    size_t const n = basis.size();
    auto G_data = calculate_extended_metric_matrix(config);
    std::mdspan G{G_data.data(), n, n};

    prd_rules r_weight_dual_rules;

    // For each basis element
    for (size_t i = 0; i < n; ++i) {
        const std::string& basis_element = basis[i];

        // Step 1: Apply r_cmpl to input basis element
        auto it_rcmpl = r_cmpl_rules.find(basis_element);
        if (it_rcmpl == r_cmpl_rules.end()) {
            r_weight_dual_rules[basis_element] = zero_str();
            continue;
        }
        std::string r_cmpl_u = it_rcmpl->second;

        if (r_cmpl_u == zero_str()) {
            r_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 2: Find the index of r_cmpl(u) in the basis (strip sign if
        // present)
        bool r_cmpl_has_minus = (r_cmpl_u.find(minus_str()) == 0);
        std::string r_cmpl_u_unsigned = r_cmpl_has_minus ? r_cmpl_u.substr(1) : r_cmpl_u;

        size_t r_cmpl_index = n; // invalid index
        for (size_t j = 0; j < n; ++j) {
            if (basis[j] == r_cmpl_u_unsigned) {
                r_cmpl_index = j;
                break;
            }
        }

        if (r_cmpl_index >= n) {
            r_weight_dual_rules[basis_element] = zero_str();
            continue;
        }

        // Step 3: Multiply by extended metric G[r_cmpl_index,
        // r_cmpl_index] Result is G × r_cmpl(u)
        int metric_value = G[r_cmpl_index, r_cmpl_index];

        if (metric_value == 1) {
            // metric × r_cmpl = r_cmpl (preserve sign)
            r_weight_dual_rules[basis_element] = r_cmpl_u;
        }
        else if (metric_value == -1) {
            // metric × r_cmpl = -r_cmpl (flip sign)
            if (r_cmpl_has_minus) {
                r_weight_dual_rules[basis_element] = r_cmpl_u_unsigned; // remove minus
            }
            else {
                r_weight_dual_rules[basis_element] = minus_str() + r_cmpl_u; // add minus
            }
        }
        else if (metric_value == 0) {
            r_weight_dual_rules[basis_element] = zero_str();
        }
        else {
            // General case: metric_value × r_cmpl
            r_weight_dual_rules[basis_element] =
                std::to_string(metric_value) + mul_str() + r_cmpl_u;
        }
    }

    return r_weight_dual_rules;
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
    if (config.has_metric_matrix()) {
        // non-orthogonal metric: the geometric product of basis blades is
        // multi-term — generated (and cross-checked) via the mt path; the
        // single-term slot stays empty
        result.geometric_product_mt = generate_geometric_product_rules_mt(config);
    }
    else {
        result.geometric_product = generate_geometric_product_rules(config);
    }
    result.wedge_product = generate_wedge_product_rules(config);
    result.dot_product = generate_dot_product_rules(config);

    // Generate complement rules based on algebra dimensionality
    // Even-dimensional algebras (EGA2D: 2D, STA4D: 4D, PGA3DP: 4D) have left and right
    // complements Odd-dimensional algebras (EGA3D: 3D, PGA2DP: 3D) have a single
    // complement
    size_t num_basis_vectors = config.basis_vectors.size();
    bool is_even_dimensional = (num_basis_vectors % 2 == 0);

    // Detect PGA algebras (have degenerate dimensions with metric signature 0).
    // Only meaningful for orthogonal metrics: a full metric matrix may carry
    // zeros on its diagonal (null vectors) while being non-degenerate overall
    // (e.g. a conformal null pair, det != 0) — those algebras get the regular
    // (non-PGA) dual path. validate_metric_matrix_config() REJECTS degenerate
    // full-matrix metrics (exact integer determinant), so this branch cannot
    // misclassify one.
    bool is_pga = false;
    if (!config.has_metric_matrix()) {
        for (int metric_val : config.metric_signature) {
            if (metric_val == 0) {
                is_pga = true;
                break;
            }
        }
    }

    // Antidual support for non-orthogonal (matrix) algebras: the antidual is
    // the dual w.r.t. the metric anti-exomorphism 𝔾 = C·G·Cᵀ. For the
    // conformal null-basis metrics G is a symmetric signed permutation with
    // G·G = I and det G = -1, which makes 𝔾 = -G and hence
    // antidual(u) = -dual(u). The identity is VERIFIED here (throws on
    // violation) rather than assumed, so the shortcut cannot silently be
    // applied to an algebra where it does not hold.
    auto verify_antimetric_is_negated_metric = [&](prd_rules const& cmpl_rules) {
        auto const G_ext = calculate_extended_metric_matrix_full(config);
        auto const G_anti =
            calculate_regressive_extended_metric_matrix_full(config, cmpl_rules);
        if (G_ext.size() != G_anti.size()) {
            throw std::runtime_error("antidual generation: metric matrix size mismatch");
        }
        for (size_t k = 0; k < G_ext.size(); ++k) {
            if (G_anti[k] != -G_ext[k]) {
                throw std::runtime_error(
                    "antidual generation: anti-exomorphism is not the negated "
                    "extended metric for this algebra -- the antidual(u) = -dual(u) "
                    "shortcut does not apply");
            }
        }
    };
    auto negate_rules = [](prd_rules const& rules) {
        prd_rules negated;
        for (auto const& [key, value] : rules) {
            if (value == zero_str()) {
                negated[key] = value;
            }
            else if (value.starts_with(minus_str())) {
                negated[key] = value.substr(1, value.size() - 1);
            }
            else {
                negated[key] = minus_str() + value;
            }
        }
        return negated;
    };

    if (is_even_dimensional) {
        // Generate both left and right complements for even algebras
        result.l_cmpl = generate_l_cmpl_rules(config, result.wedge_product);
        result.r_cmpl = generate_r_cmpl_rules(config, result.wedge_product);

        // Generate left and right dual rules from complements (only for non-PGA)
        if (!is_pga) {
            result.l_dual = generate_l_dual_rules(config, result.l_cmpl);
            result.r_dual = generate_r_dual_rules(config, result.r_cmpl);
        }

        // Antiduals for non-orthogonal metrics (verified 𝔾 = -G, see above)
        if (config.has_metric_matrix()) {
            verify_antimetric_is_negated_metric(result.r_cmpl);
            result.l_antidual = negate_rules(result.l_dual);
            result.r_antidual = negate_rules(result.r_dual);
        }

        // For even-dimensional PGA, generate bulk and weight duals
        if (is_pga) {
            result.l_bulk_dual = generate_l_bulk_dual_rules(config, result.l_cmpl);
            result.r_bulk_dual = generate_r_bulk_dual_rules(config, result.r_cmpl);
            result.l_weight_dual = generate_l_weight_dual_rules(config, result.l_cmpl);
            result.r_weight_dual = generate_r_weight_dual_rules(config, result.r_cmpl);
        }
    }
    else {
        // Generate single complement for odd algebras
        result.complement = generate_cmpl_rules(config, result.wedge_product);

        // Generate dual rules from complement (only for non-PGA)
        if (!is_pga) {
            result.dual = generate_dual_rules(config, result.complement);
        }

        // Antidual for non-orthogonal metrics (verified 𝔾 = -G, see above)
        if (config.has_metric_matrix()) {
            verify_antimetric_is_negated_metric(result.complement);
            result.antidual = negate_rules(result.dual);
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