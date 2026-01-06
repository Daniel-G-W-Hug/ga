// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_prdxpr_ega2d.hpp"
#include "ga_prdxpr_ega3d.hpp"
#include "ga_prdxpr_pga2dp.hpp"
#include "ga_prdxpr_pga3dp.hpp"
#include "ga_prdxpr_rule_generator.hpp"
#include "ga_prdxpr_sta4d.hpp"

#include <iostream>
#include <mdspan>

// Helper function to print all rules in grade-ordered way
void print_all_rules(const prd_rules& rules, const std::string& title,
                     const mvec_coeff& basis_order, const std::string& algebra_name = "")
{
    if (!algebra_name.empty()) {
        fmt::println("\n{} {}", algebra_name, title);
    }
    else {
        fmt::println("\n=== {} ===", title);
    }

    // Print in grade order by iterating through basis in order
    for (const auto& a : basis_order) {
        for (const auto& b : basis_order) {
            std::string key = a + " * " + b;
            if (title.find("wedge") != std::string::npos) {
                key = a + " ^ " + b; // Use wedge operator for wedge products
            }

            auto it = rules.find(key);
            if (it != rules.end()) {
                fmt::println("    {{\"{}\", \"{}\"}},", key, it->second);
            }
        }
    }
}

// Print complement rules in grade order (scalars, vectors, bivectors, etc.)
void print_complement_rules(const prd_rules& rules, const std::string& title,
                            const mvec_coeff& basis_order,
                            const std::string& algebra_name = "")
{
    if (!algebra_name.empty()) {
        fmt::println("\n{} {}", algebra_name, title);
    }
    else {
        fmt::println("\n=== {} ===", title);
    }

    // Print complement rules in grade order by iterating through basis in order
    for (const auto& basis_element : basis_order) {
        auto it = rules.find(basis_element);
        if (it != rules.end()) {
            fmt::println("    {{\"{}\", \"{}\"}},", basis_element, it->second);
        }
    }
}

// Helper function to do complete rule comparison
bool compare_all_rules(const prd_rules& generated, const prd_rules& reference,
                       const std::string& product_name)
{
    fmt::println("\n=== Complete {} Validation ===", product_name);

    bool all_match = true;
    int total_rules = 0;
    int matching_rules = 0;

    // Check all reference rules exist and match in generated
    for (const auto& [key, ref_value] : reference) {
        total_rules++;
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
        else {
            matching_rules++;
        }
    }

    // Check for extra rules in generated that aren't in reference
    for (const auto& [key, gen_value] : generated) {
        if (reference.find(key) == reference.end()) {
            fmt::println("  EXTRA: {} -> {}", key, gen_value);
            all_match = false;
        }
    }

    fmt::println("  rules: {}/{} match", matching_rules, total_rules);
    fmt::println("  generated rules count: {}", generated.size());
    fmt::println("  reference rules count: {}", reference.size());

    if (all_match && generated.size() == reference.size()) {
        fmt::println("  ✓ PERFECT MATCH: All {} rules match character-by-character",
                     product_name);
    }
    else {
        fmt::println("  ✗ DIFFERENCES: Some {} rules differ from reference",
                     product_name);
    }

    return all_match && generated.size() == reference.size();
}

// Helper function to build antimetric matrix from complement rules
std::vector<int> build_antimetric_matrix(mvec_coeff const& basis,
                                         prd_rules const& complement_rules)
{
    size_t const n = basis.size();
    std::vector<int> antimetric(n * n, 0);
    std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> G_anti{
        antimetric.data(), n, n};

    // For each basis element, find its complement and encode in matrix
    for (size_t i = 0; i < n; ++i) {
        auto it = complement_rules.find(basis[i]);
        if (it != complement_rules.end()) {
            std::string cmpl_value = it->second;

            // Parse sign and basis element
            int sign = 1;
            if (!cmpl_value.empty() && cmpl_value[0] == '-') {
                sign = -1;
                cmpl_value = cmpl_value.substr(1);
            }

            // Find which basis element this is
            for (size_t j = 0; j < n; ++j) {
                if (basis[j] == cmpl_value) {
                    G_anti[i, j] = sign;
                    break;
                }
            }
        }
    }

    return antimetric;
}

// Matrix multiplication for extended metric verification
std::vector<int> multiply_matrices(std::vector<int> const& A, std::vector<int> const& B,
                                    size_t n)
{
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> A_view{
        A.data(), n, n};
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> B_view{
        B.data(), n, n};

    std::vector<int> result(n * n, 0);
    std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> C{
        result.data(), n, n};

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                C[i, j] += A_view[i, k] * B_view[k, j];
            }
        }
    }

    return result;
}

// Investigate antimetric relationship: G · 𝔾 = det(g) · I
void investigate_antimetric_relationship(AlgebraConfig const& config,
                                         std::string const& algebra_name)
{
    fmt::println("\n=== Antimetric Relationship Investigation: {} ===", algebra_name);

    // Generate rules
    auto rules = generate_algebra_rules(config);

    // Get extended metric matrix
    auto G_data = calculate_extended_metric_matrix(config);
    size_t n = rules.basis.size();

    // Calculate determinant of base metric
    int det_g = 1;
    for (int m : config.metric_signature) {
        det_g *= m;
    }

    fmt::println("Metric signature: {}", fmt::join(config.metric_signature, ", "));
    fmt::println("det(g) = {}", det_g);

    // Get complement rules (for odd/even dimensional algebras)
    prd_rules complement_rules;
    bool is_even = (config.basis_vectors.size() % 2 == 0);

    if (is_even) {
        // Use right complement for even algebras
        complement_rules = rules.right_complement;
        fmt::println("Algebra dimension: {} (even) - using right complement", config.basis_vectors.size());
    } else {
        complement_rules = rules.complement;
        fmt::println("Algebra dimension: {} (odd) - using complement", config.basis_vectors.size());
    }

    // Build antimetric matrix from complement rules
    auto G_anti_data = build_antimetric_matrix(rules.basis, complement_rules);

    // Compute G · 𝔾
    auto product = multiply_matrices(G_data, G_anti_data, n);

    // Check if product equals det(g) · I
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> product_view{
        product.data(), n, n};

    fmt::println("\nVerifying: G · 𝔾 = det(g) · I");
    fmt::println("Expected diagonal value: {}", det_g);

    bool matches = true;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            int expected = (i == j) ? det_g : 0;
            if (product_view[i, j] != expected) {
                if (matches) {
                    fmt::println("Mismatches found:");
                }
                fmt::println("  [{},{}]: got {}, expected {}", i, j, product_view[i, j], expected);
                matches = false;
            }
        }
    }

    if (matches) {
        fmt::println("✓ VERIFIED: G · 𝔾 = det(g) · I holds for {}", algebra_name);
    } else {
        fmt::println("✗ FAILED: G · 𝔾 ≠ det(g) · I for {}", algebra_name);
    }
}

// Function to display generated rules for a specific algebra
void display_algebra_rules(const AlgebraConfig& config, const std::string& algebra_name)
{
    std::string separator(80, '=');
    fmt::println("\n{}", separator);
    fmt::println("=== {} algebra rules ===", algebra_name);
    fmt::println("algebra: {} - {}", algebra_name,
                 algebra_name.find("ega2d") != std::string::npos    ? "Euclidean 2d"
                 : algebra_name.find("ega3d") != std::string::npos  ? "Euclidean 3d"
                 : algebra_name.find("pga2dp") != std::string::npos ? "projective 2d"
                 : algebra_name.find("pga3dp") != std::string::npos ? "projective 3d"
                 : algebra_name.find("sta4d") != std::string::npos  ? "space-time 3d"
                                                                    : "UNKNOWN");
    fmt::println("basis vectors: {}", fmt::join(config.basis_vectors, ", "));
    fmt::println("metric signature: {}", fmt::join(config.metric_signature, ", "));

    // Generate all rules
    auto generated_rules = generate_algebra_rules(config);

    fmt::println("\n=== generated basis (grade-ordered) ===");
    fmt::println("basis: {}", fmt::join(generated_rules.basis, ", "));

    if (generated_rules.basis.size() == 4) {
        fmt::println("  grade 0 (scalar):  {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors): {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2]);
        fmt::println("  grade 2 (bivector): {}", generated_rules.basis[3]);
    }
    else if (generated_rules.basis.size() == 8) {
        fmt::println("  grade 0 (scalar):   {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors):  {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3]);
        fmt::println("  grade 2 (bivectors): {}, {}, {}", generated_rules.basis[4],
                     generated_rules.basis[5], generated_rules.basis[6]);
        std::string trivector_name = generated_rules.basis[7];
        fmt::println("  grade 3 (trivector): {} {}", trivector_name,
                     trivector_name == "e123"   ? "(EGA)"
                     : trivector_name == "e321" ? "(PGA)"
                                                : "");
    }
    else if (generated_rules.basis.size() == 16) {
        fmt::println("  grade 0 (scalar):     {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors):    {}, {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3],
                     generated_rules.basis[4]);
        fmt::println("  grade 2 (bivectors):  {}, {}, {}, {}, {}, {}",
                     generated_rules.basis[5], generated_rules.basis[6],
                     generated_rules.basis[7], generated_rules.basis[8],
                     generated_rules.basis[9], generated_rules.basis[10]);
        fmt::println("  grade 3 (trivectors): {}, {}, {}, {}", generated_rules.basis[11],
                     generated_rules.basis[12], generated_rules.basis[13],
                     generated_rules.basis[14]);
        fmt::println("  grade 4 (pseudoscalar): {}", generated_rules.basis[15]);
    }

    // Calculate and display extended metric (diagonal only - for reference)
    auto extended_metric = calculate_extended_metric(config);
    fmt::println("\n{} extended metric (diagonal)", algebra_name);
    fmt::println("extended metric diagonal: {}", fmt::join(extended_metric, ", "));

    // Show mapping of basis elements to their extended metric values
    fmt::println("basis → extended metric diagonal mapping:");
    for (size_t i = 0; i < std::min(generated_rules.basis.size(), extended_metric.size());
         ++i) {
        fmt::println("  {} → {}", generated_rules.basis[i], extended_metric[i]);
    }

    // Detect PGA algebras (have degenerate dimensions with metric signature 0)
    bool is_pga = false;
    for (int metric_val : config.metric_signature) {
        if (metric_val == 0) {
            is_pga = true;
            break;
        }
    }

    // Calculate and display regressive extended metric for PGA algebras
    if (is_pga) {
        auto regressive_metric = calculate_regressive_extended_metric(config);
        fmt::println("\n{} regressive extended metric (diagonal)", algebra_name);
        fmt::println("regressive extended metric diagonal: {}",
                     fmt::join(regressive_metric, ", "));

        // Show mapping of basis elements to their regressive extended metric values
        fmt::println("basis → regressive extended metric diagonal mapping:");
        for (size_t i = 0;
             i < std::min(generated_rules.basis.size(), regressive_metric.size()); ++i) {
            fmt::println("  {} → {}", generated_rules.basis[i], regressive_metric[i]);
        }
    }

    // Calculate and display full extended metric matrix
    auto matrix_data = calculate_extended_metric_matrix(config);
    const size_t n = generated_rules.basis.size();
    std::mdspan G{matrix_data.data(), n, n};

    fmt::println("\n{} extended metric matrix (full)", algebra_name);
    fmt::println("Extended metric matrix G ({}×{}):", n, n);

    // Print column headers
    fmt::print("      ");
    for (size_t j = 0; j < n; ++j) {
        fmt::print("{:>6}", generated_rules.basis[j]);
    }
    fmt::println("");

    // Print matrix rows with row headers
    for (size_t i = 0; i < n; ++i) {
        fmt::print("{:>6}", generated_rules.basis[i]);
        for (size_t j = 0; j < n; ++j) {
            fmt::print("{:>6}", G[i, j]);
        }
        fmt::println("");
    }

    // Calculate and display full regressive extended metric matrix for PGA algebras
    if (is_pga) {
        auto regressive_matrix_data = calculate_regressive_extended_metric_matrix(config);
        std::mdspan G_reg{regressive_matrix_data.data(), n, n};

        fmt::println("\n{} regressive extended metric matrix (full)", algebra_name);
        fmt::println("Regressive extended metric matrix Ḡ ({}×{}):", n, n);

        // Print column headers
        fmt::print("      ");
        for (size_t j = 0; j < n; ++j) {
            fmt::print("{:>6}", generated_rules.basis[j]);
        }
        fmt::println("");

        // Print matrix rows with row headers
        for (size_t i = 0; i < n; ++i) {
            fmt::print("{:>6}", generated_rules.basis[i]);
            for (size_t j = 0; j < n; ++j) {
                fmt::print("{:>6}", G_reg[i, j]);
            }
            fmt::println("");
        }
    }

    // Print all generated rules clearly with algebra prefix
    print_all_rules(generated_rules.geometric_product, "geometric product rules",
                    generated_rules.basis, algebra_name);
    print_all_rules(generated_rules.wedge_product, "wedge product rules",
                    generated_rules.basis, algebra_name);
    print_all_rules(generated_rules.dot_product, "dot product rules",
                    generated_rules.basis, algebra_name);

    // Print complement rules based on algebra dimensionality in grade order
    size_t num_basis_vectors = config.basis_vectors.size();
    bool is_even_dimensional = (num_basis_vectors % 2 == 0);

    if (is_even_dimensional) {
        // Even-dimensional algebras have left and right complements
        print_complement_rules(generated_rules.left_complement, "left complement rules",
                               generated_rules.basis, algebra_name);
        print_complement_rules(generated_rules.right_complement, "right complement rules",
                               generated_rules.basis, algebra_name);

        // Print dual rules (only for non-PGA algebras)
        if (is_pga) {
            // For even-dimensional PGA, only print bulk and weight duals
            print_complement_rules(generated_rules.left_bulk_dual, "left bulk dual rules",
                                   generated_rules.basis, algebra_name);
            print_complement_rules(generated_rules.right_bulk_dual,
                                   "right bulk dual rules", generated_rules.basis,
                                   algebra_name);
            print_complement_rules(generated_rules.left_weight_dual,
                                   "left weight dual rules", generated_rules.basis,
                                   algebra_name);
            print_complement_rules(generated_rules.right_weight_dual,
                                   "right weight dual rules", generated_rules.basis,
                                   algebra_name);
        }
        else {
            // For non-PGA even-dimensional algebras (EGA, STA), print standard duals
            print_complement_rules(generated_rules.left_dual, "left dual rules",
                                   generated_rules.basis, algebra_name);
            print_complement_rules(generated_rules.right_dual, "right dual rules",
                                   generated_rules.basis, algebra_name);
        }
    }
    else {
        // Odd-dimensional algebras have a single complement
        print_complement_rules(generated_rules.complement, "complement rules",
                               generated_rules.basis, algebra_name);

        // Print dual rules (only for non-PGA algebras)
        if (is_pga) {
            // For odd-dimensional PGA, only print bulk and weight duals
            print_complement_rules(generated_rules.bulk_dual, "bulk dual rules",
                                   generated_rules.basis, algebra_name);
            print_complement_rules(generated_rules.weight_dual, "weight dual rules",
                                   generated_rules.basis, algebra_name);
        }
        else {
            // For non-PGA odd-dimensional algebras (EGA), print standard dual
            print_complement_rules(generated_rules.dual, "dual rules",
                                   generated_rules.basis, algebra_name);
        }
    }
}

// Function to test a specific algebra with complement validation
bool test_algebra_with_complements(
    const AlgebraConfig& config, const std::string& algebra_name,
    const mvec_coeff& reference_basis, const prd_rules& reference_gpr,
    const prd_rules& reference_wdg, const prd_rules& reference_dot,
    const prd_rules* reference_lcmpl = nullptr,
    const prd_rules* reference_rcmpl = nullptr, const prd_rules* reference_cmpl = nullptr)
{
    std::string separator(80, '=');
    fmt::println("\n{}", separator);
    fmt::println("=== testing {} ===", algebra_name);
    fmt::println("algebra: {} - {}", algebra_name,
                 algebra_name.find("ega2d") != std::string::npos    ? "Euclidean 2d"
                 : algebra_name.find("ega3d") != std::string::npos  ? "Euclidean 3d"
                 : algebra_name.find("pga2dp") != std::string::npos ? "projective 2d"
                 : algebra_name.find("pga3dp") != std::string::npos ? "projective 3d"
                 : algebra_name.find("sta4d") != std::string::npos  ? "space-time 3d"
                                                                    : "UNKNOWN");
    fmt::println("basis vectors: {}", fmt::join(config.basis_vectors, ", "));
    fmt::println("metric signature: {}", fmt::join(config.metric_signature, ", "));

    // Generate all rules
    auto generated_rules = generate_algebra_rules(config);

    fmt::println("\n=== generated basis (grade-ordered) ===");
    fmt::println("generated: {}", fmt::join(generated_rules.basis, ", "));

    if (generated_rules.basis.size() == 4) {
        fmt::println("  grade 0 (scalar):  {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors): {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2]);
        fmt::println("  grade 2 (bivector): {}", generated_rules.basis[3]);
    }
    else if (generated_rules.basis.size() == 8) {
        fmt::println("  grade 0 (scalar):   {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors):  {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3]);
        fmt::println("  grade 2 (bivectors): {}, {}, {}", generated_rules.basis[4],
                     generated_rules.basis[5], generated_rules.basis[6]);
        std::string trivector_name = generated_rules.basis[7];
        fmt::println("  grade 3 (trivector): {} {}", trivector_name,
                     trivector_name == "e123"   ? "(EGA)"
                     : trivector_name == "e321" ? "(PGA)"
                                                : "");
    }
    else if (generated_rules.basis.size() == 16) {
        fmt::println("  grade 0 (scalar):     {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors):    {}, {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3],
                     generated_rules.basis[4]);
        fmt::println("  grade 2 (bivectors):  {}, {}, {}, {}, {}, {}",
                     generated_rules.basis[5], generated_rules.basis[6],
                     generated_rules.basis[7], generated_rules.basis[8],
                     generated_rules.basis[9], generated_rules.basis[10]);
        fmt::println("  grade 3 (trivectors): {}, {}, {}, {}", generated_rules.basis[11],
                     generated_rules.basis[12], generated_rules.basis[13],
                     generated_rules.basis[14]);
        fmt::println("  grade 4 (pseudoscalar): {}", generated_rules.basis[15]);
    }

    fmt::println("\nreference: {}", fmt::join(reference_basis, ", "));

    // Calculate and display extended metric (diagonal only - for reference)
    auto extended_metric = calculate_extended_metric(config);
    fmt::println("\n{} extended metric (diagonal)", algebra_name);
    fmt::println("extended metric diagonal: {}", fmt::join(extended_metric, ", "));

    // show mapping of basis elements to their extended metric values
    fmt::println("basis → extended metric diagonal mapping:");
    for (size_t i = 0; i < std::min(generated_rules.basis.size(), extended_metric.size());
         ++i) {
        fmt::println("  {} → {}", generated_rules.basis[i], extended_metric[i]);
    }

    // Detect PGA algebras (have degenerate dimensions with metric signature 0)
    bool is_pga_test = false;
    for (int metric_val : config.metric_signature) {
        if (metric_val == 0) {
            is_pga_test = true;
            break;
        }
    }

    // Calculate and display regressive extended metric for PGA algebras
    if (is_pga_test) {
        auto regressive_metric = calculate_regressive_extended_metric(config);
        fmt::println("\n{} regressive extended metric (diagonal)", algebra_name);
        fmt::println("regressive extended metric diagonal: {}",
                     fmt::join(regressive_metric, ", "));

        // Show mapping of basis elements to their regressive extended metric values
        fmt::println("basis → regressive extended metric diagonal mapping:");
        for (size_t i = 0;
             i < std::min(generated_rules.basis.size(), regressive_metric.size()); ++i) {
            fmt::println("  {} → {}", generated_rules.basis[i], regressive_metric[i]);
        }
    }

    // Calculate and display full extended metric matrix
    auto matrix_data = calculate_extended_metric_matrix(config);
    const size_t n = generated_rules.basis.size();
    std::mdspan G{matrix_data.data(), n, n};

    fmt::println("\n{} extended metric matrix (full)", algebra_name);
    fmt::println("Extended metric matrix G ({}×{}):", n, n);

    // Print column headers
    fmt::print("      ");
    for (size_t j = 0; j < n; ++j) {
        fmt::print("{:>6}", generated_rules.basis[j]);
    }
    fmt::println("");

    // Print matrix rows with row headers
    for (size_t i = 0; i < n; ++i) {
        fmt::print("{:>6}", generated_rules.basis[i]);
        for (size_t j = 0; j < n; ++j) {
            fmt::print("{:>6}", G[i, j]);
        }
        fmt::println("");
    }

    // Calculate and display full regressive extended metric matrix for PGA algebras
    if (is_pga_test) {
        auto regressive_matrix_data = calculate_regressive_extended_metric_matrix(config);
        std::mdspan G_reg{regressive_matrix_data.data(), n, n};

        fmt::println("\n{} regressive extended metric matrix (full)", algebra_name);
        fmt::println("Regressive extended metric matrix Ḡ ({}×{}):", n, n);

        // Print column headers
        fmt::print("      ");
        for (size_t j = 0; j < n; ++j) {
            fmt::print("{:>6}", generated_rules.basis[j]);
        }
        fmt::println("");

        // Print matrix rows with row headers
        for (size_t i = 0; i < n; ++i) {
            fmt::print("{:>6}", generated_rules.basis[i]);
            for (size_t j = 0; j < n; ++j) {
                fmt::print("{:>6}", G_reg[i, j]);
            }
            fmt::println("");
        }
    }

    // validate basis
    bool basis_match = (generated_rules.basis == reference_basis);
    fmt::println("basis match: {}", basis_match ? "✓ PERFECT" : "✗ DIFFERENT");

    // Print complete generated rules
    print_all_rules(generated_rules.geometric_product, "geometric product rules",
                    generated_rules.basis);
    print_all_rules(reference_gpr, "reference geometric product rules", reference_basis);

    print_all_rules(generated_rules.wedge_product, "wedge product rules",
                    generated_rules.basis);
    print_all_rules(reference_wdg, "reference wedge product rules", reference_basis);

    print_all_rules(generated_rules.dot_product, "dot product rules",
                    generated_rules.basis);
    print_all_rules(reference_dot, "reference dot product rules", reference_basis);

    // Do complete character-by-character validation
    fmt::println("\n{}", separator);
    fmt::println("COMPREHENSIVE VALIDATION - CHARACTER-BY-CHARACTER COMPARISON");
    fmt::println("{}", separator);

    bool gpr_perfect = compare_all_rules(generated_rules.geometric_product, reference_gpr,
                                         "geometric product");
    bool wdg_perfect =
        compare_all_rules(generated_rules.wedge_product, reference_wdg, "wedge product");
    bool dot_perfect =
        compare_all_rules(generated_rules.dot_product, reference_dot, "dot product");

    // validate complements if provided
    bool lcmpl_perfect = true, rcmpl_perfect = true, cmpl_perfect = true;

    if (reference_lcmpl) {
        lcmpl_perfect = compare_all_rules(generated_rules.left_complement,
                                          *reference_lcmpl, "left complement");
    }
    if (reference_rcmpl) {
        rcmpl_perfect = compare_all_rules(generated_rules.right_complement,
                                          *reference_rcmpl, "right complement");
    }
    if (reference_cmpl) {
        cmpl_perfect =
            compare_all_rules(generated_rules.complement, *reference_cmpl, "complement");
    }

    // Final summary
    fmt::println("\n{}", separator);
    fmt::println("FINAL VALIDATION SUMMARY - {}", algebra_name);
    fmt::println("{}", separator);
    fmt::println("basis generation:     {}", basis_match ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("geometric product:    {}", gpr_perfect ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("wedge product:        {}", wdg_perfect ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("dot product:          {}", dot_perfect ? "✓ PERFECT" : "✗ FAILED");

    if (reference_lcmpl) {
        fmt::println("left complement:      {}",
                     lcmpl_perfect ? "✓ PERFECT" : "✗ FAILED");
    }
    if (reference_rcmpl) {
        fmt::println("right complement:     {}",
                     rcmpl_perfect ? "✓ PERFECT" : "✗ FAILED");
    }
    if (reference_cmpl) {
        fmt::println("complement:           {}", cmpl_perfect ? "✓ PERFECT" : "✗ FAILED");
    }

    bool overall_success = basis_match && gpr_perfect && wdg_perfect && dot_perfect &&
                           lcmpl_perfect && rcmpl_perfect && cmpl_perfect;
    fmt::println("\nOVERALL RESULT: {}",
                 overall_success ? "✓ 100% PERFECT MATCH" : "✗ DIFFERENCES FOUND");

    return overall_success;
}

// Function to test a specific algebra (keep for compatibility)
bool test_algebra(const AlgebraConfig& config, const std::string& algebra_name,
                  const mvec_coeff& reference_basis, const prd_rules& reference_gpr,
                  const prd_rules& reference_wdg, const prd_rules& reference_dot)
{

    std::string separator(80, '=');
    fmt::println("\n{}", separator);
    fmt::println("=== testing {} ===", algebra_name);
    fmt::println("algebra: {} - {}", algebra_name,
                 algebra_name.find("ega2d") != std::string::npos    ? "Euclidean 2d"
                 : algebra_name.find("ega3d") != std::string::npos  ? "Euclidean 3d"
                 : algebra_name.find("pga2dp") != std::string::npos ? "projective 2d"
                 : algebra_name.find("pga3dp") != std::string::npos ? "projective 3d"
                 : algebra_name.find("sta4d") != std::string::npos  ? "space-time 3d"
                                                                    : "UNKNOWN");
    fmt::println("basis vectors: {}", fmt::join(config.basis_vectors, ", "));
    fmt::println("metric signature: {}", fmt::join(config.metric_signature, ", "));

    // Generate all rules
    auto generated_rules = generate_algebra_rules(config);

    fmt::println("\n=== generated basis (grade-ordered) ===");
    fmt::println("generated: {}", fmt::join(generated_rules.basis, ", "));

    if (generated_rules.basis.size() == 4) {
        fmt::println("  grade 0 (scalar):  {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors): {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2]);
        fmt::println("  grade 2 (bivector): {}", generated_rules.basis[3]);
    }
    else if (generated_rules.basis.size() == 8) {
        fmt::println("  grade 0 (scalar):   {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors):  {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3]);
        fmt::println("  grade 2 (bivectors): {}, {}, {}", generated_rules.basis[4],
                     generated_rules.basis[5], generated_rules.basis[6]);
        std::string trivector_name = generated_rules.basis[7];
        fmt::println("  grade 3 (trivector): {} {}", trivector_name,
                     trivector_name == "e123"   ? "(EGA)"
                     : trivector_name == "e321" ? "(PGA)"
                                                : "");
    }
    else if (generated_rules.basis.size() == 16) {
        fmt::println("  grade 0 (scalar):     {}", generated_rules.basis[0]);
        fmt::println("  grade 1 (vectors):    {}, {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3],
                     generated_rules.basis[4]);
        fmt::println("  grade 2 (bivectors):  {}, {}, {}, {}, {}, {}",
                     generated_rules.basis[5], generated_rules.basis[6],
                     generated_rules.basis[7], generated_rules.basis[8],
                     generated_rules.basis[9], generated_rules.basis[10]);
        fmt::println("  grade 3 (trivectors): {}, {}, {}, {}", generated_rules.basis[11],
                     generated_rules.basis[12], generated_rules.basis[13],
                     generated_rules.basis[14]);
        fmt::println("  grade 4 (pseudoscalar): {}", generated_rules.basis[15]);
    }

    fmt::println("\nreference: {}", fmt::join(reference_basis, ", "));

    // calculate and display extended metric
    auto extended_metric = calculate_extended_metric(config);
    fmt::println("\n=== extended metric ===");
    fmt::println("extended metric: {}", fmt::join(extended_metric, ", "));

    // show mapping of basis elements to their extended metric values
    fmt::println("basis → extended metric mapping:");
    for (size_t i = 0; i < std::min(generated_rules.basis.size(), extended_metric.size());
         ++i) {
        fmt::println("  {} → {}", generated_rules.basis[i], extended_metric[i]);
    }

    // Validate basis
    bool basis_match = (generated_rules.basis == reference_basis);
    fmt::println("basis match: {}", basis_match ? "✓ PERFECT" : "✗ DIFFERENT");

    // print complete generated rules
    print_all_rules(generated_rules.geometric_product, "geometric product rules",
                    generated_rules.basis);
    print_all_rules(reference_gpr, "reference geometric product rules", reference_basis);

    print_all_rules(generated_rules.wedge_product, "wedge product rules",
                    generated_rules.basis);
    print_all_rules(reference_wdg, "reference wedge product rules", reference_basis);

    print_all_rules(generated_rules.dot_product, "dot product rules",
                    generated_rules.basis);
    print_all_rules(reference_dot, "reference dot product rules", reference_basis);

    // Do complete character-by-character validation
    fmt::println("\n{}", separator);
    fmt::println("COMPREHENSIVE VALIDATION - CHARACTER-BY-CHARACTER COMPARISON");
    fmt::println("{}", separator);

    bool gpr_perfect = compare_all_rules(generated_rules.geometric_product, reference_gpr,
                                         "geometric product");
    bool wdg_perfect =
        compare_all_rules(generated_rules.wedge_product, reference_wdg, "wedge product");
    bool dot_perfect =
        compare_all_rules(generated_rules.dot_product, reference_dot, "dot product");

    // Final summary
    fmt::println("\n{}", separator);
    fmt::println("FINAL VALIDATION SUMMARY - {}", algebra_name);
    fmt::println("{}", separator);
    fmt::println("basis generation:     {}", basis_match ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("geometric product:    {}", gpr_perfect ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("wedge product:        {}", wdg_perfect ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("dot product:          {}", dot_perfect ? "✓ PERFECT" : "✗ FAILED");

    bool overall_success = basis_match && gpr_perfect && wdg_perfect && dot_perfect;
    fmt::println("\nOVERALL RESULT: {}",
                 overall_success ? "✓ 100% PERFECT MATCH" : "✗ DIFFERENCES FOUND");

    return overall_success;
}

// Test recursive extraction of extended metric
// Compares recursive extraction against the current hardcoded implementation
void test_extended_metric_recursive_extraction(AlgebraConfig const& config,
                                                std::string const& algebra_name)
{
    std::string separator(80, '-');
    fmt::println("\n{}", separator);
    fmt::println("RECURSIVE EXTENDED METRIC EXTRACTION TEST - {}", algebra_name);
    fmt::println("{}", separator);

    // Generate wedge and dot product rules (independent of extended metric)
    auto wedge_rules = generate_wedge_product_rules(config);
    auto dot_rules = generate_dot_product_rules(config);

    // Get extended metric from current hardcoded implementation
    auto metric_hardcoded = calculate_extended_metric(config);

    // Get extended metric from new recursive extraction
    auto metric_recursive = calculate_extended_metric_recursive(config, wedge_rules, dot_rules);

    // Compare the results
    bool perfect_match = (metric_hardcoded == metric_recursive);

    fmt::println("\nComparison Results:");
    fmt::println("  Basis size: {}", config.multivector_basis.size());
    fmt::println("  Match status: {}", perfect_match ? "✓ PERFECT" : "✗ DIFFERENCES FOUND");

    if (!perfect_match) {
        fmt::println("\nDifferences found:");
        for (size_t i = 0; i < config.multivector_basis.size(); ++i) {
            if (metric_hardcoded[i] != metric_recursive[i]) {
                fmt::println("  [{}] {}: hardcoded={}, recursive={}",
                             i, config.multivector_basis[i],
                             metric_hardcoded[i], metric_recursive[i]);
            }
        }
    }
    else {
        fmt::println("\n✓ Recursive extraction produces identical results!");
        fmt::println("  All {} extended metric values match.", metric_hardcoded.size());
    }
}

int main(int argc, char* argv[])
{
    try {
        // Check if user wants consistency test mode (display is now default)
        //
        // NOTE: The --test_consistency mode validates DETERMINISM and CONSISTENCY
        // of rule generation by comparing freshly generated rules against rules
        // generated at static initialization time. It does NOT compare against
        // hand-crafted mathematical references.
        //
        // What it validates:
        // - Determinism: Generator produces identical results with same configuration
        // - Static initialization stability: No initialization order issues
        // - Configuration consistency: Test config matches production config
        // - Regression detection: Algorithm changes are caught immediately
        //
        // This serves as a safety net during refactoring and ensures the rule
        // generation system remains stable across code changes.
        bool test_consistency =
            (argc > 1 && std::string(argv[1]) == "--test_consistency");

        if (test_consistency) {
            fmt::println("=== Automatic GA Rule Generation - Multi-Algebra Test ===");
            fmt::println("Running determinism and consistency validation...");
        }
        else {
            fmt::println("=== Automatic GA Rule Generation - Rule Display ===");
            fmt::println(
                "Tip: Use '--test_consistency' argument to run determinism testing");
        }

        // Test results (only used in test mode)
        std::vector<bool> test_results;

        // Configure and test EGA2D algebra
        // Extract basis prefix from vector basis and validate consistency
        std::string const ega2d_prefix = extract_basis_prefix(mv2d_basis_kvec[1]);
        validate_basis_consistency(mv2d_basis, mv2d_basis_kvec, ega2d_prefix, one_str());

        AlgebraConfig ega2d_config = {
            .basis_vectors = mv2d_basis_kvec[1],       // Use vector basis from header
            .metric_signature = mv2d_metric_signature, // Use metric from header
            .multivector_basis = mv2d_basis,           // Use mv2d_basis from header
            .scalar_name = one_str(),
            .basis_prefix = ega2d_prefix}; // Use extracted and validated prefix

        if (test_consistency) {
            bool ega2d_success = test_algebra_with_complements(
                ega2d_config, "ega2d", mv2d_basis, gpr_ega2d_rules, wdg_ega2d_rules,
                dot_ega2d_rules, &lcmpl_ega2d_rules, &rcmpl_ega2d_rules);
            test_results.push_back(ega2d_success);
        }
        else {
            display_algebra_rules(ega2d_config, "ega2d");
        }

        // Configure and test EGA3D algebra
        // Extract basis prefix from vector basis and validate consistency
        std::string const ega3d_prefix = extract_basis_prefix(mv3d_basis_kvec[1]);
        validate_basis_consistency(mv3d_basis, mv3d_basis_kvec, ega3d_prefix, one_str());

        AlgebraConfig ega3d_config = {
            .basis_vectors = mv3d_basis_kvec[1],       // Use vector basis from header
            .metric_signature = mv3d_metric_signature, // Use metric from header
            .multivector_basis = mv3d_basis,           // Use mv3d_basis from header
            .scalar_name = one_str(),
            .basis_prefix = ega3d_prefix}; // Use extracted and validated prefix

        if (test_consistency) {
            bool ega3d_success = test_algebra_with_complements(
                ega3d_config, "ega3d", mv3d_basis, gpr_ega3d_rules, wdg_ega3d_rules,
                dot_ega3d_rules, nullptr, nullptr, &cmpl_ega3d_rules);
            test_results.push_back(ega3d_success);
        }
        else {
            display_algebra_rules(ega3d_config, "ega3d");
        }

        // Configure and test PGA2DP algebra
        // Extract basis prefix from vector basis and validate consistency
        std::string const pga2dp_prefix = extract_basis_prefix(mv2dp_basis_kvec[1]);
        validate_basis_consistency(mv2dp_basis, mv2dp_basis_kvec, pga2dp_prefix,
                                   one_str());

        AlgebraConfig pga2dp_config = {
            .basis_vectors = mv2dp_basis_kvec[1],       // Use vector basis from header
            .metric_signature = mv2dp_metric_signature, // Use metric from header
            .multivector_basis = mv2dp_basis,           // Use mv2dp_basis from header
            .scalar_name = one_str(),
            .basis_prefix = pga2dp_prefix}; // Use extracted and validated prefix

        if (test_consistency) {
            bool pga2dp_success = test_algebra_with_complements(
                pga2dp_config, "pga2dp", mv2dp_basis, gpr_pga2dp_rules, wdg_pga2dp_rules,
                dot_pga2dp_rules, nullptr, nullptr, &cmpl_pga2dp_rules);
            test_results.push_back(pga2dp_success);
        }
        else {
            display_algebra_rules(pga2dp_config, "pga2dp");
        }

        // Configure and test PGA3DP algebra
        // Extract basis prefix from vector basis and validate consistency
        std::string const pga3dp_prefix = extract_basis_prefix(mv3dp_basis_kvec[1]);
        validate_basis_consistency(mv3dp_basis, mv3dp_basis_kvec, pga3dp_prefix,
                                   one_str());

        AlgebraConfig pga3dp_config = {
            .basis_vectors = mv3dp_basis_kvec[1],       // Use vector basis from header
            .metric_signature = mv3dp_metric_signature, // Use metric from header
            .multivector_basis = mv3dp_basis,           // Use mv3dp_basis from header
            .scalar_name = one_str(),
            .basis_prefix = pga3dp_prefix}; // Use extracted and validated prefix

        if (test_consistency) {
            bool pga3dp_success = test_algebra_with_complements(
                pga3dp_config, "pga3dp", mv3dp_basis, gpr_pga3dp_rules, wdg_pga3dp_rules,
                dot_pga3dp_rules, &lcmpl_pga3dp_rules, &rcmpl_pga3dp_rules);
            test_results.push_back(pga3dp_success);
        }
        else {
            display_algebra_rules(pga3dp_config, "pga3dp");
        }

        // Configure and test STA4D algebra (Space-Time Algebra)
        // Extract basis prefix from vector basis and validate consistency
        std::string const sta4d_prefix = extract_basis_prefix(mvsta4d_basis_kvec[1]);
        validate_basis_consistency(mvsta4d_basis, mvsta4d_basis_kvec, sta4d_prefix,
                                   one_str());

        AlgebraConfig sta4d_config = {
            .basis_vectors = mvsta4d_basis_kvec[1],       // Use vector basis from header
            .metric_signature = mvsta4d_metric_signature, // Use metric from header
            .multivector_basis = mvsta4d_basis,           // Use mvsta4d_basis from header
            .scalar_name = one_str(),
            .basis_prefix = sta4d_prefix}; // Use extracted and validated prefix

        if (test_consistency) {
            bool sta4d_success = test_algebra(
                sta4d_config, "sta4d", mvsta4d_basis, gpr_sta4d_rules, wdg_sta4d_rules,
                dot_sta4d_rules); // testing all three product types - complete manual
                                  // implementation
            test_results.push_back(sta4d_success);
        }
        else {
            display_algebra_rules(sta4d_config, "sta4d");
        }

        // Test recursive extended metric extraction for all algebras
        fmt::println("\n");
        std::string separator_long(80, '=');
        fmt::println("{}", separator_long);
        fmt::println("RECURSIVE EXTENDED METRIC EXTRACTION TESTS");
        fmt::println("{}", separator_long);

        test_extended_metric_recursive_extraction(ega2d_config, "ega2d");
        test_extended_metric_recursive_extraction(ega3d_config, "ega3d");
        test_extended_metric_recursive_extraction(pga2dp_config, "pga2dp");
        test_extended_metric_recursive_extraction(pga3dp_config, "pga3dp");
        test_extended_metric_recursive_extraction(sta4d_config, "sta4d");

        // Investigate antimetric relationship for EGA3D and STA4D
        fmt::println("\n");
        fmt::println("{}", separator_long);
        fmt::println("ANTIMETRIC RELATIONSHIP INVESTIGATION");
        fmt::println("{}", separator_long);

        investigate_antimetric_relationship(ega3d_config, "ega3d");
        investigate_antimetric_relationship(sta4d_config, "sta4d");

        if (test_consistency) {
            // Overall summary for test mode
            std::string separator(80, '=');
            fmt::println("\n{}", separator);
            fmt::println("MULTI-ALGEBRA TEST SUMMARY");
            fmt::println("{}", separator);
            fmt::println("ega2d (G(2,0,0)): {}",
                         test_results[0] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("ega3d (G(3,0,0)): {}",
                         test_results[1] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("pga2dp (G(2,0,1)): {}",
                         test_results[2] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("pga3dp (G(3,0,1)): {}",
                         test_results[3] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("sta4d (G(1,3,0)): {}",
                         test_results[4] ? "✓ PERFECT" : "✗ FAILED");

            bool all_success = std::all_of(test_results.begin(), test_results.end(),
                                           [](bool b) { return b; });
            fmt::println("\nOVERALL RESULT: {}", all_success ? "✓ ALL ALGEBRAS PERFECT"
                                                             : "✗ SOME ALGEBRAS FAILED");

            if (all_success) {
                fmt::println("\n🎉 SUCCESS: Generated rules are character-identical to "
                             "reference implementations!");
                fmt::println(
                    "The automatic rule generation system is working perfectly for "
                    "all tested algebras.");
            }
            else {
                fmt::println("\n❌ Issues found in rule generation. Please check the "
                             "differences above.");
                return 1;
            }
        }
        else {
            // Simple completion message for display mode
            std::string separator(80, '=');
            fmt::println("\n{}", separator);
            fmt::println("RULE DISPLAY COMPLETE");
            fmt::println("{}", separator);
            fmt::println(
                "All generated rules have been displayed for visual inspection.");
            fmt::println(
                "The rules are generated algorithmically and available for use.");
        }

        return 0;
    }
    catch (const std::exception& e) {
        fmt::println("Error: {}", e.what());
        return 1;
    }
}