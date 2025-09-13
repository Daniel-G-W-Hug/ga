// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_ega2d.hpp"
#include "ga_prdxpr_ega3d.hpp"
#include "ga_prdxpr_pga2dp.hpp"
#include "ga_prdxpr_pga3dp.hpp"
#include "ga_prdxpr_rule_generator.hpp"
#include "ga_prdxpr_sta3d.hpp"

#include <iostream>

// Helper function to print all rules in grade-ordered way
void print_all_rules(const prd_rules& rules, const std::string& title,
                     const mvec_coeff& basis_order, const std::string& algebra_name = "")
{
    if (!algebra_name.empty()) {
        fmt::println("\n{} {}", algebra_name, title);
    } else {
        fmt::println("\n=== {} ===", title);
    }

    // Print in grade order by iterating through basis in order
    for (const auto& a : basis_order) {
        for (const auto& b : basis_order) {
            std::string key = a + " * " + b;
            if (title.find("wedge") != std::string::npos || title.find("Wedge") != std::string::npos) {
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
                            const mvec_coeff& basis_order, const std::string& algebra_name = "")
{
    if (!algebra_name.empty()) {
        fmt::println("\n{} {}", algebra_name, title);
    } else {
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

// Function to display generated rules for a specific algebra
void display_algebra_rules(const AlgebraConfig& config, const std::string& algebra_name)
{
    std::string separator(80, '=');
    fmt::println("\n{}", separator);
    fmt::println("=== {} algebra rules ===", algebra_name);
    fmt::println("algebra: {} - {}", algebra_name,
                 algebra_name.find("ega2d") != std::string::npos    ? "euclidean 2d"
                 : algebra_name.find("ega3d") != std::string::npos  ? "euclidean 3d"
                 : algebra_name.find("pga2dp") != std::string::npos ? "projective 2d"
                 : algebra_name.find("pga3dp") != std::string::npos ? "projective 3d"
                 : algebra_name.find("sta3d") != std::string::npos  ? "space-time 3d"
                                                                    : "Unknown");
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

    // Calculate and display extended metric
    auto extended_metric = calculate_extended_metric(config);
    fmt::println("\n=== extended metric ===");
    fmt::println("extended metric: {}", fmt::join(extended_metric, ", "));

    // Show mapping of basis elements to their extended metric values
    fmt::println("basis → extended metric mapping:");
    for (size_t i = 0; i < std::min(generated_rules.basis.size(), extended_metric.size());
         ++i) {
        fmt::println("  {} → {}", generated_rules.basis[i], extended_metric[i]);
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
    }
    else {
        // Odd-dimensional algebras have a single complement
        print_complement_rules(generated_rules.complement, "complement rules",
                               generated_rules.basis, algebra_name);
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
                 algebra_name.find("ega2d") != std::string::npos    ? "euclidean 2d"
                 : algebra_name.find("ega3d") != std::string::npos  ? "euclidean 3d"
                 : algebra_name.find("pga2dp") != std::string::npos ? "projective 2d"
                 : algebra_name.find("pga3dp") != std::string::npos ? "projective 3d"
                 : algebra_name.find("sta3d") != std::string::npos  ? "space-time 3d"
                                                                    : "Unknown");
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

    // Calculate and display extended metric
    auto extended_metric = calculate_extended_metric(config);
    fmt::println("\n=== extended metric ===");
    fmt::println("extended metric: {}", fmt::join(extended_metric, ", "));

    // show mapping of basis elements to their extended metric values
    fmt::println("basis → extended metric mapping:");
    for (size_t i = 0; i < std::min(generated_rules.basis.size(), extended_metric.size());
         ++i) {
        fmt::println("  {} → {}", generated_rules.basis[i], extended_metric[i]);
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
                 algebra_name.find("ega2d") != std::string::npos    ? "euclidean 2d"
                 : algebra_name.find("ega3d") != std::string::npos  ? "euclidean 3d"
                 : algebra_name.find("pga2dp") != std::string::npos ? "projective 2d"
                 : algebra_name.find("pga3dp") != std::string::npos ? "projective 3d"
                 : algebra_name.find("sta3d") != std::string::npos  ? "space-time 3d"
                                                                    : "Unknown");
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

int main(int argc, char* argv[])
{
    try {
        // Check if user wants consistency test mode (display is now default)
        bool test_consistency =
            (argc > 1 && std::string(argv[1]) == "--test_consistency");

        if (test_consistency) {
            fmt::println("=== Automatic GA Rule Generation - Multi-Algebra Test ===");
            fmt::println("Running internal consistency tests...");
        }
        else {
            fmt::println("=== Automatic GA Rule Generation - Rule Display ===");
            fmt::println(
                "Tip: Use '--test_consistency' argument to run comparison testing");
        }

        // Test results (only used in test mode)
        std::vector<bool> test_results;

        // Configure and test EGA2D algebra
        AlgebraConfig ega2d_config = {
            .basis_vectors = {"e1", "e2"},
            .metric_signature = {+1, +1},
            .multivector_basis = {"1", "e1", "e2",
                                  "e12"}, // Copy from ga_prdxpr_ega2d.hpp line 12
            .scalar_name = one_str(),
            .basis_prefix = "e"};

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
        AlgebraConfig ega3d_config = {
            .basis_vectors = {"e1", "e2", "e3"},
            .metric_signature = {+1, +1, +1},
            .multivector_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12",
                                  "e123"}, // Copy from ga_prdxpr_ega3d.hpp line 12
            .scalar_name = one_str(),
            .basis_prefix = "e"};

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
        AlgebraConfig pga2dp_config = {
            .basis_vectors = {"e1", "e2", "e3"},
            .metric_signature = {+1, +1, 0}, // G(2,0,1) - e1²=+1, e2²=+1, e3²=0
            .multivector_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12",
                                  "e321"}, // Copy from ga_prdxpr_pga2dp.hpp line 12
            .scalar_name = one_str(),
            .basis_prefix = "e"};

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
        AlgebraConfig pga3dp_config = {
            .basis_vectors = {"e1", "e2", "e3", "e4"},
            .metric_signature = {+1, +1, +1,
                                 0}, // G(3,0,1) - e1²=+1, e2²=+1, e3²=+1, e4²=0
            .multivector_basis = {"1", "e1", "e2", "e3", "e4", "e41", "e42", "e43", "e23",
                                  "e31", "e12", "e423", "e431", "e412", "e321",
                                  "e1234"}, // Copy from ga_prdxpr_pga3dp.hpp line 12-14
            .scalar_name = one_str(),
            .basis_prefix = "e"};

        if (test_consistency) {
            bool pga3dp_success = test_algebra_with_complements(
                pga3dp_config, "pga3dp", mv3dp_basis, gpr_pga3dp_rules, wdg_pga3dp_rules,
                dot_pga3dp_rules, &lcmpl_pga3dp_rules, &rcmpl_pga3dp_rules);
            test_results.push_back(pga3dp_success);
        }
        else {
            display_algebra_rules(pga3dp_config, "pga3dp");
        }

        // Configure and test STA3D algebra (Space-Time Algebra)
        AlgebraConfig sta3d_config = {
            .basis_vectors = {"g0", "g1", "g2", "g3"},
            .metric_signature = {+1, -1, -1,
                                 -1}, // G(1,3,0) - g0²=+1, g1²=-1, g2²=-1, g3²=-1
            .multivector_basis = {"1", "g0", "g1", "g2", "g3", "g01", "g02", "g03", "g23",
                                  "g31", "g12", "g023", "g031", "g012", "g123",
                                  "g0123"}, // Copy from ga_prdxpr_sta3d.hpp line 16-18
            .scalar_name = one_str(),
            .basis_prefix = "g" // Uses gamma notation instead of e
        };

        if (test_consistency) {
            bool sta3d_success = test_algebra(
                sta3d_config, "sta3d", mvsta3d_basis, gpr_sta3d_rules, wdg_sta3d_rules,
                dot_sta3d_rules); // testing all three product types - complete manual
                                  // implementation
            test_results.push_back(sta3d_success);
        }
        else {
            display_algebra_rules(sta3d_config, "sta3d");
        }

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
            fmt::println("sta3d (G(1,3,0)): {}",
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