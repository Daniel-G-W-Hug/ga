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
                     const mvec_coeff& basis_order)
{
    fmt::println("\n=== {} ===", title);

    // Print in grade order by iterating through basis in order
    for (const auto& a : basis_order) {
        for (const auto& b : basis_order) {
            std::string key = a + " * " + b;
            if (title.find("Wedge") != std::string::npos) {
                key = a + " ^ " + b; // Use wedge operator for wedge products
            }

            auto it = rules.find(key);
            if (it != rules.end()) {
                fmt::println("    {{\"{}\", \"{}\"}},", key, it->second);
            }
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

    fmt::println("  Rules: {}/{} match", matching_rules, total_rules);
    fmt::println("  Generated rules count: {}", generated.size());
    fmt::println("  Reference rules count: {}", reference.size());

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
    fmt::println("=== {} ALGEBRA RULES ===", algebra_name);
    fmt::println("Algebra: {} - {}", algebra_name,
                 algebra_name.find("EGA2D") != std::string::npos    ? "Euclidean 2D"
                 : algebra_name.find("EGA3D") != std::string::npos  ? "Euclidean 3D"
                 : algebra_name.find("PGA2DP") != std::string::npos ? "Projective 2D"
                 : algebra_name.find("PGA3DP") != std::string::npos ? "Projective 3D"
                 : algebra_name.find("STA3D") != std::string::npos  ? "Space-Time 3D"
                                                                   : "Unknown");
    fmt::println("Basis vectors: {}", fmt::join(config.basis_vectors, ", "));
    fmt::println("Metric signature: {}", fmt::join(config.metric_signature, ", "));

    // Generate all rules
    auto generated_rules = generate_algebra_rules(config);

    fmt::println("\n=== Generated Basis (Grade-Ordered) ===");
    fmt::println("Basis: {}", fmt::join(generated_rules.basis, ", "));

    if (generated_rules.basis.size() == 4) {
        fmt::println("  Grade 0 (scalar):  {}", generated_rules.basis[0]);
        fmt::println("  Grade 1 (vectors): {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2]);
        fmt::println("  Grade 2 (bivector): {}", generated_rules.basis[3]);
    }
    else if (generated_rules.basis.size() == 8) {
        fmt::println("  Grade 0 (scalar):   {}", generated_rules.basis[0]);
        fmt::println("  Grade 1 (vectors):  {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3]);
        fmt::println("  Grade 2 (bivectors): {}, {}, {}", generated_rules.basis[4],
                     generated_rules.basis[5], generated_rules.basis[6]);
        std::string trivector_name = generated_rules.basis[7];
        fmt::println("  Grade 3 (trivector): {} {}", trivector_name,
                     trivector_name == "e123"   ? "(EGA)"
                     : trivector_name == "e321" ? "(PGA)"
                                                : "");
    }
    else if (generated_rules.basis.size() == 16) {
        fmt::println("  Grade 0 (scalar):     {}", generated_rules.basis[0]);
        fmt::println("  Grade 1 (vectors):    {}, {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3],
                     generated_rules.basis[4]);
        fmt::println("  Grade 2 (bivectors):  {}, {}, {}, {}, {}, {}",
                     generated_rules.basis[5], generated_rules.basis[6],
                     generated_rules.basis[7], generated_rules.basis[8],
                     generated_rules.basis[9], generated_rules.basis[10]);
        fmt::println("  Grade 3 (trivectors): {}, {}, {}, {}", generated_rules.basis[11],
                     generated_rules.basis[12], generated_rules.basis[13],
                     generated_rules.basis[14]);
        fmt::println("  Grade 4 (pseudoscalar): {}", generated_rules.basis[15]);
    }

    // Calculate and display extended metric
    auto extended_metric = calculate_extended_metric(config);
    fmt::println("\n=== Extended Metric ===");
    fmt::println("Extended Metric: {}", fmt::join(extended_metric, ", "));

    // Show mapping of basis elements to their extended metric values
    fmt::println("Basis → Extended Metric mapping:");
    for (size_t i = 0; i < std::min(generated_rules.basis.size(), extended_metric.size());
         ++i) {
        fmt::println("  {} → {}", generated_rules.basis[i], extended_metric[i]);
    }

    // Print all generated rules clearly
    print_all_rules(generated_rules.geometric_product,
                    "Geometric Product Rules", generated_rules.basis);
    print_all_rules(generated_rules.wedge_product, "Wedge Product Rules",
                    generated_rules.basis);
    print_all_rules(generated_rules.dot_product, "Dot Product Rules",
                    generated_rules.basis);
}

// Function to test a specific algebra (keep for compatibility)
bool test_algebra(const AlgebraConfig& config, const std::string& algebra_name,
                  const mvec_coeff& reference_basis, const prd_rules& reference_gpr,
                  const prd_rules& reference_wdg, const prd_rules& reference_dot)
{

    std::string separator(80, '=');
    fmt::println("\n{}", separator);
    fmt::println("=== Testing {} ===", algebra_name);
    fmt::println("Algebra: {} - {}", algebra_name,
                 algebra_name.find("EGA2D") != std::string::npos    ? "Euclidean 2D"
                 : algebra_name.find("EGA3D") != std::string::npos  ? "Euclidean 3D"
                 : algebra_name.find("PGA2DP") != std::string::npos ? "Projective 2D"
                 : algebra_name.find("PGA3DP") != std::string::npos ? "Projective 3D"
                 : algebra_name.find("STA3D") != std::string::npos  ? "Space-Time 3D"
                                                                    : "Unknown");
    fmt::println("Basis vectors: {}", fmt::join(config.basis_vectors, ", "));
    fmt::println("Metric signature: {}", fmt::join(config.metric_signature, ", "));

    // Generate all rules
    auto generated_rules = generate_algebra_rules(config);

    fmt::println("\n=== Generated Basis (Grade-Ordered) ===");
    fmt::println("Generated: {}", fmt::join(generated_rules.basis, ", "));

    if (generated_rules.basis.size() == 4) {
        fmt::println("  Grade 0 (scalar):  {}", generated_rules.basis[0]);
        fmt::println("  Grade 1 (vectors): {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2]);
        fmt::println("  Grade 2 (bivector): {}", generated_rules.basis[3]);
    }
    else if (generated_rules.basis.size() == 8) {
        fmt::println("  Grade 0 (scalar):   {}", generated_rules.basis[0]);
        fmt::println("  Grade 1 (vectors):  {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3]);
        fmt::println("  Grade 2 (bivectors): {}, {}, {}", generated_rules.basis[4],
                     generated_rules.basis[5], generated_rules.basis[6]);
        std::string trivector_name = generated_rules.basis[7];
        fmt::println("  Grade 3 (trivector): {} {}", trivector_name,
                     trivector_name == "e123"   ? "(EGA)"
                     : trivector_name == "e321" ? "(PGA)"
                                                : "");
    }
    else if (generated_rules.basis.size() == 16) {
        fmt::println("  Grade 0 (scalar):     {}", generated_rules.basis[0]);
        fmt::println("  Grade 1 (vectors):    {}, {}, {}, {}", generated_rules.basis[1],
                     generated_rules.basis[2], generated_rules.basis[3],
                     generated_rules.basis[4]);
        fmt::println("  Grade 2 (bivectors):  {}, {}, {}, {}, {}, {}",
                     generated_rules.basis[5], generated_rules.basis[6],
                     generated_rules.basis[7], generated_rules.basis[8],
                     generated_rules.basis[9], generated_rules.basis[10]);
        fmt::println("  Grade 3 (trivectors): {}, {}, {}, {}", generated_rules.basis[11],
                     generated_rules.basis[12], generated_rules.basis[13],
                     generated_rules.basis[14]);
        fmt::println("  Grade 4 (pseudoscalar): {}", generated_rules.basis[15]);
    }

    fmt::println("\nReference: {}", fmt::join(reference_basis, ", "));

    // Calculate and display extended metric
    auto extended_metric = calculate_extended_metric(config);
    fmt::println("\n=== Extended Metric ===");
    fmt::println("Extended Metric: {}", fmt::join(extended_metric, ", "));

    // Show mapping of basis elements to their extended metric values
    fmt::println("Basis → Extended Metric mapping:");
    for (size_t i = 0; i < std::min(generated_rules.basis.size(), extended_metric.size());
         ++i) {
        fmt::println("  {} → {}", generated_rules.basis[i], extended_metric[i]);
    }

    // Validate basis
    bool basis_match = (generated_rules.basis == reference_basis);
    fmt::println("Basis match: {}", basis_match ? "✓ PERFECT" : "✗ DIFFERENT");

    // Print complete generated rules
    print_all_rules(generated_rules.geometric_product,
                    "Generated Geometric Product Rules", generated_rules.basis);
    print_all_rules(reference_gpr, "Reference Geometric Product Rules", reference_basis);

    print_all_rules(generated_rules.wedge_product, "Generated Wedge Product Rules",
                    generated_rules.basis);
    print_all_rules(reference_wdg, "Reference Wedge Product Rules", reference_basis);

    print_all_rules(generated_rules.dot_product, "Generated Dot Product Rules",
                    generated_rules.basis);
    print_all_rules(reference_dot, "Reference Dot Product Rules", reference_basis);

    // Do complete character-by-character validation
    fmt::println("\n{}", separator);
    fmt::println("COMPREHENSIVE VALIDATION - CHARACTER-BY-CHARACTER COMPARISON");
    fmt::println("{}", separator);

    bool gpr_perfect = compare_all_rules(generated_rules.geometric_product, reference_gpr,
                                         "Geometric Product");
    bool wdg_perfect =
        compare_all_rules(generated_rules.wedge_product, reference_wdg, "Wedge Product");
    bool dot_perfect =
        compare_all_rules(generated_rules.dot_product, reference_dot, "Dot Product");

    // Final summary
    fmt::println("\n{}", separator);
    fmt::println("FINAL VALIDATION SUMMARY - {}", algebra_name);
    fmt::println("{}", separator);
    fmt::println("Basis Generation:     {}", basis_match ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("Geometric Product:    {}", gpr_perfect ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("Wedge Product:        {}", wdg_perfect ? "✓ PERFECT" : "✗ FAILED");
    fmt::println("Dot Product:          {}", dot_perfect ? "✓ PERFECT" : "✗ FAILED");

    bool overall_success = basis_match && gpr_perfect && wdg_perfect && dot_perfect;
    fmt::println("\nOVERALL RESULT: {}",
                 overall_success ? "✓ 100% PERFECT MATCH" : "✗ DIFFERENCES FOUND");

    return overall_success;
}

int main(int argc, char* argv[])
{
    try {
        // Check if user wants consistency test mode (display is now default)
        bool test_consistency = (argc > 1 && std::string(argv[1]) == "--test_consistency");
        
        if (test_consistency) {
            fmt::println("=== Automatic GA Rule Generation - Multi-Algebra Test ===");
            fmt::println("Running internal consistency tests...");
        } else {
            fmt::println("=== Automatic GA Rule Generation - Rule Display ===");
            fmt::println("Tip: Use '--test_consistency' argument to run comparison testing");
        }

        // Test results (only used in test mode)
        std::vector<bool> test_results;

        // Configure and test EGA2D algebra
        AlgebraConfig ega2d_config = {
            .basis_vectors = {"e1", "e2"},
            .metric_signature = {+1, +1},
            .multivector_basis = {"1", "e1", "e2",
                                  "e12"}, // Copy from ga_prdxpr_ega2d.hpp line 12
            .scalar_name = one_str,
            .basis_prefix = "e"};

        if (test_consistency) {
            bool ega2d_success =
                test_algebra(ega2d_config, "EGA2D", mv2d_basis, gpr_ega2d_rules,
                             wdg_ega2d_rules, dot_ega2d_rules);
            test_results.push_back(ega2d_success);
        } else {
            display_algebra_rules(ega2d_config, "EGA2D");
        }

        // Configure and test EGA3D algebra
        AlgebraConfig ega3d_config = {
            .basis_vectors = {"e1", "e2", "e3"},
            .metric_signature = {+1, +1, +1},
            .multivector_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12",
                                  "e123"}, // Copy from ga_prdxpr_ega3d.hpp line 12
            .scalar_name = one_str,
            .basis_prefix = "e"};

        if (test_consistency) {
            bool ega3d_success =
                test_algebra(ega3d_config, "EGA3D", mv3d_basis, gpr_ega3d_rules,
                             wdg_ega3d_rules, dot_ega3d_rules);
            test_results.push_back(ega3d_success);
        } else {
            display_algebra_rules(ega3d_config, "EGA3D");
        }

        // Configure and test PGA2DP algebra
        AlgebraConfig pga2dp_config = {
            .basis_vectors = {"e1", "e2", "e3"},
            .metric_signature = {+1, +1, 0}, // G(2,0,1) - e1²=+1, e2²=+1, e3²=0
            .multivector_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12",
                                  "e321"}, // Copy from ga_prdxpr_pga2dp.hpp line 12
            .scalar_name = one_str,
            .basis_prefix = "e"};

        if (test_consistency) {
            bool pga2dp_success =
                test_algebra(pga2dp_config, "PGA2DP", mv2dp_basis, gpr_pga2dp_rules,
                             wdg_pga2dp_rules, dot_pga2dp_rules);
            test_results.push_back(pga2dp_success);
        } else {
            display_algebra_rules(pga2dp_config, "PGA2DP");
        }

        // Configure and test PGA3DP algebra
        AlgebraConfig pga3dp_config = {
            .basis_vectors = {"e1", "e2", "e3", "e4"},
            .metric_signature = {+1, +1, +1,
                                 0}, // G(3,0,1) - e1²=+1, e2²=+1, e3²=+1, e4²=0
            .multivector_basis = {"1", "e1", "e2", "e3", "e4", "e41", "e42", "e43", "e23",
                                  "e31", "e12", "e423", "e431", "e412", "e321",
                                  "e1234"}, // Copy from ga_prdxpr_pga3dp.hpp line 12-14
            .scalar_name = one_str,
            .basis_prefix = "e"};

        if (test_consistency) {
            bool pga3dp_success =
                test_algebra(pga3dp_config, "PGA3DP", mv3dp_basis, gpr_pga3dp_rules,
                             wdg_pga3dp_rules, dot_pga3dp_rules);
            test_results.push_back(pga3dp_success);
        } else {
            display_algebra_rules(pga3dp_config, "PGA3DP");
        }

        // Configure and test STA3D algebra (Space-Time Algebra)
        AlgebraConfig sta3d_config = {
            .basis_vectors = {"g0", "g1", "g2", "g3"},
            .metric_signature = {+1, -1, -1,
                                 -1}, // G(1,3,0) - g0²=+1, g1²=-1, g2²=-1, g3²=-1
            .multivector_basis = {"1", "g0", "g1", "g2", "g3", "g01", "g02", "g03", "g23",
                                  "g31", "g12", "g023", "g031", "g012", "g123",
                                  "g0123"}, // Copy from ga_prdxpr_sta3d.hpp line 16-18
            .scalar_name = one_str,
            .basis_prefix = "g" // Uses gamma notation instead of e
        };

        if (test_consistency) {
            bool sta3d_success = test_algebra(
                sta3d_config, "STA3D", mvsta3d_basis, gpr_sta3d_rules, wdg_sta3d_rules,
                dot_sta3d_rules); // Testing all three product types - complete manual
                                  // implementation
            test_results.push_back(sta3d_success);
        } else {
            display_algebra_rules(sta3d_config, "STA3D");
        }

        if (test_consistency) {
            // Overall summary for test mode
            std::string separator(80, '=');
            fmt::println("\n{}", separator);
            fmt::println("MULTI-ALGEBRA TEST SUMMARY");
            fmt::println("{}", separator);
            fmt::println("EGA2D (G(2,0,0)): {}", test_results[0] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("EGA3D (G(3,0,0)): {}", test_results[1] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("PGA2DP (G(2,0,1)): {}", test_results[2] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("PGA3DP (G(3,0,1)): {}", test_results[3] ? "✓ PERFECT" : "✗ FAILED");
            fmt::println("STA3D (G(1,3,0)): {}", test_results[4] ? "✓ PERFECT" : "✗ FAILED");

            bool all_success = std::all_of(test_results.begin(), test_results.end(),
                                           [](bool b) { return b; });
            fmt::println("\nOVERALL RESULT: {}",
                         all_success ? "✓ ALL ALGEBRAS PERFECT" : "✗ SOME ALGEBRAS FAILED");

            if (all_success) {
                fmt::println("\n🎉 SUCCESS: Generated rules are character-identical to "
                             "reference implementations!");
                fmt::println("The automatic rule generation system is working perfectly for "
                             "all tested algebras.");
            }
            else {
                fmt::println("\n❌ Issues found in rule generation. Please check the "
                             "differences above.");
                return 1;
            }
        } else {
            // Simple completion message for display mode
            std::string separator(80, '=');
            fmt::println("\n{}", separator);
            fmt::println("RULE DISPLAY COMPLETE");
            fmt::println("{}", separator);
            fmt::println("All generated rules have been displayed for visual inspection.");
            fmt::println("The rules are generated algorithmically and available for use.");
        }

        return 0;
    }
    catch (const std::exception& e) {
        fmt::println("Error: {}", e.what());
        return 1;
    }
}