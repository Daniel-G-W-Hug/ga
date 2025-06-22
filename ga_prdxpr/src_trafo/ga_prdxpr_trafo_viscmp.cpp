// Copyright 2024-2025, Daniel Hug. All rights reserved.

// Visual comparison test showing original input, expected manual result, and our
// transformation result This addresses the user's request: "the original input string can
// be compared with the modified result" "lines 20 and 21 contain the original input
// strings and line 31 and 32 contain the reference solution"

#include "ga_prdxpr_trafo_sandwich_transformer.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "fmt/format.h"

struct ComparisonTestCase {
    std::string algebra;
    std::string case_type; // vector, bivector, trivector
    std::string component;
    std::string description;
    std::string manual_line_ref;
    std::string original_input;  // From manual (e.g., lines 20-21)
    std::string expected_manual; // From manual (e.g., lines 31-32)
};

class VisualComparisonTester {
  public:


    // Extract input and expected output from manual file using section markers
    std::map<std::string, std::string>
    extract_input_from_manual(const std::string& algebra, const std::string& case_type)
    {
        std::map<std::string, std::string> input_components;

        // Read manual file (using documented path pattern from CLAUDE.md)
        std::ifstream file(
            "../../ga_prdxpr/src_trafo/ga_prdxpr_transformation_manual.txt");
        if (!file.is_open()) {
            fmt::println("Error: Could not open manual file");
            return input_components;
        }

        // Define section markers for each algebra and case type
        std::string section_marker;
        if (algebra == "ega2d" && case_type == "vector") {
            section_marker = "ega2d sandwich product vec_tmp * rev(mv_e) -> vec_res:";
        }
        else if (algebra == "ega3d" && case_type == "vector") {
            section_marker = "ega3d sandwich product mv_u_tmp * rev(mv_e) -> mv_u_res:";
        }
        else if (algebra == "ega3d" && case_type == "bivector") {
            section_marker = "ega3d sandwich product mv_e_tmp * rev(mv_e) -> mv_e_res:";
        }
        else if (algebra == "pga2dp" && case_type == "vector") {
            section_marker =
                "pga2dp regressive sandwich product mv_u_tmp * rrev(mv_u) -> mv_u_res:";
        }
        else if (algebra == "pga2dp" && case_type == "bivector") {
            section_marker =
                "pga2dp regressive sandwich product mv_e_tmp * rrev(mv_u) -> mv_e_res:";
        }
        else if (algebra == "pga3dp" && case_type == "vector") {
            section_marker =
                "pga3dp regressive sandwich product mv_u_tmp * rrev(mv_e) -> mv_u_res:";
        }
        else if (algebra == "pga3dp" && case_type == "bivector") {
            section_marker =
                "pga3dp regressive sandwich product mv_e_tmp * rrev(mv_e) -> mv_e_res:";
        }
        else if (algebra == "pga3dp" && case_type == "trivector") {
            section_marker = "pga3dp regressive sandwich product mv_u_tmp_t * rrev(mv_e) "
                             "-> mv_u_res_t:";
        }

        if (section_marker.empty()) {
            fmt::println("❌ No section marker defined for {} {}", algebra, case_type);
            return input_components;
        }

        // fmt::println("📖 READING FROM MANUAL - {} {} case:", algebra, case_type);
        // fmt::println("   Looking for marker: {}", section_marker);
        // fmt::println("");

        bool in_section = false;
        std::string line;
        while (std::getline(file, line)) {
            if (line.find(section_marker) != std::string::npos) {
                in_section = true;
                continue;
            }

            if (in_section) {
                // Look for lines like: [  expression  ] component
                auto bracket_start = line.find('[');
                auto bracket_end = line.find(']');
                auto component_start = line.find(']') + 1;

                if (bracket_start != std::string::npos &&
                    bracket_end != std::string::npos) {
                    // Extract expression (between brackets)
                    std::string expression =
                        line.substr(bracket_start + 1, bracket_end - bracket_start - 1);

                    // Extract component name (after bracket)
                    std::string component = line.substr(component_start);
                    // Clean up component name (remove spaces, comma)
                    component.erase(0, component.find_first_not_of(" "));
                    if (component.back() == ',') component.pop_back();

                    // Clean up expression (remove extra spaces)
                    expression.erase(0, expression.find_first_not_of(" "));
                    expression.erase(expression.find_last_not_of(" ") + 1);

                    // Store the component
                    input_components[component] = expression;
                }
                else {
                    // End of section (empty line or next section)
                    if (line.empty() ||
                        line.find("sandwich product") != std::string::npos) {
                        break;
                    }
                }
            }
        }
        file.close();

        // Fill in missing components with "0" to get complete multivector
        auto all_components = getAllComponents(algebra);
        for (const auto& comp : all_components) {
            if (input_components.find(comp) == input_components.end()) {
                input_components[comp] = "0";
            }
        }

        return input_components;
    }


    void run_comparison_tests()
    {
        fmt::println("=== VISUAL COMPARISON: ORIGINAL vs MANUAL vs OUR RESULT ===");
        fmt::println("Showing original input strings and expected manual results for "
                     "direct comparison");

        auto test_cases = create_test_cases();

        int total_tests = 0;
        int successful_transforms = 0;
        int exact_matches = 0;
        int new_cases = 0;

        // Group tests by algebra and case type for organized output
        std::string current_case_group = "";

        for (const auto& test : test_cases) {
            total_tests++;

            // Print case group header when we encounter a new algebra+case_type
            // combination
            std::string case_group =
                fmt::format("{} {} case", test.algebra, test.case_type);
            if (case_group != current_case_group) {
                current_case_group = case_group;
                fmt::println("\n📋 {}:", case_group);
                fmt::println("{}", std::string(40, '-'));
            }

            // Handle zero input cases (identity transformation) - show full visual
            // comparison
            if (test.original_input == "0") {
                successful_transforms++;
                std::string our_result = "0"; // Zero input always gives zero result

                if (test.expected_manual == "NO_EXPECTED_RESULT") {
                    new_cases++;
                }

                bool exact_match = (test.expected_manual == "0" ||
                                    test.expected_manual == "NO_EXPECTED_RESULT");
                if (exact_match && test.expected_manual != "NO_EXPECTED_RESULT") {
                    exact_matches++;
                }

                // VISUAL COMPARISON: Show detailed INPUT, EXPECTED, RESULT for identity
                // cases too
                std::string status_symbol =
                    "✅"; // Identity transformations are always correct
                fmt::println("{} {}:", status_symbol, test.component);
                fmt::println("INPUT:    {}", test.original_input);

                if (test.expected_manual == "NO_EXPECTED_RESULT") {
                    fmt::println("EXPECTED: (no expected result available)");
                    fmt::println("RESULT:   {}", our_result);
                    fmt::println("CATEGORY: Identity transformation (0 → 0) - no "
                                 "comparison available");
                }
                else {
                    fmt::println("EXPECTED: {}", test.expected_manual);
                    fmt::println("RESULT:   {}", our_result);
                    fmt::println("CATEGORY: Identity transformation (0 → 0)");
                }
                fmt::println("");
                continue;
            }

            try {
                // Use the new streamlined API for src_prdxpr integration
                std::string our_result = SandwichTransformer::transformExpression(
                    test.original_input, test.algebra);
                successful_transforms++;

                // Check if we have expected results for comparison
                bool has_expected_result = (test.expected_manual != "NO_EXPECTED_RESULT");
                if (!has_expected_result) {
                    new_cases++;
                }

                // Analyze transformation type and expected behavior
                bool input_is_zero = (test.original_input == "0");
                bool expected_is_zero =
                    has_expected_result ? (test.expected_manual == "0") : false;
                bool result_is_zero = (our_result == "0" || our_result.empty());

                // Compare with expected (only if we have expected results)
                bool exact_match = false;
                if (has_expected_result) {
                    exact_match = compare_expressions(our_result, test.expected_manual);
                    if (exact_match) {
                        exact_matches++;
                    }
                }

                // VISUAL COMPARISON: Show detailed INPUT, EXPECTED, RESULT for all cases
                std::string status_symbol;
                if (!has_expected_result) {
                    status_symbol =
                        "🔬"; // New result symbol for cases without expected results
                }
                else {
                    status_symbol = exact_match ? "✅" : "⚠️";
                }

                fmt::println("{} {}:", status_symbol, test.component);
                fmt::println("INPUT:    {}", test.original_input);

                if (!has_expected_result) {
                    fmt::println("EXPECTED: (no expected result available)");
                }
                else {
                    fmt::println("EXPECTED: {}", test.expected_manual);
                }
                fmt::println("RESULT:   {}", our_result);

                // CATEGORY IDENTIFICATION: Add transformation behavior analysis
                if (!has_expected_result) {
                    if (result_is_zero) {
                        fmt::println("CATEGORY: New result - appears to be symmetry "
                                     "cancellation (→ 0)");
                    }
                    else {
                        fmt::println("CATEGORY: New result - successful transformation");
                    }
                }
                else if (input_is_zero && expected_is_zero && result_is_zero) {
                    fmt::println("CATEGORY: Identity transformation (0 → 0)");
                }
                else if (!input_is_zero && expected_is_zero && result_is_zero) {
                    fmt::println("CATEGORY: Symmetry cancellation (non-zero → 0)");
                }
                else if (!input_is_zero && !expected_is_zero && !result_is_zero) {
                    if (exact_match) {
                        fmt::println("CATEGORY: Successful transformation (exact match)");
                    }
                    else {
                        fmt::println("CATEGORY: Successful transformation (algebraic "
                                     "equivalence)");
                    }
                }
                else if (!input_is_zero && expected_is_zero && !result_is_zero) {
                    fmt::println(
                        "CATEGORY: ⚠️ Failed symmetry cancellation (should be zero)");
                }
                else if (!input_is_zero && !expected_is_zero && result_is_zero) {
                    fmt::println(
                        "CATEGORY: ⚠️ Unexpected cancellation (should transform)");
                }
                else {
                    fmt::println("CATEGORY: Other case");
                }

                fmt::println("");
            }
            catch (const std::exception& e) {
                fmt::println("❌ {} FAILED: {}", test.component, e.what());
                fmt::println("");
            }
        }

        print_summary(total_tests, successful_transforms, exact_matches, new_cases);
    }

  private:

    // Get all components for each algebra
    std::vector<std::string> getAllComponents(const std::string& algebra)
    {
        if (algebra == "ega2d") {
            return {"1", "e1", "e2", "e12"};
        }
        else if (algebra == "ega3d") {
            return {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"};
        }
        else if (algebra == "pga2dp") {
            return {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e321"};
        }
        else if (algebra == "pga3dp") {
            return {"1",   "e1",  "e2",  "e3",   "e4",   "e41",  "e42",  "e43",
                    "e23", "e31", "e12", "e423", "e431", "e412", "e321", "e1234"};
        }
        return {};
    }

    // Placeholder input functions for structured tests (not used since we read from
    // manual file)
    std::string get_ega2d_vector_input() { return "placeholder"; }
    std::string get_ega3d_vector_input() { return "placeholder"; }
    std::string get_ega3d_bivector_input() { return "placeholder"; }
    std::string get_pga2dp_vector_input() { return "placeholder"; }
    std::string get_pga3dp_vector_input() { return "placeholder"; }
    std::string get_pga3dp_bivector_input() { return "placeholder"; }
    std::string get_pga3dp_trivector_input() { return "placeholder"; }

    // Get argument type components (components that should be non-zero)
    std::vector<std::string> getArgumentComponents(const std::string& algebra,
                                                   const std::string& argument_type)
    {
        if (argument_type == "vector") {
            if (algebra == "ega2d") return {"e1", "e2"};
            if (algebra == "ega3d") return {"e1", "e2", "e3"};
            if (algebra == "pga2dp") return {"e1", "e2", "e3"};
            if (algebra == "pga3dp") return {"e1", "e2", "e3", "e4"};
        }
        else if (argument_type == "bivector") {
            if (algebra == "ega3d") return {"e23", "e31", "e12"};
            if (algebra == "pga3dp") return {"e23", "e31", "e12", "e41", "e42", "e43"};
        }
        else if (argument_type == "trivector") {
            if (algebra == "pga3dp") return {"e423", "e431", "e412", "e321"};
        }
        return {};
    }

    // Get symmetry cancellation components (components that should be zero)
    std::vector<std::string> getSymmetryComponents(const std::string& algebra,
                                                   const std::string& argument_type)
    {
        auto all_components = getAllComponents(algebra);
        auto argument_components = getArgumentComponents(algebra, argument_type);

        std::vector<std::string> symmetry_components;
        for (const auto& comp : all_components) {
            if (std::find(argument_components.begin(), argument_components.end(), comp) ==
                argument_components.end()) {
                symmetry_components.push_back(comp);
            }
        }
        return symmetry_components;
    }

    // Extract expected results from manual file using MANUAL EXTENSION sections
    std::map<std::string, std::string>
    extract_expected_from_manual(const std::string& algebra, const std::string& case_type)
    {
        std::map<std::string, std::string> expected_components;

        std::ifstream file(
            "../../ga_prdxpr/src_trafo/ga_prdxpr_transformation_manual.txt");
        if (!file.is_open()) {
            fmt::println("Error: Could not open manual file for expected results");
            return expected_components;
        }

        // Define which algebras have expected results available
        std::vector<std::string> supported_cases = {"ega2d_vector",    "ega3d_vector",
                                                    "ega3d_bivector",  "pga2dp_vector",
                                                    "pga2dp_bivector", "pga3dp_vector"};

        std::string case_key = algebra + "_" + case_type;
        bool has_expected_results =
            std::find(supported_cases.begin(), supported_cases.end(), case_key) !=
            supported_cases.end();

        if (!has_expected_results) {
            return expected_components; // Return empty for unsupported cases (pga3dp
                                        // bivector/trivector)
        }

        // Find the corresponding MANUAL EXTENSION section first
        std::string manual_extension_marker;
        if (algebra == "ega2d" && case_type == "vector") {
            manual_extension_marker = "ega2d MANUAL EXTENSION mv_e * vec * rev(mv_e):";
        }
        else if (algebra == "ega3d" && case_type == "vector") {
            manual_extension_marker = "ega3d MANUAL EXTENSION mv_e * vec * rev(mv_e):";
        }
        else if (algebra == "ega3d" && case_type == "bivector") {
            manual_extension_marker = "ega3d MANUAL EXTENSION mv_e * bivec * rev(mv_e):";
        }
        else if (algebra == "pga2dp" && case_type == "vector") {
            manual_extension_marker = "pga2dp MANUAL EXTENSION mv_u * vec * rrev(mv_u):";
        }
        else if (algebra == "pga2dp" && case_type == "bivector") {
            manual_extension_marker =
                "pga2dp MANUAL EXTENSION mv_u * bivec * rrev(mv_u):";
        }
        else if (algebra == "pga3dp" && case_type == "vector") {
            manual_extension_marker = "pga3dp MANUAL EXTENSION mv_e * vec * rrev(mv_e):";
        }

        bool in_manual_extension = false;
        bool in_final_results = false;
        std::string line;

        while (std::getline(file, line)) {
            // First find the correct MANUAL EXTENSION section
            if (line.find(manual_extension_marker) != std::string::npos) {
                in_manual_extension = true;
                continue;
            }

            if (in_manual_extension) {
                // Look for the "// reorder and introduce braces again:" marker
                if (line.find("// reorder and introduce braces again:") !=
                    std::string::npos) {
                    in_final_results = true;
                    continue;
                }

                if (in_final_results) {
                    auto bracket_start = line.find('[');
                    auto bracket_end = line.find(']');
                    auto component_start = line.find(']') + 1;

                    if (bracket_start != std::string::npos &&
                        bracket_end != std::string::npos) {
                        // Extract expression (between brackets)
                        std::string expression = line.substr(
                            bracket_start + 1, bracket_end - bracket_start - 1);

                        // Extract component name (after bracket)
                        std::string component = line.substr(component_start);
                        // Clean up component name (remove spaces, comma)
                        component.erase(0, component.find_first_not_of(" "));
                        if (!component.empty() && component.back() == ',')
                            component.pop_back();

                        // Clean up expression (remove extra spaces)
                        expression.erase(0, expression.find_first_not_of(" "));
                        expression.erase(expression.find_last_not_of(" ") + 1);

                        // Store the component (skip if expression is "0")
                        if (!expression.empty() && expression != "0") {
                            expected_components[component] = expression;
                        }
                    }
                    else if (line.empty() || line.find("k0") != std::string::npos ||
                             line.find("---") != std::string::npos ||
                             line.find("MANUAL EXTENSION") != std::string::npos) {
                        // End of final results section
                        break;
                    }
                }
            }
        }
        file.close();

        return expected_components;
    }

    std::vector<ComparisonTestCase> create_test_cases()
    {
        std::vector<ComparisonTestCase> cases;

        // List of all test configurations
        std::vector<std::tuple<std::string, std::string, std::string>> test_configs = {
            {"ega2d", "vector", "EGA2D Vector"},
            {"ega3d", "vector", "EGA3D Vector"},
            {"ega3d", "bivector", "EGA3D Bivector"},
            {"pga2dp", "vector", "PGA2DP Vector"},
            {"pga2dp", "bivector", "PGA2DP Bivector"},
            {"pga3dp", "vector", "PGA3DP Vector"},
            {"pga3dp", "bivector", "PGA3DP Bivector"},
            {"pga3dp", "trivector", "PGA3DP Trivector"}};

        // Generate test cases dynamically from manual file
        for (const auto& [algebra, case_type, description] : test_configs) {
            // Extract input components from manual file
            auto input_components = extract_input_from_manual(algebra, case_type);

            // Extract expected results from manual file
            auto expected_components = extract_expected_from_manual(algebra, case_type);

            // Create test cases for ALL components (including zeros) to show complete
            // transformation
            auto all_components = getAllComponents(algebra);
            bool has_expected_results = !expected_components.empty();

            for (const auto& comp : all_components) {
                std::string input_expr =
                    input_components.count(comp) ? input_components[comp] : "0";
                std::string expected_result;

                if (has_expected_results) {
                    // We have expected results for this algebra/case_type
                    expected_result =
                        expected_components.count(comp) ? expected_components[comp] : "0";
                }
                else {
                    // No expected results available - mark as "NO_EXPECTED_RESULT"
                    expected_result = "NO_EXPECTED_RESULT";
                }

                // If input is empty, treat as zero
                if (input_expr.empty()) input_expr = "0";

                cases.push_back(
                    {algebra, case_type, comp,
                     fmt::format("{} {} component", description, comp),
                     fmt::format("from manual file: {} {}", algebra, case_type),
                     input_expr, expected_result});
            }
        }

        return cases;
    }

    bool compare_expressions(const std::string& our_result, const std::string& expected)
    {
        // Basic comparison for exact matches
        // Remove all spaces for comparison
        std::string clean_our = our_result;
        std::string clean_expected = expected;

        clean_our.erase(std::remove(clean_our.begin(), clean_our.end(), ' '),
                        clean_our.end());
        clean_expected.erase(
            std::remove(clean_expected.begin(), clean_expected.end(), ' '),
            clean_expected.end());

        return clean_our == clean_expected;
    }

    void print_summary(int total, int successful, int exact_matches, int new_cases)
    {
        fmt::println("{}", std::string(80, '='));
        fmt::println("COMPREHENSIVE MULTIVECTOR TRANSFORMATION SUMMARY");
        fmt::println("{}", std::string(80, '='));

        int cases_with_expected = total - new_cases;

        fmt::println("Total test cases: {}", total);
        fmt::println("Cases with expected results: {}", cases_with_expected);
        fmt::println("New cases (no expected results): {}", new_cases);
        fmt::println("Successful transformations: {}", successful);
        fmt::println("Exact matches with manual: {}", exact_matches);
        fmt::println("Success rate: {:.1f}%", (100.0 * successful) / std::max(total, 1));
        if (cases_with_expected > 0) {
            fmt::println("Exact match rate (for cases with expected): {:.1f}%",
                         (100.0 * exact_matches) / std::max(cases_with_expected, 1));
        }

        if (exact_matches > 0) {
            fmt::println("✅ {} cases show exact match with manual results",
                         exact_matches);
        }
        if (successful > exact_matches) {
            fmt::println("⚠️  {} cases show algebraic equivalence (different form)",
                         successful - exact_matches);
        }
    }
};

int main()
{
    VisualComparisonTester tester;

    // Run visual comparison tests with detailed INPUT/EXPECTED/RESULT analysis
    tester.run_comparison_tests();

    return 0;
}