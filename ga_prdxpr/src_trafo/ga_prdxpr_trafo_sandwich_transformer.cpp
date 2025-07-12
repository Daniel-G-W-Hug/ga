// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo_sandwich_transformer.hpp"
#include "ga_prdxpr_trafo_nary_expression.hpp"

#include <algorithm>
#include <numeric>
#include <regex>
#include <sstream>

#include "fmt/format.h"

///////////////////////////////////////////////////////////////////////////////
// Main transformation interface implementation
///////////////////////////////////////////////////////////////////////////////

SandwichTransformer::MatrixTransformation
SandwichTransformer::transformSandwichProduct(const std::string& sandwich_expression,
                                              const std::string& algebra_type)
{

    MatrixTransformation result;
    result.input_expression = sandwich_expression;

    // Get algebra configuration
    auto config = AlgebraRegistry::getConfig(algebra_type);
    result.result_dimension = config.result_components.size();
    result.input_dimension = config.geometric_variables.size();

    // Core transformation pipeline
    auto [simplified_terms, simplified_ast] = parseAndSimplify(sandwich_expression);
    extractMatrixCoefficients(simplified_terms, result, algebra_type);

    // Store the N-ary AST result for exact format matching - capture before binary AST
    // conversion Parse and get the N-ary factorized result directly
    Parser parser(sandwich_expression);
    auto ast = parser.parse();
    auto expanded_ast = ExpressionSimplifier::expandProducts(ast);
    auto nary_expr = NAryConverter::fromBinaryAST(expanded_ast);
    nary_expr.normalizeSignsAndCommutativity();
    nary_expr.combineTermsAndRegroup();
    nary_expr.removeZeroTerms();
    result.nary_ast_result = nary_expr.toString();

    // Use the traditional approach for matrix form
    generateSimplifiedExpressions(result, algebra_type);

    // Validation
    result.coefficients_isolated = !result.matrix_coefficients.empty();
    result.symmetric_terms_cancelled = simplified_terms.size() < 10; // heuristic

    return result;
}

std::vector<SandwichTransformer::MatrixTransformation>
SandwichTransformer::transformMultipleExpressions(
    const std::vector<std::string>& expressions, const std::string& algebra_type)
{

    std::vector<MatrixTransformation> results;
    results.reserve(expressions.size());

    for (const auto& expr : expressions) {
        results.push_back(transformSandwichProduct(expr, algebra_type));
    }

    return results;
}

std::string
SandwichTransformer::generateMatrixFormExpression(const MatrixTransformation& transform)
{
    std::string result = "Matrix form transformation:\n";

    for (const auto& [result_comp, simplified_expr] : transform.simplified_expressions) {
        result += fmt::format("[{}] = {}\n", result_comp, simplified_expr);
    }

    return result;
}

std::string
SandwichTransformer::generateOptimizedCoefficients(const MatrixTransformation& transform)
{
    std::string result = "Optimized coefficients:\n";

    // Extract common sub-expressions like R.c0*R.c0, R.c0*R.c1, etc.
    std::map<std::string, std::string> k_coefficients;

    for (const auto& [result_comp, var_coeffs] : transform.matrix_coefficients) {
        for (const auto& [var, coeff] : var_coeffs) {

            // Find patterns like R.c0*R.c0, R.c0*R.c1
            std::regex pattern(R"(R\.c(\d+)\s*\*\s*R\.c(\d+))");
            std::smatch match;
            std::string optimized_coeff = coeff;

            auto iter = std::sregex_iterator(coeff.begin(), coeff.end(), pattern);
            auto end = std::sregex_iterator();

            for (; iter != end; ++iter) {
                match = *iter;
                std::string full_match = match.str();
                std::string k_name = fmt::format("k{}{}", match[1].str(), match[2].str());

                if (k_coefficients.find(k_name) == k_coefficients.end()) {
                    k_coefficients[k_name] = full_match;
                }
            }
        }
    }

    // Generate k-coefficient definitions
    for (const auto& [k_name, k_expr] : k_coefficients) {
        result += fmt::format("{} = {};\n", k_name, k_expr);
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////
// Core pipeline implementation
///////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<SimplifiedTerm>, std::shared_ptr<ast_node>>
SandwichTransformer::parseAndSimplify(const std::string& expression)
{
    // Parse expression using existing parser
    Parser parser(expression);
    auto ast = parser.parse();

    // Use the new N-ary simplification approach
    auto simplified_ast = NAryExpressionSimplifier::simplify(ast);

    // Convert final result to simplified terms (for compatibility with existing matrix
    // extraction)
    auto simplified_terms = ExpressionSimplifier::astToTerms(simplified_ast);

    return {simplified_terms, simplified_ast};
}

void SandwichTransformer::extractMatrixCoefficients(
    const std::vector<SimplifiedTerm>& terms, MatrixTransformation& result,
    const std::string& algebra_type)
{

    auto config = AlgebraRegistry::getConfig(algebra_type);

    // Group terms by geometric variables
    std::map<std::string, std::vector<SimplifiedTerm>> terms_by_variable;

    for (const auto& term : terms) {
        for (const auto& [var, power] : term.factors) {
            if (std::find(config.geometric_variables.begin(),
                          config.geometric_variables.end(),
                          var) != config.geometric_variables.end()) {
                terms_by_variable[var].push_back(term);
            }
        }
    }

    // For each result component, extract coefficients for each geometric variable
    for (const auto& result_comp : config.result_components) {
        for (const auto& geom_var : config.geometric_variables) {

            // Combine all terms for this geometric variable
            std::vector<std::string> coefficients;

            if (terms_by_variable.count(geom_var) > 0) {
                for (const auto& term : terms_by_variable[geom_var]) {
                    std::string coeff = extractRotorCoefficient(term, geom_var);
                    if (!coeff.empty()) {
                        coefficients.push_back(coeff);
                    }
                }
            }

            if (!coefficients.empty()) {
                result.matrix_coefficients[result_comp][geom_var] =
                    combineCoefficients(coefficients);
            }
            else {
                result.matrix_coefficients[result_comp][geom_var] = "0";
            }
        }
    }
}

void SandwichTransformer::generateSimplifiedExpressions(MatrixTransformation& result,
                                                        const std::string& algebra_type)
{
    auto config = AlgebraRegistry::getConfig(algebra_type);

    for (const auto& result_comp : config.result_components) {
        std::vector<std::string> terms;

        for (const auto& geom_var : config.geometric_variables) {
            auto coeff = result.matrix_coefficients[result_comp][geom_var];
            if (coeff != "0") {
                if (coeff == "1") {
                    terms.push_back(geom_var);
                }
                else {
                    terms.push_back(fmt::format("({}) * {}", coeff, geom_var));
                }
            }
        }

        if (terms.empty()) {
            result.simplified_expressions[result_comp] = "0";
        }
        else {
            // Combine terms with proper sign handling
            std::string expression = terms[0];
            for (size_t i = 1; i < terms.size(); ++i) {
                const std::string& term = terms[i];

                // Simply add the term as-is - parentheses are already correct from matrix
                // coefficients
                expression += " + " + term;
            }
            result.simplified_expressions[result_comp] = expression;
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

std::string SandwichTransformer::extractRotorCoefficient(const SimplifiedTerm& term,
                                                         const std::string& geometric_var)
{

    // Calculate the signed coefficient
    double signed_coeff = term.sign * term.numeric_coefficient;

    // Add rotor factors (everything except geometric_var)
    std::vector<std::string> rotor_factors;
    for (const auto& [var, power] : term.factors) {
        if (var != geometric_var) {
            // Generate C++-compatible multiplication instead of exponentiation
            for (int p = 0; p < power; ++p) {
                rotor_factors.push_back(var);
            }
        }
    }

    // Build the coefficient string with proper sign handling
    if (rotor_factors.empty()) {
        // Only numeric coefficient
        if (std::abs(signed_coeff - 1.0) < 1e-10) {
            return "1";
        }
        else if (std::abs(signed_coeff + 1.0) < 1e-10) {
            return "-1";
        }
        else {
            return fmt::format("{}", signed_coeff);
        }
    }
    else {
        // Combine numeric coefficient with rotor factors
        std::string factors_str = std::accumulate(
            rotor_factors.begin() + 1, rotor_factors.end(), rotor_factors[0],
            [](const std::string& a, const std::string& b) { return a + " * " + b; });

        if (std::abs(signed_coeff - 1.0) < 1e-10) {
            return factors_str;
        }
        else if (std::abs(signed_coeff + 1.0) < 1e-10) {
            return "-" + factors_str;
        }
        else {
            return fmt::format("{} * {}", signed_coeff, factors_str);
        }
    }
}

std::string
SandwichTransformer::combineCoefficients(const std::vector<std::string>& coeffs)
{
    if (coeffs.empty()) return "0";
    if (coeffs.size() == 1) return coeffs[0];

    // Combine terms with proper sign handling
    std::string result = coeffs[0];

    for (size_t i = 1; i < coeffs.size(); ++i) {
        const std::string& coeff = coeffs[i];

        // Check if the coefficient starts with a minus sign
        if (coeff.starts_with("-")) {
            // Already has minus sign, use subtraction
            result += " - " + coeff.substr(1);
        }
        else {
            // Positive coefficient, use addition
            result += " + " + coeff;
        }
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////
// New API functions for src_prdxpr integration
///////////////////////////////////////////////////////////////////////////////

std::string SandwichTransformer::transformExpression(const std::string& expression,
                                                     const std::string& algebra_type)
{

    if (expression == "0" || expression.empty()) {
        return "0";
    }

    // Use existing transformation pipeline with minimal overhead
    auto result = transformSandwichProduct(expression, algebra_type);
    return result.nary_ast_result;
}

std::vector<std::string> SandwichTransformer::transformSandwichMultivector(
    const std::vector<std::string>& component_expressions,
    const std::string& algebra_type)
{

    std::vector<std::string> transformed_results;
    transformed_results.reserve(component_expressions.size());

    for (const auto& expr : component_expressions) {
        transformed_results.push_back(transformExpression(expr, algebra_type));
    }

    return transformed_results;
}

///////////////////////////////////////////////////////////////////////////////
// Algebra configuration implementation
///////////////////////////////////////////////////////////////////////////////

AlgebraConfig AlgebraRegistry::getConfig(const std::string& algebra_type)
{
    if (algebra_type == "ega2d") {
        return createEGA2DConfig();
    }
    else if (algebra_type == "ega3d") {
        return createEGA3DConfig();
    }
    else if (algebra_type == "pga2dp") {
        return createPGA2DPConfig();
    }
    else if (algebra_type == "pga3dp") {
        return createPGA3DPConfig();
    }
    else {
        throw std::runtime_error("Unknown algebra type: " + algebra_type);
    }
}

AlgebraConfig AlgebraRegistry::createEGA2DConfig()
{
    return {.name = "ega2d",
            .geometric_variables = {"v.x", "v.y"},
            .result_components = {"e1", "e2"},
            .rotor_coefficients = {"R.c0", "R.c1"},
            .matrix_size = 2};
}

AlgebraConfig AlgebraRegistry::createEGA3DConfig()
{
    return {.name = "ega3d",
            .geometric_variables = {"v.x", "v.y", "v.z", "B.x", "B.y", "B.z"},
            .result_components = {"e1", "e2", "e3", "e23", "e31", "e12"},
            .rotor_coefficients = {"R.c0", "R.c1", "R.c2", "R.c3"},
            .matrix_size = 6};
}

AlgebraConfig AlgebraRegistry::createPGA2DPConfig()
{
    return {
        .name = "pga2dp",
        .geometric_variables = {"v.x", "v.y", "v.z"}, // Projective 2D has 3 coordinates
        .result_components = {"e1", "e2", "e3"},
        .rotor_coefficients = {"R.c0", "R.c1", "R.c2", "R.c3"},
        .matrix_size = 3};
}

AlgebraConfig AlgebraRegistry::createPGA3DPConfig()
{
    return {.name = "pga3dp",
            .geometric_variables = {"v.x", "v.y", "v.z",
                                    "v.w"}, // Projective 3D has 4 coordinates
            .result_components = {"e1", "e2", "e3", "e4"},
            .rotor_coefficients = {"R.c0", "R.c1", "R.c2", "R.c3", "R.c4", "R.c5", "R.c6",
                                   "R.c7"},
            .matrix_size = 4};
}

///////////////////////////////////////////////////////////////////////////////
// Test implementation
///////////////////////////////////////////////////////////////////////////////

void TransformationTests::testEGA2DVectorTransformation()
{
    fmt::println("=== Testing EGA2D Vector Transformation ===");

    // Test case from manual examples
    std::string input =
        "(R.c0 * v.x + R.c1 * v.y) * R.c0 + (R.c0 * v.y - R.c1 * v.x) * R.c1";

    auto result = SandwichTransformer::transformSandwichProduct(input, "ega2d");

    fmt::println("Input: {}", input);
    fmt::println("");

    // Show matrix coefficients
    fmt::println("Matrix coefficients:");
    for (const auto& [result_comp, var_coeffs] : result.matrix_coefficients) {
        for (const auto& [var, coeff] : var_coeffs) {
            fmt::println("M[{}][{}] = {}", result_comp, var, coeff);
        }
    }
    fmt::println("");

    // Show simplified expressions
    fmt::println("{}", SandwichTransformer::generateMatrixFormExpression(result));

    // Expected results verification
    fmt::println("Expected vs Actual:");
    fmt::println("Expected e1: (R.c0*R.c0 - R.c1*R.c1) * v.x + (2*R.c0*R.c1) * v.y");
    fmt::println("Actual   e1: {}", result.simplified_expressions["e1"]);
    fmt::println("");
}

void TransformationTests::testCompletePipeline()
{
    fmt::println("=== Testing Complete Pipeline ===");

    // Simulate multiple ga_prdxpr outputs
    std::vector<std::string> prdxpr_outputs = {
        "(R.c0 * v.x + R.c1 * v.y) * R.c0 + (R.c0 * v.y - R.c1 * v.x) * R.c1",
        "-(R.c0 * v.x + R.c1 * v.y) * R.c1 + (R.c0 * v.y - R.c1 * v.x) * R.c0"};

    auto results =
        SandwichTransformer::transformMultipleExpressions(prdxpr_outputs, "ega2d");

    for (size_t i = 0; i < results.size(); ++i) {
        fmt::println("Expression {}: {}", i + 1, prdxpr_outputs[i]);
        fmt::println("{}", SandwichTransformer::generateMatrixFormExpression(results[i]));
    }
}

void TransformationTests::testEGA3DVectorTransformation()
{
    fmt::println("=== Testing EGA3D Vector Transformation ===");

    // Test case from manual examples - e1 component
    // This is the expanded form from line 94 of the manual
    std::string input_e1 =
        "R.c0 * R.c0 * v.x + R.c1 * R.c1 * v.x - R.c2 * R.c2 * v.x - R.c3 * R.c3 * v.x + "
        "2 * R.c0 * R.c3 * v.y + 2 * R.c1 * R.c2 * v.y - "
        "2 * R.c0 * R.c2 * v.z + 2 * R.c1 * R.c3 * v.z";

    auto result = SandwichTransformer::transformSandwichProduct(input_e1, "ega3d");

    fmt::println("Input (e1 component): {}", input_e1);
    fmt::println("");

    // Show matrix coefficients
    fmt::println("Matrix coefficients extracted:");
    for (const auto& [result_comp, var_coeffs] : result.matrix_coefficients) {
        for (const auto& [var, coeff] : var_coeffs) {
            fmt::println("  M[{}][{}] = {}", result_comp, var, coeff);
        }
    }
    fmt::println("");

    // Show simplified expressions
    fmt::println("{}", SandwichTransformer::generateMatrixFormExpression(result));

    // Expected results verification
    fmt::println("Expected vs Actual:");
    fmt::println(
        "Expected e1: (R.c0*R.c0 + R.c1*R.c1 - R.c2*R.c2 - R.c3*R.c3) * v.x + "
        "(2*R.c0*R.c3 + 2*R.c1*R.c2) * v.y + (-2*R.c0*R.c2 + 2*R.c1*R.c3) * v.z");
    fmt::println("Actual   e1: {}", result.simplified_expressions["e1"]);
    fmt::println("");
}

void TransformationTests::runAllTests()
{
    testEGA2DVectorTransformation();
    testEGA3DVectorTransformation();
    testCompletePipeline();
}