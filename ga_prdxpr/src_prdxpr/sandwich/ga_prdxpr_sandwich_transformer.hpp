#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "sandwich/ga_prdxpr_sandwich_simplifier.hpp"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// Main transformation interface for ga_prdxpr integration
///////////////////////////////////////////////////////////////////////////////

class SandwichTransformer {
  public:

    struct MatrixTransformation {
        // Original expression from ga_prdxpr
        std::string input_expression;

        // Simplified matrix coefficients
        // Key: result_component (e.g., "e1", "e2", "e23", etc.)
        // Value: map of geometric_variable -> coefficient_expression
        std::map<std::string, std::map<std::string, std::string>> matrix_coefficients;

        // Generated matrix form expressions
        std::map<std::string, std::string> simplified_expressions;

        // N-ary AST result (exact format matching)
        std::string nary_ast_result;

        // Matrix dimensions determined from algebra
        size_t result_dimension;
        size_t input_dimension;

        // Verification: can be used to check correctness
        bool coefficients_isolated;
        bool symmetric_terms_cancelled;
    };

    // Main transformation interface
    // Input: Raw sandwich product expression from ga_prdxpr
    // Output: Matrix transformation with simplified coefficients
    static MatrixTransformation transformSandwichProduct(
        std::string const& sandwich_expression,
        std::string const& algebra_type = "ega2d" // ega2d, ega3d, pga2dp, pga3dp, sta4ds
    );

    // Main transformation interface with custom patterns
    static MatrixTransformation
    transformSandwichProduct(std::string const& sandwich_expression,
                             std::string const& algebra_type,
                             GeometricVariablePatterns const& patterns);

    // Simple string-to-string interface for src_prdxpr integration
    static std::string transformExpression(std::string const& expression,
                                           std::string const& algebra_type = "ega2d");

    // Simple string-to-string interface with custom patterns
    static std::string transformExpression(std::string const& expression,
                                           std::string const& algebra_type,
                                           GeometricVariablePatterns const& patterns);

    // Transform complete multivector result (vector of component expressions)
    static std::vector<std::string>
    transformSandwichMultivector(std::vector<std::string> const& component_expressions,
                                 std::string const& algebra_type = "ega2d");

    // Batch processing for multiple expressions (existing interface)
    static std::vector<MatrixTransformation>
    transformMultipleExpressions(std::vector<std::string> const& expressions,
                                 std::string const& algebra_type = "ega2d");

    // Generate final matrix form string suitable for code generation
    static std::string
    generateMatrixFormExpression(MatrixTransformation const& transform);

    // Generate optimized coefficient expressions (with k-coefficients)
    static std::string
    generateOptimizedCoefficients(MatrixTransformation const& transform);

  private:

    // Core transformation pipeline
    static std::pair<std::vector<SimplifiedTerm>, std::shared_ptr<ast_node>>
    parseAndSimplify(std::string const& expression);
    static void extractMatrixCoefficients(std::vector<SimplifiedTerm> const& terms,
                                          MatrixTransformation& result,
                                          std::string const& algebra_type);
    static void generateSimplifiedExpressions(MatrixTransformation& result,
                                              std::string const& algebra_type);

    // Algebra-specific patterns
    static std::vector<std::string>
    getGeometricVariables(std::string const& algebra_type);
    static std::vector<std::string> getResultComponents(std::string const& algebra_type);
    static std::string extractRotorCoefficient(SimplifiedTerm const& term,
                                               std::string const& geometric_var);

    // Coefficient manipulation
    static std::string combineCoefficients(std::vector<std::string> const& coeffs);
    static std::string optimizeCoefficient(std::string const& coeff_expr);
};

///////////////////////////////////////////////////////////////////////////////
// Algebra-specific configuration
///////////////////////////////////////////////////////////////////////////////

struct SandwichAlgebraConfig {
    std::string name;
    std::vector<std::string> geometric_variables; // e.g., {"v.x", "v.y"} for vectors
    std::vector<std::string> result_components;   // e.g., {"e1", "e2"} for vector results
    std::vector<std::string> rotor_coefficients;  // e.g., {"R.c0", "R.c1"} for EGA2D
    size_t matrix_size;
};

class AlgebraRegistry {
  public:

    static SandwichAlgebraConfig getConfig(std::string const& algebra_type);

  private:

    static SandwichAlgebraConfig createEGA2DConfig();
    static SandwichAlgebraConfig createEGA3DConfig();
    static SandwichAlgebraConfig createPGA2DPConfig();
    static SandwichAlgebraConfig createPGA3DPConfig();
    static SandwichAlgebraConfig createSTA4DSConfig();
    static SandwichAlgebraConfig createCGA2DCConfig();
    static SandwichAlgebraConfig createCGA3DCConfig();
};

///////////////////////////////////////////////////////////////////////////////
// Integration test cases
///////////////////////////////////////////////////////////////////////////////

class TransformationTests {
  public:

    // Test EGA2D vector sandwich (matches manual example)
    static void testEGA2DVectorTransformation();

    // Test EGA3D vector sandwich
    static void testEGA3DVectorTransformation();

    // Test coefficient isolation and matrix generation
    static void testCoefficientExtraction();

    // Test complete pipeline: ga_prdxpr output -> simplified matrix form
    static void testCompletePipeline();

    static void runAllTests();
};
