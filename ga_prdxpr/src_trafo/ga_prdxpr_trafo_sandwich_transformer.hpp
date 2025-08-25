#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo_expr_simplifier.hpp"

#include <map>
#include <string>
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
        const std::string& sandwich_expression,
        const std::string& algebra_type = "ega2d" // ega2d, ega3d, pga2dp, pga3dp
    );
    
    // Main transformation interface with custom patterns
    static MatrixTransformation transformSandwichProduct(
        const std::string& sandwich_expression,
        const std::string& algebra_type,
        const GeometricVariablePatterns& patterns
    );

    // Simple string-to-string interface for src_prdxpr integration
    static std::string transformExpression(const std::string& expression,
                                           const std::string& algebra_type = "ega2d");
                                           
    // Simple string-to-string interface with custom patterns
    static std::string transformExpression(const std::string& expression,
                                           const std::string& algebra_type,
                                           const GeometricVariablePatterns& patterns);

    // Transform complete multivector result (vector of component expressions)
    static std::vector<std::string>
    transformSandwichMultivector(const std::vector<std::string>& component_expressions,
                                 const std::string& algebra_type = "ega2d");

    // Batch processing for multiple expressions (existing interface)
    static std::vector<MatrixTransformation>
    transformMultipleExpressions(const std::vector<std::string>& expressions,
                                 const std::string& algebra_type = "ega2d");

    // Generate final matrix form string suitable for code generation
    static std::string
    generateMatrixFormExpression(const MatrixTransformation& transform);

    // Generate optimized coefficient expressions (with k-coefficients)
    static std::string
    generateOptimizedCoefficients(const MatrixTransformation& transform);

  private:

    // Core transformation pipeline
    static std::pair<std::vector<SimplifiedTerm>, std::shared_ptr<ast_node>>
    parseAndSimplify(const std::string& expression);
    static void extractMatrixCoefficients(const std::vector<SimplifiedTerm>& terms,
                                          MatrixTransformation& result,
                                          const std::string& algebra_type);
    static void generateSimplifiedExpressions(MatrixTransformation& result,
                                              const std::string& algebra_type);

    // Algebra-specific patterns
    static std::vector<std::string>
    getGeometricVariables(const std::string& algebra_type);
    static std::vector<std::string> getResultComponents(const std::string& algebra_type);
    static std::string extractRotorCoefficient(const SimplifiedTerm& term,
                                               const std::string& geometric_var);

    // Coefficient manipulation
    static std::string combineCoefficients(const std::vector<std::string>& coeffs);
    static std::string optimizeCoefficient(const std::string& coeff_expr);
};

///////////////////////////////////////////////////////////////////////////////
// Algebra-specific configuration
///////////////////////////////////////////////////////////////////////////////

struct AlgebraConfig {
    std::string name;
    std::vector<std::string> geometric_variables; // e.g., {"v.x", "v.y"} for vectors
    std::vector<std::string> result_components;   // e.g., {"e1", "e2"} for vector results
    std::vector<std::string> rotor_coefficients;  // e.g., {"R.c0", "R.c1"} for EGA2D
    size_t matrix_size;
};

class AlgebraRegistry {
  public:

    static AlgebraConfig getConfig(const std::string& algebra_type);

  private:

    static AlgebraConfig createEGA2DConfig();
    static AlgebraConfig createEGA3DConfig();
    static AlgebraConfig createPGA2DPConfig();
    static AlgebraConfig createPGA3DPConfig();
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