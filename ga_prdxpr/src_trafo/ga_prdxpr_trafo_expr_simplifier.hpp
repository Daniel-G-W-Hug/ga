#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

// Forward declarations
class SimplifiedTerm;
class ExpressionSimplifier;

///////////////////////////////////////////////////////////////////////////////
// Simplified representation for algebraic manipulation
///////////////////////////////////////////////////////////////////////////////

// Represents a single term: coefficient * product_of_variables
// Example: "2 * R.c0 * R.c1 * v.x" or "-R.c2 * R.c2 * v.z"
class SimplifiedTerm {
  public:

    double numeric_coefficient;         // Numeric part (e.g., 2, -1, 0.5)
    std::map<std::string, int> factors; // Variable factors with powers
    int sign;                           // +1 or -1

    SimplifiedTerm();
    SimplifiedTerm(double coeff, const std::map<std::string, int>& vars, int s = 1);

    // Create from string representation
    static SimplifiedTerm fromString(const std::string& term_str);

    // Operations
    SimplifiedTerm operator*(const SimplifiedTerm& other) const;
    SimplifiedTerm operator+(const SimplifiedTerm& other) const; // Only if same variables
    bool canCombineWith(const SimplifiedTerm& other) const;
    bool isZero() const;
    bool isEquivalentTo(const SimplifiedTerm& other) const;

    // Convert back to string
    std::string toString() const;

  private:

    void normalize(); // Ensure canonical form
};

///////////////////////////////////////////////////////////////////////////////
// Expression simplification engine
///////////////////////////////////////////////////////////////////////////////

class ExpressionSimplifier {
  public:

    // Main simplification pipeline
    static std::shared_ptr<ast_node> simplify(std::shared_ptr<ast_node> ast);

    // Core simplification operations
    static std::shared_ptr<ast_node> expandProducts(std::shared_ptr<ast_node> ast);
    static std::shared_ptr<ast_node> combineTerms(std::shared_ptr<ast_node> ast);
    static std::shared_ptr<ast_node> applyCommutativity(std::shared_ptr<ast_node> ast);
    static std::shared_ptr<ast_node> cancelZeros(std::shared_ptr<ast_node> ast);

    // 3-step manual transformation process
    static std::shared_ptr<ast_node>
    normalizeSignsAndApplyCommutativity(std::shared_ptr<ast_node> ast);
    static std::shared_ptr<ast_node>
    combineTermsAndRegroup(std::shared_ptr<ast_node> ast);

    // Convert AST to simplified term representation
    static std::vector<SimplifiedTerm> astToTerms(std::shared_ptr<ast_node> ast);

    // Convert simplified terms back to AST
    static std::shared_ptr<ast_node> termsToAst(const std::vector<SimplifiedTerm>& terms);

    // Test functions for specific GA patterns
    static void testEGA2DExpansion();
    static void testSymmetryCancellation();
    static void testCommutativeReordering();

    // Public for testing
    static std::vector<SimplifiedTerm>
    combineTermsVector(const std::vector<SimplifiedTerm>& terms);

  private:

    // Helper functions for AST traversal and manipulation
    static void extractTermsFromExpression(std::shared_ptr<Expression> expr,
                                           std::vector<SimplifiedTerm>& terms,
                                           int sign = 1);
    static void extractTermsFromTerm(std::shared_ptr<Term> term,
                                     std::vector<SimplifiedTerm>& terms, int sign = 1);
    static void extractTermsFromFactor(std::shared_ptr<Factor> factor,
                                       SimplifiedTerm& current_term, int sign = 1);
    static std::vector<SimplifiedTerm> extractFactorTerms(std::shared_ptr<Factor> factor,
                                                          int sign = 1);

    // Algebraic simplification rules (moved to public)
    static SimplifiedTerm expandSingleTerm(const SimplifiedTerm& term);
    static bool
    isCommutativeVariable(const std::string& var); // R.c0, R.c1, etc. are commutative
    static SimplifiedTerm applyCommutativeReordering(const SimplifiedTerm& term);
};

///////////////////////////////////////////////////////////////////////////////
// Specific GA algebraic rules
///////////////////////////////////////////////////////////////////////////////

class GAAlgebraRules {
  public:

    // Comprehensive canonical ordering for all GA variable types
    // Returns priority value where lower numbers come first in ordering
    static int getCanonicalOrderPriority(const std::string& var);

    // Get sorted vector of variable pairs in canonical order
    static std::vector<std::pair<std::string, int>>
    getSortedVariablePairs(const std::map<std::string, int>& factors);

    // Legacy function - kept for backward compatibility but uses new canonical ordering
    static std::map<std::string, int>
    reorderCommutativeFactors(const std::map<std::string, int>& factors);

    // Detect and apply GA-specific simplifications
    // Example: (R.c0^2 + R.c1^2 - R.c2^2 - R.c3^2) patterns
    static SimplifiedTerm applyGASimplifications(const SimplifiedTerm& term);

    // Symmetry-based cancellation detection
    // Example: terms that cancel due to GA bivector antisymmetry
    static bool cancelsToZero(const SimplifiedTerm& term1, const SimplifiedTerm& term2);

  private:

    static bool isRotorCoefficient(const std::string& var);
    static bool isVectorComponent(const std::string& var);
    static bool isBivectorComponent(const std::string& var);
};

///////////////////////////////////////////////////////////////////////////////
// Test cases for validation
///////////////////////////////////////////////////////////////////////////////

class SimplificationTests {
  public:

    // Test EGA2D vector sandwich expansion
    // Input: "(R.c0 * v.x + R.c1 * v.y) * R.c0 + (R.c0 * v.y - R.c1 * v.x) * R.c1"
    // Expected: "R.c0*R.c0*v.x + R.c1*R.c0*v.y + R.c0*R.c1*v.y - R.c1*R.c1*v.x"
    static void testEGA2DExpansion();

    // Test commutativity application
    // Input: "R.c1*R.c0*v.y + R.c0*R.c1*v.y"
    // Expected: "2*R.c0*R.c1*v.y"
    static void testCommutativeCollection();

    // Test symmetry cancellation
    // Input: "R.c0*v.x*R.c1 - R.c1*v.x*R.c0"
    // Expected: "0" (should cancel due to commutativity)
    static void testSymmetryCancellation();

    // Test complete EGA2D transformation
    // Verify against manual example results
    static void testCompleteEGA2DTransformation();

    static void runAllTests();
};

///////////////////////////////////////////////////////////////////////////////
// Helper function declarations
///////////////////////////////////////////////////////////////////////////////

// Convert single simplified term back to AST representation
std::shared_ptr<ast_node> convertSingleTermToAst(const SimplifiedTerm& term);