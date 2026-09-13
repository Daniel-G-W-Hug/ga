#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "sandwich/ga_prdxpr_sandwich_parser.hpp"

#include <map>
#include <memory>
#include <string>
#include <utility>
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
    SimplifiedTerm(double coeff, std::map<std::string, int> const& vars, int s = 1);

    // Create from string representation
    static SimplifiedTerm fromString(std::string const& term_str);

    // Operations
    SimplifiedTerm operator*(SimplifiedTerm const& other) const;
    SimplifiedTerm operator+(SimplifiedTerm const& other) const; // Only if same variables
    bool canCombineWith(SimplifiedTerm const& other) const;
    bool isZero() const;
    bool isEquivalentTo(SimplifiedTerm const& other) const;

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
    static std::shared_ptr<ast_node> termsToAst(std::vector<SimplifiedTerm> const& terms);

    // Test functions for specific GA patterns
    static void testEGA2DExpansion();
    static void testSymmetryCancellation();
    static void testCommutativeReordering();

    // Public for testing
    static std::vector<SimplifiedTerm>
    combineTermsVector(std::vector<SimplifiedTerm> const& terms);

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
    static SimplifiedTerm expandSingleTerm(SimplifiedTerm const& term);
    static bool
    isCommutativeVariable(std::string const& var); // R.c0, R.c1, etc. are commutative
    static SimplifiedTerm applyCommutativeReordering(SimplifiedTerm const& term);
};

///////////////////////////////////////////////////////////////////////////////
// Geometric Variable Pattern Configuration
///////////////////////////////////////////////////////////////////////////////

// Configuration structure for all geometric algebra variable patterns
struct GeometricVariablePatterns {
    std::string coeff_prefix = "R.c"; // Coefficient prefix: "R.c" or "M.c"
    std::vector<std::string> vectors = {"v.x", "v.y", "v.z", "v.w"}; // Vector components
    std::vector<std::string> bivectors = {"B.x", "B.y", "B.z"}; // Bivector components
    std::vector<std::string> trivectors = {"t.x", "t.y", "t.z",
                                           "t.w"}; // Trivector components
    std::vector<std::string> quadvectors = {};     // Quadvector components (5D algebras)

    // Predefined pattern factories
    static GeometricVariablePatterns createEGA3DPatterns();
    static GeometricVariablePatterns createPGA3DPPatterns();
    static GeometricVariablePatterns createEGA2DPatterns();
    static GeometricVariablePatterns createPGA2DPPatterns();
    static GeometricVariablePatterns createSTA4DSPatterns();
    static GeometricVariablePatterns createCGA2DCPatterns();
    static GeometricVariablePatterns createCGA3DCPatterns();
};

///////////////////////////////////////////////////////////////////////////////
// Specific GA algebraic rules
///////////////////////////////////////////////////////////////////////////////

class GAAlgebraRules {
  public:

    // Comprehensive canonical ordering for all GA variable types
    // Returns priority value where lower numbers come first in ordering
    static int getCanonicalOrderPriority(std::string const& var,
                                         GeometricVariablePatterns const& patterns = {});

    // Get sorted vector of variable pairs in canonical order
    static std::vector<std::pair<std::string, int>>
    getSortedVariablePairs(std::map<std::string, int> const& factors,
                           GeometricVariablePatterns const& patterns = {});

    // Legacy overload for simple coefficient prefix (backward compatibility)
    static std::vector<std::pair<std::string, int>>
    getSortedVariablePairs(std::map<std::string, int> const& factors,
                           std::string const& coeff_prefix);

    // Legacy function - kept for backward compatibility but uses new canonical ordering
    static std::map<std::string, int>
    reorderCommutativeFactors(std::map<std::string, int> const& factors,
                              GeometricVariablePatterns const& patterns = {});

    // Legacy overload for simple coefficient prefix (backward compatibility)
    static std::map<std::string, int>
    reorderCommutativeFactors(std::map<std::string, int> const& factors,
                              std::string const& coeff_prefix);

    // Detect and apply GA-specific simplifications
    // Example: (R.c0^2 + R.c1^2 - R.c2^2 - R.c3^2) patterns
    static SimplifiedTerm applyGASimplifications(SimplifiedTerm const& term);

    // Symmetry-based cancellation detection
    // Example: terms that cancel due to GA bivector antisymmetry
    static bool cancelsToZero(SimplifiedTerm const& term1, SimplifiedTerm const& term2);

  private:

    // Private helper functions removed - were unused legacy code
};

// Test classes have been moved to ga_prdxpr_trafo_tests.hpp for better organization

///////////////////////////////////////////////////////////////////////////////
// Helper function declarations
///////////////////////////////////////////////////////////////////////////////

// Convert single simplified term back to AST representation
std::shared_ptr<ast_node> convertSingleTermToAst(SimplifiedTerm const& term);
