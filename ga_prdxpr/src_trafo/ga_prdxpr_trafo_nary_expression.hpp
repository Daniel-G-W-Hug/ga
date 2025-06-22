#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// N-ary Expression Tree - Clean structure for manipulation
///////////////////////////////////////////////////////////////////////////////

// Represents a single multiplicative term: coefficient * var1 * var2 * ...
class NAryTerm {
  public:

    double coefficient;                   // e.g., -2.0
    std::map<std::string, int> variables; // e.g., {"R.c0": 1, "R.c1": 1, "v.x": 1}

    NAryTerm(double coeff = 1.0) : coefficient(coeff) {}

    // Operations
    NAryTerm operator*(const NAryTerm& other) const;
    bool canCombineWith(const NAryTerm& other) const; // Same variables
    NAryTerm operator+(const NAryTerm& other) const;  // Combine coefficients

    // Utilities
    std::string toString() const;
    bool isZero() const;
    void applyCommutativity(); // Sort variables canonically

  private:

    std::map<std::string, int> getCanonicalVariables() const;
};

// Represents an additive expression: term1 + term2 + term3 + ...
class NAryExpression {
  public:

    std::vector<NAryTerm> terms;

    NAryExpression() = default;
    NAryExpression(const std::vector<NAryTerm>& t) : terms(t) {}

    // Core 3-step transformation process (following manual approach)
    void expandDistributiveProducts(); // Step 1: Resolve braces by multiplication
    void
    normalizeSignsAndCommutativity(); // Step 2: Sign normalization and variable ordering
    void combineTermsAndRegroup();    // Step 3: Combine like terms
    void factorCommonVariables();     // Step 4: Factor out common variables for cleaner
                                      // parenthesization

    // Combined transformation
    void simplify()
    {
        expandDistributiveProducts();
        normalizeSignsAndCommutativity();
        combineTermsAndRegroup();
        factorCommonVariables();
        removeZeroTerms();
    }

    // Conversion and utilities
    std::string toString() const;
    void removeZeroTerms();
    size_t termCount() const { return terms.size(); }

    // Debug support
    void printDebug(const std::string& step_name) const;

  private:

    // Helper for factorization
    double getRotorCoefficientValue(const std::map<std::string, int>& rotor_vars);
};

///////////////////////////////////////////////////////////////////////////////
// Conversion functions: Binary AST ↔ N-ary Tree
///////////////////////////////////////////////////////////////////////////////

class NAryConverter {
  public:

    // Convert binary AST to N-ary expression tree (do this once)
    static NAryExpression fromBinaryAST(std::shared_ptr<ast_node> ast);

    // Convert N-ary expression tree back to binary AST (do this once)
    static std::shared_ptr<ast_node> toBinaryAST(const NAryExpression& expr);

  private:

    // Helper functions for AST → N-ary conversion
    static void extractTermsFromExpression(std::shared_ptr<Expression> expr,
                                           std::vector<NAryTerm>& terms,
                                           double sign = 1.0);
    static void extractTermsFromTerm(std::shared_ptr<Term> term,
                                     std::vector<NAryTerm>& terms, double sign = 1.0);
    static void extractFactorsFromTerm(std::shared_ptr<ast_node> node,
                                       NAryTerm& current_term, double sign = 1.0);

    // Helper functions for N-ary → AST conversion
    static std::shared_ptr<ast_node> buildTermAST(const NAryTerm& term);
    static std::shared_ptr<ast_node>
    buildExpressionAST(const std::vector<NAryTerm>& terms);
};

///////////////////////////////////////////////////////////////////////////////
// New Expression Simplifier using N-ary trees
///////////////////////////////////////////////////////////////////////////////

class NAryExpressionSimplifier {
  public:

    // Main entry point - replaces ExpressionSimplifier::simplify
    static std::shared_ptr<ast_node> simplify(std::shared_ptr<ast_node> ast);

    // For testing and debugging
    static NAryExpression simplifyToNAry(std::shared_ptr<ast_node> ast);
};