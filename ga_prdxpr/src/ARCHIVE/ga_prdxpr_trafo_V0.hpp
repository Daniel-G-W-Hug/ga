#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.
//

// transformation of simple product expressions

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Custom exceptions
class ParseError : public std::runtime_error {
  public:

    explicit ParseError(const std::string& msg) : std::runtime_error(msg) {}
};

class Term {
  public:

    double coefficient;
    std::map<std::string, int> variables; // variable name -> power

    Term();
    std::string toString() const;
    bool operator==(const Term& other) const;
    bool hasSameVariables(const Term& other) const;
};

class Expression {
  private:

    std::vector<Term> terms;
    static std::vector<std::string> tokenize(const std::string& expr);
    void multiplyTerms(Term& result, const Term& t1, const Term& t2);
    Expression parseSubExpression(const std::vector<std::string>& tokens, size_t& pos);
    void multiplyByCoefficient(double coeff);
    void addTerms(const std::vector<Term>& newTerms);

    static void validateParentheses(const std::string& expr);
    static void validateOperators(const std::string& expr);
    static void validateVariables(const std::string& expr);
    static bool isValidVariableName(const std::string& var);

    void combineLikeTerms();
    void removeZeroTerms();
    void sortTerms();

  public:

    void parse(const std::string& expr);
    void groupByVariable(const std::string& varPrefix);
    std::vector<Term> getTerms() const { return terms; }
    void simplify();
};
