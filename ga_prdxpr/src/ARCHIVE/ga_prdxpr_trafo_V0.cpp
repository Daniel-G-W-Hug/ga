// Copyright 2024-2025, Daniel Hug. All rights reserved.
//

#include "ga_prdxpr_trafo.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <regex>
#include <sstream>
#include <stack>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

Term::Term() : coefficient(1.0) {}

bool Term::operator==(const Term& other) const
{
    return std::abs(coefficient - other.coefficient) < 1e-10 &&
           variables == other.variables;
}

bool Term::hasSameVariables(const Term& other) const
{
    return variables == other.variables;
}

std::string Term::toString() const
{
    std::stringstream ss;
    if (coefficient != 1.0 || variables.empty()) {
        if (coefficient == -1.0 && !variables.empty()) {
            ss << "-";
        }
        else {
            ss << coefficient;
        }
    }
    bool first = true;
    for (const auto& var : variables) {
        if (!first || coefficient != 1.0) {
            ss << "*";
        }
        ss << var.first;
        // if (var.second > 1) {
        //     ss << "^" << var.second;
        // }
        if (var.second > 1) {
            ss << "*" << var.first;
        }
        first = false;
    }
    return ss.str();
}

void Expression::validateParentheses(const std::string& expr)
{
    std::stack<char> parentheses;
    for (char c : expr) {
        if (c == '(') {
            parentheses.push(c);
        }
        else if (c == ')') {
            if (parentheses.empty()) {
                throw ParseError("Unmatched closing parenthesis");
            }
            parentheses.pop();
        }
    }
    if (!parentheses.empty()) {
        throw ParseError("Unmatched opening parenthesis");
    }
}

void Expression::validateOperators(const std::string& expr)
{
    std::regex validOperators(R"([+\-*/])");
    std::regex invalidOperators(R"([/\\%^&|!@#$])");

    if (std::regex_search(expr, invalidOperators)) {
        throw ParseError("Invalid operators found. Only +, -, * are allowed");
    }

    std::regex consecutiveOperators(R"([+\-*/]\s*[+\-*/])");
    if (std::regex_search(expr, consecutiveOperators)) {
        throw ParseError("Invalid consecutive operators");
    }
}

bool Expression::isValidVariableName(const std::string& var)
{
    std::regex validPattern(R"(^(R\.c[0-3]|v\.[xyz])$)");
    return std::regex_match(var, validPattern);
}

void Expression::validateVariables(const std::string& expr)
{
    std::regex varPattern(R"(\b\w+\.\w+\b)");
    auto vars_begin = std::sregex_iterator(expr.begin(), expr.end(), varPattern);
    auto vars_end = std::sregex_iterator();

    for (std::sregex_iterator i = vars_begin; i != vars_end; ++i) {
        std::string var = i->str();
        if (!isValidVariableName(var)) {
            throw ParseError("Invalid variable name: " + var);
        }
    }
}

void Expression::combineLikeTerms()
{
    std::vector<Term> simplifiedTerms;

    for (const Term& term : terms) {
        bool combined = false;
        for (Term& existingTerm : simplifiedTerms) {
            if (term.hasSameVariables(existingTerm)) {
                existingTerm.coefficient += term.coefficient;
                combined = true;
                break;
            }
        }
        if (!combined) {
            simplifiedTerms.push_back(term);
        }
    }
    terms = simplifiedTerms;
}

void Expression::removeZeroTerms()
{
    terms.erase(std::remove_if(
                    terms.begin(), terms.end(),
                    [](const Term& term) { return std::abs(term.coefficient) < 1e-10; }),
                terms.end());
}

void Expression::sortTerms()
{
    std::sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
        if (a.variables.size() != b.variables.size()) {
            return a.variables.size() > b.variables.size();
        }
        if (!a.variables.empty() && !b.variables.empty()) {
            return a.variables.begin()->first < b.variables.begin()->first;
        }
        return std::abs(a.coefficient) > std::abs(b.coefficient);
    });
}

std::vector<std::string> Expression::tokenize(const std::string& expr)
{
    std::vector<std::string> tokens;
    std::regex pattern(R"([\+\-\*\(\)]|\w+(?:\.\w+)?)");

    auto words_begin = std::sregex_iterator(expr.begin(), expr.end(), pattern);
    auto words_end = std::sregex_iterator();

    for (std::sregex_iterator i = words_begin; i != words_end; ++i) {
        tokens.push_back(i->str());
    }
    return tokens;
}

void Expression::multiplyTerms(Term& result, const Term& t1, const Term& t2)
{
    result.coefficient = t1.coefficient * t2.coefficient;
    result.variables.clear();

    for (const auto& var : t1.variables) {
        result.variables[var.first] += var.second;
    }
    for (const auto& var : t2.variables) {
        result.variables[var.first] += var.second;
    }
}

void Expression::multiplyByCoefficient(double coeff)
{
    for (Term& term : terms) {
        term.coefficient *= coeff;
    }
}

void Expression::addTerms(const std::vector<Term>& newTerms)
{
    terms.insert(terms.end(), newTerms.begin(), newTerms.end());
}

Expression Expression::parseSubExpression(const std::vector<std::string>& tokens,
                                          size_t& pos)
{
    Expression result;
    std::vector<Term> currentTerms;
    Term currentTerm;
    bool expectOperator = false;
    bool unaryMinus = false;

    while (pos < tokens.size()) {
        const std::string& token = tokens[pos];

        // added for printing
        fmt::println("currentToken: '{}'\n", token);

        if (token == "(") {
            pos++;
            Expression subExpr = parseSubExpression(tokens, pos);

            // added for printing
            std::vector<Term> subTerms = subExpr.getTerms();
            for (auto const& st : subTerms) {
                fmt::println("subTerms: '{}'", st.toString());
            }

            size_t nextPos = pos + 1;
            while (nextPos < tokens.size() && tokens[nextPos] == "*") {
                nextPos++;
                if (nextPos < tokens.size()) {
                    Term coefTerm;
                    coefTerm.variables[tokens[nextPos]] = 1;

                    std::vector<Term> multipliedTerms;
                    for (const Term& term : subExpr.terms) {
                        Term resultTerm;
                        multiplyTerms(resultTerm, term, coefTerm);
                        multipliedTerms.push_back(resultTerm);
                    }
                    // added for printing
                    for (auto const& mt : multipliedTerms) {
                        fmt::println("        multipliedTerms: '{}'", mt.toString());
                    }
                    subExpr.terms = multipliedTerms;
                    pos = nextPos;
                }
                nextPos++;
            }

            if (unaryMinus) {
                for (Term& term : subExpr.terms) {
                    term.coefficient *= -1;
                }
                unaryMinus = false;
            }

            currentTerms.insert(currentTerms.end(), subExpr.terms.begin(),
                                subExpr.terms.end());
            expectOperator = true;

            // added for printing
            for (auto const& ct : currentTerms) {
                fmt::println("        currentTerms: '{}'", ct.toString());
            }
        }
        else if (token == ")") {
            if (!currentTerm.variables.empty() || currentTerm.coefficient != 1.0) {
                currentTerms.push_back(currentTerm);
            }
            terms = currentTerms;
            return *this;
        }
        else if (token == "*") {
            expectOperator = false;
        }
        else if (token == "+") {
            if (!currentTerm.variables.empty() || currentTerm.coefficient != 1.0) {
                currentTerms.push_back(currentTerm);
            }
            currentTerm = Term();
            expectOperator = false;
        }
        else if (token == "-") {
            if (!expectOperator) {
                unaryMinus = true;
            }
            else {
                if (!currentTerm.variables.empty() || currentTerm.coefficient != 1.0) {
                    currentTerms.push_back(currentTerm);
                }
                currentTerm = Term();
                currentTerm.coefficient = -1.0;
            }
            expectOperator = false;
        }
        else {
            if (unaryMinus) {
                currentTerm.coefficient *= -1;
                unaryMinus = false;
            }
            if (std::isdigit(token[0])) {
                currentTerm.coefficient *= std::stod(token);
            }
            else {
                currentTerm.variables[token]++;
            }
            expectOperator = true;
        }
        pos++;

        // added for printing
        fmt::println("        currentTerm: '{}'", currentTerm.toString());
    }

    // added for printing
    fmt::println("    currentTerm: '{}'", currentTerm.toString());

    if (!currentTerm.variables.empty() || currentTerm.coefficient != 1.0) {
        currentTerms.push_back(currentTerm);
    }
    terms = currentTerms;

    // added for printing
    for (auto const& ct : currentTerms) {
        fmt::println("    currentTerms: '{}'", ct.toString());
    }

    return *this;
}

void Expression::parse(const std::string& expr)
{
    try {
        validateParentheses(expr);
        validateOperators(expr);
        validateVariables(expr);

        std::vector<std::string> tokens = tokenize(expr);

        fmt::println("");
        fmt::println("tokens: '{}'", fmt::join(tokens, "' '"));
        fmt::println("");

        size_t pos = 0;
        *this = parseSubExpression(tokens, pos);

        // simplify();
    }
    catch (const ParseError& e) {
        std::cerr << "Parse error: " << e.what() << std::endl;
        throw;
    }
    catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
        throw;
    }
}

void Expression::simplify()
{
    combineLikeTerms();
    removeZeroTerms();
    sortTerms();
}

void Expression::groupByVariable(const std::string& varPrefix)
{
    std::map<std::string, std::vector<Term>> groupedTerms;

    for (const Term& term : terms) {
        std::string key;
        for (const auto& var : term.variables) {
            if (var.first.substr(0, varPrefix.length()) == varPrefix) {
                key = var.first;
                break;
            }
        }
        if (!key.empty()) {
            groupedTerms[key].push_back(term);
        }
    }

    for (const auto& group : groupedTerms) {
        std::cout << group.first << " * (";
        bool first = true;
        for (const Term& term : group.second) {
            if (!first) {
                if (term.coefficient >= 0) {
                    std::cout << " + ";
                }
                else {
                    std::cout << " ";
                }
            }
            else if (term.coefficient < 0) {
                std::cout << "-";
            }
            first = false;
            Term simplified = term;
            simplified.coefficient = std::abs(term.coefficient);
            simplified.variables.erase(group.first);
            std::cout << simplified.toString();
        }
        std::cout << ")" << std::endl;
    }
}
