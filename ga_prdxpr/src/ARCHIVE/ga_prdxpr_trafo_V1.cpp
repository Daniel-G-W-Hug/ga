// author: Daniel Hug, 2025
//

#include "ga_prdxpr_trafo.hpp"

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

#include <algorithm>
#include <cctype>
#include <iterator>
#include <sstream>
#include <stack>
#include <stdexcept>

// maximum number of recursive call levels for the parser
// and maximum number of operands in expressions
constexpr inline size_t max_lvl = 32;

// stack to transfer sign information between parsing expression, term and primary on
// different levels
//
// opening a brace with minus in front needs to be put onto and removed from the stack in
// parsePrimary(); new level also starts with new terms in expressions and reading new
// primaries for the right side of new terms
//
// parseExpression() and parseTerm() use the information to push minus signs to their
// operands as long as minus_in_flight indicates that they should push a resulting minus
// sign to their operands. (which is true if the stack size is odd)
// push value on stack that shows the level of the negative brace starting with one
static std::vector<size_t> mobr_on_lvl(max_lvl, 0); // minus before opening brace on lvl
static std::vector<size_t> mcvr_on_lvl(max_lvl, 0); // minus carry over from lower lvl

// TODO: refactor from global static object to at least static sign tracking object
// bundling the tracking fields and the tracking routines

void reset_sign_tracking_vectors()
{
    // reset the static fields for sign tracking to zero,
    // relevant in case of consecutive calls to parse()
    for (size_t i = 0; i < max_lvl; ++i) {
        mobr_on_lvl[i] = 0;
        mcvr_on_lvl[i] = 0;
    }
}

void print_sign_tracking_vectors()
{
    fmt::println("mobr_on_lvl  : {}", fmt::join(mobr_on_lvl, ", "));
    fmt::println("mcvr_on_lvl  : {}", fmt::join(mcvr_on_lvl, ", "));
    fmt::print("'-' in flight: ");
    for (size_t i = 0; i < max_lvl; ++i) {
        fmt::print("{:d}, ", minus_in_flight_on_level(i));
    }
}

void mark_is_negative_on_level(size_t lvl)
{
    mobr_on_lvl[lvl] = 1; // indicates '-' in front of opening brace
    mcvr_on_lvl[lvl] = 0; // no need for carry over from previous braces
}

void check_for_carry_over_on_level(size_t lvl)
{
    // requires a call to mark_is_negative_on_level prior to check for carry over
    if ((mobr_on_lvl[lvl - 1] && !mobr_on_lvl[lvl]) ||
        (mcvr_on_lvl[lvl - 1] && !mobr_on_lvl[lvl])) {
        // we have a carry over from previous opening braces
        mcvr_on_lvl[lvl] = 1;
    }
}

bool minus_in_flight_on_level(size_t lvl)
{

    if (lvl >= max_lvl) {
        throw std::runtime_error("Level too large. Increase max_lvl. lvl = " +
                                 std::to_string(lvl));
    }

    size_t cnt{0};
    for (size_t i = 0; i <= lvl; ++i) {
        if (mobr_on_lvl[i] == 1) cnt++;
    }
    auto number_of_minus_signs_is_odd = ((cnt % 2) == 1);
    // fmt::println("l={}, cnt={}, is_odd={:d}", lvl, cnt, number_of_minus_signs_is_odd);


    if ((mobr_on_lvl[lvl] || mcvr_on_lvl[lvl]) && number_of_minus_signs_is_odd) {
        // odd number of signs given and either sign directly or carry over
        return true;
    }
    else {
        return false;
    }
}

// for counting components in expressions
// (start with 1 in each call of expression on the respective level)
static std::vector<size_t> expr_cnt_components_on_level(max_lvl, 0);


// ExpressionValidator implementation
bool ExpressionValidator::validateParentheses(std::string const& expr)
{
    std::stack<char> parentheses;
    for (size_t i = 0; i < expr.length(); i++) {
        char c = expr[i];
        if (c == '(') {
            parentheses.push(c);
        }
        else if (c == ')') {
            if (parentheses.empty()) {
                throw std::runtime_error("Unmatched closing parenthesis at position " +
                                         std::to_string(i));
            }
            parentheses.pop();
        }
    }
    if (!parentheses.empty()) {
        throw std::runtime_error("Unmatched opening parenthesis found");
    }
    return true;
}

bool ExpressionValidator::validateVariableName(std::string const& name)
{
    if (name.empty()) return false;

    if (!isalpha(name[0]) && name[0] != '_') {
        throw std::runtime_error("Invalid variable name '" + name +
                                 "': must start with letter or underscore");
    }

    size_t dotPos = name.find('.');
    if (dotPos != std::string::npos) {
        std::string baseName = name.substr(0, dotPos);
        if (!validateIdentifier(baseName)) {
            throw std::runtime_error("Invalid base variable name '" + baseName + "'");
        }

        if (dotPos + 1 >= name.length()) {
            throw std::runtime_error("Missing component after dot in '" + name + "'");
        }
        std::string componentName = name.substr(dotPos + 1);
        if (!validateIdentifier(componentName)) {
            throw std::runtime_error("Invalid component name '" + componentName + "'");
        }
    }
    else {
        if (!validateIdentifier(name)) {
            throw std::runtime_error("Invalid variable name '" + name + "'");
        }
    }

    return true;
}

bool ExpressionValidator::validateOperatorSequence(std::string const& expr)
{
    char lastChar = '\0';
    bool wasOperator = false;

    for (size_t i = 0; i < expr.length(); i++) {
        char c = expr[i];

        if (isspace(c)) continue;

        bool isOperator = (c == '+' || c == '-' || c == '*');

        if (isOperator && wasOperator) {
            throw std::runtime_error("Consecutive operators found near position " +
                                     std::to_string(i));
        }

        if (isOperator && i == expr.length() - 1) {
            throw std::runtime_error("Expression cannot end with an operator");
        }

        if (isOperator && (lastChar == '\0' || lastChar == '(') && c != '-') {
            throw std::runtime_error("Invalid unary operator '" + std::string(1, c) +
                                     "' at position " + std::to_string(i));
        }

        lastChar = c;
        wasOperator = isOperator;
    }

    return true;
}

bool ExpressionValidator::validateIdentifier(std::string const& name)
{
    if (name.empty()) return false;

    if (!isalpha(name[0]) && name[0] != '_') return false;

    for (size_t i = 1; i < name.length(); i++) {
        if (!isalnum(name[i]) && name[i] != '_') return false;
    }

    return true;
}

// Token implementation
Token::Token(Token_t t, std::string const& str_val, size_t pos) :
    type(t), str_value(str_val), position(pos)
{
}

// Lexer implementation
Lexer::Lexer(std::string const& input) : input(input), position(0)
{
    ExpressionValidator::validateParentheses(input);
    ExpressionValidator::validateOperatorSequence(input);
}

char Lexer::peek() const
{
    if (position >= input.length()) return '\0';
    return input[position];
}

char Lexer::advance()
{
    if (position >= input.length()) return '\0';
    return input[position++];
}

void Lexer::skipWhitespace()
{
    while (isspace(peek()))
        advance();
}

Token Lexer::readNumber()
{
    size_t startPos = position;
    std::string num;
    bool hasDecimal = false;

    if (peek() == '-' || peek() == '+') {
        num += advance();
    }

    while (isdigit(peek()) || peek() == '.' || peek() == 'd' || peek() == 'D') {
        char c = advance();
        if (c == '.') {
            if (hasDecimal) {
                throw std::runtime_error(
                    "Invalid number format: multiple decimal points at position " +
                    std::to_string(position));
            }
            hasDecimal = true;
        }
        num += c;
    }

    if (num == "." || num == "-." || num == "+.") {
        throw std::runtime_error("Invalid number format at position " +
                                 std::to_string(startPos));
    }

    return Token(Token_t::NUMBER, num, startPos);
}

Token Lexer::readIdentifier()
{
    size_t startPos = position;
    std::string id;
    while (isalnum(peek()) || peek() == '_' || peek() == '.') {
        id += advance();
    }

    try {
        ExpressionValidator::validateVariableName(id);
    }
    catch (const std::runtime_error& e) {
        throw std::runtime_error(std::string(e.what()) + " at position " +
                                 std::to_string(startPos));
    }

    return Token(Token_t::IDENTIFIER, id, startPos);
}

Token Lexer::nextToken()
{
    skipWhitespace();

    char c = peek();
    if (c == '\0') return Token(Token_t::END, "", position);

    if (isdigit(c)) {
        return readNumber();
    }

    if (isalpha(c) || c == '_') {
        return readIdentifier();
    }

    size_t currentPos = position;
    advance(); // consume the character
    switch (c) {
        case '+':
            return Token(Token_t::PLUS, "+", currentPos);
        case '-':
            return Token(Token_t::MINUS, "-", currentPos);
        case '*':
            return Token(Token_t::MULTIPLY, "*", currentPos);
        case '(':
            return Token(Token_t::LPAREN, "(", currentPos);
        case ')':
            return Token(Token_t::RPAREN, ")", currentPos);
        case '.':
            return Token(Token_t::DOT, ".", currentPos);
        default:
            throw std::runtime_error("Invalid character '" + std::string(1, c) +
                                     "' at position " + std::to_string(currentPos));
    }
}

// // debug version of Lexer::nextToken()
// Token Lexer::nextToken()
// {
//     skipWhitespace();

//     char c = peek();
//     if (c == '\0') {
//         std::println("DEBUG: Lexer returning END token at position {}", position);
//         return Token(Token_t::END, "", position);
//     }

//     std::println("DEBUG: Lexer c = {}", c);

//     Token token(Token_t::END);

//     if (isdigit(c)) {
//         token = readNumber();
//     }
//     else if (isalpha(c) || c == '_') {
//         token = readIdentifier();
//     }
//     else {
//         size_t currentPos = position;
//         advance(); // consume the character
//         switch (c) {
//             case '+':
//                 token = Token(Token_t::PLUS, "+", currentPos);
//                 break;
//             case '-':
//                 token = Token(Token_t::MINUS, "-", currentPos);
//                 break;
//             case '*':
//                 token = Token(Token_t::MULTIPLY, "*", currentPos);
//                 break;
//             case '(':
//                 token = Token(Token_t::LPAREN, "(", currentPos);
//                 break;
//             case ')':
//                 token = Token(Token_t::RPAREN, ")", currentPos);
//                 break;
//             case '.':
//                 token = Token(Token_t::DOT, ".", currentPos);
//                 break;
//             default:
//                 throw std::runtime_error("Invalid character '" + std::string(1, c)
//                 +
//                                          "' at position " +
//                                          std::to_string(currentPos));
//         }
//     }
//     // Debug output for each token
//     fmt::println("DEBUG: Token_t = {}, Value = '{}', Position = {}",
//     Token_t_toString(token.type), token.str_value, token.position);
//
//     return token;
// }


// Expression implementation
Expression::Expression(std::shared_ptr<ExprNode> l, char o, std::shared_ptr<ExprNode> r) :
    left(l), op(o), right(r)
{
}

std::string Expression::toString(bool use_braces) const
{

    std::string prefix{};
    std::string postfix{};
    std::string r_prefix{};
    std::string r_postfix{};

    if (str_value_starts_with_minus) {
        prefix = "-("s;
        postfix = ")"s;
    }
    else if (use_braces) {
        prefix = "("s;
        postfix = ")"s;
    }

    if (op == '-' && right->nodeType() == ExprNode_t::EXPRESSION) {
        r_prefix = "("s;
        r_postfix = ")"s;
    }

    return prefix + left->toString() + " " + op + " " + r_prefix + right->toString() +
           r_postfix + postfix;

    /////

    // if (str_value_starts_with_minus) {
    //     return "-(" + left->toString() + " " + op + " " + right->toString() + ")";
    // }

    // if (op == '+') {
    //     return left->toString() + " " + op + " " + right->toString();
    // }
    // else if (op == '-' && right->nodeType() == ExprNode_t::EXPRESSION) {
    //     return left->toString() + " " + op + " (" + right->toString() + ")";
    // }
    // else {
    //     return left->toString() + " " + op + " " + right->toString();
    // }
}

ExprNode_t Expression::nodeType() const { return ExprNode_t::EXPRESSION; }

std::shared_ptr<ExprNode> Expression::simplify()
{
    auto simplifiedLeft = left->simplify();
    auto simplifiedRight = right->simplify();
    return std::make_shared<Expression>(simplifiedLeft, op, simplifiedRight);
}

// Term implementation
Term::Term(std::shared_ptr<ExprNode> l, std::shared_ptr<ExprNode> r) : left(l), right(r)
{
}

std::string Term::toString(bool use_braces) const
{

    std::string prefix{};
    std::string postfix{};
    std::string l_prefix{};
    std::string l_postfix{};
    std::string r_prefix{};
    std::string r_postfix{};

    if (str_value_starts_with_minus) {
        prefix = "-"s;
        postfix = ""s;
    }
    else if (use_braces) {
        prefix = "("s;
        postfix = ")"s;
    }

    if (left->nodeType() == ExprNode_t::EXPRESSION) {
        l_prefix = "("s;
        l_postfix = ")"s;
    }
    if (right->nodeType() == ExprNode_t::EXPRESSION) {
        r_prefix = "("s;
        r_postfix = ")"s;
    }

    return prefix + l_prefix + left->toString() + l_postfix + " * " + r_prefix +
           right->toString() + r_postfix + postfix;
}

ExprNode_t Term::nodeType() const { return ExprNode_t::TERM; }

std::shared_ptr<ExprNode> Term::simplify()
{
    auto simplifiedLeft = left->simplify();
    auto simplifiedRight = right->simplify();
    return std::make_shared<Term>(simplifiedLeft, simplifiedRight);
}

// Primary implementation
Primary::Primary(std::string const& str_val, Primary_t t) :
    type(t), str_value(str_val), str_value_starts_with_minus(str_val[0] == '-')
{
    if (t == Primary_t::NUMBER) {
        num_value = std::stod(str_value);
    }
}


std::string Primary::toString(bool use_braces) const
{
    switch (type) {

        case Primary_t::NUMBER:
        case Primary_t::VARIABLE:

            return str_value;
            break;

        case Primary_t::EXPRESSION:

            return "("s + str_value + ")"s;
            break;

        default:
            return "ERROR"s;
    }
}

ExprNode_t Primary::nodeType() const { return ExprNode_t::PRIMARY; }

std::shared_ptr<ExprNode> Primary::simplify()
{
    // if (type == Primary_t::EXPRESSION) {
    //     auto simplifiedExpr = expr->simplify();
    //     return std::make_shared<Primary>(simplifiedExpr,
    //     str_value_starts_with_minus);
    // }
    return std::make_shared<Primary>(*this);
}

void toggle_sign_of_primary(std::shared_ptr<Primary> ptr)
{

    if (ptr->type == Primary_t::NUMBER) {
        if (ptr->str_value_starts_with_minus) {
            // remove the minus
            ptr->str_value = ptr->str_value.substr(1, ptr->str_value.size());
            ptr->num_value *= -1.0;
            ptr->str_value_starts_with_minus = false;
        }
        else {
            // add a minus
            ptr->str_value = '-' + ptr->str_value;
            ptr->num_value *= -1.0;
            ptr->str_value_starts_with_minus = true;
        }
    }
    if (ptr->type == Primary_t::VARIABLE) {
        // if (ptr->str_value_sarts_with_minus) {
        if (ptr->str_value_starts_with_minus) {
            // remove the minus
            ptr->str_value = ptr->str_value.substr(1, ptr->str_value.size());
            ptr->str_value_starts_with_minus = false;
        }
        else {
            // add a minus
            ptr->str_value = '-' + ptr->str_value;
            ptr->str_value_starts_with_minus = true;
        }
    }
}

void resolve_signs_in_expression_with_primaries(std::shared_ptr<Expression> expr_ptr)
{

    if (expr_ptr->left->nodeType() == ExprNode_t::PRIMARY) {

        auto ptr = std::dynamic_pointer_cast<Primary>(expr_ptr->left);
        toggle_sign_of_primary(ptr);
    }

    // toggle the operator
    if (expr_ptr->op == '-') expr_ptr->op = '+';
    if (expr_ptr->op == '+') expr_ptr->op = '-';

    // toggle the sign in connection with the binary operator
    if (expr_ptr->right->nodeType() == ExprNode_t::PRIMARY) {

        auto ptr = std::dynamic_pointer_cast<Primary>(expr_ptr->right);

        if (ptr->str_value_starts_with_minus) {
            toggle_sign_of_primary(ptr);
        }
    }
}

void resolve_signs_in_expression_with_term(std::shared_ptr<Term> term_ptr)
{

    if (term_ptr->left->nodeType() == ExprNode_t::PRIMARY) {

        auto ptr = std::dynamic_pointer_cast<Primary>(term_ptr->left);

        if (ptr->str_value_starts_with_minus) {
            toggle_sign_of_primary(ptr);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Parser implementation
///////////////////////////////////////////////////////////////////////////////

Parser::Parser(std::string const& input) : lexer(input), currentToken() { advance(); }

void Parser::advance() { currentToken = lexer.nextToken(); }

// // debug version of Parser::advance()
// void Parser::advance()
// {
//     currentToken = lexer.nextToken();
//     fmt::println(
//         "DEBUG: Parser advanced to token: {} at position {}",
//         Token_t_toString(currentToken.type), currentToken.position);
// }

void Parser::validateBinaryOperation(Token const& op,
                                     std::shared_ptr<ExprNode> const& left,
                                     std::shared_ptr<ExprNode> const& right)
{
    if (!left || !right) {
        throw std::runtime_error("Missing operand for operator '" + op.str_value +
                                 "' at position " + std::to_string(op.position));
    }
}

std::shared_ptr<ExprNode> Parser::parsePrimary(size_t lvl)
{

    // fmt::println("Parser::parsePrimary():    lvl = {0}, "
    //              "minus_in_flight_on_level({0}) = {1}, currentToken.str_value =
    //              '{2}', " "currentToken type = {3}", lvl,
    //              minus_in_flight_on_level(lvl), currentToken.str_value,
    //              Token_t_toString(currentToken.type));

    // local analysis only (minus found directly in front of number or variable)
    bool is_negative_in_primary = false;

    if (currentToken.type == Token_t::MINUS) {

        // fmt::println("Parser::parsePrimary():    local minus");
        is_negative_in_primary = true; // now we know that we have a minus sign
        advance();                     // consume '-'
    }

    if (currentToken.type == Token_t::NUMBER) {

        // fmt::println("Parser::parsePrimary():    read NUMBER");
        std::string str_value = currentToken.str_value;

        if ((is_negative_in_primary && !minus_in_flight_on_level(lvl)) ||
            (!is_negative_in_primary && minus_in_flight_on_level(lvl))) {
            str_value = '-' + currentToken.str_value;
        }
        advance(); // consume the parsed number

        auto ptr = std::make_shared<Primary>(str_value, Primary_t::NUMBER);
        return std::dynamic_pointer_cast<ExprNode>(ptr);
    }

    if (currentToken.type == Token_t::IDENTIFIER) {
        // fmt::println("Parser::parsePrimary():    read VARIABLE");
        std::string str_value = currentToken.str_value;

        if ((is_negative_in_primary && !minus_in_flight_on_level(lvl)) ||
            (!is_negative_in_primary && minus_in_flight_on_level(lvl))) {
            str_value = '-' + currentToken.str_value;
        }
        advance(); // consume the parsed variable

        auto ptr = std::make_shared<Primary>(str_value, Primary_t::VARIABLE);
        return std::dynamic_pointer_cast<ExprNode>(ptr);
    }

    if (currentToken.type == Token_t::LPAREN) {

        // tracking of minus signs in front of opening parenthesis
        if (is_negative_in_primary) {

            // now we have a negative sign which needs to be handled in upcoming calls
            // of parseExpression(), parseTerm() and parsePrimary() accordingly
            mark_is_negative_on_level(lvl + 1);

            // for debug printing only
            // fmt::println("Parser::parsePrimary():    opening brace found for lvl {}",
            //              lvl + 1);
        }
        check_for_carry_over_on_level(lvl + 1);
    }

    size_t openParenPos = currentToken.position;
    advance(); // consume '('

    auto expr = parseExpression(lvl + 1);

    // fmt::println("currentToken.type = {}", currentToken.value);
    if (currentToken.type != Token_t::RPAREN) {
        throw std::runtime_error("Expected ')' after expression starting at position " +
                                 std::to_string(openParenPos));
    }

    advance(); // consume ')'

    return std::dynamic_pointer_cast<ExprNode>(expr);


    // if (currentToken.type == Token_t::MINUS) {
    //     last_call_was_minus = true; // mark that to assign the right sign to
    //     numbers Token minus = currentToken; advance(); auto primary =
    //     parsePrimary(); if (!primary) {
    //         throw std::runtime_error("Missing operand after unary '-' at position "
    //         +
    //                                  std::to_string(minus.position));
    //     }
    //     return std::make_shared<Primary>(primary);
    // }

    throw std::runtime_error("Unexpected token '" + currentToken.str_value +
                             "' at position " + std::to_string(currentToken.position));
}

std::shared_ptr<ExprNode> Parser::parseTerm(size_t lvl)
{
    // fmt::println("Parser::parseTerm():       before calling parsePrimary()");

    auto left = parsePrimary(lvl);

    // fmt::println("Parser::parseTerm():       after calling parsePrimary()");
    // fmt::println("Parser::parseTerm():       lvl = {0}, "
    //              "minus_in_flight_on_level({0}) = {1}, "
    //              " currentToken.str_value = '{2}', currentToken.type = {3}",
    //              lvl, minus_in_flight_on_level(lvl), currentToken.str_value,
    //              Token_t_toString(currentToken.type));

    // fmt::println("Parser::parseTerm():       should potentially remove stuff
    // here");

    bool first_term = true;
    while (currentToken.type == Token_t::MULTIPLY) {

        // fmt::println("Parser::parseTerm():       should remove stuff here");
        // fmt::println("Parser::parseTerm():       left->nodeType() = {}",
        //              ExprNode_t_toString(left->nodeType()));

        if (first_term && minus_in_flight_on_level(lvl) &&
            expr_cnt_components_on_level[lvl] > 1) {

            // push minus sign into left operand of term, if there is a minus in
            // flight and we are at least at the second component of the expression
            auto ptr = std::dynamic_pointer_cast<Primary>(left);

            toggle_sign_of_primary(ptr);
        }

        first_term = false;

        Token op = currentToken;
        advance(); // consume '*'

        // fmt::println(
        //     "Parser::parseTerm():       parsing next primary calling
        //     parsePrimary(lvl+1)");

        auto right = parsePrimary(lvl + 1);
        validateBinaryOperation(op, left, right);
        left = std::make_shared<Term>(left, right);
    }

    return left;
}

std::shared_ptr<ExprNode> Parser::parseExpression(size_t lvl)
{

    // fmt::println("Parser::parseExpression(): before calling parseTerm()");

    if (lvl >= max_lvl) {
        throw std::runtime_error("Level too large. Increase max_lvl. Level = " +
                                 std::to_string(lvl));
    }
    expr_cnt_components_on_level[lvl] = 1;

    auto left = parseTerm(lvl);

    // fmt::println("Parser::parseExpression(): left->nodeType() = {}",
    //              ExprNode_t_toString(left->nodeType()));

    // fmt::println("Parser::parseExpression(): after calling parsePrimary()");
    // fmt::println("Parser::parseExpression(): lvl = {0}, "
    //              "minus_in_flight_on_level({0}) = {1}, "
    //              " currentToken.str_value = '{2}', currentToken.type = {3}",
    //              lvl, minus_in_flight_on_level(lvl), currentToken.str_value,
    //              Token_t_toString(currentToken.type));

    // fmt::println("lvl: {0}, minus_in_flight_on_level\({0}\): {1}, "
    //              "opening_brace_with_minus: {2}",
    //              lvl, minus_in_flight_on_level(lvl),
    //              fmt::join(opening_brace_with_minus, ", "));

    // resolve expressions with primaries on the lhs
    // if (minus_in_flight_on_level(lvl) && left->nodeType() == ExprNode_t::EXPRESSION) {

    //     // push minus sign into left operand of current expression, if there is a
    //     // minus in flight
    //     auto expr_ptr = std::dynamic_pointer_cast<Expression>(left);

    //     resolve_signs_in_expression_with_primaries(expr_ptr);
    // }

    while (currentToken.type == Token_t::PLUS || currentToken.type == Token_t::MINUS) {

        // now we get to the next components, thus increase the counter
        // to be used in parseTerm() for distinguishing between the first and other
        // terms
        expr_cnt_components_on_level[lvl] += 1;

        Token op = currentToken;
        advance(); // consume '+' or '-'
        auto right = parseTerm(lvl);
        validateBinaryOperation(op, left, right);

        char expression_operator = op.type == Token_t::PLUS ? '+' : '-';

        if (minus_in_flight_on_level(lvl)) {
            expression_operator = op.type == Token_t::PLUS ? '-' : '+';
        }

        // toggle the sign in connection with the binary operator
        if (minus_in_flight_on_level(lvl) && right->nodeType() == ExprNode_t::PRIMARY) {

            auto ptr = std::dynamic_pointer_cast<Primary>(right);

            if (ptr->str_value_starts_with_minus) {
                toggle_sign_of_primary(ptr);
            }
        }

        // toggle the sign in connection with the binary operator
        if (minus_in_flight_on_level(lvl) && right->nodeType() == ExprNode_t::TERM) {

            auto term_ptr = std::dynamic_pointer_cast<Term>(right);
            resolve_signs_in_expression_with_term(term_ptr);
        }

        // // resolve signs in expressions with primaries on the rhs
        // if (minus_in_flight_on_level(lvl) &&
        //     right->nodeType() == ExprNode_t::EXPRESSION) {
        //     auto expr_ptr = std::dynamic_pointer_cast<Expression>(right);
        //     resolve_signs_in_expression_with_primaries(expr_ptr);
        // }

        left = std::make_shared<Expression>(left, expression_operator, right);
    }

    return left;
}

std::shared_ptr<ExprNode> Parser::parse()
{
    // reset the static fields for sign tracking to zero,
    // essential in case of more than one call to parse()
    reset_sign_tracking_vectors();

    auto result = parseExpression(0); // start on level 0
    if (currentToken.type != Token_t::END) {
        throw std::runtime_error("Unexpected tokens after expression at position " +
                                 std::to_string(currentToken.position));
    }
    return result;
}

// Helper function implementation
std::string parse_and_print_ast(std::string const& input)
{
    try {
        Parser parser(input);
        auto ast = parser.parse();

        // fmt::println("");
        // fmt::println("ast: {}", ast->toString());
        // fmt::println("");

        fmt::println("");
        fmt::println("ast:");
        print_parse_tree(ast);
        fmt::println("");

        auto simplified = ast->simplify();
        return simplified->toString();
    }
    catch (std::exception const& e) {
        return std::string("Error: ") + e.what();
    }
}

// Helper function implementation
std::string parse_no_print_ast(std::string const& input)
{
    try {
        Parser parser(input);
        auto ast = parser.parse();

        auto simplified = ast->simplify();
        return simplified->toString();
    }
    catch (std::exception const& e) {
        return std::string("Error: ") + e.what();
    }
}

void print_parse_tree(std::shared_ptr<ExprNode> const& ast)
{

    if (ast.get() == nullptr) return;

    ExprNode_t type = ast->nodeType();

    switch (type) {

        case ExprNode_t::EXPRESSION:

        {
            // downcast to derived
            auto ptr = std::dynamic_pointer_cast<Expression>(ast);
            print_expression_node(ptr);

            // upcast back to base after getting the left and right sides
            print_parse_tree(std::dynamic_pointer_cast<ExprNode>(ptr->left));
            print_parse_tree(std::dynamic_pointer_cast<ExprNode>(ptr->right));
        }

        break;

        case ExprNode_t::TERM:

        {
            // downcast to derived
            auto ptr = std::dynamic_pointer_cast<Term>(ast);
            print_term_node(ptr);

            // upcast back to base after getting the left and right sides
            print_parse_tree(std::dynamic_pointer_cast<ExprNode>(ptr->left));
            print_parse_tree(std::dynamic_pointer_cast<ExprNode>(ptr->right));
        }

        break;

        case ExprNode_t::PRIMARY:

            print_primary_node(std::dynamic_pointer_cast<Primary>(ast));
            break;

        default:
            fmt::println("printParseTree: Unknown type.");
    }
}

void print_expression_node(std::shared_ptr<Expression> const& ptr)
{

    if (!ptr) return;

    fmt::println("node type                     : {}",
                 ExprNode_t_toString(ptr->nodeType()));
    fmt::println("    operation symbol          : {}", ptr->op);
    fmt::println("    expression to_String()    : {}", ptr->toString());
    fmt::println("    left  node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->left->nodeType()), ptr->left->toString());
    fmt::println("    right node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->right->nodeType()), ptr->right->toString());
    fmt::println("    expression starts with '-': {}", ptr->str_value_starts_with_minus);
    fmt::println("");
}

void print_term_node(std::shared_ptr<Term> const& ptr)
{

    if (!ptr) return;

    fmt::println("node type                     : {}",
                 ExprNode_t_toString(ptr->nodeType()));
    fmt::println("    term to_String()          : {}", ptr->toString());
    fmt::println("    left  node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->left->nodeType()), ptr->left->toString());
    fmt::println("    right node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->right->nodeType()), ptr->right->toString());
    fmt::println("    term starts with '-'      : {}", ptr->str_value_starts_with_minus);
    fmt::println("");

    return;
}

void print_primary_node(std::shared_ptr<Primary> const& ptr)
{

    if (!ptr) return;

    fmt::println("node type                     : {}",
                 ExprNode_t_toString(ptr->nodeType()));
    fmt::println("    primary to_String()       : {}", ptr->toString());
    fmt::println("    primary type              : {}", Primary_t_toString(ptr->type));
    fmt::println("    primary str_value         : {}", ptr->str_value);
    fmt::println("    primary num_value         : {}", ptr->num_value);
    fmt::println("    primary starts with '-'   : {}", ptr->str_value_starts_with_minus);
    fmt::println("");

    return;
}