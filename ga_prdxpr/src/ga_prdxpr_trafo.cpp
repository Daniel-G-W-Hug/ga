// author: Daniel Hug, 2025
//

#include "ga_prdxpr_trafo.hpp"

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

#include <cctype>
#include <sstream>
#include <stdexcept>

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
//                 throw std::runtime_error("Invalid character '" + std::string(1, c) +
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

std::string Expression::toString(bt_t bt) const
{
    if (str_value_starts_with_minus) {
        return "-(" + left->toString(bt_t::no_braces) + " " + op + " " +
               right->toString(bt_t::no_braces) + ")";
    }

    switch (bt) {
        case bt_t::use_braces:
            return "(" + left->toString(bt_t::no_braces) + " " + op + " " +
                   right->toString(bt_t::no_braces) + ")";
            break;
        case bt_t::no_braces:
            return left->toString(bt_t::no_braces) + " " + op + " " +
                   right->toString(bt_t::no_braces);
            break;
        default:
            return "UNKNOWN";
    }
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

std::string Term::toString(bt_t bt) const
{
    switch (bt) {
        case bt_t::use_braces:
            return "(" + left->toString(bt_t::no_braces) + ") * (" +
                   right->toString(bt_t::no_braces) + ")";
            break;
        case bt_t::no_braces:
            return left->toString(bt_t::no_braces) + " * " +
                   right->toString(bt_t::no_braces);
            break;
        default:
            return "UNKNOWN";
    }
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

// Primary::Primary(std::shared_ptr<ExprNode> e, bool neg) :
//     type(Primary_t::EXPRESSION), expr(e), str_value_starts_with_minus(neg)
// {
// }

std::string Primary::toString(bt_t bt) const
{
    std::string unary_minus = str_value_starts_with_minus ? "-" : "";
    switch (type) {

        case Primary_t::NUMBER:
        case Primary_t::VARIABLE:

            switch (bt) {

                case bt_t::use_braces:
                    // return "(" + unary_minus + str_value + ")";
                    return "(" + str_value + ")";
                    break;
                case bt_t::no_braces:
                    // return unary_minus + str_value;
                    return str_value;
                    break;
                default:
                    return "UNKNOWN";
            }
        case Primary_t::EXPRESSION:

            switch (bt) {
                case bt_t::use_braces:
                    return unary_minus + "(" + toString(bt_t::use_braces) + ")";
                    break;
                case bt_t::no_braces:
                    return unary_minus + toString(bt_t::no_braces);
                    break;
                default:
                    return "UNKNOWN";
            }

        default:
            return "ERROR";
    }
}

ExprNode_t Primary::nodeType() const { return ExprNode_t::PRIMARY; }

std::shared_ptr<ExprNode> Primary::simplify()
{
    // if (type == Primary_t::EXPRESSION) {
    //     auto simplifiedExpr = expr->simplify();
    //     return std::make_shared<Primary>(simplifiedExpr, str_value_starts_with_minus);
    // }
    return std::make_shared<Primary>(*this);
}

// Parser implementation
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

std::shared_ptr<ExprNode> Parser::parsePrimary()
{

    bool is_negative = false;

    // fmt::println(
    //     "Parser::parsePrimary():    currentToken.str_value = '{}', currentToken type "
    //     "= {}",
    //     currentToken.str_value, Token_t_toString(currentToken.type));

    if (currentToken.type == Token_t::MINUS) {
        is_negative = true; // now we know that we have a minus sign in front
        advance();          // consume '-'
    }

    if (currentToken.type == Token_t::NUMBER) {
        std::string str_value = currentToken.str_value;
        if (is_negative) {
            str_value = '-' + currentToken.str_value;
        }
        advance(); // consume the parsed number

        auto ptr = std::make_shared<Primary>(str_value, Primary_t::NUMBER);
        return std::dynamic_pointer_cast<ExprNode>(ptr);
    }

    if (currentToken.type == Token_t::IDENTIFIER) {
        std::string str_value = currentToken.str_value;
        if (is_negative) {
            str_value = '-' + currentToken.str_value;
        }
        advance(); // consume the parsed variable

        auto ptr = std::make_shared<Primary>(str_value, Primary_t::VARIABLE);
        return std::dynamic_pointer_cast<ExprNode>(ptr);
    }

    if (currentToken.type == Token_t::LPAREN) {
        size_t openParenPos = currentToken.position;
        advance(); // consume '('

        auto expr = parseExpression();

        // fmt::println("currentToken.type = {}", currentToken.value);
        if (currentToken.type != Token_t::RPAREN) {
            throw std::runtime_error(
                "Expected ')' after expression starting at position " +
                std::to_string(openParenPos));
        }

        advance(); // consume ')'


        if (is_negative) {

            auto expr_ptr = std::dynamic_pointer_cast<Expression>(expr);

            // push unary minus to the operands for primaries

            // lhs change via pushing the sign into the primary
            if (expr_ptr->left->nodeType() == ExprNode_t::PRIMARY) {

                // work directly on the primary on the left
                auto ptr = std::dynamic_pointer_cast<Primary>(expr_ptr->left);

                // unary minus has to change sign
                if ((ptr->type == Primary_t::VARIABLE ||
                     ptr->type == Primary_t::NUMBER) &&
                    !ptr->str_value_starts_with_minus) {

                    // add minus sign
                    ptr->str_value = '-' + ptr->str_value;
                    if (ptr->type == Primary_t::NUMBER) ptr->num_value *= -1.0;
                    ptr->str_value_starts_with_minus = true;
                }
                else if ((ptr->type == Primary_t::VARIABLE ||
                          ptr->type == Primary_t::NUMBER) &&
                         ptr->str_value_starts_with_minus) {

                    // remove minus sign
                    ptr->str_value = ptr->str_value.substr(1, ptr->str_value.size());
                    if (ptr->type == Primary_t::NUMBER) ptr->num_value *= -1.0;
                    ptr->str_value_starts_with_minus = false;
                }
            }

            // lhs change via pushing the sign into the expression
            if (expr_ptr->left->nodeType() == ExprNode_t::EXPRESSION) {

                auto ptr = std::dynamic_pointer_cast<Expression>(expr_ptr->left);

                // push unary minus to the operands for primaries

                // lhs change via pushing the sign into the primary
                if (ptr->left->nodeType() == ExprNode_t::PRIMARY) {

                    // work directly on the primary on the left
                    auto e_ptr = std::dynamic_pointer_cast<Primary>(ptr->left);

                    // unary minus has to change sign
                    if ((e_ptr->type == Primary_t::VARIABLE ||
                         e_ptr->type == Primary_t::NUMBER) &&
                        !e_ptr->str_value_starts_with_minus) {

                        // add minus sign
                        e_ptr->str_value = '-' + e_ptr->str_value;
                        if (e_ptr->type == Primary_t::NUMBER) e_ptr->num_value *= -1.0;
                        e_ptr->str_value_starts_with_minus = true;
                    }
                    else if ((e_ptr->type == Primary_t::VARIABLE ||
                              e_ptr->type == Primary_t::NUMBER) &&
                             e_ptr->str_value_starts_with_minus) {

                        // remove minus sign
                        e_ptr->str_value =
                            e_ptr->str_value.substr(1, e_ptr->str_value.size());
                        if (e_ptr->type == Primary_t::NUMBER) e_ptr->num_value *= -1.0;
                        e_ptr->str_value_starts_with_minus = false;
                    }
                }

                // rhs change via operator '+' or '-' (there is no unary minus on the rhs)
                if (ptr->right->nodeType() == ExprNode_t::PRIMARY && ptr->op == '+') {
                    ptr->op = '-';
                }
                else if (ptr->right->nodeType() == ExprNode_t::PRIMARY &&
                         ptr->op == '-') {
                    ptr->op = '+';
                }
            }

            // rhs change via operator '+' or '-' (there is no unary minus on the rhs)
            if (expr_ptr->right->nodeType() == ExprNode_t::PRIMARY &&
                expr_ptr->op == '+') {
                expr_ptr->op = '-';
            }
            else if (expr_ptr->right->nodeType() == ExprNode_t::PRIMARY &&
                     expr_ptr->op == '-') {
                expr_ptr->op = '+';
            }
        }

        return std::dynamic_pointer_cast<ExprNode>(expr);
    }


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

std::shared_ptr<ExprNode> Parser::parseTerm()
{
    auto left = parsePrimary();

    // fmt::println("Parser::parseTerm():       currentToken.str_value = '{}', "
    //              "currentToken.type = {}",
    //              currentToken.str_value, Token_t_toString(currentToken.type));

    while (currentToken.type == Token_t::MULTIPLY) {
        Token op = currentToken;
        advance(); // consume '*'
        auto right = parsePrimary();
        validateBinaryOperation(op, left, right);
        left = std::make_shared<Term>(left, right);
    }

    return left;
}

std::shared_ptr<ExprNode> Parser::parseExpression()
{
    auto left = parseTerm();

    // fmt::println("Parser::parseExpression(): currentToken.str_value = '{}', "
    //              "currentToken.type = {}",
    //              currentToken.str_value, Token_t_toString(currentToken.type));

    while (currentToken.type == Token_t::PLUS || currentToken.type == Token_t::MINUS) {
        Token op = currentToken;
        advance(); // consume '+' or '-'
        auto right = parseTerm();
        validateBinaryOperation(op, left, right);
        left = std::make_shared<Expression>(left, op.type == Token_t::PLUS ? '+' : '-',
                                            right);
    }

    return left;
}

std::shared_ptr<ExprNode> Parser::parse()
{
    auto result = parseExpression();
    if (currentToken.type != Token_t::END) {
        throw std::runtime_error("Unexpected tokens after expression at position " +
                                 std::to_string(currentToken.position));
    }
    return result;
}

// Helper function implementation
std::string parseAndSimplify(std::string const& input)
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

    if (ptr.get() == nullptr) return;

    fmt::println("node type                     : {}",
                 ExprNode_t_toString(ptr->nodeType()));
    fmt::println("    operation symbol          : {}", ptr->op);
    // fmt::println("    expression to_String()    : {}",
    // ptr->toString(bt_t::no_braces));
    fmt::println("    expression to_String()    : {}", ptr->toString());
    fmt::println("    left  node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->left->nodeType()),
                 ptr->left->toString(bt_t::no_braces));
    fmt::println("    right node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->right->nodeType()),
                 ptr->right->toString(bt_t::no_braces));
    fmt::println("    expression starts with '-': {}", ptr->str_value_starts_with_minus);
    fmt::println("");
}

void print_term_node(std::shared_ptr<Term> const& ptr)
{

    if (ptr.get() == nullptr) return;

    fmt::println("node type                     : {}",
                 ExprNode_t_toString(ptr->nodeType()));
    fmt::println("    term to_String()          : {}", ptr->toString());
    fmt::println("    left  node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->left->nodeType()),
                 ptr->left->toString(bt_t::no_braces));
    fmt::println("    right node type -> value  : {} -> {}",
                 ExprNode_t_toString(ptr->right->nodeType()),
                 ptr->right->toString(bt_t::no_braces));
    fmt::println("    term starts with '-'      : {}", ptr->str_value_starts_with_minus);
    fmt::println("");

    return;
}

void print_primary_node(std::shared_ptr<Primary> const& ptr)
{

    if (ptr.get() == nullptr) return;

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