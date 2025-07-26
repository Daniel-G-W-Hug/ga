// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo.hpp"

#include <algorithm>
#include <cctype>
#include <iterator>
#include <sstream>
#include <stack>
#include <stdexcept>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

///////////////////////////////////////////////////////////////////////////////
// Validation classes
///////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////
// Lexer implementation
///////////////////////////////////////////////////////////////////////////////

Lexer::Lexer(std::string const& input_in) : input(input_in), position(0)
{
    ExpressionValidator::validateParentheses(input);
    ExpressionValidator::validateOperatorSequence(input);
    advance(); // make currentToken available

    auto token = getCurrentToken();
    // fmt::println("Lexer::Lexer()  : token.str_value = {}", token.str_value);
    // fmt::println("Lexer::Lexer()  : token.type = {}", Token_t_toString(token.type));
    // fmt::println("");
}

char Lexer::peek() const
{
    if (position >= input.length()) return '\0';
    return input[position];
}

char Lexer::advance_char()
{
    if (position >= input.length()) return '\0';
    return input[position++];
}

void Lexer::skipWhitespace()
{
    while (position < input.length() && std::isspace(input[position])) {
        position++;
    }
}

bool Lexer::isDigit(char c) const { return std::isdigit(c) != 0; }

bool Lexer::isLetter(char c) const { return std::isalpha(c) != 0 || c == '_'; }

Token Lexer::parseNumber()
{
    std::string number;
    size_t start_pos = position;

    // Handle sign
    if (peek() == '+' || peek() == '-') {
        number += advance_char();
    }

    // Parse integer part
    while (isDigit(peek())) {
        number += advance_char();
    }

    // Parse decimal part (beginning with '.')
    if (peek() == '.') {
        number += advance_char();
        while (isDigit(peek())) {
            number += advance_char();
        }
    }

    // Handle scientific notation
    if (peek() == 'd' || peek() == 'D') {
        number += advance_char();
    }

    if (number.empty() ||
        (number.length() == 1 && (number[0] == '+' || number[0] == '-'))) {
        throw std::runtime_error("Invalid number at position " +
                                 std::to_string(start_pos));
    }

    return Token(Token_t::NUMBER, number, start_pos);
}

Token Lexer::parseIdentifier()
{
    std::string identifier;
    size_t start_pos = position;

    // First character must be a letter or underscore
    if (!isLetter(peek())) {
        throw std::runtime_error("Invalid identifier at position " +
                                 std::to_string(position));
    }

    while (isLetter(peek()) || isDigit(peek())) {
        identifier += advance_char();
    }

    return Token(Token_t::IDENTIFIER, identifier, start_pos);
}

Token Lexer::getNextToken()
{
    skipWhitespace();

    if (position >= input.length()) {
        return Token(Token_t::END, "", position);
    }

    char current = peek();

    if (isDigit(current)) {
        return parseNumber();
    }

    if (isLetter(current)) {
        return parseIdentifier();
    }

    size_t currentPos = position;
    advance_char(); // consume character
    switch (current) {
        case '.':
            return Token(Token_t::DOT, ".", currentPos);
        case '+':
            return Token(Token_t::PLUS, "+", currentPos);
        case '-':
            return Token(Token_t::MINUS, "-", currentPos);
        case '*':
            return Token(Token_t::MULTIPLY, "*", currentPos);
        case '/':
            return Token(Token_t::DEVIDE, "/", currentPos);
        case '(':
            return Token(Token_t::LPAREN, "(", currentPos);
        case ')':
            return Token(Token_t::RPAREN, ")", currentPos);
        default:
            throw std::runtime_error("Invalid character at position " +
                                     std::to_string(currentPos));
    }
}


///////////////////////////////////////////////////////////////////////////////
// Parser implementation
///////////////////////////////////////////////////////////////////////////////

Parser::Parser(const std::string& input) : lexer(input) {}

std::shared_ptr<ast_node> Parser::parse()
{

    auto result = Expression::parse(lexer);

    auto currentToken = lexer.getCurrentToken();
    if (currentToken.type != Token_t::END) {
        throw std::runtime_error("Unexpected tokens after expression at position " +
                                 std::to_string(currentToken.position));
    }

    return result;
}

void Parser::validateBinaryOperation(Token const& op,
                                     std::shared_ptr<ast_node> const& left,
                                     std::shared_ptr<ast_node> const& right)
{
    if (!left || !right) {
        throw std::runtime_error("Missing operand for operator '" + op.str_value +
                                 "' at position " + std::to_string(op.position));
    }
}

// Expression ctor
Expression::Expression(std::shared_ptr<ast_node> l, char o, std::shared_ptr<ast_node> r) :
    left(l), op(o), right(r)
{
}

// Expression parsing
std::shared_ptr<ast_node> Expression::parse(Lexer& lexer)
{
    auto left = Term::parse(lexer);

    while (true) {
        auto token = lexer.getCurrentToken();

        // fmt::println("Expression::parse() 2nd: token.str_value = {}", token.str_value);
        // fmt::println("Expression::parse() 2nd: token.type = {}",
        //              Token_t_toString(token.type));
        // fmt::println("");

        if (token.type != Token_t::PLUS && token.type != Token_t::MINUS) {
            break;
        }

        char op = token.str_value[0];
        lexer.advance(); // consume '+' or '-'

        auto right = Term::parse(lexer);
        Parser::validateBinaryOperation(token, left, right);

        left = std::make_shared<Expression>(left, op, right);
    }

    return left;
}

// Term ctor
Term::Term(std::shared_ptr<ast_node> l, char o, std::shared_ptr<ast_node> r) :
    left(l), op(o), right(r)
{
}

// Term parsing
std::shared_ptr<ast_node> Term::parse(Lexer& lexer)
{
    auto left = Factor::parse(lexer);

    while (true) {
        auto token = lexer.getCurrentToken();

        // fmt::println("Term::parse() 2nd: token.str_value = {}", token.str_value);
        // fmt::println("Term::parse() 2nd: token.type = {}",
        // Token_t_toString(token.type)); fmt::println("");

        if (token.type != Token_t::MULTIPLY && token.type != Token_t::DEVIDE) {
            break;
        }

        auto op = token.str_value[0];
        lexer.advance(); // consume '*' or '/'

        auto right = Factor::parse(lexer);
        Parser::validateBinaryOperation(token, left, right);
        left = std::make_shared<Term>(left, op, right);
    }

    return left;
}

// Factor ctor
Factor::Factor(char s, std::shared_ptr<ast_node> v) : sign(s), prim_val(v) {}


// Factor parsing
std::shared_ptr<ast_node> Factor::parse(Lexer& lexer)
{

    auto token = lexer.getCurrentToken();
    // fmt::println("Factor::parse() : token.str_value = {}", token.str_value);
    // fmt::println("Factor::parse() : token.type = {}", Token_t_toString(token.type));
    // fmt::println("");

    // work on the sign, if present
    char sign = '\0'; // '\0' indicates no sign
    if (token.type == Token_t::PLUS || token.type == Token_t::MINUS) {
        sign = token.str_value[0];
        lexer.advance(); // consume sign
    }

    // read the value
    auto primary = Primary::parse(lexer);

    return make_shared<Factor>(sign, primary);
}


// Primary parsing
std::shared_ptr<ast_node> Primary::parse(Lexer& lexer)
{
    auto primary = std::make_shared<Primary>();

    auto token = lexer.getCurrentToken();
    // fmt::println("Primary::parse(): token.str_value = {}", token.str_value);
    // fmt::println("Primary::parse(): token.type = {}", Token_t_toString(token.type));
    // fmt::println("");

    switch (token.type) {

        case Token_t::NUMBER: {
            primary->type = Primary_t::NUMBER;
            primary->str_value = token.str_value;
            primary->num_value = std::stod(token.str_value);
            lexer.advance(); // consume the number
            break;
        }

        case Token_t::IDENTIFIER: {
            primary->type = Primary_t::VARIABLE;
            primary->str_value = token.str_value;
            lexer.advance(); // consume the variable

            // Check for component access (e.g., "v.x")
            if (lexer.getCurrentToken().type == Token_t::DOT) {
                lexer.advance();
                if (lexer.getCurrentToken().type != Token_t::IDENTIFIER) {
                    throw std::runtime_error(
                        "Expected identifier after dot at position " +
                        std::to_string(lexer.getCurrentToken().position));
                }
                primary->str_value += "." + lexer.getCurrentToken().str_value;
                lexer.advance(); // consume the component
            }
            break;
        }

        case Token_t::LPAREN: {

            primary->type = Primary_t::EXPRESSION;
            lexer.advance(); // consume the opening parenthesis

            primary->expr = Expression::parse(lexer);
            primary->str_value = primary->expr->to_string();

            // check for closing parenthesis
            if (lexer.getCurrentToken().type != Token_t::RPAREN) {
                throw std::runtime_error(
                    "Expected closing parenthesis at position " +
                    std::to_string(lexer.getCurrentToken().position));
            }
            lexer.advance(); // consume the closing parenthesis
            break;
        }

        default:
            throw std::runtime_error("Unexpected token at position " +
                                     std::to_string(token.position));
    }

    return primary;
}

ast_node::ast_node_t Expression::nodeType() const
{
    return ast_node::ast_node_t::EXPRESSION;
}
ast_node::ast_node_t Term::nodeType() const { return ast_node::ast_node_t::TERM; }
ast_node::ast_node_t Factor::nodeType() const { return ast_node::ast_node_t::FACTOR; }
ast_node::ast_node_t Primary::nodeType() const { return ast_node::ast_node_t::PRIMARY; }

std::string Expression::nodeType_to_string() const { return "EXPRESSION"s; }
std::string Term::nodeType_to_string() const { return "TERM"s; }
std::string Factor::nodeType_to_string() const { return "FACTOR"s; }
std::string Primary::nodeType_to_string() const { return "PRIMARY"s; }


std::string Primary::primaryType_to_string() const
{
    switch (type) {
        case Primary_t::NUMBER:
            return "NUMBER";
            break;
        case Primary_t::VARIABLE:
            return "VARIABLE";
            break;
        case Primary_t::EXPRESSION:
            return "EXPRESSION";
            break;
        default:
            return "UNKNOWN";
    }
};

///////////////////////////////////////////////////////////////////////////////
// Helper function implementation
///////////////////////////////////////////////////////////////////////////////

std::string parse_and_print_ast(std::string const& input)
{
    try {
        Parser parser(input);
        auto ast = parser.parse();

        // fmt::println("");
        // fmt::println("ast: {}", ast->to_string());
        // fmt::println("");

        {
            fmt::println("");
            fmt::println("ast (pre_order traversal):\n");
            int node_cnt = -1; // start with a non-existing node number -1
            print_parse_tree_pre_order(ast, 0, node_cnt); // start at lvl 0
            fmt::println("");
        }

        {
            fmt::println("");
            fmt::println("ast (level_order traversal):\n");
            std::queue<lvl_node> node_queue;
            lvl_node root{0, ast};
            node_queue.push(root); // push root node into queue
            int node_cnt = -1;     // start with a non-existing node number -1
            print_parse_tree_level_order(node_queue, node_cnt);
            fmt::println("");
        }


        return ast->to_string();
    }
    catch (std::exception const& e) {
        return std::string("Error: ") + e.what();
    }
}

std::string parse_only(std::string const& input)
{
    try {
        Parser parser(input);
        auto ast = parser.parse();

        return ast->to_string();
    }
    catch (std::exception const& e) {
        return std::string("Error: ") + e.what();
    }
}

///////////////////////////////////////////////////////////////////////////////
// Printing functions
///////////////////////////////////////////////////////////////////////////////

void print_parse_tree_pre_order(std::shared_ptr<ast_node> const& ast, size_t lvl,
                                int& node_cnt)
{

    if (!ast) return;

    ++node_cnt;

    // fmt::println("print_parse_tree(): ast node type: {}",
    // ast->nodeType_to_string());

    switch (ast->nodeType()) {

        case ast_node::ast_node_t::EXPRESSION:

        {
            // downcast to derived
            auto ptr = std::dynamic_pointer_cast<Expression>(ast);
            print_expression_node(ptr, lvl, node_cnt);

            // upcast back to base after getting the left and right sides
            if (ptr->left) {
                print_parse_tree_pre_order(std::dynamic_pointer_cast<ast_node>(ptr->left),
                                           lvl + 1, node_cnt);
            }
            if (ptr->right) {
                print_parse_tree_pre_order(
                    std::dynamic_pointer_cast<ast_node>(ptr->right), lvl + 1, node_cnt);
            }

        } break;

        case ast_node::ast_node_t::TERM:

        {
            // downcast to derived
            auto ptr = std::dynamic_pointer_cast<Term>(ast);
            print_term_node(ptr, lvl, node_cnt);

            // upcast back to base after getting the left and right sides
            if (ptr->left) {
                print_parse_tree_pre_order(std::dynamic_pointer_cast<ast_node>(ptr->left),
                                           lvl + 1, node_cnt);
            }
            if (ptr->right) {
                print_parse_tree_pre_order(
                    std::dynamic_pointer_cast<ast_node>(ptr->right), lvl + 1, node_cnt);
            }

        } break;

        case ast_node::ast_node_t::FACTOR:

        {
            auto ptr = std::dynamic_pointer_cast<Factor>(ast);
            if (ptr) {
                print_factor_node(ptr, lvl, node_cnt);
            }
        } break;


        case ast_node::ast_node_t::PRIMARY:

        {
            auto ptr = std::dynamic_pointer_cast<Primary>(ast);
            if (ptr) {
                print_primary_node(ptr, lvl, node_cnt);
            }
        } break;

        default:
            fmt::println("printParseTree: Unknown type.");
    }

    // fmt::println("print_parse_tree_pre_order(): printing finished.");
}

void print_parse_tree_level_order(std::queue<lvl_node>& node_queue, int& node_cnt)
{

    while (!node_queue.empty()) {

        auto node = node_queue.front(); // get the next node in the queue
        node_queue.pop();               // and remove the node from the queue

        if (!node.ptr) break;

        auto lvl = node.lvl;
        auto ast = node.ptr;

        ++node_cnt;

        switch (ast->nodeType()) {

            case ast_node::ast_node_t::EXPRESSION:

            {
                // downcast to derived
                auto ptr = std::dynamic_pointer_cast<Expression>(ast);
                print_expression_node(ptr, lvl, node_cnt);

                if (ptr->left) {

                    lvl_node left{node.lvl + 1, ptr->left};
                    node_queue.push(left);
                }
                if (ptr->right) {

                    lvl_node right{node.lvl + 1, ptr->right};
                    node_queue.push(right);
                }

            } break;

            case ast_node::ast_node_t::TERM:

            {
                // downcast to derived
                auto ptr = std::dynamic_pointer_cast<Term>(ast);
                print_term_node(ptr, lvl, node_cnt);

                if (ptr->left) {

                    lvl_node left{node.lvl + 1, ptr->left};
                    node_queue.push(left);
                }
                if (ptr->right) {

                    lvl_node right{node.lvl + 1, ptr->right};
                    node_queue.push(right);
                }

            } break;

            case ast_node::ast_node_t::FACTOR:

            {
                auto ptr = std::dynamic_pointer_cast<Factor>(ast);
                if (ptr) {
                    print_factor_node_level_order(node_queue, ptr, lvl, node_cnt);
                }
            } break;


            case ast_node::ast_node_t::PRIMARY:

            {
                auto ptr = std::dynamic_pointer_cast<Primary>(ast);
                if (ptr) {
                    print_primary_node_level_order(node_queue, ptr, lvl, node_cnt);
                }
            } break;

            default:
                fmt::println("printParseTree: Unknown type.");
        }
    }

    // fmt::println("print_parse_tree_level_order(): printing finished.");

    return;
}

void print_expression_node(std::shared_ptr<Expression> const& ptr, size_t lvl,
                           int& node_cnt)
{

    if (!ptr) return;

    fmt::println("lvl {:3}, node# {:3}, {}, op: '{}', string: {}", lvl, node_cnt,
                 ptr->nodeType_to_string(), ptr->op, ptr->to_string());
    if (ptr->left) {
        fmt::println("         left  {} -> {}", ptr->left->nodeType_to_string(),
                     ptr->left->to_string());
    }
    if (ptr->right) {
        fmt::println("         right {} -> {}", ptr->right->nodeType_to_string(),
                     ptr->right->to_string());
    }
    fmt::println("");
}

void print_term_node(std::shared_ptr<Term> const& ptr, size_t lvl, int& node_cnt)
{

    if (!ptr) return;

    fmt::println("lvl {:3}, node# {:3}, {}, op: '{}', string: {}", lvl, node_cnt,
                 ptr->nodeType_to_string(), ptr->op, ptr->to_string());
    if (ptr->left) {
        fmt::println("         left  {} -> {}", ptr->left->nodeType_to_string(),
                     ptr->left->to_string());
    }
    if (ptr->right) {
        fmt::println("         right {} -> {}", ptr->right->nodeType_to_string(),
                     ptr->right->to_string());
    }
    fmt::println("");

    return;
}

void print_factor_node(std::shared_ptr<Factor> const& ptr, size_t lvl, int& node_cnt)
{

    if (!ptr) return;

    fmt::println("lvl {:3}, node# {:3}, {}, sign: '{}'", lvl, node_cnt,
                 ptr->nodeType_to_string(), ptr->sign);

    if (ptr->prim_val) {
        auto primary_ptr = std::dynamic_pointer_cast<Primary>(ptr->prim_val);
        print_primary_node(primary_ptr, lvl, node_cnt);
    }

    return;
}

void print_primary_node(std::shared_ptr<Primary> const& ptr, size_t lvl, int& node_cnt)
{

    if (!ptr) return;

    ++node_cnt;

    fmt::println("         node# {:3}, {}_{} -> {}", node_cnt, ptr->nodeType_to_string(),
                 ptr->primaryType_to_string(), ptr->str_value);
    if (ptr->primaryType_to_string() == "NUMBER") {
        fmt::println("         num_val -> {}", ptr->num_value);
    }
    fmt::println("");

    if (ptr->expr) {
        auto expr_ptr = std::dynamic_pointer_cast<Expression>(ptr->expr);
        print_parse_tree_pre_order(std::dynamic_pointer_cast<ast_node>(expr_ptr), lvl + 1,
                                   node_cnt);
    }

    return;
}

void print_factor_node_level_order(std::queue<lvl_node>& node_queue,
                                   std::shared_ptr<Factor> const& ptr, size_t lvl,
                                   int& node_cnt)
{

    if (!ptr) return;

    fmt::println("lvl {:3}, node# {:3}, {}, sign: '{}'", lvl, node_cnt,
                 ptr->nodeType_to_string(), ptr->sign);

    if (ptr->prim_val) {
        auto primary_ptr = std::dynamic_pointer_cast<Primary>(ptr->prim_val);
        print_primary_node_level_order(node_queue, primary_ptr, lvl, node_cnt);
    }

    return;
}

void print_primary_node_level_order(std::queue<lvl_node>& node_queue,
                                    std::shared_ptr<Primary> const& ptr, size_t lvl,
                                    int& node_cnt)
{

    if (!ptr) return;

    ++node_cnt;

    fmt::println("         node# {:3}, {}_{} -> {}", node_cnt, ptr->nodeType_to_string(),
                 ptr->primaryType_to_string(), ptr->str_value);
    if (ptr->primaryType_to_string() == "NUMBER") {
        fmt::println("         num_val -> {}", ptr->num_value);
    }
    fmt::println("");

    if (ptr->expr) {
        lvl_node expr_node{lvl + 1, ptr->expr};
        node_queue.push(expr_node);
    }

    return;
}


///////////////////////////////////////////////////////////////////////////////
// xxx::to_string() implementations
///////////////////////////////////////////////////////////////////////////////


// ToString implementations for AST nodes
std::string Expression::to_string() const
{
    std::stringstream ss;
    if (left) {
        ss << left->to_string();
    }
    if (right) {
        ss << " " << op << " " << right->to_string();
    }
    return ss.str();
}

std::string Term::to_string() const
{
    std::stringstream ss;
    if (left) {
        ss << left->to_string();
    }
    if (right) {
        ss << " " << op << " " << right->to_string();
    }
    return ss.str();
}

std::string Factor::to_string() const
{
    std::stringstream ss;
    if (sign != '\0') {
        ss << sign;
    }
    if (prim_val) {
        ss << prim_val->to_string();
    }
    return ss.str();
}

std::string Primary::to_string() const
{
    switch (type) {
        case Primary_t::NUMBER:
        case Primary_t::VARIABLE:
            return str_value;
        case Primary_t::EXPRESSION:
            if (expr) {
                return "(" + expr->to_string() + ")";
            }
            else {
                return "EXPRESSION EMPTY";
            }
        default:
            // return "";
            return "EXPRESSION EMPTY";
    }
}
