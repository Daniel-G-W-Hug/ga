#include "prdxpr_parser.hpp"
#include <cctype>
#include <sstream>

// Lexer implementation
Lexer::Lexer(const std::string& input) :
    input(input), position(0), currentToken(TokenType::END, "", 0)
{
    advance();
}

char Lexer::peek() const
{
    if (position >= input.length()) {
        return '\0';
    }
    return input[position];
}

char Lexer::advance_char()
{
    if (position >= input.length()) {
        return '\0';
    }
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
    bool hasDot = false;

    // Handle sign
    if (peek() == '+' || peek() == '-') {
        number += advance_char();
    }

    // Parse integer part
    while (isDigit(peek())) {
        number += advance_char();
    }

    // Parse decimal part
    if (peek() == '.') {
        hasDot = true;
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
        throw ParserError("Invalid number at position " + std::to_string(start_pos));
    }

    return Token(TokenType::NUMBER, number, start_pos);
}

Token Lexer::parseIdentifier()
{
    std::string identifier;
    size_t start_pos = position;

    // First character must be a letter or underscore
    if (!isLetter(peek())) {
        throw ParserError("Invalid identifier at position " + std::to_string(position));
    }

    while (isLetter(peek()) || isDigit(peek())) {
        identifier += advance_char();
    }

    return Token(TokenType::IDENTIFIER, identifier, start_pos);
}

Token Lexer::getNextToken()
{
    skipWhitespace();

    if (position >= input.length()) {
        return Token(TokenType::END, "", position);
    }

    char current = peek();

    if (isDigit(current) || current == '+' || current == '-' || current == '.') {
        return parseNumber();
    }

    if (isLetter(current)) {
        return parseIdentifier();
    }

    position++;
    switch (current) {
        case '.':
            return Token(TokenType::DOT, ".", position - 1);
        case '+':
            return Token(TokenType::PLUS, "+", position - 1);
        case '-':
            return Token(TokenType::MINUS, "-", position - 1);
        case '*':
            return Token(TokenType::MULTIPLY, "*", position - 1);
        case '/':
            return Token(TokenType::DIVIDE, "/", position - 1);
        case '(':
            return Token(TokenType::LPAREN, "(", position - 1);
        case ')':
            return Token(TokenType::RPAREN, ")", position - 1);
        default:
            throw ParserError("Invalid character at position " +
                              std::to_string(position - 1));
    }
}

// Parser implementation
Parser::Parser(const std::string& input) : lexer(input) {}

std::unique_ptr<Expression> Parser::parse() { return Expression::parse(lexer); }

// Expression parsing
std::unique_ptr<Expression> Expression::parse(Lexer& lexer)
{
    auto expr = std::make_unique<Expression>();
    expr->left = Term::parse(lexer);

    while (true) {
        auto token = lexer.getCurrentToken();
        if (token.type != TokenType::PLUS && token.type != TokenType::MINUS) {
            break;
        }

        expr->op = token.value[0];
        lexer.advance();

        auto right = Term::parse(lexer);
        auto new_expr = std::make_unique<Expression>();
        new_expr->left = std::move(expr->left);
        new_expr->op = expr->op;
        new_expr->right = std::move(expr);
        expr = std::move(new_expr);
    }

    return expr;
}

// Term parsing
std::unique_ptr<Term> Term::parse(Lexer& lexer)
{
    auto term = std::make_unique<Term>();
    term->left = Factor::parse(lexer);

    while (true) {
        auto token = lexer.getCurrentToken();
        if (token.type != TokenType::MULTIPLY && token.type != TokenType::DIVIDE) {
            break;
        }

        term->op = token.value[0];
        lexer.advance();

        auto right = Factor::parse(lexer);
        auto new_term = std::make_unique<Term>();
        new_term->left = std::move(term->left);
        new_term->op = term->op;
        new_term->right = std::move(term);
        term = std::move(new_term);
    }

    return term;
}

// Factor parsing
std::unique_ptr<Factor> Factor::parse(Lexer& lexer)
{
    auto factor = std::make_unique<Factor>();

    auto token = lexer.getCurrentToken();
    if (token.type == TokenType::PLUS || token.type == TokenType::MINUS) {
        factor->sign = token.value[0];
        lexer.advance();
        factor->value = Primary::parse(lexer);
    }
    else {
        factor->sign = '\0';
        factor->value = Primary::parse(lexer);
    }

    return factor;
}

// Primary parsing
std::unique_ptr<Primary> Primary::parse(Lexer& lexer)
{
    auto primary = std::make_unique<Primary>();
    auto token = lexer.getCurrentToken();

    switch (token.type) {
        case TokenType::NUMBER:
            primary->type = Type::NUMBER;
            primary->value = token.value;
            lexer.advance();
            break;

        case TokenType::IDENTIFIER: {
            primary->type = Type::VARIABLE;
            primary->value = token.value;
            lexer.advance();

            // Check for component access (e.g., "v.x")
            if (lexer.getCurrentToken().type == TokenType::DOT) {
                lexer.advance();
                if (lexer.getCurrentToken().type != TokenType::IDENTIFIER) {
                    throw ParserError("Expected identifier after dot at position " +
                                      std::to_string(lexer.getCurrentToken().position));
                }
                primary->value += "." + lexer.getCurrentToken().value;
                lexer.advance();
            }
            break;
        }

        case TokenType::LPAREN:
            primary->type = Type::EXPRESSION;
            lexer.advance();
            primary->expr = Expression::parse(lexer);
            if (lexer.getCurrentToken().type != TokenType::RPAREN) {
                throw ParserError("Expected closing parenthesis at position " +
                                  std::to_string(lexer.getCurrentToken().position));
            }
            lexer.advance();
            break;

        default:
            throw ParserError("Unexpected token at position " +
                              std::to_string(token.position));
    }

    return primary;
}

// ToString implementations for AST nodes
std::string Expression::toString() const
{
    std::stringstream ss;
    ss << left->toString();
    if (right) {
        ss << " " << op << " " << right->toString();
    }
    return ss.str();
}

std::string Term::toString() const
{
    std::stringstream ss;
    ss << left->toString();
    if (right) {
        ss << " " << op << " " << right->toString();
    }
    return ss.str();
}

std::string Factor::toString() const
{
    std::stringstream ss;
    if (sign != '\0') {
        ss << sign;
    }
    ss << value->toString();
    return ss.str();
}

std::string Primary::toString() const
{
    switch (type) {
        case Type::NUMBER:
        case Type::VARIABLE:
            return value;
        case Type::EXPRESSION:
            return "(" + expr->toString() + ")";
        default:
            return "";
    }
}