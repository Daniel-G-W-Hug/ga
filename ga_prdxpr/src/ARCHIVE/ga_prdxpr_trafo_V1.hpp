#pragma once

// author: Daniel Hug, 2025
//

// transformation of simple product expressions

#include <limits>
#include <map>
#include <memory>
#include <stack>
#include <string>
#include <vector>

using namespace std::literals::string_literals;

// helpers to handle unary minus in front of expressions
void reset_sign_tracking_vectors();
void print_sign_tracking_vectors();
void mark_is_negative_on_level(size_t lvl);
void check_for_carry_over_on_level(size_t lvl);
bool minus_in_flight_on_level(size_t lvl);

// Forward declarations
class ExprNode;
class Expression;
class Term;
class Factor;
class Primary;

// Expression Validator class declaration
class ExpressionValidator {
  public:

    static bool validateParentheses(std::string const& expr);
    static bool validateVariableName(std::string const& name);
    static bool validateOperatorSequence(std::string const& expr);

  private:

    static bool validateIdentifier(std::string const& name);
};

// Token types
enum class Token_t {
    NUMBER,
    IDENTIFIER,
    PLUS,
    MINUS,
    MULTIPLY,
    DEVIDE,
    LPAREN,
    RPAREN,
    DOT,
    END
};

inline constexpr std::string Token_t_toString(Token_t t)
{
    switch (t) {
        case Token_t::NUMBER:
            return "NUMBER";
            break;
        case Token_t::IDENTIFIER:
            return "IDENTIFIER";
            break;
        case Token_t::PLUS:
            return "PLUS";
            break;
        case Token_t::MINUS:
            return "MINUS";
            break;
        case Token_t::MULTIPLY:
            return "MULTIPLY";
            break;
        case Token_t::DEVIDE:
            return "DEVIDE";
            break;
        case Token_t::LPAREN:
            return "LPAREN";
            break;
        case Token_t::RPAREN:
            return "RPAREN";
            break;
        case Token_t::DOT:
            return "DOT";
            break;
        case Token_t::END:
            return "END";
            break;
        default:
            return "UNKNOWN";
    }
}

// Token structure
struct Token {
    Token_t type{Token_t::END};
    std::string str_value{""s};
    size_t position{0};
    Token() = default;
    Token(Token_t t, std::string const& v = "", size_t pos = 0);
};

// Lexer class
class Lexer {
  private:

    std::string input;
    size_t position;

    char peek() const;
    char advance();
    void skipWhitespace();
    Token readNumber();
    Token readIdentifier();

  public:

    explicit Lexer(std::string const& input);
    Token nextToken();
};

enum class ExprNode_t { EXPRESSION, TERM, PRIMARY };

inline constexpr std::string ExprNode_t_toString(ExprNode_t t)
{
    switch (t) {
        case ExprNode_t::EXPRESSION:
            return "EXPRESSION";
            break;
        case ExprNode_t::TERM:
            return "TERM";
            break;
        case ExprNode_t::PRIMARY:
            return "PRIMARY";
            break;
        default:
            return "UNKNOWN";
    }
};

// Abstract base class for expression nodes
class ExprNode {
  public:

    virtual ~ExprNode() = default;
    virtual std::string toString(bool use_braces = false) const = 0;
    virtual ExprNode_t nodeType() const = 0;
    virtual std::shared_ptr<ExprNode> simplify() = 0;
};

// Expression class for addition and subtraction
class Expression : public ExprNode {
  public:

    std::shared_ptr<ExprNode> left;
    char op;
    std::shared_ptr<ExprNode> right;
    bool str_value_starts_with_minus{false}; // for starting with unary minus

    Expression(std::shared_ptr<ExprNode> l, char o, std::shared_ptr<ExprNode> r);
    std::string toString(bool use_braces = false) const override;
    ExprNode_t nodeType() const override;
    std::shared_ptr<ExprNode> simplify() override;
};

// Term class for multiplication
class Term : public ExprNode {
  public:

    std::shared_ptr<ExprNode> left;
    std::shared_ptr<ExprNode> right;
    bool str_value_starts_with_minus{false}; // for starting with unary minus

    Term(std::shared_ptr<ExprNode> l, std::shared_ptr<ExprNode> r);
    std::string toString(bool use_braces = false) const override;
    ExprNode_t nodeType() const override;
    std::shared_ptr<ExprNode> simplify() override;
};

enum class Primary_t { NUMBER, VARIABLE, EXPRESSION };

inline constexpr std::string Primary_t_toString(Primary_t t)
{
    switch (t) {
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

// Primary class for numbers, variables, and parenthesized expressions
class Primary : public ExprNode {
  public:

    Primary_t type;
    std::string str_value;
    double num_value{std::numeric_limits<double>::min()}; // indicates not assigned to
    bool str_value_starts_with_minus{false};

    Primary(std::string const& str_val, Primary_t t);
    std::string toString(bool use_braces = false) const override;
    ExprNode_t nodeType() const override;
    std::shared_ptr<ExprNode> simplify() override;
};

void toggle_sign_of_primary(std::shared_ptr<Primary> ptr);
void resolve_signs_in_expression_with_primaries(std::shared_ptr<Expression> ptr);

// Parser class
class Parser {
  private:

    Lexer lexer;
    Token currentToken;
    std::stack<Token> operatorStack;

    void advance();
    void validateBinaryOperation(Token const& op, std::shared_ptr<ExprNode> const& left,
                                 std::shared_ptr<ExprNode> const& right);
    std::shared_ptr<ExprNode> parsePrimary(size_t lvl);
    std::shared_ptr<ExprNode> parseTerm(size_t lvl);
    std::shared_ptr<ExprNode> parseExpression(size_t lvl);

  public:

    explicit Parser(std::string const& input);
    std::shared_ptr<ExprNode> parse();
};

// Helper function
std::string parse_and_print_ast(std::string const& input);
std::string parse_no_print_ast(std::string const& input);

void print_parse_tree(std::shared_ptr<ExprNode> const& ast);
void print_expression_node(std::shared_ptr<Expression> const& ptr);
void print_term_node(std::shared_ptr<Term> const& ptr);
void print_primary_node(std::shared_ptr<Primary> const& ptr);
