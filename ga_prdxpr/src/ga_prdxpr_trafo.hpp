#pragma once

// author: Daniel Hug, 2025
//

// transformation of simple product expressions

#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <string>
#include <vector>

using namespace std::literals::string_literals;

// Forward declarations
class ast_node;
class Expression;
class Term;
class Factor;
class Primary;

///////////////////////////////////////////////////////////////////////////////
// Expression Validator class declaration
///////////////////////////////////////////////////////////////////////////////
class ExpressionValidator {
  public:

    static bool validateParentheses(std::string const& expr);
    static bool validateVariableName(std::string const& name);
    static bool validateOperatorSequence(std::string const& expr);

  private:

    static bool validateIdentifier(std::string const& name);
};


///////////////////////////////////////////////////////////////////////////////
// Token types
///////////////////////////////////////////////////////////////////////////////

enum class Token_t {
    NUMBER,
    IDENTIFIER,
    DOT,
    PLUS,
    MINUS,
    MULTIPLY,
    DEVIDE,
    LPAREN,
    RPAREN,
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
        case Token_t::DOT:
            return "DOT";
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
        case Token_t::END:
            return "END";
            break;
        default:
            return "UNKNOWN";
    }
}

///////////////////////////////////////////////////////////////////////////////
// Token structure
///////////////////////////////////////////////////////////////////////////////

struct Token {
    Token_t type{Token_t::END};
    std::string str_value{""s};
    size_t position{0};
    Token() = default;
    Token(Token_t t, std::string const& str_val = "", size_t pos = 0) :
        type(t), str_value(str_val), position(pos)
    {
    }
};


///////////////////////////////////////////////////////////////////////////////
// Lexer class
///////////////////////////////////////////////////////////////////////////////
class Lexer {
  public:

    explicit Lexer(const std::string& input);
    Token getNextToken();
    const Token& getCurrentToken() const { return currentToken; }
    void advance() { currentToken = getNextToken(); }

  private:

    std::string input;
    size_t position;
    Token currentToken;

    char peek() const;
    char advance_char();
    void skipWhitespace();
    Token parseNumber();
    Token parseIdentifier();
    bool isDigit(char c) const;
    bool isLetter(char c) const;
};


///////////////////////////////////////////////////////////////////////////////
// Abstract base class for expression nodes
///////////////////////////////////////////////////////////////////////////////

class ast_node {
  public:

    enum class ast_node_t { EXPRESSION, TERM, FACTOR, PRIMARY };

    virtual ~ast_node() = default;
    virtual std::string to_string() const = 0;
    virtual ast_node_t nodeType() const = 0;
    virtual std::string nodeType_to_string() const = 0;


    //   private:
};


///////////////////////////////////////////////////////////////////////////////
// Expression class for addition and subtraction
///////////////////////////////////////////////////////////////////////////////

class Expression : public ast_node {
  public:

    Expression(std::shared_ptr<ast_node> l, char o, std::shared_ptr<ast_node> r);

    static std::shared_ptr<ast_node> parse(Lexer& lexer);

    std::string to_string() const override;
    ast_node_t nodeType() const override;
    std::string nodeType_to_string() const override;

    //   private:

    std::shared_ptr<ast_node> left;
    char op; // '+' or '-'
    std::shared_ptr<ast_node> right;
};


///////////////////////////////////////////////////////////////////////////////
// Term class for multiplication
///////////////////////////////////////////////////////////////////////////////

class Term : public ast_node {
  public:

    Term(std::shared_ptr<ast_node> l, char o, std::shared_ptr<ast_node> r);

    static std::shared_ptr<ast_node> parse(Lexer& lexer);

    std::string to_string() const override;
    ast_node_t nodeType() const override;
    std::string nodeType_to_string() const override;

    //   private:

    std::shared_ptr<ast_node> left;
    char op; // '*' or '/'
    std::shared_ptr<ast_node> right;
};


///////////////////////////////////////////////////////////////////////////////
// Factor class for tight coupling of unary signs
///////////////////////////////////////////////////////////////////////////////

class Factor : public ast_node {
  public:

    Factor(char sign, std::shared_ptr<ast_node> v);

    static std::shared_ptr<ast_node> parse(Lexer& lexer);

    std::string to_string() const override;
    ast_node_t nodeType() const override;
    std::string nodeType_to_string() const override;

    //   private:

    char sign;                          // '+' or '-' or '\0'
    std::shared_ptr<ast_node> prim_val; // ptr to primary that holds the value
};


///////////////////////////////////////////////////////////////////////////////
// Primary class for numbers, variables, and expressions in parenthesis
///////////////////////////////////////////////////////////////////////////////

class Primary : public ast_node {
  public:

    Primary() = default;

    static std::shared_ptr<ast_node> parse(Lexer& lexer);

    std::string to_string() const override;
    ast_node_t nodeType() const override;
    std::string nodeType_to_string() const override;
    std::string primaryType_to_string() const;

    //   private:

    enum class Primary_t { NUMBER, VARIABLE, EXPRESSION };

    Primary_t type;
    std::string str_value; // for NUMBER and VARIABLE as string

    // inital value indicates not assigned to
    double num_value{std::numeric_limits<double>::min()}; // for NUMBER as value

    std::shared_ptr<ast_node> expr; // for EXPRESSION
};


///////////////////////////////////////////////////////////////////////////////
// Parser class
///////////////////////////////////////////////////////////////////////////////
class Parser {
  public:

    explicit Parser(std::string const& input);
    std::shared_ptr<ast_node> parse();

    static void validateBinaryOperation(Token const& op,
                                        std::shared_ptr<ast_node> const& left,
                                        std::shared_ptr<ast_node> const& right);

    //   private:

    Lexer lexer;
};


///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

std::string parse_and_print_ast(std::string const& input);
std::string parse_only(std::string const& input);
std::string parse_and_analyse(std::string const& input);


///////////////////////////////////////////////////////////////////////////////
// Printing functions
///////////////////////////////////////////////////////////////////////////////

void print_parse_tree(std::shared_ptr<ast_node> const& ast);
void print_expression_node(std::shared_ptr<Expression> const& ptr);
void print_term_node(std::shared_ptr<Term> const& ptr);
void print_factor_node(std::shared_ptr<Factor> const& ptr);
void print_primary_node(std::shared_ptr<Primary> const& ptr);


///////////////////////////////////////////////////////////////////////////////
// Expression analyzer class
///////////////////////////////////////////////////////////////////////////////

struct ExpressionAnalyzerTerm {
    std::string coefficient;              // The non-variable part
    std::map<std::string, int> variables; // Variables and their powers

    bool operator<(const ExpressionAnalyzerTerm& other) const
    {
        return variables < other.variables;
    }
};

class ExpressionAnalyzer {
  public:

    explicit ExpressionAnalyzer(std::shared_ptr<ast_node> const& ast);

    // Get terms grouped by variable combinations
    std::map<std::set<std::string>, std::vector<ExpressionAnalyzerTerm>>
    getGroupedTerms() const;

    // Get coefficients for specific variables
    std::map<std::string, double> getCoefficients(const std::string& varName) const;

    // Convert back to string with terms grouped by variables
    std::string toGroupedString() const;

  private:

    std::vector<ExpressionAnalyzerTerm> terms;

    void analyzeExpression(std::shared_ptr<Expression> const& expr, int sign = 1);
    void analyzeTerm(std::shared_ptr<Term> const& term, const std::string& coefficient,
                     int sign);
    void analyzeFactor(std::shared_ptr<Factor> const& factor, std::string& coefficient,
                       std::map<std::string, int>& variables, int sign);
    void analyzePrimary(std::shared_ptr<Primary> const& primary, std::string& coefficient,
                        std::map<std::string, int>& variables);

    static bool isVariable(const std::string& str, const std::string& varPrefix);
    static std::string combineCoefficients(const std::string& a, const std::string& b,
                                           char op);
};
