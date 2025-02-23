#ifndef PRDXPR_PARSER_HPP
#define PRDXPR_PARSER_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Forward declarations
class Expression;
class Term;
class Factor;
class Primary;

// Token types
enum class TokenType {
    NUMBER,
    IDENTIFIER,
    DOT,
    PLUS,
    MINUS,
    MULTIPLY,
    DIVIDE,
    LPAREN,
    RPAREN,
    END
};

// Token structure
struct Token {
    TokenType type;
    std::string value;
    size_t position;

    Token(TokenType t, const std::string& v, size_t pos) :
        type(t), value(v), position(pos)
    {
    }
};

// Lexer class
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

// Abstract Syntax Tree node classes
class ASTNode {
  public:

    virtual ~ASTNode() = default;
    virtual std::string toString() const = 0;
};

class Expression : public ASTNode {
  public:

    std::string toString() const override;
    static std::unique_ptr<Expression> parse(Lexer& lexer);

  private:

    std::unique_ptr<Term> left;
    char op; // '+' or '-'
    std::unique_ptr<Expression> right;
};

class Term : public ASTNode {
  public:

    std::string toString() const override;
    static std::unique_ptr<Term> parse(Lexer& lexer);

  private:

    std::unique_ptr<Factor> left;
    char op; // '*' or '/'
    std::unique_ptr<Term> right;
};

class Factor : public ASTNode {
  public:

    std::string toString() const override;
    static std::unique_ptr<Factor> parse(Lexer& lexer);

  private:

    char sign; // '+' or '-' or '\0'
    std::unique_ptr<Primary> value;
};

class Primary : public ASTNode {
  public:

    std::string toString() const override;
    static std::unique_ptr<Primary> parse(Lexer& lexer);

  private:

    enum class Type { NUMBER, VARIABLE, EXPRESSION };
    Type type;
    std::string value;                // for NUMBER and VARIABLE
    std::unique_ptr<Expression> expr; // for EXPRESSION
};

// Parser class
class Parser {
  public:

    explicit Parser(const std::string& input);
    std::unique_ptr<Expression> parse();

  private:

    Lexer lexer;
};

// Custom exception for parser errors
class ParserError : public std::runtime_error {
  public:

    explicit ParserError(const std::string& message) : std::runtime_error(message) {}
};

#endif // PRDXPR_PARSER_HPP