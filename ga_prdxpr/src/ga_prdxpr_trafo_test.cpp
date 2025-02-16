//
// author: Daniel Hug, 2025
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

#include <string>

#include "ga_prdxpr_trafo.hpp"

TEST_SUITE("expression transformation")
{

    // Expressions to test
    // 1.0 -> 1.0
    // a -> a
    // a*a*a + 1 -> a*a*a + 1
    // a + b -> a + b
    // a + a -> 2.0 * a
    // +a + a -> 2.0 * a
    // -a + a -> 0.0
    // a - a -> 0.0
    // 2*a + 2.0*a -> 4.0 * a
    // -a + 2*a -> a
    // 2*a - 3*a -> -a
    // -2*a + 4*a -> 2*a
    // -2*a + b + 4*a -> 2*a + b
    // -(2*a + b) + 4*a -> 2*a - b
    // +(2*a + b) + 4*a -> 6*a + b
    // 4*a + (2*a + b)  -> 6*a + b
    // 4*a - (2*a + b)  -> 2*a - b
    // +4*a - (2*a + b)  -> 2*a - b
    // (a + b)*(a - b) -> a*a - b*b
    // (a + b)*(a + b) -> a*a + a*b + b*a + b*b -> (for scalar a & b) -> a*a + 2*a*b + b*b

    TEST_CASE("basic expressions")
    {

        // multi-line strings must be R"( .... )";

        // tested from here
        // std::string s = R"(2.0)";
        // std::string s = R"(-2.0)";
        // std::string s = R"(a)";
        // std::string s = R"(-a)";
        // std::string s = R"(R.c2)";
        // std::string s = R"(-R.c2)";
        // std::string s = R"(-R.c2 + R.c3)";
        // std::string s = R"(-R.c2 - R.c3)";
        // std::string s = R"(-(R.c2 + R.c3))";
        // std::string s = R"(-(-R.c2 + R.c3))";
        // std::string s = R"(-(-R.c2 - R.c3))";

        // std::string s = R"(2.0 + 3.0)";
        // std::string s = R"(-2.0 + 3.0)";
        // std::string s = R"(-(-2.0 - 3.0))";
        // std::string s = R"(-(2.0 + 3.0))";
        // std::string s = R"(-(-2.0 + 3.0))";
        // std::string s = R"(-(-2.0 - 3.0))";

        // std::string s = R"(R.c2 + R.c3 + R.c0)";
        // std::string s = R"(-(R.c2 + R.c3 + R.c0))";
        // std::string s = R"(-(R.c2 - R.c3 + R.c0))";
        // std::string s = R"(-(R.c2 + R.c3 - R.c0))";

        // std::string s = R"(2.0 + 3.0 + 1.0)";
        // std::string s = R"(-2.0 + 3.0 + 1.0)";
        // std::string s = R"(-(-2.0 - 3.0 + 1.0))";
        // std::string s = R"(-(2.0 + 3.0) + 1.0)";
        // std::string s = R"(-(2.0 + 3.0) - 1.0)";
        // std::string s = R"(-(-2.0 + 3.0) + 1.0)";
        // std::string s = R"(-(-2.0 - 3.0) + 1.0)";


        // std::string s = R"(2.0 - 1.0)";
        // std::string s = R"(-2.0 + 1.0)";
        // std::string s = R"(a*a*a + 1)";
        std::string s = R"(a*b)";
        // std::string s = R"(-(a*b))";
        // std::string s = R"(-(a*a*a - 1))";
        // std::string s = R"(-(-1 + a*a*a - 1))";
        // std::string s = R"((a + b)*(a + b))";
        // std::string s = R"(-(a + b)*(a + b))";
        // std::string s = R"((a + b)*(a - b))";

        // still needs correct handling of unary minus (expression and term)
        // std::string s = R"(-(2.0 + 1.0))";
        // std::string s = R"(-(R.c2 + R.c3 + R.c0))";
        // std::string s = R"((a + b)*(a - b))";
        // std::string s = R"(-(a + b)*(a - b))";

        // std::string s = R"(-2.0 + 1.0)";
        // std::string s = R"(2.0 * 1.0)";
        // std::string s = R"(-2.0 * 1.0)";
        // std::string s = R"(2.0 * -1.0)";
        // std::string s = R"(-2.0 * -1.0)";
        // std::string s = R"(a + a)";
        // std::string s = R"(-a + a)";
        // std::string s = R"(a + 1)";
        // std::string s = R"(-a + 1)";
        // std::string s = R"(2.0*a + 1)";
        // std::string s = R"(-2.0*a + 1)";
        // std::string s = R"(R.c2 + R.c3 + R.c0)";
        // std::string s = R"(_a + b)";

        // std::string s = R"(( - R.c2 * v.z ) * R.c0)";
        // std::string s = R"(( - R.c2 * v.z ) * R.c0)";
        // std::string s = R"(( - R.c2 * v.z ) * R.c0)";

        // std::string s = R"((R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c0 +
        //                    (R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c3 -
        //                    (R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c2 +
        //                    (R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c1)";
        // std::string s = R"(( - R.c2 * v.z + R.c3 * v.y) * R.c0)";
        // std::string s = R"(( R.c2 * v.z + R.c3 * v.y) * R.c0)";


        std::string r = parseAndSimplify(s);

        fmt::println("");
        fmt::println("s: '{}'", s);
        fmt::println("");
        fmt::println("r: '{}'", r);
        fmt::println("");

        // std::vector<Term> terms = expression.getTerms();
        // for (auto const& t : terms) {
        //     fmt::println("term: '{}'", t.toString());
        // }
        // expression.groupByVariable("v.");
    }

} // TEST_SUITE("expression transformation")