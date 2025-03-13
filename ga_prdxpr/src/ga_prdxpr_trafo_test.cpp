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

    // Expressions to test (partly for case with simplification)
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
    // -4*a - (2*a + b)  -> -6*a - b
    // (a + b)*(a - b) -> a*a - b*b
    // (a + b)*(a + b) -> a*a + a*b + b*a + b*b -> (for scalar a & b) -> a*a + 2*a*b + b*b

    // multi-line strings must be R"( .... )";

    TEST_CASE("basic expressions - parse only, no simplifications")
    {

        CHECK(parse_only(R"(2.0)") == "2.0"s);
        CHECK(parse_only(R"(-2.0)") == "-2.0"s);
        CHECK(parse_only(R"(-(2.0))") == "-(2.0)"s);
        CHECK(parse_only(R"(-(-2.0))") == "-(-2.0)"s);

        CHECK(parse_only(R"(a)") == "a"s);
        CHECK(parse_only(R"((a))") == "(a)"s);
        CHECK(parse_only(R"(-a)") == "-a"s);
        CHECK(parse_only(R"(-(a))") == "-(a)"s);
        CHECK(parse_only(R"(-(-a))") == "-(-a)"s);

        CHECK(parse_only(R"(2 + 1)") == "2 + 1"s);
        CHECK(parse_only(R"(-(2 + 1))") == "-(2 + 1)"s);
        CHECK(parse_only(R"(-(2 - 1))") == "-(2 - 1)"s);
        CHECK(parse_only(R"(-(-2 + 1))") == "-(-2 + 1)"s);
        CHECK(parse_only(R"(-(-2 - 1))") == "-(-2 - 1)"s);
        CHECK(parse_only(R"(-2 - 1 + 3)") == "-2 - 1 + 3"s);
        CHECK(parse_only(R"(-(-2 - 1 + 3))") == "-(-2 - 1 + 3)"s);

        CHECK(parse_only(R"(a + 1)") == "a + 1"s);
        CHECK(parse_only(R"((a + 1))") == "(a + 1)"s);
        CHECK(parse_only(R"(((a + 1)))") == "((a + 1))"s);
        CHECK(parse_only(R"(-(a + 1))") == "-(a + 1)"s);
        CHECK(parse_only(R"(-(a - 1))") == "-(a - 1)"s);
        CHECK(parse_only(R"(-(-a - 1))") == "-(-a - 1)"s);

        CHECK(parse_only(R"(1 - a + 1)") == "1 - a + 1"s);
        CHECK(parse_only(R"((1 - a + 1))") == "(1 - a + 1)"s);
        CHECK(parse_only(R"((1 - (a) + 1))") == "(1 - (a) + 1)"s);
        CHECK(parse_only(R"(-(1 - a + 1))") == "-(1 - a + 1)"s);
        CHECK(parse_only(R"(-(-1 - a + 1))") == "-(-1 - a + 1)"s);
        CHECK(parse_only(R"(-(1 - a + 1 - 1))") == "-(1 - a + 1 - 1)"s);

        CHECK(parse_only(R"(R.c2)") == "R.c2"s);
        CHECK(parse_only(R"(-R.c2 + R.c3)") == "-R.c2 + R.c3"s);
        CHECK(parse_only(R"(-R.c2 - R.c3)") == "-R.c2 - R.c3"s);
        CHECK(parse_only(R"(-(R.c2 + R.c3))") == "-(R.c2 + R.c3)"s);
        CHECK(parse_only(R"(-(-R.c2 + R.c3))") == "-(-R.c2 + R.c3)"s);
        CHECK(parse_only(R"(-(-R.c2 - R.c3))") == "-(-R.c2 - R.c3)"s);

        CHECK(parse_only(R"(-(a*a - 1 + 1))") == "-(a * a - 1 + 1)"s);
        CHECK(parse_only(R"(-(-a*a + 1))") == "-(-a * a + 1)"s);
        CHECK(parse_only(R"(-(1 - a*a))") == "-(1 - a * a)"s);
        CHECK(parse_only(R"(-(1 + a*a))") == "-(1 + a * a)"s);
        CHECK(parse_only(R"(-(1 - a*a + 1))") == "-(1 - a * a + 1)"s);
        CHECK(parse_only(R"(-(1 - (a*a) + 1))") == "-(1 - (a * a) + 1)"s);
        CHECK(parse_only(R"(-(1 - a*a + 1 + a*a + 1))") ==
              "-(1 - a * a + 1 + a * a + 1)"s);
        CHECK(parse_only(R"(-(1 - a*a*a + 1 - a*a + 1))") ==
              "-(1 - a * a * a + 1 - a * a + 1)"s);

        CHECK(parse_only(R"(-(-(a) + 1))") == "-(-(a) + 1)"s);
        CHECK(parse_only(R"(-((-a) + 1))") == "-((-a) + 1)"s);
        CHECK(parse_only(R"(-((-a + 1) + 1))") == "-((-a + 1) + 1)"s);
        CHECK(parse_only(R"(-(-(a) + 1))") == "-(-(a) + 1)"s);
        CHECK(parse_only(R"(-(1 + (-a + 1)))") == "-(1 + (-a + 1))"s);

        CHECK(parse_only(R"(((a)))") == "((a))"s);
        CHECK(parse_only(R"((((a))))") == "(((a)))"s);
        CHECK(parse_only(R"(((((a)))))") == "((((a))))"s);
        CHECK(parse_only(R"(-((-((a)))))") == "-((-((a))))"s);
        CHECK(parse_only(R"(-(-(-((a)))))") == "-(-(-((a))))"s);
        CHECK(parse_only(R"(-(-(-(1+a))))") == "-(-(-(1 + a)))"s);
        CHECK(parse_only(R"(-(-(-(1+(a)))))") == "-(-(-(1 + (a))))"s);
        CHECK(parse_only(R"(-(((a))))") == "-(((a)))"s);
        CHECK(parse_only(R"(-((-(a))))") == "-((-(a)))"s);
        CHECK(parse_only(R"(-((-(-a))))") == "-((-(-a)))"s);

        CHECK(parse_only(R"(-((-2 - 1) + 3 - 1))") == "-((-2 - 1) + 3 - 1)"s);
        CHECK(parse_only(R"(-(-2 - 1 + (3 - 1)))") == "-(-2 - 1 + (3 - 1))"s);

        CHECK(parse_only(R"(-(-2 - 1 + (3 - 1)))") == "-(-2 - 1 + (3 - 1))"s);
        CHECK(parse_only(R"(-(a*a*a - 1))") == "-(a * a * a - 1)"s);
        CHECK(parse_only(R"(-(-1 + a*a*a - 1))") == "-(-1 + a * a * a - 1)"s);

        CHECK(parse_only(R"((a + b)*(a + b))") == "(a + b) * (a + b)"s);
        CHECK(parse_only(R"((a + b)*(a - b))") == "(a + b) * (a - b)"s);
        CHECK(parse_only(R"(-(a + b)*(a + b))") == "-(a + b) * (a + b)"s);
        CHECK(parse_only(R"(-(a + b)*a)") == "-(a + b) * a"s);
        CHECK(parse_only(R"((a + (-1 - 1 + 3)))") == "(a + (-1 - 1 + 3))"s);

        /// NOT validated down from here


        CHECK(parse_only(R"(2.0 - 1.0)") == "2.0 - 1.0"s);
        CHECK(parse_only(R"(-2.0 + 1.0)") == "-2.0 + 1.0"s);
        CHECK(parse_only(R"(a*a*a + 1)") == "a * a * a + 1"s);
        CHECK(parse_only(R"(2.0 + 3.0 + 1.0)") == "2.0 + 3.0 + 1.0"s);
        CHECK(parse_only(R"(-2.0 + 3.0 + 1.0)") == "-2.0 + 3.0 + 1.0"s);
        CHECK(parse_only(R"(-(-2.0 - 3.0) + 1.0)") == "-(-2.0 - 3.0) + 1.0"s);
        CHECK(parse_only(R"(-(-2.0 - 3.0 + 1.0 + 2.0))") == "-(-2.0 - 3.0 + 1.0 + 2.0)"s);
        CHECK(parse_only(R"(-(-2.0 - 3.0 + 1.0 + a*b))") ==
              "-(-2.0 - 3.0 + 1.0 + a * b)"s);
        CHECK(parse_only(R"(-(-2.0 - 3.0) + 1.0 + 2*(a*b))") ==
              "-(-2.0 - 3.0) + 1.0 + 2 * (a * b)"s);


        CHECK(parse_only(R"(a*b)") == "a * b"s);
        CHECK(parse_only(R"((a*b))") == "(a * b)"s);
        CHECK(parse_only(R"(-(a*b))") == "-(a * b)"s);
        CHECK(parse_only(R"((-a*b))") == "(-a * b)"s);
        CHECK(parse_only(R"(-a*(-b))") == "-a * (-b)"s);
        CHECK(parse_only(R"(-(-a*(-b)))") == "-(-a * (-b))"s);

        CHECK(parse_only(R"(-(2.0 + 1.0))") == "-(2.0 + 1.0)"s);
        CHECK(parse_only(R"(-(R.c2 + R.c3 + R.c0))") == "-(R.c2 + R.c3 + R.c0)"s);
        CHECK(parse_only(R"((a + b)*(a - b))") == "(a + b) * (a - b)"s);
        CHECK(parse_only(R"(-(a + b)*(a - b))") == "-(a + b) * (a - b)"s);

        CHECK(parse_only(R"(2.0 * 2.0)") == "2.0 * 2.0"s);
        CHECK(parse_only(R"(-2.0 * 2.0)") == "-2.0 * 2.0"s);
        CHECK(parse_only(R"(2.0 * (-3.0))") == "2.0 * (-3.0)"s);
        CHECK(parse_only(R"(-2.0 * (-6.0))") == "-2.0 * (-6.0)"s);
        CHECK(parse_only(R"(a + a)") == "a + a"s);
        CHECK(parse_only(R"(-a + a)") == "-a + a"s);
        CHECK(parse_only(R"(a + 1)") == "a + 1"s);
        CHECK(parse_only(R"(-a + 1)") == "-a + 1"s);
        CHECK(parse_only(R"(2.0*a + 1)") == "2.0 * a + 1"s);
        CHECK(parse_only(R"(-2.0*a + 1)") == "-2.0 * a + 1"s);
        CHECK(parse_only(R"(_a + b)") == "_a + b"s);

        CHECK(parse_only(R"(( - R.c2 * v.z ) * R.c0)") == "(-R.c2 * v.z) * R.c0"s);

        CHECK(parse_only("(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c0 + "
                         "(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c3 - "
                         "(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c2 + "
                         "(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c1"s) ==
              "(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c0 + "
              "(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c3 - "
              "(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c2 + "
              "(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c1"s);

        CHECK(parse_only("-(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c3 + "
                         "(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c0 + "
                         "(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c1 + "
                         "(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c2"s) ==
              "-(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c3 + "
              "(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c0 + "
              "(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c1 + "
              "(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c2"s);

        CHECK(parse_only("-(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c1 - "
                         "(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c2 - "
                         "(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c3 + "
                         "(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c0"s) ==
              "-(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c1 - "
              "(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c2 - "
              "(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c3 + "
              "(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c0"s);


        CHECK(parse_only(R"(( - R.c2 * v.z + R.c3 * v.y) * R.c0)") ==
              "(-R.c2 * v.z + R.c3 * v.y) * R.c0"s);
        CHECK(parse_only(R"(( R.c2 * v.z + R.c3 * v.y) * R.c0)") ==
              "(R.c2 * v.z + R.c3 * v.y) * R.c0"s);

        // std::string s = "-(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c1 - "
        //                 "(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c2 - "
        //                 "(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c3 + "
        //                 "(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c0"s;
        // std::string s = R"((-R.c2 * v.z + R.c3 * v.y) * R.c0)";

        // std::string s = R"(1.0)";
        // std::string s = R"(1.0 * 2.0)";
        //   std::string s = R"(1.0 + 2.0)";
        // std::string s = R"(-1.0)";
        // std::string s = R"(-(1.0 + 2))";
        // std::string s = R"(a)";
        // std::string s = R"(a * b)";
        // std::string s = R"(a + b)";
        // std::string s = R"(a*(a + b))";
        // std::string s = R"(-a*(a + b))";
        // std::string s = R"(-(a*(a + b)))";
        // std::string s = R"(-(a*(a + b)) + 3*(a+b))";
        // std::string s = R"((a + b)*a)";
        // std::string s = R"((a + b)*(a + b))";
        // std::string s = "R.c0 * v.x - R.c2 * v.z + R.c3 * v.y + R.c1 * v.x"s;

        //   std::string s = "(R.c0 * v.x + R.c1 * v.y) * R.c0";
        // std::string s = "-(R.c0 * v.x + R.c1 * v.y) * R.c0";
        // std::string s = "(R.c0 * v.x + R.c1 * v.y + R.c2 * v.z) * R.c0";

        // ega2d:
        std::string s = "(R.c0 * v.x + R.c1 * v.y) * R.c0 + "
                        "(R.c0 * v.y - R.c1 * v.x) * R.c1"s;
        //   std::string s = "-(R.c0 * v.x + R.c1 * v.y) * R.c1 + "
        //                   "(R.c0 * v.y - R.c1 * v.x) * R.c0 "s;

        // ega3d:
        //   std::string s = "(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c0"
        //                   " + (R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c3"
        //                   " - (R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c2"
        //                   " + (R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c1"s;
        //   std::string s = "-(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c3"
        //                   " + (R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c0"
        //                   " + (R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c1"
        //                   " + (R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c2"s;
        //   std::string s = "(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c2"
        //                   " - (R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c1"
        //                   " + (R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c0"
        //                   " + (R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c3"s;
        //   std::string s = "-(R.c0 * v.x - R.c2 * v.z + R.c3 * v.y) * R.c1"
        //                   " -(R.c0 * v.y + R.c1 * v.z - R.c3 * v.x) * R.c2"
        //                   " -(R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c3"
        //                   " +(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c0"s;

        // pga2dp regressive:
        //   std::string s = "-(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c1"
        //                   " + (-R.c2 * v.z) * R.c0 - (R.c1 * v.z) * R.c3"
        //                   " + (-R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c2"s;
        //   std::string s = "-(R.c1 * v.x + R.c2 * v.y + R.c3 * v.z) * R.c2"
        //                   " + (-R.c2 * v.z) * R.c3 + (R.c1 * v.z) * R.c0"
        //                   " - (-R.c0 * v.z - R.c1 * v.y + R.c2 * v.x) * R.c1"s;
        //     std::string s = "-(-R.c2 * v.z) * R.c2 + (R.c1 * v.z) * R.c1"s;
        //   std::string s = "-(-R.c2 * v.z) * R.c1 - (R.c1 * v.z) * R.c2"s;


        fmt::println("Initial string to parse: s: '{}'\n", s);

        std::string r1 = parse_and_print_ast(s);
        fmt::println("s: '{}'\n", s);
        fmt::println("r: '{}'", r1);
        // std::string r2 = parse_and_print_ast(r1);
        // fmt::println("r: '{}'", r2);
        fmt::println("");

        std::string r = parse_and_analyse(s);
        fmt::println("s: '{}'\n", s);
        fmt::println("r: '{}'", r);
        fmt::println("");
    }

} // TEST_SUITE("expression transformation")