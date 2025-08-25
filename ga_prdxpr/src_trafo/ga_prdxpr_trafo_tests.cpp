// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo_tests.hpp"
#include "ga_prdxpr_trafo.hpp"

#include <fmt/format.h>

///////////////////////////////////////////////////////////////////////////////
// Test implementation for validation
///////////////////////////////////////////////////////////////////////////////

void SimplificationTests::testEGA2DExpansion()
{
    fmt::println("=== Testing EGA2D Expression Expansion ===");

    std::string input =
        "(R.c0 * v.x + R.c1 * v.y) * R.c0 + (R.c0 * v.y - R.c1 * v.x) * R.c1";
    fmt::println("Input: {}", input);

    Parser parser(input);
    auto ast = parser.parse();

    auto simplified = ExpressionSimplifier::simplify(ast);

    fmt::println("Simplified: {}", simplified->to_string());
    fmt::println("");
}

void SimplificationTests::testMotorCoefficientOrdering()
{
    fmt::println("=== Testing Pattern-Based Variable Ordering ===");

    // Test comprehensive pattern-based sorting
    std::map<std::string, int> test_factors = {{"v.x", 1}, {"M.c2", 1}, {"M.c0", 1},
                                               {"B.y", 1}, {"M.c1", 1}, {"B.vx", 1},
                                               {"t.z", 1}, {"v.w", 1}};

    fmt::println("Original factors: v.x, M.c2, M.c0, B.y, M.c1, B.vx, t.z, v.w");

    // Test EGA3D patterns (R.c coefficients, B.x/B.y/B.z bivectors)
    auto ega3d_patterns = GeometricVariablePatterns::createEGA3DPatterns();
    auto sorted_ega3d =
        GAAlgebraRules::getSortedVariablePairs(test_factors, ega3d_patterns);
    fmt::print("EGA3D ordering: ");
    for (const auto& [var, power] : sorted_ega3d) {
        fmt::print("{} ", var);
    }
    fmt::println("");

    // Test PGA3DP patterns (M.c coefficients, B.vx/B.vy/etc. bivectors)
    auto pga3dp_patterns = GeometricVariablePatterns::createPGA3DPPatterns();
    auto sorted_pga3dp =
        GAAlgebraRules::getSortedVariablePairs(test_factors, pga3dp_patterns);
    fmt::print("PGA3DP ordering: ");
    for (const auto& [var, power] : sorted_pga3dp) {
        fmt::print("{} ", var);
    }
    fmt::println("");

    // Test backward compatibility with simple coefficient prefix
    GeometricVariablePatterns simple_motor_patterns;
    simple_motor_patterns.coeff_prefix = "M.c";
    auto sorted_motor_compat =
        GAAlgebraRules::getSortedVariablePairs(test_factors, simple_motor_patterns);
    fmt::print("Motor compat ordering: ");
    for (const auto& [var, power] : sorted_motor_compat) {
        fmt::print("{} ", var);
    }
    fmt::println("");

    fmt::println(
        "Expected: Coefficients first, then vectors, then bivectors, then trivectors");
    fmt::println("");
}

void SimplificationTests::testPatternAwareTransformation()
{
    fmt::println("=== Testing Pattern-Aware Transformation Pipeline ===");

    // Test expression with both R.c and M.c coefficients
    std::string test_expr = "M.c0 * v.x + R.c1 * B.vx - M.c2 * t.z";
    fmt::println("Test expression: {}", test_expr);

    // Test with different pattern configurations
    fmt::println("\nDefault patterns (R.c priority):");
    Parser parser1(test_expr);
    auto ast1 = parser1.parse();
    auto expanded1 = ExpressionSimplifier::expandProducts(ast1);
    auto nary1 = NAryConverter::fromBinaryAST(expanded1);
    nary1.normalizeSignsAndCommutativity(); // Default patterns
    fmt::println("Result: {}",
                 nary1.toString(GeometricVariablePatterns{})); // Use default patterns

    fmt::println("\nMotor patterns (M.c priority):");
    GeometricVariablePatterns motor_patterns;
    motor_patterns.coeff_prefix = "M.c";
    Parser parser2(test_expr);
    auto ast2 = parser2.parse();
    auto expanded2 = ExpressionSimplifier::expandProducts(ast2);
    auto nary2 = NAryConverter::fromBinaryAST(expanded2);
    nary2.normalizeSignsAndCommutativity(motor_patterns);       // Motor patterns
    fmt::println("Result: {}", nary2.toString(motor_patterns)); // Use motor patterns

    fmt::println("\nPGA3DP patterns (complete motor configuration):");
    auto pga3dp_patterns = GeometricVariablePatterns::createPGA3DPPatterns();
    Parser parser3(test_expr);
    auto ast3 = parser3.parse();
    auto expanded3 = ExpressionSimplifier::expandProducts(ast3);
    auto nary3 = NAryConverter::fromBinaryAST(expanded3);
    nary3.normalizeSignsAndCommutativity(pga3dp_patterns);       // PGA3DP patterns
    fmt::println("Result: {}", nary3.toString(pga3dp_patterns)); // Use PGA3DP patterns

    fmt::println("Expected: Motor coefficients should sort differently based on pattern "
                 "configuration");

    // Add realistic test cases to demonstrate practical usage
    fmt::println("\n=== Realistic Usage Examples ===");

    // EGA2D realistic test (rotors only) - deliberately wrong order in input
    fmt::println("\n1. EGA2D Rotor Expression (R.c coefficients only):");
    std::string ega_expr =
        "v.x * R.c2 + B.x * R.c0 - v.y * R.c1"; // Wrong order: geometric vars first
    fmt::println("Input (wrong order): {}", ega_expr);

    Parser ega_parser(ega_expr);
    auto ega_ast = ega_parser.parse();
    auto ega_expanded = ExpressionSimplifier::expandProducts(ega_ast);
    auto ega_nary = NAryConverter::fromBinaryAST(ega_expanded);

    // Test with EGA2D patterns
    auto ega2d_patterns = GeometricVariablePatterns::createEGA2DPatterns();
    ega_nary.normalizeSignsAndCommutativity(ega2d_patterns);
    fmt::println("EGA2D result (corrected): {}", ega_nary.toString(ega2d_patterns));

    // PGA2DP realistic test (motors only) - deliberately wrong order in input
    fmt::println("\n2. PGA2DP Motor Expression (M.c coefficients only):");
    std::string pga_expr =
        "v.z * M.c3 + B.x * M.c1 + v.x * M.c2 - v.y * M.c0"; // Wrong order: geometric
                                                             // vars first, coeffs mixed
    fmt::println("Input (wrong order): {}", pga_expr);

    Parser pga_parser(pga_expr);
    auto pga_ast = pga_parser.parse();
    auto pga_expanded = ExpressionSimplifier::expandProducts(pga_ast);
    auto pga_nary = NAryConverter::fromBinaryAST(pga_expanded);

    // Test with PGA2DP patterns
    auto pga2dp_patterns = GeometricVariablePatterns::createPGA2DPPatterns();
    pga_nary.normalizeSignsAndCommutativity(pga2dp_patterns);
    fmt::println("PGA2DP result (corrected): {}", pga_nary.toString(pga2dp_patterns));

    fmt::println(
        "\nNote: These examples demonstrate the pattern system automatically reordering");
    fmt::println("expressions from wrong order (geometric vars first) to correct order "
                 "(coefficients first).");
    fmt::println("");
}

void SimplificationTests::testPredefinedPatterns()
{
    fmt::println("=== Testing Predefined Algebra Patterns ===");

    // Test factors covering all algebra types
    std::map<std::string, int> comprehensive_factors = {// Coefficients
                                                        {"R.c0", 1},
                                                        {"R.c1", 1},
                                                        {"M.c0", 1},
                                                        {"M.c2", 1},
                                                        // Vectors (various algebras)
                                                        {"v.x", 1},
                                                        {"v.y", 1},
                                                        {"v.z", 1},
                                                        {"v.w", 1},
                                                        // Bivectors (EGA and PGA styles)
                                                        {"B.x", 1},
                                                        {"B.y", 1},
                                                        {"B.z", 1},
                                                        {"B.vx", 1},
                                                        {"B.vy", 1},
                                                        {"B.mx", 1},
                                                        // Trivectors
                                                        {"t.x", 1},
                                                        {"t.y", 1},
                                                        {"t.w", 1}};

    fmt::println("Input factors: R.c0, R.c1, M.c0, M.c2, v.x, v.y, v.z, v.w, B.x, B.y, "
                 "B.z, B.vx, B.vy, B.mx, t.x, t.y, t.w");
    fmt::println("");

    // Test each predefined algebra pattern
    fmt::println("EGA2D patterns:");
    auto ega2d = GeometricVariablePatterns::createEGA2DPatterns();
    auto sorted_ega2d =
        GAAlgebraRules::getSortedVariablePairs(comprehensive_factors, ega2d);
    fmt::print("  Coeff prefix: {} | Vectors: ", ega2d.coeff_prefix);
    for (const auto& v : ega2d.vectors)
        fmt::print("{} ", v);
    fmt::print("| Ordering: ");
    for (const auto& [var, power] : sorted_ega2d)
        fmt::print("{} ", var);
    fmt::println("");

    fmt::println("EGA3D patterns:");
    auto ega3d = GeometricVariablePatterns::createEGA3DPatterns();
    auto sorted_ega3d =
        GAAlgebraRules::getSortedVariablePairs(comprehensive_factors, ega3d);
    fmt::print("  Coeff prefix: {} | Vectors: ", ega3d.coeff_prefix);
    for (const auto& v : ega3d.vectors)
        fmt::print("{} ", v);
    fmt::print("| Bivectors: ");
    for (const auto& b : ega3d.bivectors)
        fmt::print("{} ", b);
    fmt::print("| Ordering: ");
    for (const auto& [var, power] : sorted_ega3d)
        fmt::print("{} ", var);
    fmt::println("");

    fmt::println("PGA2DP patterns:");
    auto pga2dp = GeometricVariablePatterns::createPGA2DPPatterns();
    auto sorted_pga2dp =
        GAAlgebraRules::getSortedVariablePairs(comprehensive_factors, pga2dp);
    fmt::print("  Coeff prefix: {} | Vectors: ", pga2dp.coeff_prefix);
    for (const auto& v : pga2dp.vectors)
        fmt::print("{} ", v);
    fmt::print("| Bivectors: ");
    for (const auto& b : pga2dp.bivectors)
        fmt::print("{} ", b);
    fmt::print("| Ordering: ");
    for (const auto& [var, power] : sorted_pga2dp)
        fmt::print("{} ", var);
    fmt::println("");

    fmt::println("PGA3DP patterns:");
    auto pga3dp = GeometricVariablePatterns::createPGA3DPPatterns();
    auto sorted_pga3dp =
        GAAlgebraRules::getSortedVariablePairs(comprehensive_factors, pga3dp);
    fmt::print("  Coeff prefix: {} | Vectors: ", pga3dp.coeff_prefix);
    for (const auto& v : pga3dp.vectors)
        fmt::print("{} ", v);
    fmt::print("| Bivectors: ");
    for (const auto& b : pga3dp.bivectors)
        fmt::print("{} ", b);
    fmt::print("| Trivectors: ");
    for (const auto& t : pga3dp.trivectors)
        fmt::print("{} ", t);
    fmt::print("| Ordering: ");
    for (const auto& [var, power] : sorted_pga3dp)
        fmt::print("{} ", var);
    fmt::println("");

    fmt::println("Expected: Each algebra prioritizes its own coefficient prefix and "
                 "recognizes its geometric variable patterns");
    fmt::println("");
}

void SimplificationTests::testCommutativeCollection()
{
    // TODO: Implement commutative collection test
    fmt::println("=== Testing Commutative Collection ===");
    fmt::println(
        "TODO: Implement test for R.c1*R.c0*v.y + R.c0*R.c1*v.y -> 2*R.c0*R.c1*v.y");
    fmt::println("");
}

void SimplificationTests::testSymmetryCancellation()
{
    // TODO: Implement symmetry cancellation test
    fmt::println("=== Testing Symmetry Cancellation ===");
    fmt::println("TODO: Implement test for R.c0*v.x*R.c1 - R.c1*v.x*R.c0 -> 0");
    fmt::println("");
}

void SimplificationTests::testCompleteEGA2DTransformation()
{
    // TODO: Implement complete EGA2D transformation test
    fmt::println("=== Testing Complete EGA2D Transformation ===");
    fmt::println("TODO: Implement complete pipeline test");
    fmt::println("");
}

void SimplificationTests::runAllTests()
{
    fmt::println("=== Running All Simplification Tests ===");
    testEGA2DExpansion();
    testMotorCoefficientOrdering();
    testPatternAwareTransformation();
    testPredefinedPatterns();
    testCommutativeCollection();
    testSymmetryCancellation();
    testCompleteEGA2DTransformation();
    fmt::println("=== All Simplification Tests Complete ===\n");
}

///////////////////////////////////////////////////////////////////////////////
// SandwichTransformationTests implementation
///////////////////////////////////////////////////////////////////////////////

void SandwichTransformationTests::testEGA2DVectorTransformation()
{
    // TODO: Implement EGA2D vector transformation test
    fmt::println("=== Testing EGA2D Vector Transformation ===");
    fmt::println("TODO: Implement sandwich transformation test");
    fmt::println("");
}

void SandwichTransformationTests::testEGA3DVectorTransformation()
{
    // TODO: Implement EGA3D vector transformation test
    fmt::println("=== Testing EGA3D Vector Transformation ===");
    fmt::println("TODO: Implement sandwich transformation test");
    fmt::println("");
}

void SandwichTransformationTests::testCoefficientExtraction()
{
    // TODO: Implement coefficient extraction test
    fmt::println("=== Testing Coefficient Extraction ===");
    fmt::println("TODO: Implement matrix coefficient extraction test");
    fmt::println("");
}

void SandwichTransformationTests::testCompletePipeline()
{
    // TODO: Implement complete pipeline test
    fmt::println("=== Testing Complete Pipeline ===");
    fmt::println("TODO: Implement ga_prdxpr output -> simplified matrix form test");
    fmt::println("");
}

void SandwichTransformationTests::runAllTests()
{
    fmt::println("=== Running All Transformation Tests ===");
    testEGA2DVectorTransformation();
    testEGA3DVectorTransformation();
    testCoefficientExtraction();
    testCompletePipeline();
    fmt::println("=== All Transformation Tests Complete ===\n");
}