#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo_expr_simplifier.hpp"
#include "ga_prdxpr_trafo_nary_expression.hpp"
#include "ga_prdxpr_trafo_sandwich_transformer.hpp"

///////////////////////////////////////////////////////////////////////////////
// Test functions for transformation system development and validation
///////////////////////////////////////////////////////////////////////////////

class SimplificationTests {
  public:

    // Test EGA2D vector sandwich expansion
    // Input: "(R.c0 * v.x + R.c1 * v.y) * R.c0 + (R.c0 * v.y - R.c1 * v.x) * R.c1"
    // Expected: "R.c0*R.c0*v.x + R.c1*R.c0*v.y + R.c0*R.c1*v.y - R.c1*R.c1*v.x"
    static void testEGA2DExpansion();

    // Test pattern-based variable ordering system
    static void testMotorCoefficientOrdering();

    // Test complete transformation pipeline with patterns
    static void testPatternAwareTransformation();

    // Test all predefined algebra patterns
    static void testPredefinedPatterns();

    // Test commutativity application
    // Input: "R.c1*R.c0*v.y + R.c0*R.c1*v.y"
    // Expected: "2*R.c0*R.c1*v.y"
    static void testCommutativeCollection();

    // Test symmetry cancellation
    // Input: "R.c0*v.x*R.c1 - R.c1*v.x*R.c0"
    // Expected: "0" (should cancel due to commutativity)
    static void testSymmetryCancellation();

    // Test complete EGA2D transformation
    // Verify against manual example results
    static void testCompleteEGA2DTransformation();

    // Run all available tests
    static void runAllTests();
};

///////////////////////////////////////////////////////////////////////////////
// Integration test cases for sandwich transformer
///////////////////////////////////////////////////////////////////////////////

class SandwichTransformationTests {
  public:

    // Test EGA2D vector sandwich (matches manual example)
    static void testEGA2DVectorTransformation();

    // Test EGA3D vector sandwich
    static void testEGA3DVectorTransformation();

    // Test coefficient isolation and matrix generation
    static void testCoefficientExtraction();

    // Test complete pipeline: ga_prdxpr output -> simplified matrix form
    static void testCompletePipeline();

    // Run all transformation tests
    static void runAllTests();
};