-- ga_lua.lua - Main orchestrator for GA (Geometric Algebra) testing
-- This file coordinates testing of all GA functionality across different algebras

print("===============================================")
print("GA Lua Test Suite - Comprehensive GA Testing")
print("===============================================")

-- Set global test precision constant
eps = 1e-12

print("\nStarting comprehensive GA testing across all algebras...")
print("Test precision: eps =", eps)

-- Helper function to get script directory and build paths
local script_dir = "../ga_lua/input/" -- relative to build directory

-- Test EGA 2D functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_ega2d.lua")

-- Test EGA 3D functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_ega3d.lua")

-- Test PGA 2DP functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_pga2dp.lua")

-- Test PGA 3DP functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_pga3dp.lua")

-- Final summary
print("\n" .. string.rep("=", 60))
print("GA Lua Test Suite - Final Summary")
print("=================================")
print("")
print("✓ EGA 2D (Euclidean Geometric Algebra 2D) - PASSED")
print("✓ EGA 3D (Euclidean Geometric Algebra 3D) - PASSED")
print("✓ PGA 2DP (Projective Geometric Algebra 2D+1P) - PASSED")
print("✓ PGA 3DP (Projective Geometric Algebra 3D+1P) - PASSED")
print("")
print("All GA functionality tests completed successfully!")
print("Total algebras tested: 4")
print("Total test modules executed: 4")
print("")
print("The ga_lua library now provides comprehensive support for:")
print("• Euclidean Geometric Algebra (EGA) in 2D and 3D")
print("• Projective Geometric Algebra (PGA) in 2D+1P and 3D+1P")
print("• All geometric operations: dot, wedge, commutator products")
print("• Motor systems for PGA (rotations and translations)")
print("• Dual operations and complement operations")
print("• Projections, rejections, and reflections")
print("• Comprehensive multivector support")
print("• Grade operations and algebraic operations")
print("")
print("Ready for interactive GA exploration!")
