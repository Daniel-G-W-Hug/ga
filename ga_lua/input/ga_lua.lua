-- ga_lua.lua - Main orchestrator for GA (Geometric Algebra) testing
-- This file coordinates testing of GA functionality across different algebras

-- Set global test precision constant
local eps = 1e-12

print("\nStarting comprehensive GA testing across all algebras...")
print("Test precision: eps =", eps)

-- Helper function to get script directory and build paths
-- script_dir is set by the C++ code to handle platform-specific paths
-- (Windows MSVC uses Debug/Release subdirs, others don't)

-- Test EGA 2D functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_ega2d.lua")

-- Demonstrate reference vs value semantics in GA objects
print("\n" .. string.rep("=", 60))
print("=== Demonstrating reference vs value semantics ===")
print("")

print("1. Reference semantics (default assignment behavior):")
print("----------------------------------------------------")
local u = vec2d.new(2,3)
local v = u  -- This creates an alias/reference
print("u = ", u)
print("v = ", v)
v.x = 2.3
print("after v.x = 2.3:")
print("u = ", u, " (changed - they share same object)")
print("v = ", v, " (changed - they share same object)")

print("")
print("2. Value semantics (using copy method):")
print("---------------------------------------")
local u2 = vec2d.new(2,3)
local v2 = u2:copy()  -- This creates an independent copy
print("u2 = ", u2)
print("v2 = ", v2)
v2.x = 2.3
print("after v2.x = 2.3:")
print("u2 = ", u2, " (unchanged - independent object)")
print("v2 = ", v2, " (changed - independent object)")

print("")
print("=== Summary ===")
print("Use local v = u:copy() for independent copies")
print("Use local v = u for shared references")

-- Test EGA 3D functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_ega3d.lua")

-- Test PGA 2DP functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_pga2dp.lua")

-- Test PGA 3DP functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_pga3dp.lua")

-- Test MIXED functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_application.lua")
