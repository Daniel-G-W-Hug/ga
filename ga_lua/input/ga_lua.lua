-- ga_lua.lua - Main orchestrator for GA (Geometric Algebra) testing
-- This file coordinates testing of GA functionality across different algebras

-- Set global test precision constant
local eps = 1e-12

print("\nStarting comprehensive GA testing across all algebras...")
print("Test precision: eps =", eps)

-- Helper function to get script directory and build paths
-- local script_dir = "../ga_lua/input/" -- relative to build directory
local script_dir = "../../ga_lua/input/" -- relative to build/ga_lua directory

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

-- Test MIXED functionality
print("\n" .. string.rep("=", 60))
dofile(script_dir .. "test_application.lua")
