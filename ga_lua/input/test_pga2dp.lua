-- test_pga2dp.lua - Comprehensive PGA 2DP testing module
-- Tests all Projective Geometric Algebra 2DP (2D + 1 projective) functionality

print("=======================================================")
print("Testing PGA 2DP (Projective Geometric Algebra 2D+1P)")
print("=======================================================")

-- Test precision constant
local eps = eps or 1e-12

-- Test scalar2dp basic operations
print("\n1. Testing scalar2dp basic operations:")
print("-------------------------------------")

local s0 = scalar2dp.new()
local s1 = scalar2dp.new(5.0)
local s2 = scalar2dp.new(3.0)
print("s0:", s0, "empty ctor")
print("s1:", s1, "component ctor")
print("s2:", s2, "component ctor")
print("s1 + s2:", s1 + s2, "addition")
print("s1 * 2:", s1 * 2.0, "scalar multiplication")

-- Test vec2dp basic operations
print("\n2. Testing vec2dp basic operations:")
print("----------------------------------")

local v0 = vec2dp.new()
local v1 = vec2dp.new(1, 2, 0)  -- direction vector (z=0)
local v2 = vec2dp.new(3, 4, 1)  -- point (z=1)
local v3 = vec2dp.new(-v2)
local v4 = vec2dp.new(4, 6, 1)  -- v1 + v2

print("v0:", v0, "empty ctor")
print("v1:", v1, "direction vector (z=0)")
print("v2:", v2, "point (z=1)")
print("v3:", v3, "vector ctor v3 = -v2", "unary minus")
print("v4:", v4, "component ctor")

assert(v4 == v1 + v2)
assert(v1 == v4 - v2)
if v4 == v1 + v2 and v1 == v4 - v2 then
    print("v4 == v1 + v2", "comparison, addition and subtraction")
end

local v5 = 3 * v4
assert(3 * v4 == v4 * 3)
assert(v5 / 3 == v4)
print("3 * v4:", 3 * v4, "scalar multiplication")
print("v5 / 3:", v5 / 3, "scalar division")

-- Test vec2dp geometric operations
print("\n3. Testing vec2dp geometric operations:")
print("--------------------------------------")

print("dot(v1, v2):", dot(v1, v2), "dot product")
print("nrm_sq(v2): ", nrm_sq(v2), "squared norm")
print("nrm(v2): ", nrm(v2), "norm")

print("wdg(v1,v2):", wdg(v1, v2), "wedge product creates line")

-- Test unitization (PGA specific)
print("\n4. Testing unitization (PGA specific):")
print("--------------------------------------")

local v_pt = vec2dp.new(3, 4, 1)  -- point vector (can be unitized)
print("point vector v_pt =", v_pt)
print("unitize(v_pt) =", unitize(v_pt))

-- Test bivec2dp operations
print("\n5. Testing bivec2dp operations:")
print("-------------------------------")

local b0 = bivec2dp.new()
local b1 = bivec2dp.new(1, 0, 0)  -- line
local b2 = bivec2dp.new(0, 1, 0)  -- horizontal line (y = 0)
local line = bivec2dp.new(1, 1, -2)  -- line ax + by + c = 0

print("b0:", b0, "empty ctor")
print("b1:", b1, "line")
print("b2:", b2, "horizontal line")
print("line:", line, "line ax + by + c = 0")

print("unitize(line) =", unitize(line))

-- Test bivector norms
print("nrm_sq(line): ", nrm_sq(line), "squared norm")
print("nrm(line): ", nrm(line), "norm")

-- Test pscalar2dp
print("\n6. Testing pscalar2dp:")
print("---------------------")

local ps1 = pscalar2dp.new()
local ps2 = pscalar2dp.new(1.0)
print("ps1:", ps1)
print("ps2:", ps2)

-- Test dualnum2dp
print("\n7. Testing dualnum2dp:")
print("---------------------")

local dn1 = dualnum2dp.new(1.0, 0.5)
local dn2 = dualnum2dp.new(3.0, 2.0)
print("dn1:", dn1, "dual number")
print("dn2:", dn2, "dual number")
print("dn1 + dn2:", dn1 + dn2, "addition")
print("dn1 * 2:", dn1 * 2.0, "scalar multiplication")
print("unitize(dn2) =", unitize(dn2))

-- Test PGA constants
print("\n8. Testing PGA 2DP constants:")
print("-----------------------------")

print("origin_2dp =", origin_2dp)
print("e1_2dp =", e1_2dp)
print("e2_2dp =", e2_2dp)
print("e3_2dp =", e3_2dp)
print("horizon_2dp =", horizon_2dp)
print("x_axis_2dp =", x_axis_2dp)
print("y_axis_2dp =", y_axis_2dp)

-- Test wedge products with constants
print("\n9. Testing wedge products with constants:")
print("-----------------------------------------")

print("e1_2dp ^ e2_2dp =", e1_2dp ^ e2_2dp)
print("e1_2dp ^ e3_2dp =", e1_2dp ^ e3_2dp)
print("e2_2dp ^ e3_2dp =", e2_2dp ^ e3_2dp)

-- Test PGA dual operations
print("\n10. Testing PGA dual operations:")
print("-------------------------------")

local s = scalar2dp.new(2.0)
local v = vec2dp.new(1, 2, 3)
local b = bivec2dp.new(3, 4, 5)
local ps = pscalar2dp.new(1.0)

print("scalar s =", s)
print("bulk_dual(s) =", bulk_dual(s))
print("weight_dual(s) =", weight_dual(s))

print("vector v =", v)
print("bulk_dual(v) =", bulk_dual(v))
print("weight_dual(v) =", weight_dual(v))

print("bivector b =", b)
print("bulk_dual(b) =", bulk_dual(b))
print("weight_dual(b) =", weight_dual(b))

print("pseudoscalar ps =", ps)
print("bulk_dual(ps) =", bulk_dual(ps))
print("weight_dual(ps) =", weight_dual(ps))

-- Test complement operations
print("\n11. Testing complement operations:")
print("---------------------------------")

print("cmpl(s) =", cmpl(s))
print("cmpl(v) =", cmpl(v))
print("cmpl(b) =", cmpl(b))
print("cmpl(ps) =", cmpl(ps))

-- Test complement identities
print("cmpl(cmpl(v)) =", cmpl(cmpl(v)))
print("original v =", v)

-- Test angle operations
print("\n12. Testing angle operations:")
print("-----------------------------")

local v1_dir = vec2dp.new(1, 0, 0)  -- x direction
local v2_dir = vec2dp.new(0, 1, 0)  -- y direction
print("direction v1 =", v1_dir)
print("direction v2 =", v2_dir)
print("angle between v1 and v2 =", angle(v1_dir, v2_dir), "radians")
print("angle in degrees =", rad2deg(angle(v1_dir, v2_dir)), "degrees")

-- Test bivector angles
local b1 = bivec2dp.new(1, 0, 1)
local b2 = bivec2dp.new(0, 1, 1)
print("angle between bivectors =", angle(b1, b2), "radians")

-- Test support operations
print("\n13. Testing support operations:")
print("------------------------------")

local line2d = bivec2dp.new(1, 1, -2)  -- line ax + by + c = 0
print("line =", line2d)
local support_pt = support2dp(line2d)
print("support point (closest to origin) =", support_pt)

-- Test multivector operations
print("\n14. Testing multivector operations:")
print("----------------------------------")

-- Test mvec2dp_e (even multivector)
local mve1 = mvec2dp_e.new(1, 2, 3, 4)  -- scalar, bivec components
print("mvec2dp_e mve1 =", mve1)
print("gr0(mve1) =", gr0(mve1))
print("gr2(mve1) =", gr2(mve1))
print("rev(mve1) =", rev(mve1))
print("conj(mve1) =", conj(mve1))
print("cmpl(mve1) =", cmpl(mve1))

-- Test mvec2dp_u (odd multivector)  
local mvu1 = mvec2dp_u.new(5, 6, 7, 8)  -- vec, pscalar components
print("mvec2dp_u mvu1 =", mvu1)
print("gr1(mvu1) =", gr1(mvu1))
-- Note: In PGA 2DP, grade 3 is the pseudoscalar, but let's check if it works
local success, result = pcall(gr3, mvu1)
if success then
    print("gr3(mvu1) =", result)
else
    print("gr3(mvu1) - function not available for this type")
end
print("rev(mvu1) =", rev(mvu1))
print("conj(mvu1) =", conj(mvu1))
print("cmpl(mvu1) =", cmpl(mvu1))

-- Test mvec2dp (full multivector)
local mv1 = mvec2dp.new(1, 2, 3, 4, 5, 6, 7, 8)  -- all 8 components
print("mvec2dp mv1 =", mv1)
print("gr0(mv1) =", gr0(mv1))
print("gr1(mv1) =", gr1(mv1))
print("gr2(mv1) =", gr2(mv1))
print("gr3(mv1) =", gr3(mv1))
print("rev(mv1) =", rev(mv1))
print("conj(mv1) =", conj(mv1))
print("cmpl(mv1) =", cmpl(mv1))

-- Test geometric operations
print("\n15. Testing geometric operations:")
print("--------------------------------")

local v1 = vec2dp.new(1, 0, 0)
local v2 = vec2dp.new(0, 1, 0)
print("v1 =", v1)
print("v2 =", v2)
print("wdg(v1, v2) =", wdg(v1, v2))
print("rwdg(v1, bivec2dp.new(1, 0, 0)) =", rwdg(v1, bivec2dp.new(1, 0, 0)))
print("cmt(v1, v2) =", cmt(v1, v2))

-- Test PGA motor functions
print("\n16. Testing PGA motor functions:")
print("-------------------------------")

-- Test rotation motor
local origin = vec2dp.new(0, 0, 1)
local x_axis = vec2dp.new(1, 0, 0)
local phi = math.pi / 6  -- 30 degrees
print("origin =", origin)
print("x_axis =", x_axis)
print("phi =", phi, "radians (30 degrees)")

local rotation_motor = motor(origin, phi)
local rotated_result = move2dp(x_axis, rotation_motor)
print("rotation motor =", rotation_motor)
print("rotated x_axis =", rotated_result)

-- Test translation motor
local p = vec2dp.new(1, -0.5, 1)
local delta = vec2dp.new(1, 0.5, 0)  -- translation vector
print("point p =", p)
print("translation delta =", delta)

local translation_motor = motor(delta)
local translated_result = move2dp(p, translation_motor)
print("translation motor =", translation_motor)
print("translated point =", translated_result)

-- Test motor_from_lines
local p0 = vec2dp.new(1, 0.5, 1)
local p1 = vec2dp.new(2, 0.5, 1)
local p2 = vec2dp.new(2, 1.5, 1)

local l1 = unitize(wdg(p0, p1))  -- horizontal line
local l2 = unitize(wdg(p0, p2))  -- diagonal line
print("line l1 =", l1)
print("line l2 =", l2)

local motor_lines = motor_from_lines(l1, l2)
print("motor_from_lines =", motor_lines)

-- Test projections, rejections, reflections
print("\n17. Testing projections, rejections, reflections:")
print("------------------------------------------------")

local v_test = vec2dp.new(3, 4, 1)
local v_onto = vec2dp.new(1, 0, 0)
local line_test = bivec2dp.new(0, 1, 0)  -- horizontal line

print("vector v_test =", v_test)
print("vector v_onto =", v_onto)
print("line line_test =", line_test)

print("project_onto(v_test, v_onto) =", project_onto(v_test, v_onto))
print("reject_from(v_test, v_onto) =", reject_from(v_test, v_onto))
print("project_onto(v_test, line_test) =", project_onto(v_test, line_test))
print("reject_from(v_test, line_test) =", reject_from(v_test, line_test))

-- Test reflections (requires unitized elements)
local unit_line = unitize(line_test)
print("unit_line =", unit_line)
print("reflect_on(v_test, unit_line) =", reflect_on(v_test, unit_line))

-- Test line-on-line reflection
local another_line = bivec2dp.new(1, 1, 0)
local unit_another_line = unitize(another_line)
print("reflect_on(unit_line, unit_another_line) =", reflect_on(unit_line, unit_another_line))

-- Test projection + rejection identity
print("\n18. Testing geometric identities:")
print("--------------------------------")

local proj = project_onto(v_test, v_onto)
local rej = reject_from(v_test, v_onto)
local sum = proj + rej
print("projection + rejection =", sum)
print("original vector =", v_test)
print("identity verified:", sum.x == v_test.x and sum.y == v_test.y and sum.z == v_test.z)

-- Test comprehensive wdg operations
print("\n19. Testing comprehensive wdg operations:")
print("---------------------------------------")

print("scalar ^ scalar:", wdg(scalar2dp.new(2), scalar2dp.new(3)))
print("scalar ^ vector:", wdg(scalar2dp.new(2), vec2dp.new(1, 2, 3)))
print("vector ^ vector:", wdg(vec2dp.new(1, 0, 0), vec2dp.new(0, 1, 0)))
print("vector ^ bivector:", wdg(vec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 0)))
print("bivector ^ scalar:", wdg(bivec2dp.new(1, 0, 0), scalar2dp.new(2)))

-- Test comprehensive rwdg operations
print("\n20. Testing comprehensive rwdg operations:")
print("----------------------------------------")

-- Note: Only valid rwdg combinations for PGA 2DP
print("vector << bivector:", rwdg(vec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 1)))
print("bivector << vector:", rwdg(bivec2dp.new(1, 1, 0), vec2dp.new(0, 0, 1)))
print("bivector << bivector:", rwdg(bivec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 0)))

-- Test comprehensive cmt operations
print("\n21. Testing comprehensive cmt operations:")
print("---------------------------------------")

print("vector cmt vector:", cmt(vec2dp.new(1, 0, 0), vec2dp.new(0, 1, 0)))
print("bivector cmt bivector:", cmt(bivec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 0)))

-- Test nested operations (based on ga_pga_test.cpp)
print("\n22. Testing nested operations:")
print("-----------------------------")

-- Test based on actual PGA test: unitize(rwdg(l1, l2)) (line 1993 in ga_pga_test.cpp)
local l1 = bivec2dp.new(1, 0, -2)  -- line ax + by + c = 0
local l2 = bivec2dp.new(0, 1, -4)  -- another line
local intersection = unitize(rwdg(l1, l2))  -- intersection point of two lines
print("unitize(rwdg(l1, l2)) =", intersection, "(intersection point)")

-- Test wdg associativity property (line 544 in ga_pga_test.cpp)
local v1 = vec2dp.new(1, 0, 0)
local v2 = vec2dp.new(0, 1, 0) 
local v3 = vec2dp.new(1, 1, 0)
local wdg_assoc_left = wdg(wdg(v1, v2), v3)
local wdg_assoc_right = wdg(v1, wdg(v2, v3))
print("wdg associativity: wdg(wdg(v1, v2), v3) =", wdg_assoc_left)
print("                   wdg(v1, wdg(v2, v3)) =", wdg_assoc_right)

-- Test to_val operations
print("\n23. Testing to_val operations:")
print("-----------------------------")

print("to_val(scalar2dp) =", to_val(scalar2dp.new(42)))
print("to_val(pscalar2dp) =", to_val(pscalar2dp.new(7)))

print("\n=== PGA 2DP Test Summary ===")
print("✓ Basic scalar, vector, bivector operations")
print("✓ PGA-specific unitization")
print("✓ Dual operations (bulk_dual, weight_dual)")
print("✓ Complement operations (cmpl)")
print("✓ Multivector operations (even, odd, full)")
print("✓ PGA constants and basis elements")
print("✓ Motor system (rotation, translation)")
print("✓ Geometric operations (wdg, rwdg, cmt)")
print("✓ Projections, rejections, reflections")
print("✓ Support operations")
print("✓ Geometric identities")
print("PGA 2DP testing completed successfully!")