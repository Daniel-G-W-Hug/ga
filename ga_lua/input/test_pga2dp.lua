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
local v1 = vec2dp.new(1, 2, 0) -- direction vector (z=0)
local v2 = vec2dp.new(3, 4, 1) -- point (z=1)
local v3 = vec2dp.new(-v2)
local v4 = vec2dp.new(4, 6, 1) -- v1 + v2

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

print("pga.dot(v1, v2):", pga.dot(v1, v2), "dot product")
print("pga.bulk_nrm_sq(v2): ", pga.bulk_nrm_sq(v2), "squared bulk norm")
print("pga.bulk_nrm(v2): ", pga.bulk_nrm(v2), "bulk norm")
print("pga.weight_nrm(v2): ", pga.weight_nrm(v2), "weight norm")
print("pga.geom_nrm(v2): ", pga.geom_nrm(v2), "geometric norm")

print("pga.wdg(v1,v2):", pga.wdg(v1, v2), "the wedge product creates a line")

-- Test unitization (PGA specific)
print("\n4. Testing unitization (PGA specific):")
print("--------------------------------------")

local v_pt = vec2dp.new(3, 4, 1) -- point vector (can be unitized)
print("point vector v_pt =", v_pt)
print("pga.unitize(v_pt) =", pga.unitize(v_pt))

-- Test bivec2dp operations
print("\n5. Testing bivec2dp operations:")
print("-------------------------------")

local b0 = bivec2dp.new()
local b1 = bivec2dp.new(1, 0, 0) -- line
local b2 = bivec2dp.new(0, 1, 0) -- horizontal line (y = 0)
local line = bivec2dp.new(1, 1, -2) -- line ax + by + c = 0

print("b0:", b0, "empty ctor")
print("b1:", b1, "line")
print("b2:", b2, "horizontal line")
print("line:", line, "line ax + by + c = 0")

print("pga.unitize(line) =", pga.unitize(line))

-- Test bivector norms
print("pga.bulk_nrm_sq(line): ", pga.bulk_nrm_sq(line), "squared bulk norm")
print("pga.bulk_nrm(line): ", pga.bulk_nrm(line), "bulk norm")
print("pga.weight_nrm_sq(line): ", pga.weight_nrm_sq(line), "squared weight norm")
print("pga.weight_nrm(line): ", pga.weight_nrm(line), "weight norm")
print("pga.geom_nrm_sq(line): ", pga.geom_nrm_sq(line), "squared geometric norm")
print("pga.geom_nrm(line): ", pga.geom_nrm(line), "geometric norm")

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

-- Test basic constructors
local dn1 = dualnum2dp.new(1.0, 0.5)
local dn2 = dualnum2dp.new(3.0, 2.0)
print("dn1:", dn1, "dual number (from two components)")
print("dn2:", dn2, "dual number (from two components)")
print("dn1 + dn2:", dn1 + dn2, "addition")
print("dn1 * 2:", dn1 * 2.0, "scalar multiplication")
print("pga.unitize(dn2) =", pga.unitize(dn2))

-- Test constructors from scalar and pseudoscalar
local s_dn = scalar2dp.new(2.0)
local ps_dn = pscalar2dp.new(3.0)
local dn3 = dualnum2dp.new(s_dn, ps_dn)
print("dn3 from scalar+pseudoscalar:", dn3)

-- Test grade extraction
print("pga.gr0(dn3) =", pga.gr0(dn3), "(scalar component)")
print("pga.gr3(dn3) =", pga.gr3(dn3), "(pseudoscalar component)")

-- Test operations with scalar and pseudoscalar
local dn4 = dn1 + s_dn
print("dn1 + scalar =", dn4)
local dn5 = dn1 + ps_dn
print("dn1 + pseudoscalar =", dn5)
local dn6 = dn1 - s_dn
print("dn1 - scalar =", dn6)
local dn7 = dn1 - ps_dn
print("dn1 - pseudoscalar =", dn7)

-- Test PGA constants
print("\n8. Testing PGA 2DP constants:")
print("-----------------------------")

print("O_2dp =", O_2dp)
print("e1_2dp =", e1_2dp)
print("e2_2dp =", e2_2dp)
print("e3_2dp =", e3_2dp)
print("H_2dp =", H_2dp)
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
print("pga.bulk_dual(s) =", pga.bulk_dual(s))
print("pga.weight_dual(s) =", pga.weight_dual(s))

print("vector v =", v)
print("pga.bulk_dual(v) =", pga.bulk_dual(v))
print("pga.weight_dual(v) =", pga.weight_dual(v))

print("bivector b =", b)
print("pga.bulk_dual(b) =", pga.bulk_dual(b))
print("pga.weight_dual(b) =", pga.weight_dual(b))

print("pseudoscalar ps =", ps)
print("pga.bulk_dual(ps) =", pga.bulk_dual(ps))
print("pga.weight_dual(ps) =", pga.weight_dual(ps))

-- Test complement operations
print("\n11. Testing complement operations:")
print("---------------------------------")

print("pga.cmpl(s) =", pga.cmpl(s))
print("pga.cmpl(v) =", pga.cmpl(v))
print("pga.cmpl(b) =", pga.cmpl(b))
print("pga.cmpl(ps) =", pga.cmpl(ps))

-- Test complement identities
print("pga.cmpl(pga.cmpl(v)) =", pga.cmpl(pga.cmpl(v)))
print("original v =", v)

-- Test angle operations
print("\n12. Testing angle operations:")
print("-----------------------------")

local v1_dir = vec2dp.new(1, 0, 0) -- x direction
local v2_dir = vec2dp.new(0, 1, 0) -- y direction
print("direction v1 =", v1_dir)
print("direction v2 =", v2_dir)
print("angle between v1 and v2 =", pga.angle(v1_dir, v2_dir), "radians")
print("angle in degrees =", rad2deg(pga.angle(v1_dir, v2_dir)), "degrees")

-- Test bivector angles
local b1 = bivec2dp.new(1, 0, 1)
local b2 = bivec2dp.new(0, 1, 1)
print("angle between bivectors =", pga.angle(b1, b2), "radians")

-- Test support operations
print("\n13. Testing support operations:")
print("------------------------------")

local line2d = bivec2dp.new(1, 1, -2) -- line ax + by + c = 0
print("line =", line2d)
local support_pt = pga.sup(line2d)
print("support point (closest to origin) =", support_pt)

-- Test multivector operations
print("\n14. Testing multivector operations:")
print("----------------------------------")

-- Test mvec2dp_e (even multivector)
local mve1 = mvec2dp_e.new(1, 2, 3, 4) -- scalar, bivec components
print("mvec2dp_e mve1 =", mve1)
print("pga.gr0(mve1) =", pga.gr0(mve1))
print("pga.gr2(mve1) =", pga.gr2(mve1))
print("pga.rev(mve1) =", pga.rev(mve1))
print("pga.conj(mve1) =", pga.conj(mve1))
print("pga.cmpl(mve1) =", pga.cmpl(mve1))

-- Test mvec2dp_u (odd multivector)
local mvu1 = mvec2dp_u.new(5, 6, 7, 8) -- vec, pscalar components
print("mvec2dp_u mvu1 =", mvu1)
print("pga.gr1(mvu1) =", pga.gr1(mvu1))
-- Note: In PGA 2DP, grade 3 is the pseudoscalar, but let's check if it works
local success, result = pcall(gr3, mvu1)
if success then
	print("pga.gr3(mvu1) =", result)
else
	print("pga.gr3(mvu1) - function not available for this type")
end
print("pga.rev(mvu1) =", pga.rev(mvu1))
print("pga.conj(mvu1) =", pga.conj(mvu1))
print("pga.cmpl(mvu1) =", pga.cmpl(mvu1))

-- Test mvec2dp (full multivector)
local mv1 = mvec2dp.new(1, 2, 3, 4, 5, 6, 7, 8) -- all 8 components
print("mvec2dp mv1 =", mv1)
print("pga.gr0(mv1) =", pga.gr0(mv1))
print("pga.gr1(mv1) =", pga.gr1(mv1))
print("pga.gr2(mv1) =", pga.gr2(mv1))
print("pga.gr3(mv1) =", pga.gr3(mv1))
print("pga.rev(mv1) =", pga.rev(mv1))
print("pga.conj(mv1) =", pga.conj(mv1))
print("pga.cmpl(mv1) =", pga.cmpl(mv1))

-- Test geometric operations
print("\n15. Testing geometric operations:")
print("--------------------------------")

local v1 = vec2dp.new(1, 0, 0)
local v2 = vec2dp.new(0, 1, 0)
print("v1 =", v1)
print("v2 =", v2)
print("pga.wdg(v1, v2) =", pga.wdg(v1, v2))
print("pga.rwdg(v1, bivec2dp.new(1, 0, 0)) =", pga.rwdg(v1, bivec2dp.new(1, 0, 0)))
print("pga.cmt(v1, v2) =", pga.cmt(v1, v2))

-- Test PGA motor functions
print("\n16. Testing PGA motor functions:")
print("-------------------------------")

-- Test rotation motor
local origin = vec2dp.new(0, 0, 1)
local x_axis = vec2dp.new(1, 0, 0)
local phi = math.pi / 6 -- 30 degrees
print("origin =", origin)
print("x_axis =", x_axis)
print("phi =", phi, "radians (30 degrees)")

local rotation_motor = pga.get_motor(origin, phi)
local rotated_result = pga.move2dp(x_axis, rotation_motor)
print("rotation motor =", rotation_motor)
print("rotated x_axis =", rotated_result)

-- Test translation motor
local p = vec2dp.new(1, -0.5, 1)
local delta = vec2dp.new(1, 0.5, 0) -- translation vector
print("point p =", p)
print("translation delta =", delta)

local translation_motor = pga.get_motor(delta)
local translated_result = pga.move2dp(p, translation_motor)
print("translation motor =", translation_motor)
print("translated point =", translated_result)

-- Test motor_from_lines
local p0 = vec2dp.new(1, 0.5, 1)
local p1 = vec2dp.new(2, 0.5, 1)
local p2 = vec2dp.new(2, 1.5, 1)

local l1 = pga.unitize(pga.wdg(p0, p1)) -- horizontal line
local l2 = pga.unitize(pga.wdg(p0, p2)) -- diagonal line
print("line l1 =", l1)
print("line l2 =", l2)

local motor_lines = pga.get_motor_from_lines(l1, l2)
print("motor_from_lines =", motor_lines)

-- Test projections, rejections, reflections
print("\n17. Testing projections, rejections, reflections:")
print("------------------------------------------------")

local v_test = vec2dp.new(3, 4, 1)
local v_onto = vec2dp.new(1, 0, 0)
local line_test = bivec2dp.new(0, 1, 0) -- horizontal line

print("vector v_test =", v_test)
print("vector v_onto =", v_onto)
print("line line_test =", line_test)

print("pga.project_onto(v_test, v_onto) =", pga.project_onto(v_test, v_onto))
print("pga.reject_from(v_test, v_onto) =", pga.reject_from(v_test, v_onto))
print("pga.project_onto(v_test, line_test) =", pga.project_onto(v_test, line_test))
print("pga.reject_from(v_test, line_test) =", pga.reject_from(v_test, line_test))

-- Test reflections (requires unitized elements)
local unit_line = pga.unitize(line_test)
print("unit_line =", unit_line)
print("pga.reflect_on(v_test, unit_line) =", pga.reflect_on(v_test, unit_line))

-- Test line-on-line reflection
local another_line = bivec2dp.new(1, 1, 0)
local unit_another_line = pga.unitize(another_line)
print("pga.reflect_on(unit_line, unit_another_line) =", pga.reflect_on(unit_line, unit_another_line))

-- Test projection + rejection identity
print("\n18. Testing geometric identities:")
print("--------------------------------")

local proj = pga.project_onto(v_test, v_onto)
local rej = pga.reject_from(v_test, v_onto)
local sum = proj + rej
print("projection + rejection =", sum)
print("original vector =", v_test)
print("identity verified:", sum.x == v_test.x and sum.y == v_test.y and sum.z == v_test.z)

-- Test comprehensive wdg operations
print("\n19. Testing comprehensive wdg operations:")
print("---------------------------------------")

print("scalar ^ scalar:", pga.wdg(scalar2dp.new(2), scalar2dp.new(3)))
print("scalar ^ vector:", pga.wdg(scalar2dp.new(2), vec2dp.new(1, 2, 3)))
print("vector ^ vector:", pga.wdg(vec2dp.new(1, 0, 0), vec2dp.new(0, 1, 0)))
print("vector ^ bivector:", pga.wdg(vec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 0)))
print("bivector ^ scalar:", pga.wdg(bivec2dp.new(1, 0, 0), scalar2dp.new(2)))

-- Test comprehensive rwdg operations
print("\n20. Testing comprehensive rwdg operations:")
print("----------------------------------------")

-- Note: Only valid rwdg combinations for PGA 2DP
print("vector << bivector:", pga.rwdg(vec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 1)))
print("bivector << vector:", pga.rwdg(bivec2dp.new(1, 1, 0), vec2dp.new(0, 0, 1)))
print("bivector << bivector:", pga.rwdg(bivec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 0)))

-- Test comprehensive cmt operations
print("\n21. Testing comprehensive cmt operations:")
print("---------------------------------------")

print("vector cmt vector:", pga.cmt(vec2dp.new(1, 0, 0), vec2dp.new(0, 1, 0)))
print("bivector cmt bivector:", pga.cmt(bivec2dp.new(1, 0, 0), bivec2dp.new(0, 1, 0)))

-- Test nested operations (based on ga_pga_test.cpp)
print("\n22. Testing nested operations:")
print("-----------------------------")

-- Test based on actual PGA test: pga.unitize(pga.rwdg(l1, l2)) (line 1993 in ga_pga_test.cpp)
local l1 = bivec2dp.new(1, 0, -2) -- line ax + by + c = 0
local l2 = bivec2dp.new(0, 1, -4) -- another line
local intersection = pga.unitize(pga.rwdg(l1, l2)) -- intersection point of two lines
print("pga.unitize(pga.rwdg(l1, l2)) =", intersection, "(intersection point)")

-- Test wdg associativity property (line 544 in ga_pga_test.cpp)
local v1 = vec2dp.new(1, 0, 0)
local v2 = vec2dp.new(0, 1, 0)
local v3 = vec2dp.new(1, 1, 0)
local wdg_assoc_left = pga.wdg(pga.wdg(v1, v2), v3)
local wdg_assoc_right = pga.wdg(v1, pga.wdg(v2, v3))
print("wdg associativity: pga.wdg(pga.wdg(v1, v2), v3) =", wdg_assoc_left)
print("                   pga.wdg(v1, pga.wdg(v2, v3)) =", wdg_assoc_right)

-- Test to_val operations
print("\n23. Testing to_val operations:")
print("-----------------------------")

print("pga.to_val(scalar2dp) =", pga.to_val(scalar2dp.new(42)))
print("pga.to_val(pscalar2dp) =", pga.to_val(pscalar2dp.new(7)))
