-- test_pga3dp.lua - Comprehensive PGA 3DP testing module
-- Tests all Projective Geometric Algebra 3DP (3D + 1 projective) functionality

print("=======================================================")
print("Testing PGA 3DP (Projective Geometric Algebra 3D+1P)")
print("=======================================================")

-- Test precision constant
local eps = eps or 1e-12

-- Test scalar3dp basic operations
print("\n1. Testing scalar3dp basic operations:")
print("-------------------------------------")

local s0 = scalar3dp.new()
local s1 = scalar3dp.new(7.0)
local s2 = scalar3dp.new(2.0)
print("s0:", s0, "empty ctor")
print("s1:", s1, "component ctor")
print("s2:", s2, "component ctor")
print("s1 + s2:", s1 + s2, "addition")
print("s1 * 3:", s1 * 3.0, "scalar multiplication")

-- Test vec3dp basic operations
print("\n2. Testing vec3dp basic operations:")
print("----------------------------------")

local v0 = vec3dp.new()
local v1 = vec3dp.new(1, 2, 3, 0) -- direction vector (w=0)
local v2 = vec3dp.new(3, 4, 5, 1) -- point (w=1)
local v3 = vec3dp.new(-v2)
local v4 = vec3dp.new(4, 6, 8, 1) -- v1 + v2

print("v0:", v0, "empty ctor")
print("v1:", v1, "direction vector (w=0)")
print("v2:", v2, "point (w=1)")
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

-- Test vec3dp geometric operations
print("\n3. Testing vec3dp geometric operations:")
print("--------------------------------------")

print("dot(v1, v2):", dot(v1, v2), "dot product")
print("nrm_sq(v2): ", nrm_sq(v2), "squared norm")
print("nrm(v2): ", nrm(v2), "norm")

print("wdg(v1,v2):", wdg(v1, v2), "wedge product creates line")

-- Test unitization (PGA specific)
print("\n4. Testing unitization (PGA specific):")
print("--------------------------------------")

local v_pt = vec3dp.new(3, 4, 5, 1) -- point vector (can be unitized)
print("point vector v_pt =", v_pt)
print("unitize(v_pt) =", unitize(v_pt))

-- Test bivec3dp operations
print("\n5. Testing bivec3dp operations:")
print("-------------------------------")

local b0 = bivec3dp.new()
local b1 = bivec3dp.new(1, 0, 0, 2, 3, 4) -- line
local b2 = bivec3dp.new(0, 0, 1, 0, 0, 0) -- z-axis line
local line = bivec3dp.new(1, 0, 0, 0, 0, 1) -- 3d line

print("b0:", b0, "empty ctor")
print("b1:", b1, "line")
print("b2:", b2, "z-axis line")
print("line:", line, "3d line")

print("unitize(line) =", unitize(line))

-- Test bivector norms
print("nrm_sq(line): ", nrm_sq(line), "squared norm")
print("nrm(line): ", nrm(line), "norm")

-- Test trivec3dp operations
print("\n6. Testing trivec3dp operations:")
print("-------------------------------")

local tv0 = trivec3dp.new()
local tv1 = trivec3dp.new(1, 2, 3, 4) -- plane
local tv2 = trivec3dp.new(0, 0, 1, 0) -- xy-plane (z = 0)

print("tv0:", tv0, "empty ctor")
print("tv1:", tv1, "plane")
print("tv2:", tv2, "xy-plane")

print("unitize(tv1) =", unitize(tv1))

-- Test trivector norms
print("nrm_sq(tv1): ", nrm_sq(tv1), "squared norm")
print("nrm(tv1): ", nrm(tv1), "norm")

-- Test pscalar3dp
print("\n7. Testing pscalar3dp:")
print("---------------------")

local ps1 = pscalar3dp.new()
local ps2 = pscalar3dp.new(1.0)
print("ps1:", ps1)
print("ps2:", ps2)

-- Test dualnum3dp
print("\n8. Testing dualnum3dp:")
print("---------------------")

local dn1 = dualnum3dp.new(2.0, 1.5)
local dn2 = dualnum3dp.new(4.0, 3.0)
print("dn1:", dn1, "dual number")
print("dn2:", dn2, "dual number")
print("dn1 + dn2:", dn1 + dn2, "addition")
print("dn1 * 2:", dn1 * 2.0, "scalar multiplication")
print("unitize(dn2) =", unitize(dn2))

-- Test PGA 3DP constants
print("\n9. Testing PGA 3DP constants:")
print("-----------------------------")

print("e1_3dp =", e1_3dp)
print("e2_3dp =", e2_3dp)
print("e3_3dp =", e3_3dp)
print("e4_3dp =", e4_3dp)
print("origin_3dp =", origin_3dp)
print("horizon_3dp =", horizon_3dp)
print("x_axis_3dp =", x_axis_3dp)
print("xy_plane_3dp =", xy_plane_3dp)
print("I_3dp =", I_3dp)

-- Test wedge products with constants
print("\n10. Testing wedge products with constants:")
print("-----------------------------------------")

print("e1_3dp ^ e2_3dp =", e1_3dp ^ e2_3dp)
print("e1_3dp ^ e3_3dp =", e1_3dp ^ e3_3dp)
print("e1_3dp ^ e4_3dp =", e1_3dp ^ e4_3dp)
print("e2_3dp ^ e3_3dp =", e2_3dp ^ e3_3dp)

-- Test PGA dual operations
print("\n11. Testing PGA dual operations:")
print("-------------------------------")

local s = scalar3dp.new(2.0)
local v = vec3dp.new(1, 2, 3, 4)
local b = bivec3dp.new(3, 4, 5, 6, 7, 8)
local t = trivec3dp.new(1, 2, 3, 4)
local ps = pscalar3dp.new(1.0)

print("scalar s =", s)
print("left_bulk_dual(s) =", left_bulk_dual(s))
print("left_weight_dual(s) =", left_weight_dual(s))
print("right_bulk_dual(s) =", right_bulk_dual(s))
print("right_weight_dual(s) =", right_weight_dual(s))

print("vector v =", v)
print("left_bulk_dual(v) =", left_bulk_dual(v))
print("left_weight_dual(v) =", left_weight_dual(v))
print("right_bulk_dual(v) =", right_bulk_dual(v))
print("right_weight_dual(v) =", right_weight_dual(v))

print("bivector b =", b)
print("left_bulk_dual(b) =", left_bulk_dual(b))
print("left_weight_dual(b) =", left_weight_dual(b))
print("right_bulk_dual(b) =", right_bulk_dual(b))
print("right_weight_dual(b) =", right_weight_dual(b))

print("trivector t =", t)
print("left_bulk_dual(t) =", left_bulk_dual(t))
print("left_weight_dual(t) =", left_weight_dual(t))
print("right_bulk_dual(t) =", right_bulk_dual(t))
print("right_weight_dual(t) =", right_weight_dual(t))

print("pseudoscalar ps =", ps)
print("left_bulk_dual(ps) =", left_bulk_dual(ps))
print("left_weight_dual(ps) =", left_weight_dual(ps))
print("right_bulk_dual(ps) =", right_bulk_dual(ps))
print("right_weight_dual(ps) =", right_weight_dual(ps))

-- Test complement operations
print("\n12. Testing complement operations:")
print("---------------------------------")

print("lcmpl(s) =", lcmpl(s))
print("rcmpl(s) =", rcmpl(s))
print("lcmpl(v) =", lcmpl(v))
print("rcmpl(v) =", rcmpl(v))
print("lcmpl(b) =", lcmpl(b))
print("rcmpl(b) =", rcmpl(b))
print("lcmpl(t) =", lcmpl(t))
print("rcmpl(t) =", rcmpl(t))
print("lcmpl(ps) =", lcmpl(ps))
print("rcmpl(ps) =", rcmpl(ps))

-- Test complement identities
print("lcmpl(rcmpl(v)) =", lcmpl(rcmpl(v)))
print("rcmpl(lcmpl(v)) =", rcmpl(lcmpl(v)))
print("original v =", v)

-- Test angle operations
print("\n13. Testing angle operations:")
print("-----------------------------")

local v1_dir = vec3dp.new(1, 0, 0, 0) -- x direction
local v2_dir = vec3dp.new(0, 1, 0, 0) -- y direction
print("direction v1 =", v1_dir)
print("direction v2 =", v2_dir)
print("angle between v1 and v2 =", angle(v1_dir, v2_dir), "radians")
print("angle in degrees =", rad2deg(angle(v1_dir, v2_dir)), "degrees")

-- Test bivector angles
local b1 = bivec3dp.new(1, 0, 0, 0, 0, 1)
local b2 = bivec3dp.new(0, 1, 0, 0, 1, 0)
print("angle between bivectors =", angle(b1, b2), "radians")

-- Test trivector angles
local t1 = trivec3dp.new(1, 0, 0, 1)
local t2 = trivec3dp.new(0, 1, 0, 1)
print("angle between trivectors =", angle(t1, t2), "radians")

-- Test support operations
print("\n14. Testing support operations:")
print("------------------------------")

local line3d = bivec3dp.new(1, 0, 0, 0, 0, 1) -- 3d line
print("line =", line3d)
local support_pt = support3dp(line3d)
print("support point for 3D line =", support_pt)

-- Test multivector operations
print("\n15. Testing multivector operations:")
print("----------------------------------")

-- Test mvec3dp_e (even multivector)
local mve1 = mvec3dp_e.new(1, 2, 3, 4, 5, 6, 7, 8) -- scalar, bivec, pscalar
print("mvec3dp_e mve1 =", mve1)
print("gr0(mve1) =", gr0(mve1))
print("gr2(mve1) =", gr2(mve1))
print("gr4(mve1) =", gr4(mve1))
print("rev(mve1) =", rev(mve1))
print("conj(mve1) =", conj(mve1))
print("lcmpl(mve1) =", lcmpl(mve1))
print("rcmpl(mve1) =", rcmpl(mve1))

-- Test mvec3dp_u (odd multivector)
local mvu1 = mvec3dp_u.new(9, 10, 11, 12, 13, 14, 15, 16) -- vec, trivec
print("mvec3dp_u mvu1 =", mvu1)
print("gr1(mvu1) =", gr1(mvu1))
print("gr3(mvu1) =", gr3(mvu1))
print("rev(mvu1) =", rev(mvu1))
print("conj(mvu1) =", conj(mvu1))
print("lcmpl(mvu1) =", lcmpl(mvu1))
print("rcmpl(mvu1) =", rcmpl(mvu1))

-- Test mvec3dp (full multivector)
local mv1 = mvec3dp.new(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16) -- all 16 components
print("mvec3dp mv1 =", mv1)
print("gr0(mv1) =", gr0(mv1))
print("gr1(mv1) =", gr1(mv1))
print("gr2(mv1) =", gr2(mv1))
print("gr3(mv1) =", gr3(mv1))
print("gr4(mv1) =", gr4(mv1))
print("rev(mv1) =", rev(mv1))
print("conj(mv1) =", conj(mv1))
print("lcmpl(mv1) =", lcmpl(mv1))
print("rcmpl(mv1) =", rcmpl(mv1))

-- Test geometric operations
print("\n16. Testing geometric operations:")
print("--------------------------------")

local v1 = vec3dp.new(1, 0, 0, 0)
local v2 = vec3dp.new(0, 1, 0, 0)
print("v1 =", v1)
print("v2 =", v2)
print("wdg(v1, v2) =", wdg(v1, v2))
print("rwdg(v1, bivec3dp.new(1, 0, 0, 0, 0, 0)) =", rwdg(v1, bivec3dp.new(1, 0, 0, 0, 0, 0)))
print("cmt(v1, v2) =", cmt(v1, v2))

-- Test PGA motor functions
print("\n17. Testing PGA motor functions:")
print("-------------------------------")

-- Test rotation motor
local x_axis = vec3dp.new(1, 0, 0, 0)
local z_axis = bivec3dp.new(0, 0, 1, 0, 0, 0) -- z-axis bivector
local phi = math.pi / 6 -- 30 degrees
print("x_axis =", x_axis)
print("z_axis =", z_axis)
print("phi =", phi, "radians (30 degrees)")

local rotation_motor = motor(z_axis, phi)
local rotated_result = move3dp(x_axis, rotation_motor)
print("rotation motor =", rotation_motor)
print("rotated x_axis =", rotated_result)

-- Test translation motor
local p = vec3dp.new(1, -0.5, 0, 1)
local delta = vec3dp.new(1, 0.5, 0, 0) -- translation vector
print("point p =", p)
print("translation delta =", delta)

local translation_motor = motor(delta)
local translated_result = move3dp(p, translation_motor)
print("translation motor =", translation_motor)
print("translated point =", translated_result)

-- Test motor_from_planes
local p0_3d = vec3dp.new(1, 0.5, 0, 1)
local p1_3d = vec3dp.new(2, 0.5, 0, 1)
local p2_3d = vec3dp.new(2, 1.5, 0, 1)
local p3_3d = vec3dp.new(1, 0.5, 1, 1)

local pl1 = unitize(wdg(wdg(p0_3d, p3_3d), p1_3d)) -- plane
local pl2 = unitize(wdg(wdg(p0_3d, p3_3d), p2_3d)) -- plane
print("plane pl1 =", pl1)
print("plane pl2 =", pl2)

local motor_planes = motor_from_planes(pl1, pl2)
print("motor_from_planes =", motor_planes)

-- Test projections, rejections, reflections
print("\n18. Testing projections, rejections, reflections:")
print("------------------------------------------------")

local v_test = vec3dp.new(2, 3, 4, 1)
local v_onto = vec3dp.new(1, 0, 0, 0)
local line_test = bivec3dp.new(0, 0, 1, 0, 0, 0) -- z-axis line
local plane_test = trivec3dp.new(0, 0, 1, 0) -- xy-plane

print("vector v_test =", v_test)
print("vector v_onto =", v_onto)
print("line line_test =", line_test)
print("plane plane_test =", plane_test)

print("project_onto(v_test, v_onto) =", project_onto(v_test, v_onto))
print("reject_from(v_test, v_onto) =", reject_from(v_test, v_onto))
print("project_onto(v_test, line_test) =", project_onto(v_test, line_test))
print("reject_from(v_test, line_test) =", reject_from(v_test, line_test))
print("project_onto(v_test, plane_test) =", project_onto(v_test, plane_test))
print("reject_from(v_test, plane_test) =", reject_from(v_test, plane_test))

-- Test reflections (requires unitized elements)
local unit_plane = unitize(plane_test)
print("unit_plane =", unit_plane)
print("reflect_on(v_test, unit_plane) =", reflect_on(v_test, unit_plane))

-- Test line-on-plane reflection
local unit_line = unitize(line_test)
print("reflect_on(unit_line, unit_plane) =", reflect_on(unit_line, unit_plane))

-- Test plane-on-plane reflection
local another_plane = trivec3dp.new(1, 0, 0, 0) -- yz-plane
local unit_another_plane = unitize(another_plane)
print("reflect_on(unit_plane, unit_another_plane) =", reflect_on(unit_plane, unit_another_plane))

-- Test projection + rejection identity
print("\n19. Testing geometric identities:")
print("--------------------------------")

local proj = project_onto(v_test, v_onto)
local rej = reject_from(v_test, v_onto)
local sum = proj + rej
print("projection + rejection =", sum)
print("original vector =", v_test)
print("identity verified:", sum.x == v_test.x and sum.y == v_test.y and sum.z == v_test.z and sum.w == v_test.w)

-- Test comprehensive wdg operations
print("\n20. Testing comprehensive wdg operations:")
print("---------------------------------------")

print("scalar ^ scalar:", wdg(scalar3dp.new(2), scalar3dp.new(3)))
print("scalar ^ vector:", wdg(scalar3dp.new(2), vec3dp.new(1, 2, 3, 4)))
print("vector ^ vector:", wdg(vec3dp.new(1, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))
print("vector ^ bivector:", wdg(vec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))
print("bivector ^ vector:", wdg(bivec3dp.new(0, 1, 0, 0, 0, 0), vec3dp.new(0, 0, 1, 0)))
print("bivector ^ bivector:", wdg(bivec3dp.new(1, 0, 0, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))

-- Test comprehensive rwdg operations
print("\n21. Testing comprehensive rwdg operations:")
print("----------------------------------------")

-- Note: Only valid rwdg combinations for PGA 3DP
print("vector << bivector:", rwdg(vec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 1, 0, 0, 0)))
print("bivector << vector:", rwdg(bivec3dp.new(1, 1, 0, 0, 0, 0), vec3dp.new(0, 0, 1, 0)))
print("bivector << trivector:", rwdg(bivec3dp.new(1, 0, 0, 0, 0, 0), trivec3dp.new(0, 1, 0, 0)))
print("trivector << vector:", rwdg(trivec3dp.new(1, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))
print("trivector << bivector:", rwdg(trivec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))

-- Test comprehensive cmt operations
print("\n22. Testing comprehensive cmt operations:")
print("---------------------------------------")

print("vector cmt vector:", cmt(vec3dp.new(1, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))
print("bivector cmt bivector:", cmt(bivec3dp.new(1, 0, 0, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))
print("vector cmt bivector:", cmt(vec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))
print("bivector cmt vector:", cmt(bivec3dp.new(1, 0, 0, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))

-- Test nested operations
print("\n23. Testing nested operations:")
print("-----------------------------")

local p0_nest = vec3dp.new(1, 0, 0, 1)
local p1_nest = vec3dp.new(2, 0, 0, 1)
local p2_nest = vec3dp.new(1, 1, 0, 1)

-- Create simpler nested wdg calls
local line_nest = wdg(p0_nest, p1_nest) -- line between two points
local plane_nest = wdg(line_nest, p2_nest) -- plane through line and point
print("nested wdg result =", plane_nest)

-- Test to_val operations
print("\n24. Testing to_val operations:")
print("-----------------------------")

print("to_val(scalar3dp) =", to_val(scalar3dp.new(42)))
print("to_val(pscalar3dp) =", to_val(pscalar3dp.new(7)))
