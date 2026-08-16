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

print("pga.dot(v1, v2):", pga.dot(v1, v2), "dot product")
print("pga.bulk_nrm_sq(v2): ", pga.bulk_nrm_sq(v2), "squared bulk norm")
print("pga.bulk_nrm(v2): ", pga.bulk_nrm(v2), "bulk norm")
print("pga.weight_nrm(v2): ", pga.weight_nrm(v2), "weight norm")
print("pga.geom_nrm(v2): ", pga.geom_nrm(v2), "geometric norm")

print("pga.wdg(v1,v2):", pga.wdg(v1, v2), "the wedge product creates a line")

-- Test unitization (PGA specific)
print("\n4. Testing unitization (PGA specific):")
print("--------------------------------------")

local v_pt = vec3dp.new(3, 4, 5, 1) -- point vector (can be unitized)
print("point vector v_pt =", v_pt)
print("pga.unitize(v_pt) =", pga.unitize(v_pt))

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

print("pga.unitize(line) =", pga.unitize(line))

-- Test bivector norms
print("pga.bulk_nrm_sq(line): ", pga.bulk_nrm_sq(line), "squared bulk norm")
print("pga.bulk_nrm(line): ", pga.bulk_nrm(line), "bulk norm")
print("pga.weight_nrm_sq(line): ", pga.weight_nrm_sq(line), "squared weight norm")
print("pga.weight_nrm(line): ", pga.weight_nrm(line), "weight norm")
print("pga.geom_nrm_sq(line): ", pga.geom_nrm_sq(line), "squared geometric norm")
print("pga.geom_nrm(line): ", pga.geom_nrm(line), "geometric norm")

-- Test trivec3dp operations
print("\n6. Testing trivec3dp operations:")
print("-------------------------------")

local tv0 = trivec3dp.new()
local tv1 = trivec3dp.new(1, 2, 3, 4) -- plane
local tv2 = trivec3dp.new(0, 0, 1, 0) -- xy-plane (z = 0)

print("tv0:", tv0, "empty ctor")
print("tv1:", tv1, "plane")
print("tv2:", tv2, "xy-plane")

print("pga.unitize(tv1) =", pga.unitize(tv1))

-- Test trivector norms
print("pga.bulk_nrm_sq(tv1): ", pga.bulk_nrm_sq(tv1), "squared bulk norm")
print("pga.bulk_nrm(tv1): ", pga.bulk_nrm(tv1), "bulk norm")
print("pga.weight_nrm_sq(tv1): ", pga.weight_nrm_sq(tv1), "squared weight norm")
print("pga.weight_nrm(tv1): ", pga.weight_nrm(tv1), "weight norm")
print("pga.geom_nrm_sq(tv1): ", pga.geom_nrm_sq(tv1), "squared geometric norm")
print("pga.geom_nrm(tv1): ", pga.geom_nrm(tv1), "geometric norm")

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

-- Test basic constructors
local dn1 = dualnum3dp.new(2.0, 1.5)
local dn2 = dualnum3dp.new(4.0, 3.0)
print("dn1:", dn1, "dual number (from two components)")
print("dn2:", dn2, "dual number (from two components)")
print("dn1 + dn2:", dn1 + dn2, "addition")
print("dn1 * 2:", dn1 * 2.0, "scalar multiplication")
print("pga.unitize(dn2) =", pga.unitize(dn2))

-- Test constructors from scalar and pseudoscalar
local s_dn = scalar3dp.new(5.0)
local ps_dn = pscalar3dp.new(7.0)
local dn3 = dualnum3dp.new(s_dn, ps_dn)
print("dn3 from scalar+pseudoscalar:", dn3)

-- Test grade extraction
print("pga.gr0(dn3) =", pga.gr0(dn3), "(scalar component)")
print("pga.gr4(dn3) =", pga.gr4(dn3), "(pseudoscalar component)")

-- Test operations with scalar and pseudoscalar
local dn4 = dn1 + s_dn
print("dn1 + scalar =", dn4)
local dn5 = dn1 + ps_dn
print("dn1 + pseudoscalar =", dn5)
local dn6 = dn1 - s_dn
print("dn1 - scalar =", dn6)
local dn7 = dn1 - ps_dn
print("dn1 - pseudoscalar =", dn7)

-- Test PGA 3DP constants
print("\n9. Testing PGA 3DP constants:")
print("-----------------------------")

print("e1_3dp =", e1_3dp)
print("e2_3dp =", e2_3dp)
print("e3_3dp =", e3_3dp)
print("e4_3dp =", e4_3dp)
print("O_3dp =", O_3dp)
print("H_3dp =", H_3dp)
print("x_axis_3dp =", x_axis_3dp)
print("xy_3dp =", xy_3dp)
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
print("pga.l_bulk_dual(s) =", pga.l_bulk_dual(s))
print("pga.l_weight_dual(s) =", pga.l_weight_dual(s))
print("pga.r_bulk_dual(s) =", pga.r_bulk_dual(s))
print("pga.r_weight_dual(s) =", pga.r_weight_dual(s))

print("vector v =", v)
print("pga.l_bulk_dual(v) =", pga.l_bulk_dual(v))
print("pga.l_weight_dual(v) =", pga.l_weight_dual(v))
print("pga.r_bulk_dual(v) =", pga.r_bulk_dual(v))
print("pga.r_weight_dual(v) =", pga.r_weight_dual(v))

print("bivector b =", b)
print("pga.l_bulk_dual(b) =", pga.l_bulk_dual(b))
print("pga.l_weight_dual(b) =", pga.l_weight_dual(b))
print("pga.r_bulk_dual(b) =", pga.r_bulk_dual(b))
print("pga.r_weight_dual(b) =", pga.r_weight_dual(b))

print("trivector t =", t)
print("pga.l_bulk_dual(t) =", pga.l_bulk_dual(t))
print("pga.l_weight_dual(t) =", pga.l_weight_dual(t))
print("pga.r_bulk_dual(t) =", pga.r_bulk_dual(t))
print("pga.r_weight_dual(t) =", pga.r_weight_dual(t))

print("pseudoscalar ps =", ps)
print("pga.l_bulk_dual(ps) =", pga.l_bulk_dual(ps))
print("pga.l_weight_dual(ps) =", pga.l_weight_dual(ps))
print("pga.r_bulk_dual(ps) =", pga.r_bulk_dual(ps))
print("pga.r_weight_dual(ps) =", pga.r_weight_dual(ps))

-- Test complement operations
print("\n12. Testing complement operations:")
print("---------------------------------")

print("pga.l_cmpl(s) =", pga.l_cmpl(s))
print("pga.r_cmpl(s) =", pga.r_cmpl(s))
print("pga.l_cmpl(v) =", pga.l_cmpl(v))
print("pga.r_cmpl(v) =", pga.r_cmpl(v))
print("pga.l_cmpl(b) =", pga.l_cmpl(b))
print("pga.r_cmpl(b) =", pga.r_cmpl(b))
print("pga.l_cmpl(t) =", pga.l_cmpl(t))
print("pga.r_cmpl(t) =", pga.r_cmpl(t))
print("pga.l_cmpl(ps) =", pga.l_cmpl(ps))
print("pga.r_cmpl(ps) =", pga.r_cmpl(ps))

-- Test complement identities
print("pga.l_cmpl(pga.r_cmpl(v)) =", pga.l_cmpl(pga.r_cmpl(v)))
print("pga.r_cmpl(pga.l_cmpl(v)) =", pga.r_cmpl(pga.l_cmpl(v)))
print("original v =", v)

-- Test angle operations
print("\n13. Testing angle operations:")
print("-----------------------------")

local v1_dir = vec3dp.new(1, 0, 0, 0) -- x direction
local v2_dir = vec3dp.new(0, 1, 0, 0) -- y direction
print("direction v1 =", v1_dir)
print("direction v2 =", v2_dir)
print("angle between v1 and v2 =", pga.angle(v1_dir, v2_dir), "radians")
print("angle in degrees =", rad2deg(pga.angle(v1_dir, v2_dir)), "degrees")

-- Test bivector angles
local b1 = bivec3dp.new(1, 0, 0, 0, 0, 1)
local b2 = bivec3dp.new(0, 1, 0, 0, 1, 0)
print("angle between bivectors =", pga.angle(b1, b2), "radians")

-- Test trivector angles
local t1 = trivec3dp.new(1, 0, 0, 1)
local t2 = trivec3dp.new(0, 1, 0, 1)
print("angle between trivectors =", pga.angle(t1, t2), "radians")

-- Test support operations
print("\n14. Testing support operations:")
print("------------------------------")

local line3d = bivec3dp.new(1, 0, 0, 0, 0, 1) -- 3d line
print("line =", line3d)
local support_pt = pga.sup(line3d)
print("support point for 3D line =", support_pt)

-- Test multivector operations
print("\n15. Testing multivector operations:")
print("----------------------------------")

-- Test mvec3dp_e (even multivector)
local mve1 = mvec3dp_e.new(1, 2, 3, 4, 5, 6, 7, 8) -- scalar, bivec, pscalar
print("mvec3dp_e mve1 =", mve1)
print("pga.gr0(mve1) =", pga.gr0(mve1))
print("pga.gr2(mve1) =", pga.gr2(mve1))
print("pga.gr4(mve1) =", pga.gr4(mve1))
print("pga.rev(mve1) =", pga.rev(mve1))
print("pga.conj(mve1) =", pga.conj(mve1))
print("pga.l_cmpl(mve1) =", pga.l_cmpl(mve1))
print("pga.r_cmpl(mve1) =", pga.r_cmpl(mve1))

-- Test mvec3dp_u (odd multivector)
local mvu1 = mvec3dp_u.new(9, 10, 11, 12, 13, 14, 15, 16) -- vec, trivec
print("mvec3dp_u mvu1 =", mvu1)
print("pga.gr1(mvu1) =", pga.gr1(mvu1))
print("pga.gr3(mvu1) =", pga.gr3(mvu1))
print("pga.rev(mvu1) =", pga.rev(mvu1))
print("pga.conj(mvu1) =", pga.conj(mvu1))
print("pga.l_cmpl(mvu1) =", pga.l_cmpl(mvu1))
print("pga.r_cmpl(mvu1) =", pga.r_cmpl(mvu1))

-- Test mvec3dp (full multivector)
local mv1 = mvec3dp.new(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16) -- all 16 components
print("mvec3dp mv1 =", mv1)
print("pga.gr0(mv1) =", pga.gr0(mv1))
print("pga.gr1(mv1) =", pga.gr1(mv1))
print("pga.gr2(mv1) =", pga.gr2(mv1))
print("pga.gr3(mv1) =", pga.gr3(mv1))
print("pga.gr4(mv1) =", pga.gr4(mv1))
print("pga.rev(mv1) =", pga.rev(mv1))
print("pga.conj(mv1) =", pga.conj(mv1))
print("pga.l_cmpl(mv1) =", pga.l_cmpl(mv1))
print("pga.r_cmpl(mv1) =", pga.r_cmpl(mv1))

-- Test geometric operations
print("\n16. Testing geometric operations:")
print("--------------------------------")

local v1 = vec3dp.new(1, 0, 0, 0)
local v2 = vec3dp.new(0, 1, 0, 0)
print("v1 =", v1)
print("v2 =", v2)
print("pga.wdg(v1, v2) =", pga.wdg(v1, v2))
print("pga.rwdg(v1, bivec3dp.new(1, 0, 0, 0, 0, 0)) =", pga.rwdg(v1, bivec3dp.new(1, 0, 0, 0, 0, 0)))
print("pga.cmt(v1, v2) =", pga.cmt(v1, v2))

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

local rotation_motor = pga.get_motor(z_axis, phi)
local rotated_result = pga.move3dp(x_axis, rotation_motor)
print("rotation motor =", rotation_motor)
print("rotated x_axis =", rotated_result)

-- Test translation motor
local p = vec3dp.new(1, -0.5, 0, 1)
local delta = vec3dp.new(1, 0.5, 0, 0) -- translation vector
print("point p =", p)
print("translation delta =", delta)

local translation_motor = pga.get_motor(delta)
local translated_result = pga.move3dp(p, translation_motor)
print("translation motor =", translation_motor)
print("translated point =", translated_result)

-- Test motor_from_planes
local p0_3d = vec3dp.new(1, 0.5, 0, 1)
local p1_3d = vec3dp.new(2, 0.5, 0, 1)
local p2_3d = vec3dp.new(2, 1.5, 0, 1)
local p3_3d = vec3dp.new(1, 0.5, 1, 1)

local pl1 = pga.unitize(pga.wdg(pga.wdg(p0_3d, p3_3d), p1_3d)) -- plane
local pl2 = pga.unitize(pga.wdg(pga.wdg(p0_3d, p3_3d), p2_3d)) -- plane
print("plane pl1 =", pl1)
print("plane pl2 =", pl2)

local motor_planes = pga.get_motor_from_planes(pl1, pl2)
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

print("pga.project_onto(v_test, v_onto) =", pga.project_onto(v_test, v_onto))
print("pga.reject_from(v_test, v_onto) =", pga.reject_from(v_test, v_onto))
print("pga.project_onto(v_test, line_test) =", pga.project_onto(v_test, line_test))
print("pga.reject_from(v_test, line_test) =", pga.reject_from(v_test, line_test))
print("pga.project_onto(v_test, plane_test) =", pga.project_onto(v_test, plane_test))
print("pga.reject_from(v_test, plane_test) =", pga.reject_from(v_test, plane_test))

-- Test reflections (requires unitized elements)
local unit_plane = pga.unitize(plane_test)
print("unit_plane =", unit_plane)
print("pga.reflect_on(v_test, unit_plane) =", pga.reflect_on(v_test, unit_plane))

-- Test line-on-plane reflection
local unit_line = pga.unitize(line_test)
print("pga.reflect_on(unit_line, unit_plane) =", pga.reflect_on(unit_line, unit_plane))

-- Test plane-on-plane reflection
local another_plane = trivec3dp.new(1, 0, 0, 0) -- yz-plane
local unit_another_plane = pga.unitize(another_plane)
print("pga.reflect_on(unit_plane, unit_another_plane) =", pga.reflect_on(unit_plane, unit_another_plane))

-- Test projection + rejection identity
print("\n19. Testing geometric identities:")
print("--------------------------------")

local proj = pga.project_onto(v_test, v_onto)
local rej = pga.reject_from(v_test, v_onto)
local sum = proj + rej
print("projection + rejection =", sum)
print("original vector =", v_test)
print("identity verified:", sum.x == v_test.x and sum.y == v_test.y and sum.z == v_test.z and sum.w == v_test.w)

-- Test comprehensive wdg operations
print("\n20. Testing comprehensive wdg operations:")
print("---------------------------------------")

print("scalar ^ scalar:", pga.wdg(scalar3dp.new(2), scalar3dp.new(3)))
print("scalar ^ vector:", pga.wdg(scalar3dp.new(2), vec3dp.new(1, 2, 3, 4)))
print("vector ^ vector:", pga.wdg(vec3dp.new(1, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))
print("vector ^ bivector:", pga.wdg(vec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))
print("bivector ^ vector:", pga.wdg(bivec3dp.new(0, 1, 0, 0, 0, 0), vec3dp.new(0, 0, 1, 0)))
print("bivector ^ bivector:", pga.wdg(bivec3dp.new(1, 0, 0, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))

-- Test comprehensive rwdg operations
print("\n21. Testing comprehensive rwdg operations:")
print("----------------------------------------")

-- Note: Only valid rwdg combinations for PGA 3DP
print("vector << bivector:", pga.rwdg(vec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 1, 0, 0, 0)))
print("bivector << vector:", pga.rwdg(bivec3dp.new(1, 1, 0, 0, 0, 0), vec3dp.new(0, 0, 1, 0)))
print("bivector << trivector:", pga.rwdg(bivec3dp.new(1, 0, 0, 0, 0, 0), trivec3dp.new(0, 1, 0, 0)))
print("trivector << vector:", pga.rwdg(trivec3dp.new(1, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))
print("trivector << bivector:", pga.rwdg(trivec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))

-- Test comprehensive cmt operations
print("\n22. Testing comprehensive cmt operations:")
print("---------------------------------------")

print("vector cmt vector:", pga.cmt(vec3dp.new(1, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))
print("bivector cmt bivector:", pga.cmt(bivec3dp.new(1, 0, 0, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))
print("vector cmt bivector:", pga.cmt(vec3dp.new(1, 0, 0, 0), bivec3dp.new(0, 1, 0, 0, 0, 0)))
print("bivector cmt vector:", pga.cmt(bivec3dp.new(1, 0, 0, 0, 0, 0), vec3dp.new(0, 1, 0, 0)))

-- Test nested operations
print("\n23. Testing nested operations:")
print("-----------------------------")

local p0_nest = vec3dp.new(1, 0, 0, 1)
local p1_nest = vec3dp.new(2, 0, 0, 1)
local p2_nest = vec3dp.new(1, 1, 0, 1)

-- Create simpler nested wdg calls
local line_nest = pga.wdg(p0_nest, p1_nest) -- line between two points
local plane_nest = pga.wdg(line_nest, p2_nest) -- plane through line and point
print("nested wdg result =", plane_nest)

-- Test to_val operations
print("\n24. Testing to_val operations:")
print("-----------------------------")

print("pga.to_val(scalar3dp) =", pga.to_val(scalar3dp.new(42)))
print("pga.to_val(pscalar3dp) =", pga.to_val(pscalar3dp.new(7)))
