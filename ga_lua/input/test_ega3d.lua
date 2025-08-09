-- test_ega3d.lua - Comprehensive EGA 3D testing module
-- Tests all Euclidean Geometric Algebra 3D functionality

print("=================================================")
print("Testing EGA 3D (Euclidean Geometric Algebra 3D)")
print("=================================================")

-- Test precision constant
local eps = eps or 1e-12

-- Test vec3d basic operations
print("\n1. Testing vec3d basic operations:")
print("---------------------------------")

local v0 = vec3d.new()
local v1 = vec3d.new(1, 1, -1)
local v2 = vec3d.new(1.5, -2, 2)
local v3 = vec3d.new(-v2)
local v4 = vec3d.new(2.5, -1, 1)

print("v0:", v0, "empty ctor")
print("v1:", v1, "component ctor")
print("v2:", v2, "component ctor")
print("v3:", v3, "vector ctor v3 = -v2", "unary minus")
print("v4:", v4, "component ctor v4 = v1 + v3 (by hand)")

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

-- Test vec3d geometric operations
print("\n2. Testing vec3d geometric operations:")
print("-------------------------------------")

print("dot(v1, v2):", dot(v1, v2), "dot product")

print("nrm_sq(v1): ", nrm_sq(v1), "squared norm")
print("nrm(v1): ", nrm(v1), "norm")
assert(math.abs(nrm_sq(v1) - 3) < eps)
assert(math.abs(nrm(v1) - math.sqrt(3)) < eps)

print("normalize(v1): ", normalize(v1), "normalized")
print("math.abs( nrm(normalize(v1))-1.0 ):", math.abs(nrm(normalize(v1)) - 1.0))
assert(math.abs(nrm(normalize(v1)) - 1.0) < eps)

print("inv(v1):", inv(v1), "inverse")
print("math.abs(nrm(v1)*nrm(inv(v1))-1):", math.abs(nrm(v1) * nrm(inv(v1)) - 1))
assert(math.abs(nrm(v1) * nrm(inv(v1)) - 1) < eps)

print("wdg(v1,v2):", wdg(v1, v2), "dual(wdg(v1,v2)):", dual(wdg(v1, v2)))

-- Test bivector operations
print("\n3. Testing bivector operations:")
print("------------------------------")

local bv1 = wdg(vec3d.new(0, 1, 0), vec3d.new(0, 0, 1))
print("bv1 = wdg(vec3d.new(0,1,0),vec3d.new(0,0,1)):", bv1, "dual(bv1):", dual(bv1))

local bv2 = wdg(vec3d.new(0, 0, 1), vec3d.new(1, 0, 0))
print("bv2 = wdg(vec3d.new(0,0,1),vec3d.new(1,0,0)):", bv2, "dual(bv2):", dual(bv2))

local bv3 = wdg(vec3d.new(1, 0, 0), vec3d.new(0, 1, 0))
print("bv3 = wdg(vec3d.new(1,0,0),vec3d.new(0,1,0)):", bv3, "dual(bv3):", dual(bv3))

-- Test cross product vs dual
local vc = cross(vec3d.new(1, 0, 0), vec3d.new(0, 1, 0))
print("vc = (1,0,0)x(0,1,0):", vc)
print("(1,0,0)x(0,1,0) = -dual(bv3):", -dual(bv3))
print("vc == -dual(bv3) --> dual definition of Hestenes, Doran, Lasenby")
print("vc ==  dual(bv3) --> dual definition of Macdonald")

-- Test pscalar3d
print("\n4. Testing pscalar3d:")
print("--------------------")

local ps1 = pscalar3d.new()
local ps2 = pscalar3d.new(3.7) + pscalar3d.new(1.7)
print("ps1:", ps1)
print("ps2:", ps2)

-- Test angle calculations
print("\n5. Testing angle calculations:")
print("------------------------------")

local phi = angle(v1, v2)
print("phi (rad): ", phi)
print("phi (deg): ", rad2deg(phi))
assert(math.abs(angle(vec3d.new(1, 0, 0), vec3d.new(1, 1, 0)) - math.pi / 4) < eps)

-- Test mvec3d_e (even multivectors)
print("\n6. Testing mvec3d_e (even multivectors):")
print("----------------------------------------")

local mve0 = mvec3d_e.new()
local mve1 = mvec3d_e.new(1, 1, 1, 3)
local mve2 = mvec3d_e.new(1.5, -2, -2, -1)
local mve3 = mvec3d_e.new(-mve2)
local mve4 = mvec3d_e.new(2.5, -1, -1, 2)
local mve5 = mvec3d_e.new(scalar3d.new(-5))
local mve6 = mvec3d_e.new(bivec3d.new(-6, 1, 7))
local mve7 = mvec3d_e.new(scalar3d.new(-7), bivec3d.new(-6, 1, 7))

print("mve0:", mve0, "empty ctor")
print("mve1:", mve1, "component ctor")
print("mve2:", mve2, "component ctor")
print("mve3:", mve3, "even multivector ctor v3 = -v2", "unary minus")
print("mve4:", mve4, "component ctor v4 = v1 + v2 (by hand)")
print("mve5:", mve5, "ctor mvec3d(scalar3d(-5)")
print("mve6:", mve6, "ctor bivec3d.new(-6,1,7)")
print("mve7:", mve7, "ctor mvec3d_e(scalar3d.new(-7),bivec3d.new(-8,-6,1,7))")

assert(mve4 == mve1 + mve2)
assert(mve1 == mve4 - mve2)
if mve4 == mve1 + mve2 and mve1 == mve4 - mve2 then
	print("mve4 == mve1 + mve2", "comparison, addition and subtraction")
end

local mve5 = 3 * mve4
assert(3 * mve4 == mve4 * 3)
assert(mve5 / 3 == mve4)
print("mv5 = 3 * mve4:", 3 * mve4, "scalar multiplication")
print("mve5 / 3:", mve5 / 3, "scalar division")

print("gr0(mve2)", gr0(mve2))
print("gr2(mve2)", gr2(mve2))
assert(gr0(mve2) == scalar3d.new(1.5))
assert(gr2(mve2) == bivec3d.new(-2, -2, -1))

print("angle(bv1, v1) - (rad, deg):", angle(bv1, v1), rad2deg(angle(bv1, v1)))
print("angle(v1, bv1) - (rad, deg):", angle(v1, bv1), rad2deg(angle(v1, bv1)))

-- Test inverses
print("\n7. Testing inverses:")
print("-------------------")

local v = vec3d.new(13, 5, -4)
local vi = inv(v)
print("v:", v)
print("vi:", vi)
print("v*vi:", v * vi)
assert(math.abs(to_val(gr0(v * vi)) - 1) < eps)
assert(math.abs(nrm(gr2(v * vi)) - 0) < eps)

-- Test full mvec3d
print("\n8. Testing mvec3d (full multivectors):")
print("--------------------------------------")

local mv0 = mvec3d.new()
local mv1 = mvec3d.new(1, 1, 1, 1, 2, 2, 2, 2)
local mv2 = mvec3d.new(1.5, -2, 1, 2, 0, 0, 0, 0)
local mv3 = mvec3d.new(-mv2)
local mv4 = mvec3d.new(2.5, -1, 2, 3, 2, 2, 2, 2)
local mv5 = mvec3d.new(scalar3d.new(-5))
local mv6 = mvec3d.new(pscalar3d.new(-6))
local mv7 = mvec3d.new(vec3d.new(-2, 1, 2), pscalar3d.new(-8))
local mv8 = mvec3d.new(vec3d.new(-2, 1, 2))
local mv9 = mvec3d.new(bivec3d.new(2, 2, 2))
local mv10 = mvec3d.new(mvec3d_e.new(2.5, -1, 2, 3))
local mv11 = mvec3d.new(mvec3d_u.new(2.5, -1, 2, 3))

print("mv0 :", mv0, "empty ctor")
print("mv1 :", mv1, "component ctor")
print("mv2 :", mv2, "component ctor")
print("mv3 :", mv3, "multivector ctor v3 = -v2", "unary minus")
print("mv4 :", mv4, "component ctor v4 = v1 + v3 (by hand)")
print("mv5 :", mv5, "ctor mvec3d(scalar3d(-5))")
print("mv6 :", mv6, "ctor mvec3d(pscalar3d(-6))")
print("mv7 :", mv7, "ctor mvec3d(scalar3d(-7),pscalar3d(-8))")
print("mv8 :", mv8, "ctor mvec3d(vec3d(-2,1,2))")
print("mv9 :", mv9, "ctor mvec3d(bivec3d_e(2,2,2))")
print("mv10:", mv10, "ctor mvec3d(mvec3d_e(2.5,-1,2,3))")
print("mv11:", mv11, "ctor mvec3d(mvec3d_u(2.5,-1,2,3))")

assert(mv4 == mv1 + mv2)
assert(mv1 == mv4 - mv2)
if mv4 == mv1 + mv2 and mv1 == mv4 - mv2 then
	print("mv4 == mv1 + mv2", "comparison, addition and subtraction")
end

local mv5 = 3 * mv4
assert(3 * mv4 == mv4 * 3)
assert(mv5 / 3 == mv4)
print("mv5 = 3 * mv4:", 3 * mv4, "scalar multiplication")
print("mv5 / 3:", mv5 / 3, "scalar division")

print("mv2", mv2)
print("gr0(mv2)", gr0(mv2))
print("gr1(mv2)", gr1(mv2))
print("gr2(mv2)", gr2(mv2))
print("gr3(mv2)", gr3(mv2))
assert(gr0(mv2) == scalar3d.new(1.5))
assert(gr1(mv2) == vec3d.new(-2, 1, 2))
assert(gr2(mv2) == bivec3d.new(0, 0, 0))
assert(gr3(mv2) == pscalar3d.new(0))

print("nrm_sq(mv1), nrm(mv1)", nrm_sq(mv1), nrm(mv1))
print("rev(mv1), conj(mv1)", rev(mv1), conj(mv1))
print("inv(mv2), mv2*inv(mv2)", inv(mv2), mv2 * inv(mv2))

-- Test pscalar3d operations
print("\n9. Testing pscalar3d operations:")
print("-------------------------------")

print("I_3d*I_3d:", I_3d * I_3d, "nrm_sq(I_3d):", nrm_sq(I_3d), "nrm(I_3d):", nrm(I_3d))
local I = pscalar3d.new(3)
print("I:", I, "I*inv(I):", I * inv(I))
assert(to_val(I * inv(I)) == 1)

print("I+I:", I + I, "I-I:", I - I, "I/3:", I / 3, "to_val(I):", to_val(I))

-- Test mixed geometric operations
print("\n10. Testing mixed geometric operations:")
print("--------------------------------------")

print("bivec3d.new(0,0,1)>>vec3d.new(1,1,1))", (bivec3d.new(0, 0, 1) >> vec3d.new(1, 1, 1)))
print("vec3d.new(1,1,1)<<bivec3d.new(0,0,1))", (vec3d.new(1, 1, 1) << bivec3d.new(0, 0, 1)))

print("dot(bivec3d.new(1,1,1),bivec3d.new(0,0,1))", dot(bivec3d.new(1, 1, 1), bivec3d.new(0, 0, 1)))
print("cmt(bivec3d.new(1,1,1),bivec3d.new(0,0,1))", cmt(bivec3d.new(1, 1, 1), bivec3d.new(0, 0, 1)))

print("angle(vec3d.new(1,0,1),bivec3d.new(0,0,1))", rad2deg(angle(vec3d.new(1, 0, 1), bivec3d.new(0, 0, 1))))
print("bivec3d.new(0,0,1),angle(vec3d.new(1,0,1))", rad2deg(angle(bivec3d.new(0, 0, 1), vec3d.new(1, 0, 1))))

-- Test projections, rejections and reflections
print("\n11. Testing projections, rejections and reflections:")
print("---------------------------------------------------")

print("project_onto(vec3d.new(3,7,5),vec3d.new(1,1,1))", project_onto(vec3d.new(3, 7, 5), vec3d.new(1, 1, 1)))
print("reject_from(vec3d.new(3,7,5),vec3d.new(1,1,1))", reject_from(vec3d.new(3, 7, 5), vec3d.new(1, 1, 1)))

print("reflect_on_vec(vec3d.new(-2,2,0),vec3d.new(0,0,1))", reflect_on_vec(vec3d.new(-2, 2, 0), vec3d.new(0, 0, 1)))
print("reflect_on(vec3d.new(-2,2,0),vec3d.new(0,0,1))", reflect_on(vec3d.new(-2, 2, 0), vec3d.new(0, 0, 1)))

print("project_onto(vec3d.new(3,7,5),bivec3d.new(1,0,0)))", project_onto(vec3d.new(3, 7, 5), bivec3d.new(1, 0, 0)))
print("reject_from(vec3d.new(3,7,5),bivec3d.new(1,0,0))", reject_from(vec3d.new(3, 7, 5), bivec3d.new(1, 0, 0)))

print("reflect_on(vec3d.new(3,7,5),bivec3d.new(1,0,0))", reflect_on(vec3d.new(3, 7, 5), bivec3d.new(1, 0, 0)))

-- Test geometric products
print("\n12. Testing geometric products:")
print("------------------------------")

print(
	"mvec3d.new(2,-3,5,7,1,2,3,1)*mvec3d.new(1,2,-5,12,2,3,4,1):",
	mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1) * mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1)
)

print("nrm(vec3d.new(-3,5,7)):", nrm(vec3d.new(-3, 5, 7)))
print("nrm(vec3d.new(2,-5,12)):", nrm(vec3d.new(2, -5, 12)))

print("dot(vec3d.new(-3,5,7),vec3d.new(2,-5,12)):", dot(vec3d.new(-3, 5, 7), vec3d.new(2, -5, 12)))
print("wdg(vec3d.new(-3,5,7),vec3d.new(2,-5,12)):", wdg(vec3d.new(-3, 5, 7), vec3d.new(2, -5, 12)))
print("vec3d.new(-3,5,7)*vec3d.new(2,-5,12):", vec3d.new(-3, 5, 7) * vec3d.new(2, -5, 12))

print("wdg(vec3d.new(-3,5,7),bivec3d.new(2,3,4)):", wdg(vec3d.new(-3, 5, 7), bivec3d.new(2, 3, 4)))
print("wdg(bivec3d.new(2,3,4),vec3d.new(-3,5,7)):", wdg(bivec3d.new(2, 3, 4), vec3d.new(-3, 5, 7)))

print("vec3d.new(-3,5,7)*bivec3d.new(2,3,4)):", vec3d.new(-3, 5, 7) * bivec3d.new(2, 3, 4))
print("bivec3d.new(2,3,4)*vec3d.new(-3,5,7)):", bivec3d.new(2, 3, 4) * vec3d.new(-3, 5, 7))
print("bivec3d.new(2,3,4)*bivec3d.new(-3,5,7)):", bivec3d.new(2, 3, 4) * bivec3d.new(-3, 5, 7))

-- Test dual operations
print("\n13. Testing dual operations:")
print("---------------------------")

print("I_3d*mvec3d.new(2,-3,5,7,1,2,3,1)", I_3d * mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1))
print("dual(mvec3d.new(2,-3,5,7,1,2,3,1))", dual(mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1)))

print("I_3d*mvec3d_e.new(2,1,2,3)", I_3d * mvec3d_e.new(2, 1, 2, 3))
print("dual(mvec3d_e.new(2,1,2,3))", dual(mvec3d_e.new(2, 1, 2, 3)))

print("I_3d*mvec3d_u.new(2,1,2,3)", I_3d * mvec3d_u.new(2, 1, 2, 3))
print("dual(mvec3d_u.new(2,1,2,3))", dual(mvec3d_u.new(2, 1, 2, 3)))

print("I_3d*bivec3d.new(-3,5,7)", I_3d * bivec3d.new(-3, 5, 7))
print("dual(bivec3d.new(-3,5,7))", dual(bivec3d.new(-3, 5, 7)))

print("I_3d*vec3d.new(-3,5,7)", I_3d * vec3d.new(-3, 5, 7))
print("dual(vec3d.new(-3,5,7))", dual(vec3d.new(-3, 5, 7)))

-- Test complement operations
print("\n14. Testing complement operations:")
print("---------------------------------")

print("cmpl(scalar3d.new(5))", cmpl(scalar3d.new(5)))
print("cmpl(vec3d.new(3,4,5))", cmpl(vec3d.new(3, 4, 5)))
print("cmpl(bivec3d.new(1,2,3))", cmpl(bivec3d.new(1, 2, 3)))
print("cmpl(pscalar3d.new(7))", cmpl(pscalar3d.new(7)))

-- Test 3D rotations
print("\n15. Testing 3D rotations:")
print("-------------------------")

local bv = bivec3d.new(0, 1, 1)
print("Exponential rotations:")
for i = 0, 10 do
	local theta = math.pi / 10.0 * i
	print("theta:", theta, "exp(bv,theta):", exp(bv, theta))
end

print("\nRotor rotations:")
for i = 0, 10 do
	local theta = math.pi / 10.0 * i
	print(
		"theta:",
		theta,
		"get_rotor(bv,theta):",
		get_rotor(bv, theta),
		"rotated:",
		rotate(vec3d.new(1, 0, 0), get_rotor(bv, theta))
	)
end

-- Test wedge operator
print("\n16. Testing wedge operator:")
print("---------------------------")

local u = vec3d.new(-3, 5, 7)
local v = vec3d.new(2, -5, 12)
local uv = u ^ v
print("wdg(vec3d.new(-3,5,7),vec3d.new(2,-5,12)):", wdg(vec3d.new(-3, 5, 7), vec3d.new(2, -5, 12)))
print("u:", u, "v:", v, "u^v:", uv)

local u = vec3d.new(-3, 5, 7)
local v = bivec3d.new(2, 3, 4)
local uv = u ^ v
print("wdg(vec3d.new(-3,5,7),bivec3d.new(2,3,4)):", wdg(vec3d.new(-3, 5, 7), bivec3d.new(2, 3, 4)))
print("u:", u, "v:", v, "u^v:", uv)

local vu = v ^ u
print("wdg(bivec3d.new(2,3,4),vec3d.new(-3,5,7)):", wdg(bivec3d.new(2, 3, 4), vec3d.new(-3, 5, 7)))
print("u:", u, "v:", v, "v^u:", vu)

-- Test commutator operations
print("\n17. Testing commutator operations:")
print("---------------------------------")

local v1 = vec3d.new(1, 0, 0)
local v2 = vec3d.new(0, 1, 0)
print("v1 =", v1)
print("v2 =", v2)
print("cmt(v1, v2) =", cmt(v1, v2))

local bv = bivec3d.new(1, 0, 0)
local v = vec3d.new(0, 1, 0)
print("bivector bv =", bv)
print("vector v =", v)
print("cmt(bv, v) =", cmt(bv, v))
print("cmt(v, bv) =", cmt(v, bv))

-- Test complex multivector operations
print("\n18. Testing complex multivector operations:")
print("------------------------------------------")

print(
	"mvec3d.new(2,0,0,0,1,2,3,0)*mvec3d.new(1,2,-5,12,2,3,4,1):",
	mvec3d.new(2, 0, 0, 0, 1, 2, 3, 0) * mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1)
)
print(
	"mvec3d_e.new(2,1,2,3)*mvec3d.new(1,2,-5,12,2,3,4,1):",
	mvec3d_e.new(2, 1, 2, 3) * mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1)
)

print("mvec3d_e.new(2,1,2,3)*mvec3d_e.new(1,2,3,4):", mvec3d_e.new(2, 1, 2, 3) * mvec3d_e.new(1, 2, 3, 4))
print("mvec3d_u.new(2,1,2,3)*mvec3d_u.new(1,2,3,4):", mvec3d_u.new(2, 1, 2, 3) * mvec3d_u.new(1, 2, 3, 4))
print("mvec3d_e.new(2,1,2,3)*mvec3d_u.new(1,2,3,4):", mvec3d_e.new(2, 1, 2, 3) * mvec3d_u.new(1, 2, 3, 4))
print("mvec3d_u.new(2,1,2,3)*mvec3d_e.new(1,2,3,4):", mvec3d_u.new(2, 1, 2, 3) * mvec3d_e.new(1, 2, 3, 4))

-- Test contractions and basis elements
print("\n19. Testing contractions and basis elements:")
print("-------------------------------------------")

print("dot(e23_3d,e23_3d)=", dot(e23_3d, e23_3d))
print("dot(e31_3d,e31_3d)=", dot(e31_3d, e31_3d))
print("dot(e12_3d,e12_3d)=", dot(e12_3d, e12_3d))
print("dot(e23_3d,e12_3d)=", dot(e23_3d, e12_3d))
print("dot(e31_3d,e12_3d)=", dot(e31_3d, e12_3d))
print("dot(e12_3d,e23_3d)=", dot(e12_3d, e23_3d))

-- Test join and meet operations
print("\n20. Testing join and meet operations:")
print("------------------------------------")

local p = vec3d.new(1, 0, 1)
local q = vec3d.new(1, 2, 1)
local l = p ^ q
print("p^q = ", l)

local A = bivec3d.new(1, 0, 0)
local B = bivec3d.new(0, 1, 0)
local AB = A * B
print("AB = ", AB)

-- Test contraction patterns
print("\n21. Testing contraction patterns:")
print("--------------------------------")

local u = vec3d.new(3, 3, 0)
local v = vec3d.new(-0.2, 2, 0)
local uv = u ^ v
print("nrm(uv)", nrm(uv))
print("u << uv", normalize(u << uv))
print("uv >> u", normalize(uv >> u))
print("v << uv", normalize(v << uv))
print("uv >> v", normalize(uv >> v))

-- Test mixed product decomposition
print("\n22. Testing mixed product decomposition:")
print("---------------------------------------")

local A = mvec3d.new(1, -2, 3, 1, 1, 2, 3, -1)
local B = bivec3d.new(1, 2.5, -1)
print("A = ", A)
print("B = ", B)
print("B*A                       ", B * A)
print("B<<A                      ", B << A)
print("-0.5*(A*B-B*A)            ", -0.5 * (A * B - B * A))
print("(B^A)                     ", B ^ A)
print("B<<A -0.5*(A*B-B*A) + B^A", (B << A) - 0.5 * (A * B - B * A) + B ^ A)

-- Test pseudoscalar operations
print("\n23. Testing pseudoscalar operations:")
print("----------------------------------")

print("pscalar3d.new(1)*pscalar3d.new(1):", pscalar3d.new(1) * pscalar3d.new(1))
print("pscalar3d.new(1)*rev(pscalar3d.new(1)):", pscalar3d.new(1) * rev(pscalar3d.new(1)))

-- Test dual with multivector
local mv = mvec3d.new(0, 0, 0, 0, 1, 2, 3, 0)
print("mv:", mv, "dual(mv):", dual(mv))

-- cross product and equivalence to outer product (wdg)
print("\n24. Testing cross product and equivalence to wedge product:")
print("----------------------------------")
local a = vec3d.new(0, 0, 1)
local b = vec3d.new(1, 0, 0)
local c = cross(a, b)
local wab = wdg(a, b)

print("a = ", a)
print("b = ", b)
print("c = a x b = ", c)
print("a ^ b = ", wab)
print("dual(a ^ b) = ", dual(wab)) -- here they are simply duals of each other
print("dual(a x b) = ", dual(c))
print()
print("-I_3d * (a ^ b) = ", -I_3d * wab) -- Hestenes form of duality
print("I_3d * (a x b) = ", I_3d * c)

-- cross product and equivalence to outer product (wdg)
print("\n25. Testing basic vector product symmetries:")
print("----------------------------------")
local a = vec3d.new(1, 0, 2)
local b = vec3d.new(2, 0, 1)
local dot_ab = dot(a, b)
local wdg_ab = wdg(a, b)
local ab = a * b
local ba = b * a
local sym_ab = 0.5 * (a * b + b * a)
local asym_ab = 0.5 * (a * b - b * a)

print("a = ", a)
print("b = ", b)
print("dot_ab = ", dot_ab)
print("wdg_ab = ", wdg_ab)
print("a*b = ", ab)
print("b*a = ", ba)
print("sym(a,b)  = 0.5*(a*b + b*a) = ", sym_ab, ", gr0 = ", gr0(sym_ab))
print("asym(a,b) = 0.5*(a*b - b*a) = ", asym_ab, ", gr2 = ", gr2(asym_ab))
