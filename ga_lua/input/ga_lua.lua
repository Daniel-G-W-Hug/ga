-- lua script applying available C++ functions

-- vec2d
print("\nvec2d related:\n")

local v0 = vec2d.new()
local v1 = vec2d.new(1, 1)
local v2 = vec2d.new(1.5, -2)
local v3 = vec2d.new(-v2)
local v4 = vec2d.new(2.5, -1)

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

print("dot(v1, v2):", dot(v1, v2), "dot product")

print("nrm_sq(v1): ", nrm_sq(v1), "squared norm")
print("nrm(v1): ", nrm(v1), "norm")
assert(math.abs(nrm_sq(v1) - 2) < eps)
assert(math.abs(nrm(v1) - math.sqrt(2)) < eps)

print("normalize(v1): ", normalize(v1), "normalized")
print("math.abs( nrm(normalize(v1))-1.0 ):", math.abs(nrm(normalize(v1)) - 1.0))
assert(math.abs(nrm(normalize(v1)) - 1.0) < eps)

print("inv(v1):", inv(v1), "inverse")
print("math.abs(nrm(v1)*nrm(inv(v1))-1):", math.abs(nrm(v1) * nrm(inv(v1)) - 1))
assert(math.abs(nrm(v1) * nrm(inv(v1)) - 1) < eps)

print("wdg(v1,v2):", wdg(v1, v2))

local ps1 = pscalar2d.new()
local ps2 = pscalar2d.new(3.7)
assert(ps1 == pscalar2d.new())
print("ps1:", ps1)
print("ps2:", ps2)

local phi = angle(v1, v2)
print("phi (rad): ", phi)
print("phi (deg): ", rad2deg(phi))
assert(math.abs(angle(vec2d.new(1, 0), vec2d.new(1, 1)) - math.pi / 4) < eps)

-- mvec2d_e
print("\nmvec2d_e related:\n")

local mve0 = mvec2d_e.new()
local mve1 = mvec2d_e.new(1, 1)
local mve2 = mvec2d_e.new(1.5, -2)
local mve3 = mvec2d_e.new(-mve2)
local mve4 = mvec2d_e.new(2.5, -1)
local mve5 = mvec2d_e.new(scalar2d.new(-5))
local mve6 = mvec2d_e.new(pscalar2d.new(-6))
local mve7 = mvec2d_e.new(scalar2d.new(-7), pscalar2d.new(-8))

print("mve0:", mve0, "empty ctor")
print("mve1:", mve1, "component ctor")
print("mve2:", mve2, "component ctor")
print("mve3:", mve3, "even multivector ctor v3 = -v2", "unary minus")
print("mve4:", mve4, "component ctor v4 = v1 + v3 (by hand)")
print("mve5:", mve5, "ctor mvec2d(scalar2d(-5))")
print("mve6:", mve6, "ctor mvec2d(pscalar2d(-6))")
print("mve7:", mve7, "ctor mvec2d(scalar2d(-7),pscalar2d(-8))")

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
assert(gr0(mve2) == scalar2d.new(1.5))
assert(gr2(mve2) == pscalar2d.new(-2))

print("angle_to_re(mve1) - (rad, deg):", angle_to_re(mve1), rad2deg(angle_to_re(mve1)))
print("angle_to_re(mve2) - (rad, deg):", angle_to_re(mve2), rad2deg(angle_to_re(mve2)))

-- inverses
print("\ninverses:\n")
local v = vec2d.new(13, 5)
local vi = inv(v)
print("v:", v)
print("vi:", vi)
print("v*vi:", v * vi)
assert(math.abs(to_val(gr0(v * vi)) - 1) < eps)
assert(math.abs(to_val(gr2(v * vi)) - 0) < eps)
print()

local vc = mvec2d_e.new(13, 5)
local vci = inv(vc)
print("vc:", vc)
print("vci:", vci)
print("vc*vci:", vc * vci)
assert(math.abs(to_val(gr0(vc * vci)) - 1) < eps)
assert(math.abs(to_val(gr2(vc * vci)) - 0) < eps)
print()

local vm = mvec2d.new(13, -27, 3, 5)
local vmi = inv(vm)
print("vm:", vm)
print("vmi:", vmi)
print("vm*vmi:", vm * vmi)
assert(math.abs(to_val(gr0(vm * vmi)) - 1) < eps)
assert(nrm(gr1(vm * vmi)) < eps)
assert(math.abs(to_val(gr2(vm * vmi)) - 0) < eps)
print()

-- mvec2d
print("\nmvec2d related:\n")

local mv0 = mvec2d.new()
local mv1 = mvec2d.new(1, 1, 1, 1)
local mv2 = mvec2d.new(1.5, -2, 1, 2)
local mv3 = mvec2d.new(-mv2)
local mv4 = mvec2d.new(2.5, -1, 2, 3)
local mv5 = mvec2d.new(scalar2d.new(-5))
local mv6 = mvec2d.new(pscalar2d.new(-6))
local mv7 = mvec2d.new(scalar2d.new(-7), pscalar2d.new(-8))
local mv8 = mvec2d.new(vec2d.new(-2, 1))
local mv9 = mvec2d.new(mvec2d_e.new(2.5, 3))

print("mv0:", mv0, "empty ctor")
print("mv1:", mv1, "component ctor")
print("mv2:", mv2, "component ctor")
print("mv3:", mv3, "multivector ctor v3 = -v2", "unary minus")
print("mv4:", mv4, "component ctor v4 = v1 + v3 (by hand)")
print("mv5:", mv5, "ctor mvec2d(scalar2d(-5))")
print("mv6:", mv6, "ctor mvec2d(pscalar2d(-6))")
print("mv7:", mv7, "ctor mvec2d(scalar2d(-7),pscalar2d(-8))")
print("mv8:", mv8, "ctor mvec2d(vec2d(-2,1))")
print("mv9:", mv9, "ctor mvec2d(mvec2d_e(2.5,3))")

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

print("gr0(mv2)", gr0(mv2))
print("gr1(mv2)", gr1(mv2))
print("gr2(mv2)", gr2(mv2))
assert(gr0(mv2) == scalar2d.new(1.5))
assert(gr1(mv2) == vec2d.new(-2, 1))
assert(gr2(mv2) == pscalar2d.new(2))

print("nrm_sq(mv1), nrm(mv1)", nrm_sq(mv1), nrm(mv1))
print("rev(mv1), conj(mv1)", rev(mv1), conj(mv1))

print("inv(mv2), mv2*inv(mv2)", inv(mv2), mv2 * inv(mv2))

-- handling of pscalar2d
print("I_2d*I_2d:", I_2d * I_2d, "nrm_sq(I_2d):", nrm_sq(I_2d), "nrm(I_2d):", nrm(I_2d))
local I = pscalar2d.new(3)
print("I:", I, "I*inv(I):", I * inv(I))
assert(to_val(I * inv(I)) == 1)

print("I+I:", I + I, "I-I:", I - I, "I/3:", I / 3, "to_val(I):", to_val(I))

-- 2d_ops

-- projections, rejections and reflections
print("project_onto(vec2d.new(3,7),vec2d.new(1,1))", project_onto(vec2d.new(3, 7), vec2d.new(1, 1)))
print("reject_from(vec2d.new(3,7),vec2d.new(1,1))", reject_from(vec2d.new(3, 7), vec2d.new(1, 1)))

print("reflect_on_vec(vec2d.new(-2,2),vec2d.new(0,1))", reflect_on_vec(vec2d.new(-2, 2), vec2d.new(0, 1)))
print("reflect_on_hyp(vec2d.new(-2,2),vec2d.new(0,1))", reflect_on_hyp(vec2d.new(-2, 2), vec2d.new(0, 1)))

-- geometric products
print("mvec2d.new(2,-3,5,7)*mvec2d.new(1,2,-5,12):", mvec2d.new(2, -3, 5, 7) * mvec2d.new(1, 2, -5, 12))
print()
print("dot(vec2d.new(2,-3),vec2d.new(1,2)):", dot(vec2d.new(2, -3), vec2d.new(1, 2)))
print("wdg(vec2d.new(2,-3),vec2d.new(1,2)):", to_val(wdg(vec2d.new(2, -3), vec2d.new(1, 2))))
print("angle(vec2d.new(2,-3),vec2d.new(1,2)) [deg]:", rad2deg(angle(vec2d.new(2, -3), vec2d.new(1, 2))))
print("nrm(vec2d.new(2,-3)):", nrm(vec2d.new(2, -3)))
print("nrm(vec2d.new(1,2)):", nrm(vec2d.new(1, 2)))
print("nrm(vec2d.new(2,-3))*nrm(vec2d.new(1,2))", nrm(vec2d.new(2, -3)) * nrm(vec2d.new(1, 2)))
print(
	"sin(angle)*nrm(vec2d.new(2,-3))*nrm(vec2d.new(1,2)))",
	math.sin(angle(vec2d.new(2, -3), vec2d.new(1, 2))) * nrm(vec2d.new(2, -3)) * nrm(vec2d.new(1, 2))
)
print()
print("vec2d.new(2,-3)*vec2d.new(1,2):", vec2d.new(2, -3) * vec2d.new(1, 2))
print()
print("I_2d*mvec2d.new(2,-3,5,7)", I_2d * mvec2d.new(2, -3, 5, 7))
print("dual(mvec2d.new(2,-3,5,7))", dual(mvec2d.new(2, -3, 5, 7)))
print()
print("I_2d*mvec2d_e.new(2,7)", I_2d * mvec2d_e.new(2, 7))
print("dual(mvec2d_e.new(2,7))", dual(mvec2d_e.new(2, 7)))
print()
print("I_2d*vec2d.new(2,7)", I_2d * vec2d.new(2, 7))
print("dual(vec2d.new(2,7))", dual(vec2d.new(2, 7)))
print()
print("mvec2d.new(2,-3,5,7)*I_2d", mvec2d.new(2, -3, 5, 7) * I_2d)
print("dual(mvec2d.new(2,-3,5,7)*I_2d)", dual(mvec2d.new(2, -3, 5, 7)))
print()
print("mvec2d_e.new(2,7)*I_2d", mvec2d_e.new(2, 7) * I_2d)
print("dual(mvec2d_e.new(2,7))", dual(mvec2d_e.new(2, 7)))
print()
print("vec2d.new(2,7)*I_2d", vec2d.new(2, 7) * I_2d)
print("dual(vec2d.new(2,7))", dual(vec2d.new(2, 7)))
print()
print("dual(scalar2d.new(5))", dual(scalar2d.new(5)))
print("dual(pscalar2d.new(3))", dual(pscalar2d.new(3)))
print()
print("vec2d.new(2,7)*mvec2d.new(1,2,-5,12):", vec2d.new(2, 7) * mvec2d.new(1, 2, -5, 12))
print("mvec2d.new(0,2,7,0)*mvec2d.new(1,2,-5,12):", mvec2d.new(0, 2, 7, 0) * mvec2d.new(1, 2, -5, 12))
print()
print("mvec2d.new(1,2,-5,12)*vec2d.new(2,7):", mvec2d.new(1, 2, -5, 12) * vec2d.new(2, 7))
print("mvec2d.new(1,2,-5,12)*mvec2d.new(0,2,7,0):", mvec2d.new(1, 2, -5, 12) * mvec2d.new(0, 2, 7, 0))
print()
print("mvec2d_e.new(2,7)*mvec2d.new(1,2,-5,12):", mvec2d_e.new(2, 7) * mvec2d.new(1, 2, -5, 12))
print("mvec2d_e.new(2,7)*vec2d.new(2,-5):", mvec2d_e.new(2, 7) * vec2d.new(2, -5))
print("mvec2d.new(1,2,-5,12)*mvec2d_e.new(2,7):", mvec2d.new(1, 2, -5, 12) * mvec2d_e.new(2, 7))
print("vec2d.new(2,-5)*mvec2d_e.new(2,7):", vec2d.new(2, -5) * mvec2d_e.new(2, 7))
print()
print("mvec2d_e.new(1,12)*mvec2d_e.new(2,7):", mvec2d_e.new(1, 12) * mvec2d_e.new(2, 7))
print()
-- 2d rotations

for i = 0, 10 do
	local theta = math.pi / 10.0 * i
	print("theta:", theta, "exp(I_2d,theta):", exp(I_2d, theta))
end
print()

for i = 0, 10 do
	local theta = math.pi / 10.0 * i
	print(
		"theta:",
		theta,
		"rotor(I_2d,theta):",
		rotor(I_2d, theta),
		"rotated:",
		rotate(vec2d.new(1, 0), rotor(I_2d, theta))
	)
end
print()

-- vec3d
print("\nvec3d related:\n")

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

print()
local bv1 = wdg(vec3d.new(0, 1, 0), vec3d.new(0, 0, 1))
print("bv1 = wdg(vec3d.new(0,1,0),vec3d.new(0,0,1)):", bv1, "dual(bv1):", dual(bv1))
print()
local bv2 = wdg(vec3d.new(0, 0, 1), vec3d.new(1, 0, 0))
print("bv2 = wdg(vec3d.new(0,0,1),vec3d.new(1,0,0)):", bv2, "dual(bv2):", dual(bv2))
print()
local bv3 = wdg(vec3d.new(1, 0, 0), vec3d.new(0, 1, 0))
print("bv3 = wdg(vec3d.new(1,0,0),vec3d.new(0,1,0)):", bv3, "dual(bv3):", dual(bv3))
print()

local vc = cross(vec3d.new(1, 0, 0), vec3d.new(0, 1, 0))
print("vc = (1,0,0)x(0,1,0):", vc)
print("(1,0,0)x(0,1,0) = -dual(bv3):", -dual(bv3))
print("vc == -dual(bv3) --> dual definition of Hestenes, Doran, Lasenby")
print("vc ==  dual(bv3) --> dual definition of Macdonald")
print()

local ps1 = pscalar3d.new()
local ps2 = pscalar3d.new(3.7) + pscalar3d.new(1.7)
print("ps1:", ps1)
print("ps2:", ps2)

local phi = angle(v1, v2)
print("phi (rad): ", phi)
print("phi (deg): ", rad2deg(phi))
assert(math.abs(angle(vec3d.new(1, 0, 0), vec3d.new(1, 1, 0)) - math.pi / 4) < eps)

-- mvec3d_e
print("\nmvec3d_e related:\n")

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
print()

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

-- inverses
print("\ninverses:\n")
local v = vec3d.new(13, 5, -4)
local vi = inv(v)
print("v:", v)
print("vi:", vi)
print("v*vi:", v * vi)
assert(math.abs(to_val(gr0(v * vi)) - 1) < eps)
assert(math.abs(nrm(gr2(v * vi)) - 0) < eps)
print()

-- mvec3d
print("\nmvec3d related:\n")

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
print()
print("mv2", mv2)
print("gr0(mv2)", gr0(mv2))
print("gr1(mv2)", gr1(mv2))
print("gr2(mv2)", gr2(mv2))
print("gr3(mv2)", gr3(mv2))
assert(gr0(mv2) == scalar3d.new(1.5))
assert(gr1(mv2) == vec3d.new(-2, 1, 2))
assert(gr2(mv2) == bivec3d.new(0, 0, 0))
assert(gr3(mv2) == pscalar3d.new(0))
print()

print("nrm_sq(mv1), nrm(mv1)", nrm_sq(mv1), nrm(mv1))
print("rev(mv1), conj(mv1)", rev(mv1), conj(mv1))
print()
print("inv(mv2), mv2*inv(mv2)", inv(mv2), mv2 * inv(mv2))
print()

-- handling of pscalar3d
print("I_3d*I_3d:", I_3d * I_3d, "nrm_sq(I_3d):", nrm_sq(I_3d), "nrm(I_3d):", nrm(I_3d))
local I = pscalar3d.new(3)
print("I:", I, "I*inv(I):", I * inv(I))
assert(to_val(I * inv(I)) == 1)

print("I+I:", I + I, "I-I:", I - I, "I/3:", I / 3, "to_val(I):", to_val(I))

-- 3d_ops

-- mixed geometric operations
print()
print("bivec3d.new(0,0,1)>>vec3d.new(1,1,1))", (bivec3d.new(0, 0, 1) >> vec3d.new(1, 1, 1)))
print("vec3d.new(1,1,1)<<bivec3d.new(0,0,1))", (vec3d.new(1, 1, 1) << bivec3d.new(0, 0, 1)))
print()
print("dot(bivec3d.new(1,1,1),bivec3d.new(0,0,1))", dot(bivec3d.new(1, 1, 1), bivec3d.new(0, 0, 1)))
print("cmt(bivec3d.new(1,1,1),bivec3d.new(0,0,1))", cmt(bivec3d.new(1, 1, 1), bivec3d.new(0, 0, 1)))
print()
print("angle(vec3d.new(1,0,1),bivec3d.new(0,0,1))", rad2deg(angle(vec3d.new(1, 0, 1), bivec3d.new(0, 0, 1))))
print("bivec3d.new(0,0,1),angle(vec3d.new(1,0,1))", rad2deg(angle(bivec3d.new(0, 0, 1), vec3d.new(1, 0, 1))))
print()

-- projections, rejections and reflections
print("project_onto(vec3d.new(3,7,5),vec3d.new(1,1,1))", project_onto(vec3d.new(3, 7, 5), vec3d.new(1, 1, 1)))
print("reject_from(vec3d.new(3,7,5),vec3d.new(1,1,1))", reject_from(vec3d.new(3, 7, 5), vec3d.new(1, 1, 1)))
print()
print("reflect_on_vec(vec3d.new(-2,2,0),vec3d.new(0,0,1))", reflect_on_vec(vec3d.new(-2, 2, 0), vec3d.new(0, 0, 1)))
print("reflect_on_hyp(vec3d.new(-2,2,0),vec3d.new(0,0,1))", reflect_on_hyp(vec3d.new(-2, 2, 0), vec3d.new(0, 0, 1)))
print()

print("project_onto(vec3d.new(3,7,5),bivec3d.new(1,0,0)))", project_onto(vec3d.new(3, 7, 5), bivec3d.new(1, 0, 0)))
print("reject_from(vec3d.new(3,7,5),bivec3d.new(1,0,0))", reject_from(vec3d.new(3, 7, 5), bivec3d.new(1, 0, 0)))
print()
print("reflect_on(vec3d.new(3,7,5),bivec3d.new(1,0,0))", reflect_on(vec3d.new(3, 7, 5), bivec3d.new(1, 0, 0)))
print()

-- geometric products
print(
	"mvec3d.new(2,-3,5,7,1,2,3,1)*mvec3d.new(1,2,-5,12,2,3,4,1):",
	mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1) * mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1)
)
print()
print("nrm(vec3d.new(-3,5,7)):", nrm(vec3d.new(-3, 5, 7)))
print("nrm(vec3d.new(2,-5,12)):", nrm(vec3d.new(2, -5, 12)))
print()
print("dot(vec3d.new(-3,5,7),vec3d.new(2,-5,12)):", dot(vec3d.new(-3, 5, 7), vec3d.new(2, -5, 12)))
print("wdg(vec3d.new(-3,5,7),vec3d.new(2,-5,12)):", wdg(vec3d.new(-3, 5, 7), vec3d.new(2, -5, 12)))
print("vec3d.new(-3,5,7)*vec3d.new(2,-5,12):", vec3d.new(-3, 5, 7) * vec3d.new(2, -5, 12))
print()
print("wdg(vec3d.new(-3,5,7),bivec3d.new(2,3,4)):", wdg(vec3d.new(-3, 5, 7), bivec3d.new(2, 3, 4)))
print()
print("wdg(bivec3d.new(2,3,4),vec3d.new(-3,5,7)):", wdg(bivec3d.new(2, 3, 4), vec3d.new(-3, 5, 7)))
print()
print("vec3d.new(-3,5,7)*bivec3d.new(2,3,4)):", vec3d.new(-3, 5, 7) * bivec3d.new(2, 3, 4))
print()
print("bivec3d.new(2,3,4)*vec3d.new(-3,5,7)):", bivec3d.new(2, 3, 4) * vec3d.new(-3, 5, 7))
print()
print("bivec3d.new(2,3,4)*bivec3d.new(-3,5,7)):", bivec3d.new(2, 3, 4) * bivec3d.new(-3, 5, 7))
print()

print("I_3d*mvec3d.new(2,-3,5,7,1,2,3,1)", I_3d * mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1))
print("dual(mvec3d.new(2,-3,5,7,1,2,3,1))", dual(mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1)))
print()
print("I_3d*mvec3d_e.new(2,1,2,3)", I_3d * mvec3d_e.new(2, 1, 2, 3))
print("dual(mvec3d_e.new(2,1,2,3))", dual(mvec3d_e.new(2, 1, 2, 3)))
print()
print("I_3d*mvec3d_u.new(2,1,2,3)", I_3d * mvec3d_u.new(2, 1, 2, 3))
print("dual(mvec3d_u.new(2,1,2,3))", dual(mvec3d_u.new(2, 1, 2, 3)))
print()
print("I_3d*bivec3d.new(-3,5,7)", I_3d * bivec3d.new(-3, 5, 7))
print("dual(bivec3d.new(-3,5,7))", dual(bivec3d.new(-3, 5, 7)))
print()
print("I_3d*vec3d.new(-3,5,7)", I_3d * vec3d.new(-3, 5, 7))
print("dual(vec3d.new(-3,5,7))", dual(vec3d.new(-3, 5, 7)))
print()

print("mvec3d.new(2,-3,5,7,1,2,3,1)*I_3d", mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1) * I_3d)
print("dual(mvec3d.new(2,-3,5,7,1,2,3,1))", dual(mvec3d.new(2, -3, 5, 7, 1, 2, 3, 1)))
print()
print("mvec3d_e.new(2,1,2,3)*I_3d", mvec3d_e.new(2, 1, 2, 3) * I_3d)
print("dual(mvec3d_e.new(2,1,2,3))", dual(mvec3d_e.new(2, 1, 2, 3)))
print()
print("mvec3d_u.new(2,1,2,3)*I_3d", mvec3d_u.new(2, 1, 2, 3) * I_3d)
print("dual(mvec3d_u.new(2,1,2,3))", dual(mvec3d_u.new(2, 1, 2, 3)))
print()
print("bivec3d.new(-3,5,7)*I_3d", bivec3d.new(-3, 5, 7) * I_3d)
print("dual(bivec3d.new(-3,5,7))", dual(bivec3d.new(-3, 5, 7)))
print()
print("vec3d.new(-3,5,7)*I_3d", vec3d.new(-3, 5, 7) * I_3d)
print("dual(vec3d.new(-3,5,7))", dual(vec3d.new(-3, 5, 7)))
print()

print(
	"mvec3d.new(2,0,0,0,1,2,3,0)*mvec3d.new(1,2,-5,12,2,3,4,1):",
	mvec3d.new(2, 0, 0, 0, 1, 2, 3, 0) * mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1)
)
print(
	"mvec3d_e.new(2,1,2,3)*mvec3d.new(1,2,-5,12,2,3,4,1):",
	mvec3d_e.new(2, 1, 2, 3) * mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1)
)
print()
print("mvec3d_e.new(2,1,2,3)*mvec3d_e.new(1,2,3,4):", mvec3d_e.new(2, 1, 2, 3) * mvec3d_e.new(1, 2, 3, 4))
print()
print("mvec3d_u.new(2,1,2,3)*mvec3d_u.new(1,2,3,4):", mvec3d_u.new(2, 1, 2, 3) * mvec3d_u.new(1, 2, 3, 4))
print()
print("mvec3d_e.new(2,1,2,3)*mvec3d_u.new(1,2,3,4):", mvec3d_e.new(2, 1, 2, 3) * mvec3d_u.new(1, 2, 3, 4))
print()
print("mvec3d_u.new(2,1,2,3)*mvec3d_e.new(1,2,3,4):", mvec3d_u.new(2, 1, 2, 3) * mvec3d_e.new(1, 2, 3, 4))
print()
print(
	"mvec3d.new(1,2,-5,12,2,3,4,1)*mvec3d.new(2,0,0,0,1,2,3,0):",
	mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1) * mvec3d.new(2, 0, 0, 0, 1, 2, 3, 0)
)
print(
	"mvec3d.new(1,2,-5,12,2,3,4,1)*mvec3d_e.new(2,1,2,3):",
	mvec3d.new(1, 2, -5, 12, 2, 3, 4, 1) * mvec3d_e.new(2, 1, 2, 3)
)
print()
print("mvec3d_e.new(2,1,2,3)*bivec3d.new(2,3,4):", mvec3d_e.new(2, 1, 2, 3) * bivec3d.new(2, 3, 4))
print()
print("bivec3d.new(2,3,4)*mvec3d_e.new(2,1,2,3):", bivec3d.new(2, 3, 4) * mvec3d_e.new(2, 1, 2, 3))
print()
print("pscalar3d.new(1)*pscalar3d.new(1):", pscalar3d.new(1) * pscalar3d.new(1))
print("pscalar3d.new(1)*rev(pscalar3d.new(1)):", pscalar3d.new(1) * rev(pscalar3d.new(1)))
print()
print()

-- 3d rotations

local bv = bivec3d.new(0, 1, 1)
for i = 0, 10 do
	local theta = math.pi / 10.0 * i
	print("theta:", theta, "exp(bv,theta):", exp(bv, theta))
end
print()

for i = 0, 10 do
	local theta = math.pi / 10.0 * i
	print(
		"theta:",
		theta,
		"rotor(bv,theta):",
		rotor(bv, theta),
		"rotated:",
		rotate(vec3d.new(1, 0, 0), rotor(bv, theta))
	)
end
print()

local u = vec2d.new(4.5, 3)
local v = vec2d.new(2.5, 4)
local v_par = project_onto(v, u)
local v_perp = reject_from(v, u)
print("v =  ", v)
print("u = ", u)
print("v_par =  ", v_par)
print("v_perp = ", v_perp)
print("v = v_par + v_perp ", v_par + v_perp)

local mv = mvec3d.new(0, 0, 0, 0, 1, 2, 3, 0)
print("mv:", mv, "dual(mv):", dual(mv))

-- wegde as operator
local u = vec2d.new(-3, 5)
local v = vec2d.new(2, -5)
local uv = u ^ v
print("wdg(vec2d.new(-3,5),vec2d.new(2,-5)):", wdg(vec2d.new(-3, 5), vec2d.new(2, -5)))
print("u:", u, "v:", v, "u^v:", uv)
print()

local u = vec3d.new(-3, 5, 7)
local v = vec3d.new(2, -5, 12)
local uv = u ^ v
print("wdg(vec3d.new(-3,5,7),vec3d.new(2,-5,12)):", wdg(vec3d.new(-3, 5, 7), vec3d.new(2, -5, 12)))
print("u:", u, "v:", v, "u^v:", uv)
print()

local u = vec3d.new(-3, 5, 7)
local v = bivec3d.new(2, 3, 4)
local uv = u ^ v
print("wdg(vec3d.new(-3,5,7),bivec3d.new(2,3,4)):", wdg(vec3d.new(-3, 5, 7), bivec3d.new(2, 3, 4)))
print("u:", u, "v:", v, "u^v:", uv)
print()

local vu = v ^ u
print("wdg(bivec3d.new(2,3,4),vec3d.new(-3,5,7)):", wdg(bivec3d.new(2, 3, 4), vec3d.new(-3, 5, 7)))
print("u:", u, "v:", v, "v^u:", vu)
print()

-- small test for projection and rejections
local u = vec2d.new(3, 0)
local v = vec2d.new(2, 1)

local v_par = dot(v, u) * inv(u)
local v_perp = (v ^ u) * inv(u)

print("u=", u, "v=", v)
print("dot(v, u)=", dot(v, u), "wdg(v, u)=", v ^ u)
print("inv(u)", inv(u))
print("v_par=", v_par, "v_perp=", v_perp)
print("project_onto(v,u)=", project_onto(v, u), "reject_from(v,u)=", reject_from(v, u))
print()

--print("dot(v,I_2d)*inv(I_2d)=", dot(v, I_2d) * inv(I_2d))
print("(v<<I_2d)*inv(I_2d)=", (v << I_2d) * inv(I_2d))
print()

print("v*I_2d=", v * I_2d, "I_2d*v=", I_2d * v)
--print("dot(v,I_2d)=", dot(v, I_2d), "dot(I_2d,v)=", dot(I_2d, v))
print("(v<<I_2d)=", (v << I_2d), "(I_2d>>v)=", (I_2d >> v))
print("0.5*(v*I_2d-I_2d*v)=", 0.5 * (v * I_2d - I_2d * v))
print("0.5*(I_2d*v-v*I_2d)=", 0.5 * (I_2d * v - v * I_2d))
print("0.5*(v*I_2d+I_2d*v)=", 0.5 * (v * I_2d + I_2d * v))
print("0.5*(I_2d*v+v*I_2d)=", 0.5 * (I_2d * v + v * I_2d))

-- inverses
v = vec2d.new(1, 1)
print("v=", v, "inv(v)=", inv(v), "v*inv(v)=", v * inv(v))

local c = mvec2d_e.new(1, 1)
print("c=", c, "inv(c)=", inv(c), "c*inv(c)=", c * inv(c))

-- tests for join and meet:

local p = vec3d.new(1, 0, 1)
local q = vec3d.new(1, 2, 1)
local l = p ^ q
print("p^q = ", l)

local A = bivec3d.new(1, 0, 0)
local B = bivec3d.new(0, 1, 0)
local AB = A * B
print("AB = ", AB)

print("dot(e23_3d,e23_3d)=", dot(e23_3d, e23_3d))
print("dot(e31_3d,e31_3d)=", dot(e31_3d, e31_3d))
print("dot(e12_3d,e12_3d)=", dot(e12_3d, e12_3d))
print("dot(e23_3d,e12_3d)=", dot(e23_3d, e12_3d))
print("dot(e31_3d,e12_3d)=", dot(e31_3d, e12_3d))
print("dot(e12_3d,e23_3d)=", dot(e12_3d, e23_3d))
print("")
-- print("wdg(e23_3d,e23_3d)=", wdg(e23_3d, e23_3d))
-- print("wdg(e31_3d,e31_3d)=", wdg(e31_3d, e31_3d))
-- print("wdg(e12_3d,e12_3d)=", wdg(e12_3d, e12_3d))
-- print("wdg(e23_3d,e12_3d)=", wdg(e23_3d, e12_3d))
-- print("wdg(e31_3d,e12_3d)=", wdg(e31_3d, e12_3d))
-- print("wdg(e12_3d,e23_3d)=", wdg(e12_3d, e23_3d))

local s = scalar2d.new(-6)
print("s = ", s)

local u = vec2d.new(3, 3)
local v = vec2d.new(-0.2, 2)
local uv = u ^ v
print("nrm(uv)", nrm(uv))
print("u << uv", normalize(u << uv))
print("uv >> u", normalize(uv >> u))
print("v << uv", normalize(v << uv))
print("uv >> v", normalize(uv >> v))
print()

local A = mvec3d.new(1, -2, 3, 1, 1, 2, 3, -1)
local B = bivec3d.new(1, 2.5, -1)
print("A = ", A)
print("B = ", B)
print()
print("B*A                       ", B * A)
print()
print("B<<A                      ", B << A)
print("-0.5*(A*B-B*A)            ", -0.5 * (A * B - B * A))
print("(B^A)                     ", B ^ A)
print("B<<A -0.5*(A*B-B*A) + B^A", (B << A) - 0.5 * (A * B - B * A) + B ^ A)

local s = scalar2d.new(2)
local ps = pscalar2d.new(-3)
print("s   ", s)
print("ps  ", ps)
print("s+ps", s + ps)
print("s-ps", s - ps)
