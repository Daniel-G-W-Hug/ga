-- test_sta.lua - STA (sta4ds) testing module
-- Tests the Space-Time Algebra bindings: G(1,3,0), where the time direction g4
-- squares to +1 and the three spatial directions to -1. The signature carries
-- the physics: the sign of a vector's square is its causal character, and a
-- Lorentz boost is the same rotor sandwich as a rotation, with a time-space
-- generator instead of a purely spatial one.

print("=======================================================")
print("Testing STA (Space-Time Algebra, G(1,3,0))")
print("=======================================================")

local eps = 1e-12
local function approx(a, b) return math.abs(a - b) < 1e-9 end

-- Test causal character
print("\n1. Testing causal character:")
print("----------------------------")

local u_rest = vec4ds.new(0.0, 0.0, 0.0, 1.0) -- observer at rest (4-velocity g4)
local x_axis = vec4ds.new(1.0, 0.0, 0.0, 0.0) -- purely spatial direction
local photon = vec4ds.new(1.0, 0.0, 0.0, 1.0) -- light ray: |dx| = |dt|
print("u_rest:", u_rest, "nrm_sq:", sta.nrm_sq(u_rest), "timelike:", sta.is_timelike(u_rest))
print("x_axis:", x_axis, "nrm_sq:", sta.nrm_sq(x_axis), "spacelike:", sta.is_spacelike(x_axis))
print("photon:", photon, "nrm_sq:", sta.nrm_sq(photon), "lightlike:", sta.is_lightlike(photon))
assert(sta.is_timelike(u_rest) and sta.is_spacelike(x_axis) and sta.is_lightlike(photon))
assert(approx(sta.nrm_sq(u_rest), 1.0) and approx(sta.nrm_sq(x_axis), -1.0))

-- Test the observer split
print("\n2. Testing the spacetime split:")
print("-------------------------------")

local v = vec4ds.new(1.0, 2.0, 3.0, 4.0)
print("event v:", v)
print("sta.time_split(v, u_rest):", sta.time_split(v, u_rest), "the observer's clock")
print("sta.space_split(v, u_rest):", sta.space_split(v, u_rest), "relative 3-vector (bivector)")
assert(sta.time_split(v, u_rest) == scalar4ds.new(4.0))

-- Test boosts
print("\n3. Testing boosts (the rotor sandwich):")
print("---------------------------------------")

local phi = 0.5
local B = sta.get_boost(g14_4ds, phi) -- generator: the time-space plane g14
local u = sta.transform(u_rest, B)
print("boost by rapidity", phi, "in the g14 (t,x) plane")
print("boosted 4-velocity:", u)
print("gamma:", u.w, "= math.cosh(phi) =", math.cosh(phi))
print("beta*gamma:", u.x, "= math.sinh(phi) =", math.sinh(phi))
print("sta.rapidity(u_rest, u):", sta.rapidity(u_rest, u), "recovers phi")
print("nrm_sq stays +1:", sta.nrm_sq(u), "a boost preserves the metric")
assert(approx(u.w, math.cosh(phi)) and approx(u.x, math.sinh(phi)))
assert(approx(sta.rapidity(u_rest, u), phi))
assert(approx(sta.nrm_sq(u), 1.0))

-- collinear boosts add sta.rapidity(velocities would need the relativistic formula)
local B2 = sta.get_boost(g14_4ds, 0.3)
local u12 = sta.transform(u_rest, B * B2)
print("0.5 then 0.3 -> rapidity", sta.rapidity(u_rest, u12), "= 0.8 (additive)")
assert(approx(sta.rapidity(u_rest, u12), 0.8))

-- Test rotations
print("\n4. Testing rotations (same sandwich, spatial generator):")
print("--------------------------------------------------------")

local R = sta.get_rotor(g12_4ds, math.pi / 2.0)
local g1 = vec4ds.new(1.0, 0.0, 0.0, 0.0)
local g2 = vec4ds.new(0.0, 1.0, 0.0, 0.0)
print("rotor for theta = pi/2 in the g12 plane")
print("g1 ->", sta.transform(g1, R))
print("g2 ->", sta.transform(g2, R))
print("i.e. g1 -> cos(theta) g1 + sin(theta) g2 --- the same sense as ega3d's e1 in")
print("the e12 plane. Both STA builders use the POSITIVE half-argument sta.exp(+p/2 B):")
print("ega3d needs an explicit -theta/2 for that sense, STA does not, because the")
print("spatial vectors squaring to -1 already supply the flip")
assert(approx(sta.transform(g1, R).y, 1.0))
assert(sta.is_spacelike(sta.transform(g1, R))) -- a rotation keeps the causal character
print("boost then rotate:", sta.transform(u_rest, R * B), "(rotors compose by product)")

-- Test exp / log / sqrt and the double cover
print("\n5. Testing exp / log / sqrt and the double cover:")
print("-------------------------------------------------")

print("sta.log(B):", sta.log(B), "the generator carries phi/2 = 0.25")
print("sta.exp(sta.log(B)) == B:", sta.is_same_transform(sta.exp(sta.log(B)), B, eps))
local half = sta.sqrt(B)
print("sta.sqrt(B):", half)
print("sta.sqrt(B)^2 == B:", sta.is_same_transform(half * half, B, eps))
print("half the rapidity:", sta.rapidity(u_rest, sta.transform(u_rest, half)), "= 0.25")
assert(sta.is_same_transform(sta.exp(sta.log(B)), B, eps))
assert(sta.is_same_transform(half * half, B, eps))
assert(approx(sta.rapidity(u_rest, sta.transform(u_rest, half)), phi / 2.0))
print("sta.is_same_transform(B, -B):", sta.is_same_transform(B, -B, eps),
      "rotors double cover the transformations")
assert(sta.is_same_transform(B, -B, eps))

-- Test the constants and the wedge
print("\n6. Testing constants and products:")
print("----------------------------------")

print("g1_4ds ^ g2_4ds:", g1_4ds ^ g2_4ds, "= g12_4ds:", g12_4ds)
print("one_4ds:", one_4ds, "I_4ds:", I_4ds)
print("I_4ds_mv:", I_4ds_mv, "pseudoscalar as a multivector")
print("sta.rev(g12_4ds):", sta.rev(g12_4ds), "reversion flips a bivector's sign")
assert(g1_4ds ^ g2_4ds == g12_4ds)
assert(sta.rev(g12_4ds) == -g12_4ds)

print("\n=======================================================")
print("STA tests completed successfully")
print("=======================================================")
