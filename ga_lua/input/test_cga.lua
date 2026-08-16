-- test_cga.lua - CGA (cga2dc / cga3dc) testing module
-- Tests the Conformal Geometric Algebra bindings: the null embedding, the round
-- and flat objects as joins, reading the geometry back (center / radius /
-- carrier / container), meets, and the conformal versors.

print("=======================================================")
print("Testing CGA (Conformal Geometric Algebra 2DC / 3DC)")
print("=======================================================")

local eps = 1e-12
-- graded results are userdata (compare with ==); radius_sq / center_nrm and the
-- component accessors return plain Lua numbers
local function approx(a, b) return math.abs(a - b) < 1e-9 end

-- Test the null embedding
print("\n1. Testing the null embedding (cga3dc):")
print("---------------------------------------")

local a = cga.round_point3dc(0.0, 0.0, 0.0, 0.0)
local b = cga.round_point3dc(3.0, 4.0, 0.0, 0.0)
print("a:", a, "origin as a round point")
print("b:", b, "the point (3,4,0)")
print("cga.dot(a, b):", cga.dot(a, b), "= -1/2 * distance^2 = -12.5")
assert(cga.dot(a, b) == scalar3dc.new(-12.5))

-- Test the round objects as joins of points
print("\n2. Testing rounds as joins of points:")
print("-------------------------------------")

local p1 = cga.round_point3dc(1.0, 0.0, 0.0, 0.0)
local p2 = cga.round_point3dc(0.0, 1.0, 0.0, 0.0)
local p3 = cga.round_point3dc(0.0, 0.0, 1.0, 0.0)
local p4 = cga.round_point3dc(-1.0, 0.0, 0.0, 0.0)
print("p1 ^ p2:", cga.wdg(p1, p2), "dipole (point pair)")
print("p1 ^ p2 ^ p3:", cga.wdg(cga.wdg(p1, p2), p3), "circle")
local s = cga.wdg(cga.wdg(cga.wdg(p1, p2), p3), p4)
print("p1 ^ p2 ^ p3 ^ p4:", s, "sphere")
assert(cga.is_congruent(s, cga.sphere3dc(0.0, 0.0, 0.0, 1.0), eps))

-- Test reading the geometry back off an object
print("\n3. Testing center, radius and the properties:")
print("---------------------------------------------")

local sph = cga.unitize(cga.sphere3dc(3.0, 4.0, 0.0, 2.0))
print("sphere:", sph, "center (3,4,0), radius 2")
print("radius_sq:", cga.radius_sq(sph), "= 4")
print("center_nrm:", cga.center_nrm(sph), "= |(3,4,0)| = 5")
print("cen:", cga.cen(sph), "center as a round point")
local circ = cga.circle3dc(1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 1.0)
print("cga.car(circle):", cga.car(circ), "carrier (its plane)")
print("cga.con(circle):", cga.con(circ), "container (its sphere)")
print("cga.att(circle):", cga.att(circ), "attitude")
assert(approx(cga.radius_sq(sph), 4.0))
assert(approx(cga.center_nrm(sph), 5.0))

-- Test the meet --- one chain, no case distinction
print("\n4. Testing meets (line/plane cutting a sphere):")
print("-----------------------------------------------")

local s2 = cga.sphere3dc(0.0, 0.0, 0.0, 2.0)
local dip = cga.rwdg(s2, cga.line3dc(0.0, 1.0, 0.0, 1.0, 0.0, 0.0))
print("sphere ^ line:", dip, "off-center chord")
print("radius_sq:", cga.radius_sq(dip), "= 4 - 1 = 3 (real intersection)")
print("center:", cga.cen(cga.unitize(dip)), "midpoint of the chord")
assert(approx(cga.radius_sq(dip), 3.0))
local missing = cga.rwdg(s2, cga.line3dc(0.0, 3.0, 0.0, 1.0, 0.0, 0.0))
print("no intersection:", cga.radius_sq(missing), "< 0 --- the sign is the answer")
assert(cga.radius_sq(missing) < 0.0)
print("sphere ^ plane:", cga.radius_sq(cga.rwdg(s2, cga.plane3dc(0.0, 0.0, 1.0, 1.0))),
      "circle of r^2 = 3")

-- Test the conformal versors
print("\n5. Testing the conformal versors:")
print("---------------------------------")

local M = cga.get_translation(1.0, 2.0, 3.0)
local moved = cga.unitize(cga.transform(s2, M))
print("translation:", M)
print("moved center:", cga.cen(moved), "radius_sq:", cga.radius_sq(moved), "(unchanged)")
assert(approx(cga.radius_sq(moved), 4.0))
local D = cga.get_dilation(0.0, 0.0, 0.0, 3.0)
print("dilation by 3:", cga.radius_sq(cga.unitize(cga.transform(s2, D))), "= 4 * 9 = 36")
assert(approx(cga.radius_sq(cga.unitize(cga.transform(s2, D))), 36.0))
local R = cga.get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.5)
print("rotation:", R)
print("cga.rexp(cga.rlog(R)) == R:", cga.is_same_transform(cga.rexp(cga.rlog(R)), R, eps))
print("versors double cover:", cga.is_same_transform(R, -R, eps), "(M and -M act alike)")
assert(cga.is_same_transform(cga.rexp(cga.rlog(R)), R, eps))
assert(cga.is_same_transform(R, -R, eps))

-- Test cga2dc (the planar algebra) and the constants
print("\n6. Testing cga2dc and the constants:")
print("------------------------------------")

local c1 = cga.circle2dc(0.0, 0.0, 2.0)
local c2 = cga.circle2dc(3.0, 0.0, 2.0)
print("c1:", c1, "circle at the origin, radius 2")
print("c1 ^ c2:", cga.rwdg(c1, c2), "the two intersection points")
print("radius_sq:", cga.radius_sq(cga.rwdg(c1, c2)), "= 1.75")
assert(approx(cga.radius_sq(cga.rwdg(c1, c2)), 1.75))
print("e1_2dc ^ e2_2dc:", e1_2dc ^ e2_2dc, "basis bivector e12")
print("one_2dc:", one_2dc, "I_2dc:", I_2dc)
print("I_3dc_mv_u:", I_3dc_mv_u, "pseudoscalar as an odd multivector")

print("\n=======================================================")
print("CGA tests completed successfully")
print("=======================================================")
