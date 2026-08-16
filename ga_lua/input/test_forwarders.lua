-- test_forwarders.lua - the Lua-prelude forwarders (register_forwarders)
--
-- These 26 functions are defined in Lua, not C++, so nothing type-checks them at
-- build time: a wrong name inside one is a runtime nil call in that function only.
-- No other input script calls them, which is how all 26 came to be broken at once
-- when free functions moved into per-algebra tables (the bodies still called the
-- old global names) while the eight other scripts stayed green.
--
-- So this script CALLS every one of them. It checks results where a closed form is
-- obvious and otherwise just pins that the call returns a value of the right type;
-- the point is coverage of the call path.

print("=======================================================")
print("Testing the Lua-prelude forwarders (register_forwarders)")
print("=======================================================")

local function is_userdata(x) return type(x) == "userdata" end
local function approx(a, b) return math.abs(a - b) < 1e-9 end

local checked = 0
local function check(name, value)
    assert(value ~= nil, name .. " returned nil")
    assert(is_userdata(value), name .. " did not return a GA object")
    checked = checked + 1
    print(string.format("  %-24s -> %s", name, tostring(value)))
end

-- PGA 2dp: contractions and expansions --------------------------------------
print("\n1. pga 2dp contractions / expansions:")
print("-------------------------------------")

local p2 = point2dp.new(1.0, 2.0, 1.0)
local l2 = pga.wdg(point2dp.new(0.0, 0.0, 1.0), point2dp.new(1.0, 0.0, 1.0)) -- x-axis

check("l_bulk_contract2dp", pga.l_bulk_contract2dp(l2, p2))
check("l_weight_contract2dp", pga.l_weight_contract2dp(l2, p2))
check("r_bulk_contract2dp", pga.r_bulk_contract2dp(l2, p2))
check("r_weight_contract2dp", pga.r_weight_contract2dp(l2, p2))
check("l_bulk_expand2dp", pga.l_bulk_expand2dp(p2, l2))
check("l_weight_expand2dp", pga.l_weight_expand2dp(p2, l2))
check("r_bulk_expand2dp", pga.r_bulk_expand2dp(p2, l2))
check("r_weight_expand2dp", pga.r_weight_expand2dp(p2, l2))

-- PGA 3dp: contractions and expansions --------------------------------------
print("\n2. pga 3dp contractions / expansions:")
print("-------------------------------------")

local p3 = point3dp.new(1.0, 2.0, 3.0, 1.0)
local l3 = pga.wdg(point3dp.new(0.0, 0.0, 0.0, 1.0), point3dp.new(1.0, 0.0, 0.0, 1.0))

check("l_bulk_contract3dp", pga.l_bulk_contract3dp(l3, p3))
check("l_weight_contract3dp", pga.l_weight_contract3dp(l3, p3))
check("r_bulk_contract3dp", pga.r_bulk_contract3dp(l3, p3))
check("r_weight_contract3dp", pga.r_weight_contract3dp(l3, p3))
check("l_bulk_expand3dp", pga.l_bulk_expand3dp(p3, l3))
check("l_weight_expand3dp", pga.l_weight_expand3dp(p3, l3))
check("r_bulk_expand3dp", pga.r_bulk_expand3dp(p3, l3))
check("r_weight_expand3dp", pga.r_weight_expand3dp(p3, l3))

-- projections ----------------------------------------------------------------
print("\n3. projections (the _unitize_if_needed path):")
print("---------------------------------------------")

check("ortho_proj2dp", pga.ortho_proj2dp(p2, l2))
check("central_proj2dp", pga.central_proj2dp(p2, l2))
check("ortho_antiproj2dp", pga.ortho_antiproj2dp(p2, l2))
check("ortho_proj3dp", pga.ortho_proj3dp(p3, l3))
check("central_proj3dp", pga.central_proj3dp(p3, l3))
check("ortho_antiproj3dp", pga.ortho_antiproj3dp(p3, l3))

-- the orthogonal projection of (1,2) onto the x-axis is (1,0): a closed form
local proj = pga.unitize(pga.ortho_proj2dp(p2, l2))
print(string.format("  projection of (1,2) onto the x-axis = (%g, %g)", proj.x, proj.y))
assert(approx(proj.x, 1.0) and approx(proj.y, 0.0))

-- distances (both gr(a)+gr(b) branches) --------------------------------------
print("\n4. distances (exercises pga.gr and both branches):")
print("--------------------------------------------------")

-- point-to-point: gr 1 + 1 = 2, so the "else" branch (bulk_nrm of the attitude)
local d2_pp = pga.dist2dp(point2dp.new(0.0, 0.0, 1.0), point2dp.new(3.0, 4.0, 1.0))
print(string.format("  dist2dp(point, point)   = %s -> %g", tostring(d2_pp), d2_pp.c0 / d2_pp.c1))
assert(approx(d2_pp.c0 / d2_pp.c1, 5.0))
checked = checked + 1

-- point-to-line: gr 1 + 2 = 3 in 2dp, so the rwdg branch
local d2_pl = pga.dist2dp(point2dp.new(0.0, 2.0, 1.0), l2)
print(string.format("  dist2dp(point, line)    = %s -> %g", tostring(d2_pl), d2_pl.c0 / d2_pl.c1))
assert(approx(math.abs(d2_pl.c0 / d2_pl.c1), 2.0))

local d3_pp = pga.dist3dp(point3dp.new(0.0, 0.0, 0.0, 1.0), point3dp.new(3.0, 4.0, 0.0, 1.0))
print(string.format("  dist3dp(point, point)   = %s -> %g", tostring(d3_pp), d3_pp.c0 / d3_pp.c1))
assert(approx(d3_pp.c0 / d3_pp.c1, 5.0))
checked = checked + 1

local pl3 = pga.wdg(l3, point3dp.new(0.0, 1.0, 0.0, 1.0)) -- xy-plane
local d3_pp2 = pga.dist3dp(point3dp.new(0.0, 0.0, 2.0, 1.0), pl3)
print(string.format("  dist3dp(point, plane)   = %s -> %g", tostring(d3_pp2), d3_pp2.c0 / d3_pp2.c1))
assert(approx(math.abs(d3_pp2.c0 / d3_pp2.c1), 2.0))

-- STA expansions --------------------------------------------------------------
print("\n5. sta expansions:")
print("------------------")

check("l_expand4ds", sta.l_expand4ds(g1_4ds, g12_4ds))
check("r_expand4ds", sta.r_expand4ds(g1_4ds, g12_4ds))

print(string.format("\n%d forwarder calls exercised", checked))
print("=======================================================")
print("Forwarder tests completed successfully")
print("=======================================================")
