-- test_functions.lua - the free functions bound most recently (2026-08-30)
--
-- moment_about (pga), is_simple / boost_part / rot_part (sta), and
-- reciprocal_frame (ega + sta). These verify the BINDING, which is a different
-- question from the mathematics: the C++ suites already pin the maths (see
-- ga_ega{2,3}d_test.hpp and ga_sta4ds_test.hpp for reciprocal_frame, including
-- the Lorentzian sign structure a Euclidean check cannot see). What can only go
-- wrong at this layer is overload selection, container conversion, and whether a
-- C++ exception reaches Lua as an error instead of as garbage.

print("=======================================================")
print("Testing recently bound free functions")
print("=======================================================")

local function approx(a, b) return math.abs(a - b) < 1e-12 end
local function num(v) return tonumber(tostring(v):match("%(([-%d%.e+]+)%)")) end

print("\n1. pga.moment_about - the torque of a force line about a point")
print("--------------------------------------------------------------")
-- a 3 N force in +y applied at (2, 0): torque about the origin is +6
local P = vec2dp.new(2.0, 0.0, 1.0)  -- point of application (z = 1)
local f = vec2dp.new(0.0, 3.0, 0.0)  -- force direction  (z = 0)
local F = pga.wdg(P, f)              -- the force LINE
local O = vec2dp.new(0.0, 0.0, 1.0)
local M = pga.moment_about(O, F)
print("  F                  = " .. tostring(F))
print("  moment_about(O, F) = " .. tostring(M))
-- r in "bulk(F) - r ^ att(F)" is the DIRECTION O -> R, which vanishes at the
-- origin, so the moment there is exactly bulk(F)
assert(tostring(M) == tostring(pga.bulk(F)), "moment about the origin != bulk(F)")
-- and about the point of application the moment must vanish
assert(tostring(pga.moment_about(P, F)) == tostring(bivec2dp.new(0.0, 0.0, 0.0)),
       "moment about the application point is not zero")
-- the 3dp overload must resolve to its own signature, not the 2dp one
local P3 = vec3dp.new(1.0, 0.0, 0.0, 1.0)
local F3 = pga.wdg(P3, vec3dp.new(0.0, 0.0, 2.0, 0.0))
assert(pcall(function() return pga.moment_about(vec3dp.new(0, 0, 0, 1), F3) end),
       "the 3dp moment_about overload does not resolve")
print("  2dp and 3dp overloads both correct")

print("\n2. sta bivector split - is_simple / boost_part / rot_part")
print("---------------------------------------------------------")
local Brot = bivec4ds.new(0, 0, 0, 1, 0, 0)   -- a spatial plane: pure rotation
local Bboost = bivec4ds.new(1, 0, 0, 0, 0, 0) -- a time-space plane: pure boost
local zero = bivec4ds.new(0, 0, 0, 0, 0, 0)
assert(sta.is_simple(Brot), "a pure rotation bivector is simple")
assert(sta.is_simple(Bboost), "a pure boost bivector is simple")
assert(tostring(sta.rot_part(Brot)) == tostring(Brot), "rot_part of a rotation")
assert(tostring(sta.boost_part(Brot)) == tostring(zero), "boost_part of a rotation")
assert(tostring(sta.boost_part(Bboost)) == tostring(Bboost), "boost_part of a boost")
local Bmix = Brot + Bboost
assert(not sta.is_simple(Bmix), "a mixed bivector is NOT simple")
assert(tostring(sta.boost_part(Bmix) + sta.rot_part(Bmix)) == tostring(Bmix),
       "boost_part + rot_part must reconstruct the bivector")
print("  pure rotation, pure boost and the mixed case all split correctly")

print("\n3. reciprocal_frame - dot(a_i, r_j) == delta_ij")
print("-----------------------------------------------")
local function check_frame(label, frame, recip, dotf)
    assert(#recip == #frame, label .. ": wrong number of reciprocal vectors")
    for i = 1, #frame do
        for j = 1, #recip do
            local want = (i == j) and 1.0 or 0.0
            assert(approx(num(dotf(frame[i], recip[j])), want),
                   string.format("%s: dot(a_%d, r_%d) != %g", label, i, j, want))
        end
    end
    print(string.format("  %-8s %d x %d deltas correct", label, #frame, #recip))
end

-- SKEWED frames on purpose: an orthonormal one is its own reciprocal and would
-- pass even if the function simply returned its argument
local a2 = {vec2d.new(2.0, 0.0), vec2d.new(1.0, 3.0)}
check_frame("ega 2d", a2, ega.reciprocal_frame(a2[1], a2[2]), ega.dot)

local a3 = {vec3d.new(1.0, 0.0, 0.0), vec3d.new(1.0, 2.0, 0.0), vec3d.new(0.0, 1.0, 3.0)}
check_frame("ega 3d", a3, ega.reciprocal_frame(a3[1], a3[2], a3[3]), ega.dot)

local a4 = {vec4ds.new(1, 0, 0, 0), vec4ds.new(0, 1, 0, 0),
            vec4ds.new(0, 0, 1, 0), vec4ds.new(0, 0, 0, 1)}
check_frame("sta 4d", a4, sta.reciprocal_frame(a4[1], a4[2], a4[3], a4[4]), sta.dot)

-- the orthonormal control: here the reciprocal IS the frame
local e = {vec3d.new(1, 0, 0), vec3d.new(0, 1, 0), vec3d.new(0, 0, 1)}
local re = ega.reciprocal_frame(e[1], e[2], e[3])
assert(tostring(re[1]) == tostring(e[1]) and tostring(re[2]) == tostring(e[2])
       and tostring(re[3]) == tostring(e[3]), "orthonormal frame is self-reciprocal")
print("  orthonormal frame is its own reciprocal")

-- a degenerate frame throws std::invalid_argument in C++; that must arrive as a
-- Lua error rather than as infinities. sol prints the exception as it crosses the
-- boundary, so the next two "[sol2] An exception occurred" lines are EXPECTED --
-- they are the check succeeding, not a failure. Announced here so that neither a
-- reader nor a log scan mistakes them for one.
print("  (expect two '[sol2] An exception occurred: ... degenerate frame' lines)")
assert(not pcall(function()
    return ega.reciprocal_frame(vec2d.new(1, 1), vec2d.new(2, 2)) end),
    "a degenerate 2d frame must raise")
assert(not pcall(function()
    return ega.reciprocal_frame(vec3d.new(1, 0, 0), vec3d.new(2, 0, 0),
                                vec3d.new(3, 0, 0)) end),
    "a degenerate 3d frame must raise")
print("  degenerate frames raise a lua error (the C++ throw crosses sol)")

print("\nOK - moment_about, sta bivector split, reciprocal_frame")
