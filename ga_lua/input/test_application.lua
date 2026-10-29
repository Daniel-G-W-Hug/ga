-- test_application.lua - application testing module
-- Tests application functionality

print("=================================================")
print("Testing application specific usage")
print("=================================================")

-- equilibrium of moments
print("\n1. Testing equilibrium of moments (ega2d vs. ega2dp):")
print("-----------------------------")
print("\nega2d:")
local F = vec2d.new(0.5, 0)
local lever1 = vec2d.new(0, 0.5)
local lever2 = vec2d.new(0, -0.5)
local mlever1 = ega.wdg(lever1, F) -- moment1 relative to origin
local mlever2 = ega.wdg(lever2, F) -- moment2 relative to origin
local mres = mlever1 + mlever2
print("F = ", F)
print("lever1 = ", lever1)
print("lever2 = ", lever2)
print("mlever1 = ", mlever1)
print("mlever2 = ", mlever2)
print("mres = mlever1 + mlever2 = ", mres)
print("")

print("\npga2dp:")
local f = vec2dp.new(0.5, 0, 0)
local P1 = vec2dp.new(0, 0.5, 1) -- point P1 the force acts at
local P2 = vec2dp.new(0, -0.5, 1) -- point P2 the force acts at
local m1 = pga.wdg(P1, f) -- resulting moment line from f acting at P1
local m2 = pga.wdg(P2, f) -- resulting moment line from f acting at P2
mres = m1 + m2
print("f = ", f)
print("P1 = ", P1)
print("P2 = ", P2)
print("m1 = ", m1)
print("m2 = ", m2)
print("mres = m1 + m2 = ", mres)
print("pga.att(m1) = ", pga.att(m1))
print("pga.att(m2) = ", pga.att(m2))
print("pga.att(mres) = ", pga.att(mres))
print("pga.sup(m1) = ", pga.unitize(pga.sup(m1)))
print("pga.sup(m2) = ", pga.unitize(pga.sup(m2)))
print("pga.sup(mres) = ", pga.unitize(pga.sup(mres)))

-- a joint's own specification: where it may go (range) and what its actuator can
-- deliver (drive). Both default to the NEUTRAL value -- the unrestricted joint, the ideal
-- actuator -- so a joint that states neither behaves as it always did.
print("\n2. A joint's own specification: range and drive (pga2dp, pga3dp):")
print("-----------------------------")
local r0 = joint_range2dp()
assert(r0.lo == -math.huge and r0.hi == math.huge and pga.is_unrestricted(r0))
local r = joint_range2dp(-math.pi, -0.1745)
assert(r.lo == -math.pi and r.hi == -0.1745 and not pga.is_unrestricted(r))
assert(r.k_stop == 0 and r.c_stop == 0) -- stated, not enforced
local rs = joint_range3dp(-0.5, 0.5, 100.0, 2.0)
assert(rs.k_stop == 100.0 and rs.c_stop == 2.0 and not pga.is_unrestricted(rs))
local d = joint_drive2dp()
assert(d.tau_max == math.huge and d.qd_max == math.huge and d.armature == 0 and d.actuated)
d.tau_max = 25.0
local j = joint_state2dp()
assert(pga.is_unrestricted(j.range) and j.drive.tau_max == math.huge)
j.range = r
j.drive = d
assert(j.range.lo == -math.pi and j.range.hi == -0.1745 and j.drive.tau_max == 25.0)
local j3 = joint_state3dp()
assert(pga.is_unrestricted(j3.range) and joint_drive3dp().tau_max == math.huge)
j3.range = rs
assert(j3.range.k_stop == 100.0)
print("range: ", r)
print("drive: ", d)

-- sizing a stop's damper for an explicit integrator: critical for the link it stops,
-- and the stiffest spring that link carries at the step
local I = 2.65e-4
assert(math.abs(critical_damping(4.0, I) - 2.0 * math.sqrt(4.0 * I)) < 1e-15)
assert(math.abs(max_explicit_stiffness(I, 1e-3) - I * (1.0 / 2e-3) ^ 2) < 1e-9)
assert(math.abs(max_explicit_stiffness(I, 1e-3, 2.0) - I * (2.0 / 2e-3) ^ 2) < 1e-9)
print("critical_damping(4, I) = ", critical_damping(4.0, I))
print("max_explicit_stiffness(I, 1 ms) = ", max_explicit_stiffness(I, 1e-3))
print("OK - joint range and drive round-trip through joint_state, is_unrestricted, the stop sizing")

-- the screw records and the contact kinds (bound since 2026-10-29), and the motor joints
-- of joint3dp, which Lua listed only in part until then
print("\nscrew axes, screw systems, contact kinds, joint families:")
local L = bivec3dp.new(0.0, 0.0, 1.0, 1.0, -2.0, 0.0)   -- through (2, 1, 0), along z
local s3 = pga.screw_axis(pga.get_motor(L, 0.7, 1.5))  -- the screw motor
assert(math.abs(s3.angle - 0.7) < 1e-12 and math.abs(s3.dist - 1.5) < 1e-12)
assert(math.abs(s3:pitch() - 1.5 / 0.7) < 1e-12)
local s2 = pga.screw_axis(pga.get_motor(vec2dp.new(1.0, 2.0, 1.0), 0.4))
assert(math.abs(s2.angle - 0.4) < 1e-12)
local sys = pga.screw_system({bivec3dp.new(0, 0, 1, 0, 0, 0), bivec3dp.new(0, 0, 1, 1, 0, 0)})
assert(sys.span == 2)
print("screw axis:", s3, "  screw system:", sys)
assert(contact_kind3dp.flat ~= contact_kind3dp.point and contact_kind2dp.point ~= nil)
for _, k in ipairs({"helical", "cylindrical", "spherical", "planar"}) do
    assert(joint3dp[k] ~= nil, "joint3dp." .. k .. " is not bound")
end

-- the viability bound: (sqrt(2 a_max d) - rate) / dt
assert(math.abs(stopping_accel_bound(0.5, 1.0, 4.0, 0.01) - (2.0 - 1.0) / 0.01) < 1e-12)
assert(stopping_accel_bound(-0.1, 0.5, 4.0, 0.01) == -50.0) -- past the wall: brake only
print("OK - screw records, contact kinds, the joint families, the viability bound")
