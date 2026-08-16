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
