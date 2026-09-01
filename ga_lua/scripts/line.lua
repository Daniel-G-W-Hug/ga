--demo line as a bivector and its attitude
local P = vec2dp.new(1, 1, 1)
local Q = vec2dp.new(3, 2, 1)
local ln_PQ = P ^ Q

print("P = ", P)
print("Q = ", Q)
print("ln_PQ = ", ln_PQ)
print()
print("att(ln_PQ)            = ", pga.att(ln_PQ))
print("rwdg(ln_PQ, cmpl(e3)) = ", pga.rwdg(ln_PQ, pga.cmpl(e3_2dp)))
-- attitude returns the direction vector of the line
