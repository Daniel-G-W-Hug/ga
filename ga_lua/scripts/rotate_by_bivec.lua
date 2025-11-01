-- demo of rotation by bivector multiplication
local u = vec2d.new(3,1)
local e12 = pscalar2d.new(1) -- a bivector
--
print("u   = ", u)
print("e12 = ", e12)
print()
print("e12*u = ",e12*u)
print("u*e12 = ", u*e12)
