-- script to demonstrate behavior of contractions
local u = vec2d.new(1,0)   -- e1
local v = vec2d.new(0,1)   -- e2
local B = u^v -- e12
--
print("contractions in 2D:")
print("u =", u)
print("v =", v)
print("B = ", B)
print()
print("u*B = ", u*B)
print("B>>u = ",B>>u)
print("B>>v = ",B>>v)
print("u<<B = ",u<<B)
print("v<<B = ",v<<B)
print("u^B = ", u^B)
print()
print("u*B = ", u*B)
print("B*u = ", B*u)
print("0.5*(u*B - B*u) = ", 0.5*(u*B - B*u))
print("0.5*(u*B + B*u) = ", 0.5*(u*B + B*u))
print()
--
local u = vec3d.new(1.5,-0.25,1.5)
local v = vec3d.new(1.5,0,-1.5)
local w = vec3d.new(0.75,2,0.1)
local B = u^v
--
print("contractions in 3D:")
print("u =", u)
print("v =", v)
print("w =", w)
print("B = ", B)
print()
print("project_onto(w,B) = ",project_onto(w,B))
print("B>>w = ", B>>w)
print("w<<B = ", w<<B)
print()