-- script to check result of expansion in ega3d
local u = vec3d.new(1,0,0)
local v = vec3d.new(0,1,0)
local w = vec3d.new(1,1,1)
local B = u^v
--
print("expansions in 3D:")
print("u =", u)
print("v =", v)
print("w =", w)
print("B = ", B)
print()
local e = w ^ dual(B)
print("e = expand(w,B) = w ^ dual(B) = ", w ^ dual(B))
print("dual(e) = ", dual(e)) -- normal vector of plane 
print()