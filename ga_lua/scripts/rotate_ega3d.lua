-- rotation demo in EGA3D
local a = vec3d.new(1,0,0)
local b = vec3d.new(math.cos(deg2rad(15)),math.sin(deg2rad(15)),0)
local c = vec3d.new(math.cos(deg2rad(30)),math.sin(deg2rad(30)),0)
local R = b*a -- rotate from a to b
--
print("a = ", a)
print("b = ", b)
print("R = ", R)
print()
print("c = ", c)
print()
print("ega.gr1(R*a*ega.rev(R)) = ", ega.gr1(R*a*ega.rev(R)))
print()