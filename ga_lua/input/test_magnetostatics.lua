-- test_magnetostatics.lua - the magnetostatics of ga_usr_magnetostatics.hpp
--
-- These verify the BINDING: the C++ suites pin the mathematics (ga_appl2dp_test,
-- "magnetostatics 2d", and ga_appl3dp_test, "magnetostatics 3d"). What can only go
-- wrong at this layer is the 2D / 3D overload selection, the table conversions, and
-- stress_wrench taking its field as a Lua function. References are closed forms written
-- out here.

print("=======================================================")
print("Testing magnetostatics (ga_usr_magnetostatics.hpp)")
print("=======================================================")

local function rel(a, b) return math.abs(a - b) / math.abs(b) end

print("\n1. mu0 and the field of a line current (2D)")
assert(rel(mu0, 4.0e-7 * math.pi) < 1e-15, "mu0")
local c = line_current.new(vec2d.new(0.3, -0.2), pscalar2d.new(250.0))
local b = b_line_current(vec2d.new(1.0, 0.0), c)
local dx, dy = 0.7, 0.2
local k = mu0 * 250.0 / (2.0 * math.pi * (dx * dx + dy * dy))
assert(rel(b.x, -k * dy) < 1e-14 and rel(b.y, k * dx) < 1e-14, "b_line_current")
print("  " .. tostring(c) .. " -> b = " .. tostring(b))

print("\n2. traction: the 2D and the 3D overload, Maxwell's tensor")
local t2 = traction(vec2d.new(0.9, -0.4), vec2d.new(0.6, 0.8))
local bn, b2 = 0.9 * 0.6 - 0.4 * 0.8, 0.81 + 0.16
assert(rel(t2.x, (bn * 0.9 - 0.5 * b2 * 0.6) / mu0) < 1e-14, "traction 2D")
local t3 = traction(vec3d.new(0.9, -0.4, 0.3), vec3d.new(0.0, 0.6, 0.8))
bn, b2 = -0.4 * 0.6 + 0.3 * 0.8, 0.81 + 0.16 + 0.09
assert(rel(t3.z, (bn * 0.3 - 0.5 * b2 * 0.8) / mu0) < 1e-14, "traction 3D")
print("  2D " .. tostring(t2) .. ", 3D " .. tostring(t3))

print("\n3. a long segment is the line current; tables of segments")
local s = current_segment.new(vec3d.new(0, 0, -1.0e4), vec3d.new(0, 0, 1.0e4), 100.0)
local bs = b_segment(vec3d.new(0.5, 0, 0), s)
assert(rel(bs.y, mu0 * 100.0 / (2.0 * math.pi * 0.5)) < 1e-8, "b_segment")
local loop = circle_loop(vec3d.new(0, 0, 0), vec3d.new(0, 0, 1), 0.01, 256, 40.0)
assert(#loop == 256, "circle_loop returns a table of 256 segments")
local far = vec3d.new(0, 0, 1.0)
local mu = vec3d.new(0, 0, 40.0 * 0.5 * 256 * 1e-4 * math.sin(2.0 * math.pi / 256))
assert(rel(b_segments(far, loop).z, b_dipole(far, mu).z) < 1e-3, "a loop is a dipole")
print("  b_segment " .. tostring(bs) .. "; loop on axis at 100 R ~ the dipole")

print("\n4. stress_wrench with a Lua field: two parallel line currents")
local w1 = line_current.new(vec2d.new(0.10, 0.05), pscalar2d.new(300.0))
local w2 = line_current.new(vec2d.new(0.45, -0.25), pscalar2d.new(-180.0))
local both = { w1, w2 }
local W = stress_wrench(function(r) return b_line_currents(r, both) end,
                        circle_contour(w1.p, 0.1, 720))
local F1 = lorentz_force(b_line_current(w1.p, w2), w1.I)
-- the force on its line of action through p1: W = wdg(p1, F1)
local Wref = pga.wdg(vec2dp.new(0.10, 0.05, 1.0), vec2dp.new(F1.x, F1.y, 0.0))
assert(rel(W.x, Wref.x) < 1e-9 and rel(W.y, Wref.y) < 1e-9 and rel(W.z, Wref.z) < 1e-9,
       "the 2D stress wrench")
print("  W = " .. tostring(W) .. ", wdg(p1, F1) = " .. tostring(Wref))

print("\n5. stress_wrench in 3D: a dipole in another's field")
local pj, mi, mj = vec3d.new(0.06, -0.02, 0.03), vec3d.new(0, 0, 2.0), vec3d.new(0.5, 1.2, -0.4)
local W3 = stress_wrench(function(r)
    return b_dipole(r, mi) + b_dipole(r - pj, mj)
end, sphere_surface(pj, 0.02, 8, 96))
-- the C++ gate's reference: the force line through pj plus the couple
-- T = mu_j x b_i(pj) (the force taken from W3 itself: this checks the binding and
-- the couple, the C++ suite checks the force against Yung et al.)
local T = ega.dual(ega.wdg(mj, b_dipole(pj, mi)))
local Wref3 = pga.wdg(vec3dp.new(pj.x, pj.y, pj.z, 1.0), vec3dp.new(W3.vx, W3.vy, W3.vz, 0.0))
    + bivec3dp.new(0.0, 0.0, 0.0, T.x, T.y, T.z)
-- against the moment's norm: its z component vanishes here (1e-17), so a per-component
-- relative error would measure rounding
local function mdiff(A, B)
    local d = math.max(math.abs(A.mx - B.mx), math.abs(A.my - B.my), math.abs(A.mz - B.mz))
    return d / math.sqrt(B.mx * B.mx + B.my * B.my + B.mz * B.mz)
end
assert(mdiff(W3, Wref3) < 1e-8, "the 3D wrench carries the couple")
-- the falsifier: without the couple the moment is wrong
local Wnoc = pga.wdg(vec3dp.new(pj.x, pj.y, pj.z, 1.0), vec3dp.new(W3.vx, W3.vy, W3.vz, 0.0))
assert(mdiff(W3, Wnoc) > 0.1, "the couple is not negligible")
print("  W = " .. tostring(W3))

print("\n6. a magnetised box: its field, its orientation, the wrench on it")
local one = mvec3d_e.new(1, 0, 0, 0)
local cube = magnetised_box.new(vec3d.new(0, 0, 0), vec3d.new(0.005, 0.005, 0.005), one,
                                vec3d.new(0, 0, 1.0e6))
-- the default box: the identity rotor, as in C++
assert(magnetised_box().R.c0 == 1.0, "the default orientation is the identity")
-- inside a cube the demagnetising field is -m/3, so b = (2/3) mu0 m at its centre
assert(rel(b_box(vec3d.new(0, 0, 0), cube).z, 2.0 / 3.0 * mu0 * 1.0e6) < 1e-14,
       "b_box inside")
-- outside it approaches its dipole mu = 8 h_x h_y h_z m at (a/r)^4: 16x per doubling
local mu_c = vec3d.new(0, 0, 8 * 0.005 ^ 3 * 1.0e6)
local function dev(z)
    local r = vec3d.new(0, 0, z)
    return rel(b_box(r, cube).z, b_dipole(r, mu_c).z)
end
local q = dev(0.1) / dev(0.2)
assert(q > 15.5 and q < 16.5, "the cube differs from its dipole at (a/r)^4")
-- R reaches the code: a flat box turned 90 degrees about z is the box with its x and y
-- extents swapped (the magnetisation is in world coordinates), and not the unturned one
local Rz = ega.get_rotor(bivec3d.new(0, 0, 1), math.pi / 2)
local mt = vec3d.new(0.3e6, 0.8e6, 0.5e6)
local turned = magnetised_box.new(vec3d.new(0.01, 0, 0), vec3d.new(0.01, 0.004, 0.002),
                                  Rz, mt)
local swapped = magnetised_box.new(vec3d.new(0.01, 0, 0), vec3d.new(0.004, 0.01, 0.002),
                                   one, mt)
local plain = magnetised_box.new(vec3d.new(0.01, 0, 0), vec3d.new(0.01, 0.004, 0.002),
                                 one, mt)
local rp = vec3d.new(0.02, 0.013, 0.007)
local function vdiff(A, B)
    return math.max(math.abs(A.x - B.x), math.abs(A.y - B.y), math.abs(A.z - B.z))
        / math.sqrt(B.x * B.x + B.y * B.y + B.z * B.z)
end
assert(vdiff(b_box(rp, turned), b_box(rp, swapped)) < 1e-12, "the orientation is used")
assert(vdiff(b_box(rp, turned), b_box(rp, plain)) > 0.1, "and it matters here")
-- in a uniform field the wrench is a pure couple, T = mu x B = dual(mu ^ B)
local B0 = vec3d.new(0.2, -0.1, 0.4)
local Wu = box_wrench(function(p) return B0 end, turned)
local Tu = ega.dual(ega.wdg(8 * 0.01 * 0.004 * 0.002 * mt, B0))
assert(math.max(math.abs(Wu.vx), math.abs(Wu.vy), math.abs(Wu.vz)) < 1e-12,
       "no force in a uniform field")
assert(vdiff(vec3d.new(Wu.mx, Wu.my, Wu.mz), Tu) < 1e-12, "the couple mu x B")
-- panels reaches the code: in the field of a cube 12.6 mm below (the gap between them
-- 2.6 mm) one panel per face side and four differ, by little
local src = magnetised_box.new(vec3d.new(0, 0, -0.0126), vec3d.new(0.005, 0.005, 0.005),
                               one, vec3d.new(0, 0, 1.0e6))
local near_b = function(p) return b_box(p, src) end
local W1, W4 = box_wrench(near_b, cube), box_wrench(near_b, cube, 4)
assert(W1.vz ~= W4.vz and rel(W1.vz, W4.vz) < 1e-8, "the panels argument")
assert(W4.vz < -19.0, "the cubes attract")
print("  " .. tostring(cube))
print("  force between the two cubes: " .. W4.vz .. " N")

print("\nAll magnetostatics binding checks passed.")
