-- test_geodesics.lua - the geodesy layer from Lua: angles in degree/minute/second text,
-- the positions built from them, and the Helmert datum transformation with its exact
-- inverse (all bound in Lua; the dms text and the inverse since 2026-10-29)

print("=======================================================")
print("Testing the geodesy layer")
print("=======================================================")

local function close(a, b, tol) return math.abs(a - b) < tol end
local deg = math.pi / 180.0

-- degree/minute/second text, both ways
assert(close(dms2deg("52°31'12\"N", geo_angle.latitude), 52.52, 1e-12))
assert(close(dms2deg("13°24'36\"W", geo_angle.longitude), -13.41, 1e-12))
assert(close(dms2rad("52°31'12\"N", geo_angle.latitude), 52.52 * deg, 1e-15))
assert(deg2dms(52.52, geo_angle.latitude) == "52°31'12.0\"N")
assert(rad2dms(-13.41 * deg, geo_angle.longitude, 0) == "13°24'36\"W")
print("dms2deg / deg2dms / dms2rad / rad2dms: OK")

-- positions quoted in dms text
local p = to_geo_pos(geo_pos_dms("52°31'12\"N", "13°24'36\"E", 34.0))
assert(close(p.lat, 52.52 * deg, 1e-15) and close(p.lon, 13.41 * deg, 1e-15))
local q = to_geo_pos(geo_pos_dms2dp("52°31'12\"N", 34.0))
assert(close(q.lat, 52.52 * deg, 1e-15) and q.height == 34.0)
print("to_geo_pos(geo_pos_dms), to_geo_pos(geo_pos_dms2dp):", p, q)

-- the Helmert transformation (the German DHDN set) and its exact inverse: each inverse
-- undoes its namesake to the rounding of the coordinates; the other convention's
-- inverse is about 110 m off
local h = pga.make_helmert7(598.1, 73.7, 418.2, 0.202, 0.045, -2.455, 6.7)
local X = vec3dp.new(4157222.543, 664789.307, 4774952.099, 1.0)
local function dist(A, B)
    return math.sqrt((A.x - B.x) ^ 2 + (A.y - B.y) ^ 2 + (A.z - B.z) ^ 2)
end
local pv = pga.helmert_position_vector(X, h)
local cf = pga.helmert_coordinate_frame(X, h)
local e_pv = dist(pga.helmert_position_vector_inverse(pv, h), X)
local e_cf = dist(pga.helmert_coordinate_frame_inverse(cf, h), X)
local e_mix = dist(pga.helmert_coordinate_frame_inverse(pv, h), X)
print(string.format("Helmert round trips %.1e m / %.1e m; the other convention %.1f m",
                    e_pv, e_cf, e_mix))
assert(e_pv < 1e-8 and e_cf < 1e-8 and e_mix > 100.0)

print("\n=======================================================")
print("Geodesy tests completed successfully")
print("=======================================================")
