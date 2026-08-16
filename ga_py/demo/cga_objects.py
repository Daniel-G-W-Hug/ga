"""CGA in Python — conformal objects, their properties, and the versors.

The point of conformal GA: a circle, a sphere, a line and a plane are all just
multivectors, so ONE operation intersects any pair of them and the result is
another object of the same kind. There is no case analysis — the squared radius
of the result carries it (positive: two real intersections, zero: tangency,
negative: none).

Run from the project root (see ga_py/README.md §3.3 for PYTHONPATH):

    PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python" \\
        ga_py/.venv/bin/python ga_py/demo/cga_objects.py

The numbers printed here are the ones pinned in ga_py/tests/test_cga.py and in
the C++ suite (ga_cga_test), so a divergence is a failure, not a variant.
"""
import math

from ga_py import cga


print("--- 1. points are null vectors -------------------------------------")
a = cga.round_point3dc(0.0, 0.0, 0.0, 0.0)  # the origin, radius 0
b = cga.round_point3dc(3.0, 4.0, 0.0, 0.0)  # the point (3,4,0)
print("a           =", a)
print("b           =", b)
print("a . a       =", float((a * a).gr0()), "(null: a point squares to zero)")
print("a . b       =", float(cga.dot(a, b)), "= -1/2 * distance^2 = -12.5")


print("\n--- 2. rounds are joins of points ----------------------------------")
p1 = cga.round_point3dc(1.0, 0.0, 0.0, 0.0)
p2 = cga.round_point3dc(0.0, 1.0, 0.0, 0.0)
p3 = cga.round_point3dc(0.0, 0.0, 1.0, 0.0)
p4 = cga.round_point3dc(-1.0, 0.0, 0.0, 0.0)
print("p1 ^ p2         =", cga.wdg(p1, p2), "-> dipole (point pair)")
print("p1 ^ p2 ^ p3    =", cga.wdg(cga.wdg(p1, p2), p3), "-> circle")
sphere = cga.wdg(cga.wdg(cga.wdg(p1, p2), p3), p4)
print("p1^p2^p3^p4     =", sphere, "-> sphere")
print("== unit sphere  :", cga.is_congruent(sphere, cga.sphere3dc(0.0, 0.0, 0.0, 1.0),
                                            1e-12))
# a flat is the same join with the point at infinity appended
print("p1 ^ e5         =", cga.wdg(p1, cga.e5_3dc), "-> flat point")


print("\n--- 3. an object knows its own geometry ----------------------------")
s = cga.unitize(cga.sphere3dc(3.0, 4.0, 0.0, 2.0))
print("sphere at (3,4,0), r = 2")
print("  radius_sq   =", cga.radius_sq(s))
print("  center_nrm  =", cga.center_nrm(s), "= |(3,4,0)|")
print("  cen         =", cga.cen(s), "(the center, as a round point of radius 2)")
circle = cga.circle3dc(1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 1.0)
print("circle at (1,2,3), r = 2, normal +z")
print("  car         =", cga.car(circle), "(carrier: its plane)")
print("  con         =", cga.con(circle), "(container: its sphere)")
print("  att         =", cga.att(circle), "(attitude)")
print("  car ^ con   ==  the circle itself:",
      cga.is_congruent(cga.rwdg(cga.car(circle), cga.con(circle)), circle, 1e-12))


print("\n--- 4. one intersection, no case analysis --------------------------")
ball = cga.sphere3dc(0.0, 0.0, 0.0, 2.0)


def chord(sphere_, line_):
    """Intersect and, if the meet is real, return the two points."""
    d = cga.rwdg(sphere_, line_)
    r_sq = cga.radius_sq(d)
    if r_sq < 0.0:
        return r_sq, None
    du = cga.unitize(d)
    c, n, r = cga.unitize(cga.cen(du)), cga.att(du), math.sqrt(cga.radius_sq(du))
    return r_sq, ((c.x + r * n.x, c.y + r * n.y, c.z + r * n.z),
                  (c.x - r * n.x, c.y - r * n.y, c.z - r * n.z))


for y, label in ((1.0, "off-center secant"), (2.0, "tangent"), (3.0, "no intersection")):
    r_sq, pts = chord(ball, cga.line3dc(0.0, y, 0.0, 1.0, 0.0, 0.0))
    print(f"  line at y={y}  r^2 = {r_sq:+.1f}   {label}")
    if pts is not None:
        print("                hits:", tuple(round(c, 7) for c in pts[0]),
              tuple(round(c, 7) for c in pts[1]))

cut = cga.rwdg(ball, cga.plane3dc(0.0, 0.0, 1.0, 1.0))
print("  sphere ^ plane at z=1  ->  circle of r^2 =", cga.radius_sq(cut),
      "= 4 - 1  (centered on the plane)")

# the same question in the plane: two circles
c1, c2 = cga.circle2dc(0.0, 0.0, 2.0), cga.circle2dc(3.0, 0.0, 2.0)
print("  circle ^ circle (2dc)  ->  point pair of r^2 =",
      cga.radius_sq(cga.rwdg(c1, c2)))


print("\n--- 5. conformal versors -------------------------------------------")
T = cga.get_translation(1.0, 2.0, 3.0)
moved = cga.unitize(cga.transform(ball, T))
print("translate:  center ->", cga.cen(moved), " r^2 =", cga.radius_sq(moved),
      "(unchanged)")
D = cga.get_dilation(0.0, 0.0, 0.0, 3.0)
print("dilate x3:  r^2 ->", cga.radius_sq(cga.unitize(cga.transform(ball, D))),
      "= 4 * 3^2")
R = cga.get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.5)
print("rotate:     rexp(rlog(R)) == R:", cga.is_same_transform(cga.rexp(cga.rlog(R)),
                                                               R, 1e-12))
print("            M and -M act alike:", cga.is_same_transform(R, -R, 1e-12),
      "(versors double cover the transformations)")

# inversion in a sphere is a flector, and its own inverse
unit = cga.sphere3dc(0.0, 0.0, 0.0, 1.0)
q = cga.round_point3dc(2.0, 0.0, 0.0, 0.0)
once = cga.unitize(cga.invert_on(q, unit))
twice = cga.unitize(cga.invert_on(cga.invert_on(q, unit), unit))
print("invert in the unit sphere:  x = 2 ->", round(once.x, 7),
      "-> back to", round(twice.x, 7))
