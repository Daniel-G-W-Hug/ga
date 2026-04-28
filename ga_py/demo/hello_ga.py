"""Hello-world for ga_py — verifies the install and demonstrates the two
core GA operations: the wedge product and the geometric product.

Run from the project root (see ga_py/README.md §3.4):

    PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python" \\
        ga_py/.venv/bin/python ga_py/demo/hello_ga.py

If this prints results without `ImportError`, the install is good.
"""
import ga_py


# EGA: build a vector and a bivector, combine them
v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
B = ga_py.ega.bivec3d(0.0, 0.0, 1.0)
m = v + B
print("m =", m)                # cross-grade multivector
print("|v| =", ga_py.ega.nrm(v))


# Wedge and geometric products --- the two GA fundamentals
v1 = ga_py.ega.vec3d(1.0, 2.0, 1.0)
v2 = ga_py.ega.vec3d(0.5, 3.0, 2.0)
print("v1 ^ v2 =", v1 ^ v2)    # wedge -> bivec3d (oriented area spanned by v1, v2)
print("v1 * v2 =", v1 * v2)    # geometric product -> mvec3d_e (scalar + bivector parts)


# Fundamental GA identity for vectors: v1 * v2 = dot(v1, v2) + wdg(v1, v2)
# Verify by extracting the grade-0 and grade-2 parts of the geometric product.
prod = v1 * v2
print("gr0(v1*v2) =", float(prod.gr0()), "  dot(v1, v2) =", float(ga_py.ega.dot(v1, v2)))
print("gr2(v1*v2) =", prod.gr2(), "   wdg(v1, v2) =", ga_py.ega.wdg(v1, v2))


# PGA: a point and a line, compute incidence
p = ga_py.pga.point2d(1.0, 2.0)
q = ga_py.pga.point2d(4.0, 6.0)
line = ga_py.pga.wdg(p, q)     # line joining two points
print("line =", line)
