"""PGA3DP demo: reference plane and tumbling tool plane (grinding application).

Python port of the first test case in
``ga_test/src/ga_appl3dp_appl_test.hpp``::

    TEST_CASE("pga3dp: reference and tumbling plane (for grinding application)")

Side-by-side reading is intended — every line below has a direct
counterpart in the C++ source. The math is identical; only the host
language differs. See ``ga_py/README.md`` §3 for setup.

Run from the project root:

    PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python" \\
        ga_py/.venv/bin/python ga_py/demo/tumbling_plane.py

Geometric scenario
------------------
A reference plane ``p_ref`` is reconstructed from three measured points
A, B, C on a real (slightly tilted) workpiece. A grinding-tool axis
``z`` should be perpendicular to that plane, but in practice the tip
tumbles in a small circle of radius ``r``. For each tumbling angle
``phi`` we build the current tool plane ``p_tool`` and intersect it
with ``p_ref`` to get the contact line ``l``. The line direction
``att(l)`` and its angle to ``e1`` characterise the contact geometry.

Notes for Python users
----------------------
* ``ga_py.pga.bulk_nrm(...)`` returns a typed ``scalar3dp`` wrapper —
  wrap in ``float(...)`` before doing Python arithmetic with it
  (see README §4.3 "Typed scalars vs. plain floats").
* ``pi`` is not re-exported as a top-level constant; use ``math.pi``.
* GA objects support f-string format specs (``f"{v:>-8.5f}"``) just like
  built-in numeric types. The spec is dispatched through ``fmt::format``
  on the C++ side and applied **per component** — the typical PEP 3101
  fields (``f``/``e``/``g``, width, precision, sign, alignment) work
  the same as for Python ``float``. Bad specs raise ``ValueError``.
"""
import math

import ga_py
from ga_py import pga  # `from ga_py import ega/pga` works (attr access)


# ---------------------------------------------------------------------
# Reference plane from three measured points (close to the e1-e2 plane,
# slightly tilted by a 0.01 z-offset on B to mimic a real measurement).
# C++: vec3dp(1.0, 0.0, 0.0, 1)  → Python: pga.vec3dp(1.0, 0.0, 0.0, 1)
# The fourth coordinate is the projective weight: w=1 → point, w=0 → direction.
# ---------------------------------------------------------------------
A = pga.vec3dp(1.0, 0.0, 0.0, 1)
B = pga.vec3dp(1.0, 1.0, 0.01, 1)
C = pga.vec3dp(0.0, 1.0, 0.0, 1)

# Plane through A, B, C as a wedge of three points, then unitized so the
# weight-norm is 1 (canonical representative of the plane subspace).
p_ref = pga.unitize(pga.wdg(pga.wdg(A, B), C))     # trivec3dp
# Plane normal as a (direction) vector — left weight dual of the plane.
p_n = pga.l_weight_dual(p_ref)                     # vec3dp

# f-string format specs work directly on GA types. The spec is applied
# per component (matches the C++ `fmt::println("{:>-8.5f}", ...)` output).
print("pga3dp: reference and tumbling plane (for grinding application)\n")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print(f"           p_ref  = {p_ref:>-8.5f}")
print(f"           p_n    =    {p_n:>-8.5f}, |p_n| = {float(pga.bulk_nrm(p_n))}")
print()


# ---------------------------------------------------------------------
# Tool tip T (centre of the rotating grinding tool) close to p_ref,
# nominal tool axis z in +e3 direction.
# ---------------------------------------------------------------------
T = pga.vec3dp(3.0, 1.0, 0.0, 1)
z = pga.vec3dp(0.0, 0.0, 1.0, 0)


# ---------------------------------------------------------------------
# Sweep the tumbling angle phi from 0 to 2*pi in 15-degree steps. At
# each step compute the current tool plane and its intersection with
# p_ref (the contact line).
# ---------------------------------------------------------------------
r_dim = 0.01                                       # tumbling radius (1% of |z|)

phi = 0.0
step = ga_py.deg2rad(15)                           # `deg2rad` is at ga_py.* top level
while phi < 2.0 * math.pi:
    # Tumbling offset r at the tip — a direction vector (w=0) of length r_dim.
    r = r_dim * pga.vec3dp(math.cos(phi), math.sin(phi), 0.0, 0)

    # Point R is the tip of the (slightly tumbled) tool axis. The
    # normalized vector T → R is the current tool-plane normal.
    R = T + z + r
    z_n = (R - T) / float(pga.bulk_nrm(R - T))

    print(f"phi = {math.degrees(phi):>3.0f}: z_n         =    {z_n:>-8.5f}")

    # Build the current tool plane in two steps:
    #   1) line through T in direction z_n
    #   2) plane containing T, perpendicular to that line
    l_h = pga.unitize(pga.wdg(z_n, T))
    p_tool = pga.unitize(pga.wdg(T, pga.r_weight_dual(l_h)))
    print(f"           p_tool      = {p_tool:>-8.5f}")

    # Intersection of the two planes — a line (bivec3dp). `att(l)` is
    # the line's direction (an ideal point at infinity).
    l = pga.rwdg(p_tool, p_ref)
    print(f"           att(l)      =    {pga.att(l):>-8.5f}")
    print(f"           angle(l,e1) =    "
          f"{math.degrees(pga.angle(pga.att(l), pga.e1_3dp)):>5.2f}°")
    print()

    phi += step
