"""PGA3DP demo: two intersecting discs (wafer + grinding tool).

Python port of the second test case in
``ga_test/src/ga_appl3dp_appl_test.hpp``::

    TEST_CASE("pga3dp: intersecting discs (for grinding application)")

Side-by-side reading is intended — every line below has a direct
counterpart in the C++ source. Setup: see ``ga_py/README.md`` §3.

Run from the project root:

    PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python" \\
        ga_py/.venv/bin/python ga_py/demo/intersecting_discs.py

Geometric scenario
------------------
Two discs (wafer on a chuck and a ring-shaped grinding tool) sit
slightly overlapping in z. We integrate over the tool's annulus
(``r_min`` … ``r_max``, full revolution), test each ring-segment
centroid for being above the wafer, and accumulate contact area,
material-overlap volume, and relative speed at the contact point.

Notes for Python users
----------------------
* ``ga_py.pga.bulk_nrm(...)`` returns a typed ``scalar3dp`` —
  use ``float(...)`` before mixing with Python numbers
  (see README §4.3).
* ``pi`` and ``to_val`` are not exposed at ``ga_py.*``; use ``math.pi``
  and ``float(...)``.
* ``float * bivec3dp`` and ``bivec3dp >> vec3dp`` are bound, mirroring
  the C++ scalar product and right contraction. The right contraction
  ``>>`` has lower precedence than ``*`` in Python (same as C++), so
  ``rs * B >> v`` parses as ``(rs * B) >> v`` — what we want.
* GA types support f-string format specs (``f"{v:>-8.5f}"``); the spec
  is applied per component on the C++ side via ``fmt::format``.
"""
import math

import ga_py
from ga_py import pga


# ---------------------------------------------------------------------
# A disc as a Python class — direct mirror of the C++ `struct disc`.
# Validation in __init__ matches the C++ ctor's invariant checks.
# ---------------------------------------------------------------------
class Disc:
    """A flat disc (no thickness modelled directly).

    Parameters
    ----------
    cp     : centre point as ``vec3dp`` with w == 1.
    d      : diameter [mm] (>= 0).
    r_min  : inner radius [mm] (for ring-shaped tools); 0 for a full disc.
    r_max  : outer radius [mm].
    rs     : rotation speed [rad/s].
    u, v   : orthonormal in-plane direction vectors (w == 0, |u|=|v|=1, u·v=0).
    """

    def __init__(self, cp, d, r_min, r_max, rs, u, v):
        if cp.w != 1.0:
            raise ValueError("disc: cp must be a point (w=1).")
        if d <= 0.0:
            raise ValueError("disc: invalid diameter.")
        r = d / 2.0
        if r_min < 0.0 or r_max < 0.0 or r_min > r_max or r > r_max:
            raise ValueError("disc: invalid radius.")
        if float(pga.dot(u, v)) != 0.0:
            raise ValueError("disc: u and v must be perpendicular.")
        if float(pga.bulk_nrm(u)) != 1.0 or u.w != 0.0:
            raise ValueError("disc: u must be a unit-vector.")
        if float(pga.bulk_nrm(v)) != 1.0 or v.w != 0.0:
            raise ValueError("disc: v must be a unit-vector.")

        self.cp = cp
        self.r = r
        self.r_min = r_min
        self.r_max = r_max
        self.rs = rs
        self.u = u
        self.v = v

        # Derived quantities — same as the C++ initialiser list.
        self.B_uv = pga.wdg(u, v)                    # bivec3dp (disc plane)
        self.pl = pga.wdg(cp, self.B_uv)             # trivec3dp (plane equation)
        self.pl_normal = pga.l_weight_dual(self.pl)  # vec3dp


# ---------------------------------------------------------------------
# Wafer: 200 mm, 300 rpm, sitting on the e1-e2 plane through the origin.
# u, v are the global axes — wafer is not tilted.
# ---------------------------------------------------------------------
u = pga.vec3dp(1, 0, 0, 0)
v = pga.vec3dp(0, 1, 0, 0)

wafer = Disc(pga.O_3dp, 200.0, 0.0, 100.0, ga_py.rpm2radps(300.0), u, v)

print("pga3dp: intersecting discs (for grinding application)\n")
print(f"wafer.cp         = {wafer.cp}")
print(f"wafer.r [mm]     = {wafer.r}")
print(f"wafer.r_min [mm] = {wafer.r_min}")
print(f"wafer.r_max [mm] = {wafer.r_max}")
print(f"wafer.rs [rad/s] = {wafer.rs}")
print(f"wafer.pl         = {wafer.pl:>-8.5f}")
print(f"wafer.pl_normal  = {wafer.pl_normal:>-8.5f}")
print(f"u                = {u:>-8.5f}")
print(f"v                = {v:>-8.5f}")
print()


# ---------------------------------------------------------------------
# Tool: 200 mm OD, 95 mm ID (ring), 2000 rpm, centred 150 mm out in +e1
# and 0.5 mm below the wafer (overlapping in -e3 by 0.5 mm).
# Note: the tool uses (u, -v) so its B_uv has the opposite orientation
# of the wafer — same physical plane, opposite spin direction.
# ---------------------------------------------------------------------
tc = pga.vec3dp(150, 0, -0.5, 1)
tool = Disc(tc, 200.0, 95.0, 100.0, ga_py.rpm2radps(2000), u, -v)

print(f"tool.cp         = {tool.cp}")
print(f"tool.r [mm]     = {tool.r}")
print(f"tool.r_min [mm] = {tool.r_min}")
print(f"tool.r_max [mm] = {tool.r_max}")
print(f"tool.rs [rad/s] = {tool.rs}")
print(f"tool.pl         = {tool.pl:>-8.5f}")
print(f"tool.pl_normal  = {tool.pl_normal:>-8.5f}")
print(f"u               = {u:>-8.5f}")
print(f"v               = {-v:>-8.5f}")
print()


# ---------------------------------------------------------------------
# Detour from the C++ source: shifting a plane by x while keeping its
# orientation unchanged. P ^ u ^ v + x ^ u ^ v == (P + x) ^ (u ^ v).
# ---------------------------------------------------------------------
x = pga.vec3dp(0, 0, 0.05, 0)
pl_shifted = wafer.pl + pga.wdg(x, wafer.B_uv)
assert pl_shifted == pga.wdg(wafer.cp + x, wafer.B_uv)


# ---------------------------------------------------------------------
# Integrate over the tool annulus: for every ring-segment centroid,
# rotate it around the tool axis, project onto the wafer plane, and
# accumulate area / volume / contact speed when there's a real overlap.
# ---------------------------------------------------------------------
volume = 0.0
area = 0.0
min_tot_speed = float("inf")
max_tot_speed = -float("inf")

if float(pga.bulk_nrm(tool.cp - wafer.cp)) < tool.r + wafer.r:

    nr = 25                                    # radial steps
    dr = (tool.r_max - tool.r_min) / nr

    nphi = 720                                 # angular steps
    dphi = 2.0 * math.pi / nphi

    print(f"nr = {nr}, dr = {dr:>-8.5f}, nphi = {nphi}, dphi = {dphi:>-8.5f}\n")

    for j in range(nr):
        r_inner = tool.r_min + j * dr          # current inner radius
        rm = r_inner + 0.5 * dr                # mean radius (for area)

        for i in range(nphi):
            phi = i * dphi
            phim = phi + 0.5 * dphi            # mean angle (centroid)

            dA = rm * dphi * dr

            # Geometric centroid of the annular ring segment.
            r_i_sq = r_inner * r_inner
            r_i_tr = r_i_sq * r_inner
            r_o = r_inner + dr
            r_o_sq = r_o * r_o
            r_o_tr = r_o_sq * r_o

            r_gc = (2.0 * (r_o_tr - r_i_tr) * math.sin(0.5 * dphi) /
                    (3.0 * (r_o_sq - r_i_sq) * 0.5 * dphi))
            ref_pos = tool.cp + r_gc * u

            # Rotation line (axis of the tool) and motor for angle phim.
            rot_line = pga.wdg(tool.cp, tool.pl_normal)
            cur_pos = pga.move3dp(ref_pos, pga.get_motor(rot_line, phim))

            # Orthogonal projection of cur_pos onto the wafer plane.
            cur_proj_pos = pga.rwdg(
                wafer.pl, pga.r_weight_expand3dp(cur_pos, wafer.pl)
            )

            # Inside the wafer footprint?
            if float(pga.bulk_nrm(cur_proj_pos - wafer.cp)) <= wafer.r:
                delta_vec = cur_pos - cur_proj_pos

                # Real overlap only when delta_vec points along the tool
                # normal (i.e. the tool material reaches the wafer).
                if ga_py.sign(float(pga.dot(delta_vec, tool.pl_normal))) == 1:
                    area += dA
                    volume += dA * float(pga.bulk_nrm(delta_vec))

                    # Surface speeds at the contact point. Note Python's
                    # `>>` precedence: `rs * B >> v` is `(rs * B) >> v`.
                    wafer_speed = (wafer.rs * wafer.B_uv) >> (cur_proj_pos - wafer.cp)
                    tool_speed = (tool.rs * tool.B_uv) >> (cur_pos - tool.cp)
                    tot_speed_vec = tool_speed + wafer_speed
                    tot_speed = float(pga.bulk_nrm(tot_speed_vec))

                    if tot_speed != 0.0:
                        min_tot_speed = min(min_tot_speed, tot_speed)
                        max_tot_speed = max(max_tot_speed, tot_speed)


print(f"area = {area:>-8.5f} mm^2")
print(f"volume = {volume:>-8.5f} mm^3")
print(f"min_tot_speed = {min_tot_speed:>-8.5f} mm/s")
print(f"max_tot_speed = {max_tot_speed:>-8.5f} mm/s")
print()
print(f"wafer speed at r = {wafer.r_max} mm = "
      f"{wafer.rs * wafer.r_max:>-8.5f} mm/s")
print(f"tool  speed at r = {tool.r_max} mm = "
      f"{tool.rs * tool.r_max:>-8.5f} mm/s")
print()


# The loop is intentionally pure-Python so the data flow is readable.
# With nphi=720, nr=25 this is ~18 000 iterations of a few PGA ops
# each — expect a runtime of a few seconds. For production use,
# vectorising over raw component arrays (numpy buffer protocol) is
# the right next step; see README §6.3.
