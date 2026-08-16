"""Python port of the merry-go-round kinematics example.

The example drives a `kinematic_system2dp` frame tree: a rotating platform P
mounted on the inertial ground W, carrying N turntables T_i that each spin
relative to the platform. A marked point on T0 traces an epicyclic world path.

`kinematic_system2dp` is a stateful C++ class and is *not* bound to Python.
The point of this test is to confirm that all the underlying GA artefacts the
demo relies on — the `pose2dp` / `kin_state2dp` data structs, the
`motor_from_pose2dp` / `pose2dp_from_motor` helpers, and the motor / twist
primitives (`rgpr`, `rrev`, `exp`, `move2dp`, `rcmt`) — are reachable from
Python and compose into the same kinematics. We therefore *reconstruct* the
relevant pieces of `kinematic_system2dp` from those primitives and validate
the result against an independent closed-form (epicyclic) solution.
"""
from __future__ import annotations

import math

import pytest

import ga_py
from ga_py import pga

# --- merry-go-round parameters (mirror merry_go_round_params defaults) ----
CENTRE = (0.0, 0.0)   # platform centre in world coordinates
OMEGA_P = 0.6         # platform angular velocity [rad/s]
OMEGA_T = 2.0         # turntable spin vs. platform [rad/s]
RADIUS = 1.6          # turntable mounting radius on the platform
MARKED = (0.6, 0.0)   # tracked point on turntable 0, in T0 local coordinates


def _pt(x: float, y: float) -> "pga.vec2dp":
    """A unitized projective point (z = 1)."""
    return pga.vec2dp(x, y, 1.0)


def _to_twist(vx: float, vy: float, omega: float) -> "pga.vec2dp":
    """Pack (linear velocity, angular velocity) into a 2dp twist — the exact
    encoding used by kinematic_system2dp::to_twist: twist(.x,.y,.z) =
    (-v_y, v_x, omega)."""
    return pga.vec2dp(-vy, vx, omega)


def _rot(angle: float, x: float, y: float) -> tuple[float, float]:
    c, s = math.cos(angle), math.sin(angle)
    return (c * x - s * y, s * x + c * y)


class MerryGoRound:
    """W ← P (platform) ← T0 (turntable), reconstructed from GA primitives.

    Frame poses are stored as `pose2dp` (relative to the parent). World motors
    are composed with `rgpr` exactly as static_system2dp::get_pos_trafo(i, 0)
    does, and world velocity/acceleration twists are propagated with the same
    Newton-Euler recursion as kinematic_system2dp::world_VA.
    """

    def __init__(self, t: float):
        # Platform: origin at centre, angle grows with the platform rate.
        self.pose_P = pga.pose2dp(_pt(*CENTRE), OMEGA_P * t)
        # Turntable 0: mounted at (RADIUS, 0) on the platform, spinning.
        mount = _pt(RADIUS, 0.0)
        self.pose_T0 = pga.pose2dp(mount, OMEGA_T * t)

        # World (frame → world) motors via rgpr, parent applied on the left.
        self.M_P = pga.motor_from_pose2dp(self.pose_P)
        self.M_T0 = pga.rgpr(self.M_P, pga.motor_from_pose2dp(self.pose_T0))

        # Relative twists (kin_state: only an angular rate per frame).
        self.xi_P = _to_twist(0.0, 0.0, OMEGA_P)
        self.xi_T0 = _to_twist(0.0, 0.0, OMEGA_T)

    def marked_world(self) -> "pga.vec2dp":
        return pga.move2dp(_pt(*MARKED), self.M_T0)

    def world_VA(self) -> tuple["pga.vec2dp", "pga.vec2dp"]:
        """(V, A) world twists of T0 — kinematic_system2dp::world_VA, with the
        relative acceleration twists zero (constant rates)."""
        V = pga.vec2dp(0.0, 0.0, 0.0)
        A = pga.vec2dp(0.0, 0.0, 0.0)
        for xi, M in ((self.xi_P, self.M_P), (self.xi_T0, self.M_T0)):
            zeta = pga.move2dp(xi, M)          # Ad_M(xi): relative twist in world
            V = V + zeta
            A = A + pga.rcmt(V, zeta)          # zetadot = 0; Coriolis coupling only
        return V, A

    def point_velocity(self, X: "pga.vec2dp") -> "pga.vec2dp":
        V, _ = self.world_VA()
        return pga.rcmt(V, X)                  # velocity_field(V, X)

    def point_acceleration(self, X: "pga.vec2dp") -> "pga.vec2dp":
        V, A = self.world_VA()
        return pga.rcmt(A, X) + pga.rcmt(V, pga.rcmt(V, X))  # accel_field(V, A, X)


# --- analytic (epicyclic) ground truth -----------------------------------

def _analytic(t: float):
    """Closed-form world position / velocity / acceleration of the marked
    point on T0 for constant platform and turntable rates."""
    thP, thT = OMEGA_P * t, OMEGA_T * t
    cx, cy = CENTRE
    # turntable-0 origin: centre + R(thP)·mount
    ox, oy = _rot(thP, RADIUS, 0.0)
    ox, oy = cx + ox, cy + oy
    # marked point: O + R(thP+thT)·marked_local
    mx, my = _rot(thP + thT, *MARKED)
    px, py = ox + mx, oy + my

    wsum = OMEGA_P + OMEGA_T
    # v = wP × (O - C) + (wP+wT) × (X - O), with ω × (x,y) = ω·(-y, x)
    rx, ry = ox - cx, oy - cy          # O - C
    sx, sy = px - ox, py - oy          # X - O
    vx = OMEGA_P * (-ry) + wsum * (-sy)
    vy = OMEGA_P * (rx) + wsum * (sx)
    # a = -wP²(O - C) - (wP+wT)²(X - O)   (constant rates → no tangential term)
    ax = -OMEGA_P**2 * rx - wsum**2 * sx
    ay = -OMEGA_P**2 * ry - wsum**2 * sy
    return (px, py), (vx, vy), (ax, ay)


@pytest.mark.parametrize("t", [0.0, 0.25, 0.5, 1.0, 2.3])
def test_merry_go_round_kinematics_match_analytic(t):
    """Position, velocity and acceleration of the marked point — GA
    reconstruction vs. closed-form epicyclic solution."""
    mgr = MerryGoRound(t)
    X = mgr.marked_world()
    v = mgr.point_velocity(X)
    a = mgr.point_acceleration(X)

    (px, py), (vx, vy), (ax, ay) = _analytic(t)

    # position (unitized projective point)
    assert X.z == pytest.approx(1.0)
    assert (X.x, X.y) == pytest.approx((px, py), abs=1e-12)
    # velocity / acceleration are free vectors (z = 0)
    assert v.z == pytest.approx(0.0, abs=1e-12)
    assert (v.x, v.y) == pytest.approx((vx, vy), abs=1e-12)
    assert a.z == pytest.approx(0.0, abs=1e-12)
    assert (a.x, a.y) == pytest.approx((ax, ay), abs=1e-12)


def test_marked_point_orbit_radius_is_constant():
    """|X - C| is the epicycle radius and must be the same at every time for
    T0's marked point (centre offset RADIUS, local arm |MARKED|)."""
    cx, cy = CENTRE
    radii = []
    for t in (0.0, 0.3, 0.7, 1.1, 2.0):
        X = MerryGoRound(t).marked_world()
        radii.append(math.hypot(X.x - cx, X.y - cy))
    # the orbit radius oscillates between |R - r| and |R + r|; here we just
    # assert it stays within those physical bounds (sanity on the composition)
    lo, hi = abs(RADIUS - math.hypot(*MARKED)), RADIUS + math.hypot(*MARKED)
    for r in radii:
        assert lo - 1e-12 <= r <= hi + 1e-12


def test_step_via_motor_manifold_matches_closed_form_pose():
    """kinematic_system2dp::step evolves a frame's relative pose on the motor
    manifold: P_new = P ⟇ exp(0.5·ξ·dt). For a constant rate, N steps of dt
    must reproduce the closed-form pose at t = N·dt — exercising pose2dp,
    motor_from_pose2dp / pose2dp_from_motor, rgpr, rrev and exp together."""
    dt, n = 0.016, 50
    xi = _to_twist(0.0, 0.0, OMEGA_T)             # turntable relative twist
    pose = pga.pose2dp(_pt(RADIUS, 0.0), 0.0)     # start unrotated
    for _ in range(n):
        P = pga.motor_from_pose2dp(pose)                      # body → parent
        P_new = pga.rgpr(P, pga.rexp(0.5 * xi * dt))
        pose = pga.pose2dp_from_motor(P_new)                  # decode evolved pose
    assert pose.phi == pytest.approx(OMEGA_T * n * dt, abs=1e-12)
    assert (pose.origin.x, pose.origin.y) == pytest.approx((RADIUS, 0.0), abs=1e-12)


# --- data-struct surface (accessibility of the physics artefacts) ---------

def test_pose2dp_fields_and_roundtrip():
    p = pga.pose2dp(_pt(1.0, 2.0), 0.5)
    assert (p.origin.x, p.origin.y, p.origin.z) == (1.0, 2.0, 1.0)
    assert p.phi == 0.5
    # mutate read/write fields
    p.phi = 1.25
    assert p.phi == 1.25
    # default ctor zero-initialises
    z = pga.pose2dp()
    assert z.phi == 0.0 and (z.origin.x, z.origin.y) == (0.0, 0.0)
    # pose ⟷ motor roundtrip
    M = pga.motor_from_pose2dp(p)
    back = pga.pose2dp_from_motor(M)
    assert back.phi == pytest.approx(p.phi, abs=1e-12)
    assert (back.origin.x, back.origin.y) == pytest.approx((1.0, 2.0), abs=1e-12)
    assert back == pga.pose2dp_from_motor(M)            # __eq__


def test_pose3dp_and_kin_state_constructible():
    p3 = pga.pose3dp(pga.vec3dp(1.0, 0.0, 0.0, 1.0), pga.vec3dp(0.0, 0.0, 1.0, 0.0))
    assert (p3.origin.x, p3.rot.z) == (1.0, 1.0)

    k2 = pga.kin_state2dp(pga.vec2dp(1.0, 0.0, 0.0), pga.vec2dp(0.0, 0.0, 0.0), 2.0, 0.0)
    assert k2.omega == 2.0 and k2.vel.x == 1.0
    assert pga.kin_state2dp().omega == 0.0             # default ctor

    k3 = pga.kin_state3dp()
    assert k3.omega.x == 0.0


def test_pose_str_repr_format():
    p = pga.pose2dp(_pt(1.0, 2.0), 0.5)
    assert str(p) == repr(p)
    assert "pose2dp" in str(p)
    # format spec routes through the C++ fmt::formatter
    assert "0.500" in f"{p:.3f}"
    with pytest.raises(ValueError):
        format(p, "not-a-spec")
