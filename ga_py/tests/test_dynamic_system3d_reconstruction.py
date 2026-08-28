"""Validate the 3dp dynamic reconstruction (ga_py.systems.DynamicSystem3dp) against the
two headline correctness metrics of the dynamic layer, generic (no application content):

  * a torque-free free body conserves BOTH kinetic energy and world angular momentum
    (the Poinsot / Dzhanibekov check) -- exercising step_free_body, compute_omega_dot,
    the inertia map and momentum_world;
  * a coupled articulated chain (a planar double pendulum, z-axis hinges) satisfies the
    mass-matrix / kinetic-energy identity  1/2 q̇ᵀ M(q) q̇ == KE, and conserves total
    energy across the swing -- the 3D lift of test_dynamic_system_reconstruction.py.

Runnable under pytest, or standalone:
    PYTHONPATH=build/external/ga/ga_py:external/ga/ga_py/python \
        python3 external/ga/ga_py/tests/test_dynamic_system3d_reconstruction.py
"""
from __future__ import annotations

import math

import ga_py
from ga_py import pga
from ga_py.systems import DynamicSystem3dp, make_cuboid_body


def _kin_state(omega: "pga.vec3dp") -> "pga.kin_state3dp":
    z = pga.vec3dp(0.0, 0.0, 0.0, 0.0)
    return pga.kin_state3dp(z, z, omega, z)  # (vel, acc, omega, alpha)


def _bivec_comps(b: "pga.bivec3dp") -> tuple:
    """The six line components, read by name exactly like vec3dp.x/.y/.z."""
    return (b.vx, b.vy, b.vz, b.mx, b.my, b.mz)


def _bivec_dist(a: "pga.bivec3dp", b: "pga.bivec3dp") -> float:
    ca, cb = _bivec_comps(a), _bivec_comps(b)
    return math.sqrt(sum((ca[i] - cb[i]) ** 2 for i in range(6)))


# --- torque-free tumble (free body) ---------------------------------------------------

def _make_tumbler(wx: float, wy: float, wz: float) -> DynamicSystem3dp:
    s = DynamicSystem3dp()
    s.set_gravity(pga.vec3dp(0.0, 0.0, 0.0, 0.0))  # torque-free
    s.add_frame("W")                                # inertial root
    # an asymmetric top: distinct extents -> three distinct principal moments
    s.add_body("top", make_cuboid_body(1.0, 1.0, 2.0, 3.0),
               state=_kin_state(pga.vec3dp(wx, wy, wz, 0.0)), parent_idx=0)
    return s


def test_free_body_conserves_energy_and_momentum():
    # spin mostly about the intermediate axis + small perturbation (Dzhanibekov regime)
    s = _make_tumbler(0.05, 2.0, 0.03)
    e0 = s.kinetic_energy()                     # g = 0 -> total energy == KE
    l0 = s.momentum_world("top")                # world momentum bivector
    emin = emax = e0
    lerr = 0.0
    dt = 1.0e-3
    for _ in range(2000):                       # 2 s
        s.step(dt)
        e = s.kinetic_energy()
        emin, emax = min(emin, e), max(emax, e)
        lerr = max(lerr, _bivec_dist(s.momentum_world("top"), l0))
    escale = abs(e0)
    lscale = math.sqrt(sum(c * c for c in _bivec_comps(l0)))
    assert (emax - emin) / escale < 1e-6         # energy conserved
    assert lerr / lscale < 1e-6                  # angular momentum conserved


# --- coupled double pendulum (z-axis hinges) ------------------------------------------

def _make_double_pendulum(phi1, w1, phi2, w2) -> DynamicSystem3dp:
    """Two thin square plates hinged about parallel z-axes (planar motion computed via the
    full 3D dynamics). Mirrors the 2D double-pendulum test: plate hinged to the world at
    its corner Q = (1, 1); plate 2 hinged to plate 1's opposite corner; rest cm poses
    (0, 0, 0) and (-2, -2, 0) relative to the parent.
    """
    body = make_cuboid_body(1.0, 2.0, 2.0, 0.1)  # thin square plate, m = 1, 2x2
    q = pga.vec3dp(1.0, 1.0, 0.0, 1.0)           # body-frame hinge point
    z = pga.vec3dp(0.0, 0.0, 1.0, 0.0)           # hinge axis
    zero_rot = pga.vec3dp(0.0, 0.0, 0.0, 0.0)
    s = DynamicSystem3dp()
    s.add_frame("W")
    s.add_revolute_body("p1", body, q, z, phi1, w1,
                        pose=pga.pose3dp(pga.vec3dp(0.0, 0.0, 0.0, 1.0), zero_rot),
                        parent_idx=0)
    s.add_revolute_body("p2", body, q, z, phi2, w2,
                        pose=pga.pose3dp(pga.vec3dp(-2.0, -2.0, 0.0, 1.0), zero_rot),
                        parent_idx=s.index_of("p1"))
    return s


def _states():
    return [(2.0, 1.0, -1.5, -2.0), (0.4, -0.7, 1.1, 0.5),
            (-1.2, 2.3, 0.8, -1.4), (0.0, 0.0, 0.0, 0.0)]


def test_double_pendulum_mass_matrix_matches_ke():
    for phi1, w1, phi2, w2 in _states():
        s = _make_double_pendulum(phi1, w1, phi2, w2)
        M = s.mass_matrix()  # 2x2 row-major
        qd = [w1, w2]
        ke_M = 0.5 * sum(qd[a] * M[a * 2 + b] * qd[b] for a in range(2) for b in range(2))
        assert abs(M[1] - M[2]) < 1e-12                   # symmetric
        assert abs(ke_M - s.kinetic_energy()) < 1e-10     # 1/2 q̇ᵀ M q̇ == KE


def test_double_pendulum_energy_conserved():
    s = _make_double_pendulum(2.0, 1.0, -1.5, -2.0)
    e0 = s.total_energy()
    emin = emax = e0
    kemax = s.kinetic_energy()
    dt = 5.0e-4
    for _ in range(4000):  # 2 s
        s.step(dt)
        e = s.total_energy()
        emin, emax = min(emin, e), max(emax, e)
        kemax = max(kemax, s.kinetic_energy())
    scale = max(kemax, abs(e0))
    assert (emax - emin) / scale < 1e-4


# --- standalone runner (no pytest needed) ---------------------------------------------

if __name__ == "__main__":
    checks = [
        ("free-body energy + momentum conservation (Dzhanibekov)",
         test_free_body_conserves_energy_and_momentum),
        ("double pendulum: 1/2 q̇ᵀ M q̇ == KE identity",
         test_double_pendulum_mass_matrix_matches_ke),
        ("double pendulum: total energy conserved over the swing",
         test_double_pendulum_energy_conserved),
    ]
    ok = True
    print("\nDynamicSystem3dp reconstruction — Module B validation\n")
    for name, fn in checks:
        try:
            fn()
            print(f"  [PASS] {name}")
        except AssertionError as e:
            ok = False
            print(f"  [FAIL] {name}: {e}")
    print(f"\n{'ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED'}\n")
    raise SystemExit(0 if ok else 1)


# --- an actuator torque at a joint ----------------------------------------------------

def test_joint_torque_pendulum_equilibrium():
    """set_joint_torque mirror: a constant torque against gravity on a damped pendulum
    settles at m g l sin(q) = tau; at rest the angular acceleration is tau / I_pivot."""
    m, l, g, w = 2.0, 0.5, 9.81, 0.05
    body = make_cuboid_body(m, w, 2.0 * l, w)
    zero_rot = pga.vec3dp(0.0, 0.0, 0.0, 0.0)
    z = pga.vec3dp(0.0, 0.0, 1.0, 0.0)

    def make(gg):
        s = DynamicSystem3dp()
        s.set_gravity(pga.vec3dp(0.0, -gg, 0.0, 0.0))
        s.add_frame("W")
        s.add_revolute_body("P", body, pga.vec3dp(0.0, l, 0.0, 1.0), z,
                            pose=pga.pose3dp(pga.vec3dp(0.0, -l, 0.0, 1.0), zero_rot),
                            parent_idx=0)
        return s

    i_piv = m * (w * w + 4.0 * l * l) / 12.0 + m * l * l
    s = make(0.0)
    s.set_joint_torque("P", lambda t: 0.7)
    assert abs(s.joint_accel("P") - 0.7 / i_piv) < 1e-12

    s = make(g)
    tau = 0.3 * m * g * l
    s.set_joint_torque("P", lambda t: tau)
    s.set_joint_spring_damper("P", 0.0, 2.0)
    for _ in range(8000):
        s.step(2.5e-3)
    assert abs(s.joint_phi("P") - math.asin(0.3)) < 1e-6
