"""Reconstruct the joint-space forward dynamics of dynamic_system2dp in Python,
from the bound GA primitives -- the dynamic-layer analogue of
test_merry_go_round.py (which reconstructs the kinematic layer).

The stateful dynamic_system2dp class is NOT bound; this test rebuilds its
coupled forward dynamics for a double pendulum out of the primitives that ARE
bound, exercising the joint descriptors the closed-loop work added to ga_py
(joint2dp / joint_state2dp) together with the inertia map, the spatial pairing
(rwdg), the motor algebra (rgpr / exp / move2dp / rrev) and the shared solvers
(lu_solve).

Mirrors dynamic_system2dp exactly:
  * a joint is a joint_state2dp {type, screw_b, rest, phi, omega}; its
    body->parent motor is  M(q) = rest (x) exp(0.5 q screw_b)  (rgpr + exp);
  * world motors compose with rgpr, parent on the left (= get_pos_trafo(i, 0));
  * the joint-space mass matrix is  M_jk = sum_i <S_j^(i), I_i[S_k^(i)]>  with
    the spatial pairing  <a, P> = -rwdg(a, P)  and the screws transported into
    each body frame by the adjoint sandwich;
  * the right-hand side is gravity minus the velocity-product (Coriolis) bias,
    read from the same world velocity/acceleration recursion (world_VA);
  * the state (phi, omega) per joint is advanced by RK4.

Validated by the two checks the C++ suite uses: the mass-matrix / kinetic-energy
identity  1/2 qdot^T M(q) qdot == KE, and energy conservation across the swing.
"""
from __future__ import annotations

import math

import pytest

import ga_py
from ga_py import pga

GRAV = pga.vec2dp(0.0, -9.81, 0.0)  # uniform gravity (a direction, z = 0)


def _pt(x: float, y: float) -> "pga.vec2dp":
    return pga.vec2dp(x, y, 1.0)


def _spatial_dot(xi: "pga.vec2dp", mom: "pga.bivec2dp") -> float:
    """Reciprocal (spatial) pairing twist x wrench -> scalar: <xi, mom> =
    -rwdg(xi, mom)  (dynamic_system2dp::spatial_dot)."""
    return -float(pga.rwdg(xi, mom))


class DoublePendulum:
    """Two revolute-jointed square plates, reconstructed from primitives.

    Geometry matches the C++ double-pendulum test: square plates (w = h = 2,
    m = 1); plate 1 hinged to the world at its corner Q = (1, 1); plate 2 hinged
    to plate 1's opposite corner, its own corner also Q = (1, 1). Rest poses
    (relative to the parent): cm1 at (0, 0), cm2 at (-2, -2).
    """

    M = 1.0
    W = H = 2.0
    Q = _pt(1.0, 1.0)  # the body-frame hinge / revolute screw generator

    def __init__(self, phi1: float, omega1: float, phi2: float, omega2: float):
        I = pga.get_plate_inertia(self.M, self.W, self.H)  # about each cm
        # rest = body->parent motor at q = 0 (motor_from_pose2dp of the rest pose)
        rest1 = pga.motor_from_pose2dp(pga.pose2dp(_pt(0.0, 0.0), 0.0))   # identity
        rest2 = pga.motor_from_pose2dp(pga.pose2dp(_pt(-2.0, -2.0), 0.0)) # translation

        # the joint descriptors the closed-loop work unlocked in ga_py
        self.joints = [
            pga.joint_state2dp(pga.joint2dp.revolute, self.Q, rest1, phi1, omega1),
            pga.joint_state2dp(pga.joint2dp.revolute, self.Q, rest2, phi2, omega2),
        ]
        self.inertia = [I, I]
        self._refresh()

    # --- kinematics rebuilt from the joint state -----------------------------

    def _joint_motor(self, j: "pga.joint_state2dp") -> "pga.mvec2dp_u":
        # M(q) = rest (x) exp(0.5 q screw_b)
        return pga.rgpr(j.rest, pga.exp(0.5 * j.phi * j.screw_b))

    def _refresh(self) -> None:
        # world (body i -> world) motors: parent applied on the left via rgpr
        m0 = self._joint_motor(self.joints[0])
        m1 = pga.rgpr(m0, self._joint_motor(self.joints[1]))
        self.Mw = [m0, m1]

    def _world_VA(self, upto: int) -> tuple["pga.vec2dp", "pga.vec2dp"]:
        """World velocity & acceleration twists of body `upto`, via the
        Newton-Euler recursion with the relative acceleration twists zeroed
        (the q-ddot = 0 velocity-product bias) -- as in world_VA."""
        V = pga.vec2dp(0.0, 0.0, 0.0)
        A = pga.vec2dp(0.0, 0.0, 0.0)
        for k in range(upto + 1):
            xi = self.joints[k].omega * self.joints[k].screw_b   # rel twist = q-dot screw
            zeta = pga.move2dp(xi, self.Mw[k])                   # Ad_M(xi) into world
            V = V + zeta
            A = A + pga.rcmt(V, zeta)                            # zetadot = 0; Coriolis
        return V, A

    def _world_screw(self, j: int) -> "pga.vec2dp":
        return pga.move2dp(self.joints[j].screw_b, self.Mw[j])

    def _cm(self, i: int) -> "pga.vec2dp":
        return pga.move2dp(_pt(0.0, 0.0), self.Mw[i])

    # --- dynamic-layer quantities -------------------------------------------

    def kinetic_energy(self) -> float:
        ke = 0.0
        for i in range(2):
            Vb = pga.move2dp(self._world_VA(i)[0], pga.rrev(self.Mw[i]))  # body twist
            ke += 0.5 * _spatial_dot(Vb, self.inertia[i](Vb))
        return ke

    def potential_energy(self) -> float:
        pe = 0.0
        for i in range(2):
            cm = self._cm(i)
            pe += -self.M * (GRAV.x * cm.x + GRAV.y * cm.y)
        return pe

    def total_energy(self) -> float:
        return self.kinetic_energy() + self.potential_energy()

    def mass_matrix(self) -> list[list[float]]:
        n = len(self.joints)
        S = [self._world_screw(j) for j in range(n)]
        Mmat = [[0.0] * n for _ in range(n)]
        for i in range(n):  # body i has joints 0..i as ancestors (a chain)
            Iinv_to_body = pga.rrev(self.Mw[i])
            I = self.inertia[i]
            for j in range(i + 1):
                xj = pga.move2dp(S[j], Iinv_to_body)  # joint-j screw in body i frame
                for k in range(i + 1):
                    xk = pga.move2dp(S[k], Iinv_to_body)
                    Mmat[j][k] += _spatial_dot(xj, I(xk))
        return Mmat

    def rhs(self) -> list[float]:
        """Generalised force tau_j = sum_i m_i ( vcm_i(S_j).g - vcm_i(S_j).b_cm_i ),
        the gravity force minus the velocity-product (Coriolis) bias."""
        n = len(self.joints)
        S = [self._world_screw(j) for j in range(n)]
        rhs = [0.0] * n
        for i in range(n):
            cm = self._cm(i)
            V, A = self._world_VA(i)
            b = pga.rcmt(A, cm) + pga.rcmt(V, pga.rcmt(V, cm))  # accel field, q-ddot = 0
            for j in range(i + 1):
                vj = pga.rcmt(S[j], cm)  # cm velocity per unit rate of joint j
                rhs[j] += self.M * (vj.x * GRAV.x + vj.y * GRAV.y)
                rhs[j] -= self.M * (vj.x * b.x + vj.y * b.y)
        return rhs

    def joint_accel(self) -> list[float]:
        n = len(self.joints)
        M = self.mass_matrix()
        flat = [M[a][b] for a in range(n) for b in range(n)]
        return ga_py.lu_solve(flat, self.rhs(), n)

    def step(self, dt: float) -> None:
        """RK4 on the joint state u = [phi0, phi1, omega0, omega1]."""
        n = len(self.joints)

        def deriv(u):
            for k in range(n):
                self.joints[k].phi = u[k]
                self.joints[k].omega = u[n + k]
            self._refresh()
            qdd = self.joint_accel()
            return [u[n + k] for k in range(n)] + qdd  # (dphi, domega) = (omega, qddot)

        u0 = [self.joints[k].phi for k in range(n)] + [
            self.joints[k].omega for k in range(n)
        ]
        k1 = deriv(u0)
        k2 = deriv([u0[i] + 0.5 * dt * k1[i] for i in range(2 * n)])
        k3 = deriv([u0[i] + 0.5 * dt * k2[i] for i in range(2 * n)])
        k4 = deriv([u0[i] + dt * k3[i] for i in range(2 * n)])
        u = [u0[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(2 * n)]
        for k in range(n):
            self.joints[k].phi = u[k]
            self.joints[k].omega = u[n + k]
        self._refresh()


# --- tests ----------------------------------------------------------------

def test_joint_state_descriptor_drives_the_motor():
    """The joint_state2dp descriptor (type / screw_b / rest / phi / omega) is
    what the reconstruction reads -- a smoke check that it round-trips."""
    dp = DoublePendulum(0.7, 0.0, -0.3, 0.0)
    assert dp.joints[0].type == pga.joint2dp.revolute
    assert dp.joints[0].screw_b == pga.vec2dp(1.0, 1.0, 1.0)
    assert dp.joints[0].phi == pytest.approx(0.7)


@pytest.mark.parametrize(
    "phi1,w1,phi2,w2",
    [
        (2.0, 1.0, -1.5, -2.0),  # the C++ test's energetic initial state
        (0.4, -0.7, 1.1, 0.5),
        (-1.2, 2.3, 0.8, -1.4),
        (0.0, 0.0, 0.0, 0.0),
    ],
)
def test_mass_matrix_matches_kinetic_energy(phi1, w1, phi2, w2):
    """The built-in consistency check: 1/2 qdot^T M(q) qdot == KE, the same
    identity dynamic_system2dp's tests assert. Validates the reconstructed
    mass matrix, world twists and inertia-map pairing."""
    dp = DoublePendulum(phi1, w1, phi2, w2)
    M = dp.mass_matrix()
    qd = [w1, w2]
    ke_M = 0.5 * sum(qd[a] * M[a][b] * qd[b] for a in range(2) for b in range(2))
    assert M[0][1] == pytest.approx(M[1][0], abs=1e-12)  # symmetric
    assert ke_M == pytest.approx(dp.kinetic_energy(), rel=1e-12, abs=1e-12)


def test_energy_conserved_over_the_swing():
    """Constrained-free coupled dynamics conserve energy to integrator
    tolerance -- the headline correctness metric of the dynamic layer,
    reproduced here from the bound primitives."""
    dp = DoublePendulum(2.0, 1.0, -1.5, -2.0)
    e0 = dp.total_energy()
    emin = emax = e0
    kemax = dp.kinetic_energy()
    dt = 0.0005
    for _ in range(4000):  # 2 s
        dp.step(dt)
        e = dp.total_energy()
        emin = min(emin, e)
        emax = max(emax, e)
        kemax = max(kemax, dp.kinetic_energy())
    scale = max(kemax, abs(e0))
    assert (emax - emin) / scale < 1e-4
