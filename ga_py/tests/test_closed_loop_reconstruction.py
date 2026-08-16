"""Reconstruct the closed-loop solver of closed_loop_system2dp in Python, from
the bound GA primitives -- the closed-loop analogue of
test_dynamic_system_reconstruction.py (open-chain dynamics) and
test_merry_go_round.py (kinematics).

The stateful closed_loop_system2dp class is NOT bound (like every system class);
this test rebuilds its three solve tiers for a planar FOUR-BAR linkage out of the
primitives that ARE bound -- the joint descriptors (joint2dp / joint_state2dp),
the loop-closure descriptor (loop_constraint2dp / constraint2dp), the motor
algebra (rgpr / exp / move2dp / rrev), the velocity bracket (rcmt), the spatial
pairing (rwdg), the inertia map, and the shared solvers (lu_solve / lstsq_solve).

Mirrors closed_loop_system2dp exactly:
  * a spanning TREE of revolute bodies (here a genuine branching tree: crank ->
    coupler is one branch, the rocker a second branch from the ground), each body
    a joint_state2dp whose body->parent motor is  M(q) = rest (x) exp(0.5 q S_b);
  * the closure residual g(q) = P_a - P_b for one point-coincidence, and its
    Jacobian G whose columns are the relative anchor partial velocities
    rcmt(S_j, P) over the ancestor joints (= constraint_jacobian);
  * POSITION: assemble() Newton-solves g(q) = 0  (lstsq_solve on G_dep);
  * VELOCITY: solve_velocities() distributes a driver rate via  G_dep qdot_dep =
    -G_drv qdot_drv;
  * DYNAMICS: the bordered acceleration-level KKT system
        | M  Gᵀ | |qddot |   | tau        |
        | G  0  | |lambda| = | -Gdot qdot |
    (M, tau from the open-loop assembly; -Gdot qdot the anchor velocity-product
    bias) solved through the bound lu_solve, integrated by RK4 with post-step GGL
    projection of position (assemble) and velocity (qdot <- qdot - G⁺(G qdot)).

Validated by the same checks the C++ four-bar tests use: the assembled
configuration matches the analytic four-bar angles (Phase 1), the distributed
velocities satisfy the velocity-level closure G qdot = 0 (Phase 2), and the
released mechanism conserves energy while ‖g‖ stays bounded (Phase 3).
"""
from __future__ import annotations

import math
from typing import NamedTuple

import pytest

import ga_py
from ga_py import pga

GRAV = pga.vec2dp(0.0, -9.81, 0.0)  # uniform gravity (a direction, z = 0)

_ORIGIN_SCREW = pga.vec2dp(0.0, 0.0, 1.0)  # revolute screw: rotation about the frame origin


def _pt(x: float, y: float) -> "pga.vec2dp":
    return pga.vec2dp(x, y, 1.0)


def _spatial_dot(xi: "pga.vec2dp", mom: "pga.bivec2dp") -> float:
    """Reciprocal (spatial) pairing twist x wrench -> scalar (spatial_dot)."""
    return -float(pga.rwdg(xi, mom))


class Loop(NamedTuple):
    """A point-coincidence loop closure: anchor_a (in frame_a) == anchor_b (in
    frame_b), in world coordinates (loop_constraint2dp, constraint2dp.coincidence)."""
    frame_a: int
    anchor_a: "pga.vec2dp"
    frame_b: int
    anchor_b: "pga.vec2dp"


class ClosedLoopSystem:
    """closed_loop_system2dp reconstructed from primitives: a spanning tree of
    revolute bodies + point-coincidence loop closures.

    Each body i carries: a parent index (None = the world root), a joint_state2dp
    descriptor, a mass, and an inertia map (about the cm = the frame origin)."""

    def __init__(self):
        self.parent: list[int | None] = []
        self.joints: list["pga.joint_state2dp"] = []
        self.mass: list[float] = []
        self.inertia: list = []
        self.loops: list[Loop] = []
        self.Mw: list["pga.mvec2dp_u"] = []
        self._identity = pga.motor_from_pose2dp(pga.pose2dp(_pt(0.0, 0.0), 0.0))

    # --- build the spanning tree -------------------------------------------------

    def add_revolute_body(self, parent, rest_pose, mass, w, h, phi0=0.0, omega0=0.0) -> int:
        rest = pga.motor_from_pose2dp(rest_pose)
        self.parent.append(parent)
        self.joints.append(
            # trailing 0.0s = no spring/damper force elements (stiffness, damping, q_rest)
            pga.joint_state2dp(pga.joint2dp.revolute, _ORIGIN_SCREW, rest, phi0, omega0,
                               0.0, 0.0, 0.0)
        )
        self.mass.append(mass)
        self.inertia.append(pga.get_plate_inertia(mass, w, h))
        self._refresh()
        return len(self.joints) - 1

    def add_loop(self, loop: Loop) -> None:
        self.loops.append(loop)

    @property
    def dof(self) -> list[int]:
        return list(range(len(self.joints)))

    # --- kinematics rebuilt from the joint state ---------------------------------

    def _joint_motor(self, j: "pga.joint_state2dp") -> "pga.mvec2dp_u":
        # M(q) = rest (x) exp(0.5 q screw_b)
        return pga.rgpr(j.rest, pga.rexp(0.5 * j.phi * j.screw_b))

    def _refresh(self) -> None:
        """World (body i -> world) motors: parent applied on the left via rgpr.
        Parents are added before children, so a single forward pass suffices."""
        self.Mw = [self._identity] * len(self.joints)
        for i, j in enumerate(self.joints):
            par = self.parent[i]
            base = self._identity if par is None else self.Mw[par]
            self.Mw[i] = pga.rgpr(base, self._joint_motor(j))

    def _anc_path(self, i: int) -> list[int]:
        """Ancestor joints of body i, root-first (the path root..i)."""
        path = []
        k: int | None = i
        while k is not None:
            path.append(k)
            k = self.parent[k]
        path.reverse()
        return path

    def _world_screw(self, j: int) -> "pga.vec2dp":
        return pga.move2dp(self.joints[j].screw_b, self.Mw[j])

    def _cm(self, i: int) -> "pga.vec2dp":
        return pga.move2dp(_pt(0.0, 0.0), self.Mw[i])

    def _world_VA(self, frame: int) -> tuple["pga.vec2dp", "pga.vec2dp"]:
        """World velocity & acceleration twists of `frame`, via the Newton-Euler
        recursion along its ancestor path with the relative accel twists zeroed
        (the q-ddot = 0 velocity-product bias) -- as in world_VA."""
        V = pga.vec2dp(0.0, 0.0, 0.0)
        A = pga.vec2dp(0.0, 0.0, 0.0)
        for k in self._anc_path(frame):
            xi = self.joints[k].omega * self.joints[k].screw_b  # rel twist = q-dot screw
            zeta = pga.move2dp(xi, self.Mw[k])                  # Ad_M(xi) into world
            V = V + zeta
            A = A + pga.rcmt(V, zeta)                           # zetadot = 0; Coriolis
        return V, A

    def _point_accel(self, P: "pga.vec2dp", frame: int) -> "pga.vec2dp":
        """Acceleration of point P rigidly attached to `frame`, at q-ddot = 0
        (point_acceleration): rcmt(A, P) + rcmt(V, rcmt(V, P))."""
        V, A = self._world_VA(frame)
        return pga.rcmt(A, P) + pga.rcmt(V, pga.rcmt(V, P))

    # --- closure residual & Jacobian ---------------------------------------------

    def _anchor_world(self, frame: int, anchor: "pga.vec2dp") -> "pga.vec2dp":
        P = pga.move2dp(anchor, self.Mw[frame])
        return pga.vec2dp(P.x / P.z, P.y / P.z, 1.0)

    def residual(self) -> list[float]:
        g = []
        for lc in self.loops:
            Pa = self._anchor_world(lc.frame_a, lc.anchor_a)
            Pb = self._anchor_world(lc.frame_b, lc.anchor_b)
            g += [Pa.x - Pb.x, Pa.y - Pb.y]
        return g

    def residual_norm(self) -> float:
        g = self.residual()
        return max((abs(x) for x in g), default=0.0)

    def constraint_jacobian(self, cols: list[int]) -> list[list[float]]:
        """G (m x len(cols), m = 2 * loop_count): column k for constraint c is the
        relative anchor partial velocity rcmt(S_j, P_a) - rcmt(S_j, P_b), each term
        present only when joint cols[k] is an ancestor of that anchor's frame."""
        m = 2 * len(self.loops)
        G = [[0.0] * len(cols) for _ in range(m)]
        Sw = [self._world_screw(j) for j in cols]
        for c, lc in enumerate(self.loops):
            Pa = self._anchor_world(lc.frame_a, lc.anchor_a)
            Pb = self._anchor_world(lc.frame_b, lc.anchor_b)
            anc_a = set(self._anc_path(lc.frame_a))
            anc_b = set(self._anc_path(lc.frame_b))
            for k, j in enumerate(cols):
                dPa = pga.rcmt(Sw[k], Pa) if j in anc_a else pga.vec2dp(0.0, 0.0, 0.0)
                dPb = pga.rcmt(Sw[k], Pb) if j in anc_b else pga.vec2dp(0.0, 0.0, 0.0)
                G[2 * c + 0][k] = dPa.x - dPb.x
                G[2 * c + 1][k] = dPa.y - dPb.y
        return G

    # --- POSITION: assemble (Newton on g = 0) ------------------------------------

    def assemble(self, driven: tuple[int, ...] = (), tol: float = 1e-12,
                 max_iter: int = 50) -> float:
        dep = [j for j in self.dof if j not in driven]
        for it in range(max_iter + 1):
            g = self.residual()
            gnorm = max((abs(x) for x in g), default=0.0)
            if gnorm < tol:
                return gnorm
            if it == max_iter:
                break
            G = self.constraint_jacobian(dep)
            flat = [G[i][k] for i in range(len(g)) for k in range(len(dep))]
            b = [-x for x in g]
            delta = ga_py.lstsq_solve(flat, b, len(dep))
            for k, j in enumerate(dep):
                self.joints[j].phi += delta[k]
            self._refresh()
        raise RuntimeError(f"assemble: Newton did not converge (residual {self.residual_norm()})")

    # --- VELOCITY: distribute driver rate (G qdot = 0) ---------------------------

    def solve_velocities(self, driven: tuple[int, ...]) -> list[float]:
        dep = [j for j in self.dof if j not in driven]
        m = 2 * len(self.loops)
        Gdrv = self.constraint_jacobian(list(driven))
        rhs = [0.0] * m
        for i in range(m):
            for k, j in enumerate(driven):
                rhs[i] -= Gdrv[i][k] * self.joints[j].omega
        Gdep = self.constraint_jacobian(dep)
        flat = [Gdep[i][k] for i in range(m) for k in range(len(dep))]
        qdot_dep = ga_py.lstsq_solve(flat, rhs, len(dep))
        for k, j in enumerate(dep):
            self.joints[j].omega = qdot_dep[k]
        self._refresh()
        return qdot_dep

    def project_velocities(self) -> None:
        """Velocity projection qdot <- qdot - G⁺(G qdot) (GGL stabilisation)."""
        rj = self.dof
        n = len(rj)
        m = 2 * len(self.loops)
        G = self.constraint_jacobian(rj)
        qd = [self.joints[j].omega for j in rj]
        Gv = [sum(G[i][k] * qd[k] for k in range(n)) for i in range(m)]
        flat = [G[i][k] for i in range(m) for k in range(n)]
        dqd = ga_py.lstsq_solve(flat, Gv, n)
        for k, j in enumerate(rj):
            self.joints[j].omega -= dqd[k]
        self._refresh()

    # --- DYNAMICS: open-loop M / tau + the bordered KKT solve --------------------

    def _mass_matrix(self, rj: list[int]) -> list[list[float]]:
        n = len(rj)
        Sw = {j: self._world_screw(j) for j in rj}
        M = [[0.0] * n for _ in range(n)]
        for bi in self.dof:                     # bodies and joints are 1:1
            anc = set(self._anc_path(bi))
            Iinv_to_body = pga.rrev(self.Mw[bi])
            I = self.inertia[bi]
            xb = {j: pga.move2dp(Sw[j], Iinv_to_body) for j in rj if j in anc}
            for a, ja in enumerate(rj):
                if ja not in anc:
                    continue
                for b, jb in enumerate(rj):
                    if jb not in anc:
                        continue
                    M[a][b] += _spatial_dot(xb[ja], I(xb[jb]))
        return M

    def _generalised_forces(self, rj: list[int]) -> list[float]:
        """tau_j = sum_i m_i ( vcm_i(S_j).g - vcm_i(S_j).b_cm_i ): gravity force
        minus the velocity-product (Coriolis) bias, summed over the bodies that
        have joint j as an ancestor."""
        idx = {j: k for k, j in enumerate(rj)}
        Sw = {j: self._world_screw(j) for j in rj}
        tau = [0.0] * len(rj)
        for bi in self.dof:
            anc = set(self._anc_path(bi))
            cm = self._cm(bi)
            V, A = self._world_VA(bi)
            b = pga.rcmt(A, cm) + pga.rcmt(V, pga.rcmt(V, cm))  # accel field, q-ddot = 0
            m = self.mass[bi]
            for j in rj:
                if j not in anc:
                    continue
                vj = pga.rcmt(Sw[j], cm)  # cm velocity per unit rate of joint j
                tau[idx[j]] += m * (vj.x * GRAV.x + vj.y * GRAV.y)
                tau[idx[j]] -= m * (vj.x * b.x + vj.y * b.y)
        return tau

    def _constraint_bias(self) -> list[float]:
        """G-dot q-dot: the relative anchor acceleration evaluated at q-ddot = 0."""
        gd = []
        for lc in self.loops:
            Pa = self._anchor_world(lc.frame_a, lc.anchor_a)
            Pb = self._anchor_world(lc.frame_b, lc.anchor_b)
            aa = self._point_accel(Pa, lc.frame_a)
            ab = self._point_accel(Pb, lc.frame_b)
            gd += [aa.x - ab.x, aa.y - ab.y]
        return gd

    def joint_accelerations(self) -> list[float]:
        """Solve the bordered KKT system through the bound lu_solve, returning the
        joint accelerations (the lambda block is discarded here)."""
        rj = self.dof
        n = len(rj)
        m = 2 * len(self.loops)
        M = self._mass_matrix(rj)
        G = self.constraint_jacobian(rj)
        tau = self._generalised_forces(rj)
        gd = self._constraint_bias()

        N = n + m
        K = [[0.0] * N for _ in range(N)]
        for a in range(n):
            for b in range(n):
                K[a][b] = M[a][b]
        for c in range(m):
            for a in range(n):
                K[a][n + c] = G[c][a]  # Gᵀ
                K[n + c][a] = G[c][a]  # G
        r = [tau[a] for a in range(n)] + [-gd[c] for c in range(m)]
        flat = [K[a][b] for a in range(N) for b in range(N)]
        sol = ga_py.lu_solve(flat, r, N)
        return sol[:n]

    def step(self, dt: float) -> None:
        """RK4 on the coupled joint state, then GGL projection of (position,
        velocity) back onto the constraint manifold."""
        rj = self.dof
        n = len(rj)

        def setu(u):
            for k, j in enumerate(rj):
                self.joints[j].phi = u[k]
                self.joints[j].omega = u[n + k]
            self._refresh()

        def deriv(u):
            setu(u)
            qdd = self.joint_accelerations()
            return [u[n + k] for k in range(n)] + qdd

        u0 = [self.joints[j].phi for j in rj] + [self.joints[j].omega for j in rj]
        k1 = deriv(u0)
        k2 = deriv([u0[i] + 0.5 * dt * k1[i] for i in range(2 * n)])
        k3 = deriv([u0[i] + 0.5 * dt * k2[i] for i in range(2 * n)])
        k4 = deriv([u0[i] + dt * k3[i] for i in range(2 * n)])
        u = [u0[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(2 * n)]
        setu(u)

        self.assemble(driven=())   # position projection (min-norm Newton)
        self.project_velocities()  # velocity projection

    # --- energy ------------------------------------------------------------------

    def kinetic_energy(self) -> float:
        ke = 0.0
        for i in self.dof:
            V, _ = self._world_VA(i)
            Vb = pga.move2dp(V, pga.rrev(self.Mw[i]))  # body twist
            ke += 0.5 * _spatial_dot(Vb, self.inertia[i](Vb))
        return ke

    def potential_energy(self) -> float:
        pe = 0.0
        for i in self.dof:
            cm = self._cm(i)
            pe += -self.mass[i] * (GRAV.x * cm.x + GRAV.y * cm.y)
        return pe

    def total_energy(self) -> float:
        return self.kinetic_energy() + self.potential_energy()


# --- the planar four-bar (matches ga_appl2dp_mechanics_test.hpp) ------------------
#
#   a (crank) = 2, b (coupler) = 3, c (rocker) = sqrt(5), d (ground) = 4
#   O2 = (0,0), O4 = (d,0); spanning tree crank -> coupler + rocker branch, closed
#   by coupler tip (b,0 in CO) == rocker tip (c,0 in RO). Each link is a unit plate
#   whose cm sits at its frame origin (the joint).

_A, _B, _C, _D = 2.0, 3.0, math.sqrt(5.0), 4.0


def _four_bar(theta2: float, phi3_guess: float, phi4_guess: float) -> tuple[ClosedLoopSystem, int]:
    cl = ClosedLoopSystem()
    CR = cl.add_revolute_body(None, pga.pose2dp(_pt(0.0, 0.0), 0.0), 1.0, 1.0, 1.0, theta2)
    CO = cl.add_revolute_body(CR, pga.pose2dp(_pt(_A, 0.0), 0.0), 1.0, 1.0, 1.0, phi3_guess)
    RO = cl.add_revolute_body(None, pga.pose2dp(_pt(_D, 0.0), 0.0), 1.0, 1.0, 1.0, phi4_guess)
    cl.add_loop(Loop(CO, _pt(_B, 0.0), RO, _pt(_C, 0.0)))
    return cl, CR


# --- tests ----------------------------------------------------------------------

def test_loop_descriptor_round_trips():
    """A smoke check that the bound closure descriptor (loop_constraint2dp /
    constraint2dp) carries the data the reconstruction reads."""
    lc = pga.loop_constraint2dp(1, _pt(3.0, 0.0), 2, _pt(math.sqrt(5.0), 0.0),
                                pga.constraint2dp.coincidence)
    assert lc.frame_a == 1 and lc.frame_b == 2
    assert lc.anchor_a == _pt(3.0, 0.0)
    assert lc.type == pga.constraint2dp.coincidence


def test_assembles_to_analytic_four_bar_config():
    """Phase 1: drive the crank to pi/2 and Newton-solve the dependent coupler /
    rocker angles. At theta2 = pi/2: A = (0,2), B = (3,2) -> theta3 = -pi/2,
    theta4 = atan2(2,-1), the analytic closed configuration."""
    cl, CR = _four_bar(math.pi / 2.0, phi3_guess=-1.4, phi4_guess=2.2)

    assert cl.residual_norm() > 0.1  # the loop starts open
    gfin = cl.assemble(driven=(CR,))
    assert gfin < 1e-12
    assert cl.joints[1].phi == pytest.approx(-math.pi / 2.0, abs=1e-9)      # coupler
    assert cl.joints[2].phi == pytest.approx(math.atan2(2.0, -1.0), abs=1e-9)  # rocker
    assert cl.joints[0].phi == pytest.approx(math.pi / 2.0, abs=1e-15)      # driven crank


def test_velocity_distribution_closes_the_loop():
    """Phase 2: with the crank driven at a unit rate, the distributed coupler /
    rocker rates satisfy the velocity-level closure G qdot = 0."""
    cl, CR = _four_bar(1.0, phi3_guess=-1.4, phi4_guess=2.2)
    cl.assemble(driven=(CR,))
    cl.joints[CR].omega = 1.3  # prescribe the driver rate
    cl._refresh()
    cl.solve_velocities(driven=(CR,))

    # G qdot must vanish (the loop neither opens nor closes at velocity level)
    qd = [cl.joints[j].omega for j in cl.dof]
    G = cl.constraint_jacobian(cl.dof)
    for row in G:
        assert sum(row[k] * qd[k] for k in range(len(qd))) == pytest.approx(0.0, abs=1e-10)


def test_energy_conserved_over_the_swing():
    """Phase 3: released from rest under gravity, the four-bar's constraint forces
    do no work, so energy is conserved to integrator tolerance while the closure
    error ‖g‖ stays bounded -- reconstructed here entirely from bound primitives."""
    cl, CR = _four_bar(1.2, phi3_guess=-1.4, phi4_guess=2.2)
    cl.assemble(driven=(CR,))   # start from a consistent closed configuration
    assert cl.residual_norm() < 1e-12

    e0 = cl.total_energy()
    emin = emax = e0
    kemax = cl.kinetic_energy()
    gmax = 0.0
    dt = 0.001
    for _ in range(1000):  # 1 s
        cl.step(dt)
        e = cl.total_energy()
        emin = min(emin, e)
        emax = max(emax, e)
        kemax = max(kemax, cl.kinetic_energy())
        gmax = max(gmax, cl.residual_norm())

    scale = max(kemax, abs(e0))
    assert (emax - emin) / scale < 1e-4  # energy conserved (constraint forces do no work)
    assert gmax < 1e-9                   # the loop stays closed (projection controls drift)
    assert kemax > 0.05                  # the mechanism actually moved
