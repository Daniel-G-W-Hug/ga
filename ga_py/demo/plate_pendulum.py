"""PGA2DP demo: rigid 2D plate pendulum pinned at its top-right corner.

Python port of the test case at
``ga_test/src/ga_appl2dp_physics_test.hpp:1490-1768``::

    TEST_CASE("pga2dp: combined motion (rigid 1-body system)")

Side-by-side reading is intended — the demo mirrors the C++ source line
for line, swapping mdspan plumbing for the std::vector ``rk4_step``
overload bound for Python. See ``ga_py/README.md`` §3 for setup.

Run from the project root:

    PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python" \\
        ga_py/.venv/bin/python ga_py/demo/plate_pendulum.py

What this demonstrates
----------------------
The full physics workflow in projective GA from Python:

* body-frame motor formulation: ``M(t) = M0 ⟇ exp(½ B_b(t))``;
* inertia construction with parallel-axis (Steiner) correction
  (``get_plate_inertia`` with a non-default pivot);
* LU-based inertia inversion (``get_inertia_inverse``);
* body-frame Euler ODE rhs (``compute_omega_dot``);
* RK4 time integration via ``ga_py.rk4_step`` (4 sub-steps per time step).

Energy conservation and the pivot constraint are checked at every step.

Notes for Python users
----------------------
* ``rk4_step`` returns ``(u_new, uh_new)`` (a 2-tuple of lists) — the
  C++ overload mutates in place; Python takes by value and returns. See
  README §6.5 for the rationale.
* PGA's ``vec2dp`` exposes ``.x .y .z`` (z is the projective coord).
  Bivectors expose ``.x .y .z`` as the three e-products.
* GA types support f-string format specs (``f"{v:>-7.3f}"``); the spec
  is applied per component on the C++ side via ``fmt::format``.
"""
import math

import ga_py
from ga_py import pga


# ---------------------------------------------------------------------
# 1-body simulation class — direct mirror of C++ `sim_ode_plate_pga2dp`.
# ---------------------------------------------------------------------
class SimOdePlatePga2dp:
    """2D plate pendulum: rigid body pinned at its top-right corner Q_b.

    State (mirrors the C++ struct):
      u[0] = B_b   (position bivector;       B_b = phi   * Q_b)
      u[1] = Omega (velocity / rate bivector;Omega = omega * Q_b)
    Both are stored as `vec2dp` in PGA2DP (the dual of a "rate of change"
    bivector in 2dp lives in vec2dp, see ga/ga_pga2dp_ops_physics.hpp).
    """

    def __init__(self, m, w, h, cm_pos, cm_phi, cm_omega):
        self.m = m
        self.width = w
        self.height = h
        self.cm_w_pos0 = cm_pos
        self.cm_w_phi0 = cm_phi
        self.cm_w_omega0 = cm_omega

        # Inertia about the body-frame pivot Q_b = (hw, hh, 1) — the
        # parallel-axis correction is applied inside get_plate_inertia.
        hw, hh = w / 2.0, h / 2.0
        Q_b = pga.vec2dp(hw, hh, 1.0)
        print("inertia of plate in body system (about pivot Q_b = TR corner):")
        self.I = pga.get_plate_inertia(m, w, h, Q_b)
        self.I_inv = pga.get_inertia_inverse(self.I)
        print(f"system inertia (about pivot Q_b):")
        print(f"I     = {self.I:>-7.3f}")
        print(f"I_inv = {self.I_inv:>-7.3f}")

        # rk4 working state — initialised in set_initial_values().
        self.u = [pga.vec2dp(0, 0, 0), pga.vec2dp(0, 0, 0)]
        self.uh = ([pga.vec2dp(0, 0, 0), pga.vec2dp(0, 0, 0)],
                   [pga.vec2dp(0, 0, 0), pga.vec2dp(0, 0, 0)])
        self.M0 = None       # initial translation motor (mvec2dp_u)
        self.pivot_w = None  # world position of the pivot (vec2dp)

    def set_initial_values(self):
        hw, hh = self.width / 2.0, self.height / 2.0

        # M0: pure translation placing the body origin (= cm) at cm_w_pos0.
        # Encoding: translation by (tx, ty) → motor exp(0.5*(-ty, tx, 0)).
        self.M0 = pga.exp(0.5 * pga.vec2dp(
            -self.cm_w_pos0.y, self.cm_w_pos0.x, 0.0))

        # World position of the body-frame pivot Q_b — fixed throughout the
        # rotation: move2dp(Q_b, M0 ⟇ exp(½ phi Q_b)) = move2dp(Q_b, M0).
        Q_b = pga.vec2dp(hw, hh, 1.0)
        self.pivot_w = pga.move2dp(Q_b, self.M0)
        print(f"pivot_w = {self.pivot_w:>-7.3f}  "
              f"(world position of pivot Q_b via M0)")

        # B_b(0) = phi0 * Q_b ; Omega(0) = omega0 * Q_b (pivot constraint).
        self.u[0] = self.cm_w_phi0 * Q_b
        self.u[1] = self.cm_w_omega0 * Q_b

    def calc_rhs(self):
        """Body-frame ODE right-hand side.

        Force couple in world frame: gravity at cm + reaction at pivot.
        The reaction is split out so we can transform the *net torque*
        into the body frame for compute_omega_dot.
        """
        hw, hh = self.width / 2.0, self.height / 2.0
        Q_b = pga.vec2dp(hw, hh, 1.0)

        B = self.u[0]      # B_b = phi * Q_b
        Omega = self.u[1]  # Omega_b = omega * Q_b

        # Current motor M(t) = M0 ⟇ exp(½ B_b(t)) and the cm world position.
        M = pga.rgpr(self.M0, pga.exp(0.5 * B))
        cm_w = pga.move2dp(pga.vec2dp(0.0, 0.0, 1.0), M)

        # Force couple. Gravity g acts at cm; the reaction holds the pivot in
        # place. The centrifugal term (cm - pivot) * omega² accounts for the
        # acceleration of the cm frame relative to the inertial frame.
        omega = Omega.z
        a_cf = (cm_w - self.pivot_w) * omega * omega
        g = pga.vec2dp(0.0, -9.81, 0.0)
        f_cm = self.m * (a_cf + g)
        F_cm = pga.wdg(cm_w, f_cm)
        F_piv = pga.wdg(self.pivot_w, -f_cm)
        F_b = pga.move2dp(F_cm + F_piv, pga.rrev(M))  # pull into body frame

        # Pivot constraint: Omega stays along Q_b. Project compute_omega_dot
        # onto Q_b via its scalar component (Q_b.z = 1).
        dOmega = pga.compute_omega_dot(self.I_inv, F_b, Omega, self.I)
        alpha = dOmega.z

        return [Omega, alpha * Q_b]   # [dB/dt, dOmega/dt]

    def calc_rkstep(self, dt):
        """One full RK4 step = 4 sub-steps. ga_py.rk4_step is the
        std::vector overload bound for Python; it takes the state by
        value and returns the updated state."""
        for rk_step in (1, 2, 3, 4):
            rhs = self.calc_rhs()
            self.u, self.uh = ga_py.rk4_step(self.u, self.uh, rhs, dt, rk_step)

    def print_sim(self, t):
        B = self.u[0]
        M = pga.rgpr(self.M0, pga.exp(0.5 * B))
        cm_w = pga.move2dp(pga.vec2dp(0.0, 0.0, 1.0), M)
        # B_b = phi * Q_b, so phi = B.z. World-frame B via M0 — the rotation
        # angle is frame-invariant, only the moment terms change.
        B_w = pga.move2dp(B, self.M0)
        print(f"t = {t:>-7.3f}: "
              f"cm_w = {cm_w:>-6.3f}, "
              f"phi_b = {B.z:>-6.3f}, "
              f"phi_w = {B_w.z:>-6.3f}, "
              f"B_w = {B_w:>-6.3f}, "
              f"B = {B:>-6.3f}")

    # Validation helpers ---------------------------------------------------
    def get_cm_world(self):
        M = pga.rgpr(self.M0, pga.exp(0.5 * self.u[0]))
        return pga.move2dp(pga.vec2dp(0.0, 0.0, 1.0), M)

    def get_omega(self):
        return self.u[1].z   # Omega = omega * Q_b, Q_b.z = 1

    def get_I_zz_pivot(self):
        return self.I[2, 2]


# ---------------------------------------------------------------------
# Run a full simulation: 2x2 plate, mass 1 kg, pinned at TR corner held
# fixed at the world origin. cm starts at (-1, -1) so the TR pivot
# coincides with world (0, 0). 60 RK4 steps over 2 seconds.
# ---------------------------------------------------------------------
print("pga2dp: combined motion (rigid 1-body system)\n")

m, w, h = 1.0, 2.0, 2.0
cm_pos = pga.vec2dp(-1.0, -1.0, 1.0)   # cm position (unitized point)
cm_phi = 0.0                            # initial rotation angle
cm_omega = 0.0                          # initial angular velocity

sim = SimOdePlatePga2dp(m, w, h, cm_pos, cm_phi, cm_omega)
sim.set_initial_values()
sim.print_sim(0.0)

# Static validation (matches the C++ test):
hw, hh = w / 2.0, h / 2.0
expected_I_pivot = 8.0 / 3.0
assert abs(sim.get_I_zz_pivot() - expected_I_pivot) < 1e-12, \
    f"I[2,2] = {sim.get_I_zz_pivot()} (expected 8/3 = {expected_I_pivot})"

# Initial mechanical energy (KE = 0, PE = m*g*cm.y).
g = 9.81
E_0 = m * g * cm_pos.y               # = -9.81 J
pivot_dist_sq = hw * hw + hh * hh    # = 2.0


# Time integration: 60 steps over [0, 2 s]. At every step verify
# (a) energy conservation and (b) the pivot constraint |cm|² = 2.
n_steps = 60
T_max = 2.0
dt = T_max / n_steps

for n in range(1, n_steps + 1):
    t = n * dt
    sim.calc_rkstep(dt)
    sim.print_sim(t)

    cm_w = sim.get_cm_world()
    omega = sim.get_omega()

    T_kin = 0.5 * sim.get_I_zz_pivot() * omega * omega
    V_pot = m * g * cm_w.y
    energy_drift = abs((T_kin + V_pot) - E_0)
    r_sq = cm_w.x * cm_w.x + cm_w.y * cm_w.y
    pivot_drift = abs(r_sq - pivot_dist_sq)

    # RK4 drift tolerance (matches the C++ test's 1e-3 epsilon).
    assert energy_drift < 1e-2 * abs(E_0), \
        f"energy not conserved at t={t}: T+V = {T_kin + V_pot}, E_0 = {E_0}"
    assert pivot_drift < 1e-2 * pivot_dist_sq, \
        f"pivot constraint broken at t={t}: |cm|² = {r_sq}, expected {pivot_dist_sq}"

print()
print(f"final energy drift     = {energy_drift:.3e} J  (limit 1e-2 * |E_0|)")
print(f"final pivot constraint = {pivot_drift:.3e} m² (limit 1e-2 * 2.0)")
