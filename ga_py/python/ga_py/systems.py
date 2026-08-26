"""ga_py.systems — pure-Python reconstruction of the PGA3DP multibody system classes.

The stateful C++ classes `static_system3dp` / `kinematic_system3dp` /
`dynamic_system3dp` (ga/ga_pga3dp_ops_mechanics.hpp) are NOT bound to Python by design.
This module rebuilds them from the bound GA primitives in `ga_py.pga`, so the same
frame-tree kinematics that the C++ classes provide is reachable from Python.

The reconstruction is faithful because the primitives carry the semantics: the motor
algebra (`motor_from_pose3dp` / `pose3dp_from_motor` / `rgpr` / `rrev` / `move3dp`), the
twist velocity field (`rcmt`) and `exp` do the geometric work. The "system" classes are
thin orchestration — a frame tree and a Newton-Euler recursion — so they port 1:1 to any
language that has the primitives. There is NO Python-side GA math here; every geometric
step dispatches into a bound C++ function.

C++ ↔ Python map (this module):

    static_system3dp::add_frame          -> StaticSystem3dp.add_frame
    static_system3dp::get_pos_trafo      -> StaticSystem3dp.get_pos_trafo   (LCA walk)
    static_system3dp::step_pos_trafo     -> StaticSystem3dp.step_pos_trafo
    kinematic_system3dp::to_twist        -> KinematicSystem3dp.to_twist
    kinematic_system3dp::twist_world     -> KinematicSystem3dp.twist_world
    kinematic_system3dp::velocity_field  -> KinematicSystem3dp.velocity_field
    kinematic_system3dp::accel_field     -> KinematicSystem3dp.accel_field
    kinematic_system3dp::point_velocity  -> KinematicSystem3dp.point_velocity
    kinematic_system3dp::point_accel...  -> KinematicSystem3dp.point_acceleration

`dynamic_system3dp` (the inertia / force-element tier, incl. the `extra_wrenches()`
subclass seam) is reconstructed in a later addition to this module.

Every method that takes a frame reference accepts either the integer frame index or the
frame name (string), unifying the C++ index/name overloads.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import pga

# Sentinel: default parent is the previously added frame (mirrors
# static_system3dp::prev_frame). A distinct object so a real index 0 is never mistaken
# for "use the default".
PREV_FRAME = object()

_ZERO = pga.vec3dp(0.0, 0.0, 0.0, 0.0)  # a zero DIRECTION (w = 0)


def _identity_pose() -> "pga.pose3dp":
    """The identity pose: origin at O (w = 1), no rotation."""
    return pga.pose3dp(pga.O_3dp, _ZERO)


def _identity_motor() -> "pga.mvec3dp_e":
    """The rgpr neutral element (the unit pseudoscalar I_3dp, as an even multivector).

    Built as the motor of the identity pose so no named-constant binding is required;
    mirrors static_system3dp using `I_3dp_mv_e` as the neutral element of `rgpr`.
    """
    return pga.motor_from_pose3dp(_identity_pose())


class StaticSystem3dp:
    """A tree of rigid coordinate frames, each posed relative to its parent.

    Reconstruction of `static_system3dp`: frames are added with a pose (relative to a
    parent), and `get_pos_trafo(from, to)` returns the motor that carries a point from
    one frame to another by walking up to the lowest common ancestor and back down —
    the identical motor-composition algorithm the C++ class uses.
    """

    def __init__(self) -> None:
        self._pose: list[pga.pose3dp] = []      # pose per frame (relative to parent)
        self._name: list[str] = []              # display name per frame
        self._parent: list[int] = []            # parent index per frame (root: itself)
        self._name_to_idx: dict[str, int] = {}

    # -- construction -----------------------------------------------------------------

    def add_frame(self, name: str, pose: "pga.pose3dp | None" = None,
                  parent_idx=PREV_FRAME) -> int:
        """Add a frame `name` with `pose` (relative to its parent) and return its index.

        `parent_idx` defaults to the previously added frame (a plain sequence of calls
        builds a linear chain); pass an explicit earlier index to branch. The first
        frame is the root (its parent is itself). `pose` defaults to the identity.
        """
        new_idx = len(self._pose)
        if pose is None:
            pose = _identity_pose()
        if pose.origin.w != 1.0:
            raise ValueError(
                "StaticSystem3dp: unitized projective origin (w == 1) expected, "
                f"got w == {pose.origin.w}")
        if parent_idx is PREV_FRAME:
            parent_idx = 0 if new_idx == 0 else new_idx - 1
        elif parent_idx >= new_idx:
            raise ValueError(
                "StaticSystem3dp: parent_idx must refer to an already added frame in "
                f"[0,{new_idx}), but got parent_idx == {parent_idx}")
        if name in self._name_to_idx:
            raise ValueError(
                f"StaticSystem3dp: frame name must be unique, but '{name}' exists")
        self._name_to_idx[name] = new_idx
        self._pose.append(pose)
        self._name.append(name)
        self._parent.append(parent_idx)
        return new_idx

    # -- lookup -----------------------------------------------------------------------

    def index_of(self, name: str) -> int:
        try:
            return self._name_to_idx[name]
        except KeyError:
            raise ValueError(f"StaticSystem3dp: no frame named '{name}'") from None

    def _resolve(self, ref) -> int:
        """Accept an int index or a frame name and return the index."""
        return ref if isinstance(ref, int) else self.index_of(ref)

    def size(self) -> int:
        return len(self._pose)

    def empty(self) -> bool:
        return not self._pose

    def parent(self, idx: int) -> int:
        return self._parent[idx]

    def name(self, idx: int) -> str:
        return self._name[idx]

    def pose(self, ref) -> "pga.pose3dp":
        return self._pose[self._resolve(ref)]

    def set_pose(self, ref, pose: "pga.pose3dp") -> None:
        self._pose[self._resolve(ref)] = pose

    # -- transforms -------------------------------------------------------------------

    def step_pos_trafo(self, child_idx: int) -> "pga.mvec3dp_e":
        """Single parent(child) -> child point transform, built from the child's pose.

        `motor_from_pose3dp` gives the body->parent motor; its regressive reverse is the
        parent->child transform (mirrors static_system3dp::step_pos_trafo).
        """
        return pga.rrev(pga.motor_from_pose3dp(self._pose[child_idx]))

    def _ancestor_chain(self, idx: int) -> list[int]:
        """[idx, parent(idx), ..., root]."""
        chain = [idx]
        while self._parent[idx] != idx:  # the root is its own parent
            idx = self._parent[idx]
            chain.append(idx)
        return chain

    def get_pos_trafo(self, from_ref, to_ref) -> "pga.mvec3dp_e":
        """Motor carrying a point from frame `from_ref` to frame `to_ref`.

        Walks up from `from` to the lowest common ancestor (LCA), then down to `to` —
        the identical algorithm to static_system3dp::get_pos_trafo. Use as
        `p_to = move3dp(p_from, M)`.
        """
        from_idx = self._resolve(from_ref)
        to_idx = self._resolve(to_ref)
        if from_idx >= self.size() or to_idx >= self.size():
            raise ValueError("StaticSystem3dp: frame index out of range")

        ident = _identity_motor()
        if from_idx == to_idx:
            return ident  # identity (the pseudoscalar, rgpr's neutral element)

        to_chain = self._ancestor_chain(to_idx)  # [to, ..., root]

        # M_up: from -> LCA. Each child->parent step is rrev(step_pos_trafo(child)) =
        # the body->parent motor; a further-up step multiplies on the LEFT.
        m_up = ident
        node = from_idx
        while node not in to_chain:
            m_up = pga.rgpr(pga.rrev(self.step_pos_trafo(node)), m_up)
            node = self._parent[node]
        lca = node

        # M_down: LCA -> to. Each parent->child step is step_pos_trafo(child); the
        # deepest child (to) ends up on the far LEFT.
        m_down = ident
        j = 0
        while to_chain[j] != lca:
            m_down = pga.rgpr(m_down, self.step_pos_trafo(to_chain[j]))
            j += 1

        return pga.rgpr(m_down, m_up)  # up-segment first, then down-segment


class KinematicSystem3dp(StaticSystem3dp):
    """A `StaticSystem3dp` with a momentary relative velocity/acceleration twist per
    frame; reconstruction of `kinematic_system3dp`.

    World-frame twists are propagated root -> frame by the Newton-Euler adjoint +
    commutator recursion (`world_VA`), and point velocities/accelerations follow from
    the twist velocity field `rcmt`.
    """

    def __init__(self) -> None:
        super().__init__()
        self._rel_vtwist: list[pga.bivec3dp] = []
        self._rel_atwist: list[pga.bivec3dp] = []

    @staticmethod
    def to_twist(v: "pga.vec3dp", omega: "pga.vec3dp") -> "pga.bivec3dp":
        """Pack (linear velocity, angular velocity) into a bivec3dp twist.

        Weight slots (e41,e42,e43) carry the angular velocity; bulk slots (e23,e31,e12)
        carry the linear velocity — the exact encoding of kinematic_system3dp::to_twist.
        """
        return pga.bivec3dp(omega.x, omega.y, omega.z, v.x, v.y, v.z)

    @staticmethod
    def _zero_twist() -> "pga.bivec3dp":
        return pga.bivec3dp(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # -- construction (adds the twist bookkeeping in sync with the frame list) ---------

    def add_frame(self, name: str, pose: "pga.pose3dp | None" = None,
                  state: "pga.kin_state3dp | None" = None,
                  parent_idx=PREV_FRAME) -> int:
        """Add a frame with `pose` and an optional momentary kinematic `state`
        (a `kin_state3dp`, relative to the parent). A frame without a state is at rest.
        """
        idx = super().add_frame(name, pose, parent_idx)
        if state is None:
            self._rel_vtwist.append(self._zero_twist())
            self._rel_atwist.append(self._zero_twist())
        else:
            self._rel_vtwist.append(self.to_twist(state.vel, state.omega))
            self._rel_atwist.append(self.to_twist(state.acc, state.alpha))
        return idx

    def set_state(self, ref, state: "pga.kin_state3dp") -> None:
        idx = self._resolve(ref)
        self._rel_vtwist[idx] = self.to_twist(state.vel, state.omega)
        self._rel_atwist[idx] = self.to_twist(state.acc, state.alpha)

    def set_twist(self, ref, twist: "pga.bivec3dp") -> None:
        self._rel_vtwist[self._resolve(ref)] = twist

    def relative_twist(self, ref) -> "pga.bivec3dp":
        return self._rel_vtwist[self._resolve(ref)]

    # -- world twists & point kinematics -----------------------------------------------

    def twist_world(self, ref) -> "pga.bivec3dp":
        """World-frame velocity twist of a frame: the sum along the root -> frame path
        of each relative twist transported to world by that frame's world motor (the
        adjoint `move3dp`). Spatial twists add: V_i = V_parent + Ad(xi_i).
        """
        idx = self._resolve(ref)
        v = self._zero_twist()
        n = idx
        while self._parent[n] != n:
            v = v + pga.move3dp(self._rel_vtwist[n], self.get_pos_trafo(n, 0))
            n = self._parent[n]
        return v

    @staticmethod
    def velocity_field(twist: "pga.bivec3dp", x: "pga.vec3dp") -> "pga.vec3dp":
        """PGA rate of change of a point under a twist: Xdot = rcmt(V, X). Argument
        order matters — rcmt(V, X) == -rcmt(X, V).
        """
        return pga.rcmt(twist, x)

    @staticmethod
    def accel_field(vel_twist: "pga.bivec3dp", acc_twist: "pga.bivec3dp",
                    x: "pga.vec3dp") -> "pga.vec3dp":
        """Rigid-body acceleration field at X: rcmt(A, X) + rcmt(V, rcmt(V, X))
        (frame/Euler + centripetal).
        """
        return pga.rcmt(acc_twist, x) + pga.rcmt(vel_twist, pga.rcmt(vel_twist, x))

    def point_velocity(self, x_world: "pga.vec3dp", ref) -> "pga.vec3dp":
        """Velocity of a world-space point X rigidly attached to a frame."""
        return self.velocity_field(self.twist_world(ref), x_world)

    def point_acceleration(self, x_world: "pga.vec3dp", ref) -> "pga.vec3dp":
        """Acceleration of a world-space point X rigidly attached to a frame."""
        v, a = self._world_va(self._resolve(ref))
        return self.accel_field(v, a, x_world)

    def accel_twist_world(self, ref) -> "pga.bivec3dp":
        """World-frame acceleration twist A = dV/dt of a frame."""
        _, a = self._world_va(self._resolve(ref))
        return a

    def _world_va(self, idx: int) -> "tuple[pga.bivec3dp, pga.bivec3dp]":
        """World velocity & acceleration twists of frame idx by the Newton-Euler
        recursion (mirrors kinematic_system3dp::world_VA):

            V_i = V_parent + Ad(xi_i)
            A_i = A_parent + Ad(xidot_i) + [V_i, Ad(xi_i)]

        The se(3) twist Lie bracket [.,.] is the regressive commutator rcmt.
        """
        path: list[int] = []
        n = idx
        while True:
            path.append(n)
            if self._parent[n] == n:
                break
            n = self._parent[n]
        path.reverse()  # root -> idx

        v = self._zero_twist()
        a = self._zero_twist()
        for n in path:
            if self._parent[n] == n:
                continue  # the root contributes nothing
            m = self.get_pos_trafo(n, 0)
            zeta = pga.move3dp(self._rel_vtwist[n], m)     # Ad(xi_n)
            zetadot = pga.move3dp(self._rel_atwist[n], m)  # Ad(xidot_n)
            v = v + zeta
            a = a + zetadot + pga.rcmt(v, zeta)            # [V, zeta]
        return v, a


# ---------------------------------------------------------------------------
# Dynamic tier — inertia, force elements and forward dynamics.
# Reconstruction of `dynamic_system3dp`.
# ---------------------------------------------------------------------------

# Joint kinds (mirror of the C++ `joint3dp` enum). Both 1-DOF kinds run through
# the same exponential M(q) = rest ⟇ rexp(½ q·screw) — only the screw differs.
FREE = "free"
REVOLUTE = "revolute"
PRISMATIC = "prismatic"


@dataclass
class Body3dp:
    """Rigid-body inertial properties of a frame (mirror of `body3dp`): the inertia map
    `I` (body twist -> body momentum), its cached inverse `I_inv`, and the total mass.
    """

    I: "pga.inertia3dp"
    I_inv: "pga.inertia3dp"
    mass: float = 0.0


@dataclass
class JointState3dp:
    """Per-frame joint state (mirror of `joint_state3dp`). `screw_b` is the body-frame
    twist generator (a rotation-axis line for a revolute, an ideal line for a prismatic);
    `rest` is the body->parent motor at q = 0. The optional spring/damper acts on the
    generalised coordinate: tau += -stiffness·(q - q_rest) - damping·q̇.
    """

    type: str = FREE
    screw_b: "pga.bivec3dp" = field(default=None)
    rest: "pga.mvec3dp_e" = field(default=None)
    phi: float = 0.0
    omega: float = 0.0
    stiffness: float = 0.0
    damping: float = 0.0
    q_rest: float = 0.0


@dataclass
class GroundedSpring3dp:
    """A grounded spatial spring + damper (mirror of `grounded_spring3dp`): a body-fixed
    point `anchor_b` tied to the inertial anchor `p0_world` by anisotropic world-axis
    stiffness `k` + isotropic damping `c`. One element yields both translational and (via
    its lever arm) tilt stiffness.
    """

    anchor_b: "pga.vec3dp"
    p0_world: "pga.vec3dp"
    k: "pga.vec3dp"
    c: float = 0.0


def _empty_inertia() -> "pga.inertia3dp":
    return pga.inertia3dp()


def make_cuboid_body(m: float, w: float, h: float, d: float) -> Body3dp:
    """Body for a uniform cuboid (extents w,h,d along e1,e2,e3), origin at the cm.
    m = 0 builds a massless carrier link (a joint holder without inertia); its inverse
    inertia stays zero, as in the C++ builder -- the inverse only serves the free-body
    integrator, which never runs for a massless body."""
    I = pga.get_cuboid_inertia(m, w, h, d)
    return Body3dp(I, pga.get_inertia_inverse(I) if m > 0.0 else _empty_inertia(), m)


class DynamicSystem3dp(KinematicSystem3dp):
    """Forces/inertia tier on top of `KinematicSystem3dp`; reconstruction of
    `dynamic_system3dp`.

    1-DOF revolute/prismatic joints form a coupled chain integrated in reduced (joint)
    coordinates by joint-space forward dynamics `M(q) q̈ = τ` (the shared `lu_solve`);
    free 6-DOF bodies integrate independently by the se(3) Euler equation. Force elements
    — joint spring/damper, time-varying applied wrench, grounded spatial spring, driven
    (moving-base) joints — all feed the generalised force additively. A subclass injects
    application-specific wrenches by overriding `extra_wrenches()`.
    """

    def __init__(self) -> None:
        super().__init__()
        self._body: list[Body3dp] = []
        self._joint: list[JointState3dp] = []
        self._grav = pga.vec3dp(0.0, -9.81, 0.0, 0.0)  # world direction (w = 0)
        self._wrench: dict[int, "callable"] = {}       # frame -> fn(t) -> bivec3dp
        self._driven: dict[int, tuple[float, float]] = {}  # frame -> (rate, q0)
        self._springs: dict[int, list[GroundedSpring3dp]] = {}
        self._time = 0.0

    # -- extension seam (mirror of the protected virtual extra_wrenches) ---------------

    def extra_wrenches(self) -> "list[tuple[int, pga.bivec3dp]]":
        """Injection point for a subclass: return (frame_idx, world-wrench) pairs folded
        onto each frame's supporting joints, evaluated every sub-step. The generic base
        returns none.
        """
        return []

    # -- construction ------------------------------------------------------------------

    def add_frame(self, name: str, pose: "pga.pose3dp | None" = None,
                  state: "pga.kin_state3dp | None" = None,
                  parent_idx=PREV_FRAME) -> int:
        idx = super().add_frame(name, pose, state, parent_idx)
        self._body.append(Body3dp(_empty_inertia(), _empty_inertia(), 0.0))
        self._joint.append(JointState3dp())
        return idx

    def add_body(self, name: str, body: Body3dp,
                 state: "pga.kin_state3dp | None" = None,
                 pose: "pga.pose3dp | None" = None, parent_idx=PREV_FRAME) -> int:
        """Add a free (6-DOF) rigid body with inertial properties + initial state."""
        idx = self.add_frame(name, pose, state, parent_idx)
        self._body[idx] = body
        return idx

    def add_revolute_body(self, name: str, body: Body3dp, pivot_b: "pga.vec3dp",
                          axis_b: "pga.vec3dp", phi0: float = 0.0, omega0: float = 0.0,
                          pose: "pga.pose3dp | None" = None, parent_idx=PREV_FRAME) -> int:
        """Add a 1-DOF hinge about the body-fixed axis line through `pivot_b` along unit
        `axis_b`; the joint screw is the line wdg(pivot_b, axis_b).
        """
        return self._add_screw_joint(name, body, REVOLUTE, pga.wdg(pivot_b, axis_b),
                                     phi0, omega0, pose, parent_idx)

    def add_prismatic_body(self, name: str, body: Body3dp, direction: "pga.vec3dp",
                           s0: float = 0.0, v0: float = 0.0,
                           pose: "pga.pose3dp | None" = None, parent_idx=PREV_FRAME) -> int:
        """Add a 1-DOF slider along body-fixed unit `direction`; the joint screw is the
        translation generator (an ideal line).
        """
        screw = pga.bivec3dp(0.0, 0.0, 0.0, direction.x, direction.y, direction.z)
        return self._add_screw_joint(name, body, PRISMATIC, screw, s0, v0, pose, parent_idx)

    def _add_screw_joint(self, name, body, jtype, screw_b, q0, qdot0, pose,
                         parent_idx) -> int:
        idx = self.add_frame(name, pose, None, parent_idx)
        self._body[idx] = body
        js = JointState3dp(type=jtype, screw_b=screw_b,
                           rest=pga.rrev(self.step_pos_trafo(idx)), phi=q0, omega=qdot0)
        self._joint[idx] = js
        self._apply_joint_state(idx)
        return idx

    # -- force elements & configuration ------------------------------------------------

    def set_gravity(self, g: "pga.vec3dp") -> None:
        self._grav = g

    def gravity(self) -> "pga.vec3dp":
        return self._grav

    def set_joint_spring_damper(self, ref, k: float, c: float, q0: float = 0.0) -> None:
        js = self._joint[self._resolve(ref)]
        js.stiffness, js.damping, js.q_rest = k, c, q0

    def set_applied_wrench(self, ref, fn) -> None:
        """Attach a time-varying world wrench fn(t) -> bivec3dp to a frame (None clears)."""
        idx = self._resolve(ref)
        if fn is None:
            self._wrench.pop(idx, None)
        else:
            self._wrench[idx] = fn

    def set_driven_rate(self, ref, rate: float, q0: float = 0.0) -> None:
        """Prescribe a 1-DOF joint at constant rate q(t) = q0 + rate·t (a moving base)."""
        idx = self._resolve(ref)
        self._driven[idx] = (rate, q0)
        self._apply_driven_joints()

    def clear_driven_joint(self, ref) -> None:
        self._driven.pop(self._resolve(ref), None)

    def is_driven(self, ref) -> bool:
        return self._resolve(ref) in self._driven

    def add_grounded_spring(self, ref, anchor_b: "pga.vec3dp", k: "pga.vec3dp",
                            c: float = 0.0, p0_world: "pga.vec3dp | None" = None) -> None:
        """Tie a body-fixed point `anchor_b` to an inertial anchor by anisotropic world
        stiffness `k` + isotropic damping `c`. If `p0_world` is omitted the current world
        position of the point is taken as the rest position.
        """
        idx = self._resolve(ref)
        if p0_world is None:
            p0_world = pga.unitize(pga.move3dp(anchor_b, self.get_pos_trafo(idx, 0)))
        self._springs.setdefault(idx, []).append(
            GroundedSpring3dp(anchor_b, p0_world, k, c))

    def clear_grounded_springs(self, ref) -> None:
        self._springs.pop(self._resolve(ref), None)

    def time(self) -> float:
        return self._time

    def set_time(self, t: float) -> None:
        self._time = t

    def body_props(self, ref) -> Body3dp:
        return self._body[self._resolve(ref)]

    def joint_phi(self, ref) -> float:
        return self._joint[self._resolve(ref)].phi

    def joint_omega(self, ref) -> float:
        return self._joint[self._resolve(ref)].omega

    # -- forward dynamics & integration ------------------------------------------------

    def joint_accel(self, ref) -> float:
        """Current acceleration of a joint from the coupled forward dynamics (no step)."""
        rj = self._dof_joints()
        qdd = self._forward_dynamics(rj)
        idx = self._resolve(ref)
        for k, j in enumerate(rj):
            if j == idx:
                return qdd[k]
        return 0.0

    def sync_accelerations(self) -> None:
        """Write the dynamic joint accelerations into the relative accel twists, so the
        world accel queries return the actual accelerations (not just the bias)."""
        rj = self._dof_joints()
        qdd = self._forward_dynamics(rj)
        for k, j in enumerate(rj):
            self._rel_atwist[j] = qdd[k] * self._joint[j].screw_b

    def step(self, dt: float) -> None:
        """Advance the system by dt: the 1-DOF joints as a coupled RK4 chain, free bodies
        independently. Threads the sub-step clock into any time-varying wrench / driven
        joint, then prescribes the driven joints at t + dt.
        """
        rj = self._dof_joints()
        if rj:
            self._coupled_step(rj, dt)
        for i in range(1, self.size()):
            if self._joint[i].type == FREE and self._body[i].mass > 0.0:
                self._step_free_body(i, dt)
        self._time += dt
        self._apply_driven_joints()

    # -- energy / momentum diagnostics -------------------------------------------------

    def kinetic_energy(self) -> float:
        ke = 0.0
        for i in range(1, self.size()):
            minv = pga.rrev(self.get_pos_trafo(i, 0))          # world -> body i
            vb = pga.move3dp(self.twist_world(i), minv)        # body twist
            ke += 0.5 * self._spatial_dot(vb, self._body[i].I(vb))
        return ke

    def potential_energy(self) -> float:
        pe = 0.0
        g = self._grav
        for i in range(1, self.size()):
            cm = pga.move3dp(pga.O_3dp, self.get_pos_trafo(i, 0))
            pe += -self._body[i].mass * (g.x * cm.x + g.y * cm.y + g.z * cm.z)
            dq = self._joint[i].phi - self._joint[i].q_rest
            pe += 0.5 * self._joint[i].stiffness * dq * dq
        for fi, sps in self._springs.items():
            m = self.get_pos_trafo(fi, 0)
            for sp in sps:
                p = pga.unitize(pga.move3dp(sp.anchor_b, m))
                dx, dy, dz = p.x - sp.p0_world.x, p.y - sp.p0_world.y, p.z - sp.p0_world.z
                pe += 0.5 * (sp.k.x * dx * dx + sp.k.y * dy * dy + sp.k.z * dz * dz)
        return pe

    def total_energy(self) -> float:
        return self.kinetic_energy() + self.potential_energy()

    def momentum_world(self, ref) -> "pga.bivec3dp":
        """World-frame momentum bivector of a body (conserved for a torque-free body)."""
        idx = self._resolve(ref)
        m = self.get_pos_trafo(idx, 0)                       # body -> world
        vb = pga.move3dp(self.twist_world(idx), pga.rrev(m))  # body twist
        return pga.move3dp(self._body[idx].I(vb), m)          # momentum back to world

    def mass_matrix(self) -> "list[float]":
        """Joint-space mass matrix M(q) (n*n row-major) at the current configuration."""
        rj = self._dof_joints()
        n = len(rj)
        S = [pga.move3dp(self._joint[rj[i]].screw_b, self.get_pos_trafo(rj[i], 0))
             for i in range(n)]
        minv = [pga.rrev(self.get_pos_trafo(rj[i], 0)) for i in range(n)]
        mmat = [0.0] * (n * n)
        for i in range(n):
            I = self._body[rj[i]].I
            for j in range(n):
                if not self._is_ancestor(rj[j], rj[i]):
                    continue
                xj = pga.move3dp(S[j], minv[i])
                for k in range(n):
                    if not self._is_ancestor(rj[k], rj[i]):
                        continue
                    xk = pga.move3dp(S[k], minv[i])
                    mmat[j * n + k] += self._spatial_dot(xj, I(xk))
        return mmat

    # -- internals ---------------------------------------------------------------------

    @staticmethod
    def _spatial_dot(xi: "pga.bivec3dp", mom: "pga.bivec3dp") -> float:
        """Spatial (reciprocal) pairing twist x wrench -> scalar: -rwdg(xi, mom)."""
        return -float(pga.rwdg(xi, mom))

    def _joint_motor(self, idx: int, q: float) -> "pga.mvec3dp_e":
        js = self._joint[idx]
        return pga.rgpr(js.rest, pga.rexp(0.5 * q * js.screw_b))

    def _apply_joint_state(self, idx: int) -> None:
        js = self._joint[idx]
        self.set_pose(idx, pga.pose3dp_from_motor(self._joint_motor(idx, js.phi)))
        self.set_twist(idx, js.omega * js.screw_b)  # rel twist = q̇·screw

    def _apply_driven_joints(self) -> None:
        for idx, (rate, q0) in self._driven.items():
            self._joint[idx].phi = q0 + rate * self._time
            self._joint[idx].omega = rate
            self._apply_joint_state(idx)

    def _dof_joints(self) -> "list[int]":
        return [i for i in range(1, self.size())
                if self._joint[i].type in (REVOLUTE, PRISMATIC) and i not in self._driven]

    def _is_ancestor(self, jf: int, bf: int) -> bool:
        n = bf
        while True:
            if n == jf:
                return True
            if self.parent(n) == n:
                return False
            n = self.parent(n)

    def _assemble_mass_bias(self, rj: "list[int]") -> "tuple[list[float], list[float]]":
        """Assemble M(q) and the generalised-force RHS for the chain `rj` (mirror of
        assemble_mass_bias). SIDE EFFECT: runs the bias pass, zeroing the chain's relative
        accel twists so the world accel queries return only the velocity-product bias.
        """
        n = len(rj)
        z = self._zero_twist()
        for j in rj:
            self._rel_atwist[j] = z  # velocity-product (bias) pass

        S = [pga.move3dp(self._joint[rj[j]].screw_b, self.get_pos_trafo(rj[j], 0))
             for j in range(n)]

        # inertia-bearing bodies: dof joints AND driven (moving-base) joints
        bl = list(rj) + list(self._driven.keys())
        nb = len(bl)
        cm = [None] * nb
        minv = [None] * nb
        fbias = [None] * nb
        mass = [0.0] * nb
        for i in range(nb):
            fb = bl[i]
            m = self.get_pos_trafo(fb, 0)
            cm[i] = pga.move3dp(pga.O_3dp, m)
            minv[i] = pga.rrev(m)
            mass[i] = self._body[fb].mass
            vb = pga.move3dp(self.twist_world(fb), minv[i])
            ab = pga.move3dp(self.accel_twist_world(fb), minv[i])  # rel_atwist = 0
            I = self._body[fb].I
            fbias[i] = I(ab) + pga.rcmt(vb, I(vb))  # I·A_bias + gyroscopic V ×* (I V)

        g = self._grav
        mmat = [0.0] * (n * n)
        rhs = [0.0] * n
        for j in range(n):
            for i in range(nb):
                if not self._is_ancestor(rj[j], bl[i]):
                    continue
                vj = self.velocity_field(S[j], cm[i])          # cm velocity, unit rate j
                xj = pga.move3dp(S[j], minv[i])                # joint-j screw in body i
                rhs[j] += mass[i] * (vj.x * g.x + vj.y * g.y + vj.z * g.z) \
                    - self._spatial_dot(xj, fbias[i])
                I = self._body[bl[i]].I
                for k in range(n):
                    if not self._is_ancestor(rj[k], bl[i]):
                        continue
                    xk = pga.move3dp(S[k], minv[i])
                    mmat[j * n + k] += self._spatial_dot(xj, I(xk))

        # joint spring/damper (additive, diagonal): τ_j += -k(q-q0) - c·q̇
        for j in range(n):
            js = self._joint[rj[j]]
            rhs[j] += -js.stiffness * (js.phi - js.q_rest) - js.damping * js.omega

        # applied wrenches (world, at the current clock)
        for fi, fn in self._wrench.items():
            if fn is None:
                continue
            w = fn(self._time)
            for j in range(n):
                if self._is_ancestor(rj[j], fi):
                    rhs[j] += self._spatial_dot(S[j], w)

        # grounded spatial springs/dampers (configuration-dependent)
        for fi, sps in self._springs.items():
            m = self.get_pos_trafo(fi, 0)
            vw = self.twist_world(fi)
            for sp in sps:
                p = pga.unitize(pga.move3dp(sp.anchor_b, m))
                v = self.velocity_field(vw, p)
                fx = -sp.k.x * (p.x - sp.p0_world.x) - sp.c * v.x
                fy = -sp.k.y * (p.y - sp.p0_world.y) - sp.c * v.y
                fz = -sp.k.z * (p.z - sp.p0_world.z) - sp.c * v.z
                w = pga.wdg(p, pga.vec3dp(fx, fy, fz, 0.0))
                for j in range(n):
                    if self._is_ancestor(rj[j], fi):
                        rhs[j] += self._spatial_dot(S[j], w)

        # subclass-provided wrenches (empty in the generic base)
        for fi, w in self.extra_wrenches():
            for j in range(n):
                if self._is_ancestor(rj[j], fi):
                    rhs[j] += self._spatial_dot(S[j], w)

        return mmat, rhs

    def _forward_dynamics(self, rj: "list[int]") -> "list[float]":
        mmat, rhs = self._assemble_mass_bias(rj)
        import ga_py
        return ga_py.lu_solve(mmat, rhs, len(rj))

    def _write_u(self, rj: "list[int]", u: "list[float]") -> None:
        n = len(rj)
        for k in range(n):
            self._joint[rj[k]].phi = u[k]
            self._joint[rj[k]].omega = u[n + k]
            self._apply_joint_state(rj[k])

    def _coupled_step(self, rj: "list[int]", dt: float) -> None:
        n = len(rj)
        dim = 2 * n
        u0 = [self._joint[rj[k]].phi for k in range(n)] + \
             [self._joint[rj[k]].omega for k in range(n)]
        t0 = self._time

        def f(t, u):
            self._time = t
            self._apply_driven_joints()
            self._write_u(rj, u)
            qdd = self._forward_dynamics(rj)
            return [u[n + k] for k in range(n)] + qdd

        k1 = f(t0, u0)
        k2 = f(t0 + 0.5 * dt, [u0[i] + 0.5 * dt * k1[i] for i in range(dim)])
        k3 = f(t0 + 0.5 * dt, [u0[i] + 0.5 * dt * k2[i] for i in range(dim)])
        k4 = f(t0 + dt, [u0[i] + dt * k3[i] for i in range(dim)])
        u = [u0[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(dim)]
        self._time = t0
        self._write_u(rj, u)

    def _step_free_body(self, idx: int, dt: float) -> None:
        """RK4-integrate one free rigid body over dt on the Lie-algebra pair (B, Omega):
        dB/dt = dexp⁻¹(B) Omega = Omega + ½[B, Omega] + 1/12 [B, [B, Omega]] (the bracket
        is rcmt; the series is truncated after the double bracket, enough for fourth
        order since B = O(dt) restarts at 0 each step), dOmega/dt = I⁻¹[W_body -
        rcmt(Omega, I(Omega))] (compute_omega_dot). Pose evolves on the motor manifold
        M(t) = M0 ⟇ rexp(½ B). Bracket terms included -- without them the scheme is only
        second order on the pose.

        Mirror scope: the C++ library integrates a massive free body as a 6-coordinate
        motor joint INSIDE the coupled solve (a floating base: a chain hung on it loads
        it). This mirror keeps the older separate route, which coincides with the
        library for an isolated free body -- the case its cross-check covers -- but
        does not reproduce a floating base carrying joints.
        """
        m0 = pga.rrev(self.step_pos_trafo(idx))  # current body -> parent motor
        bd = self._body[idx]
        I, I_inv, mass = bd.I, bd.I_inv, bd.mass
        g = self._grav

        def omega_dot(B, Om):
            m = pga.rgpr(m0, pga.rexp(0.5 * B))
            f_world = pga.vec3dp(mass * g.x, mass * g.y, mass * g.z, 0.0)
            w_w = pga.wdg(pga.move3dp(pga.O_3dp, m), f_world)  # gravity wrench (world)
            w_b = pga.move3dp(w_w, pga.rrev(m))                # into the body frame
            return pga.compute_omega_dot(I_inv, w_b, Om, I)

        def deriv(B, Om):
            c1 = pga.rcmt(B, Om)
            b_dot = Om + 0.5 * c1 + (1.0 / 12.0) * pga.rcmt(B, c1)  # dexp⁻¹(B) Omega
            return b_dot, omega_dot(B, Om)  # (dB/dt, dOmega/dt)

        B0 = self._zero_twist()
        Om0 = self.relative_twist(idx)
        kB1, kO1 = deriv(B0, Om0)
        kB2, kO2 = deriv(B0 + 0.5 * dt * kB1, Om0 + 0.5 * dt * kO1)
        kB3, kO3 = deriv(B0 + 0.5 * dt * kB2, Om0 + 0.5 * dt * kO2)
        kB4, kO4 = deriv(B0 + dt * kB3, Om0 + dt * kO3)
        B = B0 + dt / 6.0 * (kB1 + 2 * kB2 + 2 * kB3 + kB4)
        Om = Om0 + dt / 6.0 * (kO1 + 2 * kO2 + 2 * kO3 + kO4)

        self.set_pose(idx, pga.pose3dp_from_motor(pga.rgpr(m0, pga.rexp(0.5 * B))))
        self.set_twist(idx, Om)
