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
