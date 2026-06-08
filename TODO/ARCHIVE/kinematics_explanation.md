# Frame trees, kinematic chains & articulated dynamics in PGA (documentation collection)

This file collects the explanatory material for the **multi-frame system layer** built on top
of the single-body PGA mechanics already developed in `ga_docu`: the classes `static_system`,
`kinematic_system`, `dynamic_system` (in `ga/ga_pga{2,3}dp_ops_physics.hpp`). It is a
**staging area**: once reviewed and refined it folds into `ga_docu/5_ga_modelling_physics.tex`,
directly **after** the subsubsection "Moving coordinate systems"
(`subsubsection:moving_coordinate_systems`), as the extension from *one* moving body to a
*tree* of coupled frames.

**This builds on existing chapters — it does not re-derive them.** Chapters 3 and 5 already
establish, for a single rigid body, everything this layer composes: motors and `exp`
(`3_ga`, `subsection:motor_generators_PGA`), the generators-add-vs-motors-multiply distinction
and Chasles/screw (`3_ga`, "Independent generators, coupled finite displacements"), the
body↔world transport of twists by the sandwich (`5_ga`, `eq:body_to_world_pga`), the force
line / wrench `F = Q∧f` (`5_ga`, `subsubsection:modelling_force_and_torque`), the momentum
wrench `P = I[Ω]` (`subsubsection:modelling_momentum`), the inertia map
(`subsubsection:structure_inertia_maps`, `eq:inertia_map_components_{2D,3D}`), and the body-frame
Euler ODE `Ω̇_b = I_b⁻¹[F_b − rcmt(Ω_b, I_b[Ω_b])]` (`eq:ode_system_body_frame`). The new
material below **references** those results and concentrates on what they do *not* cover: an
arbitrary **frame tree**, velocity/acceleration **transport along a chain**, and **articulated
(jointed) coupled dynamics** in reduced coordinates.

**Scope and emphasis.** Not an API reference — no exhaustive argument lists. The goal is that a
reader understands the concepts and the basic mathematics of the system layer, with the reading
list (§8) pointing to primary sources. Code names are anchors back into the library, not the
subject.

**LaTeX transfer plan (agreed):**

- §1–6 become the formal capability extension in the normal text flow, with the new
  forward-dynamics / mass-matrix **equations spelled out**; overlapping foundations stay as
  `\ref`/`\cite` cross-references into chapters 3/5.
- §7 (the honest GA assessment) and §8 (further reading) form a separate **addendum**; the
  literature references are defined there and `\cite`d from the main text.

Presented for `pga2dp`; everything lifts unchanged to `pga3dp` (only grades change — §6).

**Notation** (project convention): `⟑` geometric product (`gpr`), `⟇` regressive geometric
product (`rgpr`), `rrev(M)` regressive reverse. In PGA the motor sandwich is
`move2dp(v, M) = M ⟇ v ⟇ rrev(M)` (`5_ga`, `eq:pga_sandwich_tranformation`); for a unit motor
`M ⟇ rrev(M) = 1`.

**Primal vs. dual PGA.** This library is the **primal** `G(2,0,1)`/`G(3,0,1)` (Lengyel: points/
lines primary, motions via the regressive product `⟇`), *not* the **dual** plane-based
`G*(2,0,1)`/`G*(3,0,1)` of Gunn and Dorst & De Keninck (planes primary, motions via the
geometric product). The two are dual-equivalent via `a ⟇ b ≡ lcmpl(rcmpl(a) ∧ rcmpl(b))`. This
is already stated in `5_ga` (`subsection:modelling_physics`, "inspired by … but uses Lengyel's
direct representation") and `1_ga_intro.tex`. Consequence for the reading list: concepts from
plane-based sources carry over, but their *formulae* are the dual of ours — §8 tags each source
`[primal]` / `[dual]` / `[metric-free]`.

---

## 1. Architecture: three layers on one frame tree

The implementation is a strict three-layer stack, each layer *inheriting* the one below and
adding exactly one new kind of state. The geometry — a tree of coordinate frames — is shared;
the layers differ only in what they attach to each frame.

| Layer                 | Adds per frame                      | Answers                              |
| --------------------- | ----------------------------------- | ------------------------------------ |
| `static_system2dp`    | a **pose** (relative to parent)     | *where* is each frame?               |
| `kinematic_system2dp` | a **velocity + acceleration twist** | *how is it moving?* (no forces)      |
| `dynamic_system2dp`   | **inertia + joints + gravity**      | *how do forces move it?* (Newton)    |

```text
static_system2dp  ──▶  kinematic_system2dp  ──▶  dynamic_system2dp
   pose tree              + twists                  + inertia / joints / forces
```

Each capability is *additive* and testable on its own: the transformation algebra at the static
layer, velocity/acceleration transport at the kinematic layer, energy conservation at the
dynamic layer. A user who only needs coordinate transforms never pays for the dynamics. The
same three classes exist in `pga3dp` with identical structure and algorithms — only the
representations change (§6).

---

## 2. Pose, frame, and the frame tree (`static_system`)

The existing chapters move a *single* body by a motor `M(t)`. The static system generalises
this to a **tree of frames**, each posed relative to its parent, with rigid transforms between
*any* two frames.

### 2.1 Pose — a relative rigid placement

A **pose** is an origin plus an orientation, expressing one frame relative to another (in 2D:
`{vec2dp origin; value_t phi;}`). It is the human-readable form of a **motor** `M`; the two
interconvert via `motor_from_pose2dp` (build `M = translate(origin) ⟇ rotate(phi)`) and
`pose2dp_from_motor` (read the origin off `move2dp(O, M)` and the angle off the rotated `e1`).
The motor is built exactly as in `3_ga` (`eq:motor_pga_body_frame`); the pose decode is the 2D
special case of a constrained motor `log`, exact on round-trip.

> **3D note.** `pose3dp` replaces the scalar `phi` with an *axis·angle* vector. The 2D scalar
> angle is the only dimension-specific piece of the pose.

### 2.2 Frame — a named pose, always relative to its parent

`static_frame2dp` is a *named* pose. The invariant: **a frame's pose is always stored relative
to its parent**, never in world coordinates. The world (root) frame is its own parent.

### 2.3 The frame tree and its capabilities

`static_system2dp` holds a **tree** of frames (each records its parent). Adding a frame without
naming a parent links it to the previously added one — a plain sequence builds a linear chain
`root → f1 → f2 → …`; naming an earlier frame branches the tree. This buys, for free:

- **Branching** — a body frame with two arms: both arm chains parent off the body, the body off
  the (possibly moving) world. Any sub-assembly is a subtree.
- **Name addressing** — frames carry unique names; operations have by-name overloads
  (`get_pos_trafo("hand", "world")`).
- **Reposing** — `set_pose(idx, …)` moves a frame relative to its parent; the whole subtree
  follows automatically, since all descendant poses are relative.

> **Worked usage.** The frame tree is exercised end-to-end in two places. In `ga_view`, the
> *frame-transformation* scene (`ga_view/src/active_frame_trafo.{hpp,cpp}`) drives a
> `static_system2dp` directly — building a parent/child tree and reading transforms between
> frames — while the *merry-go-round* scene
> (`ga_view/src/active_merry_go_round.{hpp,cpp}`) shows the same tree carrying motion one layer
> up (a `kinematic_system2dp` rotating platform, advanced via `step(dt)`). For non-graphical,
> assertion-checked examples see the system-class test cases in
> `ga_test/src/ga_appl2dp_physics_test.hpp` (and the 3D analogue
> `ga_appl3dp_physics_test.hpp`).

### 2.4 Transformation logic — one sandwich for points *and* vectors

`get_pos_trafo(from, to)` returns a single motor `M`; **any** geometric object is transformed by
the *same* `move2dp(X, M)`. This is worth stressing because classical mechanics hands you several
different rules: a **point** transforms affinely (rotate *and* translate); a **free vector**
linearly (rotate only); and **basis vectors** transform *covariantly* while a vector's
*coordinates* transform *contravariantly* (the inverse-transpose).

In PGA there is one operation, and the object's **homogeneous weight** selects the behaviour
automatically: a point (weight `1`) is rotated *and* translated; a direction (weight `0`, a
point at infinity) is rotated *only*, because translating a point at infinity does nothing. No
contravariant/covariant split survives — the motor sandwich is a single self-consistent
conjugation transporting every grade (points, lines, directions, and — §3 — twists and momenta)
correctly. The classical "which transformation law?" decision disappears into the algebra.

### 2.5 How `get_pos_trafo` walks the tree (LCA)

Because poses are relative, the transform between two arbitrary frames is composed along the tree
path. `get_pos_trafo(from, to)`:

1. walks **up** from `from` to the **lowest common ancestor** (LCA), composing each child→parent
   step (the regressive reverse of that frame's pose motor);
2. walks **down** from the LCA to `to`, composing each parent→child step;
3. returns `rgpr(M_down, M_up)`.

Special cases fall out of the same code: `from == to` → identity (the unit pseudoscalar, neutral
for `⟇`); a pure ancestor/descendant relation makes one segment empty (forward/backward
transform); a linear chain is the degenerate tree where the LCA is the shallower frame.
Composition is motor multiplication `⟇`, **dimension-neutral** — the tree machinery lifts to 3D
unchanged.

---

## 3. Transport of motion along a chain (`kinematic_system`)

Chapter 5 transports the twist of a *single* body between its body and world frames
(`eq:body_to_world_pga`, `eq:rcmt_pga_world`). The kinematic system generalises this to a
**chain of frames**: each frame carries a relative velocity and acceleration twist (w.r.t. its
parent), and the layer answers velocity/acceleration queries for any point of any frame.

### 3.1 Twists and the adjoint (recap + the new ingredient)

A **twist** is an instantaneous rigid motion (angular + linear velocity) — the Lie-algebra
element with `M = exp(½B)`; twists **add** while motors **multiply** (`3_ga`,
`subsection:motor_generators_PGA`). The library tags it with a weak alias
(`twist2dp = vec2dp`, encoding `(−v_y, v_x, ω)`; in 3D a genuine `BiVec3dp`) to mark the
*semantic role* — "this is a twist, transported by the adjoint" — without a new type.

The one fact needed below: re-expressing a twist in another frame is the **adjoint**, and in GA
the adjoint *is the same sandwich* used on a point,
`Ad_M(B) = M ⟇ B ⟇ rrev(M) = move2dp(B, M)` (this is `eq:body_to_world_pga` read as a twist
transport). So `move2dp` serves double duty — point action and twist adjoint — because both are
conjugation `M(·)rrev(M)`.

### 3.2 Velocities add along the path — `twist_world` (new)

This is the genuinely new step beyond the single-body chapter. In a chain, each frame `k`
contributes a relative twist `ξ_k` in *its own* frame. Transport each to the world frame by the
adjoint, then — because twists live in a vector space — **superpose**:

```text
V_world(i) = Σ over the path root→i   Ad_{M_k}(ξ_k)
           = Σ   move2dp(rel_twist[k], get_pos_trafo(k, 0))
```

They add (not compose) because they are infinitesimal: finite motions must be multiplied, but
to first order simultaneous infinitesimal motions sum. This is the GA form of the robotics
**spatial Jacobian** — and it is exactly what a single moving body does not exhibit.

### 3.3 Velocity and acceleration fields (recap)

Given a frame's world twist, the velocity of any point `X` of that frame is the regressive
commutator `velocity_field(V, X) = rcmt(V, X)`, and
`accel_field(V, A, X) = rcmt(A, X) + rcmt(V, rcmt(V, X))` — i.e. the frame (Euler) and
centripetal terms of the "Moving coordinate systems" subsubsection (`5_ga`, table there;
`eq:rcmt_pga_world`), now applied per chain frame. `point_velocity` / `point_acceleration` build
the world twists for a frame and evaluate the field at `X`.

### 3.4 Acceleration transport along the chain — the Coriolis bracket (new)

Differentiating the adjoint along a moving chain gives the recursive Newton–Euler relations:

```text
V_i = V_parent + Ad(ξ_i)
A_i = A_parent + Ad(ξ̇_i) + [ V_i , Ad(ξ_i) ]
```

The second acceleration term is the **Lie bracket** of twists — the **Coriolis / centrifugal
coupling** — and in this primal PGA it is a native product, `[A, B] = rcmt(A, B)` (the se(2)/
se(3) twist bracket; the implementation routes through `rcmt` directly, component-identical to
the hand-derived planar formula). The single-frame transport theorem is in chapter 5; the
**recursion over a chain** is what the kinematic system adds.

---

## 4. Generators add, motors multiply — pointer

The velocity-vs-finite distinction (twists add in the Lie algebra; motors multiply
non-commutatively, `exp(½(A+B)) ≠ exp(½A)⟇exp(½B)` unless `rcmt(A,B)=0`; Chasles; the screw) is
developed fully in `3_ga`, "Independent generators, coupled finite displacements". It is the
reason `step_pos_trafo` composes the parent→child pose as `rgpr(M_rot, M_tra)` rather than a
single `exp` of summed generators — *generators are for velocities, motors are for poses*. Only
this implementation consequence belongs in the new section; the theory is cross-referenced.

---

## 5. Articulated coupled dynamics (`dynamic_system`)

Chapter 5 solves a *single* free rigid body (`eq:ode_system_body_frame`) and even a *single*
revolute pivot (`eq:pivot_invariance`, `eq:constraint_projection_2d`). The dynamic system adds
the genuinely new dynamics: **several bodies coupled by joints**, integrated together in
**reduced (joint) coordinates**.

> **Worked usage.** The canonical example is the **double pendulum** — two revolute-jointed
> plates swinging chaotically under gravity while conserving energy. It is visualised in
> `ga_view` (`ga_view/src/active_double_pendulum.{hpp,cpp}`, a `dynamic_system2dp` advanced via
> `step(dt)`) and pinned by assertion-checked tests: the planar case in
> `ga_test/src/ga_appl2dp_physics_test.hpp` (TEST_CASE *"pga2dp: dynamic_system2dp - double
> pendulum energy conservation (M3)"*) and the genuinely-3D, non-parallel-axes case in
> `ga_test/src/ga_appl3dp_physics_test.hpp` (TEST_SUITE *"PGA3DP: dynamic_system3dp (M3)"*,
> TEST_CASE *"pga3dp: spatial double pendulum - non-parallel axes (M3)"*). Both check energy
> conservation to RK4 tolerance and the mass-matrix-vs-kinetic-energy identity (§5.3).

### 5.1 Inertia, energy and the spatial pairing (recap + helper)

Each body carries the inertia map `I` (twist → momentum wrench `P = I(Ω)`) of
`subsubsection:structure_inertia_maps`, with mass and rotational inertia in one bivector→bivector
operator. The implementation pairs a twist with a wrench through the dimension-agnostic
**spatial / reciprocal product** `spatial_dot(ξ, P) = −rwdg(ξ, P)`, so that
`spatial_dot(ξ, I(ξ)) = 2·KE` carries mass and angular terms uniformly (no `½m|v|²+½Iω²` split).
This one pairing drives kinetic energy, the mass matrix, and forward dynamics, and lifts to 3D
by swapping the map for `Inertia3dp` and `rwdg` for the `BiVec3dp` reciprocal.

### 5.2 Joints — the PGA unification of rotation and translation

A **joint** is the reduced coordinate linking a body to its parent: `free` (unconstrained
3-DOF), `revolute` (1-DOF hinge), or `prismatic` (1-DOF slider). The unification: **both 1-DOF
joints run identical code.** The body→parent motor is `M(q) = rest ⟇ exp(½ q · screw)`, the
relative twist is `q̇ · screw`, the Jacobian column is `velocity_field(screw, ·)` — only the
**screw generator** differs: a **finite point** `Q_b` (revolute → `exp` is a rotation about
`Q_b`) versus an **ideal point / direction** (prismatic → `exp` is a translation). Chapter 5
treats the single revolute pivot (the invariance `eq:pivot_invariance` is the `rest ⟇ exp`
constraint for one joint); the *generic* screw machinery, the prismatic case, and Steiner
inertia emerging automatically from the spatial Jacobian are new here.

### 5.3 Forward dynamics — reduced coordinates over a chain (new)

The new dynamics is the coupling. In reduced coordinates the joint *is* the constraint, so there
are no explicit constraint forces, no drift, and clean energy — at the cost of a small
joint-space solve. `forward_dynamics` assembles, by virtual work over the bodies, the
joint-space equation of motion `M(q) q̈ = RHS(q, q̇)`:

```text
M[j][k] = Σ_i  spatial_dot( S_j^{body i}, I_i( S_k^{body i} ) )     (joint-space mass matrix)
RHS[j]  = Σ_i  m_i ( v_cm,i(S_j)·g  −  v_cm,i(S_j)·b_cm,i )         (gravity − Coriolis bias)
```

`S_j` is the world joint screw; `S_j^{body i} = move2dp(S_j, rrev(M_i))` transports it into body
`i`'s frame (where its inertia map lives); `v_cm,i(S_j) = velocity_field(S_j, cm_i)` is the
spatial-Jacobian column; `b_cm,i` is the velocity-product (Coriolis/centripetal) cm-acceleration
at `q̈ = 0`. The sums run over bodies having joint `j` (and `k`) as an ancestor; the system is
solved by the shared LU solver. The identity `½ q̇ᵀ M(q) q̇ == kinetic_energy()` is a built-in
consistency check. Each body's contribution still obeys the per-body Euler equation
`eq:ode_system_body_frame`; the mass matrix is its projection onto the joint coordinates.

> **3D note.** The 3D RHS uses the full spatial bias **wrench**
> `F_bias = I(A_bias) + rcmt(V, I(V))` (the `rcmt` term is the gyroscopic torque), which
> vanishes in the planar case to recover the 2D `RHS`. The per-body Euler line
> `Ω̇ = I⁻¹(W − rcmt(Ω, I(Ω)))` is identical in 2D and 3D up to the grade of `Ω`.

### 5.4 Time integration and diagnostics

`step(dt)` advances by RK4: a coupled joint chain is integrated together in reduced coordinates
(`coupled_step`), free bodies independently (`step_free_body`). The energy diagnostics
(`kinetic_energy`, `potential_energy`, `total_energy`) are the correctness metric the design is
tuned for: **energy conserved to RK4 tolerance** — validated for the free body, compound
pendulum, prismatic slider, and chaotic double pendulum, in both 2D and 3D.

---

## 6. The 2D → 3D lift in one table

The `pga3dp` stack is the same classes, methods and algorithms; only the representations change:

| concept           | `pga2dp`                  | `pga3dp`                         |
| ----------------- | ------------------------- | -------------------------------- |
| motor (pose)      | `mvec2dp_u` (odd)         | `mvec3dp_e` (even)               |
| pose angle        | scalar `phi`              | axis·angle `rot` (a vector)      |
| twist (velocity)  | `vec2dp` (gr1 shadow)     | `BiVec3dp` (gr2)                 |
| point             | `vec2dp`, z = 1           | `Vec3dp`, w = 1                  |
| joint screw       | a finite/ideal point      | a line / screw axis (`BiVec3dp`) |
| inertia           | 3×3 map (`Inertia2dp`)    | 6×6 map (`Inertia3dp`)           |
| twist bracket     | `rcmt(vec2dp, vec2dp)`    | `rcmt(BiVec3dp, BiVec3dp)`       |
| spatial pairing   | `−rwdg(vec2dp, bivec2dp)` | `−rwdg(BiVec3dp, BiVec3dp)`      |

The sandwich, the adjoint, the tree walk, the joint exponential, the inertia-map pairing and the
Euler equation are **identical formulae** — the lift is a change of representation, not a
rewrite. That is itself the headline result of doing 2D fully first.

---

> **ADDENDUM (not part of the main text flow).** In the LaTeX document the following two sections
> form a separate addendum; the literature references are defined here and `\cite`d normally
> from the main text where relevant.

## 7. What GA actually buys for rigid-body dynamics (and what it doesn't)

Deliberately even-handed. The GA literature has a habit of overselling; the goal here is to state
plainly where the projective-GA (PGA) formulation is genuinely advantageous and where a classical
method is the better tool. (This library is the **primal** PGA — see the notation section; the
advantages below are of the geometric formulation as such and hold for either representation.)

**The baseline** is not "GA vs. nothing" but GA vs. 4×4 `SE(3)` matrices or Featherstone's
spatial vector algebra (6-vectors with a hand-defined `×*`) — both mature and fast. We are
reorganising a solved problem, so the question is whether the reorganisation pays.

### 7.1 Where GA is genuinely advantageous

1. **Pose is a motor, not a 4×4 matrix.** No orthonormality constraint, no renormalisation
   drift; one `exp`/`log` covers rotation + translation as a screw; the pose evolves on the
   manifold exactly via `M ⟇ exp(½ Ω dt)`.
2. **Twists, momenta and wrenches are bivectors, and the "spatial cross product" is the
   commutator.** Featherstone's 6-vectors are Plücker screw coordinates by another name and his
   `×*` is hand-defined; here a twist *is* a screw, the momentum wrench is `I(Ω)`, the bias
   wrench `I(A_bias)+rcmt(V,I(V))` is a bivector, and `×*` is `rcmt` — native, not bespoke.
3. **The adjoint equals the point sandwich.** Transporting a twist is the same `move2dp(Ω, M)`
   as moving a point — classical needs a separate 6×6 `Ad_M`. One operation, double duty (§3.1).
4. **A force is its line of action — a wrench.** `F = wdg(Q, f) = Q∧f` carries force (weight) and
   torque-about-origin (bulk) in one bivector — the literature's **wrench** (Dorst's "forque").
   Moments are not a separate pseudovector with `r × F`. Developed in `5_ga`,
   `subsubsection:modelling_force_and_torque`.
5. **The Euler equation is one dimension-agnostic line** (`eq:ode_system_body_frame`),
   `Ω̇ = I⁻¹(W − rcmt(Ω, I(Ω)))`, identical in 2D and 3D up to grade — the lift is the same
   formula, not a rewrite.

### 7.2 The showcase: Coriolis/centrifugal coupling = one commutator

The messy part of a classical double-pendulum Lagrangian is the Coriolis/centrifugal coupling —
Christoffel symbols from derivatives of `M(q)`. Here it collapses to **one term: the Lie bracket
of twists** `[V_i, Ad(ξ_i)] = rcmt(…)`, already present in the chain acceleration recursion
(§3.4). GA's message, *"the Coriolis force is the geometric fact that twists do not commute"*, is
one sentence backed by running code.

These objects are also *drawable* (unlike 6-vectors): each body's **instantaneous centre of
rotation** is a point read off the world twist, `(B.x/B.z, B.y/B.z)` (`3_ga` already gives this
for the single body); the **joint screw is literally the joint point**; the **momentum wrench**
is a geometric bivector. The dynamical state *is* geometry — the visualisation is not a bolt-on.

### 7.3 Where GA does not help (honest caveats)

- **Inertia is still a linear operator** (3×3 / 6×6); for an anisotropic body it does not
  collapse to a single multivector.
- **The mass matrix `M(q)` is still a small numeric solve.** Reducing to generalised coordinates
  is a robotics move: the geometry computes the entries, but `q̈ = M⁻¹(τ − h)` is ordinary linear
  algebra.
- **Performance is not the selling point.** Mature spatial-vector / matrix codes are highly
  optimised; PGA wins on clarity, unification and visualisability, not speed.

### 7.4 The trade-off we made: reduced vs. maximal coordinates

A purer-GA formulation we deliberately did *not* choose is **maximal coordinates** — every body a
full motor + twist, joints as *constraint wrenches* — more GA-native, but its constraint
stabilisation injects/removes energy, fighting an energy-conservation test. We chose **reduced
(joint) coordinates**: the revolute joint *is* the constraint, so no explicit constraint forces,
no drift, clean energy — at the cost of the small `M(q)` solve. ~90% of the showcase value while
actually conserving energy. Documenting *why* we passed on the purer option is part of the honest
account.

### 7.5 Are we just redoing the classical method?

Partly. The *algorithm* (recursive Newton–Euler) is classical. What is GA is the **representation
and concepts** — screws/motors/bivectors/commutator vs. matrices/6-vectors/hand-defined operators
— plus three things the standard tooling does not give for free: (a) the same primitives scale
unchanged from one body to a coupled tree; (b) the state is *visualisable geometry*; (c) 2D → 3D
is the same algebra. A legitimate showcase: not "GA solves the unsolved", but "GA expresses
rigid-body dynamics more honestly, unifiedly and visually than the standard tooling."

---

## 8. Further reading — understanding the concepts behind the library

Each idea (motors, the sandwich, the adjoint, twists as a Lie algebra, the screw unification, the
inertia map, the Coriolis commutator) is a standard result with a primary source. The grouping
goes from the algebra, through plane-based PGA, to rigid-body mechanics, and finally the
robotics / Lie-theory and classical viewpoints the GA formulation translates.

**Standing caveat — primal vs. dual.** As set out in the notation section, `5_ga`
(`subsection:modelling_physics`) and `1_ga_intro.tex`, this library is the **primal**
`G(2,0,1)`/`G(3,0,1)` (Lengyel), while Gunn and Dorst & De Keninck are **dual** plane-based.
Read those for the **concepts** (representation-independent: motors, bivector twists, the
adjoint, the screw, the commutator), but translate their **formulae** through
`a ⟇ b ≡ lcmpl(rcmpl(a) ∧ rcmpl(b))` before matching code. Tags: `[primal]`, `[dual]`,
`[metric-free]`.

**1 — Geometric algebra foundations:**

- `[metric-free]` **Macdonald — "Linear and Geometric Algebra"** (`Macdonald:2010`). Gentle,
  rigorous first contact with the geometric product, blades, rotors.
- `[metric-free]` **Hestenes — "New Foundations for Classical Mechanics"**
  (`David_Hestenes_NFCM:1986`). The origin of mechanics in GA; his ch. 7 on rigid-body mechanics
  already points toward PGA (as `5_ga` notes). The conceptual ancestor of §5.
- `[metric-free]` **Browne — "Grassmann Algebra, Vol. 1"** (`Browne_Grassmann-Algebra_Vol1:2012`).
  Exterior-algebra background for the wedge / regressive products, complements, and the force/
  momentum lines (`5_ga` cites it for exactly this).

**2 — Plane-based PGA and the primal/dual choice:**

- `[primal]` **Lengyel — "Projective Geometric Algebra Illuminated"**
  (`Lengyel_pga-illuminated:2024`; projectivegeometricalgebra.org). **The convention this library
  follows.** Argument for the primal choice: `Lengyel_poor-foundations_GA:2024` /
  `Lengyel_transwedge_product:2025`.
- `[dual]` **Gunn — "Projective geometric algebra: a new framework for doing Euclidean geometry"**
  (`Gunn_PGA:2020`). The foundational plane-based paper — read for structure, dualize formulae.
- `[dual]` **Dorst & De Keninck — "A Guided Tour to the Plane-based Geometric Algebra PGA"**
  (`Plane-based_PGA_Dorst-DeKennik:2022`). The most accessible PGA tour; excellent intuition,
  dual formulae.

**3 — Rigid-body mechanics in PGA** (closest to §3–§5; both plane-based/dual):

- `[dual]` **Dorst & De Keninck — "May the Forque Be with You — Dynamics in PGA"**
  (`Dynamics_in_PGA_plane-based_Dorst-DeKennik:2023`). The direct analogue of the dynamics layer:
  momentum bivector, the "forque", the motor equation of motion — in the dual model
  (geometric where ours is regressive). `5_ga` explicitly adapts its inertia-map results "to our
  approach using regressive products". The single most relevant comparison for §5/§7.
- `[dual]` **Charles Gunn — "Geometry, Kinematics, and Rigid Body Mechanics in Cayley–Klein
  Geometries"** (PhD thesis, TU Berlin, 2011; *not yet in the bib*). The most complete treatment
  of rigid-body kinematics/dynamics with motors and bivector twists in (dual) projective GA.

**4 — Lie-theory / robotics translation** `[metric-free]` (same groups SE(2)/SE(3) as matrices;
map matrix↔motor, `Ad`↔sandwich, bracket↔`rcmt`). *Not yet in the bib — add on transfer.*

- **Solà, Deray, Atchuthan — "A micro Lie theory for state estimation in robotics"**
  (arXiv:1812.01537). Compact: group ↔ algebra, `exp`/`log`, adjoint, bracket, Jacobians. Best
  quick read for the Lie picture.
- **Lynch & Park — "Modern Robotics", Ch. 3** (free PDF + course). Twists, adjoint, screws; the
  recursion `A_i = A_p + Ad(ξ̇) + [V, Ad(ξ)]` (§3.4) is straight from here.
- **Murray, Li, Sastry — "A Mathematical Introduction to Robotic Manipulation", Ch. 2–3** (free
  PDF). Canonical thorough treatment of twists, `Ad`, the Lie bracket, product of exponentials.
- **Selig — "Geometric Fundamentals of Robotics"** (2005). Screw theory and SE(3) as a Lie group.

**5 — The classical baseline §7 compares against** `[metric-free]` *(not yet in the bib)*:

- **Featherstone — "Rigid Body Dynamics Algorithms"** (2008). Spatial vector algebra: 6-component
  twists/wrenches and `×*` — the Plücker cousin of our bivector twists and `rcmt`, and the
  reference for the forward-dynamics algorithm (§5.3).

**Suggested path.** Macdonald (GA) → Lengyel *PGA Illuminated* (our primal PGA) → Dorst & De
Keninck tour + *Forque* dynamics (PGA mechanics, mind the dual) → Solà and Lynch & Park Ch. 3
(the Lie/robotics translation of twists and the adjoint).
