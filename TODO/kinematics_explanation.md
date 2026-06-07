# Kinematics documentation (collection for later)

This file collects explanatory material produced while building the kinematics layer
(`static_system2dp`, `kinematic_system2dp`, twists, velocity/acceleration). It is a staging
area — the pieces here are to be folded into the proper documentation (the `ga_docu` LaTeX
and/or the library headers) once the kinematic layer is complete.

**Notation** follows the project convention:

- `⟑` — geometric product (`gpr` in code)
- `⟇` — regressive geometric product (`rgpr` in code)
- `rrev(M)` — regressive reverse of `M` (`\undertilde{M}` / `\utilde{M}`)

In PGA, motors, the motor sandwich, and `exp()` are all built on the **regressive** product
`⟇`. In particular the rigid-motion sandwich is

```text
move2dp(v, M) = gr1(rgpr(rgpr(M, v), rrev(M)))   i.e.   M ⟇ v ⟇ rrev(M)
```

For a unit motor, `rrev(M)` is the inverse under `⟇`: `M ⟇ rrev(M) = 1` (the rgpr identity).

The sandwich is algebra-dependent — the product and the reverse must match. **EGA and STA**
use the *geometric* product with the ordinary reverse for their rotor sandwich,
`R ⟑ X ⟑ rev(R)` (`rotate(X, R)` / `transform(X, R)`). **PGA** uses the *regressive* product
with the regressive reverse, `M ⟇ X ⟇ rrev(M)` (`move2dp` / `move3dp`) — this is the form
used throughout the rest of this document, since the kinematics live in pga2dp.

## The adjoint and twist transport

### Two kinds of "change of frame"

There are two different objects you re-express from one frame to another, and they
transform by *different* rules:

- A **point** `X` (a location). Re-expressed by the motor sandwich:
  `X_world = move2dp(X_body, M) = M ⟇ X_body ⟇ rrev(M)`.
- A **twist** `B` (an instantaneous motion — an angular + linear velocity). Re-expressed by
  the **adjoint**, written `Ad_M(B)`.

The geometric-algebra punchline: the adjoint is the *same* sandwich.

```text
Ad_M(B) = M ⟇ B ⟇ rrev(M) = move2dp(B, M)
```

So `move2dp` serves double duty — point action *and* adjoint — because both are just
*conjugation* `M (·) rrev(M)`. That is why `twist_world` calls the same `move2dp` you would
use on a point. Conceptually they are still different operations; only the formula coincides.

### What a twist actually is

Two layers:

- **Motors** `M` — finite rigid motions. They form a **group**: compose by the product
  `M1 ⟇ M2` (`rgpr`); the identity is "no motion".
- **Twists** `B` — *infinitesimal* motions (velocities), the Lie algebra (tangent space at
  the identity). The bridge is the exponential `M = exp(½ B)` (w.r.t. `⟇`).

A twist is "the velocity you would integrate to get a motor". It lives in a **vector
space**, which is why twists *add*, while motors *multiply*.

### Why conjugation re-expresses a twist (intuition)

If `M` maps body coordinates to world coordinates and `B_body` is a motion described in the
body frame, then the *same physical motion* in world coordinates is obtained by

1. `rrev(M)` — step from world into the body frame,
2. apply the body motion `B_body`,
3. `M` — step back out to the world,

which is exactly `M ⟇ B_body ⟇ rrev(M) = Ad_M(B_body)`. The adjoint is "do the motion, but
in the other frame's language".

The clean identity behind it:

```text
M ⟇ exp(½ B) ⟇ rrev(M) = exp(½ · Ad_M(B))
```

— conjugating the *motor* equals building a motor from the *adjoint-transported twist*.
One-line proof: conjugation distributes through the power series, because
`M ⟇ Bⁿ ⟇ rrev(M) = (M ⟇ B ⟇ rrev(M))ⁿ` (the inner `rrev(M) ⟇ M = 1` telescopes).

### Why the velocities just add (`twist_world`)

In a chain, each joint `k` contributes a relative twist `ξ_k` expressed in *its own* frame.
To combine them, first put them all in a *common* (world) frame via the adjoint, then —
because twists live in a vector space — they superpose linearly:

```text
V_world(i) = Σ over the path   Ad_{M_k}(ξ_k)
           = Σ   move2dp(rel_twist[k], get_pos_trafo(k, 0))
```

Why add and not compose? Because they are *infinitesimal*. Finite motions do not commute
and must be multiplied; to first order, simultaneous infinitesimal motions just sum. This is
the GA version of the robotics "spatial Jacobian".

### The 2D picture (concrete)

For a planar motion `M` = rotate by θ, translate by `t`, the adjoint does two things to a
twist `(v, ω)`:

- the **angular part ω is unchanged** — a rotation rate is the same in every frame (hence
  `twist_world.z` is just the sum of the per-joint `ω`);
- the **linear part gets a lever-arm term**: `v ↦ R(θ)·v + ω·(t rotated 90°)`. That `ω × t`
  piece is the extra velocity a point picks up because the rotation acts at an offset.

So the adjoint *is* the rotating-frame transport theorem, encoded as one conjugation.

### The acceleration bracket

Acceleration is the time-derivative of velocity, so differentiate the adjoint. Because
`M_i(t)` is itself moving with spatial twist `V_i`:

```text
d/dt [ Ad_{M_i}(ξ_i) ] = Ad_{M_i}(ξ̇_i) + [ V_i , Ad_{M_i}(ξ_i) ]
```

The first term is the transported relative acceleration; the second — the **Lie bracket**
`[·,·]` (the `bracket(...)` in `world_VA`) — is the **Coriolis / centrifugal coupling**. The
bracket measures how much two twists fail to commute, which is what produces those cross
terms. That is the entire content of the Newton–Euler line
`A_i = A_parent + Ad(ξ̇_i) + [V_i, Ad(ξ_i)]`.

### Why a twist is a distinct semantic role (justifies `twist2dp`)

A point and a twist both ride the *same* `move2dp` sandwich, so in code they look
interchangeable — yet they transform under **different representations** of the motor group
(the *standard* action vs. the *adjoint* action) and live in different spaces (a point is a
location; a twist is a Lie-algebra element). They are genuinely different *kinds* of thing
that happen to share a storage type (`vec2dp`) and a transformation formula.

That is exactly why the weak `twist2dp` alias is appropriate: it labels the **semantic role**
("this `vec2dp` is a Lie-algebra twist, transformed via the adjoint") without pretending the
math needs a separate type — because the GA sandwich already handles both uniformly. In
pga3dp the twist is a genuine `BiVec3dp` (grade 2); the `twist2dp` alias keeps the 2D API
structurally parallel.

## Generators add, motors multiply (and do not commute)

This is the single most error-prone point in the whole framework, so it is worth stating
sharply. There are two levels, and additivity is valid on exactly one of them.

**Velocity level (the Lie algebra) — things add.** A twist splits additively into its
rotational and translational parts, `Ω = Ω_rot + Ω_tra` (the "Euclidean split"); the
velocity field is *linear* in the twist, `rcmt(Ω_rot + Ω_tra, X) = rcmt(Ω_rot, X) +
rcmt(Ω_tra, X)`; and a single generator accumulates component-wise, `B ← B + Ω·dt`. All of
this is correct because twists live in a vector space.

**Finite level (the Lie group) — things multiply, and order matters.** A finite motion is
applied by its motor `M = exp(½ B)`, and two *separate* motions compose by the regressive
geometric product (motor multiplication), which is **non-commutative**:

```text
M2 ⟇ M1  ≠  M1 ⟇ M2          (in general)
exp(½(A+B))  ≠  exp(½A) ⟇ exp(½B)     unless  rcmt(A, B) = 0
```

Adding generators equals multiplying motors *only* when the generators commute. For a pure
translation and a pure rotation that happens in exactly one situation: the **3D screw**
(translation **along** the rotation axis). In 2D it never happens for a nonzero rotation —
Chasles' theorem says every planar displacement is a *single* rotation about some center (or
a pure translation), so a rotation and a translation always collapse to one rotation about a
*shifted* center.

**The subtle trap.** Writing `B = B_rot + B_tra` and then `exp(½ B)` is *not* wrong — but it
produces the **single** screw (3D) / rotation-about-a-shifted-pivot (2D) generated by that
*combined twist*. It is **not** the composition "translate, then rotate". If you actually
want "translate by `t`, then rotate about `Q`", you must build and multiply two motors,
`M = M_rot ⟇ M_tra`, in the intended order — adding the generators gives a different motion.

This is exactly why `step_pos_trafo` composes the one-step parent→child pose as
`rgpr(M_rot, M_tra)` rather than `exp(−½(t + r))`: the single-exp form is a screw about the
wrong center and fails to map the child origin to `(0,0)`. *Generators are for velocities;
motors are for poses.*

(See `ga_docu/3_ga_modelling_motion.tex`, the subsection "Independent generators, coupled
finite displacements", eqn. `motors_do_not_commute`, for the same statement in the formal
documentation.)

## What GA actually buys for rigid-body dynamics (and what it doesn't)

This section is deliberately even-handed. The GA literature has a long habit of overselling
— sweeping claims of elegance and universality that don't survive contact with a real
implementation. The goal here is the opposite: state plainly where the plane-based GA (PGA)
formulation of rigid-body dynamics is genuinely advantageous, and where a classical method
is simply the better tool. The honest case is the more convincing one, and being clear
about the limits is the part the field has historically done badly.

**The baseline.** The fair comparison is not "GA vs. nothing" but GA vs. the standard
tooling: 4×4 `SE(3)` homogeneous matrices, or Featherstone's *spatial vector algebra*
(6-component twists/wrenches with a hand-defined spatial cross product `×*`). Both are
mature and fast. We are reorganising a solved problem, not solving an open one — so the
question is whether the reorganisation pays for itself.

### Where GA is genuinely advantageous

1. **Pose is a motor, not a 4×4 matrix.** No orthonormality constraint to maintain, no
   matrix renormalisation drift; one `exp`/`log` covers rotation + translation together as a
   screw. The pose evolves on the manifold exactly via `M ⟇ exp(½ Ω dt)`.
2. **Twists and wrenches are bivectors, and the "spatial cross product" is the algebra's
   commutator.** This is the strongest structural point. Featherstone's 6-vectors are
   Plücker coordinates of screws *by another name*, and his `×*` is defined by hand with
   sign conventions. In PGA a twist literally *is* a screw, the momentum/wrench is its dual
   bivector, and `×*` is `rcmt` (the regressive commutator) — a native product, not a
   bespoke operator. The geometry Featherstone encodes by convention, PGA encodes
   structurally.
3. **The adjoint equals the point sandwich.** Transporting a twist between frames is the
   *same* operation `move2dp(Ω, M) = M ⟇ Ω ⟇ rrev(M)` as moving a point — classical needs a
   separate 6×6 `Ad_M`. One operation, double duty (see "The adjoint and twist transport").
4. **Force-about-a-point is the wedge.** `wdg(Q, f)` is the applied wrench; momentum is the
   inertia map to a bivector. Moments are not a separate pseudovector with an `r × F` cross
   product.
5. **The Euler equation is one dimension-agnostic line:**
   `Ω̇ = I⁻¹(W − rcmt(Ω, I(Ω)))` (`compute_omega_dot`), identical in 2D and 3D up to which
   grades the twist and wrench occupy. The 2D→3D lift is the *same formula*, not a rewrite.

### The showcase: Coriolis/centrifugal coupling = a single commutator

The genuinely messy part of a classical double-pendulum Lagrangian is the
Coriolis/centrifugal coupling — the Christoffel symbols built from derivatives of the mass
matrix `M(q)`. In the geometric formulation this collapses to **one term: the Lie bracket
of twists** `[V_i, Ad(ξ_i)]`, already built in `world_VA`'s acceleration recursion (see "The
acceleration bracket"). GA's message is literally *"the Coriolis force is the geometric fact
that twists do not commute"*, and the commutator is a native product. That is a real
conceptual win — statable in one sentence and backed by running code.

These objects are also *drawable*, which 6-vectors are not:

- each body's **instantaneous centre of rotation** is a point read straight off the world
  twist, `pivot = (B.x/B.z, B.y/B.z)`; for the double pendulum it is a moving point tracing
  a curve;
- the **joint screw is literally the joint point** (`Ω_rot = ω·Q`): joint 2 is the moving
  corner of plate 1, and plate 2's relative twist is `ω₂` times that point;
- the **momentum / angular-momentum bivector** is a geometric object, not a column of
  numbers.

So the dynamical state *is* geometry — the visualisation is not a bolt-on.

### Where GA does not help (honest caveats)

- **Inertia is still a linear operator** (3×3 in 2D, 6×6 in 3D). GA represents it as a
  bivector→bivector map, but for an anisotropic body it does not collapse to a single
  multivector. No free lunch.
- **The joint-space mass matrix `M(q)` is still a small numeric solve.** Reducing to
  generalised coordinates is a *robotics* move, slightly against the GA grain: the geometry
  computes the entries, but the final `q̈ = M⁻¹(τ − h)` is ordinary linear algebra.
- **Performance is not the selling point.** Mature spatial-vector / matrix codes are highly
  optimised. PGA here wins on clarity, unification and visualisability — not speed.

### Are we just redoing the classical method?

Partly, and it is worth saying so. The *algorithm* (recursive Newton–Euler) is classical and
not ours to claim. What is GA is the **representation and the concepts** —
screws/motors/bivectors/commutator in place of matrices/6-vectors/hand-defined operators —
plus three things the standard tooling does not give for free: (a) the same primitives scale
unchanged from one body to a coupled tree; (b) the state is *visualisable geometry*; (c)
2D→3D is the same algebra, not a rewrite. That is a legitimate showcase: not "GA solves the
unsolved", but "GA expresses rigid-body dynamics more honestly, more unifiedly and more
visually than the standard tooling."

### The trade-off we actually made: reduced vs. maximal coordinates

There is a purer-GA formulation we deliberately did *not* choose: **maximal coordinates**,
where every body stays a full motor + twist and joints become *constraint wrenches*
(bivectors) — the "everything is a multivector, nothing reduces to a joint angle"
formulation. It is more GA-native, but its constraint stabilisation injects/removes energy,
which fights an energy-conservation test directly. We chose **reduced (joint) coordinates**:
the revolute joint *is* the constraint (relative twist forced to `ω·Q_joint`), so there are
no explicit constraint forces, no drift, and energy conservation is clean — at the cost of
the small `M(q)` solve above. That is roughly 90% of the GA showcase value while actually
conserving energy. Documenting *why* we passed on the purer option is part of the honest
account: the purer formulation is more elegant on paper and worse on the metric we care
about.

## Further reading (Lie groups vs. Lie algebras)

Orientation: the robotics literature treats these as *matrix* Lie groups (SE(2)/SE(3)) — the
adjoint is a matrix, the bracket a matrix commutator. GA/PGA represents the *same* groups as
motors, where the adjoint is the conjugation sandwich `M ⟇ (·) ⟇ rrev(M)` and the bracket is
`rcmt`. Read a robotics source for the *structure* and a PGA source for the *translation*.
Wikipedia is fine for isolated definitions (adjoint representation, exponential map, Lie
bracket) but too scattered to build the integrated picture quickly.

- **Solà, Deray, Atchuthan — "A micro Lie theory for state estimation in robotics"**
  (arXiv:1812.01537, ~17 pp). Compact, modern, built around SE(2)/SE(3): the group ↔ algebra
  correspondence, `exp`/`log`, the adjoint, the bracket, Jacobians. Best single quick read.
- **Lynch & Park — "Modern Robotics", Ch. 3** (free PDF + free online course). Twists, the
  adjoint, screws, exponential coordinates — the clearest intuition for velocity propagation;
  the Newton–Euler line `A_i = A_p + Ad(ξ̇) + [V, Ad(ξ)]` comes straight from here.
- **Murray, Li, Sastry — "A Mathematical Introduction to Robotic Manipulation", Ch. 2–3**
  (free PDF). The canonical, thorough treatment of twists, `Ad`, the Lie bracket, and the
  product of exponentials.
- **De Keninck & Dorst — "A Guided Tour to the Plane-based Geometric Algebra PGA"**
  (bivectors.net). Motors as `exp` of bivectors, bivectors as the Lie algebra — the GA/PGA
  bridge that connects the robotics structure to our sandwich / regressive products.
- **Ethan Eade — "Lie Groups for 2D and 3D Transformations"** (ethaneade.com). A terse
  `exp`/`log`/adjoint formula cheat-sheet for SO/SE.

Suggested path: Solà (concepts) → Lynch & Park Ch. 3 (twist/adjoint intuition) →
De Keninck & Dorst (see it as motors / bivectors).
