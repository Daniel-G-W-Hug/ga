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
