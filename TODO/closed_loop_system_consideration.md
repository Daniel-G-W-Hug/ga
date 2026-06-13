# Closed-loop systems — findings and implementation plan

Status: design note / proposal. **Phase 0 (reuse-seam prep) landed**; Phases 1-5 not yet
implemented. Captures (a) why the current `static_/kinematic_/dynamic_system{2,3}dp` tier
is **open-chain only**, and (b) a plan to add **closed-loop** (parallel mechanism) support
as a *separate, additive* layer that reuses the open-loop code without complicating it.

Origin: question "can the three-stage tree approach model closed-chain systems like a
Stewart-Gough platform?" — answer below.

**Progress (2026-06-13).**

- Phase 0 done: `forward_dynamics` split into `assemble_mass_bias` (the reuse seam) + a
  thin solve, in both `ga/ga_pga2dp_ops_physics.hpp` and `ga/ga_pga3dp_ops_physics.hpp`.
  Behaviour-preserving, validated by the existing physics suites (details in Phase 0).
- Prerequisite fixed: a clang-22 doctest `__LINE__` symbol collision that broke the
  `ga_appl3dp_test` build (so the 3D physics tests could not even run) — resolved in
  `ga_test/CMakeLists.txt` (`DOCTEST_COUNTER=__COUNTER__` for all test targets).

---

## 1. Findings — the current approach is structurally open-chain only

The three-stage classes live in [ga/ga_pga2dp_ops_physics.hpp](../ga/ga_pga2dp_ops_physics.hpp)
and [ga/ga_pga3dp_ops_physics.hpp](../ga/ga_pga3dp_ops_physics.hpp)
(`static_system → kinematic_system → dynamic_system`, each deriving from the previous).
Three independent properties make them a **rooted tree with reduced (joint) coordinates**,
which excludes kinematic loops by construction:

1. **One parent per frame.** `std::vector<size_t> parent_of` stores exactly one parent
   index per frame (the root's parent is itself). `add_frame(rf, parent_idx)` takes a
   single parent. A Stewart-Gough platform is reached through *six* legs at once — its
   frame would need six parents. The data model cannot express that.

2. **Transforms assume a unique path.** `get_pos_trafo(from, to)` walks *up to the lowest
   common ancestor, then down*. A unique LCA / unique path exists only in a tree; a loop
   gives multiple paths and no unique LCA.

3. **Reduced coordinates + no constraint machinery.** Joints are only
   `{free, revolute, prismatic}`; each is a 1-DOF (or free) relation *to its parent*.
   `step()` integrates the tree in joint-angle coordinates via `forward_dynamics`
   (assembles `M(q) q̈ = RHS` and `lu_solve`s it). A search for
   `loop | closure | constraint | lagrange | multiplier | baumgarte` finds nothing; the
   only "loop" in the code is the root self-loop that terminates the upward walk.

This is a **deliberate** design choice, already documented in
[ga_docu/5_ga_modelling_physics.tex](../ga_docu/5_ga_modelling_physics.tex) ("Reduced vs.
maximal coordinates"): maximal coordinates with explicit constraint wrenches were *not*
chosen because their constraint stabilization injects/removes energy and fights the
energy-conservation test the design is tuned for. The reduced (joint) coordinates make the
joint *itself* the constraint, so there are no constraint forces, no drift, clean energy —
at the cost of being tree-only.

### What this means for a Stewart-Gough platform

Cannot be modelled with the current tier as-is. The same holds for any closed loop:
four-bar linkage, delta robot, parallelogram linkage, etc.

### What already works / partial capability

- **Each leg as an independent open chain** — fine.
- **Inverse geometry with the platform pose given as input** — the static layer can
  compute leg geometry (anchor points as `vec`, leg lines/lengths via `wdg`/join and
  norms). PGA is elegant here. But the platform pose must be *supplied*, not solved.
- The hard part it cannot do: **forward kinematics / dynamics** (legs → platform pose),
  which *is* the loop-closure problem.

---

## 2. Requirement / boundary condition for the extension

Keep the existing open-loop tier exactly as it is — *simple and separately usable*. The
closed-loop capability must be a **separate, additional** part of the library, built *on
top of* the open-loop primitives by **composition** (uses what is there, draws in what is
needed), not by folding constraint logic into the existing classes. A user who only needs
open chains never includes the new header and pays no complexity. This mirrors how the
open-loop tier stands on its own today.

---

## 3. Recommended approach — reduced tree + loop-closure constraints (composition)

Build the open-chain part with the existing `dynamic_system{2,3}dp` as a **spanning tree**
(unchanged). Then close loops by registering **loop-closure constraints** between existing
tree frames, and solve the resulting *constrained* dynamics in a new layer. The user keeps
building a tree as today, then calls `add_loop_constraint(...)` to close loops — no manual
"cutting" needed; topology = the existing tree + extra constraint edges.

This is **coordinate partitioning / reduced-coordinates-with-constraints**, not maximal
coordinates.

### 3.1 Why this over maximal coordinates

| | reduced tree + loop constraints (recommended) | maximal coordinates |
| --- | --- | --- |
| Reuse of existing code | **high** — reuses `M(q)`, the spatial Jacobian, `forward_dynamics`, RK4 | low — every body free, the tree reduction unused |
| System size | small (n joint DOF + m constraints) | large (6·#bodies + all joint constraints) |
| Energy behaviour | good — few constraints to stabilise | needs heavier stabilisation (the energy concern the docs raise) |
| GA-nativeness | constraints are PGA meet/coincidence | also GA-native (joints as constraint wrenches) |
| Topology generality | any loop expressible as tree + closure edges | any topology uniformly |

Recommend the reduced/hybrid path because it directly satisfies the "reuse, stay separate,
preserve energy" boundary condition. Maximal coordinates is noted as a future alternative
(a cleaner fit if fully general non-tree topologies are ever needed) but is **out of scope
for the first implementation**.

### 3.2 Constrained equations of motion (the math)

Tree dynamics from the existing assembly: `M(q) q̈ = τ(q, q̇)` (with `τ` = gravity +
Coriolis/centripetal bias, both produced by `assemble_mass_bias`, the Phase 0 split). Loop
closure adds `m`
scalar constraints `g(q) = 0` with Jacobian `G = ∂g/∂q`. Differentiating twice gives the
acceleration-level KKT system solved each step:

```text
| M   Gᵀ | | q̈ |   |  τ                |
|        | |    | = |                   |
| G   0  | | λ  |   | -Ġ q̇  (+ stab.)   |
```

`λ` are the Lagrange multipliers (the constraint/leg forces). Solve via the shared LU
solver (`detail/ga_solver.hpp`) on the bordered system, or a Schur complement
`G M⁻¹ Gᵀ λ = G M⁻¹ τ + Ġ q̇` reusing `M⁻¹` from the existing path.

### 3.3 Constraint primitives (PGA-native)

Start with the one primitive that covers four-bar and Stewart-Gough leg attachments and
grow as needed:

- **point coincidence** (a planar pin / a spatial spherical joint): two anchor points,
  `P_a` on one branch and `P_b` on another, must coincide. Residual = the unitized
  difference `P_a − P_b` (2 scalar eqns in 2D, 3 in 3D), or equivalently the vanishing of
  their join. This is the loop-closure condition for Stewart-Gough (each leg top pinned to
  the platform).
- later: **point-on-line**, **distance / rigid-link length**, **parallel/perpendicular** —
  all expressible as vanishing `wdg`/`rwdg`/`dot` expressions.

A `loop_constraint{2,3}dp` descriptor: `{ frame_a, anchor_a, frame_b, anchor_b, type }`.

### 3.4 Constraint Jacobian = reuse of the spatial Jacobian columns

Key reuse insight: the constraint Jacobian rows are **relative partial velocities** of the
two anchor points, and the partial velocity of a point w.r.t. a joint rate is already
computed during `assemble_mass_bias` as `velocity_field(S_j, P) = rcmt(S_j, P)` (a spatial
Jacobian column). So for each constraint and each joint `j` on the path to its anchors:

```text
G[c][j] = (∂P_a/∂q_j − ∂P_b/∂q_j)·basis = velocity_field(S_j, P_a) − velocity_field(S_j, P_b)
```

No new geometric machinery — the same columns that build `M` and `τ` build `G`. `Ġ q̇` is
the constraint-level velocity-product term, obtainable from the same bias pass that
already yields `b_cm` (relative accel twists zeroed).

### 3.5 Constraint stabilisation — projection (recommended)

To respect the energy-conservation property the open-loop design protects:

- **Recommended: post-step projection (GGL-style).** After each RK4 step, project the
  state back onto the manifold by a few Gauss-Newton iterations: position-level
  `q ← q − G⁺ g(q)` until `g≈0`, then velocity-level `q̇ ← q̇ − G⁺ (G q̇)`. This removes
  drift without continuously injecting energy.
- **Alternative: Baumgarte** (`-Ġq̇ + 2α(Gq̇) + β²g` on the RHS) — simpler, one extra term,
  but injects/removes energy (the concern the docs raise). Keep as a fallback.

---

## 4. Module and API layout (separation)

New, self-contained headers — the open-loop files are **not modified** beyond the minimal
additive exposure in §5:

- `ga/ga_pga2dp_ops_constraints.hpp` — `closed_loop_system2dp`
- `ga/ga_pga3dp_ops_constraints.hpp` — `closed_loop_system3dp`

`closed_loop_system2dp` **has-a** `dynamic_system2dp` (composition, not inheritance) as
its spanning tree, plus a `std::vector<loop_constraint2dp>`:

```cpp
class closed_loop_system2dp {
    dynamic_system2dp tree;                  // the open-chain spanning tree (reused as-is)
    std::vector<loop_constraint2dp> loops;   // the extra closure edges
  public:
    // build the tree through the SAME familiar API (forwarded):
    void add_revolute_body(...);  void add_prismatic_body(...);  void add_body(...);
    // then close loops:
    void add_loop_constraint(loop_constraint2dp const&);
    // three-stage mirror of the open-loop tier:
    void assemble();              // static  : solve g(q)=0 (Newton)  -> consistent config
    void solve_velocities();      // kinematic: dependent rates from G q̇ = 0
    void step(value_t dt);        // dynamic : constrained KKT + projection + RK4
    // diagnostics reuse the tree's: kinetic_energy(), total_energy(), get_pos_trafo(), ...
};
```

Open-loop users include only `..._ops_physics.hpp`; closed-loop users include
`..._ops_constraints.hpp`. The split mirrors the existing optional-capability pattern.

---

## 5. Reuse surface — minimal, additive exposure from the base classes

The new layer needs the tree's mass matrix, bias forces, and partial velocities. The
*behaviour-preserving* refactor (does **not** change the open-loop logic or results):

- **DONE.** `forward_dynamics` is split into `assemble_mass_bias(rj) -> {M, rhs}` + a thin
  solve, in both 2dp and 3dp; the existing `step()` path is unchanged.
  `assemble_mass_bias` is **private for now** (open-loop API unchanged) and carries the
  bias-pass side effect (zeroes `rel_atwist`); a comment marks it as the reuse seam.
- **TODO when the layer lands.** Expose `assemble_mass_bias` (make public, or add a
  `friend class closed_loop_system{2,3}dp`); expose `velocity_field(S, P)` (partial
  velocity), `spatial_dot`, and the joint-screw builder
  (`move2dp(joint.screw_b, get_pos_trafo(j,0))`). These stay internal until a caller
  exists, to avoid widening the open-loop API speculatively.

If exposing the base internals is later judged undesirable, the new layer can instead
**recompute** `M`/`τ` from the already-public primitives (`get_pos_trafo`,
`point_acceleration`, `body.I`, `spatial_dot`) at the cost of some duplicated assembly.
The split makes the reuse path the natural choice (single source of truth for the
assembly).

---

## 6. Implementation phases

- **Phase 0 — expose (additive). [DONE 2026-06-13]** `forward_dynamics` split into
  `assemble_mass_bias` + a thin solve in both 2dp and 3dp (`<utility>` added,
  clang-formatted), documented as the closed-loop reuse seam. Validated unchanged:
  `ga_appl2dp_test` (42 cases / 497 assertions) and `ga_appl3dp_test` (28 / 7601) pass
  with the same energy tolerances (double pendulum `dE ~ 1e-13..1e-14`). Helper exposure
  (`assemble_mass_bias` public/friend; `velocity_field`/`spatial_dot`/screw builder) is
  deferred to Phase 1, when a caller exists. Prerequisite also fixed: a clang-22 doctest
  `__LINE__` collision that blocked the `ga_appl3dp_test` build
  (`DOCTEST_COUNTER=__COUNTER__` for all test targets in `ga_test/CMakeLists.txt`). No
  ga_py impact (verified): the split is a private method of an unbound stateful class —
  absent from `manifest.json` and `ga_py/src/generated/`.
- **Phase 1 — constraints + assembly (position level).** `loop_constraint2dp`, residual
  `g(q)`, Jacobian `G` (point-coincidence first). `assemble()` solves `g(q)=0` by Newton.
  Validate on a **planar four-bar** (assemble a consistent closed configuration).
  *ga_py note:* the new pure-data `loop_constraint{2,3}dp` struct (public-field aggregate)
  auto-binds on a plain bindgen regeneration — add a dedicated ga_py test for it (the
  coverage lists do not assert completeness). The stateful `closed_loop_system` class is
  NOT bound (reconstructed in Python from primitives, see Phase 5). Give the struct an fmt
  formatter in `ga/detail/fmt/` so the generated `__str__`/`__format__` work.
- **Phase 2 — kinematic closed loop.** Distribute driver joint rate to dependent rates
  (`G q̇ = 0`) and accelerations (`G q̈ = -Ġ q̇`). Validate four-bar velocity/accel ratios
  vs. analytic.
- **Phase 3 — dynamic closed loop.** KKT (or Schur) solve for `q̈, λ`; post-step projection
  stabilisation; RK4 via the shared `rk4_step`. Validate four-bar **energy conservation**
  to integrator tolerance (the headline correctness metric) and bounded constraint drift.
- **Phase 4 — 3D lift.** Mirror onto `closed_loop_system3dp` over the `pga3dp` base
  (point-coincidence → spatial spherical joints). Validate a **Stewart-Gough / delta**
  analogue: forward dynamics of a 6-leg platform; energy check.
- **Phase 5 — docs, demo, wrappers.** Extend `ga_docu` with the new facility as a
  first-class deliverable (scope + mandatory style match in §9); add a `ga_view` four-bar
  (2D) and Stewart/delta (3D) scene; optional Python exposure (the class is stateful, so
  reconstructed in Python from the bound primitives, as done for `kinematic_system2dp` in
  `ga_py/tests/test_merry_go_round.py`).

---

## 7. Testing and validation

Mirror the existing test files
([ga_test/src/ga_appl2dp_physics_test.hpp](../ga_test/src/ga_appl2dp_physics_test.hpp) and
the 3dp analogue):

- four-bar linkage (planar, 1-DOF closed loop): assembly, velocity ratios, energy
  conservation, constraint drift `‖g‖` bounded over a long run;
- a planar slider-crank (mixed revolute/prismatic loop);
- 3D: a minimal Stewart-Gough (or delta) — forward dynamics under gravity, energy check,
  `‖g‖` bound;
- regression: all existing open-loop tests still pass byte-identical after the Phase 0
  refactor.

---

## 8. Open design decisions (where input would refine the plan)

1. **Stabilisation default** — projection (recommended, energy-clean) vs. Baumgarte
   (simpler). Plan assumes projection; confirm.
2. **Solver form** — bordered KKT via `lu_solve` vs. Schur complement reusing `M⁻¹`.
   Either works; Schur reuses the existing path more directly.
3. **Maximal coordinates** — left as a future alternative, not in the first cut. Promote
   it only if fully general (non-tree-reducible) topologies are required.
4. **Constraint vocabulary** — start with point-coincidence only, or seed point-on-line /
   distance from the outset?

These are refinements, not blockers — the plan can proceed on the defaults above.

---

## 9. Documentation (ga_docu) — scope and required style

The closed-loop facility must be documented in `ga_docu` as a first-class part of the
modelling chapter, not an afterthought. The write-up is a deliverable of equal weight to
the code, and **must read as if written by the same hand as the existing documentation**.

**Where.** A new subsubsection in
[ga_docu/5_ga_modelling_physics.tex](../ga_docu/5_ga_modelling_physics.tex), placed after
"Articulated coupled dynamics" (the open-loop dynamics tier), cross-referenced from the
"Frame trees" overview and from the existing "Reduced vs. maximal coordinates" addendum
(which already foreshadows this work). Title in the descriptive register of the
surrounding sections, e.g. "Closed-loop and parallel mechanisms".

**What.** Mirror the development arc of the open-loop sections:

- the loop-closure problem and why the tree tier cannot express it (condense §1);
- the constrained equations of motion (the KKT system of §3.2) in the doc's own notation,
  with the constraint Jacobian shown as the reuse of the spatial-Jacobian columns (§3.4);
- the PGA-native constraint primitives (coincidence as a vanishing meet/join), §3.3;
- stabilisation (projection vs Baumgarte) and the energy argument, tied back to the
  reduced-vs-maximal discussion;
- a worked example carried through (planar four-bar, then the Stewart-Gough / delta lift),
  matching how the double pendulum is used in the open-loop sections.

**Style (mandatory — match the existing documentation).** The new prose must be
indistinguishable in voice from the hand-written sections:

- follow the **"Documentation prose voice"** rules in `CLAUDE.md` (terse, declarative,
  equation-first; no marketing/persuasive register; sparse `\emph{}`);
- use the project **notation conventions** (`CLAUDE.md` table: ⟑/⟇, `rev`/`rrev`,
  `move{2,3}dp`, `rcmt`, `rwdg`, bulk/weight) and the existing LaTeX macros (`\veedot`,
  `\rrev`, `\bulk`, `\weight`, …); cite Lengyel where the constraint geometry draws on it;
- sub-topics use **run-in `\paragraph{Label:}`** labels (not titled blocks), and never end
  the preceding paragraph with `\\` before a `\paragraph`;
- wrap prose at column 90; keep tables and equations on their own lines; reference
  equations with `\ref`/`\eqref` as the rest of the chapter does.

**Done when** the chapter builds (PDF recompiles, no LaTeX errors) and a read-through
shows the new subsection sitting seamlessly beside the open-loop material.

---

## 10. References (in-repo)

- open-loop tier: [ga/ga_pga2dp_ops_physics.hpp](../ga/ga_pga2dp_ops_physics.hpp),
  [ga/ga_pga3dp_ops_physics.hpp](../ga/ga_pga3dp_ops_physics.hpp)
- reduced-vs-maximal rationale + dynamics derivation:
  [ga_docu/5_ga_modelling_physics.tex](../ga_docu/5_ga_modelling_physics.tex)
- shared numerics: `detail/ga_solver.hpp` (`lu_solve`); `ga_usr_utilities.hpp`
  (`rk4_step`)
- moving-frame kinematics (`rcmt` velocity/accel fields):
  [ga_docu/3_ga_modelling_motion.tex](../ga_docu/3_ga_modelling_motion.tex)
