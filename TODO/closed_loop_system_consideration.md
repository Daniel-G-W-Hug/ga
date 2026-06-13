# Closed-loop systems — findings and implementation plan

Status: design note / proposal. **The full 2D AND 3D stack landed — Phase 0 (reuse-seam
prep), Phase 1 (C++ core + ga_py), Phase 2 (kinematic), Phase 3 (dynamic, energy-conserving
KKT) and Phase 4 (3D lift over `pga3dp`, with the numeric kernels shared via
`detail/ga_solver.hpp`)**; only Phase 5 (docs/demo/wrappers) remains. Captures
(a) why the current `static_/kinematic_/dynamic_system{2,3}dp` tier
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
- Phase 1 C++ core done: new additive header
  [ga/ga_pga2dp_ops_constraints.hpp](../ga/ga_pga2dp_ops_constraints.hpp) with
  `constraint2dp`, the pure-data `loop_constraint2dp` descriptor, and
  `closed_loop_system2dp` (composition: has-a `dynamic_system2dp` tree). Position level:
  `residual()` g(q), `constraint_jacobian()` G (reusing the `velocity_field(S_j, P)`
  spatial-Jacobian columns), and `assemble(driven)` — Newton solve for a consistent closed
  configuration (square LU on the dependent joints; min-norm / least-squares fallbacks for
  redundant / over-constrained loops). Reuse seam exposed via `friend class
  closed_loop_system2dp` on `dynamic_system2dp` (open-loop public API byte-unchanged);
  fmt formatter added in `ga/detail/fmt/ga_fmt_physics.hpp`; wired into `ga_pga.hpp`.
  Validated: a planar four-bar (`ga_appl2dp_physics_test.hpp`) assembles to an
  analytically-known closed configuration (θ3, θ4 to 1e-9; closure residual ‖g‖ ~1e-14;
  bounded ~1e-13 over a crank sweep). Open-loop suites unchanged: `ga_appl2dp_test`
  43/508, `ga_pga_test` 169/2758, `ga_appl3dp_test` 28/7601 all pass.
- Phase 1 ga_py done: `loop_constraint2dp` does NOT auto-bind on a plain regeneration — its
  `constraint2dp` enum field was filtered by the emitter (same reason `joint_state2dp`/
  `body2dp` were unbound). Fixed by adding **scoped-enum support to ga_bindgen**: `model.py`
  gains an `Enum` schema; `scan.py` collects namespace-scope `ENUM_DECL`s; `emit_nanobind.py`
  emits `nb::enum_` modules (`register_enums_<submod>`, called first in `register_all`),
  injects enum names into `type_map`, and switches the data-struct field ctor to C++20
  parenthesized aggregate init (so an `int`-bound `size_t` field like `frame_a` does not trip
  brace-init narrowing). This binds `constraint2dp`/`joint2dp`/`joint3dp` and the structs
  that use them (`loop_constraint2dp`, plus the bonus `joint_state2dp`/`joint_state3dp` — the
  latter needed a `twist3dp`→`bivec3dp` `type_map` alias and new fmt formatters in
  `ga/detail/fmt/ga_fmt_physics.hpp`). New `ga_py/tests/test_constraints.py` (13 cases);
  `test_stubs.py` updated to exclude `enum.Enum` classes from the `__format__` completeness
  check. Full chain regenerated (manifest, generated cpp, `.pyi`, cross-check JSON);
  `ga_py` suite 711 passed. The stateful `closed_loop_system2dp` stays unbound by design
  (reconstructed in Python from primitives, à la `dynamic_system2dp`).

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
- **Phase 1 — constraints + assembly (position level). [DONE 2026-06-13]**
  `loop_constraint2dp`, residual `g(q)`, Jacobian `G` (point-coincidence first).
  `assemble()` solves `g(q)=0` by Newton. Validated on a **planar four-bar** (assembles a
  consistent closed configuration against an analytic ground truth). ga_py exposure done
  (see progress note above — required adding scoped-enum support to ga_bindgen, since the
  `constraint2dp` enum field blocked the plain-regeneration auto-bind originally assumed
  here). The stateful `closed_loop_system` class is NOT bound (reconstructed in Python from
  primitives, see Phase 5).
- **Phase 2 — kinematic closed loop. [DONE 2026-06-13]** `closed_loop_system2dp` gains
  `set_joint_rate`/`joint_rate`, `solve_velocities(driven)` (`G_dep q̇_dep = -G_drv q̇_drv`)
  and `solve_accelerations(driven, driven_accels)` (`G_dep q̈_dep = -G_drv q̈_drv - Ġ q̇`,
  with the `Ġ q̇` velocity-product term read off as the relative anchor acceleration at
  dependent `q̈=0` via the kinematic layer's `point_acceleration` bias trick — no extra
  geometry). The position Newton, velocity and acceleration solves all route through one
  refactored `solve_jacobian_system(G,b)` (square LU / min-norm / least-squares). Validated
  on the four-bar at θ2=π/2 against analytic ratios (ω: coupler −4/3, rocker 1; α: coupler
  −2/9, rocker 2/3 for ω2=1, α2=0 — all exact to ~1e-12) plus velocity/acceleration closure
  of the shared coupler tip (`‖dv‖,‖da‖ ~1e-16`). Test:
  `ga_appl2dp_physics_test.hpp` (now 44 cases / 517 assertions). No ga_py impact (only
  member functions added to the unbound class).
- **Phase 3 — dynamic closed loop. [DONE 2026-06-13]** `closed_loop_system2dp` gains
  `joint_accelerations()` (the bordered acceleration-level KKT solve for `q̈, λ` — `M`/`τ`
  reused from `assemble_mass_bias`, `G` from `constraint_jacobian`, `-Ġq̇` from
  `constraint_bias` = the relative anchor accel at `q̈=0` via the `point_acceleration` bias
  trick; solved with the shared LU) and `step(dt)` (RK4 via the shared `rk4_step`, then GGL
  post-step projection: position by a min-norm Newton `assemble({})`, velocity by
  `q̇ ← q̇ − G⁺(Gq̇)` — drift control with no energy injection). Validated on the four-bar
  released under gravity from θ2=1.2 (π/2 is a radial-gravity zero-torque equilibrium —
  started off it): **energy conserved `dE/scale ≈ 1.3e-6`** over 3 s (KEmax 37.9, the
  KE↔PE exchange balances E0=18.29 exactly), closure drift `max‖g‖ ≈ 1e-12`. Test:
  `ga_appl2dp_physics_test.hpp` (now 45 cases / 521 assertions). No ga_py impact.
  (Slider-crank mixed revolute/prismatic loop from §7 still open as a secondary test.)
- **Phase 4 — 3D lift.**
  - *Step 0 — shared-kernel extraction. [DONE 2026-06-13]* Before mirroring, the two
    domain- and dimension-agnostic numeric kernels were promoted to
    [ga/detail/ga_solver.hpp](../ga/detail/ga_solver.hpp) (beside `lu_solve`, templated on
    T): `lstsq_solve(A, b, ncols)` (rectangular square / min-norm / least-squares solve) and
    `kkt_solve(M, G, f, g, n, m, λ*)` (bordered saddle-point solve). `closed_loop_system2dp`
    was refactored onto them (dropping its private `newton_step`/`solve_jacobian_system` and
    the inline KKT assembly); behaviour byte-identical (all 24 closed-loop assertions, 521
    total, unchanged). Decision (user-confirmed): generic numerics live in the **general
    numeric layer**, NOT a physics-shared header, so the planned electromagnetics extension
    reaches them without depending on mechanics (same rationale as `lu_solve`/`rk4_step`).
    The GA-mechanics machinery (`constraint_jacobian` via `velocity_field`, `constraint_bias`
    via `point_acceleration`, the orchestration) stays in `_ops_constraints.hpp`, mirrored
    2dp/3dp (parallel-class convention, not templated sharing).
  - *Main lift. [DONE 2026-06-13]* `closed_loop_system3dp` in
    [ga/ga_pga3dp_ops_constraints.hpp](../ga/ga_pga3dp_ops_constraints.hpp) — a faithful
    mirror of the 2D class over the `pga3dp` base (point-coincidence → spatial spherical
    joint; residual = 3-component unitized anchor difference; `G` columns = `velocity_field`
    on `vec3dp`; `twist3dp`/`move3dp`/`mvec3dp_e` swaps; the `(pivot,axis)` revolute
    signature), consuming the shared `lstsq_solve`/`kkt_solve` verbatim. Friend hook +
    forward-decl added to `dynamic_system3dp`; `loop_constraint3dp` + `constraint3dp` fmt
    formatter; wired into `ga_pga.hpp`. **Validated** on a genuinely-spatial 1-DOF two-arm
    reacher (two 2R arms, shoulder-yaw e3 + elbow-pitch e2, forearm tips pinned — full-rank
    3D constraint; a planar embedding would make the KKT singular): assembly `g 0.385 →
    5e-13`, velocity closure of the shared hand `~1e-12`, and **energy conservation
    `dE/scale ≈ 9e-10`** with closure drift `~1e-12` over a 2 s spatial swing.
    `ga_appl3dp_physics_test.hpp` (now 29 cases / 7609 assertions). ga_py: the
    `loop_constraint3dp` struct and `constraint3dp` enum auto-bound by the Phase-1 enum
    support (pga now 34 types / 4 enums);
    `ga_py/tests/test_constraints.py` extended (715 ga_py tests). The Stewart-Gough / delta
    6-leg platform remains a richer future demo (Phase 5), but the 3D layer + energy metric
    are proven here.
- **Phase 5 — docs, demo, wrappers.** Extend `ga_docu` with the new facility as a
  first-class deliverable (scope + mandatory style match in §9); add a `ga_view` four-bar
  (2D) and Stewart/delta (3D) scene; optional Python exposure (the class is stateful, so
  reconstructed in Python from the bound primitives, as done for `kinematic_system2dp` in
  `ga_py/tests/test_merry_go_round.py`).
  *ga_view side-by-side open- vs. closed-loop demo:* beyond the per-mechanism scenes above,
  add ONE `ga_view` 2D scene that shows an **open-chain** tree system and a **closed-loop**
  system next to each other, animated together, to make the structural difference legible.
  Left panel: an open chain built straight from `dynamic_system2dp` — e.g. a 2-3 link robot
  arm, or a simplified stick-figure person (torso → upper/lower-arm and upper/lower-leg
  branches off a common root, the existing branching tree). Its end-effector / hand swings
  freely (the chain has no closure). Right panel: a `closed_loop_system2dp` whose spanning
  tree is the *same* open chain plus an `add_loop_constraint` that pins the end-effector
  (e.g. a four-bar, or the arm with its hand pinned to a fixed world point / to the other
  hand), so the same driver motion now propagates through the loop and the dependent joints
  are solved. Drive both from one shared parameter (a crank angle / a joint sweep) so the
  viewer sees, in lock-step: the open chain free vs. the closed loop's coupled, constrained
  motion; optionally overlay the closure point and the constraint residual `‖g‖`. This is
  the headline visual for the "tree vs. tree+closure" distinction the whole note is about
  (§1 vs. §3) — keep ga_view's left-handed screen convention local per
  `[[feedback_library_docs_right_handed]]`.
  *ga_py reconstruction test (enabled by the Phase 1 enum work):* the Phase 1 ga_py change
  also bound the joint enums (`joint2dp`/`joint3dp`) and the joint descriptors
  (`joint_state2dp`/`joint_state3dp`) — currently with **no test consumer**. Add a
  `dynamic_system2dp` reconstruction test that exercises them: reproduce the joint-space
  forward dynamics (`assemble_mass_bias` → `M(q) q̈ = τ` → RK4) in Python from the bound
  primitives, modelling each joint faithfully via `joint2dp.revolute`/`prismatic`, the
  body-frame screw `screw_b`, and the rest motor — exactly the way `test_merry_go_round.py`
  uses the bound `pose2dp` rather than a hand-rolled pose tuple. Validate against the C++
  open-loop result (e.g. the double-pendulum energy / mass-matrix identity, or the four-bar
  closed-loop assembly). This is the payoff that turns the now-bound joint types from
  available-but-unused into a tested, faithful Python mirror of the dynamic layer, and is
  the natural stepping-stone to a Python `closed_loop_system2dp` reconstruction.

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
