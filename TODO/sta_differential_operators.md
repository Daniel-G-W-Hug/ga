# STA Differential Operators and Electrodynamics (idea collection)

Status: **idea / not started** (2026-07-04). This note collects the plan for the missing
STA physics tier so it can be refined step by step before implementation. Origin: the
ga/ga_docu structure review identified this as the one substantive modelling gap — STA
has the full kinematic layer (boosts, general rotor `exp`/`log`/`sqrt`, causal character,
observer splits, `transform_opt`) but nothing that *uses* it dynamically, the way
chapter 5 (mechanics) uses chapter 3's PGA motion layer.

---

## Goal

A `ga_sta4ds_ops_electrodynamics.hpp` layer (name follows the mechanics precedent set by
`ga_pga{2,3}dp_ops_mechanics.hpp`) providing:

1. the spacetime **vector derivative** $\nabla = g^\mu \partial_\mu$ acting on
   multivector fields,
2. **Maxwell's equation** in GA form, $\nabla F = J$ (one equation replacing all four),
3. the **Lorentz force law** and relativistic particle motion, incl. the rotor equation
   of motion, integrated in proper time.

This gives ga_docu's "Motion in STA" (3.5) a physics counterpart and completes the
symmetry of the document (EGA → PGA mechanics; STA → electrodynamics).

## What already exists (foundation is complete)

- **Numeric substrate** (all ported/shared, zero new dependencies needed):
  - `ga/detail/ga_stencil.hpp` — FD stencil generator (`stencil_t`): weights, order,
    truncation error for explicit *and* compact schemes; built for the discretized
    $\nabla$ (added 2026-07-04, incl. two fixes vs the personal-lib source).
  - `ga/detail/ga_solver.hpp` — `lu_solve` / `lstsq_solve` / `kkt_solve`.
  - `ga/ga_usr_utilities.hpp` — `rk4_step` (mdspan + vector forms) for the proper-time
    ODE; step/easing helpers.
- **STA kinematics** (`ga_sta4ds_ops.hpp`): rotor sandwich `transform`/`transform_opt`,
  `get_rotor`/`get_boost` (mind the documented half-angle sign split!), general
  (non-simple) rotor `exp`/`log`/`sqrt`, causal character via the geometric square
  `B² = gr0(X⟑X)` (NOT `nrm_sq` — see CLAUDE.md "STA4D rotor operations").
- **Faraday-bivector observer split** (`ga_sta4ds_ops.hpp`): `rel_vec_split(F, u)` /
  `rel_bivec_split(F, u)` — electric-type / magnetic-type parts relative to a timelike
  observer u; under a boost the parts mix (the EM field transformation).

## Design sketch

### Two derivative paths — do NOT couple them

- **Closed-form path** (analytic fields): fields given as callables returning
  multivectors; $\nabla$ applied symbolically per case (plane wave, Coulomb, constant
  F). No stencils involved. This path drives the analytic test gates.
- **Discretized path** (grid fields): fields sampled on a (regular, later maybe
  stretched) grid; $\nabla$ assembled from `stencil_t` weights per axis. This is the
  consumer the stencil port prepared. Compute stencils in normalized coordinates (x/h)
  — the order-detection threshold in `stencil_t` assumes O(1) spacing (documented in
  the header).

### The vector derivative

$\nabla = g^\mu \partial_\mu$ with $g^\mu$ the reciprocal basis (STA metric
{-1,-1,-1,+1}, so reciprocal vectors carry signs). Acting on a multivector field A:

- $\nabla A = \nabla \cdot A + \nabla \wedge A$ (grade-lowering + grade-raising part);
- on a bivector field F: $\nabla F$ has vector + trivector parts →
  $\nabla \cdot F = J$ (sources) and $\nabla \wedge F = 0$ (no magnetic monopoles)
  are the two grade-parts of the single equation $\nabla F = J$.

Implementation idea: a small `field4ds<T>` view type (mdspan over grid + spacing +
per-axis stencil choice), plus free functions `nabla(field)`, `nabla_dot(field)`,
`nabla_wdg(field)` returning fields of the appropriate grade. Start with vector and
bivector fields only (what Maxwell needs); generalize later if needed.

### Maxwell in GA form

- $F$ = Faraday bivector (6 components = E (g_k4 part) + B (g_jk part) — exactly the
  `rel_vec_split`/`rel_bivec_split` decomposition relative to u = g4).
- $J$ = four-current vector (charge density + current density).
- $\nabla F = J$: verify grade-wise against the classical Maxwell equations for the
  standard observer. Gauge potential option: $F = \nabla \wedge A$ with vector
  potential A; then $\nabla^2 A = J$ in Lorenz gauge ($\nabla \cdot A = 0$) — the
  d'Alembertian $\nabla^2$ is a good first stencil consumer (pure second derivatives).

### Lorentz force and particle motion

- Force law: $\dot p = q \, (F \cdot u)$ — grade-1 contraction of the Faraday bivector
  with the four-velocity; $u \cdot u = 1$ preserved automatically (F·u ⊥ u).
- **Rotor equation of motion** (the elegant STA formulation):
  $\dot R = \tfrac{q}{2m} F R$, with $u(\tau) = R \, g_4 \, \tilde R$. Integrates with
  `rk4_step` in proper time; renormalize R via the versor norm (NOT `nrm_sq` — same
  trap as `sqrt(rotor)`, see CLAUDE.md).
- For constant F this has the closed form $R(\tau) = \exp(\tfrac{q\tau}{2m} F) R_0$ —
  directly reuses the general rotor `exp` (a boost·rotation product for generic F) and
  is the primary analytic gate.

## Analytic test gates (in the repo's gate-first style)

1. **Constant magnetic field**: cyclotron motion; frequency $\omega = qB/(\gamma m)$,
   circular orbit radius, $\gamma$ constant. RK4 vs closed-form rotor exp.
2. **Constant electric field**: hyperbolic motion (constant proper acceleration);
   closed form in rapidity, $\gamma(\tau) = \cosh(a\tau)$.
3. **Crossed fields E×B** (|E| < |B|): drift velocity $v_d = E/B$; boost to the drift
   frame kills E (`rel_vec_split` after `transform` — ties the EM field transformation
   to the particle dynamics).
4. **Plane wave**: $F = f \, \exp(I k \cdot x)$ null field; $\nabla F = 0$ closed-form;
   discretized $\nabla F \to 0$ at the stencil's reported order under grid refinement
   (reuses the convergence-rate test pattern from `ga_stencil_test`).
5. **Coulomb field**: $\nabla F = 0$ off the source; total flux ↔ charge (discrete).
6. **Wave equation**: $\nabla^2 A = 0$ plane-wave dispersion on the grid — pins the
   d'Alembertian stencil assembly (compact vs explicit stencil comparison possible).

## Open questions (to refine step by step)

- **Field storage**: mdspan-of-multivectors (AoS) vs per-component grids (SoA)? AoS is
  simpler and matches `rk4_step` usage; SoA vectorizes better. Start AoS, benchmark
  later (`ga_test/utilities/` bench pattern).
- **Boundary handling** for the discretized $\nabla$: one-sided stencils from
  `stencil_t` at grid edges (it generates them naturally) vs ghost cells vs periodic.
  Start periodic (plane-wave gates) + one-sided as the second option.
- **Time integration of fields** (FDTD-style $\partial_t F$ from $\nabla F = J$):
  needed for wave propagation demos; leapfrog/Yee vs RK4 — decide when the demo is
  chosen. Not needed for the particle-motion tier.
- **Units/scaling**: natural units (c = 1) throughout, matching the STA metric
  convention; document in the header like the sign conventions.
- **Binding scope**: closed-form particle motion (rotor ODE, gates 1-3) is a natural
  ga_py surface; grid fields probably stay C++-only initially (mdspan not bindable —
  same reason `rk4_step` has the vector-form twin).
- **ga_view**: a 2D projection demo (cyclotron / E×B drift trajectories) would fit the
  existing scene pattern; a real field visualization is a 3D-viz question →
  [vtk_3d_visualization_context.md](vtk_3d_visualization_context.md).
- **ga_prdxpr**: no new products needed (gpr/wdg/dot cover $\nabla$'s algebra); the
  reciprocal-basis signs are just the STA metric — but double-check the contraction
  identities used by `nabla_dot` against the sta4ds `.cases` when implementing.

## Suggested phasing

- **Phase E.1 — particle tier (no grids)**: Lorentz force + rotor EOM + proper-time
  RK4; gates 1-3. Smallest useful slice, zero new infrastructure.
- **Phase E.2 — closed-form field tier**: $\nabla$ on analytic fields, $\nabla F = J$
  grade-split, gauge potential; gates 4-5 (closed-form halves).
- **Phase E.3 — discretized tier**: `field4ds` grid view + stencil-assembled $\nabla$;
  gates 4-6 (discrete halves), convergence-order checks.
- **Phase E.4 — docu**: "Electrodynamics in STA" section in ga_docu (equation-first,
  matching the house voice), glossary entries (Faraday bivector, four-current, proper
  time), literature (Hestenes STA papers, Doran/Lasenby ch. 7).

(Phase letter "E" continues the wafer-grinding plan's lettering only by coincidence —
renumber freely when this becomes a real plan.)
