# Wafer-grinding dynamics — modelling & optimization plan

Status: **planning + initial port done (2026-06-14).** Long-horizon project; this file
evolves and is the living plan across sessions. The kinematics/dynamics infrastructure
(`static_/kinematic_/dynamic_system{2,3}dp`, the closed-loop layer) was built as the
runway for exactly this. Companion memory: `project_wafer_grinding`.

**DONE so far:** the two grinding mini-tests translated onto `static_/kinematic_system3dp`
(pedagogical, in `ga_appl3dp_appl_test.hpp`, right after the originals); the
`is_congruent` refinement (intermediate step below); the Fig. 1 **machine geometry &
kinematic frame tree** specified (new section + Phase 0); **Phase 0.a + 0.b DONE &
visually confirmed** (grain- trajectory app-test passes; `active_grinding_marks` ga_view
scene with both views, contact masking, and C-key ratio cycling — see Phase 0; user
committing this state).

**Phase B.1 DONE 2026-06-14.** Sommerfeld explicit-forcing warm-up landed as app-test
`"pga3dp: Sommerfeld unbalanced-rotor warm-up (Phase B.1)"` in `ga_appl3dp_appl_test.hpp`
(right after the Phase-0 grain-trajectory case). Reproduces the steady-state analysis of
Bisoi et al. (2020): the spindle CM is an unbalanced rotor of total mass `m+M` on an
anisotropic foundation (springs Kx,Ky + dampers Rx,Ry along e2/e3) with PRESCRIBED spin;
the rotating unbalance `m·e·Ω²(cosΩt e2 + sinΩt e3)` decouples the two radial axes into the
paper's **Eq. (2)** forced oscillators. Net radial force built as a GA `vec3dp` direction
(the wrench bookkeeping the library force elements will reuse), RK4 via the shared
`rk4_step`/`rk4_get_time`. **Gate met:** swept Ω∈[6,40], the RK4 steady-state amplitude
matches the **closed-form Eq. (4)** `A=m·e·Ω²/√((ΩR)²+(K−(m+M)Ω²)²)` to **<0.03 %** on both
axes (1 % epsilon gate); resonance peaks emerge at `√(K/(m+M))` = 16.50/28.57 rad/s (paper
16.51/28.57). Uses the **actual Table-1 params** (m=4.9, M=2.45, e=0.008336, Kx=2000,
Ky=6000, Rx=5, Ry=10) → reproduces **Fig. 2** (peak Ax≈0.135 m, Ay≈0.117 m). 33 cases/9143
assertions pass. Scope: pure-prototype (no library change yet); GA content thin by design.
The Sommerfeld **jump** (Figs. 3–10) needs the two-way non-ideal DC-motor coupling (power
balance Eq. 10 / bond graph) — DEFERRED; prescribed spin gives the resonance backbone only.
Restore point before this: `d25ebe6`.

**RESUME HERE (next session) →** pick one of:

1. **Phase B.2** — replace B.1's *explicit* `m·e·Ω²` forcing with an **emergent** spinning
   offset-cm rigid body (the unbalance arises from a body whose cm is offset from the spin
   axis), and reproduce B.1's swept amplitudes. Validates "`meω²` for free" + the
   prescribed-spin bias path; the natural next step on the warm-up. Still a `ga_appl3dp`
   prototype.
2. **Phase A (library promotion)** — fold the proven spring/damper generalized force +
   applied wrench into `dynamic_system3dp::tau` (additive; gravity path byte-unchanged), so
   `step()` can integrate forced systems. Do this once B.1/B.2 have proven the primitive.
3. **Phase 0.c** — the two-pane top (−e1, `e423`) + side (−e2, `e431`) CS-relative-motion
   view (deferred). Useful mainly once dynamics add radial/axial/tilt; can wait.

**Recommendation: B.2** (closes out the warm-up's GA payoff), then A (promote), then C.

**RESOLVED (2026-06-14): static-vs-feed = PRESCRIBED INFEED.** The spindle carries a
steady prescribed translation (axial infeed `x_a` along −e1 and/or radial traverse) ON TOP
of the 5 vibration DOFs — so `x_a` (and possibly the radial offset) is a time-varying
prescribed input, not a constant. Phase C must drive the nominal pose with this prescribed
motion; the 5 vibration DOFs perturb around it. Phase 0 / B.1 still run at the static
nominal pose (the warm-up validates the spring/damper primitive, independent of feed).
**Minor/optional:** silence the app-wide Qt "Sans Serif" font warning via a
default font in `ga_view/src/main.cpp` (`QApplication::setFont(QFont("Helvetica", …))`).

## Goal

Model the geometric **and** dynamic aspects of **wafer self-rotational grinding** to
understand the grinding process and — visionary goal — **optimize process parameters**
(feed rate of the grinding wheel relative to the chuck, speed ratio) to dynamically
control the grinding force and **prevent wafer damage**.

Target system (Tao et al., *Int. J. Mech. Sci.* 232 (2022) 107620): a cup grinding wheel
on an **aerostatic-bearing spindle** grinds a wafer held on a rotating **chuck**. The
spindle vibrates in radial (x,y), axial (z) and **tilting** (θ,φ) DOFs (5-DOF model, mass
eccentricity `e` → centrifugal `meω²`, gyroscopic tilt coupling). The axial runout `z_b =
z − R_w·φ` (Eq. 14) modulates grain cutting depth → surface waviness (WMD along the mark
direction, WCD circumferential; wavelengths λ_m, λ_c). A **real machine exists** with
**measured time- and frequency-domain vibration data** at several points, for calibration
and validation.

## Starting point (already in the repo)

Two prototype tests in `ga_test/src/ga_appl3dp_appl_test.hpp` (pre-infra, raw PGA):

- *"reference and tumbling plane"* — wheel **tilt kinematics**: a wobbling tool normal,
  the tool∩wafer intersection line and its angle. = Tao's (θ,φ); `z_b` is just a rim
  point's axial coordinate under a motor sandwich, no separate equation.
- *"intersecting discs"* — **contact/removal kinematics**: tool ring vs. wafer disc,
  overlap MRR volume, and the relative **surface speed** `rs·B_uv >> (P − cp)`. = Tao's
  Eq. (26) grain velocity.

The geometry/observation half is ~prototyped; the missing half is the **dynamics** driving
it. **Translated** versions of both (expressing the same idea via
`static_/kinematic_system3dp`) sit right after the originals — pedagogical, no new
insight, asserted congruent to the raw result.

## Machine geometry & kinematic frame tree (Tao Fig. 1)

The geometry of Fig. 1 maps directly onto a `static_/kinematic_system3dp` frame tree: a
stationary inertial root with **two independent chains** branching off it — the
chuck/wafer chain and the spindle/wheel chain. Each frame is posed RELATIVE to its parent
(`static_frame3dp(name, origin_in_parent, rot_axis*angle)`), and spin is the
parent-relative angular velocity carried in `kin_state3dp{.omega = ...}`. This is the
geometric foundation for everything downstream: the grain-trajectory curves, mini-test 1
(tilt), and the Sommerfeld warm-up all build on this exact tree.

**Root frame `chuck_ctr_stat`** (inertial reference, identity in the library): sits on the
chuck axis at the chuck surface that holds the wafer. Physical axis labels: `e1` (x) = up
= chuck-surface normal (the chuck/wafer spin axis); `e3` (z) = to the right; `e2` (y) = to
the front, toward the observer. Right-handed: `e1 x e2 = e3` (up x front = right). All
other frames are expressed relative to this root.

**Frame tree** (origin/rot are RELATIVE to the parent; `e1,e2,e3` are the parent's axes):

| Frame                   | Parent             | Origin in parent                                                              | Orientation `rot` (axis·angle)  | Spin `omega`    | Role                                                                         |
|-------------------------|--------------------|-------------------------------------------------------------------------------|---------------------------------|-----------------|------------------------------------------------------------------------------|
| `chuck_ctr_stat`        | (self / root)      | `O`                                                                           | `0`                             | —               | inertial reference; chuck surface                                            |
| `chuck_ctr_rot`         | `chuck_ctr_stat`   | `O` (coincident at t=0)                                                       | `0`                             | `n_s · e1`      | rotates with the chuck about its axis (e1)                                   |
| `wafer_top_avg_rot`     | `chuck_ctr_rot`    | `(tw_avg, 0, 0)` (+e1)                                                        | `0` (co-rotates)                | —               | average top surface of the wafer; X-deviation = local profile                |
| `spindle_cm_stat`       | `chuck_ctr_stat`   | `(x_a, −R/√2, +R/√2)` (radial offset \|·\|=R; `x_a` = axial infeed along −e1) | `(0, -pi/2, 0)` (−90° about e2) | —               | spindle centre of mass; nominal Fig.-1 placement + 5 vibration DOFs vs. root |
| `tool_top_avg_rot`      | `spindle_cm_stat`  | `(0, 0, l3)` (+ local z)                                                      | `0` (co-rotates)                | `n_w · z_local` | average grinding-tool surface; rotates with the wheel about its axis         |
| `tool_surface_avg_at_R` | `tool_top_avg_rot` | `(R_w, 0, 0)` (+ local x)                                                     | `0` (co-rotates)                | —               | a point on the wheel's outer rim (radius R_w)                                |

Notes / parameters to pin down:

- **Spin naming (per Fig. 1, counterintuitive):** Fig. 1 labels the **chuck** spin `n_s`
  (~200–300 rpm) and the **grinding-wheel** spin `n_w` (~3000 rpm) — i.e. the `_s`
  subscript sits on the chuck and `_w` on the wheel, opposite to what the letters suggest.
  This doc follows Fig. 1. (The paper text may differ; if so, reconcile then.)
- **Nominal radial placement (Fig. 1, fixed for now):** `spindle_cm_stat` is offset from
  the chuck centre by `−R/√2` along root e2 (Y) and `+R/√2` along root e3 (Z), with `R =
  R_wheel = R_wafer` (assumed equal for now). The horizontal offset magnitude is then
  exactly **R**, so the wheel rim (radius R) reaches from the **wafer centre** (R−R) out
  to 2R and sweeps one full wafer radius — the standard self-rotational grinding layout:
  the rim just touches the wafer centre.
- **Axial infeed `x_a` (OPEN — needs paper check):** along the spindle's +z_local (= −e1,
  toward the wafer); with no tilt/shift the wheel first contacts the wafer top at the
  geometric event `tool-surface reaches the wafer` (≈ when the +l3 tool surface meets
  x=0). **Open question (user to clarify against Tao):** is the setup *static* (fixed
  nominal pose, only the dynamic vibration disturbances), or is there an *additional
  prescribed translational motion* of the wheel — a steady axial infeed and/or a radial
  traverse — beyond the disturbances? This decides whether `x_a` (and possibly the radial
  offset) is a constant or a time-varying prescribed input. **Resolve before Phase B/C
  dynamics; Phase 0 (grain trajectories) can run with the static nominal pose either
  way.**
- The 5 vibration DOFs (radial x,y; axial z; tilt θ,φ) are perturbations ON TOP of
  whichever nominal (static or fed) pose the above resolves to.
- `spindle_cm_stat`'s −90° about `e2` makes its local axes `x→e3`, `y→e2`, `z→−e1`: the
  spindle/tool axis (`+z_local`) points DOWN toward the wafer (−up), as for a cup wheel
  grinding the top surface. l1, l2 (spring/damper stations along the axis) are positions
  in this frame, deferred until the dynamics.
- Time-varying parameters (`tw_avg` as the wafer thins, `l3` as the tool wears) change
  slowly — treat as constant within ~10 chuck revolutions for now.

**Projection planes for ga_view** (results are 3DP; ga_view is 2D, so project):

- **e423_3dp** — the plane ⊥`e1` (carrier `e2∧e3`): the plane the chuck/wafer rotates in.
  Grain-trajectory curves (the Fig. 1 grey-dashed / Fig. 7 red-dashed marks) are shown
  here.
- **e431_3dp** — the plane ⊥`e2` (the x–z side elevation): shows the relative orientation
  of wafer and grinding tool.

**Grain-trajectory reproduction (precursor, pure kinematics — see Phase 0).** A grain
fixed on the wheel rim is `tool_surface_avg_at_R`'s origin; its curve in the **rotating
wafer frame** is `get_pos_trafo("tool_surface_avg_at_R", "wafer_top_avg_rot")` applied to
`O` sampled over time, then projected onto `e423_3dp`. For given `(n_w, n_s, R_w, wafer
radius)` this reproduces the epicyclic grinding-mark patterns of Fig. 1 / Fig. 7 —
achievable before any dynamics, and a direct functional check of the whole frame tree.

## The one capability gap: FORCE ELEMENTS

`dynamic_system{2,3}dp` has only gravity + the Coriolis bias as generalized forces, and
only **rigid** joints. Both the target and the warm-up need **springs + dampers to ground,
rotating unbalance `meω²`, and an applied wrench** (grinding force; motor torque).

Clean mapping onto the existing infra: model the spindle as **one rigid rotor body
connected to ground by a serial stack of single-DOF joints** — 3 prismatic (x,y,z), 2
revolute (θ,φ), 1 driven revolute (spin Ω) — **inertia on the final body**, a linear `−k·q
− c·q̇` force on each mount joint, and eccentricity `e` = the rotor's **cm offset**. Then
`assemble_mass_bias` already yields `M(q)` + the gyroscopic/Coriolis bias; we only add the
spring/damper and applied-wrench terms to `tau` and integrate with RK4.

**Two payoffs for free:** the centrifugal `meω²` forcing emerges from a spinning offset-cm
body; the gyroscopic tilt coupling `(J_y−J_z)ω φ̇` (Tao hand-derives in Eq. 13) emerges
from rigid-body Euler dynamics. **Honest caveats** (not geometry — must be added as force
elements): bearing reactions, unbalance reaction, the grinding-force law; the DC-motor
power balance is outside GA; the stiff bearing springs (Kz ~ 1e6 N/m) make the ODE
**stiff** → RK4 needs a small dt or an implicit integrator.

## Intermediate library step: refine `is_congruent` [DONE 2026-06-14]

Surfaced while porting mini-test 1: comparing two **independently computed** planes/lines
(results of wedge / dual / motor chains) with `is_congruent` fails at the default
`tolerance = eps`, even when the operands are identical to ~5 decimals. Two compounding
causes:

1. **Default tolerance is single-op machine eps**, but congruence is almost always applied
   to *computed* geometry whose error has accumulated over many GA ops. (CLAUDE.md's
   congruence note already says "absolute eps is too strict after GA ops that accumulate
   error" — that reasoning should extend to the **default**, not only the relative
   scaling.)
2. The algorithm fixes the scale factor `k` from the **first non-zero** component; when
   that component is small (e.g. a 0.01 tilt term next to a ~1 term), its relative
   rounding is amplified onto the dominant components and lands right at the eps boundary.

**Refinement [DONE]:** all 11 graded overloads (EGA2D/3D, PGA2DP/3DP, STA4D) now delegate
to a shared `hd::ga::detail::coeffs_congruent<N>` in `ga/ga_value_t.hpp` that (a) pivots
`k` on the **largest-magnitude** component of `b`, and (b) defaults to the new
`eps_congruent` (`1e-12`) relative tolerance reflecting accumulated GA error
(scalars/pscalars keep their trivial overloads + `eps`). The translated grinding test
reverted to plain `is_congruent(...)`. All 5 C++ suites + 725 ga_py tests pass; CLAUDE.md
congruence section updated.

## Plan (phased) — decisions taken 2026-06-14

- **0 — geometry & kinematic frame tree [precursor, no dynamics].** Build the Fig. 1 frame
  tree above as a `kinematic_system3dp` (two chains off `chuck_ctr_stat`). First
  functional payoff: reproduce the **grain-trajectory** curves (Fig. 1 grey-dashed / Fig.
  7 red-dashed) — the rim grain's path in the rotating wafer frame, projected onto
  `e423_3dp`, for a few `(n_s, n_w, R)` combinations. The Sommerfeld warm-up (B) and the
  5-DOF spindle (C) reuse THIS tree, so it is the on-ramp to the complex system.
  - **0.a [DONE 2026-06-14]:** app-test `"pga3dp: grinding grain trajectory in the wafer
    frame (Phase 0)"` in `ga_appl3dp_appl_test.hpp` builds the full six-frame tree and
    samples the rim grain over one chuck revolution. Gates pass: wheel axis at `|C|=R`,
    grain on the rim (`|grain−C|=R`), radial sweep `r ∈ [0, 2R]` (rim through wafer
    centre), `r` frame-invariant (root == wafer), rim speed `= n_w·R`. Emits the
    `e423_3dp` (y,z) projection column for ga_view. Used `n_s=300, n_w=3000 rpm` (ratio 10
    → pattern closes in one chuck rev), `R=150 mm`.
  - **0.b [DONE 2026-06-14, pending visual check]:** ga_view scene `active_grinding_marks`
    (the 5-place hand-sync surface) runs the six-frame `kinematic_system3dp` and draws the
    rim grain traced in the rotating wafer frame, projected onto `e423_3dp` (top view down
    −e1; e2 horizontal, e3 vertical). Fixed wafer disk + orbiting wheel disk + orange
    grinding-mark rosette; SPACE/R/T = pause/reset/trace. Normalized `R=1`, `n_s:n_w =
    1:10`. Builds and smoke-launches; **user to page to the scenes for the visual-layout
    check** (legend placement, colours, pace). Two scenes share one
    `active_grinding_marks` (a `gm_view` param picks the display frame). Marks are stored
    in **wafer-local coords + a contact flag** and drawn only while the grain is on the
    wafer (`r ≤ R`) — short strokes, not full loops. Display transform =
    `get_pos_trafo(wafer → display-frame)`: (1) **wafer frame** — marks sit still (the
    carved pattern), wheel orbits; (2) **global frame** — the same marks rotate rigidly
    with the wafer (wheel fixed, new strokes cut at the contact region). **C** cycles the
    n_w/n_s ratio (10.0 / 10.3 / 9.7) so off-integer ratios precess instead of repeating
    (shown bottom-left). Restore point before this refinement: commit `6a074e1`.
  - **0.c [proposed, deferred — decide after 0.b check]:** a second ga_view scene showing
    the relative motion of the *coordinate systems* (frame glyphs, not just the grain) in
    two projected views side by side — a **top view** along root −e1 (the `e423_3dp` chuck
    plane, as in the merry-go-round) and a **side view** along root −e2 (the `e431_3dp`
    elevation). Purpose: visualize the radial / axial / tilting spindle motions once the
    dynamics add them. Needs normalized rpm + sensible scaling, and is also practice at
    drawing projected views. Shift to a later session if 0.b suffices for now.
- **A — force-element extension.** Linear spring + damper generalized forces on joint
  coordinates, plus an applied external wrench at a body-fixed point, folded into `tau`.
  Additive; the gravity-only path stays byte-unchanged. (Mirrors how the closed-loop layer
  was added.)
- **B — Sommerfeld warm-up [CHOSEN to do first].** Unbalanced rotor on anisotropic x,y
  spring/damper supports, **prescribed spin**; validate steady-state amplitudes against
  the **closed-form Eq. (4)** of the Sommerfeld paper (Bisoi et al., *Advances in Rotor
  Dynamics*, 2020) and the resonances `ω₁,₂ = √(K/(m+M))` (= 16.51 / 28.57 rad/s for the
  paper's Table-1 params). The analytic correctness gate, exactly like the four-bar
  angles.
  - **B.1 [DONE 2026-06-14]:** explicit rotating unbalance forcing `m·e·Ω²(cos,sin)` +
    anisotropic radial springs/dampers, prescribed spin → matches the closed-form forced
    response to <0.05 % over an Ω sweep, resonances at `√(K/mt)` = 16.51 / 28.57 rad/s.
    App-test `"pga3dp: Sommerfeld unbalanced-rotor warm-up (Phase B.1)"`; force built as a
    GA `vec3dp`, integrated with the shared `rk4_step`. Isolates the spring/damper
    primitive. (See the DONE note at the top.)
  - **B.2 (next concrete step):** replace the explicit forcing with an **emergent**
    spinning offset-cm body → must reproduce B.1 (validates "`meω²` for free" + the
    prescribed-spin bias path).
- **C — Tao 5-DOF spindle.** The 6-joint stack above; reproduce Eq. (13) responses and the
  Fig. 4 characteristic frequencies `f_x, f_z, f_θ, f_φ`; derive `z_b` directly from a rim
  point; **calibrate/validate against the machine's measured data**. Port mini-test 1
  (tilt now driven by the dynamics, not a prescribed φ-sweep).
- **D — grinding contact + two-way coupling [the actual goal].** Promote mini-test 2 to a
  contact/observation layer; add a grinding-force law (MRR × relative-speed, calibrated to
  measured data) → wrench on the spindle; close the loop vibration → `z_b` → depth → MRR →
  force → vibration. Observe surface waviness (λ_c, λ_m) vs. measured topography.
- **E — optimization [visionary].** With the calibrated closed-loop model, optimize feed
  rate / speed ratio to bound the grinding force and `z_b` (damage prevention).

### Decisions

- **Non-ideal DC motor / Sommerfeld speed-jump (two-way drive↔structure): DEFERRED** until
  the grinding-force feedback (D) needs a reacting drive. Warm-up uses prescribed spin.
- **Code home: PROTOTYPE the force elements in a `ga_appl3dp` app-test first**, validate
  vs. Eq. (4), THEN promote the proven primitive into the library. (Matches how the two
  mini-tests started.)
- **Warm-up params:** use the Sommerfeld paper's Table 1 (they match the closed-form
  gate); bring the machine's measured stiffness/damping/eccentricity in at Phase C.

## Validation gates (per phase)

- 0: grain-trajectory curves match the Fig. 1 / Fig. 7 patterns for the paper's `(n_w,
  n_s, R_w, wafer radius)`; the frame tree round-trips (`get_pos_trafo` LCA walk) and a
  rim point's surface speed equals `n_s · R_w` (sanity, as in the translated mini-test).
- B: simulated amplitudes == Eq. (4) across an ω sweep; resonances at 16.51 / 28.57 rad/s.
- C: characteristic frequencies `f_x, f_z, f_θ, f_φ` and the `z_b` waveform vs. measured
  data.
- D: waviness wavelengths λ_c, λ_m vs. measured wafer topography.

## Out of scope (for now)

The full DC-motor/bond-graph drive model; the grain-statistics surface-topography
simulation as a standalone renderer (we observe `z_b`/trajectories, not a full grain
ensemble); any 3D `ga_view` rendering of the spindle (ga_view is a 2D viewer).
