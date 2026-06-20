# Wafer-grinding dynamics — modelling & optimization plan

> **PROJECT CLOSED 2026-06-20 (by user decision).** The work is concluded here; the
> remaining parts are **ABANDONED, not deferred** — they will not be picked up unless the
> project is explicitly reopened. **Done & committed:** Phases 0, A, B.1, C (Tao 5-DOF
> spindle) + the Eq.13 paper-error analysis, D.1 (surface formation), the whole D.2
> contact/force-loop tier (D.2a feed-forward → D.2b closed loop → D.2c thinning + chuck
> reaction → D.2d RK4/ABM2 integrators + stiffness verdict + adaptive dt → D.2e feed-control
> force limiting), `get_disc_inertia`/`make_disc_body` (2D+3D), and the ga_docu Wafer Grinding
> subsubsection (`6_ga_applications_pga.tex`). **Abandoned:** E (feed/speed optimization), the
> chip-thickness force-law refinement, D.2d-4 (implicit/stiff solver), Phase 0.c (CS-view
> ga_view scene), the non-ideal DC-motor / Sommerfeld speed-jump coupling, and the 2D twin of
> the contact element (the contact element stays 3D-only by design). The library
> infrastructure that landed (force elements, the swappable grinding-force-law callback, the
> integrator selector, the feed controller) remains in place and usable; only the further
> modelling/optimization work is dropped.

Status: **CLOSED 2026-06-20 — see the box above.** (Historical: planning + initial port done
2026-06-14.) The kinematics/dynamics infrastructure
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

**Phase A.1 DONE 2026-06-14 (spring/damper force elements, UNCOMMITTED).** Folded a linear
spring + damper on each 1-DOF joint coordinate into `dynamic_system3dp::tau`
(`ga/ga_pga3dp_ops_physics.hpp`): `tau_j += -k_j(q_j - q0_j) - c_j q̇_j`. Changes:
`joint_state3dp` gains `stiffness`/`damping`/`q_rest` (default 0 → gravity/bias path
byte-unchanged); public `set_joint_spring_damper(idx,k,c,q0=0)`; additive diagonal RHS term
in `assemble_mass_bias` (so coupled forward dynamics + `step()` pick it up automatically);
spring potential `½k(q-q0)²` added to `potential_energy()`. **Gate met:** new app-tests in
`ga_appl3dp_physics_test.hpp` + `ga_appl2dp_physics_test.hpp` (M2 suites) — a prismatic
spring/damper **damped harmonic oscillator** matches the closed form to **1.8e-12**, and the
**undamped** rig conserves KE+spring-PE to **2.7e-14** (identical in 2D and 3D). **2D twin
mirrored** (same 4 edits in `ga_pga2dp_ops_physics.hpp`). **ga_py regenerated** (scan →
emit_nanobind → build → stubs → test-data): `joint_state{2,3}dp` now expose
`stiffness`/`damping`/`q_rest`; the positional-ctor hand-sync (3 reconstruction/roundtrip
call sites + roundtrip field asserts) updated. **All green:** appl2dp 47/532, appl3dp
35/9148, pga core 169/2758, **ga_py 725**.

**Phase A.2 DONE 2026-06-14 (applied wrench + time threading, UNCOMMITTED).** Added a
time-varying applied external wrench folded into `tau` + RK4 sub-step time threading, in
both `dynamic_system{2,3}dp`. Per-frame sparse `wrench_` map (`std::function<bivecNdp(value_t
t)>`, world frame) + a `time_` clock; public `set_applied_wrench(idx, fn)`, `time()`,
`set_time()`. Fold-in in `assemble_mass_bias`: for each applied wrench W on frame fi,
`RHS[j] += spatial_dot(S[j], W)` for joints j supporting fi (the same reciprocal pairing as
gravity; sign verified +). `coupled_step` sets `time_ = rk4_get_time(t0,dt,s-1)` per stage
(restores on exit); `step()` advances `time_ += dt`. **Gate met:** new app-tests reproduce
B.1's Sommerfeld forced response INSIDE the library — single prismatic joint (mass m+M) +
`set_joint_spring_damper` + rotating `m·e·Ω²` wrench `wdg(O,F(t))`; steady amplitude ==
closed-form Eq.(4) to **≤0.05 %** (x@om1 0.13475 vs 0.13476; y@om2 0.11670 exact), plus a
constant-force `joint_accel == F0/mt` sign check. 2D twin mirrored (forced oscillator, same
numbers). **NO ga_py impact** (dynamic_system is unbound; no bound-struct change) — confirmed
ga_py 725 still green. **All green:** appl2dp 48/535, appl3dp 36/9152, pga core 169/2758,
ga_py 725. Phase A (force-element tier) is now COMPLETE: spring/damper + applied wrench +
time, 2D+3D, validated against the analytic gate.

**Phase A.3 DONE 2026-06-14 (driven/prescribed-rate joint, UNCOMMITTED).** A 1-DOF joint
can now be KINEMATICALLY DRIVEN at constant rate `q(t)=q0+rate·t` (motor spin / steady
feed), in both `dynamic_system{2,3}dp`. Design: stored in a side `driven_` map (NOT in the
bound `joint_state`, so **no ga_py regen**); a joint is driven iff it's in the map — keeps
its revolute/prismatic screw, is excluded from `dof_joints`, and is re-evaluated from `q(t)`
at each RK4 sub-step (via the A.2 clock). It acts as a **moving base**: its velocity
propagates into the dynamic sub-chain's Newton-Euler bias, so centrifugal/Coriolis emerge.
API: `set_driven_rate(idx,rate,q0=0)`, `clear_driven_joint`, `is_driven`. (Constant rate →
zero joint accel; general `q(t)` with `q̈≠0` is a future extension.) **Gate met:** 3D
kinematics test (`phi==Ω·t`, `|v_P|==Ω·L` exact) + 3D & 2D **spinning radial slider**: a
mass on a spring spun by a driven base settles to the centrifugal equilibrium
`q_eq=mΩ²L/(k−mΩ²)` → **0.24999 vs 0.25000** (centrifugal emerges, nothing added by hand).
All green: appl2dp 49/536, appl3dp 38/9156, pga 169/2758, ga_py 725. CLAUDE.md gained a note
on the bound-struct field-addition ga_py hand-sync trap (+ the side-map avoidance pattern).

**PHASES A.1 / A.2 / A.3 COMPLETE — the force-element + driven-joint tier is done (2D+3D).**

**Phase C DONE 2026-06-16 (Tao 5-DOF spindle, UNCOMMITTED).** Built the 6-joint stack
(3 prismatic x,y,z + 2 revolute θ,φ vibration DOFs on massless frames + a revolute SPIN
joint carrying the rotor body with cm offset `e` from the spin axis). Bearings modelled with
a NEW force element (see C.1) on the NON-spinning housing (the φ frame): two radial springs
at axial `±Lb` give the full radial stiffness AND — via the lever arm — the emergent tilt
stiffness; one axial spring at the cm gives the axial stiffness. Spin is a **driven**
(motor-clamped, Ω̈=0) joint. Tests in `ga_appl3dp_physics_test.hpp`
(`TEST_SUITE("PGA3DP: Tao wheel-spindle (Phase C)")`); 2D parity of the C.1 primitive in
`ga_appl2dp_physics_test.hpp`. Suites green: appl2dp 50/539, appl3dp 42/9168, pga 169/2758.
**No ga_py change** (the new element + the library fix are internal to `dynamic_system`,
which is unbound).

- **C.1 — grounded spatial spring/damper element** (the headline capability). New
  `grounded_spring{2,3}dp` + `add_grounded_spring(idx, anchor_b, [p0_world,] k, c)` on
  `dynamic_system{2,3}dp`: a body-fixed point tied to a ground anchor by anisotropic
  (world-axis) stiffness `k` + isotropic damping `c`. CONFIGURATION-dependent — its restoring
  wrench `wdg(P, F)` is recomputed from the live pose/velocity each RK4 sub-step (folded into
  `assemble_mass_bias`, NOT a function of time like `set_applied_wrench`); contributes its
  potential to `potential_energy()`. Because the wrench acts at the physical point, ONE spring
  yields both a translational stiffness and (via its lever arm) a tilt stiffness `k·l²`.
  Gate: emergent translational `k` + emergent tilt `k·l²`, 3D+2D, exact.
- **C.2 — characteristic frequencies.** f_x=2387, f_z=5485 (FIRST-PRINCIPLES, table
  stiffnesses, Eqs 16/17), f_θ0=6702 (calibrated, see deviation #1). All emerge from the
  assembled dynamics.
- **C.3 — centrifugal + gyroscopic.** Centrifugal `ẍ = −e·ω² = −0.19739 m/s²` EXACT (offset
  cm). Gyroscopic coupling `φ̈/θ̇ = (J_z/Jₓ)·ω` EXACT (textbook polar-inertia rotordynamics).
- **C.4 — z_b = z − R_w·φ** (Eq.14), the surface-error driver: the GA rim-point axial position
  equals `z − R_w·φ` exactly.

**LIBRARY FIX (2026-06-16, part of Phase C): driven-joint inertia as a moving base.**
`assemble_mass_bias` summed inertia only over the dof joints, so a body on a *driven* joint
had its inertia SILENTLY DROPPED (→ singular mass matrix when the driven joint carried the
only inertia below a dof joint — exactly the spindle's clamped spin). Fixed in both
`ga_pga{2,3}dp_ops_physics.hpp`: the body list is now `dof joints ∪ driven joints`, each
projected onto its ANCESTOR dof joints. Purely ADDITIVE (a driven *base* with no dof joints
above it — e.g. the A.3 spinning slider — is unaffected); all suites stayed green. This is a
genuine latent-bug fix, INDEPENDENT of the modelling deviations below.

**THREE DEVIATIONS from Tao's lumped Eq.13 — RESOLVED 2026-06-17, full derivation in
[tao_eq13_derivation.md](tao_eq13_derivation.md). Both are PAPER ERRORS, not coordinate
artifacts nor GA-model gaps; the GA model gives the physically correct coefficients.
Summary: (1) tilt stiffness is 4× too large = 2 [Tao adds k_y, but a tilt about y is
resisted only by the x-springs] × 2 [k_x is the TOTAL of two parallel bearings (pinned by
f_x), so per-bearing is k_x/2] — so Lb≈0.231 m just absorbs the 4×, and the real f_θ is
~34% lower. (2) Gyroscopic (J_y−J_z)ω is Tao dropping the kinematic frame-rotation term
+J_t ω φ̇ from d/dt(J ω_c); the correct coupling is the polar J_z ω (matches the C.3 gate
φ̈/θ̇=(J_z/J_x)ω and standard rotordynamics). Original flagged text retained below for
context:**

1. **Tilt stiffness** — Tao's `(kₓ+k_y)(l1²+l2²)` is ~4× a consistent two-spring lever-arm
   model (where translation must equal kₓ). To match Fig.4's f_θ we CALIBRATE the radial-
   bearing offset to `Lb = sqrt(K_tilt/kₓ) ≈ 0.231 m` (≈2× the table's l1=0.1). Confirmed NOT
   caused by units (strict SI) or by `e` (2 µm; the ±e in l1,l2 shifts l1²+l2² by <1e-7).
2. **Gyroscopic coupling** — the GA model gives the textbook polar-inertia coupling `J_z·ω`
   (forward/backward whirl split `(J_z/Jₓ)ω ≈ 1.76 Hz`); Tao's Eq.13 uses `(J_y−J_z)ω ≈ Jₓ·ω`
   (the TRANSVERSE inertia → his ±f_s/2 ≈ 25 Hz). The GA result matches standard
   rotordynamics; Tao's does not.
3. Both #1 and #2 point at the same theme: Tao's hand-lumped SCALAR Eq.13 coefficients do not
   cleanly reduce from consistent rigid-body mechanics — possibly artifacts of his
   inertial-frame tilt-angle definitions (Eqs 1/2), possibly errors.

**Eq.13 ANALYSIS DONE 2026-06-17** ([tao_eq13_derivation.md](tao_eq13_derivation.md)) —
both deviations are paper errors, the GA model is correct, and `Lb` is now documented as
absorbing Tao's 4× over-count (not a physical offset). Code comments updated to point at
the derivation (no numeric code change needed). The doc also records the **falsifiable
predictions**: corrections are tilt-only, so Figs. 4(c,d)/6(b)/15 and Eqs. 18/19/22 change
(single tilt peak ~3.7–4.4 kHz, no spin-dependent ±f_s/2 split — which matches Tao's own
"peaks same at 2250/2650 r/min" §4.2.1), while the x/y/z dynamics and the entire
experimentally-validated topography half (§3/§4, Table 3) are unchanged (they use the
*measured* f_b = 6253 Hz, not Eq.13's predicted frequency).

**Phase D.1 DONE 2026-06-17 (UNCOMMITTED): grinding-mark wavelengths.** App-test
`"pga3dp: grinding-mark wavelengths lambda_m / lambda_c (Phase D.1)"` in
`ga_appl3dp_appl_test.hpp` (right after the Phase-0 grain-trajectory case). Reproduces Tao's
Table-3 SIMULATED wavelengths from the GA twist field: `λ = (relative surface speed)/f_b`,
where the wheel surface speed `v_w=n_w·R_w` (→ WMD `λ_m`) and the wafer surface speed
`v_s=n_s·r` (→ WCD `λ_c`) are read straight off `kinematic_system3dp::point_velocity` of the
spinning frames (same twist field as the disc-surface-speeds test; same Fig-1 tree as
Phase 0). f_b = measured 6253 Hz (Tao §4.2.2), i.e. the surface half consumes the MEASURED
tilt freq, independent of the Eq.13 controversy. Gate (all <1 µm vs Table 3 sim): λ_m=5652
(n_w=2250), 6657 (2650); λ_c=133 (n_s=265), 163 (325). + Fig-12 trends (λ_c∝r, λ_m indep of
n_s). Full appl3dp suite green: 43 cases / 9174 assertions. Pure kinematics, NO library
change, NO ga_py impact.

This is the LOW-RISK feed-forward slice of the contact/force loop (no feedback).

**Phase D.1b DONE 2026-06-17 (UNCOMMITTED): surface-generation model.** App-test
`"pga3dp: simulated wafer topography profiles (Phase D.1b)"` (right after D.1). Implements
Tao's grain-cutting model Eqs. 23-28 in 1D: Gaussian grain sizes (Eqs 23/24), spacing
`l_g = μ·(π/6N_g)^⅓` (Eq.25), each grain centred at arc length `s_k` cuts to axial depth
`z_b(t_k)=A_b·sin(2πf_b t_k)` (Eq.26 G^z; feed ramp negligible) and leaves a spherical-cap
furrow `0.5d−√(0.25d²−δ²)` (Eq.27); surface = lower envelope over all grains (Eq.28, min).
The waviness EMERGES from the carve (not the formula) and its measured period (box-smooth +
zero-crossings) matches the D.1 prediction `v/f_b`: WMD λ_m carved 5720 vs 5652 µm, WCD λ_c
carved 134.6 vs 133.1 µm (~1%, gate 5%). Sparklines in the log show ~5 clean sine periods.
Pure kinematics, NO library/ga_py change; appl3dp 44/9176 green. (Full 2D raster of Figs
9/10/16/17 with circumferential jaggedness = a ga_view scene, D.1c, deferred — visualization,
no new gate.)

**Phase D.1c DONE 2026-06-17 (COMMITTED 5639a8f + 879ac0d): two ga_view scenes.** Both
build the rendered geometry from `kinematic_system3dp::point_velocity` (the GA twist field),
visually confirmed by the user.

- **`active_grinding_topo`** (5639a8f) — Tao Figs. 9/10 ground-surface topography as a jet
  heatmap; C toggles WCD (square) / WMD (near-square box, height exaggerated). Shows live
  λ_m/λ_c + a height colorbar + quantified patch axes + an animated feed-scroll readout. The
  7.6° grain-mark inclination is GA-computed at the true contact point with the co-rotating
  (down-grinding) wheel — confirmed to match Tao's 7.6° (counter-rotating gives ~11.5°).
- **Fig. 12 wavelength trends** — Tao Fig. 12 as **four separate NATIVE ga_view line
  charts** (`get_grinding_wavelength_panels()` in `w_mainwindow.cpp`), one static
  `Coordsys_model` per panel (a/b/c/d), each with its own paper axis range/labels. Curves
  built inline by the same `kinematic_system3dp` `surface_speeds` code as the D.1 app-test
  (`λ = point_velocity / f_b`); families drawn with `add_ln` polylines + per-family colour
  & marker (□/○/+/×) + a `diagram_legend`. Reproduces the paper AND the GA-derived
  independences (λ_c indep of n_w; λ_m indep of r and of n_s). **PAPER ERROR documented in
  code:** Fig. 12 caption says f_b=4000 for panels (b)/(d), but the plotted curves are at
  f_b=6000 (clean 1/f_b factor; user confirmed (d) tops ~7200 µm) — we use 6000 to match.
  **REWORKED 2026-06-18:** the original `active_grinding_wavelength` (879ac0d) hand-painted
  a custom 4-panel chart `QGraphicsItem` (axes/ticks/legend all hand-drawn, ~355 lines) —
  overkill for a GA viewer. Replaced by the four native charts above, enabled by a small
  reusable ga_view core feature: a per-model axis override (`Coordsys_model::axis_cfg` +
  `plot_axis_cfg`; `Coordsys::set_axes`; `w_Coordsys::switch_to_model` retargets the shared
  aspect-locked axes per chart model and restores the geometric axes on exit, re-emitting
  `viewResized` so ctor-anchored marker items reposition). Old scene + its 5 wiring points
  removed.

**PHASE D.1 COMPLETE** (D.1a wavelengths + D.1b surface model + D.1c two ga_view scenes; all
driven from the GA twist field, all committed). The whole surface-formation half of Tao is
now reproduced and validated against the paper from the GA kinematics, and shown independent
of the Eq.13 controversy.

**RESUME HERE (next session) → D.2 (the risky/open research part): the two-way force loop** —
grinding-force law → wrench → vibration → z_b → depth → MRR → force. Tao has NO force model
(assumes F=25/25/100 N const), so D.2 is genuinely uncertain modeling + may stiffen the ODE
(RK4 dt≈1e-8 already). **Detailed, clarified plan now in the dedicated section
[Phase D.2 detailed plan](#phase-d2-detailed-plan-two-way-force-loop) below** (user
clarifications 2026-06-20: build flexible infrastructure for a student to
experiment with force laws — prepare BOTH the depth-based and MRR-based laws, BOTH loop
topologies, BOTH integrators; control strategy as a separate opt-in top layer). Then E —
feed/speed optimization. (Off critical path: **B.2** —
CLOSED 2026-06-18, subsumed by C.3 (the emergent centrifugal `m·e·Ω²` from a spinning
offset-cm body was validated EXACTLY there, `ẍ=−e·ω²`; no separate B.2 test needed);
**Phase 0.c** CS-view scene — still deferred.)

**RESOLVED (2026-06-14): static-vs-feed = PRESCRIBED INFEED.** The spindle carries a
steady prescribed translation (axial infeed `x_a` along −e1 and/or radial traverse) ON TOP
of the 5 vibration DOFs — so `x_a` (and possibly the radial offset) is a time-varying
prescribed input, not a constant. Phase C must drive the nominal pose with this prescribed
motion; the 5 vibration DOFs perturb around it. Phase 0 / B.1 still run at the static
nominal pose (the warm-up validates the spring/damper primitive, independent of feed).
**Qt "Sans Serif" font warning — FIXED 2026-06-18:** `ga_view/src/main.cpp` now sets a
default app font (`QApplication::setFont(QFont("Helvetica", 10))`); the startup warning is
gone.

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
  - **B.2 [CLOSED 2026-06-18 — subsumed by C.3]:** the goal was to replace the explicit
    forcing with an **emergent** spinning offset-cm body and reproduce B.1 (validating
    "`meω²` for free" + the prescribed-spin bias path). Phase C.3 already does exactly this
    in the full 5-DOF spindle: the centrifugal acceleration `ẍ = −e·ω²` emerges EXACTLY
    from the spinning offset-cm rotor body (driven-joint moving base), so a separate B.2
    prototype adds nothing. Closed without new code.
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

## Phase D.2 detailed plan (two-way force loop)

Clarified with the user 2026-06-20. **Design principle:** build *flexible infrastructure a
student can experiment on*, not one hard-wired model. Keep the experiment-specific pieces
(the **force law**, the **control schedule**) cleanly separated from the infrastructure, and
add each new capability as an **opt-in layer that is only paid for when used** — the same
side-map / wrapper pattern as the driven joints, applied wrenches and the closed-loop tier.
Every choice below is prepared for BOTH alternatives the user named, with a stated default.

### Decisions captured (2026-06-20)

- **Force law — prepare BOTH, default (a).** (a) engagement-depth / chip-thickness law
  [DEFAULT]; (b) MRR × specific-energy law. Interchangeable behind one registration hook; a
  student swaps laws without touching the infrastructure.
- **Loop topology — prepare BOTH, default spindle-only.** (1) spindle-side only [start];
  (2) + chuck/wafer reaction, needed because the wafer thinning makes engagement (hence
  force) sensitive. The wafer-side reaction + the wafer-thinning state are opt-in.
- **Integrator — BOTH, RK4 default.** RK4 [default] + an Adams–Bashforth multistep with
  adjusted (adaptive) stepping for speed. A benchmark test reports RK4-vs-AAB accuracy
  (gate: ≤1 % delta) and wall-clock. (AAB read as an Adams–Bashforth predictor / Adams–
  Moulton corrector pair — "alternating" = the predictor-corrector alternation — with
  step-size control from the P-C error estimate; CONFIRM exact scheme/order with the user.)
- **Home — app-test prototype first**, then promote the proven primitive into the library
  (standing project decision). Validation ladder: feed-forward force → close loop → thinning
  → control.
- **Control strategy — a SEPARATE top layer**, only used/paid for when attached (mirrors
  `closed_loop_system` on top of `dynamic_system`). Options: (a) hold force ≤ 70 % of the
  constant-feed reference by controlling feed; (b) constant feed until wafer = 150 µm, then
  taper the force target linearly 100 %→70 % between 150 µm and 80 µm. Both are starting
  values to be tuned later.

### Infrastructure to add (mirrors the existing force-element side-maps)

1. **`contact_state` (POD)** — everything any force law could read, computed by the infra
   from the live kinematics each RK4 sub-step: engagement depth `delta` (rim axial position
   `z_b = z − R_w·φ` vs the nominal wafer-top, so vibration modulates it), relative surface
   speed `v_rel` (from `kinematic_system::point_velocity`, the Phase-0/D.1 twist field),
   nominal feed velocity `v_feed`, current wafer thickness `tw`, contact point `P` (world)
   and feed-normal `n`, swept area / MRR, clock `t`.
2. **Force-law hook** — `using grinding_law = std::function<vec3dp(contact_state const&)>;`
   registered per spindle frame: `set_contact_force(idx, law, contact_params)`. The infra
   builds the wrench `wdg(P, F)` and folds it into `assemble_mass_bias` exactly like the
   grounded spring (config-dependent, NOT a function of time). Two laws shipped as free
   functions a student picks or replaces: `grinding_force_depth` (a, default) and
   `grinding_force_mrr` (b). `contact_params` = geometry/material constants (R_w, nominal
   wafer-top, contact width/arc, specific energy `k_e`, …).
3. **Wafer-thinning state + chuck reaction (opt-in)** — `wafer_thickness_` integrated from
   removed volume (`tw_dot = −MRR / area`); feeds back into `delta`.
   `set_contact_reaction(chuck_idx)` enables the equal-and-opposite wrench `−wdg(P,F)` on the
   wafer/chuck frame (topology 2). Absent ⇒ topology 1, gravity/bias path byte-unchanged.
4. **Integrator selector** — `set_integrator(rk4 | adams_bashforth)`, RK4 default. AAB stepper
   added to `ga_usr_utilities.hpp` beside `rk4_step` (history buffer, RK4 self-start, adaptive
   dt from the P-C error estimate). `coupled_step` dispatches on the selection.
5. **Feed-control layer (separate, opt-in)** — a `grinding_controller` holding the schedule
   (option a/b params) + `feed_rate(measured_force, tw) -> rate`; wired by the test into
   `set_driven_rate` each macro-step (or a `set_feed_controller` side hook). Absent ⇒
   constant feed; the layer is inert and unpriced when not attached.

### 3D-only by design (no 2D twin) — decided 2026-06-20

The contact / grinding-force element (D.2b/c) is added to `ga_pga3dp_ops_physics.hpp` ONLY,
breaking the otherwise-strict 2D+3D force-element parity **deliberately**. The earlier force
elements (spring/damper, applied wrench, driven joint, grounded spring) are general mechanics
primitives with real 2D consumers (2D closed-loop scenes, the merry-go-round), so parity
earned its keep. This element is the first **application-specific** one: wafer
self-rotational grinding is inherently 3D (`z_b = z − R_w·φ`, tilt, radial) and its vocabulary
(`wafer_thickness`, `removal_rate`, feed) is grinding-flavoured with no 2D consumer — so a 2D
mirror would be parity-for-its-own-sake against the "only paid for if used" principle. The
general kernel buried in it (a configuration-dependent applied wrench from a user callback on
the live state — `set_applied_wrench` generalised from time-only to state-dependent) can be
factored into 2D IF a genuine 2D contact/penalty use ever arises. Recorded in a header
comment at `contact_state3dp` so a future session does not "fix" the asymmetry by reflex.

### Sub-phases (the validation ladder)

- **D.2a — contact kinematics + force laws, FEED-FORWARD (no back-reaction). DONE
  2026-06-20 (UNCOMMITTED).** App-test `"pga3dp: grinding contact force, feed-forward (Phase
  D.2a)"` in `ga_appl3dp_appl_test.hpp` (right after D.1b). Prototypes the `contact_state`
  POD + a swappable depth force law (`F_normal = k·δ`, tangential `μ·F_normal` opposing the
  GA-derived sliding direction) locally in the test — to be promoted into `dynamic_system` in
  D.2b. The relative sliding velocity is `point_velocity(wheel grain) − point_velocity(wafer
  point)` (same Phase-0/D.1 twist field); the wafer normal is the chuck axis e1 and the wheel
  spins about −e1, so the tangential force carries the in-plane (F_x,F_y) and the axial F_z
  comes purely from the depth term. **Calibration (the documented anchor):** with nominal
  engagement `δ0 = 1 µm`, `k = F_z_ref/δ0 = 1.0e8 N/m` reproduces Tao's axial `F_z = 100 N`,
  and `μ = √(25²+25²)/100 = 0.3536` matches the in-plane reference `35.36 N`. **Gates (6
  assertions):** axial == 100 N + in-plane == 35.36 N (calibration round-trip); sliding ⟂ e1;
  rim grain GA surface speed == `n_w·R_w` (47123.89 mm/s exact); and the feed-forward `z_b =
  A_b·sin(2πf_b t)` modulation ripples F_z over [90, 110] N == `k·(δ0∓A_b)` (clean 3-period
  waveform in the log). Pure kinematics, NO library/ga_py change; appl3dp 45/9182 green.
- **D.2b — close the loop, spindle side (topology 1). DONE 2026-06-20 (UNCOMMITTED).**
  PROMOTED the contact element into the LIBRARY `dynamic_system3dp`
  (`ga/ga_pga3dp_ops_physics.hpp`): `contact_state3dp` POD + `grinding_law3dp` swappable
  callback + free law `grinding_force_depth(c,k,mu)` (default) + `set_contact_force(idx,
  wafer_idx, law, n_hat, delta0, v_feed, tw, contact_b)` / `clear_/has_contact_force` /
  diagnostics `contact_force(idx)` `contact_engagement(idx)`. A private `eval_contact_state`
  helper builds the live state (δ = δ0 − (P·n̂ − z_ref); v_rel from the twist field) shared by
  the `assemble_mass_bias` fold-in (same `wdg(P,F)` projection as the grounded spring,
  config-dependent NOT time) and the diagnostic. New suite
  `TEST_SUITE("PGA3DP: grinding force loop (Phase D.2b)")` in `ga_appl3dp_physics_test.hpp`:
  a SIMPLIFIED rig (axial prismatic DOF + bearing spring/damper k_z + a driven massless
  flywheel for the real sliding speed; full Tao 5-DOF spindle = next refinement). The loop
  `m z̈ + c ż + k_spring z = k_grind(δ0 − z)` has the analytic equilibrium
  `z_eq = k_grind·δ0/(k_spring+k_grind)`, `F_eq = k_grind(δ0 − z_eq)`. **Gates (8 assertions):**
  M[0]==m + initial accel == F_ff/m; rim sliding speed == n_w·R_w; **closed loop settles to
  z_eq/F_eq exactly** (RK4 dt=1e-6, 20 ms); and the HEADLINE — spindle compliance pulls the
  force BELOW the feed-forward nominal (100 N → **90.476 N, −9.5%**), with the analytic limit
  `F_eq = F_ff·k_spring/(k_spring+k_grind)` → F_ff as k_spring→∞ (recovers D.2a). Library
  change purely ADDITIVE: pga core 169/2758 UNCHANGED, appl3dp 46/9190; no ga_py impact
  (dynamic_system unbound). 2D twin of the contact element = deferred decision (grinding is
  inherently 3D; mirror only if a 2D analog is wanted).
- **D.2c — wafer thinning + chuck reaction (topology 2). DONE 2026-06-20 (UNCOMMITTED).**
  Library (`ga/ga_pga3dp_ops_physics.hpp`): `contact_spec` += `k_mrr` (removal const) + `react`
  (toggle); `set_contact_force` += optional `k_mrr` arg; new `set_contact_reaction(idx,on)`,
  `update_wafer_thinning(dt_macro)`, `wafer_thickness(idx)`, `removal_rate(idx)`; reaction
  fold-in applies `−wdg(P,F)` on the wafer frame's joints when `react`. Two cases in new
  `TEST_SUITE("PGA3DP: grinding wafer thinning + reaction (Phase D.2c)")`: (1) thinning at
  `k_mrr=0.01` thins 0.4264 µm/step monotonically, decrement == `MRR·dt_macro` exact
  (`MRR=k_mrr·δ_eq·v_rel=4.26e-7 m/s` ~ Tao feed scale); (2) third-law — two prismatic DOFs
  (m1=0.8, m2=2.0), reaction ON → accel ratio −0.400 = −m1/m2, OFF → wafer accel 0. Additive:
  pga core 169/2758 unchanged, appl3dp 48/9208; no ga_py impact. — original sub-plan: (1)
  **Wafer thinning** — a
  quasi-static thickness decremented between macro-steps by a SIMPLE removal law `MRR =
  k_mrr·δ·v_rel` with ONE tunable constant `k_mrr` (lumps contact width / removal
  efficiency; documented, like the force calibration). The live `tw` is exposed in
  `contact_state` so a force law / the D.2e control can read it. (2) **Chuck-side reaction**
  — an OPTIONAL user-configurable toggle `set_contact_reaction(idx, on)`; **default = rigid
  wafer (off)**. When on, the equal/opposite wrench `−wdg(P,F)` is folded onto the wafer
  frame's supporting joints. The realistic chuck-side reaction physics is a LATER refinement;
  D.2c just builds + tests the infrastructure for a simple case. Gate: `tw` decreases
  monotonically at the analytic rate `k_mrr·δ·v_rel`; with the reaction ON, a compliant
  wafer DOF accelerates opposite the spindle (Newton's third law, accel ratio = −m_spindle/
  m_wafer) and with it OFF the wafer stays put.
- **D.2d — Adams–Bashforth integrator + adaptive stepping + comparison.**
  **D.2d-1 (constant dt) DONE 2026-06-20 (UNCOMMITTED).** Two integrators added to
  `ga/ga_usr_utilities.hpp` behind a uniform `step(rhs, u, t, dt)` callable interface:
  `rk4_integrator` (WRAPS the canonical mdspan `rk4_step` — not a 2nd RK4) and
  `abm2_integrator` (Adams–Bashforth-2 predictor / Adams–Moulton-2 corrector, PECE, RK4
  self-start; 2 rhs evals/step vs RK4's 4). NEW dedicated test target **`ga_integrator_test`**
  (6th doctest executable; NORMAL build so it RESPECTS `CMAKE_BUILD_TYPE` — unlike the
  forced-O3 benches in `utilities/` — so its timing is meaningful in Debug AND Release, and
  it rebuilds/reruns on Windows). Damped-oscillator analytic cases (3): accuracy vs closed
  form (RK4 6.1e-10, ABM2 6.1e-5; agree <1 %), convergence order (RK4 16.02 = 4th, ABM2 4.01
  = 2nd), and a timing case printing the build mode. **KEY DEBUG-vs-RELEASE FINDING:** Debug
  (-O0) RK4 813 ms / ABM2 39 ms (ratio 0.05, RK4 ~20× slower — mdspan indexing not inlined,
  so the comparison is MEANINGLESS in Debug); Release (-O3) RK4 49 ms / ABM2 21 ms (ratio
  0.43 = the true 2-vs-4-eval ratio). Documented in a comment in the timing case. All green
  (8/8; ega/pga unchanged). **User directive honoured: simple analytic cases UPFRONT, separate
  test file, BEFORE the grinding loop / full Tao.**
  - **D.2d-2 DONE 2026-06-20 (UNCOMMITTED): integrator selector on the grinding loop +
    stiffness verdict.** `dynamic_system3dp::set_integrator(integrator_kind::rk4 | abm2)` +
    `get_integrator()`; `coupled_step` refactored to a SINGLE shared forward-dynamics rhs
    dispatched to `rk4_integrator` (default) or a persistent `abm2_integrator` member (the
    multistep history; reset on switch). **RK4 stays BYTE-IDENTICAL** — full physics suite
    unchanged (appl2dp 50/539, pga 169/2758, appl3dp pre-D.2d counts). New suite
    `TEST_SUITE("PGA3DP: grinding loop integrator (Phase D.2d-2)")`: both integrators reach the
    analytic equilibrium F_eq=90.4762 N and **agree to 1.57e-16** (the steady state is a fixed
    point, order-independent). **STIFFNESS VERDICT (the purpose of this step): NOT classically
    stiff** at k_grind=1e8 — w_eff≈36229 rad/s, working dt=1e-6 is **28× below the ABM2
    stability limit (2.76e-5)**, so dt is ACCURACY-bound (resolving the ~6 kHz vibration), not
    stability-bound. → an implicit solver (D.2d-4) is NOT needed now; only if a near-rigid
    contact penalty injects a fast mode. appl3dp 49/9212.
  - **D.2d-3 (add-on) DONE 2026-06-20 (UNCOMMITTED): variable/adaptive dt.** New
    `abm2_adaptive_integrator` in `ga/ga_usr_utilities.hpp`: the predictor-corrector
    difference is the Milne local-error estimate `err≈(1/6)|u_c−u_p|`; a standard controller
    `h_new=h·clamp(0.9·err^(−1/3),…)` sets the step (reject + retry on `err>1`). Non-uniform
    steps use the VARIABLE-step AB2 predictor `u_p=u_n+h[(1+r/2)f_n−(r/2)f_{n−1}]`,
    `r=h_n/h_{n−1}` (the trapezoidal corrector is ratio-invariant); RK4 self-start. Two
    analytic cases in `ga_integrator_test` (new suite): (1) error control on the damped
    oscillator — tighter rtol → smaller error + more steps (1e-4: err 3.3e-3, 1579 steps;
    1e-7: err 3.1e-5, 16443 steps); (2) step ADAPTATION on `y=e^{−t²}` under absolute control
    — dt shrinks where `|y'''|` peaks (steep middle 2.84e-3) and grows at the flat tail
    (4.21e-3), range 14.8× (`dtmax/dtmin`), err 1.3e-6. **(3) SPEEDUP vs fixed dt (step count
    at MATCHED global accuracy; both ~2 evals/step, so step ratio = cost ratio): on `e^{−t²}`
    over [0,6] adaptive 689 vs fixed 1984 → ~2.9× fewer (the flat tail cost adaptive +44 steps
    vs +992 for fixed; grows with a longer flat region). On a UNIFORM-timescale problem (the
    oscillator) there is nothing to exploit → adaptive ≈ fixed, the honest caveat.** 15/15;
    physics suites unchanged. **NOTE: adaptive dt is NOT a stiffness remedy** — for an EXPLICIT method the step on a stiff system
    is STABILITY-bound, so the controller just auto-shrinks `dt` to that limit and crawls; it
    buys error control / auto-tuning + efficiency on VARYING dynamics, not stiffness. (Not
    wired into `dynamic_system` — the grinding loop is accuracy-bound & non-stiff per D.2d-2,
    so a variable-dt `step()` there adds nothing now; wire it only if a use needs it.)
  - **D.2d-4 (the actual stiffness lever, IF D.2d-2 shows stiffness): an IMPLICIT solver** —
    A/L-stable: implicit trapezoidal / Adams–Moulton solved with a Newton iteration, or a BDF.
    Only these take large steps on a genuinely stiff system (limited by accuracy, not
    stability). Likely triggered by a near-rigid CONTACT penalty (large `k_grind`) injecting a
    mode far above the ~6 kHz physics band. Its own analytic stiff-ODE test case upfront.
- **D.2e — control layer (force limiting). DONE 2026-06-20 (UNCOMMITTED).** New
  `grinding_controller3dp` (library, `ga/ga_pga3dp_ops_physics.hpp`): an OPT-IN, self-contained
  process-level controller — `target_fraction(tw)` gives the force target ((a) constant `cap`
  when `tw_hi<=tw_lo`; (b) linear taper `1→cap` between `tw_hi` and `tw_lo`), `feed_command(tw)
  = v_feed_nominal·target_fraction` (feed-forward; F∝feed in the removal balance). Nothing in
  `dynamic_system` depends on it. New suite `TEST_SUITE("PGA3DP: grinding feed control (Phase
  D.2e)")` runs the quasi-static PROCESS loop (`v_feed=MRR` balance ⇒ steady `F=gain·feed`,
  `tw -= feed·dt`). **THE KEY RESULT (answers "is 70 % force just feed→70 %?" — NO):** baseline
  100 N/240 s; option (a) 70 N but **343 s = 1.43× longer** (the throughput cost: removal also
  drops to 70 %); option (b) taper keeps F=100 N above 150 µm, 267 s = 1.11× (cheaper, only
  throttles the thin-wafer phase). Gates: F_a==70 N, t_a/t_base==1/0.7; taper schedule
  (full≥150 µm, 0.85 at 115 µm, 0.7≤80 µm). Documented secondary effects (controller header):
  feed↔force not a static gain (compliance/transients); waviness `z_b`/`λ=v/f_b` UNAFFECTED by
  feed (benefit = damage avoidance, not finish); thinner wafer more fragile (option b
  rationale); F∝feed linearity is law-specific. Additive: pga 169/2758 unchanged, appl3dp
  50/9221. No ga_py impact (controller is a plain struct, dynamic_system unbound).

Then **E — optimization** on the calibrated closed-loop model (feed/speed ratio to bound
force and `z_b`).

### Open items — RESOLVED 2026-06-20

1. **AAB scheme = AB2 predictor + AM2 (trapezoidal) corrector**, a PECE predictor-corrector
   (explicit Adams–Bashforth alternated with the implicit Adams–Moulton step), with adaptive
   dt from the P-C error estimate. Start with this simplest 2nd-order form — "something
   working that can be extended" to higher order later.
2. **Force law: simplest first, chip-thickness as the first refinement.** D.2a ships the
   linear depth law `F_normal = k·δ`, tangential `F_tan = μ·F_normal` along `−v_rel` — get
   the loop running at all. The chip-thickness form (force per active grain × engagement
   area) is the first refinement once the simple law closes the loop.
3. **Calibration anchor [MUST DOCUMENT — do not forget]:** tune the depth-law constant `k`
   so the *constant-feed* nominal force reproduces Tao's reference `F = (F_x,F_y,F_z) =
   (25, 25, 100) N` (Tao §2.2, from Zhu [70]). This fixes a concrete number for D.2e's
   "≤ 70 % of constant-feed force" target and ties the model back to the paper. Record the
   calibrated `k` (and the contact geometry it assumes) in the app-test header comment AND
   here when D.2a lands.
4. **Wafer thinning = quasi-static, updated slowly between MACRO-steps** (not integrated
   in-loop). The thinning timescale (µm/min) is ~10⁹× slower than the kHz vibration, so the
   thickness is held constant within a macro-step and decremented by `MRR·Δt_macro / area`
   between steps.

## Parked markers (future capabilities)

- **`get_disc_inertia` (round geometry) — DONE 2026-06-20.** Added in 2D and 3D, mirroring
  `get_plate_inertia`/`make_plate_body` (2D) and `get_cuboid_inertia`/`make_cuboid_body` (3D)
  incl. the optional parallel-axis (Steiner) pivot. **2D** `get_disc_inertia(m, r)` +
  `make_disc_body(m, r)`: polar moment `J = ½ m r²`. **3D** `get_disc_inertia(m, r, t)` +
  `make_disc_body(m, r, t)`: solid cylinder, SYMMETRY/spin axis along **e3** (disc in the
  e1-e2 plane, thickness t along e3) — axial `I_zz = ½ m r²`, transverse
  `I_xx = I_yy = m(r²/4 + t²/12)`. Unit tests (appl2dp/appl3dp "physics tests prep"): moments,
  `I·I_inv = identity`, Steiner (rim pivot adds `m r²`), the thin-disc perpendicular-axis
  theorem `I_zz = I_xx + I_yy`. No ga_py impact (the sibling inertia helpers aren't bound).
  Ready for the wheel/wafer bodies when their true inertia matters (replaces the cuboid
  approximation in the rigs).

## Documentation plan (grinding-specific subsection) — DONE 2026-06-20

**DONE:** new file [ga_docu/6_ga_applications_pga.tex](ga_docu/6_ga_applications_pga.tex) —
`\subsection{Applications of PGA}` → `\subsubsection{Wafer Grinding}`, wired into
`0_ga_docu.tex` (after the physics chapter) + Tao bib entry (`Tao_wafer_grinding:2022`).
Compiles clean (PDF §1.5.1). User-approved ("good for a student to get started with the
code"). Covers, concept/result-focused (NOT phase/step history): the frame tree + `λ=v/f_b`;
the additive opt-in force-element tier (spring/damper, applied wrench, driven moving base,
grounded spring `k·l²`); the grinding contact + closed loop (`F=k·δ`, `wdg(P,F)`,
`δ=δ0−z_b`, compliance equilibrium) with the swappable-law callback called out as the
experiment surface; the 5-DOF spindle + the two Eq.13 corrections (gyroscopic `J_z`, tilt
stiffness 4×) as findings + implications (tilt-only; topography half unaffected); the
RK4/ABM2 selector + stability-vs-accuracy (accuracy-bound, not stiff) + the integrator-limits
DOC TASK below FOLDED IN (incl. debug-vs-release → inline-control); and the feed controller +
the 70 %→1.43× throughput trade-off. Extensibility seams surfaced throughout; only
built features documented. The fuller Eq.13 algebra stays in
[tao_eq13_derivation.md](tao_eq13_derivation.md).

--- original note (for context) ---

What we have documented so far in `ga_docu` (§5 multi-body motion + the §8 glossary) is
**generic** rigid-body / multi-body dynamics — force/torque, momentum, inertia maps,
Newton–Euler, frame trees, articulated + closed-loop dynamics. The wafer-grinding work is now
getting **application-specific** (the contact/force element, the grinding law, wafer thinning,
the closed force loop, the Tao spindle, the integrator trade-off, the two Tao paper errors),
so it warrants its **own dedicated subsection** in the documentation — NOT folded into the
generic multi-body sections. Write it as the capstone case study (after D.2/E land, per the
agreed order: D.2 first, then docs), covering: the Fig.1 frame tree, the force-element +
contact-loop tier, the Eq.13 paper-error findings ([tao_eq13_derivation.md](tao_eq13_derivation.md)),
and the surface-formation half (D.1). Keep the generic-vs-grinding split explicit.

**DOC TASK — integrator stability & accuracy limits (worked on the grinding model), noted
2026-06-20.** Document, with the grinding model as the running example, the two distinct
upper bounds on the time step and the lower bound from compute:

- **Stability limit (explicit methods):** RK4 `dt ≲ 2.78/ω`, ABM2 (PECE) `dt ≲ ~1/ω`, with
  `ω = √(k_tot/m)` the fastest structural mode (grinding loop: `ω≈36229 rad/s` →
  `dt_RK4≈7.7e-5`, `dt_ABM2≈2.8e-5`). Derive the boundaries from each method's amplification
  factor on `ẏ=λy` (the imaginary-axis interval for an oscillatory mode).
- **Accuracy limit:** `dt` small enough to resolve the physics of interest — the `~6 kHz`
  vibration `f_b` (`dt ≲ 1/(20 f_b) ≈ 8e-6`). On the grinding loop this is the BINDING
  constraint (well below the stability limits) — i.e. the loop is accuracy-bound, not stiff
  (the D.2d-2 verdict). Show the convergence-order evidence (RK4 4th, ABM2 2nd) from
  `ga_integrator_test`.
- **Compute lower bound (the real machine question):** the durable end-point — for INLINE /
  real-time control on the actual machine, how small can `dt` be? Set by the per-step compute
  cost (ns/step, Debug vs Release, from `ga_integrator_test`: e.g. Release ~21–49 ns/step for
  this 2-DOF rig) versus the control-loop period and the model size. The feasibility window is
  `[dt_min(compute, real-time budget), dt_max(stability ∧ accuracy)]`; document when ABM2's
  ~2× fewer evals/step buys real-time headroom and when a smaller/implicit step is forced.
  This ties the Debug-vs-Release timing finding directly to a deployment constraint.

## Deferred documentation tasks — DONE 2026-06-18

Glossary entries added to [ga_docu/8_ga_glossary.tex](ga_docu/8_ga_glossary.tex) (compiles
clean, all cross-refs resolve). Eleven entries, in the file's existing register (terse,
declarative, equation-first), each with `\ref{}` back-pointers into the full §3/§5
treatments:

- The six originally deferred: **pose**, **pose vs. twist**, **twist**, **wrench**,
  **Lie group vs. Lie algebra**, **infinitesimal vs. finite movements**.
- Group A (connective tissue the above reference): **motor and rotor**, **screw (Chasles'
  theorem)**, **inertia map**.
- Group B (Lie-cluster stubs): **adjoint**, **commutator (Lie bracket)**.

Skipped by decision (Group C — §2-basics territory, defined there already): bulk/weight,
unitization, congruence, attitude, forque.
