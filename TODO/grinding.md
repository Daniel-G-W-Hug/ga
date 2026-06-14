# Wafer-grinding dynamics — modelling & optimization plan

Status: **planning + initial port done (2026-06-14).** Long-horizon project; this file
evolves and is the living plan across sessions. The kinematics/dynamics infrastructure
(`static_/kinematic_/dynamic_system{2,3}dp`, the closed-loop layer) was built as the
runway for exactly this. Companion memory: `project_wafer_grinding`.

**DONE so far:** the two grinding mini-tests translated onto `static_/kinematic_system3dp`
(pedagogical, in `ga_appl3dp_appl_test.hpp`, right after the originals); the `is_congruent`
refinement (intermediate step below). **RESUME HERE → Phase B.1** (Sommerfeld force-element
prototype): see "First chunk (next): Phase B.1" at the end of Phase B, plus the design
already agreed in chat — explicit rotating unbalance wrench `meω²(cos,sin)` + x,y spring/
damper on a `(m+M)` body, prescribed spin, validated against the closed-form Eq. (4)
amplitudes over an ω sweep (resonances 16.51 / 28.57 rad/s, Sommerfeld Table-1 params).
Prototype in a `ga_appl3dp` app-test first; promote to the library only after it matches.

## Goal

Model the geometric **and** dynamic aspects of **wafer self-rotational grinding** to
understand the grinding process and — visionary goal — **optimize process parameters**
(feed rate of the grinding wheel relative to the chuck, speed ratio) to dynamically control
the grinding force and **prevent wafer damage**.

Target system (Tao et al., *Int. J. Mech. Sci.* 232 (2022) 107620): a cup grinding wheel on
an **aerostatic-bearing spindle** grinds a wafer held on a rotating **chuck**. The spindle
vibrates in radial (x,y), axial (z) and **tilting** (θ,φ) DOFs (5-DOF model, mass
eccentricity `e` → centrifugal `meω²`, gyroscopic tilt coupling). The axial runout
`z_b = z − R_w·φ` (Eq. 14) modulates grain cutting depth → surface waviness (WMD along the
mark direction, WCD circumferential; wavelengths λ_m, λ_c). A **real machine exists** with
**measured time- and frequency-domain vibration data** at several points, for calibration
and validation.

## Starting point (already in the repo)

Two prototype tests in `ga_test/src/ga_appl3dp_appl_test.hpp` (pre-infra, raw PGA):

- *"reference and tumbling plane"* — wheel **tilt kinematics**: a wobbling tool normal, the
  tool∩wafer intersection line and its angle. = Tao's (θ,φ); `z_b` is just a rim point's
  axial coordinate under a motor sandwich, no separate equation.
- *"intersecting discs"* — **contact/removal kinematics**: tool ring vs. wafer disc, overlap
  MRR volume, and the relative **surface speed** `rs·B_uv >> (P − cp)`. = Tao's Eq. (26)
  grain velocity.

The geometry/observation half is ~prototyped; the missing half is the **dynamics** driving
it. **Translated** versions of both (expressing the same idea via `static_/kinematic_system3dp`)
sit right after the originals — pedagogical, no new insight, asserted congruent to the raw
result.

## The one capability gap: FORCE ELEMENTS

`dynamic_system{2,3}dp` has only gravity + the Coriolis bias as generalized forces, and only
**rigid** joints. Both the target and the warm-up need **springs + dampers to ground,
rotating unbalance `meω²`, and an applied wrench** (grinding force; motor torque).

Clean mapping onto the existing infra: model the spindle as **one rigid rotor body connected
to ground by a serial stack of single-DOF joints** — 3 prismatic (x,y,z), 2 revolute (θ,φ),
1 driven revolute (spin Ω) — **inertia on the final body**, a linear `−k·q − c·q̇` force on
each mount joint, and eccentricity `e` = the rotor's **cm offset**. Then `assemble_mass_bias`
already yields `M(q)` + the gyroscopic/Coriolis bias; we only add the spring/damper and
applied-wrench terms to `tau` and integrate with RK4.

**Two payoffs for free:** the centrifugal `meω²` forcing emerges from a spinning offset-cm
body; the gyroscopic tilt coupling `(J_y−J_z)ω φ̇` (Tao hand-derives in Eq. 13) emerges from
rigid-body Euler dynamics. **Honest caveats** (not geometry — must be added as force
elements): bearing reactions, unbalance reaction, the grinding-force law; the DC-motor power
balance is outside GA; the stiff bearing springs (Kz ~ 1e6 N/m) make the ODE **stiff** → RK4
needs a small dt or an implicit integrator.

## Intermediate library step: refine `is_congruent` [DONE 2026-06-14]

Surfaced while porting mini-test 1: comparing two **independently computed** planes/lines
(results of wedge / dual / motor chains) with `is_congruent` fails at the default
`tolerance = eps`, even when the operands are identical to ~5 decimals. Two compounding
causes:

1. **Default tolerance is single-op machine eps**, but congruence is almost always applied
   to *computed* geometry whose error has accumulated over many GA ops. (CLAUDE.md's
   congruence note already says "absolute eps is too strict after GA ops that accumulate
   error" — that reasoning should extend to the **default**, not only the relative scaling.)
2. The algorithm fixes the scale factor `k` from the **first non-zero** component; when that
   component is small (e.g. a 0.01 tilt term next to a ~1 term), its relative rounding is
   amplified onto the dominant components and lands right at the eps boundary.

**Refinement [DONE]:** all 11 graded overloads (EGA2D/3D, PGA2DP/3DP, STA4D) now delegate to
a shared `hd::ga::detail::coeffs_congruent<N>` in `ga/ga_value_t.hpp` that (a) pivots `k` on
the **largest-magnitude** component of `b`, and (b) defaults to the new `eps_congruent`
(`1e-12`) relative tolerance reflecting accumulated GA error (scalars/pscalars keep their
trivial overloads + `eps`). The translated grinding test reverted to plain `is_congruent(...)`.
All 5 C++ suites + 725 ga_py tests pass; CLAUDE.md congruence section updated.

## Plan (phased) — decisions taken 2026-06-14

- **A — force-element extension.** Linear spring + damper generalized forces on joint
  coordinates, plus an applied external wrench at a body-fixed point, folded into `tau`.
  Additive; the gravity-only path stays byte-unchanged. (Mirrors how the closed-loop layer
  was added.)
- **B — Sommerfeld warm-up [CHOSEN to do first].** Unbalanced rotor on anisotropic x,y
  spring/damper supports, **prescribed spin**; validate steady-state amplitudes against the
  **closed-form Eq. (4)** of the Sommerfeld paper (Bisoi et al., *Advances in Rotor
  Dynamics*, 2020) and the resonances `ω₁,₂ = √(K/(m+M))` (= 16.51 / 28.57 rad/s for the
  paper's Table-1 params). The analytic correctness gate, exactly like the four-bar angles.
  - **B.1 (next concrete step):** explicit rotating unbalance wrench `meω²(cos,sin)` +
    springs/dampers → match Eq. (4) over an ω sweep. Isolates the spring/damper primitive.
  - **B.2:** replace the explicit forcing with an **emergent** spinning offset-cm body →
    must reproduce B.1 (validates "`meω²` for free" + the prescribed-spin bias path).
- **C — Tao 5-DOF spindle.** The 6-joint stack above; reproduce Eq. (13) responses and the
  Fig. 4 characteristic frequencies `f_x, f_z, f_θ, f_φ`; derive `z_b` directly from a rim
  point; **calibrate/validate against the machine's measured data**. Port mini-test 1 (tilt
  now driven by the dynamics, not a prescribed φ-sweep).
- **D — grinding contact + two-way coupling [the actual goal].** Promote mini-test 2 to a
  contact/observation layer; add a grinding-force law (MRR × relative-speed, calibrated to
  measured data) → wrench on the spindle; close the loop vibration → `z_b` → depth → MRR →
  force → vibration. Observe surface waviness (λ_c, λ_m) vs. measured topography.
- **E — optimization [visionary].** With the calibrated closed-loop model, optimize feed
  rate / speed ratio to bound the grinding force and `z_b` (damage prevention).

### Decisions

- **Non-ideal DC motor / Sommerfeld speed-jump (two-way drive↔structure): DEFERRED** until
  the grinding-force feedback (D) needs a reacting drive. Warm-up uses prescribed spin.
- **Code home: PROTOTYPE the force elements in a `ga_appl3dp` app-test first**, validate vs.
  Eq. (4), THEN promote the proven primitive into the library. (Matches how the two
  mini-tests started.)
- **Warm-up params:** use the Sommerfeld paper's Table 1 (they match the closed-form gate);
  bring the machine's measured stiffness/damping/eccentricity in at Phase C.

## Validation gates (per phase)

- B: simulated amplitudes == Eq. (4) across an ω sweep; resonances at 16.51 / 28.57 rad/s.
- C: characteristic frequencies `f_x, f_z, f_θ, f_φ` and the `z_b` waveform vs. measured data.
- D: waviness wavelengths λ_c, λ_m vs. measured wafer topography.

## Out of scope (for now)

The full DC-motor/bond-graph drive model; the grain-statistics surface-topography simulation
as a standalone renderer (we observe `z_b`/trajectories, not a full grain ensemble); any 3D
`ga_view` rendering of the spindle (ga_view is a 2D viewer).
