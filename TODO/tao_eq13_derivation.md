# Deriving Tao's Eq. (13) from the GA rigid-body model

Resolves the two deviations flagged in Phase C of the wafer-grinding work (see
[grinding.md](grinding.md) and the `PGA3DP: Tao wheel-spindle (Phase C)` suite in
`ga_test/src/ga_appl3dp_physics_test.hpp`). Reference: Tao et al., *Int. J. Mech.
Sci.* 232 (2022) 107620, §2.1, Eqs. (1)–(3), (8), (12), (13), (16)–(19), Table 1.

Both deviations are confirmed **paper errors / over-lumpings**, not coordinate-choice
artifacts and not gaps in the GA model. The GA `dynamic_system3dp` spindle reproduces
the physically correct coefficients; Tao's hand-lumped scalar Eq. (13) does not reduce
from consistent rigid-body mechanics.

## Setup (Tao's own kinematics)

Tao's angle parameterization is fine; the errors are downstream of it. The body-frame
angular velocity (his Eq. 2, with `Ω = ωt`, `ω̇ = 0`) is

    ω_x = θ̇ cosΩ + φ̇ sinΩ
    ω_y = −θ̇ sinΩ + φ̇ cosΩ
    ω_z = ω

which inverts to

    θ̇ = ω_x cosΩ − ω_y sinΩ
    φ̇ = ω_x sinΩ + ω_y cosΩ                                         (★)

Symmetric rotor: `J_x = J_y ≡ J_t` (transverse), `J_z` (polar). Ignore the
off-diagonal `−me(l1−l2)` of Eq. (3) for the gyroscopic argument — it produces only
the `∓meω²(l1−l2)` unbalance forcing on the RHS, which the GA model reproduces exactly
(Phase C.3 centrifugal gate: `ẍ = eω²`).

The body-frame Euler equations `J·ω̇_c + ω_c × (J·ω_c) = M` are, with `J_x = J_y = J_t`:

    J_t ω̇_x + (J_z − J_t) ω ω_y = M_x
    J_t ω̇_y + (J_t − J_z) ω ω_x = M_y

## Deviation 2 — gyroscopic coupling: `J_z`, not `(J_y − J_z)`

To get the θ-equation, project the Euler equations the same way Tao's angles are
defined: `(θ-eqn) = (x-eqn)·cosΩ − (y-eqn)·sinΩ`.

**Inertia term.** Differentiating `ω_x, ω_y` (note `Ω̇ = ω`):

    ω̇_x = θ̈ cosΩ − θ̇ ω sinΩ + φ̈ sinΩ + φ̇ ω cosΩ
    ω̇_y = −θ̈ sinΩ − θ̇ ω cosΩ + φ̈ cosΩ − φ̇ ω sinΩ

so

    ω̇_x cosΩ − ω̇_y sinΩ = θ̈ + ω φ̇                                  (†)

The `θ̈` is expected; the extra `+ω φ̇` is the **kinematic frame-rotation term** —
it comes from differentiating the rotating cosΩ/sinΩ projection, NOT from the cross
product. This is the term Tao drops.

**Gyroscopic (cross-product) term**, using (★):

    (J_z − J_t) ω ω_y cosΩ − (J_t − J_z) ω ω_x sinΩ
      = (J_z − J_t) ω (ω_y cosΩ + ω_x sinΩ)
      = (J_z − J_t) ω φ̇

**Sum** (LHS of the θ-equation):

    J_t (θ̈ + ω φ̇) + (J_z − J_t) ω φ̇  =  J_t θ̈ + J_z ω φ̇

The two `ωφ̇` contributions combine: the transverse `J_t` from (†) cancels the `−J_t`
in the cross-product, leaving the **polar moment `J_z`**:

    J_t θ̈ + d_x θ̇ + K_tilt·θ = −J_z ω φ̇ + (forcing)                 (CORRECT)

This is the textbook gyroscopic rotor equation — coupling coefficient = polar moment
of inertia `J_z`. The φ-equation gives `+J_z ω θ̇` symmetrically.

**Where Tao's `(J_y − J_z)` comes from.** He used `J_t(ω̇_x cosΩ − ω̇_y sinΩ) = J_t θ̈`
— i.e. dropped the `+J_t ω φ̇` of (†) — keeping only the cross-product gyroscopic term
`(J_z − J_t)ωφ̇` and moving it to the RHS:

    J_t θ̈ + ... = (J_t − J_z) ω φ̇ + ... = (J_y − J_z) ω φ̇ + ...     (Tao Eq. 13)

So Tao's coefficient is short by exactly the dropped `J_t ω φ̇`:
`(J_y − J_z) = J_t − J_z` instead of the correct `J_z`.

**Numerical consequence.** Forward/backward tilt-whirl split:

- Correct (`J_z`):  `(J_z / 2J_t) ω` ⇒ peaks split by `(J_z/J_t) f_s ≈ 1.76 Hz`
  at `f_s = 50 Hz` (essentially coincident — Tao's 6727/6680 would be ~6701/6702).
- Tao (`≈ J_t`):    `(J_t / 2J_t) ω = ω/2` ⇒ split `= ±f_s/2 ≈ ±25 Hz` (50 Hz apart),
  giving his Eqs. (18)/(19) `± f_s/2` terms and the resolved twin peaks.

GA model (Phase C.3 gate): `φ̈/θ̇ = (J_z/J_x)ω = 11.0538 s⁻¹` exactly → 1.76 Hz split,
matching standard rotordynamics, contradicting Tao.

**Verdict: paper error.** A dropped kinematic term in `d/dt(J ω_c)`, not a coordinate
artifact (Tao's Eq. 1/2 are self-consistent; the error is in carrying the time
derivative through to `θ̈, φ̈`).

## Deviation 1 — tilt stiffness: factor 4 over-count

Tao's tilt stiffness (Eqs. 12, 18, 19) is

    K_tilt = (k_x + k_y)(l1² + l2²) + k_z R_m²

**Radial part, from a consistent two-bearing model.** Two radial bearings at axial
stations `±l_i`. The translational equation (Eq. 13, line 1) reads `m ẍ + k_x x`, so
`k_x` is the **total** radial stiffness = sum of the two parallel bearings; each bearing
carries `k_x/2` per direction. (The `f_x` cross-check pins this: `f_x = (1/2π)√(k_x/m)`
matches 2385 Hz only with `k_x = 1.8e8 N/m` as the *total* — if `k_x` were per-bearing,
`f_x` would be a factor √2 high.)

A tilt `φ` about the y-axis displaces a bearing point at axial `l_i` by `l_i φ` **in x
only** (rotation about y maps z→x), engaging only the x-springs. Restoring torque about
y from each bearing: `(k_x/2) l_i² φ`. Summing the two bearings:

    K_tilt^radial (consistent) = (k_x/2)(l1² + l2²)

Compare Tao's radial part `(k_x + k_y)(l1² + l2²) = 2 k_x (l1² + l2²)` (since
`k_x = k_y`). The ratio is **4**, decomposing as two independent factors of 2:

1. `(k_x + k_y)` vs `k_x` — Tao adds the orthogonal-direction stiffness `k_y`, but a
   tilt about y is resisted only by the x-springs. (×2)
2. `k_x` vs `k_x/2` per bearing — Tao uses the *total* radial stiffness as each
   bearing's contribution rather than distributing it over the two parallel bearings.
   (×2)

Equivalently: Tao's coefficient treats each of the two bearings as having the full
radial stiffness in **both** directions, where the geometry justifies only (full/2)
per bearing in the **one** direction that resists the tilt.

**Axial part.** `k_z R_m²` likely carries a similar ×2 (a distributed axial bearing at
radius `R_m` gives `(k_z/2)R_m²` from the cos² average), but the bearing layout for the
axial film is less explicit in the paper, so this term is left as-is in the model.

**Numerical consequence / the `Lb` calibration.** At the table's `l1 = l2 = 0.1 m`, a
consistent radial model gives `1/4` of Tao's radial tilt stiffness; keeping `k_z R_m²`
as-is, `f_θ` comes out ≈ 4430 Hz (~34% below Tao's ≈ 6700). The Phase C model absorbs
the gap by inflating the radial-bearing lever to

    Lb = √(K_tilt / k_x) ≈ 0.231 m  (≈ 2.3 × the physical l1 = 0.1 m)

with two springs of `k_x/2` at `±Lb`, so the emergent tilt stiffness
`2·(k_x/2)·Lb² = k_x Lb² = K_tilt` reproduces Tao's `f_θ` by construction.

**Verdict: paper over-count.** `Lb` is therefore *not* a physical bearing offset — it
is the value that reproduces Tao's published `f_θ`, which is itself 4× too stiff in its
radial part. A first-principles model at the physical `l1 = l2 = 0.1 m` predicts the
real tilt frequency ≈ 34% lower. The `Lb` calibration is documented (not "mysterious"):
it equals `√(4·(physical l1²+l2²)/2 + k_zR_m²/k_x)` ≈ how the 4× shows up as a ≈2.3×
lever.

## Where the corrected model predicts something different from the paper

Both errors live entirely in the **tilt (θ, φ)** dynamics; the translational (x, y, z)
equations of Eq. (13) are untouched. So the corrected model agrees with the paper
everywhere except the figures/quantities that come from *solving the tilt equations*.

**Differs (tilt-only):**

- **Fig. 4(c),(d) — θ/φ PSD.** Tao shows twin peaks at f_θ = 6727 and f_φ = 6680 Hz
  (47 Hz apart). Corrected: (i) the center drops — a consistent two-bearing tilt
  stiffness is 4× smaller in its radial part, so with the table's parameters f_tilt ≈
  3.7–4.4 kHz, not ~6.7 kHz; (ii) the 47 Hz separation collapses — the genuine
  gyroscopic split is (J_z/J_t) f_s ≈ 1.76 Hz, unresolvable, so a *single* tilt peak,
  not a pair.
- **Fig. 6(b) — z_b PSD** (z_b = z − R_w φ carries the tilt frequencies): same — single
  lower tilt peak instead of the twin 6680/6727.
- **Fig. 15 — simulated acceleration** ("acquired via Eqs. 13–14"): the plotted beat is
  the f_θ–f_φ pair beating at ~42 Hz (~24 ms envelope). Corrected, the two tilt modes are
  ~1 Hz apart, so that beat is essentially gone on the plotted scale.
- **Eqs. (18),(19),(22):** `(k_x+k_y)(l1²+l2²) → (k_x/2)(l1²+l2²)` (radial part 4× smaller)
  and `±f_s/2 → ±(J_z/2J_t) f_s` (negligible).

**A consistency check the correction passes and the paper's Eq. (13) fails.** Tao's
±f_s/2 makes the tilt-peak *separation* grow with spin (37.5 / 44.2 / 50 Hz at
n_w = 2250 / 2650 / 3000 r/min). His own measurement (§4.2.1) says the two tilting peaks
(6103, 6403 Hz) are "the same under n_w = 2250 and 2650 r/min" — i.e. **spin-independent**.
A real gyroscopic split scales with spin; a structural/anisotropic (k_x ≠ k_y, asymmetric
inertia) split does not. So the measured spin-independence is consistent with the
corrected picture (gyroscopic split negligible; the ~300 Hz measured separation is a
bearing/structure anisotropy, not whirl) and *inconsistent* with Eq. (13)'s gyroscopic
±f_s/2 origin. Tao does not flag this.

**An honest tension the correction exposes (not a win for our model).** The corrected
rigid-body tilt frequency (~3.7–4.4 kHz at the table's stiffnesses) sits ~30–40% *below*
both Tao's ~6.7 kHz and the measured ~6.25 kHz. Tao's 4× radial over-count fortuitously
lifts his theoretical f_θ up near the measured value. The lesson: the published
*translational* bearing stiffnesses cannot be naively reduced to a tilt stiffness — the
real aerostatic film is stiffer in tilt than `(k_x/2)(l1²+l2²) + (k_z/2)R_m²` implies. A
faithful model should take the tilt stiffness from a tilt measurement (or a film model),
not inherit it from k_x, k_z; Eq. (12)'s lumped coefficient hides this behind a wrong
factor.

**Does NOT differ (so the paper's headline results stand):**

- **Figs. 3(a), 4(a),(b), 5 — translational x, y, z** and the radial centrifugal/`f_x`,
  `f_z` responses: unaffected (the errors are tilt-only). The GA model reproduces
  f_x = 2387, f_z = 5485 Hz from first principles.
- **§3 topography (Figs. 9–12) and §4 experimental wavelength validation (Table 3,
  λ_c/λ_m, the 11.7%/8.8% agreement):** the surface-formation model takes the vibration
  as an *input* sine `z_b = A_b sin(2π f_b t + φ_b)`, and the experimental comparison uses
  the **measured** tilting frequency (mean 6253 Hz, §4.2.2) as f_b — never Eq. (22)'s
  predicted value. So the topography results and their agreement with experiment do not
  depend on the Eq. (13) errors and are unchanged.

Net: the corrections rewrite Tao's *predicted spindle-dynamics spectra* (Figs. 4c/d, 6,
15; Eqs. 18/19/22), reconcile the model with his own spin-independence observation, and
leave the experimentally-validated surface-topography half of the paper untouched.

## Status of the GA model

No code change is required for correctness — the GA `dynamic_system3dp` already gives
the physically correct gyroscopic coefficient (`J_z`) and the correct emergent tilt
stiffness `k·l²` per spring. The only modeling concession is the calibrated `Lb`, which
is now explained as absorbing Tao's 4× radial over-count rather than being unexplained.

Open follow-up: continue the roadmap — contact/force loop (grinding force ↔ vibration ↔
`z_b` ↔ cutting depth ↔ MRR, Eqs. 14/26) and feed/speed optimization.
