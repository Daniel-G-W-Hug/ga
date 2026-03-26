# Coupled Newton-Euler Equations for Rigid Body Dynamics

## Why Coupling Arises

If you express everything in an **inertial frame**, translation and rotation *seem*
decoupled:

$$m\ddot{\mathbf{r}} = \mathbf{F}, \qquad \dot{\mathbf{L}} = \boldsymbol{\tau}$$

But the moment you work in a **body-fixed frame** (which is almost always necessary for
non-spherical bodies, since the inertia tensor $\mathbf{I}$ is then constant), the
rotation of the frame itself introduces coupling terms.

---

## The Full Coupled System

**Translational equation (body frame):**

$$m(\dot{\mathbf{v}} + \boldsymbol{\omega} \times \mathbf{v}) = \mathbf{F}^b$$

**Rotational equation (Euler, body frame):**

$$\mathbf{I}\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times
(\mathbf{I}\boldsymbol{\omega}) = \boldsymbol{\tau}^b$$

where:
- $\mathbf{v}$ — translational velocity in the body frame
- $\boldsymbol{\omega}$ — angular velocity in the body frame
- $\mathbf{I}$ — (constant) inertia tensor in the body frame
- $\mathbf{F}^b, \boldsymbol{\tau}^b$ — external force/torque expressed in the body frame

The **coupling term** $\boldsymbol{\omega} \times \mathbf{v}$ in the translational
equation arises from the rotating reference frame.

---

## Solving for the Derivatives (state-space form)

$$\dot{\mathbf{v}} = \frac{\mathbf{F}^b}{m} - \boldsymbol{\omega} \times \mathbf{v}$$

$$\dot{\boldsymbol{\omega}} = \mathbf{I}^{-1}\left(\boldsymbol{\tau}^b -
\boldsymbol{\omega} \times \mathbf{I}\boldsymbol{\omega}\right)$$

The state vector is $[\mathbf{v},\, \boldsymbol{\omega}]^T \in \mathbb{R}^6$.

---

## Also Integrate the Orientation

To recover position/orientation in the inertial frame:

$$\dot{\mathbf{r}} = \mathbf{R}\,\mathbf{v}$$

$$\dot{\mathbf{R}} = \mathbf{R}\,[\boldsymbol{\omega}]_\times \quad \text{(or use
quaternions: } \dot{\mathbf{q}} = \tfrac{1}{2}\mathbf{q} \otimes
\boldsymbol{\omega}\text{)}$$

**Quaternions are strongly recommended** to avoid gimbal lock.

---

## Complete State Vector for Integration

| State | Size | Equation |
|---|---|---|
| $\mathbf{r}$ | 3 | $\dot{\mathbf{r}} = \mathbf{R}\mathbf{v}$ |
| $\mathbf{q}$ (quaternion) | 4 | $\dot{\mathbf{q}} = \frac{1}{2}\mathbf{q}\otimes\boldsymbol{\omega}$ |
| $\mathbf{v}$ (body frame) | 3 | $\dot{\mathbf{v}} = \mathbf{F}^b/m - \boldsymbol{\omega}\times\mathbf{v}$ |
| $\boldsymbol{\omega}$ (body frame) | 3 | $\dot{\boldsymbol{\omega}} = \mathbf{I}^{-1}(\boldsymbol{\tau}^b - \boldsymbol{\omega}\times\mathbf{I}\boldsymbol{\omega})$ |

**Total: 13 coupled 1st-order ODEs** (with quaternion normalization constraint
$|\mathbf{q}|=1$).

---

## Key Differences from a "Parallel" Approach

| | Parallel (inertial frame) | Coupled (body frame) |
|---|---|---|
| Inertia tensor $\mathbf{I}$ | Time-varying $\mathbf{I}(t) = \mathbf{R}\mathbf{I}_0\mathbf{R}^T$ | Constant |
| Coupling term | None explicitly | $\boldsymbol{\omega}\times\mathbf{v}$ and $\boldsymbol{\omega}\times\mathbf{I}\boldsymbol{\omega}$ |
| Practical for non-spherical bodies | Expensive | Standard choice |

The $\boldsymbol{\omega}\times\mathbf{I}\boldsymbol{\omega}$ term is already present in
the Euler equations. The only *new* coupling when adding translation is the
$\boldsymbol{\omega}\times\mathbf{v}$ term, which accounts for the rotating body frame.

---

## Component-wise Matrix Form

The coupled equation component-wise as vector equation with matrices shows how the state
vector components are actually coupled.

### Full Component-wise Newton-Euler System

Using **principal body axes** (diagonal $\mathbf{I} = \text{diag}(I_x, I_y, I_z)$), the
system in state-space form $\dot{\mathbf{x}} = \mathbf{A}(\boldsymbol{\omega})\,\mathbf{x}+ \mathbf{b}$ is:

$$\begin{pmatrix} \dot{v}_x \\ \dot{v}_y \\ \dot{v}_z \\ \dot{\omega}_x \\ \dot{\omega}_y
\\ \dot{\omega}_z \end{pmatrix} = \left(\begin{array}{ccc|ccc} 0 & +\omega_z & -\omega_y &
0 & 0 & 0 \\ -\omega_z & 0 & +\omega_x & 0 & 0 & 0 \\ +\omega_y & -\omega_x & 0 & 0 & 0 &
0 \\ \hline 0 & 0 & 0 & 0 & +\frac{I_y}{I_x}\omega_z & -\frac{I_z}{I_x}\omega_y \\ 0 & 0 &
0 & -\frac{I_x}{I_y}\omega_z & 0 & +\frac{I_z}{I_y}\omega_x \\ 0 & 0 & 0 &
+\frac{I_x}{I_z}\omega_y & -\frac{I_y}{I_z}\omega_x & 0 \end{array}\right) \begin{pmatrix}
v_x \\ v_y \\ v_z \\ \omega_x \\ \omega_y \\ \omega_z \end{pmatrix} + \begin{pmatrix}
F_x^b/m \\ F_y^b/m \\ F_z^b/m \\ \tau_x^b/I_x \\ \tau_y^b/I_y \\ \tau_z^b/I_z
\end{pmatrix}$$

---

### Block Structure of $\mathbf{A}(\boldsymbol{\omega})$

$$\mathbf{A}(\boldsymbol{\omega}) = \begin{pmatrix}
\underbrace{-[\boldsymbol{\omega}]_\times}_{\text{frame rotation}} & \mathbf{0} \\
\mathbf{0} &
\underbrace{-\mathbf{I}^{-1}[\boldsymbol{\omega}]_\times\mathbf{I}}_{\text{gyroscopic}}
\end{pmatrix}$$

| Block | Meaning |
|---|---|
| **Top-left** $-[\boldsymbol{\omega}]_\times$ | $v_x, v_y, v_z$ couple to each other — rotating frame drags translational velocity |
| **Bottom-right** $-\mathbf{I}^{-1}[\boldsymbol{\omega}]_\times\mathbf{I}$ | $\omega_x, \omega_y, \omega_z$ couple to each other — gyroscopic/Euler torque |
| **Off-diagonal** $\mathbf{0}$ | No direct $v_i \leftrightarrow \omega_j$ coupling in the state |

---

### Key Observations

- The two subsystems ($\mathbf{v}$ and $\boldsymbol{\omega}$) **do not appear in each
  other's state slots** — off-diagonal blocks are zero.
- The coupling is **parametric**: $\boldsymbol{\omega}$ appears as *coefficients* in the
  translational block.
- $\boldsymbol{\omega}$ drives $\mathbf{v}$ indirectly via the rotating frame;
  $\mathbf{v}$ does **not** influence $\dot{\boldsymbol{\omega}}$.
- The system is **nonlinear** because $\mathbf{A}$ itself depends on the state
  $\boldsymbol{\omega}$.
- Since $\boldsymbol{\omega}$ must be known to evaluate $\mathbf{A}$, both subsystems
  **must be integrated together** with a shared $\boldsymbol{\omega}$ at each timestep —
  they cannot be decoupled.

---

## 2D Planar Case (2 Translational + 1 Rotational DOF)

In 2D, only $\omega_z = \omega$ survives (rotation about the out-of-plane axis), and
$\omega_x = \omega_y = 0$. The cross product $\boldsymbol{\omega} \times \mathbf{v}$
reduces to a simple 2D rotation, and the gyroscopic term $\boldsymbol{\omega} \times
\mathbf{I}\boldsymbol{\omega}$ vanishes entirely.

### Component-wise Matrix Form

$$\begin{pmatrix} \dot{v}_x \\ \dot{v}_y \\ \dot{\omega} \end{pmatrix} = \begin{pmatrix} 0
& +\omega & 0 \\ -\omega & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} v_x \\ v_y \\
\omega \end{pmatrix} + \begin{pmatrix} F_x^b/m \\ F_y^b/m \\ \tau^b/I_z \end{pmatrix}$$

---

### Reading the Structure

| Block | Meaning |
|---|---|
| Top-left $2\times2$: $\begin{pmatrix}0 & +\omega \\ -\omega & 0\end{pmatrix}$ | $v_x, v_y$ couple to each other via frame rotation — same $-[\boldsymbol{\omega}]_\times$ as 3D, just reduced |
| Bottom-right scalar: $0$ | **No gyroscopic term** — $\boldsymbol{\omega}\times\mathbf{I}\boldsymbol{\omega} = 0$ since there is only one rotation axis |
| Off-diagonal column 3, rows 1–2: $0$ | $\omega$ still does not appear as a *state* in the $\dot{\mathbf{v}}$ equation — it is a *coefficient* only |
| Off-diagonal row 3, cols 1–2: $0$ | $v_x, v_y$ have **no influence** on $\dot{\omega}$ — same as in 3D |

The rotational equation degenerates to the trivially simple:

$$\dot{\omega} = \frac{\tau^b}{I_z}$$

---

### Kinematic Equations (inertial frame)

$$\begin{pmatrix} \dot{x} \\ \dot{y} \\ \dot{\theta} \end{pmatrix} = \begin{pmatrix}
\cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{pmatrix}
\begin{pmatrix} v_x \\ v_y \\ \omega \end{pmatrix}$$

---

### Complete 2D State Vector

| State | Size | Equation |
|---|---|---|
| $x, y$ | 2 | $\dot{\mathbf{r}} = \mathbf{R}(\theta)\,\mathbf{v}$ |
| $\theta$ | 1 | $\dot{\theta} = \omega$ |
| $v_x, v_y$ | 2 | $\dot{v}_x = +\omega v_y + F_x^b/m$, $\;\dot{v}_y = -\omega v_x + F_y^b/m$ |
| $\omega$ | 1 | $\dot{\omega} = \tau^b/I_z$ |

**Total: 6 coupled 1st-order ODEs** — compared to 13 in 3D.

---

### Comparison: 2D vs 3D

| Property | 3D | 2D |
|---|---|---|
| State size | 13 (with quaternion) | 6 |
| Gyroscopic term $\boldsymbol{\omega}\times\mathbf{I}\boldsymbol{\omega}$ | Present, nonlinear | **Zero** |
| Rotational DOF | 3 ($\omega_x, \omega_y, \omega_z$ coupled) | 1 ($\omega$ scalar, decoupled) |
| Orientation representation | Quaternion $\mathbf{q}$ (4 states) | Scalar angle $\theta$ (1 state) |
| $\dot{\omega}$ equation | Nonlinear Euler equations | Simple $\tau^b/I_z$ |

---

## Constrained Motion: Rotation About a Fixed Pivot (not the CM)

### The Core Confusion: Frame vs. Rotation Axis

When a body is constrained so that a specific body-frame point $Q_b$ (the pivot) remains
fixed in the world, a subtle but important correction is required. It is easy to conflate
two independent questions:

1. **In which frame are the equations expressed?** (body frame vs. world frame)
2. **About which axis does the body actually rotate?** (cm vs. pivot)

These are completely independent. The correction for the rotation axis (parallel axis
theorem) must be applied regardless of which frame you use to express the equations.

---

### Why the Inertia About the CM Is Not Sufficient

The standard Euler equation `compute_omega_dot` computes:

$$\dot{\boldsymbol{\Omega}} = \mathbf{I}_{cm}^{-1}\bigl(\mathbf{F} -
\text{rcmt}(\boldsymbol{\Omega},\,\mathbf{I}_{cm}\boldsymbol{\Omega})\bigr)$$

This is the correct equation for a **free** rigid body rotating about its own center of
mass. It gives:

$$\alpha_{\text{free}} = \frac{\tau_{\text{eff}}}{I_{cm,\text{rot}}}$$

But when the pivot is fixed (not the cm), the body is **not free**. The cm follows a
circular arc of radius $r$ around the pivot. This circular motion requires centripetal
acceleration that is provided by the constraint reaction at the pivot. The constraint
effectively increases the inertia the applied torque must overcome.

The correct equation for a pivot-constrained body is:

$$\alpha_{\text{pivot}} = \frac{\tau_{\text{eff}}}{I_{\text{pivot}}} =
\frac{\tau_{\text{eff}}}{I_{cm,\text{rot}} + m r^2}$$

where $r$ is the distance from the cm to the pivot (parallel axis theorem / Steiner's
theorem). Since $\tau_{\text{eff}}$ is the same in both (a pure force couple has equal
torque about any reference point), the two are related by a simple scalar correction:

$$\alpha_{\text{pivot}} = \alpha_{\text{free}} \cdot \frac{I_{cm,\text{rot}}}{I_{\text{pivot}}}$$

---

### Why the Integration Frame Is Irrelevant

The integration runs entirely in the body frame: `B_b = phi * Q_b` and `Omega_b = omega *
Q_b` are body-frame quantities. But `B_b` and `Omega_b` are not describing motion *within*
the body frame — they are body-frame *coordinates* of the motor's generator, which
describes how the body frame moves relative to the world.

The correction `I_{cm}/I_{\text{pivot}}` is **not** a frame transformation. It is a
correction for the physical rotation axis. Whether you express the equations in the body
frame or the world frame, the body is rotating about an axis that passes through $Q_b$, at
distance $r$ from the cm. The kinetic energy is $\frac{1}{2} I_{\text{pivot}} \omega^2$
regardless of frame; therefore the equation of motion requires $I_{\text{pivot}}$,
regardless of frame.

**Summary:** The frame choice affects how quantities are *expressed*; the rotation axis
choice affects the *physics*. The parallel-axis correction addresses the physics.

---

### Implementation in PGA2DP Body-Frame ODE (2D Plate Pendulum)

The parallel-axis correction is applied once at construction time by passing the pivot
point directly to `get_plate_inertia`. In `calc_rhs` no manual scaling is needed:

```
// Constructor (once):
Q_b  = Vec2dp{hw, hh, 1}                         (pivot = TR corner in body frame)
I    = get_plate_inertia(m, w, h, Q_b)            (inertia already about pivot Q_b)
I_inv = get_inertia_inverse(I)

// calc_rhs (every step):
Omega_b = omega * Q_b         (constraint: B_b ∝ Q_b, rotation about pivot Q_b)

dOmega = compute_omega_dot(I_inv, F_b, Omega, I)
         // I is about Q_b, so dOmega.z = tau_eff / I_pivot[2,2] directly

alpha  = dOmega.z             // no further scaling needed
rhs[1] = alpha * Q_b          // enforce constraint: dOmega_b/dt ∝ Q_b
```

The force computation must briefly pass through the world frame to locate the force
application points (`pivot_w` and `cm_w`), but the torque bivector is immediately
transformed back to the body frame via `F_b = move2dp(F_up_w + F_dn_w, rrev(M))`.
All ODE state (`B_b`, `Omega_b`, `rhs`) is body-frame throughout.

---

### Parallel-Axis Correction: Pivot Point (2D) vs. Pivot Line (3D)

The correction is built directly into the body inertia constructors via an optional pivot
parameter. Only the **scalar rotational entries** of the inertia matrix are updated — the
mass block (upper-left) and off-diagonal coupling terms are left unchanged:

**PGA2DP — pivot is a point** `Vec2dp Q_b`:

The pivot point $Q_b = (P_x, P_y, 1)$ gives the scalar Steiner correction directly:

$$I_{\text{pivot}}[2,2] = J_{cm} + m\,(P_x^2 + P_y^2)$$

All other entries of $\mathbf{I}$ remain equal to $\mathbf{I}_{cm}$.

```cpp
get_plate_inertia(m, w, h, Q_b)   // Q_b = Vec2dp{hw, hh, 1} — TR corner
```

Default: `Q_b = Vec2dp{0, 0, 1}` (body origin = cm) → no correction.

**Why off-diagonal coupling terms must not be added.** The body frame always has the cm
at its origin, so the general off-diagonal formula for a body whose cm is displaced does
not apply. Adding coupling terms via `get_point_inertia(m, Q_b)` would cause
$\mathbf{I}_{\text{pivot}} \cdot (\omega Q_b) = (0,\,0,\,J_{cm}\,\omega)$ instead of
$(0,\,0,\,I_{\text{pivot},zz}\,\omega)$: the Steiner term $mr^2$ is silently cancelled,
and `I_inv[2,2] = 1/J_cm` rather than `1/I_pivot_zz` — an error of factor
$I_{\text{pivot},zz}/J_{cm}$ (about $4$ for a square plate rotating about its corner).

**PGA3DP — pivot is a line** `BiVec3dp L_pivot`:

In 3D the pivot entity is a line (grade 2). The foot of the perpendicular from the body
origin $O_b$ to $L_{\text{pivot}}$ gives the correction point:

$$\mathbf{n} = (L.vx,\, L.vy,\, L.vz), \quad
  \text{mom} = (L.mx,\, L.my,\, L.mz), \quad
  \mathbf{P}_{\text{foot}} = \frac{\mathbf{n} \times \text{mom}}{|\mathbf{n}|^2}$$

Only the diagonal lower-left block of the $6\times 6$ inertia matrix is updated:

$$I[3,0] \mathrel{+}= m(P_y^2+P_z^2), \quad
  I[4,1] \mathrel{+}= m(P_x^2+P_z^2), \quad
  I[5,2] \mathrel{+}= m(P_x^2+P_y^2)$$

```cpp
get_cuboid_inertia(m, w, h, d, L_pivot)   // L_pivot = BiVec3dp pivot line
```

Default: `L_pivot = BiVec3dp{}` (zero ideal part $\mathbf{n} = 0$, line through origin)
→ no correction. This correctly handles both the zero-bivector sentinel and any valid
pivot line that passes through the cm.

**Constraint projection is always required.** The ODE state is constrained to
$\boldsymbol{\Omega}_b \propto$ pivot entity, but `compute_omega_dot` returns a full
unconstrained result (it knows nothing about the kinematic constraint). The projection
extracts the scalar angular acceleration $\alpha$ and re-aligns the result:

- 2D: `alpha = dOmega.z` (since $Q_b.z = 1$, $\omega = \Omega_b.z$)
- 3D: `alpha = dOmega · L̂` (project result onto normalized pivot line direction)

Then `rhs[1] = alpha * Q_b` (2D) or `rhs[1] = alpha * L_pivot` (3D) enforces the
constraint $\boldsymbol{\Omega}_b \propto$ pivot entity.

---

## When the Pivot Is Not Fixed: Multi-DOF Systems and the Instantaneous Center of Rotation

### Motivation: Thread-Connected Pivot

Consider replacing the fixed pivot with a thread of length $l$ attached to a fixed anchor
$A$ above. The TR corner $Q_b$ now hangs from $A$ and can swing freely. The constraint
changes fundamentally:

- **Fixed pivot**: $Q_b$ has **zero velocity** in the world. The reaction force is
  determined implicitly — it takes whatever value keeps $Q_b$ stationary. The system has
  **1 DOF** ($\varphi$, the plate rotation angle), and the pivot inertia $I_{\text{pivot}}$
  fully captures the dynamics.
- **Thread pivot**: $Q_b$ moves on a circle of radius $l$ around $A$. The thread provides
  a constraint force **along** the thread direction, with a magnitude that depends on
  centripetal acceleration. The system now has **2 DOF**: $\varphi_1$ (thread swing angle)
  and $\varphi_2$ (plate rotation about $Q_b$).

The single-pivot ODE formulation with a corrected $I_{\text{pivot}}$ does **not** extend
to the thread case.

---

### The Instantaneous Center of Rotation (ICR)

For any 2D rigid body in general planar motion, there exists an instantaneous center of
zero velocity at every moment. In PGA2DP this is particularly elegant: the world-frame
velocity bivector $\boldsymbol{\Omega}_{\text{world}}$ **is** the ICR as an unnormalized
point.

For the thread + plate compound pendulum, the world-frame velocity is a superposition of
the two rotation contributions:

$$\boldsymbol{\Omega}_{\text{world}} = \omega_1 \cdot A + \omega_2 \cdot Q_{b,\text{world}}$$

where $A$ is the world-fixed thread anchor point, $Q_{b,\text{world}} =
\text{move2dp}(Q_b, M_{\text{thread}}(\varphi_1))$ is the current world position of the
TR corner, and $\omega_1$, $\omega_2$ are the angular velocities of the two DOF. Since
both are normalized points ($.z = 1$):

$$\text{ICR} = \frac{\boldsymbol{\Omega}_{\text{world}}}{\omega_1 + \omega_2} =
\frac{\omega_1 \cdot A + \omega_2 \cdot Q_{b,\text{world}}}{\omega_1 + \omega_2}$$

The ICR is the **weighted average** of the two rotation centers, weighted by their
respective angular velocities. In PGA3DP the same relation holds with the ICR replaced by
an instantaneous rotation **axis** (a `BiVec3dp` line).

**Critical limitation**: computing the ICR requires knowing $\omega_1$ and $\omega_2$
separately. The ICR is a kinematic consequence of the dynamics — it cannot substitute for
a proper multi-DOF formulation. Using the ICR as a substitute pivot in the single-DOF
equation would be circular.

---

### Correct Formulation: Motor Composition and 2-DOF ODE

The natural PGA formulation decomposes the total motor into a chain:

$$M(\varphi_1, \varphi_2) = M_{\text{thread}}(\varphi_1) \mathbin{⟇} \exp\!\bigl(\tfrac{1}{2}\,\varphi_2 \cdot Q_b\bigr)$$

where $M_{\text{thread}}(\varphi_1)$ encodes the swing of $Q_b$ about anchor $A$
(rotation by $\varphi_1$ about $A$ in the world frame), and $\exp(\tfrac{1}{2}\,\varphi_2
\cdot Q_b)$ is the plate's own rotation about $Q_b$ in the body frame. The ODE state
becomes:

$$\bigl(B_1,\,\Omega_1,\,B_2,\,\Omega_2\bigr) \;\in\; \mathbb{R}^4$$

two coupled second-order ODEs (four first-order). The thread tension must be computed
explicitly from the centripetal constraint and appears as a coupling force between the two
equations.

---

### Body Frame: Still Useful, but Role Changes

| | Fixed pivot | Thread pivot |
|---|---|---|
| DOF | 1 | 2 |
| Pivot constraint | Zero velocity at $Q_b$ | $Q_b$ constrained to circle of radius $l$ |
| Reaction force | Free (implicitly absorbed) | Along thread, magnitude from centripetal eq. |
| Body-frame advantage | $\mathbf{I}$ constant, single $\alpha$ | $\mathbf{I}$ still constant for plate rotation |
| Motor structure | $M_0 \mathbin{⟇} \exp(\tfrac{1}{2} B_b)$ | $M_{\text{thread}}(\varphi_1) \mathbin{⟇} \exp(\tfrac{1}{2} \varphi_2 \cdot Q_b)$ |

The plate's inertia tensor $\mathbf{I}$ (parallel-axis corrected about $Q_b$) remains
constant in the body frame and reusable in `calc_rhs` for the plate's own rotational
equation. What changes is that `calc_rhs` must also produce $\ddot{\varphi}_1$ for the
thread, and the two equations are coupled through the thread tension — which requires
solving the constraint equations simultaneously at each step.

The motor composition structure $M_{\text{thread}} \mathbin{⟇} M_{\text{plate}}$ scales
naturally to 3D and to chains of additional bodies, making PGA the natural language for
this class of multi-body problem.
