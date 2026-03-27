# Modelling Rigid Body Motion in PGA

This document collects the transformation laws for rigid body motion modelling in
PGA2DP and PGA3DP as used in this codebase. It covers how geometric objects transform
under motors, how to construct motors, how to formulate body-frame ODEs, and how forces
and velocities map between frames.

---

## The Motor and the Sandwich Product

A **motor** $M$ encodes a rigid body configuration (rotation + translation). It acts on
geometric objects via the **sandwich product**:

$$X_w = M \mathbin{⟇} X_b \mathbin{⟇} \utilde{M}$$

In code this is written as:

```cpp
X_w = move2dp(X_b, M);   // PGA2DP
X_w = move3dp(X_b, M);   // PGA3DP
```

The notation follows the documentation convention:
$\tilde{M}$ denotes the **reverse** (`rev(M)`) and $\utilde{M}$ denotes the
**regressive reverse** (`rrev(M)`). Since this codebase uses the regressive product
convention throughout (`rgpr` for motor composition, `rrev` for the sandwich product),
all motor sandwich products use $\utilde{M}$. For unit motors the two reverses
coincide as group inverses, but the function names and symbols are kept distinct.

**Notation source by document type:** in Markdown files `\utilde{}` is used (KaTeX
native, no package needed); in the LaTeX documentation (`ga_docu`) `\undertilde{}`
from the `accents` package (`\usepackage{accents}`) is used instead.

**Inverse transformation** (world → body):

$$X_b = \utilde{M} \mathbin{⟇} X_w \mathbin{⟇} M \quad \Longleftrightarrow \quad
\texttt{move2dp}(X_w,\, \texttt{rrev}(M))$$

---

## PGA Objects and How They Transform

### PGA2DP — G(2,0,1)

| Object | Type | Homogeneous coord | Effect of $z=0$ |
|---|---|---|---|
| Finite point $(x,y)$ | `Vec2dp{x, y, 1}` | $z = 1$ | Motor translates + rotates |
| Ideal point / direction | `Vec2dp{x, y, 0}` | $z = 0$ | Motor rotates **only** (translation has no effect) |
| Line / force bivector | `BiVec2dp{x, y, z}` | — | Motor translates + rotates |

All three transform via `move2dp(X, M)`.

**Key consequence:** body-frame corner offsets stored as direction vectors
`{bx, by, 0}` are rotated but not translated by the motor — the cm translation is
added separately:

```cpp
vec2dp const offset_b{corner.x, corner.y, 0.0};     // direction (z=0)
vec2dp const offset_w = move2dp(offset_b, M);        // pure rotation
vec2dp const corner_w{cm_w.x + offset_w.x, cm_w.y + offset_w.y, 1.0};
```

### PGA3DP — G(3,0,1)

| Object | Type | Homogeneous coord |
|---|---|---|
| Finite point $(x,y,z)$ | `Vec3dp{x, y, z, 1}` | $w = 1$ |
| Ideal point / direction | `Vec3dp{x, y, z, 0}` | $w = 0$ |
| Line | `BiVec3dp{vx,vy,vz, mx,my,mz}` | ideal part $\mathbf{v}$, real part $\mathbf{m}$ |
| Plane | `TriVec3dp{...}` | — |

All transform via `move3dp(X, M)` with the same sandwich product.

---

## Motor Construction

### Pure Translation by $(t_x, t_y)$ — PGA2DP

$$M_t = \exp\!\Bigl(\tfrac{1}{2}\,\text{vec2dp}\{-t_y,\; t_x,\; 0\}\Bigr)$$

```cpp
mvec2dp_u M_t = exp(0.5 * vec2dp{-ty, tx, 0.0});
```

This places the body origin at $(t_x, t_y)$ in the world.

### Pure Rotation by $\varphi$ About Point $P = (p_x, p_y, 1)$ — PGA2DP

$$M_r = \exp\!\Bigl(\tfrac{1}{2}\,\varphi \cdot P\Bigr) = \exp\!\Bigl(\tfrac{\varphi}{2}\,
\text{vec2dp}\{p_x,\; p_y,\; 1\}\Bigr)$$

```cpp
auto M_r = exp(0.5 * phi * vec2dp{px, py, 1.0});
```

### Combined: Body-Frame Motor $M_0 \mathbin{⟇} M_r$

The motor for a body initially placed at rest by $M_0$, then rotated by $\varphi$ about
body-frame point $Q_b$ (= the pivot), is:

$$M(t) = M_0 \mathbin{⟇} \exp\!\Bigl(\tfrac{1}{2}\,B_b(t)\Bigr)$$

```cpp
auto M = rgpr(M0, exp(0.5 * B_b));
```

where `rgpr` is the regressive geometric product (motor composition in this convention).

**$M_0$ construction:** given the initial world cm position `cm_pos0 = (cx, cy, 1)`:

```cpp
mvec2dp_u M0 = exp(0.5 * vec2dp{-cm_pos0.y, cm_pos0.x, 0.0});
```

**Recovering the cm world position from the motor:**

```cpp
vec2dp const O_b{0.0, 0.0, 1.0};          // body origin = cm in body frame
auto const cm_w = move2dp(O_b, M);        // cm in world frame
```

---

## Body-Frame ODE Formulation

The ODE state is $(B_b,\,\Omega_b)$ where both are body-frame vectors (`Vec2dp` in 2D,
`BiVec3dp` in 3D):

$$\frac{d B_b}{dt} = \Omega_b, \qquad \frac{d \Omega_b}{dt} = \mathbf{I}_{\text{pivot}}^{-1}
\Bigl(\mathbf{F}_b - \text{rcmt}(\Omega_b,\,\mathbf{I}_{\text{pivot}}\,\Omega_b)\Bigr)$$

The second equation is `compute_omega_dot(I_inv, F_b, Omega, I)`.

**For pivot-constrained motion** the constraint $B_b \propto Q_b$ (2D) or $B_b \propto
L_b$ (3D) must be enforced explicitly after calling `compute_omega_dot`:

```cpp
// 2D: extract scalar alpha, project onto constraint direction
auto const dOmega = compute_omega_dot(I_inv, F_b, Omega, I);
double const alpha = dOmega.z;   // since Q_b.z = 1
rhs[1] = alpha * Q_b;
```

**Why body-frame is preferred:**
The inertia tensor $\mathbf{I}$ is *constant* in the body frame (it does not rotate with
the body). In the world frame, $\mathbf{I}(t) = M \mathbf{I}_0 \utilde{M}$ must be
recomputed at every step — expensive and numerically delicate for 3D.

---

## Force and Torque Transformation

### Constructing a Force Bivector

A force $\mathbf{f}$ applied at point $P$ produces the PGA bivector (wrench):

$$\mathbf{F} = P \wedge \mathbf{f} \quad \Longleftrightarrow \quad
\texttt{wdg}(P,\, \mathbf{f})$$

The `BiVec2dp` basis is $\{\mathbf{e}_{31},\, \mathbf{e}_{32},\, \mathbf{e}_{12}\}$, so for
$\mathbf{F} = \texttt{wdg}(\{p_x, p_y, 1\},\, \{f_x, f_y, 0\})$ the components are:

$$\mathbf{F} = \begin{pmatrix} f_x \\ f_y \\ p_x f_y - p_y f_x \end{pmatrix}
\begin{bmatrix} \mathbf{e}_{31} \\ \mathbf{e}_{32} \\ \mathbf{e}_{12} \end{bmatrix}$$

Two parts of $\mathbf{F}$ can be extracted afterwards if needed:

- **Weight** — the ideal (metric-free) part; the $\{\mathbf{e}_{31}, \mathbf{e}_{32}\}$
  components, encoding the line of action as an ideal line:
  ```cpp
  auto w = weight(F);   // BiVec2dp{fx, fy, 0} — ideal line (force encoded)
  ```

- **Bulk** — the metric-dependent part; the `e12` component, encoding the scalar torque
  about the origin:
  ```cpp
  auto tau = bulk(F);   // BiVec2dp{0, 0, px*fy - py*fx} — torque scalar
  ```

To recover the force as a direction vector, apply the **attitude function** `att()` directly
to $\mathbf{F}$. It returns the $\{\mathbf{e}_{31}, \mathbf{e}_{32}\}$ components as an
ideal vector:

```cpp
auto f_dir = att(F);   // Vec2dp{fx, fy, 0} — force direction vector
```

In 3D the same applies with six bivector components
$\{\mathbf{e}_{41}, \mathbf{e}_{42}, \mathbf{e}_{43}, \mathbf{e}_{23}, \mathbf{e}_{31}, \mathbf{e}_{12}\}$:
the first three encode the force and `att(F)` returns `Vec3dp{fx, fy, fz, 0}`; the last
three encode the torque and are returned by `bulk(F)`.

### World Frame → Body Frame

```cpp
auto const F_b = move2dp(F_w, rrev(M));   // 2D
auto const F_b = move3dp(F_w, rrev(M));   // 3D
```

### Why Force Application Points Are in the World Frame

The application points (`pivot_w`, `cm_w`) are naturally world-frame quantities:
`pivot_w` is stored as a fixed world coordinate; `cm_w = move2dp(O_b, M)` comes from
the motor. The wedge product `wdg(P, f)` requires the **actual geometric location** of
$P$ to produce the correct torque. Forces are therefore computed in the world frame and
the resulting bivector is transformed to body frame in one step:

```cpp
auto const F_up_w = wdg(pivot_w, vec2dp{0.0, m*g, 0.0});
auto const F_dn_w = wdg(cm_w,    vec2dp{0.0, -m*g, 0.0});
auto const F_b    = move2dp(F_up_w + F_dn_w, rrev(M));
```

This is equivalent to computing the forces directly in the body frame (using constant
body-frame application points and a rotated gravity direction), but simpler because
gravity direction is fixed in the world frame.

---

## Pivot Invariance Identity

For any finite point $P$ (with $P.z = 1$) and any scalar $\varphi$:

$$\text{move2dp}(P,\, \exp(\tfrac{1}{2}\,\varphi \cdot P)) = P$$

A point is invariant under rotation about itself. This is why encoding $B_b = \varphi \cdot
Q_b$ keeps the world image of $Q_b$ fixed:

$$\text{move2dp}(Q_b,\, M_0 \mathbin{⟇} \exp(\tfrac{1}{2}\,\varphi \cdot Q_b)) =
\text{move2dp}(Q_b,\, M_0) = \text{pivot}_w = \text{const}$$

The pivot stays at its initial world position for all $\varphi$ — no explicit constraint
enforcement is needed for the *position* of the pivot, only for the *velocity direction*
(the $B_b \propto Q_b$ constraint on the ODE state).

---

## Velocity: Body Frame vs. World Frame

The motor velocity bivector can be expressed in either frame:

$$\boldsymbol{\Omega}_b = \utilde{M}\, \frac{dM}{dt} \qquad \text{(body frame)}$$
$$\boldsymbol{\Omega}_w = \frac{dM}{dt}\, \utilde{M} \qquad \text{(world frame)}$$

They are related by the adjoint action of $M$:
$\boldsymbol{\Omega}_w = \text{move2dp}(\boldsymbol{\Omega}_b, M)$.

For the body-frame ODE, $\boldsymbol{\Omega}_b$ is the ODE state variable and satisfies:

$$\frac{d B_b}{dt} = \boldsymbol{\Omega}_b \quad \Longrightarrow \quad
\frac{dM}{dt} = M \mathbin{⟇} \boldsymbol{\Omega}_b$$

For pivot-constrained motion with pivot $Q_b$:
$\boldsymbol{\Omega}_b = \omega \cdot Q_b$, so $\omega = \boldsymbol{\Omega}_b.z$ (since
$Q_b.z = 1$).

---

## Body Frame vs. World Frame: Advantages

The body-frame formulation retains significant advantages as system complexity grows.

### Constant Inertia Tensor

The single most important advantage:

- **Body frame**: $\mathbf{I}$ is constant — invert once at initialization.
- **World frame**: $\mathbf{I}_w(t) = M\,\mathbf{I}_0\,\utilde{M}$ must be recomputed
  and re-inverted at every integration step ($O(n^3)$ per step in 3D).

For a single rigid body this overhead is manageable. For a system of $N$ coupled bodies
it becomes $N$ matrix inversions per step.

### Gyroscopic Terms

The Euler equations in the body frame:

$$\dot{\boldsymbol{\Omega}}_b = \mathbf{I}^{-1}
\bigl(\mathbf{F}_b - \boldsymbol{\Omega}_b \times \mathbf{I}\,\boldsymbol{\Omega}_b\bigr)$$

The gyroscopic term $\boldsymbol{\Omega}_b \times \mathbf{I}\,\boldsymbol{\Omega}_b$
(the `rcmt` term in the PGA code) is straightforward because $\mathbf{I}$ is constant.
In the world frame this term expands into time-derivatives of $\mathbf{I}_w(t)$, which
are messier to compute and accumulate numerical drift.

### Multi-Body Systems (Joints, Chains, Trees)

For articulated bodies (robot arms, vehicle suspensions, character rigs):

- Each body has its own body frame. Inter-body constraints (joints) are expressed as
  **constant** offsets in the respective body frames.
- In the world frame, every joint attachment point must be recomputed each step via
  the parent body's motor.
- The recursive Newton-Euler algorithm (RNEA) and articulated body algorithm (ABA) —
  the standard $O(N)$ algorithms for rigid body chains — are formulated entirely in body
  frames precisely because of this.

### Numerical Integration Stability

Body-frame ODEs integrate a **constant-coefficient** linear system locally (the
linearized dynamics around any state). World-frame ODEs have time-varying coefficients
tied to the current orientation — this can increase stiffness and worsen step-size
constraints for explicit integrators.

### When the World Frame Is Equally Good

- **2D single-body, no rotation**: $\mathbf{I}$ does not rotate → both frames identical.
- **Forces and constraints that are naturally world-frame** (gravity, ground contact,
  fixed pivots): computed in the world frame and transformed with one `move2dp`/`move3dp`
  call — cheap regardless of which frame the ODE lives in.

### Summary

| Criterion | Body frame | World frame |
| --- | --- | --- |
| Inertia tensor | Constant, invert once | Recomputed each step |
| Gyroscopic term | Simple (`rcmt`) | Requires $\dot{\mathbf{I}}_w$ |
| Multi-body chains | $O(N)$ algorithms native | Recompute all attachments |
| Force application | One transform per force | Natural |
| Integration stiffness | Lower (constant coefficients) | Higher (time-varying) |
| Code complexity (single body) | Slightly more setup | Slightly simpler setup |

For a single body with a fixed pivot the difference is small. As soon as a second body,
a joint, or free 3D rotation is introduced, the body-frame advantages compound quickly.
The PGA body-frame formulation $M(t) = M_0 \mathbin{⟇} \exp(\tfrac{1}{2}B_b(t))$ scales
to those cases without reformulation.

---

## 2D vs. 3D Reference

| Concept | PGA2DP | PGA3DP |
|---|---|---|
| Point | `Vec2dp{x, y, 1}` | `Vec3dp{x, y, z, 1}` |
| Direction | `Vec2dp{x, y, 0}` | `Vec3dp{x, y, z, 0}` |
| Line / force | `BiVec2dp{x, y, z}` | `BiVec3dp{vx,vy,vz,mx,my,mz}` |
| Motor type | `mvec2dp_u` | `mvec3dp_e` |
| Motor action | `move2dp(X, M)` | `move3dp(X, M)` |
| Reverse | `rrev(M)` | `rrev(M)` |
| Motor composition | `rgpr(M0, M1)` | `rgpr(M0, M1)` |
| Translation motor | `exp(0.5 * vec2dp{-ty, tx, 0})` | `exp(0.5 * bivec3dp{...})` |
| Rotation motor | `exp(0.5 * phi * Vec2dp point)` | `exp(0.5 * phi * BiVec3dp line)` |
| ODE velocity type | `Vec2dp` (`vec2dp`) | `BiVec3dp` (`bivec3dp`) |
| Pivot entity | `Vec2dp` point $Q_b$ | `BiVec3dp` line $L_b$ |
| Inertia matrix | `Inertia2dp` (3×3) | `Inertia3dp` (6×6) |
| Inertia about pivot | `get_plate_inertia(m,w,h, Q_b)` | `get_cuboid_inertia(m,w,h,d, L_b)` |

**Pivot entity grade matches velocity grade:** in PGA2DP the velocity $\Omega_b$ is a
`Vec2dp` (grade 1) and the pivot $Q_b$ is a point (grade 1); in PGA3DP the velocity is a
`BiVec3dp` (grade 2) and the pivot $L_b$ is a line (grade 2). The constraint
$\Omega_b \propto \text{pivot}$ is grade-consistent in both cases.

---

## Addendum: Lie Algebra vs. Lie Group — Independent Generators, Coupled Finite Displacements

### The Two Levels

The motion generator $\boldsymbol{\Omega}_b$ lives in the **Lie algebra** $\mathfrak{se}(2)$
(2D) or $\mathfrak{se}(3)$ (3D) — the tangent space to the identity of the motor group.
The motor $M$ itself lives in the **Lie group** $SE(2)$ or $SE(3)$. These two levels have
fundamentally different structure.

### PGA2DP: Independent at Velocity Level

The Lie algebra $\mathfrak{se}(2)$ is three-dimensional with basis generators:

| Generator | `vec2dp` encoding | Effect |
|---|---|---|
| Translation in $x$ | `{0, 1, 0}` | body drifts in $+x$ |
| Translation in $y$ | `{-1, 0, 0}` | body drifts in $+y$ |
| Rotation about origin | `{0, 0, 1}` | body spins about world origin |

A general instantaneous velocity $\boldsymbol{\Omega}_b = \{-v_y,\; v_x,\; \omega\}$
means **simultaneously**: translate at $(v_x, v_y)$ **and** rotate at $\omega$ — with
**no coupling** between the three components. Integrating velocity over an infinitesimal
step $dt$ gives:

$$B_b \leftarrow B_b + \boldsymbol{\Omega}_b \cdot dt$$

All three components of $B_b$ accumulate independently.

### The Coupling Enters Through the Exponential Map (Chasles' Theorem, 2D)

When computing a **finite** displacement via $M = \exp(\tfrac{1}{2} B_b)$, the structure
of the 2D plane imposes a hard constraint:

> **Chasles' theorem (2D):** Any rigid body displacement in 2D is either a pure
> rotation about some center, or a pure translation. There is no screw motion.

Consequently, for any $B_b$ with $B_b.z \neq 0$, the motor $\exp(\tfrac{1}{2} B_b)$
is **always a pure rotation** — the translational components merely determine the
rotation center:

$$\text{pivot} = \left(\frac{B_b.x}{B_b.z},\; \frac{B_b.y}{B_b.z}\right)$$

The "translation" encoded in $B_b.x$ and $B_b.y$ does not produce an independent
translational displacement; it shifts the rotation center away from the world origin.

**At the velocity (Lie algebra) level the components are independent.
At the displacement (Lie group) level they are not.**

### PGA3DP: True Screw Motion

In 3D, the Lie algebra $\mathfrak{se}(3)$ is six-dimensional: three for rotation axis
direction, three for translation velocity. The Chasles theorem for 3D states:

> Any rigid body displacement in 3D is a **screw motion**: rotation about an axis
> combined with a translation **along** that same axis (helical motion).

Exponentiating a combined 3D generator therefore yields a genuine screw — rotation and
translation are simultaneously present in the finite displacement. The rotation axis
direction and the translational pitch are independently adjustable because there exists a
direction (along the axis) that is geometrically orthogonal to the plane of rotation.

In 2D no such direction exists (the axis is perpendicular to the entire plane), so
simultaneous rotation and translation always collapse to rotation about a different center.

### Summary: 2D vs. 3D Behaviour of the Exponential Map

| | PGA2DP | PGA3DP |
|---|---|---|
| Lie algebra dimension | 3 | 6 |
| Velocity components | $(-v_y, v_x, \omega)$ — independent | $(v_x,v_y,v_z,\omega_x,\omega_y,\omega_z)$ — independent |
| $\exp(B)$ with rotation+translation | Pure rotation about shifted center | Screw motion (rotation + translation along axis) |
| Finite pivot / screw axis | Point $(B.x/B.z,\; B.y/B.z)$ | Line (Plücker coordinates of screw axis) |
| "Translation" in combined $B$ | Shifts rotation center | Contributes independent pitch along axis |

### Practical Consequence for the Accumulated Generator $B_b$

The displayed $B_b$ (and $B_w$) are accumulated Lie algebra elements — their components
are additive and independent **as generators**. The pivot decoded from $B_b$ as
$(B_b.x/B_b.z,\; B_b.y/B_b.z)$ is meaningful as the instantaneous rotation center of
the **finite** displacement $\exp(\tfrac{1}{2} B_b)$, not of the instantaneous velocity
$\boldsymbol{\Omega}_b$.

When $B_b = B_{b,0} + t \cdot \boldsymbol{\Omega}_b$ (uniform motion from initial
condition $B_{b,0}$):

- If $B_{b,0} = \mathbf{0}$: the decoded pivot equals the body-frame pivot of
  $\boldsymbol{\Omega}_b$ for all $t$, i.e. $(B_b.x/B_b.z, B_b.y/B_b.z) =
  (\Omega_b.x/\Omega_b.z, \Omega_b.y/\Omega_b.z) = Q_b$ — constant and
  geometrically meaningful.

- If $B_{b,0} \neq \mathbf{0}$ with $B_{b,0}.z = 0$ (a translational initial offset):
  the decoded pivot starts at infinity (pure translation at $t=0$) and converges
  asymptotically to $Q_b$ as $t \to \infty$. This makes the decoded pivot hard to
  interpret and should be avoided in visualizations intended to explain the encoding.

**Design rule for kinematic visualizations:** use $B_{b,0} = \mathbf{0}$ (body starts
at the $M_0$ reference pose) so that the accumulated generator retains a clean geometric
meaning throughout the animation.
