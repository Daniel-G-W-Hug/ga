# Structure of Inertia and Inertia Maps

## Overview

The inertia map $I[\Omega]$ is the linear map that takes the velocity twist $\Omega$ of a
rigid body and returns its momentum wrench $P$:

$$P = I[\Omega]$$

This is the PGA formulation of the combined Newton–Euler law $\dot{P} = F$, where the
momentum wrench $P$ unifies linear and angular momentum into a single algebraic object,
and the applied forque $F$ unifies force and torque similarly (see the preceding
subsections on modelling force and motion). The structure of the inertia map — where mass
and moment-of-inertia values appear in its matrix representation — is not arbitrary. It
is entirely determined by two things: the choice of basis for the domain and codomain, and
the Hodge complement structure of the respective algebra. Understanding this structure
makes the connection to classical mechanics transparent.

The derivation of the matrix elements (as a function of a single point mass) follows
directly from $I_{2D}[\Omega] = X \wedge m\,\mathrm{rcmt}(\Omega, X)$ in 2D and
$I_{3D}[\Omega] = X \wedge m\,\mathrm{rcmt}(\Omega, X)$ in 3D (see preceding subsection
on the inertia map). What is discussed here is how to read the resulting matrices and how
their entries relate to classical quantities.

---

## Part 1 — PGA2DP: 2D Inertia Map

### 1.1 Domain and Codomain

In PGA2DP the velocity twist $\Omega$ is a **vector** and the momentum wrench $P$ is a
**bivector**. The inertia map therefore acts between two different grades:

$$I_{2D}[\Omega] : \mathrm{Vec2dp} \to \mathrm{BiVec2dp}$$

The basis elements and their physical interpretation are:

| Index | Vec2dp — velocity twist $\Omega$ | BiVec2dp — momentum wrench $P$ |
|-------|----------------------------------|-------------------------------|
| 0 | `e1` — translational velocity, x-direction | `e31` — x-directed line, x-momentum |
| 1 | `e2` — translational velocity, y-direction | `e32` — y-directed line, y-momentum |
| 2 | `e3` — projective direction, angular velocity | `e12` — ideal line, angular momentum |

Note that `e3` is the projective basis vector (squaring to zero) and `e12` is the ideal
line (the line at infinity), which carries angular momentum in the PGA2DP encoding.

### 1.2 Single-Point Inertia Matrix

For a unitized point at position $(X_x, X_y)$ with $X_z = 1$ (projective coordinate)
and mass $m$, the inertia matrix is:

$$I_{2D} = m \begin{bmatrix}
0        &  X_z^2   & -X_y X_z      \\
-X_z^2   &  0       &  X_x X_z      \\
-X_x X_z & -X_y X_z & X_x^2 + X_y^2
\end{bmatrix}$$

Substituting $X_z = 1$ for a unitized point, and writing $r^2 = X_x^2 + X_y^2$:

$$I_{2D} = m \begin{bmatrix}
0    &  1   & -X_y      \\
-1   &  0   &  X_x      \\
-X_x & -X_y & r^2
\end{bmatrix}$$

where rows correspond to output components (e31, e32, e12) and columns to input
components ($\Omega_x$, $\Omega_y$, $\Omega_z$).

### 1.3 The Hodge Complement and Why Mass Is Off-Diagonal

The dominant feature of the matrix is that the mass $m$ appears at off-diagonal positions
$I[0,1]$ and $I[1,0]$, not on the diagonal. This is a direct consequence of the **Hodge
complement** structure of PGA2DP, not a coincidence or a sign of an unusual encoding.

The mass contribution to the inertia map is $I_\mathrm{mass}(\Omega) = -m \cdot
\mathbf{!}\Omega$, where $\mathbf{!}$ denotes the Hodge complement in PGA2DP. The
complement maps translational basis vectors as follows:

$$\mathbf{!}e_1 = e_{32}, \qquad \mathbf{!}e_2 = -e_{31}$$

Applying $I_\mathrm{mass} = -m \cdot \mathbf{!}$:

- $\Omega_x$ (component along `e1`, index 0) maps to $-m \cdot e_{32}$, which is at row
  index 1 → the entry $I[1,0] = -m$.
- $\Omega_y$ (component along `e2`, index 1) maps to $-m \cdot (-e_{31}) = +m \cdot
  e_{31}$, which is at row index 0 → the entry $I[0,1] = +m$.

The diagonal entries $I[0,0]$ and $I[1,1]$ are zero because the complement of `e1` is
`e32` (not `e31`), and the complement of `e2` is `-e31` (not `e32`): there is no
self-coupling. This is a consequence of the degenerate metric, where $e_3^2 = 0$.

The resulting mass block is **antisymmetric** with values $+m$ and $-m$. This
antisymmetry reflects the fact that the Hodge complement swaps the x- and y-directions
while introducing a sign change. It is the PGA algebraic encoding of the classical
momentum $\mathbf{p} = m\mathbf{v}$, correctly recovering $p_x = m v_x$ and $p_y = m
v_y$ once the kinematic encoding of the velocity twist is taken into account.

> **Note on kinematics**: The velocity twist $\Omega$ for a pure translation with
> classical velocity $(v_x, v_y)$ is encoded in PGA2DP as $\Omega = (-v_y, v_x, 0)$,
> i.e.\ the x- and y-components are swapped and one changes sign. This is a consequence
> of the motor-exponential encoding (translation as two reflections on parallel lines).
> With this encoding the two sign flips in the mass block cancel: $I[1,0] \cdot
> \Omega_x = (-m)(-v_y) = +m v_y$, correctly recovering $p_y = m v_y$.

### 1.4 Moment of Inertia

The entry $I[2,2] = m r^2$ maps angular velocity (index 2, `e3`) to angular momentum
(index 2, `e12`). Both the input and output carry the same rotational index, so this is
the diagonal of the rotational sub-block. It recovers the classical scalar moment of
inertia $J = m r^2$ for a point mass at distance $r$ from the rotation centre.

The Hodge complement of `e3` is $\mathbf{!}e_3 = -e_{12}$ (the ideal line), connecting
the angular velocity generator to the angular momentum carrier. The negative sign is
absorbed into the definition so that $J > 0$.

### 1.5 Cross-Coupling: Angular Velocity and Linear Momentum

The entries $I[0,2] = -m X_y$ and $I[1,2] = +m X_x$ couple the angular velocity
component to the linear momentum components. They encode the classical relation

$$\mathbf{L}_\mathrm{linear} = m(\boldsymbol{\omega} \times \mathbf{r})$$

i.e.\ a rotating body at offset $(X_x, X_y)$ from the rotation centre carries linear
momentum in addition to angular momentum. These entries vanish when the point is at the
origin ($X_x = X_y = 0$). By the transpose symmetry of the inertia map (which holds in
the body frame), the entries $I[2,0] = -m X_x$ and $I[2,1] = -m X_y$ encode the
conjugate coupling: translational velocity of an off-axis point contributes to angular
momentum.

### 1.6 Inertia Matrix in the Body Frame (CoM at Origin)

Integrating the single-point inertia matrix $I_{2D}(m, X)$ over a rigid body with centre
of mass at the origin (body frame) gives:

$$I_\mathrm{body} = \begin{bmatrix}
0  &  m  &  0 \\
-m &  0  &  0 \\
0  &  0  &  J
\end{bmatrix}$$

where $J = \int (x^2 + y^2)\,dm$ is the classical scalar moment of inertia about the
centre of mass. The mass block (upper-left 2×2) and the moment of inertia entry
(lower-right 1×1) are cleanly separated: translational and rotational degrees of freedom
decouple at the centre of mass, exactly as in classical mechanics.

### 1.7 General Frame: Centre of Mass Not at Origin

For a body with centre of mass at $(c_x, c_y)$ relative to the chosen reference point,
the full integrated inertia matrix is:

$$I = \begin{bmatrix}
0     &  m     & -m c_y         \\
-m    &  0     &  m c_x         \\
-m c_x & -m c_y & J_\mathrm{cm} + m(c_x^2 + c_y^2)
\end{bmatrix}$$

where $J_\mathrm{cm}$ is the classical scalar moment of inertia about the centre of mass.
The matrix has the following block structure:

| Sub-block | Entries | Physical meaning |
|-----------|---------|-----------------|
| Upper-left 2×2 | $[0, m;\;{-m}, 0]$ | Pure antisymmetric mass block — **frame-invariant**, unchanged regardless of CoM position |
| Upper-right 2×1 | $[-m c_y;\; m c_x]$ | Angular velocity → linear momentum coupling. Zero only when CoM is at origin. Encodes $\mathbf{L}_\mathrm{linear} = m(\boldsymbol{\omega} \times \mathbf{r}_\mathrm{cm})$ |
| Lower-left 1×2 | $[-m c_x,\; -m c_y]$ | Translational velocity → angular momentum coupling. Transpose of upper-right (up to sign convention) |
| Lower-right 1×1 | $J_\mathrm{cm} + m r^2$ | Classical parallel axis theorem: $J_\mathrm{cm}$ increased by $m(c_x^2 + c_y^2)$ |

The CoM offset only populates the cross-coupling blocks. The mass block $[0, m; -m, 0]$
is rigid and frame-invariant, a direct consequence of mass being a Hodge-dual coupling
that does not depend on position.

### 1.8 Continuous Limit: Uniform Rectangular Plate

As a reference example, consider a uniform rectangular plate of width $W$ (along `e1`),
height $H$ (along `e2`), total mass $M$, centred at the origin. Integrating the
single-point formula over the plate gives:

$$I_\mathrm{plate} = M \begin{bmatrix}
0  &  1  &  0 \\
-1 &  0  &  0 \\
0  &  0  & \frac{W^2 + H^2}{12}
\end{bmatrix}$$

The entry $I[2,2] = M(W^2 + H^2)/12$ is the standard formula for the second moment of
area of a rectangle. The cross-coupling terms vanish because the centre of mass is at the
origin, and the off-diagonal mass terms take their canonical body-frame values $\pm M$.

For a discrete grid of $n \times n$ equally spaced point masses (endpoint-inclusive), the
relative error in $I[2,2]$ compared to the continuous limit is $2/(n-1)$, independent of
the plate dimensions. At least $n = 202$ points per dimension are needed to keep the
error below $1\%$; the discrete sum always converges from above as $n \to \infty$.

### 1.9 Classical Quantities and PGA2DP Entries

| Classical quantity | PGA2DP entry | Condition |
|--------------------|-------------|-----------|
| Total mass $m$ | $I[0,1] = +m$, $I[1,0] = -m$ | Always |
| Scalar moment of inertia $J$ about CoM | $I[2,2] = J$ | CoM at origin |
| CoM offset $(c_x, c_y)$ | $I[0,2] = -m c_y$, $I[1,2] = +m c_x$ (and transposed) | General frame |
| Parallel axis correction | $I[2,2] = J_\mathrm{cm} + m(c_x^2 + c_y^2)$ | General frame |

---

## Part 2 — PGA3DP: 3D Inertia Map

### 2.1 Domain and Codomain

In PGA3DP the velocity twist $\Omega$ is a **bivector** and the momentum wrench $P$ is
also a **bivector**. In contrast to the 2D case, the inertia map acts within the same
grade:

$$I_{3D}[\Omega] : \mathrm{BiVec3dp} \to \mathrm{BiVec3dp}$$

The six basis bivectors split naturally into two groups, corresponding to translational
and rotational degrees of freedom:

| Index | BiVec3dp — velocity twist $\Omega$ | BiVec3dp — momentum wrench $P$ |
|-------|-----------------------------------|-----------------------------|
| 0 | `e41` — ideal line, translational velocity x | `e41` — ideal line, linear momentum x |
| 1 | `e42` — ideal line, translational velocity y | `e42` — ideal line, linear momentum y |
| 2 | `e43` — ideal line, translational velocity z | `e43` — ideal line, linear momentum z |
| 3 | `e23` — real line, angular velocity about `e1` | `e23` — real line, angular momentum about `e1` |
| 4 | `e31` — real line, angular velocity about `e2` | `e31` — real line, angular momentum about `e2` |
| 5 | `e12` — real line, angular velocity about `e3` | `e12` — real line, angular momentum about `e3` |

The **ideal lines** (`e41`, `e42`, `e43`) are the bivectors containing the projective
direction `e4` (which squares to zero). They carry translational information. The **real
lines** (`e23`, `e31`, `e12`) are the bivectors of the Euclidean subspace and carry
rotational information.

### 2.2 Single-Point Inertia Matrix

For a unitized point at position $(X_x, X_y, X_z)$ with $X_w = 1$ and mass $m$, the
$6 \times 6$ inertia matrix is:

$$I_{3D} = m \begin{bmatrix}
0              &  X_z X_w       & -X_y X_w       & X_w^2    &  0       &  0       \\
-X_z X_w       &  0             &  X_x X_w       & 0        &  X_w^2   &  0       \\
 X_y X_w       & -X_x X_w       &  0             & 0        &  0       &  X_w^2   \\
 X_y^2 + X_z^2 & -X_x X_y       & -X_x X_z       & 0        & -X_z X_w &  X_y X_w \\
-X_x X_y       &  X_x^2 + X_z^2 & -X_y X_z       &  X_z X_w &  0       & -X_x X_w \\
-X_x X_z       & -X_y X_z       &  X_x^2 + X_y^2 & -X_y X_w &  X_x X_w &  0
\end{bmatrix}$$

Substituting $X_w = 1$:

$$I_{3D} = m \begin{bmatrix}
0              &  X_z       & -X_y       & 1    &  0   &  0   \\
-X_z           &  0         &  X_x       & 0    &  1   &  0   \\
 X_y           & -X_x       &  0         & 0    &  0   &  1   \\
 X_y^2 + X_z^2 & -X_x X_y   & -X_x X_z   & 0    & -X_z &  X_y \\
-X_x X_y       &  X_x^2+X_z^2 & -X_y X_z &  X_z &  0   & -X_x \\
-X_x X_z       & -X_y X_z   &  X_x^2+X_y^2 & -X_y &  X_x &  0
\end{bmatrix}$$

where rows correspond to output components (`e41, e42, e43, e23, e31, e12`) and columns
to input components ($\Omega_{vx}, \Omega_{vy}, \Omega_{vz}, \Omega_{mx}, \Omega_{my},
\Omega_{mz}$).

### 2.3 Block Structure and the Hodge Complement in 3D

The $6 \times 6$ matrix splits into four $3 \times 3$ blocks. For the single-point matrix
with $X_w = 1$ the block structure at an arbitrary position is:

```
         [ vx  vy  vz | mx  my  mz ]     ← Ω input
(e41) vx [  0  Xz -Xy |  1   0   0 ]
(e42) vy [ -Xz  0  Xx |  0   1   0 ]   upper-right = I₃  (mass, Hodge-dual coupling)
(e43) vz [  Xy -Xx  0 |  0   0   1 ]
         [------------|-------------]
(e23) mx [ Xy²+Xz²  -XxXy  -XxXz |  0  -Xz  Xy ]
(e31) my [ -XxXy  Xx²+Xz² -XyXz  |  Xz   0  -Xx ]   lower-left = J  (moments)
(e12) mz [ -XxXz  -XyXz  Xx²+Xy² | -Xy  Xx   0  ]
```

The four blocks play distinct physical roles:

- **Upper-left [0:3, 0:3]**: translational velocity → linear momentum. For a point at the
  origin this block is zero (no mass self-coupling in the translational subspace due to
  the degenerate metric). For a point at $(X_x, X_y, X_z)$ the upper-left becomes an
  antisymmetric matrix encoding the coupling $\mathbf{p} = m(X_w \mathbf{v} +
  \boldsymbol{\omega} \times \mathbf{X})$ when CoM is not at origin.
- **Upper-right [0:3, 3:6]**: angular velocity → linear momentum. At $X_w = 1$ this
  block equals the identity $I_3$, scaled by $m$. This is the 3D analogue of the
  off-diagonal mass entries in the 2D case: mass couples angular velocity generators
  (real lines) to linear momentum carriers (ideal lines) via the Hodge complement.
- **Lower-left [3:6, 0:3]**: translational velocity → angular momentum. This block
  carries the classical moments and products of inertia — the PGA form of the 3D inertia
  tensor.
- **Lower-right [3:6, 3:6]**: angular velocity → angular momentum. For a point at the
  origin this block is zero. For a displaced point it becomes antisymmetric, encoding the
  parallel axis contribution to angular momentum.

The placement of mass in the **upper-right** block (connecting real-line velocity to
ideal-line momentum) follows from the Hodge complement in PGA3DP:

$$\mathbf{!}e_{23} = e_{41}, \quad \mathbf{!}e_{31} = e_{42}, \quad \mathbf{!}e_{12} =
e_{43}$$

The complement maps each real line (rotational generator) to the corresponding ideal line
(linear momentum carrier). The mass contribution $I_\mathrm{mass}(\Omega) = -m \cdot
\mathbf{!}\Omega$ therefore places $m$ on the diagonal of the upper-right block. Diagonal
entries in the upper-right mean that each angular velocity component $\Omega_{mx}$,
$\Omega_{my}$, $\Omega_{mz}$ couples directly to the corresponding linear momentum
component.

### 2.4 Inertia Matrix in the Body Frame (CoM at Origin)

Integrating the single-point inertia matrix over a rigid body with centre of mass at the
origin gives the body-frame inertia matrix:

$$I_\mathrm{body} = \begin{bmatrix}
0    & 0    & 0    &  m   & 0    & 0   \\
0    & 0    & 0    &  0   & m    & 0   \\
0    & 0    & 0    &  0   & 0    & m   \\
I_{xx} & -I_{xy} & -I_{xz} & 0 & 0 & 0 \\
-I_{xy} & I_{yy} & -I_{yz} & 0 & 0 & 0 \\
-I_{xz} & -I_{yz} & I_{zz} & 0 & 0 & 0
\end{bmatrix}$$

written schematically as:

$$I_\mathrm{body} = \begin{bmatrix}
\mathbf{0} & m\, I_3 \\
J_\mathrm{cm} & \mathbf{0}
\end{bmatrix}$$

where $\mathbf{0}$ denotes the $3 \times 3$ zero block, $m\,I_3$ is $m$ times the
identity matrix, and $J_\mathrm{cm}$ is the classical $3 \times 3$ inertia tensor about
the centre of mass:

$$J_\mathrm{cm} = \begin{bmatrix}
I_{xx} & -I_{xy} & -I_{xz} \\
-I_{xy} & I_{yy}  & -I_{yz} \\
-I_{xz} & -I_{yz} & I_{zz}
\end{bmatrix}$$

with the classical integrals:

$$I_{xx} = \int(y^2+z^2)\,dm, \quad I_{yy} = \int(x^2+z^2)\,dm, \quad I_{zz} =
\int(x^2+y^2)\,dm$$
$$I_{xy} = \int xy\,dm, \quad I_{xz} = \int xz\,dm, \quad I_{yz} = \int yz\,dm$$

Again, translational and rotational degrees of freedom are cleanly separated at the
centre of mass: the diagonal blocks are zero and all physical content sits in the
off-diagonal blocks.

### 2.5 General Frame: Centre of Mass Not at Origin

For a body with centre of mass at $\mathbf{r} = (c_x, c_y, c_z)$ the full $6 \times 6$
inertia matrix is:

$$I = \begin{bmatrix}
0       &  m c_z   & -m c_y   &  m   &  0   &  0  \\
-m c_z  &  0       &  m c_x   &  0   &  m   &  0  \\
 m c_y  & -m c_x   &  0       &  0   &  0   &  m  \\
I_{xx}  & -I_{xy}  & -I_{xz}  &  0   & -m c_z &  m c_y \\
-I_{xy} &  I_{yy}  & -I_{yz}  &  m c_z &  0   & -m c_x \\
-I_{xz} & -I_{yz}  &  I_{zz}  & -m c_y &  m c_x &  0
\end{bmatrix}$$

where the lower-left inertia tensor entries already include the parallel axis corrections:

$$I_{xx} = I_{xx,\mathrm{cm}} + m(c_y^2+c_z^2), \quad I_{xy} = I_{xy,\mathrm{cm}} + m
c_x c_y, \quad \text{etc.}$$

The block structure in the general frame is:

| Block | Body frame (CoM at origin) | General frame | Physical meaning |
|-------|--------------------------|---------------|-----------------|
| Upper-left [0:3, 0:3] | $\mathbf{0}$ | $m[\mathbf{r}\times]^\top$ | Trans.\ vel.\ → linear mom.\ (CoM coupling) |
| Upper-right [0:3, 3:6] | $m\,I_3$ | $m\,I_3$ | Angular vel.\ → linear mom.\ (**always** $m\,I_3$, **frame-invariant**) |
| Lower-left [3:6, 0:3] | $J_\mathrm{cm}$ | $J_\mathrm{cm}$ + parallel axis | Trans.\ vel.\ → angular mom. |
| Lower-right [3:6, 3:6] | $\mathbf{0}$ | $m[\mathbf{r}\times]$ | Angular vel.\ → angular mom.\ (CoM coupling) |

Here $[\mathbf{r}\times]$ is the $3\times 3$ antisymmetric cross-product matrix of the
CoM position vector:

$$[\mathbf{r}\times] = \begin{bmatrix}
0    & -c_z &  c_y \\
c_z  &  0   & -c_x \\
-c_y &  c_x &  0
\end{bmatrix}$$

Two important observations carry over from the 2D case:

1. **The upper-right block $m\,I_3$ is frame-invariant.** It equals $m\,I_3$ regardless
   of the CoM position. Mass is a Hodge-dual coupling between ideal and real lines that
   does not depend on position — this is the 3D counterpart of the invariant mass block
   in the 2D case.
2. **Upper-left and lower-right are conjugate:** upper-left $= -\text{lower-right}^\top$.
   Both are antisymmetric and driven by the same CoM offset $\mathbf{r}$, reflecting the
   conjugate nature of the translational-to-translational and angular-to-angular couplings
   induced by the offset.

This $6\times 6$ structure is the PGA form of the classical **spatial inertia matrix**
from rigid body dynamics (cf.\ Featherstone), expressed in the (translational, angular)
velocity ordering of PGA3DP.

### 2.6 Continuous Limit: Uniform Rectangular Cuboid

As a reference example, consider a uniform rectangular cuboid of width $w$ (along `e1`),
height $h$ (along `e2`), depth $d$ (along `e3`), total mass $m$, centred at the origin.
All cross-inertia and centre-of-mass coupling terms vanish, giving:

$$I_\mathrm{cuboid} = m \begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
\frac{h^2+d^2}{12} & 0 & 0 & 0 & 0 & 0 \\
0 & \frac{w^2+d^2}{12} & 0 & 0 & 0 & 0 \\
0 & 0 & \frac{w^2+h^2}{12} & 0 & 0 & 0
\end{bmatrix}$$

The three diagonal entries of the lower-left block are the classical rectangular moments
of inertia:

| Entry | Value | Rotation axis | Perpendicular extents |
|-------|-------|---------------|-----------------------|
| $I[3,0]$ | $m(h^2+d^2)/12$ | `e1` (x-axis) | $h$ and $d$ |
| $I[4,1]$ | $m(w^2+d^2)/12$ | `e2` (y-axis) | $w$ and $d$ |
| $I[5,2]$ | $m(w^2+h^2)/12$ | `e3` (z-axis) | $w$ and $h$ |

### 2.7 Classical Quantities and PGA3DP Entries

In the body frame (CoM at origin):

| Classical quantity | PGA3DP block | Entries |
|-------------------|-------------|---------|
| Total mass $m$ | Upper-right diagonal ($m\,I_3$) | $I[0,3]=I[1,4]=I[2,5]=m$ |
| Diagonal moments $I_{xx}, I_{yy}, I_{zz}$ | Lower-left diagonal | $I[3,0], I[4,1], I[5,2]$ |
| Products of inertia $I_{xy}, I_{xz}, I_{yz}$ | Lower-left off-diagonal | $I[3,1]=I[4,0]=-I_{xy}$, $I[3,2]=I[5,0]=-I_{xz}$, $I[4,2]=I[5,1]=-I_{yz}$ |
| CoM offset (general frame) | All four $3\times 3$ blocks | via parallel axis theorem |

In the **principal-axis body frame** all products of inertia vanish ($I_{xy} = I_{xz} =
I_{yz} = 0$), leaving only the three diagonal moments — exactly the cuboid formula above.

---

## Part 3 — Unified Picture: 2D and 3D

### 3.1 Common Algebraic Origin: the Hodge Complement

In both PGA2DP and PGA3DP the mass contribution to the inertia map takes the same
algebraic form:

$$I_\mathrm{mass}(\Omega) = -m \cdot \mathbf{!}\Omega$$

The Hodge complement $\mathbf{!}$ maps each velocity generator to its momentum carrier
in the dual subspace. In PGA the translational and rotational generators live in
Hodge-dual subspaces (ideal lines ↔ real lines), so the mass term always places $m$ in
the off-diagonal block connecting these two subspaces. Where exactly that block appears in
the matrix depends on the basis ordering.

### 3.2 Basis Ordering Determines Block Position

The difference in where mass appears (upper-left 2×2 in 2D vs upper-right 3×3 in 3D) is
entirely due to the basis ordering convention:

| Algebra | Translational indices | Rotational indices | Mass block position |
|---------|-----------------------|--------------------|---------------------|
| PGA2DP | 0, 1 (lower) | 2 (higher) | **Upper-left** 2×2: trans.\ vel.\ → lin.\ mom. |
| PGA3DP | 0, 1, 2 (lower) | 3, 4, 5 (higher) | **Upper-right** 3×3: angular vel.\ → lin.\ mom. |

In PGA2DP the translational generators (`e1`, `e2`) sit at lower indices and the
rotational generator (`e3`) at the highest index. The Hodge complement maps `e1` → `e32`
and `e2` → `−e31`, which are the momentum carriers in the **translational** output
rows — hence mass lands in the upper-left sub-block.

In PGA3DP the translational generators (`e41`, `e42`, `e43`) are ideal lines and sit at
lower indices (0–2), while the rotational generators (`e23`, `e31`, `e12`) are real lines
at higher indices (3–5). The Hodge complement maps `e23` → `e41`, `e31` → `e42`, `e12`
→ `e43`, which are the momentum carriers in the **translational** output rows (indices
0–2) — hence mass lands in the upper-right sub-block (column range 3–5, row range 0–2).

### 3.3 Frame Invariance of the Mass Block

In both dimensions the mass block is **frame-invariant**: it retains the same values
regardless of where the centre of mass is located. The mass term is a Hodge-dual coupling
between velocity generators and momentum carriers that does not depend on position. This
is why the mass block does not participate in the parallel axis theorem — only the moment
block and the coupling blocks change when the reference frame is shifted.

### 3.4 Unified Summary Table

| Classical quantity | PGA2DP ($3\times 3$) | PGA3DP ($6\times 6$) |
|-------------------|---------------------|---------------------|
| Total mass $m$ | Upper-left 2×2: $I[0,1]=+m$, $I[1,0]=-m$ | Upper-right 3×3 diagonal: $I[0,3]=I[1,4]=I[2,5]=m$ |
| Scalar/diagonal moments | Lower-right 1×1: $I[2,2]=J$ | Lower-left 3×3 diagonal: $I[3,0], I[4,1], I[5,2]$ |
| Products of inertia | n/a (scalar $J$ in 2D) | Lower-left 3×3 off-diagonal |
| CoM coupling | $I[0,2], I[1,2]$ and transposed | All four 3×3 blocks non-zero |
| Parallel axis | $I[2,2] = J_\mathrm{cm} + m r^2$ | Lower-left 3×3 (standard 3D parallel axis) |

### 3.5 Decoupling at the Centre of Mass

In both 2D and 3D, choosing the body frame with origin at the centre of mass ($c_x = c_y
= 0$ in 2D, $c_x = c_y = c_z = 0$ in 3D) leads to a clean decoupling:

- The mass block retains its invariant form.
- All cross-coupling blocks (translational velocity to angular momentum and vice versa)
  become zero.
- The moment-of-inertia block takes its simplest form (the classical tensor $J_\mathrm{cm}$
  without any parallel axis corrections).

This matches the classical result that Newton's and Euler's equations decouple at the
centre of mass, and it is the reason the body frame with CoM at the origin is the natural
choice for integrating the equations of motion (see preceding subsection on solving the
Newton–Euler equations).
