# Coupled Newton-Euler Equations for Rigid Body Dynamics

### Why Coupling Arises

If you express everything in an **inertial frame**, translation and rotation *seem*
decoupled:

$$m\ddot{\mathbf{r}} = \mathbf{F}, \qquad \dot{\mathbf{L}} = \boldsymbol{\tau}$$

But the moment you work in a **body-fixed frame** (which is almost always necessary for
non-spherical bodies, since the inertia tensor $\mathbf{I}$ is then constant), the
rotation of the frame itself introduces coupling terms.

---

### The Full Coupled System

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

### Solving for the Derivatives (state-space form)

$$\dot{\mathbf{v}} = \frac{\mathbf{F}^b}{m} - \boldsymbol{\omega} \times \mathbf{v}$$

$$\dot{\boldsymbol{\omega}} = \mathbf{I}^{-1}\left(\boldsymbol{\tau}^b -
\boldsymbol{\omega} \times \mathbf{I}\boldsymbol{\omega}\right)$$

The state vector is $[\mathbf{v},\, \boldsymbol{\omega}]^T \in \mathbb{R}^6$.

---

### Also Integrate the Orientation

To recover position/orientation in the inertial frame:

$$\dot{\mathbf{r}} = \mathbf{R}\,\mathbf{v}$$

$$\dot{\mathbf{R}} = \mathbf{R}\,[\boldsymbol{\omega}]_\times \quad \text{(or use
quaternions: } \dot{\mathbf{q}} = \tfrac{1}{2}\mathbf{q} \otimes
\boldsymbol{\omega}\text{)}$$

**Quaternions are strongly recommended** to avoid gimbal lock.

---

### Complete State Vector for Integration

| State | Size | Equation |
|---|---|---|
| $\mathbf{r}$ | 3 | $\dot{\mathbf{r}} = \mathbf{R}\mathbf{v}$ |
| $\mathbf{q}$ (quaternion) | 4 | $\dot{\mathbf{q}} = \frac{1}{2}\mathbf{q}\otimes\boldsymbol{\omega}$ |
| $\mathbf{v}$ (body frame) | 3 | $\dot{\mathbf{v}} = \mathbf{F}^b/m - \boldsymbol{\omega}\times\mathbf{v}$ |
| $\boldsymbol{\omega}$ (body frame) | 3 | $\dot{\boldsymbol{\omega}} = \mathbf{I}^{-1}(\boldsymbol{\tau}^b - \boldsymbol{\omega}\times\mathbf{I}\boldsymbol{\omega})$ |

**Total: 13 coupled 1st-order ODEs** (with quaternion normalization constraint
$|\mathbf{q}|=1$).

---

### Key Differences from a "Parallel" Approach

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
