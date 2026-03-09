# Inertia Map Analysis

## PGA Inertia Matrix Structure — Why mass is off-diagonal

### Core insight: the map crosses two dual spaces

Both the PGA2DP and PGA3DP inertia maps connect **two geometrically dual spaces**, not
the same space like the classical inertia tensor (which maps ω → L within one vector
space). The apparent paradox of mass appearing off-diagonal is entirely a consequence of
the Hodge complement structure of PGA.

---

## PGA2DP: 3×3 Inertia Map  `I[Ω]: Vec2dp → BiVec2dp`

### Basis correspondence (Vec2dp)

| Index | Vec2dp input (velocity twist Omega) | BiVec2dp output (momentum wrench L) |
| --- | --- | --- |
| 0 | `e1` - x-component of velocity twist | `e31` - x-directed line, x-momentum |
| 1 | `e2` - y-component of velocity twist | `e32` - y-directed line, y-momentum |
| 2 | `e3` - projective component, angular velocity | `e12` - ideal line, angular momentum |

### Point mass matrix (unitized point at (Xx, Xy), Xz = 1)

```text
I = m * [  0    1   -Xy ]     row 0 → e31 (x-momentum)
        [ -1    0    Xx ]     row 1 → e32 (y-momentum)
        [-Xx  -Xy    r² ]     row 2 → e12 (angular momentum),  r² = Xx² + Xy²
```

### Why mass appears off-diagonal: the Hodge complement crossing

The mass terms of the inertia map follow `I_mass(Ω) = −m · !Ω`, where `!` is the Hodge
complement in PGA2DP. The complement of the translational basis elements crosses indices:

```text
!e1 = e32    (x-component maps to the y-momentum line)
!e2 = −e31   (y-component maps to the −x-momentum line)
```

Applying `I_mass = −m · !`:

- `Ω.x` (e1, index 0) → `I_mass(e1) = −m · !e1 = −m · e32` → lands at `I[1,0] = −m`
- `Ω.y` (e2, index 1) → `I_mass(e2) = −m · !e2 = −m · (−e31) = +m · e31` → lands at `I[0,1] = +m`

The Hodge complement swaps x↔y and adds a sign, which is why mass appears **off-diagonal**
and **antisymmetric** (`+m` / `−m`).

`I[0,0] = I[1,1] = 0`: the complement of e1 is e32 (not e31), and of e2 is −e31 (not e32),
so there is no self-coupling on the diagonal for the mass terms.

**Note on kinematics**: The above is a direct consequence of the Hodge structure. The
connection to classical Newton's law `p = mv` requires the additional kinematic fact that
in PGA2DP the motor-generator (velocity twist Omega) for a translation with classical
velocity `(vx, vy)` is `Omega = att(bulk_dual(v)) = (−vy, vx, 0)`, not `(vx, vy, 0)`.
This is specific to the motor-exponential encoding (translation as two reflections on
parallel lines) and does **not** mean that e1 generally encodes y-velocity: in any normal
Vec2dp (position, force, acceleration), e1 is always the x-component and e2 the y-component.
With the twisted Omega encoding, the two sign flips cancel: `I[1,0] · Omega.x = (−m)(−vy) = +m·vy`,
recovering `p_y = m·vy`.

### Why the moment of inertia is at (2,2)

`I[2,2] = m·r²` maps angular velocity (index 2, `e3`, origin point) to angular momentum
(index 2, `e12`, ideal line). Both carry the same rotational index → this IS the
diagonal of the rotational sub-block, recovering the classical scalar moment of inertia.

### Cross-coupling terms I[0,2], I[1,2]

`I[0,2] = -m·Xy`, `I[1,2] = m·Xx`: angular velocity of a body offset from the rotation
axis carries linear momentum (PGA encoding of `L_linear = m·(ω × r)`). Vanish when
center of mass is at the origin.

### Continuous limit: uniform rectangular plate

For a plate of width W (e1), height H (e2), mass M, centered at origin:

```text
I_plate = M * [  0    1         0        ]
              [ -1    0         0        ]
              [  0    0    (W²+H²)/12    ]
```

`I[2,2] = M·(W²+H²)/12` is the standard rectangle formula.

### Discrete approximation (endpoint-inclusive n×n grid)

Relative error in `I[2,2]` for an n-point endpoint-inclusive grid: **2/(n−1)**, independent
of plate dimensions.

- n = 25:  error ≈ 8.3%  (I[2,2] = 65/36 ≈ 1.806 vs. continuous 5/3 ≈ 1.667)
- n = 202: error < 1%
- n → ∞:  converges from above to M·(W²+H²)/12

### PGA2DP entry summary

| Entry | Value | Physical meaning | Why there |
| --- | --- | --- | --- |
| `I[0,1]` | `+m` | Omega.y (e2) -> e31 (x-momentum) | Hodge: !e2=-e31, I_mass=-m*! gives +m |
| `I[1,0]` | `-m` | Omega.x (e1) -> e32 (y-momentum) | Hodge: !e1=+e32, I_mass=-m*! gives -m |
| `I[0,0]`, `I[1,1]` | `0` | no self-coupling | degenerate metric e3^2=0 |
| `I[0,2]`, `I[1,2]` | `~ position` | angular vel -> linear momentum (offset) | omega x r coupling |
| `I[2,0]`, `I[2,1]` | `~ position` | translational vel -> angular momentum | same coupling, transposed |
| `I[2,2]` | `m*r^2` | angular vel -> angular momentum | rotational diagonal, classical m*r^2 |

---

## PGA3DP: 6×6 Inertia Map  `I[Ω]: BiVec3dp → BiVec3dp`

### Basis correspondence (BiVec3dp)

| Index | BiVec3dp input (velocity twist Omega) | BiVec3dp output (momentum wrench L) |
| --- | --- | --- |
| 0 | `e41` (vx) - ideal line, translational velocity x | `e41` (vx) - linear momentum x |
| 1 | `e42` (vy) - ideal line, translational velocity y | `e42` (vy) - linear momentum y |
| 2 | `e43` (vz) - ideal line, translational velocity z | `e43` (vz) - linear momentum z |
| 3 | `e23` (mx) - real line, angular velocity about e1 | `e23` (mx) - angular momentum about e1 |
| 4 | `e31` (my) - real line, angular velocity about e2 | `e31` (my) - angular momentum about e2 |
| 5 | `e12` (mz) - real line, angular velocity about e3 | `e12` (mz) - angular momentum about e3 |

### Point mass matrix (unitized point at (Xx, Xy, Xz), Xw = 1)

```text
I = m *
  row 0: [  0          Xz·Xw      -Xy·Xw       Xw²    0        0     ]
  row 1: [ -Xz·Xw       0          Xx·Xw        0      Xw²      0     ]
  row 2: [  Xy·Xw      -Xx·Xw       0           0      0        Xw²   ]
  row 3: [  Xy²+Xz²   -Xx·Xy      -Xx·Xz        0     -Xz·Xw    Xy·Xw ]
  row 4: [ -Xx·Xy       Xx²+Xz²   -Xy·Xz        Xz·Xw  0       -Xx·Xw ]
  row 5: [ -Xx·Xz      -Xy·Xz       Xx²+Xy²    -Xy·Xw  Xx·Xw    0     ]
```

### Block structure and the Hodge crossing in 3D

The 6×6 matrix splits into four 3×3 blocks:

```text
        [ vx  vy  vz ]   [ mx  my  mz ]       ← Ω input
   vx   [             |               ]
   vy   [    0        |    m·I₃       ]   upper-right: angular vel → linear momentum
   vz   [             |               ]
        [ ------------|-------------- ]
   mx   [             |               ]
   my   [    m·J_rot  |    0          ]   lower-left: translational vel → angular momentum
   mz   [             |               ]
```

- **Upper-left [0:3, 0:3] = 0**: no translational-to-translational coupling
- **Upper-right [0:3, 3:6] = m·I₃**: angular velocity (e23/e31/e12, real lines)
  → linear momentum (e41/e42/e43, ideal lines) via mass m
- **Lower-left [3:6, 0:3] = m·J_rot**: translational velocity (e41/e42/e43)
  → angular momentum (e23/e31/e12) via classical moments of inertia
- **Lower-right [3:6, 3:6] = 0**: no rotational-to-rotational coupling

The Hodge complement in PGA3DP swaps ideal lines (e41,e42,e43) ↔ real lines
(e23,e31,e12). This places mass in the upper-right block (angular→linear) and moments of
inertia in the lower-left block (translational→angular) — the exact 3D analogue of mass
appearing off-diagonal in the 2D case.

### Continuous limit: uniform rectangular cuboid

For a cuboid of width w (e1), height h (e2), depth d (e3), mass m, centered at origin,
all cross-inertia and center-of-mass terms vanish:

```text
I_cuboid = m *
  [  0              0              0             1    0    0  ]
  [  0              0              0             0    1    0  ]
  [  0              0              0             0    0    1  ]
  [ (h²+d²)/12      0              0             0    0    0  ]
  [  0             (w²+d²)/12      0             0    0    0  ]
  [  0              0             (w²+h²)/12     0    0    0  ]
```

J_rot diagonal entries are the classical rectangle-rule moments of inertia:

| Entry | Value | Rotation axis | Perpendicular extents |
| --- | --- | --- | --- |
| `I[3,0]` | `m*(h^2+d^2)/12` | e1 (x-axis) | h (e2) and d (e3) |
| `I[4,1]` | `m*(w^2+d^2)/12` | e2 (y-axis) | w (e1) and d (e3) |
| `I[5,2]` | `m*(w^2+h^2)/12` | e3 (z-axis) | w (e1) and h (e2) |

### Unified picture across 2D and 3D

In both cases, mass enters where **Hodge-dual** spaces are connected:

| Dimension | Mass coupling | Direction | Reason |
| --- | --- | --- | --- |
| PGA2DP | `I[0,1]=+m`, `I[1,0]=-m` | translational vel <-> linear momentum | e1<->e31, e2<->e23 are complements |
| PGA3DP | `I[0,3]=I[1,4]=I[2,5]=m` | angular vel -> linear momentum | e23<->e41, e31<->e42, e12<->e43 are complements |

Moments of inertia always appear where a component maps **to its own Hodge-dual in the
output space**, i.e., at the rotational sub-block diagonal (2D: `I[2,2]`;
3D: `I[3,0]`, `I[4,1]`, `I[5,2]`).
