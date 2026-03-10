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

### Mapping classical quantities to PGA2DP

Integrating `get_point_inertia(dm, X)` over a rigid body gives for an arbitrary body frame
with center of mass at `(cx = (1/m)·∫x dm, cy = (1/m)·∫y dm)`:

```text
I = [  0      m       -m*cy        ]
    [ -m      0        m*cx        ]
    [ -m*cx  -m*cy    J_cm+m*r^2  ]
```

where `r^2 = cx^2+cy^2` and `J_cm` is the classical scalar moment of inertia about the CoM
J_cm = ∫(x²+y²) dm.

**In the body frame (CoM at origin, cx = cy = 0):**

```text
I_body = [  0    m    0 ]
         [ -m    0    0 ]
         [  0    0    J ]
```

| Classical quantity | PGA2DP entry | Condition |
| --- | --- | --- |
| Total mass `m` | `I[0,1]=+m`, `I[1,0]=-m` | always |
| Scalar moment of inertia `J` about CoM | `I[2,2]=J` | CoM at origin |
| CoM offset `(cx,cy)` | `I[0,2]=-m*cy`, `I[1,2]=+m*cx` (and transposed) | general frame |
| Parallel axis correction | `I[2,2]=J_cm+m*(cx^2+cy^2)` | general frame |

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

### Mapping classical quantities to PGA3DP

Integrating `get_point_inertia(dm, X)` over a rigid body gives for the **body frame
(CoM at origin)**:

```text
         [ vx    vy    vz  |  mx    my    mz  ]
  (e41)  [  0     0     0  |   m     0     0  ]
  (e42)  [  0     0     0  |   0     m     0  ]   upper-right = m*I3
  (e43)  [  0     0     0  |   0     0     m  ]
         [-----------------|------------------]
  (e23)  [ Ixx  -Ixy  -Ixz |   0     0     0  ]
  (e31)  [-Ixy   Iyy  -Iyz |   0     0     0  ]   lower-left = J_cm
  (e12)  [-Ixz  -Iyz   Izz |   0     0     0  ]
```

Classical 3D inertia tensor components:

```text
I_xx = ∫(y²+z²) dm,   I_yy = ∫(x²+z²) dm,   I_zz = ∫(x²+y²) dm
I_xy = ∫xy dm,        I_xz = ∫xz dm,        I_yz = ∫yz dm
```

| Classical quantity | PGA3DP block | Entries |
| --- | --- | --- |
| Total mass `m` | upper-right diagonal (m*I3) | `I[0,3]=I[1,4]=I[2,5]=m` |
| Diagonal moments `I_xx, I_yy, I_zz` | lower-left diagonal | `I[3,0], I[4,1], I[5,2]` |
| Products of inertia `I_xy, I_xz, I_yz` | lower-left off-diagonal | `I[3,1]=I[4,0]=-I_xy`, `I[3,2]=I[5,0]=-I_xz`, `I[4,2]=I[5,1]=-I_yz` |
| CoM offset (general frame) | coupling in all four 3x3 blocks | via parallel axis theorem |

Note: in the **principal-axis body frame** all products of inertia vanish (`I_xy=I_xz=I_yz=0`),
leaving only the three diagonal moments — exactly the cuboid formula below.

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

The basis ordering differs between 2D and 3D, which shifts where mass and moments appear:

**PGA2DP** — Vec2dp input (indices 0,1 = translational e1,e2; index 2 = angular e3),
BiVec2dp output (indices 0,1 = linear momentum e31,e32; index 2 = angular momentum e12):

- Mass occupies the **upper-left 2×2** (translational vel → linear momentum):
  `I[0,1]=+m`, `I[1,0]=-m` — antisymmetric due to Hodge crossing `!e1=e32`, `!e2=-e31`
- Moments occupy the **lower-right 1×1** (angular vel → angular momentum):
  `I[2,2]=J` — Hodge crossing `!e3=-e12` connects angular velocity e3 to angular momentum e12

**PGA3DP** — BiVec3dp input (indices 0-2 = translational e41,e42,e43; indices 3-5 = angular e23,e31,e12),
BiVec3dp output (indices 0-2 = linear momentum e41,e42,e43; indices 3-5 = angular momentum e23,e31,e12):

- Mass occupies the **upper-right 3×3** (angular vel → linear momentum):
  `I[0,3]=I[1,4]=I[2,5]=m` — Hodge crossing `!e23=e41`, `!e31=e42`, `!e12=e43`
- Moments occupy the **lower-left 3×3** (translational vel → angular momentum):
  `I[3,0], I[4,1], I[5,2]` — classical J_cm diagonal

In both cases mass is placed by `I_mass = -m * !` (negative Hodge complement) in the
sub-block connecting the Hodge-dual pair of velocity and momentum spaces. The different
block positions (upper-left in 2D vs upper-right in 3D) arise from the different basis
orderings: in 2D translational generators are at lower indices (0,1) and angular at
index 2; in 3D translational generators (e41,e42,e43) are at lower indices (0-2) and
angular generators (e23,e31,e12) at higher indices (3-5).

### Unified mapping: classical quantities to PGA (body frame, CoM at origin)

| Classical quantity | PGA2DP (3x3) | PGA3DP (6x6) |
| --- | --- | --- |
| Total mass `m` | upper-left 2×2: `I[0,1]=+m`, `I[1,0]=-m` | upper-right 3×3: `I[0,3]=I[1,4]=I[2,5]=m` |
| Scalar/diagonal moments `J` / `I_xx,I_yy,I_zz` | lower-right 1×1: `I[2,2]=J` | lower-left 3×3 diagonal: `I[3,0], I[4,1], I[5,2]` |
| Products of inertia | n/a (scalar J in 2D) | lower-left 3×3 off-diagonal: `I[3,1]=I[4,0]=-I_xy`, `I[3,2]=I[5,0]=-I_xz`, `I[4,2]=I[5,1]=-I_yz` |
| CoM offset (general frame) | `I[0,2], I[1,2], I[2,0], I[2,1]` non-zero | all four 3×3 blocks non-zero |

---

## General frame (CoM not at origin)

### PGA2DP: general frame with CoM at (cx, cy)

The 3×3 matrix splits into a (2+1)×(2+1) block structure:

```text
          [ trans (e1,e2) | angular (e3) ]
trans     [   0    +m     |   -m·cy      ]   upper-left: mass (invariant)
(e31,e32) [  -m     0     |   +m·cx      ]   upper-right: angular→linear coupling
          [---------------|-------------]
angular   [  -m·cx -m·cy  | J_cm+m·r²   ]   lower-left: trans→angular coupling
(e12)                                        lower-right: moment (with parallel axis)
```

- **Upper-left 2×2** `[0, m; -m, 0]`: pure antisymmetric mass block — **frame-invariant**, always the same regardless of CoM position.
- **Upper-right 2×1** `[-m·cy; m·cx]`: angular velocity → linear momentum coupling. Zero only when CoM is at origin. Encodes `L_linear = m·(ω × r_cm)`.
- **Lower-left 1×2** `[-m·cx, -m·cy]`: translational velocity → angular momentum coupling. Transpose of upper-right (up to sign).
- **Lower-right 1×1** `J_cm + m·(cx²+cy²)`: classical parallel axis theorem.

The CoM offset `(cx, cy)` only populates the cross-coupling blocks; the mass block is rigid.

### PGA3DP: general frame with CoM at r = (cx, cy, cz)

The 6×6 matrix splits into four 3×3 blocks. In the body frame the upper-left and
lower-right are zero; in a general frame they become antisymmetric cross-product matrices:

```text
         [ vx      vy      vz    |  mx      my      mz    ]
(e41)    [   0    +m·cz  -m·cy   |   m       0       0    ]
(e42)    [ -m·cz    0    +m·cx   |   0       m       0    ]   upper-right = m·I₃ (invariant)
(e43)    [ +m·cy  -m·cx    0     |   0       0       m    ]
         [---------------------- | ----------------------- ]
(e23)    [ I_xx  -I_xy  -I_xz   |   0    -m·cz  +m·cy   ]
(e31)    [-I_xy   I_yy  -I_yz   | +m·cz    0    -m·cx   ]   lower-right = -upper-left
(e12)    [-I_xz  -I_yz   I_zz   | -m·cy  +m·cx    0     ]
```

where the lower-left inertia tensor entries include the parallel axis correction:

```text
I_xx = I_xx_cm + m·(cy²+cz²)       I_xy = I_xy_cm + m·cx·cy
I_yy = I_yy_cm + m·(cx²+cz²)       I_xz = I_xz_cm + m·cx·cz
I_zz = I_zz_cm + m·(cx²+cy²)       I_yz = I_yz_cm + m·cy·cz
```

Block summary for the general frame:

| Block | Body frame | General frame | Physical meaning |
| --- | --- | --- | --- |
| Upper-left [0:3, 0:3] | `0` | `m·[r×]ᵀ` | translational vel → linear mom (CoM coupling) |
| Upper-right [0:3, 3:6] | `m·I₃` | `m·I₃` | angular vel → linear mom (always m, invariant) |
| Lower-left [3:6, 0:3] | `J_cm` | `J_cm + parallel axis` | translational vel → angular mom |
| Lower-right [3:6, 3:6] | `0` | `m·[r×]` | angular vel → angular mom (CoM coupling) |

where `[r×]` is the 3×3 antisymmetric cross-product matrix of the CoM position vector:

```text
[r×] = [  0   -cz   cy ]
       [  cz    0  -cx ]
       [ -cy   cx    0 ]
```

Note that **upper-left = −lower-right = lower-right^T**: the two coupling blocks are the
same antisymmetric matrix with opposite sign. This reflects the fact that the coupling
from translational velocity to linear momentum and from angular velocity to angular
momentum are both driven by the same CoM offset — they are conjugate couplings.

The upper-right `m·I₃` block is **frame-invariant**: it always equals `m·I₃` regardless
of where the CoM is, because mass is a Hodge-dual coupling between ideal and real lines
that does not depend on position.

This 6×6 structure is the PGA form of the classical **spatial inertia matrix** from
rigid body dynamics (Featherstone notation), expressed in the (translational, angular)
velocity ordering of PGA3DP.

### Considerations for comparison of classical inertia and the pga version

## PGA2DP — 3×3 matrix

Integrating get_point_inertia(dm, X) over a body gives:

I = [ 0       m        -m·cy  ]
    [ -m      0         m·cx  ]
    [ -m·cx  -m·cy   J_origin ]
where cx = (1/m)·∫x dm, cy = (1/m)·∫y dm (center of mass), and J_origin = ∫(x²+y²) dm.

By the parallel axis theorem: J_origin = J_cm + m·(cx²+cy²).

In the body frame (CoM at origin): cx = cy = 0, so:

I_body = [ 0    m    0  ]
         [ -m   0    0  ]
         [ 0    0    J  ]
The mapping from classical to PGA2DP is then simply:

Classical quantity PGA2DP location
mass m I[0,1] = +m, I[1,0] = -m
scalar moment of inertia J about CoM I[2,2] = J
CoM offset if not at origin I[0,2] = -m·cy, I[1,2] = m·cx, and transposed; I[2,2] = J_cm +
m·r²

## PGA3DP — 6×6 matrix

Integrating get_point_inertia(dm, X) over a body, for a body with CoM at origin:

         [ vx  vy  vz | mx  my  mz ]
   (e41) [  0   0   0 |  m   0   0  ]   upper-right = m·I₃
   (e42) [  0   0   0 |  0   m   0  ]
   (e43) [  0   0   0 |  0   0   m  ]
         [-----------+------------  ]
   (e23) [ Ixx -Ixy -Ixz|  0   0   0  ]  lower-left = classical J_cm
   (e31) [-Ixy  Iyy -Iyz|  0   0   0  ]
   (e12) [-Ixz -Iyz  Izz|  0   0   0  ]
where J_cm is the classical 3D inertia tensor about the CoM:

J_cm = [ I_xx  -I_xy  -I_xz ]      with  I_xx = ∫(y²+z²) dm
       [-I_xy   I_yy  -I_yz ]            I_yy = ∫(x²+z²) dm
       [-I_xz  -I_yz   I_zz ]            I_zz = ∫(x²+y²) dm
                                          I_xy = ∫xy dm,  etc.
The mapping from classical to PGA3DP (body frame, CoM at origin):

Classical quantity PGA3DP location
mass m upper-right diagonal: I[0,3]=I[1,4]=I[2,5]=m
diagonal moments I_xx, I_yy, I_zz I[3,0], I[4,1], I[5,2]
products of inertia I_xy, I_xz, I_yz off-diagonal of lower-left: I[3,1]=I[4,0]=-I_xy, etc.
CoM offset if not at origin coupling terms appear in all four 3×3 blocks
Unified pattern across 2D and 3D:

The PGA inertia matrix in the body frame (CoM at origin) always decomposes cleanly:

Mass occupies the Hodge-dual block (off-diagonal 2D, upper-right 3D) — because Newton's p
= m·v crosses dual spaces Moments of inertia occupy the same-grade block (diagonal 2D at
[2,2], lower-left 3×3 in 3D) — these connect translational velocity to angular momentum
All products of inertia (off-diagonal in J_cm) vanish in the body principal-axis frame —
just as classically
