# Inertia Map Analysis

## Matrix from Document (Section 5)

The inertia tensor for rotation about point **X** = (Xₓ, Xᵧ, Xᵤ):

```text
I[Ω] = m · [ Xz²       0         -Xx·Xz     ]
           [ 0         Xz²       -Xy·Xz     ]
           [ -Xx·Xz    -Xy·Xz    Xx² + Xy²  ]
```

## Symbolic Invertibility

**Status**: ✅ Can be inverted symbolically

### Determinant

```text
det(I[Ω]) = m³ · Xz² · (Xx² + Xy² + Xz²) = m³ · Xz² · |X|²
```

where `|X|² = Xx² + Xy² + Xz²` is the squared magnitude of position vector **X**.

### Existence Conditions

The inverse exists when `det ≠ 0`, requiring:

1. **m ≠ 0** (non-zero mass)
2. **Xz ≠ 0** (position not in xy-plane)
3. **|X| ≠ 0** (position not at origin)

### Symbolic Inverse Formula

```text
I[Ω]⁻¹ = 1/(m · Xz² · |X|²) ·

    [ Xz²·(Xy² + Xz²)    Xx·Xy·Xz         Xx·Xz·|X|²  ]
    [ Xx·Xy·Xz           Xz²·(Xx² + Xz²)  Xy·Xz·|X|²  ]
    [ Xx·Xz·|X|²         Xy·Xz·|X|²       Xz²·|X|²    ]
```

## Singularity Analysis

### Critical Singularity: Xz = 0

When Xz → 0, the matrix becomes **singular**. This has physical meaning:

- The inertia tensor represents rotation about point **X** = (Xₓ, Xᵧ, Xᵤ)
- When Xz = 0, the point lies in the xy-plane
- The tensor degenerates because the z-component constraint vanishes

### Numerical Considerations

**Condition Number**: Deteriorates as Xz → 0

```text
κ(I[Ω]) ≈ |X|² / Xz²    (for small Xz)
```

The matrix becomes increasingly ill-conditioned when the rotation point approaches the
xy-plane.

## Implementation Recommendations

### When to Use Symbolic Inversion

✅ Use symbolic formula when:

- Xz is guaranteed bounded away from zero (`|Xz| > ε` for some threshold ε)
- Modeling rotations about points with significant z-coordinate
- Need exact algebraic expressions for further derivations

### When to Use Numerical Inversion

⚠️ Use numerical methods when:

- Xz might approach zero (rotation points near xy-plane)
- Need robust handling across full configuration space
- Working with measured/noisy data where Xz uncertainty is significant

### Suggested Strategy

For robust implementation:

```cpp
if (std::abs(X.z) > threshold) {
    // Use symbolic inversion (faster, exact)
    return symbolic_inverse(m, X);
} else {
    // Use numerical inversion with regularization
    return numerical_inverse_with_svd(I_matrix);
}
```

Threshold recommendation: `threshold ≈ 0.1 · |X|` to maintain condition number κ < 100.

## Geometric Interpretation

This inertia tensor structure suggests:

- **Parallel axis theorem** application with special z-axis treatment
- Rotation about arbitrary point **X** decomposed into z-direction and xy-plane components
- The block structure reflects separation of in-plane (xy) and out-of-plane (z) dynamics

## Further Analysis Needed

- [ ] Verify physical origin of this specific tensor form
- [ ] Check if this relates to 2D rigid body dynamics embedded in 3D
- [ ] Investigate whether GA formulation can eliminate the singularity
- [ ] Compare with standard parallel axis theorem formulation

---

## 6×6 Inertia Matrix Analysis (PGA3DP Form)

### Matrix Definition

From section 5 of the documentation (equations 612-620):

```text
I[Ω] = m ·
    [ 0              X_z·X_w       -X_y·X_w       X_w²           0              0            ]
    [-X_z·X_w        0              X_x·X_w       0              X_w²           0            ]
    [ X_y·X_w       -X_x·X_w        0             0              0              X_w²         ]
    [ X_y² + X_z²   -X_x·X_y       -X_x·X_z       0             -X_z·X_w        X_y·X_w      ]
    [-X_x·X_y        X_x² + X_z²   -X_y·X_z       X_z·X_w        0             -X_x·X_w      ]
    [-X_x·X_z       -X_y·X_z        X_x² + X_y²  -X_y·X_w        X_x·X_w        0            ]
```

where:

- `(X_x, X_y, X_z, X_w)` are homogeneous coordinates in PGA3DP
- `X_w` is the homogeneous coordinate (0 for ideal points/directions, ≠0 for finite
  points)
- `m` is the mass

### Invertibility Analysis

**Status**: ✅ Invertible under specific conditions

#### Determinant Analysis

For the matrix M (without scalar m), invertibility depends on:

1. **m ≠ 0** (non-zero mass - physically required)
2. **X_w ≠ 0** (finite point, not at infinity in projective space)

#### Critical Singularity: X_w = 0

When `X_w = 0` (point at infinity), the matrix becomes:

```text
M = [ 0            0           0          0    0    0  ]
    [ 0            0           0          0    0    0  ]
    [ 0            0           0          0    0    0  ]
    [ X_y² + X_z² -X_x·X_y    -X_x·X_z    0    0    0  ]
    [-X_x·X_y      X_x² + X_z² -X_y·X_z   0    0    0  ]
    [-X_x·X_z     -X_y·X_z     X_x² + X_y² 0    0    0  ]
```

This is **singular** (determinant = 0) due to the 3×3 zero block in upper-left and
upper-right.

**Physical Interpretation**: The inertia map is undefined for ideal points (directions at
infinity) in PGA, which makes physical sense - you cannot define inertia about a
direction, only about a position.

### Matrix Structure

The 6×6 matrix has special block structure:

- **Upper-left 3×3**: Skew-symmetric block (antisymmetric)
- **Upper-right 3×3**: Diagonal-like structure with X_w
- **Lower-left 3×3**: Symmetric block (moment of inertia components)
- **Lower-right 3×3**: Skew-symmetric block with X_w

This structure reflects the coupling between:

- Rotational components (bivector part of angular velocity)
- Translational components (vector part related to screw motion)

### Implementation Strategy

#### Recommendation: **Numerical Inversion**

For this 6×6 inertia matrix in physics simulation context, **numerical inversion is
strongly recommended**.

**Rationale:**

1. **Matrix Size**: 6×6 symbolic inverse would produce 36 expressions, each with
   potentially hundreds of terms - unwieldy and difficult to maintain

2. **Physics Application**: Runtime calculations in rigid body dynamics naturally use
   numerical methods

3. **Computational Efficiency**: Numerical inversion of 6×6 is extremely fast with
   optimized libraries (Eigen, LAPACK)

4. **Conditioning**: For physical configurations with `X_w ≠ 0`, the matrix should be
   reasonably well-conditioned

5. **Maintainability**: Easier to extend/modify for different physics models

6. **Integration**: Natural fit with physics simulation pipelines

#### Suggested Implementation

```cpp
// Using Eigen library
Eigen::Matrix<value_t, 6, 6> I_omega = /* construct from X_x, X_y, X_z, X_w, m */;

// Option 1: Direct inversion (fast, for well-conditioned cases)
Eigen::Matrix<value_t, 6, 6> I_omega_inv = I_omega.inverse();

// Option 2: With condition number check (robust)
Eigen::JacobiSVD<Eigen::Matrix<value_t, 6, 6>> svd(I_omega);
value_t cond = svd.singularValues()(0) / svd.singularValues()(5);

if (cond < condition_threshold) {
    I_omega_inv = I_omega.inverse();
} else {
    // Handle ill-conditioned case (e.g., regularization or warning)
    throw std::runtime_error("Inertia matrix is ill-conditioned");
}

// Option 3: SVD-based pseudo-inverse (most robust)
I_omega_inv = I_omega.completeOrthogonalDecomposition().pseudoInverse();
```

#### When to Consider Symbolic Inversion

Only consider symbolic inversion if:

- Theoretical analysis/derivation for papers
- Need exact closed-form solutions for equations of motion
- Formal verification of mathematical properties
- Educational purposes (demonstrating structure)

**Not recommended for**: Production physics code, runtime simulation, interactive
applications

### Numerical Stability Considerations

**Condition Number Behavior**:

- Well-conditioned when `X_w` is bounded away from zero
- Deteriorates as `X_w → 0` (approaching ideal points)
- Mass scaling (`m`) factors out and doesn't affect conditioning

**Recommended Safeguards**:

```cpp
// Check for degenerate cases before inversion
const value_t X_w_threshold = 1e-6; // Adjust based on physical units

if (std::abs(X_w) < X_w_threshold) {
    throw std::runtime_error("Cannot compute inertia map for ideal points (X_w ≈ 0)");
}

if (std::abs(m) < mass_threshold) {
    throw std::runtime_error("Cannot compute inertia map for massless objects");
}
```

### Geometric Interpretation (PGA3DP Context)

This 6×6 structure represents:

- **Motor algebra**: Coupling rotations and translations in rigid body motion
- **Screw theory**: Inertia properties for general screw motions (not just pure rotations)
- **Projective geometry**: Natural handling of points at finite positions (X_w ≠ 0) vs.
  directions at infinity (X_w = 0)

The block structure separates:

- Pure rotational inertia effects
- Translation-rotation coupling (Coriolis-like terms)
- Translational inertia effects

### Further Investigation

- [ ] Verify connection to motor algebra and screw theory formulations
- [ ] Determine physical meaning of each 3×3 block
- [ ] Check if there's a relationship to spatial inertia matrices in robotics
- [ ] Investigate whether PGA formulation provides advantages over traditional 6×6 spatial
  inertia
- [ ] Compare numerical stability with other inertia representations
