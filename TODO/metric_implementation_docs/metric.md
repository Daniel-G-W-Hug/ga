# Extended Metric Computation for Geometric Algebras

## Table of Contents

1. [Overview](#overview)
2. [Theoretical Background](#theoretical-background)
3. [Algebra Specifications](#algebra-specifications)
4. [Algorithm Description](#algorithm-description)
5. [C++ Implementation](#c-implementation)
6. [Expected Results](#expected-results)
7. [Usage Examples](#usage-examples)

-----

## Overview

This document provides a complete specification for computing the extended metric (outermorphism) for various geometric algebras including:

- **EGA** (Euclidean Geometric Algebra)
- **STA** (Spacetime Algebra)
- **PGA** (Projective Geometric Algebra)
- **CGA** (Conformal Geometric Algebra)

The extended metric extends the base metric from vectors to all multivectors (blades of any grade).

-----

## Theoretical Background

### Metric Definition

In geometric algebra, a **metric** is a symmetric bilinear form that defines the geometric product between basis vectors. For basis vectors **e**ᵢ and **e**ⱼ, the metric is specified by their inner products:

```text
eᵢ · eⱼ = gᵢⱼ
```

where gᵢⱼ are the components of the metric tensor. This determines how basis vectors square:

- **e**ᵢ² = gᵢᵢ (can be +1, -1, or 0)

The geometric product follows:

```text
eᵢeⱼ = eᵢ · eⱼ + eᵢ ∧ eⱼ
```

### Extended Metric (Outermorphism)

The **extended metric** or **outermorphism** extends the metric from vectors to all multivectors (blades of any grade). It is denoted by angle brackets ⟨ ⟩.

For basis blades **B**ᵢ and **B**ⱼ, the extended metric matrix elements are:

```text
Gᵢⱼ = ⟨Bᵢ, Bⱼ⟩
```

For blades **A** = **a**₁ ∧ **a**₂ ∧ … ∧ **aₖ** and **B** = **b**₁ ∧ **b**₂ ∧ … ∧ **bₖ** of the same grade k:

```text
⟨A, B⟩ = det(G)
```

where G is the k×k Gram matrix with entries:

```text
G[p,q] = aₚ · bᵩ
```

### Key Properties

1. **Different grades**: If blades have different grades, ⟨**A**, **B**⟩ = 0
1. **Scalars**: ⟨1, 1⟩ = 1
1. **Vectors**: ⟨**v**, **w**⟩ = **v** · **w**
1. **Block-diagonal structure**: The extended metric is block-diagonal by grade
1. **Null vectors**: Within each grade block, off-diagonal elements can appear when null vectors are present (as in CGA with null basis)

-----

## Algebra Specifications

### EGA 2D - Euclidean Geometric Algebra 2D

- **Signature**: ℝ²’⁰’⁰
- **Basis**: {e₁, e₂}
- **Metric tensor**: diag(1, 1)
- **Properties**: All basis vectors square to +1
- **Extended metric**: Fully diagonal, all elements ±1
- **Basis ordering**: e₁, e₂

### EGA 3D - Euclidean Geometric Algebra 3D

- **Signature**: ℝ³’⁰’⁰
- **Basis**: {e₁, e₂, e₃}
- **Metric tensor**: diag(1, 1, 1)
- **Properties**: All basis vectors square to +1
- **Extended metric**: Fully diagonal, all elements ±1
- **Basis ordering**: e₁, e₂, e₃

### STA 3D - Spacetime Algebra 3D

- **Signature**: ℝ¹’³’⁰
- **Basis**: {e₀, e₁, e₂, e₃} (time, x, y, z)
- **Metric tensor**: diag(1, -1, -1, -1)
- **Properties**: Minkowski metric for spacetime
- **Extended metric**: Fully diagonal, elements ±1
- **Basis ordering**: e₀, e₁, e₂, e₃

### PGA 2D - Projective Geometric Algebra 2D

- **Signature**: ℝ²’⁰’¹
- **Basis**: {e₁, e₂, e₀} (e₀ is last, degenerate)
- **Metric tensor**: diag(1, 1, 0)
- **Properties**: e₀² = 0 (null/degenerate dimension)
- **Extended metric**: Diagonal with zeros for blades containing e₀
- **Basis ordering**: e₁, e₂, e₀

### PGA 3D - Projective Geometric Algebra 3D

- **Signature**: ℝ³’⁰’¹
- **Basis**: {e₁, e₂, e₃, e₀} (e₀ is last, degenerate)
- **Metric tensor**: diag(1, 1, 1, 0)
- **Properties**: e₀² = 0 (null/degenerate dimension)
- **Extended metric**: Diagonal with zeros for blades containing e₀
- **Basis ordering**: e₁, e₂, e₃, e₀

### CGA 2D - Conformal Geometric Algebra 2D (Diagonal Basis)

- **Signature**: ℝ³’¹’⁰
- **Basis**: {e₁, e₂, e₊, e₋}
- **Metric tensor**: diag(1, 1, 1, -1)
- **Properties**: Mixed signature with one negative square
- **Extended metric**: Fully diagonal, elements ±1
- **Basis ordering**: e₁, e₂, e₊, e₋

### CGA 3D - Conformal Geometric Algebra 3D (Diagonal Basis)

- **Signature**: ℝ⁴’¹’⁰
- **Basis**: {e₁, e₂, e₃, e₊, e₋}
- **Metric tensor**: diag(1, 1, 1, 1, -1)
- **Properties**: Mixed signature with one negative square
- **Extended metric**: Fully diagonal, elements ±1
- **Basis ordering**: e₁, e₂, e₃, e₊, e₋

### CGA 2D - Conformal Geometric Algebra 2D (Null Basis)

- **Signature**: ℝ³’¹’⁰
- **Basis**: {e₁, e₂, n₀, n∞}
- **Metric tensor**:
  
  ```text
  [1   0   0    0 ]
  [0   1   0    0 ]
  [0   0   0   -1 ]
  [0   0  -1    0 ]
  ```

- **Properties**:
  - n₀² = 0, n∞² = 0
  - n₀ · n∞ = -1 (off-diagonal coupling)
- **Extended metric**: Block-diagonal by grade with off-diagonal elements within each grade block
- **Basis ordering**: e₁, e₂, n₀, n∞

### CGA 3D - Conformal Geometric Algebra 3D (Null Basis)

- **Signature**: ℝ⁴’¹’⁰
- **Basis**: {e₁, e₂, e₃, n₀, n∞}
- **Metric tensor**:
  
  ```text
  [1   0   0   0    0 ]
  [0   1   0   0    0 ]
  [0   0   1   0    0 ]
  [0   0   0   0   -1 ]
  [0   0   0  -1    0 ]
  ```

- **Properties**:
  - n₀² = 0, n∞² = 0
  - n₀ · n∞ = -1 (off-diagonal coupling)
- **Extended metric**: Block-diagonal by grade with off-diagonal elements within each grade block
- **Basis ordering**: e₁, e₂, e₃, n₀, n∞

-----

## Algorithm Description

### Step 1: Generate All Basis Blades

For an n-dimensional algebra, generate all 2ⁿ basis blades ordered by:

1. Grade (0 to n)
1. Lexicographically within each grade

**Example for n=3**:

- Grade 0: {1}
- Grade 1: {e₁, e₂, e₃}
- Grade 2: {e₁₂, e₁₃, e₂₃}
- Grade 3: {e₁₂₃}

Total: 1 + 3 + 3 + 1 = 8 = 2³ blades

### Step 2: Compute Blade Inner Products

For each pair of basis blades Bᵢ and Bⱼ:

```text
function blade_inner_product(blade1, blade2, metric):
    // Different grades -> zero
    if len(blade1) ≠ len(blade2):
        return 0
    
    // Both scalars
    if len(blade1) = 0:
        return 1
    
    // Same grade k: compute Gram matrix determinant
    k = len(blade1)
    G = k×k matrix
    
    for p = 1 to k:
        for q = 1 to k:
            G[p,q] = metric[blade1[p], blade2[q]]
    
    return det(G)
```

### Step 3: Build Extended Metric Matrix

The extended metric is a 2ⁿ × 2ⁿ matrix where:

```text
ExtendedMetric[i,j] = blade_inner_product(blade_i, blade_j, metric)
```

-----

## C++ Implementation

### Complete Code

```cpp
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

// ============================================================================
// Matrix utilities
// ============================================================================

class Matrix {
public:
    std::vector<std::vector<double>> data;
    size_t rows, cols;
    
    Matrix(size_t r, size_t c) : rows(r), cols(c) {
        data.resize(r, std::vector<double>(c, 0.0));
    }
    
    double& operator()(size_t i, size_t j) {
        return data[i][j];
    }
    
    const double& operator()(size_t i, size_t j) const {
        return data[i][j];
    }
    
    static Matrix diagonal(const std::vector<double>& diag) {
        size_t n = diag.size();
        Matrix m(n, n);
        for (size_t i = 0; i < n; ++i) {
            m(i, i) = diag[i];
        }
        return m;
    }
    
    // Compute determinant using LU decomposition
    double determinant() const {
        if (rows != cols) return 0.0;
        if (rows == 0) return 1.0;
        if (rows == 1) return data[0][0];
        if (rows == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }
        
        // Create a copy for LU decomposition
        std::vector<std::vector<double>> a = data;
        int n = rows;
        double det = 1.0;
        
        for (int i = 0; i < n; ++i) {
            // Find pivot
            int pivot = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(a[j][i]) > std::abs(a[pivot][i])) {
                    pivot = j;
                }
            }
            
            if (std::abs(a[pivot][i]) < 1e-12) {
                return 0.0; // Singular matrix
            }
            
            // Swap rows if needed
            if (pivot != i) {
                std::swap(a[i], a[pivot]);
                det = -det;
            }
            
            det *= a[i][i];
            
            // Eliminate column
            for (int j = i + 1; j < n; ++j) {
                double factor = a[j][i] / a[i][i];
                for (int k = i + 1; k < n; ++k) {
                    a[j][k] -= factor * a[i][k];
                }
            }
        }
        
        return det;
    }
};

// ============================================================================
// Basis blade representation
// ============================================================================

struct Blade {
    std::vector<int> indices; // Sorted list of basis vector indices
    
    Blade() {}
    Blade(const std::vector<int>& idx) : indices(idx) {}
    
    int grade() const {
        return indices.size();
    }
    
    std::string toString(const std::vector<std::string>& basis_names) const {
        if (indices.empty()) return "1";
        std::string result = "e_";
        for (int idx : indices) {
            result += basis_names[idx];
        }
        return result;
    }
};

// ============================================================================
// Generate all basis blades
// ============================================================================

void generateBladesRecursive(int n, int start, int grade, 
                             std::vector<int>& current,
                             std::vector<Blade>& blades) {
    if (current.size() == grade) {
        blades.push_back(Blade(current));
        return;
    }
    
    for (int i = start; i < n; ++i) {
        current.push_back(i);
        generateBladesRecursive(n, i + 1, grade, current, blades);
        current.pop_back();
    }
}

std::vector<Blade> generateBasisBlades(int n) {
    std::vector<Blade> blades;
    
    for (int grade = 0; grade <= n; ++grade) {
        std::vector<int> current;
        generateBladesRecursive(n, 0, grade, current, blades);
    }
    
    return blades;
}

// ============================================================================
// Extended metric computation
// ============================================================================

double bladeInnerProduct(const Blade& blade1, const Blade& blade2, 
                         const Matrix& metric) {
    // Different grades -> zero
    if (blade1.grade() != blade2.grade()) {
        return 0.0;
    }
    
    // Both scalars
    if (blade1.grade() == 0) {
        return 1.0;
    }
    
    // Same grade k: compute determinant of Gram matrix
    int k = blade1.grade();
    Matrix gram(k, k);
    
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            gram(i, j) = metric(blade1.indices[i], blade2.indices[j]);
        }
    }
    
    return gram.determinant();
}

Matrix computeExtendedMetric(const Matrix& metric, 
                             std::vector<Blade>& blades) {
    int n = metric.rows;
    blades = generateBasisBlades(n);
    int size = blades.size();
    
    Matrix extended(size, size);
    
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            extended(i, j) = bladeInnerProduct(blades[i], blades[j], metric);
        }
    }
    
    return extended;
}

// ============================================================================
// Algebra metric definitions
// ============================================================================

Matrix ega2d_metric() {
    return Matrix::diagonal({1.0, 1.0});
}

Matrix ega3d_metric() {
    return Matrix::diagonal({1.0, 1.0, 1.0});
}

Matrix sta4d_metric() {
    return Matrix::diagonal({1.0, -1.0, -1.0, -1.0});
}

Matrix pga2d_metric() {
    // e1, e2, e0 where e0^2 = 0
    return Matrix::diagonal({1.0, 1.0, 0.0});
}

Matrix pga3d_metric() {
    // e1, e2, e3, e0 where e0^2 = 0
    return Matrix::diagonal({1.0, 1.0, 1.0, 0.0});
}

Matrix cga2d_metric() {
    return Matrix::diagonal({1.0, 1.0, 1.0, -1.0});
}

Matrix cga3d_metric() {
    return Matrix::diagonal({1.0, 1.0, 1.0, 1.0, -1.0});
}

Matrix cga2d_metric_null() {
    Matrix m(4, 4);
    m(0, 0) = 1.0;   // e1
    m(1, 1) = 1.0;   // e2
    m(2, 3) = -1.0;  // n0 · n∞
    m(3, 2) = -1.0;  // n∞ · n0
    return m;
}

Matrix cga3d_metric_null() {
    Matrix m(5, 5);
    m(0, 0) = 1.0;   // e1
    m(1, 1) = 1.0;   // e2
    m(2, 2) = 1.0;   // e3
    m(3, 4) = -1.0;  // n0 · n∞
    m(4, 3) = -1.0;  // n∞ · n0
    return m;
}

// ============================================================================
// Display utilities
// ============================================================================

bool isDiagonal(const Matrix& m, double tol = 1e-10) {
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            if (i != j && std::abs(m(i, j)) > tol) {
                return false;
            }
        }
    }
    return true;
}

void printExtendedMetric(const Matrix& extended, 
                        const std::vector<Blade>& blades,
                        const std::vector<std::string>& basis_names,
                        const std::string& name) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << name << "\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Matrix shape: " << extended.rows << "x" << extended.cols << "\n";
    
    bool is_diagonal = isDiagonal(extended);
    std::cout << "Diagonal: " << (is_diagonal ? "Yes" : "No") << "\n";
    
    if (is_diagonal) {
        std::cout << "\nDiagonal elements:\n";
        for (size_t i = 0; i < blades.size(); ++i) {
            std::string blade_str = blades[i].toString(basis_names);
            std::cout << "  " << std::setw(10) << std::left << blade_str 
                     << ": " << std::setw(7) << std::fixed << std::setprecision(3)
                     << extended(i, i) << "\n";
        }
    } else {
        std::cout << "\nFull matrix (showing non-zero elements):\n";
        
        // Group by grade
        std::vector<std::vector<int>> grades;
        int max_grade = 0;
        for (const auto& blade : blades) {
            max_grade = std::max(max_grade, blade.grade());
        }
        grades.resize(max_grade + 1);
        
        for (size_t i = 0; i < blades.size(); ++i) {
            grades[blades[i].grade()].push_back(i);
        }
        
        for (int g = 0; g <= max_grade; ++g) {
            if (grades[g].empty()) continue;
            
            std::cout << "\n  Grade " << g << ":\n";
            
            for (int idx_i : grades[g]) {
                std::string blade_str = blades[idx_i].toString(basis_names);
                std::cout << "    " << std::setw(8) << std::left << blade_str << ": [";
                
                bool first = true;
                for (int idx_j : grades[g]) {
                    if (!first) std::cout << ", ";
                    first = false;
                    
                    double val = extended(idx_i, idx_j);
                    if (std::abs(val) > 1e-10) {
                        std::cout << std::setw(7) << std::fixed 
                                 << std::setprecision(3) << val;
                    } else {
                        std::cout << "   .   ";
                    }
                }
                std::cout << "]\n";
            }
        }
    }
}

// ============================================================================
// Main
// ============================================================================

int main() {
    struct AlgebraConfig {
        std::string name;
        Matrix metric;
        std::vector<std::string> basis_names;
        
        AlgebraConfig(const std::string& n, const Matrix& m, 
                     const std::vector<std::string>& b)
            : name(n), metric(m), basis_names(b) {}
    };
    
    std::vector<AlgebraConfig> algebras;
    algebras.push_back(AlgebraConfig("EGA 2D", ega2d_metric(), {"1", "2"}));
    algebras.push_back(AlgebraConfig("EGA 3D", ega3d_metric(), {"1", "2", "3"}));
    algebras.push_back(AlgebraConfig("STA 3D", sta4d_metric(), {"0", "1", "2", "3"}));
    algebras.push_back(AlgebraConfig("PGA 2D", pga2d_metric(), {"1", "2", "0"}));
    algebras.push_back(AlgebraConfig("PGA 3D", pga3d_metric(), {"1", "2", "3", "0"}));
    algebras.push_back(AlgebraConfig("CGA 2D (diagonal basis)", cga2d_metric(), 
                                    {"1", "2", "+", "-"}));
    algebras.push_back(AlgebraConfig("CGA 3D (diagonal basis)", cga3d_metric(), 
                                    {"1", "2", "3", "+", "-"}));
    algebras.push_back(AlgebraConfig("CGA 2D (null basis)", cga2d_metric_null(), 
                                    {"1", "2", "0", "∞"}));
    algebras.push_back(AlgebraConfig("CGA 3D (null basis)", cga3d_metric_null(), 
                                    {"1", "2", "3", "0", "∞"}));
    
    for (const auto& algebra : algebras) {
        std::vector<Blade> blades;
        Matrix extended = computeExtendedMetric(algebra.metric, blades);
        printExtendedMetric(extended, blades, algebra.basis_names, algebra.name);
    }
    
    return 0;
}
```

### Compilation

```bash
g++ -std=c++11 -O2 extended_metric.cpp -o extended_metric
./extended_metric
```

-----

## Expected Results

### Summary by Algebra Type

1. **EGA 2D/3D**:

- Fully diagonal extended metric
- All diagonal elements are ±1
- 4 elements (2D) or 8 elements (3D)

1. **STA 3D**:

- Fully diagonal extended metric
- Mix of +1 and -1 elements
- 16 elements total

1. **PGA 2D/3D**:

- Diagonal extended metric
- Contains zeros for any blade containing e₀
- Example PGA 2D: blades like e₀, e₁₀, e₂₀, e₁₂₀ have zero metric

1. **CGA (diagonal basis)**:

- Fully diagonal extended metric
- All diagonal elements are ±1
- 16 elements (2D) or 32 elements (3D)

1. **CGA (null basis)**:

- Block-diagonal by grade
- **Off-diagonal elements appear** within grade blocks
- These arise from n₀ · n∞ = -1 coupling
- Example: In grade 1, the metric block is not diagonal

### Detailed Example: PGA 2D

Basis blades (8 total):

```text
Grade 0: 1
Grade 1: e₁, e₂, e₀
Grade 2: e₁₂, e₁₀, e₂₀
Grade 3: e₁₂₀
```

Extended metric diagonal elements:

```text
⟨1⟩     = 1.000
⟨e₁⟩    = 1.000
⟨e₂⟩    = 1.000
⟨e₀⟩    = 0.000  ← zero because e₀² = 0
⟨e₁₂⟩   = 1.000
⟨e₁₀⟩   = 0.000  ← contains e₀
⟨e₂₀⟩   = 0.000  ← contains e₀
⟨e₁₂₀⟩  = 0.000  ← contains e₀
```

### Detailed Example: CGA 2D (Null Basis)

Grade 1 metric block (4×4):

```text
         e₁      e₂      n₀      n∞
e₁   [  1.000   0.000   0.000   0.000 ]
e₂   [  0.000   1.000   0.000   0.000 ]
n₀   [  0.000   0.000   0.000  -1.000 ]  ← off-diagonal!
n∞   [  0.000   0.000  -1.000   0.000 ]  ← off-diagonal!
```

Notice the off-diagonal elements at positions (2,3) and (3,2) with value -1, arising from n₀ · n∞ = -1.

-----

## Usage Examples

### Example 1: Computing Extended Metric for Custom Algebra

```cpp
// Define a custom metric (e.g., 2+1 dimensional spacetime)
Matrix custom_metric(3, 3);
custom_metric(0, 0) = 1.0;   // time
custom_metric(1, 1) = -1.0;  // x
custom_metric(2, 2) = -1.0;  // y

// Compute extended metric
std::vector<Blade> blades;
Matrix extended = computeExtendedMetric(custom_metric, blades);

// Display
std::vector<std::string> names = {"t", "x", "y"};
printExtendedMetric(extended, blades, names, "Custom 2+1 Spacetime");
```

### Example 2: Accessing Individual Elements

```cpp
// Get extended metric for a specific algebra
std::vector<Blade> blades;
Matrix extended = computeExtendedMetric(pga2d_metric(), blades);

// Access specific element
int blade_index_e12 = 4;  // depends on blade ordering
double metric_e12 = extended(blade_index_e12, blade_index_e12);

std::cout << "Metric of e₁₂: " << metric_e12 << "\n";
```

### Example 3: Finding Blade Index

```cpp
// Find the index of blade e₁₂ in the blade list
std::vector<Blade> blades = generateBasisBlades(3);

for (size_t i = 0; i < blades.size(); ++i) {
    if (blades[i].indices == std::vector<int>{0, 1}) {  // e₁₂
        std::cout << "Blade e₁₂ is at index " << i << "\n";
        break;
    }
}
```

-----

## Mathematical Verification

### Example Calculation: Bivector in EGA 3D

For blade e₁₂ = e₁ ∧ e₂ in EGA 3D:

Gram matrix:

```text
G = [ e₁·e₁  e₁·e₂ ]   [ 1  0 ]
    [ e₂·e₁  e₂·e₂ ] = [ 0  1 ]
```

Extended metric:

```text
⟨e₁₂, e₁₂⟩ = det(G) = 1×1 - 0×0 = 1
```

### Example Calculation: Bivector in CGA (Null Basis)

For blade n₀ ∧ e₁ in CGA:

Gram matrix (using n₀ as first vector, e₁ as second):

```text
G = [ n₀·n₀  n₀·e₁ ]   [ 0  0 ]
    [ e₁·n₀  e₁·e₁ ] = [ 0  1 ]
```

Extended metric:

```text
⟨n₀∧e₁, n₀∧e₁⟩ = det(G) = 0×1 - 0×0 = 0
```

For blade n₀ ∧ n∞:

Gram matrix:

```text
G = [ n₀·n₀  n₀·n∞ ]   [  0  -1 ]
    [ n∞·n₀  n∞·n∞ ] = [ -1   0 ]
```

Extended metric:

```text
⟨n₀∧n∞, n₀∧n∞⟩ = det(G) = 0×0 - (-1)×(-1) = -1
```

This shows the effect of null vectors on the extended metric.

-----

## Notes and Important Considerations

1. **Blade Ordering**: The code uses lexicographic ordering within each grade. This is important for consistent indexing.
1. **Numerical Precision**: The determinant computation uses a tolerance of 1e-12 for singularity detection.
1. **Memory Efficiency**: For large algebras (n > 5), the 2ⁿ × 2ⁿ matrix can become very large. Consider sparse matrix representations.
1. **Alternative Bases**: CGA can be represented in either diagonal basis (e₊, e₋) or null basis (n₀, n∞). The null basis shows off-diagonal structure more clearly.
1. **Signature Notation**: The notation ℝᵖ’ᵍ’ʳ means:

- p basis vectors with e² = +1
- q basis vectors with e² = -1
- r basis vectors with e² = 0 (degenerate)

-----

## References and Further Reading

- **Geometric Algebra**: Doran & Lasenby, “Geometric Algebra for Physicists”
- **Conformal GA**: Dorst et al., “Geometric Algebra for Computer Science”
- **Projective GA**: Gunn, “Geometric Algebras for Euclidean Geometry”
- **Implementation**: <https://bivector.net> for practical GA implementations

-----

## License

This code and documentation are provided for educational and research purposes.

-----

**Document Version**: 1.0  
**Last Updated**: 2024  
**Author**: Generated via Claude (Anthropic)
