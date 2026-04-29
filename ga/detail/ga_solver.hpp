#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// Small dense linear-system solver used by GA physics modelling.
//
// LU decomposition with partial pivoting plus back-substitution, and a
// matrix determinant built on top. Sized for the small inertia matrices
// produced by PGA rigid-body dynamics (3x3 in 2D, 6x6 in 3D), where
// pulling in a full BLAS/LAPACK stack would be excessive.
//
// Usage:
//
//   1.) LU decomposition of matrix A:
//       hd::ga::lu_decomp(A, perm);
//
//   2.) Solve A*x = b by back-substitution; the solution overwrites b
//       (can be repeated cheaply for many right-hand sides on the same A):
//       hd::ga::lu_backsubs(A, perm, b);
//
//   3.) Matrix determinant via LU factorization:
//       T d = hd::ga::det(A);
//
// Adapted from the hd utility library and made internal to the ga library
// so the physics ops carry no external dependency.
/////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>     // std::abs
#include <mdspan>    // std::mdspan, std::dextents, std::extents
#include <stdexcept> // std::runtime_error, std::invalid_argument
#include <string>    // std::string
#include <vector>    // std::vector (scratch storage in det / lu_decomp)

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// Solver_error: thrown by lu_decomp on a singular or malformed input.
//
// Caught locally in det() (singular matrix → returns 0). Higher-level
// physics callers (e.g. get_inertia_inverse) detect singularity by
// checking the determinant first and never let this exception escape.
/////////////////////////////////////////////////////////////////////////////////////////
struct Solver_error : std::runtime_error {
    explicit Solver_error(char const* msg) : std::runtime_error(msg) {}
    explicit Solver_error(std::string const& msg) : std::runtime_error(msg) {}
};


/////////////////////////////////////////////////////////////////////////////////////////
// LU decomposition with partial pivoting (Crout's algorithm with row pivoting).
//
//   a    - n×n matrix; on exit, contains L (below diag, unit diagonal implicit)
//          and U (on/above diag) packed in place.
//   perm - length-n permutation vector recording row swaps.
/////////////////////////////////////////////////////////////////////////////////////////
inline void lu_decomp(std::mdspan<double, std::dextents<size_t, 2>> a,
                      std::mdspan<int, std::dextents<size_t, 1>> perm)
{
    if (a.extent(0) != a.extent(1) || a.extent(0) != perm.extent(0)) {
        throw Solver_error("hd::ga::lu_decomp: matrix is non-square or "
                           "permutation-vector size incompatible.");
    }

    constexpr double TINY = 1.e-20;
    int const ubound = static_cast<int>(a.extent(0)) - 1;

    // Per-row scaling factors (for implicit pivoting).
    std::vector<double> vv(a.extent(0));
    for (int i = 0; i <= ubound; ++i) {
        double aamax = 0.0;
        for (int j = 0; j <= ubound; ++j) {
            if (std::abs(a[i, j]) > aamax) aamax = std::abs(a[i, j]);
        }
        if (aamax == 0.0) {
            throw Solver_error("hd::ga::lu_decomp: singular matrix.");
        }
        vv[i] = 1.0 / aamax;
    }

    int imax = 0;
    for (int j = 0; j <= ubound; ++j) {
        if (j > 0) {
            for (int i = 0; i <= j - 1; ++i) {
                double sum = a[i, j];
                if (i > 0) {
                    for (int k = 0; k <= i - 1; ++k) sum -= a[i, k] * a[k, j];
                    a[i, j] = sum;
                }
            }
        }
        double aamax = 0.0;
        for (int i = j; i <= ubound; ++i) {
            double sum = a[i, j];
            if (j > 0) {
                for (int k = 0; k <= j - 1; ++k) sum -= a[i, k] * a[k, j];
                a[i, j] = sum;
            }
            double const dum = vv[i] * std::abs(sum);
            if (dum >= aamax) {
                imax = i;
                aamax = dum;
            }
        }
        if (j != imax) {
            for (int k = 0; k <= ubound; ++k) {
                double const dum = a[imax, k];
                a[imax, k] = a[j, k];
                a[j, k] = dum;
            }
            vv[imax] = vv[j];
        }
        perm[j] = imax;
        if (j != ubound) {
            if (a[j, j] == 0.0) a[j, j] = TINY;
            double const dum = 1.0 / a[j, j];
            for (int i = j + 1; i <= ubound; ++i) a[i, j] *= dum;
        }
    }
    if (a[ubound, ubound] == 0.0) a[ubound, ubound] = TINY;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Back-substitution for the LU factorization produced by lu_decomp().
//
//   a    - LU-decomposed matrix from lu_decomp().
//   perm - permutation vector from lu_decomp().
//   b    - right-hand side; on exit, contains the solution x of A*x = b.
//
// Cheap to repeat for many right-hand sides against the same factorization.
/////////////////////////////////////////////////////////////////////////////////////////
inline void lu_backsubs(std::mdspan<double const, std::dextents<size_t, 2>> a,
                        std::mdspan<int const, std::dextents<size_t, 1>> perm,
                        std::mdspan<double, std::dextents<size_t, 1>> b)
{
    if (a.extent(0) != a.extent(1) || a.extent(0) != perm.extent(0) ||
        a.extent(0) != b.extent(0)) {
        throw Solver_error("hd::ga::lu_backsubs: matrix non-square, or "
                           "permutation / right-hand side size incompatible.");
    }

    int const ubound = static_cast<int>(a.extent(0)) - 1;

    int ii = -1; // sentinel: index of first non-zero element seen so far
    for (int i = 0; i <= ubound; ++i) {
        int const ll = perm[i];
        double sum = b[ll];
        b[ll] = b[i];
        if (ii != -1) {
            for (int j = ii; j <= i - 1; ++j) sum -= a[i, j] * b[j];
        }
        else if (sum != 0.0) {
            ii = i;
        }
        b[i] = sum;
    }

    for (int i = ubound; i >= 0; --i) {
        double sum = b[i];
        if (i < ubound) {
            for (int j = i + 1; j <= ubound; ++j) sum -= a[i, j] * b[j];
        }
        b[i] = sum / a[i, i];
    }
}


/////////////////////////////////////////////////////////////////////////////////////////
// Determinant of a square matrix via the LU factorization.
//
// The input span is not modified; data is copied into a local scratch
// buffer first. Returns T(0) for a singular matrix (rather than
// propagating Solver_error), which is the convention physics callers
// rely on for "is this inertia tensor invertible?" checks.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
T det(std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy> A)
{
    size_t const n = A.extent(0);

    if (n != A.extent(1)) {
        throw std::invalid_argument("hd::ga::det: matrix must be square.");
    }

    std::vector<double> data(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            data[i * n + j] = static_cast<double>(A[i, j]);
        }
    }

    std::mdspan<double, std::dextents<size_t, 2>> a(data.data(), n, n);
    std::vector<int> perm_data(n);
    std::mdspan<int, std::dextents<size_t, 1>> perm(perm_data.data(), n);

    try {
        lu_decomp(a, perm);
    }
    catch (Solver_error const&) {
        return T(0); // Singular matrix has determinant 0.
    }

    double result = 1.0;
    for (size_t i = 0; i < n; ++i) result *= a[i, i];

    int swaps = 0;
    for (size_t i = 0; i < n; ++i) {
        if (perm[i] != static_cast<int>(i)) ++swaps;
    }

    return static_cast<T>((swaps % 2 == 0) ? result : -result);
}

} // namespace hd::ga
