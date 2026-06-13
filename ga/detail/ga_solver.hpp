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
                    for (int k = 0; k <= i - 1; ++k)
                        sum -= a[i, k] * a[k, j];
                    a[i, j] = sum;
                }
            }
        }
        double aamax = 0.0;
        for (int i = j; i <= ubound; ++i) {
            double sum = a[i, j];
            if (j > 0) {
                for (int k = 0; k <= j - 1; ++k)
                    sum -= a[i, k] * a[k, j];
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
            for (int i = j + 1; i <= ubound; ++i)
                a[i, j] *= dum;
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
            for (int j = ii; j <= i - 1; ++j)
                sum -= a[i, j] * b[j];
        }
        else if (sum != 0.0) {
            ii = i;
        }
        b[i] = sum;
    }

    for (int i = ubound; i >= 0; --i) {
        double sum = b[i];
        if (i < ubound) {
            for (int j = i + 1; j <= ubound; ++j)
                sum -= a[i, j] * b[j];
        }
        b[i] = sum / a[i, i];
    }
}


/////////////////////////////////////////////////////////////////////////////////////////
// Convenience: solve a small dense system A x = b given as flat ROW-MAJOR std::vector
// (A has n*n entries, b has n). Wraps lu_decomp + lu_backsubs -- the factorization is
// done in double, then cast back to T -- so the solver lives in ONE place for both 2D and
// 3D rigid-body physics (the joint-space mass matrix in forward dynamics). Returns x.
// Throws Solver_error on a singular matrix.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> lu_solve(std::vector<T> const& A_in, std::vector<T> const& b_in, size_t n)
{
    std::vector<double> a(n * n), b(n);
    for (size_t i = 0; i < n * n; ++i)
        a[i] = static_cast<double>(A_in[i]);
    for (size_t i = 0; i < n; ++i)
        b[i] = static_cast<double>(b_in[i]);

    std::vector<int> perm(n);
    std::mdspan<double, std::dextents<size_t, 2>> am(a.data(), n, n);
    std::mdspan<int, std::dextents<size_t, 1>> pm(perm.data(), n);
    lu_decomp(am, pm);

    std::mdspan<double const, std::dextents<size_t, 2>> ac(a.data(), n, n);
    std::mdspan<int const, std::dextents<size_t, 1>> pc(perm.data(), n);
    std::mdspan<double, std::dextents<size_t, 1>> bm(b.data(), n);
    lu_backsubs(ac, pc, bm);

    std::vector<T> x(n);
    for (size_t i = 0; i < n; ++i)
        x[i] = static_cast<T>(b[i]);
    return x;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Least-squares / minimum-norm dense solve of a (possibly non-square) system A x = b,
// where A is m x ncols (flat ROW-MAJOR, m = b.size()). Returns x (length ncols). Routed
// through the shared square lu_solve via the normal equations, so the whole library
// shares ONE dense LU. Three regimes (the Moore-Penrose pseudo-inverse solution x = A^+
// b):
//
//   ncols == m : square          -> solve A x = b directly
//   ncols >  m : underdetermined  -> minimum-norm   x = A^T (A A^T)^-1 b
//   ncols <  m : overdetermined   -> least-squares  x = (A^T A)^-1 A^T b
//
// Domain- and dimension-agnostic (pure linear algebra over T): used by the closed-loop
// constraint solver for the position Newton step (A = constraint Jacobian G, b = -g), the
// velocity / acceleration distribution, and the velocity projection -- in both 2D and 3D,
// and available to any other caller (it carries no GA or physics knowledge).
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> lstsq_solve(std::vector<T> const& A, std::vector<T> const& b, size_t ncols)
{
    size_t const m = b.size();

    if (ncols == m) {
        return lu_solve(A, b, m); // square A
    }
    if (ncols > m) {
        // minimum-norm: solve (A A^T) y = b (m x m), then x = A^T y
        std::vector<T> AAt(m * m, T(0));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j) {
                T s = T(0);
                for (size_t k = 0; k < ncols; ++k)
                    s += A[i * ncols + k] * A[j * ncols + k];
                AAt[i * m + j] = s;
            }
        std::vector<T> const y = lu_solve(AAt, b, m);
        std::vector<T> x(ncols, T(0));
        for (size_t k = 0; k < ncols; ++k) {
            T s = T(0);
            for (size_t i = 0; i < m; ++i)
                s += A[i * ncols + k] * y[i];
            x[k] = s;
        }
        return x;
    }
    // overdetermined: normal equations (A^T A) x = A^T b (ncols x ncols)
    std::vector<T> AtA(ncols * ncols, T(0)), Atb(ncols, T(0));
    for (size_t a = 0; a < ncols; ++a) {
        for (size_t i = 0; i < m; ++i)
            Atb[a] += A[i * ncols + a] * b[i];
        for (size_t bb = 0; bb < ncols; ++bb) {
            T s = T(0);
            for (size_t i = 0; i < m; ++i)
                s += A[i * ncols + a] * A[i * ncols + bb];
            AtA[a * ncols + bb] = s;
        }
    }
    return lu_solve(AtA, Atb, ncols);
}


/////////////////////////////////////////////////////////////////////////////////////////
// Bordered (saddle-point / KKT) dense solve. For the equality-constrained system
//
//   | M   G^T | | x |   | f |
//   |         | |   | = |   |
//   | G   0   | | l |   | g |
//
// with M (n x n) symmetric, G (m x n) the constraint Jacobian (both flat ROW-MAJOR), f
// (length n) and g (length m) the right-hand sides. Builds the (n+m) x (n+m) bordered
// matrix and solves it through the shared square lu_solve. Returns x (length n); writes
// the Lagrange multipliers l (length m) into `lambda_out` if non-null.
//
// Domain- and dimension-agnostic (pure linear algebra over T): the closed-loop dynamics
// (2D and 3D) use it with M = joint-space mass matrix, G = constraint Jacobian, f =
// generalised force tau, g = -G-dot q-dot, to get the joint accelerations x = q-ddot and
// the constraint forces l. Carries no GA or physics knowledge.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> kkt_solve(std::vector<T> const& M, std::vector<T> const& G,
                         std::vector<T> const& f, std::vector<T> const& g, size_t n,
                         size_t m, std::vector<T>* lambda_out = nullptr)
{
    size_t const N = n + m;
    std::vector<T> K(N * N, T(0)), r(N, T(0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            K[i * N + j] = M[i * n + j]; // M block (top-left)
    for (size_t c = 0; c < m; ++c)
        for (size_t j = 0; j < n; ++j) {
            K[j * N + (n + c)] = G[c * n + j]; // G^T block (top-right)
            K[(n + c) * N + j] = G[c * n + j]; // G  block (bottom-left)
        }
    for (size_t i = 0; i < n; ++i)
        r[i] = f[i];
    for (size_t c = 0; c < m; ++c)
        r[n + c] = g[c];

    std::vector<T> const sol = lu_solve(K, r, N);
    if (lambda_out) lambda_out->assign(sol.begin() + n, sol.end());
    return std::vector<T>(sol.begin(), sol.begin() + n);
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
    for (size_t i = 0; i < n; ++i)
        result *= a[i, i];

    int swaps = 0;
    for (size_t i = 0; i < n; ++i) {
        if (perm[i] != static_cast<int>(i)) ++swaps;
    }

    return static_cast<T>((swaps % 2 == 0) ? result : -result);
}

} // namespace hd::ga
