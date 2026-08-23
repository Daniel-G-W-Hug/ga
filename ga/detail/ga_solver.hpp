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
//   4.) Tridiagonal systems (compact/Pade finite differences, 1d diffusion) do NOT
//       need any of the above -- they are solved in O(n) without forming a matrix:
//       auto x = hd::ga::tridiag_solve(a, b, c, d);
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
// Tridiagonal solve (Thomas algorithm): forward elimination + back substitution.
//
//     a[i] x_{i-1} + b[i] x_i + c[i] x_{i+1} = d[i]
//
// with a[0] and c[n-1] outside the matrix and ignored. The three diagonals are given
// as separate length-n vectors -- the n*n matrix is never formed.
//
// This is the solver the COMPACT (Pade) finite-difference schemes ask for. A compact
// scheme states the derivative implicitly, coupling its values at i-1, i and i+1, so
// evaluating one over a grid IS a tridiagonal system: one row per node, three entries
// each. Routing that through lu_decomp/lu_backsubs would build the full n*n matrix and
// factor it in O(n^3) to recover an O(n) answer -- for a 10^3-node grid, ~10^6 stored
// doubles and ~10^9 operations in place of ~10^4. Hence this one, which allocates two
// length-n work vectors and touches each row twice.
//
// The rhs and the solution are a GENERIC type: the coefficients are scalars while the
// unknowns may be multivectors (the derivative of a multivector field at every node).
// Only T + T and double * T are required of it -- which is why the elimination below
// adds (-a[i]) * x rather than subtracting, T needing no operator-.
//
// NO PIVOTING is done. That is what makes the algorithm O(n), and it is safe here
// because the systems this serves are DIAGONALLY DOMINANT -- for the compact first
// derivative |2/3| > |1/6| + |1/6|, for the second |5/6| > 2*|1/12|, and a
// discretized diffusion operator is dominant by construction. Dominance also bounds
// the multipliers, so no growth factor appears. A pivot that collapses anyway means
// the system is not of that kind, and throws rather than returning noise.
//
// -------------------------------------------------------------------------------------
// BOUNDARY ROWS -- how the first and last rows are built
// -------------------------------------------------------------------------------------
//
// The caller assembles a, b, c and d; this routine applies no boundary conditions of
// its own. The ONE rule it imposes is the band: row i may reference columns i-1, i and
// i+1 and nothing else. Interior rows of a 3-point stencil satisfy that automatically.
// The first and last rows do not, because a centred stencil has no neighbour to one
// side, so each needs a deliberate choice. The four that come up:
//
// 1.) DIRICHLET, as an identity row (simplest; the boundary node stays an unknown).
//     Give the row the trivial equation u_0 = g_0:
//
//         a[0] = 0      b[0] = 1      c[0] = 0      d[0] = g_0
//
//     Indexing stays 1:1 with the grid and the solution comes back including the
//     boundary values. It also makes that row strictly dominant, which is what braces
//     an interior that is only weakly so. This is what fd_derivative() does, and the
//     route to reach for unless the extra unknown is a problem.
//
// 2.) DIRICHLET, by ELIMINATION (the boundary node is not an unknown at all).
//     Solve for the interior only. Row 1's coefficient a[1] multiplies the KNOWN u_0,
//     so move that product to the right-hand side and zero the coefficient:
//
//         d[1] -= a[1] * g_0        a[1] = 0
//         d[n-2] -= c[n-2] * g_1    c[n-2] = 0
//
//     (with the system then assembled over the interior nodes alone). The minus sign
//     is the whole content of "bringing a known term to the other side": the equation
//     reads sum_k coeff_k u_k = rhs, so a term whose u is known leaves the left side
//     as MINUS its contribution. Both routes give the same interior values -- a test
//     in ga_stencil_test.cpp solves one problem both ways and compares.
//
// 3.) NEUMANN -- a derivative is prescribed, so the row IS a one-sided fd stencil,
//     and this is where the band rule bites. The 2-point one-sided stencil for u'
//     touches columns 0 and 1 and fits as it stands. The 2nd-order one-sided stencil
//     touches columns 0, 1 AND 2 -- one column too many. Restore the band with a
//     single step of Gaussian elimination against row 1, the only other row that
//     touches column 2:
//
//         f = w_2 / c[1]                          [w_2 = the out-of-band weight]
//         b[0] = w_0 - f * a[1]
//         c[0] = w_1 - f * b[1]
//         d[0] = g_0 - f * d[1]                   [g_0 = the prescribed derivative]
//
//     and symmetrically at the far end against row n-2. This is the general answer
//     whenever ONE row overhangs the band by a column, Neumann or not.
//
//     What it does NOT buy is accuracy, and it is worth being clear about that: a
//     boundary row one order below the interior scheme does not generally cost global
//     order. Measured on the advection-diffusion problem in ga_stencil_test.cpp, the
//     folded 2nd-order row and the plain 1st-order 2-point row BOTH converge at rate
//     2.00, and the 2-point row is the more accurate of the two by ~1.6x. So reach
//     for the elimination when the wide row is what you have -- not in the
//     expectation that a higher-order flux row will improve the solution. A ghost
//     node is the other classical answer and needs no elimination, at the cost of one
//     extra unknown.
//
// 4.) IMPLICITLY, by a BC-AWARE STENCIL -- the boundary condition built into the
//     boundary row's stencil rather than stated as a row of its own. A prescribed
//     derivative IS supplied derivative data, which is what make_scheme()'s aux_nodes
//     argument takes (ga_usr_fd.hpp), so the stencil is generated knowing it:
//
//         auto sc = make_scheme(x[0], 2, {x[0], x[1]}, {x[0]}, {x[0]});
//         //                            ^ f samples    ^ lhs   ^ f' supplied here
//
//     The prescribed value g_0 then enters the right-hand side through sc.aux_weights
//     and the u-columns stay on the left:
//
//         u''(x_0) = sc.weights[0] u_0 + sc.weights[1] u_1 + sc.aux_weights[0] g_0
//
//     What this buys is FOOTPRINT: one supplied derivative value is worth one node.
//     u'' from {x_0, x_1} plus a prescribed u'(x_0) is 1st order and touches columns
//     0 and 1; the plain one-sided u'' needs {x_0, x_1, x_2} for the same order and
//     overhangs the band. Supplying the condition is therefore the way to make a
//     boundary row fit the three diagonals WITHOUT the elimination of recipe 3 --
//     and it holds one order up too: 3 f-nodes plus u'(x_0) reach 2nd order where the
//     plain form needs 4. A homogeneous condition (u' = 0, the reflective boundary)
//     is the same construction with g_0 = 0, where the rhs term simply drops out.
//
// The same four apply to any operator assembled this way -- the shape of the boundary
// row is a property of the boundary condition, not of the equation being discretized.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> tridiag_solve(std::vector<double> const& a, std::vector<double> const& b,
                             std::vector<double> const& c, std::vector<T> const& d)
{
    size_t const n = b.size();
    if (n < 2 || a.size() != n || c.size() != n || d.size() != n) {
        throw Solver_error("hd::ga::tridiag_solve: diagonals and rhs must have the "
                           "same length n >= 2.");
    }

    // a pivot is judged against the size of the row it came from, so the test means
    // the same thing whatever the system is scaled in
    auto checked_pivot = [&](double piv, size_t i) {
        double const row = std::abs(a[i]) + std::abs(b[i]) + std::abs(c[i]);
        if (std::abs(piv) <= 1.0e-14 * (row > 0.0 ? row : 1.0)) {
            throw Solver_error("hd::ga::tridiag_solve: vanishing pivot -- the system "
                               "is singular or not diagonally dominant.");
        }
        return piv;
    };

    std::vector<double> cp(n, 0.0); // eliminated super-diagonal
    std::vector<T> dp(n, T{});      // eliminated right-hand side

    // forward elimination -- sweep down, removing the sub-diagonal
    double piv = checked_pivot(b[0], 0);
    cp[0] = c[0] / piv;
    dp[0] = (1.0 / piv) * d[0];
    for (size_t i = 1; i < n; ++i) {
        piv = checked_pivot(b[i] - a[i] * cp[i - 1], i);
        cp[i] = c[i] / piv;
        dp[i] = (1.0 / piv) * (d[i] + (-a[i]) * dp[i - 1]);
    }

    // back substitution -- sweep up
    std::vector<T> x(n, T{});
    x[n - 1] = dp[n - 1];
    for (size_t i = n - 1; i-- > 0;) {
        x[i] = dp[i] + (-cp[i]) * x[i + 1];
    }
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
