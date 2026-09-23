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

#include <algorithm> // std::min, std::max, std::swap
#include <cmath>     // std::abs
#include <mdspan>    // std::mdspan, std::dextents, std::extents
#include <stdexcept> // std::runtime_error, std::invalid_argument
#include <string>    // std::string
#include <utility>   // std::move
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
// Guarded dense solve: lu_decomp + lu_backsubs, but REFUSING a numerically rank-deficient
// matrix instead of returning a silently inflated answer.
//
// lu_decomp substitutes TINY (1e-20) for a vanishing pivot, so a singular system does not
// raise -- it produces a solution scaled by that pivot's reciprocal. That is the right
// behaviour for a solver that must never abort mid-iteration, and the wrong behaviour for
// a caller that believes it asked for a pseudo-inverse. This wrapper compares the
// smallest and largest |U_ii| after factorization (a reciprocal condition estimate) and
// throws when the ratio falls below `rcond_min`.
//
// Used by lstsq_solve for the normal-equation solves; not part of the public solver API.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> lu_solve_guarded(std::vector<T> const& A_in, std::vector<T> const& b_in,
                                size_t n, double rcond_min, char const* what)
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

    // reciprocal condition estimate from the U diagonal
    double umin = std::abs(a[0]), umax = std::abs(a[0]);
    for (size_t i = 1; i < n; ++i) {
        double const u = std::abs(a[i * n + i]);
        if (u < umin) umin = u;
        if (u > umax) umax = u;
    }
    if (umax == 0.0 || umin <= rcond_min * umax) {
        throw Solver_error(std::string("hd::ga::") + what +
                           ": matrix is numerically rank-deficient (reciprocal condition "
                           "estimate below tolerance); the normal-equation route cannot "
                           "form a pseudo-inverse here -- pass a damping factor, or "
                           "reformulate.");
    }

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


// The factorization primitives live in hd::ga::detail: they are the machinery behind
// matrix_rank / minnorm_solve / nullspace_basis, not API (and a struct of public
// fields would otherwise be picked up by the binding generator).
namespace detail {

/////////////////////////////////////////////////////////////////////////////////////////
// Rank-revealing QR: Householder reflections with column pivoting on an m x n matrix
// (flat ROW-MAJOR, double). On return `qr` holds R on and above the diagonal and the
// Householder vectors below it (unit leading entry implicit), `tau` the reflector
// scalars, `perm` the column permutation (A P = Q R, column k of A P is column perm[k]
// of A), and `rank` the numerical rank: the first k with |R_kk| <= rtol * |R_00|,
// R_00 being the largest-norm column's pivot. Column pivoting is what makes the
// diagonal of R reveal the rank -- it moves the largest remaining column to the front
// at every step, so a negligible pivot means every remaining column is negligible.
//
// This is the factorization the library reaches for wherever the rank of a map is in
// question: a constraint Jacobian at a singular configuration (a knee lock, a contact
// transition -- Featherstone's varying rank), a redundant task, a screw system's
// dimension. It replaces "throw on rank deficiency" with a rank DETERMINATION and a
// minimum-norm answer (see lstsq_solve). Householder rather than Gram-Schmidt because it
// is backward stable; QR rather than SVD because rank, a null-space basis and a
// minimum-norm solve need no singular values, and the SVD's extra work would buy only
// the manipulability ellipsoid, which belongs to a separate step.
/////////////////////////////////////////////////////////////////////////////////////////
struct qr_factor {
    std::vector<double> qr;   // m x n row-major: R above the diagonal, reflectors below
    std::vector<double> tau;  // reflector scalars, length min(m, n)
    std::vector<size_t> perm; // column permutation, length n
    size_t m{0}, n{0};
    size_t rank{0};
};

inline qr_factor qr_decomp(std::vector<double> A, size_t m, size_t n,
                           double rtol = 1.0e-12)
{
    if (A.size() != m * n) {
        throw Solver_error("hd::ga::qr_decomp: matrix size does not match m x n.");
    }
    qr_factor f;
    f.m = m;
    f.n = n;
    f.qr = std::move(A);
    size_t const kmax = std::min(m, n);
    f.tau.assign(kmax, 0.0);
    f.perm.resize(n);
    for (size_t j = 0; j < n; ++j)
        f.perm[j] = j;
    auto& a = f.qr;

    // squared column norms of the remaining submatrix, for the pivot choice
    std::vector<double> cn(n, 0.0);
    for (size_t j = 0; j < n; ++j)
        for (size_t i = 0; i < m; ++i)
            cn[j] += a[i * n + j] * a[i * n + j];

    double r00 = 0.0;
    f.rank = kmax;
    for (size_t k = 0; k < kmax; ++k) {
        // pivot: the column of largest remaining norm (recomputed from scratch -- these
        // matrices are small, and downdating norms is where pivoted QR goes wrong)
        size_t p = k;
        double best = -1.0;
        for (size_t j = k; j < n; ++j) {
            double s = 0.0;
            for (size_t i = k; i < m; ++i)
                s += a[i * n + j] * a[i * n + j];
            cn[j] = s;
            if (s > best) {
                best = s;
                p = j;
            }
        }
        if (p != k) {
            for (size_t i = 0; i < m; ++i)
                std::swap(a[i * n + k], a[i * n + p]);
            std::swap(f.perm[k], f.perm[p]);
            std::swap(cn[k], cn[p]);
        }
        // Householder reflector for column k, rows k..m-1
        double alpha = std::sqrt(best);
        if (k == 0) r00 = alpha;
        if (alpha <= rtol * r00) {
            // every remaining column is negligible: the rank is k
            f.rank = k;
            for (size_t kk = k; kk < kmax; ++kk)
                f.tau[kk] = 0.0;
            break;
        }
        if (a[k * n + k] > 0.0) alpha = -alpha; // v = x - alpha e_1 without cancellation
        double const v0 = a[k * n + k] - alpha;
        // scale the reflector so its leading entry is 1: v = (1, a[k+1..]/v0)
        for (size_t i = k + 1; i < m; ++i)
            a[i * n + k] /= v0;
        f.tau[k] = -v0 / alpha; // H = I - tau v v^T with v_0 = 1
        a[k * n + k] = alpha;   // R_kk
        // apply H to the remaining columns
        for (size_t j = k + 1; j < n; ++j) {
            double s = a[k * n + j];
            for (size_t i = k + 1; i < m; ++i)
                s += a[i * n + k] * a[i * n + j];
            s *= f.tau[k];
            a[k * n + j] -= s;
            for (size_t i = k + 1; i < m; ++i)
                a[i * n + j] -= s * a[i * n + k];
        }
    }
    return f;
}

// apply Q^T (transpose = true) or Q (false) of a factorization to a length-m vector
inline void qr_apply_q(qr_factor const& f, std::vector<double>& v, bool transpose)
{
    size_t const m = f.m, n = f.n;
    size_t const kmax = std::min(m, n);
    auto reflect = [&](size_t k) {
        if (f.tau[k] == 0.0) return;
        double s = v[k];
        for (size_t i = k + 1; i < m; ++i)
            s += f.qr[i * n + k] * v[i];
        s *= f.tau[k];
        v[k] -= s;
        for (size_t i = k + 1; i < m; ++i)
            v[i] -= s * f.qr[i * n + k];
    };
    if (transpose) { // Q^T = H_{kmax-1} ... H_0
        for (size_t k = 0; k < kmax; ++k)
            reflect(k);
    }
    else { // Q = H_0 ... H_{kmax-1} (each H is its own inverse)
        for (size_t k = kmax; k-- > 0;)
            reflect(k);
    }
}

} // namespace detail


/////////////////////////////////////////////////////////////////////////////////////////
// Numerical rank of an m x ncols matrix (flat ROW-MAJOR): the rank revealed by the
// pivoted QR at relative tolerance rtol. The library's answer to "how many independent
// constraints / how much mobility here" -- Featherstone's r in mobility = n - r.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
size_t matrix_rank(std::vector<T> const& A, size_t rows, size_t ncols,
                   double rtol = 1.0e-12)
{
    if (rows == 0 || ncols == 0) return 0;
    std::vector<double> a(rows * ncols);
    for (size_t i = 0; i < rows * ncols; ++i)
        a[i] = static_cast<double>(A[i]);
    return detail::qr_decomp(std::move(a), rows, ncols, rtol).rank;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Minimum-norm least-squares solve of A x = b for ANY shape and rank, through a complete
// orthogonal decomposition: pivoted QR of A gives A P = Q [R11 R12; 0 0] with R11 the
// r x r leading block (r = rank); a second, unpivoted QR of [R11 R12]^T gives
// [R11 R12] = L Z^T with L lower-triangular r x r and Z orthonormal n x r, so
//
//     A P = Q [L Z^T; 0]      x = P Z L^-1 (Q^T b)_{1:r}
//
// is the solution of least residual and, among those, of least norm -- the Moore-Penrose
// answer A^+ b, at every rank. Where r = n the second QR is skipped (the solution is
// then plain back-substitution). This is what makes a mechanism survive a singular
// configuration: the solve returns a finite, continuous answer instead of throwing or
// inflating. Returns x (length ncols); writes the rank if `rank_out` is non-null.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> minnorm_solve(std::vector<T> const& A, std::vector<T> const& b,
                             size_t ncols, size_t* rank_out = nullptr,
                             double rtol = 1.0e-12)
{
    size_t const m = b.size();
    if (A.size() != m * ncols) {
        throw Solver_error(
            "hd::ga::minnorm_solve: A must have b.size() * ncols entries.");
    }
    std::vector<double> a(m * ncols), bb(m);
    for (size_t i = 0; i < m * ncols; ++i)
        a[i] = static_cast<double>(A[i]);
    for (size_t i = 0; i < m; ++i)
        bb[i] = static_cast<double>(b[i]);

    detail::qr_factor const f = detail::qr_decomp(std::move(a), m, ncols, rtol);
    size_t const r = f.rank;
    if (rank_out) *rank_out = r;
    std::vector<T> x(ncols, T(0));
    if (r == 0) return x;

    detail::qr_apply_q(f, bb, true);   // c = Q^T b; only the first r entries matter
    std::vector<double> y(ncols, 0.0); // the solution in the permuted column order

    if (r == ncols) {
        // full column rank: back-substitute R11 y = c
        for (size_t i = r; i-- > 0;) {
            double s = bb[i];
            for (size_t j = i + 1; j < r; ++j)
                s -= f.qr[i * ncols + j] * y[j];
            y[i] = s / f.qr[i * ncols + i];
        }
    }
    else {
        // rank-deficient / underdetermined: T = [R11 R12] (r x ncols); factor T^T = Z L^T
        std::vector<double> Tt(ncols * r, 0.0); // ncols x r row-major
        for (size_t i = 0; i < r; ++i)
            for (size_t j = i; j < ncols; ++j)
                Tt[j * r + i] = f.qr[i * ncols + j];
        detail::qr_factor const g = detail::qr_decomp(std::move(Tt), ncols, r, 0.0);
        // T has full row rank r, so the threshold 0 keeps every column -- but qr_decomp
        // still PIVOTS them: T^T P = Q R, column i of T^T P being column perm[i] of T^T.
        // Hence T = P R^T Q^T, and T y = c becomes R^T w = P^T c with w = (Q^T y)_{1:r},
        // where (P^T c)_i = c_{perm[i]}: the right-hand side is read in the PIVOTED order
        // (forward substitution on the transposed upper factor). The inverse reading,
        // c_{perm[i]} = bb_i, agrees only when the permutation is its own inverse (the
        // identity, one swap) and otherwise returns a y that does not solve T y = c.
        std::vector<double> c(r);
        for (size_t i = 0; i < r; ++i)
            c[i] = bb[g.perm[i]];
        std::vector<double> w(r, 0.0);
        for (size_t i = 0; i < r; ++i) {
            double s = c[i];
            for (size_t j = 0; j < i; ++j)
                s -= g.qr[j * r + i] * w[j];
            w[i] = s / g.qr[i * r + i];
        }
        // y = Z w = Q_g [w; 0]  (Q_g is ncols x ncols, Z its first r columns)
        std::vector<double> yy(ncols, 0.0);
        for (size_t i = 0; i < r; ++i)
            yy[i] = w[i];
        detail::qr_apply_q(g, yy, false);
        y = yy;
    }
    // undo the column permutation of the first factorization
    for (size_t k = 0; k < ncols; ++k)
        x[f.perm[k]] = static_cast<T>(y[k]);
    return x;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Orthonormal basis of the null space of an m x ncols matrix A: the ncols - rank
// vectors z with A z = 0, returned as an ncols x (ncols - rank) ROW-MAJOR matrix (one
// basis vector per column). From the complete orthogonal decomposition above, they are
// the trailing columns of Z, i.e. Q of the second factorization applied to the unit
// vectors e_r .. e_{ncols-1}, permuted back. Empty for full column rank. This is the
// redundancy of a task (the motions a Jacobian cannot see) and, for a screw system, the
// reciprocal system's dimension count.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> nullspace_basis(std::vector<T> const& A, size_t rows, size_t ncols,
                               size_t* rank_out = nullptr, double rtol = 1.0e-12)
{
    std::vector<double> a(rows * ncols);
    for (size_t i = 0; i < rows * ncols; ++i)
        a[i] = static_cast<double>(A[i]);
    detail::qr_factor const f = detail::qr_decomp(std::move(a), rows, ncols, rtol);
    size_t const r = f.rank;
    if (rank_out) *rank_out = r;
    size_t const k = ncols - r;
    std::vector<T> N(ncols * k, T(0));
    if (k == 0) return N;
    // T = [R11 R12] (r x ncols), factor T^T; the null space of T is spanned by the last
    // ncols - r columns of Q_g (r == 0: every unit vector)
    std::vector<double> Tt(ncols * std::max<size_t>(r, 1), 0.0);
    for (size_t i = 0; i < r; ++i)
        for (size_t j = i; j < ncols; ++j)
            Tt[j * r + i] = f.qr[i * ncols + j];
    for (size_t c = 0; c < k; ++c) {
        std::vector<double> e(ncols, 0.0);
        e[r + c] = 1.0;
        if (r > 0) {
            detail::qr_factor const g = detail::qr_decomp(Tt, ncols, r, 0.0);
            detail::qr_apply_q(g, e, false);
        }
        for (size_t j = 0; j < ncols; ++j)
            N[f.perm[j] * k + c] = static_cast<T>(e[j]);
    }
    return N;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Least-squares / minimum-norm dense solve of a (possibly non-square) system A x = b,
// where A is m x ncols (flat ROW-MAJOR, m = b.size()). Returns x (length ncols).
//
// Without damping this IS the Moore-Penrose solution A^+ b at every shape and rank,
// through the rank-revealing QR (minnorm_solve): square, underdetermined
// (minimum-norm), overdetermined (least-squares) and rank-deficient inputs are one
// case. A rank-deficient input no longer throws -- it gets the finite, continuous
// minimum-norm answer, which is what a mechanism passing through a singular
// configuration needs; ask matrix_rank() when the rank itself is the question.
//
// `damping` (default 0 = off) selects the OTHER route: Tikhonov regularization of the
// normal equations, adding lambda^2 = damping * trace(Gram)/rows to the diagonal before
// the shared LU. That trades exactness for a damped, continuous answer -- the standard
// damped-least-squares / Levenberg step -- and is kept because it is a different
// estimator, not a worse implementation of the same one: near a singularity the damped
// step stays bounded where the pseudo-inverse step (finite, but the inverse of a small
// pivot) grows. The scaling is relative to the Gram trace so one damping value works
// across problem scales. Ignored for the square regime, which forms no Gram matrix.
//
// Domain- and dimension-agnostic (pure linear algebra over T): used by the closed-loop
// constraint solver for the position Newton step (A = constraint Jacobian G, b = -g), the
// velocity / acceleration distribution, and the velocity projection -- in both 2D and 3D,
// and available to any other caller (it carries no GA or physics knowledge).
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> lstsq_solve(std::vector<T> const& A, std::vector<T> const& b, size_t ncols,
                           double damping = 0.0)
{
    size_t const m = b.size();
    if (damping <= 0.0) return minnorm_solve(A, b, ncols);

    // add lambda^2 = damping * trace(G)/n to the diagonal of an n x n Gram matrix
    auto const regularize = [damping](std::vector<T>& G, size_t n) {
        T tr = T(0);
        for (size_t i = 0; i < n; ++i)
            tr += G[i * n + i];
        T const lam2 = T(damping) * tr / T(n);
        for (size_t i = 0; i < n; ++i)
            G[i * n + i] += lam2;
    };
    double const rcond_min = 0.0; // regularizing IS the answer to rank deficiency

    if (ncols == m) {
        // square A: no Gram matrix is formed, so damping does not apply
        return minnorm_solve(A, b, ncols);
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
        regularize(AAt, m);
        std::vector<T> const y = lu_solve_guarded(AAt, b, m, rcond_min, "lstsq_solve");
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
    regularize(AtA, ncols);
    return lu_solve_guarded(AtA, Atb, ncols, rcond_min, "lstsq_solve");
}


/////////////////////////////////////////////////////////////////////////////////////////
// BOUNDED-VARIABLE least squares: the minimiser of ||A x - b|| subject to a BOX on the
// decision variables, lo <= x <= hi. A (m x ncols, flat ROW-MAJOR), b (length m), and
// lo / hi (length ncols; +-infinity where a variable is unbounded, which is the neutral
// value and reduces the answer to lstsq_solve's).
//
// WHY THIS AND NOT A CLAMP. Solving unconstrained and then cutting each variable to its
// bound is one comparison per variable and it answers a DIFFERENT question: the cut
// vector is no longer the least-squares solution of anything, so whatever the residual
// was arranged to express is silently lost. This returns the constrained optimum
// instead, and for a BOX it is not an approximation of the quadratic program -- it IS
// the quadratic program's answer, reached in finitely many steps.
//
// HOW (the classical active-set method, Stark & Parker): hold each variable either FREE
// or AT a bound; solve the free columns alone against the residual the bound-held ones
// leave (dropping columns -- which is why minnorm_solve taking the matrix by column
// count is all that was needed); if a free variable lands outside its box, walk along
// the segment to the first bound it crosses and hold it there; when every free variable
// is inside, read the gradient g = A^T (A x - b) and RELEASE any held variable whose
// gradient says the residual would fall if it came off its bound. Each step either
// holds or releases one variable and strictly decreases the residual, so it terminates;
// `max_iter` is a guard against a pathological cycle, not part of the method, and
// `iters_out` reports what it actually took.
//
// WHAT IT DOES NOT DO, so the boundary is known in advance: the box is on the VARIABLES.
// A constraint that couples them -- a friction cone, a centre-of-pressure polygon, a
// bound on some A x rather than on x -- is a general polytope, where an active-set
// iteration can cycle and is not guaranteed to land on the optimum. That is where a
// real quadratic program earns its place.
//
// Domain- and dimension-agnostic (pure linear algebra over T), like the rest of this
// file: it carries no GA and no physics knowledge.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> bvls_solve(std::vector<T> const& A, std::vector<T> const& b, size_t ncols,
                          std::vector<T> const& lo, std::vector<T> const& hi,
                          size_t* iters_out = nullptr, size_t max_iter = 0)
{
    size_t const m = b.size();
    if (A.size() != m * ncols) {
        throw Solver_error("hd::ga::bvls_solve: A must have b.size() * ncols entries.");
    }
    if (lo.size() != ncols || hi.size() != ncols) {
        throw Solver_error("hd::ga::bvls_solve: lo and hi must have ncols entries.");
    }
    for (size_t j = 0; j < ncols; ++j)
        if (!(lo[j] <= hi[j]))
            throw Solver_error("hd::ga::bvls_solve: every bound must have lo <= hi.");

    if (max_iter == 0) max_iter = 3 * ncols + 10;

    enum class state { free_var, at_lo, at_hi };
    std::vector<state> st(ncols, state::free_var);
    std::vector<T> x(ncols, T(0));

    // start feasible: a variable whose box excludes 0 begins held at the near bound
    for (size_t j = 0; j < ncols; ++j) {
        if (lo[j] > T(0)) {
            x[j] = lo[j];
            st[j] = state::at_lo;
        }
        else if (hi[j] < T(0)) {
            x[j] = hi[j];
            st[j] = state::at_hi;
        }
    }

    auto residual = [&](std::vector<T> const& v) { // r = A v - b
        std::vector<T> r(m);
        for (size_t i = 0; i < m; ++i) {
            T acc = -b[i];
            for (size_t j = 0; j < ncols; ++j)
                acc += A[i * ncols + j] * v[j];
            r[i] = acc;
        }
        return r;
    };

    size_t used = 0;
    for (; used < max_iter; ++used) {
        // the free columns, and what the held ones leave of b
        std::vector<size_t> fr;
        for (size_t j = 0; j < ncols; ++j)
            if (st[j] == state::free_var) fr.push_back(j);

        std::vector<T> cand = x;
        if (!fr.empty()) {
            std::vector<T> Af(m * fr.size()), bf(m);
            for (size_t i = 0; i < m; ++i) {
                T held = T(0);
                for (size_t j = 0; j < ncols; ++j)
                    if (st[j] != state::free_var) held += A[i * ncols + j] * x[j];
                bf[i] = b[i] - held;
                for (size_t c = 0; c < fr.size(); ++c)
                    Af[i * fr.size() + c] = A[i * ncols + fr[c]];
            }
            std::vector<T> const xf = minnorm_solve(Af, bf, fr.size());
            for (size_t c = 0; c < fr.size(); ++c)
                cand[fr[c]] = xf[c];
        }

        // does the free solve stay inside its box?
        T step = T(1);
        size_t hit = ncols;
        bool to_lo = false;
        for (size_t const j : fr) {
            if (cand[j] < lo[j] || cand[j] > hi[j]) {
                T const target = (cand[j] < lo[j]) ? lo[j] : hi[j];
                T const denom = cand[j] - x[j];
                T const t = (denom == T(0)) ? T(0) : (target - x[j]) / denom;
                if (t < step) {
                    step = t;
                    hit = j;
                    to_lo = (cand[j] < lo[j]);
                }
            }
        }
        if (hit != ncols) { // walk to the first bound crossed and hold that variable
            for (size_t j = 0; j < ncols; ++j)
                x[j] += step * (cand[j] - x[j]);
            x[hit] = to_lo ? lo[hit] : hi[hit];
            st[hit] = to_lo ? state::at_lo : state::at_hi;
            continue;
        }
        x = cand;

        // every free variable is inside: release a held one if the gradient wants it
        // -- but only if it wants it by more than rounding. Without that floor a
        // gradient of order the residual's own error releases a variable the next solve
        // puts straight back, and the iteration walks to its guard instead of stopping
        // (measured: every system in the seeded sweep hit `max_iter`, while still
        // returning the optimum, before this floor was added).
        std::vector<T> const r = residual(x);
        T scale = T(0); // the gradient's own magnitude, so the floor is relative
        for (size_t i = 0; i < m; ++i)
            scale = std::max(scale, std::abs(r[i]));
        T gscale = T(0);
        for (size_t i = 0; i < m * ncols; ++i)
            gscale = std::max(gscale, std::abs(A[i]));
        T const gtol = T(1.0e-12) * std::max(T(1), scale * gscale * T(m));
        size_t best = ncols;
        T best_g = gtol;
        for (size_t j = 0; j < ncols; ++j) {
            if (st[j] == state::free_var) continue;
            T g = T(0); // g_j = (A^T r)_j
            for (size_t i = 0; i < m; ++i)
                g += A[i * ncols + j] * r[i];
            // at the LOWER bound x may only rise, which helps when g < 0; at the UPPER
            // bound it may only fall, which helps when g > 0
            T const gain = (st[j] == state::at_lo) ? -g : g;
            if (gain > best_g) {
                best_g = gain;
                best = j;
            }
        }
        if (best == ncols) break; // optimal: no held variable wants to move
        st[best] = state::free_var;
    }
    if (iters_out) *iters_out = used;
    return x;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Least-squares QUADRATIC PROGRAM with equalities, general inequalities and a box:
//
//     minimize  1/2 |A x - b|^2   subject to   E x = e,   C x >= d,   lo <= x <= hi ,
//
// A (p x ncols), E (q x ncols) and C (r x ncols) flat ROW-MAJOR, b / e / d their
// right-hand sides (q or r may be zero: pass empty vectors), lo / hi of length ncols
// (+-infinity where a variable is unbounded). `x` holds a FEASIBLE starting point on
// entry -- every constraint satisfied, to `tol` -- and the answer on return. A caller
// always has one: the box's point nearest zero when E and C are empty, or the previous
// solution of a hierarchy whose next level only adds constraints that solution meets.
// Throws Solver_error if the start is not feasible.
//
// WHY THIS BESIDE bvls_solve. That one's box is on the VARIABLES, and its active set is
// searched one variable at a time. A constraint that COUPLES the variables -- a friction
// cone, a centre-of-pressure polygon, a bound on some rows of C x, a limit stated in one
// set of coordinates while the objective lives in another -- is a general polytope, the
// case bvls_solve's own comment names as where a real quadratic program earns its place.
// This is that program, sized for the small dense problems this library poses (tens of
// variables, tens of rows). With E and C empty it is bvls_solve's problem and returns
// its answer; with the box loose and C empty it is the equality-constrained least
// squares, whose solution satisfies the KKT conditions A^T (A x - b) = E^T lambda.
//
// HOW: the primal active-set method with a null-space step (Nocedal & Wright, ch. 16).
// A WORKING SET W holds constraints as equalities: every row of E, plus the inequality
// rows and box faces currently active. From the feasible x the step to the minimizer on
// W's affine set is
//
//     p = N z,   N = nullspace_basis(E_W),   z = argmin |(A N) z - (b - A x)| ,
//
// with z the minimum-norm solution, so a rank-deficient objective picks one minimizer
// rather than failing. If p is blocked by an inactive inequality, x moves to the first
// one hit and that row joins W; if p is zero, the multipliers
// lambda = argmin |E_W^T lambda - A^T (A x - b)| are read: an inequality whose
// multiplier is negative is holding the objective up and is RELEASED (the most negative
// first); when none is, the KKT conditions hold and x is the constrained optimum --
// global, since the objective is convex. Every iterate is feasible. A blocking
// constraint met at zero step length is added without moving, so the working set
// grows and a degenerate vertex is left in finitely many steps; `max_iter` (default
// 20 * (ncols + q + r + 1)) guards a pathological cycle, and `iters_out` reports what
// it took.
//
// WHAT IT DOES NOT DO: it is not a general QP with an arbitrary Hessian (the objective
// is a least-squares residual, which is what every use here has), it does not scale
// the data (the caller's rows are already in comparable units), and its tolerance is
// one number, `tol`, used relative to the size of x for "zero step", "active" and
// "negative multiplier" alike.
//
// Domain- and dimension-agnostic, like the rest of this file.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T>&
qp_ls_solve(std::vector<T> const& A, std::vector<T> const& b, size_t ncols,
            std::vector<T> const& E, std::vector<T> const& e, std::vector<T> const& C,
            std::vector<T> const& d, std::vector<T> const& lo, std::vector<T> const& hi,
            std::vector<T>& x, size_t* iters_out = nullptr, size_t max_iter = 0,
            double tol = 1.0e-10)
{
    size_t const n = ncols, p = b.size(), q = e.size(), r = d.size();
    if (A.size() != p * n || E.size() != q * n || C.size() != r * n) {
        throw Solver_error("hd::ga::qp_ls_solve: A, E and C must each have rows * ncols "
                           "entries.");
    }
    if (lo.size() != n || hi.size() != n || x.size() != n) {
        throw Solver_error("hd::ga::qp_ls_solve: lo, hi and x must have ncols entries.");
    }
    T const inf = std::numeric_limits<T>::infinity();
    if (max_iter == 0) max_iter = 20 * (n + q + r + 1);

    // every inequality as a "≥" row a^T x >= rhs: the r rows of C, then the lower and
    // the upper face of each bounded variable
    struct ineq {
        std::vector<T> a;
        T rhs;
    };
    std::vector<ineq> in;
    for (size_t i = 0; i < r; ++i)
        in.push_back({std::vector<T>(C.begin() + std::ptrdiff_t(i * n),
                                     C.begin() + std::ptrdiff_t((i + 1) * n)),
                      d[i]});
    for (size_t j = 0; j < n; ++j) {
        if (lo[j] > -inf) {
            std::vector<T> a(n, T(0));
            a[j] = T(1);
            in.push_back({a, lo[j]});
        }
        if (hi[j] < inf) {
            std::vector<T> a(n, T(0));
            a[j] = T(-1);
            in.push_back({a, -hi[j]});
        }
    }
    size_t const ni = in.size();
    auto scale = [&]() {
        T s = T(1);
        for (size_t j = 0; j < n; ++j)
            s = std::max(s, std::abs(x[j]));
        return s;
    };
    auto dotx = [&](std::vector<T> const& a, std::vector<T> const& v) {
        T s = T(0);
        for (size_t j = 0; j < n; ++j)
            s += a[j] * v[j];
        return s;
    };
    // the start must be feasible
    {
        T const feas = T(tol) * scale();
        for (size_t i = 0; i < q; ++i) {
            T s = -e[i];
            for (size_t j = 0; j < n; ++j)
                s += E[i * n + j] * x[j];
            if (std::abs(s) > feas * T(100))
                throw Solver_error("hd::ga::qp_ls_solve: the starting point violates an "
                                   "equality constraint.");
        }
        for (auto const& c : in)
            if (dotx(c.a, x) < c.rhs - feas * T(100))
                throw Solver_error("hd::ga::qp_ls_solve: the starting point violates an "
                                   "inequality constraint or the box.");
    }
    std::vector<bool> active(ni, false);
    for (size_t i = 0; i < ni; ++i)
        active[i] = std::abs(dotx(in[i].a, x) - in[i].rhs) <= T(tol) * scale();

    size_t used = 0;
    for (; used < max_iter; ++used) {
        // the working set as an equality block
        std::vector<size_t> wi;
        for (size_t i = 0; i < ni; ++i)
            if (active[i]) wi.push_back(i);
        size_t const kw = q + wi.size();
        std::vector<T> Ew(kw * n, T(0));
        for (size_t i = 0; i < q; ++i)
            for (size_t j = 0; j < n; ++j)
                Ew[i * n + j] = E[i * n + j];
        for (size_t w = 0; w < wi.size(); ++w)
            for (size_t j = 0; j < n; ++j)
                Ew[(q + w) * n + j] = in[wi[w]].a[j];

        // the step on W's affine set: p = N z
        std::vector<T> resid(p);
        for (size_t i = 0; i < p; ++i) {
            T s = b[i];
            for (size_t j = 0; j < n; ++j)
                s -= A[i * n + j] * x[j];
            resid[i] = s;
        }
        std::vector<T> step(n, T(0));
        size_t k = n;
        std::vector<T> N;
        if (kw > 0) {
            size_t rk = 0;
            N = nullspace_basis(Ew, kw, n, &rk);
            k = n - rk;
        }
        if (k > 0) {
            std::vector<T> AN(p * k, T(0));
            for (size_t i = 0; i < p; ++i)
                for (size_t c = 0; c < k; ++c) {
                    T s = T(0);
                    for (size_t j = 0; j < n; ++j)
                        s += A[i * n + j] *
                             (kw > 0 ? N[j * k + c] : (j == c ? T(1) : T(0)));
                    AN[i * k + c] = s;
                }
            auto const z = lstsq_solve(AN, resid, k);
            for (size_t j = 0; j < n; ++j) {
                T s = T(0);
                for (size_t c = 0; c < k; ++c)
                    s += (kw > 0 ? N[j * k + c] : (j == c ? T(1) : T(0))) * z[c];
                step[j] = s;
            }
        }
        T pn = T(0);
        for (size_t j = 0; j < n; ++j)
            pn = std::max(pn, std::abs(step[j]));

        if (pn <= T(tol) * scale()) {
            // stationary on W: the multipliers decide whether to release or to stop
            if (wi.empty()) break;
            std::vector<T> g(n, T(0)); // A^T (A x - b) = -A^T resid
            for (size_t j = 0; j < n; ++j) {
                T s = T(0);
                for (size_t i = 0; i < p; ++i)
                    s -= A[i * n + j] * resid[i];
                g[j] = s;
            }
            std::vector<T> EwT(n * kw);
            for (size_t i = 0; i < kw; ++i)
                for (size_t j = 0; j < n; ++j)
                    EwT[j * kw + i] = Ew[i * n + j];
            auto const lambda = lstsq_solve(EwT, g, kw);
            T gscale = T(tol);
            for (size_t j = 0; j < n; ++j)
                gscale = std::max(gscale, std::abs(g[j]));
            size_t worst = ni;
            T most = -T(tol) * gscale;
            for (size_t w = 0; w < wi.size(); ++w)
                if (lambda[q + w] < most) {
                    most = lambda[q + w];
                    worst = wi[w];
                }
            if (worst == ni) break; // KKT: every active inequality pushes the right way
            active[worst] = false;
            continue;
        }

        // the ratio test over the inactive inequalities
        T alpha = T(1);
        size_t block = ni;
        for (size_t i = 0; i < ni; ++i) {
            if (active[i]) continue;
            T const ap = dotx(in[i].a, step);
            if (ap >= -T(tol)) continue; // the step does not approach this one
            T const room = dotx(in[i].a, x) - in[i].rhs; // >= 0 by feasibility
            T const a_i = std::max(T(0), room) / (-ap);
            if (a_i < alpha) {
                alpha = a_i;
                block = i;
            }
        }
        for (size_t j = 0; j < n; ++j)
            x[j] += alpha * step[j];
        if (block < ni) {
            active[block] = true;
            // land exactly on the face, so the working-set solves see it as tight
            T const viol = dotx(in[block].a, x) - in[block].rhs;
            if (std::abs(viol) > T(0)) {
                T aa = T(0);
                for (size_t j = 0; j < n; ++j)
                    aa += in[block].a[j] * in[block].a[j];
                if (aa > T(0))
                    for (size_t j = 0; j < n; ++j)
                        x[j] -= viol * in[block].a[j] / aa;
            }
        }
    }
    if (iters_out) *iters_out = used;
    return x;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Null-space projection of a secondary objective: given an m x ncols map A and a desired
// rate v (length ncols), return the part of v that A cannot see,
//
//     v_ns = (I - A^+ A) v = v - A^+ (A v) ,
//
// so that A v_ns == 0 to the accuracy of the solve. This is the standard redundancy
// resolution step: a redundant system's primary task fixes only m of its ncols freedoms,
// and the remaining ncols - m can serve a secondary objective WITHOUT disturbing the
// task -- posture centering, joint-limit avoidance, obstacle repulsion.
//
// `damping` is forwarded to lstsq_solve; pass a small positive value when A may approach
// rank deficiency, which is exactly where an unprojected secondary objective would
// otherwise be amplified.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> nullspace_project(std::vector<T> const& A, std::vector<T> const& v,
                                 size_t rows, size_t ncols, double damping = 0.0)
{
    if (v.size() != ncols) {
        throw Solver_error("hd::ga::nullspace_project: v must have ncols entries.");
    }
    std::vector<T> Av(rows, T(0)); // the task component of v
    for (size_t i = 0; i < rows; ++i) {
        T s = T(0);
        for (size_t k = 0; k < ncols; ++k)
            s += A[i * ncols + k] * v[k];
        Av[i] = s;
    }
    std::vector<T> const corr = lstsq_solve(A, Av, ncols, damping); // A^+ (A v)
    std::vector<T> out(ncols);
    for (size_t k = 0; k < ncols; ++k)
        out[k] = v[k] - corr[k];
    return out;
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
// the Lagrange multipliers l (length m) into `lambda_out` if non-null, and the rank of G
// into `rank_out` if non-null.
//
// G MAY BE RANK-DEFICIENT, and it is as a matter of course: the rank of a closed-loop
// constraint matrix varies with the configuration (Featherstone 8.5, 8.10 -- a knee lock,
// a contact transition, an over-constrained loop). The bordered matrix is then singular,
// x stays uniquely determined while l is only determined up to the null space of G^T,
// and the constraint rows may be slightly inconsistent -- so the system is solved in the
// least-squares sense with the minimum-norm l (minnorm_solve). With G of full row rank
// the LU path is taken unchanged. Before this, a singular bordered matrix went through
// the LU's TINY pivot substitution and returned multipliers of order 1e16 without
// complaint (measured at a biped's knee lock).
//
// Domain- and dimension-agnostic (pure linear algebra over T): the closed-loop dynamics
// (2D and 3D) use it with M = joint-space mass matrix, G = constraint Jacobian, f =
// generalised force tau, g = -G-dot q-dot, to get the joint accelerations x = q-ddot and
// the constraint forces l. Carries no GA or physics knowledge.
/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> kkt_solve(std::vector<T> const& M, std::vector<T> const& G,
                         std::vector<T> const& f, std::vector<T> const& g, size_t n,
                         size_t m, std::vector<T>* lambda_out = nullptr,
                         size_t* rank_out = nullptr)
{
    size_t const N = n + m;
    size_t const r = (m > 0) ? matrix_rank(G, m, n) : 0;
    if (rank_out) *rank_out = r;
    std::vector<T> K(N * N, T(0)), rhs(N, T(0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            K[i * N + j] = M[i * n + j]; // M block (top-left)
    for (size_t c = 0; c < m; ++c)
        for (size_t j = 0; j < n; ++j) {
            K[j * N + (n + c)] = G[c * n + j]; // G^T block (top-right)
            K[(n + c) * N + j] = G[c * n + j]; // G  block (bottom-left)
        }
    for (size_t i = 0; i < n; ++i)
        rhs[i] = f[i];
    for (size_t c = 0; c < m; ++c)
        rhs[n + c] = g[c];

    // full row rank: the bordered matrix is regular, the shared LU is exact and cheap
    if (r == m) {
        std::vector<T> const sol = lu_solve(K, rhs, N);
        if (lambda_out) lambda_out->assign(sol.begin() + n, sol.end());
        return std::vector<T>(sol.begin(), sol.begin() + n);
    }

    // A rank drop: some rows of G are dependent (a knee lock, a body resting on two
    // pins, an over-constrained loop). q-ddot is still unique when the independent
    // rows constrain the coordinates consistently, so solve the REDUCED bordered
    // system built from a maximal independent subset of the rows (greedy: a row is
    // kept when it raises the rank) by LU; then recover the multipliers of ALL rows
    // as the minimum-norm solution of Gᵀ λ = f - M q-ddot -- the constraint force is
    // determined, its split over dependent rows is not, and the least-norm split is
    // the convention. (A minimum-norm solve of the whole singular bordered matrix,
    // the earlier route, is WRONG here: it minimises the joint norm of (q-ddot,
    // lambda) and returns accelerations for a system that cannot move.)
    std::vector<size_t> keep;
    std::vector<T> Gk; // the kept rows, r x n
    for (size_t c = 0; c < m && keep.size() < r; ++c) {
        std::vector<T> trial = Gk;
        trial.insert(trial.end(), G.begin() + c * n, G.begin() + (c + 1) * n);
        if (matrix_rank(trial, keep.size() + 1, n) > keep.size()) {
            keep.push_back(c);
            Gk = std::move(trial);
        }
    }
    size_t const rk = keep.size();
    size_t const Nk = n + rk;
    std::vector<T> Kk(Nk * Nk, T(0)), rk_rhs(Nk, T(0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            Kk[i * Nk + j] = M[i * n + j];
    for (size_t c = 0; c < rk; ++c)
        for (size_t j = 0; j < n; ++j) {
            Kk[j * Nk + (n + c)] = Gk[c * n + j];
            Kk[(n + c) * Nk + j] = Gk[c * n + j];
        }
    for (size_t i = 0; i < n; ++i)
        rk_rhs[i] = f[i];
    for (size_t c = 0; c < rk; ++c)
        rk_rhs[n + c] = g[keep[c]];
    std::vector<T> const solk =
        (rk + n > 0) ? lu_solve(Kk, rk_rhs, Nk) : std::vector<T>{};
    std::vector<T> qdd(solk.begin(), solk.begin() + n);
    if (lambda_out) {
        // Gᵀ (n x m) λ = f - M q-ddot, minimum-norm λ over all m rows
        std::vector<T> GT(n * m), rhs_l(n);
        for (size_t c = 0; c < m; ++c)
            for (size_t j = 0; j < n; ++j)
                GT[j * m + c] = G[c * n + j];
        for (size_t i = 0; i < n; ++i) {
            T acc = f[i];
            for (size_t j = 0; j < n; ++j)
                acc -= M[i * n + j] * qdd[j];
            rhs_l[i] = acc;
        }
        *lambda_out = minnorm_solve(GT, rhs_l, m);
    }
    return qdd;
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
