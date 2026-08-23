#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// Finite-difference stencil generator (Fornberg-style, via a small linear solve).
//
// Given a development point x0 and point sets for f, f' and f'', stencil_t computes the
// weights of a finite-difference approximation by solving the Taylor-matching system
// with hd::ga::lu_decomp / lu_backsubs, plus the achieved order and the coefficient of
// the leading truncation term. Both explicit stencils (all derivative information on one
// side, e.g. central differences) and compact/implicit stencils (several f' or f''
// points on the left-hand side, e.g. Pade schemes) are supported.
//
// Usage (weights are available after construction), for the central second
// derivative f''(x0) ~ wf0[0]*f(x0-h) + wf0[1]*f(x0) + wf0[2]*f(x0+h):
//
//     stencil_t s(x0, stencil_lhs::f2, {x0 - h, x0, x0 + h}, {}, {x0});
//     // -> s.wf0 = {1/h^2, -2/h^2, 1/h^2}, s.order = 2, s.trunc_err = h^2/12
//
// The lhs derivative points (xf1 for stencil_lhs::f1, xf2 for stencil_lhs::f2) are part
// of their point sets; the corresponding weights (wf1 / wf2) are normalized so that the
// lhs coefficients sum to 1.
//
// Intended as internal infrastructure for discretized field derivatives (numerical
// nabla on grids, e.g. for the planned STA electrodynamics tier). Closed-form
// derivatives of analytic fields do not need stencils and should not route through
// this header.
//
// The order/truncation-error detection compares each residual against the magnitude of
// the contributions it is made of, so it is independent of how the node set is scaled:
// a small physical spacing and a wide stencil normalized into [-1, 1] are read alike.
// If no leading term is found the stencil is not characterizable and the ctor throws,
// rather than reporting order 0.
//
// Adapted from the hd utility library (hd_stencil.hpp) and made internal to the ga
// library so the field-operator layers carry no external dependency. Three fixes
// relative to the source version: the leading truncation term is the FIRST
// non-vanishing residual order (the source kept overwriting, reporting the last), lhs
// terms enter the residual with the same sign convention used to build the matrix, and
// the order test is relative rather than against a fixed absolute threshold.
/////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>     // std::pow, std::abs
#include <limits>    // std::numeric_limits (fact beyond 170!)
#include <mdspan>    // std::mdspan (views onto the scratch storage)
#include <stdexcept> // std::invalid_argument
#include <vector>    // std::vector

#include "ga_solver.hpp" // lu_decomp, lu_backsubs

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// n! as a double -- the stencil matrix entries are Taylor coefficients 1/i!, so this
// is called from the innermost loops of calc_stencil().
//
// Adapted from the hd utility library (hd_functions.hpp), whose table lookup replaces
// the multiply-n-times loop that used to stand here. Measured over the call pattern
// that matters (runtime arguments, n < 20):
//
//     plain loop (was here)        3.89 ns/call      rel err 5e-16 at n = 100
//     hd fact(), lazy table        0.60 ns/call      rel err 1e-10 at n = 100
//     this, compile-time table     0.50 ns/call      rel err 5e-16 at n = 100
//
// Two deliberate departures from the hd version, both behind those numbers:
//
//   - the table is built at COMPILE time and never written to. hd fills its table
//     lazily into a function-local static, which is a data race the moment two
//     threads build a stencil at once. A constexpr table also keeps fact() usable in
//     a constant expression.
//
//   - no log_gamma branch. hd switches to exp(log_gamma(n+1)) above n = 32, which
//     costs five significant digits (1e-10 against 5e-16) because exp amplifies the
//     error of its argument. 171 entries cover every factorial a double can hold, so
//     the branch has nothing left to do: 171! is infinity in double, and returning
//     that is more honest than returning an approximation of it.
/////////////////////////////////////////////////////////////////////////////////////////

namespace detail {

// 170! is the last factorial below DBL_MAX, so this spans the representable domain
struct fact_table_t {
    double v[171]{};

    constexpr fact_table_t()
    {
        v[0] = 1.0;
        for (int i = 1; i < 171; ++i) {
            v[i] = v[i - 1] * static_cast<double>(i);
        }
    }
};

inline constexpr fact_table_t fact_table{};

} // namespace detail

constexpr double fact(int n)
{
    if (n < 0) {
        throw std::invalid_argument("hd::ga::fact: negative argument.");
    }
    if (n > 170) {
        return std::numeric_limits<double>::infinity(); // n! exceeds DBL_MAX
    }
    return detail::fact_table.v[n];
}

enum class stencil_lhs {
    f1, // f'  terms considered to be on lhs of the finite difference
    f2  // f'' terms considered to be on lhs of the finite difference
};

struct stencil_t {

    // after calling the ctor, the "output values" below can be used
    stencil_t(double x0, stencil_lhs lhs_t, std::vector<double> xf0,
              std::vector<double> xf1, std::vector<double> xf2) :
        x0{x0}, lhs_t{lhs_t}, xf0{std::move(xf0)}, xf1{std::move(xf1)},
        xf2{std::move(xf2)}
    {
        // consistency checks (nf0() >= 2: fewer function samples cannot pin a
        // derivative, and nf0() == 0 would underflow the j0e index below)
        if ((nf1() == 0 && nf2() == 0) || n() < 3 || nf0() < 2 ||
            (nf1() == 0 && lhs_t == stencil_lhs::f1) ||
            (nf2() == 0 && lhs_t == stencil_lhs::f2)) {
            throw std::invalid_argument(
                "Inconsistent specification of stencil in ctor of hd::ga::stencil_t.");
        }

        wf0.reserve(this->xf0.size());
        wf1.reserve(this->xf1.size());
        wf2.reserve(this->xf2.size());

        // provide the weights, order and truncation error for further use
        calc_stencil();
    }

    // input values provided by ctor

    double const x0; // development point of stencil
                     // (should be within or at least close to the point coordinates)

    stencil_lhs const lhs_t; // f1 or f2 terms on lhs; all other terms on rhs

    std::vector<double> const xf0; // coordinates of points for f
    std::vector<double> const xf1; // coordinates of points for f'
    std::vector<double> const xf2; // coordinates of points for f''

    // output values filled in by calc_stencil()

    std::vector<double> wf0; // weights of points for f
    std::vector<double> wf1; // weights of points for f'
    std::vector<double> wf2; // weights of points for f''

    int order{0};          // order of the fd stencil
    double trunc_err{0.0}; // coefficient of the leading neglected Taylor term

    // helpers for number of points
    size_t nf0() const { return xf0.size(); }          // number of points for f
    size_t nf1() const { return xf1.size(); }          // number of points for f'
    size_t nf2() const { return xf2.size(); }          // number of points for f''
    size_t n() const { return nf0() + nf1() + nf2(); } // total number of points

  private:

    void calc_stencil();
};

inline void stencil_t::calc_stencil()
{
    // scratch storage for matrix, permutation and rhs vector, initialized with 0
    std::vector<double> mem_matrix(n() * n(), 0.0);
    std::vector<int> mem_perm(n(), 0);
    std::vector<double> mem_rhs(n(), 0.0);

    // create views
    std::mdspan matrix{mem_matrix.data(), n(), n()};
    std::mdspan perm{mem_perm.data(), n()};
    std::mdspan rhs{mem_rhs.data(), n()};

    // setup column index ranges (begin/end indices for the f, f', f'' columns);
    // the f'/f'' ranges are set and used only under the nf1()/nf2() guards below,
    // but MSVC cannot prove the guards match, so zero-initialize to silence C4701
    size_t j0b, j0e, col;
    size_t j1b = 0, j1e = 0, j2b = 0, j2e = 0;

    // f
    j0b = 0;
    j0e = nf0() - 1;
    col = j0e;

    // f'
    if (nf1() > 0) {
        j1b = col + 1;
        j1e = col + nf1();
        col += nf1();
    }

    // f''
    if (nf2() > 0) {
        j2b = col + 1;
        j2e = col + nf2();
        // col += nf2(); // only needed for a further extension with f'''
    }

    // setup the Taylor-matching matrix columnwise;
    // unwanted terms of the series expansion of the lhs are moved to the rhs (sfact)

    // f
    for (size_t j = j0b; j <= j0e; ++j) {
        matrix[0, j] = 1.0;
        for (size_t i = 1; i < n(); ++i) {
            matrix[i, j] = std::pow(xf0[j - j0b] - x0, i) / fact(static_cast<int>(i));
        }
    }

    double sfact; // locates terms on lhs (-1.0) or on rhs (+1.0)

    // f'
    double sfact1 = 1.0;
    if (nf1() > 0) {
        sfact1 = (lhs_t == stencil_lhs::f1) ? -1.0 : 1.0;
        sfact = sfact1;
        for (size_t j = j1b; j <= j1e; ++j) {
            matrix[0, j] = 0.0;
            matrix[1, j] = 1.0;
            for (size_t i = 2; i < n(); ++i) {
                matrix[i, j] = sfact * std::pow(xf1[j - j1b] - x0, i - 1) /
                               fact(static_cast<int>(i) - 1);
            }
        }
    }

    // f''
    double sfact2 = 1.0;
    if (nf2() > 0) {
        sfact2 = (lhs_t == stencil_lhs::f2) ? -1.0 : 1.0;
        sfact = sfact2;
        for (size_t j = j2b; j <= j2e; ++j) {
            matrix[0, j] = 0.0;
            matrix[1, j] = 0.0;
            matrix[2, j] = 1.0;
            for (size_t i = 3; i < n(); ++i) {
                matrix[i, j] = sfact * std::pow(xf2[j - j2b] - x0, i - 2) /
                               fact(static_cast<int>(i) - 2);
            }
        }
    }

    // setup rhs
    if (lhs_t == stencil_lhs::f1) {
        rhs[1] = 1.0;
    }
    if (lhs_t == stencil_lhs::f2) {
        rhs[2] = 1.0;
    }

    // normalization: replace the last equation with the normalization condition (sum of
    // coefficients on lhs = 1), i.e. set the coefficients of the primary derivative to
    // 1.0 in the last equation and remove them from the matrix line they came from
    for (size_t j = 0; j < n(); ++j) {
        matrix[n() - 1, j] = 0.0;
    }
    rhs[n() - 1] = 1.0;

    if (lhs_t == stencil_lhs::f1) {
        for (size_t j = j1b; j <= j1e; ++j) {
            matrix[n() - 1, j] = 1.0;
            matrix[1, j] = 0.0;
        }
    }
    if (lhs_t == stencil_lhs::f2) {
        for (size_t j = j2b; j <= j2e; ++j) {
            matrix[n() - 1, j] = 1.0;
            matrix[2, j] = 0.0;
        }
    }

    // solve system
    lu_decomp(matrix, perm);
    lu_backsubs(matrix, perm, rhs); // weights are now on rhs

    // assign weights to the output vectors
    for (size_t j = j0b; j <= j0e; ++j) {
        wf0.push_back(rhs[j]);
    }
    if (nf1() > 0) {
        for (size_t j = j1b; j <= j1e; ++j) {
            wf1.push_back(rhs[j]);
        }
    }
    if (nf2() > 0) {
        for (size_t j = j2b; j <= j2e; ++j) {
            wf2.push_back(rhs[j]);
        }
    }

    // compute order and truncation error: continue the Taylor expansion beyond the
    // matched orders and report the FIRST non-vanishing residual term (lhs terms enter
    // with the same sign convention used to build the matrix above).
    //
    // The vanishing test is RELATIVE to the magnitude of the contributions the sum is
    // built from. A matched order cancels to zero only up to rounding, so its residual
    // sits ~1e-16 below that magnitude, while the leading term is a finite fraction of
    // it -- the two are separated by many orders of magnitude no matter how the node
    // set is scaled. An ABSOLUTE threshold is not: it reads every term as vanishing
    // for a small physical spacing, and equally for a wide stencil normalized into
    // [-1, 1], where the leading coefficient is itself tiny (a 9-point central stencil
    // so normalized has |trunc_err| ~ 1e-7).
    //
    // The scan runs two orders past n() because a symmetric node set annihilates whole
    // orders by parity, which can push the leading term beyond the matched range.

    double const rel_tol = 1.0e-8; // >> rounding (~1e-16), << any genuine term

    for (size_t i = nf0(); i <= n() + 2; ++i) {

        double sumte = 0.0;  // the residual of Taylor order i
        double sumabs = 0.0; // magnitude of the contributions it cancels between

        auto add = [&](double term) {
            sumte += term;
            sumabs += std::abs(term);
        };

        // f
        for (size_t j = j0b; j <= j0e; ++j) {
            add(std::pow(xf0[j - j0b] - x0, i) / fact(static_cast<int>(i)) * rhs[j]);
        }

        // f'
        if (nf1() > 0) {
            for (size_t j = j1b; j <= j1e; ++j) {
                add(sfact1 * std::pow(xf1[j - j1b] - x0, i - 1) /
                    fact(static_cast<int>(i) - 1) * rhs[j]);
            }
        }

        // f''
        if (nf2() > 0) {
            for (size_t j = j2b; j <= j2e; ++j) {
                add(sfact2 * std::pow(xf2[j - j2b] - x0, i - 2) /
                    fact(static_cast<int>(i) - 2) * rhs[j]);
            }
        }

        if (sumabs > 0.0 && std::abs(sumte) > rel_tol * sumabs) {
            trunc_err = sumte;
            if (lhs_t == stencil_lhs::f1) {
                order = static_cast<int>(i) - 1;
            }
            if (lhs_t == stencil_lhs::f2) {
                order = static_cast<int>(i) - 2;
            }
            break; // leading term found
        }
    }

    // No leading term within the scan means the stencil is not characterizable -- an
    // ill-posed node set, or one so nearly singular that the solve produced noise.
    // Reporting order 0 would look like a valid answer and silently poison anything
    // deriving a step size from it, so fail here instead.
    if (order == 0) {
        throw std::runtime_error("hd::ga::stencil_t: could not determine the order of "
                                 "the stencil (degenerate or ill-conditioned node set).");
    }
}

} // namespace hd::ga
