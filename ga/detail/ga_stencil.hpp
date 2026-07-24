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
// NOTE: the order/truncation-error detection compares the residual of successive Taylor
// orders against an absolute threshold (1e-6), which assumes O(1) node spacing. For
// very small h compute the stencil in normalized coordinates (x/h) instead.
//
// Adapted from the hd utility library (hd_stencil.hpp) and made internal to the ga
// library so the field-operator layers carry no external dependency. Two fixes relative
// to the source version: the leading truncation term is the FIRST non-vanishing
// residual order (the source kept overwriting, reporting the last), and lhs terms enter
// the residual with the same sign convention used to build the matrix.
/////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>     // std::pow, std::abs
#include <mdspan>    // std::mdspan (views onto the scratch storage)
#include <stdexcept> // std::invalid_argument
#include <vector>    // std::vector

#include "ga_solver.hpp" // lu_decomp, lu_backsubs

namespace hd::ga {

// factorial as double (stencil matrix entries are Taylor coefficients 1/i!)
constexpr double factorial(size_t n)
{
    double f = 1.0;
    for (size_t i = 2; i <= n; ++i) {
        f *= static_cast<double>(i);
    }
    return f;
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
        // consistency checks
        if ((nf1() == 0 && nf2() == 0) || n() < 3 ||
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
            matrix[i, j] = std::pow(xf0[j - j0b] - x0, i) / factorial(i);
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
                matrix[i, j] =
                    sfact * std::pow(xf1[j - j1b] - x0, i - 1) / factorial(i - 1);
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
                matrix[i, j] =
                    sfact * std::pow(xf2[j - j2b] - x0, i - 2) / factorial(i - 2);
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
    // with the same sign convention used to build the matrix above)
    for (size_t i = nf0(); i <= n(); ++i) {

        double sumte = 0.0;

        // f
        for (size_t j = j0b; j <= j0e; ++j) {
            sumte += std::pow(xf0[j - j0b] - x0, i) / factorial(i) * rhs[j];
        }

        // f'
        if (nf1() > 0) {
            for (size_t j = j1b; j <= j1e; ++j) {
                sumte += sfact1 * std::pow(xf1[j - j1b] - x0, i - 1) / factorial(i - 1) *
                         rhs[j];
            }
        }

        // f''
        if (nf2() > 0) {
            for (size_t j = j2b; j <= j2e; ++j) {
                sumte += sfact2 * std::pow(xf2[j - j2b] - x0, i - 2) / factorial(i - 2) *
                         rhs[j];
            }
        }

        double const order_tol = 1.0e-6;
        if (std::abs(sumte) > order_tol) {
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
}

} // namespace hd::ga
