#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// Finite-difference schemes and grid derivatives.
//
// Algebra-INDEPENDENT numerics: nothing here knows about a geometric algebra. The
// sample type is generic, so scalars and multivectors are differentiated alike, and
// the per-algebra field-operator layers (ga_ega{2,3}d_ops_calculus.hpp) build on it.
//
// - fd_kind, fd_scheme          the scheme description (explicit or compact/Pade)
// - central_scheme(deriv,order) equidistant explicit scheme, normalised offsets
// - compact_scheme(deriv)       equidistant O4 Pade scheme on i-1, i, i+1
// - scheme_at(x, i, ...)        scheme from arbitrary node coordinates
// - fd_derivative(x, f, ...)    derivative at every node of a 1d grid
// - fd_step(scheme, scale)      step balancing truncation against round-off
/////////////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // std::lower_bound, std::min, std::max
#include <cmath>     // std::pow, std::abs
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::invalid_argument
#include <vector>    // std::vector

#include "detail/ga_stencil.hpp" // stencil_t (finite-difference weight generator)
#include "ga_value_t.hpp"        // value_t

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// finite-difference schemes for the field operators
////////////////////////////////////////////////////////////////////////////////
//
// A scheme is a set of sample offsets in units of the step h together with the
// weights stencil_t generates for them. The stencil is built in NORMALIZED
// coordinates (spacing exactly 1), which makes the weights independent of h --
// applying the scheme scales them by 1/h^deriv -- and side-steps the absolute
// order-detection threshold documented in ga_stencil.hpp, which assumes O(1)
// node spacing and would misreport the order for a small physical h.
//
// Building a scheme costs a small LU solve, so build it once and reuse it; the
// operators keep a cached default and take a scheme by reference otherwise.
////////////////////////////////////////////////////////////////////////////////

enum class fd_kind {
    explicit_fd, // all derivative information on one point: f'(x0) = sum w_k f(x_k)
    compact_fd   // Pade: derivative values at several points on the lhs (implicit)
};

struct fd_scheme {
    std::vector<double> nodes;       // sample points for f
    std::vector<double> weights;     // rhs weights for those samples
    std::vector<double> lhs_nodes;   // points carrying the derivative (lhs)
    std::vector<double> lhs_weights; // lhs weights, normalised to sum to 1
    int deriv{1};                    // which derivative (1 or 2)
    int order{0};                    // achieved order of accuracy
    double trunc_err{0.0};           // coefficient of the leading neglected term

    // an EXPLICIT scheme carries the derivative at a single point, so it evaluates
    // directly; a COMPACT one couples neighbouring derivative values and needs a solve
    bool is_compact() const { return lhs_nodes.size() > 1; }
};

namespace detail {

// build a scheme from arbitrary node coordinates -- the general case behind both the
// equidistant helpers and the non-uniform grid path
inline fd_scheme make_scheme(double x0, int deriv, std::vector<double> nodes,
                             std::vector<double> lhs_nodes)
{
    if (deriv != 1 && deriv != 2) {
        throw std::invalid_argument("fd scheme: deriv must be 1 or 2");
    }

    // Build the stencil in NORMALIZED coordinates, (x - x0)/h_loc. ga_stencil.hpp
    // detects the order by comparing Taylor residuals against an ABSOLUTE threshold,
    // which assumes O(1) node spacing: on a grid with a small physical spacing it
    // would both misreport the order and solve a badly conditioned system. Working in
    // units of the local scale removes both problems; the weights are scaled back into
    // the caller's units on return.
    double h_loc = 0.0;
    for (double xn : nodes)
        h_loc = std::max(h_loc, std::abs(xn - x0));
    if (h_loc <= 0.0) throw std::invalid_argument("fd scheme: degenerate node set");

    auto norm = [&](std::vector<double> const& v) {
        std::vector<double> out;
        out.reserve(v.size());
        for (double xn : v)
            out.push_back((xn - x0) / h_loc);
        return out;
    };

    stencil_t const s(0.0, (deriv == 1) ? stencil_lhs::f1 : stencil_lhs::f2, norm(nodes),
                      (deriv == 1) ? norm(lhs_nodes) : std::vector<double>{},
                      (deriv == 2) ? norm(lhs_nodes) : std::vector<double>{});

    // rhs weights carry the derivative's units; the lhs weights are dimensionless
    // (normalised to sum to 1) and need no rescaling
    double const inv = 1.0 / std::pow(h_loc, double(deriv));
    std::vector<double> w = s.wf0;
    for (double& wi : w)
        wi *= inv;

    return {std::move(nodes),
            std::move(w),
            std::move(lhs_nodes),
            (deriv == 1) ? s.wf1 : s.wf2,
            deriv,
            s.order,
            s.trunc_err * std::pow(h_loc, double(s.order))};
}

} // namespace detail

// central EXPLICIT scheme of the requested even accuracy order, equidistant/normalised
inline fd_scheme central_scheme(int deriv, int order = 2)
{
    if (order < 2 || order % 2 != 0) {
        throw std::invalid_argument("central_scheme: order must be even and >= 2");
    }
    int const m = order / 2;
    std::vector<double> nodes;
    for (int i = -m; i <= m; ++i)
        nodes.push_back(double(i));
    return detail::make_scheme(0.0, deriv, std::move(nodes), {0.0});
}

// COMPACT (Pade) scheme on the three points i-1, i, i+1, equidistant/normalised.
// It carries f AND the derivative at all three points, which buys 4th order from a
// three-point footprint where an explicit stencil of the same width reaches only 2nd:
//
//   deriv 1:  f'_{i-1} + 4 f'_i + f'_{i+1} = 3 (f_{i+1} - f_{i-1}) / h
//   deriv 2:  f''_{i-1} + 10 f''_i + f''_{i+1} = 12 (f_{i-1} - 2 f_i + f_{i+1}) / h^2
//
// (the weights below are those relations normalised so the lhs sums to 1). Being
// implicit, a compact scheme cannot be evaluated at a single point -- it defines a
// banded system over the whole grid, which fd_derivative() solves.
inline fd_scheme compact_scheme(int deriv)
{
    return detail::make_scheme(0.0, deriv, {-1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0});
}

// scheme at node i of an ARBITRARY (possibly non-uniform) grid, from the coordinates
// themselves. The weights then differ from node to node, which is exactly what a
// stretched grid requires.
inline fd_scheme scheme_at(std::vector<double> const& x, std::size_t i, int deriv,
                           fd_kind kind = fd_kind::explicit_fd, int order = 2)
{
    if (x.size() < 3) throw std::invalid_argument("scheme_at: need at least 3 nodes");

    if (kind == fd_kind::compact_fd) {
        if (i == 0 || i + 1 >= x.size()) {
            throw std::invalid_argument("scheme_at: compact needs an interior node");
        }
        std::vector<double> const nb{x[i - 1], x[i], x[i + 1]};
        return detail::make_scheme(x[i], deriv, nb, nb);
    }

    // explicit: take the order+1 nodes centred on i as far as the ends allow, so the
    // same call also produces the one-sided closures at the boundaries
    std::size_t const want = std::size_t(order) + std::size_t(deriv == 2 ? 1 : 1);
    std::size_t const n = std::min(x.size(), want + 1);
    std::size_t begin = (i < n / 2) ? 0 : i - n / 2;
    if (begin + n > x.size()) begin = x.size() - n;
    std::vector<double> nodes(x.begin() + long(begin), x.begin() + long(begin + n));
    return detail::make_scheme(x[i], deriv, std::move(nodes), {x[i]});
}

////////////////////////////////////////////////////////////////////////////////
// derivatives of sampled data on a 1d grid (equidistant or not)
////////////////////////////////////////////////////////////////////////////////
//
// Given node coordinates x and samples f (any type supporting addition and scaling
// by a double -- a scalar or any multivector), return the derivative at EVERY node.
//
// The grid may be uniform or arbitrarily stretched: the weights are generated from
// the coordinates themselves, so a non-uniform spacing costs nothing but a stencil
// per node.
//
// Explicit schemes evaluate node by node. A COMPACT scheme cannot: its lhs couples
// the derivative at i-1, i and i+1, so the interior rows form a tridiagonal system
// which is solved here in one pass (Thomas algorithm, O(N)). The two boundary rows
// are closed with one-sided explicit stencils of the requested order, which keeps
// the matrix tridiagonal.
//
// Note that a compact scheme reaches 4th order from a three-point footprint where an
// explicit stencil of the same width reaches only 2nd -- that is the point of it. On
// a NON-uniform grid the first-derivative scheme keeps 4th order, while the
// second-derivative one drops to 3rd; scheme_at() reports the order it achieved, so
// check it rather than assuming.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline std::vector<T> fd_derivative(std::vector<double> const& x, std::vector<T> const& f,
                                    int deriv, fd_kind kind = fd_kind::explicit_fd,
                                    int order = 2)
{
    std::size_t const n = x.size();
    if (n != f.size()) throw std::invalid_argument("fd_derivative: size mismatch");
    if (n < 3) throw std::invalid_argument("fd_derivative: need at least 3 nodes");

    std::vector<T> d(n, T{});

    if (kind == fd_kind::explicit_fd) {
        for (std::size_t i = 0; i < n; ++i) {
            fd_scheme const sc = scheme_at(x, i, deriv, kind, order);
            T acc{};
            for (std::size_t k = 0; k < sc.nodes.size(); ++k) {
                if (sc.weights[k] == 0.0) continue;
                std::size_t const j = std::size_t(
                    std::lower_bound(x.begin(), x.end(), sc.nodes[k] - 1.0e-12) -
                    x.begin());
                acc = acc + sc.weights[k] * f[j];
            }
            d[i] = acc;
        }
        return d;
    }

    // compact: tridiagonal system  a_i d_{i-1} + b_i d_i + c_i d_{i+1} = rhs_i
    std::vector<double> a(n, 0.0), b(n, 1.0), c(n, 0.0);
    std::vector<T> rhs(n, T{});

    for (std::size_t i = 1; i + 1 < n; ++i) {
        fd_scheme const sc = scheme_at(x, i, deriv, fd_kind::compact_fd);
        a[i] = sc.lhs_weights[0];
        b[i] = sc.lhs_weights[1];
        c[i] = sc.lhs_weights[2];
        T acc{};
        for (std::size_t k = 0; k < 3; ++k) {
            if (sc.weights[k] == 0.0) continue;
            acc = acc + sc.weights[k] * f[i - 1 + k];
        }
        rhs[i] = acc;
    }
    // boundaries: one-sided EXPLICIT closure, so the matrix stays tridiagonal
    for (std::size_t i : {std::size_t(0), n - 1}) {
        // a compact solve couples every node, so an under-resolved boundary pollutes
        // the whole interior: close at 4th order to match the interior scheme
        fd_scheme const sc =
            scheme_at(x, i, deriv, fd_kind::explicit_fd, std::max(order, 4));
        T acc{};
        for (std::size_t k = 0; k < sc.nodes.size(); ++k) {
            if (sc.weights[k] == 0.0) continue;
            std::size_t const j = std::size_t(
                std::lower_bound(x.begin(), x.end(), sc.nodes[k] - 1.0e-12) - x.begin());
            acc = acc + sc.weights[k] * f[j];
        }
        a[i] = 0.0;
        b[i] = 1.0;
        c[i] = 0.0;
        rhs[i] = acc;
    }

    // Thomas algorithm (tridiagonal, O(N))
    std::vector<double> cp(n, 0.0);
    std::vector<T> dp(n, T{});
    cp[0] = c[0] / b[0];
    dp[0] = (1.0 / b[0]) * rhs[0];
    for (std::size_t i = 1; i < n; ++i) {
        double const m = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / m;
        dp[i] = (1.0 / m) * (rhs[i] + (-a[i]) * dp[i - 1]);
    }
    d[n - 1] = dp[n - 1];
    for (std::size_t i = n - 1; i-- > 0;) {
        d[i] = dp[i] + (-cp[i]) * d[i + 1];
    }
    return d;
}

// Step size balancing truncation against round-off: an order-p scheme for the d-th
// derivative has total error ~ h^p + eps/h^d, minimised at h ~ eps^(1/(p+d)). Scaled
// by the magnitude of the evaluation point so the choice tracks the field's own scale.
inline double fd_step(fd_scheme const& sc, double scale)
{
    double const e = std::numeric_limits<double>::epsilon();
    double const h0 = std::pow(e, 1.0 / double(sc.order + sc.deriv));
    return h0 * std::max(1.0, std::abs(scale));
}

} // namespace hd::ga
