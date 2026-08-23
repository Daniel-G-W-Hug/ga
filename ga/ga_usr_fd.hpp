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
// - fd_scheme                     the scheme description, in stencil_t's full model
// - make_scheme(x0, deriv, ...)   the general builder: any node set, x0 need NOT be a
//                                 node, and f' / f'' data may be supplied alongside f
// - central_scheme(deriv,order)   equidistant explicit scheme, normalised offsets
// - staggered_scheme(deriv,order) explicit scheme at the MIDPOINT between its nodes
// - compact_scheme(deriv)         equidistant O4 Pade scheme on i-1, i, i+1
// - scheme_at(x, i, ...)          scheme from arbitrary node coordinates
// - apply_scheme<T>(sc, ...)      evaluate an offset scheme against a sampled callable
// - fd_derivative(x, f, ...)      derivative at every node of a 1d grid
// - fd_step(scheme, scale)        step balancing truncation against round-off
/////////////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // std::min, std::max
#include <cmath>     // std::pow, std::abs
#include <cstddef>   // std::size_t
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::invalid_argument
#include <utility>   // std::move
#include <vector>    // std::vector

#include "detail/ga_stencil.hpp" // stencil_t (finite-difference weight generator)
#include "ga_value_t.hpp"        // value_t

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// finite-difference schemes for the field operators
////////////////////////////////////////////////////////////////////////////////
//
// A scheme is a set of sample points together with the weights stencil_t generates
// for them. It carries stencil_t's full model, which is wider than "f at a few
// points" in three ways worth knowing about:
//
//   - the development point x0 need NOT coincide with a node. Sampling on the
//     half-offsets +-1/2, +-3/2 and developing at 0 gives the midpoint (staggered)
//     derivative a Yee-type lattice needs -- see staggered_scheme().
//   - the sought derivative may live at SEVERAL points, which makes the scheme
//     implicit (compact / Pade) -- see compact_scheme() and is_compact().
//   - values of the OTHER derivative order may be supplied as rhs data (Hermite),
//     which buys accuracy from the same footprint wherever the caller already knows
//     them analytically -- see the aux_nodes argument of make_scheme().
//
// The equidistant builders lay their nodes out at UNIT spacing, so their weights are
// independent of the step h and applying the scheme scales them by 1/h^deriv. Nothing
// requires that: stencil_t solves on whatever coordinates it is given, which is what
// scheme_at() feeds it from a stretched grid. Building a scheme costs a small LU
// solve, so build it once and reuse it; the operators keep a cached default and take
// a scheme by reference otherwise.
////////////////////////////////////////////////////////////////////////////////

enum class fd_kind {
    explicit_fd, // all derivative information on one point: f'(x0) = sum w_k f(x_k)
    compact_fd   // Pade: derivative values at several points on the lhs (implicit)
};

struct fd_scheme {
    std::vector<double> nodes;       // sample points for f
    std::vector<double> weights;     // rhs weights for those samples
    std::vector<double> lhs_nodes;   // points carrying the SOUGHT derivative (lhs)
    std::vector<double> lhs_weights; // lhs weights, normalised to sum to 1
    std::vector<double> aux_nodes;   // points carrying the OTHER derivative order (rhs)
    std::vector<double> aux_weights; // weights for that supplied derivative data
    double x0{0.0};                  // development point; need NOT be one of the nodes
    int deriv{1};                    // which derivative (1 or 2)
    int order{0};                    // achieved order of accuracy
    double trunc_err{0.0};           // coefficient of the leading neglected term

    // index of nodes[0] in the grid scheme_at() was built from; meaningless (and
    // unused) for the offset-style schemes, which carry no grid
    std::size_t node0{0};

    bool offset{false}; // nodes are step multiples, not absolute grid coordinates

    // an EXPLICIT scheme carries the derivative at a single point, so it evaluates
    // directly; a COMPACT one couples neighbouring derivative values and needs a solve
    bool is_compact() const { return lhs_nodes.size() > 1; }

    // a HERMITE scheme additionally consumes supplied values of the other derivative
    // order, so a plain sampler of f is not enough to evaluate it
    bool is_hermite() const { return !aux_nodes.empty(); }

    // the derivative order the aux_nodes carry (the one that is not sought)
    int aux_deriv() const { return (deriv == 1) ? 2 : 1; }

    // an OFFSET scheme has its nodes in units of a step still to be chosen, so it can
    // be applied anywhere; a scheme from scheme_at() carries absolute coordinates and
    // is tied to the grid it was built from. Recorded rather than inferred from
    // x0 == 0.0, which a grid with a node at the origin would satisfy by accident.
    bool is_offset() const { return offset; }
};

// The general builder -- stencil_t's full model, with the node coordinates in the
// caller's own units.
//
//   x0         development point; need NOT be one of the nodes
//   deriv      the derivative sought (1 or 2); selects which set is the implicit lhs
//   nodes      points where f is sampled
//   lhs_nodes  points carrying the sought derivative (one -> explicit, several ->
//              compact/implicit)
//   aux_nodes  points where values of the OTHER derivative order are supplied as rhs
//              data (Hermite); empty for an ordinary scheme
//
// The coordinates are handed to stencil_t as they are. An earlier version normalised
// them to unit scale and scaled the weights back afterwards, to work around an
// ABSOLUTE threshold in stencil_t's order detection; that threshold is relative now,
// and measurement across fourteen decades of spacing shows the raw solve loses
// nothing (the Taylor matrix scales column-wise, which the pivoting absorbs). So the
// normalisation bought only the bookkeeping needed to undo it -- one rescaling per
// output quantity, each a place to put the exponent wrong.
inline fd_scheme make_scheme(double x0, int deriv, std::vector<double> nodes,
                             std::vector<double> lhs_nodes,
                             std::vector<double> aux_nodes = {})
{
    if (deriv != 1 && deriv != 2) {
        throw std::invalid_argument("fd scheme: deriv must be 1 or 2");
    }
    if (lhs_nodes.empty()) {
        throw std::invalid_argument("fd scheme: need at least one lhs node");
    }

    // route the two derivative point sets: the sought order goes on the lhs, the
    // supplied one (if any) stays on the rhs as data
    stencil_t const s(x0, (deriv == 1) ? stencil_lhs::f1 : stencil_lhs::f2, nodes,
                      (deriv == 1) ? lhs_nodes : aux_nodes,
                      (deriv == 1) ? aux_nodes : lhs_nodes);

    fd_scheme out;
    out.nodes = std::move(nodes);
    out.weights = s.wf0;
    out.lhs_nodes = std::move(lhs_nodes);
    out.lhs_weights = (deriv == 1) ? s.wf1 : s.wf2;
    out.aux_nodes = std::move(aux_nodes);
    out.aux_weights = (deriv == 1) ? s.wf2 : s.wf1;
    out.x0 = x0;
    out.deriv = deriv;
    out.order = s.order;
    out.trunc_err = s.trunc_err;
    return out;
}

// central EXPLICIT scheme of the requested even accuracy order, equidistant/normalised.
// A centred node set annihilates one Taylor order by parity, so order+1 nodes reach
// accuracy `order` for the first AND the second derivative alike.
inline fd_scheme central_scheme(int deriv, int order = 2)
{
    if (order < 2 || order % 2 != 0) {
        throw std::invalid_argument("central_scheme: order must be even and >= 2");
    }
    int const m = order / 2;
    std::vector<double> nodes;
    for (int i = -m; i <= m; ++i)
        nodes.push_back(double(i));
    fd_scheme sc = make_scheme(0.0, deriv, std::move(nodes), {0.0});
    sc.offset = true;
    return sc;
}

// STAGGERED explicit scheme: the derivative at the MIDPOINT between its samples, taken
// on the half-offsets +-1/2, +-3/2, ... This is the arrangement a Yee-type staggered
// lattice asks for, where a field component lives at a cell face and its derivative is
// wanted at the cell centre (or the other way round). The development point is 0, so
// none of the nodes coincides with it.
inline fd_scheme staggered_scheme(int deriv, int order = 2)
{
    if (order < 2 || order % 2 != 0) {
        throw std::invalid_argument("staggered_scheme: order must be even and >= 2");
    }
    // a first derivative reaches `order` from `order` half-offset nodes; a second
    // derivative needs one node pair more, having no sample at the centre to lean on
    int const m = order / 2 + deriv - 1;
    std::vector<double> nodes;
    for (int i = -m; i <= m; ++i) {
        if (i == 0) continue;
        nodes.push_back((i < 0) ? (double(i) + 0.5) : (double(i) - 0.5));
    }
    fd_scheme sc = make_scheme(0.0, deriv, std::move(nodes), {0.0});
    sc.offset = true;
    return sc;
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
    fd_scheme sc = make_scheme(0.0, deriv, {-1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0});
    sc.offset = true;
    return sc;
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
        fd_scheme sc = make_scheme(x[i], deriv, nb, nb);
        sc.node0 = i - 1;
        return sc;
    }

    // Explicit. A CENTRED stencil gains one order from parity, so order+1 nodes reach
    // accuracy `order` for either derivative. A one-sided closure at a grid end has no
    // such symmetry and needs order+deriv nodes for the same accuracy -- widen only
    // there, so an interior node keeps the classical minimal centred footprint (three
    // points at order 2, not four).
    std::size_t const n_centred = std::size_t(order) + 1;
    std::size_t const n_sided = std::size_t(order) + std::size_t(deriv);
    bool const centred = (i >= n_centred / 2) && (i + n_centred / 2 < x.size());

    std::size_t const n = std::min(x.size(), centred ? n_centred : n_sided);
    std::size_t begin = (i < n / 2) ? 0 : i - n / 2;
    if (begin + n > x.size()) begin = x.size() - n;

    std::vector<double> nodes(x.begin() + long(begin), x.begin() + long(begin + n));
    fd_scheme sc = make_scheme(x[i], deriv, std::move(nodes), {x[i]});
    sc.node0 = begin;
    return sc;
}

////////////////////////////////////////////////////////////////////////////////
// applying a scheme
////////////////////////////////////////////////////////////////////////////////

// Evaluate an OFFSET scheme (x0 == 0, nodes in units of the step) against a callable
// that samples the field at a given offset from the evaluation point. The accumulator
// type T is explicit because the sampled type is usually narrower than the type the
// result is accumulated in -- a vector field summed into a multivector, say.
//
// This is the one place the weights are applied; the per-algebra field operators all
// route through it so the loop, its zero-weight skip and its 1/h^deriv scaling exist
// once rather than once per algebra and operator.
template <typename T, typename F>
inline T apply_scheme(fd_scheme const& sc, F&& sample, double h)
{
    if (sc.is_compact()) {
        throw std::invalid_argument(
            "apply_scheme: compact schemes are implicit and need a grid; "
            "use fd_derivative() or pass an explicit scheme");
    }
    if (sc.is_hermite()) {
        throw std::invalid_argument(
            "apply_scheme: a Hermite scheme also consumes supplied derivative values, "
            "which a sampler of f alone cannot provide");
    }
    if (!sc.is_offset()) {
        throw std::invalid_argument(
            "apply_scheme: needs an offset scheme (x0 == 0); a scheme from scheme_at() "
            "is tied to the grid it was built from");
    }
    if (h <= 0.0) throw std::invalid_argument("apply_scheme: step must be positive");

    T acc{};
    for (std::size_t k = 0; k < sc.nodes.size(); ++k) {
        if (sc.weights[k] == 0.0) continue; // the centre drops out of odd stencils
        acc = acc + sc.weights[k] * T(sample(sc.nodes[k] * h));
    }
    return acc * ((sc.deriv == 1) ? (1.0 / h) : (1.0 / (h * h)));
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
// the derivative at i-1, i and i+1, so the interior rows form a tridiagonal system,
// handed to tridiag_solve() (detail/ga_solver.hpp) for one O(N) pass. The two
// boundary rows are closed with one-sided explicit stencils of the requested order,
// which is what keeps the matrix tridiagonal -- a wider closure would spill outside
// the three diagonals and force the general solver.
//
// Note that a compact scheme reaches 4th order from a three-point footprint where an
// explicit stencil of the same width reaches only 2nd -- that is the point of it. On
// a NON-uniform grid the first-derivative scheme keeps 4th order, while the
// second-derivative one drops to 3rd; scheme_at() reports the order it achieved, so
// check it rather than assuming.
////////////////////////////////////////////////////////////////////////////////

namespace detail {

// gather one grid scheme against the samples it was built from. scheme_at() records
// where its node window starts, so the samples are indexed directly -- no search back
// from the node coordinates, whose fuzzy comparison would be scale-dependent.
template <typename T>
inline T gather_grid_scheme(fd_scheme const& sc, std::vector<T> const& f)
{
    T acc{};
    for (std::size_t k = 0; k < sc.nodes.size(); ++k) {
        if (sc.weights[k] == 0.0) continue;
        acc = acc + sc.weights[k] * f[sc.node0 + k];
    }
    return acc;
}

} // namespace detail

template <typename T>
inline std::vector<T> fd_derivative(std::vector<double> const& x, std::vector<T> const& f,
                                    int deriv, fd_kind kind = fd_kind::explicit_fd,
                                    int order = 2)
{
    std::size_t const n = x.size();
    if (n != f.size()) throw std::invalid_argument("fd_derivative: size mismatch");
    if (n < 3) throw std::invalid_argument("fd_derivative: need at least 3 nodes");

    if (kind == fd_kind::explicit_fd) {
        std::vector<T> d(n, T{});
        for (std::size_t i = 0; i < n; ++i) {
            d[i] = detail::gather_grid_scheme(scheme_at(x, i, deriv, kind, order), f);
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
        rhs[i] = detail::gather_grid_scheme(sc, f);
    }
    // boundaries: one-sided EXPLICIT closure, so the matrix stays tridiagonal
    for (std::size_t i : {std::size_t(0), n - 1}) {
        // a compact solve couples every node, so an under-resolved boundary pollutes
        // the whole interior: close at 4th order to match the interior scheme
        fd_scheme const sc =
            scheme_at(x, i, deriv, fd_kind::explicit_fd, std::max(order, 4));
        a[i] = 0.0;
        b[i] = 1.0;
        c[i] = 0.0;
        rhs[i] = detail::gather_grid_scheme(sc, f);
    }

    // one O(N) tridiagonal solve -- see tridiag_solve() in detail/ga_solver.hpp, which
    // is why a compact scheme costs no more than an explicit one per node
    return tridiag_solve(a, b, c, rhs);
}

// Step size balancing truncation against round-off. For an order-p scheme of the d-th
// derivative the two error sources are
//
//     truncation  ~ |trunc_err| * h^p        round-off ~ eps * sum|w_k| / h^d
//
// and their sum is smallest where the derivatives balance, at
//
//     h = ( d * sum|w_k| * eps / (p * |trunc_err|) )^(1/(p+d))
//
// Both stencil-dependent factors are known: the leading truncation COEFFICIENT is what
// stencil_t computes alongside the order, and sum|w_k| is the factor by which a
// relative sampling error is amplified into the estimate. Using them instead of the
// exponent alone reproduces the textbook optimum -- for the central first difference
// it gives cbrt(3 eps), where the exponent-only form is off by a factor ~1.4.
//
// The result is scaled by the magnitude of the evaluation point so the choice tracks
// the field's own scale. Meant for OFFSET schemes; a grid scheme from scheme_at()
// already carries its spacing in its weights.
inline double fd_step(fd_scheme const& sc, double scale)
{
    double const e = std::numeric_limits<double>::epsilon();
    double const p = double(sc.order);
    double const d = double(sc.deriv);

    double w = 0.0;
    for (double wi : sc.weights)
        w += std::abs(wi);

    // fall back to the exponent-only estimate if the scheme carries no leading
    // coefficient (it cannot, as built here -- stencil_t throws rather than report
    // order 0 -- but a hand-filled fd_scheme might)
    double const h0 =
        (p > 0.0 && w > 0.0 && sc.trunc_err != 0.0)
            ? std::pow(d * w * e / (p * std::abs(sc.trunc_err)), 1.0 / (p + d))
            : std::pow(e, 1.0 / (std::max(p, 1.0) + d));

    return h0 * std::max(1.0, std::abs(scale));
}

} // namespace hd::ga
