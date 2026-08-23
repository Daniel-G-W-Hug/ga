#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_stencil.hpp" // stencil_t (finite-difference weight generator)
#include "ga_ega3d_ops.hpp"      // ega3d ops (products, basics, exp/log, rotations)

#include <algorithm>   // std::lower_bound, std::min, std::max
#include <cmath>       // std::sqrt, std::cbrt
#include <limits>      // std::numeric_limits
#include <type_traits> // std::invoke_result_t
#include <utility>     // std::forward


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


namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega3d differential operations on multivector FIELDS:
//
// - nabla()                        -> the vector derivative, sum_i e_i * d_i f
// - nabla_dot()                    -> inner derivative  (grade-lowering part of nabla)
// - nabla_wdg()                    -> outer derivative  (grade-raising part of nabla)
// - laplacian()                    -> the scalar (grade-preserving) operator sum_i d_i^2
// - d_dt()                         -> the time derivative of a time-dependent field
//
// A FIELD here is any callable, f(r) or f(r, t), returning a multivector type. Nothing
// is stored: the operators sample the callable, so analytic and tabulated fields are
// treated alike and no grid, no ownership and no allocation are involved.
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// the vector derivative (nabla operator)
////////////////////////////////////////////////////////////////////////////////
//
// nabla acts on a field by the geometric product with the basis vectors:
//
//     nabla f = sum_i e_i * d_i f
//
// which is a SINGLE object carrying both classical derivatives at once. For a field
// of definite grade k they are its two grade parts,
//
//     inner derivative   nabla . f = gr_{k-1}(nabla f)     (grade lowering)
//     outer derivative   nabla ^ f = gr_{k+1}(nabla f)     (grade raising)
//
// so gradient, divergence and curl are not three operators but three projections of
// one. In 3d, for the four field grades:
//
//     scalar field  phi : nabla phi          = grad phi                (a vector)
//     vector field  A   : gr0(nabla A)       = div A                   (a scalar)
//                         gr2(nabla A)       = I (curl A)              (a bivector)
//     bivector field B  : gr1(nabla B)       = -curl of its dual       (a vector)
//                         gr3(nabla B)       = I div of its dual       (a trivector)
//     trivector field T : gr2(nabla T)       = I grad of its coeff     (a bivector)
//
// The identities nabla ^ (nabla ^ M) = 0 and nabla . (nabla . M) = 0 hold, and contain
// curl(grad) = 0 and div(curl) = 0 as their classical cases. The Laplacian factors as
//
//     nabla . (nabla ^ X) + nabla ^ (nabla . X) = laplacian(X)
//
// so the vector derivative is a square root of the Laplacian -- a property no operator
// of classical vector calculus has.
//
// NUMERICS: the default step is a central difference, whose optimal step for a first
// derivative scales as cbrt(eps) and for a second as sqrt(sqrt(eps)); both are computed
// from the field's own scale rather than fixed, so the operators behave over a wide
// range of magnitudes. Pass an explicit h to override.
////////////////////////////////////////////////////////////////////////////////

namespace detail {

// i-th cartesian basis vector of ega3d
template <typename T> inline Vec3d<T> ega3d_basis(int i)
{
    return (i == 0) ? Vec3d<T>(T(1.0), T(0.0), T(0.0))
                    : ((i == 1) ? Vec3d<T>(T(0.0), T(1.0), T(0.0))
                                : Vec3d<T>(T(0.0), T(0.0), T(1.0)));
}

// apply a scheme along axis i to get the directional derivative of a field at r
template <typename F>
inline MVec3d<value_t> ega3d_axis_deriv(F&& f, Vec3d<value_t> const& r, int i,
                                        fd_scheme const& sc, value_t h)
{
    // a compact scheme is implicit: it couples derivative values at neighbouring nodes
    // and cannot be evaluated at a single point. Use fd_derivative() on a grid instead.
    if (sc.is_compact()) {
        throw std::invalid_argument(
            "nabla/laplacian: compact schemes are implicit and need a grid; "
            "use fd_derivative() or pass an explicit scheme");
    }
    Vec3d<value_t> const e = ega3d_basis<value_t>(i);
    MVec3d<value_t> acc{};
    for (std::size_t k = 0; k < sc.nodes.size(); ++k) {
        if (sc.weights[k] == 0.0) continue; // the centre point drops out of odd stencils
        acc = acc + sc.weights[k] * MVec3d<value_t>(f(r + (sc.nodes[k] * h) * e));
    }
    value_t const scale = (sc.deriv == 1) ? (value_t(1.0) / h) : (value_t(1.0) / (h * h));
    return acc * scale;
}

// the cached default schemes (2nd-order central), built once
inline fd_scheme const& ega3d_default_scheme(int deriv)
{
    static fd_scheme const d1 = central_scheme(1, 2);
    static fd_scheme const d2 = central_scheme(2, 2);
    return (deriv == 1) ? d1 : d2;
}

} // namespace detail

// nabla f at r -- the full vector derivative, as a multivector
template <typename F>
inline auto nabla(F&& f, Vec3d<value_t> const& r, fd_scheme const& sc,
                  value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, nrm(r));

    MVec3d<value_t> acc{};
    for (int i = 0; i < 3; ++i) {
        acc = acc + MVec3d<value_t>(detail::ega3d_basis<value_t>(i)) *
                        detail::ega3d_axis_deriv(f, r, i, sc, h);
    }
    return acc;
}

template <typename F>
inline auto nabla(F&& f, Vec3d<value_t> const& r, value_t h = value_t(0.0))
{
    return nabla(std::forward<F>(f), r, detail::ega3d_default_scheme(1), h);
}

// nabla f at (r, t) for a time-dependent field
template <typename F>
inline auto nabla(F&& f, Vec3d<value_t> const& r, value_t t, value_t h)
{
    return nabla([&](Vec3d<value_t> const& q) { return f(q, t); }, r, h);
}

// The inner and outer derivatives are the two grade projections of nabla f. Which
// projection that is depends on the grade of the FIELD, so it is deduced from the
// callable's return type: a grade-k field gives gr_{k-1} and gr_{k+1}. The extreme
// cases are identically zero (a scalar field has no inner derivative, a trivector
// field no outer one) and are reported as a zero of the natural type.

// inner derivative: nabla . f  (grade lowering)
template <typename F>
inline auto nabla_dot(F&& f, Vec3d<value_t> const& r, value_t h = value_t(0.0))
{
    using ret_t = std::invoke_result_t<F, Vec3d<value_t> const&>;
    auto const N = nabla(std::forward<F>(f), r, h);

    if constexpr (std::is_same_v<ret_t, Scalar3d<value_t>>) {
        return Scalar3d<value_t>(0.0); // a scalar field has no inner derivative
    }
    else if constexpr (std::is_same_v<ret_t, Vec3d<value_t>>) {
        return gr0(N); // div
    }
    else if constexpr (std::is_same_v<ret_t, BiVec3d<value_t>>) {
        return gr1(N);
    }
    else if constexpr (std::is_same_v<ret_t, PScalar3d<value_t>>) {
        return gr2(N);
    }
    else {
        return N; // mixed-grade field: the caller projects
    }
}

// outer derivative: nabla ^ f  (grade raising)
template <typename F>
inline auto nabla_wdg(F&& f, Vec3d<value_t> const& r, value_t h = value_t(0.0))
{
    using ret_t = std::invoke_result_t<F, Vec3d<value_t> const&>;
    auto const N = nabla(std::forward<F>(f), r, h);

    if constexpr (std::is_same_v<ret_t, Scalar3d<value_t>>) {
        return gr1(N); // grad
    }
    else if constexpr (std::is_same_v<ret_t, Vec3d<value_t>>) {
        return gr2(N); // I * curl
    }
    else if constexpr (std::is_same_v<ret_t, BiVec3d<value_t>>) {
        return gr3(N);
    }
    else if constexpr (std::is_same_v<ret_t, PScalar3d<value_t>>) {
        return PScalar3d<value_t>(0.0); // no grade 4 in 3d
    }
    else {
        return N; // mixed-grade field: the caller projects
    }
}

// the Laplacian -- a SCALAR operator, so it preserves the grade of the field
template <typename F>
inline auto laplacian(F&& f, Vec3d<value_t> const& r, fd_scheme const& sc,
                      value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, nrm(r));

    MVec3d<value_t> acc{};
    for (int i = 0; i < 3; ++i) {
        acc = acc + detail::ega3d_axis_deriv(f, r, i, sc, h);
    }
    return acc;
}

template <typename F>
inline auto laplacian(F&& f, Vec3d<value_t> const& r, value_t h = value_t(0.0))
{
    return laplacian(std::forward<F>(f), r, detail::ega3d_default_scheme(2), h);
}

// d f / dt of a time-dependent field
template <typename F>
inline auto d_dt(F&& f, Vec3d<value_t> const& r, value_t t, fd_scheme const& sc,
                 value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, t);

    MVec3d<value_t> acc{};
    for (std::size_t k = 0; k < sc.nodes.size(); ++k) {
        if (sc.weights[k] == 0.0) continue;
        acc = acc + sc.weights[k] * MVec3d<value_t>(f(r, t + sc.nodes[k] * h));
    }
    value_t const scale = (sc.deriv == 1) ? (value_t(1.0) / h) : (value_t(1.0) / (h * h));
    return acc * scale;
}

template <typename F>
inline auto d_dt(F&& f, Vec3d<value_t> const& r, value_t t, value_t h = value_t(0.0))
{
    return d_dt(std::forward<F>(f), r, t, detail::ega3d_default_scheme(1), h);
}

} // namespace hd::ga::ega
