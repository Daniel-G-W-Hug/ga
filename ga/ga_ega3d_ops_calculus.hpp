#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_stencil.hpp" // stencil_t (finite-difference weight generator)
#include "ga_ega3d_ops.hpp"      // ega3d ops (products, basics, exp/log, rotations)

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

struct fd_scheme {
    std::vector<double> offsets; // sample points in units of h
    std::vector<double> weights; // stencil weights for those points
    int deriv{1};                // which derivative the weights approximate (1 or 2)
    int order{0};                // achieved order of accuracy
    double trunc_err{0.0};       // coefficient of the leading neglected Taylor term
};

// central scheme for the 1st or 2nd derivative at the requested even accuracy order
inline fd_scheme central_scheme(int deriv, int order = 2)
{
    if (deriv != 1 && deriv != 2) {
        throw std::invalid_argument("central_scheme: deriv must be 1 or 2");
    }
    if (order < 2 || order % 2 != 0) {
        throw std::invalid_argument("central_scheme: order must be even and >= 2");
    }

    int const m = order / 2;
    std::vector<double> offsets;
    for (int i = -m; i <= m; ++i)
        offsets.push_back(double(i));

    stencil_t const s(0.0, (deriv == 1) ? stencil_lhs::f1 : stencil_lhs::f2, offsets,
                      (deriv == 1) ? std::vector<double>{0.0} : std::vector<double>{},
                      (deriv == 2) ? std::vector<double>{0.0} : std::vector<double>{});

    return {offsets, s.wf0, deriv, s.order, s.trunc_err};
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
    Vec3d<value_t> const e = ega3d_basis<value_t>(i);
    MVec3d<value_t> acc{};
    for (std::size_t k = 0; k < sc.offsets.size(); ++k) {
        if (sc.weights[k] == 0.0) continue; // the centre point drops out of odd stencils
        acc = acc + sc.weights[k] * MVec3d<value_t>(f(r + (sc.offsets[k] * h) * e));
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
    for (std::size_t k = 0; k < sc.offsets.size(); ++k) {
        if (sc.weights[k] == 0.0) continue;
        acc = acc + sc.weights[k] * MVec3d<value_t>(f(r, t + sc.offsets[k] * h));
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
