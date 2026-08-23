#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_ega2d_ops.hpp" // ega2d ops (products, basics, exp/log, rotations)
#include "ga_usr_fd.hpp"    // fd_scheme, central/compact schemes, fd_derivative

#include <cmath>       // std::sqrt
#include <type_traits> // std::invoke_result_t, std::is_same_v
#include <utility>     // std::forward


namespace hd::ga::detail {

// Helpers for the ega2d field operators. Kept in hd::ga::detail (NOT
// hd::ga::ega::detail): opening a detail namespace inside hd::ga::ega shadows
// the hd::ga::detail names that the algebra headers use unqualified, and whether
// that breaks the build then depends on include order.

// i-th cartesian basis vector of ega2d
template <typename T> inline Vec2d<T> ega2d_basis(int i)
{
    return (i == 0) ? Vec2d<T>(T(1.0), T(0.0)) : Vec2d<T>(T(0.0), T(1.0));
}

// directional derivative of a field along axis i -- apply_scheme() owns the weight
// loop and rejects a scheme that cannot be evaluated at a single point
template <typename F>
inline MVec2d<value_t> ega2d_axis_deriv(F&& f, Vec2d<value_t> const& r, int i,
                                        fd_scheme const& sc, value_t h)
{
    Vec2d<value_t> const e = ega2d_basis<value_t>(i);
    return apply_scheme<MVec2d<value_t>>(
        sc, [&](value_t off) { return f(r + off * e); }, h);
}


} // namespace hd::ga::detail


namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega2d differential operations on multivector FIELDS:
//
// - nabla()                        -> the vector derivative, sum_i e_i * d_i f
// - nabla_dot()                    -> inner derivative  (grade-lowering part of nabla)
// - nabla_wdg()                    -> outer derivative  (grade-raising part of nabla)
// - laplacian()                    -> the scalar (grade-preserving) operator sum_i d_i^2
// - d_dt()                         -> the time derivative of a time-dependent field
//
// A FIELD is any callable, f(r) or f(r, t), returning a multivector type. Nothing is
// stored: the operators sample the callable, so analytic and tabulated fields are
// treated alike and no grid, no ownership and no allocation are involved.
//
// The 2d counterpart of ga_ega3d_ops_calculus.hpp; see that header for the reasoning.
// In 2d the grade ladder is shorter, so:
//
//     scalar field  phi : nabla ^ phi = grad phi              (a vector)
//     vector field  A   : nabla . A   = div A                 (a scalar)
//                         nabla ^ A   = the scalar curl       (a pseudoscalar)
//     pseudoscalar  p   : nabla . p                           (a vector)
//
// NAMING -- why these are not called div and curl
//
// The classical trio are single GRADES of the pair above, not operators in their own
// right, so naming the operators after them would be wrong almost everywhere:
//
//     grad = the OUTER derivative of a SCALAR field
//     div  = the INNER derivative of a VECTOR field
//     curl = the DUAL of the outer derivative of a vector field, in 3d only
//
// Three consequences:
//
//   - nabla_dot is a divergence at exactly one grade. On a bivector field it returns a
//     vector, on a trivector field a bivector; neither is a divergence in any classical
//     sense. Likewise nabla_wdg is a gradient only on a scalar field.
//   - nabla_wdg(A) is not curl A but its DUAL: the two carry the same components and
//     different meaning, a bivector (the plane) versus a vector (the axis). The
//     classical vector is dual(nabla_wdg(A)) -- the DUAL, which carries the metric,
//     not the complement (ega cannot distinguish them, sta4ds can) and not I * (...),
//     which is a different map again away from grades 0 and 1.
//   - curl is 3d-only. The outer derivative of a vector field is a pseudoscalar in 2d
//     and a trivector in sta4ds, so the name would mean something different, or
//     nothing, per algebra. nabla_dot / nabla_wdg mean the same in all of them.
//
// A COMMON MISREADING: "gradient = divergence + curl". The identity is
//
//     nabla A = nabla . A + nabla ^ A          for a VECTOR field A
//
// whose left side is the vector derivative, not the gradient, and whose second term is
// the dual of the curl. For a SCALAR field the inner part vanishes identically, so
// there nabla phi IS the gradient and nothing is added to it.
//
// The 2d curl is a pseudoscalar rather than a bivector-valued object, which is the
// usual "curl is a scalar in the plane" statement read off the grade structure.
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// METRIC ASSUMPTION -- read before porting this to another algebra
////////////////////////////////////////////////////////////////////////////////
//
// The vector derivative is properly written with the RECIPROCAL basis,
//
//     nabla = e^mu d_mu ,      e^mu . e_nu = delta^mu_nu
//
// and the implementation below uses the basis vectors e_i directly. That is
// correct here only because the ega metric is the identity, where e^i = e_i.
//
// It does NOT carry over. In sta4ds the metric is {-1,-1,-1,+1}, so
// e^1 = -g1, e^2 = -g2, e^3 = -g3 and e^4 = +g4: writing the sum with the basis
// instead of the reciprocal basis flips the sign of every spacelike term, and
// does so SILENTLY, because the result is still a well-formed multivector. A
// port to STA needs the reciprocal basis first. In PGA it cannot be ported at
// all -- a degenerate direction has no reciprocal vector.
////////////////////////////////////////////////////////////////////////////////


// nabla f at r -- the full vector derivative, as a multivector
template <typename F>
inline auto nabla(F&& f, Vec2d<value_t> const& r, fd_scheme const& sc,
                  value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, nrm(r));

    MVec2d<value_t> acc{};
    for (int i = 0; i < 2; ++i) {
        acc = acc + MVec2d<value_t>(detail::ega2d_basis<value_t>(i)) *
                        detail::ega2d_axis_deriv(f, r, i, sc, h);
    }
    return acc;
}

template <typename F>
inline auto nabla(F&& f, Vec2d<value_t> const& r, value_t h = value_t(0.0))
{
    static fd_scheme const d1 = central_scheme(1, 2);
    return nabla(std::forward<F>(f), r, d1, h);
}

template <typename F>
inline auto nabla(F&& f, Vec2d<value_t> const& r, value_t t, value_t h)
{
    return nabla([&](Vec2d<value_t> const& q) { return f(q, t); }, r, h);
}

// The inner and outer derivatives are the two grade projections of nabla f, deduced
// from the callable's return type. The extreme cases are identically zero (a scalar
// field has no inner derivative, a pseudoscalar field no outer one in 2d).

// inner derivative: nabla . f  (grade lowering)
template <typename F>
inline auto nabla_dot(F&& f, Vec2d<value_t> const& r, value_t h = value_t(0.0))
{
    using ret_t = std::invoke_result_t<F, Vec2d<value_t> const&>;
    auto const N = nabla(std::forward<F>(f), r, h);

    if constexpr (std::is_same_v<ret_t, Scalar2d<value_t>>) {
        return Scalar2d<value_t>(0.0);
    }
    else if constexpr (std::is_same_v<ret_t, Vec2d<value_t>>) {
        return gr0(N); // div
    }
    else if constexpr (std::is_same_v<ret_t, PScalar2d<value_t>>) {
        return gr1(N);
    }
    else {
        return N; // mixed-grade field: the caller projects
    }
}

// outer derivative: nabla ^ f  (grade raising)
template <typename F>
inline auto nabla_wdg(F&& f, Vec2d<value_t> const& r, value_t h = value_t(0.0))
{
    using ret_t = std::invoke_result_t<F, Vec2d<value_t> const&>;
    auto const N = nabla(std::forward<F>(f), r, h);

    if constexpr (std::is_same_v<ret_t, Scalar2d<value_t>>) {
        return gr1(N); // grad
    }
    else if constexpr (std::is_same_v<ret_t, Vec2d<value_t>>) {
        return gr2(N); // the scalar curl, as a pseudoscalar
    }
    else if constexpr (std::is_same_v<ret_t, PScalar2d<value_t>>) {
        return PScalar2d<value_t>(0.0); // no grade 3 in 2d
    }
    else {
        return N;
    }
}

// the Laplacian -- a SCALAR operator, so it preserves the grade of the field
template <typename F>
inline auto laplacian(F&& f, Vec2d<value_t> const& r, fd_scheme const& sc,
                      value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, nrm(r));

    MVec2d<value_t> acc{};
    for (int i = 0; i < 2; ++i) {
        acc = acc + detail::ega2d_axis_deriv(f, r, i, sc, h);
    }
    return acc;
}

template <typename F>
inline auto laplacian(F&& f, Vec2d<value_t> const& r, value_t h = value_t(0.0))
{
    static fd_scheme const d2 = central_scheme(2, 2);
    return laplacian(std::forward<F>(f), r, d2, h);
}

// d f / dt of a time-dependent field
template <typename F>
inline auto d_dt(F&& f, Vec2d<value_t> const& r, value_t t, fd_scheme const& sc,
                 value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, t);

    return apply_scheme<MVec2d<value_t>>(
        sc, [&](value_t off) { return f(r, t + off); }, h);
}

template <typename F>
inline auto d_dt(F&& f, Vec2d<value_t> const& r, value_t t, value_t h = value_t(0.0))
{
    static fd_scheme const d1 = central_scheme(1, 2);
    return d_dt(std::forward<F>(f), r, t, d1, h);
}

} // namespace hd::ga::ega
