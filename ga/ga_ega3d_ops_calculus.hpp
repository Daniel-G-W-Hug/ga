#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_ega3d_ops.hpp" // ega3d ops (products, basics, exp/log, rotations)
#include "ga_usr_fd.hpp"    // fd_scheme, central/compact schemes, fd_derivative

#include <cmath>       // std::sqrt, std::cbrt
#include <limits>      // std::numeric_limits
#include <type_traits> // std::invoke_result_t
#include <utility>     // std::forward


namespace hd::ga::detail {

// Helpers for the ega3d field operators. Kept in hd::ga::detail (NOT
// hd::ga::ega::detail): opening a detail namespace inside hd::ga::ega shadows
// the hd::ga::detail names that the algebra headers use unqualified, and whether
// that breaks the build then depends on include order.

// i-th cartesian basis vector of ega3d
template <typename T> inline Vec3d<T> ega3d_basis(int i)
{
    return (i == 0) ? Vec3d<T>(T(1.0), T(0.0), T(0.0))
                    : ((i == 1) ? Vec3d<T>(T(0.0), T(1.0), T(0.0))
                                : Vec3d<T>(T(0.0), T(0.0), T(1.0)));
}

// directional derivative of a field along axis i -- apply_scheme() owns the weight
// loop and rejects a scheme that cannot be evaluated at a single point
template <typename F>
inline MVec3d<value_t> ega3d_axis_deriv(F&& f, Vec3d<value_t> const& r, int i,
                                        fd_scheme const& sc, value_t h)
{
    Vec3d<value_t> const e = ega3d_basis<value_t>(i);
    return apply_scheme<MVec3d<value_t>>(
        sc, [&](value_t off) { return f(r + off * e); }, h);
}

// the cached default schemes (2nd-order central), built once
inline fd_scheme const& ega3d_default_scheme(int deriv)
{
    static fd_scheme const d1 = central_scheme(1, 2);
    static fd_scheme const d2 = central_scheme(2, 2);
    return (deriv == 1) ? d1 : d2;
}


} // namespace hd::ga::detail


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
//     scalar field  phi : nabla phi     = grad phi                 (a vector)
//     vector field  A   : gr0(nabla A)  = div A                    (a scalar)
//                         gr2(nabla A)  = dual(curl A)             (a bivector)
//     bivector field B  : gr1(nabla B)  = -curl of its dual        (a vector)
//                         gr3(nabla B)  = dual(div of its dual)    (a trivector)
//     trivector field T : gr2(nabla T)  = dual(grad of its coeff)  (a bivector)
//
// The map is the DUAL -- dual(), which carries the metric -- and NOT the complement,
// even though ega cannot tell the two apart: with an identity metric cmpl(X) == dual(X)
// at every grade, so either name reads correctly here and only one of them survives a
// port. sta4ds settles it, its metric being {-1,-1,-1,+1}:
//
//     l_cmpl(v) = TriVec4ds(-0.3, 0.5, -0.2, 1.4)
//     l_dual(v) = TriVec4ds( 0.3,-0.5,  0.2, 1.4)      -- different maps
//
// and the identity that relates the pseudoscalar to them, I * X = l_dual(rev(X)),
// holds for the DUAL at every grade and for the complement at none.
//
// Nor is the dual multiplication by the pseudoscalar. I * X = l_dual(rev(X)) carries a
// reversion, so the two agree at grades 0 and 1 and differ in sign from grade 2. Writing
// "I (curl A)" above would happen to be right -- curl A is a vector -- while teaching a
// rule that is wrong one grade over. dual() is an involution here, so the relation reads
// both ways: dual(gr2(nabla A)) = curl A.
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

    return apply_scheme<MVec3d<value_t>>(
        sc, [&](value_t off) { return f(r, t + off); }, h);
}

template <typename F>
inline auto d_dt(F&& f, Vec3d<value_t> const& r, value_t t, value_t h = value_t(0.0))
{
    return d_dt(std::forward<F>(f), r, t, detail::ega3d_default_scheme(1), h);
}

} // namespace hd::ga::ega
