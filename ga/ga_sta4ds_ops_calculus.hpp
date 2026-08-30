#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_sta4ds_ops.hpp" // sta4ds ops (products, basics, exp/log, transform)
#include "ga_usr_consts.hpp" // named basis-blade constants
#include "ga_usr_fd.hpp"     // fd_scheme, central/compact schemes, fd_derivative

#include <cmath>       // std::sqrt
#include <type_traits> // std::invoke_result_t, std::is_same_v
#include <utility>     // std::forward


namespace hd::ga::detail {

// Helpers for the sta4ds field operators. Kept in hd::ga::detail (NOT
// hd::ga::sta::detail): opening a detail namespace inside hd::ga::sta shadows the
// hd::ga::detail names that the algebra headers use unqualified, and whether that
// breaks the build then depends on include order.

// mu-th coordinate direction of sta4ds, in the library's order (g1, g2, g3 space,
// g4 time)
template <typename T> inline Vec4ds<T> sta4ds_basis(int mu)
{
    switch (mu) {
        case 0:
            return Vec4ds<T>(T(1.0), T(0.0), T(0.0), T(0.0));
        case 1:
            return Vec4ds<T>(T(0.0), T(1.0), T(0.0), T(0.0));
        case 2:
            return Vec4ds<T>(T(0.0), T(0.0), T(1.0), T(0.0));
        default:
            return Vec4ds<T>(T(0.0), T(0.0), T(0.0), T(1.0));
    }
}

// mu-th RECIPROCAL direction: e^mu . e_nu = delta. For the orthogonal standard basis
// this is inv(e_mu) = e_mu / (e_mu . e_mu), i.e. -g1, -g2, -g3, +g4. Using e_mu here
// instead would silently flip every spacelike term -- see the block in the ops header.
template <typename T> inline Vec4ds<T> sta4ds_rec_basis(int mu)
{
    return (mu == 3) ? sta4ds_basis<T>(mu) : Vec4ds<T>(-sta4ds_basis<T>(mu));
}

// directional derivative of a field along coordinate mu -- apply_scheme() owns the
// weight loop and rejects a scheme that cannot be evaluated at a single point
template <typename F>
inline MVec4ds<value_t> sta4ds_axis_deriv(F&& f, Vec4ds<value_t> const& x, int mu,
                                          fd_scheme const& sc, value_t h)
{
    Vec4ds<value_t> const e = sta4ds_basis<value_t>(mu);
    return apply_scheme<MVec4ds<value_t>>(
        sc, [&](value_t off) { return f(x + off * e); }, h);
}

} // namespace hd::ga::detail


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// provides sta4ds differential operations on multivector FIELDS:
//
// - nabla()                        -> the vector derivative, sum_mu e^mu * d_mu f
// - nabla_dot()                    -> inner derivative  (grade-lowering part of nabla)
// - nabla_wdg()                    -> outer derivative  (grade-raising part of nabla)
// - dalembertian()                 -> the wave operator, nabla . nabla
//
// A FIELD is any callable f(x) of a spacetime position, returning a multivector type.
// Nothing is stored: the operators sample the callable, so analytic and tabulated
// fields are treated alike.
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// the spacetime vector derivative
////////////////////////////////////////////////////////////////////////////////
//
// Unlike its ega counterparts this operator is written over the RECIPROCAL basis,
//
//     nabla = e^mu d_mu ,      e^mu . e_nu = delta^mu_nu
//
// which in the metric {-1,-1,-1,+1} means e^1 = -g1, e^2 = -g2, e^3 = -g3 and
// e^4 = +g4. Writing the sum over the basis instead -- as the ega layers legitimately
// do, since there e^i = e_i -- would flip the sign of every spacelike term, and would
// do so SILENTLY, because the result is still a well-formed multivector. The general
// frame is available as reciprocal_frame(); for this orthogonal basis it reduces to
// inv(e_mu), which is what sta4ds_rec_basis() uses.
//
// TIME IS A COORDINATE HERE. Unlike ega, where time had to be carried as a separate
// argument, x is a spacetime position and d/dt is just d_4. So there is no d_dt():
// the fourth term of the sum is the time derivative.
//
// The grade split is the same as everywhere: for a grade-k field the inner derivative
// is gr_{k-1}(nabla f) and the outer one gr_{k+1}(nabla f). For the Faraday bivector
// F those two parts are the two halves of Maxwell's equation, nabla F = J, with
// nabla . F = J (sources) and nabla ^ F = 0 (no magnetic monopoles).
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
// TESTING A RESULT FOR ZERO: nrm() CANNOT DO IT HERE. In a Minkowski metric the
// norm is null-capable -- nrm(X) == 0 holds for a whole cone of non-zero X -- so
// "nrm(nabla f) < tol" is NOT a check that the field solves its equation. It will
// report success for a field that does not. Compare componentwise instead. This is
// not hypothetical: while writing the tests for this header a plane wave with the
// wrong magnetic pairing gave nabla F with components of order 0.1 and nrm exactly
// 0, and the nrm-based check passed.
//
// The d'Alembertian is nabla . nabla = d_4^2 - d_1^2 - d_2^2 - d_3^2, a scalar
// operator, so it preserves the grade of the field. Note the sign convention: it
// follows the metric, so a wave equation reads box f = 0 with no extra minus.
////////////////////////////////////////////////////////////////////////////////

// nabla f at the spacetime point x -- the full vector derivative
template <typename F>
inline auto nabla(F&& f, Vec4ds<value_t> const& x, fd_scheme const& sc,
                  value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, value_t(nrm(x)));

    MVec4ds<value_t> acc{};
    for (int mu = 0; mu < 4; ++mu) {
        acc = acc + MVec4ds<value_t>(detail::sta4ds_rec_basis<value_t>(mu)) *
                        detail::sta4ds_axis_deriv(f, x, mu, sc, h);
    }
    return acc;
}

template <typename F>
inline auto nabla(F&& f, Vec4ds<value_t> const& x, value_t h = value_t(0.0))
{
    static fd_scheme const d1 = central_scheme(1, 2);
    return nabla(std::forward<F>(f), x, d1, h);
}

// The inner and outer derivatives are the two grade projections of nabla f, deduced
// from the callable's return type. The extreme cases are identically zero: a scalar
// field has no inner derivative, a pseudoscalar field no outer one.

template <typename F>
inline auto nabla_dot(F&& f, Vec4ds<value_t> const& x, value_t h = value_t(0.0))
{
    using ret_t = std::invoke_result_t<F, Vec4ds<value_t> const&>;
    auto const N = nabla(std::forward<F>(f), x, h);

    if constexpr (std::is_same_v<ret_t, Scalar4ds<value_t>>) {
        return Scalar4ds<value_t>(0.0);
    }
    else if constexpr (std::is_same_v<ret_t, Vec4ds<value_t>>) {
        return gr0(N);
    }
    else if constexpr (std::is_same_v<ret_t, BiVec4ds<value_t>>) {
        return gr1(N); // for the Faraday bivector: the source term of nabla F = J
    }
    else if constexpr (std::is_same_v<ret_t, TriVec4ds<value_t>>) {
        return gr2(N);
    }
    else if constexpr (std::is_same_v<ret_t, PScalar4ds<value_t>>) {
        return gr3(N);
    }
    else {
        return N; // mixed-grade field: the caller projects
    }
}

template <typename F>
inline auto nabla_wdg(F&& f, Vec4ds<value_t> const& x, value_t h = value_t(0.0))
{
    using ret_t = std::invoke_result_t<F, Vec4ds<value_t> const&>;
    auto const N = nabla(std::forward<F>(f), x, h);

    if constexpr (std::is_same_v<ret_t, Scalar4ds<value_t>>) {
        return gr1(N);
    }
    else if constexpr (std::is_same_v<ret_t, Vec4ds<value_t>>) {
        return gr2(N);
    }
    else if constexpr (std::is_same_v<ret_t, BiVec4ds<value_t>>) {
        return gr3(N); // for the Faraday bivector: vanishes (no magnetic monopoles)
    }
    else if constexpr (std::is_same_v<ret_t, TriVec4ds<value_t>>) {
        return gr4(N);
    }
    else if constexpr (std::is_same_v<ret_t, PScalar4ds<value_t>>) {
        return PScalar4ds<value_t>(0.0); // no grade 5 in 4d
    }
    else {
        return N;
    }
}

// the d'Alembertian, nabla . nabla = d_4^2 - d_1^2 - d_2^2 - d_3^2. A scalar
// operator, so it preserves the grade of the field.
template <typename F>
inline auto dalembertian(F&& f, Vec4ds<value_t> const& x, fd_scheme const& sc,
                         value_t h = value_t(0.0))
{
    if (h <= value_t(0.0)) h = fd_step(sc, value_t(nrm(x)));

    MVec4ds<value_t> acc{};
    for (int mu = 0; mu < 4; ++mu) {
        // the metric sign of the coordinate: +1 for time, -1 for the space directions
        value_t const g = (mu == 3) ? value_t(1.0) : value_t(-1.0);
        acc = acc + g * detail::sta4ds_axis_deriv(f, x, mu, sc, h);
    }
    return acc;
}

template <typename F>
inline auto dalembertian(F&& f, Vec4ds<value_t> const& x, value_t h = value_t(0.0))
{
    static fd_scheme const d2 = central_scheme(2, 2);
    return dalembertian(std::forward<F>(f), x, d2, h);
}

} // namespace hd::ga::sta
