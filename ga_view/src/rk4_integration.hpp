#ifndef HD_RK4_INTEGRATION_H
#define HD_RK4_INTEGRATION_H
////////////////////////////////////////////////////////////////////////////////
// time integration with Runge-Kutta 4-step method of 4th order
//
// du/dt = rhs(t,u(t)), u(t0)=u0
//
// with:
//  k1 = rhs(t_i       , u_i)             for Pred. 1: Euler forward to t+0.5*dt
//  k2 = rhs(t_i + dt/2, u_i + dt/2*k1)   for Corr. 1: Euler backward to t+0.5*dt
//  k3 = rhs(t_i + dt/2, u_i + dt/2*k2)   for Pred. 2: Midpoint rule to t+dt
//  k4 = rhs(t_i + dt  , u_i + dt*k3)     for Corr. 2: Simpson rule to t+dt
//
// u_i+1 = u_i + dt/6*(k1 + 2*k2 + 2*k3 + k4)
//
// Adapted for ga_view from experiments/ode_int
////////////////////////////////////////////////////////////////////////////////

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "mdspan/mdspan.hpp"
// make mdspan less verbose
using namespace Kokkos;

#include "ga/ga_pga.hpp"
using namespace hd::ga;

#include <stdexcept> // std::invalid_argument

namespace hd {

// Get time at RK sub-step (inline in header)
inline double get_time_rkstep(double ti, double dt, size_t rk)
{
    switch (rk) {
        case 0:
            return ti;
        case 1:
            [[fallthrough]];
        case 2:
            return ti + dt / 2.;
        case 3:
            return ti + dt;
        default:
            throw std::invalid_argument(
                "hd::get_time_rkstep: rk out of range (0 <= rk < 4).");
    }
}

// Templatized RK4 step to support vec2d, vec2dp, vec3d, vec3dp, etc.
// This allows modeling dynamic systems in different geometric algebras (EGA, PGA)
template <typename VecType>
void rk4_step(mdspan<VecType, dextents<size_t, 1>> u,
              mdspan<VecType, dextents<size_t, 2>> uh,
              mdspan<VecType const, dextents<size_t, 1>> rhs, double const dt,
              size_t rk_step)
{
    double rk1 = 1. / 6. * dt;
    double rk2 = 1. / 3. * dt;
    double rk3 = 1. / 2. * dt;
    double rk4 = dt;

    size_t n = u.extent(0);

    switch (rk_step) {
        case 1: // predictor 1: Euler forward to t + 0.5*dt
            for (size_t i = 0; i < n; ++i) {
                uh[0, i] = u[i];
            }
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0, i] + rk3 * rhs[i];
                uh[1, i] = rk1 * rhs[i];
            }
            break;

        case 2: // corrector 1: Euler backward to t + 0.5*dt
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0, i] + rk3 * rhs[i];
                uh[1, i] += rk2 * rhs[i];
            }
            break;

        case 3: // predictor 2: midpoint rule to t + dt
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0, i] + rk4 * rhs[i];
                uh[1, i] += rk2 * rhs[i];
            }
            break;

        case 4: // corrector 2: Simpson rule to t + dt
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0, i] + uh[1, i] + rk1 * rhs[i];
            }
            break;
    }
}

} // namespace hd

#endif // HD_RK4_INTEGRATION_H
