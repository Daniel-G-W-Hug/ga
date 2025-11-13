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

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "mdspan/mdspan.hpp"
// make mdspan less verbose
using namespace Kokkos;

#include "ga/ga_pga.hpp"
using namespace hd::ga;

namespace hd {

double get_time_rkstep(double ti, double dt, size_t rk);

// Templatized to support both vec2d, vec2dp, vec3dp, etc.
template<typename VecType>
void rk4_step(mdspan<VecType, dextents<size_t, 1>> u,
              mdspan<VecType, dextents<size_t, 2>> uh,
              mdspan<VecType const, dextents<size_t, 1>> rhs,
              double const dt,
              size_t rk_step);

} // namespace hd

// Include implementation since it's a template
#include "rk4_integration_impl.hpp"

#endif // HD_RK4_INTEGRATION_H
