#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <cmath>     // std::cos, std::sin
#include <mdspan>    // std::mdspan, std::dextents (used by rk4_step)
#include <numbers>   // math constants like pi
#include <stdexcept> // std::invalid_argument

#include "detail/type_t/ga_scalar_t.hpp"
#include "ga_value_t.hpp"

// Bring std::mdspan / std::dextents into the global namespace so that
// existing call sites (`mdspan<...> u(...)`) keep compiling unchanged
// after the Kokkos-fallback shim was removed.
using std::dextents;
using std::mdspan;

namespace hd::ga {

using std::numbers::pi; // make pi available for users

// conversion between degrees and radians [360° = 2*pi rad]
constexpr value_t deg2rad(value_t value) { return value / 180.0 * pi; }
constexpr value_t rad2deg(value_t value) { return value * 180.0 / pi; }

// conversion between revolutions per minute [rpm] and radians per second [rad/s = radps]
constexpr value_t rpm2radps(value_t value) { return value * 2.0 * pi / 60.0; };
constexpr value_t radps2rpm(value_t value) { return value * 60.0 / (2.0 * pi); };

// conversion between revolutions per second (rps) = Hertz [1 rps = 1 Hz = 1 1/s^(-1)] and
// radians per second [rad/s = radps]
constexpr value_t Hz2radps(value_t value) { return value * 2.0 * pi; };
constexpr value_t radps2Hz(value_t value) { return value / (2.0 * pi); };

// sign function for floating point types
// returns +1.0 for value >= 0.0 and -1.0 for value < 0.0
// HINT: this is intentionally NOT the typical signum function
//       which returns 0.0 for value == 0.0
template <typename T>
    requires(numeric_type<T>)
constexpr T sign(T value)
{
    if (value >= 0.0) {
        return 1.0;
    }
    return -1.0;
}

// sign function overload for Scalar_t types
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr T sign(Scalar_t<T, Tag> s)
{
    return sign(T(s));
}

// Templates for is_even and is_odd work with any integer type
template <typename T>
    requires(std::is_integral_v<T>) // is_even: argument must be an integral type
bool is_even(T arg)
{
    return arg % 2 == 0;
}

template <typename T>
    requires(std::is_integral_v<T>) // is_odd: argument must be an integral type
bool is_odd(T arg)
{
    return arg % 2 != 0;
}

template <typename T>
    requires(numeric_type<T>)
class CosSin {

  public:

    // ctor
    CosSin() = delete; // always require an argument for initialization
    CosSin(T const& x) : m_c(std::cos(x)), m_s(std::sin(x)) {}

    // simple getter (no setter, this must be done via the ctor)
    T const& c() const { return m_c; }
    T const& s() const { return m_s; }

  private:

    T m_c;
    T m_s;
};

class discrete_range {
  public:

    discrete_range(double min, double max, size_t steps) :
        m_min(min), m_max(max), m_steps(steps), m_delta((m_max - m_min) / m_steps)
    {
        if (max < min)
            throw std::invalid_argument("discrete_range: requires min <= max.");
    }

    double min() const { return m_min; }
    double max() const { return m_max; }
    double delta() const { return m_delta; }
    size_t steps() const { return m_steps; }

  private:

    // user input in ctor
    double m_min;
    double m_max;
    size_t m_steps;

    // computed value
    double m_delta;
};

////////////////////////////////////////////////////////////////////////////////
// RK4 Integration (Runge-Kutta 4-step method of 4th order)
//
// Solves: du/dt = rhs(t, u(t)), u(t0) = u0
//
// with:
//  k1 = rhs(t_i       , u_i)             for Pred. 1: Euler forward to t+0.5*dt
//  k2 = rhs(t_i + dt/2, u_i + dt/2*k1)   for Corr. 1: Euler backward to t+0.5*dt
//  k3 = rhs(t_i + dt/2, u_i + dt/2*k2)   for Pred. 2: Midpoint rule to t+dt
//  k4 = rhs(t_i + dt  , u_i + dt*k3)     for Corr. 2: Simpson rule to t+dt
//
// u_i+1 = u_i + dt/6*(k1 + 2*k2 + 2*k3 + k4)
//
// The template parameter VecType supports any GA vector type (vec2d, vec2dp,
// vec3d, vec3dp, etc.) allowing physics simulations in different algebras.
////////////////////////////////////////////////////////////////////////////////

// Get time at RK sub-step
// rk_step: 0 -> t_i, 1 or 2 -> t_i + dt/2, 3 -> t_i + dt
inline value_t rk4_get_time(value_t ti, value_t dt, size_t rk_step)
{
    switch (rk_step) {
        case 0:
            return ti;
        case 1:
            [[fallthrough]];
        case 2:
            return ti + dt / 2.0;
        case 3:
            return ti + dt;
        default:
            throw std::invalid_argument(
                "hd::ga::rk4_get_time: rk_step out of range (0 <= rk_step < 4).");
    }
}

// RK4 integration step - templatized to support any GA vector type
// Parameters:
//   u   - state vector (modified in place)
//   uh  - helper storage (2 x n array for intermediate values)
//   rhs - right-hand side values (computed externally based on current state)
//   dt  - time step size
//   rk_step - RK sub-step (1, 2, 3, or 4)
template <typename VecType>
void rk4_step(
    std::mdspan<VecType, std::dextents<size_t, 1>> u,
    std::mdspan<VecType, std::dextents<size_t, 2>> uh,
    std::mdspan<VecType const, std::dextents<size_t, 1>> rhs,
    value_t const dt, size_t rk_step)
{

    if (rk_step < 1 || rk_step > 4)
        throw std::invalid_argument("rk4_step requires argument: 1 <= rk_step <= 4.");

    value_t const rk1 = 1.0 / 6.0 * dt;
    value_t const rk2 = 1.0 / 3.0 * dt;
    value_t const rk3 = 1.0 / 2.0 * dt;
    value_t const rk4 = dt;

    size_t const n = u.extent(0);

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

} // namespace hd::ga