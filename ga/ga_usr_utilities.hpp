#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <algorithm> // std::clamp, std::max, std::min (adaptive step controller)
#include <array>     // std::array (rk4_step vector overload)
#include <cmath>     // std::cos, std::sin
#include <mdspan>    // std::mdspan, std::dextents (used by rk4_step)
#include <numbers>   // math constants like pi
#include <stdexcept> // std::invalid_argument
#include <utility>   // std::pair, std::move (rk4_step vector overload)
#include <vector>    // std::vector (rk4_step vector overload)

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

// step functions mapping x to the range [0.0, 1.0] (e.g. for blending/easing):
// each normalizes x over [low_x, high_x], clamps to [0.0, 1.0] and applies its shape
//
//                  0.                               for x < low_x
// linear_step(x) = (x - low_x)/(high_x - low_x)     for low_x <= x <= high_x
//                  1.                               for x > high_x
//
constexpr value_t linear_step(value_t low_x, value_t high_x, value_t x)
{
    return std::clamp((x - low_x) / (high_x - low_x), 0.0, 1.0);
}

//                  0.                for x < low_x
// smooth_step(x) = 3*x^2 - 2*x^3     for low_x <= x <= high_x
//                  1.                for x > high_x
//
// origin: 3rd order polynomial with df/dx = 0 at low_x and high_x
//
constexpr value_t smooth_step(value_t low_x, value_t high_x, value_t x)
{
    value_t const t = std::clamp((x - low_x) / (high_x - low_x), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

//                    0.                          for x < low_x
// smoother_step(x) = 6*x^5 - 15*x^4 + 10*x^3     for low_x <= x <= high_x
//                    1.                          for x > high_x
//
// origin: 5th order polynomial with df/dx = 0 and d²f/dx² = 0 at low_x and high_x
//
constexpr value_t smoother_step(value_t low_x, value_t high_x, value_t x)
{
    value_t const t = std::clamp((x - low_x) / (high_x - low_x), 0.0, 1.0);
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

// sign function for floating point types
// returns +1.0 for value >= 0.0 and -1.0 for value < 0.0
// HINT: this is intentionally NOT the typical signum function
//       which returns 0.0 for value == 0.0 (that one is signum() below);
//       use sign() where a sign factor must never vanish
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

// signum function for floating point types (the classical three-valued one)
// returns +1.0 for value > 0.0, 0.0 for value == 0.0 and -1.0 for value < 0.0
// HINT: use this where zero is a case of its own, e.g. to classify a signed
//       squared quantity as positive, zero or negative; use sign() above when a
//       sign factor is required that never becomes zero
template <typename T>
    requires(numeric_type<T>)
constexpr T signum(T value)
{
    if (value > 0.0) {
        return 1.0;
    }
    if (value < 0.0) {
        return -1.0;
    }
    return 0.0;
}

// signum function overload for Scalar_t types
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr T signum(Scalar_t<T, Tag> s)
{
    return signum(T(s));
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
void rk4_step(std::mdspan<VecType, std::dextents<size_t, 1>> u,
              std::mdspan<VecType, std::dextents<size_t, 2>> uh,
              std::mdspan<VecType const, std::dextents<size_t, 1>> rhs, value_t const dt,
              size_t rk_step)
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


////////////////////////////////////////////////////////////////////////////////
// rk4_step — std::vector overload (Python-bindable companion to the mdspan
// version above).
//
// Same algorithm and identical numerical results as the mdspan version;
// the differences are purely API ergonomics:
//
//   - Inputs (u, uh) are taken by value and the updated state is returned
//     as `std::pair<u, uh>`. In C++ this is essentially zero-cost via move
//     semantics + NRVO; from Python it matches the natural pattern
//     `u, uh = rk4_step(u, uh, rhs, dt, rk_step)` (nanobind has no caster
//     for std::mdspan, so the original signature is unreachable from
//     Python — see ga_py/README.md).
//
//   - `uh` is shaped as `std::array<std::vector<VecType>, 2>` mirroring the
//     `[2 x n]` layout of the mdspan version. From Python this surfaces as
//     a 2-tuple of lists.
//
// Use the mdspan overload for in-place hot loops in C++; use this overload
// when calling from Python or when a value-style API reads cleaner.
////////////////////////////////////////////////////////////////////////////////

template <typename VecType>
std::pair<std::vector<VecType>, std::array<std::vector<VecType>, 2>>
rk4_step(std::vector<VecType> u, std::array<std::vector<VecType>, 2> uh,
         std::vector<VecType> const& rhs, value_t const dt, size_t rk_step)
{
    if (rk_step < 1 || rk_step > 4) {
        throw std::invalid_argument("rk4_step requires argument: 1 <= rk_step <= 4.");
    }
    size_t const n = u.size();
    if (uh[0].size() != n || uh[1].size() != n || rhs.size() != n) {
        throw std::invalid_argument(
            "rk4_step: u, uh[0], uh[1] and rhs must all have the same length.");
    }

    value_t const rk1 = 1.0 / 6.0 * dt;
    value_t const rk2 = 1.0 / 3.0 * dt;
    value_t const rk3 = 1.0 / 2.0 * dt;
    value_t const rk4 = dt;

    switch (rk_step) {
        case 1: // predictor 1: Euler forward to t + 0.5*dt
            for (size_t i = 0; i < n; ++i) {
                uh[0][i] = u[i];
            }
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0][i] + rk3 * rhs[i];
                uh[1][i] = rk1 * rhs[i];
            }
            break;

        case 2: // corrector 1: Euler backward to t + 0.5*dt
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0][i] + rk3 * rhs[i];
                uh[1][i] += rk2 * rhs[i];
            }
            break;

        case 3: // predictor 2: midpoint rule to t + dt
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0][i] + rk4 * rhs[i];
                uh[1][i] += rk2 * rhs[i];
            }
            break;

        case 4: // corrector 2: Simpson rule to t + dt
            for (size_t i = 0; i < n; ++i) {
                u[i] = uh[0][i] + uh[1][i] + rk1 * rhs[i];
            }
            break;
    }

    return {std::move(u), std::move(uh)};
}

////////////////////////////////////////////////////////////////////////////////
// Uniform ODE integrators with a callable right-hand side
//
// Both classes integrate  du/dt = f(t, u)  for a flat state std::vector<double> via a
// single  step(f, u, t, dt) -> t + dt  call that advances `u` in place. `f` is any
// callable
//
//     void f(double t, std::vector<double> const& u, std::vector<double>& dudt)
//
// writing the derivative into the provided buffer (no allocation in the hot loop). The
// uniform interface lets a caller swap integrators and benchmark them head to head; it is
// also the seam dynamic_system uses to select RK4 vs the multistep method.
//
// rk4_integrator wraps the canonical substage-based rk4_step (above) -- it is NOT a
// second RK4 implementation -- so the classic 4th-order, 4-evaluations-per-step method
// stays the single source of truth.
////////////////////////////////////////////////////////////////////////////////
class rk4_integrator {

  public:

    explicit rk4_integrator(size_t n) : uh_(2 * n), rhs_(n) {}

    // advance u from t to t + dt in place (4 rhs evaluations); returns t + dt.
    template <typename RHS>
    double step(RHS&& f, std::vector<double>& u, double t, double dt)
    {
        size_t const n = u.size();
        auto us = std::mdspan<double, std::dextents<size_t, 1>>(u.data(), n);
        auto uhs = std::mdspan<double, std::dextents<size_t, 2>>(uh_.data(), 2, n);
        auto rs = std::mdspan<double const, std::dextents<size_t, 1>>(rhs_.data(), n);
        for (size_t s = 1; s <= 4; ++s) {
            // rhs at the current substage state, evaluated at the substage time
            f(rk4_get_time(t, dt, s - 1), u, rhs_);
            rk4_step(us, uhs, rs, dt, s);
        }
        return t + dt;
    }

  private:

    std::vector<double> uh_;  // [2 x n] RK4 scratch (flattened)
    std::vector<double> rhs_; // per-substage derivative buffer
};

////////////////////////////////////////////////////////////////////////////////
// Adams-Bashforth-Moulton 2nd-order predictor-corrector, fixed step (PECE)
//
//   Predict  (Adams-Bashforth 2):  u_p     = u_n + dt (3/2 f_n - 1/2 f_{n-1})
//   Evaluate                       f_p     = f(t_{n+1}, u_p)
//   Correct  (Adams-Moulton 2):    u_{n+1} = u_n + dt (1/2 f_p + 1/2 f_n)
//   Evaluate                       f_{n+1} = f(t_{n+1}, u_{n+1})   [reused as f_n next
//   step]
//
// A linear MULTISTEP method: 2 rhs evaluations per step (vs RK4's 4) but 2nd order (vs
// RK4's 4th), keeping one past derivative f_{n-1}. The first step self-starts with one
// RK4 step (reusing rk4_integrator). This is the fixed-dt baseline; an adaptive-dt
// variant layers the predictor-corrector difference as a local error estimate on top of
// the same formulas.
////////////////////////////////////////////////////////////////////////////////
class abm2_integrator {

  public:

    explicit abm2_integrator(size_t n) : rk4_(n), f_n_(n), f_nm1_(n), u_p_(n), f_p_(n) {}

    // advance u from t to t + dt in place (2 rhs evaluations after self-start); returns
    // t + dt. The first call self-starts with a single RK4 step.
    template <typename RHS>
    double step(RHS&& f, std::vector<double>& u, double t, double dt)
    {
        size_t const n = u.size();
        if (!started_) {
            f(t, u, f_nm1_);        // f_0, kept as f_{n-1} for the next step
            rk4_.step(f, u, t, dt); // self-start: u_0 -> u_1 via the canonical RK4
            f(t + dt, u, f_n_);     // f_1
            started_ = true;
            return t + dt;
        }
        for (size_t i = 0; i < n; ++i) // P: Adams-Bashforth 2 predictor
            u_p_[i] = u[i] + dt * (1.5 * f_n_[i] - 0.5 * f_nm1_[i]);
        f(t + dt, u_p_, f_p_);         // E
        for (size_t i = 0; i < n; ++i) // C: Adams-Moulton 2 (trapezoidal) corrector
            u[i] += dt * (0.5 * f_p_[i] + 0.5 * f_n_[i]);
        f_nm1_.swap(f_n_);  // f_{n-1} <- f_n
        f(t + dt, u, f_n_); // E: f_{n+1}, reused as f_n on the next step
        return t + dt;
    }

    void reset() { started_ = false; } // restart the multistep history
    size_t dim() const { return f_n_.size(); }

  private:

    rk4_integrator rk4_; // self-start for the first step
    bool started_{false};
    std::vector<double> f_n_, f_nm1_; // f(t_n, u_n), f(t_{n-1}, u_{n-1})
    std::vector<double> u_p_, f_p_;   // predictor state + its derivative
};

////////////////////////////////////////////////////////////////////////////////
// Adams-Bashforth-Moulton 2nd-order with ADAPTIVE (variable) step size
//
// Same predictor-corrector as abm2_integrator, but the step is chosen automatically from
// a LOCAL ERROR ESTIMATE -- the Milne device: the difference between the (explicit)
// predictor and the (implicit) corrector estimates the local truncation error,
//
//     err ~ (1/6) |u_corrector - u_predictor|.
//
// Because the steps are non-uniform, the predictor uses the VARIABLE-step
// Adams-Bashforth-2 coefficients with r = h_n / h_{n-1}:
//
//     u_p = u_n + h_n [ (1 + r/2) f_n - (r/2) f_{n-1} ]    (-> 3/2, -1/2 when r = 1),
//
// while the trapezoidal corrector  u_{n+1} = u_n + (h_n/2)(f_p + f_n)  is step-ratio
// invariant. The scaled error (per-component  atol + rtol*|u|, infinity norm) drives a
// standard step controller  h_new = h * clamp(0.9 * err^(-1/3), ...); a step with err > 1
// is rejected and retried smaller. First step self-starts with RK4.
//
// NOTE: this is an EXPLICIT method -- adaptivity gives error control + efficiency on
// problems with VARYING dynamics, it is NOT a stiff-system remedy (a stiff problem forces
// the step to the stability limit regardless). For stiffness use an implicit method.
////////////////////////////////////////////////////////////////////////////////
class abm2_adaptive_integrator {

  public:

    explicit abm2_adaptive_integrator(size_t n, double dt_min = 1.0e-12,
                                      double dt_max = 1.0e30) :
        rk4_(n), f_n_(n), f_nm1_(n), u_p_(n), u_c_(n), f_tmp_(n), dt_min_(dt_min),
        dt_max_(dt_max)
    {
    }

    // Take one ACCEPTED adaptive step from t; u is updated in place. `dt` is in/out: the
    // step to attempt on entry, the suggested next step on exit. atol/rtol set the
    // per-component tolerance. Internally rejects + retries with a smaller step until err
    // <= 1 (or dt_min). Returns t + (accepted step).
    template <typename RHS>
    double step(RHS&& f, std::vector<double>& u, double t, double& dt, double atol,
                double rtol)
    {
        size_t const n = u.size();
        if (!started_) {
            f(t, u, f_nm1_);        // f_0
            rk4_.step(f, u, t, dt); // self-start u_0 -> u_1 (fixed dt)
            f(t + dt, u, f_n_);     // f_1
            h_prev_ = dt;
            last_dt_ = dt;
            started_ = true;
            ++accepted_;
            return t + dt;
        }
        for (;;) {
            double const r = dt / h_prev_;
            for (size_t i = 0; i < n; ++i) // variable-step AB2 predictor
                u_p_[i] = u[i] + dt * ((1.0 + 0.5 * r) * f_n_[i] - 0.5 * r * f_nm1_[i]);
            f(t + dt, u_p_, f_tmp_);       // E (f at predictor)
            for (size_t i = 0; i < n; ++i) // trapezoidal corrector (ratio-invariant)
                u_c_[i] = u[i] + 0.5 * dt * (f_tmp_[i] + f_n_[i]);

            double err = 0.0; // scaled local error estimate (Milne 1/6), infinity norm
            for (size_t i = 0; i < n; ++i) {
                double const sc =
                    atol + rtol * std::max(std::abs(u[i]), std::abs(u_c_[i]));
                err = std::max(err, std::abs(u_c_[i] - u_p_[i]) / (6.0 * sc));
            }

            if (err <= 1.0 || dt <= dt_min_ * (1.0 + 1.0e-9)) { // accept
                for (size_t i = 0; i < n; ++i)
                    u[i] = u_c_[i];
                f(t + dt, u, f_tmp_); // f_{n+1}, reused as f_n next step
                f_nm1_.swap(f_n_);
                f_n_.swap(f_tmp_);
                h_prev_ = dt;
                last_dt_ = dt;
                double const t_new = t + dt;
                double const fac =
                    (err > 0.0) ? 0.9 * std::pow(1.0 / err, 1.0 / 3.0) : 5.0;
                dt = std::clamp(dt * std::clamp(fac, 0.2, 5.0), dt_min_, dt_max_);
                ++accepted_;
                return t_new;
            }
            ++rejected_; // reject: shrink and retry
            double const fac = 0.9 * std::pow(1.0 / err, 1.0 / 3.0);
            dt = std::max(dt * std::clamp(fac, 0.1, 1.0), dt_min_);
        }
    }

    // Integrate from t0 to t_end starting with dt0; the final step is clamped to land
    // exactly on t_end. Returns t_end. accepted()/rejected()/last_dt() report the run.
    template <typename RHS>
    double integrate(RHS&& f, std::vector<double>& u, double t0, double t_end, double dt0,
                     double atol, double rtol)
    {
        double t = t0, dt = dt0;
        while (t < t_end - 1.0e-12 * std::max(1.0, std::abs(t_end))) {
            double dt_io = std::min(dt, t_end - t); // clamp to land on t_end
            t = step(f, u, t, dt_io, atol, rtol);
            dt = dt_io;
        }
        return t;
    }

    void reset()
    {
        started_ = false;
        accepted_ = rejected_ = 0;
    }
    size_t accepted() const { return accepted_; }
    size_t rejected() const { return rejected_; }
    double last_dt() const { return last_dt_; }

  private:

    rk4_integrator rk4_; // self-start for the first step
    bool started_{false};
    std::vector<double> f_n_, f_nm1_;       // f(t_n, u_n), f(t_{n-1}, u_{n-1})
    std::vector<double> u_p_, u_c_, f_tmp_; // predictor / corrector / scratch derivative
    double h_prev_{0.0};                    // last accepted step (for the step ratio r)
    double last_dt_{0.0};                   // last accepted step (diagnostic)
    double dt_min_, dt_max_;
    size_t accepted_{0}, rejected_{0};
};

} // namespace hd::ga