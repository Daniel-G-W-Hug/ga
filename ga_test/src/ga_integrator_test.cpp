// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Standalone test/benchmark target for the ODE integrators (ga/ga_usr_utilities.hpp):
// RK4 (canonical 4th-order) vs ABM2 (Adams-Bashforth-Moulton 2nd-order
// predictor-corrector).
//
// This is a NORMAL doctest executable -- it respects CMAKE_BUILD_TYPE, so the timing
// cases are meaningful in Debug AND Release (unlike the forced-O3 benchmarks in
// utilities/). Build and run it in both configs to gauge the debug-vs-release impact; the
// timing case prints the active build mode. The same binary builds and runs on Windows
// for cross-checking.
//
// Phase D.2d of the wafer-grinding plan (TODO/grinding.md): SIMPLE analytic cases UPFRONT
// (the damped harmonic oscillator with a closed form) document the integrators' accuracy
// and cost BEFORE the method is applied to the grinding loop and the full Tao model.
// FIXED dt here; an adaptive-dt variant is layered on top later.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "fmt/format.h"

#include "ga/ga_ega.hpp" // brings in value_t + ga/ga_usr_utilities.hpp (the integrators)

using hd::ga::abm2_adaptive_integrator;
using hd::ga::abm2_integrator;
using hd::ga::rk4_integrator;

// ---------------------------------------------------------------------------------
// Reference problem: the damped harmonic oscillator  m x'' + c x' + k x = 0,
// state u = [x, v], with the closed-form underdamped solution for x(0)=x0, v(0)=0.
// ---------------------------------------------------------------------------------
namespace {

struct oscillator {
    double m, k, c;
    double wn() const { return std::sqrt(k / m); }
    double zeta() const { return c / (2.0 * std::sqrt(k * m)); }
    double wd() const { return wn() * std::sqrt(1.0 - zeta() * zeta()); }

    // rhs: du/dt = [v, -(c/m) v - (k/m) x]
    void operator()(double, std::vector<double> const& u, std::vector<double>& du) const
    {
        du[0] = u[1];
        du[1] = -(c / m) * u[1] - (k / m) * u[0];
    }

    // closed-form x(t) for x0, v0 = 0
    double x_exact(double t, double x0) const
    {
        double const z = zeta(), w = wn(), d = wd();
        return std::exp(-z * w * t) *
               (x0 * std::cos(d * t) + (z * w * x0 / d) * std::sin(d * t));
    }
};

// integrate to time T at fixed dt, returning the max |x - x_exact| over the trajectory
// and the final x. Integ is any integrator with step(f, u, t, dt).
template <typename Integ>
std::pair<double, double> run(oscillator const& osc, Integ& integ, double dt, double T,
                              double x0)
{
    std::vector<double> u{x0, 0.0};
    double t = 0.0, maxerr = 0.0;
    size_t const N = size_t(T / dt + 0.5);
    for (size_t i = 0; i < N; ++i) {
        t = integ.step(osc, u, t, dt);
        maxerr = std::max(maxerr, std::abs(u[0] - osc.x_exact(t, x0)));
    }
    return {maxerr, u[0]};
}

} // namespace

TEST_SUITE("integrators: RK4 vs ABM2 (fixed dt)")
{

    TEST_CASE("damped oscillator: accuracy vs closed form + RK4/ABM2 agreement")
    {
        fmt::println("");
        fmt::println("integrators: damped oscillator accuracy (fixed dt)");
        fmt::println("");

        oscillator const osc{1.0, 100.0, 1.0}; // wn=10, zeta=0.05 (underdamped)
        double const x0 = 1.0, dt = 1.0e-3, T = 10.0;

        rk4_integrator rk(2);
        abm2_integrator ab(2);
        auto const [err_rk, xf_rk] = run(osc, rk, dt, T, x0);
        auto const [err_ab, xf_ab] = run(osc, ab, dt, T, x0);

        fmt::println("  wn = {:.3f}, zeta = {:.3f}, dt = {:.0e}, T = {:.1f}", osc.wn(),
                     osc.zeta(), dt, T);
        fmt::println("  max|x - x_exact|:  RK4 = {:.3e},  ABM2 = {:.3e}", err_rk, err_ab);
        fmt::println(
            "  final x:           RK4 = {:.6f},  ABM2 = {:.6f}  (|diff| = {:.3e})", xf_rk,
            xf_ab, std::abs(xf_rk - xf_ab));

        // RK4 (4th order) is accurate to roundoff-ish at this dt; ABM2 (2nd order) is
        // less accurate but bounded. Both track the closed form.
        CHECK(err_rk < 1.0e-5);
        CHECK(err_ab < 1.0e-2);
        // the two solutions agree to <= 1% of the amplitude (the Phase D.2d gate)
        CHECK(std::abs(xf_rk - xf_ab) < 0.01 * x0);
        fmt::println("");
    }

    TEST_CASE("convergence order: RK4 ~ 4th, ABM2 ~ 2nd")
    {
        fmt::println("");
        fmt::println("integrators: convergence order (error ratio on dt halving)");
        fmt::println("");

        oscillator const osc{1.0, 100.0, 1.0};
        double const x0 = 1.0, T = 1.0; // smooth span, little decay

        auto err_at = [&](double dt) {
            rk4_integrator rk(2);
            abm2_integrator ab(2);
            return std::pair{run(osc, rk, dt, T, x0).first,
                             run(osc, ab, dt, T, x0).first};
        };
        auto const [e_rk_1, e_ab_1] = err_at(2.0e-3);
        auto const [e_rk_2, e_ab_2] = err_at(1.0e-3);

        double const r_rk = e_rk_1 / e_rk_2; // ~16 for 4th order
        double const r_ab = e_ab_1 / e_ab_2; // ~4 for 2nd order
        fmt::println("  error ratio on dt halving:  RK4 = {:.2f} (expect ~16),  "
                     "ABM2 = {:.2f} (expect ~4)",
                     r_rk, r_ab);

        // robust lower bounds: ABM2 is clearly > 1st order, RK4 clearly > 2nd order
        CHECK(r_ab > 3.0);
        CHECK(r_rk > 8.0);
        fmt::println("");
    }

    // OBSERVED (Apple M-series, clang): the Debug-vs-Release contrast is the point.
    //   Debug   (-O0): RK4 ~813 ms, ABM2 ~39 ms  -> ratio ~0.05 (RK4 ~20x slower)
    //   Release (-O3): RK4  ~49 ms, ABM2 ~21 ms  -> ratio ~0.43 (RK4  ~2x slower)
    // RK4 wraps the mdspan-based rk4_step, whose per-element indexing is NOT inlined at
    // -O0, so in Debug the mdspan overhead dwarfs everything and the comparison is
    // meaningless. At -O3 it inlines away and RK4 lands at ~2x ABM2 -- the true
    // eval-count ratio (4 vs 2 rhs evaluations per step). LESSON: only trust the
    // integrator timing in Release; ABM2 wins ~2x per step there, but is 2nd order (vs
    // RK4's 4th), so for a target accuracy it needs a smaller dt -- the real trade-off
    // the grinding loop will weigh. (Numbers are machine-specific; the qualitative
    // Debug/Release split is the durable result.)
    TEST_CASE("timing: RK4 vs ABM2 wall-clock (Debug vs Release)")
    {
        fmt::println("");
        fmt::println("integrators: timing (build mode = {})",
#ifdef NDEBUG
                     "Release / NDEBUG"
#else
                     "Debug"
#endif
        );
        fmt::println("");

        oscillator const osc{1.0, 100.0, 1.0};
        double const x0 = 1.0, dt = 1.0e-4, T = 100.0; // ~1e6 steps
        size_t const N = size_t(T / dt + 0.5);

        auto time_run = [&](auto make_integ) {
            auto integ = make_integ();
            std::vector<double> u{x0, 0.0};
            double t = 0.0;
            auto const t0 = std::chrono::steady_clock::now();
            for (size_t i = 0; i < N; ++i)
                t = integ.step(osc, u, t, dt);
            auto const t1 = std::chrono::steady_clock::now();
            double const ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            return std::pair{ms, u[0]};
        };

        auto const [ms_rk, xf_rk] = time_run([] { return rk4_integrator(2); });
        auto const [ms_ab, xf_ab] = time_run([] { return abm2_integrator(2); });

        fmt::println("  {} steps, dt = {:.0e}:", N, dt);
        fmt::println("    RK4  (4 evals/step): {:>8.2f} ms  ({:.2f} ns/step)", ms_rk,
                     ms_rk * 1e6 / double(N));
        fmt::println("    ABM2 (2 evals/step): {:>8.2f} ms  ({:.2f} ns/step)", ms_ab,
                     ms_ab * 1e6 / double(N));
        fmt::println("    ABM2/RK4 time ratio: {:.3f}", ms_ab / ms_rk);
        fmt::println("  (run this target in Debug and Release to compare the impact)");

        // sanity: both integrators actually advanced to a finite, decayed state
        CHECK(std::isfinite(xf_rk));
        CHECK(std::isfinite(xf_ab));
        CHECK(std::abs(xf_rk) < x0); // damped -> below the initial amplitude
        fmt::println("");
    }

} // TEST_SUITE("integrators: RK4 vs ABM2 (fixed dt)")

// =================================================================================
// Adaptive (variable-dt) ABM2 -- Phase D.2d-3.
//
// The step size is chosen from the predictor-corrector difference (Milne local error
// estimate). Two analytic cases, increasing in complexity: (1) error control on the
// damped oscillator -- tighter tolerance => smaller error + more steps; (2) step
// ADAPTATION on a localized feature y = e^{-t^2} (slow-fast-slow) -- the step shrinks
// through the steep middle and grows at the flat ends. Adaptivity is error control +
// efficiency on varying dynamics; it is NOT a stiffness remedy (see D.2d-4).
// =================================================================================
TEST_SUITE("integrators: ABM2 adaptive (variable dt)")
{

    TEST_CASE("adaptive ABM2: error control on the damped oscillator")
    {
        fmt::println("");
        fmt::println("integrators: adaptive ABM2 error control (damped oscillator)");
        fmt::println("");

        oscillator const osc{1.0, 100.0, 1.0};
        double const x0 = 1.0, T = 10.0;

        // local adaptive driver: integrate 0->T, recording max error, step counts, dt
        // range
        auto run = [&](double atol, double rtol) {
            abm2_adaptive_integrator integ(2);
            std::vector<double> u{x0, 0.0};
            double t = 0.0, dt = 1.0e-3, maxerr = 0.0, dtmin = 1.0e30, dtmax = 0.0;
            while (t < T - 1.0e-12) {
                double dt_io = std::min(dt, T - t);
                t = integ.step(osc, u, t, dt_io, atol, rtol);
                dt = dt_io;
                maxerr = std::max(maxerr, std::abs(u[0] - osc.x_exact(t, x0)));
                dtmin = std::min(dtmin, integ.last_dt());
                dtmax = std::max(dtmax, integ.last_dt());
            }
            return std::tuple{maxerr, integ.accepted(), integ.rejected(), dtmin, dtmax};
        };
        auto const [e4, n4, r4, mn4, mx4] = run(1.0e-9, 1.0e-4);
        auto const [e7, n7, r7, mn7, mx7] = run(1.0e-12, 1.0e-7);

        fmt::println("  rtol 1e-4: max err = {:.2e}, steps = {} (+{} rej), dt in "
                     "[{:.2e}, {:.2e}]",
                     e4, n4, r4, mn4, mx4);
        fmt::println("  rtol 1e-7: max err = {:.2e}, steps = {} (+{} rej), dt in "
                     "[{:.2e}, {:.2e}]",
                     e7, n7, r7, mn7, mx7);

        CHECK(e4 < 1.0e-2); // error is controlled / bounded
        CHECK(e7 < e4);     // tighter tolerance -> smaller error
        CHECK(n7 > n4);     // tighter tolerance -> more (smaller) steps
        fmt::println("");
    }

    TEST_CASE("adaptive ABM2: step size adapts to a localized feature (e^{-t^2})")
    {
        fmt::println("");
        fmt::println("integrators: adaptive ABM2 step adaptation (y = e^{{-t^2}})");
        fmt::println("");

        // y' = -2 t y -> y = e^{-t^2}: a slow-fast-slow profile. Under ABSOLUTE error
        // control (atol, rtol = 0) the step shrinks where |y'''| peaks -- the visually
        // steep middle
        // (~t in [0.6, 0.9]) -- and grows at the flat ends where |y'''| -> 0.
        // (Relative/rtol control would instead demand the tail, where y'''/y ~ t^3 grows
        // even as y -> 0.)
        auto gauss = [](double t, std::vector<double> const& u, std::vector<double>& du) {
            du[0] = -2.0 * t * u[0];
        };
        abm2_adaptive_integrator integ(1);
        std::vector<double> u{1.0};
        double t = 0.0, dt = 1.0e-3, T = 3.0;
        double maxerr = 0.0, dtmin = 1.0e30, dtmax = 0.0, dt_steep = 0.0, dt_tail = 0.0;
        while (t < T - 1.0e-12) {
            double dt_io = std::min(dt, T - t);
            t = integ.step(gauss, u, t, dt_io, 1.0e-8, 0.0); // absolute error control
            dt = dt_io;
            maxerr = std::max(maxerr, std::abs(u[0] - std::exp(-t * t)));
            dtmin = std::min(dtmin, integ.last_dt());
            dtmax = std::max(dtmax, integ.last_dt());
            if (t > 0.6 && t < 0.9 && dt_steep == 0.0) dt_steep = integ.last_dt();
            if (t > 2.0) dt_tail = integ.last_dt();
        }
        fmt::println("  max err = {:.2e}, steps = {} (+{} rej)", maxerr, integ.accepted(),
                     integ.rejected());
        fmt::println("  dt: steep-middle (t~0.75) = {:.2e}, flat-tail (t~2.5) = {:.2e}, "
                     "range [{:.2e}, {:.2e}]",
                     dt_steep, dt_tail, dtmin, dtmax);

        CHECK(maxerr < 1.0e-4);     // matches the closed form (error controlled)
        CHECK(dt_steep < dt_tail);  // smaller step through the steep feature
        CHECK(dtmax / dtmin > 5.0); // the step genuinely varies over the run
        fmt::println("");
    }

} // TEST_SUITE("integrators: ABM2 adaptive (variable dt)")
