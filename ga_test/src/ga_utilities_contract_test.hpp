// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// ga_usr_utilities: the unit conversions, the two step-size bounds, and the integrators'
// bookkeeping (reset, dim, the adaptive integrate loop) -- each against its stated
// formula or an exact solution
/////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <vector>

TEST_SUITE("utilities: contracts")
{
    using hd::ga::abm2_adaptive_integrator;
    using hd::ga::abm2_integrator;
    using hd::ga::critical_damping;
    using hd::ga::Hz2radps;
    using hd::ga::max_explicit_stiffness;
    using hd::ga::radps2Hz;
    using hd::ga::radps2rpm;
    using hd::ga::rpm2radps;
    using hd::ga::stopping_accel_bound;
    using hd::ga::value_t;

    TEST_CASE("unit conversions: rpm, Hz and rad/s, both ways")
    {
        value_t const tau = 2.0 * std::numbers::pi;
        CHECK(rpm2radps(60.0) == doctest::Approx(tau));
        CHECK(radps2rpm(tau) == doctest::Approx(60.0));
        CHECK(Hz2radps(1.0) == doctest::Approx(tau));
        CHECK(radps2Hz(tau) == doctest::Approx(1.0));
        CHECK(radps2rpm(rpm2radps(1234.5)) == doctest::Approx(1234.5));
        CHECK(radps2Hz(Hz2radps(50.0)) == doctest::Approx(50.0));
        CHECK(rpm2radps(radps2rpm(3.0)) == doctest::Approx(Hz2radps(radps2Hz(3.0))));
    }

    TEST_CASE("max_explicit_stiffness: the critical damper then spends exactly the ratio")
    {
        // c = 2 sqrt(k I) with k the bound gives c dt / I == ratio, whatever I and dt
        for (value_t const inertia : {1.0e-3, 0.5, 40.0}) {
            for (value_t const dt : {1.0e-4, 1.0e-3}) {
                value_t const k = max_explicit_stiffness(inertia, dt, 0.8);
                CHECK(critical_damping(k, inertia) * dt / inertia ==
                      doctest::Approx(0.8));
            }
        }
        // the bound falls with the square of the step and in proportion to the inertia
        CHECK(max_explicit_stiffness(2.0, 1.0e-3) ==
              doctest::Approx(4.0 * max_explicit_stiffness(2.0, 2.0e-3)));
        CHECK(max_explicit_stiffness(2.0, 1.0e-3) ==
              doctest::Approx(2.0 * max_explicit_stiffness(1.0, 1.0e-3)));
        // a non-positive inertia admits nothing
        CHECK(max_explicit_stiffness(-1.0, 1.0e-3) == 0.0);
    }

    TEST_CASE("stopping_accel_bound: the braking (viability) bound towards a wall")
    {
        value_t const a_max = 2.0, dt = 0.01;
        // a coordinate moving at exactly the rate the distance can still absorb may not
        // accelerate further towards the wall
        value_t const d = 4.0;
        value_t const v_allow = std::sqrt(2.0 * a_max * d);
        CHECK(stopping_accel_bound(d, v_allow, a_max, dt) == doctest::Approx(0.0));
        // slower than that, it may: by the rate deficit over one step
        CHECK(stopping_accel_bound(d, v_allow - 0.5, a_max, dt) ==
              doctest::Approx(0.5 / dt));
        // at or past the wall only braking is allowed: the bound removes the whole rate
        CHECK(stopping_accel_bound(0.0, 1.5, a_max, dt) == doctest::Approx(-1.5 / dt));
        CHECK(stopping_accel_bound(-0.2, 1.5, a_max, dt) == doctest::Approx(-1.5 / dt));
        // moving away from the wall, the bound only grows
        CHECK(stopping_accel_bound(d, -1.0, a_max, dt) >
              stopping_accel_bound(d, 1.0, a_max, dt));
    }

    TEST_CASE("abm2 integrators: dim, reset, and the adaptive integrate loop")
    {
        // u' = -u from u(0) = 1: u(t) = exp(-t)
        auto const decay = [](double, std::vector<double> const& u,
                              std::vector<double>& du) { du[0] = -u[0]; };

        abm2_integrator ab(1);
        CHECK(ab.dim() == 1);
        std::vector<double> u1{1.0};
        double const t1 = ab.step(decay, u1, 0.0, 0.1); // self-start (one RK4 step)
        ab.step(decay, u1, t1, 0.1);
        // reset restarts the multistep history: the next step is a self-start again,
        // identical to a fresh integrator's first step
        ab.reset();
        std::vector<double> ua{1.0}, ub{1.0};
        ab.step(decay, ua, 0.0, 0.1);
        abm2_integrator fresh(1);
        fresh.step(decay, ub, 0.0, 0.1);
        CHECK(ua[0] == ub[0]);

        abm2_adaptive_integrator ad(1);
        std::vector<double> u{1.0};
        double const t_end = ad.integrate(decay, u, 0.0, 2.0, 1.0e-3, 1.0e-9, 1.0e-9);
        CHECK(t_end == 2.0); // the last step is clamped onto t_end exactly
        CHECK(u[0] == doctest::Approx(std::exp(-2.0)).epsilon(1.0e-6));
        CHECK(ad.accepted() > 0);
        CHECK(ad.last_dt() > 0.0);
        ad.reset();
        CHECK(ad.accepted() == 0);
        CHECK(ad.rejected() == 0);
    }

    TEST_CASE("the extended-test gate reaches this target")
    {
        // the library's safety gates are delivered by linking the ga target; this pins
        // that the build's answer and the macro the target was compiled with agree
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        CHECK(hd::ga::detail::extended_testing_enabled());
#else
        CHECK_FALSE(hd::ga::detail::extended_testing_enabled());
#endif
    }

    TEST_CASE("safe_epsilon: an integral type contributes no epsilon of its own")
    {
        using hd::ga::detail::safe_epsilon;
        CHECK(safe_epsilon<int, double>() == safe_epsilon<double>());
        CHECK(safe_epsilon<double, int>() == safe_epsilon<double>());
        CHECK(safe_epsilon<float, double>() == safe_epsilon<float, float>());
    }
}
