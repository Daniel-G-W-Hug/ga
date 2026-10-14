// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// ga_usr_utilities: the RK4 stage helpers and discrete_range -- their stage times and
// accessors, and the arguments they refuse
/////////////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <stdexcept>
#include <vector>

TEST_SUITE("utilities: preconditions")
{
    using hd::ga::discrete_range;
    using hd::ga::rk4_get_time;
    using hd::ga::rk4_step;
    using hd::ga::value_t;

    TEST_CASE("rk4_get_time: the four stage times, and nothing beyond")
    {
        value_t const t = 2.0, dt = 0.5;
        CHECK(rk4_get_time(t, dt, 0) == t);
        CHECK(rk4_get_time(t, dt, 1) == t + dt / 2);
        CHECK(rk4_get_time(t, dt, 2) == t + dt / 2);
        CHECK(rk4_get_time(t, dt, 3) == t + dt);
        CHECK_THROWS_WITH_AS((void)rk4_get_time(t, dt, 4), doctest::Contains("rk_step"),
                             std::invalid_argument);
    }

    TEST_CASE("rk4_step: stages are numbered 1..4")
    {
        std::vector<value_t> const u{1.0};
        std::array<std::vector<value_t>, 2> const uh{u, u};
        std::vector<value_t> const rhs{0.0};
        CHECK_THROWS_AS((void)rk4_step(u, uh, rhs, 0.1, 0), std::invalid_argument);
        CHECK_THROWS_AS((void)rk4_step(u, uh, rhs, 0.1, 5), std::invalid_argument);
        CHECK_NOTHROW((void)rk4_step(u, uh, rhs, 0.1, 1));

        // the mdspan form refuses the same stage numbers
        std::vector<value_t> um{1.0}, uhm{1.0, 1.0};
        std::mdspan<value_t, std::dextents<size_t, 1>> us(um.data(), 1);
        std::mdspan<value_t, std::dextents<size_t, 2>> uhs(uhm.data(), 2, 1);
        std::mdspan<value_t const, std::dextents<size_t, 1>> rs(rhs.data(), 1);
        CHECK_THROWS_AS(rk4_step(us, uhs, rs, value_t(0.1), 0), std::invalid_argument);
        CHECK_THROWS_AS(rk4_step(us, uhs, rs, value_t(0.1), 5), std::invalid_argument);
    }

    TEST_CASE("discrete_range: accessors, and an inverted range is refused")
    {
        discrete_range const r(-1.0, 3.0, 8);
        CHECK(r.min() == -1.0);
        CHECK(r.max() == 3.0);
        CHECK(r.steps() == 8);
        CHECK(r.delta() == 0.5);
        CHECK_THROWS_WITH_AS(discrete_range(1.0, 0.0, 4), doctest::Contains("min <= max"),
                             std::invalid_argument);
    }
}
