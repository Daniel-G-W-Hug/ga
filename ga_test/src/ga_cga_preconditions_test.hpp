// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// CGA radius accessors: an object without a round weight (a flat, or a degenerate meet)
// has no radius. radius() and radius_sq() each refuse it with their own message rather
// than divide by zero; the message is pinned because this build's division check would
// throw anyway, so a bare CHECK_THROWS could not see the guard.
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("CGA preconditions")
{
    TEST_CASE("cga2dc: radius needs a round weight (vector and trivector)")
    {
        // a vector whose round-weight component (e3) is zero
        vec2dc const v(1.0, 0.0, 0.0, 2.0);
        CHECK(round_weight_nrm_sq(v) == 0.0);
        CHECK_THROWS_WITH_AS((void)radius_sq(v), doctest::Contains("no round weight"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)radius(v), doctest::Contains("no round weight"),
                             std::runtime_error);

        // a trivector whose round-weight component is zero
        trivec2dc const t(1.0, 2.0, 0.5, 0.0);
        CHECK(round_weight_nrm_sq(t) == 0.0);
        CHECK_THROWS_WITH_AS((void)radius_sq(t), doctest::Contains("no round weight"),
                             std::runtime_error);
    }

    TEST_CASE("cga3dc: radius needs a round weight (vector and bivector)")
    {
        // a vector whose round-weight component (w) is zero
        vec3dc const v(0.0, 0.0, 1.0, 0.0, 3.0);
        CHECK(round_weight_nrm_sq(v) == 0.0);
        CHECK_THROWS_WITH_AS((void)radius_sq(v), doctest::Contains("no round weight"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)radius(v), doctest::Contains("no round weight"),
                             std::runtime_error);

        // a bivector whose round-weight part (vx, vy, vz) is zero
        bivec3dc const B(0.0, 0.0, 0.0, 1.0, 2.0, 0.5, 0.25, 1.0, 0.0, 3.0);
        CHECK(round_weight_nrm_sq(B) == 0.0);
        CHECK_THROWS_WITH_AS((void)radius_sq(B), doctest::Contains("no round weight"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)radius(B), doctest::Contains("no round weight"),
                             std::runtime_error);
    }
}
