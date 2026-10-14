// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// PGA operations that need a finite object: a motor about a fixed point or a fixed line
// is refused for an IDEAL one, which has no position to turn about
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA preconditions")
{
    TEST_CASE("pga2dp: get_motor needs a finite fixed point")
    {
        vec2dp const ideal(1.0, 2.0, 0.0); // a direction, not a point
        CHECK_THROWS_WITH_AS((void)get_motor(ideal, pi / 2), doctest::Contains("ideal"),
                             std::invalid_argument);

        // a non-unitized finite point is accepted and acts as its unitized self
        vec2dp const P(2.0, 4.0, 2.0);
        auto const M = get_motor(P, pi / 2);
        CHECK(move2dp(unitize(P), M) == unitize(P));
    }

    TEST_CASE("pga3dp: get_motor needs a finite rotation axis")
    {
        // a line at infinity (zero weight) has no axis to turn about
        bivec3dp const ideal(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        CHECK_THROWS_WITH_AS((void)get_motor(ideal, pi / 2), doctest::Contains("ideal"),
                             std::invalid_argument);

        // the axis itself stays fixed under the motor
        bivec3dp const L = join(vec3dp(1.0, 0.0, 0.0, 1.0), vec3dp(1.0, 0.0, 1.0, 1.0));
        auto const M = get_motor(L, pi / 3);
        CHECK(move3dp(vec3dp(1.0, 0.0, 0.5, 1.0), M) == vec3dp(1.0, 0.0, 0.5, 1.0));
    }
}
