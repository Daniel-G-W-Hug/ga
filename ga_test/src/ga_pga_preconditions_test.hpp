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

    TEST_CASE("pga2dp / pga3dp: an ideal line or plane has no angle")
    {
        // the angle is that of the unitized objects; one at infinity (zero weight) has
        // no attitude to measure, so it is refused rather than reported as pi/2
        bivec2dp const l2 = join(vec2dp(0.0, 0.0, 1.0), vec2dp(1.0, 1.0, 1.0));
        bivec2dp const ideal2(0.0, 0.0, 3.0);
        CHECK_THROWS_WITH_AS((void)angle(ideal2, l2), doctest::Contains("angle: line"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(l2, ideal2), doctest::Contains("angle: line"),
                             std::runtime_error);
        CHECK(angle(l2, 2.0 * join(vec2dp(0.0, 0.0, 1.0), vec2dp(1.0, 0.0, 1.0))) ==
              doctest::Approx(pi / 4));

        bivec3dp const l3 = join(vec3dp(0.0, 0.0, 0.0, 1.0), vec3dp(1.0, 0.0, 0.0, 1.0));
        bivec3dp const ideal_l(0.0, 0.0, 0.0, 1.0, 2.0, 3.0);
        trivec3dp const pl(0.0, 0.0, 1.0, -3.0); // z == 3
        trivec3dp const ideal_p(0.0, 0.0, 0.0, 2.0);
        CHECK_THROWS_WITH_AS((void)angle(ideal_l, l3), doctest::Contains("angle: line"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(l3, ideal_l), doctest::Contains("angle: line"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(ideal_p, l3), doctest::Contains("angle: plane"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(pl, ideal_l), doctest::Contains("angle: line"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(l3, ideal_p), doctest::Contains("angle: plane"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(ideal_l, pl), doctest::Contains("angle: line"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(ideal_p, pl), doctest::Contains("angle: plane"),
                             std::runtime_error);
        CHECK_THROWS_WITH_AS((void)angle(pl, ideal_p), doctest::Contains("angle: plane"),
                             std::runtime_error);

        // finite objects still measure, whatever their scale: the x axis lies in z == 3's
        // attitude (angle 0 to the plane), and two parallel planes are at 0
        CHECK(angle(3.0 * pl, l3) == doctest::Approx(0.0));
        CHECK(angle(l3, 0.5 * pl) == doctest::Approx(0.0));
        CHECK(angle(pl, 2.0 * pl) == doctest::Approx(0.0));
        CHECK(angle(l3, 4.0 * l3) == doctest::Approx(0.0));
    }
}
