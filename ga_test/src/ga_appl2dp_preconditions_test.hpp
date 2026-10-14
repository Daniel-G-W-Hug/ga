// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga2dp mechanics and geodesics: the inputs each entry point refuses, with the message
// pinned so that an unrelated error cannot pass for the guard under test
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA2DP mechanics preconditions")
{
    TEST_CASE("pga2dp: a frame tree refuses malformed frames and indices")
    {
        // a frame's origin is a unitized point
        CHECK_THROWS_WITH_AS(static_frame2dp("A", vec2dp(1.0, 2.0, 2.0)),
                             doctest::Contains("Unitized projective point expected"),
                             std::runtime_error);

        static_system2dp sys;
        // an empty system has no transformation, not even the identity
        CHECK_THROWS_WITH_AS((void)sys.get_pos_trafo(0, 0),
                             doctest::Contains("Index range"), std::runtime_error);

        sys.add_frame(static_frame2dp("W"));
        // a parent must already exist
        CHECK_THROWS_WITH_AS(sys.add_frame(static_frame2dp("B"), 3),
                             doctest::Contains("parent_idx must refer"),
                             std::runtime_error);
        sys.add_frame(static_frame2dp("B", vec2dp(1.0, 0.0, 1.0)), 0);
        CHECK_NOTHROW((void)sys.get_pos_trafo(1, 0));
        CHECK_THROWS_WITH_AS((void)sys.get_pos_trafo(0, 2),
                             doctest::Contains("Index range"), std::runtime_error);
    }

    TEST_CASE("pga2dp: a single point mass at the origin has no inverse inertia")
    {
        // every row of the point-mass map vanishes but the two rotational couplings, so
        // the determinant is exactly zero
        auto const I = get_point_inertia(2.0, vec2dp(0.0, 0.0, 1.0));
        CHECK(det(I.view()) == 0.0);
        CHECK_THROWS_WITH_AS((void)get_inertia_inverse(I), doctest::Contains("singular"),
                             std::invalid_argument);
    }

    TEST_CASE("pga2dp: get_inertia_inverse refuses a singular map relative to its scale")
    {
        // one point mass OFF the origin: no zero row, and det() returns a rounding-level
        // ~1e-20 instead of 0 -- an exact-zero test let this through with entries of 2e20
        auto const single = get_point_inertia(3.0, vec2dp(1.0, 2.0, 1.0));
        CHECK(det(single.view()) != 0.0);
        CHECK_THROWS_WITH_AS((void)get_inertia_inverse(single),
                             doctest::Contains("singular"), std::invalid_argument);
        // ... at any mass: the test is relative
        CHECK_THROWS_AS(
            (void)get_inertia_inverse(get_point_inertia(1.0e6, vec2dp(1.0, 2.0, 1.0))),
            std::invalid_argument);

        // two point masses span the plane's motions: invertible, at any mass unit
        for (value_t const m : {1.0e-9, 1.0, 1.0e9}) {
            auto I = get_point_inertia(m, vec2dp(1.0, 2.0, 1.0));
            I += get_point_inertia(m, vec2dp(-0.5, 0.25, 1.0));
            CHECK_NOTHROW((void)get_inertia_inverse(I));
        }
        // a thin body is not a singular one (a diagonal inertia passes however unequal)
        CHECK_NOTHROW((void)make_plate_body(2.0, 1.0, 1.0e-4));
    }

    TEST_CASE("pga2dp: inverse_dynamics needs one acceleration per dof coordinate")
    {
        dynamic_system2dp s;
        s.add_frame(static_frame2dp("W"));
        s.add_body(static_frame2dp("B"), make_disc_body(1.0, 0.5)); // a free body
        REQUIRE(s.dof_coords().size() == 3);
        CHECK_NOTHROW((void)s.inverse_dynamics(std::vector<value_t>(3, 0.0)));
        CHECK_THROWS_WITH_AS((void)s.inverse_dynamics(std::vector<value_t>(2, 0.0)),
                             doctest::Contains("one entry per dof coordinate (3), got 2"),
                             std::invalid_argument);
    }

    TEST_CASE("geodesics: dms2deg and deg2dms refuse what they cannot read or write")
    {
        // a lone decimal point is not a number
        CHECK_THROWS_WITH_AS((void)dms2deg(".N", geo_angle::latitude),
                             doctest::Contains("malformed number"),
                             std::invalid_argument);
        // nothing may follow the hemisphere letter
        CHECK_THROWS_WITH_AS((void)dms2deg("52.5N x", geo_angle::latitude),
                             doctest::Contains("trailing characters"),
                             std::invalid_argument);
        // a longitude beyond a full turn
        CHECK_THROWS_WITH_AS((void)dms2deg("400E", geo_angle::longitude),
                             doctest::Contains("longitude outside"),
                             std::invalid_argument);
        CHECK(dms2deg("359.5E", geo_angle::longitude) == 359.5);

        // the seconds' resolution is 0..9 decimals
        CHECK_THROWS_AS((void)deg2dms(52.5, geo_angle::latitude, 10),
                        std::invalid_argument);
        CHECK_THROWS_AS((void)deg2dms(52.5, geo_angle::latitude, -1),
                        std::invalid_argument);
        CHECK(deg2dms(52.5, geo_angle::latitude, 9) == "52°30'00.000000000\"N");
    }
}
