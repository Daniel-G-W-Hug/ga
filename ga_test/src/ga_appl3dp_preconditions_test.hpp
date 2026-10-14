// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga3dp mechanics: the inputs each entry point refuses, with the message pinned so that
// an unrelated error cannot pass for the guard under test
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP mechanics preconditions")
{
    TEST_CASE("pga3dp: a frame tree refuses malformed frames, names and indices")
    {
        CHECK_THROWS_WITH_AS(static_frame3dp("A", vec3dp(1.0, 2.0, 3.0, 2.0)),
                             doctest::Contains("Unitized projective point expected"),
                             std::runtime_error);

        static_system3dp sys;
        CHECK_THROWS_WITH_AS((void)sys.get_pos_trafo(0, 0),
                             doctest::Contains("Index range"), std::runtime_error);

        sys.add_frame(static_frame3dp("W"));
        CHECK_THROWS_WITH_AS(sys.add_frame(static_frame3dp("B"), 1),
                             doctest::Contains("parent_idx must refer"),
                             std::runtime_error);
        // names are unique: they are how a frame is found again
        CHECK_THROWS_WITH_AS(sys.add_frame(static_frame3dp("W")),
                             doctest::Contains("must be unique"), std::runtime_error);
        sys.add_frame(static_frame3dp("B", vec3dp(0.0, 1.0, 0.0, 1.0)), 0);
        CHECK(sys.index_of("B") == 1);
        CHECK_THROWS_WITH_AS((void)sys.index_of("C"),
                             doctest::Contains("no frame named 'C'"), std::runtime_error);
        CHECK_NOTHROW((void)sys.get_pos_trafo(1, 0));
        CHECK_THROWS_WITH_AS((void)sys.get_pos_trafo(2, 0),
                             doctest::Contains("Index range"), std::runtime_error);
    }

    TEST_CASE("pga3dp: a single point mass at the origin has no inverse inertia")
    {
        auto const I = get_point_inertia(2.0, vec3dp(0.0, 0.0, 0.0, 1.0));
        CHECK(det(I.view()) == 0.0);
        CHECK_THROWS_WITH_AS((void)get_inertia_inverse(I), doctest::Contains("singular"),
                             std::invalid_argument);
    }

    TEST_CASE("pga3dp: get_inertia_inverse refuses a singular map relative to its scale")
    {
        // one point mass OFF the origin: no zero row, det() ~1e-75 instead of 0
        auto const single = get_point_inertia(3.0, vec3dp(1.0, 2.0, 0.5, 1.0));
        CHECK(det(single.view()) != 0.0);
        CHECK_THROWS_WITH_AS((void)get_inertia_inverse(single),
                             doctest::Contains("singular"), std::invalid_argument);
        // two point masses still leave the rotation about their joining line free
        auto pair = get_point_inertia(1.0, vec3dp(1.0, 2.0, 0.5, 1.0));
        pair += get_point_inertia(1.0, vec3dp(-0.5, 0.25, 1.0, 1.0));
        CHECK_THROWS_AS((void)get_inertia_inverse(pair), std::invalid_argument);

        // three non-collinear point masses span every motion: invertible at any mass unit
        for (value_t const m : {1.0e-9, 1.0, 1.0e9}) {
            auto I = get_point_inertia(m, vec3dp(1.0, 2.0, 0.5, 1.0));
            I += get_point_inertia(m, vec3dp(-0.5, 0.25, 1.0, 1.0));
            I += get_point_inertia(m, vec3dp(0.3, -1.0, -0.7, 1.0));
            CHECK_NOTHROW((void)get_inertia_inverse(I));
        }
        // a thin body is not a singular one
        CHECK_NOTHROW((void)make_cuboid_body(2.0, 1.0, 1.0e-4, 1.0e-4));
    }

    TEST_CASE("pga3dp: inverse_dynamics needs one acceleration per dof coordinate")
    {
        dynamic_system3dp s;
        s.add_frame(static_frame3dp("W"));
        s.add_body(static_frame3dp("B"), make_cuboid_body(1.0, 1.0, 1.0, 1.0));
        REQUIRE(s.dof_coords().size() == 6);
        CHECK_NOTHROW((void)s.inverse_dynamics(std::vector<value_t>(6, 0.0)));
        CHECK_THROWS_WITH_AS((void)s.inverse_dynamics(std::vector<value_t>(7, 0.0)),
                             doctest::Contains("one entry per dof coordinate (6), got 7"),
                             std::invalid_argument);
    }
}
