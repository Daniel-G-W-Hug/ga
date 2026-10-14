// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga2dp mechanics API: the accessors and switches of the kinematic, dynamic,
// closed-loop and contact layers that the topical tests did not reach -- each against
// its index twin, an exact answer, or its documented default
/////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <vector>

TEST_SUITE("PGA2DP mechanics API")
{
    TEST_CASE("pga2dp: Inertia2dp -= removes what += added")
    {
        auto const I1 = get_point_inertia(2.0, vec2dp(1.0, 0.5, 1.0));
        auto const I2 = get_point_inertia(3.0, vec2dp(-0.5, 2.0, 1.0));
        auto I = I1;
        I += I2;
        I -= I2;
        auto const a = I.view();
        auto const b = I1.view();
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
                CHECK(a[i, j] == doctest::Approx(b[i, j]));
    }

    TEST_CASE("kinematic_system2dp: by-name accessors are the by-index ones")
    {
        kinematic_system2dp ks;
        CHECK_FALSE(ks.is_valid_for_transformations());
        ks.add_frame(static_frame2dp("W"));
        CHECK_FALSE(ks.is_valid_for_transformations());
        ks.add_frame(static_frame2dp("A"), kin_state2dp{.omega = 2.0, .alpha = 0.5});
        ks.add_frame(
            static_frame2dp("B", vec2dp(1.0, 0.0, 1.0)),
            kin_state2dp{.vel = vec2dp(0.0, 1.0, 0.0), .acc = vec2dp(0.5, 0.0, 0.0)});
        CHECK(ks.is_valid_for_transformations());
        size_t const B = ks.index_of("B");
        CHECK(ks.relative_accel_twist("B") == ks.relative_accel_twist(B));
        CHECK(ks.accel_twist_world("B") == ks.accel_twist_world(B));
        CHECK(ks.accel_twist_world("B") != ks.accel_twist_world(ks.index_of("A")));
    }

    TEST_CASE("dynamic_system2dp: driven joints, motor joints, springs, torques")
    {
        dynamic_system2dp sys;
        sys.set_gravity(vec2dp(0.0, -9.81, 0.0));
        sys.add_frame(static_frame2dp("W"));
        sys.add_revolute_body(static_frame2dp("R"), make_plate_body(1.0, 0.1, 0.1),
                              vec2dp(0.0, 0.0, 1.0));
        sys.add_body(static_frame2dp("F", vec2dp(2.0, 0.0, 1.0)),
                     make_disc_body(1.0, 0.2), {},
                     0); // a free body: a motor joint with three screws
        size_t const R = sys.index_of("R"), F = sys.index_of("F");
        CHECK(sys.dof_coords().size() == 4);

        CHECK_FALSE(sys.is_driven(R));
        sys.set_driven_rate(R, 1.5);
        CHECK(sys.is_driven(R));
        CHECK(sys.dof_coords().size() == 3);
        sys.clear_driven_joint(R);
        CHECK_FALSE(sys.is_driven(R));
        CHECK(sys.dof_coords().size() == 4);

        // a motor joint's configuration is its motor: setting it poses the body
        auto const M = get_motor(vec2dp(0.0, 0.0, 1.0), deg2rad(90.0));
        sys.set_joint_motor(F, M);
        CHECK(sys.joint_motor(F) == M);
        vec2dp const tip_b(0.5, 0.0, 1.0);
        CHECK(unitize(move2dp(tip_b, sys.get_pos_trafo(F, 0))) ==
              unitize(move2dp(tip_b, M)));

        sys.set_joint_torque(R, [](value_t t) { return 2.0 + t; });
        sys.set_time(0.5);
        auto const tau = sys.joint_torques();
        auto const rc = sys.dof_coords();
        REQUIRE(tau.size() == rc.size());
        for (size_t c = 0; c < rc.size(); ++c)
            CHECK(tau[c] == doctest::Approx(rc[c].frame == R ? 2.5 : 0.0));

        value_t const E0 = sys.potential_energy();
        sys.add_grounded_spring(F, tip_b, vec2dp(0.0, 0.0, 1.0),
                                vec2dp(100.0, 100.0, 0.0), 0.0);
        CHECK(sys.potential_energy() > E0);
        sys.clear_grounded_springs(F);
        CHECK(sys.potential_energy() == doctest::Approx(E0));
    }

    TEST_CASE("closed_loop_system2dp: the constraint rank drops at a singular pose")
    {
        // two unit links from the origin with the tip pinned to the ground: generically
        // the pin's two rows are independent; stretched straight, they are not
        // (Featherstone 8.10: the rank varies with the configuration)
        auto build = [](value_t q1, value_t q2, vec2dp const& anchor) {
            closed_loop_system2dp cl;
            auto const link = make_plate_body(1.0, 1.0, 0.1);
            cl.add_frame(static_frame2dp("W"));
            cl.add_revolute_body(static_frame2dp("L1"), link, vec2dp(0.0, 0.0, 1.0), q1);
            cl.add_revolute_body(static_frame2dp("L2", vec2dp(1.0, 0.0, 1.0)), link,
                                 vec2dp(0.0, 0.0, 1.0), q2);
            cl.add_loop_constraint(loop_constraint2dp{cl.index_of("L2"),
                                                      vec2dp(1.0, 0.0, 1.0), 0, anchor,
                                                      constraint2dp::coincidence});
            return cl;
        };
        // 60 deg and -120 deg relative: the tip is at (1, 0)
        auto generic = build(deg2rad(60.0), deg2rad(-120.0), vec2dp(1.0, 0.0, 1.0));
        CHECK(generic.residual_norm() < 1.0e-12);
        CHECK(generic.constraint_rows() == 2);
        CHECK(generic.constraint_rank() == 2);
        auto straight = build(0.0, 0.0, vec2dp(2.0, 0.0, 1.0));
        CHECK(straight.residual_norm() < 1.0e-12);
        CHECK(straight.constraint_rank() == 1);

        // the read-only view of the tree reports the same state
        closed_loop_system2dp const& cg = generic;
        CHECK(cg.system().joint_phi(cg.index_of("L1")) == doctest::Approx(deg2rad(60.0)));
    }

    TEST_CASE("closed_loop_system2dp: a prismatic body")
    {
        closed_loop_system2dp cl;
        cl.add_frame(static_frame2dp("W"));
        cl.add_prismatic_body(static_frame2dp("slide"), make_plate_body(1.0, 0.1, 0.1),
                              vec2dp(1.0, 0.0, 0.0), 0.75);
        size_t const P = cl.index_of("slide");
        CHECK(cl.system().joint_phi(P) == doctest::Approx(0.75));
        CHECK(unitize(move2dp(O_2dp, cl.system().get_pos_trafo(P, 0))) ==
              vec2dp(0.75, 0.0, 1.0));
    }

    TEST_CASE("ground_contact2dp: defaults and the run-time switches")
    {
        closed_loop_system2dp cl;
        cl.add_frame(static_frame2dp("W"));
        cl.add_body(static_frame2dp("rod", vec2dp(0.0, 0.3, 1.0)),
                    make_plate_body(2.0, 1.0, 0.02));
        ground_contact2dp gc(cl, ground_contact2dp::ground_line(vec2dp(0.0, 0.0, 1.0),
                                                                vec2dp(1.0, 0.0, 1.0)));
        CHECK(gc.simultaneity() == doctest::Approx(1.0e-6));
        gc.set_simultaneity(1.0e-4);
        CHECK(gc.simultaneity() == doctest::Approx(1.0e-4));
        CHECK(gc.release_threshold() == doctest::Approx(0.5));
        gc.set_release_threshold(0.1);
        CHECK(gc.release_threshold() == doctest::Approx(0.1));

        size_t const c = gc.add({cl.index_of("rod"), vec2dp(-0.5, 0.0, 1.0)});
        CHECK(gc.unilateral(c));
        gc.set_unilateral(c, false);
        CHECK_FALSE(gc.unilateral(c));
    }
}
