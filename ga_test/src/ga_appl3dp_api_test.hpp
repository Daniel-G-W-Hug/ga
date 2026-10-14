// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga3dp mechanics API: the accessors and switches of the kinematic, dynamic,
// closed-loop and contact layers that the topical tests did not reach -- each against
// its index twin, an exact answer, or its documented default
/////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <vector>

TEST_SUITE("PGA3DP mechanics API")
{
    TEST_CASE("pga3dp: Inertia3dp -= removes what += added")
    {
        auto const I1 = get_point_inertia(2.0, vec3dp(1.0, 0.5, -1.0, 1.0));
        auto const I2 = get_point_inertia(3.0, vec3dp(-0.5, 2.0, 0.25, 1.0));
        auto I = I1;
        I += I2;
        I -= I2;
        auto const a = I.view();
        auto const b = I1.view();
        for (size_t i = 0; i < 6; ++i)
            for (size_t j = 0; j < 6; ++j)
                CHECK(a[i, j] == doctest::Approx(b[i, j]));
    }

    TEST_CASE("kinematic_system3dp: by-name accessors, and the centripetal acceleration")
    {
        value_t const w = 2.0;
        kinematic_system3dp ks;
        CHECK_FALSE(ks.is_valid_for_transformations());
        ks.add_frame(static_frame3dp("W"));
        CHECK_FALSE(ks.is_valid_for_transformations()); // one frame relates nothing
        // A spins about the z axis through the world origin; B rides on A, offset along
        // x, and moves and accelerates relative to it
        ks.add_frame(static_frame3dp("A"),
                     kin_state3dp{.omega = vec3dp(0.0, 0.0, w, 0.0)});
        ks.add_frame(static_frame3dp("B", vec3dp(1.0, 0.0, 0.0, 1.0)),
                     kin_state3dp{.vel = vec3dp(0.0, 1.0, 0.0, 0.0),
                                  .acc = vec3dp(0.5, 0.0, 0.0, 0.0)});
        CHECK(ks.is_valid_for_transformations());
        CHECK(ks.is_linear_chain());

        size_t const A = ks.index_of("A"), B = ks.index_of("B");
        vec3dp const X(1.5, 0.0, 0.0, 1.0);
        CHECK(ks.twist_world("B") == ks.twist_world(B));
        CHECK(ks.relative_accel_twist("B") == ks.relative_accel_twist(B));
        CHECK(ks.accel_twist_world("B") == ks.accel_twist_world(B));
        CHECK(ks.point_acceleration(X, "B") == ks.point_acceleration(X, B));

        // a point at radius r on a frame spinning at w about z: speed w r tangentially,
        // acceleration w^2 r towards the axis
        CHECK(ks.point_velocity(X, "A") == vec3dp(0.0, w * 1.5, 0.0, 0.0));
        CHECK(ks.point_acceleration(X, "A") == vec3dp(-w * w * 1.5, 0.0, 0.0, 0.0));

        // the by-name setters write what the by-index ones do
        kinematic_system3dp ks2 = ks;
        kin_state3dp const k{.vel = vec3dp(0.0, 0.0, 2.0, 0.0),
                             .omega = vec3dp(1.0, 0.0, 0.0, 0.0)};
        ks.set_state("B", k);
        ks2.set_state(B, k);
        CHECK(ks.relative_twist(B) == ks2.relative_twist(B));
        twist3dp const T(0.5, 0.0, 0.0, 0.0, 1.0, 0.0);
        ks.set_twist("A", T);
        CHECK(ks.relative_twist(A) == T);

        // a second child of the root branches the tree
        ks.add_frame(static_frame3dp("C"), kin_state3dp{}, 0);
        CHECK_FALSE(ks.is_linear_chain());
    }

    TEST_CASE("dynamic_system3dp: driven joints, motor joints, integrator, springs, "
              "torques")
    {
        dynamic_system3dp sys;
        sys.set_gravity(vec3dp(0.0, 0.0, -9.81, 0.0));
        sys.add_frame(static_frame3dp("W"));
        sys.add_revolute_body(static_frame3dp("R"), make_cuboid_body(1.0, 0.1, 0.1, 0.1),
                              O_3dp, vec3dp(0.0, 0.0, 1.0, 0.0));
        sys.add_spherical_body(static_frame3dp("S", vec3dp(1.0, 0.0, 0.0, 1.0)),
                               make_cuboid_body(1.0, 0.2, 0.2, 0.2), O_3dp);
        size_t const R = sys.index_of("R"), S = sys.index_of("S");
        CHECK(sys.dof_coords().size() == 4); // 1 + 3

        // a driven joint leaves the dynamic coordinates, and clearing restores it
        CHECK_FALSE(sys.is_driven(R));
        sys.set_driven_rate(R, 1.5);
        CHECK(sys.is_driven(R));
        CHECK(sys.dof_coords().size() == 3);
        sys.clear_driven_joint(R);
        CHECK_FALSE(sys.is_driven(R));
        CHECK(sys.dof_coords().size() == 4);

        // a motor joint's configuration is its motor: setting it poses the body
        auto const M = get_motor(join(O_3dp, vec3dp(0.0, 0.0, 1.0, 1.0)), deg2rad(90.0));
        sys.set_joint_motor(S, M);
        CHECK(sys.joint_motor(S) == M);
        vec3dp const tip_b(0.5, 0.0, 0.0, 1.0);
        auto const tip_world = move3dp(tip_b, sys.get_pos_trafo(S, R));
        CHECK(unitize(tip_world) == unitize(move3dp(tip_b, M)));

        // the registered torques at the current clock, one per dof coordinate
        sys.set_joint_torque(R, [](value_t t) { return 2.0 + t; });
        sys.set_time(0.5);
        auto const tau = sys.joint_torques();
        auto const rc = sys.dof_coords();
        REQUIRE(tau.size() == rc.size());
        for (size_t c = 0; c < rc.size(); ++c)
            CHECK(tau[c] == doctest::Approx(rc[c].frame == R ? 2.5 : 0.0));

        // the integrator is a switch: RK4 by default, ABM2 on request
        CHECK(sys.get_integrator() == integrator_kind::rk4);
        sys.set_integrator(integrator_kind::abm2);
        CHECK(sys.get_integrator() == integrator_kind::abm2);
        sys.set_integrator(integrator_kind::rk4);

        // a grounded spring adds potential energy; clearing it takes exactly that back
        value_t const E0 = sys.potential_energy();
        sys.add_grounded_spring(S, tip_b, vec3dp(0.0, 0.0, 0.0, 1.0),
                                vec3dp(100.0, 100.0, 100.0, 0.0), 0.0);
        CHECK(sys.potential_energy() > E0);
        sys.clear_grounded_springs(S);
        CHECK(sys.potential_energy() == doctest::Approx(E0));
    }

    TEST_CASE("dynamic_system3dp: RK4 and ABM2 integrate the same pendulum")
    {
        auto run = [](integrator_kind k) {
            dynamic_system3dp sys;
            sys.set_gravity(vec3dp(0.0, 0.0, -9.81, 0.0));
            sys.add_frame(static_frame3dp("W"));
            sys.add_revolute_body(static_frame3dp("P", vec3dp(0.0, 0.0, -1.0, 1.0)),
                                  make_cuboid_body(1.0, 0.05, 0.05, 0.05),
                                  vec3dp(0.0, 0.0, 1.0, 1.0), vec3dp(0.0, 1.0, 0.0, 0.0),
                                  0.3);
            sys.set_integrator(k);
            for (int n = 0; n < 1000; ++n)
                sys.step(1.0e-4);
            return sys.joint_phi(sys.index_of("P"));
        };
        // a 2nd-order and a 4th-order scheme at dt = 1e-4 over 0.1 s: agreement far
        // below the swing, and a swing that actually happened
        value_t const q4 = run(integrator_kind::rk4);
        value_t const q2 = run(integrator_kind::abm2);
        CHECK(std::abs(q4 - 0.3) > 1.0e-3);
        CHECK(std::abs(q4 - q2) < 1.0e-6);
    }

    TEST_CASE("closed_loop_system3dp: the planar four-bar in space")
    {
        // the 2D chapter's four-bar (a = 2, b = 3, c = sqrt 5, d = 4, crank at 90 deg),
        // built in the xy plane of space with every hinge along z: the same exact
        // velocities and accelerations as in the plane
        value_t const a = 2.0, b = 3.0, c = std::sqrt(5.0), d = 4.0;
        auto const link = make_cuboid_body(1.0, 1.0, 1.0, 1.0);
        vec3dp const z_axis(0.0, 0.0, 1.0, 0.0);
        closed_loop_system3dp cl;
        cl.add_frame(static_frame3dp("W"));
        cl.add_revolute_body(static_frame3dp("CR"), link, O_3dp, z_axis, 0.0, 0.0,
                             cl.index_of("W"));
        cl.add_revolute_body(static_frame3dp("CO", vec3dp(a, 0.0, 0.0, 1.0)), link, O_3dp,
                             z_axis, -1.4, 0.0, cl.index_of("CR"));
        cl.add_revolute_body(static_frame3dp("RO", vec3dp(d, 0.0, 0.0, 1.0)), link, O_3dp,
                             z_axis, 2.2, 0.0, cl.index_of("W"));
        size_t const CR = cl.index_of("CR"), CO = cl.index_of("CO"),
                     RO = cl.index_of("RO");
        cl.add_loop_constraint(loop_constraint3dp{CO, vec3dp(b, 0.0, 0.0, 1.0), RO,
                                                  vec3dp(c, 0.0, 0.0, 1.0),
                                                  constraint3dp::coincidence});
        cl.set_joint(CR, pi / 2.0);
        CHECK(cl.system().joint_phi(CR) == doctest::Approx(pi / 2.0));
        CHECK(cl.assemble({CR}) < 1.0e-12);

        // a coincidence in space is 3 rows, but a planar loop moves in the plane only:
        // the out-of-plane row is dependent
        CHECK(cl.constraint_rows() == 3);
        CHECK(cl.constraint_rank() == 2);

        cl.set_joint_rate(CR, 1.0);
        cl.solve_velocities({CR});
        CHECK(cl.joint_rate(CO) == doctest::Approx(-4.0 / 3.0).epsilon(1e-9));
        CHECK(cl.joint_rate(RO) == doctest::Approx(1.0).epsilon(1e-9));
        auto const adep = cl.solve_accelerations({CR}, {0.0});
        REQUIRE(adep.size() == 2);
        CHECK(adep[0] == doctest::Approx(-2.0 / 9.0).epsilon(1e-9));
        CHECK(adep[1] == doctest::Approx(2.0 / 3.0).epsilon(1e-9));

        // the read-only view of the tree reports the same state
        closed_loop_system3dp const& ccl = cl;
        CHECK(ccl.system().joint_phi(CO) == cl.system().joint_phi(CO));
    }

    TEST_CASE("closed_loop_system3dp: prismatic and spherical bodies")
    {
        closed_loop_system3dp cl;
        cl.add_frame(static_frame3dp("W"));
        cl.add_prismatic_body(static_frame3dp("slide"),
                              make_cuboid_body(1.0, 0.1, 0.1, 0.1),
                              vec3dp(1.0, 0.0, 0.0, 0.0), 0.75);
        cl.add_spherical_body(static_frame3dp("ball", vec3dp(0.0, 0.0, 1.0, 1.0)),
                              make_cuboid_body(1.0, 0.1, 0.1, 0.1), O_3dp);
        size_t const P = cl.index_of("slide");
        CHECK(cl.system().joint_phi(P) == doctest::Approx(0.75));
        CHECK(cl.system().dof_coords().size() == 4); // 1 + 3
        // the slider's origin sits 0.75 along x
        CHECK(unitize(move3dp(O_3dp, cl.system().get_pos_trafo(P, 0))) ==
              vec3dp(0.75, 0.0, 0.0, 1.0));
    }

    TEST_CASE("ground_contact3dp: accessors, defaults, release and the separation test")
    {
        value_t const g = 9.81;
        auto const floor = ground_contact3dp::ground_plane(vec3dp(0.0, 0.0, 0.0, 1.0),
                                                           vec3dp(1.0, 0.0, 0.0, 1.0),
                                                           vec3dp(0.0, 1.0, 0.0, 1.0));
        closed_loop_system3dp cl;
        cl.system().set_gravity(vec3dp(0.0, 0.0, -g, 0.0));
        cl.add_frame(static_frame3dp("W"));
        cl.add_body(static_frame3dp("rod", vec3dp(0.0, 0.0, 0.3, 1.0),
                                    vec3dp(0.0, -0.1, 0.0, 0.0)),
                    make_cuboid_body(2.0, 1.0, 0.02, 0.02));
        ground_contact3dp gc(cl, floor);
        CHECK(gc.ground() == floor);

        // the documented defaults
        CHECK(gc.release_threshold() == doctest::Approx(0.5));
        CHECK(gc.simultaneity() == doctest::Approx(1.0e-6));
        gc.set_release_threshold(0.2);
        CHECK(gc.release_threshold() == doctest::Approx(0.2));

        size_t const c = gc.add({cl.index_of("rod"), vec3dp(-0.5, 0.0, 0.0, 1.0)});
        CHECK(gc.count() == 1);
        CHECK(gc.kind(c) == contact_kind3dp::point);
        CHECK(gc.unilateral(c));
        gc.set_unilateral(c, false);
        CHECK_FALSE(gc.unilateral(c));
        gc.set_unilateral(c, true);

        // fall onto the low end: it lands and pins
        for (int i = 0; i < 400 && gc.events().empty(); ++i)
            gc.step(1.0e-3);
        REQUIRE(gc.active(c));
        CHECK(cl.constraint_rows() == 3);

        // the lift-off rule (a pull beyond the threshold AND a trial step that
        // separates): under gravity the pinned end presses into the floor and stays
        // engaged ...
        size_t const n_events = gc.events().size();
        for (int i = 0; i < 50; ++i)
            gc.step(1.0e-3);
        CHECK(gc.active(c));
        CHECK(gc.events().size() == n_events);
        // ... with gravity reversed it is pulled off, and the step releases it
        cl.system().set_gravity(vec3dp(0.0, 0.0, g, 0.0));
        for (int i = 0; i < 50 && gc.active(c); ++i)
            gc.step(1.0e-3);
        CHECK_FALSE(gc.active(c));
        REQUIRE(gc.events().size() > n_events);
        CHECK_FALSE(gc.events().back().touchdown);
        cl.system().set_gravity(vec3dp(0.0, 0.0, -g, 0.0));

        // released by hand: land again first, then release
        for (int i = 0; i < 2000 && !gc.active(c); ++i)
            gc.step(1.0e-3);
        REQUIRE(gc.active(c));
        // released by hand: the constraint goes, the contact stays registered
        gc.release(c);
        CHECK_FALSE(gc.active(c));
        CHECK(cl.constraint_rows() == 0);
        CHECK(gc.count() == 1);
    }
}
