// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga3dp mechanics API: the accessors and switches of the kinematic, dynamic,
// closed-loop and contact layers that the topical tests did not reach -- each against
// its index twin, an exact answer, or its documented default
/////////////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <functional>
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

    TEST_CASE("pga3dp: a joint stop bounces the joint back elastically, and never pulls")
    {
        // the lower stop: a hinge spinning at -omega0 runs into lo; without damping it
        // penetrates by omega0 sqrt(I / k) and leaves at +omega0
        value_t const lo = -0.3, k = 400.0, omega0 = 1.0, dt = 1.0e-4;
        auto build = [&](value_t c_stop) {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{}); // the default gravity is on
            s.add_frame(static_frame3dp("W"));
            s.add_revolute_body(
                static_frame3dp("P"), make_cuboid_body(1.0, 1.0, 0.1, 0.1),
                vec3dp(-0.5, 0.0, 0.0, 1.0), vec3dp(0.0, 0.0, 1.0, 0.0), 0.0, -omega0);
            s.set_joint_range(s.index_of("P"), joint_range3dp{lo, -lo, k, c_stop});
            return s;
        };
        auto bounce = [&](dynamic_system3dp& s, value_t& min_q) {
            size_t const P = s.index_of("P");
            min_q = s.joint_phi(P);
            for (int n = 0; n < 40000 && !(s.joint_phi(P) > lo && s.joint_omega(P) > 0.0);
                 ++n) {
                s.step(dt);
                min_q = std::min(min_q, s.joint_phi(P));
            }
            return s.joint_omega(P);
        };
        auto elastic = build(0.0);
        value_t const I = elastic.mass_matrix()[0];
        value_t min_q = 0.0;
        CHECK(bounce(elastic, min_q) == doctest::Approx(omega0).epsilon(1.0e-3));
        CHECK(lo - min_q == doctest::Approx(omega0 * std::sqrt(I / k)).epsilon(1.0e-2));

        auto damped = build(3.0 * 2.0 * std::sqrt(k * I));
        value_t min_qd = 0.0;
        value_t const omega_d = bounce(damped, min_qd);
        CHECK(damped.joint_phi(damped.index_of("P")) > lo); // it leaves: never pulled
        CHECK(omega_d > 0.0);
        CHECK(omega_d < omega0);
    }

    TEST_CASE("pga3dp: a disc's inertia about an offset pivot line follows Steiner")
    {
        // a pivot line parallel to e3 at distance d along e1: the moments about e2 and e3
        // gain m d^2, the one about e1 nothing -- the parallel-axis theorem
        value_t const m = 2.0, r = 0.5, t = 0.1, d = 0.75;
        auto const I0 = get_disc_inertia(m, r, t);
        auto const Id = get_disc_inertia(
            m, r, t, join(vec3dp(d, 0.0, 0.0, 1.0), vec3dp(d, 0.0, 1.0, 1.0)));
        auto const a = I0.view();
        auto const b = Id.view();
        CHECK(b[3, 0] == doctest::Approx(a[3, 0]));
        CHECK(b[4, 1] == doctest::Approx(a[4, 1] + m * d * d));
        CHECK(b[5, 2] == doctest::Approx(a[5, 2] + m * d * d));
        // the linear block is the mass however the pivot lies
        CHECK(b[0, 3] == doctest::Approx(m));
    }

    TEST_CASE("pga3dp: a subclass adds its own force elements through extra_wrenches()")
    {
        // the extension seam: a subclass contributes world wrenches the base knows
        // nothing about. A constant force F through a free body's centre of mass: x = F
        // t^2 / 2m
        struct pushed : dynamic_system3dp {
            bivec3dp W{};

          protected:

            std::vector<std::pair<size_t, bivec3dp>> extra_wrenches() override
            {
                return {{1, W}};
            }
        };
        value_t const m = 2.0, F = 3.0, T = 0.5;
        pushed s;
        s.set_gravity(vec3dp{});
        s.add_frame(static_frame3dp("W"));
        s.add_body(static_frame3dp("B"), make_cuboid_body(m, 0.2, 0.2, 0.2));
        s.W = wdg(vec3dp(0.0, 0.0, 0.0, 1.0), vec3dp(F, 0.0, 0.0, 0.0)); // through the cm
        for (int n = 0; n < 500; ++n)
            s.step(T / 500.0);
        vec3dp const C = unitize(move3dp(O_3dp, s.get_pos_trafo(1, 0)));
        CHECK(C.x == doctest::Approx(F * T * T / (2.0 * m)).epsilon(1.0e-9));
        CHECK(std::abs(C.y) < 1.0e-12);
        CHECK(std::abs(C.z) < 1.0e-12);
    }

    TEST_CASE("pga3dp: frames without a coordinate, and clearing a motor joint's torque")
    {
        dynamic_system3dp s;
        s.add_frame(static_frame3dp("W"));
        s.add_body(static_frame3dp("B"), make_cuboid_body(1.0, 0.2, 0.2, 0.2));
        s.add_frame(static_frame3dp("tip", vec3dp(1.0, 0.0, 0.0, 1.0))); // rigid on B
        CHECK(s.joint_accel(s.index_of("W")) == 0.0);
        CHECK(s.joint_accel(s.index_of("tip")) == 0.0);

        size_t const B = s.index_of("B");
        s.set_joint_torque(B, [](value_t) {
            return std::vector<value_t>{1.0, 0.0, 0.0, 0.0, 0.0, -1.0};
        });
        value_t sum = 0.0;
        for (auto x : s.joint_torques())
            sum += std::abs(x);
        CHECK(sum == doctest::Approx(2.0));
        s.set_joint_torque(B, std::function<std::vector<value_t>(value_t)>{});
        for (auto x : s.joint_torques())
            CHECK(x == 0.0);
    }

    TEST_CASE("dynamic_system3dp: an empty function clears a wrench or a joint torque")
    {
        // a hinge about e3 without gravity: a torque or a wrench speeds it up, and once
        // cleared it spins on at a constant rate
        auto build = [] {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{}); // the default gravity is on
            s.add_frame(static_frame3dp("W"));
            s.add_revolute_body(
                static_frame3dp("P"), make_cuboid_body(1.0, 1.0, 0.1, 0.1),
                vec3dp(-0.5, 0.0, 0.0, 1.0), vec3dp(0.0, 0.0, 1.0, 0.0), 0.0, 0.0);
            return s;
        };
        auto gain = [](dynamic_system3dp& s) {
            size_t const P = s.index_of("P");
            value_t const before = s.joint_omega(P);
            for (int i = 0; i < 10; ++i)
                s.step(1.0e-3);
            return s.joint_omega(P) - before;
        };
        {
            auto s = build();
            size_t const P = s.index_of("P");
            s.set_joint_torque(P, [](value_t) { return 2.0; });
            CHECK(gain(s) > 0.0);
            s.set_joint_torque(P, std::function<value_t(value_t)>{});
            CHECK(std::abs(gain(s)) < 1.0e-12);
        }
        {
            auto s = build();
            size_t const P = s.index_of("P");
            // a force along e2 at the body's far end: a moment about the hinge
            s.set_applied_wrench(P, [](value_t) {
                return wdg(vec3dp(1.0, 0.0, 0.0, 1.0), vec3dp(0.0, 1.0, 0.0, 0.0));
            });
            CHECK(std::abs(gain(s)) > 1.0e-3);
            s.set_applied_wrench(P, std::function<bivec3dp(value_t)>{});
            CHECK(std::abs(gain(s)) < 1.0e-12);
        }
    }

    TEST_CASE("ground_contact3dp: a point contact promoted to a flat foot")
    {
        value_t const g = 9.81;
        auto const floor = ground_contact3dp::ground_plane(vec3dp(0.0, 0.0, 0.0, 1.0),
                                                           vec3dp(1.0, 0.0, 0.0, 1.0),
                                                           vec3dp(0.0, 1.0, 0.0, 1.0));
        auto build = [&](closed_loop_system3dp& cl, value_t th) {
            cl.system().set_gravity(vec3dp(0.0, 0.0, -g, 0.0));
            cl.add_frame(static_frame3dp("W"));
            cl.add_body(static_frame3dp("box", vec3dp(0.0, 0.0, 0.2, 1.0),
                                        vec3dp(0.0, th, 0.0, 0.0)),
                        make_cuboid_body(2.0, 1.0, 0.4, 0.02));
        };
        // the mechanism's whole momentum: zero once the weld has absorbed the motion
        auto at_rest = [](closed_loop_system3dp& cl) {
            bivec3dp const P = cl.system().total_momentum();
            return std::max({std::abs(P.vx), std::abs(P.vy), std::abs(P.vz),
                             std::abs(P.mx), std::abs(P.my), std::abs(P.mz)});
        };

        { // a flat foot placed where it stands: engage() welds it without an impact,
          // and the body stays put
            closed_loop_system3dp cl;
            build(cl, 0.0);
            ground_contact3dp gc(cl, floor);
            ground_contact3dp::spec s{cl.index_of("box"), vec3dp(0.0, 0.0, -0.2, 1.0)};
            s.kind = contact_kind3dp::flat;
            s.half_length = 0.5;
            s.half_width = 0.2;
            size_t const c = gc.add(s);
            gc.engage(c);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 6);
            REQUIRE(gc.events().size() == 1);
            CHECK(gc.events()[0].impulse.empty());
            for (int i = 0; i < 200; ++i)
                gc.step(1.0e-3);
            CHECK(gc.active(c));
            CHECK(at_rest(cl) < 1.0e-9);
        }
        { // promoted in the air: the weld is made on demand and engages at touchdown
            closed_loop_system3dp cl;
            build(cl, 0.0);
            ground_contact3dp gc(cl, floor);
            size_t const c = gc.add({cl.index_of("box"), vec3dp(0.0, 0.0, 0.0, 1.0)});
            gc.set_kind(c, contact_kind3dp::flat);
            CHECK(gc.kind(c) == contact_kind3dp::flat);
            for (int i = 0; i < 1000 && !gc.active(c); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 6);
            CHECK(at_rest(cl) < 1.0e-9);
        }
        { // promoted after landing as a pin: the switch stops all motion; back to a point
          // leaves only the pin's rows
            closed_loop_system3dp cl;
            build(cl, -0.1);
            ground_contact3dp gc(cl, floor);
            size_t const c = gc.add({cl.index_of("box"), vec3dp(-0.5, 0.0, 0.0, 1.0)});
            for (int i = 0; i < 1000 && !gc.active(c); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 3);
            gc.set_kind(c, contact_kind3dp::flat);
            CHECK(cl.constraint_rows() == 6);
            CHECK(at_rest(cl) < 1.0e-9);
            gc.set_kind(c, contact_kind3dp::point);
            CHECK(cl.constraint_rows() == 3);
        }
    }
}
