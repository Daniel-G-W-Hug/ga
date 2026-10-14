// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga2dp mechanics API: the accessors and switches of the kinematic, dynamic,
// closed-loop and contact layers that the topical tests did not reach -- each against
// its index twin, an exact answer, or its documented default
/////////////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <functional>
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

    TEST_CASE("pga2dp: a joint stop bounces the joint back elastically, and never pulls")
    {
        // no gravity: a plate spinning at omega0 about its end reaches the top of its
        // range. Without damping the stop is a spring: the joint penetrates by
        // omega0 sqrt(I / k) and comes back out at -omega0 (energy conserved).
        value_t const hi = 0.3, k = 400.0, omega0 = 1.0, dt = 1.0e-4;
        auto build = [&](value_t c_stop) {
            dynamic_system2dp s;
            s.set_gravity(vec2dp(0.0, 0.0, 0.0)); // the default gravity is on
            s.add_frame(static_frame2dp("W"));
            s.add_revolute_body(static_frame2dp("P"), make_plate_body(1.0, 1.0, 0.1),
                                vec2dp(-0.5, 0.0, 1.0), 0.0, omega0);
            s.set_joint_range(s.index_of("P"), joint_range2dp{-hi, hi, k, c_stop});
            return s;
        };
        auto bounce = [&](dynamic_system2dp& s, value_t& max_q) {
            size_t const P = s.index_of("P");
            max_q = s.joint_phi(P);
            for (int n = 0; n < 40000 && !(s.joint_phi(P) < hi && s.joint_omega(P) < 0.0);
                 ++n) {
                s.step(dt);
                max_q = std::max(max_q, s.joint_phi(P));
            }
            return s.joint_omega(P);
        };

        auto elastic = build(0.0);
        value_t const I = elastic.mass_matrix()[0];
        value_t max_q = 0.0;
        value_t const omega_out = bounce(elastic, max_q);
        CHECK(omega_out == doctest::Approx(-omega0).epsilon(1.0e-3));
        CHECK(max_q - hi == doctest::Approx(omega0 * std::sqrt(I / k)).epsilon(1.0e-2));

        // heavily over-damped (c = 3x critical): a damper that could PULL would hold the
        // joint at the boundary on its way out; the clamped stop only pushes, so the
        // joint still leaves -- slower, with energy lost
        auto damped = build(3.0 * 2.0 * std::sqrt(k * I));
        value_t max_qd = 0.0;
        value_t const omega_d = bounce(damped, max_qd);
        CHECK(damped.joint_phi(damped.index_of("P")) < hi);
        CHECK(omega_d < 0.0);
        CHECK(std::abs(omega_d) < omega0);
        CHECK(max_qd < max_q);
    }

    TEST_CASE("pga2dp: an empty function clears a joint torque or an applied wrench")
    {
        dynamic_system2dp s;
        s.set_gravity(vec2dp(0.0, -9.81, 0.0));
        s.add_frame(static_frame2dp("W"));
        s.add_revolute_body(static_frame2dp("R"), make_plate_body(1.0, 1.0, 0.1),
                            vec2dp(-0.5, 0.0, 1.0));
        s.add_body(static_frame2dp("F", vec2dp(3.0, 0.0, 1.0)), make_disc_body(1.0, 0.2),
                   {}, 0);
        size_t const R = s.index_of("R"), F = s.index_of("F");
        auto const sum = [](std::vector<value_t> const& v) {
            value_t a = 0.0;
            for (auto x : v)
                a += std::abs(x);
            return a;
        };

        s.set_joint_torque(R, [](value_t) { return 2.0; });
        s.set_joint_torque(F,
                           [](value_t) { return std::vector<value_t>{1.0, -1.0, 0.5}; });
        CHECK(sum(s.joint_torques()) == doctest::Approx(4.5));
        s.set_joint_torque(R, std::function<value_t(value_t)>{});
        s.set_joint_torque(F, std::function<std::vector<value_t>(value_t)>{});
        CHECK(sum(s.joint_torques()) == 0.0);

        auto const rhs0 = s.mass_bias().second;
        s.set_applied_wrench(R, [](value_t) { return bivec2dp(0.0, 0.0, 5.0); });
        CHECK(s.mass_bias().second != rhs0);
        s.set_applied_wrench(R, std::function<bivec2dp(value_t)>{});
        CHECK(s.mass_bias().second == rhs0);
    }

    TEST_CASE("pga2dp: frames without a coordinate, and a linear chain")
    {
        dynamic_system2dp s;
        s.add_frame(static_frame2dp("W"));
        s.add_revolute_body(static_frame2dp("A"), make_plate_body(1.0, 1.0, 0.1),
                            vec2dp(0.0, 0.0, 1.0));
        s.add_frame(static_frame2dp("tip", vec2dp(1.0, 0.0, 1.0))); // a rigid frame
        CHECK(s.is_linear_chain());
        // a frame that carries no coordinate (the root, a rigid frame) has no
        // acceleration
        CHECK(s.joint_accel(s.index_of("W")) == 0.0);
        CHECK(s.joint_accel(s.index_of("tip")) == 0.0);
    }

    TEST_CASE("ground_contact2dp: a point contact promoted to a flat foot")
    {
        value_t const g = 9.81;
        auto const floor =
            ground_contact2dp::ground_line(vec2dp(0.0, 0.0, 1.0), vec2dp(1.0, 0.0, 1.0));
        auto build = [&](closed_loop_system2dp& cl, value_t th) {
            cl.system().set_gravity(vec2dp(0.0, -g, 0.0));
            cl.add_frame(static_frame2dp("W"));
            cl.add_body(static_frame2dp("rod", vec2dp(0.0, 0.2, 1.0), th),
                        make_plate_body(2.0, 1.0, 0.02));
        };
        // the mechanism's whole momentum: zero once the weld has absorbed the motion
        auto at_rest = [](closed_loop_system2dp& cl) {
            bivec2dp const P = cl.system().total_momentum();
            return std::max({std::abs(P.x), std::abs(P.y), std::abs(P.z)});
        };

        { // a flat foot placed where it stands: engage() welds it without an impact,
          // and the body stays put
            closed_loop_system2dp cl;
            build(cl, 0.0);
            ground_contact2dp gc(cl, floor);
            ground_contact2dp::spec s{cl.index_of("rod"), vec2dp(0.0, -0.2, 1.0)};
            s.kind = contact_kind2dp::flat;
            s.half_length = 0.5;
            size_t const c = gc.add(s);
            gc.engage(c);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 3);
            REQUIRE(gc.events().size() == 1);
            CHECK(gc.events()[0].touchdown);
            CHECK(gc.events()[0].impulse.empty());
            gc.engage(c); // already active: nothing happens
            CHECK(gc.events().size() == 1);
            for (int i = 0; i < 200; ++i)
                gc.step(1.0e-3);
            CHECK(gc.active(c));
            CHECK(at_rest(cl) < 1.0e-9);
        }
        { // promoted while still in the air: the weld is made on demand and engages at
          // touchdown -- the impact absorbs the whole motion, nothing moves afterwards
            closed_loop_system2dp cl;
            build(cl, 0.0);
            ground_contact2dp gc(cl, floor);
            size_t const c = gc.add({cl.index_of("rod"), vec2dp(0.0, 0.0, 1.0)});
            gc.set_kind(c, contact_kind2dp::flat);
            CHECK(gc.kind(c) == contact_kind2dp::flat);
            CHECK_FALSE(gc.active(c));
            for (int i = 0; i < 1000 && !gc.active(c); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 3);
            CHECK(at_rest(cl) < 1.0e-9);
        }
        { // promoted after landing as a pin: the switch itself is an impact that stops
          // the rotation about the pin; and back to a point leaves only the pin's rows
            closed_loop_system2dp cl;
            build(cl, 0.1);
            ground_contact2dp gc(cl, floor);
            size_t const c = gc.add({cl.index_of("rod"), vec2dp(-0.5, 0.0, 1.0)});
            for (int i = 0; i < 1000 && !gc.active(c); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 2);
            gc.set_kind(c, contact_kind2dp::flat);
            CHECK(cl.constraint_rows() == 3);
            CHECK(at_rest(cl) < 1.0e-9);
            gc.set_kind(c, contact_kind2dp::point);
            CHECK(cl.constraint_rows() == 2);
        }
    }

    TEST_CASE("ground_contact2dp: a segment released by hand, and lifted off by a pull")
    {
        value_t const g = 9.81;
        auto const floor =
            ground_contact2dp::ground_line(vec2dp(0.0, 0.0, 1.0), vec2dp(1.0, 0.0, 1.0));
        auto build = [&](closed_loop_system2dp& cl) {
            cl.system().set_gravity(vec2dp(0.0, -g, 0.0));
            cl.add_frame(static_frame2dp("W"));
            cl.add_body(static_frame2dp("rod", vec2dp(0.0, 0.3, 1.0), 0.1),
                        make_plate_body(2.0, 1.0, 0.02));
        };
        auto segment = [](size_t frame) {
            ground_contact2dp::spec sp;
            sp.frame = frame;
            sp.point_b = vec2dp(-0.5, 0.0, 1.0);
            sp.segment = true;
            sp.end_b = vec2dp(0.5, 0.0, 1.0);
            return sp;
        };

        { // landed on its lower end, then released by hand: no end is pinned any more
            closed_loop_system2dp cl;
            build(cl);
            ground_contact2dp gc(cl, floor);
            size_t const c = gc.add(segment(cl.index_of("rod")));
            for (int i = 0; i < 1000 && !gc.active(c); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            REQUIRE(gc.pinned_end(c) >= 0);
            gc.release(c);
            CHECK_FALSE(gc.active(c));
            CHECK(gc.pinned_end(c) == -1);
            CHECK(gc.kind(c) == contact_kind2dp::point);
        }
        { // pinned on an end, then pulled UP (gravity reversed): the lift-off rule lets
          // go
            closed_loop_system2dp cl;
            build(cl);
            ground_contact2dp gc(cl, floor);
            size_t const c = gc.add(segment(cl.index_of("rod")));
            for (int i = 0; i < 1000 && !gc.active(c); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.pinned_end(c) >= 0);
            cl.system().set_gravity(vec2dp(0.0, g, 0.0));
            for (int i = 0; i < 100 && gc.active(c); ++i)
                gc.step(1.0e-3);
            CHECK_FALSE(gc.active(c));
            CHECK_FALSE(gc.events().back().touchdown);
        }
        { // flat on both ends, then pulled up: it lifts off as a whole
            closed_loop_system2dp cl;
            build(cl);
            ground_contact2dp gc(cl, floor);
            size_t const c = gc.add(segment(cl.index_of("rod")));
            for (int i = 0; i < 3000 && !(gc.active(c) && gc.pinned_end(c) == -1); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            REQUIRE(gc.pinned_end(c) == -1); // welded flat
            cl.system().set_gravity(vec2dp(0.0, g, 0.0));
            for (int i = 0; i < 100 && gc.active(c); ++i)
                gc.step(1.0e-3);
            CHECK_FALSE(gc.active(c));
        }
    }

    TEST_CASE("geodesics: dms2deg reads the Unicode prime and double prime")
    {
        // U+2032 PRIME for minutes and U+2033 DOUBLE PRIME for seconds read like ' and "
        CHECK(dms2deg("52°31′12″N", geo_angle::latitude) ==
              doctest::Approx(dms2deg("52°31'12\"N", geo_angle::latitude)));
        CHECK(dms2deg("13° 24′ 36.5″ E", geo_angle::longitude) ==
              doctest::Approx(13.0 + 24.0 / 60.0 + 36.5 / 3600.0));
    }
}
