// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga3dp mechanics formatters: the enums, the value structs, the frame tree and the
// kinematic system (whose formatter decodes each stored twist back into omega and v).
// Every spec is also given at RUN time: a compile-time checked format string runs the
// formatter's parse() only during compilation, so without fmt::runtime that code never
// executes.
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP mechanics formatting")
{
    // the same text whether the spec is checked at compile time or parsed at run time
    auto const same_at_runtime = [](auto const& x) {
        CHECK(fmt::format(fmt::runtime("{}"), x) == fmt::format("{}", x));
        CHECK(fmt::format(fmt::runtime("{:.2f}"), x) == fmt::format("{:.2f}", x));
    };

    TEST_CASE("pga3dp: the mechanics enums print their enumerator names")
    {
        CHECK(fmt::format("{}", joint3dp::free) == "free");
        CHECK(fmt::format("{}", joint3dp::revolute) == "revolute");
        CHECK(fmt::format("{}", joint3dp::prismatic) == "prismatic");
        CHECK(fmt::format("{}", joint3dp::helical) == "helical");
        CHECK(fmt::format("{}", joint3dp::cylindrical) == "cylindrical");
        CHECK(fmt::format("{}", joint3dp::spherical) == "spherical");
        CHECK(fmt::format("{}", joint3dp::planar) == "planar");
        CHECK(fmt::format("{}", constraint3dp::coincidence) == "coincidence");
        CHECK(fmt::format("{}", constraint3dp::distance) == "distance");
        CHECK(fmt::format("{}", constraint3dp::frame) == "frame");
        CHECK(fmt::format("{}", contact_kind3dp::point) == "point");
        CHECK(fmt::format("{}", contact_kind3dp::flat) == "flat");
        CHECK(fmt::format("{}", integrator_kind::rk4) == "rk4");
        CHECK(fmt::format("{}", integrator_kind::abm2) == "abm2");

        CHECK(fmt::format("[{:^9}]", joint3dp::helical) == "[ helical ]");
        CHECK(fmt::format(fmt::runtime("[{:>6}]"), integrator_kind::abm2) == "[  abm2]");
        CHECK_THROWS_AS((void)fmt::format(fmt::runtime("{:.2f}"), joint3dp::free),
                        fmt::format_error);
    }

    TEST_CASE("pga3dp: formatting of the mechanics value structs")
    {
        pose3dp const p{vec3dp(1, -2, 0.5, 1), vec3dp(0, 0, 0.25, 0)};
        CHECK(fmt::format("{}", p) ==
              "pose3dp(origin = Vec3dp(1,-2,0.5,1), rot = Vec3dp(0,0,0.25,0))");
        CHECK(
            fmt::format("{:.1f}", p) ==
            "pose3dp(origin = Vec3dp(1.0,-2.0,0.5,1.0), rot = Vec3dp(0.0,0.0,0.2,0.0))");
        same_at_runtime(p);

        kin_state3dp const k{.vel = vec3dp(1, 0, 0, 0),
                             .acc = vec3dp(0, -2, 0, 0),
                             .omega = vec3dp(0, 0, 3, 0),
                             .alpha = vec3dp(0.5, 0, 0, 0)};
        CHECK(fmt::format("{}", k) ==
              "kin_state3dp(vel = Vec3dp(1,0,0,0), acc = Vec3dp(0,-2,0,0), "
              "omega = Vec3dp(0,0,3,0), alpha = Vec3dp(0.5,0,0,0))");
        same_at_runtime(k);

        joint_state3dp j;
        j.type = joint3dp::helical;
        j.screw_b = bivec3dp(0, 0, 1, 0, 0, 0.5);
        j.rest = mvec3dp_e(1, 0, 0, 0, 0, 0, 0, 0);
        j.phi = -0.5;
        j.omega = 2;
        CHECK(fmt::format("{}", j) ==
              "joint_state3dp(type = helical, screw_b = BiVec3dp(0,0,1,0,0,0.5), "
              "rest = MVec3dp_E(1,0,0,0,0,0,0,0), phi = -0.5, omega = 2)");
        same_at_runtime(j);

        screw_axis3dp const s{bivec3dp(0, 0, 1, 0, 0, 0), 0.5, 2};
        CHECK(fmt::format("{}", s) ==
              "screw_axis3dp(axis = BiVec3dp(0,0,1,0,0,0), angle = 0.5, dist = 2)");
        same_at_runtime(s);

        screw_system3dp const sys{6, 6};
        CHECK(fmt::format("{:.3f}", sys) == "screw_system3dp(span = 6, closure = 6)");
        same_at_runtime(sys);

        loop_constraint3dp const c{
            1, vec3dp(1, 0, 0, 1), 3, vec3dp(0, 1, 0, 1), constraint3dp::frame, 0, true};
        CHECK(fmt::format("{}", c) ==
              "loop_constraint3dp(frame_a = 1, anchor_a = Vec3dp(1,0,0,1), frame_b = 3, "
              "anchor_b = Vec3dp(0,1,0,1), type = frame, length = 0, active = true)");
        same_at_runtime(c);

        body3dp const b = make_cuboid_body(2.0, 1.0, 2.0, 3.0);
        std::string const bs = fmt::format("{:.3f}", b);
        CHECK(bs.starts_with("body3dp(mass = 2.000, I = "));
        CHECK(bs.find(", I_inv = ") != std::string::npos);
        same_at_runtime(b);
    }

    TEST_CASE("pga3dp: formatting of a frame tree")
    {
        static_frame3dp const f("A", vec3dp(1, 2, 3, 1), vec3dp(0, 0, 0.5, 0));
        CHECK(fmt::format("{}", f) ==
              "static_frame3dp(name = 'A', origin = Vec3dp(1,2,3,1), "
              "rot = Vec3dp(0,0,0.5,0))");
        same_at_runtime(f);

        static_system3dp sys;
        CHECK(fmt::format("{}", sys) == "static_system3dp()");

        sys.add_frame(static_frame3dp("W"));
        sys.add_frame(static_frame3dp("B", vec3dp(1, 0, 0, 1)));
        // one line per frame; the root names itself as its parent
        CHECK(fmt::format("{:.1f}", sys) ==
              "static_system3dp(\n"
              "    rf[0] (parent 0) = static_frame3dp(name = 'W', "
              "origin = Vec3dp(0.0,0.0,0.0,1.0), rot = Vec3dp(0.0,0.0,0.0,0.0))\n"
              "    rf[1] (parent 0) = static_frame3dp(name = 'B', "
              "origin = Vec3dp(1.0,0.0,0.0,1.0), rot = Vec3dp(0.0,0.0,0.0,0.0))\n"
              ")");
        same_at_runtime(sys);
    }

    TEST_CASE("pga3dp: formatting of a kinematic system decodes omega and v")
    {
        kinematic_system3dp ks;
        CHECK(fmt::format("{}", ks) == "kinematic_system3dp()");

        // a frame AT its parent's origin, so the relative twist carries the state as
        // given: omega about e3, v along e1, alpha about e1, a along -e2
        ks.add_frame(static_frame3dp("W"));
        ks.add_frame(static_frame3dp("B"), kin_state3dp{.vel = vec3dp(1, 0, 0, 0),
                                                        .acc = vec3dp(0, -2, 0, 0),
                                                        .omega = vec3dp(0, 0, 3, 0),
                                                        .alpha = vec3dp(0.5, 0, 0, 0)});
        std::string const s = fmt::format("{}", ks);
        CHECK(s.starts_with("kinematic_system3dp(\n    rf[0] (parent 0) = "));
        CHECK(s.find("  | rel omega = (0, 0, 0), v = (0, 0, 0); alpha = (0, 0, 0), "
                     "a = (0, 0, 0)\n") != std::string::npos);
        CHECK(s.find("  | rel omega = (0, 0, 3), v = (1, 0, 0); alpha = (0.5, 0, 0), "
                     "a = (0, -2, 0)\n") != std::string::npos);
        CHECK(s.ends_with("\n)"));
        same_at_runtime(ks);
    }

    TEST_CASE("pga3dp: a numeric spec reaches every field of a joint's range and drive")
    {
        joint_range3dp r;
        r.lo = -1.5;
        r.hi = 2.25;
        r.k_stop = 100.0;
        r.c_stop = 0.5;
        CHECK(fmt::format("{:.2f}", r) == "joint_range3dp(lo = -1.50, hi = 2.25, "
                                          "k_stop = 100.00, c_stop = 0.50)");
        same_at_runtime(r);

        joint_drive3dp dr;
        dr.tau_max = 3.0;
        dr.qd_max = 4.5;
        dr.qdd_max = 6.0;
        dr.armature = 0.125;
        dr.actuated = false;
        CHECK(fmt::format("{:.1f}", dr) ==
              "joint_drive3dp(tau_max = 3.0, qd_max = 4.5, qdd_max = 6.0, "
              "armature = 0.1, actuated = false)");
        same_at_runtime(dr);
    }
}
