// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// pga2dp mechanics formatters: the enums, the value structs, and the spec forwarding of
// the composite formatters (one captured spec handed down to every numeric member).
// Every spec is also given at RUN time: a compile-time checked format string runs the
// formatter's parse() only during compilation, so without fmt::runtime that code never
// executes.
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA2DP mechanics formatting")
{
    // the same text whether the spec is checked at compile time or parsed at run time
    auto const same_at_runtime = [](auto const& x) {
        CHECK(fmt::format(fmt::runtime("{}"), x) == fmt::format("{}", x));
        CHECK(fmt::format(fmt::runtime("{:.2f}"), x) == fmt::format("{:.2f}", x));
    };

    TEST_CASE("pga2dp: the mechanics enums print their enumerator names")
    {
        CHECK(fmt::format("{}", joint2dp::free) == "free");
        CHECK(fmt::format("{}", joint2dp::revolute) == "revolute");
        CHECK(fmt::format("{}", joint2dp::prismatic) == "prismatic");
        CHECK(fmt::format("{}", constraint2dp::coincidence) == "coincidence");
        CHECK(fmt::format("{}", constraint2dp::distance) == "distance");
        CHECK(fmt::format("{}", constraint2dp::frame) == "frame");
        CHECK(fmt::format("{}", contact_kind2dp::point) == "point");
        CHECK(fmt::format("{}", contact_kind2dp::flat) == "flat");

        // a string spec applies to the name: width, fill and alignment
        CHECK(fmt::format("[{:>9}]", joint2dp::revolute) == "[ revolute]");
        CHECK(fmt::format(fmt::runtime("[{:*<7}]"), contact_kind2dp::flat) ==
              "[flat***]");

        // a numeric spec is rejected, not ignored
        CHECK_THROWS_AS((void)fmt::format(fmt::runtime("{:.2f}"), joint2dp::free),
                        fmt::format_error);
    }

    TEST_CASE("pga2dp: formatting of the mechanics value structs")
    {
        pose2dp const p{vec2dp(1, -2, 1), 0.5};
        CHECK(fmt::format("{}", p) == "pose2dp(origin = Vec2dp(1,-2,1), phi = 0.5)");
        CHECK(fmt::format("{:.2f}", p) ==
              "pose2dp(origin = Vec2dp(1.00,-2.00,1.00), phi = 0.50)");
        same_at_runtime(p);
        CHECK(fmt::format("{}", pose2dp{vec2dp(-0.0, 0, 1), -0.0}) ==
              "pose2dp(origin = Vec2dp(0,0,1), phi = 0)");

        joint_state2dp j;
        j.type = joint2dp::revolute;
        j.rest = mvec2dp_u(0, 0, 0.5, 1);
        j.phi = 0.25;
        j.omega = -1;
        CHECK(fmt::format("{}", j) ==
              "joint_state2dp(type = revolute, screw_b = Vec2dp(0,0,1), "
              "rest = MVec2dp_U(0,0,0.5,1), phi = 0.25, omega = -1)");
        CHECK(fmt::format("{:.1f}", j) ==
              "joint_state2dp(type = revolute, screw_b = Vec2dp(0.0,0.0,1.0), "
              "rest = MVec2dp_U(0.0,0.0,0.5,1.0), phi = 0.2, omega = -1.0)");
        same_at_runtime(j);

        screw_axis2dp const s{vec2dp(1, 2, 1), 0.5, vec2dp{}};
        CHECK(fmt::format("{}", s) ==
              "screw_axis2dp(centre = Vec2dp(1,2,1), angle = 0.5, "
              "translation = Vec2dp(0,0,0))");
        same_at_runtime(s);

        // the dimensions are counts: a numeric spec is accepted and does not apply
        screw_system2dp const sys{2, 3};
        CHECK(fmt::format("{}", sys) == "screw_system2dp(span = 2, closure = 3)");
        CHECK(fmt::format("{:.2f}", sys) == "screw_system2dp(span = 2, closure = 3)");
        same_at_runtime(sys);

        loop_constraint2dp const c{
            0, vec2dp(1, 0, 1), 2, vec2dp(0, 1, 1), constraint2dp::distance, 1.5, false};
        CHECK(fmt::format("{}", c) ==
              "loop_constraint2dp(frame_a = 0, anchor_a = Vec2dp(1,0,1), frame_b = 2, "
              "anchor_b = Vec2dp(0,1,1), type = distance, length = 1.5, active = false)");
        CHECK(fmt::format("{:.1f}", c) ==
              "loop_constraint2dp(frame_a = 0, anchor_a = Vec2dp(1.0,0.0,1.0), "
              "frame_b = 2, anchor_b = Vec2dp(0.0,1.0,1.0), type = distance, "
              "length = 1.5, active = false)");
        same_at_runtime(c);

        // a body prints its mass and both inertia maps, the spec reaching all three
        body2dp const b = make_disc_body(2.0, 0.5);
        std::string const bs = fmt::format("{:.3f}", b);
        CHECK(bs.starts_with("body2dp(mass = 2.000, I = "));
        CHECK(bs.find(", I_inv = ") != std::string::npos);
        CHECK(bs.ends_with(")"));
        CHECK(bs.find("2.000") != std::string::npos);
        same_at_runtime(b);
    }

    TEST_CASE("pga2dp: a numeric spec reaches every field of a joint's range and drive")
    {
        joint_range2dp r;
        r.lo = -1.5;
        r.hi = 2.25;
        r.k_stop = 100.0;
        r.c_stop = 0.5;
        CHECK(fmt::format("{:.2f}", r) == "joint_range2dp(lo = -1.50, hi = 2.25, "
                                          "k_stop = 100.00, c_stop = 0.50)");
        same_at_runtime(r);

        joint_drive2dp dr;
        dr.tau_max = 3.0;
        dr.qd_max = 4.5;
        dr.qdd_max = 6.0;
        dr.armature = 0.125;
        dr.actuated = false;
        CHECK(fmt::format("{:.1f}", dr) ==
              "joint_drive2dp(tau_max = 3.0, qd_max = 4.5, qdd_max = 6.0, "
              "armature = 0.1, actuated = false)");
        same_at_runtime(dr);

        geo_pos_dms const g{"52°31'12\"N", "13°24'36\"E", 34.5};
        CHECK(fmt::format("{:.2f}", g) ==
              "geo_pos_dms(lat = 52°31'12\"N, lon = 13°24'36\"E, height = 34.50 m)");
        same_at_runtime(g);
        geo_pos_dms2dp const g2{"52°31'12\"N", 34.5};
        CHECK(fmt::format("{:.1f}", g2).ends_with("height = 34.5 m)"));
        same_at_runtime(g2);
    }
}
