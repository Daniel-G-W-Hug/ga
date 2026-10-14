// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_fmt_check.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// PGA formatters: every value type and every named geometric object through fmt
// (compile-time and run-time spec) and through an ostream -- see ga_test_fmt_check.hpp
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA formatting")
{
    TEST_CASE("pga2dp: formatting contract of every value type")
    {
        fmt_check::check<scalar2dp>("Scalar2dp");
        fmt_check::check<vec2dp>("Vec2dp");
        fmt_check::check<bivec2dp>("BiVec2dp");
        fmt_check::check<pscalar2dp>("PScalar2dp");
        fmt_check::check<mvec2dp_e>("MVec2dp_E");
        fmt_check::check<mvec2dp_u>("MVec2dp_U");
        fmt_check::check<mvec2dp>("MVec2dp");
        fmt_check::check<dualnum2dp>("DualNum2dp");
    }

    TEST_CASE("pga2dp: formatting contract of the geometric objects")
    {
        // each prints under its own name, not under the name of the type it wraps
        fmt_check::check_value(vector2d(vec2d(1, -2)), vector2d(vec2d(-0.0, -0.0)),
                               "Vector2d");
        fmt_check::check_value(point2d(1, -2), point2d(-0.0, -0.0), "Point2d");
        fmt_check::check_value(point2dp(vec2dp(1, -2, 0.5)),
                               point2dp(vec2dp(-0.0, -0.0, -0.0)), "Point2dp");
        fmt_check::check<line2d>("Line2d");
    }

    TEST_CASE("pga3dp: formatting contract of every value type")
    {
        fmt_check::check<scalar3dp>("Scalar3dp");
        fmt_check::check<vec3dp>("Vec3dp");
        fmt_check::check<bivec3dp>("BiVec3dp");
        fmt_check::check<trivec3dp>("TriVec3dp");
        fmt_check::check<pscalar3dp>("PScalar3dp");
        fmt_check::check<mvec3dp_e>("MVec3dp_E");
        fmt_check::check<mvec3dp_u>("MVec3dp_U");
        fmt_check::check<mvec3dp>("MVec3dp");
        fmt_check::check<dualnum3dp>("DualNum3dp");
    }

    TEST_CASE("pga3dp: formatting contract of the geometric objects")
    {
        fmt_check::check_value(vector3d(vec3d(1, -2, 0.5)),
                               vector3d(vec3d(-0.0, -0.0, -0.0)), "Vector3d");
        fmt_check::check<point3d>("Point3d");
        fmt_check::check_value(point3dp(vec3dp(1, -2, 0.5, 3)),
                               point3dp(vec3dp(-0.0, -0.0, -0.0, -0.0)), "Point3dp");
        fmt_check::check<line3d>("Line3d");
        fmt_check::check_value(plane3d(trivec3dp(1, -2, 0.5, 3)),
                               plane3d(trivec3dp(-0.0, -0.0, -0.0, -0.0)), "Plane3d");
    }
}
