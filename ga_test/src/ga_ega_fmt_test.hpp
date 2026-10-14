// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_fmt_check.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// EGA formatters: every value type through fmt (compile-time and run-time spec) and
// through an ostream -- see ga_test_fmt_check.hpp for the contract
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA formatting")
{
    TEST_CASE("ega2d: formatting contract of every value type")
    {
        fmt_check::check<scalar2d>("Scalar2d");
        fmt_check::check<vec2d>("Vec2d");
        fmt_check::check<pscalar2d>("PScalar2d");
        fmt_check::check<mvec2d_e>("MVec2d_E");
        fmt_check::check<mvec2d>("MVec2d");
    }

    TEST_CASE("ega3d: formatting contract of every value type")
    {
        fmt_check::check<scalar3d>("Scalar3d");
        fmt_check::check<vec3d>("Vec3d");
        fmt_check::check<bivec3d>("BiVec3d");
        fmt_check::check<pscalar3d>("PScalar3d");
        fmt_check::check<mvec3d_e>("MVec3d_E");
        fmt_check::check<mvec3d_u>("MVec3d_U");
        fmt_check::check<mvec3d>("MVec3d");
    }
}
