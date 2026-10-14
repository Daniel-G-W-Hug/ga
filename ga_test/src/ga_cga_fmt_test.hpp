// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_fmt_check.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// CGA formatters: every value type through fmt (compile-time and run-time spec) and
// through an ostream -- see ga_test_fmt_check.hpp for the contract
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("CGA formatting")
{
    TEST_CASE("cga2dc: formatting contract of every value type")
    {
        fmt_check::check<scalar2dc>("Scalar2dc");
        fmt_check::check<vec2dc>("Vec2dc");
        fmt_check::check<bivec2dc>("BiVec2dc");
        fmt_check::check<trivec2dc>("TriVec2dc");
        fmt_check::check<pscalar2dc>("PScalar2dc");
        fmt_check::check<mvec2dc_e>("MVec2dc_E");
        fmt_check::check<mvec2dc_u>("MVec2dc_U");
        fmt_check::check<mvec2dc>("MVec2dc");
    }

    TEST_CASE("cga3dc: formatting contract of every value type")
    {
        fmt_check::check<scalar3dc>("Scalar3dc");
        fmt_check::check<vec3dc>("Vec3dc");
        fmt_check::check<bivec3dc>("BiVec3dc");
        fmt_check::check<trivec3dc>("TriVec3dc");
        fmt_check::check<quadvec3dc>("QuadVec3dc");
        fmt_check::check<pscalar3dc>("PScalar3dc");
        fmt_check::check<mvec3dc_e>("MVec3dc_E");
        fmt_check::check<mvec3dc_u>("MVec3dc_U");
        fmt_check::check<mvec3dc>("MVec3dc");
    }
}
