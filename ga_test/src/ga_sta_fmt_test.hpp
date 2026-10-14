// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_fmt_check.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// STA formatters: every value type through fmt (compile-time and run-time spec) and
// through an ostream -- see ga_test_fmt_check.hpp for the contract
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("STA formatting")
{
    TEST_CASE("sta4ds: formatting contract of every value type")
    {
        fmt_check::check<scalar4ds>("Scalar4ds");
        fmt_check::check<vec4ds>("Vec4ds");
        fmt_check::check<bivec4ds>("BiVec4ds");
        fmt_check::check<trivec4ds>("TriVec4ds");
        fmt_check::check<pscalar4ds>("PScalar4ds");
        fmt_check::check<mvec4ds_e>("MVec4ds_E");
        fmt_check::check<mvec4ds_u>("MVec4ds_U");
        fmt_check::check<mvec4ds>("MVec4ds");
    }
}
