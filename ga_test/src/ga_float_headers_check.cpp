// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Compile-only check: every public header of the library compiles with
// value_t = float (the build sets _HD_GA_VALUE_T_FLOAT for this translation unit, see
// ga/ga_value_t.hpp). Templates are checked where this includes instantiate them; the
// float variants of the base test suites instantiate the algebras' operations.

#if !defined(_HD_GA_VALUE_T_FLOAT)
#error "this check is meant to be compiled with _HD_GA_VALUE_T_FLOAT defined"
#endif

#include "ga/ga_cga.hpp"
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"
#include "ga/ga_sta.hpp"

#include "ga/ga_usr_fd.hpp"
#include "ga/ga_usr_geodesics.hpp"
#include "ga/ga_usr_magnetostatics.hpp"
#include "ga/ga_usr_types_mechanics.hpp"

#include "ga/detail/ga_quadrature.hpp"
#include "ga/detail/ga_segment_distance.hpp"
#include "ga/detail/ga_solver.hpp"
#include "ga/detail/ga_stencil.hpp"

#include <span>
#include <type_traits>
#include <vector>

static_assert(std::is_same_v<hd::ga::value_t, float>);
static_assert(std::is_same_v<decltype(hd::ga::eps), float const>);
static_assert(std::is_same_v<decltype(hd::ga::eps_congruent), float const>);

// the public templates no base suite instantiates, instantiated here under float
namespace {

using namespace hd::ga;

[[maybe_unused]] void instantiate_templates()
{
    auto const b2 = [](vec2d const& p) { return p; };
    auto const b3 = [](vec3d const& p) { return p; };
    std::vector<contour_point> const contour(1);
    std::vector<surface_point> const surface(1);
    (void)stress_wrench(b2, std::span<contour_point const>(contour));
    (void)stress_wrench(b3, std::span<surface_point const>(surface));
    (void)box_wrench(b3, magnetised_box{});

    auto const rhs = [](value_t, std::vector<value_t> const& u,
                        std::vector<value_t>& du) { du = u; };
    std::vector<value_t> u(2, value_t(1.0));
    (void)rk4_integrator(2).step(rhs, u, value_t(0.0), value_t(0.1));
    (void)abm2_integrator(2).step(rhs, u, value_t(0.0), value_t(0.1));
    (void)abm2_adaptive_integrator(2).integrate(rhs, u, value_t(0.0), value_t(1.0),
                                                value_t(0.1), value_t(1.0e-3),
                                                value_t(1.0e-3));
    std::vector<vec3dp> v(2), r(2);
    std::array<std::vector<vec3dp>, 2> uh{v, v};
    (void)rk4_step(v, uh, r, value_t(0.1), 1);
}

} // namespace
