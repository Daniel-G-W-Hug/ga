// Hand-written bindings for hd::ga::rk4_step.
//
// Why hand-written: rk4_step is templated on `VecType`. The libclang scanner
// records a single overload per template, with `VecType` as a literal type
// parameter — the bindgen's resolver only enumerates *concrete* types from
// the manifest's type_map, so it cannot synthesise per-VecType bindings.
// Enumerating explicitly here is the simplest path; it also serves as the
// project's first hand-written nanobind binding (see ga_py/README.md §6.5).
//
// What's bound: the std::vector / std::array overload of rk4_step (the
// mdspan overload is unreachable from Python — nanobind has no mdspan
// caster). VecType is enumerated for the four vector-shaped GA types
// commonly used in physics simulations: vec2d, vec3d, vec2dp, vec3dp.
// Bivectors / trivectors are mechanical follow-ups (one m.def per type).
//
// Surface: each binding lives at top level (`ga_py.rk4_step`). Overload
// resolution is by argument type — calling `rk4_step([vec3d(...)], ...)`
// dispatches to the vec3d instantiation automatically.

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <utility>
#include <vector>

#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

namespace nb = nanobind;
using namespace hd::ga;

namespace {

// Single-template helper that emits the m.def(...) lambda for one VecType.
// Defined as a function template so each VecType gets its own concrete
// instantiation — nanobind requires concrete signatures at the binding site.
template <typename VecType>
void bind_rk4_step_one(nb::module_& m)
{
    m.def("rk4_step",
          [](std::vector<VecType> u, std::array<std::vector<VecType>, 2> uh,
             std::vector<VecType> const& rhs, value_t dt, size_t rk_step) {
              return rk4_step<VecType>(std::move(u), std::move(uh), rhs, dt,
                                       rk_step);
          },
          nb::arg("u"), nb::arg("uh"), nb::arg("rhs"), nb::arg("dt"),
          nb::arg("rk_step"),
          "Runge-Kutta 4 sub-step (vector overload). Returns "
          "(u_new, (uh0_new, uh1_new)). VecType is dispatched on the element "
          "type of `u`.");
}

} // anonymous namespace

void bind_rk4_step(nb::module_& m)
{
    bind_rk4_step_one<vec2d>(m);
    bind_rk4_step_one<vec3d>(m);
    bind_rk4_step_one<vec2dp>(m);
    bind_rk4_step_one<vec3dp>(m);
}
