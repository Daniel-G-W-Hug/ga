// Module entrypoint for ga_py — hand-written, small, stable.
// Everything else under src/generated/ is produced by ga_bindgen.

#include <nanobind/nanobind.h>

namespace nb = nanobind;

// Implemented in src/generated/register_all.cpp (generated).
// `top` receives free functions from `hd::ga`; `ega` / `pga` receive
// types and free functions from their respective sub-namespaces.
void register_all(nb::module_& top, nb::module_& ega, nb::module_& pga);

// Hand-written bindings — see src/bindings_rk4_step.cpp for the rationale.
void bind_rk4_step(nb::module_& m);

NB_MODULE(_ga_py, m) {
    m.doc() = "Python bindings for the hd::ga geometric algebra library.\n"
              "Generated from C++ headers by ga_bindgen.";

    auto ega = m.def_submodule("ega", "Euclidean geometric algebra (2D, 3D, 4D)");
    auto pga = m.def_submodule("pga", "Projective geometric algebra (2dp, 3dp)");

    register_all(m, ega, pga);
    bind_rk4_step(m);
}
