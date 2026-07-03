// Module entrypoint for ga_py — hand-written, small, stable.
// Everything else under src/generated/ is produced by ga_bindgen.

#include <nanobind/nanobind.h>

namespace nb = nanobind;

// Implemented in src/generated/register_all.cpp (generated).
// `top` receives free functions from `hd::ga`; `ega` / `pga` / `sta` receive
// types and free functions from their respective sub-namespaces.
void register_all(nb::module_& top, nb::module_& ega, nb::module_& pga,
                  nb::module_& sta);

// Hand-written bindings — see the corresponding bindings_*.cpp for rationale.
void bind_rk4_step(nb::module_& m);
void bind_mechanics(nb::module_& pga);

NB_MODULE(_ga_py, m) {
    m.doc() = "Python bindings for the hd::ga geometric algebra library.\n"
              "Generated from C++ headers by ga_bindgen.";

    auto ega = m.def_submodule("ega", "Euclidean geometric algebra (2D, 3D)");
    auto pga = m.def_submodule("pga", "Projective geometric algebra (2dp, 3dp)");
    auto sta = m.def_submodule("sta", "Space-time algebra (G(1,3,0), types *_4ds)");

    register_all(m, ega, pga, sta);
    bind_rk4_step(m);
    bind_mechanics(pga);
}
