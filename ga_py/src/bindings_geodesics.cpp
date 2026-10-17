// Hand-written bindings for the geodesy surface of ga_usr_geodesics.hpp:
//
// - the reference ellipsoid. It carries member functions (flattening(), N(), ...), so it
//   is no pure-data struct and the bindgen does not list it -- and every geodesy
//   function takes an `ellipsoid const&`. ga_bindgen/src/emit_nanobind.py registers the
//   name as hand-bound, so with this class the generated layer binds geo_to_ecef,
//   enu_at, the projections onto the tangent plane, tm_forward / tm_inverse, ... and the
//   generated constants wgs84 / grs80. Those constants are converted at import, so this
//   class must be registered BEFORE register_all (see module.cpp).
//
// - unproject(m, x): the C++ form reports a miss through a `bool*` out-parameter, which
//   has no Python mapping; bound here returning the position, or None when the map
//   coordinates are outside the image of the ellipsoid (beyond the limb).
//
// The ellipsoid lives in ga_py.ega beside geo_pos and transverse_mercator, where the
// bindgen routes the plain value types of hd::ga.

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>

#include <fmt/format.h>

#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include "ga/ga_pga.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace hd::ga;
using namespace hd::ga::pga;

void bind_geodesics_types(nb::module_& ega)
{
    nb::class_<ellipsoid>(
        ega, "ellipsoid",
        "Reference ellipsoid of revolution, given by its two radii [m].")
        .def(
            "__init__",
            [](ellipsoid* e, double r_equator, double r_pole) {
                new (e) ellipsoid{r_equator, r_pole};
            },
            "r_equator"_a, "r_pole"_a)
        .def_rw("r_equator", &ellipsoid::r_equator, "semi-major axis a [m]")
        .def_rw("r_pole", &ellipsoid::r_pole, "semi-minor axis b [m]")
        .def("flattening", &ellipsoid::flattening, "f = (a - b) / a")
        .def("e_sq", &ellipsoid::e_sq, "first eccentricity squared (a^2 - b^2) / a^2")
        .def("radius", &ellipsoid::radius, "lat"_a,
             "distance of the surface from the geocenter at geodetic latitude lat [rad]")
        .def("N", &ellipsoid::N, "lat"_a,
             "radius of curvature in the prime vertical at geodetic latitude lat [rad]")
        .def("geocentric_lat", &ellipsoid::geocentric_lat, "lat"_a,
             "the geocentric latitude for a geodetic one [rad]")
        .def("__eq__",
             [](ellipsoid const& a, ellipsoid const& b) {
                 return a.r_equator == b.r_equator && a.r_pole == b.r_pole;
             })
        .def("__repr__", [](ellipsoid const& e) { return fmt::format("{}", e); })
        .def("__str__", [](ellipsoid const& e) { return fmt::format("{}", e); })
        .def("__format__", [](ellipsoid const& e, std::string_view spec) {
            try {
                if (spec.empty()) return fmt::format("{}", e);
                return fmt::format(fmt::runtime("{:" + std::string(spec) + "}"), e);
            }
            catch (fmt::format_error const& err) {
                throw std::invalid_argument(err.what());
            }
        });
}

void bind_geodesics(nb::module_& pga)
{
    pga.def(
        "unproject",
        [](azimuthal_map const& m, vec2d const& x) -> std::optional<geo_pos> {
            bool hit = false;
            geo_pos const p = unproject(m, x, &hit);
            if (!hit) return std::nullopt;
            return p;
        },
        "m"_a, "x"_a,
        "The position on the ellipsoid that maps to the map coordinates x, or None when "
        "x lies outside the image of the ellipsoid (beyond the limb).");
}
