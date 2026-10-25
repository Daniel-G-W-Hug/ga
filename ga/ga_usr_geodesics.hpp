#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <algorithm> // std::max
#include <array>     // std::array (transverse Mercator coefficients)
#include <cctype>    // std::isspace, std::toupper
#include <cmath>     // std::sin, std::cos, std::sqrt, std::atan2, std::abs
#include <cstdio>    // std::snprintf (deg2dms)
#include <cstdlib>   // std::strtod
#include <stdexcept> // std::invalid_argument
#include <string>    // std::string
#include <utility>   // std::swap
#include <vector>    // std::vector (densify)

#include "detail/ga_quadrature.hpp" // gauss_legendre12 (the geodesic integrals)
#include "ega/ga_ega2d_ops.hpp"  // normalize, operator*(vec, I_2d), ... for the 2D frame
#include "ega/ga_ega3d_ops.hpp"  // wdg, dual, normalize, ... for the direction frame
#include "ga_usr_consts.hpp"     // e3_3d, I_3d, O_3dp, x_axis_3dp, z_axis_3dp, ...
#include "ga_usr_types.hpp"      // vec3d, vec3dp, mvec3dp_e
#include "ga_usr_utilities.hpp"  // deg2rad, rad2deg, pi
#include "ga_value_t.hpp"        // value_t
#include "pga/ga_pga2dp_ops.hpp" // get_motor, move2dp, unitize, ... for the 2D case
#include "pga/ga_pga3dp_ops.hpp" // get_motor, move3dp, unitize, ... for points/motors

/////////////////////////////////////////////////////////////////////////////////////////
// Geodetic coordinates on a reference ellipsoid.
//
// The two coordinate systems related here (both right-handed):
//
// ECEF: Earth-Centered, Earth-Fixed
// -> cs rotates with earth, a point fixed on the earth's surface has fixed coordinates
// -> origin of cs in center of mass of earth (geocenter) at (0,0,0,1)
// -> e1 axis (X) in equatorial plane through prime meridian (0° longitude)
// -> e2 axis (Y) in equatorial plane 90° east of e1 (90° East longitude)
// -> e3 axis (Z) aligns with the earth's mean rotational axis, pointing north
//
// ENU: East-North-Up
// -> local Cartesian cs attached to one point on the earth's surface
// -> e1 axis (X) points locally east
// -> e2 axis (Y) points locally north
// -> e3 axis (Z) points locally upwards
//
// The user-facing input is what one reads off an encyclopedia entry for a town: two
// angles in degree/minute/second notation and an elevation in m above sea level. That
// input form (geo_pos_dms) is kept as given; the calculation form (geo_pos, angles in
// rad) is obtained from it by an explicit conversion, so the parsing cost is paid once,
// before the calculation, and never inside it.
//
// Conventions committed to throughout (see geo_pos and enu_at for the detail): the
// latitude is the GEODETIC one, the height is measured along the ELLIPSOID NORMAL, and
// ENU's up is that same normal -- the local vertical, as used in practice. None of these
// is the geocentric ray, which sits up to ~11.5' away from the normal.
//
// provides in namespace hd::ga:
//
// - ellipsoid                 : reference ellipsoid, given by its two radii; carries
//                               radius(), N(), flattening(), e_sq(), geocentric_lat()
// - wgs84 / grs80             : the two reference ellipsoids provided
// - geo_angle                 : which angle is being parsed (latitude / longitude)
// - geo_pos_dms               : position as the user provides it (two strings +
//                               elevation), left unparsed
// - geo_pos                   : position ready to calculate with (lat/lon in rad +
//                               height along the ellipsoid normal)
// - dms2deg()                 : parse one dms-notation angle into decimal degrees
// - dms2rad()                 : the same, in rad
// - deg2dms() / rad2dms()     : the inverse -- write an angle back out in dms
//                               notation, closing the parser's round trip
// - to_geo_pos()              : geo_pos_dms -> geo_pos (the explicit conversion above)
// - densify()                 : intermediate positions on an edge straight in (lat, lon)
// - transverse_mercator       : the parameters of a transverse Mercator grid
// - utm_zone()                : those of a UTM zone
// - tm_forward() / tm_inverse(): geodetic position <-> grid coordinates (easting,
//                               northing) -- e.g. to read geodata delivered in UTM
// - lambert_azimuthal         : the parameters of a Lambert azimuthal equal-area map
// - laea_at() / etrs89_laea() : that map centred at a position / the EU's EPSG:3035 grid
// - laea_forward() / laea_inverse(): geodetic position <-> map coordinates, AREA kept --
//                               maps for statistics and densities
// - geodesic_inverse()        : distance and both azimuths of the shortest path between
//                               two positions (geodesic_solution)
// - geodesic_direct()         : where a path of given start, azimuth and length ends
//                               (geodesic_destination)
// - geodesic_path()           : a shortest path as positions, to draw a route
// - azimuthal_equidistant     : the parameters of an azimuthal equidistant map; aeqd_at()
// - aeqd_forward() / aeqd_inverse(): geodetic position <-> map coordinates, DISTANCE and
//                               DIRECTION from the centre kept -- the flight map of a hub
// - distance_from_geocenter() : the position's "total radius" [m], derived
// - geo_pos_dms2dp/geo_pos2dp : the same pair for the MERIDIAN SECTION (no longitude);
//                               to_geo_pos() and distance_from_geocenter() overload on
//                               them
// - detail::meridian_from_geodetic() / geodetic_from_meridian() : the shared core both
//                               dimensions are built on -- 3D is this map plus one
//                               rotation by the longitude
//
// provides in namespace hd::ga::pga:
//
// - geo_to_ecef()        : geo_pos -> ECEF, as a pga3dp point
// - ecef_to_geo()        : ECEF point -> geo_pos (Bowring's closed form)
// - enu_frame            : the local east/north/up directions, expressed in ECEF
// - enu_at()             : that frame from the GEOMETRY -- meridian plane, contraction
//                          and dual; no transcendental function is evaluated. Overloads
//                          for an ECEF point and for a geo_pos.
// - enu_basis_at()       : the same frame written out from the ANGLES -- the conventional
//                          closed form, kept alongside enu_at() so the two cross-check
// - enu_to_ecef_motor()  : motor taking ENU coordinates to ECEF coordinates (its rrev is
//                          the ECEF -> ENU direction)
// - enu_motor_at()       : the same motor for a station given as an ECEF point
// - helmert7             : the seven parameters of a datum transformation, built in the
//                          published units (m, arc seconds, ppm) by make_helmert7()
// - helmert_position_vector() / helmert_coordinate_frame() : ECEF -> ECEF between two
//                          datums, one motor sandwich plus a scale; the two EPSG
//                          conventions differ in the sign of the rotation, so each has
//                          its own name
// - azimuthal            : the projections onto the tangent plane -- orthographic,
//                          gnomonic, stereographic, perspective (what each preserves and
//                          what it is used for is stated at the enum)
// - azimuthal_options    : the per-kind settings (viewpoint height, drawable radius)
// - azimuthal_map        : one projection at one station, built by make_azimuthal_map()
// - to_local()           : an ECEF point in the map's local frame
// - is_drawable()        : inside the map's drawable region (visible, for the horizon
// kinds)
// - project()            : local point (or geo_pos) -> map coordinates [m], the meet of
//                          the line from the projection centre with the tangent plane
// - clip_at()            : where a segment crosses a plane -- cutting an outline at the
//                          map's limit
// - unproject()          : map coordinates -> the position on the ellipsoid
// - limit_outline()      : points around the edge of the drawable region (the horizon)
//
// and the two-dimensional counterparts, for a position in the meridian section:
//
// - geo_to_ecef() / ecef_to_geo()  : overloads on geo_pos2dp / vec2dp
// - un_frame / un_at()   : the local frame, in the order (UP, NORTH) -- forced, since
//                          (north, up) is negatively oriented and no motor could carry
//                          the reference axes onto it. No degenerate case but the
//                          geocenter: the 3D pole degeneracy comes from the longitude.
// - un_basis_at()        : the same frame from the angle, as the cross-check
// - un_to_ecef_motor()   : ONE rotation -- and it is the latitude itself
// - un_motor_at()        : the same motor for a station given as a point
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// reference ellipsoid
/////////////////////////////////////////////////////////////////////////////////////////

// The earth as an ellipsoid of revolution: wider at the equator than at the poles. It is
// given by its two radii, and its radius at an intermediate latitude follows from them.
//
// With the geodetic latitude lat, the radius measured from the geocenter is
//
//     R(lat) = sqrt( ((a^2 cos(lat))^2 + (b^2 sin(lat))^2)
//                  / (( a  cos(lat))^2 + ( b  sin(lat))^2) )
//
// which runs from R(0) = a at the equator down to R(pi/2) = b at the pole.
struct ellipsoid {

    value_t r_equator; // semi-major axis a [m]
    value_t r_pole;    // semi-minor axis b [m]

    // flattening f = (a - b)/a
    constexpr value_t flattening() const { return (r_equator - r_pole) / r_equator; }

    // first eccentricity squared e^2 = (a^2 - b^2)/a^2
    constexpr value_t e_sq() const
    {
        return (r_equator * r_equator - r_pole * r_pole) / (r_equator * r_equator);
    }

    // radius from the geocenter at geodetic latitude lat [rad]
    value_t radius(value_t lat) const
    {
        value_t const a = r_equator;
        value_t const b = r_pole;
        value_t const c = std::cos(lat);
        value_t const s = std::sin(lat);

        value_t const num = (a * a * c) * (a * a * c) + (b * b * s) * (b * b * s);
        value_t const den = (a * c) * (a * c) + (b * s) * (b * s);

        return std::sqrt(num / den);
    }

    // The geodetic latitude is the angle of the ellipsoid NORMAL against the equatorial
    // plane, the geocentric one the angle of the ray from the geocenter. On a sphere the
    // two coincide; on the earth they differ by up to ~11.5' (at 45°), so a position
    // placed along a ray must use the geocentric angle:
    //
    //     tan(lat_geocentric) = (1 - e^2) tan(lat_geodetic)
    //
    value_t geocentric_lat(value_t lat) const
    {
        return std::atan2((1.0 - e_sq()) * std::sin(lat), std::cos(lat));
    }

    // Radius of curvature in the prime vertical at geodetic latitude lat [rad]: the
    // distance from the position along the ellipsoid NORMAL down to the polar axis.
    // This -- not radius() -- is what places a position in ECEF, because a geodetic
    // height is measured along that normal.
    value_t N(value_t lat) const
    {
        value_t const s = std::sin(lat);
        return r_equator / std::sqrt(1.0 - e_sq() * s * s);
    }
};

// Exact equality of the two radii. Free, not a member: the binding generator drops a
// struct carrying a member function, and the generated equality of every struct with an
// ellipsoid field (azimuthal_map) compares it with this.
inline bool operator==(ellipsoid const& a, ellipsoid const& b)
{
    return a.r_equator == b.r_equator && a.r_pole == b.r_pole;
}

// WGS84 (the GPS datum) and GRS80 (the ITRS/ETRS89 datum); the two differ by well below
// a millimetre on the earth's surface
inline ellipsoid const wgs84{6378137.0, 6356752.314245179};
inline ellipsoid const grs80{6378137.0, 6356752.314140356};


/////////////////////////////////////////////////////////////////////////////////////////
// the meridian section -- the core both dimensions are built on
/////////////////////////////////////////////////////////////////////////////////////////
//
// All of the ellipsoid's geometry lives in ONE plane. Cut the ellipsoid of revolution
// along a meridian and what remains is an ellipse with the same two radii, in which a
// position is described by the geodetic latitude and a height along the normal:
//
//     r = (N + h) cos(lat)              distance from the polar axis
//     z = (N (1 - e^2) + h) sin(lat)    distance along it
//
// The three-dimensional case is this map plus a rotation by the longitude about the
// polar axis -- geo_to_ecef() below is literally these two lines with r resolved into
// (r cos(lon), r sin(lon)). The two-dimensional case IS the meridian section, so it uses
// them unchanged. Sharing them is therefore not code reuse for its own sake: it is the
// statement that 3D geodesy is 2D geodesy plus one rotation.

namespace detail {

struct meridian_rz {
    value_t r; // distance from the polar axis [m]
    value_t z; // distance along the polar axis [m]
};

struct meridian_lat_h {
    value_t lat;    // geodetic latitude [rad]
    value_t height; // height above the ellipsoid, along its normal [m]
};

// geodetic (lat, height) -> the meridian-plane coordinates (r, z)
inline meridian_rz meridian_from_geodetic(ellipsoid const& el, value_t lat,
                                          value_t height)
{
    value_t const N = el.N(lat);

    return meridian_rz{(N + height) * std::cos(lat),
                       (N * (1.0 - el.e_sq()) + height) * std::sin(lat)};
}

// (r, z) -> geodetic (lat, height), via Bowring's closed-form solution (1976): the
// auxiliary (parametric) latitude theta = atan2(z a, r b) turns the implicit relation
// between geodetic latitude and height into an explicit one.
//
// It is a SINGLE pass, not an iteration, so its accuracy depends on the height -- the
// approximation is built around a point near the surface. Measured round-trip error at
// latitude 45 deg:
//
//     height        latitude error      height error
//     0             2e-11"              4e-11 m       (exact, to double precision)
//     1 km          3e-10"              1e-08 m
//     100 km        3e-06"              9e-05 m
//     1000 km       2e-04"              7e-03 m
//     36000 km      1e-03"              0.26 m        (geostationary)
//
// So it is exact for anything on or near the ground and still good to a metre out at
// geostationary altitude -- but the sub-millimetre claim holds only up to roughly
// 100 km. Iterate (feed the result back through) if a satellite altitude is needed to
// better than that.
//
// r < 0 is not meaningful (it is a distance from the axis); on the axis itself, r == 0,
// the latitude is +/- 90 deg and the height follows from z alone.
inline meridian_lat_h geodetic_from_meridian(ellipsoid const& el, value_t r, value_t z)
{
    value_t const a = el.r_equator;
    value_t const b = el.r_pole;
    value_t const e2 = el.e_sq();
    value_t const ep2 = e2 / ((1.0 - el.flattening()) * (1.0 - el.flattening()));

    if (r < eps * a) { // on the polar axis
        return meridian_lat_h{(z >= 0.0) ? 0.5 * pi : -0.5 * pi, std::abs(z) - b};
    }

    value_t const theta = std::atan2(z * a, r * b);
    value_t const st = std::sin(theta);
    value_t const ct = std::cos(theta);

    value_t const lat = std::atan2(z + ep2 * b * st * st * st, r - e2 * a * ct * ct * ct);

    value_t const N = el.N(lat);
    value_t const sp = std::sin(lat);

    // near the poles r -> 0 makes h = r/cos(lat) - N ill-conditioned; switch to the
    // z-based form there (each is well-conditioned where the other is not)
    value_t const h =
        (std::abs(sp) > 0.5) ? z / sp - N * (1.0 - e2) : r / std::cos(lat) - N;

    return meridian_lat_h{lat, h};
}

} // namespace detail


/////////////////////////////////////////////////////////////////////////////////////////
// parsing degree/minute/second notation
/////////////////////////////////////////////////////////////////////////////////////////

// which of the two angles is being parsed -- fixes the admissible hemisphere letters
// and the admissible range
enum class geo_angle { latitude, longitude };

namespace detail {

// is s[i] the start of the degree sign, as UTF-8 (C2 B0) or as bare Latin-1 (B0)?
inline size_t degree_sign_len(std::string const& s, size_t i)
{
    auto const u = [&](size_t k) { return static_cast<unsigned char>(s[k]); };
    if (i + 1 < s.size() && u(i) == 0xC2u && u(i + 1) == 0xB0u) return 2;
    if (i < s.size() && u(i) == 0xB0u) return 1;
    return 0;
}

// is s[i] the start of the given U+20xx prime mark (E2 80 xx in UTF-8)?
inline size_t prime_len(std::string const& s, size_t i, unsigned char last)
{
    auto const u = [&](size_t k) { return static_cast<unsigned char>(s[k]); };
    if (i + 2 < s.size() && u(i) == 0xE2u && u(i + 1) == 0x80u && u(i + 2) == last) {
        return 3;
    }
    return 0;
}

inline void skip_space(std::string const& s, size_t& i)
{
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i])))
        ++i;
}

} // namespace detail

// Parse an angle given in degree/minute/second notation into decimal degrees.
//
// ACCEPTED INPUT FORMAT
//
//     [sign] number [unit] [number [unit] [number [unit]]] [hemisphere]
//
// - sign       : an optional leading '+' or '-'. Mutually exclusive with a hemisphere
//                letter -- giving both is an error, not a silent preference.
//
// - unit       : the delimiter FOLLOWING a number sets that number's unit
//                   '°' (U+00B0)          -> degrees
//                   '\'' or '′' (U+2032)  -> minutes
//                   '"'  or '″' (U+2033)  -> seconds
//                A number with NO delimiter is read as decimal degrees. The units must
//                appear in the order degrees, minutes, seconds; each at most once.
//                Minutes and seconds must lie in [0, 60).
//
// - hemisphere : an optional trailing letter, upper or lower case
//                   'N' / 'S'  for a latitude   (geo_angle::latitude)
//                   'E' / 'W'  for a longitude  (geo_angle::longitude)
//                Only these four are accepted -- deliberately NOT the German 'O' for
//                "Ost" and no other localized spelling, since 'O' reads as "Ouest"
//                (= West) in French and would silently flip the sign of the position.
//                A letter that does not fit the angle being parsed is an error.
//
// - whitespace : ignored everywhere.
//
// Accepted examples (all of them valid):
//
//     "52°31'12.0\"N"   "52°31'12\"N"   "52°31'N"   "52°31' N"   "52.52N"   "-33.8688"
//
// Rejected examples (each throws std::invalid_argument):
//
//     "13°24'36\"O"  (localized hemisphere letter)   "-52°31'N"  (sign AND hemisphere)
//     "52°31'E"      (letter does not match a latitude)
//     "52'31°N"      (units out of order)            "52°70'N"   (minutes >= 60)
//
// Degrees are returned positive towards north / east. Anything not covered above throws
// std::invalid_argument naming the offending string -- the parser never guesses.
inline value_t dms2deg(std::string const& s, geo_angle which)
{
    using namespace detail;

    size_t i = 0;
    skip_space(s, i);

    value_t sign = 1.0;
    bool signed_explicitly = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
        sign = (s[i] == '-') ? -1.0 : 1.0;
        signed_explicitly = true;
        ++i;
    }

    value_t deg = 0.0, min = 0.0, sec = 0.0;
    int unit_seen = 0; // 0: none yet, 1: degrees, 2: minutes, 3: seconds

    while (true) {
        skip_space(s, i);
        if (i >= s.size()) break;
        if (!(std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.')) break;

        char const* first = s.c_str() + i;
        char* last = nullptr;
        value_t const val = std::strtod(first, &last);
        if (last == first) {
            throw std::invalid_argument("dms2deg: malformed number in '" + s + "'");
        }
        i += static_cast<size_t>(last - first);

        skip_space(s, i);

        size_t const n_deg = degree_sign_len(s, i);
        size_t const n_min = prime_len(s, i, 0xB2u); // U+2032
        size_t const n_sec = prime_len(s, i, 0xB3u); // U+2033

        int unit = 0;
        if (n_deg > 0) {
            unit = 1;
            i += n_deg;
        }
        else if (i < s.size() && s[i] == '\'') {
            unit = 2;
            ++i;
        }
        else if (n_min > 0) {
            unit = 2;
            i += n_min;
        }
        else if (i < s.size() && s[i] == '"') {
            unit = 3;
            ++i;
        }
        else if (n_sec > 0) {
            unit = 3;
            i += n_sec;
        }
        else {
            unit = 1; // a bare number is decimal degrees
        }

        if (unit <= unit_seen) {
            throw std::invalid_argument("dms2deg: repeated or out-of-order unit in '" +
                                        s + "'");
        }
        unit_seen = unit;

        switch (unit) {
            case 1:
                deg = val;
                break;
            case 2:
                min = val;
                break;
            case 3:
                sec = val;
                break;
        }
    }

    if (unit_seen == 0) {
        throw std::invalid_argument("dms2deg: no angle found in '" + s + "'");
    }
    if (min < 0.0 || min >= 60.0 || sec < 0.0 || sec >= 60.0) {
        throw std::invalid_argument("dms2deg: minutes/seconds out of [0,60) in '" + s +
                                    "'");
    }

    skip_space(s, i);

    if (i < s.size()) {
        char const h = static_cast<char>(std::toupper(static_cast<unsigned char>(s[i])));
        ++i;
        skip_space(s, i);
        if (i != s.size()) {
            throw std::invalid_argument("dms2deg: trailing characters in '" + s + "'");
        }
        if (signed_explicitly) {
            throw std::invalid_argument("dms2deg: both a sign and a hemisphere letter"
                                        " in '" +
                                        s + "'");
        }
        bool const is_lat = (which == geo_angle::latitude);
        if (is_lat && (h == 'N' || h == 'S')) {
            sign = (h == 'S') ? -1.0 : 1.0;
        }
        else if (!is_lat && (h == 'E' || h == 'W')) {
            sign = (h == 'W') ? -1.0 : 1.0;
        }
        else {
            // fail loudly and name the accepted set: a localized letter ('O' for the
            // German "Ost", say) is rejected rather than guessed at, because the same
            // letter means the opposite hemisphere in other languages
            throw std::invalid_argument(
                std::string("dms2deg: unusable hemisphere letter '") + h + "' in '" + s +
                "' -- expected " +
                (which == geo_angle::latitude ? "'N' or 'S' for a latitude"
                                              : "'E' or 'W' for a longitude") +
                " (localized spellings are not accepted)");
        }
    }

    value_t const angle = sign * (deg + min / 60.0 + sec / 3600.0);

    if (which == geo_angle::latitude && std::abs(angle) > 90.0) {
        throw std::invalid_argument("dms2deg: latitude outside [-90,90] in '" + s + "'");
    }
    if (which == geo_angle::longitude && std::abs(angle) > 360.0) {
        throw std::invalid_argument("dms2deg: longitude outside [-360,360] in '" + s +
                                    "'");
    }

    return angle;
}

// the same angle in rad (see dms2deg for the accepted format)
inline value_t dms2rad(std::string const& s, geo_angle which)
{
    return deg2rad(dms2deg(s, which));
}

// Emit an angle in degree/minute/second notation -- the inverse of dms2deg(), so that a
// position can be written back in the form it was read in. Round-trips exactly:
// dms2deg(deg2dms(x, w), w) == x to the precision requested.
//
// sec_decimals fixes the resolution: 0 gives whole seconds (~31 m of latitude), the
// default 1 gives ~3 m, 3 gives ~3 cm. Rounding is applied to the TOTAL seconds before
// the value is split, so a value that rounds up through 60 carries into the minutes (and
// on into the degrees) instead of printing an impossible 60".
inline std::string deg2dms(value_t deg, geo_angle which, int sec_decimals = 1)
{
    if (sec_decimals < 0 || sec_decimals > 9) {
        throw std::invalid_argument("deg2dms: sec_decimals must lie in [0, 9]");
    }

    bool const negative = (deg < 0.0);

    value_t const scale = std::pow(value_t(10.0), value_t(sec_decimals));
    value_t total_sec = std::round(std::abs(deg) * 3600.0 * scale) / scale;

    int const d = static_cast<int>(total_sec / 3600.0);
    total_sec -= value_t(d) * 3600.0;
    int const m = static_cast<int>(total_sec / 60.0);
    value_t const s = total_sec - value_t(m) * 60.0;

    char const hemi =
        (which == geo_angle::latitude) ? (negative ? 'S' : 'N') : (negative ? 'W' : 'E');

    int const width = (sec_decimals > 0) ? sec_decimals + 3 : 2;

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%d°%02d'%0*.*f\"%c", d, m, width, sec_decimals, s,
                  hemi);
    return std::string(buf);
}

// the same, from radians
inline std::string rad2dms(value_t rad, geo_angle which, int sec_decimals = 1)
{
    return deg2dms(rad2deg(rad), which, sec_decimals);
}


/////////////////////////////////////////////////////////////////////////////////////////
// geodetic position
/////////////////////////////////////////////////////////////////////////////////////////

// A position as the user provides it -- copied straight out of an encyclopedia entry,
// with the two angles left in the notation they were written in:
//
//     geo_pos_dms Berlin{"52°31'12\"N", "13°24'36\"E", 35.0};
//
// Nothing is parsed here; to_geo_pos() below turns this into the calculation form.
struct geo_pos_dms {

    std::string lat;  // latitude in dms notation, e.g. "52°31'12\"N"
    std::string lon;  // longitude in dms notation, e.g. "13°24'36\"E"
    value_t height{}; // elevation [m] above sea level, as quoted
};

// The same position ready to calculate with: both angles in rad.
//
// WHAT YOU GET -- the two conventions this type commits to:
//
// - The latitude is the GEODETIC one, i.e. the angle of the ellipsoid NORMAL at the
//   position against the equatorial plane. That is the latitude an encyclopedia entry,
//   a map or a GNSS receiver quotes. It is not the angle of the ray from the geocenter
//   (the geocentric latitude), which is smaller by up to ~11.5' at mid-latitudes.
//
// - The height is measured ALONG THAT NORMAL, from the reference ellipsoid, which is
//   likewise the geodetic convention. It is therefore not simply added to a radius:
//   distance_from_geocenter() below derives the radius rather than storing it, because
//   at non-zero height the position no longer sits on the geocentric ray.
//
// NOTE on the height datum -- the one place where easy-to-look-up data is not exact.
// A quoted town elevation is an ORTHOMETRIC height H, measured against mean sea level,
// i.e. against the GEOID. The height stored here is an ELLIPSOIDAL height h, measured
// against the reference ellipsoid. They differ by the geoid undulation N:
//
//     h = H + N
//
// N is the separation of the two surfaces: a MEASURED field, not derivable from any
// formula, and -- unlike every other quantity in this header -- a function of BOTH
// latitude AND longitude. Globally it runs from about -106 m at the Indian Ocean geoid
// low, south of India, to about +85 m near New Guinea; Europe sits a few tens of metres
// positive.
//
// That longitude dependence is exactly why N is supplied rather than computed, and why
// it is not a member of ellipsoid: the reference ellipsoid is a surface of REVOLUTION,
// so everything derived from it -- radius(), N(), geocentric_lat() -- is a function of
// latitude alone, and a longitude argument would have nothing to do. Modelling the geoid
// instead means a spherical-harmonic expansion (EGM2008 runs to degree 2190), which is a
// data table, a different kind of object from the two radii this header is built on.
//
// See to_geo_pos(),
// which takes N as an optional argument and defaults it to zero: with the default, an
// encyclopedia elevation is passed through unchanged and the position is right relative
// to the ELLIPSOID, which is what keeps the type usable with data anyone can look up.
// Supply N when the vertical datum matters.
struct geo_pos {

    value_t lat{};    // geodetic latitude [rad], positive towards north
    value_t lon{};    // longitude [rad], positive towards east
    value_t height{}; // [m] above the reference ellipsoid, along the ellipsoid normal
};

// Exact equality, field by field (free, for the reason given at the ellipsoid's).
inline bool operator==(geo_pos const& a, geo_pos const& b)
{
    return a.lat == b.lat && a.lon == b.lon && a.height == b.height;
}

// distance of the position from the geocenter [m] -- the "total radius". It equals
// el.radius(lat) + height only on the ellipsoid surface (height == 0); off the surface
// the height runs along the normal, which is not the geocentric ray.
inline value_t distance_from_geocenter(geo_pos const& p, ellipsoid const& el = wgs84)
{
    value_t const N = el.N(p.lat);
    value_t const s = std::sin(p.lat);
    value_t const c = std::cos(p.lat);

    value_t const r_eq = (N + p.height) * c;                     // equatorial component
    value_t const r_ax = (N * (1.0 - el.e_sq()) + p.height) * s; // axial component

    return std::sqrt(r_eq * r_eq + r_ax * r_ax);
}

// The explicit conversion between the two, done once before the calculation.
//
// geoid_undulation is the separation N between the geoid (mean sea level, what a quoted
// elevation is measured against) and the reference ellipsoid (what geo_pos measures
// against), so the stored ellipsoidal height becomes
//
//     h = H_quoted + geoid_undulation
//
// Leaving it at zero treats the quoted elevation as if it were ellipsoidal: the position
// is then right relative to the ellipsoid and off by N vertically -- some tens of metres
// -- relative to sea level. That is the deliberate default, because N cannot be computed
// from lat/lon and looking it up is a separate step; pass it when the vertical datum
// matters. It does NOT affect latitude or longitude, and so does not affect a local
// frame's orientation -- only how far out along the normal the position sits.
inline geo_pos to_geo_pos(geo_pos_dms const& p, value_t geoid_undulation = 0.0)
{
    return geo_pos{dms2rad(p.lat, geo_angle::latitude),
                   dms2rad(p.lon, geo_angle::longitude), p.height + geoid_undulation};
}

// Intermediate positions on an edge that is straight in (lat, lon) -- the convention of
// most border data sets, whose vertices are joined by straight lines in a plate carree.
// Projected onto a plane such an edge is a curve, and drawing only its end points
// replaces it by a chord, visibly so for a long edge such as a border along a parallel.
//
// Returns the positions strictly between a and b, spaced so that neither angle changes
// by more than max_step [rad] between neighbours; the height is interpolated alike. The
// longitude difference is taken the short way round, so an edge across the 180-degree
// meridian is not walked around the globe; the longitudes returned may then run past
// +-180 deg, which geo_to_ecef accepts.
inline std::vector<geo_pos> densify(geo_pos const& a, geo_pos const& b, value_t max_step)
{
    if (!(max_step > 0.0)) {
        throw std::invalid_argument("densify: max_step must be positive");
    }

    value_t dlon = b.lon - a.lon;
    if (dlon > pi) dlon -= 2.0 * pi;
    else if (dlon < -pi) dlon += 2.0 * pi;
    value_t const dlat = b.lat - a.lat;
    value_t const dh = b.height - a.height;

    value_t const span = std::max(std::abs(dlat), std::abs(dlon));
    auto const n = static_cast<size_t>(std::ceil(span / max_step)); // segments

    std::vector<geo_pos> out;
    if (n > 1) {
        out.reserve(n - 1);
        for (size_t i = 1; i < n; ++i) {
            value_t const f = static_cast<value_t>(i) / static_cast<value_t>(n);
            out.push_back(geo_pos{a.lat + f * dlat, a.lon + f * dlon, a.height + f * dh});
        }
    }
    return out;
}


/////////////////////////////////////////////////////////////////////////////////////////
// transverse Mercator: geodetic position <-> grid coordinates (UTM, national grids)
/////////////////////////////////////////////////////////////////////////////////////////
//
// The conformal map onto a cylinder that touches the ellipsoid along a central meridian:
// the basis of UTM and of most national grids. Official geodata is often delivered only
// in such GRID coordinates -- easting and northing in metres -- so reading it needs the
// inverse, tm_inverse(), and writing it the forward map, tm_forward().
//
// This is classical numerics with no geometric-algebra form, and is written as such:
// Krueger's series in the third flattening n = f / (2 - f), to order n^6, in the form
// given by Karney ("Transverse Mercator with an accuracy of a few nanometers", J. Geodesy
// 85, 2011), on the EXACT conformal latitude -- the inverse solves for the geodetic
// latitude by Newton's method instead of a second series. Within a UTM zone (3 deg either
// side of the central meridian) it agrees with PROJ to well below a millimetre (tested);
// like every series form it degrades far outside the zone.
//
// Grid coordinates are returned as a vec2d: x the easting, y the northing [m].

// the parameters of one transverse Mercator grid
struct transverse_mercator {
    value_t lat0{0.0};           // latitude of the origin [rad]
    value_t lon0{0.0};           // the central meridian [rad]
    value_t k0{1.0};             // the scale on the central meridian
    value_t false_easting{0.0};  // [m]
    value_t false_northing{0.0}; // [m]
};

// UTM zone 1..60: the central meridian at 6 zone - 183 deg, k0 = 0.9996, a false easting
// of 500 km, and a false northing of 0 (northern) or 10 000 km (southern hemisphere)
inline transverse_mercator utm_zone(int zone, bool north = true)
{
    if (zone < 1 || zone > 60) {
        throw std::invalid_argument("utm_zone: the zone must be in 1..60");
    }
    return transverse_mercator{0.0, deg2rad(6.0 * zone - 183.0), 0.9996, 500000.0,
                               north ? 0.0 : 10000000.0};
}

namespace detail {

// the coefficients of Krueger's series for one ellipsoid (Karney 2011, eqs. 14 and 15)
struct tm_series {
    value_t e;                    // first eccentricity
    value_t A;                    // the rectifying radius
    std::array<value_t, 6> alpha; // forward: conformal sphere -> rectifying
    std::array<value_t, 6> beta;  // inverse: rectifying -> conformal sphere
};

inline tm_series tm_coefficients(ellipsoid const& el)
{
    value_t const f = el.flattening();
    value_t const n = f / (2.0 - f);
    value_t const n2 = n * n, n3 = n2 * n, n4 = n3 * n, n5 = n4 * n, n6 = n5 * n;

    tm_series s{};
    s.e = std::sqrt(el.e_sq());
    s.A = el.r_equator / (1.0 + n) * (1.0 + n2 / 4.0 + n4 / 64.0 + n6 / 256.0);
    s.alpha = {n / 2.0 - 2.0 * n2 / 3.0 + 5.0 * n3 / 16.0 + 41.0 * n4 / 180.0 -
                   127.0 * n5 / 288.0 + 7891.0 * n6 / 37800.0,
               13.0 * n2 / 48.0 - 3.0 * n3 / 5.0 + 557.0 * n4 / 1440.0 +
                   281.0 * n5 / 630.0 - 1983433.0 * n6 / 1935360.0,
               61.0 * n3 / 240.0 - 103.0 * n4 / 140.0 + 15061.0 * n5 / 26880.0 +
                   167603.0 * n6 / 181440.0,
               49561.0 * n4 / 161280.0 - 179.0 * n5 / 168.0 + 6601661.0 * n6 / 7257600.0,
               34729.0 * n5 / 80640.0 - 3418889.0 * n6 / 1995840.0,
               212378941.0 * n6 / 319334400.0};
    s.beta = {n / 2.0 - 2.0 * n2 / 3.0 + 37.0 * n3 / 96.0 - n4 / 360.0 -
                  81.0 * n5 / 512.0 + 96199.0 * n6 / 604800.0,
              n2 / 48.0 + n3 / 15.0 - 437.0 * n4 / 1440.0 + 46.0 * n5 / 105.0 -
                  1118711.0 * n6 / 3870720.0,
              17.0 * n3 / 480.0 - 37.0 * n4 / 840.0 - 209.0 * n5 / 4480.0 +
                  5569.0 * n6 / 90720.0,
              4397.0 * n4 / 161280.0 - 11.0 * n5 / 504.0 - 830251.0 * n6 / 7257600.0,
              4583.0 * n5 / 161280.0 - 108847.0 * n6 / 3991680.0,
              20648693.0 * n6 / 638668800.0};
    return s;
}

// tan of the conformal latitude, from tan of the geodetic one:
// tau' = sinh(asinh(tau) - e atanh(e sin(lat)))
inline value_t tm_taup(value_t tau, value_t e)
{
    value_t const s = tau / std::sqrt(1.0 + tau * tau); // sin(lat)
    return std::sinh(std::asinh(tau) - e * std::atanh(e * s));
}

// the inverse of tm_taup, by Newton's method on
// d tau' / d tau = (1 - e^2) sqrt(1 + tau'^2) sqrt(1 + tau^2) / (1 + (1 - e^2) tau^2)
inline value_t tm_tau(value_t taup, value_t e)
{
    value_t const e2m = 1.0 - e * e;
    value_t tau = taup / e2m; // close enough to start from everywhere
    for (int i = 0; i < 10; ++i) {
        value_t const tp = tm_taup(tau, e);
        value_t const d = (taup - tp) * (1.0 + e2m * tau * tau) /
                          (e2m * std::sqrt(1.0 + tp * tp) * std::sqrt(1.0 + tau * tau));
        tau += d;
        if (std::abs(d) <= 1.0e-15 * std::max(value_t(1.0), std::abs(tau))) break;
    }
    return tau;
}

// the rectifying-sphere coordinate xi of the central meridian at latitude lat
inline value_t tm_xi0(tm_series const& s, value_t lat)
{
    value_t const xip = std::atan(tm_taup(std::tan(lat), s.e));
    value_t xi = xip;
    for (int j = 1; j <= 6; ++j)
        xi += s.alpha[size_t(j - 1)] * std::sin(2.0 * j * xip);
    return xi;
}

} // namespace detail

// geodetic position -> grid coordinates (x easting, y northing) [m]; the height is
// ignored, as in every map projection
inline vec2d tm_forward(geo_pos const& p, transverse_mercator const& tm,
                        ellipsoid const& el = wgs84)
{
    auto const s = detail::tm_coefficients(el);

    value_t lam = p.lon - tm.lon0; // the longitude from the central meridian, (-pi, pi]
    lam = std::remainder(lam, 2.0 * pi);

    value_t const taup = detail::tm_taup(std::tan(p.lat), s.e);
    value_t const xip = std::atan2(taup, std::cos(lam));
    value_t const etap = std::asinh(std::sin(lam) / std::hypot(taup, std::cos(lam)));

    value_t xi = xip, eta = etap;
    for (int j = 1; j <= 6; ++j) {
        value_t const a = s.alpha[size_t(j - 1)];
        xi += a * std::sin(2.0 * j * xip) * std::cosh(2.0 * j * etap);
        eta += a * std::cos(2.0 * j * xip) * std::sinh(2.0 * j * etap);
    }

    return vec2d{tm.false_easting + tm.k0 * s.A * eta,
                 tm.false_northing + tm.k0 * s.A * (xi - detail::tm_xi0(s, tm.lat0))};
}

// grid coordinates (x easting, y northing) [m] -> geodetic position, height 0
inline geo_pos tm_inverse(vec2d const& en, transverse_mercator const& tm,
                          ellipsoid const& el = wgs84)
{
    auto const s = detail::tm_coefficients(el);

    value_t const xi =
        (en.y - tm.false_northing) / (tm.k0 * s.A) + detail::tm_xi0(s, tm.lat0);
    value_t const eta = (en.x - tm.false_easting) / (tm.k0 * s.A);

    value_t xip = xi, etap = eta;
    for (int j = 1; j <= 6; ++j) {
        value_t const b = s.beta[size_t(j - 1)];
        xip -= b * std::sin(2.0 * j * xi) * std::cosh(2.0 * j * eta);
        etap -= b * std::cos(2.0 * j * xi) * std::sinh(2.0 * j * eta);
    }

    value_t const taup = std::sin(xip) / std::hypot(std::sinh(etap), std::cos(xip));
    value_t const lam = std::atan2(std::sinh(etap), std::cos(xip));

    return geo_pos{std::atan(detail::tm_tau(taup, s.e)), tm.lon0 + lam, 0.0};
}


/////////////////////////////////////////////////////////////////////////////////////////
// Lambert azimuthal equal-area: geodetic position <-> map coordinates, area kept
/////////////////////////////////////////////////////////////////////////////////////////
//
// The map that keeps AREA: a region has the same area on the map as on the ellipsoid, at
// every distance from the centre and in every direction. That is what a map is for when
// it compares quantities per area -- population density, land use, statistics -- and why
// the EU's statistical grid is this projection: ETRS89-LAEA, EPSG:3035, centred at 52 deg
// N 10 deg E (etrs89_laea()). The price is shape: away from the centre a region is
// squeezed along the radius and stretched across it by the same factor, so its area stays
// and its outline does not.
//
// Unlike the projections onto the tangent plane in hd::ga::pga, this map is NOT a central
// projection: no centre exists from which the lines through the points meet a plane at
// the map positions, so the meet X = rwdg(wdg(C, P), plane) has nothing to offer here. It
// is classical numerics, written as such (EPSG method 9820; Snyder, "Map Projections -- A
// Working Manual", USGS 1987, section 24). The ellipsoid is replaced by the AUTHALIC
// sphere -- the sphere of the same surface area, radius Rq -- with each geodetic latitude
// moved to the authalic latitude beta, which divides the ellipsoid's area between the
// equator and the latitude in the same ratio as the sphere's. On that sphere the
// equal-area azimuthal map has a closed form, and a factor D restores the true scale at
// the centre. The inverse returns from beta to the geodetic latitude by Newton's method
// on the authalic function q, not by the usual truncated series. Agreement with PROJ is
// a few micrometres (tested).
//
// The whole ellipsoid is mapped into a disc of radius 2 Rq, except the ANTIPODE of the
// centre, which has no image (it would be the whole rim). Map coordinates are returned
// as a vec2d: x towards east, y towards north at the centre [m]. NOTE: EPSG lists the
// axes of EPSG:3035 as (northing, easting), so tools that honour EPSG's axis order print
// y before x.

// the parameters of one Lambert azimuthal equal-area map
struct lambert_azimuthal {
    value_t lat0{0.0};           // latitude of the centre [rad]
    value_t lon0{0.0};           // longitude of the centre [rad]
    value_t false_easting{0.0};  // [m]
    value_t false_northing{0.0}; // [m]
};

// the map centred at a position, the centre at map coordinates (0, 0)
inline lambert_azimuthal laea_at(geo_pos const& centre)
{
    return lambert_azimuthal{centre.lat, centre.lon, 0.0, 0.0};
}

// ETRS89-LAEA (EPSG:3035), the EU's equal-area grid: centred at 52 deg N 10 deg E, false
// easting 4321 km, false northing 3210 km. Defined on GRS80 -- pass grs80 to the
// functions below; WGS84 differs from it by about 0.1 mm in the pole radius
inline lambert_azimuthal etrs89_laea()
{
    return lambert_azimuthal{deg2rad(52.0), deg2rad(10.0), 4321000.0, 3210000.0};
}

namespace detail {

// atanh(e s) / e, which is s on a sphere
inline value_t laea_atanh_e(value_t s, value_t e)
{
    return (e == 0.0) ? s : std::atanh(e * s) / e;
}

// the authalic function (Snyder eq. 3-12): the area of the ellipsoid between the equator
// and the latitude whose sine is s is pi a^2 q(s)
inline value_t laea_q(value_t s, value_t e)
{
    value_t const e2 = e * e;
    return (1.0 - e2) * (s / (1.0 - e2 * s * s) + laea_atanh_e(s, e));
}

// the inverse of laea_q: the geodetic latitude for a value of q, by Newton's method
// (Snyder eq. 3-16, dq/dlat = 2 (1 - e^2) cos(lat) / (1 - e^2 sin^2(lat))^2)
inline value_t laea_lat_from_q(value_t q, value_t qp, value_t e)
{
    if (std::abs(q) >= qp) return std::copysign(0.5 * pi, q);
    value_t const e2 = e * e;
    value_t lat = std::asin(std::clamp(q / 2.0, value_t(-1.0), value_t(1.0)));
    for (int i = 0; i < 20; ++i) {
        value_t const s = std::sin(lat), c = std::cos(lat);
        if (c < 1.0e-12) break; // at a pole, where q already is +-qp
        value_t const w = 1.0 - e2 * s * s;
        value_t const d =
            w * w / (2.0 * c) * (q / (1.0 - e2) - s / w - laea_atanh_e(s, e));
        lat += d;
        if (std::abs(d) <= 1.0e-15) break;
    }
    return lat;
}

// what one map needs from its ellipsoid and centre
struct laea_constants {
    value_t e;      // first eccentricity
    value_t qp;     // q at the pole
    value_t Rq;     // radius of the authalic sphere
    value_t sin_b0; // sine and cosine of the centre's authalic latitude
    value_t cos_b0;
    value_t D; // restores the true scale at the centre (1 at a pole)
};

inline laea_constants laea_setup(lambert_azimuthal const& m, ellipsoid const& el)
{
    laea_constants c{};
    c.e = std::sqrt(el.e_sq());
    c.qp = laea_q(1.0, c.e);
    c.Rq = el.r_equator * std::sqrt(0.5 * c.qp);
    value_t const s0 = std::sin(m.lat0);
    c.sin_b0 = std::clamp(laea_q(s0, c.e) / c.qp, value_t(-1.0), value_t(1.0));
    c.cos_b0 = std::sqrt(std::max(value_t(0.0), 1.0 - c.sin_b0 * c.sin_b0));
    // D = a m0 / (Rq cos(beta0)), m0 = cos(lat0) / sqrt(1 - e^2 sin^2(lat0)); its limit
    // at a pole is 1, where the oblique form below becomes the polar one
    c.D = (c.cos_b0 < 1.0e-12)
              ? 1.0
              : el.r_equator * std::cos(m.lat0) / std::sqrt(1.0 - el.e_sq() * s0 * s0) /
                    (c.Rq * c.cos_b0);
    return c;
}

} // namespace detail

// geodetic position -> map coordinates (x towards east, y towards north) [m]; the
// height is ignored, as in every map projection. Throws for the antipode of the centre.
inline vec2d laea_forward(geo_pos const& p, lambert_azimuthal const& m,
                          ellipsoid const& el = wgs84)
{
    auto const c = detail::laea_setup(m, el);

    value_t const sin_b = std::clamp(detail::laea_q(std::sin(p.lat), c.e) / c.qp,
                                     value_t(-1.0), value_t(1.0));
    value_t const cos_b = std::sqrt(std::max(value_t(0.0), 1.0 - sin_b * sin_b));
    value_t const dlon = std::remainder(p.lon - m.lon0, 2.0 * pi);

    // 1 + cos of the angle from the centre on the authalic sphere
    value_t const den = 1.0 + c.sin_b0 * sin_b + c.cos_b0 * cos_b * std::cos(dlon);
    if (!(den > 1.0e-14)) {
        throw std::invalid_argument("laea_forward: the antipode of the centre has no "
                                    "image");
    }
    value_t const B = c.Rq * std::sqrt(2.0 / den);

    return vec2d{m.false_easting + B * c.D * cos_b * std::sin(dlon),
                 m.false_northing +
                     B / c.D * (c.cos_b0 * sin_b - c.sin_b0 * cos_b * std::cos(dlon))};
}

// map coordinates (x towards east, y towards north) [m] -> geodetic position, height 0.
// Throws for a point outside the disc of radius 2 Rq the map fills.
inline geo_pos laea_inverse(vec2d const& xy, lambert_azimuthal const& m,
                            ellipsoid const& el = wgs84)
{
    auto const c = detail::laea_setup(m, el);

    value_t const x = xy.x - m.false_easting;
    value_t const y = xy.y - m.false_northing;
    value_t const rho = std::hypot(x / c.D, c.D * y);
    if (rho == 0.0) return geo_pos{m.lat0, m.lon0, 0.0};
    if (rho > 2.0 * c.Rq * (1.0 + 1.0e-12)) {
        throw std::invalid_argument("laea_inverse: the point lies outside the map");
    }

    // the angle from the centre on the authalic sphere, and the authalic latitude there
    value_t const C = 2.0 * std::asin(std::min(value_t(1.0), rho / (2.0 * c.Rq)));
    value_t const sC = std::sin(C), cC = std::cos(C);
    value_t const sin_b = std::clamp(cC * c.sin_b0 + c.D * y * sC * c.cos_b0 / rho,
                                     value_t(-1.0), value_t(1.0));
    value_t const lon = m.lon0 + std::atan2(x * sC, c.D * rho * c.cos_b0 * cC -
                                                        c.D * c.D * y * c.sin_b0 * sC);

    return geo_pos{detail::laea_lat_from_q(c.qp * sin_b, c.qp, c.e), lon, 0.0};
}


/////////////////////////////////////////////////////////////////////////////////////////
// geodesics on the ellipsoid: the inverse and the direct problem
/////////////////////////////////////////////////////////////////////////////////////////
//
// The shortest path on the ellipsoid between two positions -- what a long-haul flight
// follows, and what a "great circle" means once the earth is not a sphere. Two problems:
//
// - INVERSE: given two positions, the distance s12 between them and the azimuth of the
//   path at each end (azi1 departing, azi2 arriving; clockwise from north)
// - DIRECT:  given a position, an azimuth and a distance, where the path ends
//
// The path is not a plane section of the ellipsoid, so no product expression of the
// algebra yields it; this is classical numerics, written as such. It follows Bessel's
// method in the form of Karney ("Algorithms for geodesics", J. Geodesy 87, 2013): the
// geodesic is mapped onto an auxiliary sphere through the reduced latitude beta,
// tan(beta) = (1 - f) tan(lat), where it is a great circle with the arc length sigma
// and the equatorial azimuth alpha0 (sin(alpha0) = sin(alpha) cos(beta), Clairaut's
// constant), and distance and longitude follow from two integrals over sigma,
//
//     s      = b  int sqrt(1 + k^2 sin^2(sigma)) dsigma ,        k^2 = e'^2 cos^2(alpha0)
//     lambda = omega - f sin(alpha0) int (2 - f) / (1 + (1 - f) sqrt(1 + k^2 sin^2))
//     dsigma
//
// (Karney's eqs. 7 and 8; omega the longitude on the auxiliary sphere). Where Karney
// expands both integrals in series, they are evaluated here by Gauss-Legendre
// quadrature -- the integrands are smooth and nearly constant (k^2 < 0.007 on the earth),
// so twelve nodes per quarter circle reach machine precision. The inverse problem is
// solved by bracketing the departure azimuth: with the two ends arranged canonically
// (the first the one farther from the equator, in the southern hemisphere, the second
// east of it), the longitude difference reached grows monotonically with the azimuth
// from 0 (due north) to pi (due south, across the pole) -- Karney's section 4 -- so a
// regula falsi on [0, pi] converges for every pair, nearly antipodal ones included,
// where Vincenty's iteration fails. Agreement with GeographicLib (PROJ's geod) is a
// fraction of a micrometre over the whole earth (tested).
//
// Degenerate cases, handled as Karney does: an end at a pole, where the azimuth is the
// meridian's longitude; both ends on the equator, where the equator itself is the
// shortest path up to (1 - f) pi of longitude, and beyond that a path that leaves it and
// meets it again half a circle on (the northern of the two mirror images).

// the solution of the inverse problem
struct geodesic_solution {
    value_t s12{0.0};  // distance along the ellipsoid [m]
    value_t azi1{0.0}; // azimuth at the start, departing [rad], clockwise from north
    value_t azi2{0.0}; // azimuth at the end, arriving (the direction of travel there)
};

// the solution of the direct problem
struct geodesic_destination {
    geo_pos pos{};     // where the path ends (height 0)
    value_t azi2{0.0}; // azimuth there, arriving [rad]
};

namespace detail {

// the 12-point Gauss-Legendre rule: gauss_legendre12() in detail/ga_quadrature.hpp

// int_a^b g(sigma) dsigma, twelve nodes per panel of at most a quarter circle
template <typename F> inline value_t geod_quadrature(F const& g, value_t a, value_t b)
{
    auto const& q = gauss_legendre12();
    int const n = std::max(1, static_cast<int>(std::ceil(std::abs(b - a) / (0.5 * pi))));
    value_t const h = (b - a) / n;
    value_t sum = 0.0;
    for (int p = 0; p < n; ++p) {
        value_t const c = a + (p + 0.5) * h;
        for (size_t i = 0; i < 12; ++i)
            sum += q.w[i] * g(c + 0.5 * h * q.x[i]);
    }
    return 0.5 * h * sum;
}

// one geodesic on the auxiliary sphere: the two integrals of Karney's eqs. 7 and 8
struct geod_line {
    value_t f;      // flattening
    value_t sin_a0; // Clairaut's constant sin(alpha0)
    value_t k2;     // e'^2 cos^2(alpha0)

    value_t g1(value_t s) const
    {
        value_t const sn = std::sin(s);
        return std::sqrt(1.0 + k2 * sn * sn);
    }
    // the distance integral, in units of the pole radius
    value_t i1(value_t a, value_t b) const
    {
        return geod_quadrature([this](value_t s) { return g1(s); }, a, b);
    }
    // the longitude correction lambda - omega
    value_t dlam(value_t a, value_t b) const
    {
        if (sin_a0 == 0.0) return 0.0;
        return -f * sin_a0 *
               geod_quadrature(
                   [this](value_t s) { return (2.0 - f) / (1.0 + (1.0 - f) * g1(s)); }, a,
                   b);
    }
};

inline geod_line geod_line_from(ellipsoid const& el, value_t sin_a0)
{
    value_t const e2 = el.e_sq();
    value_t const ep2 = e2 / (1.0 - e2);
    return geod_line{el.flattening(), sin_a0, ep2 * (1.0 - sin_a0 * sin_a0)};
}

// the reduced latitude, as its sine and cosine
inline void reduced_latitude(ellipsoid const& el, value_t lat, value_t& sb, value_t& cb)
{
    value_t const f = el.flattening();
    value_t const y = (1.0 - f) * std::sin(lat), x = std::cos(lat);
    value_t const r = std::hypot(x, y);
    sb = y / r;
    cb = x / r;
}

// the canonical inverse problem: beta1 <= 0, |beta2| <= |beta1|, 0 <= lam12 <= pi;
// the departure azimuth in [0, pi], the path followed to its first northward crossing of
// beta2. Returns lam12 reached for the azimuth a1, and the line's quantities.
struct geod_shot {
    value_t lam12;
    value_t sig1, sig2;
    value_t sin_a0;
    value_t ca2cb2; // cos(alpha2) cos(beta2) at the end
};

inline geod_shot geod_shoot(ellipsoid const& el, value_t sb1, value_t cb1, value_t sb2,
                            value_t cb2, value_t a1)
{
    value_t const sa1 = std::sin(a1), ca1 = std::cos(a1);
    value_t const sin_a0 = sa1 * cb1;
    value_t const sig1 = std::atan2(sb1, ca1 * cb1);
    // Karney eq. 45: cos(alpha2) cos(beta2), northward, written to avoid cancellation
    value_t const ca2cb2 = std::sqrt(
        std::max(value_t(0.0), ca1 * ca1 * cb1 * cb1 + (cb2 - cb1) * (cb2 + cb1)));
    value_t const sig2 = std::atan2(sb2, ca2cb2);
    // omega on the auxiliary sphere: omega1 in (-pi, 0], omega2 in [-pi/2, pi/2], so the
    // difference is continuous in a1 and needs no unwrapping
    value_t const om1 = std::atan2(sin_a0 * std::sin(sig1), std::cos(sig1));
    value_t const om2 = std::atan2(sin_a0 * std::sin(sig2), std::cos(sig2));
    auto const L = geod_line_from(el, sin_a0);
    return geod_shot{om2 - om1 + L.dlam(sig1, sig2), sig1, sig2, sin_a0, ca2cb2};
}

// the root of an increasing function on [lo, hi] (F(lo) <= 0 <= F(hi)), by regula
// falsi with the Illinois modification: superlinear, and never leaves the bracket
template <typename F> inline value_t geod_root(F const& fn, value_t lo, value_t hi)
{
    value_t f_lo = fn(lo), f_hi = fn(hi);
    if (f_lo >= 0.0) return lo;
    if (f_hi <= 0.0) return hi;
    value_t x = lo;
    int side = 0;
    for (int it = 0; it < 200; ++it) {
        x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);
        value_t const fx = fn(x);
        if (fx == 0.0 || hi - lo < 1.0e-15 || std::abs(fx) < 1.0e-15) break;
        if (fx < 0.0) {
            lo = x;
            f_lo = fx;
            if (side == -1) f_hi *= 0.5;
            side = -1;
        }
        else {
            hi = x;
            f_hi = fx;
            if (side == 1) f_lo *= 0.5;
            side = 1;
        }
    }
    return x;
}

} // namespace detail

// the inverse problem: distance and azimuths of the shortest path from a to b (heights
// ignored)
inline geodesic_solution geodesic_inverse(geo_pos const& a, geo_pos const& b,
                                          ellipsoid const& el = wgs84)
{
    value_t lat1 = a.lat, lat2 = b.lat;
    value_t lam12 = std::remainder(b.lon - a.lon, 2.0 * pi);

    // canonical arrangement: the end farther from the equator first, in the south, the
    // other east of it
    bool const swapped = std::abs(lat1) < std::abs(lat2);
    if (swapped) {
        std::swap(lat1, lat2);
        lam12 = -lam12;
    }
    value_t const lat_sign = (lat1 > 0.0) ? -1.0 : 1.0;
    lat1 *= lat_sign;
    lat2 *= lat_sign;
    value_t const lon_sign = (lam12 < 0.0) ? -1.0 : 1.0;
    lam12 *= lon_sign;

    value_t sb1, cb1, sb2, cb2;
    detail::reduced_latitude(el, lat1, sb1, cb1);
    detail::reduced_latitude(el, lat2, sb2, cb2);

    value_t a1 = 0.0, a2 = 0.0, s12 = 0.0;
    value_t const f = el.flattening();
    if (cb1 < 1.0e-15) {
        // the first end at the (south) pole: the path is the meridian of the other end,
        // the azimuth at the pole that meridian's longitude
        auto const L = detail::geod_line_from(el, 0.0);
        s12 = el.r_pole * L.i1(-0.5 * pi, std::atan2(sb2, cb2));
        a1 = lam12;
        a2 = 0.0;
    }
    else if (sb1 == 0.0 && lam12 <= (1.0 - f) * pi) {
        // both on the equator, close enough for the equator to be the shortest path
        s12 = el.r_equator * lam12;
        a1 = a2 = 0.5 * pi;
    }
    else if (sb1 == 0.0) {
        // both on the equator, farther apart than (1 - f) pi: the shortest path leaves
        // the equator northward and meets it again at sigma = pi, where
        // lam12(a1) = pi - f sin(a1) int_0^pi ... falls from pi (a1 = 0, over the pole)
        // to (1 - f) pi (a1 = pi/2, along the equator). Of the two mirror images the
        // northern one is taken, as GeographicLib does
        auto const lam_eq = [&](value_t az) {
            return pi + detail::geod_line_from(el, std::sin(az)).dlam(0.0, pi);
        };
        a1 = detail::geod_root([&](value_t az) { return lam12 - lam_eq(az); }, 0.0,
                               0.5 * pi);
        s12 = el.r_pole * detail::geod_line_from(el, std::sin(a1)).i1(0.0, pi);
        a2 = pi - a1;
    }
    else {
        // the general case: lam12(a1) increases monotonically over [0, pi]
        a1 = detail::geod_root(
            [&](value_t az) {
                return detail::geod_shoot(el, sb1, cb1, sb2, cb2, az).lam12 - lam12;
            },
            0.0, pi);
        auto const shot = detail::geod_shoot(el, sb1, cb1, sb2, cb2, a1);
        auto const L = detail::geod_line_from(el, shot.sin_a0);
        s12 = el.r_pole * L.i1(shot.sig1, shot.sig2);
        a2 = std::atan2(shot.sin_a0, shot.ca2cb2);
    }

    // undo the arrangement: the mirrors first, then the exchange of the ends
    if (lon_sign < 0.0) {
        a1 = -a1;
        a2 = -a2;
    }
    if (lat_sign < 0.0) {
        a1 = pi - a1;
        a2 = pi - a2;
    }
    if (swapped) {
        value_t const t = a1;
        a1 = a2 + pi;
        a2 = t + pi;
    }
    return geodesic_solution{s12, std::remainder(a1, 2.0 * pi),
                             std::remainder(a2, 2.0 * pi)};
}

// the direct problem: where the path from a, leaving at azimuth azi1 [rad, clockwise from
// north], ends after the distance s12 [m]
inline geodesic_destination geodesic_direct(geo_pos const& a, value_t azi1, value_t s12,
                                            ellipsoid const& el = wgs84)
{
    value_t sb1, cb1;
    detail::reduced_latitude(el, a.lat, sb1, cb1);
    value_t const sa1 = std::sin(azi1), ca1 = std::cos(azi1);
    value_t sin_a0 = sa1 * cb1;
    value_t cos_a0 = std::hypot(sb1, ca1 * cb1);
    value_t sig1 = std::atan2(sb1, ca1 * cb1);
    value_t lon1 = a.lon;
    if (cb1 < 1.0e-15) {
        // a start at a pole: the azimuth names the meridian taken -- its longitude is
        // lon1 + azi1 from the south pole, lon1 + pi - azi1 from the north pole (the
        // convention geodesic_inverse returns). On the auxiliary sphere the north pole
        // is sigma = pi/2, past which omega = pi: so there the base meridian is the
        // opposite one, lon1 - azi1
        sin_a0 = 0.0;
        cos_a0 = 1.0;
        sig1 = std::copysign(0.5 * pi, sb1);
        lon1 += (sb1 < 0.0) ? azi1 : -azi1;
    }
    auto const L = detail::geod_line_from(el, sin_a0);

    // the arc sigma2 with b * I1(sigma1, sigma2) = s12, by Newton (I1 increases)
    value_t sig2 = sig1 + s12 / el.r_pole;
    for (int it = 0; it < 30; ++it) {
        value_t const r = el.r_pole * L.i1(sig1, sig2) - s12;
        sig2 -= r / (el.r_pole * L.g1(sig2));
        if (std::abs(r) < 1.0e-9) break;
    }

    value_t const ss2 = std::sin(sig2), cs2 = std::cos(sig2);
    value_t const lat2 = std::atan2(cos_a0 * ss2, (1.0 - el.flattening()) *
                                                      std::hypot(sin_a0, cos_a0 * cs2));
    // omega12 matters only modulo 2 pi here: the longitude is wrapped at the end
    value_t const om1 = std::atan2(sin_a0 * std::sin(sig1), std::cos(sig1));
    value_t const om2 = std::atan2(sin_a0 * ss2, cs2);
    value_t const lon2 =
        std::remainder(lon1 + (om2 - om1) + L.dlam(sig1, sig2), 2.0 * pi);
    return geodesic_destination{geo_pos{lat2, lon2, 0.0},
                                std::atan2(sin_a0, cos_a0 * cs2)};
}

// the shortest path from a to b as positions, both ends included, spaced at most
// max_step [m] apart -- to draw a route
inline std::vector<geo_pos> geodesic_path(geo_pos const& a, geo_pos const& b,
                                          value_t max_step, ellipsoid const& el = wgs84)
{
    if (!(max_step > 0.0)) {
        throw std::invalid_argument("geodesic_path: max_step must be positive");
    }
    auto const g = geodesic_inverse(a, b, el);
    auto const n = std::max(size_t(1), static_cast<size_t>(std::ceil(g.s12 / max_step)));
    std::vector<geo_pos> out;
    out.reserve(n + 1);
    out.push_back(geo_pos{a.lat, a.lon, 0.0});
    for (size_t i = 1; i < n; ++i) {
        out.push_back(
            geodesic_direct(a, g.azi1, g.s12 * static_cast<value_t>(i) / n, el).pos);
    }
    out.push_back(geo_pos{b.lat, b.lon, 0.0});
    return out;
}


/////////////////////////////////////////////////////////////////////////////////////////
// azimuthal equidistant: true distance and direction from one point
/////////////////////////////////////////////////////////////////////////////////////////
//
// The map that keeps DISTANCE and DIRECTION from its centre: every point is drawn at its
// geodesic distance from the centre, in the direction of the geodesic's azimuth there,
//
//     x = s12 sin(azi1) ,   y = s12 cos(azi1)
//
// So the shortest route from the centre to anywhere is a STRAIGHT line on the map, of
// its true length, and the circles around the centre are the places at one distance. It
// is the map of everything seen from one place: a flight map from one airport (the
// emblem of the United Nations is its polar aspect), radio and seismic ranges, the reach
// of a missile. Distances and directions between two OTHER points are not kept, and
// shape and area grow worse with the distance; the whole earth fits into a disc of
// radius about 20 000 km, its rim the antipode of the centre.
//
// Like the equal-area map it is no central projection, so the tangent-plane meet does not
// produce it; it is the inverse geodesic problem in polar coordinates. On the ellipsoid
// this is the rigorous form (PROJ's aeqd, GeographicLib's AzimuthalEquidistant), not the
// series approximation of EPSG method 1125 / 9832, which is valid only near the centre.

// the parameters of one azimuthal equidistant map
struct azimuthal_equidistant {
    value_t lat0{0.0};           // latitude of the centre [rad]
    value_t lon0{0.0};           // longitude of the centre [rad]
    value_t false_easting{0.0};  // [m]
    value_t false_northing{0.0}; // [m]
};

// the map centred at a position, the centre at map coordinates (0, 0)
inline azimuthal_equidistant aeqd_at(geo_pos const& centre)
{
    return azimuthal_equidistant{centre.lat, centre.lon, 0.0, 0.0};
}

// geodetic position -> map coordinates (x towards east, y towards north at the centre)
inline vec2d aeqd_forward(geo_pos const& p, azimuthal_equidistant const& m,
                          ellipsoid const& el = wgs84)
{
    auto const g = geodesic_inverse(geo_pos{m.lat0, m.lon0, 0.0}, p, el);
    return vec2d{m.false_easting + g.s12 * std::sin(g.azi1),
                 m.false_northing + g.s12 * std::cos(g.azi1)};
}

// map coordinates -> geodetic position, height 0 (a point beyond the antipode's distance
// continues along its geodesic)
inline geo_pos aeqd_inverse(vec2d const& xy, azimuthal_equidistant const& m,
                            ellipsoid const& el = wgs84)
{
    value_t const x = xy.x - m.false_easting, y = xy.y - m.false_northing;
    value_t const s = std::hypot(x, y);
    if (s == 0.0) return geo_pos{m.lat0, m.lon0, 0.0};
    return geodesic_direct(geo_pos{m.lat0, m.lon0, 0.0}, std::atan2(x, y), s, el).pos;
}


/////////////////////////////////////////////////////////////////////////////////////////
// the two-dimensional case: a position in the meridian section
/////////////////////////////////////////////////////////////////////////////////////////
//
// Same ellipse, same geodetic latitude, same height along the normal -- only the
// longitude is gone, because a meridian section has no second angle. The position types
// below are therefore the three-dimensional ones minus one field, and they feed the same
// detail::meridian_* core.
//
// Coverage: a geodetic latitude runs over [-90, +90] deg, which sweeps the HALF ellipse
// with r >= 0 -- one meridian. The other half (r < 0) is the 180-degree meridian, which
// in three dimensions is reached by the longitude and here has no representation.

// a position as the user provides it: one angle in dms notation plus an elevation
struct geo_pos_dms2dp {

    std::string lat;  // latitude in dms notation, e.g. "52°31'12\"N"
    value_t height{}; // elevation [m] above sea level, as quoted
};

// the same position ready to calculate with (see geo_pos for what the two conventions
// commit to -- they are identical here, longitude aside)
struct geo_pos2dp {

    value_t lat{};    // geodetic latitude [rad], positive towards north
    value_t height{}; // [m] above the reference ellipsoid, along the ellipsoid normal
};

// the explicit conversion, with the same optional geoid undulation as its 3D counterpart
inline geo_pos2dp to_geo_pos(geo_pos_dms2dp const& p, value_t geoid_undulation = 0.0)
{
    return geo_pos2dp{dms2rad(p.lat, geo_angle::latitude), p.height + geoid_undulation};
}

// distance of the position from the geocenter [m] -- see the 3D overload
inline value_t distance_from_geocenter(geo_pos2dp const& p, ellipsoid const& el = wgs84)
{
    auto const m = detail::meridian_from_geodetic(el, p.lat, p.height);
    return std::sqrt(m.r * m.r + m.z * m.z);
}

} // namespace hd::ga


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// geodetic <-> ECEF <-> ENU
/////////////////////////////////////////////////////////////////////////////////////////

// geodetic position -> ECEF, returned as a unitized pga3dp point [m]:
//
//     X = (N + h) cos(lat) cos(lon)
//     Y = (N + h) cos(lat) sin(lon)
//     Z = (N (1 - e^2) + h) sin(lat)
//
// with N = el.N(lat) the radius of curvature in the prime vertical. The (1 - e^2) factor
// on Z is what makes lat the geodetic latitude rather than the geocentric one, and the
// height h enters all three components alike because it runs along the ellipsoid normal.
inline vec3dp geo_to_ecef(geo_pos const& p, ellipsoid const& el = wgs84)
{
    auto const m = detail::meridian_from_geodetic(el, p.lat, p.height);

    // the meridian section, turned to the position's longitude about the polar axis
    return vec3dp{m.r * std::cos(p.lon), m.r * std::sin(p.lon), m.z, 1.0};
}

// ECEF point -> geodetic position: undo the longitude rotation, then solve the meridian
// section (see detail::geodetic_from_meridian, which carries Bowring's method and its
// accuracy-versus-height table).
//
// On the polar axis the longitude is undefined and is taken as 0 by convention -- the
// same choice enu_at(vec3dp) makes for the frame there.
inline geo_pos ecef_to_geo(vec3dp const& P, ellipsoid const& el = wgs84)
{
    auto const p = unitize(P);

    value_t const r = std::sqrt(p.x * p.x + p.y * p.y); // distance from the polar axis
    value_t const lon = (r < eps * el.r_equator) ? 0.0 : std::atan2(p.y, p.x);

    auto const g = detail::geodetic_from_meridian(el, r, p.z);

    return geo_pos{g.lat, lon, g.height};
}


// the local ENU basis directions at a position, expressed in ECEF (unit directions,
// i.e. points at infinity with w = 0)
struct enu_frame {

    vec3dp east;
    vec3dp north;
    vec3dp up;
};

namespace detail {

// The frame construction shared by both enu_at() overloads, from the local vertical u
// and the horizontal meridian direction r (both unit, r towards the position):
//
//     N     = wdg(e3, r)   the meridian plane: the earth's axis joined with r
//     north = (u << N)     u contracted into that plane -- the in-plane direction
//                          perpendicular to up, pointing towards the pole
//     east  = dual(N)      the normal of the meridian plane
//
// Building N from r rather than from u is what keeps this non-degenerate: wdg(e3, u)
// equals cos(lat) * wdg(e3, r), so the two agree in direction everywhere (cos(lat) > 0)
// but only the r form survives at the poles, where u becomes parallel to e3. The frame
// is orthonormal and right-handed in the order (east, north, up).
//
// The ORIENTATION of N is the whole content of its line: taking wdg(r, e3) instead
// reverses N, and with it both east and north -- which leaves the frame orthonormal and
// still right-handed (it is a half turn about up), so only a pinned direction catches it.
inline enu_frame enu_from(vec3d const& u, vec3d const& r)
{
    auto const N = ega::wdg(ega::e3_3d, r);

    auto const n = ega::normalize(ega::operator<<(u, N)); // north: in the meridian plane
    auto const e = ega::normalize(ega::dual(N));          // east:  normal of that plane

    return enu_frame{vec3dp{e.x, e.y, e.z, 0.0}, vec3dp{n.x, n.y, n.z, 0.0},
                     vec3dp{u.x, u.y, u.z, 0.0}};
}

} // namespace detail

// WHAT UP IS HERE: the ELLIPSOID NORMAL at the position -- the local vertical, which is
// the convention used in practice (GNSS, surveying, and every "ENU"/"local tangent
// plane" frame in the literature). It is NOT the ray from the geocenter; the two differ
// by up to ~11.5' at mid-latitudes, which would tilt every locally computed direction by
// that much.
//
// The frame is built from the geometry of the position alone -- no angles are taken
// apart and no rotation matrix is written out (cf. A. Arsenovic, "What is Geometric
// Algebra", PyData Virginia 2025, which does this with a radial up):
//
//     g = grad of the ellipsoid through the position, i.e. the local vertical
//     N = wdg(e3, g)   the meridian plane: spanned by the earth's axis and the vertical
//     up    = g        the ellipsoid normal
//     north = (g << N) the direction in the meridian plane perpendicular to up, i.e. up
//                      contracted into that plane, pointing towards the pole
//     east  = dual(N)  the normal of the meridian plane, i.e. perpendicular to both the
//                      earth's axis and the vertical
//
// The three come out orthonormal after normalization, and right-handed in the order
// (east, north, up). Only directions are involved, so this is done in ega3d (the pga3dp
// directions carry the same three components with w = 0), and no transcendental function
// is evaluated at all.
//
// ACCURACY: taking the gradient AT the position rather than at its foot point on the
// ellipsoid is exact on the surface and degrades slowly with height -- below 0.11" at
// 1 km, below 1" at 9 km. Use the geo_pos overload (or enu_basis_at) if the position is
// far above the surface, since those read the geodetic latitude directly and are exact
// at any height.
inline enu_frame enu_at(vec3dp const& P, ellipsoid const& el = wgs84)
{
    // the frame is made of directions only, so it is built in ega3d (HINT 2). The ega
    // names must be qualified: unqualified lookup from inside hd::ga::pga finds the pga
    // overload set first and stops there, so a using-directive would not reach them.
    auto const p = unitize(P);

    value_t const aa = el.r_equator * el.r_equator;
    value_t const bb = el.r_pole * el.r_pole;

    // Both degeneracy tests below must be SCALE-FREE. The gradient carries a factor
    // 1/a^2, so its length is ~1e-7 for a point on the earth and nrm_sq(g) ~ 1e-14 --
    // already close to an absolute eps, and the meridian plane's ~1e-14 sin^2 falls
    // under it well before the pole. So the position is tested at its own scale, and the
    // normal is unitized before the plane is built.
    if (p.x * p.x + p.y * p.y + p.z * p.z < eps) {
        throw std::invalid_argument("enu_at: no local frame at the geocenter");
    }

    // the outward normal of x^2/a^2 + y^2/a^2 + z^2/b^2 = const through the position
    auto const u = ega::normalize(vec3d{p.x / aa, p.y / aa, p.z / bb}); // up

    // The meridian direction: the horizontal unit vector towards the position. Away from
    // the polar axis this is the position's own horizontal part; ON the axis a point
    // carries no longitude at all, and the frame is genuinely free to rotate about up --
    // every horizontal direction is south at the north pole. The prime meridian is then
    // taken by CONVENTION, which is the same choice ecef_to_geo makes when it reports
    // lon == 0 for a point on the axis. Use the geo_pos overload to name a meridian.
    value_t const r_xy = std::sqrt(p.x * p.x + p.y * p.y);
    auto const r = (r_xy > eps * std::sqrt(aa)) ? vec3d{p.x / r_xy, p.y / r_xy, 0.0}
                                                : vec3d{1.0, 0.0, 0.0};

    return detail::enu_from(u, r);
}

// The frame at a geodetic position. This overload has NO degenerate case: the meridian
// plane is built from the LONGITUDE, which is well defined at every latitude including
// the poles, rather than from the position, which stops carrying one there.
inline enu_frame enu_at(geo_pos const& p, ellipsoid const& el = wgs84)
{
    // on the surface the gradient IS the normal, so evaluate it there and stay exact
    auto const S = geo_to_ecef(geo_pos{p.lat, p.lon, 0.0}, el);

    value_t const aa = el.r_equator * el.r_equator;
    value_t const bb = el.r_pole * el.r_pole;

    auto const u = ega::normalize(vec3d{S.x / aa, S.y / aa, S.z / bb});
    auto const r = vec3d{std::cos(p.lon), std::sin(p.lon), 0.0};

    return detail::enu_from(u, r);
}

// The same frame written out from the angles instead -- the conventional form, kept
// alongside enu_at() so the two can be compared against each other. It evaluates one
// sine and one cosine per angle where enu_at() takes a wedge, a contraction, a dual and
// three normalizations.
//
// Up is the ellipsoid normal here too, and directly so: the geodetic latitude IS the
// angle of that normal, so no correction is applied to it. enu_at() is the primary
// definition and this one is the independent check on it.
inline enu_frame enu_basis_at(geo_pos const& p,
                              [[maybe_unused]] ellipsoid const& el = wgs84)
{
    value_t const sp = std::sin(p.lat);
    value_t const cp = std::cos(p.lat);
    value_t const sl = std::sin(p.lon);
    value_t const cl = std::cos(p.lon);

    return enu_frame{vec3dp{-sl, cl, 0.0, 0.0}, vec3dp{-sp * cl, -sp * sl, cp, 0.0},
                     vec3dp{cp * cl, cp * sl, sp, 0.0}};
}

// Motor taking ENU coordinates at the station to ECEF coordinates (body -> world):
//
//     M = T(P0) (x) Rz(lon + 90°) (x) Rx(90° - lat)
//
// read right to left -- Rx tilts the equatorial frame up to the station's latitude, Rz
// swings it onto the local meridian and puts e1 onto east, T moves it out to the station
// point. The two rotations ARE the longitude and the co-latitude, i.e. the same two
// angles that define the position; the latitude enters un-corrected because it already
// is the angle of the ellipsoid normal, which is what up is.
//
// The inverse map ECEF -> ENU is the regressive reverse rrev(M), so an ECEF point X is
// read in the local frame as move3dp(X, rrev(M)).
inline mvec3dp_e enu_to_ecef_motor(geo_pos const& p, ellipsoid const& el = wgs84)
{
    auto const P0 = geo_to_ecef(p, el);

    auto const M_tra = get_motor(vec3dp{P0.x, P0.y, P0.z, 0.0});
    auto const M_lon = get_motor(z_axis_3dp, p.lon + 0.5 * pi);
    auto const M_lat = get_motor(x_axis_3dp, 0.5 * pi - p.lat);

    // concatenation order is rgpr(applied_last, applied_first)
    return rgpr(M_tra, rgpr(M_lon, M_lat));
}

// Motor for a station given as an ECEF POINT rather than as a geo_pos -- the companion to
// enu_at(vec3dp), so a caller holding only a fix does not have to convert by hand.
//
// It deliberately does NOT build the motor from the frame enu_at() already computed,
// which looks like the natural unification. That needs a rotation carrying the reference
// frame onto the local one, and every "align a onto b" construction -- sqrt(rgpr(l2,
// rrev(l1))) via get_motor_from_lines, or the closed form 1 + sum(b_i a_i) -- is singular
// when that rotation is by pi, where its axis is undetermined. Measured on a 5-degree
// grid, the two-step line-alignment form fails across the WHOLE lon = 90 degree meridian
// (20 of 888 points, axes reversed), because there the second step must send e1 to -e1.
//
// enu_to_ecef_motor has no such inversion: longitude and co-latitude are the position's
// own coordinates, so it is correct at every latitude and longitude, poles included. That
// is why the frame and the motor are built by two different routes and cross-checked
// against each other in the tests, rather than one being derived from the other.
inline mvec3dp_e enu_motor_at(vec3dp const& P, ellipsoid const& el = wgs84)
{
    return enu_to_ecef_motor(ecef_to_geo(P, el), el);
}


/////////////////////////////////////////////////////////////////////////////////////////
// datum transformation: the 7-parameter Helmert transformation in ECEF
/////////////////////////////////////////////////////////////////////////////////////////
//
// Two geodetic datums -- an old national one on the Bessel ellipsoid, say, and ETRS89 --
// place their ECEF axes slightly differently: shifted by up to some hundred metres,
// turned by a few arc seconds, and with a scale differing by some parts per million. The
// Helmert transformation carries ECEF coordinates from one datum to the other,
//
//     X' = T + (1 + ds) R X
//
// with the translation T = (tx, ty, tz), a rotation R by the small angles (rx, ry, rz)
// and the scale correction ds. To convert geodetic positions between datums: geo_to_ecef
// on the source datum's ellipsoid, the Helmert transformation, ecef_to_geo on the target
// datum's ellipsoid.
//
// In PGA the rigid part IS a motor, M = T (x) R, and the scale about the geocentre is a
// change of the point's weight, so the transformation is one sandwich:
//
//     X' = move3dp(vec3dp(x, y, z, w / (1 + ds)), M)
//
// The rotation is the motor of the LINE through the geocentre with direction
// (rx, ry, rz), whose weight is the angle: rexp(0.5 * bivec3dp(rx, ry, rz, 0, 0, 0)) --
// a right-handed turn by |r| about that axis.
//
// TWO CONVENTIONS, published parameter sets use both, and they differ only in the sign of
// the rotation. Applying a parameter set with the wrong one is a silent error of about
// 2 |r| |X| -- 110 m for the German DHDN set (2.5 arc seconds, measured) -- so each has a
// function named after it:
//
// - helmert_position_vector()  : EPSG method 1033 (and 1037 geog3D), PROJ's towgs84; the
//                                rotation turns the POINT by +r
// - helmert_coordinate_frame() : EPSG method 1032 (and 1038), common in Europe; the
//                                rotation turns the AXES by +r, i.e. the point by -r
//
// The parameter set says which one it is meant for. EPSG writes R in its linearized form
// I + [r]x, valid to first order in the angles; the motor here is the EXACT rotation. The
// two differ by second-order terms, at most |r|^2 |X| / 2: 0.33 mm for the DHDN set, and
// 0.04 mm from PROJ's exact form, which composes three axis rotations instead of turning
// about one axis -- far below the accuracy of any published parameter set (both measured
// and tested). The inverse transformation is, to first order, the same function with
// every parameter negated (EPSG's "reversible" methods) -- only to first order: for the
// DHDN set the round trip misses by 9 mm (measured), the rotation acting on a
// translation of 740 m.

// the seven parameters: translation [m], rotation angles [rad], scale correction (ppm
// times 1e-6). Published sets give the angles in arc seconds and the scale in ppm --
// make_helmert7() takes those units
struct helmert7 {
    value_t tx{0.0}, ty{0.0}, tz{0.0}; // [m]
    value_t rx{0.0}, ry{0.0}, rz{0.0}; // [rad]
    value_t ds{0.0};                   // scale correction, dimensionless
};

// the parameters in the units they are published in: metres, arc seconds and ppm
inline helmert7 make_helmert7(value_t tx, value_t ty, value_t tz, value_t rx_arcsec,
                              value_t ry_arcsec, value_t rz_arcsec, value_t ds_ppm)
{
    value_t const as = deg2rad(1.0 / 3600.0);
    return helmert7{
        tx, ty, tz, rx_arcsec * as, ry_arcsec * as, rz_arcsec * as, ds_ppm * 1.0e-6};
}

namespace detail {

// X' = T + (1 + ds) R(sign * r) X, as one motor sandwich of the re-weighted point
inline vec3dp helmert_apply(vec3dp const& X, helmert7 const& h, value_t sign)
{
    auto const M_rot = rexp(0.5 * sign * bivec3dp{h.rx, h.ry, h.rz, 0.0, 0.0, 0.0});
    auto const M_tra = get_motor(vec3dp{h.tx, h.ty, h.tz, 0.0});
    // concatenation order is rgpr(applied_last, applied_first)
    auto const M = rgpr(M_tra, M_rot);
    return unitize(move3dp(vec3dp{X.x, X.y, X.z, X.w / (1.0 + h.ds)}, M));
}

} // namespace detail

// position vector convention (EPSG 1033, PROJ towgs84): ECEF point -> ECEF point
inline vec3dp helmert_position_vector(vec3dp const& X, helmert7 const& h)
{
    return detail::helmert_apply(X, h, 1.0);
}

// coordinate frame convention (EPSG 1032): ECEF point -> ECEF point
inline vec3dp helmert_coordinate_frame(vec3dp const& X, helmert7 const& h)
{
    return detail::helmert_apply(X, h, -1.0);
}


/////////////////////////////////////////////////////////////////////////////////////////
// projections onto the tangent plane (the azimuthal projections)
/////////////////////////////////////////////////////////////////////////////////////////
//
// A map of a region of the ellipsoid onto the plane that touches it at one point, the
// STATION. In the station's local frame (east, north, up -- see enu_to_ecef_motor) that
// plane is U = 0, i.e. xy_3dp, and every projection of this family is ONE construction:
// the line from a projection centre C through the point P, met with the plane,
//
//     X = rwdg(wdg(C, P), xy_3dp)
//
// The kinds below differ only in where C sits, and that choice decides what the map
// preserves. Cartography calls them the azimuthal projections: on a sphere each of them
// shows the true direction (azimuth) from the station to every point, and they differ in
// how distance from the station is scaled. Map coordinates are metres in the plane,
// x towards east and y towards north at the station.
//
// They are taken on the REAL ellipsoid: the plane touches the ellipsoid, and the centres
// are placed relative to it. On a sphere they reduce to the textbook forms (Snyder, "Map
// Projections -- A Working Manual", USGS 1987, sections 20-23), which the tests check.
// On the earth two of them ARE the EPSG ellipsoidal methods, which define them by this
// construction: the orthographic (EPSG 9840) and the vertical perspective (EPSG 9838)
// agree with it to 2e-9 m (tested). The EPSG oblique stereographic (9809) does not: it
// routes through a conformal sphere to stay exactly conformal on the ellipsoid, and the
// stereographic map here, conformal only on a sphere, differs from it by up to 3.2e-4 of
// the distance from the station -- 1.4 km across Europe seen from Berlin (measured). EPSG
// has no gnomonic method.
//
// Using one: make_azimuthal_map() once per station, then per point to_local(),
// is_drawable() and project(). The two steps are separate because a drawing pipeline
// needs the local point twice -- to decide whether it is drawable, and to project it --
// and because an outline that leaves the drawable region is cut there with clip_at(),
// which needs the local points too.

// which projection: where the centre C sits, and what that gives the map
enum class azimuthal {

    // C at infinity along up: every point drops perpendicularly onto the plane.
    //
    // The earth as it looks from very far away. Near the station the scale is close to
    // true; towards the limb features are compressed radially, to zero at the limb.
    // Neither
    // conformal nor equal-area. Only the near half is drawable, bounded by the limb.
    //
    // Used for pictorial views of the globe (the earth or a planet seen from space), for
    // showing which hemisphere faces a point, and wherever a map should look like a globe
    // rather than measure on it.
    orthographic,

    // C at the geocentre: every point is projected from the earth's centre.
    //
    // Every great circle lies in a plane through the centre, so it maps to a STRAIGHT
    // LINE -- the defining property. Scale grows fast away from the station (on a sphere
    // radially as 1/cos^2 of the angular distance: x4 at 60 deg), and less than a
    // hemisphere can be shown at all: the edge of the hemisphere maps to infinity.
    //
    // Used for great-circle routes in navigation (a straight line on the chart is the
    // shortest route, to be transferred to a conformal chart for steering), seismic and
    // radio wave paths, meteor trails, and the faces of cube-shaped global grids.
    //
    // On the ellipsoid the straight lines are the images of the plane sections through
    // the geocentre (great ellipses), and the geocentre lies OFF the station's normal
    // everywhere between the equator and the poles (21 km at 49 deg). That adds a
    // north-south tilt the three kinds with their centre on the normal do not have: over
    // a
    // region centred on the station it cancels, over one extending unevenly around it the
    // area excess against the orthographic map is no longer twice the stereographic one,
    // as it is for a centre on the normal at the same depth (tested).
    gnomonic,

    // C at the other end of the normal chord through the station (on a sphere: the
    // antipode).
    //
    // On a sphere CONFORMAL: angles and small shapes are true and circles map to circles
    // (on the ellipsoid nearly so; see the EPSG note above).
    // Scale grows slowly away from the station (on a sphere radially as 1/cos^2 of half
    // the angular distance: x2 at 90 deg), so a large region keeps recognisable shapes.
    //
    // Used for maps of a region around a point where shapes must be right: polar maps,
    // national grids built on the oblique stereographic (the Netherlands, Romania),
    // weather charts, and outside geodesy the celestial sphere and crystallography.
    stereographic,

    // C at a finite height above the station (azimuthal_options::height).
    //
    // What a camera at that height sees, drawn in the tangent plane: the vertical
    // perspective. Only the part of the earth visible from C is drawable, bounded by the
    // horizon as seen from there; the higher C, the closer the map comes to the
    // orthographic one, which is its limit at infinite height.
    //
    // Used for views from orbit or from an aircraft -- the earth as seen from a
    // geostationary satellite, the ISS or a balloon -- and for matching satellite images.
    perspective
};

// the per-kind settings; each applies to some kinds only and is ignored by the others
struct azimuthal_options {

    // perspective: height of the viewpoint above the station [m], must be positive. The
    // default is the geostationary altitude.
    value_t height{35786.0e3};

    // gnomonic, stereographic: the radius of the drawable region around the station, as
    // an angular distance [rad]; 0 selects the kind's default (60 deg / 90 deg). These
    // two have no horizon of their own -- the gnomonic map runs to infinity at 90 deg,
    // the stereographic one at 180 deg -- so where to stop is a choice, not a property.
    // The orthographic and perspective maps stop at their horizon.
    value_t max_angle{0.0};
};

// A projection onto the tangent plane at one station, ready to use. Built by
// make_azimuthal_map(); all its geometry is in the station's local frame.
struct azimuthal_map {

    geo_pos station; // the tangent point, on the ellipsoid (height 0)
    ellipsoid el;    // the ellipsoid the plane touches
    azimuthal kind;  // which projection
    mvec3dp_e M;     // local (east, north, up) -> ECEF; rrev(M) is ECEF -> local
    vec3dp C;        // the projection centre, local (an ideal point for orthographic)
    trivec3dp limit; // bounds the drawable region, local; wdg(P, limit) > 0 inside
};

namespace detail {

// The parameters t at which the line p + t d meets the ellipsoid surface
//
//     x^2/a^2 + y^2/a^2 + z^2/b^2 = 1
//
// for a point p and a direction d, both ECEF. A quadric is not a blade, so the algebra
// has no meet for it; this is the plain quadratic. Lengths are taken in units of a, so
// every term is of order one. Returns false when the line misses the ellipsoid, and
// otherwise the two roots with t1 <= t2 (equal for a tangent line).
inline bool line_meets_ellipsoid(vec3d const& p, vec3d const& d, ellipsoid const& el,
                                 value_t& t1, value_t& t2)
{
    value_t const a = el.r_equator;
    value_t const q = (a / el.r_pole) * (a / el.r_pole); // a^2/b^2 on the z terms

    vec3d const ps{p.x / a, p.y / a, p.z / a};
    vec3d const ds{d.x / a, d.y / a, d.z / a};

    // A t^2 + 2 B t + c = 0
    value_t const A = ds.x * ds.x + ds.y * ds.y + q * ds.z * ds.z;
    value_t const B = ps.x * ds.x + ps.y * ds.y + q * ps.z * ds.z;
    value_t const c = ps.x * ps.x + ps.y * ps.y + q * ps.z * ps.z - 1.0;

    value_t const disc = B * B - A * c;
    if (!(A > 0.0) || disc < 0.0) return false;

    // the cancellation-free pair of roots: the larger one from -(B + sign(B) sqrt), the
    // other from the product of the roots, c/A
    value_t const s = std::sqrt(disc);
    value_t const k = -(B + ((B >= 0.0) ? s : -s));

    value_t r1 = k / A;
    value_t r2 = (k != 0.0) ? c / k : r1;
    if (r1 > r2) std::swap(r1, r2);

    t1 = r1;
    t2 = r2;
    return true;
}

} // namespace detail

// The projection of the given kind onto the plane tangent to the ellipsoid at the
// station. The station's height is ignored: the plane touches the ellipsoid.
//
// The drawable region is bounded by a plane, `limit`:
//
// - orthographic and perspective: the POLAR PLANE of C with respect to the ellipsoid,
//   which meets the ellipsoid exactly in the horizon seen from C. For a point C = (c, w)
//   it is the set of x with x^T G c = w, G = diag(1/a^2, 1/a^2, 1/b^2), and a point P on
//   the ellipsoid is visible from C (C lies on the outer side of P's tangent plane)
//   exactly when P^T G c > w. One rule covers both kinds: for the orthographic map C is
//   ideal, w = 0, and the polar plane is the diametral plane conjugate to up.
//
// - gnomonic and stereographic: the plane U = -R (1 - cos(max_angle)) with R the
//   ellipsoid's radius at the station, i.e. the depth below the tangent plane of a point
//   at that angular distance on a sphere of that radius. It is a cut at a chosen radius,
//   approximate by design.
inline azimuthal_map make_azimuthal_map(geo_pos const& station, azimuthal kind,
                                        ellipsoid const& el = wgs84,
                                        azimuthal_options const& opt = {})
{
    geo_pos const s{station.lat, station.lon, 0.0};
    auto const M = enu_to_ecef_motor(s, el);
    auto const Mr = rrev(M);

    vec3dp C{};
    switch (kind) {
        case azimuthal::orthographic:
            C = e3_3dp; // the ideal point along up
            break;
        case azimuthal::gnomonic:
            C = move3dp(O_3dp, Mr); // the geocentre, read in the local frame
            break;
        case azimuthal::stereographic: {
            // the second point where the normal through the station meets the ellipsoid
            // (the first is the station itself, t = 0)
            auto const S = geo_to_ecef(s, el);
            auto const u = enu_at(s, el).up;
            value_t t1{}, t2{};
            detail::line_meets_ellipsoid(vec3d{S.x, S.y, S.z}, vec3d{u.x, u.y, u.z}, el,
                                         t1, t2);
            C = vec3dp{0.0, 0.0, (std::abs(t1) > std::abs(t2)) ? t1 : t2, 1.0};
            break;
        }
        case azimuthal::perspective:
            if (!(opt.height > 0.0)) {
                throw std::invalid_argument("make_azimuthal_map: the perspective "
                                            "viewpoint needs a positive height");
            }
            C = vec3dp{0.0, 0.0, opt.height, 1.0};
            break;
    }

    trivec3dp limit{};
    if (kind == azimuthal::orthographic || kind == azimuthal::perspective) {
        // the polar plane of C, built in ECEF and read in the local frame. Multiplied by
        // a^2, so that its normal carries the size of c and its offset the size of a^2 w
        auto const c = move3dp(C, M);
        value_t const a = el.r_equator;
        value_t const q = (a / el.r_pole) * (a / el.r_pole);
        limit = move3dp(trivec3dp{c.x, c.y, q * c.z, -a * a * c.w}, Mr);
    }
    else {
        bool const gn = (kind == azimuthal::gnomonic);
        value_t const alpha =
            (opt.max_angle > 0.0) ? opt.max_angle : deg2rad(gn ? 60.0 : 90.0);
        if (alpha >= (gn ? 0.5 * pi : pi)) {
            throw std::invalid_argument(
                gn ? "make_azimuthal_map: the gnomonic max_angle must stay below 90 deg"
                   : "make_azimuthal_map: the stereographic max_angle must stay below "
                     "180 deg");
        }
        value_t const depth = el.radius(s.lat) * (1.0 - std::cos(alpha));
        limit = trivec3dp{0.0, 0.0, 1.0, depth}; // the plane U = -depth
    }

    // orient the limit so that the station, the local origin, is on its positive side
    if (value_t(wdg(O_3dp, limit)) < 0.0) limit = -limit;

    return azimuthal_map{s, el, kind, M, C, limit};
}

// an ECEF point read in the map's local frame (east, north, up at the station)
inline vec3dp to_local(azimuthal_map const& m, vec3dp const& P)
{
    return move3dp(P, rrev(m.M));
}

// Is the local point inside the drawable region, i.e. on the station's side of the
// limit? For the orthographic and perspective maps and a point on the ellipsoid this is
// exactly "visible from C". The point must carry a positive weight, as every point from
// geo_to_ecef does.
inline bool is_drawable(azimuthal_map const& m, vec3dp const& P_local)
{
    return value_t(wdg(P_local, m.limit)) > 0.0;
}

// The map coordinates [m] of a local point: x towards east, y towards north at the
// station. The point must be drawable (is_drawable); outside the region the line from C
// may meet the plane on the far side, or not at all, and then unitize() throws.
inline vec2d project(azimuthal_map const& m, vec3dp const& P_local)
{
    auto const X = unitize(rwdg(wdg(m.C, P_local), xy_3dp));
    return vec2d{X.x, X.y};
}

// The same for a geodetic position, for single points (labels, markers). A pipeline
// over many points should call to_local() and is_drawable() itself, see above.
inline vec2d project(azimuthal_map const& m, geo_pos const& p)
{
    return project(m, to_local(m, geo_to_ecef(p, m.el)));
}

// The point where the segment A-B crosses the plane t -- the cut, when an outline leaves
// the drawable region (t = the map's limit). Only the line through A and B enters, so
// call it for a segment whose ends lie on opposite sides of t.
inline vec3dp clip_at(vec3dp const& A, vec3dp const& B, trivec3dp const& t)
{
    return unitize(rwdg(wdg(A, B), t));
}

// The inverse: the position on the ellipsoid that maps to the map coordinates x -- for a
// cursor readout, or to re-centre a map on a point picked in it.
//
// It is the line from C through the image point, met with the ellipsoid. Of the two
// crossings the one the projection uses is taken: the nearer to C when C is outside the
// ellipsoid (orthographic, perspective -- the visible side), the one on the image side of
// C when C is inside or on it (gnomonic, stereographic).
//
// Map coordinates outside the image of the ellipsoid (beyond the limb of an orthographic
// or perspective map) have no position. With hit_out given, *hit_out reports whether
// there is one and a default geo_pos is returned when not; without it, a miss throws.
inline geo_pos unproject(azimuthal_map const& m, vec2d const& x, bool* hit_out = nullptr)
{
    auto const X = vec3dp{x.x, x.y, 0.0, 1.0};
    // the direction from C towards X, as an ideal point: w_C (X - c) for a finite C,
    // -C for the ideal one
    auto const d = X * m.C.w - m.C;

    auto const Xe = move3dp(X, m.M);
    auto const de = move3dp(vec3dp{d.x, d.y, d.z, 0.0}, m.M);

    value_t t1{}, t2{};
    bool const hit = detail::line_meets_ellipsoid(vec3d{Xe.x, Xe.y, Xe.z},
                                                  vec3d{de.x, de.y, de.z}, m.el, t1, t2);
    if (hit_out) *hit_out = hit;
    if (!hit) {
        if (hit_out) return geo_pos{};
        throw std::domain_error(
            "unproject: these map coordinates are outside the image of the ellipsoid");
    }

    bool const c_outside =
        (m.kind == azimuthal::orthographic || m.kind == azimuthal::perspective);
    value_t const t = c_outside ? t1 : t2;

    return ecef_to_geo(vec3dp{Xe.x + t * de.x, Xe.y + t * de.y, Xe.z + t * de.z, 1.0},
                       m.el);
}

// The boundary of the drawable region: n points (ECEF, on the ellipsoid) around the
// curve where the map's limit plane cuts the ellipsoid -- the horizon of an orthographic
// or perspective map, the chosen cut of a gnomonic or stereographic one. Projected, it is
// the outline of the map: the edge of the globe in an orthographic view.
//
// A plane cuts an ellipsoid in an ellipse whose centre K lies on the diameter CONJUGATE
// to the plane: for the plane n.x = d and the ellipsoid x^T G x = 1,
//
//     K = d / (n^T G^-1 n) * G^-1 n
//
// and around K the curve is K + r(theta) v(theta) for the unit directions v in the
// plane, with r^2 = (1 - K^T G K) / (v^T G v): the cross term vanishes because K is the
// centre. (For the polar plane of C, G^-1 n is C itself, so K lies on the line from the
// geocentre towards C.) Lengths are taken in units of a, as in line_meets_ellipsoid.
inline std::vector<vec3dp> limit_outline(azimuthal_map const& m, size_t n = 360)
{
    value_t const a = m.el.r_equator;
    value_t const q = (a / m.el.r_pole) * (a / m.el.r_pole); // G = diag(1, 1, q) / a^2

    // the limit in ECEF, n.x + d0 = 0; in units of a: (a n).xi = -d0
    auto const L = move3dp(m.limit, m.M);
    vec3d const nn{a * L.x, a * L.y, a * L.z};
    value_t const d = -L.w;

    vec3d const gin{nn.x, nn.y, nn.z / q}; // G^-1 n, up to the common 1/a^2
    vec3d const K = (d / value_t(ega::dot(nn, gin))) * gin;
    value_t const r2_num = 1.0 - (K.x * K.x + K.y * K.y + q * K.z * K.z);
    if (!(r2_num > 0.0)) {
        throw std::domain_error("limit_outline: the limit plane misses the ellipsoid");
    }

    // an orthonormal basis of the plane's directions
    auto const nh = ega::normalize(nn);
    vec3d const seed =
        (std::abs(nh.x) < 0.9) ? vec3d{1.0, 0.0, 0.0} : vec3d{0.0, 1.0, 0.0};
    auto const e1 = ega::normalize(seed - value_t(ega::dot(seed, nh)) * nh);
    auto const e2 = ega::normalize(ega::dual(ega::wdg(nh, e1)));

    std::vector<vec3dp> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        value_t const th = 2.0 * pi * static_cast<value_t>(i) / static_cast<value_t>(n);
        vec3d const v = std::cos(th) * e1 + std::sin(th) * e2;
        value_t const r = std::sqrt(r2_num / (v.x * v.x + v.y * v.y + q * v.z * v.z));
        vec3d const xi = K + r * v;
        out.push_back(vec3dp{a * xi.x, a * xi.y, a * xi.z, 1.0});
    }
    return out;
}


/////////////////////////////////////////////////////////////////////////////////////////
// the two-dimensional case: ECEF section <-> geodetic, and the local frame
/////////////////////////////////////////////////////////////////////////////////////////
//
// The frame here is the meridian-plane section of ECEF, drawn the way one draws it:
//
//     e1 -> towards the equator (the prime meridian's direction)
//     e2 -> along the polar axis, towards the north pole
//
// which is exactly what geo_to_ecef() produces in three dimensions once the longitude
// rotation is stripped off: (r, z) becomes (e1, e2).

// geodetic position -> the meridian section, as a unitized pga2dp point [m]
inline vec2dp geo_to_ecef(geo_pos2dp const& p, ellipsoid const& el = wgs84)
{
    auto const m = hd::ga::detail::meridian_from_geodetic(el, p.lat, p.height);

    return vec2dp{m.r, m.z, 1.0};
}

// the meridian section -> geodetic position (see detail::geodetic_from_meridian for the
// method and its accuracy against height).
//
// (The 3D projections onto the tangent plane are further up, after enu_motor_at.)
//
// A point with e1 < 0 lies on the far half of the ellipse -- the 180-degree meridian --
// which has no geodetic latitude in [-90, +90] and so no representation here; in three
// dimensions that half is reached by the longitude instead. It is refused rather than
// folded onto a wrong latitude.
inline geo_pos2dp ecef_to_geo(vec2dp const& P, ellipsoid const& el = wgs84)
{
    auto const p = unitize(P);

    if (p.x < -eps * el.r_equator) {
        throw std::invalid_argument(
            "ecef_to_geo: a point with e1 < 0 is on the 180-degree meridian, which a "
            "single meridian section does not represent (use the 3D overload)");
    }

    auto const g =
        hd::ga::detail::geodetic_from_meridian(el, std::max(p.x, value_t(0.0)), p.y);

    return geo_pos2dp{g.lat, g.height};
}


// The local frame in the meridian plane, expressed in the section's own coordinates
// (unit directions, i.e. points at infinity with w = 0).
//
// The order is (UP, NORTH), vertical first, and that is forced rather than chosen: with
// e1 towards the equator and e2 towards the pole, north = (-sin lat, cos lat) and
// up = (cos lat, sin lat) give
//
//     north ^ up = -e12      but      up ^ north = +e12
//
// so the pair (north, up) is NEGATIVELY oriented. A motor is a rotation, so no motor
// could ever carry (e1, e2) onto (north, up) -- that map is a reflection. Taken in the
// order (up, north) the frame is a proper rotation, and a particularly simple one: it is
// the rotation by the latitude itself, since at latitude zero the vertical IS e1.
struct un_frame {

    vec2dp up;
    vec2dp north;
};

// The frame from the geometry of the position alone. In two dimensions this needs no
// meridian plane -- there is no longitude to lose -- so it is simply the ellipse normal
// and its quarter turn, and unlike the 3D case it has NO degenerate configuration except
// the geocenter. The pole is unremarkable here: the 3D pole degeneracy comes from the
// longitude, which does not exist in a meridian section.
//
//     up    = grad of the ellipse through the position, normalized
//     north = up * I_2d      a quarter turn towards the pole (the ega2d pseudoscalar
//                            rotates a vector by +90 degrees)
inline un_frame un_at(vec2dp const& P, ellipsoid const& el = wgs84)
{
    auto const p = unitize(P);

    if (p.x * p.x + p.y * p.y < eps) {
        throw std::invalid_argument("un_at: no local frame at the geocenter");
    }

    value_t const aa = el.r_equator * el.r_equator;
    value_t const bb = el.r_pole * el.r_pole;

    // ega names must be qualified: unqualified lookup from inside hd::ga::pga finds the
    // pga overload set first and stops there
    auto const u = ega::normalize(vec2d{p.x / aa, p.y / bb});
    auto const n = ega::operator*(u, ega::I_2d);

    return un_frame{vec2dp{u.x, u.y, 0.0}, vec2dp{n.x, n.y, 0.0}};
}

inline un_frame un_at(geo_pos2dp const& p, ellipsoid const& el = wgs84)
{
    // on the surface the gradient IS the normal, so evaluate it there and stay exact
    return un_at(geo_to_ecef(geo_pos2dp{p.lat, 0.0}, el), el);
}

// The same frame written out from the angle instead -- the conventional form, kept
// alongside un_at() so the two cross-check, exactly as enu_basis_at() does in 3D.
inline un_frame un_basis_at(geo_pos2dp const& p,
                            [[maybe_unused]] ellipsoid const& el = wgs84)
{
    value_t const s = std::sin(p.lat);
    value_t const c = std::cos(p.lat);

    return un_frame{vec2dp{c, s, 0.0}, vec2dp{-s, c, 0.0}};
}

// Motor taking local (up, north) coordinates to the meridian section (body -> world):
//
//     M = T(P0) (x) R(lat)
//
// ONE rotation, where the 3D motor needs two -- the longitude rotation has nothing to
// turn. And the remaining rotation is the latitude itself, unmodified: it is already the
// angle of the ellipse normal, which is what up is.
//
// The inverse map, section -> local, is the regressive reverse rrev(M).
inline mvec2dp_u un_to_ecef_motor(geo_pos2dp const& p, ellipsoid const& el = wgs84)
{
    auto const P0 = geo_to_ecef(p, el);

    auto const M_tra = get_motor(vec2dp{P0.x, P0.y, 0.0});
    auto const M_rot = get_motor(O_2dp, p.lat); // about the geocenter

    return rgpr(M_tra, M_rot); // rotate first, then translate
}

// the same motor for a station given as a point -- companion to un_at(vec2dp)
inline mvec2dp_u un_motor_at(vec2dp const& P, ellipsoid const& el = wgs84)
{
    return un_to_ecef_motor(ecef_to_geo(P, el), el);
}

} // namespace hd::ga::pga

// Printing support for the types above -- and it has to stay HERE, at the very end.
// ga_fmt_geodesics.hpp includes this header back (so it compiles on its own), but
// that include is a no-op by the time it runs: the #pragma once guard is already
// set. Hoisted to the top, the formatter specialisations are therefore written
// before the types exist -- measured: "no member named 'ellipsoid' in namespace
// 'hd::ga'", 7 errors, both standalone and through ga_pga.hpp.
#include "detail/fmt/ga_fmt_geodesics.hpp"
