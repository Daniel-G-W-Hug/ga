#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <cctype>    // std::isspace, std::toupper
#include <cmath>     // std::sin, std::cos, std::sqrt, std::atan2, std::abs
#include <cstdio>    // std::snprintf (deg2dms)
#include <cstdlib>   // std::strtod
#include <stdexcept> // std::invalid_argument
#include <string>    // std::string

#include "ga_ega2d_ops.hpp"     // normalize, operator*(vec, I_2d), ... for the 2D frame
#include "ga_ega3d_ops.hpp"     // wdg, dual, normalize, ... for the direction frame
#include "ga_pga2dp_ops.hpp"    // get_motor, move2dp, unitize, ... for the 2D case
#include "ga_pga3dp_ops.hpp"    // get_motor, move3dp, unitize, ... for points/motors
#include "ga_usr_consts.hpp"    // e3_3d, I_3d, O_3dp, x_axis_3dp, z_axis_3dp, ...
#include "ga_usr_types.hpp"     // vec3d, vec3dp, mvec3dp_e
#include "ga_usr_utilities.hpp" // deg2rad, rad2deg, pi
#include "ga_value_t.hpp"       // value_t

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

#include "detail/fmt/ga_fmt_geodesics.hpp" // printing support for the types above
