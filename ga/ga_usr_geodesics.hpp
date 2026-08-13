#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <cctype>    // std::isspace, std::toupper
#include <cmath>     // std::sin, std::cos, std::sqrt, std::atan2, std::abs
#include <cstdlib>   // std::strtod
#include <stdexcept> // std::invalid_argument
#include <string>    // std::string

#include "ga_ega3d_ops.hpp"     // wdg, dual, normalize, ... for the direction frame
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
// - to_geo_pos()              : geo_pos_dms -> geo_pos (the explicit conversion above)
// - distance_from_geocenter() : the position's "total radius" [m], derived
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
// N is the separation of the two surfaces. It is a measured field, not derivable from
// lat/lon by any formula, and reaches some tens of metres -- about +45 m around Berlin,
// about +50 m over Madrid, and roughly -30 m over the Indian Ocean. See to_geo_pos(),
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
    value_t const sp = std::sin(p.lat);
    value_t const cp = std::cos(p.lat);
    value_t const sl = std::sin(p.lon);
    value_t const cl = std::cos(p.lon);

    value_t const N = el.N(p.lat);

    return vec3dp{(N + p.height) * cp * cl, (N + p.height) * cp * sl,
                  (N * (1.0 - el.e_sq()) + p.height) * sp, 1.0};
}

// ECEF point -> geodetic position, via Bowring's closed-form solution (1976): the
// auxiliary (parametric) latitude theta = atan2(Z a, r b) turns the implicit relation
// between geodetic latitude and height into an explicit one.
//
// It is a SINGLE pass, not an iteration, so its accuracy depends on the height -- the
// approximation is built around a point near the surface. Measured round-trip error of
// ecef_to_geo(geo_to_ecef(p)) at latitude 45 deg:
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
inline geo_pos ecef_to_geo(vec3dp const& P, ellipsoid const& el = wgs84)
{
    auto const p = unitize(P);

    value_t const a = el.r_equator;
    value_t const b = el.r_pole;
    value_t const e2 = el.e_sq();
    value_t const ep2 = e2 / ((1.0 - el.flattening()) * (1.0 - el.flattening()));

    value_t const r = std::sqrt(p.x * p.x + p.y * p.y); // distance from the polar axis

    if (r < eps * a) { // on the polar axis: the longitude is undefined, take 0
        value_t const lat = (p.z >= 0.0) ? 0.5 * pi : -0.5 * pi;
        return geo_pos{lat, 0.0, std::abs(p.z) - b};
    }

    value_t const theta = std::atan2(p.z * a, r * b);
    value_t const st = std::sin(theta);
    value_t const ct = std::cos(theta);

    value_t const lat =
        std::atan2(p.z + ep2 * b * st * st * st, r - e2 * a * ct * ct * ct);
    value_t const lon = std::atan2(p.y, p.x);

    value_t const N = el.N(lat);
    value_t const sp = std::sin(lat);

    // near the poles r -> 0 makes h = r/cos(lat) - N ill-conditioned; switch to the
    // Z-based form there (each is well-conditioned where the other is not)
    value_t const h =
        (std::abs(sp) > 0.5) ? p.z / sp - N * (1.0 - e2) : r / std::cos(lat) - N;

    return geo_pos{lat, lon, h};
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

} // namespace hd::ga::pga

#include "detail/fmt/ga_fmt_geodesics.hpp" // printing support for the types above
