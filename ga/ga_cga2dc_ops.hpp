#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_cga2dc_ops_basics.hpp"
#include "ga_cga2dc_ops_products.hpp"

#include <cmath>     // dipole_points
#include <complex>   // exp/log/sqrt of regressive versors (central subalgebra = C)
#include <stdexcept> // dipole_points: the imaginary dipole
#include <utility>   // std::pair


namespace hd::ga::cga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides cga2dc geometric operations (layer under construction):
//
// - rexp()                  -> versor from a bivector generator (w.r.t. rgpr,
//                              closed form)
// - rlog()                  -> principal generator of an even unit versor
//                              (w.r.t. rgpr, inverse of exp)
// - rsqrt()                 -> square root of an even unit versor (w.r.t. rgpr)
// - get_translation()       -> motor translating by (dx, dy)
// - get_rotation()          -> motor rotating by theta about a point
// - get_dilation()          -> motor scaling by sigma about a point
// - transform()             -> apply motor: sandwich M (v) u (v) rrev(M)
// - invert_on()             -> inversion in a circle or line (flector sandwich)
//
// Object construction (grade encodes the object kind: round point = vector,
// dipole / flat point = bivector, circle / line = trivector):
//
// - round_point2dc()        -> round point from center and radius (r = 0: the
//                              null-point embedding of a Euclidean point)
// - flat_point2dc()         -> flat point from a Euclidean position
// - dipole2dc()             -> dipole from center, radius and normal
// - circle2dc()             -> circle from center and radius
//
// (the bulk/weight split, the part/center norms, unitize and cconj live in
// ga_cga2dc_ops_basics.hpp with the other algebras' norm layer; the radius
// norms are HERE because they need the product-layer rdot)
//
// - radius_nrm_sq()         -> the antidot square u (o) u = r^2 (signed:
//                              negative for imaginary radii)
// - radius_sq()             -> squared radius (radius_nrm_sq / round_weight^2)
//
// Geometric properties of the round objects:
//
// - car()                   -> carrier: the lowest-dimensional flat containing u
// - ccr()                   -> cocarrier: the carrier of the antidual
// - cen()                   -> center as a round point (same position + radius)
// - con()                   -> container: the smallest circle containing u
// - par()                   -> partner: same center/carrier, r^2 negated
// - att()                   -> attitude (direction content; removes the origin)
//
// Worked constructions:
//
// - dipole_points()         -> the two points of a dipole (the end of a meet chain)
// - two_link_ik2dc()        -> planar two-link inverse kinematics as a circle meet
//
// - is_congruent()          -> same subspace up to a non-zero scalar factor
// - is_close()              -> same value within a RELATIVE tolerance
// - is_same_transform()     -> do two motors act as the same conformal map?
//                              (versors double-cover the transformations)
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// exponential, logarithm and square root w.r.t. the regressive geometric
// product -- the versor side of the algebra: conformal transformations are
// applied as rgpr sandwiches, so their generators and roots live here.
//
// Under rgpr the roles of scalar and pseudoscalar are exchanged: the
// PSEUDOSCALAR I is the identity element and the SCALAR acts as an imaginary
// unit (rgpr(1,1) = -I). The span {I, 1} is central in the even subalgebra and
// isomorphic to the complex numbers via
//
//     a*I + b*1  <->  a + i b
//
// For a bivector B the regressive square rgpr(B,B) = a*I + b*1 lies in this
// center (the grade-2 part of any square vanishes), so the exponential series
// sums in closed form with complex arithmetic:
//
//     rexp(B) = cosh(w) + rgpr(B, sinh(w)/w),    w = rsqrt(a + i b)
//
// with cosh(w) and sinh(w)/w read back through the isomorphism (both are even
// in w, so the branch of the square root does not matter). This is the
// non-degenerate analog of the dual-number closed form of the degenerate pga
// motors. The 6-dimensional bivector space generates the full 6-dimensional
// conformal group of the plane: rotations, translations, dilations and
// transversions (translations are the parabolic case rgpr(B,B) = 0, where the
// series terminates as rexp(B) = I + B).
//
// rlog() returns the principal generator: complex acosh of the central part,
// then B = rgpr(V, S^-1) with V the bivector part and S = sinh(w)/w. It
// expects an even UNIT versor (exp always produces one) and throws where the
// logarithm is singular (S = 0, i.e. a full-turn versor with central part
// -I). rsqrt(M) = rexp(0.5 * rlog(M)) inherits the principal branch and the
// same guard.
////////////////////////////////////////////////////////////////////////////////
// HINT: exp/log/sqrt cannot be constexpr (std::complex math; log may throw)

// helper: map a + i b back to the central element a*I + b*1 and multiply a
// bivector by a central element: rgpr(B, a*I + b*1) = a*B + b*rgpr(B, 1)
template <typename T>
    requires(numeric_type<T>)
inline BiVec2dc<T> rgpr_central(BiVec2dc<T> const& B, std::complex<T> const& s)
{
    return s.real() * B + s.imag() * rgpr(B, Scalar2dc<T>(T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_E<T> rexp(BiVec2dc<T> const& B)
{
    auto const sq = rgpr(B, B); // a*I + b*1, central
    std::complex<T> const z(T(gr4(sq)), T(gr0(sq)));
    std::complex<T> const w = std::sqrt(z);
    std::complex<T> const ch = std::cosh(w);
    // sinh(w)/w, series for small |z| (even in w, so express it in z):
    //     sinh(w)/w = 1 + z/6 + z^2/120 + ...
    std::complex<T> const sh_ow = (std::abs(z) > T(1e-4))
                                      ? std::sinh(w) / w
                                      : T(1.0) + z / T(6.0) + z * z / T(120.0);
    return MVec2dc_E<T>(Scalar2dc<T>(ch.imag()), rgpr_central(B, sh_ow),
                        PScalar2dc<T>(ch.real()));
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec2dc<T> rlog(MVec2dc_E<T> const& M)
{
    std::complex<T> const c(T(gr4(M)), T(gr0(M))); // central part
    BiVec2dc<T> const V = gr2(M);                  // V = rgpr(B, sinh(w)/w)
    std::complex<T> const w = std::acosh(c);       // principal branch
    std::complex<T> const w2 = w * w;
    std::complex<T> const sh_ow = (std::abs(w) > T(1e-2))
                                      ? std::sinh(w) / w
                                      : T(1.0) + w2 / T(6.0) + w2 * w2 / T(120.0);
    // singular exactly at the full-turn versors (central part -I): any
    // generator with rgpr(B,B) = -pi^2 * I exponentiates to -I
    hd::ga::detail::check_normalization<T>(std::abs(sh_ow),
                                           "regressive logarithm (full-turn versor)");
    return rgpr_central(V, T(1.0) / sh_ow);
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_E<T> rsqrt(MVec2dc_E<T> const& M)
{
    return rexp(T(0.5) * rlog(M));
}


////////////////////////////////////////////////////////////////////////////////
// conformal transformations
//
// An even unit versor (motor) M transforms any object with the regressive
// sandwich
//
//     u' = transform(u, M) = M (v) u (v) rrev(M)
//
// Reflections generalize to INVERSIONS in circles (and lines, which are
// circles through infinity): an odd sandwich with the trivector c,
//
//     u' = invert_on(u, c) = -( c (v) u (v) rrev(c) )
//
// (the leading minus is the flector-sandwich convention of the reference).
// Composing two inversions gives a motor: parallel lines -> translation,
// intersecting lines -> rotation, concentric circles -> dilation.
//
// The named builders produce unit motors directly from the exponential. All of
// them are generated by FLAT POINTS: a finite flat point at m has the
// regressive square -I (elliptic), so
//
//     rotation about m by theta:  M = rexp( (theta/2) * flat_point(m) )
//
// its dual is the "imaginary circle in the horizon centered at m" with
// regressive square +I (hyperbolic), so with delta = -1/2 ln(sigma)
//
//     dilation about m by sigma:  M = rexp( delta * r_dual(flat_point(m)) )
//
// and a flat point IN the horizon (zero weight) is parabolic (regressive
// square 0), so the translation exponential terminates:
//
//     translation by (dx, dy):    M = I + B,  B = (dy/2) e14 - (dx/2) e24
//
// Weights under sandwiches are homogeneous: a non-unit versor scales all
// image weights uniformly (odd sandwiches may flip their sign); normalize
// results with unitize() where a weight-one representative is needed.
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> transform(Vec2dc<T> const& v,
                                                     MVec2dc_E<U> const& M)
{
    return gr1(rgpr(rgpr(M, v), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> transform(BiVec2dc<T> const& B,
                                                       MVec2dc_E<U> const& M)
{
    return gr2(rgpr(rgpr(M, B), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> transform(TriVec2dc<T> const& t,
                                                        MVec2dc_E<U> const& M)
{
    return gr3(rgpr(rgpr(M, t), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec2dc<std::common_type_t<T, U>> invert_on(Vec2dc<T> const& v,
                                                     TriVec2dc<U> const& c)
{
    return -gr1(rgpr(rgpr(c, v), rrev(c)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec2dc<std::common_type_t<T, U>> invert_on(BiVec2dc<T> const& B,
                                                       TriVec2dc<U> const& c)
{
    return -gr2(rgpr(rgpr(c, B), rrev(c)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec2dc<std::common_type_t<T, U>> invert_on(TriVec2dc<T> const& t,
                                                        TriVec2dc<U> const& c)
{
    return -gr3(rgpr(rgpr(c, t), rrev(c)));
}

// translation by (dx, dy): parabolic, the exponential series terminates
template <typename T>
    requires(numeric_type<T>)
constexpr MVec2dc_E<T> get_translation(T dx, T dy)
{
    return MVec2dc_E<T>(
        Scalar2dc<T>(T(0.0)),
        BiVec2dc<T>(T(0.0), T(0.0), T(0.0), T(0.5) * dy, -T(0.5) * dx, T(0.0)),
        PScalar2dc<T>(T(1.0)));
}

// rotation about (mx, my) by theta (counterclockwise for theta > 0); the
// generator is the unit flat point at (mx, my) -- written out inline, see
// flat_point2dc() below
template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_E<T> get_rotation(T mx, T my, T theta)
{
    auto const fp = BiVec2dc<T>(T(0.0), T(0.0), T(0.0), -mx, -my, T(-1.0));
    return rexp(BiVec2dc<T>(T(0.5) * theta * fp));
}

// dilation about (mx, my) by the scale factor sigma > 0; the generator is the
// dual of the unit flat point at (mx, my) -- written out inline, see
// flat_point2dc() below
template <typename T>
    requires(numeric_type<T>)
inline MVec2dc_E<T> get_dilation(T mx, T my, T sigma)
{
    hd::ga::detail::check_normalization<T>(sigma, "dilation scale factor");
    T const delta = -T(0.5) * std::log(sigma);
    auto const fp = BiVec2dc<T>(T(0.0), T(0.0), T(0.0), -mx, -my, T(-1.0));
    return rexp(BiVec2dc<T>(delta * r_dual(fp)));
}


////////////////////////////////////////////////////////////////////////////////
// geometric object construction
//
// The round objects and their component semantics (fields as stored; the
// object kind is encoded in the grade):
//
//     round point (vector):   a = (ax, ay, az, aw), center (ax,ay)/az,
//                             aw = (p^2 + r^2)/2 * az; null for r = 0
//     dipole (bivector):      carrier-line part (vx, vy, vz), cocarrier
//                             position part (mx, my, mz); a dipole with zero
//                             v-part is a FLAT POINT
//     circle (trivector):     flat-line part (x, y, z), carrier weight w; a
//                             circle with w = 0 is a LINE (through infinity)
//
// The circle's weight component is intentionally NEGATIVE for the standard
// construction so the center can be read directly off the flat-line part.
////////////////////////////////////////////////////////////////////////////////

// round point with center (px, py) and radius r (unitized: origin weight 1);
// r = 0 gives the null-point embedding q with q * q = 0
template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> round_point2dc(T px, T py, T r)
{
    return Vec2dc<T>(px, py, T(1.0), T(0.5) * (px * px + py * py + r * r));
}

// flat point at Euclidean position (px, py) (unitized)
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> flat_point2dc(T px, T py)
{
    return BiVec2dc<T>(T(0.0), T(0.0), T(0.0), -px, -py, T(-1.0));
}

// dipole with center (px, py), radius r and normal (nx, ny); its two points
// lie at distance r from the center, PERPENDICULAR to the normal (in 2D the
// normal is perpendicular to the carrier line, as for lines)
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> dipole2dc(T px, T py, T r, T nx, T ny)
{
    T const pn = px * nx + py * ny; // p . n
    T const hs = T(0.5) * (px * px + py * py + r * r);
    T const pxn = px * ny - py * nx; // p x n (2d cross product)
    T const gx = nx;
    T const gy = ny;
    T const gz = -pn;
    T const ppx = hs * ny - pxn * px;
    T const ppy = -hs * nx - pxn * py;
    T const ppz = -pxn;
    return BiVec2dc<T>(gy, -gx, gz, -ppx, -ppy, -ppz);
}

// circle with center (px, py) and radius r (unitized, carrier weight -1)
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> circle2dc(T px, T py, T r)
{
    return TriVec2dc<T>(py, -px, T(-0.5) * (px * px + py * py - r * r), T(-1.0));
}


////////////////////////////////////////////////////////////////////////////////
// radius norm and squared radius
//
// The radius norm square is the ANTIDOT square u (o) u = rdot(u, u) = r^2 (for
// a unitized object; the dot square gives dot(u, u) = -r^2). It is SIGNED:
// negative values mean an imaginary radius (e.g. the meet of two
// non-intersecting circles), so no square root is taken here. These live in
// the geometric layer (not in ops_basics with the other norms) because rdot is
// a product-layer operation.
////////////////////////////////////////////////////////////////////////////////

// the antidot square (SIGNED; = r^2 for a unitized object)
template <typename T>
    requires(numeric_type<T>)
constexpr T radius_nrm_sq(Vec2dc<T> const& v)
{
    return T(rdot(v, v));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T radius_nrm_sq(BiVec2dc<T> const& B)
{
    return T(rdot(B, B));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T radius_nrm_sq(TriVec2dc<T> const& t)
{
    return T(rdot(t, t));
}

// squared radius of a round object of any weight (throws for flat objects,
// whose round weight is zero)
// HINT: radius_sq()/unitize() cannot be constexpr due to the checks for
//       division by zero which might throw

template <typename T>
    requires(numeric_type<T>)
inline T radius_sq(Vec2dc<T> const& v)
{
    T const wn_sq = round_weight_nrm_sq(v);
    hd::ga::detail::check_normalization<T>(wn_sq, "round point (round weight)");
    return radius_nrm_sq(v) / wn_sq;
}

template <typename T>
    requires(numeric_type<T>)
inline T radius_sq(BiVec2dc<T> const& B)
{
    T const wn_sq = round_weight_nrm_sq(B);
    hd::ga::detail::check_normalization<T>(wn_sq, "dipole (round weight)");
    return radius_nrm_sq(B) / wn_sq;
}

template <typename T>
    requires(numeric_type<T>)
inline T radius_sq(TriVec2dc<T> const& t)
{
    T const wn_sq = round_weight_nrm_sq(t);
    hd::ga::detail::check_normalization<T>(wn_sq, "circle (round weight)");
    return radius_nrm_sq(t) / wn_sq;
}


////////////////////////////////////////////////////////////////////////////////
// geometric properties of the round objects
//
// All are defined coordinate-free over the products already in place (the
// component tables they produce were verified against the reference
// literature entry by entry):
//
//     car(u) = u ^ e4                    (grade k -> k+1; annihilates the flat
//                                         part, lifts the round part)
//     ccr(u) = l_antidual(u) ^ e4        (= car of the antidual)
//     cen(u) = rwdg( ccr(u), u )         (center as a ROUND POINT: position
//                                         AND radius of u)
//     con(u) = u ^ l_antidual( car(u) )  (smallest circle containing u)
//     par(u) = (-1)^(gr(u)+1) rwdg( con(l_antidual(u)), car(u) )
//                                        (same center + carrier, r^2 negated;
//                                         u . par(u) = 0 always)
//     att(u) = rwdg( u, r_cmpl(e3) )     (removes the origin content; for a
//                                         dipole: carrier direction + weight)
//
// The wedge with e4 is taken on the RIGHT, matching this library's trailing
// round-dimension blade ordering (the same convention the reference uses for
// its own trailing ordering).
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> car(Vec2dc<T> const& v)
{
    return wdg(v, Vec2dc<T>(T(0.0), T(0.0), T(0.0), T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> car(BiVec2dc<T> const& B)
{
    return wdg(B, Vec2dc<T>(T(0.0), T(0.0), T(0.0), T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> car(TriVec2dc<T> const& t)
{
    return wdg(t, Vec2dc<T>(T(0.0), T(0.0), T(0.0), T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> ccr(Vec2dc<T> const& v)
{
    return car(l_antidual(v));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> ccr(BiVec2dc<T> const& B)
{
    return car(l_antidual(B));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> ccr(TriVec2dc<T> const& t)
{
    return car(l_antidual(t));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> cen(Vec2dc<T> const& v)
{
    return rwdg(ccr(v), v);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> cen(BiVec2dc<T> const& B)
{
    return rwdg(ccr(B), B);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> cen(TriVec2dc<T> const& t)
{
    return rwdg(ccr(t), t);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> con(Vec2dc<T> const& v)
{
    return wdg(v, l_antidual(car(v)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> con(BiVec2dc<T> const& B)
{
    return wdg(B, l_antidual(car(B)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> con(TriVec2dc<T> const& t)
{
    return wdg(t, l_antidual(car(t)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> par(Vec2dc<T> const& v)
{
    // gr = 1: sign +1
    return rwdg(con(l_antidual(v)), car(v));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> par(BiVec2dc<T> const& B)
{
    // gr = 2: sign -1
    return T(-1.0) * rwdg(con(l_antidual(B)), car(B));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> par(TriVec2dc<T> const& t)
{
    // gr = 3: sign +1
    return rwdg(con(l_antidual(t)), car(t));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> att(Vec2dc<T> const& v)
{
    return rwdg(v, r_cmpl(Vec2dc<T>(T(0.0), T(0.0), T(1.0), T(0.0))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> att(BiVec2dc<T> const& B)
{
    return rwdg(B, r_cmpl(Vec2dc<T>(T(0.0), T(0.0), T(1.0), T(0.0))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> att(TriVec2dc<T> const& t)
{
    return rwdg(t, r_cmpl(Vec2dc<T>(T(0.0), T(0.0), T(1.0), T(0.0))));
}


////////////////////////////////////////////////////////////////////////////////
// worked constructions: the points of a dipole, two-link inverse kinematics
////////////////////////////////////////////////////////////////////////////////

// the two points of a dipole d, with d unitized:
//
//     p_pm = cen(d) +/- sqrt(radius_sq(d)) * att(d)
//
// returned as null points of unit weight (round_point2dc(x, y, 0)). This is the end
// of the meet chain (circle v circle, circle v line): the sign of radius_sq(d)
// decides beforehand whether real points exist, so an imaginary dipole (radius_sq
// < 0) throws -- test the sign first where the meet may miss. A flat object has no
// round weight and throws in radius_sq.
template <typename T>
    requires(numeric_type<T>)
inline std::pair<Vec2dc<T>, Vec2dc<T>> dipole_points(BiVec2dc<T> const& d)
{
    T const r_sq = radius_sq(d);
    if (r_sq < T(0.0))
        throw std::runtime_error("dipole_points: imaginary dipole, no real points");
    auto const du = unitize(d);
    auto const c = cen(du);
    auto const a = att(du);
    T const r = std::sqrt(r_sq);
    T const cx = c.x / c.z, cy = c.y / c.z; // center at weight one
    return {round_point2dc(cx + r * a.x, cy + r * a.y, T(0.0)),
            round_point2dc(cx - r * a.x, cy - r * a.y, T(0.0))};
}

// planar two-link inverse kinematics as a circle meet (the 2R arm, Lynch & Park
// 6.1; a leg from hip to ankle): the joint between a link of length L1 from the
// base point b and a link of length L2 to the tip point t lies on both circles, so
//
//     j_pm = dipole_points( circle2dc(b, L1) v circle2dc(t, L2) )
//
// are the two configurations (elbow up / elbow down, knee forward / back). No
// discriminant is formed and no branch is taken: the meet is one product, and the
// tip is out of reach (|t - b| > L1 + L2 or < |L1 - L2|) when that dipole is
// imaginary, which throws; test radius_sq of the meet first where reach is in
// question. Coincident base and tip (concentric circles) throw as a flat meet.
template <typename T>
    requires(numeric_type<T>)
inline std::pair<Vec2dc<T>, Vec2dc<T>> two_link_ik2dc(T bx, T by, T tx, T ty, T L1, T L2)
{
    return dipole_points(rwdg(circle2dc(bx, by, L1), circle2dc(tx, ty, L2)));
}


////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar2dc<T> a, Scalar2dc<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero scalars are congruent (represent the same 0-dimensional subspace)
    return true;
}

// For vectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Vec2dc<T> const& a, Vec2dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, tolerance);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec2dc<T> const& a, BiVec2dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<6>({value_t(a.vx), value_t(a.vy), value_t(a.vz),
                                        value_t(a.mx), value_t(a.my), value_t(a.mz)},
                                       {value_t(b.vx), value_t(b.vy), value_t(b.vz),
                                        value_t(b.mx), value_t(b.my), value_t(b.mz)},
                                       tolerance);
}

// For trivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(TriVec2dc<T> const& a, TriVec2dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, tolerance);
}

// For pseudoscalars: all non-zero pseudoscalars represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar2dc<T> a, PScalar2dc<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars are congruent (represent the full 4D space)
    return true;
}


////////////////////////////////////////////////////////////////////////////////
// is_close(): equality within a RELATIVE tolerance
////////////////////////////////////////////////////////////////////////////////
//
// Same question as operator==, but with the threshold scaled by the operands (see
// detail::coeffs_close). Reach for it whenever the values carry a physical scale:
// operator== measures against an absolute eps, which cannot resolve anything once
// coordinates grow large, where a single ulp already exceeds it. In a conformal
// algebra that regime is the default: the round (u) component of an embedded
// point grows QUADRATICALLY with the coordinates.
//
// Distinct from is_congruent(), which allows an arbitrary scale factor between the
// operands and so answers "same subspace", not "same value".

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Scalar2dc<T> a, Scalar2dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Vec2dc<T> const& a, Vec2dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(BiVec2dc<T> const& a, BiVec2dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<6>({value_t(a.vx), value_t(a.vy), value_t(a.vz),
                                    value_t(a.mx), value_t(a.my), value_t(a.mz)},
                                   {value_t(b.vx), value_t(b.vy), value_t(b.vz),
                                    value_t(b.mx), value_t(b.my), value_t(b.mz)},
                                   rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(TriVec2dc<T> const& a, TriVec2dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(PScalar2dc<T> a, PScalar2dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec2dc_E<T> const& a, MVec2dc_E<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<8>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec2dc_U<T> const& a, MVec2dc_U<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<8>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec2dc<T> const& a, MVec2dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<16>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7), value_t(a.c8), value_t(a.c9),
         value_t(a.c10), value_t(a.c11), value_t(a.c12), value_t(a.c13), value_t(a.c14),
         value_t(a.c15)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7), value_t(b.c8), value_t(b.c9),
         value_t(b.c10), value_t(b.c11), value_t(b.c12), value_t(b.c13), value_t(b.c14),
         value_t(b.c15)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(DualNum2dc<T> const& a, DualNum2dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<2>({value_t(a.c0), value_t(a.c1)},
                                   {value_t(b.c0), value_t(b.c1)}, rel_tol);
}


// Do two motors describe the same conformal transformation? NOT the same
// question as M1 == M2: versors double-cover the transformations (M and -M
// act identically), and any uniform rescaling of a versor changes only the
// homogeneous weights of the images, not the geometry. Decided by ACTION: a
// Moebius transformation is pinned by the images of three points, so the
// images of three spread round points are compared up to weight.
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_same_transform(MVec2dc_E<T> const& M1, MVec2dc_E<U> const& M2,
                       value_t rel_tol = eps_congruent)
{
    auto const p0 = Vec2dc<value_t>(0.0, 0.0, 1.0, 0.0); // embedded origin
    auto const p1 = Vec2dc<value_t>(1.0, 0.0, 1.0, 0.5); // embedded (1,0)
    auto const p2 = Vec2dc<value_t>(0.0, 1.0, 1.0, 0.5); // embedded (0,1)
    return is_congruent(transform(p0, M1), transform(p0, M2), rel_tol) &&
           is_congruent(transform(p1, M1), transform(p1, M2), rel_tol) &&
           is_congruent(transform(p2, M1), transform(p2, M2), rel_tol);
}

} // namespace hd::ga::cga
