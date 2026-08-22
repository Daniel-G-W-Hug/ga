#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_cga3dc_ops_basics.hpp"
#include "ga_cga3dc_ops_products.hpp"

#include <complex> // exp/log/sqrt of regressive versors (central subalgebra = C)
#include <limits>  // get_translation zero-direction guard


namespace hd::ga::cga {

////////////////////////////////////////////////////////////////////////////////
// provides cga3dc geometric operations (layer under construction):
//
// - rexp()                  -> versor from a simple trivector generator
//                              (w.r.t. rgpr, closed form; throws non-simple)
// - rlog()                  -> principal generator of a simple-image motor
// - rsqrt()                 -> square root of such a motor
// - get_translation()       -> motor translating by (dx, dy, dz)
// - get_rotation()          -> motor rotating about a line by theta
// - get_dilation()          -> motor scaling by sigma about a point
// - get_transversion()      -> special conformal transformation motor
// - get_loxodromic()        -> two-fixed-point motor from a dipole
// - transform()             -> apply motor: sandwich M (v) u (v) rrev(M)
// - invert_on()             -> inversion in a sphere or plane (flector
//                              sandwich)
//
// Object construction (grade encodes the object kind: round point = vector,
// dipole / flat point = bivector, circle / line = trivector, sphere / plane =
// quadvector):
//
// - round_point3dc()        -> round point from center and radius (r = 0: the
//                              null-point embedding of a Euclidean point)
// - flat_point3dc()         -> flat point from a Euclidean position
// - sphere3dc()             -> sphere from center and radius
// - circle3dc()             -> circle from center, radius and unit normal
// - dipole3dc()             -> dipole from center, radius and unit axis
// - line3dc()               -> line through a point with a direction
// - plane3dc()              -> plane from a unit normal and origin distance
//
// (the bulk/weight split, the part/center norms, unitize and cconj live in
// ga_cga3dc_ops_basics.hpp with the other algebras' norm layer; the radius
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
// - con()                   -> container: the smallest sphere containing u
// - par()                   -> partner: same center/carrier, r^2 negated
// - att()                   -> attitude (direction content; removes the origin)
//
// - is_congruent()          -> same subspace up to a non-zero scalar factor
// - is_close()              -> same value within a RELATIVE tolerance
// - is_same_transform()     -> do two motors act as the same conformal map?
////////////////////////////////////////////////////////////////////////////////


// forward declarations of the object constructors used by the transformation
// builders (defined in the object-construction section below, which follows
// the canonical per-algebra ordering: transformations first, objects after)
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> flat_point3dc(T px, T py, T pz);
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> line3dc(T px, T py, T pz, T vx, T vy, T vz);
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> plane3dc(T nx, T ny, T nz, T d);
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> car(BiVec3dc<T> const& B);


////////////////////////////////////////////////////////////////////////////////
// exponential, logarithm and square root w.r.t. the regressive geometric
// product -- the versor side of the algebra: conformal transformations are
// applied as rgpr sandwiches, so their generators and roots live here.
//
// In the odd-dimensional cga3dc the regressive versors (motors) are ODD
// multivectors (the rgpr identity is the odd pseudoscalar I), and the
// exponential generators are TRIVECTORS (antigrade 2). The regressive square
// of a generator t is
//
//     rgpr(t, t) = alpha * I + rwdg(t, t)
//
// where the grade-1 part rwdg(t, t) measures NON-SIMPLICITY. For a SIMPLE
// generator (a blade: every line, circle and dual flat point is one) the
// square is purely central and the series sums in closed form,
//
//     rexp(t) = cosh(w) * I + (sinh(w)/w) * t,   w = rsqrt(alpha)
//
// evaluated with complex arithmetic so the elliptic (alpha < 0), hyperbolic
// (alpha > 0) and parabolic (alpha = 0, translations: rexp(t) = I + t)
// branches come out of ONE formula (cosh(w) and sinh(w)/w are even in w, so
// the square-root branch is irrelevant). Non-simple generators (screws,
// loxodromics) are NOT covered yet -- exp throws for them; compose their
// motors as rgpr products of simple ones instead.
//
// rlog() returns the principal generator: complex acosh of the pseudoscalar
// part, then t = V / (sinh(w)/w) with V the trivector part. It is restricted
// to motors in the image of the simple exp (vanishing vector part) and throws
// otherwise, and at the true singularity (full-turn versor, central part -I).
// rsqrt(M) = rexp(rlog(M)/2).
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> rexp(TriVec3dc<T> const& t)
{
    auto const sq = rgpr(t, t); // alpha * I + rwdg(t,t), central iff simple
    auto const nsp = gr1(sq);   // non-simplicity (vector) part
    T const alpha = T(gr5(sq));
    // component norm (the metric dot is not definite on the null directions)
    T const nsp_nrm = std::sqrt(nsp.x * nsp.x + nsp.y * nsp.y + nsp.z * nsp.z +
                                nsp.w * nsp.w + nsp.u * nsp.u);
    if (nsp_nrm >
        T(1e3) * std::numeric_limits<T>::epsilon() * (std::abs(alpha) + T(1.0))) {
        throw std::invalid_argument(
            "cga3dc exp: non-simple generator (rwdg(t,t) != 0) not supported -- "
            "compose motors of simple generators via rgpr instead");
    }
    std::complex<T> const z(alpha, T(0.0));
    std::complex<T> const w = std::sqrt(z);
    T const ch = std::cosh(w).real();
    T const sh_ow = (std::abs(alpha) > T(1e-4))
                        ? (std::sinh(w) / w).real()
                        : T(1.0) + alpha / T(6.0) + alpha * alpha / T(120.0);
    return MVec3dc_U<T>(Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)), sh_ow * t,
                        PScalar3dc<T>(ch));
}

template <typename T>
    requires(numeric_type<T>)
inline TriVec3dc<T> rlog(MVec3dc_U<T> const& M)
{
    // restricted to the image of the simple exp: vector part must vanish
    auto const v = gr1(M);
    T const v_nrm = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w + v.u * v.u);
    if (v_nrm >
        T(1e3) * std::numeric_limits<T>::epsilon() * (std::abs(T(gr5(M))) + T(1.0))) {
        throw std::invalid_argument(
            "cga3dc log: motor outside the simple-generator image (non-zero "
            "vector part)");
    }
    T const c = T(gr5(M));
    TriVec3dc<T> const V = gr3(M); // V = (sinh(w)/w) * t
    std::complex<T> const w = std::acosh(std::complex<T>(c, T(0.0)));
    std::complex<T> const w2 = w * w;
    std::complex<T> const sh_ow = (std::abs(w) > T(1e-2))
                                      ? std::sinh(w) / w
                                      : T(1.0) + w2 / T(6.0) + w2 * w2 / T(120.0);
    // singular exactly at the full-turn versors (central part -I)
    hd::ga::detail::check_normalization<T>(std::abs(sh_ow),
                                           "regressive logarithm (full-turn versor)");
    return T(1.0) / sh_ow.real() * V;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> rsqrt(MVec3dc_U<T> const& M)
{
    return rexp(T(0.5) * rlog(M));
}


////////////////////////////////////////////////////////////////////////////////
// conformal transformations
//
// An odd unit versor (motor) M transforms any object with the regressive
// sandwich
//
//     u' = transform(u, M) = M (v) u (v) rrev(M)
//
// Reflections generalize to INVERSIONS in spheres (and planes, which are
// spheres through infinity): an anti-odd sandwich with the quadvector s,
//
//     u' = invert_on(u, s) = -( s (v) u (v) rrev(s) )
//
// (the leading minus is the flector-sandwich convention of the reference).
// Composing two inversions gives a motor: parallel planes -> translation,
// intersecting planes -> rotation about their meet line, concentric spheres
// -> dilation about their center.
//
// The named builders produce unit motors from the exponential (translation
// directly in the terminated parabolic form):
//
//     translation by d:            M = I + B,  B from the parabolic generator
//     rotation about a line by theta:  M = rexp( (theta/2) * unit line )
//     dilation about m by sigma:   M = rexp( delta * dual flat point at m ),
//                                  delta = -1/2 ln(sigma)
//
// Weights under sandwiches are homogeneous: a non-unit versor scales all
// image weights uniformly; normalize results with unitize() where a
// weight-one representative is needed.
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> transform(Vec3dc<T> const& v,
                                                     MVec3dc_U<U> const& M)
{
    return gr1(rgpr(rgpr(M, v), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> transform(BiVec3dc<T> const& B,
                                                       MVec3dc_U<U> const& M)
{
    return gr2(rgpr(rgpr(M, B), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> transform(TriVec3dc<T> const& t,
                                                        MVec3dc_U<U> const& M)
{
    return gr3(rgpr(rgpr(M, t), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> transform(QuadVec3dc<T> const& Q,
                                                         MVec3dc_U<U> const& M)
{
    return gr4(rgpr(rgpr(M, Q), rrev(M)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr Vec3dc<std::common_type_t<T, U>> invert_on(Vec3dc<T> const& v,
                                                     QuadVec3dc<U> const& s)
{
    return -gr1(rgpr(rgpr(s, v), rrev(s)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BiVec3dc<std::common_type_t<T, U>> invert_on(BiVec3dc<T> const& B,
                                                       QuadVec3dc<U> const& s)
{
    return -gr2(rgpr(rgpr(s, B), rrev(s)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr TriVec3dc<std::common_type_t<T, U>> invert_on(TriVec3dc<T> const& t,
                                                        QuadVec3dc<U> const& s)
{
    return -gr3(rgpr(rgpr(s, t), rrev(s)));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr QuadVec3dc<std::common_type_t<T, U>> invert_on(QuadVec3dc<T> const& Q,
                                                         QuadVec3dc<U> const& s)
{
    return -gr4(rgpr(rgpr(s, Q), rrev(s)));
}

// translation by (dx, dy, dz): parabolic, the exponential series terminates
// (generator = half the line at infinity dual to the direction; SIGN pinned
// by the tests)
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> get_translation(T dx, T dy, T dz)
{
    // exact two-plane form: reflection on the plane through the origin
    // followed by reflection on the parallel plane at half the distance
    T const n2 = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (n2 < std::numeric_limits<T>::epsilon()) {
        return MVec3dc_U<T>(Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)),
                            TriVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                                         T(0.0), T(0.0), T(0.0), T(0.0)),
                            PScalar3dc<T>(T(1.0)));
    }
    T const nx = dx / n2, ny = dy / n2, nz = dz / n2;
    auto const M = rgpr(plane3dc(nx, ny, nz, T(0.5) * n2), plane3dc(nx, ny, nz, T(0.0)));
    return gr1(M) + gr3(M) + gr5(M); // already grades {3, 5}; assembled as mv_u
}

// rotation about the line through (px, py, pz) with direction (vx, vy, vz)
// by theta (right-handed about the direction for theta > 0); the generator is
// the unit line
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> get_rotation(T px, T py, T pz, T vx, T vy, T vz, T theta)
{
    auto const l = line3dc(px, py, pz, vx, vy, vz);
    auto const sq = rgpr(l, l); // = alpha * I, alpha = -|v|^2 (elliptic)
    T const alpha = T(gr5(sq));
    hd::ga::detail::check_normalization<T>(std::abs(alpha), "rotation axis line");
    auto const lu = T(1.0) / std::sqrt(-alpha) * l;
    return rexp(TriVec3dc<T>(T(0.5) * theta * lu));
}

// dilation about (mx, my, mz) by the scale factor sigma > 0; the generator is
// the ANTIDUAL of the flat point at m (hyperbolic, regressive square +I; the
// antidual carries the star role in this library, as for the properties)
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> get_dilation(T mx, T my, T mz, T sigma)
{
    hd::ga::detail::check_normalization<T>(sigma, "dilation scale factor");
    T const delta = -T(0.5) * std::log(sigma);
    auto const g = antidual(flat_point3dc(mx, my, mz));
    return rexp(TriVec3dc<T>(delta * g));
}


// transversion (special conformal transformation) by (tx, ty, tz): the
// conjugate of a translation by the unit-sphere inversion,
//
//     x' = (x + x^2 t) / (1 + 2 t.x + t^2 x^2)
//
// built directly as the versor product s (v) T (v) s with the unit sphere s
// at the origin (an even flector conjugating the odd translation motor back
// to a motor)
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> get_transversion(T tx, T ty, T tz)
{
    auto const s = QuadVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.5), T(-1.0)); // unit sphere
    auto const Tm = get_translation(tx, ty, tz);
    auto const M = rgpr(rgpr(s, Tm), s);
    return gr1(M) + gr3(M) + gr5(M);
}

// loxodromic motor about the two surface points of the (unitized) dipole d:
// rotation by phi about the carrier line combined with the hyperbolic motion
// of rapidity delta that flows from one dipole point to the other (both fix
// the two points; delta > 0 flows towards the point at center + r * axis). The generator
// phi/2 * l + delta * h is NON-simple, but its two parts commute, so the motor is the
// product of the two SIMPLE exponentials:
//
//     M = rexp( (phi/2) * l ) (v) rexp( delta * h )
//
// with l the unit carrier line, h = antidual(d) normalized to rgpr(h,h) = +I.
template <typename T>
    requires(numeric_type<T>)
inline MVec3dc_U<T> get_loxodromic(BiVec3dc<T> const& d, T phi, T delta)
{
    auto const du = unitize(d);
    // elliptic part: the unit carrier line of the dipole
    auto const l = car(du);
    auto const sql = rgpr(l, l);
    T const al = T(gr5(sql));
    hd::ga::detail::check_normalization<T>(std::abs(al), "dipole carrier line");
    auto const lu = T(1.0) / std::sqrt(-al) * l;
    // hyperbolic part: the antidual of the dipole (fixes the same two points)
    auto const h = antidual(du);
    auto const sqh = rgpr(h, h);
    T const ah = T(gr5(sqh));
    hd::ga::detail::check_normalization<T>(std::abs(ah), "dipole antidual");
    auto const hu = T(1.0) / std::sqrt(ah) * h;
    return rgpr(rexp(TriVec3dc<T>(T(0.5) * phi * lu)), rexp(TriVec3dc<T>(-delta * hu)));
}


////////////////////////////////////////////////////////////////////////////////
// geometric object construction
//
// The round objects and their component semantics (the object kind is
// encoded in the grade): round point (vector), dipole (bivector), circle
// (trivector), sphere (quadvector); the flats are the rounds through the
// point at infinity: flat point (bivector), line (trivector), plane
// (quadvector).
//
// The derived constructors are product chains over round points (wedge =
// join, rwdg = meet), so their component forms follow from the generated
// products; the orientation conventions are pinned in the tests.
////////////////////////////////////////////////////////////////////////////////

// round point with center (px, py, pz) and radius r (unitized: origin weight
// 1); r = 0 gives the null-point embedding q with q * q = 0
template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> round_point3dc(T px, T py, T pz, T r)
{
    return Vec3dc<T>(px, py, pz, T(1.0), T(0.5) * (px * px + py * py + pz * pz + r * r));
}

// flat point at Euclidean position (px, py, pz) (unitized), reference (4.25):
// p = px e15 + py e25 + pz e35 + e45 -- POSITIVE components (the 3d reference
// blade order matches this library directly, unlike the 2d case where the
// origin-leading reference blades flip the sign); equals the carrier of the
// round point at the same position, car(round_point3dc(p, r)) == flat_point3dc(p)
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> flat_point3dc(T px, T py, T pz)
{
    return BiVec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), px, py, pz,
                       T(1.0));
}

// sphere with center (px, py, pz) and radius r (unitized, round weight -1 as
// for circle2dc): the codim-1 round object is the negated antidual of the
// grade-1 carrier with the (p^2 - r^2)/2 infinity coefficient (r^2 enters the
// codim-1 form with the OPPOSITE sign vs. the round point, as in circle2dc);
// the orientation matches the container: sphere3dc(p, r) == con(round point)
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> sphere3dc(T px, T py, T pz, T r)
{
    return -antidual(
        Vec3dc<T>(px, py, pz, T(1.0), T(0.5) * (px * px + py * py + pz * pz - r * r)));
}

// line through (px, py, pz) with direction (vx, vy, vz): the join of two round
// points on it with the point at infinity
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> line3dc(T px, T py, T pz, T vx, T vy, T vz)
{
    auto const a1 = round_point3dc(px, py, pz, T(0.0));
    auto const a2 = round_point3dc(px + vx, py + vy, pz + vz, T(0.0));
    return wdg(wdg(a1, a2), Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(1.0)));
}

// plane with unit normal (nx, ny, nz) and signed origin distance d: the
// carrier of any circle in it -- built as the flat of a great circle
template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> plane3dc(T nx, T ny, T nz, T d)
{
    // three points spanning the plane: p0 = d*n and p0 +/- two orthogonal
    // in-plane directions; the join of their circle with infinity is the plane
    // (normalized in the tests; orientation pinned there)
    T const ax = (std::abs(nx) < T(0.9)) ? T(1.0) : T(0.0);
    T const ay = (std::abs(nx) < T(0.9)) ? T(0.0) : T(1.0);
    // u1 = a x n, u2 = n x u1 (both in-plane)
    T const u1x = ay * nz;
    T const u1y = -ax * nz;
    T const u1z = ax * ny - ay * nx;
    T const u2x = ny * u1z - nz * u1y;
    T const u2y = nz * u1x - nx * u1z;
    T const u2z = nx * u1y - ny * u1x;
    auto const p0 = round_point3dc(d * nx, d * ny, d * nz, T(0.0));
    auto const p1 = round_point3dc(d * nx + u1x, d * ny + u1y, d * nz + u1z, T(0.0));
    auto const p2 = round_point3dc(d * nx + u2x, d * ny + u2y, d * nz + u2z, T(0.0));
    return wdg(wdg(wdg(p0, p1), p2), Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(1.0)));
}

// circle with center (px, py, pz), radius r and unit normal (nx, ny, nz): the
// meet of the sphere with the plane through the center; negated so that
// att() recovers the requested normal with POSITIVE sign
template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> circle3dc(T px, T py, T pz, T r, T nx, T ny, T nz)
{
    T const d = px * nx + py * ny + pz * nz;
    return -rwdg(sphere3dc(px, py, pz, r), plane3dc(nx, ny, nz, d));
}

// dipole with center (px, py, pz), radius r and unit axis (nx, ny, nz) (the
// two points lie at center +/- r * axis): the meet of the sphere with the
// line through the center ALONG the axis (in 3D the normal points along the
// carrier line, unlike the 2D convention); negated so that att() recovers the
// requested axis with POSITIVE sign
template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> dipole3dc(T px, T py, T pz, T r, T nx, T ny, T nz)
{
    return -rwdg(sphere3dc(px, py, pz, r), line3dc(px, py, pz, nx, ny, nz));
}


////////////////////////////////////////////////////////////////////////////////
// radius norm and squared radius
//
// The radius norm square is the ANTIDOT square u (o) u = rdot(u, u) = r^2 (for
// a unitized object; the dot square gives dot(u, u) = -r^2). It is SIGNED:
// negative values mean an imaginary radius (e.g. the meet of two
// non-intersecting spheres), so no square root is taken here. These live in
// the geometric layer (not in ops_basics with the other norms) because rdot is
// a product-layer operation.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr T radius_nrm_sq(Vec3dc<T> const& v)
{
    return T(rdot(v, v));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T radius_nrm_sq(BiVec3dc<T> const& B)
{
    return T(rdot(B, B));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T radius_nrm_sq(TriVec3dc<T> const& t)
{
    return T(rdot(t, t));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T radius_nrm_sq(QuadVec3dc<T> const& Q)
{
    return T(rdot(Q, Q));
}

// squared radius of a round object of any weight (throws for flat objects,
// whose round weight is zero)
// HINT: radius_sq()/unitize() cannot be constexpr due to the checks for
//       division by zero which might throw

template <typename T>
    requires(numeric_type<T>)
inline T radius_sq(Vec3dc<T> const& v)
{
    T const wn_sq = round_weight_nrm_sq(v);
    hd::ga::detail::check_normalization<T>(wn_sq, "round point (round weight)");
    return radius_nrm_sq(v) / wn_sq;
}

template <typename T>
    requires(numeric_type<T>)
inline T radius_sq(BiVec3dc<T> const& B)
{
    T const wn_sq = round_weight_nrm_sq(B);
    hd::ga::detail::check_normalization<T>(wn_sq, "dipole (round weight)");
    return radius_nrm_sq(B) / wn_sq;
}

template <typename T>
    requires(numeric_type<T>)
inline T radius_sq(TriVec3dc<T> const& t)
{
    T const wn_sq = round_weight_nrm_sq(t);
    hd::ga::detail::check_normalization<T>(wn_sq, "circle (round weight)");
    return radius_nrm_sq(t) / wn_sq;
}

template <typename T>
    requires(numeric_type<T>)
inline T radius_sq(QuadVec3dc<T> const& Q)
{
    T const wn_sq = round_weight_nrm_sq(Q);
    hd::ga::detail::check_normalization<T>(wn_sq, "sphere (round weight)");
    return radius_nrm_sq(Q) / wn_sq;
}

// scale a round object so its round weight norm becomes one (sign preserved;
// throws for flat objects)


////////////////////////////////////////////////////////////////////////////////
// geometric properties of the round objects
//
// All are defined coordinate-free over the products already in place (the
// dimension-free defining equations of the reference; the odd-dimensional
// cga3dc has the singular antidual in the star role):
//
//     car(u) = u ^ e5                    (grade k -> k+1; annihilates the flat
//                                         part, lifts the round part)
//     ccr(u) = antidual(u) ^ e5          (= car of the antidual)
//     cen(u) = rwdg( ccr(u), u )         (center as a ROUND POINT: position
//                                         AND radius of u)
//     con(u) = u ^ antidual( car(u) )    (smallest sphere containing u)
//     par(u) = (-1)^(gr(u)+1) rwdg( con(antidual(u)), car(u) )
//                                        (same center + carrier, r^2 negated;
//                                         u . par(u) = 0 always)
//     att(u) = rwdg( u, cmpl(e4) )       (removes the origin content; for a
//                                         dipole: carrier direction + weight)
//
// The wedge with e5 is taken on the RIGHT, matching this library's trailing
// infinity-blade ordering (the same convention the reference uses for its
// own trailing ordering).
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> car(Vec3dc<T> const& v)
{
    return wdg(v, Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> car(BiVec3dc<T> const& B)
{
    return wdg(B, Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> car(TriVec3dc<T> const& t)
{
    return wdg(t, Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> car(QuadVec3dc<T> const& Q)
{
    return wdg(Q, Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(0.0), T(1.0)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> ccr(Vec3dc<T> const& v)
{
    return car(antidual(v));
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> ccr(BiVec3dc<T> const& B)
{
    return car(antidual(B));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> ccr(TriVec3dc<T> const& t)
{
    return car(antidual(t));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> ccr(QuadVec3dc<T> const& Q)
{
    return car(antidual(Q));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> cen(Vec3dc<T> const& v)
{
    return rwdg(ccr(v), v);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> cen(BiVec3dc<T> const& B)
{
    return rwdg(ccr(B), B);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> cen(TriVec3dc<T> const& t)
{
    return rwdg(ccr(t), t);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> cen(QuadVec3dc<T> const& Q)
{
    return rwdg(ccr(Q), Q);
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> con(Vec3dc<T> const& v)
{
    return wdg(v, antidual(car(v)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> con(BiVec3dc<T> const& B)
{
    return wdg(B, antidual(car(B)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> con(TriVec3dc<T> const& t)
{
    return wdg(t, antidual(car(t)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> con(QuadVec3dc<T> const& Q)
{
    return wdg(Q, antidual(car(Q)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> par(Vec3dc<T> const& v)
{
    // gr = 1: sign +1
    return rwdg(con(antidual(v)), car(v));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> par(BiVec3dc<T> const& B)
{
    // gr = 2: sign -1
    return T(-1.0) * rwdg(con(antidual(B)), car(B));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> par(TriVec3dc<T> const& t)
{
    // gr = 3: sign +1
    return rwdg(con(antidual(t)), car(t));
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> par(QuadVec3dc<T> const& Q)
{
    // gr = 4: sign -1
    return T(-1.0) * rwdg(con(antidual(Q)), car(Q));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> att(Vec3dc<T> const& v)
{
    return rwdg(v, cmpl(Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(1.0), T(0.0))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> att(BiVec3dc<T> const& B)
{
    return rwdg(B, cmpl(Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(1.0), T(0.0))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> att(TriVec3dc<T> const& t)
{
    return rwdg(t, cmpl(Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(1.0), T(0.0))));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> att(QuadVec3dc<T> const& Q)
{
    return rwdg(Q, cmpl(Vec3dc<T>(T(0.0), T(0.0), T(0.0), T(1.0), T(0.0))));
}


////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar3dc<T> a, Scalar3dc<U> b, value_t tolerance = eps)
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
bool is_congruent(Vec3dc<T> const& a, Vec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)},
        tolerance);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec3dc<T> const& a, BiVec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        tolerance);
}

// For trivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(TriVec3dc<T> const& a, TriVec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        tolerance);
}

// For quadvectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(QuadVec3dc<T> const& a, QuadVec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)},
        tolerance);
}

// For pseudoscalars: all non-zero pseudoscalars represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar3dc<T> a, PScalar3dc<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars are congruent
    return true;
}


////////////////////////////////////////////////////////////////////////////////
// is_close(): equality within a RELATIVE tolerance
//
// operator== compares against an absolute eps and cannot resolve anything once
// coordinates are large; is_close scales the per-component threshold by the
// operands (see detail::coeffs_close). Reach for it whenever the values carry a
// physical scale. Distinct from is_congruent(), which allows an arbitrary scale
// factor between the operands.
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Scalar3dc<T> a, Scalar3dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Vec3dc<T> const& a, Vec3dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(BiVec3dc<T> const& a, BiVec3dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(TriVec3dc<T> const& a, TriVec3dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(QuadVec3dc<T> const& a, QuadVec3dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(PScalar3dc<T> a, PScalar3dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec3dc_E<T> const& a, MVec3dc_E<U> const& b,
              value_t rel_tol = eps_congruent)
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
bool is_close(MVec3dc_U<T> const& a, MVec3dc_U<U> const& b,
              value_t rel_tol = eps_congruent)
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
bool is_close(MVec3dc<T> const& a, MVec3dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<32>(
        {value_t(a.c0),  value_t(a.c1),  value_t(a.c2),  value_t(a.c3),  value_t(a.c4),
         value_t(a.c5),  value_t(a.c6),  value_t(a.c7),  value_t(a.c8),  value_t(a.c9),
         value_t(a.c10), value_t(a.c11), value_t(a.c12), value_t(a.c13), value_t(a.c14),
         value_t(a.c15), value_t(a.c16), value_t(a.c17), value_t(a.c18), value_t(a.c19),
         value_t(a.c20), value_t(a.c21), value_t(a.c22), value_t(a.c23), value_t(a.c24),
         value_t(a.c25), value_t(a.c26), value_t(a.c27), value_t(a.c28), value_t(a.c29),
         value_t(a.c30), value_t(a.c31)},
        {value_t(b.c0),  value_t(b.c1),  value_t(b.c2),  value_t(b.c3),  value_t(b.c4),
         value_t(b.c5),  value_t(b.c6),  value_t(b.c7),  value_t(b.c8),  value_t(b.c9),
         value_t(b.c10), value_t(b.c11), value_t(b.c12), value_t(b.c13), value_t(b.c14),
         value_t(b.c15), value_t(b.c16), value_t(b.c17), value_t(b.c18), value_t(b.c19),
         value_t(b.c20), value_t(b.c21), value_t(b.c22), value_t(b.c23), value_t(b.c24),
         value_t(b.c25), value_t(b.c26), value_t(b.c27), value_t(b.c28), value_t(b.c29),
         value_t(b.c30), value_t(b.c31)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(DualNum3dc<T> const& a, DualNum3dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<2>({value_t(a.c0), value_t(a.c1)},
                                   {value_t(b.c0), value_t(b.c1)}, rel_tol);
}

// Do two motors describe the same conformal transformation? NOT the same
// question as M1 == M2: versors double-cover the transformations (M and -M
// act identically), and any uniform rescaling of a versor changes only the
// homogeneous weights of the images, not the geometry. Decided by ACTION: a
// 3d Moebius transformation is pinned by the images of a spread point set,
// so the images of five embedded points in general position are compared up
// to weight.
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_same_transform(MVec3dc_U<T> const& M1, MVec3dc_U<U> const& M2,
                       value_t rel_tol = eps_congruent)
{
    auto const p0 = Vec3dc<value_t>(0.0, 0.0, 0.0, 1.0, 0.0); // embedded origin
    auto const p1 = Vec3dc<value_t>(1.0, 0.0, 0.0, 1.0, 0.5); // embedded (1,0,0)
    auto const p2 = Vec3dc<value_t>(0.0, 1.0, 0.0, 1.0, 0.5); // embedded (0,1,0)
    auto const p3 = Vec3dc<value_t>(0.0, 0.0, 1.0, 1.0, 0.5); // embedded (0,0,1)
    auto const p4 = Vec3dc<value_t>(1.0, 1.0, 1.0, 1.0, 1.5); // embedded (1,1,1)
    return is_congruent(transform(p0, M1), transform(p0, M2), rel_tol) &&
           is_congruent(transform(p1, M1), transform(p1, M2), rel_tol) &&
           is_congruent(transform(p2, M1), transform(p2, M2), rel_tol) &&
           is_congruent(transform(p3, M1), transform(p3, M2), rel_tol) &&
           is_congruent(transform(p4, M1), transform(p4, M2), rel_tol);
}

} // namespace hd::ga::cga
