#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <cmath>   // std::abs, std::cos, std::log, std::sin, std::sqrt
#include <cstddef> // size_t
#include <span>    // std::span
#include <utility> // std::pair
#include <vector>  // std::vector

#include "detail/ga_quadrature.hpp" // gauss_legendre12 (the closed-surface rule)
#include "ega/ga_ega2d_ops.hpp"     // the 2D field algebra
#include "ega/ga_ega3d_ops.hpp"     // the 3D field algebra
#include "ga_usr_types.hpp"         // vec2d, vec3d, bivec2dp, bivec3dp, ...
#include "ga_usr_utilities.hpp"     // pi
#include "ga_value_t.hpp"           // value_t
#include "pga/ga_pga2dp_ops.hpp"    // the 2D wrench: wdg of a point and a direction
#include "pga/ga_pga3dp_ops.hpp"    // the 3D wrench

/////////////////////////////////////////////////////////////////////////////////////////
// Magnetostatics in GA, in two and three dimensions: the field of the elementary sources,
// the force on a current, the Maxwell stress as one sandwich, and the wrench of a field
// on what a closed contour or surface encloses. The FIELD quantities live in EGA (where
// a value is evaluated is not part of what it is); the WRENCH lives in PGA (a force's
// line of action is part of what it is), the EGA -> PGA lift being the point Vec(P, 1)
// and the direction Vec(t, 0).
//
// THE GRADES. The flux density is carried as a VECTOR b in both dimensions. In 2D the
// field of a z-independent current distribution lies in the plane, and the out-of-plane
// quantities are pseudoscalars: a line current I = I_z e12, the vector potential
// A = A_z e12, a magnetisation's surface current K = m ^ n. In 3D a current element is
// the vector I dl, the vector potential a vector, a dipole moment a vector mu; the
// induction BIVECTOR B = I3 b = dual(b) is the other common convention, and each formula
// below says whether it changes sign under it.
//
// THE 2D KERNEL. The field of a line current is the vector inverse,
//
//     b(r) = mu0/(2 pi) (r - p)^-1 I         ((r - p)^-1 = (r - p) / |r - p|^2)
//
// the 2D Green's function of the vector derivative (nabla b = mu0 J, one equation for
// div b = 0 and curl b = mu0 J_z), with b = -nabla A and A = -(mu0/(2 pi)) ln|r - p| I.
// The flux per unit length through a path from P to Q is A(Q) - A(P).
//
// THE 3D KERNEL. A straight segment from a to b carrying I (unit direction u, rho the
// rejection of r - a from the segment's line) has the field
//
//     b(r) = (mu0 I / (4 pi)) (c_a - c_b) cmpl(u ^ rho^-1)
//     c_a = u . (r - a) / |r - a|,   c_b = u . (r - b) / |r - b|
//
// and the 2D kernel is its limit for an infinitely long segment (c_a - c_b -> 2). Its
// vector potential lies along the segment:
//
//     A(r) = (mu0 I / (4 pi)) ln((|r - b| - u.(r - b)) / (|r - a| - u.(r - a))) u
//
// THE DIPOLE. A small loop of moment mu = I (area) n has, far away, the field
//
//     b(r) = mu0 / (8 pi |r|^3) (mu + 3 r_hat mu r_hat)
//
// Because I3 is central in 3D, the bivector form is the SAME expression with mu ->
// dual(mu) and b -> B: no sign change. The coupling of a second moment to that field is
// the scalar <m B> (the energy is its negative), the torque on it cmt(B, m). In 2D the
// line dipole's field is the pure sandwich (mu0 / (2 pi |r|^2)) r_hat mu r_hat.
//
// THE STRESS. The traction on a surface with unit normal n is one sandwich in either
// dimension,
//
//     t(n) = (1/(2 mu0)) b n b                (its grade-1 part in 3D)
//
// equal to Maxwell's tensor (1/mu0)((b.n) b - (1/2) b^2 n) because b n b = 2 (b.n) b -
// b^2 n for vectors. SIGN: with the bivector B = I3 b the same traction is
// -(1/(2 mu0)) B n B, since I3^2 = -1.
//
// THE WRENCH. The force and torque of the field on everything inside a closed contour
// (2D) or surface (3D) is the sum of the tractions on their points of action,
//
//     W = sum_i wdg(P_i, t_i) ds_i    (2D: dA_i in 3D)
//
// the same expression in both dimensions. In 2D W is a line of action; in 3D the force
// lines through different points sum to a SCREW, so a couple (a dipole's torque) appears
// with no couple term written anywhere.
//
// provides in namespace hd::ga:
//
// - mu0                       : the vacuum permeability [Vs/Am]
//
// 2D (ega2d field, pga2dp wrench):
//
// - line_current              : a z-directed line current at a point
// - b_line_current()          : its field (the 2D kernel)
// - A_line_current()          : its vector potential, -(mu0/(2 pi)) ln|r - p| I
// - b_line_currents()         : the field of a set of line currents
// - lorentz_force(vec2d, pscalar2d) : the force per length on a line current, b I
// - traction(vec2d, vec2d)    : the stress sandwich t(n)
// - contour_point             : a quadrature point of a closed contour (P, outward n, ds)
// - circle_contour()          : the midpoint rule on a circle
// - stress_wrench(b, contour) : W = sum wdg(P, t(n)) ds, a bivec2dp
// - surface_current(vec2d, vec2d)   : K = m ^ n, a magnetisation's surface current
// - magnetised_disc_currents(): a uniformly magnetised disc as boundary line currents
//
// 3D (ega3d field, pga3dp wrench):
//
// - current_segment           : a straight current from a to b
// - b_segment(), A_segment()  : its field and vector potential (Biot-Savart, closed)
// - b_segments()              : the field of a set of segments
// - circle_loop()             : a circular loop as N chords
// - lorentz_force(vec3d, vec3d) : the force on a current element I dl, cmpl(I dl ^ b)
// - b_dipole(r, vec3d mu) / b_dipole(r, bivec3d m) : a point dipole's field, vector or
//                               bivector convention
// - dipole_coupling()         : <m_j B_i(r)>, the second moment's coupling to the first's
//                               field (bivector convention)
// - traction(vec3d, vec3d)    : the stress sandwich t(n)
// - surface_point             : a quadrature point of a closed surface (P, n, dA)
// - sphere_surface()          : Gauss-Legendre in cos(theta) x the midpoint rule in phi
// - stress_wrench(b, surface) : W = sum wdg(P, t(n)) dA, a bivec3dp
// - surface_current(vec3d, vec3d)   : K = m x n = cmpl(m ^ n)
// - magnetised_sphere_currents(): a uniformly magnetised sphere as current rings
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

inline value_t const mu0 = 4.0e-7 * pi; // vacuum permeability [Vs/Am]

/////////////////////////////////////////////////////////////////////////////////////////
// 2D
/////////////////////////////////////////////////////////////////////////////////////////

// a line current I = I_z e12 [A] along z through the point p
struct line_current {
    vec2d p{};
    pscalar2d I{};
};

// the field of a line current: b(r) = mu0/(2 pi) (r - p)^-1 I
inline vec2d b_line_current(vec2d const& r, line_current const& c)
{
    using namespace ega;
    return (mu0 / (2.0 * pi)) * (inv(r - c.p) * c.I);
}

// its vector potential, A = -(mu0/(2 pi)) ln|r - p| I (zero on the unit circle about p)
inline pscalar2d A_line_current(vec2d const& r, line_current const& c)
{
    using namespace ega;
    return pscalar2d(-(mu0 / (2.0 * pi)) * std::log(nrm(r - c.p)) * value_t(c.I));
}

inline vec2d b_line_currents(vec2d const& r, std::span<line_current const> cs)
{
    vec2d b{};
    for (auto const& c : cs)
        b += b_line_current(r, c);
    return b;
}

// the force per unit length on a line current in a field: F = b I (= I e_z x B)
inline vec2d lorentz_force(vec2d const& b, pscalar2d I)
{
    using namespace ega;
    return b * I;
}

// the magnetic traction on a surface with unit normal n: t = (1/(2 mu0)) b n b
inline vec2d traction(vec2d const& b, vec2d const& n)
{
    using namespace ega;
    return (1.0 / (2.0 * mu0)) * (b * n * b); // a vector in ega2d: no other grade
}

// one quadrature point of a closed contour: the point, its OUTWARD unit normal, the
// length it stands for
struct contour_point {
    vec2d P{};
    vec2d n{};
    value_t ds{0.0};
};

// the midpoint rule on a circle -- spectrally accurate for a smooth periodic integrand
inline std::vector<contour_point> circle_contour(vec2d const& centre, value_t radius,
                                                 size_t N)
{
    std::vector<contour_point> c;
    c.reserve(N);
    value_t const ds = 2.0 * pi * radius / value_t(N);
    for (size_t i = 0; i < N; ++i) {
        value_t const phi = 2.0 * pi * (value_t(i) + 0.5) / value_t(N);
        vec2d const n{std::cos(phi), std::sin(phi)};
        c.push_back({centre + radius * n, n, ds});
    }
    return c;
}

// the wrench of a field on what a closed contour encloses: W = sum wdg(P, t(n)) ds
template <typename Field>
inline bivec2dp stress_wrench(Field&& b, std::span<contour_point const> contour)
{
    using namespace pga;
    bivec2dp W{};
    for (auto const& q : contour) {
        vec2d const t = traction(b(q.P), q.n);
        W += wdg(vec2dp(q.P.x, q.P.y, 1.0), vec2dp(t.x, t.y, 0.0)) * q.ds;
    }
    return W;
}

// the equivalent surface current of a magnetisation m at an edge with outward unit
// normal n: K = m ^ n [A/m] (the classical M x n, z component)
inline pscalar2d surface_current(vec2d const& m, vec2d const& n)
{
    using namespace ega;
    return wdg(m, n);
}

// a disc of radius R about c, uniformly magnetised m, as N boundary line currents
// K ds at the midpoints of N equal arcs
inline std::vector<line_current> magnetised_disc_currents(vec2d const& c, value_t R,
                                                          vec2d const& m, size_t N)
{
    std::vector<line_current> cs;
    cs.reserve(N);
    for (auto const& q : circle_contour(c, R, N))
        cs.push_back({q.P, surface_current(m, q.n) * q.ds});
    return cs;
}

/////////////////////////////////////////////////////////////////////////////////////////
// 3D
/////////////////////////////////////////////////////////////////////////////////////////

// a straight current I [A] from a to b
struct current_segment {
    vec3d a{};
    vec3d b{};
    value_t I{0.0};
};

// the field of a segment: (mu0 I / (4 pi)) (c_a - c_b) cmpl(u ^ rho^-1)
inline vec3d b_segment(vec3d const& r, current_segment const& s)
{
    using namespace ega;
    vec3d const L = s.b - s.a;
    vec3d const u = L / nrm(L);
    vec3d const Ra = r - s.a, Rb = r - s.b;
    vec3d const rho = Ra - value_t(dot(Ra, u)) * u;
    value_t const ca = value_t(dot(u, Ra)) / nrm(Ra);
    value_t const cb = value_t(dot(u, Rb)) / nrm(Rb);
    return (mu0 * s.I / (4.0 * pi) * (ca - cb)) * cmpl(wdg(u, inv(rho)));
}

// its vector potential, along the segment
inline vec3d A_segment(vec3d const& r, current_segment const& s)
{
    using namespace ega;
    vec3d const L = s.b - s.a;
    vec3d const u = L / nrm(L);
    vec3d const Ra = r - s.a, Rb = r - s.b;
    value_t const la = nrm(Ra) - value_t(dot(u, Ra));
    value_t const lb = nrm(Rb) - value_t(dot(u, Rb));
    return (mu0 * s.I / (4.0 * pi) * std::log(lb / la)) * u;
}

inline vec3d b_segments(vec3d const& r, std::span<current_segment const> ss)
{
    vec3d b{};
    for (auto const& s : ss)
        b += b_segment(r, s);
    return b;
}

namespace detail {

// two unit vectors completing the unit n to a right-handed frame (e_a, e_b, n)
inline std::pair<vec3d, vec3d> frame_about(vec3d const& n)
{
    using namespace ega;
    vec3d const h = std::abs(n.x) < 0.9 ? vec3d{1.0, 0.0, 0.0} : vec3d{0.0, 1.0, 0.0};
    vec3d ea = h - value_t(dot(h, n)) * n;
    ea = ea / nrm(ea);
    return {ea, cmpl(wdg(n, ea))};
}

} // namespace detail

// a circular loop about c, radius R, in the plane with unit normal n, as N chords
// carrying I counter-clockwise about n: its moment is I (N/2) R^2 sin(2 pi/N) n
inline std::vector<current_segment> circle_loop(vec3d const& c, vec3d const& n, value_t R,
                                                size_t N, value_t I)
{
    auto const [ea, eb] = detail::frame_about(n);
    std::vector<current_segment> s;
    s.reserve(N);
    auto P = [&](size_t k) {
        value_t const phi = 2.0 * pi * value_t(k) / value_t(N);
        return c + R * (std::cos(phi) * ea + std::sin(phi) * eb);
    };
    for (size_t k = 0; k < N; ++k)
        s.push_back({P(k), P(k + 1), I});
    return s;
}

// the force on a current element I dl in a field: F = cmpl(I dl ^ b) (= I dl x b)
inline vec3d lorentz_force(vec3d const& b, vec3d const& I_dl)
{
    using namespace ega;
    return cmpl(wdg(I_dl, b));
}

// the field of a point dipole mu at the origin: mu0 / (8 pi r^3) (mu + 3 r_hat mu r_hat)
inline vec3d b_dipole(vec3d const& r, vec3d const& mu)
{
    using namespace ega;
    value_t const rn = nrm(r);
    vec3d const rh = r / rn;
    value_t const c = mu0 / (8.0 * pi * rn * rn * rn);
    return c * (mu + 3.0 * gr1(rh * mu * rh));
}

// the same in the bivector convention, m = dual(mu), B = dual(b): the same expression,
// since I3 is central
inline bivec3d b_dipole(vec3d const& r, bivec3d const& m)
{
    using namespace ega;
    value_t const rn = nrm(r);
    vec3d const rh = r / rn;
    value_t const c = mu0 / (8.0 * pi * rn * rn * rn);
    return c * (m + 3.0 * gr2(rh * m * rh));
}

// the coupling of a dipole m_j at r (relative to m_i) to the field of m_i: <m_j B_i(r)>
// (= mu_j . b_i(r); the interaction energy is its negative)
inline value_t dipole_coupling(vec3d const& r, bivec3d const& m_i, bivec3d const& m_j)
{
    using namespace ega;
    return value_t(dot(m_j, b_dipole(r, m_i)));
}

// the magnetic traction on a surface with unit normal n: t = (1/(2 mu0)) b n b
inline vec3d traction(vec3d const& b, vec3d const& n)
{
    using namespace ega;
    return (1.0 / (2.0 * mu0)) * gr1(b * n * b); // the trivector part is b ^ n ^ b = 0
}

// one quadrature point of a closed surface: the point, its OUTWARD unit normal, the
// area it stands for
struct surface_point {
    vec3d P{};
    vec3d n{};
    value_t dA{0.0};
};

// a sphere about c: Gauss-Legendre in u = cos(theta) on `panels` panels of 12 points,
// the midpoint rule in phi on N_phi points (spectral for a smooth periodic integrand)
inline std::vector<surface_point> sphere_surface(vec3d const& c, value_t R, size_t panels,
                                                 size_t N_phi)
{
    auto const& g = detail::gauss_legendre12();
    std::vector<surface_point> s;
    s.reserve(panels * 12 * N_phi);
    value_t const h = 2.0 / value_t(panels), dphi = 2.0 * pi / value_t(N_phi);
    for (size_t p = 0; p < panels; ++p)
        for (size_t i = 0; i < 12; ++i) {
            value_t const u = -1.0 + h * (value_t(p) + 0.5 * (g.x[i] + 1.0));
            value_t const w = 0.5 * h * g.w[i], sn = std::sqrt(1.0 - u * u);
            for (size_t j = 0; j < N_phi; ++j) {
                value_t const phi = dphi * (value_t(j) + 0.5);
                vec3d const n{sn * std::cos(phi), sn * std::sin(phi), u};
                s.push_back({c + R * n, n, R * R * w * dphi});
            }
        }
    return s;
}

// the wrench of a field on what a closed surface encloses: W = sum wdg(P, t(n)) dA
template <typename Field>
inline bivec3dp stress_wrench(Field&& b, std::span<surface_point const> surface)
{
    using namespace pga;
    bivec3dp W{};
    for (auto const& q : surface) {
        vec3d const t = traction(b(q.P), q.n);
        W += wdg(vec3dp(q.P.x, q.P.y, q.P.z, 1.0), vec3dp(t.x, t.y, t.z, 0.0)) * q.dA;
    }
    return W;
}

// the equivalent surface current of a magnetisation m at a surface with outward unit
// normal n: K = m x n = cmpl(m ^ n) [A/m] (2D: K = m ^ n, its z component)
inline vec3d surface_current(vec3d const& m, vec3d const& n)
{
    using namespace ega;
    return cmpl(wdg(m, n));
}

// a sphere of radius R about c, uniformly magnetised m, as current rings: K = m x n
// circulates about m's axis with |K| = |m| sin(theta), so a band du = d cos(theta) of the
// surface carries |m| R du -- Gauss-Legendre in u on `panels` panels, each ring as
// N_seg chords
inline std::vector<current_segment> magnetised_sphere_currents(vec3d const& c, value_t R,
                                                               vec3d const& m,
                                                               size_t panels,
                                                               size_t N_seg)
{
    using namespace ega;
    auto const& g = detail::gauss_legendre12();
    value_t const M = nrm(m);
    vec3d const mh = m / M;
    std::vector<current_segment> s;
    value_t const h = 2.0 / value_t(panels);
    for (size_t p = 0; p < panels; ++p)
        for (size_t i = 0; i < 12; ++i) {
            value_t const u = -1.0 + h * (value_t(p) + 0.5 * (g.x[i] + 1.0));
            value_t const w = 0.5 * h * g.w[i];
            auto ring = circle_loop(c + R * u * mh, mh, R * std::sqrt(1.0 - u * u), N_seg,
                                    M * R * w);
            s.insert(s.end(), ring.begin(), ring.end());
        }
    return s;
}

} // namespace hd::ga

// Printing support for the types above -- at the very end, as in ga_usr_geodesics.hpp:
// ga_fmt_magnetostatics.hpp includes this header back, so hoisted to the top its
// specialisations would be written before the types exist.
#include "detail/fmt/ga_fmt_magnetostatics.hpp"
