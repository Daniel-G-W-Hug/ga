// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <cmath>
#include <iostream>
#include <memory>

#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

// hd::ga::det / lu_decomp / lu_backsubs come transitively via ga_pga.hpp.

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP dynamic_system3dp - Phase C: the Tao wheel-spindle vibration model
//   (Tao et al., Int. J. Mech. Sci. 232 (2022) 107620). A 6-joint serial stack on the
//   inertial root reproduces the paper's 5-DOF spindle: 3 prismatic (x,y,z) + 2 revolute
//   (theta,phi) vibration DOFs (massless frames) carrying a driven revolute SPIN joint
//   with the rotor body (mass m, principal inertias Jx=Jy, Jz, cm offset by e from the
//   spin axis). Aerostatic bearings = GROUNDED spatial springs on the NON-spinning
//   housing (the phi frame): two radial springs at axial +-Lb give the full radial
//   stiffness AND
//   -- via the lever arm -- the emergent tilt stiffness; one axial spring at the cm gives
//   the axial stiffness. All parameters are Table-1 values in strict SI. Validated
//   against the closed-form characteristic frequencies Eqs (15)-(19) and the centrifugal
//   / gyroscopic structure of Eq. (13).
//
// Lb is a CALIBRATION, not a physical bearing offset -- and the reason is now understood
// (derived in TODO/tao_eq13_derivation.md): Tao's lumped radial tilt stiffness
// (kx+ky)(l1^2+l2^2) is a factor of 4 too large vs a consistent two-bearing model
// (kx/2)(l1^2+l2^2). The 4x = 2 [Tao adds ky, but a tilt about y is resisted only by the
// x-springs] x 2 [kx is the TOTAL radial stiffness of the two parallel bearings (pinned
// by the f_x cross-check), so each bearing carries kx/2, not kx]. At the table's
// l1=l2=100 mm a consistent model gives f_theta ~34% low (~4430 vs ~6700 Hz). We absorb
// the over-count by inflating the lever to Lb = sqrt(K_tilt/kx) ~ 231 mm (~2.3x physical)
// with two kx/2 springs at +-Lb, so the emergent tilt stiffness kx Lb^2 = K_tilt
// reproduces Fig. 4. So Lb encodes Tao's 4x; the physical tilt frequency is really ~34%
// lower. Confirmed NOT caused by units (strict SI throughout) or by e (the +-e in l1,l2
// shifts l1^2+l2^2 by <1e-7). The GA model itself is correct -- it is Tao's coefficient
// that does not reduce from rigid-body mechanics.
//
// FALSIFIABLE PREDICTIONS (where the corrected model parts from the paper -- both errors
// are tilt-only, so x/y/z and all of section 3/4 are untouched; full note in
// TODO/tao_eq13_derivation.md): (a) Fig. 4(c,d) and Fig. 6(b) -- a SINGLE tilt peak, not
// the twin 6680/6727: the gyroscopic split is (Jz/Jt)*f_s ~ 1.76 Hz (unresolvable), and
// the center sits ~30-40% lower (~3.7-4.4 kHz) for a consistent tilt stiffness. (b) Tao's
// +-f_s/2 makes the peak SEPARATION grow with spin; his own measurement (sec 4.2.1: peaks
// "the same under 2250 and 2650 r/min") says it does NOT -- which supports the corrected
// (negligible-gyroscopic) picture and refutes Eq.13's gyroscopic origin. (c) Topography
// results / experimental wavelengths (Table 3) are UNCHANGED: that half feeds the surface
// model with the MEASURED f_b = 6253 Hz, never Eq.13's predicted frequency.
/////////////////////////////////////////////////////////////////////////////////////////

namespace {

// Tao Table-1 spindle performance parameters, in strict SI (m, kg, s, N). The table's
// stiffness rows read "N/um" (180, 950); a compliance of 180 um/N would be absurd and the
// frequency cross-checks only close as 180 N/um = 1.8e8 N/m. e is 2 um (not 2 mm).
struct tao_spindle_params {
    value_t m = 0.8;     // rotor mass [kg]
    value_t kx = 1.8e8;  // radial stiffness kx=ky [N/m]  (180 N/um)
    value_t kz = 9.5e8;  // axial stiffness   [N/m]       (950 N/um)
    value_t Jx = 5.4e-3; // transverse inertia Jx=Jy [kg m^2]  (5.4 g m^2)
    value_t Jz = 1.9e-4; // spin-axis inertia [kg m^2]         (0.19 g m^2)
    value_t e = 2.0e-6;  // mass eccentricity [m]              (2 um)
    value_t l1 = 0.1;    // bearing axial stations [m] (the +-e correction is negligible:
    value_t l2 = 0.1;    //   l1^2+l2^2 = 0.02 to 8 sig figs, so dropped here)
    value_t Rm = 0.05;   // rotor radius (axial-bearing lever) [m]
    value_t Rw = 0.15;   // grinding-wheel radius [m]
    value_t l3 = 0.18;   // wheel-surface-to-rotor-center axial distance [m]

    // Tao's lumped tilt stiffness (Eq.12/18/19) and the CALIBRATED radial-bearing offset
    // that makes the emergent tilt stiffness (2 * (kx/2) * Lb^2 = kx Lb^2) equal it. Lb
    // is a fit to Fig. 4, NOT a physical offset: at the table's l1=l2 a consistent
    // two-bearing model gives 1/4 of Tao's radial K_tilt -- the factor-4 over-count is
    // derived in
    // TODO/tao_eq13_derivation.md (Tao adds ky and uses the total kx per bearing).
    value_t K_tilt() const { return (kx + kx) * (l1 * l1 + l2 * l2) + kz * Rm * Rm; }
    value_t Lb() const { return std::sqrt(K_tilt() / kx); }

    // closed-form characteristic frequencies [Hz], Eqs (16)-(19) (no-spin tilt = f_th0)
    value_t f_x() const { return std::sqrt(kx / m) / (2.0 * pi); }
    value_t f_z() const { return std::sqrt(kz / m) / (2.0 * pi); }
    value_t f_th0() const { return std::sqrt(K_tilt() / Jx) / (2.0 * pi); }
};

// rotor body3dp with principal moments (Jx=Jy, Jz) about the cm, via a cuboid whose
// extents reproduce them: w=h with Izz=m(w^2+h^2)/12=Jz, Ixx=m(h^2+d^2)/12=Jx.
inline body3dp make_rotor_body(value_t m, value_t Jx, value_t Jz)
{
    value_t const w = std::sqrt(6.0 * Jz / m);
    value_t const d = std::sqrt(12.0 * Jx / m - 6.0 * Jz / m);
    return make_cuboid_body(m, w, w, d);
}

struct spindle_idx {
    size_t jx, jy, jz, jth, jph, spin, housing;
};

// Build the 6-joint Tao spindle into `sys`. Initial vibration offsets q0 = {x,y,z,th,ph};
// if spin_rate != 0 the spin joint is kinematically driven at that constant rate.
inline spindle_idx build_tao_spindle(dynamic_system3dp& sys, tao_spindle_params const& p,
                                     std::array<value_t, 5> const& q0 = {},
                                     value_t spin_rate = 0.0,
                                     std::array<value_t, 5> const& qdot0 = {},
                                     std::array<value_t, 2> const& tilt_rest = {},
                                     value_t c_bearing = 0.0)
{
    auto const massless = body3dp{};
    auto const rotor = make_rotor_body(p.m, p.Jx, p.Jz);
    sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
    sys.add_frame(static_frame3dp("W"));
    sys.add_prismatic_body(static_frame3dp("Jx"), massless, vec3dp{1, 0, 0, 0}, q0[0],
                           qdot0[0]);
    sys.add_prismatic_body(static_frame3dp("Jy"), massless, vec3dp{0, 1, 0, 0}, q0[1],
                           qdot0[1]);
    sys.add_prismatic_body(static_frame3dp("Jz"), massless, vec3dp{0, 0, 1, 0}, q0[2],
                           qdot0[2]);
    sys.add_revolute_body(static_frame3dp("Jth"), massless, O_3dp, vec3dp{1, 0, 0, 0},
                          q0[3], qdot0[3]);
    sys.add_revolute_body(static_frame3dp("Jph"), massless, O_3dp, vec3dp{0, 1, 0, 0},
                          q0[4], qdot0[4]);
    // Spin joint: axis e3 through (e,0,0) so the rotor cm (body origin) orbits at radius
    // e. If spin_rate != 0 the spin is kinematically DRIVEN (motor-clamped, Omega-ddot =
    // 0, matching Tao) -- see set_driven_rate below; assemble_mass_bias now retains a
    // driven joint's inertia as a moving base, so the clamped gyroscopic dynamics are
    // produced.
    sys.add_revolute_body(static_frame3dp("rotor"), rotor, vec3dp{p.e, 0, 0, 1},
                          vec3dp{0, 0, 1, 0}, 0.0, spin_rate);

    spindle_idx ix;
    ix.jx = sys.index_of("Jx");
    ix.jy = sys.index_of("Jy");
    ix.jz = sys.index_of("Jz");
    ix.jth = sys.index_of("Jth");
    ix.jph = sys.index_of("Jph");
    ix.spin = sys.index_of("rotor");
    ix.housing = ix.jph; // bearings on the non-spinning phi frame

    value_t const Lb = p.Lb(), kr = p.kx / 2.0;
    // Optional STATIC bearing-anchor misalignment (Zhou-flatness bridge, Phase Dyn.1):
    // offset the two radial ground anchors so the unstressed (zero-force) configuration
    // is tilted by (alpha about e1 -> theta, beta about e2 -> phi). Under that small tilt
    // the
    // +-Lb body points tip (linearized) to (+-Lb*beta, -+Lb*alpha, +-Lb); placing each
    // ground anchor there makes the tilted pose the spring rest, so the spindle SETTLES
    // to (theta,phi) = (alpha,beta) -- the static alignment deviation, now produced by
    // the dynamics rather than hand-set. Default {0,0} keeps the anchors at +-Lb (Phase C
    // byte-unchanged). c_bearing adds isotropic point damping so the settling is
    // non-oscillatory (default 0 -> Phase C byte-unchanged).
    value_t const ta = tilt_rest[0], tb = tilt_rest[1];
    sys.add_grounded_spring(ix.housing, vec3dp{0, 0, Lb, 1},
                            vec3dp{Lb * tb, -Lb * ta, Lb, 1}, vec3dp{kr, kr, 0, 0},
                            c_bearing);
    sys.add_grounded_spring(ix.housing, vec3dp{0, 0, -Lb, 1},
                            vec3dp{-Lb * tb, Lb * ta, -Lb, 1}, vec3dp{kr, kr, 0, 0},
                            c_bearing);
    sys.add_grounded_spring(ix.housing, O_3dp, O_3dp, vec3dp{0, 0, p.kz, 0}, c_bearing);
    if (spin_rate != 0.0) sys.set_driven_rate(ix.spin, spin_rate); // motor-clamped spin
    return ix;
}

// Measure the oscillation frequency [Hz] of joint `jidx` by downward zero-crossings (with
// linear sub-step interpolation) of its coordinate over `nsteps` of size `dt`.
inline value_t measure_freq(dynamic_system3dp& sys, size_t jidx, value_t dt, int nsteps)
{
    std::vector<value_t> cross;
    value_t qprev = sys.joint_phi(jidx), t = 0.0;
    for (int n = 0; n < nsteps; ++n) {
        sys.step(dt);
        t += dt;
        value_t const q = sys.joint_phi(jidx);
        if (qprev > 0.0 && q <= 0.0) // downward crossing, interpolated
            cross.push_back(t - dt * q / (q - qprev));
        qprev = q;
    }
    value_t T = 0.0;
    for (size_t i = 1; i < cross.size(); ++i)
        T += cross[i] - cross[i - 1];
    return (cross.size() - 1.0) / T; // mean period -> frequency
}

} // namespace

TEST_SUITE("PGA3DP: Tao wheel-spindle (Phase C)")
{

    TEST_CASE("pga3dp: Tao spindle - radial/axial/tilt free frequencies (Phase C.2)")
    {
        fmt::println("pga3dp: Tao spindle - characteristic frequencies (Phase C.2)");

        // No-spin free-vibration frequencies of the assembled spindle. Each vibration DOF
        // is released from a small offset (all others at rest, no spin, no damping) and
        // its frequency measured. f_x, f_z are FIRST-PRINCIPLES (no calibration): the
        // radial/axial bearing stiffnesses are the table values, so matching Eqs
        // (16),(17) genuinely validates the assembled model. f_th0 instead checks that
        // the CALIBRATED Lb reproduces Tao's tilt stiffness -- matched by construction,
        // not derived (the ~factor-4 gap to a consistent l1=l2 two-spring model is
        // unresolved; see the suite header). With spin off the DOFs decouple (gyroscopic
        // coupling vanishes).
        tao_spindle_params const p;
        fmt::println("  Lb (calibrated radial-bearing offset) = {:.4f} m (l1 = {:.3f})",
                     p.Lb(), p.l1);

        // radial f_x: release x
        {
            dynamic_system3dp sys;
            build_tao_spindle(sys, p, {1.0e-7, 0, 0, 0, 0});
            value_t const f = measure_freq(sys, sys.index_of("Jx"), 2.0e-8, 60000);
            fmt::println("  f_x = {:.1f} Hz  (Eq.16: {:.1f}, paper 2385)", f, p.f_x());
            CHECK(f == doctest::Approx(p.f_x()).epsilon(2e-3));
        }
        // axial f_z: release z
        {
            dynamic_system3dp sys;
            build_tao_spindle(sys, p, {0, 0, 1.0e-7, 0, 0});
            value_t const f = measure_freq(sys, sys.index_of("Jz"), 1.0e-8, 60000);
            fmt::println("  f_z = {:.1f} Hz  (Eq.17: {:.1f}, paper 5490)", f, p.f_z());
            CHECK(f == doctest::Approx(p.f_z()).epsilon(2e-3));
        }
        // no-spin tilt f_th0: release theta (validates the emergent calibrated tilt
        // stiff.)
        {
            dynamic_system3dp sys;
            build_tao_spindle(sys, p, {0, 0, 0, 1.0e-7, 0});
            value_t const f = measure_freq(sys, sys.index_of("Jth"), 1.0e-8, 50000);
            fmt::println("  f_th0 = {:.1f} Hz  (Eq.18/19 at f_s=0: {:.1f})", f,
                         p.f_th0());
            CHECK(f == doctest::Approx(p.f_th0()).epsilon(3e-3));
        }

        fmt::println("");
    }

    TEST_CASE("pga3dp: Tao spindle - gyroscopic + centrifugal emergence (Phase C.3)")
    {
        fmt::println("pga3dp: Tao spindle - gyroscopic + centrifugal (Phase C.3)");

        // The two "for free" payoffs of building on rigid-body GA dynamics: with the spin
        // joint DRIVEN at omega, BOTH the gyroscopic tilt coupling (J_y-J_z)omega phi-dot
        // (Eq.13) and the centrifugal unbalance forcing m e omega^2 EMERGE from the
        // Newton-Euler velocity-product bias of the spinning offset-cm rotor -- nothing
        // is added by hand. Validated instantaneously (forward dynamics at a chosen
        // state) against the closed-form coefficients, the same altitude as the A.3
        // centrifugal slider.
        tao_spindle_params const p;
        value_t const f_s = 50.0;             // spin frequency [Hz] (N_w = 3000 r/min)
        value_t const omega = 2.0 * pi * f_s; // spin rate [rad/s]

        // --- Gyroscopic coupling == textbook polar-inertia value (clamped spin) -------
        // With the spin motor-CLAMPED (driven, Omega-ddot = 0, as Tao assumes), a tilt
        // rate theta-dot induces phi-dd through the gyroscopic velocity-product bias. The
        // standard rotordynamics result is that the coupling moment is the POLAR angular
        // momentum J_z*omega, so  phi-dd = (J_z / J_x) * omega * theta-dot  (forward /
        // backward whirl split = (J_z/J_x)*omega). The GA rigid-body model reproduces
        // this textbook value exactly. (assemble_mass_bias now retains the driven rotor's
        // inertia, so the clamped case runs; the value equals the free-spin case here
        // because the spin-DOF coupling to phi vanishes at this config.)
        //
        // DEVIATION (now resolved -- derived in TODO/tao_eq13_derivation.md): Tao's Eq.13
        // uses (J_y - J_z)*omega ~ J_x*omega (the TRANSVERSE inertia), giving his +-f_s/2
        // ~ 25 Hz split. That disagrees with both standard rotordynamics and this model,
        // which give the polar-inertia whirl split (J_z/J_x)*omega ~ 1.76 Hz. It is a
        // PAPER ERROR, not an artefact of his Eq.1/2 angles: projecting the body-frame
        // Euler equations onto theta,phi, the time derivative d/dt(J*omega_c) contributes
        // a kinematic frame-rotation term +J_t*omega*phi-dot that combines with the
        // cross-product term (J_z - J_t)*omega*phi-dot to leave the correct coupling
        // J_z*omega. Tao dropped the +J_t*omega*phi-dot term, leaving (J_t - J_z) = (J_y
        // - J_z). The GA model keeps the full derivative, so it is correct.
        {
            value_t const thd0 = 1.0; // theta-dot [rad/s]
            dynamic_system3dp sys;
            auto const ix = build_tao_spindle(sys, p, {}, omega, {0, 0, 0, thd0, 0});
            value_t const phidd = sys.joint_accel(ix.jph);
            value_t const coupling =
                p.Jz / p.Jx * omega;                     // textbook (J_z/J_t)*omega [1/s]
            value_t const split = coupling / (2.0 * pi); // whirl split [Hz]
            fmt::println(
                "  gyroscopic: phi-dd/theta-dot = {:.4f} (textbook (Jz/Jx)*omega "
                "= {:.4f}); whirl split = {:.2f} Hz  [Tao Eq.13 -> ~{:.0f} Hz; "
                "deviation flagged]",
                phidd / thd0, coupling, split, f_s);
            CHECK(std::abs(phidd / thd0) == doctest::Approx(coupling).epsilon(1e-3));
        }

        // --- Centrifugal unbalance forcing: x-dd at rest, spin on ---------------------
        // Eq.13 x-equation: m x-dd + ... = m e omega^2 sin(omega t) + F_x. At t=0 (cm
        // offset along the spin-axis radial, x=0, x-dot=0, no grinding force) the radial
        // acceleration the lib returns must be the centrifugal forcing / mass = e
        // omega^2.
        {
            dynamic_system3dp sys;
            auto const ix = build_tao_spindle(sys, p, {}, omega); // only spin, all q=0
            value_t const xdd = sys.joint_accel(ix.jx);
            value_t const ydd = sys.joint_accel(ix.jy);
            value_t const pred = p.e * omega * omega; // m e omega^2 / m
            fmt::println("  centrifugal: x-dd = {:.5f} m/s^2 (predicted e*omega^2 = "
                         "{:.5f}); y-dd = {:.2e}",
                         xdd, pred, ydd);
            CHECK(std::abs(xdd) == doctest::Approx(pred).epsilon(1e-3));
            CHECK(std::abs(ydd) < 1e-3 * pred); // forcing is along x at this spin phase
        }

        fmt::println("");
    }

    TEST_CASE("pga3dp: Tao spindle - axial runout z_b = z - R_w*phi (Phase C.4)")
    {
        fmt::println("pga3dp: Tao spindle - axial runout z_b (Eq.14, Phase C.4)");

        // Eq.14: the total wheel-wafer axial displacement at the grinding rim is
        //   z_b = z - R_w * phi
        // -- the axial DOF z plus the rim's axial swing under tilt phi (a point at radius
        // R_w tips out of the grinding plane by R_w*phi). z_b is THE surface-error driver
        // in Tao's model (it sets the grinding-mark depth, Eq.26). We put the spindle in
        // a known (z, phi) configuration and confirm the world axial position of a rim
        // point at radius R_w, computed by the GA model, equals z - R_w*phi.
        tao_spindle_params const p;
        value_t const z0 = 1.0e-6, phi0 = 5.0e-6; // small axial offset + tilt [m, rad]
        dynamic_system3dp sys;
        build_tao_spindle(sys, p, {0, 0, z0, 0, phi0}); // z on Jz, phi on Jph
        size_t const jz = sys.index_of("Jz"), jph = sys.index_of("Jph");

        // world axial position of a rim point at radius R_w on the (non-spinning) housing
        vec3dp const rim_w =
            unitize(move3dp(vec3dp{p.Rw, 0.0, 0.0, 1.0}, sys.get_pos_trafo("Jph", "W")));
        value_t const zb_model = rim_w.z; // from the GA pose
        value_t const zb_eq14 = sys.joint_phi(jz) - p.Rw * std::sin(sys.joint_phi(jph));
        value_t const zb_lin = z0 - p.Rw * phi0; // Eq.14 (small angle)
        fmt::println("  z = {:.3e}, phi = {:.3e} -> z_b(model) = {:.4e}, z_b(Eq.14) = "
                     "{:.4e} (linear {:.4e})",
                     z0, phi0, zb_model, zb_eq14, zb_lin);
        CHECK(zb_model ==
              doctest::Approx(zb_eq14).epsilon(1e-9)); // GA pose == z - R_w sin(phi)
        CHECK(zb_model ==
              doctest::Approx(zb_lin).epsilon(1e-3)); // == Eq.14 to small angle

        fmt::println("");
    }

    // -------------------------------------------------------------------------------------
    // Phase Dyn.1 -- the static-deviation run: the Zhou-geometry <-> Tao-dynamics BRIDGE.
    //
    // Until now the two error sources lived in separate worlds: the global FLATNESS came
    // from a HAND-SET kinematic wafer tilt (Zhou Figs 5-7, ga_appl3dp_appl_test.hpp), and
    // the spindle DYNAMICS (Phase C) only produced free vibrations about an ideal pose.
    // This case fuses them. The alignment misalignment (alpha,beta) is injected as a
    // STATIC bearing-anchor offset (build_tao_spindle's tilt_rest), so under its own
    // springs/dampers the spindle SETTLES to a tilted equilibrium
    // (theta,phi)=(alpha,beta). That settled tilt -- now an OUTPUT of the dynamics, not
    // an input -- is then fed into the Zhou single-grain cutting-path sampler and shown
    // to reproduce the same cone/dome/bowl global-flatness shapes.
    //
    // Axis relabeling: the spindle's two independent tilt DOFs are theta (about e1) and
    // phi (about e2); the Zhou sampler tilts the wafer by (a about e2, b about e3). Each
    // spindle tilt DOF maps onto one of Zhou's two perpendicular tilt slots -- theta -> a
    // (the SYMMETRIC slot, along the grain's radial line -> CONE), phi -> b (the
    // perpendicular slot -> sign-dependent DOME/BOWL). The bridge is the
    // magnitude-and-shape correspondence: a dynamically-settled tilt of magnitude alpha
    // carves Zhou's cone, one of magnitude beta carves the dome.
    //
    // A constant infeed (prescribed axial translation on Jz) is a UNIFORM z offset of the
    // whole path -- it sets the engagement depth but not the flatness SHAPE -- so it is
    // omitted from this shape gate and enters in the Step-2 dynamic-vibration run (where
    // the time-varying z_b adds the waviness on top of this mean tilt).
    TEST_CASE(
        "pga3dp: Tao spindle - static tilt deviation -> Zhou flatness (Phase Dyn.1)")
    {
        fmt::println("pga3dp: Tao spindle - static tilt deviation -> Zhou flatness "
                     "(Phase Dyn.1)");

        tao_spindle_params const p;
        value_t const deg = pi / 180.0;
        value_t const alpha = 0.1 * deg, beta = 0.1 * deg; // target misalignment [rad]
        value_t const c_bearing = 4000.0; // ~critical tilt damping -> monotone settle
        value_t const dt = 2.0e-8;
        int const nsteps = 60000; // 1.2 ms ~ 8 tilt periods: fully settled

        // Release the spindle from the IDEAL (theta=phi=0) pose with the anchors offset
        // to the misaligned rest; integrate until it settles, and report the equilibrium
        // tilt.
        auto settle = [&](value_t a, value_t b) {
            dynamic_system3dp sys;
            auto const ix = build_tao_spindle(sys, p, {}, 0.0, {}, {a, b}, c_bearing);
            for (int n = 0; n < nsteps; ++n)
                sys.step(dt);
            struct {
                value_t th, ph, rad;
            } r;
            r.th = sys.joint_phi(ix.jth);
            r.ph = sys.joint_phi(ix.jph);
            // the radial/axial DOFs are NOT driven by a pure tilt offset (the two anchor
            // forces cancel) -- confirm they stay at the ideal
            r.rad =
                std::max({std::abs(sys.joint_phi(ix.jx)), std::abs(sys.joint_phi(ix.jy)),
                          std::abs(sys.joint_phi(ix.jz))});
            return r;
        };

        // (theta,phi) = (alpha,0): settles to the alpha-tilt equilibrium
        auto const ce = settle(alpha, 0.0);
        fmt::println("  settle(alpha,0): theta_eq = {:.5e} (target {:.5e}), phi_eq = "
                     "{:.2e}, |x,y,z| <= {:.1e} m",
                     ce.th, alpha, ce.ph, ce.rad);
        CHECK(ce.th == doctest::Approx(alpha).epsilon(1e-3)); // dynamics reproduce alpha
        CHECK(std::abs(ce.ph) < 1e-3 * alpha);                // no phi cross-coupling
        CHECK(ce.rad < 1e-7);                                 // radial/axial undisturbed

        // (theta,phi) = (0,beta): settles to the beta-tilt equilibrium
        auto const de = settle(0.0, beta);
        fmt::println("  settle(0,beta): phi_eq   = {:.5e} (target {:.5e}), theta_eq = "
                     "{:.2e}",
                     de.ph, beta, de.th);
        CHECK(de.ph == doctest::Approx(beta).epsilon(1e-3)); // dynamics reproduce beta
        CHECK(std::abs(de.th) < 1e-3 * beta);

        // Confirm (alpha,beta) is a genuine EQUILIBRIUM: with the spindle placed exactly
        // at the misaligned rest (no spin, no velocity), the tilt accelerations vanish --
        // so the settled state above is the static fixed point, not a slow transient.
        // Compare to the restoring acceleration from the ideal (theta=phi=0) pose.
        {
            dynamic_system3dp eq, off;
            auto const ie =
                build_tao_spindle(eq, p, {0, 0, 0, alpha, beta}, 0.0, {}, {alpha, beta});
            auto const io = build_tao_spindle(off, p, {}, 0.0, {}, {alpha, beta});
            value_t const a_eq = std::abs(eq.joint_accel(ie.jth));
            value_t const a_off =
                std::abs(off.joint_accel(io.jth)); // restoring from ideal
            fmt::println(
                "  equilibrium: |theta-ddot| at rest = {:.3e} rad/s^2 (restoring "
                "from ideal = {:.3e}); ratio {:.1e}",
                a_eq, a_off, a_eq / a_off);
            CHECK(a_eq <
                  1e-5 * a_off); // negligible vs the restoring scale -> equilibrium
        }

        // ---- Flatness consequence: feed the settled tilt into the Zhou sampler
        // ---------- One wheel grain (radius r1=R, axis offset L=R, half-overlap) sampled
        // in the tilted, spinning wafer frame over one wafer revolution; z = wafer-axis
        // component = profile height. Identical to the committed Zhou Figs 5-7 case, but
        // the tilt now comes from the dynamics (theta_eq -> a, phi_eq -> b). Geometry in
        // mm; only the dimensionless settled ANGLES cross over from the SI dynamics.
        double const R = 150.0, L = R, r1 = R;
        double const w1 = rpm2radps(1500.0), w2 = rpm2radps(50.0); // n2/n1 = 1/30
        auto profile = [&](double a, double b) {
            kinematic_system3dp sys;
            sys.add_frame(static_frame3dp("W"));
            sys.add_frame(static_frame3dp("wheel_rot", vec3dp{0.0, L, 0.0, 1.0}),
                          kin_state3dp{.omega = vec3dp{w1, 0.0, 0.0, 0.0}},
                          sys.index_of("W"));
            sys.add_frame(static_frame3dp("grain", vec3dp{0.0, -r1, 0.0, 1.0}),
                          kin_state3dp{}, sys.index_of("wheel_rot"));
            sys.add_frame(static_frame3dp("wafer_tilt", O_3dp, vec3dp{0.0, a, b, 0.0}),
                          kin_state3dp{}, sys.index_of("W"));
            sys.add_frame(static_frame3dp("wafer_rot"),
                          kin_state3dp{.omega = vec3dp{w2, 0.0, 0.0, 0.0}},
                          sys.index_of("wafer_tilt"));
            size_t const wheel = sys.index_of("wheel_rot"),
                         wafer = sys.index_of("wafer_rot");
            double const T = 2.0 * pi / w2;
            int const N = 6000;
            struct {
                double zmin = 1e9, zmax = -1e9, z_ctr = 0.0, rmax = 0.0;
            } pr;
            for (int k = 0; k <= N; ++k) {
                double const t = T * double(k) / double(N);
                sys.set_pose(wheel, pose3dp{vec3dp{0.0, L, 0.0, 1.0},
                                            vec3dp{w1 * t, 0.0, 0.0, 0.0}});
                sys.set_pose(wafer, pose3dp{O_3dp, vec3dp{w2 * t, 0.0, 0.0, 0.0}});
                auto const g =
                    unitize(move3dp(O_3dp, sys.get_pos_trafo("grain", "wafer_rot")));
                double const z = g.x, r = std::sqrt(g.y * g.y + g.z * g.z);
                if (r > R + 1e-6) continue;
                pr.zmin = std::min(pr.zmin, z);
                pr.zmax = std::max(pr.zmax, z);
                if (r < 3.0) pr.z_ctr = z;
                pr.rmax = std::max(pr.rmax, r);
            }
            return pr;
        };

        double const amp = R * std::tan(0.1 * deg); // 0.2618 mm unclipped rim amplitude

        auto const flat = profile(0.0, 0.0); // ideal alignment -> FLAT
        CHECK(std::abs(flat.zmin) < 1e-9);
        CHECK(std::abs(flat.zmax) < 1e-9);

        auto const cone = profile(ce.th, 0.0); // settled alpha -> CONE (symmetric)
        CHECK(std::abs(cone.zmax + cone.zmin) < 1e-3 * amp);
        CHECK(std::abs(cone.z_ctr) < 1e-3);
        CHECK(cone.zmax > 0.5 * amp);
        CHECK(cone.zmax < 1.01 * amp);
        CHECK(cone.rmax > 0.99 * R);

        auto const dome = profile(0.0, de.ph); // settled +beta -> CONVEX dome
        CHECK(dome.zmin == doctest::Approx(0.0).epsilon(1e-4));
        CHECK(dome.zmax > 0.3 * amp);
        CHECK(dome.zmax < 0.75 * cone.zmax); // ~half the cone (Zhou: beta ~ half alpha)

        auto const bowl = profile(0.0, -de.ph); // settled -beta -> CONCAVE bowl (mirror)
        CHECK(bowl.zmax == doctest::Approx(0.0).epsilon(1e-4));
        CHECK(bowl.zmin == doctest::Approx(-dome.zmax).epsilon(1e-4));

        fmt::println("  flat z in [{:+.4f},{:+.4f}] | cone [{:+.4f},{:+.4f}] | dome "
                     "[{:+.4f},{:+.4f}] | bowl [{:+.4f},{:+.4f}] mm",
                     flat.zmin, flat.zmax, cone.zmin, cone.zmax, dome.zmin, dome.zmax,
                     bowl.zmin, bowl.zmax);
        fmt::println("  -> dynamically-settled tilt reproduces Zhou's cone/dome/bowl "
                     "flatness (dome/cone = {:.2f})",
                     dome.zmax / cone.zmax);
        fmt::println("");
    }

    // -------------------------------------------------------------------------------------
    // Phase Dyn.2 -- the dynamic-vibration run: ONE live spindle run yields BOTH error
    // sources at once. On top of the Phase-Dyn.1 static misalignment (the MEAN tilt ->
    // global flatness), the spindle now spins (driven Omega, mass unbalance e) and its
    // lightly-damped tilt DOFs RING at their natural frequency, so the axial runout
    //
    //     z_b = z - R_w * sin(phi)      (Eq.14)
    //
    // becomes time-varying. The live z_b then carries both signatures of Tao's model: its
    // MEAN is the static tilt (flatness, Zhou cone/dome of Dyn.1), and its OSCILLATION at
    // the tilt natural frequency f_b (Tao Eq.22 = f_th0) carves the WAVINESS whose
    // spatial wavelength is (surface speed)/f_b (Phase D.1). The headline: f_b is now an
    // OUTPUT of the live dynamics, not the measured 6253 Hz fed into the Phase-D.1
    // surface model.
    //
    // The constant grinding infeed (feed ~ 30 um/min) is a UNIFORM slow z-ramp -- over
    // the ~1.6 ms vibration window it advances < 1 nm, i.e. a constant offset of the mean
    // cutting depth, invisible to the ripple. It sets the engagement depth (and combines
    // with z_b for the carve), not the waviness frequency, so the axial DOF is left free
    // here (z ~ 0) and the feed enters the depth in the Step-3 visualization.
    TEST_CASE(
        "pga3dp: Tao spindle - dynamic vibration: flatness + waviness (Phase Dyn.2)")
    {
        fmt::println("pga3dp: Tao spindle - dynamic vibration: flatness + waviness "
                     "(Phase Dyn.2)");

        tao_spindle_params const p;
        value_t const deg = pi / 180.0;
        value_t const beta = 0.1 * deg; // static mean tilt -> global flatness
        value_t const delta = 1.0e-5;   // tilt perturbation -> z_b ripple [rad]
        value_t const f_s = 50.0;       // 3000 r/min driven spin
        value_t const Omega = 2.0 * pi * f_s;
        value_t const dt = 2.0e-8;
        int const nsteps = 80000; // 1.6 ms ~ 11 tilt periods -> ~10 z_b crossings

        // Static misalignment via the anchor offset (mean flatness); UNDAMPED so the tilt
        // rings cleanly at its natural freq over the window (the table damping ratio
        // ~5e-5 is negligible here); driven spin + unbalance on; phi released `delta`
        // above its tilted equilibrium so z_b ripples about the mean tilt.
        dynamic_system3dp sys;
        auto const ix = build_tao_spindle(sys, p, {0, 0, 0, 0.0, beta + delta}, Omega, {},
                                          {0.0, beta}, 0.0);

        std::vector<value_t> zb, ts;
        zb.reserve(nsteps + 1);
        ts.reserve(nsteps + 1);
        value_t phi_sum = 0.0, t = 0.0;
        for (int n = 0; n <= nsteps; ++n) {
            value_t const z = sys.joint_phi(ix.jz);
            value_t const phi = sys.joint_phi(ix.jph);
            zb.push_back(z - p.Rw * std::sin(phi));
            ts.push_back(t);
            phi_sum += phi;
            if (n < nsteps) {
                sys.step(dt);
                t += dt;
            }
        }
        value_t const phi_mean = phi_sum / value_t(nsteps + 1);

        // (1) MEAN tilt == the static beta -> the GLOBAL FLATNESS survives the vibration:
        // averaging the live phi recovers the prescribed alignment tilt (Dyn.1
        // cone/dome).
        fmt::println("  mean tilt phi = {:.5e} rad (static beta = {:.5e}) -> flatness",
                     phi_mean, beta);
        CHECK(phi_mean == doctest::Approx(beta).epsilon(2e-2));

        // (2) z_b OSCILLATION frequency == the tilt natural freq f_b (Tao Eq.22 = f_th0),
        // measured from the live run by downward zero-crossings of (z_b - mean), linearly
        // interpolated.
        value_t zb_mean = 0.0;
        for (auto v : zb)
            zb_mean += v;
        zb_mean /= value_t(zb.size());
        std::vector<value_t> cross;
        for (size_t i = 1; i < zb.size(); ++i) {
            value_t const a = zb[i - 1] - zb_mean, b = zb[i] - zb_mean;
            if (a > 0.0 && b <= 0.0) cross.push_back(ts[i - 1] + dt * a / (a - b));
        }
        value_t Tsum = 0.0;
        for (size_t i = 1; i < cross.size(); ++i)
            Tsum += cross[i] - cross[i - 1];
        value_t const f_b = (cross.size() - 1.0) / Tsum;
        value_t const A_b = 0.5 * (*std::max_element(zb.begin(), zb.end()) -
                                   *std::min_element(zb.begin(), zb.end())); // [m]
        fmt::println(
            "  z_b ripple: f_b = {:.1f} Hz (Eq.22 f_th0 = {:.1f}), A_b = {:.3f} um", f_b,
            p.f_th0(), A_b * 1e6);
        CHECK(f_b == doctest::Approx(p.f_th0()).epsilon(5e-3));

        // (3) WAVINESS wavelength lambda = (surface speed)/f_b (Phase D.1), now with f_b
        // from the LIVE dynamics. WMD lambda_m uses the wheel surface speed v_w =
        // n_w*R_w; WCD lambda_c uses the wafer surface speed v_s = n_s*r at a
        // representative radius.
        value_t const v_w = rpm2radps(3000.0) * p.Rw; // wheel surface speed [m/s]
        value_t const v_s = rpm2radps(265.0) * 0.030; // wafer @ r=30 mm [m/s]
        value_t const lambda_m = v_w / f_b * 1e3;     // WMD wavelength [mm]
        value_t const lambda_c = v_s / f_b * 1e3;     // WCD wavelength [mm]
        fmt::println("  waviness from live f_b: lambda_m = {:.3f} mm (WMD), lambda_c = "
                     "{:.1f} um (WCD)",
                     lambda_m, lambda_c * 1e3);
        CHECK(lambda_m > 5.0); // WMD ~ mm scale (Tao Table 3: 5-7 mm)
        CHECK(lambda_m < 10.0);
        CHECK(lambda_c * 1e3 > 50.0); // WCD ~ 100 um scale (Tao Table 3: 127-163 um)
        CHECK(lambda_c * 1e3 < 300.0);

        fmt::println(
            "  -> one live run: MEAN tilt = flatness (Zhou cone/dome), z_b ripple "
            "at f_b = waviness (lambda = v/f_b)");
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: Tao wheel-spindle (Phase C)")

// =================================================================================
// PGA3DP dynamic_system3dp - Phase D.2b: closing the grinding force loop
//
// D.2a (ga_appl3dp_appl_test.hpp) computed the grinding force FEED-FORWARD from a fixed
// engagement. D.2b promotes the contact element into dynamic_system3dp (set_contact_force
// + the swappable grinding_law3dp + the assemble_mass_bias fold-in) and CLOSES the loop:
// the spindle's own axial compliance sets the engagement, the force law turns it into a
// wrench, the wrench feeds back. The closed loop has an analytic equilibrium, which is
// the gate.
//
// SIMPLIFIED RIG (the full Tao 5-DOF spindle of Phase C is the next refinement): a single
// axial DOF -- a prismatic joint along the wafer normal e1 carrying the rotor mass m,
// with an axial bearing spring k_spring + damper c to ground -- plus a DRIVEN wheel spin
// about e1 (massless, purely kinematic) so the contact has a real relative-sliding speed
// n_w*R_w. A static wafer frame supplies the wafer side of the relative velocity. The
// in-plane (tangential) part of the force is perpendicular to the axial DOF, so only the
// axial F_z drives the loop -- the equation of motion is
//
//     m z'' + c z' + k_spring z = F_z = k_grind (delta0 - z)        [z_ref = 0 at start]
//
// i.e. a damped oscillator with effective stiffness (k_spring + k_grind) and an offset
// drive k_grind*delta0, settling to
//
//     z_eq = k_grind delta0 / (k_spring + k_grind),   F_eq = k_grind (delta0 - z_eq).
//
// The HEADLINE: closing the loop pulls the force BELOW the feed-forward nominal
// (k_grind*delta0 = 100 N) because the spindle deflects and backs out of the cut; in the
// weak-coupling / stiff-spindle limit (k_grind << k_spring) F_eq -> the feed-forward
// value, recovering D.2a. Reaction on the wafer (chuck side) stays OFF here -- topology 1
// (D.2c).
// =================================================================================
TEST_SUITE("PGA3DP: grinding force loop (Phase D.2b)")
{

    TEST_CASE("pga3dp: closed grinding-force loop, axial spindle (Phase D.2b)")
    {
        fmt::println("");
        fmt::println("pga3dp: closed grinding-force loop, axial spindle (Phase D.2b)");
        fmt::println("");

        // SI units. Calibration from D.2a: delta0 = 1 um, k_grind = 1e8 N/m -> 100 N
        // nominal; mu = 0.3536 (Tao in-plane). Axial bearing k_spring = Tao k_z ~ 950
        // N/um.
        value_t const m = 0.8;                 // rotor mass [kg] (Tao Table 1)
        value_t const R_w = 0.15;              // grinding-wheel radius [m]
        value_t const n_w = rpm2radps(3000.0); // wheel spin [rad/s]
        value_t const delta0 = 1.0e-6;         // nominal engagement [m]
        value_t const k_grind = 1.0e8; // depth-law stiffness [N/m] (-> 100 N at delta0)
        value_t const mu = std::sqrt(2.0) * 25.0 / 100.0; // in-plane / axial = 0.3536
        value_t const k_spring = 9.5e8; // axial bearing stiffness [N/m] (Tao k_z)
        value_t const k_tot = k_spring + k_grind;
        value_t const c = 2.0 * 0.7 * std::sqrt(k_tot * m); // damping, zeta_eff = 0.7

        // analytic closed-loop equilibrium
        value_t const z_eq = k_grind * delta0 / k_tot;
        value_t const F_eq = k_grind * (delta0 - z_eq);
        value_t const F_ff = k_grind * delta0; // feed-forward nominal (D.2a) = 100 N

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0}); // isolate the axial loop
        sys.add_frame(static_frame3dp("W"));         // inertial root
        // static wafer surface (wafer side of the relative-sliding velocity)
        sys.add_frame(static_frame3dp("wafer"), sys.index_of("W"));
        // axial DOF: prismatic along e1, carries the rotor mass, bearing spring + damper
        sys.add_prismatic_body(static_frame3dp("spindle_z"),
                               make_cuboid_body(m, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        size_t const zj = sys.index_of("spindle_z");
        sys.set_joint_spring_damper(zj, k_spring, c);
        // driven wheel spin about e1: a negligible-mass kinematic flywheel (1e-6 kg, just
        // enough for a non-singular inertia) -- it supplies the relative-sliding speed
        // n_w*R_w without materially loading the axial DOF (M[0] stays ~ m).
        sys.add_revolute_body(static_frame3dp("wheel"),
                              make_cuboid_body(1.0e-6, R_w, R_w, 0.01), O_3dp,
                              vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, zj);
        size_t const wj = sys.index_of("wheel");
        sys.set_driven_rate(wj, n_w);

        // register the grinding contact: rim point on the wheel vs the wafer, depth law
        vec3dp const rim_b{0.0, R_w, 0.0, 1.0}; // rim point in the wheel body frame
        sys.set_contact_force(
            wj, sys.index_of("wafer"),
            [=](contact_state3dp const& s) {
                return grinding_force_depth(s, k_grind, mu);
            },
            vec3dp{1.0, 0.0, 0.0, 0.0}, delta0, 30.0e-6 / 60.0, 0.8e-3, rim_b);

        // GATE 1: at the nominal pose (z = 0) the only force is the contact normal force
        // k_grind*delta0 = 100 N on the mass m -> initial axial accel 125 m/s^2
        CHECK(sys.mass_matrix()[0] == doctest::Approx(m));
        CHECK(sys.joint_accel(zj) == doctest::Approx(F_ff / m));
        // GATE 2: the contact has a real relative-sliding speed (wheel-dominated) ==
        // n_w*R_w
        vec3dp const P0 = unitize(move3dp(rim_b, sys.get_pos_trafo(wj, 0)));
        CHECK(to_val(bulk_nrm(sys.point_velocity(P0, wj))) == doctest::Approx(n_w * R_w));

        // integrate to steady state (zeta_eff = 0.7 -> settles in a few periods)
        value_t const dt = 1.0e-6;
        for (size_t nstep = 1; nstep <= 20000; ++nstep) // 20 ms
            sys.step(dt);

        value_t const z_ss = sys.joint_phi(zj);
        value_t const F_ss = sys.contact_force(wj).x; // axial component (n_hat = e1)
        value_t const d_ss = sys.contact_engagement(wj);
        fmt::println("  feed-forward nominal force : {:>7.3f} N  (k_grind*delta0)", F_ff);
        fmt::println("  closed-loop equilibrium    : z_eq = {:.4f} um, F_eq = {:>7.3f} N "
                     "(analytic)",
                     z_eq * 1.0e6, F_eq);
        fmt::println(
            "  simulated steady state     : z_ss = {:.4f} um, F_ss = {:>7.3f} N, "
            "delta_ss = {:.4f} um",
            z_ss * 1.0e6, F_ss, d_ss * 1.0e6);
        fmt::println("  -> spindle compliance backs the wheel out of the cut: "
                     "{:.1f} N -> {:.1f} N ({:+.1f}%)",
                     F_ff, F_ss, 100.0 * (F_ss - F_ff) / F_ff);

        // GATE 3: the closed loop settles to the analytic equilibrium
        CHECK(z_ss == doctest::Approx(z_eq).epsilon(1e-4));
        CHECK(F_ss == doctest::Approx(F_eq).epsilon(1e-4));
        CHECK(d_ss == doctest::Approx(delta0 - z_eq).epsilon(1e-4));
        // GATE 4: closing the loop reduced the force below the feed-forward nominal
        CHECK(F_ss < F_ff);
        CHECK(F_eq == doctest::Approx(F_ff * k_spring / k_tot)); // == feed-forward as
                                                                 // k_spring -> infinity

        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding force loop (Phase D.2b)")

// =================================================================================
// PGA3DP dynamic_system3dp - Phase D.2c: wafer thinning + optional chuck-side reaction
//
// Two pieces of infrastructure on top of the D.2b loop:
//  (1) WAFER THINNING -- a quasi-static wafer thickness decremented between feed
//  macro-steps
//      by a SIMPLE removal law MRR = k_mrr * delta * v_rel (one tunable constant k_mrr,
//      lumping contact width / removal efficiency). The thickness is exposed in
//      contact_state so a thickness-dependent force law or the D.2e feed control can read
//      it.
//  (2) CHUCK-SIDE REACTION -- the optional, user-configurable equal/opposite force line
//      -wdg(P, F) on the wafer frame (set_contact_reaction). DEFAULT off = rigid wafer;
//      the realistic reaction physics is a later refinement, so this just builds + tests
//      the infrastructure (Newton's third law for a simple two-DOF case).
// =================================================================================
TEST_SUITE("PGA3DP: grinding wafer thinning + reaction (Phase D.2c)")
{

    TEST_CASE("pga3dp: quasi-static wafer thinning from MRR (Phase D.2c)")
    {
        fmt::println("");
        fmt::println("pga3dp: quasi-static wafer thinning from MRR (Phase D.2c)");
        fmt::println("");

        // the D.2b axial loop + a removal constant k_mrr. The wafer thins
        // quasi-statically (held constant within the kHz dynamics, decremented between
        // feed macro-steps).
        value_t const m = 0.8, R_w = 0.15, n_w = rpm2radps(3000.0);
        value_t const delta0 = 1.0e-6, k_grind = 1.0e8;
        value_t const mu = std::sqrt(2.0) * 25.0 / 100.0;
        value_t const k_spring = 9.5e8, k_tot = k_spring + k_grind;
        value_t const c = 2.0 * 0.7 * std::sqrt(k_tot * m);
        value_t const tw0 = 0.8e-3; // initial wafer thickness [m]
        // tunable removal constant [1/m]: MRR = k_mrr*delta*v_rel is the
        // thickness-reduction RATE [m/s]; chosen so the thinning is ~ Tao's feed scale
        // (30 um/min ~ 0.4 um/s).
        value_t const k_mrr = 0.01;

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
        sys.add_frame(static_frame3dp("W"));
        sys.add_frame(static_frame3dp("wafer"), sys.index_of("W"));
        sys.add_prismatic_body(static_frame3dp("spindle_z"),
                               make_cuboid_body(m, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        size_t const zj = sys.index_of("spindle_z");
        sys.set_joint_spring_damper(zj, k_spring, c);
        sys.add_revolute_body(static_frame3dp("wheel"),
                              make_cuboid_body(1.0e-6, R_w, R_w, 0.01), O_3dp,
                              vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, zj);
        size_t const wj = sys.index_of("wheel");
        sys.set_driven_rate(wj, n_w);
        vec3dp const rim_b{0.0, R_w, 0.0, 1.0};
        sys.set_contact_force(
            wj, sys.index_of("wafer"),
            [=](contact_state3dp const& s) {
                return grinding_force_depth(s, k_grind, mu);
            },
            vec3dp{1.0, 0.0, 0.0, 0.0}, delta0, 30.0e-6 / 60.0, tw0, rim_b, k_mrr);

        // settle to the axial equilibrium, then thin quasi-statically
        value_t const dt = 1.0e-6;
        for (size_t nstep = 1; nstep <= 20000; ++nstep)
            sys.step(dt);

        value_t const d_eq = sys.contact_engagement(wj);
        value_t const rate = sys.removal_rate(wj); // MRR = k_mrr * delta_eq * v_rel
        CHECK(rate == doctest::Approx(k_mrr * d_eq * n_w *
                                      R_w)); // v_rel = n_w*R_w (static wafer)
        fmt::println("  delta_eq = {:.4f} um, v_rel = {:.2f} m/s, MRR = {:.4e} m/s",
                     d_eq * 1e6, n_w * R_w, rate);

        value_t const dt_macro = 1.0; // feed macro-step [s]
        value_t prev = sys.wafer_thickness(wj);
        CHECK(prev == doctest::Approx(tw0));
        fmt::println("  wafer thickness over 5 feed macro-steps (dt = {:.2f} s):",
                     dt_macro);
        for (int i = 0; i < 5; ++i) {
            sys.update_wafer_thinning(dt_macro);
            value_t const tw = sys.wafer_thickness(wj);
            CHECK(tw < prev);                                     // monotone decrease
            CHECK(prev - tw == doctest::Approx(rate * dt_macro)); // exact decrement
            fmt::println("    step {}: tw = {:.3f} um  (removed {:.4f} um)", i + 1,
                         tw * 1e6, (prev - tw) * 1e6);
            prev = tw;
        }
        fmt::println("");
    }

    TEST_CASE("pga3dp: optional chuck-side reaction wrench (Phase D.2c)")
    {
        fmt::println("");
        fmt::println("pga3dp: optional chuck-side reaction wrench (Phase D.2c)");
        fmt::println("");

        // Newton's third law for the contact: with the reaction ON, the equal/opposite
        // force line -wdg(P, F) drives a COMPLIANT wafer DOF opposite the spindle; with
        // it OFF (default rigid wafer) the wafer stays put. Two prismatic DOFs along e1
        // (spindle mass m1, wafer mass m2) facing across the contact; at t = 0 (both at
        // rest, mu = 0, delta = delta0) the contact normal force is F = k_grind*delta0 =
        // 100 N.
        value_t const m1 = 0.8, m2 = 2.0;
        value_t const delta0 = 1.0e-6, k_grind = 1.0e8,
                      mu = 0.0;             // pure normal force here
        value_t const F = k_grind * delta0; // = 100 N

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
        sys.add_frame(static_frame3dp("W"));
        sys.add_prismatic_body(static_frame3dp("spindle"),
                               make_cuboid_body(m1, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        sys.add_prismatic_body(static_frame3dp("wafer"),
                               make_cuboid_body(m2, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        size_t const si = sys.index_of("spindle"), wi = sys.index_of("wafer");
        sys.set_contact_force(
            si, wi,
            [=](contact_state3dp const& s) {
                return grinding_force_depth(s, k_grind, mu);
            },
            vec3dp{1.0, 0.0, 0.0, 0.0}, delta0);

        // reaction OFF (default rigid wafer): only the spindle accelerates
        CHECK(sys.joint_accel(si) == doctest::Approx(F / m1));
        CHECK(sys.joint_accel(wi) == doctest::Approx(0.0));

        // reaction ON: the wafer accelerates opposite, ratio -m1/m2 (Newton's third law)
        sys.set_contact_reaction(si, true);
        value_t const a_s = sys.joint_accel(si), a_w = sys.joint_accel(wi);
        CHECK(a_s == doctest::Approx(F / m1));
        CHECK(a_w == doctest::Approx(-F / m2));
        CHECK(a_w / a_s == doctest::Approx(-m1 / m2));
        fmt::println(
            "  reaction ON : spindle accel = {:>8.2f}, wafer accel = {:>8.2f} m/s^2"
            "  (ratio {:.3f} = -m1/m2)",
            a_s, a_w, a_w / a_s);

        // toggling back OFF restores the rigid-wafer behaviour
        sys.set_contact_reaction(si, false);
        CHECK(sys.joint_accel(wi) == doctest::Approx(0.0));
        fmt::println("  reaction OFF: wafer accel = 0 (rigid wafer)");
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding wafer thinning + reaction (Phase D.2c)")

// =================================================================================
// PGA3DP dynamic_system3dp - Phase D.2d-2: integrator selector on the grinding loop
//
// dynamic_system3dp::set_integrator(rk4 | abm2) picks the time integrator for
// coupled_step (RK4 default, byte-identical to the old loop; ABM2 =
// Adams-Bashforth-Moulton 2nd order, a drop-in sharing the same forward-dynamics rhs).
// Run the D.2b axial grinding loop with both and (1) confirm both reach the analytic
// equilibrium and agree, and (2) MEASURE whether the loop is actually stiff -- the reason
// this step comes before any adaptive/implicit work (TODO/grinding.md Phase D.2d).
// Verdict: the working dt sits far inside both integrators' stability regions, so dt is
// ACCURACY-bound, not stability-bound -> NOT classically stiff at the current contact
// stiffness; an implicit solver (D.2d-4) is only needed if a near-rigid contact penalty
// injects a mode far above the ~6 kHz physics band.
// =================================================================================
TEST_SUITE("PGA3DP: grinding loop integrator (Phase D.2d-2)")
{

    TEST_CASE(
        "pga3dp: RK4 vs ABM2 on the grinding loop + stiffness verdict (Phase D.2d-2)")
    {
        fmt::println("");
        fmt::println("pga3dp: RK4 vs ABM2 on the grinding loop (Phase D.2d-2)");
        fmt::println("");

        value_t const m = 0.8, R_w = 0.15, n_w = rpm2radps(3000.0);
        value_t const delta0 = 1.0e-6, k_grind = 1.0e8;
        value_t const mu = std::sqrt(2.0) * 25.0 / 100.0;
        value_t const k_spring = 9.5e8, k_tot = k_spring + k_grind;
        value_t const c = 2.0 * 0.7 * std::sqrt(k_tot * m);
        value_t const z_eq = k_grind * delta0 / k_tot;
        value_t const F_eq = k_grind * (delta0 - z_eq); // analytic equilibrium force

        auto build = [&](integrator_kind ik) {
            auto sys = std::make_unique<dynamic_system3dp>();
            sys->set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            sys->add_frame(static_frame3dp("W"));
            sys->add_frame(static_frame3dp("wafer"), sys->index_of("W"));
            sys->add_prismatic_body(
                static_frame3dp("spindle_z"), make_cuboid_body(m, 0.05, 0.05, 0.05),
                vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys->index_of("W"));
            size_t const zj = sys->index_of("spindle_z");
            sys->set_joint_spring_damper(zj, k_spring, c);
            sys->add_revolute_body(static_frame3dp("wheel"),
                                   make_cuboid_body(1.0e-6, R_w, R_w, 0.01), O_3dp,
                                   vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, zj);
            size_t const wj = sys->index_of("wheel");
            sys->set_driven_rate(wj, n_w);
            sys->set_contact_force(
                wj, sys->index_of("wafer"),
                [=](contact_state3dp const& s) {
                    return grinding_force_depth(s, k_grind, mu);
                },
                vec3dp{1.0, 0.0, 0.0, 0.0}, delta0, 30.0e-6 / 60.0, 0.8e-3,
                vec3dp{0.0, R_w, 0.0, 1.0});
            sys->set_integrator(ik);
            return sys;
        };

        value_t const dt = 1.0e-6;
        size_t const N = 20000;
        auto run = [&](integrator_kind ik) {
            auto sys = build(ik);
            for (size_t i = 0; i < N; ++i)
                sys->step(dt);
            return sys->contact_force(sys->index_of("wheel")).x; // axial force
        };
        value_t const F_rk = run(integrator_kind::rk4);
        value_t const F_ab = run(integrator_kind::abm2);

        fmt::println("  analytic equilibrium F_eq = {:.4f} N", F_eq);
        fmt::println("  RK4  : F = {:.4f} N", F_rk);
        fmt::println("  ABM2 : F = {:.4f} N   (|F_rk - F_ab| / F_eq = {:.2e})", F_ab,
                     std::abs(F_rk - F_ab) / F_eq);

        // both integrators reach the analytic equilibrium and agree to <= 1% (D.2d gate)
        CHECK(F_rk == doctest::Approx(F_eq).epsilon(1e-4));
        CHECK(F_ab == doctest::Approx(F_eq).epsilon(1e-3)); // ABM2 lower order, looser
        CHECK(std::abs(F_rk - F_ab) / F_eq < 0.01);

        // stiffness verdict: compare the working dt to the explicit-stability dt limits.
        // If dt is far below them, stability is NOT the binding constraint (dt is
        // accuracy-set to resolve the ~6 kHz vibration) -> the loop is not classically
        // stiff.
        value_t const w_eff = std::sqrt(k_tot / m); // ~36229 rad/s
        value_t const dt_rk = 2.78 / w_eff;         // RK4 stability limit (order of mag.)
        value_t const dt_ab = 1.0 / w_eff;          // ABM2 PECE stability limit (smaller)
        fmt::println("  w_eff = {:.0f} rad/s; stability dt: RK4 ~{:.2e}, ABM2 ~{:.2e}; "
                     "working dt = {:.0e}",
                     w_eff, dt_rk, dt_ab, dt);
        fmt::println(
            "  -> working dt is {:.0f}x below the ABM2 stability limit: dt is "
            "ACCURACY-bound, the loop is NOT classically stiff (k_grind = {:.0e}).",
            dt_ab / dt, k_grind);
        CHECK(dt < 0.1 * dt_ab); // stability is not binding -> not stiff at this k_grind
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding loop integrator (Phase D.2d-2)")

// =================================================================================
// PGA3DP - Phase D.2e: grinding feed control (force limiting) + the throughput cost
//
// grinding_controller3dp is the optional control tier on top of the force loop: it holds
// the grinding force within a target envelope by commanding the infeed rate. This runs
// the PROCESS-level (macro-timescale) loop -- the quasi-static reduction of the fast
// force loop, in which the steady force is ~proportional to the feed (removal balance
// v_feed = MRR = k_mrr*delta*v_rel -> delta proportional to feed, F = k_grind*delta).
//
// The point of the test (and the answer to "is reducing the force to 70% just feed ->
// 70%?"): cutting the force target to a fraction cuts the feed to the same fraction AND
// the removal rate to the same fraction, so the wafer thins slower and the process takes
// ~1/fraction longer -- the throughput cost of force-limiting. Two schedules: (a)
// constant 70% cap; (b) thickness-scheduled taper 100% -> 70% between 150 and 80 um.
// =================================================================================
TEST_SUITE("PGA3DP: grinding feed control (Phase D.2e)")
{

    TEST_CASE("pga3dp: feed control -- force limiting + throughput cost (Phase D.2e)")
    {
        fmt::println("");
        fmt::println("pga3dp: grinding feed control -- force limiting + throughput cost "
                     "(Phase D.2e)");
        fmt::println("");

        value_t const F_ref = 100.0;          // reference grinding force [N] (Tao F_z)
        value_t const v_nom = 30.0e-6 / 60.0; // nominal feed 30 um/min -> [m/s]
        value_t const gain =
            F_ref / v_nom; // quasi-static F <-> feed gain (removal balance)
        value_t const tw0 = 200.0e-6, tw_end = 80.0e-6; // grind 200 um -> 80 um

        // step the slow process: feed = ctrl.feed_command(tw); steady F = gain*feed; the
        // removal rate equals the feed (steady balance) so tw -= feed*dt. Returns
        // {process time, peak force}.
        auto grind = [&](grinding_controller3dp const& ctrl) {
            value_t tw = tw0, t = 0.0, F_max = 0.0;
            value_t const dt = 1.0; // macro-step [s]
            while (tw > tw_end) {
                value_t const feed = ctrl.feed_command(tw);
                F_max = std::max(F_max, gain * feed);
                tw -= feed * dt; // removal rate = feed (steady-state)
                t += dt;
            }
            return std::pair{t, F_max};
        };

        grinding_controller3dp base{v_nom, F_ref, 1.0, 0.0, 0.0};      // no cap (100%)
        grinding_controller3dp capA{v_nom, F_ref, 0.7, 0.0, 0.0};      // (a) constant 70%
        grinding_controller3dp capB{v_nom, F_ref, 0.7, 150e-6, 80e-6}; // (b) taper
        auto const [t_base, F_base] = grind(base);
        auto const [t_a, F_a] = grind(capA);
        auto const [t_b, F_b] = grind(capB);

        fmt::println("  baseline (100%):    F = {:>5.1f} N, process = {:>4.0f} s", F_base,
                     t_base);
        fmt::println("  option a (70% cap): F = {:>5.1f} N, process = {:>4.0f} s  "
                     "({:.2f}x longer -- the throughput cost)",
                     F_a, t_a, t_a / t_base);
        fmt::println("  option b (taper):   F_max = {:>5.1f} N (full above 150um), "
                     "process = {:>4.0f} s ({:.2f}x)",
                     F_b, t_b, t_b / t_base);

        // option (a): force held at 70%, and the process is ~1/0.7 = 1.43x longer
        CHECK(F_a == doctest::Approx(0.7 * F_ref));
        CHECK(t_a / t_base == doctest::Approx(1.0 / 0.7).epsilon(0.01));
        // option (b): peak force is still 100% (above 150 um); the taper slows the late
        // phase -> longer than baseline but cheaper than the constant cap
        CHECK(F_b == doctest::Approx(F_ref));
        CHECK(t_b > t_base);
        CHECK(t_b < t_a);

        // the taper schedule itself: full above tw_hi, capped at/below tw_lo, linear
        // between
        CHECK(capB.target_fraction(200e-6) == doctest::Approx(1.0));  // above 150
        CHECK(capB.target_fraction(150e-6) == doctest::Approx(1.0));  // at 150
        CHECK(capB.target_fraction(115e-6) == doctest::Approx(0.85)); // midpoint -> 85%
        CHECK(capB.target_fraction(80e-6) == doctest::Approx(0.7));   // at 80 -> capped
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding feed control (Phase D.2e)")
