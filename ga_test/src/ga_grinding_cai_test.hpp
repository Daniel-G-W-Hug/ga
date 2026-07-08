// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <algorithm> // std::min, std::max, std::clamp
#include <array>     // std::array (Cai HTM 4x4 reference matrices, Phase F)
#include <cmath>     // std::sin, std::cos, std::atan
#include <limits>    // std::numeric_limits
#include <random>    // std::mt19937, distributions (grain-size randomness)
#include <utility>   // std::pair
#include <vector>    // std::vector

// include functions to be tested
#include "ga/ga_ega.hpp" // for cross product
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)
using namespace hd::ga::ega; // use specific operations of EGA (for cross product)

TEST_SUITE("PGA3DP: Cai volumetric-error machine")
{

    /////////////////////////////////////////////////////////////////////////////////////
    // Phase F of the wafer-grinding plan (TODO/grinding.md): the machine geometric /
    // VOLUMETRIC error model of Cai et al., Measurement 234 (2024) 114825 -- a third,
    // distinct error source from Tao (spindle dynamics) and Zhou (axis-tilt kinematics):
    // the deterministic machine-tool geometric error budget and the resulting volumetric
    // error at the functional point.
    //
    // F.0 -- the geometric on-ramp: the Cai two-branch kinematic chain (Fig. 2). The Cai
    // machine has two TRANSLATIONAL axes (X, Z) and two ROTATIONAL axes (C1 = grinding
    // wheel, angle bt; C2 = wafer, angle bw). Mapping machine axes to the GA basis:
    // machine X -> e1, Y -> e2, Z -> e3. The C1/C2 rotations are about e3; the wheel-tilt
    // angles ax, ay are rotations about e1 and e2 (Cai Eq.10, Tx and Ty). The functional
    // point is a grain on the wheel rim, Pt = (Rt cos(th), Rt sin(th), 0) in the tilted
    // wheel frame.
    //
    // The IDEAL transform from the wheel rim to the wafer frame is Cai Eq.(1):
    //
    //     Pi_t = Rz(-bw) . T(-X,0,0) . T(0,0,Z) . Rz(bt) . [Tx(ax).Ty(ay)] . Pt
    //
    // i.e. (wafer chain)^-1 . (wheel chain) . Pt, where the wafer chain is
    // bed -> X-translation -> C2-rotation and the wheel chain is
    // bed -> Z-translation -> C1-rotation -> wheel-tilt. This is exactly the LCA walk
    // get_pos_trafo(wheel_rim, wafer) on the frame tree below. The point of F.0: the GA
    // motor chain reproduces the paper's HTM matrix chain. The GATE computes Eq.(1) as an
    // INDEPENDENT 4x4 HTM product (the paper's own method) and checks it against the GA
    // result component-wise. (Cai Eq.1 linearises the tilt; we use the EXACT Tx.Ty of
    // Eq.10 on BOTH sides, so the gate is exact -- the small-angle tilt linearisation is
    // a Phase-F.1/F.2 concern, not here.)
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: Cai volumetric-error machine, ideal transform Eq.1 (Phase F.0)")
    {
        fmt::println("");
        fmt::println(
            "pga3dp: Cai volumetric-error machine, ideal transform Eq.1 (Phase F.0)");
        fmt::println("");

        // ---- independent 4x4 HTM reference (Cai's own method) -------------------------
        using mat4 = std::array<std::array<double, 4>, 4>;

        auto ident = []() {
            mat4 M{};
            for (int i = 0; i < 4; ++i)
                M[i][i] = 1.0;
            return M;
        };
        auto matmul = [](mat4 const& A, mat4 const& B) {
            mat4 C{};
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 4; ++k)
                        s += A[i][k] * B[k][j];
                    C[i][j] = s;
                }
            return C;
        };
        auto transl = [&](double dx, double dy, double dz) {
            mat4 M = ident();
            M[0][3] = dx;
            M[1][3] = dy;
            M[2][3] = dz;
            return M;
        };
        auto Rz = [&](double a) { // rotation about e3 (machine Z; C1/C2 axes)
            mat4 M = ident();
            M[0][0] = std::cos(a);
            M[0][1] = -std::sin(a);
            M[1][0] = std::sin(a);
            M[1][1] = std::cos(a);
            return M;
        };
        auto Rx = [&](double a) { // Cai Eq.10 Tx: rotation about e1 (machine X)
            mat4 M = ident();
            M[1][1] = std::cos(a);
            M[1][2] = -std::sin(a);
            M[2][1] = std::sin(a);
            M[2][2] = std::cos(a);
            return M;
        };
        auto Ry = [&](double a) { // Cai Eq.10 Ty: rotation about e2 (machine Y)
            mat4 M = ident();
            M[0][0] = std::cos(a);
            M[0][2] = std::sin(a);
            M[2][0] = -std::sin(a);
            M[2][2] = std::cos(a);
            return M;
        };

        // ---- machine parameters (non-trivial sample; Rt from Table 2) -----------------
        double const Rt = 72.5;            // grinding-wheel radius [mm]
        double const X = 100.0;            // X-axis position [mm]
        double const Z = 50.0;             // Z-axis position [mm]
        double const bt = 0.30;            // C1 (wheel) angle [rad]
        double const bw = 0.70;            // C2 (wafer) angle [rad]
        double const ax = 0.01 * pi / 180; // wheel tilt about e1 (Table 2: 0.01 deg)
        double const ay = 0.01 * pi / 180; // wheel tilt about e2 (Table 2: 0.01 deg)

        // ---- GA frame tree: two branches off the machine bed --------------------------
        kinematic_system3dp sys;
        sys.add_frame(static_frame3dp("bed")); // inertial root (machine bed P)

        // wafer chain: bed -> X-axis translation -> C2 (wafer) rotation about e3
        sys.add_frame(static_frame3dp("x_axis", vec3dp{X, 0.0, 0.0, 1.0}),
                      sys.index_of("bed"));
        sys.add_frame(static_frame3dp("c2_wafer", O_3dp, vec3dp{0.0, 0.0, bw, 0.0}),
                      sys.index_of("x_axis"));

        // wheel chain: bed -> Z-axis translation -> C1 (wheel) rotation about e3 ->
        // wheel tilt Tx(ax) then Ty(ay) (two frames so each is a single axis*angle)
        sys.add_frame(static_frame3dp("z_axis", vec3dp{0.0, 0.0, Z, 1.0}),
                      sys.index_of("bed"));
        sys.add_frame(static_frame3dp("c1_wheel", O_3dp, vec3dp{0.0, 0.0, bt, 0.0}),
                      sys.index_of("z_axis"));
        sys.add_frame(static_frame3dp("tilt_x", O_3dp, vec3dp{ax, 0.0, 0.0, 0.0}),
                      sys.index_of("c1_wheel"));
        sys.add_frame(static_frame3dp("tilt_y", O_3dp, vec3dp{0.0, ay, 0.0, 0.0}),
                      sys.index_of("tilt_x"));

        // tilt = Tx(ax) . Ty(ay) (Cai Eq.10), shared by the HTM reference
        mat4 const tilt = matmul(Rx(ax), Ry(ay));

        // sweep the phase angle th of the rim grain over one wheel revolution
        int const nth = 12;
        for (int k = 0; k <= nth; ++k) {
            double const th = 2.0 * pi * double(k) / double(nth);

            // GA: the rim grain (in the tilted wheel frame) expressed in the wafer frame
            auto const Pt_local = vec3dp{Rt * std::cos(th), Rt * std::sin(th), 0.0, 1.0};
            auto const Pi_ga =
                unitize(move3dp(Pt_local, sys.get_pos_trafo("tilt_y", "c2_wafer")));

            // HTM reference (Cai Eq.1): Rz(-bw) T(-X) T(Z) Rz(bt) tilt Pt
            mat4 const chain = matmul(
                Rz(-bw), matmul(transl(-X, 0.0, 0.0),
                                matmul(transl(0.0, 0.0, Z), matmul(Rz(bt), tilt))));
            double const px = Rt * std::cos(th), py = Rt * std::sin(th), pz = 0.0;
            double const rx =
                chain[0][0] * px + chain[0][1] * py + chain[0][2] * pz + chain[0][3];
            double const ry =
                chain[1][0] * px + chain[1][1] * py + chain[1][2] * pz + chain[1][3];
            double const rz =
                chain[2][0] * px + chain[2][1] * py + chain[2][2] * pz + chain[2][3];

            CHECK(Pi_ga.x == doctest::Approx(rx).epsilon(1e-9));
            CHECK(Pi_ga.y == doctest::Approx(ry).epsilon(1e-9));
            CHECK(Pi_ga.z == doctest::Approx(rz).epsilon(1e-9));

            if (k % (nth / 6) == 0)
                fmt::println("th = {:>5.2f}: Pi_ga = ({:>9.4f},{:>9.4f},{:>9.4f}) | "
                             "HTM = ({:>9.4f},{:>9.4f},{:>9.4f})",
                             th, Pi_ga.x, Pi_ga.y, Pi_ga.z, rx, ry, rz);
        }

        // ---- smoke check: trivial machine (no offsets, no tilt, bt == bw) is identity -
        {
            kinematic_system3dp s0;
            s0.add_frame(static_frame3dp("bed"));
            s0.add_frame(static_frame3dp("x0", O_3dp), s0.index_of("bed"));
            s0.add_frame(static_frame3dp("c2_0", O_3dp, vec3dp{0.0, 0.0, 0.4, 0.0}),
                         s0.index_of("x0"));
            s0.add_frame(static_frame3dp("z0", O_3dp), s0.index_of("bed"));
            s0.add_frame(static_frame3dp("c1_0", O_3dp, vec3dp{0.0, 0.0, 0.4, 0.0}),
                         s0.index_of("z0"));
            auto const Pt0 = vec3dp{Rt, 0.0, 0.0, 1.0};
            auto const P0 = unitize(move3dp(Pt0, s0.get_pos_trafo("c1_0", "c2_0")));
            CHECK(P0.x == doctest::Approx(Rt)); // bt == bw, no offset -> Pt unchanged
            CHECK(P0.y == doctest::Approx(0.0));
            CHECK(P0.z == doctest::Approx(0.0));
        }

        fmt::println("");
        fmt::println(
            "Eq.(1) ideal transform: GA motor chain == HTM 4x4 product (Phase F.0)");
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // F.1 -- the error-motor primitive (the one genuinely new piece of Phase F). Each
    // geometric error of the Cai budget (Table 1) is a small rigid perturbation, given in
    // the paper as a LINEARISED 4x4 matrix -- a small rotation eps = (ex, ey, ez) plus a
    // translation del = (dx, dy, dz):
    //
    //     M_lin = [ 1   -ez   ey   dx ]
    //             [ ez   1   -ex   dy ]
    //             [-ey   ex   1    dz ]
    //             [ 0    0    0    1  ]
    //
    // In GA this is just a MOTOR: rotation about the local origin by eps, then
    // translation by del -- exactly motor_from_pose3dp(pose3dp{origin = del, rot = eps})
    // (no new library symbol needed; it reuses the pose<->motor converters). Three things
    // are shown:
    //
    //   (a) the error motor reproduces the paper's linearised error matrix (to O(eps^2),
    //       the regime where the linearisation is valid);
    //   (b) the error motor is the EXACT rigid map -- it preserves distances to machine
    //   eps,
    //       whereas the paper's linearised matrix does not (a GA advantage, for free);
    //   (c) the relative error motor M_err = M_actual (rgpr) rrev(M_ideal) isolates the
    //       inserted error, and pose3dp_from_motor / log() recover it as a SCREW carrying
    //       BOTH the orientation error eps AND the position error del in one object --
    //       where the paper's E = Pe - Pi is position-only.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: Cai geometric error as a motor (Phase F.1)")
    {
        fmt::println("");
        fmt::println("pga3dp: Cai geometric error as a motor (Phase F.1)");
        fmt::println("");

        double const as = pi / 180.0 / 3600.0; // 1 arcsec in radians

        // a representative PDGE (Z-axis-scale: del ~ um in mm, eps ~ arcsec)
        auto const del = vec3dp{-3.6e-3, -2.5e-3, 1.0e-3, 1.0};    // (dx, dy, dz) [mm]
        double const ex = 4.0 * as, ey = -3.0 * as, ez = 2.0 * as; // (ex, ey, ez) [rad]
        auto const eps_rot = vec3dp{ex, ey, ez, 0.0};

        // the error motor: rotate about the origin by eps, then translate by del
        auto const M_e = motor_from_pose3dp(pose3dp{del, eps_rot});

        // paper's linearised 4x4 error applied to a point (M_lin . p)
        auto lin_apply = [&](vec3dp const& p) {
            return vec3dp{p.x - ez * p.y + ey * p.z + del.x,
                          ez * p.x + p.y - ex * p.z + del.y,
                          -ey * p.x + ex * p.y + p.z + del.z, 1.0};
        };

        // test points: a few rim points plus off-axis points
        std::vector<vec3dp> const pts = {
            vec3dp{72.5, 0.0, 0.0, 1.0}, vec3dp{0.0, 72.5, 0.0, 1.0},
            vec3dp{50.0, -30.0, 20.0, 1.0}, vec3dp{-100.0, 40.0, 10.0, 1.0}};

        // (a) error motor == paper linearised matrix, to O(eps^2). With eps ~ arcsec the
        // gap is ~1e-9 mm over a ~150 mm part -- far below any geometric signal.
        double max_lin_gap = 0.0;
        for (auto const& p : pts) {
            auto const q_ga = unitize(move3dp(p, M_e));
            auto const q_lin = lin_apply(p);
            double const g = std::sqrt((q_ga.x - q_lin.x) * (q_ga.x - q_lin.x) +
                                       (q_ga.y - q_lin.y) * (q_ga.y - q_lin.y) +
                                       (q_ga.z - q_lin.z) * (q_ga.z - q_lin.z));
            max_lin_gap = std::max(max_lin_gap, g);
            CHECK(q_ga.x == doctest::Approx(q_lin.x).epsilon(1e-6));
            CHECK(q_ga.y == doctest::Approx(q_lin.y).epsilon(1e-6));
            CHECK(q_ga.z == doctest::Approx(q_lin.z).epsilon(1e-6));
        }
        fmt::println(
            "(a) max |motor - linearised| over test points = {:.3e} mm  (O(eps^2))",
            max_lin_gap);

        // (b) the motor is EXACTLY rigid (preserves distances); the linearised matrix is
        // not
        double motor_dist_err = 0.0, lin_dist_err = 0.0;
        for (size_t i = 0; i < pts.size(); ++i)
            for (size_t j = i + 1; j < pts.size(); ++j) {
                auto const d0 = bulk_nrm(pts[i] - pts[j]);
                auto const dm = bulk_nrm(unitize(move3dp(pts[i], M_e)) -
                                         unitize(move3dp(pts[j], M_e)));
                auto const dl = bulk_nrm(lin_apply(pts[i]) - lin_apply(pts[j]));
                motor_dist_err = std::max(motor_dist_err, std::abs(to_val(dm - d0)));
                lin_dist_err = std::max(lin_dist_err, std::abs(to_val(dl - d0)));
            }
        CHECK(motor_dist_err < 1e-9); // rigid to machine eps
        fmt::println("(b) max distance error: motor = {:.3e} mm (rigid), linearised = "
                     "{:.3e} mm (not)",
                     motor_dist_err, lin_dist_err);

        // (c) isolate + recover the error from a chain: M_err = M_actual (rgpr)
        // rrev(M_ideal)
        auto const M_ideal = motor_from_pose3dp(
            pose3dp{vec3dp{30.0, -10.0, 5.0, 1.0}, vec3dp{0.0, 0.0, 0.5, 0.0}});
        auto const M_actual =
            rgpr(M_e, M_ideal); // error applied in the ideal parent frame
        auto const M_err = rgpr(M_actual, rrev(M_ideal));

        // M_err recovers the inserted error motor (equal up to rgpr round-trip rounding)
        CHECK(std::abs(M_err.c0 - M_e.c0) < 1e-12);
        CHECK(std::abs(M_err.c1 - M_e.c1) < 1e-12);
        CHECK(std::abs(M_err.c2 - M_e.c2) < 1e-12);
        CHECK(std::abs(M_err.c3 - M_e.c3) < 1e-12);
        CHECK(std::abs(M_err.c4 - M_e.c4) < 1e-12);
        CHECK(std::abs(M_err.c5 - M_e.c5) < 1e-12);
        CHECK(std::abs(M_err.c6 - M_e.c6) < 1e-12);
        CHECK(std::abs(M_err.c7 - M_e.c7) < 1e-12);

        // decode the screw: pose3dp_from_motor gives (origin = del, rot = eps); log()
        // gives the full screw generator whose weight (rotation) part is eps -- ONE
        // object carrying BOTH position and orientation error
        auto const rec = pose3dp_from_motor(M_err);
        CHECK(rec.origin.x == doctest::Approx(del.x).epsilon(1e-6));
        CHECK(rec.origin.y == doctest::Approx(del.y).epsilon(1e-6));
        CHECK(rec.origin.z == doctest::Approx(del.z).epsilon(1e-6));
        CHECK(rec.rot.x == doctest::Approx(ex).epsilon(1e-6));
        CHECK(rec.rot.y == doctest::Approx(ey).epsilon(1e-6));
        CHECK(rec.rot.z == doctest::Approx(ez).epsilon(1e-6));

        auto const B = log(M_err); // screw generator; weight slots carry the rotation eps
        CHECK(B.vx == doctest::Approx(ex).epsilon(1e-5));
        CHECK(B.vy == doctest::Approx(ey).epsilon(1e-5));
        CHECK(B.vz == doctest::Approx(ez).epsilon(1e-5));
        fmt::println("(c) recovered error screw: rot = ({:.3e},{:.3e},{:.3e}) rad, "
                     "origin = ({:.3e},{:.3e},{:.3e}) mm",
                     rec.rot.x, rec.rot.y, rec.rot.z, rec.origin.x, rec.origin.y,
                     rec.origin.z);
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // F.2 -- the volumetric error map (Cai Eq.2/Eq.3). The full geometric-error budget
    // (Table 1: 24 PDGEs + 5 PIGEs) is composed into the two-branch chain. Cai Eq.(2)
    // inserts, per axis, a PIGE matrix (position-INdependent: squareness / parallelism)
    // and a PDGE matrix (position-dependent: per-axis translation + angular errors) into
    // the ideal chain of F.0:
    //
    //   wheel chain  M(bed<-wheel) = eta . T(Z) . PDGE_Z . PIGE_C1 . PDGE_C1 . Rz(bt) .
    //   tilt wafer chain  M(bed<-wafer) = T(X) . PDGE_X . PIGE_C2 . PDGE_C2 . Rz(bw)
    //   actual point Pe = M(bed<-wafer)^-1 . M(bed<-wheel) . Pt   (= get_pos_trafo walk)
    //
    // In GA each error is the F.1 error motor (motor_from_pose3dp{origin = del, rot =
    // eps}); the chain is just extra frames in the tree. The volumetric error is then E =
    // Pe - Pi (Cai Eq.3) -- which the paper hand-expands into a page-long closed form
    // that GA makes unnecessary: we compose motors and subtract the points.
    //
    // GATES: (1) the GA error chain reproduces an INDEPENDENT 4x4 HTM product of Eq.(2)
    // to O(eps^2) over a th sweep (full error set); (2) isolated-error spot checks
    // against the clean leading terms of Eq.(3): the relative yaw eps_zz/eps_zx gives the
    // in-plane E ~ Rt(eps_zz - eps_zx) term, an axial PDGE delta_zc1 gives E_z =
    // delta_zc1, etc. (The full Eq.3 closed form is the hand-expansion GA replaces; its
    // known Eq.3<->Eq.8 sign discrepancies are exactly why hand-expansion is
    // error-prone.)
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: Cai volumetric error map Eq.2/Eq.3 (Phase F.2)")
    {
        fmt::println("");
        fmt::println("pga3dp: Cai volumetric error map Eq.2/Eq.3 (Phase F.2)");
        fmt::println("");

        // machine parameters (as in F.0)
        double const Rt = 72.5, X = 100.0, Z = 50.0, bt = 0.30, bw = 0.70;
        double const ax = 0.01 * pi / 180, ay = 0.01 * pi / 180;
        double const as = pi / 180.0 / 3600.0; // 1 arcsec [rad]

        // the full Cai error budget (Table 1). Angles in arcsec, translations in um (mm).
        struct cai_err {
            // X-axis PDGE
            double dxx = 0, dyx = 0, dzx = 0, exx = 0, eyx = 0, ezx = 0;
            // Z-axis PDGE
            double dxz = 0, dyz = 0, dzz = 0, exz = 0, eyz = 0, ezz = 0;
            // C1-axis PDGE (ezc1 absorbed into bt)
            double dxc1 = 0, dyc1 = 0, dzc1 = 0, exc1 = 0, eyc1 = 0;
            // C2-axis PDGE (ezc2 absorbed into bw)
            double dxc2 = 0, dyc2 = 0, dzc2 = 0, exc2 = 0, eyc2 = 0;
            // PIGE: squareness X-Z, parallelism C1, parallelism C2
            double eta = 0, azoc1 = 0, bzoc1 = 0, azoc2 = 0, bzoc2 = 0;
        };

        // --- GA functional point Pe for a given error set + phase angle th
        // --------------
        auto Pe_ga = [&](cai_err const& e, double th) {
            kinematic_system3dp s;
            s.add_frame(static_frame3dp("bed"));
            // wheel chain: eta . T(Z) . PDGE_Z . PIGE_C1 . PDGE_C1 . Rz(bt) . tilt
            s.add_frame(static_frame3dp("z_pige", O_3dp, vec3dp{0, e.eta, 0, 0}),
                        s.index_of("bed"));
            s.add_frame(static_frame3dp("z_axis", vec3dp{0, 0, Z, 1}),
                        s.index_of("z_pige"));
            s.add_frame(static_frame3dp("z_pdge", vec3dp{e.dxz, e.dyz, e.dzz, 1},
                                        vec3dp{e.exz, e.eyz, e.ezz, 0}),
                        s.index_of("z_axis"));
            s.add_frame(static_frame3dp("c1_pige", O_3dp, vec3dp{e.azoc1, e.bzoc1, 0, 0}),
                        s.index_of("z_pdge"));
            s.add_frame(static_frame3dp("c1_pdge", vec3dp{e.dxc1, e.dyc1, e.dzc1, 1},
                                        vec3dp{e.exc1, e.eyc1, 0, 0}),
                        s.index_of("c1_pige"));
            s.add_frame(static_frame3dp("c1_rot", O_3dp, vec3dp{0, 0, bt, 0}),
                        s.index_of("c1_pdge"));
            s.add_frame(static_frame3dp("tilt_x", O_3dp, vec3dp{ax, 0, 0, 0}),
                        s.index_of("c1_rot"));
            s.add_frame(static_frame3dp("tilt_y", O_3dp, vec3dp{0, ay, 0, 0}),
                        s.index_of("tilt_x"));
            // wafer chain: T(X) . PDGE_X . PIGE_C2 . PDGE_C2 . Rz(bw). The C2 errors are
            // entered in NATURAL (non-negated) form: get_pos_trafo inverts the whole
            // wafer branch automatically, which reproduces Eq.(3)'s wafer-side signs
            // (e.g. -delta_zc2, -delta_zx). (Eq.2 prints the C2 matrices pre-negated
            // because the paper distributed the inverse by hand; we must not double-apply
            // it.)
            s.add_frame(static_frame3dp("x_axis", vec3dp{X, 0, 0, 1}), s.index_of("bed"));
            s.add_frame(static_frame3dp("x_pdge", vec3dp{e.dxx, e.dyx, e.dzx, 1},
                                        vec3dp{e.exx, e.eyx, e.ezx, 0}),
                        s.index_of("x_axis"));
            s.add_frame(static_frame3dp("c2_pige", O_3dp, vec3dp{e.azoc2, e.bzoc2, 0, 0}),
                        s.index_of("x_pdge"));
            s.add_frame(static_frame3dp("c2_pdge", vec3dp{e.dxc2, e.dyc2, e.dzc2, 1},
                                        vec3dp{e.exc2, e.eyc2, 0, 0}),
                        s.index_of("c2_pige"));
            s.add_frame(static_frame3dp("c2_rot", O_3dp, vec3dp{0, 0, bw, 0}),
                        s.index_of("c2_pdge"));
            auto const Pt = vec3dp{Rt * std::cos(th), Rt * std::sin(th), 0.0, 1.0};
            return unitize(move3dp(Pt, s.get_pos_trafo("tilt_y", "c2_rot")));
        };

        // --- independent 4x4 HTM reference (Cai's own method)
        // ---------------------------
        using mat4 = std::array<std::array<double, 4>, 4>;
        auto ident = []() {
            mat4 M{};
            for (int i = 0; i < 4; ++i)
                M[i][i] = 1.0;
            return M;
        };
        auto mul = [](mat4 const& A, mat4 const& B) {
            mat4 C{};
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 4; ++k)
                        s += A[i][k] * B[k][j];
                    C[i][j] = s;
                }
            return C;
        };
        auto transl = [&](double dx, double dy, double dz) {
            mat4 M = ident();
            M[0][3] = dx, M[1][3] = dy, M[2][3] = dz;
            return M;
        };
        auto Rz = [&](double a) {
            mat4 M = ident();
            M[0][0] = std::cos(a), M[0][1] = -std::sin(a), M[1][0] = std::sin(a),
            M[1][1] = std::cos(a);
            return M;
        };
        auto Rx = [&](double a) {
            mat4 M = ident();
            M[1][1] = std::cos(a), M[1][2] = -std::sin(a), M[2][1] = std::sin(a),
            M[2][2] = std::cos(a);
            return M;
        };
        auto Ry = [&](double a) {
            mat4 M = ident();
            M[0][0] = std::cos(a), M[0][2] = std::sin(a), M[2][0] = -std::sin(a),
            M[2][2] = std::cos(a);
            return M;
        };
        // linearised geometric-error matrix (eps = rotation, del = translation)
        auto emat = [&](double dx, double dy, double dz, double ex, double ey,
                        double ez) {
            mat4 M = ident();
            M[0][1] = -ez, M[0][2] = ey, M[1][0] = ez, M[1][2] = -ex, M[2][0] = -ey,
            M[2][1] = ex;
            M[0][3] = dx, M[1][3] = dy, M[2][3] = dz;
            return M;
        };
        auto rigid_inv = [&](mat4 const& M) { // [R|t]^-1 = [R^T | -R^T t]
            mat4 I = ident();
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    I[i][j] = M[j][i];
            for (int i = 0; i < 3; ++i)
                I[i][3] = -(I[i][0] * M[0][3] + I[i][1] * M[1][3] + I[i][2] * M[2][3]);
            return I;
        };
        auto Pe_htm = [&](cai_err const& e, double th) {
            mat4 const tilt = mul(Rx(ax), Ry(ay));
            mat4 K = mul(emat(0, 0, 0, 0, e.eta, 0),
                         mul(transl(0, 0, Z),
                             mul(emat(e.dxz, e.dyz, e.dzz, e.exz, e.eyz, e.ezz),
                                 mul(emat(0, 0, 0, e.azoc1, e.bzoc1, 0),
                                     mul(emat(e.dxc1, e.dyc1, e.dzc1, e.exc1, e.eyc1, 0),
                                         mul(Rz(bt), tilt))))));
            mat4 W = mul(
                transl(X, 0, 0),
                mul(emat(e.dxx, e.dyx, e.dzx, e.exx, e.eyx, e.ezx),
                    mul(emat(0, 0, 0, e.azoc2, e.bzoc2, 0),
                        mul(emat(e.dxc2, e.dyc2, e.dzc2, e.exc2, e.eyc2, 0), Rz(bw)))));
            mat4 const T = mul(rigid_inv(W), K);
            double const px = Rt * std::cos(th), py = Rt * std::sin(th);
            return std::array<double, 3>{T[0][0] * px + T[0][1] * py + T[0][3],
                                         T[1][0] * px + T[1][1] * py + T[1][3],
                                         T[2][0] * px + T[2][1] * py + T[2][3]};
        };

        // representative full error set
        cai_err full{.dxx = 2e-3,
                     .dyx = 1e-3,
                     .dzx = 1.5e-3,
                     .exx = 3 * as,
                     .eyx = 2 * as,
                     .ezx = 1 * as,
                     .dxz = -3.6e-3,
                     .dyz = -2.5e-3,
                     .dzz = 1e-3,
                     .exz = 2 * as,
                     .eyz = -1.5 * as,
                     .ezz = 1.5 * as,
                     .dxc1 = 1e-3,
                     .dyc1 = -1e-3,
                     .dzc1 = 2e-3,
                     .exc1 = 1 * as,
                     .eyc1 = 1.5 * as,
                     .dxc2 = -1e-3,
                     .dyc2 = 2e-3,
                     .dzc2 = -1.5e-3,
                     .exc2 = 2 * as,
                     .eyc2 = -1 * as,
                     .eta = 2 * as,
                     .azoc1 = 1.5 * as,
                     .bzoc1 = 1 * as,
                     .azoc2 = -1 * as,
                     .bzoc2 = 2 * as};

        // (1) GA error chain == independent HTM product of Eq.(2), to O(eps^2)
        double max_gap = 0.0;
        int const nth = 12;
        for (int k = 0; k <= nth; ++k) {
            double const th = 2.0 * pi * double(k) / double(nth);
            auto const g = Pe_ga(full, th);
            auto const h = Pe_htm(full, th);
            max_gap = std::max({max_gap, std::abs(g.x - h[0]), std::abs(g.y - h[1]),
                                std::abs(g.z - h[2])});
            CHECK(g.x == doctest::Approx(h[0]).epsilon(1e-6));
            CHECK(g.y == doctest::Approx(h[1]).epsilon(1e-6));
            CHECK(g.z == doctest::Approx(h[2]).epsilon(1e-6));
        }
        fmt::println("(1) GA chain vs HTM Eq.(2): max gap = {:.3e} mm  (O(eps^2))",
                     max_gap);

        // the volumetric error E = Pe - Pi (Cai Eq.3); print its magnitude over th
        cai_err const none{};
        double e_max = 0.0;
        for (int k = 0; k <= nth; ++k) {
            double const th = 2.0 * pi * double(k) / double(nth);
            auto const pe = Pe_ga(full, th);
            auto const pi_ = Pe_ga(none, th);
            double const ex = pe.x - pi_.x, ey = pe.y - pi_.y, ez = pe.z - pi_.z;
            e_max = std::max(e_max, std::sqrt(ex * ex + ey * ey + ez * ez));
            if (k % (nth / 4) == 0)
                fmt::println("th={:>5.2f}: E = ({:>9.4f}, {:>9.4f}, {:>9.4f}) um", th,
                             1e3 * ex, 1e3 * ey, 1e3 * ez);
        }
        CHECK(e_max < 0.05); // volumetric error stays sub-50-um for arcsec/um inputs
        fmt::println("    |E|_max = {:.4f} um", 1e3 * e_max);

        // (2a) Eq.(3) leading term: isolated yaw eps_zz (Z-axis) & eps_zx (X-axis) give
        // the in-plane E ~ -Rt(eps_zz - eps_zx) sin(th+bt-bw), +Rt(...) cos(th+bt-bw)
        {
            cai_err iso{};
            iso.ezz = 5 * as;
            iso.ezx = 2 * as;
            double const drel = (iso.ezz - iso.ezx);
            for (int k = 0; k <= nth; ++k) {
                double const th = 2.0 * pi * double(k) / double(nth);
                auto const pe = Pe_ga(iso, th);
                auto const pi_ = Pe_ga(none, th);
                double const ex = pe.x - pi_.x, ey = pe.y - pi_.y;
                CHECK(ex ==
                      doctest::Approx(-Rt * drel * std::sin(th + bt - bw)).epsilon(2e-3));
                CHECK(ey ==
                      doctest::Approx(Rt * drel * std::cos(th + bt - bw)).epsilon(2e-3));
            }
            fmt::println(
                "(2a) isolated yaw eps_zz/eps_zx: E_xy == Rt(eps_zz-eps_zx) term "
                "of Eq.(3)  [OK]");
        }

        // (2b) Eq.(3): isolated axial PDGEs give constant E_z = dzc1 + dzz - dzx - dzc2
        {
            cai_err iso{};
            iso.dzc1 = 2e-3, iso.dzz = 1e-3, iso.dzx = 1.5e-3, iso.dzc2 = -1e-3;
            double const ez_pred = iso.dzc1 + iso.dzz - iso.dzx - iso.dzc2;
            for (int k = 0; k <= nth; ++k) {
                double const th = 2.0 * pi * double(k) / double(nth);
                double const ez = Pe_ga(iso, th).z - Pe_ga(none, th).z;
                CHECK(ez == doctest::Approx(ez_pred).epsilon(1e-3));
            }
            fmt::println(
                "(2b) isolated axial PDGEs: E_z == dzc1+dzz-dzx-dzc2 = {:.4f} um "
                "(const) [OK]",
                1e3 * ez_pred);
        }

        fmt::println("");
        fmt::println("Volumetric map: GA motor composition reproduces Cai Eq.(2)/Eq.(3) "
                     "(Phase F.2)");
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // F.3 -- the Abbe/Bryan-optimized volumetric error (Cai Eq.5-8). The traditional
    // model (F.2) applies each axis's geometric error AS IF measured at the functional
    // point Pt. In reality each error is measured at a point OFFSET from Pt by L(j) (Fig.
    // 3/4): an Abbe offset for the positioning error, a Bryan offset for the straightness
    // errors. Cai's optimization adds, per axis j, the extra volumetric error that this
    // offset introduces:
    //
    //     E(j) = Tc2 . EM(j) . L(j)                                       (Cai Eq.5-7)
    //
    // where Tc2 = Rz(bw) (the common wafer rotation), L(j) is the measurement-point
    // position (Pt's nominal position plus the offset L), and EM(j) is the axis-j error
    // as a DIFFERENTIAL (skew) matrix -- the linearised error matrix with the identity
    // removed (zero diagonal: the eps skew + the del column). EM(j).L(j) is therefore the
    // FIRST-ORDER DISPLACEMENT of the point L(j) under axis j's error.
    //
    // GA reading (the showcase): EM(j).L(j) is exactly the error twist's VELOCITY FIELD
    // at L(j). The error motor of F.1 (M_e = motor_from_pose3dp{del, eps}) has screw
    // generator twist = 2*log(M_e) (the Omega = 2 Mdot (rcmt) rrev(M) convention; log
    // returns the half-angle generator because motor_from_pose3dp builds exp(0.5*B)). The
    // rate of change of a point under that twist is Xdot = rcmt(twist, X) -- so
    //
    //     E(j) = move3dp( rcmt(2*log(M_e_j), L(j)) , Rz(bw)-motor )
    //
    // evaluates the error twist at the OFFSET point and rotates it into the wafer frame.
    // The Abbe/Bryan cross terms that Eq.(8) carries and Eq.(3) lacks (eps.L lever arms
    // like eyz.Lzz, -ezz.Lyz, eyc1.Lzc1) are exactly the omega x L part of that velocity
    // field -- ONE operation per axis, no per-axis bespoke equation. The full optimized
    // error is E' = E(F.2) + sum_j E(j) over j = {z, x, c1, c2} (PDGE offsets) and {xoz,
    // zoc1, zoc2} (the PIGE-derived squareness/parallelism, evaluated at Pt's nominal
    // position with no extra offset).
    //
    // GATES: (1) the GA Abbe/Bryan correction sum reproduces an INDEPENDENT 4x4 HTM of
    // Eqs.(5-7) over a th sweep; (2) the FULL optimized error E' (GA base F.2 + GA
    // correction) == the assembled HTM of Eq.(8) to O(eps^2); (3) an isolated-offset spot
    // check (only eyz + Lzz) reproduces the eyz.Lzz cross term of Eq.(8)'s E'_x.
    //
    // PAPER-ERROR FLAG (resolved by the GA): Eq.(8)'s printed E'_x cbw-brace has the
    // straightness combination delta_xz - delta_xx, but its E'_y sbw-brace prints
    // delta_xx - delta_xz -- the SAME brace with the sign flipped (Eq.(3) prints
    // delta_xz - delta_xx in BOTH). In the GA both components come from rotating ONE
    // displacement vector by ONE Rz(bw), so they MUST share the brace: the GA confirms
    // Eq.(3)'s (delta_xz - delta_xx) and flags Eq.(8)'s E'_y print as a transcription
    // typo (check 4 below).
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: Cai Abbe/Bryan-optimized volumetric error Eq.5-8 (Phase F.3)")
    {
        fmt::println("");
        fmt::println(
            "pga3dp: Cai Abbe/Bryan-optimized volumetric error Eq.5-8 (Phase F.3)");
        fmt::println("");

        // machine parameters (as in F.0/F.2)
        double const Rt = 72.5, X = 100.0, Z = 50.0, bt = 0.30, bw = 0.70;
        double const ax = 0.01 * pi / 180, ay = 0.01 * pi / 180;
        double const as = pi / 180.0 / 3600.0; // 1 arcsec [rad]

        // the full Cai error budget (Table 1) PLUS the Abbe/Bryan measurement offsets
        // L(j) (Eq.4-7): Lsj is the s-directional distance from the j-axis measurement
        // point to the grinding-wheel center (mm). PIGE-derived terms (xoz, zoc1, zoc2)
        // carry no extra offset.
        struct cai_err {
            // X-axis PDGE
            double dxx = 0, dyx = 0, dzx = 0, exx = 0, eyx = 0, ezx = 0;
            // Z-axis PDGE
            double dxz = 0, dyz = 0, dzz = 0, exz = 0, eyz = 0, ezz = 0;
            // C1-axis PDGE (ezc1 absorbed into bt)
            double dxc1 = 0, dyc1 = 0, dzc1 = 0, exc1 = 0, eyc1 = 0;
            // C2-axis PDGE (ezc2 absorbed into bw)
            double dxc2 = 0, dyc2 = 0, dzc2 = 0, exc2 = 0, eyc2 = 0;
            // PIGE: squareness X-Z, parallelism C1, parallelism C2
            double eta = 0, azoc1 = 0, bzoc1 = 0, azoc2 = 0, bzoc2 = 0;
            // Abbe/Bryan measurement offsets L(j) (Eq.4-7), in mm
            double Lxz = 0, Lyz = 0, Lzz = 0;
            double Lxx = 0, Lyx = 0, Lzx = 0;
            double Lxc1 = 0, Lyc1 = 0, Lzc1 = 0;
            double Lxc2 = 0, Lyc2 = 0, Lzc2 = 0;
        };

        // ---- GA base functional point Pe (Cai Eq.2, identical to F.2)
        // -------------------
        auto Pe_ga = [&](cai_err const& e, double th) {
            kinematic_system3dp s;
            s.add_frame(static_frame3dp("bed"));
            s.add_frame(static_frame3dp("z_pige", O_3dp, vec3dp{0, e.eta, 0, 0}),
                        s.index_of("bed"));
            s.add_frame(static_frame3dp("z_axis", vec3dp{0, 0, Z, 1}),
                        s.index_of("z_pige"));
            s.add_frame(static_frame3dp("z_pdge", vec3dp{e.dxz, e.dyz, e.dzz, 1},
                                        vec3dp{e.exz, e.eyz, e.ezz, 0}),
                        s.index_of("z_axis"));
            s.add_frame(static_frame3dp("c1_pige", O_3dp, vec3dp{e.azoc1, e.bzoc1, 0, 0}),
                        s.index_of("z_pdge"));
            s.add_frame(static_frame3dp("c1_pdge", vec3dp{e.dxc1, e.dyc1, e.dzc1, 1},
                                        vec3dp{e.exc1, e.eyc1, 0, 0}),
                        s.index_of("c1_pige"));
            s.add_frame(static_frame3dp("c1_rot", O_3dp, vec3dp{0, 0, bt, 0}),
                        s.index_of("c1_pdge"));
            s.add_frame(static_frame3dp("tilt_x", O_3dp, vec3dp{ax, 0, 0, 0}),
                        s.index_of("c1_rot"));
            s.add_frame(static_frame3dp("tilt_y", O_3dp, vec3dp{0, ay, 0, 0}),
                        s.index_of("tilt_x"));
            s.add_frame(static_frame3dp("x_axis", vec3dp{X, 0, 0, 1}), s.index_of("bed"));
            s.add_frame(static_frame3dp("x_pdge", vec3dp{e.dxx, e.dyx, e.dzx, 1},
                                        vec3dp{e.exx, e.eyx, e.ezx, 0}),
                        s.index_of("x_axis"));
            s.add_frame(static_frame3dp("c2_pige", O_3dp, vec3dp{e.azoc2, e.bzoc2, 0, 0}),
                        s.index_of("x_pdge"));
            s.add_frame(static_frame3dp("c2_pdge", vec3dp{e.dxc2, e.dyc2, e.dzc2, 1},
                                        vec3dp{e.exc2, e.eyc2, 0, 0}),
                        s.index_of("c2_pige"));
            s.add_frame(static_frame3dp("c2_rot", O_3dp, vec3dp{0, 0, bw, 0}),
                        s.index_of("c2_pdge"));
            auto const Pt = vec3dp{Rt * std::cos(th), Rt * std::sin(th), 0.0, 1.0};
            return unitize(move3dp(Pt, s.get_pos_trafo("tilt_y", "c2_rot")));
        };

        // ---- GA Abbe/Bryan correction sum (Cai Eq.5-7)
        // --------------------------------- Each axis correction is the error twist's
        // velocity field at the offset point L(j), rotated into the wafer frame by
        // Rz(bw). twist = 2*log(M_e) is the full screw generator (log returns the
        // half-angle generator). The displacement is a direction (w=0); a pure Rz
        // rotation acts on its x,y,z exactly like Tc2 on the homogeneous EM.L point.
        auto const Rz_bw = motor_from_pose3dp(pose3dp{O_3dp, vec3dp{0, 0, bw, 0}});
        auto E_corr_ga = [&](vec3dp const& del, vec3dp const& eps_rot,
                             vec3dp const& Lpt) {
            auto const M_e =
                motor_from_pose3dp(pose3dp{vec3dp{del.x, del.y, del.z, 1.0}, eps_rot});
            auto const twist = 2.0 * log(M_e); // full screw generator Omega
            auto const d = rcmt(twist, Lpt);   // velocity field at the offset point
            return move3dp(vec3dp{d.x, d.y, d.z, 0.0}, Rz_bw); // rotate into wafer frame
        };
        auto E_ab_ga = [&](cai_err const& e, double th) {
            double const bc = Rt * std::cos(th + bw), bs = Rt * std::sin(th + bw);
            // measurement points L(j) (Eq.4-7): Pt's nominal position + offset L. The X-
            // and C2-axes carry the carriage offset (-X in x, +Z in z); the Z-, C1-axes
            // and the PIGE terms do not.
            vec3dp const sum =
                E_corr_ga({e.dxz, e.dyz, e.dzz, 0}, {e.exz, e.eyz, e.ezz, 0},
                          {bc + e.Lxz, bs + e.Lyz, e.Lzz, 1}) + // E(z)
                E_corr_ga({e.dxx, e.dyx, e.dzx, 0}, {e.exx, e.eyx, e.ezx, 0},
                          {bc - X + e.Lxx, bs + e.Lyx, Z + e.Lzx, 1}) + // E(x)
                E_corr_ga({e.dxc1, e.dyc1, e.dzc1, 0}, {e.exc1, e.eyc1, 0, 0},
                          {bc + e.Lxc1, bs + e.Lyc1, e.Lzc1, 1}) + // E(c1)
                E_corr_ga({e.dxc2, e.dyc2, e.dzc2, 0}, {e.exc2, e.eyc2, 0, 0},
                          {bc - X + e.Lxc2, bs + e.Lyc2, Z + e.Lzc2, 1}) +  // E(c2)
                E_corr_ga({0, 0, 0, 0}, {0, e.eta, 0, 0}, {bc, bs, Z, 1}) + // E(xoz)
                E_corr_ga({0, 0, 0, 0}, {e.azoc1, e.bzoc1, 0, 0},
                          {bc, bs, Z, 1}) + // E(zoc1)
                E_corr_ga({0, 0, 0, 0}, {e.azoc2, e.bzoc2, 0, 0},
                          {bc - X, bs, Z, 1}); // E(zoc2)
            return sum;
        };

        // ---- independent 4x4 HTM reference (Cai's own method)
        // --------------------------
        using mat4 = std::array<std::array<double, 4>, 4>;
        auto ident = []() {
            mat4 M{};
            for (int i = 0; i < 4; ++i)
                M[i][i] = 1.0;
            return M;
        };
        auto mul = [](mat4 const& A, mat4 const& B) {
            mat4 C{};
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 4; ++k)
                        s += A[i][k] * B[k][j];
                    C[i][j] = s;
                }
            return C;
        };
        auto transl = [&](double dx, double dy, double dz) {
            mat4 M = ident();
            M[0][3] = dx, M[1][3] = dy, M[2][3] = dz;
            return M;
        };
        auto Rz = [&](double a) {
            mat4 M = ident();
            M[0][0] = std::cos(a), M[0][1] = -std::sin(a), M[1][0] = std::sin(a),
            M[1][1] = std::cos(a);
            return M;
        };
        auto Rx = [&](double a) {
            mat4 M = ident();
            M[1][1] = std::cos(a), M[1][2] = -std::sin(a), M[2][1] = std::sin(a),
            M[2][2] = std::cos(a);
            return M;
        };
        auto Ry = [&](double a) {
            mat4 M = ident();
            M[0][0] = std::cos(a), M[0][2] = std::sin(a), M[2][0] = -std::sin(a),
            M[2][2] = std::cos(a);
            return M;
        };
        auto emat = [&](double dx, double dy, double dz, double ex, double ey,
                        double ez) {
            mat4 M = ident();
            M[0][1] = -ez, M[0][2] = ey, M[1][0] = ez, M[1][2] = -ex, M[2][0] = -ey,
            M[2][1] = ex;
            M[0][3] = dx, M[1][3] = dy, M[2][3] = dz;
            return M;
        };
        auto rigid_inv = [&](mat4 const& M) {
            mat4 I = ident();
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    I[i][j] = M[j][i];
            for (int i = 0; i < 3; ++i)
                I[i][3] = -(I[i][0] * M[0][3] + I[i][1] * M[1][3] + I[i][2] * M[2][3]);
            return I;
        };
        // base Pe (Eq.2) -- same chain as F.2's HTM reference
        auto Pe_htm = [&](cai_err const& e, double th) {
            mat4 const tilt = mul(Rx(ax), Ry(ay));
            mat4 K = mul(emat(0, 0, 0, 0, e.eta, 0),
                         mul(transl(0, 0, Z),
                             mul(emat(e.dxz, e.dyz, e.dzz, e.exz, e.eyz, e.ezz),
                                 mul(emat(0, 0, 0, e.azoc1, e.bzoc1, 0),
                                     mul(emat(e.dxc1, e.dyc1, e.dzc1, e.exc1, e.eyc1, 0),
                                         mul(Rz(bt), tilt))))));
            mat4 W = mul(
                transl(X, 0, 0),
                mul(emat(e.dxx, e.dyx, e.dzx, e.exx, e.eyx, e.ezx),
                    mul(emat(0, 0, 0, e.azoc2, e.bzoc2, 0),
                        mul(emat(e.dxc2, e.dyc2, e.dzc2, e.exc2, e.eyc2, 0), Rz(bw)))));
            mat4 const T = mul(rigid_inv(W), K);
            double const px = Rt * std::cos(th), py = Rt * std::sin(th);
            return std::array<double, 3>{T[0][0] * px + T[0][1] * py + T[0][3],
                                         T[1][0] * px + T[1][1] * py + T[1][3],
                                         T[2][0] * px + T[2][1] * py + T[2][3]};
        };
        // EM(j) (zero-diagonal differential matrix) applied to L(j), then Tc2 = Rz(bw)
        auto Tc2_EM_L = [&](double dx, double dy, double dz, double ex, double ey,
                            double ez, std::array<double, 3> const& L) {
            double const vx =
                -ez * L[1] + ey * L[2] + dx; // EM.L (linearised displacement)
            double const vy = ez * L[0] - ex * L[2] + dy;
            double const vz = -ey * L[0] + ex * L[1] + dz;
            return std::array<double, 3>{std::cos(bw) * vx - std::sin(bw) * vy,
                                         std::sin(bw) * vx + std::cos(bw) * vy, vz};
        };
        auto E_ab_htm = [&](cai_err const& e, double th) {
            double const bc = Rt * std::cos(th + bw), bs = Rt * std::sin(th + bw);
            std::array<double, 3> s{0, 0, 0};
            auto add = [&](std::array<double, 3> const& a) {
                s[0] += a[0], s[1] += a[1], s[2] += a[2];
            };
            add(Tc2_EM_L(e.dxz, e.dyz, e.dzz, e.exz, e.eyz, e.ezz,
                         {bc + e.Lxz, bs + e.Lyz, e.Lzz}));
            add(Tc2_EM_L(e.dxx, e.dyx, e.dzx, e.exx, e.eyx, e.ezx,
                         {bc - X + e.Lxx, bs + e.Lyx, Z + e.Lzx}));
            add(Tc2_EM_L(e.dxc1, e.dyc1, e.dzc1, e.exc1, e.eyc1, 0,
                         {bc + e.Lxc1, bs + e.Lyc1, e.Lzc1}));
            add(Tc2_EM_L(e.dxc2, e.dyc2, e.dzc2, e.exc2, e.eyc2, 0,
                         {bc - X + e.Lxc2, bs + e.Lyc2, Z + e.Lzc2}));
            add(Tc2_EM_L(0, 0, 0, 0, e.eta, 0, {bc, bs, Z}));             // E(xoz)
            add(Tc2_EM_L(0, 0, 0, e.azoc1, e.bzoc1, 0, {bc, bs, Z}));     // E(zoc1)
            add(Tc2_EM_L(0, 0, 0, e.azoc2, e.bzoc2, 0, {bc - X, bs, Z})); // E(zoc2)
            return s;
        };

        // representative full error set + Abbe/Bryan offsets
        cai_err full{.dxx = 2e-3,
                     .dyx = 1e-3,
                     .dzx = 1.5e-3,
                     .exx = 3 * as,
                     .eyx = 2 * as,
                     .ezx = 1 * as,
                     .dxz = -3.6e-3,
                     .dyz = -2.5e-3,
                     .dzz = 1e-3,
                     .exz = 2 * as,
                     .eyz = -1.5 * as,
                     .ezz = 1.5 * as,
                     .dxc1 = 1e-3,
                     .dyc1 = -1e-3,
                     .dzc1 = 2e-3,
                     .exc1 = 1 * as,
                     .eyc1 = 1.5 * as,
                     .dxc2 = -1e-3,
                     .dyc2 = 2e-3,
                     .dzc2 = -1.5e-3,
                     .exc2 = 2 * as,
                     .eyc2 = -1 * as,
                     .eta = 2 * as,
                     .azoc1 = 1.5 * as,
                     .bzoc1 = 1 * as,
                     .azoc2 = -1 * as,
                     .bzoc2 = 2 * as,
                     .Lxz = 10.0,
                     .Lyz = -8.0,
                     .Lzz = 15.0,
                     .Lxx = 12.0,
                     .Lyx = 5.0,
                     .Lzx = -10.0,
                     .Lxc1 = 8.0,
                     .Lyc1 = -6.0,
                     .Lzc1 = 20.0,
                     .Lxc2 = -7.0,
                     .Lyc2 = 9.0,
                     .Lzc2 = 12.0};
        cai_err const none{};

        int const nth = 12;

        // (1) GA Abbe/Bryan correction == independent HTM of Eqs.(5-7), over a th sweep
        double max_ab_gap = 0.0, ab_max = 0.0;
        for (int k = 0; k <= nth; ++k) {
            double const th = 2.0 * pi * double(k) / double(nth);
            auto const g = E_ab_ga(full, th);
            auto const h = E_ab_htm(full, th);
            max_ab_gap = std::max({max_ab_gap, std::abs(g.x - h[0]), std::abs(g.y - h[1]),
                                   std::abs(g.z - h[2])});
            ab_max = std::max(ab_max, std::sqrt(g.x * g.x + g.y * g.y + g.z * g.z));
            CHECK(g.x == doctest::Approx(h[0]).epsilon(1e-5));
            CHECK(g.y == doctest::Approx(h[1]).epsilon(1e-5));
            CHECK(g.z == doctest::Approx(h[2]).epsilon(1e-5));
        }
        fmt::println("(1) GA Abbe/Bryan sum vs HTM Eq.(5-7): max gap = {:.3e} mm, "
                     "|E_AB|_max = {:.4f} um",
                     max_ab_gap, 1e3 * ab_max);

        // (2) full optimized error E' = E(F.2) + sum_j E(j): fully-GA vs assembled HTM of
        // Eq.(8). E(F.2) = Pe - Pi (base volumetric error); E_AB the correction.
        double max_eprime_gap = 0.0, eprime_max = 0.0;
        for (int k = 0; k <= nth; ++k) {
            double const th = 2.0 * pi * double(k) / double(nth);
            auto const pe = Pe_ga(full, th);
            auto const pi_ = Pe_ga(none, th);
            auto const ab = E_ab_ga(full, th);
            double const gx = (pe.x - pi_.x) + ab.x, gy = (pe.y - pi_.y) + ab.y,
                         gz = (pe.z - pi_.z) + ab.z;

            auto const he = Pe_htm(full, th);
            auto const hi = Pe_htm(none, th);
            auto const hab = E_ab_htm(full, th);
            double const hx = (he[0] - hi[0]) + hab[0], hy = (he[1] - hi[1]) + hab[1],
                         hz = (he[2] - hi[2]) + hab[2];

            max_eprime_gap = std::max({max_eprime_gap, std::abs(gx - hx),
                                       std::abs(gy - hy), std::abs(gz - hz)});
            eprime_max = std::max(eprime_max, std::sqrt(gx * gx + gy * gy + gz * gz));
            CHECK(gx == doctest::Approx(hx).epsilon(1e-5));
            CHECK(gy == doctest::Approx(hy).epsilon(1e-5));
            CHECK(gz == doctest::Approx(hz).epsilon(1e-5));
            if (k % (nth / 4) == 0)
                fmt::println("th={:>5.2f}: E' = ({:>9.4f}, {:>9.4f}, {:>9.4f}) um", th,
                             1e3 * gx, 1e3 * gy, 1e3 * gz);
        }
        fmt::println("(2) full E' (GA base+corr) vs HTM Eq.(8): max gap = {:.3e} mm, "
                     "|E'|_max = {:.4f} um",
                     max_eprime_gap, 1e3 * eprime_max);

        // (3) isolated-offset spot check: only eyz (Z-axis pitch) + its z-offset Lzz. The
        // Abbe/Bryan E_AB then has the lone cross term eyz.Lzz in EM(z).L(z)'s x-row, so
        // E_AB_x == cbw * eyz * Lzz -- exactly the eyz.Lzz term of Eq.(8)'s E'_x
        // cbw-brace.
        {
            cai_err iso{};
            iso.eyz = 4 * as;
            iso.Lzz = 25.0;
            for (int k = 0; k <= nth; ++k) {
                double const th = 2.0 * pi * double(k) / double(nth);
                auto const ab = E_ab_ga(iso, th);
                double const term = std::cos(bw) * iso.eyz * iso.Lzz;
                CHECK(ab.x == doctest::Approx(term).epsilon(1e-3));
            }
            fmt::println("(3) isolated eyz/Lzz: E_AB_x == cbw*eyz*Lzz = {:.4f} um, "
                         "the eyz.Lzz term of Eq.(8) E'_x  [OK]",
                         1e3 * std::cos(bw) * iso.eyz * iso.Lzz);
        }

        // (4) PAPER-ERROR flag (the GA is the tie-break). Isolate the X/Z straightness
        // delta_xx, delta_xz in the BASE error E(F.2) (= Eq.3). Their brace enters E_x as
        // (delta_xz - delta_xx).cbw and E_y as -(delta_xz - delta_xx).sbw -- ONE brace,
        // because both come from rotating one displacement by one Rz(bw). The GA thus
        // confirms Eq.(3)'s (delta_xz - delta_xx) in BOTH rows and flags Eq.(8)'s printed
        // E'_y brace (delta_xx - delta_xz) as a sign typo.
        {
            cai_err iso{};
            iso.dxx = 5e-3, iso.dxz = 2e-3;
            double const brace = iso.dxz - iso.dxx; // Eq.(3) form
            for (int k = 0; k <= nth; ++k) {
                double const th = 2.0 * pi * double(k) / double(nth);
                auto const pe = Pe_ga(iso, th);
                auto const pi_ = Pe_ga(none, th);
                double const ex = pe.x - pi_.x, ey = pe.y - pi_.y;
                CHECK(ex == doctest::Approx(brace * std::cos(bw)).epsilon(2e-2));
                CHECK(ey == doctest::Approx(-brace * std::sin(bw)).epsilon(2e-2));
            }
            fmt::println("(4) GA confirms Eq.(3) brace (dxz-dxx) in BOTH E_x and E_y; "
                         "Eq.(8)'s printed E'_y (dxx-dxz) is a sign typo  [flagged]");
        }

        fmt::println("");
        fmt::println("Abbe/Bryan: GA error-twist velocity field at the offset point "
                     "reproduces Cai Eq.(5-8) (Phase F.3)");
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // F.4 -- wafer surface topography + TTV from the volumetric error (Cai Eq.9-18). The
    // surface-formation half: a wheel-rim grain (functional point Pt) traces a path in
    // the rotating wafer frame (Eq.9), the volumetric error perturbs it (Eq.12 xp1 = xp0
    // + E'), and the ground surface is the per-radius MIN-envelope of the grain height
    // (Eq.13 Z = zp1-min | (X,Y)) -- rotationally symmetric, so it collapses to a 1D
    // radial profile. Wafer thickness t = zp1-min + t0 - zc-min (Eq.14), TTV = max(t) -
    // min(t) (Eq.16).
    //
    // GA reading (the showcase): the paper builds the ideal trajectory (Eq.9), then a
    // page-long LINEARISED error map E' (Eq.8), then ADDS them (Eq.12). In GA we inject
    // the error as a MOTOR into the kinematic chain and trace the ACTUAL rim point --
    // Eq.12 done exactly, no linearisation, no E' assembly. The chain POSITION of the
    // error is the physics: a Z-axis angular error (eps_xz, eps_yz) sits BEFORE the wheel
    // spin (a fixed, non-rotating wheel-axis tilt) -> it breaks the radial symmetry -> a
    // CONE; the nominal adjustment tilt (ax, ay) sits AFTER the spin (rotates WITH the
    // wheel) -> it averages out under the min-envelope -> FLAT. This is exactly why Cai
    // Fig.8(a) is flat despite ax = ay = 0.01 deg, and it falls out of where the frames
    // sit in the tree.
    //
    // GATES: (a) eps_xz = eps_yz = 0 -> flat (Fig.8a); (b) eps_xz = 5" -> a single-signed
    // sin-cone of amplitude Rt*eps_xz (the tilt lever arm); (c) eps_yz = 5" -> a
    // symmetric cos-cone spanning center-to-rim, TTV ~ 2*Rt*eps_yz; (d) both -> warped,
    // TTV ~ the larger; the TTV magnitudes (1.8 / 3.5 um) bracket Fig.8's colorbars (1.4
    // / 3.4 um). Plus: the exact z-displacement == Cai Eq.8's E'_z term Rt(eps_xz s(bt) -
    // eps_yz c(bt)) (ties F.4 to F.3); the Eq.14 thickness with a supplied chuck height
    // field zc(r); and the ZHOU<->CAI cross-check (the cone obeys the same lever-arm*tilt
    // law as the Zhou wafer-tilt cone above).
    //
    // PAPER-ERROR WATCH: Fig.8 labels (b) eps_xz with the LARGER (3.4 um) colorbar and
    // (c) eps_yz with the smaller (1.4 um), but the GA gives eps_xz the SMALLER
    // (sin-cone, ~1.8 um) and eps_yz the LARGER (cos-cone, ~3.5 um) -- i.e. the
    // eps_xz<->eps_yz roles are swapped vs Fig.8's print. This is consistent with the
    // malformed RMt in Cai Eq.10 (its second column [1,0,0]^T is not a rotation); the GA
    // uses a proper rotation, so it is the tie-break. Flag, do not silently reconcile.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE(
        "pga3dp: Cai wafer topography + TTV from volumetric error Eq.9-18 (Phase F.4)")
    {
        fmt::println("");
        fmt::println("pga3dp: Cai wafer topography + TTV from volumetric error Eq.9-18 "
                     "(Phase F.4)");
        fmt::println("");

        // Table-2 parameters
        double const Rt = 72.5, Rw = 150.0; // wheel / wafer radius [mm]
        double const X = Rt; // wheel-wafer centre offset: rim reaches the wafer centre
        double const Z = 50.0;
        double const ax = 0.01 * pi / 180,
                     ay = 0.01 * pi / 180;     // nominal adjustment tilt
        double const as = pi / 180.0 / 3600.0; // 1 arcsec [rad]
        double const wt = rpm2radps(2400.0), ww = rpm2radps(80.0); // wheel / wafer spin

        int const NR = 75; // radial Z-map bins over [0, Rw]
        double const dr = Rw / NR;

        // ---- min-envelope radial topography (Eq.13): trace the wheel-rim grain through
        // the Cai machine with a Z-axis angular error (exz,eyz) injected BEFORE the wheel
        // spin; take the lowest grain height per wafer-radius bin over one wafer
        // revolution.
        auto topo = [&](double exz, double eyz) {
            std::vector<double> zmin(NR, 1e9);
            std::vector<int> hit(NR, 0);
            kinematic_system3dp s;
            s.add_frame(static_frame3dp("bed"));
            s.add_frame(static_frame3dp("z_axis", vec3dp{0, 0, Z, 1}), s.index_of("bed"));
            s.add_frame(static_frame3dp("z_err", O_3dp, vec3dp{exz, eyz, 0, 0}),
                        s.index_of("z_axis")); // fixed wheel-axis tilt (before the spin)
            s.add_frame(static_frame3dp("c1_rot", O_3dp, vec3dp{0, 0, 0, 0}),
                        s.index_of("z_err"));
            s.add_frame(static_frame3dp("tilt_x", O_3dp, vec3dp{ax, 0, 0, 0}),
                        s.index_of("c1_rot")); // nominal tilt (after the spin)
            s.add_frame(static_frame3dp("tilt_y", O_3dp, vec3dp{0, ay, 0, 0}),
                        s.index_of("tilt_x"));
            s.add_frame(static_frame3dp("x_axis", vec3dp{X, 0, 0, 1}), s.index_of("bed"));
            s.add_frame(static_frame3dp("c2_rot", O_3dp, vec3dp{0, 0, 0, 0}),
                        s.index_of("x_axis"));
            size_t const c1 = s.index_of("c1_rot"), c2 = s.index_of("c2_rot");
            double const Tw = 2.0 * pi / ww; // one wafer revolution [s]
            int const N = 20000;
            auto const Pt = vec3dp{Rt, 0, 0, 1};
            for (int k = 0; k <= N; ++k) {
                double const t = Tw * double(k) / double(N);
                s.set_pose(c1, pose3dp{O_3dp, vec3dp{0, 0, wt * t, 0}});
                s.set_pose(c2, pose3dp{O_3dp, vec3dp{0, 0, ww * t, 0}});
                auto const g = unitize(move3dp(Pt, s.get_pos_trafo("tilt_y", "c2_rot")));
                double const r = std::sqrt(g.x * g.x + g.y * g.y); // wafer radius (e1-e2)
                if (r > Rw) continue;
                int const b = std::min(NR - 1, int(r / dr));
                // *** GROUND-SURFACE HEIGHT PROFILE Z(r) (Cai Eq.13) ***: the wafer upper
                // surface at radius r is the LOWEST grain pass there (the deepest cut
                // wins, Z = zp1-min). zmin[b] accumulates that per-radius minimum over
                // the whole grain trajectory -> the raw height profile that the cone is
                // read from.
                if (g.z < zmin[b]) zmin[b] = g.z;
                ++hit[b];
            }
            return std::make_pair(zmin, hit);
        };

        // deviation of an error topography vs the ideal one -> the volumetric-error
        // contribution (Fig.8), with the feed / nominal-tilt baseline cancelled. Returns
        // {lo, hi, center, ttv} in mm over the hit radii.
        auto const ideal = topo(0.0, 0.0);
        auto dev = [&](std::pair<std::vector<double>, std::vector<int>> const& topo_e) {
            double lo = 1e9, hi = -1e9, ctr = 0.0;
            for (int i = 0; i < NR; ++i) {
                if (topo_e.second[i] == 0 || ideal.second[i] == 0) continue;
                // *** THE WAFER HEIGHT PROFILE (the Fig.8 cone) ***: the deviation of the
                // error topography from the IDEAL (no-error, flat) min-envelope at radius
                // r. Subtracting the ideal cancels the feed / nominal-tilt baseline, so d
                // is purely the volumetric-error contribution -- the
                // convex/concave/warped cone the gates (b)/(c)/(d) below report.
                double const d = topo_e.first[i] - ideal.first[i];
                lo = std::min(lo, d);
                hi = std::max(hi, d);
                if (i == 0) ctr = d; // innermost bin ~ wafer centre
            }
            return std::array<double, 4>{lo, hi, ctr, hi - lo};
        };

        // (a) no Z-axis angular error -> FLAT (Fig.8a). The nominal ax,ay tilt averages
        // out under the min-envelope because it rotates with the wheel spin.
        {
            double spread = 0.0, lo = 1e9, hi = -1e9;
            for (int i = 0; i < NR; ++i) {
                if (ideal.second[i] == 0) continue;
                lo = std::min(lo, ideal.first[i]);
                hi = std::max(hi, ideal.first[i]);
            }
            spread = hi - lo;
            CHECK(spread < 1e-5); // flat to < 0.01 um despite ax = ay = 0.01 deg
            fmt::println("(a) eps_xz=eps_yz=0: ideal surface spread = {:.4f} um -> FLAT "
                         "(Fig.8a)",
                         1e3 * spread);
        }

        double const lever =
            Rt * 5.0 * as; // 72.5 mm * 5 arcsec = 1.7574 um tilt lever arm

        // (b) eps_xz = 5" -> single-signed sin-cone of amplitude Rt*eps_xz
        auto const db = dev(topo(5.0 * as, 0.0));
        {
            CHECK(db[1] < 0.05e-3); // single-signed (no positive lobe)
            CHECK(std::abs(db[0]) ==
                  doctest::Approx(lever).epsilon(0.03)); // amplitude == Rt*eps_xz
            CHECK(std::abs(db[2]) < 0.1e-3);             // ~0 at the wafer centre
            fmt::println("(b) eps_xz=5\": dev in [{:+.3f},{:+.3f}] um, TTV={:.3f} um  "
                         "-> sin-cone (amp=Rt*eps_xz={:.3f})",
                         1e3 * db[0], 1e3 * db[1], 1e3 * db[3], 1e3 * lever);
        }

        // (c) eps_yz = 5" -> symmetric cos-cone center-to-rim, TTV ~ 2*Rt*eps_yz
        auto const dc = dev(topo(0.0, 5.0 * as));
        {
            CHECK(dc[0] < 0.0); // spans both signs
            CHECK(dc[1] > 0.0);
            CHECK(dc[3] == doctest::Approx(2.0 * lever).epsilon(0.05)); // TTV ~ 2*lever
            fmt::println("(c) eps_yz=5\": dev in [{:+.3f},{:+.3f}] um, TTV={:.3f} um  "
                         "-> cos-cone (TTV~2*Rt*eps_yz={:.3f})",
                         1e3 * dc[0], 1e3 * dc[1], 1e3 * dc[3], 2e3 * lever);
        }

        // (d) both -> warped, TTV ~ the larger (cos + sin superpose)
        auto const dd = dev(topo(5.0 * as, 5.0 * as));
        {
            CHECK(dd[3] > dc[3] * 0.9);   // at least as large as the eps_yz cone
            CHECK(dd[3] < dc[3] + db[3]); // sub-additive (not a simple sum)
            fmt::println("(d) both 5\": dev in [{:+.3f},{:+.3f}] um, TTV={:.3f} um  "
                         "-> warped cone (Fig.8d)",
                         1e3 * dd[0], 1e3 * dd[1], 1e3 * dd[3]);
        }

        // (e) tie to F.3: the EXACT grain z-displacement at the functional point equals
        // Cai Eq.8's E'_z term Rt(eps_xz*s(bt) - eps_yz*c(bt)) over a wheel-phase sweep
        // -- the same error-twist z-velocity field F.3 validated, now driving the
        // surface.
        {
            double const exz = 5.0 * as, eyz = 3.0 * as;
            auto zsingle = [&](double ex, double ey, double bt) {
                kinematic_system3dp s;
                s.add_frame(static_frame3dp("bed"));
                s.add_frame(static_frame3dp("z_axis", vec3dp{0, 0, Z, 1}),
                            s.index_of("bed"));
                s.add_frame(static_frame3dp("z_err", O_3dp, vec3dp{ex, ey, 0, 0}),
                            s.index_of("z_axis"));
                s.add_frame(static_frame3dp("c1_rot", O_3dp, vec3dp{0, 0, bt, 0}),
                            s.index_of("z_err"));
                s.add_frame(static_frame3dp("tilt_x", O_3dp, vec3dp{ax, 0, 0, 0}),
                            s.index_of("c1_rot"));
                s.add_frame(static_frame3dp("tilt_y", O_3dp, vec3dp{0, ay, 0, 0}),
                            s.index_of("tilt_x"));
                s.add_frame(static_frame3dp("x_axis", vec3dp{X, 0, 0, 1}),
                            s.index_of("bed"));
                s.add_frame(static_frame3dp("c2_rot", O_3dp, vec3dp{0, 0, 0.3, 0}),
                            s.index_of("x_axis"));
                return double(unitize(move3dp(vec3dp{Rt, 0, 0, 1},
                                              s.get_pos_trafo("tilt_y", "c2_rot")))
                                  .z);
            };
            double max_gap = 0.0;
            for (int k = 0; k < 8; ++k) {
                double const bt = 2.0 * pi * double(k) / 8.0;
                double const dz = zsingle(exz, eyz, bt) - zsingle(0, 0, bt);
                double const eq8 = Rt * (exz * std::sin(bt) - eyz * std::cos(bt));
                max_gap = std::max(max_gap, std::abs(dz - eq8));
                CHECK(dz == doctest::Approx(eq8).epsilon(1e-4));
            }
            fmt::println("(e) exact dZ == Eq.8 E'_z = Rt(eps_xz s(bt)-eps_yz c(bt)) to "
                         "{:.2e} mm  [ties F.4->F.3]",
                         max_gap);
        }

        // (f) Eq.14 thickness + Eq.16 TTV with a supplied chuck height field zc(r) (the
        // DESIGN PRINCIPLE: chuck = rigid carrier + residual scalar field; here the
        // carrier is identity and zc(r) is the 1D radial dressing profile). Flat chuck ->
        // TTV is the topography spread; a chuck cone matching the wafer cone cancels it
        // (the wafer lower surface conforms to the chuck, t = zp1-min + t0 - zc).
        {
            auto const topo_e = topo(0.0, 5.0 * as);
            double const t0 = 0.8; // commanded wafer thickness [mm]
            auto ttv = [&](auto const& zc) {
                double lo = 1e9, hi = -1e9;
                for (int i = 0; i < NR; ++i) {
                    if (topo_e.second[i] == 0) continue;
                    double const r = (i + 0.5) * dr;
                    // *** WAFER THICKNESS vs THE CHUCK (Cai Eq.14) ***: t =
                    // (upper-surface height profile zp1-min) + (commanded thickness t0) -
                    // (chuck profile zc(r)). zc IS the chuck height field; for a FLAT
                    // chuck zc==0 (below), so the thickness profile == the height profile
                    // -- the wafer's deviation referenced to the flat chuck.
                    double const t = topo_e.first[i] + t0 - zc(r); // Eq.14
                    lo = std::min(lo, t);
                    hi = std::max(hi, t);
                }
                return hi - lo; // Eq.16 TTV = max(t) - min(t)
            };
            // flat chuck (zc == 0): the lower surface is planar, so the thickness profile
            // reduces to the upper-surface height profile -> TTV == the topography spread
            double const ttv_flat = ttv([](double) { return 0.0; });
            CHECK(ttv_flat == doctest::Approx(dc[3]).epsilon(0.02)); // == topography TTV
            // a chuck dressed to the SAME error cone (zc = zp1-min) makes t = t0 + const
            double const ttv_conform = ttv(
                [&](double r) { return topo_e.first[std::min(NR - 1, int(r / dr))]; });
            CHECK(ttv_conform < 1e-9); // the chuck height field cancels the wafer cone
            fmt::println("(f) Eq.14/16 TTV: flat chuck = {:.3f} um (= topography); chuck "
                         "dressed to the cone -> TTV = {:.3e} um (cancels)",
                         1e3 * ttv_flat, 1e3 * ttv_conform);
        }

        // (g) ZHOU<->CAI cross-check: the Cai eps_xz/eps_yz cone obeys the SAME
        // lever-arm * tilt law as the Zhou wafer-tilt cone (R*tan(alpha), validated in
        // the "wafer-tilt flatness profile" test above) -- a fixed axis tilt of magnitude
        // eta gives a radial cone of amplitude (lever radius)*eta. Cai's lever is the
        // wheel radius Rt; Zhou's is the wafer radius Rw. Same mechanism, two independent
        // models.
        {
            CHECK(std::abs(db[0]) == doctest::Approx(Rt * 5.0 * as).epsilon(0.03));
            fmt::println(
                "(g) Zhou<->Cai: Cai cone amp = Rt*eps = {:.3f} um obeys the same "
                "lever*tilt cone law as Zhou R*tan(alpha)",
                1e3 * std::abs(db[0]));
        }

        fmt::println("");
        fmt::println("Topography: GA trajectory + error-motor injection + min-envelope "
                     "reproduces Cai Fig.8 / TTV (Phase F.4)");
        fmt::println("");
    }
}
