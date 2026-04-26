// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <chrono>
#include <iostream>
#include <tuple>
#include <vector>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// Euclidean geometric algebra of three-dimensional space: ega3d, pga3dp comparison
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: comparison tests")
{

    TEST_CASE("MVec3dp: rotated motors (for robotics applications)")
    {
        fmt::println("MVec3dp: rotated motors (for robotics applications)");

        // get motor that rotates around z_axis (e3 direction) by +45°
        auto M1 = get_motor(z_axis_3dp, deg2rad(45));
        CHECK(move3dp(e1_3dp, M1) == vec3dp{1, 1, 0, 0} / sqrt(2.));

        // now modify the rotor to rotate around x_axis (e1 direction) by +90° in addition
        auto M2 = get_motor(x_axis_3dp, deg2rad(90));
        auto M1_rotated = move3dp(M1, M2);
        CHECK(move3dp(e1_3dp, M1_rotated) == vec3dp{1, 0, 1, 0} / sqrt(2.));
        fmt::println("");
    }

    TEST_CASE("pga3dp: speed of rotation (via bivector)")
    {

        fmt::println("");
        fmt::println("pga3dp: speed of rotation (via bivector)");
        fmt::println("");


        // EGA3D:
        //        for omega = const. (a ega3d vector representing a rotational axis)
        //
        //        ds/dt = v = omega x r    (with x symbolizing the cross product)
        //        dv/dt = a = omega x v
        //
        //        with Omega = cmpl(omega)  (a bivector)
        //        ds/dt = v = Omega >> r = cmt(r, Omega)
        //        dv/dt = a = Omega >> v = cmt(v, Omega)
        //
        // PGA3DP:
        //        for omega = const. (a ega3dp vector representing a rotational axis)
        //
        //        with Omega = att(r_bulk_dual(omega)) a bivector)
        //        ds/dt = v = Omega >> r = cmt(r, Omega)
        //        dv/dt = a = Omega >> v = cmt(v, Omega)

        // rotating plane defined by point P and bivector B

        // defined in EGA3D
        auto omega_ega = cross(vec3d{1, 0.5, 0}, vec3d{0, 1, 0}); // normal of rot.-plane
                                                                  // = rotation axis

        omega_ega /= nrm(omega_ega); // normalize to length 1

        auto r_ega = vec3d{1, 1, 0} / nrm(vec3d{1, 1, 0}); // relative vector to point
                                                           // (direction only, length 1)

        // defined in PGA3DP
        auto B = wdg(vec3dp{1, 0.5, 0, 0}, vec3dp{0, 1, 0, 0}); // bivector in rot.-plane

        // bulk_normalize B to prepare for multiplication with turning rate, if required
        // (otherwise speed would be scaled with bulk_nrm(B))
        B /= to_val(bulk_nrm(B));


        double omega = rpm2radps(300); // angular speed in rad/s

        omega_ega *= omega;     // angular velocity [rad/s] (encoded in axis of rotation)
        auto Omega = omega * B; // bivector that encodes turing rate and plane of rotation


        // speed v of rotation at point r relative to center of rotation
        auto r = vec3dp{1, 1, 0, 0};
        r /= to_val(bulk_nrm(r)); // normalize r to get direction vector of length 1

        auto v1_ega = cross(omega_ega, r_ega); // speed calculated in EGA3D
        auto Omega_ega = cmpl(omega_ega);      // express omega as bivector
        CHECK(v1_ega == Omega_ega >> r_ega);
        CHECK(v1_ega == cmt(r_ega, Omega_ega));

        auto v2 = Omega >> r;  // right contraction turns in direction of B (or Omega)
        auto a2 = Omega >> v2; // acceleration vector (not used here)

        // convert v1_ega to pga3dp type for comparison
        auto v1 = vec3dp{v1_ega.x, v1_ega.y, v1_ega.z, 0};

        CHECK(v1 == v2);
        CHECK(v2 == cmt(r, Omega));  // right contraction Omega >> r and cmt(r, Omega)
                                     // are equivalent
        CHECK(a2 == cmt(v2, Omega)); // right contraction Omega >> v and cmt(v, Omega)
        // are equivalent

        auto omega3dp = vec3dp{omega_ega.x, omega_ega.y, omega_ega.z, 1};
        CHECK(Omega == att(r_bulk_dual(omega3dp)));
        CHECK(Omega == att(r_cmpl(omega3dp)));

        fmt::println("B = {}", B);
        fmt::println("bulk(B)   = {}, bulk_nrm(B)   = {}", bulk(B), bulk_nrm(B));
        fmt::println("weight(B) = {}, weight_nrm(B) = {}", weight(B), weight_nrm(B));
        fmt::println("");
        fmt::println("r = {}", r);
        fmt::println("omega [rad/s]     = {}", omega);
        fmt::println("omega_ega         = {}", omega_ega);
        fmt::println("Omega = omega * B = {}", Omega);
        fmt::println("    att(r_bulk_dual(omega3dp))) = {}", att(r_bulk_dual(omega3dp)));
        fmt::println("    att(r_cmpl(omega3dp)) = {}", att(r_cmpl(omega3dp)));
        fmt::println("");
        fmt::println("v1_ega = omega * (omega_ega x r_ega) = {}", v1_ega);
        fmt::println("v1     = omega * (v1_ega, 1.0)       = {}, bulk_nrm(v1) = {}", v1,
                     to_val(bulk_nrm(v1)));
        fmt::println("v2 = Omega >> r                      = {}, bulk_nrm(v2) = {}", v2,
                     to_val(bulk_nrm(v2)));
        fmt::println("a2 = Omega >> v2                     = {}, bulk_nrm(a2) = {}", a2,
                     to_val(bulk_nrm(a2)));

        fmt::println("");

        fmt::println("");
    }

    TEST_CASE("pga3dp: application tests - force and moment (I/II)")
    {
        fmt::println("pga3dp: application tests - force and moment (I/II)");


        // PGA implementation
        auto P = vec3dp{2.25, 1, 0, 1};
        auto f = vec3dp{0, 0.75, 0, 0};
        auto F = wdg(P, f); // force line resulting from f acting at P (=forque)

        fmt::println("P = {}", P);
        fmt::println("f = {}", f);
        fmt::println("F = P^f = {}, att(f) = {}, M_O = bulk(F) = {}", F, att(F), bulk(F));
        fmt::println("");
        fmt::println("O_3dp^f = {}, att(O_3dp^f) = {} (=force)", wdg(O_3dp, f),
                     att(wdg(O_3dp, f)));
        fmt::println("p^f = {}, bulk(p^f) = {} (=torque)", wdg(P - O_3dp, f),
                     bulk(wdg(P - O_3dp, f)));
        fmt::println("");
        CHECK(sup(F) == vec3dp{2.25, 0, 0, 1});
        CHECK(wdg(P, f) == wdg(O_3dp, f) + wdg(P - O_3dp, f));

        auto R1 = vec3dp{1.5, 2, 0, 1};
        auto R2 = vec3dp{3.0, 2, 0, 1};
        auto R3 = vec3dp{3.0, 0, 0, 1};
        auto R4 = vec3dp{1.5, 0, 0, 1};
        auto F_R1 = wdg(P - R1, f);
        auto F_R2 = wdg(P - R2, f);
        auto F_R3 = wdg(P - R3, f);
        auto F_R4 = wdg(P - R4, f);

        fmt::println("R1 = {}", R1);
        fmt::println("R2 = {}", R2);
        fmt::println("R3 = {}", R3);
        fmt::println("R4 = {}", R4);
        fmt::println("");
        fmt::println("F_R1 = (P - R1)^f = {}, att(f) = {}, M_R1 = bulk(F_R1) = {}", F_R1,
                     att(F_R1), bulk(F_R1));
        fmt::println("F_R2 = (P - R2)^f = {}, att(f) = {}, M_R2 = bulk(F_R2) = {}", F_R2,
                     att(F_R2), bulk(F_R2));
        fmt::println("F_R3 = (P - R3)^f = {}, att(f) = {}, M_R3 = bulk(F_R3) = {}", F_R3,
                     att(F_R3), bulk(F_R3));
        fmt::println("F_R4 = (P - R4)^f = {}, att(f) = {}, M_R4 = bulk(F_R4) = {}", F_R4,
                     att(F_R4), bulk(F_R4));
        fmt::println("");

        CHECK(F_R1 == F - wdg(R1, f));
        CHECK(F_R2 == F - wdg(R2, f));
        CHECK(F_R3 == F - wdg(R3, f));
        CHECK(F_R4 == F - wdg(R4, f));

        fmt::println("");
    }

    TEST_CASE("pga3dp: application tests - force and moment (II/II)")
    {
        fmt::println("pga3dp: application tests - force and moment (II/II)");

        // EGA implementation
        auto p_ega = vec3d{1, 0, 0};
        auto f_ega = vec3d{0, 0, -1}; // yields a positive moment in a right-handed system

        auto p_cross_f_ega = cross(p_ega, f_ega);
        auto p_wdg_f_ega = wdg(p_ega, f_ega);

        fmt::println("p_ega = {}", p_ega);
        fmt::println("f_ega = {}", f_ega);
        fmt::println("pxf = {}, nrm(pxf) = {}", p_cross_f_ega, nrm(p_cross_f_ega));
        fmt::println("p^f = {}, nrm(p^f) = {}", p_wdg_f_ega, nrm(p_wdg_f_ega));
        fmt::println("");

        CHECK(nrm(p_cross_f_ega) == 1.0);

        // PGA implementation
        auto P = vec3dp{1, 0, 0, 1};
        auto f = vec3dp{0, 0, 1, 0}; // yields a positive moment in a right-handed system
        auto F = wdg(P, f);

        fmt::println("P = {}", P);
        fmt::println("f = {}", f);
        fmt::println("F = P^f = {},  att(f) = {}, bulk(F) = {}", F, att(F), bulk(F));
        fmt::println("");
        fmt::println("O_3dp^f = {}, att(O_3dp^f) = {} (=force)", wdg(O_3dp, f),
                     att(wdg(O_3dp, f)));
        fmt::println("p^f = {}, bulk(p^f) = {} (=torque)", wdg(P - O_3dp, f),
                     bulk(wdg(P - O_3dp, f)));
        fmt::println("");

        CHECK(sup(F) == P);

        auto R1 = vec3dp{0.5, 0, -1, 1};
        auto R2 = vec3dp{1.5, 0, -1, 1};

        auto F_R1 = wdg(P - R1, f);
        auto F_R2 = wdg(P - R2, f);
        auto d_perp1 = ortho_proj3dp(R1, F) - R1;
        auto d_perp2 = ortho_proj3dp(R2, F) - R2;

        fmt::println("R1 = {}", R1);
        fmt::println("R2 = {}", R2);
        fmt::println("F_R1 = wdg(P - R1, f) = {}, att(F_R1) = {}, bulk(F_R1) = {}", F_R1,
                     att(F_R1), bulk(F_R1));
        fmt::println("F_R2 = wdg(P - R2, f) = {}, att(F_R1) = {}, bulk(F_R2) = {}", F_R2,
                     att(F_R2), bulk(F_R2));
        fmt::println("");
        fmt::println("ortho_proj3dp(R1, F) = {}", ortho_proj3dp(R1, F));
        fmt::println("ortho_proj3dp(R2, F) = {}", ortho_proj3dp(R2, F));
        fmt::println("d_perp1 = {}", d_perp1);
        fmt::println("d_perp2 = {}", d_perp2);
        fmt::println("d_perp1 ^ f = {}", wdg(d_perp1, f));
        fmt::println("d_perp2 ^ f = {}", wdg(d_perp2, f));

        CHECK(wdg(P - R1, f) == wdg(d_perp1, f));
        CHECK(wdg(P - R2, f) == wdg(d_perp2, f));

        fmt::println("");
    }


    TEST_CASE("pga3dp: pure rotation (manually, ega3d, pga3dp)")
    {

        fmt::println("");
        fmt::println("pga3dp: pure rotation (manually, ega3d, pga3dp)");

        constexpr size_t tsteps = 12;
        constexpr double omega = Hz2radps(1);

        auto constexpr t0 = 0.0;
        auto constexpr tmax = 1.0;
        auto constexpr dt = (tmax - t0) / tsteps;

        // assume constant element and time-dependent transformation in world frame
        // bivector rate is constant, i.e. constant angular velocity

        fmt::println("");
        fmt::println("manual description:");

        for (size_t nt = 0; nt <= tsteps; nt++) {
            auto const t = t0 + nt * dt;
            auto const pm = vec3d{std::cos(omega * t), 0, -std::sin(omega * t)};
            // vm = d(pm)/dt
            auto const vm =
                vec3d{-omega * std::sin(omega * t), 0, -omega * std::cos(omega * t)};
            fmt::println("t = {:>-5.2f},  pm = {:>-7.4f},  vm = {:>-7.4f}", t, pm, vm);
        }

        fmt::println("");
        fmt::println("ega description:");

        auto constexpr p0e = vec3d{1, 0, 0};      // initial vector
        auto constexpr Be_rot = bivec3d{0, 1, 0}; // rotate in plane e31
        fmt::println("Be_rot = {:>-7.4f}", Be_rot);


        for (size_t nt = 0; nt <= tsteps; nt++) {
            auto const t = t0 + nt * dt;
            auto const Me_rot = get_rotor(Be_rot, omega * t); // time dependent rot. angle
            auto pe = rotate(p0e, Me_rot);
            auto ve = cmt(pe, Be_rot * omega);
            fmt::println("t = {:>-5.2f},  pe = {:>-7.4f},  "
                         "ve = {:>-7.4f}",
                         t, pe, ve);
            // fmt::println("t = {:>-5.2f},  pe = {:>-7.4f},  "
            //              "ve = {:>-7.4f},  Me_rot = {:>-7.4f}",
            //              t, pe, ve, Me_rot);
        }

        fmt::println("");
        fmt::println("pga description:");

        auto constexpr P0 = vec3dp{1, 0, 0, 1};            // initial vector
        auto constexpr B_rot = bivec3dp{0, 1, 0, 0, 0, 0}; // rotate around e2-axis
        fmt::println("B_rot = {:>-7.4f}", B_rot);

        for (size_t nt = 0; nt <= tsteps; nt++) {
            auto const t = t0 + nt * dt;
            auto const M_rot = get_motor(B_rot, omega * t); // time dependent rot. angle
            auto const P = move3dp(P0, M_rot);
            // auto const v = -rcmt(P, B_rot * omega);
            auto const v = rcmt(B_rot * omega, P); // get rid of minus sign
            fmt::println("t = {:>-5.2f},  P = {:>-7.4f},  "
                         "v = {:>-7.4f}",
                         t, P, v);
            // fmt::println("t = {:>-5.2f},  P = {:>-7.4f},  "
            //              "v = {:>-7.4f},  M_rot = {:>-7.4f}",
            //              t, P, v, M_rot);
        }

        // fmt::println("");
        // for (size_t nt = 0; nt <= tsteps; nt++) {
        //     auto t = t0 + nt * dt;
        //     auto M_rot = get_motor(B_rot, omega * t);
        //     fmt::println("t = {:>-5.2f}, M_rot = {:>-7.4f},  "
        //                  "rrev(M_rot) = {:>-7.4f}",
        //                  t, M_rot, rrev(M_rot));
        //     CHECK(gr4(rgpr(rrev(M_rot), M_rot)) == pscalar3dp{1.0});
        //     CHECK(gr4(rgpr(M_rot, rrev(M_rot))) == pscalar3dp{1.0});
        // }


        fmt::println("");
    }


    TEST_CASE("pga3dp: rotation and translation (manually, ega3d, pga3dp)")
    {

        fmt::println("");
        fmt::println("pga3dp: rotation and translation (manually, ega3d, pga3dp)");

        constexpr size_t tsteps = 12;
        constexpr double omega = Hz2radps(1);

        auto constexpr t0 = 0.0;
        auto constexpr tmax = 1.0;
        auto constexpr dt = (tmax - t0) / tsteps;

        // assume constant element and time-dependent transformation in world frame
        // bivector rate is constant, i.e. constant angular velocity

        fmt::println("");
        fmt::println("manual description:");

        auto const v0e = vec3d{1., 1., 0.} / std::sqrt(2.);

        for (size_t nt = 0; nt <= tsteps; nt++) {
            auto const t = t0 + nt * dt;
            auto const pm = vec3d{std::cos(omega * t), 0, -std::sin(omega * t)} + v0e * t;
            // vm = d(pm)/dt
            auto const vm =
                vec3d{-omega * std::sin(omega * t), 0, -omega * std::cos(omega * t)} +
                v0e;
            fmt::println("t = {:>-5.2f},  pm = {:>-7.4f},  vm = {:>-7.4f}", t, pm, vm);
        }

        fmt::println("");
        fmt::println("ega description:");

        auto constexpr p0e = vec3d{1, 0, 0};      // initial vector
        auto constexpr Be_rot = bivec3d{0, 1, 0}; // rotate in plane e31
        fmt::println("Be_rot = {:>-7.4f}", Be_rot);


        for (size_t nt = 0; nt <= tsteps; nt++) {
            auto const t = t0 + nt * dt;

            // compute in body system
            auto const Me_rot = get_rotor(Be_rot, omega * t); // time dependent rotation
            auto pe = rotate(p0e, Me_rot);
            auto ve = cmt(pe, Be_rot * omega);

            // compute in world system
            pe += v0e * t; // add displacement of system
            ve += v0e;     // add velocity of system
            fmt::println("t = {:>-5.2f},  pe = {:>-7.4f},  "
                         "ve = {:>-7.4f}",
                         t, pe, ve);
        }

        fmt::println("");
        fmt::println("pga description:");

        auto const v0 = vec3dp{1., 1., 0., 0.} / std::sqrt(2.);

        auto constexpr P0 = vec3dp{1, 0, 0, 1};            // initial vector
        auto constexpr B_rot = bivec3dp{0, 1, 0, 0, 0, 0}; // rotate around e2-axis
        fmt::println("B_rot = {:>-7.4f}", B_rot);

        for (size_t nt = 0; nt <= tsteps; nt++) {
            auto const t = t0 + nt * dt;

            // compute in body system
            auto const M_rot = get_motor(B_rot, omega * t); // time dependent rot.
            auto P = move3dp(P0, M_rot);
            auto v = -rcmt(P, B_rot * omega);

            // compute in world system
            auto const M_tra = get_motor(v0 * t);
            P = move3dp(P, M_tra);
            v += v0;

            fmt::println("t = {:>-5.2f},  P = {:>-7.4f},  "
                         "v = {:>-7.4f}",
                         t, P, v);
        }

        fmt::println("");
    }

    TEST_CASE("pga3dp: elementary movement")
    {
        fmt::println("");
        fmt::println("pga3dp: elementary movement");

        // rotational movement with axis through origin
        // (rotational plane and rotational speed):
        //
        // - the weight contains only components which have a factor of e4
        // - directional information of B is contained in its weight
        // - multiplying a vector with e4 only keeps the e4x components
        // - zero weight means B is contained in the horizon
        //
        // - if the vector is a direction of the rotation axis,
        //   B_rot contains this info after the wedge product with e4
        //
        //   B_rot = e4 ^ rot_axis_vec = O_3dp ^ rot_axis_vec
        //                             = O_3dp ^ l_weight_dual(plane of rotation plr)
        //                             = l_cmpl(rwdg(H_3dp, r_bulk_dual(plr)))
        //
        auto n = vec3dp{1, 1, 1, 0}; // only components e1, e2, e3 relevant
        auto plr = r_bulk_dual(n);   // plane with normal n, containing origin
        //
        auto B_rot = wdg(O_3dp, n); // calc from normal (only keeps e4x components)
                                    // only keeps e41, e42 and e43, since e4 ^ e4 == 0
        auto B_rot2 = wdg(O_3dp,
                          l_weight_dual(plr)); // calc from plane directly
        auto B_rot3 = l_cmpl(rwdg(H_3dp, r_cmpl(n)));
        auto B_rot4 = l_cmpl(rwdg(H_3dp, weight(plr)));

        fmt::println("n        -> normal vector = {}", n);
        fmt::println("plr      -> rot. plane    = {}", plr);
        fmt::println("B_rot                     = {}", B_rot);
        fmt::println("");
        fmt::println("B_rot2                    = {}", B_rot2);
        fmt::println("B_rot3                    = {}", B_rot3);
        fmt::println("B_rot4                    = {}", B_rot4);
        fmt::println("");


        // translational movement (for position, speed, acceleration, etc.)
        //
        // - the bulk only contains components which do not contain a factor of e4
        // - positional information of B is contained in its bulk
        // - zero bulk means that B contains the origin
        //
        //   B_tra = att(r_bulk_dual(v)) = att(plane with normal v)
        //
        auto const v = vec3dp{1, 2, 3, 0};
        auto B_tra = att(r_bulk_dual(v));
        // creates a plane that has a normal v (a trivec)
        // att(plane) is a bivector with respective orientation

        fmt::println("v                          = {}", v);
        fmt::println("r_bulk_dual(v)             = {}", r_bulk_dual(v));
        fmt::println("att(r_bulk_dual(v))        = {}", att(r_bulk_dual(v)));
        fmt::println("rwdg(r_bulk_dual(v),H_3dp) = {}", rwdg(r_bulk_dual(v), H_3dp));
        fmt::println("B_tra                      = {}", B_tra);
        fmt::println("");

        // split a bivector rate into the parts carrying info on rotation and translation
        auto B_full = bivec3dp{1, 2, 3, 4, 5, 6};
        fmt::println("B_full = {}", B_full);
        fmt::println("B_full_bulk   = bulk(B_full)   = {}", bulk(B_full));
        fmt::println("B_full_weight = weight(B_full) = {}", weight(B_full));
        fmt::println("");

        auto const B_tot = B_rot + B_tra;
        auto const B_tot_alt = wdg(O_3dp, n) - r_cmpl(wdg(O_3dp, v));

        fmt::println("B_tot = B_rot + B_tra                           = {}", B_tot);
        fmt::println("B_tot_alt = wdg(O_3dp, n) - r_cmpl(wdg(O_3dp, v) = {}", B_tot_alt);
        fmt::println("");

        fmt::println("Derivatives at P_i (must be normlized points with p.w == 1):");
        fmt::println("");
        auto const P_1 = vec3dp{1, 0, 0, 1};
        auto const P_2 = vec3dp{2, 0, 0, 1};
        auto const P_3 = vec3dp{0, 1, 0, 1};
        auto const P_4 = vec3dp{0, 2, 0, 1};
        auto const P_5 = vec3dp{1, 1, 0, 1};
        fmt::println("P_1 = {}", P_1);
        fmt::println("P_2 = {}", P_2);
        fmt::println("P_3 = {}", P_3);
        fmt::println("P_4 = {}", P_4);
        fmt::println("P_5 = {}", P_5);
        fmt::println("");
        fmt::println("Rotational speed:");
        fmt::println("v_rot = rcmt(B_rot,P_1) = {}", rcmt(B_rot, P_1));
        fmt::println("v_rot = rcmt(B_rot,P_2) = {}", rcmt(B_rot, P_2));
        fmt::println("v_rot = rcmt(B_rot,P_3) = {}", rcmt(B_rot, P_3));
        fmt::println("v_rot = rcmt(B_rot,P_4) = {}", rcmt(B_rot, P_4));
        fmt::println("v_rot = rcmt(B_rot,P_5) = {}", rcmt(B_rot, P_5));
        fmt::println("");
        fmt::println("Translational speed:");
        fmt::println("v_tra = rcmt(B_tra,P_1) = {}", rcmt(B_tra, P_1));
        fmt::println("v_tra = rcmt(B_tra,P_2) = {}", rcmt(B_tra, P_2));
        fmt::println("v_tra = rcmt(B_tra,P_3) = {}", rcmt(B_tra, P_3));
        fmt::println("v_tra = rcmt(B_tra,P_4) = {}", rcmt(B_tra, P_4));
        fmt::println("v_tra = rcmt(B_tra,P_5) = {}", rcmt(B_tra, P_5));
        fmt::println("");
        fmt::println("Total speed:");
        fmt::println("v_tot = rcmt(B_tot,P_1) = {}", rcmt(B_tot, P_1));
        fmt::println("v_tot = rcmt(B_tot,P_2) = {}", rcmt(B_tot, P_2));
        fmt::println("v_tot = rcmt(B_tot,P_3) = {}", rcmt(B_tot, P_3));
        fmt::println("v_tot = rcmt(B_tot,P_4) = {}", rcmt(B_tot, P_4));
        fmt::println("v_tot = rcmt(B_tot,P_5) = {}", rcmt(B_tot, P_5));
        fmt::println("");
        fmt::println("Pure rotational speed around reference point P1:");
        fmt::println("v_rot = rcmt(B_rot,P_1) = {}", rcmt(B_rot, P_1 - P_1));
        fmt::println("v_rot = rcmt(B_rot,P_2) = {}", rcmt(B_rot, P_2 - P_1));
        fmt::println("v_rot = rcmt(B_rot,P_3) = {}", rcmt(B_rot, P_3 - P_1));
        fmt::println("v_rot = rcmt(B_rot,P_4) = {}", rcmt(B_rot, P_4 - P_1));
        fmt::println("v_rot = rcmt(B_rot,P_5) = {}", rcmt(B_rot, P_5 - P_1));
        fmt::println("");
        fmt::println("rgpr(B_tra, B_tra) = {}", rgpr(B_tra, B_tra));
        fmt::println("rgpr(B_rot, B_rot) = {}", rgpr(B_rot, B_rot));

        CHECK(rgpr(B_tra, B_tra) == mvec3dp_e{});      // should be zero
                                                       // (stops expansion of exp(B_tra)
                                                       //  after second term)
        CHECK(gr0(rgpr(B_rot, B_rot)) == scalar3dp{}); // should be zero
        CHECK(gr2(rgpr(B_rot, B_rot)) == bivec3dp{});  // should be zero
        CHECK(gr4(rgpr(B_rot, B_rot)) ==
              -pscalar3dp{bulk_nrm_sq(n)}); // should be negative
                                            // and nrm_sq should
                                            // be bulk_nrm_sq(n)

        fmt::println("");
    }

    TEST_CASE("pga3dp: momentum bivector")
    {
        fmt::println("");
        fmt::println("pga3dp: momentum bivector");

        // a mass point at position X (=typically a distance to the center of force)
        // has mass m and speed dX/dt = Xv
        // => total momentum P ist defined as bivector mX^Xv
        //
        // P = X ^ m*Xv = (O + x) ^ m*X_v
        //              = (O ^ m*X_v) + (x ^ m*X_v)
        // where att(P) is the classical linear momentum p = m*v = m*Xv
        // and bulk(P) is the classical angular momentum l = m*(r x v) (as vector)
        //                                               L = X^m*Xv    (as bivector)
        //
        // angular momentum w.r.t. onother point R is given by bivector equilibrium:
        // P_O(X) = X ^ m*Xv
        // P_Y(X) = (X-Y) ^ m*Xv
        // Thus the momentum at any point can be calculated from the given momentum line!
        //
        // => this is completely analogous to the force line and force and torque relation

        // momentum at point X
        auto m = 2.5;
        auto X = vec3dp{2, 0, 0, 1};
        auto Xv = vec3dp{0, 2, 0, 0};

        // other reference points
        auto R1 = vec3dp{1, 0, 0, 1};
        auto R2 = vec3dp{3, 0, 0, 1};

        // total momentum (as bivector)
        auto P = m * wdg(X, Xv);

        fmt::println("m                                    = {}", m);
        fmt::println("X                                    = {}", X);
        fmt::println("Xv                                   = {}", Xv);
        fmt::println("P                                    = {}", P);
        fmt::println("");
        fmt::println("att(P)  = linear momentum            = {}", att(P));
        fmt::println("bulk(P) = angular momentum vs. O_3dp = {}", bulk(P));
        fmt::println("simple blade test: wdg(P,P) == 0:    = {}", wdg(P, P));
        fmt::println("");
        fmt::println("R1                                                = {}", R1);
        fmt::println("bulk(wdg(X - R1,att(P)) = angular momentum vs. R1 = {}",
                     bulk(wdg(X - R1, att(P))));
        fmt::println("R2                                                = {}", R2);
        fmt::println("bulk(wdg(X - R2,att(P)) = angular momentum vs. R2 = {}",
                     bulk(wdg(X - R2, att(P))));

        fmt::println("");
    }

    TEST_CASE("pga3dp: rotating duck - VA, EGA, PGA")
    {
        fmt::println("");
        fmt::println("pga3dp: rotating duck - VA, EGA, PGA");
        fmt::println("");

        // two coordinated systems in relative rotation to each other with
        // a duck moving in the rotating system

        double omega_s = 0.4; // rad/s
        double alpha_s = 0.1; // rad/s^2

        vec3d v_cm = vec3d{0, 0, 0};
        vec3d a_cm = vec3d{0, 0, 0};

        double v_b_s = 0.2; // m/s (speed of duck in rotating system)
        double a_b_s = 0.1; // m/s^2 (acceleration of duck in rotating system)
        double r_b_s = 1.0; // m (position of duck)

        fmt::println("  1st: rotating duck in VA");
        fmt::println("");

        auto omega_v = omega_s * e3_3d;
        auto alpha_v = alpha_s * e3_3d;

        auto v_b_v = v_b_s * e1_3d;
        auto a_b_v = a_b_s * e1_3d;
        auto r_b_v = r_b_s * e1_3d;

        auto v_duck_va = v_cm + cross(omega_v, r_b_v) + v_b_v;
        auto a_duck_va = a_cm + cross(alpha_v, r_b_v) +
                         cross(omega_v, cross(omega_v, r_b_v)) +
                         2.0 * cross(omega_v, v_b_v) + a_b_v;

        fmt::println("VA:   v_duck_va = {}", v_duck_va);
        fmt::println("VA:   a_duck_va = {}", a_duck_va);
        fmt::println("VA:   a_cm  = {}", a_cm);
        fmt::println("VA:   a_t   = {}", cross(alpha_v, r_b_v));
        fmt::println("VA:   a_n   = {}", cross(omega_v, cross(omega_v, r_b_v)));
        fmt::println("VA:   a_cor = {}", 2.0 * cross(omega_v, v_b_v));
        fmt::println("VA:   a_rel = {}", a_b_v);


        fmt::println("");
        fmt::println("  2nd: rotating duck in EGA");
        fmt::println("");

        auto omega_B = dual(omega_v);
        auto alpha_B = dual(alpha_v);

        auto v_duck_ega = v_cm + cmt(r_b_v, omega_B) + v_b_v;
        auto a_duck_ega = a_cm + cmt(r_b_v, alpha_B) + cmt(cmt(r_b_v, omega_B), omega_B) +
                          2.0 * cmt(v_b_v, omega_B) + a_b_v;

        fmt::println("EGA:  v_duck_ega = {}", v_duck_ega);
        fmt::println("EGA:  a_duck_ega = {}", a_duck_ega);
        fmt::println("EGA:  a_cm  = {}", a_cm);
        fmt::println("EGA:  a_t   = {}", cmt(r_b_v, alpha_B));
        fmt::println("EGA:  a_n   = {}", cmt(cmt(r_b_v, omega_B), omega_B));
        fmt::println("EGA:  a_cor = {}", 2.0 * cmt(v_b_v, omega_B));
        fmt::println("EGA:  a_rel = {}", a_b_v);

        fmt::println("");
        fmt::println("  3rd: rotating duck in PGA");
        fmt::println("");

        // vec3dp v_cm_pga = vec3dp{0, 0, 0, 0};
        // vec3dp a_cm_pga = vec3dp{0, 0, 0, 0};
        vec3dp v_cm_pga = vec3dp(v_cm, 0);
        vec3dp a_cm_pga = vec3dp(a_cm, 0);

        auto omega_v_pga = omega_s * e3_3dp;
        auto alpha_v_pga = alpha_s * e3_3dp;

        auto v_b_v_pga = v_b_s * e1_3dp;
        auto a_b_v_pga = a_b_s * e1_3dp;
        auto r_b_v_pga = r_b_s * e1_3dp;

        auto omega_B_pga = wdg(O_3dp, omega_v_pga);
        auto alpha_B_pga = wdg(O_3dp, alpha_v_pga);

        auto v_duck_pga = v_cm_pga + rcmt(omega_B_pga, r_b_v_pga) + v_b_v_pga;
        auto a_duck_pga = a_cm_pga + rcmt(alpha_B_pga, r_b_v_pga) +
                          rcmt(omega_B_pga, rcmt(omega_B_pga, r_b_v_pga)) +
                          2.0 * rcmt(omega_B_pga, v_b_v_pga) + a_b_v_pga;

        fmt::println("PGA:  v_duck_pga = {}", v_duck_pga);
        fmt::println("PGA:  a_duck_pga = {}", a_duck_pga);
        fmt::println("PGA:  a_cm  = {}", a_cm_pga);
        fmt::println("PGA:  a_t   = {}", rcmt(alpha_B_pga, r_b_v_pga));
        fmt::println("PGA:  a_n   = {}", rcmt(omega_B_pga, rcmt(omega_B_pga, r_b_v_pga)));
        fmt::println("PGA:  a_cor = {}", 2.0 * rcmt(omega_B_pga, v_b_v_pga));
        fmt::println("PGA:  a_rel = {}", a_b_v_pga);

        fmt::println("");
    }

    TEST_CASE("VA - ega3d - pga3dp: rotational speed")
    {
        fmt::println("");
        fmt::println("VA - ega3d - pga3dp: rotational speed");

        auto omega3d = vec3d{0, 0, 2 * pi};
        auto r3d = vec3d{2, 0, 0};
        auto v3d = vec3d{-2, 0, 0};
        auto Omega3d = dual(omega3d);

        auto omega3dp = vec3dp{0, 0, 2 * pi, 0};
        auto Q = O_3dp;
        // auto Q = vec3dp{1, -1, 0, 1};
        auto R3dp = vec3dp{2, 0, 0, 1};
        auto v3dp = vec3dp{-2, 0, 0, 0};
        auto Omega3dp = wdg(Q, omega3dp);

        fmt::println("VA:  omega3d = {}", omega3d);
        fmt::println("VA:  v = omega3d x r3d = {}", cross(omega3d, r3d));
        fmt::println("VA:  a_cp = omega3d x omega3d x r3d = {}",
                     cross(omega3d, cross(omega3d, r3d)));
        fmt::println("VA:  a_co = 2*omega3d x v3d = {}", 2 * cross(omega3d, v3d));
        fmt::println("");
        fmt::println("EGA: Omega3d = {}", Omega3d);
        fmt::println("EGA: v = Omega3d >> r3d = {}", Omega3d >> r3d);
        fmt::println("EGA: v = cmt(r3d,Omega3d) = {}", cmt(r3d, Omega3d));
        fmt::println("EGA: a_cp = Omega3d >> (Omega3d >> r3d) = {}",
                     Omega3d >> (Omega3d >> r3d));
        fmt::println("EGA: a_cp = cmt(cmt(r3d,Omega3d),Omega3d) = {}",
                     cmt(cmt(r3d, Omega3d), Omega3d));
        fmt::println("EGA: a_co = 2*cmt(v3d,Omega3d) = {}", 2 * cmt(v3d, Omega3d));
        fmt::println("EGA: a_co = 2*(Omega3d >> v3d) = {}", 2 * (Omega3d >> v3d));
        fmt::println("");
        fmt::println("PGA: Omega3dp = {}", Omega3dp);
        fmt::println("PGA: v = rcmt(Omega3dp, R3dp) = {}", rcmt(Omega3dp, R3dp));
        fmt::println("PGA: a_cp = rcmt(Omega3dp,rcmt(Omega3dp, R3dp)) = {}",
                     rcmt(Omega3dp, rcmt(Omega3dp, R3dp)));
        fmt::println("PGA: a_co = 2*rcmt(Omega3dp, v3dp) = {}", 2 * rcmt(Omega3dp, v3dp));

        fmt::println("");
    }

} // EGA3D PGA3DP Comparison Tests
