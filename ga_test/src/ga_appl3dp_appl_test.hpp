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


/////////////////////////////////////////////////////////////////////////////////////////
// projective geometric algebra pga3dp application tests
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: application tests")
{

    TEST_CASE("pga3dp: reference and tumbling plane (for grinding application)")
    {

        fmt::println("");
        fmt::println("pga3dp: reference and tumbling plane (for grinding application)");
        fmt::println("");

        // define a plane by 3 measuring points A,B,C in coordinate system (e1,e2,e3,e4)
        // with e4 as the projective dimension (e4 = 0 for directions, e4 = 1 for points)
        // (here idealized plane close to e1, e2 plane, derived from fictive measurement)
        auto A = vec3dp(1.0, 0.0, 0.0, 1);
        auto B = vec3dp(1.0, 1.0, 0.01, 1);
        auto C = vec3dp(0.0, 1.0, 0.0, 1);

        // compute the unitized plane and its normal:
        auto p_ref = unitize(wdg(wdg(A, B), C)); // plane equation
        auto p_n = l_weight_dual(p_ref);         // normal of plane

        fmt::println("A = {}", A);
        fmt::println("B = {}", B);
        fmt::println("C = {}", C);
        fmt::println("           p_ref  = {:>-8.5f}", p_ref);
        fmt::println("           p_n    =    {:>-8.5f}, |p_n| = {}", p_n,
                     to_val(bulk_nrm(p_n)));
        fmt::println("");

        // define a reference point T (=center of rotating grinding tool)
        // located close to our measured reference plane
        auto T = vec3dp(3.0, 1.0, 0.0, 1);

        // normal of the tool plane is vector z in e3 direction
        // z is potentially tumbling with radius r at tip in e1, e2 plane
        auto z = vec3dp(0.0, 0.0, 1.0, 0);

        for (auto phi = 0.0; phi < 2 * pi; phi += deg2rad(15)) {

            auto r_dim = 0.01; // radius of tumbling (1% of normal vector)
            auto r = r_dim * vec3dp(cos(phi), sin(phi), 0.0, 0);

            // Point R is at tip of tumbling vector, z_n is normalized vector towards R
            auto R = T + z + r;
            auto z_n = (R - T) / to_val(bulk_nrm(R - T)); // normalized vector from T to
                                                          // R which is the normal of
                                                          // the current plane p_tool

            fmt::println("phi = {:>3.0f}: z_n         =    {:>-8.5f}", rad2deg(phi), z_n);

            // T and z_n form a line from which we can calculate the plane
            // 1st step: line through T in direction of z_n
            auto l_h = unitize(wdg(z_n, T));
            // 2nd step: plane containing T, perpendicular to l
            auto p_tool = unitize(wdg(T, r_weight_dual(l_h)));
            fmt::println("           p_tool      = {:>-8.5f}", p_tool);

            // both planes intersect in a line l
            auto l = rwdg(p_tool, p_ref);
            // print direction vector of line l -> att(l)
            fmt::println("           att(l)      =    {:>-8.5f}", att(l));
            fmt::println("           angle(l,e1) =    {:>-5.2f}°",
                         rad2deg(angle(att(l), e1_3dp)));

            fmt::println("");
        }
    }

    TEST_CASE("pga3dp: intersecting discs (for grinding application)")
    {

        fmt::println("");
        fmt::println("pga3dp: intersecting discs (for grinding application)");
        fmt::println("");

        // to get started simple: no thickness of disc modelled, but assumed to have
        // material thickness in direction against normal vector of modelled plane
        // ctor interface needs to be improved for practical use (distinguishable types)
        struct disc {

            disc(vec3dp const& cp_in, double d_in, double r_min_in, double r_max_in,
                 double rs_in, vec3dp const& u_in, vec3dp const& v_in) :
                cp(cp_in), r(d_in / 2.0), r_min(r_min_in), r_max(r_max_in), rs(rs_in),
                u(u_in), v(v_in), B_uv(wdg(u, v)), pl(wdg(cp, B_uv)),
                pl_normal(l_weight_dual(pl))
            {
                if (cp.w != 1.0)
                    throw std::runtime_error("disc: cp must be a point (w=1).");
                if (d_in <= 0.0) throw std::runtime_error("disc: invalid diameter.");
                if (r_min < 0.0 || r_max < 0.0 || r_min > r_max || r > r_max)
                    throw std::runtime_error("disc: invalid radius.");
                if (dot(u, v) != 0.0)
                    throw std::runtime_error("disc: u and v must be perpendicular.");
                if (bulk_nrm(u) != 1.0 || u.w != 0.0)
                    throw std::runtime_error("disc: u must be a unit-vector.");
                if (bulk_nrm(v) != 1.0 || v.w != 0.0)
                    throw std::runtime_error("disc: v must be a unit-vector.");
            }

            vec3dp cp;    // center point [mm]
            double r;     // radius [mm]
            double r_min; // min radius [mm] (for tool with ring shape)
            double r_max; // max radius [mm] (for tool with ring shape)
            double rs;    // rotating speed [rad/s]
            vec3dp u;     // unit vector in plane of disc (local cs, u-direction)
            vec3dp v; // unit vector in plane of disc (local cs, v-direction, perp. to u)

            bivec3dp B_uv;    // bivector (calculated)
            trivec3dp pl;     // plane (calculated)
            vec3dp pl_normal; // plane normal (calculated)
        };

        auto u = vec3dp{1, 0, 0, 0}; // unit vector in plane of disc
        auto v = vec3dp{0, 1, 0, 0}; // unit vector in plane of disc

        // wafer is located on chuck
        // normal of plane of wafer showing upwards (towards e3)
        // center point (=reference point) is on top of wafer -> here: cp at origin
        // (wafer thickness would be going downwards, i.e. against direction of normal)
        // wafer plane is identical with e12 plane, rotation axis is e3
        // wafer diameter is 200 mm
        // wafer rotation speed is 300 rpm (defined by rotation speed of chuck)
        // u and v define right-handed local coordinate system of disc
        // (here identical with global coordinate system e1,e2 for untilted disc)
        auto wafer = disc(O_3dp, 200.0, 0.0, 100.0, rpm2radps(300.0), u, v);

        fmt::println("wafer.cp         = {}", wafer.cp);
        fmt::println("wafer.r [mm]     = {}", wafer.r);
        fmt::println("wafer.r_min [mm] = {}", wafer.r_min);
        fmt::println("wafer.r_max [mm] = {}", wafer.r_max);
        fmt::println("wafer.rs [rad/s] = {}", wafer.rs);
        fmt::println("wafer.pl         = {:>-8.5f}", wafer.pl);
        fmt::println("wafer.pl_normal  = {:>-8.5f}", wafer.pl_normal);
        fmt::println("u                = {:>-8.5f}", u);
        fmt::println("v                = {:>-8.5f}", v);
        fmt::println("");

        // tool center
        auto tc = vec3dp{150, 0, -0.5, 1}; // overlap is 0.5 mm in z direction
        // modify u to tilt the tool-plane (+e3-comp. => volume becomes larger)
        //                                 (-e3-comp. => volume becomes smaller)
        // u += vec3dp{0.0, 0.0, 0.005, 0};
        // u = u / to_val(bulk_nrm(u));
        //
        // wafer plane is identical with -e12 plane
        // tool is located at tc
        // tool diameter is 200 mm (r=100 mm)
        // tool min radius is 95 mm (ring shape)
        // tool max radius is 100 mm
        // tool rotation speed is 2000 rpm (typical value for grinding is 1500-3000 rpm)
        // u and v define right-handed local coordinate system of disc
        auto tool = disc(tc, 200.0, 95.0, 100.0, rpm2radps(2000), u, -v);

        fmt::println("tool.cp         = {}", tool.cp);
        fmt::println("tool.r [mm]     = {}", tool.r);
        fmt::println("tool.r_min [mm] = {}", tool.r_min);
        fmt::println("tool.r_max [mm] = {}", tool.r_max);
        fmt::println("tool.rs [rad/s] = {}", tool.rs);
        fmt::println("tool.pl         = {:>-8.5f}", tool.pl);
        fmt::println("tool.pl_normal  = {:>-8.5f}", tool.pl_normal);
        fmt::println("u               = {:>-8.5f}", u);
        fmt::println("v               = {:>-8.5f}", -v);
        fmt::println("");

        // begin detour: just to show how it works, not directly needed here
        //
        // shift a plane with its reference point in a direction x, while keeping the same
        // orientation:
        // 1.) pl = P ^ u ^ v
        // 2.) pl_shifted = P ^ u ^ v + x ^ u ^ v = (P + x) ^ (u ^ v)
        auto x = vec3dp{0, 0, 0.05, 0}; // vector to shift the plane
        auto pl_shifted = wafer.pl + wdg(x, wafer.B_uv);
        CHECK(pl_shifted == wdg(wafer.cp + x, wafer.B_uv));
        // end detour

        auto volume = 0.0;
        auto area = 0.0;
        double min_tot_speed = std::numeric_limits<double>::max();
        double max_tot_speed = std::numeric_limits<double>::min();

        if (bulk_nrm(tool.cp - wafer.cp) < tool.r + wafer.r) {

            // calculate overlapping volume, if both discs are in reach of each other
            // simple version (not fully efficent): always use full tool area
            // r- and phi-range could be reduced to actually needed range with a little
            // more effort

            size_t nr = 25; // number of steps in radial direction
            double dr = (tool.r_max - tool.r_min) / nr;

            size_t nphi = 720;             // number of steps in angular direction
            double dphi = 2.0 * pi / nphi; // angle in rad

            fmt::println("nr = {}, dr = {:>-8.5f}, nphi = {}, dphi = {:>-8.5f}", nr, dr,
                         nphi, dphi);
            fmt::println("");

            for (size_t j = 0; j < nr; ++j) {

                // integrate from r_min to r_max only
                double r = tool.r_min + j * dr; // current inner radius
                double rm = r + dr * 0.5;       // mean radius for calculation of area

                for (size_t i = 0; i < nphi; ++i) {

                    double phi = i * dphi;          // current angle
                    double phim = phi + dphi * 0.5; // mean angle f. calc. geom. centroid

                    double dA = rm * dphi * dr;
                    // fmt::println("r = {}, phi = {:>-8.5f}, da = {:>-8.5f}", r, phi,
                    // dA);

                    // calculate position of current geometric centroid
                    auto r_i_sq = r * r;               // inner radius r_i^2
                    auto r_i_tr = r_i_sq * r;          // inner radius r_i^3
                    auto r_o_sq = (r + dr) * (r + dr); // outer radius r_o^2
                    auto r_o_tr = r_o_sq * (r + dr);   // outer radius r_o^3

                    // distance to geometric centroid of ring segment (annulus)
                    auto r_gc = 2.0 * (r_o_tr - r_i_tr) * sin(dphi * 0.5) /
                                (3.0 * (r_o_sq - r_i_sq) * dphi * 0.5);
                    auto ref_pos = tool.cp + r_gc * u;

                    // rotation line is defined by tool center and tool normal
                    auto rot_line = wdg(tool.cp, tool.pl_normal);

                    // move from reference position to position defined by rot_line & phim
                    auto cur_pos = move3dp(ref_pos, get_motor(rot_line, phim));

                    // fmt::println("ref_pos = {:>-8.5f}", ref_pos);
                    // fmt::println("cur_pos = {:>-8.5f}", cur_pos);

                    // orthogonal projection of current position onto wafer plane
                    auto cur_proj_pos =
                        rwdg(wafer.pl, r_weight_expand3dp(cur_pos, wafer.pl));
                    // fmt::println("cur_proj_pos = {:>-8.5f}", cur_proj_pos);

                    // if cur_proj_pos is within the wafer, potentially add to area and
                    // volume
                    if (bulk_nrm(cur_proj_pos - wafer.cp) <= wafer.r) {
                        // now the two planes are in reach of each other,
                        // but we want to make sure that area or volume is only counted,
                        // if the planes really overlap

                        auto delta_vec = cur_pos - cur_proj_pos;
                        // fmt::println("delta_vec = {}", delta_vec);

                        if (sign(dot(delta_vec, tool.pl_normal)) == 1) {
                            // only if delta vec is aligned with the direction of tool
                            // normal vector there is a real overlap of the materials of
                            // wafer and tool

                            area += dA;
                            volume += dA * to_val(bulk_nrm(delta_vec));

                            // calculate total speed at current position
                            auto wafer_speed =
                                wafer.rs * wafer.B_uv >> (cur_proj_pos - wafer.cp);

                            auto tool_speed = tool.rs * tool.B_uv >> (cur_pos - tool.cp);
                            auto tot_speed_vec = tool_speed + wafer_speed;
                            auto tot_speed = to_val(bulk_nrm(tot_speed_vec));
                            // fmt::println("wafer_speed = {:>-8.5f}, nrm = {:>-8.5f}",
                            //              wafer_speed, to_val(bulk_nrm(wafer_speed)));
                            // fmt::println("tool_speed  = {:>-8.5f}, nrm = {:>-8.5f}",
                            //              tool_speed, to_val(bulk_nrm(tool_speed)));
                            // fmt::println("tot_speed_vec = {:>-8.5f}, nrm = {:>-8.5f}",
                            //              tot_speed_vec, tot_speed);
                            // fmt::println("");

                            if (tot_speed != 0.0) {
                                min_tot_speed = std::min(min_tot_speed, tot_speed);
                                max_tot_speed = std::max(max_tot_speed, tot_speed);
                            }
                        }
                    }
                }
            }
        }

        // double tool_area = pi * tool.r * tool.r;
        // fmt::println("area = {:>-8.5f} mm^2, tool_area = {:>-8.5f} mm^2", area,
        //              tool_area);
        fmt::println("area = {:>-8.5f} mm^2", area);
        fmt::println("volume = {:>-8.5f} mm^3", volume);
        fmt::println("min_tot_speed = {:>-8.5f} mm/s", min_tot_speed);
        fmt::println("max_tot_speed = {:>-8.5f} mm/s", max_tot_speed);
        fmt::println("");
        fmt::println("wafer speed at r = {} mm = {:>-8.5f} mm/s", wafer.r_max,
                     wafer.rs * wafer.r_max); // v = omega * r
        fmt::println("tool  speed at r = {} mm = {:>-8.5f} mm/s", tool.r_max,
                     tool.rs * tool.r_max); // v = omega * r
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // The two tests below re-express the two grinding mini-tests above using the
    // static_/kinematic_system3dp frame-tree infrastructure. They bring NO new geometric
    // insight -- the point is purely pedagogical: to show how the same very simple idea
    // reads once it is phrased through the framework (define geometry in a local frame,
    // let the system carry it to world coordinates / ask the system for a point's
    // velocity) rather than with hand-built duals and bivector velocity fields. Each
    // asserts its result against the raw-PGA computation so the equivalence is explicit.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: tumbling tool plane via static_system3dp (grinding application)")
    {

        fmt::println("");
        fmt::println("pga3dp: tumbling tool plane via static_system3dp (grinding "
                     "application)");
        fmt::println("");

        // Same question as "reference and tumbling plane": where does the tumbling tool
        // plane cut the wafer reference plane? Instead of rebuilding the tilted tool
        // normal and its plane via duals at every angle, we attach the tool's cutting
        // plane to a tool FRAME once, and let static_system3dp carry it to world
        // coordinates as the frame is re-posed. Tumbling = a precessing tilt of that
        // frame.

        // reference (wafer) plane from 3 measured points (unchanged from the raw test)
        auto A = vec3dp(1.0, 0.0, 0.0, 1);
        auto B = vec3dp(1.0, 1.0, 0.01, 1);
        auto C = vec3dp(0.0, 1.0, 0.0, 1);
        auto p_ref = unitize(wdg(wdg(A, B), C));

        // tool reference point (center of the rotating grinding tool)
        auto T = vec3dp(3.0, 1.0, 0.0, 1);
        double const r_dim = 0.01; // tumbling radius (1% of the unit normal)

        // frame tree: world/wafer root W, tool frame T_f posed relative to it
        static_system3dp sys;
        sys.add_frame(static_frame3dp("W"));      // inertial root (wafer side)
        sys.add_frame(static_frame3dp("T_f", T)); // tool frame, initially upright
        size_t const W_idx = sys.index_of("W");
        size_t const T_idx = sys.index_of("T_f");

        // the tool cutting plane in its OWN frame: the horizontal plane through the frame
        // origin spanned by local e1, e2 (normal = local e3). Defined once.
        auto const pl_local = wdg(O_3dp, wdg(vec3dp{1, 0, 0, 0}, vec3dp{0, 1, 0, 0}));

        for (auto phi = 0.0; phi < 2 * pi; phi += deg2rad(15)) {

            // tumbling = tilt the tool frame by atan(r_dim) about the axis perpendicular
            // to the momentary tilt direction (cos phi, sin phi); the axis precesses with
            // phi.
            auto const tilt_axis = vec3dp{-std::sin(phi), std::cos(phi), 0.0, 0.0};
            auto const rotvec =
                std::atan(r_dim) * tilt_axis; // axis-angle rotation vector
            sys.set_pose(T_idx, pose3dp(T, rotvec));

            // tool plane in world: carry the local plane through the frame tree
            auto const M = sys.get_pos_trafo(T_idx, W_idx);
            auto const p_tool = unitize(move3dp(pl_local, M));

            // intersection line of the two planes
            auto const l = rwdg(p_tool, p_ref);

            // raw-PGA reference computation (as in the original test) for cross-check
            auto const r = r_dim * vec3dp(std::cos(phi), std::sin(phi), 0.0, 0);
            auto const z = vec3dp(0.0, 0.0, 1.0, 0);
            auto const R = T + z + r;
            auto const z_n = (R - T) / to_val(bulk_nrm(R - T));
            auto const l_h = unitize(wdg(z_n, T));
            auto const p_tool_raw = unitize(wdg(T, r_weight_dual(l_h)));
            auto const l_raw = rwdg(p_tool_raw, p_ref);

            // same tool plane and same intersection line, up to scale / orientation
            CHECK(is_congruent(p_tool, p_tool_raw));
            CHECK(is_congruent(l, l_raw));

            fmt::println("phi = {:>3.0f}: angle(l,e1) = {:>5.2f}° (framework) == "
                         "{:>5.2f}° (raw)",
                         rad2deg(phi), rad2deg(angle(att(l), e1_3dp)),
                         rad2deg(angle(att(l_raw), e1_3dp)));
        }
        fmt::println("");
    }

    TEST_CASE(
        "pga3dp: disc surface speeds via kinematic_system3dp (grinding application)")
    {

        fmt::println("");
        fmt::println("pga3dp: disc surface speeds via kinematic_system3dp (grinding "
                     "application)");
        fmt::println("");

        // The "intersecting discs" test computes the relative grinding speed as a raw
        // bivector velocity field  rs * B_uv >> (P - cp). Here the same surface speeds
        // come from kinematic_system3dp::point_velocity of two spinning frames -- the
        // idiomatic way to ask "how fast does this world point move because its frame
        // rotates". The overlap / MRR integration of the raw test is unchanged and not
        // repeated; only the speed query is re-expressed and cross-checked against the
        // raw formula.

        double const rs_wafer = rpm2radps(300.0);      // chuck/wafer spin [rad/s]
        double const rs_tool = rpm2radps(2000.0);      // grinding-wheel spin [rad/s]
        auto const tc = vec3dp{150.0, 0.0, -0.5, 1.0}; // tool center (0.5 mm overlap)
        double const rm = 97.5; // mean ring radius of the tool [mm]

        // spinning frames: wafer about +e3 through the origin; tool about +e3 through tc.
        kinematic_system3dp ksys;
        ksys.add_frame(static_frame3dp("W")); // inertial root
        ksys.add_frame(static_frame3dp("wafer", O_3dp),
                       kin_state3dp{.omega = vec3dp{0.0, 0.0, rs_wafer, 0.0}},
                       ksys.index_of("W"));
        ksys.add_frame(static_frame3dp("tool", tc),
                       kin_state3dp{.omega = vec3dp{0.0, 0.0, rs_tool, 0.0}},
                       ksys.index_of("W"));

        // raw bivector spin rates for the cross-check (untilted discs, same +e3 sense)
        auto const B_uv = wdg(vec3dp{1, 0, 0, 0}, vec3dp{0, 1, 0, 0});
        auto const Omega_wafer = rs_wafer * B_uv;
        auto const Omega_tool = rs_tool * B_uv;

        double min_tot = std::numeric_limits<double>::max();
        double max_tot = std::numeric_limits<double>::lowest();

        for (auto a = 0.0; a < 2 * pi; a += deg2rad(45)) {

            // a sample point on the tool ring, and its projection onto the wafer (z = 0)
            auto const P_tool = tc + rm * vec3dp{std::cos(a), std::sin(a), 0.0, 0.0};
            auto const P_proj = vec3dp{P_tool.x, P_tool.y, 0.0, 1.0};

            // framework: surface speeds from the spinning frames
            auto const v_tool = ksys.point_velocity(P_tool, "tool");
            auto const v_wafer = ksys.point_velocity(P_proj, "wafer");
            auto const tot_fw = to_val(bulk_nrm(v_tool + v_wafer));

            // raw: same speeds as the bivector velocity field of the original test
            auto const v_tool_raw = (Omega_tool >> (P_tool - tc));
            auto const v_wafer_raw = (Omega_wafer >> (P_proj - O_3dp));
            auto const tot_raw = to_val(bulk_nrm(v_tool_raw + v_wafer_raw));

            // tool surface speed magnitude is omega * (distance from spin axis) = rs * rm
            CHECK(to_val(bulk_nrm(v_tool)) == doctest::Approx(rs_tool * rm));
            // framework and raw relative-speed magnitudes agree (same physics)
            CHECK(tot_fw == doctest::Approx(tot_raw));

            min_tot = std::min(min_tot, tot_fw);
            max_tot = std::max(max_tot, tot_fw);

            fmt::println("a = {:>3.0f}: |v_tool| = {:>8.2f}, |v_wafer| = {:>8.2f}, "
                         "|v_rel| = {:>8.2f} mm/s",
                         rad2deg(a), to_val(bulk_nrm(v_tool)), to_val(bulk_nrm(v_wafer)),
                         tot_fw);
        }
        fmt::println("min_tot_speed = {:>-8.2f} mm/s", min_tot);
        fmt::println("max_tot_speed = {:>-8.2f} mm/s", max_tot);
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Phase 0 of the wafer-grinding plan (TODO/grinding.md): the geometric on-ramp. Build
    // the Fig.-1 two-chain frame tree (chuck/wafer + spindle/wheel) and reproduce the
    // GRAIN TRAJECTORY -- the path a wheel-rim grain traces in the ROTATING WAFER frame
    // (Fig. 1 grey-dashed / Fig. 7 red-dashed marks). Pure kinematics, no dynamics: it
    // validates the frame tree end-to-end before any force elements are added, and emits
    // the e423_3dp (chuck-plane) projection that ga_view will draw.
    //
    // Spin naming follows Fig. 1 (counterintuitive): n_s = chuck (slow), n_w = wheel
    // (fast). Nominal radial placement: the spindle CM is offset (-R/v2, +R/v2) in the
    // root e2-e3 plane so its distance from the chuck axis is exactly R; with wheel
    // radius R the rim reaches the wafer centre and sweeps one wafer radius (self-rot
    // layout).
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: grinding grain trajectory in the wafer frame (Phase 0)")
    {

        fmt::println("");
        fmt::println("pga3dp: grinding grain trajectory in the wafer frame (Phase 0)");
        fmt::println("");

        // machine geometry (Tao Fig. 1; see TODO/grinding.md)
        double const R = 150.0;    // wheel radius == wafer radius [mm] (equal for now)
        double const l3 = 30.0;    // tool surface distance from the spindle CM [mm]
        double const tw_avg = 0.8; // average wafer thickness [mm]
        double const x_a = l3;     // axial infeed so the grain sits at the chuck plane

        // spin rates per Fig. 1: chuck n_s (slow), wheel n_w (fast). n_w / n_s = 10 here,
        // so the grinding-mark pattern closes after one chuck revolution.
        double const ns = rpm2radps(300.0);  // chuck spin [rad/s]
        double const nw = rpm2radps(3000.0); // wheel spin [rad/s]

        // nominal radial placement of the spindle CM: |offset|_yz == R (see header note)
        double const off = R / std::sqrt(2.0);
        auto const spindle_origin = vec3dp{x_a, -off, off, 1.0}; // (x_a, -R/v2, +R/v2)

        // frame tree
        kinematic_system3dp sys;
        sys.add_frame(static_frame3dp("chuck_ctr_stat")); // inertial root

        // chain 1 (wafer): chuck spin about e1, then the average wafer top surface
        sys.add_frame(static_frame3dp("chuck_ctr_rot"),
                      kin_state3dp{.omega = vec3dp{ns, 0.0, 0.0, 0.0}},
                      sys.index_of("chuck_ctr_stat"));
        sys.add_frame(static_frame3dp("wafer_top_avg_rot", vec3dp{tw_avg, 0.0, 0.0, 1.0}),
                      kin_state3dp{}, sys.index_of("chuck_ctr_rot"));

        // chain 2 (spindle/wheel): CM reoriented -90 deg about e2 (tool axis -> -e1), the
        // tool surface at +l3 along the spindle z spinning at n_w, rim grain at +R along
        // x
        sys.add_frame(static_frame3dp("spindle_cm_stat", spindle_origin,
                                      vec3dp{0.0, -pi / 2.0, 0.0, 0.0}),
                      kin_state3dp{}, sys.index_of("chuck_ctr_stat"));
        sys.add_frame(static_frame3dp("tool_top_avg_rot", vec3dp{0.0, 0.0, l3, 1.0}),
                      kin_state3dp{.omega = vec3dp{0.0, 0.0, nw, 0.0}},
                      sys.index_of("spindle_cm_stat"));
        sys.add_frame(static_frame3dp("tool_surface_avg_at_R", vec3dp{R, 0.0, 0.0, 1.0}),
                      kin_state3dp{}, sys.index_of("tool_top_avg_rot"));

        size_t const chuck_rot = sys.index_of("chuck_ctr_rot");
        size_t const tool_top = sys.index_of("tool_top_avg_rot");

        // wheel-axis location in the root e2-e3 (chuck) plane: fixed (independent of both
        // spins). Its distance from the chuck axis must be exactly R -- the key geometry.
        auto const C = unitize(
            move3dp(O_3dp, sys.get_pos_trafo("tool_top_avg_rot", "chuck_ctr_stat")));
        CHECK(std::sqrt(C.y * C.y + C.z * C.z) == doctest::Approx(R));

        // sample the grain over one chuck revolution
        double const T = 2.0 * pi / ns; // one chuck revolution [s]
        int const nsteps = 720;
        double r_min = std::numeric_limits<double>::max();
        double r_max = std::numeric_limits<double>::lowest();

        for (int k = 0; k <= nsteps; ++k) {
            double const t = T * double(k) / double(nsteps);

            // advance the two spins (every other frame co-rotates with its parent)
            sys.set_pose(chuck_rot, pose3dp{O_3dp, vec3dp{ns * t, 0.0, 0.0, 0.0}});
            sys.set_pose(tool_top, pose3dp{vec3dp{0.0, 0.0, l3, 1.0},
                                           vec3dp{0.0, 0.0, nw * t, 0.0}});

            // the rim grain (origin of the rim frame) in the rotating wafer frame and
            // root
            auto const g_wafer = unitize(move3dp(
                O_3dp, sys.get_pos_trafo("tool_surface_avg_at_R", "wafer_top_avg_rot")));
            auto const g_root = unitize(move3dp(
                O_3dp, sys.get_pos_trafo("tool_surface_avg_at_R", "chuck_ctr_stat")));

            // radial distance from the chuck axis (e1) is invariant under the e1 spin, so
            // it must agree in the root and wafer frames
            double const r_wafer =
                std::sqrt(g_wafer.y * g_wafer.y + g_wafer.z * g_wafer.z);
            double const r_root = std::sqrt(g_root.y * g_root.y + g_root.z * g_root.z);
            CHECK(r_wafer == doctest::Approx(r_root));

            // the grain stays on the wheel rim: distance R from the fixed wheel axis C
            double const dy = g_root.y - C.y;
            double const dz = g_root.z - C.z;
            CHECK(std::sqrt(dy * dy + dz * dz) == doctest::Approx(R));

            r_min = std::min(r_min, r_wafer);
            r_max = std::max(r_max, r_wafer);

            // e423_3dp projection (drop the e1 component): the chuck-plane (y,z) the mark
            // is drawn in -- this is the column ga_view will plot
            if (k % (nsteps / 12) == 0)
                fmt::println("t/T = {:>4.2f}: grain_wafer (e423 y,z) = ({:>8.2f}, "
                             "{:>8.2f}) mm,  r = {:>7.2f} mm",
                             t / T, g_wafer.y, g_wafer.z, r_wafer);
        }

        // the rim sweeps from the wafer centre (r ~ 0) out to 2R (self-rot grinding
        // layout)
        CHECK(r_min < 0.1 * R);
        CHECK(r_max == doctest::Approx(2.0 * R).epsilon(0.001));
        fmt::println("r range over one chuck rev: [{:>5.2f}, {:>6.2f}] mm  (expect ~[0, "
                     "{:.0f}])",
                     r_min, r_max, 2.0 * R);

        // rim grain speed from the wheel spin: |v| = n_w * R (distance from the wheel
        // axis)
        auto const g_now = unitize(
            move3dp(O_3dp, sys.get_pos_trafo("tool_surface_avg_at_R", "chuck_ctr_stat")));
        auto const v_grain = sys.point_velocity(g_now, "tool_surface_avg_at_R");
        CHECK(to_val(bulk_nrm(v_grain)) == doctest::Approx(nw * R));
        fmt::println("rim grain speed |v| = {:>10.2f} mm/s  (expect n_w*R = {:>10.2f})",
                     to_val(bulk_nrm(v_grain)), nw * R);
        fmt::println("");
    }

    // Wafer-axis TILT -> global flatness (Zhou et al., Precision Eng. 27 (2003) 175):
    // dropping the parallel-axis idealization. Our flat-surface assumption holds only
    // when the wheel axis is exactly parallel to the wafer axis. In reality the wafer
    // axis is slightly misaligned by two small alignment angles -- alpha about X, beta
    // about Y (Zhou Eq.1, matrices B and C) -- and that tilt is what sets the GLOBAL
    // FLATNESS of the ground wafer (distinct from the WAVINESS the Phase-D.1 axial-runout
    // model carves). This reproduces Zhou Figs 5-7: a single wheel grain, sampled in the
    // (tilted, spinning) wafer frame over a common period, traces a 3D cutting path whose
    // height z (the wafer-axis component) is a CONE for an alpha-tilt, a CONVEX dome for
    // +beta and a CONCAVE bowl for -beta -- flat only when alpha = beta = 0.
    //
    // Zhou Eq.1 is r2 = A B C D [r1, L, -f]^T (matrix chain). The same kinematics in PGA:
    // the wheel grain (radius r1, spin omega1, axis offset L) expressed in the wafer
    // frame (spin omega2, static tilt alpha,beta) via a motor sandwich -- our
    // kinematic_system3dp get_pos_trafo. The cutting PATH (not yet the removed-material
    // envelope or the density-driven centre concavity of Zhou Figs 8-16) is what Figs 5-7
    // plot.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: wafer-tilt flatness profile, Zhou Figs 5-7 (cone / dome / bowl)")
    {
        fmt::println("");
        fmt::println(
            "pga3dp: wafer-tilt flatness profile, Zhou Figs 5-7 (cone/dome/bowl)");
        fmt::println("");

        // Zhou geometry: R1 = R2 = 150 mm, L = R1 (half-overlap, rim through wafer
        // centre), grain at (r1, theta1) = (150, 0). Speed ratio n2/n1 = 1/30 (wafer /
        // wheel), as in the Fig.5 rosette. Root: e1 = wafer/wheel axis (Zhou Z), e2-e3 =
        // grinding plane (Zhou X-Y); the wheel-axis offset L is along +e2 (Zhou X).
        double const R = 150.0;              // wheel radius == wafer radius [mm]
        double const L = R;                  // wheel-axis offset (half-overlap)
        double const r1 = R;                 // grain radius on the wheel
        double const w1 = rpm2radps(1500.0); // wheel spin n1 [rad/s]
        double const w2 = rpm2radps(50.0);   // wafer spin n2 [rad/s]  (n2/n1 = 1/30)

        // Sample one grain in the wafer frame over one wafer revolution (the single-grain
        // cutting path of Zhou Figs 5-7) for a given pair of alignment tilts and return
        // the profile-height (z = wafer-axis component) statistics + the radial extent.
        auto profile = [&](double a_deg, double b_deg) {
            double const a = a_deg * pi / 180.0;
            double const b = b_deg * pi / 180.0;

            kinematic_system3dp sys;
            sys.add_frame(static_frame3dp("W")); // inertial root (wheel-axis aligned)
            // wheel: centre offset L along +e2, spins about its own axis e1 at n1
            sys.add_frame(static_frame3dp("wheel_rot", vec3dp{0.0, L, 0.0, 1.0}),
                          kin_state3dp{.omega = vec3dp{w1, 0.0, 0.0, 0.0}},
                          sys.index_of("W"));
            // grain at radius r1 from the wheel centre; placed so it reaches the wafer
            // centre (L - r1 = 0) at t = 0
            sys.add_frame(static_frame3dp("grain", vec3dp{0.0, -r1, 0.0, 1.0}),
                          kin_state3dp{}, sys.index_of("wheel_rot"));
            // wafer: a STATIC misalignment tilt (alpha about e2, beta about e3) carrying
            // a frame that spins about the tilted axis at n2 -- Zhou's B*C then A
            sys.add_frame(static_frame3dp("wafer_tilt", O_3dp, vec3dp{0.0, a, b, 0.0}),
                          kin_state3dp{}, sys.index_of("W"));
            sys.add_frame(static_frame3dp("wafer_rot"),
                          kin_state3dp{.omega = vec3dp{w2, 0.0, 0.0, 0.0}},
                          sys.index_of("wafer_tilt"));

            size_t const wheel = sys.index_of("wheel_rot");
            size_t const wafer = sys.index_of("wafer_rot");
            double const T = 2.0 * pi / w2; // one wafer revolution
            int const N = 6000;
            struct {
                double zmin = 1e9, zmax = -1e9, z_ctr = 0.0, rmax = 0.0;
            } p;
            for (int k = 0; k <= N; ++k) {
                double const t = T * double(k) / double(N);
                sys.set_pose(wheel, pose3dp{vec3dp{0.0, L, 0.0, 1.0},
                                            vec3dp{w1 * t, 0.0, 0.0, 0.0}});
                sys.set_pose(wafer, pose3dp{O_3dp, vec3dp{w2 * t, 0.0, 0.0, 0.0}});
                // grain in the wafer frame: z = e1 (wafer-axis) component = profile
                // height
                auto const g =
                    unitize(move3dp(O_3dp, sys.get_pos_trafo("grain", "wafer_rot")));
                double const z = g.x;
                double const r = std::sqrt(g.y * g.y + g.z * g.z);
                // the grain orbits out to 2R, but only marks the wafer where it is ON it
                // (r2 <= R2) -- Zhou clips the cutting path to the wafer in Figs 5-7
                if (r > R + 1e-6) continue;
                p.zmin = std::min(p.zmin, z);
                p.zmax = std::max(p.zmax, z);
                if (r < 3.0) p.z_ctr = z;
                p.rmax = std::max(p.rmax, r);
            }
            return p;
        };

        double const tol = 1.0e-9; // the chain is exact -> tight tolerances
        double const amp = R * std::tan(0.1 * pi / 180.0); // 150*tan(0.1 deg) = 0.2618 mm

        // (1) perfect alignment alpha = beta = 0 -> FLAT (z == 0 everywhere)
        auto const flat = profile(0.0, 0.0);
        CHECK(std::abs(flat.zmin) < tol);
        CHECK(std::abs(flat.zmax) < tol);
        fmt::println("alpha=0.0 beta= 0.0 : z in [{:+.4f},{:+.4f}] mm  -> FLAT",
                     flat.zmin, flat.zmax);

        // The exact clipped rim amplitude is not a closed form (it depends on which wafer
        // points the discrete cutting path reaches), so gate the SHAPE -- symmetric cone
        // / one-sided dome / mirror bowl -- and a sensible amplitude band, not an exact
        // value. The reference R*tan(alpha) = 0.2618 mm is the unclipped rim amplitude.

        // (2) alpha = 0.1, beta = 0 -> CONE: symmetric about z = 0, apex at the centre,
        // reaching the wafer rim with an amplitude of order R*tan(alpha)
        auto const cone = profile(0.1, 0.0);
        CHECK(std::abs(cone.zmax + cone.zmin) < 1e-3 * amp); // symmetric about z = 0
        CHECK(std::abs(cone.z_ctr) < 1e-3);                  // apex at the centre
        CHECK(cone.zmax > 0.5 * amp);                        // a real cone (~0.2 mm)
        CHECK(cone.zmax < 1.01 * amp);                       // bounded by R*tan(alpha)
        CHECK(cone.rmax > 0.99 * R);                         // path reaches the wafer rim
        fmt::println(
            "alpha=0.1 beta= 0.0 : z in [{:+.4f},{:+.4f}] mm  -> CONE (symmetric)",
            cone.zmin, cone.zmax);

        // (3) beta = +0.1 -> CONVEX dome: z stays >= 0 (one-sided, bulges up)
        auto const dome = profile(0.0, 0.1);
        CHECK(dome.zmin == doctest::Approx(0.0).epsilon(1e-4));
        CHECK(dome.zmax > 0.3 * amp); // a real dome
        // the on-wafer beta-dome amplitude comes out at ABOUT HALF the alpha-cone
        // amplitude -- reproducing Zhou's observation that "the effect of beta on the
        // global flatness is about half that of alpha". It emerges from the clipping
        // geometry, not put in by hand.
        CHECK(dome.zmax < 0.75 * cone.zmax);
        fmt::println("alpha=0.0 beta=+0.1 : z in [{:+.4f},{:+.4f}] mm  -> CONVEX (dome); "
                     "dome/cone = {:.2f} (Zhou: beta ~ half alpha)",
                     dome.zmin, dome.zmax, dome.zmax / cone.zmax);

        // (4) beta = -0.1 -> CONCAVE bowl: z stays <= 0 (the mirror of the dome)
        auto const bowl = profile(0.0, -0.1);
        CHECK(bowl.zmax == doctest::Approx(0.0).epsilon(1e-4));
        CHECK(bowl.zmin < -0.3 * amp);
        CHECK(bowl.zmin == doctest::Approx(-dome.zmax).epsilon(1e-4)); // mirror of +beta
        fmt::println("alpha=0.0 beta=-0.1 : z in [{:+.4f},{:+.4f}] mm  -> CONCAVE (bowl)",
                     bowl.zmin, bowl.zmax);

        // global-flatness magnitude: 0.1 deg of alpha gives a rim amplitude of
        // R*tan(0.1 deg) = 0.262 mm over the 150 mm radius, matching Zhou's "0.1 deg in
        // alpha creates a variation of about 0.2 mm over 300 mm".
        CHECK(amp == doctest::Approx(0.2618).epsilon(1e-3));
        fmt::println("");
        fmt::println("global-flatness rim amplitude at 0.1 deg: R*tan(a) = {:.4f} mm "
                     "(Zhou: ~0.2 mm over 300 mm)",
                     amp);
        fmt::println("");
    }

    // Phase D.1 of the wafer-grinding plan (TODO/grinding.md): the surface-formation
    // on-ramp. The wheel-spindle vibration -- the axial runout z_b (Eq.14) -- oscillates
    // at frequency f_b and modulates the grain cutting depth, carving a waviness whose
    // spatial wavelength is (relative surface speed) / f_b. Tao observes two patterns
    // (Fig.8): a waviness along the mark direction (WMD, wavelength lambda_m) set by the
    // WHEEL surface speed v_w = n_w*R_w, and a waviness along the circumferential
    // direction (WCD, lambda_c) set by the WAFER surface speed v_s = n_s*r. This
    // reproduces Tao's Table-3 SIMULATED wavelengths directly from the GA kinematics: the
    // two surface speeds are read off kinematic_system3dp::point_velocity of the spinning
    // frames (the same twist field as the "disc surface speeds" test), and lambda =
    // |point_velocity| / f_b.
    //
    // f_b is the MEASURED tilting frequency (6253 Hz = mean of the two measured peaks,
    // Tao sec 4.2.2) -- exactly the value the paper feeds into its surface model for the
    // experimental comparison. Per the Eq.13 analysis (TODO/tao_eq13_derivation.md) this
    // surface-formation half consumes the measured f_b, never the (erroneous) PREDICTED
    // tilt frequency, so it is independent of that controversy. The measured/simulated
    // wavelengths in Table 3 still differ by 5-12% (grain wear, elastoplasticity,
    // grinding heat -- error sources Tao does not model); we gate against the SIMULATED
    // column.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: grinding-mark wavelengths lambda_m / lambda_c (Phase D.1)")
    {
        fmt::println("");
        fmt::println("pga3dp: grinding-mark wavelengths lambda_m / lambda_c (Phase D.1)");
        fmt::println("");

        double const R = 150.0;    // wheel radius == wafer radius [mm]
        double const l3 = 30.0;    // tool-surface distance from the spindle CM [mm]
        double const tw_avg = 0.8; // average wafer thickness [mm]
        double const x_a = l3;     // axial infeed (grain at the chuck plane)
        double const f_b = 6253.0; // measured tilting frequency [Hz] (Tao sec 4.2.2)
        double const r_wcd = 30.0; // radial distance for the WCD measurement [mm]

        struct speeds {
            double v_w; // wheel surface speed (WMD driver) [mm/s]
            double v_s; // wafer surface speed at r_wcd (WCD driver) [mm/s]
        };

        // Build the Fig.-1 two-chain tree (same as the Phase-0 grain-trajectory test)
        // with the chuck spinning at ns and the wheel at nw, then read the two surface
        // speeds straight from the GA twist field via point_velocity (v_s sampled at
        // radius r).
        auto surface_speeds = [&](double ns, double nw, double r) -> speeds {
            kinematic_system3dp sys;
            sys.add_frame(static_frame3dp("chuck_ctr_stat")); // inertial root
            sys.add_frame(static_frame3dp("chuck_ctr_rot"),
                          kin_state3dp{.omega = vec3dp{ns, 0.0, 0.0, 0.0}},
                          sys.index_of("chuck_ctr_stat"));
            sys.add_frame(
                static_frame3dp("wafer_top_avg_rot", vec3dp{tw_avg, 0.0, 0.0, 1.0}),
                kin_state3dp{}, sys.index_of("chuck_ctr_rot"));
            double const off = R / std::sqrt(2.0);
            sys.add_frame(static_frame3dp("spindle_cm_stat", vec3dp{x_a, -off, off, 1.0},
                                          vec3dp{0.0, -pi / 2.0, 0.0, 0.0}),
                          kin_state3dp{}, sys.index_of("chuck_ctr_stat"));
            sys.add_frame(static_frame3dp("tool_top_avg_rot", vec3dp{0.0, 0.0, l3, 1.0}),
                          kin_state3dp{.omega = vec3dp{0.0, 0.0, nw, 0.0}},
                          sys.index_of("spindle_cm_stat"));
            sys.add_frame(
                static_frame3dp("tool_surface_avg_at_R", vec3dp{R, 0.0, 0.0, 1.0}),
                kin_state3dp{}, sys.index_of("tool_top_avg_rot"));

            // wheel-rim grain (root coords): its speed = n_w * R_w
            auto const g = unitize(move3dp(
                O_3dp, sys.get_pos_trafo("tool_surface_avg_at_R", "chuck_ctr_stat")));
            double const v_w =
                to_val(bulk_nrm(sys.point_velocity(g, "tool_surface_avg_at_R")));
            // wafer material point at radius r (root coords): its speed = n_s * r
            auto const P = vec3dp{tw_avg, r, 0.0, 1.0};
            double const v_s =
                to_val(bulk_nrm(sys.point_velocity(P, "wafer_top_avg_rot")));
            return {v_w, v_s};
        };

        // --- WMD: lambda_m = v_w / f_b, vs Tao Table 3 (simulated column) -------------
        // (Table-3 lambda_m rows are at n_s = 265 r/min; lambda_m is independent of n_s.)
        struct wmd_case {
            double nw;      // wheel speed [rpm]
            double lam_sim; // Tao simulated lambda_m [um]
        };
        for (auto const& c : {wmd_case{2250.0, 5652.0}, wmd_case{2650.0, 6656.0}}) {
            auto const s = surface_speeds(rpm2radps(265.0), rpm2radps(c.nw), r_wcd);
            double const lam_m = s.v_w / f_b * 1000.0; // mm -> um
            fmt::println(
                "  n_w = {:>4.0f} rpm: v_w = {:>8.1f} mm/s -> lambda_m = {:>5.0f} "
                "um  (Tao sim {:.0f})",
                c.nw, s.v_w, lam_m, c.lam_sim);
            CHECK(lam_m == doctest::Approx(c.lam_sim).epsilon(2e-3));
        }

        // --- WCD: lambda_c = v_s / f_b, vs Tao Table 3 (simulated column) -------------
        // (Table-3 lambda_c rows are at n_w = 3000 r/min, r = 30 mm; lambda_c is
        // independent of n_w.)
        struct wcd_case {
            double ns;      // wafer speed [rpm]
            double lam_sim; // Tao simulated lambda_c [um]
        };
        for (auto const& c : {wcd_case{265.0, 133.0}, wcd_case{325.0, 163.0}}) {
            auto const s = surface_speeds(rpm2radps(c.ns), rpm2radps(3000.0), r_wcd);
            double const lam_c = s.v_s / f_b * 1000.0; // mm -> um
            fmt::println(
                "  n_s = {:>4.0f} rpm: v_s = {:>8.1f} mm/s -> lambda_c = {:>5.0f} "
                "um  (Tao sim {:.0f})",
                c.ns, s.v_s, lam_c, c.lam_sim);
            CHECK(lam_c == doctest::Approx(c.lam_sim).epsilon(5e-3));
        }

        // Fig.-12 trends (qualitative): lambda_c grows linearly with radial distance r
        // (v_s = n_s*r), while lambda_m does not depend on the wafer speed n_s. One spot
        // check of each: doubling r_wcd doubles lambda_c; changing n_s leaves v_w fixed.
        auto const s30 = surface_speeds(rpm2radps(265.0), rpm2radps(3000.0), 30.0);
        auto const s60 = surface_speeds(rpm2radps(265.0), rpm2radps(3000.0), 60.0);
        CHECK(s60.v_s == doctest::Approx(2.0 * s30.v_s).epsilon(1e-9)); // lambda_c ~ r
        CHECK(
            surface_speeds(rpm2radps(325.0), rpm2radps(2250.0), 30.0).v_w ==
            doctest::Approx(s30.v_w * 2250.0 / 3000.0).epsilon(1e-9)); // v_w indep of n_s
        fmt::println(
            "  Fig.12 trends: lambda_c proportional to r; lambda_m independent of "
            "wafer speed -- confirmed");
        fmt::println("");
    }

    // Phase D.1b of the wafer-grinding plan (TODO/grinding.md): the surface-generation
    // model (Tao Eqs. 23-28). Where D.1 took the wavelength law lambda = v/f_b as given,
    // here the waviness EMERGES from actually carving a surface: many diamond grains
    // (Gaussian size distribution, Eqs. 23-25) sweep the cut at the relative surface
    // speed; a grain centred at arc length s_k cuts to axial depth z_b(t_k) =
    // A_b*sin(2*pi*f_b*t_k) (Eq.26 G^z, the slow feed ramp is sub-femtometre over one
    // window and dropped) and leaves a spherical-cap furrow cap(d) = 0.5*dg -
    // sqrt(0.25*dg^2 - d^2) (Eq.27); the ground surface is the lower envelope of all
    // grains (Eq.28, min). The resulting 1D profile is Tao's Fig.-11 waviness; its
    // measured period must equal the D.1 prediction v/f_b -- for the wheel speed v_w
    // (mark direction, lambda_m) and the wafer speed v_s (circumferential, lambda_c).
    // This validates the surface model itself, not just the formula. (The full 2D
    // topography raster of Figs. 9/10/16/17 -- with the circumferential grain-mark
    // jaggedness -- is a visualization left to a ga_view scene, D.1c.)
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga3dp: simulated wafer topography profiles (Phase D.1b)")
    {
        fmt::println("");
        fmt::println("pga3dp: simulated wafer topography profiles (Phase D.1b)");
        fmt::println("");

        // surface-generation parameters (Tao secs 3.1/3.2). All lengths in micrometres.
        double const f_b = 6253.0;               // vibration (tilting) frequency [Hz]
        double const A_b = 0.1;                  // vibration amplitude [um] (Tao sec 3.2)
        double const d_min = 4.0, d_max = 6.0;   // grain diameter range [um]
        double const Ng = 0.25;                  // grain volume fraction
        double const mu = 0.5 * (d_max + d_min); // Eq.23 mean grain size
        double const sigma = (d_max - d_min) / 6.0;         // Eq.24 std deviation
        double const l_g = mu * std::cbrt(pi / (6.0 * Ng)); // Eq.25 grain spacing [um]

        // relative surface speeds [um/s] (= |point_velocity|, gated in D.1): wheel rim
        // n_w*R_w (mark direction) and wafer n_s*r (circumferential).
        double const Rw = 150.0e3, r_wcd = 30.0e3;         // [um]
        double const v_w = rpm2radps(2250.0) * Rw;         // WMD speed [um/s]
        double const v_s = rpm2radps(265.0) * r_wcd;       // WCD speed [um/s]
        double const lam_m = v_w / f_b, lam_c = v_s / f_b; // predicted wavelengths [um]

        // Carve a 1D axial surface profile h(s) over 5 wavelengths (step lam/250) by the
        // grain model above; sizes ~ N(mu, sigma) (clamped to [d_min,d_max]), positions
        // uniform random with dense overlap so every grid point is covered, fixed seed
        // for a deterministic gate.
        auto carve = [&](double v, double lam) -> std::pair<std::vector<double>, double> {
            double const span = 5.0 * lam, ds = lam / 250.0;
            int const N = int(span / ds) + 1;
            std::vector<double> h(N, A_b * 10.0); // start above any possible cut
            std::mt19937 rng(12345);
            std::normal_distribution<double> grain_d(mu, sigma);
            std::uniform_real_distribution<double> pos(0.0, span);
            int const M = int(20.0 * span / mu); // dense overlap -> full coverage
            for (int k = 0; k < M; ++k) {
                double const sk = pos(rng);
                double const dk = std::clamp(grain_d(rng), d_min, d_max);
                double const zc = A_b * std::sin(2.0 * pi * f_b * sk / v); // z_b at t_k
                int const i0 = std::max(0, int((sk - 0.5 * dk) / ds));
                int const i1 = std::min(N - 1, int((sk + 0.5 * dk) / ds));
                for (int i = i0; i <= i1; ++i) {
                    double const delta = i * ds - sk;
                    double const rad = 0.25 * dk * dk - delta * delta;
                    if (rad < 0.0) continue;
                    double const cap = 0.5 * dk - std::sqrt(rad); // Eq.27 furrow profile
                    h[i] = std::min(h[i], zc + cap);              // Eq.28 lower envelope
                }
            }
            return {h, ds};
        };

        // Measure the waviness wavelength: box-smooth (window << lambda but >> grain, so
        // the sub-grain jaggedness averages out and the z_b waviness survives), subtract
        // the mean, take the mean spacing of upward zero-crossings (cf. measure_freq for
        // time).
        auto wavelength = [&](std::vector<double> const& h, double ds) {
            int const N = int(h.size());
            int const W = std::max(1, N / 50); // smoothing window [points] ~ lam/10
            std::vector<double> g(N);
            double acc = 0.0;
            for (int i = 0; i < N; ++i) {
                acc += h[i];
                if (i >= W) acc -= h[i - W];
                g[i] = acc / std::min(i + 1, W);
            }
            double mean = 0.0;
            for (double x : g)
                mean += x;
            mean /= N;
            std::vector<double> cross;
            for (int i = 1; i < N; ++i) {
                double const a = g[i - 1] - mean, b = g[i] - mean;
                if (a <= 0.0 && b > 0.0) // upward zero-crossing, interpolated
                    cross.push_back((i - 1 + (-a) / (b - a)) * ds);
            }
            double T = 0.0;
            for (size_t i = 1; i < cross.size(); ++i)
                T += cross[i] - cross[i - 1];
            return cross.size() >= 2 ? T / double(cross.size() - 1) : 0.0;
        };

        // compact sparkline of a profile (smoothed, mean-removed, normalized) for the log
        auto sparkline = [&](std::vector<double> const& h) {
            static char const* lv = " .:-=+*#";
            int const N = int(h.size()), cols = 64, W = std::max(1, N / 50);
            std::string out;
            double lo = 1e30, hi = -1e30;
            std::vector<double> s(cols);
            for (int c = 0; c < cols; ++c) {
                int const i = c * (N - 1) / (cols - 1);
                double acc = 0.0;
                int cnt = 0;
                for (int j = std::max(0, i - W); j <= std::min(N - 1, i + W); ++j, ++cnt)
                    acc += h[j];
                s[c] = acc / cnt;
                lo = std::min(lo, s[c]);
                hi = std::max(hi, s[c]);
            }
            for (int c = 0; c < cols; ++c) {
                int const idx = hi > lo ? int((s[c] - lo) / (hi - lo) * 7.0 + 0.5) : 0;
                out += lv[std::clamp(idx, 0, 7)];
            }
            return out;
        };

        // --- WMD (mark direction): carved lambda_m == v_w / f_b -----------------------
        {
            auto const [h, ds] = carve(v_w, lam_m);
            double const meas = wavelength(h, ds);
            fmt::println("  WMD: lambda_m carved = {:>6.0f} um (predicted v_w/f_b = "
                         "{:>6.0f})",
                         meas, lam_m);
            fmt::println("       |{}|  (5 wavelengths, waviness ~{:.2f} um)",
                         sparkline(h), A_b);
            CHECK(meas == doctest::Approx(lam_m).epsilon(0.05));
        }
        // --- WCD (circumferential): carved lambda_c == v_s / f_b ----------------------
        {
            auto const [h, ds] = carve(v_s, lam_c);
            double const meas = wavelength(h, ds);
            fmt::println("  WCD: lambda_c carved = {:>6.1f} um (predicted v_s/f_b = "
                         "{:>6.1f})",
                         meas, lam_c);
            fmt::println("       |{}|  (5 wavelengths, waviness ~{:.2f} um)",
                         sparkline(h), A_b);
            CHECK(meas == doctest::Approx(lam_c).epsilon(0.05));
        }
        fmt::println("  grain spacing l_g = {:.2f} um (Eq.25); sizes ~ N({:.1f}, {:.2f}) "
                     "um",
                     l_g, mu, sigma);
        fmt::println("");
    }

    // Phase D.2a of the wafer-grinding plan (TODO/grinding.md): the contact / grinding-
    // force on-ramp -- FEED-FORWARD only (no back-reaction yet; the loop closes in D.2b).
    // The wheel rim engages the wafer to a depth delta; a grinding-force LAW maps the
    // live contact state to a force which -- once the loop is closed -- becomes a wrench
    // on the spindle. Here we only COMPUTE and report the force at constant feed:
    // "calculate the resulting forces with constant feed velocity first" (user,
    // 2026-06-20).
    //
    // INFRASTRUCTURE PROTOTYPE (promoted into dynamic_system in D.2b): a contact_state
    // POD bundles everything any force law could read; the law is a swappable callback so
    // a student can experiment. Two laws are planned -- (a) a depth law [DEFAULT] and (b)
    // an MRR x specific-energy law. D.2a ships the simplest depth law, F_normal =
    // k*delta, tangential F_tan = mu*F_normal opposing the relative sliding (the
    // chip-thickness form is the first refinement in a later step).
    //
    // CALIBRATION ANCHOR [documented so it is not forgotten -- TODO/grinding.md D.2 #3]:
    // Tao gives no force model and just assumes a constant grinding force
    // F = (F_x, F_y, F_z) = (25, 25, 100) N (Tao sec 2.2, after Zhu [70]). We calibrate
    // the depth-law constant k so the NOMINAL (zero-vibration) engagement delta0
    // reproduces the axial reference F_z = 100 N, and mu so the in-plane reference
    // magnitude sqrt(F_x^2 + F_y^2) = 35.36 N is matched. The in-plane SPLIT (F_x vs F_y)
    // then falls out of the GA contact kinematics (the relative-sliding direction), not a
    // hand-set vector. delta0 is the tunable calibration handle.
    //
    // The relative sliding velocity at the contact is read from the GA twist field
    // (kinematic_system3dp::point_velocity of the wheel grain minus the coincident wafer
    // point) -- the same field as Phase 0 / D.1, so the force direction is GA-derived.
    // The wafer normal is the chuck axis e1; the wheel spins about (root) -e1, so BOTH
    // surface velocities are in-plane -> the tangential force carries the in-plane
    // (F_x,F_y) and the axial F_z comes purely from the depth term.
    //
    // GATES: (1) calibration round-trip -- zero vibration -> axial force == 100 N and
    // in-plane magnitude == 35.36 N; (2) the sliding direction is perpendicular to the
    // wafer normal e1; (3) the rim grain's GA surface speed == n_w*R_w (wheel-dominated);
    // (4) feed-forward modulation -- with z_b(t) = A_b*sin(2*pi*f_b*t) the engagement
    // delta = delta0 - z_b ripples the axial force about 100 N with amplitude k*A_b (the
    // force the spindle will feel once D.2b closes the loop).
    TEST_CASE("pga3dp: grinding contact force, feed-forward (Phase D.2a)")
    {
        fmt::println("");
        fmt::println("pga3dp: grinding contact force, feed-forward (Phase D.2a)");
        fmt::println("");

        // --- machine geometry & spins (Tao Fig.1; identical frame tree to Phase 0)
        // ------
        double const R = 150.0;              // wheel radius == wafer radius [mm]
        double const l3 = 30.0;              // tool surface distance from spindle CM [mm]
        double const tw_avg = 0.8;           // average wafer thickness [mm]
        double const x_a = l3;               // axial infeed -> grain at the chuck plane
        double const ns = rpm2radps(265.0);  // chuck spin [rad/s] (Tao Table 2)
        double const nw = rpm2radps(3000.0); // wheel spin [rad/s]
        double const off = R / std::sqrt(2.0);
        auto const spindle_origin = vec3dp{x_a, -off, off, 1.0};

        kinematic_system3dp sys;
        sys.add_frame(static_frame3dp("chuck_ctr_stat"));
        sys.add_frame(static_frame3dp("chuck_ctr_rot"),
                      kin_state3dp{.omega = vec3dp{ns, 0.0, 0.0, 0.0}},
                      sys.index_of("chuck_ctr_stat"));
        sys.add_frame(static_frame3dp("wafer_top_avg_rot", vec3dp{tw_avg, 0.0, 0.0, 1.0}),
                      kin_state3dp{}, sys.index_of("chuck_ctr_rot"));
        sys.add_frame(static_frame3dp("spindle_cm_stat", spindle_origin,
                                      vec3dp{0.0, -pi / 2.0, 0.0, 0.0}),
                      kin_state3dp{}, sys.index_of("chuck_ctr_stat"));
        sys.add_frame(static_frame3dp("tool_top_avg_rot", vec3dp{0.0, 0.0, l3, 1.0}),
                      kin_state3dp{.omega = vec3dp{0.0, 0.0, nw, 0.0}},
                      sys.index_of("spindle_cm_stat"));
        sys.add_frame(static_frame3dp("tool_surface_avg_at_R", vec3dp{R, 0.0, 0.0, 1.0}),
                      kin_state3dp{}, sys.index_of("tool_top_avg_rot"));

        // --- contact-state POD + swappable force law (infra prototype)
        // ------------------ Promoted into dynamic_system in D.2b; here all quantities
        // are local. SI units (metres / newtons) for the force law, independent of the mm
        // frame tree.
        struct contact_state {
            double delta;  // engagement depth [m] (rim into the nominal wafer surface)
            double v_rel;  // relative sliding speed [m/s] (GA twist field)
            vec3dp t_hat;  // unit sliding direction (world)
            vec3dp n_hat;  // unit wafer-normal / axial direction (world) = +F_z axis
            double v_feed; // nominal feed velocity [m/s]
            double tw;     // current wafer thickness [m]
            double t;      // clock [s]
        };

        // depth law (DEFAULT, option a): normal force k*delta along the wafer normal,
        // plus a tangential friction-like force mu*F_normal opposing the sliding
        // direction. Zero when disengaged (delta <= 0).
        struct depth_law {
            double k, mu;
        };
        auto force_depth = [](contact_state const& c, depth_law const& p) -> vec3dp {
            double const Fn = p.k * std::max(c.delta, 0.0); // axial (normal) magnitude
            double const Ft = p.mu * Fn;                    // tangential magnitude
            return vec3dp{Fn * c.n_hat.x - Ft * c.t_hat.x,
                          Fn * c.n_hat.y - Ft * c.t_hat.y,
                          Fn * c.n_hat.z - Ft * c.t_hat.z, 0.0};
        };

        // --- contact kinematics: relative sliding velocity from the GA twist field
        // ------- The rim grain's velocity has magnitude n_w*R_w for ANY rim point (wheel
        // spin dominates), so the sliding speed/direction does not depend on the engaged
        // phase; delta is the separate modelling quantity below.
        auto const P = unitize(
            move3dp(O_3dp, sys.get_pos_trafo("tool_surface_avg_at_R", "chuck_ctr_stat")));
        auto const v_wheel = sys.point_velocity(P, "tool_surface_avg_at_R");
        auto const v_wafer = sys.point_velocity(P, "wafer_top_avg_rot");
        vec3dp const v_rel_vec = v_wheel - v_wafer; // [mm/s]
        double const v_rel_mm = to_val(bulk_nrm(v_rel_vec));
        vec3dp const t_hat = vec3dp{v_rel_vec.x / v_rel_mm, v_rel_vec.y / v_rel_mm,
                                    v_rel_vec.z / v_rel_mm, 0.0};
        vec3dp const n_hat = vec3dp{1.0, 0.0, 0.0, 0.0}; // chuck axis e1 = wafer normal

        // --- calibration to Tao's constant reference force
        // ------------------------------
        double const delta0 =
            1.0e-6; // nominal engagement [m] (tunable calibration handle)
        double const Fz_ref = 100.0; // Tao axial reference [N]
        double const Fxy_ref = std::sqrt(25.0 * 25.0 + 25.0 * 25.0); // in-plane ref [N]
        depth_law const law{.k = Fz_ref / delta0, .mu = Fxy_ref / Fz_ref};
        double const v_feed = 30.0e-6 / 60.0; // Tao feed 30 um/min -> [m/s]

        contact_state const cs{.delta = delta0,
                               .v_rel = v_rel_mm * 1.0e-3,
                               .t_hat = t_hat,
                               .n_hat = n_hat,
                               .v_feed = v_feed,
                               .tw = tw_avg * 1.0e-3,
                               .t = 0.0};
        vec3dp const F0 = force_depth(cs, law);
        double const F_axial = F0.x; // along e1 = n_hat
        double const F_inplane =
            std::sqrt(F0.y * F0.y + F0.z * F0.z); // F_x,F_y magnitude

        fmt::println("  calibrated depth law: k = {:.3e} N/m  (delta0 = {:.2f} um), "
                     "mu = {:.4f}",
                     law.k, delta0 * 1.0e6, law.mu);
        fmt::println(
            "  nominal force: axial F_z = {:>7.3f} N, in-plane |F_xy| = {:>7.3f} N "
            "(Tao ref 100 / 35.36)",
            F_axial, F_inplane);
        fmt::println("  rim grain GA surface speed |v_wheel| = {:>10.2f} mm/s "
                     "(n_w*R_w = {:>10.2f}); |v_rel| = {:>10.2f}",
                     to_val(bulk_nrm(v_wheel)), nw * R, v_rel_mm);

        // GATE 1: calibration round-trip (axial 100 N, in-plane 35.36 N)
        CHECK(F_axial == doctest::Approx(Fz_ref));
        CHECK(F_inplane == doctest::Approx(Fxy_ref));
        // GATE 2: sliding is in-plane (perpendicular to the wafer normal e1)
        CHECK(std::abs(t_hat.x) < 1.0e-9);
        // GATE 3: rim grain GA surface speed == n_w*R_w (wheel-dominated twist field)
        CHECK(to_val(bulk_nrm(v_wheel)) == doctest::Approx(nw * R));

        // --- GATE 4: feed-forward force modulation by the axial runout z_b
        // --------------- z_b(t) = A_b sin(2 pi f_b t) modulates the engagement delta =
        // delta0 - z_b, so the axial force ripples about 100 N between k*(delta0 -+ A_b).
        // This is the force the spindle will feel once D.2b closes the loop.
        double const f_b = 6253.0; // measured tilting frequency [Hz] (Tao sec 4.2.2)
        double const A_b = 0.1e-6; // vibration amplitude [m] (Tao, 0.1 um)
        double Fmax = -1.0e30, Fmin = 1.0e30;
        int const NS = 400;
        double const Tper = 1.0 / f_b;
        std::string wave;
        static char const* lv = " .:-=+*#";
        for (int i = 0; i <= NS; ++i) {
            double const t = 3.0 * Tper * double(i) / double(NS); // 3 vibration periods
            double const zb = A_b * std::sin(2.0 * pi * f_b * t);
            contact_state c = cs;
            c.delta = delta0 - zb;
            c.t = t;
            double const Fz = force_depth(c, law).x;
            Fmax = std::max(Fmax, Fz);
            Fmin = std::min(Fmin, Fz);
            if (i % (NS / 64) == 0) {
                int const k = int((Fz - 90.0) / 20.0 * 7.0 + 0.5);
                wave += lv[std::clamp(k, 0, 7)];
            }
        }
        fmt::println("  z_b modulation: F_z in [{:>6.2f}, {:>6.2f}] N "
                     "(expect k*(delta0 -+ A_b) = [{:>6.2f}, {:>6.2f}])",
                     Fmin, Fmax, law.k * (delta0 - A_b), law.k * (delta0 + A_b));
        fmt::println("  F_z(t) |{}|  (3 periods of f_b = {:.0f} Hz)", wave, f_b);
        CHECK(Fmax == doctest::Approx(law.k * (delta0 + A_b)));
        CHECK(Fmin == doctest::Approx(law.k * (delta0 - A_b)));
        fmt::println("");
    }

    TEST_CASE("pga3dp: Sommerfeld unbalanced-rotor warm-up (Phase B.1)")
    {

        fmt::println("");
        fmt::println("pga3dp: Sommerfeld unbalanced-rotor warm-up (Phase B.1)");
        fmt::println("");

        // FORCE-ELEMENT WARM-UP (TODO/grinding.md Phase B.1). Prototype the spring/damper
        // + applied-wrench force elements on the Phase-0 spindle's radial DOFs, before
        // promoting them into dynamic_system3dp::tau (Phase A). Reproduces the
        // steady-state analysis of Bisoi et al., "Sommerfeld Effect Characterization in
        // Anisotropic Non-ideal Rotor System", Advances in Rotor Dynamics 2020. The
        // spindle centre of mass is modelled as an unbalanced rotor of total mass (m + M)
        // on an anisotropic foundation in the root radial plane (perpendicular to the
        // chuck axis e1): springs Kx,Ky and dampers Rx,Ry along e2 (x) and e3 (y), with a
        // PRESCRIBED spin omega. The mass eccentricity drives the centrifugal forcing
        // m*e*omega^2, rotating in the radial plane: F_u(t) = m*e*omega^2 * (cos(omega t)
        // e2 + sin(omega t) e3).
        //
        // With prescribed spin the two radial axes DECOUPLE into independent forced
        // oscillators -- the paper's steady-state equation of motion Eq. (2):
        //
        //     (m+M) x'' + Rx x' + Kx x = m*e*omega^2 cos(omega t)        [x along e2]
        //     (m+M) y'' + Ry y' + Ky y = m*e*omega^2 sin(omega t)        [y along e3]
        //
        // GATE: the RK4 steady-state amplitude matches the closed-form Eq. (4)
        //
        //     A(omega) = m*e*omega^2 / sqrt((omega R)^2 + (K - (m+M)omega^2)^2)
        //
        // across an omega sweep. Plotting A vs omega reproduces the paper's Fig. 2 (the
        // forced-response backbone); the resonances are the undamped naturals
        // omega_n = sqrt(K/(m+M)) = 16.51 (x) and 28.57 rad/s (y). (The Sommerfeld speed
        // JUMP itself -- Figs. 3-10 -- needs the two-way non-ideal DC-motor coupling,
        // which is DEFERRED; with prescribed spin we get the resonance backbone, not the
        // jump.)

        // Bisoi et al. (2020) Table 1 parameters
        double const m_rotor = 4.9;          // rotor (disk) mass m [kg]
        double const M_motor = 2.45;         // motor mass M [kg]
        double const mt = m_rotor + M_motor; // total mass m + M = 7.35 kg
        double const e_ecc = 0.008336;       // eccentricity e [m]
        double const me = m_rotor * e_ecc;   // unbalance m*e (forcing m*e*omega^2) [kg m]
        double const Kx = 2000.0;            // foundation stiffness, e2 (x) [N/m]
        double const Ky = 6000.0;            // foundation stiffness, e3 (y) [N/m]
        double const Rx = 5.0;               // foundation damping, e2 (x) [N s/m]
        double const Ry = 10.0;              // foundation damping, e3 (y) [N s/m]

        // undamped natural frequencies, paper: omega1 = 16.51, omega2 = 28.57 rad/s
        double const om1 = std::sqrt(Kx / mt);
        double const om2 = std::sqrt(Ky / mt);
        CHECK(om1 == doctest::Approx(16.51).epsilon(0.005));
        CHECK(om2 == doctest::Approx(28.57).epsilon(0.005));

        // Eq. (4): closed-form steady-state amplitude of one damped SDOF, m*e*omega^2
        // forced
        auto amp_cf = [&](double Om, double K, double R) {
            double const d = K - mt * Om * Om;
            return me * Om * Om / std::sqrt((Om * R) * (Om * R) + d * d);
        };

        // RK4-integrate the decoupled radial dynamics to steady state and return the
        // (x, y) = (e2, e3) amplitudes. State u = [x, y, x', y']; the net radial force is
        // built as a GA direction (vec3dp, w = 0) in the root frame, then Newton's law
        // gives the acceleration -- the same wrench bookkeeping the library force
        // elements will use.
        auto simulate = [&](double Om) -> std::pair<double, double> {
            double const dt = 5.0e-4;  // step [s]; ~125 steps/rev at the top of the sweep
            double const t_end = 30.0; // light damping (Rx/2mt ~ 0.34/s) -> long settle
            double const t_meas = 24.0; // measure the amplitude over the last 6 s

            std::array<value_t, 4> u_mem{0.0, 0.0, 0.0, 0.0};
            std::array<value_t, 8> uh_mem{};
            std::array<value_t, 4> rhs_mem{};
            auto u = std::mdspan<value_t, std::dextents<size_t, 1>>(u_mem.data(), 4);
            auto uh = std::mdspan<value_t, std::dextents<size_t, 2>>(uh_mem.data(), 2, 4);
            auto const rhs =
                std::mdspan<value_t const, std::dextents<size_t, 1>>(rhs_mem.data(), 4);

            double q2_max = std::numeric_limits<double>::lowest();
            double q2_min = std::numeric_limits<double>::max();
            double q3_max = std::numeric_limits<double>::lowest();
            double q3_min = std::numeric_limits<double>::max();

            for (double t = 0.0; t < t_end; t += dt) {
                for (size_t s = 1; s <= 4; ++s) {
                    double const ts = rk4_get_time(t, dt, s - 1);

                    // rotating unbalance force m*e*omega^2 in the root radial plane
                    // (e2,e3)
                    auto const Fu =
                        me * Om * Om *
                        vec3dp{0.0, std::cos(Om * ts), std::sin(Om * ts), 0.0};
                    // anisotropic spring + damper reactions on the CM (the force
                    // elements)
                    auto const Fs = vec3dp{0.0, -Kx * u[0], -Ky * u[1], 0.0};
                    auto const Fd = vec3dp{0.0, -Rx * u[2], -Ry * u[3], 0.0};
                    auto const a = (1.0 / mt) * (Fu + Fs + Fd); // Newton (decoupled axes)

                    rhs_mem[0] = u[2]; // dx/dt = x'
                    rhs_mem[1] = u[3]; // dy/dt = y'
                    rhs_mem[2] = a.y;  // dx'/dt along e2
                    rhs_mem[3] = a.z;  // dy'/dt along e3
                    rk4_step(u, uh, rhs, dt, s);
                }
                if (t >= t_meas) {
                    q2_max = std::max(q2_max, value_t(u[0]));
                    q2_min = std::min(q2_min, value_t(u[0]));
                    q3_max = std::max(q3_max, value_t(u[1]));
                    q3_min = std::min(q3_min, value_t(u[1]));
                }
            }
            return {0.5 * (q2_max - q2_min), 0.5 * (q3_max - q3_min)};
        };

        // sweep across both resonances (Fig. 2 abscissa)
        std::array<double, 8> const sweep{6.0, 12.0, om1, 20.0, 24.0, om2, 34.0, 40.0};
        std::array<double, 8> sim2{}, sim3{};

        fmt::println(
            "   omega      Ax_sim     Ax_cf   err%      Ay_sim     Ay_cf   err%");
        for (size_t i = 0; i < sweep.size(); ++i) {
            double const Om = sweep[i];
            auto const [a2, a3] = simulate(Om);
            double const c2 = amp_cf(Om, Kx, Rx);
            double const c3 = amp_cf(Om, Ky, Ry);
            sim2[i] = a2;
            sim3[i] = a3;
            fmt::println(
                "{:>7.2f}  {:>10.3e} {:>9.3e} {:>5.2f}  {:>10.3e} {:>9.3e} {:>5.2f}", Om,
                a2, c2, 100.0 * std::abs(a2 - c2) / c2, a3, c3,
                100.0 * std::abs(a3 - c3) / c3);

            // GATE: simulated steady-state amplitude matches the closed-form Eq. (4)
            // (both axes; observed error <0.1 %, so 1 % is a comfortable, strong gate)
            CHECK(a2 == doctest::Approx(c2).epsilon(0.01));
            CHECK(a3 == doctest::Approx(c3).epsilon(0.01));
        }

        // resonance peaks emerge from the simulation: the e2 amplitude at om1 exceeds its
        // neighbours, and likewise the e3 amplitude at om2 (sweep indices 2 and 5)
        CHECK(sim2[2] > sim2[1]);
        CHECK(sim2[2] > sim2[3]);
        CHECK(sim3[5] > sim3[4]);
        CHECK(sim3[5] > sim3[6]);

        fmt::println("");
        fmt::println("resonances sqrt(K/mt): e2(x) = {:.2f}, e3(y) = {:.2f} rad/s "
                     "(Bisoi Table-1: 16.51, 28.57)",
                     om1, om2);
        fmt::println("");
    }

    TEST_CASE("pga3dp: line perpendicular to non-intersecting lines (rcmt)")
    {

        fmt::println("");
        fmt::println("pga3dp: line perpendicular to non-intersecting lines (rcmt)");
        fmt::println("");

        auto l1 = x_axis_3dp;
        auto l2 = wdg(vec3dp{0, 1, 0, 1}, vec3dp{0, 1, -1, 1});

        auto l3 = rcmt(l2, rrev(l1));
        CHECK(l3 == y_axis_3dp);

        fmt::println("l1 = {}", l1);
        fmt::println("l2 = {}", l2);
        fmt::println("l3 = {}", l3);
        fmt::println("");

        // use the opportunity to test the sign function
        fmt::println("sign(val) yields +1.0 for val >= 0.0 and -1.0 for val < 0.0).");
        fmt::println("val==0 is always regarded as positive sign:");
        fmt::println("sign(l2.vx) = {} (postive zero)", sign(l2.vx));
        fmt::println("sign(l2.vz) = {} (negative one)", sign(l2.vz));
        fmt::println("sign(l3.vx) = {} (negative zero)", sign(l3.vx));
        fmt::println("sign(l3.vy) = {} (postive one)", sign(l3.vy));
        CHECK(sign(l2.vx) == 1.0);
        CHECK(sign(l2.vz) == -1.0);
        CHECK(sign(l3.vx) == 1.0);
        CHECK(sign(l3.vy) == 1.0);
        fmt::println("");
    }

    TEST_CASE("pga: bivector rate of change")
    {

        fmt::println("");
        fmt::println("pga: bivector rate of change - translation:");

        auto v2dp = vec2dp{1, 2, 0};    // translation vector 2d
        auto v3dp = vec3dp{1, 2, 3, 0}; // translation vector 3d

        auto omega_tra_2dp = att(bulk_dual(v2dp));   // rate of change is a vec (2D)
        auto omega_tra_3dp = att(r_bulk_dual(v3dp)); // rate of change is a bivec (3D)

        CHECK(omega_tra_2dp == -cmpl(wdg(e3_2dp, v2dp)));   // -cmpl(bivec)  = vec
        CHECK(omega_tra_3dp == -r_cmpl(wdg(e4_3dp, v3dp))); // -r_cmpl(bivec) = bivec

        fmt::println("v2dp          = {}", v2dp);
        fmt::println("v3dp          = {}", v3dp);
        fmt::println("omega_tra_2dp = {}", omega_tra_2dp);
        fmt::println("omega_tra_3dp = {}", omega_tra_3dp);
        fmt::println("");
        // fmt::println("att(bulk_dual(v2dp)) = {}", att(bulk_dual(v2dp)));
        // fmt::println("");

        fmt::println("pga: bivector rate of change - rotation:");

        auto omeg2dp = scalar2dp{3};       // omega as angular speed which is a scalar
        auto omeg3dp = vec3dp{1, 2, 3, 0}; // rotation rate omega as vec

        auto omega_rot_2dp = wdg(e3_2dp, omeg2dp); // rate of change is a vec (2D)
        auto omega_rot_3dp = wdg(e4_3dp, omeg3dp); // rate of change is a bivec (3D)

        CHECK(omega_rot_2dp == wdg(e3_2dp, weight_dual(bulk_dual(omeg2dp))));
        CHECK(omega_rot_3dp == wdg(e4_3dp, l_weight_dual(r_bulk_dual(omeg3dp))));

        fmt::println("omeg2dp       = {}", omeg2dp);
        fmt::println("omeg3dp       = {}", omeg3dp);
        fmt::println("omega_rot_2dp = {}", omega_rot_2dp);
        // fmt::println("plr_2dp       = {}", bulk_dual(omeg2dp));
        fmt::println("omega_rot_3dp = {}", omega_rot_3dp);
        // fmt::println("plr_3dp       = {}", r_bulk_dual(omeg3dp));
        fmt::println("");

        fmt::println("pga: bivector rate of change - checking translational speed.");

        auto X1_2dp = vec2dp{1, 0, 1};
        auto X2_2dp = vec2dp{2, 0, 1};
        auto X3_2dp = vec2dp{0, 1, 1};

        CHECK(rcmt(omega_tra_2dp, X1_2dp) == v2dp);
        CHECK(rcmt(omega_tra_2dp, X2_2dp) == v2dp);
        CHECK(rcmt(omega_tra_2dp, X3_2dp) == v2dp);

        auto X1_3dp = vec3dp{1, 0, 0, 1};
        auto X2_3dp = vec3dp{2, 0, 0, 1};
        auto X3_3dp = vec3dp{0, 1, 0, 1};

        CHECK(rcmt(omega_tra_3dp, X1_3dp) == v3dp);
        CHECK(rcmt(omega_tra_3dp, X2_3dp) == v3dp);
        CHECK(rcmt(omega_tra_3dp, X3_3dp) == v3dp);

        fmt::println("pga: bivector rate of change - checking rotational speed.");

        CHECK(bulk_nrm(rcmt(omega_rot_2dp, X1_2dp)) == omeg2dp * bulk_nrm(X1_2dp));
        CHECK(bulk_nrm(rcmt(omega_rot_2dp, X2_2dp)) == omeg2dp * bulk_nrm(X2_2dp));
        CHECK(bulk_nrm(rcmt(omega_rot_2dp, X3_2dp)) == omeg2dp * bulk_nrm(X3_2dp));

        auto L = wdg(O_3dp, omeg3dp); // rotational axis
        CHECK(bulk_nrm(rcmt(omega_rot_3dp, X1_3dp)) ==
              bulk_nrm(omeg3dp) * bulk_nrm(ortho_proj3dp(X1_3dp, L) - X1_3dp));
        CHECK(bulk_nrm(rcmt(omega_rot_3dp, X2_3dp)) ==
              bulk_nrm(omeg3dp) * bulk_nrm(ortho_proj3dp(X2_3dp, L) - X2_3dp));
        CHECK(bulk_nrm(rcmt(omega_rot_3dp, X3_3dp)) ==
              bulk_nrm(omeg3dp) * bulk_nrm(ortho_proj3dp(X3_3dp, L) - X3_3dp));

        fmt::println("");
    }

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
}
