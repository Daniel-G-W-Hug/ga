// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <algorithm> // std::min, std::max, std::clamp
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
}
