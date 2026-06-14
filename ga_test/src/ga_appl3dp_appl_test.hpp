// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <algorithm> // std::min, std::max
#include <cmath>     // std::sin, std::cos, std::atan
#include <limits>    // std::numeric_limits

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
