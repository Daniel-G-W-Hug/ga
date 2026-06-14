#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <QColor>
#include <QPen>

#include <cassert> // attribute [[maybe_unused]]
#include <optional>
#include <string>
#include <vector>

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

// use the ga framework
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)

using namespace fmt::literals; // just to make the format literals visible

// required literals for string handling
using namespace std::literals::string_literals;

enum Symbol { plus, cross, circle, square };

using pt2d = hd::ga::Vec2d<value_t>;
using pt2dp = hd::ga::Vec2dp<value_t>;
using bivt2dp = hd::ga::BiVec2dp<value_t>;

// this struct should be used by the user to mark points
struct pt2d_mark {

    QPen pen{QPen(Qt::blue, 2, Qt::SolidLine)};
    Symbol symbol{plus}; // define the marking symbol
    int nsize{4};        // characteristic size of mark symbol in pixels
};

static pt2d_mark const pt2d_mark_default; // for default arguments

// this struct should be used by the user to mark lines
struct ln2d_mark {

    QPen pen{QPen(Qt::black, 2, Qt::SolidLine)};
    QColor area_col{QColor(0, 128, 0, 128)};

    pt2d_mark pm{};
    size_t delta{1}; // 1 shows every point, 2 every second, ...

    bool mark_pts{false};
    bool mark_area{false};
};

static ln2d_mark const ln2d_mark_default; // for default arguments;

// vector as directed line between two points
struct vt2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points later on
    pt2d beg; // beg = (0,0) => position vector, otherwise => free vector
    pt2d end;

    vt2d(pt2d const& p2) : beg{pt2d{0.0, 0.0}}, end{p2} {}     // position vector
    vt2d(pt2d const& p1, pt2d const& p2) : beg{p1}, end{p2} {} // vector from beg to end

    vt2d(vt2d const&) = default;
    vt2d& operator=(vt2d const&) = default;
    vt2d(vt2d&&) = default;
    vt2d& operator=(vt2d&&) = default;
    vt2d() = default;
    ~vt2d() noexcept = default;
};

// this struct should be used by the user to mark vectors
struct vt2d_mark {

    QPen pen{QPen(Qt::black, 2, Qt::SolidLine)};
    std::string display_name = ""s; // default empty string, i.e. don't draw anyting
};

const vt2d_mark vt2d_mark_default; // for default arguments;

// this struct should be used by the user to mark projective lines (=bivectors)
struct bivt2dp_mark {

    QPen pen{QPen(Qt::gray, 1, Qt::SolidLine)};
    std::string display_name = ""s; // default empty string, i.e. don't draw anyting
};

const bivt2dp_mark bivt2dp_mark_default; // for default arguments;


// active vector as directed line between two active points
struct avt2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t beg_idx; // beg = (0,0) => position vector, otherwise => free vector
    size_t end_idx;
};


// active bivector between three active points (one common to both vectors)
struct abivt2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t uend_idx;
    size_t vend_idx;
};

// active reflection using three active points (one common to both vectors)
struct aproj2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t beg_idx;
    size_t uend_idx;
    size_t vend_idx;
};

// active reflection using three active points (one common to both vectors)
struct arefl2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t n1end_idx;
    size_t n2end_idx;
};

// active projective vector as directed line between two active points
struct abivt2dp {

    // reference to points used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t beg_idx;
    size_t end_idx;
};

// active projective reflection using four active points
struct arefl2dp {

    size_t p1_idx; // bvt1 = wdg(p1,p2)
    size_t p2_idx;
    size_t p3_idx; // bvt2 = wdg(p2,p3)
    size_t p4_idx;
};

// active open kinematics chain using an arbitrary number of active points
struct akinematic2dp {

    // first point is base, last point is target
    // open kinematic chain is defined until one point before target
    std::vector<size_t> apt; // index of active point
};

// Physical parameters for spring-mass system
struct spring_params {
    double m = 2.0;  // mass [kg]
    double c = 0.01; // damping constant [Ns/m]
    double k = 40.0; // spring constant [N/m]
    double l0 = 0.1; // initial spring length (w/o force being created by spring)
    double g = 9.81; // gravitational acceleration [m/s²] in -y direction
};

// active ODE spring-mass system using fixation point
struct aode_spring2d {

    size_t fixation_idx;  // index to active_pt2d (fixation point)
    spring_params params; // physical parameters
};

// Physical parameters for rigid plate pendulum (PGA2D)
struct plate_params {
    double m = 1.0;          // mass [kg]
    double w = 2.0;          // width [m] (e1 direction in body frame)
    double h = 2.0;          // height [m] (e2 direction in body frame)
    double g = 9.81;         // gravitational acceleration [m/s²]
    double c = 0.0;          // angular damping coefficient [Nms/rad]
    double phi_init = 0.0;   // initial tilt angle [rad] (rotation about pivot)
    double omega_init = 0.0; // initial angular velocity [rad/s]
};

// active plate pendulum: rigid plate fixed at top-right corner (pivot = active_pt2d)
struct aode_plate_pga2dp {
    size_t fixation_idx; // index to active_pt2d (pivot = top-right corner)
    plate_params params;
};

// Kinematic case selector for body-frame / world-frame transformation demo
enum class kin_case_t { Translation, Rotation, Combined };

// Parameters for the body-frame / world-frame transformation visualization.
// One item per sub-scene (M0=identity upper half, M0 != identity lower half).
struct frame_trafo_params {
    kin_case_t kin_case{kin_case_t::Translation};
    bool m0_is_identity{true}; // informational only (for legend text)

    // M0: initial motor mapping body frame to world frame at t=0.
    // Identity: exp(0.5 * vec2dp{0,0,0}) = ps.
    // Translation to (tx,ty): exp(0.5 * vec2dp{-ty, tx, 0}).
    mvec2dp_u M0;

    // Initial B_b accumulator at t=0 (pre-positions body; B_w always starts at 0).
    // Translation to (tx,ty) from body-frame origin: vec2dp{-ty, tx, 0}.
    vec2dp B_b_init{0.0, 0.0, 0.0};

    // Constant velocity generator in body frame (Omega_b).
    // Translation vx=1, vy=0: vec2dp{0, 1, 0}.
    // Rotation omega about Q_b: omega * Q_b.
    vec2dp Omega_b{0.0, 0.0, 0.0};

    // Animation duration [s] — animation loops at end.
    double duration{6.0};

    // World frame W marker position in world coordinates (fixed reference axes).
    double wx{0.0}, wy{0.0};

    // If false, skip drawing the W marker and axes (use for sub-scenes sharing
    // one world frame — only the first sub-scene should draw it).
    bool draw_world_frame{true};
    bool draw_gen_vectors{false}; // Ω_b / Ω_w overlay arrows (off by default)

    // Text overlay anchor y-coordinate in world coordinates.
    // If left at the sentinel value (−999.0) the overlay defaults to wy + 0.25.
    double text_wy{-999.0};

    // World-frame-drive mode: M(t) = T(O(t)) ⟇ R(cm_omega*t)
    // where O(t) = (cm_ox + cm_vx*t, cm_oy + cm_vy*t).
    // When true, Omega_b is ignored for M(t); B_b/B_w are Euler-integrated
    // observations from the time-varying body/world generators.
    bool world_frame_drive{false};
    double cm_ox{0.0};    // initial CM x (world frame)
    double cm_oy{0.0};    // initial CM y (world frame)
    double cm_vx{0.0};    // CM velocity x (world frame)
    double cm_vy{0.0};    // CM velocity y (world frame)
    double cm_omega{0.0}; // spin rate about CM (rad/s, CCW positive)
};

// Active item: body-frame / world-frame transformation demo (no active points)
struct aframe_trafo {
    frame_trafo_params params;
};

// Parameters for the merry-go-round demo (kinematic_system2dp: rotating platform with
// N turntables mounted at radius r, 120 deg apart, each spinning relative to the
// platform).
// Layout selector for the merry-go-round item.
enum class mgr_layout {
    scene,    // full merry-go-round drawn in world coordinates
    dashboard // mini merry-go-round (corner) + magnitude plots vs platform angle theta
};

struct merry_go_round_params {
    vec2dp centre{0.0, 0.0, 1.0}; // platform centre, in world-frame W coordinates
    double platform_omega{0.6};   // platform angular velocity [rad/s]
    double turntable_omega{2.0};  // turntable spin vs. platform [rad/s]
    double radius{1.6};           // mounting radius of turntables on the platform
    int n_turntables{3};          // number of turntables (>= 1)
    vec2dp marked_point{0.6, 0.0,
                        1.0}; // tracked point on turntable 0 (local coordinates)
    double dt{0.016};         // integration step per tick [s]
    // Canvas anchor of the world frame W: the diagram-grid position at which W's origin
    // is drawn. A pure rigid translation between world coords and canvas coords; the
    // kinematics are unaffected (the tree is purely relative). Default (0,0) draws W at
    // the grid origin (world coords == canvas coords).
    vec2dp world_origin{0.0, 0.0, 1.0};
    mgr_layout layout{mgr_layout::scene}; // scene (full) vs. dashboard (mini + plots)
};

// Active item: merry-go-round demo (no active points)
struct amerry_go_round {
    merry_go_round_params params;
};

// Parameters for the double-pendulum demo (dynamic_system2dp: two revolute-jointed rigid
// plates swinging under gravity -- a chaotic, energy-conserving system). Plate 1 is
// hinged to the world at `pivot` (its top-right corner); plate 2 is hinged to plate 1's
// opposite (bottom-left) corner. Released from (phi1_0, phi2_0) with rates (omega1_0,
// omega2_0).
struct double_pendulum_params {
    vec2dp pivot{0.0, 1.5, 1.0};         // fixed world hinge of plate 1
    double w1{1.4}, h1{1.0}, m1{1.0};    // plate 1: width, height, mass
    double w2{1.4}, h2{1.0}, m2{1.0};    // plate 2: width, height, mass
    double phi1_0{1.9}, phi2_0{-1.1};    // initial joint angles [rad]
    double omega1_0{0.0}, omega2_0{0.0}; // initial joint rates [rad/s]
    double dt{0.004};                    // RK4 sub-step [s]
    int substeps{4};                     // sub-steps per tick (real-time at ~60 fps)
};

// Active item: double-pendulum demo (no active points)
struct adouble_pendulum {
    double_pendulum_params params;
};

// Parameters of the four-bar linkage demo (a closed_loop_system2dp). A Grashof
// crank-rocker: the ground link (crank pivot O2 -> rocker pivot O4) is the longest and
// the crank the shortest, so the crank rotates fully while the rocker oscillates. The
// crank pivot O2 sits at `base`; O4 at base + (ground, 0). `base` centers the linkage in
// the view. The crank is motor-driven at crank_rate.
struct four_bar_params {
    vec2dp base{-1.6, 0.0,
                1.0};    // world position of the crank pivot O2 (centres the figure)
    double ground{3.2};  // O2 -> O4 (fixed) link length
    double crank{0.8};   // O2 -> A   (driver) link length
    double coupler{2.8}; // A  -> B   link length
    double rocker{2.4};  // O4 -> B   link length
    vec2dp coupler_pt{1.4, 0.7, 1.0}; // point on the coupler (origin at A, e1 toward B)
                                      // whose world path is the classic coupler curve
    double crank_rate{1.2};           // driven crank angular rate [rad/s]
    double dt{0.004};                 // time step [s]
    int substeps{4};                  // steps per tick (real-time at ~60 fps)
};

// Active item: four-bar linkage demo (no active points)
struct afour_bar {
    four_bar_params params;
};

// Parameters of the open-vs-closed side-by-side demo. TWO copies of the same 3-link arm
// (shoulder -> elbow -> wrist, each link `link` long), driven by one shared oscillating
// shoulder angle. LEFT (base_open): an open chain -- elbow/wrist held fixed, so the arm
// is rigid and its hand sweeps a free arc. RIGHT (base_closed): a closed loop -- the hand
// is pinned to a fixed world point, so the elbow/wrist are re-solved each frame and the
// arm morphs around the planted hand. The nominal pose places the hand 2*link straight
// above the shoulder, so the shoulder can oscillate +/- amp while the IK stays
// well-conditioned.
struct open_vs_closed_params {
    vec2dp base_open{-1.8, -0.5, 1.0};  // shoulder (ground pivot) of the open arm
    vec2dp base_closed{1.5, -0.5, 1.0}; // shoulder (ground pivot) of the closed arm
    double link{1.0};                   // length of each of the 3 links
    double shoulder0{1.5707963};        // nominal shoulder angle (pi/2, points up)
    double elbow0{-1.0471976};          // nominal elbow angle (-60 deg)
    double wrist0{2.0943951};           // nominal wrist angle (120 deg)
    double amp{0.5};                    // shoulder oscillation amplitude [rad]
    double period{4.0};                 // oscillation period [s]
    double dt{0.004};                   // time step [s]
    int substeps{4};                    // steps per tick
};

// Active item: open-vs-closed side-by-side demo (no active points)
struct aopen_vs_closed {
    open_vs_closed_params params;
};

// Parameters of the planar 5-bar (2-RRR) parallel manipulator demo (a
// closed_loop_system2dp). It is the 2D analogue of a delta / Stewart-Gough robot: a
// moving end-effector positioned IN PARALLEL by two 2-link arms rising from two fixed,
// actuated shoulders. Each arm is shoulder -> upper link -> elbow -> forearm -> tip; both
// tips are pinned to the SAME end-effector by one point-coincidence. 4 revolute joints (2
// shoulders
// + 2 elbows), 2 constraint equations -> 2 DOF. The shoulders SA, SB sit at base_a,
// base_b. The two SHOULDERS are the actuated inputs (the natural parallel-robot drive):
// each frame set_joint(SA/SB, .) + assemble({SA, SB}) solves the two elbows so the tips
// stay coincident (forward kinematics of the closed loop). Phase-shifted shoulder
// sinusoids make the effector trace a Lissajous-like closed curve. Nominal angles point
// the upper arms up-and-outward so the elbows splay and the forearms bend back to the
// central tip -- the well-conditioned interior of the workspace, away from the
// full-extension singularity.
struct planar_delta_params {
    vec2dp base_a{-1.0, -0.6, 1.0}; // world pivot of the left shoulder SA
    vec2dp base_b{1.0, -0.6, 1.0};  // world pivot of the right shoulder SB
    double upper{1.5};              // length of each upper-arm link (shoulder -> elbow)
    double fore{1.8};               // length of each forearm link (elbow -> tip)
    double shoulder_a0{1.9198622};  // nominal left  shoulder angle (110 deg, up-left)
    double shoulder_b0{1.2217305};  // nominal right shoulder angle (70 deg, up-right)
    double amp{0.25};               // shoulder oscillation amplitude [rad] (keeps the
                                    // forearm "knee" ~26 deg clear of the parallel
                                    // singularity across the whole driven sweep)
    double period{5.0};             // oscillation period [s]
    double phase{1.5707963};        // phase lead of shoulder B (pi/2 -> Lissajous figure)
    double dt{0.004};               // time step [s]
    int substeps{4};                // steps per tick (real-time at ~60 fps)
};

// Active item: planar 5-bar parallel manipulator demo (no active points)
struct aplanar_delta {
    planar_delta_params params;
};

// Parameters of the wafer-grinding grain-trajectory demo (Phase 0 of TODO/grinding.md).
// It runs the six-frame kinematic_system3dp of Tao Fig. 1 -- two chains off the
// stationary chuck centre: chuck_ctr_rot -> wafer_top_avg_rot (the wafer, spun at n_s)
// and spindle_cm_stat -> tool_top_avg_rot -> tool_surface_avg_at_R (the wheel, spun at
// n_w) -- and traces a wheel-rim grain IN THE ROTATING WAFER FRAME, projected onto the
// chuck plane e423_3dp (drop e1; e2 horizontal, e3 vertical). That path is the
// grinding-mark pattern of Fig. 1 / Fig. 7. Everything is normalized to R = wheel radius
// = wafer radius = 1, and the spin rates are normalized for a watchable pace -- only the
// ratio n_w / n_s shapes the pattern (here 10, matching n_s = 300 / n_w = 3000 rpm).
// Which frame the grain trajectory is drawn in. wafer_frame: the rotating wafer frame --
// the grain sweeps the grinding-mark rosette (Fig. 1 / Fig. 7). root_frame: the
// stationary global frame -- the wheel is fixed and the grain just circles the wheel rim,
// while the wafer spins beneath. Same computed curve, two viewpoints.
enum class gm_view { wafer_frame, root_frame };

struct grinding_marks_params {
    double R{1.0};      // wheel radius == wafer radius (normalized; the view spans ~2 R)
    double ns{0.8};     // chuck (wafer) spin rate [rad/s, normalized for the animation]
    double ratio{10.0}; // wheel:chuck spin ratio n_w/n_s (the pattern parameter; the
                        // wheel rate is ns * ratio). Cycle it live with the 'C' key --
                        // non-integer ratios (e.g. 10.3) make the marks precess instead
                        // of repeating.
    double dt{0.016};   // sim time advanced per tick [s]
    gm_view view{gm_view::wafer_frame}; // frame the marks are drawn in
};

// Active item: wafer-grinding grain-trajectory demo (no active points)
struct agrinding_marks {
    grinding_marks_params params;
};

// One row in the key-binding table of a legend
struct key_legend_entry {
    std::string key;         // label shown in "Key" column, e.g. "U", "SPACE"
    std::string description; // label shown in "Function" column
};

// Legend box overlaid on the diagram at a fixed viewport position.
// The anchor is the top-left corner of the box.
// x_pct / y_pct are the distances from the left / top edge of the
// drawing area (inside the axes) expressed as a fraction of the area
// width / height. size_pct is the box width as a fraction of area width.
// Text that exceeds the available column width wraps automatically.
struct diagram_legend {
    std::string heading{};                   // bold title row (always shown)
    std::vector<key_legend_entry> entries{}; // optional key-binding rows

    double x_pct{0.02};    // left margin of anchor as fraction of area width
    double y_pct{0.02};    // top  margin of anchor as fraction of area height
    double size_pct{0.30}; // box width as fraction of area width
};

// ----------------------------------------------------------------------------
// convenience alias to make pt2d and ln2d look similar
// convenience alias to make pt2dp and ln2dp look similar
// ----------------------------------------------------------------------------
using ln2d = std::vector<pt2d>;
using cln2dp = std::vector<pt2dp>;
// ----------------------------------------------------------------------------


class Coordsys_model {
  public:

    // add passive point
    [[maybe_unused]] size_t add_pt(pt2d const& p_in,
                                   pt2d_mark const& m = pt2d_mark_default);
    [[maybe_unused]] size_t add_pt(pt2dp const& p_in,
                                   pt2d_mark const& m = pt2d_mark_default);

    // add passive line
    [[maybe_unused]] size_t add_ln(std::vector<pt2d> const& vp_in,
                                   ln2d_mark const& m = ln2d_mark_default);
    [[maybe_unused]] size_t add_ln(std::vector<pt2dp> const& vp_in,
                                   ln2d_mark const& m = ln2d_mark_default);

    // add passive vector
    [[maybe_unused]] size_t add_vt(vt2d const& vt_in,
                                   vt2d_mark const& m = vt2d_mark_default);

    // add passive bivector, i.e. a projective line
    [[maybe_unused]] size_t add_bivtp(bivt2dp const& bivtp_in,
                                      bivt2dp_mark const& m = bivt2dp_mark_default);

    // add active point
    [[maybe_unused]] size_t add_apt(pt2d const& pt_in);

    // add active vector
    [[maybe_unused]] size_t add_avt(avt2d const& vt_in);

    // add active bivector
    [[maybe_unused]] size_t add_abivt(abivt2d const& bivt_in);

    // add active projection
    [[maybe_unused]] size_t add_aproj(aproj2d const& proj_in);

    // add active reflection
    [[maybe_unused]] size_t add_arefl(arefl2d const& refl_in);

    // add active projective bivector (= a projective line)
    [[maybe_unused]] size_t add_abivtp(abivt2dp const& abivtp_in);

    // add active projective reflection
    [[maybe_unused]] size_t add_areflp(arefl2dp const& reflp_in);

    // add active open kinematic chaing
    [[maybe_unused]] size_t add_akinem(akinematic2dp const& akinem_in);

    // add active ODE spring-mass system
    [[maybe_unused]] size_t add_aode(aode_spring2d const& aode_in);

    // add active plate pendulum ODE system
    [[maybe_unused]] size_t add_aode_plate(aode_plate_pga2dp const& aode_plate_in);

    // add body-frame / world-frame transformation demo item
    [[maybe_unused]] size_t add_aframe_trafo(aframe_trafo const& aft_in);

    // add merry-go-round demo item
    [[maybe_unused]] size_t add_merry_go_round(amerry_go_round const& amgr_in);

    // add double-pendulum demo item
    [[maybe_unused]] size_t add_double_pendulum(adouble_pendulum const& adp_in);

    // add four-bar linkage demo item
    [[maybe_unused]] size_t add_four_bar(afour_bar const& afb_in);

    // add open-vs-closed side-by-side demo item
    [[maybe_unused]] size_t add_open_vs_closed(aopen_vs_closed const& aovc_in);

    // add planar 5-bar parallel manipulator demo item
    [[maybe_unused]] size_t add_planar_delta(aplanar_delta const& apld_in);

    [[maybe_unused]] size_t add_grinding_marks(agrinding_marks const& agm_in);

    void set_label(std::string new_label) { m_label = std::move(new_label); };
    std::string label() const { return m_label; }

    void set_legend(diagram_legend leg) { legend = std::move(leg); };

    // reset model to empty state, e.g. for reuse in new model
    void clear();

    // data for points (same index is for same point)
    std::vector<pt2d> pt;
    std::vector<pt2d_mark> pt_mark;

    std::vector<pt2dp> ptp;
    std::vector<pt2d_mark> ptp_mark;

    // data for lines consisting of points (same index is for same line)
    std::vector<ln2d> ln;
    std::vector<ln2d_mark> ln_mark;

    std::vector<cln2dp> clnp;
    std::vector<ln2d_mark> clnp_mark;

    // data for vectors (same index is for same vector)
    std::vector<vt2d> vt;
    std::vector<vt2d_mark> vt_mark;

    // data for projective lines (=bivectors)
    std::vector<bivt2dp> bivtp;
    std::vector<bivt2dp_mark> bivtp_mark;

    // data for active points (same index is for same point)
    std::vector<pt2d> apt;

    // data for active vectors using active points (same index is for same point)
    std::vector<avt2d> avt;

    // data for active bivectors using active points (same index is for same point)
    std::vector<abivt2d> abivt;

    // data for active reflections using active points (same index is for same point)
    std::vector<aproj2d> aproj;

    // data for active reflections using active points (same index is for same point)
    std::vector<arefl2d> arefl;

    // data for active projective bivectors using active points
    std::vector<abivt2dp> abivtp;

    // data for active projective reflections using active points
    std::vector<arefl2dp> areflp;

    // data for active projective kinematics chains using active points
    std::vector<akinematic2dp> akinemp;

    // data for active ODE spring-mass systems using active points
    std::vector<aode_spring2d> aode;

    // data for active plate pendulum ODE systems
    std::vector<aode_plate_pga2dp> aode_plate;

    // data for body-frame / world-frame transformation demo items
    std::vector<aframe_trafo> aft;

    // data for merry-go-round demo items
    std::vector<amerry_go_round> amgr;

    // data for double-pendulum demo items
    std::vector<adouble_pendulum> adp;

    // data for four-bar linkage demo items
    std::vector<afour_bar> afb;

    // data for open-vs-closed side-by-side demo items
    std::vector<aopen_vs_closed> aovc;

    // data for planar 5-bar parallel manipulator demo items
    std::vector<aplanar_delta> apld;

    // data for wafer-grinding grain-trajectory demo items
    std::vector<agrinding_marks> agm;

    // optional legend overlay (key bindings or plain heading)
    std::optional<diagram_legend> legend{};

    // model label, e.g. time stamp description of current Coordsys_model
    std::string m_label{};
};

// ----------------------------------------------------------------------------
// printing support via fmt library
// ----------------------------------------------------------------------------

// formating for user defined types (pt2d)
template <> struct fmt::formatter<pt2d> : fmt::nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(pt2d const& pt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<pt2d>::parse(ParseContext& ctx)
{
    return nested_formatter<double>::parse(ctx);
}

template <typename FormatContext>
auto fmt::formatter<pt2d>::format(pt2d const& pt, FormatContext& ctx) const
{
    auto out = fmt::format_to(ctx.out(), "pt2d(");
    out = fmt::format_to(out, "{}", nested(pt.x));
    out = fmt::format_to(out, ", ");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(pt.y));
    return fmt::format_to(out, ")");
}

// formating for user defined types (pt2dp)
template <> struct fmt::formatter<pt2dp> : fmt::nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(const pt2dp& pt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<pt2dp>::parse(ParseContext& ctx)
{
    return nested_formatter<double>::parse(ctx);
}

template <typename FormatContext>
auto fmt::formatter<pt2dp>::format(pt2dp const& pt, FormatContext& ctx) const
{
    auto out = fmt::format_to(ctx.out(), "pt2dp(");
    out = fmt::format_to(out, "{}", nested(pt.x));
    out = fmt::format_to(out, ", ");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(pt.y));
    out = fmt::format_to(out, ", ");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(pt.z));
    return fmt::format_to(out, ")");
}


// formating for user defined types (vt2d)
template <> struct fmt::formatter<vt2d> : fmt::nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(vt2d const& vt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<vt2d>::parse(ParseContext& ctx)
{
    return nested_formatter<double>::parse(ctx);
}

template <typename FormatContext>
auto fmt::formatter<vt2d>::format(vt2d const& vt, FormatContext& ctx) const
{
    auto out = fmt::format_to(ctx.out(), "vt2d(pt2d(");
    out = fmt::format_to(out, "{}", nested(vt.beg.x));
    out = fmt::format_to(out, ", ");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(vt.beg.y));
    out = fmt::format_to(out, "), pt2d(");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(vt.end.x));
    out = fmt::format_to(out, ", ");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(vt.end.y));
    return fmt::format_to(out, "))");
}

// formating for user defined types (bivt2dp)
template <> struct fmt::formatter<bivt2dp> : fmt::nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(bivt2dp const& bivt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<bivt2dp>::parse(ParseContext& ctx)
{
    return nested_formatter<double>::parse(ctx);
}

template <typename FormatContext>
auto fmt::formatter<bivt2dp>::format(bivt2dp const& bivt, FormatContext& ctx) const
{
    auto out = fmt::format_to(ctx.out(), "bivt2dp(");
    out = fmt::format_to(out, "{}", nested(bivt.x));
    out = fmt::format_to(out, ", ");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(bivt.y));
    out = fmt::format_to(out, ", ");
    ctx.advance_to(out);
    out = fmt::format_to(out, "{}", nested(bivt.z));
    return fmt::format_to(out, ")");
}

// Bsp. für Anwendung
// std::vector<pt2d> vp1{{1.0, 1.0}, {1.5, 2.0}};
// pt2d p{1.0, 2.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));