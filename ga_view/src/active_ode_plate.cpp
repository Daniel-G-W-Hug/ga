// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_common.hpp"
#include "active_ode_plate.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <algorithm> // std::min, std::max
#include <cmath>
#include <fmt/core.h>

using namespace hd::ga;
using namespace hd::ga::pga;

active_ode_plate::active_ode_plate(Coordsys* cs, w_Coordsys* wcs, active_pt2d* pivot_pt,
                                   plate_params const& params, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_pivot_pt{pivot_pt}, m_params{params},
    u_mem(SYS_SIZE), uh_mem(2 * SYS_SIZE), rhs_mem(SYS_SIZE)
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    // Body-frame corners (plate symmetric about cm = body origin O_b)
    double const hw = m_params.w / 2.0;
    double const hh = m_params.h / 2.0;
    m_corners_b[0] = vec2dp{-hw, -hh, 1.0}; // bottom-left
    m_corners_b[1] = vec2dp{hw, -hh, 1.0};  // bottom-right
    m_corners_b[2] = vec2dp{hw, hh, 1.0};   // top-right = pivot corner Q_b
    m_corners_b[3] = vec2dp{-hw, hh, 1.0};  // top-left

    // Pivot Q_b = TR corner in body frame; inertia about Q_b via parallel-axis correction
    vec2dp const Q_b{hw, hh, 1.0};
    I = get_plate_inertia(m_params.m, m_params.w, m_params.h, Q_b);
    I_inv = get_inertia_inverse(I);

    // Store initial pivot position for reset
    pt2d const fix_pt = pivot_pt->scenePos();
    m_initial_pivot = fix_pt;

    // M0: pure translation placing body origin (= cm) at initial world cm position
    // cm_pos0 = pivot_world - (hw, hh) offset
    vec2dp const cm_pos0{fix_pt.x - hw, fix_pt.y - hh, 1.0};
    m_M0 = exp(0.5 * vec2dp{-cm_pos0.y, cm_pos0.x, 0.0});

    // World position of pivot Q_b via M0 (stays fixed during simulation)
    m_pivot_w = move2dp(Q_b, m_M0);

    // Initial body-frame state: B_b = phi_init * Q_b, Omega_b = omega_init * Q_b
    u_mem[0] = m_params.phi_init * Q_b;
    u_mem[1] = m_params.omega_init * Q_b;

    // Connect signals
    connect(pivot_pt, &active_pt2d::pointMoved, this, &active_ode_plate::pivotMoved);
    connect(wcs, &w_Coordsys::viewResized, this, &active_ode_plate::viewChanged);

    // Create and start integration timer
    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_ode_plate::integrationStep);
    m_timer->start(UPDATE_INTERVAL_MS);
}

QPointF active_ode_plate::toScreen(vec2dp const& pt) const
{
    return QPointF(cs->x.au_to_w(pt.x), cs->y.au_to_w(pt.y));
}

void active_ode_plate::calculateRHS()
{
    auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), SYS_SIZE);
    auto rhs = mdspan<vec2dp, dextents<size_t, 1>>(rhs_mem.data(), SYS_SIZE);

    vec2dp const B = u[0];     // B_b = phi * Q_b  (body-frame rotation bivector)
    vec2dp const Omega = u[1]; // Omega_b = omega * Q_b  (body-frame angular velocity)

    // Current motor: M(t) = M0 ⟇ exp(½ B_b(t))  [body-frame formulation]
    auto const M = rgpr(m_M0, exp(0.5 * B));

    // cm world position: body origin O_b = (0,0,1) mapped to world via M
    vec2dp const O_b{0.0, 0.0, 1.0};
    auto const cm_w = move2dp(O_b, M);

    // Force at cm: gravity + centrifugal (pivot reaction = -f_cm)
    vec2dp const g_vec{0.0, -m_params.g, 0.0};
    value_t const omega_val = Omega.z;
    auto const a_cf = (cm_w - m_pivot_w) * (omega_val * omega_val);
    auto const f_cm = value_t(m_params.m) * (a_cf + g_vec);

    // Cache centrifugal component for visualization (static gravity always = (0,-mg,0))
    m_f_cf_w = value_t(m_params.m) * a_cf;

    // Force couple: f_cm at cm, -f_cm at pivot (equal and opposite → pure torque)
    auto const F_cm_w = wdg(cm_w, f_cm);
    auto const F_piv_w = wdg(m_pivot_w, -f_cm);
    auto const F_b = move2dp(F_cm_w + F_piv_w, rrev(M));

    // Cache net torque (e12-component of body-frame force bivector) for visualization
    m_torque_b = F_b.z;

    // Angular damping torque in body frame (opposes angular velocity)
    BiVec2dp<value_t> const F_damping_b{0.0, 0.0, -m_params.c * Omega.z};

    // Enforce pivot constraint explicitly: extract scalar alpha from dOmega.z and
    // project back onto Q_b so that Omega_b = omega*Q_b is preserved exactly,
    // preventing numerical drift of off-constraint components (e31/e32) during
    // long-running RK4 integration.
    // In principle correct without projection: rhs[1] = compute_omega_dot(I_inv, F_b +
    // F_damping_b, Omega, I);
    double const hw = m_params.w / 2.0;
    double const hh = m_params.h / 2.0;
    vec2dp const Q_b{hw, hh, 1.0};
    auto const dOmega = compute_omega_dot(I_inv, F_b + F_damping_b, Omega, I);
    double const alpha = dOmega.z; // = tau_eff / I_pivot[2,2]

    rhs[0] = Omega;       // dB_b/dt = Omega_b
    rhs[1] = alpha * Q_b; // dOmega_b/dt = alpha * Q_b  (pivot constraint)
}

void active_ode_plate::updateState()
{
    auto const M = rgpr(m_M0, exp(0.5 * u_mem[0]));
    vec2dp const O_b{0.0, 0.0, 1.0};
    vec2dp const cm_w = move2dp(O_b, M);
    vec2dp const e2_w = move2dp(m_e2_b, M); // e2 body direction in world frame

    m_trajectory.push_back({cm_w, e2_w});
    if (m_trajectory.size() > MAX_TRAJ_PTS) {
        m_trajectory.pop_front();
    }

    m_time += DT;
}

void active_ode_plate::integrationStep()
{
    if (m_paused) return;

    auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), SYS_SIZE);
    auto uh = mdspan<vec2dp, dextents<size_t, 2>>(uh_mem.data(), 2, SYS_SIZE);
    auto rhs = mdspan<vec2dp const, dextents<size_t, 1>>(rhs_mem.data(), SYS_SIZE);

    // RK4 integration (4 sub-steps)
    for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
        calculateRHS();
        rk4_step(u, uh, rhs, DT, rk_step);
    }

    updateState();
    calculateRHS(); // refresh final state for paint()

    update(); // trigger repaint
}

void active_ode_plate::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                             QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    // Clip to coordsys active area
    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));

    qp->save();

    // Compute current motor and world-frame quantities
    auto const M = rgpr(m_M0, exp(0.5 * u_mem[0]));
    vec2dp const O_b{0.0, 0.0, 1.0};
    vec2dp const cm_w = move2dp(O_b, M);

    // e1/e2 body-frame directions in world frame (direction vecs, z=0 → not translated)
    vec2dp const e1_w = move2dp(m_e1_b, M);
    vec2dp const e2_w = move2dp(m_e2_b, M);

    // Transform plate corners to world frame → screen
    // Corners are stored as body-frame offsets from cm; rotate only (z=0 direction
    // vectors)
    std::array<QPointF, NUM_CORNERS> corners_screen;
    for (size_t i = 0; i < NUM_CORNERS; ++i) {
        vec2dp const offset_b{m_corners_b[i].x, m_corners_b[i].y, 0.0}; // direction (z=0)
        vec2dp const offset_w = move2dp(offset_b, M);                   // rotate only
        vec2dp const corner_w{cm_w.x + offset_w.x, cm_w.y + offset_w.y, 1.0};
        corners_screen[i] = toScreen(corner_w);
    }

    QPointF const cm_screen = toScreen(cm_w);

    // 1. Draw cm trajectory trace (light gray, thin)
    if (m_trajectory.size() > 1) {
        qp->setPen(QPen(QColor(180, 180, 180, 120), 1, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_trajectory.size(); ++i) {
            qp->drawLine(toScreen(m_trajectory[i - 1].cm_w),
                         toScreen(m_trajectory[i].cm_w));
        }
    }

    // // 2. Draw e2 orientation markers along trajectory (every ORIENT_STEP points)
    // constexpr size_t ORIENT_STEP = 10;
    // constexpr double ORIENT_LEN = 0.08; // model coordinate units
    // qp->setPen(QPen(QColor(80, 80, 180, 160), 1, Qt::SolidLine));
    // for (size_t i = 0; i < m_trajectory.size(); i += ORIENT_STEP) {
    //     TrajPoint const& tp = m_trajectory[i];
    //     vec2dp const tip_w{tp.cm_w.x + tp.e2_w.x * ORIENT_LEN,
    //                        tp.cm_w.y + tp.e2_w.y * ORIENT_LEN, 1.0};
    //     qp->drawLine(toScreen(tp.cm_w), toScreen(tip_w));
    // }

    // 3. Draw semi-transparent filled plate polygon
    QPolygonF plate_poly;
    for (QPointF const& p : corners_screen) {
        plate_poly << p;
    }
    qp->setPen(QPen(Qt::darkBlue, 2, Qt::SolidLine));
    qp->setBrush(QBrush(QColor(70, 130, 180, 80))); // steel blue, semi-transparent
    qp->drawPolygon(plate_poly);

    // 4. Draw body-frame axes at cm
    double const axis_len = std::min(m_params.w, m_params.h) * 0.35;

    // e1 axis: red arrow
    vec2dp const e1_tip_w{cm_w.x + e1_w.x * axis_len, cm_w.y + e1_w.y * axis_len, 1.0};
    qp->setPen(QPen(Qt::red, 2, Qt::SolidLine));
    qp->setBrush(Qt::NoBrush);
    qp->drawPath(arrowLine(cm_screen, toScreen(e1_tip_w)));
    qp->drawPath(arrowHead(cm_screen, toScreen(e1_tip_w)));

    // e2 axis: dark green arrow
    vec2dp const e2_tip_w{cm_w.x + e2_w.x * axis_len, cm_w.y + e2_w.y * axis_len, 1.0};
    qp->setPen(QPen(Qt::darkGreen, 2, Qt::SolidLine));
    qp->drawPath(arrowLine(cm_screen, toScreen(e2_tip_w)));
    qp->drawPath(arrowHead(cm_screen, toScreen(e2_tip_w)));

    // 5. Draw cm marker (black cross)
    constexpr double CM_MARK_PX = 5.0;
    qp->setPen(QPen(Qt::black, 2));
    qp->drawLine(cm_screen + QPointF(-CM_MARK_PX, 0.0),
                 cm_screen + QPointF(CM_MARK_PX, 0.0));
    qp->drawLine(cm_screen + QPointF(0.0, -CM_MARK_PX),
                 cm_screen + QPointF(0.0, CM_MARK_PX));

    // 6. Draw pivot marker (white-filled circle)
    QPointF const pivot_screen = toScreen(m_pivot_w);
    constexpr double PIVOT_RADIUS_PX = 5.0;
    qp->setPen(QPen(Qt::black, 2));
    qp->setBrush(QBrush(Qt::white));
    qp->drawEllipse(pivot_screen, PIVOT_RADIUS_PX, PIVOT_RADIUS_PX);

    // 7. Draw forces at cm and pivot: static (gravity/reaction) and dynamic (centrifugal)
    constexpr double FORCE_SCALE = 0.06; // model units per Newton
    qp->setBrush(Qt::NoBrush);

    QColor const col_grav{180, 0, 0};       // dark red:    gravity at cm (static)
    QColor const col_react{0, 140, 60};     // forest green: static reaction at pivot
    QColor const col_cf{230, 115, 0};       // orange:      centrifugal at cm (dynamic)
    QColor const col_cf_react{60, 60, 200}; // medium blue: centrifugal reaction at pivot

    double const mg = m_params.m * m_params.g;

    // Static: gravity downward at cm
    {
        vec2dp const tip{cm_w.x, cm_w.y - mg * FORCE_SCALE, 1.0};
        qp->setPen(QPen(col_grav, 2, Qt::SolidLine));
        qp->drawPath(arrowLine(cm_screen, toScreen(tip)));
        qp->drawPath(arrowHead(cm_screen, toScreen(tip)));
    }
    // Static: reaction upward at pivot
    {
        vec2dp const tip{m_pivot_w.x, m_pivot_w.y + mg * FORCE_SCALE, 1.0};
        qp->setPen(QPen(col_react, 2, Qt::SolidLine));
        qp->drawPath(arrowLine(pivot_screen, toScreen(tip)));
        qp->drawPath(arrowHead(pivot_screen, toScreen(tip)));
    }

    // Dynamic: centrifugal force at cm — threshold in screen pixels to avoid
    // atan2 instability when the tip is within numerical noise of the base.
    {
        QPointF const cf_tip = toScreen(vec2dp{cm_w.x + m_f_cf_w.x * FORCE_SCALE,
                                               cm_w.y + m_f_cf_w.y * FORCE_SCALE, 1.0});
        QPointF const d = cf_tip - cm_screen;
        if (d.x() * d.x() + d.y() * d.y() >= ARROWSIZE * ARROWSIZE) {
            qp->setPen(QPen(col_cf, 2, Qt::SolidLine));
            qp->drawPath(arrowLine(cm_screen, cf_tip));
            qp->drawPath(arrowHead(cm_screen, cf_tip));
        }
    }
    // Dynamic: centrifugal reaction at pivot (equal and opposite)
    {
        QPointF const cfr_tip =
            toScreen(vec2dp{m_pivot_w.x - m_f_cf_w.x * FORCE_SCALE,
                            m_pivot_w.y - m_f_cf_w.y * FORCE_SCALE, 1.0});
        QPointF const d = cfr_tip - pivot_screen;
        if (d.x() * d.x() + d.y() * d.y() >= ARROWSIZE * ARROWSIZE) {
            qp->setPen(QPen(col_cf_react, 2, Qt::SolidLine));
            qp->drawPath(arrowLine(pivot_screen, cfr_tip));
            qp->drawPath(arrowHead(pivot_screen, cfr_tip));
        }
    }

    // 7b. Torque visualization: semi-transparent circle around pivot, area ∝ |torque|,
    //     four tangential arrows at 3/6/9/12 o'clock indicating rotation direction.
    if (std::abs(m_torque_b) > 1e-4) {
        constexpr double TORQUE_AREA_SCALE = 0.05; // circle area per N·m [model units²]
        double const r_model =
            std::sqrt(std::abs(m_torque_b) * TORQUE_AREA_SCALE / std::acos(-1.0));
        double const r_px =
            std::abs(cs->x.au_to_w(m_pivot_w.x + r_model) - cs->x.au_to_w(m_pivot_w.x));

        // Color: teal for CCW (positive torque), orange for CW (negative torque)
        QColor const torque_color =
            (m_torque_b >= 0.0) ? QColor(0, 170, 130) : QColor(210, 90, 0);
        QColor torque_fill = torque_color;
        torque_fill.setAlpha(35);

        qp->setPen(QPen(torque_color, 1.5, Qt::SolidLine));
        qp->setBrush(QBrush(torque_fill));
        qp->drawEllipse(pivot_screen, r_px, r_px);

        // Four tangential arrows: position offsets and CCW tangent directions
        // [dx, dy, tx_ccw, ty_ccw]: offset from pivot, tangent direction for CCW
        struct ClockArrow {
            double dx, dy, tx, ty;
        };
        std::array<ClockArrow, 4> const clocks = {{
            {r_model, 0.0, 0.0, 1.0},   // 3 o'clock:  right, CCW tangent = up
            {0.0, r_model, -1.0, 0.0},  // 12 o'clock: top,   CCW tangent = left
            {-r_model, 0.0, 0.0, -1.0}, // 9 o'clock:  left,  CCW tangent = down
            {0.0, -r_model, 1.0, 0.0},  // 6 o'clock:  bottom,CCW tangent = right
        }};
        double const sgn = (m_torque_b >= 0.0) ? 1.0 : -1.0;
        double const arrow_half = 0.25 * r_model;
        qp->setPen(QPen(torque_color, 2, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (ClockArrow const& c : clocks) {
            double const px = m_pivot_w.x + c.dx;
            double const py = m_pivot_w.y + c.dy;
            double const tx = sgn * c.tx;
            double const ty = sgn * c.ty;
            QPointF const arw_start =
                toScreen(vec2dp{px - tx * arrow_half, py - ty * arrow_half, 1.0});
            QPointF const arw_end =
                toScreen(vec2dp{px + tx * arrow_half, py + ty * arrow_half, 1.0});
            qp->drawPath(arrowHead(arw_start, arw_end));
        }
    }

    // 8. Draw simulation time near pivot
    qp->setPen(QPen(Qt::darkBlue, 1));
    QFont mono_font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    mono_font.setPointSize(9);
    qp->setFont(mono_font);
    qp->drawText(pivot_screen + QPointF(8.0, -8.0),
                 QString("t = %1 s").arg(m_time, 0, 'f', 2));

    qp->restore();
}

QRectF active_ode_plate::boundingRect() const
{
    auto const M = rgpr(m_M0, exp(0.5 * u_mem[0]));
    vec2dp const O_b{0.0, 0.0, 1.0};
    vec2dp const cm_w = move2dp(O_b, M);

    double min_x = std::min(cm_w.x, m_pivot_w.x);
    double max_x = std::max(cm_w.x, m_pivot_w.x);
    double min_y = std::min(cm_w.y, m_pivot_w.y);
    double max_y = std::max(cm_w.y, m_pivot_w.y);

    // Include all plate corners
    for (size_t i = 0; i < NUM_CORNERS; ++i) {
        vec2dp const offset_b{m_corners_b[i].x, m_corners_b[i].y, 0.0};
        vec2dp const offset_w = move2dp(offset_b, M);
        vec2dp const corner_w{cm_w.x + offset_w.x, cm_w.y + offset_w.y, 1.0};
        min_x = std::min(min_x, corner_w.x);
        max_x = std::max(max_x, corner_w.x);
        min_y = std::min(min_y, corner_w.y);
        max_y = std::max(max_y, corner_w.y);
    }

    // Include trajectory
    for (TrajPoint const& tp : m_trajectory) {
        min_x = std::min(min_x, tp.cm_w.x);
        max_x = std::max(max_x, tp.cm_w.x);
        min_y = std::min(min_y, tp.cm_w.y);
        max_y = std::max(max_y, tp.cm_w.y);
    }

    constexpr double MARGIN = 1.0;
    min_x -= MARGIN;
    max_x += MARGIN;
    min_y -= MARGIN;
    max_y += MARGIN;

    return QRectF(QPointF(cs->x.au_to_w(min_x), cs->y.au_to_w(max_y)),
                  QPointF(cs->x.au_to_w(max_x), cs->y.au_to_w(min_y)));
}

QPainterPath active_ode_plate::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void active_ode_plate::viewChanged() { update(); }

void active_ode_plate::pivotMoved()
{
    // Stop integration timer so all state changes are applied atomically before
    // the next integration step starts.
    m_timer->stop();

    // Preserve accumulated rotation angle and angular velocity (body-frame scalars).
    // B_b = phi * Q_b and Omega_b = omega * Q_b: Q_b is constant in the body frame,
    // so only the scalar multiples change — and they don't change here (pivot move
    // only relocates the body in the world, not its current rotation state).
    double const theta = u_mem[0].z; // current rotation angle [rad]  (B_b.z = phi)
    double const omega = u_mem[1].z; // current angular velocity [rad/s]

    pt2d const fix_pt = m_pivot_pt->scenePos();
    double const hw = m_params.w / 2.0;
    double const hh = m_params.h / 2.0;

    // Recompute M0 for the new pivot position:
    // cm rest position = pivot - (hw, hh) body-frame offset
    vec2dp const cm_pos0_new{fix_pt.x - hw, fix_pt.y - hh, 1.0};
    m_M0 = exp(0.5 * vec2dp{-cm_pos0_new.y, cm_pos0_new.x, 0.0});

    // Update stored world pivot (= move2dp(Q_b, m_M0) = (fix_pt.x, fix_pt.y, 1))
    vec2dp const Q_b{hw, hh, 1.0};
    m_pivot_w = move2dp(Q_b, m_M0);

    // Body-frame state unchanged in form: B_b ∝ Q_b, Omega_b ∝ Q_b
    u_mem[0] = theta * Q_b;
    u_mem[1] = omega * Q_b;

    // Refresh RHS so the first integration sub-step and paint() use consistent values
    calculateRHS();

    m_timer->start(UPDATE_INTERVAL_MS);
    prepareGeometryChange();
    update();
}

void active_ode_plate::resetSimulation()
{
    // Restore pivot to initial position
    m_pivot_pt->setScenePos(m_initial_pivot);

    double const hw = m_params.w / 2.0;
    double const hh = m_params.h / 2.0;

    // Recompute M0 from initial pivot position
    vec2dp const cm_pos0_reset{m_initial_pivot.x - hw, m_initial_pivot.y - hh, 1.0};
    m_M0 = exp(0.5 * vec2dp{-cm_pos0_reset.y, cm_pos0_reset.x, 0.0});

    // Restore world pivot
    vec2dp const Q_b{hw, hh, 1.0};
    m_pivot_w = move2dp(Q_b, m_M0);

    // Reset body-frame integration state
    u_mem[0] = m_params.phi_init * Q_b;
    u_mem[1] = m_params.omega_init * Q_b;

    // Reset time and trajectory
    m_time = 0.0;
    m_trajectory.clear();

    update();
}

void active_ode_plate::togglePause()
{
    m_paused = !m_paused;
    update();
}
