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

    // Compute inertia matrix in body frame from plate geometry
    I = get_plate_inertia(m_params.m, m_params.w, m_params.h);
    I_inv = get_inertia_inverse(I);

    // Body-frame corners (plate symmetric about cm = body origin O_b)
    double const hw = m_params.w / 2.0;
    double const hh = m_params.h / 2.0;
    m_corners_b[0] = vec2dp{-hw, -hh, 1.0}; // bottom-left
    m_corners_b[1] = vec2dp{hw, -hh, 1.0};  // bottom-right
    m_corners_b[2] = vec2dp{hw, hh, 1.0};   // top-right = pivot corner
    m_corners_b[3] = vec2dp{-hw, hh, 1.0};  // top-left

    // Store initial pivot position for reset
    pt2d const fix_pt = pivot_pt->scenePos();
    m_initial_pivot = fix_pt;

    // Initial cm world position (rest position: pivot is TR corner)
    m_cm_pos0 = vec2dp{fix_pt.x - hw, fix_pt.y - hh, 1.0};

    // Initial B encodes rotation by phi_init about pivot: B = phi*(px*e23 + py*e31 + e12)
    // For phi_init=0 this is zero (identity motor), plate starts in rest position.
    u_mem[0] = vec2dp{fix_pt.x * m_params.phi_init, fix_pt.y * m_params.phi_init,
                      m_params.phi_init};
    // Initial angular velocity: Omega = omega_init*(px*e23 + py*e31 + e12)
    u_mem[1] = vec2dp{fix_pt.x * m_params.omega_init, fix_pt.y * m_params.omega_init,
                      m_params.omega_init};

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

    vec2dp const B = u[0];
    vec2dp const Omega = u[1];

    // Motor from current bivector B (half-angle convention)
    auto const M = exp(0.5 * B);

    // cm position in world frame: apply motor to initial cm position
    vec2dp const cm_w = move2dp(m_cm_pos0, M);

    // Pivot position in world frame (user-draggable, updated every step)
    pt2d const fix_pt = m_pivot_pt->scenePos();
    vec2dp const pivot_w{fix_pt.x, fix_pt.y, 1.0};

    // Forces in world frame as PGA2D bivectors: F = wdg(application_point, force_vector)
    // Reaction at pivot: +mg upward (constraint force balancing gravity for net-zero
    // translation)
    auto const F_up_w = wdg(pivot_w, vec2dp(0.0, m_params.m * m_params.g, 0.0));
    // Gravity at cm: -mg downward
    auto const F_dn_w = wdg(cm_w, vec2dp(0.0, -m_params.m * m_params.g, 0.0));

    // Transform total force bivector from world to body frame using reverse motor
    auto const F_b = move2dp(F_up_w + F_dn_w, rrev(M));

    // Angular damping torque in body frame (e12 component opposes angular velocity)
    BiVec2dp<value_t> const F_damping_b{0.0, 0.0, -m_params.c * Omega.z};

    // ODE right-hand side:
    //   dB/dt     = Omega
    //   dOmega/dt = I_inv[ F - rcmt(Omega, I[Omega]) ]
    rhs[0] = Omega;
    rhs[1] = compute_omega_dot(I_inv, F_b + F_damping_b, Omega, I);
}

void active_ode_plate::updateState()
{
    auto const M = exp(0.5 * u_mem[0]);
    vec2dp const cm_w = move2dp(m_cm_pos0, M);
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
    auto const M = exp(0.5 * u_mem[0]);
    vec2dp const cm_w = move2dp(m_cm_pos0, M);

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
    pt2d const fix_pt = m_pivot_pt->scenePos();
    QPointF const pivot_screen = toScreen(vec2dp{fix_pt.x, fix_pt.y, 1.0});
    constexpr double PIVOT_RADIUS_PX = 5.0;
    qp->setPen(QPen(Qt::black, 2));
    qp->setBrush(QBrush(Qt::white));
    qp->drawEllipse(pivot_screen, PIVOT_RADIUS_PX, PIVOT_RADIUS_PX);

    // 7. Draw force vectors (scaled for visibility)
    constexpr double FORCE_SCALE = 0.02; // model units per Newton
    double const mg = m_params.m * m_params.g;

    // Gravity at cm: red arrow pointing downward (-y model direction)
    vec2dp const grav_tip_w{cm_w.x, cm_w.y - mg * FORCE_SCALE, 1.0};
    qp->setPen(QPen(Qt::red, 1, Qt::SolidLine));
    qp->setBrush(Qt::NoBrush);
    qp->drawPath(arrowLine(cm_screen, toScreen(grav_tip_w)));
    qp->drawPath(arrowHead(cm_screen, toScreen(grav_tip_w)));

    // Reaction at pivot: dark green arrow pointing upward (+y model direction)
    vec2dp const react_tip_w{fix_pt.x, fix_pt.y + mg * FORCE_SCALE, 1.0};
    qp->setPen(QPen(Qt::darkGreen, 1, Qt::SolidLine));
    qp->drawPath(arrowLine(pivot_screen, toScreen(react_tip_w)));
    qp->drawPath(arrowHead(pivot_screen, toScreen(react_tip_w)));

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
    auto const M = exp(0.5 * u_mem[0]);
    vec2dp const cm_w = move2dp(m_cm_pos0, M);
    pt2d const fix_pt = m_pivot_pt->scenePos();

    double min_x = std::min(cm_w.x, fix_pt.x);
    double max_x = std::max(cm_w.x, fix_pt.x);
    double min_y = std::min(cm_w.y, fix_pt.y);
    double max_y = std::max(cm_w.y, fix_pt.y);

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

    // The world-frame bivectors B and Omega both encode the pivot coordinates:
    //   B     = theta * (px*e23 + py*e31 + e12)  — rotation by theta about (px,py)
    //   Omega = omega * (px*e23 + py*e31 + e12)  — rotation rate omega about (px,py)
    // When the pivot moves, only the scalar parts (B.z = theta, Omega.z = omega)
    // are pivot-independent and must be preserved; the x/y components must be
    // re-expressed for the new pivot coordinates.
    double const theta = u_mem[0].z; // accumulated rotation angle [rad]
    double const omega = u_mem[1].z; // angular velocity [rad/s]

    pt2d const fix_pt = m_pivot_pt->scenePos();
    double const hw = m_params.w / 2.0;
    double const hh = m_params.h / 2.0;

    // New cm rest-position: pivot is the TR corner when motor = identity
    m_cm_pos0 = vec2dp{fix_pt.x - hw, fix_pt.y - hh, 1.0};

    // Re-express B and Omega as rotation about the new pivot
    u_mem[0] = vec2dp{fix_pt.x * theta, fix_pt.y * theta, theta};
    u_mem[1] = vec2dp{fix_pt.x * omega, fix_pt.y * omega, omega};

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

    // Recompute initial cm world position from original pivot
    double const hw = m_params.w / 2.0;
    double const hh = m_params.h / 2.0;
    m_cm_pos0 = vec2dp{m_initial_pivot.x - hw, m_initial_pivot.y - hh, 1.0};

    // Reset integration state (B=phi_init rotation bivector, Omega from omega_init)
    u_mem[0] = vec2dp{m_initial_pivot.x * m_params.phi_init,
                      m_initial_pivot.y * m_params.phi_init, m_params.phi_init};
    u_mem[1] = vec2dp{m_initial_pivot.x * m_params.omega_init,
                      m_initial_pivot.y * m_params.omega_init, m_params.omega_init};

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
