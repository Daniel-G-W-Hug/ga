// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_common.hpp"
#include "active_ode.hpp"

#include "ga/ga_pga.hpp" // includes ga_usr_utilities.hpp with RK4 integration

#include <algorithm>  // std::min, std::max
#include <cmath>      // std::sqrt
#include <fmt/core.h> // fmt::print

using namespace hd::ga;
using namespace hd::ga::pga;

active_ode::active_ode(Coordsys* cs, w_Coordsys* wcs, active_pt2d* fixation_pt,
                       spring_params const& params, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_fixation_pt{fixation_pt}, m_params{params},
    m_paused{false}
{
    // It is essential to NOT use QGraphicsItem::ItemIsMovable, but keep the other flags
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    // Calculate equilibrium position: y_eq = m*g/k + m_params.l0
    double y_eq = m_params.m * m_params.g / m_params.k + m_params.l0;

    // Store initial conditions
    pt2d fix_pt = fixation_pt->scenePos();
    m_initial_fixation = vec2dp(fix_pt.x, fix_pt.y, 1.0);
    vec2dp downward{0.0, y_eq, 0.0};
    m_initial_mass = m_initial_fixation - downward;

    // Initialize velocity to zero
    m_velocity = vec2dp(0.0, 0.0, 0.0);

    // Allocate RK4 memory
    u_mem.resize(FST_ODR_SYS_SIZE);
    uh_mem.resize(2 * FST_ODR_SYS_SIZE);
    rhs_mem.resize(FST_ODR_SYS_SIZE);

    // Set initial state: u[0] = position, u[1] = velocity
    u_mem[0] = m_initial_mass;        // Store position as point (w=1.0)
    u_mem[1] = vec2dp{0.0, 0.0, 0.0}; // Zero velocity as vector (w=0.0)

    // Initialize forces
    m_spring_force = vec2dp(0.0, m_params.m * m_params.g, 0.0);
    m_damping_force = vec2dp(0.0, 0.0, 0.0);
    m_gravity_force = vec2dp(0.0, -m_params.m * m_params.g, 0.0);

    // Connect signals
    connect(m_fixation_pt, &active_pt2d::pointMoved, this, &active_ode::fixationMoved);
    connect(wcs, &w_Coordsys::viewResized, this, &active_ode::viewChanged);

    // Create and start timer
    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_ode::integrationStep);
    m_timer->start(UPDATE_INTERVAL_MS);
}

void active_ode::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                       QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    // clipping area is active area of coordsys
    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));

    qp->save();

    // 1. Draw trajectory trace (light gray, thin line)
    if (m_trajectory.size() > 1) {
        QPen trace_pen(QColor(200, 200, 200, 100), 1, Qt::SolidLine);
        qp->setPen(trace_pen);

        for (size_t i = 1; i < m_trajectory.size(); ++i) {
            QPointF p1(cs->x.au_to_w(m_trajectory[i - 1].x),
                       cs->y.au_to_w(m_trajectory[i - 1].y));
            QPointF p2(cs->x.au_to_w(m_trajectory[i].x),
                       cs->y.au_to_w(m_trajectory[i].y));
            qp->drawLine(p1, p2);
        }
    }

    // 2. Draw mass point (filled circle, dark gray)
    vec2dp mass_model_pos = u_mem[0]; // Current mass position from integration
    QPointF mass_widget_pos(cs->x.au_to_w(mass_model_pos.x),
                            cs->y.au_to_w(mass_model_pos.y));
    constexpr double mass_radius = 6.0; // pixels
    qp->setPen(QPen(Qt::darkGray, 2));
    qp->setBrush(QBrush(Qt::darkGray));
    qp->drawEllipse(mass_widget_pos, mass_radius, mass_radius);

    // 3. Draw spring line (blue, medium thickness)
    QPen spring_pen(Qt::blue, 2, Qt::SolidLine);
    qp->setPen(spring_pen);
    qp->setBrush(Qt::NoBrush);

    pt2d fix_pt = m_fixation_pt->scenePos();
    QPointF fix_pos(cs->x.au_to_w(fix_pt.x), cs->y.au_to_w(fix_pt.y));
    qp->drawLine(fix_pos, mass_widget_pos);

    // 4. Draw equilibrium position marker (orange horizontal line)
    // Calculate equilibrium position: y_eq below fixation point
    double y_eq = m_params.m * m_params.g / m_params.k + m_params.l0;
    vec2dp fixation_model_pos(fix_pt.x, fix_pt.y, 1.0);
    vec2dp downward{0.0, y_eq, 0.0};
    vec2dp equilibrium_pos = fixation_model_pos - downward;

    // Draw horizontal line at equilibrium position along the spring line
    constexpr double eq_marker_width = 0.1; // width in model coordinates
    QPointF eq_left(cs->x.au_to_w(equilibrium_pos.x - eq_marker_width),
                    cs->y.au_to_w(equilibrium_pos.y));
    QPointF eq_right(cs->x.au_to_w(equilibrium_pos.x + eq_marker_width),
                     cs->y.au_to_w(equilibrium_pos.y));

    qp->setPen(QPen(QColor(255, 165, 0), 2, Qt::SolidLine)); // orange
    qp->drawLine(eq_left, eq_right);

    // 5. Draw force vectors (scaled for visibility)
    // Scale factor to make forces visible (adjust as needed)
    constexpr double force_scale = 0.02; // scale forces to coordinate units

    // Spring force (red)
    if (bulk_nrm(m_spring_force) > 1e-3) {
        vec2dp spring_vec = m_spring_force * force_scale;
        QPointF spring_start(cs->x.au_to_w(mass_model_pos.x),
                             cs->y.au_to_w(mass_model_pos.y));
        QPointF spring_end(cs->x.au_to_w(mass_model_pos.x + spring_vec.x),
                           cs->y.au_to_w(mass_model_pos.y + spring_vec.y));

        qp->setPen(QPen(Qt::red, 1, Qt::SolidLine));
        qp->drawPath(arrowLine(spring_start, spring_end));
        qp->drawPath(arrowHead(spring_start, spring_end));
    }

    // Damping force (magenta/purple)
    // Note: Damping force is always along spring direction (by physics calculation)
    if (bulk_nrm(m_damping_force) > 1e-3 && m_damping_force != vec2dp(0.0, 0.0, 0.0)) {
        vec2dp damping_vec = m_damping_force * force_scale;
        QPointF damping_start(cs->x.au_to_w(mass_model_pos.x),
                              cs->y.au_to_w(mass_model_pos.y));
        QPointF damping_end(cs->x.au_to_w(mass_model_pos.x + damping_vec.x),
                            cs->y.au_to_w(mass_model_pos.y + damping_vec.y));

        // qp->setPen(QPen(Qt::magenta, 1, Qt::SolidLine));
        // qp->drawPath(arrowLine(damping_start, damping_end));
        // qp->drawPath(arrowHead(damping_start, damping_end));
    }

    // Gravity force (green)
    if (bulk_nrm(m_gravity_force) > 1e-3) {
        vec2dp gravity_vec = m_gravity_force * force_scale;
        QPointF gravity_start(cs->x.au_to_w(mass_model_pos.x),
                              cs->y.au_to_w(mass_model_pos.y));
        QPointF gravity_end(cs->x.au_to_w(mass_model_pos.x + gravity_vec.x),
                            cs->y.au_to_w(mass_model_pos.y + gravity_vec.y));

        qp->setPen(QPen(Qt::green, 1, Qt::SolidLine));
        qp->drawPath(arrowLine(gravity_start, gravity_end));
        qp->drawPath(arrowHead(gravity_start, gravity_end));
    }

    qp->restore();
}

QRectF active_ode::boundingRect() const
{
    // Bounding box should include fixation, mass, and trajectory
    vec2dp mass_pos = u_mem[0];
    pt2d fix_pt = m_fixation_pt->scenePos();
    double min_x = std::min(fix_pt.x, mass_pos.x);
    double max_x = std::max(fix_pt.x, mass_pos.x);
    double min_y = std::min(fix_pt.y, mass_pos.y);
    double max_y = std::max(fix_pt.y, mass_pos.y);

    // Extend to include trajectory
    for (const auto& pt : m_trajectory) {
        min_x = std::min(min_x, pt.x);
        max_x = std::max(max_x, pt.x);
        min_y = std::min(min_y, pt.y);
        max_y = std::max(max_y, pt.y);
    }

    // Add margin for force vectors
    constexpr double margin = 0.5; // model coordinate units
    min_x -= margin;
    max_x += margin;
    min_y -= margin;
    max_y += margin;

    return QRectF(QPointF(cs->x.au_to_w(min_x), cs->y.au_to_w(max_y)),
                  QPointF(cs->x.au_to_w(max_x), cs->y.au_to_w(min_y)));
}

QPainterPath active_ode::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void active_ode::viewChanged()
{
    // View changed by external influence (e.g., window resize)
    update();
}

void active_ode::fixationMoved()
{
    // Fixation point was moved by user - physics continues with new fixation position
    update();
}

void active_ode::calculateRHS()
{
    // Get current state
    pt2d fix_pt = m_fixation_pt->scenePos();
    vec2dp fixation_pos(fix_pt.x, fix_pt.y, 1.0);
    vec2dp mass_pos = u_mem[0]; // position is stored as point (w=1.0) in u[0]
    vec2dp velocity = u_mem[1]; // velocity is stored as vector (w=0.0) in u[1]

    // Calculate spring displacement from fixation to mass minus initial spring length l0
    vec2dp displacement_vec = mass_pos - fixation_pos;
    vec2dp displacement_unit_vec = bulk_normalize(displacement_vec);

    double elongated_spring_length =
        std::max(to_val(bulk_nrm(displacement_vec)) - m_params.l0, 0.0);
    vec2dp elongation_vec = displacement_unit_vec * elongated_spring_length;

    // Spring force: F = -k * elongation_vec (Hooke's law)
    // At equilibrium (displacement = -y_eq where y_eq = m*g/k + m_params.l0):
    //   F_spring = -k * (0, -y_eq) = (0, +k*y_eq) = (0, +m*g) (upward)
    //   F_gravity = (0, -m*g) (downward)
    //   F_total = 0 ✓
    m_spring_force = -m_params.k * elongation_vec;

    // Damping force: only damps velocity component along spring direction
    // F_damping = -c * (v · ŝ) * ŝ, where ŝ is unit vector along spring
    if (elongated_spring_length > 1e-10) {
        double velocity_along_spring =
            to_val(dot(velocity, displacement_unit_vec)); // v · ŝ
        m_damping_force = -m_params.c * velocity_along_spring * displacement_unit_vec;
    }
    else {
        m_damping_force = vec2dp(0.0, 0.0, 0.0); // avoid division by zero
    }

    m_gravity_force = vec2dp(0.0, -m_params.m * m_params.g, 0.0);

    // Debug output when damping force is significant
    // double vel_mag = std::sqrt(velocity.x * velocity.x + velocity.y * velocity.y);
    // double damp_mag = std::sqrt(m_damping_force.x * m_damping_force.x +
    //                             m_damping_force.y * m_damping_force.y);
    // if (vel_mag > 0.01 || damp_mag > 0.01) {
    //     fmt::print("velocity: ({:.6f}, {:.6f}), mag={:.6f}\n", velocity.x, velocity.y,
    //     vel_mag); fmt::print("m_damping_force: ({:.6f}, {:.6f}), mag={:.6f}\n",
    //     m_damping_force.x, m_damping_force.y, damp_mag);
    // }

    vec2dp F_total = m_spring_force + m_damping_force + m_gravity_force;
    vec2dp acceleration = F_total / m_params.m;

    // Set right-hand side for ODE system
    // u[0]' = velocity
    // u[1]' = acceleration
    rhs_mem[0] = velocity;
    rhs_mem[1] = acceleration;
}

void active_ode::updateMassPosition()
{
    // Update velocity
    m_velocity = u_mem[1];

    // Add to trajectory
    vec2dp new_mass_pos = u_mem[0];
    m_trajectory.push_back(new_mass_pos);
    if (m_trajectory.size() > MAX_TRAJECTORY_POINTS) {
        m_trajectory.pop_front();
    }
}

void active_ode::integrationStep()
{
    if (m_paused) {
        return; // Skip integration if paused
    }

    // Create mdspan views for RK4 integration
    auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), FST_ODR_SYS_SIZE);
    auto uh = mdspan<vec2dp, dextents<size_t, 2>>(uh_mem.data(), 2, FST_ODR_SYS_SIZE);
    auto rhs =
        mdspan<vec2dp const, dextents<size_t, 1>>(rhs_mem.data(), FST_ODR_SYS_SIZE);

    // Perform RK4 integration (4 sub-steps)
    for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
        calculateRHS();
        rk4_step(u, uh, rhs, DT, rk_step);
    }

    // Update mass position and trajectory
    updateMassPosition();

    // Recalculate forces with final integrated state for visualization
    calculateRHS();

    update(); // Trigger repaint
}

void active_ode::resetSimulation()
{
    // Restore initial fixation position
    m_fixation_pt->setScenePos(pt2d(m_initial_fixation.x, m_initial_fixation.y));

    // Reset integration state
    u_mem[0] = m_initial_mass;        // Store position as point (w=1.0)
    u_mem[1] = vec2dp{0.0, 0.0, 0.0}; // Zero velocity as vector (w=0.0)

    // Reset velocity
    m_velocity = vec2dp{0.0, 0.0, 0.0};

    // Reset forces
    m_spring_force = vec2dp(0.0, m_params.m * m_params.g, 0.0);
    m_damping_force = vec2dp(0.0, 0.0, 0.0);
    m_gravity_force = vec2dp(0.0, -m_params.m * m_params.g, 0.0);

    // Clear trajectory
    m_trajectory.clear();

    update();
}

void active_ode::togglePause()
{
    m_paused = !m_paused;
    update(); // Repaint to show paused state (could add visual indicator later)
}
