#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_pt2d.hpp"
#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include "ga/ga_pga.hpp"

#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QMarginsF>
#include <QPainter>
#include <QPointF>
#include <QRectF>
#include <QTimer>
#include <QWidget>

#include <deque>
#include <vector>

// Note: spring_params is defined in coordsys_model.hpp (included above)

// active_ode models a damped vertical spring-mass system with continuous ODE integration
// - Fixation point can be moved by mouse (active_pt2d)
// - Mass point is driven by physics (passive visualization)
// - Uses RK4 time integration at 60 Hz
// - Visualizes: spring line, trajectory trace, force vectors
// - Interactive reset and pause/resume via keyboard

class active_ode : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_ode)

  public:

    active_ode(Coordsys* cs, w_Coordsys* wcs, active_pt2d* fixation_pt,
               spring_params const& params, QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void viewChanged();
    void fixationMoved();
    void integrationStep(); // Called by timer at 60 Hz
    void resetSimulation(); // Reset to initial state
    void togglePause();     // Pause/resume integration

  private:

    void calculateRHS();
    void updateMassPosition();

    Coordsys* cs;
    w_Coordsys* wcs;

    active_pt2d* m_fixation_pt; // fixation point (user-controllable)

    spring_params m_params;

    // State variables (in model coordinates, not widget pixels)
    vec2dp m_initial_fixation;  // For reset (point: w=1.0)
    vec2dp m_initial_mass;      // For reset (point: w=1.0)
    vec2dp m_velocity;          // Current velocity [m/s] (vector: w=0.0)

    // Force vectors for visualization (in Newton [N])
    vec2dp m_spring_force;   // (vector: w=0.0)
    vec2dp m_damping_force;  // (vector: w=0.0)
    vec2dp m_gravity_force;  // (vector: w=0.0)

    // RK4 integration state (FST_ODR_SYS_SIZE = 2: [position, velocity])
    std::vector<vec2dp> u_mem;   // [position(w=1.0), velocity(w=0.0)]
    std::vector<vec2dp> uh_mem;  // helper for integration
    std::vector<vec2dp> rhs_mem; // right-hand side values

    // Trajectory trace
    std::deque<vec2dp> m_trajectory; // Limited-length trajectory (points: w=1.0)
    static constexpr size_t MAX_TRAJECTORY_POINTS = 250;

    // Timer for continuous integration
    QTimer* m_timer;
    bool m_paused;
    static constexpr int UPDATE_INTERVAL_MS = 16; // ~60 Hz (1000/60 ≈ 16.67)
    static constexpr double DT = 0.016;           // Time step [s]

    // System size for 1st order ODE system (position + velocity)
    static constexpr size_t FST_ODR_SYS_SIZE = 2;
};
