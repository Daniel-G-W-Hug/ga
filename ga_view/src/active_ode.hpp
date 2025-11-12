#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_pt2d.hpp"
#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include "ga/ga_ega.hpp"

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
    pt2d m_initial_fixation;  // For reset
    pt2d m_initial_mass;      // For reset
    hd::ga::vec2d m_velocity; // Current velocity [m/s]

    // Force vectors for visualization (in Newton [N])
    hd::ga::vec2d m_spring_force;
    hd::ga::vec2d m_damping_force;
    hd::ga::vec2d m_gravity_force;

    // RK4 integration state (FST_ODR_SYS_SIZE = 2: [position, velocity])
    std::vector<hd::ga::vec2d> u_mem;   // [position, velocity]
    std::vector<hd::ga::vec2d> uh_mem;  // helper for integration
    std::vector<hd::ga::vec2d> rhs_mem; // right-hand side values

    // Trajectory trace
    std::deque<pt2d> m_trajectory; // Limited-length trajectory
    static constexpr size_t MAX_TRAJECTORY_POINTS = 250;

    // Timer for continuous integration
    QTimer* m_timer;
    bool m_paused;
    static constexpr int UPDATE_INTERVAL_MS = 16; // ~60 Hz (1000/60 ≈ 16.67)
    static constexpr double DT = 0.016;           // Time step [s]

    // System size for 1st order ODE system (position + velocity)
    static constexpr size_t FST_ODR_SYS_SIZE = 2;
};
