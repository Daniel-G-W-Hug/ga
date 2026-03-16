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

#include <array>
#include <deque>
#include <vector>

// active_ode_plate models a rigid flat plate acting as a pendulum fixed at its
// top-right corner (the pivot point, user-draggable).
// - PGA2D formulation: motor M = exp(0.5*B), ODE state = (B, Omega), both vec2dp
// - cm world position: move2dp(m_cm_pos0, M) where m_cm_pos0 is the initial cm position
// - Forces: reaction at pivot (+mg upward) + gravity at cm (-mg downward) = pure torque
// - Visualizes: semi-transparent plate polygon, body-frame axes e1/e2, cm trajectory
//   with e2 orientation markers, force vectors at cm and pivot

class active_ode_plate : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_ode_plate)

  public:

    active_ode_plate(Coordsys* cs, w_Coordsys* wcs, active_pt2d* pivot_pt,
                     plate_params const& params, QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void viewChanged();
    void pivotMoved();
    void integrationStep();
    void resetSimulation();
    void togglePause();

  private:

    void calculateRHS();
    void updateState();
    QPointF toScreen(vec2dp const& pt) const;

    Coordsys* cs;
    w_Coordsys* wcs;
    active_pt2d* m_pivot_pt;
    plate_params m_params;

    // Inertia in body frame (computed once from m, w, h)
    inertia2dp I;
    inertia2dp I_inv;

    // Initial conditions (stored for reset)
    vec2dp m_cm_pos0;        // initial cm world position: cm_w(t) = move2dp(m_cm_pos0, M)
    double m_time{0.0};      // simulation elapsed time [s]
    pt2d m_initial_pivot;

    // RK4 integration state: u[0]=B (position bivector), u[1]=Omega (velocity bivector)
    std::vector<vec2dp> u_mem;   // [B, Omega]
    std::vector<vec2dp> uh_mem;  // RK4 helper (2 * SYS_SIZE entries)
    std::vector<vec2dp> rhs_mem; // right-hand side [dB/dt, dOmega/dt]

    // Body-frame plate corners: [BL, BR, TR, TL] symmetric about cm = body origin
    // TR = pivot corner = (hw, hh, 1)
    static constexpr size_t NUM_CORNERS = 4;
    std::array<vec2dp, NUM_CORNERS> m_corners_b;

    // Body-frame axes (constant unit direction vectors, z=0)
    vec2dp m_e1_b{1.0, 0.0, 0.0};
    vec2dp m_e2_b{0.0, 1.0, 0.0};

    // Trajectory of cm with e2 orientation markers
    struct TrajPoint {
        vec2dp cm_w; // cm position in world frame (z=1)
        vec2dp e2_w; // e2 body axis direction in world frame (z=0)
    };
    std::deque<TrajPoint> m_trajectory;
    static constexpr size_t MAX_TRAJ_PTS = 500;

    // Timer for continuous integration (~60 Hz)
    QTimer* m_timer;
    bool m_paused{false};
    static constexpr int UPDATE_INTERVAL_MS = 16;
    static constexpr double DT = 0.016; // time step [s]
    static constexpr size_t SYS_SIZE = 2;
};
