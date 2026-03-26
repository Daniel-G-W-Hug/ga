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
// - PGA2D body-frame formulation: M(t) = M0 ⟇ exp(½ B_b(t)), ODE state = (B_b, Omega_b)
// - M0: pure translation motor placing body origin (= cm) at initial world position
// - cm world position: move2dp(O_b, M) where O_b = (0,0,1) is the body origin
// - Pivot Q_b = TR corner = (hw, hh, 1) in body frame; m_pivot_w is its world image via M0
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

    // Inertia in body frame about pivot Q_b (parallel-axis corrected via get_plate_inertia)
    inertia2dp I;
    inertia2dp I_inv;

    // Initial conditions (stored for reset)
    mvec2dp_u m_M0;          // translation motor: body origin (= cm) to initial world pos
    vec2dp m_pivot_w;        // world position of pivot Q_b = TR corner (updated on drag)
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

    // Cached forces and torque for visualization (updated in calculateRHS)
    vec2dp m_f_cf_w{0.0, 0.0, 0.0}; // centrifugal force at cm in world frame (z=0)
    double m_torque_b{0.0};          // e12-component of F_b in body frame [N·m]

    // Timer for continuous integration (~60 Hz)
    QTimer* m_timer;
    bool m_paused{false};
    static constexpr int UPDATE_INTERVAL_MS = 16;
    static constexpr double DT = 0.016; // time step [s]
    static constexpr size_t SYS_SIZE = 2;
};
