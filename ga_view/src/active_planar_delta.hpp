#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include "ga/ga_pga.hpp"

#include <QGraphicsItem>
#include <QPainter>
#include <QPointF>
#include <QRectF>
#include <QTimer>
#include <QWidget>

#include <deque>
#include <string>

// active_planar_delta visualizes a closed_loop_system2dp: a planar 5-bar (2-RRR) parallel
// manipulator, the 2D analogue of a delta / Stewart-Gough robot. A moving end-effector is
// positioned IN PARALLEL by two 2-link arms rising from two fixed, actuated shoulders.
// The spanning tree is two independent open chains -- shoulder SA -> elbow EA and
// shoulder SB -> elbow EB; the loop is closed by ONE point-coincidence pinning the two
// forearm tips to the shared end-effector. 4 revolute joints, 2 constraint equations -> 2
// DOF.
//
// The two SHOULDERS are the actuated inputs (the natural parallel-robot drive). Each tick
// sets the two shoulder angles from phase-shifted sinusoids and calls
// closed_loop_system2dp::assemble() to solve the dependent elbow angles (Newton on the
// constraint Jacobian), so the tips stay coincident -- forward kinematics of the closed
// loop. The live closure error ||g|| (~1e-13) shows it does. The end-effector then traces
// a Lissajous-like closed curve.
//
// SPACE toggles pause/resume; R resets to the start configuration; T toggles the
// effector-path trace.

class active_planar_delta : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_planar_delta)

  public:

    active_planar_delta(Coordsys* cs, w_Coordsys* wcs, planar_delta_params const& params,
                        QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:

    void tick();           // advance the shoulders, re-close the loop
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void toggleTrace();    // T - show/hide the effector-path trace
    void viewChanged();

  private:

    QPointF toScreen(vec2dp const& pt) const;
    void build_system(); // (re)build m_cl from m_params at the start configuration

    // current driven shoulder angles at animation time m_t
    double shoulderA() const;
    double shoulderB() const;

    // world point of a body-frame point `p` carried by tree frame `name`. Non-const
    // because get_pos_trafo (the tree walk) is a non-const query on the system.
    vec2dp worldPoint(std::string const& name, vec2dp const& p);

    // the joint points in world coordinates (SA, SB are fixed; elbows + tip move)
    vec2dp pivotSA() const; // left  shoulder ground pivot
    vec2dp pivotSB() const; // right shoulder ground pivot
    vec2dp elbowA();        // left  arm elbow (upper tip / forearm base)
    vec2dp elbowB();        // right arm elbow
    vec2dp effector();      // shared end-effector (forearm tips coincide here)

    void drawBar(QPainter* qp, vec2dp const& a, vec2dp const& b, QColor const& col,
                 double width) const;
    void drawPivot(QPainter* qp, vec2dp const& w, QColor const& col, bool grounded) const;

    Coordsys* cs;
    w_Coordsys* wcs;
    planar_delta_params m_params;
    closed_loop_system2dp m_cl;

    // cached frame indices (set in build_system)
    size_t m_SA{0}, m_EA{0}, m_SB{0}, m_EB{0};
    double m_t{0.0};
    bool m_paused{false};
    bool m_show_trace{true};
    QTimer* m_timer;

    std::deque<vec2dp> m_path; // world-frame trace of the end-effector

    static constexpr size_t MAX_PATH = 4000;
    static constexpr int UPDATE_MS = 16;
};
