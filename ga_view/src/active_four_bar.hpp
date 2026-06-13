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

// active_four_bar visualizes a closed_loop_system2dp: a planar four-bar linkage, the
// canonical closed kinematic loop. The spanning tree is crank -> coupler (an open chain),
// with the rocker a second branch from the ground; the loop is closed by ONE
// point-coincidence constraint pinning the coupler tip to the rocker tip. The mechanism
// is a Grashof crank-rocker, so the crank rotates fully while the rocker oscillates.
//
// The crank is motor-driven at a constant rate. Each tick advances the crank angle and
// calls closed_loop_system2dp::assemble() to solve the dependent coupler and rocker
// angles (Newton on the constraint Jacobian), so the loop stays closed; the live closure
// error
// ||g|| (~1e-13) shows it does. A point rigidly attached to the coupler traces the
// classic coupler curve.
//
// SPACE toggles pause/resume; R resets to the crank's start angle; T toggles the
// coupler-curve trace.

class active_four_bar : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_four_bar)

  public:

    active_four_bar(Coordsys* cs, w_Coordsys* wcs, four_bar_params const& params,
                    QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:

    void tick();           // advance the crank, re-close the loop
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void toggleTrace();    // T - show/hide the coupler-curve trace
    void viewChanged();

  private:

    QPointF toScreen(vec2dp const& pt) const;
    void build_system(); // (re)build m_cl from m_params at the crank start angle

    // world point of a body-frame point `p` carried by tree frame `name`. Non-const
    // because get_pos_trafo (the tree walk) is a non-const query on the system.
    vec2dp worldPoint(std::string const& name, vec2dp const& p);

    // the joint points in world coordinates (O2, O4 are fixed; A, B, P move)
    vec2dp pivotO2() const; // crank ground pivot
    vec2dp pivotO4() const; // rocker ground pivot
    vec2dp jointA();        // crank tip / coupler base
    vec2dp jointB();        // coupler tip == rocker tip (the closure pin)
    vec2dp couplerPoint();

    void drawBar(QPainter* qp, vec2dp const& a, vec2dp const& b, QColor const& col,
                 double width) const;
    void drawPivot(QPainter* qp, vec2dp const& w, QColor const& col, bool grounded) const;

    Coordsys* cs;
    w_Coordsys* wcs;
    four_bar_params m_params;
    closed_loop_system2dp m_cl;

    // cached frame indices (set in build_system)
    size_t m_CR{0}, m_CO{0}, m_RO{0};
    double m_crank0{0.0}; // crank start angle (for reset)
    double m_crank{0.0};  // current crank angle [rad]
    double m_t{0.0};
    bool m_paused{false};
    bool m_show_trace{true};
    QTimer* m_timer;

    std::deque<vec2dp> m_path; // world-frame trace of the coupler point (coupler curve)

    static constexpr size_t MAX_PATH = 4000;
    static constexpr int UPDATE_MS = 16;
};
