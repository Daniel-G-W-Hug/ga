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

// active_frame_trafo visualizes the body-frame / world-frame transformation
// framework (PGA2D) for three kinematic cases:
//   Translation, Rotation, Combined
//
// Each scene contains two items: one for M0 = identity (upper half) and one
// for M0 != identity (lower half).  Both animate in parallel so the user can
// compare B_b and B_w side-by-side.
//
// State variables (all in PGA2D):
//   B_b  - accumulated position generator in body  frame (integrated)
//   B_w  - accumulated position generator in world frame (integrated)
//   M(t) = M0 ⟇ exp(½ B_b(t))   [body-frame motor formula]
//
// The animation loops every params.duration seconds.
// SPACE toggles pause/resume; R resets to t=0.

class active_frame_trafo : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_frame_trafo)

  public:

    active_frame_trafo(Coordsys* cs, w_Coordsys* wcs,
                       frame_trafo_params const& params,
                       QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void tick();             // advance time, update state
    void togglePause();      // SPACE
    void resetAnimation();   // R
    void viewChanged();

  private:

    QPointF toScreen(vec2dp const& pt) const;

    void drawWorldFrame(QPainter* qp, vec2dp const& O_w,
                        vec2dp const& e1_w, vec2dp const& e2_w);
    void drawTrajectory(QPainter* qp);
    void drawBody(QPainter* qp, vec2dp const& O_w,
                  vec2dp const& e1_w, vec2dp const& e2_w);
    void drawHomeMarker(QPainter* qp); // M0 reference: body position when B_b=0
    void drawGeneratorVectors(QPainter* qp, vec2dp const& O_w,
                              vec2dp const& e1_w, vec2dp const& e2_w);
    void drawTextOverlay(QPainter* qp);

    Coordsys* cs;
    w_Coordsys* wcs;
    frame_trafo_params m_params;

    double m_t{0.0};
    bool   m_paused{false};
    QTimer* m_timer;

    vec2dp    m_B_b{0.0, 0.0, 0.0};
    vec2dp    m_B_w{0.0, 0.0, 0.0};
    mvec2dp_u m_M;

    struct TrajPt { vec2dp pos; };
    std::deque<TrajPt> m_trajectory;

    static constexpr size_t MAX_TRAJ   = 350;
    static constexpr int    UPDATE_MS  = 16;
    static constexpr double DT         = 0.016; // [s] per tick

    // Visual scales (world-coordinate units)
    static constexpr double AXIS_LEN    = 0.50; // frame-axis arrow length
    static constexpr double POLE_OFFSET = 0.20; // pole base offset along e1_b from body origin
    static constexpr double POLE_LEN    = 0.45; // flagpole height along e2_b
    static constexpr double FLAG_W      = 0.32; // flag width  (e1_b direction)
    static constexpr double FLAG_H      = 0.24; // flag height (–e2_b direction from pole tip)

    // Velocity-vector overlay scales
    static constexpr double VEC_SCALE   = 0.25; // model-units per unit of Ω magnitude
    static constexpr double VEC_MAX     = 0.90; // max arrow length (model units)
    static constexpr double ARC_R       = 0.14; // rotation-arc radius (model units)
};
