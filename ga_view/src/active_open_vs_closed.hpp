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

// active_open_vs_closed contrasts an OPEN chain with a CLOSED loop, side by side, driven
// by one shared oscillating shoulder angle. Both panels are the SAME 3-link arm (shoulder
// -> elbow -> wrist):
//
//   LEFT  (open chain):  a closed_loop_system2dp with NO loop constraint -- the elbow and
//                        wrist are held at their nominal angles, so the arm is rigid and
//                        its hand sweeps a free arc as the shoulder oscillates.
//   RIGHT (closed loop): the SAME arm whose hand is pinned to a fixed world point by one
//                        point-coincidence constraint. Driving the same shoulder now
//                        forces the elbow and wrist (the dependent joints) to be
//                        re-solved each frame so the hand stays planted -- the arm morphs
//                        around it.
//
// Using closed_loop_system2dp for both keeps the code uniform: the open arm's assemble()
// is a no-op (no constraints). The headline visual for "tree vs. tree + closure": one
// input, two behaviours.
//
// SPACE toggles pause/resume; R resets; T toggles the free-hand trace.

class active_open_vs_closed : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_open_vs_closed)

  public:

    active_open_vs_closed(Coordsys* cs, w_Coordsys* wcs,
                          open_vs_closed_params const& params,
                          QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:

    void tick();
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void toggleTrace();    // T - show/hide the free-hand trace
    void viewChanged();

  private:

    QPointF toScreen(vec2dp const& pt) const;
    void build_system();
    void build_arm(closed_loop_system2dp& cl,
                   vec2dp const& base); // 3-link arm at nominal

    // world point of a body-frame point `p` carried by frame `name` of system `cl`
    vec2dp worldPt(closed_loop_system2dp& cl, std::string const& name, vec2dp const& p);
    vec2dp handPt(closed_loop_system2dp& cl); // tip of the wrist link (the hand)

    void drawBar(QPainter* qp, vec2dp const& a, vec2dp const& b, QColor const& col,
                 double width) const;
    void drawPivot(QPainter* qp, vec2dp const& w, QColor const& col, bool grounded) const;
    void drawArm(QPainter* qp, closed_loop_system2dp& cl);

    Coordsys* cs;
    w_Coordsys* wcs;
    open_vs_closed_params m_params;

    closed_loop_system2dp m_open;   // left: free hand (no constraint)
    closed_loop_system2dp m_closed; // right: pinned hand (one coincidence)
    size_t m_shoulder{1};           // shoulder joint index (frame "L1") in both systems
    vec2dp m_pin{0.0, 0.0, 1.0};    // fixed world pin point of the closed arm's hand

    double m_t{0.0};
    bool m_paused{false};
    bool m_show_trace{true};
    QTimer* m_timer;

    std::deque<vec2dp> m_open_hand_path; // world trace of the OPEN arm's free hand

    static constexpr size_t MAX_PATH = 4000;
    static constexpr int UPDATE_MS = 16;
};
