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
#include <functional>

// active_merry_go_round visualizes a kinematic_system2dp: a rotating platform carrying
// N turntables (mounted at radius r, evenly spaced 360/N apart), each spinning relative
// to the platform. One marked point on turntable 0 traces its world path (an epicycloid).
// The system is advanced in place via kinematic_system2dp::step(dt).
//
// Two layouts (selected via merry_go_round_params::layout):
//   - mgr_layout::scene     : the full merry-go-round drawn in world coordinates.
//   - mgr_layout::dashboard : a small merry-go-round in the top-right corner plus three
//     stacked panels plotting |position|, |velocity|, |acceleration| of the marked point
//     against the platform turning angle theta.
//
// SPACE toggles pause/resume; R resets to t=0; C shows/hides the magnitude circles.

class active_merry_go_round : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_merry_go_round)

  public:

    active_merry_go_round(Coordsys* cs, w_Coordsys* wcs,
                          merry_go_round_params const& params,
                          QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void tick();           // advance time, update state
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void toggleCircles();  // C - show/hide momentary-magnitude circles
    void viewChanged();

  private:

    // one recorded sample for the dashboard plots
    struct plot_sample {
        double theta; // platform turning angle [rad] (monotone, cumulative)
        double r;     // |position| of marked point w.r.t. platform centre [world units]
        double v;     // |velocity| of marked point [world units / s]
        double a;     // |acceleration| of marked point [world units / s^2]
    };

    QPointF toScreen(vec2dp const& pt) const;
    void build_system();              // (re)build m_ds from m_params at t=0
    void record_sample();             // append marked-point path + magnitude sample
    std::string tt_name(int i) const; // turntable frame name "T<i>"

    // draw the merry-go-round geometry through an arbitrary world->screen mapping `map`
    // (px/py = pixels per world unit for that mapping). `full` enables the W frame glyph
    // and arrow-tip labels (used full-size; omitted in the dashboard thumbnail).
    void drawSystem(QPainter* qp, std::function<QPointF(vec2dp const&)> const& map,
                    double px, double py, bool full);
    void paintDashboard(QPainter* qp); // mini system + magnitude plots vs theta

    Coordsys* cs;
    w_Coordsys* wcs;
    merry_go_round_params m_params;
    kinematic_system2dp m_ds;

    double m_t{0.0};
    bool m_paused{false};
    bool m_show_circles{true}; // draw momentary-magnitude circles (toggled by 'C')
    QTimer* m_timer;

    std::deque<vec2dp> m_path;      // world-frame trace of the marked point
    std::deque<plot_sample> m_plot; // magnitude history for the dashboard plots

    static constexpr size_t MAX_PATH = 4000;
    static constexpr int UPDATE_MS = 16;
    static constexpr double PLOT_WINDOW =
        2.0 * hd::ga::pi; // plotted theta span [rad]: one platform revolution
    static constexpr double TT_RADIUS = 0.45; // drawn turntable disk radius (world units)
    static constexpr double AXIS_LEN_W = 1.00; // ground-frame axis length (world units)
    static constexpr double AXIS_LEN_P = 0.80; // platform-frame axis length
    static constexpr double AXIS_LEN_T = 0.35; // turntable-frame axis length
    static constexpr double VEL_SCALE =
        0.40; // world units drawn per (unit/s) of velocity
    static constexpr double ACC_SCALE = 0.15; // world units drawn per (unit/s^2) of accel
};
