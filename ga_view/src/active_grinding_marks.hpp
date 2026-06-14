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

// active_grinding_marks visualizes Phase 0 of the wafer-grinding plan (TODO/grinding.md):
// the GRAIN TRAJECTORY a wheel-rim grain traces in the ROTATING WAFER frame -- the
// grinding-mark pattern of Tao Fig. 1 / Fig. 7. It runs the six-frame kinematic_system3dp
// of the geometry test (two chains off the stationary chuck centre) and projects each
// world point onto the chuck plane e423_3dp -- a top view looking down the chuck axis
// -e1, with e2 horizontal and e3 vertical (matching the test's printed (y,z) column).
//
// In the wafer frame the wafer is fixed (the unit circle at the origin) and the wheel,
// fixed in the root frame, appears to orbit -- so the rim grain sweeps out the rosette.
// Everything is normalized to R = wheel radius = wafer radius = 1; only the ratio
// n_w / n_s shapes the pattern. SPACE pause/resume; R reset; T toggle the trace.

class active_grinding_marks : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_grinding_marks)

  public:

    active_grinding_marks(Coordsys* cs, w_Coordsys* wcs,
                          grinding_marks_params const& params,
                          QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void tick();           // advance time, update state
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void toggleTrace();    // T - show/hide the grain-mark trace
    void cycleRatio();     // C - cycle the wheel:chuck spin ratio (10.0 / 10.3 / 9.7)
    void viewChanged();

  private:

    // one recorded grain position, stored in WAFER-LOCAL coordinates plus whether the
    // grain was on the wafer (a mark is only cut while it is). Display transforms it by
    // the wafer->display-frame motor: identity in the wafer view, the live wafer rotation
    // in the global view -- so the SAME carved marks are seen from either frame.
    struct mark_sample {
        vec3dp p_wafer; // grain position in the wafer frame
        bool contact;   // grain within the wafer disk (r <= R) -> a mark is cut
    };

    // project a 3dp world point onto the chuck plane e423_3dp (drop the e1 component) and
    // map to canvas: e2 -> horizontal, e3 -> vertical
    QPointF toScreen(vec3dp const& p) const;

    char const* display_frame() const; // frame the scene is drawn in (per m_params.view)
    double nw() const;                 // wheel spin rate = ns * current ratio
    void build_system();               // (re)build m_sys from m_params at t=0
    void advance_to(double t);         // set the chuck + wheel poses for sim time t
    void record_sample(); // append the rim grain's wafer-frame position + flag

    // draw a circle of world radius r (in the e423 projection) centred at world point c
    void drawCircle(QPainter* qp, vec3dp const& c, double r, QColor const& edge,
                    QColor const& fill);

    Coordsys* cs;
    w_Coordsys* wcs;
    grinding_marks_params m_params;
    kinematic_system3dp m_sys;

    double m_t{0.0};
    double m_ratio{10.0}; // current wheel:chuck ratio (init from params; cycled by 'C')
    bool m_paused{false};
    bool m_show_trace{true};
    QTimer* m_timer;

    std::deque<mark_sample> m_path; // grain marks in wafer-local coords (+ contact flag)

    static constexpr size_t MAX_PATH = 12000;
    static constexpr int UPDATE_MS = 16;
};
