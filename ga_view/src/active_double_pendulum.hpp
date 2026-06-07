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

// active_double_pendulum visualizes a dynamic_system2dp: two revolute-jointed rigid
// plates swinging under gravity -- a chaotic, energy-conserving double pendulum. Plate 1
// is hinged to the world at a fixed pivot (its top-right corner); plate 2 is hinged to
// plate 1's opposite (bottom-left) corner. The marked tip of plate 2 traces its (chaotic)
// world path. The system is advanced in place via dynamic_system2dp::step(dt) (RK4,
// sub-stepped).
//
// SPACE toggles pause/resume; R resets to t=0.

class active_double_pendulum : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_double_pendulum)

  public:

    active_double_pendulum(Coordsys* cs, w_Coordsys* wcs,
                           double_pendulum_params const& params,
                           QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void tick();           // advance time, update state
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void toggleTrace();    // T - show/hide the cg / path trace
    void toggleICR();      // I - show/hide the instantaneous centre of rotation
    void toggleForces();   // F - show/hide the force / moment overlay
    void viewChanged();

  private:

    QPointF toScreen(vec2dp const& pt) const;
    void build_system(); // (re)build m_ds from m_params at t=0
    void recordICR();    // append plate 2's ICR (or a NaN gap) to m_icr2_path

    // draw one plate (frame `name`, size w x h) as a filled polygon + body-frame axes
    void drawPlate(QPainter* qp, std::string const& name, double w, double h,
                   QColor const& fill, QColor const& edge, char const* label);

    // draw the force overlay (gravity at each cm + hinge reaction forces). Toggled by
    // 'F'.
    void drawForces(QPainter* qp);

    // draw the instantaneous-centre-of-rotation overlay (light-grey trace of plate 2's
    // ICR
    // + a crosshair at each plate's current ICR). Toggled by 'I'.
    void drawICR(QPainter* qp);

    Coordsys* cs;
    w_Coordsys* wcs;
    double_pendulum_params m_params;
    dynamic_system2dp m_ds;

    double m_t{0.0};
    double m_E0{0.0}; // initial total energy (for the live conservation readout)
    bool m_paused{false};
    bool m_show_trace{true};   // cg / path trace          (toggled by 'T', default on)
    bool m_show_icr{false};    // instantaneous centre of rotation (toggled by 'I')
    bool m_show_forces{false}; // gravity + hinge reactions (toggled by 'F')
    QTimer* m_timer;

    std::deque<vec2dp> m_path;      // world-frame trace of plate 2's swinging tip
    std::deque<vec2dp> m_icr2_path; // trace of plate 2's instantaneous centre of rotation
                                    // (gap markers carry x = NaN where the ICR is at
                                    // infinity, i.e. when plate 2 is momentarily not
                                    // rotating)

    static constexpr size_t MAX_PATH = 2000;
    static constexpr int UPDATE_MS = 16;
    static constexpr double FORCE_SCALE = 0.02; // world units drawn per unit of force
    static constexpr double ICR_MIN_OMEGA =
        0.01; // below this |omega| the ICR is "at inf"
    static constexpr double ICR_MAX_REACH = 1000.0; // drop ICR samples farther than this
};
