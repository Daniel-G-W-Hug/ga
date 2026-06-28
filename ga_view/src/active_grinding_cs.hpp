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

// active_grinding_cs visualizes Phase 0.c of the wafer-grinding plan (TODO/grinding.md):
// the RELATIVE MOTION OF THE COORDINATE SYSTEMS of Tao's wheel-spindle vs. the wafer/
// chuck, to make the spindle's RADIAL / AXIAL / TILTING vibration DOFs legible before the
// dynamic model adds them for real. The spindle placement is a prescribed (illustrative)
// wobble built from PGA3DP motors; the wafer/chuck is fixed.
//
// Two projected views are drawn side by side in the shared aspect-locked Coordsys:
//   - LEFT  "top view"  (look down root -e1, the e423 chuck plane): e2 -> horizontal,
//            e3 -> vertical. RADIAL motion (spindle CM drift in the chuck plane) shows
//            here.
//   - RIGHT "side view" (look along root -e2, the e431 elevation): e3 -> horizontal,
//            e1 -> vertical. AXIAL motion (along the chuck axis e1) and TILTING (the tool
//            axis tipping) show here.
//
// A grey dashed "ghost" marks the ideal (unperturbed) spindle placement, so the live
// deviation is visible. Everything is normalized to R = wheel radius = wafer radius = 1.
// SPACE pause/resume; R reset; C cycles the isolated DOF group (radial / axial / tilt /
// all) so each effect can be seen on its own.

class active_grinding_cs : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_grinding_cs)

  public:

    active_grinding_cs(Coordsys* cs, w_Coordsys* wcs, grinding_cs_params const& params,
                       QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:

    void tick();           // advance time, update state
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void cycleDof();       // C - radial -> axial -> tilt -> all (combined last)
    void toggleIncl();     // T - static equilibrium inclination on / off
    void viewChanged();

  private:

    // which DOF group is currently driven (the rest held at their equilibrium value)
    enum class dof_group { radial, axial, tilt, all };

    // the spindle placement motor at the current time. `ideal` returns the perfect square
    // placement (the grey ghost); otherwise the imperfect equilibrium (optional static
    // inclination) plus the prescribed wobble of the active DOF group.
    mvec3dp_e spindle_motor(bool ideal) const;

    // The two stacked panels of the shared Coordsys share the e3 horizontal axis (so a
    // feature lines up vertically between top and side -- engineering-drawing style).
    // top view: look down -e1 (e423 plane), e3 -> horizontal, e2 -> vertical (depth).
    // side view: look along -e2 (e431 elevation), e3 -> horizontal, e1 -> vertical.
    QPointF toTop(vec3dp const& p) const;
    QPointF toSide(vec3dp const& p) const;

    // draw one spindle rendering (cup-wheel rim + tool axis + origin) into a panel via
    // the supplied projection; `ghost` selects the ideal grey dashed styling.
    void drawSpindle(QPainter* qp,
                     QPointF (active_grinding_cs::*proj)(vec3dp const&) const,
                     bool ghost) const;
    void drawWafer(QPainter* qp, bool side) const;
    void drawPanelFrame(QPainter* qp, double cx, double cy, char const* title) const;

    grinding_cs_params m_params;
    Coordsys* cs;
    w_Coordsys* wcs;

    QTimer* m_timer{nullptr};
    double m_t{0.0};
    bool m_paused{false};
    dof_group m_dof{dof_group::radial}; // start with the easy single-DOF case
    bool m_incl_on{false};              // start from perfect (uninclined) equilibrium
};
