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

#include <vector>

// active_grinding_flatness reproduces Zhou et al. (Precision Eng. 27 (2003) 175) Figs
// 5-7: the GLOBAL FLATNESS error a small wafer-axis misalignment imprints on the ground
// wafer. Dropping the parallel-axis idealization, the wafer axis is tilted by two
// alignment angles -- alpha about X, beta about Y. A single wheel grain (radius r1,
// offset L), sampled in the tilted, spinning wafer frame via PGA motors
// (kinematic_system3dp), traces a 3D cutting path whose height z (wafer-axis component)
// is a CONE for alpha, a CONVEX dome for +beta and a CONCAVE bowl for -beta -- flat only
// when alpha = beta = 0.
//
// Three projected panels (as in the paper): X-Y (the cutting-path rosette), X-Z and Y-Z
// (the profile cross-sections). The in-plane X,Y use a spatial scale; the height Z is
// drawn on a hugely exaggerated scale (z is ~0.2 mm against a 300 mm diameter), exactly
// as Zhou plots it. The path is clipped to the wafer (r2 <= R2).
//
// SPACE pause/resume the trace build-up; R restart it; C steps through the 3x3
// (alpha,beta) grid of Fig. 7 (cone -> dome -> bowl and the combinations).

class active_grinding_flatness : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_grinding_flatness)

  public:

    active_grinding_flatness(Coordsys* cs, w_Coordsys* wcs,
                             grinding_flatness_params const& params,
                             QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:

    void tick();           // grow the revealed trace
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void cycleTilt();      // C - step the (alpha, beta) grid of Fig. 7
    void viewChanged();

  private:

    // one sampled cutting-path point on the wafer, in mm: (x, y) in the wafer plane and
    // the profile height z along the wafer axis. `pen` is true when a line should be
    // drawn from the previous point (false where the grain just re-entered the wafer).
    struct path_pt {
        double x, y, z;
        bool pen;
    };

    void recompute(); // rebuild the cutting path for the current (alpha, beta)

    // panel projections (au), with a spatial scale for the plane coords and an
    // exaggerated scale for z: X-Y rosette, X-Z profile, Y-Z profile.
    QPointF toXY(path_pt const& p) const;
    QPointF toXZ(path_pt const& p) const;
    QPointF toYZ(path_pt const& p) const;

    void drawPath(QPainter* qp,
                  QPointF (active_grinding_flatness::*proj)(path_pt const&) const) const;
    void drawPanel(QPainter* qp, double cx, double cy, double hw, double hh,
                   char const* title) const;

    grinding_flatness_params m_params;
    Coordsys* cs;
    w_Coordsys* wcs;

    QTimer* m_timer{nullptr};
    std::vector<path_pt> m_path; // the full cutting path (mm)
    size_t m_shown{0};           // how much of m_path is revealed (trace build-up)
    bool m_paused{false};
    int m_tilt{4}; // index into the 9-case (alpha,beta) grid; start at (0,0) = flat
};
