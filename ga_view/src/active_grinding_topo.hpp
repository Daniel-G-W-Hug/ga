#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include "ga/ga_pga.hpp"

#include <QGraphicsItem>
#include <QImage>
#include <QPainter>
#include <QPointF>
#include <QRectF>
#include <QTimer>
#include <QWidget>

// active_grinding_topo visualizes Phase D.1c of the wafer-grinding plan
// (TODO/grinding.md): the ground-wafer surface TOPOGRAPHY of Tao 2022, Figs. 9/10. The
// wheel-spindle vibration (axial runout z_b, Eq.14) at frequency f_b modulates the grain
// cutting depth, carving a waviness whose spatial wavelength is (relative surface speed)
// / f_b: a fine ripple along the grinding-mark direction (WMD, lambda_m = v_w / f_b, the
// wheel rim speed) and a coarse modulation around the circumference (WCD, lambda_c = v_s
// / f_b, the wafer speed at radius r). Both wavelengths are computed live from the GA
// kinematics -- kinematic_system3dp::point_velocity of the spinning Fig.-1 frames -- so
// the rendered height field is built straight from the library's twist field, and its
// wavelengths match the paper (validated to <1 um in the Phase-D.1 tests).
//
// Two patch presets (toggle with C) mirror the paper's two figures: a near-square window
// (Fig. 9, the circumferential waviness lambda_c) and a wide strip (Fig. 10, the
// mark-direction waviness lambda_m, height exaggerated). SPACE pause/resume; R reset; the
// pattern slowly scrolls as the feed advances the surface under the wheel.

class active_grinding_topo : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_grinding_topo)

  public:

    active_grinding_topo(Coordsys* cs, w_Coordsys* wcs,
                         grinding_topo_params const& params,
                         QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void tick();           // advance the feed-scroll phase, redraw
    void togglePause();    // SPACE
    void resetAnimation(); // R
    void cyclePatch();     // C - toggle the WCD (square) / WMD (wide strip) window
    void viewChanged();

  private:

    // Build the Fig.-1 frame tree from the (real, SI) Tao parameters and read the two
    // relative surface speeds off the GA twist field: v_w = |point_velocity| of the wheel
    // rim, v_s = |point_velocity| of a wafer point at radius r. Sets m_lam_m, m_lam_c,
    // the mark inclination m_beta, then rebuilds the height-field image.
    void build_system();
    void rebuild_image(); // recompute the topography QImage for the current patch + phase

    // the patch window [um] for the current preset, and its placement on the canvas [au]:
    // [m_x0,m_x1] x [m_y0,m_y1], to the right of the text legend, with a colour legend to
    // its right. Set in build_system / cyclePatch.
    double patch_w_um() const;
    double patch_h_um() const;
    void set_layout(); // fill m_x0..m_y1 for the current patch preset

    // jet-style colormap: t in [0,1] -> blue .. cyan .. green .. yellow .. red
    static QRgb jet(double t);

    QPointF auPt(double ax, double ay) const; // abstract-unit point -> widget coords

    // vertical colour legend mapping the jet ramp to surface height [um], drawn in the
    // abstract-unit box [cbx0,cbx1] x [cby0,cby1]
    void drawColorbar(QPainter* qp, double cbx0, double cbx1, double cby0,
                      double cby1) const;

    Coordsys* cs;
    w_Coordsys* wcs;
    grinding_topo_params m_params;
    kinematic_system3dp m_sys;

    double m_lam_m{0.0}; // mark-direction wavelength [um] = v_w / f_b
    double m_lam_c{0.0}; // circumferential wavelength [um] = v_s / f_b
    double m_v_w{0.0};   // wheel rim speed [um/s]
    double m_v_s{0.0};   // wafer surface speed at r [um/s]
    double m_beta{0.0};  // grinding-mark inclination [rad] (GA-computed)
    double m_phase{0.0}; // feed-scroll phase [rad]

    // patch placement on the canvas [abstract units], to the right of the text legend
    double m_x0{0.0}, m_x1{0.0}, m_y0{0.0}, m_y1{0.0};

    bool m_paused{false};
    gt_patch m_patch; // current preset (WCD square / WMD strip)
    QImage m_image;   // cached topography heatmap
    QTimer* m_timer;

    static constexpr int UPDATE_MS = 33;
};
