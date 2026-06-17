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
#include <QWidget>

#include <array>
#include <string>
#include <vector>

// active_grinding_wavelength visualizes Tao 2022, Fig. 12: how the grinding-waviness
// wavelengths depend on the process parameters. Every point is lambda = (relative surface
// speed)/f_b, with the wheel speed v_w and wafer speed v_s read from
// kinematic_system3dp::point_velocity of the Fig.-1 frames -- so the four panels are line
// plots built entirely from the GA twist field. They reproduce the paper's trends AND
// make the GA-derived independence relations visible:
//   (a) lambda_c vs f_b, family in radius r            -> ~1/f_b, spaced by r
//   (b) lambda_c vs wafer speed, family in n_w         -> linear; the n_w curves OVERLAP
//   (c) lambda_m vs f_b, family in r                   -> ~1/f_b; the r curves OVERLAP
//   (d) lambda_m vs wafer speed, family in n_w         -> FLAT lines, spaced by n_w
// C cycles the four panels.

class active_grinding_wavelength : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    Q_DISABLE_COPY_MOVE(active_grinding_wavelength)

  public:

    active_grinding_wavelength(Coordsys* cs, w_Coordsys* wcs,
                               grinding_wave_params const& params,
                               QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void cyclePanel();     // C - cycle Fig.-12 panels (a)->(b)->(c)->(d)
    void resetAnimation(); // R - back to panel (a)
    void viewChanged();

  private:

    struct speeds {
        double v_w; // wheel rim speed n_w*R_w [um/s]
        double v_s; // wafer speed n_s*r       [um/s]
    };

    // surface speeds from the GA twist field for the given (n_s, n_w, r)
    speeds surface_speeds(double ns_rpm, double nw_rpm, double r_um) const;

    // wavelength [um] at swept value x with curve-family value cval, for the current
    // panel
    double wavelength_at(double x, double cval) const;

    void
    build_curves(); // (re)compute the four curves + axis ranges for the current panel

    QPointF dataPt(double x, double y) const; // data coords -> widget coords

    Coordsys* cs;
    w_Coordsys* wcs;
    grinding_wave_params m_params;

    int m_panel{0}; // 0=a, 1=b, 2=c, 3=d

    // computed curves (data coords) + axis ranges for the current panel
    std::array<std::vector<QPointF>, 4> m_curves;
    std::array<double, 4> m_cvals{}; // family values (r [mm] or n_w [rpm])
    double m_xmin{0.0}, m_xmax{1.0}, m_ymin{0.0}, m_ymax{1.0};

    // plot box on the canvas [abstract units], right of the text legend
    static constexpr double BX0 = 0.25, BX1 = 2.7, BY0 = -1.05, BY1 = 1.55;
};
