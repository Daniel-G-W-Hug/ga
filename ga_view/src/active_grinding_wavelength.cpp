// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_grinding_wavelength.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <algorithm>
#include <cmath>

using namespace hd::ga;
using namespace hd::ga::pga;

namespace {
// curve colours (up to four families)
QColor const kCurveCol[4] = {QColor(200, 40, 40), QColor(30, 150, 40),
                             QColor(40, 90, 200), QColor(220, 140, 0)};
} // namespace

active_grinding_wavelength::active_grinding_wavelength(Coordsys* cs, w_Coordsys* wcs,
                                                       grinding_wave_params const& params,
                                                       QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}, m_panel{params.panel}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);
    build_curves();
    connect(wcs, &w_Coordsys::viewResized, this,
            &active_grinding_wavelength::viewChanged);
}

// ---------------------------------------------------------------------------
// Surface speeds from the GA twist field (same Fig.-1 tree as the topo scene)
// ---------------------------------------------------------------------------

active_grinding_wavelength::speeds
active_grinding_wavelength::surface_speeds(double ns_rpm, double nw_rpm,
                                           double r_um) const
{
    double const ns = ns_rpm * pi / 30.0, nw = nw_rpm * pi / 30.0;
    double const Rw = m_params.Rw_um, off = Rw / std::sqrt(2.0);
    double const l3 = 0.3 * Rw, x_a = l3;

    kinematic_system3dp s;
    s.add_frame(static_frame3dp("chuck_ctr_stat"));
    s.add_frame(static_frame3dp("chuck_ctr_rot"),
                kin_state3dp{.omega = vec3dp{ns, 0.0, 0.0, 0.0}},
                s.index_of("chuck_ctr_stat"));
    s.add_frame(static_frame3dp("wafer_top_avg_rot", vec3dp{0.0, 0.0, 0.0, 1.0}),
                kin_state3dp{}, s.index_of("chuck_ctr_rot"));
    s.add_frame(static_frame3dp("spindle_cm_stat", vec3dp{x_a, -off, off, 1.0},
                                vec3dp{0.0, -pi / 2.0, 0.0, 0.0}),
                kin_state3dp{}, s.index_of("chuck_ctr_stat"));
    s.add_frame(static_frame3dp("tool_top_avg_rot", vec3dp{0.0, 0.0, l3, 1.0}),
                kin_state3dp{.omega = vec3dp{0.0, 0.0, -nw, 0.0}},
                s.index_of("spindle_cm_stat"));
    s.add_frame(static_frame3dp("tool_surface_avg_at_R", vec3dp{Rw, 0.0, 0.0, 1.0}),
                kin_state3dp{}, s.index_of("tool_top_avg_rot"));

    vec3dp const g = unitize(
        move3dp(O_3dp, s.get_pos_trafo("tool_surface_avg_at_R", "chuck_ctr_stat")));
    double const vw = to_val(bulk_nrm(s.point_velocity(g, "tool_surface_avg_at_R")));
    double const vs = to_val(
        bulk_nrm(s.point_velocity(vec3dp{0.0, r_um, 0.0, 1.0}, "wafer_top_avg_rot")));
    return {vw, vs};
}

// lambda [um] for the current panel: x is the swept variable, cval the curve-family
// value.
//
// Fixed vibration frequency for the fixed-f_b panels (b) and (d). Tao's Fig.-12 caption
// states f_b = 4000 Hz here, but that is INCONSISTENT with the plotted curves: at 4000 Hz
// lambda = v/f_b gives lambda_m up to ~11000 um (n_w = 2800) and lambda_c up to ~314 um,
// whereas the published Fig. 12(b)/(d) top out near ~7200 um and ~210 um. Those match
// f_b = 6000 Hz (= 11000/1.5, 314/1.5) -- the same f_b Tao uses in Figs. 9/10. So the
// caption's "4000 Hz" is a paper error; we use 6000 Hz to reproduce the figure. (Panels
// (a)/(c) sweep f_b directly and are unaffected.)
static constexpr double FB_FIXED =
    6000.0; // [Hz] (Fig.-12 caption says 4000 -- see above)
double active_grinding_wavelength::wavelength_at(double x, double cval) const
{
    switch (m_panel) {
        case 0: // (a) lambda_c vs f_b, family r [mm]
            return surface_speeds(200.0, 2000.0, cval * 1.0e3).v_s / x;
        case 1: // (b) lambda_c vs wafer speed n_s, family n_w
            return surface_speeds(x, cval, 40.0e3).v_s / FB_FIXED;
        case 2: // (c) lambda_m vs f_b, family r [mm]  (v_w independent of r -> overlap)
            return surface_speeds(200.0, 2000.0, cval * 1.0e3).v_w / x;
        default: // (d) lambda_m vs wafer speed n_s, family n_w  (v_w independent of n_s)
            return surface_speeds(x, cval, 40.0e3).v_w / FB_FIXED;
    }
}

void active_grinding_wavelength::build_curves()
{
    bool const fb_axis = (m_panel == 0 || m_panel == 2);
    m_xmin = fb_axis ? 4000.0 : 100.0;
    m_xmax = fb_axis ? 8000.0 : 300.0;
    bool const r_family = (m_panel == 0 || m_panel == 2);
    m_cvals = r_family ? std::array<double, 4>{30.0, 40.0, 50.0, 60.0}
                       : std::array<double, 4>{1600.0, 2000.0, 2400.0, 2800.0};

    int const N = 81;
    m_ymin = 1e30;
    m_ymax = -1e30;
    for (int c = 0; c < 4; ++c) {
        m_curves[c].clear();
        for (int i = 0; i < N; ++i) {
            double const x = m_xmin + (m_xmax - m_xmin) * double(i) / double(N - 1);
            double const y = wavelength_at(x, m_cvals[c]);
            m_curves[c].push_back(QPointF(x, y));
            m_ymin = std::min(m_ymin, y);
            m_ymax = std::max(m_ymax, y);
        }
    }
    double const pad = 0.08 * (m_ymax - m_ymin);
    m_ymin -= pad;
    m_ymax += pad;
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_grinding_wavelength::cyclePanel()
{
    m_panel = (m_panel + 1) % 4;
    build_curves();
    prepareGeometryChange();
    update();
}

void active_grinding_wavelength::resetAnimation()
{
    m_panel = m_params.panel;
    build_curves();
    prepareGeometryChange();
    update();
}

void active_grinding_wavelength::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Painting
// ---------------------------------------------------------------------------

QPointF active_grinding_wavelength::dataPt(double x, double y) const
{
    double const ax = BX0 + (x - m_xmin) / (m_xmax - m_xmin) * (BX1 - BX0);
    double const ay = BY0 + (y - m_ymin) / (m_ymax - m_ymin) * (BY1 - BY0);
    return QPointF(cs->x.au_to_w(ax), cs->y.au_to_w(ay));
}

void active_grinding_wavelength::paint(QPainter* qp,
                                       QStyleOptionGraphicsItem const* option,
                                       QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    bool const fb_axis = (m_panel == 0 || m_panel == 2);
    bool const lam_c = (m_panel == 0 || m_panel == 1);
    char const* const titles[4] = {
        "(a) lambda_c vs f_b   (family: radius r)",
        "(b) lambda_c vs n_s, f_b=6000   (n_w curves COINCIDE)",
        "(c) lambda_m vs f_b   (r curves COINCIDE)",
        "(d) lambda_m vs n_s, f_b=6000   (FLAT: indep. of n_s)"};

    // plot frame
    QPointF const tl(cs->x.au_to_w(BX0), cs->y.au_to_w(BY1));
    QPointF const br(cs->x.au_to_w(BX1), cs->y.au_to_w(BY0));
    QRectF const box(tl, br);
    qp->setPen(QPen(QColor(40, 40, 40), 1.4));
    qp->setBrush(QBrush(QColor(252, 252, 252)));
    qp->drawRect(box);

    // gridlines + tick labels
    QFont tk = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    tk.setPointSize(9);
    qp->setFont(tk);
    for (int k = 0; k <= 4; ++k) {
        double const f = double(k) / 4.0;
        double const xv = m_xmin + f * (m_xmax - m_xmin);
        double const yv = m_ymin + f * (m_ymax - m_ymin);
        // vertical gridline at xv
        qp->setPen(QPen(QColor(220, 220, 220), 1.0, Qt::DashLine));
        qp->drawLine(dataPt(xv, m_ymin), dataPt(xv, m_ymax));
        qp->drawLine(dataPt(m_xmin, yv), dataPt(m_xmax, yv));
        qp->setPen(QPen(QColor(20, 20, 20), 1));
        qp->drawText(dataPt(xv, m_ymin) + QPointF(-18.0, 16.0),
                     QString("%1").arg(xv, 0, 'f', 0));
        qp->drawText(dataPt(m_xmin, yv) + QPointF(-46.0, 4.0),
                     QString("%1").arg(yv, 0, 'f', 0));
    }

    // curves (uniform width). On the coincident-family panels (b, c) the four curves lie
    // exactly on top of each other -- the title flags this ("... curves COINCIDE") and
    // the legend lists all four families.
    for (int c = 0; c < 4; ++c) {
        QPolygonF poly;
        for (auto const& p : m_curves[c])
            poly << dataPt(p.x(), p.y());
        qp->setPen(QPen(kCurveCol[c], 2.2));
        qp->setBrush(Qt::NoBrush);
        qp->drawPolyline(poly);
    }

    // titles + axis labels
    QFont lab = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lab.setPointSize(10);
    lab.setBold(true);
    qp->setFont(lab);
    qp->setPen(QPen(QColor(20, 20, 20), 1));
    qp->drawText(dataPt(m_xmin, m_ymax) + QPointF(-10.0, -14.0), titles[m_panel]);
    qp->drawText(QPointF(0.5 * (tl.x() + br.x()) - 110.0, br.y() + 34.0),
                 fb_axis ? "vibration frequency f_b [Hz]" : "wafer speed n_s [rpm]");
    qp->save();
    qp->translate(tl.x() - 56.0, 0.5 * (tl.y() + br.y()));
    qp->rotate(-90.0);
    qp->drawText(QPointF(-44.0, 0.0), lam_c ? "lambda_c [um]" : "lambda_m [um]");
    qp->restore();

    // curve legend (top-right inside the box)
    QFont lg = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lg.setPointSize(10);
    qp->setFont(lg);
    bool const r_family = (m_panel == 0 || m_panel == 2);
    double const lx = BX1 - 0.75;
    for (int c = 0; c < 4; ++c) {
        double const ly = BY1 - 0.12 - 0.13 * c;
        QPointF const p = QPointF(cs->x.au_to_w(lx), cs->y.au_to_w(ly));
        qp->setPen(QPen(kCurveCol[c], 3.0));
        qp->drawLine(p + QPointF(0.0, -4.0), p + QPointF(22.0, -4.0));
        qp->setPen(QPen(QColor(20, 20, 20), 1));
        qp->drawText(p + QPointF(28.0, 0.0),
                     r_family ? QString("r = %1 mm").arg(m_cvals[c], 0, 'f', 0)
                              : QString("n_w = %1 rpm").arg(m_cvals[c], 0, 'f', 0));
    }

    // footer: the relation, and that it comes from the GA twist field
    qp->setFont(lab);
    qp->drawText(
        QPointF(cs->x.au_to_w(BX0), cs->y.au_to_w(BY0 - 0.55)),
        "lambda = surface_speed / f_b   (speeds: point_velocity;  C: cycle a-d)");

    qp->restore();
}

QRectF active_grinding_wavelength::boundingRect() const
{
    return QRectF(QPointF(cs->x.au_to_w(BX0 - 0.9), cs->y.au_to_w(BY1 + 0.6)),
                  QPointF(cs->x.au_to_w(BX1 + 0.4), cs->y.au_to_w(BY0 - 0.9)));
}

QPainterPath active_grinding_wavelength::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
