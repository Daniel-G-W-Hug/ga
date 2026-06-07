// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_double_pendulum.hpp"
#include "active_common.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <QPolygonF>
#include <QTimer>
#include <algorithm>
#include <cmath>
#include <string>

using namespace hd::ga;
using namespace hd::ga::pga;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_double_pendulum::active_double_pendulum(Coordsys* cs, w_Coordsys* wcs,
                                               double_pendulum_params const& params,
                                               QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    build_system();

    connect(wcs, &w_Coordsys::viewResized, this, &active_double_pendulum::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_double_pendulum::tick);
    m_timer->start(UPDATE_MS);
}

// (Re)build the dynamic system at t=0: world -> plate 1 (revolute) -> plate 2 (revolute).
// Geometry: plate 1 hinged to the world at its top-right corner Q1 = (w1/2, h1/2), placed
// at m_params.pivot; plate 2 hinged to plate 1's opposite (bottom-left) corner via its
// own top-right corner Q2 = (w2/2, h2/2).
void active_double_pendulum::build_system()
{
    auto const& p = m_params;
    m_ds = dynamic_system2dp{};
    m_ds.add_frame(static_frame2dp("W")); // inertial root

    vec2dp const Q1{p.w1 / 2.0, p.h1 / 2.0, 1.0};
    vec2dp const cm1_rest{p.pivot.x - Q1.x, p.pivot.y - Q1.y, 1.0}; // rel. to world
    m_ds.add_revolute_body(static_frame2dp("B1", cm1_rest, 0.0),
                           make_plate_body(p.m1, p.w1, p.h1), Q1, p.phi1_0, p.omega1_0);

    vec2dp const Q2{p.w2 / 2.0, p.h2 / 2.0, 1.0};
    // plate 2 cm relative to plate 1: plate 1's BL corner minus plate 2's TR corner
    vec2dp const cm2_rest{-(p.w1 / 2.0 + p.w2 / 2.0), -(p.h1 / 2.0 + p.h2 / 2.0), 1.0};
    m_ds.add_revolute_body(static_frame2dp("B2", cm2_rest, 0.0),
                           make_plate_body(p.m2, p.w2, p.h2), Q2, p.phi2_0, p.omega2_0,
                           m_ds.index_of("B1"));

    m_path.clear();
    m_icr2_path.clear();
    m_t = 0.0;
    m_E0 = m_ds.total_energy();
    m_ds.sync_accelerations(); // make accel queries reflect the forces (for the overlay)
    m_path.push_back(
        move2dp(vec2dp{-p.w2 / 2.0, -p.h2 / 2.0, 1.0},
                m_ds.get_pos_trafo("B2", "W"))); // tip = far corner of plate 2
    recordICR();
}

// Append plate 2's instantaneous centre of rotation to the trace, or a NaN gap marker
// when it is at (or near) infinity -- i.e. when plate 2 is momentarily not rotating, or
// when the ICR has run far from the system. ICR = (V.x/V.z, V.y/V.z) of the world twist
// V.
void active_double_pendulum::recordICR()
{
    vec2dp const gap{std::nan(""), std::nan(""), 0.0};
    twist2dp const V = m_ds.twist_world("B2");
    if (std::abs(V.z) > ICR_MIN_OMEGA) {
        vec2dp const icr{V.x / V.z, V.y / V.z, 1.0};
        bool const near = std::hypot(icr.x - m_params.pivot.x, icr.y - m_params.pivot.y) <
                          ICR_MAX_REACH;
        m_icr2_path.push_back(near ? icr : gap);
    }
    else {
        m_icr2_path.push_back(gap);
    }
    if (m_icr2_path.size() > MAX_PATH) m_icr2_path.pop_front();
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_double_pendulum::tick()
{
    if (m_paused) return;

    int const n = std::max(1, m_params.substeps);
    for (int i = 0; i < n; ++i)
        m_ds.step(m_params.dt);
    m_t += n * m_params.dt;
    m_ds.sync_accelerations(); // refresh accelerations for the force overlay

    // record a trace only while its layer is enabled, so toggling it on starts it "from
    // now on" (the toggle slot clears it) rather than revealing the whole past history.
    if (m_show_trace) {
        m_path.push_back(move2dp(vec2dp{-m_params.w2 / 2.0, -m_params.h2 / 2.0, 1.0},
                                 m_ds.get_pos_trafo("B2", "W")));
        if (m_path.size() > MAX_PATH) m_path.pop_front();
    }
    if (m_show_icr) recordICR();

    prepareGeometryChange();
    update();
}

void active_double_pendulum::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_double_pendulum::toggleTrace()
{
    m_show_trace = !m_show_trace;
    if (m_show_trace) m_path.clear(); // start the trace "from now on"
    update();
}

void active_double_pendulum::toggleICR()
{
    m_show_icr = !m_show_icr;
    if (m_show_icr) m_icr2_path.clear(); // start the ICR trace "from now on"
    update();
}

void active_double_pendulum::toggleForces()
{
    m_show_forces = !m_show_forces;
    update();
}

void active_double_pendulum::resetAnimation()
{
    build_system();
    prepareGeometryChange();
    update();
}

void active_double_pendulum::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Coordinate helper
// ---------------------------------------------------------------------------

QPointF active_double_pendulum::toScreen(vec2dp const& pt) const
{
    return QPointF(cs->x.au_to_w(pt.x), cs->y.au_to_w(pt.y));
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

void active_double_pendulum::drawPlate(QPainter* qp, std::string const& name, double w,
                                       double h, QColor const& fill, QColor const& edge,
                                       char const* label)
{
    auto const M = m_ds.get_pos_trafo(name, "W");
    double const hw = w / 2.0, hh = h / 2.0;
    vec2dp const corners[4] = {
        {-hw, -hh, 1.0}, {hw, -hh, 1.0}, {hw, hh, 1.0}, {-hw, hh, 1.0}};

    QPolygonF poly;
    for (auto const& c : corners)
        poly << toScreen(move2dp(c, M));
    qp->setPen(QPen(edge, 2.0));
    qp->setBrush(QBrush(fill));
    qp->drawPolygon(poly);

    // body-frame axes at the cm (= body origin): e1 red, e2 green
    vec2dp const cm = move2dp(O_2dp, M);
    vec2dp const d1 = move2dp(vec2dp{1.0, 0.0, 0.0}, M);
    vec2dp const d2 = move2dp(vec2dp{0.0, 1.0, 0.0}, M);
    double const len = std::min(w, h) * 0.40;
    QPointF const cm_s = toScreen(cm);
    qp->setBrush(Qt::NoBrush);
    qp->setPen(QPen(QColor(200, 0, 0), 1.6));
    qp->drawPath(
        arrowLine(cm_s, toScreen(vec2dp{cm.x + d1.x * len, cm.y + d1.y * len, 1.0})));
    qp->drawPath(
        arrowHead(cm_s, toScreen(vec2dp{cm.x + d1.x * len, cm.y + d1.y * len, 1.0})));
    qp->setPen(QPen(QColor(0, 140, 0), 1.6));
    qp->drawPath(
        arrowLine(cm_s, toScreen(vec2dp{cm.x + d2.x * len, cm.y + d2.y * len, 1.0})));
    qp->drawPath(
        arrowHead(cm_s, toScreen(vec2dp{cm.x + d2.x * len, cm.y + d2.y * len, 1.0})));

    // cm dot + label
    QFont lbl_font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl_font.setPointSize(9);
    lbl_font.setBold(true);
    qp->setPen(QPen(Qt::black, 2));
    qp->setBrush(QBrush(Qt::black));
    qp->drawEllipse(cm_s, 2.5, 2.5);
    qp->setFont(lbl_font);
    qp->setPen(QPen(Qt::black, 1));
    qp->drawText(cm_s + QPointF(6.0, -6.0), label);
}

// Force overlay (toggled by 'F'): gravity m*g at each cm, and the hinge reaction forces
// reconstructed by a backward Newton-Euler pass (R2 = m2 a2 - m2 g; R1 = m1 a1 - m1 g +
// R2) using the full cm accelerations (m_ds.sync_accelerations() refreshed them after
// step()).
void active_double_pendulum::drawForces(QPainter* qp)
{
    vec2dp const g{0.0, -9.81, 0.0};
    auto const M1 = m_ds.get_pos_trafo("B1", "W");
    auto const M2 = m_ds.get_pos_trafo("B2", "W");
    vec2dp const cm1 = move2dp(O_2dp, M1);
    vec2dp const cm2 = move2dp(O_2dp, M2);
    vec2dp const hinge2 =
        move2dp(vec2dp{-m_params.w1 / 2.0, -m_params.h1 / 2.0, 1.0}, M1);

    // full cm accelerations + backward Newton-Euler hinge reactions
    vec2dp const a1 = m_ds.point_acceleration(cm1, "B1");
    vec2dp const a2 = m_ds.point_acceleration(cm2, "B2");
    vec2dp const fg1{m_params.m1 * g.x, m_params.m1 * g.y, 0.0};
    vec2dp const fg2{m_params.m2 * g.x, m_params.m2 * g.y, 0.0};
    vec2dp const R2{m_params.m2 * a2.x - fg2.x, m_params.m2 * a2.y - fg2.y, 0.0};
    vec2dp const R1{m_params.m1 * a1.x - fg1.x + R2.x, m_params.m1 * a1.y - fg1.y + R2.y,
                    0.0};

    QFont lbl = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl.setPointSize(8);
    lbl.setBold(true);

    auto draw_force = [&](vec2dp const& base, vec2dp const& f, QColor const& col,
                          char const* text) {
        QPointF const b = toScreen(base);
        QPointF const t =
            toScreen(vec2dp{base.x + f.x * FORCE_SCALE, base.y + f.y * FORCE_SCALE, 1.0});
        qp->setBrush(Qt::NoBrush);
        qp->setPen(QPen(col, 2.0));
        qp->drawPath(arrowLine(b, t));
        qp->drawPath(arrowHead(b, t));
        qp->setFont(lbl);
        qp->setPen(QPen(col, 1));
        qp->drawText(t + QPointF(3.0, -3.0), text);
    };

    QColor const col_grav{230, 115, 0};  // orange: gravity at the cm
    QColor const col_react{150, 0, 160}; // purple: hinge reaction force
    draw_force(cm1, fg1, col_grav, "m1 g");
    draw_force(cm2, fg2, col_grav, "m2 g");
    draw_force(m_params.pivot, R1, col_react, "R1");
    draw_force(hinge2, R2, col_react, "R2");
}

// Instantaneous-centre-of-rotation overlay (toggled by 'I'): each plate's ICR is decoded
// from its world twist V as (V.x/V.z, V.y/V.z) -- plate 1's is the fixed pivot, plate 2's
// is a moving point. Drawn as a teal crosshair, with plate 2's ICR leaving a light-grey
// trace (the GA showcase).
void active_double_pendulum::drawICR(QPainter* qp)
{
    // plate 2's ICR trace (light grey, drawn first so the crosshairs sit on top). Skip
    // segments touching a NaN gap or jumping farther than a plate span (the ICR can
    // teleport across the system near the omega -> 0 singularities).
    if (m_icr2_path.size() >= 2) {
        double const maxseg = std::hypot(m_params.w2, m_params.h2) + 1.0;
        qp->setPen(QPen(QColor(150, 150, 150, 170), 1.2, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_icr2_path.size(); ++i) {
            vec2dp const& a = m_icr2_path[i - 1];
            vec2dp const& b = m_icr2_path[i];
            if (std::isnan(a.x) || std::isnan(b.x)) continue;
            if (std::hypot(b.x - a.x, b.y - a.y) > maxseg) continue;
            qp->drawLine(toScreen(a), toScreen(b));
        }
    }

    QFont lbl = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl.setPointSize(8);
    lbl.setBold(true);

    auto draw_icr = [&](char const* frame, char const* text) {
        twist2dp const V = m_ds.twist_world(frame);
        if (std::abs(V.z) < 1.0e-6) return; // ~pure translation -> ICR at infinity
        QPointF const p = toScreen(vec2dp{V.x / V.z, V.y / V.z, 1.0});
        QColor const col{0, 130, 130}; // teal
        qp->setPen(QPen(col, 1.5));
        qp->setBrush(Qt::NoBrush);
        qp->drawEllipse(p, 5.0, 5.0);
        qp->drawLine(p + QPointF(-7.0, 0.0), p + QPointF(7.0, 0.0));
        qp->drawLine(p + QPointF(0.0, -7.0), p + QPointF(0.0, 7.0));
        qp->setFont(lbl);
        qp->setPen(QPen(col, 1));
        qp->drawText(p + QPointF(8.0, -8.0), text);
    };
    draw_icr("B1", "ICR1");
    draw_icr("B2", "ICR2");
}

void active_double_pendulum::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                                   QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    // --- cg / path trace (drawn first, underneath; toggled by 'T') ---
    if (m_show_trace && m_path.size() >= 2) {
        qp->setPen(QPen(QColor(40, 90, 200, 200), 1.4, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_path.size(); ++i)
            qp->drawLine(toScreen(m_path[i - 1]), toScreen(m_path[i]));
    }

    // --- the two plates (with body-frame axes) ---
    drawPlate(qp, "B1", m_params.w1, m_params.h1, QColor(210, 215, 230, 170),
              QColor(70, 90, 150), "B1");
    drawPlate(qp, "B2", m_params.w2, m_params.h2, QColor(230, 215, 210, 170),
              QColor(150, 90, 70), "B2");

    // --- hinges: fixed world pivot (plate 1) + inter-plate joint (plate 1 BL corner) ---
    auto const hinge_dot = [&](vec2dp const& w, QColor const& col) {
        qp->setPen(QPen(col, 2));
        qp->setBrush(QBrush(QColor(255, 255, 255)));
        qp->drawEllipse(toScreen(w), 4.0, 4.0);
        qp->setPen(QPen(col, 2));
        qp->setBrush(QBrush(col));
        qp->drawEllipse(toScreen(w), 1.6, 1.6);
    };
    vec2dp const hinge2_w = move2dp(vec2dp{-m_params.w1 / 2.0, -m_params.h1 / 2.0, 1.0},
                                    m_ds.get_pos_trafo("B1", "W"));
    hinge_dot(m_params.pivot, QColor(20, 20, 20)); // fixed pivot (black)
    hinge_dot(hinge2_w, QColor(150, 90, 70));      // inter-plate hinge

    // --- marked tip of plate 2 ---
    vec2dp const tip_w = move2dp(vec2dp{-m_params.w2 / 2.0, -m_params.h2 / 2.0, 1.0},
                                 m_ds.get_pos_trafo("B2", "W"));
    qp->setPen(QPen(QColor(200, 60, 0), 2));
    qp->setBrush(QBrush(QColor(230, 90, 20)));
    qp->drawEllipse(toScreen(tip_w), 4.0, 4.0);

    // --- overlays (each independently toggleable) ---
    if (m_show_icr) drawICR(qp);       // 'I': instantaneous centre of rotation
    if (m_show_forces) drawForces(qp); // 'F': gravity + hinge reactions

    // --- ground frame W (fixed, solid red/green axes at the world origin) ---
    {
        double const len = 0.8;
        QPointF const o = toScreen(O_2dp);
        qp->setBrush(Qt::NoBrush);
        qp->setPen(QPen(QColor(200, 0, 0), 1.6));
        qp->drawPath(arrowLine(o, toScreen(vec2dp{len, 0.0, 1.0})));
        qp->drawPath(arrowHead(o, toScreen(vec2dp{len, 0.0, 1.0})));
        qp->setPen(QPen(QColor(0, 140, 0), 1.6));
        qp->drawPath(arrowLine(o, toScreen(vec2dp{0.0, len, 1.0})));
        qp->drawPath(arrowHead(o, toScreen(vec2dp{0.0, len, 1.0})));
        QFont lbl_font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        lbl_font.setPointSize(9);
        lbl_font.setBold(true);
        qp->setFont(lbl_font);
        qp->setPen(QPen(Qt::black, 1));
        qp->drawText(o + QPointF(-16.0, 14.0), "W");
    }

    // --- live energy readout (top-left, in pixel coordinates): conservation on screen
    // ---
    {
        double const ke = m_ds.kinetic_energy();
        double const pe = m_ds.potential_energy();
        double const e = ke + pe;
        double const scale = std::max(std::abs(m_E0), 1.0);
        QFont mono = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        mono.setPointSize(12);
        qp->setFont(mono);
        qp->setPen(QPen(QColor(40, 40, 40), 1));
        double const x = cs->x.nmin() + 12.0;
        double const y = cs->y.nmax() + 20.0;
        double const dy = 19.0;
        qp->drawText(QPointF(x, y), QString("t = %1 s").arg(m_t, 0, 'f', 1));
        qp->drawText(QPointF(x, y + dy),
                     QString("KE = %1   PE = %2").arg(ke, 0, 'f', 2).arg(pe, 0, 'f', 2));
        qp->drawText(QPointF(x, y + 2.0 * dy), QString("E = %1   (drift %2)")
                                                   .arg(e, 0, 'f', 3)
                                                   .arg((e - m_E0) / scale, 0, 'e', 1));
    }

    qp->restore();
}

QRectF active_double_pendulum::boundingRect() const
{
    double const reach =
        std::hypot(m_params.w1, m_params.h1) + std::hypot(m_params.w2, m_params.h2) + 1.0;
    double xmin = m_params.pivot.x - reach;
    double xmax = m_params.pivot.x + reach;
    double ymin = m_params.pivot.y - reach;
    double ymax = m_params.pivot.y + reach;

    for (vec2dp const& p : m_path) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }
    xmin = std::min(xmin, 0.0);
    ymin = std::min(ymin, 0.0);

    return QRectF(QPointF(cs->x.au_to_w(xmin), cs->y.au_to_w(ymax)),
                  QPointF(cs->x.au_to_w(xmax), cs->y.au_to_w(ymin)));
}

QPainterPath active_double_pendulum::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
