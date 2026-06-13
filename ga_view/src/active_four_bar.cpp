// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_four_bar.hpp"
#include "active_common.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <QtMath>
#include <algorithm>
#include <cmath>
#include <string>

using namespace hd::ga;
using namespace hd::ga::pga;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_four_bar::active_four_bar(Coordsys* cs, w_Coordsys* wcs,
                                 four_bar_params const& params, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    build_system();

    connect(wcs, &w_Coordsys::viewResized, this, &active_four_bar::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_four_bar::tick);
    m_timer->start(UPDATE_MS);
}

// (Re)build the closed loop at the crank start angle. Spanning tree:
//   world --(revolute O2)--> crank --(revolute A)--> coupler   [branch 1]
//   world --(revolute O4)--> rocker                            [branch 2]
// closed by pinning the coupler tip (coupler length along the coupler frame) to the
// rocker tip (rocker length along the rocker frame). The dependent coupler / rocker start
// angles come from a closed-form circle-circle intersection of the two reachable tips;
// assemble() then polishes the configuration to g ~ 0.
void active_four_bar::build_system()
{
    auto const& p = m_params;
    m_cl = closed_loop_system2dp{};
    m_cl.add_frame(static_frame2dp("W")); // inertial root / ground

    // upper intersection of circle(centre c0, radius r0) and circle(c1, r1)
    auto circ2 = [](vec2dp const& c0, double r0, vec2dp const& c1, double r1) -> vec2dp {
        double const dx = c1.x - c0.x, dy = c1.y - c0.y;
        double const d = std::hypot(dx, dy);
        double const a = (r0 * r0 - r1 * r1 + d * d) / (2.0 * d);
        double const h2 = r0 * r0 - a * a;
        double const h = h2 > 0.0 ? std::sqrt(h2) : 0.0;
        double const px = c0.x + a * dx / d, py = c0.y + a * dy / d;
        return vec2dp{px - h * dy / d, py + h * dx / d, 1.0}; // perp = (-dy, dx)/d
    };

    vec2dp const O2{p.base.x, p.base.y, 1.0};
    vec2dp const O4{p.base.x + p.ground, p.base.y, 1.0};
    vec2dp const A0{p.base.x + p.crank, p.base.y, 1.0}; // crank tip at crank angle 0
    vec2dp const B0 = circ2(A0, p.coupler, O4, p.rocker);
    double const theta3_0 =
        std::atan2(B0.y - A0.y, B0.x - A0.x); // coupler vs crank (==0)
    double const theta4_0 = std::atan2(B0.y - O4.y, B0.x - O4.x); // rocker

    // branch 1: crank (revolute about O2) then coupler (revolute about A)
    m_cl.add_revolute_body(static_frame2dp("CR", O2, 0.0), make_plate_body(1.0, 1.0, 1.0),
                           vec2dp{0.0, 0.0, 1.0}, /*phi0*/ 0.0, /*omega0*/ 0.0,
                           m_cl.index_of("W"));
    m_cl.add_revolute_body(static_frame2dp("CO", vec2dp{p.crank, 0.0, 1.0}, 0.0),
                           make_plate_body(1.0, 1.0, 1.0), vec2dp{0.0, 0.0, 1.0},
                           /*phi0*/ theta3_0, /*omega0*/ 0.0, m_cl.index_of("CR"));
    // branch 2: rocker (revolute about O4)
    m_cl.add_revolute_body(static_frame2dp("RO", O4, 0.0), make_plate_body(1.0, 1.0, 1.0),
                           vec2dp{0.0, 0.0, 1.0}, /*phi0*/ theta4_0, /*omega0*/ 0.0,
                           m_cl.index_of("W"));

    m_CR = m_cl.index_of("CR");
    m_CO = m_cl.index_of("CO");
    m_RO = m_cl.index_of("RO");

    // close the loop: coupler tip == rocker tip
    m_cl.add_loop_constraint(loop_constraint2dp{m_CO, vec2dp{p.coupler, 0.0, 1.0}, m_RO,
                                                vec2dp{p.rocker, 0.0, 1.0},
                                                constraint2dp::coincidence});

    m_crank0 = 0.0;
    m_crank = m_crank0;
    m_t = 0.0;
    m_cl.assemble(/*driven*/ {m_CR}); // snap to a consistent closed configuration

    m_path.clear();
    m_path.push_back(couplerPoint());
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_four_bar::tick()
{
    if (m_paused) return;

    int const n = std::max(1, m_params.substeps);
    m_crank += m_params.crank_rate * m_params.dt * n;
    m_t += m_params.dt * n;

    m_cl.set_joint(m_CR, m_crank);    // advance the driven crank
    m_cl.assemble(/*driven*/ {m_CR}); // re-solve coupler + rocker (loop stays closed)

    if (m_show_trace) {
        m_path.push_back(couplerPoint());
        if (m_path.size() > MAX_PATH) m_path.pop_front();
    }

    prepareGeometryChange();
    update();
}

void active_four_bar::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_four_bar::toggleTrace()
{
    m_show_trace = !m_show_trace;
    if (m_show_trace) m_path.clear(); // start the trace "from now on"
    update();
}

void active_four_bar::resetAnimation()
{
    build_system();
    prepareGeometryChange();
    update();
}

void active_four_bar::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

QPointF active_four_bar::toScreen(vec2dp const& pt) const
{
    return QPointF(cs->x.au_to_w(pt.x), cs->y.au_to_w(pt.y));
}

vec2dp active_four_bar::worldPoint(std::string const& name, vec2dp const& p)
{
    return move2dp(p, m_cl.system().get_pos_trafo(name, "W"));
}

vec2dp active_four_bar::pivotO2() const
{
    return vec2dp{m_params.base.x, m_params.base.y, 1.0};
}
vec2dp active_four_bar::pivotO4() const
{
    return vec2dp{m_params.base.x + m_params.ground, m_params.base.y, 1.0};
}
vec2dp active_four_bar::jointA() { return worldPoint("CO", vec2dp{0.0, 0.0, 1.0}); }
vec2dp active_four_bar::jointB()
{
    return worldPoint("CO", vec2dp{m_params.coupler, 0.0, 1.0});
}
vec2dp active_four_bar::couplerPoint() { return worldPoint("CO", m_params.coupler_pt); }

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

void active_four_bar::drawBar(QPainter* qp, vec2dp const& a, vec2dp const& b,
                              QColor const& col, double width) const
{
    qp->setPen(QPen(col, width, Qt::SolidLine, Qt::RoundCap));
    qp->drawLine(toScreen(a), toScreen(b));
}

void active_four_bar::drawPivot(QPainter* qp, vec2dp const& w, QColor const& col,
                                bool grounded) const
{
    QPointF const s = toScreen(w);
    if (grounded) {
        // small ground hatch triangle under the fixed pivot
        QPolygonF tri;
        tri << s + QPointF(0.0, 0.0) << s + QPointF(-7.0, 12.0) << s + QPointF(7.0, 12.0);
        qp->setPen(QPen(QColor(90, 90, 90), 1.5));
        qp->setBrush(QBrush(QColor(200, 200, 200)));
        qp->drawPolygon(tri);
    }
    qp->setPen(QPen(col, 2.0));
    qp->setBrush(QBrush(Qt::white));
    qp->drawEllipse(s, 4.0, 4.0);
    qp->setBrush(QBrush(col));
    qp->drawEllipse(s, 1.6, 1.6);
}

void active_four_bar::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                            QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    vec2dp const O2 = pivotO2();
    vec2dp const O4 = pivotO4();
    vec2dp const A = jointA();
    vec2dp const B = jointB();
    vec2dp const P = couplerPoint();

    QColor const col_ground{120, 120, 120};
    QColor const col_crank{70, 90, 150};   // blue
    QColor const col_coupler{150, 90, 70}; // brown
    QColor const col_rocker{0, 130, 90};   // teal-green
    QColor const col_pin{210, 90, 20};     // orange: the closure pin

    // --- coupler-curve trace (drawn first, underneath; toggled by 'T') ---
    if (m_show_trace && m_path.size() >= 2) {
        qp->setPen(QPen(QColor(40, 90, 200, 200), 1.4, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_path.size(); ++i)
            qp->drawLine(toScreen(m_path[i - 1]), toScreen(m_path[i]));
    }

    // --- the four bars ---
    drawBar(qp, O2, O4, col_ground, 3.0); // ground (fixed)
    drawBar(qp, O2, A, col_crank, 4.0);   // crank (driver)
    drawBar(qp, A, B, col_coupler, 4.0);  // coupler
    drawBar(qp, O4, B, col_rocker, 4.0);  // rocker

    // --- the rigid coupler triangle A-B-P (dashed), carrying the traced point P ---
    qp->setPen(QPen(QColor(150, 90, 70, 140), 1.2, Qt::DashLine));
    qp->drawLine(toScreen(A), toScreen(P));
    qp->drawLine(toScreen(B), toScreen(P));

    // --- pivots and joints ---
    drawPivot(qp, O2, QColor(20, 20, 20), /*grounded*/ true);
    drawPivot(qp, O4, QColor(20, 20, 20), /*grounded*/ true);
    drawPivot(qp, A, col_crank, /*grounded*/ false);
    drawPivot(qp, B, col_pin, /*grounded*/ false); // the closure pin

    // --- the traced coupler point ---
    qp->setPen(QPen(QColor(40, 90, 200), 2));
    qp->setBrush(QBrush(QColor(60, 110, 220)));
    qp->drawEllipse(toScreen(P), 3.5, 3.5);

    // --- frame labels ---
    QFont lbl = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl.setPointSize(9);
    lbl.setBold(true);
    qp->setFont(lbl);
    qp->setPen(QPen(Qt::black, 1));
    qp->drawText(toScreen(O2) + QPointF(-22.0, 16.0), "O2");
    qp->drawText(toScreen(O4) + QPointF(8.0, 16.0), "O4");
    qp->drawText(toScreen(A) + QPointF(6.0, -6.0), "A");
    qp->drawText(toScreen(B) + QPointF(6.0, -6.0), "B");
    qp->drawText(toScreen(P) + QPointF(6.0, -6.0), "P");

    // --- live readout (top-left, in pixel coordinates) ---
    {
        QFont mono = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        mono.setPointSize(12);
        qp->setFont(mono);
        qp->setPen(QPen(QColor(40, 40, 40), 1));
        double const x = cs->x.nmin() + 12.0;
        double const y = cs->y.nmax() + 20.0;
        double const dy = 19.0;
        double const crank_deg = std::fmod(qRadiansToDegrees(m_crank), 360.0);
        qp->drawText(QPointF(x, y), QString("t = %1 s").arg(m_t, 0, 'f', 1));
        qp->drawText(QPointF(x, y + dy),
                     QString("crank = %1 deg").arg(crank_deg, 0, 'f', 0));
        qp->drawText(QPointF(x, y + 2.0 * dy),
                     QString("closure ||g|| = %1").arg(m_cl.residual_norm(), 0, 'e', 1));
    }

    qp->restore();
}

QRectF active_four_bar::boundingRect() const
{
    double const reach = m_params.crank + m_params.coupler + 1.0;
    double const bx = m_params.base.x, by = m_params.base.y;
    double xmin = bx - reach;
    double xmax = bx + m_params.ground + reach;
    double ymin = by - reach;
    double ymax = by + reach;

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

QPainterPath active_four_bar::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
