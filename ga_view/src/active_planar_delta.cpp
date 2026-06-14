// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_planar_delta.hpp"
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

active_planar_delta::active_planar_delta(Coordsys* cs, w_Coordsys* wcs,
                                         planar_delta_params const& params,
                                         QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    build_system();

    connect(wcs, &w_Coordsys::viewResized, this, &active_planar_delta::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_planar_delta::tick);
    m_timer->start(UPDATE_MS);
}

// (Re)build the closed loop at the start configuration. Spanning tree (two open chains):
//   world --(revolute SA)--> upper_a --(revolute EA)--> fore_a   [arm A]
//   world --(revolute SB)--> upper_b --(revolute EB)--> fore_b   [arm B]
// closed by pinning the two forearm tips together (one point-coincidence). The dependent
// elbow start angles come from a closed-form circle-circle intersection of the two
// forearm reaches; assemble() then polishes the configuration to g ~ 0.
void active_planar_delta::build_system()
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

    // seed at the t = 0 *driven* config (not the bare nominal): shoulder B carries the
    // phase lead, so shoulderB(0) = shoulder_b0 + amp*sin(phase) != shoulder_b0. Using the
    // driven values here makes the start configuration coincide with the t = 0 point of the
    // Lissajous, so the trace begins continuously (no first-frame jump).
    m_t = 0.0;
    double const phiA0 = shoulderA();
    double const phiB0 = shoulderB();

    // elbow positions reachable from each shoulder at the start shoulder angles
    vec2dp const EA0{p.base_a.x + p.upper * std::cos(phiA0),
                     p.base_a.y + p.upper * std::sin(phiA0), 1.0};
    vec2dp const EB0{p.base_b.x + p.upper * std::cos(phiB0),
                     p.base_b.y + p.upper * std::sin(phiB0), 1.0};
    vec2dp const tip = circ2(EA0, p.fore, EB0, p.fore); // shared effector start position

    // elbow joint angles (forearm direction relative to the upper-arm frame)
    double const thetaA0 = std::atan2(tip.y - EA0.y, tip.x - EA0.x) - phiA0;
    double const thetaB0 = std::atan2(tip.y - EB0.y, tip.x - EB0.x) - phiB0;

    // arm A: shoulder (revolute about base_a) then elbow (revolute about its upper tip)
    m_cl.add_revolute_body(static_frame2dp("SA", p.base_a, 0.0),
                           make_plate_body(1.0, 1.0, 1.0), vec2dp{0.0, 0.0, 1.0},
                           /*phi0*/ phiA0, /*omega0*/ 0.0, m_cl.index_of("W"));
    m_cl.add_revolute_body(static_frame2dp("EA", vec2dp{p.upper, 0.0, 1.0}, 0.0),
                           make_plate_body(1.0, 1.0, 1.0), vec2dp{0.0, 0.0, 1.0},
                           /*phi0*/ thetaA0, /*omega0*/ 0.0, m_cl.index_of("SA"));
    // arm B: shoulder (revolute about base_b) then elbow
    m_cl.add_revolute_body(static_frame2dp("SB", p.base_b, 0.0),
                           make_plate_body(1.0, 1.0, 1.0), vec2dp{0.0, 0.0, 1.0},
                           /*phi0*/ phiB0, /*omega0*/ 0.0, m_cl.index_of("W"));
    m_cl.add_revolute_body(static_frame2dp("EB", vec2dp{p.upper, 0.0, 1.0}, 0.0),
                           make_plate_body(1.0, 1.0, 1.0), vec2dp{0.0, 0.0, 1.0},
                           /*phi0*/ thetaB0, /*omega0*/ 0.0, m_cl.index_of("SB"));

    m_SA = m_cl.index_of("SA");
    m_EA = m_cl.index_of("EA");
    m_SB = m_cl.index_of("SB");
    m_EB = m_cl.index_of("EB");

    // close the loop: forearm tip of arm A == forearm tip of arm B (shared effector)
    m_cl.add_loop_constraint(loop_constraint2dp{m_EA, vec2dp{p.fore, 0.0, 1.0}, m_EB,
                                                vec2dp{p.fore, 0.0, 1.0},
                                                constraint2dp::coincidence});

    m_cl.assemble(/*driven*/ {m_SA, m_SB}); // snap to a consistent closed configuration

    m_path.clear();
    m_path.push_back(effector());
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_planar_delta::tick()
{
    if (m_paused) return;

    int const n = std::max(1, m_params.substeps);
    m_t += m_params.dt * n;

    m_cl.set_joint(m_SA, shoulderA());      // advance the two driven shoulders
    m_cl.set_joint(m_SB, shoulderB());      //
    m_cl.assemble(/*driven*/ {m_SA, m_SB}); // re-solve elbows (loop stays closed)

    if (m_show_trace) {
        m_path.push_back(effector());
        if (m_path.size() > MAX_PATH) m_path.pop_front();
    }

    prepareGeometryChange();
    update();
}

void active_planar_delta::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_planar_delta::toggleTrace()
{
    m_show_trace = !m_show_trace;
    if (m_show_trace) m_path.clear(); // start the trace "from now on"
    update();
}

void active_planar_delta::resetAnimation()
{
    build_system();
    prepareGeometryChange();
    update();
}

void active_planar_delta::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

double active_planar_delta::shoulderA() const
{
    return m_params.shoulder_a0 +
           m_params.amp * std::sin(2.0 * M_PI * m_t / m_params.period);
}

double active_planar_delta::shoulderB() const
{
    return m_params.shoulder_b0 +
           m_params.amp * std::sin(2.0 * M_PI * m_t / m_params.period + m_params.phase);
}

QPointF active_planar_delta::toScreen(vec2dp const& pt) const
{
    return QPointF(cs->x.au_to_w(pt.x), cs->y.au_to_w(pt.y));
}

vec2dp active_planar_delta::worldPoint(std::string const& name, vec2dp const& p)
{
    return move2dp(p, m_cl.system().get_pos_trafo(name, "W"));
}

vec2dp active_planar_delta::pivotSA() const { return m_params.base_a; }
vec2dp active_planar_delta::pivotSB() const { return m_params.base_b; }
vec2dp active_planar_delta::elbowA() { return worldPoint("EA", vec2dp{0.0, 0.0, 1.0}); }
vec2dp active_planar_delta::elbowB() { return worldPoint("EB", vec2dp{0.0, 0.0, 1.0}); }
vec2dp active_planar_delta::effector()
{
    return worldPoint("EA", vec2dp{m_params.fore, 0.0, 1.0});
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

void active_planar_delta::drawBar(QPainter* qp, vec2dp const& a, vec2dp const& b,
                                  QColor const& col, double width) const
{
    qp->setPen(QPen(col, width, Qt::SolidLine, Qt::RoundCap));
    qp->drawLine(toScreen(a), toScreen(b));
}

void active_planar_delta::drawPivot(QPainter* qp, vec2dp const& w, QColor const& col,
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

void active_planar_delta::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                                QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    vec2dp const SA = pivotSA();
    vec2dp const SB = pivotSB();
    vec2dp const EA = elbowA();
    vec2dp const EB = elbowB();
    vec2dp const P = effector();

    QColor const col_ground{120, 120, 120};
    QColor const col_arm_a{70, 90, 150}; // blue   : left arm
    QColor const col_arm_b{0, 130, 90};  // green  : right arm
    QColor const col_pin{210, 90, 20};   // orange : the shared end-effector pin

    // --- effector-path trace (drawn first, underneath; toggled by 'T') ---
    if (m_show_trace && m_path.size() >= 2) {
        qp->setPen(QPen(QColor(180, 70, 20, 200), 1.4, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_path.size(); ++i)
            qp->drawLine(toScreen(m_path[i - 1]), toScreen(m_path[i]));
    }

    // --- the fixed base line between the two grounded shoulders ---
    drawBar(qp, SA, SB, col_ground, 3.0);

    // --- the two 2-link arms ---
    drawBar(qp, SA, EA, col_arm_a, 4.0); // arm A upper
    drawBar(qp, EA, P, col_arm_a, 4.0);  // arm A forearm
    drawBar(qp, SB, EB, col_arm_b, 4.0); // arm B upper
    drawBar(qp, EB, P, col_arm_b, 4.0);  // arm B forearm

    // --- pivots and joints ---
    drawPivot(qp, SA, QColor(20, 20, 20), /*grounded*/ true);
    drawPivot(qp, SB, QColor(20, 20, 20), /*grounded*/ true);
    drawPivot(qp, EA, col_arm_a, /*grounded*/ false);
    drawPivot(qp, EB, col_arm_b, /*grounded*/ false);
    drawPivot(qp, P, col_pin, /*grounded*/ false); // the shared effector pin

    // --- the end-effector ---
    qp->setPen(QPen(QColor(180, 70, 20), 2));
    qp->setBrush(QBrush(QColor(230, 110, 40)));
    qp->drawEllipse(toScreen(P), 3.5, 3.5);

    // --- frame labels ---
    QFont lbl = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl.setPointSize(9);
    lbl.setBold(true);
    qp->setFont(lbl);
    qp->setPen(QPen(Qt::black, 1));
    qp->drawText(toScreen(SA) + QPointF(-26.0, 16.0), "SA");
    qp->drawText(toScreen(SB) + QPointF(10.0, 16.0), "SB");
    qp->drawText(toScreen(EA) + QPointF(-26.0, -6.0), "EA");
    qp->drawText(toScreen(EB) + QPointF(10.0, -6.0), "EB");
    qp->drawText(toScreen(P) + QPointF(8.0, -8.0), "P");

    // --- live readout (top-left, in pixel coordinates) ---
    {
        QFont mono = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        mono.setPointSize(12);
        qp->setFont(mono);
        qp->setPen(QPen(QColor(40, 40, 40), 1));
        double const x = cs->x.nmin() + 12.0;
        double const y = cs->y.nmax() + 20.0;
        double const dy = 19.0;
        double const a_deg = std::fmod(qRadiansToDegrees(shoulderA()), 360.0);
        double const b_deg = std::fmod(qRadiansToDegrees(shoulderB()), 360.0);
        qp->drawText(QPointF(x, y), QString("t = %1 s").arg(m_t, 0, 'f', 1));
        qp->drawText(QPointF(x, y + dy), QString("shoulders = %1 / %2 deg")
                                             .arg(a_deg, 0, 'f', 0)
                                             .arg(b_deg, 0, 'f', 0));
        qp->drawText(QPointF(x, y + 2.0 * dy),
                     QString("closure ||g|| = %1").arg(m_cl.residual_norm(), 0, 'e', 1));
    }

    qp->restore();
}

QRectF active_planar_delta::boundingRect() const
{
    double const reach = m_params.upper + m_params.fore + 1.0;
    double xmin = std::min(m_params.base_a.x, m_params.base_b.x) - reach;
    double xmax = std::max(m_params.base_a.x, m_params.base_b.x) + reach;
    double ymin = std::min(m_params.base_a.y, m_params.base_b.y) - reach;
    double ymax = std::max(m_params.base_a.y, m_params.base_b.y) + reach;

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

QPainterPath active_planar_delta::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
