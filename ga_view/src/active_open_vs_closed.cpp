// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_open_vs_closed.hpp"
#include "active_common.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

using namespace hd::ga;
using namespace hd::ga::pga;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_open_vs_closed::active_open_vs_closed(Coordsys* cs, w_Coordsys* wcs,
                                             open_vs_closed_params const& params,
                                             QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    build_system();

    connect(wcs, &w_Coordsys::viewResized, this, &active_open_vs_closed::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_open_vs_closed::tick);
    m_timer->start(UPDATE_MS);
}

// Build a 3-link arm (world -> L1 shoulder -> L2 elbow -> L3 wrist) at its nominal pose,
// rooted at `base`. Each joint is a revolute about its own frame origin.
void active_open_vs_closed::build_arm(closed_loop_system2dp& cl, vec2dp const& base)
{
    auto const& p = m_params;
    auto const body =
        make_plate_body(1.0, 1.0, 1.0); // inertia irrelevant (kinematic drive)
    vec2dp const pivot{0.0, 0.0, 1.0};
    vec2dp const along{p.link, 0.0, 1.0}; // next joint, one link out along the local e1

    cl.add_frame(static_frame2dp("W")); // root / ground
    cl.add_revolute_body(static_frame2dp("L1", base, 0.0), body, pivot, p.shoulder0, 0.0,
                         cl.index_of("W"));
    cl.add_revolute_body(static_frame2dp("L2", along, 0.0), body, pivot, p.elbow0, 0.0,
                         cl.index_of("L1"));
    cl.add_revolute_body(static_frame2dp("L3", along, 0.0), body, pivot, p.wrist0, 0.0,
                         cl.index_of("L2"));
}

void active_open_vs_closed::build_system()
{
    // LEFT: open arm -- no loop constraint, so it stays rigid at the nominal elbow/wrist
    m_open = closed_loop_system2dp{};
    build_arm(m_open, m_params.base_open);

    // RIGHT: closed arm -- pin the hand to the world point where it nominally sits, so
    // the elbow/wrist become the dependent joints solved each frame.
    m_closed = closed_loop_system2dp{};
    build_arm(m_closed, m_params.base_closed);
    m_pin = handPt(m_closed); // the hand's nominal world position
    m_closed.add_loop_constraint(
        loop_constraint2dp{m_closed.index_of("L3"), vec2dp{m_params.link, 0.0, 1.0},
                           m_closed.index_of("W"), m_pin, constraint2dp::coincidence});

    m_shoulder = m_closed.index_of("L1");
    m_closed.assemble(/*driven*/ {m_shoulder}); // already consistent (a no-op snap)

    m_t = 0.0;
    m_open_hand_path.clear();
    m_open_hand_path.push_back(handPt(m_open));
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_open_vs_closed::tick()
{
    if (m_paused) return;

    int const n = std::max(1, m_params.substeps);
    m_t += m_params.dt * n;

    // one shared driver: oscillate the shoulder of BOTH arms identically
    double const theta =
        m_params.shoulder0 +
        m_params.amp * std::sin(2.0 * std::numbers::pi * m_t / m_params.period);

    m_open.set_joint(m_shoulder, theta);   // open: elbow/wrist fixed -> rigid swing
    m_closed.set_joint(m_shoulder, theta); // closed: drive the shoulder ...
    m_closed.assemble(
        /*driven*/ {m_shoulder}); // ... and re-solve elbow/wrist (hand pinned)

    if (m_show_trace) {
        m_open_hand_path.push_back(handPt(m_open));
        if (m_open_hand_path.size() > MAX_PATH) m_open_hand_path.pop_front();
    }

    prepareGeometryChange();
    update();
}

void active_open_vs_closed::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_open_vs_closed::toggleTrace()
{
    m_show_trace = !m_show_trace;
    if (m_show_trace) m_open_hand_path.clear();
    update();
}

void active_open_vs_closed::resetAnimation()
{
    build_system();
    prepareGeometryChange();
    update();
}

void active_open_vs_closed::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

QPointF active_open_vs_closed::toScreen(vec2dp const& pt) const
{
    return QPointF(cs->x.au_to_w(pt.x), cs->y.au_to_w(pt.y));
}

vec2dp active_open_vs_closed::worldPt(closed_loop_system2dp& cl, std::string const& name,
                                      vec2dp const& p)
{
    return move2dp(p, cl.system().get_pos_trafo(name, "W"));
}

vec2dp active_open_vs_closed::handPt(closed_loop_system2dp& cl)
{
    return worldPt(cl, "L3", vec2dp{m_params.link, 0.0, 1.0});
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

void active_open_vs_closed::drawBar(QPainter* qp, vec2dp const& a, vec2dp const& b,
                                    QColor const& col, double width) const
{
    qp->setPen(QPen(col, width, Qt::SolidLine, Qt::RoundCap));
    qp->drawLine(toScreen(a), toScreen(b));
}

void active_open_vs_closed::drawPivot(QPainter* qp, vec2dp const& w, QColor const& col,
                                      bool grounded) const
{
    QPointF const s = toScreen(w);
    if (grounded) {
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

void active_open_vs_closed::drawArm(QPainter* qp, closed_loop_system2dp& cl)
{
    vec2dp const sh = worldPt(cl, "L1", vec2dp{0.0, 0.0, 1.0}); // shoulder (base)
    vec2dp const el = worldPt(cl, "L2", vec2dp{0.0, 0.0, 1.0}); // elbow
    vec2dp const wr = worldPt(cl, "L3", vec2dp{0.0, 0.0, 1.0}); // wrist
    vec2dp const hd = handPt(cl);                               // hand (wrist tip)

    drawBar(qp, sh, el, QColor(70, 90, 150), 4.0); // upper arm (blue)
    drawBar(qp, el, wr, QColor(150, 90, 70), 4.0); // forearm (brown)
    drawBar(qp, wr, hd, QColor(0, 130, 90), 4.0);  // hand link (green)

    drawPivot(qp, sh, QColor(20, 20, 20), /*grounded*/ true);
    drawPivot(qp, el, QColor(70, 90, 150), false);
    drawPivot(qp, wr, QColor(150, 90, 70), false);
}

void active_open_vs_closed::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                                  QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    QFont lbl = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl.setPointSize(11);
    lbl.setBold(true);

    // --- vertical divider between the two panels ---
    {
        double const xmid = 0.5 * (m_params.base_open.x + m_params.base_closed.x);
        qp->setPen(QPen(QColor(200, 200, 200), 1.0, Qt::DashLine));
        qp->drawLine(toScreen(vec2dp{xmid, cs->y.min(), 1.0}),
                     toScreen(vec2dp{xmid, cs->y.max(), 1.0}));
    }

    // --- OPEN arm: free-hand trace (under the arm), then the arm ---
    if (m_show_trace && m_open_hand_path.size() >= 2) {
        qp->setPen(QPen(QColor(40, 90, 200, 200), 1.4, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_open_hand_path.size(); ++i)
            qp->drawLine(toScreen(m_open_hand_path[i - 1]),
                         toScreen(m_open_hand_path[i]));
    }
    drawArm(qp, m_open);
    {
        vec2dp const hd = handPt(m_open);
        qp->setPen(QPen(QColor(40, 90, 200), 2));
        qp->setBrush(QBrush(QColor(60, 110, 220)));
        qp->drawEllipse(toScreen(hd), 4.5, 4.5); // free hand (blue, moving)
    }

    // --- CLOSED arm + the fixed pin ---
    drawArm(qp, m_closed);
    {
        // the pin: a fixed target (cross) with the planted hand on top
        QPointF const pin = toScreen(m_pin);
        qp->setPen(QPen(QColor(210, 90, 20), 2.0));
        qp->setBrush(Qt::NoBrush);
        qp->drawEllipse(pin, 8.0, 8.0);
        qp->drawLine(pin + QPointF(-11.0, 0.0), pin + QPointF(11.0, 0.0));
        qp->drawLine(pin + QPointF(0.0, -11.0), pin + QPointF(0.0, 11.0));
        qp->setBrush(QBrush(QColor(230, 110, 30)));
        qp->drawEllipse(pin, 4.5, 4.5); // pinned hand (orange, fixed)
    }

    // --- panel labels (world coords, above each shoulder) ---
    qp->setFont(lbl);
    qp->setPen(QPen(QColor(40, 40, 40), 1));
    qp->drawText(toScreen(vec2dp{m_params.base_open.x - 1.0, 1.9, 1.0}),
                 "OPEN CHAIN  (free hand)");
    qp->drawText(toScreen(vec2dp{m_params.base_closed.x - 1.0, 1.9, 1.0}),
                 "CLOSED LOOP  (hand pinned)");

    // --- live readout (top-left, in pixel coordinates) ---
    {
        QFont mono = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        mono.setPointSize(12);
        qp->setFont(mono);
        qp->setPen(QPen(QColor(40, 40, 40), 1));
        double const x = cs->x.nmin() + 12.0;
        double const y = cs->y.nmax() + 20.0;
        double const dy = 19.0;
        qp->drawText(QPointF(x, y), QString("t = %1 s").arg(m_t, 0, 'f', 1));
        qp->drawText(QPointF(x, y + dy), "one shared shoulder drive");
        qp->drawText(
            QPointF(x, y + 2.0 * dy),
            QString("closed-loop ||g|| = %1").arg(m_closed.residual_norm(), 0, 'e', 1));
    }

    qp->restore();
}

QRectF active_open_vs_closed::boundingRect() const
{
    double const reach = 3.0 * m_params.link + 1.0;
    double const xl = m_params.base_open.x, xr = m_params.base_closed.x;
    double const yb = std::min(m_params.base_open.y, m_params.base_closed.y);
    double xmin = std::min(xl, xr) - reach;
    double xmax = std::max(xl, xr) + reach;
    double ymin = yb - reach;
    double ymax = yb + reach;

    for (vec2dp const& p : m_open_hand_path) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }

    return QRectF(QPointF(cs->x.au_to_w(xmin), cs->y.au_to_w(ymax)),
                  QPointF(cs->x.au_to_w(xmax), cs->y.au_to_w(ymin)));
}

QPainterPath active_open_vs_closed::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
