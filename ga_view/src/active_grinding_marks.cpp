// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_grinding_marks.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <QTimer>
#include <cmath>

using namespace hd::ga;
using namespace hd::ga::pga;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_grinding_marks::active_grinding_marks(Coordsys* cs, w_Coordsys* wcs,
                                             grinding_marks_params const& params,
                                             QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    build_system();

    connect(wcs, &w_Coordsys::viewResized, this, &active_grinding_marks::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_grinding_marks::tick);
    m_timer->start(UPDATE_MS);
}

// ---------------------------------------------------------------------------
// System (re)build and time advance
// ---------------------------------------------------------------------------

// The six-frame Tao Fig.-1 tree, normalized to R = wheel radius = wafer radius = 1. The
// axial quantities (l3, the spindle's axial offset x_a) lie along e1 and so drop out of
// the e423 projection -- they are kept only so the tree matches the geometry test.
void active_grinding_marks::build_system()
{
    double const R = m_params.R;
    double const off = R / std::sqrt(2.0); // |spindle offset|_yz == R (rim hits centre)
    double const l3 = 0.3 * R; // tool surface vs. spindle CM (along e1, hidden)
    double const x_a = l3;     // axial infeed (along e1, hidden)

    m_sys = kinematic_system3dp{};
    m_sys.add_frame(static_frame3dp("chuck_ctr_stat")); // inertial root

    // chain 1 (wafer): chuck spin about e1, then the wafer top surface (at the chuck
    // axis)
    m_sys.add_frame(static_frame3dp("chuck_ctr_rot"),
                    kin_state3dp{.omega = vec3dp{m_params.ns, 0.0, 0.0, 0.0}},
                    m_sys.index_of("chuck_ctr_stat"));
    m_sys.add_frame(static_frame3dp("wafer_top_avg_rot", vec3dp{0.0, 0.0, 0.0, 1.0}),
                    kin_state3dp{}, m_sys.index_of("chuck_ctr_rot"));

    // chain 2 (spindle/wheel): CM reoriented -90 deg about e2, tool surface at +l3 along
    // the spindle z spinning at n_w, rim grain at +R along the local x
    m_sys.add_frame(static_frame3dp("spindle_cm_stat", vec3dp{x_a, -off, off, 1.0},
                                    vec3dp{0.0, -pi / 2.0, 0.0, 0.0}),
                    kin_state3dp{}, m_sys.index_of("chuck_ctr_stat"));
    m_sys.add_frame(static_frame3dp("tool_top_avg_rot", vec3dp{0.0, 0.0, l3, 1.0}),
                    kin_state3dp{.omega = vec3dp{0.0, 0.0, m_params.nw, 0.0}},
                    m_sys.index_of("spindle_cm_stat"));
    m_sys.add_frame(static_frame3dp("tool_surface_avg_at_R", vec3dp{R, 0.0, 0.0, 1.0}),
                    kin_state3dp{}, m_sys.index_of("tool_top_avg_rot"));

    m_path.clear();
    m_t = 0.0;
    advance_to(0.0);
    record_sample();
}

// Set the two spin poses for sim time t (every other frame co-rotates with its parent).
void active_grinding_marks::advance_to(double t)
{
    double const l3 = 0.3 * m_params.R;
    m_sys.set_pose(m_sys.index_of("chuck_ctr_rot"),
                   pose3dp{O_3dp, vec3dp{m_params.ns * t, 0.0, 0.0, 0.0}});
    m_sys.set_pose(
        m_sys.index_of("tool_top_avg_rot"),
        pose3dp{vec3dp{0.0, 0.0, l3, 1.0}, vec3dp{0.0, 0.0, m_params.nw * t, 0.0}});
}

// Append the rim grain's position IN THE WAFER FRAME (the grinding mark).
void active_grinding_marks::record_sample()
{
    vec3dp const g = unitize(move3dp(
        O_3dp, m_sys.get_pos_trafo("tool_surface_avg_at_R", "wafer_top_avg_rot")));
    m_path.push_back(g);
    if (m_path.size() > MAX_PATH) m_path.pop_front();
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_grinding_marks::tick()
{
    if (m_paused) return;

    m_t += m_params.dt;
    advance_to(m_t);
    record_sample();

    prepareGeometryChange();
    update();
}

void active_grinding_marks::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_grinding_marks::toggleTrace()
{
    m_show_trace = !m_show_trace;
    update();
}

void active_grinding_marks::resetAnimation()
{
    build_system();
    prepareGeometryChange();
    update();
}

void active_grinding_marks::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Coordinate helpers
// ---------------------------------------------------------------------------

// Project a 3dp point onto the chuck plane e423_3dp (drop e1) and map to canvas:
// e2 (root y) -> horizontal, e3 (root z) -> vertical.
QPointF active_grinding_marks::toScreen(vec3dp const& p) const
{
    vec3dp const pu = unitize(p);
    return QPointF(cs->x.au_to_w(pu.y), cs->y.au_to_w(pu.z));
}

void active_grinding_marks::drawCircle(QPainter* qp, vec3dp const& c, double r,
                                       QColor const& edge, QColor const& fill)
{
    double const px = std::abs(cs->x.px_density_rng());
    double const py = std::abs(cs->y.px_density_rng());
    qp->setPen(QPen(edge, 1.5));
    qp->setBrush(QBrush(fill));
    qp->drawEllipse(toScreen(c), r * px, r * py);
}

// ---------------------------------------------------------------------------
// Painting
// ---------------------------------------------------------------------------

void active_grinding_marks::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                                  QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    double const R = m_params.R;

    QFont lbl_font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl_font.setPointSize(9);
    lbl_font.setBold(true);

    // --- wafer: fixed unit disk at the wafer-frame origin (the wafer centre) ---
    drawCircle(qp, O_3dp, R, QColor(50, 100, 170), QColor(180, 205, 235, 110));

    // wafer orientation marker: a thin reference spoke centre -> +e2 rim (fixed in the
    // wafer frame, so its slow turn shows the chuck rotation)
    qp->setPen(QPen(QColor(50, 100, 170, 160), 1.2, Qt::SolidLine));
    qp->drawLine(toScreen(O_3dp), toScreen(vec3dp{0.0, R, 0.0, 1.0}));

    // --- grinding-mark trace of the rim grain (drawn under the live geometry) ---
    if (m_show_trace && m_path.size() >= 2) {
        qp->setPen(QPen(QColor(200, 60, 0, 200), 1.6, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_path.size(); ++i)
            qp->drawLine(toScreen(m_path[i - 1]), toScreen(m_path[i]));
    }

    // --- wheel: current disk (its centre, fixed in the root frame, orbits in the wafer
    // frame) + spoke to the grain ---
    vec3dp const wheel_c = unitize(
        move3dp(O_3dp, m_sys.get_pos_trafo("tool_top_avg_rot", "wafer_top_avg_rot")));
    drawCircle(qp, wheel_c, R, QColor(150, 90, 0), QColor(245, 200, 130, 70));
    qp->setPen(QPen(QColor(150, 90, 0), 2));
    qp->setBrush(QBrush(QColor(150, 90, 0)));
    qp->drawEllipse(toScreen(wheel_c), 3.0, 3.0);

    vec3dp const grain = m_path.empty() ? O_3dp : m_path.back();
    qp->setPen(QPen(QColor(150, 90, 0, 170), 1.4, Qt::DashLine));
    qp->drawLine(toScreen(wheel_c), toScreen(grain)); // wheel radius to the grain

    // --- the grain itself ---
    qp->setPen(QPen(QColor(200, 60, 0), 2));
    qp->setBrush(QBrush(QColor(230, 90, 20)));
    qp->drawEllipse(toScreen(grain), 4.0, 4.0);

    // --- wafer centre dot + labels ---
    qp->setPen(QPen(Qt::black, 2));
    qp->setBrush(QBrush(Qt::black));
    qp->drawEllipse(toScreen(O_3dp), 2.5, 2.5);

    qp->setFont(lbl_font);
    qp->setPen(QPen(QColor(30, 60, 110), 1));
    qp->drawText(toScreen(vec3dp{0.0, R, 0.0, 1.0}) + QPointF(6.0, -6.0), "wafer");
    qp->setPen(QPen(QColor(150, 90, 0), 1));
    qp->drawText(toScreen(wheel_c) + QPointF(6.0, -6.0), "wheel");

    qp->restore();
}

QRectF active_grinding_marks::boundingRect() const
{
    double const reach = 2.4 * m_params.R; // grain sweeps r in [0, 2 R]; small margin
    return QRectF(toScreen(vec3dp{0.0, -reach, reach, 1.0}),
                  toScreen(vec3dp{0.0, reach, -reach, 1.0}));
}

QPainterPath active_grinding_marks::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
