// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_grinding_cs.hpp"

#include "ga/ga_pga.hpp"

#include <QTimer>
#include <cmath>

using namespace hd::ga;
using namespace hd::ga::pga;

namespace {

constexpr int UPDATE_MS = 16; // animation tick [ms]
constexpr double S = 0.55;    // au per normalized unit (panel scale)

// The two panels are STACKED on the right (the left of the view is reserved for the
// legend) and share the e3 horizontal axis, so a feature lines up vertically between the
// top view (above) and the side view (below) -- engineering-drawing convention.
constexpr double CX = 1.35;       // shared horizontal centre of both panels [au]
constexpr double CY_TOP = 1.25;   // top-view panel centre (upper) [au]
constexpr double CY_SIDE = -1.25; // side-view panel centre (lower) [au]
constexpr double PANEL_HW = 1.1;  // panel half-width  [au]
constexpr double PANEL_HH = 1.05; // panel half-height [au]

// content centres, so the wafer + offset wheel sit centred in each panel
constexpr double E3C = 0.5; // horizontal (e3): wafer [-1,1] + wheel [0,2] -> centre 0.5
constexpr double E2C = 0.0; // top-view vertical (e2)
constexpr double E1C = 0.3; // side-view vertical (e1)

constexpr double TWO_PI = 2.0 * pi;

// nominal placement, normalized to R = wheel = wafer = 1: the wheel centre is offset by R
// along +e3 (to the RIGHT of the chuck, so the rim passes through the chuck centre), the
// tool face lifted by X_A along the chuck axis e1, reoriented -pi/2 about e2 so the tool
// axis (local z) points down toward the wafer.
constexpr double OFF_R = 1.0; // wheel-centre offset along +e3
constexpr double X_A = 0.30;  // tool-face lift along the chuck axis e1

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_grinding_cs::active_grinding_cs(Coordsys* cs, w_Coordsys* wcs,
                                       grinding_cs_params const& params,
                                       QGraphicsItem* parent) :
    QGraphicsItem(parent), m_params{params}, cs{cs}, wcs{wcs}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    connect(wcs, &w_Coordsys::viewResized, this, &active_grinding_cs::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_grinding_cs::tick);
    m_timer->start(UPDATE_MS);
}

// ---------------------------------------------------------------------------
// The spindle placement motor (perfect ghost, or equilibrium + prescribed wobble)
// ---------------------------------------------------------------------------

mvec3dp_e active_grinding_cs::spindle_motor(bool ideal) const
{
    double th = 0.0, ph = 0.0;                    // tool-axis tilt about e2, e3
    double rad_y = 0.0, rad_z = 0.0, axial = 0.0; // CM offsets (e2, e3, e1)

    if (!ideal) {
        // the imperfect EQUILIBRIUM: a small static inclination (toggleable with 'T')
        if (m_incl_on) {
            th = m_params.incl_e2;
            ph = m_params.incl_e3;
        }
        // the prescribed illustrative wobble of the active DOF group, about equilibrium
        bool const r = m_dof == dof_group::radial || m_dof == dof_group::all;
        bool const a = m_dof == dof_group::axial || m_dof == dof_group::all;
        bool const t = m_dof == dof_group::tilt || m_dof == dof_group::all;
        double const t_ = m_t;
        if (r) {
            rad_y = m_params.radial_amp * std::sin(TWO_PI * m_params.f_radial * t_);
            rad_z = m_params.radial_amp * std::cos(TWO_PI * m_params.f_radial * t_);
        }
        if (a) axial = m_params.axial_amp * std::sin(TWO_PI * m_params.f_axial * t_);
        if (t) {
            th += m_params.tilt_amp * std::sin(TWO_PI * m_params.f_tilt * t_);
            ph += m_params.tilt_amp * std::cos(TWO_PI * m_params.f_tilt * t_);
        }
    }

    // perturbed origin: nominal +e3 offset + (axial along e1, radial in the e2,e3 plane)
    vec3dp const O_w{X_A + axial, rad_y, OFF_R + rad_z, 1.0};

    // M = T(O_w)  <rgpr>  Rtilt(about e2,e3)  <rgpr>  R0(-pi/2 about e2)
    // applied right-to-left: nominal reorient, then tip, then translate into place.
    mvec3dp_e const R0 =
        motor_from_pose3dp(pose3dp{O_3dp, vec3dp{0.0, -pi / 2.0, 0.0, 0.0}});
    mvec3dp_e const Rtilt = motor_from_pose3dp(pose3dp{O_3dp, vec3dp{0.0, th, ph, 0.0}});
    mvec3dp_e const T = motor_from_pose3dp(pose3dp{O_w, vec3dp{0.0, 0.0, 0.0, 0.0}});
    return rgpr(T, rgpr(Rtilt, R0));
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_grinding_cs::tick()
{
    if (m_paused) return;
    m_t += m_params.dt;
    prepareGeometryChange();
    update();
}

void active_grinding_cs::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_grinding_cs::resetAnimation()
{
    m_t = 0.0;
    prepareGeometryChange();
    update();
}

// C: radial -> axial -> tilt -> all, so the simple single-DOF cases come first and the
// combined case is shown last.
void active_grinding_cs::cycleDof()
{
    switch (m_dof) {
        case dof_group::radial:
            m_dof = dof_group::axial;
            break;
        case dof_group::axial:
            m_dof = dof_group::tilt;
            break;
        case dof_group::tilt:
            m_dof = dof_group::all;
            break;
        case dof_group::all:
            m_dof = dof_group::radial;
            break;
    }
    update();
}

void active_grinding_cs::toggleIncl()
{
    m_incl_on = !m_incl_on;
    update();
}

void active_grinding_cs::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Projections (the two stacked panels; both share the e3 horizontal axis)
// ---------------------------------------------------------------------------

QPointF active_grinding_cs::toTop(vec3dp const& p) const
{
    vec3dp const pu = unitize(p); // e3 -> horizontal, e2 -> vertical
    return QPointF(cs->x.au_to_w(CX + (pu.z - E3C) * S),
                   cs->y.au_to_w(CY_TOP + (pu.y - E2C) * S));
}

QPointF active_grinding_cs::toSide(vec3dp const& p) const
{
    vec3dp const pu = unitize(p); // e3 -> horizontal, e1 (chuck axis) -> vertical
    return QPointF(cs->x.au_to_w(CX + (pu.z - E3C) * S),
                   cs->y.au_to_w(CY_SIDE + (pu.x - E1C) * S));
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

// The wafer disk lies in the chuck plane (perp e1) at e1 == 0. In the top view it is a
// circle (face-on); in the side view it is edge-on -> a horizontal segment.
void active_grinding_cs::drawWafer(QPainter* qp, bool side) const
{
    qp->setPen(QPen(QColor(70, 110, 200), 2.0));
    double const R = m_params.R;
    if (!side) {
        QPolygonF poly;
        for (int i = 0; i <= 64; ++i) {
            double const a = TWO_PI * i / 64.0;
            poly << toTop(vec3dp{0.0, R * std::cos(a), R * std::sin(a), 1.0});
        }
        qp->drawPolyline(poly);
    }
    else {
        qp->drawLine(toSide(vec3dp{0.0, 0.0, -R, 1.0}), toSide(vec3dp{0.0, 0.0, R, 1.0}));
    }
}

// One spindle rendering through the supplied panel projection: the cup-wheel rim circle
// (in the tool's local x-y plane), the tool / spin axis (local z) and the origin marker.
// `ghost` draws the ideal square placement in grey dashed.
void active_grinding_cs::drawSpindle(QPainter* qp,
                                     QPointF (active_grinding_cs::*proj)(vec3dp const&)
                                         const,
                                     bool ghost) const
{
    mvec3dp_e const M = spindle_motor(ghost);
    double const R = m_params.R;

    QPolygonF rim;
    for (int i = 0; i <= 64; ++i) {
        double const a = TWO_PI * i / 64.0;
        vec3dp const w = move3dp(vec3dp{R * std::cos(a), R * std::sin(a), 0.0, 1.0}, M);
        rim << (this->*proj)(w);
    }

    QPointF const ax0 = (this->*proj)(move3dp(vec3dp{0.0, 0.0, -0.7 * R, 1.0}, M));
    QPointF const ax1 = (this->*proj)(move3dp(vec3dp{0.0, 0.0, 0.35 * R, 1.0}, M));
    QPointF const O = (this->*proj)(move3dp(O_3dp, M));

    if (ghost) {
        qp->setPen(QPen(QColor(150, 150, 150, 170), 1.5, Qt::DashLine));
        qp->drawPolyline(rim);
        qp->drawLine(ax0, ax1);
    }
    else {
        qp->setPen(QPen(QColor(230, 120, 30), 2.4)); // cup wheel: orange
        qp->drawPolyline(rim);
        qp->setPen(QPen(QColor(150, 30, 30), 2.8)); // tool / spin axis: dark red
        qp->drawLine(ax0, ax1);
        qp->setBrush(QColor(150, 30, 30));
        qp->drawEllipse(O, 3.0, 3.0);
        qp->setBrush(Qt::NoBrush);
    }
}

// panel chrome: a light bounding box and a title at the top-left of the box.
void active_grinding_cs::drawPanelFrame(QPainter* qp, double cx, double cy,
                                        char const* title) const
{
    QPointF const tl(cs->x.au_to_w(cx - PANEL_HW), cs->y.au_to_w(cy + PANEL_HH));
    QPointF const br(cs->x.au_to_w(cx + PANEL_HW), cs->y.au_to_w(cy - PANEL_HH));
    qp->setPen(QPen(QColor(200, 200, 200), 1.0));
    qp->drawRect(QRectF(tl, br));
    qp->setPen(QPen(Qt::black, 1.0));
    qp->drawText(QPointF(tl.x() + 6.0, tl.y() + 16.0), title);
}

void active_grinding_cs::paint(QPainter* qp, QStyleOptionGraphicsItem const*, QWidget*)
{
    qp->setRenderHint(QPainter::Antialiasing, true);
    qp->setBrush(Qt::NoBrush);

    drawPanelFrame(qp, CX, CY_TOP, "top view (down -e1): radial");
    drawPanelFrame(qp, CX, CY_SIDE, "side view (along -e2): axial + tilt");

    // top panel
    drawWafer(qp, /*side=*/false);
    drawSpindle(qp, &active_grinding_cs::toTop, /*ghost=*/true);
    drawSpindle(qp, &active_grinding_cs::toTop, /*ghost=*/false);

    // side panel
    drawWafer(qp, /*side=*/true);
    drawSpindle(qp, &active_grinding_cs::toSide, /*ghost=*/true);
    drawSpindle(qp, &active_grinding_cs::toSide, /*ghost=*/false);

    // status readout (lower-left, under the legend)
    char const* mode = m_dof == dof_group::radial  ? "radial"
                       : m_dof == dof_group::axial ? "axial"
                       : m_dof == dof_group::tilt  ? "tilt"
                                                   : "all (radial + axial + tilt)";
    qp->setPen(QPen(Qt::black, 1.0));
    qp->drawText(QPointF(cs->x.au_to_w(-3.3), cs->y.au_to_w(-2.1)),
                 QString("driven DOF (C): ") + mode);
    qp->drawText(QPointF(cs->x.au_to_w(-3.3), cs->y.au_to_w(-2.35)),
                 QString("equilibrium (T): ") +
                     (m_incl_on ? "inclined (imperfect)" : "ideal square"));
}

QRectF active_grinding_cs::boundingRect() const
{
    return QRectF(QPointF(cs->x.au_to_w(-3.5), cs->y.au_to_w(2.5)),
                  QPointF(cs->x.au_to_w(3.5), cs->y.au_to_w(-2.5)));
}

QPainterPath active_grinding_cs::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
