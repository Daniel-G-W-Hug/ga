// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_grinding_flatness.hpp"

#include "ga/ga_pga.hpp"

#include <QTimer>
#include <array>
#include <cmath>

using namespace hd::ga;
using namespace hd::ga::pga;

namespace {

constexpr int UPDATE_MS = 16;

// panels stacked on the right (legend reserved on the left): the X-Y rosette above its
// two profile cross-sections X-Z and Y-Z. All share the spatial scale S_XY for the
// in-plane coords; the height z is drawn on the hugely exaggerated S_Z (z ~ 0.2 mm vs a
// 300 mm diameter), exactly as Zhou plots it.
constexpr double CX = 1.15;      // shared horizontal centre of the panels [au]
constexpr double CY_XY = 1.25;   // X-Y rosette centre [au]
constexpr double CY_XZ = -0.45;  // X-Z profile centre [au]
constexpr double CY_YZ = -1.65;  // Y-Z profile centre [au]
constexpr double HW = 1.05;      // panel half-width [au]
constexpr double HH_XY = 1.05;   // rosette half-height [au]
constexpr double HH_PROF = 0.55; // profile half-height [au]
constexpr double S_XY =
    0.9 / 150.0;            // au per mm in the wafer plane (R = 150 mm -> 0.9 au)
constexpr double S_Z = 1.4; // au per mm of height (exaggerated)

// Fig. 7's 3x3 grid of (alpha, beta) in units of the tilt magnitude, row-major from the
// top-left: alpha = +,+,+, 0,0,0, -,-,- ; beta = -,0,+ repeating. Index 4 = (0,0) flat.
constexpr std::array<std::pair<int, int>, 9> GRID = {
    {{+1, -1}, {+1, 0}, {+1, +1}, {0, -1}, {0, 0}, {0, +1}, {-1, -1}, {-1, 0}, {-1, +1}}};

} // namespace

active_grinding_flatness::active_grinding_flatness(Coordsys* cs, w_Coordsys* wcs,
                                                   grinding_flatness_params const& params,
                                                   QGraphicsItem* parent) :
    QGraphicsItem(parent), m_params{params}, cs{cs}, wcs{wcs}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    recompute();

    connect(wcs, &w_Coordsys::viewResized, this, &active_grinding_flatness::viewChanged);
    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_grinding_flatness::tick);
    m_timer->start(UPDATE_MS);
}

// ---------------------------------------------------------------------------
// Cutting-path computation (Zhou Eq.1 in PGA, the app-test geometry)
// ---------------------------------------------------------------------------

void active_grinding_flatness::recompute()
{
    double const R = m_params.R, L = m_params.L, r1 = m_params.r1;
    double const w1 = rpm2radps(m_params.n1), w2 = rpm2radps(m_params.n2);
    double const a = double(GRID[m_tilt].first) * m_params.tilt * pi / 180.0;  // alpha
    double const b = double(GRID[m_tilt].second) * m_params.tilt * pi / 180.0; // beta

    kinematic_system3dp sys;
    sys.add_frame(static_frame3dp("W"));
    sys.add_frame(static_frame3dp("wheel_rot", vec3dp{0.0, L, 0.0, 1.0}),
                  kin_state3dp{.omega = vec3dp{w1, 0.0, 0.0, 0.0}}, sys.index_of("W"));
    sys.add_frame(static_frame3dp("grain", vec3dp{0.0, -r1, 0.0, 1.0}), kin_state3dp{},
                  sys.index_of("wheel_rot"));
    sys.add_frame(static_frame3dp("wafer_tilt", O_3dp, vec3dp{0.0, a, b, 0.0}),
                  kin_state3dp{}, sys.index_of("W"));
    sys.add_frame(static_frame3dp("wafer_rot"),
                  kin_state3dp{.omega = vec3dp{w2, 0.0, 0.0, 0.0}},
                  sys.index_of("wafer_tilt"));

    size_t const wheel = sys.index_of("wheel_rot"), wafer = sys.index_of("wafer_rot");
    double const T = 2.0 * pi / w2; // one wafer revolution
    int const N = m_params.samples;

    m_path.clear();
    m_path.reserve(N + 1);
    bool prev_on = false;
    for (int k = 0; k <= N; ++k) {
        double const t = T * double(k) / double(N);
        sys.set_pose(wheel,
                     pose3dp{vec3dp{0.0, L, 0.0, 1.0}, vec3dp{w1 * t, 0.0, 0.0, 0.0}});
        sys.set_pose(wafer, pose3dp{O_3dp, vec3dp{w2 * t, 0.0, 0.0, 0.0}});
        auto const g = unitize(move3dp(O_3dp, sys.get_pos_trafo("grain", "wafer_rot")));
        double const r = std::sqrt(g.y * g.y + g.z * g.z);
        bool const on = r <= R + 1e-6; // only mark where the grain is on the wafer
        if (on) {
            // (e2, e3) = wafer-plane (x, y); e1 = profile height z
            m_path.push_back({g.y, g.z, g.x, /*pen=*/prev_on});
        }
        prev_on = on;
    }
    m_shown = 0;
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_grinding_flatness::tick()
{
    if (m_paused || m_path.empty()) return;
    if (m_shown < m_path.size()) {
        m_shown = std::min(m_path.size(), m_shown + m_path.size() / 240 + 1);
        update();
    }
}

void active_grinding_flatness::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_grinding_flatness::resetAnimation()
{
    m_shown = 0;
    update();
}

void active_grinding_flatness::cycleTilt()
{
    m_tilt = (m_tilt + 1) % 9;
    recompute();
    update();
}

void active_grinding_flatness::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Projections + drawing
// ---------------------------------------------------------------------------

QPointF active_grinding_flatness::toXY(path_pt const& p) const
{
    return QPointF(cs->x.au_to_w(CX + p.x * S_XY), cs->y.au_to_w(CY_XY + p.y * S_XY));
}
QPointF active_grinding_flatness::toXZ(path_pt const& p) const
{
    return QPointF(cs->x.au_to_w(CX + p.x * S_XY), cs->y.au_to_w(CY_XZ + p.z * S_Z));
}
QPointF active_grinding_flatness::toYZ(path_pt const& p) const
{
    return QPointF(cs->x.au_to_w(CX + p.y * S_XY), cs->y.au_to_w(CY_YZ + p.z * S_Z));
}

void active_grinding_flatness::drawPath(
    QPainter* qp, QPointF (active_grinding_flatness::*proj)(path_pt const&) const) const
{
    qp->setPen(QPen(QColor(150, 30, 30), 1.0)); // cutting path: dark red
    size_t const n = std::min(m_shown, m_path.size());
    for (size_t i = 1; i < n; ++i) {
        if (m_path[i].pen)
            qp->drawLine((this->*proj)(m_path[i - 1]), (this->*proj)(m_path[i]));
    }
}

void active_grinding_flatness::drawPanel(QPainter* qp, double cx, double cy, double hw,
                                         double hh, char const* title) const
{
    QPointF const tl(cs->x.au_to_w(cx - hw), cs->y.au_to_w(cy + hh));
    QPointF const br(cs->x.au_to_w(cx + hw), cs->y.au_to_w(cy - hh));
    qp->setPen(QPen(QColor(200, 200, 200), 1.0));
    qp->drawRect(QRectF(tl, br));
    // a faint z = 0 reference line through the profile panels
    if (std::abs(cy - CY_XY) > 0.01) {
        qp->setPen(QPen(QColor(210, 210, 210), 1.0, Qt::DashLine));
        qp->drawLine(QPointF(tl.x(), cs->y.au_to_w(cy)),
                     QPointF(br.x(), cs->y.au_to_w(cy)));
    }
    qp->setPen(QPen(Qt::black, 1.0));
    qp->drawText(QPointF(tl.x() + 6.0, tl.y() + 16.0), title);
}

void active_grinding_flatness::paint(QPainter* qp, QStyleOptionGraphicsItem const*,
                                     QWidget*)
{
    qp->setRenderHint(QPainter::Antialiasing, true);
    qp->setBrush(Qt::NoBrush);

    drawPanel(qp, CX, CY_XY, HW, HH_XY, "X-Y cutting path (rosette)");
    drawPanel(qp, CX, CY_XZ, HW, HH_PROF, "X-Z profile (height exaggerated)");
    drawPanel(qp, CX, CY_YZ, HW, HH_PROF, "Y-Z profile (height exaggerated)");

    drawPath(qp, &active_grinding_flatness::toXY);
    drawPath(qp, &active_grinding_flatness::toXZ);
    drawPath(qp, &active_grinding_flatness::toYZ);

    // --- scale key (lower-left): the radial dimension vs the (exaggerated) height, so
    // the two very different scales are legible. The 150 mm radial bar drawn next to the
    // 0.2 mm height bar makes the exaggeration self-evident.
    auto hbar = [&](double x, double y, double mm, char const* lbl) {
        double const x0 = cs->x.au_to_w(x), x1 = cs->x.au_to_w(x + mm * S_XY);
        double const yw = cs->y.au_to_w(y);
        qp->setPen(QPen(Qt::black, 1.4));
        qp->drawLine(QPointF(x0, yw), QPointF(x1, yw));
        qp->drawLine(QPointF(x0, yw - 4.0), QPointF(x0, yw + 4.0));
        qp->drawLine(QPointF(x1, yw - 4.0), QPointF(x1, yw + 4.0));
        qp->drawText(QPointF(x1 + 8.0, yw + 4.0), lbl);
    };
    auto vbar = [&](double x, double y, double mm, char const* lbl) {
        double const xw = cs->x.au_to_w(x);
        double const y0 = cs->y.au_to_w(y), y1 = cs->y.au_to_w(y + mm * S_Z);
        qp->setPen(QPen(Qt::black, 1.4));
        qp->drawLine(QPointF(xw, y0), QPointF(xw, y1));
        qp->drawLine(QPointF(xw - 4.0, y0), QPointF(xw + 4.0, y0));
        qp->drawLine(QPointF(xw - 4.0, y1), QPointF(xw + 4.0, y1));
        qp->drawText(QPointF(xw + 8.0, 0.5 * (y0 + y1) + 4.0), lbl);
    };
    qp->setPen(QPen(Qt::black, 1.0));
    qp->drawText(QPointF(cs->x.au_to_w(-3.3), cs->y.au_to_w(-0.95)), "Scale:");
    hbar(-3.2, -1.2, 150.0, "150 mm  (radial X / Y)");
    vbar(-3.2, -1.85, 0.2, "0.2 mm  (height z)");
    qp->drawText(QPointF(cs->x.au_to_w(-3.3), cs->y.au_to_w(-2.0)),
                 QString::asprintf("height drawn ~%.0fx exaggerated", S_Z / S_XY));

    // (alpha, beta) + shape readout
    double const a = double(GRID[m_tilt].first) * m_params.tilt;
    double const b = double(GRID[m_tilt].second) * m_params.tilt;
    char const* shape = (a == 0.0 && b == 0.0) ? "flat"
                        : (b == 0.0)           ? "cone"
                        : (a == 0.0 && b > 0)  ? "convex (dome)"
                        : (a == 0.0 && b < 0)  ? "concave (bowl)"
                                               : "combined";
    qp->setPen(QPen(Qt::black, 1.0));
    qp->drawText(
        QPointF(cs->x.au_to_w(-3.3), cs->y.au_to_w(-2.25)),
        QString::asprintf("alpha = %+.1f deg, beta = %+.1f deg  ->  %s", a, b, shape));
    qp->drawText(QPointF(cs->x.au_to_w(-3.3), cs->y.au_to_w(-2.45)),
                 "C: step the (alpha, beta) grid (Fig. 7)");
}

QRectF active_grinding_flatness::boundingRect() const
{
    return QRectF(QPointF(cs->x.au_to_w(-3.5), cs->y.au_to_w(2.5)),
                  QPointF(cs->x.au_to_w(3.5), cs->y.au_to_w(-2.5)));
}

QPainterPath active_grinding_flatness::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
