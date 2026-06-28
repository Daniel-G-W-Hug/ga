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

double active_grinding_flatness::lambda_m() const
{
    // waviness wavelength along the mark direction = wheel surface speed / f_b [mm].
    // v_w = omega_w * R [mm/s] (the relative grinding speed that sweeps the ripple).
    double const v_w = rpm2radps(m_params.n1) * m_params.R; // [mm/s]
    return v_w / m_params.f_b;                              // [mm]
}

// The dynamic-mode bottom panel: the spindle-vibration waviness z_b(s) = A_b sin(2 pi s /
// lambda_m) along the mark direction. The waviness is ~0.1 um against a ~0.2 mm flatness,
// so it gets its OWN large exaggeration here (it would be invisible in the flatness
// panels). About four wavelengths are shown, with the amplitude A_b and the wavelength
// lambda_m each quantified by a labelled dimension bar so both are readable off the
// panel.
void active_grinding_flatness::drawWaviness(QPainter* qp, double cy, double hh) const
{
    double const lam = lambda_m();              // [mm]
    double const nwave = 4.0;                   // wavelengths shown across the panel
    double const span = nwave * lam;            // mark-direction extent shown [mm]
    double const Ab_mm = m_params.A_b * 1.0e-3; // waviness amplitude [mm]
    double const sx = 2.0 * HW / span;          // au per mm along the mark
    double const amp = 0.40 * hh;               // curve half-height in the panel [au]
    double const sz = amp / std::max(Ab_mm, 1.0e-9); // au per mm of z_b (own exag.)
    double const xL = CX - HW;                       // panel left edge [au]

    auto X = [&](double mm) { return cs->x.au_to_w(xL + mm * sx); };
    auto Z = [&](double zmm) { return cs->y.au_to_w(cy + zmm * sz); };

    // dashed z_b = 0 midline (the cutting-depth reference the ripple modulates)
    qp->setPen(QPen(QColor(205, 205, 205), 1.0, Qt::DashLine));
    qp->drawLine(QPointF(X(0.0), cs->y.au_to_w(cy)), QPointF(X(span), cs->y.au_to_w(cy)));

    // the waviness curve, growing with the trace build-up
    size_t const n = std::min(m_shown, m_path.size());
    double const frac = m_path.empty() ? 0.0 : double(n) / double(m_path.size());
    qp->setPen(QPen(QColor(30, 90, 170), 1.5)); // waviness: blue (distinct from flatness)
    QPointF prev;
    bool have_prev = false;
    int const NW = 600;
    for (int i = 0; i <= NW; ++i) {
        double const u = double(i) / double(NW); // 0..1 across the panel
        if (u > frac) break;                     // grow with the trace build-up
        double const s = u * span;               // arc length along the mark [mm]
        QPointF const pt(X(s), Z(Ab_mm * std::sin(2.0 * pi * s / lam)));
        if (have_prev) qp->drawLine(prev, pt);
        prev = pt;
        have_prev = true;
    }

    // --- amplitude callout: a vertical 0 -> A_b bar at the first crest, the value placed
    // ABOVE the crest tip (in the clear strip between the crests and the panel top --
    // well below the title row, and clear of the descending curve to either side).
    double const xa = X(0.25 * lam);
    double const ya0 = cs->y.au_to_w(cy), yaA = Z(Ab_mm);
    qp->setPen(QPen(Qt::black, 1.3));
    qp->drawLine(QPointF(xa, ya0), QPointF(xa, yaA));
    qp->drawLine(QPointF(xa - 4.0, yaA), QPointF(xa + 4.0, yaA)); // crest tick
    qp->drawLine(QPointF(xa - 4.0, ya0), QPointF(xa + 4.0, ya0)); // midline tick
    qp->drawText(QPointF(xa - 26.0, yaA - 7.0),
                 QString::asprintf("A_b = %.2f um", m_params.A_b));

    // --- wavelength callout: a horizontal one-period bar between the two MIDDLE crests
    // (centred in the panel), below the troughs, with the lambda_m label centred under it
    // over the FULL panel width so the text never clips (the lambda = v_w / f_b relation
    // is spelled out in the legend).
    double const yw = cs->y.au_to_w(cy - 0.70 * hh);
    double const xw0 = X(1.25 * lam),
                 xw1 = X(2.25 * lam); // middle crest-to-crest = one lam
    qp->setPen(QPen(Qt::black, 1.3));
    qp->drawLine(QPointF(xw0, yw), QPointF(xw1, yw));
    qp->drawLine(QPointF(xw0, yw - 4.0), QPointF(xw0, yw + 4.0));
    qp->drawLine(QPointF(xw1, yw - 4.0), QPointF(xw1, yw + 4.0));
    qp->drawText(QRectF(QPointF(X(0.0), yw + 6.0), QPointF(X(span), yw + 24.0)),
                 Qt::AlignHCenter | Qt::AlignTop,
                 QString::asprintf("lambda_m = %.2f mm", lam));
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

    bool const dyn = m_params.dynamic;

    drawPanel(qp, CX, CY_XY, HW, HH_XY, "X-Y cutting path (rosette)");
    drawPanel(qp, CX, CY_XZ, HW, HH_PROF,
              dyn ? "X-Z flatness profile (mean tilt; exaggerated)"
                  : "X-Z profile (height exaggerated)");
    drawPath(qp, &active_grinding_flatness::toXY);
    drawPath(qp, &active_grinding_flatness::toXZ);

    if (dyn) {
        // dynamic view: the bottom panel is the spindle-vibration WAVINESS (z_b ripple),
        // not the Y-Z flatness cross-section -- the new, second error source. f_b rides
        // in the title (keeps the in-panel callouts short enough not to clip).
        QByteArray const wtitle =
            QString::asprintf("z_b waviness, f_b = %.0f Hz (Tao Dyn.2)", m_params.f_b)
                .toLatin1();
        drawPanel(qp, CX, CY_YZ, HW, HH_PROF, wtitle.constData());
        drawWaviness(qp, CY_YZ, HH_PROF);
    }
    else {
        drawPanel(qp, CX, CY_YZ, HW, HH_PROF, "Y-Z profile (height exaggerated)");
        drawPath(qp, &active_grinding_flatness::toYZ);
    }

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
