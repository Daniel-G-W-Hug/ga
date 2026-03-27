// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_common.hpp"
#include "active_frame_trafo.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <QPolygonF>
#include <QTimer>
#include <algorithm>
#include <cmath>
#include <fmt/core.h>

using namespace hd::ga;
using namespace hd::ga::pga;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_frame_trafo::active_frame_trafo(Coordsys* cs, w_Coordsys* wcs,
                                       frame_trafo_params const& params,
                                       QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    // Initialise state: both accumulators start at B_b_init so they measure
    // from the same reference pose and are directly comparable.
    // For pure translation Omega_w = Omega_b, so B_w ≡ B_b at all times.
    // For rotation they diverge, making the frame-dependence visible.
    m_B_b = m_params.B_b_init;
    m_B_w = m_params.B_b_init;
    m_M   = rgpr(m_params.M0, exp(0.5 * m_B_b));

    connect(wcs, &w_Coordsys::viewResized, this, &active_frame_trafo::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_frame_trafo::tick);
    m_timer->start(UPDATE_MS);
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_frame_trafo::tick()
{
    if (m_paused) return;

    if (m_params.world_frame_drive) {
        // Analytical M(t) = T(O(t)) ⟇ R(cm_omega*t)
        // O(t) moves in a straight line → exact CM path regardless of rotation.
        double const ox    = m_params.cm_ox + m_params.cm_vx * m_t;
        double const oy    = m_params.cm_oy + m_params.cm_vy * m_t;
        double const theta = m_params.cm_omega * m_t;
        mvec2dp_u const R_t = exp(0.5 * vec2dp{0.0, 0.0, theta});
        mvec2dp_u const T_t = exp(0.5 * vec2dp{-oy, ox, 0.0});
        m_M = rgpr(T_t, R_t);

        // Instantaneous world-frame generator: Momentanpol at IC = O + (-vy/ω, vx/ω)
        // L_w = (ω*ox - vy)*e31 + (ω*oy + vx)*e32 + ω*e12
        double const o_z = m_params.cm_omega;
        vec2dp const Omega_w{o_z * ox - m_params.cm_vy,
                              o_z * oy + m_params.cm_vx,
                              o_z};

        // Body-frame generator = inverse sandwich of world-frame generator
        vec2dp const Omega_b_inst = move2dp(Omega_w, rrev(m_M));

        // Euler-integrate both accumulators for informational text display
        m_B_w = m_B_w + Omega_w * DT;
        m_B_b = m_B_b + Omega_b_inst * DT;

    } else {
        // Body-frame driven: constant Omega_b, M from accumulated B_b
        m_M = rgpr(m_params.M0, exp(0.5 * m_B_b));

        // World-frame velocity: sandwich Omega_b with M
        vec2dp const Omega_w = move2dp(m_params.Omega_b, m_M);

        // Euler integration of both accumulators
        m_B_b = m_B_b + m_params.Omega_b * DT;
        m_B_w = m_B_w + Omega_w * DT;
    }

    m_t += DT;

    // Record body-origin trajectory
    vec2dp const O_b{0.0, 0.0, 1.0};
    vec2dp const O_w = move2dp(O_b, m_M);
    m_trajectory.push_back({O_w});
    if (m_trajectory.size() > MAX_TRAJ) m_trajectory.pop_front();

    if (m_t >= m_params.duration) resetAnimation(); // loop

    prepareGeometryChange();
    update();
}

void active_frame_trafo::resetAnimation()
{
    m_t   = 0.0;
    m_B_b = m_params.B_b_init;
    m_B_w = m_params.B_b_init;

    if (m_params.world_frame_drive) {
        mvec2dp_u const T_0 = exp(0.5 * vec2dp{-m_params.cm_oy, m_params.cm_ox, 0.0});
        mvec2dp_u const R_0 = exp(0.5 * vec2dp{0.0, 0.0, 0.0}); // identity rotation
        m_M = rgpr(T_0, R_0);
    } else {
        m_M = rgpr(m_params.M0, exp(0.5 * m_B_b));
    }

    m_trajectory.clear();
    update();
}

void active_frame_trafo::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_frame_trafo::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Coordinate helper
// ---------------------------------------------------------------------------

QPointF active_frame_trafo::toScreen(vec2dp const& pt) const
{
    return QPointF(cs->x.au_to_w(pt.x), cs->y.au_to_w(pt.y));
}

// ---------------------------------------------------------------------------
// Drawing helpers
// ---------------------------------------------------------------------------

// World frame W: two fixed reference arrows at (wx, wy).
// Skipped entirely when m_params.draw_world_frame == false.
void active_frame_trafo::drawWorldFrame(QPainter* qp, vec2dp const& /*O_w*/,
                                        vec2dp const& /*e1_w*/, vec2dp const& /*e2_w*/)
{
    if (!m_params.draw_world_frame) return;
    vec2dp const W_orig{m_params.wx, m_params.wy, 1.0};
    vec2dp const W_e1_tip{m_params.wx + AXIS_LEN, m_params.wy, 1.0};
    vec2dp const W_e2_tip{m_params.wx, m_params.wy + AXIS_LEN, 1.0};

    QPointF const W_s  = toScreen(W_orig);
    QPointF const e1_s = toScreen(W_e1_tip);
    QPointF const e2_s = toScreen(W_e2_tip);

    // e1 axis: solid red, thick
    qp->setPen(QPen(Qt::red, 2.5, Qt::SolidLine));
    qp->setBrush(Qt::NoBrush);
    qp->drawPath(arrowLine(W_s, e1_s));
    qp->drawPath(arrowHead(W_s, e1_s));

    // e2 axis: solid green, thick
    qp->setPen(QPen(Qt::darkGreen, 2.5, Qt::SolidLine));
    qp->drawPath(arrowLine(W_s, e2_s));
    qp->drawPath(arrowHead(W_s, e2_s));

    // Origin dot
    qp->setPen(QPen(Qt::black, 2));
    qp->setBrush(QBrush(Qt::black));
    qp->drawEllipse(W_s, 3.0, 3.0);

    // Label "W"
    QFont lbl_font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl_font.setPointSize(9);
    lbl_font.setBold(true);
    qp->setFont(lbl_font);
    qp->setPen(QPen(Qt::black, 1));
    qp->drawText(W_s + QPointF(-14.0, -6.0), "W");

    // (M0 indicator removed — M0 is displayed in the text overlay instead)
}

// cm trajectory (light gray polyline)
void active_frame_trafo::drawTrajectory(QPainter* qp)
{
    if (m_trajectory.size() < 2) return;

    qp->setPen(QPen(QColor(160, 160, 160, 140), 1, Qt::SolidLine));
    qp->setBrush(Qt::NoBrush);
    for (size_t i = 1; i < m_trajectory.size(); ++i) {
        qp->drawLine(toScreen(m_trajectory[i - 1].pos), toScreen(m_trajectory[i].pos));
    }
}

// Body frame B: axes + flagpole + rectangular flag with Bosch logo
void active_frame_trafo::drawBody(QPainter* qp, vec2dp const& O_w,
                                  vec2dp const& e1_w, vec2dp const& e2_w)
{
    QPointF const O_s = toScreen(O_w);

    // --- Body-frame axes ---

    vec2dp const B_e1_tip{O_w.x + e1_w.x * AXIS_LEN, O_w.y + e1_w.y * AXIS_LEN, 1.0};
    vec2dp const B_e2_tip{O_w.x + e2_w.x * AXIS_LEN, O_w.y + e2_w.y * AXIS_LEN, 1.0};
    QPointF const Be1_s = toScreen(B_e1_tip);
    QPointF const Be2_s = toScreen(B_e2_tip);

    // e1_b: dark red dashed shaft, solid head
    qp->setPen(QPen(QColor(180, 0, 0), 2.0, Qt::DashLine));
    qp->setBrush(Qt::NoBrush);
    qp->drawPath(arrowLine(O_s, Be1_s));
    qp->setPen(QPen(QColor(180, 0, 0), 2.0, Qt::SolidLine));
    qp->drawPath(arrowHead(O_s, Be1_s));

    // e2_b: dark green dashed shaft, solid head
    qp->setPen(QPen(QColor(0, 120, 0), 2.0, Qt::DashLine));
    qp->drawPath(arrowLine(O_s, Be2_s));
    qp->setPen(QPen(QColor(0, 120, 0), 2.0, Qt::SolidLine));
    qp->drawPath(arrowHead(O_s, Be2_s));

    // Label "B"
    QFont lbl_font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl_font.setPointSize(9);
    lbl_font.setBold(true);
    qp->setFont(lbl_font);
    qp->setPen(QPen(QColor(100, 0, 0), 1));
    qp->drawText(O_s + QPointF(6.0, -6.0), "B");

    // --- Flagpole: base at body-frame (POLE_OFFSET, 0), along e2_b ---

    vec2dp const pole_base_w{O_w.x + e1_w.x * POLE_OFFSET,
                              O_w.y + e1_w.y * POLE_OFFSET, 1.0};
    vec2dp const pole_tip_w{pole_base_w.x + e2_w.x * POLE_LEN,
                             pole_base_w.y + e2_w.y * POLE_LEN, 1.0};
    QPointF const pole_base_s = toScreen(pole_base_w);
    QPointF const pole_tip_s  = toScreen(pole_tip_w);

    qp->setPen(QPen(QColor(60, 60, 60), 1.8, Qt::SolidLine));
    qp->setBrush(Qt::NoBrush);
    qp->drawLine(pole_base_s, pole_tip_s);

    // --- Rectangular flag: hangs to +e1_b and -e2_b from pole tip ---
    // Corners (body-frame offsets from pole_tip_w):
    //   P0 = pole_tip                 (top-left, where pole meets flag)
    //   P1 = pole_tip + FLAG_W*e1_b  (top-right)
    //   P2 = P1 - FLAG_H*e2_b        (bottom-right)
    //   P3 = pole_tip - FLAG_H*e2_b  (bottom-left)

    auto flag_pt = [&](double u1, double u2) -> QPointF {
        vec2dp const p{pole_tip_w.x + u1 * e1_w.x + u2 * e2_w.x,
                       pole_tip_w.y + u1 * e1_w.y + u2 * e2_w.y, 1.0};
        return toScreen(p);
    };

    QPolygonF flag_poly;
    flag_poly << flag_pt(0.0,    0.0)
              << flag_pt(FLAG_W, 0.0)
              << flag_pt(FLAG_W, -FLAG_H)
              << flag_pt(0.0,    -FLAG_H);

    qp->setPen(QPen(QColor(50, 100, 200), 1.2));
    qp->setBrush(QBrush(QColor(70, 130, 180, 200)));
    qp->drawPolygon(flag_poly);

    // --- Bosch anchor logo: render PNG pixmap into the flag rectangle ---
    // The pixmap is embedded via Qt resources (ga_view.qrc) and cached statically.
    static QPixmap const bosch_logo = QPixmap(":/bosch_anchor.png");

    if (!bosch_logo.isNull()) {
        double const px_per_x = std::abs(cs->x.px_density_rng());
        double const px_per_y = std::abs(cs->y.px_density_rng());
        double const flag_w_px = FLAG_W * px_per_x;
        double const flag_h_px = FLAG_H * px_per_y;

        // Body rotation angle in screen space.
        // Model e1_w = (cos θ, sin θ); screen y is flipped → screen angle = −θ.
        // QPainter::rotate is clockwise-positive.
        double const angle_deg = std::atan2(-e1_w.y, e1_w.x) * 180.0 / M_PI;

        // Fit the PNG at its native aspect ratio inside the flag rectangle
        // so the logo circle appears undistorted.
        double const png_ar = static_cast<double>(bosch_logo.width()) / bosch_logo.height();
        double dest_w = flag_h_px * png_ar; // try fitting by height first
        double dest_h = flag_h_px;
        if (dest_w > flag_w_px) {           // too wide → fit by width instead
            dest_w = flag_w_px;
            dest_h = dest_w / png_ar;
        }
        double const x_off = 0.5 * (flag_w_px - dest_w);
        double const y_off = 0.5 * (flag_h_px - dest_h);

        qp->save();
        qp->translate(pole_tip_s);
        qp->rotate(angle_deg);
        qp->drawPixmap(QRectF(x_off, y_off, dest_w, dest_h),
                       bosch_logo, QRectF(bosch_logo.rect()));
        qp->restore();
    }

    // --- Body-origin marker (small filled circle) ---

    qp->setPen(QPen(Qt::black, 1));
    qp->setBrush(QBrush(QColor(60, 60, 60)));
    qp->drawEllipse(O_s, 3.5, 3.5);
}

// Home marker: shows where M0 places the body origin when B_b = 0.
// Drawn as short orange dotted axes + "M₀·O" label.
// For M0=id this coincides with W; for M0≠id it is visibly offset, making
// the effect of M0 on the initial body placement immediately apparent.
void active_frame_trafo::drawHomeMarker(QPainter* qp)
{
    // For world_frame_drive: home = initial config T(cm_ox,cm_oy) ⟇ R(0)
    // For body-frame drive:  home = M0 (the preset initial motor)
    mvec2dp_u const M_home = m_params.world_frame_drive
        ? exp(0.5 * vec2dp{-m_params.cm_oy, m_params.cm_ox, 0.0})
        : m_params.M0;

    // Body origin and directions at t=0: apply M_home
    vec2dp const O_b {0.0, 0.0, 1.0};
    vec2dp const e1_b{1.0, 0.0, 0.0};
    vec2dp const e2_b{0.0, 1.0, 0.0};

    vec2dp const H_o  = move2dp(O_b,  M_home);
    vec2dp const H_e1 = move2dp(e1_b, M_home); // unit direction
    vec2dp const H_e2 = move2dp(e2_b, M_home);

    double const HL = AXIS_LEN * 0.55; // shorter than body-frame arrows
    vec2dp const H_e1_tip{H_o.x + H_e1.x * HL, H_o.y + H_e1.y * HL, 1.0};
    vec2dp const H_e2_tip{H_o.x + H_e2.x * HL, H_o.y + H_e2.y * HL, 1.0};

    QPointF const H_s   = toScreen(H_o);
    QPointF const He1_s = toScreen(H_e1_tip);
    QPointF const He2_s = toScreen(H_e2_tip);

    QColor const col(210, 110, 0); // orange

    // e1 at home: dotted orange shaft + solid head
    qp->setPen(QPen(col, 1.8, Qt::DotLine));
    qp->setBrush(Qt::NoBrush);
    qp->drawPath(arrowLine(H_s, He1_s));
    qp->setPen(QPen(col, 1.8, Qt::SolidLine));
    qp->drawPath(arrowHead(H_s, He1_s));

    // e2 at home: dotted orange shaft + solid head
    qp->setPen(QPen(col, 1.8, Qt::DotLine));
    qp->drawPath(arrowLine(H_s, He2_s));
    qp->setPen(QPen(col, 1.8, Qt::SolidLine));
    qp->drawPath(arrowHead(H_s, He2_s));

    // Home origin dot
    qp->setPen(QPen(col, 1));
    qp->setBrush(QBrush(col));
    qp->drawEllipse(H_s, 3.5, 3.5);

    // Label "M₀·O"
    QFont lbl = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    lbl.setPointSize(8);
    lbl.setBold(true);
    qp->setFont(lbl);
    qp->setPen(QPen(col, 1));
    qp->drawText(H_s + QPointF(6.0, -6.0), "M\u2080\u00B7O");
}

// Velocity-generator overlay: draws Ω_b (blue) and Ω_w (magenta) at the body
// origin, showing:
//   • translation part: arrow in direction (Ω.y, -Ω.x) [world-frame coords]
//   • rotation part: ↺/↻ arc whose span direction matches the rotation sense
//
// For Ω_b the direction is first expressed in body-frame coords and then
// rotated into world frame via e1_w / e2_w so the arrow is always body-aligned.
// For Ω_w the direction is already in world frame.
//
// Key insights made visible:
//   – Pure translation : Ω_b ≡ Ω_w (arrows always coincide)
//   – Rotation / Combined, M0=id    : arrows coincide at t=0, may diverge
//   – Rotation / Combined, M0≠id   : arrows differ immediately (M0 couples frames)
void active_frame_trafo::drawGeneratorVectors(QPainter* qp,
                                              vec2dp const& O_w,
                                              vec2dp const& e1_w,
                                              vec2dp const& e2_w)
{
    if (!m_params.draw_gen_vectors) return;
    vec2dp const& Omega_b = m_params.Omega_b;
    vec2dp const  Omega_w = move2dp(Omega_b, m_M);

    QPointF const O_s      = toScreen(O_w);
    double const  px_per_x = std::abs(cs->x.px_density_rng());

    // Draw one generator: translation arrow + rotation-sense arc, both at body origin.
    // dir_x, dir_y: velocity direction in world-frame coords.
    auto draw_gen = [&](QColor const& col, double dir_x, double dir_y, double omega_z) {

        // Translation component: arrow from body origin
        double const tmag = std::sqrt(dir_x * dir_x + dir_y * dir_y);
        if (tmag > 1e-6) {
            double const len = std::min(tmag * VEC_SCALE, VEC_MAX);
            double const ux  = dir_x / tmag;
            double const uy  = dir_y / tmag;
            vec2dp const tip_w{O_w.x + ux * len, O_w.y + uy * len, 1.0};
            QPointF const tip_s = toScreen(tip_w);
            qp->setPen(QPen(col, 2.0, Qt::SolidLine));
            qp->setBrush(Qt::NoBrush);
            qp->drawPath(arrowLine(O_s, tip_s));
            qp->drawPath(arrowHead(O_s, tip_s));
        }

        // Rotation component: open arc (↺ / ↻) centred at body origin.
        // omega_z > 0  →  CCW in math coords  →  CW on Qt screen (y flipped)
        //              →  negative spanAngle in Qt's 1/16-degree units
        if (std::abs(omega_z) > 1e-4) {
            double const r_px = ARC_R * px_per_x;
            QRectF const rect(O_s.x() - r_px, O_s.y() - r_px, 2.0 * r_px, 2.0 * r_px);
            int const start_qt = 45 * 16;
            int const span_qt  = (omega_z > 0) ? -270 * 16 : +270 * 16;
            qp->setPen(QPen(col, 2.0, Qt::SolidLine));
            qp->setBrush(Qt::NoBrush);
            qp->drawArc(rect, start_qt, span_qt);
        }
    };

    // Ω_b: body-frame direction (Omega_b.y, -Omega_b.x) rotated to world frame
    double const Ob_dir_x = Omega_b.y * e1_w.x + (-Omega_b.x) * e2_w.x;
    double const Ob_dir_y = Omega_b.y * e1_w.y + (-Omega_b.x) * e2_w.y;

    // Ω_w: world-frame direction directly from generator components
    double const Ow_dir_x =  Omega_w.y;
    double const Ow_dir_y = -Omega_w.x;

    // Draw Ω_w (magenta) first so Ω_b (blue) appears on top when they overlap
    draw_gen(QColor(180, 0, 180), Ow_dir_x, Ow_dir_y, Omega_w.z);
    draw_gen(QColor(0,  60, 200), Ob_dir_x, Ob_dir_y, Omega_b.z);
}

// Text overlay: B_b, B_w (with decoded pivot when rotating), M0, M(t), time.
// Uses a running line counter so conditional pivot sub-lines don't shift
// the fixed-layout lines in the pure-translation scene (B.z == 0 → no sub-line).
void active_frame_trafo::drawTextOverlay(QPainter* qp)
{
    QFont mono = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    mono.setPointSize(9);
    qp->setFont(mono);

    // Anchor: fixed model-space position left of world frame marker.
    // text_wy sentinel −999.0 falls back to wy + 0.25.
    double const twy = (m_params.text_wy > -990.0) ? m_params.text_wy
                                                    : m_params.wy + 0.25;
    QPointF const anchor = toScreen(vec2dp{m_params.wx - 3.3, twy, 1.0});

    QFontMetrics const fm(mono);
    int const lh = fm.height() + 2;
    int line = 0;

    // Pivot sub-lines are shown for Rotation and Combined scenes (B.z grows
    // during animation).  Layout is fixed by kin_case, not by the runtime value
    // of B.z, so reset never causes lines to appear/disappear mid-animation.
    // Translation scenes never show pivot lines (B.z == 0 always).
    bool const show_pivot = (m_params.kin_case != kin_case_t::Translation);
    double const pivot_eps = 1e-4;

    auto at = [&](int l) { return anchor + QPointF(0.0, l * lh); };

    // Helper: one pivot decoded sub-line; shows "n/a" at t=0 when B.z≈0.
    auto draw_pivot_line = [&](vec2dp const& B, QString const& label, QColor const& col) {
        qp->setPen(QPen(col, 1));
        if (std::abs(B.z) > pivot_eps) {
            qp->drawText(at(line++),
                QString("     \u2192 %1=(%2, %3),  \u03C6=%4 rad")
                    .arg(label)
                    .arg(B.x / B.z, 6, 'f', 3)
                    .arg(B.y / B.z, 6, 'f', 3)
                    .arg(B.z,       6, 'f', 3));
        } else {
            qp->drawText(at(line++),
                QString("     \u2192 %1=(  n/a,   n/a),  \u03C6=%2 rad")
                    .arg(label)
                    .arg(B.z, 6, 'f', 3));
        }
    };

    // --- B_b (blue) ---
    qp->setPen(QPen(QColor(0, 60, 200), 1));
    qp->drawText(at(line++),
        QString("B_b=(%1, %2, %3)")
            .arg(m_B_b.x, 7, 'f', 3)
            .arg(m_B_b.y, 7, 'f', 3)
            .arg(m_B_b.z, 7, 'f', 3));
    if (show_pivot) draw_pivot_line(m_B_b, "pivot_b", QColor(0, 80, 160));

    // --- B_w (magenta) ---
    qp->setPen(QPen(QColor(180, 0, 180), 1));
    qp->drawText(at(line++),
        QString("B_w=(%1, %2, %3)")
            .arg(m_B_w.x, 7, 'f', 3)
            .arg(m_B_w.y, 7, 'f', 3)
            .arg(m_B_w.z, 7, 'f', 3));
    if (show_pivot) draw_pivot_line(m_B_w, "pivot_w", QColor(140, 0, 140));

    // --- M0 / initial-config block (dark olive, constant) ---
    if (m_params.world_frame_drive) {
        // Show world-frame drive parameters instead of M0 motor components
        qp->setPen(QPen(QColor(100, 80, 0), 1));
        qp->drawText(at(line++),
            QString("v_w =(%1, %2),  \u03C9=%3 rad/s")
                .arg(m_params.cm_vx,    6, 'f', 3)
                .arg(m_params.cm_vy,    6, 'f', 3)
                .arg(m_params.cm_omega, 6, 'f', 3));
        qp->setPen(QPen(QColor(130, 105, 0), 1));
        qp->drawText(at(line++),
            QString("     \u2192 home: O_w=(%1, %2),  \u03C6\u2080=%3 rad")
                .arg(m_params.cm_ox, 6, 'f', 3)
                .arg(m_params.cm_oy, 6, 'f', 3)
                .arg(0.0,            6, 'f', 3));
    } else {
        qp->setPen(QPen(QColor(100, 80, 0), 1));
        qp->drawText(at(line++),
            QString("M\u2080 =(e1=%1, e2=%2, e3=%3, e321=%4)")
                .arg(m_params.M0.c0, 6, 'f', 3)
                .arg(m_params.M0.c1, 6, 'f', 3)
                .arg(m_params.M0.c2, 6, 'f', 3)
                .arg(m_params.M0.c3, 6, 'f', 3));
        // M0 geometric decoding: body-origin position and frame orientation in world
        {
            vec2dp const O_h  = move2dp(vec2dp{0.0, 0.0, 1.0}, m_params.M0);
            vec2dp const e1_h = move2dp(vec2dp{1.0, 0.0, 0.0}, m_params.M0);
            double const phi  = std::atan2(e1_h.y, e1_h.x);
            qp->setPen(QPen(QColor(130, 105, 0), 1));
            qp->drawText(at(line++),
                QString("     \u2192 home: O_w=(%1, %2),  \u03C6\u2080=%3 rad")
                    .arg(O_h.x, 6, 'f', 3)
                    .arg(O_h.y, 6, 'f', 3)
                    .arg(phi,   6, 'f', 3));
        }
    }

    // --- M(t) (dark gray) ---
    qp->setPen(QPen(QColor(80, 80, 80), 1));
    qp->drawText(at(line++),
        QString("M  =(e1=%1, e2=%2, e3=%3, e321=%4)")
            .arg(m_M.c0, 6, 'f', 3)
            .arg(m_M.c1, 6, 'f', 3)
            .arg(m_M.c2, 6, 'f', 3)
            .arg(m_M.c3, 6, 'f', 3));

    // --- time (dark blue) ---
    qp->setPen(QPen(Qt::darkBlue, 1));
    qp->drawText(at(line),
        QString("t  = %1 s%2")
            .arg(m_t, 5, 'f', 2)
            .arg(m_paused ? "  [PAUSED]" : ""));
}

// ---------------------------------------------------------------------------
// Paint
// ---------------------------------------------------------------------------

void active_frame_trafo::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                                QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    // Current motor and body-frame quantities in world frame
    if (m_params.world_frame_drive) {
        double const ox    = m_params.cm_ox + m_params.cm_vx * m_t;
        double const oy    = m_params.cm_oy + m_params.cm_vy * m_t;
        double const theta = m_params.cm_omega * m_t;
        mvec2dp_u const R_t = exp(0.5 * vec2dp{0.0, 0.0, theta});
        mvec2dp_u const T_t = exp(0.5 * vec2dp{-oy, ox, 0.0});
        m_M = rgpr(T_t, R_t);
    } else {
        m_M = rgpr(m_params.M0, exp(0.5 * m_B_b));
    }

    vec2dp const O_b{0.0, 0.0, 1.0};
    vec2dp const e1_b{1.0, 0.0, 0.0};
    vec2dp const e2_b{0.0, 1.0, 0.0};

    vec2dp const O_w  = move2dp(O_b, m_M);
    vec2dp const e1_w = move2dp(e1_b, m_M); // direction (z=0): not translated
    vec2dp const e2_w = move2dp(e2_b, m_M);

    drawTrajectory(qp);
    drawHomeMarker(qp);
    drawWorldFrame(qp, O_w, e1_w, e2_w);
    drawGeneratorVectors(qp, O_w, e1_w, e2_w);
    drawBody(qp, O_w, e1_w, e2_w);
    drawTextOverlay(qp);

    qp->restore();
}

// ---------------------------------------------------------------------------
// Bounding rect / shape
// ---------------------------------------------------------------------------

QRectF active_frame_trafo::boundingRect() const
{
    // Conservative bounding rect covering both W marker and animation range
    double const margin = 1.5;
    double xmin = m_params.wx - 4.0 - margin;
    double xmax = m_params.wx + 4.0 + margin;
    double ymin = m_params.wy - 2.0 - margin;
    double ymax = m_params.wy + 2.0 + margin;

    // Include trajectory
    for (TrajPt const& tp : m_trajectory) {
        xmin = std::min(xmin, tp.pos.x);
        xmax = std::max(xmax, tp.pos.x);
        ymin = std::min(ymin, tp.pos.y);
        ymax = std::max(ymax, tp.pos.y);
    }

    return QRectF(QPointF(cs->x.au_to_w(xmin), cs->y.au_to_w(ymax)),
                  QPointF(cs->x.au_to_w(xmax), cs->y.au_to_w(ymin)));
}

QPainterPath active_frame_trafo::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
