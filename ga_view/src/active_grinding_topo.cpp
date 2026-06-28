// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_grinding_topo.hpp"
#include "view_fonts.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <QTimer>
#include <algorithm>
#include <cmath>

using namespace hd::ga;
using namespace hd::ga::pga;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_grinding_topo::active_grinding_topo(Coordsys* cs, w_Coordsys* wcs,
                                           grinding_topo_params const& params,
                                           QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}, m_patch{params.patch}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    build_system();

    connect(wcs, &w_Coordsys::viewResized, this, &active_grinding_topo::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_grinding_topo::tick);
    m_timer->start(UPDATE_MS);
}

// ---------------------------------------------------------------------------
// Build the GA frame tree and read the surface speeds off the twist field
// ---------------------------------------------------------------------------

// The Fig.-1 two-chain tree, built in REAL units (micrometres, rad/s) so the surface
// speeds -- and hence the wavelengths -- come out in micrometres and match the paper. The
// wheel and chuck spin (kin_state omega) feed kinematic_system3dp::point_velocity, which
// supplies the two relative surface speeds the wavelengths are built from. Axial offsets
// lie along e1 and drop out of the chuck-plane (e2,e3) quantities used here.
void active_grinding_topo::build_system()
{
    double const ns = m_params.ns_rpm * pi / 30.0; // rpm -> rad/s
    double const nw = m_params.nw_rpm * pi / 30.0;
    double const Rw = m_params.Rw_um;
    double const r = m_params.r_um;
    double const off = Rw / std::sqrt(2.0);
    double const l3 = 0.3 * Rw, x_a = l3; // along e1 (hidden in the e2-e3 plane)

    m_sys = kinematic_system3dp{};
    m_sys.add_frame(static_frame3dp("chuck_ctr_stat")); // inertial root
    m_sys.add_frame(static_frame3dp("chuck_ctr_rot"),
                    kin_state3dp{.omega = vec3dp{ns, 0.0, 0.0, 0.0}},
                    m_sys.index_of("chuck_ctr_stat"));
    m_sys.add_frame(static_frame3dp("wafer_top_avg_rot", vec3dp{0.0, 0.0, 0.0, 1.0}),
                    kin_state3dp{}, m_sys.index_of("chuck_ctr_rot"));
    m_sys.add_frame(static_frame3dp("spindle_cm_stat", vec3dp{x_a, -off, off, 1.0},
                                    vec3dp{0.0, -pi / 2.0, 0.0, 0.0}),
                    kin_state3dp{}, m_sys.index_of("chuck_ctr_stat"));
    // wheel spins so it CO-ROTATES with the wafer about the chuck axis (Tao's
    // down-grinding mode, Table 2): the spindle local z points along -e1, so a -nw spin
    // about it gives +nw about e1 -- the same sense as the wafer's +ns. This is what
    // reproduces the paper's 7.6 deg grinding-mark inclination (the counter-rotating
    // sense gives ~11.5 deg).
    m_sys.add_frame(static_frame3dp("tool_top_avg_rot", vec3dp{0.0, 0.0, l3, 1.0}),
                    kin_state3dp{.omega = vec3dp{0.0, 0.0, -nw, 0.0}},
                    m_sys.index_of("spindle_cm_stat"));
    m_sys.add_frame(static_frame3dp("tool_surface_avg_at_R", vec3dp{Rw, 0.0, 0.0, 1.0}),
                    kin_state3dp{}, m_sys.index_of("tool_top_avg_rot"));

    // wheel rim speed v_w = |point_velocity| of the rim grain (= n_w * R_w)
    vec3dp const grain = unitize(
        move3dp(O_3dp, m_sys.get_pos_trafo("tool_surface_avg_at_R", "chuck_ctr_stat")));
    m_v_w = to_val(bulk_nrm(m_sys.point_velocity(grain, "tool_surface_avg_at_R")));

    // wafer surface speed v_s = |point_velocity| of a wafer point at radius r (= n_s * r)
    vec3dp const P_r = vec3dp{0.0, r, 0.0, 1.0};
    m_v_s = to_val(bulk_nrm(m_sys.point_velocity(P_r, "wafer_top_avg_rot")));

    m_lam_m = m_v_w / m_params.f_b; // mark-direction wavelength [um]
    m_lam_c = m_v_s / m_params.f_b; // circumferential wavelength [um]

    // grinding-mark inclination: the grain velocity RELATIVE to the wafer, evaluated at
    // the TRUE contact point Q (a wafer point at radius r that also lies on the wheel
    // rim), then measured against the local radial direction. With the co-rotating wheel
    // this reproduces Tao's 7.6 deg mark inclination (Fig. 9) from the GA twist field.
    //
    // Q solves |Q| = r and |Q - C| = R_w (C = wheel centre in the chuck plane). With
    // |C| = R_w this gives Q.C = r^2/2, a line intersected with the radius-r circle.
    vec3dp const Croot = unitize(
        move3dp(O_3dp, m_sys.get_pos_trafo("tool_top_avg_rot", "chuck_ctr_stat")));
    double const a = Croot.y, b = Croot.z, rhs = 0.5 * r * r;
    double Qy = 0.0, Qz = r;
    if (std::abs(b) > 1e-9) {
        double const A = 1.0 + (a / b) * (a / b);
        double const B = -2.0 * a * rhs / (b * b);
        double const Cc = (rhs / b) * (rhs / b) - r * r;
        double const disc = B * B - 4.0 * A * Cc;
        if (disc >= 0.0) {
            Qy = (-B + std::sqrt(disc)) / (2.0 * A);
            Qz = (rhs - a * Qy) / b;
        }
    }
    vec3dp const Q = vec3dp{0.0, Qy, Qz, 1.0};
    vec3dp const rel = m_sys.point_velocity(Q, "tool_surface_avg_at_R") -
                       m_sys.point_velocity(Q, "wafer_top_avg_rot");
    double const rn = std::sqrt(Qy * Qy + Qz * Qz);
    double const rel_r = (rel.y * Qy + rel.z * Qz) / rn;  // radial component
    double const rel_c = (-rel.y * Qz + rel.z * Qy) / rn; // circumferential component
    m_beta = std::atan2(rel_c, rel_r);                    // mark angle from radial [rad]
    if (m_beta > 0.5 * pi) m_beta -= pi;                  // reduce to a line tilt
    if (m_beta <= -0.5 * pi) m_beta += pi;

    set_layout();
    rebuild_image();
}

// ---------------------------------------------------------------------------
// Patch geometry
// ---------------------------------------------------------------------------

// physical window of the rendered patch [um]: a near-square crop for the circumferential
// waviness (Fig. 9), or a wide thin strip for the mark-direction waviness (Fig. 10).
double active_grinding_topo::patch_w_um() const
{
    return (m_patch == gt_patch::wcd) ? 600.0 : 6000.0;
}
double active_grinding_topo::patch_h_um() const
{
    return (m_patch == gt_patch::wcd) ? 600.0 : 30.0;
}

// placement on the canvas [abstract units]: the patch sits to the RIGHT of the text
// legend (which occupies the top-left, out to ~x=-0.6), leaving room for a colour legend
// on its right. WCD is a centred square; WMD is a wide strip (its height is exaggerated
// -- the true patch is 6000 x 30 um, a 200:1 sliver).
void active_grinding_topo::set_layout()
{
    if (m_patch == gt_patch::wcd) {
        m_x0 = -0.1;
        m_x1 = 2.5;
        m_y0 = -1.3;
        m_y1 = 1.3; // 2.6 x 2.6 square
    }
    else {
        m_x0 = -0.1;
        m_x1 = 2.5;
        m_y0 = -1.0;
        m_y1 = 1.0; // near-square box (the 30 um height is exaggerated, like Tao Fig. 10)
    }
}

// ---------------------------------------------------------------------------
// Topography height field -> heatmap image
// ---------------------------------------------------------------------------

// Classic "jet" colormap (blue -> cyan -> green -> yellow -> red), t in [0,1].
QRgb active_grinding_topo::jet(double t)
{
    auto cl = [](double x) { return std::clamp(x, 0.0, 1.0); };
    double const r = cl(1.5 - std::abs(4.0 * t - 3.0));
    double const g = cl(1.5 - std::abs(4.0 * t - 2.0));
    double const b = cl(1.5 - std::abs(4.0 * t - 1.0));
    return qRgb(int(255.0 * r), int(255.0 * g), int(255.0 * b));
}

// Carve the topography height field z(s_m, s_c) = A_b*[sin(2pi*mark/lambda_m)
// + sin(2pi*s_c/lambda_c)] over the patch: horizontal = mark/radial direction (carries
// lambda_m), vertical = circumferential (carries lambda_c). The mark family is tilted by
// the GA-computed inclination m_beta; m_phase scrolls the pattern as the feed advances.
void active_grinding_topo::rebuild_image()
{
    int const Wpx = (m_patch == gt_patch::wcd) ? 480 : 960;
    int const Hpx = (m_patch == gt_patch::wcd) ? 480 : 300;
    double const W = patch_w_um(), H = patch_h_um();

    m_image = QImage(Wpx, Hpx, QImage::Format_RGB32);
    for (int j = 0; j < Hpx; ++j) {
        double const sc = (Hpx > 1) ? H * double(j) / double(Hpx - 1) : 0.0;
        for (int i = 0; i < Wpx; ++i) {
            double const sm = (Wpx > 1) ? W * double(i) / double(Wpx - 1) : 0.0;
            double const arg_m =
                2.0 * pi * (sm * std::cos(m_beta) + sc * std::sin(m_beta)) / m_lam_m;
            double const arg_c = 2.0 * pi * sc / m_lam_c + m_phase;
            double const z = std::sin(arg_m + m_phase) + std::sin(arg_c); // in [-2, 2]
            m_image.setPixel(i, j, jet((z + 2.0) / 4.0));
        }
    }
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_grinding_topo::tick()
{
    if (m_paused) return;
    m_phase += 0.04; // slow feed-scroll
    rebuild_image();
    update();
}

void active_grinding_topo::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_grinding_topo::cyclePatch()
{
    m_patch = (m_patch == gt_patch::wcd) ? gt_patch::wmd : gt_patch::wcd;
    set_layout();
    rebuild_image();
    prepareGeometryChange();
    update();
}

void active_grinding_topo::resetAnimation()
{
    m_phase = 0.0;
    m_patch = m_params.patch;
    build_system();
    prepareGeometryChange();
    update();
}

void active_grinding_topo::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Painting
// ---------------------------------------------------------------------------

QPointF active_grinding_topo::auPt(double ax, double ay) const
{
    return QPointF(cs->x.au_to_w(ax), cs->y.au_to_w(ay));
}

// Vertical colour legend: the jet ramp from the lowest surface height (blue, bottom) to
// the highest (red, top), with tick labels in micrometres. The field z = A_b*(sin + sin)
// spans
// +-2*A_b, so the legend runs from -2*A_b to +2*A_b.
void active_grinding_topo::drawColorbar(QPainter* qp, double cbx0, double cbx1,
                                        double cby0, double cby1) const
{
    int const nseg = 64;
    for (int k = 0; k < nseg; ++k) {
        double const t0 = double(k) / nseg, t1 = double(k + 1) / nseg;
        double const y0 = cby0 + t0 * (cby1 - cby0), y1 = cby0 + t1 * (cby1 - cby0);
        QRectF seg(auPt(cbx0, y1), auPt(cbx1, y0)); // higher t -> higher (top)
        qp->fillRect(seg, QColor(jet(0.5 * (t0 + t1))));
    }
    qp->setPen(QPen(QColor(40, 40, 40), 1.2));
    qp->setBrush(Qt::NoBrush);
    qp->drawRect(QRectF(auPt(cbx0, cby1), auPt(cbx1, cby0)));

    // tick labels (surface height [um]) at t = 0, 0.25, .. 1
    QFont f = mono_font(9);
    qp->setFont(f);
    qp->setPen(QPen(QColor(20, 20, 20), 1));
    double const hmax = 2.0 * m_params.A_b; // +-2*A_b [um]
    for (int k = 0; k <= 4; ++k) {
        double const t = 0.25 * k;
        double const yy = cby0 + t * (cby1 - cby0);
        double const h = (2.0 * t - 1.0) * hmax; // map t in [0,1] -> [-hmax, +hmax]
        qp->drawText(auPt(cbx1 + 0.05, yy) + QPointF(0.0, 4.0),
                     QString("%1").arg(h, 5, 'f', 2));
    }
    f.setBold(true);
    qp->setFont(f);
    qp->drawText(auPt(cbx0 - 0.05, cby1 + 0.12), "height [um]");
}

void active_grinding_topo::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                                 QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    // blit the height-field image into the patch box (to the right of the text legend)
    QRectF const target(auPt(m_x0, m_y1), auPt(m_x1, m_y0)); // top-left, bottom-right
    qp->drawImage(target, m_image);
    qp->setPen(QPen(QColor(40, 40, 40), 1.2));
    qp->setBrush(Qt::NoBrush);
    qp->drawRect(target);

    // colour legend, just right of the patch (taller than a thin WMD strip so it reads)
    double const yc = 0.5 * (m_y0 + m_y1);
    double const cbh = std::max(0.5 * (m_y1 - m_y0), 0.9);
    drawColorbar(qp, m_x1 + 0.18, m_x1 + 0.34, yc - cbh, yc + cbh);

    // patch-direction axis labels
    QFont lbl = mono_font(10, /*bold=*/true);
    qp->setFont(lbl);
    qp->setPen(QPen(QColor(30, 30, 30), 1));
    qp->drawText(auPt(0.5 * (m_x0 + m_x1) - 0.5, m_y0 - 0.30), "mark direction (radial)");
    qp->save();
    QPointF const cl = auPt(m_x0 - 0.20, yc);
    qp->translate(cl.x(), cl.y());
    qp->rotate(-90.0);
    qp->drawText(QPointF(-58.0, 0.0), "circumferential");
    qp->restore();

    // numeric extent of the patch [um]: mark direction 0..W along the bottom edge,
    // circumferential 0..H along the left edge (s_c = 0 is the top image row).
    QFont tick = mono_font(9);
    qp->setFont(tick);
    qp->setPen(QPen(QColor(30, 30, 30), 1));
    double const W = patch_w_um(), H = patch_h_um();
    qp->drawText(auPt(m_x0, m_y0) + QPointF(2.0, 15.0), "0");
    qp->drawText(auPt(m_x1, m_y0) + QPointF(-52.0, 15.0),
                 QString("%1 um").arg(W, 0, 'f', 0));
    qp->drawText(auPt(m_x0, m_y1) + QPointF(-26.0, 0.0), "0");
    qp->drawText(auPt(m_x0, m_y0) + QPointF(-58.0, 4.0),
                 QString("%1 um").arg(H, 0, 'f', 0));

    // info block, anchored BELOW the colour legend so it never runs into it (the legend
    // can extend lower than a thin WMD patch): the GA-derived wavelengths, the preset,
    // and the animated quantity.
    QFont info = mono_font(10, /*bold=*/true);
    qp->setFont(info);
    qp->setPen(QPen(QColor(20, 20, 20), 1));
    // first line sits below BOTH the colour legend and the "mark direction" label
    double const iy = std::min(m_y0 - 0.55, (yc - cbh) - 0.22);
    char const* name = (m_patch == gt_patch::wcd) ? "WCD (Fig. 9)" : "WMD (Fig. 10)";
    qp->drawText(auPt(m_x0, iy), QString("%1   patch %2 x %3 um   (C: toggle WCD/WMD)")
                                     .arg(name)
                                     .arg(patch_w_um(), 0, 'f', 0)
                                     .arg(patch_h_um(), 0, 'f', 0));
    qp->drawText(
        auPt(m_x0, iy - 0.24),
        QString("f_b = %1 Hz   lambda_m = %2 um   lambda_c = %3 um  (point_vel.)")
            .arg(m_params.f_b, 0, 'f', 0)
            .arg(m_lam_m, 0, 'f', 0)
            .arg(m_lam_c, 0, 'f', 0));
    // the animated quantity: the feed advancing the surface under the wheel scrolls the
    // pattern; show the circumferential scroll offset (wraps every lambda_c).
    double const scroll = std::fmod(m_phase / (2.0 * pi) * m_lam_c, m_lam_c);
    qp->drawText(auPt(m_x0, iy - 0.48),
                 QString("ANIMATED feed scroll = %1 / %2 um (lambda_c; SPACE to pause)")
                     .arg(scroll, 0, 'f', 0)
                     .arg(m_lam_c, 0, 'f', 0));

    qp->restore();
}

QRectF active_grinding_topo::boundingRect() const
{
    // generous: covers the patch, the colour legend to its right, and the info block,
    // which is anchored below the (possibly taller-than-patch) colour legend.
    double const yc = 0.5 * (m_y0 + m_y1);
    double const cbh = std::max(0.5 * (m_y1 - m_y0), 0.9);
    double const iy = std::min(m_y0 - 0.55, (yc - cbh) - 0.22);
    double const ybot = iy - 0.7; // below the lowest info line (iy - 0.48)
    return QRectF(auPt(m_x0 - 0.7, m_y1 + 0.5), auPt(m_x1 + 1.2, ybot));
}

QPainterPath active_grinding_topo::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
