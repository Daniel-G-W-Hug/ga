// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_merry_go_round.hpp"
#include "active_common.hpp"

#include "ga/ga_pga.hpp"

#include <QFontDatabase>
#include <QTimer>
#include <algorithm>
#include <cmath>
#include <string>

using namespace hd::ga;
using namespace hd::ga::pga;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

active_merry_go_round::active_merry_go_round(Coordsys* cs, w_Coordsys* wcs,
                                             merry_go_round_params const& params,
                                             QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_params{params}
{
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    build_system();

    connect(wcs, &w_Coordsys::viewResized, this, &active_merry_go_round::viewChanged);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &active_merry_go_round::tick);
    m_timer->start(UPDATE_MS);
}

std::string active_merry_go_round::tt_name(int i) const
{
    return "T" + std::to_string(i);
}

// (Re)build the kinematic system at t=0: ground -> platform -> N turntables.
void active_merry_go_round::build_system()
{
    m_ds = kinematic_system2dp{};
    m_ds.add_frame(static_frame2dp("W"));                      // ground (inertial root)
    m_ds.add_frame(static_frame2dp("P", m_params.centre, 0.0), // rotating platform
                   kin_state2dp{.omega = m_params.platform_omega});
    int const n = std::max(1, m_params.n_turntables);
    for (int i = 0; i < n; ++i) {
        double const a = 2.0 * M_PI * i / n;
        vec2dp const mount{m_params.radius * std::cos(a), m_params.radius * std::sin(a),
                           1.0};
        m_ds.add_frame(static_frame2dp(tt_name(i), mount, 0.0),
                       kin_state2dp{.omega = m_params.turntable_omega},
                       m_ds.index_of("P"));
    }

    m_path.clear();
    m_plot.clear();
    m_t = 0.0;
    record_sample();
}

// Append the current marked-point world position (path) and its magnitude sample (plots),
// then trim both histories. Called once per tick and at (re)build.
void active_merry_go_round::record_sample()
{
    vec2dp const Pm = move2dp(m_params.marked_point, m_ds.get_pos_trafo("T0", "W"));
    m_path.push_back(Pm);
    if (m_path.size() > MAX_PATH) m_path.pop_front();

    vec2dp const vv = m_ds.point_velocity(Pm, "T0");
    vec2dp const aa = m_ds.point_acceleration(Pm, "T0");
    double const theta = m_t * m_params.platform_omega; // platform angle (constant omega)
    double const r =
        std::hypot(Pm.x - m_params.centre.x, Pm.y - m_params.centre.y); // orbit radius
    m_plot.push_back({theta, r, std::hypot(vv.x, vv.y), std::hypot(aa.x, aa.y)});
    while (m_plot.size() > 1 && m_plot.front().theta < theta - PLOT_WINDOW)
        m_plot.pop_front();
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void active_merry_go_round::tick()
{
    if (m_paused) return;

    m_ds.step(m_params.dt);
    m_t += m_params.dt;
    record_sample();

    prepareGeometryChange();
    update();
}

void active_merry_go_round::togglePause()
{
    m_paused = !m_paused;
    update();
}

void active_merry_go_round::toggleCircles()
{
    m_show_circles = !m_show_circles;
    update();
}

void active_merry_go_round::resetAnimation()
{
    build_system();
    prepareGeometryChange();
    update();
}

void active_merry_go_round::viewChanged() { update(); }

// ---------------------------------------------------------------------------
// Coordinate helper
// ---------------------------------------------------------------------------

QPointF active_merry_go_round::toScreen(vec2dp const& pt) const
{
    // world coords -> canvas coords is a constant rigid translation by world_origin
    // (where the world frame W is anchored on the diagram grid).
    return QPointF(cs->x.au_to_w(pt.x + m_params.world_origin.x),
                   cs->y.au_to_w(pt.y + m_params.world_origin.y));
}

// ---------------------------------------------------------------------------
// Painting
// ---------------------------------------------------------------------------

void active_merry_go_round::paint(QPainter* qp, QStyleOptionGraphicsItem const* option,
                                  QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));
    qp->save();

    if (m_params.layout == mgr_layout::dashboard) {
        paintDashboard(qp);
    }
    else {
        // full merry-go-round in world coordinates (toScreen mapping, includes labels)
        double const px = std::abs(cs->x.px_density_rng());
        double const py = std::abs(cs->y.px_density_rng());
        drawSystem(qp, [this](vec2dp const& p) { return toScreen(p); }, px, py, true);
    }

    qp->restore();
}

// Draw the merry-go-round geometry through the supplied world->screen mapping `map`
// (px/py = pixels per world unit for that mapping). Used both for the full-size scene
// (map = toScreen) and the dashboard thumbnail (map = a fit-to-box scaling). `full`
// enables the W-frame glyph and the P/Ti/v/a text labels (omitted in the thumbnail).
void active_merry_go_round::drawSystem(QPainter* qp,
                                       std::function<QPointF(vec2dp const&)> const& map,
                                       double px, double py, bool full)
{
    auto draw_disk = [&](vec2dp const& c_world, double r_world, QColor const& fill,
                         QColor const& edge) {
        qp->setPen(QPen(edge, 1.5));
        qp->setBrush(QBrush(fill));
        qp->drawEllipse(map(c_world), r_world * px, r_world * py);
    };
    auto draw_spoke = [&](vec2dp const& a_world, vec2dp const& b_world,
                          QColor const& col) {
        qp->setPen(QPen(col, 1.8, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        qp->drawLine(map(a_world), map(b_world));
    };

    // Free-vector arrow anchored at base_w, pointing along vec_w (a velocity/acceleration
    // field value), with vec_w scaled to world units for legibility. A thin circle of
    // radius = scaled vector length is drawn around base_w so the momentary magnitude
    // (and its fluctuation over the orbit) is readable independently of the arrow's
    // direction.
    auto draw_arrow = [&](vec2dp const& base_w, vec2dp const& vec_w, double scale,
                          QColor const& col) {
        double const mag = std::hypot(vec_w.x, vec_w.y) * scale; // arrow length [world u]
        QPointF const b_s = map(base_w);
        QPointF const t_s =
            map(vec2dp{base_w.x + vec_w.x * scale, base_w.y + vec_w.y * scale, 1.0});
        qp->setBrush(Qt::NoBrush);
        if (m_show_circles) {
            qp->setPen(
                QPen(col, 0.6, Qt::SolidLine)); // momentary-magnitude circle (thin)
            qp->drawEllipse(b_s, mag * px, mag * py);
        }
        qp->setPen(QPen(col, 2.0, Qt::SolidLine)); // arrow
        qp->drawPath(arrowLine(b_s, t_s));
        qp->drawPath(arrowHead(b_s, t_s));
    };

    // Small coordinate system of a frame: e1 (red) / e2 (green) unit axes obtained by
    // moving the body directions through the frame's motor M. Solid for the fixed ground
    // frame, dashed for moving frames (matches the world/body convention in
    // active_frame_trafo).
    auto draw_axes = [&](auto const& M, double len, Qt::PenStyle style) {
        vec2dp const o = move2dp(O_2dp, M);
        vec2dp const d1 = move2dp(vec2dp{1.0, 0.0, 0.0}, M); // e1 unit direction (world)
        vec2dp const d2 = move2dp(vec2dp{0.0, 1.0, 0.0}, M); // e2 unit direction (world)
        QPointF const o_s = map(o);
        QPointF const e1_s = map(vec2dp{o.x + d1.x * len, o.y + d1.y * len, 1.0});
        QPointF const e2_s = map(vec2dp{o.x + d2.x * len, o.y + d2.y * len, 1.0});
        qp->setBrush(Qt::NoBrush);
        qp->setPen(QPen(QColor(200, 0, 0), 1.6, style)); // e1: red
        qp->drawPath(arrowLine(o_s, e1_s));
        qp->setPen(QPen(QColor(200, 0, 0), 1.6, Qt::SolidLine));
        qp->drawPath(arrowHead(o_s, e1_s));
        qp->setPen(QPen(QColor(0, 140, 0), 1.6, style)); // e2: green
        qp->drawPath(arrowLine(o_s, e2_s));
        qp->setPen(QPen(QColor(0, 140, 0), 1.6, Qt::SolidLine));
        qp->drawPath(arrowHead(o_s, e2_s));
    };

    QFont lbl_font = mono_font(9, /*bold=*/true);

    // --- traced path of the marked point (drawn first, underneath) ---
    if (m_path.size() >= 2) {
        qp->setPen(QPen(QColor(40, 90, 200, 200), 1.5, Qt::SolidLine));
        qp->setBrush(Qt::NoBrush);
        for (size_t i = 1; i < m_path.size(); ++i)
            qp->drawLine(map(m_path[i - 1]), map(m_path[i]));
    }

    // --- platform: disk + structural arms to each turntable + axes + centre ---
    auto const Mp = m_ds.get_pos_trafo("P", "W");
    vec2dp const P_o = move2dp(O_2dp, Mp);
    int const n = std::max(1, m_params.n_turntables);
    draw_disk(P_o, m_params.radius + TT_RADIUS, QColor(210, 210, 210, 90),
              QColor(120, 120, 120));
    for (int i = 0; i < n; ++i) {
        double const a = 2.0 * M_PI * i / n;
        vec2dp const mount_w = move2dp(
            vec2dp{m_params.radius * std::cos(a), m_params.radius * std::sin(a), 1.0},
            Mp);
        draw_spoke(P_o, mount_w, QColor(120, 120, 120));
    }
    draw_axes(Mp, AXIS_LEN_P, Qt::DashLine); // platform frame (moving)
    qp->setPen(QPen(Qt::black, 2));
    qp->setBrush(QBrush(Qt::black));
    qp->drawEllipse(map(P_o), 3.0, 3.0);
    if (full) {
        qp->setFont(lbl_font);
        qp->setPen(QPen(Qt::black, 1));
        qp->drawText(map(P_o) + QPointF(8.0, -8.0), "P");
    }

    // --- turntables: disk + frame axes + origin dot + label ---
    for (int i = 0; i < n; ++i) {
        auto const Mt = m_ds.get_pos_trafo(tt_name(i), "W");
        vec2dp const T_o = move2dp(O_2dp, Mt);
        draw_disk(T_o, TT_RADIUS, QColor(180, 205, 235, 150), QColor(50, 100, 170));
        draw_axes(Mt, AXIS_LEN_T, Qt::SolidLine); // turntable frame (moving)
        qp->setPen(QPen(QColor(30, 60, 110), 2));
        qp->setBrush(QBrush(QColor(30, 60, 110)));
        qp->drawEllipse(map(T_o), 2.5, 2.5);
        if (full) {
            qp->setFont(lbl_font);
            qp->setPen(QPen(QColor(30, 60, 110), 1));
            qp->drawText(map(T_o) + QPointF(6.0, -6.0),
                         QString::fromStdString(tt_name(i)));
        }
    }

    // --- ground frame W (fixed at world origin, solid axes) ---
    if (full) {
        draw_axes(m_ds.get_pos_trafo("W", "W"), AXIS_LEN_W, Qt::SolidLine);
        qp->setFont(lbl_font);
        qp->setPen(QPen(Qt::black, 1));
        qp->drawText(map(O_2dp) + QPointF(-16.0, 14.0), "W");
    }

    // --- marked point on turntable 0 + its velocity / acceleration arrows ---
    vec2dp const P_marked = move2dp(m_params.marked_point, m_ds.get_pos_trafo("T0", "W"));

    // velocity (magenta) and acceleration (black) of the point fixed in T0, world frame
    vec2dp const v_marked = m_ds.point_velocity(P_marked, "T0");
    vec2dp const a_marked = m_ds.point_acceleration(P_marked, "T0");
    draw_arrow(P_marked, a_marked, ACC_SCALE, QColor(0, 0, 0));     // accel (under)
    draw_arrow(P_marked, v_marked, VEL_SCALE, QColor(190, 0, 190)); // velocity (over)

    qp->setPen(QPen(QColor(200, 60, 0), 2));
    qp->setBrush(QBrush(QColor(230, 90, 20)));
    qp->drawEllipse(map(P_marked), 4.0, 4.0);

    if (full) {
        // arrow-tip labels
        qp->setFont(lbl_font);
        qp->setPen(QPen(QColor(150, 0, 150), 1));
        qp->drawText(map(vec2dp{P_marked.x + v_marked.x * VEL_SCALE,
                                P_marked.y + v_marked.y * VEL_SCALE, 1.0}) +
                         QPointF(4.0, -4.0),
                     "v");
        qp->setPen(QPen(Qt::black, 1));
        qp->drawText(map(vec2dp{P_marked.x + a_marked.x * ACC_SCALE,
                                P_marked.y + a_marked.y * ACC_SCALE, 1.0}) +
                         QPointF(4.0, -4.0),
                     "a");
    }
}

// Dashboard layout: a mini merry-go-round in the top-right corner plus three stacked
// panels plotting |position| / |velocity| / |acceleration| of the marked point against
// the platform turning angle theta (windowed to one platform revolution).
void active_merry_go_round::paintDashboard(QPainter* qp)
{
    double const L = cs->x.nmin();
    double const R = cs->x.nmax();
    double const T = cs->y.nmax(); // top    (smaller pixel y)
    double const B = cs->y.nmin(); // bottom (larger pixel y)
    double const W = R - L;
    double const H = B - T;
    double const mgn = 12.0;

    QFont title_font = mono_font(9, /*bold=*/true);
    QFont small_font = mono_font(8);

    // --- top band: legend (left, overlay) + theta caption + mini (right) ---
    // capH reserves a row at the very top of the band for the shared theta caption, so it
    // stays inside the drawing area (it was previously placed above T and clipped).
    double const topH = std::min(H * 0.34, W * 0.34);
    double const capH = 18.0;
    double const miniSide = topH - capH - mgn;
    QRectF const miniBox(R - miniSide - mgn, T + capH, miniSide, miniSide);

    qp->setBrush(QColor(252, 252, 252, 235));
    qp->setPen(QPen(QColor(120, 120, 120), 1));
    qp->drawRect(miniBox);

    // mini mapping: fit (centre +/- reach) into miniBox, world-y up -> screen-y down
    double const reach = m_params.radius + TT_RADIUS + 0.3;
    double const s = (miniSide * 0.92) / (2.0 * reach);
    QPointF const miniCtr = miniBox.center();
    vec2dp const C = m_params.centre;
    auto miniMap = [=](vec2dp const& p) {
        return QPointF(miniCtr.x() + (p.x - C.x) * s, miniCtr.y() - (p.y - C.y) * s);
    };
    qp->save();
    qp->setClipRect(miniBox);
    drawSystem(qp, miniMap, s, s, false);
    qp->restore();

    // Shared current state for ALL three plots: the platform turning angle theta (the
    // common "now" on the x-axis). Drawn once in the reserved caption row at the top of
    // the band, right-aligned over the mini (inset from the right edge).
    double const theta_now = m_plot.empty() ? 0.0 : m_plot.back().theta;
    qp->setFont(title_font);
    qp->setPen(QColor(50, 50, 50));
    qp->drawText(
        QRectF(miniBox.left() - 90.0, T + 2.0, miniBox.width() + 74.0, capH - 2.0),
        Qt::AlignRight | Qt::AlignVCenter,
        QString("now:  theta = %1 rad  (%2 rev)")
            .arg(theta_now, 0, 'f', 2)
            .arg(theta_now / (2.0 * pi), 0, 'f', 2));

    // Heading / key-hints / colour key live in the diagram_legend (top-left, draggable),
    // set on the dashboard scene in w_mainwindow.cpp — not drawn here.
    //
    // BACKLOG (deferred, see kinematics roadmap): the plot panels below currently roll
    // their own X/Y pixel mapping. A nicer version would draw a real, independent
    // Coordsys per panel (theta on x, magnitude on y) so the axes/ticks/grid are rendered
    // by the shared cs machinery rather than hand-laid-out here. Worth it only if more
    // plots get added; skipped for now to keep the dashboard self-contained.

    // --- three stacked plot panels below the top band ---
    double const panelsTop = T + topH;
    double const panelsBot = B - mgn;
    double const gap = 10.0;
    double const panelW = W - 2.0 * mgn;
    double const panelH = (panelsBot - panelsTop - 2.0 * gap) / 3.0;
    double const panelX = L + mgn;

    auto draw_plot = [&](QRectF const& box, QString const& title, QString const& sym,
                         QColor const& col,
                         std::function<double(plot_sample const&)> const& acc,
                         double cur) {
        qp->setBrush(QColor(252, 252, 252, 235));
        qp->setPen(QPen(QColor(120, 120, 120), 1));
        qp->drawRect(box);

        double const pL = 50.0, pR = 12.0, pT = 18.0, pB = 16.0;
        QRectF const plot(box.left() + pL, box.top() + pT, box.width() - pL - pR,
                          box.height() - pT - pB);

        // axes
        qp->setPen(QPen(QColor(90, 90, 90), 1));
        qp->drawLine(plot.bottomLeft(), plot.bottomRight()); // theta axis
        qp->drawLine(plot.topLeft(), plot.bottomLeft());     // magnitude axis

        // title
        qp->setFont(title_font);
        qp->setPen(Qt::black);
        qp->drawText(QPointF(box.left() + 6.0, box.top() + 13.0), title);

        if (m_plot.size() >= 2) {
            double const th_max = m_plot.back().theta;
            double const th_min = th_max - PLOT_WINDOW;
            double vmax = 1.0e-9;
            for (plot_sample const& sm : m_plot)
                vmax = std::max(vmax, acc(sm));
            vmax *= 1.15;

            auto X = [&](double th) {
                return plot.left() + (th - th_min) / PLOT_WINDOW * plot.width();
            };
            auto Y = [&](double v) { return plot.bottom() - (v / vmax) * plot.height(); };

            // horizontal grid + y-scale labels at 0, 1/4, 1/2, 3/4, max
            qp->setFont(small_font);
            for (int k = 0; k <= 4; ++k) {
                double const yv = vmax * k / 4.0;
                double const yy = Y(yv);
                if (k > 0 && k < 4) {
                    qp->setPen(QPen(QColor(228, 228, 228), 1)); // light gridline
                    qp->drawLine(QPointF(plot.left(), yy), QPointF(plot.right(), yy));
                }
                qp->setPen(QColor(70, 70, 70));
                qp->drawText(QPointF(box.left() + 6.0, yy + 3.0),
                             QString::number(yv, 'f', 2));
            }
            qp->setPen(col);
            qp->drawText(QPointF(plot.right() - 150.0, box.top() + 13.0),
                         QString("%1 = %2  (at theta now)").arg(sym).arg(cur, 0, 'f', 3));

            // curve
            qp->setPen(QPen(col, 1.6));
            QPointF prev;
            bool have = false;
            for (plot_sample const& sm : m_plot) {
                if (sm.theta < th_min) continue;
                QPointF const p(X(sm.theta), Y(acc(sm)));
                if (have) qp->drawLine(prev, p);
                prev = p;
                have = true;
            }
        }

        qp->setFont(small_font);
        qp->setPen(QColor(70, 70, 70));
        qp->drawText(QPointF(plot.right() - 14.0, box.bottom() - 4.0), "theta");
    };

    double const r_now = m_plot.empty() ? 0.0 : m_plot.back().r;
    double const v_now = m_plot.empty() ? 0.0 : m_plot.back().v;
    double const a_now = m_plot.empty() ? 0.0 : m_plot.back().a;

    draw_plot(
        QRectF(panelX, panelsTop, panelW, panelH),
        "|r|  orbit radius about platform centre", "|r|", QColor(40, 90, 200),
        [](plot_sample const& sm) { return sm.r; }, r_now);
    draw_plot(
        QRectF(panelX, panelsTop + (panelH + gap), panelW, panelH), "|v|  speed", "|v|",
        QColor(190, 0, 190), [](plot_sample const& sm) { return sm.v; }, v_now);
    draw_plot(
        QRectF(panelX, panelsTop + 2.0 * (panelH + gap), panelW, panelH),
        "|a|  acceleration", "|a|", QColor(0, 0, 0),
        [](plot_sample const& sm) { return sm.a; }, a_now);
}

QRectF active_merry_go_round::boundingRect() const
{
    // dashboard fills the whole drawing area (it lays out in pixel coordinates)
    if (m_params.layout == mgr_layout::dashboard) {
        return QRectF(QPointF(cs->x.nmin(), cs->y.nmax()),
                      QPointF(cs->x.nmax(), cs->y.nmin()));
    }

    double const reach = m_params.radius + TT_RADIUS + 1.0;
    double xmin = m_params.centre.x - reach;
    double xmax = m_params.centre.x + reach;
    double ymin = m_params.centre.y - reach;
    double ymax = m_params.centre.y + reach;

    for (vec2dp const& p : m_path) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }

    // include the world-frame W glyph drawn at world origin (0,0), axes length AXIS_LEN_W
    xmin = std::min(xmin, 0.0);
    ymin = std::min(ymin, 0.0);
    xmax = std::max(xmax, AXIS_LEN_W);
    ymax = std::max(ymax, AXIS_LEN_W);

    // map world-coordinate corners to canvas via toScreen (applies the world_origin
    // shift)
    return QRectF(toScreen(vec2dp{xmin, ymax, 1.0}), toScreen(vec2dp{xmax, ymin, 1.0}));
}

QPainterPath active_merry_go_round::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}
