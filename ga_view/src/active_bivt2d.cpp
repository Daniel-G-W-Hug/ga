// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_bivt2d.hpp"
#include "active_common.hpp"

#include "ga/ga_ega.hpp"

active_bivt2d::active_bivt2d(Coordsys* cs, w_Coordsys* wcs, active_pt2d* uend,
                             active_pt2d* vend, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_uend{uend}, m_vend{vend}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    // setZValue(18);

    connect(wcs, &w_Coordsys::viewResized, m_uend, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_vend, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::rotateModeChanged, this, &active_bivt2d::rotateModeChanged);
}

void active_bivt2d::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                          QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    // clipping area is active area of coordsys
    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));

    // draw in item coordinate system
    qp->save();

    // beg_pos is always at origin
    QPointF beg_pos = QPointF(cs->x.au_to_w(0.0), cs->y.au_to_w(0.0));
    QPointF end_upos =
        QPointF(cs->x.au_to_w(m_uend->scenePos().x), cs->y.au_to_w(m_uend->scenePos().y));
    QPointF end_vpos =
        QPointF(cs->x.au_to_w(m_vend->scenePos().x), cs->y.au_to_w(m_vend->scenePos().y));

    QPointF tip_pos = end_upos + end_vpos - beg_pos;

    // draw parallelogram corresponding to u^v
    QPolygonF polygon;
    polygon << beg_pos << end_upos << tip_pos << end_vpos << beg_pos;

    // The sign has to be reversed here, since device coordinates are in a left-handed
    // system. The angle calculation itself is done in a classical right-handed system.
    //
    // Be aware: Depending on aspect ratio of x- vs. y-axis calculated angles will change,
    // if calculated from device coordinate input. However, they are visually consistent.
    //
    // Actual angle calculations for mathematical/physical purposes must use logical
    // coordinates as input for the calculation exclusively!
    qreal angle_rel = -angle_between_lines(beg_pos, end_upos, end_vpos);
    // qDebug() << "active_bivt2d::paint: angle_rel = " << angle_rel;

    if (angle_rel >= 0.0) {
        qp->setPen(QPen(QBrush(col_lgreen), 1, Qt::SolidLine));
        qp->setBrush(col_lgreen);
    }
    else {
        qp->setPen(QPen(QBrush(col_lblue), 1, Qt::SolidLine));
        qp->setBrush(col_lblue);
    }
    qp->drawPolygon(polygon);

    if (angle_rel >= 0.0) {
        qp->setPen(QPen(QBrush(Qt::darkGreen), 4, Qt::SolidLine));
    }
    else {
        qp->setPen(QPen(QBrush(Qt::darkBlue), 4, Qt::SolidLine));
    }
    qp->setBrush(Qt::NoBrush);
    qp->drawPath(anglePath(beg_pos, end_upos, end_vpos));


    // draw vector u
    qp->setPen(QPen(QBrush(Qt::black), 2, Qt::SolidLine));
    qp->setBrush(Qt::black);

    if (m_mouse_hover && !m_mouse_l_pressed) {
        qp->setPen(QPen(QBrush(col_green), 2, Qt::SolidLine)); // hover: green
        qp->setBrush(col_green);
    }
    if (m_mouse_hover && m_mouse_l_pressed) {
        qp->setPen(QPen(QBrush(col_red), 2, Qt::SolidLine)); // pressed: red
        qp->setBrush(col_red);
    }
    if (m_mouse_hover && !m_mouse_l_pressed &&
        (m_rotate_pt_mode == rotate_pt_mode::rotate_u ||
         m_rotate_pt_mode == rotate_pt_mode::rotate_both_uv)) {
        qp->setPen(QPen(QBrush(col_blue), 2, Qt::SolidLine)); // rotate: blue
        qp->setBrush(col_blue);
    }
    if (m_mouse_hover && m_mouse_l_pressed &&
        (m_rotate_pt_mode == rotate_pt_mode::rotate_u ||
         m_rotate_pt_mode == rotate_pt_mode::rotate_both_uv)) {
        qp->setPen(QPen(QBrush(col_blue), 2, Qt::SolidLine)); // rotate: blue
        qp->setBrush(col_blue);
    }

    qp->drawPath(arrowLine(beg_pos, end_upos));

    // from here on we want to draw with a small pen to get a pointy vector head
    QPen pen = qp->pen();
    pen.setWidth(1);
    qp->setPen(pen);
    qp->drawPath(arrowHead(beg_pos, end_upos));

    qp->drawPath(arrowLine(beg_pos, end_upos));


    // draw vector v
    qp->setPen(QPen(QBrush(Qt::black), 2, Qt::SolidLine));
    qp->setBrush(Qt::black);

    if (m_mouse_hover && !m_mouse_l_pressed) {
        qp->setPen(QPen(QBrush(col_green), 2, Qt::SolidLine)); // hover: green
        qp->setBrush(col_green);
    }
    if (m_mouse_hover && m_mouse_l_pressed) {
        qp->setPen(QPen(QBrush(col_red), 2, Qt::SolidLine)); // pressed: red
        qp->setBrush(col_red);
    }
    if (m_mouse_hover && !m_mouse_l_pressed &&
        (m_rotate_pt_mode == rotate_pt_mode::rotate_v ||
         m_rotate_pt_mode == rotate_pt_mode::rotate_both_uv)) {
        qp->setPen(QPen(QBrush(col_blue), 2, Qt::SolidLine)); // rotate: blue
        qp->setBrush(col_blue);
    }
    if (m_mouse_hover && m_mouse_l_pressed &&
        (m_rotate_pt_mode == rotate_pt_mode::rotate_v ||
         m_rotate_pt_mode == rotate_pt_mode::rotate_both_uv)) {
        qp->setPen(QPen(QBrush(col_blue), 2, Qt::SolidLine)); // rotate: blue
        qp->setBrush(col_blue);
    }

    qp->drawPath(arrowLine(beg_pos, end_vpos));

    // from here on we want to draw with a small pen to get a pointy vector head
    pen = qp->pen();
    pen.setWidth(1);
    qp->setPen(pen);
    qp->drawPath(arrowHead(beg_pos, end_vpos));


    // draw name, length, optionally unit of vectors and bivector

    qp->setFont(QFont("Helvetica", 14, QFont::Normal));
    QFontMetrics fm = qp->fontMetrics();
    // qp->setPen(Qt::white);  // set font color

    try {
        // make all available that is needed for GA
        using namespace hd::ga;
        using namespace hd::ga::ega;

        auto u = vec2d{scenePos_uend().x, scenePos_uend().y};
        auto v = vec2d{scenePos_vend().x, scenePos_vend().y};
        auto v_perp = normalize(reject_from(v, u)) * nrm(v);


        auto wdg_uv = wdg(u, v);
        // orientation relative to I_2d (= right-handed system of e1 and e2)
        auto orientation = sign(wdg_uv / I_2d);
        // qDebug() << "orientation " << orientation;

        auto dot_uv = dot(u, v);
        auto dot_height = dot_uv / nrm(u);
        // the factor orientiation is required to show consistent direction for
        // positive dot product "above u" within right-handed system
        auto v_perp_dot = normalize(reject_from(v, u)) * dot_height * orientation;

        auto gpr_uv = u * v;

        // fmt::println("u = {}", u);
        // fmt::println("v = {}", v);
        // fmt::println("v_perp = {}", v_perp);
        // fmt::println("nrm(u) = {}", nrm(u));
        // fmt::println("nrm(v) = {}", nrm(v));
        // fmt::println("nrm(dot_uv) = {}", nrm(dot_uv));
        // fmt::println("nrm(wdg_uv) = {}", nrm(wdg_uv));
        // fmt::println("nrm(gpr_uv) = {}", nrm(gpr_uv));
        // fmt::println("sum_dot_wdg = {}", nrm(dot_uv) + nrm(wdg_uv));
        // fmt::println("nrm_sq(dot_uv) = {}", nrm_sq(dot_uv));
        // fmt::println("nrm_sq(wdg_uv) = {}", nrm_sq(wdg_uv));
        // fmt::println("nrm_sq(gpr_uv) = {}", nrm_sq(gpr_uv));
        // fmt::println("sum_dot_wdg_sq = {}", nrm_sq(dot_uv) + nrm_sq(wdg_uv));

        auto v_perp_pos = QPointF(cs->x.au_to_w(v_perp.x), cs->y.au_to_w(v_perp.y));
        auto v_perp_dot_pos =
            QPointF(cs->x.au_to_w(v_perp_dot.x), cs->y.au_to_w(v_perp_dot.y));

        QPointF u_name_pos = beg_pos + 0.5 * (end_upos - beg_pos);
        QPointF v_name_pos = beg_pos + 0.5 * (end_vpos - beg_pos);

        QPointF dot_uv_name_pos = end_upos + 0.5 * (v_perp_dot_pos - beg_pos);
        QPointF wdg_uv_name_pos = u_name_pos + 0.5 * (end_vpos - beg_pos);
        QPointF gpr_uv_name_pos = v_perp_pos + 0.5 * (end_upos - beg_pos);

        // draw rectangle corresponding to u·v
        QPolygonF dot_rect;
        dot_rect << beg_pos << v_perp_dot_pos << v_perp_dot_pos + (end_upos - beg_pos)
                 << end_upos << beg_pos;
        qp->setPen(QPen(QBrush(col_lyel), 1, Qt::SolidLine));
        qp->setBrush(col_lyel);
        qp->drawPolygon(dot_rect);

        // draw rectangle corresponding to u * v
        QPolygonF gpr_rect;
        gpr_rect << beg_pos << v_perp_pos << v_perp_pos + (end_upos - beg_pos) << end_upos
                 << beg_pos;
        qp->setPen(QPen(QBrush(col_lblue), 1, Qt::SolidLine));
        qp->setBrush(col_lred);
        qp->drawPolygon(gpr_rect);

        // for shift of text to the right or to the left of the vector by using
        // contraction
        //
        // original version: not well-behaved if u, v or uv become too small)
        // auto dur = normalize(uv >> u);
        // auto dvl = normalize(v << uv);
        //
        // improved version: not well-behaved, only if either u or v become too small)
        auto dur = normalize(u << orientation * I_2d);
        auto dvl = normalize(orientation * I_2d >> v);
        auto dvpr = normalize(v_perp << orientation * I_2d);

        QString u_nrm = QString::number(nrm(u), 'f', 2);
        QString v_nrm = QString::number(nrm(v), 'f', 2);

        QString dot_uv_nrm = QString::number(nrm(dot_uv), 'f', 2);
        QString wdg_uv_nrm = QString::number(nrm(wdg_uv), 'f', 2);
        QString gpr_uv_nrm = QString::number(nrm(gpr_uv), 'f', 2);

        QString dot_uv_nrm_sq = QString::number(nrm_sq(dot_uv), 'f', 2);
        QString wdg_uv_nrm_sq = QString::number(nrm_sq(wdg_uv), 'f', 2);
        QString gpr_uv_nrm_sq = QString::number(nrm_sq(gpr_uv), 'f', 2);

        QString su = QString("|u| = ") + u_nrm;
        QString sv = QString("|v| = ") + v_nrm;

        QString sdot_uv = QString("u·v = s(");
        if (sign(dot_uv) < 0) {
            sdot_uv += QString("-");
        }
        sdot_uv += dot_uv_nrm + ")";

        QString suv = QString("u^v = ps(");
        if (orientation < 0) {
            suv += QString("-");
        }
        suv += wdg_uv_nrm + ")";

        QString sgpr_uv_mv = QString("u*v = mv_e(") + QString::number(gpr_uv.c0, 'f', 2) +
                             "," + QString::number(gpr_uv.c1, 'f', 2) + ")";
        QString sgpr_uv = QString("|u*v| = ");
        sgpr_uv += gpr_uv_nrm + "; ";
        sgpr_uv += QString("|u*v|² = |u·v|² + |u^v|² = ");
        sgpr_uv += dot_uv_nrm_sq + " + " + wdg_uv_nrm_sq + " = " + gpr_uv_nrm_sq;

        qp->setPen(QPen(QBrush(Qt::black), 1, Qt::SolidLine));

        const int offset = 25;
        qp->drawText(u_name_pos.x() - fm.horizontalAdvance(su) / 2 + int(offset * dur.x),
                     u_name_pos.y() + fm.height() / 2 - int(offset * dur.y), su);

        qp->drawText(v_name_pos.x() - fm.horizontalAdvance(sv) / 2 +
                         int(1.4 * offset * dvl.x),
                     v_name_pos.y() + fm.height() / 2 - int(1.4 * offset * dvl.y), sv);

        if (orientation >= 0) {
            qp->setPen(QPen(QBrush(Qt::darkRed), 1, Qt::SolidLine));
            qp->drawText(gpr_uv_name_pos.x() - fm.horizontalAdvance(sgpr_uv_mv) / 2 +
                             int(offset * dur.x),
                         gpr_uv_name_pos.y() - 1.3 * fm.height() + int(offset * dur.y),
                         sgpr_uv_mv);

            qp->setPen(QPen(QBrush(Qt::darkRed), 1, Qt::SolidLine));
            qp->drawText(gpr_uv_name_pos.x() - fm.horizontalAdvance(sgpr_uv) / 2 +
                             int(offset * dur.x),
                         gpr_uv_name_pos.y() + fm.height() / 2 + int(offset * dur.y),
                         sgpr_uv);
        }
        else {
            qp->setPen(QPen(QBrush(Qt::darkRed), 1, Qt::SolidLine));
            qp->drawText(gpr_uv_name_pos.x() - fm.horizontalAdvance(sgpr_uv_mv) / 2 +
                             int(offset * dur.x),
                         gpr_uv_name_pos.y() + fm.height() / 2 + int(offset * dur.y),
                         sgpr_uv_mv);

            qp->setPen(QPen(QBrush(Qt::darkRed), 1, Qt::SolidLine));
            qp->drawText(gpr_uv_name_pos.x() - fm.horizontalAdvance(sgpr_uv) / 2 +
                             int(offset * dur.x),
                         gpr_uv_name_pos.y() + 2.3 * fm.height() + int(offset * dur.y),
                         sgpr_uv);
        }

        if (orientation >= 0) {
            qp->setPen(QPen(QBrush(Qt::darkGreen), 1, Qt::SolidLine));
        }
        else {
            qp->setPen(QPen(QBrush(Qt::darkBlue), 1, Qt::SolidLine));
        }
        qp->drawText(wdg_uv_name_pos.x() - fm.horizontalAdvance(suv) / 2,
                     wdg_uv_name_pos.y() + fm.height() / 2, suv);

        qp->setPen(QPen(QBrush(Qt::darkYellow), 1, Qt::SolidLine));
        qp->drawText(dot_uv_name_pos.x() - fm.horizontalAdvance(sdot_uv) / 2 +
                         int(1.8 * offset * dvpr.x),
                     dot_uv_name_pos.y() + fm.height() / 2 - int(1.8 * offset * dvpr.y),
                     sdot_uv);
    }
    catch (std::exception const& e) {
        qDebug() << e.what();
    }

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    // draw shape (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawPath(shape());

    qp->restore();
}

QRectF active_bivt2d::boundingRect() const
{
    // give bounding box in item coordinate system
    QPointF beg_pos = QPointF(cs->x.a_to_w(0.0), cs->y.a_to_w(0.0));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(m_uend->scenePos().x), cs->y.a_to_w(m_uend->scenePos().y));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(m_vend->scenePos().x), cs->y.a_to_w(m_vend->scenePos().y));

    QPointF tip_pos = end_upos + end_vpos - beg_pos;

    // parallelogram corresponding to u^v
    QPolygonF polygon;

    polygon << beg_pos << end_upos << tip_pos << end_vpos << beg_pos;

    return polygon.boundingRect();
}

QPainterPath active_bivt2d::shape() const
{

    QPointF beg_pos = QPointF(cs->x.a_to_w(0.0), cs->y.a_to_w(0.0));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(m_uend->scenePos().x), cs->y.a_to_w(m_uend->scenePos().y));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(m_vend->scenePos().x), cs->y.a_to_w(m_vend->scenePos().y));

    QPainterPath path = vectorShape(beg_pos, end_upos);
    path += vectorShape(beg_pos, end_vpos);

    return path;
}

void active_bivt2d::setScenePos_uend(pt2d const& pos)
{
    if (pos != m_uend->scenePos()) {
        prepareGeometryChange();
        m_uend->setScenePos(pos);
    }
}

void active_bivt2d::setScenePos_vend(pt2d const& pos)
{
    if (pos != m_vend->scenePos()) {
        prepareGeometryChange();
        m_vend->setScenePos(pos);
    }
}

pt2d active_bivt2d::scenePos_uend() const { return m_uend->scenePos(); }
pt2d active_bivt2d::scenePos_vend() const { return m_vend->scenePos(); }


void active_bivt2d::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_bivt2d::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_bivt2d::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_bivt2d::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_bivt2d::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_bivt2d::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_bivt2d: Qt::LeftButton.";
        m_mouse_l_pressed = true;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_bivt2d: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }

    update();
    QGraphicsItem::mousePressEvent(event); // call default implementation
}

void active_bivt2d::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_bivt2d::mouseReleaseEvent.";

    // qDebug() << "active_bivt2d::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_bivt2d::scenePos_end():" << scenePos_end();

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt2d: Qt::LeftButton.";
        m_mouse_l_pressed = false;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt2d: Qt::RightButton.";
        m_mouse_r_pressed = false;
    }

    update();
    QGraphicsItem::mouseReleaseEvent(event); // call default implementation
}

void active_bivt2d::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_bivt2d::mouseMoveEvent.";

    if (m_mouse_l_pressed) {

        QPointF delta = event->scenePos() - event->lastScenePos();

        if (delta != QPointF(0, 0)) {

            // qDebug() << "scenePos():" << event->scenePos();
            // qDebug() << "lastScenePos():" << event->lastScenePos();
            // qDebug() << "delta:" << delta;

            switch (m_rotate_pt_mode) {
                case rotate_pt_mode::rotate_none:

                    m_uend->moveBy(delta.x(), delta.y());
                    m_vend->moveBy(delta.x(), delta.y());

                    m_uend->posChanged();
                    m_vend->posChanged();
                    break;

                case rotate_pt_mode::rotate_u: {
                    // get current turning angle
                    auto scnPos = vec2d(cs->x.w_to_au(event->scenePos().x()),
                                        cs->y.w_to_au(event->scenePos().y()));
                    auto lscnPos = vec2d(cs->x.w_to_au(event->lastScenePos().x()),
                                         cs->y.w_to_au(event->lastScenePos().y()));
                    auto cur_ang = angle(lscnPos, scnPos);
                    // fmt::println("ang = {}", rad2deg(cur_ang));

                    auto cur_u = vec2d(m_uend->scenePos().x, m_uend->scenePos().y);
                    auto rot = get_rotor(I_2d, cur_ang);
                    auto new_u = rotate(cur_u, rot);
                    m_uend->setScenePos(new_u);

                } break;

                case rotate_pt_mode::rotate_v: {
                    // get current turning angle
                    auto scnPos = vec2d(cs->x.w_to_au(event->scenePos().x()),
                                        cs->y.w_to_au(event->scenePos().y()));
                    auto lscnPos = vec2d(cs->x.w_to_au(event->lastScenePos().x()),
                                         cs->y.w_to_au(event->lastScenePos().y()));
                    auto cur_ang = angle(lscnPos, scnPos);
                    // fmt::println("ang = {}", rad2deg(cur_ang));

                    auto cur_v = vec2d(m_vend->scenePos().x, m_vend->scenePos().y);
                    auto rot = get_rotor(I_2d, cur_ang);
                    auto new_v = rotate(cur_v, rot);
                    m_vend->setScenePos(new_v);

                } break;

                case rotate_pt_mode::rotate_both_uv: {
                    // get current turning angle
                    auto scnPos = vec2d(cs->x.w_to_au(event->scenePos().x()),
                                        cs->y.w_to_au(event->scenePos().y()));
                    auto lscnPos = vec2d(cs->x.w_to_au(event->lastScenePos().x()),
                                         cs->y.w_to_au(event->lastScenePos().y()));
                    auto cur_ang = angle(lscnPos, scnPos);
                    // fmt::println("ang = {}", rad2deg(cur_ang));

                    auto cur_u = vec2d(m_uend->scenePos().x, m_uend->scenePos().y);
                    auto cur_v = vec2d(m_vend->scenePos().x, m_vend->scenePos().y);

                    auto rot = get_rotor(I_2d, cur_ang);

                    auto new_u = rotate(cur_u, rot);
                    auto new_v = rotate(cur_v, rot);

                    // Block all signals during update (avoid recursive updates)
                    m_uend->blockSignals(true);
                    m_vend->blockSignals(true);

                    m_uend->setScenePos(new_u);
                    m_vend->setScenePos(new_v);

                    // Re-enable signals - Qt will update automatically
                    m_uend->blockSignals(false);
                    m_vend->blockSignals(false);

                } break;

                default:
                    std::unreachable();
            }
        }
    }
}

void active_bivt2d::rotateModeChanged(rotate_pt_mode mode)
{

    if (mode != m_rotate_pt_mode) { // only adjust if actually changed
        m_rotate_pt_mode = mode;
    }
    update();
}
