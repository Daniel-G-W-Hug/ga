// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_reflection.hpp"
#include "active_common.hpp"

#include "ga/ga_ega.hpp"

#include <algorithm> // std::max
#include <array>

active_reflection::active_reflection(Coordsys* cs, w_Coordsys* wcs, active_pt2d* n1end,
                                     active_pt2d* n2end, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_n1end{n1end}, m_n2end{n2end}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    // setZValue(18);

    connect(wcs, &w_Coordsys::viewResized, m_n1end, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_n2end, &active_pt2d::viewChanged);

    connect(this, &active_reflection::pointsMoved, m_n1end, &active_pt2d::viewChanged);
    connect(this, &active_reflection::pointsMoved, m_n2end, &active_pt2d::viewChanged);

    connect(this, &active_reflection::viewMoved, m_n1end, &active_pt2d::posChanged);
    connect(this, &active_reflection::viewMoved, m_n2end, &active_pt2d::posChanged);
}

void active_reflection::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    // draw first normal vector and corresponding hyperplane

    QPointF beg_pos = QPointF(cs->x.au_to_w(0.0), cs->y.au_to_w(0.0));

    QPointF end_n1pos_nrm_scene = QPointF(scenePos_n1end().x, scenePos_n1end().y) /
                                  nrm(QPointF(scenePos_n1end().x, scenePos_n1end().y));
    QPointF end_n1pos =
        QPointF(cs->x.au_to_w(scenePos_n1end().x), cs->y.au_to_w(scenePos_n1end().y));
    QPointF end_n1pos_nrm = QPointF(cs->x.au_to_w(end_n1pos_nrm_scene.x()),
                                    cs->y.au_to_w(end_n1pos_nrm_scene.y()));

    qp->setPen(QPen(QBrush(col_lgreen), 2, Qt::SolidLine));
    qp->drawPath(arrowLine(beg_pos, end_n1pos_nrm));
    qp->setPen(Qt::gray);
    qp->drawPath(arrowLine(end_n1pos_nrm, end_n1pos));
    qp->setPen(QPen(QBrush(col_lgreen), 1, Qt::SolidLine));
    qp->setBrush(col_lgreen);
    // from here on we want to draw with a small pen to get a pointy vector head
    qp->drawPath(arrowHead(beg_pos, end_n1pos_nrm));
    qp->setPen(QPen(QBrush(col_lgreen), 2, Qt::SolidLine));
    qp->drawPath(planeLine(beg_pos, end_n1pos,
                           std::max(cs->x.widget_size(), cs->y.widget_size())));


    // draw second normal vector and corresponding hyperplane

    QPointF end_n2pos_nrm_scene = QPointF(scenePos_n2end().x, scenePos_n2end().y) /
                                  nrm(QPointF(scenePos_n2end().x, scenePos_n2end().y));
    QPointF end_n2pos =
        QPointF(cs->x.au_to_w(scenePos_n2end().x), cs->y.au_to_w(scenePos_n2end().y));
    QPointF end_n2pos_nrm = QPointF(cs->x.au_to_w(end_n2pos_nrm_scene.x()),
                                    cs->y.au_to_w(end_n2pos_nrm_scene.y()));

    qp->setPen(QPen(QBrush(col_lred), 2, Qt::SolidLine));
    qp->drawPath(arrowLine(beg_pos, end_n2pos_nrm));
    qp->setPen(Qt::gray);
    qp->drawPath(arrowLine(end_n2pos_nrm, end_n2pos));
    qp->setPen(QPen(QBrush(col_lred), 2, Qt::SolidLine));
    // from here on we want to draw with a small pen to get a pointy vector head
    qp->setPen(QPen(QBrush(col_lred), 1, Qt::SolidLine));
    qp->setBrush(col_lred);
    qp->drawPath(arrowHead(beg_pos, end_n2pos_nrm));
    qp->setPen(QPen(QBrush(col_lred), 2, Qt::SolidLine));
    qp->drawPath(planeLine(beg_pos, end_n2pos,
                           std::max(cs->x.widget_size(), cs->y.widget_size())));


    // let's do a little ga
    using hd::ga::vec2d;
    std::array<vec2d, 3> const v = {vec2d{0.5, 0.5}, vec2d{2.0, 0.5}, vec2d{2.5, 1.0}};
    vec2d n1(end_n1pos_nrm_scene.x(), end_n1pos_nrm_scene.y());
    vec2d n2(end_n2pos_nrm_scene.x(), end_n2pos_nrm_scene.y());
    std::vector<vec2d> vr, vrr;

    using hd::ga::ega::reflect_on;
    for (auto const& e : v) {
        vec2d v_tmp = reflect_on(e, n1);
        vr.push_back(v_tmp);
        vrr.push_back(reflect_on(v_tmp, n2));
    }

    QPainterPath org, refl1, refl2;
    org.moveTo(QPointF(cs->x.au_to_w(v[0].x), cs->y.au_to_w(v[0].y)));
    refl1.moveTo(QPointF(cs->x.au_to_w(vr[0].x), cs->y.au_to_w(vr[0].y)));
    refl2.moveTo(QPointF(cs->x.au_to_w(vrr[0].x), cs->y.au_to_w(vrr[0].y)));
    for (size_t i = 1; i < v.size(); ++i) {
        org.lineTo(QPointF(cs->x.au_to_w(v[i].x), cs->y.au_to_w(v[i].y)));
        refl1.lineTo(QPointF(cs->x.au_to_w(vr[i].x), cs->y.au_to_w(vr[i].y)));
        refl2.lineTo(QPointF(cs->x.au_to_w(vrr[i].x), cs->y.au_to_w(vrr[i].y)));
    }
    org.closeSubpath();
    refl1.closeSubpath();
    refl2.closeSubpath();

    qp->setPen(QPen(QBrush(col_lblue), 1, Qt::SolidLine)); // pressed: red
    qp->setBrush(col_lblue);
    qp->drawPath(org);

    qp->setPen(QPen(QBrush(col_lgreen), 1, Qt::SolidLine)); // pressed: red
    qp->setBrush(col_lgreen);
    qp->drawPath(refl1);

    qp->setPen(QPen(QBrush(col_lred), 1, Qt::SolidLine)); // pressed: red
    qp->setBrush(col_lred);
    qp->drawPath(refl2);

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    // draw shape (optional for testing)
    qp->setPen(col_yel);
    qp->setBrush(col_yel);
    qp->drawPath(shape());

    qp->restore();
}

QRectF active_reflection::boundingRect() const
{
    // give bounding box in item coordinate system

    // choose full drawing/clipping area here, because of complex alignment of reflected
    // shapes
    return QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                 cs->y.nmin() - cs->y.nmax());
}

QPainterPath active_reflection::shape() const
{
    QPointF beg_pos = QPointF(cs->x.au_to_w(0.0), cs->y.au_to_w(0.0));
    QPointF end_n1pos =
        QPointF(cs->x.au_to_w(scenePos_n1end().x), cs->y.au_to_w(scenePos_n1end().y));
    QPointF end_n2pos =
        QPointF(cs->x.au_to_w(scenePos_n2end().x), cs->y.au_to_w(scenePos_n2end().y));

    QPainterPath path = vectorShape(beg_pos, end_n1pos);
    path += vectorShape(beg_pos, end_n2pos);

    return path;
}

void active_reflection::setScenePos_n1end(pt2d const& pos)
{
    if (pos != m_n1end->scenePos()) {
        prepareGeometryChange();
        m_n1end->setScenePos(pos);
    }
}

void active_reflection::setScenePos_n2end(pt2d const& pos)
{
    if (pos != m_n2end->scenePos()) {
        prepareGeometryChange();
        m_n2end->setScenePos(pos);
    }
}

void active_reflection::setScenePos_n1end_wo_update(pt2d const& pos)
{
    if (pos != m_n1end->scenePos()) {
        m_n1end->setScenePos(pos);
    }
}

void active_reflection::setScenePos_n2end_wo_update(pt2d const& pos)
{
    if (pos != m_n2end->scenePos()) {
        m_n2end->setScenePos(pos);
    }
}

pt2d active_reflection::scenePos_n1end() const { return m_n1end->scenePos(); }

pt2d active_reflection::scenePos_n2end() const { return m_n2end->scenePos(); }


void active_reflection::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_reflection::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_reflection::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_reflection::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_reflection::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_reflection::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_reflection: Qt::LeftButton.";
        m_mouse_l_pressed = true;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_reflection: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }
    update();
}

void active_reflection::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_reflection::mouseReleaseEvent.";

    // qDebug() << "active_reflection::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_reflection::scenePos_n1end():" << scenePos_n1end();

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt2d: Qt::LeftButton.";
        m_mouse_l_pressed = false;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt2d: Qt::RightButton.";
        m_mouse_r_pressed = false;
    }
    update();
}

void active_reflection::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{

    // qDebug() << "active_reflection::mouseMoveEvent.";

    if (m_mouse_l_pressed) {

        QPointF delta = event->scenePos() - event->lastScenePos();

        if (delta != QPointF(0, 0)) {

            // qDebug() << "scenePos():" << event->scenePos();
            // qDebug() << "lastScenePos():" << event->lastScenePos();
            // qDebug() << "delta:" << delta;

            // get current turning angle
            auto scnPos = vec2d(cs->x.w_to_au(event->scenePos().x()),
                                cs->y.w_to_au(event->scenePos().y()));
            auto lscnPos = vec2d(cs->x.w_to_au(event->lastScenePos().x()),
                                 cs->y.w_to_au(event->lastScenePos().y()));
            auto cur_ang = angle(lscnPos, scnPos);
            // fmt::println("ang = {}", rad2deg(cur_ang));

            // get normal vectors and transform them with current turning angle
            auto cur_n1 = vec2d(scenePos_n1end());
            auto cur_n2 = vec2d(scenePos_n2end());

            auto new_n1 = rotate(cur_n1, get_rotor(I_2d, cur_ang));
            auto new_n2 = rotate(cur_n2, get_rotor(I_2d, cur_ang));

            setScenePos_n1end_wo_update(new_n1);
            setScenePos_n2end_wo_update(new_n2);

            emit pointsMoved();
        }
    }
}
