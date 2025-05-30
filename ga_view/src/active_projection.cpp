// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_projection.hpp"
#include "active_common.hpp"

#include "ga/ga_ega.hpp"

active_projection::active_projection(Coordsys* cs, w_Coordsys* wcs, active_pt2d* beg,
                                     active_pt2d* uend, active_pt2d* vend,
                                     QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_beg{beg}, m_uend{uend}, m_vend{vend}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    // setZValue(18);

    connect(wcs, &w_Coordsys::viewResized, m_beg, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_uend, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_vend, &active_pt2d::viewChanged);

    connect(this, &active_projection::viewMoved, m_beg, &active_pt2d::posChanged);
    connect(this, &active_projection::viewMoved, m_uend, &active_pt2d::posChanged);
    connect(this, &active_projection::viewMoved, m_vend, &active_pt2d::posChanged);
}

void active_projection::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x), cs->y.a_to_w(scenePos_beg().y));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(scenePos_uend().x), cs->y.a_to_w(scenePos_uend().y));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(scenePos_vend().x), cs->y.a_to_w(scenePos_vend().y));

    try {

        // make all available that is needed for GA
        using namespace hd::ga;
        using namespace hd::ga::ega;

        auto u = vec2d{scenePos_uend().x - scenePos_beg().x,
                       scenePos_uend().y - scenePos_beg().y};
        auto v = vec2d{scenePos_vend().x - scenePos_beg().x,
                       scenePos_vend().y - scenePos_beg().y};

        auto u_inv = vec2d{};
        if (nrm_sq(u) != 0.0) {
            u_inv = inv(u); // this could throw if nrm_sq(u) == 0.0
        }
        auto v_par = dot(u, v) * u_inv;

        QPointF vpar = QPointF(scenePos_beg().x + v_par.x, scenePos_beg().y + v_par.y);

        // qDebug() << "vpar:  " << vpar;
        // qDebug() << "scenePos_beg(): " << scenePos_beg();

        QPointF end_vpar_pos = QPointF(cs->x.a_to_w(vpar.x()), cs->y.a_to_w(vpar.y()));

        qreal angle_rel = angle(u, v);
        // qDebug() << "active_projection::paint: angle_rel = " << angle_rel;

        if (angle_rel >= 0.0) {
            qp->setPen(QPen(QBrush(Qt::darkGreen), 4, Qt::SolidLine));
        }
        else {
            qp->setPen(QPen(QBrush(Qt::darkBlue), 4, Qt::SolidLine));
        }
        qp->setBrush(Qt::NoBrush);
        qp->drawPath(anglePath(beg_pos, end_upos, end_vpos));


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

        // draw vectors u and v
        qp->drawPath(arrowLine(beg_pos, end_upos));

        // from here on we want to draw with a small pen to get a pointy vector head
        QPen pen = qp->pen();
        pen.setWidth(1);
        qp->setPen(pen);
        qp->drawPath(arrowHead(beg_pos, end_upos));

        pen.setWidth(2);
        qp->setPen(pen);
        qp->drawPath(arrowLine(beg_pos, end_vpos));

        // from here on we want to draw with a small pen to get a pointy vector head
        pen.setWidth(1);
        qp->setPen(pen);
        qp->drawPath(arrowHead(beg_pos, end_vpos));


        // draw vectors v_par and v_perp
        qp->setPen(QPen(QBrush(col_blue), 3, Qt::SolidLine));
        qp->drawPath(arrowLine(beg_pos, end_vpar_pos));
        pen = qp->pen();
        pen.setWidth(1);
        qp->setPen(pen);
        qp->setBrush(col_blue);
        qp->drawPath(arrowHead(beg_pos, end_vpar_pos));

        qp->setPen(QPen(QBrush(col_blue), 3, Qt::SolidLine));
        qp->drawPath(arrowLine(end_vpar_pos, end_vpos));
        pen = qp->pen();
        pen.setWidth(1);
        qp->setPen(pen);
        qp->setBrush(col_blue);
        qp->drawPath(arrowHead(end_vpar_pos, end_vpos));
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

QRectF active_projection::boundingRect() const
{
    // ToDo: check for right bounding box (taken from bivec; tip_pos to be replaced)

    // give bounding box in item coordinate system
    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x), cs->y.a_to_w(scenePos_beg().y));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(scenePos_uend().x), cs->y.a_to_w(scenePos_uend().y));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(scenePos_vend().x), cs->y.a_to_w(scenePos_vend().y));

    QPointF tip_pos = end_upos + end_vpos - beg_pos;

    // parallelogram corresponding to u^v
    QPolygonF polygon;

    polygon << beg_pos << end_upos << tip_pos << end_vpos << beg_pos;

    return polygon.boundingRect();
}

QPainterPath active_projection::shape() const
{

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x), cs->y.a_to_w(scenePos_beg().y));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(scenePos_uend().x), cs->y.a_to_w(scenePos_uend().y));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(scenePos_vend().x), cs->y.a_to_w(scenePos_vend().y));

    QPainterPath path = vectorShape(beg_pos, end_upos);
    path += vectorShape(beg_pos, end_vpos);

    return path;
}

void active_projection::setScenePos_beg(pt2d const& pos)
{
    if (pos != m_beg->scenePos()) {
        prepareGeometryChange();
        m_beg->setScenePos(pos);
    }
}

void active_projection::setScenePos_uend(pt2d const& pos)
{
    if (pos != m_uend->scenePos()) {
        prepareGeometryChange();
        m_uend->setScenePos(pos);
    }
}

void active_projection::setScenePos_vend(pt2d const& pos)
{
    if (pos != m_vend->scenePos()) {
        prepareGeometryChange();
        m_vend->setScenePos(pos);
    }
}

pt2d active_projection::scenePos_beg() const { return m_beg->scenePos(); }
pt2d active_projection::scenePos_uend() const { return m_uend->scenePos(); }
pt2d active_projection::scenePos_vend() const { return m_vend->scenePos(); }


void active_projection::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_projection::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_projection::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_projection::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_projection::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_projection::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_projection: Qt::LeftButton.";
        m_mouse_l_pressed = true;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_projection: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }

    update();
    QGraphicsItem::mousePressEvent(event); // call default implementation
}

void active_projection::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_projection::mouseReleaseEvent.";

    // qDebug() << "active_projection::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_projection::scenePos_end():" << scenePos_end();

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

void active_projection::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_projection::mouseMoveEvent.";

    if (m_mouse_l_pressed) {

        QPointF delta = event->scenePos() - event->lastScenePos();

        if (delta != QPointF(0, 0)) {

            // qDebug() << "scenePos():" << event->scenePos();
            // qDebug() << "lastScenePos():" << event->lastScenePos();
            // qDebug() << "delta:" << delta;

            m_beg->moveBy(delta.x(), delta.y());
            m_uend->moveBy(delta.x(), delta.y());
            m_vend->moveBy(delta.x(), delta.y());

            emit viewMoved();
        }
    }
}
