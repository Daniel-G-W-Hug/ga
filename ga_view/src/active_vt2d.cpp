// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_vt2d.hpp"
#include "active_common.hpp"

active_vt2d::active_vt2d(Coordsys* cs, w_Coordsys* wcs, active_pt2d* beg,
                         active_pt2d* end, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_beg{beg}, m_end{end}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    // setZValue(18);

    connect(wcs, &w_Coordsys::viewResized, m_beg, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_end, &active_pt2d::viewChanged);

    connect(this, &active_vt2d::viewMoved, m_beg, &active_pt2d::posChanged);
    connect(this, &active_vt2d::viewMoved, m_end, &active_pt2d::posChanged);
}

void active_vt2d::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                        QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    // clipping area is active area of coordsys
    qp->setClipRect(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                          cs->y.nmin() - cs->y.nmax()));


    // draw in item coordinate system
    qp->save();

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

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x), cs->y.a_to_w(scenePos_beg().y));
    QPointF end_pos =
        QPointF(cs->x.a_to_w(scenePos_end().x), cs->y.a_to_w(scenePos_end().y));

    QPen pen = qp->pen();
    pen.setWidth(2);
    qp->drawPath(arrowLine(beg_pos, end_pos));

    // from here on we want to draw with a small pen to get a pointy vector head
    pen.setWidth(1);
    qp->drawPath(arrowHead(beg_pos, end_pos));

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

QRectF active_vt2d::boundingRect() const
{
    // give bounding box in item coordinate system
    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x), cs->y.a_to_w(scenePos_beg().y));
    QPointF end_pos =
        QPointF(cs->x.a_to_w(scenePos_end().x), cs->y.a_to_w(scenePos_end().y));
    return QRectF(beg_pos, end_pos).normalized();
}

QPainterPath active_vt2d::shape() const
{
    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x), cs->y.a_to_w(scenePos_beg().y));
    QPointF end_pos =
        QPointF(cs->x.a_to_w(scenePos_end().x), cs->y.a_to_w(scenePos_end().y));
    return vectorShape(beg_pos, end_pos);
}

void active_vt2d::setScenePos_beg(pt2d const& pos)
{
    if (pos != m_beg->scenePos()) {
        prepareGeometryChange();
        m_beg->setScenePos(pos);
    }
}

void active_vt2d::setScenePos_end(pt2d const& pos)
{
    if (pos != m_end->scenePos()) {
        prepareGeometryChange();
        m_end->setScenePos(pos);
    }
}

pt2d active_vt2d::scenePos_beg() const { return m_beg->scenePos(); }
pt2d active_vt2d::scenePos_end() const { return m_end->scenePos(); }


void active_vt2d::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_vt2d::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_vt2d::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_vt2d::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_vt2d::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_vt2d::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_vt2d: Qt::LeftButton.";
        m_mouse_l_pressed = true;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_vt2d: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }

    update();
    QGraphicsItem::mousePressEvent(event); // call default implementation
}

void active_vt2d::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_vt2d::mouseReleaseEvent.";

    // qDebug() << "active_vt2d::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_vt2d::scenePos_end():" << scenePos_end();

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

void active_vt2d::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_vt2d::mouseMoveEvent.";

    if (m_mouse_l_pressed) {

        QPointF delta = event->scenePos() - event->lastScenePos();

        if (delta != QPointF(0, 0)) {

            // qDebug() << "scenePos():" << event->scenePos();
            // qDebug() << "lastScenePos():" << event->lastScenePos();
            // qDebug() << "delta:" << delta;

            m_beg->moveBy(delta.x(), delta.y());
            m_end->moveBy(delta.x(), delta.y());

            emit viewMoved();
        }
    }
}
