// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_pt2d.hpp"
#include "active_common.hpp"

active_pt2d::active_pt2d(Coordsys* cs, w_Coordsys* wcs, pt2d& pos,
                         QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, m_pos{pos}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    connect(wcs, &w_Coordsys::viewResized, this, &active_pt2d::viewChanged);

    setPos(cs->x.a_to_w(m_pos.x),
           cs->y.a_to_w(m_pos.y)); // set item to scene coordinates

    setZValue(100); // active points should always be on top
}

void active_pt2d::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    qp->setPen(col_blue);
    qp->setBrush(col_blue); // selectable: blue (default)

    if (m_mouse_hover && !m_mouse_l_pressed) {
        qp->setPen(col_green);
        qp->setBrush(col_green); // hover: green
    }
    if (m_mouse_hover && m_mouse_l_pressed) {
        qp->setPen(col_red);
        qp->setBrush(col_red); // selected: red
    }

    qp->drawEllipse(QRectF(QPointF(-RADIUS, -RADIUS), QPointF(RADIUS, RADIUS)));

    qp->restore();
}

QRectF active_pt2d::boundingRect() const
{
    return QRectF(QPointF(-RADIUS, -RADIUS), QPointF(RADIUS, RADIUS));
}

QPainterPath active_pt2d::shape() const
{
    QPainterPath path;
    path.addEllipse(boundingRect());
    return path;
}

void active_pt2d::setScenePos(pt2d const& pos)
{
    // qDebug() << "active_pt2d::setScenePos called.";

    if (m_pos != pos) {

        // qDebug() << "active_pt2d::setScenePos changed.";

        prepareGeometryChange();
        m_pos = pos;
        setPos(cs->x.au_to_w(m_pos.x), cs->y.au_to_w(m_pos.y));
        emit pointMoved();
    }
}

void active_pt2d::setScenePos_wo_update(pt2d const& pos)
{

    // just store, no update (for use when several points need to be updated)

    // qDebug() << "active_pt2d::setScenePos called.";
    if (m_pos != pos) {
        // qDebug() << "active_pt2d::setScenePos changed.";
        m_pos = pos;
    }
}

pt2d active_pt2d::scenePos() const { return m_pos; }

void active_pt2d::viewChanged()
{
    // qDebug() << "active_pt2d: viewChanged() received.";

    // view changed by external influence, set to m_pos
    setPos(cs->x.au_to_w(m_pos.x), cs->y.au_to_w(m_pos.y));
    update();
}

void active_pt2d::posChanged()
{
    // qDebug() << "active_pt2d: posChanged() received.";

    // position changed by external influence, update m_pos
    QPointF npos = pos();
    m_pos = pt2d(cs->x.w_to_au(npos.x()), cs->y.w_to_au(npos.y()));
    emit pointMoved();
}

void active_pt2d::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_pt2d::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_pt2d::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_pt2d::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_pt2d::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_pt2d::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt2d: Qt::LeftButton.";
        m_mouse_l_pressed = true;
        update();
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt2d: Qt::RightButton.";
        m_mouse_r_pressed = true;
        update();
    }
}

void active_pt2d::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_pt2d::mouseReleaseEvent.";
    // qDebug() << "active_vt2d::scenePos:" << scenePos();

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt2d: Qt::LeftButton.";
        m_mouse_l_pressed = false;
        update();
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt2d: Qt::RightButton.";
        m_mouse_r_pressed = false;
        update();
    }
}

void active_pt2d::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_pt2d::mouseMoveEvent.";

    if (m_mouse_l_pressed) {

        QPointF delta = event->scenePos() - event->lastScenePos();

        // qDebug() << "scenePos():" << event->scenePos();
        // qDebug() << "lastScenePos():" << event->lastScenePos();
        // qDebug() << "delta:" << delta;
        // qDebug() << "pos():" << event->pos();

        if (delta != QPointF(0, 0)) {
            moveBy(delta.x(), delta.y());
            posChanged();
            update();
        }
    }
}