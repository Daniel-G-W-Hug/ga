//
// author: Daniel Hug, 2024
//

#include "active_pt.hpp"
#include "active_common.hpp"

active_pt::active_pt(Coordsys* cs, w_Coordsys* wcs, QPointF const& pos,
                     QGraphicsItem* parent) : QGraphicsItem(parent), cs{cs}, m_pos{pos}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    connect(wcs, &w_Coordsys::viewResized, this, &active_pt::viewChanged);

    setPos(cs->x.a_to_w(m_pos.x()),
           cs->y.a_to_w(m_pos.y())); // set item to scene coordinates

    setZValue(100); // active points should always be on top
}

void active_pt::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                      QWidget* widget)
{

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

QRectF active_pt::boundingRect() const
{
    return QRectF(QPointF(-RADIUS, -RADIUS), QPointF(RADIUS, RADIUS));
}

QPainterPath active_pt::shape() const
{
    QPainterPath path;
    path.addEllipse(boundingRect());
    return path;
}

void active_pt::setScenePos(QPointF const& pos)
{
    // qDebug() << "active_pt::setScenePos called.";

    if (m_pos != pos) {

        // qDebug() << "active_pt::setScenePos changed.";

        prepareGeometryChange();
        m_pos = pos;
    }
}

QPointF active_pt::scenePos() { return m_pos; }


void active_pt::viewChanged()
{
    // qDebug() << "active_pt: viewChanged() received.";

    // view changed by external influence, set to m_pos
    setPos(cs->x.a_to_w(m_pos.x()), cs->y.a_to_w(m_pos.y()));
}

void active_pt::posChanged()
{
    // qDebug() << "active_pt: posChanged() received.";

    // position changed by external influence, update m_pos
    QPointF npos = pos();
    m_pos = QPointF(cs->x.w_to_a(npos.x()), cs->y.w_to_a(npos.y()));
}

void active_pt::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    // qDebug() << "active_pt::hoverEnterEvent.";
    m_mouse_hover = true;

    update();
    QGraphicsItem::hoverEnterEvent(event);
}

void active_pt::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    // qDebug() << "active_pt::hoverLeaveEvent.";
    m_mouse_hover = false;

    update();
    QGraphicsItem::hoverLeaveEvent(event);
}

void active_pt::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_pt::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt: Qt::LeftButton.";
        m_mouse_l_pressed = true;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }

    update();
    QGraphicsItem::mousePressEvent(event);
}

void active_pt::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_pt::mouseReleaseEvent.";
    // qDebug() << "active_vec::scenePos:" << scenePos();

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt: Qt::LeftButton.";
        m_mouse_l_pressed = false;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt: Qt::RightButton.";
        m_mouse_r_pressed = false;
    }

    update();
    QGraphicsItem::mouseReleaseEvent(event);
}

void active_pt::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_pt::mouseMoveEvent.";

    if (m_mouse_l_pressed) {

        QPointF delta = event->scenePos() - event->lastScenePos();

        // qDebug() << "scenePos():" << event->scenePos();
        // qDebug() << "lastScenePos():" << event->lastScenePos();
        // qDebug() << "delta:" << delta;
        // qDebug() << "pos():" << event->pos();

        if (delta != QPointF(0, 0)) {
            moveBy(delta.x(), delta.y());
            posChanged();
        }
    }

    update();
    QGraphicsItem::mouseMoveEvent(event); // move the item normally
}