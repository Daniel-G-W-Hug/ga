//
// author: Daniel Hug, 2024
//

#include "active_vec.hpp"
#include "active_common.hpp"

active_vec::active_vec(Coordsys* cs, w_Coordsys* wcs, active_pt* beg, active_pt* end,
                       QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_beg{beg}, m_end{end}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    // setZValue(18);

    connect(wcs, &w_Coordsys::viewResized, m_beg, &active_pt::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_end, &active_pt::viewChanged);

    connect(this, &active_vec::viewMoved, m_beg, &active_pt::posChanged);
    connect(this, &active_vec::viewMoved, m_end, &active_pt::posChanged);
}

void active_vec::paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                       QWidget* widget)
{

    // clipping area is active area of coordsys
    painter->setClipRect(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax()));


    // draw in item coordinate system
    painter->save();

    painter->setPen(QPen(QBrush(Qt::black), 2, Qt::SolidLine));
    painter->setBrush(Qt::black);

    if (m_mouse_hover && !m_mouse_l_pressed) {
        painter->setPen(QPen(QBrush(col_green), 2, Qt::SolidLine)); // hover: green
        painter->setBrush(col_green);
    }
    if (m_mouse_hover && m_mouse_l_pressed) {
        painter->setPen(QPen(QBrush(col_red), 2, Qt::SolidLine)); // pressed: red
        painter->setBrush(col_red);
    }

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x()), cs->y.a_to_w(scenePos_beg().y()));
    QPointF end_pos =
        QPointF(cs->x.a_to_w(scenePos_end().x()), cs->y.a_to_w(scenePos_end().y()));

    QPen pen = painter->pen();
    pen.setWidth(2);
    painter->drawPath(arrowLine(beg_pos, end_pos));

    // from here on we want to draw with a small pen to get a pointy vector head
    pen.setWidth(1);
    painter->drawPath(arrowHead(beg_pos, end_pos));

    // draw bounding box (optional for testing)
    // painter->setPen(col_yel);
    // painter->setBrush(col_yel);
    // painter->drawRect(boundingRect());

    // draw shape (optional for testing)
    // painter->setPen(col_yel);
    // painter->setBrush(col_yel);
    // painter->drawPath(shape());

    painter->restore();
}

QRectF active_vec::boundingRect() const
{
    // give bounding box in item coordinate system
    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x()), cs->y.a_to_w(scenePos_beg().y()));
    QPointF end_pos =
        QPointF(cs->x.a_to_w(scenePos_end().x()), cs->y.a_to_w(scenePos_end().y()));
    return QRectF(beg_pos, end_pos).normalized();
}

QPainterPath active_vec::shape() const
{
    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x()), cs->y.a_to_w(scenePos_beg().y()));
    QPointF end_pos =
        QPointF(cs->x.a_to_w(scenePos_end().x()), cs->y.a_to_w(scenePos_end().y()));
    return vectorShape(beg_pos, end_pos);
}

void active_vec::setScenePos_beg(QPointF const& pos)
{
    if (pos != m_beg->scenePos()) {
        prepareGeometryChange();
        m_beg->setScenePos(pos);
    }
}

void active_vec::setScenePos_end(QPointF const& pos)
{
    if (pos != m_end->scenePos()) {
        prepareGeometryChange();
        m_end->setScenePos(pos);
    }
}

QPointF active_vec::scenePos_beg() const { return m_beg->scenePos(); }
QPointF active_vec::scenePos_end() const { return m_end->scenePos(); }


void active_vec::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_vec::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_vec::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_vec::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_vec::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_vec::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_vec: Qt::LeftButton.";
        m_mouse_l_pressed = true;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_vec: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }

    update();
    QGraphicsItem::mousePressEvent(event); // call default implementation
}

void active_vec::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_vec::mouseReleaseEvent.";

    // qDebug() << "active_vec::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_vec::scenePos_end():" << scenePos_end();

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt: Qt::LeftButton.";
        m_mouse_l_pressed = false;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt: Qt::RightButton.";
        m_mouse_r_pressed = false;
    }

    update();
    QGraphicsItem::mouseReleaseEvent(event); // call default implementation
}

void active_vec::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_vec::mouseMoveEvent.";

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
