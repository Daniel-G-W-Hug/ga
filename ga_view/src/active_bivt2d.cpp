//
// author: Daniel Hug, 2024
//

#include "active_bivt2d.hpp"
#include "active_common.hpp"


active_bivt2d::active_bivt2d(Coordsys* cs, w_Coordsys* wcs, active_pt2d* beg,
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
}

void active_bivt2d::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                          QWidget* widget)
{

    // clipping area is active area of coordsys
    qp->setClipRect(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                          cs->y.nmin() - cs->y.nmax()));

    // draw in item coordinate system
    qp->save();

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(m_beg->scenePos().x), cs->y.a_to_w(m_beg->scenePos().y));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(m_uend->scenePos().x), cs->y.a_to_w(m_uend->scenePos().y));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(m_vend->scenePos().x), cs->y.a_to_w(m_vend->scenePos().y));

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
    QPointF beg_pos =
        QPointF(cs->x.a_to_w(m_beg->scenePos().x), cs->y.a_to_w(m_beg->scenePos().y));
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

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(m_beg->scenePos().x), cs->y.a_to_w(m_beg->scenePos().y));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(m_uend->scenePos().x), cs->y.a_to_w(m_uend->scenePos().y));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(m_vend->scenePos().x), cs->y.a_to_w(m_vend->scenePos().y));

    QPainterPath path = vectorShape(beg_pos, end_upos);
    path += vectorShape(beg_pos, end_vpos);

    return path;
}

void active_bivt2d::setScenePos_beg(pt2d const& pos)
{
    if (pos != m_beg->scenePos()) {
        prepareGeometryChange();
        m_beg->setScenePos(pos);
    }
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

pt2d active_bivt2d::scenePos_beg() const { return m_beg->scenePos(); }
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

            m_beg->moveBy(delta.x(), delta.y());
            m_uend->moveBy(delta.x(), delta.y());
            m_vend->moveBy(delta.x(), delta.y());

            m_beg->posChanged();
            m_uend->posChanged();
            m_vend->posChanged();
        }
    }
}
