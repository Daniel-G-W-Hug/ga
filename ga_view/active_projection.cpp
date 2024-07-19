//
// author: Daniel Hug, 2024
//

#include "active_projection.hpp"
#include "active_common.hpp"


active_projection::active_projection(Coordsys* cs, w_Coordsys* wcs, active_pt* beg,
                                     active_pt* uend, active_pt* vend,
                                     QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_beg{beg}, m_uend{uend}, m_vend{vend}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    // setZValue(18);

    connect(wcs, &w_Coordsys::viewResized, m_beg, &active_pt::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_uend, &active_pt::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_vend, &active_pt::viewChanged);

    connect(this, &active_projection::viewMoved, m_beg, &active_pt::posChanged);
    connect(this, &active_projection::viewMoved, m_uend, &active_pt::posChanged);
    connect(this, &active_projection::viewMoved, m_vend, &active_pt::posChanged);
}

void active_projection::paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                              QWidget* widget)
{

    // clipping area is active area of coordsys
    painter->setClipRect(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax()));

    // draw in item coordinate system
    painter->save();

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x()), cs->y.a_to_w(scenePos_beg().y()));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(scenePos_uend().x()), cs->y.a_to_w(scenePos_uend().y()));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(scenePos_vend().x()), cs->y.a_to_w(scenePos_vend().y()));

    QPointF u = scenePos_uend() - scenePos_beg();
    QPointF v = scenePos_vend() - scenePos_beg();

    QPointF u_inv = u / (u.x() * u.x() + u.y() * u.y());

    qreal nrm_prod =
        sqrt(u.x() * u.x() + u.y() * u.y()) * sqrt(v.x() * v.x() + v.y() * v.y());

    qreal dot = u.x() * v.x() + u.y() * v.y();    // dot(v,u)
    qreal wdg = -(u.x() * v.y() - u.y() * v.x()); // wdg(v,u)

    QPointF vpar = scenePos_beg() + dot * u_inv;
    QPointF vperp = QPointF(wdg * u_inv.y(), -wdg * u_inv.x()); // wdg calc. manually

    // qDebug() << "vpar:  " << vpar;
    // qDebug() << "vperp: " << vperp;
    // qDebug() << "scenePos_beg(): " << scenePos_beg();

    QPointF end_vpar_pos = QPointF(cs->x.a_to_w(vpar.x()), cs->y.a_to_w(vpar.y()));
    QPointF end_vperp_pos = QPointF(cs->x.a_to_w(vperp.x()), cs->y.a_to_w(vperp.y()));

    // The sign has to be reversed here, since device coordinates are in a left-handed
    // system. The angle calculation itself is done in a classical right-handed system.
    //
    // Be aware: Depending on aspect ratio of x- vs. y-axis calculated angles will change,
    // if calculated from device coordinate input. However, they are visually consistent.
    //
    // Actual angle calculations for mathematical/physical purposes must use logical
    // coordinates as input for the calculation exclusively!
    qreal angle_rel = -angle_between_lines(beg_pos, end_upos, end_vpos);
    // qDebug() << "active_projection::paint: angle_rel = " << angle_rel;

    if (angle_rel >= 0.0) {
        painter->setPen(QPen(QBrush(Qt::darkGreen), 4, Qt::SolidLine));
    }
    else {
        painter->setPen(QPen(QBrush(Qt::darkBlue), 4, Qt::SolidLine));
    }
    painter->setBrush(Qt::NoBrush);
    painter->drawPath(anglePath(beg_pos, end_upos, end_vpos));

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

    // draw vectors u and v
    painter->drawPath(arrowLine(beg_pos, end_upos));

    // from here on we want to draw with a small pen to get a pointy vector head
    QPen pen = painter->pen();
    pen.setWidth(1);
    painter->setPen(pen);
    painter->drawPath(arrowHead(beg_pos, end_upos));

    pen.setWidth(2);
    painter->setPen(pen);
    painter->drawPath(arrowLine(beg_pos, end_vpos));

    // from here on we want to draw with a small pen to get a pointy vector head
    pen.setWidth(1);
    painter->setPen(pen);
    painter->drawPath(arrowHead(beg_pos, end_vpos));


    // draw vectors v_par and v_perp
    painter->setPen(QPen(QBrush(col_blue), 3, Qt::SolidLine));
    painter->drawPath(arrowLine(beg_pos, end_vpar_pos));
    pen = painter->pen();
    pen.setWidth(1);
    painter->setPen(pen);
    painter->setBrush(col_blue);
    painter->drawPath(arrowHead(beg_pos, end_vpar_pos));

    painter->setPen(QPen(QBrush(col_blue), 3, Qt::SolidLine));
    painter->drawPath(arrowLine(end_vpar_pos, end_vpos));
    pen = painter->pen();
    pen.setWidth(1);
    painter->setPen(pen);
    painter->setBrush(col_blue);
    painter->drawPath(arrowHead(end_vpar_pos, end_vpos));

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

QRectF active_projection::boundingRect() const
{
    // ToDo: check for right bounding box (taken from bivec; tip_pos to be replaced)

    // give bounding box in item coordinate system
    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x()), cs->y.a_to_w(scenePos_beg().y()));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(scenePos_uend().x()), cs->y.a_to_w(scenePos_uend().y()));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(scenePos_vend().x()), cs->y.a_to_w(scenePos_vend().y()));

    QPointF tip_pos = end_upos + end_vpos - beg_pos;

    // parallelogram corresponding to u^v
    QPolygonF polygon;

    polygon << beg_pos << end_upos << tip_pos << end_vpos << beg_pos;

    return polygon.boundingRect();
}

QPainterPath active_projection::shape() const
{

    QPointF beg_pos =
        QPointF(cs->x.a_to_w(scenePos_beg().x()), cs->y.a_to_w(scenePos_beg().y()));
    QPointF end_upos =
        QPointF(cs->x.a_to_w(scenePos_uend().x()), cs->y.a_to_w(scenePos_uend().y()));
    QPointF end_vpos =
        QPointF(cs->x.a_to_w(scenePos_vend().x()), cs->y.a_to_w(scenePos_vend().y()));

    QPainterPath path = vectorShape(beg_pos, end_upos);
    path += vectorShape(beg_pos, end_vpos);

    return path;
}

void active_projection::setScenePos_beg(QPointF const& pos)
{
    if (pos != m_beg->scenePos()) {
        prepareGeometryChange();
        m_beg->setScenePos(pos);
    }
}

void active_projection::setScenePos_uend(QPointF const& pos)
{
    if (pos != m_uend->scenePos()) {
        prepareGeometryChange();
        m_uend->setScenePos(pos);
    }
}

void active_projection::setScenePos_vend(QPointF const& pos)
{
    if (pos != m_vend->scenePos()) {
        prepareGeometryChange();
        m_vend->setScenePos(pos);
    }
}

QPointF active_projection::scenePos_beg() const { return m_beg->scenePos(); }
QPointF active_projection::scenePos_uend() const { return m_uend->scenePos(); }
QPointF active_projection::scenePos_vend() const { return m_vend->scenePos(); }


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
