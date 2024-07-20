//
// author: Daniel Hug, 2024
//

#include "active_reflection.hpp"
#include "active_common.hpp"

#include "ga/ga.hpp"

#include <algorithm> // std::max

active_reflection::active_reflection(Coordsys* cs, w_Coordsys* wcs, active_pt* n1end,
                                     active_pt* n2end, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_n1end{n1end}, m_n2end{n2end}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    // setZValue(18);

    connect(wcs, &w_Coordsys::viewResized, m_n1end, &active_pt::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_n2end, &active_pt::viewChanged);

    connect(this, &active_reflection::viewMoved, m_n1end, &active_pt::posChanged);
    connect(this, &active_reflection::viewMoved, m_n2end, &active_pt::posChanged);
}

void active_reflection::paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                              QWidget* widget)
{

    // clipping area is active area of coordsys
    painter->setClipRect(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax()));

    // draw in item coordinate system
    painter->save();

    // draw first normal vector and corresponding hyperplane

    QPointF beg_pos = QPointF(cs->x.a_to_w(0.0), cs->y.a_to_w(0.0));

    QPointF end_n1pos_nrm_scene = scenePos_n1end() / nrm(scenePos_n1end());
    QPointF end_n1pos =
        QPointF(cs->x.a_to_w(scenePos_n1end().x()), cs->y.a_to_w(scenePos_n1end().y()));
    QPointF end_n1pos_nrm = QPointF(cs->x.a_to_w(end_n1pos_nrm_scene.x()),
                                    cs->y.a_to_w(end_n1pos_nrm_scene.y()));

    painter->setPen(QPen(QBrush(col_lgreen), 2, Qt::SolidLine));
    painter->setBrush(col_lgreen);
    painter->drawPath(arrowLine(beg_pos, end_n1pos_nrm));
    painter->setPen(Qt::gray);
    painter->drawPath(arrowLine(end_n1pos_nrm, end_n1pos));
    painter->setPen(QPen(QBrush(col_lgreen), 2, Qt::SolidLine));
    // from here on we want to draw with a small pen to get a pointy vector head
    QPen pen = painter->pen();
    pen.setWidth(1);
    painter->drawPath(arrowHead(beg_pos, end_n1pos_nrm));
    pen.setWidth(2);
    painter->drawPath(planeLine(beg_pos, end_n1pos,
                                std::max(cs->x.widget_size(), cs->y.widget_size())));


    // draw second normal vector and corresponding hyperplane

    QPointF end_n2pos_nrm_scene = scenePos_n2end() / nrm(scenePos_n2end());
    QPointF end_n2pos =
        QPointF(cs->x.a_to_w(scenePos_n2end().x()), cs->y.a_to_w(scenePos_n2end().y()));
    QPointF end_n2pos_nrm = QPointF(cs->x.a_to_w(end_n2pos_nrm_scene.x()),
                                    cs->y.a_to_w(end_n2pos_nrm_scene.y()));

    painter->setPen(QPen(QBrush(col_lred), 2, Qt::SolidLine));
    painter->setBrush(col_lred);
    painter->drawPath(arrowLine(beg_pos, end_n2pos_nrm));
    painter->setPen(Qt::gray);
    painter->drawPath(arrowLine(end_n2pos_nrm, end_n2pos));
    painter->setPen(QPen(QBrush(col_lred), 2, Qt::SolidLine));
    painter->setBrush(col_lred);
    // from here on we want to draw with a small pen to get a pointy vector head
    pen.setWidth(1);
    painter->drawPath(arrowHead(beg_pos, end_n2pos_nrm));
    pen.setWidth(2);
    painter->drawPath(planeLine(beg_pos, end_n2pos,
                                std::max(cs->x.widget_size(), cs->y.widget_size())));


    // let's do a little ga
    using hd::ga::vec2d;
    std::vector<vec2d> v = {{0.5, 0.5}, {2.0, 0.5}, {2.5, 1.0}};
    vec2d n1{end_n1pos_nrm_scene.x(), end_n1pos_nrm_scene.y()};
    vec2d n2{end_n2pos_nrm_scene.x(), end_n2pos_nrm_scene.y()};
    std::vector<vec2d> vr, vrr;

    using hd::ga::reflect_on_hyp;
    for (auto const& e : v) {
        vec2d v_tmp = reflect_on_hyp(e, n1);
        vr.push_back(v_tmp);
        vrr.push_back(reflect_on_hyp(v_tmp, n2));
    }

    QPainterPath org, refl1, refl2;
    org.moveTo(QPointF(cs->x.a_to_w(v[0].x), cs->y.a_to_w(v[0].y)));
    refl1.moveTo(QPointF(cs->x.a_to_w(vr[0].x), cs->y.a_to_w(vr[0].y)));
    refl2.moveTo(QPointF(cs->x.a_to_w(vrr[0].x), cs->y.a_to_w(vrr[0].y)));
    for (size_t i = 1; i < v.size(); ++i) {
        org.lineTo(QPointF(cs->x.a_to_w(v[i].x), cs->y.a_to_w(v[i].y)));
        refl1.lineTo(QPointF(cs->x.a_to_w(vr[i].x), cs->y.a_to_w(vr[i].y)));
        refl2.lineTo(QPointF(cs->x.a_to_w(vrr[i].x), cs->y.a_to_w(vrr[i].y)));
    }
    org.closeSubpath();
    refl1.closeSubpath();
    refl2.closeSubpath();

    painter->setPen(QPen(QBrush(col_lblue), 1, Qt::SolidLine)); // pressed: red
    painter->setBrush(col_lblue);
    painter->drawPath(org);

    painter->setPen(QPen(QBrush(col_lgreen), 1, Qt::SolidLine)); // pressed: red
    painter->setBrush(col_lgreen);
    painter->drawPath(refl1);

    painter->setPen(QPen(QBrush(col_lred), 1, Qt::SolidLine)); // pressed: red
    painter->setBrush(col_lred);
    painter->drawPath(refl2);

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
    QPointF beg_pos = QPointF(cs->x.a_to_w(0.0), cs->y.a_to_w(0.0));
    QPointF end_n1pos =
        QPointF(cs->x.a_to_w(scenePos_n1end().x()), cs->y.a_to_w(scenePos_n1end().y()));
    QPointF end_n2pos =
        QPointF(cs->x.a_to_w(scenePos_n2end().x()), cs->y.a_to_w(scenePos_n2end().y()));

    QPainterPath path = vectorShape(beg_pos, end_n1pos);
    path += vectorShape(beg_pos, end_n2pos);

    return path;
}

void active_reflection::setScenePos_n1end(QPointF const& pos)
{
    if (pos != m_n1end->scenePos()) {
        prepareGeometryChange();
        m_n1end->setScenePos(pos);
    }
}

void active_reflection::setScenePos_n2end(QPointF const& pos)
{
    if (pos != m_n2end->scenePos()) {
        prepareGeometryChange();
        m_n2end->setScenePos(pos);
    }
}

QPointF active_reflection::scenePos_n1end() const { return m_n1end->scenePos(); }

QPointF active_reflection::scenePos_n2end() const { return m_n2end->scenePos(); }


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
    QGraphicsItem::mousePressEvent(event); // call default implementation
}

void active_reflection::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_reflection::mouseReleaseEvent.";

    // qDebug() << "active_reflection::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_reflection::scenePos_n1end():" << scenePos_n1end();

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

void active_reflection::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_reflection::mouseMoveEvent.";

    // if (m_mouse_l_pressed) {

    //     QPointF delta = event->scenePos() - event->lastScenePos();

    //     if (delta != QPointF(0, 0)) {

    //         // qDebug() << "scenePos():" << event->scenePos();
    //         // qDebug() << "lastScenePos():" << event->lastScenePos();
    //         // qDebug() << "delta:" << delta;

    //         m_n1end->moveBy(delta.x(), delta.y());
    //         m_n2end->moveBy(delta.x(), delta.y());

    //         emit viewMoved();
    //     }
    // }
}
