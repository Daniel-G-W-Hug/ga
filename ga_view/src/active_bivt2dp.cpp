// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "active_bivt2dp.hpp"
#include "active_common.hpp"

active_bivt2dp::active_bivt2dp(Coordsys* cs, w_Coordsys* wcs, active_pt2d* beg,
                               active_pt2d* end, QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_beg{beg}, m_end{end}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);

    connect(wcs, &w_Coordsys::viewResized, m_beg, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_end, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, this, &active_bivt2dp::viewChanged);

    connect(this, &active_bivt2dp::viewMoved, m_beg, &active_pt2d::posChanged);
    connect(this, &active_bivt2dp::viewMoved, m_end, &active_pt2d::posChanged);

    connect(m_beg, &active_pt2d::pointMoved, this, &active_bivt2dp::viewChanged);
    connect(m_end, &active_pt2d::pointMoved, this, &active_bivt2dp::viewChanged);


    reset_item_data(); // update geometry (beg_pos, end_pos, and min- & max-values)
    setZValue(0);
}

void active_bivt2dp::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    qp->drawPath(arrowLine(beg_pos, end_pos));

    // from here on we want to draw with a small pen to get a pointy vector head
    QPen pen = qp->pen();
    pen.setWidth(1);
    qp->setBrush(pen.color());
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

QRectF active_bivt2dp::boundingRect() const
{
    return QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x), cs->y.au_to_w(max_y))),
                  mapFromScene(QPointF(cs->x.au_to_w(max_x), cs->y.au_to_w(min_y))));
}

QPainterPath active_bivt2dp::shape() const { return vectorShape(beg_pos, end_pos); }

void active_bivt2dp::setScenePos_beg(pt2d const& pos)
{
    if (pos != m_beg->scenePos()) {
        prepareGeometryChange();
        m_beg->setScenePos(pos);
    }
}

void active_bivt2dp::setScenePos_end(pt2d const& pos)
{
    if (pos != m_end->scenePos()) {
        prepareGeometryChange();
        m_end->setScenePos(pos);
    }
}

pt2d active_bivt2dp::scenePos_beg() const { return m_beg->scenePos(); }
pt2d active_bivt2dp::scenePos_end() const { return m_end->scenePos(); }


void active_bivt2dp::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_bivt2dp::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_bivt2dp::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_bivt2dp::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_bivt2dp::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_bivt2dp::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_bivt2dp: Qt::LeftButton.";
        m_mouse_l_pressed = true;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_bivt2dp: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }

    update();
    QGraphicsItem::mousePressEvent(event); // call default implementation
}

void active_bivt2dp::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_bivt2dp::mouseReleaseEvent.";

    // qDebug() << "active_bivt2dp::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_bivt2dp::scenePos_end():" << scenePos_end();

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

void active_bivt2dp::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_bivt2dp::mouseMoveEvent.";

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

void active_bivt2dp::viewChanged()
{
    // qDebug() << "active_bivt2dp: viewChanged() received.";

    // view changed by external influence, set to position in model
    reset_item_data(); // update geometry (beg_pos, end_pos, and min- & max-values)
    update();
}

void active_bivt2dp::reset_item_data()
{

    // update geometry, for initialization of if view changed

    try {

        // create the projective bivector using positions of both active points
        auto p = pt2dp(m_beg->scenePos().x, m_beg->scenePos().y, 1.0);
        auto q = pt2dp(m_end->scenePos().x, m_end->scenePos().y, 1.0);
        auto bvt = wdg(p, q);

        // Determine the angle of the projective line using atan2
        // For bivector (x,y,z), the line direction is (-y, x)
        // This gives us the full orientation in range [-π, π]
        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // qDebug() << "line_angle = " << rad2deg(line_angle);

        auto p1 = vec2dp(cs->x.min(), cs->y.max(), 1.0);
        auto p2 = vec2dp(cs->x.max(), cs->y.max(), 1.0);
        auto p3 = vec2dp(cs->x.min(), cs->y.min(), 1.0);
        auto p4 = vec2dp(cs->x.max(), cs->y.min(), 1.0);

        auto l1 = unitize(wdg(p1, p2)); // from min.x to max.x @ max.y
        auto l2 = unitize(wdg(p3, p1)); // from min.y to max.y @ min.x
        auto l3 = unitize(wdg(p3, p4)); // from min.x to max.x @ min.y
        auto l4 = unitize(wdg(p4, p2)); // from min.y to max.y @ max.x

        vec2dp p_from;
        vec2dp p_to;
        vec2dp p_from_24;
        vec2dp p_to_24;
        vec2dp p_from_13;
        vec2dp p_to_13;

        // Calculate intersection candidates for horizontal (left/right) boundaries
        // Skip when on ±y-axis: line_angle ≈ ±π (pos y-axis) or line_angle ≈ 0 (neg
        // y-axis)
        if (!(std::abs(std::abs(line_angle) - pi) < eps || std::abs(line_angle) < eps)) {
            p_from_24 = unitize(rwdg(l2, bvt)); // left boundary (lower x)
            p_to_24 = unitize(rwdg(l4, bvt));   // right boundary (higher x)
        }

        // Calculate intersection candidates for vertical (top/bottom) boundaries
        // Skip when on pos x-axis (line_angle ≈ -π/2) or neg x-axis (line_angle ≈ π/2)
        if (!(std::abs(line_angle - pi / 2.0) < eps ||
              std::abs(line_angle + pi / 2.0) < eps)) {
            p_from_13 = unitize(rwdg(l3, bvt)); // bottom boundary (lower y)
            p_to_13 = unitize(rwdg(l1, bvt));   // top boundary (higher y)
        }

        // Select final endpoints based on line_angle using atan2 ranges
        // pos. x-axis: line_angle ≈ -π/2
        if (std::abs(line_angle + pi / 2.0) < eps) {
            p_from = p_from_24;
            p_to = p_to_24;
        }
        // first quadrant: line_angle ∈ (-π, -π/2)
        else if (line_angle < -pi / 2.0 - eps && line_angle > -pi + eps) {
            if (p_to_13.x < cs->x.max()) {
                p_to = p_to_13;
            }
            else {
                p_to = p_to_24;
            }
            if (p_from_13.x > cs->x.min()) {
                p_from = p_from_13;
            }
            else {
                p_from = p_from_24;
            }
        }
        // pos. y-axis: line_angle ≈ ±π
        else if (std::abs(std::abs(line_angle) - pi) < eps) {
            p_from = p_from_13;
            p_to = p_to_13;
        }
        // second quadrant: line_angle ∈ (π/2, π)
        else if (line_angle > pi / 2.0 + eps && line_angle < pi - eps) {
            if (p_from_24.y < cs->y.max()) {
                p_to = p_from_24;
            }
            else {
                p_to = p_to_13;
            }
            if (p_to_24.y > cs->y.min()) {
                p_from = p_to_24;
            }
            else {
                p_from = p_from_13;
            }
        }
        // neg. x-axis: line_angle ≈ π/2
        else if (std::abs(line_angle - pi / 2.0) < eps) {
            p_from = p_to_24;
            p_to = p_from_24;
        }
        // third quadrant: line_angle ∈ (0, π/2)
        else if (line_angle > eps && line_angle < pi / 2.0 - eps) {
            if (p_from_13.x > cs->x.min()) {
                p_to = p_from_13;
            }
            else {
                p_to = p_from_24;
            }
            if (p_to_13.x < cs->x.max()) {
                p_from = p_to_13;
            }
            else {
                p_from = p_to_24;
            }
        }
        // neg. y-axis: line_angle ≈ 0
        else if (std::abs(line_angle) < eps) {
            p_from = p_to_13;
            p_to = p_from_13;
        }
        // fourth quadrant: line_angle ∈ (-π/2, 0)
        else if (line_angle > -pi / 2.0 + eps && line_angle < -eps) {
            if (p_to_24.y > cs->y.min()) {
                p_to = p_to_24;
            }
            else {
                p_to = p_from_13;
            }
            if (p_from_24.y < cs->y.max()) {
                p_from = p_from_24;
            }
            else {
                p_from = p_to_13;
            }
        }

        // fmt::println("idx={},       B.x={:.4g}, B.y={:.4g}, B.z={:.4g}", idx,
        //              bvt.x, bvt.y, bvt.z);
        // fmt::println("phi_x={:.3g}, phi_y={:.3g}, att={:.4g}", rad2deg(phi_x),
        //              rad2deg(phi_y), att(bvt));
        // fmt::println("p_from_24 = {:.2g}, p_to_24 = {:.2g}", p_from_24, p_to_24);
        // fmt::println("p_from_13 = {:.2g}, p_to_13 = {:.2g}", p_from_13, p_to_13);
        // fmt::println("p_from    = {:.2g}, p_to    = {:.2g}", p_from, p_to);

        // create the lines and map them to the coordinate system
        beg_pos = mapFromScene(QPointF(cs->x.au_to_w(p_from.x), cs->y.au_to_w(p_from.y)));
        end_pos = mapFromScene(QPointF(cs->x.au_to_w(p_to.x), cs->y.au_to_w(p_to.y)));

        // find min and max values for bounding rectangle:
        // QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x),cs->y.au_to_w(max_y))),
        //        mapFromScene(QPointF(cs->x.au_to_w(max_x),cs->y.au_to_w(min_y))))
        min_x = std::min(min_x, p_from.x);
        max_x = std::max(max_x, p_from.x);
        min_y = std::min(min_y, p_from.y);
        max_y = std::max(max_y, p_from.y);
        min_x = std::min(min_x, p_to.x);
        max_x = std::max(max_x, p_to.x);
        min_y = std::min(min_y, p_to.y);
        max_y = std::max(max_y, p_to.y);

        // qDebug() << "min_x: " << min_x;
        // qDebug() << "max_x: " << max_x;
        // qDebug() << "min_y: " << min_y;
        // qDebug() << "max_y: " << max_y;

        // fmt::println("");
    }
    catch (std::exception const& e) {
        qDebug() << e.what();
    }
}
