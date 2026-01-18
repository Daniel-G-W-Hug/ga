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

        // Two-angle approach: measure angle from x-axis and (-y)-axis Using (-y)-axis
        // gives consistent angle ranges that work correctly with the angle() function
        // which returns values in [0, π]
        auto const x_axis = bivec2dp{1, 0, 0};
        auto const neg_y_axis = bivec2dp{0, -1, 0};

        auto const phi_x = angle(x_axis, bvt);     // angle from x-axis [0, π]
        auto const phi_y = angle(neg_y_axis, bvt); // angle from (-y)-axis [0, π]

        // Define corner points of the coordinate system boundary
        auto p1 = vec2dp(cs->x.min(), cs->y.max(), 1.0);
        auto p2 = vec2dp(cs->x.max(), cs->y.max(), 1.0);
        auto p3 = vec2dp(cs->x.min(), cs->y.min(), 1.0);
        auto p4 = vec2dp(cs->x.max(), cs->y.min(), 1.0);

        // Define boundary lines (edges of the coordinate system)
        auto l1 = unitize(wdg(p1, p2)); // top:    from min.x to max.x @ max.y
        auto l2 = unitize(wdg(p3, p1)); // left:   from min.y to max.y @ min.x
        auto l3 = unitize(wdg(p3, p4)); // bottom: from min.x to max.x @ min.y
        auto l4 = unitize(wdg(p4, p2)); // right:  from min.y to max.y @ max.x

        vec2dp p_from;
        vec2dp p_to;
        vec2dp p_from_24; // intersection with left boundary (l2)
        vec2dp p_to_24;   // intersection with right boundary (l4)
        vec2dp p_from_13; // intersection with bottom boundary (l3)
        vec2dp p_to_13;   // intersection with top boundary (l1)

        // Calculate intersection with left/right boundaries (l2, l4)
        // Skip when line is near horizontal (phi_x ≈ π/2), parallel to these boundaries
        if (std::abs(phi_x - pi / 2.0) > eps) {
            p_from_24 = unitize(rwdg(l2, bvt)); // left boundary (lower x)
            p_to_24 = unitize(rwdg(l4, bvt));   // right boundary (higher x)
        }

        // Calculate intersection with top/bottom boundaries (l1, l3)
        // Skip when line is near vertical (phi_y ≈ π/2), parallel to these boundaries
        if (std::abs(phi_y - pi / 2.0) > eps) {
            p_from_13 = unitize(rwdg(l3, bvt)); // bottom boundary (lower y)
            p_to_13 = unitize(rwdg(l1, bvt));   // top boundary (higher y)
        }

        // Select final endpoints based on phi_x and phi_y angles
        // With (-y)-axis reference:
        // phi_x ≈ 0:   line parallel to +x axis
        // phi_x ≈ π:   line parallel to -x axis
        // phi_y ≈ 0:   line parallel to -y axis (pointing down)
        // phi_y ≈ π:   line parallel to +y axis (pointing up)
        //
        // Quadrant angle ranges:
        // Q1 (upper-right): phi_x ∈ (0, π/2), phi_y ∈ (π/2, π)
        // Q2 (upper-left):  phi_x ∈ (π/2, π), phi_y ∈ (π/2, π)
        // Q3 (lower-left):  phi_x ∈ (π/2, π), phi_y ∈ (0, π/2)
        // Q4 (lower-right): phi_x ∈ (0, π/2), phi_y ∈ (0, π/2)

        // positive x-axis: phi_x ≈ 0, phi_y ≈ π/2
        if ((std::abs(phi_x) < eps) && (std::abs(phi_y - pi / 2.0) < eps)) {
            p_from = p_from_24;
            p_to = p_to_24;
        }

        // first quadrant (upper-right): phi_x ∈ (0, π/2), phi_y ∈ (π/2, π)
        else if (phi_x > eps && phi_x < (pi / 2.0 - eps) && phi_y > (pi / 2.0 + eps) &&
                 phi_y < (pi - eps)) {
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

        // positive y-axis: phi_x ≈ π/2, phi_y ≈ π
        else if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y - pi) < eps)) {
            p_from = p_from_13;
            p_to = p_to_13;
        }

        // second quadrant (upper-left): phi_x ∈ (π/2, π), phi_y ∈ (π/2, π)
        else if (phi_x > (pi / 2.0 + eps) && phi_x < (pi - eps) &&
                 phi_y > (pi / 2.0 + eps) && phi_y < (pi - eps)) {
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

        // negative x-axis: phi_x ≈ π, phi_y ≈ π/2
        else if ((std::abs(phi_x - pi) < eps) && (std::abs(phi_y - pi / 2.0) < eps)) {
            p_from = p_to_24;
            p_to = p_from_24;
        }

        // third quadrant (lower-left): phi_x ∈ (π/2, π), phi_y ∈ (0, π/2)
        else if (phi_x > (pi / 2.0 + eps) && phi_x < (pi - eps) && phi_y > eps &&
                 phi_y < (pi / 2.0 - eps)) {
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

        // negative y-axis: phi_x ≈ π/2, phi_y ≈ 0
        else if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y) < eps)) {
            p_from = p_to_13;
            p_to = p_from_13;
        }

        // fourth quadrant (lower-right): phi_x ∈ (0, π/2), phi_y ∈ (0, π/2)
        else if (phi_x > eps && phi_x < (pi / 2.0 - eps) && phi_y > eps &&
                 phi_y < (pi / 2.0 - eps)) {
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
