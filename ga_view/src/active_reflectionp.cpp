// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_reflectionp.hpp"
#include "active_common.hpp"

#include "ga/ga_ega.hpp"

#include <algorithm> // std::max
#include <array>

active_reflectionp::active_reflectionp(Coordsys* cs, w_Coordsys* wcs, active_pt2d* p1,
                                       active_pt2d* p2, active_pt2d* p3, active_pt2d* p4,
                                       QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_p1{p1}, m_p2{p2}, m_p3{p3}, m_p4{p4}
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
             QGraphicsItem::ItemIsFocusable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(true);


    connect(wcs, &w_Coordsys::viewResized, m_p1, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_p2, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_p3, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, m_p4, &active_pt2d::viewChanged);
    connect(wcs, &w_Coordsys::viewResized, this, &active_reflectionp::viewChanged);

    connect(wcs, &w_Coordsys::moveModeChanged, this,
            &active_reflectionp::moveModeChanged);

    connect(this, &active_reflectionp::viewMoved, m_p1, &active_pt2d::posChanged);
    connect(this, &active_reflectionp::viewMoved, m_p2, &active_pt2d::posChanged);
    connect(this, &active_reflectionp::viewMoved, m_p3, &active_pt2d::posChanged);
    connect(this, &active_reflectionp::viewMoved, m_p4, &active_pt2d::posChanged);

    connect(this, &active_reflectionp::pointsMoved, m_p1, &active_pt2d::viewChanged);
    connect(this, &active_reflectionp::pointsMoved, m_p2, &active_pt2d::viewChanged);
    connect(this, &active_reflectionp::pointsMoved, m_p3, &active_pt2d::viewChanged);
    connect(this, &active_reflectionp::pointsMoved, m_p4, &active_pt2d::viewChanged);

    connect(m_p1, &active_pt2d::pointMoved, this, &active_reflectionp::viewChanged);
    connect(m_p2, &active_pt2d::pointMoved, this, &active_reflectionp::viewChanged);
    connect(m_p3, &active_pt2d::pointMoved, this, &active_reflectionp::viewChanged);
    connect(m_p4, &active_pt2d::pointMoved, this, &active_reflectionp::viewChanged);

    reset_item_data(); // update geometry (beg_pos, end_pos, and min- & max-values)
    setZValue(0);
}

void active_reflectionp::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    // draw first hyperplane (=line)

    qp->drawPath(arrowLine(bvt1_beg_pos, bvt1_end_pos));
    // from here on we want to draw with a small pen to get a pointy vector head
    QPen pen = qp->pen();
    pen.setWidth(1);
    qp->setBrush(pen.color());
    qp->drawPath(arrowHead(bvt1_beg_pos, bvt1_end_pos));


    // draw second hyperplane (=line)

    qp->drawPath(arrowLine(bvt2_beg_pos, bvt2_end_pos));
    // from here on we want to draw with a small pen to get a pointy vector head
    pen = qp->pen();
    pen.setWidth(1);
    qp->setBrush(pen.color());
    qp->drawPath(arrowHead(bvt2_beg_pos, bvt2_end_pos));


    // let's do a little ga
    std::array<vec2dp, 3> const v = {vec2dp{0.5, 0.5, 1.0}, vec2dp{2.0, 0.5, 1.0},
                                     vec2dp{2.5, 1.0, 1.0}};

    auto p1 = vec2dp(m_p1->scenePos().x, m_p1->scenePos().y, 1.0);
    auto q1 = vec2dp(m_p2->scenePos().x, m_p2->scenePos().y, 1.0);
    auto bvt1 = unitize(wdg(p1, q1));
    // auto bvt1 = wdg(p1, q1);

    auto p2 = vec2dp(m_p3->scenePos().x, m_p3->scenePos().y, 1.0);
    auto q2 = vec2dp(m_p4->scenePos().x, m_p4->scenePos().y, 1.0);
    auto bvt2 = unitize(wdg(p2, q2));

    std::vector<vec2dp> vr, vrr;
    vr.reserve(v.size());
    vrr.reserve(v.size());

    for (auto const& e : v) {
        vec2dp v_tmp = reflect_on(e, bvt1);
        vr.push_back(v_tmp);
        vrr.emplace_back(reflect_on(v_tmp, bvt2));
    }

    QPainterPath org, refl1, refl2;
    org.moveTo(QPointF(cs->x.au_to_w(v[0].x), cs->y.au_to_w(v[0].y)));
    refl1.moveTo(QPointF(cs->x.au_to_w(vr[0].x), cs->y.au_to_w(vr[0].y)));
    refl2.moveTo(QPointF(cs->x.au_to_w(vrr[0].x), cs->y.au_to_w(vrr[0].y)));
    for (size_t i = 1; i < v.size(); ++i) {
        org.lineTo(QPointF(cs->x.au_to_w(v[i].x), cs->y.au_to_w(v[i].y)));
        refl1.lineTo(QPointF(cs->x.au_to_w(vr[i].x), cs->y.au_to_w(vr[i].y)));
        refl2.lineTo(QPointF(cs->x.au_to_w(vrr[i].x), cs->y.au_to_w(vrr[i].y)));
    }
    org.closeSubpath();
    refl1.closeSubpath();
    refl2.closeSubpath();

    qp->setPen(QPen(QBrush(col_lblue), 1, Qt::SolidLine)); // pressed: red
    qp->setBrush(col_lblue);
    qp->drawPath(org);

    qp->setPen(QPen(QBrush(col_lgreen), 1, Qt::SolidLine)); // pressed: red
    qp->setBrush(col_lgreen);
    qp->drawPath(refl1);

    qp->setPen(QPen(QBrush(col_lred), 1, Qt::SolidLine)); // pressed: red
    qp->setBrush(col_lred);
    qp->drawPath(refl2);

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    // draw shape (optional for testing)
    qp->setPen(col_yel);
    qp->setBrush(col_yel);
    qp->drawPath(shape());

    qp->restore();
}

QRectF active_reflectionp::boundingRect() const
{
    return QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x), cs->y.au_to_w(max_y))),
                  mapFromScene(QPointF(cs->x.au_to_w(max_x), cs->y.au_to_w(min_y))))
        .normalized();
}

QPainterPath active_reflectionp::shape() const
{
    QPainterPath path = vectorShape(bvt1_beg_pos, bvt1_end_pos);
    path += vectorShape(bvt2_beg_pos, bvt2_end_pos);

    return path;
}

void active_reflectionp::setScenePos_p1(pt2d const& pos)
{
    if (pos != m_p1->scenePos()) {
        prepareGeometryChange();
        m_p1->setScenePos(pos);
        reset_item_data();
    }
}

void active_reflectionp::setScenePos_p2(pt2d const& pos)
{
    if (pos != m_p2->scenePos()) {
        prepareGeometryChange();
        m_p2->setScenePos(pos);
        reset_item_data();
    }
}

void active_reflectionp::setScenePos_p3(pt2d const& pos)
{
    if (pos != m_p3->scenePos()) {
        prepareGeometryChange();
        m_p3->setScenePos(pos);
        reset_item_data();
    }
}

void active_reflectionp::setScenePos_p4(pt2d const& pos)
{
    if (pos != m_p4->scenePos()) {
        prepareGeometryChange();
        m_p4->setScenePos(pos);
        reset_item_data();
    }
}

void active_reflectionp::setScenePos_p1_wo_update(pt2d const& pos)
{
    if (pos != m_p1->scenePos()) {
        m_p1->setScenePos_wo_update(pos);
    }
}

void active_reflectionp::setScenePos_p2_wo_update(pt2d const& pos)
{
    if (pos != m_p2->scenePos()) {
        m_p2->setScenePos_wo_update(pos);
    }
}

void active_reflectionp::setScenePos_p3_wo_update(pt2d const& pos)
{
    if (pos != m_p3->scenePos()) {
        m_p3->setScenePos_wo_update(pos);
    }
}

void active_reflectionp::setScenePos_p4_wo_update(pt2d const& pos)
{
    if (pos != m_p4->scenePos()) {
        m_p4->setScenePos_wo_update(pos);
    }
}

pt2d active_reflectionp::scenePos_p1() const { return m_p1->scenePos(); }

pt2d active_reflectionp::scenePos_p2() const { return m_p2->scenePos(); }

pt2d active_reflectionp::scenePos_p3() const { return m_p3->scenePos(); }

pt2d active_reflectionp::scenePos_p4() const { return m_p4->scenePos(); }


void active_reflectionp::hoverEnterEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_reflectionp::hoverEnterEvent.";
    m_mouse_hover = true;
    update();
}

void active_reflectionp::hoverLeaveEvent(QGraphicsSceneHoverEvent* event)
{
    Q_UNUSED(event)

    // qDebug() << "active_reflectionp::hoverLeaveEvent.";
    m_mouse_hover = false;
    update();
}

void active_reflectionp::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_reflectionp::mousePressEvent.";

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_reflectionp: Qt::LeftButton.";
        m_mouse_l_pressed = true;

        if (m_move_mode == move_mode::rotate_both_lines) {
            // calculate and store the current intersection point
            auto p1 = vec2dp(m_p1->scenePos().x, m_p1->scenePos().y, 1.0);
            auto q1 = vec2dp(m_p2->scenePos().x, m_p2->scenePos().y, 1.0);
            auto bvt1 = unitize(wdg(p1, q1));

            auto p2 = vec2dp(m_p3->scenePos().x, m_p3->scenePos().y, 1.0);
            auto q2 = vec2dp(m_p4->scenePos().x, m_p4->scenePos().y, 1.0);
            auto bvt2 = unitize(wdg(p2, q2));

            // set private variable tp to current intersection point
            // "turning point" tp is used in mouseMoveEvent
            tp = rwdg(bvt1, bvt2); // must not be unitized, can be at infinity
        }
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_reflectionp: Qt::RightButton.";
        m_mouse_r_pressed = true;
    }
    update();
}

void active_reflectionp::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    // qDebug() << "active_reflectionp::mouseReleaseEvent.";
    // qDebug() << "active_reflectionp::scenePos_beg():" << scenePos_beg();
    // qDebug() << "active_reflectionp::scenePos_p1():" << scenePos_p1();

    if (event->button() == Qt::LeftButton) {
        // qDebug() << "active_pt2d: Qt::LeftButton.";
        m_mouse_l_pressed = false;
    }
    if (event->button() == Qt::RightButton) {
        // qDebug() << "active_pt2d: Qt::RightButton.";
        m_mouse_r_pressed = false;
    }
    update();
}

void active_reflectionp::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // Q_UNUSED(event)

    // qDebug() << "active_reflectionp::mouseMoveEvent.";

    if (m_mouse_l_pressed) {

        QPointF delta = event->scenePos() - event->lastScenePos();

        if (delta != QPointF(0, 0)) {

            // qDebug() << "scenePos():" << event->scenePos();
            // qDebug() << "lastScenePos():" << event->lastScenePos();
            // qDebug() << "delta:" << delta;d

            switch (m_move_mode) {
                case move_mode::shift_both_lines:
                    m_p1->moveBy(delta.x(), delta.y());
                    m_p2->moveBy(delta.x(), delta.y());
                    m_p3->moveBy(delta.x(), delta.y());
                    m_p4->moveBy(delta.x(), delta.y());
                    emit viewMoved();
                    break;

                case move_mode::shift_line12:
                    m_p1->moveBy(delta.x(), delta.y());
                    m_p2->moveBy(delta.x(), delta.y());
                    emit viewMoved();
                    break;

                case move_mode::shift_line34:
                    m_p3->moveBy(delta.x(), delta.y());
                    m_p4->moveBy(delta.x(), delta.y());
                    emit viewMoved();
                    break;

                case move_mode::rotate_both_lines: {

                    // get current turning angle
                    auto scnPos = vec2d(cs->x.w_to_au(event->scenePos().x()),
                                        cs->y.w_to_au(event->scenePos().y()));
                    auto lscnPos = vec2d(cs->x.w_to_au(event->lastScenePos().x()),
                                         cs->y.w_to_au(event->lastScenePos().y()));
                    auto cur_ang = angle(scnPos, lscnPos);
                    // fmt::println("ang = {}", rad2deg(cur_ang));


                    auto cur_p1p = vec2dp(m_p1->scenePos().x, m_p1->scenePos().y, 1.0);
                    auto cur_p2p = vec2dp(m_p2->scenePos().x, m_p2->scenePos().y, 1.0);
                    auto cur_p3p = vec2dp(m_p3->scenePos().x, m_p3->scenePos().y, 1.0);
                    auto cur_p4p = vec2dp(m_p4->scenePos().x, m_p4->scenePos().y, 1.0);

                    // move all 4 projective points around the turning point pt
                    // auto mot = motor(cur_p1p, cur_ang);
                    auto mot = motor(tp, cur_ang);
                    // fmt::println("tp = {}", tp);
                    // fmt::println("mot = {}", mot);

                    auto new_p1p = unitize(move2dp(cur_p1p, mot));
                    auto new_p2p = unitize(move2dp(cur_p2p, mot));
                    auto new_p3p = unitize(move2dp(cur_p3p, mot));
                    auto new_p4p = unitize(move2dp(cur_p4p, mot));

                    setScenePos_p1_wo_update(vec2d(new_p1p.x, new_p1p.y));
                    setScenePos_p2_wo_update(vec2d(new_p2p.x, new_p2p.y));
                    setScenePos_p3_wo_update(vec2d(new_p3p.x, new_p3p.y));
                    setScenePos_p4_wo_update(vec2d(new_p4p.x, new_p4p.y));
                }
                    emit pointsMoved();
                    break;
                default:
                    std::unreachable();
            }
            viewChanged();
        }
    }
}

void active_reflectionp::viewChanged()
{
    // qDebug() << "active_bivt2dp: viewChanged() received.";

    // view changed by external influence, set to position in model
    reset_item_data(); // update geometry (beg_pos, end_pos, and min- & max-values)
    update();
}

void active_reflectionp::moveModeChanged(move_mode mode)
{

    if (mode != m_move_mode) { // only adjust, if actually changed
        m_move_mode = mode;

        if (m_move_mode == move_mode::rotate_both_lines) {
            // calculate and store the current intersection point
            auto p1 = vec2dp(m_p1->scenePos().x, m_p1->scenePos().y, 1.0);
            auto q1 = vec2dp(m_p2->scenePos().x, m_p2->scenePos().y, 1.0);
            auto bvt1 = unitize(wdg(p1, q1));

            auto p2 = vec2dp(m_p3->scenePos().x, m_p3->scenePos().y, 1.0);
            auto q2 = vec2dp(m_p4->scenePos().x, m_p4->scenePos().y, 1.0);
            auto bvt2 = unitize(wdg(p2, q2));

            // set private variable tp to current intersection point
            // "turning point" tp is used in mouseMoveEvent
            tp = rwdg(bvt1, bvt2); // must not be unitized, can be at infinity
        }
        update();
    }
}

void active_reflectionp::reset_item_data()
{

    struct res_vec2dp {
        vec2dp p_from;
        vec2dp p_to;
    };

    // get the actual start and end points of the line in the coordsystem
    auto reset_bivecp = [this](pt2dp p, pt2dp q) -> res_vec2dp {
        // determine the angle of the projective line
        auto const x_axis = bivec2dp{0, 1, 0};
        auto const y_axis = bivec2dp{1, 0, 0}; // really this is -y_axis_2dp

        auto bvt = wdg(p, q);
        auto const phi_x = angle(x_axis, bvt);
        auto const phi_y = angle(y_axis, bvt);

        // qDebug() << "phi_x = " << rad2deg(phi_x);
        // qDebug() << "phi_y = " << rad2deg(phi_y);

        auto p1 = vec2dp(cs->x.min(), cs->y.max(), 1.0);
        auto p2 = vec2dp(cs->x.max(), cs->y.max(), 1.0);
        auto p3 = vec2dp(cs->x.min(), cs->y.min(), 1.0);
        auto p4 = vec2dp(cs->x.max(), cs->y.min(), 1.0);

        auto l1 = unitize(wdg(p1, p2)); // from min.x to max.x @ max.y
        auto l2 = unitize(wdg(p3, p1)); // from min.y to max.y @ min.x
        auto l3 = unitize(wdg(p3, p4)); // from min.x to max.x @ min.y
        auto l4 = unitize(wdg(p4, p2)); // from min.y to max.y @ max.x

        res_vec2dp res;
        vec2dp p_from_24;
        vec2dp p_to_24;
        vec2dp p_from_13;
        vec2dp p_to_13;

        if (std::abs(phi_x - pi / 2.0) > eps) {
            // if not on pos. or neg. y-axis
            p_from_24 = unitize(rwdg(l2, bvt)); // from lower x values
            p_to_24 = unitize(rwdg(l4, bvt));   // to higher x values
        }

        if (std::abs(phi_y - pi / 2.0) > eps) {
            // if not on pos.or neg.x - axis
            p_from_13 = unitize(rwdg(l3, bvt)); // from lower y values
            p_to_13 = unitize(rwdg(l1, bvt));   // to higher y values
        }

        // pos. x-axis
        if ((std::abs(phi_x) < eps) && (std::abs(phi_y - pi / 2.0) < eps)) {

            // qDebug() << "On positive x-axis.";

            res.p_from = p_from_24;
            res.p_to = p_to_24;
        }

        // first quadrant
        if (phi_x > eps && (phi_x < (pi / 2. - eps)) && (phi_y > (pi / 2. + eps)) &&
            (phi_y < (pi - eps))) {

            // qDebug() << "In first quadrant.";

            if (p_to_13.x < cs->x.max()) {
                res.p_to = p_to_13;
            }
            else {
                res.p_to = p_to_24;
            }

            if (p_from_13.x > cs->x.min()) {
                res.p_from = p_from_13;
            }
            else {
                res.p_from = p_from_24;
            }
        }

        // pos. y-axis
        if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y - pi) < eps)) {

            // qDebug() << "On positive y-axis.";

            res.p_from = p_from_13;
            res.p_to = p_to_13;
        }

        // second quadrant
        if ((phi_x > (pi / 2.0 + eps)) && (phi_x < (pi - eps)) && (phi_y < (pi - eps)) &&
            (phi_y > (pi / 2.0 + eps))) {

            // qDebug() << "In second quadrant.";

            if (p_from_24.y < cs->y.max()) {
                res.p_to = p_from_24;
            }
            else {
                res.p_to = p_to_13;
            }

            if (p_to_24.y > cs->y.min()) {
                res.p_from = p_to_24;
            }
            else {
                res.p_from = p_from_13;
            }
        }

        // neg. x-axis
        if ((std::abs(phi_x - pi) < eps) && (std::abs(phi_y - pi / 2.0)) < eps) {

            // qDebug() << "On negative x-axis.";

            res.p_from = p_to_24;
            res.p_to = p_from_24;
        }

        // third quadrant
        if ((phi_x < pi - eps) && (phi_x > pi / 2.0 + eps) && (phi_y < pi / 2.0 - eps) &&
            (phi_y > eps)) {

            // qDebug() << "in third quadrant.";

            if (p_from_13.x > cs->x.min()) {
                res.p_to = p_from_13;
            }
            else {
                res.p_to = p_from_24;
            }

            if (p_to_13.x < cs->x.max()) {
                res.p_from = p_to_13;
            }
            else {
                res.p_from = p_to_24;
            }
        }

        // neg. y-axis
        if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y) < eps)) {

            // qDebug() << "On negative y-axis.";

            res.p_from = p_to_13;
            res.p_to = p_from_13;
        }

        // fourth quadrant
        if ((phi_x < pi / 2.0 - eps) && (phi_x > eps) && (phi_y > eps) &&
            (phi_y < pi / 2. - eps)) {

            // qDebug() << "In fourth quadrant.";

            if (p_to_24.y > cs->y.min()) {
                res.p_to = p_to_24;
            }
            else {
                res.p_to = p_from_13;
            }

            if (p_from_24.y < cs->y.max()) {
                res.p_from = p_from_24;
            }
            else {
                res.p_from = p_to_13;
            }
        }

        return res;
    };

    // update geometry, for initialization of if view changed
    try {

        // create the first projective bivector using positions of both active points
        auto p = pt2dp(m_p1->scenePos().x, m_p1->scenePos().y, 1.0);
        auto q = pt2dp(m_p2->scenePos().x, m_p2->scenePos().y, 1.0);

        auto res = reset_bivecp(p, q);

        // create the lines and map them to the coordinate system
        bvt1_beg_pos = mapFromScene(
            QPointF(cs->x.au_to_w(res.p_from.x), cs->y.au_to_w(res.p_from.y)));
        bvt1_end_pos =
            mapFromScene(QPointF(cs->x.au_to_w(res.p_to.x), cs->y.au_to_w(res.p_to.y)));


        // find min and max values for bounding rectangle:
        auto min_x1 = std::min(std::numeric_limits<double>::max(), res.p_from.x);
        auto max_x1 = std::max(std::numeric_limits<double>::lowest(), res.p_from.x);
        auto min_y1 = std::min(std::numeric_limits<double>::max(), res.p_from.y);
        auto max_y1 = std::max(std::numeric_limits<double>::lowest(), res.p_from.y);
        min_x1 = std::min(min_x1, res.p_to.x);
        max_x1 = std::max(max_x1, res.p_to.x);
        min_y1 = std::min(min_y1, res.p_to.y);
        max_y1 = std::max(max_y1, res.p_to.y);

        // qDebug() << "min_x1: " << min_x1;
        // qDebug() << "max_x1: " << max_x1;
        // qDebug() << "min_y1: " << min_y1;
        // qDebug() << "max_y1: " << max_y1;

        // fmt::println("");

        // create the second projective bivector using positions of both active points
        p = pt2dp(m_p3->scenePos().x, m_p3->scenePos().y, 1.0);
        q = pt2dp(m_p4->scenePos().x, m_p4->scenePos().y, 1.0);

        res = reset_bivecp(p, q);

        // create the lines and map them to the coordinate system
        bvt2_beg_pos = mapFromScene(
            QPointF(cs->x.au_to_w(res.p_from.x), cs->y.au_to_w(res.p_from.y)));
        bvt2_end_pos =
            mapFromScene(QPointF(cs->x.au_to_w(res.p_to.x), cs->y.au_to_w(res.p_to.y)));


        // find min and max values for bounding rectangle:
        auto min_x2 = std::min(std::numeric_limits<double>::max(), res.p_from.x);
        auto max_x2 = std::max(std::numeric_limits<double>::lowest(), res.p_from.x);
        auto min_y2 = std::min(std::numeric_limits<double>::max(), res.p_from.y);
        auto max_y2 = std::max(std::numeric_limits<double>::lowest(), res.p_from.y);
        min_x2 = std::min(min_x2, res.p_to.x);
        max_x2 = std::max(max_x2, res.p_to.x);
        min_y2 = std::min(min_y2, res.p_to.y);
        max_y2 = std::max(max_y2, res.p_to.y);

        // qDebug() << "min_x2: " << min_x2;
        // qDebug() << "max_x2: " << max_x2;
        // qDebug() << "min_y2: " << min_y2;
        // qDebug() << "max_y2: " << max_y2;

        min_x = std::min(min_x1, min_x2);
        max_x = std::max(max_x1, max_x2);
        min_y = std::min(min_y1, min_y2);
        max_y = std::max(max_y1, max_y2);

        // qDebug() << "min_x: " << min_x;
        // qDebug() << "max_x: " << max_x;
        // qDebug() << "min_y: " << min_y;
        // qDebug() << "max_y: " << max_y;
    }
    catch (std::exception const& e) {
        qDebug() << e.what();
    }
}
