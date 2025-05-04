// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "item_bivt2dp.hpp"
#include "active_common.hpp"

#include "ga/ga_pga.hpp"
using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)

item_bivt2dp::item_bivt2dp(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                           QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}

{
    Q_UNUSED(wcs)

    connect(wcs, &w_Coordsys::viewResized, this, &item_bivt2dp::viewChanged);

    reset_item_data(); // update geometry (beg_pos, end_pos, and min- & max-values)
    setZValue(0);
}

void item_bivt2dp::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                         QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    // clipping area is active area of coordsys
    qp->setClipRect(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                          cs->y.nmin() - cs->y.nmax()));

    // draw in item coordinate system
    qp->save();

    qp->setPen(cm->bivtp_mark[idx].pen);
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

QRectF item_bivt2dp::boundingRect() const
{
    return QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x), cs->y.au_to_w(max_y))),
                  mapFromScene(QPointF(cs->x.au_to_w(max_x), cs->y.au_to_w(min_y))));
}

QPainterPath item_bivt2dp::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void item_bivt2dp::viewChanged()
{
    // qDebug() << "item_bivt2dp: viewChanged() received.";

    // view changed by external influence, set to position in model
    reset_item_data(); // update geometry (beg_pos, end_pos, and min- & max-values)
}

void item_bivt2dp::reset_item_data()
{

    // update geometry, for initialization of if view changed

    try {
        // determine the angle of the projective line
        auto const x_axis = bivec2dp{0, 1, 0};
        auto const y_axis = bivec2dp{1, 0, 0}; // really this is -y_axis_2dp

        auto phi_x = angle(x_axis, cm->bivtp[idx]);
        auto phi_y = angle(y_axis, cm->bivtp[idx]);

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

        vec2dp p_from;
        vec2dp p_to;
        vec2dp p_from_24;
        vec2dp p_to_24;
        vec2dp p_from_13;
        vec2dp p_to_13;

        // constexpr double eps = 1.e-6;  // not needed here, would shadow the global eps

        if (std::abs(phi_x - pi / 2.0) > eps) {
            // if not on pos. or neg. y-axis
            p_from_24 = unitize(rwdg(l2, cm->bivtp[idx])); // from lower x values
            p_to_24 = unitize(rwdg(l4, cm->bivtp[idx]));   // to higher x values
        }

        if (std::abs(phi_y - pi / 2.0) > eps) {
            // if not on pos.or neg.x - axis
            p_from_13 = unitize(rwdg(l3, cm->bivtp[idx])); // from lower y values
            p_to_13 = unitize(rwdg(l1, cm->bivtp[idx]));   // to higher y values
        }

        // pos. x-axis
        if ((std::abs(phi_x) < eps) && (std::abs(phi_y - pi / 2.0) < eps)) {

            // qDebug() << "On positive x-axis.";

            p_from = p_from_24;
            p_to = p_to_24;
        }

        // first quadrant
        if (phi_x > eps && (phi_x < (pi / 2. - eps)) && (phi_y > (pi / 2. + eps)) &&
            (phi_y < (pi - eps))) {

            // qDebug() << "In first quadrant.";

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

        // pos. y-axis
        if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y - pi) < eps)) {

            // qDebug() << "On positive y-axis.";

            p_from = p_from_13;
            p_to = p_to_13;
        }

        // second quadrant
        if ((phi_x > (pi / 2.0 + eps)) && (phi_x < (pi - eps)) && (phi_y < (pi - eps)) &&
            (phi_y > (pi / 2.0 + eps))) {

            // qDebug() << "In second quadrant.";

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

        // neg. x-axis
        if ((std::abs(phi_x - pi) < eps) && (std::abs(phi_y - pi / 2.0)) < eps) {

            // qDebug() << "On negative x-axis.";

            p_from = p_to_24;
            p_to = p_from_24;
        }

        // third quadrant
        if ((phi_x < pi - eps) && (phi_x > pi / 2.0 + eps) && (phi_y < pi / 2.0 - eps) &&
            (phi_y > eps)) {

            // qDebug() << "in third quadrant.";

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

        // neg. y-axis
        if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y) < eps)) {

            // qDebug() << "On negative y-axis.";

            p_from = p_to_13;
            p_to = p_from_13;
        }

        // fourth quadrant
        if ((phi_x < pi / 2.0 - eps) && (phi_x > eps) && (phi_y > eps) &&
            (phi_y < pi / 2. - eps)) {

            // qDebug() << "In fourth quadrant.";

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
        //              cm->bivtp[idx].x, cm->bivtp[idx].y, cm->bivtp[idx].z);
        // fmt::println("phi_x={:.3g}, phi_y={:.3g}, att={:.4g}", rad2deg(phi_x),
        //              rad2deg(phi_y), att(cm->bivtp[idx]));
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
