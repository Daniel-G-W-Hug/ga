// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "item_cln2dp.hpp"
#include "active_common.hpp"

#include <algorithm> // std::min, std::max
#include <limits>    // std::numeric_limits

#include "QPointF"
#include "QRectF"
#include "QSizeF"

#include "ga/ga_pga.hpp"
using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


item_cln2dp::item_cln2dp(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                         QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}
{
    Q_UNUSED(wcs)

    // setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
    //          QGraphicsItem::ItemSendsGeometryChanges |
    //          QGraphicsItem::ItemSendsScenePositionChanges);
    // setAcceptHoverEvents(true);

    // connect(wcs, &w_Coordsys::viewResized, this, &item_cln2dp::viewChanged);

    // make sure the points are normalized
    for (size_t j = 0; j < cm->clnp[idx].size(); ++j) {
        if (cm->clnp[idx][j].z == 0.0) {
            throw std::runtime_error("pt2de in ln2 requires z-component z == 1.0");
        }
        if (cm->clnp[idx][j].z != 1.0) {
            // if required, unitze the point as first step
            cm->clnp[idx][j].x /= cm->clnp[idx][j].z;
            cm->clnp[idx][j].y /= cm->clnp[idx][j].z;
            cm->clnp[idx][j].z = 1.0;
        }
    }

    // // create corresponding lines as bivectors
    // std::vector<bivec2dp> bve;
    // for (size_t j = 0; j < cm->clnp[idx].size() - 1; ++j) {
    //     bve.emplace_back(wdg(cm->clnp[idx][j], cm->clnp[idx][j + 1]));
    //     fmt::println("bve.[{}]          = {}", j, bve[j]);
    // }
    // fmt::println("");
    // fmt::println("meet(bve.[0],bve[1]          = {}", meet(bve[0], bve[1]));
    // fmt::println("unitize(meet(bve.[0],bve[1]) = {}", unitize(meet(bve[0], bve[1])));

    // auto lo = vec2dp(cs->x.min(), cs->y.max(), 1.0);
    // auto ro = vec2dp(cs->x.max(), cs->y.max(), 1.0);
    // auto upper_line = wdg(lo, ro);
    // auto cut_at_upper = unitize(meet(upper_line, bve[1]));
    // fmt::println("cut_at_upper = {:g}", cut_at_upper);

    // find min and max values for bounding rectangle:
    // QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x),cs->y.au_to_w(max_y))),
    //        mapFromScene(QPointF(cs->x.au_to_w(max_x),cs->y.au_to_w(min_y))))
    for (size_t j = 0; j < cm->clnp[idx].size(); ++j) {
        min_x = std::min(min_x, cm->clnp[idx][j].x);
        max_x = std::max(max_x, cm->clnp[idx][j].x);
        min_y = std::min(min_y, cm->clnp[idx][j].y);
        max_y = std::max(max_y, cm->clnp[idx][j].y);
    }

    if (cm->clnp_mark[idx].mark_area) {
        // in case of area marking also consider y=0 (all else is taken care of above)
        min_y = std::min(min_y, 0.0);
        max_y = std::max(max_y, 0.0);
    }
    // qDebug() << "min_x: " << min_x;
    // qDebug() << "max_x: " << max_x;
    // qDebug() << "min_y: " << min_y;
    // qDebug() << "max_y: " << max_y;

    setZValue(0); // passive lines should be far down in the stack
}

void item_cln2dp::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    qp->setPen(cm->clnp_mark[idx].pen);

    // connect all points on each line
    for (size_t j = 0; j < cm->clnp[idx].size() - 1; ++j) {
        int nx1 = cs->x.au_to_w(cm->clnp[idx][j].x);
        int ny1 = cs->y.au_to_w(cm->clnp[idx][j].y);
        int nx2 = cs->x.au_to_w(cm->clnp[idx][j + 1].x);
        int ny2 = cs->y.au_to_w(cm->clnp[idx][j + 1].y);
        qp->drawLine(mapFromScene(QPointF(nx1, ny1)), mapFromScene(QPointF(nx2, ny2)));
    }

    if (cm->clnp_mark[idx].mark_area) {
        qp->setPen(cm->clnp_mark[idx].pen);
        qp->setBrush(cm->clnp_mark[idx].area_col);

        QPainterPath polyPath;
        polyPath.moveTo(
            mapFromScene(QPointF(cs->x.au_to_w(cm->clnp[idx][0].x), cs->y.au_to_w(0.0))));

        // connect all points on each line
        for (size_t j = 0; j < cm->clnp[idx].size(); ++j) {
            polyPath.lineTo(mapFromScene(QPointF(cs->x.au_to_w(cm->clnp[idx][j].x),
                                                 cs->y.au_to_w(cm->clnp[idx][j].y))));
        }

        polyPath.lineTo(
            mapFromScene(QPointF(cs->x.au_to_w(cm->clnp[idx][cm->clnp[idx].size() - 1].x),
                                 cs->y.au_to_w(0.0))));
        polyPath.closeSubpath();
        qp->drawPath(polyPath);
    }

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    qp->restore();
}

QRectF item_cln2dp::boundingRect() const
{
    return QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x), cs->y.au_to_w(max_y))),
                  mapFromScene(QPointF(cs->x.au_to_w(max_x), cs->y.au_to_w(min_y))));
}

QPainterPath item_cln2dp::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}