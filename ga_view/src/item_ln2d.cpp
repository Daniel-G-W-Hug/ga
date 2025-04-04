// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "item_ln2d.hpp"
#include "active_common.hpp"

#include <algorithm> // std::min, std::max
#include <limits>    // std::numeric_limits

#include "QPointF"
#include "QRectF"
#include "QSizeF"

item_ln2d::item_ln2d(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                     QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}
{
    Q_UNUSED(wcs)

    // setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
    //          QGraphicsItem::ItemSendsGeometryChanges |
    //          QGraphicsItem::ItemSendsScenePositionChanges);
    // setAcceptHoverEvents(true);

    // connect(wcs, &w_Coordsys::viewResized, this, &item_ln2d::viewChanged);

    // find min and max values for bounding rectangle:
    // QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x),cs->y.au_to_w(max_y))),
    //        mapFromScene(QPointF(cs->x.au_to_w(max_x),cs->y.au_to_w(min_y))))
    for (size_t j = 0; j < cm->ln[idx].size(); ++j) {
        min_x = std::min(min_x, cm->ln[idx][j].x);
        max_x = std::max(max_x, cm->ln[idx][j].x);
        min_y = std::min(min_y, cm->ln[idx][j].y);
        max_y = std::max(max_y, cm->ln[idx][j].y);
    }

    if (cm->ln_mark[idx].mark_area) {
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

void item_ln2d::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    qp->setPen(cm->ln_mark[idx].pen);

    // connect all points on each line
    for (size_t j = 0; j < cm->ln[idx].size() - 1; ++j) {
        int nx1 = cs->x.au_to_w(cm->ln[idx][j].x);
        int ny1 = cs->y.au_to_w(cm->ln[idx][j].y);
        int nx2 = cs->x.au_to_w(cm->ln[idx][j + 1].x);
        int ny2 = cs->y.au_to_w(cm->ln[idx][j + 1].y);
        qp->drawLine(mapFromScene(QPointF(nx1, ny1)), mapFromScene(QPointF(nx2, ny2)));
    }

    if (cm->ln_mark[idx].mark_area) {
        qp->setPen(cm->ln_mark[idx].pen);
        qp->setBrush(cm->ln_mark[idx].area_col);

        QPainterPath polyPath;
        polyPath.moveTo(
            mapFromScene(QPointF(cs->x.au_to_w(cm->ln[idx][0].x), cs->y.au_to_w(0.0))));

        // connect all points on each line
        for (size_t j = 0; j < cm->ln[idx].size(); ++j) {
            polyPath.lineTo(mapFromScene(QPointF(cs->x.au_to_w(cm->ln[idx][j].x),
                                                 cs->y.au_to_w(cm->ln[idx][j].y))));
        }

        polyPath.lineTo(mapFromScene(QPointF(
            cs->x.au_to_w(cm->ln[idx][cm->ln[idx].size() - 1].x), cs->y.au_to_w(0.0))));
        polyPath.closeSubpath();
        qp->drawPath(polyPath);
    }

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    qp->restore();
}

QRectF item_ln2d::boundingRect() const
{
    return QRectF(mapFromScene(QPointF(cs->x.au_to_w(min_x), cs->y.au_to_w(max_y))),
                  mapFromScene(QPointF(cs->x.au_to_w(max_x), cs->y.au_to_w(min_y))));
}

QPainterPath item_ln2d::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}