// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "item_pt2dp.hpp"
#include "active_common.hpp"

#include "QPointF"
#include "QRectF"
#include "QSizeF"

#include <stdexcept> // std::runtime_error

item_pt2dp::item_pt2dp(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                       QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}
{
    // setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
    //          QGraphicsItem::ItemSendsGeometryChanges |
    //          QGraphicsItem::ItemSendsScenePositionChanges);
    // setAcceptHoverEvents(true);

    connect(wcs, &w_Coordsys::viewResized, this, &item_pt2dp::viewChanged);


    if (cm->ptp[idx].z == 0.0) {
        throw std::runtime_error(
            "Point at infinity. pt2dp requires z-component z == 1.0");
    };
    if (cm->ptp[idx].z != 1.0) {
        // if required, unitze the point as first step
        cm->ptp[idx].x /= cm->ptp[idx].z;
        cm->ptp[idx].y /= cm->ptp[idx].z;
        cm->ptp[idx].z = 1.0;
    };
    setPos(cs->x.au_to_w(cm->ptp[idx].x),
           cs->y.au_to_w(cm->ptp[idx].y)); // set item to scene coordinates

    setZValue(50); // passive points should be on top, but below active points
}

void item_pt2dp::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    qp->setPen(cm->ptp_mark[idx].pen);

    switch (cm->ptp_mark[idx].symbol) {
        case Symbol::plus: {
            qp->drawLine(-cm->ptp_mark[idx].nsize, 0, +cm->ptp_mark[idx].nsize, 0);
            qp->drawLine(0, -cm->ptp_mark[idx].nsize, 0, cm->ptp_mark[idx].nsize);
            break;
        }
        case Symbol::cross: {
            qp->drawLine(-cm->ptp_mark[idx].nsize, -cm->ptp_mark[idx].nsize,
                         cm->ptp_mark[idx].nsize, cm->ptp_mark[idx].nsize);
            qp->drawLine(-cm->ptp_mark[idx].nsize, cm->ptp_mark[idx].nsize,
                         cm->ptp_mark[idx].nsize, -cm->ptp_mark[idx].nsize);
            break;
        }
        case Symbol::circle: {
            qp->drawEllipse(QPoint(0, 0), cm->ptp_mark[idx].nsize,
                            cm->ptp_mark[idx].nsize);
            break;
        }
        case Symbol::square: {
            qp->drawLine(-cm->ptp_mark[idx].nsize, -cm->ptp_mark[idx].nsize,
                         cm->ptp_mark[idx].nsize, -cm->ptp_mark[idx].nsize);
            qp->drawLine(cm->ptp_mark[idx].nsize, -cm->ptp_mark[idx].nsize,
                         cm->ptp_mark[idx].nsize, cm->ptp_mark[idx].nsize);
            qp->drawLine(cm->ptp_mark[idx].nsize, cm->ptp_mark[idx].nsize,
                         -cm->ptp_mark[idx].nsize, cm->ptp_mark[idx].nsize);
            qp->drawLine(-cm->ptp_mark[idx].nsize, cm->ptp_mark[idx].nsize,
                         -cm->ptp_mark[idx].nsize, -cm->ptp_mark[idx].nsize);
            break;
        }
    }

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    qp->restore();
}

QRectF item_pt2dp::boundingRect() const
{
    return QRectF(QPointF(-cm->ptp_mark[idx].nsize, -cm->ptp_mark[idx].nsize),
                  QPointF(cm->ptp_mark[idx].nsize, cm->ptp_mark[idx].nsize));
}

QPainterPath item_pt2dp::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void item_pt2dp::viewChanged()
{
    // qDebug() << "item_pt2dp: viewChanged() received.";

    // view changed by external influence, set to position in model
    setPos(cs->x.au_to_w(cm->ptp[idx].x), cs->y.au_to_w(cm->ptp[idx].y));
}