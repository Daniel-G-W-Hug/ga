// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "item_pt2de.hpp"
#include "active_common.hpp"

#include "QPointF"
#include "QRectF"
#include "QSizeF"

#include <stdexcept> // std::runtime_error

item_pt2de::item_pt2de(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                       QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}
{
    // setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
    //          QGraphicsItem::ItemSendsGeometryChanges |
    //          QGraphicsItem::ItemSendsScenePositionChanges);
    // setAcceptHoverEvents(true);

    connect(wcs, &w_Coordsys::viewResized, this, &item_pt2de::viewChanged);


    if (cm->pte[idx].z == 0.0) {
        throw std::runtime_error(
            "Point at infinity. pt2de requires z-component z == 1.0");
    };
    if (cm->pte[idx].z != 1.0) {
        // if required, unitze the point as first step
        cm->pte[idx].x /= cm->pte[idx].z;
        cm->pte[idx].y /= cm->pte[idx].z;
        cm->pte[idx].z = 1.0;
    };
    setPos(cs->x.au_to_w(cm->pte[idx].x),
           cs->y.au_to_w(cm->pte[idx].y)); // set item to scene coordinates

    setZValue(50); // passive points should be on top, but below active points
}

void item_pt2de::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    qp->setPen(cm->pte_mark[idx].pen);

    switch (cm->pte_mark[idx].symbol) {
        case Symbol::plus: {
            qp->drawLine(-cm->pte_mark[idx].nsize, 0, +cm->pte_mark[idx].nsize, 0);
            qp->drawLine(0, -cm->pte_mark[idx].nsize, 0, cm->pte_mark[idx].nsize);
            break;
        }
        case Symbol::cross: {
            qp->drawLine(-cm->pte_mark[idx].nsize, -cm->pte_mark[idx].nsize,
                         cm->pte_mark[idx].nsize, cm->pte_mark[idx].nsize);
            qp->drawLine(-cm->pte_mark[idx].nsize, cm->pte_mark[idx].nsize,
                         cm->pte_mark[idx].nsize, -cm->pte_mark[idx].nsize);
            break;
        }
        case Symbol::circle: {
            qp->drawEllipse(QPoint(0, 0), cm->pte_mark[idx].nsize,
                            cm->pte_mark[idx].nsize);
            break;
        }
        case Symbol::square: {
            qp->drawLine(-cm->pte_mark[idx].nsize, -cm->pte_mark[idx].nsize,
                         cm->pte_mark[idx].nsize, -cm->pte_mark[idx].nsize);
            qp->drawLine(cm->pte_mark[idx].nsize, -cm->pte_mark[idx].nsize,
                         cm->pte_mark[idx].nsize, cm->pte_mark[idx].nsize);
            qp->drawLine(cm->pte_mark[idx].nsize, cm->pte_mark[idx].nsize,
                         -cm->pte_mark[idx].nsize, cm->pte_mark[idx].nsize);
            qp->drawLine(-cm->pte_mark[idx].nsize, cm->pte_mark[idx].nsize,
                         -cm->pte_mark[idx].nsize, -cm->pte_mark[idx].nsize);
            break;
        }
    }

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    qp->restore();
}

QRectF item_pt2de::boundingRect() const
{
    return QRectF(QPointF(-cm->pte_mark[idx].nsize, -cm->pte_mark[idx].nsize),
                  QPointF(cm->pte_mark[idx].nsize, cm->pte_mark[idx].nsize));
}

QPainterPath item_pt2de::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void item_pt2de::viewChanged()
{
    // qDebug() << "item_pt2de: viewChanged() received.";

    // view changed by external influence, set to position in model
    setPos(cs->x.au_to_w(cm->pte[idx].x), cs->y.au_to_w(cm->pte[idx].y));
}