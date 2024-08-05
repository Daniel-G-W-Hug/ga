//
// author: Daniel Hug, 2024
//

#include "item_pt2d.hpp"
#include "active_common.hpp"

#include "QPointF"
#include "QRectF"
#include "QSizeF"

item_pt2d::item_pt2d(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                     QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}
{
    // setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
    //          QGraphicsItem::ItemSendsGeometryChanges |
    //          QGraphicsItem::ItemSendsScenePositionChanges);
    // setAcceptHoverEvents(true);

    connect(wcs, &w_Coordsys::viewResized, this, &item_pt2d::viewChanged);

    setPos(cs->x.au_to_w(cm->pt[idx].x),
           cs->y.au_to_w(cm->pt[idx].y)); // set item to scene coordinates

    setZValue(50); // passive points should be on top, but below active points
}

void item_pt2d::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                      QWidget* widget)
{

    // clipping area is active area of coordsys
    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));

    // draw in item coordinate system
    qp->save();

    qp->setPen(cm->pt_mark[idx].pen);

    switch (cm->pt_mark[idx].symbol) {
        case Symbol::plus: {
            qp->drawLine(-cm->pt_mark[idx].nsize, 0, +cm->pt_mark[idx].nsize, 0);
            qp->drawLine(0, -cm->pt_mark[idx].nsize, 0, cm->pt_mark[idx].nsize);
            break;
        }
        case Symbol::cross: {
            qp->drawLine(-cm->pt_mark[idx].nsize, -cm->pt_mark[idx].nsize,
                         cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize);
            qp->drawLine(-cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize,
                         cm->pt_mark[idx].nsize, -cm->pt_mark[idx].nsize);
            break;
        }
        case Symbol::circle: {
            qp->drawEllipse(QPoint(0, 0), cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize);
            break;
        }
        case Symbol::square: {
            qp->drawLine(-cm->pt_mark[idx].nsize, -cm->pt_mark[idx].nsize,
                         cm->pt_mark[idx].nsize, -cm->pt_mark[idx].nsize);
            qp->drawLine(cm->pt_mark[idx].nsize, -cm->pt_mark[idx].nsize,
                         cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize);
            qp->drawLine(cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize,
                         -cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize);
            qp->drawLine(-cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize,
                         -cm->pt_mark[idx].nsize, -cm->pt_mark[idx].nsize);
            break;
        }
    }

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    qp->restore();
}

QRectF item_pt2d::boundingRect() const
{
    return QRectF(QPointF(-cm->pt_mark[idx].nsize, -cm->pt_mark[idx].nsize),
                  QPointF(cm->pt_mark[idx].nsize, cm->pt_mark[idx].nsize));
}

QPainterPath item_pt2d::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void item_pt2d::viewChanged()
{
    // qDebug() << "item_pt2d: viewChanged() received.";

    // view changed by external influence, set to position in model
    setPos(cs->x.au_to_w(cm->pt[idx].x), cs->y.au_to_w(cm->pt[idx].y));
}