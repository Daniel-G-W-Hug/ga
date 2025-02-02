//
// author: Daniel Hug, 2024 & 2025
//

#include "item_vt2d.hpp"
#include "active_common.hpp"

#include <algorithm> // std::min, std::max
#include <limits>    // std::numeric_limits

#include "QPointF"
#include "QRectF"

item_vt2d::item_vt2d(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                     QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}
{

    // setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable |
    //          QGraphicsItem::ItemSendsGeometryChanges |
    //          QGraphicsItem::ItemSendsScenePositionChanges);
    // setAcceptHoverEvents(true);

    connect(wcs, &w_Coordsys::viewResized, this, &item_vt2d::viewChanged);

    // setZValue(0); // passive vectors should be far down in the stack

    beg_pos = mapFromScene(
        QPointF(cs->x.au_to_w(cm->vt[idx].beg.x), cs->y.au_to_w(cm->vt[idx].beg.y)));
    end_pos = mapFromScene(
        QPointF(cs->x.au_to_w(cm->vt[idx].end.x), cs->y.au_to_w(cm->vt[idx].end.y)));
}

void item_vt2d::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
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

    // draw vector from beg to end
    qp->setPen(cm->vt_mark[idx].pen);
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

QRectF item_vt2d::boundingRect() const { return QRectF(beg_pos, end_pos).normalized(); }

QPainterPath item_vt2d::shape() const { return vectorShape(beg_pos, end_pos); }

void item_vt2d::viewChanged()
{
    // qDebug() << "item_vt2d: viewChanged() received.";

    // view changed by external influence, set to position in model
    beg_pos = mapFromScene(
        QPointF(cs->x.au_to_w(cm->vt[idx].beg.x), cs->y.au_to_w(cm->vt[idx].beg.y)));
    end_pos = mapFromScene(
        QPointF(cs->x.au_to_w(cm->vt[idx].end.x), cs->y.au_to_w(cm->vt[idx].end.y)));
}