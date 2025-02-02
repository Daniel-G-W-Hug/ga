//
// author: Daniel Hug, 2024 & 2025
//


// TODO: active bivector in EGA should be at origin exclusively
//       move away from origion should only be possible in PGA

#include "item_bivt2de.hpp"
#include "active_common.hpp"

#include "ga/ga_ega.hpp"

item_bivt2de::item_bivt2de(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                           QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, cm{cm}, idx{idx}

{
    Q_UNUSED(wcs)

    setZValue(0);
}

void item_bivt2de::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                         QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    // clipping area is active area of coordsys
    qp->setClipRect(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                          cs->y.nmin() - cs->y.nmax()));

    // draw in item coordinate system
    qp->save();


    try {

        // make all available that is needed for GA
        using namespace hd::ga;
        using namespace hd::ga::ega;
    }
    catch (std::exception const& e) {
        qDebug() << e.what();
    }

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

QRectF item_bivt2de::boundingRect() const
{
    // give bounding box in item coordinate system

    QPolygonF polygon;

    // polygon << beg_pos << end_upos << tip_pos << end_vpos << beg_pos;

    return polygon.boundingRect();
}

QPainterPath item_bivt2de::shape() const
{


    // QPainterPath path = vectorShape(beg_pos, end_upos);
    // path += vectorShape(beg_pos, end_vpos);

    QPainterPath path;

    return path;
}
