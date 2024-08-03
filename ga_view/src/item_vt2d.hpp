#pragma once
//
// author: Daniel Hug, 2024
//

#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsView>
#include <QPainter>
#include <QPointF>
#include <QWidget>

// active_pt can be moved by mouse

class item_vt2d : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    item_vt2d(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
              QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void viewChanged(); // the view changed, e.g. resize

  private:

    Coordsys* cs;       // coordsys for drawing
    Coordsys_model* cm; // model to draw from
    size_t idx;         // index in line vector of model to draw from

    QPointF beg_pos; // beg_pos of vector (already mapped to cs)
    QPointF end_pos; // beg_pos of vector (already mapped to cs)
};