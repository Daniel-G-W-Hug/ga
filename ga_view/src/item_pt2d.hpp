#pragma once
//
// author: Daniel Hug, 2024 & 2025
//

#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsView>
#include <QPainter>
#include <QWidget>

// passive pt2d

class item_pt2d : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    item_pt2d(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
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
    size_t idx;         // index in pt vector of model to draw from
};