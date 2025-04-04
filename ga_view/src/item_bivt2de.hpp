#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_pt2d.hpp"
#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QPointF>
#include <QPolygonF>
#include <QRectF>
#include <QWidget>

// item_bivt2de has three active points. Can be manipulated and moved by mouse.
// shows wedge product of two vector u and v: u^v

class item_bivt2de : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    item_bivt2de(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                 QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  private:

    Coordsys* cs;
    w_Coordsys* wcs;

    Coordsys_model* cm; // model to draw from
    size_t idx;         // index in vector of bivectors of model to draw from
};