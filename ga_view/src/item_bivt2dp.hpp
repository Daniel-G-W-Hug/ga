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

// item_bivt2dp has three active points. Can be manipulated and moved by mouse.
// it is intended to show projective lines with a support point
// (initially closest point to origin), a rotation option and a translation option

class item_bivt2dp : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    item_bivt2dp(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                 QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;


  public slots:
    void viewChanged(); // the view changed, e.g. resize

  private:

    void reset_item_data(); // update geometry, for initialization of if view changed

    Coordsys* cs;
    w_Coordsys* wcs;

    Coordsys_model* cm; // model to draw from
    size_t idx;         // index in vector of bivectors of model to draw from

    QPointF beg_pos; // beg_pos of line (already mapped to cs)
    QPointF end_pos; // beg_pos of line (already mapped to cs)

    // for bounding box of line
    double min_x{std::numeric_limits<double>::max()};
    double max_x{std::numeric_limits<double>::lowest()};
    double min_y{std::numeric_limits<double>::max()};
    double max_y{std::numeric_limits<double>::lowest()};
};