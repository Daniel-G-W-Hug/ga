#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsView>
#include <QPainter>
#include <QWidget>

// connection line cln2dp - connecting passive projective points

class item_cln2dp : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    item_cln2dp(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm, size_t idx,
                QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  private:

    Coordsys* cs;       // coordsys for drawing
    Coordsys_model* cm; // model to draw from
    size_t idx;         // index in line vector of model to draw from

    double min_x{std::numeric_limits<double>::max()};
    double max_x{std::numeric_limits<double>::lowest()};
    double min_y{std::numeric_limits<double>::max()};
    double max_y{std::numeric_limits<double>::lowest()};
};