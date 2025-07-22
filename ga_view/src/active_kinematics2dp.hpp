#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_pt2d.hpp"
#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QMarginsF>
#include <QPainter>
#include <QPointF>
#include <QRectF>
#include <QWidget>

// active_kinematics2dp has a vector of active points which can be moved by mouse.
// The last of the n points serves as target, while n-1 points constitute the open
// kinematics chain.

class active_kinematics2dp : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    active_kinematics2dp(Coordsys* cs, w_Coordsys* wcs, std::vector<active_pt2d*> pts,
                         QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

  public slots:
    void viewChanged();
    void pointMoved();

  private:

    void solve_kinematic_chain();

    Coordsys* cs;
    w_Coordsys* wcs;

    std::vector<active_pt2d*> m_pts;

    const double segment_length{0.5}; // length of each kinematic segment
};