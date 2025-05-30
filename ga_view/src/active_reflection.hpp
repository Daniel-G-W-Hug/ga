#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_pt2d.hpp"
#include "coordsys.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QMarginsF>
#include <QPainter>
#include <QPointF>
#include <QRectF>
#include <QWidget>

// active_reflection has two active points. Can be manipulated and moved by mouse.

class active_reflection : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 5 }; // for use in w_common.cpp
    int type() const override { return Type; }

    // active points are used here as endpoints of normal vectors from (0,0) to endpoint
    // describing a normal to a plane used for reflection
    active_reflection(Coordsys* cs, w_Coordsys* wcs, active_pt2d* n1end,
                      active_pt2d* n2end, QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_n1end(pt2d const& pos);
    void setScenePos_n2end(pt2d const& pos);
    void setScenePos_n1end_wo_update(pt2d const& pos);
    void setScenePos_n2end_wo_update(pt2d const& pos);

    pt2d scenePos_n1end() const;
    pt2d scenePos_n2end() const;

    bool isHovered() { return m_mouse_hover; }

  signals:
    void viewMoved();
    void pointsMoved();

  protected:

    void hoverEnterEvent(QGraphicsSceneHoverEvent* event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* event) override;

    void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event) override;

    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;

  private:

    Coordsys* cs;
    w_Coordsys* wcs;

    active_pt2d* m_n1end; // active_pt2d at end position n1 (implicit m_beg at (0,0))
    active_pt2d* m_n2end; // active_pt2d at end position n2 (implicit m_beg at (0,0))

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed
};