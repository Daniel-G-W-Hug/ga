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

// active_bivt2dp has two active points. Can be manipulated and moved by mouse.

class active_bivt2dp : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 6 }; // for use in w_common.cpp
    int type() const override { return Type; }

    active_bivt2dp(Coordsys* cs, w_Coordsys* wcs, active_pt2d* beg, active_pt2d* end,
                   QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_beg(pt2d const& pos);
    void setScenePos_end(pt2d const& pos);
    pt2d scenePos_beg() const;
    pt2d scenePos_end() const;

    bool isHovered() const { return m_mouse_hover; }

  public slots:
    void viewChanged(); // the view changed, e.g. resize

  signals:
    void viewMoved();

  protected:

    void hoverEnterEvent(QGraphicsSceneHoverEvent* event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* event) override;

    void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event) override;

    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;

  private:

    void reset_item_data(); // update geometry, for initialization of if view changed

    Coordsys* cs;
    w_Coordsys* wcs;

    active_pt2d* m_beg; // active_pt2d at beginning position (active point)
    active_pt2d* m_end; // active_pt2d at end position (active point)

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed

    QPointF beg_pos; // beg_pos of line to be drawn (already mapped to cs)
    QPointF end_pos; // end_pos of line to be drawn (already mapped to cs)

    // for bounding box of line
    double min_x{std::numeric_limits<double>::max()};
    double max_x{std::numeric_limits<double>::lowest()};
    double min_y{std::numeric_limits<double>::max()};
    double max_y{std::numeric_limits<double>::lowest()};
};