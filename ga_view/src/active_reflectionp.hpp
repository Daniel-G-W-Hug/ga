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

// active_reflectionp has four active points. Can be manipulated and moved by mouse.
// active projective reflection

class active_reflectionp : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 7 }; // for use in w_common.cpp
    int type() const override { return Type; }

    // active points are used here as to calculate two bivectors
    // bvt1 = wdg(p1,p2) and bvt2 = wdg(p3,p4)
    active_reflectionp(Coordsys* cs, w_Coordsys* wcs, active_pt2d* p1, active_pt2d* p2,
                       active_pt2d* p3, active_pt2d* p4, QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_p1(pt2d const& pos);
    void setScenePos_p2(pt2d const& pos);
    void setScenePos_p3(pt2d const& pos);
    void setScenePos_p4(pt2d const& pos);

    void setScenePos_p1_wo_update(pt2d const& pos);
    void setScenePos_p2_wo_update(pt2d const& pos);
    void setScenePos_p3_wo_update(pt2d const& pos);
    void setScenePos_p4_wo_update(pt2d const& pos);

    pt2d scenePos_p1() const;
    pt2d scenePos_p2() const;
    pt2d scenePos_p3() const;
    pt2d scenePos_p4() const;

    bool isHovered() { return m_mouse_hover; }

  public slots:
    void viewChanged(); // the view changed, e.g. resize
    void moveModeChanged(move_mode mode);

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

    void reset_item_data(); // update geometry, for initialization of if view changed

    Coordsys* cs;
    w_Coordsys* wcs;

    active_pt2d* m_p1; // active_pt2d at position p1
    active_pt2d* m_p2; // active_pt2d at position p2
    active_pt2d* m_p3; // active_pt2d at position p3
    active_pt2d* m_p4; // active_pt2d at position p4

    bool m_mouse_hover{false};                          // mouse is hovering over the item
    bool m_mouse_l_pressed{false};                      // left button mouse is pressed
    bool m_mouse_r_pressed{false};                      // right button mouse is pressed
    move_mode m_move_mode{move_mode::shift_both_lines}; // how to move lines

    QPointF bvt1_beg_pos; // beg_pos of bvt1 line to be drawn (already mapped to cs)
    QPointF bvt1_end_pos; // end_pos of bvt1 line to be drawn (already mapped to cs)
    QPointF bvt2_beg_pos; // beg_pos of bvt2 line to be drawn (already mapped to cs)
    QPointF bvt2_end_pos; // end_pos of bvt2 line to be drawn (already mapped to cs)

    vec2dp tp; // current intersection point (for rotational case w/ mouseMoveEvent)

    // for bounding box of line
    double min_x;
    double max_x;
    double min_y;
    double max_y;
};