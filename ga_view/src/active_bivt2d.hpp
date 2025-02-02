#pragma once
//
// author: Daniel Hug, 2024 & 2025
//

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

// active_bivt2d has three active points. Can be manipulated and moved by mouse.
// shows wedge product of two vector u and v: u^v

class active_bivt2d : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 3 };
    int type() const override { return Type; }

    active_bivt2d(Coordsys* cs, w_Coordsys* wcs, active_pt2d* beg, active_pt2d* uend,
                  active_pt2d* vend, QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_beg(pt2d const& pos);
    void setScenePos_uend(pt2d const& pos);
    void setScenePos_vend(pt2d const& pos);
    pt2d scenePos_beg() const;
    pt2d scenePos_uend() const;
    pt2d scenePos_vend() const;

    bool isHovered() { return m_mouse_hover; }

  protected:

    void hoverEnterEvent(QGraphicsSceneHoverEvent* event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* event) override;

    void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event) override;

    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;

  private:

    Coordsys* cs;
    w_Coordsys* wcs;

    active_pt2d* m_beg;  // active_pt2d at beginning position (common starting point)
    active_pt2d* m_uend; // active_pt2d at end position of u vector
    active_pt2d* m_vend; // active_pt2d at end position of v vector

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed
};