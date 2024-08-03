#pragma once
//
// author: Daniel Hug, 2024
//

#include "active_pt.hpp"
#include "coordsys.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QPointF>
#include <QPolygonF>
#include <QRectF>
#include <QWidget>

// active_projection has three active points. Can be manipulated and moved by mouse.
// shows wedge product of two vector u and v: u^v

class active_projection : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 4 };
    int type() const override { return Type; }

    active_projection(Coordsys* cs, w_Coordsys* wcs, active_pt* beg, active_pt* uend,
                      active_pt* vend, QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_beg(QPointF const& pos);
    void setScenePos_uend(QPointF const& pos);
    void setScenePos_vend(QPointF const& pos);
    QPointF scenePos_beg() const;
    QPointF scenePos_uend() const;
    QPointF scenePos_vend() const;

    bool isHovered() { return m_mouse_hover; }

  signals:
    void viewMoved();

  protected:

    void hoverEnterEvent(QGraphicsSceneHoverEvent* event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* event) override;

    void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event) override;

    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;

  private:

    Coordsys* cs;
    w_Coordsys* wcs;

    active_pt* m_beg;  // active_pt at beginning position (common starting point)
    active_pt* m_uend; // active_pt at end position of u vector
    active_pt* m_vend; // active_pt at end position of v vector

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed
};