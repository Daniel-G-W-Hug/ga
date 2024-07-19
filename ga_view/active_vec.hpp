#pragma once
//
// author: Daniel Hug, 2024
//

#include "active_pt.hpp"
#include "coordsys.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QMarginsF>
#include <QPainter>
#include <QPointF>
#include <QRectF>
#include <QWidget>

// active_vec has two active points. Can be manipulated and moved by mouse.

class active_vec : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 2 };
    int type() const override { return Type; }

    active_vec(Coordsys* cs, w_Coordsys* wcs, active_pt* beg, active_pt* end,
               QGraphicsItem* parent = nullptr);
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_beg(QPointF const& pos);
    void setScenePos_end(QPointF const& pos);
    QPointF scenePos_beg() const;
    QPointF scenePos_end() const;

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

    active_pt* m_beg; // active_pt at beginning position
    active_pt* m_end; // active_pt at end position

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed
};