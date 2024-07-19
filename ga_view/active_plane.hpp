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

// active_plane has two active points. Can be manipulated and moved by mouse.

class active_plane : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 5 };
    int type() const override { return Type; }

    active_plane(Coordsys* cs, w_Coordsys* wcs, active_pt* n1end, active_pt* n2end,
                 QGraphicsItem* parent = nullptr);
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_n1end(QPointF const& pos);
    void setScenePos_n2end(QPointF const& pos);
    QPointF scenePos_n1end() const;
    QPointF scenePos_n2end() const;

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

    active_pt* m_n1end; // active_pt at end position n1 (implicit m_beg at (0,0))
    active_pt* m_n2end; // active_pt at end position n2 (implicit m_beg at (0,0))

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed
};