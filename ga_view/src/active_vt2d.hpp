#pragma once
//
// author: Daniel Hug, 2024
//

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

// active_vt2d has two active points. Can be manipulated and moved by mouse.

class active_vt2d : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 2 };
    int type() const override { return Type; }

    active_vt2d(Coordsys* cs, w_Coordsys* wcs, active_pt2d* beg, active_pt2d* end,
                QGraphicsItem* parent = nullptr);
    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos_beg(pt2d const& pos);
    void setScenePos_end(pt2d const& pos);
    pt2d scenePos_beg() const;
    pt2d scenePos_end() const;

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

    active_pt2d* m_beg; // active_pt2d at beginning position
    active_pt2d* m_end; // active_pt2d at end position

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed
};