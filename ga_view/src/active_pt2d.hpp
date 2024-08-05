#pragma once
//
// author: Daniel Hug, 2024
//

#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_coordsys.hpp"

#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>
#include <QPainter>
#include <QWidget>

// active_pt2d can be moved by mouse

class active_pt2d : public QObject, public QGraphicsItem {

    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

  public:

    enum { Type = UserType + 1 };
    int type() const override { return Type; }

    // interactive changes of pos will be reflected in the model
    active_pt2d(Coordsys* cs, w_Coordsys* wcs, pt2d& pos,
                QGraphicsItem* parent = nullptr);

    void paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

    void setScenePos(pt2d const& pos);
    pt2d scenePos() const;

    bool isHovered() { return m_mouse_hover; }

  public slots:
    void viewChanged(); // the view changed, e.g. resize
    void posChanged();  // the position was changed, e.g. shifted in viewport

  protected:

    void hoverEnterEvent(QGraphicsSceneHoverEvent* event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* event) override;

    void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event) override;

    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;

  private:

    Coordsys* cs; // coordsys for drawing
    pt2d& m_pos;  // position of the center point of the item (as scene position)

    bool m_mouse_hover{false};     // mouse is hovering over the item
    bool m_mouse_l_pressed{false}; // left button mouse is pressed
    bool m_mouse_r_pressed{false}; // right button mouse is pressed
};