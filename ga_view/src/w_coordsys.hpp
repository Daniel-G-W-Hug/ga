#pragma once
//
// author: Daniel Hug, 2024
//

#include "coordsys.hpp"
#include "coordsys_model.hpp"
#include "w_common.hpp"

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPainter>
#include <QWidget>

class w_Coordsys : public QGraphicsView {
    Q_OBJECT

  public:
    w_Coordsys(Coordsys* cs, std::vector<Coordsys_model*>& vm, QGraphicsScene* scene,
               QWidget* parent = nullptr);

  protected:
    void resizeEvent(QResizeEvent* event);
    void paintEvent(QPaintEvent* event);
    void drawBackground(QPainter* qp, const QRectF& rect);
    void drawForeground(QPainter* qp, const QRectF& rect);

    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);

    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);

    void push_to_history();  // for undo
    void pop_from_history(); // undo

  public slots:
    void switch_to_model(size_t);

  signals:
    void viewResized();
    void mouseMoved(bool hot, mouse_pos_t mouse_pos);
    void modeChanged(pz_action action, pz_mode mode);
    void undoChanged(size_t undo_steps);
    void scalingChanged(axis_scal xscal, axis_scal yscal);

  private:
    Coordsys* cs;          // coordinate system
    QGraphicsScene* scene; // graphics scene

    Coordsys_model* cm;              // current modell
    std::vector<Coordsys_model*> vm; // vector of models
                                     // that might be switched between
                                     // in case of several models

    std::vector<Coordsys> cs_history; // history of coordinate-systems (for undo)

    // mouse status
    int m_nx{0};                         // x-position of mouse in widget
    int m_ny{0};                         // y-position of mouse in widget
    int m_nx_hot{0};                     // x-position of mouse in hot area
    int m_ny_hot{0};                     // y-position of mouse in hot area
    bool m_hot{false};                   // mouse within coordsys area
    bool m_leftButton{false};            // left button pressed in hot area
    bool m_rightButton{false};           // right button pressed in hot area
    pz_action m_action{pz_action::none}; // normal movement, no pan or zoom
    pz_mode m_mode{pz_mode::x_and_y};    // no pan or zoom restrictions
    int m_nx_leftPress{0};               // x-position of leftButtonPress-Event
    int m_ny_leftPress{0};               // y-position of leftButtonPress-Event
};
