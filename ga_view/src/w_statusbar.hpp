#pragma once
//
// author: Daniel Hug, 2024
//

#include "coordsys.hpp"
#include "w_coordsys.hpp"

#include <QPainter>
#include <QWidget>
#include <QtWidgets>

class w_Statusbar : public QWidget {
    Q_OBJECT

  public:

    w_Statusbar(int width, QWidget* parent = nullptr);

  protected:

    void resizeEvent(QResizeEvent* event);
    void paintEvent(QPaintEvent* event);
    void draw(QPainter* qp);

  public slots:

    // mouse within (min...max) (true) or not (false) cs area, current position is
    // x, y
    void on_mouseMoved(bool hot, mouse_pos_t mouse_pos);
    void on_currentModelChanged(int currentModelIndex);
    void on_maximumModelChanged(int maximumModelIndex);
    void on_modeChanged(pz_action action, pz_mode mode);
    void on_undoChanged(int undo_steps);
    void on_labelChanged(std::string label);
    void on_scalingChanged(axis_scal xscal, axis_scal yscal);

  private:

    int w_width;
    const int w_height{20};

    // data to be displayed in status bar

    // mouse position within coordsys
    bool m_hot;      // mouse is within cs area
    int m_nx, m_ny;  // mouse position in device coordinates
    double m_x, m_y; // mouse position in cs


    int m_currentModel{0}; // model index (default: show first step)
    int m_maximumModel{0}; // maximum index of models

    // model label
    std::string m_label{};

    // pan and zoom action
    pz_action m_action{pz_action::none};

    // pan and zoom mode restrictions
    // (x_and_y: no restriction)
    pz_mode m_mode{pz_mode::x_and_y};

    // number of undo steps available
    int m_undo_steps{0};

    // axis scaling
    axis_scal m_xscaling{axis_scal::linear};
    axis_scal m_yscaling{axis_scal::linear};
};
