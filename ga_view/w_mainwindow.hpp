#pragma once
//
// author: Daniel Hug, 2024
//

#include "coordsys.hpp"
#include "w_coordsys.hpp"
#include "w_statusbar.hpp"

#include <QGraphicsScene>
#include <QMainWindow>
#include <QWidget>

class w_MainWindow : public QMainWindow {

    Q_OBJECT

  public:

    explicit w_MainWindow(QWidget* parent = nullptr);

  private:

    Coordsys* cs;
    QGraphicsScene* scene;

    w_Coordsys* wcs;
    w_Statusbar* wsb;
};