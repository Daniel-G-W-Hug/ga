#pragma once
//
// author: Daniel Hug, 2024
//

#include "coordsys.hpp"
#include "coordsys_model.hpp"

#include "w_coordsys.hpp"
#include "w_statusbar.hpp"

#include <QGraphicsScene>
#include <QMainWindow>
#include <QWidget>

class w_MainWindow : public QMainWindow {

    Q_OBJECT

  public:

    explicit w_MainWindow(QWidget* parent = nullptr);

  public slots:
    void changeModel(int);

  signals:
    void updateModel(int);
    void labelChanged(std::string new_label);

  private:

    Coordsys* cs;

    w_Coordsys* wcs;
    QSlider* slider;
    w_Statusbar* wsb;

    QGraphicsScene* scene;

    std::vector<Coordsys_model> models; // model storage
    std::vector<Coordsys_model*> vm;    // vector of pointer to models
                                        // that might be switched between
                                        // in case of several models
};