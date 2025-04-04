#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "coordsys.hpp"
#include "coordsys_model.hpp"

#include "w_coordsys.hpp"
#include "w_statusbar.hpp"

#include <QGraphicsScene>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QMainWindow>
#include <QSlider>
#include <QSpinBox>
#include <QWidget>

class w_MainWindow : public QMainWindow {

    Q_OBJECT

  public:

    explicit w_MainWindow(QWidget* parent = nullptr);

  public slots:
    void changeModel(size_t);

  signals:
    void updateCurrentModel(size_t);
    void updateMaximumModel(size_t);
    void labelChanged(std::string new_label);

  private:

    Coordsys* cs;

    w_Coordsys* wcs;
    QSlider* slider;
    QSpinBox* spinbox;
    QGroupBox* sliderSpinGroup;
    w_Statusbar* wsb;

    QGraphicsScene* scene;

    std::vector<Coordsys_model> models; // model storage
    std::vector<Coordsys_model*> vm;    // vector of pointer to models
                                        // that might be switched between
                                        // in case of several models
};