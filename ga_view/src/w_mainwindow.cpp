// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "w_mainwindow.hpp"
#include "w_coordsys.hpp"
#include "w_statusbar.hpp"

#include "scenes.hpp"

#include <QVBoxLayout>

#include <stdexcept> // std::runtime_error

w_MainWindow::w_MainWindow(QWidget* parent) : QMainWindow(parent)
{
    models = build_models();

    for (size_t i = 0; i < models.size(); ++i) {
        vm.push_back(&models[i]);
        // fmt::println("Model: {}, Label: {}", i, vm[i]->label());
    }
    // fmt::println("vm.size() = {}", vm.size());

    if (vm.size() < 1) throw std::runtime_error("w_MainWindow requires model size >= 1.");

    cs = get_initial_cs();

    scene = new QGraphicsScene(this);

    wcs = new w_Coordsys(cs, vm, scene, this);
    wsb = new w_Statusbar(cs->x.widget_size(), this);

    sliderSpinGroup = new QGroupBox(this);

    slider = new QSlider(Qt::Horizontal, sliderSpinGroup);
    spinbox = new QSpinBox(sliderSpinGroup);

    slider->setRange(0, int(vm.size() - 1));
    slider->setFocusPolicy(Qt::StrongFocus); // allows movement with keys when in focus

    spinbox->setRange(0, int(vm.size() - 1));
    spinbox->setFocusPolicy(Qt::StrongFocus); // allows movement with keys when in focus

    auto hb = new QHBoxLayout;
    hb->setContentsMargins(0, 0, 0, 0);
    hb->addWidget(slider);
    hb->addSpacing(5);
    hb->addWidget(spinbox);
    hb->addSpacing(10);

    sliderSpinGroup->setLayout(hb);

    auto layout = new QVBoxLayout;
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addWidget(wcs);
    layout->addSpacing(-1);
    layout->addWidget(sliderSpinGroup);
    layout->addSpacing(-3);
    layout->addWidget(wsb);
    // setLayout(layout);

    auto widget = new QWidget;
    widget->setLayout(layout);
    setCentralWidget(widget);
    setWindowTitle("ga_view");

    // make the group of slider and spinbox work as a unit
    connect(slider, &QSlider::valueChanged, spinbox, &QSpinBox::setValue);
    connect(spinbox, &QSpinBox::valueChanged, slider, &QSlider::setValue);

    // link mainwindow model changes to coordinate system and status bar
    connect(slider, &QSlider::valueChanged, this, &w_MainWindow::changeModel);
    connect(this, &w_MainWindow::updateCurrentModel, wcs, &w_Coordsys::switch_to_model);
    connect(this, &w_MainWindow::updateCurrentModel, wsb,
            &w_Statusbar::on_currentModelChanged);
    connect(this, &w_MainWindow::updateMaximumModel, wsb,
            &w_Statusbar::on_maximumModelChanged);
    connect(this, &w_MainWindow::labelChanged, wsb, &w_Statusbar::on_labelChanged);

    // link coordsys to statusbar
    connect(wcs, &w_Coordsys::mouseMoved, wsb, &w_Statusbar::on_mouseMoved);
    connect(wcs, &w_Coordsys::modeChanged, wsb, &w_Statusbar::on_modeChanged);
    connect(wcs, &w_Coordsys::undoChanged, wsb, &w_Statusbar::on_undoChanged);
    connect(wcs, &w_Coordsys::scalingChanged, wsb, &w_Statusbar::on_scalingChanged);

    // update status bar with current axis scaling
    emit wcs->scalingChanged(cs->x.scaling(), cs->y.scaling());

    // update status bar with maximum model number
    emit updateMaximumModel(vm.size() - 1);

    // start with the first model
    changeModel(size_t(0));
}

void w_MainWindow::changeModel(size_t idx)
{
    if (idx < vm.size()) {
        // fmt::println("w_MainWindow got signal model_changed with value {}", idx);
        populate_scene(cs, wcs, vm[idx], scene); // set new scene with model content
        emit labelChanged(vm[idx]->label());     // inform status bar
        emit updateCurrentModel(idx);            // inform coordinate system widget
    }
}
