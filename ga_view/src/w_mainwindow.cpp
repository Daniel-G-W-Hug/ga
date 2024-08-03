//
// author: Daniel Hug, 2024
//

#include "w_mainwindow.hpp"
#include "w_coordsys.hpp"
#include "w_statusbar.hpp"

#include "active_bivec.hpp"
#include "active_projection.hpp"
#include "active_pt.hpp"
#include "active_reflection.hpp"
#include "active_vec.hpp"

#include "item_ln2d.hpp"
#include "item_pt2d.hpp"
#include "item_vt2d.hpp"

#include "coordsys_model.hpp"
#include "hd/hd_functions.hpp"

#include <QPainter>
#include <QPointF>
#include <QVBoxLayout>

#include <stdexcept> // std::runtime_error

std::vector<Coordsys_model> get_single_model_with_lots_of_stuff()
{

    Coordsys_model cm;
    std::vector<Coordsys_model> vm;

    // pt2d p0;
    // p0.y += 0.1;
    // cm.add_pt(p0);

    // pt2d p1;
    // p1.x = 3.0;
    // p1.y = 2.0;
    // cm.add_pt(p1);

    // pt2d p2(4.0, 3.0);
    // pt2d_mark p2m;
    // p2m.symbol = Symbol::circle;
    // p2m.pen = QPen(Qt::green, 2, Qt::SolidLine);
    // cm.add_pt(p2, p2m);

    // pt2d p3{4.0, 3.0};
    // pt2d_mark p3m;
    // p3m.symbol = Symbol::square;
    // p3m.pen = QPen(Qt::red, 2, Qt::SolidLine);
    // cm.add_pt(p3, p3m);

    // ln2d l1;
    // l1.push_back(p0);
    // l1.push_back(p3);
    // cm.add_ln(l1);

    // {
    //     ln2d l2;
    //     double x = -0.5;
    //     double dx = 0.01;
    //     double x_eps = 0.1 * dx;

    //     while (x <= 1.5 + x_eps) {
    //         l2.push_back(pt2d(x, hd::linear_step(0.0, 1.0, x)));
    //         x += dx;
    //     }
    //     ln2d_mark l2m;
    //     l2m.mark_pts = true;
    //     // l2m.delta = 10; // show every tenth point only
    //     l2m.delta = 5;
    //     l2m.pm.symbol = circle;
    //     l2m.pm.pen = QPen(Qt::green, 2, Qt::SolidLine);
    //     cm.add_ln(l2, l2m);
    // }

    // {

    //     ln2d l2;
    //     double x = -0.5;
    //     double dx = 0.01;
    //     double x_eps = 0.1 * dx;

    //     while (x <= 1.5 + x_eps) {
    //         l2.push_back(pt2d(x, hd::smooth_step(0.0, 1.0, x)));
    //         x += dx;
    //     }
    //     ln2d_mark l2m;
    //     l2m.mark_pts = true;
    //     l2m.delta = 10; // show every tenth point only
    //     l2m.pm.symbol = plus;
    //     l2m.pm.pen = QPen(Qt::cyan, 2, Qt::SolidLine);

    //     cm.add_ln(l2, l2m);
    // }

    // {

    //     ln2d l2;
    //     double x = -0.5;
    //     double dx = 0.01;
    //     double x_eps = 0.1 * dx;

    //     while (x <= 1.5 + x_eps) {
    //         l2.push_back(pt2d(x, hd::smoother_step(0.0, 1.0, x)));
    //         x += dx;
    //     }
    //     ln2d_mark l2m;
    //     l2m.mark_pts = true;
    //     l2m.delta = 10; // show every tenth point only
    //     l2m.pm.symbol = square;
    //     l2m.pm.pen = QPen(Qt::blue, 2, Qt::SolidLine);

    //     cm.add_ln(l2, l2m);
    // }

    // {

    //     ln2d l2;
    //     double x = 0.25;
    //     double dx = 0.01;
    //     double x_eps = 0.1 * dx;

    //     while (x <= .75 + x_eps) {
    //         l2.push_back(pt2d(x, hd::smoother_step(0.0, 1.0, x)));
    //         x += dx;
    //     }
    //     ln2d_mark l2m;
    //     l2m.mark_area = true;
    //     l2m.pen = QPen(QColor(0, 0, 128, 50), 2, Qt::SolidLine);
    //     l2m.area_col = QColor(0, 0, 128, 50);

    //     cm.add_ln(l2, l2m);
    // }

    {
        pt2d p1(1, 1);
        pt2d p2(2, 1);

        vt2d v1{p1}, v2{p2}, v3{p1, p2};

        cm.add_vt(v1);
        cm.add_vt(v2);
        cm.add_vt(v3);
    }

    cm.set_label("init label");

    vm.push_back(cm);

    return vm;
}

std::vector<Coordsys_model> get_moving_line()
{

    std::vector<Coordsys_model> vm;

    ln2d_mark lm;
    lm.pen = QPen(Qt::red, 2, Qt::SolidLine);
    lm.mark_pts = true;
    lm.pm.symbol = Symbol::circle;
    lm.pm.nsize = 4;
    lm.pm.pen = QPen(Qt::green, 1, Qt::SolidLine);

    double tmin = 0.0;
    double tmax = 2.0;
    double dt = 0.01;
    double t_eps = 0.01 * dt;

    double T = 1.0;
    double lambda = 2.0;
    //
    double omega = 2. * M_PI / T;
    double k = 2. * M_PI / lambda;

    double xmin = -2.0;
    double xmax = 2.0;
    double dx = 0.05;
    double x_eps = 0.01 * dx;

    // rechtslaufende Welle (Schwingungsdauer T, Wellenlänge lambda)
    // y(x,t) = y_max * sin(2*pi*(t/T - x/lambda))
    // y(x,t) = y_max * sin(omega*t - k*x)
    //
    // es gelten: omega = 2*pi/T; f = 1/T; Wellenzahl k = 2*pi/lambda
    // Ausbreitungsgeschwindigkeit: c = lambda*f => c*T = lambda

    for (double t = tmin; t <= tmax + t_eps; t += dt) {

        Coordsys_model cm;

        // add an active point
        if (t == tmin) {
            cm.add_apt(pt2d{1, 1});
        }

        // define a poly line and a model
        ln2d l;
        for (double x = xmin; x <= xmax + x_eps; x += dx) {
            l.push_back(pt2d(x, std::sin(omega * t - k * x)));
        }

        // add the poly line to the model
        cm.add_ln(l, lm);

        cm.set_label(fmt::format("t={:.3f}", t));

        // store the model in the vector
        vm.push_back(cm);
    }

    return vm;
}

void populate_scene(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm,
                    QGraphicsScene* scene)
{

    ///////////////////////////////////////////////////////////////////////////
    // vector
    ///////////////////////////////////////////////////////////////////////////

    // active_pt* pt0 = new active_pt(cs, wcs, QPointF(0, 1));
    // scene->addItem(pt0);
    // active_pt* pt1 = new active_pt(cs, wcs, QPointF(3, 1));
    // scene->addItem(pt1);
    // active_pt* pt2 = new active_pt(cs, wcs, QPointF(0, 2));
    // scene->addItem(pt2);

    // active_vec* v0 = new active_vec(cs, wcs, pt0, pt1);
    // scene->addItem(v0);
    // active_vec* v1 = new active_vec(cs, wcs, pt1, pt2);
    // scene->addItem(v1);

    ///////////////////////////////////////////////////////////////////////////
    // bivector
    ///////////////////////////////////////////////////////////////////////////

    // active_pt* pt0 = new active_pt(cs, wcs, QPointF(0, 0));
    // scene->addItem(pt0);
    // active_pt* pt1 = new active_pt(cs, wcs, QPointF(3, 0));
    // scene->addItem(pt1);
    // active_pt* pt2 = new active_pt(cs, wcs, QPointF(0, 2));
    // scene->addItem(pt2);
    // active_bivec* bv1 = new active_bivec(cs, wcs, pt0, pt1, pt2);
    // scene->addItem(bv1);

    ///////////////////////////////////////////////////////////////////////////
    // projection
    ///////////////////////////////////////////////////////////////////////////

    // active_pt* pt0 = new active_pt(cs, wcs, QPointF(0, 0));
    // scene->addItem(pt0);
    // active_pt* pt1 = new active_pt(cs, wcs, QPointF(3, 0));
    // scene->addItem(pt1);
    // active_pt* pt2 = new active_pt(cs, wcs, QPointF(2, 2));
    // scene->addItem(pt2);
    // active_projection* proj1 = new active_projection(cs, wcs, pt0, pt1, pt2);
    // scene->addItem(proj1);

    ///////////////////////////////////////////////////////////////////////////
    // reflection
    ///////////////////////////////////////////////////////////////////////////

    // active_pt* pt1 = new active_pt(cs, wcs, QPointF(0, 1));
    // scene->addItem(pt1);
    // active_pt* pt2 = new active_pt(cs, wcs, QPointF(1, 0));
    // scene->addItem(pt2);
    // active_reflection* pl = new active_reflection(cs, wcs, pt1, pt2);
    // scene->addItem(pl);

    // remove content from previous scene
    scene->clear();

    // register everything from current model with scene

    // register all vectors with scene
    for (int idx = 0; idx < cm->vt.size(); ++idx) {
        scene->addItem(new item_vt2d(cs, wcs, cm, idx));
    }

    // register all lines with scene
    for (int idx = 0; idx < cm->ln.size(); ++idx) {
        scene->addItem(new item_ln2d(cs, wcs, cm, idx));
    }

    // register all points with scene
    for (int idx = 0; idx < cm->pt.size(); ++idx) {
        scene->addItem(new item_pt2d(cs, wcs, cm, idx));
    }

    // register all active points with scene
    // for (int idx = 0; idx < cm->apt.size(); ++idx) {
    //     scene->addItem(new active_pt(cs, wcs, cm->apt[idx]));
    // }
}

Coordsys* get_initial_cs()
{
    // TODO: read default parameters from config file (e.g. .lua)

    axis_data ax(axis_rng(-3.5, 3.5), axis_dir::x, axis_scal::linear, "x label",
                 axis_ticks(0.0, 0.5, 5));

    axis_data ay(axis_rng(-2.5, 2.5), axis_dir::y, axis_scal::linear, "y label",
                 axis_ticks(0.0, 0.5, 5));

    widget_axis_data wx(900, 65, 820);
    widget_axis_data wy(700, 50, 620);

    Axis x(wx, ax);
    Axis y(wy, ay, x.px_density_rng()); // enable aspect_ratio = 1.0

    coordsys_data cd("Coordsys Title");

    return new Coordsys(x, y, cd, keep_aspect_ratio::yes);
}


w_MainWindow::w_MainWindow(QWidget* parent) : QMainWindow(parent)
{
    // models = get_single_model_with_lots_of_stuff();
    models = get_moving_line();
    for (int i = 0; i < models.size(); ++i) {
        vm.push_back(&models[i]);
        // fmt::println("Model: {}, Label: {}", i, vm[i]->label());
    }

    cs = get_initial_cs();

    scene = new QGraphicsScene(this);

    wcs = new w_Coordsys(cs, vm, scene, this);
    wsb = new w_Statusbar(cs->x.widget_size(), this);

    slider = new QSlider(Qt::Horizontal, this);
    slider->setRange(0,
                     vm.size() - 1); // only allow to switch between existing models

    auto layout = new QVBoxLayout;
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addWidget(wcs);
    layout->addSpacing(2);
    layout->addWidget(slider);
    layout->addSpacing(2);
    layout->addWidget(wsb);
    // setLayout(layout);

    auto widget = new QWidget;
    widget->setLayout(layout);
    setCentralWidget(widget);
    setWindowTitle("ga_view");

    // link mainwindow model changes to coordinate system and status bar
    connect(slider, &QAbstractSlider::valueChanged, this, &w_MainWindow::changeModel);
    connect(this, &w_MainWindow::updateModel, wcs, &w_Coordsys::switch_to_model);
    connect(this, &w_MainWindow::updateModel, wsb, &w_Statusbar::on_modelChanged);
    connect(this, &w_MainWindow::labelChanged, wsb, &w_Statusbar::on_labelChanged);

    // link coordsys to statusbar
    connect(wcs, &w_Coordsys::mouseMoved, wsb, &w_Statusbar::on_mouseMoved);
    connect(wcs, &w_Coordsys::modeChanged, wsb, &w_Statusbar::on_modeChanged);
    connect(wcs, &w_Coordsys::undoChanged, wsb, &w_Statusbar::on_undoChanged);
    connect(wcs, &w_Coordsys::scalingChanged, wsb, &w_Statusbar::on_scalingChanged);

    // update status bar with current axis scaling
    emit wcs->scalingChanged(cs->x.scaling(), cs->y.scaling());

    // update with first model
    if (vm.size() < 1) throw std::runtime_error("w_MainWindow requires model size >= 1.");
    changeModel(0);
}

void w_MainWindow::changeModel(int idx)
{
    if (idx >= 0 && idx < vm.size()) {
        // fmt::println("w_MainWindow got signal model_changed with value {}", idx);
        populate_scene(cs, wcs, vm[idx], scene); // set new scene with model content
        emit labelChanged(vm[idx]->label());     // inform status bar
        emit updateModel(idx);                   // inform coordinate system widget
    }
}
