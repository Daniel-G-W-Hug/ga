//
// author: Daniel Hug, 2024
//

#include "w_mainwindow.hpp"
#include "active_bivec.hpp"
#include "active_plane.hpp"
#include "active_projection.hpp"
#include "active_pt.hpp"
#include "active_vec.hpp"
#include "w_coordsys.hpp"
#include "w_statusbar.hpp"

#include <QPainter>
#include <QPointF>
#include <QVBoxLayout>

Coordsys* get_initial_cs()
{
    // TODO: read default parameters from config file (e.g. .lua)

    axis_data ax(axis_rng(-3.5, 3.5), axis_dir::x, axis_scal::linear, "x label",
                 axis_ticks(0.0, 0.5, 5));

    axis_data ay(axis_rng(-2.5, 2.5), axis_dir::y, axis_scal::linear, "y label",
                 axis_ticks(0.0, 0.5, 5));

    widget_axis_data wx(600, 65, 520);
    widget_axis_data wy(500, 50, 420);

    Axis x(wx, ax);
    Axis y(wy, ay, x.px_density_rng()); // enable aspect_ratio = 1.0

    coordsys_data cd("Coordsys Title");
    Coordsys* cs = new Coordsys(x, y, cd, keep_aspect_ratio::yes);

    return cs;
}

void populate_scene(Coordsys* cs, w_Coordsys* wcs, QGraphicsScene* scene)
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

    // active_pt* pt0 = new active_pt(cs, wcs, QPointF(0, 1));
    // scene->addItem(pt0);
    // active_pt* pt1 = new active_pt(cs, wcs, QPointF(3, 1));
    // scene->addItem(pt1);
    // active_pt* pt2 = new active_pt(cs, wcs, QPointF(0, 2));
    // scene->addItem(pt2);
    // active_bivec* bv1 = new active_bivec(cs, wcs, pt0, pt1, pt2);
    // scene->addItem(bv1);

    ///////////////////////////////////////////////////////////////////////////
    // projection
    ///////////////////////////////////////////////////////////////////////////

    // active_pt* pt0 = new active_pt(cs, wcs, QPointF(0, 1));
    // scene->addItem(pt0);
    // active_pt* pt1 = new active_pt(cs, wcs, QPointF(3, 1));
    // scene->addItem(pt1);
    // active_pt* pt2 = new active_pt(cs, wcs, QPointF(2, 2));
    // scene->addItem(pt2);
    // active_projection* proj1 = new active_projection(cs, wcs, pt0, pt1, pt2);
    // scene->addItem(proj1);

    ///////////////////////////////////////////////////////////////////////////
    // reflection
    ///////////////////////////////////////////////////////////////////////////

    active_pt* pt1 = new active_pt(cs, wcs, QPointF(0, 1));
    scene->addItem(pt1);
    active_pt* pt2 = new active_pt(cs, wcs, QPointF(1, 0));
    scene->addItem(pt2);
    active_plane* pl = new active_plane(cs, wcs, pt1, pt2);
    scene->addItem(pl);
}


w_MainWindow::w_MainWindow(QWidget* parent) : QMainWindow(parent)
{

    cs = get_initial_cs();

    scene = new QGraphicsScene(this);
    // scene->setSceneRect(QRectF(0, 0, cs->x.widget_size(), cs->y.widget_size()));

    wcs = new w_Coordsys(cs, scene, this);
    wsb = new w_Statusbar(cs->x.widget_size(), this);

    auto layout = new QVBoxLayout;
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addWidget(wcs);
    layout->addWidget(wsb);
    // setLayout(layout);

    auto widget = new QWidget;
    widget->setLayout(layout);
    setCentralWidget(widget);
    setWindowTitle("qt-cs");

    // link coordsys to statusbar
    connect(wcs, &w_Coordsys::mouseMoved, wsb, &w_Statusbar::on_mouseMoved);
    connect(wcs, &w_Coordsys::modeChanged, wsb, &w_Statusbar::on_modeChanged);
    connect(wcs, &w_Coordsys::undoChanged, wsb, &w_Statusbar::on_undoChanged);
    connect(wcs, &w_Coordsys::labelChanged, wsb, &w_Statusbar::on_labelChanged);
    connect(wcs, &w_Coordsys::scalingChanged, wsb, &w_Statusbar::on_scalingChanged);

    // // update status bar with label of first model
    // emit wcs->labelChanged(cm->label());
    // emit wcs->labelChanged(vm[0]->label());

    // update status bar with axis scaling
    emit wcs->scalingChanged(cs->x.scaling(), cs->y.scaling());

    populate_scene(cs, wcs, scene);
}
