// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "w_coordsys.hpp"
#include "w_mainwindow.hpp"
#include "w_statusbar.hpp"


#include "active_bivt2d.hpp"
#include "active_bivt2dp.hpp"
#include "active_projection.hpp"
#include "active_pt2d.hpp"
#include "active_reflection.hpp"
#include "active_reflectionp.hpp"
#include "active_vt2d.hpp"

#include "item_bivt2dp.hpp"
#include "item_cln2dp.hpp"
#include "item_ln2d.hpp"
#include "item_pt2d.hpp"
#include "item_pt2dp.hpp"
#include "item_vt2d.hpp"

#include "coordsys_model.hpp"
#include "hd/hd_functions.hpp"

#include <QColor>
#include <QPainter>
#include <QPointF>
#include <QVBoxLayout>

#include <map>       // std::map
#include <stdexcept> // std::runtime_error

std::vector<Coordsys_model> get_model_with_lots_of_stuff()
{

    std::vector<Coordsys_model> vm;

    {
        Coordsys_model cm;

        pt2d p0;
        p0.y += 0.1;
        cm.add_pt(p0);

        pt2d p1;
        p1.x = 3.0;
        p1.y = 2.0;
        cm.add_pt(p1);

        pt2d p2(4.0, 3.0);
        pt2d_mark p2m;
        p2m.symbol = Symbol::circle;
        p2m.pen = QPen(Qt::green, 2, Qt::SolidLine);
        cm.add_pt(p2, p2m);

        pt2d p3{4.0, 3.0};
        pt2d_mark p3m;
        p3m.symbol = Symbol::square;
        p3m.pen = QPen(Qt::red, 2, Qt::SolidLine);
        cm.add_pt(p3, p3m);

        ln2d l1;
        l1.push_back(p0);
        l1.push_back(p3);
        cm.add_ln(l1);

        {
            ln2d l2;
            double x = -0.5;
            double dx = 0.01;
            double x_eps = 0.1 * dx;

            while (x <= 1.5 + x_eps) {
                l2.emplace_back(pt2d(x, hd::linear_step(0.0, 1.0, x)));
                x += dx;
            }
            ln2d_mark l2m;
            l2m.mark_pts = true;
            // l2m.delta = 10; // show every tenth point only
            l2m.delta = 5;
            l2m.pm.symbol = circle;
            l2m.pm.pen = QPen(Qt::green, 2, Qt::SolidLine);
            cm.add_ln(l2, l2m);
        }

        {

            ln2d l2;
            double x = -0.5;
            double dx = 0.01;
            double x_eps = 0.1 * dx;

            while (x <= 1.5 + x_eps) {
                l2.emplace_back(pt2d(x, hd::smooth_step(0.0, 1.0, x)));
                x += dx;
            }
            ln2d_mark l2m;
            l2m.mark_pts = true;
            l2m.delta = 10; // show every tenth point only
            l2m.pm.symbol = plus;
            l2m.pm.pen = QPen(Qt::cyan, 2, Qt::SolidLine);

            cm.add_ln(l2, l2m);
        }

        {

            ln2d l2;
            double x = -0.5;
            double dx = 0.01;
            double x_eps = 0.1 * dx;

            while (x <= 1.5 + x_eps) {
                l2.emplace_back(pt2d(x, hd::smoother_step(0.0, 1.0, x)));
                x += dx;
            }
            ln2d_mark l2m;
            l2m.mark_pts = true;
            l2m.delta = 10; // show every tenth point only
            l2m.pm.symbol = square;
            l2m.pm.pen = QPen(Qt::blue, 2, Qt::SolidLine);

            cm.add_ln(l2, l2m);
        }

        {

            ln2d l2;
            double x = 0.25;
            double dx = 0.01;
            double x_eps = 0.1 * dx;

            while (x <= .75 + x_eps) {
                l2.emplace_back(pt2d(x, hd::smoother_step(0.0, 1.0, x)));
                x += dx;
            }
            ln2d_mark l2m;
            l2m.mark_area = true;
            l2m.pen = QPen(QColor(0, 0, 128, 50), 2, Qt::SolidLine);
            l2m.area_col = QColor(0, 0, 128, 50);

            cm.add_ln(l2, l2m);
        }

        cm.set_label("initial model");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2d p0(0, 0);
        pt2d p1(1, 1);
        pt2d p2(2, 1);

        cm.add_pt(p0);
        cm.add_pt(p1);
        cm.add_pt(p2);

        vt2d v1{p0, p1}, v2{p0, p2}, v3{p1, p2};

        cm.add_vt(v1);
        cm.add_vt(v2);
        cm.add_vt(v3);

        cm.set_label("vector model 1");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2dp p0(0, 0, 1);
        pt2dp p1(1, 1.5, 1);
        pt2dp p2(2, 1, 1);

        pt2dp p3(2, 1, 2); // projective point with z != 1
        pt2d_mark p3m;
        p3m.symbol = Symbol::square;
        p3m.pen = QPen(Qt::red, 2, Qt::SolidLine);

        cm.add_pt(p0);
        cm.add_pt(p1);
        cm.add_pt(p2);
        cm.add_pt(p3, p3m);

        cln2dp l1; // connecting line between projective points
        l1.push_back(p0);
        l1.push_back(p2);
        l1.push_back(p1);

        ln2d_mark lm;
        lm.mark_pts = true;
        lm.pen = QPen(Qt::green, 1, Qt::SolidLine);
        lm.pm.symbol = circle;
        lm.pm.pen = QPen(Qt::green, 2, Qt::SolidLine);

        cm.add_ln(l1, lm);

        cm.set_label("projective model 1");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        for (int i = 0; i <= 23; ++i) {
            double phi = i * pi / 12 + pi / 2.;
            bivt2dp b(cos(phi), sin(phi), 0); // lines through origin
            cm.add_bivtp(b);
        }

        cm.set_label("proj. model 2 - lines through orgin");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        for (int i = 0; i <= 23; ++i) {
            double phi = i * pi / 12 + pi / 2.;
            bivt2dp b(cos(phi), sin(phi), 1.5); // tangent lines to circle with r = 1.5
            cm.add_bivtp(b);
        }

        cm.set_label("proj. model 3 - lines tangent to circle r = 1.5");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        // lines and reflected lines parallel to x-axis
        auto p1x = pt2dp{-2, 0.5, 1};
        auto p2x = pt2dp{2, 0.5, 1};
        auto p1xr = reflect_on(p1x, x_axis_2dp);
        auto p2xr = reflect_on(p2x, x_axis_2dp);
        pt2d_mark pxm;
        pxm.symbol = Symbol::circle;
        pxm.pen = QPen(Qt::green, 2, Qt::SolidLine);
        cm.add_pt(p1x, pxm);
        cm.add_pt(p2x, pxm);
        pxm.pen = QPen(Qt::cyan, 2, Qt::SolidLine);
        cm.add_pt(p1xr, pxm);
        cm.add_pt(p2xr, pxm);

        auto b12x = wdg(p1x, p2x);
        auto b12xr = reflect_on(b12x, x_axis_2dp);
        cm.add_bivtp(b12x);
        cm.add_bivtp(b12xr);

        // lines and reflected lines
        auto p = pt2dp{-1.5, -1.5, 1};
        auto q = pt2dp{0.5, 1, 1};
        auto b = wdg(p, q);
        auto brx = reflect_on(b, x_axis_2dp); // reflect on x-axis
        auto prx = reflect_on(p, x_axis_2dp); // reflect on x-axis
        auto qrx = reflect_on(q, x_axis_2dp); // reflect on x-axis

        auto bry = reflect_on(b, y_axis_2dp); // reflect on y-axis
        auto pry = reflect_on(p, y_axis_2dp); // reflect on y-axis
        auto qry = reflect_on(q, y_axis_2dp); // reflect on y-axis


        cm.add_pt(p);
        cm.add_pt(q);

        pt2d_mark qrm;
        qrm.pen = QPen(Qt::darkBlue, 2, Qt::SolidLine);
        cm.add_pt(prx, qrm);
        cm.add_pt(qrx, qrm);

        qrm.symbol = Symbol::square;
        qrm.pen = QPen(Qt::darkBlue, 2, Qt::SolidLine);
        cm.add_pt(pry, qrm);
        cm.add_pt(qry, qrm);

        cm.add_bivtp(b);
        cm.add_bivtp(brx);
        cm.add_bivtp(bry);
        cm.add_bivtp(x_axis_2dp);
        cm.add_bivtp(y_axis_2dp);

        // fmt::println("p = {}", p);
        // fmt::println("q = {}", q);
        // fmt::println("qrx = {}", qrx);
        // fmt::println("b = {}, att(b) = {}", b, att(b));
        // fmt::println("brx = {}, att(brx) = {}", brx, att(brx));
        // fmt::println("bry = {}, att(bry) = {}", bry, att(bry));

        // fmt::println("wdg(vec2dp{{0, 0, 1}}, vec2dp{{1, 0, 1}}) = {}, att = {}",
        //              wdg(vec2dp{0, 0, 1}, vec2dp{1, 0, 1}),
        //              att(wdg(vec2dp{0, 0, 1}, vec2dp{1, 0, 1})));
        // fmt::println("wdg(vec2dp{{0, 0, 1}}, vec2dp{{0, 1, 1}}) = {}, att = {}",
        //              wdg(vec2dp{0, 0, 1}, vec2dp{0, 1, 1}),
        //              att(wdg(vec2dp{0, 0, 1}, vec2dp{0, 1, 1})));

        // fmt::println("wdg(vec3dp{{0, 0, 0, 1}}, vec3dp{{1, 0, 0, 1}}) = {}, att = {}",
        //              wdg(vec3dp{0, 0, 0, 1}, vec3dp{1, 0, 0, 1}),
        //              att(wdg(vec3dp{0, 0, 0, 1}, vec3dp{1, 0, 0, 1})));
        // fmt::println("wdg(vec3dp{{0, 0, 0, 1}}, vec3dp{{0, 1, 0, 1}}) = {}, att = {}",
        //              wdg(vec3dp{0, 0, 0, 1}, vec3dp{0, 1, 0, 1}),
        //              att(wdg(vec3dp{0, 0, 0, 1}, vec3dp{0, 1, 0, 1})));
        // fmt::println("wdg(vec3dp{{0, 0, 0, 1}}, vec3dp{{0, 0, 1, 1}}) = {}, att = {}",
        //              wdg(vec3dp{0, 0, 0, 1}, vec3dp{0, 0, 1, 1}),
        //              att(wdg(vec3dp{0, 0, 0, 1}, vec3dp{0, 0, 1, 1})));

        cm.set_label("proj. model 4 - line through p, q & reflected");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p_id = cm.add_apt(pt2d{-1.5, -1});
        size_t q_id = cm.add_apt(pt2d{0.5, 1});

        cm.add_abivtp(abivt2dp{p_id, q_id});

        // fmt::println("p = {}", p);
        // fmt::println("q = {}", q);

        cm.set_label("proj. model 5 - line through active p, q");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2d p0(0, 0);
        pt2d p1(-1, 1);
        pt2d p2(-2, 1);

        cm.add_vt(vt2d{p0, p1});
        cm.add_vt(vt2d{p0, p2});
        cm.add_vt(vt2d{p1, p2});

        cm.set_label("vector model 2");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 1});
        size_t p1a_id = cm.add_apt(pt2d{-1.5, 1});
        size_t p1b_id = cm.add_apt(pt2d{-1.5, 1});
        size_t p2_id = cm.add_apt(pt2d{0, 2});

        cm.add_avt(avt2d{p0_id, p1a_id});
        cm.add_avt(avt2d{p1b_id, p2_id});

        cm.set_label("avts with separate apts");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 1});
        size_t p1_id = cm.add_apt(pt2d{3, 1});
        size_t p2_id = cm.add_apt(pt2d{0, 2});

        cm.add_avt(avt2d{p0_id, p1_id});
        cm.add_avt(avt2d{p1_id, p2_id});

        cm.set_label("avts with common apt");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p1_id = cm.add_apt(pt2d{3, 0});
        size_t p2_id = cm.add_apt(pt2d{0, 2});

        cm.add_abivt(abivt2d{p1_id, p2_id});

        cm.set_label("active bivector");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 0});
        size_t p1_id = cm.add_apt(pt2d{3, 0});
        size_t p2_id = cm.add_apt(pt2d{2, 2});

        cm.add_aproj(aproj2d{p0_id, p1_id, p2_id});

        cm.set_label("active projection");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 1});
        size_t p1_id = cm.add_apt(pt2d{1, 0});

        cm.add_arefl(arefl2d{p0_id, p1_id});

        cm.set_label("active reflection (2x)");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 0});
        size_t p1_id = cm.add_apt(pt2d{1, 0});
        size_t p2_id = cm.add_apt(pt2d{0, 1});

        cm.add_areflp(arefl2dp{p0_id, p1_id, p0_id, p2_id});

        cm.set_label("active reflection projective, common pt (2x)");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 0});
        size_t p1_id = cm.add_apt(pt2d{1, 0});
        size_t p2_id = cm.add_apt(pt2d{0, 0});
        size_t p3_id = cm.add_apt(pt2d{0, 1});

        cm.add_areflp(arefl2dp{p0_id, p1_id, p2_id, p3_id});

        cm.set_label("active reflection projective (2x)");

        vm.push_back(cm);
    }

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
        // add an active point
        if (t == tmin + dt) {
            cm.add_apt(pt2d{1, 2});
        }

        // define a poly line and a model
        ln2d l;
        for (double x = xmin; x <= xmax + x_eps; x += dx) {
            l.emplace_back(pt2d(x, std::sin(omega * t - k * x)));
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

    // remove content from previous scene
    scene->clear();

    // register everything from current model with scene

    // register all vectors
    for (size_t idx = 0; idx < cm->vt.size(); ++idx) {
        scene->addItem(new item_vt2d(cs, wcs, cm, idx));
    }

    // register all projective bivectors (=lines)
    for (size_t idx = 0; idx < cm->bivtp.size(); ++idx) {
        scene->addItem(new item_bivt2dp(cs, wcs, cm, idx));
    }

    // register all lines
    for (size_t idx = 0; idx < cm->ln.size(); ++idx) {
        scene->addItem(new item_ln2d(cs, wcs, cm, idx));
    }

    // register all projective passive lines
    for (size_t idx = 0; idx < cm->clnp.size(); ++idx) {
        scene->addItem(new item_cln2dp(cs, wcs, cm, idx));
    }

    // register all passive points
    for (size_t idx = 0; idx < cm->pt.size(); ++idx) {
        scene->addItem(new item_pt2d(cs, wcs, cm, idx));
    }

    // register all projective passive points
    for (size_t idx = 0; idx < cm->ptp.size(); ++idx) {
        scene->addItem(new item_pt2dp(cs, wcs, cm, idx));
    }

    ///////////////////////////////////////////////////////////////////////////
    // active points
    ///////////////////////////////////////////////////////////////////////////

    // for mapping of idx of active points to users of active points
    std::map<size_t, active_pt2d*> apt2d_map;

    // register all active points
    for (size_t idx = 0; idx < cm->apt.size(); ++idx) {
        // register singular active pts to scene
        active_pt2d* apt = new active_pt2d(cs, wcs, cm->apt[idx]);
        apt2d_map[idx] = apt;
        scene->addItem(apt);
    }

    ///////////////////////////////////////////////////////////////////////////
    // active vectors
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->avt.size(); ++idx) {
        scene->addItem(new active_vt2d(cs, wcs, apt2d_map[cm->avt[idx].beg_idx],
                                       apt2d_map[cm->avt[idx].end_idx]));
    }

    ///////////////////////////////////////////////////////////////////////////
    // active bivectors
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->abivt.size(); ++idx) {
        scene->addItem(new active_bivt2d(cs, wcs, apt2d_map[cm->abivt[idx].uend_idx],
                                         apt2d_map[cm->abivt[idx].vend_idx]));
    }

    ///////////////////////////////////////////////////////////////////////////
    // active projections
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->aproj.size(); ++idx) {
        scene->addItem(new active_projection(cs, wcs, apt2d_map[cm->aproj[idx].beg_idx],
                                             apt2d_map[cm->aproj[idx].uend_idx],
                                             apt2d_map[cm->aproj[idx].vend_idx]));
    }

    ///////////////////////////////////////////////////////////////////////////
    // active reflections
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->arefl.size(); ++idx) {
        scene->addItem(new active_reflection(cs, wcs, apt2d_map[cm->arefl[idx].n1end_idx],
                                             apt2d_map[cm->arefl[idx].n2end_idx]));
    }

    ///////////////////////////////////////////////////////////////////////////
    // active projective bivectors
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->abivtp.size(); ++idx) {
        scene->addItem(new active_bivt2dp(cs, wcs, apt2d_map[cm->abivtp[idx].beg_idx],
                                          apt2d_map[cm->abivtp[idx].end_idx]));
    }

    ///////////////////////////////////////////////////////////////////////////
    // active projective reflections
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->areflp.size(); ++idx) {
        scene->addItem(new active_reflectionp(
            cs, wcs, apt2d_map[cm->areflp[idx].p1_idx], apt2d_map[cm->areflp[idx].p2_idx],
            apt2d_map[cm->areflp[idx].p3_idx], apt2d_map[cm->areflp[idx].p4_idx]));
    }
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
    models = get_model_with_lots_of_stuff();
    // models = get_moving_line();
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
    emit updateMaximumModel(int(vm.size() - 1));

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
