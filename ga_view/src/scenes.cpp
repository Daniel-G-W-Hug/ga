// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

using namespace hd::ga;
using namespace hd::ga::ega;
using namespace hd::ga::pga;

#include "w_coordsys.hpp"
#include "w_mainwindow.hpp"
#include "w_statusbar.hpp"


#include "active_bivt2d.hpp"
#include "active_bivt2dp.hpp"
#include "active_double_pendulum.hpp"
#include "active_four_bar.hpp"
#include "active_frame_trafo.hpp"
#include "active_grinding_marks.hpp"
#include "active_grinding_topo.hpp"
#include "active_kinematics2dp.hpp"
#include "active_merry_go_round.hpp"
#include "active_ode.hpp"
#include "active_ode_plate.hpp"
#include "active_open_vs_closed.hpp"
#include "active_planar_delta.hpp"
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

#include <cmath>     // std::sqrt
#include <map>       // std::map
#include <stdexcept> // std::runtime_error
#include <utility>   // std::move
#include <vector>    // std::vector

#include "scenes.hpp"

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

        diagram_legend leg;
        leg.heading = "Showing visualization options.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

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

        p0 = pt2d(0, 0);
        p1 = pt2d(1, 1);
        p2 = pt2d(2, 1);

        cm.add_pt(p0);
        cm.add_pt(p1);
        cm.add_pt(p2);

        vt2d v1{p0, p1}, v2{p0, p2}, v3{p1, p2};

        cm.add_vt(v1);
        cm.add_vt(v2);
        cm.add_vt(v3);

        // independent vecs
        size_t p0_id = cm.add_apt(pt2d{-1, -1});
        size_t p1a_id = cm.add_apt(pt2d{-2.5, -1});
        size_t p1b_id = cm.add_apt(pt2d{-2.5, -1});
        size_t p2_id = cm.add_apt(pt2d{-1, -2});

        cm.add_avt(avt2d{p0_id, p1a_id});
        cm.add_avt(avt2d{p1b_id, p2_id});

        // coupled vecs
        p0_id = cm.add_apt(pt2d{1, -1});
        p1a_id = cm.add_apt(pt2d{2.5, -1});
        p2_id = cm.add_apt(pt2d{1, -2});

        cm.add_avt(avt2d{p0_id, p1a_id});
        cm.add_avt(avt2d{p1a_id, p2_id});

        cm.set_label("vector model");

        diagram_legend leg;
        leg.heading = "Drag active points to modify vectors.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }


    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 0});
        size_t p1_id = cm.add_apt(pt2d{3, 0});
        size_t p2_id = cm.add_apt(pt2d{2, 2});

        cm.add_aproj(aproj2d{p0_id, p1_id, p2_id});

        cm.set_label("active projection");

        diagram_legend leg;
        leg.heading = "Drag active points to modify vectors, projection and rejection.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p1_id = cm.add_apt(pt2d{3, 0});
        size_t p2_id = cm.add_apt(pt2d{0, 2});

        cm.add_abivt(abivt2d{p1_id, p2_id});

        cm.set_label("visualized geometric product");

        diagram_legend leg;
        leg.heading = "Key bindings (hover on vectors to activate)";
        leg.entries = {{"U", "rotate vector u"},
                       {"V", "rotate vector v"},
                       {"B", "rotate both u and v"}};
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

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

        cm.set_label("proj. 1 - points");

        diagram_legend leg;
        leg.heading = "Red square shows P(2,1,2) after unitization.\nYields P'(1,0.5,1).";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        for (int i = 0; i <= 23; ++i) {
            double phi = i * pi / 12 + pi / 2.;
            bivt2dp b(cos(phi), sin(phi), 0.0); // lines through origin
            cm.add_bivtp(b);
        }

        cm.set_label("proj. 2 - lines through orgin");

        diagram_legend leg;
        leg.heading = "Lines through origin: \nbivt2dp b(cos(phi), sin(phi), 0.0)";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        for (int i = 0; i <= 23; ++i) {
            double phi = i * pi / 12 + pi / 2.;
            bivt2dp b(cos(phi), sin(phi), 1.5); // tangent lines to circle with r = 1.5
            cm.add_bivtp(b);
        }

        cm.set_label("proj. 3 - lines tangent to circle r = 1.5");

        diagram_legend leg;
        leg.heading =
            "Lines with const. distance to origin: \nbivt2dp b(cos(phi), sin(phi), 1.5)";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p_id = cm.add_apt(pt2d{-1.5, -1});
        size_t q_id = cm.add_apt(pt2d{0.5, 1});

        cm.add_abivtp(abivt2dp{p_id, q_id});

        // fmt::println("p = {}", p);
        // fmt::println("q = {}", q);

        cm.set_label("proj. 4 - join line p ^ q");

        diagram_legend leg;
        leg.heading = "Move active points to drag join line from resulting bivector.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2d p1(1.5, -1); // point p1 (f1 is attached to p1)
        pt2d p2(1, -0.5); // point p2 (f2 is attached to p2)

        auto p1p = vec2dp{p1.x, p1.y, 1}; // point p1 as projective point
        auto p2p = vec2dp{p2.x, p2.y, 1}; // point p2 as projective point

        double phi = deg2rad(15); // direction of force f1 vs. e1
        double f1_magn = 1.5;     // magnitude of force f1
        auto f1_vec = f1_magn * vec2dp{cos(phi), sin(phi), 0.0}; // force vector f1
        auto f1 = wdg(p1p, f1_vec);                              // force f1 acting at p1

        phi = deg2rad(75);    // direction of force f2 vs. e1
        double f2_magn = 3.0; // magnitude of force f2
        auto f2_vec = f2_magn * vec2dp{cos(phi), sin(phi), 0.0}; // force vector f2
        auto f2 = wdg(p2p, f2_vec);                              // force f2 acting at p2

        auto pip = unitize(rwdg(f1, f2)); // intersection point of force lines
        pt2d p_intsec{pip.x, pip.y};
        cm.add_pt(p_intsec);

        // -> variant with force vectors starting at p1, p2
        // cm.add_pt(p1);
        // cm.add_pt(p2);
        // auto p1pf = p1p + f1_vec;
        // auto p2pf = p2p + f2_vec;
        // pt2d p1f{p1pf.x, p1pf.y};
        // pt2d p2f{p2pf.x, p2pf.y};
        // vt2d v1{p1, p1f}, v2{p2, p2f};

        // -> variant with force vectors starting at pi
        auto p1pf = pip + f1_vec;
        auto p2pf = pip + f2_vec;
        pt2d p1f{p1pf.x, p1pf.y};
        pt2d p2f{p2pf.x, p2pf.y};
        vt2d v1(p_intsec, p1f);
        vt2d v2(p_intsec, p2f);

        cm.add_pt(p1f);
        cm.add_pt(p2f);

        vt2d_mark vec_mark;
        vec_mark.display_name = fmt::format("f1 = {:.2f}", f1);
        cm.add_vt(v1, vec_mark);
        vec_mark.display_name = fmt::format("f2 = {:.2f}", f2);
        cm.add_vt(v2, vec_mark);

        bivt2dp_mark force_mark;

        force_mark.pen = QPen(Qt::green, 2, Qt::SolidLine);
        cm.add_bivtp(f1, force_mark);

        auto v1_beg = vec2d{p1.x, p1.y};
        auto v1_end = vec2d{p1f.x, p1f.y};


        force_mark.pen = QPen(Qt::darkGreen, 2, Qt::SolidLine);
        cm.add_bivtp(f2, force_mark);

        force_mark.pen = QPen(Qt::darkRed, 2, Qt::SolidLine);
        auto fres = f1 + f2;
        // fmt::println("inclined: fres = {}", fres);
        auto pfresp = pip + att(fres);
        pt2d pfres{pfresp.x, pfresp.y};
        cm.add_pt(pfres);

        vt2d vres(p_intsec, pfres);
        vec_mark.display_name = fmt::format("fres = {:.2f}", fres);
        cm.add_vt(vres, vec_mark);

        cm.add_bivtp(fres);

        cm.set_label("proj. 5 - adding force lines (inclined)");

        diagram_legend leg;
        leg.heading =
            "Adding force lines (bivectors). \nResult is a force line that contains "
            "resulting line of action, force and moment w.r.t. to origin.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2d p1(0, -0.5); // point p1 (f1 is attached to p1)
        pt2d p2(0, 0.5);  // point p2 (f2 is attached to p2)

        auto p1p = vec2dp{p1.x, p1.y, 1}; // point p1 as projective point
        auto p2p = vec2dp{p2.x, p2.y, 1}; // point p2 as projective point

        double phi = deg2rad(0); // direction of force f1 vs. e1
        double f1_magn = 1.0;    // magnitude of force f1
        auto f1_vec = f1_magn * vec2dp{cos(phi), sin(phi), 0.0}; // force vector f1
        auto f1 = wdg(p1p, f1_vec);                              // force f1 acting at p1

        phi = deg2rad(0);     // direction of force f2 vs. e1
        double f2_magn = 0.5; // magnitude of force f2
        auto f2_vec = f2_magn * vec2dp{cos(phi), sin(phi), 0.0}; // force vector f2
        auto f2 = wdg(p2p, f2_vec);                              // force f2 acting at p2

        // -> variant with force vectors starting at p1, p2
        cm.add_pt(p1);
        cm.add_pt(p2);
        auto p1pf = p1p + f1_vec;
        auto p2pf = p2p + f2_vec;
        pt2d p1f{p1pf.x, p1pf.y};
        pt2d p2f{p2pf.x, p2pf.y};
        vt2d v1{p1, p1f}, v2{p2, p2f};

        cm.add_pt(p1f);
        cm.add_pt(p2f);

        vt2d_mark vec_mark;
        vec_mark.display_name = "f1";
        cm.add_vt(v1, vec_mark);
        vec_mark.display_name = "f2";
        cm.add_vt(v2, vec_mark);

        bivt2dp_mark force_mark;

        force_mark.pen = QPen(Qt::green, 2, Qt::SolidLine);
        force_mark.display_name = fmt::format("bivec1 = {:.2f}", f1);
        cm.add_bivtp(f1, force_mark);

        auto v1_beg = vec2d{p1.x, p1.y};
        auto v1_end = vec2d{p1f.x, p1f.y};


        force_mark.pen = QPen(Qt::darkGreen, 2, Qt::SolidLine);
        force_mark.display_name = fmt::format("bivec2 = {:.2f}", f2);
        cm.add_bivtp(f2, force_mark);


        auto fres = f1 + f2;
        // fmt::println("parallel: f1 = {}", f1);
        // fmt::println("parallel: f2 = {}", f2);
        // fmt::println("parallel: fres = {}", fres);
        auto psupp = unitize(sup(fres)); // unitize is essential or resulting force
                                         // vector might not be located on line
        auto pfresp = psupp + att(fres);
        pt2d psup = pt2d{psupp.x, psupp.y};
        pt2d pfres{pfresp.x, pfresp.y};

        cm.add_pt(psup);
        cm.add_pt(pfres);
        vt2d vres{psup, pfres};
        vec_mark.display_name = "f_res";
        cm.add_vt(vres, vec_mark);

        force_mark.pen = QPen(Qt::darkRed, 2, Qt::SolidLine);
        force_mark.display_name = fmt::format("bivec_res = {:.2f}", fres);
        cm.add_bivtp(fres, force_mark);

        cm.set_label("proj. 6 - adding force lines (parallel)");

        diagram_legend leg;
        leg.heading =
            "Adding force lines (bivectors). \nResult is a force line that contains "
            "resulting line of action, force and moment w.r.t. to origin.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2d p1(0, -1); // point p1 (f1 is attached to p1)
        pt2d p2(0, 1);  // point p2 (f2 is attached to p2)

        auto p1p = vec2dp{p1.x, p1.y, 1}; // point p1 as projective point
        auto p2p = vec2dp{p2.x, p2.y, 1}; // point p2 as projective point

        double phi = deg2rad(180); // direction of force f1 vs. e1
        double f1_magn = 1.5;      // magnitude of force f1
        auto f1_vec = f1_magn * vec2dp{cos(phi), sin(phi), 0.0}; // force vector f1
        auto f1 = wdg(p1p, f1_vec);                              // force f1 acting at p1

        phi = deg2rad(0);     // direction of force f2 vs. e1
        double f2_magn = 1.5; // magnitude of force f2
        auto f2_vec = f2_magn * vec2dp{cos(phi), sin(phi), 0.0}; // force vector f2
        auto f2 = wdg(p2p, f2_vec);                              // force f2 acting at p2

        // -> variant with force vectors starting at p1, p2
        cm.add_pt(p1);
        cm.add_pt(p2);
        auto p1pf = p1p + f1_vec;
        auto p2pf = p2p + f2_vec;
        pt2d p1f{p1pf.x, p1pf.y};
        pt2d p2f{p2pf.x, p2pf.y};
        vt2d v1{p1, p1f}, v2{p2, p2f};

        cm.add_pt(p1f);
        cm.add_pt(p2f);

        vt2d_mark vec_mark;
        vec_mark.display_name = "f1";
        cm.add_vt(v1, vec_mark);
        vec_mark.display_name = "f2";
        cm.add_vt(v2, vec_mark);

        bivt2dp_mark force_mark;

        force_mark.pen = QPen(Qt::green, 2, Qt::SolidLine);
        force_mark.display_name = fmt::format("bivec1 = {:.2f}", f1);
        cm.add_bivtp(f1, force_mark);

        auto v1_beg = vec2d{p1.x, p1.y};
        auto v1_end = vec2d{p1f.x, p1f.y};


        force_mark.pen = QPen(Qt::darkGreen, 2, Qt::SolidLine);
        force_mark.display_name = fmt::format("bivec2 = {:.2f}", f2);
        cm.add_bivtp(f2, force_mark);


        auto fres = f1 + f2;
        // fmt::println("anti-parallel: fres = {}", fres);
        auto psupp = sup(fres);
        auto pfresp = psupp + att(fres);
        pt2d psup = pt2d{psupp.x, psupp.y};
        pt2d pfres{pfresp.x, pfresp.y};

        cm.add_pt(psup);
        cm.add_pt(pfres);
        vt2d vres{psup, pfres};
        vec_mark.display_name = "f_res";
        cm.add_vt(vres, vec_mark);

        force_mark.pen = QPen(Qt::darkRed, 2, Qt::SolidLine);
        force_mark.display_name = fmt::format("bivec_res = {:.2f}", fres);
        cm.add_bivtp(fres, force_mark);

        cm.set_label("proj. 7 - adding force lines (anti-parallel)");

        diagram_legend leg;
        leg.heading =
            "Adding force lines (bivectors). \nResult is a force line that contains "
            "resulting line of action, force and moment w.r.t. to origin.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

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

        bivt2dp_mark bvm;
        bvm.pen = pxm.pen;

        cm.add_pt(p1x, pxm);
        cm.add_pt(p2x, pxm);
        auto b12x = wdg(p1x, p2x);
        cm.add_bivtp(b12x, bvm);

        pxm.pen = QPen(Qt::cyan, 2, Qt::SolidLine);
        cm.add_pt(p1xr, pxm);
        cm.add_pt(p2xr, pxm);
        bvm.pen = pxm.pen;
        auto b12xr = reflect_on(b12x, x_axis_2dp);
        cm.add_bivtp(b12xr, bvm);

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

        // default pt2d_mark (blue, 2pt, solid; plus; size=4)
        pt2d_mark ptm;
        ptm.symbol = Symbol::circle;
        ptm.pen = QPen(Qt::blue, 2, Qt::SolidLine);

        cm.add_pt(p, ptm);
        cm.add_pt(q, ptm);
        bvm.pen = ptm.pen;
        cm.add_bivtp(b, bvm);

        ptm.pen = QPen(Qt::darkBlue, 2, Qt::SolidLine);
        cm.add_pt(prx, ptm);
        cm.add_pt(qrx, ptm);
        bvm.pen = ptm.pen;
        cm.add_bivtp(brx, bvm);

        ptm.pen = QPen(Qt::red, 2, Qt::SolidLine);
        cm.add_pt(pry, ptm);
        cm.add_pt(qry, ptm);
        bvm.pen = ptm.pen;
        cm.add_bivtp(bry, bvm);

        bvm.pen = QPen(Qt::black, 2, Qt::SolidLine);
        cm.add_bivtp(x_axis_2dp, bvm);
        cm.add_bivtp(y_axis_2dp, bvm);

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

        cm.set_label("proj. 8 - various reflections");

        diagram_legend leg;
        leg.heading = "Various reflections on x-axis and y-axis.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        // reference line l and reference point p
        auto l = wdg(pt2dp{-0.5, 1, 1}, pt2dp{1, 1.5, 1});
        auto p = pt2dp{-0.5, 0.5, 1};

        // default pt2d_mark (blue, 2pt, solid; plus; size=4)
        pt2d_mark ptm;
        ptm.symbol = Symbol::circle;
        ptm.pen = QPen(Qt::red, 2, Qt::SolidLine);
        bivt2dp_mark bvm;
        bvm.pen = ptm.pen;
        cm.add_pt(p, ptm);
        cm.add_bivtp(l, bvm);

        // projection perpendiclar to line through p
        auto bv_perp = r_weight_expand2dp(p, l); // line perpendicular to l through p
        auto p_perp = rwdg(l, bv_perp);          // intersection point bv_perp and l
        ptm.pen = QPen(Qt::gray, 2, Qt::SolidLine);
        bvm.pen = ptm.pen;
        cm.add_pt(p_perp, ptm);
        cm.add_bivtp(bv_perp, bvm);

        // orthogonal antiprojection of line onto p (i.e. line parallel to l through p)
        auto dir_antiproj = r_weight_contract2dp(l, p); // direction of line
                                                        // (=same attidude as l)
        auto bv_antiproj = wdg(p, dir_antiproj);        // line through p in direction
        bvm.pen = QPen(Qt::darkBlue, 2, Qt::SolidLine);
        cm.add_bivtp(bv_antiproj, bvm);

        // reflect point p on line l
        auto pr = reflect_on(p, l);
        ptm.pen = QPen(Qt::darkRed, 2, Qt::SolidLine);
        cm.add_pt(pr, ptm);

        // (point-)reflect line l on point p
        auto lpr = invert_on(l, p);
        bvm.pen = ptm.pen;
        cm.add_bivtp(lpr, bvm);

        // show coord-axes as bivectors
        bvm.pen = QPen(Qt::black, 2, Qt::SolidLine);
        cm.add_bivtp(x_axis_2dp, bvm);
        cm.add_bivtp(y_axis_2dp, bvm);

        cm.set_label("proj. 9 - products, reflection, orthogonal proj.");

        diagram_legend leg;
        leg.heading =
            "Reference line l and reference point p (both red) in various operations:";
        leg.entries = {{"gray:", "projection perpendiclar to l through p"},
                       {"dark blue:", "orthogonal antiprojection of l onto p (i.e. "
                                      "line parallel to l through p)"},
                       {"dark red line:", "(point-)reflect line l on point p"},
                       {"dark red point:", "reflect point p on l"}};
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        // lines and reflected lines parallel to x-axis
        auto l = wdg(pt2dp{-0.5, 1, 1}, pt2dp{1, 1.5, 1});
        auto p = pt2dp{-0.5, 0.5, 1};

        // default pt2d_mark (blue, 2pt, solid; plus; size=4)
        pt2d_mark ptm;
        ptm.symbol = Symbol::circle;
        ptm.pen = QPen(Qt::red, 2, Qt::SolidLine);
        cm.add_pt(p, ptm);

        bivt2dp_mark bvm;
        bvm.pen = ptm.pen;
        cm.add_bivtp(l, bvm);

        // central projection through p towards origin
        auto bv_central = r_bulk_expand2dp(p, l); // proj. line
        auto p_central = rwdg(l, bv_central);     // central proj. point
        ptm.pen = QPen(Qt::cyan, 2, Qt::SolidLine);
        bvm.pen = QPen(Qt::cyan, 1, Qt::SolidLine);
        cm.add_pt(p_central, ptm);
        cm.add_bivtp(bv_central, bvm);

        // central antiprojection of line onto p (i.e. line parallel to l through p)
        auto dir_cproj = r_bulk_contract2dp(l, p); // direction of line
        auto bv_cproj = wdg(p, dir_cproj);         // line through p in direction
        auto p_cproj = rwdg(l, bv_cproj);
        ptm.pen = QPen(Qt::blue, 2, Qt::SolidLine);
        bvm.pen = QPen(Qt::blue, 1, Qt::SolidLine);
        cm.add_pt(p_cproj, ptm);
        cm.add_bivtp(bv_cproj, bvm);

        // orthogonal projection through blue point onto cyan line
        auto lp = r_weight_expand2dp(p_cproj, bv_central);
        bvm.pen = QPen(Qt::green, 1, Qt::SolidLine);
        cm.add_bivtp(lp, bvm);


        // show coord-axes as bivectors
        bvm.pen = QPen(Qt::black, 2, Qt::SolidLine);
        cm.add_bivtp(x_axis_2dp, bvm);
        cm.add_bivtp(y_axis_2dp, bvm);

        cm.set_label("proj. 10 - products, central projections");

        diagram_legend leg;
        leg.heading =
            "Reference line l and reference point p (both red) in various operations:";
        leg.entries = {
            {"cyan:", "central projection through p towards origin"},
            {"blue:", "central antiprojection of line onto p (i.e. line "
                      "'projectively parallel' to l through p)"},
            {"green:", "orthogonal projection through blue point onto cyan line"}};
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 2.5});
        size_t p1_id = cm.add_apt(pt2d{2.5, 0});

        cm.add_arefl(arefl2d{p0_id, p1_id});

        cm.set_label("Rotations from reflections (EGA)");

        diagram_legend leg;
        leg.heading = "EGA: Two consequtive reflections yield a rotation. \nEither move "
                      "active points or active line segments between vector tip and "
                      "points (green on hover) by dragging.";
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 0});
        size_t p1_id = cm.add_apt(pt2d{1, 0});
        size_t p2_id = cm.add_apt(pt2d{0, 1});

        cm.add_areflp(arefl2dp{p0_id, p1_id, p0_id, p2_id});

        cm.set_label("Rotations from reflections (PGA)");

        diagram_legend leg;
        leg.heading = "PGA: Two consequtive reflections yield a rotation.\nCommon "
                      "intersection point.\nEither move active points or active lines "
                      "(green on hover) by dragging.";
        leg.entries = {{"A:", "move 1st line"},
                       {"S:", "move 2nd line"},
                       {"D:", "rotate both lines"}};
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        size_t p0_id = cm.add_apt(pt2d{0, 0});
        size_t p1_id = cm.add_apt(pt2d{1, 0});
        size_t p2_id = cm.add_apt(pt2d{0, 0});
        size_t p3_id = cm.add_apt(pt2d{0, 1});

        cm.add_areflp(arefl2dp{p0_id, p1_id, p2_id, p3_id});

        cm.set_label("Rotations from reflections (PGA)");

        diagram_legend leg;
        leg.heading = "PGA: Two consequtive reflections yield a rotation.\nIndependent "
                      "points.\nEither move active points or active lines (green on "
                      "hover) by dragging.";
        leg.entries = {{"A:", "move 1st line"},
                       {"S:", "move 2nd line"},
                       {"D:", "rotate both lines\n(only for intersecting lines)"}};
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2dp pc(-1, -1, 1); // rotation center

        pt2d_mark m;
        m.symbol = Symbol::circle;
        m.pen = QPen(Qt::red, 2, Qt::SolidLine);
        cm.add_pt(pc, m);

        // reference lines
        pt2dp p0(0, 0, 1);
        pt2dp p1(1, 0, 1);
        pt2dp p2(0, 1, 1);

        m.symbol = Symbol::square;
        m.pen = QPen(Qt::magenta, 2, Qt::SolidLine);

        cm.add_pt(p0, m);
        cm.add_pt(p1, m);
        cm.add_pt(p2, m);

        cln2dp l1;
        l1.push_back(p0);
        l1.push_back(p1);
        cln2dp l2;
        l2.push_back(p0);
        l2.push_back(p2);
        cm.add_ln(l1);
        cm.add_ln(l2);

        // first rotation by 7.5 degrees counter-clockwise
        auto mot = get_motor(pc, deg2rad(7.5));

        m.symbol = Symbol::square;
        m.pen = QPen(Qt::green, 2, Qt::SolidLine);

        auto p0r = move2dp(p0, mot);
        auto p1r = move2dp(p1, mot);
        auto p2r = move2dp(p2, mot);

        cm.add_pt(p0r, m);
        cm.add_pt(p1r, m);
        cm.add_pt(p2r, m);

        cln2dp l1r;
        l1r.push_back(p0r);
        l1r.push_back(p1r);
        cln2dp l2r;
        l2r.push_back(p0r);
        l2r.push_back(p2r);
        cm.add_ln(l1r);
        cm.add_ln(l2r);

        // second rotation by 15 degrees counter-clockwise
        mot = get_motor(pc, deg2rad(15));

        m.symbol = Symbol::square;
        m.pen = QPen(Qt::cyan, 2, Qt::SolidLine);

        p0r = move2dp(p0, mot);
        p1r = move2dp(p1, mot);
        p2r = move2dp(p2, mot);

        cm.add_pt(p0r, m);
        cm.add_pt(p1r, m);
        cm.add_pt(p2r, m);

        cln2dp l1rr;
        l1rr.push_back(p0r);
        l1rr.push_back(p1r);
        cln2dp l2rr;
        l2rr.push_back(p0r);
        l2rr.push_back(p2r);
        cm.add_ln(l1rr);
        cm.add_ln(l2rr);

        // after rotation additionally translate the last points by (-2,1,0)
        mot = get_motor(vec2dp{-2, 1, 0});

        m.symbol = Symbol::square;
        m.pen = QPen(Qt::blue, 2, Qt::SolidLine);

        auto p0rt = move2dp(p0r, mot);
        auto p1rt = move2dp(p1r, mot);
        auto p2rt = move2dp(p2r, mot);

        cm.add_pt(p0rt, m);
        cm.add_pt(p1rt, m);
        cm.add_pt(p2rt, m);

        cln2dp l1rt;
        l1rt.push_back(p0rt);
        l1rt.push_back(p1rt);
        cln2dp l2rt;
        l2rt.push_back(p0rt);
        l2rt.push_back(p2rt);
        cm.add_ln(l1rt);
        cm.add_ln(l2rt);

        // first translate the original points by (-2,1,0)
        mot = get_motor(vec2dp{-2, 1, 0});


        auto p0t = move2dp(p0, mot);
        auto p1t = move2dp(p1, mot);
        auto p2t = move2dp(p2, mot);

        // then rotate the translated points by 15°
        mot = get_motor(pc, deg2rad(15));

        auto p0tr = move2dp(p0t, mot);
        auto p1tr = move2dp(p1t, mot);
        auto p2tr = move2dp(p2t, mot);

        m.symbol = Symbol::square;
        m.pen = QPen(Qt::darkRed, 2, Qt::SolidLine);
        cm.add_pt(p0tr, m);
        cm.add_pt(p1tr, m);
        cm.add_pt(p2tr, m);

        m.pen = QPen(Qt::red, 2, Qt::SolidLine);
        cm.add_pt(p0t, m);

        cln2dp l1tr;
        l1tr.push_back(p0tr);
        l1tr.push_back(p1tr);
        cln2dp l2tr;
        l2tr.push_back(p0tr);
        l2tr.push_back(p2tr);

        ln2d_mark l;
        l.pen = QPen(Qt::gray, 2, Qt::DashLine);
        cm.add_ln(l1tr, l);
        cm.add_ln(l2tr, l);

        // now draw the connecting lines
        cln2dp cln_r0;
        cln_r0.push_back(pc);
        cln_r0.push_back(p0);

        cln2dp cln_r0r;
        cln_r0r.push_back(pc);
        cln_r0r.push_back(p0r);

        cln2dp cln_rt;
        cln_rt.push_back(p0r);
        cln_rt.push_back(p0rt);

        l.pen = QPen(Qt::blue, 1, Qt::DashDotLine);
        cm.add_ln(cln_r0, l);
        cm.add_ln(cln_r0r, l);
        cm.add_ln(cln_rt, l);


        cln2dp cln_tr;
        cln_tr.push_back(p0);
        cln_tr.push_back(p0t);

        cln2dp cln_t0t;
        cln_t0t.push_back(pc);
        cln_t0t.push_back(p0t);

        cln2dp cln_t0tr;
        cln_t0tr.push_back(pc);
        cln_t0tr.push_back(p0tr);

        l.pen = QPen(Qt::red, 1, Qt::DashDotLine);
        cm.add_ln(cln_tr, l);
        cm.add_ln(cln_t0t, l);
        cm.add_ln(cln_t0tr, l);

        cm.set_label("Sequence of rotation & translation matters");

        diagram_legend leg;
        leg.heading = "Sequence of rotation and translation matters:";
        leg.entries = {{"blue lines:", "rotate 1st, then translate"},
                       {"red lines:", "translate 1st, then rotate"}};
        leg.x_pct = 0.68;
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        akinematic2dp pt_chain;
        size_t p0_id = cm.add_apt(pt2d{-3.5, 0.5});
        pt_chain.apt.push_back(p0_id);
        size_t p1_id = cm.add_apt(pt2d{-3.0, 0.5});
        pt_chain.apt.push_back(p1_id);
        size_t p2_id = cm.add_apt(pt2d{-2.5, 0.5});
        pt_chain.apt.push_back(p2_id);
        size_t p3_id = cm.add_apt(pt2d{-2.0, 0.5});
        pt_chain.apt.push_back(p3_id);
        size_t p4_id = cm.add_apt(pt2d{-1.5, 0.5});
        pt_chain.apt.push_back(p4_id);
        size_t p5_id = cm.add_apt(pt2d{-1.0, 0.5});
        pt_chain.apt.push_back(p5_id);
        size_t p6_id = cm.add_apt(pt2d{-0.5, 0.5});
        pt_chain.apt.push_back(p6_id);
        size_t p7_id = cm.add_apt(pt2d{0, 0.5});
        pt_chain.apt.push_back(p7_id);
        size_t p8_id = cm.add_apt(pt2d{0.5, 0.5});
        pt_chain.apt.push_back(p8_id);
        size_t p9_id = cm.add_apt(pt2d{1.0, 0.5});
        pt_chain.apt.push_back(p9_id);
        size_t p10_id = cm.add_apt(pt2d{1.5, 0.5});
        pt_chain.apt.push_back(p10_id);
        size_t p11_id = cm.add_apt(pt2d{2.0, 0.5});
        pt_chain.apt.push_back(p11_id);
        size_t p12_id = cm.add_apt(pt2d{2.5, 0.5});
        pt_chain.apt.push_back(p12_id);
        size_t p13_id = cm.add_apt(pt2d{3.0, 0.5});
        pt_chain.apt.push_back(p13_id);
        size_t p14_id = cm.add_apt(pt2d{3.5, 0.5});
        pt_chain.apt.push_back(p14_id);

        cm.add_akinem(pt_chain);

        akinematic2dp pt_chain2;
        size_t p20_id = cm.add_apt(pt2d{-3.5, 2.5});
        pt_chain2.apt.push_back(p20_id);
        size_t p21_id = cm.add_apt(pt2d{-3.0, 2.5});
        pt_chain2.apt.push_back(p21_id);
        size_t p22_id = cm.add_apt(pt2d{-2.5, 2.5});
        pt_chain2.apt.push_back(p22_id);
        size_t p23_id = cm.add_apt(pt2d{-2.0, 2.5});
        pt_chain2.apt.push_back(p23_id);

        cm.add_akinem(pt_chain2);

        cm.set_label("Active kinematics");

        diagram_legend leg;
        leg.heading = "Active kinematic chains. Drag free end point by mouse.";
        leg.x_pct = 0.68;
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        // Create spring-mass ODE system
        spring_params params; // Use default values:
                              // m=2.0, c=0.01, k=40.0, g=9.81, l0=0.1, g=9.81
        // params.c = 1.0;
        // params.l0 = 0.0;
        // params.m = 5.0;

        // Fixation point (user-controllable)
        size_t fix_id = cm.add_apt(pt2d{0.0, 0.0});

        // Create ODE system (mass point is drawn internally, not an active point)
        aode_spring2d spring_system;
        spring_system.fixation_idx = fix_id;
        spring_system.params = params;

        cm.add_aode(spring_system);

        cm.set_label("Spring-mass ODE system");

        diagram_legend leg;
        leg.heading = "Spring-mass system (2nd order ODE).\nDrag the fixation point to "
                      "initiate/change movement.";
        leg.entries = {{"SPACE:", "stop/re-start integration"},
                       {"R:", "reset simulation"}};
        // Defaults (no need to set explicitly, shown here for reference):
        //   leg.x_pct    = 0.02;  // 2 % from left edge of drawing area
        //   leg.y_pct    = 0.02;  // 2 % from top  edge of drawing area
        //   leg.size_pct = 0.30;  // box width = 30 % of drawing area width
        // Example overrides to move / resize the legend:
        //   leg.x_pct    = 0.75;  // near right edge
        //   leg.y_pct    = 0.80;  // near bottom edge
        //   leg.size_pct = 0.20;  // narrower box (20 %)
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        // Rigid flat plate pendulum using PGA2D formulation
        // Plate is fixed at its top-right corner (pivot = user-draggable active point)
        // Gravity creates a torque that drives pendulum-like oscillation
        plate_params params;
        // params.m = 10.0;         // mass [kg]
        // params.w = 2.0;          // width [m]
        // params.h = 2.0;          // height [m]
        // params.g = 9.81;         // gravity [m/s²]
        // params.c = 0.0;          // angular damping (set > 0 to add dissipation)
        // params.phi_init = 0.0;   // initial tilt angle [rad]
        // params.omega_init = 0.0; // initial angular velocity [rad/s]

        // Pivot point (user-draggable, top-right corner of plate)
        size_t pivot_id = cm.add_apt(pt2d{0.0, 0.0});

        aode_plate_pga2dp plate_system;
        plate_system.fixation_idx = pivot_id;
        plate_system.params = params;

        cm.add_aode_plate(plate_system);

        cm.set_label("Plate pendulum (PGA2D)");

        diagram_legend leg;
        leg.heading =
            "Plate pendulum (PGA2D) (2nd order ODE).\nDrag the fixation point to "
            "initiate/change movement.";
        leg.entries = {{"SPACE:", "stop/re-start integration"},
                       {"R:", "reset simulation"},
                       {"-----", "----------"},
                       {"orange:", "gravity at cm (static)"},
                       {"orange:", "static reaction at pivot"},
                       {"orange:", "centrifugal force at cm (dynamic)"},
                       {"orange:", "centrifugal reaction at pivot (dynamic)"},
                       {"orange:", "torque moment (circle area = |torque|)"}};
        leg.size_pct = 0.48; // wider box to fit force-colour descriptions
        cm.set_legend(leg);

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

    ///////////////////////////////////////////////////////////////////////////
    // active projective kinematics chains
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->akinemp.size(); ++idx) {

        std::vector<active_pt2d*> kinematic_chain;
        for (size_t pt_idx = 0; pt_idx < cm->akinemp[idx].apt.size(); ++pt_idx) {
            kinematic_chain.push_back(apt2d_map[cm->akinemp[idx].apt[pt_idx]]);
        }
        scene->addItem(new active_kinematics2dp(cs, wcs, kinematic_chain));
    }

    ///////////////////////////////////////////////////////////////////////////
    // active ODE spring-mass systems
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->aode.size(); ++idx) {
        active_ode* ode_system = new active_ode(
            cs, wcs, apt2d_map[cm->aode[idx].fixation_idx], cm->aode[idx].params);
        // Connect keyboard signals for reset and pause/resume
        QObject::connect(wcs, &w_Coordsys::resetRequested, ode_system,
                         &active_ode::resetSimulation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, ode_system,
                         &active_ode::togglePause);
        scene->addItem(ode_system);
    }

    ///////////////////////////////////////////////////////////////////////////
    // active plate pendulum ODE systems (PGA2D rigid body dynamics)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->aode_plate.size(); ++idx) {
        active_ode_plate* plate_system =
            new active_ode_plate(cs, wcs, apt2d_map[cm->aode_plate[idx].fixation_idx],
                                 cm->aode_plate[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, plate_system,
                         &active_ode_plate::resetSimulation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, plate_system,
                         &active_ode_plate::togglePause);
        scene->addItem(plate_system);
    }

    ///////////////////////////////////////////////////////////////////////////
    // body-frame / world-frame transformation demo items (PGA2D kinematics)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->aft.size(); ++idx) {
        active_frame_trafo* ft = new active_frame_trafo(cs, wcs, cm->aft[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, ft,
                         &active_frame_trafo::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, ft,
                         &active_frame_trafo::togglePause);
        scene->addItem(ft);
    }

    ///////////////////////////////////////////////////////////////////////////
    // merry-go-round demo items (PGA2D kinematic_system2dp)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->amgr.size(); ++idx) {
        active_merry_go_round* mgr =
            new active_merry_go_round(cs, wcs, cm->amgr[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, mgr,
                         &active_merry_go_round::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, mgr,
                         &active_merry_go_round::togglePause);
        QObject::connect(wcs, &w_Coordsys::circlesToggleRequested, mgr,
                         &active_merry_go_round::toggleCircles);
        scene->addItem(mgr);
    }

    ///////////////////////////////////////////////////////////////////////////
    // double-pendulum demo items (PGA2D dynamic_system2dp)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->adp.size(); ++idx) {
        active_double_pendulum* dp =
            new active_double_pendulum(cs, wcs, cm->adp[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, dp,
                         &active_double_pendulum::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, dp,
                         &active_double_pendulum::togglePause);
        QObject::connect(wcs, &w_Coordsys::forcesToggleRequested, dp,
                         &active_double_pendulum::toggleForces);
        QObject::connect(wcs, &w_Coordsys::traceToggleRequested, dp,
                         &active_double_pendulum::toggleTrace);
        QObject::connect(wcs, &w_Coordsys::icrToggleRequested, dp,
                         &active_double_pendulum::toggleICR);
        scene->addItem(dp);
    }

    ///////////////////////////////////////////////////////////////////////////
    // four-bar linkage demo items (PGA2D closed_loop_system2dp)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->afb.size(); ++idx) {
        active_four_bar* fb = new active_four_bar(cs, wcs, cm->afb[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, fb,
                         &active_four_bar::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, fb,
                         &active_four_bar::togglePause);
        QObject::connect(wcs, &w_Coordsys::traceToggleRequested, fb,
                         &active_four_bar::toggleTrace);
        scene->addItem(fb);
    }

    ///////////////////////////////////////////////////////////////////////////
    // open-vs-closed side-by-side demo (PGA2D dynamic_system2dp vs
    // closed_loop_system2dp)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->aovc.size(); ++idx) {
        active_open_vs_closed* ovc =
            new active_open_vs_closed(cs, wcs, cm->aovc[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, ovc,
                         &active_open_vs_closed::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, ovc,
                         &active_open_vs_closed::togglePause);
        QObject::connect(wcs, &w_Coordsys::traceToggleRequested, ovc,
                         &active_open_vs_closed::toggleTrace);
        scene->addItem(ovc);
    }

    ///////////////////////////////////////////////////////////////////////////
    // planar 5-bar parallel manipulator demo items (PGA2D closed_loop_system2dp)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->apld.size(); ++idx) {
        active_planar_delta* pd = new active_planar_delta(cs, wcs, cm->apld[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, pd,
                         &active_planar_delta::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, pd,
                         &active_planar_delta::togglePause);
        QObject::connect(wcs, &w_Coordsys::traceToggleRequested, pd,
                         &active_planar_delta::toggleTrace);
        scene->addItem(pd);
    }

    ///////////////////////////////////////////////////////////////////////////
    // wafer-grinding grain-trajectory demo items (PGA3DP kinematic_system3dp)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->agm.size(); ++idx) {
        active_grinding_marks* gm =
            new active_grinding_marks(cs, wcs, cm->agm[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, gm,
                         &active_grinding_marks::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, gm,
                         &active_grinding_marks::togglePause);
        QObject::connect(wcs, &w_Coordsys::traceToggleRequested, gm,
                         &active_grinding_marks::toggleTrace);
        QObject::connect(wcs, &w_Coordsys::circlesToggleRequested, gm,
                         &active_grinding_marks::cycleRatio);
        scene->addItem(gm);
    }

    ///////////////////////////////////////////////////////////////////////////
    // wafer-grinding surface-topography demo items (PGA3DP kinematic_system3dp)
    ///////////////////////////////////////////////////////////////////////////
    for (size_t idx = 0; idx < cm->agto.size(); ++idx) {
        active_grinding_topo* gt =
            new active_grinding_topo(cs, wcs, cm->agto[idx].params);
        QObject::connect(wcs, &w_Coordsys::resetRequested, gt,
                         &active_grinding_topo::resetAnimation);
        QObject::connect(wcs, &w_Coordsys::pauseToggleRequested, gt,
                         &active_grinding_topo::togglePause);
        QObject::connect(wcs, &w_Coordsys::circlesToggleRequested, gt,
                         &active_grinding_topo::cyclePatch);
        scene->addItem(gt);
    }

    // Set focus to wcs widget so that key presses are received immediately
    wcs->setFocus();
}


// ---------------------------------------------------------------------------
// Body-frame / world-frame transformation scenes (3 kinematic cases)
// ---------------------------------------------------------------------------
//
// Each Coordsys_model holds two active_frame_trafo items that animate in
// parallel: the UPPER item uses M0 = identity, the LOWER item uses M0 ≠ id.
// The user can observe B_b, B_w, and M simultaneously for both sub-cases.
//
// PGA2D encoding conventions used here:
//   Translation to (tx, ty):  B = vec2dp{-ty, tx, 0}
//   Motor:                    M = exp(0.5 * B)
//   Velocity generator:       Omega_b = (-vy, vx, 0) for translation vx,vy
//                             Omega_b = omega * Q_b   for rotation about Q_b
// ---------------------------------------------------------------------------

std::vector<Coordsys_model> get_frame_trafo_scenes()
{
    std::vector<Coordsys_model> vm;

    // Identity motor (reused as M0 for all "M0=id" items)
    mvec2dp_u const M0_id = exp(0.5 * vec2dp{0.0, 0.0, 0.0});

    // Helper lambda: motor for pure translation to world position (tx, ty)
    auto make_M0_trans = [](double tx, double ty) -> mvec2dp_u {
        return exp(0.5 * vec2dp{-ty, tx, 0.0});
    };

    // -----------------------------------------------------------------------
    // Scene 1: Pure translation
    //   Omega_b = (0, 1, 0)  =>  vx = 1 unit/s, vy = 0
    //   6 s animation: body slides 6 units in +x direction.
    //   Key insight: B_w = B_b always for translation, regardless of M0.
    // -----------------------------------------------------------------------
    {
        Coordsys_model cm;

        // Upper half – M0 = identity; W frame at world origin (shared reference).
        // Body slides through y=0: at t=3 s it passes exactly through the W origin,
        // making M0=id self-evident. B_b_init=(0,-3,0) → body starts at (-3, 0).
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Translation;
            p.m0_is_identity = true;
            p.M0 = M0_id;
            p.B_b_init = vec2dp{0.0, -3.0, 0.0}; // body starts at (-3, 0)
            p.Omega_b = vec2dp{0.0, 1.0, 0.0};   // vx = 1 unit/s
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = 0.0; // W marker at world origin
            p.draw_world_frame = true;
            p.text_wy = 1.10; // text 1.1 unit above animation lane (y=0)
            cm.add_aframe_trafo({p});
        }

        // Lower half – M0 = translation to (0, −2); body slides through y=−2.
        // At t=3 s body passes through the orange M₀·O home marker at (0,−2),
        // making the M0 offset immediately visible against the upper lane.
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Translation;
            p.m0_is_identity = false;
            p.M0 = make_M0_trans(0.0, -2.0);     // body home at (0,-2)
            p.B_b_init = vec2dp{0.0, -3.0, 0.0}; // body starts at (-3,-2)
            p.Omega_b = vec2dp{0.0, 1.0, 0.0};   // vx = 1 unit/s
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = -2.0;
            p.draw_world_frame = false; // share the single W at origin
            p.text_wy = -0.9;           // 1.1 unit above lower animation lane (y=-2)
            cm.add_aframe_trafo({p});
        }

        cm.set_label("Frame transformation: pure translation");

        diagram_legend leg;
        leg.heading = "PGA2D: body-frame / world-frame  (pure translation)";
        leg.entries = {{"SPACE:", "pause / resume animation"},
                       {"R:", "reset animation to t=0"},
                       {"─────", "──────────────────────"},
                       {"red/grn (solid):", "world frame W (fixed)"},
                       {"drk-red/grn (dash):", "body frame B (animated)"},
                       {"orange (dot):", "M\u2080\u00B7O \u2014 body home at B_b=0"},
                       {"blue:", "B_b  \u2014 body-frame generator"},
                       {"magenta:", "B_w  \u2014 world-frame generator"},
                       {"─────", "──────────────────────"},
                       {"Key insight:", "B_w = B_b for pure translation"},
                       {"", "translation velocity is frame-indep."}};
        leg.size_pct = 0.36;
        leg.x_pct = 0.62; // upper-right corner (1.0 - size_pct - 0.02)
        leg.y_pct = 0.02;
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    // -----------------------------------------------------------------------
    // Scene 2: Pure rotation
    //   Q_b = (0.8, 0, 1)  (pivot in body frame, distance 0.8 from body origin)
    //   omega = pi rad/s  =>  3 full rotations in 6 s
    //   Omega_b = omega * Q_b = (0.8*pi, 0, pi)
    //   Key insight: B_b.z = B_w.z always (rotation angle is frame-independent);
    //                B_b.xy ≠ B_w.xy when M0 ≠ id (pivot encodes differently).
    // -----------------------------------------------------------------------
    {
        Coordsys_model cm;

        double const omega = M_PI;          // rad/s
        vec2dp const Q_b{0.8, 0.0, 1.0};    // pivot in body frame
        vec2dp const Omega_b = omega * Q_b; // = (0.8π, 0, π)

        // Upper half – M0 = identity, body at (0,0), pivot world at (0.8, 0)
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Rotation;
            p.m0_is_identity = true;
            p.M0 = M0_id;
            p.B_b_init = vec2dp{0.0, 0.0, 0.0}; // body at world origin
            p.Omega_b = Omega_b;
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = 0.0; // W at world origin
            p.draw_world_frame = true;
            p.text_wy = 1.10;
            cm.add_aframe_trafo({p});
        }

        // Lower half – M0 = translation to (0, -2), pivot world at (0.8, -2)
        // Consistent with Scene 1: lower lane 2 units below upper.
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Rotation;
            p.m0_is_identity = false;
            p.M0 = make_M0_trans(0.0, -2.0);
            p.B_b_init = vec2dp{0.0, 0.0, 0.0}; // body at M0 = (0,-2)
            p.Omega_b = Omega_b;
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = -2.0;
            p.draw_world_frame = false; // share single W at world origin
            p.text_wy = -0.90;
            cm.add_aframe_trafo({p});
        }

        cm.set_label("Frame transformation: pure rotation");

        diagram_legend leg;
        leg.heading = "PGA2D: body-frame / world-frame  (pure rotation)";
        leg.entries = {{"SPACE:", "pause / resume animation"},
                       {"R:", "reset animation to t=0"},
                       {"─────", "──────────────────────"},
                       {"red/grn (solid):", "world frame W (fixed)"},
                       {"drk-red/grn (dash):", "body frame B (animated)"},
                       {"orange (dot):", "M\u2080\u00B7O \u2014 body home at B_b=0"},
                       {"blue:", "B_b  \u2014 body-frame generator"},
                       {"magenta:", "B_w  \u2014 world-frame generator"},
                       {"─────", "──────────────────────"},
                       {"Key insight:", "B_b.z = B_w.z always"},
                       {"", "(rotation angle: frame-indep.)"},
                       {"", "B_b.xy \u2260 B_w.xy when M0\u2260id"},
                       {"", "(pivot encodes differently)"}};
        leg.size_pct = 0.36;
        leg.x_pct = 0.62;
        leg.y_pct = 0.02;
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    // -----------------------------------------------------------------------
    // Scene 3: Same Omega_b as Scene 2 (pure rotation, Q_b, omega).
    // Upper: M0 = identity  (reference — identical to Scene 2 upper).
    // Lower: M0 = translation(0,-2) composed with rotation(π/4).
    //   → body home at (0,-2), body frame oriented at 45° to world frame.
    //
    // This isolates the effect of a *combined* M0 on the world-frame generator:
    //   Scene 2 lower: M0 = translation only  → pivot_w shifted but not rotated
    //   Scene 3 lower: M0 = trans + rot       → pivot_w rotated AND shifted
    //                                            ⇒ all B_w components differ from B_b
    // -----------------------------------------------------------------------
    {
        Coordsys_model cm;

        // Same dynamics as Scene 2 — only M0 differs
        double const omega = M_PI;
        vec2dp const Q_b{0.8, 0.0, 1.0};
        vec2dp const Omega_b = omega * Q_b;

        auto make_M0_rot = [](double phi) -> mvec2dp_u {
            return exp(0.5 * vec2dp{0.0, 0.0, phi});
        };

        // Upper half – M0 = identity (reference)
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Combined;
            p.m0_is_identity = true;
            p.M0 = M0_id;
            p.B_b_init = vec2dp{0.0, 0.0, 0.0};
            p.Omega_b = Omega_b;
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = 0.0;
            p.draw_world_frame = true;
            p.text_wy = 1.10;
            cm.add_aframe_trafo({p});
        }

        // Lower half – M0 = translation(0,-2) ⟇ rotation(π/4)
        //   body home at (0,-2), body frame at 45° to world frame
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Combined;
            p.m0_is_identity = false;
            p.M0 = rgpr(make_M0_trans(0.0, -2.0), make_M0_rot(M_PI / 4.0));
            p.B_b_init = vec2dp{0.0, 0.0, 0.0};
            p.Omega_b = Omega_b;
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = -2.0;
            p.draw_world_frame = false;
            p.text_wy = -0.90;
            cm.add_aframe_trafo({p});
        }

        cm.set_label("Frame transformation: M0 = translation + rotation");

        diagram_legend leg;
        leg.heading = "PGA2D: body-frame / world-frame  (M0 combined)";
        leg.entries = {{"SPACE:", "pause / resume animation"},
                       {"R:", "reset animation to t=0"},
                       {"─────", "──────────────────────"},
                       {"red/grn (solid):", "world frame W (fixed)"},
                       {"drk-red/grn (dash):", "body frame B (animated)"},
                       {"orange (dot):", "M\u2080\u00B7O \u2014 body home at B_b=0"},
                       {"blue:", "B_b  \u2014 body-frame generator"},
                       {"magenta:", "B_w  \u2014 world-frame generator"},
                       {"─────", "──────────────────────"},
                       {"Key insight:", "M\u2080 encodes rot+trans:"},
                       {"", "pivot_w = R(\u03C6\u2080)\u00B7Q_b + t\u2080"},
                       {"", "\u21d2 pivot_b \u2260 pivot_w (rotated+shifted)"}};
        leg.size_pct = 0.36;
        leg.x_pct = 0.62;
        leg.y_pct = 0.02;
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    // -----------------------------------------------------------------------
    // Scene 4: World-frame-driven motion — straight-line CM + spin
    //   Upper: pure translation  (-3, +1.5) → (+3, +1.5), vx=1, vy=0
    //   Lower: same translation  (-3, -1.5) → (+3, -1.5) + 2 full rotations
    //          M(t) = T(O(t)) ⟇ R(ω·t)  — analytical, exact straight-line path
    //
    //   Both M0=id.  W marker at world origin (between the two tracks).
    //   Key insight: with constant Omega_w (not Omega_b), the CM moves in a
    //                straight line; body-frame velocity Omega_b(t) varies in time.
    // -----------------------------------------------------------------------
    {
        Coordsys_model cm;

        double const vx_wf = 1.0;               // CM velocity x (world)
        double const vy_wf = 0.0;               // CM velocity y (world)
        double const omega4 = 4.0 * M_PI / 6.0; // 2 full rotations in 6 s

        // Upper half – pure horizontal translation at y = +1.5
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Translation;
            p.m0_is_identity = false;
            p.M0 = make_M0_trans(-3.0, 1.0);
            // B_b_init: body starts at (-3, +1)  →  tx=-3, ty=1  →  B={-1,-3,0}
            p.B_b_init = vec2dp{0.0, 0.0, 0.0};
            p.Omega_b = vec2dp{0.0, 1.0, 0.0}; // vx=1, vy=0
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = 0.0; // W marker at world origin
            p.draw_world_frame = true;
            p.cm_ox = -3.0;
            p.cm_oy = -1.5;
            p.text_wy = 2.30; // top-left text area
            cm.add_aframe_trafo({p});
        }

        // Lower half – translation + 2 full spins at y = −1.5, world-frame driven
        {
            frame_trafo_params p;
            p.kin_case = kin_case_t::Combined;
            p.m0_is_identity = true;
            p.M0 = M0_id;
            // B_b_init initial display value: matches body at (-3, -1.5)
            p.B_b_init = vec2dp{1.5, -3.0, 0.0};
            p.Omega_b = vec2dp{0.0, 0.0, 0.0}; // unused when world_frame_drive
            p.duration = 6.0;
            p.wx = 0.0;
            p.wy = 0.0;
            p.draw_world_frame = false;
            p.world_frame_drive = true;
            p.cm_ox = -3.0;
            p.cm_oy = -1.5;
            p.cm_vx = vx_wf;
            p.cm_vy = vy_wf;
            p.cm_omega = omega4;
            p.text_wy = 0.50; // between tracks, text extends downward
            cm.add_aframe_trafo({p});
        }

        cm.set_label("Frame transformation: world-frame driven (screw)");

        diagram_legend leg;
        leg.heading = "PGA2D: world-frame driven  (straight path + spin)";
        leg.entries = {{"SPACE:", "pause / resume animation"},
                       {"R:", "reset animation to t=0"},
                       {"─────", "──────────────────────"},
                       {"upper:", "pure translation  (Omega_b = const)"},
                       {"lower:", "translation + 2 spins"},
                       {"", "M(t)=T(O(t))\u27C7R(\u03C9t),  exact straight path"},
                       {"─────", "──────────────────────"},
                       {"red/grn (solid):", "world frame W (fixed)"},
                       {"drk-red/grn (dash):", "body frame B (animated)"},
                       {"orange (dot):", "home position at t=0"},
                       {"blue:", "B_b  \u2014 body-frame generator"},
                       {"magenta:", "B_w  \u2014 world-frame generator"},
                       {"─────", "──────────────────────"},
                       {"Key insight:", "const Omega_w \u21d2 straight CM path;"},
                       {"", "Omega_b(t) varies with body angle"}};
        leg.size_pct = 0.40;
        leg.x_pct = 0.59;
        leg.y_pct = 0.02;
        cm.set_legend(leg);

        vm.push_back(cm);
    }

    return vm;
}

// ---------------------------------------------------------------------------
// Wafer-grinding waviness-wavelength panels (Tao 2022, Fig. 12)
//
// Four NATIVE ga_view line charts -- one static Coordsys_model each, with its own axis
// range / labels supplied through Coordsys_model::axis_cfg (the shared Coordsys is
// retargeted on model switch; aspect-ratio lock off). Every plotted point is the
// grinding-waviness wavelength lambda = (relative surface speed) / f_b, with the surface
// speed read off kinematic_system3dp::point_velocity of the Tao Fig.-1 frame tree -- the
// SAME GA twist-field computation validated in the app-test
// "pga3dp: grinding-mark wavelengths lambda_m / lambda_c (Phase D.1)". No active item, no
// hand-drawn axes -- native ga_view line charts.
// ---------------------------------------------------------------------------
namespace {

// {wheel-rim, wafer} surface speeds [um/s] from the GA twist field, for chuck speed
// ns_rpm, wheel speed nw_rpm and wafer radius r_um (mirrors the app-test surface_speeds
// lambda and the Phase-0 Fig.-1 frame tree).
struct gw_speeds {
    double v_w;
    double v_s;
};

gw_speeds gw_surface_speeds(double ns_rpm, double nw_rpm, double r_um, double Rw_um)
{
    double const ns = ns_rpm * pi / 30.0, nw = nw_rpm * pi / 30.0;
    double const off = Rw_um / std::sqrt(2.0);
    double const l3 = 0.3 * Rw_um, x_a = l3;

    kinematic_system3dp s;
    s.add_frame(static_frame3dp("chuck_ctr_stat"));
    s.add_frame(static_frame3dp("chuck_ctr_rot"),
                kin_state3dp{.omega = vec3dp{ns, 0.0, 0.0, 0.0}},
                s.index_of("chuck_ctr_stat"));
    s.add_frame(static_frame3dp("wafer_top_avg_rot", vec3dp{0.0, 0.0, 0.0, 1.0}),
                kin_state3dp{}, s.index_of("chuck_ctr_rot"));
    s.add_frame(static_frame3dp("spindle_cm_stat", vec3dp{x_a, -off, off, 1.0},
                                vec3dp{0.0, -pi / 2.0, 0.0, 0.0}),
                kin_state3dp{}, s.index_of("chuck_ctr_stat"));
    s.add_frame(static_frame3dp("tool_top_avg_rot", vec3dp{0.0, 0.0, l3, 1.0}),
                kin_state3dp{.omega = vec3dp{0.0, 0.0, -nw, 0.0}},
                s.index_of("spindle_cm_stat"));
    s.add_frame(static_frame3dp("tool_surface_avg_at_R", vec3dp{Rw_um, 0.0, 0.0, 1.0}),
                kin_state3dp{}, s.index_of("tool_top_avg_rot"));

    vec3dp const g = unitize(
        move3dp(O_3dp, s.get_pos_trafo("tool_surface_avg_at_R", "chuck_ctr_stat")));
    double const vw = to_val(bulk_nrm(s.point_velocity(g, "tool_surface_avg_at_R")));
    double const vs = to_val(
        bulk_nrm(s.point_velocity(vec3dp{0.0, r_um, 0.0, 1.0}, "wafer_top_avg_rot")));
    return {vw, vs};
}

} // namespace

std::vector<Coordsys_model> get_grinding_wavelength_panels()
{
    double const Rw_um = 150.0e3; // grinding-wheel radius [um]
    // Fixed vibration frequency for the fixed-f_b panels (b)/(d). Tao's Fig.-12 caption
    // states f_b = 4000 Hz, but the plotted curves match f_b = 6000 Hz (lambda = v/f_b;
    // at 4000 Hz lambda_m would top ~11000 um, lambda_c ~314 um, vs. the figure's ~7200 /
    // ~210 um == 6000 Hz). The caption "4000 Hz" is a paper error; we use 6000 Hz to
    // reproduce the figure. Panels (a)/(c) sweep f_b directly and are unaffected.
    double const fb_fixed = 6000.0; // [Hz]
    int const N = 81;

    // four curve families: distinct colour + marker symbol (ga_view has plus/cross/
    // circle/square), matching Tao's per-family markers.
    QColor const col[4] = {QColor(200, 40, 40), QColor(30, 150, 40), QColor(40, 90, 200),
                           QColor(220, 140, 0)};
    Symbol const sym[4] = {Symbol::square, Symbol::circle, Symbol::plus, Symbol::cross};

    auto add_family = [&](Coordsys_model& cm, int fam, double xmin, double xmax,
                          auto&& lambda_of_x) {
        ln2d curve;
        for (int i = 0; i < N; ++i) {
            double const x = xmin + (xmax - xmin) * double(i) / double(N - 1);
            curve.emplace_back(pt2d(x, lambda_of_x(x)));
        }
        ln2d_mark m;
        m.pen = QPen(col[fam], 2, Qt::SolidLine);
        m.mark_pts = true;
        m.delta = 10;
        m.pm.symbol = sym[fam];
        m.pm.pen = QPen(col[fam], 2, Qt::SolidLine);
        m.pm.nsize = 5;
        cm.add_ln(curve, m);
    };

    std::vector<Coordsys_model> vm;

    // --- (a) lambda_c vs f_b, family radius r [mm] -----------------------------------
    {
        Coordsys_model cm;
        double const r_mm[4] = {30.0, 40.0, 50.0, 60.0};
        for (int f = 0; f < 4; ++f)
            add_family(cm, f, 4000.0, 8000.0, [&, f](double fb) {
                return gw_surface_speeds(200.0, 2000.0, r_mm[f] * 1.0e3, Rw_um).v_s / fb;
            });
        cm.set_label("Wafer grinding: lambda_c vs f_b (Tao Fig. 12a)");
        cm.set_axis_cfg(plot_axis_cfg{.xmin = 4000.0,
                                      .xmax = 8000.0,
                                      .ymin = 70.0,
                                      .ymax = 350.0,
                                      .xlabel = "vibration frequency f_b [Hz]",
                                      .ylabel = "lambda_c [um]",
                                      .x_major_delta = 1000.0,
                                      .y_major_delta = 70.0});
        diagram_legend leg;
        leg.heading = "Tao Fig. 12(a): circumferential-waviness wavelength lambda_c = "
                      "v_s/f_b vs vibration frequency. lambda_c falls with f_b and grows "
                      "with radial distance r. v_s from kinematic_system3dp.";
        leg.entries = {{"red []", "r = 30 mm"},
                       {"green ()", "r = 40 mm"},
                       {"blue +", "r = 50 mm"},
                       {"orange x", "r = 60 mm"}};
        leg.x_pct = 0.62;
        leg.y_pct = 0.04;
        leg.size_pct = 0.34;
        cm.set_legend(leg);
        vm.push_back(std::move(cm));
    }

    // --- (b) lambda_c vs wafer speed n_s, family n_w [rpm] (curves COINCIDE) ----------
    {
        Coordsys_model cm;
        double const nw[4] = {1600.0, 2000.0, 2400.0, 2800.0};
        for (int f = 0; f < 4; ++f)
            add_family(cm, f, 100.0, 300.0, [&, f](double ns) {
                return gw_surface_speeds(ns, nw[f], 40.0e3, Rw_um).v_s / fb_fixed;
            });
        cm.set_label("Wafer grinding: lambda_c vs n_s (Tao Fig. 12b)");
        cm.set_axis_cfg(plot_axis_cfg{.xmin = 100.0,
                                      .xmax = 300.0,
                                      .ymin = 70.0,
                                      .ymax = 210.0,
                                      .xlabel = "wafer speed n_s [rpm]",
                                      .ylabel = "lambda_c [um]",
                                      .x_major_delta = 50.0,
                                      .y_major_delta = 35.0});
        diagram_legend leg;
        leg.heading = "Tao Fig. 12(b): lambda_c = v_s/f_b vs wafer speed (f_b = 6000 Hz; "
                      "caption says 4000 -- paper error). lambda_c rises linearly with "
                      "n_s and is INDEPENDENT of wheel speed n_w -> the four curves "
                      "coincide.";
        leg.entries = {{"red []", "n_w = 1600 rpm"},
                       {"green ()", "n_w = 2000 rpm"},
                       {"blue +", "n_w = 2400 rpm"},
                       {"orange x", "n_w = 2800 rpm"}};
        leg.x_pct = 0.04;
        leg.y_pct = 0.04;
        leg.size_pct = 0.36;
        cm.set_legend(leg);
        vm.push_back(std::move(cm));
    }

    // --- (c) lambda_m vs f_b, family radius r [mm] (curves COINCIDE) ------------------
    {
        Coordsys_model cm;
        double const r_mm[4] = {30.0, 40.0, 50.0, 60.0};
        for (int f = 0; f < 4; ++f)
            add_family(cm, f, 4000.0, 8000.0, [&, f](double fb) {
                return gw_surface_speeds(200.0, 2000.0, r_mm[f] * 1.0e3, Rw_um).v_w / fb;
            });
        cm.set_label("Wafer grinding: lambda_m vs f_b (Tao Fig. 12c)");
        cm.set_axis_cfg(plot_axis_cfg{.xmin = 4000.0,
                                      .xmax = 8000.0,
                                      .ymin = 3500.0,
                                      .ymax = 8500.0,
                                      .xlabel = "vibration frequency f_b [Hz]",
                                      .ylabel = "lambda_m [um]",
                                      .x_major_delta = 1000.0,
                                      .y_major_delta = 1000.0});
        diagram_legend leg;
        leg.heading =
            "Tao Fig. 12(c): mark-direction wavelength lambda_m = v_w/f_b vs "
            "vibration frequency. lambda_m falls with f_b and is INDEPENDENT of "
            "radial distance r -> the four curves coincide. v_w from "
            "kinematic_system3dp.";
        leg.entries = {{"red []", "r = 30 mm"},
                       {"green ()", "r = 40 mm"},
                       {"blue +", "r = 50 mm"},
                       {"orange x", "r = 60 mm"}};
        leg.x_pct = 0.62;
        leg.y_pct = 0.04;
        leg.size_pct = 0.34;
        cm.set_legend(leg);
        vm.push_back(std::move(cm));
    }

    // --- (d) lambda_m vs wafer speed n_s, family n_w [rpm] (FLAT, distinct) -----------
    {
        Coordsys_model cm;
        double const nw[4] = {1600.0, 2000.0, 2400.0, 2800.0};
        for (int f = 0; f < 4; ++f)
            add_family(cm, f, 100.0, 300.0, [&, f](double ns) {
                return gw_surface_speeds(ns, nw[f], 40.0e3, Rw_um).v_w / fb_fixed;
            });
        cm.set_label("Wafer grinding: lambda_m vs n_s (Tao Fig. 12d)");
        cm.set_axis_cfg(plot_axis_cfg{.xmin = 100.0,
                                      .xmax = 300.0,
                                      .ymin = 4000.0,
                                      .ymax = 7600.0,
                                      .xlabel = "wafer speed n_s [rpm]",
                                      .ylabel = "lambda_m [um]",
                                      .x_major_delta = 50.0,
                                      .y_major_delta = 600.0});
        diagram_legend leg;
        leg.heading = "Tao Fig. 12(d): lambda_m = v_w/f_b vs wafer speed (f_b = 6000 Hz; "
                      "caption says 4000 -- paper error). lambda_m is INDEPENDENT of n_s "
                      "(flat lines) and rises with wheel speed n_w.";
        leg.entries = {{"red []", "n_w = 1600 rpm"},
                       {"green ()", "n_w = 2000 rpm"},
                       {"blue +", "n_w = 2400 rpm"},
                       {"orange x", "n_w = 2800 rpm"}};
        leg.x_pct = 0.04;
        leg.y_pct = 0.04;
        leg.size_pct = 0.36;
        cm.set_legend(leg);
        vm.push_back(std::move(cm));
    }

    return vm;
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

    coordsys_data cd("Geometric Algebra (GA) demos");

    return new Coordsys(x, y, cd, keep_aspect_ratio::yes);
}

std::vector<Coordsys_model> build_models()
{
    std::vector<Coordsys_model> models = get_model_with_lots_of_stuff();
    // models = get_moving_line();

    // Append body-frame / world-frame transformation scenes
    {
        auto ft_scenes = get_frame_trafo_scenes();
        for (auto& s : ft_scenes)
            models.push_back(std::move(s));
    }

    // Append merry-go-round scene (kinematic_system2dp)
    {
        Coordsys_model cm;
        // Anchor the world frame W at canvas (-2,-2) to move its axes clear of the busy
        // centre. The kinematic tree is purely relative, so re-anchoring W is just a
        // constant draw-time translation. The platform then sits at (3,2) relative to W,
        // which keeps it drawn at the same canvas spot (1,0) it had before: (1,0) - W's
        // (-2,-2) = (3,2).
        merry_go_round_params mgp;
        mgp.world_origin = vec2dp{-2.0, -2.0, 1.0};
        mgp.centre = vec2dp{3.0, 2.0, 1.0};
        cm.add_merry_go_round(amerry_go_round{mgp});
        cm.set_label("Merry-go-round (kinematic_system2dp)");

        diagram_legend leg;
        leg.heading = "PGA2D kinematics: rotating platform P carrying 3 spinning "
                      "turntables T0..T2.\nThe marked point on T0 traces an epicycloid.";
        leg.entries = {{"SPACE:", "pause / resume animation"},
                       {"R:", "reset animation to t=0"},
                       {"C:", "show / hide magnitude circles"},
                       {"─────", "──────────"},
                       {"red / green:", "frame axes e1 / e2"},
                       {"solid axes:", "ground frame W (fixed)"},
                       {"dashed axes:", "platform frame P (rotating)"},
                       {"grey arms:", "platform structure to turntables"},
                       {"orange dot:", "marked point on T0"},
                       {"blue curve:", "world path of the marked point"},
                       {"magenta arrow:", "marked-point velocity v"},
                       {"black arrow:", "marked-point acceleration a"},
                       {"thin circles:", "momentary |v| / |a| magnitude"}};
        leg.size_pct = 0.34;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append merry-go-round dashboard scene: same simulation, but laid out as a mini
    // merry-go-round in the corner plus |r| / |v| / |a| plots vs the platform angle.
    // The heading / key hints / colour key use the shared diagram_legend (top-left,
    // draggable); the mini sits top-right and the plots below, so they do not collide.
    {
        Coordsys_model cm;
        merry_go_round_params mgp;
        mgp.layout = mgr_layout::dashboard;
        cm.add_merry_go_round(amerry_go_round{mgp});
        cm.set_label("Merry-go-round dashboard (magnitudes vs angle)");

        diagram_legend leg;
        leg.heading = "PGA2D kinematics: marked-point magnitudes over one platform "
                      "revolution (theta = 0..2 pi; newest sample on the right).";
        leg.entries = {{"SPACE:", "pause / resume animation"},
                       {"R:", "reset animation to t=0"},
                       {"C:", "show / hide magnitude circles"},
                       {"─────", "──────────"},
                       {"top-right:", "live merry-go-round (mini)"},
                       {"blue:", "|r|  orbit radius about platform centre"},
                       {"magenta:", "|v|  speed"},
                       {"black:", "|a|  acceleration"}};
        leg.x_pct = 0.02;
        leg.y_pct = 0.02;
        leg.size_pct = 0.32;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append double-pendulum scene (dynamic_system2dp) -- the last example: two coupled
    // revolute plates swinging chaotically under gravity, energy conserved.
    {
        Coordsys_model cm;
        cm.add_double_pendulum(adouble_pendulum{});
        cm.set_label("Double pendulum (dynamic_system2dp)");

        diagram_legend leg;
        leg.heading = "PGA2D dynamics: a chaotic double pendulum of two rigid plates "
                      "(coupled revolute joints). Total energy is conserved.";
        leg.entries = {{"SPACE:", "pause / resume animation"},
                       {"R:", "reset animation to t=0"},
                       {"T:", "toggle path trace (default on)"},
                       {"I:", "toggle centre of rotation (default off)"},
                       {"F:", "toggle forces & reactions (default off)"},
                       {"─────", "──────────"},
                       {"red / green:", "body-frame axes e1 / e2"},
                       {"black dot:", "fixed world pivot"},
                       {"brown dot:", "inter-plate hinge"},
                       {"orange dot:", "tracked tip of plate 2"},
                       {"blue curve:", "world path of the tip (T)"},
                       {"grey curve:", "trace of plate-2 ICR (I)"},
                       {"teal cross:", "instantaneous centre of rotation (I)"},
                       {"orange arrows:", "gravity m g at each cm (F)"},
                       {"purple arrows:", "hinge reaction forces R1 / R2 (F)"},
                       {"top-left:", "live KE / PE / E (energy check)"}};
        leg.x_pct = 0.62;
        leg.y_pct = 0.02;
        leg.size_pct = 0.36;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append four-bar linkage scene (closed_loop_system2dp): a motor-driven Grashof
    // crank-rocker. The spanning tree (crank -> coupler, plus the rocker as a second
    // branch) is closed by one point-coincidence; the dependent joints are re-solved each
    // frame so the loop stays closed, and a coupler point traces the classic coupler
    // curve.
    {
        Coordsys_model cm;
        cm.add_four_bar(afour_bar{});
        cm.set_label("Four-bar linkage (closed_loop_system2dp)");

        diagram_legend leg;
        leg.heading =
            "PGA2D closed loop: a motor-driven four-bar (Grashof crank-rocker). "
            "The dependent coupler / rocker angles are re-solved each frame, so "
            "the loop stays closed (live ||g|| top-left).";
        leg.entries = {{"SPACE / R:", "pause-resume / reset"},
                       {"T:", "toggle coupler-curve trace"},
                       {"─────", "──────────"},
                       {"bars:", "crank (blue) / coupler (brown) / rocker (green)"},
                       {"O2, O4:", "fixed ground pivots (grey link between)"},
                       {"B (orange):", "loop-closure pin (coupler tip = rocker tip)"},
                       {"P (blue):", "coupler point tracing the coupler curve"}};
        leg.x_pct = 0.43; // right-aligned (left margin = 1 - size - small gap)
        leg.y_pct = 0.64; // lower part of the view, a touch above the bottom
        leg.size_pct = 0.55;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append open-vs-closed side-by-side scene: the same 3-link arm driven by one shared
    // shoulder oscillation. LEFT (dynamic_system2dp tree) -- free hand sweeps an arc;
    // RIGHT (closed_loop_system2dp) -- hand pinned, so the elbow/wrist re-solve and the
    // arm morphs around the planted hand. The headline "tree vs. tree + closure" visual.
    {
        Coordsys_model cm;
        cm.add_open_vs_closed(aopen_vs_closed{});
        cm.set_label("Open chain vs. closed loop (same arm, one shoulder drive)");

        diagram_legend leg;
        leg.heading =
            "PGA2D tree vs. tree + closure: ONE shared shoulder drive, two "
            "behaviours. Left: open chain, the hand sweeps free. Right: the hand "
            "is pinned, so the elbow / wrist re-solve to keep it planted.";
        leg.entries = {{"SPACE / R:", "pause-resume / reset"},
                       {"T:", "toggle free-hand trace"},
                       {"─────", "──────────"},
                       {"left arm:", "open chain -- elbow/wrist fixed, hand free"},
                       {"blue dot + curve:", "the free hand and its swept arc"},
                       {"right arm:", "closed loop -- same shoulder drive"},
                       {"orange cross:", "fixed pin; the hand stays on it (live ||g||)"}};
        leg.x_pct = 0.20; // centred horizontally (left margin = (1 - size) / 2)
        leg.y_pct = 0.65; // a touch higher than before
        leg.size_pct = 0.60;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append planar 5-bar parallel manipulator scene (closed_loop_system2dp): the 2D
    // analogue of a delta / Stewart-Gough robot. Two 2-link arms rise from two fixed,
    // actuated shoulders; their forearm tips are pinned to a shared end-effector by one
    // point-coincidence. The two shoulders are driven by phase-shifted sinusoids; the
    // dependent elbows are re-solved each frame so the tips stay coincident (forward
    // kinematics of the closed loop) and the effector traces a Lissajous-like curve.
    {
        Coordsys_model cm;
        cm.add_planar_delta(aplanar_delta{});
        cm.set_label("Planar 5-bar parallel manipulator (closed_loop_system2dp)");

        diagram_legend leg;
        leg.heading =
            "PGA2D parallel robot (2D delta analogue): two arms position one effector IN "
            "PARALLEL. The two shoulders are driven; the elbows are re-solved each frame "
            "so the tips stay pinned (live ||g|| top-left), tracing a Lissajous curve.";
        leg.entries = {{"SPACE / R:", "pause-resume / reset"},
                       {"T:", "toggle effector-path trace"},
                       {"─────", "──────────"},
                       {"arms:", "left (blue) / right (green), each shoulder->elbow"},
                       {"SA, SB:", "fixed actuated shoulder pivots (grey base link)"},
                       {"EA, EB:", "elbows (re-solved by the loop closure)"},
                       {"P (orange):", "shared end-effector tracing the curve"}};
        leg.x_pct = 0.43; // right-aligned (left margin = 1 - size - small gap)
        leg.y_pct = 0.64; // lower part of the view, a touch above the bottom
        leg.size_pct = 0.55;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append wafer-grinding grain-trajectory scene (kinematic_system3dp, Phase 0 of
    // TODO/grinding.md): the six-frame Tao Fig.-1 tree, with a wheel-rim grain traced in
    // the rotating wafer frame and projected onto the chuck plane e423_3dp -- the
    // grinding-mark rosette of Fig. 1 / Fig. 7. Normalized to R = wheel = wafer = 1.
    {
        Coordsys_model cm;
        cm.add_grinding_marks(agrinding_marks{});
        cm.set_label(
            "Wafer grinding: grain marks in the wafer frame (kinematic_system3dp)");

        diagram_legend leg;
        leg.heading =
            "PGA3D kinematics (Tao 2022, Fig. 1): a wheel-rim grain seen in the ROTATING "
            "WAFER frame, projected onto the chuck plane e423 (top view down -e1). Marks "
            "are cut only while the grain is on the wafer (r <= R), so they sit still "
            "here "
            "as the wheel orbits. C cycles the n_w/n_s ratio (shown bottom-left).";
        leg.entries = {
            {"SPACE:", "pause / resume animation"},
            {"R:", "reset animation to t=0"},
            {"T:", "show / hide the grinding marks"},
            {"C:", "cycle spin ratio (10.0 / 10.3 / 9.7)"},
            {"─────", "──────────"},
            {"blue disk:", "wafer (fixed in this frame), unit radius R"},
            {"blue spoke:", "wafer orientation mark (shows chuck spin)"},
            {"orange disk:", "grinding wheel (orbits as the wafer turns)"},
            {"orange dot:", "the tracked rim grain"},
            {"orange marks:", "grain strokes cut into the wafer (r <= R only)"},
            {"axes:", "horizontal e2 (root y), vertical e3 (root z)"}};
        leg.x_pct = 0.02;
        leg.y_pct = 0.02;
        leg.size_pct = 0.24;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append the SAME grinding kinematics drawn in the stationary GLOBAL/root frame: the
    // wheel sits at a fixed position and the grain just circles the wheel rim, while the
    // wafer (its orientation spoke) spins beneath. The rosette of the previous scene is
    // the SAME motion seen in the rotating wafer frame -- a direct frame-of-reference
    // contrast.
    {
        Coordsys_model cm;
        grinding_marks_params gmp;
        gmp.view = gm_view::root_frame;
        cm.add_grinding_marks(agrinding_marks{gmp});
        cm.set_label("Wafer grinding: grain path in the global frame (wheel fixed)");

        diagram_legend leg;
        leg.heading =
            "PGA3D kinematics (Tao 2022, Fig. 1): the SAME grinding marks, seen in the "
            "STATIONARY global frame (projected onto e423). The wheel is fixed; the "
            "marked "
            "wafer rotates, carrying its strokes with it, while new strokes are cut at "
            "the "
            "contact region. C cycles the n_w/n_s ratio (shown bottom-left).";
        leg.entries = {
            {"SPACE:", "pause / resume animation"},
            {"R:", "reset animation to t=0"},
            {"T:", "show / hide the grinding marks"},
            {"C:", "cycle spin ratio (10.0 / 10.3 / 9.7)"},
            {"─────", "──────────"},
            {"blue disk:", "wafer (spins; its spoke shows the chuck rotation)"},
            {"orange disk:", "grinding wheel (fixed in the global frame)"},
            {"orange dot:", "the tracked rim grain"},
            {"orange marks:", "strokes cut into the wafer (rotate with it)"},
            {"axes:", "horizontal e2 (root y), vertical e3 (root z)"}};
        leg.x_pct = 0.02;
        leg.y_pct = 0.02;
        leg.size_pct = 0.24;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append the wafer-grinding surface-topography scene (Phase D.1c of TODO/grinding.md;
    // Tao 2022, Figs. 9/10): a heatmap of the ground-wafer height field, whose waviness
    // wavelengths lambda_m / lambda_c are computed live from kinematic_system3dp::
    // point_velocity (the GA twist field). C toggles the WCD (Fig. 9) and WMD (Fig. 10)
    // crops.
    {
        Coordsys_model cm;
        cm.add_grinding_topo(agrinding_topo{});
        cm.set_label("Wafer grinding: ground-surface topography (Tao Figs. 9/10)");

        diagram_legend leg;
        leg.heading =
            "PGA3D surface topography (Tao 2022, Figs. 9/10): the wheel-spindle "
            "vibration "
            "at f_b modulates the grain cut depth, carving a waviness of wavelength "
            "lambda = (surface speed)/f_b. The wheel speed v_w and wafer speed v_s are "
            "read "
            "from kinematic_system3dp::point_velocity, so lambda_m / lambda_c (shown "
            "bottom-left) are built straight from the GA kinematics and match the paper.";
        leg.entries = {{"SPACE:", "pause / resume the feed-scroll"},
                       {"R:", "reset"},
                       {"C:", "toggle WCD (Fig. 9) / WMD (Fig. 10)"},
                       {"─────", "──────────"},
                       {"colour:", "ground-surface height (jet: blue low .. red high)"},
                       {"lambda_c:", "circumferential waviness = v_s / f_b (WCD)"},
                       {"lambda_m:", "mark-direction waviness = v_w / f_b (WMD)"},
                       {"horizontal:", "mark / radial direction"},
                       {"vertical:", "circumferential direction"}};
        leg.x_pct = 0.02;
        leg.y_pct = 0.02;
        leg.size_pct = 0.22;
        cm.set_legend(leg);

        models.push_back(std::move(cm));
    }

    // Append the four wafer-grinding wavelength panels (Tao Fig. 12 a/b/c/d) -- one
    // ga_view line chart each, with its own axis range via Coordsys_model::axis_cfg.
    {
        auto gw_panels = get_grinding_wavelength_panels();
        for (auto& p : gw_panels)
            models.push_back(std::move(p));
    }

    return models;
}
