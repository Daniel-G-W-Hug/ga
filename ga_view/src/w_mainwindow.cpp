// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "w_mainwindow.hpp"
#include "w_coordsys.hpp"
#include "w_statusbar.hpp"


#include "active_bivt2d.hpp"
#include "active_bivt2dp.hpp"
#include "active_kinematics2dp.hpp"
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

        // first rotation by 5 degrees counter-clockwise
        auto mot = get_motor(pc, deg2rad(5));

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

        // second rotation by 10 degrees counter-clockwise
        mot = get_motor(pc, deg2rad(10));

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

        cm.set_label("proj. 0 - rotated lines (not origin)");

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

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        for (int i = 0; i <= 23; ++i) {
            double phi = i * pi / 12 + pi / 2.;
            bivt2dp b(cos(phi), sin(phi), 0.0); // lines through origin
            cm.add_bivtp(b);
        }

        cm.set_label("proj. - lines through orgin");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        for (int i = 0; i <= 23; ++i) {
            double phi = i * pi / 12 + pi / 2.;
            bivt2dp b(cos(phi), sin(phi), 1.5); // tangent lines to circle with r = 1.5
            cm.add_bivtp(b);
        }

        cm.set_label("proj. - lines tangent to circle r = 1.5");

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
        pt2d pi{pip.x, pip.y};
        cm.add_pt(pi);

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
        vt2d v1{pi, p1f}, v2{pi, p2f};

        cm.add_pt(p1f);
        cm.add_pt(p2f);
        cm.add_vt(v1);
        cm.add_vt(v2);

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

        vt2d vres{pi, pfres};
        cm.add_vt(vres);

        cm.add_bivtp(fres, force_mark);

        cm.set_label("proj. - adding force lines (inclined)");

        vm.push_back(cm);
    }

    {
        Coordsys_model cm;

        pt2d p1(0, -0.5); // point p1 (f1 is attached to p1)
        pt2d p2(0, 0.5);  // point p2 (f2 is attached to p2)

        auto p1p = vec2dp{p1.x, p1.y, 1}; // point p1 as projective point
        auto p2p = vec2dp{p2.x, p2.y, 1}; // point p2 as projective point

        double phi = deg2rad(0); // direction of force f1 vs. e1
        double f1_magn = 1.5;    // magnitude of force f1
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
        cm.add_vt(v1);
        cm.add_vt(v2);

        bivt2dp_mark force_mark;

        force_mark.pen = QPen(Qt::green, 2, Qt::SolidLine);
        cm.add_bivtp(f1, force_mark);

        auto v1_beg = vec2d{p1.x, p1.y};
        auto v1_end = vec2d{p1f.x, p1f.y};


        force_mark.pen = QPen(Qt::darkGreen, 2, Qt::SolidLine);
        cm.add_bivtp(f2, force_mark);

        force_mark.pen = QPen(Qt::darkRed, 2, Qt::SolidLine);
        auto fres = f1 + f2;
        // fmt::println("parallel: fres = {}", fres);
        auto psupp = support2dp(fres);
        auto pfresp = psupp + att(fres);
        pt2d psup = pt2d{psupp.x, psupp.y};
        pt2d pfres{pfresp.x, pfresp.y};

        cm.add_pt(psup);
        cm.add_pt(pfres);
        vt2d vres{psup, pfres};
        cm.add_vt(vres);

        cm.add_bivtp(fres, force_mark);

        cm.set_label("proj. - adding force lines (parallel)");

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
        cm.add_vt(v1);
        cm.add_vt(v2);

        bivt2dp_mark force_mark;

        force_mark.pen = QPen(Qt::green, 2, Qt::SolidLine);
        cm.add_bivtp(f1, force_mark);

        auto v1_beg = vec2d{p1.x, p1.y};
        auto v1_end = vec2d{p1f.x, p1f.y};


        force_mark.pen = QPen(Qt::darkGreen, 2, Qt::SolidLine);
        cm.add_bivtp(f2, force_mark);

        force_mark.pen = QPen(Qt::darkRed, 2, Qt::SolidLine);
        auto fres = f1 + f2;
        // fmt::println("anti-parallel: fres = {}", fres);
        auto psupp = support2dp(fres);
        auto pfresp = psupp + att(fres);
        pt2d psup = pt2d{psupp.x, psupp.y};
        pt2d pfres{pfresp.x, pfresp.y};

        cm.add_pt(psup);
        cm.add_pt(pfres);
        vt2d vres{psup, pfres};
        cm.add_vt(vres);

        cm.add_bivtp(fres, force_mark);

        cm.set_label("proj. - adding force lines (anti-parallel)");

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

        cm.set_label("proj. - various reflections");

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
        auto bv_perp = right_weight_expand2dp(p, l); // line perpendicular to l through p
        auto p_perp = rwdg(l, bv_perp);              // intersection point bv_perp and l
        ptm.pen = QPen(Qt::gray, 2, Qt::SolidLine);
        bvm.pen = ptm.pen;
        cm.add_pt(p_perp, ptm);
        cm.add_bivtp(bv_perp, bvm);

        // orthogonal antiprojection of line onto p (i.e. line parallel to l through p)
        auto dir_antiproj = right_weight_contract2dp(l, p); // direction of line
                                                            // (=same attidude as l)
        auto bv_antiproj = wdg(p, dir_antiproj);            // line through p in direction
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

        cm.set_label("proj. - products, refl./ortho. proj.");

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
        auto bv_central = right_bulk_expand2dp(p, l); // proj. line
        auto p_central = rwdg(l, bv_central);         // central proj. point
        ptm.pen = QPen(Qt::cyan, 2, Qt::SolidLine);
        bvm.pen = QPen(Qt::cyan, 1, Qt::SolidLine);
        cm.add_pt(p_central, ptm);
        cm.add_bivtp(bv_central, bvm);

        // central antiprojection of line onto p (i.e. line parallel to l through p)
        auto dir_cproj = right_bulk_contract2dp(l, p); // direction of line
        auto bv_cproj = wdg(p, dir_cproj);             // line through p in direction
        auto p_cproj = rwdg(l, bv_cproj);
        ptm.pen = QPen(Qt::cyan, 2, Qt::SolidLine);
        bvm.pen = QPen(Qt::cyan, 1, Qt::SolidLine);
        cm.add_pt(p_cproj, ptm);
        cm.add_bivtp(bv_cproj, bvm);

        auto lp = right_weight_expand2dp(p_cproj, bv_central);
        bvm.pen = QPen(Qt::gray, 1, Qt::SolidLine);
        cm.add_bivtp(lp, bvm);


        // show coord-axes as bivectors
        bvm.pen = QPen(Qt::black, 2, Qt::SolidLine);
        cm.add_bivtp(x_axis_2dp, bvm);
        cm.add_bivtp(y_axis_2dp, bvm);

        cm.set_label("proj. - products, central proj.");

        vm.push_back(cm);
    }


    {
        Coordsys_model cm;

        size_t p_id = cm.add_apt(pt2d{-1.5, -1});
        size_t q_id = cm.add_apt(pt2d{0.5, 1});

        cm.add_abivtp(abivt2dp{p_id, q_id});

        // fmt::println("p = {}", p);
        // fmt::println("q = {}", q);

        cm.set_label("proj. - join line p ^ q");

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

        cm.set_label("active kinematics");

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

        size_t p0_id = cm.add_apt(pt2d{0, 2.5});
        size_t p1_id = cm.add_apt(pt2d{2.5, 0});

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

    // Set focus to wcs widget so that key presses are received immediately
    wcs->setFocus();
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
