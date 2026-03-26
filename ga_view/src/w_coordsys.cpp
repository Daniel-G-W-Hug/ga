// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "w_coordsys.hpp"
#include "coordsys_model.hpp"

#include <QCursor>
#include <QFont>
#include <QFontMetrics>
#include <QPainter>
#include <QPalette>
#include <QPen>
#include <QPoint>
#include <QRect>
#include <QString>
#include <QWheelEvent>

#include <algorithm> // for std::min and std::max
#include <cmath>     // for axis scaling (and mathematical functions)


w_Coordsys::w_Coordsys(Coordsys* cs, std::vector<Coordsys_model*>& vm,
                       QGraphicsScene* scene, QWidget* parent) :
    QGraphicsView(parent), cs(cs), scene(scene), vm(vm)
{
    scene->setItemIndexMethod(QGraphicsScene::NoIndex);
    scene->setSceneRect(0, 0, cs->x.widget_size(), cs->y.widget_size());
    setScene(scene);

    setCacheMode(QGraphicsView::CacheBackground);
    // setCacheMode(QGraphicsView::CacheNone);
    setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
    setRenderHint(QPainter::Antialiasing);
    setTransformationAnchor(AnchorUnderMouse);
    setMinimumSize(cs->x.widget_size(), cs->y.widget_size());
    setAlignment(Qt::AlignTop | Qt::AlignLeft);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    // setSizePolicy(QSizePolicy::MinimumExpanding,
    // QSizePolicy::MinimumExpanding);
    updateGeometry();

    // receive Mouse Move Events even when no button is pressed
    // (default is false) required to inform the status bar
    // about the current mouse position
    setMouseTracking(true);
    // Accept KeyPress and KeyRelease Events
    setFocusPolicy(Qt::StrongFocus);
}

void w_Coordsys::resizeEvent(QResizeEvent* event)
{
    QSize oldSize = event->oldSize();
    QSize currentSize = event->size(); // new widget size after resize
    if (oldSize != currentSize) {
        // qDebug() << "oldSize:     " << oldSize;
        // qDebug() << "currentSize: " << currentSize;
        // qDebug() << "cs->x.widget_size(): " << cs->x.widget_size();
        // qDebug() << "cs->y.widget_size(): " << cs->y.widget_size();
        cs->adjust_to_resized_widget(currentSize.width(), currentSize.height());
        // qDebug() << "cs->x.widget_size(): " << cs->x.widget_size();
        // qDebug() << "cs->y.widget_size(): " << cs->y.widget_size();
        scene->setSceneRect(0, 0, cs->x.widget_size(), cs->y.widget_size());

        emit viewResized(); // to update items
    }
}

void w_Coordsys::paintEvent(QPaintEvent* event)
{
    QPainter painter(viewport());


    QGraphicsView::paintEvent(event);
}


// Returns the legend box in viewport pixel coordinates (empty QRect if no legend).
// The anchor is the top-left corner of the box; x_pct / y_pct measure the
// distance from the left / top edge of the drawing area.
// Heights account for word-wrapped text.
QRect w_Coordsys::legendRect() const
{
    if (!cm || !cm->legend.has_value()) return QRect{};

    diagram_legend const& leg = *cm->legend;

    int const area_w  = cs->x.nmax() - cs->x.nmin();
    int const area_h  = cs->y.nmin() - cs->y.nmax(); // nmin > nmax (y-axis inverted)
    int const box_w   = static_cast<int>(leg.size_pct * area_w);
    int const padding = 8;
    int const inner_w = box_w - 2 * padding;    // available text width

    QFontMetrics const fm_h{QFont("Helvetica", 14, QFont::Bold)};
    QFontMetrics const fm_t{QFont("Helvetica", 12, QFont::Normal)};
    int const key_col_w = box_w / 4;
    int const desc_w    = inner_w - key_col_w;

    // Heading height respects word-wrap within inner_w
    int const heading_h = fm_h.boundingRect(0, 0, inner_w, 0,
                              Qt::TextWordWrap | Qt::AlignLeft,
                              QString::fromStdString(leg.heading)).height();
    int box_h = 2 * padding + heading_h + 6;

    if (!leg.entries.empty()) {
        box_h += 9;                          // separator
        box_h += fm_t.height() + 4;         // column header row
        for (key_legend_entry const& entry : leg.entries) {
            int const key_h  = fm_t.boundingRect(0, 0, key_col_w, 0,
                                   Qt::TextWordWrap | Qt::AlignLeft,
                                   QString::fromStdString(entry.key)).height();
            int const desc_h = fm_t.boundingRect(0, 0, desc_w, 0,
                                   Qt::TextWordWrap | Qt::AlignLeft,
                                   QString::fromStdString(entry.description)).height();
            box_h += std::max({fm_t.height(), key_h, desc_h}) + 4;
        }
    }

    // Anchor: top-left corner of the legend box
    int const nx_topleft = cs->x.nmin() + static_cast<int>(leg.x_pct * area_w);
    int const ny_topleft = cs->y.nmax() + static_cast<int>(leg.y_pct * area_h);

    return QRect{nx_topleft, ny_topleft, box_w, box_h};
}


void w_Coordsys::drawBackground(QPainter* qp, const QRectF& rect)
{
    Q_UNUSED(rect);

    cs->draw(qp);
    // cm->draw(qp, cs);

    // fmt::print("w_Coorsys::drawBackground (0,0,{},{})\n", cs->x.widget_size(),
    //            cs->y.widget_size());
}


void w_Coordsys::drawForeground(QPainter* qp, const QRectF& rect)
{
    Q_UNUSED(rect);

    // fmt::print("w_Coorsys::drawForeground\n");

    if (m_leftButton) {

        // fmt::print("w_Coorsys::drawForeground m_leftButton pressed\n");

        qp->save();
        qp->setPen(QPen(Qt::blue, 2, Qt::SolidLine));
        qp->setBrush(QColor(240, 230, 50, 128)); // transparent yellow

        switch (m_mode) {
            case pz_mode::x_and_y:
                // draw zoom rectangle

                qp->drawRect(m_nx_leftPress, m_ny_leftPress, m_nx_hot - m_nx_leftPress,
                             m_ny_hot - m_ny_leftPress);

                break;
            case pz_mode::x_only:
                qp->drawRect(m_nx_leftPress, cs->y.nmax(), m_nx_hot - m_nx_leftPress,
                             cs->y.nmin() - cs->y.nmax());
                break;
            case pz_mode::y_only:
                qp->drawRect(cs->x.nmin(), m_ny_leftPress, cs->x.nmax() - cs->x.nmin(),
                             m_ny_hot - m_ny_leftPress);
                break;
        }
        qp->restore();
    }

    // Draw legend overlay (always on top, fixed in viewport, independent of zoom/pan)
    if (cm && cm->legend.has_value()) {
        QRect const box = legendRect();
        if (box.isValid()) {
            diagram_legend const& leg = *cm->legend;

            qp->save();

            // Semi-transparent rounded background
            qp->setBrush(QColor(248, 248, 248, 210));
            qp->setPen(QPen(QColor(Qt::darkGray), 1, Qt::SolidLine));
            qp->drawRoundedRect(box, 4, 4);

            int const padding = 8;
            int const inner_w = box.width() - 2 * padding;
            int       x       = box.left() + padding;
            int       y       = box.top()  + padding;

            // Heading: Helvetica 14 Bold (matches axis label style), word-wrapped
            qp->setFont(QFont("Helvetica", 14, QFont::Bold));
            QFontMetrics const fm_h{qp->font()};
            qp->setPen(Qt::black);
            QString const heading_str = QString::fromStdString(leg.heading);
            int const heading_h = fm_h.boundingRect(0, 0, inner_w, 0,
                                      Qt::TextWordWrap | Qt::AlignLeft,
                                      heading_str).height();
            qp->drawText(QRect(x, y, inner_w, heading_h),
                         Qt::TextWordWrap | Qt::AlignLeft,
                         heading_str);
            y += heading_h + 6;

            if (!leg.entries.empty()) {
                // Horizontal separator
                qp->setPen(QPen(Qt::gray, 1, Qt::SolidLine));
                qp->drawLine(box.left() + 4, y + 4, box.right() - 4, y + 4);
                y += 9;

                qp->setFont(QFont("Helvetica", 12, QFont::Normal));
                QFontMetrics const fm_t{qp->font()};
                int const key_col_w = box.width() / 4; // "Key" column: ~25% of box
                int const desc_w    = inner_w - key_col_w;

                // Column headers (italic)
                qp->setFont(QFont("Helvetica", 12, QFont::Normal, true /*italic*/));
                qp->setPen(Qt::darkGray);
                qp->drawText(x,             y + fm_t.ascent(), "Key");
                qp->drawText(x + key_col_w, y + fm_t.ascent(), "Function");
                y += fm_t.height() + 4;

                // Entry rows — description wraps if wider than its column
                qp->setFont(QFont("Helvetica", 12, QFont::Normal));
                qp->setPen(Qt::black);
                for (key_legend_entry const& entry : leg.entries) {
                    QString const key_str  = QString::fromStdString(entry.key);
                    QString const desc_str = QString::fromStdString(entry.description);
                    int const key_h  = fm_t.boundingRect(0, 0, key_col_w, 0,
                                           Qt::TextWordWrap | Qt::AlignLeft,
                                           key_str).height();
                    int const desc_h = fm_t.boundingRect(0, 0, desc_w, 0,
                                           Qt::TextWordWrap | Qt::AlignLeft,
                                           desc_str).height();
                    int const row_h = std::max({fm_t.height(), key_h, desc_h}) + 4;
                    qp->drawText(QRect(x,             y, key_col_w, key_h),
                                 Qt::TextWordWrap | Qt::AlignLeft,
                                 key_str);
                    qp->drawText(QRect(x + key_col_w, y, desc_w,    desc_h),
                                 Qt::TextWordWrap | Qt::AlignLeft,
                                 desc_str);
                    y += row_h;
                }
            }

            qp->restore();
        }
    }
}

void w_Coordsys::keyPressEvent(QKeyEvent* event)
{
    // ignore key repetition, just change the mode if required
    if (event->key() == Qt::Key_X && m_mode != pz_mode::x_only) {
        m_mode = pz_mode::x_only;
        // fmt::print("X pressed\n");
        emit modeChanged(m_action, m_mode);
    }
    if (event->key() == Qt::Key_Y && m_mode != pz_mode::y_only) {
        m_mode = pz_mode::y_only;
        // fmt::print("Y pressed\n");
        emit modeChanged(m_action, m_mode);
    }
    if (event->key() == Qt::Key_Z && (event->modifiers() & Qt::ControlModifier)) {
        // call undo function to reinstate last coordsys
        pop_from_history();
    }

    if (event->key() == Qt::Key_A && m_move_mode != move_mode::shift_line12) {
        m_move_mode = move_mode::shift_line12;
        emit moveModeChanged(m_move_mode);
        // fmt::print("A pressed\n");
    }
    if (event->key() == Qt::Key_S && m_move_mode != move_mode::shift_line34) {
        m_move_mode = move_mode::shift_line34;
        emit moveModeChanged(m_move_mode);
        // fmt::print("S pressed\n");
    }
    if (event->key() == Qt::Key_D && m_move_mode != move_mode::rotate_both_lines) {
        m_move_mode = move_mode::rotate_both_lines;
        emit moveModeChanged(m_move_mode);
        // fmt::print("S pressed\n");
    }

    if (event->key() == Qt::Key_U && m_rotate_pt_mode != rotate_pt_mode::rotate_u) {
        m_rotate_pt_mode = rotate_pt_mode::rotate_u;
        emit rotateModeChanged(m_rotate_pt_mode);
        // fmt::print("U pressed\n");
    }
    if (event->key() == Qt::Key_V && m_rotate_pt_mode != rotate_pt_mode::rotate_v) {
        m_rotate_pt_mode = rotate_pt_mode::rotate_v;
        emit rotateModeChanged(m_rotate_pt_mode);
        // fmt::print("V pressed\n");
    }
    if (event->key() == Qt::Key_B && m_rotate_pt_mode != rotate_pt_mode::rotate_both_uv) {
        m_rotate_pt_mode = rotate_pt_mode::rotate_both_uv;
        emit rotateModeChanged(m_rotate_pt_mode);
        // fmt::print("B pressed\n");
    }

    // Reset active simulations (ODE systems)
    if (event->key() == Qt::Key_R && !event->isAutoRepeat()) {
        emit resetRequested();
        // fmt::print("R pressed - reset requested\n");
    }

    // Toggle pause/resume for active simulations
    if (event->key() == Qt::Key_Space && !event->isAutoRepeat()) {
        emit pauseToggleRequested();
        // fmt::print("Space pressed - pause toggle requested\n");
    }
}

void w_Coordsys::keyReleaseEvent(QKeyEvent* event)
{
    if (event->key() == Qt::Key_X) {
        m_mode = pz_mode::x_and_y;
        // fmt::print("X released\n");
        emit modeChanged(m_action, m_mode);
    }
    if (event->key() == Qt::Key_Y) {
        m_mode = pz_mode::x_and_y;
        // fmt::print("Y released\n");
        emit modeChanged(m_action, m_mode);
    }

    if (event->key() == Qt::Key_A) {
        m_move_mode = move_mode::shift_both_lines;
        emit moveModeChanged(m_move_mode);
        // fmt::print("A released\n");
    }
    if (event->key() == Qt::Key_S) {
        m_move_mode = move_mode::shift_both_lines;
        emit moveModeChanged(m_move_mode);
        // fmt::print("S released\n");
    }
    if (event->key() == Qt::Key_D) {
        m_move_mode = move_mode::shift_both_lines;
        emit moveModeChanged(m_move_mode);
        // fmt::print("D released\n");
    }

    if (event->key() == Qt::Key_U) {
        m_rotate_pt_mode = rotate_pt_mode::rotate_none;
        emit rotateModeChanged(m_rotate_pt_mode);
        // fmt::print("U released\n");
    }
    if (event->key() == Qt::Key_V) {
        m_rotate_pt_mode = rotate_pt_mode::rotate_none;
        emit rotateModeChanged(m_rotate_pt_mode);
        // fmt::print("V released\n");
    }
    if (event->key() == Qt::Key_B) {
        m_rotate_pt_mode = rotate_pt_mode::rotate_none;
        emit rotateModeChanged(m_rotate_pt_mode);
        // fmt::print("B released\n");
    }
}

void w_Coordsys::mousePressEvent(QMouseEvent* event)
{
    // Legend drag: takes priority over zoom/pan; check before hot-area logic
    if (event->button() == Qt::LeftButton && cm && cm->legend.has_value()) {
        QRect const lbox = legendRect();
        if (lbox.isValid() && lbox.contains(event->pos())) {
            m_legend_dragging    = true;
            // Store offset from legend's top-left anchor to the click position
            m_legend_drag_offset = QPoint(event->pos().x() - lbox.left(),
                                          event->pos().y() - lbox.top());
            setCursor(QCursor(Qt::SizeAllCursor));
            update();
            return; // do not start a zoom rectangle
        }
    }

    // accept mouse presses only in hot area
    if (m_hot) {

        bool anyHover = anyHovered(items());

        if (event->button() == Qt::LeftButton && !anyHover) {
            // fmt::print("w_Coordsys::mousePressEvent() left button (zoom)\n");

            // zoom
            m_leftButton = true;
            m_action = pz_action::zoom;
            m_nx_leftPress = event->pos().x();
            m_ny_leftPress = event->pos().y();
            emit modeChanged(m_action, m_mode);
        }

        if (event->button() == Qt::RightButton) {
            // fmt::print("w_Coordsys::mousePressEvent() right button (pan)\n");

            // pan
            m_rightButton = true;
            m_action = pz_action::pan;
            setCursor(QCursor(Qt::OpenHandCursor));

            // store undo info (before pan starts, don't store intermediate steps)
            // call update to statusbar only on keyRelease, to not confuse the user
            push_to_history();
            emit modeChanged(m_action, m_mode);
        }
    }

    update();
    QGraphicsView::mousePressEvent(event);
}

void w_Coordsys::mouseReleaseEvent(QMouseEvent* event)
{
    // Release legend drag
    if (event->button() == Qt::LeftButton && m_legend_dragging) {
        m_legend_dragging = false;
        setCursor(m_hot ? QCursor(Qt::CrossCursor) : QCursor());
        update();
        return;
    }

    // end of zoom event triggered by release of left mouse button
    if (event->button() == Qt::LeftButton && m_leftButton) {
        // fmt::print("w_Coordsys::mouseReleaseEvent() left button\n");

        // reset action, if active action was zoom
        if (m_action == pz_action::zoom) {
            m_leftButton = false;
            m_action = pz_action::none;
            emit modeChanged(m_action, m_mode);
        }

        // return scaled values as new min and max
        double new_xmin = std::min(cs->x.w_to_a(m_nx_leftPress), cs->x.w_to_a(m_nx_hot));
        double new_xmax = std::max(cs->x.w_to_a(m_nx_leftPress), cs->x.w_to_a(m_nx_hot));
        double new_ymin = std::min(cs->y.w_to_a(m_ny_leftPress), cs->y.w_to_a(m_ny_hot));
        double new_ymax = std::max(cs->y.w_to_a(m_ny_leftPress), cs->y.w_to_a(m_ny_hot));

        // fmt::print("xmin={}, xmax={}, ymin={}, ymax={}\n", new_xmin, new_xmax,
        //            new_ymin, new_ymax);

        // only adjust if remaining new x- and y-axis are larger than zero
        // otherwise ignore zoom request
        if (new_xmin != new_xmax && new_ymin != new_ymax) {

            // fmt::print("w_Coordsys::mouseReleaseEvent() left button with change\n\n ");

            // store undo info (only if zoom actually is performed)
            // and update statusbar emediately
            push_to_history();
            emit undoChanged(cs_history.size()); // update undo info in status bar

            // adjust and update
            switch (m_mode) {
                case pz_mode::x_and_y:
                    cs->adjust_to_zoom(new_xmin, new_xmax, new_ymin, new_ymax);
                    break;
                case pz_mode::x_only:
                    cs->adjust_to_zoom(new_xmin, new_xmax, cs->y.min(), cs->y.max());
                    break;
                case pz_mode::y_only:
                    cs->adjust_to_zoom(cs->x.min(), cs->x.max(), new_ymin, new_ymax);
                    break;
            }
            invalidateScene();
            emit viewResized(); // to update items
            updateSceneRect(scene->sceneRect());
        }
    }

    // end of pan event triggered by release of right mouse button
    if (event->button() == Qt::RightButton && m_rightButton) {
        // fmt::print("w_Coordsys::mouseReleaseEvent() right button\n");

        // reset action, if active action was pan
        if (m_action == pz_action::pan) {

            // fmt::print("w_Coordsys::mouseReleaseEvent() end of pan.\n\n");

            // update statusbar (cs_history.size() was changed on keyPressEvent)
            emit undoChanged(cs_history.size());

            m_rightButton = false;
            m_action = pz_action::none;
            emit modeChanged(m_action, m_mode);
        }

        if (m_hot) {
            setCursor(QCursor(Qt::CrossCursor));
        }
        else {
            setCursor(QCursor());
        }
        invalidateScene();
        emit viewResized(); // to update items
        updateSceneRect(scene->sceneRect());
    }

    update();
    QGraphicsView::mouseReleaseEvent(event);
}

void w_Coordsys::mouseMoveEvent(QMouseEvent* event)
{
    // Handle legend drag first; bypasses all normal pan/zoom processing
    if (m_legend_dragging && cm && cm->legend.has_value()) {
        int const nx          = event->pos().x();
        int const ny          = event->pos().y();
        int const nx_topleft  = nx - m_legend_drag_offset.x();
        int const ny_topleft  = ny - m_legend_drag_offset.y();
        int const area_w      = cs->x.nmax() - cs->x.nmin();
        int const area_h      = cs->y.nmin() - cs->y.nmax();
        if (area_w > 0 && area_h > 0) {
            cm->legend->x_pct = std::clamp(
                static_cast<double>(nx_topleft - cs->x.nmin()) / area_w, 0.0, 0.98);
            cm->legend->y_pct = std::clamp(
                static_cast<double>(ny_topleft - cs->y.nmax()) / area_h, 0.0, 0.98);
        }
        m_nx = nx;
        m_ny = ny;
        scene->update(); // triggers drawForeground repaint (same as zoom-rect update)
        return;
    }

    // current mouse position in widget
    int nx = event->pos().x();
    int ny = event->pos().y();

    if (nx != m_nx || ny != m_ny) {
        // mouse moved to a new postion

        // convert coordinates to scaled values
        double x_pos = cs->x.w_to_a(nx);
        double y_pos = cs->y.w_to_a(ny);

        // determine if mouse is in active cs area
        bool hot = false;
        if (x_pos >= cs->x.min() && x_pos <= cs->x.max() && y_pos >= cs->y.min() &&
            y_pos <= cs->y.max()) {
            hot = true;
        }
        m_hot = hot; // store whether mouse is still in hot area

        // send current mouse position to status bar
        mouse_pos_t mouse_pos{nx, ny, x_pos, y_pos};
        emit mouseMoved(hot, mouse_pos);

        // current mouse position in hot area (needed for zoom rectangle)
        if (nx < cs->x.nmin()) {
            m_nx_hot = cs->x.nmin();
        }
        else if (nx > cs->x.nmax()) {
            m_nx_hot = cs->x.nmax();
        }
        else {
            m_nx_hot = nx;
        }

        // different coordinate growth direction
        if (ny > cs->y.nmin()) {
            m_ny_hot = cs->y.nmin();
        }
        else if (ny < cs->y.nmax()) {
            m_ny_hot = cs->y.nmax();
        }
        else {
            m_ny_hot = ny;
        }

        // switch to crosshair cursor in hot area if right button is not pressed
        if (m_hot && cursor() != Qt::CrossCursor && m_rightButton == false) {
            setCursor(QCursor(Qt::CrossCursor));
        }
        // switch to open hand cursor in hot area if right button is pressed
        if (m_hot && cursor() != Qt::OpenHandCursor && m_rightButton == true) {
            setCursor(QCursor(Qt::OpenHandCursor));
        }

        // switch back to default cursor outside of hot area
        if (!m_hot) {
            setCursor(QCursor());
        }

        // Override cursor when hovering over the legend box (signals it is draggable)
        if (cm && cm->legend.has_value()) {
            QRect const lbox = legendRect();
            if (lbox.isValid() && lbox.contains(QPoint(nx, ny))) {
                setCursor(QCursor(Qt::SizeAllCursor));
            }
        }

        // pan (only in hot area)
        if (m_rightButton && m_hot) {
            double dx = x_pos - cs->x.w_to_a(m_nx);
            double dy = y_pos - cs->y.w_to_a(m_ny);
            // fmt::print("dx = {}, dy = {}\n", dx, dy);
            switch (m_mode) {
                case pz_mode::x_and_y:
                    cs->adjust_to_pan(dx, dy);
                    break;
                case pz_mode::x_only:
                    cs->adjust_to_pan(dx, 0.0);
                    break;
                case pz_mode::y_only:
                    cs->adjust_to_pan(0.0, dy);
                    break;
            }
            invalidateScene();
            emit viewResized(); // to update items
            updateSceneRect(scene->sceneRect());
        }

        // store current position
        m_nx = nx;
        m_ny = ny;

        // zoom (initiate update for drawing the zoom frame)
        // allow update outside of hot area as well to include boundaries of hot
        // area easily
        if (m_leftButton) {
            scene->update();
        }

        // fmt::print("m_nx = {}, m_ny = {}\n", m_nx, m_ny);
        // fmt::print("nx = {}, ny = {}\n", nx, ny);
        // fmt::print("delta_nx = {}, delta_ny = {}\n", delta_nx, delta_ny);

        // fmt::print("hot = {}, left = {}, right = {}\n", m_hot, m_leftButton,
        //            m_rightButton);

        // fmt::print("nx = {}, ny = {}, x_pos = {}, y_pos = {}\n", nx, ny, x_pos,
        //            y_pos);
    }

    update();
    QGraphicsView::mouseMoveEvent(event);
}


void w_Coordsys::wheelEvent(QWheelEvent* event)
{
    // qDebug() << "w_coordsys.cpp: wheelEvent";

    // one tick corresponds to 1/8°
    // scroll wheel ticks towards user are > 0; ticks away from user are negative
    // on MacOS two finger movement without pressing is also registered as
    // wheelEvent

    if (event->isEndEvent()) // this currently works reliably only on Mac platform
    {
        // qDebug() << "w_Coordsys::wheelEvent->isEndEvent()";

        // reset action, if active action was wheel_zoom
        if (m_action == pz_action::wheel_zoom) {
            m_action = pz_action::none;
            emit modeChanged(m_action, m_mode);
        }
        return;
    }

    // just use the y-direction ticks as indication for requested zoom factor
    int numTicks = event->angleDelta().y();

    if (m_hot && (numTicks != 0)) {

        if (m_action == pz_action::none) {
            // for new scrollwheel moves, store undo information

            // qDebug() << "w_Coordsys::wheelEvent Begin.";

            // store undo info only if wheel zoom actually is started
            // and update statusbar immediately
            m_action = pz_action::wheel_zoom;
            emit modeChanged(m_action, m_mode);

            push_to_history();
            emit undoChanged(cs_history.size()); // update undo info in status bar
        }

        // fmt::print("w_Coordsys::wheelEvent moving.\n");

        // numTicks is used for scaling: positive values for zoom in, negative for
        // zoom out center point for the scaling is the mouse cursor position
        // m_nx_hot, m_ny_hot
        // the distance between the x- and y-values at the current mouse position
        // and the current min- and max-positions are scaled

        // fmt::print(
        //     "Received a wheelEvent at m_nx_hot = {}, m_ny_hot = {}, m_hot = {}, "
        //     "numTicks = {}\n",
        //     m_nx_hot, m_ny_hot, m_hot, numTicks);

        double x = cs->x.w_to_a(m_nx_hot);
        double y = cs->y.w_to_a(m_ny_hot);

        double dx_min = x - cs->x.min();
        double dx_max = cs->x.max() - x;

        double dy_min = y - cs->y.min();
        double dy_max = cs->y.max() - y;

        double scale_fact = 1.0 - 0.01 * 0.25 * numTicks; // 4 numTicks = 1% scaling

        double new_xmin = x - scale_fact * dx_min;
        double new_xmax = x + scale_fact * dx_max;

        double new_ymin = y - scale_fact * dy_min;
        double new_ymax = y + scale_fact * dy_max;

        // fmt::print("new_xmin={}, new_xmax={}, new_ymin={}, new_ymax={}\n\n",
        //            new_xmin, new_xmax, new_ymin, new_ymax);

        // fmt::print("w_coordsys::wheelEvent: xtarget_ratio = {}, ytarget_ratio =
        // {}\n\n",
        //            cs->get_xtarget_ratio(), cs->get_ytarget_ratio());

        // adjust and update
        switch (m_mode) {
            case pz_mode::x_and_y:
                cs->adjust_to_wheel_zoom(new_xmin, new_xmax, new_ymin, new_ymax,
                                         cs->get_xtarget_ratio(),
                                         cs->get_ytarget_ratio());
                break;
            case pz_mode::x_only:
                cs->adjust_to_wheel_zoom(new_xmin, new_xmax, cs->y.min(), cs->y.max(),
                                         cs->get_xtarget_ratio(),
                                         cs->get_ytarget_ratio());
                break;
            case pz_mode::y_only:
                cs->adjust_to_wheel_zoom(cs->x.min(), cs->x.max(), new_ymin, new_ymax,
                                         cs->get_xtarget_ratio(),
                                         cs->get_ytarget_ratio());
                break;
        }
        invalidateScene();
        emit viewResized(); // to update items
        updateSceneRect(scene->sceneRect());
    }

    update();
    QGraphicsView::wheelEvent(event);
}

void w_Coordsys::push_to_history()
{
    // push a copy of the current cs into cs_history
    cs_history.emplace_back(Coordsys(cs->x, cs->y, cs->get_coordsys_data()));
}

void w_Coordsys::pop_from_history()
{
    // undo function to restore older version of cs
    if (cs_history.size() > 0) {
        // fmt::print("got undo signal: cs_history.size()={}\n",
        // cs_history.size());
        *cs = cs_history[cs_history.size() - 1];
        cs_history.pop_back();
        // just in case the widet was resized => adjust to current size
        cs->adjust_to_resized_widget(width(), height());
        emit undoChanged(cs_history.size()); // update undo info in status bar
        invalidateScene();
        emit viewResized(); // to update items
        updateSceneRect(scene->sceneRect());
    }
}

void w_Coordsys::switch_to_model(size_t idx)
{
    if (idx < vm.size()) {
        // fmt::println("w_Coordsys got signal switch_to_model with value {}", idx);
        cm = vm[idx];
        invalidateScene();
        updateSceneRect(scene->sceneRect());
    }
}