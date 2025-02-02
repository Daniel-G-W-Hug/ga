#pragma once
//
// author: Daniel Hug, 2024 & 2025
//

#include <QColor>
#include <QDebug>
#include <QPainterPath>
#include <QPointF>

// common operations used by several active items

static constexpr float RADIUS = 7.0;     // active_pt2d radius in pixels
static constexpr float MARGIN = 4.0;     // margin for shape of vector
static constexpr float ARROWSIZE = 12.0; // arrow size of vector

static constexpr float ARCRADIUS = 45.0; // radius of arc to show direction of bivector
static constexpr float ARCMARK = 3.0;    // radius of arc mark to show turning direction
static constexpr float ARCDELTA = 22.0;  // max. delta angle of arc mark (begin & end)

// the inline keyword makes sure that the consts are created only once and not in every
// translation unit
inline QColor const col_blue = QColor(0, 0, 255, 127);
inline QColor const col_lblue = QColor(0, 0, 127, 50);
inline QColor const col_green = QColor(0, 255, 0, 255);
inline QColor const col_lgreen = QColor(0, 127, 0, 50);
inline QColor const col_red = QColor(255, 0, 0, 255);
inline QColor const col_lred = QColor(127, 0, 0, 50);
inline QColor const col_yel = QColor(255, 255, 0, 127);

qreal nrm(QPointF const& u);
qreal dot(QPointF const& u, QPointF const& v);
qreal wdg(QPointF const& u, QPointF const& v);

qreal angle_of_line(QPointF const& beg_pos, QPointF const& end_pos);
qreal angle_between_lines(QPointF const& beg_pos, QPointF const& end_upos,
                          QPointF const& end_vpos);
QPainterPath arrowLine(QPointF const& beg_pos, QPointF const& end_pos);
QPainterPath arrowHead(QPointF const& beg_pos, QPointF const& end_pos);
QPainterPath vectorShape(QPointF const& beg_pos, QPointF const& end_pos);
QPainterPath anglePath(QPointF const& beg_pos, QPointF const& end_upos,
                       QPointF const& end_vpos);
QPainterPath planeLine(QPointF const& beg_pos, QPointF const& end_pos, qreal max_size);
