#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <QGraphicsItem>
#include <QGraphicsView>
#include <QList>

// this is for enums used for communication across widgets
//
// this is functions that create dependencies between otherwise independent widgets

// pan, zoom and wheel_zoom actions
enum class pz_action { none, pan, zoom, wheel_zoom };

// mode restriction for pan and zoom handling
// x_and_y: no restriction
// x_only:  restrict pan/zoom to x axis
// y_only:  restrict pan/zoom to y axis
enum class pz_mode { x_and_y, x_only, y_only };

// mouse move mode for active_reflectionp
// shifts relative to mouse move
// rotation around intersection point of lines
enum class move_mode { shift_both_lines, shift_line12, shift_line34, rotate_both_lines };

// mouse move mode for active_bivt2d
// rotate endpoint(s) of vectors instead of translation
enum class rotate_pt_mode { rotate_none, rotate_u, rotate_v, rotate_both_uv };

bool anyHovered(QList<QGraphicsItem*> const& items);