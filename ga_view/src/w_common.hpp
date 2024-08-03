#pragma once
//
// author: Daniel Hug, 2024
//

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

bool anyHovered(QList<QGraphicsItem*> const& items);