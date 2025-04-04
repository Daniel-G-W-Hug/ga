// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "w_common.hpp"

#include "active_bivt2d.hpp"
#include "active_projection.hpp"
#include "active_pt2d.hpp"
#include "active_reflection.hpp"
#include "active_vt2d.hpp"

bool anyHovered(QList<QGraphicsItem*> const& items)
{

    // return true, if the mouse in hovering on any of the items in the list
    //
    // HINT FOR EXTENSION: ADD ANY NEW ACTIVE TYPES CAPABALE OF HOVERING HERE!
    //
    bool anyHover{false};
    foreach (QGraphicsItem* item, items) {
        if (item->type() == active_pt2d::Type)
            anyHover = anyHover || qgraphicsitem_cast<active_pt2d*>(item)->isHovered();
        if (item->type() == active_vt2d::Type)
            anyHover = anyHover || qgraphicsitem_cast<active_vt2d*>(item)->isHovered();
        if (item->type() == active_bivt2d::Type)
            anyHover = anyHover || qgraphicsitem_cast<active_bivt2d*>(item)->isHovered();
        if (item->type() == active_projection::Type)
            anyHover =
                anyHover || qgraphicsitem_cast<active_projection*>(item)->isHovered();
        if (item->type() == active_reflection::Type)
            anyHover =
                anyHover || qgraphicsitem_cast<active_reflection*>(item)->isHovered();
    }
    return anyHover;
}