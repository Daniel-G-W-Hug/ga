// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "w_common.hpp"

#include "active_bivt2d.hpp"
#include "active_bivt2dp.hpp"
#include "active_projection.hpp"
#include "active_pt2d.hpp"
#include "active_reflection.hpp"
#include "active_reflectionp.hpp"
#include "active_vt2d.hpp"

bool anyHovered(QList<QGraphicsItem*> const& items)
{

    // return true, if the mouse in hovering on any of the items in the list
    //
    // TODO: for extensions ADD ANY NEW ACTIVE TYPES CAPABALE OF HOVERING HERE!
    //
    bool anyHover{false};
    foreach (QGraphicsItem* item, items) {
        if (item->type() == active_pt2d::Type) // type() == 1
            anyHover = anyHover || qgraphicsitem_cast<active_pt2d*>(item)->isHovered();
        if (item->type() == active_vt2d::Type) // type() == 2
            anyHover = anyHover || qgraphicsitem_cast<active_vt2d*>(item)->isHovered();
        if (item->type() == active_bivt2d::Type) // type() == 3
            anyHover = anyHover || qgraphicsitem_cast<active_bivt2d*>(item)->isHovered();
        if (item->type() == active_projection::Type) // type() == 4
            anyHover =
                anyHover || qgraphicsitem_cast<active_projection*>(item)->isHovered();
        if (item->type() == active_reflection::Type) // type() == 5
            anyHover =
                anyHover || qgraphicsitem_cast<active_reflection*>(item)->isHovered();
        if (item->type() == active_bivt2dp::Type) // type() == 6
            anyHover = anyHover || qgraphicsitem_cast<active_bivt2dp*>(item)->isHovered();
        if (item->type() == active_reflectionp::Type) // type() == 7
            anyHover =
                anyHover || qgraphicsitem_cast<active_reflectionp*>(item)->isHovered();
    }
    return anyHover;
}