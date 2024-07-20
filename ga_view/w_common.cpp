#include "w_common.hpp"

#include "active_bivec.hpp"
#include "active_projection.hpp"
#include "active_pt.hpp"
#include "active_reflection.hpp"
#include "active_vec.hpp"

bool anyHovered(QList<QGraphicsItem*> const& items)
{

    // return true, if the mouse in hovering on any of the items in the list
    bool anyHover{false};
    foreach (QGraphicsItem* item, items) {
        if (item->type() == active_pt::Type)
            anyHover = anyHover || qgraphicsitem_cast<active_pt*>(item)->isHovered();
        if (item->type() == active_vec::Type)
            anyHover = anyHover || qgraphicsitem_cast<active_vec*>(item)->isHovered();
        if (item->type() == active_bivec::Type)
            anyHover = anyHover || qgraphicsitem_cast<active_bivec*>(item)->isHovered();
        if (item->type() == active_projection::Type)
            anyHover =
                anyHover || qgraphicsitem_cast<active_projection*>(item)->isHovered();
        if (item->type() == active_reflection::Type)
            anyHover =
                anyHover || qgraphicsitem_cast<active_reflection*>(item)->isHovered();
    }
    return anyHover;
}