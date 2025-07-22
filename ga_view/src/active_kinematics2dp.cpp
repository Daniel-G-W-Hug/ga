// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "active_kinematics2dp.hpp"
#include "active_common.hpp"

#include <algorithm> // std::min, std::max


active_kinematics2dp::active_kinematics2dp(Coordsys* cs, w_Coordsys* wcs,
                                           std::vector<active_pt2d*> pts,
                                           QGraphicsItem* parent) :
    QGraphicsItem(parent), cs{cs}, wcs{wcs}, m_pts{pts}
{
    // It is essential to NOT use QGraphicsItem::ItemIsMovable, but keep the other flags
    // they are important for window resizing to work correctly
    setFlags(QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges |
             QGraphicsItem::ItemSendsScenePositionChanges);
    setAcceptHoverEvents(false);

    // setZValue(18);

    // connect active points to relevant signals
    for (size_t i = 0; i < m_pts.size(); ++i) {
        connect(wcs, &w_Coordsys::viewResized, m_pts[i], &active_pt2d::viewChanged);
        connect(m_pts[i], &active_pt2d::pointMoved, this,
                &active_kinematics2dp::pointMoved);
    }
    connect(wcs, &w_Coordsys::viewResized, this, &active_kinematics2dp::viewChanged);

    // TODO: make sure that pts vector is large enough (base + at least 2 points + target)
}

void active_kinematics2dp::paint(QPainter* qp, const QStyleOptionGraphicsItem* option,
                                 QWidget* widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)

    // clipping area is active area of coordsys
    qp->setClipRect(
        mapRectFromScene(QRect(cs->x.nmin(), cs->y.nmax(), cs->x.nmax() - cs->x.nmin(),
                               cs->y.nmin() - cs->y.nmax())));


    // draw in item coordinate system
    qp->save();

    qp->setPen(QPen(QBrush(Qt::black), 2, Qt::SolidLine));
    qp->setBrush(Qt::black);

    // draw arrows between active points (w/o the target point)
    for (size_t i = 0; i < m_pts.size() - 2; ++i) {

        QPointF beg_pos = QPointF(cs->x.au_to_w(m_pts[i]->scenePos().x),
                                  cs->y.au_to_w(m_pts[i]->scenePos().y));
        QPointF end_pos = QPointF(cs->x.au_to_w(m_pts[i + 1]->scenePos().x),
                                  cs->y.au_to_w(m_pts[i + 1]->scenePos().y));

        QPen pen = qp->pen();
        pen.setWidth(2);
        qp->drawPath(arrowLine(beg_pos, end_pos));

        // from here on we want to draw with a small pen to get a pointy vector head
        pen.setWidth(1);
        qp->drawPath(arrowHead(beg_pos, end_pos));
    }

    // draw bounding box (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawRect(boundingRect());

    // draw shape (optional for testing)
    // qp->setPen(col_yel);
    // qp->setBrush(col_yel);
    // qp->drawPath(shape());

    qp->restore();
}

QRectF active_kinematics2dp::boundingRect() const
{
    // initialize with extreme values to be displayed in cs
    double min_x = cs->x.max();
    double max_x = cs->x.min();
    double min_y = cs->y.max();
    double max_y = cs->y.min();

    // find actual extreme values
    for (size_t i = 0; i < m_pts.size(); ++i) {
        min_x = std::min(min_x, m_pts[i]->scenePos().x);
        min_y = std::min(min_y, m_pts[i]->scenePos().y);
        max_x = std::max(max_x, m_pts[i]->scenePos().x);
        max_y = std::max(max_y, m_pts[i]->scenePos().y);
    }

    // give bounding box in item coordinate system
    return QRectF(QPointF(cs->x.au_to_w(min_x) - RADIUS, cs->y.au_to_w(max_y) - RADIUS),
                  QPointF(cs->x.au_to_w(max_x) + RADIUS, cs->y.au_to_w(min_y) + RADIUS));
}

QPainterPath active_kinematics2dp::shape() const
{
    QPainterPath path;
    path.addRect(active_kinematics2dp::boundingRect());
    return path;
}

void active_kinematics2dp::viewChanged()
{
    // qDebug() << "active_kinematics2dp: viewChanged() received.";

    // view changed by external influence, set to position in model
    update();
}

void active_kinematics2dp::pointMoved()
{
    // qDebug() << "active_kinematics2dp: pointMoved() received.";

    // redo the kinematics calculation and redraw

    // Block all signals during solving (avoid too many recursive updates)
    for (auto* pt : m_pts) {
        pt->blockSignals(true);
    }

    solve_kinematic_chain();

    // Re-enable signals - Qt will update automatically
    for (auto* pt : m_pts) {
        pt->blockSignals(false);
    }

    update();
}

void active_kinematics2dp::solve_kinematic_chain()
{
    // qDebug() << "active_kinematics2dp::solve_kinematic_chain().";

    // create unitized projective points from pt2d for GA calculations
    std::vector<vec2dp> pts;
    pts.reserve(m_pts.size());
    for (size_t i = 0; i < m_pts.size(); ++i) {
        pts.push_back(vec2dp(m_pts[i]->scenePos().x, m_pts[i]->scenePos().y, 1.0));
    }

    // indices of important elements
    const size_t target = m_pts.size() - 1; // point representing the target
    const size_t last_elem = target - 1;    // last element of the kinematic chain
    const size_t base_elem = 0;             // base element of the kinematic chain


    // set end of chain to target
    pts[last_elem] = pts[target];

    // go backwards along kinematic chain and restore length of elements
    for (size_t i = last_elem - 1; i > base_elem; --i) {

        // move against current direction between elements
        auto direction = -normalize(pts[i + 1] - pts[i]) * segment_length;
        pts[i] = move2dp(pts[i + 1], motor(direction));
        // fmt::println("backwards: pts[{}] = {}", i, pts[i]);
    }

    // go forwards along kinematic chain and restore length of elements again
    for (size_t i = base_elem + 1; i < target; ++i) {

        // move in current direction between elements
        auto direction = normalize(pts[i] - pts[i - 1]) * segment_length;
        pts[i] = move2dp(pts[i - 1], motor(direction));
        // fmt::println("forwards: pts[{}] = {}", i, pts[i]);
    }

    // write changes back to the visualization (keep base and target stable)
    for (size_t i = 1; i < m_pts.size() - 1; ++i) {
        m_pts[i]->setScenePos(pt2d(pts[i].x, pts[i].y));
    }
}
