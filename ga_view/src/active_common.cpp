//
// author: Daniel Hug, 2024 & 2025
//

#include "active_common.hpp"

#include <cmath>
#include <numbers>

using std::atan2;
using std::cos;
using std::sin;
using std::sqrt;

qreal nrm(QPointF const& u) { return sqrt(u.x() * u.x() + u.y() * u.y()); }
qreal dot(QPointF const& u, QPointF const& v) { return u.x() * v.x() + u.y() * v.y(); }
qreal wdg(QPointF const& u, QPointF const& v) { return u.x() * v.y() - u.y() * v.x(); }

qreal angle_of_line(QPointF const& beg_pos, QPointF const& end_pos)
{
    // return angle of line between beg_pos and end_pos
    // in left-handed qt device coordinate system
    qreal dy = (end_pos - beg_pos).y();
    qreal dx = (end_pos - beg_pos).x();
    return atan2(-dy, dx);
}

qreal angle_between_lines(QPointF const& beg_pos, QPointF const& end_upos,
                          QPointF const& end_vpos)
{
    // Calculation of angle here is done in a right-handed system. If using device
    // coordinates, i.e. a left-handed system, the angle requires a minus sign vs. the
    // caclulation. However, the advantage is, that the calculation itself can be used for
    // for right-handed systems using locical coordinates as input as well.
    //
    // Returns the relative angle between lines line are defined by beg_pos to
    // end_upos (line 1) and by beg_pos to end_vpos (line 2)

    using std::numbers::pi;

    QPointF u = end_upos - beg_pos;
    QPointF v = end_vpos - beg_pos;

    qreal nrm_prod = nrm(u) * nrm(v);

    qreal cos_angle = std::clamp(dot(u, v) / nrm_prod, -1.0, 1.0);
    qreal sin_angle = std::clamp(wdg(u, v) / nrm_prod, -1.0, 1.0);

    if (cos_angle >= 0.0) {

        // quadrant I or IV
        return std::asin(sin_angle);
    }
    else if (cos_angle < 0.0 && sin_angle >= 0.0) {
        // quadrant II
        return pi - std::asin(sin_angle);
    }
    else {
        // cos_angle < 0.0 && sin_angle < 0.0)
        // quadrant III
        return -pi - std::asin(sin_angle);
    }
}

QPainterPath arrowLine(QPointF const& beg_pos, QPointF const& end_pos)
{
    // return line between beg_pos and end_pos to draw a vector

    using std::numbers::pi;

    qreal angle = angle_of_line(beg_pos, end_pos);

    QPainterPath path;
    path.moveTo(beg_pos);
    path.lineTo(end_pos - QPointF(sin(angle + pi / 2) * 0.7 * ARROWSIZE,
                                  cos(angle + pi / 2) * 0.7 * ARROWSIZE));
    return path;
}

QPainterPath arrowHead(QPointF const& beg_pos, QPointF const& end_pos)
{
    // return arrow head at end_pos to draw a vector
    // direction of vector head is defined by beg_pos and end_pos

    using std::numbers::pi;

    qreal angle = angle_of_line(beg_pos, end_pos);

    QPainterPath path;
    path.moveTo(end_pos);
    path.lineTo(end_pos - QPointF(sin(angle + pi / 2.5) * ARROWSIZE,
                                  cos(angle + pi / 2.5) * ARROWSIZE));
    path.lineTo(end_pos - QPointF(sin(angle + pi / 2) * 0.7 * ARROWSIZE,
                                  cos(angle + pi / 2) * 0.7 * ARROWSIZE));
    path.lineTo(end_pos - QPointF(sin(angle + pi - pi / 2.5) * ARROWSIZE,
                                  cos(angle + pi - pi / 2.5) * ARROWSIZE));
    path.closeSubpath();
    return path;
}

QPainterPath vectorShape(QPointF const& beg_pos, QPointF const& end_pos)
{
    using std::numbers::pi;

    // return shape of vector incl. line and head between beg_pos and end_pos
    qreal angle = angle_of_line(beg_pos, end_pos);

    QPainterPath path;
    path.moveTo(beg_pos -
                QPointF(sin(angle + pi / 2) * MARGIN, cos(angle + pi / 2) * MARGIN) -
                QPointF(-cos(angle + pi / 2) * MARGIN, sin(angle + pi / 2) * MARGIN));
    path.lineTo(end_pos +
                QPointF(sin(angle + pi / 2) * MARGIN, cos(angle + pi / 2) * MARGIN) -
                QPointF(-cos(angle + pi / 2) * MARGIN, sin(angle + pi / 2) * MARGIN));
    path.lineTo(end_pos +
                QPointF(sin(angle + pi / 2) * MARGIN, cos(angle + pi / 2) * MARGIN) +
                QPointF(-cos(angle + pi / 2) * MARGIN, sin(angle + pi / 2) * MARGIN));
    path.lineTo(beg_pos -
                QPointF(sin(angle + pi / 2) * MARGIN, cos(angle + pi / 2) * MARGIN) +
                QPointF(-cos(angle + pi / 2) * MARGIN, sin(angle + pi / 2) * MARGIN));
    path.closeSubpath();
    return path;
}

QPainterPath anglePath(QPointF const& beg_pos, QPointF const& end_upos,
                       QPointF const& end_vpos)
{
    // Must use device coordinates here for the angles to visually make sense.
    // Difference to real angles is caused by aspect ratio of selected coordinate system
    // on screen. Angle sign differences compared to physical angles are caused by
    // switching between left- and right-handed coordinate systems

    using std::numbers::pi;

    QPointF tl = QPointF(beg_pos.x() - ARCRADIUS, beg_pos.y() - ARCRADIUS);
    QRectF br = QRectF(tl, QSizeF(2 * ARCRADIUS, 2 * ARCRADIUS));

    // sign relative to x-axis (angle > 0 is counterclockwise vs. x-axis)
    qreal angle_u = angle_of_line(beg_pos, end_upos) * 180 / pi;
    // sign has to be reversed here since device coordinates are in left-handed system.
    // (calculation itself is done in classical right-handed system)
    qreal angle_sweep = -angle_between_lines(beg_pos, end_upos, end_vpos) * 180 / pi;
    qreal angle_delta = ARCDELTA * (angle_sweep / 180); // scale to 0 for small angles

    // qDebug() << "angle_u: " << angle_u << ", angle_sweep: " << angle_sweep
    //          << ", angle_delta: " << angle_delta;

    QPainterPath path;
    path.moveTo(beg_pos);
    path.arcMoveTo(br, angle_u + angle_delta);
    path.arcTo(br, angle_u + angle_delta, angle_sweep - 2 * angle_delta);
    path.addEllipse(path.currentPosition(), ARCMARK, ARCMARK);
    return path;
}

QPainterPath planeLine(QPointF const& beg_pos, QPointF const& end_pos, qreal max_size)
{
    // return line normal to vector between beg_pos and end_pos

    using std::numbers::pi;

    qreal angle = angle_of_line(beg_pos, end_pos) - pi / 2;

    // qDebug() << "planeLine: angle = " << angle;

    QPainterPath path;
    path.moveTo(beg_pos + QPointF(cos(angle) * max_size, -sin(angle) * max_size));
    path.lineTo(beg_pos - QPointF(cos(angle) * max_size, -sin(angle) * max_size));
    return path;
}