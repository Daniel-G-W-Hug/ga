#pragma once

#include <QColor>
#include <QPen>

#include <cassert> // attribute [[maybe_unused]]
#include <vector>

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

using namespace fmt::literals; // just to make the format literals visible

enum Symbol { plus, cross, circle, square };

struct pt2d { // coordinates of point on x and y axis
    double x{0.0}, y{0.0};
    pt2d(double x_in, double y_in) : x(x_in), y(y_in) {}
    pt2d(pt2d const&) = default;
    pt2d& operator=(pt2d const&) = default;
    pt2d(pt2d&&) = default;
    pt2d& operator=(pt2d&&) = default;
    pt2d() = default;
    ~pt2d() noexcept = default;

    bool operator==(pt2d const& rhs) const
    {
        auto abs_delta_x = std::abs(rhs.x - x);
        auto abs_delta_y = std::abs(rhs.y - y);
        if (abs_delta_x == 0.0 && abs_delta_y == 0.0) return true;
        return false;
    }
};

// this struct should be used by the user to mark points
struct pt2d_mark {

    Symbol symbol{plus}; // define the marking symbol
    int nsize{4};        // characteristic size of mark symbol in pixels
    QPen pen{QPen(Qt::blue, 2, Qt::SolidLine)};

    int grp{0}; // user provided group the
                // item shall belong to (for selection)
};

const pt2d_mark pt2d_mark_default; // for default arguments

// this struct should be used by the user to mark lines
struct ln2d_mark {

    QPen pen{QPen(Qt::black, 2, Qt::SolidLine)};

    bool mark_pts{false};
    int delta{1}; // 1 shows every point, 2 every second, ...
    pt2d_mark pm{};

    bool mark_area{false};
    QColor area_col{QColor(0, 128, 0, 128)};

    int grp{0}; // user provided group the
                // item shall belong to (for selection)
};

const ln2d_mark ln2d_mark_default; // for default arguments;

// vector as directed line
struct vt2d {
    pt2d beg; // beg = (0,0) => position vector
              // else => free vector
    pt2d end;

    vt2d(pt2d const& p2) : beg{pt2d{0.0, 0.0}}, end{p2} {}     // position vector
    vt2d(pt2d const& p1, pt2d const& p2) : beg{p1}, end{p2} {} // vector from beg to end

    vt2d(vt2d const&) = default;
    vt2d& operator=(vt2d const&) = default;
    vt2d(vt2d&&) = default;
    vt2d& operator=(vt2d&&) = default;
    vt2d() = default;
    ~vt2d() noexcept = default;
};

// this struct should be used by the user to mark vectors
struct vt2d_mark {

    QPen pen{QPen(Qt::black, 2, Qt::SolidLine)};

    int grp{0}; // user provided group the
                // item shall belong to (for selection)
};

const vt2d_mark vt2d_mark_default; // for default arguments;

// ----------------------------------------------------------------------------
// this is used internally, not by the user directly
// ----------------------------------------------------------------------------

struct mark_id {

    // system provided unique_id to identify each item in the model
    int id{-1};           // valid unique ids are positive values
    int linked_to_id{-1}; // associated with other id (e.g. pt_marks with line)
    bool active{true};    // active elements are diplayed (on by default)
};

// ----------------------------------------------------------------------------
// this is used internally up to here, not by the user directly
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// convenience alias to make pt2d and ln2d look similar
// ----------------------------------------------------------------------------
using ln2d = std::vector<pt2d>;
// ----------------------------------------------------------------------------

class Coordsys_model {
  public:

    // add passive point
    [[maybe_unused]] int add_pt(const pt2d& p_in, const pt2d_mark m = pt2d_mark_default);

    // add passive line
    [[maybe_unused]] int add_ln(const std::vector<pt2d>& vp_in,
                                const ln2d_mark m = ln2d_mark_default);
    // add passive vector
    [[maybe_unused]] int add_vt(const vt2d& v_in, const vt2d_mark m = vt2d_mark_default);

    // add active point
    [[maybe_unused]] int add_apt(const pt2d& p_in);

    void set_label(std::string new_label) { m_label = new_label; };
    std::string label() { return m_label; }

    // reset model to empty state, e.g. for reuse in new model
    void clear();

    int unique_id{-1}; // id = unique id with rng: [0..) to identify each item in model.
                       // Assigned when model is setup using push_back calls.

    // data for points (same index is for same point)
    std::vector<pt2d> pt;
    std::vector<pt2d_mark> pt_mark;
    std::vector<mark_id> pt_id;

    // data for lines (same index is for same line)
    std::vector<ln2d> ln;
    std::vector<ln2d_mark> ln_mark;
    std::vector<mark_id> ln_id;

    // data for vectors (same index is for same vector)
    std::vector<vt2d> vt;
    std::vector<vt2d_mark> vt_mark;
    std::vector<mark_id> vt_id;

    // data for active points (same index is for same point)
    std::vector<pt2d> apt;
    std::vector<mark_id> apt_id;

    // model label (e.g. time stamp description)
    std::string m_label{};
};

// ----------------------------------------------------------------------------
// printing support via fmt library
// ----------------------------------------------------------------------------

// formating for user defined types (pt2d)
template <> struct fmt::formatter<pt2d> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext> auto format(const pt2d& pt, FormatContext& ctx);
};

template <typename ParseContext>
constexpr auto fmt::formatter<pt2d>::parse(ParseContext& ctx)
{
    return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<pt2d>::format(const pt2d& pt, FormatContext& ctx)
{
    // return fmt::format_to(ctx.out(), "{{{0}, {1}}}", pt.x, pt.y);
    return fmt::format_to(ctx.out(), "({0}, {1})", pt.x, pt.y);
}

// Bsp. für Anwendung
// std::vector<pt2d> vp1{{1.0, 1.0}, {1.5, 2.0}};
// pt2d p{1.0, 2.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));