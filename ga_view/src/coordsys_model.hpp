#pragma once
//
// author: Daniel Hug, 2024
//

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
    // pt2d(pt2d const&) = default;
    // pt2d& operator=(pt2d const&) = default;
    // pt2d(pt2d&&) = default;
    // pt2d& operator=(pt2d&&) = default;
    pt2d() = default;
    // ~pt2d() noexcept = default;

    bool operator==(pt2d const& rhs) const
    {
        if (x == rhs.x && y == rhs.y) return true;
        return false;
    }
};

// this struct should be used by the user to mark points
struct pt2d_mark {

    QPen pen{QPen(Qt::blue, 2, Qt::SolidLine)};
    Symbol symbol{plus}; // define the marking symbol
    int nsize{4};        // characteristic size of mark symbol in pixels
};

static pt2d_mark const pt2d_mark_default; // for default arguments

// this struct should be used by the user to mark lines
struct ln2d_mark {

    QPen pen{QPen(Qt::black, 2, Qt::SolidLine)};
    QColor area_col{QColor(0, 128, 0, 128)};

    pt2d_mark pm{};
    size_t delta{1}; // 1 shows every point, 2 every second, ...

    bool mark_pts{false};
    bool mark_area{false};
};

static ln2d_mark const ln2d_mark_default; // for default arguments;

// vector as directed line between two points
struct vt2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points later on
    pt2d beg; // beg = (0,0) => position vector, otherwise => free vector
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
};

const vt2d_mark vt2d_mark_default; // for default arguments;


// active vector as directed line between two active points
struct avt2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t beg_idx; // beg = (0,0) => position vector, otherwise => free vector
    size_t end_idx;
};


// active bivector between three active points (one common to both vectors)
struct abivt2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t beg_idx;
    size_t uend_idx;
    size_t vend_idx;
};

// active reflection using three active points (one common to both vectors)
struct aproj2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t beg_idx;
    size_t uend_idx;
    size_t vend_idx;
};

// active reflection using three active points (one common to both vectors)
struct arefl2d {

    // reference to point used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t n1end_idx;
    size_t n2end_idx;
};

// ----------------------------------------------------------------------------
// convenience alias to make pt2d and ln2d look similar
// ----------------------------------------------------------------------------
using ln2d = std::vector<pt2d>;
// ----------------------------------------------------------------------------


class Coordsys_model {
  public:

    // add passive point
    [[maybe_unused]] size_t add_pt(pt2d const& p_in,
                                   pt2d_mark const& m = pt2d_mark_default);

    // add passive line
    [[maybe_unused]] size_t add_ln(std::vector<pt2d> const& vp_in,
                                   ln2d_mark const& m = ln2d_mark_default);
    // add passive vector
    [[maybe_unused]] size_t add_vt(vt2d const& vt_in,
                                   vt2d_mark const& m = vt2d_mark_default);

    // add active point
    [[maybe_unused]] size_t add_apt(pt2d const& pt_in);

    // add active vector
    [[maybe_unused]] size_t add_avt(avt2d const& vt_in);

    // add active bivector
    [[maybe_unused]] size_t add_abivt(abivt2d const& bivt_in);

    // add active projection
    [[maybe_unused]] size_t add_aproj(aproj2d const& proj_in);

    // add active reflection
    [[maybe_unused]] size_t add_arefl(arefl2d const& refl_in);

    void set_label(std::string new_label) { m_label = std::move(new_label); };
    std::string label() { return m_label; }

    // reset model to empty state, e.g. for reuse in new model
    void clear();

    // data for points (same index is for same point)
    std::vector<pt2d> pt;
    std::vector<pt2d_mark> pt_mark;

    // data for lines consisting of points (same index is for same line)
    std::vector<ln2d> ln;
    std::vector<ln2d_mark> ln_mark;

    // data for vectors (same index is for same vector)
    std::vector<vt2d> vt;
    std::vector<vt2d_mark> vt_mark;

    // data for active points (same index is for same point)
    std::vector<pt2d> apt;

    // data for active vectors using active points (same index is for same point)
    std::vector<avt2d> avt;

    // data for active vectors using active points (same index is for same point)
    std::vector<abivt2d> abivt;

    // data for active reflections using active points (same index is for same point)
    std::vector<aproj2d> aproj;

    // data for active reflections using active points (same index is for same point)
    std::vector<arefl2d> arefl;

    // model label, e.g. time stamp description of current Coordsys_model
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
    return fmt::format_to(ctx.out(), "({}, {})", pt.x, pt.y);
}


// formating for user defined types (vt2d)
template <> struct fmt::formatter<vt2d> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext> auto format(const vt2d& vt, FormatContext& ctx);
};

template <typename ParseContext>
constexpr auto fmt::formatter<vt2d>::parse(ParseContext& ctx)
{
    return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<vt2d>::format(const vt2d& vt, FormatContext& ctx)
{
    return fmt::format_to(ctx.out(), "(({}, {}),({}, {}))", vt.beg.x, vt.beg.y, vt.end.x,
                          vt.end.y);
}

// Bsp. für Anwendung
// std::vector<pt2d> vp1{{1.0, 1.0}, {1.5, 2.0}};
// pt2d p{1.0, 2.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));