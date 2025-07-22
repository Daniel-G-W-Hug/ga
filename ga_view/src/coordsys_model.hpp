#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <QColor>
#include <QPen>

#include <cassert> // attribute [[maybe_unused]]
#include <vector>

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

// use the ga framework
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"
using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)

using namespace fmt::literals; // just to make the format literals visible

enum Symbol { plus, cross, circle, square };

using pt2d = hd::ga::Vec2d<value_t>;
using pt2dp = hd::ga::Vec2dp<value_t>;
using bivt2dp = hd::ga::BiVec2dp<value_t>;

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

// this struct should be used by the user to mark projective lines (=bivectors)
struct bivt2dp_mark {

    QPen pen{QPen(Qt::gray, 1, Qt::SolidLine)};
};

const bivt2dp_mark bivt2dp_mark_default; // for default arguments;


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

// active projective vector as directed line between two active points
struct abivt2dp {

    // reference to points used here, because this enables
    // use of unique address of pt2d for mapping to users of active points in scene
    size_t beg_idx;
    size_t end_idx;
};

// active projective reflection using four active points
struct arefl2dp {

    size_t p1_idx; // bvt1 = wdg(p1,p2)
    size_t p2_idx;
    size_t p3_idx; // bvt2 = wdg(p2,p3)
    size_t p4_idx;
};

// active open kinematics chain using an arbitrary number of active points
struct akinematic2dp {

    // first point is base, last point is target
    // open kinematic chain is defined until one point before target
    std::vector<size_t> apt; // index of active point
};

// ----------------------------------------------------------------------------
// convenience alias to make pt2d and ln2d look similar
// convenience alias to make pt2dp and ln2dp look similar
// ----------------------------------------------------------------------------
using ln2d = std::vector<pt2d>;
using cln2dp = std::vector<pt2dp>;
// ----------------------------------------------------------------------------


class Coordsys_model {
  public:

    // add passive point
    [[maybe_unused]] size_t add_pt(pt2d const& p_in,
                                   pt2d_mark const& m = pt2d_mark_default);
    [[maybe_unused]] size_t add_pt(pt2dp const& p_in,
                                   pt2d_mark const& m = pt2d_mark_default);

    // add passive line
    [[maybe_unused]] size_t add_ln(std::vector<pt2d> const& vp_in,
                                   ln2d_mark const& m = ln2d_mark_default);
    [[maybe_unused]] size_t add_ln(std::vector<pt2dp> const& vp_in,
                                   ln2d_mark const& m = ln2d_mark_default);

    // add passive vector
    [[maybe_unused]] size_t add_vt(vt2d const& vt_in,
                                   vt2d_mark const& m = vt2d_mark_default);

    // add passive bivector, i.e. a projective line
    [[maybe_unused]] size_t add_bivtp(bivt2dp const& bivtp_in,
                                      bivt2dp_mark const& m = bivt2dp_mark_default);

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

    // add active projective bivector (= a projective line)
    [[maybe_unused]] size_t add_abivtp(abivt2dp const& abivtp_in);

    // add active projective reflection
    [[maybe_unused]] size_t add_areflp(arefl2dp const& reflp_in);

    // add active open kinematic chaing
    [[maybe_unused]] size_t add_akinem(akinematic2dp const& akinem_in);

    void set_label(std::string new_label) { m_label = std::move(new_label); };
    std::string label() const { return m_label; }

    // reset model to empty state, e.g. for reuse in new model
    void clear();

    // data for points (same index is for same point)
    std::vector<pt2d> pt;
    std::vector<pt2d_mark> pt_mark;

    std::vector<pt2dp> ptp;
    std::vector<pt2d_mark> ptp_mark;

    // data for lines consisting of points (same index is for same line)
    std::vector<ln2d> ln;
    std::vector<ln2d_mark> ln_mark;

    std::vector<cln2dp> clnp;
    std::vector<ln2d_mark> clnp_mark;

    // data for vectors (same index is for same vector)
    std::vector<vt2d> vt;
    std::vector<vt2d_mark> vt_mark;

    // data for projective lines (=bivectors)
    std::vector<bivt2dp> bivtp;
    std::vector<bivt2dp_mark> bivtp_mark;

    // data for active points (same index is for same point)
    std::vector<pt2d> apt;

    // data for active vectors using active points (same index is for same point)
    std::vector<avt2d> avt;

    // data for active bivectors using active points (same index is for same point)
    std::vector<abivt2d> abivt;

    // data for active reflections using active points (same index is for same point)
    std::vector<aproj2d> aproj;

    // data for active reflections using active points (same index is for same point)
    std::vector<arefl2d> arefl;

    // data for active projective bivectors using active points
    std::vector<abivt2dp> abivtp;

    // data for active projective reflections using active points
    std::vector<arefl2dp> areflp;

    // data for active projective reflections using active points
    std::vector<akinematic2dp> akinemp;

    // model label, e.g. time stamp description of current Coordsys_model
    std::string m_label{};
};

// ----------------------------------------------------------------------------
// printing support via fmt library
// ----------------------------------------------------------------------------

// formating for user defined types (pt2d)
template <> struct fmt::formatter<pt2d> : nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(pt2d const& pt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<pt2d>::parse(ParseContext& ctx)
{
    return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<pt2d>::format(pt2d const& pt, FormatContext& ctx) const
{
    return fmt::format_to(ctx.out(), "pt2d({}, {})", pt.x, pt.y);
}

// formating for user defined types (pt2dp)
template <> struct fmt::formatter<pt2dp> : nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(const pt2dp& pt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<pt2dp>::parse(ParseContext& ctx)
{
    return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<pt2dp>::format(pt2dp const& pt, FormatContext& ctx) const
{
    return fmt::format_to(ctx.out(), "pt2dp({}, {}, {})", pt.x, pt.y, pt.z);
}


// formating for user defined types (vt2d)
template <> struct fmt::formatter<vt2d> : nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(vt2d const& vt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<vt2d>::parse(ParseContext& ctx)
{
    return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<vt2d>::format(vt2d const& vt, FormatContext& ctx) const
{
    return fmt::format_to(ctx.out(), "vt2d(pt2d({}, {}), pt2d({}, {}))", vt.beg.x,
                          vt.beg.y, vt.end.x, vt.end.y);
}

// formating for user defined types (bivt2dp)
template <> struct fmt::formatter<bivt2dp> : nested_formatter<double> {
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx);
    template <typename FormatContext>
    auto format(bivt2dp const& bivt, FormatContext& ctx) const;
};

template <typename ParseContext>
constexpr auto fmt::formatter<bivt2dp>::parse(ParseContext& ctx)
{
    return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<bivt2dp>::format(bivt2dp const& bivt, FormatContext& ctx) const
{
    return fmt::format_to(ctx.out(), "bivt2dp({}, {}, {})", bivt.x, bivt.y, bivt.z);
}

// Bsp. für Anwendung
// std::vector<pt2d> vp1{{1.0, 1.0}, {1.5, 2.0}};
// pt2d p{1.0, 2.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));