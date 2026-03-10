// Verify edge cases: lines exactly on axes

#include "ga/ga_pga.hpp"
#include <fmt/core.h>
#include <numbers>

using namespace hd::ga;
using namespace hd::ga::pga;
using std::numbers::pi;

int main()
{
    value_t x_min = -5.0, x_max = 5.0;
    value_t y_min = -5.0, y_max = 5.0;

    auto p_origin = vec2dp(0.0, 0.0, 1.0);
    value_t radius = 3.0;

    // Define viewport boundaries
    auto p1 = vec2dp(x_min, y_max, 1.0);
    auto p2 = vec2dp(x_max, y_max, 1.0);
    auto p3 = vec2dp(x_min, y_min, 1.0);
    auto p4 = vec2dp(x_max, y_min, 1.0);

    auto l1 = unitize(wdg(p1, p2)); // top
    auto l2 = unitize(wdg(p3, p1)); // left
    auto l3 = unitize(wdg(p3, p4)); // bottom
    auto l4 = unitize(wdg(p4, p2)); // right

    struct TestCase {
        double angle_deg;
        std::string axis_name;
    };

    TestCase cases[] = {
        {0, "pos x-axis"}, {90, "pos y-axis"}, {180, "neg x-axis"}, {270, "neg y-axis"}};

    fmt::println("Testing edge cases (lines exactly on axes):");
    fmt::println("{:>5} {:>15} {:>10} {:>15} {:>15} {:>25} {:>25}", "angle", "axis",
                 "line_ang", "calc p_24?", "calc p_13?", "p_from", "p_to");
    fmt::println("{:-^130}", "");

    for (auto& tc : cases) {
        value_t angle_rad = tc.angle_deg * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // Check intersection calculation conditions (corrected)
        bool calc_24 =
            !(std::abs(std::abs(line_angle) - pi) < eps || std::abs(line_angle) < eps);
        bool calc_13 = !(std::abs(line_angle - pi / 2.0) < eps ||
                         std::abs(line_angle + pi / 2.0) < eps);

        vec2dp p_from_24, p_to_24, p_from_13, p_to_13;
        vec2dp p_from, p_to;

        // Calculate intersections based on conditions
        if (calc_24) {
            p_from_24 = unitize(rwdg(l2, bvt));
            p_to_24 = unitize(rwdg(l4, bvt));
        }

        if (calc_13) {
            p_from_13 = unitize(rwdg(l3, bvt));
            p_to_13 = unitize(rwdg(l1, bvt));
        }

        // Determine final endpoints based on axis
        if (std::abs(line_angle + pi / 2.0) < eps) { // pos x-axis
            p_from = p_from_24;
            p_to = p_to_24;
        }
        else if (std::abs(std::abs(line_angle) - pi) < eps) { // pos y-axis
            p_from = p_from_13;
            p_to = p_to_13;
        }
        else if (std::abs(line_angle - pi / 2.0) < eps) { // neg x-axis
            p_from = p_to_24;
            p_to = p_from_24;
        }
        else if (std::abs(line_angle) < eps) { // neg y-axis
            p_from = p_to_13;
            p_to = p_from_13;
        }

        fmt::println("{:5.0f} {:>15} {:+10.1f} {:>15} {:>15} ({:+6.2f},{:+6.2f}) "
                     "({:+6.2f},{:+6.2f})",
                     tc.angle_deg, tc.axis_name, line_angle * 180.0 / pi,
                     calc_24 ? "YES" : "NO", calc_13 ? "YES" : "NO", p_from.x, p_from.y,
                     p_to.x, p_to.y);
    }

    fmt::println("\n\nExpected results:");
    fmt::println("  0° (pos x-axis): p_from=(-5, 0), p_to=(+5, 0)  [horizontal line, "
                 "calc p_24 only]");
    fmt::println(" 90° (pos y-axis): p_from=(0, -5), p_to=(0, +5)  [vertical line, calc "
                 "p_13 only]");
    fmt::println("180° (neg x-axis): p_from=(+5, 0), p_to=(-5, 0)  [horizontal line, "
                 "calc p_24 skipped!]");
    fmt::println("270° (neg y-axis): p_from=(0, +5), p_to=(0, -5)  [vertical line, calc "
                 "p_13 only]");

    return 0;
}
