// Test to understand the original bivector orientation logic

#include "ga/ga_pga.hpp"
#include <algorithm>
#include <cmath>
#include <fmt/core.h>
#include <numbers>
#include <vector>

using namespace hd::ga;
using namespace hd::ga::pga;
using std::numbers::pi;

int main()
{
    // Fixed viewport
    value_t x_min = -5.0, x_max = 5.0;
    value_t y_min = -5.0, y_max = 5.0;

    // Test 1: Fixed p at origin, rotate q around it
    fmt::println("=== Test 1: p fixed at (0,0), q rotating every 5° ===");
    auto p_fixed = vec2dp(0.0, 0.0, 1.0);
    value_t radius = 3.0;

    // Define viewport boundaries for simulation
    auto p1 = vec2dp(x_min, y_max, 1.0);
    auto p2 = vec2dp(x_max, y_max, 1.0);
    auto p3 = vec2dp(x_min, y_min, 1.0);
    auto p4 = vec2dp(x_max, y_min, 1.0);

    auto l_top = unitize(wdg(p1, p2));
    auto l_left = unitize(wdg(p3, p1));
    auto l_bottom = unitize(wdg(p3, p4));
    auto l_right = unitize(wdg(p4, p2));

    for (int angle_deg = 0; angle_deg <= 360; angle_deg += 5) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto q_rotating =
            vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_fixed, q_rotating);

        // Calculate line angle using atan2
        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // Simulate intersection calculation
        vec2dp p_left = unitize(rwdg(l_left, bvt));
        vec2dp p_right = unitize(rwdg(l_right, bvt));
        vec2dp p_bottom = unitize(rwdg(l_bottom, bvt));
        vec2dp p_top = unitize(rwdg(l_top, bvt));

        auto is_in_bounds = [&](vec2dp const& pt) {
            return pt.x >= x_min - eps && pt.x <= x_max + eps && pt.y >= y_min - eps &&
                   pt.y <= y_max + eps;
        };

        std::vector<vec2dp> valid_points;
        if (is_in_bounds(p_left)) valid_points.push_back(p_left);
        if (is_in_bounds(p_right)) valid_points.push_back(p_right);
        if (is_in_bounds(p_bottom)) valid_points.push_back(p_bottom);
        if (is_in_bounds(p_top)) valid_points.push_back(p_top);

        vec2dp p_from, p_to;
        if (valid_points.size() >= 2) {
            value_t cos_line = std::cos(line_angle);
            value_t sin_line = std::sin(line_angle);

            auto projection = [&](vec2dp const& pt) {
                return pt.x * cos_line + pt.y * sin_line;
            };

            auto min_it = std::min_element(valid_points.begin(), valid_points.end(),
                                           [&](vec2dp const& a, vec2dp const& b) {
                                               return projection(a) < projection(b);
                                           });
            auto max_it = std::max_element(valid_points.begin(), valid_points.end(),
                                           [&](vec2dp const& a, vec2dp const& b) {
                                               return projection(a) < projection(b);
                                           });

            p_from = *min_it;
            p_to = *max_it;

            fmt::println("{:3}°: line_ang={:+7.2f}°, p_from=({:+.2f},{:+.2f}), "
                         "p_to=({:+.2f},{:+.2f}), q=({:+.2f},{:+.2f})",
                         angle_deg, line_angle * 180.0 / pi, p_from.x, p_from.y, p_to.x,
                         p_to.y, q_rotating.x, q_rotating.y);
        }
    }

    fmt::println("\n=== Test 2: q fixed at (0,0), p rotating every 5° ===");
    auto q_fixed = vec2dp(0.0, 0.0, 1.0);

    for (int angle_deg = 0; angle_deg <= 360; angle_deg += 5) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto p_rotating =
            vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_rotating, q_fixed);

        // Calculate line angle using atan2
        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // Simulate intersection calculation
        vec2dp p_left = unitize(rwdg(l_left, bvt));
        vec2dp p_right = unitize(rwdg(l_right, bvt));
        vec2dp p_bottom = unitize(rwdg(l_bottom, bvt));
        vec2dp p_top = unitize(rwdg(l_top, bvt));

        auto is_in_bounds = [&](vec2dp const& pt) {
            return pt.x >= x_min - eps && pt.x <= x_max + eps && pt.y >= y_min - eps &&
                   pt.y <= y_max + eps;
        };

        std::vector<vec2dp> valid_points;
        if (is_in_bounds(p_left)) valid_points.push_back(p_left);
        if (is_in_bounds(p_right)) valid_points.push_back(p_right);
        if (is_in_bounds(p_bottom)) valid_points.push_back(p_bottom);
        if (is_in_bounds(p_top)) valid_points.push_back(p_top);

        vec2dp p_from, p_to;
        if (valid_points.size() >= 2) {
            value_t cos_line = std::cos(line_angle);
            value_t sin_line = std::sin(line_angle);

            auto projection = [&](vec2dp const& pt) {
                return pt.x * cos_line + pt.y * sin_line;
            };

            auto min_it = std::min_element(valid_points.begin(), valid_points.end(),
                                           [&](vec2dp const& a, vec2dp const& b) {
                                               return projection(a) < projection(b);
                                           });
            auto max_it = std::max_element(valid_points.begin(), valid_points.end(),
                                           [&](vec2dp const& a, vec2dp const& b) {
                                               return projection(a) < projection(b);
                                           });

            p_from = *min_it;
            p_to = *max_it;

            fmt::println("{:3}°: line_ang={:+7.2f}°, p_from=({:+.2f},{:+.2f}), "
                         "p_to=({:+.2f},{:+.2f}), p=({:+.2f},{:+.2f})",
                         angle_deg, line_angle * 180.0 / pi, p_from.x, p_from.y, p_to.x,
                         p_to.y, p_rotating.x, p_rotating.y);
        }
    }

    return 0;
}
