// Verify complete equivalence: compare actual p_from/p_to between original and new

#include "ga/ga_pga.hpp"
#include <fmt/core.h>
#include <numbers>

using namespace hd::ga;
using namespace hd::ga::pga;
using std::numbers::pi;

int main()
{
    // Fixed viewport
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

    fmt::println("COMPLETE EQUIVALENCE TEST");
    fmt::println("=========================\n");
    fmt::println("Comparing ACTUAL p_from/p_to outputs between original and new "
                 "implementations.\n");

    int mismatch_count = 0;
    int total_tests = 0;

    fmt::println("{:>10} {:>12} {:>30} {:>30} {:>10}", "q angle", "line_angle",
                 "ORIGINAL p_from → p_to", "NEW p_from → p_to", "Match?");
    fmt::println("{:-^105}", "");

    for (int angle_deg = 0; angle_deg <= 360; angle_deg += 5) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        total_tests++;

        // ============================================================
        // ORIGINAL IMPLEMENTATION
        // ============================================================
        auto const x_axis = bivec2dp{0, 1, 0};
        auto const y_axis = bivec2dp{1, 0, 0};
        auto const phi_x = angle(x_axis, bvt);
        auto const phi_y = angle(y_axis, bvt);

        vec2dp orig_p_from, orig_p_to;
        vec2dp orig_p_from_24, orig_p_to_24, orig_p_from_13, orig_p_to_13;

        // Calculate intersections based on original conditions
        if (std::abs(phi_x - pi / 2.0) > eps) {
            orig_p_from_24 = unitize(rwdg(l2, bvt));
            orig_p_to_24 = unitize(rwdg(l4, bvt));
        }
        if (std::abs(phi_y - pi / 2.0) > eps) {
            orig_p_from_13 = unitize(rwdg(l3, bvt));
            orig_p_to_13 = unitize(rwdg(l1, bvt));
        }

        // Original case logic
        if ((std::abs(phi_x) < eps) && (std::abs(phi_y - pi / 2.0) < eps)) {
            orig_p_from = orig_p_from_24;
            orig_p_to = orig_p_to_24;
        }
        else if (phi_x > eps && (phi_x < (pi / 2. - eps)) && (phi_y > (pi / 2. + eps)) &&
                 (phi_y < (pi - eps))) {
            if (orig_p_to_13.x < x_max) {
                orig_p_to = orig_p_to_13;
            }
            else {
                orig_p_to = orig_p_to_24;
            }
            if (orig_p_from_13.x > x_min) {
                orig_p_from = orig_p_from_13;
            }
            else {
                orig_p_from = orig_p_from_24;
            }
        }
        else if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y - pi) < eps)) {
            orig_p_from = orig_p_from_13;
            orig_p_to = orig_p_to_13;
        }
        else if ((phi_x > (pi / 2.0 + eps)) && (phi_x < (pi - eps)) &&
                 (phi_y < (pi - eps)) && (phi_y > (pi / 2.0 + eps))) {
            if (orig_p_from_24.y < y_max) {
                orig_p_to = orig_p_from_24;
            }
            else {
                orig_p_to = orig_p_to_13;
            }
            if (orig_p_to_24.y > y_min) {
                orig_p_from = orig_p_to_24;
            }
            else {
                orig_p_from = orig_p_from_13;
            }
        }
        else if ((std::abs(phi_x - pi) < eps) && (std::abs(phi_y - pi / 2.0)) < eps) {
            orig_p_from = orig_p_to_24;
            orig_p_to = orig_p_from_24;
        }
        else if ((phi_x < pi - eps) && (phi_x > pi / 2.0 + eps) &&
                 (phi_y < pi / 2.0 - eps) && (phi_y > eps)) {
            if (orig_p_from_13.x > x_min) {
                orig_p_to = orig_p_from_13;
            }
            else {
                orig_p_to = orig_p_from_24;
            }
            if (orig_p_to_13.x < x_max) {
                orig_p_from = orig_p_to_13;
            }
            else {
                orig_p_from = orig_p_to_24;
            }
        }
        else if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y) < eps)) {
            orig_p_from = orig_p_to_13;
            orig_p_to = orig_p_from_13;
        }
        else if ((phi_x < pi / 2.0 - eps) && (phi_x > eps) && (phi_y > eps) &&
                 (phi_y < pi / 2. - eps)) {
            if (orig_p_to_24.y > y_min) {
                orig_p_to = orig_p_to_24;
            }
            else {
                orig_p_to = orig_p_from_13;
            }
            if (orig_p_from_24.y < y_max) {
                orig_p_from = orig_p_from_24;
            }
            else {
                orig_p_from = orig_p_to_13;
            }
        }

        // ============================================================
        // NEW IMPLEMENTATION
        // ============================================================
        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        vec2dp new_p_from, new_p_to;
        vec2dp new_p_from_24, new_p_to_24, new_p_from_13, new_p_to_13;

        // Calculate intersections based on new conditions
        if (!(std::abs(std::abs(line_angle) - pi) < eps || std::abs(line_angle) < eps)) {
            new_p_from_24 = unitize(rwdg(l2, bvt));
            new_p_to_24 = unitize(rwdg(l4, bvt));
        }
        if (!(std::abs(line_angle - pi / 2.0) < eps ||
              std::abs(line_angle + pi / 2.0) < eps)) {
            new_p_from_13 = unitize(rwdg(l3, bvt));
            new_p_to_13 = unitize(rwdg(l1, bvt));
        }

        // New case logic
        if (std::abs(line_angle + pi / 2.0) < eps) {
            new_p_from = new_p_from_24;
            new_p_to = new_p_to_24;
        }
        else if (line_angle < -pi / 2.0 - eps && line_angle > -pi + eps) {
            if (new_p_to_13.x < x_max) {
                new_p_to = new_p_to_13;
            }
            else {
                new_p_to = new_p_to_24;
            }
            if (new_p_from_13.x > x_min) {
                new_p_from = new_p_from_13;
            }
            else {
                new_p_from = new_p_from_24;
            }
        }
        else if (std::abs(std::abs(line_angle) - pi) < eps) {
            new_p_from = new_p_from_13;
            new_p_to = new_p_to_13;
        }
        else if (line_angle > pi / 2.0 + eps && line_angle < pi - eps) {
            if (new_p_from_24.y < y_max) {
                new_p_to = new_p_from_24;
            }
            else {
                new_p_to = new_p_to_13;
            }
            if (new_p_to_24.y > y_min) {
                new_p_from = new_p_to_24;
            }
            else {
                new_p_from = new_p_from_13;
            }
        }
        else if (std::abs(line_angle - pi / 2.0) < eps) {
            new_p_from = new_p_to_24;
            new_p_to = new_p_from_24;
        }
        else if (line_angle > eps && line_angle < pi / 2.0 - eps) {
            if (new_p_from_13.x > x_min) {
                new_p_to = new_p_from_13;
            }
            else {
                new_p_to = new_p_from_24;
            }
            if (new_p_to_13.x < x_max) {
                new_p_from = new_p_to_13;
            }
            else {
                new_p_from = new_p_to_24;
            }
        }
        else if (std::abs(line_angle) < eps) {
            new_p_from = new_p_to_13;
            new_p_to = new_p_from_13;
        }
        else if (line_angle > -pi / 2.0 + eps && line_angle < -eps) {
            if (new_p_to_24.y > y_min) {
                new_p_to = new_p_to_24;
            }
            else {
                new_p_to = new_p_from_13;
            }
            if (new_p_from_24.y < y_max) {
                new_p_from = new_p_from_24;
            }
            else {
                new_p_from = new_p_to_13;
            }
        }

        // ============================================================
        // COMPARE RESULTS
        // ============================================================
        bool p_from_match = (std::abs(orig_p_from.x - new_p_from.x) < 1e-10) &&
                            (std::abs(orig_p_from.y - new_p_from.y) < 1e-10);
        bool p_to_match = (std::abs(orig_p_to.x - new_p_to.x) < 1e-10) &&
                          (std::abs(orig_p_to.y - new_p_to.y) < 1e-10);
        bool full_match = p_from_match && p_to_match;

        if (!full_match) {
            mismatch_count++;
        }

        fmt::println("{:10.0f}° {:+12.1f}° ({:+.1f},{:+.1f})→({:+.1f},{:+.1f}) "
                     "({:+.1f},{:+.1f})→({:+.1f},{:+.1f}) {:>10}",
                     (double)angle_deg, line_angle * 180.0 / pi, orig_p_from.x,
                     orig_p_from.y, orig_p_to.x, orig_p_to.y, new_p_from.x, new_p_from.y,
                     new_p_to.x, new_p_to.y, full_match ? "✓" : "✗ MISMATCH");
    }

    fmt::println("\n\nFINAL RESULTS:");
    fmt::println("==============");
    fmt::println("Total tests: {}", total_tests);
    fmt::println("Mismatches: {}", mismatch_count);
    fmt::println("Matches: {}", total_tests - mismatch_count);
    fmt::println("");

    if (mismatch_count == 0) {
        fmt::println("✓✓✓ PERFECT! Both implementations produce IDENTICAL results! ✓✓✓");
    }
    else {
        fmt::println(
            "✗✗✗ FAILURE! Found {} mismatches - implementations are NOT equivalent! ✗✗✗",
            mismatch_count);
    }

    return mismatch_count > 0 ? 1 : 0;
}
