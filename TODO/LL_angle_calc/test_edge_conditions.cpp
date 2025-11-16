// Test to verify edge case conditions for intersection calculation

#include "ga/ga_pga.hpp"
#include <fmt/core.h>
#include <numbers>

using namespace hd::ga;
using namespace hd::ga::pga;
using std::numbers::pi;

int main()
{
    auto p_origin = vec2dp(0.0, 0.0, 1.0);
    value_t radius = 3.0;

    fmt::println("Testing edge case conditions:");
    fmt::println("{:>5} {:>10} {:>10} {:>15} {:>15} {:>25} {:>25}", "angle", "phi_x",
                 "phi_y", "line_angle", "|phi_x-π/2|", "orig: calc p_from/to_24?",
                 "orig: calc p_from/to_13?");
    fmt::println("{:-^120}", "");

    for (int angle_deg = 0; angle_deg <= 360; angle_deg += 15) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        // Original approach
        auto const x_axis = bivec2dp{0, 1, 0};
        auto const y_axis = bivec2dp{1, 0, 0};
        auto const phi_x = angle(x_axis, bvt);
        auto const phi_y = angle(y_axis, bvt);

        // New approach
        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // Original conditions for calculating intersections
        bool calc_24_orig = (std::abs(phi_x - pi / 2.0) > eps); // not on ±y-axis
        bool calc_13_orig = (std::abs(phi_y - pi / 2.0) > eps); // not on ±x-axis

        fmt::println("{:5.0f} {:10.1f} {:10.1f} {:+15.1f} {:15.6f} {:>25} {:>25}",
                     (double)angle_deg, phi_x * 180.0 / pi, phi_y * 180.0 / pi,
                     line_angle * 180.0 / pi, std::abs(phi_x - pi / 2.0),
                     calc_24_orig ? "YES (calc p_24)" : "NO (skip p_24)",
                     calc_13_orig ? "YES (calc p_13)" : "NO (skip p_13)");
    }

    fmt::println("\n\nKey observations:");
    fmt::println("Original condition for p_from_24/p_to_24: |phi_x - π/2| > eps");
    fmt::println("  - This is FALSE when phi_x ≈ π/2 (i.e., on ±y-axis)");
    fmt::println("  - phi_x = π/2 occurs when line_angle = ±π or π/2");
    fmt::println("\nOriginal condition for p_from_13/p_to_13: |phi_y - π/2| > eps");
    fmt::println("  - This is FALSE when phi_y ≈ π/2 (i.e., on ±x-axis)");
    fmt::println("  - phi_y = π/2 occurs when line_angle = -π/2 or π/2");

    return 0;
}
