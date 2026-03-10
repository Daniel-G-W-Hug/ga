// Systematic walkthrough of all angles for line orientations in PGA

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

    fmt::println("Systematic angle mapping for PGA lines:");
    fmt::println("========================================\n");
    fmt::println(
        "Key insight: A LINE has no unique direction, but the BIVECTOR p∧q does!");
    fmt::println(
        "The bivector orientation depends on the ORDER of points: p∧q = -(q∧p)\n");
    fmt::println("When we form bvt = p∧q with p at origin and q rotating CCW:");
    fmt::println("  - The bivector normal rotates with q");
    fmt::println(
        "  - atan2(-bvt.y, bvt.x) gives the perpendicular direction to the line");
    fmt::println(
        "  - This angle uniquely identifies which boundary-based 'quadrant' we're in\n");

    fmt::println("{:>10} {:>15} {:>12} {:>12} {:>30}", "q angle", "line_angle", "bvt.x",
                 "bvt.y", "Case (boundary-based)");
    fmt::println("{:-^85}", "");

    for (int angle_deg = 0; angle_deg <= 360; angle_deg += 15) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // Determine case based on line_angle
        std::string case_name;
        if (std::abs(line_angle + pi / 2.0) < eps) {
            case_name = "pos x-axis";
        }
        else if (line_angle < -pi / 2.0 - eps && line_angle > -pi + eps) {
            case_name = "first quadrant";
        }
        else if (std::abs(std::abs(line_angle) - pi) < eps) {
            case_name = "pos y-axis";
        }
        else if (line_angle > pi / 2.0 + eps && line_angle < pi - eps) {
            case_name = "second quadrant";
        }
        else if (std::abs(line_angle - pi / 2.0) < eps) {
            case_name = "neg x-axis";
        }
        else if (line_angle > eps && line_angle < pi / 2.0 - eps) {
            case_name = "third quadrant";
        }
        else if (std::abs(line_angle) < eps) {
            case_name = "neg y-axis";
        }
        else if (line_angle > -pi / 2.0 + eps && line_angle < -eps) {
            case_name = "fourth quadrant";
        }
        else {
            case_name = "UNKNOWN";
        }

        fmt::println("{:10.0f}° {:+15.1f}° {:+12.2f} {:+12.2f}   {:30}",
                     (double)angle_deg, line_angle * 180.0 / pi, bvt.x, bvt.y, case_name);
    }

    fmt::println("\n\nSYSTEMATIC MAPPING:");
    fmt::println("===================\n");

    fmt::println("Axes (4 distinct line directions, but 2 bivector orientations each):");
    fmt::println("---------------------------------------------------------------------");
    fmt::println("  Horizontal line (x-axis):");
    fmt::println(
        "    • line_angle = -90° (-π/2) → bivector points along +x → 'pos x-axis'");
    fmt::println(
        "    • line_angle = +90° (+π/2) → bivector points along -x → 'neg x-axis'");
    fmt::println("  Vertical line (y-axis):");
    fmt::println(
        "    • line_angle = ±180° (±π)  → bivector points along +y → 'pos y-axis'");
    fmt::println(
        "    • line_angle = 0°          → bivector points along -y → 'neg y-axis'\n");

    fmt::println("Quadrants (4 distinct line orientations):");
    fmt::println("------------------------------------------");
    fmt::println("  'first quadrant':  line_angle ∈ (-180°, -90°) = (-π, -π/2)");
    fmt::println("    → Line crosses bottom-left to top-right");
    fmt::println("    → Bivector normal points into 2nd classical quadrant\n");

    fmt::println("  'second quadrant': line_angle ∈ (+90°, +180°) = (+π/2, +π)");
    fmt::println("    → Line crosses bottom-right to top-left");
    fmt::println("    → Bivector normal points into 2nd classical quadrant\n");

    fmt::println("  'third quadrant':  line_angle ∈ (0°, +90°) = (0, +π/2)");
    fmt::println("    → Line crosses top-left to bottom-right");
    fmt::println("    → Bivector normal points into 1st classical quadrant\n");

    fmt::println("  'fourth quadrant': line_angle ∈ (-90°, 0°) = (-π/2, 0)");
    fmt::println("    → Line crosses top-right to bottom-left");
    fmt::println("    → Bivector normal points into 4th classical quadrant\n");

    fmt::println("\nCRITICAL UNDERSTANDING:");
    fmt::println("=======================");
    fmt::println("• The GEOMETRIC LINE (as a set of points) has no direction");
    fmt::println("• The BIVECTOR p∧q encodes both the line AND an orientation");
    fmt::println("• Swapping p and q flips the bivector: p∧q = -(q∧p)");
    fmt::println(
        "• atan2(-bvt.y, bvt.x) gives a unique angle in [-π, π] for the bivector");
    fmt::println("• This angle determines which viewport boundaries the line crosses");
    fmt::println("• The 'quadrants' are NOT classical quadrants - they're boundary "
                 "crossing patterns!");

    return 0;
}
