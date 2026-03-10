// Complete table of phi_x, phi_y, and line_angle throughout full rotation

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

    fmt::println("COMPLETE ANGLE TABLE (p at origin, q rotating CCW)");
    fmt::println("==================================================\n");
    fmt::println(
        "Testing if line_angle provides a UNIQUE angle for each line direction.\n");

    fmt::println("{:>10} {:>10} {:>10} {:>12} {:>12} {:>12} {:>30}", "q angle", "phi_x",
                 "phi_y", "line_angle", "bvt.x", "bvt.y", "Case");
    fmt::println("{:-^110}", "");

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

        // Determine case
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

        fmt::println(
            "{:10.0f}° {:10.1f}° {:10.1f}° {:+12.1f}° {:+12.2f} {:+12.2f}   {:30}",
            (double)angle_deg, phi_x * 180.0 / pi, phi_y * 180.0 / pi,
            line_angle * 180.0 / pi, bvt.x, bvt.y, case_name);
    }

    fmt::println("\n\nANALYSIS:");
    fmt::println("=========\n");

    fmt::println("Question: Does a geometric line have a UNIQUE angle?\n");

    fmt::println("GEOMETRIC LINE DIRECTIONS (ignoring orientation):");
    fmt::println("  • Horizontal line:  appears at q=0° and q=180°");
    fmt::println("    - At q=0°:   phi_x=0°,   phi_y=90°,  line_angle=-90°");
    fmt::println("    - At q=180°: phi_x=180°, phi_y=90°,  line_angle=+90°");
    fmt::println("    → Same LINE, different BIVECTOR orientations!\n");

    fmt::println("  • Vertical line:    appears at q=90° and q=270°");
    fmt::println("    - At q=90°:  phi_x=90°,  phi_y=180°, line_angle=±180°");
    fmt::println("    - At q=270°: phi_x=90°,  phi_y=0°,   line_angle=0°");
    fmt::println("    → Same LINE, different BIVECTOR orientations!\n");

    fmt::println("  • Diagonal 45°:     appears at q=45° and q=225°");
    fmt::println("    - At q=45°:  phi_x=45°,  phi_y=135°, line_angle=-135°");
    fmt::println("    - At q=225°: phi_x=135°, phi_y=45°,  line_angle=+45°");
    fmt::println("    → Same LINE, different BIVECTOR orientations!\n");

    fmt::println("CONCLUSION:");
    fmt::println("-----------");
    fmt::println("• A GEOMETRIC LINE does NOT have a unique angle (it's undirected)");
    fmt::println("• A BIVECTOR p∧q DOES have a unique angle (it encodes orientation)");
    fmt::println("• The original code uses phi_x and phi_y (both in [0,π])");
    fmt::println("• The new code uses line_angle (in [-π,π])");
    fmt::println("• Both encode the SAME information: bivector orientation");
    fmt::println("• The mapping IS unique and correct!");

    fmt::println("\n\nQUESTION: Can we use a simpler angle representation?");
    fmt::println("======================================================\n");

    fmt::println("Option 1: Use line_angle ∈ [-π, π] (current implementation)");
    fmt::println("  ✓ Single atan2 call");
    fmt::println("  ✓ Unique for each bivector orientation");
    fmt::println("  ✓ Maps cleanly to boundary crossing patterns");
    fmt::println("  ✓ Already implemented and verified\n");

    fmt::println("Option 2: Use a single angle ∈ [0, π] for undirected lines");
    fmt::println("  × Loses orientation information needed for arrow direction");
    fmt::println("  × Cannot distinguish p∧q from q∧p");
    fmt::println("  × Would require additional orientation flag\n");

    fmt::println("RECOMMENDATION: Keep line_angle ∈ [-π, π] - it's the right choice!");

    return 0;
}
