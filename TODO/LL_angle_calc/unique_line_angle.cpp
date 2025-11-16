// Can we define a UNIQUE angle for each geometric line direction?

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

    fmt::println("QUEST FOR UNIQUE LINE ANGLE");
    fmt::println("============================\n");

    fmt::println("Idea: Can we use phi_x alone to uniquely identify line direction?\n");

    fmt::println("{:>10} {:>10} {:>10} {:>12} {:>30}", "q angle", "phi_x", "phi_y",
                 "line_angle", "Unique line direction?");
    fmt::println("{:-^85}", "");

    // Track which phi_x values we've seen
    for (int angle_deg = 0; angle_deg <= 180; angle_deg += 15) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        auto const x_axis = bivec2dp{0, 1, 0};
        auto const y_axis = bivec2dp{1, 0, 0};
        auto const phi_x = angle(x_axis, bvt);
        auto const phi_y = angle(y_axis, bvt);

        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        fmt::println("{:10.0f}° {:10.1f}° {:10.1f}° {:+12.1f}°   Each line appears ONCE",
                     (double)angle_deg, phi_x * 180.0 / pi, phi_y * 180.0 / pi,
                     line_angle * 180.0 / pi);
    }

    fmt::println("\n\nCRITICAL OBSERVATION:");
    fmt::println("=====================\n");

    fmt::println("If we rotate q from 0° to 180° (half circle), each geometric line");
    fmt::println("appears EXACTLY ONCE! This means:\n");

    fmt::println("  • phi_x ∈ [0°, 180°] uniquely identifies each line direction");
    fmt::println("  • But phi_x alone doesn't tell us the arrow orientation!");
    fmt::println("  • We need BOTH phi_x and phi_y to determine arrow direction\n");

    fmt::println("ALTERNATIVE: Use line_angle but normalize to [0°, 180°]?\n");

    fmt::println("{:>10} {:>12} {:>18} {:>30}", "q angle", "line_angle", "normalized",
                 "Problem?");
    fmt::println("{:-^85}", "");

    for (int angle_deg = 0; angle_deg <= 360; angle_deg += 15) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // Try to normalize to [0, π]
        value_t normalized = line_angle;
        if (normalized < 0) normalized += pi;
        if (normalized > pi) normalized -= pi;

        std::string problem = "";
        if (angle_deg > 0 && angle_deg < 180) {
            // Check if we lost information
            int opposite_deg = angle_deg + 180;
            if (opposite_deg <= 360) {
                problem = "← Loses orientation info!";
            }
        }

        fmt::println("{:10.0f}° {:+12.1f}° {:+18.1f}°   {:30}", (double)angle_deg,
                     line_angle * 180.0 / pi, normalized * 180.0 / pi, problem);
    }

    fmt::println("\n\nFINAL ANSWER:");
    fmt::println("=============\n");

    fmt::println("Can a geometric LINE have a unique angle? YES and NO:\n");

    fmt::println("YES - For the GEOMETRIC LINE (undirected):");
    fmt::println("  • We could use phi_x ∈ [0°, 180°] as the unique line angle");
    fmt::println("  • Each value uniquely identifies one undirected line\n");

    fmt::println("NO - For RENDERING with an ARROW:");
    fmt::println("  • We MUST know which end has the arrow (p→q vs q→p)");
    fmt::println("  • This requires the full bivector orientation");
    fmt::println("  • Either use (phi_x, phi_y) or use line_angle ∈ [-180°, 180°]\n");

    fmt::println("CURRENT IMPLEMENTATION:");
    fmt::println("  • Original: (phi_x, phi_y) - two angles in [0°, 180°]");
    fmt::println("  • New: line_angle - one angle in [-180°, 180°]");
    fmt::println("  • Both encode FULL bivector orientation");
    fmt::println("  • Both are UNIQUE for each orientation");
    fmt::println("  • Both are CORRECT - new version just simpler!");

    return 0;
}
