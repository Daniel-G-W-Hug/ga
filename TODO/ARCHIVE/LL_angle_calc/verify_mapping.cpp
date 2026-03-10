// Verify that atan2 mapping matches original phi_x/phi_y logic

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

    fmt::println("Verifying quadrant mapping:");
    fmt::println("{:>5} {:>10} {:>10} {:>25} {:>25}", "angle", "phi_x", "phi_y",
                 "Original Logic", "New atan2 Logic");
    fmt::println("{:-^90}", "");

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

        // Determine which case in ORIGINAL logic
        std::string orig_case = "UNKNOWN";

        // pos. x-axis
        if ((std::abs(phi_x) < eps) && (std::abs(phi_y - pi / 2.0) < eps)) {
            orig_case = "pos x-axis";
        }
        // first quadrant
        else if (phi_x > eps && (phi_x < (pi / 2. - eps)) && (phi_y > (pi / 2. + eps)) &&
                 (phi_y < (pi - eps))) {
            orig_case = "first quadrant";
        }
        // pos. y-axis
        else if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y - pi) < eps)) {
            orig_case = "pos y-axis";
        }
        // second quadrant
        else if ((phi_x > (pi / 2.0 + eps)) && (phi_x < (pi - eps)) &&
                 (phi_y < (pi - eps)) && (phi_y > (pi / 2.0 + eps))) {
            orig_case = "second quadrant";
        }
        // neg. x-axis
        else if ((std::abs(phi_x - pi) < eps) && (std::abs(phi_y - pi / 2.0)) < eps) {
            orig_case = "neg x-axis";
        }
        // third quadrant
        else if ((phi_x < pi - eps) && (phi_x > pi / 2.0 + eps) &&
                 (phi_y < pi / 2.0 - eps) && (phi_y > eps)) {
            orig_case = "third quadrant";
        }
        // neg. y-axis
        else if ((std::abs(phi_x - pi / 2.0) < eps) && (std::abs(phi_y) < eps)) {
            orig_case = "neg y-axis";
        }
        // fourth quadrant
        else if ((phi_x < pi / 2.0 - eps) && (phi_x > eps) && (phi_y > eps) &&
                 (phi_y < pi / 2. - eps)) {
            orig_case = "fourth quadrant";
        }

        // Determine which case in NEW logic
        std::string new_case = "UNKNOWN";

        if (std::abs(line_angle + pi / 2.0) < eps) {
            new_case = "pos x-axis";
        }
        else if (line_angle < -pi / 2.0 - eps && line_angle > -pi + eps) {
            new_case = "first quadrant";
        }
        else if (std::abs(std::abs(line_angle) - pi) < eps) {
            new_case = "pos y-axis";
        }
        else if (line_angle > pi / 2.0 + eps && line_angle < pi - eps) {
            new_case = "second quadrant";
        }
        else if (std::abs(line_angle - pi / 2.0) < eps) {
            new_case = "neg x-axis";
        }
        else if (line_angle > eps && line_angle < pi / 2.0 - eps) {
            new_case = "third quadrant";
        }
        else if (std::abs(line_angle) < eps) {
            new_case = "neg y-axis";
        }
        else if (line_angle > -pi / 2.0 + eps && line_angle < -eps) {
            new_case = "fourth quadrant";
        }

        bool match = (orig_case == new_case);
        std::string status = match ? "✓" : "✗ MISMATCH!";
        fmt::println("{:5.0f} {:10.1f} {:10.1f} {:>25} {:>25} {:>15}", (double)angle_deg,
                     phi_x * 180.0 / pi, phi_y * 180.0 / pi, orig_case, new_case, status);
    }

    return 0;
}
