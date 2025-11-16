// Map original phi_x/phi_y logic to atan2 ranges

#include "ga/ga_pga.hpp"
#include <fmt/core.h>
#include <numbers>

using namespace hd::ga;
using namespace hd::ga::pga;
using std::numbers::pi;

int main()
{
    // Test various line orientations to understand mapping
    auto p_origin = vec2dp(0.0, 0.0, 1.0);

    struct TestCase {
        double angle_deg;
        std::string expected_case;
    };

    TestCase cases[] = {{0, "pos x-axis"},   {45, "first quadrant"},
                        {90, "pos y-axis"},  {135, "second quadrant"},
                        {180, "neg x-axis"}, {225, "third quadrant"},
                        {270, "neg y-axis"}, {315, "fourth quadrant"}};

    fmt::println("Mapping original quadrant logic to atan2:");
    fmt::println("{:>5} {:>20} {:>10} {:>10} {:>10}", "angle", "expected", "phi_x",
                 "phi_y", "atan2");
    fmt::println("{:-^60}", "");

    for (auto& tc : cases) {
        double angle_rad = tc.angle_deg * pi / 180.0;
        auto q = vec2dp(3.0 * std::cos(angle_rad), 3.0 * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        // Original approach
        auto const x_axis = bivec2dp{0, 1, 0};
        auto const y_axis = bivec2dp{1, 0, 0};
        auto const phi_x = angle(x_axis, bvt);
        auto const phi_y = angle(y_axis, bvt);

        // New approach
        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        fmt::println("{:5.0f} {:>20} {:10.1f} {:10.1f} {:+10.1f}", tc.angle_deg,
                     tc.expected_case, phi_x * 180.0 / pi, phi_y * 180.0 / pi,
                     line_angle * 180.0 / pi);
    }

    return 0;
}
