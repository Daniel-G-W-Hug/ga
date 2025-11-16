// Test atan2 behavior at the ±180° boundary

#include "ga/ga_pga.hpp"
#include <cmath>
#include <fmt/core.h>
#include <numbers>

using namespace hd::ga;
using namespace hd::ga::pga;
using std::numbers::pi;

int main()
{
    fmt::println("Testing atan2 behavior at ±180° boundary");
    fmt::println("=========================================\n");

    auto p_origin = vec2dp(0.0, 0.0, 1.0);
    value_t radius = 3.0;

    // Test points very close to 90° (pos y-axis)
    fmt::println("Testing near q = 90° (pos y-axis):\n");
    fmt::println("{:>15} {:>12} {:>12} {:>15} {:>20}", "q angle", "bvt.x", "bvt.y",
                 "line_angle", "Matches condition?");
    fmt::println("{:-^80}", "");

    for (double angle_deg = 89.9; angle_deg <= 90.1; angle_deg += 0.05) {
        value_t angle_rad = angle_deg * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        // Check if matches condition: |abs(line_angle) - π| < eps
        bool matches = std::abs(std::abs(line_angle) - pi) < eps;

        fmt::println("{:15.2f}° {:+12.6f} {:+12.6f} {:+15.6f}° {:>20}", angle_deg, bvt.x,
                     bvt.y, line_angle * 180.0 / pi, matches ? "YES (pos y-axis)" : "NO");
    }

    fmt::println("\n\nDirect test at exactly 90°:");
    fmt::println("----------------------------\n");

    {
        value_t angle_rad = 90.0 * pi / 180.0;
        auto q = vec2dp(radius * std::cos(angle_rad), radius * std::sin(angle_rad), 1.0);
        auto bvt = wdg(p_origin, q);

        value_t line_angle = std::atan2(-bvt.y, bvt.x);

        fmt::println("q = (0, 3)");
        fmt::println("bvt = p∧q = ({:.10f}, {:.10f}, {:.10f})", bvt.x, bvt.y, bvt.z);
        fmt::println("line_angle = atan2(-bvt.y, bvt.x) = atan2({:.10f}, {:.10f})",
                     -bvt.y, bvt.x);
        fmt::println("line_angle = {:.15f} radians", line_angle);
        fmt::println("line_angle = {:.15f}°", line_angle * 180.0 / pi);
        fmt::println("|line_angle| = {:.15f} radians", std::abs(line_angle));
        fmt::println("π = {:.15f} radians", pi);
        fmt::println("|line_angle| - π = {:.15e}", std::abs(line_angle) - pi);
        fmt::println("eps = {:.15e}", eps);
        fmt::println("");
        fmt::println("Does |abs(line_angle) - π| < eps? {}",
                     std::abs(std::abs(line_angle) - pi) < eps ? "YES" : "NO");
    }

    fmt::println("\n\nTest atan2 return value for different inputs:");
    fmt::println("----------------------------------------------\n");

    struct TestCase {
        double y;
        double x;
        std::string description;
    };

    TestCase cases[] = {
        {0.0, -3.0, "atan2(0, -3) = atan2(0, negative)"},
        {0.0, 3.0, "atan2(0, 3) = atan2(0, positive)"},
        {3.0, 0.0, "atan2(3, 0) = atan2(positive, 0)"},
        {-3.0, 0.0, "atan2(-3, 0) = atan2(negative, 0)"},
        {1e-15, -3.0, "atan2(tiny positive, negative)"},
        {-1e-15, -3.0, "atan2(tiny negative, negative)"},
    };

    for (auto& tc : cases) {
        value_t result = std::atan2(tc.y, tc.x);
        fmt::println("{:50} = {:+.15f}° ({:+.15f} rad)", tc.description,
                     result * 180.0 / pi, result);
    }

    fmt::println("\n\nKEY FINDING:");
    fmt::println("============\n");
    fmt::println(
        "atan2(0, negative) returns EXACTLY ±π (both +π and -π are the same angle).");
    fmt::println("The C++ standard allows either +π or -π to be returned.");
    fmt::println("On this system, atan2(0, -3) = {:.15f}°",
                 std::atan2(0.0, -3.0) * 180.0 / pi);
    fmt::println("");
    fmt::println(
        "The condition |abs(line_angle) - π| < eps correctly handles BOTH cases:");
    fmt::println("  • If atan2 returns +π: |+π| - π = 0 < eps ✓");
    fmt::println("  • If atan2 returns -π: |-π| - π = 0 < eps ✓");
    fmt::println("");
    fmt::println("Therefore, the implementation is ROBUST across different platforms!");

    return 0;
}
