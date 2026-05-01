// Cross-check exporter for the Python wrapper.
//
// Builds a JSON file of (op, inputs, expected output) tuples by calling the
// C++ library directly. The Python test suite (ga_py/tests/test_cross_check.py)
// loads the JSON, reconstructs the inputs as Python objects, calls the same
// operation through the bindings, and verifies bit-for-bit (within tolerance)
// equivalence.
//
// Output: written to <project_root>/ga_py/tests/data/ga_test_cases.json
//   (path can be overridden by the first command-line argument).
//
// The case set is **representative**, not exhaustive: one case per
// (algebra, op, type-pair) combination, drawn from the existing
// hand-written values in ga_test/src/ga_*_test.hpp. By construction the
// "expected" field is what the C++ library currently produces; if a
// future C++ change alters semantics, regenerate this JSON and the
// Python tests will follow.
//
// Build target: ga_export_python_cases
// Regenerate via:
//   cd build && cmake --build . --target ga_export_python_cases \
//     && ./ga_test/ga_export_python_cases

#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

#include <cstdio>
#include <fmt/format.h>
#include <fstream>
#include <string>
#include <vector>

using namespace hd::ga;
using namespace hd::ga::ega;
using namespace hd::ga::pga;

namespace {

// --------------------------------------------------------------------------- //
// JSON serialization helpers
// --------------------------------------------------------------------------- //

// One field in a Case: a typed value with a list of doubles in canonical
// field order (matching the order the Python wrapper exposes them).
struct JsonValue {
    std::string type;
    std::vector<double> fields;
};

std::string emit_value(JsonValue const& v)
{
    std::string fields;
    for (size_t i = 0; i < v.fields.size(); ++i) {
        if (i) fields += ", ";
        // Use the highest-precision repr that round-trips through double.
        fields += fmt::format("{:.17g}", v.fields[i]);
    }
    return fmt::format(R"({{"type": "{}", "fields": [{}]}})", v.type, fields);
}

// --------------------------------------------------------------------------- //
// to_json overloads --- one per bound user type
// --------------------------------------------------------------------------- //
// Field order MUST match what the Python wrapper exposes. For struct-shaped
// types (vec*, bivec*, trivec*) this is the C++ struct field order. For
// scalar/pscalar it is just the value. For multivectors it is c0..cN.

// Plain numeric scalars (returned by nrm_sq, bulk_nrm_sq, etc.).
// In Python these become a `float`, not a typed scalar wrapper.
[[maybe_unused]] JsonValue to_json(double v) { return {"double", {v}}; }

[[maybe_unused]] JsonValue to_json(scalar2d v) { return {"scalar2d", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(scalar3d v) { return {"scalar3d", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(scalar4d v) { return {"scalar4d", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar2d v) { return {"pscalar2d", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar3d v) { return {"pscalar3d", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar4d v) { return {"pscalar4d", {value_t(v)}}; }

[[maybe_unused]] JsonValue to_json(vec2d const& v) { return {"vec2d", {v.x, v.y}}; }
[[maybe_unused]] JsonValue to_json(vec3d const& v) { return {"vec3d", {v.x, v.y, v.z}}; }
[[maybe_unused]] JsonValue to_json(vec4d const& v) { return {"vec4d", {v.x, v.y, v.z, v.w}}; }
[[maybe_unused]] JsonValue to_json(bivec3d const& B) { return {"bivec3d", {B.x, B.y, B.z}}; }
[[maybe_unused]] JsonValue to_json(bivec4d const& B)
{
    return {"bivec4d", {B.vx, B.vy, B.vz, B.mx, B.my, B.mz}};
}
[[maybe_unused]] JsonValue to_json(trivec4d const& t) { return {"trivec4d", {t.x, t.y, t.z, t.w}}; }

[[maybe_unused]] JsonValue to_json(scalar2dp v) { return {"scalar2dp", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(scalar3dp v) { return {"scalar3dp", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar2dp v) { return {"pscalar2dp", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar3dp v) { return {"pscalar3dp", {value_t(v)}}; }

[[maybe_unused]] JsonValue to_json(vec2dp const& v) { return {"vec2dp", {v.x, v.y, v.z}}; }
[[maybe_unused]] JsonValue to_json(vec3dp const& v) { return {"vec3dp", {v.x, v.y, v.z, v.w}}; }
[[maybe_unused]] JsonValue to_json(bivec2dp const& B) { return {"bivec2dp", {B.x, B.y, B.z}}; }
[[maybe_unused]] JsonValue to_json(bivec3dp const& B)
{
    return {"bivec3dp", {B.vx, B.vy, B.vz, B.mx, B.my, B.mz}};
}
[[maybe_unused]] JsonValue to_json(trivec3dp const& t) { return {"trivec3dp", {t.x, t.y, t.z, t.w}}; }

// Multivectors: emit c0..cN in component order.
[[maybe_unused]] JsonValue to_json(MVec2d<value_t> const& M)
{
    return {"mvec2d", {M.c0, M.c1, M.c2, M.c3}};
}
[[maybe_unused]] JsonValue to_json(MVec2d_E<value_t> const& M) { return {"mvec2d_e", {M.c0, M.c1}}; }
[[maybe_unused]] JsonValue to_json(MVec3d<value_t> const& M)
{
    return {"mvec3d", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}
[[maybe_unused]] JsonValue to_json(MVec3d_E<value_t> const& M)
{
    return {"mvec3d_e", {M.c0, M.c1, M.c2, M.c3}};
}
[[maybe_unused]] JsonValue to_json(MVec3d_U<value_t> const& M)
{
    return {"mvec3d_u", {M.c0, M.c1, M.c2, M.c3}};
}
[[maybe_unused]] JsonValue to_json(MVec2dp<value_t> const& M)
{
    return {"mvec2dp", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}
[[maybe_unused]] JsonValue to_json(MVec2dp_E<value_t> const& M)
{
    return {"mvec2dp_e", {M.c0, M.c1, M.c2, M.c3}};
}
[[maybe_unused]] JsonValue to_json(MVec2dp_U<value_t> const& M)
{
    return {"mvec2dp_u", {M.c0, M.c1, M.c2, M.c3}};
}
[[maybe_unused]] JsonValue to_json(MVec3dp<value_t> const& M)
{
    return {"mvec3dp",
            {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11,
             M.c12, M.c13, M.c14, M.c15}};
}
[[maybe_unused]] JsonValue to_json(MVec3dp_E<value_t> const& M)
{
    return {"mvec3dp_e", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}
[[maybe_unused]] JsonValue to_json(MVec3dp_U<value_t> const& M)
{
    return {"mvec3dp_u", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}

[[maybe_unused]] JsonValue to_json(DualNum2dp<value_t> const& d) { return {"dualnum2dp", {d.c0, d.c1}}; }
[[maybe_unused]] JsonValue to_json(DualNum3dp<value_t> const& d) { return {"dualnum3dp", {d.c0, d.c1}}; }

// --------------------------------------------------------------------------- //
// Case collection
// --------------------------------------------------------------------------- //

struct Case {
    std::string id;
    std::string submodule; // "ega" | "pga" | "top"
    std::string op;
    std::vector<JsonValue> args;
    JsonValue expected;
    double tolerance = 1e-12;
};

std::vector<Case> g_cases;

// Helper: emit a case, computing the result via a lambda.
template <typename ResultT, typename... Args>
void add_case(std::string const& id, std::string const& submodule, std::string const& op,
              ResultT const& result, double tolerance, Args const&... args)
{
    Case c;
    c.id = id;
    c.submodule = submodule;
    c.op = op;
    c.tolerance = tolerance;
    (c.args.push_back(to_json(args)), ...);
    c.expected = to_json(result);
    g_cases.push_back(std::move(c));
}

// Convenience wrapper with default tolerance.
template <typename ResultT, typename... Args>
void add(std::string const& id, std::string const& submodule, std::string const& op,
         ResultT const& result, Args const&... args)
{
    add_case(id, submodule, op, result, 1e-12, args...);
}

// --------------------------------------------------------------------------- //
// EGA cases
// --------------------------------------------------------------------------- //

void emit_ega_cases()
{
    // --- vec2d ---
    {
        vec2d a{1.0, 2.0}, b{3.0, -4.0};
        add("ega_vec2d_dot", "ega", "dot", dot(a, b), a, b);
        add("ega_vec2d_wdg", "ega", "wdg", wdg(a, b), a, b);
        add("ega_vec2d_add", "ega", "+", a + b, a, b);
        add("ega_vec2d_sub", "ega", "-", a - b, a, b);
        add("ega_vec2d_neg", "ega", "neg", -a, a);
        add("ega_vec2d_nrm_sq", "ega", "nrm_sq", nrm_sq(a), a);
    }

    // --- vec3d ---
    {
        vec3d a{1.0, 2.0, 3.0}, b{4.0, 5.0, 6.0};
        add("ega_vec3d_dot", "ega", "dot", dot(a, b), a, b);
        add("ega_vec3d_wdg", "ega", "wdg", wdg(a, b), a, b);
        add("ega_vec3d_cross", "ega", "cross", cross(a, b), a, b);
        add("ega_vec3d_add", "ega", "+", a + b, a, b);
        add("ega_vec3d_sub", "ega", "-", a - b, a, b);
        add("ega_vec3d_neg", "ega", "neg", -a, a);
        add("ega_vec3d_nrm_sq", "ega", "nrm_sq", nrm_sq(a), a);
        add("ega_vec3d_rev", "ega", "rev", rev(a), a);
        // Geometric product of two vectors -> mvec3d_e (scalar + bivector parts)
        add("ega_vec3d_gpr", "ega", "*", a * b, a, b);
        // Inverse: vec3d -> vec3d
        vec3d c{2.0, 6.0, -4.0};
        add("ega_vec3d_inv", "ega", "inv", inv(c), c);
    }

    // --- bivec3d ---
    {
        bivec3d B1{1.0, -2.0, 3.0}, B2{0.5, 1.5, -0.5};
        add("ega_bivec3d_dot", "ega", "dot", dot(B1, B2), B1, B2);
        add("ega_bivec3d_add", "ega", "+", B1 + B2, B1, B2);
        add("ega_bivec3d_neg", "ega", "neg", -B1, B1);
        add("ega_bivec3d_nrm_sq", "ega", "nrm_sq", nrm_sq(B1), B1);
        add("ega_bivec3d_rev", "ega", "rev", rev(B1), B1);
        // bivec * bivec -> mvec3d_e (scalar + bivector)
        add("ega_bivec3d_gpr", "ega", "*", B1 * B2, B1, B2);
        add("ega_bivec3d_cmpl", "ega", "cmpl", cmpl(B1), B1);
    }

    // --- complement / dual ---
    {
        vec3d v{1.0, 2.0, 3.0};
        add("ega_vec3d_cmpl", "ega", "cmpl", cmpl(v), v);
        add("ega_vec3d_dual", "ega", "dual", dual(v), v);
    }

    // --- rotor / rotation ---
    // Rotation of e1 by 90° in the e12 plane around the origin -> e2.
    {
        auto R = get_rotor(e12_3d, deg2rad(90.0));
        vec3d v{1.0, 0.0, 0.0};
        add("ega_rotate_e1_90deg", "ega", "rotate", rotate(v, R), v, R);
    }
    // Rotor itself (mvec3d_e) for a known angle.
    {
        auto R = get_rotor(e12_3d, deg2rad(45.0));
        // get_rotor takes (BiVec3d, value_t) — but value_t isn't a JSON-emittable
        // GA type; emit the rotor's components only as the expected output of a
        // single-arg rotor-info case. Bind to a representative call: rev(R).
        add("ega_rotor_45deg_rev", "ega", "rev", rev(R), R);
    }

    // --- mvec3d operations ---
    {
        mvec3d M1{scalar3d{1.0}, vec3d{2.0, 0.0, 0.0}, bivec3d{0.0, 0.0, 3.0},
                  pscalar3d{4.0}};
        add("ega_mvec3d_rev", "ega", "rev", rev(M1), M1);
        add("ega_mvec3d_neg", "ega", "neg", -M1, M1);
    }
}

// --------------------------------------------------------------------------- //
// PGA cases
// --------------------------------------------------------------------------- //

void emit_pga_cases()
{
    // --- vec2dp / vec3dp basic algebra ---
    {
        vec2dp a{1.0, 2.0, 1.0}, b{3.0, 4.0, 1.0};
        add("pga_vec2dp_dot", "pga", "dot", dot(a, b), a, b);
        add("pga_vec2dp_wdg", "pga", "wdg", wdg(a, b), a, b);
        add("pga_vec2dp_add", "pga", "+", a + b, a, b);
        add("pga_vec2dp_neg", "pga", "neg", -a, a);
        add("pga_vec2dp_unitize", "pga", "unitize", unitize(a), a);
        add("pga_vec2dp_bulk", "pga", "bulk", bulk(a), a);
        add("pga_vec2dp_weight", "pga", "weight", weight(a), a);
        // PGA2dp is odd-dimensional, so plain bulk_dual / weight_dual are
        // defined (no l_/r_ prefix needed).
        add("pga_vec2dp_bulk_dual", "pga", "bulk_dual", bulk_dual(a), a);
        add("pga_vec2dp_weight_dual", "pga", "weight_dual", weight_dual(a), a);
        // PGA2dp is odd-dimensional, so plain cmpl is defined (no l_/r_ split).
        add("pga_vec2dp_cmpl", "pga", "cmpl", cmpl(a), a);
    }

    {
        vec3dp p{1.0, 2.0, 3.0, 1.0}, q{4.0, 5.0, 6.0, 1.0};
        add("pga_vec3dp_dot", "pga", "dot", dot(p, q), p, q);
        add("pga_vec3dp_wdg", "pga", "wdg", wdg(p, q), p, q);
        add("pga_vec3dp_join", "pga", "join", join(p, q), p, q);
        add("pga_vec3dp_add", "pga", "+", p + q, p, q);
        add("pga_vec3dp_neg", "pga", "neg", -p, p);
        add("pga_vec3dp_unitize", "pga", "unitize", unitize(p), p);
        add("pga_vec3dp_bulk", "pga", "bulk", bulk(p), p);
        add("pga_vec3dp_weight", "pga", "weight", weight(p), p);
        add("pga_vec3dp_att", "pga", "att", att(p), p);
        add("pga_vec3dp_rev", "pga", "rev", rev(p), p);
        add("pga_vec3dp_r_cmpl", "pga", "r_cmpl", r_cmpl(p), p);
        add("pga_vec3dp_l_cmpl", "pga", "l_cmpl", l_cmpl(p), p);
        // PGA3dp is even-dimensional, so the dual splits and there is no plain
        // `bulk_dual` / `weight_dual` --- you must pick l_ or r_ explicitly.
        add("pga_vec3dp_r_bulk_dual", "pga", "r_bulk_dual", r_bulk_dual(p), p);
        add("pga_vec3dp_l_bulk_dual", "pga", "l_bulk_dual", l_bulk_dual(p), p);
        add("pga_vec3dp_r_weight_dual", "pga", "r_weight_dual", r_weight_dual(p), p);
        add("pga_vec3dp_l_weight_dual", "pga", "l_weight_dual", l_weight_dual(p), p);
        add("pga_vec3dp_bulk_nrm_sq", "pga", "bulk_nrm_sq", bulk_nrm_sq(p), p);
        add("pga_vec3dp_weight_nrm_sq", "pga", "weight_nrm_sq", weight_nrm_sq(p), p);
    }

    // --- bivec3dp (line) ---
    {
        bivec3dp L{1.0, 2.0, -1.0, 3.0, -2.0, 4.0};
        add("pga_bivec3dp_neg", "pga", "neg", -L, L);
        add("pga_bivec3dp_rev", "pga", "rev", rev(L), L);
        add("pga_bivec3dp_rrev", "pga", "rrev", rrev(L), L);
        add("pga_bivec3dp_bulk", "pga", "bulk", bulk(L), L);
        add("pga_bivec3dp_weight", "pga", "weight", weight(L), L);
        add("pga_bivec3dp_att", "pga", "att", att(L), L);
        add("pga_bivec3dp_r_cmpl", "pga", "r_cmpl", r_cmpl(L), L);
        add("pga_bivec3dp_bulk_nrm_sq", "pga", "bulk_nrm_sq", bulk_nrm_sq(L), L);
        add("pga_bivec3dp_weight_nrm_sq", "pga", "weight_nrm_sq", weight_nrm_sq(L), L);
    }

    // --- trivec3dp (plane) ---
    {
        trivec3dp T{2.0, -1.0, 3.0, 4.0};
        add("pga_trivec3dp_neg", "pga", "neg", -T, T);
        add("pga_trivec3dp_rev", "pga", "rev", rev(T), T);
        add("pga_trivec3dp_bulk", "pga", "bulk", bulk(T), T);
        add("pga_trivec3dp_weight", "pga", "weight", weight(T), T);
        add("pga_trivec3dp_att", "pga", "att", att(T), T);
        add("pga_trivec3dp_r_cmpl", "pga", "r_cmpl", r_cmpl(T), T);
    }

    // --- distance + dualnum ---
    {
        vec3dp p1{0.0, 0.0, 0.0, 1.0};
        vec3dp p2{1.0, 0.0, 0.0, 1.0}; // unit distance along x
        add("pga_dist3dp_unit", "pga", "dist3dp", dist3dp(p1, p2), p1, p2);
    }

    // --- meet / join chain (line ∩ plane = point) ---
    {
        vec3dp p1{2.0, 0.0, 1.0, 1.0};
        vec3dp p2{4.0, 3.0, 1.0, 1.0};
        auto l12 = join(p1, p2);
        vec3dp p3{0.0, 6.0, 1.0, 1.0};
        auto plane = unitize(join(l12, p3));
        add("pga_join_points_to_line", "pga", "join", l12, p1, p2);
        add("pga_meet_line_with_plane", "pga", "meet", meet(l12, plane), l12, plane);
    }

    // --- motor build + sandwich (uses rgpr; emits the produced multivector) ---
    {
        auto p0 = O_3dp;
        auto p1 = vec3dp{1.0, 0.0, 0.0, 1.0};
        auto p2 = vec3dp{1.0, 1.0, 0.0, 1.0};
        auto p3 = vec3dp{0.0, 0.0, 1.0, 1.0};
        auto pl1 = unitize(wdg(wdg(p0, p3), p1));
        auto pl2 = unitize(wdg(wdg(p0, p3), p2));
        auto motor = rgpr(pl2, pl1);
        add("pga_motor_from_two_planes", "pga", "rgpr", motor, pl2, pl1);
        // move3dp(point, motor) — sandwich
        vec3dp pt{1.0, -0.5, 0.0, 1.0};
        add_case("pga_move3dp_point", "pga", "move3dp", move3dp(pt, motor), 1e-9, pt,
                 motor);
    }
}

// --------------------------------------------------------------------------- //
// Top-level cases (functions in hd::ga)
// --------------------------------------------------------------------------- //

void emit_top_cases()
{
    // Top-level numeric helpers take and return `value_t` (plain double).
    add("top_deg2rad_180", "top", "deg2rad", deg2rad(180.0), 180.0);
    add("top_rad2deg_pi", "top", "rad2deg", rad2deg(3.141592653589793),
        3.141592653589793);
    add("top_Hz2radps_1", "top", "Hz2radps", Hz2radps(1.0), 1.0);
    add("top_radps2Hz_1", "top", "radps2Hz", radps2Hz(1.0), 1.0);
    add("top_sign_neg",   "top", "sign",     static_cast<double>(sign(-3.5)), -3.5);
}

// --------------------------------------------------------------------------- //
// JSON file emission
// --------------------------------------------------------------------------- //

void write_json(std::string const& path)
{
    std::ofstream out(path);
    if (!out) {
        fmt::print(stderr, "FAILED to open {} for writing\n", path);
        std::exit(1);
    }
    out << "{\n";
    out << R"(  "version": 1,)" << "\n";
    out << R"(  "note": "AUTO-GENERATED by ga_test/src/export_python_cases.cpp -- do not edit by hand. Regenerate via: cmake --build . --target ga_export_python_cases && ./ga_test/ga_export_python_cases",)"
        << "\n";
    out << R"(  "cases": [)" << "\n";
    for (size_t i = 0; i < g_cases.size(); ++i) {
        auto const& c = g_cases[i];
        out << "    {\n";
        out << fmt::format(R"(      "id": "{}",)", c.id) << "\n";
        out << fmt::format(R"(      "submodule": "{}",)", c.submodule) << "\n";
        out << fmt::format(R"(      "op": "{}",)", c.op) << "\n";
        out << R"(      "args": [)" << "\n";
        for (size_t j = 0; j < c.args.size(); ++j) {
            out << "        " << emit_value(c.args[j]);
            out << (j + 1 == c.args.size() ? "\n" : ",\n");
        }
        out << "      ],\n";
        out << "      \"expected\": " << emit_value(c.expected) << ",\n";
        out << fmt::format(R"(      "tolerance": {:.0e})", c.tolerance) << "\n";
        out << "    }" << (i + 1 == g_cases.size() ? "\n" : ",\n");
    }
    out << "  ]\n";
    out << "}\n";
}

} // anonymous namespace

int main(int argc, char** argv)
{
    std::string path = "ga_py/tests/data/ga_test_cases.json";
    if (argc > 1) path = argv[1];

    emit_ega_cases();
    emit_pga_cases();
    emit_top_cases();

    write_json(path);
    fmt::print("Wrote {} cases to {}\n", g_cases.size(), path);
    return 0;
}
