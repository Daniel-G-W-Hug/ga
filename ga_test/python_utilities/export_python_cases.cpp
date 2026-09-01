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

#include "ga/ga_cga.hpp"
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"
#include "ga/ga_sta.hpp"

#include <cstdio>
#include <fmt/format.h>
#include <fstream>
#include <string>
#include <vector>

using namespace hd::ga;
using namespace hd::ga::ega;
using namespace hd::ga::pga;
using namespace hd::ga::cga;
using namespace hd::ga::sta;

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
[[maybe_unused]] JsonValue to_json(pscalar2d v) { return {"pscalar2d", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar3d v) { return {"pscalar3d", {value_t(v)}}; }

[[maybe_unused]] JsonValue to_json(vec2d const& v) { return {"vec2d", {v.x, v.y}}; }
[[maybe_unused]] JsonValue to_json(vec3d const& v) { return {"vec3d", {v.x, v.y, v.z}}; }

[[maybe_unused]] JsonValue to_json(bivec3d const& B)
{
    return {"bivec3d", {B.x, B.y, B.z}};
}

[[maybe_unused]] JsonValue to_json(scalar2dp v) { return {"scalar2dp", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(scalar3dp v) { return {"scalar3dp", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar2dp v) { return {"pscalar2dp", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar3dp v) { return {"pscalar3dp", {value_t(v)}}; }

[[maybe_unused]] JsonValue to_json(scalar4ds v) { return {"scalar4ds", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar4ds v) { return {"pscalar4ds", {value_t(v)}}; }

[[maybe_unused]] JsonValue to_json(vec2dp const& v)
{
    return {"vec2dp", {v.x, v.y, v.z}};
}
[[maybe_unused]] JsonValue to_json(vec3dp const& v)
{
    return {"vec3dp", {v.x, v.y, v.z, v.w}};
}
[[maybe_unused]] JsonValue to_json(vec4ds const& v)
{
    return {"vec4ds", {v.x, v.y, v.z, v.w}};
}
[[maybe_unused]] JsonValue to_json(bivec2dp const& B)
{
    return {"bivec2dp", {B.x, B.y, B.z}};
}
[[maybe_unused]] JsonValue to_json(bivec3dp const& B)
{
    return {"bivec3dp", {B.vx, B.vy, B.vz, B.mx, B.my, B.mz}};
}
[[maybe_unused]] JsonValue to_json(bivec4ds const& B)
{
    return {"bivec4ds", {B.vx, B.vy, B.vz, B.mx, B.my, B.mz}};
}
[[maybe_unused]] JsonValue to_json(trivec3dp const& t)
{
    return {"trivec3dp", {t.x, t.y, t.z, t.w}};
}
[[maybe_unused]] JsonValue to_json(trivec4ds const& t)
{
    return {"trivec4ds", {t.x, t.y, t.z, t.w}};
}

// Multivectors: emit c0..cN in component order.
[[maybe_unused]] JsonValue to_json(MVec2d<value_t> const& M)
{
    return {"mvec2d", {M.c0, M.c1, M.c2, M.c3}};
}
[[maybe_unused]] JsonValue to_json(MVec2d_E<value_t> const& M)
{
    return {"mvec2d_e", {M.c0, M.c1}};
}
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
[[maybe_unused]] JsonValue to_json(MVec4ds<value_t> const& M)
{
    return {"mvec4ds",
            {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11,
             M.c12, M.c13, M.c14, M.c15}};
}
[[maybe_unused]] JsonValue to_json(MVec4ds_E<value_t> const& M)
{
    return {"mvec4ds_e", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}
[[maybe_unused]] JsonValue to_json(MVec4ds_U<value_t> const& M)
{
    return {"mvec4ds_u", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}

[[maybe_unused]] JsonValue to_json(scalar2dc v) { return {"scalar2dc", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(scalar3dc v) { return {"scalar3dc", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar2dc v) { return {"pscalar2dc", {value_t(v)}}; }
[[maybe_unused]] JsonValue to_json(pscalar3dc v) { return {"pscalar3dc", {value_t(v)}}; }

[[maybe_unused]] JsonValue to_json(vec2dc const& v)
{
    return {"vec2dc", {v.x, v.y, v.z, v.w}};
}
[[maybe_unused]] JsonValue to_json(trivec2dc const& t)
{
    return {"trivec2dc", {t.x, t.y, t.z, t.w}};
}
[[maybe_unused]] JsonValue to_json(bivec2dc const& B)
{
    return {"bivec2dc", {B.vx, B.vy, B.vz, B.mx, B.my, B.mz}};
}
[[maybe_unused]] JsonValue to_json(vec3dc const& v)
{
    return {"vec3dc", {v.x, v.y, v.z, v.w, v.u}};
}
[[maybe_unused]] JsonValue to_json(quadvec3dc const& q)
{
    return {"quadvec3dc", {q.x, q.y, q.z, q.w, q.u}};
}
[[maybe_unused]] JsonValue to_json(bivec3dc const& B)
{
    return {"bivec3dc", {B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py, B.pz, B.pw}};
}
[[maybe_unused]] JsonValue to_json(trivec3dc const& T)
{
    return {"trivec3dc", {T.vx, T.vy, T.vz, T.mx, T.my, T.mz, T.px, T.py, T.pz, T.pw}};
}

[[maybe_unused]] JsonValue to_json(MVec2dc<value_t> const& M)
{
    return {"mvec2dc",
            {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11,
             M.c12, M.c13, M.c14, M.c15}};
}
[[maybe_unused]] JsonValue to_json(MVec2dc_E<value_t> const& M)
{
    return {"mvec2dc_e", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}
[[maybe_unused]] JsonValue to_json(MVec2dc_U<value_t> const& M)
{
    return {"mvec2dc_u", {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7}};
}
[[maybe_unused]] JsonValue to_json(MVec3dc<value_t> const& M)
{
    return {"mvec3dc",
            {M.c0,  M.c1,  M.c2,  M.c3,  M.c4,  M.c5,  M.c6,  M.c7,  M.c8,  M.c9,  M.c10,
             M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17, M.c18, M.c19, M.c20, M.c21,
             M.c22, M.c23, M.c24, M.c25, M.c26, M.c27, M.c28, M.c29, M.c30, M.c31}};
}
[[maybe_unused]] JsonValue to_json(MVec3dc_E<value_t> const& M)
{
    return {"mvec3dc_e",
            {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11,
             M.c12, M.c13, M.c14, M.c15}};
}
[[maybe_unused]] JsonValue to_json(MVec3dc_U<value_t> const& M)
{
    return {"mvec3dc_u",
            {M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11,
             M.c12, M.c13, M.c14, M.c15}};
}

[[maybe_unused]] JsonValue to_json(DualNum2dp<value_t> const& d)
{
    return {"dualnum2dp", {d.c0, d.c1}};
}
[[maybe_unused]] JsonValue to_json(DualNum3dp<value_t> const& d)
{
    return {"dualnum3dp", {d.c0, d.c1}};
}

// --------------------------------------------------------------------------- //
// Case collection
// --------------------------------------------------------------------------- //

struct Case {
    std::string id;
    std::string submodule; // "ega" | "pga" | "cga" | "sta" | "top"
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
// CGA cases (cga2dc = G(3,1,0), cga3dc = G(4,1,0) --- types *_2dc / *_3dc)
// --------------------------------------------------------------------------- //
// Conformal points are null vectors, so the round objects are joins of points
// and the flats are the same joins with the point at infinity appended. The
// cases below cover both directions (build by join, read back by cen / car /
// con / par / att / radius_sq) plus the conformal versors.

void emit_cga_cases()
{
    // --- cga3dc: null embedding and the point-pair / circle / sphere joins ---
    {
        vec3dc a = round_point3dc(1.0, 0.0, 0.0, 0.0);
        vec3dc b = round_point3dc(0.0, 1.0, 0.0, 0.0);
        vec3dc c = round_point3dc(0.0, 0.0, 1.0, 0.0);
        vec3dc d = round_point3dc(-1.0, 0.0, 0.0, 0.0);
        add("cga_vec3dc_dot", "cga", "dot", dot(a, b), a, b);
        add("cga_vec3dc_add", "cga", "+", a + b, a, b);
        add("cga_vec3dc_sub", "cga", "-", a - b, a, b);
        add("cga_vec3dc_neg", "cga", "neg", -a, a);
        add("cga_vec3dc_wdg_dipole", "cga", "wdg", wdg(a, b), a, b);
        add("cga_vec3dc_wdg_circle", "cga", "wdg", wdg(wdg(a, b), c), wdg(a, b), c);
        add("cga_vec3dc_wdg_sphere", "cga", "wdg", wdg(wdg(wdg(a, b), c), d),
            wdg(wdg(a, b), c), d);
        add("cga_vec3dc_gpr", "cga", "*", a * b, a, b);
        add("cga_vec3dc_rev", "cga", "rev", rev(wdg(a, b)), wdg(a, b));
        add("cga_vec3dc_conj", "cga", "conj", conj(wdg(a, b)), wdg(a, b));
        add("cga_vec3dc_cconj", "cga", "cconj", cconj(wdg(a, b)), wdg(a, b));
        add("cga_vec3dc_gr_inv", "cga", "gr_inv", gr_inv(wdg(a, b)), wdg(a, b));
    }

    // --- cga3dc: reading the geometry back off an object ---
    {
        quadvec3dc s = sphere3dc(3.0, 4.0, 0.0, 2.0);
        add("cga_sphere3dc_radius_sq", "cga", "radius_sq", radius_sq(s), s);
        add("cga_sphere3dc_center_nrm_sq", "cga", "center_nrm_sq", center_nrm_sq(s), s);
        add("cga_sphere3dc_cen", "cga", "cen", cen(s), s);
        add("cga_sphere3dc_par", "cga", "par", par(s), s);
        add("cga_sphere3dc_unitize", "cga", "unitize", unitize(s), s);
        add("cga_sphere3dc_dot_self", "cga", "dot", dot(unitize(s), unitize(s)),
            unitize(s), unitize(s));

        trivec3dc circ = circle3dc(1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 1.0);
        add("cga_circle3dc_car", "cga", "car", car(circ), circ);
        add("cga_circle3dc_con", "cga", "con", con(circ), circ);
        add("cga_circle3dc_ccr", "cga", "ccr", ccr(circ), circ);
        add("cga_circle3dc_att", "cga", "att", att(circ), circ);

        vec3dc rp = round_point3dc(1.0, 2.0, 3.0, 2.0);
        add("cga_round_point3dc_con", "cga", "con", con(rp), rp);
        add("cga_round_point3dc_radius_sq", "cga", "radius_sq", radius_sq(rp), rp);
        add("cga_round_point3dc_round_bulk", "cga", "round_bulk", round_bulk(rp), rp);
        add("cga_round_point3dc_round_weight", "cga", "round_weight", round_weight(rp),
            rp);

        bivec3dc fp = flat_point3dc(1.0, 2.0, 3.0);
        add("cga_flat_point3dc_flat_bulk", "cga", "flat_bulk", flat_bulk(fp), fp);
        add("cga_flat_point3dc_flat_weight", "cga", "flat_weight", flat_weight(fp), fp);
    }

    // --- cga3dc: meet --- a line cutting a sphere off-center gives a point pair ---
    {
        quadvec3dc s = sphere3dc(0.0, 0.0, 0.0, 2.0);
        trivec3dc l = line3dc(0.0, 1.0, 0.0, 1.0, 0.0, 0.0);
        bivec3dc dip = rwdg(s, l);
        add("cga_sphere_line_rwdg", "cga", "rwdg", dip, s, l);
        add("cga_dipole3dc_radius_sq", "cga", "radius_sq", radius_sq(dip), dip);
        add("cga_dipole3dc_cen", "cga", "cen", cen(unitize(dip)), unitize(dip));
        add("cga_dipole3dc_att", "cga", "att", att(unitize(dip)), unitize(dip));

        quadvec3dc pl = plane3dc(0.0, 0.0, 1.0, 1.0);
        add("cga_sphere_plane_rwdg", "cga", "rwdg", rwdg(s, pl), s, pl);
    }

    // --- cga3dc: complements and duals (odd-dimensional --- single variant) ---
    {
        vec3dc a = round_point3dc(1.0, 2.0, 3.0, 1.0);
        add("cga_vec3dc_cmpl", "cga", "cmpl", cmpl(a), a);
        add("cga_vec3dc_dual", "cga", "dual", dual(a), a);
        add("cga_vec3dc_antidual", "cga", "antidual", antidual(a), a);
        add("cga_vec3dc_inv", "cga", "inv", inv(a), a);
    }

    // --- cga3dc: the conformal versors (odd-grade motors in 3dc) ---
    {
        quadvec3dc s = sphere3dc(0.0, 0.0, 0.0, 2.0);
        mvec3dc_u T = get_translation(1.0, 2.0, 3.0);
        mvec3dc_u R = get_rotation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.5);
        mvec3dc_u D = get_dilation(0.0, 0.0, 0.0, 3.0);
        add("cga_get_translation", "cga", "get_translation", T, 1.0, 2.0, 3.0);
        add("cga_get_rotation", "cga", "get_rotation", R, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            0.5);
        add("cga_get_dilation", "cga", "get_dilation", D, 0.0, 0.0, 0.0, 3.0);
        add("cga_transform_sphere_translate", "cga", "transform", transform(s, T), s, T);
        add("cga_transform_sphere_dilate", "cga", "transform", transform(s, D), s, D);
        add("cga_rgpr_versor_compose", "cga", "rgpr", rgpr(T, R), T, R);
        add("cga_rlog_rotation", "cga", "rlog", rlog(R), R);
        add("cga_rexp_rlog_rotation", "cga", "rexp", rexp(rlog(R)), rlog(R));
        add("cga_rsqrt_rotation", "cga", "rsqrt", rsqrt(R), R);
        add("cga_rrev_versor", "cga", "rrev", rrev(R), R);

        // inversion in a sphere: the point at 2 maps to r^2 / 2
        vec3dc p = round_point3dc(2.0, 0.0, 0.0, 0.0);
        quadvec3dc us = sphere3dc(0.0, 0.0, 0.0, 1.0);
        add("cga_invert_on_sphere", "cga", "invert_on", invert_on(p, us), p, us);
    }

    // --- cga2dc: the planar algebra (circles are trivectors, dipoles bivectors) ---
    {
        vec2dc a = round_point2dc(1.0, 0.0, 0.0);
        vec2dc b = round_point2dc(0.0, 1.0, 0.0);
        vec2dc c = round_point2dc(-1.0, 0.0, 0.0);
        add("cga_vec2dc_dot", "cga", "dot", dot(a, b), a, b);
        add("cga_vec2dc_wdg_dipole", "cga", "wdg", wdg(a, b), a, b);
        add("cga_vec2dc_wdg_circle", "cga", "wdg", wdg(wdg(a, b), c), wdg(a, b), c);
        add("cga_vec2dc_gpr", "cga", "*", a * b, a, b);

        trivec2dc c1 = circle2dc(0.0, 0.0, 2.0);
        trivec2dc c2 = circle2dc(3.0, 0.0, 2.0);
        add("cga_circle2dc_rwdg", "cga", "rwdg", rwdg(c1, c2), c1, c2);
        add("cga_circle2dc_radius_sq", "cga", "radius_sq", radius_sq(c1), c1);
        add("cga_circle2dc_cen", "cga", "cen", cen(c2), c2);
        add("cga_circle2dc_att", "cga", "att", att(c1), c1);
        add("cga_circle2dc_unitize", "cga", "unitize", unitize(c2), c2);

        // dipole2dc(center_x, center_y, radius, direction_x, direction_y)
        // (the same dipole as before the argument became a direction: the
        // former normal (1, 0) is the direction (0, -1))
        bivec2dc dip = dipole2dc(0.0, 0.0, 1.0, 0.0, -1.0);
        add("cga_dipole2dc_radius_sq", "cga", "radius_sq", radius_sq(dip), dip);
        add("cga_dipole2dc_car", "cga", "car", car(dip), dip);

        // versors are EVEN-grade in cga2dc (odd in cga3dc)
        mvec2dc_e T2 = get_translation(1.0, 2.0);
        mvec2dc_e R2 = get_rotation(0.0, 0.0, 0.6);
        add("cga_get_translation2dc", "cga", "get_translation", T2, 1.0, 2.0);
        add("cga_get_rotation2dc", "cga", "get_rotation", R2, 0.0, 0.0, 0.6);
        add("cga_transform_circle2dc", "cga", "transform", transform(c1, T2), c1, T2);
        add("cga_rgpr_versor_compose2dc", "cga", "rgpr", rgpr(T2, R2), T2, R2);
        add("cga_rlog_rotation2dc", "cga", "rlog", rlog(R2), R2);
    }
}

// --------------------------------------------------------------------------- //
// STA cases (G(1,3,0) --- Space-Time Algebra, types *_4ds)
// --------------------------------------------------------------------------- //
// Metric is (-,-,-,+): for a vector nrm_sq = -x^2 - y^2 - z^2 + w^2, so a vector
// is spacelike when nrm_sq < 0 and timelike when nrm_sq > 0 (w the time component).
// STA is even-dimensional, so duals / complements split into l_/r_ variants
// (no plain cmpl / dual), exactly as for pga3dp.

void emit_sta_cases()
{
    // --- vec4ds basic algebra ---
    {
        vec4ds a{1.0, 2.0, 3.0, 4.0}, b{5.0, -6.0, 7.0, -8.0};
        add("sta_vec4ds_dot", "sta", "dot", dot(a, b), a, b);
        add("sta_vec4ds_wdg", "sta", "wdg", wdg(a, b), a, b);
        add("sta_vec4ds_add", "sta", "+", a + b, a, b);
        add("sta_vec4ds_sub", "sta", "-", a - b, a, b);
        add("sta_vec4ds_neg", "sta", "neg", -a, a);
        add("sta_vec4ds_nrm_sq", "sta", "nrm_sq", nrm_sq(a), a);
        add("sta_vec4ds_rev", "sta", "rev", rev(a), a);
        // geometric product of two vectors -> mvec4ds_e (scalar + bivector parts)
        add("sta_vec4ds_gpr", "sta", "*", a * b, a, b);
        // complement / dual of a vector -> trivec4ds (left and right variants differ)
        add("sta_vec4ds_l_cmpl", "sta", "l_cmpl", l_cmpl(a), a);
        add("sta_vec4ds_r_cmpl", "sta", "r_cmpl", r_cmpl(a), a);
        add("sta_vec4ds_l_dual", "sta", "l_dual", l_dual(a), a);
        add("sta_vec4ds_r_dual", "sta", "r_dual", r_dual(a), a);
        // inverse of a non-null vector -> vec4ds (nrm_sq != 0; here spacelike)
        vec4ds c{2.0, 6.0, -4.0, 1.0}; // nrm_sq = -55, invertible
        add("sta_vec4ds_inv", "sta", "inv", inv(c), c);
    }

    // --- bivec4ds basic algebra ---
    {
        bivec4ds B1{1.0, -2.0, 3.0, -4.0, 5.0, -6.0};
        bivec4ds B2{0.5, 1.5, -0.5, 2.0, -1.0, 0.5};
        add("sta_bivec4ds_add", "sta", "+", B1 + B2, B1, B2);
        add("sta_bivec4ds_sub", "sta", "-", B1 - B2, B1, B2);
        add("sta_bivec4ds_neg", "sta", "neg", -B1, B1);
        add("sta_bivec4ds_nrm_sq", "sta", "nrm_sq", nrm_sq(B1), B1);
        add("sta_bivec4ds_rev", "sta", "rev", rev(B1), B1);
        add("sta_bivec4ds_rrev", "sta", "rrev", rrev(B1), B1);
        add("sta_bivec4ds_conj", "sta", "conj", conj(B1), B1);
        // bivec * bivec -> mvec4ds_e (scalar + bivector + pseudoscalar)
        add("sta_bivec4ds_gpr", "sta", "*", B1 * B2, B1, B2);
        // complement of a bivector -> bivec4ds
        add("sta_bivec4ds_l_cmpl", "sta", "l_cmpl", l_cmpl(B1), B1);
        add("sta_bivec4ds_r_cmpl", "sta", "r_cmpl", r_cmpl(B1), B1);
    }

    // --- trivec4ds basic algebra ---
    {
        trivec4ds T{2.0, -1.0, 3.0, 4.0};
        add("sta_trivec4ds_neg", "sta", "neg", -T, T);
        add("sta_trivec4ds_nrm_sq", "sta", "nrm_sq", nrm_sq(T), T);
        add("sta_trivec4ds_rev", "sta", "rev", rev(T), T);
        // complement of a trivector -> vec4ds
        add("sta_trivec4ds_l_cmpl", "sta", "l_cmpl", l_cmpl(T), T);
        add("sta_trivec4ds_r_cmpl", "sta", "r_cmpl", r_cmpl(T), T);
    }

    // --- spacetime split of a vector relative to the standard observer g4 ---
    {
        vec4ds x{2.0, 3.0, 5.0, 7.0};
        vec4ds u{0.0, 0.0, 0.0, 1.0}; // g4, unit timelike observer
        add("sta_time_split", "sta", "time_split", time_split(x, u), x, u);
        add("sta_space_split", "sta", "space_split", space_split(x, u), x, u);
    }

    // --- rotation (spatial) and boost (Lorentz) via the rotor sandwich ---
    {
        // rotate g1 by 90 deg in the g12 plane -> -g2  (cos90 g1 - sin90 g2)
        bivec4ds g12{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
        auto Rrot = get_rotor(g12, deg2rad(90.0));
        vec4ds g1{1.0, 0.0, 0.0, 0.0};
        add_case("sta_transform_rotate_g1_90deg", "sta", "transform", transform(g1, Rrot),
                 1e-9, g1, Rrot);

        // boost the rest 4-velocity g4 by rapidity phi in the g14 plane
        bivec4ds g14{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        auto Rboost = get_boost(g14, 0.5);
        vec4ds g4{0.0, 0.0, 0.0, 1.0};
        add_case("sta_transform_boost_g4", "sta", "transform", transform(g4, Rboost),
                 1e-9, g4, Rboost);

        // sqrt(rotor) halves the angle/rapidity: sqrt(R(x)) == R(x/2).
        // Use a generic rotor (boost composed with rotation) exercising all 8 coeffs.
        auto Rgen = get_boost(g14, 0.5) * get_rotor(g12, 0.6);
        add_case("sta_sqrt_rotor", "sta", "sqrt", sqrt(Rgen), 1e-9, Rgen);

        // closed-form transform_opt (vec / bivec / trivec) -- must match transform()
        vec4ds v{2.0, 3.0, 5.0, 7.0};
        bivec4ds B{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        trivec4ds t{1.0, 2.0, 3.0, 4.0};
        add_case("sta_transform_opt_vec", "sta", "transform_opt", transform_opt(v, Rgen),
                 1e-9, v, Rgen);
        add_case("sta_transform_opt_bivec", "sta", "transform_opt",
                 transform_opt(B, Rgen), 1e-9, B, Rgen);
        add_case("sta_transform_opt_trivec", "sta", "transform_opt",
                 transform_opt(t, Rgen), 1e-9, t, Rgen);
    }

    // --- projection / rejection / reflection (vector onto vector) ---
    {
        vec4ds v{2.0, 3.0, 5.0, 7.0};
        vec4ds t{0.0, 0.0, 0.0, 1.0}; // non-null target (timelike)
        add("sta_project_onto", "sta", "project_onto", project_onto(v, t), v, t);
        add("sta_reject_from", "sta", "reject_from", reject_from(v, t), v, t);
        add("sta_reflect_on_vec", "sta", "reflect_on_vec", reflect_on_vec(v, t), v, t);
    }

    // --- angle (spacelike) / rapidity (timelike) --- both return plain double ---
    {
        // two purely spatial (spacelike) vectors -> ordinary Euclidean angle
        vec4ds a{1.0, 2.0, 3.0, 0.0}, b{2.0, 0.0, 1.0, 0.0};
        add("sta_angle_spacelike", "sta", "angle", angle(a, b), a, b);
        // two timelike 4-velocities -> relative rapidity
        vec4ds u1{0.0, 0.0, 0.0, 1.0}, u2{0.5, 0.0, 0.0, 1.0};
        add("sta_rapidity_timelike", "sta", "rapidity", rapidity(u1, u2), u1, u2);
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
    add("top_sign_neg", "top", "sign", static_cast<double>(sign(-3.5)), -3.5);
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
    out << R"(  "note": "AUTO-GENERATED by ga_test/python_utilities/export_python_cases.cpp -- do not edit by hand. Regenerate via: cmake --build . --target ga_export_python_cases && ./ga_test/ga_export_python_cases",)"
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
    // Default output path is baked in at build time (absolute, so the
    // exporter works regardless of the working directory it is launched
    // from). A command-line argument still overrides it.
#ifdef _GA_PY_CASES_JSON_DEFAULT
    std::string path = _GA_PY_CASES_JSON_DEFAULT;
#else
    std::string path = "ga_py/tests/data/ga_test_cases.json";
#endif
    if (argc > 1) path = argv[1];

    emit_ega_cases();
    emit_pga_cases();
    emit_cga_cases();
    emit_sta_cases();
    emit_top_cases();

    write_json(path);
    fmt::print("Wrote {} cases to {}\n", g_cases.size(), path);
    return 0;
}
