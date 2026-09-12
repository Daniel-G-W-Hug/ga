// Hand-written bindings for the generic two-argument templates: the projection /
// antiprojection families, the euclidean distances, try_unitize, and the defaulted
// tolerance of is_simple.
//
// Why hand-written: these are declared as fully generic forwarding templates,
//
//     template <typename arg1, typename arg2> decltype(auto) ortho_proj3dp(arg1&&,
//     arg2&&)
//
// so the manifest records them with `param_types: ["arg1 &&", "arg2 &&"]` and a
// `decltype(auto)` return, which the emitter's resolver has nothing to map onto a
// nanobind signature. The scanner DOES see them (FUNCTION_TEMPLATE is collected beside
// FUNCTION_DECL) -- the obstacle is the generic form, not the template. Contrast
// project_onto, which is bound automatically: it is templated on the scalar type only and
// names concrete types in its signature.
//
// Why not reimplement them in Python (which is what this file replaces): the projections
// carry two things a Python composition cannot reproduce without copying constants -- the
// ideal-target threshold in detail::by_weight_sq (divide only when
// weight_nrm_sq > safe_epsilon^2) and the blade-target check under
// _HD_GA_EXTENDED_TEST_BLADE_TARGET -- and dist*dp carries an `if constexpr` dispatch on
// the grade sum, which a Python shim can only mirror through a hand-kept table of class
// names. All three had drifted. Binding the library function makes drift impossible.
//
// What is NOT here, deliberately: the bulk/weight contractions and expansions are generic
// templates too, but each is a one-line composition of already-bound primitives
// (l_bulk_contract2dp(a, b) = rwdg(bulk_dual(a), b)), so the Python forwarders in
// ga_py/python/ga_py/__init__.py ARE that composition and cover every grade pair the
// bound rwdg / wdg / dual accept. Enumerating pairs here would narrow them.
//
// The operand pairs are hand-kept: one per meaningful grade combination, with
// gr(a) < gr(b) for a projection and gr(a) > gr(b) for an antiprojection (the library's
// requires, which does not compile otherwise). A pair omitted here is a missing Python
// overload -- a TypeError -- rather than a wrong answer. The derived geometric primitives
// (point*, line*, plane*) need no entries of their own: their bindings declare the graded
// type as a base class, so nanobind upcasts them into these overloads.

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>

#include <utility>

#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"
#include "ga/ga_sta.hpp"

namespace nb = nanobind;
using namespace hd::ga;

namespace {

// One concrete instantiation of a generic two-argument template. `op` is a generic
// lambda forwarding to the library function, so a family states its name, its docstring
// and its operand pairs once each.
template <typename A, typename B, typename Op>
void bind_pair(nb::module_& m, char const* name, Op op, char const* doc)
{
    m.def(
        name, [op](A const& a, B const& b) { return op(a, b); }, nb::arg("a"),
        nb::arg("b"), doc);
}

// try_unitize reports whether it divided through a trailing bool* out-parameter, which
// has no nanobind mapping -- that is why the generator skipped it entirely. In Python it
// returns the pair instead: (object, unitized).
template <typename T> void bind_try_unitize(nb::module_& m, char const* doc)
{
    m.def(
        "try_unitize",
        [](T const& x) {
            bool unitized = false;
            auto const r = pga::try_unitize(x, &unitized);
            return std::pair{r, unitized};
        },
        nb::arg("x"), doc);
}

char const* const DOC_ORTHO_PROJ =
    "Orthogonal projection of a onto the larger-grade b. The target's scale is "
    "divided out and the source's is kept, so f(a, k*b) == f(a, b) and "
    "f(k*a, b) == k*f(a, b). The result is NOT unitized -- call try_unitize() on "
    "it for a canonical representative, since a projection may legitimately land "
    "at infinity. REQUIRES gr(a) < gr(b), and a target that is a blade.";

char const* const DOC_CENTRAL_PROJ =
    "Central projection of a onto the larger-grade b (towards the origin): the "
    "same expression with the bulk expansion in place of the weight one. Same "
    "scaling contract as ortho_proj. REQUIRES gr(a) < gr(b).";

char const* const DOC_ORTHO_ANTIPROJ =
    "Orthogonal antiprojection: the smallest object of a's grade that contains b. "
    "REQUIRES gr(a) > gr(b).";

char const* const DOC_CENTRAL_ANTIPROJ =
    "Central antiprojection: the antiprojection counterpart of central_proj, "
    "completing the four-way family. REQUIRES gr(a) > gr(b).";

char const* const DOC_DIST2DP =
    "Euclidean distance between two 2dp objects as a homogeneous magnitude, "
    "returned as dualnum2dp(distance_numerator, weight). Divide the first by the "
    "second for the metric distance.";

char const* const DOC_DIST3DP =
    "Euclidean distance between two 3dp objects as a homogeneous magnitude, "
    "returned as dualnum3dp(distance_numerator, weight). Divide the first by the "
    "second for the metric distance.";

char const* const DOC_TRY_UNITIZE =
    "Unitize if the object carries weight, else return it unchanged. Returns "
    "(object, unitized): the flag is False for an ideal object, which unitize() "
    "would refuse instead.";

void bind_pga2dp(nb::module_& m)
{
    auto const op = [](auto const& a, auto const& b) { return pga::ortho_proj2dp(a, b); };
    bind_pair<vec2dp, bivec2dp>(m, "ortho_proj2dp", op, DOC_ORTHO_PROJ);

    auto const cp = [](auto const& a, auto const& b) {
        return pga::central_proj2dp(a, b);
    };
    bind_pair<vec2dp, bivec2dp>(m, "central_proj2dp", cp, DOC_CENTRAL_PROJ);

    auto const oa = [](auto const& a, auto const& b) {
        return pga::ortho_antiproj2dp(a, b);
    };
    bind_pair<bivec2dp, vec2dp>(m, "ortho_antiproj2dp", oa, DOC_ORTHO_ANTIPROJ);

    auto const ca = [](auto const& a, auto const& b) {
        return pga::central_antiproj2dp(a, b);
    };
    bind_pair<bivec2dp, vec2dp>(m, "central_antiproj2dp", ca, DOC_CENTRAL_ANTIPROJ);

    auto const d = [](auto const& a, auto const& b) { return pga::dist2dp(a, b); };
    // line-line is absent because the library's expression is: the grade sum is 4,
    // past the pseudoscalar, so neither branch of dist2dp is well-formed for it
    bind_pair<vec2dp, vec2dp>(m, "dist2dp", d, DOC_DIST2DP);
    bind_pair<vec2dp, bivec2dp>(m, "dist2dp", d, DOC_DIST2DP);
    bind_pair<bivec2dp, vec2dp>(m, "dist2dp", d, DOC_DIST2DP);
}

void bind_pga3dp(nb::module_& m)
{
    auto const op = [](auto const& a, auto const& b) { return pga::ortho_proj3dp(a, b); };
    bind_pair<vec3dp, bivec3dp>(m, "ortho_proj3dp", op, DOC_ORTHO_PROJ);
    bind_pair<vec3dp, trivec3dp>(m, "ortho_proj3dp", op, DOC_ORTHO_PROJ);
    bind_pair<bivec3dp, trivec3dp>(m, "ortho_proj3dp", op, DOC_ORTHO_PROJ);

    auto const cp = [](auto const& a, auto const& b) {
        return pga::central_proj3dp(a, b);
    };
    bind_pair<vec3dp, bivec3dp>(m, "central_proj3dp", cp, DOC_CENTRAL_PROJ);
    bind_pair<vec3dp, trivec3dp>(m, "central_proj3dp", cp, DOC_CENTRAL_PROJ);
    bind_pair<bivec3dp, trivec3dp>(m, "central_proj3dp", cp, DOC_CENTRAL_PROJ);

    auto const oa = [](auto const& a, auto const& b) {
        return pga::ortho_antiproj3dp(a, b);
    };
    bind_pair<bivec3dp, vec3dp>(m, "ortho_antiproj3dp", oa, DOC_ORTHO_ANTIPROJ);
    bind_pair<trivec3dp, vec3dp>(m, "ortho_antiproj3dp", oa, DOC_ORTHO_ANTIPROJ);
    bind_pair<trivec3dp, bivec3dp>(m, "ortho_antiproj3dp", oa, DOC_ORTHO_ANTIPROJ);

    auto const ca = [](auto const& a, auto const& b) {
        return pga::central_antiproj3dp(a, b);
    };
    bind_pair<bivec3dp, vec3dp>(m, "central_antiproj3dp", ca, DOC_CENTRAL_ANTIPROJ);
    bind_pair<trivec3dp, vec3dp>(m, "central_antiproj3dp", ca, DOC_CENTRAL_ANTIPROJ);
    bind_pair<trivec3dp, bivec3dp>(m, "central_antiproj3dp", ca, DOC_CENTRAL_ANTIPROJ);

    auto const d = [](auto const& a, auto const& b) { return pga::dist3dp(a, b); };
    // line-plane and plane-plane are absent for the same reason as line-line in 2dp:
    // the grade sum exceeds 4, so the wedge the else-branch needs does not exist
    bind_pair<vec3dp, vec3dp>(m, "dist3dp", d, DOC_DIST3DP);
    bind_pair<vec3dp, bivec3dp>(m, "dist3dp", d, DOC_DIST3DP);
    bind_pair<vec3dp, trivec3dp>(m, "dist3dp", d, DOC_DIST3DP);
    bind_pair<bivec3dp, vec3dp>(m, "dist3dp", d, DOC_DIST3DP);
    bind_pair<bivec3dp, bivec3dp>(m, "dist3dp", d, DOC_DIST3DP);
    bind_pair<trivec3dp, vec3dp>(m, "dist3dp", d, DOC_DIST3DP);
}

} // anonymous namespace

void bind_projections(nb::module_& m_ega, nb::module_& m_pga, nb::module_& m_sta)
{
    bind_pga2dp(m_pga);
    bind_pga3dp(m_pga);

    // EGA 3d: the only higher-grade target below the pseudoscalar is the bivector, so
    // there is one pair. (ega2d has none at all -- projecting onto its pseudoscalar is
    // the identity, which is why no ortho_proj2d exists.)
    auto const ep = [](auto const& a, auto const& b) { return ega::ortho_proj3d(a, b); };
    bind_pair<vec3d, bivec3d>(m_ega, "ortho_proj3d", ep, DOC_ORTHO_PROJ);

    auto const sp = [](auto const& a, auto const& b) { return sta::ortho_proj4ds(a, b); };
    bind_pair<vec4ds, bivec4ds>(m_sta, "ortho_proj4ds", sp, DOC_ORTHO_PROJ);
    bind_pair<vec4ds, trivec4ds>(m_sta, "ortho_proj4ds", sp, DOC_ORTHO_PROJ);
    bind_pair<bivec4ds, trivec4ds>(m_sta, "ortho_proj4ds", sp, DOC_ORTHO_PROJ);

    // The derived primitives come FIRST: nanobind resolves overloads in registration
    // order, so binding vec3dp before point3dp would upcast a point and return a vec.
    bind_try_unitize<point2dp>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<line2d>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<point3dp>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<line3d>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<plane3d>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<vec2dp>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<bivec2dp>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<vec3dp>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<bivec3dp>(m_pga, DOC_TRY_UNITIZE);
    bind_try_unitize<trivec3dp>(m_pga, DOC_TRY_UNITIZE);

    // is_simple's tolerance defaults to eps_congruent in C++; the generator drops
    // default arguments, so the one-argument form is added here. Registered AFTER the
    // generated two-argument overload, which keeps taking an explicit tolerance.
    m_pga.def(
        "is_simple", [](bivec3dp const& B) { return pga::is_simple(B); }, nb::arg("B"),
        "Is this bivector a blade, i.e. does it represent a line? The Pluecker "
        "condition wdg(B, B) == 0, measured relative to the sum of the squared "
        "components. Tolerance defaults to eps_congruent.");
    m_sta.def(
        "is_simple", [](bivec4ds const& B) { return sta::is_simple(B); }, nb::arg("B"),
        "Is this bivector a blade, i.e. does it represent a 2-plane? The Pluecker "
        "condition wdg(B, B) == 0, measured relative to the sum of the squared "
        "components. Tolerance defaults to eps_congruent.");
}
