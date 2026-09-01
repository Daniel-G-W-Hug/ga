# ga_lua coverage utilities

Helpers for keeping the **Lua** wrapper (`ga_lua`) in step with the **Python**
wrapper (`ga_py`). The Python binding is the parity target: it is the
manifest-driven (`ga_bindgen/manifest.json`), POD-filtered realisation of the public
GA API, so "ga_lua is as complete as ga_py" is exactly "`lua_coverage.py` reports
nothing missing".

## `lua_coverage.py`

Introspects the **built** `ga_py` extension and diffs its exposed surface against
what `ga_lua/src/ga_lua.hpp` binds, reporting the gaps **per algebra**
(`ega` / `pga` / `sta`) plus the top-level free functions.

It compares four kinds:

| kind        | ga_py source                 | ga_lua source (parsed from `ga_lua.hpp`) |
| ----------- | ---------------------------- | ---------------------------------------- |
| `types`     | classes on the submodule     | `new_usertype<NAME>`                     |
| `functions` | callables on the submodule   | `set_function("NAME", ...)`              |
| `constants` | data attributes              | `lua["NAME"] = ...`                      |
| `operators` | ordered pairs where `a + b` / `a - b` succeeds | `sol::resolve<R(A, B)>(operator+)`, placed where Lua can actually reach it |

### The `operators` kind

Types-only checking says nothing about the operator surface, and that blind spot
was worth 340 of 428 missing type pairs until 2026-08-30. This kind closes it, and
it models Lua's dispatch rather than merely grepping for the resolve:

- for `a + b` Lua takes `__add` from **a**'s metatable and only falls back to
  **b**'s if a has none at all;
- so a resolve for `A + B` sitting in B's block, when A binds `operator+` itself,
  is **dead code** -- it compiles, reads correctly, and can never fire. Those are
  listed under `=== UNREACHABLE ===` (there were 10 such at the time this was
  written, e.g. `pscalar2dp + dualnum2dp`);
- `sol::base_classes` widens reachability: a resolve taking a base accepts a
  derived argument, so `point3dp + vec3dp` is reachable through `vec3dp`'s set.
  The converse bites too -- two derived operands with no operator between them
  fail with "attempt to perform arithmetic", which base classes do not fix.

The target set is obtained by TRYING every ordered pair in ga_py rather than from
a table of type names, so cross-family pairs (a 2d type with a 3d one, sharing the
`ega` submodule but no algebra) drop out by themselves and newly bound types are
picked up without editing the script.

### Usage

Build `ga_py` first, then run with the `ga_py` venv (so `import ga_py` works):

```bash
cmake -S . -B build -D_GA_BUILD_PYTHON=ON
cmake --build build --target _ga_py

ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py            # full report
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --summary  # totals only
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --algebra=sta
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --kind=types
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --kind=operators
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --present  # also list bound
```

The script adds `build/ga_py` and `ga_py/python` to `sys.path` itself, so it works
from any cwd as long as the in-tree `build/` is current.

### What it does and does NOT check

- **Per-table, but not overload-level.** Since free functions moved into the
  per-algebra tables (`ega.wdg`, `pga.rgpr`, `cga.dot`, `sta.wdg`), the tool matches
  the TABLE too: a function bound only in `ega` no longer counts as present for
  `pga`. That distinction found a real mis-filing the moment it was added
  (`expand` sat in `ega` while `ga_py` has it as `pga.expand`). It still does not
  verify that every per-type overload within a table exists. The dominant gaps (a whole missing algebra, or an unbound
  function) still surface clearly through the type/constant lists. Audit per-overload
  completeness by reading the relevant `register_*` block.
- **Operators are out of scope.** `ga_py` exposes them as dunder methods, `ga_lua` as
  `sol::meta_function` entries inside each `new_usertype`. They ride along with the
  type, so they are not listed as free functions.
- **`ga_py` must be current.** The tool reads the built module, so rebuild `_ga_py`
  after any binding/regeneration change before trusting the numbers.

### Two things it cannot see

- **A forwarder that is defined but broken.** The `register_forwarders` prelude is
  Lua source inside a C++ string: the tool (and the compiler) see the definitions,
  not whether their bodies resolve. All 26 were dead at once after the per-algebra
  move — their bodies still called the old global names — while every input script
  stayed green because none of them called a forwarder. `ga_lua/input/test_forwarders.lua`
  now calls all 26, so that path has a gate.
- **Placement that differs from `ga_py` on purpose.** `gr` / `rgr` are top-level in
  `ga_py` and per-algebra here (they dispatch on argument type either way), so a
  top-level target found in any table counts as bound.

## `gen_lua_overloads.py`

Drafts the `sol::resolve` overload blocks that close the gaps `lua_coverage.py` finds
(the *find* vs *apply* split mirrors ga_prdxpr's validation_utilities vs utilities).
From `ga_bindgen/manifest.json` it emits one `lua.set_function("name", sol::overload(
... ))` per requested name, in `ga_lua`'s exact `sol::resolve<Ret(Args)>(name)` style.

```bash
ga_py/.venv/bin/python ga_lua/utilities/gen_lua_overloads.py \
    project_onto,rotate,angle                 # default --algebra=ega
ga_py/.venv/bin/python ga_lua/utilities/gen_lua_overloads.py \
    sqrt,log,twdg1 --algebra=ega,pga          # merge families into one block
```

Mechanics: overloads are selected by the manifest `namespace` (`hd::ga::<family>`);
template types are concretised to Lua aliases by lowercasing the outer name
(`Vec2d`->`vec2d`, `MVec2d_E`->`mvec2d_e`), with `value_t`/`bool`/`std::vector<...>`
and bare scalar template params handled specially; reference vs value passing is
preserved. Duplicate signatures (a function listed once per header) are de-duplicated.

It is a **drafting aid** — review each block before pasting:

- a name spanning families must bind once, so pass every family it covers in a single
  `--algebra=` list (a second `set_function` for the same name silently replaces the
  first); the generated overloads for not-yet-registered type families (e.g. STA before
  Phase 3) must wait until those `new_usertype`s exist;
- decide per overload whether it belongs in the REPL;
- after pasting, `clang-format` the file, build `ga_lua`, and re-run `lua_coverage.py`.

## `splice_overloads.py`

The *apply* counterpart to `gen_lua_overloads.py`. That one drafts a whole new
`set_function` block for a name nothing binds yet; this one **extends an existing
block** with one family's overloads, in place. Use it when a name is already bound for
a family and newly added library overloads have to join the same block.

A name is bound **per algebra table** — `ega.wdg`, `pga.wdg`, `cga.wdg` and `sta.wdg`
are four separate bindings of one name, and 53 names sit on more than one table. So
`--algebra` selects both the manifest namespace to read *and* the table to edit.

```bash
python3 ga_lua/utilities/splice_overloads.py --algebra=cga --dry-run   # report only
python3 ga_lua/utilities/splice_overloads.py --algebra=cga             # all bound names
python3 ga_lua/utilities/splice_overloads.py --algebra=sta nrm,dot,inv # a subset
```

It reads the family's overloads from the manifest (sharing the concretiser with the
generator), paren-matches the enclosing `sol::overload(`, and inserts the missing
`sol::resolve` lines before its close under a `// <family>` comment. It creates no new
blocks and skips any name that is not already bound.

**`--exclude REGEX` expresses a per-CLASS decision.** The tool works per *name*, but a
name's overloads rarely all belong to one class -- 4 of pga `wdg`'s 70 gaps take the
named `point`/`line`/`plane` types and the other 66 do not. Any generated signature
matching the regex is skipped and counted, so a class decision becomes a recorded,
repeatable invocation instead of a hand-edit afterwards:

```bash
EXCL='\b(point|line|plane)[0-9]|std::vector|dualnum'     # named types, batch, dualnum
python3 splice_overloads.py --algebra=pga --exclude "$EXCL" "wdg,rwdg,rgpr,dot"
# cmt/rcmt additionally: a commutator with a SCALAR is zero in every algebra, and with
# a PSEUDOSCALAR wherever I is central (ega3d, pga2dp) -- binding those adds surface
# that can only answer 0
python3 splice_overloads.py --algebra=pga \
    --exclude "$EXCL"'|\bscalar|\bpscalar3d\b|\bpscalar2dp\b' "cmt,rcmt"
```

**It is idempotent, and that has to survive clang-format.** The presence test compares
signatures with all whitespace removed, because the file it edits has already been
formatted and a long signature is wrapped at whatever column the formatter chose:

```cpp
sol::resolve<std::vector<vec4ds>(std::vector<vec4ds> const&,
                                 mvec4ds_e const&)>(transform_opt),
```

A line-based substring test does not find that, and appends a **second, identical**
`sol::resolve` — which compiles, since `sol::overload` accepts duplicates, so nothing
downstream complains. Measured on the tree as it stood: 7 such duplicates across four
families, all in the longest signatures (the `std::vector` batch overloads of
`transform_opt`/`rotate_opt`, the 7-argument `circle3dc`/`dipole3dc`). Fixed 2026-09-01;
if the tool ever grows another comparison against the file's text, normalise it the
same way.

After splicing: `clang-format` the file, build `ga_lua`, re-run `lua_coverage.py`, run a
smoke script, and review the diff — decide per overload whether it belongs in the REPL.

## Closure roadmap (bring `ga_lua` to `ga_py` parity)

Two scope decisions are fixed: **bind only the pure-data PODs** (the stateful
`*_system` classes stay unbound, exactly as in `ga_py`), and **take orientation from
`ga_py`** (parity where reasonably possible). Work the phases in order; end each by
re-running `lua_coverage.py` and a smoke `.lua` script, and watch the MISSING count
fall.

- **Phase 1 — EGA2D / EGA3D — DONE** (100%: types 12/12, functions 37/37, constants
  43/43). Bound the 9 missing functions (`gs_orthogonal`/`gs_orthonormal`,
  `l_expand`/`r_expand`, `rotate_opt`, `log`, `sqrt`, `twdg1`, `rtwdg1`) and the 6
  `one_*` constants. Most of the larger geometric set (`project_onto`, `reflect_on`,
  `rotate`, `angle`, `is_congruent`, …) was already bound.
- **Phase 2 — PGA2DP / PGA3DP — DONE** (effectively complete: types 38/40, functions
  83/91, constants 80/80; the rest are intentional skips). Split into: **2a** the
  concrete C++ functions (`rinv`, `move*_opt`, `invert_on`) + all PGA constants; the
  Lua-prelude **forwarders** (`register_forwarders`: the 24 contraction / expansion /
  projection / `dist` templates, mirroring ga_py's Python one-liners over bound
  primitives) + `gr()`; **2b-i** the convenience **types** (`point/vector/line/plane`,
  via `register_convenience_types`, subclasses bound with `sol::base_classes`) +
  `expand`; **2b-ii** the physics **PODs** (`pose`, `kin_state`, `joint_state`,
  `loop_constraint`, via `register_physics_pods`) + the scoped enums
  (`joint*`/`constraint*`) + `motor_from_pose`/`pose_from_motor`.
  **Skipped by decision:** the `inertia2dp/3dp` matrix type and the dynamics functions
  (`make_*_body`, `get_*_inertia`, `compute_omega_dot`, `get_inertia_inverse`).
- **Phase 3 — STA4D — DONE** (100%: types 9/9, functions 39/39, constants 45/45). Added
  `#include "ga/ga_sta.hpp"` and `register_4ds_types` (adapted from `register_3dp_types`
  by `3dp→4ds` / `pga→sta`, since STA shares pga3dp's grade structure and underlying
  `_t` templates). Added the 45 STA constants, the 13 STA-only functions
  (`transform`/`transform_opt`, the space-time splits, `rapidity`, `get_boost`,
  `is_{timelike,spacelike,lightlike}`, `l_undual`/`r_undual`), and STA overloads on the
  26 shared `set_function`s
  (`nrm`/`dot`/`inv`/`exp`/`log`/`sqrt`/duals/complements/`wdg`/`rwdg`/…) via
  `splice_overloads.py`. `register_functions` gained `using namespace hd::ga::sta`.
  Smoke script: `ga_lua/input/test_sta.lua` (added 2026-08-16 alongside the CGA one).
- **Phase 5 — CGA2DC / CGA3DC — DONE** (100%: types 19/19, functions 70/70, constants
  137/137). Added `#include "ga/ga_cga.hpp"` and `register_2dc_types` /
  `register_3dc_types` (mirroring `register_4ds_types`; two metric differences from STA:
  the pseudoscalar's self-product is a *scalar*, and the CGA dual numbers are the plain
  `MVec2_t` without the scalar/pseudoscalar mixing overloads). The 43 CGA-only functions
  (`cen`/`car`/`con`/`par`/`ccr`, `radius_sq`, the round/flat bulk+weight norms, the
  object constructors, `get_translation`/`get_rotation`/`get_dilation`/
  `get_transversion`/`get_loxodromic`, `antidual`/`l_antidual`/`r_antidual`, `cconj`)
  came from `gen_lua_overloads.py --algebra=cga`; CGA overloads on the 27 shared
  `set_function`s from `splice_overloads.py --algebra=cga`; the 137 constants generated
  from `ga_usr_consts.hpp`'s `hd::ga::cga` block. `register_functions` and
  `register_constants` gained `using namespace hd::ga::cga`. Also bound in the same pass,
  since they were missing for *every* algebra: `is_close` (ega,pga,cga,sta),
  `is_same_rotation` (ega), `is_same_motion` (pga), `is_same_transform` (cga,sta), and
  the top-level `signum`. Smoke script: `ga_lua/input/test_cga.lua`.

- **Phase 4 — top-level helpers — DONE** (modulo deliberate skips): bound the unit
  conversions (`Hz2radps`, `radps2Hz`, `radps2rpm`, `rpm2radps`), `rgr`, and `sign`
  (`gr` was bound in Phase 2). The solver / integrator free functions (`lu_solve`,
  `lstsq_solve`, `rk4_step`, `rk4_get_time`) operate on matrix/vector spans and don't map
  to a Lua REPL — **skipped by decision**.

## Status

`ga_lua` is at `ga_py` parity except for deliberate skips — **630/653** bound
(`lua_coverage.py --summary`). EGA, CGA and STA are 100%; PGA is complete bar the
`inertia2dp/3dp` matrix type and the 8 dynamics functions; the top-level remainder is
the solver/integrator helpers plus `linear_step`/`smooth_step`/`smoother_step`,
`to_geo_pos` and the `geo_angle` type. Re-run `lua_coverage.py` after any library API change
to catch new gaps.
