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

It compares three kinds:

| kind        | ga_py source              | ga_lua source (regex on `ga_lua.hpp`) |
| ----------- | ------------------------- | ------------------------------------- |
| `types`     | classes on the submodule  | `new_usertype<NAME>`                  |
| `functions` | callables on the submodule| `set_function("NAME", ...)`           |
| `constants` | data attributes           | `lua["NAME"] = ...`                   |

### Usage

Build `ga_py` first, then run with the `ga_py` venv (so `import ga_py` works):

```bash
cmake -S . -B build -D_GA_BUILD_PYTHON=ON
cmake --build build --target _ga_py

ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py            # full report
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --summary  # totals only
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --algebra=sta
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --kind=types
ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py --present  # also list bound
```

The script adds `build/ga_py` and `ga_py/python` to `sys.path` itself, so it works
from any cwd as long as the in-tree `build/` is current.

### What it does and does NOT check

- **Name-level, not overload-level.** `ga_lua` binds one global `nrm` carrying
  per-type `sol::overload`s; `ga_py` exposes `ega.nrm`, `pga.nrm`, … A name reported
  "bound" is present *somewhere* in `ga_lua` — the tool does not verify that every
  per-type overload exists. The dominant gaps (a whole missing algebra, or an unbound
  function) still surface clearly through the type/constant lists. Audit per-overload
  completeness by reading the relevant `register_*` block.
- **Operators are out of scope.** `ga_py` exposes them as dunder methods, `ga_lua` as
  `sol::meta_function` entries inside each `new_usertype`. They ride along with the
  type, so they are not listed as free functions.
- **`ga_py` must be current.** The tool reads the built module, so rebuild `_ga_py`
  after any binding/regeneration change before trusting the numbers.

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
  (`transform`/`transform_opt`, the space-time splits, `is_{timelike,spacelike,lightlike}`,
  `rapidity`, `get_boost`, `l_undual`/`r_undual`), and STA overloads on the 26 shared
  `set_function`s (`nrm`/`dot`/`inv`/`exp`/`log`/`sqrt`/duals/complements/`wdg`/`rwdg`/…)
  via `splice_overloads.py`. `register_functions` gained `using namespace hd::ga::sta`.
- **Phase 4 — top-level helpers — DONE** (modulo deliberate skips): bound the unit
  conversions (`Hz2radps`, `radps2Hz`, `radps2rpm`, `rpm2radps`), `rgr`, and `sign`
  (`gr` was bound in Phase 2). The solver / integrator free functions (`lu_solve`,
  `lstsq_solve`, `rk4_step`, `rk4_get_time`) operate on matrix/vector spans and don't map
  to a Lua REPL — **skipped by decision**.

## Status

`ga_lua` is at `ga_py` parity except for deliberate skips — **395/409** bound
(`lua_coverage.py --summary`). EGA and STA are 100%; PGA is complete bar the
`inertia2dp/3dp` matrix type and the 8 dynamics functions; the top-level remainder is
the 4 solver/integrator helpers. Re-run `lua_coverage.py` after any library API change
to catch new gaps.
