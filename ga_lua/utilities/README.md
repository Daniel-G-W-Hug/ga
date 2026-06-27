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

## Closure roadmap (bring `ga_lua` to `ga_py` parity)

Two scope decisions are fixed: **bind only the pure-data PODs** (the stateful
`*_system` classes stay unbound, exactly as in `ga_py`), and **take orientation from
`ga_py`** (parity where reasonably possible). Work the phases in order; end each by
re-running `lua_coverage.py` and a smoke `.lua` script, and watch the MISSING count
fall.

- **Phase 1 — EGA2D / EGA3D** (types already bound; close function + constant gaps):
  projections/reflections/rotations (`project_onto`, `reject_from`, `reflect_on`,
  `rotate`/`rotate_opt`), `angle`, `gs_orthogonal`/`gs_orthonormal`, `det`,
  `is_congruent`, generic `gr`, `expand`/`l_expand`/`r_expand`, `l_undual`/`r_undual`,
  `abs`, `sign`; the missing `one_*` / metric constants.
- **Phase 2 — PGA2DP / PGA3DP**: regressive products `rgpr`/`rwdg`/`rdot`/`rinv`; the
  projective contractions and expansions; motor sandwiches `move2dp`/`move3dp`(`_opt`);
  `ortho_proj`/`ortho_antiproj`/`central_proj`, `dist2dp`/`dist3dp`, `is_congruent`.
  Add the geometric convenience **types** (`line2d/3d`, `point2d/2dp/3d/3dp`,
  `plane3d`, `vector2d/3d`), the PGA POD types (`pose`, `kin_state`, `joint_state`,
  `loop_constraint`, `inertia`), the simple scoped enums (`joint2dp`, `constraint2dp`),
  and the missing PGA multivector-form / metric constants.
- **Phase 3 — STA4D** (new algebra, the largest single gap — currently 0 types, 0
  constants bound): add `#include "ga/ga_sta.hpp"`, a `register_sta_types` block
  (`scalar4ds … mvec4ds`, `dualnum4ds`) mirroring an existing `register_*_types`;
  STA functions (`transform`/`transform_opt`, the space-time splits, causal predicates
  `is_{timelike,spacelike,lightlike}`, `rapidity`, `get_boost`/`get_rotor`,
  `exp`/`log`/`sqrt`, `l_expand4ds`/`r_expand4ds`); the STA basis constants.
- **Phase 4 — top-level helpers**: unit conversions (`Hz2radps`, `rpm2radps`, …) and
  any remaining free functions (`gr`, `rgr`, `sign`). The solver / integrator
  free functions (`lu_solve`, `lstsq_solve`, `rk4_step`, …) operate on
  matrix/vector spans and are a judgment call for a Lua REPL — bind last, or skip.
