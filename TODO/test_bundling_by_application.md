# Bundling tests by application case (intermediate step)

Reorganize the test suite so **each application case is one self-contained test
executable**, mirroring how `ga_pga_test.cpp` bundles all PGA tests. This is the
preparatory step for the public/private repo split
([public_private_repo.md](public_private_repo.md) §4a): once an application's tests live in
dedicated files, moving that application to the private repo is a whole-file move.

Deliverables: `ga_grinding_test` (new) and `ga_maglev_test` (consolidate the existing
maglev + dipole files) — plus a **naming convention** every future application follows.

---

## 1. Target pattern (mirror `ga_pga_test.cpp` / `ga_appl3dp_test.cpp`)

One executable per application, split into a thin entry `.cpp` + chapter `.hpp` files:

- **Entry `.cpp`** carries the doctest main and shared includes/usings, then `#include`s
  the chapter headers. Verbatim shape of `ga_appl3dp_test.cpp`:

  ```cpp
  #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
  #include "doctest/doctest.h"
  #include "ga_<app>_<chapter1>_test.hpp"
  #include "ga_<app>_<chapter2>_test.hpp"
  ```

- **Chapter `.hpp`** files hold the `TEST_SUITE`/`TEST_CASE`s (each keeps its own
  `#include`s + `using namespace` block — idempotent, matches the current
  `ga_appl3dp_*_test.hpp` convention).

**Ordering caveat:** doctest runs suites ordered by **filename** (`--order-by=file`
default; see the comment in `ga_appl3dp_test.cpp`). Chapter headers must sort in the
intended run order. Keeping the `_appl_` (`a`) chapter before `_mechanics_` (`m`) preserves
the current relative order. The tests are self-contained (no cross-suite state), so
splitting into a separate executable does not affect pass/fail.

---

## 2. Current state — three application clusters are mixed into shared files

Investigation of the two shared 3dp headers (case titles are the stable anchors; line
numbers drift as edits land):

### `ga_appl3dp_appl_test.hpp` (2742 lines, one `TEST_SUITE("PGA3DP: application tests")`)

| Cluster | Cases (title anchors) |
| ------- | --------------------- |
| **GRINDING** (contiguous, first in the suite) | reference & tumbling plane; intersecting discs; tumbling tool plane via `static_system3dp`; disc surface speeds via `kinematic_system3dp`; grain trajectory (Phase 0); wafer-tilt flatness Zhou Figs 5-7; grinding-mark wavelengths (D.1); wafer topography (D.1b); grinding contact force feed-forward (D.2a); Sommerfeld warm-up (B.1) |
| **generic PGA** | line perpendicular to non-intersecting lines (rcmt); bivector rate of change |
| **CAI machine (Phase F)** | Cai volumetric-error machine Eq.1 (F.0); Cai geometric error as a motor (F.1); Cai volumetric error map Eq.2/3 (F.2); Cai Abbe/Bryan Eq.5-8 (F.3); Cai wafer topography + TTV Eq.9-18 (F.4) |

The grinding cases are **contiguous at the top** (before the first generic case), so they
cut cleanly. No file-scope helpers exist in this file (each case defines its own locals,
e.g. the `disc` struct) — verified, so nothing straddles the boundary.

### `ga_appl3dp_mechanics_test.hpp` (2674 lines)

| Cluster | Suites |
| ------- | ------ |
| **generic mechanics** (top, lines ~1-1601) | physics tests prep (Inertia/rk4/disc-body); physics tests application (ODE cuboid pendulum); `dynamic_system3dp` M1 / M2 / M3; `closed_loop_system3dp`. **Includes the generic force-element tests** (spring/damper, grounded-spring stiffness, driven-joint) — these test *library features* and have **no grinding keyword** in their titles, so they stay generic. |
| **GRINDING** (contiguous tail, ~1603-EOF) | a Phase-C comment block; an **anonymous-namespace helper set** (`tao_spindle_params`, `make_rotor_body`, `spindle_idx`, `build_tao_spindle`, `measure_freq`); `TEST_SUITE`s: Tao wheel-spindle (C.2-C.4 + Dyn.1/Dyn.2); grinding force loop (D.2b); wafer thinning (D.2c); loop integrator (D.2d); feed control (D.2e) |

The grinding helpers sit **after** every generic suite (defined at ~1644, generic suites
end at ~1601), so no generic test can reference them — the grinding tail (comment +
helpers + suites) is one contiguous, self-contained block. Clean cut.

### Maglev — already isolated (separate executables)

- `ga_maglev_test.cpp` (692 lines): `TEST_SUITE`s "maglev: rotating Halbach actuator" +
  "maglev: hover digital twin".
- `ga_dipole_test.cpp` (195 lines): `TEST_SUITE("magnetostatic dipoles in GA")`.
- `ga_dipole_model.hpp`: shared physics core (`using namespace dipole_model;`).

No mixing with generic tests — but they are **two** executables, not one bundle.

---

## 3. Grinding extraction → `ga_grinding_test`

Create three new files under `ga_test/src/`:

1. **`ga_grinding_test.cpp`** — entry. Copy `ga_appl3dp_test.cpp` verbatim; replace the
   three `#include`s with:

   ```cpp
   #include "ga_grinding_appl_test.hpp"
   #include "ga_grinding_mechanics_test.hpp"
   ```

2. **`ga_grinding_appl_test.hpp`** — move the 10 GRINDING cases out of
   `ga_appl3dp_appl_test.hpp`. Give the head the same includes/usings block as the source
   file; wrap the cases in `TEST_SUITE("PGA3DP: grinding application tests")`.

3. **`ga_grinding_mechanics_test.hpp`** — move the GRINDING tail out of
   `ga_appl3dp_mechanics_test.hpp`: the Phase-C comment, the **anonymous-namespace helper
   set** (must travel with the suites), and the five grinding `TEST_SUITE`s. Same
   includes/usings head.

Then **delete** the moved content from the two `ga_appl3dp_*` files:

- `ga_appl3dp_appl_test.hpp` → retains the 2 generic cases (rcmt, bivector rate) + the Cai
  cluster (until §5 moves Cai out too).
- `ga_appl3dp_mechanics_test.hpp` → retains the generic mechanics suites (1-1601), which is
  the bulk of the file. Unchanged behaviour.

**CMake** (`ga_test/CMakeLists.txt`): add a `ga_grinding_test` executable mirroring the
existing `ga_appl3dp_test` target block (new `EXEC_NAME`, `set(GRINDING_FILES
src/ga_grinding_test.cpp)`, `add_executable`, `target_include_directories(... PRIVATE
${CMAKE_SOURCE_DIR})`, `target_link_libraries(... doctest::doctest ga)`).

---

## 4. Maglev consolidation → one `ga_maglev_test` bundle

The user's goal is *one file collecting everything for maglev*. Fold the dipole tests into
the maglev executable as a chapter:

1. **`ga_maglev_test.cpp`** becomes the thin entry: doctest main + the shared includes
   (`ga/ga_ega.hpp`/`ga_pga.hpp` as needed, `ga/detail/ga_solver.hpp`,
   `ga_test/src/ga_dipole_model.hpp`), then:

   ```cpp
   #include "ga_maglev_dipole_test.hpp"   // magnetostatics foundation
   #include "ga_maglev_twin_test.hpp"     // rotating actuator + hover twin
   ```

2. **`ga_maglev_dipole_test.hpp`** ← the body of today's `ga_dipole_test.cpp`
   (`TEST_SUITE("magnetostatic dipoles in GA")`).
3. **`ga_maglev_twin_test.hpp`** ← the two current maglev `TEST_SUITE`s.
4. **Retire** the standalone `ga_dipole_test` executable (its cases now run under
   `ga_maglev_test`); keep `ga_dipole_model.hpp` as the shared physics header. Update
   `ga_test/CMakeLists.txt` (drop `EXEC_NAME8`/`ga_dipole_test`; `ga_maglev_test` now builds
   from `ga_maglev_test.cpp` including the two chapters).

> Alternative (lighter): leave `ga_maglev_test.cpp` + `ga_dipole_test.cpp` as-is — they are
> already isolated and private-ready. Choose the fold only if a single maglev executable is
> wanted. Recommended: fold, so every application follows the same one-executable pattern.

---

## 5. Cai machine (Phase F) — a third application, same treatment

The Cai volumetric-error machine (5 Phase-F cases) is mixed into
`ga_appl3dp_appl_test.hpp` exactly like grinding was. While the file is open, extract it
too — otherwise the "generic" appl file is left holding a full second application:

- New `ga_cai_test.cpp` (entry) + `ga_cai_test.hpp` (`TEST_SUITE("PGA3DP: Cai
  volumetric-error machine")`), new CMake target `ga_cai_test`.

> **Open question (does not block bundling):** is Cai a public generic demo or a private
> application? Bundling it into `ga_cai_test` is worthwhile either way (it de-mixes the
> file); the public/private assignment is a separate call for the main split plan. Flagged
> in [public_private_repo.md](public_private_repo.md) §9.

After §3 + §5, `ga_appl3dp_appl_test.hpp` holds only the 2 generic cases (rcmt, bivector
rate) — thin enough to consider merging into `ga_appl3dp_3d_cmp_test.hpp` or renaming to a
generic-pga3dp file. Minor; decide during the edit.

---

## 6. Verification (pure relocation — behaviour must be identical)

- **Assertion-count parity.** Before the change, run `./ga_appl3dp_test`, `./ga_dipole_test`,
  `./ga_maglev_test` and record each doctest `assertions:`/`test cases:` total. After, the
  sum over `./ga_appl3dp_test` (generic-only) + `./ga_grinding_test` + `./ga_cai_test` +
  `./ga_maglev_test` (folded) must equal the original sum, all green. No assertion is added
  or removed — cases only move.
- **Public build unaffected:** `cmake --build .` succeeds; `ga_appl3dp_test` still passes
  (now generic-only). `ga_appl3dp_test` is the long pole (~3 min Debug); give it ≥300 s.
- **doctest line-collision guard:** the new entry `.cpp`s inherit
  `DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN`; if a `DOCTEST_CLANG >= 22` line-number collision
  surfaces (as guarded in `ga_pga_test.cpp` with `#define DOCTEST_COUNTER __COUNTER__`),
  add the same guard to the affected entry.
- **clang-format** every new/edited file (`clang-format -i`, in-tree) before done.

---

## 7. Convention for future application cases (the durable outcome)

Each application case = one bundle, name-prefixed `ga_<app>_`:

- `ga_<app>_test.cpp` — entry (doctest main + shared includes + chapter `#include`s).
- `ga_<app>_<chapter>_test.hpp` — one per natural grouping (`appl` for
  kinematics/geometry, `mechanics` for dynamics, …), named to sort into the intended
  doctest run order.
- matching `active_<app>_*` view scenes and a `ga_docu` chapter (per the main split plan).

A new application never lands in a shared `ga_appl3dp_*` file again — it gets its own
bundle from the start. This is what makes the eventual private-repo move a whole-file
operation.

---

## 8. Suggested order of work

1. Grinding first (largest, clearest): create the 3 files (§3), cut from the two shared
   headers, add the CMake target, verify assertion parity.
2. Cai next (§5): same de-mixing, verify the generic appl file is now trivial.
3. Maglev fold (§4): consolidate to one executable, retire `ga_dipole_test`.
4. Full verification (§6); commit as a self-contained public refactor **before** any
   repo-split work begins.
