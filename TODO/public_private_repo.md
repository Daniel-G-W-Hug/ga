# Public / Private repository separation

Split the codebase along one simple, durable rule: the **generic GA library and its
generic demonstrations stay public**; **specific application cases move to a private repo,
code and docs together**. The two initial application cases are **wafer grinding** and
**magnetic levitation**. A strong, mirrored structural link ties the two repos together.

**Why a code-and-docs-together rule (not case-by-case confidentiality):** deciding
file-by-file whether something is "confidential enough" is unworkable long-term. Instead:
an *application* (a concrete real-world use case — a machine, a process) lives in the
private repo in full; the *library* and its *generic demonstrations* (textbook mechanisms,
feature demos) stay public. This is easy to apply going forward and needs no per-file
judgement.

Two application cases differ only in **git-history treatment** (see §7):

- **Magnetic levitation** — contains confidential physics → move **and purge from public
  history** (earlier commits must become undownloadable from the remote).
- **Wafer grinding** — nothing confidential yet → **move only** (plain removal); no history
  purge needed. It moves because it is an application case, not because it is secret.

> **Phase 0 — DONE (2026-07-07): `ga` is now private.** Verified via the GitHub API (an
> unauthenticated request to `Daniel-G-W-Hug/ga` returns HTTP 404, vs. 200 for a public
> control repo). This stops all ongoing public exposure and removes the time pressure from
> the whole split. The consequence for the plan: the entire reorganization and the maglev
> **history purge run while the repo is private**, and the repo is flipped **back to
> public only at the very end**, once the tree and history are verified clean. That closes
> the SHA-fetch window §7 otherwise could not (see the reordered §6/§7).

---

## 1. Guiding principles

1. **Public = generic library + generic demos.** The header-only `ga/` library (incl. the
   general-purpose force elements and `dynamic_system`/`closed_loop_system`), `ga_prdxpr`,
   `ga_lua`, `ga_py`, the core algebra tests, the **generic mechanism demos** (four-bar,
   double pendulum, merry-go-round, planar delta, closed-loop 5-bar, open-vs-closed) and
   **feature demos** (projection, reflection, bivectors), and the generic `ga_docu`
   chapters.
2. **Private = application cases (code + docs).** Wafer grinding and magnetic levitation,
   each as a self-contained bundle: tests, benches, view scenes, and doc chapter.
3. **Dependency direction is one-way.** Private code depends on the public library; public
   code never depends on private code. Already true for maglev/dipole (verified: no refs
   in `ga/`, `ga_prdxpr/`, `ga_lua/`, `ga_py/`).
4. **Mirror the structure.** A private file sits at the same relative path as its public
   sibling (`ga_test/src/…`, `ga_view/src/…`, `ga_docu/…`), so moving is mechanical.
5. **Bundle each application behind a name prefix.** Everything for an application shares a
   prefix (`ga_maglev_*`, `ga_grinding_*`, `active_maglev*`, `active_grinding_*`) so a
   whole case is one grep and one move. This drives the extraction refactor in §4a.
6. **Confidential content leaves public history; non-confidential application code only
   leaves `HEAD`.** Two-tier, see §7.

---

## 2. What moves vs. what stays

### 2a. The boundary (decided 2026-07-07)

- Generic mechanism/feature demos **stay public** as library documentation — they are not
  "applications".
- Everything **grinding-labeled** (`grinding|wafer|tao|sommerfeld|zhou`) moves private as
  one bundle, including the generic building-block cases ("tumbling plane", "intersecting
  discs") and the Sommerfeld rotordynamics warm-up. One rule, no per-case exceptions.
- All maglev/dipole moves private.

### 2b. MAGLEV bundle → private (move + purge history, §7)

| Public path today | Kind |
| ----------------- | ---- |
| `ga_test/src/ga_dipole_model.hpp` | private physics core (dipole field/force/torque, Halbach group, `maglev_rig`, `mover_body`, `solve_alpha`) |
| `ga_test/src/ga_dipole_test.cpp` | test target `ga_dipole_test` |
| `ga_test/src/ga_maglev_test.cpp` | test target `ga_maglev_test` |
| `ga_test/utilities/bench_dipole_ga_vs_classical.cpp` | bench `ga_dipole_bench` |
| `ga_test/utilities/bench_maglev_stability.cpp` | bench `ga_maglev_stability` |
| `ga_test/utilities/maglev_scene_replica.cpp` | validator `maglev_scene_replica` |
| `ga_view/src/active_maglev.{hpp,cpp}` | view scene (seam, §4c/§5) |
| `TODO/magnetic_levitation.md`, `TODO/maglev_S1_dipole_formulation.md`, `TODO/CtrlXFlow6d.pptx` | docs (already git-excluded) |

### 2c. GRINDING bundle → private (move only, no history purge, §7)

Currently **interleaved** in shared 3dp test headers and must be **extracted into dedicated
files first** (§4a). The grinding-labeled cases (verified list):

- In [`ga_appl3dp_appl_test.hpp`](../ga_test/src/ga_appl3dp_appl_test.hpp): tumbling plane
  & intersecting discs (grinding building blocks), tumbling tool plane via
  `static_system3dp`, grain trajectory (Phase 0), wafer-tilt flatness (Zhou Figs 5–7),
  grinding-mark wavelengths (D.1), wafer topography (D.1b), contact force feed-forward
  (D.2a), Sommerfeld warm-up (B.1).
- In [`ga_appl3dp_mechanics_test.hpp`](../ga_test/src/ga_appl3dp_mechanics_test.hpp), whole
  `TEST_SUITE`s: Tao wheel-spindle (C.2–C.4), grinding force loop (D.2b), wafer thinning
  (D.2c), loop integrator (D.2d), feed control (D.2e).
- View scenes: `ga_view/src/active_grinding_{cs,flatness,marks,topo}.{hpp,cpp}` (4 scenes,
  same hand-sync surface as maglev — see §5).
- Python demos: `ga_py/demo/{tumbling_plane,intersecting_discs}.py` (grinding-labeled
  building blocks).
- Docs: `TODO/grinding.md`, `TODO/tao_eq13_derivation.md`,
  `TODO/Grinding_setup_CS_draft.md`; and the grinding half of the `ga_docu` applications
  chapter (moves whole, §4d).

> Note: the **force-control work (Phase D.2a–e) is already written and public** — it is not
> future work. It moves now with the rest of the grinding bundle.

### 2d. Stays PUBLIC

- The whole `ga/` library, incl. force elements (`set_joint_spring_damper`,
  `set_applied_wrench`, `add_grounded_spring`, `grounded_spring{2,3}dp`), `dynamic_system`,
  `closed_loop_system` — general-purpose GA mechanics.
- `ga_prdxpr`, `ga_lua`, `ga_py` (minus the two grinding demos), core algebra tests,
  `bench_sta4ds_transform`.
- The generic mechanism/feature view scenes and their tests; the 2D application tests
  (`ga_appl2dp_*` — verified grinding-free).
- The generic `ga_docu` chapters (intro, basics, modelling-motion, mechanics, comparison,
  glossary, literature).

---

## 3. Link model (DECIDED): submodule superset

The private repo is a **superset that embeds the public repo as a git submodule**, pinned
to an exact public commit. Location: `/Users/hud3bh/prg/cpp/pj/ga_private` (sibling of
`ga`, mirroring the `industrial_ai` layout).

```text
ga/                              PUBLIC repo (github.com/Daniel-G-W-Hug/ga)
ga_private/                      PRIVATE repo (superset)
├── external/ga/                 ← git submodule → public repo, pinned to a commit
├── ga_test/src/     ga_maglev_*, ga_dipole_*, ga_grinding_*  (application tests)
├── ga_test/utilities/           maglev benches + replica
├── ga_test/private_targets.cmake        declares the private test/bench targets
├── ga_view/src/     active_maglev.*, active_grinding_*.*  (+ registration TUs)
├── ga_view/private_scenes.cmake         adds the private scene sources to ga_view
├── ga_docu/         6_ga_applications_pga.tex  (applications chapter)
├── TODO/            confidential + application docs
└── CMakeLists.txt   → add_subdirectory(external/ga); pass GA_PRIVATE_DIR down
```

The library-consuming pieces (tests, benches, replica) move trivially — they only
`target_link_libraries(... ga)` (maglev additionally includes the private
`ga_dipole_model.hpp`). The structural work is the `ga_view` registration seam (§4c/§5)
and the grinding extraction refactor (§4a).

### 3a. Workflows after the split

**Public user** — `git clone …/ga`. Builds the library, `ga_prdxpr`, `ga_lua`, `ga_py`,
public `ga_view` (generic demos only), and all public tests. `GA_PRIVATE_DIR` unset → no
application target is declared. The tree has zero maglev content in `HEAD` **and** history,
and no grinding at `HEAD` (grinding remains in older public history — not confidential).

**Private user** — `git clone --recurse-submodules …/ga_private`. CMake configures the
public submodule via `add_subdirectory(external/ga)` and passes
`-DGA_PRIVATE_DIR=${CMAKE_SOURCE_DIR}` down, lighting up the application test/bench targets
and folding the maglev + grinding scenes into `ga_view`. One `ga_view` shows generic demos
**plus** the application scenes.

**Bumping the public dependency** — advance the pin explicitly:

```bash
cd external/ga && git fetch && git checkout <new-public-commit> && cd -
git add external/ga && git commit -m "bump public ga to <short-sha>"
```

records exactly which public commit each private state builds against (the "strong link").

**Documentation** — generic GA docs stay in public `ga_docu`. The applications chapter and
confidential write-ups live in `ga_private/`. A private doc references a public one by
pinned path (`external/ga/ga_docu/…`).

**Implementing private add-ons on public parts** — one working copy holds both: edit public
code inside `external/ga` (a normal public checkout on its own branch), private code in the
parent tree. Commits land in whichever repo owns the file; the public change is later
pinned in via a submodule bump.

---

## 4. Prerequisite public-repo refactors (do FIRST, in public)

Behaviour-preserving changes committed to public `develop` so the application bundles can
leave without breaking any public build.

### 4a. Extract the interleaved grinding cases into dedicated files (the new big item)

Grinding is currently woven through `ga_appl3dp_appl_test.hpp` and
`ga_appl3dp_mechanics_test.hpp`. To move it as whole files:

1. Create a grinding test executable mirroring the existing pattern: `ga_grinding_test`
   built from `ga_test/src/ga_grinding_test.cpp` including
   `ga_grinding_appl_test.hpp` (kinematics/geometry/topology cases) and
   `ga_grinding_mechanics_test.hpp` (Tao spindle + D.2 force-loop suites).
2. **Cut** every grinding-labeled `TEST_CASE`/`TEST_SUITE` (the §2c list) out of the two
   `ga_appl3dp_*` headers into the new grinding headers. Leave the generic PGA application
   tests behind (audit the remainder to confirm none is grinding-labeled).
3. Verify `ga_appl3dp_test` still passes (now generic-only) and the new `ga_grinding_test`
   passes with the same assertions — this is a pure relocation, no assertion changes.
4. Commit this **in the public repo** first (it is a legitimate public refactor). Only then
   move `ga_grinding_*` to private.

This establishes the **per-application file convention** (principle 5): each application is
its own `ga_<app>_test` executable + `active_<app>_*` scenes + doc chapter.

### 4b. `ga_test` — optional private target overlay

```cmake
# Opt-in private overlay: extra test/bench targets living outside this repo.
if(GA_PRIVATE_DIR AND EXISTS "${GA_PRIVATE_DIR}/ga_test/private_targets.cmake")
    include("${GA_PRIVATE_DIR}/ga_test/private_targets.cmake")
endif()
```

`private_targets.cmake` (in `ga_private`) declares `ga_maglev_test`, `ga_dipole_test`,
`ga_grinding_test`, and the maglev benches/replica (linking `doctest::doctest ga`, include
dir = private root). Unset `GA_PRIVATE_DIR` → nothing added.

### 4c. `ga_view` — generic scene registry seam (one binary, add-on compiled in)

Five scenes now leave public (maglev + 4 grinding), each with the **same six-point
hand-sync** (see §5). That tips the design decisively to a **generic active-item registry**
rather than per-scene special-casing:

- Replace the per-scene `std::vector<a…>` members in `coordsys_model` and the per-scene
  `add_<name>` / `clear()` lines with a **type-erased active-item registry**, so
  `coordsys_model` carries no application-specific type at all. This removes the whole
  hand-sync class of problems for every current and future scene — public scene files
  become application-free.
- Add a registration hook the app iterates in `populate_scene`/`build_models`, e.g.
  `register_scene(scene_descriptor)` where the descriptor carries the model-builder and the
  active-item factory (the `populate_scene` block).
- CMake seam in `ga_view/CMakeLists.txt`:

  ```cmake
  if(GA_PRIVATE_DIR AND EXISTS "${GA_PRIVATE_DIR}/ga_view/private_scenes.cmake")
      include("${GA_PRIVATE_DIR}/ga_view/private_scenes.cmake")  # adds sources
      target_compile_definitions(ga_view PRIVATE GA_VIEW_PRIVATE_SCENES)
  endif()
  ```

  `private_scenes.cmake` appends `active_maglev.cpp`, `active_grinding_*.cpp`, and their
  `register_*_scene.cpp` TUs (each self-registers) to the `ga_view` target.

> Fallback if the full registry is too big for one pass: move the application types into
> private headers and gate their registration behind a single `#ifdef
> GA_VIEW_PRIVATE_SCENES` call site — still one binary, still no application code in public,
> at the cost of one guarded stub in `scenes.cpp`. With **five** scenes moving, the generic
> registry is strongly preferred; the `#ifdef` fallback scales worse.

### 4d. `ga_docu` — move the whole applications chapter private

The applications chapter [`6_ga_applications_pga.tex`](../ga_docu/6_ga_applications_pga.tex)
(Wafer Grinding + Magnetic Levitation) is one self-contained file, `\input` once from
[`0_ga_docu.tex`](../ga_docu/0_ga_docu.tex) (l.186). Move it whole:

1. Move `6_ga_applications_pga.tex` to `ga_private/ga_docu/`.
2. Guard the include in public `0_ga_docu.tex`:

   ```latex
   % Application chapter lives in the private repo; included only in the private build.
   \IfFileExists{6_ga_applications_pga.tex}{\input{6_ga_applications_pga.tex}}{}
   ```

   Public: file absent → chapter omitted, builds clean. Private: TEXINPUTS/symlink resolves
   it (public `ga_docu/6_ga_applications_pga.tex` is git-ignored so the private file drops
   in), chapter appears.
3. **Fix the one dangling cross-ref:** [`8_ga_glossary.tex`](../ga_docu/8_ga_glossary.tex)
   l.210 `\ref{subsubsection:wafer_grinding}` → reword to drop the `\ref{}` (verified: the
   only public cross-ref into the chapter).

Future private chapters attach the same way — one guarded `\IfFileExists{…}{\input}{}` per
chapter in `0_ga_docu.tex`, no public `.tex` edits. History: whole-file → clean path purge
(§7).

### 4e. Confirm the clean cut

```bash
grep -rln -iE 'maglev|dipole|grind|wafer|tao|sommerfeld|zhou' \
    ga/ ga_prdxpr/ ga_lua/ ga_py/ga_py/ ga_docu/
# after 4a–4d + moves: expect no matches (grinding demos & app chapter gone)
```

---

## 5. `ga_view` hand-sync surface (what 4c must absorb)

Each moving scene touches the same six shared points. Maglev, for reference:

- `scenes.cpp`: `#include` (l.26); item-creation block in `populate_scene` (l.1688–1701,
  connects reset/pause/trace/forces/magnets slots); registration block (l.2854–2865,
  `add_maglev` + label + legend).
- `coordsys_model.hpp`: `maglev_params`, `amaglev`, `add_maglev` decl, `std::vector<amaglev>
  amag`.
- `coordsys_model.cpp`: `add_maglev` impl, `amag.clear()` in `clear()`.
- `ga_view/CMakeLists.txt`: the `.cpp` source.

The 4 grinding scenes (`active_grinding_{cs,flatness,marks,topo}`) each replicate this:
`grinding_*_params` + `agrinding_*` structs, `add_grinding_*` methods, `*_mark`/`agr*`
vectors + their `clear()` lines, plus scenes.cpp blocks. The generic registry (4c) collapses
all five into self-registering descriptors and removes every one of these hand-sync points
from public files.

---

## 6. Migration steps (order of operations)

**The whole sequence runs while `ga` is private (Phase 0 done); the repo is flipped back to
public only in the final step.**

0. **`ga` is private** (done, 2026-07-07). No public exposure during any step below.
1. **(Optional but recommended) Bundle tests by application first** — the intermediate step
   in [test_bundling_by_application.md](test_bundling_by_application.md): extract the
   scattered grinding cases into `ga_grinding_test`, fold dipole into a single
   `ga_maglev_test`, de-mix Cai. Verify assertion-count parity. This makes step 4 a
   whole-file move.
2. **Public refactors (§4a–§4d), committed to `develop`.** Verify the build + full test
   suite pass with `GA_PRIVATE_DIR` unset (grinding still builds from `ga_grinding_test`,
   maglev from its files at this point).
3. **Create `ga_private`** with the submodule superset (§3); add `external/ga`.
4. **Move** the maglev bundle (§2b) and grinding bundle (§2c) into mirrored private paths.
   Add `private_targets.cmake` (§4b), `private_scenes.cmake` + registration TUs (§4c), and
   the applications chapter (§4d). Get the private superset building; all application
   tests, benches, and view scenes pass **from the private repo**.
5. **Remove from `ga`** (plain `git rm`): the maglev files, the grinding files
   (`ga_grinding_*`, `active_grinding_*`, the two Python demos), and the application doc
   chapter; delete their references from public CMake / `coordsys_model` (now behind the
   registry). Verify build + tests green; `ga_view` shows only generic demos.
6. **History: purge MAGLEV (§7) — while private.** Grinding is move-only, no purge.
7. **Move the git-excluded TODO docs** into `ga_private/TODO/`; drop the `.git/info/exclude`
   entries.
8. **Update docs**: `CLAUDE.md` (the split + per-app file convention + where each
   application lives), `README.md`, `TODO/github_setup.md` (branch model for both repos).
9. **Verify the clean state (§8)** on `ga`: grep gate, history gate, public doc build — all
   green — then **flip `ga` back to public**. Because the purge happened while private,
   there is never a public window in which the removed commits are fetchable by SHA.

---

## 7. Git history — two-tier, executed while private

Because `ga` is private (Phase 0), the purge is done **before the repo is made public
again** — so the rewritten-away commits are never publicly served, not even transiently by
SHA. This is strictly cleaner than purging an already-public repo, where unreachable objects
stay fetchable until GitHub GC.

**Maglev (confidential) — purge from history so the eventual public repo never carried it.**
Scope is narrow (verified): the maglev/dipole files live **only on `develop`**; `main`,
`sta4ds`, and the coworker `develop_becker` carry **zero**. Procedure:

1. **Backup** a full mirror (`git clone --mirror`) before any rewrite.
2. **`git filter-repo`** on `develop` to purge these paths from every reachable commit:
   `ga_test/src/ga_dipole_model.hpp`, `ga_test/src/ga_dipole_test.cpp`,
   `ga_test/src/ga_maglev_test.cpp`, the two maglev benches, `maglev_scene_replica.cpp`,
   `ga_view/src/active_maglev.{hpp,cpp}`, and `ga_docu/6_ga_applications_pga.tex` (the app
   chapter contains the confidential maglev section; purging it also drops the
   non-confidential grinding doc from history — harmless).
   Do the §4 refactors + the plain `git rm` of these files as a *normal* commit first, so
   only historical revisions need scrubbing.
3. **Force-push** `develop`: `git push --force-with-lease origin develop`.
4. **Clear residual reachability on GitHub** (still private, no public exposure risk):
   close/delete any PR built on the affected history and delete stale branches/tags pointing
   into it. Optionally ask GitHub Support to `git gc`; while private this is precautionary,
   not exposure-critical.
5. **Coordinate the coworker re-clone** (`develop_becker`, per `TODO/github_setup.md`): after
   the force-push they must re-clone or hard-reset; a plain `git pull` conflicts on the
   rewritten history.
6. **Verify, then go public** (§6 step 9): the grep + history gates (§8) must be green on the
   private `ga` before flipping visibility back to public.

**Grinding (not confidential) — move only, NO history purge.** Remove the grinding files
from the `develop` tip in a normal commit. Older commits retain grinding; that is fine —
nothing there is secret. (Should any grinding content ever become confidential, a targeted
purge can be run then.)

> The git-excluded `TODO/` docs were never committed — confirm with
> `git log --all --oneline -- TODO/magnetic_levitation.md` (expect empty).
>
> Residual limitation (only relevant if `ga` was ever public before Phase 0): copies already
> cloned/forked/cached externally cannot be retracted — going private stops *future* access
> but is not retroactive. Check the repo's fork list once. For access gained *after* Phase 0
> there is no such window.

---

## 8. Verification / acceptance criteria

- **Public-only build** (`GA_PRIVATE_DIR` unset): configures/builds with no application
  targets; `ga_ega_test`, `ga_pga_test`, `ga_appl{2,3}dp_test` (now generic-only),
  `ga_sta_test`, etc. pass; `ga_view` launches and shows only the generic demos.
- **Grep gate**: the §4e grep returns nothing at `HEAD`.
- **History gate**: `git log --all --oneline -- ga_view/src/active_maglev.cpp` **and**
  `-- ga_docu/6_ga_applications_pga.tex` are empty after §7. (Grinding paths may still
  appear in history — expected, move-only.)
- **Public doc build**: `latexmk -pdf 0_ga_docu.tex` in a public checkout compiles with no
  undefined references (chapter omitted, glossary reworded).
- **Private build** (submodule superset): builds the public library + view app + all
  application targets; `ga_maglev_test`, `ga_dipole_test`, `ga_grinding_test`,
  `maglev_scene_replica` pass; maglev + grinding scenes appear in `ga_view` and match
  pre-split behaviour (compare against `maglev_scene_replica`); private doc build includes
  the applications chapter with references resolved.
- **No behaviour change from the grinding extraction (§4a)**: the relocated grinding
  assertions run identically before (public, pre-move) and after (private).
- **Structural link**: a private file and its public sibling share a relative path; the
  submodule resolves.

---

## 9. Decisions (resolved) and remaining open points

**Resolved (2026-07-07):**

- **Sequencing** — `ga` taken **private now** (Phase 0 done); the reorganization and the
  maglev history purge run while private, and `ga` is flipped **back to public only at the
  end**, after the clean-state gates pass. No public-exposure clock during the work.
- **Split rule** — public = generic library + generic demos; private = application cases
  (code + docs together). Generic mechanism/feature demos stay public.
- **Grinding scope** — everything grinding-labeled (incl. tumbling-plane / intersecting-
  discs building blocks + Sommerfeld) moves private as one bundle.
- **Link model** — submodule superset; public embedded at `external/ga`, pinned;
  `ga_private` sibling of `ga`.
- **`ga_view`** — one binary, add-on compiled in via a generic active-item registry
  (5 scenes moving makes this the clear choice over per-scene `#ifdef`).
- **Docs** — whole applications chapter moves private; guarded `\input`.
- **History** — two-tier: maglev purged from `develop`; grinding move-only.

**Still open, to settle when implementing:**

1. **Registry vs `#ifdef` for `ga_view`** — confirm the full generic active-item registry
   is in scope for the first pass (preferred), or start with the `#ifdef` fallback.
2. **Python demos** — `ga_py/demo/{tumbling_plane,intersecting_discs}.py` are grinding-
   labeled but pure generic geometry; move with grinding (per the one-rule decision) or keep
   as public geometry demos? Default: move.
3. **Private `ga_py`/`ga_lua` bindings** — application work is C++-only today; expose
   application types to Python/Lua in the private repo, or keep C++-only?
