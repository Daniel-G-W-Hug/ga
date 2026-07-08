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
>
> **FINAL PARTITION — MODULE-LEVEL (decided 2026-07-08).** The split is now cut mostly by
> whole module, which retires the delicate `ga_view` scene-by-scene seam:
>
> | Public | Private |
> | ------ | ------- |
> | `ga/` (library — ships the generated product code), `ga_bindgen/`, `ga_py/`, `ga_lua/`, **`ga_prdxpr/`** (kept public: regeneration + "missing products easily added" are documented features) | **`ga_view/` (whole module)** — viz explains GA but is not self-standing; moving it whole dissolves the scene seam and takes the maglev + grinding scenes with it |
> | `ga_test/` **minus** the maglev + grinding bundles | `ga_test/` maglev + grinding bundles + `maglev_utilities/`; `ga/detail`-free `ga_dipole_model.hpp` |
> | `ga_docu/` **generated PDF only** (§4d) | `ga_docu/` LaTeX sources |
> | top-level CMake, README, LICENSE, `.clang-format` | `TODO/` (whole, principle 7) |
>
> Within-module splits remaining: `ga_test` (generic vs. maglev/grinding — done via the test
> bundling) and `ga_docu` (LaTeX private / PDF public). Everything else is whole-directory.
> **History purge stays maglev-only** — `ga_prdxpr`/`ga_view` history is uncritical and is
> **not** rewritten (§7).

---

## 1. Guiding principles

1. **Public = generic library + generic demos.** The header-only `ga/` library (incl. the
   general-purpose force elements and `dynamic_system`/`closed_loop_system`), `ga_prdxpr`,
   `ga_lua`, `ga_py`, the core algebra tests, the **generic mechanism demos** (four-bar,
   double pendulum, merry-go-round, planar delta, closed-loop 5-bar, open-vs-closed) and
   **feature demos** (projection, reflection, bivectors). Documentation is a **generated
   PDF only** — the `ga_docu` LaTeX sources are private; public ships just the generic-
   chapters PDF (§4d).
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
7. **The entire `TODO/` directory is private (decided 2026-07-08).** All plans and working
   notes — tracked or git-excluded, application-specific or generic — move to the private
   repo. Plans and old plans do not belong in the public part regardless of the character
   of the task they describe. The public repo carries no `TODO/`. This removes any
   per-file plan triage (see the simplified doc handling in §2b/§2c). One consequence:
   public `CLAUDE.md` must not point at `TODO/*` files — those references get genericised
   or dropped in the public-repo docs pass (§8/§6 step 8).

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

Post-bundling layout (the test-bundling intermediate step is DONE — see
[test_bundling_by_application.md](test_bundling_by_application.md)):

| Public path today | Kind |
| ----------------- | ---- |
| `ga_test/src/ga_dipole_model.hpp` | private physics core (dipole field/force/torque, Halbach group, `maglev_rig`, `mover_body`, `solve_alpha`) |
| `ga_test/src/ga_maglev_test.cpp` | entry of the `ga_maglev_test` bundle (doctest main) |
| `ga_test/src/ga_maglev_dipole_test.hpp` | dipole chapter (was `ga_dipole_test.cpp`, folded in) |
| `ga_test/src/ga_maglev_twin_test.hpp` | actuator + hover-twin chapter (was the body of `ga_maglev_test.cpp`) |
| `ga_test/maglev_utilities/` (whole dir) | `bench_dipole_ga_vs_classical.cpp`, `bench_maglev_stability.cpp`, `maglev_scene_replica.cpp` + its `CMakeLists.txt` (targets `ga_dipole_bench`, `ga_maglev_stability`, `maglev_scene_replica`) |
| `ga_view/src/active_maglev.{hpp,cpp}` | in `ga_view/` → moves with the whole module (§4c) |

Maglev docs move with the whole `TODO/` (principle 7) — no per-file listing needed
(`magnetic_levitation.md`, `maglev_S1_dipole_formulation.md`, `CtrlXFlow6d.pptx` are among
them; the first three are already git-excluded).

### 2c. GRINDING bundle → private (move only, no history purge, §7)

The grinding cases were **extracted from the shared 3dp test headers into a dedicated
`ga_grinding_test` bundle** (the intermediate step, now DONE). Files that move:

- `ga_test/src/ga_grinding_test.cpp` — entry (doctest main).
- `ga_test/src/ga_grinding_appl_test.hpp` — kinematics/geometry/topology cases (tumbling
  plane, intersecting discs, tumbling tool plane, grain trajectory Phase 0, wafer-tilt
  flatness Zhou, grinding-mark wavelengths D.1, wafer topography D.1b, contact force
  feed-forward D.2a, Sommerfeld warm-up B.1).
- `ga_test/src/ga_grinding_mechanics_test.hpp` — Tao wheel-spindle (C.2–C.4) + grinding
  force loop (D.2b), wafer thinning (D.2c), loop integrator (D.2d), feed control (D.2e).
- `ga_test/src/ga_grinding_cai_test.hpp` — the Cai volumetric-error machine (Phase F);
  folded into the grinding package as its own chapter (see [test_bundling_by_application.md](test_bundling_by_application.md)).
- View scenes: `ga_view/src/active_grinding_{cs,flatness,marks,topo}.{hpp,cpp}` (4 scenes) —
  these move with the **whole `ga_view/` module** (§4c), not individually.
- Python demos: `ga_py/demo/{tumbling_plane,intersecting_discs}.py` (grinding-labeled
  building blocks).
- Docs: the grinding write-ups (`grinding.md`, `tao_eq13_derivation.md`,
  `Grinding_setup_CS_draft.md`) move with the whole `TODO/` (principle 7); the grinding
  half of the `ga_docu` applications chapter moves whole (§4d).

> Note: the **force-control work (Phase D.2a–e) is already written and public** — it is not
> future work. It moves with the rest of the grinding bundle. **Cai now shares grinding's
> fate** (private) since it is a chapter of `ga_grinding_test`.

### 2d. Stays PUBLIC

- The whole `ga/` library, incl. force elements (`set_joint_spring_damper`,
  `set_applied_wrench`, `add_grounded_spring`, `grounded_spring{2,3}dp`), `dynamic_system`,
  `closed_loop_system` — general-purpose GA mechanics.
- `ga_prdxpr`, `ga_lua`, `ga_py` (minus the two grinding demos), core algebra tests,
  `bench_sta4ds_transform`.
- The generic mechanism/feature view scenes and their tests; the 2D application tests
  (`ga_appl2dp_*` — verified grinding-free).
- A **generated PDF** of the generic `ga_docu` chapters (intro, basics, modelling-motion,
  mechanics, comparison, glossary, literature). The LaTeX sources themselves are private
  and the PDF is built there and copied over (§4d) — public carries no `.tex`.

---

## 3. Link model (DECIDED): submodule superset

The private repo is a **superset that embeds the public repo as a git submodule**, pinned
to an exact public commit. Location: `/Users/hud3bh/prg/cpp/pj/ga_private` (sibling of
`ga`, mirroring the `industrial_ai` layout).

```text
ga/                              PUBLIC repo (github.com/Daniel-G-W-Hug/ga)
ga_private/                      PRIVATE repo (superset)
├── external/ga/                 ← git submodule → public repo, pinned to a commit
├── .clang-format                ← own copy at the root (its files search up to HERE, not
│                                  into external/ga); keep byte-identical to the public one
├── ga_test/src/     ga_maglev_* (dipole+twin+model), ga_grinding_* (appl/cai/mechanics)
├── ga_test/maglev_utilities/    bench_dipole, bench_maglev_stability, maglev_scene_replica
├── ga_test/private_targets.cmake        declares the private test/bench targets
├── ga_view/src/     active_maglev.*, active_grinding_*.*  (+ registration TUs)
├── ga_view/private_scenes.cmake         adds the private scene sources to ga_view
├── ga_docu/         6_ga_applications_pga.tex  (applications chapter)
├── TODO/            the WHOLE TODO/ directory (all plans, principle 7)
└── CMakeLists.txt   → add_subdirectory(external/ga); pass GA_PRIVATE_DIR down
```

**`.clang-format` in both repos (added 2026-07-08):** the repo style is now a committed
`.clang-format` at the public repo root (a verbatim copy of the maintainer's former global
config), so co-workers format identically regardless of their own global. The private repo
needs its **own** copy at `ga_private/` root — clang-format searches upward from a file and
would stop at `ga_private/` without descending into the `external/ga` submodule, so a
private-repo file would otherwise fall back to LLVM defaults. Keep the two copies
byte-identical.

The library-consuming pieces (tests, benches, replica) move trivially — they only
`target_link_libraries(... ga)` (maglev additionally includes the private
`ga_dipole_model.hpp`). With `ga_view` moving **wholesale** (§4c) there is no scene seam to
build; the remaining pieces are the `ga_test` overlay (§4b) and the grinding test extraction
(§4a, already done).

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

**Documentation** — all `ga_docu` LaTeX lives in `ga_private/`; public carries only a
generated generic-chapters PDF, refreshed from private when docs change (§4d). So editing
docs is always a private-repo activity; a doc update ends with rebuilding + copying the
public PDF into the `external/ga` submodule checkout and bumping it in public.

**Implementing private add-ons on public parts** — one working copy holds both: edit public
code inside `external/ga` (a normal public checkout on its own branch), private code in the
parent tree. Commits land in whichever repo owns the file; the public change is later
pinned in via a submodule bump.

---

## 4. Prerequisite public-repo refactors (do FIRST, in public)

Behaviour-preserving changes committed to public `develop` so the application bundles can
leave without breaking any public build.

### 4a. Extract the interleaved grinding cases into dedicated files — DONE

**Completed as the test-bundling intermediate step** (public refactor, committed to
`develop`; details in [test_bundling_by_application.md](test_bundling_by_application.md)).
Grinding was cut out of `ga_appl3dp_{appl,mechanics}_test.hpp` into the `ga_grinding_test`
bundle (`ga_grinding_{appl,cai,mechanics}_test.hpp`); the Cai machine folded in as its own
chapter; the maglev bundle consolidated (`ga_maglev_test` = dipole + twin chapters); the
maglev tools moved to `ga_test/maglev_utilities/` and the exporter to
`ga_test/python_utilities/`. All verified by assertion-count parity (pure relocation).

This established the **per-application file convention** (principle 5): each application is
its own `ga_<app>_test` executable (+ chapters) + `active_<app>_*` scenes + doc chapter.
Grinding/maglev are now whole-file/whole-directory moves for the repo split.

### 4b. `ga_test` — optional private target overlay

```cmake
# Opt-in private overlay: extra test/bench targets living outside this repo.
if(GA_PRIVATE_DIR AND EXISTS "${GA_PRIVATE_DIR}/ga_test/private_targets.cmake")
    include("${GA_PRIVATE_DIR}/ga_test/private_targets.cmake")
endif()
```

`private_targets.cmake` (in `ga_private`) declares `ga_maglev_test` and `ga_grinding_test`,
and includes the `maglev_utilities/` subdirectory (benches + replica) (linking
`doctest::doctest ga`, include dir = private root). Unset `GA_PRIVATE_DIR` → nothing added.

### 4c. `ga_view` — moved WHOLESALE to private (decided 2026-07-08; no seam)

**Superseding decision:** the entire `ga_view/` directory moves to the private repo — no
scene-by-scene split, no `coordsys_model`/`scenes.cpp` surgery, no registry, no `#ifdef`.
Rationale (user): the visualization explains GA aspects but is not self-explanatory on its
own, so the public repo does without it; and moving the whole module dissolves the delicate
scene seam entirely. The grinding *scenes* (`active_grinding_*`) and the maglev scene
(`active_maglev`) go private automatically as part of `ga_view` — nothing app-specific has
to be teased out of shared view files.

Mechanics:

- **Public top-level `CMakeLists.txt`:** guard the viewer so a public checkout (no `ga_view/`)
  skips it: `if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ga_view" AND Qt6_FOUND)
  add_subdirectory(ga_view) endif()`.
- **Private superset:** its top-level CMake adds `add_subdirectory(external/ga)` (the public
  library etc.) and then `add_subdirectory(ga_view)` from `ga_private/ga_view/`; `ga_view`
  links the `ga` target across the submodule boundary (added first, so the target exists).
- Nothing public build-depends on `ga_view` (verified) — removing it leaves the public build
  green. The only public reference is a docstring in
  `ga_py/tests/test_merry_go_round.py` (it reconstructs the kinematics independently) →
  reword the comment.

This retires the former scene-registry seam design entirely (the wholesale move is simpler
and lower-risk).

### 4d. `ga_docu` — all LaTeX private; public gets a generated PDF only (decided 2026-07-08)

**Simplification (user, 2026-07-08):** the public repo keeps **no LaTeX at all** — the
entire `ga_docu/` `.tex` sources + build machinery move to the private repo, and the public
repo carries a **single generated PDF** of the *generic* chapters (the public library
documentation), refreshed from private whenever the docs change. This removes the whole
LaTeX seam from public (no `\IfFileExists` chapter guard, no glossary-`\ref` reword, no
cross-ref management, no LaTeX toolchain for public users).

Layout:

- **Private** owns `ga_docu/*.tex` (all chapters, incl. the applications chapter) + the
  build. It builds **two** PDFs:
  - the *full* internal PDF (all chapters, incl. applications — private only), and
  - the *public* PDF = generic chapters only (applications chapter omitted). The omit is
    a private-side build switch — e.g. keep the `\IfFileExists{6_ga_applications_pga.tex}`
    guard in `0_ga_docu.tex` and build the public variant with that file absent, or a
    `\if@public` flag. The one dangling `8_ga_glossary.tex` `\ref{subsubsection:wafer_grinding}`
    is reworded once (private-side) so the public variant has no undefined reference.
- **Public** carries only `ga_docu/ga_docu.pdf` (the generic-chapters PDF) — no `.tex`.
- **Sync workflow:** a small private script (`make public-pdf`) rebuilds the public PDF and
  copies it to `external/ga/ga_docu/ga_docu.pdf` (the submodule checkout), committed to
  public as a normal artifact bump when docs change.

**Public-repo removal (in §6):** `git rm ga_docu/*.tex` and the LaTeX aux/build files; keep
(or replace with) the generic-only PDF. **History:** the current public `.tex` (applications
chapter) *and* the current committed PDF already render maglev, so both the applications
`.tex` and the old `0_ga_docu.pdf` blobs fall under the maglev history purge (§7) — add
`ga_docu/6_ga_applications_pga.tex` and `ga_docu/0_ga_docu.pdf` to the `filter-repo --path`
set; the public repo then re-adds a fresh generic-only PDF at a clean name (e.g.
`ga_docu/ga_docu.pdf`).

### 4e. Confirm the clean cut

```bash
grep -rln -iE 'maglev|dipole|grind|wafer|tao|sommerfeld|zhou' \
    ga/ ga_prdxpr/ ga_lua/ ga_py/
# after 4a–4d + moves: expect no matches (ga_docu is gone from public entirely)
```

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
   the **whole `ga_docu/` LaTeX** (§4d). Get the private superset building; all application
   tests, benches, and view scenes pass **from the private repo**; the private doc build
   produces both the full and the public (generic-only) PDF.
5. **Remove from `ga`** (plain `git rm`): the maglev files, the grinding files
   (`ga_grinding_*`, `active_grinding_*`, the two Python demos), and **all `ga_docu/*.tex` +
   LaTeX build files** (replace with the generic-only `ga_docu/ga_docu.pdf` from private);
   delete their references from public CMake / `coordsys_model` (now behind the registry).
   Verify build + tests green; `ga_view` shows only generic demos.
6. **History: purge MAGLEV (§7) — while private.** Grinding is move-only, no purge.
7. **Move the WHOLE `TODO/` directory** into `ga_private/TODO/` (principle 7 — all plans go
   private) and `git rm -r TODO/` from public; drop the now-moot `.git/info/exclude`
   entries. (Decide separately whether the tracked plan docs also need a history purge, or
   just leave older history — see §7.)
8. **Update docs**: public `CLAUDE.md` (the split + per-app file convention; **remove/
   genericise every `TODO/*` reference** since `TODO/` is now private — principle 7; drop
   the `ga_view` build/scene sections). Public **top-level `README.md`: remove the 3
   `ga_view` lines** (module list, the `ga_view: fmt, qt6` dependency, the `brew install
   qt6` line) — audited 2026-07-08, that is the *only* README cross-reference to a
   private-bound module; the `ga_prdxpr` README section **stays** (public), and all
   per-module READMEs are already clean. The branch-model note (`github_setup.md`) is
   itself a `TODO/` doc → private. Also reword the `ga_view` docstring in
   `ga_py/tests/test_merry_go_round.py`.
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
   `ga_test/src/ga_dipole_model.hpp`, `ga_test/src/ga_maglev_test.cpp`,
   `ga_test/src/ga_maglev_dipole_test.hpp`, `ga_test/src/ga_maglev_twin_test.hpp`, the whole
   `ga_test/maglev_utilities/` directory (the two benches + the replica),
   `ga_view/src/active_maglev.{hpp,cpp}`, `ga_docu/6_ga_applications_pga.tex` (the app
   chapter contains the confidential maglev section; purging it also drops the
   non-confidential grinding doc from history — harmless), and **`ga_docu/0_ga_docu.pdf`**
   (every committed build of it *renders* the maglev chapter, so the binary blobs carry the
   confidential content too — purge the PDF path, then the public repo re-adds a fresh
   generic-only PDF at a clean name, §4d). Note the pre-bundling `ga_dipole_test.cpp` is also
   present in older history and must be purged too (git-log the path to confirm and add it to
   the filter-repo `--path` set).
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

> **`TODO/` (whole directory, principle 7):** the git-excluded docs
> (`magnetic_levitation.md`, `maglev_S1_dipole_formulation.md`, `CtrlXFlow6d.pptx`) were
> never committed (confirm with `git log --all --oneline -- TODO/magnetic_levitation.md` →
> empty), so they need no purge. The *tracked* plan docs (`grinding.md`,
> `public_private_repo.md`, `github_setup.md`, …) are in public history; `git rm -r TODO/`
> removes them from `HEAD`, and — since they are plans, not confidential physics — leaving
> them in older history is acceptable. Purge them too only if you want the public repo to
> carry no plan history at all (they can be added to the maglev `filter-repo --path` set).
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
  application targets; `ga_maglev_test`, `ga_grinding_test`, and the `maglev_utilities/`
  tools (`maglev_scene_replica`) pass; maglev + grinding scenes appear in `ga_view` and match
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
