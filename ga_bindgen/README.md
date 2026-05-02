# ga_bindgen — Binding Generator for `ga/`

> **For contributors regenerating Python bindings.** End users only need [`ga_py/README.md`](../ga_py/README.md) — generated bindings are committed to git, so building or using `ga_py` does not require this tool.

Tool that scans the C++ headers in `ga/` with libclang and emits Python (nanobind) bindings for `ga_py/`. Single source of truth: the C++ headers themselves.

## Components

```text
ga_bindgen/
├── src/
│   ├── clang_setup.py     # find_libclang() + compile flags
│   ├── model.py           # manifest dataclasses (schema_version=1)
│   ├── scan.py            # ga/*.hpp → manifest.json
│   └── emit_nanobind.py   # manifest.json → ga_py/src/generated/*.cpp
├── manifest.json          # last-known API snapshot (committed)
└── README.md
```

## Usage

The generator runs in its own venv (`build/spike_libclang/.venv`) — separate from the wrapper venv. See [`ga_py/README.md §5.1`](../ga_py/README.md#51-the-two-virtual-environments--dont-mix-them) for the rationale and side-by-side activation commands.

```bash
source build/spike_libclang/.venv/bin/activate

# Step 1 — scan headers, write manifest
python3 ga_bindgen/src/scan.py

# Step 2 — emit nanobind glue from manifest (full regeneration is the default)
python3 ga_bindgen/src/emit_nanobind.py

# Step 3 — recompile (back in the wrapper venv or directly)
cmake --build build
```

Both scripts take `--help`. `emit_nanobind.py` regenerates every binding by
default; pass `--type <name>` (e.g. `--type vec3d`) for single-type emission
during binding-script development. The legacy `--all` flag is still accepted
as a no-op alias for the default.

## Dependencies

- Python 3.10+
- `libclang` Python package — wraps the native libclang shared library. Install with `pip install libclang`. The pip package ships the shared library on Linux and Windows; on macOS the bundled resource dir is broken, so a system/Homebrew LLVM is required (see below).
- `nanobind` — only needed at C++ compile time for the generated bindings, not by `ga_bindgen` itself.

### Platform notes

**macOS:** Install LLVM via Homebrew so `clang_setup.py` can find a complete resource dir:

```bash
brew install llvm
pip install libclang
```

**Linux:**

```bash
apt install libclang-dev   # or the distro equivalent
pip install libclang
```

**Windows:** Install [LLVM](https://releases.llvm.org/download.html) (the pre-built Windows installer); `clang_setup.py` probes `C:\Program Files\LLVM\lib\libclang.dll` automatically. Then create the venv:

```bat
python -m venv build\spike_libclang\.venv
build\spike_libclang\.venv\Scripts\pip install libclang
```

## Design notes

- **Generated files committed to git**: regeneration is an explicit step run by contributors, not an implicit build dependency. PR diffs show exactly what new Python surface a C++ change created.
- **Filter:** only symbols declared inside `ga/` are collected; system headers, `fmt`, etc. are excluded. The allow-list in `scan.py` is `{"hd::ga", "hd::ga::ega", "hd::ga::pga"}` — `hd::ga::detail` is excluded by design (see [`ga_py/README.md §6.8`](../ga_py/README.md#68-excluded-namespace-hdgadetail)).
- **Value type frozen at `value_t = double`** in v1.

## Open items / future work

These are tooling features that have not landed yet. Wrapper-side open work lives in [`ga_py/README.md §6`](../ga_py/README.md#6-open-items--future-work).

### `audit_lua.py` — per-overload Lua coverage audit

`ga_lua/src/ga_lua.hpp` is hand-maintained and known incomplete. The plan is to invert the relationship — let `ga_bindgen` produce the authoritative manifest and audit Lua against it.

A 2026-04 string-search spike confirmed gaps (`rk4_step`, `rk4_get_time`, `Hz2radps`, `radps2Hz`, `rpm2radps`, `radps2rpm`, `*_metric_view` helpers). A per-overload audit driven by the manifest is the next step, with output landing in `ga_bindgen/reports/lua_coverage_*.md`. The Python wrapper does not depend on it.

### Stub generation `emit_stubs.py`

Emit `.pyi` files alongside the `bindings_*.cpp` so IDEs and `mypy --strict` can see the Python surface. See [`ga_py/README.md §6.2`](../ga_py/README.md#62-stub-generation-pyi-and-mypy) for the user-visible motivation.

### Manifest delta report

A `reports/delta_YYYY-MM-DD.md` showing "new since last run", "removed since last run", "changed signature" diffs would help PR reviewers. Low priority — `git diff` on `manifest.json` mostly serves the same purpose.
