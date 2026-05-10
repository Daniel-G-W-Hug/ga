# ga_bindgen — Binding Generator for `ga/`

> **For contributors regenerating Python bindings.** End users only need
> [`ga_py/README.md`](../ga_py/README.md) — generated bindings are committed to git, so
> building or using `ga_py` does not require this tool.

Tool that scans the C++ headers in `ga/` with libclang and emits Python (nanobind)
bindings for `ga_py/`. Single source of truth: the C++ headers themselves.

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

## First-time setup

The generator runs in its own venv (`ga_bindgen/.venv`) — separate from the wrapper venv
at `ga_py/.venv`. See [`ga_py/README.md
§5.1`](../ga_py/README.md#51-the-two-virtual-environments--dont-mix-them) for the
rationale and side-by-side activation commands.

> **Working directory throughout: the project root** (the directory containing `ga/`,
> `ga_py/`, `ga_bindgen/`, `CMakeLists.txt`). Every command below uses paths relative to
> that root.

### 1. Install LLVM / libclang on the system

**macOS:**

```bash
brew install llvm
```

**Linux (Debian / Ubuntu):**

```bash
apt install libclang-dev    # or the distro equivalent
```

**Windows:** No separate LLVM install is required. The pip `libclang` package (installed
in step 2) bundles `libclang.dll` and is sufficient for scanning the `ga/` headers.
Optionally install the [LLVM release](https://releases.llvm.org/download.html) pre-built
installer if you want the full clang toolchain; `clang_setup.py` will prefer the system
install over the bundled DLL when it finds one at `C:\Program Files\LLVM\lib\libclang.dll`.

### 2. Create the bindgen venv and install Python deps

**macOS / Linux** — cwd: project root:

```bash
python3 -m venv ga_bindgen/.venv
ga_bindgen/.venv/bin/pip install libclang
```

**Windows** (Developer Command Prompt or PowerShell) — cwd: project root:

```bat
python -m venv ga_bindgen\.venv
ga_bindgen\.venv\Scripts\pip install libclang
```

`ga_bindgen/.venv` is gitignored. The `nanobind` package is **not** needed here — it's
only required at C++ compile time for the generated bindings.

### 3. Verify the setup

A quick sanity check that libclang loads and `clang_setup.py` can find a usable resource
dir:

**macOS / Linux** — cwd: project root:

```bash
ga_bindgen/.venv/bin/python ga_bindgen/src/clang_setup.py
```

**Windows** (Developer Command Prompt or PowerShell) — cwd: project root:

```bat
ga_bindgen\.venv\Scripts\python.exe ga_bindgen/src/clang_setup.py
```

Expected output: a path to a `libclang.{dylib,so,dll}` followed by the assembled compile
flags. If this raises `RuntimeError: Could not locate a usable libclang`, recheck step 1
and 2 — the candidate paths the script probes are listed under
[Platform notes](#platform-notes) below.

## Regenerating the manifest and bindings

Once the venv is set up, regeneration is three steps run from the project root.

**macOS / Linux** — cwd: project root:

```bash
source ga_bindgen/.venv/bin/activate

# Step 1 — scan headers, write manifest
python3 ga_bindgen/src/scan.py

# Step 2 — emit nanobind glue from manifest (full regeneration is the default)
python3 ga_bindgen/src/emit_nanobind.py

# Step 3 — recompile (back in the wrapper venv or directly)
cmake --build build
```

**Windows (Developer Command Prompt)** — recommended Windows shell; cwd: project root.
Open via Start → "Developer Command Prompt for VS 2022".

```bat
ga_bindgen\.venv\Scripts\activate.bat

rem Step 1 — scan headers, write manifest
python ga_bindgen\src\scan.py

rem Step 2 — emit nanobind glue from manifest
python ga_bindgen\src\emit_nanobind.py

rem Step 3 — recompile
cmake --build build --target _ga_py --config Debug
```

**Windows (PowerShell)** — cwd: project root:

```powershell
& ga_bindgen\.venv\Scripts\Activate.ps1

# Step 1 — scan headers, write manifest
python ga_bindgen/src/scan.py

# Step 2 — emit nanobind glue from manifest (full regeneration is the default)
python ga_bindgen/src/emit_nanobind.py

# Step 3 — recompile
cmake --build build --target _ga_py --config Debug
```

> **PowerShell note:** if running `Activate.ps1` opens a file-association dialog,
> prefix it with `&` as shown above. If the execution policy blocks it, run
> `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once first.

Both scripts take `--help`. `emit_nanobind.py` regenerates every binding by default; pass
`--type <name>` (e.g. `--type vec3d`) for single-type emission during binding-script
development. The legacy `--all` flag is still accepted as a no-op alias for the default.

`ga_py/CMakeLists.txt` also exposes a `regenerate_python_bindings` umbrella target that
drives both scripts; it auto-detects `ga_bindgen/.venv` if present, otherwise falls back
to whatever `find_package(Python)` resolved.

## Platform notes

`clang_setup.py` probes these paths in order; the first hit wins:

- **macOS:** `/opt/homebrew/opt/llvm/lib/libclang.dylib` (Apple Silicon) or
  `/usr/local/opt/llvm/lib/libclang.dylib` (Intel)
- **Linux:** `/usr/lib/x86_64-linux-gnu/libclang-{18,17,16,15,14}.so.1`,
  `/usr/lib/llvm-{18,17,16,15,14}/lib/libclang.so.1`, then
  `/usr/lib/x86_64-linux-gnu/libclang.so`
- **Windows:** `C:\Program Files\LLVM\lib\libclang.dll`,
  `C:\Program Files (x86)\LLVM\lib\libclang.dll`, then the bundled DLL inside the pip
  `libclang` package at `ga_bindgen\.venv\Lib\site-packages\clang\native\libclang.dll`.

On Windows, MSVC and Windows SDK include paths are discovered automatically via
`vswhere.exe` — no manual path configuration needed when Visual Studio 2022 is installed.

If your install lives elsewhere (macOS / Linux), either symlink it into one of those
paths or extend the candidate list in [`src/clang_setup.py`](src/clang_setup.py).

## Dependencies

- Python 3.10+
- `libclang` Python package (installed into `ga_bindgen/.venv` per step 2 above)
- **macOS / Linux:** A real LLVM install on the system (per step 1) — the pip package's
  bundled resource dir is not sufficient on those platforms
- **Windows:** The pip `libclang` package's bundled `libclang.dll` is sufficient; Visual
  Studio 2022 (for MSVC headers) must be installed
- `nanobind` is only needed at C++ compile time for the generated bindings, not by
  `ga_bindgen` itself

## Design notes

- **Generated files committed to git**: regeneration is an explicit step run by
  contributors, not an implicit build dependency. PR diffs show exactly what new Python
  surface a C++ change created.
- **Filter:** only symbols declared inside `ga/` are collected; system headers, `fmt`,
  etc. are excluded. The allow-list in `scan.py` is `{"hd::ga", "hd::ga::ega",
  "hd::ga::pga"}` — `hd::ga::detail` is excluded by design (see [`ga_py/README.md
  §6.8`](../ga_py/README.md#68-excluded-namespace-hdgadetail)).
- **Value type frozen at `value_t = double`** in v1.

## Open items / future work

These are tooling features that have not landed yet. Wrapper-side open work lives in
[`ga_py/README.md §6`](../ga_py/README.md#6-open-items--future-work).

### `audit_lua.py` — per-overload Lua coverage audit

`ga_lua/src/ga_lua.hpp` is hand-maintained and known incomplete. The plan is to invert the
relationship — let `ga_bindgen` produce the authoritative manifest and audit Lua against
it.

A 2026-04 string-search spike confirmed gaps (`rk4_step`, `rk4_get_time`, `Hz2radps`,
`radps2Hz`, `rpm2radps`, `radps2rpm`, `*_metric_view` helpers). A per-overload audit
driven by the manifest is the next step, with output landing in
`ga_bindgen/reports/lua_coverage_*.md`. The Python wrapper does not depend on it.

### Stub generation `emit_stubs.py`

Emit `.pyi` files alongside the `bindings_*.cpp` so IDEs and `mypy --strict` can see the
Python surface. See [`ga_py/README.md
§6.2`](../ga_py/README.md#62-stub-generation-pyi-and-mypy) for the user-visible
motivation.

### Manifest delta report

A `reports/delta_YYYY-MM-DD.md` showing "new since last run", "removed since last run",
"changed signature" diffs would help PR reviewers. Low priority — `git diff` on
`manifest.json` mostly serves the same purpose.
