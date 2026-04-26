# ga_bindgen — Binding Generator for `ga/`

Tool that scans the C++ headers in `ga/` with libclang and emits Python (nanobind) bindings for `ga_py/`. Single source of truth: the C++ headers themselves.

Status: **end-to-end pipeline proof** — currently emits bindings only for `vec3d` to validate the full scanner → emitter → CMake → Python import flow. Expanding to all 39 user types is mechanical from here.

See `TODO/considerations_python_wrapper.md` for the full design context.

## Two virtual environments — don't mix them

This project keeps the binding generator and the wrapper in separate venvs. They serve different roles, have different dependencies, and must not be conflated:

| Venv path | Used for | Dependencies |
| --- | --- | --- |
| `build/spike_libclang/.venv` | **Regenerating bindings via `ga_bindgen`** (this tool). Only contributors touching `ga/*.hpp` or this generator need it. | `libclang` |
| `ga_py/.venv` | **Running the wrapper and its tests.** Everyday development of `ga_py/`. | `pytest`, `hypothesis`, `numpy` |

If you run `python ga_bindgen/src/scan.py` from the wrong venv, libclang won't be importable. If you run `pytest ga_py/tests/` from the wrong venv, `pytest` won't be on PATH (or `hypothesis` will be missing).

**macOS / Linux:**

```bash
# For regenerating bindings (this tool):
source build/spike_libclang/.venv/bin/activate

# For wrapper development / running tests (separate session):
source ga_py/.venv/bin/activate
```

**Windows** (`Scripts\` instead of `bin/`):

```bat
build\spike_libclang\.venv\Scripts\activate

rem separate session for wrapper tests:
ga_py\.venv\Scripts\activate
```

## Components

```text
ga_bindgen/
├── src/
│   ├── clang_setup.py     # find_libclang() + compile flags
│   ├── model.py           # manifest dataclasses (schema_version=1)
│   ├── scan.py            # ga/*.hpp → manifest.json
│   └── emit_nanobind.py   # manifest.json → ga_py/src/generated/*.cpp
├── manifest.json          # last-known API snapshot (committed once stable)
└── README.md
```

## Usage

Run from the project root with a venv that has `libclang` and `nanobind` installed:

```bash
# Step 1 — scan headers, write manifest
python3 ga_bindgen/src/scan.py

# Step 2 — emit nanobind glue from manifest
python3 ga_bindgen/src/emit_nanobind.py
```

Both scripts take `--help`.

## Dependencies

- Python 3.10+
- `libclang` Python package — wraps the native libclang shared library. Install with `pip install libclang`. The pip package ships the shared library on Linux and Windows; on macOS the bundled resource dir is broken, so a system/Homebrew LLVM is required (see below).
- `nanobind` — only needed at C++ compile time for the generated bindings, not by ga_bindgen itself.

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

**Windows:** Install [LLVM](https://releases.llvm.org/download.html) (the pre-built Windows installer); `clang_setup.py` probes `C:\Program Files\LLVM\lib\libclang.dll` automatically. Then:

```bat
pip install libclang
```

Create the venv and install:

```bat
python -m venv build\spike_libclang\.venv
build\spike_libclang\.venv\Scripts\pip install libclang
```

## Design notes

- **Generated files committed to git**: regeneration is an explicit step run by contributors, not an implicit build dependency. PR diffs show exactly what new Python surface a C++ change created.
- **Filter:** only symbols declared inside `ga/` are collected; system headers, `fmt`, etc. are excluded.
- **Value type frozen at `value_t = double`** in v1.
