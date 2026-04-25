# ga_bindgen — Binding Generator for `ga/`

Tool that scans the C++ headers in `ga/` with libclang and emits Python (nanobind) bindings for `ga_py/`. Single source of truth: the C++ headers themselves.

Status: **end-to-end pipeline proof** — currently emits bindings only for `vec3d` to validate the full scanner → emitter → CMake → Python import flow. Expanding to all 39 user types is mechanical from here.

See `TODO/considerations_python_wrapper.md` for the full design context.

## Components

```
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
- `libclang` Python package (probes for system/homebrew libclang at runtime; pip-bundled libclang's resource dir is broken on macOS)
- `nanobind` (only needed by the generated code at compile time, not by ga_bindgen itself)

## Design notes

- **Generated files committed to git**: regeneration is an explicit step run by contributors, not an implicit build dependency. PR diffs show exactly what new Python surface a C++ change created.
- **Filter:** only symbols declared inside `ga/` are collected; system headers, `fmt`, etc. are excluded.
- **Value type frozen at `value_t = double`** in v1.
