# ga_py — Python wrapper around the `ga/` library

`ga_py` exposes the geometric algebra library as a Python module. Types (`vec3d`, `bivec3d`, `mvec3d`, …), operators (`+`, `-`, `*`, `^`, `<<`, `>>`), and free functions (`gpr`, `wdg`, `dot`, `nrm`, …) all map 1:1 to their C++ counterparts.

```python
import ga_py

v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
B = ga_py.ega.bivec3d(0.0, 0.0, 1.0)
m = v + B                                    # cross-grade addition → mvec3d
print(m, m.gr1(), m.gr2())                   # extract grade components
```

The package is **generated** from the C++ headers by `ga_bindgen/` — see [TODO/considerations_python_wrapper.md](../TODO/considerations_python_wrapper.md) for the design and [ga_bindgen/README.md](../ga_bindgen/README.md) for the toolchain.

## Performance — what's compiled, what's interpreted

A common worry about Python wrappers is that the math will be re-implemented in interpreted Python and pay an interpreter tax per operation. **That is not the case here.**

When Python evaluates `v1 + v2`:

```
Python:  v1 + v2
            │
            ▼  nanobind dispatch — type check + unwrap (one C call, ~few hundred ns)
C++:     [](vec3d const& a, vec3d const& b) { return a + b; }
            │
            ▼  compiled C++ — ADL resolves to operator+ for Vec3_t<double, vec3d_tag>
C++:     hd::ga::operator+(...)
            │
            ▼  inlined and optimized exactly like a normal C++ call
            
Result wrapped back into a Python `vec3d` object.
```

The lambdas in [src/generated/](src/generated/) are real C++ code, compiled with full optimization, instantiating the GA templates for `value_t = double`. There is **no Python re-implementation of the math** anywhere. The only per-call overhead is the Python ↔ C nanobind dispatch, which is unavoidable for any Python binding regardless of how the C++ is packaged.

### "Why isn't there a separate `libga_compiled.a` that the wrapper links against?"

It would not change runtime speed. The same template instantiations would run with the same machine code whether they live inside `_ga_py.so` directly (today) or in a separate static archive (alternative). The only thing such a library can improve is **compile time**: today every `bindings_*.cpp` re-instantiates the same `Vec3_t<double, vec3d_tag>` etc., and a precompiled archive could do the work once. For the current ~50 binding files compiling in ~30s, this is not yet painful, so the precompiled-archive optimization is deferred. See [TODO/considerations_python_wrapper.md §10b](../TODO/considerations_python_wrapper.md) for the location-of-such-a-library question (which would land in `ga_py/` first if it ever does).

### Where Python *would* be slower than raw C++

Two places — both unrelated to whether the math itself is "compiled":

1. **Per-call dispatch overhead.** Each Python → C++ crossing costs a few hundred ns. For million-call inner loops on tiny types (e.g., adding millions of `vec3d`s in a Python `for` loop), this dominates. The fix is **batching** (e.g., numpy-array-of-multivectors operations), which is the Tier 2 numpy work deferred per the design doc.
2. **Object lifetime overhead.** Every result is a new Python object with a refcounted handle. Negligible for normal use; adds up in tight loops.

Neither is fixed by precompiled libraries; both are fixed (when needed) by vectorized array APIs.

## Layout

```
ga_py/
├── CMakeLists.txt              # nanobind + scikit-build target
├── README.md                   # this file
├── src/
│   ├── module.cpp              # hand-written entry point (NB_MODULE)
│   └── generated/              # produced by ga_bindgen — do not edit
│       ├── bindings_*.cpp      # one file per bound user type
│       ├── register_all.cpp    # dispatcher in inheritance-topological order
│       └── bindings_list.cmake # generated source list, included by CMakeLists.txt
├── python/
│   └── ga_py/__init__.py       # re-exports ega/pga submodules
└── tests/                      # pytest
```

## Building

From the project root, with a venv that has `nanobind` and `pytest` installed:

```bash
NANOBIND_DIR=$(python3 -c 'import nanobind; print(nanobind.cmake_dir())')
mkdir -p build/ga_py
cmake -S ga_py -B build/ga_py \
      -DPython_EXECUTABLE=$(which python3) \
      -DCMAKE_PREFIX_PATH="$NANOBIND_DIR" \
      -GNinja
cmake --build build/ga_py
```

Run tests:

```bash
PYTHONPATH="build/ga_py:ga_py/python" pytest ga_py/tests/
```

## Regenerating bindings after C++ changes

When `ga/` gains new types, operators, or functions, regenerate from the headers:

```bash
python3 ga_bindgen/src/scan.py        # re-scan headers → manifest.json
python3 ga_bindgen/src/emit_nanobind.py --all
cmake --build build/ga_py             # recompile
```

Generated files are committed to git so PR diffs show exactly what new Python surface a C++ change exposes.
