# py-gromos — stage contract

## Job
L4+ user-facing Python package. The importable `gromos` package: Python wrappers, high-level
simulation runners, analysis helpers, notebooks, and examples. Built with `maturin` on top of
`pyo3-gromos`. Contains no Rust code — only Python that calls the compiled extension.

## Inputs (consumes from)
`pyo3-gromos` (compiled `.so` extension, via `maturin develop` or installed wheel).

## Outputs
- `python/gromos/` — importable Python package
  - `__init__.py` — top-level API surface
  - `md_runners.py` — high-level simulation runner wrappers
  - `analysis.py` — Python-side analysis helpers (calls gromos-analysis via pyo3)
  - `gromos.pyi` — type stubs for IDE support
- `notebooks/` — Jupyter notebooks (education + demonstration)
- `examples/` — standalone Python scripts (17 scripts)
- `tests/` — Python-level integration tests

## Status
- Basic wrappers ✓; `maturin develop` builds ✓
- Remaining P3 items:
  - [ ] `md_runners.py` simplify; `analysis.py` expose gromos-analysis to Python
  - [ ] Method chaining: `sim.run(steps=1000).energies().plot()`
  - [ ] Energy timeseries as DataFrame (Polars/pandas interop)
  - [ ] Rich `__repr__` / `_repr_html_` for Jupyter (Topology, Configuration, Energy)
  - [ ] Rewrite `notebooks/`: 01 inspect + single-point energy; 02 short MD + energy conservation; 03 NVE/NVT/NPT comparison
  - [ ] Rewrite `examples/` (17 scripts) on the new API
  - [ ] Fix `test_basic.py`, `test_advanced_features.py`

## Key files
```
python/gromos/__init__.py    — top-level API
python/gromos/md_runners.py  — simulation runner wrappers
python/gromos/analysis.py    — analysis helpers
python/gromos/gromos.pyi     — type stubs
notebooks/                   — Jupyter notebooks
examples/                    — standalone scripts
pyproject.toml               — package metadata (maturin)
```

## Crate-specific rules
- **No physics, no data structures** — this layer only wraps what pyo3-gromos exposes.
- **API design reference:** model method chaining and DataFrame interop on Polars patterns (`.local/polars`).
- Changes here require `maturin develop` in this directory to rebuild the extension before testing.
- `gromos.abi3.so` is a compiled artifact — never edit it directly; rebuild via maturin.
