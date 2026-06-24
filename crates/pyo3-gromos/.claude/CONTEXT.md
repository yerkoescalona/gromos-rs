# pyo3-gromos — stage contract

## Job
L4 Python bindings (PyO3). A thin wrapper over the Rust core — no physics, no data duplication.

## Inputs (consumes from)
All Rust crates.

## Outputs
Python-callable API for running simulations and analysing trajectories.

## Status
- Compositional Simulation API ✓, AlgorithmSequence API ✓
- Python reference tests: 62 passed ✓
- `.pyi` stubs ✓
- Remaining P3 items:
  - [ ] Expose ForceField evaluation (single-point energy/force)
  - [ ] Expose SHAKE / constraint info
  - [ ] Expose energy decomposition (bonded, LJ, CRF, kinetic, pressure)
  - [ ] Method chaining: `sim.run(steps=1000).energies().plot()`
  - [ ] Energy timeseries as DataFrame (Polars/pandas interop)
  - [ ] Rich `__repr__` / `_repr_html_` for Jupyter (Topology, Configuration, Energy)

## Crate-specific rules
- **Thin wrapper only.** Zero physics, zero data structures that duplicate the Rust core.
- **API design reference:** study Polars pyo3 patterns (`PyDataFrame` / `PyExpr` / `PyLazyFrame`) at `.local/polars/py-polars/src/`.
- `py-gromos` (the Python package, separate maturin build) sits alongside this crate; `maturin develop` must pass after any change here.
