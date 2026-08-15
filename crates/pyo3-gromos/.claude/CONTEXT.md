# pyo3-gromos — stage contract

## Job
L4 Python bindings (PyO3). A thin wrapper over the Rust core — no physics, no data duplication.

## Inputs (consumes from)
All Rust crates.

## Outputs
Python-callable API for running simulations and analysing trajectories.

## Status
- Compositional Simulation API ✓, AlgorithmSequence API ✓
- `System` + `InputParameters` factories + 2-arg `Simulation(system, params)` ✓ (P3.1/3.2)
- `sim.run(steps, ene_freq)` → numpy array; energy decomposition (bond/angle/dihedral/improper/lj/coulomb) ✓ (P3.3)
- Python tests: 121 (118 passed, 2 documented position-mismatch skips) ✓ (`py-gromos/tests/`)
- `.pyi` stubs ✓
- **PLAN.md 3.6 ✓ done** — `build_simulation()` (`simulation.rs`) now dispatches every gromosXX
  feature that was already implemented lower in the stack but silently unwired from Python: NTIVEL=1
  velocity generation, SETTLE/LINCS (`ConstraintSelection::from_imd`), Nosé-Hoover (+chain)
  thermostat (`push_thermostat()`), distance/position restraints (`distrest=`/`posresspec=`/
  `refpos=` kwargs), triclinic/truncated-octahedron box (`NTB=-1`, plus `sim.forces` frame
  rotation). `REFERENCE_SYSTEMS` (`py-gromos/tests/test_gromosXX_references.py`) now covers all 37
  active Rust reference systems except the 2 explicitly-deferred FEP ones. **Known, tracked
  divergence:** `gromos-rs`'s own `.trc` position output disagrees with gromosXX for
  triclinic/EM-shift-frame systems (pre-existing, Rust-side, not a Python binding bug) — see
  `POSITION_MISMATCH_SYSTEMS` in the test file.
- **PLAN.md 3.5 ✓ done (M1, M2, SD)**:
  - **M1**: `InputParameters.nve/nvt/npt` take `constraints="none"|"hbonds"|"allbonds"` (default
    `"none"`) → `ntc=1|2|3` (`ntc_from_constraints()`, `gromos-io/imd.rs`); readable `.constraints`
    getter. No builder change needed — the existing `shake_enabled = imd.ntc > 1 || …` guard picks it
    up. Verified empirically: factory `nvt(..., constraints="hbonds")` on `aladip_solvated` stays
    stable at `dt=0.002`; `constraints="none"` on the same system still diverges (documented
    contrast, not a regression). Test: `test_constrained_system_stable_with_factory_params`.
  - **M2**: `temperature` getter now uses `PySimulation.total_dof` (constraint-aware, computed once
    via the shared `compute_total_dof()` in `algorithm_sequence.rs` — same value both builders and
    the thermostat use) instead of bare `n_atoms*3`.
  - **SD**: `build_simulation()` branches on `imd.ntem > 0` — pushes `SteepestDescentAlgorithm`
    instead of leap-frog, skips thermostat/barostat/`TemperatureCalculation`/COM removal (mirrors
    the `md` binary's own EM sequence). New `SteepestDescent` algorithm-sequence building block +
    `AlgorithmSequence.minimize()` preset; `from_parameters()` dispatches to it when `ntem > 0` — the
    composable path and `Simulation(system, params)` now agree. `shake_algorithm_from_imd()`
    extracted so EM and standard-MD SHAKE construction can't drift apart. Verified: `aladip_vacuum`
    under `steepest_descent(steps=50)` actually drops -13.5→-47.98 kJ/mol (previously a silent no-op
    at `dt=0`). Tests: `test_steepest_descent_via_simulation`, `test_steepest_descent_via_algorithm_sequence`.
- `sim.volume`/`sim.pressure` getters (mirror `temperature`'s shape). **`pressure` is only
  physically meaningful under NPT** — the virial term is only populated by `PressureCalculation`,
  which only NPT's sequence includes; under NVE/NVT it returns the misleading kinetic-only term, not
  zero. Documented in the getter, `.pyi`, API reference, and `02_short_md.ipynb`. Test:
  `test_volume_and_pressure_getters`.
- Remaining P3 items:
  - [ ] Expose ForceField evaluation (single-point energy/force)
  - [ ] Rich `__repr__` / `_repr_html_` for Jupyter (Topology, Configuration, Energy)
  - [ ] PLAN.md 3.7 — ML potential binding (`PySchNetPotential`, `Simulation.add_ml_potential`,
    selection-string zones). Blocked on `gromos-forces` P2.8-6 (no orchestrator-aware step path
    exists to bind yet) and on `--features ml` needing `libtorch` (not available here).
- Deferred (tech debt, not scheduled — full audit in `~/.claude/plans/golden-baking-liskov.md`):
  unify the two imd→sequence builders (`build_simulation` vs `resolve_algorithm_sequence`),
  consolidate scattered defaults, constraints-as-`System`-attribute, symmetric `InputParameters`
  construction (kwargs ctor + setters + `write_imd_file`/`.save()`).

## Crate-specific rules
- **Thin wrapper only.** Zero physics, zero data structures that duplicate the Rust core.
- **API design reference:** study Polars pyo3 patterns (`PyDataFrame` / `PyExpr` / `PyLazyFrame`) at `.local/polars/py-polars/src/`.
- `py-gromos` (the Python package, separate maturin build) sits alongside this crate; `maturin develop` must pass after any change here.
