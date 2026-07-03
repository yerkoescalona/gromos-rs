# Learning Guide

A guide to understanding the py-gromos codebase and where things are headed.

## What exists and what is stub

The cleanest way to understand the current state is to look at `__init__.py`:
it imports only what actually works from the Rust extension. Everything that
raises `NotImplementedError` or shells out to a binary lives in sub-modules
(`md_runners`, `analysis`, `system_builder`) and is not in the default namespace.

| Name | Backed by | Status |
|------|-----------|--------|
| `System`, `Topology`, `Configuration` | Rust | ✅ working |
| `InputParameters` + factories | Rust | ✅ working |
| `Simulation`, `AlgorithmSequence` | Rust | ✅ working |
| `Vec3`, `Energy`, `Frame`, `rmsd`, `rdf` | Rust | ✅ working |
| `md_runners.*` (MDSimulation, run_standard_md, …) | Python subprocess → `md` binary | ⚠️ legacy, requires binary in PATH |
| `analysis.*` | Python subprocess → gromos++ programs | 🔜 programs not yet ported |
| `system_builder.*` (ForceField, molecule, …) | Python stub | 🔜 design sketch, raises NotImplementedError |

## Source layout

```
py-gromos/
├── python/gromos/
│   ├── __init__.py        ← re-exports working names; the module contract
│   ├── gromos.pyi         ← type stubs for the Rust extension
│   ├── md_runners.py      ← legacy subprocess wrappers (deprecated path)
│   ├── analysis.py        ← future analysis wrappers (mostly stub)
│   └── system_builder.py  ← future system-builder design sketch
├── tests/
│   ├── test_basic.py                  ← unit tests for Vec3, Energy, Frame, rmsd, rdf
│   └── test_gromosXX_references.py    ← 21-system reference validation suite
└── docs/                  ← this documentation
```

The Rust side lives in `crates/pyo3-gromos/src/`:

| File | Exposes |
|------|---------|
| `system.rs` | `System` |
| `topology.rs` | `Topology` |
| `py_conf.rs` | `Configuration` |
| `parameters.rs` | `InputParameters` + factories |
| `simulation.rs` | `Simulation` |
| `algorithm_sequence.rs` | `AlgorithmSequence` |
| `lib.rs` | `Vec3`, `Energy`, `Frame`, `rmsd`, `rdf`, module assembly |

## Running the test suite

```bash
make test-python           # build + full suite (89 pass, 11 skip, 0 fail)

# Just the reference energy/force/position tests
.venv/bin/pytest py-gromos/tests/test_gromosXX_references.py -v

# Just the P3.2 workflow tests
.venv/bin/pytest py-gromos/tests/test_gromosXX_references.py \
    -k "system_constructor or factory_workflow or factory_nvt" -v
```

The reference suite validates against double-precision gromosXX output for
21 systems: vacuum pairs, single molecules, small solvated boxes, bulk water
(NVE/NVT/NPT), and solvated alanine dipeptide.

## Roadmap

### P3.3 — Energy reporters ✓ done

`sim.run(steps, ene_freq)` streams energies into a NumPy array without writing
a `.tre` file; `EnergyTimeseries` (`gromos.timeseries`) wraps it with
`block_average()`, `to_dataframe()` (polars/pandas/dict), and `plot()`
(plotly/matplotlib). The `angle`/`dihedral`/`improper` energy components are
wired end to end — including a fix in the underlying Rust bonded-force
combiner (`gromos-forces/src/bonded/mod.rs`), which lumped all four bonded
terms into `bond_total` even outside the Python bindings.

### P3.4 — Working notebooks ✓ done

`notebooks/01_load_and_inspect.ipynb` and `02_short_md.ipynb` use the real
`from_files → nvt → run` path against reference systems, including the
`Topology.solvate()` + `System(topo, conf)` path for systems `from_files()`
can't load directly (unsolvated topology vs. solvated configuration).

### P3.5 — Constraints and energy minimization ✓ done

Two gaps found while writing the P3.4 notebooks, both closed:

- **Constraints on the factories.** `InputParameters.nve/nvt/npt` take
  `constraints="none"|"hbonds"|"allbonds"` (default `"none"`, matching GROMOS's
  own `NTC=1` default). Previously factory-built params had no way to turn SHAKE
  on, so a constrained system (e.g. one with solute H-bonds) would silently run
  unconstrained and diverge. `02_short_md.ipynb` §2a demonstrates the failure
  mode and the fix side by side, live.
- **Energy minimization via `Simulation`.** `InputParameters.steepest_descent(steps)`
  now actually minimizes when run through `Simulation` — previously it silently
  fell back to plain leap-frog with `dt=0` (a no-op). Also available as
  `AlgorithmSequence.minimize(topo, params)` / the `SteepestDescent` building
  block for custom sequences.

### FUTURE — System builder algebra

The full design is in `FUTURE.md`. Key idea: replace the traditional 8-binary
gromos++ pipeline with a composable Python API:

```python
ff     = ForceField.load("54A7")
system = molecule("ALA", ff) * 10 + solvent("SPC", n=2000)
system.neutralize(ion="CL")
```

Open decisions (tracked in `FUTURE.md` and `system_builder.py` comments) include:
native topology building vs. subprocess `make_top`, the `.solvate()` / `.pack()`
coordinate pipeline, and the `+` / `*` assembly algebra.

### FUTURE — SoA memory layout

The Rust core stores atoms as `Vec<Vec3>` (array-of-structs). A migration to
structure-of-arrays (`struct Soa { x: Vec<f64>, y: Vec<f64>, z: Vec<f64> }`)
would allow true zero-copy NumPy sharing and cleaner SIMD in the nonbonded loop.
Tracked in `FUTURE.md` §Dimension 1.

## References

- GROMOS force fields and theory: [www.gromos.net](https://www.gromos.net)
- PyO3 (Rust-Python bindings): [pyo3.rs](https://pyo3.rs)
- Maturin (build tool): [github.com/PyO3/maturin](https://github.com/PyO3/maturin)
