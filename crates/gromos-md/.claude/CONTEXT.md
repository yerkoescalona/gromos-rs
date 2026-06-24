# gromos-md — stage contract

## Job
L3 orchestration. The `md` binary: composes gromos-core + gromos-forces + gromos-integrators +
gromos-io into a full MD simulation loop matching gromosXX output bit-for-bit.

## Inputs (consumes from)
All lower crates (gromos-core, gromos-forces, gromos-integrators, gromos-io).

## Outputs
`md` binary: runs GROMOS MD simulations; writes trajectories, energies, force trajectories
matching gromosXX to reference tolerances.

## Status
- Full MD loop ✓: NVE / NVT (Berendsen + NHC) / NPT
- Force trajectory output (@trf) ✓: FREEFORCERED + CONSFORCERED, atom-by-atom tolerance 1e-6
- **34/34 reference tests passing** (see PLAN.md Reference Test Status for the full matrix)
- `CellListPairlistAlgorithm` exists in gromos-core but is **not yet wired here** — `StandardPairlistAlgorithm` is still the default; system-size heuristic is the remaining step (P1.6 deferred)

## Key files
```
src/bin/md.rs                               — main MD driver, CLI, simulation setup
tests/test_gromosXX_references.rs           — integration tests vs gromosXX
tests/run_references.py                     — generate gromosXX reference data
tests/gromosXX_references/                  — reference input + expected output
```

## Crate-specific rules
- **Reference test oracle = gromosXX.** A change that alters wired, tested output is a regression until proven a deliberate, documented decision.
- **DON'T edit `gromosXX_references/*/expected/`** — those are ground truth.
- **All simulation parameters from `@input` .imd/.in file** (not CLI flags).
- **CLI:** clap `#[derive(Parser)]` with `gromos_args()` pre-processor (`@key → --key`, `@f argfile` expansion). No custom arg parsers.
- **Tolerances:** force_abs=1e-6, energy_rel=1e-8, position_abs=1e-9.

## How to add a reference test
1. Add to `SYSTEMS` list in `tests/run_references.py`
2. Run `python3 tests/run_references.py --md-binary .local/gromosXX/md++/build/program/md`
3. Add `ref_test!(name, "dir")` in `tests/test_gromosXX_references.rs`
