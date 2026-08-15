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
- **37/39 reference tests passing** (2 ignored: `aladip_vacuum_fep` known FEP mismatch,
  `aladip_vacuum_em` EM frame-count off-by-one — see PLAN.md Reference Test Status for the full matrix)
- `PairlistAlgorithm::from_imd(imd.algorithm, …)` dispatches Standard/CellList here (PLAN.md 9a-0/
  9a-1); CellList only activates on explicit `ALGORITHM grid_cell`, never via a size heuristic — all
  37 active reference systems still use `standard`.
- `test_provider_reference.rs` / `test_orchestrator_reference.rs` — gromosXX-reference tests for the
  `PotentialProvider`/`ProviderOrchestrator` seam (PLAN.md 2.6/P2.8-2), not the main `md` binary path.
  `schnet_nve_loop.rs` (feature `ml`) — real leapfrog loop driving `SchNetInteraction`, also not
  wired into `md` itself (PLAN.md P2.8-6, still open).

## Key files
```
src/bin/md.rs                               — main MD driver, CLI, simulation setup
tests/test_gromosXX_references.rs           — integration tests vs gromosXX
tests/run_references.py                     — generate gromosXX reference data
tests/gromosXX_references/                  — reference input + expected output
tests/test_provider_reference.rs            — LjCrfInteraction vs gromosXX (provider seam)
tests/test_orchestrator_reference.rs        — ProviderOrchestrator vs gromosXX (multi-provider seam)
tests/schnet_nve_loop.rs                    — SchNetInteraction + LeapFrog NVE loop (feature `ml`)
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
