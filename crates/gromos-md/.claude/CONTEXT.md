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
- **37/40 reference tests passing** (3 ignored: `aladip_vacuum_fep` known FEP mismatch,
  `aladip_vacuum_em` EM frame-count off-by-one, `water_216_nve_nobath` — a *correct* reference the
  engine fails because the IMD parser turns an absent MULTIBATH block into a Berendsen bath, PLAN.md
  3.9 A18 — see PLAN.md Reference Test Status for the full matrix)
- `PairlistAlgorithm::from_imd(imd.algorithm, …)` dispatches Standard/CellList here (PLAN.md 9a-0/
  9a-1); CellList only activates on explicit `ALGORITHM grid_cell`, never via a size heuristic — all
  37 active reference systems still use `standard`.
- `test_provider_reference.rs` / `test_orchestrator_reference.rs` — gromosXX-reference tests for the
  `PotentialProvider`/`ProviderOrchestrator` seam (PLAN.md 2.6/P2.8-2), not the main `md` binary path.
  `schnet_nve_loop.rs` (feature `ml`) — real leapfrog loop driving `SchNetInteraction`;
  `xtb_orchestrator_sequence.rs` — `ProviderOrchestratorAlgorithm` inside a real `AlgorithmSequence`
  (PLAN.md P2.8-6 ✓). Neither is wired into the `md` CLI yet.
- **PLAN.md 3.9 (next): `md.rs`'s run assembly moves out of this crate.** `src/lib.rs` exports
  nothing, which is why `pyo3-gromos` re-implements `md.rs`'s IMD→`AlgorithmSequence` assembly (and
  has drifted: `four_pi_eps_i`, `parallel_nonbonded`, NSM-from-coordinates, DOF, FEP). After a
  five-framework review the target is a new library crate `crates/gromos-run` (no clap/env_logger/
  mpi/cudarc; serde non-optional) with `RunRecipe` (+ `from_imd`/`to_imd`, `Diagnostics`),
  `prepare_system`, `build_plan` → `Vec<AlgorithmSpec>`, `validate_plan`, `instantiate` (reads only
  the plan), `start` (init + step 0); this crate's binaries and `pyo3-gromos` both call it. Do not add
  a second sequence builder or a second IMD reader anywhere. Known parser hazard recorded there: an
  `.imd` without MULTIBATH silently runs a Berendsen bath (τ = 0.1) — fixed by 3.9 step 2.

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
