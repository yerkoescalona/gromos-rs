# gromos-analysis — stage contract

## Job
L4 facade. GROMOS-style analysis programs as a thin layer over the shared engine core.
Zero physics re-implementation — all physics comes from `gromos-forces`.

## Inputs (consumes from)
`gromos-core` + `gromos-io` + `gromos-forces` (for `single_point_energy`).

## Outputs
Analysis binaries producing results comparable to GROMOS analysis tools.

## Status — implemented
- `atominfo` ✓ — reads topology + AtomSpecifier, prints TITLE+ATOMS; verified against GROMOS atominfo
- `rmsd` ✓ — Kabsch fit (`fit.rs`), @atomspec, @ref, @pbc, @nofit
- `nhoparam` ✓ — N-H order parameters S², rotational fit, window averaging, `ee()`
- `frameout` ✓ — full GROMOS feature parity: @pbc gather, @include SOLUTE/SOLVENT/ALL,
  @ref/@atomsfit rotational fit, @spec ALL/EVERY/SPEC, @time range, cnf/pdb/trc output, @single;
  all I/O through gromos-io; 8 integration tests
- `ener` ✓ — energy recalculation from trajectory via `gromos-forces::energy::single_point_energy`
- `ext_ti_ana` ✓ — trapezoidal ΔG ± ee() from `.trg` files
- `bar` ✓ — BAR iteration (numerically stable log-sum-exp), bootstrap error
- `ext_ti_merge` ✓ — λ interpolation, trapezoidal ΔG
- `fit.rs` ✓ — Kabsch rotational fit (Horn 1987 quaternion); `superimpose()`, `rmsd()`, 7 unit tests

## Status — stubs / parked
- `visco`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy` — stubs; parked

## Key files
```
src/fit.rs                          — Kabsch fit, superimpose(), rmsd()
src/bin/structural/rmsd.rs          — rmsd binary
src/bin/structural/atominfo.rs      — atominfo binary
src/bin/trajectory/frameout.rs      — frameout binary (full GROMOS feature parity)
src/bin/noe/nhoparam.rs             — N-H order parameters
src/bin/energy/ener.rs              — energy recalculation
src/bin/free_energy/ext_ti_ana.rs   — TI integration
src/bin/free_energy/bar.rs          — BAR estimator
src/bin/free_energy/ext_ti_merge.rs — λ window interpolation
tests/test_frameout.rs              — 8 frameout integration tests
```

## Crate-specific rules
- **Zero physics re-implementation.** All force/energy calls go through `gromos-forces`.
- **All I/O through gromos-io.** No hand-rolled `BufWriter`/`File` in binary entry points;
  use `write_g96`, `write_pdb_positions`, `TrajectoryWriter::write_trc_frame`.
- **Kabsch fit lives here** (`src/fit.rs`), not in gromos-core — it is analysis-only.
- **PBC gathering lives in gromos-core** (`src/gather.rs`) — shared with MD engine.
