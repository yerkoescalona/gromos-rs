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

- **Ports of gromos++ programs with reference tests (2026-08-30, `tests/gromospp_references.rs`,
  fixtures + exact gromos++ commands in `tests/data/README.md`):** `tser` (rewritten: property
  specifiers, distributions), `mdf`, `dg_ener`, `dfmult`, `matrix_overlap`. In gromos-tools:
  `explode`, `duplicate`, `pt_top`.

## Shared machinery (use it; do not re-implement per binary)
- `gromos_io::args::Arguments` — gromos++ `Arguments` (`@key` multi-valued, `count`, `@f` files).
- `gromos_io::pbc::Pbc` — `@pbc <r|v> [cog|nog]`, `periodicity(box)`, `gather()` = gromos++ `coggather`.
- `gromos_analysis::time::Time` — `@time t0 dt` or the frame's time; `Time::fmt` = gromos++ layout.
- `gromos_analysis::property` — `d/a/t/tp/o/op/pr/pa` specifiers, `PropertyContainer`,
  `atoms_to_string` (gromos++ `AtomSpecifier::toString`), `ordered_atoms` (order-preserving).
- `gromos_analysis::distribution` — gromos++ `Distribution`; `cpp_g` = C++ default float format.
- `gromos_analysis::lnexp` — log-exponential averages / covariance / statistical inefficiency.
- `gromos_io::table::read_columns` — numeric column files; `gromos_io::coordinate::format_g96`.
- `gromos_core::Stat::ee_strict` — gromos++'s error estimate (`nan` below 200 values).
- gromos-core `nearest_image(ri, rj)` returns the minimum-image **vector** ri − rj; gromos++'s
  `nearestImage(ref, pos)` is the image **position** `ref − that vector`. Every port must mind this.

## Status — stubs / parked
- `visco`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy` — stubs; parked
- Ported with reference tests (batch 2): `bilayer_dist`, `bilayer_oparam`, `jval`, `edyn`
  (eigenvector projections — gromos++'s are a bug, see CHANGELOG 0.0.40), `gca`.
- Ported with reference tests (batch 3): `eds_update_1`, `eds_update_2` (`src/eds.rs` shared),
  `jepot` (`gromos_io::jvalue` readers, `time::TimeSpec` shared with `jval`), `pocket` (no
  gathering — gromos++ never applies `@pbc` there), `dfgrid` (gromos++'s open-set Dijkstra kept
  for identical tie-breaking; `cnf` output, real atoms only); `make_sasa_top` lives in
  gromos-tools.
- Not ported, with reasons (CHANGELOG 0.0.41): `cos_dipole`, `cos_epsilon` (COS special
  trajectories), `dGslv_pbsolv` (PB solver), `prep_bb` (interactive), `prep_xray_le` (X-ray).

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
