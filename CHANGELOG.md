# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [0.0.10] (2026-05-17)

### Refactor

- **workspace:** restructure crates into focused responsibilities
  - `gromos-md`: 8 simulation engine binaries (md, md_mpi, md_mpi_cuda, mdf, remd, repex_mpi, eds, gamd) + integration tests
  - `gromos-tools`: 30 system construction binaries organized in subdirectories (topology/, box/, conversion/, utilities/)
  - `gromos-analysis`: 66 analysis binaries organized in subdirectories (structural/, energy/, distribution/, dynamics/, free_energy/, trajectory/, noe/, clustering/, xray/, special/) + existing library code
  - `gromos-cli`: slimmed to thin unified `gromos` multicall binary (clap only)
  - `gromos` facade: removed `gromos-analysis` re-export to avoid circular dependency; analysis bins use `gromos_core`/`gromos_io` directly
  - All 21 reference tests pass, workspace compiles clean

## [0.0.9](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.8...v0.0.9) (2026-03-29)

### Bug Fixes

- **forces:** fix CRF energy mismatch (~3.82 kJ/mol) for systems with solvent near cutoff boundary
  - Expand solvent long-range pairlist to all atom pairs (was: first-atom only with shared PBC shift)
  - Remove HEAVISIDE truncation from `lj_crf_innerloop` to match gromosXX `#undef XXHEAVISIDE` default
  - Long-range solvent now processed with `lj_crf_innerloop` (per-atom nearest_image) instead of `solvent_innerloop`

### Features

- **tests:** add water_10_box and nacl_water_box_shifted reference systems
  - water_10_box: 2 ions + 10 SPC waters (32 atoms), positions away from cutoff boundaries
  - nacl_water_box_shifted: nacl_water_box with perturbed positions near cutoff
- **tests:** promote nacl_water_box, nacl_water_box_shifted, nacl_3water_cutoff from ignored to active
  - 14 of 19 reference tests now pass (was 11)
- **io:** add verbose logging and md_output.log capture to run_references.py
- **io:** improve force trajectory writer robustness

## [0.0.8](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.7...v0.0.8) (2026-03-08)

### Features

- **core:** add solute/solvent boundary fields to Topology
  - Add `num_solute_chargegroups` field and `num_solute_atoms()` method
  - Set during `build_topology()` for solute/solvent dispatch

- **forces:** gromosXX-compatible nonbonded force architecture
  - Add HEAVISIDE truncation to `lj_crf_innerloop` (skip pairs beyond cutoff²)
  - Add `solvent_innerloop` with shared PBC shift from O-O nearest_image
  - Split `rf_excluded_interactions` into solute/solvent paths:
    - Solute: self-term + excluded pair forces + full RF energy
    - Solvent: no self-term, no forces, only distance-dependent energy
  - Add `cutoff_sq` to CRFParameters

### Refactor

- **pairlist:** separate solute/solvent CG pairs with correct distance metrics
  - Solute CGs: center-of-geometry distance, exclusion checks
  - Solvent CGs: first-atom position distance (gromosXX convention)
  - Intra-CG non-excluded pairs for solute CGs
  - Solvent-solvent stores first-atom pairs only (expanded in innerloop)
  - Add debug logging for CG positions, distances, pair classification

## [0.0.7](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.6...v0.0.7) (2026-03-02)

### Refactor

- **io:** gromosXX-compatible topology init order with solvate() ([60207d9](https://github.com/yerkoescalona/gromos-rs/commit/60207d9584e67c0ae622ed6267f19fd38965467e))
  - Follow gromosXX initialization: read_parameter → read_topology → solvate → read_configuration
  - Add Topology::solvate(nsm) method matching gromosXX topo.solvate(0, nsm)
  - Add SolventAtomTemplate/SolventConstraintTemplate types on Topology
  - build_topology() now stores solvent template without expanding
  - Remove build_topology_with_solvent() — NSM comes from IMD SYSTEM block
  - Reorder md.rs: IMD → topology+solvate → coordinates (was: IMD → coords → topo)
  - Store chargegroup_codes on Topology for later use
  - LJ matrix now includes solvent IAC types at build time
  - All 6 reference tests pass (pair_lj, pair_lj_mixed, nacl_pair, water_single, benzene_vacuum, water_3_box)

## [0.0.6](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.5...v0.0.6) (2026-02-27)

### Features

- **forces:** gromosXX-compatible bonded forces, topology parsing, and SHAKE ([5d409e7](https://github.com/yerkoescalona/gromos-rs/commit/5d409e7da1ffb9528db09d122e68515bda4df6bc))
  - Fix bonded force vector conventions to match gromosXX (v = pos(i) - pos(j))
  - Simplify quartic bond force: avoid unnecessary r division, use r^2 directly
  - Add NTF flag support to calculate_bonded_forces_ntf() for selective force terms
  - Add debug logging for individual bonded energy components
  - RF excluded corrections: add forces (not just energy), rename to rf_excluded_interactions
  - Parse dihedral/improper types and terms from topology (TORSDIHEDRALTYPE, DIHEDRAL, IMPDIHEDRAL)
  - Fix SOLUTEATOM parser: exclusions on same line (col 8+), skip INE14 line
  - Rewrite SHAKE to match gromosXX convergence (tolerance on |constraint_length - r|)
  - Add shake_algorithm.rs module
  - Add ForceStorage to nonbonded for combined force+energy tracking

## [0.0.5](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.4...v0.0.5) (2026-02-27)

### Refactor

- **io:** gromosXX-compatible CLI and IO cleanup ([c0e9ea2](https://github.com/yerkoescalona/gromos-rs/commit/c0e9ea2603942bed73944dab073cd58045966f5d))
  - Move coordinate reader from md.rs to gromos-io (read_coordinates + CoordinateData)
  - Add GENBOX, POSITIONRED, VELOCITYRED block support in coordinate reader
  - Rewrite IMD parser for gromosXX positional format (comment headers + data lines)
  - Add NTF force flags, COMTRANSROT, shake_tol to ImdParameters
  - Refactor md.rs CLI to match gromosXX @-parameters: @topo @conf @input (required), @fin @trc @tre @trf @trv (output), @pttopo @posresspec @refpos @distrest @angrest @dihrest @gamd (input), @verb @print @version @develop (control)
  - All simulation parameters now come from @input .imd/.in file
  - Fix pre-existing f64/usize type errors in energy_binary and trajectory_binary tests
  - 9 new/updated tests for coordinate and IMD parsing
  - pair_lj reference still passes at ~1e-11 precision

## [0.0.4](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.3...v0.0.4) (2026-02-24)

### Refactor

- **core:** implement gromosXX Algorithm/AlgorithmSequence pattern ([c61cd57](https://github.com/yerkoescalona/gromos-rs/commit/c61cd57c5b7b1eb29c0f067fd02038560b923ab9))
  - Add Algorithm trait and AlgorithmSequence in gromos-core with 7 unit tests
  - Implement algorithm wrappers: Forcefield, LeapFrogVelocity, LeapFrogPosition, TemperatureCalculation, EnergyCalculation
  - TemperatureCalculation uses gromosXX averaged kinetic energy formula
  - Refactor md.rs to use AlgorithmSequence instead of manual force+integrate loop
  - Energy convention follows gromosXX: write to current(), exchange_state moves to old(), read from old() for output
  - Energy output uses full f64 scientific notation for exact comparison
  - All 10 steps of pair_lj match gromosXX reference to ~1e-11 precision
  - Fix minor clippy warnings in gamd, remd, integrator, pairlist

## [0.0.3](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.2...v0.0.3) (2026-02-22)

### Refactor

- **forces:** enhance nonbonded interactions with RF corrections and update math structures ([51b8681](https://github.com/yerkoescalona/gromos-rs/commit/51b86819e74ffcbe918b30b00f180b4c5be2c9ae))
- **core:** support f64 instead of f32 ([ade5423](https://github.com/yerkoescalona/gromos-rs/commit/ade5423ca49fc8974dfdc624f71770a4b7c10a92))

### Chores

- remove dead logging module and replace black with ruff ([de003df](https://github.com/yerkoescalona/gromos-rs/commit/de003dfef3066db053255dbc58087739becf87f7))

## [0.0.2](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.1...v0.0.2) (2025-12-21)

### Refactor

- **workspace:** reorganize into crate-based workspace structure ([7d97a54](https://github.com/yerkoescalona/gromos-rs/commit/7d97a54096d902529e878374395b837717e9343f))

## [0.0.1](https://github.com/yerkoescalona/gromos-rs/commits/v0.0.1) (2025-11-22)

### Refactor

- update force writer integration tests for new API ([ec11525](https://github.com/yerkoescalona/gromos-rs/commit/ec11525f4e08aabc8d621a187cd8875afd3b90f6))

### Chores

- clean and update README files ([27f80bc](https://github.com/yerkoescalona/gromos-rs/commit/27f80bc3abcad0df3d19bd20d5a5291ce7f5fb8e))
- add GROMOS and claude entries to .gitignore ([926ff5c](https://github.com/yerkoescalona/gromos-rs/commit/926ff5c41d8658b964a09aa9d93682f7627ba264))

### Initial

- initial commit ([c0a551c](https://github.com/yerkoescalona/gromos-rs/commit/c0a551cfcac12b0a60daaad908d151f689395732))
