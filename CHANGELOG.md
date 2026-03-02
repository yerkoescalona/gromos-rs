# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [0.0.6](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.5...v0.0.6) (2026-02-22)

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

## [0.0.5](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.4...v0.0.5) (2026-02-22)

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

## [0.0.4](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.3...v0.0.4) (2026-02-22)

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
