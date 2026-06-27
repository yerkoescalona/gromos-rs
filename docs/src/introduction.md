# gromos-rs

A GROMOS-faithful molecular dynamics engine written in Rust.

gromos-rs is a ground-up reimplementation of the [GROMOS](http://www.gromos.net) MD engine
targeting bit-for-bit numerical agreement with gromosXX on all validated reference systems,
while diverging deliberately on architecture: single non-forked inner loop, owned SoA data
layout, O(N) pairlist, and a composable Python API.

## Design principles

- **Reference fidelity first.** Every physics feature lands with a GROMOS reference test.
  The expected outputs in `crates/gromos-md/tests/gromosXX_references/*/expected/` are ground
  truth and are never modified.

- **One implementation, not four.** gromosXX ships four pairlist algorithms, three constraint
  solvers, and a forked OpenMP/MPI/CUDA path for every tight loop. gromos-rs has one of each,
  chosen to be correct and fast, with parallelism added transparently via Rayon.

- **Rust's type system enforces physics invariants.** Boundary conditions, molecule roles, and
  pairlist categories are encoded as types, not integers — misuse is a compile error, not a
  silent wrong answer.

## API docs

The full API reference is generated from source:

```sh
cargo doc --workspace --no-deps --open
```

Or browse it at [docs.rs/gromos-rs](https://docs.rs/gromos-rs) once published.

## Crate map

| Crate | Role |
|-------|------|
| `gromos-core` | Math types, topology, configuration, pairlist |
| `gromos-forces` | Bonded and nonbonded force kernels |
| `gromos-integrators` | MD integrators, thermostat, barostat, constraints |
| `gromos-io` | GROMOS file format readers/writers |
| `gromos-analysis` | Analysis tools (RMSD, order parameters, FE estimators) |
| `gromos-md` | The `md` binary and reference test suite |
