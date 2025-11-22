# GROMOS-RS

Rust implementation of molecular dynamics simulation engine based on the GROMOS force field.

⚠️ **Early Development** - This is an educational/research project in early alpha stage.

## Features

- Molecular topology and force field parameters
- Integrators (Leap-Frog, Velocity Verlet, Stochastic Dynamics)
- Nonbonded interactions (LJ + Coulomb with reaction field)
- Bonded interactions (bonds, angles, dihedrals)
- Enhanced sampling methods (EDS, GaMD, REMD)
- Free energy perturbation (FEP)
- Parallel computing with Rayon
- Python bindings (see `../py-gromos`)

## Structure

```
gromos-rs/
├── src/
│   ├── topology.rs         # Molecular structure
│   ├── configuration.rs    # System state
│   ├── integrator.rs       # Time integration
│   ├── interaction/        # Force calculations
│   ├── pairlist.rs         # Neighbor lists
│   ├── algorithm/          # Constraints, thermostats
│   ├── io/                 # File I/O (G96, trajectories)
│   ├── eds.rs              # Enveloping Distribution Sampling
│   ├── gamd.rs             # Gaussian Accelerated MD
│   ├── remd.rs             # Replica Exchange MD
│   └── fep.rs              # Free Energy Perturbation
└── tests/                  # Test suite
```

## Building

### Prerequisites

- Rust 1.70+ (install from https://rustup.rs)

### Build

```bash
cargo build --release
```

### Run Tests

```bash
cargo test
```

## Optional Features

Enable optional features with cargo flags:

```bash
# FFTW for PME (requires FFTW3 installed)
cargo build --features use-fftw

# MPI for parallel replicas (requires MPI installed)
cargo build --features use-mpi

# CUDA for GPU acceleration (requires CUDA toolkit)
cargo build --features use-cuda

# Better memory allocator
cargo build --features mimalloc
```

## License

GPL-2.0 - See LICENSE file for details
