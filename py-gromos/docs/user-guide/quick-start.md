# Quick Start Guide

Python bindings for GROMOS-RS molecular dynamics simulations.

!!! warning "Educational Project - Not Functional"
    **py-gromos cannot currently be built** because gromos-rs has ~140+ compilation errors. This guide shows the intended API design for educational purposes.

## Architecture

Inspired by [Polars](https://github.com/pola-rs/polars)' design:

- **Rust core**: Rust MD engine (when functional)
- **PyO3 bindings**: Safe, zero-copy data sharing between Rust and Python
- **NumPy integration**: Direct array conversion without copying
- **Memory safety**: Guaranteed by Rust's type system

## Installation

!!! danger "Not Currently Possible"
    py-gromos cannot be built because gromos-rs doesn't compile. These instructions are aspirational.

### From source (when working)

```bash
# Install maturin (Rust-Python build tool)
pip install maturin

# Build and install in development mode
cd py-gromos
maturin develop --release

# Or build wheel for distribution
maturin build --release
pip install target/wheels/gromos-*.whl
```

### From PyPI (future)

```bash
pip install gromos  # Not available yet
```

## Quick Start Example

!!! note "API Design Only"
    This shows the intended API. The code cannot actually run yet.

```python
import gromos
import numpy as np

# Create system state
state = gromos.State(
    num_atoms=1000,
    num_temp_groups=1,
    num_energy_groups=1
)

# Set up simulation box
box = gromos.Box.rectangular(5.0, 5.0, 5.0)  # 5x5x5 nm

# Create integrator
integrator = gromos.LeapFrog(dt=0.002)  # 2 fs timestep

# Work with 3D vectors
v1 = gromos.Vec3(1.0, 2.0, 3.0)
v2 = gromos.Vec3(4.0, 5.0, 6.0)
distance = v1.distance(v2)
print(f"Distance: {distance:.3f} nm")

# Convert to/from NumPy
vec_array = v1.to_numpy()
v3 = gromos.Vec3.from_numpy(np.array([1.0, 2.0, 3.0]))

# Access simulation data as NumPy arrays (zero-copy)
positions = state.positions()  # Shape: (N, 3)
velocities = state.velocities()  # Shape: (N, 3)
forces = state.forces()  # Shape: (N, 3)
```

## Advanced Sampling

### Gaussian Accelerated MD (GaMD)

```python
# Create GaMD parameters
gamd_params = gromos.GamdParameters(
    sigma0=6.0,
    threshold_mode='lower'
)

# Initialize GaMD runner
gamd = gromos.GamdRunner(gamd_params)
```

### Enveloping Distribution Sampling (EDS)

```python
# Create EDS parameters for 4 states
eds_params = gromos.EDSParameters(
    num_states=4,
    smoothness=1.0
)

# Initialize EDS runner
eds = gromos.EDSRunner(eds_params)
```

### Replica Exchange MD (REMD)

```python
# Create REMD controller
remd = gromos.ReplicaController(
    num_replicas=8,
    exchange_interval=1000  # Attempt exchange every 1000 steps
)

print(f"Managing {remd.num_replicas()} replicas")
```

## Features

### Math Types
- `Vec3`: SIMD-accelerated 3D vectors
- `Mat3`: SIMD-accelerated 3×3 matrices

### Core Structures
- `State`: System state (positions, velocities, forces)
- `Energy`: Energy tracking (kinetic, potential, components)
- `Configuration`: Complete system configuration
- `Topology`: Molecular topology (atoms, bonds, parameters)
- `Box`: Simulation box (vacuum, rectangular, triclinic)

### Integrators
- `LeapFrog`: Fast velocity Verlet variant
- `VelocityVerlet`: Higher accuracy integrator
- `StochasticDynamics`: Langevin dynamics for implicit solvent

### Advanced Sampling
- `GamdParameters`, `GamdRunner`: Gaussian Accelerated MD
- `EDSParameters`, `EDSRunner`: Enveloping Distribution Sampling
- `ReplicaController`: Replica Exchange MD

## Design Goals

py-gromos aims to provide (when functional):

| Feature | Design Intention |
|---------|------------------|
| SIMD vectorization | Automatic via glam + wide |
| Parallel execution | Rayon multi-threading |
| Zero-copy arrays | NumPy integration via PyO3 |
| Memory safety | Rust guarantees, no segfaults |
| Zero-cost abstractions | No runtime overhead |

!!! warning "Performance Not Validated"
    Performance claims like "2-3x speedup" have NOT been measured. gromos-rs doesn't compile yet, so no benchmarks exist.

### Why Zero-Copy Matters

Like Polars, GROMOS-RS uses Apache Arrow-compatible memory layout:

```python
# No data copying - direct memory access
positions = state.positions()  # Returns view of Rust memory
positions[0] = [1.0, 2.0, 3.0]  # Would need copy for mutation

# Setting data requires a copy
new_positions = np.random.rand(1000, 3).astype(np.float32)
state.set_positions(new_positions)  # Copies into Rust
```

## Comparison with Polars

| Aspect | Polars | GROMOS-RS |
|--------|--------|-----------|
| Core language | Rust | Rust |
| Python bindings | PyO3 | PyO3 |
| Data sharing | Arrow zero-copy | NumPy zero-copy |
| Parallelism | Rayon | Rayon |
| SIMD | Yes | Yes (glam) |
| Domain | DataFrames | Molecular dynamics |

## Development

### Building

```bash
# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install maturin
pip install maturin

# Build in development mode (includes debug symbols)
maturin develop

# Build optimized release
maturin develop --release
```

### Testing

```bash
# Run Python tests
pytest tests/

# Run Rust tests
cargo test -p gromos-rs
```

### Project Structure

```
py-gromos/
├── Cargo.toml           # Rust package manifest
├── pyproject.toml       # Python package config
├── src/
│   └── lib.rs          # Rust ↔ Python bindings (PyO3)
├── python/
│   └── gromos/
│       └── __init__.py  # Python API
├── tests/              # Python tests
└── examples/           # Example scripts
```

## Next Steps

- **[API Reference](../api/reference.md)**: Complete API documentation
- **[Learning Guide](learning-guide.md)**: Learn PyO3 and Rust-Python integration
- **[Contributing](../development/contributing.md)**: Help make this functional

## License

GPL-2.0 - Same as GROMOS

## References

- **GROMOS**: [www.gromos.net](https://www.gromos.net)
- **Polars**: [github.com/pola-rs/polars](https://github.com/pola-rs/polars)
- **PyO3**: [pyo3.rs](https://pyo3.rs)
- **NumPy**: [numpy.org](https://numpy.org)

## Citation

!!! warning "Not for Research Use"
    This is an educational project. Use official GROMOS software for research.
