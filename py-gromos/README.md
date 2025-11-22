# py-gromos

Python bindings for GROMOS-RS molecular dynamics simulations.

⚠️ **Early Development** - This is an educational/research project in early alpha stage.

## Features

- Zero-copy NumPy integration for molecular coordinates
- SIMD-accelerated vector mathematics
- Pythonic API for MD simulations
- Advanced sampling methods (EDS, GaMD, REMD)
- Trajectory analysis tools
- Integration with GROMOS force field
- Built on PyO3 for performance

## Architecture

Inspired by [Polars](https://github.com/pola-rs/polars):

```
Python API (py-gromos)
         ↓
    PyO3 bindings
         ↓
   Rust Core (gromos-rs)
```

## Installation

### Prerequisites

- Python 3.8+
- Rust 1.70+ (install from https://rustup.rs)

### Build

```bash
pip install maturin
maturin develop --release
```

## Quick Start

### Basic Vector Operations

```python
import gromos
import numpy as np

# Create vectors
v1 = gromos.Vec3(1.0, 0.0, 0.0)
v2 = gromos.Vec3(0.0, 1.0, 0.0)

# Vector operations
v3 = v1 + v2
dot_product = v1.dot(v2)
cross_product = v1.cross(v2)
length = v1.length()

print(f"v1 + v2 = {v3}")
print(f"v1 · v2 = {dot_product}")
print(f"v1 × v2 = {cross_product}")
print(f"|v1| = {length}")
```

### Running MD Simulations

```python
from gromos import run_standard_md

# Run a standard MD simulation
run_standard_md(
    topology="system.top",
    coordinates="initial.cnf",
    parameters="md.imd",
    output_prefix="production"
)
```

### Working with Trajectories

```python
import gromos

# Load trajectory
traj = gromos.load_trajectory("trajectory.trc")

# Compute RMSD
rmsd = traj.rmsd(reference_frame=0)

# Analyze geometry
bonds = traj.compute_bond_lengths(atom_pairs=[(0, 1), (1, 2)])
angles = traj.compute_angles(atom_triples=[(0, 1, 2)])
```

## Examples

Comprehensive examples are available in the `examples/` directory:

- `01_basic_vectors.py` - Vector math and operations
- `02_system_setup.py` - Creating molecular systems
- `03_integrators.py` - Time integration methods
- `06_run_standard_md.py` - Running MD simulations
- `07_run_gamd.py` - Gaussian accelerated MD
- `08_run_eds.py` - Enveloping distribution sampling
- `09_run_remd.py` - Replica exchange MD
- `11_trajectory_analysis.py` - Analyzing trajectories
- `16_free_energy_perturbation.py` - FEP calculations

Run examples:

```bash
python examples/01_basic_vectors.py
```

## Documentation

Full documentation is available in the `docs/` directory:

- **[User Guide](docs/user-guide/)**: Installation, quick start, tutorials
- **[API Reference](docs/api/)**: Complete API documentation
- **[Development](docs/development/)**: Contributing and building guides

Build documentation locally:

```bash
pip install mkdocs mkdocs-material pymdown-extensions
mkdocs serve
```

Then open http://localhost:8000 in your browser.

## Testing

Run the test suite:

```bash
# Run Python tests
pytest tests/

# Run Rust tests
cargo test
```

## License

GPL-2.0 - See LICENSE file for details

## References

- **GROMOS**: https://www.gromos.net
- **gromos-rs**: Core Rust implementation (see `../gromos-rs`)
- **PyO3**: https://pyo3.rs
- **Maturin**: https://github.com/PyO3/maturin
