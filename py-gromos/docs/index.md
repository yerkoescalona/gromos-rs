# py-gromos Documentation

Python bindings for GROMOS-RS molecular dynamics simulations.

!!! warning "Educational Project - Not Production Ready"
    **py-gromos is an educational Python binding for learning purposes only.**

    - ⚠️ NOT intended for actual research or scientific work
    - 🚧 Cannot currently build due to gromos-rs compilation errors
    - 📚 Learning resource for Rust-Python integration with PyO3
    - 🔬 Requires gromos-rs to be functional first

    **For official GROMOS software**, visit: [www.gromos.net](https://www.gromos.net)

## What is py-gromos?

py-gromos is a Python binding for the GROMOS-RS library, providing:

- **High-performance**: Rust core with zero-copy NumPy integration
- **Safe**: Memory safety guaranteed by Rust's type system
- **Modern**: PyO3 bindings for seamless Rust-Python interoperability
- **Educational**: Learn how to build scientific Python packages backed by Rust

## Architecture

Inspired by [Polars](https://github.com/pola-rs/polars)' design:

```
┌─────────────────────────────────────────┐
│           Python API (py-gromos)        │
│  NumPy arrays, Pythonic interface       │
└──────────────┬──────────────────────────┘
               │ PyO3 bindings (zero-copy)
┌──────────────┴──────────────────────────┐
│         Rust Core (gromos-rs)           │
│  SIMD, parallelism, memory safety       │
└─────────────────────────────────────────┘
```

## Current Status

!!! danger "Not Functional"
    py-gromos **cannot be built** because gromos-rs has ~140+ compilation errors. This is purely educational scaffolding at this stage.

### What Exists
- ✅ PyO3 binding structure (src/lib.rs)
- ✅ API design (Math types, State, Integrators, etc.)
- ✅ Example code structure
- ✅ Testing framework setup
- ✅ CI/CD configuration

### What's Missing
- ❌ gromos-rs doesn't compile
- ❌ Bindings can't be built
- ❌ No actual functionality
- ❌ No performance validation
- ❌ No PyPI package

## Getting Started

### Prerequisites

```bash
# Install Rust toolchain (stable)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install Python 3.12+
python3 --version
```

### Installation (When Working)

!!! note "Not Currently Possible"
    These instructions are aspirational. The package cannot currently be built.

```bash
# Clone repository
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromosXX/py-gromos

# Install maturin (Rust-Python build tool)
pip install maturin

# Build and install in development mode
maturin develop --release
```

## Documentation Structure

### For Users
- **[Installation](user-guide/installation.md)**: Setup instructions
- **[Quick Start](user-guide/quick-start.md)**: First steps with py-gromos
- **[Learning Guide](user-guide/learning-guide.md)**: PyO3 and Rust-Python integration
- **[API Reference](api/reference.md)**: Complete API documentation

### For Developers
- **[Contributing](development/contributing.md)**: How to contribute
- **[Building](development/building.md)**: Build system and development setup
- **[Python Bindings Tutorial](development/python-bindings-tutorial.md)**: Deep dive into PyO3 bindings architecture

## Learning Goals

Use this project to learn:

### Rust-Python Integration
- PyO3 for Rust-Python bindings
- Zero-copy data sharing with NumPy
- Building Python packages with Maturin
- Type conversions and error handling

### Scientific Python
- NumPy integration patterns
- High-performance Python packages
- API design for scientific computing
- Testing strategies (pytest + Rust tests)

### Software Engineering
- Multi-language project structure
- CI/CD for Rust-Python packages
- Documentation with MkDocs
- Packaging and distribution

## Comparison with Polars

py-gromos is heavily inspired by Polars' architecture:

| Aspect | Polars | py-gromos |
|--------|--------|-----------|
| Core language | Rust | Rust (gromos-rs) |
| Python bindings | PyO3 | PyO3 |
| Data sharing | Apache Arrow (zero-copy) | NumPy (zero-copy) |
| Parallelism | Rayon | Rayon |
| SIMD | Yes | Yes (planned) |
| Domain | DataFrames | Molecular dynamics |
| Status | Production | Educational |

## License

py-gromos is licensed under GPL-2.0, maintaining compatibility with GROMOS.

- Copyright (C) 2025 Yerko Escalona
- Based on GROMOS software (Copyright by Biomos b.v. and GROMOS development team)
- See [LICENSE](../LICENSE) file for full text

## References

- **GROMOS**: [www.gromos.net](https://www.gromos.net)
- **Polars**: [github.com/pola-rs/polars](https://github.com/pola-rs/polars)
- **PyO3**: [pyo3.rs](https://pyo3.rs)
- **Maturin**: [github.com/PyO3/maturin](https://github.com/PyO3/maturin)

## Citation

!!! warning "Not for Research Use"
    This is an educational project and should NOT be cited in scientific publications.
    Use official GROMOS software for research work.

For educational references, cite the original GROMOS papers and acknowledge the inspiration from Polars.
