# py-gromos

Python bindings for GROMOS-RS molecular dynamics simulations.

!!! warning "Educational Project - Not Functional"
    **py-gromos cannot currently be built** because gromos-rs has ~140+ compilation errors. This is an educational project demonstrating Rust-Python integration with PyO3.

## Quick Links

- **[Documentation](docs/index.md)**: Complete documentation
- **[Installation](docs/user-guide/installation.md)**: Setup instructions (when functional)
- **[Quick Start](docs/user-guide/quick-start.md)**: Getting started guide
- **[API Reference](docs/api/reference.md)**: Complete API documentation
- **[Python Bindings Tutorial](docs/development/python-bindings-tutorial.md)**: Learn how PyO3 bindings work
- **[Contributing](docs/development/contributing.md)**: How to contribute

## Architecture

Inspired by [Polars](https://github.com/pola-rs/polars):

```
Python API (py-gromos)
         ↓
    PyO3 bindings
         ↓
   Rust Core (gromos-rs)
```

## Current Status

py-gromos is **not functional** because:
- gromos-rs has ~140+ compilation errors
- Bindings cannot be built
- No code can run

This is an educational project for learning:
- Rust-Python integration with PyO3
- Zero-copy NumPy integration
- Building scientific Python packages backed by Rust

## Installation (Aspirational)

```bash
# When functional:
pip install maturin
maturin develop --release
```

## Example API Design

```python
import gromos
import numpy as np

# API design (not functional yet)
v = gromos.Vec3(1.0, 2.0, 3.0)
print(v.length())
```

## Documentation

Full documentation is available in the `docs/` directory:

- **User Guide**: Installation, quick start, learning resources
- **API Reference**: Complete API documentation
- **Development**: Contributing and building guides

Build docs locally:

```bash
pip install mkdocs mkdocs-material
mkdocs serve
```

## License

GPL-2.0 - Same as GROMOS

Copyright (C) 2025 Yerko Escalona

## References

- **GROMOS**: [www.gromos.net](https://www.gromos.net)
- **Polars**: [github.com/pola-rs/polars](https://github.com/pola-rs/polars)
- **PyO3**: [pyo3.rs](https://pyo3.rs)
- **Maturin**: [github.com/PyO3/maturin](https://github.com/PyO3/maturin)
