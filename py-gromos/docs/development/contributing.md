# Contributing to py-gromos

Thank you for your interest in contributing to py-gromos!

!!! warning "Early Stage Project"
    py-gromos is in very early stages. gromos-rs doesn't compile yet (~140+ errors), so py-gromos cannot be built. Contributions are welcome, but understand that basic functionality doesn't exist yet.

## Current Status

### What Works
- ✅ Project structure
- ✅ PyO3 binding skeleton
- ✅ API design
- ✅ Testing framework
- ✅ CI/CD configuration

### What Doesn't Work
- ❌ gromos-rs compilation
- ❌ Building py-gromos
- ❌ Running any code
- ❌ Running tests
- ❌ Performance benchmarks

## Ways to Contribute

### 1. Fix gromos-rs Compilation Errors

The biggest blocker is that gromos-rs doesn't compile. Help fix:

- Type errors (CudaError → DriverError, etc.)
- Missing implementations
- Feature flag issues
- Dependency conflicts

See **[gromos-rs contributing guide](../../gromos-rs/docs/development/contributing.md)**

### 2. Improve Documentation

- Add examples for PyO3 patterns
- Document API design decisions
- Write tutorials for Rust-Python integration
- Improve error messages

### 3. Extend API Design

Design (but don't necessarily implement) new features:

- Additional integrators
- More sampling methods
- I/O formats
- Analysis tools

### 4. Testing Strategy

Design test strategy for when code becomes functional:

- Unit tests for bindings
- Integration tests with NumPy
- Performance benchmarks
- Memory leak detection

## Development Setup

### Prerequisites

```bash
# Install Rust (stable)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install Python 3.12+
python3 --version

# Install development tools
pip install maturin pytest ruff mypy
```

### Fork and Clone

```bash
# 1. Fork the repository on GitHub

# 2. Clone your fork
git clone https://github.com/YOUR_USERNAME/gromosXX.git
cd gromosXX/py-gromos

# 3. Add upstream remote
git remote add upstream https://github.com/yerkoescalona/gromos-rs.git

# 4. Create feature branch
git checkout -b feature/my-contribution
```

### Project Structure

```
py-gromos/
├── Cargo.toml           # Rust package manifest
├── pyproject.toml       # Python package config
├── Makefile             # Build commands
├── src/
│   └── lib.rs          # PyO3 bindings (Rust → Python)
├── python/
│   └── gromos/
│       └── __init__.py  # Python API
├── tests/              # Python tests (pytest)
├── examples/           # Example scripts
└── docs/               # MkDocs documentation
```

## Code Style

### Rust Code (src/lib.rs)

Follow gromos-rs style guide:

```rust
use pyo3::prelude::*;

/// Create Vec3 from Python
#[pyclass]
pub struct Vec3 {
    #[pyo3(get, set)]
    pub x: f64,
    #[pyo3(get, set)]
    pub y: f64,
    #[pyo3(get, set)]
    pub z: f64,
}

#[pymethods]
impl Vec3 {
    #[new]
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
}
```

**Style rules:**
- Use `#[pyclass]` for types exposed to Python
- Use `#[pymethods]` for methods
- Document with rustdoc comments (`///`)
- Follow Rust naming conventions

### Python Code (python/gromos/)

```python
"""Python API for GROMOS-RS."""

from typing import List, Optional
import numpy as np
import numpy.typing as npt

from . import _gromos  # Rust extension

def example_function(data: npt.NDArray[np.float64]) -> float:
    """Example function with type hints.

    Args:
        data: Input array

    Returns:
        Computed result
    """
    return _gromos.rust_function(data)
```

**Style rules:**
- Use type hints
- Follow PEP 8 (enforced by ruff)
- Document with Google-style docstrings
- Line length: 100 characters

### Format Code

```bash
# Rust
cargo fmt --all

# Python - ruff for linting and formatting
ruff check --fix python/ tests/ examples/
ruff format python/ tests/ examples/
```

## Testing

!!! warning "Tests Cannot Run Yet"
    These instructions are for when the project becomes functional.

### Python Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test
pytest tests/test_vec3.py::test_vector_length

# With coverage
pytest tests/ --cov=gromos --cov-report=html
```

### Rust Tests

```bash
# Test PyO3 bindings
cargo test -p py-gromos
```

### Example Tests

```python
# tests/test_vec3.py

import pytest
import gromos
import numpy as np

def test_vector_creation():
    v = gromos.Vec3(1.0, 2.0, 3.0)
    assert v.x == 1.0
    assert v.y == 2.0
    assert v.z == 3.0

def test_vector_length():
    v = gromos.Vec3(3.0, 4.0, 0.0)
    assert abs(v.length() - 5.0) < 1e-10

def test_numpy_conversion():
    v = gromos.Vec3(1.0, 2.0, 3.0)
    arr = v.to_numpy()
    assert isinstance(arr, np.ndarray)
    np.testing.assert_array_equal(arr, [1.0, 2.0, 3.0])
```

## Pull Request Process

### Before Submitting

```bash
# 1. Format code
make fmt

# 2. Lint code
make lint

# 3. Run tests (when possible)
make test

# 4. Check documentation builds
mkdocs build
```

### Creating a PR

1. **Write clear title**:
   ```
   ✓ Add Vec3.normalize() method to PyO3 bindings
   ✗ Fix stuff
   ```

2. **Write detailed description**:
   ```markdown
   ## Summary
   Adds normalize() method to Vec3 class for unit vector calculation.

   ## Changes
   - Add normalize() to Vec3 in src/lib.rs
   - Add Python tests for normalization
   - Add example to docs

   ## Testing
   - [x] Unit tests pass (when gromos-rs compiles)
   - [x] Documentation builds
   - [x] Examples updated
   ```

3. **Keep PRs focused**: One feature/fix per PR

## Documentation

### API Documentation

Document all public APIs:

```rust
/// Calculate vector length (magnitude).
///
/// Returns:
///     float: L2 norm of the vector
///
/// Example:
///     >>> v = gromos.Vec3(3.0, 4.0, 0.0)
///     >>> v.length()
///     5.0
fn length(&self) -> f64 {
    (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
}
```

### Building Docs

```bash
# Install MkDocs
pip install mkdocs mkdocs-material

# Serve locally
mkdocs serve

# Build static site
mkdocs build
```

## Getting Help

- **Questions**: Open a GitHub Discussion
- **Bugs**: Open an Issue with reproduction steps
- **Features**: Open an Issue describing the use case
- **PyO3 Help**: Check [PyO3 documentation](https://pyo3.rs)

## Code of Conduct

Be respectful, professional, and constructive in all interactions.

---

Thank you for contributing to py-gromos! 🦀🐍
