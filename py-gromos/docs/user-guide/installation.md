# Installation Guide

Installation instructions for py-gromos Python bindings.

!!! danger "Not Currently Functional"
    **py-gromos cannot be installed** because gromos-rs has ~140+ compilation errors. These instructions document the intended installation process for when the project becomes functional.

## Prerequisites

### System Requirements

- **Python 3.12 or later**
- **Rust stable (latest)**
- **C compiler** (gcc, clang, or MSVC)
- **Git** (for source installation)

### Platform Support

When functional, py-gromos should support:

- ✅ Linux (x86_64, aarch64)
- ✅ macOS (x86_64, Apple Silicon)
- ✅ Windows (x86_64)

## Installation Methods

### Method 1: From PyPI (Future)

!!! note "Not Available Yet"
    py-gromos is not published to PyPI and cannot be built yet.

```bash
# When available:
pip install gromos
```

### Method 2: From Source (Development)

!!! warning "Will Fail"
    This will fail due to gromos-rs compilation errors.

#### Step 1: Install Rust

```bash
# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Restart your shell, then verify
rustc --version
cargo --version
```

#### Step 2: Clone Repository

```bash
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromos-rs/py-gromos
```

#### Step 3: Install Maturin

Maturin is the build tool for Rust-based Python packages:

```bash
pip install maturin
```

#### Step 4: Build and Install

**Development mode** (editable install):

```bash
maturin develop --release
```

This builds the Rust extension and installs it in your current Python environment.

**Build wheel for distribution**:

```bash
maturin build --release
pip install target/wheels/gromos-*.whl
```

## Verification

!!! note "Cannot Run Yet"
    These verification steps are aspirational.

### Test Installation

```python
import gromos
print(gromos.__version__)

# Create a simple vector
v = gromos.Vec3(1.0, 2.0, 3.0)
print(f"Vector length: {v.length()}")
```

### Run Examples

```bash
cd py-gromos/examples
python 01_vectors_and_math.py
```

## Troubleshooting

### Common Issues

#### Issue: "command 'cargo' not found"

**Solution**: Rust is not installed or not in PATH.

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Add to PATH (Linux/macOS)
source $HOME/.cargo/env

# Verify
cargo --version
```

#### Issue: "maturin: command not found"

**Solution**: Maturin is not installed.

```bash
pip install --upgrade maturin
```

#### Issue: Compilation errors in gromos-rs

**Solution**: This is expected. gromos-rs has ~140+ compilation errors and doesn't build yet. The project is educational and not functional.

## Development Installation

For contributing to py-gromos:

```bash
# Clone repository
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromos-rs/py-gromos

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Install development dependencies
pip install -r requirements-dev.txt

# Install maturin
pip install maturin

# Build in debug mode (faster compilation)
maturin develop

# Run tests (will fail until gromos-rs compiles)
pytest tests/
```

## Uninstallation

```bash
pip uninstall gromos
```

## Next Steps

- **[Quick Start](quick-start.md)**: First steps with py-gromos
- **[API Reference](../api/reference.md)**: Complete API documentation
- **[Building](../development/building.md)**: Build system details
