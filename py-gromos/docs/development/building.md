# Building py-gromos

Technical details about building py-gromos from source.

!!! warning "Cannot Currently Build"
    py-gromos cannot be built because gromos-rs has ~140+ compilation errors. This document describes the intended build process.

## Build System Overview

py-gromos uses **Maturin** to build Rust extensions for Python. Maturin handles:

- Compiling Rust code with cargo
- Generating Python-compatible shared libraries
- Creating Python wheels
- Installing in development mode

## Architecture

```
┌─────────────────────────────────────────┐
│         Python Package (py-gromos)       │
│                                          │
│  ┌────────────────────────────────────┐ │
│  │  python/gromos/__init__.py         │ │  Pure Python
│  │  (imports _gromos extension)       │ │
│  └──────────────┬─────────────────────┘ │
│                 │                        │
│  ┌──────────────┴─────────────────────┐ │
│  │  _gromos.so (Rust extension)       │ │  PyO3 bindings
│  │  Compiled from src/lib.rs          │ │
│  └──────────────┬─────────────────────┘ │
└─────────────────┼───────────────────────┘
                  │
┌─────────────────┴───────────────────────┐
│         gromos-rs (Rust library)         │
│  Core MD simulation code                 │
└──────────────────────────────────────────┘
```

## Build Configurations

### Debug Build (Fast Compilation)

```bash
maturin develop
```

- **Optimization**: O0 (no optimization)
- **Compile time**: ~2-3 minutes (when working)
- **Runtime**: Slow
- **Use for**: Development, debugging

### Release Build (Optimized)

```bash
maturin develop --release
```

- **Optimization**: O3 + LTO
- **Compile time**: ~10-15 minutes (when working)
- **Runtime**: Fast
- **Use for**: Performance testing, distribution

## Build Process

### Step 1: Cargo Build

Maturin first compiles the Rust crate:

```bash
cargo build --lib --features pyo3/extension-module
```

This produces a shared library (`.so` on Linux, `.dylib` on macOS, `.dll` on Windows).

### Step 2: Python Integration

Maturin:
1. Renames the library to match Python conventions (`_gromos.so`)
2. Places it in the Python package directory
3. Creates a `.pth` file for development mode

### Step 3: Installation

**Development mode** (`maturin develop`):
- Installs package in editable mode
- Recompile with `maturin develop` after Rust changes
- No need to reinstall after pure Python changes

**Wheel mode** (`maturin build`):
- Creates distributable wheel in `target/wheels/`
- Install with `pip install target/wheels/gromos-*.whl`

## Configuration Files

### Cargo.toml

```toml
[package]
name = "py-gromos"
version = "0.1.0"
edition = "2021"

[lib]
name = "_gromos"           # Underscore prefix for Python extension
crate-type = ["cdylib"]    # Dynamic library for Python

[dependencies]
pyo3 = { version = "0.20", features = ["extension-module", "abi3-py312"] }
gromos-rs = { path = "../gromos-rs" }
numpy = "0.20"
```

**Key settings:**
- `crate-type = ["cdylib"]`: Build as dynamic library
- `extension-module`: PyO3 feature for Python extensions
- `abi3-py312`: Stable ABI (Python 3.12+)

### pyproject.toml

```toml
[build-system]
requires = ["maturin>=1.0,<2.0"]
build-backend = "maturin"

[project]
name = "gromos"
requires-python = ">=3.12"
dependencies = [
    "numpy>=1.20",
]

[tool.maturin]
features = ["pyo3/extension-module"]
```

## Build Targets

### Linux

```bash
# Build for current platform
maturin build --release

# Cross-compile for manylinux (PyPI compatibility)
docker run --rm -v $(pwd):/io konstin2/maturin build --release
```

### macOS

```bash
# Universal binary (Intel + Apple Silicon)
maturin build --release --target universal2-apple-darwin
```

### Windows

```bash
# Build with MSVC
maturin build --release
```

## Development Workflow

### Quick Iteration

For pure Python changes:
```bash
# Edit python/gromos/*.py
# No rebuild needed, changes are live
```

For Rust changes:
```bash
# Edit src/lib.rs
maturin develop  # Recompile (debug mode)
# Test your changes
```

### Testing Build

```bash
# Build in debug mode
maturin develop

# Run tests
pytest tests/

# Check imports work
python -c "import gromos; print(gromos.__version__)"
```

## Optimization Flags

### Rust Optimization

In `.cargo/config.toml`:

```toml
[profile.release]
opt-level = 3          # Maximum optimization
lto = "fat"            # Link-time optimization
codegen-units = 1      # Single codegen unit for better optimization
```

### Python Optimization

```bash
# Run with optimizations
python -O script.py

# Generate optimized bytecode
python -m compileall python/gromos/
```

## Troubleshooting

### Issue: "cannot find -lgromos_rs"

**Cause**: gromos-rs library not built or not in linker path.

**Solution**:
```bash
cd ../gromos-rs
cargo build --release --lib
```

### Issue: "undefined symbol: PyInit__gromos"

**Cause**: Wrong crate-type or missing PyO3 setup.

**Solution**: Check `Cargo.toml` has `crate-type = ["cdylib"]`.

### Issue: Compilation errors in gromos-rs

**Cause**: This is expected. gromos-rs has ~140+ compilation errors.

**Solution**: Help fix gromos-rs compilation errors first.

## Benchmarking Build Time

!!! note "When Functional"
    These are estimates for when the project compiles.

Expected build times (Intel i9-12900K, 32GB RAM):

| Configuration | Clean Build | Incremental |
|---------------|-------------|-------------|
| Debug | ~2 min | ~10 sec |
| Release | ~15 min | ~2 min |
| Release + LTO | ~25 min | ~5 min |

## CI/CD Integration

GitHub Actions workflow (`.github/workflows/test-python.yml`):

```yaml
- name: Build Python package
  run: |
    cd py-gromos
    pip install maturin
    maturin develop --release

- name: Test
  run: |
    pytest py-gromos/tests/ -v
```

## Next Steps

- **[Contributing](contributing.md)**: How to contribute
- **[Quick Start](../user-guide/quick-start.md)**: Using py-gromos
- **[API Reference](../api/reference.md)**: Complete API
