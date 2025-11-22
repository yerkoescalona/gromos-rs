# CI/CD and Testing Guide

This guide explains the comprehensive CI/CD and testing infrastructure for gromosXX, inspired by the Polars project.

## Table of Contents

- [Overview](#overview)
- [Makefiles](#makefiles)
- [GitHub Actions Workflows](#github-actions-workflows)
- [Local Development](#local-development)
- [Testing](#testing)
- [Release Process](#release-process)

## Overview

The gromosXX project consists of two main components:

1. **gromos-rs**: High-performance Rust library for molecular dynamics
2. **py-gromos**: Python bindings for gromos-rs

Our CI/CD setup ensures both components are thoroughly tested individually and together.

## Makefiles

We provide three Makefiles for streamlined development:

### Root Makefile (`./Makefile`)

The root Makefile orchestrates both projects together.

**Common commands:**

```bash
# Build both projects
make build

# Run all tests
make test

# Format, lint, test, and build everything
make all

# Quick CI check (what runs in CI)
make ci

# Format code in both projects
make fmt

# Lint code in both projects
make lint

# Build release versions
make build-release

# Clean all artifacts
make clean
```

### gromos-rs Makefile (`./gromos-rs/Makefile`)

Rust-specific development commands.

**Common commands:**

```bash
# Build in debug mode (fast compilation)
make build

# Build in release mode (optimized)
make build-release

# Run all tests
make test

# Run tests with all features
make test-all-features

# Run clippy linter
make clippy

# Format code
make fmt

# Run benchmarks
make bench

# Build documentation
make doc

# Pre-commit checks (format + lint + test)
make pre-commit

# Build MPI binaries
make build-md-mpi

# Run MPI tests (requires MPI runtime)
make test-mpi
```

### py-gromos Makefile (`./py-gromos/Makefile`)

Python bindings development commands.

**Common commands:**

```bash
# Create virtual environment
make venv

# Install dev dependencies
make install-dev

# Build Python extension (debug)
make build

# Build Python extension (release)
make build-release

# Run Python tests
make test

# Run tests with coverage
make test-coverage

# Format Python code
make fmt

# Lint Python code
make lint

# Run examples
make run-examples

# Build wheel for distribution
make build-wheel

# Clean all artifacts
make clean
```

**CPU Optimization:**

For older CPUs without AVX2 support:
```bash
LTS_CPU=1 make build-release
```

## GitHub Actions Workflows

### 1. Test Rust (`test-rust.yml`)

Tests the Rust library across platforms and feature combinations.

**Triggers:**
- Push to main or claude/* branches
- Pull requests
- Changes to gromos-rs code

**Test Matrix:**
- OS: Ubuntu, macOS, Windows
- Rust: stable, nightly (experimental)
- Features: various combinations

**Jobs:**
- Code formatting check
- Clippy linting
- Build with all features
- Run tests (debug and release)
- Documentation tests
- Feature combination tests
- Benchmarks (on main branch only)
- Security audit

### 2. Test Python (`test-python.yml`)

Tests Python bindings across platforms and Python versions.

**Triggers:**
- Push to main or claude/* branches
- Pull requests
- Changes to py-gromos or gromos-rs

**Test Matrix:**
- OS: Ubuntu, macOS, Windows
- Python: 3.9, 3.10, 3.11, 3.12, 3.13

**Jobs:**
- Build and test Python bindings
- Format check (black)
- Lint (ruff)
- Clippy for Rust bindings
- Coverage report (Python 3.12 on Ubuntu)
- Import test without dependencies
- Run benchmarks (main branch only)
- Run all examples

### 3. CI (`ci.yml`)

Complete integration testing workflow.

**Triggers:**
- Push to main or claude/* branches
- Pull requests

**Jobs:**
1. **Format Check**: Rust and Python formatting
2. **Lint**: Clippy for Rust, ruff for Python
3. **Integration Tests**: Full test suite for both projects
4. **Build Release**: Build release artifacts on all platforms
5. **Test Notebooks**: Execute Jupyter notebooks (PR only)
6. **Coverage**: Generate combined coverage report (main only)

### 4. Release (`release.yml`)

Automated release process.

**Triggers:**
- Push tags matching `v*.*.*`
- Manual workflow dispatch

**Jobs:**
- Create GitHub release
- Build wheels for all platforms and Python versions
- Build source distribution
- Publish to PyPI (on tag push)
- Publish to crates.io (on tag push)

### 5. Dependency Update (`dependency-update.yml`)

Automated dependency management.

**Triggers:**
- Schedule: Every Monday at 9:00 AM UTC
- Manual workflow dispatch

**Jobs:**
- Update Rust dependencies
- Update Python dependencies
- Run security audit
- Create PR with updates
- Create issue if vulnerabilities found

## Local Development

### First-time Setup

1. **Install Rust:**
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

2. **Install Python 3.9+:**
   ```bash
   # On Ubuntu/Debian
   sudo apt-get install python3 python3-venv python3-pip
   ```

3. **Clone and build:**
   ```bash
   git clone <repository>
   cd gromosXX
   make build
   ```

### Daily Development Workflow

```bash
# 1. Make your changes to the code

# 2. Format code
make fmt

# 3. Run linters
make lint

# 4. Run tests
make test

# 5. Before committing, run full checks
make pre-commit

# 6. Commit and push
git add .
git commit -m "Your message"
git push
```

### Testing Specific Components

```bash
# Test only Rust code
make test-rust

# Test only Python bindings
make test-python

# Test with specific features
cd gromos-rs && make test-mpi

# Run specific Python test
cd py-gromos && .venv/bin/pytest tests/test_basic.py -v

# Run tests with coverage
make test-coverage
```

## Testing

### Test Organization

**gromos-rs tests:**
- `tests/comprehensive_tests.rs` - Core functionality
- `tests/md_*.rs` - Molecular dynamics features
- `tests/io_*.rs` - File I/O operations
- `tests/mpi_tests.rs` - MPI parallel features
- `tests/validation.rs` - Validation against reference

**py-gromos tests:**
- `tests/test_basic.py` - Basic Python API
- `tests/test_advanced_features.py` - Advanced features
- `examples/*.py` - Examples that serve as integration tests
- `notebooks/*.ipynb` - Jupyter notebooks

### Running Tests Locally

```bash
# Quick test (debug mode)
make test

# Comprehensive test (all features, release mode)
make test-all

# Test with coverage
make test-coverage

# Run specific test suite
cd gromos-rs && cargo test md_integration

# Run Python tests with verbose output
cd py-gromos && make test-verbose

# Run benchmarks
cd gromos-rs && make bench
```

### Writing Tests

**Rust tests:**
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feature() {
        // Your test code
        assert_eq!(result, expected);
    }
}
```

**Python tests:**
```python
import pytest
import gromos

def test_feature():
    # Your test code
    assert result == expected
```

## Release Process

### Creating a Release

1. **Update version numbers:**
   - `gromos-rs/Cargo.toml`
   - `py-gromos/pyproject.toml`

2. **Update CHANGELOG.md** with release notes

3. **Commit changes:**
   ```bash
   git add .
   git commit -m "chore: Bump version to 0.2.0"
   git push
   ```

4. **Create and push tag:**
   ```bash
   git tag -a v0.2.0 -m "Release version 0.2.0"
   git push origin v0.2.0
   ```

5. **GitHub Actions will automatically:**
   - Create GitHub release
   - Build wheels for all platforms
   - Publish to PyPI (requires PYPI_TOKEN secret)
   - Publish to crates.io (requires CARGO_REGISTRY_TOKEN secret)

### Manual Release

If needed, you can build and publish manually:

```bash
# Build Python wheel
cd py-gromos
make build-wheel

# Publish to PyPI
make publish

# Publish Rust crate
cd ../gromos-rs
cargo publish
```

## CI/CD Best Practices

1. **Always run `make pre-commit` before pushing**
2. **Write tests for new features**
3. **Keep dependencies updated** (automated weekly)
4. **Review security audit results**
5. **Maintain code coverage above 80%**
6. **Document new features and APIs**
7. **Use conventional commits** for clear history

## Troubleshooting

### Common Issues

**Issue: Tests fail locally but pass in CI**
- Ensure you're using the same Rust/Python versions
- Clean and rebuild: `make clean && make build`

**Issue: MPI tests fail**
- Ensure MPI runtime is installed: `sudo apt-get install libopenmpi-dev`
- Run with: `mpirun -n 4 cargo test --features use-mpi`

**Issue: Python import fails**
- Rebuild Python bindings: `cd py-gromos && make build`
- Check virtual environment is activated

**Issue: Clippy errors**
- Auto-fix: `make fix`
- Manual review may be needed for some warnings

## Additional Resources

- [Cargo Book](https://doc.rust-lang.org/cargo/)
- [Maturin Documentation](https://www.maturin.rs/)
- [PyO3 Guide](https://pyo3.rs/)
- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Polars CI/CD](https://github.com/pola-rs/polars) - Inspiration for this setup

## Contributing

When contributing:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Make your changes
4. Run `make pre-commit` to ensure code quality
5. Commit your changes: `git commit -m 'feat: Add amazing feature'`
6. Push to the branch: `git push origin feature/amazing-feature`
7. Open a Pull Request

All PRs must pass CI checks before merging.
