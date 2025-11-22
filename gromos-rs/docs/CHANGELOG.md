# Changelog

All notable changes to GROMOS-RS will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

!!! warning "Educational Project"
    This is an educational project in early development. No production releases exist yet.

## [0.0.1a1] - 2025-11-22

### Added

- Initial alpha release of GROMOS-RS
- Core molecular dynamics simulation framework in Rust
- Python bindings via PyO3/Maturin (py-gromos package)
- Basic integrator implementations (MD, Leapfrog)
- Topology and configuration file parsing
- GPU acceleration support structure (CUDA)
- MPI parallelization framework
- Advanced sampling methods scaffolding (EDS, REMD, GaMD)
- Free energy perturbation (FEP) modules
- Comprehensive test suite and CI/CD pipeline
- Developer documentation and contribution guidelines

### Known Limitations

- Project in early alpha stage - not production ready
- ~140+ compilation errors with all features enabled
- Missing implementations for PME, complete MPI, EDS, and REMD
- No validation against GROMOS reference implementations
- Educational/research purposes only

## [Unreleased]

### Added (Infrastructure)
- GPL-2.0 license with proper attribution to GROMOS/Biomos
- Comprehensive CI/CD with GitHub Actions
- Developer tooling (rustfmt, clippy, black, ruff, pre-commit hooks)
- MkDocs documentation structure
- Makefile for unified build commands
- Python bindings skeleton with PyO3/Maturin
- Comprehensive .gitignore files for Rust and Python

### Added (Code Structure)
- Rust workspace with gromos-rs library
- Python package py-gromos with bindings
- Test data relocated to gromos-rs/tests/data/
- Educational disclaimer in README
- Developer documentation in gromos-rs/docs/development/

### Known Issues
- ~140+ compilation errors with all features enabled
- Missing type definitions (CudaError → DriverError needed fixing)
- Incomplete PME/MPI implementations
- Unfinished EDS/REMD functionality
- No validation against GROMOS reference implementations

### Removed
- Legacy C++ code (md++, gromosPlusPlus) - 1.5M lines removed
- 22 outdated planning/analysis markdown files - 8,800 lines removed
- Unused test_tutorial and tutorial_files directories
- Empty /docs directory (moved to gromos-rs/docs/)

### Status
**This project is NOT ready for any use - educational or otherwise.**
It serves as a learning resource for understanding the structure of a
Rust/Python molecular dynamics project.

## Version History

No releases have been made yet. This is a work in progress.

[Unreleased]: https://github.com/yerkoescalona/gromos-rs/tree/main
