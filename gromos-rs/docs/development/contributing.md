# Contributing to GROMOS-RS

Thank you for your interest in contributing to GROMOS-RS! This guide will help you get started with contributing code, documentation, and ideas to the project.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Development Setup](#development-setup)
3. [Code Style and Standards](#code-style-and-standards)
4. [Testing Guidelines](#testing-guidelines)
5. [Pull Request Process](#pull-request-process)
6. [Code Review Guidelines](#code-review-guidelines)
7. [Common Contribution Types](#common-contribution-types)

## Getting Started

### Prerequisites

Before contributing, ensure you have:

```bash
# Rust (1.70+)
rustup update stable

# Development tools
cargo install cargo-clippy
cargo install cargo-fmt
cargo install cargo-criterion

# Optional but recommended
cargo install cargo-watch    # Auto-rebuild on changes
cargo install cargo-expand   # Expand macros
cargo install cargo-tree     # Dependency tree
```

### Fork and Clone

```bash
# 1. Fork the repository on GitHub

# 2. Clone your fork
git clone https://github.com/YOUR_USERNAME/gromosXX.git
cd gromosXX/gromos-rs

# 3. Add upstream remote
git remote add upstream https://github.com/yerkoescalona/gromos-rs.git

# 4. Create feature branch
git checkout -b feature/my-contribution
```

## Development Setup

### Quick Start

```bash
# Build in debug mode
cargo build

# Run tests
cargo test

# Run with logging
RUST_LOG=debug cargo run --bin md -- --help

# Watch for changes (if installed cargo-watch)
cargo watch -x test
```

### Project Structure

Understand the layout:

```
gromos-rs/
├── src/
│   ├── lib.rs              # Library API
│   ├── bin/                # Executables (md, remd, eds, etc.)
│   ├── algorithm/          # Algorithms (constraints, thermostats)
│   ├── interaction/        # Force calculations
│   ├── io/                 # File I/O
│   └── ...
├── tests/                  # Integration tests
├── benches/                # Performance benchmarks
├── examples/               # Usage examples
└── docs/                   # Documentation
```

## Code Style and Standards

!!! info "Complete Style Guide"
    For comprehensive coding conventions, see the [Coding Style Guide](coding-style.md).
    This section provides a quick overview.

### Formatting

**Always** use `rustfmt` before committing:

```bash
# Format all code
cargo fmt

# Check formatting (CI check)
cargo fmt -- --check
```

Configuration in `.rustfmt.toml`:

```toml
max_width = 100
tab_spaces = 4
use_small_heuristics = "Default"
```

### Linting with Clippy

**Always** run Clippy and fix warnings:

```bash
# Run clippy
cargo clippy -- -D warnings

# Auto-fix some issues
cargo clippy --fix
```

**Common Clippy Rules**:

```rust
// ✗ Bad: Unnecessary clone
let v = vec.clone();
process(&v);

// ✓ Good: Borrow instead
process(&vec);

// ✗ Bad: Unnecessary collect
let sum: i32 = vec.iter().collect::<Vec<_>>().iter().sum();

// ✓ Good: Direct iteration
let sum: i32 = vec.iter().sum();

// ✗ Bad: Manual iterator
for i in 0..vec.len() {
    process(&vec[i]);
}

// ✓ Good: Use iterator
for item in &vec {
    process(item);
}
```

### Naming Conventions

Follow Rust naming conventions:

```rust
// ✓ Types: UpperCamelCase
pub struct MyStruct;
pub enum MyEnum;
pub trait MyTrait;

// ✓ Functions, variables: snake_case
pub fn calculate_energy() {}
let my_variable = 42;

// ✓ Constants: SCREAMING_SNAKE_CASE
pub const MAX_ITERATIONS: usize = 1000;

// ✓ Modules: snake_case
mod my_module;

// ✓ Lifetimes: short, lowercase
fn foo<'a, 'b>(x: &'a str, y: &'b str) {}
```

### Documentation

All public APIs **must** be documented:

```rust
/// Calculate Lennard-Jones interaction energy and forces.
///
/// Uses the 12-6 Lennard-Jones potential:
/// \[ V(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right] \]
///
/// # Arguments
///
/// * `r2` - Squared distance between atoms (nm²)
/// * `epsilon` - LJ well depth (kJ/mol)
/// * `sigma` - LJ size parameter (nm)
///
/// # Returns
///
/// Tuple of (energy in kJ/mol, force magnitude in kJ/mol/nm)
///
/// # Examples
///
/// ```
/// use gromos_rs::interaction::lj;
///
/// let r = 0.35; // nm
/// let (energy, force) = lj::calculate(r*r, 0.65, 0.32);
/// assert!(energy < 0.0); // Attractive at this distance
/// ```
///
/// # Performance
///
/// This function is designed to be #[inline] for performance.
///
/// # See Also
///
/// * [`calculate_lj_crf`] - Combined LJ + Coulomb calculation
#[inline]
pub fn calculate(r2: f64, epsilon: f64, sigma: f64) -> (f64, f64) {
    let sigma2 = sigma * sigma;
    let sigma6 = sigma2 * sigma2 * sigma2;
    let r6 = r2 * r2 * r2;

    let energy = 4.0 * epsilon * (sigma6 * sigma6 / (r6 * r6) - sigma6 / r6);
    let force = 24.0 * epsilon * (2.0 * sigma6 * sigma6 / (r6 * r6) - sigma6 / r6) / r2;

    (energy, force)
}
```

**Documentation Guidelines**:

1. **Summary**: First line should be a one-line summary
2. **Details**: Explain what the function does, not how
3. **Arguments**: Document all parameters with types and units
4. **Returns**: Explain return value with units
5. **Examples**: Include runnable examples (tested by `cargo test`)
6. **Math**: Use LaTeX for equations (rendered in docs)
7. **Performance**: Note if optimized/SIMD/parallel
8. **See Also**: Link to related functions

### Error Handling

Use `Result` for fallible operations:

```rust
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SimulationError {
    #[error("Topology file not found: {path}")]
    TopologyNotFound { path: String },

    #[error("SHAKE failed to converge after {iterations} iterations")]
    ShakeConvergence { iterations: usize },

    #[error("Invalid atom index: {index} (max: {max})")]
    InvalidAtomIndex { index: usize, max: usize },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, SimulationError>;

// Usage:
pub fn read_topology(path: &str) -> Result<Topology> {
    let file = std::fs::File::open(path)
        .map_err(|_| SimulationError::TopologyNotFound {
            path: path.to_string(),
        })?;

    // ... parse topology

    Ok(topology)
}
```

**Error Handling Rules**:

1. **Never `unwrap()` in library code** (only in tests/examples)
2. **Never `panic!()` for expected errors** (use `Result`)
3. **Use `thiserror`** for error types
4. **Provide context** in error messages
5. **Propagate with `?`** operator

### Performance-Critical Code

Mark hot paths and justify optimizations:

```rust
/// Fast inner loop for nonbonded force calculation.
///
/// # Performance
///
/// This function is expected to be the hottest path in MD simulations.
/// Optimizations applied:
/// - #[inline(always)] for zero-cost abstraction
/// - SIMD vectorization (AVX2/AVX-512)
/// - Manual loop unrolling for better instruction pipelining
/// - Careful memory layout for cache efficiency
///
/// Note: Actual performance has not been measured in gromos-rs.
#[inline(always)]
#[target_feature(enable = "avx2")]
unsafe fn lj_crf_inner_loop_simd(/* ... */) {
    // Highly optimized code with comments explaining why
}
```

**Performance Guidelines**:

1. **Profile first**: Use `cargo bench` and `perf` (when project compiles)
2. **Document optimizations**: Explain why code is complex
3. **Benchmark changes**: Include before/after numbers (when possible)
4. **Keep simple version**: For reference/testing
5. **Use `#[inline]` carefully**: Only for small, hot functions

!!! warning "Current Status"
    gromos-rs has ~140+ compilation errors and cannot currently run benchmarks. Performance optimizations are aspirational until the codebase is functional.

## Testing Guidelines

### Unit Tests

Every module should have tests:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_lj_energy() {
        let epsilon = 0.65; // kJ/mol
        let sigma = 0.32;   // nm
        let r = sigma * 2.0_f64.powf(1.0/6.0); // Minimum energy distance

        let (energy, _force) = calculate_lj(r*r, epsilon, sigma);

        assert_relative_eq!(energy, -epsilon, epsilon = 1e-10);
    }

    #[test]
    fn test_lj_force_zero_at_minimum() {
        let epsilon = 0.65;
        let sigma = 0.32;
        let r = sigma * 2.0_f64.powf(1.0/6.0);

        let (_energy, force) = calculate_lj(r*r, epsilon, sigma);

        assert_relative_eq!(force, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_lj_attractive_below_minimum() {
        let epsilon = 0.65;
        let sigma = 0.32;
        let r = sigma * 2.0_f64.powf(1.0/6.0) * 1.1; // Slightly beyond minimum

        let (energy, force) = calculate_lj(r*r, epsilon, sigma);

        assert!(energy < 0.0); // Attractive
        assert!(force > 0.0);  // Pulls together
    }

    #[test]
    fn test_lj_repulsive_above_minimum() {
        let epsilon = 0.65;
        let sigma = 0.32;
        let r = sigma * 0.9; // Inside minimum

        let (energy, force) = calculate_lj(r*r, epsilon, sigma);

        assert!(energy > 0.0); // Repulsive
        assert!(force < 0.0);  // Pushes apart
    }

    #[test]
    #[should_panic(expected = "distance cannot be zero")]
    fn test_lj_zero_distance_panics() {
        calculate_lj(0.0, 0.65, 0.32);
    }
}
```

**Testing Best Practices**:

1. **Test one thing**: Each test should verify one behavior
2. **Use `approx`**: For floating-point comparisons
3. **Test edge cases**: Zero, negative, very large, very small
4. **Test errors**: Use `#[should_panic]` or `assert!(result.is_err())`
5. **Use descriptive names**: `test_what_when_expected`

### Integration Tests

```rust
// tests/integration_md.rs

use gromos_rs::*;

#[test]
fn test_harmonic_oscillator_energy_conservation() {
    // Create 1D harmonic oscillator
    let mut sim = create_harmonic_oscillator(
        mass = 1.0,
        k = 100.0,
        x0 = 0.1,
    );

    let integrator = LeapFrog { dt: 0.001 };
    let e_initial = sim.total_energy();

    // Run 10000 steps (10 ps)
    for _ in 0..10000 {
        integrator.step(&mut sim).unwrap();
    }

    let e_final = sim.total_energy();
    let drift = (e_final - e_initial).abs() / e_initial;

    // Energy should be conserved to <0.01%
    assert!(drift < 1e-4, "Energy drift: {:.6}%", drift * 100.0);
}

#[test]
fn test_water_shake_constraints() {
    let mut water = create_water_molecule();
    let shake = Shake::new(tolerance = 1e-4, max_iter = 1000);

    // Apply random velocities
    randomize_velocities(&mut water, temp = 300.0);

    // Apply SHAKE
    shake.apply(&mut water).unwrap();

    // Check bond lengths satisfied
    let r_oh1 = distance(&water, 0, 1);
    let r_oh2 = distance(&water, 0, 2);
    let r_hh = distance(&water, 1, 2);

    assert_relative_eq!(r_oh1, 0.1, epsilon = 1e-6); // 0.1 nm = 1 Å
    assert_relative_eq!(r_oh2, 0.1, epsilon = 1e-6);
    assert_relative_eq!(r_hh, 0.163299, epsilon = 1e-6); // 104.5° angle
}
```

### Property-Based Testing

For complex algorithms, use property testing:

```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn test_pbc_nearest_image_is_shortest(
        r in -10.0..10.0,
        box_length in 1.0..20.0
    ) {
        let r_pbc = apply_pbc(r, box_length);

        // Property: PBC image should be within half box length
        prop_assert!(r_pbc.abs() <= box_length / 2.0);

        // Property: Equivalent to original modulo box_length
        let diff = (r - r_pbc).abs();
        prop_assert!(
            diff < 1e-10 ||
            (diff - box_length).abs() < 1e-10 ||
            (diff - 2.0 * box_length).abs() < 1e-10
        );
    }
}
```

### Benchmarks

Add benchmarks for performance-critical code:

```rust
// benches/nonbonded.rs

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use gromos_rs::interaction::nonbonded::*;

fn bench_nonbonded(c: &mut Criterion) {
    let mut group = c.benchmark_group("nonbonded");

    for size in [100, 1000, 10000].iter() {
        let (positions, pairlist) = create_test_system(*size);

        group.bench_with_input(
            BenchmarkId::new("serial", size),
            size,
            |b, _| {
                b.iter(|| {
                    calculate_nonbonded_serial(
                        black_box(&positions),
                        black_box(&pairlist),
                    )
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("parallel", size),
            size,
            |b, _| {
                b.iter(|| {
                    calculate_nonbonded_parallel(
                        black_box(&positions),
                        black_box(&pairlist),
                    )
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_nonbonded);
criterion_main!(benches);
```

Run benchmarks:

```bash
cargo bench

# View results
open target/criterion/report/index.html
```

## Pull Request Process

### Before Submitting

Checklist before creating a PR:

```bash
# 1. Format code
cargo fmt

# 2. Run clippy
cargo clippy -- -D warnings

# 3. Run all tests
cargo test --all-features

# 4. Run benchmarks (if you changed performance-critical code)
cargo bench

# 5. Update documentation
cargo doc --no-deps --open

# 6. Check documentation examples
cargo test --doc
```

### Creating a PR

1. **Write a clear title**:
   ```
   ✓ Good: Add Perturbed SHAKE for FEP simulations
   ✗ Bad: Fix stuff
   ```

2. **Write a detailed description**:
   ```markdown
   ## Summary
   Implements Perturbed SHAKE algorithm for λ-dependent constraints in FEP simulations.

   ## Changes
   - Add `PerturbedShake` struct with λ-dependent constraint lengths
   - Implement `apply()` method with iterative solver
   - Add unit tests for convergence and constraint satisfaction
   - Add integration test with water molecule
   - Add documentation with mathematical description

   ## Performance
   - Expected to converge in 3-5 iterations (similar to regular SHAKE)
   - Designed for minimal overhead when λ=0 or λ=1

   ## Testing
   - [x] Unit tests pass
   - [x] Integration tests pass
   - [x] Benchmarks show no regression
   - [x] Documentation builds
   - [x] Examples run successfully

   ## Related Issues
   Closes #123
   ```

3. **Keep PRs focused**: One feature/fix per PR

4. **Update CHANGELOG.md**: Add entry under "Unreleased"

5. **Add examples** (if adding new features)

### Code Review Process

Your PR will be reviewed for:

1. **Correctness**: Does it work? Are there bugs?
2. **Testing**: Are there sufficient tests?
3. **Documentation**: Is the code well-documented?
4. **Performance**: Does it maintain/improve performance?
5. **Style**: Does it follow Rust conventions?
6. **Safety**: No unsafe code without justification

**Responding to Reviews**:

- Be receptive to feedback
- Ask questions if unclear
- Make requested changes in new commits (don't force-push until approved)
- Respond to each comment

## Common Contribution Types

### Adding a New Force Term

Example: Adding harmonic dihedral potential

```rust
// 1. Add to interaction/bonded.rs

#[derive(Debug, Clone)]
pub struct HarmonicDihedral {
    pub atoms: [usize; 4],  // i-j-k-l
    pub k: f64,             // Force constant (kJ/mol/rad²)
    pub phi0: f64,          // Equilibrium angle (rad)
}

impl HarmonicDihedral {
    pub fn calculate(
        &self,
        positions: &[Vector3<f64>],
        forces: &mut [Vector3<f64>],
    ) -> f64 {
        // Get atom positions
        let r_i = positions[self.atoms[0]];
        let r_j = positions[self.atoms[1]];
        let r_k = positions[self.atoms[2]];
        let r_l = positions[self.atoms[3]];

        // Calculate dihedral angle
        let phi = calculate_dihedral_angle(r_i, r_j, r_k, r_l);

        // Energy: V = 0.5 * k * (φ - φ₀)²
        let delta_phi = phi - self.phi0;
        let energy = 0.5 * self.k * delta_phi * delta_phi;

        // Force: F = -dV/dr = -k * (φ - φ₀) * dφ/dr
        let f_magnitude = -self.k * delta_phi;
        let gradients = calculate_dihedral_gradients(r_i, r_j, r_k, r_l);

        forces[self.atoms[0]] += f_magnitude * gradients[0];
        forces[self.atoms[1]] += f_magnitude * gradients[1];
        forces[self.atoms[2]] += f_magnitude * gradients[2];
        forces[self.atoms[3]] += f_magnitude * gradients[3];

        energy
    }
}

// 2. Add tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_harmonic_dihedral_equilibrium() {
        // At equilibrium, energy should be zero and force zero
        let dihedral = HarmonicDihedral {
            atoms: [0, 1, 2, 3],
            k: 100.0,
            phi0: std::f64::consts::PI,
        };

        // Create trans configuration
        let positions = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
            Vector3::new(0.2, 0.0, 0.0),
            Vector3::new(0.3, 0.0, 0.0),
        ];

        let mut forces = vec![Vector3::zero(); 4];
        let energy = dihedral.calculate(&positions, &mut forces);

        assert_relative_eq!(energy, 0.0, epsilon = 1e-10);
        for f in &forces {
            assert_relative_eq!(f.norm(), 0.0, epsilon = 1e-10);
        }
    }
}

// 3. Add documentation
// 4. Add to BondedForces struct
// 5. Add example
// 6. Submit PR!
```

### Adding a New Integrator

See [`rust-for-gromos-rs.md`](rust-for-gromos-rs.md) for detailed examples.

### Improving Performance

!!! warning "Not Currently Applicable"
    gromos-rs has ~140+ compilation errors and cannot build in release mode. Performance work should wait until the codebase is functional.

When the project compiles:

1. **Profile first**:
   ```bash
   cargo build --release
   perf record --call-graph dwarf ./target/release/md
   perf report
   ```

2. **Benchmark before and after**:
   ```bash
   # Before
   cargo bench --bench nonbonded > before.txt

   # Make changes

   # After
   cargo bench --bench nonbonded > after.txt

   # Compare
   cargo install cargo-benchcmp
   cargo-benchcmp before.txt after.txt
   ```

3. **Document the optimization** in PR

### Adding Documentation

- Add to `docs/` directory
- Follow existing structure
- Include code examples
- Add to `mkdocs.yml` navigation

## Getting Help

- **Questions**: Open a GitHub Discussion
- **Bugs**: Open an Issue with reproduction steps
- **Features**: Open an Issue describing the use case

## Code of Conduct

Be respectful, professional, and constructive in all interactions.

---

Thank you for contributing to GROMOS-RS! 🦀
