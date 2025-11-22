# Rust API Documentation

This page provides links and information for accessing the full Rust API documentation.

## Generating Documentation

GROMOS-RS uses Rust's built-in documentation system (`rustdoc`).

### Generate Locally

```bash
cd gromos-rs

# Generate documentation for gromos-rs only
cargo doc --no-deps

# Generate documentation including dependencies
cargo doc

# Generate and open in browser
cargo doc --open
```

Documentation will be generated in `target/doc/gromos_rs/index.html`.

### Documentation Structure

The generated documentation includes:

- **Modules**: Organized by functionality
- **Structs/Enums**: Data types with field documentation
- **Functions**: Standalone functions with examples
- **Traits**: Interfaces and trait implementations
- **Examples**: Code examples for common tasks
- **Source code**: Links to implementation

## Online Documentation

!!! info "Coming Soon"
    When gromos-rs is published to crates.io, automatic documentation will be available at [https://docs.rs/gromos-rs](https://docs.rs/gromos-rs)

## Key Modules

Click through to explore each module in the generated docs:

### Core Simulation

- **`gromos_rs::integrator`**: Integration algorithms
  - Leap-Frog, Velocity Verlet
  - Stochastic Dynamics
  - Energy minimization

- **`gromos_rs::topology`**: System topology
  - Atoms, bonds, molecules
  - Force field parameters
  - Exclusions and constraints

- **`gromos_rs::config`**: Simulation configuration
  - Parameter file parsing
  - Simulation settings
  - Output control

### Forces and Interactions

- **`gromos_rs::interaction::bonded`**: Bonded forces
  - Bonds (quartic, harmonic)
  - Angles (cosine, harmonic)
  - Dihedrals (proper, improper, cross)

- **`gromos_rs::interaction::nonbonded`**: Nonbonded forces
  - Lennard-Jones
  - Coulomb
  - Soft-core (FEP)

- **`gromos_rs::interaction::electrostatics`**: Long-range electrostatics
  - Reaction Field
  - PME (Particle Mesh Ewald)
  - Ewald summation

- **`gromos_rs::interaction::restraints`**: Restraining potentials
  - Position, distance, angle, dihedral
  - Harmonic and flat-bottom

### Algorithms

- **`gromos_rs::algorithm::constraints`**: Constraint algorithms
  - SHAKE, M-SHAKE
  - SETTLE (water)
  - LINCS

- **`gromos_rs::algorithm::thermostats`**: Temperature coupling
  - Berendsen
  - Nosé-Hoover
  - Andersen

- **`gromos_rs::algorithm::barostats`**: Pressure coupling
  - Berendsen
  - Parrinello-Rahman

### Advanced Sampling

- **`gromos_rs::fep`**: Free energy perturbation
  - λ-coupling
  - Soft-core potentials
  - dH/dλ calculation

- **`gromos_rs::remd`**: Replica exchange MD
  - Temperature REMD
  - λ-REMD (Hamiltonian)
  - 2D REMD
  - Exchange statistics

- **`gromos_rs::eds`**: Enveloping distribution sampling
  - Multi-state EDS
  - Adaptive offsets (AEDS)
  - Reference state updates

- **`gromos_rs::gamd`**: Gaussian accelerated MD
  - Boost potential
  - Dual boost mode
  - Statistics collection

### I/O

- **`gromos_rs::io::topology`**: Topology file I/O
  - `.top` reader
  - GROMOS96 format

- **`gromos_rs::io::coordinate`**: Coordinate file I/O
  - `.cnf` reader/writer
  - Position, velocity, box

- **`gromos_rs::io::trajectory`**: Trajectory I/O
  - `.trc` writer (ASCII and binary)
  - Frame buffering

- **`gromos_rs::io::energy`**: Energy output
  - `.tre` writer
  - Energy component tracking

- **`gromos_rs::io::force`**: Force output
  - `.trf` writer
  - Detailed force breakdown

- **`gromos_rs::io::ptp`**: Perturbation topology
  - `.ptp` writer
  - FEP state definitions

- **`gromos_rs::io::dlg`**: Free energy output
  - `.dlg` writer (dH/dλ)
  - TI (thermodynamic integration)

### Utilities

- **`gromos_rs::pairlist`**: Neighbor list generation
  - Standard pairlist
  - Grid cell pairlist (O(N))

- **`gromos_rs::math`**: Mathematical utilities
  - Vector operations
  - Matrix operations
  - Periodic boundary conditions

- **`gromos_rs::error`**: Error types
  - `GromosError`
  - Result types

## Example: Finding Documentation

### Search for a Function

In the generated docs, use the search bar:

1. Open `cargo doc --open`
2. Type "leap_frog" in search
3. Navigate to `gromos_rs::integrator::LeapFrog`

### Browse by Module

1. Start at the main page: `gromos_rs`
2. Click "Modules" in the left sidebar
3. Navigate to module of interest
4. Click through to structs/functions

### View Source Code

All documentation includes "source" links:

1. Navigate to any function/struct
2. Click `[src]` link in top-right
3. View actual implementation

## Building Custom Documentation

### Include Private Items

```bash
cargo doc --document-private-items
```

### Exclude Dependencies

```bash
cargo doc --no-deps
```

### For a Specific Package

```bash
cargo doc -p gromos-rs --no-deps
```

## Documentation Standards

All public APIs in GROMOS-RS follow these documentation standards:

### Function Documentation

```rust
/// Brief one-line description.
///
/// More detailed explanation of what the function does,
/// including any important caveats or considerations.
///
/// # Arguments
///
/// * `param1` - Description of first parameter
/// * `param2` - Description of second parameter
///
/// # Returns
///
/// Description of return value
///
/// # Errors
///
/// Description of when this function returns an error
///
/// # Examples
///
/// ```
/// use gromos_rs::module::function;
///
/// let result = function(arg1, arg2)?;
/// assert_eq!(result, expected);
/// ```
///
/// # Panics
///
/// Description of conditions that cause panics (if any)
///
/// # Safety
///
/// For `unsafe` functions, explain safety requirements
pub fn function(param1: Type1, param2: Type2) -> Result<ReturnType> {
    // ...
}
```

### Struct Documentation

```rust
/// Brief description of the struct.
///
/// Detailed explanation of purpose and usage.
///
/// # Examples
///
/// ```
/// use gromos_rs::Struct;
///
/// let s = Struct::new();
/// s.method();
/// ```
pub struct Struct {
    /// Documentation for field1
    pub field1: Type1,

    /// Documentation for field2
    field2: Type2,  // Private field
}
```

## Testing Documentation Examples

All code examples in documentation are tested:

```bash
# Run documentation tests
cargo test --doc

# Run all tests including doc tests
cargo test
```

This ensures examples stay up-to-date!

## Contributing to Documentation

See [Contributing Guide](../development/contributing.md) for guidelines on:

- Writing clear documentation
- Adding examples
- Documenting new features
- Keeping docs up-to-date

## Additional Resources

- **The Rust Programming Language Book**: [https://doc.rust-lang.org/book/](https://doc.rust-lang.org/book/)
- **rustdoc Book**: [https://doc.rust-lang.org/rustdoc/](https://doc.rust-lang.org/rustdoc/)
- **API Guidelines**: [https://rust-lang.github.io/api-guidelines/](https://rust-lang.github.io/api-guidelines/)
