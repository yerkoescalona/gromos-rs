# Coding Style Guide

This guide documents the coding conventions and style for GROMOS-RS based on codebase analysis.

## General Principles

### 1. Scientific Domain First

We respect established scientific notation rather than forcing programming conventions:

```rust
// ✓ GOOD: Domain-standard notation
let kT = BOLTZMANN * temperature;  // Thermal energy
let pH = 7.4;                        // Hydrogen ion concentration
let deltaG = -50.0;                  // Free energy change (kJ/mol)

// ✗ AVOID: Forced programming conventions
let k_t = BOLTZMANN * temperature;   // Non-standard in physics
let p_h = 7.4;                        // Nobody writes it this way
```

See [NAMING_CONVENTIONS.md](/NAMING_CONVENTIONS.md) for full scientific naming guide.

### 2. Readability Over Cleverness

```rust
// ✓ GOOD: Clear intent
fn calculate_kinetic_energy(masses: &[f64], velocities: &[Vec3]) -> f64 {
    masses.iter()
        .zip(velocities.iter())
        .map(|(&m, v)| 0.5 * m * v.length_squared() as f64)
        .sum()
}

// ✗ AVOID: Overly clever
fn ke(m: &[f64], v: &[Vec3]) -> f64 {
    m.iter().zip(v).fold(0., |a, (&m, v)| a + 0.5 * m * (v.x * v.x + v.y * v.y + v.z * v.z) as f64)
}
```

### 3. Document the Why, Not the What

```rust
// ✓ GOOD: Explains reasoning
// Use half-step velocities for better energy conservation (Velocity Verlet)
let v_half = v + 0.5 * a * dt;

// ✗ AVOID: States the obvious
// Add 0.5 * a * dt to v
let v_half = v + 0.5 * a * dt;
```

## Naming Conventions

### Variables

#### Standard Names

```rust
// ✓ Descriptive snake_case
let atom_count = topology.num_atoms();
let force_constant = 1000.0;  // kJ/(mol·nm²)
let cutoff2_short = cutoff_short * cutoff_short;
let half_box = box_length / 2.0;
```

#### Scientific Variables

```rust
// ✓ Domain-standard notation (always allowed)
let kT = 2.479;           // kJ/mol
let pH = 7.4;             // Hydrogen ion concentration
let pKa = 4.76;           // Acid dissociation constant
let deltaG = -50.0;       // kJ/mol
let deltaH = -55.0;       // kJ/mol
let rmsd = 0.15;          // nm (already snake_case)
```

#### Standard Abbreviations

Use these consistently throughout the codebase:

```rust
conf    // configuration (not config, cfg, or configuration)
topo    // topology (not top)
pos     // position (not position, p)
vel     // velocity (not velocity, v)
coord   // coordinate (not coordinate)
dih     // dihedral (not dihedral - prefer short form)
elec    // electrostatic (not electrostatic)
crf     // charge reaction field (lowercase in code, CRF in docs)
lj      // Lennard-Jones (lowercase in code, LJ in docs)
pme     // Particle Mesh Ewald (lowercase in code, PME in docs)
```

#### Index Variables

```rust
// ✓ Single letter for simple loops
for i in 0..n_atoms {
    let atom = &atoms[i];
}

// ✓ Descriptive when multiple indices
for (atom_idx, atom) in atoms.iter().enumerate() {
    for (neighbor_idx, neighbor) in neighbors[atom_idx].iter().enumerate() {
        // ...
    }
}

// ✓ Clear context-specific names
for (atom_i, atom_j) in pairs {
    let r_ij = pos[atom_i] - pos[atom_j];
}
```

#### Temporary Variables

```rust
// ✓ GOOD: Consistent naming for same concept
let distance = (pos_i - pos_j).length();
let distance2 = distance * distance;
let inv_distance = 1.0 / distance;

// ✗ AVOID: Mixing r, dist, distance for same thing
let r = pos_i - pos_j;
let dist = r.length();
let distance2 = dist * dist;
```

### Functions

#### Verb Choice (Standardized)

```rust
// ✓ GOOD: Use "calculate_" for computations
pub fn calculate_energy(&self, positions: &[Vec3]) -> f64 { }
pub fn calculate_forces(&mut self, conf: &Configuration) { }
pub fn calculate_pressure(&self) -> f64 { }

// ✗ AVOID: "compute_" (less common in codebase)
pub fn compute_energy(&self) -> f64 { }  // Inconsistent

// ✓ GOOD: Use "update_" for mutations
pub fn update_statistics(&mut self, value: f64) { }
pub fn update_positions(&mut self, dt: f64) { }
pub fn update_velocities(&mut self, dt: f64, forces: &[Vec3]) { }

// ✗ AVOID: "set_" unless it's a simple setter
pub fn set_energy(&mut self, e: f64) { }  // Only for simple setters
```

#### Predicates

```rust
// ✓ GOOD: Use is_, has_, needs_ prefixes
pub fn is_ready(&self) -> bool { }
pub fn is_converged(&self) -> bool { }
pub fn has_constraints(&self) -> bool { }
pub fn needs_update(&self) -> bool { }
pub fn should_rebuild_pairlist(&self) -> bool { }
```

#### Getters

```rust
// ✓ GOOD: Bare name (Rust convention)
pub fn temperature(&self) -> f64 { self.temp }
pub fn energy(&self) -> f64 { self.total_energy }
pub fn positions(&self) -> &[Vec3] { &self.pos }

// ✗ AVOID: "get_" prefix (not idiomatic Rust)
pub fn get_temperature(&self) -> f64 { }
pub fn get_energy(&self) -> f64 { }
```

### Types

```rust
// ✓ PascalCase for types
pub struct SystemState { }
pub enum IntegratorType { }
pub trait BoundaryCondition { }

// ✓ Descriptive names
pub struct DistanceRestraint { }
pub struct LeapFrogIntegrator { }
pub struct VelocityVerletIntegrator { }
```

### Constants

```rust
// ✓ SCREAMING_SNAKE_CASE for constants
pub const BOLTZMANN: f64 = 0.00831446;  // kJ/(mol·K)
pub const FOUR_PI_EPS_I: f64 = 138.9354859;  // kJ mol⁻¹ nm e⁻²
pub const MAX_ITERATIONS: usize = 1000;

// ✓ Include units in comments
pub const AVOGADRO: f64 = 6.022140857e23;  // mol⁻¹
```

## Physical Constants

### Constants Module (Recommended)

Create `src/constants.rs`:

```rust
//! Physical constants used in molecular simulations
//!
//! All values in GROMOS units:
//! - Energy: kJ/mol
//! - Distance: nm
//! - Time: ps
//! - Temperature: K
//! - Angle: rad

/// Boltzmann constant (kJ/(mol·K))
pub const BOLTZMANN: f64 = 0.00831446;

/// Gas constant (kJ/(mol·K))
pub const GAS_CONSTANT: f64 = BOLTZMANN;

/// Four pi epsilon inverse (kJ mol⁻¹ nm e⁻²)
/// Used for electrostatic energy: E = FOUR_PI_EPS_I * q1 * q2 / r
pub const FOUR_PI_EPS_I: f64 = 138.9354859;

/// Avogadro's number (mol⁻¹)
pub const AVOGADRO: f64 = 6.022140857e23;

/// Geometric tolerance (nm)
/// Used to avoid division by zero in distance calculations
pub const EPSILON_GEOMETRY: f64 = 1e-10;

/// Convergence tolerance for iterative algorithms
pub const EPSILON_CONVERGENCE: f64 = 1e-6;

/// Strict convergence for energy minimization
pub const EPSILON_STRICT: f64 = 1e-15;

/// Speed of light (nm/ps)
pub const SPEED_OF_LIGHT: f64 = 299792.458;
```

### Usage

```rust
use crate::constants::*;

// ✓ GOOD: Use named constants
let thermal_energy = BOLTZMANN * temperature;
let coulomb_energy = FOUR_PI_EPS_I * q1 * q2 / distance;

if distance2 < EPSILON_GEOMETRY {
    return 0.0;  // Avoid division by zero
}

// ✗ AVOID: Magic numbers
let thermal_energy = 0.00831446 * temperature;  // What is this?
if distance2 < 1e-10 {  // Why 1e-10?
    return 0.0;
}
```

## Error Handling

### Strategy (To Be Implemented)

```rust
use thiserror::Error;

#[derive(Error, Debug)]
pub enum GromosError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("Computation failed: {0}")]
    ComputationFailed(String),

    #[error("Convergence failure after {iterations} iterations")]
    ConvergenceFailed { iterations: usize },

    #[error("File not found: {path}")]
    FileNotFound { path: String },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, GromosError>;
```

### Current Patterns

```rust
// ✓ GOOD: Use Option for missing but valid states
pub fn calculate_parameters(&self, stats: &Statistics) -> Option<(f64, f64)> {
    if !stats.is_ready() {
        return None;
    }
    // ...
}

// ✓ GOOD: Early return for degenerate cases
pub fn calculate_distance_force(&self, r2: f64) -> f64 {
    if r2 < EPSILON_GEOMETRY {
        return 0.0;  // Avoid division by zero
    }
    // ...
}

// ✗ AVOID: Unwrap in library code (only in tests/examples)
let file = std::fs::File::open(path).unwrap();  // Don't do this

// ✓ GOOD: Propagate errors
let file = std::fs::File::open(path)?;  // Better
```

## Documentation

### Module-Level

```rust
//! Pairlist (neighbor list) generation
//!
//! Implementation of Verlet pairlist algorithm for efficient
//! neighbor searching in molecular dynamics simulations.
//!
//! # Algorithm
//!
//! 1. Build pairlist when any atom moves > skin/2
//! 2. Use larger cutoff (cutoff + skin) for pairlist
//! 3. Actual interactions use smaller cutoff
//!
//! # Performance
//!
//! - O(N²) naive algorithm for small systems (N < 1000)
//! - Cell list optimization for large systems
//! - Parallel construction with rayon
//!
//! # References
//!
//! - Verlet, L. (1967). Computer "Experiments" on Classical Fluids.
//!   Physical Review, 159(1), 98-103.
//! - GROMOS manual vol. 2, section 3.2: Pairlist Algorithm
```

### Function Documentation

```rust
/// Calculate Lennard-Jones 12-6 potential energy and force
///
/// Uses the standard LJ potential:
/// $$V(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]$$
///
/// # Arguments
///
/// * `r2` - Squared distance between atoms (nm²)
/// * `c6` - Dispersion coefficient (kJ mol⁻¹ nm⁶)
/// * `c12` - Repulsion coefficient (kJ mol⁻¹ nm¹²)
///
/// # Returns
///
/// Tuple of:
/// - Energy (kJ/mol)
/// - Force magnitude (kJ/mol/nm)
///
/// # Examples
///
/// ```rust
/// use gromos_rs::interaction::nonbonded::calculate_lj;
///
/// let r = 0.35;  // nm
/// let c6 = 1.0;   // kJ mol⁻¹ nm⁶
/// let c12 = 2.0;  // kJ mol⁻¹ nm¹²
///
/// let (energy, force) = calculate_lj(r * r, c6, c12);
/// assert!(energy < 0.0);  // Attractive at this distance
/// ```
///
/// # Performance
///
/// This function is #[inline] for zero-cost abstraction.
/// Hot path in MD simulations (~40% of nonbonded time).
///
/// # Panics
///
/// Panics if `r2` is zero or negative (debug builds only).
///
/// # See Also
///
/// - [`calculate_lj_crf`] - Combined LJ + Coulomb interaction
/// - [`LennardJones`] - LJ parameter container
#[inline]
pub fn calculate_lj(r2: f64, c6: f64, c12: f64) -> (f64, f64) {
    debug_assert!(r2 > 0.0, "r2 must be positive");

    let r6 = r2 * r2 * r2;
    let r12 = r6 * r6;

    let energy = c12 / r12 - c6 / r6;
    let force = 12.0 * c12 / r12 - 6.0 * c6 / r6;

    (energy, force / r2)
}
```

### Inline Comments

```rust
// ✓ GOOD: Explain non-obvious choices
// Use max(0.0, ...) to prevent negative variance from numerical noise
let variance = (variance_sum / count - mean * mean).max(0.0);

// ✓ GOOD: Reference equations from papers
// Andersen thermostat collision probability [Andersen, 1980]
// P = 1 - exp(-dt/τ)
let collision_prob = 1.0 - (-dt / tau).exp();

// ✓ GOOD: Explain algorithm steps
// Step 1: Update velocities to half-step v(t+dt/2)
// Step 2: Update positions to full-step x(t+dt)
// Step 3: Calculate forces F(t+dt)
// Step 4: Update velocities to full-step v(t+dt)

// ✗ AVOID: Stating the obvious
// Set x to 5
let x = 5;

// Loop over atoms
for atom in &atoms {  // This is obvious from the code
    // ...
}
```

## Struct Organization

### Field Ordering

```rust
// Order fields by:
// 1. Most important/frequently accessed
// 2. Logical grouping
// 3. Configuration before state
// 4. Public before private

pub struct MDSimulation {
    // Core configuration (public, most important)
    pub timestep: f64,
    pub temperature: f64,
    pub pressure: Option<f64>,

    // State (public, frequently accessed)
    pub current_step: usize,
    pub total_energy: f64,
    pub kinetic_energy: f64,
    pub potential_energy: f64,

    // Components (private)
    integrator: Box<dyn Integrator>,
    thermostat: Option<Box<dyn Thermostat>>,
    barostat: Option<Box<dyn Barostat>>,

    // Cached/derived (private, least important)
    cached_forces: Vec<Vec3>,
    statistics: Statistics,
}
```

### Builder Pattern

```rust
pub struct DistanceRestraint {
    atom_i: usize,
    atom_j: usize,
    r0: f64,
    k: f64,
    flat_bottom: Option<f64>,
}

impl DistanceRestraint {
    /// Create new distance restraint
    pub fn new(atom_i: usize, atom_j: usize, r0: f64, k: f64) -> Self {
        Self {
            atom_i,
            atom_j,
            r0,
            k,
            flat_bottom: None,
        }
    }

    /// Add flat-bottom potential
    pub fn with_flat_bottom(mut self, half_width: f64) -> Self {
        self.flat_bottom = Some(half_width);
        self
    }
}

// Usage
let restraint = DistanceRestraint::new(0, 1, 0.5, 1000.0)
    .with_flat_bottom(0.05);
```

## Testing

### Test Naming

```rust
#[cfg(test)]
mod tests {
    use super::*;

    // ✓ GOOD: test_ prefix + descriptive name
    #[test]
    fn test_energy_conservation_in_nve() {
        // ...
    }

    #[test]
    fn test_temperature_coupling_converges_exponentially() {
        // ...
    }

    #[test]
    fn test_shake_satisfies_constraints_within_tolerance() {
        // ...
    }

    // ✗ AVOID: No test_ prefix
    #[test]
    fn energy_conservation() {  // Missing prefix
        // ...
    }

    // ✗ AVOID: Unclear names
    #[test]
    fn test_function() {  // Test what?
        // ...
    }
}
```

### Test Structure

```rust
#[test]
fn test_lj_energy_at_minimum() {
    // Arrange: Set up test data
    let epsilon = 0.65;  // kJ/mol
    let sigma = 0.32;    // nm
    let r = sigma * 2.0_f64.powf(1.0/6.0);  // Minimum energy distance

    // Act: Execute function
    let (energy, force) = calculate_lj(r * r, epsilon, sigma);

    // Assert: Check results
    assert_relative_eq!(energy, -epsilon, epsilon = 1e-10);
    assert_relative_eq!(force, 0.0, epsilon = 1e-10);
}
```

## Performance-Critical Code

### Documentation Requirements

```rust
/// Fast innerloop for nonbonded force calculation
///
/// # Performance
///
/// This is expected to be the hottest path in MD simulations.
///
/// **Optimizations applied:**
/// - #[inline(always)] for zero-cost abstraction
/// - SIMD vectorization (AVX2/AVX-512 on x86_64)
/// - Manual loop unrolling for instruction pipelining
/// - Cache-friendly memory layout (SoA vs AoS)
///
/// **Note:** Actual performance has not been measured in gromos-rs.
///
/// **Alternatives:**
/// - [`calculate_nonbonded_simple`] - Simple version for reference
/// - [`calculate_nonbonded_gpu`] - GPU-accelerated version
#[inline(always)]
#[target_feature(enable = "avx2")]
unsafe fn lj_crf_innerloop_simd(/* ... */) {
    // Complex optimized code with extensive comments
}
```

### Justification

```rust
// ✓ GOOD: Explain why optimization is needed
// Hot path: Expected to be called 10^6 times per step in 100k atom system
// Should be optimized as critical performance path in MD simulations
// Theoretical speedup from optimization: 3-4× (based on literature)

// ✗ AVOID: Optimization without explanation
// Fast version
unsafe fn fast_calc(/* ... */) { }
```

## Code Examples

### Before and After Refactoring

#### Before (Inconsistent)
```rust
// Mixed abbreviations, magic numbers, minimal docs
fn compute_e(r2: f64) -> f64 {
    if r2 < 1e-10 { return 0.0; }
    let r6 = r2 * r2 * r2;
    let c12 = 2.0;
    let c6 = 1.0;
    c12 / (r6 * r6) - c6 / r6
}
```

#### After (Following Guide)
```rust
use crate::constants::EPSILON_GEOMETRY;

/// Calculate Lennard-Jones interaction energy
///
/// E_LJ = C₁₂/r¹² - C₆/r⁶
///
/// # Arguments
/// * `r2` - Squared distance (nm²)
/// * `c6` - Dispersion coefficient (kJ mol⁻¹ nm⁶)
/// * `c12` - Repulsion coefficient (kJ mol⁻¹ nm¹²)
///
/// # Returns
/// Energy in kJ/mol
fn calculate_lj_energy(r2: f64, c6: f64, c12: f64) -> f64 {
    // Avoid division by zero
    if r2 < EPSILON_GEOMETRY {
        return 0.0;
    }

    let r6 = r2 * r2 * r2;
    let r12 = r6 * r6;

    c12 / r12 - c6 / r6
}
```

## Migration Plan

### Immediate (Apply to New Code)
- ✓ Use naming conventions from this guide
- ✓ Prefer `calculate_` over `compute_`
- ✓ Use standard abbreviations
- ✓ Document all public APIs

### Short-term (Next Release)
- Create `src/constants.rs` module
- Add module-level docs to all modules
- Unify abbreviation usage across codebase
- Complete test coverage for core modules

### Long-term (Future Releases)
- Implement `GromosError` type
- Refactor old code to match conventions
- Add performance notes to hot paths
- Create architecture documentation

## Enforcement

### Automated (CI/CD)
- ✅ rustfmt enforces formatting
- ✅ clippy catches common issues
- ✅ Naming allowances configured (kT, pH, etc.)
- ✅ Tests must pass

### Manual (Code Review)
- Check naming follows conventions
- Verify documentation completeness
- Ensure error handling strategy
- Validate test coverage

## Related Documentation

- [Contributing Guide](contributing.md) - How to contribute
- [Rust for GROMOS-RS](rust-for-gromos-rs.md) - Rust patterns
- [Tool Development](tool-development-guide.md) - Creating analysis tools
- [NAMING_CONVENTIONS.md](/NAMING_CONVENTIONS.md) - Scientific naming

---

*This guide is based on codebase analysis and documents existing patterns to maintain consistency.*
