# Rust for GROMOS-RS: A Developer's Guide

This guide introduces Rust concepts through the lens of GROMOS-RS development, helping both new and experienced developers understand how Rust is used in this molecular dynamics codebase.

## Why Rust for MD Simulations?

### The Problem with C++ in MD

Traditional MD codes (like md++) face challenges:

```cpp
// C++: Manual memory management, potential bugs
void calculate_forces(std::vector<Atom> &atoms,
                     std::vector<double> &forces) {
    // Potential issues:
    // 1. Data race if parallelized incorrectly
    // 2. Buffer overflow if indices wrong
    // 3. Use-after-free if atoms reallocated
    // 4. Iterator invalidation

    #pragma omp parallel for
    for (int i = 0; i < atoms.size(); i++) {
        forces[i] = 0.0;  // Race condition!

        for (int j = i+1; j < atoms.size(); j++) {
            double f = calc_pairwise_force(atoms[i], atoms[j]);
            forces[i] += f;   // Data race!
            forces[j] -= f;   // Data race!
        }
    }
}
```

### The Rust Solution

```rust
// Rust: Safety guaranteed at compile time
fn calculate_forces(atoms: &[Atom], forces: &mut [f64]) {
    // Compiler checks:
    // ✓ No data races (checked at compile time!)
    // ✓ No buffer overflows (bounds checking)
    // ✓ No use-after-free (ownership system)
    // ✓ No null pointers (Option<T> instead)

    forces.par_iter_mut()  // Parallel iterator
        .zip(atoms.par_iter())
        .enumerate()
        .for_each(|(i, (force, atom_i))| {
            *force = 0.0;

            for (j, atom_j) in atoms[i+1..].iter().enumerate() {
                let f = calc_pairwise_force(atom_i, atom_j);
                // This is SAFE! No data races possible
                // Compiler guarantees it!
            }
        });
}
```

**Key Benefits**:
1. **Memory safety** without garbage collection
2. **Thread safety** guaranteed by compiler
3. **Zero-cost abstractions** (no runtime overhead)
4. **Performance** comparable to C, often faster

## Rust Fundamentals for GROMOS-RS

### 1. Ownership and Borrowing

The foundation of Rust's safety guarantees.

#### Ownership Rules

```rust
// Rule 1: Each value has exactly one owner
let atoms = vec![Atom::new(), Atom::new()];  // atoms owns the Vec

// Rule 2: When owner goes out of scope, value is dropped
{
    let temp = vec![1.0, 2.0, 3.0];
}  // temp is dropped here, memory freed automatically

// Rule 3: You can transfer ownership
let atoms2 = atoms;  // Ownership moved to atoms2
// atoms is now invalid! Compile error if used

// In GROMOS-RS:
let mut sim = Simulation::new(topology, config, params);
run_simulation(sim);  // sim moved into function
// sim is now invalid! Cannot use it again
```

#### Borrowing in Force Calculations

```rust
// Real example from gromos-rs/src/interaction/nonbonded.rs

pub fn calculate_lj_crf(
    positions: &[Vector3<f64>],      // Immutable borrow
    charges: &[f64],                  // Immutable borrow
    pairlist: &[(usize, usize)],     // Immutable borrow
    lj_params: &[[f64; 2]],          // Immutable borrow
    forces: &mut [Vector3<f64>],      // Mutable borrow (exclusive!)
    energy: &mut f64,                 // Mutable borrow
) {
    // Can have:
    // - Any number of immutable borrows (&T)
    // - OR exactly one mutable borrow (&mut T)
    // - But NOT both at the same time!

    for &(i, j) in pairlist {
        let r_ij = positions[j] - positions[i];
        let r2 = r_ij.norm_squared();

        // Calculate force
        let f = lj_force(r2, lj_params[i][j]);

        // Safe: Compiler guarantees no two threads access same index
        forces[i] += f * r_ij;
        forces[j] -= f * r_ij;

        *energy += lj_energy(r2, lj_params[i][j]);
    }
}
```

**Why this matters**:
- Prevents data races at compile time
- No need for locks/mutexes (in most cases)
- Compiler optimizes better (knows about aliasing)

### 2. Pattern Matching and Error Handling

#### Result Type for Robust I/O

```rust
// gromos-rs/src/io/topology.rs

use std::fs::File;
use std::io::{self, BufRead, BufReader};

pub fn read_topology(path: &str) -> Result<Topology, IoError> {
    // Result<T, E> is Rust's way of handling errors
    // Either: Ok(value) or Err(error)

    let file = File::open(path)?;  // ? propagates errors
    // Equivalent to:
    // let file = match File::open(path) {
    //     Ok(f) => f,
    //     Err(e) => return Err(e.into()),
    // };

    let reader = BufReader::new(file);
    let mut topology = Topology::default();

    for line in reader.lines() {
        let line = line?;  // Handle I/O errors

        match parse_topology_line(&line) {
            Ok(TopologyEntry::Atom(atom)) => {
                topology.add_atom(atom);
            }
            Ok(TopologyEntry::Bond(bond)) => {
                topology.add_bond(bond);
            }
            Ok(TopologyEntry::Comment) => {
                // Skip comments
            }
            Err(e) => {
                return Err(IoError::ParseError {
                    line: line.clone(),
                    error: e,
                });
            }
        }
    }

    Ok(topology)
}

// Usage:
let topology = read_topology("system.top")?;
// If this fails, error propagates up
// No unchecked exceptions!
```

#### Pattern Matching in Integrators

```rust
// gromos-rs/src/integrator.rs

pub enum IntegratorType {
    LeapFrog { dt: f64 },
    VelocityVerlet { dt: f64 },
    StochasticDynamics { dt: f64, gamma: f64 },
    SteepestDescent { step_size: f64, tolerance: f64 },
}

impl IntegratorType {
    pub fn step(&self, system: &mut System) -> Result<(), Error> {
        match self {
            IntegratorType::LeapFrog { dt } => {
                // v(t+dt/2) = v(t-dt/2) + a(t)*dt
                for (vel, force, mass) in izip!(
                    system.velocities.iter_mut(),
                    &system.forces,
                    &system.masses
                ) {
                    *vel += force / mass * dt;
                }

                // r(t+dt) = r(t) + v(t+dt/2)*dt
                for (pos, vel) in system.positions.iter_mut()
                    .zip(&system.velocities)
                {
                    *pos += vel * dt;
                }

                Ok(())
            }

            IntegratorType::VelocityVerlet { dt } => {
                // r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt²
                // v(t+dt) = v(t) + 0.5*[a(t) + a(t+dt)]*dt

                let old_forces = system.forces.clone();

                // Update positions
                for (pos, vel, force, mass) in izip!(
                    system.positions.iter_mut(),
                    &system.velocities,
                    &old_forces,
                    &system.masses
                ) {
                    *pos += vel * dt + 0.5 * force / mass * dt * dt;
                }

                // Calculate new forces
                system.calculate_forces()?;

                // Update velocities (half old, half new forces)
                for (vel, f_old, f_new, mass) in izip!(
                    system.velocities.iter_mut(),
                    &old_forces,
                    &system.forces,
                    &system.masses
                ) {
                    *vel += 0.5 * (f_old + f_new) / mass * dt;
                }

                Ok(())
            }

            IntegratorType::StochasticDynamics { dt, gamma } => {
                // Langevin dynamics
                // v(t+dt) = v(t) + [F - γv + R]*dt/m
                // where R ~ N(0, 2γkT/m)

                use rand_distr::{Normal, Distribution};
                let mut rng = rand::thread_rng();

                for (vel, force, mass) in izip!(
                    system.velocities.iter_mut(),
                    &system.forces,
                    &system.masses
                ) {
                    // Friction term
                    let friction = -gamma * *vel;

                    // Random force (thermal fluctuations)
                    let sigma = (2.0 * gamma * K_B * system.temperature / mass).sqrt();
                    let normal = Normal::new(0.0, sigma).unwrap();
                    let random_force = normal.sample(&mut rng);

                    // Update velocity
                    *vel += (force / mass + friction + random_force) * dt;
                }

                Ok(())
            }

            IntegratorType::SteepestDescent { step_size, tolerance } => {
                // r_new = r_old - α * ∇E
                let max_force = system.forces.iter()
                    .map(|f| f.norm())
                    .fold(0.0, f64::max);

                if max_force < *tolerance {
                    return Ok(()); // Converged!
                }

                for (pos, force) in system.positions.iter_mut()
                    .zip(&system.forces)
                {
                    *pos -= step_size * force;
                }

                Ok(())
            }
        }
    }
}
```

**Key Patterns**:
- `match` is exhaustive (compiler checks all cases)
- No fall-through bugs
- Easy to add new integrators
- Type-safe dispatch

### 3. Traits: Rust's Interface System

#### Force Trait

```rust
// gromos-rs/src/interaction/mod.rs

pub trait Force {
    /// Calculate force and energy for this interaction
    fn calculate(
        &self,
        system: &System,
        forces: &mut [Vector3<f64>],
    ) -> f64;  // Returns energy

    /// Optional: calculate force and energy derivatives (for FEP)
    fn calculate_perturbed(
        &self,
        system: &System,
        lambda: f64,
        forces: &mut [Vector3<f64>],
        dhdl: &mut f64,  // dH/dλ for TI
    ) -> f64 {
        // Default implementation: non-perturbed
        self.calculate(system, forces)
    }
}

// Implement for bonded interactions
impl Force for BondedForces {
    fn calculate(
        &self,
        system: &System,
        forces: &mut [Vector3<f64>],
    ) -> f64 {
        let mut energy = 0.0;

        // Bonds
        for bond in &self.bonds {
            let i = bond.atom_i;
            let j = bond.atom_j;

            let r_ij = system.positions[j] - system.positions[i];
            let r = r_ij.norm();

            // Quartic bond: V = 1/4 * k * (r² - r0²)²
            let delta_r2 = r*r - bond.r0*bond.r0;
            let e = 0.25 * bond.k * delta_r2 * delta_r2;
            let f_scalar = -bond.k * delta_r2 * r;  // -dV/dr

            let f_vec = (f_scalar / r) * r_ij;
            forces[i] -= f_vec;
            forces[j] += f_vec;

            energy += e;
        }

        // Angles
        for angle in &self.angles {
            // ... similar
        }

        // Dihedrals
        for dihedral in &self.dihedrals {
            // ... similar
        }

        energy
    }
}

// Implement for nonbonded
impl Force for NonbondedForces {
    fn calculate(
        &self,
        system: &System,
        forces: &mut [Vector3<f64>],
    ) -> f64 {
        calculate_lj_crf(
            &system.positions,
            &system.charges,
            &self.pairlist,
            &self.lj_params,
            forces,
            &mut energy,
        );

        energy
    }
}

// Now we can use polymorphism!
fn total_energy(forces_list: &[Box<dyn Force>], system: &System) -> f64 {
    let mut total_forces = vec![Vector3::zero(); system.n_atoms];
    let mut total_energy = 0.0;

    for force in forces_list {
        total_energy += force.calculate(system, &mut total_forces);
    }

    total_energy
}
```

**Trait Benefits**:
- Zero-cost abstraction (static dispatch when possible)
- Easy to extend (new force types)
- Testable (mock implementations)

### 4. Parallel Programming with Rayon

Rayon makes parallelization trivial and **safe**.

#### Parallel Force Calculation

```rust
// gromos-rs/src/interaction/nonbonded.rs

use rayon::prelude::*;

pub fn calculate_nonbonded_parallel(
    positions: &[Vector3<f64>],
    charges: &[f64],
    pairlist: &[(usize, usize)],
    lj_params: &[[f64; 2]],
    forces: &mut [Vector3<f64>],
) -> f64 {
    // Problem: Multiple threads will write to same force array
    // Solution: Each thread gets its own force array, then reduce

    use std::sync::atomic::{AtomicU64, Ordering};

    // Atomic for energy accumulation
    let energy = AtomicU64::new(0);

    // Parallel iteration over pairlist
    let local_forces: Vec<_> = pairlist
        .par_chunks(1000)  // Process in chunks for cache efficiency
        .map(|chunk| {
            // Each thread gets its own force array
            let mut thread_forces = vec![Vector3::zero(); positions.len()];
            let mut thread_energy = 0.0;

            for &(i, j) in chunk {
                let r_ij = positions[j] - positions[i];
                let r2 = r_ij.norm_squared();

                if r2 < CUTOFF_SQUARED {
                    // LJ interaction
                    let r2_inv = 1.0 / r2;
                    let r6_inv = r2_inv * r2_inv * r2_inv;
                    let c6 = lj_params[i][0];
                    let c12 = lj_params[i][1];

                    let e_lj = c12 * r6_inv * r6_inv - c6 * r6_inv;
                    let f_lj = (12.0*c12*r6_inv*r6_inv - 6.0*c6*r6_inv) * r2_inv;

                    // Coulomb with reaction field
                    let q_ij = charges[i] * charges[j];
                    let e_crf = q_ij * (1.0/r2.sqrt() + CRF_CONST * r2);
                    let f_crf = q_ij * (1.0/r2.sqrt().powi(3) - 2.0*CRF_CONST);

                    // Total
                    let f_total = (f_lj + f_crf) * r_ij;

                    thread_forces[i] -= f_total;
                    thread_forces[j] += f_total;

                    thread_energy += e_lj + e_crf;
                }
            }

            (thread_forces, thread_energy)
        })
        .reduce(
            || (vec![Vector3::zero(); positions.len()], 0.0),
            |(mut f1, e1), (f2, e2)| {
                // Combine results from different threads
                for (force1, force2) in f1.iter_mut().zip(&f2) {
                    *force1 += force2;
                }
                (f1, e1 + e2)
            }
        );

    // Copy to output
    forces.copy_from_slice(&local_forces.0);

    local_forces.1  // Return total energy
}
```

**Rayon Magic**:
- `par_iter()` - Parallel iterator (automatic work-stealing)
- `par_chunks()` - Process in parallel chunks
- `map()` + `reduce()` - Map-reduce pattern
- **Zero data races** - Compiler guarantees safety!

#### SIMD Vectorization

```rust
// gromos-rs/src/math/simd.rs

#[cfg(target_feature = "avx2")]
use std::arch::x86_64::*;

#[target_feature(enable = "avx2")]
pub unsafe fn lj_crf_simd(
    r2: __m256,  // 8 x f32 (8 distances squared)
    c6: __m256,  // 8 x f32 (C6 parameters)
    c12: __m256, // 8 x f32 (C12 parameters)
    q_ij: __m256, // 8 x f32 (charge products)
) -> (__m256, __m256) {  // (energy, force)
    // Compute 8 LJ+CRF interactions at once!

    // r⁻²
    let r2_inv = _mm256_div_ps(_mm256_set1_ps(1.0), r2);

    // r⁻⁶ = (r⁻²)³
    let r6_inv = _mm256_mul_ps(
        _mm256_mul_ps(r2_inv, r2_inv),
        r2_inv
    );

    // r⁻¹² = (r⁻⁶)²
    let r12_inv = _mm256_mul_ps(r6_inv, r6_inv);

    // E_LJ = C₁₂/r¹² - C₆/r⁶
    let e_lj = _mm256_sub_ps(
        _mm256_mul_ps(c12, r12_inv),
        _mm256_mul_ps(c6, r6_inv)
    );

    // F_LJ = 12*C₁₂/r¹⁴ - 6*C₆/r⁸
    let f_lj = _mm256_mul_ps(
        _mm256_sub_ps(
            _mm256_mul_ps(_mm256_set1_ps(12.0), _mm256_mul_ps(c12, r12_inv)),
            _mm256_mul_ps(_mm256_set1_ps(6.0), _mm256_mul_ps(c6, r6_inv))
        ),
        r2_inv
    );

    // Coulomb (simplified)
    let r_inv = _mm256_sqrt_ps(r2_inv);
    let e_crf = _mm256_mul_ps(q_ij, r_inv);
    let f_crf = _mm256_mul_ps(e_crf, r2_inv);

    // Total
    let e_total = _mm256_add_ps(e_lj, e_crf);
    let f_total = _mm256_add_ps(f_lj, f_crf);

    (e_total, f_total)
}
```

**Performance**:
- 8 interactions per CPU cycle (AVX2)
- 16 interactions per CPU cycle (AVX-512)
- 2-3x speedup over scalar code

### 5. Type System for Physical Units

Prevent unit errors at compile time!

```rust
// gromos-rs/src/units.rs

use std::ops::{Add, Sub, Mul, Div};

/// Newtype pattern for type-safe units
#[derive(Debug, Clone, Copy)]
pub struct Length(f64);  // nm

#[derive(Debug, Clone, Copy)]
pub struct Energy(f64);  // kJ/mol

#[derive(Debug, Clone, Copy)]
pub struct Temperature(f64);  // K

#[derive(Debug, Clone, Copy)]
pub struct Pressure(f64);  // bar

impl Length {
    pub fn nm(value: f64) -> Self { Length(value) }
    pub fn angstrom(value: f64) -> Self { Length(value / 10.0) }
    pub fn as_nm(&self) -> f64 { self.0 }
}

impl Energy {
    pub fn kj_mol(value: f64) -> Self { Energy(value) }
    pub fn kcal_mol(value: f64) -> Self { Energy(value * 4.184) }
    pub fn as_kj_mol(&self) -> f64 { self.0 }
}

// Arithmetic with correct units
impl Add for Length {
    type Output = Length;
    fn add(self, other: Length) -> Length {
        Length(self.0 + other.0)
    }
}

impl Mul<f64> for Length {
    type Output = Length;
    fn mul(self, scalar: f64) -> Length {
        Length(self.0 * scalar)
    }
}

impl Div for Length {
    type Output = f64;  // Dimensionless ratio
    fn div(self, other: Length) -> f64 {
        self.0 / other.0
    }
}

// Prevent unit errors at compile time!
fn calculate_density(mass: f64, volume: Length) -> f64 {
    // mass / volume.0  // ← Would work but bypasses type system

    // Better: define volume type
    // mass / volume  // ← Won't compile! Need Volume type

    // For now:
    mass / volume.as_nm().powi(3)
}

// Usage:
let cutoff = Length::nm(1.4);
let box_size = Length::angstrom(45.0);

if cutoff > box_size / 2.0 {
    panic!("Cutoff too large for box!");
}

let temp = Temperature(300.0);
let energy = Energy::kcal_mol(150.0);

// let x = temp + energy;  // ← Won't compile! Type mismatch
```

## GROMOS-RS Code Architecture

### Directory Structure

```
gromos-rs/src/
├── lib.rs                  # Library entry point
├── main.rs                 # Binary entry points
│
├── # === CORE SYSTEM ===
├── topology.rs             # Topology (atoms, bonds, molecules)
├── configuration.rs        # Positions, velocities, box
├── simulation.rs           # Main simulation loop
├── parameters.rs           # Force field parameters
│
├── # === INTEGRATION ===
├── integrator.rs           # Integration algorithms
│   # - LeapFrog, VelocityVerlet
│   # - StochasticDynamics
│   # - SteepestDescent
│
├── # === FORCES ===
├── interaction/
│   ├── mod.rs              # Force trait
│   ├── bonded.rs           # Bonds, angles, dihedrals
│   ├── nonbonded.rs        # LJ + Coulomb
│   ├── electrostatics.rs   # RF, PME, Ewald
│   └── restraints.rs       # Position, distance restraints
│
├── # === ALGORITHMS ===
├── algorithm/
│   ├── constraints.rs      # SHAKE, SETTLE, LINCS
│   ├── thermostats.rs      # Berendsen, Nosé-Hoover
│   └── barostats.rs        # Berendsen, Parrinello-Rahman
│
├── # === ADVANCED SAMPLING ===
├── fep.rs                  # Free energy perturbation
├── remd.rs                 # Replica exchange MD
├── eds.rs                  # Enveloping distribution sampling
├── gamd.rs                 # Gaussian accelerated MD
│
├── # === I/O ===
├── io/
│   ├── topology.rs         # .top reader
│   ├── coordinate.rs       # .cnf reader/writer
│   ├── imd.rs              # .imd parameter reader
│   ├── trajectory.rs       # .trc writer
│   ├── energy.rs           # .tre writer
│   ├── force.rs            # .trf writer
│   ├── ptp.rs              # .ptp writer (FEP)
│   └── dlg.rs              # .dlg writer (dH/dλ)
│
├── # === UTILITIES ===
├── pairlist.rs             # Neighbor lists
├── math.rs                 # Vector math, PBC
├── units.rs                # Physical units
└── error.rs                # Error types
```

### Typical Development Workflow

```bash
# 1. Create feature branch
git checkout -b feature/my-new-feature

# 2. Write code (with tests!)
# gromos-rs/src/interaction/my_new_force.rs

# 3. Run tests frequently
cargo test

# 4. Run clippy (Rust linter)
cargo clippy -- -D warnings

# 5. Format code
cargo fmt

# 6. Run benchmarks
cargo bench

# 7. Commit
git commit -m "Add my new feature"

# 8. Push and create PR
git push origin feature/my-new-feature
```

## Best Practices for GROMOS-RS

### 1. Use Clippy

```bash
# Clippy is Rust's linter
cargo clippy -- -D warnings

# Common issues it catches:
# - Unnecessary clones
# - Inefficient iterations
# - Potential bugs
# - Style inconsistencies
```

### 2. Write Tests

```rust
// Every module should have tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leap_frog_energy_conservation() {
        // Create simple test system
        let mut system = create_harmonic_oscillator();
        let integrator = LeapFrog { dt: 0.001 };

        let e0 = system.total_energy();

        // Run 10000 steps
        for _ in 0..10000 {
            integrator.step(&mut system).unwrap();
        }

        let e1 = system.total_energy();
        let drift = (e1 - e0).abs() / e0;

        assert!(drift < 1e-4, "Energy drift: {}", drift);
    }

    #[test]
    fn test_shake_satisfies_constraints() {
        let mut system = create_water_molecule();
        let shake = Shake::new(1e-4, 1000);

        shake.apply(&mut system).unwrap();

        // Check bond lengths
        let r_oh1 = (system.positions[0] - system.positions[1]).norm();
        let r_oh2 = (system.positions[0] - system.positions[2]).norm();

        assert!((r_oh1 - 0.1).abs() < 1e-6);
        assert!((r_oh2 - 0.1).abs() < 1e-6);
    }
}
```

### 3. Use Type System

```rust
// Bad: Primitive obsession
fn calculate_energy(positions: &[f64]) -> f64 { /* ... */ }

// Good: Use types
fn calculate_energy(positions: &[Vector3<f64>]) -> Energy { /* ... */ }

// Better: Use newtypes
fn calculate_energy(positions: &[Position]) -> Energy { /* ... */ }
```

### 4. Avoid Premature Optimization

```rust
// First: Write clear, correct code
fn slow_but_clear(atoms: &[Atom]) -> f64 {
    atoms.iter()
        .map(|atom| atom.kinetic_energy())
        .sum()
}

// Then: Profile and optimize hot paths
fn fast_optimized(atoms: &[Atom]) -> f64 {
    atoms.par_iter()  // Parallel
        .map(|atom| atom.kinetic_energy())
        .sum()
}

// Use: cargo bench to measure
```

### 5. Document Public APIs

```rust
/// Calculate Lennard-Jones and Coulomb interactions.
///
/// # Arguments
///
/// * `positions` - Atomic positions in nm
/// * `charges` - Atomic charges in e
/// * `pairlist` - List of interacting pairs
///
/// # Returns
///
/// Tuple of (potential energy in kJ/mol, forces in kJ/mol/nm)
///
/// # Examples
///
/// ```
/// use gromos_rs::interaction::nonbonded;
///
/// let (energy, forces) = nonbonded::calculate(&positions, &charges, &pairlist);
/// ```
///
/// # Performance
///
/// This function is SIMD-vectorized and parallelized with Rayon.
/// Expected performance: ~2-3x faster than scalar C++ code.
pub fn calculate(
    positions: &[Vector3<f64>],
    charges: &[f64],
    pairlist: &[(usize, usize)],
) -> (f64, Vec<Vector3<f64>>) {
    // ...
}
```

## Common Patterns in GROMOS-RS

### Builder Pattern

```rust
pub struct Simulation {
    topology: Topology,
    configuration: Configuration,
    parameters: Parameters,
    integrator: Box<dyn Integrator>,
    forces: Vec<Box<dyn Force>>,
}

impl Simulation {
    pub fn builder() -> SimulationBuilder {
        SimulationBuilder::default()
    }
}

pub struct SimulationBuilder {
    topology: Option<Topology>,
    configuration: Option<Configuration>,
    parameters: Option<Parameters>,
    integrator: Option<Box<dyn Integrator>>,
}

impl SimulationBuilder {
    pub fn topology(mut self, topo: Topology) -> Self {
        self.topology = Some(topo);
        self
    }

    pub fn configuration(mut self, conf: Configuration) -> Self {
        self.configuration = Some(conf);
        self
    }

    pub fn build(self) -> Result<Simulation, Error> {
        Ok(Simulation {
            topology: self.topology.ok_or(Error::MissingTopology)?,
            configuration: self.configuration.ok_or(Error::MissingConfiguration)?,
            parameters: self.parameters.unwrap_or_default(),
            integrator: self.integrator.ok_or(Error::MissingIntegrator)?,
            forces: vec![],
        })
    }
}

// Usage:
let sim = Simulation::builder()
    .topology(topology)
    .configuration(config)
    .integrator(Box::new(LeapFrog { dt: 0.002 }))
    .build()?;
```

### Iterator Pattern

```rust
// Trajectory iterator
pub struct TrajectoryIterator {
    file: BufReader<File>,
    current_frame: usize,
}

impl Iterator for TrajectoryIterator {
    type Item = Result<Frame, IoError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_frame() {
            Ok(frame) => {
                self.current_frame += 1;
                Some(Ok(frame))
            }
            Err(IoError::Eof) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

// Usage:
for frame in trajectory.frames() {
    let frame = frame?;
    let rmsd = calculate_rmsd(&frame, &reference);
    println!("{}", rmsd);
}
```

## Performance Tips

### 1. Use Release Builds

```bash
# Debug (slow): cargo build
# Release (fast): cargo build --release

# 10-100x difference!
```

### 2. Profile with Criterion

```rust
// benches/nonbonded.rs
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_nonbonded(c: &mut Criterion) {
    let positions = create_test_system(10000);
    let pairlist = generate_pairlist(&positions);

    c.bench_function("nonbonded 10k atoms", |b| {
        b.iter(|| {
            calculate_nonbonded(
                black_box(&positions),
                black_box(&pairlist),
            )
        })
    });
}

criterion_group!(benches, bench_nonbonded);
criterion_main!(benches);
```

### 3. Use `inline` Judiciously

```rust
#[inline]
pub fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

#[inline(never)]  // Don't inline (for profiling)
pub fn large_function() { /* ... */ }
```

### 4. Avoid Allocations in Hot Loops

```rust
// Bad: Allocates every iteration
for _ in 0..1000000 {
    let v = vec![0.0; 100];  // Allocation!
    process(&v);
}

// Good: Allocate once, reuse
let mut v = vec![0.0; 100];
for _ in 0..1000000 {
    v.fill(0.0);  // Clear
    process(&mut v);
}
```

## Resources

### Learning Rust

1. **The Rust Book**: https://doc.rust-lang.org/book/
2. **Rust by Example**: https://doc.rust-lang.org/rust-by-example/
3. **Rustlings**: https://github.com/rust-lang/rustlings

### Advanced Topics

1. **Async Book**: https://rust-lang.github.io/async-book/
2. **Rustonomicon** (unsafe): https://doc.rust-lang.org/nomicon/
3. **Performance Book**: https://nnethercote.github.io/perf-book/

### Ecosystem

1. **Rayon** (parallelism): https://docs.rs/rayon/
2. **ndarray** (arrays): https://docs.rs/ndarray/
3. **criterion** (benchmarking): https://docs.rs/criterion/

---

**Next**: [Contributing Guide](contributing.md) | [Code Architecture](architecture.md) | [Testing Guide](testing.md)
