# API Reference Overview

GROMOS-RS provides both a library (`lib.rs`) and command-line binaries. This section documents the programmatic API.

## Rust API Documentation

The full Rust API documentation is generated using `rustdoc`:

```bash
# Generate and open documentation
cd gromos-rs
cargo doc --open
```

Documentation will be available at `target/doc/gromos_rs/index.html`.

### Online Documentation

- **docs.rs**: [https://docs.rs/gromos-rs](https://docs.rs/gromos-rs) (when published)
- **Local build**: `cargo doc --no-deps --open`

## Module Overview

### Core Modules

| Module | Description | Key Types |
|--------|-------------|-----------|
| `integrator` | Integration algorithms | `LeapFrog`, `VelocityVerlet`, `StochasticDynamics` |
| `interaction` | Force calculations | `BondedForces`, `NonbondedForces`, `Electrostatics` |
| `algorithm` | Constraints, thermostats, barostats | `Shake`, `Berendsen`, `NoseHoover` |
| `io` | File I/O | `TopologyReader`, `TrajectoryWriter`, `EnergyWriter` |
| `topology` | System topology | `Topology`, `Atom`, `Bond`, `Molecule` |
| `config` | Simulation configuration | `SimulationConfig`, `ForceFieldParams` |

### Advanced Modules

| Module | Description | Key Types |
|--------|-------------|-----------|
| `fep` | Free energy perturbation | `FEPState`, `LambdaSchedule`, `SoftCore` |
| `remd` | Replica exchange MD | `ReplicaExchange`, `ExchangeAttempt` |
| `eds` | Enveloping distribution sampling | `EDSState`, `AdaptiveOffsets` |
| `gamd` | Gaussian accelerated MD | `GaMDBoost`, `SearchMode` |
| `pairlist` | Neighbor lists | `StandardPairlist`, `GridCellPairlist` |

## Quick Examples

### Running a Simulation

```rust
use gromos_rs::{Simulation, Topology, Configuration, Parameters};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load input files
    let topology = Topology::from_file("system.top")?;
    let config = Configuration::from_file("initial.cnf")?;
    let params = Parameters::from_file("md.imd")?;

    // Create simulation
    let mut sim = Simulation::new(topology, config, params)?;

    // Run MD
    for step in 0..params.nsteps {
        sim.step()?;

        if step % 1000 == 0 {
            sim.write_trajectory()?;
            sim.write_energy()?;
        }
    }

    Ok(())
}
```

### Custom Force Calculation

```rust
use gromos_rs::interaction::{BondedForces, NonbondedForces};
use gromos_rs::topology::Topology;

fn calculate_forces(topology: &Topology, positions: &[f64]) -> Vec<f64> {
    let mut forces = vec![0.0; positions.len()];

    // Bonded forces
    BondedForces::calculate(topology, positions, &mut forces);

    // Nonbonded forces
    NonbondedForces::calculate(topology, positions, &mut forces);

    forces
}
```

### Free Energy Perturbation

```rust
use gromos_rs::fep::{FEPState, SoftCore};

fn run_fep_window(lambda: f64) -> Result<f64, Box<dyn std::error::Error>> {
    let mut sim = setup_simulation()?;

    // Set up FEP state
    let fep = FEPState::new(lambda)?;
    sim.set_fep_state(fep);

    // Run simulation
    let mut dhdl = 0.0;
    for step in 0..sim.params.nsteps {
        sim.step()?;
        dhdl += sim.calculate_dhdl()?;
    }

    Ok(dhdl / sim.params.nsteps as f64)
}
```

### Replica Exchange

```rust
use gromos_rs::remd::ReplicaExchange;

fn run_remd() -> Result<(), Box<dyn std::error::Error>> {
    // Set up replicas
    let temperatures = vec![300.0, 310.0, 320.0, 330.0, 340.0, 350.0];
    let mut remd = ReplicaExchange::new(temperatures)?;

    // Run REMD
    for cycle in 0..1000 {
        // MD steps for each replica
        remd.run_md_steps(1000)?;

        // Attempt exchanges
        remd.attempt_exchanges()?;

        // Write statistics
        remd.write_stats()?;
    }

    Ok(())
}
```

## API Design Principles

### 1. Type Safety

GROMOS-RS uses Rust's type system to prevent errors at compile time:

```rust
// Units are explicit
struct Temperature(f64);  // Kelvin
struct Pressure(f64);     // bar
struct Energy(f64);       // kJ/mol

// Can't mix incompatible types
let temp = Temperature(300.0);
let press = Pressure(1.0);
// let x = temp + press;  // Compile error!
```

### 2. Zero-Cost Abstractions

High-level APIs have no runtime overhead:

```rust
// This:
forces.par_iter_mut()
    .zip(positions.par_iter())
    .for_each(|(f, pos)| *f = calculate_force(pos));

// Compiles to the same machine code as:
for i in 0..forces.len() {
    forces[i] = calculate_force(&positions[i]);
}
```

### 3. Iterator-based Processing

```rust
// Trajectory analysis using iterators
let avg_energy = trajectory
    .frames()
    .map(|frame| frame.potential_energy())
    .sum::<f64>() / trajectory.len() as f64;
```

### 4. Error Handling

All fallible operations return `Result`:

```rust
use gromos_rs::error::GromosPSError;

fn load_system(path: &str) -> Result<Topology, GromosPSError> {
    let topology = Topology::from_file(path)?;  // ? propagates errors
    validate_topology(&topology)?;
    Ok(topology)
}
```

## Performance Considerations

### SIMD Vectorization

Many operations are automatically vectorized:

```rust
// Automatically uses AVX2/AVX-512
#[target_feature(enable = "avx2")]
unsafe fn lennard_jones_simd(r2: f32x8, c6: f32x8, c12: f32x8) -> f32x8 {
    let r6 = r2 * r2 * r2;
    let r12 = r6 * r6;
    c12 / r12 - c6 / r6
}
```

### Parallelization

Use Rayon for automatic parallelization:

```rust
use rayon::prelude::*;

// Parallel pairlist generation
let pairs: Vec<(usize, usize)> = atoms
    .par_iter()
    .enumerate()
    .flat_map(|(i, atom_i)| {
        atoms[i+1..]
            .iter()
            .enumerate()
            .filter_map(move |(j, atom_j)| {
                if distance(atom_i, atom_j) < cutoff {
                    Some((i, i + j + 1))
                } else {
                    None
                }
            })
    })
    .collect();
```

## Extending GROMOS-RS

### Custom Integrator

```rust
use gromos_rs::integrator::Integrator;

struct MyIntegrator {
    dt: f64,
}

impl Integrator for MyIntegrator {
    fn step(&mut self, system: &mut System) -> Result<(), Error> {
        // Update positions
        for (pos, vel) in system.positions.iter_mut()
            .zip(system.velocities.iter())
        {
            *pos += vel * self.dt;
        }

        // Calculate forces
        system.calculate_forces()?;

        // Update velocities
        for (vel, force, mass) in izip!(
            system.velocities.iter_mut(),
            system.forces.iter(),
            system.masses.iter()
        ) {
            *vel += force / mass * self.dt;
        }

        Ok(())
    }
}
```

### Custom Force Term

```rust
use gromos_rs::interaction::Force;

struct MyForce {
    k: f64,  // Force constant
}

impl Force for MyForce {
    fn calculate(&self, system: &System, forces: &mut [f64]) -> f64 {
        let mut energy = 0.0;

        // Custom force calculation
        for (i, pos) in system.positions.iter().enumerate() {
            let f = -self.k * pos;  // Harmonic restraint to origin
            forces[i * 3 + 0] += f * pos[0];
            forces[i * 3 + 1] += f * pos[1];
            forces[i * 3 + 2] += f * pos[2];

            energy += 0.5 * self.k * (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        }

        energy
    }
}
```

## Testing

GROMOS-RS has extensive test coverage:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leap_frog() {
        let mut sim = Simulation::new_test_system();
        let initial_energy = sim.total_energy();

        // Run NVE simulation
        for _ in 0..10000 {
            sim.leap_frog_step();
        }

        let final_energy = sim.total_energy();
        let drift = (final_energy - initial_energy).abs() / initial_energy;

        assert!(drift < 1e-4, "Energy drift too large: {}", drift);
    }
}
```

## Documentation Comments

All public APIs include documentation:

```rust
/// Calculates Lennard-Jones and Coulomb interactions.
///
/// # Arguments
///
/// * `positions` - Atomic positions [x,y,z,...] in nm
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
pub fn calculate(
    positions: &[f64],
    charges: &[f64],
    pairlist: &[(usize, usize)]
) -> (f64, Vec<f64>) {
    // ...
}
```

## Further Reading

- [Core Modules](core.md) - Detailed core API docs
- [Integration](integration.md) - Integration algorithms
- [Interactions](interactions.md) - Force field calculations
- [I/O](io.md) - File input/output

For complete API documentation, run:

```bash
cargo doc --open
```
