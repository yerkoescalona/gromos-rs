# Tool Development Guide: Writing GROMOS++ Style Programs in Rust

This guide teaches you how to write analysis and preprocessing tools for gromos-rs, following the philosophy and patterns of the 104 gromos++ programs.

## Table of Contents

1. [GROMOS++ Philosophy](#gromos-philosophy)
2. [Basic Tool Structure](#basic-tool-structure)
3. [Tutorial: Your First Tool](#tutorial-your-first-tool)
4. [Common Patterns](#common-patterns)
5. [Advanced Examples](#advanced-examples)
6. [Best Practices](#best-practices)
7. [Testing Your Tool](#testing-your-tool)
8. [Publishing Your Tool](#publishing-your-tool)

---

## GROMOS++ Philosophy

### Unix Philosophy

GROMOS++ follows the Unix philosophy:

```
1. Do one thing well
2. Work together via standard formats
3. Text-based for composability
4. Scriptable and automatable
```

Example workflow:
```bash
# Each tool does ONE thing
ene_ana @traj ener.tre @prop e_pot > energy.dat
rmsd @traj traj.trc @ref crystal.pdb > rmsd.dat
hbond @traj traj.trc > hbonds.dat

# Combine with Unix tools
paste energy.dat rmsd.dat | awk '{print $1, $2, $4}' > combined.dat
```

### Tool Design Principles

```
✓ Single Responsibility: One tool, one task
✓ Standard I/O: Read trajectories, write results
✓ Composable: Output can be input to other tools
✓ Documented: Clear help text
✓ Predictable: Same inputs → same outputs
✓ Fast: Optimized for performance
```

---

## Basic Tool Structure

### Minimal Tool Template

```rust
// gromos-rs/src/bin/my_tool.rs

use gromos_rs::{Topology, TrajectoryReader};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "my_tool")]
#[command(about = "Brief description of what this tool does", long_about = None)]
struct Args {
    /// Input topology file (.top)
    #[arg(short, long)]
    topology: String,

    /// Input trajectory file (.trc)
    #[arg(short = 'x', long)]
    trajectory: String,

    /// Output file (default: stdout)
    #[arg(short, long)]
    output: Option<String>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments
    let args = Args::parse();

    // Load topology
    let topology = Topology::from_file(&args.topology)?;

    // Open trajectory
    let traj = TrajectoryReader::open(&args.trajectory)?;

    // Process trajectory
    for (frame_idx, frame) in traj.enumerate() {
        let frame = frame?;

        // Do something with frame
        let result = process_frame(&frame, &topology);

        println!("{} {}", frame_idx, result);
    }

    Ok(())
}

fn process_frame(frame: &Frame, topology: &Topology) -> f64 {
    // Your analysis here
    0.0
}
```

### Building and Running

```bash
# Add to Cargo.toml
[[bin]]
name = "my_tool"
path = "src/bin/my_tool.rs"

# Build
cargo build --release --bin my_tool

# Run
./target/release/my_tool --topology system.top --trajectory traj.trc
```

---

## Tutorial: Your First Tool

Let's implement `rgyr` (radius of gyration calculator) step by step.

### Step 1: Define the Tool

```rust
// gromos-rs/src/bin/rgyr.rs

use gromos_rs::prelude::*;
use clap::Parser;

/// Calculate radius of gyration for selected atoms
///
/// The radius of gyration (Rg) measures the compactness of a structure:
/// Rg = sqrt(Σ m_i * r_i² / Σ m_i)
///
/// where m_i is the mass of atom i and r_i is its distance from the center of mass.
#[derive(Parser, Debug)]
#[command(name = "rgyr")]
#[command(author = "GROMOS-RS Contributors")]
#[command(version = "1.0")]
struct Args {
    /// Input topology file
    #[arg(short, long)]
    topology: String,

    /// Input trajectory file
    #[arg(short = 'x', long)]
    trajectory: String,

    /// Atom selection (e.g., "1:CA" or "1:PRO")
    #[arg(short, long, default_value = "1:all")]
    atoms: String,

    /// Use mass-weighted Rg
    #[arg(short, long, default_value_t = true)]
    mass_weighted: bool,

    /// Output file (default: stdout)
    #[arg(short, long)]
    output: Option<String>,
}
```

### Step 2: Implement Core Logic

```rust
use nalgebra::Vector3;

fn calculate_rgyr(
    positions: &[Vector3<f64>],
    masses: &[f64],
    mass_weighted: bool,
) -> f64 {
    // 1. Calculate center of mass
    let mut com = Vector3::zeros();
    let mut total_mass = 0.0;

    for (pos, &mass) in positions.iter().zip(masses) {
        let weight = if mass_weighted { mass } else { 1.0 };
        com += weight * pos;
        total_mass += weight;
    }

    com /= total_mass;

    // 2. Calculate Rg
    let mut rg2 = 0.0;

    for (pos, &mass) in positions.iter().zip(masses) {
        let r = pos - com;
        let weight = if mass_weighted { mass } else { 1.0 };
        rg2 += weight * r.norm_squared();
    }

    (rg2 / total_mass).sqrt()
}
```

### Step 3: Process Trajectory

```rust
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Load topology
    let topology = Topology::from_file(&args.topology)?;

    // Parse atom selection
    let selection = AtomSelection::from_string(&args.atoms, &topology)?;
    let atoms = selection.indices();

    // Get masses
    let masses: Vec<f64> = atoms.iter()
        .map(|&i| topology.atom(i).mass())
        .collect();

    // Open trajectory
    let traj = TrajectoryReader::open(&args.trajectory)?;

    // Setup output
    let mut output: Box<dyn std::io::Write> = match args.output {
        Some(path) => Box::new(std::fs::File::create(path)?),
        None => Box::new(std::io::stdout()),
    };

    // Write header
    writeln!(output, "# Time(ps)    Rg(nm)")?;

    // Process frames
    for frame in traj {
        let frame = frame?;

        // Extract selected atom positions
        let positions: Vec<Vector3<f64>> = atoms.iter()
            .map(|&i| frame.position(i))
            .collect();

        // Calculate Rg
        let rg = calculate_rgyr(&positions, &masses, args.mass_weighted);

        // Write result
        writeln!(output, "{:.3}      {:.6}", frame.time(), rg)?;
    }

    Ok(())
}
```

### Step 4: Add Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rgyr_sphere() {
        // Perfect sphere: Rg = sqrt(3/5) * R
        let radius = 1.0;
        let n_points = 1000;

        let mut positions = Vec::new();
        for i in 0..n_points {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n_points as f64);
            let phi = std::f64::consts::PI * (i as f64) / (n_points as f64);

            let x = radius * theta.cos() * phi.sin();
            let y = radius * theta.sin() * phi.sin();
            let z = radius * phi.cos();

            positions.push(Vector3::new(x, y, z));
        }

        let masses = vec![1.0; n_points];
        let rg = calculate_rgyr(&positions, &masses, true);

        let expected = (3.0 / 5.0).sqrt() * radius;
        assert!((rg - expected).abs() < 0.01);
    }

    #[test]
    fn test_rgyr_linear() {
        // Linear molecule along x-axis
        let positions = vec![
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        ];
        let masses = vec![1.0, 1.0, 1.0];

        let rg = calculate_rgyr(&positions, &masses, true);

        // For 3 equal masses at -1, 0, 1: Rg = sqrt(2/3)
        let expected = (2.0 / 3.0).sqrt();
        assert!((rg - expected).abs() < 1e-10);
    }
}
```

### Step 5: Add Documentation

```rust
/// Calculate radius of gyration for selected atoms
///
/// # Theory
///
/// The radius of gyration (Rg) measures the spatial distribution of atoms
/// around their center of mass:
///
/// Rg² = Σ_i m_i |r_i - r_com|² / Σ_i m_i
///
/// where:
/// - m_i is the mass of atom i
/// - r_i is the position of atom i
/// - r_com is the center of mass
///
/// # Examples
///
/// Calculate Rg for protein backbone:
/// ```bash
/// rgyr --topology protein.top --trajectory traj.trc --atoms "1:N,CA,C,O"
/// ```
///
/// Calculate unweighted Rg:
/// ```bash
/// rgyr --topology protein.top --trajectory traj.trc --no-mass-weighted
/// ```
///
/// # Output
///
/// Time series of Rg:
/// ```
/// # Time(ps)    Rg(nm)
/// 0.000      1.234
/// 2.000      1.256
/// 4.000      1.243
/// ...
/// ```
```

---

## Common Patterns

### Pattern 1: Trajectory Analysis

```rust
// Template for trajectory-based analysis

use gromos_rs::prelude::*;

pub struct MyAnalyzer {
    topology: Topology,
    selection: AtomSelection,
    // Tool-specific state
}

impl MyAnalyzer {
    pub fn new(topology: Topology, selection: AtomSelection) -> Self {
        Self { topology, selection }
    }

    pub fn analyze_frame(&mut self, frame: &Frame) -> Result<f64> {
        // Per-frame analysis
        let atoms = self.selection.indices();

        // Do calculation
        let result = 0.0; // Your calculation here

        Ok(result)
    }

    pub fn finalize(&self) -> Result<()> {
        // Post-processing, statistics, etc.
        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let topology = Topology::from_file(&args.topology)?;
    let selection = AtomSelection::from_string(&args.atoms, &topology)?;

    let mut analyzer = MyAnalyzer::new(topology, selection);

    let traj = TrajectoryReader::open(&args.trajectory)?;

    for frame in traj {
        let frame = frame?;
        let result = analyzer.analyze_frame(&frame)?;
        println!("{}", result);
    }

    analyzer.finalize()?;

    Ok(())
}
```

### Pattern 2: Pairwise Analysis (e.g., RDF, H-bonds)

```rust
// Template for analyzing pairs of atoms

pub struct PairwiseAnalyzer {
    topology: Topology,
    atoms_i: AtomSelection,
    atoms_j: AtomSelection,
    cutoff: f64,
}

impl PairwiseAnalyzer {
    pub fn analyze_frame(&self, frame: &Frame) -> Vec<(usize, usize, f64)> {
        let mut results = Vec::new();

        // Iterate over all pairs
        for &i in self.atoms_i.indices() {
            let pos_i = frame.position(i);

            for &j in self.atoms_j.indices() {
                if i == j { continue; }

                let pos_j = frame.position(j);

                // Apply PBC
                let r_ij = apply_pbc(pos_j - pos_i, frame.box_vectors());
                let distance = r_ij.norm();

                if distance < self.cutoff {
                    results.push((i, j, distance));
                }
            }
        }

        results
    }
}
```

### Pattern 3: Statistical Accumulation

```rust
// Template for accumulating statistics

pub struct StatisticsCollector {
    sum: f64,
    sum_squared: f64,
    count: usize,
    min: f64,
    max: f64,
}

impl StatisticsCollector {
    pub fn new() -> Self {
        Self {
            sum: 0.0,
            sum_squared: 0.0,
            count: 0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        }
    }

    pub fn add(&mut self, value: f64) {
        self.sum += value;
        self.sum_squared += value * value;
        self.count += 1;
        self.min = self.min.min(value);
        self.max = self.max.max(value);
    }

    pub fn mean(&self) -> f64 {
        self.sum / (self.count as f64)
    }

    pub fn std_dev(&self) -> f64 {
        let mean = self.mean();
        let variance = self.sum_squared / (self.count as f64) - mean * mean;
        variance.sqrt()
    }

    pub fn print_summary(&self) {
        println!("Statistics:");
        println!("  Count:  {}", self.count);
        println!("  Mean:   {:.6}", self.mean());
        println!("  StdDev: {:.6}", self.std_dev());
        println!("  Min:    {:.6}", self.min);
        println!("  Max:    {:.6}", self.max);
    }
}
```

---

## Advanced Examples

### Example 1: RDF (Radial Distribution Function)

Full implementation of RDF calculator:

```rust
// gromos-rs/src/bin/rdf.rs

use gromos_rs::prelude::*;
use clap::Parser;

#[derive(Parser)]
#[command(name = "rdf")]
#[command(about = "Calculate radial distribution function g(r)")]
struct Args {
    #[arg(short, long)]
    topology: String,

    #[arg(short = 'x', long)]
    trajectory: String,

    /// First atom group
    #[arg(long)]
    atoms_i: String,

    /// Second atom group
    #[arg(long)]
    atoms_j: String,

    /// Maximum distance (nm)
    #[arg(long, default_value_t = 1.5)]
    cutoff: f64,

    /// Bin width (nm)
    #[arg(long, default_value_t = 0.002)]
    bin_width: f64,
}

struct RDFCalculator {
    r_min: f64,
    r_max: f64,
    bin_width: f64,
    n_bins: usize,
    histogram: Vec<usize>,
    n_frames: usize,
    n_atoms_i: usize,
    n_atoms_j: usize,
}

impl RDFCalculator {
    fn new(r_max: f64, bin_width: f64) -> Self {
        let n_bins = (r_max / bin_width).ceil() as usize;

        Self {
            r_min: 0.0,
            r_max,
            bin_width,
            n_bins,
            histogram: vec![0; n_bins],
            n_frames: 0,
            n_atoms_i: 0,
            n_atoms_j: 0,
        }
    }

    fn analyze_frame(
        &mut self,
        frame: &Frame,
        atoms_i: &[usize],
        atoms_j: &[usize],
    ) {
        self.n_frames += 1;
        self.n_atoms_i = atoms_i.len();
        self.n_atoms_j = atoms_j.len();

        for &i in atoms_i {
            let pos_i = frame.position(i);

            for &j in atoms_j {
                if i == j { continue; }

                let pos_j = frame.position(j);
                let r_ij = apply_pbc(pos_j - pos_i, frame.box_vectors());
                let r = r_ij.norm();

                if r < self.r_max {
                    let bin = (r / self.bin_width) as usize;
                    if bin < self.n_bins {
                        self.histogram[bin] += 1;
                    }
                }
            }
        }
    }

    fn calculate_rdf(&self, box_volume: f64) -> Vec<(f64, f64)> {
        let mut rdf = Vec::with_capacity(self.n_bins);

        let rho = (self.n_atoms_j as f64) / box_volume;

        for (bin, &count) in self.histogram.iter().enumerate() {
            let r = (bin as f64 + 0.5) * self.bin_width;

            // Shell volume: 4πr²Δr
            let shell_volume = 4.0 * std::f64::consts::PI * r * r * self.bin_width;

            // Ideal count: ρ * V_shell * N_i * N_frames
            let ideal_count = rho * shell_volume *
                             (self.n_atoms_i as f64) *
                             (self.n_frames as f64);

            let g_r = if ideal_count > 0.0 {
                (count as f64) / ideal_count
            } else {
                0.0
            };

            rdf.push((r, g_r));
        }

        rdf
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Load topology
    let topology = Topology::from_file(&args.topology)?;

    // Parse selections
    let atoms_i = AtomSelection::from_string(&args.atoms_i, &topology)?;
    let atoms_j = AtomSelection::from_string(&args.atoms_j, &topology)?;

    // Initialize calculator
    let mut calculator = RDFCalculator::new(args.cutoff, args.bin_width);

    // Open trajectory
    let traj = TrajectoryReader::open(&args.trajectory)?;

    // Process frames
    let mut total_volume = 0.0;
    let mut n_frames = 0;

    for frame in traj {
        let frame = frame?;

        calculator.analyze_frame(
            &frame,
            atoms_i.indices(),
            atoms_j.indices(),
        );

        total_volume += frame.box_volume();
        n_frames += 1;
    }

    let avg_volume = total_volume / (n_frames as f64);

    // Calculate and print RDF
    let rdf = calculator.calculate_rdf(avg_volume);

    println!("# r(nm)    g(r)");
    for (r, g) in rdf {
        println!("{:.6}  {:.6}", r, g);
    }

    Ok(())
}
```

### Example 2: RMSD with Alignment

```rust
// gromos-rs/src/bin/rmsd.rs

use gromos_rs::prelude::*;
use nalgebra::{Vector3, Matrix3};

struct RMSDCalculator {
    reference: Vec<Vector3<f64>>,
    masses: Vec<f64>,
    do_fit: bool,
}

impl RMSDCalculator {
    fn new(reference: Vec<Vector3<f64>>, masses: Vec<f64>, do_fit: bool) -> Self {
        Self { reference, masses, do_fit }
    }

    fn calculate(&self, current: &[Vector3<f64>]) -> f64 {
        let mut current = current.to_vec();

        if self.do_fit {
            // Align current to reference
            self.align(&mut current);
        }

        // Calculate RMSD
        let mut sum = 0.0;
        let mut total_mass = 0.0;

        for ((&r_ref, &r_cur), &mass) in self.reference.iter()
            .zip(current.iter())
            .zip(self.masses.iter())
        {
            let diff = r_cur - r_ref;
            sum += mass * diff.norm_squared();
            total_mass += mass;
        }

        (sum / total_mass).sqrt()
    }

    fn align(&self, current: &mut [Vector3<f64>]) {
        // 1. Center both structures at origin
        let com_ref = self.center_of_mass(&self.reference);
        let com_cur = self.center_of_mass(current);

        for pos in current.iter_mut() {
            *pos -= com_cur;
        }

        // 2. Calculate optimal rotation matrix (Kabsch algorithm)
        let rotation = self.kabsch_rotation(current);

        // 3. Apply rotation
        for pos in current.iter_mut() {
            *pos = rotation * *pos;
        }
    }

    fn center_of_mass(&self, positions: &[Vector3<f64>]) -> Vector3<f64> {
        let mut com = Vector3::zeros();
        let mut total_mass = 0.0;

        for (&pos, &mass) in positions.iter().zip(self.masses.iter()) {
            com += mass * pos;
            total_mass += mass;
        }

        com / total_mass
    }

    fn kabsch_rotation(&self, current: &[Vector3<f64>]) -> Matrix3<f64> {
        // Build covariance matrix
        let mut covariance = Matrix3::zeros();

        for ((&r_ref, &r_cur), &mass) in self.reference.iter()
            .zip(current.iter())
            .zip(self.masses.iter())
        {
            let r_ref_centered = r_ref - self.center_of_mass(&self.reference);
            let r_cur_centered = r_cur - self.center_of_mass(current);

            covariance += mass * r_cur_centered * r_ref_centered.transpose();
        }

        // SVD
        let svd = covariance.svd(true, true);
        let u = svd.u.unwrap();
        let v_t = svd.v_t.unwrap();

        // Ensure right-handed coordinate system
        let mut rotation = v_t.transpose() * u.transpose();
        let det = rotation.determinant();

        if det < 0.0 {
            let mut v = v_t.transpose();
            v.column_mut(2).neg_mut();
            rotation = v * u.transpose();
        }

        rotation
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rmsd_identical() {
        let positions = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let masses = vec![1.0, 1.0, 1.0];

        let calc = RMSDCalculator::new(positions.clone(), masses, false);
        let rmsd = calc.calculate(&positions);

        assert!(rmsd < 1e-10);
    }

    #[test]
    fn test_rmsd_with_rotation() {
        // Original triangle
        let reference = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];

        // Rotated 90° around z-axis
        let rotated = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
        ];

        let masses = vec![1.0, 1.0, 1.0];

        // Without alignment: non-zero RMSD
        let calc_no_fit = RMSDCalculator::new(reference.clone(), masses.clone(), false);
        let rmsd_no_fit = calc_no_fit.calculate(&rotated);
        assert!(rmsd_no_fit > 0.1);

        // With alignment: zero RMSD (same structure)
        let calc_fit = RMSDCalculator::new(reference, masses, true);
        let rmsd_fit = calc_fit.calculate(&rotated);
        assert!(rmsd_fit < 1e-6);
    }
}
```

### Example 3: Hydrogen Bond Analysis

```rust
// gromos-rs/src/bin/hbond.rs

use gromos_rs::prelude::*;

struct HBond {
    donor_idx: usize,
    hydrogen_idx: usize,
    acceptor_idx: usize,
    distance: f64,
    angle: f64,
}

struct HBondAnalyzer {
    topology: Topology,
    donors: Vec<(usize, usize)>,      // (donor, hydrogen) pairs
    acceptors: Vec<usize>,
    max_distance: f64,
    min_angle: f64,
    hbonds: Vec<Vec<HBond>>,  // Per-frame H-bonds
}

impl HBondAnalyzer {
    fn new(
        topology: Topology,
        donors: Vec<(usize, usize)>,
        acceptors: Vec<usize>,
        max_distance: f64,
        min_angle: f64,
    ) -> Self {
        Self {
            topology,
            donors,
            acceptors,
            max_distance,
            min_angle,
            hbonds: Vec::new(),
        }
    }

    fn analyze_frame(&mut self, frame: &Frame) {
        let mut frame_hbonds = Vec::new();

        for &(donor_idx, hydrogen_idx) in &self.donors {
            let pos_d = frame.position(donor_idx);
            let pos_h = frame.position(hydrogen_idx);

            for &acceptor_idx in &self.acceptors {
                // Skip if same atom
                if donor_idx == acceptor_idx { continue; }

                let pos_a = frame.position(acceptor_idx);

                // Check distance criterion
                let r_ha = apply_pbc(pos_a - pos_h, frame.box_vectors());
                let distance = r_ha.norm();

                if distance > self.max_distance {
                    continue;
                }

                // Check angle criterion: D-H...A
                let r_dh = apply_pbc(pos_h - pos_d, frame.box_vectors());
                let angle = self.calculate_angle(&r_dh, &r_ha);

                if angle < self.min_angle {
                    continue;
                }

                // Found H-bond!
                frame_hbonds.push(HBond {
                    donor_idx,
                    hydrogen_idx,
                    acceptor_idx,
                    distance,
                    angle,
                });
            }
        }

        self.hbonds.push(frame_hbonds);
    }

    fn calculate_angle(&self, v1: &Vector3<f64>, v2: &Vector3<f64>) -> f64 {
        let cos_angle = v1.dot(v2) / (v1.norm() * v2.norm());
        cos_angle.acos().to_degrees()
    }

    fn print_statistics(&self) {
        println!("H-Bond Statistics:");
        println!("  Total frames: {}", self.hbonds.len());

        // Calculate per-pair occupancy
        let mut pair_counts: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();

        for frame_hbonds in &self.hbonds {
            for hbond in frame_hbonds {
                *pair_counts.entry((hbond.donor_idx, hbond.acceptor_idx))
                    .or_insert(0) += 1;
            }
        }

        println!("\nH-Bond Pairs:");
        println!("  Donor    Acceptor    Occupancy");

        for ((donor, acceptor), count) in pair_counts {
            let occupancy = (count as f64) / (self.hbonds.len() as f64);
            println!("  {:6}   {:6}      {:.3}", donor, acceptor, occupancy);
        }
    }
}
```

---

## Best Practices

### 1. Performance Optimization

```rust
// Use iterators (zero-cost abstractions)
let sum: f64 = values.iter().sum();  // ✓

// Not manual loops
let mut sum = 0.0;  // ✗
for &v in &values {
    sum += v;
}

// Use parallel iteration for large datasets
use rayon::prelude::*;

let results: Vec<_> = frames.par_iter()
    .map(|frame| analyze_frame(frame))
    .collect();

// Pre-allocate when size is known
let mut results = Vec::with_capacity(n_frames);  // ✓
let mut results = Vec::new();  // ✗ (will reallocate)

// Use references to avoid cloning
fn process(data: &[f64]) { /* ... */ }  // ✓
fn process(data: Vec<f64>) { /* ... */ }  // ✗ (takes ownership)
```

### 2. Error Handling

```rust
// Use Result for fallible operations
fn read_trajectory(path: &str) -> Result<Trajectory, Error> {
    // ...
}

// Provide context with anyhow
use anyhow::{Context, Result};

let traj = TrajectoryReader::open(&args.trajectory)
    .context("Failed to open trajectory file")?;

// Handle specific errors
match result {
    Ok(value) => { /* ... */ }
    Err(Error::FileNotFound) => { /* ... */ }
    Err(e) => return Err(e),
}
```

### 3. Documentation

```rust
/// Calculate property X from trajectory.
///
/// # Arguments
///
/// * `topology` - System topology
/// * `trajectory` - Input trajectory path
/// * `selection` - Atom selection string
///
/// # Returns
///
/// Time series of property X
///
/// # Example
///
/// ```bash
/// my_tool --topology system.top --trajectory traj.trc --atoms "1:CA"
/// ```
///
/// # References
///
/// [1] Author et al., Journal, Year
pub fn calculate_x(/* ... */) -> Result<Vec<f64>> {
    // ...
}
```

### 4. Testing

```rust
#[cfg(test)]
mod tests {
    use super::*;

    // Test with simple cases
    #[test]
    fn test_simple_case() {
        let result = my_function(&simple_input);
        assert_eq!(result, expected);
    }

    // Test edge cases
    #[test]
    fn test_empty_input() {
        let result = my_function(&[]);
        assert!(result.is_empty());
    }

    // Test with real data (integration test)
    #[test]
    fn test_with_real_trajectory() {
        let traj = TrajectoryReader::open("tests/data/test.trc").unwrap();
        // ...
    }

    // Benchmark critical paths
    #[bench]
    fn bench_my_function(b: &mut Bencher) {
        b.iter(|| my_function(&large_input));
    }
}
```

---

## Testing Your Tool

### Unit Tests

```bash
# Run tests
cargo test --bin my_tool

# Run with output
cargo test --bin my_tool -- --nocapture

# Run specific test
cargo test --bin my_tool test_name
```

### Integration Tests

```bash
# Test with real data
./target/release/my_tool \
    --topology tests/data/system.top \
    --trajectory tests/data/traj.trc \
    > output.dat

# Compare with reference
diff output.dat tests/data/reference.dat
```

### Performance Benchmarks

```bash
# Time execution
time ./target/release/my_tool --topology large.top --trajectory large.trc

# Profile with flamegraph
cargo install flamegraph
cargo flamegraph --bin my_tool -- --topology system.top --trajectory traj.trc
```

---

## Publishing Your Tool

### 1. Documentation

- Add comprehensive doc comments
- Include examples in documentation
- Write user guide (if complex)

### 2. Tests

- Unit tests for all functions
- Integration tests with real data
- Benchmarks for performance-critical code

### 3. Submit Pull Request

See [Contributing Guide](contributing.md) for process.

---

## Summary

You now know how to write gromos++ style analysis tools in Rust!

**Key Takeaways**:
- ✅ Follow Unix philosophy (do one thing well)
- ✅ Use standard data structures (Topology, Frame, etc.)
- ✅ Optimize with iterators and parallelization
- ✅ Handle errors properly
- ✅ Write tests and documentation
- ✅ Follow Rust best practices

**Next Steps**:
1. Read [Rust for GROMOS-RS](rust-for-gromos-rs.md) for language basics
2. Study existing tools in `gromos-rs/src/bin/`
3. Implement your own tool!
4. Share it with the community

---

**Happy coding!** 🦀
