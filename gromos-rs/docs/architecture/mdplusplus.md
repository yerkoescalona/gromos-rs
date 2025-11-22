# MD++ Architecture

MD++ is the core simulation engine of the GROMOS software suite, written in C++ with over 300,000 lines of code developed over 20+ years.

## Overview

**MD++** (Molecular Dynamics Plus Plus) is a feature-rich, high-performance molecular dynamics simulation package designed for biomolecular systems.

### Key Statistics

| Metric | Value |
|--------|-------|
| **Language** | C++11/14 |
| **Code Size** | ~300,000 lines |
| **Development** | 20+ years (2000-present) |
| **License** | GPL-2.0 |
| **Main Binary** | `md++` |
| **Primary Use** | Production MD simulations |
| **Parallelization** | OpenMP (threads), MPI (multi-node), CUDA (GPU) |

## Directory Structure

```
md++/
├── src/                          # Source code (~300K lines)
│   ├── algorithm/                # Algorithms (80K lines)
│   │   ├── integration/          # Integrators (13 algorithms)
│   │   │   ├── leapfrog.cc      # Leap-Frog integrator
│   │   │   ├── velocityverlet.cc # Velocity Verlet
│   │   │   ├── stochastic_dynamics.cc # 4 SD variants
│   │   │   ├── minimize/        # Minimization
│   │   │   │   ├── steepest_descent.cc
│   │   │   │   └── conjugate_gradient.cc
│   │   │   ├── scaled_leapfrog.cc # Multiple time-stepping
│   │   │   ├── monte_carlo.cc    # MC sampling
│   │   │   ├── multigradient.cc  # Multi-potential
│   │   │   ├── lattice_shift.cc  # FEP long-range
│   │   │   ├── eds.cc            # Enveloping distribution
│   │   │   └── gamd.cc           # Gaussian accelerated
│   │   │
│   │   ├── constraints/          # Constraint algorithms (50K lines)
│   │   │   ├── shake.cc          # Standard SHAKE
│   │   │   ├── m_shake.cc        # Mass-weighted SHAKE
│   │   │   ├── settle.cc         # Rigid water (analytical)
│   │   │   ├── lincs.cc          # Linear constraint solver
│   │   │   ├── perturbed_shake.cc # λ-dependent constraints
│   │   │   ├── flexible_constraint.cc # Time-dependent
│   │   │   ├── remove_com_motion.cc # COM removal
│   │   │   ├── angle_constraint.cc
│   │   │   └── dihedral_constraint.cc
│   │   │
│   │   ├── temperature/          # Thermostats
│   │   │   ├── berendsen.cc     # Weak coupling
│   │   │   ├── nose_hoover.cc   # Extended system
│   │   │   └── andersen.cc      # Stochastic collisions
│   │   │
│   │   ├── pressure/            # Barostats
│   │   │   ├── berendsen.cc     # Weak coupling
│   │   │   └── parrinello_rahman.cc # Extended system
│   │   │
│   │   └── virtualatoms/        # Virtual site handling
│   │       ├── virtualatom_type.cc
│   │       └── force_redistribution.cc
│   │
│   ├── interaction/             # Force calculations (120K lines)
│   │   ├── bonded/              # Bonded interactions
│   │   │   ├── quartic_bond.cc  # Quartic bond potential
│   │   │   ├── harmonic_bond.cc # Harmonic bond potential
│   │   │   ├── angle.cc         # Angle potentials (cosine/harmonic)
│   │   │   ├── dihedral.cc      # Dihedral potentials
│   │   │   ├── improper.cc      # Improper dihedrals
│   │   │   ├── crossdihedral.cc # Cross-dihedral coupling
│   │   │   └── perturbed_*.cc   # FEP variants of all above
│   │   │
│   │   ├── nonbonded/           # Nonbonded interactions
│   │   │   ├── pairlist.cc      # Neighbor list generation
│   │   │   ├── grid_cell_pairlist.cc # O(N) spatial decomp
│   │   │   ├── innerloop.cc     # LJ + Coulomb kernel
│   │   │   ├── perturbed_innerloop.cc # FEP soft-core
│   │   │   ├── rf_innerloop.cc  # Reaction field
│   │   │   └── polarisation.cc  # Polarizable force fields
│   │   │
│   │   ├── latticesum/          # Long-range electrostatics
│   │   │   ├── ewald.cc         # Direct Ewald summation
│   │   │   ├── pme.cc           # Particle Mesh Ewald
│   │   │   └── p3m.cc           # Particle-Particle Particle-Mesh
│   │   │
│   │   ├── qmmm/                # QM/MM hybrid simulations
│   │   │   ├── qm_zone.cc       # QM region definition
│   │   │   ├── qm_link.cc       # Link atoms
│   │   │   ├── qm_atom.cc       # QM atom management
│   │   │   ├── dftb_worker.cc   # DFTB+ interface (60 KB)
│   │   │   ├── gaussian_worker.cc # Gaussian interface
│   │   │   ├── mopac_worker.cc  # MOPAC interface
│   │   │   ├── mndo_worker.cc   # MNDO interface
│   │   │   ├── orca_worker.cc   # ORCA interface
│   │   │   ├── turbomole_worker.cc # Turbomole interface
│   │   │   ├── xtb_worker.cc    # xTB interface
│   │   │   ├── nn_worker.cc     # Neural network potential
│   │   │   └── nonbonded_innerloop_qmmm.cc # QM-MM interactions
│   │   │
│   │   └── special/             # Special interactions (200K lines!)
│   │       ├── distance_restraint.cc # Distance restraints
│   │       ├── position_restraint.cc # Position restraints
│   │       ├── angle_restraint.cc
│   │       ├── dihedral_restraint.cc
│   │       ├── jvalue_restraint.cc (50 KB) # NMR J-coupling
│   │       ├── rdc_restraint.cc (138 KB!) # Residual dipolar coupling
│   │       ├── order_parameter_restraint.cc # Lipid ordering
│   │       ├── xray_restraint.cc (90 KB) # X-ray refinement
│   │       ├── local_elevation.cc (75 KB) # Metadynamics
│   │       ├── distance_field_interaction.cc # Biasing
│   │       ├── electric_field_interaction.cc # External E-field
│   │       ├── nemd.cc          # Non-equilibrium MD
│   │       └── symmetry_restraint.cc # Symmetry enforcement
│   │
│   ├── topology/                # System topology (40K lines)
│   │   ├── topology.cc          # Main topology class
│   │   ├── atom.cc              # Atom definitions
│   │   ├── bond.cc              # Bond list
│   │   ├── molecule.cc          # Molecule organization
│   │   ├── virtualatom_type.cc  # Virtual sites
│   │   ├── exclusion.cc         # Exclusion lists
│   │   └── perturbed_topology.cc # FEP topology
│   │
│   ├── simulation/              # Simulation control (30K lines)
│   │   ├── simulation.cc        # Main simulation loop
│   │   ├── multistep.cc         # Multi-step integrator
│   │   └── parameter.cc         # Parameter handling
│   │
│   ├── io/                      # Input/Output (60K lines)
│   │   ├── topology/
│   │   │   ├── in_topology.cc   # .top reader
│   │   │   └── in_perturbation.cc # .ptp reader
│   │   ├── configuration/
│   │   │   ├── in_configuration.cc # .cnf reader
│   │   │   └── out_configuration.cc # .cnf writer
│   │   ├── parameter/
│   │   │   └── in_parameter.cc  # .imd reader (complex!)
│   │   ├── trajectory/
│   │   │   └── out_trajectory.cc # .trc writer
│   │   ├── energy/
│   │   │   └── out_energy.cc    # .tre writer
│   │   ├── force/
│   │   │   └── out_force.cc     # .trf writer
│   │   ├── gzstream/            # Gzip compression
│   │   └── pdb/                 # PDB format
│   │
│   ├── cukernel/                # CUDA GPU kernels (80K lines)
│   │   ├── nonbonded_kernel.cu  # GPU nonbonded (200 KB)
│   │   ├── shake_kernel.cu      # GPU SHAKE (80 KB)
│   │   ├── settle_kernel.cu     # GPU SETTLE (40 KB)
│   │   ├── rf_kernel.cu         # GPU reaction field
│   │   └── pme_kernel.cu        # GPU PME
│   │
│   ├── replicaExchange/         # Replica exchange (25K lines)
│   │   ├── replica_exchange.cc  # Main REMD
│   │   ├── replica/
│   │   │   ├── temperature_replica.cc # T-REMD
│   │   │   ├── lambda_replica.cc # λ-REMD
│   │   │   ├── eds_replica.cc   # EDS-REMD
│   │   │   └── adde_replica.cc  # ADDE reweighting
│   │   └── mpi_controller.cc    # MPI communication
│   │
│   ├── math/                    # Math utilities (15K lines)
│   │   ├── volume.cc            # Volume calculations
│   │   ├── periodicity.cc       # PBC handling
│   │   ├── transformation.cc    # Coordinate transforms
│   │   └── random.cc            # RNG
│   │
│   └── util/                    # Utilities (20K lines)
│       ├── timing.cc            # Performance timing
│       ├── error.cc             # Error handling
│       ├── debug.cc             # Debug output
│       └── template_split.cc    # Template parsing
│
├── programs/                    # Main executables
│   └── md.cc                    # md++ main program (2K lines)
│
├── config/                      # Build configuration
│   ├── Config.sh               # Configuration script
│   └── Makefile.am             # Automake files
│
└── CMakeLists.txt              # CMake build system

```

## Core Components

### 1. Integration Algorithms (13 total)

MD++ implements a comprehensive set of integrators:

#### Standard Molecular Dynamics

**Leap-Frog** (`algorithm/integration/leapfrog.cc`):
```cpp
// Leap-frog integration scheme
// v(t+dt/2) = v(t-dt/2) + a(t)*dt
// r(t+dt) = r(t) + v(t+dt/2)*dt

class Leapfrog_Algorithm {
    void apply(topology::Topology &topo,
               configuration::Configuration &conf,
               simulation::Simulation &sim);
};
```

- Most common MD integrator
- Time-reversible, symplectic
- Energy conservation in NVE
- Used in >90% of GROMOS simulations

**Velocity Verlet** (`algorithm/integration/velocityverlet.cc`):
- More accurate than leap-frog for the same timestep
- Better for Nosé-Hoover thermostat
- Synchronous positions and velocities

**Stochastic Dynamics** (`algorithm/integration/stochastic_dynamics.cc`):
- 4 different variants implemented
- Langevin dynamics with friction
- Implicit solvent simulations
- Temperature control without separate thermostat

#### Energy Minimization

**Steepest Descent** (`algorithm/integration/minimize/steepest_descent.cc`):
```cpp
// Steepest descent minimization
// r_new = r_old - α * ∇E

class Steepest_Descent {
    double step_size;
    double tolerance;
    int max_iterations;
};
```

**Conjugate Gradient** (`algorithm/integration/minimize/conjugate_gradient.cc`):
- Much faster convergence than steepest descent
- Fletcher-Reeves or Polak-Ribiere variants
- Line search with golden section or polynomial fit

#### Advanced Methods

**Scaled Leap-Frog** (`algorithm/integration/scaled_leapfrog.cc`):
- Multiple time-stepping (MTS)
- Fast forces: every step
- Slow forces: every N steps
- 2-3x speedup for large systems

**Monte Carlo** (`algorithm/integration/monte_carlo.cc`):
- Metropolis sampling
- Hybrid MC/MD
- Constant pressure MC moves

**Multi-Gradient** (`algorithm/integration/multigradient.cc`):
- Multiple potential energy surfaces
- Interpolation between states
- Advanced free energy methods

**Lattice Shift** (`algorithm/integration/lattice_shift.cc`):
- Track PBC crossings for FEP
- Required for long-range FEP calculations
- Charge transfer across periodic boundaries

### 2. Constraint Algorithms (9 total)

**SHAKE** (`algorithm/constraints/shake.cc`):
```cpp
// SHAKE algorithm for bond constraints
// Iteratively satisfy: |r_ij| = d_ij

class Shake {
    double tolerance = 1e-4;  // Convergence tolerance
    int max_iterations = 1000; // Max iterations

    void apply(std::vector<Atom> &atoms,
               std::vector<Bond> &constraints);
};
```

- Iterative constraint satisfaction
- Typically converges in 3-5 iterations
- Can constrain bonds, angles, dihedrals

**M-SHAKE** (Mass-weighted variant):
- Reduced mass formulation
- Better for heterogeneous systems
- Faster convergence for water

**SETTLE** (`algorithm/constraints/settle.cc`):
- Analytical solution for rigid water
- No iteration required
- Exact to machine precision
- 5-10x faster than SHAKE for water

**LINCS** (`algorithm/constraints/lincs.cc`):
- Linear constraint solver
- Matrix formulation
- Better parallelization than SHAKE
- Popular in GROMACS

**Perturbed SHAKE** (`algorithm/constraints/perturbed_shake.cc`):
```cpp
// λ-dependent constraint lengths for FEP
// d_ij(λ) = (1-λ)*d_ij^A + λ*d_ij^B

class Perturbed_Shake : public Shake {
    double lambda;  // FEP parameter
    std::vector<double> d_A, d_B;  // State A and B lengths
};
```

**COM Motion Removal** (`algorithm/constraints/remove_com_motion.cc`):
- Remove center-of-mass translation/rotation
- Prevent system drift
- Important for long simulations

### 3. Nonbonded Interactions

#### Pairlist Generation

**Standard Pairlist** (`interaction/nonbonded/pairlist.cc`):
- Chargegroup-based cutoff
- Update frequency: every 5-10 steps
- O(N²) scaling for small systems

**Grid Cell Pairlist** (`interaction/nonbonded/grid_cell_pairlist.cc`):
```cpp
// O(N) spatial decomposition using grid cells

class Grid_Cell_Pairlist {
    struct Cell {
        std::vector<int> atoms;  // Atoms in this cell
        std::vector<int> neighbors;  // Neighboring cells
    };

    std::vector<Cell> grid;
    double cell_size;  // ≥ cutoff distance

    void build_pairlist(const std::vector<Vector> &pos);
};
```

- O(N) scaling
- Essential for large systems (>10K atoms)
- Cell size ≥ cutoff distance

#### Force Kernels

**LJ + Coulomb Inner Loop** (`interaction/nonbonded/innerloop.cc`):
```cpp
// Lennard-Jones + Coulomb interaction
// E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] = C₁₂/r¹² - C₆/r⁶
// E_Coul = q_i*q_j / (4πε₀r)

inline void lj_crf_innerloop(
    const Vector &r_ij,
    const double q_i, const double q_j,
    const double c6, const double c12,
    double &energy, Vector &force)
{
    const double r2 = r_ij.abs2();
    const double r2_inv = 1.0 / r2;
    const double r6_inv = r2_inv * r2_inv * r2_inv;
    const double r12_inv = r6_inv * r6_inv;

    // LJ interaction
    const double e_lj = c12 * r12_inv - c6 * r6_inv;
    const double f_lj = (12*c12*r12_inv - 6*c6*r6_inv) * r2_inv;

    // Coulomb with reaction field
    const double e_crf = crf_2cut3i * r2 - crf_cut3i;
    const double f_crf = 2 * crf_2cut3i;

    energy += e_lj + q_i * q_j * e_crf;
    force += (f_lj + q_i * q_j * f_crf) * r_ij;
}
```

**Soft-Core FEP** (`interaction/nonbonded/perturbed_innerloop.cc`):
```cpp
// Soft-core potential prevents singularities at λ=0 or λ=1
// r_sc² = r² + α_LJ * σ⁶ * λᵖ * (1-λ)ᵖ
// E(λ) = (1-λ)*E_A(r_sc) + λ*E_B(r_sc)

inline double softcore_distance(
    double r2,
    double sigma6,
    double lambda,
    double alpha = 0.5,  // Soft-core parameter
    int power = 2)       // Typically 1 or 2
{
    double lambda_term = std::pow(lambda, power) *
                        std::pow(1.0 - lambda, power);
    return std::sqrt(r2 + alpha * sigma6 * lambda_term);
}
```

### 4. Long-Range Electrostatics

**Reaction Field** (`interaction/latticesum/rf.cc`):
```cpp
// Reaction field method
// Continuum dielectric beyond cutoff
// E_RF = q_i*q_j * [1/r + (ε-1)/(2ε+1) * r²/R_c³ - 3ε/(2ε+1) * 1/R_c]

class Reaction_Field {
    double cutoff;           // R_c
    double epsilon_rf;       // ε (typically 61 for water)
    double crf_2cut3i;      // 2/(2ε+1) * 1/R_c³
    double crf_cut3i;       // 3ε/(2ε+1) * 1/R_c
};
```

- GROMOS traditional method
- Fast, simple
- Good for homogeneous systems
- Typical ε_RF = 61 (water permittivity)

**Particle Mesh Ewald (PME)** (`interaction/latticesum/pme.cc`):
```cpp
// PME: Ewald summation using FFT
// E_total = E_real + E_reciprocal + E_self

class PME {
    // Real space (short-range, direct sum)
    double alpha;           // Ewald parameter (typically 0.3-0.4 nm⁻¹)
    double real_cutoff;     // Real space cutoff

    // Reciprocal space (long-range, FFT)
    int grid_x, grid_y, grid_z;  // FFT grid dimensions
    int spline_order;       // B-spline order (typically 4)
    fftw_plan forward_plan; // FFT forward transform
    fftw_plan backward_plan;// FFT backward transform

    void calculate_energy_forces();
};
```

- Exact treatment of long-range electrostatics
- O(N log N) scaling via FFT
- Essential for charged systems
- Widely used standard

**P3M** (`interaction/latticesum/p3m.cc`):
- Alternative to PME
- Different charge assignment
- Similar performance

### 5. QM/MM Implementation

MD++ has one of the most comprehensive QM/MM implementations:

**QM/MM Architecture** (`interaction/qmmm/`):

```cpp
// QM/MM hybrid simulation architecture

class QM_Zone {
    std::vector<int> qm_atoms;      // Atoms in QM region
    std::vector<int> mm_atoms;      // Atoms in MM region
    std::vector<QM_Link> link_atoms; // Link atoms at QM/MM boundary

    std::string qm_software;  // Which QM engine to use
    double charge_total;      // Total QM region charge
};

class QM_Worker {
    virtual void calculate_qm_energy_gradients(
        const std::vector<Vector> &qm_coords,
        const std::vector<double> &mm_charges,
        double &qm_energy,
        std::vector<Vector> &qm_gradients) = 0;
};

// Specific implementations:
class DFTB_Worker : public QM_Worker { /*...*/ };
class Gaussian_Worker : public QM_Worker { /*...*/ };
class ORCA_Worker : public QM_Worker { /*...*/ };
class xTB_Worker : public QM_Worker { /*...*/ };
// ... 9 total engines
```

**Supported QM Engines**:

1. **DFTB+** (60 KB code):
   - Density Functional Tight Binding
   - Fast, semi-empirical
   - Good for organic molecules
   - Popular choice

2. **Gaussian**:
   - Ab initio quantum chemistry
   - High accuracy
   - Expensive, for small QM regions

3. **ORCA**:
   - Modern QM package
   - DFT, coupled cluster, etc.
   - Good parallelization

4. **xTB**:
   - Extended tight-binding
   - Very fast
   - Good for conformational sampling

5. **MOPAC/MNDO**:
   - Semi-empirical methods (AM1, PM3, etc.)
   - Fast, moderate accuracy

6. **Turbomole**:
   - High-performance DFT
   - Good for large QM regions

7. **Neural Network Potential**:
   - ML-based QM
   - Speed of MM, accuracy approaching QM

8. **Ghost Worker**:
   - Testing/debugging

**Link Atoms** (`interaction/qmmm/qm_link.cc`):
```cpp
// Link atoms at QM/MM boundary
// Typically hydrogen atoms capping broken bonds

struct QM_Link {
    int qm_atom;       // QM atom in bond
    int mm_atom;       // MM atom in bond
    int link_atom;     // H atom capping the bond
    double scale;      // Position scaling (typically 0.7 for C-H)

    Vector position() const {
        return qm_pos + scale * (mm_pos - qm_pos);
    }
};
```

### 6. Special Interactions (>200K lines!)

**NMR Restraints** (`interaction/special/jvalue_restraint.cc`, `rdc_restraint.cc`):

```cpp
// J-value (scalar coupling) restraints
// Karplus relation: J(φ) = A*cos²(φ) + B*cos(φ) + C

class Jvalue_Restraint {
    double phi;           // Dihedral angle
    double J_exp;         // Experimental J-value
    double J_calc;        // Calculated J-value
    double A, B, C;       // Karplus parameters
    double force_const;   // Restraint force constant

    double energy() {
        J_calc = A*cos(phi)*cos(phi) + B*cos(phi) + C;
        double delta = J_calc - J_exp;
        return 0.5 * force_const * delta * delta;
    }
};

// RDC (Residual Dipolar Coupling) restraints (138 KB!)
// Very complex: alignment tensor, molecular orientation, averaging

class RDC_Restraint {
    // Alignment tensor (5 independent components)
    double Szz, Sxx_minus_Syy, Sxy, Sxz, Syz;

    // RDC calculation involves:
    // - Molecular orientation
    // - Internuclear vector
    // - Alignment tensor
    // - Averaging over ensemble

    double calculate_rdc(Vector r_ij);
    void update_alignment_tensor();
};
```

**X-ray Restraints** (`interaction/special/xray_restraint.cc` - 90 KB):
- Structure factors
- Electron density maps
- R-factor optimization
- Refinement against diffraction data

**Local Elevation** (`interaction/special/local_elevation.cc` - 75 KB):
```cpp
// Local elevation (metadynamics variant)
// Add Gaussian hills to escape energy minima

class Local_Elevation {
    struct GaussianHill {
        std::vector<double> center;  // Center in CV space
        double height;                // Hill height
        double width;                 // Gaussian width
    };

    std::vector<GaussianHill> hills;
    std::vector<CollectiveVariable> cvs;

    double buildup_energy;  // Energy from accumulated hills

    void add_hill(const std::vector<double> &position);
};
```

### 7. GPU/CUDA Acceleration

**CUDA Kernels** (`cukernel/` - 80K lines):

```cuda
// GPU nonbonded kernel
__global__ void nonbonded_kernel(
    const float3 *positions,
    const float *charges,
    const uint *iac,
    const uint2 *pairlist,
    const double *lj_params,
    float3 *forces,
    double *energies,
    int n_pairs)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    uint2 pair = pairlist[tid];
    uint i = pair.x;
    uint j = pair.y;

    // Load positions
    float3 ri = positions[i];
    float3 rj = positions[j];

    // Calculate distance
    float3 rij = {rj.x - ri.x, rj.y - ri.y, rj.z - ri.z};
    float r2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;

    // LJ + Coulomb (vectorized)
    float r2_inv = 1.0f / r2;
    float r6_inv = r2_inv * r2_inv * r2_inv;

    uint iac_i = iac[i];
    uint iac_j = iac[j];
    double c6 = lj_params[iac_i * n_types + iac_j];
    double c12 = lj_params[iac_i * n_types + iac_j + n_types*n_types];

    float e_lj = c12 * r6_inv * r6_inv - c6 * r6_inv;
    float e_coul = charges[i] * charges[j] * sqrtf(r2_inv);

    // Atomic add to forces (coalesced writes)
    atomicAdd(&forces[i].x, f * rij.x);
    atomicAdd(&forces[i].y, f * rij.y);
    atomicAdd(&forces[i].z, f * rij.z);
    atomicAdd(&forces[j].x, -f * rij.x);
    atomicAdd(&forces[j].y, -f * rij.y);
    atomicAdd(&forces[j].z, -f * rij.z);
}
```

**GPU Performance**:
- 10-50x speedup for nonbonded
- 20-30x speedup for SHAKE
- Requires CUDA-capable GPU (NVIDIA only)
- Best for large systems (>50K atoms)

### 8. Replica Exchange (MPI)

**Parallel Tempering** (`replicaExchange/replica/temperature_replica.cc`):

```cpp
// Temperature replica exchange
// Multiple replicas at different temperatures
// Exchange based on Metropolis criterion

class Temperature_Replica_Exchange {
    int n_replicas;
    std::vector<double> temperatures;
    std::vector<Simulation> replicas;

    // MPI communication
    MPI_Comm comm;
    int rank, size;

    void attempt_exchange() {
        // Calculate exchange probability
        // P = min(1, exp[(1/kT_i - 1/kT_j)(E_j - E_i)])

        for (int i = 0; i < n_replicas - 1; i += 2) {
            double T_i = temperatures[i];
            double T_j = temperatures[i+1];
            double E_i = replicas[i].energy();
            double E_j = replicas[i+1].energy();

            double delta = (1.0/(k_B*T_i) - 1.0/(k_B*T_j)) * (E_j - E_i);
            double p_exchange = std::min(1.0, std::exp(delta));

            if (random() < p_exchange) {
                swap_configurations(i, i+1);
                accepted_exchanges[i]++;
            }
            total_attempts[i]++;
        }
    }
};
```

**MPI Parallelization**:
- Each replica on different node
- Periodic exchange attempts
- Efficient communication patterns
- Scalability to 100s-1000s of replicas

## Performance Characteristics

### Typical Performance (16-core CPU)

| System Size | md++ (ns/day) | Notes |
|-------------|---------------|-------|
| 2.6K atoms (lysozyme) | ~150 | With RF |
| 23K atoms (DHFR) | ~45 | With RF |
| 85K atoms (membrane) | ~12 | With RF |
| 23K atoms (DHFR) | ~25 | With PME (slower) |

### With GPU (CUDA)

| System Size | CPU (ns/day) | GPU (ns/day) | Speedup |
|-------------|--------------|--------------|---------|
| 23K atoms | ~45 | ~600 | 13x |
| 85K atoms | ~12 | ~350 | 29x |
| 250K atoms | ~3 | ~120 | 40x |

### Scalability

**OpenMP (Shared Memory)**:
- Good scaling to 16-32 cores
- Diminishing returns beyond 32 cores
- Best for single-node workstations

**MPI (Distributed Memory)**:
- Excellent scaling to 1000s of cores
- Required for very large systems (>500K atoms)
- Domain decomposition
- HPC clusters

## Build System

### CMake Build

```bash
cd md++
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DOMP=ON \           # OpenMP parallelization
    -DMPI=ON \           # MPI support
    -DCUDA=ON \          # GPU support
    -DQM_DFTB=ON \       # DFTB+ QM/MM
    -DQM_XTBA=ON         # xTB QM/MM

make -j$(nproc)
```

### Configuration Options

Over 100 CMake options for features:
- Parallelization (OpenMP, MPI, CUDA)
- QM/MM engines (9 options)
- Precision (single/double/mixed)
- Debugging flags
- Performance profiling

## Input File Format (.imd)

MD++ uses complex input files with hundreds of parameters:

```
# Example md.imd - Full MD++ input file

TITLE
    Protein in water - NPT MD
END

SYSTEM
    NPM             1          # Number of protein molecules
    NSM             8000       # Number of solvent molecules
END

INITIALISE
    NTIVEL          0          # Initialize velocities (0=from file)
    NTISHK          0          # Initialize SHAKE (0=from file)
    NTINHT          0          # Initialize Nose-Hoover (0=from file)
    NTINHB          0          # Initialize Nose-Hoover barostat
    NTISHI          0          # Initialize stochastic integral
    NTIRTC          0          # Initialize position restraints
    NTIRDC          0          # Initialize distance restraints
    NTICOM          0          # Initialize COM motion
    NTISTI          0          # Initialize stochastic dynamics
    IG              12345      # Random seed
    TEMPI           300.0      # Initial temperature
END

STEP
    NSTLIM          500000     # Number of steps (500K = 1 ns)
    T               0.0        # Start time
    DT              0.002      # Timestep (2 fs)
END

BOUNDCOND
    NTB             1          # Boundary type (1=rectangular)
    NDFMIN          0          # Removed DOF
END

MULTIBATH
    ALGORITHM       1          # Bath coupling (1=Berendsen)
    NBATHS          2          # Number of baths

    # Bath 1: Solute
    TEMP0           300.0
    TAU             0.1
    DOFSET          1
    LAST            2641
    COMBATH         0
    IRBATH          0

    # Bath 2: Solvent
    TEMP0           300.0
    TAU             0.1
    DOFSET          2642
    LAST            26641
    COMBATH         0
    IRBATH          0
END

PRESSURESCALE
    COUPLE          1          # Pressure coupling (1=Berendsen)
    SCALE           1          # Scaling type (1=isotropic)
    COMP            4.575e-4   # Compressibility (water)
    TAUP            0.5        # Coupling time
    VIRIAL          2          # Virial calc (2=atomic)
    SEMIANISOTROPIC 0 0 0
    PRES0           0.06102 0.06102 0.06102  # Pressure (bar → gromos units)
END

FORCE
    BONDS           0          # Bond forces (0=all constrained)
    ANGLES          1          # Angle forces
    IMPROPER        1          # Improper dihedral forces
    DIHEDRAL        1          # Dihedral forces
    ELECTROSTATIC   1          # Electrostatic forces
    VDW             1          # van der Waals forces
    NEGR            50         # Energy groups
END

CONSTRAINT
    NTC             3          # Constraint algorithm (3=all bonds)
    NTCP            1          # Parallel SHAKE (1=on)
    NTCP0           0          # SHAKE reference (0=initial)
    NTCS            1          # Solvent SETTLE (1=on)
    NTCS0           0          # SETTLE reference
END

PAIRLIST
    ALGORITHM       2          # Pairlist type (2=grid cell O(N))
    NSNB            5          # Frequency (every 5 steps)
    RCUTP           0.8        # Short-range cutoff (nm)
    RCUTL           1.4        # Long-range cutoff (nm)
    SIZE            0.4        # Grid cell size
    TYPE            0          # Pairlist type
END

NONBONDED
    NLRELE          1          # Long-range electrostatics (1=RF)
    APPAK           0.0        # Reaction field κ parameter
    RCRF            1.4        # RF cutoff (nm)
    EPSRF           61.0       # RF permittivity (water)
    NSLFEXCL        1          # Exclude 1-4 LJ
    NSHEX           1          # Exclude 1-4 Coulomb
    NKX             0          # Ewald kx (0=auto)
    NKY             0          # Ewald ky
    NKZ             0          # Ewald kz
    KAPPA           0.0        # Ewald α (0=auto)
    NGX             32         # PME grid x
    NGY             32         # PME grid y
    NGZ             32         # PME grid z
    NASORD          4          # PME B-spline order
    NQEVAL          0          # QM evaluation frequency
    FACCOUL         1.0        # Coulomb scaling
    NSCALE          0          # Scaling in SHAKE
END

PRINTOUT
    NTPR            1000       # Energy every 1000 steps = 2 ps
    NTPP            0          # Block averages
END

WRITETRAJ
    NTWX            500        # Coords every 500 steps = 1 ps
    NTWSE           0          # Skip solvent in trajectory
    NTWV            0          # Write velocities
    NTWF            0          # Write forces
    NTWE            500        # Write energy every 500 steps
    NTWG            0          # Write energy groups
    NTWB            0          # Write blocks
END
```

Over 50 input blocks, 200+ parameters!

## Summary

### Strengths

1. **Feature completeness**: 20+ years of development
2. **QM/MM**: 9 QM engines, comprehensive implementation
3. **GPU acceleration**: Mature CUDA kernels
4. **MPI parallelization**: Scales to 1000s of cores
5. **Special interactions**: NMR, X-ray, metadynamics, etc.
6. **Well-tested**: Used in 1000s of publications
7. **Documentation**: Extensive manual (500+ pages)

### Weaknesses

1. **Code complexity**: 300K lines, hard to maintain
2. **Manual memory management**: Potential bugs
3. **Limited SIMD**: Not extensively vectorized
4. **Build complexity**: Many dependencies
5. **CUDA-only**: No OpenCL, no AMD/Intel GPUs
6. **Threading issues**: Potential race conditions (manual OpenMP)

### Use Cases

**Ideal for**:
- QM/MM simulations
- NMR structure refinement
- X-ray crystallography refinement
- Large-scale HPC (>1000 cores)
- GPU-accelerated MD (>50K atoms)
- Metadynamics, advanced sampling

**Less ideal for**:
- Quick prototyping
- New algorithms (code complexity)
- Non-NVIDIA GPUs
- Learning MD programming

---

**Next**: [GROMOS++ Architecture](gromosplusplus.md) | [GROMOS-RS Architecture](gromos-rs.md) | [Comparison](comparison.md)
