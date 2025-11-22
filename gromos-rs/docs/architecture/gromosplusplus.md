# GROMOS++ Architecture

GROMOS++ is the analysis and pre-processing toolkit of the GROMOS software suite, providing 104 command-line programs for molecular dynamics analysis.

## Overview

**GROMOS++** (GROMOS Plus Plus) is a collection of C++ programs for preparing, running, and analyzing molecular dynamics simulations.

### Key Statistics

| Metric | Value |
|--------|-------|
| **Language** | C++11 |
| **Code Size** | ~250,000 lines |
| **Programs** | 104 tools |
| **Development** | 20+ years (2000-present) |
| **License** | GPL-2.0 |
| **Primary Use** | MD analysis & preprocessing |
| **Dependencies** | Minimal (standard C++ library) |

### Philosophy

GROMOS++ follows the **Unix philosophy**:
- Each program does one thing well
- Programs work together via file I/O
- Text-based input/output for portability
- Scriptable for automation

## Directory Structure

```
gromosPlusPlus/
├── programs/                 # 104 executable programs
│   │
│   ├── # === SIMULATION & PREPROCESSING (10 tools) ===
│   ├── pdb2g96.cc           # PDB → GROMOS96 format
│   ├── make_top.cc          # Build topology from building blocks
│   ├── com_top.cc           # Combine topologies
│   ├── con_top.cc           # Concatenate topologies
│   ├── red_top.cc           # Reduce topology (remove atoms)
│   ├── link_top.cc          # Link topologies (covalent bonds)
│   ├── addvirt_top.cc       # Add virtual atoms to topology
│   ├── pert_top.cc          # Create perturbation topology
│   ├── check_top.cc         # Validate topology
│   ├── prep_eds.cc          # Prepare EDS simulation
│   │
│   ├── # === SYSTEM SETUP (8 tools) ===
│   ├── build_conf.cc        # Build coordinates from topology
│   ├── ran_box.cc           # Randomize solvent positions
│   ├── sim_box.cc           # Solvate system in box
│   ├── ion.cc               # Add ions to neutralize/set concentration
│   ├── amber2gromos.cc      # Convert AMBER to GROMOS
│   ├── inbox.cc             # Put atoms inside box
│   ├── unify_box.cc         # Unify box representation
│   ├── bin_box.cc           # Discretize box for analysis
│   │
│   ├── # === ENERGY ANALYSIS (10 tools) ===
│   ├── ene_ana.cc           # Energy analysis (main tool)
│   ├── int_ener.cc          # Interaction energy between groups
│   ├── dg_ener.cc           # Free energy from energy differences
│   ├── edyn.cc              # Energy dynamics
│   ├── m_widom.cc           # Widom insertion (chemical potential)
│   ├── pb_solve.cc          # Poisson-Boltzmann solver
│   ├── dgslv_pbsolv.cc      # Solvation free energy (PB)
│   ├── epsilon.cc           # Dielectric constant
│   ├── bar.cc               # Bennett Acceptance Ratio (FEP)
│   ├── ext_ti_ana.cc        # Thermodynamic integration analysis
│   │
│   ├── # === STRUCTURAL ANALYSIS (18 tools) ===
│   ├── rmsd.cc              # Root mean square deviation
│   ├── rmsf.cc              # Root mean square fluctuation
│   ├── rgyr.cc              # Radius of gyration
│   ├── dssp.cc              # Secondary structure (Kabsch-Sander)
│   ├── sasa.cc              # Solvent accessible surface area
│   ├── sasa_hasel.cc        # SASA (Hasel method)
│   ├── cog.cc               # Center of geometry
│   ├── com.cc               # Center of mass
│   ├── contactnum.cc        # Contact number
│   ├── close_pair.cc        # Close contact analysis
│   ├── cry.cc               # Crystallinity analysis
│   ├── cry_rms.cc           # Crystalline RMSD
│   ├── structurefactor.cc   # Structure factor S(q)
│   ├── bilayer_dist.cc      # Bilayer thickness distribution
│   ├── bilayer_oparam.cc    # Bilayer order parameters
│   ├── atom_info.cc         # Atom information
│   ├── follow.cc            # Follow specific atoms
│   ├── solute_entropy.cc    # Solute entropy estimation
│   │
│   ├── # === INTERACTION ANALYSIS (10 tools) ===
│   ├── hbond.cc             # Hydrogen bond analysis (complex!)
│   ├── rdf.cc               # Radial distribution function
│   ├── dipole.cc            # Dipole moment calculation
│   ├── diffus.cc            # Diffusion coefficient
│   ├── visco.cc             # Viscosity calculation
│   ├── cos_dipole.cc        # Collective dipole correlation
│   ├── cos_epsilon.cc       # Dielectric permittivity
│   ├── tcf.cc               # Time correlation functions
│   ├── ditrans.cc           # Distance-dependent transitions
│   ├── tser.cc              # Transition analysis
│   │
│   ├── # === NMR ANALYSIS (6 tools) ===
│   ├── jval.cc              # J-value calculation
│   ├── noe.cc               # NOE analysis
│   ├── prep_noe.cc          # Prepare NOE restraints
│   ├── nhoparam.cc          # N-H order parameters
│   ├── rdc_ana.cc           # RDC analysis
│   ├── rdc_prep.cc          # Prepare RDC restraints
│   │
│   ├── # === X-RAY CRYSTALLOGRAPHY (8 tools) ===
│   ├── r_factor.cc          # Calculate R-factor
│   ├── r_real_factor.cc     # Real-space R-factor
│   ├── xray2gromos.cc       # Convert X-ray structure
│   ├── prep_xray.cc         # Prepare X-ray restraints
│   ├── bfactor.cc           # B-factor calculation
│   ├── cry_rms.cc           # Crystal RMSD
│   ├── prep_bb.cc           # Prepare backbone restraints
│   ├── rot_rel.cc           # Rotational relaxation
│   │
│   ├── # === TRAJECTORY PROCESSING (12 tools) ===
│   ├── frameout.cc          # Extract trajectory frames
│   ├── filter.cc            # Filter trajectory
│   ├── tstrip.cc            # Strip solvent from trajectory
│   ├── gathtraj.cc          # Concatenate trajectories
│   ├── trs_ana.cc           # Trajectory statistics
│   ├── copy_box.cc          # Copy box vectors
│   ├── check_box.cc         # Check box validity
│   ├── inbox.cc             # Put coordinates in box
│   ├── fit_ener.cc          # Fit energy to trajectory
│   ├── epath.cc             # Energy path analysis
│   ├── distance_filter.cc   # Filter by distance
│   ├── atom_filter.cc       # Filter by atom selection
│   │
│   ├── # === CLUSTERING & CONFORMATIONS (6 tools) ===
│   ├── cluster.cc           # Conformational clustering (RMSD)
│   ├── clust_rmsd.cc        # Clustering by RMSD
│   ├── clust_dist.cc        # Clustering by distance
│   ├── confmat.cc           # Conformation matrix
│   ├── follow.cc            # Follow conformational changes
│   ├── pca.cc               # Principal component analysis
│   │
│   ├── # === FREE ENERGY TOOLS (6 tools) ===
│   ├── bar.cc               # Bennett Acceptance Ratio (main FEP tool)
│   ├── ext_ti_ana.cc        # TI analysis (comprehensive)
│   ├── ext_ti_merge.cc      # Merge TI data
│   ├── eds_update_1.cc      # EDS energy offset update (method 1)
│   ├── eds_update_2.cc      # EDS energy offset update (method 2)
│   ├── dg_ener.cc           # ΔG from energy distributions
│   │
│   └── # === UTILITIES & MISC (10 tools) ===
│       ├── mk_script.cc     # Generate simulation scripts
│       ├── prep_hybrid.cc   # Prepare hybrid topology
│       ├── prep_leap.cc     # Prepare for Leap-Frog
│       ├── prep_md.cc       # Prepare MD input
│       └── ... (6 more)
│
├── src/                     # Source code libraries
│   ├── gcore/               # Core classes (60K lines)
│   │   ├── System.cc        # System representation
│   │   ├── Molecule.cc      # Molecule class
│   │   ├── Atom.cc          # Atom class
│   │   ├── Bond.cc          # Bond, angle, dihedral
│   │   ├── GromosForceField.cc # Force field parameters
│   │   ├── Box.cc           # Periodic box
│   │   └── AtomTopology.cc  # Atom topology
│   │
│   ├── gio/                 # I/O library (40K lines)
│   │   ├── InTopology.cc    # Read .top files
│   │   ├── InG96.cc         # Read GROMOS96 coordinates
│   │   ├── InPDB.cc         # Read PDB files
│   │   ├── InIACElementNameMapping.cc # Element mapping
│   │   ├── OutG96.cc        # Write GROMOS96 format
│   │   ├── OutPDB.cc        # Write PDB format
│   │   └── Ginstream.cc     # GROMOS input stream
│   │
│   ├── gmath/               # Math library (30K lines)
│   │   ├── Vec.cc           # Vector operations
│   │   ├── Matrix.cc        # Matrix operations
│   │   ├── Stat.cc          # Statistics
│   │   ├── Distribution.cc  # Probability distributions
│   │   ├── Correlation.cc   # Correlation functions
│   │   ├── Physics.cc       # Physical constants
│   │   └── Expression.cc    # Math expression parser
│   │
│   ├── args/                # Argument parsing (10K lines)
│   │   ├── Arguments.cc     # Command-line arguments
│   │   ├── BoundaryParser.cc # Boundary condition parsing
│   │   └── GatherParser.cc  # Atom selection parser
│   │
│   ├── bound/               # Boundary conditions (5K lines)
│   │   ├── Boundary.cc      # PBC handling
│   │   ├── RectBox.cc       # Rectangular box
│   │   ├── TriclinicBox.cc  # Triclinic box
│   │   └── Vacuum.cc        # No PBC
│   │
│   ├── fit/                 # Fitting algorithms (8K lines)
│   │   ├── RotationalFit.cc # Least-squares fit
│   │   ├── TranslationalFit.cc # Translation
│   │   ├── Reference.cc     # Reference structure
│   │   └── PositionUtils.cc # Position utilities
│   │
│   ├── utils/               # Utilities (20K lines)
│   │   ├── AtomSpecifier.cc # Atom selection syntax
│   │   ├── Energy.cc        # Energy calculations
│   │   ├── groTime.cc       # Time handling
│   │   ├── Rmsd.cc          # RMSD calculation
│   │   ├── Neighbours.cc    # Neighbor search
│   │   └── Value.cc         # Generic value class
│   │
│   └── pb/                  # Poisson-Boltzmann (15K lines)
│       ├── PB_Parameters.cc # PB parameters
│       ├── FFTGridType.cc   # FFT grid
│       ├── FFTPoisson.cc    # FFT Poisson solver
│       └── PB_Solver.cc     # Main PB solver
│
└── config/                  # Build configuration
    ├── Config.sh           # Configuration script
    └── Makefile.am         # Automake configuration
```

## Key Programs (Detailed)

### 1. Energy Analysis

#### `ene_ana` - Main Energy Analysis Tool

```cpp
// Energy analysis from .tre files
// Most commonly used GROMOS++ tool

// Features:
// - Time series of energy components
// - Statistical analysis (mean, std, block averages)
// - Energy drift detection
// - Distribution analysis
// - Free energy from distributions

// Usage:
// ene_ana @traj ener.tre @prop e_pot e_kin e_tot temp pres
//         @time 0 1000 @dt 2

// Output:
// # Time series
// # time       e_pot       e_kin       e_tot       temp        pres
//   0.000    -45632.1     8234.5    -37397.6     298.3       1.013
//   2.000    -45589.3     8256.3    -37333.0     298.9       0.987
//   ...

// Statistical analysis:
// Property     Mean         Std          Min          Max
// e_pot       -45611.2     234.5       -46234.1    -45102.3
// e_kin         8245.3      89.2         8023.4      8567.2
// temp          298.6        3.2          290.1       309.4
```

**Complex Energy Properties**:
- Total energy (conservative test)
- Potential energy (bonded + nonbonded)
- Kinetic energy (temperature)
- Temperature (instantaneous, average)
- Pressure (virial, kinetic)
- Volume
- Density
- Box dimensions
- Energy components (bond, angle, dihedral, LJ, Coulomb, etc.)
- Energy groups (solute-solute, solute-solvent, etc.)

#### `bar` - Bennett Acceptance Ratio

```cpp
// BAR: Optimal free energy estimator for FEP
// Combines forward and reverse sampling

// Theory:
// ΔG = -kT ln(<exp(-ΔU/kT)>_forward * <exp(ΔU/kT)>_reverse)^(1/2)
//    ≈ ΔU_avg + kT ln(<exp(-(ΔU-ΔU_avg)/kT)>)

// Usage:
// bar @dlg lambda_*.dlg @temp 300 @nblocks 5

// Output:
// Lambda    ΔG (kJ/mol)    Error      Samples
// 0.0-0.1     12.34 ±      0.23       10000
// 0.1-0.2     14.56 ±      0.31       10000
// 0.2-0.3     15.89 ±      0.28       10000
// ...
// Total:      145.67 ±     1.23
```

**Features**:
- Optimal statistical efficiency (better than TI or FEP alone)
- Automatic overlap detection
- Block averaging for error analysis
- MBAR (Multistate BAR) extension
- Publication-ready output

#### `ext_ti_ana` - Thermodynamic Integration Analysis

```cpp
// TI: Alternative to BAR using dH/dλ
// ΔG = ∫₀¹ <dH/dλ>_λ dλ

// Numerical integration methods:
// - Trapezoidal rule
// - Simpson's rule
// - Cubic spline interpolation

// Usage:
// ext_ti_ana @dlg lambda_*.dlg @method simpson @plot

// Output:
// Lambda    <dH/dλ>      Error       ΔG(λ)
// 0.0       156.7 ±      3.2         0.0
// 0.1       142.3 ±      2.8         14.95
// 0.2       128.9 ±      2.6         28.86
// ...
// 1.0        45.2 ±      1.8        145.34

// Plus: Plot data for gnuplot/matplotlib
```

### 2. Structural Analysis

#### `rmsd` - Root Mean Square Deviation

```cpp
// RMSD calculation with various options

class RMSD_Calculator {
    // Atom selection for RMSD calculation
    utils::AtomSpecifier atoms_fit;   // Atoms for fitting
    utils::AtomSpecifier atoms_calc;  // Atoms for RMSD calculation

    // Reference structure
    gcore::System ref_system;

    // Fitting options
    bool do_fit = true;           // Fit before RMSD?
    bool mass_weighted = true;    // Mass-weighted RMSD?
    bool pbc = true;              // Consider PBC?

    double calculate_rmsd(gcore::System &sys) {
        // 1. Optional: Fit structures
        if (do_fit) {
            fit::RotationalFit fit(atoms_fit);
            fit.fit(&sys);
        }

        // 2. Calculate RMSD
        double sum = 0.0;
        double total_mass = 0.0;

        for (int i = 0; i < atoms_calc.size(); i++) {
            Vec r_ref = ref_system.atom(atoms_calc.atom(i)).pos();
            Vec r_cur = sys.atom(atoms_calc.atom(i)).pos();

            // Apply PBC if needed
            if (pbc) {
                r_cur = pbc_nearest(r_cur, r_ref);
            }

            double dist2 = (r_ref - r_cur).abs2();

            if (mass_weighted) {
                double mass = sys.atom(atoms_calc.atom(i)).mass();
                sum += mass * dist2;
                total_mass += mass;
            } else {
                sum += dist2;
                total_mass += 1.0;
            }
        }

        return sqrt(sum / total_mass);
    }
};

// Usage examples:
// rmsd @topo system.top @traj traj.trc @ref crystal.pdb
//      @atomsfit 1:CA                  # Fit on Cα atoms
//      @atomsrmsd 1:N,CA,C,O           # RMSD on backbone
//      @pbc r                          # Use PBC

// Output:
// # Time(ps)    RMSD(nm)
//   0.000       0.000
//   2.000       0.123
//   4.000       0.145
//   ...
```

**Advanced RMSD Features**:
- Per-residue RMSD
- RMSD matrix (all vs all frames)
- Cluster representatives
- Time-averaged structures
- Crystal symmetry handling

#### `dssp` - Secondary Structure (Kabsch-Sander)

```cpp
// DSSP algorithm for secondary structure assignment
// H-bond patterns determine structure

// Structure types:
// H = α-helix
// E = β-sheet (extended)
// T = turn
// S = bend
// G = 3₁₀-helix
// I = π-helix
// B = β-bridge
// C = coil (random)

// H-bond criterion:
// E_HB = 0.084 * [1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN] kcal/mol
// H-bond if E_HB < -0.5 kcal/mol

// Usage:
// dssp @topo system.top @traj traj.trc @out dssp.dat

// Output:
// # Residue    Time...
// #            0       2       4       6       8
//   1(MET)     C       C       C       C       C
//   2(GLY)     C       C       C       C       C
//   3(SER)     C       C       H       H       H
//   4(ALA)     H       H       H       H       H
//   5(LEU)     H       H       H       H       H
//   ...

// Summary:
// Time(ps)    Helix%    Sheet%    Turn%    Coil%
// 0.0         32.1      18.4      12.3     37.2
// 2.0         33.5      17.9      11.8     36.8
// ...
```

#### `sasa` - Solvent Accessible Surface Area

```cpp
// SASA calculation using probe sphere
// Two implementations: SASA (exact) and SASA_hasel (fast approximation)

// Shrake-Rupley algorithm:
// - Place probe sphere (r=0.14 nm) on molecular surface
// - Check if probe overlaps with other atoms
// - Surface points accessible = SASA

class SASA_Calculator {
    double probe_radius = 0.14;  // Water radius (nm)
    int n_points = 960;          // Surface points per atom (Fibonacci sphere)

    double calculate_sasa(gcore::System &sys,
                          utils::AtomSpecifier &atoms) {
        double total_sasa = 0.0;

        for (int i = 0; i < atoms.size(); i++) {
            Atom &atom = sys.atom(atoms.atom(i));
            double r_i = atom.radius() + probe_radius;
            Vec pos_i = atom.pos();

            int accessible_points = 0;

            // Generate Fibonacci sphere points
            for (int p = 0; p < n_points; p++) {
                Vec surface_point = pos_i + r_i * fibonacci_point(p, n_points);

                // Check if this point is buried by other atoms
                bool buried = false;
                for (int j = 0; j < sys.numAtoms(); j++) {
                    if (i == j) continue;

                    Atom &atom_j = sys.atom(j);
                    double r_j = atom_j.radius() + probe_radius;
                    Vec pos_j = atom_j.pos();

                    if ((surface_point - pos_j).abs() < r_j) {
                        buried = true;
                        break;
                    }
                }

                if (!buried) accessible_points++;
            }

            // Area of sphere: 4πr²
            double sphere_area = 4.0 * M_PI * r_i * r_i;
            double atom_sasa = sphere_area * accessible_points / n_points;
            total_sasa += atom_sasa;
        }

        return total_sasa;
    }
};

// Usage:
// sasa @topo system.top @traj traj.trc
//      @atoms 1:PRO              # Protein only
//      @probe 0.14               # Water probe radius

// Output:
// # Time(ps)    SASA(nm²)    SASA_hydrophobic    SASA_polar
//   0.0         165.3        98.2                67.1
//   2.0         167.1        99.5                67.6
//   ...

// Per-residue:
// Residue    SASA(nm²)    Buried%
// 1(MET)     0.234        78.3
// 2(GLY)     0.156        34.2
// ...
```

### 3. Interaction Analysis

#### `hbond` - Hydrogen Bond Analysis

One of the most complex GROMOS++ programs!

```cpp
// H-bond criteria:
// - Distance: r(D...A) < 0.25 nm
// - Angle: ∠(D-H...A) > 135°
// - Optionally: ∠(H...A-X) where X is atom bonded to A

class HBond_Analyzer {
    // Criteria
    double max_distance = 0.25;  // nm
    double min_angle = 135.0;    // degrees
    bool three_centered = false; // D-H...A-X angle

    // Donor/acceptor definitions
    std::vector<HBond_Donor> donors;
    std::vector<HBond_Acceptor> acceptors;

    struct HBond {
        int donor_id;
        int hydrogen_id;
        int acceptor_id;
        double distance;
        double angle_DHA;
        double angle_HAX;
        int first_frame;
        int last_frame;
        double occupancy;
    };

    void analyze_trajectory() {
        std::vector<HBond> hbonds;

        // For each frame
        for (Frame &frame : trajectory) {
            // Check all donor-acceptor pairs
            for (auto &donor : donors) {
                for (auto &acceptor : acceptors) {
                    if (is_hbond(donor, acceptor, frame)) {
                        update_hbond_statistics(hbonds, donor, acceptor);
                    }
                }
            }
        }

        // Post-processing
        calculate_occupancies(hbonds);
        identify_persistent_hbonds(hbonds);
        calculate_hbond_lifetimes(hbonds);
        calculate_hbond_dynamics(hbonds);
    }
};

// Usage:
// hbond @topo system.top @traj traj.trc
//       @donors 1:PRO@N,O        # Protein N,O as donors
//       @acceptors s:O           # Solvent O as acceptors
//       @distcut 0.25            # Distance cutoff
//       @anglecut 135            # Angle cutoff

// Output:
// # Hydrogen Bonds
// Donor           Acceptor        Occupancy    AvgDist    AvgAngle    Lifetime
// PRO:5:N-H       PRO:1:O         0.856        0.195      162.3       15.2ps
// PRO:12:N-H      PRO:8:O         0.743        0.203      158.7       8.4ps
// PRO:23:OG-H     WAT:1234:O      0.234        0.215      148.2       2.1ps
// ...

// Statistics:
// Total H-bonds:              156
// Persistent (>50%):           42
// Transient (<50%):           114
// Average per frame:           87.3
// Intramolecular:              34
// Intermolecular:             122
```

**Advanced H-bond Features**:
- Time series of individual H-bonds
- H-bond lifetime distributions
- Autocorrelation functions
- Cluster analysis of H-bond networks
- Visualization output (PyMOL, VMD)

#### `rdf` - Radial Distribution Function

```cpp
// g(r) = <ρ(r)> / <ρ_bulk>
// Probability of finding atom at distance r

class RDF_Calculator {
    double r_min = 0.0;
    double r_max = 1.5;  // nm
    double bin_width = 0.002;  // nm
    int n_bins;

    std::vector<double> histogram;
    std::vector<double> rdf;

    void calculate(AtomSpecifier &atoms_i,
                   AtomSpecifier &atoms_j,
                   Trajectory &traj) {
        n_bins = (r_max - r_min) / bin_width;
        histogram.resize(n_bins, 0.0);

        // Accumulate histogram
        for (Frame &frame : traj) {
            for (Atom &atom_i : atoms_i) {
                for (Atom &atom_j : atoms_j) {
                    if (&atom_i == &atom_j) continue;

                    Vec r_ij = pbc(atom_j.pos() - atom_i.pos());
                    double r = r_ij.abs();

                    if (r >= r_min && r < r_max) {
                        int bin = (r - r_min) / bin_width;
                        histogram[bin] += 1.0;
                    }
                }
            }
        }

        // Normalize
        double rho = atoms_j.size() / frame.box.volume();
        int n_frames = traj.size();

        for (int i = 0; i < n_bins; i++) {
            double r = r_min + (i + 0.5) * bin_width;
            double shell_volume = 4.0 * M_PI * r * r * bin_width;
            double ideal_count = rho * shell_volume * atoms_i.size() * n_frames;

            rdf[i] = histogram[i] / ideal_count;
        }
    }
};

// Usage:
// rdf @topo system.top @traj traj.trc
//     @atoms_i 1:OW              # Water oxygen
//     @atoms_j 1:OW              # Water oxygen
//     @cut 1.5                   # Maximum distance

// Output:
// # r(nm)    g(r)       n(r)      # n(r) = integrated number
//   0.001    0.000      0.000
//   0.003    0.000      0.000
//   0.005    0.002      0.001
//   ...
//   0.275    2.745      4.2       # First peak (H-bonded)
//   ...
//   0.485    0.823      12.1      # Second peak
//   ...
//   1.500    1.002      56.3      # Bulk density
```

**RDF Applications**:
- Solvation structure
- Coordination numbers
- Phase transitions
- Liquid structure

### 4. NMR Analysis

#### `jval` - J-value Calculation

```cpp
// J-value (scalar coupling) from Karplus equation
// J(φ) = A*cos²(φ) + B*cos(φ) + C

// Different Karplus parameters for different couplings:
// ³J(HN-Hα): A=6.4, B=-1.4, C=1.9 (Bystrov)
// ³J(Hα-Hβ): A=9.5, B=-1.6, C=1.8 (DeMarco)

class Jvalue_Calculator {
    struct KarplusParams {
        double A, B, C;
        std::string name;
    };

    std::map<std::string, KarplusParams> karplus_sets = {
        {"HN-HA", {6.4, -1.4, 1.9, "Bystrov"}},
        {"HA-HB", {9.5, -1.6, 1.8, "DeMarco"}},
        {"C-C",   {4.9, -0.6, 0.6, "Perez"}},
        // ... many more
    };

    double calculate_jvalue(double phi, const KarplusParams &params) {
        return params.A * cos(phi) * cos(phi) +
               params.B * cos(phi) +
               params.C;
    }

    void analyze_trajectory(Trajectory &traj,
                            std::vector<Dihedral> &dihedrals) {
        for (auto &dihedral : dihedrals) {
            std::vector<double> jvalues;

            for (Frame &frame : traj) {
                double phi = dihedral.angle(frame);
                double j = calculate_jvalue(phi, karplus_sets[dihedral.type]);
                jvalues.push_back(j);
            }

            // Statistics
            double j_avg = mean(jvalues);
            double j_std = std_dev(jvalues);
            double j_r3avg = pow(mean(pow(jvalues, -3)), -1.0/3.0);  // r⁻³ averaging

            results[dihedral.name] = {j_avg, j_std, j_r3avg};
        }
    }
};

// Usage:
// jval @topo system.top @traj traj.trc
//      @dih phi.dat              # Dihedral definitions
//      @karplus HN-HA            # Karplus parameter set

// Output:
// Residue    φ(avg)      J(avg)      J(r⁻³avg)    J(exp)      Δ
// 5(ALA)     -65.3       8.2         8.4          8.1         +0.3
// 6(LEU)     -72.1       9.1         9.3          9.5         -0.2
// ...

// RMSD from experiment: 0.42 Hz
```

### 5. Trajectory Processing

#### `frameout` - Extract Frames

```cpp
// Extract specific frames or time ranges from trajectory

// Usage modes:
// 1. Every Nth frame
// frameout @traj traj.trc @every 10 @out subset.g96

// 2. Specific time range
// frameout @traj traj.trc @time 100 500 @out range.g96

// 3. Specific frames
// frameout @traj traj.trc @frames 0,10,50,100 @out frames.g96

// 4. Nearest to reference (lowest RMSD)
// frameout @traj traj.trc @ref crystal.pdb @closest 10 @out closest.g96
```

#### `filter` - Filter Trajectories

```cpp
// Filter trajectory based on various criteria

// 1. Energy filter
// filter @traj traj.trc @ene ener.tre @prop e_pot @range -50000 -45000

// 2. Distance filter
// filter @traj traj.trc @dist "1:CA 2:CA" @range 0.3 0.5

// 3. RMSD filter
// filter @traj traj.trc @ref crystal.pdb @rmsd @range 0 0.2

// 4. Secondary structure filter
// filter @traj traj.trc @ss 5:10 @type helix
```

### 6. Free Energy Tools

All free energy methods implemented:

1. **BAR** (Bennett Acceptance Ratio) - Optimal
2. **TI** (Thermodynamic Integration) - Traditional
3. **FEP** (Forward/Reverse) - Simple but biased
4. **MBAR** (Multistate BAR) - Multiple λ windows
5. **WHAM** (Weighted Histogram) - For umbrella sampling
6. **Zwanzig** (Exponential averaging) - Fast but high variance

## Atom Selection Syntax

GROMOS++ uses a powerful atom selection language:

```
# Basic selections
1:CA                  # Atom CA in molecule 1
1:1-10                # Atoms 1-10 in molecule 1
1:ALA,GLY             # All ALA and GLY residues in molecule 1
s:O                   # All O atoms in solvent (s)

# Combinations
1:PRO@CA              # All CA atoms in protein
1:1-100&CA,CB         # CA and CB in residues 1-100
1:ALA|GLY             # ALA or GLY

# Atom types
1:a%OW                # All atoms of type OW
1:a%C*                # All atoms starting with C

# Distances
1:CA<0.5:s:O          # Protein CA within 0.5 nm of water O

# Expressions
1:(res%1-10&a%CA)     # Complex selections with grouping
```

## Common Workflows

### 1. Basic Trajectory Analysis

```bash
# Extract properties
ene_ana @traj energy.tre @prop e_pot temp pres @time 0 1000

# RMSD
rmsd @topo system.top @traj traj.trc @ref crystal.pdb \\
     @atomsfit 1:CA @atomsrmsd 1:N,CA,C,O

# RMSF
rmsf @topo system.top @traj traj.trc @atomsrmsd 1:CA

# Secondary structure
dssp @topo system.top @traj traj.trc @out dssp.dat

# H-bonds
hbond @topo system.top @traj traj.trc \\
      @donors 1:PRO@N,O @acceptors s:O
```

### 2. Free Energy Analysis

```bash
# Run FEP simulations at multiple λ
for lambda in 0.0 0.1 0.2 ... 1.0; do
    md++ @f md.imd @lambda $lambda
done

# BAR analysis
bar @dlg lambda_*.dlg @temp 300 @out free_energy.dat

# TI analysis
ext_ti_ana @dlg lambda_*.dlg @method simpson

# Compare methods
# BAR typically has lowest error
```

### 3. NMR Structure Refinement

```bash
# Calculate NOE violations
noe @topo system.top @traj traj.trc @noe noe_restraints.dat

# Calculate J-values
jval @topo system.top @traj traj.trc @dih phi_psi.dat

# RDC calculations
rdc_ana @topo system.top @traj traj.trc @rdc rdc_data.dat

# Overall quality
# R-factor, violations, RMSD to experimental
```

## Performance

GROMOS++ tools are generally very fast:

| Program | System Size | Performance |
|---------|-------------|-------------|
| `ene_ana` | Any | ~1M frames/min |
| `rmsd` | 2.6K atoms | ~10K frames/min |
| `rmsd` | 23K atoms | ~2K frames/min |
| `hbond` | 2.6K atoms | ~500 frames/min |
| `rdf` | 23K atoms | ~1K frames/min |
| `dssp` | 100 residues | ~5K frames/min |
| `bar` | 10 windows | <1 second |

## Summary

### Strengths

1. **Comprehensive**: 104 tools covering all MD analysis needs
2. **Well-tested**: 20+ years of use in publications
3. **Scriptable**: Unix philosophy, easy automation
4. **Fast**: Efficient C++ implementation
5. **Flexible**: Powerful atom selection syntax
6. **Free energy**: Complete set of FEP/TI tools
7. **NMR/X-ray**: Specialized tools for structure refinement

### Weaknesses

1. **Learning curve**: 104 programs to learn
2. **Documentation**: Scattered across 104 man pages
3. **Inconsistency**: Different programs, different interfaces
4. **No GUI**: Command-line only
5. **Limited visualization**: Text output, need external tools

### Integration with gromos-rs

**Recommended Strategy**:
- **gromos-rs**: Run simulations (2-3x faster)
- **gromos++**: Analyze trajectories (104 tools)

**Workflow**:
```bash
# Simulate with gromos-rs
gromos-rs md --traj output.trc --energy output.tre

# Analyze with gromos++
ene_ana @traj output.tre @prop e_tot
rmsd @traj output.trc @ref crystal.pdb
hbond @traj output.trc @donors 1:PRO@N,O
```

**Why not reimplement?**
- 250K lines of battle-tested code
- 20+ years of development
- Widely used, validated
- Focus gromos-rs on simulation performance

---

**Next**: [GROMOS-RS Architecture](gromos-rs.md) | [Comparison](comparison.md) | [Gap Analysis](gaps.md)
