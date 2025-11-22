# Quick Start Guide

Get started with GROMOS-RS in minutes!

## Your First Simulation

### 1. Prepare Input Files

You'll need three files:

- **Topology** (`.top`): Molecular structure
- **Coordinates** (`.cnf`): Initial positions/velocities
- **Parameters** (`.imd`): Simulation settings

Example directory structure:

```
my_simulation/
├── system.top       # Topology
├── initial.cnf      # Starting coordinates
└── md.imd          # MD parameters
```

### 2. Run MD Simulation

```bash
cd my_simulation

# Basic MD run
gromos-rs md \
    --top system.top \
    --coord initial.cnf \
    --param md.imd \
    --traj output.trc \
    --energy energy.tre
```

Or using the shorter form:

```bash
gromos-rs md -t system.top -c initial.cnf -p md.imd -o output.trc
```

### 3. Monitor Progress

GROMOS-RS writes to stdout:

```
GROMOS-RS v0.3.0 - High-Performance Molecular Dynamics
========================================================
Reading topology: system.top
  Atoms: 2,641
  Molecules: 1 protein + 832 waters
  Bonds: 2,584 (SHAKE constrained)

Reading coordinates: initial.cnf
  Initial T = 298.5 K

Reading parameters: md.imd
  Integrator: Leap-Frog
  Timestep: 2 fs
  Steps: 500,000 (1 ns total)
  Thermostat: Berendsen (τ=0.1 ps)
  Barostat: Berendsen (τ=0.5 ps)

Starting simulation...
  Step    Time(ps)    E_pot      E_kin      E_tot      Temp(K)    Press(bar)
     0       0.000  -12543.2    3214.5   -9328.7      298.5       1.02
  1000       2.000  -12487.6    3245.1   -9242.5      301.2       0.98
  2000       4.000  -12501.3    3228.7   -9272.6      299.7       1.01
  ...

Performance: 95.3 ns/day
Simulation complete!
```

### 4. Analyze Results

Use GROMOS++ tools for analysis:

```bash
# Energy analysis
ener_ana @traj energy.tre @prop e_tot

# RMSD
rmsd @topo system.top @traj output.trc @ref initial.cnf

# Hydrogen bonds
hbond @topo system.top @traj output.trc
```

## Common Workflows

### NVT Simulation (Constant Temperature)

Create `nvt.imd`:

```
# NVT ensemble
SYSTEM
    NPM      1     # 1 molecule
    NSM      832   # 832 solvent molecules
END
INITIALISE
    TEMPI    300.0  # Initial temperature
END
STEP
    NSTLIM   500000 # 500k steps = 1 ns
    DT       0.002  # 2 fs timestep
END
FORCE
    BONDS    0      # SHAKE all bonds
    ANGLES   1      # Include angles
    IMPROPER 1      # Include impropers
    DIHEDRAL 1      # Include dihedrals
    NTF      1      # Force calculation type
END
PAIRLIST
    ALGORITHM 2     # Grid cell pairlist
    RCUTP    1.4    # Short-range cutoff
    RCUTL    1.4    # Long-range cutoff
END
NONBONDED
    NLRELE   1      # Reaction field
    EPSRF    61.0   # RF permittivity
END
CONSTRAINT
    NTC      3      # SHAKE all bonds
    NTCP     1      # Parallel SHAKE
END
TEMPERATURE
    NTTB     1      # Berendsen thermostat
    TEMP0    300.0  # Target temperature
    TAUT     0.1    # Coupling time
END
```

Run:

```bash
gromos-rs md -t system.top -c initial.cnf -p nvt.imd -o nvt.trc
```

### NPT Simulation (Constant Pressure & Temperature)

Add to `npt.imd`:

```
PRESSURE
    NPTB     1      # Berendsen barostat
    PRES0    1.0    # 1 bar target
    COMP     4.575e-4  # Compressibility of water
    TAUP     0.5    # Coupling time
END
```

### Free Energy Perturbation (FEP)

#### Step 1: Create Perturbation Topology

```bash
# Generate .ptp file defining state A → state B transformation
gromos-rs make_pt_top \
    --top system.top \
    --ptp perturbation.ptp \
    --state-a "resid 1 and name CA" \
    --state-b "resid 1 and name CB"
```

#### Step 2: Run FEP Windows

```bash
# Run multiple λ windows
for lambda in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
    gromos-rs md \
        -t system.top \
        -c initial.cnf \
        -p fep.imd \
        --ptp perturbation.ptp \
        --lambda $lambda \
        --dlg fep_${lambda}.dlg \
        -o fep_${lambda}.trc
done
```

#### Step 3: Analyze Free Energy

```bash
# Use gromos++ BAR (Bennett Acceptance Ratio)
bar @dlg fep_*.dlg @temp 300 > free_energy.dat
```

### Replica Exchange MD (REMD)

```bash
# Temperature REMD with 16 replicas
gromos-rs remd \
    --top system.top \
    --coord initial.cnf \
    --param remd.imd \
    --replicas 16 \
    --temp-range 300:450 \
    --exchange-freq 1000 \
    --output-dir remd_output/
```

Analysis:

```bash
rep_ana @dir remd_output/ @temps 300:450
```

### Enhanced Sampling: EDS

```bash
# Enveloping Distribution Sampling
gromos-rs eds \
    --top system.top \
    --coord initial.cnf \
    --param eds.imd \
    --states 4 \
    --eoff-init "0.0 10.0 20.0 30.0" \
    --output eds_output.trc
```

### Enhanced Sampling: GaMD

```bash
# Gaussian Accelerated MD
gromos-rs gamd \
    --top system.top \
    --coord initial.cnf \
    --param gamd.imd \
    --mode dual \
    --search-mode boost \
    --output gamd_output.trc
```

## Performance Tips

### Use Native CPU Optimizations

Build with:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

This enables AVX2/AVX-512 for 2-3x speedup.

### Choose Efficient Pairlist Algorithm

In `.imd` file:

```
PAIRLIST
    ALGORITHM 2  # Grid cell pairlist (O(N), fastest for large systems)
END
```

Options:
- `1`: Standard (O(N²), small systems <5K atoms)
- `2`: Grid cell (O(N), large systems >5K atoms) **← Recommended**

### Use Parallel SHAKE

```
CONSTRAINT
    NTCP 1  # Parallel SHAKE
END
```

### Enable Reaction Field (Faster than PME for Most Cases)

```
NONBONDED
    NLRELE 1      # Reaction field
    EPSRF  61.0   # Permittivity (water)
END
```

Use PME only for highly charged systems:

```
NONBONDED
    NLRELE 2      # PME
    KAPPA  0.33   # Ewald parameter
    NGRID  64 64 64  # FFT grid
END
```

### Use Multiple Threads

```bash
# Set thread count
export RAYON_NUM_THREADS=16

gromos-rs md -t system.top -c initial.cnf -p md.imd -o output.trc
```

## Example Systems

### Lysozyme in Water (Small, ~2.6K atoms)

```bash
# Download example
wget https://www.gromos.net/examples/lysozyme.tar.gz
tar xzf lysozyme.tar.gz
cd lysozyme

# Run NVT equilibration (100 ps)
gromos-rs md -t lyso.top -c lyso.cnf -p nvt.imd -o nvt.trc

# Run NPT production (1 ns)
gromos-rs md -t lyso.top -c nvt_final.cnf -p npt.imd -o prod.trc
```

Expected performance: **~200 ns/day** on modern CPU (16 cores)

### DHFR in Water (Medium, ~23.5K atoms)

Expected performance: **~95 ns/day** on modern CPU (16 cores)

### Membrane Protein (Large, ~85K atoms)

Expected performance: **~30 ns/day** on modern CPU (16 cores)

## Validation

### Energy Conservation (NVE)

Run NVE (no thermostat/barostat) and check energy drift:

```bash
gromos-rs md -t system.top -c initial.cnf -p nve.imd -o nve.trc

# Analyze
ener_ana @traj nve_energy.tre @prop e_tot

# Should show <0.01% drift over 1 ns
```

### Temperature Distribution (NVT)

```bash
ener_ana @traj energy.tre @prop temp

# Should show mean ≈ 300 K, σ ≈ 5-10 K
```

### Density (NPT)

```bash
ener_ana @traj energy.tre @prop density

# Water should be ~1.0 g/cm³ at 300 K, 1 bar
```

## Next Steps

- [User Guide](../user-guide/introduction.md) - Detailed documentation
- [Input Files](../user-guide/input-files.md) - File format reference
- [Free Energy](../user-guide/free-energy.md) - Advanced FEP guide
- [Analysis](../user-guide/analysis.md) - Post-processing workflows

## Troubleshooting

### Simulation Crashes

Check:
1. Topology is valid: `check_top @topo system.top`
2. Coordinates match topology: `check_conf @topo system.top @conf initial.cnf`
3. Timestep not too large (should be ≤2 fs for constrained bonds)

### Poor Performance

1. Use grid cell pairlist (`ALGORITHM 2`)
2. Build with `target-cpu=native`
3. Set `RAYON_NUM_THREADS` to physical core count
4. Use Reaction Field instead of PME if possible

### Energy Not Conserved

1. Reduce timestep (try 1 fs)
2. Tighten SHAKE tolerance
3. Check for clashes in initial structure

## Getting Help

- **Documentation**: [Full docs](../user-guide/introduction.md)
- **Issues**: [GitHub](https://github.com/yerkoescalona/gromos-rs/issues)
- **Forum**: [GROMOS Forum](https://www.gromos.net/forum)
