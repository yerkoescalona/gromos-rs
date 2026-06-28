# API Reference

All classes and functions exported by `import gromos`.

---

## System

Pairs a `Topology` with a `Configuration`. Validates that both describe the same
number of atoms at construction time.

```python
from gromos import System
```

### Constructors

```python
System(topology: Topology, configuration: Configuration) -> System
```
Raises `ValueError` if `topology.n_atoms != configuration.n_atoms`.

```python
System.from_files(topo_file: str, conf_file: str) -> System  # staticmethod
```
Reads both files and validates in one call.

**Example**
```python
system = System.from_files("water_216.topo", "equilibrated.cnf")

topo   = Topology("water_216.topo")
conf   = Configuration("equilibrated.cnf")
system = System(topo, conf)   # identical result
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_atoms` | `int` | Total atom count |
| `charge` | `int` | Net integer charge (e) |
| `positions` | `ndarray (N,3) f64` | Positions, nm |
| `velocities` | `ndarray (N,3) f64` | Velocities, nm/ps |
| `box` | `tuple[float,float,float]` | Box dimensions (Lx, Ly, Lz), nm |
| `topology` | `Topology` | Underlying topology |
| `configuration` | `Configuration` | Underlying configuration |

### Methods

**`write(path: str)`**  
Write current coordinates to a GROMOS `.cnf` file.

---

## Topology

Reads a GROMOS topology file (`.topo`).

```python
from gromos import Topology
topo = Topology("system.topo")
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_atoms` | `int` | Total atom count (solute + solvent) |
| `n_solute_atoms` | `int` | Solute-only count |
| `n_solvent_atoms` | `int` | Solvent-only count |
| `masses` | `ndarray (N,) f64` | Atom masses, g/mol |
| `charges` | `ndarray (N,) f64` | Partial charges, e |

### Methods

**`solvate(nsm: int)`**  
Append `nsm` copies of the solvent template to the topology in-place.

---

## Configuration

Reads a GROMOS coordinate file (`.cnf` / `.g96`).

```python
from gromos import Configuration
conf = Configuration("initial.cnf")
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_atoms` | `int` | Atom count in this file |
| `positions` | `ndarray (N,3) f64` | Positions, nm |
| `velocities` | `ndarray (N,3) f64` | Velocities, nm/ps |
| `box_dimensions` | `tuple[float,float,float]` | (Lx, Ly, Lz), nm |

---

## InputParameters

MD / minimisation control parameters. Two construction paths: load an existing
`.imd` file, or use a factory that sets sensible defaults for the chosen ensemble.

```python
from gromos import InputParameters
```

### Constructors

```python
InputParameters(input_file: str)             # load from file
InputParameters.from_file(input_file: str)   # identical staticmethod alias
```

### Factory methods (staticmethod)

```python
InputParameters.nve(dt: float, steps: int) -> InputParameters
```
Microcanonical ensemble. Thermostat coupling time set to −1 (no coupling).

```python
InputParameters.nvt(dt: float, steps: int, temperature: float) -> InputParameters
```
Canonical ensemble. Berendsen thermostat, τ = 0.1 ps.

```python
InputParameters.npt(dt: float, steps: int, temperature: float, pressure: float) -> InputParameters
```
Isothermal-isobaric ensemble. Berendsen thermostat + barostat.  
`pressure` in bar. Compressibility defaults to 4.575 × 10⁻⁴ nm² kJ⁻¹ mol (water).

```python
InputParameters.steepest_descent(steps: int) -> InputParameters
```
Steepest-descent energy minimisation (`ENERGYMIN` block, `ntem=1`).

**Example**
```python
params = InputParameters.nvt(dt=0.002, steps=5000, temperature=300.0)
print(params.dt, params.nstlim, params.temperature)
```

### Properties (read-only)

| Property | Type | Description |
|----------|------|-------------|
| `dt` | `float` | Timestep, ps |
| `nstlim` | `int` | Number of steps |
| `temperature` | `float` | First bath target temperature, K |
| `cutoff` | `float` | Long-range cutoff (rcutl), nm |
| `rcutp` | `float` | Short-range pairlist cutoff, nm |
| `nsm` | `int` | Number of solvent molecules |
| `ntc` | `int` | SHAKE mode (1=none, 2=H-bonds, 3=all) |
| `ntb` | `int` | Boundary type (0=vacuum, 1=rectangular) |
| `nsnb` | `int` | Pairlist update frequency |
| `ntwx` | `int` | Trajectory write frequency |
| `ntwe` | `int` | Energy write frequency |

---

## Simulation

Builds an algorithm sequence from the parameters and runs the MD loop.

```python
from gromos import Simulation
```

### Constructors

```python
# Recommended — two-argument form
Simulation(system: System, params: InputParameters)

# Three-argument forms (legacy, still supported)
Simulation(topo: Topology, conf: Configuration, params: InputParameters)
Simulation(topo_file: str, conf_file: str, input_file: str)

# Explicit file-path staticmethod
Simulation.from_files(topo_file: str, conf_file: str, input_file: str)

# Custom algorithm sequence
Simulation.from_sequence(topo: Topology, conf: Configuration,
                         params: InputParameters, sequence: AlgorithmSequence)
```

### Running

**`step(n_steps: int)`**  
Advance the simulation by `n_steps` steps. All state properties are updated after
this call.

### State properties (read-only)

| Property | Type | Description |
|----------|------|-------------|
| `total_energy` | `float` | Total energy, kJ/mol |
| `kinetic_energy` | `float` | Kinetic energy, kJ/mol |
| `potential_energy` | `float` | Potential energy, kJ/mol |
| `temperature` | `float` | Instantaneous temperature, K |
| `positions` | `ndarray (N,3) f64` | Positions, nm |
| `velocities` | `ndarray (N,3) f64` | Velocities, nm/ps |
| `forces` | `ndarray (N,3) f64` | Forces, kJ/(mol·nm) |
| `time` | `float` | Current time, ps |
| `current_step` | `int` | Step counter |
| `dt` | `float` | Timestep, ps |
| `n_atoms` | `int` | Atom count |
| `n_solute_atoms` | `int` | Solute atom count |
| `n_solvent_atoms` | `int` | Solvent atom count |
| `algorithm_names` | `list[str]` | Names of algorithms in the sequence |
| `energies` | `Energy` | Full energy object |

---

## AlgorithmSequence

The ordered list of algorithms executed each MD step. Allows inspection and
modification of the pipeline before constructing a `Simulation`.

```python
from gromos import AlgorithmSequence
```

### Factory methods (staticmethod)

```python
AlgorithmSequence.nve(topo: Topology, params: InputParameters) -> AlgorithmSequence
AlgorithmSequence.nvt(topo: Topology, params: InputParameters) -> AlgorithmSequence
AlgorithmSequence.npt(topo: Topology, params: InputParameters) -> AlgorithmSequence
AlgorithmSequence.from_parameters(topo: Topology, params: InputParameters) -> AlgorithmSequence
```

### Modification methods

```python
seq.add(algorithm)                          # append
seq.insert_after(name: str, algorithm)      # insert after named step
seq.insert_before(name: str, algorithm)     # insert before named step
seq.remove(name: str)                       # remove by name
seq.replace(name: str, algorithm)           # swap by name
```

### Inspection

```python
seq.names          # list[str] — ordered algorithm names
len(seq)           # int
"Forcefield" in seq  # bool
```

**Example**
```python
seq = AlgorithmSequence.nvt(topo, params)
print(seq.names)
# ['RemoveCOMMotion', 'Forcefield', 'LeapFrogVelocity', 'BerendsenThermostat',
#  'LeapFrogPosition', 'TemperatureCalculation', 'EnergyCalculation']

seq.remove("RemoveCOMMotion")
sim = Simulation.from_sequence(topo, conf, params, seq)
```

---

## Algorithm building blocks

Individual algorithm objects that can be inserted into an `AlgorithmSequence`.

### `Forcefield`

Nonbonded force calculation (pairlist + LJ + reaction-field electrostatics + bonded terms).

```python
Forcefield(
    cutoff: float | None = None,
    rcutp: float | None = None,
    epsilon_rf: float | None = None,
    kappa: float | None = None,
    pairlist_update: int | None = None,
    virial: str | None = None,
    ntf_bond: bool | None = None,
    ntf_angle: bool | None = None,
    ntf_improper: bool | None = None,
    ntf_dihedral: bool | None = None,
)
```

`None` values inherit from `InputParameters`.

### `LeapFrogVelocity` / `LeapFrogPosition`

Half-step velocity update and full-step position update (leap-frog integrator).

```python
LeapFrogVelocity()
LeapFrogPosition()
```

### `BerendsenThermostat`

Velocity rescaling thermostat.

```python
BerendsenThermostat(temperature: float = 300.0, tau: float = 0.1)
```

### `BerendsenBarostat`

Isotropic pressure coupling.

```python
BerendsenBarostat(pressure: float = 1.0, tau: float = 0.5,
                  compressibility: float = 4.575e-4, virial: str | None = None)
```

### `ShakeConstraints`

Bond-length constraints via SHAKE.

```python
ShakeConstraints(tolerance: float = 1e-4, max_iterations: int = 1000, mode: str = "solute")
```

### `TemperatureCalculation` / `PressureCalculation` / `EnergyCalculation`

Accumulate kinetic/virial/potential energy into the `Energy` object each step.

```python
TemperatureCalculation()
PressureCalculation(virial: str = "atomic")
EnergyCalculation()
```

### `RemoveCOMMotion`

Remove translational (and optionally rotational) centre-of-mass velocity.

```python
RemoveCOMMotion(initial: bool = True, nscm: int = 0)
```

---

## Vec3

A 3D vector. Useful for single-atom inspection; use NumPy arrays from `Simulation`
or `System` for bulk work.

```python
from gromos import Vec3
v = Vec3(1.0, 2.0, 3.0)
```

### Properties

`x`, `y`, `z` — `float` (read-only)

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `length()` | `float` | Euclidean norm |
| `normalize()` | `Vec3` | Unit vector |
| `dot(other: Vec3)` | `float` | Dot product |
| `cross(other: Vec3)` | `Vec3` | Cross product |

### Operators

`+`, `-`, `* scalar` between `Vec3` objects.

---

## Energy

Snapshot of energy components (kJ/mol).

```python
from gromos import Energy
e = Energy()
```

### Properties

`total`, `kinetic`, `potential` — `float`

---

## Frame

A trajectory frame (time + step metadata + optional positions).

```python
from gromos import Frame
f = Frame(time=1.5, step=100)
print(f.time, f.step, f.n_atoms)
```

---

## Free functions

### `rmsd`

Root-mean-square deviation between two coordinate arrays, nm.

```python
from gromos import rmsd
import numpy as np

r = rmsd(
    positions.astype(np.float32),   # (N, 3) float32
    reference.astype(np.float32),   # (N, 3) float32
)
```

### `rdf`

Radial distribution function between two atom index groups.

```python
from gromos import rdf

r_vals, g_vals = rdf(
    positions.astype(np.float32),  # (N, 3) float32
    group_a=[0, 3, 6],             # list[int] — indices in group A
    group_b=[1, 4, 7],             # list[int] — indices in group B
    n_bins=100,
    r_max=1.5,                     # nm
)
# r_vals: (n_bins,) float64 — bin centres, nm
# g_vals: (n_bins,) float64 — g(r) values
```

---

## Planned API (not yet implemented)

### P3.3 — Energy timeseries without files

```python
# Run n steps, collect energies every ene_freq steps into a NumPy array.
# No .tre file written.
energies = sim.run(steps=5000, ene_freq=100)   # returns EnergyTimeseries
df       = energies.to_dataframe()
energies.plot("kinetic", "potential")
energies.block_average("total", block_size=50)
```

### FUTURE — System builder algebra

Design rationale in `FUTURE.md`. Nothing below is implemented.

```python
from gromos import ForceField

ff     = ForceField.load("54A7")
system = molecule("ALA", ff) * 10 + solvent("SPC", n=2000)
system.neutralize(ion="CL")
system.write("prepared.topo", "prepared.cnf")
```

---

## Legacy: `gromos.md_runners`

The `md_runners` sub-module contains Python classes that shell out to the `md`
binary (`MDSimulation`, `GaMDSimulation`, `EDSSimulation`, `REMDSimulation`,
`TISimulation`). They write temporary `.imd` files and parse output files.

Prefer `Simulation` for new code. The `md_runners` API will be deprecated once
`sim.run()` (P3.3) covers the same use cases.
