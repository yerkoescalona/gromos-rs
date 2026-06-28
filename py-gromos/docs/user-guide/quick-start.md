# Quick Start

## Installation

```bash
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromos-rs
make build-python        # creates .venv, compiles Rust, installs gromos
source .venv/bin/activate
```

See [Installation](installation.md) for details and troubleshooting.

---

## Loading a system

A `System` is a topology + coordinate file pair.  
The atom count is validated at construction — mismatched files raise `ValueError` immediately.

```python
from gromos import System

system = System.from_files("water_216.topo", "equilibrated.cnf")
print(system)
# System(n_atoms=648, charge=0, box=(3.105, 3.105, 3.105))

print(system.n_atoms)           # 648
print(system.charge)            # 0  — integer, units of e
print(system.positions.shape)   # (648, 3) — float64, nm
print(system.velocities.shape)  # (648, 3) — float64, nm/ps
print(system.box)               # (3.105, 3.105, 3.105) nm
```

You can also build a `System` from pre-loaded objects:

```python
from gromos import Topology, Configuration, System

topo   = Topology("water_216.topo")
conf   = Configuration("equilibrated.cnf")
system = System(topo, conf)   # same validation
```

### Inspecting topology and configuration directly

```python
topo = Topology("water_216.topo")
print(topo.n_atoms)          # total atoms (solute + solvent)
print(topo.n_solute_atoms)   # solute-only
print(topo.n_solvent_atoms)  # solvent-only
print(topo.masses.shape)     # (N,) float64, g/mol
print(topo.charges.shape)    # (N,) float64, e

conf = Configuration("equilibrated.cnf")
print(conf.positions.shape)   # (N, 3) float64, nm
print(conf.velocities.shape)  # (N, 3) float64, nm/ps
print(conf.box_dimensions)    # (Lx, Ly, Lz) nm
```

---

## Creating simulation parameters

Use a factory — no `.imd` file authoring required.

```python
from gromos import InputParameters

# NVT — Berendsen thermostat, τ = 0.1 ps
params = InputParameters.nvt(dt=0.002, steps=5000, temperature=300.0)

# NVE — microcanonical, no thermostat
params = InputParameters.nve(dt=0.002, steps=5000)

# NPT — Berendsen thermostat + barostat, water compressibility by default
params = InputParameters.npt(dt=0.002, steps=5000, temperature=300.0, pressure=1.0)

# Energy minimisation — steepest descent
params = InputParameters.steepest_descent(steps=500)

# Or load an existing GROMOS input file
params = InputParameters.from_file("run.imd")
params = InputParameters("run.imd")   # identical

print(params.dt, params.nstlim, params.temperature, params.cutoff)
```

---

## Running a simulation

```python
from gromos import Simulation

# Recommended: two-argument form
sim = Simulation(system, params)

# Legacy three-argument forms also work
sim = Simulation(topo, conf, params)
sim = Simulation("water_216.topo", "equilibrated.cnf", "run.imd")
```

### Advancing and reading state

```python
sim.step(100)   # advance 100 MD steps

# Energies (kJ/mol)
print(sim.total_energy)
print(sim.kinetic_energy)
print(sim.potential_energy)

# Positions and velocities ((N, 3) float64 NumPy arrays)
print(sim.positions.shape)   # nm
print(sim.velocities.shape)  # nm/ps
print(sim.forces.shape)      # kJ/(mol·nm)

# Clock
print(sim.time)          # ps
print(sim.current_step)  # step number

# Thermostat temperature
print(sim.temperature)   # K
```

### Collecting energies over a run

```python
import numpy as np

energies = []
for _ in range(50):
    sim.step(100)
    energies.append(sim.total_energy)

energies = np.array(energies)
print(f"Mean: {energies.mean():.1f} ± {energies.std():.1f} kJ/mol")
```

---

## Inspecting the algorithm sequence

Every `Simulation` runs a fixed sequence of algorithms each step. You can read
it out and even build a custom sequence before constructing the simulation.

```python
# Inspect what's running
print(sim.algorithm_names)
# ['RemoveCOMMotion', 'Forcefield', 'LeapFrogVelocity', 'BerendsenThermostat',
#  'LeapFrogPosition', 'TemperatureCalculation', 'EnergyCalculation']

# Build and modify a sequence before constructing the simulation
from gromos import AlgorithmSequence

seq = AlgorithmSequence.nvt(topo, params)
seq.remove("RemoveCOMMotion")   # disable COM motion removal
print(seq.names)
print("Forcefield" in seq)      # True

sim = Simulation.from_sequence(topo, conf, params, seq)
sim.step(100)
```

---

## Writing output

```python
# Write current coordinates to a GROMOS .cnf file
system.write("output.cnf")
```

---

## Units

| Quantity | Unit |
|----------|------|
| Length | nm |
| Time | ps |
| Energy | kJ/mol |
| Force | kJ/(mol·nm) |
| Velocity | nm/ps |
| Temperature | K |
| Pressure | bar |
| Mass | g/mol |
| Charge | elementary charge (e) |

---

## Next steps

- **[API Reference](../api/reference.md)** — every class and method documented
- **[Installation](installation.md)** — build options, troubleshooting
