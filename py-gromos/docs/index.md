# py-gromos

Python bindings for the GROMOS-RS molecular dynamics engine.

## Five lines to a running simulation

```python
from gromos import System, InputParameters, Simulation

system = System.from_files("water_216.topo", "equilibrated.cnf")
params = InputParameters.nvt(dt=0.002, steps=5000, temperature=300.0)
sim    = Simulation(system, params)
sim.step(1000)

print(sim.total_energy)     # kJ/mol
print(sim.positions.shape)  # (648, 3)
```

No `.imd` files to author. No temporary directories. All in memory.

## What works today

| Class / function | Description |
|-----------------|-------------|
| `System` | Paired topology + coordinates, atom-count validated |
| `Topology` | Load `.topo` — atoms, masses, charges, bonds |
| `Configuration` | Load `.cnf` / `.g96` — positions, velocities, box |
| `InputParameters` | Load `.imd` or build via `nve / nvt / npt / steepest_descent` factories |
| `Simulation` | Algorithm sequence + `step(n)` + full energy/position/force access |
| `AlgorithmSequence` | Inspect and modify the MD pipeline at the step level |
| `Vec3`, `Energy`, `Frame` | Utility types |
| `rmsd`, `rdf` | Analysis functions |

The engine is validated bit-for-bit against double-precision gromosXX output
(energy rel. tolerance 1 × 10⁻⁸, force abs. tolerance 1 × 10⁻⁶ kJ mol⁻¹ nm⁻¹)
across 21 reference systems.

## What is coming next

| Milestone | Description |
|-----------|-------------|
| **P3.3** | `sim.run(steps, ene_freq)` returns a NumPy energy array — no `.tre` file |
| **P3.4** | Working Jupyter notebooks for the `from_files → nvt → run` path |
| **FUTURE** | System builder algebra: `molecule("ALA", ff) * 10 + solvent("SPC", n=2000)` |

## Architecture

```
Python user code
      │
      ▼
gromos/__init__.py   — clean re-export of working names only
      │
      ▼
gromos/gromos.abi3.so   — PyO3 extension compiled from pyo3-gromos (Rust)
      │
      ▼
gromos-rs workspace
  gromos-core · gromos-io · gromos-forces · gromos-integrators
```

Positions, velocities, and forces are returned as `(N, 3)` float64 NumPy arrays.

## Documentation

- **[Quick Start](user-guide/quick-start.md)** — first simulation step by step
- **[API Reference](api/reference.md)** — every class and method
- **[Installation](user-guide/installation.md)** — build from source

## License

GPL-2.0 — same as GROMOS.
