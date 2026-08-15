"""
gromos — Python bindings for the GROMOS-RS molecular dynamics engine.

Quick start
-----------
>>> from gromos import System, InputParameters, Simulation
>>>
>>> system = System.from_files("water_216.topo", "equilibrated.cnf")
>>> params = InputParameters.nvt(dt=0.002, steps=1000, temperature=300.0)
>>> sim    = Simulation(system, params)
>>> sim.step(1000)
>>> print(sim.total_energy)   # kJ/mol

All units follow the GROMOS convention: nm, ps, kJ/mol, K.
"""

# ── Rust extension (always available) ────────────────────────────────────────
from .gromos import (
    # Math / utility
    Vec3,
    Energy,
    Frame,
    rmsd,
    rdf,
    # I/O objects
    Topology,
    Configuration,
    InputParameters,
    # Simulation objects
    System,
    Simulation,
    # Algorithm-sequence building blocks
    AlgorithmSequence,
    Forcefield,
    LeapFrogIntegrator,
    LeapFrogVelocity,
    LeapFrogPosition,
    BerendsenThermostat,
    BerendsenBarostat,
    ShakeConstraints,
    SteepestDescent,
    TemperatureCalculation,
    PressureCalculation,
    EnergyCalculation,
    RemoveCOMMotion,
    # QM potential (real xtb — always available, no --features ml needed)
    XtbPotential,
)

# ── ML potential API (only present if the extension was built --features ml) ─
# `Simulation(..., ml_potential=..., ml_region=..., ml_buffer=...)` also only accepts
# real values for these kwargs in an `ml`-enabled build; the kwargs themselves always
# exist, but passing a SchNetPotential without one raises a clear RuntimeError instead
# of ImportError-ing at import time — see `simulation.rs::resolve_ml_spec`.
try:
    from .gromos import SchNetPotential, resolve_zone_partition

    _HAS_ML = True
except ImportError:
    _HAS_ML = False

# ── Energy timeseries wrapper for Simulation.run() ────────────────────────────
from .timeseries import EnergyTimeseries

# ── Legacy subprocess runners (depend on the `md` binary being in PATH) ──────
# These wrap the GROMOS-RS command-line tool and write temporary files.
# Prefer the Simulation class for new code.
from . import md_runners

# ── Analysis (subprocess wrappers, mostly stubs) ─────────────────────────────
from . import analysis

__version__ = "0.1.0"

__all__ = [
    # Math / utility
    "Vec3",
    "Energy",
    "Frame",
    "rmsd",
    "rdf",
    # I/O
    "Topology",
    "Configuration",
    "InputParameters",
    # Simulation
    "System",
    "Simulation",
    # Algorithm sequence
    "AlgorithmSequence",
    "Forcefield",
    "LeapFrogIntegrator",
    "LeapFrogVelocity",
    "LeapFrogPosition",
    "BerendsenThermostat",
    "BerendsenBarostat",
    "ShakeConstraints",
    "SteepestDescent",
    "TemperatureCalculation",
    "PressureCalculation",
    "EnergyCalculation",
    "RemoveCOMMotion",
    # QM potential
    "XtbPotential",
    # Timeseries
    "EnergyTimeseries",
    # Sub-modules
    "md_runners",
    "analysis",
    # Meta
    "__version__",
]

if _HAS_ML:
    __all__ += ["SchNetPotential", "resolve_zone_partition"]
