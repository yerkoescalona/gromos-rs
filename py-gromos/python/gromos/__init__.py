"""
gromos — Python bindings for the GROMOS-RS molecular dynamics engine.

Quick start
-----------
>>> from gromos import System, Recipe, Simulation
>>>
>>> system = System.from_files("water_216.topo", "equilibrated.cnf")
>>> recipe = Recipe.nvt(dt=0.002, steps=1000, temperature=300.0)   # or Recipe.from_imd("run.imd")
>>> sim    = Simulation(system, recipe)
>>> sim.step(1000)
>>> print(sim.total_energy)   # kJ/mol

The recipe is the one description of a run — the same data a GROMOS ``.imd`` file carries,
grouped by concern (``recipe.control``, ``.ensemble``, ``.constraints`` …), plus additive
``terms`` (``Term("schnet", ...)``) and auxiliary ``inputs``. ``recipe.plan(system)`` is the MD
step as editable data; ``recipe.to_imd(n_atoms=...)`` is what gromosXX would run.
``gromos.terms()`` / ``gromos.algorithms()`` list what this build knows.

All units follow the GROMOS convention: nm, ps, kJ/mol, K.
"""

# ── Rust extension (always available) ────────────────────────────────────────
from .gromos import (
    Algorithm,
    # Algorithm-sequence building blocks (deprecated: Recipe.plan(system) + Plan)
    AlgorithmSequence,
    BerendsenBarostat,
    BerendsenThermostat,
    Configuration,
    Energy,
    EnergyCalculation,
    Forcefield,
    Frame,
    InputParameters,  # deprecated: Recipe.from_imd / Recipe.nve|nvt|npt|minimize
    LeapFrogIntegrator,
    LeapFrogPosition,
    LeapFrogVelocity,
    Plan,
    PressureCalculation,
    # The run as data (PLAN.md 3.9)
    Recipe,
    RemoveCOMMotion,
    ShakeConstraints,
    Simulation,
    SteepestDescent,
    # Simulation objects
    System,
    TemperatureCalculation,
    Term,
    # I/O objects
    Topology,
    # Math / utility
    Vec3,
    # QM potential (real xtb — always available, no --features ml needed)
    XtbPotential,
    algorithms,
    build_info,
    rdf,
    rmsd,
    terms,
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

# ── Exceptions (subclasses of the builtins the binding raised before) ─────────
# ── Legacy subprocess runners (depend on the `md` binary being in PATH) ──────
# These wrap the GROMOS-RS command-line tool and write temporary files.
# Prefer the Simulation class for new code.
# ── Analysis (subprocess wrappers, mostly stubs) ─────────────────────────────
from . import analysis, exceptions, md_runners
from .exceptions import MissingFeatureError, PlanError, RecipeError, RunError

# ── Energy timeseries wrapper for Simulation.run() ────────────────────────────
from .timeseries import EnergyTimeseries

try:
    from importlib.metadata import version as _dist_version

    __version__ = _dist_version("gromos")
except Exception:  # not installed as a distribution (e.g. imported from a checkout)
    __version__ = "0.0.0+unknown"

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
    # The run as data
    "Recipe",
    "Term",
    "Algorithm",
    "Plan",
    "terms",
    "algorithms",
    "build_info",
    # Exceptions
    "exceptions",
    "RecipeError",
    "PlanError",
    "MissingFeatureError",
    "RunError",
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
