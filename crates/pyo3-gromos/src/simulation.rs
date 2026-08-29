//! Interactive Simulation API for Python
//!
//! Supports both compositional and file-based construction, mirroring
//! the internal md architecture with GROMOS naming conventions.
//!
//! Since PLAN.md 3.9 step 1 this file contains **no run assembly of its own**: a
//! `Simulation` is built by the same `gromos-run` functions the `md` binary calls
//! (`prepare_system` → `build_sequence_from_imd` → `start`), so the gromosXX reference suite
//! (which drives the binary) and the Python suite guard one code path from two sides.
//!
//! # Example (Python)
//!
//! ```python
//! from gromos import Simulation, Topology, Configuration, InputParameters
//!
//! # Compositional — mirrors md internals
//! topo = Topology("system.topo")
//! conf = Configuration("initial.cnf")
//! params = InputParameters("run.imd")
//! sim = Simulation(topo, conf, params)
//!
//! # File-based (convenience, backward-compatible)
//! sim = Simulation("system.topo", "initial.cnf", "run.imd")
//!
//! # Both work identically
//! sim.step(100)
//! print(sim.energies, sim.positions)
//! print(sim.algorithm_names)  # inspect the MD sequence
//! ```

use std::path::PathBuf;

use numpy::PyArray2;
use pyo3::prelude::*;

use gromos_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::{BoxType, Configuration},
    math::{truncoct_triclinic_rotmat, Vec3},
    units::PhysicalConstants,
    Topology,
};
#[cfg(feature = "ml")]
use gromos_forces::zones::ZonePartition;
use gromos_io::{
    coordinate::read_coordinates,
    imd::{read_imd_file, ImdParameters},
    topology::{build_topology, read_topology_file},
};
use gromos_run::{
    build_sequence_from_imd, prepare_system, start, Coordinates, RunError, RunInputs, RunOptions,
};

use super::algorithm_sequence::{
    compute_total_dof, resolve_algorithm_sequence, PyAlgorithmSequence,
};
use super::parameters::PyInputParameters;
use super::py_conf::PyConfiguration;
use super::system::PySystem;
use super::topology::PyTopology;
use super::PyEnergy;

/// Map a `gromos_run::RunError` onto the builtin Python exception the binding raised for
/// the same condition before the assembly was shared (no behaviour change for callers).
pub(crate) fn run_err(e: RunError) -> PyErr {
    match e {
        RunError::Io { .. } => PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()),
        RunError::MissingInput { .. }
        | RunError::AtomCountMismatch { .. }
        | RunError::SolventCount { .. }
        | RunError::Inconsistent(_)
        | RunError::Validation { .. } => {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
        },
        RunError::Init(_) | RunError::Step { .. } => {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string())
        },
    }
}

/// Plain data extracted from a Python `SchNetPotential` + region/buffer selector strings,
/// threaded from `Simulation::new` down to `build_simulation`. Deliberately *not* itself a
/// pyo3 type (no lifetime, always compiles) so `build_simulation`'s signature doesn't need to
/// change based on the `ml` feature; only the code that actually *acts* on a `Some` value is
/// feature-gated (see `build_simulation`).
#[derive(Clone)]
#[cfg_attr(not(feature = "ml"), allow(dead_code))]
pub(crate) struct MlPotentialSpec {
    pub model_path: String,
    pub cutoff: f64,
    pub elements: Vec<i64>,
    pub region: String,
    pub buffer: Option<String>,
}

/// Extract a `MlPotentialSpec` from `Simulation.__init__`'s `ml_potential=`/`ml_region=`/
/// `ml_buffer=` kwargs. `ml_potential` is a `SchNetPotential` (a `#[cfg(feature = "ml")]`
/// pyclass) — takes `&Bound<'_, PyAny>` rather than that concrete type so this function's own
/// signature doesn't need to change across builds; the actual downcast is feature-gated inside.
fn resolve_ml_spec(
    ml_potential: Option<&Bound<'_, PyAny>>,
    ml_region: Option<String>,
    ml_buffer: Option<String>,
) -> PyResult<Option<MlPotentialSpec>> {
    match (ml_potential, ml_region) {
        (None, None) => Ok(None),
        (Some(_), None) | (None, Some(_)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "ml_potential and ml_region must be given together",
        )),
        (Some(potential), Some(region)) => {
            #[cfg(feature = "ml")]
            {
                let potential = potential
                    .extract::<PyRef<crate::ml_potential::PySchNetPotential>>()
                    .map_err(|_| {
                        PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                            "ml_potential must be a SchNetPotential",
                        )
                    })?;
                Ok(Some(MlPotentialSpec {
                    model_path: potential.model_path.clone(),
                    cutoff: potential.cutoff,
                    elements: potential.elements.clone(),
                    region,
                    buffer: ml_buffer,
                }))
            }
            #[cfg(not(feature = "ml"))]
            {
                let _ = (potential, region, ml_buffer);
                Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                    "ML support was not compiled into this build (rebuild pyo3-gromos/py-gromos \
                     with --features ml)",
                ))
            }
        },
    }
}

/// Build a simulation from raw components — the shared core of every constructor.
///
/// Prepares the system and builds the algorithm sequence through `gromos-run` (exactly what
/// the `md` binary does), optionally inserts the ML orchestrator term right after
/// `Forcefield`, then runs `init` + step 0.
fn build_simulation(
    topo: Topology,
    physical_constants: PhysicalConstants,
    positions: Vec<Vec3>,
    velocities: Vec<Vec3>,
    box_dims: Vec3,
    imd: &ImdParameters,
    inputs: &RunInputs,
    ml_spec: Option<&MlPotentialSpec>,
) -> PyResult<PySimulation> {
    let coords = Coordinates {
        positions,
        velocities,
        box_dims,
    };
    let prepared =
        prepare_system(imd, topo, physical_constants, coords, inputs).map_err(run_err)?;
    let built =
        build_sequence_from_imd(imd, &prepared, inputs, &RunOptions::default()).map_err(run_err)?;
    #[cfg(feature = "ml")]
    let periodicity_for_ml = gromos_run::periodicity_of(&prepared);

    let gromos_run::Built {
        sequence: mut md_sequence,
        summary,
    } = built;
    let gromos_run::Prepared {
        topology: topo,
        configuration: mut conf,
        ..
    } = prepared;
    let n_atoms = topo.num_atoms();

    // ML potential, if attached (PLAN.md P3.7) — inserted immediately after Forcefield,
    // matching `orchestrator_algorithm.rs`'s documented placement requirement (it adds to
    // Forcefield's already-computed force/virial). `ml_spec` is always
    // `Option<MlPotentialSpec>` regardless of the `ml` feature; only the code that *acts* on
    // `Some` is feature-gated, so `Simulation`'s constructor signature doesn't change.
    #[cfg(feature = "ml")]
    if let Some(spec) = ml_spec {
        let partition = ZonePartition::from_selections(&topo, &spec.region, spec.buffer.as_deref())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        let potential = crate::ml_potential::PySchNetPotential::from_spec(spec);
        let algorithm = crate::ml_potential::build_ml_orchestrator_algorithm(
            &potential,
            &partition,
            n_atoms,
            periodicity_for_ml,
        )?;
        let after_forcefield = md_sequence
            .algorithm_names()
            .iter()
            .position(|name| *name == "Forcefield")
            .map(|i| i + 1)
            .ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                    "built sequence has no Forcefield to attach the ML term to",
                )
            })?;
        md_sequence.insert(after_forcefield, Box::new(algorithm));
    }
    #[cfg(not(feature = "ml"))]
    {
        let _ = &ml_spec;
    }

    let mut sim_state = SimulationState::new(imd.dt, imd.nstlim);
    start(&mut md_sequence, &topo, &mut conf, &sim_state).map_err(run_err)?;
    sim_state.advance();

    Ok(PySimulation {
        topology: topo,
        configuration: conf,
        md_sequence,
        sim_state,
        dt: imd.dt,
        n_atoms,
        total_dof: summary.total_dof,
    })
}

/// Build a simulation from a user-provided algorithm sequence.
///
/// The system is prepared exactly as for `build_simulation`; the sequence descriptors are
/// then resolved into real Rust algorithms (`resolve_algorithm_sequence`, the descriptor path
/// slated for replacement by the recipe plan in PLAN.md 3.9 steps 2–4).
fn build_simulation_from_sequence(
    topo: Topology,
    physical_constants: PhysicalConstants,
    positions: Vec<Vec3>,
    velocities: Vec<Vec3>,
    box_dims: Vec3,
    imd: &ImdParameters,
    sequence: &PyAlgorithmSequence,
) -> PyResult<PySimulation> {
    let coords = Coordinates {
        positions,
        velocities,
        box_dims,
    };
    let prepared = prepare_system(imd, topo, physical_constants, coords, &RunInputs::default())
        .map_err(run_err)?;
    let gromos_run::Prepared {
        topology: topo,
        configuration: mut conf,
        box_dims,
        ..
    } = prepared;
    let n_atoms = topo.num_atoms();
    let total_dof = compute_total_dof(&topo, imd);

    let mut md_sequence =
        resolve_algorithm_sequence(sequence, &topo, &conf, imd, physical_constants, box_dims)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to resolve algorithm sequence: {}",
                    e
                ))
            })?;

    let mut sim_state = SimulationState::new(imd.dt, imd.nstlim);
    start(&mut md_sequence, &topo, &mut conf, &sim_state).map_err(run_err)?;
    sim_state.advance();

    Ok(PySimulation {
        topology: topo,
        configuration: conf,
        md_sequence,
        sim_state,
        dt: imd.dt,
        n_atoms,
        total_dof,
    })
}

// ============================================================================
// PySimulation
// ============================================================================

/// A GROMOS-RS molecular dynamics simulation.
///
/// Interactive simulation object inspired by OpenMM's Simulation class,
/// using GROMOS naming conventions and file formats.
///
/// Create from Topology + Configuration + InputParameters objects,
/// or directly from file paths. Call `step(n)` to advance and access
/// properties like `positions`, `forces`, `energies`, `temperature`.
#[pyclass(name = "Simulation", unsendable)]
pub struct PySimulation {
    topology: Topology,
    configuration: Configuration,
    md_sequence: AlgorithmSequence,
    sim_state: SimulationState,
    dt: f64,
    n_atoms: usize,
    /// Kinetic degrees of freedom (constraint- and NDFMIN-aware), computed once
    /// at build time by `gromos_run::total_dof` and reused for `temperature` —
    /// the same value the thermostat (if any) couples to.
    total_dof: f64,
}

#[pymethods]
impl PySimulation {
    /// Create a new simulation.
    ///
    /// Three-argument forms (backward-compatible):
    ///   - `Simulation(topology, configuration, parameters)` — pre-loaded objects
    ///   - `Simulation("system.topo", "initial.cnf", "run.imd")` — file paths
    ///
    /// Two-argument form (new):
    ///   - `Simulation(system, parameters)` — System object + InputParameters
    ///
    /// Optional restraint file paths (mirror the `md` binary's `@distrest`/
    /// `@posresspec`/`@refpos` flags): `distrest`, `posresspec`, `refpos`.
    /// A perturbation topology (FEP, `NTG != 0`) is applied to the `Topology`
    /// beforehand with `Topology.apply_perturbation(path)`.
    #[new]
    #[pyo3(signature = (arg1, arg2, arg3=None, *, distrest=None, posresspec=None, refpos=None, ml_potential=None, ml_region=None, ml_buffer=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        arg1: &Bound<'_, PyAny>,
        arg2: &Bound<'_, PyAny>,
        arg3: Option<&Bound<'_, PyAny>>,
        distrest: Option<String>,
        posresspec: Option<String>,
        refpos: Option<String>,
        ml_potential: Option<&Bound<'_, PyAny>>,
        ml_region: Option<String>,
        ml_buffer: Option<String>,
    ) -> PyResult<Self> {
        let inputs = RunInputs {
            pttopo: None,
            posresspec: posresspec.map(PathBuf::from),
            refpos: refpos.map(PathBuf::from),
            distrest: distrest.map(PathBuf::from),
        };
        let ml_spec = resolve_ml_spec(ml_potential, ml_region, ml_buffer)?;
        match arg3 {
            None => {
                // Two-arg form: Simulation(system, params)
                let system = arg1.extract::<PyRef<PySystem>>().map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "With two arguments, first must be a System object",
                    )
                })?;
                let params = arg2.extract::<PyRef<PyInputParameters>>().map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "With two arguments, second must be an InputParameters object",
                    )
                })?;
                build_simulation(
                    system.topology.inner.clone(),
                    system.topology.physical_constants,
                    system.configuration.pos_data.clone(),
                    system.configuration.vel_data.clone(),
                    system.configuration.box_dims,
                    &params.inner,
                    &inputs,
                    ml_spec.as_ref(),
                )
            },
            Some(arg3) => {
                // Three-arg forms: file paths or (Topology, Configuration, InputParameters)
                if let (Ok(topo_file), Ok(conf_file), Ok(input_file)) = (
                    arg1.extract::<String>(),
                    arg2.extract::<String>(),
                    arg3.extract::<String>(),
                ) {
                    return Self::_from_files_with_inputs(
                        &topo_file,
                        &conf_file,
                        &input_file,
                        &inputs,
                        ml_spec.as_ref(),
                    );
                }

                let topo = arg1.extract::<PyRef<PyTopology>>().map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "First argument must be a file path (str) or Topology object",
                    )
                })?;
                let conf = arg2.extract::<PyRef<PyConfiguration>>().map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "Second argument must be a file path (str) or Configuration object",
                    )
                })?;
                let params = arg3.extract::<PyRef<PyInputParameters>>().map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "Third argument must be a file path (str) or InputParameters object",
                    )
                })?;

                build_simulation(
                    topo.inner.clone(),
                    topo.physical_constants,
                    conf.pos_data.clone(),
                    conf.vel_data.clone(),
                    conf.box_dims,
                    &params.inner,
                    &inputs,
                    ml_spec.as_ref(),
                )
            },
        }
    }

    /// Create a simulation from file paths (alternative to constructor).
    ///
    /// Example:
    ///     sim = Simulation.from_files("system.topo", "initial.cnf", "run.imd")
    #[staticmethod]
    fn from_files(topo_file: &str, conf_file: &str, input_file: &str) -> PyResult<Self> {
        Self::_from_files_with_inputs(
            topo_file,
            conf_file,
            input_file,
            &RunInputs::default(),
            None,
        )
    }

    /// Run the simulation for `n_steps` MD steps.
    ///
    /// Example:
    ///     sim.step(1000)  # advance 1000 steps
    fn step(&mut self, n_steps: usize) -> PyResult<()> {
        for _ in 0..n_steps {
            self.md_sequence
                .run_step(&self.topology, &mut self.configuration, &self.sim_state)
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Error at step {}: {}",
                        self.sim_state.step, e
                    ))
                })?;
            self.sim_state.advance();
        }
        Ok(())
    }

    /// Run `n_steps` MD steps, sampling energies every `ene_freq` steps.
    ///
    /// Returns an `(n_frames, 12)` numpy array with columns
    /// `[time, kinetic, potential, total, volume, pressure, bond, angle, improper, dihedral, lj, coulomb]`
    /// — the same component order as the `.tre` file written by the `md` binary
    /// (see `gromos_io::energy::EnergyFrame`), so a `run()` array and a `.tre`
    /// dump of the same trajectory line up column-for-column.
    /// Frame 0 is the state before any of these steps ran; subsequent frames
    /// are sampled after every `ene_freq`-th step.
    ///
    /// Example:
    ///     energies = sim.run(1000, ene_freq=100)  # (11, 12) array
    #[pyo3(signature = (n_steps, ene_freq=100))]
    fn run<'py>(
        &mut self,
        py: Python<'py>,
        n_steps: usize,
        ene_freq: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let mut rows: Vec<Vec<f64>> = vec![self.energy_row()];
        for i in 0..n_steps {
            self.md_sequence
                .run_step(&self.topology, &mut self.configuration, &self.sim_state)
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Error at step {}: {}",
                        self.sim_state.step, e
                    ))
                })?;
            self.sim_state.advance();
            if (i + 1) % ene_freq == 0 {
                rows.push(self.energy_row());
            }
        }
        Ok(PyArray2::from_vec2_bound(py, &rows)?)
    }

    // -- State getters -------------------------------------------------------

    /// Current simulation time in picoseconds.
    #[getter]
    fn time(&self) -> f64 {
        self.sim_state.time
    }

    /// Current step number.
    #[getter]
    fn current_step(&self) -> usize {
        self.sim_state.step
    }

    /// Time step size in picoseconds.
    #[getter]
    fn dt(&self) -> f64 {
        self.dt
    }

    /// Number of atoms in the system.
    #[getter]
    fn n_atoms(&self) -> usize {
        self.n_atoms
    }

    /// Names of algorithms in the MD sequence.
    ///
    /// Example:
    ///     print(sim.algorithm_names)
    ///     # ['Forcefield', 'LeapFrogVelocity', 'LeapFrogPosition', ...]
    #[getter]
    fn algorithm_names(&self) -> Vec<String> {
        self.md_sequence
            .algorithm_names()
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    // -- Phase-space getters -------------------------------------------------

    /// Current positions as an Nx3 numpy array (nm).
    ///
    /// Returns positions from the "old" state, which holds the most
    /// recently completed step's data (GROMOS convention).
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = state.pos.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Current velocities as an Nx3 numpy array (nm/ps).
    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = state.vel.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Current forces as an Nx3 numpy array (kJ/mol/nm).
    ///
    /// For a truncated-octahedron box, GROMOS rotates FREEFORCERED output
    /// back into the original (cube) frame via
    /// `truncoct_triclinic_rotmat(false)` — positions/velocities stay in the
    /// triclinic frame, only forces are rotated back — mirrored here so this
    /// getter agrees with the `md` binary's `.trf` output (no ROTTRANS block
    /// support here, so `rmat(phi,theta,psi)` is identity, same as `md.rs`).
    #[getter]
    fn forces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = if state.box_config.box_type == BoxType::TruncatedOctahedral {
            let rot = truncoct_triclinic_rotmat(false);
            state
                .force
                .iter()
                .map(|f| {
                    let r = rot * *f;
                    vec![r.x, r.y, r.z]
                })
                .collect()
        } else {
            state.force.iter().map(|v| vec![v.x, v.y, v.z]).collect()
        };
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    // -- Energy getters ------------------------------------------------------

    /// Current energies (total, kinetic, potential, and components).
    #[getter]
    fn energies(&self) -> PyEnergy {
        let e = &self.configuration.old().energies;
        PyEnergy {
            total: e.total(),
            kinetic: e.kinetic_total,
            potential: e.potential_total,
            bond: e.bond_total,
            angle: e.angle_total,
            dihedral: e.dihedral_total,
            improper: e.improper_total,
            lj: e.lj_total,
            coulomb: e.crf_total,
        }
    }

    /// Total energy (kJ/mol).
    #[getter]
    fn total_energy(&self) -> f64 {
        self.configuration.old().energies.total()
    }

    /// Potential energy (kJ/mol).
    #[getter]
    fn potential_energy(&self) -> f64 {
        self.configuration.old().energies.potential_total
    }

    /// Kinetic energy (kJ/mol).
    #[getter]
    fn kinetic_energy(&self) -> f64 {
        self.configuration.old().energies.kinetic_total
    }

    /// Current temperature (K), computed from kinetic energy.
    ///
    /// Uses the same constraint- and NDFMIN-aware degrees-of-freedom count the
    /// thermostat couples to (`total_dof`), not a bare `3*n_atoms` — a
    /// constrained system's true DOF is lower, and using the unconstrained
    /// count here used to silently disagree with what the thermostat targets.
    #[getter]
    fn temperature(&self) -> f64 {
        let n_dof = self.total_dof.round() as usize;
        self.configuration.old().temperature(n_dof)
    }

    /// Current box volume (nm³).
    #[getter]
    fn volume(&self) -> f64 {
        self.configuration.old().box_config.volume()
    }

    /// Current instantaneous pressure (bar): `P = (2*KE - virial) / (3*V)`.
    ///
    /// The virial term is only populated by `PressureCalculation` in the
    /// algorithm sequence — automatic for NPT (`InputParameters.npt()`/
    /// `AlgorithmSequence.npt()`), absent for NVE/NVT. Without it this returns
    /// the *kinetic-only* term (`2*KE / (3*V)`, not zero) — a physically
    /// incomplete, misleadingly large number, not "no pressure to report".
    /// Only trust this getter under NPT.
    #[getter]
    fn pressure(&self) -> f64 {
        self.configuration.old().pressure()
    }

    // -- Topology getters ----------------------------------------------------

    /// Number of solute atoms.
    #[getter]
    fn n_solute_atoms(&self) -> usize {
        self.topology.num_solute_atoms()
    }

    /// Number of solvent atoms.
    #[getter]
    fn n_solvent_atoms(&self) -> usize {
        self.n_atoms - self.topology.num_solute_atoms()
    }

    /// Create a simulation with a custom algorithm sequence.
    ///
    /// This is the most flexible constructor — full control over the MD step.
    /// The sequence determines exactly what happens each step and in what order.
    ///
    /// Args:
    ///     topo: Topology object
    ///     conf: Configuration object
    ///     params: InputParameters (still needed for dt, nstlim, and defaults)
    ///     sequence: AlgorithmSequence defining the MD step
    ///
    /// Example:
    ///     seq = AlgorithmSequence.nvt(topo, params)
    ///     seq.remove("RemoveCOMMotion")  # customize
    ///     sim = Simulation.from_sequence(topo, conf, params, seq)
    ///     sim.step(1000)
    #[staticmethod]
    fn from_sequence(
        topo: &PyTopology,
        conf: &PyConfiguration,
        params: &PyInputParameters,
        sequence: &PyAlgorithmSequence,
    ) -> PyResult<Self> {
        build_simulation_from_sequence(
            topo.inner.clone(),
            topo.physical_constants,
            conf.pos_data.clone(),
            conf.vel_data.clone(),
            conf.box_dims,
            &params.inner,
            sequence,
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "Simulation(n_atoms={}, step={}, time={:.3} ps, E_tot={:.6e} kJ/mol)",
            self.n_atoms,
            self.sim_state.step,
            self.sim_state.time,
            self.configuration.old().energies.total(),
        )
    }
}

// Private helpers
impl PySimulation {
    /// `[time, kinetic, potential, total, volume, pressure, bond, angle, improper, dihedral, lj, coulomb]`
    /// for the current state — same layout as `EnergyFrame`, used by `run()` to
    /// build the energy timeseries. Temperature is left at 0.0: unlike volume
    /// and pressure it needs the degrees-of-freedom count, which this row does
    /// not carry (use the `temperature` getter).
    fn energy_row(&self) -> Vec<f64> {
        let state = self.configuration.old();
        let volume = state.box_config.volume();
        let pressure = state.pressure();
        let frame = gromos_io::energy::EnergyFrame::from_energy(
            &state.energies,
            self.sim_state.time,
            0.0,
            volume,
            pressure,
        );
        vec![
            frame.time,
            frame.kinetic,
            frame.potential,
            frame.total,
            frame.volume,
            frame.pressure,
            frame.bond,
            frame.angle,
            frame.improper,
            frame.dihedral,
            frame.lj,
            frame.coul_real,
        ]
    }

    fn _from_files_with_inputs(
        topo_file: &str,
        conf_file: &str,
        input_file: &str,
        inputs: &RunInputs,
        ml_spec: Option<&MlPotentialSpec>,
    ) -> PyResult<Self> {
        let imd = read_imd_file(input_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read input file '{}': {}",
                input_file, e
            ))
        })?;

        let parsed = read_topology_file(topo_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read topology '{}': {}",
                topo_file, e
            ))
        })?;
        let physical_constants = parsed.physical_constants;
        let topo = build_topology(parsed);

        let coord_data = read_coordinates(conf_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read coordinates '{}': {}",
                conf_file, e
            ))
        })?;

        build_simulation(
            topo,
            physical_constants,
            coord_data.positions,
            coord_data.velocities,
            coord_data.box_dims,
            &imd,
            inputs,
            ml_spec,
        )
    }
}

/// Register simulation bindings in the Python module.
pub fn register_simulation(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySimulation>()?;
    Ok(())
}
