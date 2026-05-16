//! Interactive Simulation API for Python
//!
//! Supports both compositional and file-based construction, mirroring
//! the internal md architecture with gromosXX naming conventions.
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

use pyo3::prelude::*;
use numpy::PyArray2;

use gromos_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::{Box as SimBox, Configuration},
    math::{Periodicity, Rectangular, Vacuum, Vec3},
    pairlist::{PairlistContainer, StandardPairlistAlgorithm},
    Topology,
};
use gromos_forces::nonbonded::CRFParameters;
use gromos_integrators::algorithms::{
    BerendsenBarostat, BerendsenBarostatParams, BerendsenThermostat,
    EnergyCalculation, Forcefield, LeapFrogPosition, LeapFrogVelocity,
    PressureCalculation, RemoveCOMMotion, ShakeAlgorithm, TemperatureCalculation,
    VirialType,
};
use gromos_integrators::constraints::{NtcMode, ShakeParameters};
use gromos_io::{
    coordinate::read_coordinates,
    imd::{read_imd_file, ImdParameters},
    topology::{build_topology, read_topology_file},
};

use super::PyEnergy;
use super::topology::PyTopology;
use super::py_conf::PyConfiguration;
use super::parameters::PyInputParameters;

// ============================================================================
// Shared build logic — constructs a fully initialized simulation from parts
// ============================================================================

/// Build a simulation from raw components.
///
/// This is the shared core used by both the file-path and object constructors.
/// It solvates the topology (if not already solvated), validates atom counts,
/// builds the Configuration + AlgorithmSequence, and runs step 0.
fn build_simulation(
    mut topo: Topology,
    positions: Vec<Vec3>,
    velocities: Vec<Vec3>,
    box_dims: Vec3,
    imd: &ImdParameters,
) -> PyResult<PySimulation> {
    // Solvate topology if not already solvated
    if topo.num_atoms() == topo.num_solute_atoms() && imd.nsm > 0 {
        topo.solvate(imd.nsm);
    }

    let n_atoms = topo.num_atoms();
    if positions.len() != n_atoms {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Atom count mismatch: topology={}, coordinates={}",
            n_atoms,
            positions.len()
        )));
    }

    // Build Configuration (double-buffered state)
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = positions;
    conf.current_mut().vel = if velocities.len() == n_atoms {
        velocities
    } else {
        vec![Vec3::ZERO; n_atoms]
    };
    conf.current_mut().box_config = SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z);
    conf.copy_current_to_old();

    // Extract parameters
    let dt = imd.dt;
    let n_steps = imd.nstlim;
    let cutoff = imd.rcutl;
    let rf_epsilon = imd.epsrf;
    let rf_kappa = imd.appak;
    let pairlist_update = imd.nsnb;
    let ntf = imd.ntf;

    let shake_enabled = imd.ntc > 1 || (imd.ntcs > 0 && imd.nsm > 0);
    let shake_tolerance = imd.shake_tol;

    let temperature = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].temp0.is_empty() {
        imd.temp_bath[0].temp0[0]
    } else {
        300.0
    };
    let thermostat_tau = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].tau.is_empty() {
        imd.temp_bath[0].tau[0]
    } else {
        -1.0
    };
    let thermostat_on = thermostat_tau > 0.0;

    // Nonbonded interactions
    let crf_params = CRFParameters::new(cutoff, 1.0, rf_epsilon, rf_kappa);
    let lj_params = Forcefield::convert_lj_parameters(&topo);

    // Pairlist
    let mut pairlist = PairlistContainer::new(imd.rcutp, cutoff, 0.0);
    pairlist.update_frequency = pairlist_update;

    let use_chargegroups = !topo.chargegroups.is_empty();
    let pairlist_algorithm = StandardPairlistAlgorithm::new(use_chargegroups);
    let periodicity = if box_dims.x == 0.0 && box_dims.y == 0.0 && box_dims.z == 0.0 {
        Periodicity::Vacuum(Vacuum)
    } else {
        Periodicity::Rectangular(Rectangular::new(box_dims))
    };

    pairlist_algorithm.update(&topo, &conf, &mut pairlist, &periodicity);

    // Build algorithm sequence (gromosXX Leap-Frog pattern)
    let mut md_sequence = AlgorithmSequence::new();

    // 1. COM motion removal
    if imd.nticom >= 1 || imd.nscm > 0 {
        md_sequence.push(Box::new(RemoveCOMMotion::new(imd.nticom, imd.nscm)));
    }

    // 2. Forcefield (bonded + nonbonded)
    let mut forcefield = Forcefield::new(
        lj_params, crf_params, periodicity, pairlist, pairlist_algorithm,
    );
    forcefield.ntf_bond = ntf[0] != 0;
    forcefield.ntf_angle = ntf[1] != 0;
    forcefield.ntf_improper = ntf[2] != 0;
    forcefield.ntf_dihedral = ntf[3] != 0;
    if !topo.solvent_atom_template.is_empty() {
        forcefield.atoms_per_solvent = topo.solvent_atom_template.len();
    }
    if imd.couple_pressure {
        forcefield.virial_type = match imd
            .pressure_parameters
            .as_ref()
            .map(|p| p.virial)
            .unwrap_or(0)
        {
            2 => VirialType::Molecular,
            1 => VirialType::Atomic,
            _ => VirialType::None,
        };
    }
    md_sequence.push(Box::new(forcefield));

    // 3. Leap-Frog velocity
    md_sequence.push(Box::new(LeapFrogVelocity::new()));

    // 3b. Berendsen thermostat
    if thermostat_on {
        let n_solute = topo.num_solute_atoms();
        let atoms_per_solvent = if !topo.solvent_atom_template.is_empty() {
            topo.solvent_atom_template.len()
        } else {
            1
        };
        let n_solvent_molecules = if atoms_per_solvent > 0 && n_atoms > n_solute {
            (n_atoms - n_solute) / atoms_per_solvent
        } else {
            0
        };
        let solvent_constraint_dof = if shake_enabled {
            n_solvent_molecules * topo.solvent_constraint_template.len()
        } else {
            0
        };
        let total_dof =
            (3 * n_atoms - solvent_constraint_dof) as f64 - imd.ndfmin as f64;
        md_sequence.push(Box::new(BerendsenThermostat::new_single_bath(
            temperature,
            thermostat_tau,
            total_dof,
            n_atoms,
        )));
    }

    // 4. Leap-Frog position
    md_sequence.push(Box::new(LeapFrogPosition::new()));

    // 5. SHAKE constraints
    if shake_enabled {
        let ntc_mode = match imd.ntc {
            3 => NtcMode::AllBonds,
            2 => NtcMode::HydrogenBonds,
            _ => NtcMode::SolventOnly,
        };
        let mut shake_alg = ShakeAlgorithm::new(ShakeParameters {
            tolerance: shake_tolerance,
            max_iterations: 1000,
            ntc: ntc_mode,
        });
        if imd.ntishk >= 1 {
            shake_alg.shake_initial_positions = true;
        }
        if imd.ntishk >= 2 {
            shake_alg.shake_initial_velocities = true;
        }
        md_sequence.push(Box::new(shake_alg));
    }

    // 6. Temperature calculation
    md_sequence.push(Box::new(TemperatureCalculation::new()));

    // 7. Pressure calculation and barostat
    if imd.couple_pressure {
        let virial_type = match imd
            .pressure_parameters
            .as_ref()
            .map(|p| p.virial)
            .unwrap_or(0)
        {
            2 => VirialType::Molecular,
            1 => VirialType::Atomic,
            _ => VirialType::None,
        };
        md_sequence.push(Box::new(PressureCalculation::new(virial_type)));

        let pp = imd.pressure_parameters.as_ref();
        md_sequence.push(Box::new(BerendsenBarostat::new(BerendsenBarostatParams {
            pressure0: pp.map(|p| p.pressure0[0][0]).unwrap_or(1.0),
            compressibility: pp.map(|p| p.compressibility[0][0]).unwrap_or(4.575e-4),
            tau: pp.map(|p| p.tau_p).unwrap_or(0.5),
        })));
    }

    // 8. Energy calculation
    md_sequence.push(Box::new(EnergyCalculation::new()));

    // Initialize sequence
    let mut sim_state = SimulationState::new(dt, n_steps);
    md_sequence
        .init(&topo, &mut conf, &sim_state)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to initialize algorithm sequence: {}",
                e
            ))
        })?;

    // Run step 0 (initial force evaluation, gromosXX convention)
    md_sequence
        .run_step(&topo, &mut conf, &sim_state)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Error at step 0: {}",
                e
            ))
        })?;
    sim_state.advance();

    Ok(PySimulation {
        topology: topo,
        configuration: conf,
        md_sequence,
        sim_state,
        dt,
        n_atoms,
    })
}

// ============================================================================
// PySimulation
// ============================================================================

/// A GROMOS-RS molecular dynamics simulation.
///
/// Interactive simulation object inspired by OpenMM's Simulation class,
/// using gromosXX naming conventions and file formats.
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
}

#[pymethods]
impl PySimulation {
    /// Create a new simulation.
    ///
    /// Accepts either file paths (str) or pre-loaded objects:
    ///
    ///     # From objects (compositional)
    ///     sim = Simulation(topology, configuration, parameters)
    ///
    ///     # From file paths (convenience)
    ///     sim = Simulation("system.topo", "initial.cnf", "run.imd")
    #[new]
    fn new(
        arg1: &Bound<'_, PyAny>,
        arg2: &Bound<'_, PyAny>,
        arg3: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        // String path API
        if let (Ok(topo_file), Ok(conf_file), Ok(input_file)) = (
            arg1.extract::<String>(),
            arg2.extract::<String>(),
            arg3.extract::<String>(),
        ) {
            return Self::_from_files(&topo_file, &conf_file, &input_file);
        }

        // Object API
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
            conf.pos_data.clone(),
            conf.vel_data.clone(),
            conf.box_dims,
            &params.inner,
        )
    }

    /// Create a simulation from file paths (alternative to constructor).
    ///
    /// Example:
    ///     sim = Simulation.from_files("system.topo", "initial.cnf", "run.imd")
    #[staticmethod]
    fn from_files(topo_file: &str, conf_file: &str, input_file: &str) -> PyResult<Self> {
        Self::_from_files(topo_file, conf_file, input_file)
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
    /// recently completed step's data (gromosXX convention).
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = state
            .pos
            .iter()
            .map(|v| vec![v.x, v.y, v.z])
            .collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Current velocities as an Nx3 numpy array (nm/ps).
    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = state
            .vel
            .iter()
            .map(|v| vec![v.x, v.y, v.z])
            .collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Current forces as an Nx3 numpy array (kJ/mol/nm).
    #[getter]
    fn forces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = state
            .force
            .iter()
            .map(|v| vec![v.x, v.y, v.z])
            .collect();
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
            angle: 0.0,
            dihedral: 0.0,
            improper: 0.0,
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
    #[getter]
    fn temperature(&self) -> f64 {
        let n_dof = self.n_atoms * 3;
        self.configuration.old().temperature(n_dof)
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
    fn _from_files(topo_file: &str, conf_file: &str, input_file: &str) -> PyResult<Self> {
        let imd = read_imd_file(input_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read input file '{}': {}",
                input_file, e
            ))
        })?;

        let topo_data = read_topology_file(topo_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read topology '{}': {}",
                topo_file, e
            ))
        })?;
        let topo = build_topology(topo_data);

        let coord_data = read_coordinates(conf_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read coordinates '{}': {}",
                conf_file, e
            ))
        })?;

        build_simulation(
            topo,
            coord_data.positions,
            coord_data.velocities,
            coord_data.box_dims,
            &imd,
        )
    }
}

/// Register simulation bindings in the Python module.
pub fn register_simulation(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySimulation>()?;
    Ok(())
}
