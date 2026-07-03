//! Interactive Simulation API for Python
//!
//! Supports both compositional and file-based construction, mirroring
//! the internal md architecture with GROMOS naming conventions.
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

use numpy::PyArray2;
use pyo3::prelude::*;

use gromos_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::BoxType,
    configuration::{Box as SimBox, Configuration},
    math::{
        truncoct_triclinic, truncoct_triclinic_box, truncoct_triclinic_rotmat, Mat3, Periodicity,
        Rectangular, Triclinic, Vacuum, Vec3,
    },
    pairlist::{PairlistAlgorithm, PairlistContainer},
    random::generate_velocities,
    Topology,
};
use gromos_forces::nonbonded::CRFParameters;
use gromos_forces::restraints::{
    DistanceRestraint, DistanceRestraints, PerturbedDistanceRestraint,
    PerturbedDistanceRestraints, PositionRestraint, PositionRestraints,
};
use gromos_integrators::algorithms::{
    BerendsenBarostat, BerendsenBarostatParams, BerendsenThermostat, EnergyCalculation, Forcefield,
    LeapFrogPosition, LeapFrogVelocity, LincsAlgorithm, NoseHooverThermostat, PressureCalculation,
    RemoveCOMMotion, SettleAlgorithm, ShakeAlgorithm, SteepestDescentAlgorithm,
    TemperatureCalculation, VirialType,
};
use gromos_integrators::constraints::{NtcMode, ShakeParameters};
use gromos_io::{
    coordinate::read_coordinates,
    distanceres::read_distanceres,
    imd::{read_imd_file, ImdParameters},
    posres::{build_posres_entries, read_posresspec, read_refpos},
    topology::{build_topology, read_topology_file},
};

use super::algorithm_sequence::{
    compute_total_dof, resolve_algorithm_sequence, PyAlgorithmSequence,
};
use super::parameters::PyInputParameters;
use super::py_conf::PyConfiguration;
use super::system::PySystem;
use super::topology::PyTopology;
use super::PyEnergy;

// ============================================================================
// Shared build logic — constructs a fully initialized simulation from parts
// ============================================================================

/// Build a `ShakeAlgorithm` from IMD constraint settings — the single place
/// both the standard-MD and steepest-descent branches of `build_simulation`
/// construct SHAKE from, so they can't silently drift apart.
fn shake_algorithm_from_imd(imd: &ImdParameters, solute_shake: bool, solvent_shake: bool) -> ShakeAlgorithm {
    let ntc_mode = if solute_shake {
        match imd.ntc {
            3 => NtcMode::AllBonds,
            2 => NtcMode::HydrogenBonds,
            _ => NtcMode::SolventOnly,
        }
    } else {
        NtcMode::SolventOnly
    };
    let mut shake_alg = ShakeAlgorithm::new(ShakeParameters {
        tolerance: imd.shake_tol,
        max_iterations: 1000,
        ntc: ntc_mode,
    });
    shake_alg.include_solvent = solvent_shake;
    if imd.ntishk >= 1 {
        shake_alg.shake_initial_positions = true;
    }
    if imd.ntishk >= 2 {
        shake_alg.shake_initial_velocities = true;
    }
    shake_alg
}

/// Which constraint algorithm(s) GROMOS dispatches to, mirroring `md.rs`
/// exactly: the solute algorithm is chosen by NTCP (only relevant when
/// NTC>1) and the solvent algorithm by NTCS (only relevant when solvent
/// molecules exist) — independently of each other.
struct ConstraintSelection {
    solute_shake: bool,
    solvent_shake: bool,
    settle_enabled: bool,
    lincs_enabled: bool,
    solute_lincs: bool,
    solvent_lincs: bool,
}

impl ConstraintSelection {
    /// `has_solvent` should reflect the topology's *actual* solvent atoms
    /// (`num_atoms() > num_solute_atoms()`) rather than `imd.nsm` — the
    /// compositional Python path often solvates the topology directly via
    /// `Topology.solvate()` with factory-built `InputParameters` that never
    /// set `nsm`, so relying on `nsm` alone would silently disable solvent
    /// constraints (e.g. rigid water flying apart) whenever the two diverge.
    fn from_imd(imd: &ImdParameters, has_solvent: bool) -> Self {
        let solute_constrained = imd.ntc > 1;
        let solute_lincs = solute_constrained && imd.ntcp == 2;
        let solute_shake = solute_constrained && !solute_lincs;

        let solvent_constrained = imd.ntcs > 0 && has_solvent;
        let solvent_settle = solvent_constrained && imd.ntcs == 3;
        let solvent_lincs = solvent_constrained && imd.ntcs == 2;
        let solvent_shake = solvent_constrained && !solvent_settle && !solvent_lincs;

        ConstraintSelection {
            solute_shake,
            solvent_shake,
            settle_enabled: solvent_settle,
            lincs_enabled: solute_lincs || solvent_lincs,
            solute_lincs,
            solvent_lincs,
        }
    }

    fn shake_enabled(&self) -> bool {
        self.solute_shake || self.solvent_shake
    }
}

/// Push whichever constraint algorithm(s) `ConstraintSelection` selected
/// onto the sequence, in the same SHAKE→SETTLE→LINCS order `md.rs` uses.
fn push_constraint_algorithms(
    md_sequence: &mut AlgorithmSequence,
    imd: &ImdParameters,
    sel: &ConstraintSelection,
) {
    if sel.shake_enabled() {
        md_sequence.push(Box::new(shake_algorithm_from_imd(
            imd,
            sel.solute_shake,
            sel.solvent_shake,
        )));
    }
    if sel.settle_enabled {
        md_sequence.push(Box::new(SettleAlgorithm::new()));
    }
    if sel.lincs_enabled {
        let ntc_mode = match imd.ntc {
            3 => NtcMode::AllBonds,
            2 => NtcMode::HydrogenBonds,
            _ => NtcMode::SolventOnly,
        };
        md_sequence.push(Box::new(LincsAlgorithm::new(
            ntc_mode,
            imd.lincs_order_solute,
            imd.lincs_order_solvent,
            sel.solute_lincs,
            sel.solvent_lincs,
        )));
    }
}

/// Build the thermostat algorithm GROMOS's MULTIBATH algorithm selector
/// picks: 0 = Berendsen, 1 = single-bath Nose-Hoover, n>=2 = Nose-Hoover
/// chain of length n (or `nhc_chain`, whichever is larger) — mirrors
/// `md.rs`'s `thermostat_algorithm` dispatch.
fn push_thermostat(
    md_sequence: &mut AlgorithmSequence,
    imd: &ImdParameters,
    temperature: f64,
    tau: f64,
    total_dof: f64,
    n_atoms: usize,
) {
    let algorithm = imd.temp_bath.first().map(|b| b.algorithm).unwrap_or(0);
    let nhc_chain = imd.temp_bath.first().map(|b| b.nhc_chain).unwrap_or(0);
    match algorithm {
        0 => {
            md_sequence.push(Box::new(BerendsenThermostat::new_single_bath(
                temperature,
                tau,
                total_dof,
                n_atoms,
            )));
        },
        1 => {
            md_sequence.push(Box::new(NoseHooverThermostat::new_single_bath(
                temperature,
                tau,
                total_dof,
                n_atoms,
            )));
        },
        n => {
            let chain = (n as usize).max(nhc_chain).max(2);
            md_sequence.push(Box::new(NoseHooverThermostat::new_chain_bath(
                temperature,
                tau,
                total_dof,
                n_atoms,
                chain,
            )));
        },
    }
}

/// NTB=-1 (truncated octahedron): GROMOS stores the legacy "cube edge length
/// L" in the coordinate file and converts it to the lower-triangular
/// triclinic box vectors on read, rotating positions/velocities into that
/// frame — mirrors `md.rs`'s identical `truncoct_triclinic_box`/
/// `truncoct_triclinic` calls. Returns the triclinic box matrix (`None` for
/// any other NTB, i.e. rectangular or vacuum) and mutates `positions`/
/// `velocities` in place when it applies.
fn apply_truncoct_box(
    imd: &ImdParameters,
    box_dims: Vec3,
    positions: &mut [Vec3],
    velocities: &mut [Vec3],
) -> Option<Mat3> {
    if imd.ntb == -1 {
        let cubic = Mat3::from_diagonal(box_dims);
        let triclinic_box = truncoct_triclinic_box(cubic, true);
        truncoct_triclinic(positions, true);
        truncoct_triclinic(velocities, true);
        Some(triclinic_box)
    } else {
        None
    }
}

/// Optional restraint-input file paths, threaded through from the Python
/// `Simulation`/`System` constructors down to `build_simulation` — mirrors
/// the `md` binary's `@distrest`/`@posresspec`/`@refpos` CLI flags.
#[derive(Default, Clone)]
pub(crate) struct RestraintFiles {
    pub distrest: Option<String>,
    pub posresspec: Option<String>,
    pub refpos: Option<String>,
}

/// Load and apply position/distance restraints onto a `Forcefield`, mirroring
/// `md.rs`'s `@posresspec`/`@refpos`/`@distrest` handling exactly (same file
/// formats, same NTPOR/NTPORB/NTDIR dispatch) — the only gap being FEP's
/// perturbed distance restraints' `lambda` still comes from `imd.rlam` same
/// as `md.rs`.
#[allow(clippy::too_many_arguments)]
fn apply_restraints(
    forcefield: &mut Forcefield,
    imd: &ImdParameters,
    positions: &[Vec3],
    distrest_file: Option<&str>,
    posresspec_file: Option<&str>,
    refpos_file: Option<&str>,
) -> PyResult<()> {
    fn io_err(e: impl std::fmt::Display) -> PyErr {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string())
    }

    if imd.ntpor > 0 {
        let por_file = posresspec_file.ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "NTPOR={} but no posresspec file was given",
                imd.ntpor
            ))
        })?;
        let restrained_atoms = read_posresspec(por_file).map_err(io_err)?;
        let ref_positions = if imd.ntporb >= 1 {
            let rpr_file = refpos_file.ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "NTPORB={} but no refpos file was given",
                    imd.ntporb
                ))
            })?;
            read_refpos(rpr_file).map_err(io_err)?
        } else {
            positions.to_vec()
        };
        let entries = build_posres_entries(&restrained_atoms, &ref_positions);
        let mut pr = PositionRestraints::new();
        for entry in &entries {
            pr.add_restraint(PositionRestraint::new(
                entry.atom,
                entry.reference_pos,
                imd.cpor,
            ));
        }
        forcefield.position_restraints = pr;
    }

    if imd.ntdir != 0 {
        let dr_file = distrest_file.ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "NTDIR={} but no distrest file was given",
                imd.ntdir
            ))
        })?;
        let (unperturbed, perturbed) = read_distanceres(dr_file).map_err(io_err)?;
        let mode = imd.ntdir.abs();
        let k = imd.cdir;
        let r_linear = imd.dir0;

        let mut dr = DistanceRestraints::new();
        for spec in &unperturbed {
            let mut r = DistanceRestraint::new(
                spec.atom1, spec.atom2, spec.r0, spec.w0, spec.rah, k, r_linear, mode,
            );
            if imd.ntdir < 0 {
                r = r.with_time_averaging(imd.taudir, imd.dt);
            }
            dr.add(r);
        }

        let mut pdr = PerturbedDistanceRestraints::new();
        for spec in &perturbed {
            let r = PerturbedDistanceRestraint::new(
                spec.atom1, spec.atom2, spec.n, spec.m, spec.a_r0, spec.b_r0, spec.a_w0,
                spec.b_w0, spec.rah, k, r_linear, mode,
            );
            pdr.add(r);
        }

        forcefield.distance_restraints = dr;
        forcefield.perturbed_distance_restraints = pdr;
        forcefield.lambda = imd.rlam;
    }

    Ok(())
}

/// Resolve initial velocities the same way `md.rs` does: generate a
/// Maxwell-Boltzmann distribution at `imd.tempi` (seeded by `imd.ig`) when
/// `NTIVEL=1`, otherwise use the coordinate-file velocities (or zero if
/// absent/mismatched).
fn initial_velocities(imd: &ImdParameters, velocities: &[Vec3], masses: &[f64]) -> Vec<Vec3> {
    if imd.ntivel == 1 {
        generate_velocities(imd.tempi, imd.ig as u32, masses)
    } else if velocities.len() == masses.len() {
        velocities.to_vec()
    } else {
        vec![Vec3::ZERO; masses.len()]
    }
}

/// Build a simulation from raw components.
///
/// This is the shared core used by both the file-path and object constructors.
/// It solvates the topology (if not already solvated), validates atom counts,
/// builds the Configuration + AlgorithmSequence, and runs step 0.
fn build_simulation(
    mut topo: Topology,
    mut positions: Vec<Vec3>,
    mut velocities: Vec<Vec3>,
    box_dims: Vec3,
    imd: &ImdParameters,
    restraints: &RestraintFiles,
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

    let truncoct_box_matrix =
        apply_truncoct_box(imd, box_dims, &mut positions, &mut velocities);

    // Build Configuration (double-buffered state)
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = positions;
    conf.current_mut().vel = initial_velocities(imd, &velocities, &topo.mass);
    conf.current_mut().box_config = match truncoct_box_matrix {
        Some(triclinic_box) => SimBox::truncated_octahedral(triclinic_box),
        None => SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z),
    };
    conf.copy_current_to_old();

    // Extract parameters
    let dt = imd.dt;
    let n_steps = imd.nstlim;
    let cutoff = imd.rcutl;
    let rf_epsilon = imd.epsrf;
    let rf_kappa = imd.appak;
    let pairlist_update = imd.nsnb;
    let ntf = imd.ntf;

    let has_solvent = topo.num_atoms() > topo.num_solute_atoms();
    let constraints = ConstraintSelection::from_imd(imd, has_solvent);
    let is_minimization = imd.ntem > 0;

    // Computed once, used both by the thermostat (if present) and stored on
    // PySimulation for the `temperature` getter — the two must agree.
    let total_dof = compute_total_dof(&topo, imd);

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

    let periodicity = if let Some(triclinic_box) = truncoct_box_matrix {
        Periodicity::Triclinic(Triclinic::new(triclinic_box))
    } else if box_dims.x == 0.0 && box_dims.y == 0.0 && box_dims.z == 0.0 {
        Periodicity::Vacuum(Vacuum)
    } else {
        Periodicity::Rectangular(Rectangular::new(box_dims))
    };
    let box_type = match &periodicity {
        Periodicity::Rectangular(_) => BoxType::Rectangular,
        Periodicity::Triclinic(_) => BoxType::Triclinic,
        Periodicity::Vacuum(_) => BoxType::Vacuum,
    };
    let pairlist_algorithm = PairlistAlgorithm::from_imd(
        imd.algorithm,
        topo.num_atoms(),
        box_type,
        !topo.chargegroups.is_empty(),
    );

    pairlist_algorithm.update(&topo, &conf, &mut pairlist, &periodicity);

    // Build algorithm sequence — GROMOS Leap-Frog pattern, or Steepest-Descent
    // energy minimization if NTEM > 0. EM replaces the leap-frog velocity+position
    // steps and skips COM removal / thermostat / barostat / kinetic-energy
    // calculation entirely (GROMOS convention: E_total = E_pot during EM),
    // mirroring the `md` binary's own sequence construction.
    let mut md_sequence = AlgorithmSequence::new();

    // 1. COM motion removal (not applicable during EM — no velocities)
    if !is_minimization && (imd.nticom >= 1 || imd.nscm != 0) {
        md_sequence.push(Box::new(RemoveCOMMotion::new(imd.nticom, imd.nscm)));
    }

    // 2. Forcefield (bonded + nonbonded)
    let mut forcefield = Forcefield::new(
        lj_params,
        crf_params,
        periodicity,
        pairlist,
        pairlist_algorithm,
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
    apply_restraints(
        &mut forcefield,
        imd,
        &conf.current().pos,
        restraints.distrest.as_deref(),
        restraints.posresspec.as_deref(),
        restraints.refpos.as_deref(),
    )?;
    md_sequence.push(Box::new(forcefield));

    if is_minimization {
        // Steepest-descent minimization: replaces LeapFrogVelocity + LeapFrogPosition.
        let sd = SteepestDescentAlgorithm::new()
            .with_tolerance(imd.dele)
            .with_step_sizes(imd.dx0, imd.dxm)
            .with_min_steps(imd.nmin)
            .with_force_limit(imd.flim);
        md_sequence.push(Box::new(sd));

        // GROMOS applies constraints even during minimization.
        push_constraint_algorithms(&mut md_sequence, imd, &constraints);

        // No TemperatureCalculation — EM has no velocities/kinetic energy.
        md_sequence.push(Box::new(EnergyCalculation::new()));
    } else {
        // 3. Leap-Frog velocity
        md_sequence.push(Box::new(LeapFrogVelocity::new()));

        // 3b. Thermostat (Berendsen or Nose-Hoover/chain, per MULTIBATH algorithm)
        if thermostat_on {
            push_thermostat(&mut md_sequence, imd, temperature, thermostat_tau, total_dof, n_atoms);
        }

        // 4. Leap-Frog position
        md_sequence.push(Box::new(LeapFrogPosition::new()));

        // 5. Constraints (SHAKE / SETTLE / LINCS)
        push_constraint_algorithms(&mut md_sequence, imd, &constraints);

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
    }

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

    // Run step 0 (initial force evaluation, GROMOS convention)
    md_sequence
        .run_step(&topo, &mut conf, &sim_state)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Error at step 0: {}", e))
        })?;
    sim_state.advance();

    Ok(PySimulation {
        topology: topo,
        configuration: conf,
        md_sequence,
        sim_state,
        dt,
        n_atoms,
        total_dof,
    })
}

/// Build a simulation from a user-provided algorithm sequence.
///
/// The sequence descriptors are resolved into real Rust algorithms using
/// the topology, configuration, and IMD parameters for context.
fn build_simulation_from_sequence(
    mut topo: Topology,
    positions: Vec<Vec3>,
    velocities: Vec<Vec3>,
    box_dims: Vec3,
    imd: &ImdParameters,
    sequence: &PyAlgorithmSequence,
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
    conf.current_mut().vel = initial_velocities(imd, &velocities, &topo.mass);
    conf.current_mut().box_config = SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z);
    conf.copy_current_to_old();

    let dt = imd.dt;
    let n_steps = imd.nstlim;
    let total_dof = compute_total_dof(&topo, imd);

    // Resolve descriptors into real algorithms
    let mut md_sequence = resolve_algorithm_sequence(sequence, &topo, &conf, imd, box_dims)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to resolve algorithm sequence: {}",
                e
            ))
        })?;

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

    // Run step 0 (initial force evaluation, GROMOS convention)
    md_sequence
        .run_step(&topo, &mut conf, &sim_state)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Error at step 0: {}", e))
        })?;
    sim_state.advance();

    Ok(PySimulation {
        topology: topo,
        configuration: conf,
        md_sequence,
        sim_state,
        dt,
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
    /// at build time and reused for `temperature` — must match what the
    /// thermostat (if any) is coupling to. See `compute_total_dof`.
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
    #[new]
    #[pyo3(signature = (arg1, arg2, arg3=None, *, distrest=None, posresspec=None, refpos=None))]
    fn new(
        arg1: &Bound<'_, PyAny>,
        arg2: &Bound<'_, PyAny>,
        arg3: Option<&Bound<'_, PyAny>>,
        distrest: Option<String>,
        posresspec: Option<String>,
        refpos: Option<String>,
    ) -> PyResult<Self> {
        let restraints = RestraintFiles {
            distrest,
            posresspec,
            refpos,
        };
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
                    system.configuration.pos_data.clone(),
                    system.configuration.vel_data.clone(),
                    system.configuration.box_dims,
                    &params.inner,
                    &restraints,
                )
            },
            Some(arg3) => {
                // Three-arg forms: file paths or (Topology, Configuration, InputParameters)
                if let (Ok(topo_file), Ok(conf_file), Ok(input_file)) = (
                    arg1.extract::<String>(),
                    arg2.extract::<String>(),
                    arg3.extract::<String>(),
                ) {
                    return Self::_from_files_with_restraints(
                        &topo_file,
                        &conf_file,
                        &input_file,
                        &restraints,
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
                    conf.pos_data.clone(),
                    conf.vel_data.clone(),
                    conf.box_dims,
                    &params.inner,
                    &restraints,
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
    /// and pressure it needs the degrees-of-freedom count, which `PySimulation`
    /// does not currently track (see `BerendsenThermostat` construction in
    /// `build_simulation` for that computation).
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

    fn _from_files(topo_file: &str, conf_file: &str, input_file: &str) -> PyResult<Self> {
        Self::_from_files_with_restraints(
            topo_file,
            conf_file,
            input_file,
            &RestraintFiles::default(),
        )
    }

    fn _from_files_with_restraints(
        topo_file: &str,
        conf_file: &str,
        input_file: &str,
        restraints: &RestraintFiles,
    ) -> PyResult<Self> {
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
            restraints,
        )
    }
}

/// Register simulation bindings in the Python module.
pub fn register_simulation(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySimulation>()?;
    Ok(())
}
