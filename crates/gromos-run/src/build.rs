//! Stages 2 and 3 of the builder: `instantiate` turns a validated plan into runnable
//! algorithms, `start` primes the run (PLAN.md 3.9).
//!
//! ```text
//! RunRecipe ── build_plan ──▶ Vec<AlgorithmSpec> ── validate_plan ── instantiate ──▶ AlgorithmSequence ── start
//! ```
//!
//! `instantiate` reads **only** the plan (plus the prepared topology/configuration and the
//! boundary conditions): every physical value it needs was resolved by `build_plan`, so an
//! edited plan is exactly what runs (drift guard G8). The recipe-aware wrappers below add the
//! reporting summary the binary prints.

use std::path::Path;

use gromos_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::{BoxType, Configuration},
    math::{Periodicity, Rectangular, Triclinic, Vacuum},
    pairlist::{PairlistAlgorithm, PairlistContainer},
    Topology,
};
use gromos_forces::{
    nonbonded::CRFParameters,
    restraints::{
        DistanceRestraint, DistanceRestraints, PerturbedDistanceRestraint,
        PerturbedDistanceRestraints, PositionRestraint, PositionRestraints,
    },
};
use gromos_integrators::{
    algorithms::{
        BerendsenBarostat, BerendsenBarostatParams, BerendsenThermostat, EnergyCalculation,
        Forcefield, LeapFrogPosition, LeapFrogVelocity, LincsAlgorithm, NoseHooverThermostat,
        PressureCalculation, RemoveCOMMotion, SettleAlgorithm, ShakeAlgorithm,
        SteepestDescentAlgorithm, TemperatureCalculation, VirialType,
    },
    constraints::ShakeParameters,
};
use gromos_io::{
    distanceres::read_distanceres,
    imd::ImdParameters,
    posres::{build_posres_entries, read_posresspec, read_refpos},
};

use crate::plan::{
    build_plan, validate_plan, AlgorithmSpec, DistanceRestraintsPlan, ForcefieldPlan,
    PositionRestraintsPlan,
};
use crate::recipe::{
    ConstraintAlgorithm, Diagnostics, RunRecipe, TermSpec, ThermostatAlgorithm, VirialKind,
};
use crate::{total_dof, ConstraintSelection, Prepared, RunError, RunInputs, RunOptions};

/// What the thermostat was built as — for the binary's report and the Python getters.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermostatSummary {
    /// `"Berendsen"`, `"Nose-Hoover"` or `"Nose-Hoover-Chain(n)"`.
    pub label: String,
    pub temperature: f64,
    pub tau: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BarostatSummary {
    pub pressure0: f64,
    pub tau: f64,
    pub compressibility: f64,
}

/// Facts about the built sequence that the caller may report. Never consulted by the
/// sequence itself.
#[derive(Debug, Clone)]
pub struct BuildSummary {
    pub is_minimization: bool,
    pub constraints: ConstraintSelection,
    /// SHAKE parameters when SHAKE is in the sequence.
    pub shake: Option<ShakeParameters>,
    /// LINCS expansion orders (solute, solvent).
    pub lincs_orders: (usize, usize),
    pub thermostat: Option<ThermostatSummary>,
    pub barostat: Option<BarostatSummary>,
    /// Number of position-restrained atoms.
    pub position_restraints: usize,
    /// (unperturbed, perturbed) distance restraints.
    pub distance_restraints: (usize, usize),
    /// Pairs in the initial pairlist.
    pub initial_pairs: usize,
    pub parallel_nonbonded: bool,
    /// Kinetic degrees of freedom (see [`crate::total_dof`]); what the thermostat couples to
    /// and what the temperature must be reported with.
    pub total_dof: f64,
}

/// The result of `instantiate`: the sequence plus the counts the caller may report.
pub struct Instantiated {
    pub sequence: AlgorithmSequence,
    pub initial_pairs: usize,
    pub position_restraints: usize,
    pub distance_restraints: (usize, usize),
}

/// A built run: the sequence, the data it came from, and the report.
pub struct Built {
    pub sequence: AlgorithmSequence,
    pub recipe: RunRecipe,
    pub plan: Vec<AlgorithmSpec>,
    pub diagnostics: Diagnostics,
    pub summary: BuildSummary,
}

/// The periodicity a prepared system runs under: triclinic for a truncated-octahedron box,
/// vacuum for an all-zero box, rectangular otherwise.
pub fn periodicity_of(prepared: &Prepared) -> Periodicity {
    let b = prepared.box_dims;
    if let Some(triclinic_box) = prepared.truncoct_box {
        Periodicity::Triclinic(Triclinic::new(triclinic_box))
    } else if b.x == 0.0 && b.y == 0.0 && b.z == 0.0 {
        Periodicity::Vacuum(Vacuum)
    } else {
        Periodicity::Rectangular(Rectangular::new(b))
    }
}

fn virial_type(kind: VirialKind) -> VirialType {
    match kind {
        VirialKind::Molecular => VirialType::Molecular,
        VirialKind::Atomic => VirialType::Atomic,
        VirialKind::None => VirialType::None,
    }
}

fn io_err(what: &str, path: &Path) -> impl FnOnce(gromos_io::IoError) -> RunError {
    let what = format!("{what} '{}'", path.display());
    move |source| RunError::Io { what, source }
}

fn load_position_restraints(
    plan: &PositionRestraintsPlan,
    startup_positions: &[gromos_core::math::Vec3],
) -> Result<PositionRestraints, RunError> {
    let restrained_atoms =
        read_posresspec(&plan.spec).map_err(io_err("position restraint spec", &plan.spec))?;
    // NTPORB=1: reference positions from a file; NTPORB=0: the start-up configuration.
    let ref_positions = match &plan.reference {
        Some(path) => read_refpos(path).map_err(io_err("reference positions", path))?,
        None => startup_positions.to_vec(),
    };
    let entries = build_posres_entries(&restrained_atoms, &ref_positions);
    let mut pr = PositionRestraints::new();
    for entry in &entries {
        pr.add_restraint(PositionRestraint::new(
            entry.atom,
            entry.reference_pos,
            plan.k,
        ));
    }
    Ok(pr)
}

fn load_distance_restraints(
    plan: &DistanceRestraintsPlan,
) -> Result<(DistanceRestraints, PerturbedDistanceRestraints), RunError> {
    let (unperturbed, perturbed) =
        read_distanceres(&plan.file).map_err(io_err("distance restraint file", &plan.file))?;
    let mut dr = DistanceRestraints::new();
    for spec in &unperturbed {
        let mut r = DistanceRestraint::new(
            spec.atom1,
            spec.atom2,
            spec.r0,
            spec.w0,
            spec.rah,
            plan.k,
            plan.r_linear,
            plan.mode,
        );
        if let Some(tau) = plan.time_averaging {
            r = r.with_time_averaging(tau, plan.dt);
        }
        dr.add(r);
    }
    let mut pdr = PerturbedDistanceRestraints::new();
    for spec in &perturbed {
        pdr.add(PerturbedDistanceRestraint::new(
            spec.atom1,
            spec.atom2,
            spec.n,
            spec.m,
            spec.a_r0,
            spec.b_r0,
            spec.a_w0,
            spec.b_w0,
            spec.rah,
            plan.k,
            plan.r_linear,
            plan.mode,
        ));
    }
    Ok((dr, pdr))
}

/// GROMOS: individual λ = RLAM^NLAM, dλ/dRLAM = NLAM · RLAM^(NLAM−1)
/// (same arithmetic as `ImdParameters::lambda_and_derivative`, bit for bit).
fn lambda_and_derivative(lambda: f64, nlam: i32) -> (f64, f64) {
    let n = nlam as f64;
    let l = lambda.powf(n);
    let dl = if nlam <= 0 {
        0.0
    } else {
        n * lambda.powf(n - 1.0)
    };
    (l, dl)
}

struct ForcefieldBuilt {
    forcefield: Forcefield,
    initial_pairs: usize,
    position_restraints: usize,
    distance_restraints: (usize, usize),
}

fn instantiate_forcefield(
    ff: &ForcefieldPlan,
    topo: &Topology,
    conf: &Configuration,
    periodicity: &Periodicity,
) -> Result<ForcefieldBuilt, RunError> {
    let n_atoms = topo.num_atoms();
    // CRF interior dielectric is always 1 in GROMOS.
    let crf_params = CRFParameters::new(ff.cutoff_long, 1.0, ff.epsilon_rf, ff.kappa);
    let lj_params = Forcefield::convert_lj_parameters(topo);

    let mut pairlist = PairlistContainer::new(ff.cutoff_short, ff.cutoff_long, 0.0);
    pairlist.update_frequency = ff.pairlist_every;
    pairlist.grid_size = ff.grid_size;
    let box_type = match periodicity {
        Periodicity::Rectangular(_) => BoxType::Rectangular,
        Periodicity::Triclinic(_) => BoxType::Triclinic,
        Periodicity::Vacuum(_) => BoxType::Vacuum,
    };
    let pairlist_algorithm = PairlistAlgorithm::from_imd(
        ff.pairlist_algorithm,
        n_atoms,
        box_type,
        !topo.chargegroups.is_empty(),
        ff.pairlist_type,
    );
    pairlist_algorithm.update(topo, conf, &mut pairlist, periodicity);
    let initial_pairs = pairlist.total_pairs();

    let mut forcefield = Forcefield::new(
        lj_params,
        crf_params,
        periodicity.clone(),
        pairlist,
        pairlist_algorithm,
    );
    forcefield.four_pi_eps_i = ff.four_pi_eps_i;
    forcefield.ntf_bond = ff.bonds;
    forcefield.ntf_angle = ff.angles;
    forcefield.ntf_improper = ff.impropers;
    forcefield.ntf_dihedral = ff.dihedrals;
    forcefield.parallel_nonbonded = ff.parallel;
    forcefield.atoms_per_solvent = ff.atoms_per_solvent;
    forcefield.virial_type = virial_type(ff.virial);

    let mut position_restraints = 0;
    if let Some(plan) = &ff.position_restraints {
        let pr = load_position_restraints(plan, &conf.current().pos)?;
        position_restraints = pr.restraints.len();
        forcefield.position_restraints = pr;
    }
    let mut distance_restraints = (0, 0);
    if let Some(plan) = &ff.distance_restraints {
        let (dr, pdr) = load_distance_restraints(plan)?;
        distance_restraints = (dr.restraints.len(), pdr.restraints.len());
        forcefield.distance_restraints = dr;
        forcefield.perturbed_distance_restraints = pdr;
        forcefield.lambda = plan.lambda;
    }
    if let Some(p) = &ff.perturbation {
        forcefield.lambda_and_derivative = lambda_and_derivative(p.lambda, p.nlam);
        forcefield.lambda_exp = p.nlam;
        forcefield.global_alphlj = p.alpha_lj;
        forcefield.global_alphc = p.alpha_c;
    }

    Ok(ForcefieldBuilt {
        forcefield,
        initial_pairs,
        position_restraints,
        distance_restraints,
    })
}

#[cfg(feature = "ml")]
fn instantiate_orchestrator(
    terms: &[TermSpec],
    topo: &Topology,
    periodicity: &Periodicity,
) -> Result<Box<dyn gromos_core::algorithm::Algorithm>, RunError> {
    use gromos_core::selection::AtomSelection;
    use gromos_forces::nonbonded::schnet::SchNetInteraction;
    use gromos_forces::orchestrator::ProviderOrchestrator;
    use gromos_forces::orchestrator_algorithm::ProviderOrchestratorAlgorithm;
    use gromos_forces::zones::{Zone, ZonePartition};

    let n_atoms = topo.num_atoms();
    let mut orchestrator = ProviderOrchestrator::new();
    for term in terms {
        match term {
            TermSpec::Schnet {
                model,
                cutoff,
                elements,
                region,
                buffer,
                ..
            } => {
                let partition = ZonePartition::from_selections(topo, region, buffer.as_deref())
                    .map_err(|e| RunError::Recipe(format!("schnet region: {e}")))?;
                let inner = partition.atoms_in(Zone::Inner);
                let selection = AtomSelection::from_indices(inner, n_atoms)
                    .map_err(|e| RunError::Recipe(format!("schnet region: {e}")))?;
                let provider = SchNetInteraction::load(model, *cutoff, elements.clone())
                    .map_err(|e| RunError::Recipe(format!("schnet model '{model}': {e}")))?;
                orchestrator.register(selection, Box::new(provider));
            },
        }
    }
    Ok(Box::new(ProviderOrchestratorAlgorithm::new(
        orchestrator,
        periodicity.clone(),
    )))
}

#[cfg(not(feature = "ml"))]
fn instantiate_orchestrator(
    terms: &[TermSpec],
    _topo: &Topology,
    _periodicity: &Periodicity,
) -> Result<Box<dyn gromos_core::algorithm::Algorithm>, RunError> {
    for term in terms {
        if let Some(feature) = term.feature() {
            return Err(RunError::MissingFeature {
                term: term.name().to_string(),
                feature,
            });
        }
    }
    Err(RunError::InvalidPlan(
        "orchestrator without terms".to_string(),
    ))
}

/// Stage 3: construct the algorithms of a **validated** plan.
///
/// Reads only the plan, the prepared topology/configuration and the boundary conditions —
/// no recipe, no parameter file (drift guard G8).
pub fn instantiate(
    plan: &[AlgorithmSpec],
    topo: &Topology,
    conf: &Configuration,
    periodicity: &Periodicity,
) -> Result<Instantiated, RunError> {
    let n_atoms = topo.num_atoms();
    let mut seq = AlgorithmSequence::new();
    let mut initial_pairs = 0;
    let mut position_restraints = 0;
    let mut distance_restraints = (0, 0);

    for spec in plan {
        match spec {
            AlgorithmSpec::RemoveCom { initial, every } => {
                seq.push(Box::new(RemoveCOMMotion::new(*initial, *every)));
            },
            AlgorithmSpec::Forcefield(ff) => {
                let built = instantiate_forcefield(ff, topo, conf, periodicity)?;
                initial_pairs = built.initial_pairs;
                position_restraints = built.position_restraints;
                distance_restraints = built.distance_restraints;
                seq.push(Box::new(built.forcefield));
            },
            AlgorithmSpec::Orchestrator { terms } => {
                seq.push(instantiate_orchestrator(terms, topo, periodicity)?);
            },
            AlgorithmSpec::LeapFrogVelocity => seq.push(Box::new(LeapFrogVelocity::new())),
            AlgorithmSpec::Thermostat {
                algorithm,
                temperature,
                tau,
                dof,
            } => match algorithm {
                ThermostatAlgorithm::Berendsen => seq.push(Box::new(
                    BerendsenThermostat::new_single_bath(*temperature, *tau, *dof, n_atoms),
                )),
                ThermostatAlgorithm::NoseHoover => seq.push(Box::new(
                    NoseHooverThermostat::new_single_bath(*temperature, *tau, *dof, n_atoms),
                )),
                ThermostatAlgorithm::NoseHooverChain(n) => {
                    seq.push(Box::new(NoseHooverThermostat::new_chain_bath(
                        *temperature,
                        *tau,
                        *dof,
                        n_atoms,
                        (*n).max(2),
                    )))
                },
            },
            AlgorithmSpec::LeapFrogPosition => seq.push(Box::new(LeapFrogPosition::new())),
            AlgorithmSpec::Shake {
                tolerance,
                max_iterations,
                solute,
                include_solvent,
                initial_positions,
                initial_velocities,
            } => {
                let mut alg = ShakeAlgorithm::new(ShakeParameters {
                    tolerance: *tolerance,
                    max_iterations: *max_iterations,
                    ntc: ConstraintSelection::ntc_mode_of(solute.ntc()),
                });
                alg.include_solvent = *include_solvent;
                alg.shake_initial_positions = *initial_positions;
                alg.shake_initial_velocities = *initial_velocities;
                seq.push(Box::new(alg));
            },
            AlgorithmSpec::Settle => seq.push(Box::new(SettleAlgorithm::new())),
            AlgorithmSpec::Lincs {
                solute,
                order_solute,
                order_solvent,
                constrain_solute,
                constrain_solvent,
            } => seq.push(Box::new(LincsAlgorithm::new(
                ConstraintSelection::ntc_mode_of(solute.ntc()),
                *order_solute,
                *order_solvent,
                *constrain_solute,
                *constrain_solvent,
            ))),
            AlgorithmSpec::SteepestDescent {
                tolerance,
                step0,
                step_max,
                min_steps,
                force_limit,
            } => seq.push(Box::new(
                SteepestDescentAlgorithm::new()
                    .with_tolerance(*tolerance)
                    .with_step_sizes(*step0, *step_max)
                    .with_min_steps(*min_steps)
                    .with_force_limit(*force_limit),
            )),
            AlgorithmSpec::TemperatureCalculation => {
                seq.push(Box::new(TemperatureCalculation::new()))
            },
            AlgorithmSpec::PressureCalculation { virial } => {
                seq.push(Box::new(PressureCalculation::new(virial_type(*virial))))
            },
            AlgorithmSpec::Barostat {
                pressure0,
                compressibility,
                tau,
            } => seq.push(Box::new(BerendsenBarostat::new(BerendsenBarostatParams {
                pressure0: *pressure0,
                compressibility: *compressibility,
                tau: *tau,
            }))),
            AlgorithmSpec::EnergyCalculation => seq.push(Box::new(EnergyCalculation::new())),
        }
    }

    Ok(Instantiated {
        sequence: seq,
        initial_pairs,
        position_restraints,
        distance_restraints,
    })
}

fn summarize(
    recipe: &RunRecipe,
    plan: &[AlgorithmSpec],
    topo: &Topology,
    inst: &Instantiated,
) -> BuildSummary {
    let has_solvent = topo.num_atoms() > topo.num_solute_atoms();
    let c = &recipe.constraints;
    let sel = ConstraintSelection::from_codes(
        c.solute.ntc(),
        ConstraintAlgorithm::code(c.solute_algorithm),
        ConstraintAlgorithm::code(c.solvent_algorithm),
        has_solvent,
    );
    let mut summary = BuildSummary {
        is_minimization: recipe.is_minimization(),
        constraints: sel,
        shake: None,
        lincs_orders: (c.lincs_order_solute, c.lincs_order_solvent),
        thermostat: None,
        barostat: None,
        position_restraints: inst.position_restraints,
        distance_restraints: inst.distance_restraints,
        initial_pairs: inst.initial_pairs,
        parallel_nonbonded: false,
        total_dof: total_dof(
            topo,
            &sel,
            ConstraintSelection::ntc_mode_of(c.solute.ntc()),
            recipe.boundary.ndfmin,
        ),
    };
    for spec in plan {
        match spec {
            AlgorithmSpec::Forcefield(ff) => summary.parallel_nonbonded = ff.parallel,
            AlgorithmSpec::Shake {
                tolerance,
                max_iterations,
                solute,
                ..
            } => {
                summary.shake = Some(ShakeParameters {
                    tolerance: *tolerance,
                    max_iterations: *max_iterations,
                    ntc: ConstraintSelection::ntc_mode_of(solute.ntc()),
                })
            },
            AlgorithmSpec::Thermostat {
                algorithm,
                temperature,
                tau,
                ..
            } => {
                summary.thermostat = Some(ThermostatSummary {
                    label: match algorithm {
                        ThermostatAlgorithm::Berendsen => "Berendsen".to_string(),
                        ThermostatAlgorithm::NoseHoover => "Nose-Hoover".to_string(),
                        ThermostatAlgorithm::NoseHooverChain(n) => {
                            format!("Nose-Hoover-Chain({n})")
                        },
                    },
                    temperature: *temperature,
                    tau: *tau,
                })
            },
            AlgorithmSpec::Barostat {
                pressure0,
                compressibility,
                tau,
            } => {
                summary.barostat = Some(BarostatSummary {
                    pressure0: *pressure0,
                    tau: *tau,
                    compressibility: *compressibility,
                })
            },
            _ => {},
        }
    }
    summary
}

/// Build a run from a recipe: plan → validate → instantiate, plus the report.
pub fn build_sequence_from_recipe(
    recipe: &RunRecipe,
    prepared: &Prepared,
    inputs: &RunInputs,
) -> Result<Built, RunError> {
    let topo = &prepared.topology;
    let plan = build_plan(
        recipe,
        topo,
        inputs,
        prepared.physical_constants.four_pi_eps_i,
    )?;
    validate_plan(&plan)?;
    let periodicity = periodicity_of(prepared);
    let inst = instantiate(&plan, topo, &prepared.configuration, &periodicity)?;
    let summary = summarize(recipe, &plan, topo, &inst);
    Ok(Built {
        sequence: inst.sequence,
        recipe: recipe.clone(),
        plan,
        diagnostics: Diagnostics::default(),
        summary,
    })
}

/// Build a run from GROMOS parameters: the `.imd` front-end of [`build_sequence_from_recipe`].
pub fn build_sequence_from_imd(
    imd: &ImdParameters,
    prepared: &Prepared,
    inputs: &RunInputs,
    options: &RunOptions,
) -> Result<Built, RunError> {
    let (mut recipe, diagnostics) = RunRecipe::from_imd_with(imd, &options.passthrough)?;
    recipe.execution.parallel = options.parallel;
    let mut built = build_sequence_from_recipe(&recipe, prepared, inputs)?;
    built.diagnostics = diagnostics;
    Ok(built)
}

/// Initialise the sequence and run step 0 (the initial force evaluation — GROMOS convention).
///
/// Shared so that the binary and the binding start a run identically: `init` is physics
/// (NTISHK shakes the initial positions/velocities), and step 0 primes the double-buffered
/// state. The caller advances `state` afterwards.
pub fn start(
    sequence: &mut AlgorithmSequence,
    topo: &Topology,
    conf: &mut Configuration,
    state: &SimulationState,
) -> Result<(), RunError> {
    sequence
        .init(topo, conf, state)
        .map_err(|e| RunError::Init(e.to_string()))?;
    sequence
        .run_step(topo, conf, state)
        .map_err(|e| RunError::Step {
            step: state.step,
            message: e.to_string(),
        })
}

#[cfg(test)]
mod tests {
    //! The recipe path must be **bit-identical** to the step-1 builder it replaced, on every
    //! reference system (the Rust-side A/B PLAN.md 3.9 step 2 promised).

    use std::path::{Path, PathBuf};

    use gromos_core::algorithm::SimulationState;
    use gromos_io::{
        coordinate::read_coordinates,
        imd::read_imd_file,
        topology::{build_topology, read_topology_file},
    };

    use super::*;
    use crate::{legacy_builder::legacy_build_sequence, prepare_system, Coordinates};

    fn refs() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
    }

    fn toml_str(text: &str, key: &str) -> Option<String> {
        text.lines().find_map(|l| {
            let l = l.trim();
            l.strip_prefix(&format!("{key} = "))
                .map(|v| v.trim().trim_matches('"').to_string())
        })
    }

    struct Loaded {
        imd: ImdParameters,
        inputs: RunInputs,
        prepare: Box<dyn Fn() -> Prepared>,
    }

    fn load(dir: &Path) -> Option<Loaded> {
        let toml = std::fs::read_to_string(dir.join("input.toml")).ok()?;
        let path = |key: &str| toml_str(&toml, key).map(|v| dir.join(v));
        let topo_path = path("topology")?;
        let conf_path = path("configuration")?;
        let imd = read_imd_file(path("parameters")?).ok()?;
        let inputs = RunInputs {
            pttopo: path("pttopo"),
            posresspec: path("posresspec"),
            refpos: path("refpos"),
            distrest: path("distrest"),
        };
        let imd2 = imd.clone();
        let inputs2 = inputs.clone();
        let prepare = Box::new(move || {
            let parsed = read_topology_file(&topo_path).expect("topology");
            let pc = parsed.physical_constants;
            let topo = build_topology(parsed);
            let coords: Coordinates = read_coordinates(&conf_path).expect("conf").into();
            prepare_system(&imd2, topo, pc, coords, &inputs2).expect("prepare_system")
        });
        Some(Loaded {
            imd,
            inputs,
            prepare,
        })
    }

    /// Run `steps` steps and return every number a parity comparison cares about.
    fn run(mut prepared: Prepared, mut seq: AlgorithmSequence, dt: f64, steps: usize) -> Vec<f64> {
        let topo = prepared.topology;
        let conf = &mut prepared.configuration;
        let mut state = SimulationState::new(dt, steps);
        start(&mut seq, &topo, conf, &state).expect("start");
        state.advance();
        for _ in 0..steps {
            seq.run_step(&topo, conf, &state).expect("run_step");
            state.advance();
        }
        let old = conf.old();
        let mut out = vec![
            old.energies.total(),
            old.energies.kinetic_total,
            old.energies.potential_total,
        ];
        out.extend(old.pos.iter().flat_map(|v| [v.x, v.y, v.z]));
        out.extend(old.force.iter().flat_map(|v| [v.x, v.y, v.z]));
        out
    }

    #[test]
    fn plan_matches_legacy_builder_bit_for_bit() {
        let mut dirs: Vec<PathBuf> = std::fs::read_dir(refs())
            .expect("reference dir")
            .flatten()
            .map(|e| e.path())
            .filter(|p| p.is_dir() && p.join("input.toml").exists())
            .collect();
        dirs.sort();
        assert!(dirs.len() >= 40);
        let options = RunOptions::default();
        let mut compared = 0;
        for dir in dirs {
            let name = dir.file_name().unwrap().to_string_lossy().to_string();
            let Some(loaded) = load(&dir) else { continue };
            let steps = loaded.imd.nstlim.min(10);
            let dt = loaded.imd.dt;

            let p_legacy = (loaded.prepare)();
            let legacy = legacy_build_sequence(&loaded.imd, &p_legacy, &loaded.inputs, &options)
                .unwrap_or_else(|e| panic!("{name}: legacy builder: {e}"));
            let legacy_names: Vec<String> = legacy
                .sequence
                .algorithm_names()
                .iter()
                .map(|s| s.to_string())
                .collect();
            let a = run(p_legacy, legacy.sequence, dt, steps);

            let p_new = (loaded.prepare)();
            let built = build_sequence_from_imd(&loaded.imd, &p_new, &loaded.inputs, &options)
                .unwrap_or_else(|e| panic!("{name}: recipe builder: {e}"));
            let new_names: Vec<String> = built
                .sequence
                .algorithm_names()
                .iter()
                .map(|s| s.to_string())
                .collect();
            assert_eq!(new_names, legacy_names, "{name}: sequences differ");
            let b = run(p_new, built.sequence, dt, steps);

            assert_eq!(a.len(), b.len(), "{name}");
            for (i, (x, y)) in a.iter().zip(&b).enumerate() {
                assert!(
                    x.to_bits() == y.to_bits(),
                    "{name}: value {i} differs after {steps} steps: {x:e} vs {y:e}"
                );
            }
            compared += 1;
        }
        assert!(compared >= 40, "only {compared} systems compared");
    }
}
