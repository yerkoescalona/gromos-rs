//! `AlgorithmSpec` — the MD step as data (stage 1 of the builder, PLAN.md 3.9).
//!
//! `build_plan` turns a [`RunRecipe`] into the ordered list of algorithms the GROMOS step runs,
//! with **every value resolved** (degrees of freedom, virial kind, `four_pi_eps_i`, restraint
//! file paths, the parallel decision, …). `instantiate` (in `build.rs`) reads *only* this list —
//! never the recipe — so an edited plan is exactly what runs (drift guard G8).
//!
//! `validate_plan` enforces the GROMOS ordering invariants, declared per kind in
//! [`AlgorithmSpec::rules`] (drift guard G9); every plan goes through it, edited or not.
//!
//! Registries: [`AlgorithmSpec::name`] and [`TermSpec`]'s are exhaustive `match`es, so a new
//! variant without a name is a compile error; [`AlgorithmSpec::examples`] must list one value
//! per variant (`registry_is_complete` test), which is what the Python side introspects.

use std::path::PathBuf;

use gromos_core::Topology;
use serde::{Deserialize, Serialize};

use crate::recipe::{ConstraintAlgorithm, Coupling, SoluteConstraints, ThermostatAlgorithm};
use crate::recipe::{RunRecipe, TermSpec, VirialKind};
use crate::{total_dof, ConstraintSelection, RunError, RunInputs};

/// Position restraints as the plan carries them: the input files, the force constant, and
/// whether the reference positions are the start-up configuration (`reference: None`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PositionRestraintsPlan {
    pub spec: PathBuf,
    /// `None` = the start-up positions (NTPORB=0), resolved at instantiation from the
    /// configuration, so the plan stays the same run after run.
    pub reference: Option<PathBuf>,
    pub k: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DistanceRestraintsPlan {
    pub file: PathBuf,
    /// |NTDIR|
    pub mode: i32,
    pub k: f64,
    pub r_linear: f64,
    /// `Some(tau)` for time-averaged restraints (NTDIR < 0); needs `dt`.
    pub time_averaging: Option<f64>,
    pub dt: f64,
    pub lambda: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PerturbationPlan {
    pub lambda: f64,
    pub nlam: i32,
    pub alpha_lj: f64,
    pub alpha_c: f64,
}

/// Everything `Forcefield` needs, resolved.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ForcefieldPlan {
    pub cutoff_short: f64,
    pub cutoff_long: f64,
    pub pairlist_every: usize,
    pub grid_size: f64,
    /// PAIRLIST ALGORITHM code (0 standard, 1 grid, 2 grid_cell) — `PairlistAlgorithm::from_imd`
    pub pairlist_algorithm: i32,
    /// PAIRLIST TYPE code (0 chargegroup, 1 atomic)
    pub pairlist_type: i32,
    pub epsilon_rf: f64,
    pub kappa: f64,
    pub four_pi_eps_i: f64,
    pub bonds: bool,
    pub angles: bool,
    pub impropers: bool,
    pub dihedrals: bool,
    pub parallel: bool,
    pub atoms_per_solvent: usize,
    pub virial: VirialKind,
    pub position_restraints: Option<PositionRestraintsPlan>,
    pub distance_restraints: Option<DistanceRestraintsPlan>,
    pub perturbation: Option<PerturbationPlan>,
}

/// One algorithm of the MD step, fully resolved. One spec = one `Algorithm`.
///
/// `Forcefield` is much larger than the unit variants; that is the nature of the data (it is
/// what the force field needs) and a plan has one of them, so it is not boxed.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum AlgorithmSpec {
    RemoveCom {
        initial: i32,
        every: i32,
    },
    Forcefield(ForcefieldPlan),
    Orchestrator {
        terms: Vec<TermSpec>,
    },
    LeapFrogVelocity,
    Thermostat {
        algorithm: ThermostatAlgorithm,
        temperature: f64,
        tau: f64,
        /// Kinetic degrees of freedom the bath couples to (resolved from the constraints).
        dof: f64,
    },
    LeapFrogPosition,
    Shake {
        tolerance: f64,
        max_iterations: usize,
        solute: SoluteConstraints,
        include_solvent: bool,
        initial_positions: bool,
        initial_velocities: bool,
    },
    Settle,
    Lincs {
        solute: SoluteConstraints,
        order_solute: usize,
        order_solvent: usize,
        constrain_solute: bool,
        constrain_solvent: bool,
    },
    SteepestDescent {
        tolerance: f64,
        step0: f64,
        step_max: f64,
        min_steps: usize,
        force_limit: f64,
    },
    TemperatureCalculation,
    PressureCalculation {
        virial: VirialKind,
    },
    Barostat {
        pressure0: f64,
        compressibility: f64,
        tau: f64,
    },
    EnergyCalculation,
}

/// Ordering rules of one kind, consulted by `validate_plan`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KindRules {
    /// At most one of this kind.
    pub unique: bool,
    /// Must appear exactly once.
    pub required: bool,
    /// Must be the last algorithm.
    pub last: bool,
    /// Must be the first algorithm, if present.
    pub first: bool,
    /// Must come immediately after this kind.
    pub must_follow: Option<&'static str>,
    /// Must appear somewhere after every one of these kinds.
    pub after: &'static [&'static str],
    /// Must appear somewhere before every one of these kinds.
    pub before: &'static [&'static str],
    /// Cannot coexist with these kinds.
    pub excludes: &'static [&'static str],
    /// Requires these kinds to be present.
    pub requires: &'static [&'static str],
}

const NONE: KindRules = KindRules {
    unique: false,
    required: false,
    last: false,
    first: false,
    must_follow: None,
    after: &[],
    before: &[],
    excludes: &[],
    requires: &[],
};

const DYNAMICS_ONLY: &[&str] = &[
    "leap_frog_velocity",
    "leap_frog_position",
    "thermostat",
    "temperature_calculation",
    "pressure_calculation",
    "barostat",
    "remove_com",
];

impl AlgorithmSpec {
    /// Every kind, in canonical order. `examples()` and the tests keep this list honest.
    pub const KINDS: &'static [&'static str] = &[
        "remove_com",
        "forcefield",
        "orchestrator",
        "leap_frog_velocity",
        "thermostat",
        "leap_frog_position",
        "shake",
        "settle",
        "lincs",
        "steepest_descent",
        "temperature_calculation",
        "pressure_calculation",
        "barostat",
        "energy_calculation",
    ];

    /// The serde tag of this spec — exhaustive, so a new variant cannot be nameless.
    pub fn name(&self) -> &'static str {
        match self {
            AlgorithmSpec::RemoveCom { .. } => "remove_com",
            AlgorithmSpec::Forcefield(_) => "forcefield",
            AlgorithmSpec::Orchestrator { .. } => "orchestrator",
            AlgorithmSpec::LeapFrogVelocity => "leap_frog_velocity",
            AlgorithmSpec::Thermostat { .. } => "thermostat",
            AlgorithmSpec::LeapFrogPosition => "leap_frog_position",
            AlgorithmSpec::Shake { .. } => "shake",
            AlgorithmSpec::Settle => "settle",
            AlgorithmSpec::Lincs { .. } => "lincs",
            AlgorithmSpec::SteepestDescent { .. } => "steepest_descent",
            AlgorithmSpec::TemperatureCalculation => "temperature_calculation",
            AlgorithmSpec::PressureCalculation { .. } => "pressure_calculation",
            AlgorithmSpec::Barostat { .. } => "barostat",
            AlgorithmSpec::EnergyCalculation => "energy_calculation",
        }
    }

    /// One example value per variant — the registry Python introspects (`gromos.algorithms()`)
    /// and the test that proves `KINDS` is complete.
    pub fn examples() -> Vec<AlgorithmSpec> {
        vec![
            AlgorithmSpec::RemoveCom {
                initial: 1,
                every: 10,
            },
            AlgorithmSpec::Forcefield(ForcefieldPlan {
                cutoff_short: 0.8,
                cutoff_long: 1.4,
                pairlist_every: 5,
                grid_size: 0.0,
                pairlist_algorithm: 0,
                pairlist_type: 0,
                epsilon_rf: 61.0,
                kappa: 0.0,
                four_pi_eps_i: gromos_core::units::four_pi_eps_i,
                bonds: true,
                angles: true,
                impropers: true,
                dihedrals: true,
                parallel: false,
                atoms_per_solvent: 3,
                virial: VirialKind::None,
                position_restraints: None,
                distance_restraints: None,
                perturbation: None,
            }),
            AlgorithmSpec::Orchestrator { terms: Vec::new() },
            AlgorithmSpec::LeapFrogVelocity,
            AlgorithmSpec::Thermostat {
                algorithm: ThermostatAlgorithm::Berendsen,
                temperature: 300.0,
                tau: 0.1,
                dof: 1944.0,
            },
            AlgorithmSpec::LeapFrogPosition,
            AlgorithmSpec::Shake {
                tolerance: 1e-4,
                max_iterations: 1000,
                solute: SoluteConstraints::AllBonds,
                include_solvent: true,
                initial_positions: false,
                initial_velocities: false,
            },
            AlgorithmSpec::Settle,
            AlgorithmSpec::Lincs {
                solute: SoluteConstraints::AllBonds,
                order_solute: 4,
                order_solvent: 4,
                constrain_solute: true,
                constrain_solvent: false,
            },
            AlgorithmSpec::SteepestDescent {
                tolerance: 0.1,
                step0: 0.01,
                step_max: 0.05,
                min_steps: 1,
                force_limit: 0.0,
            },
            AlgorithmSpec::TemperatureCalculation,
            AlgorithmSpec::PressureCalculation {
                virial: VirialKind::Molecular,
            },
            AlgorithmSpec::Barostat {
                pressure0: 1.0,
                compressibility: 4.575e-4,
                tau: 0.5,
            },
            AlgorithmSpec::EnergyCalculation,
        ]
    }

    /// Ordering invariants per kind (GROMOS `md.cc` step order).
    pub fn rules(kind: &str) -> KindRules {
        match kind {
            "remove_com" => KindRules {
                unique: true,
                first: true,
                excludes: &["steepest_descent"],
                ..NONE
            },
            "forcefield" => KindRules {
                unique: true,
                required: true,
                ..NONE
            },
            "orchestrator" => KindRules {
                unique: true,
                must_follow: Some("forcefield"),
                ..NONE
            },
            "leap_frog_velocity" => KindRules {
                unique: true,
                after: &["forcefield"],
                before: &["leap_frog_position"],
                requires: &["leap_frog_position", "temperature_calculation"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "thermostat" => KindRules {
                unique: true,
                after: &["leap_frog_velocity"],
                before: &["leap_frog_position"],
                requires: &["leap_frog_velocity"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "leap_frog_position" => KindRules {
                unique: true,
                after: &["leap_frog_velocity"],
                requires: &["leap_frog_velocity"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "shake" | "settle" | "lincs" => KindRules {
                unique: true,
                after: &["forcefield"],
                before: &["temperature_calculation", "energy_calculation"],
                ..NONE
            },
            "steepest_descent" => KindRules {
                unique: true,
                after: &["forcefield"],
                excludes: DYNAMICS_ONLY,
                ..NONE
            },
            "temperature_calculation" => KindRules {
                unique: true,
                after: &["leap_frog_position"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "pressure_calculation" => KindRules {
                unique: true,
                after: &["temperature_calculation"],
                before: &["barostat"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "barostat" => KindRules {
                unique: true,
                after: &["pressure_calculation"],
                requires: &["pressure_calculation"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "energy_calculation" => KindRules {
                unique: true,
                required: true,
                last: true,
                ..NONE
            },
            _ => NONE,
        }
    }
}

impl TermSpec {
    pub const KINDS: &'static [&'static str] = &["schnet"];

    pub fn name(&self) -> &'static str {
        match self {
            TermSpec::Schnet { .. } => "schnet",
        }
    }

    pub fn examples() -> Vec<TermSpec> {
        vec![TermSpec::Schnet {
            model: "model.pt".into(),
            cutoff: 0.5,
            elements: vec![8, 1, 1],
            region: "1:res(SOL:a)".into(),
            buffer: None,
            coupling: Coupling::Delta,
        }]
    }

    /// The cargo feature this term needs, if any.
    pub fn feature(&self) -> Option<&'static str> {
        match self {
            TermSpec::Schnet { .. } => Some("ml"),
        }
    }

    /// Whether this term's provider reports a virial (a barostat needs one).
    pub fn provides_virial(&self) -> bool {
        match self {
            TermSpec::Schnet { .. } => false,
        }
    }

    pub fn coupling(&self) -> Coupling {
        match self {
            TermSpec::Schnet { coupling, .. } => *coupling,
        }
    }
}

/// Stage 1: the recipe as an ordered, fully resolved list of algorithms.
///
/// Needs the topology (constraint counts → degrees of freedom, solvent shape) and the
/// auxiliary input paths (restraint files become plan values).
pub fn build_plan(
    recipe: &RunRecipe,
    topo: &Topology,
    inputs: &RunInputs,
    four_pi_eps_i: f64,
) -> Result<Vec<AlgorithmSpec>, RunError> {
    let n_atoms = topo.num_atoms();
    let has_solvent = n_atoms > topo.num_solute_atoms();
    let c = &recipe.constraints;
    let sel = ConstraintSelection::from_codes(
        c.solute.ntc(),
        ConstraintAlgorithm::code(c.solute_algorithm),
        ConstraintAlgorithm::code(c.solvent_algorithm),
        has_solvent,
    );
    let ntc_mode = ConstraintSelection::ntc_mode_of(c.solute.ntc());
    let dof = total_dof(topo, &sel, ntc_mode, recipe.boundary.ndfmin);

    let ff = &recipe.forcefield;
    let position_restraints = if ff.restraints.position.mode > 0 {
        let spec = inputs
            .posresspec
            .clone()
            .ok_or_else(|| RunError::MissingInput {
                flag: "posresspec",
                reason: format!("NTPOR={}", ff.restraints.position.mode),
            })?;
        let reference = if ff.restraints.position.reference >= 1 {
            Some(
                inputs
                    .refpos
                    .clone()
                    .ok_or_else(|| RunError::MissingInput {
                        flag: "refpos",
                        reason: format!("NTPORB={}", ff.restraints.position.reference),
                    })?,
            )
        } else {
            None
        };
        Some(PositionRestraintsPlan {
            spec,
            reference,
            k: ff.restraints.position.k,
        })
    } else {
        None
    };
    let distance_restraints = if ff.restraints.distance.mode != 0 {
        let d = &ff.restraints.distance;
        Some(DistanceRestraintsPlan {
            file: inputs
                .distrest
                .clone()
                .ok_or_else(|| RunError::MissingInput {
                    flag: "distrest",
                    reason: format!("NTDIR={}", d.mode),
                })?,
            mode: d.mode.abs(),
            k: d.k,
            r_linear: d.r_linear,
            time_averaging: (d.mode < 0).then_some(d.tau),
            dt: recipe.control.dt,
            lambda: recipe.perturbation.lambda,
        })
    } else {
        None
    };
    let perturbation = (recipe.perturbation.enabled != 0).then_some(PerturbationPlan {
        lambda: recipe.perturbation.lambda,
        nlam: recipe.perturbation.nlam,
        alpha_lj: recipe.perturbation.alpha_lj,
        alpha_c: recipe.perturbation.alpha_c,
    });

    let is_minimization = recipe.is_minimization();
    let barostat = recipe.ensemble.barostat.as_ref();
    let forcefield = AlgorithmSpec::Forcefield(ForcefieldPlan {
        cutoff_short: ff.pairlist.cutoff_short,
        cutoff_long: ff.pairlist.cutoff_long,
        pairlist_every: ff.pairlist.every,
        grid_size: ff.pairlist.grid_size,
        pairlist_algorithm: recipe.to_imd().algorithm,
        pairlist_type: recipe.to_imd().type_,
        epsilon_rf: ff.electrostatics.epsilon_rf,
        kappa: ff.electrostatics.kappa,
        four_pi_eps_i,
        bonds: ff.bonds,
        angles: ff.angles,
        impropers: ff.impropers,
        dihedrals: ff.dihedrals,
        parallel: recipe.execution.parallel.resolve(n_atoms),
        atoms_per_solvent: if topo.solvent_atom_template.is_empty() {
            3
        } else {
            topo.solvent_atom_template.len()
        },
        // The binary only configures the virial for dynamics with pressure coupling.
        virial: match (is_minimization, barostat) {
            (false, Some(b)) => b.virial,
            _ => VirialKind::None,
        },
        position_restraints,
        distance_restraints,
        perturbation,
    });

    let mut plan = Vec::new();
    let constraints = |plan: &mut Vec<AlgorithmSpec>| {
        if sel.shake_enabled() {
            plan.push(AlgorithmSpec::Shake {
                tolerance: c.shake_tolerance,
                max_iterations: 1000,
                solute: if sel.solute_shake {
                    c.solute
                } else {
                    SoluteConstraints::None
                },
                include_solvent: sel.solvent_shake,
                initial_positions: recipe.control.initial_shake >= 1,
                initial_velocities: recipe.control.initial_shake >= 2,
            });
        }
        if sel.settle_enabled {
            plan.push(AlgorithmSpec::Settle);
        }
        if sel.lincs_enabled {
            plan.push(AlgorithmSpec::Lincs {
                solute: c.solute,
                order_solute: c.lincs_order_solute,
                order_solvent: c.lincs_order_solvent,
                constrain_solute: sel.solute_lincs,
                constrain_solvent: sel.solvent_lincs,
            });
        }
    };

    if is_minimization {
        plan.push(forcefield);
        if !ff.terms.is_empty() {
            plan.push(AlgorithmSpec::Orchestrator {
                terms: ff.terms.clone(),
            });
        }
        let m = &recipe.minimisation;
        plan.push(AlgorithmSpec::SteepestDescent {
            tolerance: m.tolerance,
            step0: m.step0,
            step_max: m.step_max,
            min_steps: m.min_steps,
            force_limit: m.force_limit,
        });
        constraints(&mut plan);
        plan.push(AlgorithmSpec::EnergyCalculation);
    } else {
        let com = &recipe.control.com_motion;
        if com.initial >= 1 || com.every != 0 {
            plan.push(AlgorithmSpec::RemoveCom {
                initial: com.initial,
                every: com.every,
            });
        }
        plan.push(forcefield);
        if !ff.terms.is_empty() {
            plan.push(AlgorithmSpec::Orchestrator {
                terms: ff.terms.clone(),
            });
        }
        plan.push(AlgorithmSpec::LeapFrogVelocity);
        if let Some(t) = &recipe.ensemble.thermostat {
            if t.baths.len() > 1 {
                return Err(RunError::Recipe(format!(
                    "{} temperature baths: gromos-rs supports one bath (MULTIBATH NBATHS=1)",
                    t.baths.len()
                )));
            }
            if let Some(bath) = t.baths.first() {
                if bath.tau > 0.0 {
                    plan.push(AlgorithmSpec::Thermostat {
                        algorithm: t.algorithm,
                        temperature: bath.temperature,
                        tau: bath.tau,
                        dof,
                    });
                }
            }
        }
        plan.push(AlgorithmSpec::LeapFrogPosition);
        constraints(&mut plan);
        plan.push(AlgorithmSpec::TemperatureCalculation);
        if let Some(b) = barostat {
            plan.push(AlgorithmSpec::PressureCalculation { virial: b.virial });
            plan.push(AlgorithmSpec::Barostat {
                pressure0: b.pressure0[0][0],
                compressibility: b.compressibility[0][0],
                tau: b.tau,
            });
        }
        plan.push(AlgorithmSpec::EnergyCalculation);
    }
    Ok(plan)
}

/// Check the GROMOS ordering invariants of a plan (edited or not).
pub fn validate_plan(plan: &[AlgorithmSpec]) -> Result<(), RunError> {
    let names: Vec<&str> = plan.iter().map(|a| a.name()).collect();
    let err = |msg: String| Err(RunError::InvalidPlan(msg));
    let position = |kind: &str| names.iter().position(|n| *n == kind);
    let count = |kind: &str| names.iter().filter(|n| **n == kind).count();

    if plan.is_empty() {
        return err("the plan is empty".into());
    }
    for kind in AlgorithmSpec::KINDS {
        let rules = AlgorithmSpec::rules(kind);
        let n = count(kind);
        if rules.required && n != 1 {
            return err(format!("{kind} must appear exactly once (found {n})"));
        }
        if rules.unique && n > 1 {
            return err(format!("{kind} may appear only once (found {n})"));
        }
        let Some(i) = position(kind) else { continue };
        if rules.first && i != 0 {
            return err(format!("{kind} must be the first algorithm"));
        }
        if rules.last && i != names.len() - 1 {
            return err(format!("{kind} must be the last algorithm"));
        }
        if let Some(prev) = rules.must_follow {
            if i == 0 || names[i - 1] != prev {
                return err(format!("{kind} must come immediately after {prev}"));
            }
        }
        for other in rules.after {
            if let Some(j) = position(other) {
                if j > i {
                    return err(format!("{kind} must come after {other}"));
                }
            }
        }
        for other in rules.before {
            if let Some(j) = position(other) {
                if j < i {
                    return err(format!("{kind} must come before {other}"));
                }
            }
        }
        for other in rules.excludes {
            if position(other).is_some() {
                return err(format!("{kind} cannot be combined with {other}"));
            }
        }
        for other in rules.requires {
            if position(other).is_none() {
                return err(format!("{kind} requires {other}"));
            }
        }
    }
    // Terms: a barostat needs every provider to report a virial; `Replace` coupling needs a
    // zone-aware Forcefield (PLAN.md 2.8) — declared, never silently additive.
    for spec in plan {
        if let AlgorithmSpec::Orchestrator { terms } = spec {
            for term in terms {
                if term.coupling() == Coupling::Replace {
                    return err(format!(
                        "term {} with coupling=replace requires the zone-aware Forcefield \
                         (PLAN.md 2.8); only coupling=delta is available",
                        term.name()
                    ));
                }
                if position("barostat").is_some() && !term.provides_virial() {
                    return err(format!(
                        "term {} reports no virial and cannot be combined with a barostat",
                        term.name()
                    ));
                }
            }
        }
    }
    Ok(())
}

/// A plan as JSON (the machine form; `md @dump`).
pub fn plan_to_json(plan: &[AlgorithmSpec]) -> Result<String, RunError> {
    serde_json::to_string_pretty(plan).map_err(|e| RunError::Serde(e.to_string()))
}

pub fn plan_from_json(text: &str) -> Result<Vec<AlgorithmSpec>, RunError> {
    serde_json::from_str(text).map_err(|e| RunError::Serde(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_is_complete() {
        let names: Vec<&str> = AlgorithmSpec::examples().iter().map(|a| a.name()).collect();
        assert_eq!(names, AlgorithmSpec::KINDS, "examples() and KINDS disagree");
        let term_names: Vec<&str> = TermSpec::examples().iter().map(|t| t.name()).collect();
        assert_eq!(term_names, TermSpec::KINDS);
    }

    #[test]
    fn examples_round_trip_through_json() {
        for spec in AlgorithmSpec::examples() {
            let text = serde_json::to_string(&spec).unwrap();
            let back: AlgorithmSpec = serde_json::from_str(&text).unwrap();
            assert_eq!(back, spec);
        }
        for term in TermSpec::examples() {
            let text = serde_json::to_string(&term).unwrap();
            let back: TermSpec = serde_json::from_str(&text).unwrap();
            assert_eq!(back, term);
        }
    }

    #[test]
    fn unknown_term_parameter_is_an_error() {
        let text = r#"{"kind":"schnet","model":"m.pt","cutof":0.5,"elements":[8],"region":"a"}"#;
        assert!(serde_json::from_str::<TermSpec>(text).is_err());
    }

    fn dynamics_plan() -> Vec<AlgorithmSpec> {
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        vec![
            pick("remove_com"),
            pick("forcefield"),
            pick("leap_frog_velocity"),
            pick("thermostat"),
            pick("leap_frog_position"),
            pick("shake"),
            pick("temperature_calculation"),
            pick("pressure_calculation"),
            pick("barostat"),
            pick("energy_calculation"),
        ]
    }

    #[test]
    fn canonical_plans_validate() {
        validate_plan(&dynamics_plan()).unwrap();
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        validate_plan(&[
            pick("forcefield"),
            pick("steepest_descent"),
            pick("shake"),
            pick("energy_calculation"),
        ])
        .unwrap();
    }

    #[test]
    fn violations_are_rejected_with_a_named_reason() {
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        let mut p = dynamics_plan();
        p.remove(9); // no energy_calculation
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("energy_calculation"));

        let mut p = dynamics_plan();
        p.swap(7, 8); // barostat before pressure_calculation
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("pressure_calculation"));

        let mut p = dynamics_plan();
        p.insert(3, pick("orchestrator")); // orchestrator not right after forcefield
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("immediately after"));

        let mut p = dynamics_plan();
        p.insert(2, pick("steepest_descent")); // SD in a dynamics plan
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("cannot be combined"));

        let mut p = dynamics_plan();
        p.push(pick("forcefield")); // second forcefield, and energy_calculation no longer last
        assert!(validate_plan(&p).is_err());
    }
}
