//! Kinetic degrees of freedom — one formula for the thermostat and the reported temperature.
//!
//! gromosXX (`Multibath::calculate_degrees_of_freedom`): `3N − N_constraints − NDFMIN`, where
//! `N_constraints` counts the solvent constraints of every solvent molecule (when the solvent
//! is constrained) and the solute distance constraints (when the solute is constrained).
//!
//! Before this crate, three copies of this computation existed and none counted the solute
//! constraints (`md.rs` carried a `TODO`); they also subtracted the solvent constraints whenever
//! *any* constraint algorithm was on, even one acting on the solute only. PLAN.md 3.9 A11.

use gromos_core::Topology;
use gromos_integrators::constraints::{NtcMode, ShakeBuffers};

use crate::ConstraintSelection;

/// Total kinetic degrees of freedom of the system.
pub fn total_dof(topo: &Topology, sel: &ConstraintSelection, ntc: NtcMode, ndfmin: i32) -> f64 {
    let n_atoms = topo.num_atoms();
    let solvent_constraints = if sel.solvent_constrained() {
        topo.num_solvent_molecules() * topo.solvent_constraint_template.len()
    } else {
        0
    };
    let solute_constraints = if sel.solute_constrained() {
        ShakeBuffers::new(topo, ntc, false).solute_constraints.len()
    } else {
        0
    };
    (3 * n_atoms - solvent_constraints - solute_constraints) as f64 - ndfmin as f64
}
