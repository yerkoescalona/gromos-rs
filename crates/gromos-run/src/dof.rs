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

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::topology::{
        Bond, BondParameters, MolTypeAtom, SolventAtomTemplate, SolventConstraintTemplate,
    };

    /// Two bonded solute atoms plus one rigid SPC water: 3·5 = 15 degrees of freedom.
    fn dipeptide_and_a_water() -> Topology {
        let mut topo = Topology::new();
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0; 2];
        topo.moltypes[0].atoms = vec![
            MolTypeAtom {
                name: "C1".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
            MolTypeAtom {
                name: "C2".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
        ];
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 0.0,
            r0: 0.15,
        });
        topo.solvent_atom_template = vec![
            SolventAtomTemplate {
                iac: 0,
                name: "OW".into(),
                mass: 15.9994,
                charge: -0.82,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW1".into(),
                mass: 1.008,
                charge: 0.41,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW2".into(),
                mass: 1.008,
                charge: 0.41,
            },
        ];
        topo.solvent_constraint_template = vec![
            SolventConstraintTemplate {
                i: 0,
                j: 1,
                length: 0.1,
            },
            SolventConstraintTemplate {
                i: 0,
                j: 2,
                length: 0.1,
            },
            SolventConstraintTemplate {
                i: 1,
                j: 2,
                length: 0.1633,
            },
        ];
        topo.solvate(1);
        topo
    }

    #[test]
    fn everything_constrained() {
        let topo = dipeptide_and_a_water();
        assert_eq!(topo.num_atoms(), 5);
        let sel = ConstraintSelection::from_codes(3, 1, 1, true);
        // 15 − 3 (water) − 1 (solute bond) − 3 (NDFMIN)
        assert_eq!(total_dof(&topo, &sel, NtcMode::AllBonds, 3), 8.0);
    }

    #[test]
    fn solute_constraints_are_not_subtracted_when_only_the_solvent_is_constrained() {
        let topo = dipeptide_and_a_water();
        let sel = ConstraintSelection::from_codes(1, 1, 1, true);
        assert_eq!(total_dof(&topo, &sel, NtcMode::SolventOnly, 0), 12.0);
    }

    #[test]
    fn nothing_constrained() {
        let topo = dipeptide_and_a_water();
        let sel = ConstraintSelection::from_codes(1, 1, 0, true);
        assert_eq!(total_dof(&topo, &sel, NtcMode::SolventOnly, 0), 15.0);
    }
}
