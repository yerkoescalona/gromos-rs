//! Forcefield algorithm — wraps the complete force calculation.
//!
//! Equivalent to gromosXX `algorithm::Forcefield`:
//! 1. Update pairlist (if needed)
//! 2. Calculate bonded forces (bonds, angles, dihedrals)
//! 3. Calculate nonbonded forces (LJ + CRF)
//! 4. Assemble total forces and energies onto Configuration

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::math::Periodicity;
use gromos_core::pairlist::{PairlistContainer, StandardPairlistAlgorithm};
use gromos_core::topology::Topology;

use gromos_forces::bonded::calculate_bonded_forces_ntf;
use gromos_forces::nonbonded::{
    lj_crf_innerloop, rf_excluded_interactions, solvent_innerloop, CRFParameters, ForceStorage,
    LJParameters,
};

/// Forcefield algorithm that performs the complete force calculation each step.
///
/// Owns all the nonbonded interaction parameters and pairlist state.
/// After `apply()`, `conf.current()` has updated forces and potential energies.
pub struct Forcefield {
    /// Lennard-Jones parameter matrix [type_i][type_j]
    pub lj_params: Vec<Vec<LJParameters>>,
    /// Coulomb reaction field parameters
    pub crf_params: CRFParameters,
    /// Boundary conditions
    pub periodicity: Periodicity,
    /// Pairlist container
    pub pairlist: PairlistContainer,
    /// Pairlist construction algorithm
    pub pairlist_algorithm: StandardPairlistAlgorithm,
    /// Whether to use quartic bond potentials (vs harmonic)
    pub use_quartic_bonds: bool,
    /// Whether to run nonbonded in parallel
    pub parallel_nonbonded: bool,
    /// NTF flags: which bonded terms to include (gromosXX FORCE block)
    pub ntf_bond: bool,
    pub ntf_angle: bool,
    pub ntf_improper: bool,
    pub ntf_dihedral: bool,
    /// Number of atoms per solvent molecule (e.g. 3 for water)
    pub atoms_per_solvent: usize,
    /// Reusable nonbonded force storage (avoids allocation per step)
    nonbonded_storage: ForceStorage,
    /// Cached solute pairlist in (u32, u32) format
    pairlist_solute_u32: Vec<(u32, u32)>,
    /// Cached solvent pairlist in (u32, u32) format
    pairlist_solvent_u32: Vec<(u32, u32)>,
    /// Cached charge vector
    charges: Vec<f64>,
    /// Cached IAC vector
    iac_u32: Vec<u32>,
}

impl Forcefield {
    pub fn new(
        lj_params: Vec<Vec<LJParameters>>,
        crf_params: CRFParameters,
        periodicity: Periodicity,
        pairlist: PairlistContainer,
        pairlist_algorithm: StandardPairlistAlgorithm,
    ) -> Self {
        Self {
            lj_params,
            crf_params,
            periodicity,
            pairlist,
            pairlist_algorithm,
            use_quartic_bonds: true,
            parallel_nonbonded: false,
            ntf_bond: true,
            ntf_angle: true,
            ntf_improper: true,
            ntf_dihedral: true,
            atoms_per_solvent: 3,
            nonbonded_storage: ForceStorage::new(0),
            pairlist_solute_u32: Vec::new(),
            pairlist_solvent_u32: Vec::new(),
            charges: Vec::new(),
            iac_u32: Vec::new(),
        }
    }

    /// Convert topology LJ parameters to the nonbonded-crate format.
    pub fn convert_lj_parameters(topo: &Topology) -> Vec<Vec<LJParameters>> {
        topo.lj_parameters
            .iter()
            .map(|row| row.iter().map(Into::into).collect())
            .collect()
    }
}

impl Algorithm for Forcefield {
    fn init(
        &mut self,
        topo: &Topology,
        _conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        let n = topo.inverse_mass.len();
        self.nonbonded_storage = ForceStorage::new(n);
        self.charges = topo.charge.clone();
        self.iac_u32 = topo.iac.iter().map(|&i| i as u32).collect();
        Ok(())
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        let n_atoms = topo.inverse_mass.len();

        // Lazy init if init() was not called
        if self.nonbonded_storage.forces.len() != n_atoms {
            self.nonbonded_storage = ForceStorage::new(n_atoms);
            self.charges = topo.charge.clone();
            self.iac_u32 = topo.iac.iter().map(|&i| i as u32).collect();
        }

        // --- 1. Update pairlist if needed ---
        if self.pairlist.needs_update() {
            self.pairlist_algorithm
                .update(topo, conf, &mut self.pairlist, &self.periodicity);
        }
        self.pairlist.step();

        // Cache solute pairlist as (u32, u32) — solute_short + solute_long
        self.pairlist_solute_u32.clear();
        self.pairlist_solute_u32.extend(
            self.pairlist
                .solute_short
                .iter()
                .chain(self.pairlist.solute_long.iter())
                .map(|&(i, j)| (i as u32, j as u32)),
        );

        // Cache solvent pairlist as (u32, u32) — solvent_short + solvent_long
        self.pairlist_solvent_u32.clear();
        self.pairlist_solvent_u32.extend(
            self.pairlist
                .solvent_short
                .iter()
                .chain(self.pairlist.solvent_long.iter())
                .map(|&(i, j)| (i as u32, j as u32)),
        );

        // --- 2. Calculate bonded forces (conditional on NTF flags) ---
        let any_bonded = self.ntf_bond || self.ntf_angle || self.ntf_improper || self.ntf_dihedral;
        let bonded_result = if any_bonded {
            calculate_bonded_forces_ntf(
                topo, conf, self.use_quartic_bonds,
                self.ntf_bond, self.ntf_angle, self.ntf_dihedral, self.ntf_improper,
            )
        } else {
            gromos_forces::bonded::ForceEnergy::new(n_atoms)
        };

        // --- 3. Calculate nonbonded forces ---
        self.nonbonded_storage.clear();

        // Solute-solute and solute-solvent: lj_crf_innerloop (with HEAVISIDE truncation)
        lj_crf_innerloop(
            &conf.current().pos,
            &self.charges,
            &self.iac_u32,
            &self.pairlist_solute_u32,
            &self.lj_params,
            &self.crf_params,
            &self.periodicity,
            &mut self.nonbonded_storage,
        );

        // Solvent-solvent: solvent_innerloop (shared PBC shift, no HEAVISIDE)
        if !self.pairlist_solvent_u32.is_empty() {
            solvent_innerloop(
                &conf.current().pos,
                &self.charges,
                &self.iac_u32,
                &self.pairlist_solvent_u32,
                &self.lj_params,
                &self.crf_params,
                &self.periodicity,
                self.atoms_per_solvent,
                &mut self.nonbonded_storage,
            );
        }

        // RF self-energy and excluded-pair Coulomb (energy + forces)
        rf_excluded_interactions(
            &self.charges,
            &topo.exclusions,
            &conf.current().pos,
            &self.crf_params,
            &self.periodicity,
            &mut self.nonbonded_storage,
            topo.num_solute_atoms(),
        );

        // --- 4. Assemble forces and energies ---
        let state = conf.current_mut();
        let has_bonded = !bonded_result.forces.is_empty();
        for i in 0..n_atoms {
            let bonded_f = if has_bonded { bonded_result.forces[i] } else { gromos_core::math::Vec3::ZERO };
            state.force[i] = bonded_f + self.nonbonded_storage.forces[i];
        }

        state.energies.bond_total = bonded_result.energy;
        state.energies.lj_total = self.nonbonded_storage.e_lj;
        state.energies.crf_total = self.nonbonded_storage.e_crf;
        state.energies.update_potential_total();

        // Debug: show force magnitudes
        let f_max = state.force.iter().map(|f| f.length()).fold(0.0_f64, f64::max);
        log::debug!("  Bond: {:.10e}  LJ: {:.10e}  CRF: {:.10e}", bonded_result.energy, self.nonbonded_storage.e_lj, self.nonbonded_storage.e_crf);
        log::debug!("  Max |force|: {:.10e}, solute_pairs: {}, solvent_pairs: {}", f_max, self.pairlist_solute_u32.len(), self.pairlist_solvent_u32.len());

        Ok(())
    }

    fn name(&self) -> &str {
        "Forcefield"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::math::{Vacuum, Vec3};
    use gromos_forces::nonbonded::CRFParameters;

    /// Test that the forcefield algorithm computes correct LJ forces for the pair_lj system.
    #[test]
    fn test_forcefield_pair_lj() {
        // Build topology matching pair_lj reference:
        // 2 argon atoms, sigma=0.3405nm, epsilon=0.9961 kJ/mol
        // C6=6.2647e-3, C12=9.847e-6
        let mut topo = Topology::new();
        topo.mass = vec![39.948, 39.948];
        topo.inverse_mass = vec![1.0 / 39.948, 1.0 / 39.948];
        topo.charge = vec![0.0, 0.0];
        topo.iac = vec![0, 0];

        // Exclusions: neither atom excludes the other
        use std::collections::HashSet;
        topo.exclusions = vec![HashSet::new(), HashSet::new()];

        // LJ parameter matrix (1x1 for single atom type)
        let lj = LJParameters { c6: 6.2647e-3, c12: 9.847e-6 };
        let lj_params = vec![vec![lj]];

        // CRF params (epsilon=1, rf_epsilon=1 => pure Coulomb cutoff, but charges=0)
        let crf_params = CRFParameters::new(1.4, 1.0, 1.0, 0.0);

        let periodicity = Periodicity::Vacuum(Vacuum);

        // Pairlist: manually set to contain the single pair
        let mut pairlist = PairlistContainer::new(1.4, 1.4, 0.0);
        pairlist.solute_short = vec![(0, 1)];
        pairlist.update_frequency = 1000; // Don't auto-update in test

        let pairlist_algorithm = StandardPairlistAlgorithm::new(false);

        let mut ff = Forcefield::new(
            lj_params,
            crf_params,
            periodicity,
            pairlist,
            pairlist_algorithm,
        );

        // Configuration: 2 atoms at r=0.35 nm
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.35, 0.0, 0.0),
        ];
        conf.current_mut().vel = vec![Vec3::ZERO, Vec3::ZERO];
        conf.current_mut().force = vec![Vec3::ZERO, Vec3::ZERO];

        let sim = SimulationState::new(0.002, 10);
        ff.init(&topo, &mut conf, &sim).unwrap();
        ff.apply(&topo, &mut conf, &sim).unwrap();

        // Check forces: should be attractive/repulsive along x
        let f0 = conf.current().force[0];
        let f1 = conf.current().force[1];

        // Forces should be equal and opposite (Newton's 3rd law)
        assert!(
            (f0.x + f1.x).abs() < 1e-10,
            "Forces not equal and opposite: f0.x={}, f1.x={}",
            f0.x, f1.x
        );
        assert!(f0.y.abs() < 1e-10);
        assert!(f0.z.abs() < 1e-10);

        // At r=0.35nm for Ar-Ar: net force should be repulsive (r < sigma=0.3405*2^(1/6)=0.382)
        // Actually at r=0.35 which is just below the LJ minimum, let's check:
        // F = -(12*C12/r^13 - 6*C6/r^7)
        // For atom 0: F_x should be negative (pushed towards -x)
        let r = 0.35_f64;
        let r6 = r.powi(6);
        let r7 = r6 * r;
        let r12 = r6 * r6;
        let r13 = r12 * r;
        let c6 = 6.2647e-3;
        let c12 = 9.847e-6;
        let expected_fx = -(12.0 * c12 / r13 - 6.0 * c6 / r7);

        assert!(
            (f0.x - expected_fx).abs() < 1e-6,
            "f0.x = {:.10}, expected = {:.10}",
            f0.x,
            expected_fx
        );

        // Check potential energy
        let expected_e_lj = c12 / r12 - c6 / r6;
        let e_lj = conf.current().energies.lj_total;
        assert!(
            (e_lj - expected_e_lj).abs() < 1e-10,
            "E_LJ = {:.10e}, expected = {:.10e}",
            e_lj,
            expected_e_lj
        );
    }
}
