//! Forcefield algorithm — wraps the complete force calculation.
//!
//! Equivalent to GROMOS `algorithm::Forcefield`:
//! 1. Update pairlist (if needed)
//! 2. Calculate bonded forces (bonds, angles, dihedrals)
//! 3. Calculate nonbonded forces (LJ + CRF)
//! 4. Assemble total forces and energies onto Configuration

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::math::{Mat3, Periodicity, Vec3};
use gromos_core::pairlist::{PairlistAlgorithm, PairlistContainer};
use gromos_core::topology::Topology;

use gromos_forces::bonded::{calculate_bonded_forces_ntf, calculate_perturbed_bonded_forces};
use gromos_forces::nonbonded::{
    lj_crf_innerloop, lj_crf_innerloop_novirial,
    lj_crf_innerloop_parallel, lj_crf_innerloop_parallel_novirial,
    lj_crf_innerloop_cg_grouped, lj_crf_innerloop_cg_grouped_novirial,
    lj_crf_innerloop_cg_grouped_parallel, lj_crf_innerloop_cg_grouped_parallel_novirial,
    one_four_interaction_loop, rf_excluded_interactions,
    solvent_innerloop, solvent_innerloop_novirial,
    solvent_innerloop_parallel, solvent_innerloop_parallel_novirial,
    CGPairGroup, CRFParameters, ForceStorage, LJParamMatrix, LJParameters,
    PertAtomInfo, PertNBCorrection, PerturbedLambdaParams,
    perturbed_pairlist_correction, perturbed_self_energy_correction,
    perturbed_excluded_correction, perturbed_one_four_correction,
    perturbed_atom_pair_correction,
};
use gromos_forces::restraints::{
    DistanceRestraints, PerturbedDistanceRestraints, PositionRestraints,
};

use crate::algorithms::pressure_calculation::VirialType;

/// Forcefield algorithm that performs the complete force calculation each step.
///
/// Owns all the nonbonded interaction parameters and pairlist state.
/// After `apply()`, `conf.current()` has updated forces and potential energies.
pub struct Forcefield {
    /// Flat Lennard-Jones parameter matrix (cache-friendly contiguous layout)
    pub lj_params: LJParamMatrix,
    /// Coulomb reaction field parameters
    pub crf_params: CRFParameters,
    /// Boundary conditions
    pub periodicity: Periodicity,
    /// Pairlist container
    pub pairlist: PairlistContainer,
    /// Pairlist construction algorithm
    pub pairlist_algorithm: PairlistAlgorithm,
    /// Whether to use quartic bond potentials (vs harmonic)
    pub use_quartic_bonds: bool,
    /// Whether to run nonbonded in parallel
    pub parallel_nonbonded: bool,
    /// NTF flags: which bonded terms to include (GROMOS FORCE block)
    pub ntf_bond: bool,
    pub ntf_angle: bool,
    pub ntf_improper: bool,
    pub ntf_dihedral: bool,
    /// Number of atoms per solvent molecule (e.g. 3 for water)
    pub atoms_per_solvent: usize,
    /// Virial type: None, Atomic, or Molecular (controls KE tensor + virial correction)
    pub virial_type: VirialType,
    /// Reusable nonbonded force storage (avoids allocation per step)
    nonbonded_storage: ForceStorage,
    /// Cached solute pairlist in (u32, u32) format
    pairlist_solute_u32: Vec<(u32, u32)>,
    /// CG-pair group metadata for the solute shortrange pairlist
    solute_cg_groups: Vec<CGPairGroup>,
    /// Cached solvent pairlist in (u32, u32) format
    pairlist_solvent_u32: Vec<(u32, u32)>,
    /// Cached long-range solute pairlist (u32, u32)
    pairlist_solute_long_u32: Vec<(u32, u32)>,
    /// Cached long-range solvent pairlist (u32, u32)
    pairlist_solvent_long_u32: Vec<(u32, u32)>,
    /// Cached charge vector
    charges: Vec<f64>,
    /// Cached IAC vector
    iac_u32: Vec<u32>,
    /// Cached long-range forces (twin-range: reused between pairlist updates)
    longrange_forces: Vec<Vec3>,
    /// Cached long-range LJ energy
    longrange_e_lj: f64,
    /// Cached long-range CRF energy
    longrange_e_crf: f64,
    /// Cached long-range virial tensor (twin-range: reused between pairlist updates)
    /// Must be added to nonbonded_storage.virial each step, just like forces/energies.
    longrange_virial: [[f64; 3]; 3],
    /// Whether twin-range is active (RCUTP < RCUTL)
    twin_range_active: bool,
    /// Whether long-range forces have been computed at least once
    longrange_computed: bool,
    /// Whether CG-grouped solute innerloop is safe (box large enough relative to cutoff + CG size)
    use_cg_grouped_solute: bool,
    /// Position restraints (empty = none)
    pub position_restraints: PositionRestraints,
    /// Distance restraints (empty = none)
    pub distance_restraints: DistanceRestraints,
    /// Perturbed distance restraints (empty = none)
    pub perturbed_distance_restraints: PerturbedDistanceRestraints,
    /// Current lambda value (for perturbed restraints)
    pub lambda: f64,
    /// (lambda, lambda_derivative) for FEP/TI perturbed bonded terms.
    /// lambda_derivative = NLAM * RLAM^(NLAM-1); set from ImdParameters::lambda_and_derivative().
    pub lambda_and_derivative: (f64, f64),
    /// NLAM exponent (1 = linear, 2 = quadratic, etc.)
    pub lambda_exp: i32,
    /// Global soft-core alpha for LJ (ALPHLJ from imd PERTURBATION block).
    /// Per-atom lj_soft from .ptp is multiplied by this to get effective alpha_lj,
    /// matching GROMOS in_perturbation.cc:1308 `lj_soft *= param.perturbation.soft_vdw`.
    pub global_alphlj: f64,
    /// Global soft-core alpha for CRF (ALPHC from imd PERTURBATION block).
    pub global_alphc: f64,
    /// Per-atom perturbation data; None for non-perturbed atoms.
    /// Populated when a pttopo is loaded (topo.perturbed_solute.atoms is non-empty).
    pert_info: Vec<Option<PertAtomInfo>>,
}

impl Forcefield {
    pub fn new(
        lj_params: Vec<Vec<LJParameters>>,
        crf_params: CRFParameters,
        periodicity: Periodicity,
        pairlist: PairlistContainer,
        pairlist_algorithm: PairlistAlgorithm,
    ) -> Self {
        let twin_range_active = pairlist.short_range_cutoff < pairlist.long_range_cutoff - 1e-10;
        Self {
            lj_params: LJParamMatrix::from_nested(&lj_params),
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
            virial_type: VirialType::None,
            nonbonded_storage: ForceStorage::new(0),
            pairlist_solute_u32: Vec::new(),
            solute_cg_groups: Vec::new(),
            pairlist_solvent_u32: Vec::new(),
            pairlist_solute_long_u32: Vec::new(),
            pairlist_solvent_long_u32: Vec::new(),
            charges: Vec::new(),
            iac_u32: Vec::new(),
            longrange_forces: Vec::new(),
            longrange_e_lj: 0.0,
            longrange_e_crf: 0.0,
            longrange_virial: [[0.0; 3]; 3],
            twin_range_active,
            longrange_computed: false,
            use_cg_grouped_solute: false, // computed at init when box size is known
            position_restraints: PositionRestraints::new(),
            distance_restraints: DistanceRestraints::new(),
            perturbed_distance_restraints: PerturbedDistanceRestraints::new(),
            lambda: 0.0,
            lambda_and_derivative: (0.0, 1.0),
            lambda_exp: 1,
            global_alphlj: 1.0,
            global_alphc: 1.0,
            pert_info: Vec::new(),
        }
    }

    /// Build per-atom perturbation info from the topology's perturbed_solute.
    ///
    /// GROMOS multiplies the per-atom lj_soft/crf_soft from the .ptp file by the
    /// global ALPHLJ/ALPHC from the .imd PERTURBATION block before storing them
    /// (in_perturbation.cc:1308: `lj_soft *= param.perturbation.soft_vdw`).
    /// We replicate that here: effective alpha_lj = pa.lj_soft * global_alphlj.
    fn build_pert_info(topo: &Topology, global_alphlj: f64, global_alphc: f64) -> Vec<Option<PertAtomInfo>> {
        let n = topo.inverse_mass.len();
        let mut info: Vec<Option<PertAtomInfo>> = vec![None; n];
        for pa in &topo.perturbed_solute.atoms {
            if pa.seq < n {
                // topo.iac and PerturbedAtom.a_iac / b_iac are both 0-indexed:
                //   topo.iac   → topology.rs reader subtracts 1 from file's 1-indexed value
                //   pa.a_iac   → ptp.rs reader subtracts 1 from pttopo's 1-indexed value
                // No conversion needed here.
                info[pa.seq] = Some(PertAtomInfo {
                    a_iac:    pa.a_iac,
                    b_iac:    pa.b_iac,
                    a_charge: pa.a_charge,
                    b_charge: pa.b_charge,
                    alpha_lj:  pa.lj_soft  * global_alphlj,
                    alpha_crf: pa.crf_soft * global_alphc,
                });
            }
        }
        info
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
        self.longrange_forces = vec![Vec3::ZERO; n];
        self.pert_info = Self::build_pert_info(topo, self.global_alphlj, self.global_alphc);
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
            self.longrange_forces = vec![Vec3::ZERO; n_atoms];
            self.pert_info = Self::build_pert_info(topo, self.global_alphlj, self.global_alphc);
        }

        // IMPORTANT: Update periodicity from current box each step.
        // Under NPT, the barostat scales the box at the end of the previous step,
        // so self.periodicity (which caches box_size, half_box, inv_box) must be
        // refreshed before any force/pairlist/virial calculation uses it.
        // Without this, all PBC nearest_image calls use stale box dimensions.
        if !matches!(self.periodicity, Periodicity::Vacuum(_)) {
            use gromos_core::configuration::BoxType;
            use gromos_core::math::Rectangular;
            let box_cfg = &conf.current().box_config;
            match box_cfg.box_type {
                BoxType::Rectangular => {
                    let dims = box_cfg.dimensions();
                    if dims.x > 0.0 && dims.y > 0.0 && dims.z > 0.0 {
                        self.periodicity = Periodicity::Rectangular(Rectangular::new(dims));
                    }
                }
                BoxType::Triclinic | BoxType::TruncatedOctahedral => {
                    use gromos_core::math::Triclinic;
                    self.periodicity = Periodicity::Triclinic(Triclinic::new(box_cfg.vectors));
                }
                _ => {}
            }
        }

        // Determine if CG-grouped solute innerloop is safe for this box.
        // Safety condition: half_box_min > cutoff_long + max_cg_diameter
        // When this fails, atom pairs near the half-box boundary can get
        // different periodic image rounding than their CG reference atoms.
        let cg_data_available = topo.atom_to_chargegroup.len() == topo.num_atoms()
            && !topo.chargegroups.is_empty();
        self.use_cg_grouped_solute = cg_data_available && match &self.periodicity {
            Periodicity::Vacuum(_) => true, // no PBC → shift is always 0, safe
            Periodicity::Rectangular(rect) => {
                let half_min = rect.half_box.x.min(rect.half_box.y).min(rect.half_box.z);
                let cutoff = self.pairlist.long_range_cutoff + self.pairlist.skin;
                // Conservative CG diameter estimate (0.3nm covers water + most molecular CGs)
                half_min > cutoff + 0.3
            }
            Periodicity::Triclinic(_) => false, // be conservative for triclinic
        };

        // --- prepare_virial: compute kinetic energy tensor (GROMOS: before force calc) ---
        // Uses current velocities and positions; stored in current().kinetic_energy_tensor
        if self.virial_type != VirialType::None {
            let ke_tensor = if self.virial_type == VirialType::Molecular && !topo.pressure_groups.is_empty() {
                compute_molecular_ke_tensor(topo, conf, &self.periodicity)
            } else {
                compute_atomic_ke_tensor(topo, conf, n_atoms)
            };
            conf.current_mut().kinetic_energy_tensor = ke_tensor;
        }

        // --- 1. Update pairlist if needed ---
        let pairlist_updated = self.pairlist.needs_update();
        if pairlist_updated {
            let t_pl = std::time::Instant::now();
            self.pairlist_algorithm
                .update(topo, conf, &mut self.pairlist, &self.periodicity);
            log::trace!("    pairlist update: {:.3} ms", t_pl.elapsed().as_secs_f64() * 1000.0);
        }
        self.pairlist.step();

        // Cache short-range solute pairlist as (u32, u32) — only on pairlist update or first call
        if pairlist_updated || self.pairlist_solute_u32.is_empty() {
            self.pairlist_solute_u32.clear();
            self.pairlist_solute_u32.extend(
                self.pairlist
                    .solute_short
                    .iter()
                    .map(|&(i, j)| (i as u32, j as u32)),
            );

            // Build CG-pair group metadata only if we'll use the CG-grouped kernel.
            self.solute_cg_groups.clear();
            if self.use_cg_grouped_solute && !self.pairlist_solute_u32.is_empty() {
                let atom_to_cg = &topo.atom_to_chargegroup;
                let mut start = 0u32;
                let (first_i, first_j) = self.pairlist_solute_u32[0];
                let mut cur_cg_i = atom_to_cg[first_i as usize];
                let mut cur_cg_j = atom_to_cg[first_j as usize];
                let mut ref_atom_i = topo.chargegroups[cur_cg_i].atoms[0] as u32;
                let mut ref_atom_j = topo.chargegroups[cur_cg_j].atoms[0] as u32;

                for idx in 1..self.pairlist_solute_u32.len() {
                    let (ii, jj) = self.pairlist_solute_u32[idx];
                    let cg_i = atom_to_cg[ii as usize];
                    let cg_j = atom_to_cg[jj as usize];
                    if cg_i != cur_cg_i || cg_j != cur_cg_j {
                        // Close current group
                        self.solute_cg_groups.push(CGPairGroup {
                            ref_atom_i,
                            ref_atom_j,
                            start,
                            end: idx as u32,
                        });
                        // Start new group
                        start = idx as u32;
                        cur_cg_i = cg_i;
                        cur_cg_j = cg_j;
                        ref_atom_i = topo.chargegroups[cg_i].atoms[0] as u32;
                        ref_atom_j = topo.chargegroups[cg_j].atoms[0] as u32;
                    }
                }
                // Close final group
                self.solute_cg_groups.push(CGPairGroup {
                    ref_atom_i,
                    ref_atom_j,
                    start,
                    end: self.pairlist_solute_u32.len() as u32,
                });
            }

            // Cache short-range solvent pairlist as (u32, u32)
            self.pairlist_solvent_u32.clear();
            self.pairlist_solvent_u32.extend(
                self.pairlist
                    .solvent_short
                    .iter()
                    .map(|&(i, j)| (i as u32, j as u32)),
            );

            // Cache long-range pairlists
            self.pairlist_solute_long_u32.clear();
            self.pairlist_solute_long_u32.extend(
                self.pairlist
                    .solute_long
                    .iter()
                    .map(|&(i, j)| (i as u32, j as u32)),
            );
            self.pairlist_solvent_long_u32.clear();
            self.pairlist_solvent_long_u32.extend(
                self.pairlist
                    .solvent_long
                    .iter()
                    .map(|&(i, j)| (i as u32, j as u32)),
            );
        }

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

        // --- 2b. Perturbed bonded forces (FEP/TI, NTG != 0) ---
        // Forces are accumulated in the returned ForceEnergyLambda and applied in step 4.
        let (lambda, lambda_deriv) = self.lambda_and_derivative;
        let perturbed_bonded = if !topo.perturbed_solute.is_empty() {
            let pr = calculate_perturbed_bonded_forces(topo, conf, lambda, lambda_deriv);
            log::debug!("  Perturbed bonded: E={:.10e} dH/dλ={:.10e}", pr.energy, pr.lambda_derivative);
            Some(pr)
        } else {
            None
        };

        // --- 3. Calculate nonbonded forces ---
        self.nonbonded_storage.clear();
        let need_virial = self.virial_type != VirialType::None;

        // Short-range solute interactions
        let t_solute = std::time::Instant::now();
        if self.use_cg_grouped_solute {
            // CG-grouped: compute nearest_image once per CG pair (fast path for large boxes)
            if self.parallel_nonbonded {
                let result = if need_virial {
                    lj_crf_innerloop_cg_grouped_parallel(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solute_u32,
                        &self.solute_cg_groups,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        n_atoms,
                    )
                } else {
                    lj_crf_innerloop_cg_grouped_parallel_novirial(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solute_u32,
                        &self.solute_cg_groups,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        n_atoms,
                    )
                };
                self.nonbonded_storage.merge(&result);
            } else if need_virial {
                lj_crf_innerloop_cg_grouped(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist_solute_u32,
                    &self.solute_cg_groups,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    &mut self.nonbonded_storage,
                );
            } else {
                lj_crf_innerloop_cg_grouped_novirial(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist_solute_u32,
                    &self.solute_cg_groups,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    &mut self.nonbonded_storage,
                );
            }
        } else {
            // Per-atom nearest_image fallback (exact, for small boxes)
            if self.parallel_nonbonded {
                let result = if need_virial {
                    lj_crf_innerloop_parallel(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solute_u32,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        n_atoms,
                    )
                } else {
                    lj_crf_innerloop_parallel_novirial(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solute_u32,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        n_atoms,
                    )
                };
                self.nonbonded_storage.merge(&result);
            } else if need_virial {
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
            } else {
                lj_crf_innerloop_novirial(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist_solute_u32,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    &mut self.nonbonded_storage,
                );
            }
        }
        let e_lj_after_solute = self.nonbonded_storage.e_lj;
        let e_crf_after_solute = self.nonbonded_storage.e_crf;
        log::trace!("    solute forces: {:.3} ms ({} pairs)", t_solute.elapsed().as_secs_f64() * 1000.0, self.pairlist_solute_u32.len());
        log::debug!("  After solute innerloop: e_lj={:.10e}, e_crf={:.10e}", e_lj_after_solute, e_crf_after_solute);

        // Short-range solvent-solvent interactions (shared PBC shift)
        let t_solvent = std::time::Instant::now();
        if !self.pairlist_solvent_u32.is_empty() {
            if self.parallel_nonbonded {
                let result = if need_virial {
                    solvent_innerloop_parallel(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solvent_u32,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.atoms_per_solvent,
                        n_atoms,
                    )
                } else {
                    solvent_innerloop_parallel_novirial(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solvent_u32,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.atoms_per_solvent,
                        n_atoms,
                    )
                };
                self.nonbonded_storage.merge(&result);
            } else if need_virial {
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
            } else {
                solvent_innerloop_novirial(
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
        }
        log::trace!("    solvent forces: {:.3} ms ({} pairs)", t_solvent.elapsed().as_secs_f64() * 1000.0, self.pairlist_solvent_u32.len());
        let e_lj_after_solvent = self.nonbonded_storage.e_lj;
        let e_crf_after_solvent = self.nonbonded_storage.e_crf;
        log::debug!("  After solvent innerloop: e_lj={:.10e}, e_crf={:.10e} (delta_lj={:.10e}, delta_crf={:.10e})",
            e_lj_after_solvent, e_crf_after_solvent,
            e_lj_after_solvent - e_lj_after_solute, e_crf_after_solvent - e_crf_after_solute);

        // Twin-range long-range forces: recalculate on pairlist update, reuse otherwise
        if self.twin_range_active {
            // Recalculate long-range if pairlist was updated OR first step
            if pairlist_updated || !self.longrange_computed {
                // Recalculate long-range forces and cache them
                let mut lr_storage = ForceStorage::new(n_atoms);

                // Long-range solute interactions
                if !self.pairlist_solute_long_u32.is_empty() {
                    lj_crf_innerloop(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solute_long_u32,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        &mut lr_storage,
                    );
                }
                let lr_solute_e_lj = lr_storage.e_lj;
                let lr_solute_e_crf = lr_storage.e_crf;

                // Long-range solvent-solvent interactions
                // Pairlist stores expanded atom pairs (per-atom nearest_image)
                if !self.pairlist_solvent_long_u32.is_empty() {
                    lj_crf_innerloop(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist_solvent_long_u32,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        &mut lr_storage,
                    );
                }
                log::debug!("  LR solute: e_lj={:.10e}, e_crf={:.10e} ({} pairs)",
                    lr_solute_e_lj, lr_solute_e_crf, self.pairlist_solute_long_u32.len());
                log::debug!("  LR solvent: e_lj={:.10e}, e_crf={:.10e} ({} pairs)",
                    lr_storage.e_lj - lr_solute_e_lj, lr_storage.e_crf - lr_solute_e_crf,
                    self.pairlist_solvent_long_u32.len());

                // Cache the long-range results
                self.longrange_forces.clear();
                self.longrange_forces.extend_from_slice(&lr_storage.forces);
                self.longrange_e_lj = lr_storage.e_lj;
                self.longrange_e_crf = lr_storage.e_crf;
                self.longrange_virial = lr_storage.virial;
                self.longrange_computed = true;
                log::debug!("  Long-range (recomputed): lr_e_lj={:.10e}, lr_e_crf={:.10e}, solute_long={}, solvent_long={}",
                    lr_storage.e_lj, lr_storage.e_crf,
                    self.pairlist_solute_long_u32.len(), self.pairlist_solvent_long_u32.len());
            }

            // Add cached long-range forces, energies, and virial to current step.
            // IMPORTANT: The virial from long-range pairs MUST be included here.
            // Without it, the virial trace is ~1000 kJ/mol too small, causing
            // wildly wrong pressures and barostat box-collapse under NPT.
            for i in 0..n_atoms {
                self.nonbonded_storage.forces[i] += self.longrange_forces[i];
            }
            self.nonbonded_storage.e_lj += self.longrange_e_lj;
            self.nonbonded_storage.e_crf += self.longrange_e_crf;
            for a in 0..3 {
                for b in 0..3 {
                    self.nonbonded_storage.virial[a][b] += self.longrange_virial[a][b];
                }
            }
        } else {
            // No twin-range: evaluate long-range pairs as short-range (same as before)
            if !self.pairlist_solute_long_u32.is_empty() {
                lj_crf_innerloop(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist_solute_long_u32,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    &mut self.nonbonded_storage,
                );
            }
            if !self.pairlist_solvent_long_u32.is_empty() {
                solvent_innerloop(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist_solvent_long_u32,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.atoms_per_solvent,
                    &mut self.nonbonded_storage,
                );
            }
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

        // 1-4 interactions: LJ with cs6/cs12 + CRF with coulomb scaling
        // GROMOS: computed separately from pairlist, no cutoff
        if !topo.one_four_pairs.is_empty() {
            one_four_interaction_loop(
                &topo.one_four_pairs,
                &conf.current().pos,
                &self.charges,
                &self.iac_u32,
                &self.lj_params,
                &self.crf_params,
                &self.periodicity,
                &mut self.nonbonded_storage,
                1.0, // coulomb_scaling: 1.0 for standard GROMOS
            );
        }

        // --- 3b. Perturbed nonbonded corrections (FEP/TI) ---
        // The four steps below mirror GROMOS's perturbed_nonbonded_outerloop:
        //   perturbed_pairlist_correction   → replaces state-A contribution for perturbed pairs
        //   perturbed_self_energy_correction → corrects RF self for perturbed atoms
        //   perturbed_excluded_correction    → corrects RF excluded pairs for perturbed atoms
        //   perturbed_one_four_correction    → corrects 1-4 for perturbed atoms
        let pert_nb_dhdl = if !self.pert_info.is_empty()
            && self.pert_info.iter().any(|x| x.is_some())
        {
            let (lam, _) = self.lambda_and_derivative;
            let lp = PerturbedLambdaParams::from_lambda(
                lam, lam, lam, lam, 1.0, 1.0, 1.0, 1.0,
                self.lambda_exp,
            );
            // 3b-1: pairlist correction (delta over what state-A innerloop already computed)
            // Applied to BOTH short-range and long-range pairlists:
            // the regular innerloops computed state-A for all pairs; here we replace the
            // contribution for perturbed pairs with the full soft-core dual-topology result.
            let corr_pl = perturbed_pairlist_correction(
                &conf.current().pos,
                &self.pairlist_solute_u32,
                &self.pert_info,
                &self.charges,
                &self.iac_u32,
                &self.lj_params,
                &self.crf_params,
                &lp,
                &self.periodicity,
            );
            self.nonbonded_storage.e_lj  += corr_pl.delta_e_lj;
            self.nonbonded_storage.e_crf += corr_pl.delta_e_crf;
            let mut dhdl_nb = corr_pl.dhdl;
            for i in 0..n_atoms {
                self.nonbonded_storage.forces[i] += corr_pl.forces[i];
            }

            // 3b-1b: same correction for ALL remaining pairlists that may contain
            // perturbed pairs. Solute-solvent long-range pairs are in solute_long.
            let pls: &[(&str, &[(u32,u32)])] = &[
                ("solute_long",   &self.pairlist_solute_long_u32),
                ("solvent_short", &self.pairlist_solvent_u32),
                ("solvent_long",  &self.pairlist_solvent_long_u32),
            ];
            for (_name, pl) in pls {
                if pl.is_empty() { continue; }
                let corr = perturbed_pairlist_correction(
                    &conf.current().pos,
                    pl,
                    &self.pert_info,
                    &self.charges,
                    &self.iac_u32,
                    &self.lj_params,
                    &self.crf_params,
                    &lp,
                    &self.periodicity,
                );
                self.nonbonded_storage.e_lj  += corr.delta_e_lj;
                self.nonbonded_storage.e_crf += corr.delta_e_crf;
                dhdl_nb += corr.dhdl;
                for i in 0..n_atoms {
                    self.nonbonded_storage.forces[i] += corr.forces[i];
                }
            }
            log::debug!("  Pert corr solute_short: ({} pairs, Δe_lj={:.6e}, Δe_crf={:.6e})",
                self.pairlist_solute_u32.len(), corr_pl.delta_e_lj, corr_pl.delta_e_crf);
            log::debug!("  Perturbed NB: Δe_lj_pl={:.6e} Δe_crf_pl={:.6e} Δe_self=TBD dH/dλ_nb={:.6e}",
                corr_pl.delta_e_lj, corr_pl.delta_e_crf, dhdl_nb);

            // 3b-2: RF self-energy correction for perturbed atoms
            let (de_self, dhdl_self) = perturbed_self_energy_correction(
                &self.pert_info, &self.crf_params, &lp,
            );
            self.nonbonded_storage.e_crf += de_self;
            dhdl_nb += dhdl_self;

            // 3b-3: excluded-pair CRF correction (separate accumulator to avoid double-add)
            {
                use gromos_forces::nonbonded::PertNBCorrection;
                let mut corr_ex = PertNBCorrection::new(n_atoms);
                perturbed_excluded_correction(
                    &topo.exclusions,
                    &conf.current().pos,
                    &self.pert_info,
                    &self.charges,
                    &self.crf_params,
                    &lp,
                    &self.periodicity,
                    topo.num_solute_atoms(),
                    &mut corr_ex,
                );
                self.nonbonded_storage.e_crf += corr_ex.delta_e_crf;
                dhdl_nb += corr_ex.dhdl;
                for i in 0..n_atoms {
                    self.nonbonded_storage.forces[i] += corr_ex.forces[i];
                }
            }

            // 3b-4: 1-4 correction (separate accumulator)
            if !topo.one_four_pairs.is_empty() {
                use gromos_forces::nonbonded::PertNBCorrection;
                let mut corr14 = PertNBCorrection::new(n_atoms);
                perturbed_one_four_correction(
                    &topo.one_four_pairs,
                    &conf.current().pos,
                    &self.pert_info,
                    &self.charges,
                    &self.iac_u32,
                    &self.lj_params,
                    &self.crf_params,
                    &lp,
                    &self.periodicity,
                    &mut corr14,
                );
                log::debug!("  1-4 corr: Δe_lj={:.6e} Δe_crf={:.6e}", corr14.delta_e_lj, corr14.delta_e_crf);
                self.nonbonded_storage.e_lj  += corr14.delta_e_lj;
                self.nonbonded_storage.e_crf += corr14.delta_e_crf;
                dhdl_nb += corr14.dhdl;
                for i in 0..n_atoms {
                    self.nonbonded_storage.forces[i] += corr14.forces[i];
                }
            }

            // 3b-5: PERTATOMPAIR correction (special pair types replacing 1-4 LJ)
            if !topo.perturbed_solute.atom_pairs.is_empty() {
                let mut corr_ap = PertNBCorrection::new(n_atoms);
                perturbed_atom_pair_correction(
                    &topo.perturbed_solute.atom_pairs,
                    &topo.one_four_pairs,
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.lj_params,
                    &self.crf_params,
                    &lp,
                    &self.periodicity,
                    &mut corr_ap,
                );
                self.nonbonded_storage.e_lj += corr_ap.delta_e_lj;
                dhdl_nb += corr_ap.dhdl;
                for i in 0..n_atoms {
                    self.nonbonded_storage.forces[i] += corr_ap.forces[i];
                }
                log::debug!("  AtomPair corr: Δe_lj={:.6e} Δe_crf={:.6e}", corr_ap.delta_e_lj, corr_ap.delta_e_crf);
            }

            log::debug!("  Perturbed NB: Δe_lj_pl={:.6e} Δe_crf_pl={:.6e} Δe_self={:.6e} dH/dλ_nb={:.6e}",
                corr_pl.delta_e_lj, corr_pl.delta_e_crf, de_self, dhdl_nb);
            log::debug!("  nb_storage after all: e_lj={:.6e} e_crf={:.6e}", self.nonbonded_storage.e_lj, self.nonbonded_storage.e_crf);

            dhdl_nb
        } else {
            0.0
        };

        // --- 4. Assemble forces and energies ---
        {
            let state = conf.current_mut();
            let has_bonded = !bonded_result.forces.is_empty();
            for i in 0..n_atoms {
                let bonded_f = if has_bonded { bonded_result.forces[i] } else { gromos_core::math::Vec3::ZERO };
                let pert_f = perturbed_bonded.as_ref()
                    .and_then(|p| p.forces.get(i).copied())
                    .unwrap_or(gromos_core::math::Vec3::ZERO);
                state.force[i] = bonded_f + pert_f + self.nonbonded_storage.forces[i];
            }

            state.energies.bond_total = bonded_result.energy
                + perturbed_bonded.as_ref().map_or(0.0, |p| p.energy);
            state.energies.dhdl_total = perturbed_bonded.as_ref().map_or(0.0, |p| p.lambda_derivative)
                + pert_nb_dhdl;
            state.energies.lj_total = self.nonbonded_storage.e_lj;
            state.energies.crf_total = self.nonbonded_storage.e_crf;
            state.energies.update_potential_total();

            // Transfer virial tensor to configuration (nonbonded + bonded)
            // GROMOS convention: virial[a][b] = Σ r[a] * force[b]
            let vir = &self.nonbonded_storage.virial;
            let bvir = &bonded_result.virial;
            state.virial_tensor = Mat3::from_cols(
                Vec3::new(vir[0][0] + bvir[0][0], vir[1][0] + bvir[1][0], vir[2][0] + bvir[2][0]),
                Vec3::new(vir[0][1] + bvir[0][1], vir[1][1] + bvir[1][1], vir[2][1] + bvir[2][1]),
                Vec3::new(vir[0][2] + bvir[0][2], vir[1][2] + bvir[1][2], vir[2][2] + bvir[2][2]),
            );

            // Debug: show force magnitudes
            let f_max = state.force.iter().map(|f| f.length()).fold(0.0_f64, f64::max);
            log::debug!("  Bond: {:.10e}  LJ: {:.10e}  CRF: {:.10e}", bonded_result.energy, self.nonbonded_storage.e_lj, self.nonbonded_storage.e_crf);
            log::debug!("  Max |force|: {:.10e}, solute_pairs: {}, solvent_pairs: {}", f_max, self.pairlist_solute_u32.len(), self.pairlist_solvent_u32.len());
        }

        // --- 4b. Position restraints (special interaction) ---
        if !self.position_restraints.restraints.is_empty() {
            let e_posres = self.position_restraints.calculate_all(conf, &self.periodicity);
            conf.current_mut().energies.special_total = e_posres;
            log::debug!("  Position restraints: {:.10e} kJ/mol ({} atoms)",
                e_posres, self.position_restraints.restraints.len());
        }

        // --- 4c. Distance restraints ---
        if !self.distance_restraints.restraints.is_empty() {
            let e_dr = self.distance_restraints.calculate_all(conf, &self.periodicity);
            conf.current_mut().energies.distanceres_total = e_dr;
            log::debug!("  Distance restraints: {:.10e} kJ/mol", e_dr);
        }
        if !self.perturbed_distance_restraints.restraints.is_empty() {
            let e_pdr = self.perturbed_distance_restraints.calculate_all(conf, &self.periodicity, self.lambda);
            conf.current_mut().energies.distanceres_total += e_pdr;
            log::debug!("  Perturbed distance restraints: {:.10e} kJ/mol", e_pdr);
        }

        // --- 5. atomic_to_molecular_virial: correct virial from atomic to molecular ---
        // GROMOS: last interaction in forcefield, subtracts intramolecular virial contributions
        if self.virial_type == VirialType::Molecular && !topo.pressure_groups.is_empty() {
            apply_molecular_virial_correction(topo, conf, &self.periodicity);
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "Forcefield"
    }
}

// === Virial helper functions (GROMOS: prepare_virial.cc) ===

/// Compute atomic kinetic energy tensor: KE_ij = 0.5 * Σ_k m_k * v_k_i * v_k_j
/// Uses conf.current().vel (GROMOS convention: computed before force calc)
fn compute_atomic_ke_tensor(topo: &Topology, conf: &Configuration, n_atoms: usize) -> Mat3 {
    let vel = &conf.current().vel;
    let mut ke = [[0.0f64; 3]; 3];
    for k in 0..n_atoms {
        let m = topo.mass[k];
        let v = vel[k];
        let vv = [v.x, v.y, v.z];
        for a in 0..3 {
            for b in 0..3 {
                ke[a][b] += 0.5 * m * vv[a] * vv[b];
            }
        }
    }
    Mat3::from_cols(
        Vec3::new(ke[0][0], ke[1][0], ke[2][0]),
        Vec3::new(ke[0][1], ke[1][1], ke[2][1]),
        Vec3::new(ke[0][2], ke[1][2], ke[2][2]),
    )
}

/// Compute molecular kinetic energy tensor using COM velocities per pressure group.
/// GROMOS: prepare_virial.cc _centre_of_mass with chain-gathering.
/// Uses conf.current().vel and conf.current().pos.
fn compute_molecular_ke_tensor(topo: &Topology, conf: &Configuration, _periodicity: &Periodicity) -> Mat3 {
    let vel = &conf.current().vel;
    let mut ke = [[0.0f64; 3]; 3];

    for pg in &topo.pressure_groups {
        let mut mv = [0.0f64; 3]; // mass-weighted velocity sum
        let mut m_total = 0.0f64;

        for i in pg.clone() {
            let m = topo.mass[i];
            m_total += m;
            mv[0] += m * vel[i].x;
            mv[1] += m * vel[i].y;
            mv[2] += m * vel[i].z;
        }
        if m_total > 0.0 {
            for a in 0..3 {
                for b in 0..3 {
                    ke[a][b] += 0.5 * mv[a] * mv[b] / m_total;
                }
            }
        }
    }
    Mat3::from_cols(
        Vec3::new(ke[0][0], ke[1][0], ke[2][0]),
        Vec3::new(ke[0][1], ke[1][1], ke[2][1]),
        Vec3::new(ke[0][2], ke[1][2], ke[2][2]),
    )
}

/// Apply molecular virial correction: transform atomic virial to molecular virial.
/// GROMOS: atomic_to_molecular_virial in prepare_virial.cc
///
/// For each pressure group: compute COM position (chain-gathered), then for each atom:
///   r = nearest_image(pos_atom, com_pos)
///   corrP(b, a) += force(a) * r(b)
/// Finally: virial -= corrP
///
/// Operates on conf.current() (GROMOS convention).
fn apply_molecular_virial_correction(topo: &Topology, conf: &mut Configuration, periodicity: &Periodicity) {
    let pos = &conf.current().pos;
    let force = &conf.current().force;

    let mut corr = [[0.0f64; 3]; 3];

    for pg in &topo.pressure_groups {
        // Chain-gather COM position (GROMOS: _centre_of_mass)
        let first = pg.start;
        let mut com = [0.0f64; 3];
        let mut m_total = 0.0f64;
        let mut prev = pos[first];

        for i in pg.clone() {
            let m = topo.mass[i];
            m_total += m;
            // nearest_image returns (pos[i] - prev) with PBC
            let p = periodicity.nearest_image(pos[i], prev);
            // gathered position = p + prev
            let gathered = p + prev;
            com[0] += m * gathered.x;
            com[1] += m * gathered.y;
            com[2] += m * gathered.z;
            prev = gathered; // chain: advance reference to gathered position
        }
        if m_total > 0.0 {
            com[0] /= m_total;
            com[1] /= m_total;
            com[2] /= m_total;
        }

        let com_vec = Vec3::new(com[0], com[1], com[2]);

        // Accumulate correction: corrP(b,a) += f(a) * (pos - com)(b)
        for i in pg.clone() {
            let r = periodicity.nearest_image(pos[i], com_vec);
            let f = force[i];
            let rv = [r.x, r.y, r.z];
            let fv = [f.x, f.y, f.z];
            for a in 0..3 {
                for b in 0..3 {
                    corr[b][a] += fv[a] * rv[b];
                }
            }
        }
    }

    // virial -= corrP
    let corr_mat = Mat3::from_cols(
        Vec3::new(corr[0][0], corr[1][0], corr[2][0]),
        Vec3::new(corr[0][1], corr[1][1], corr[2][1]),
        Vec3::new(corr[0][2], corr[1][2], corr[2][2]),
    );
    let state = conf.current_mut();
    state.virial_tensor = Mat3::from_cols(
        state.virial_tensor.x_axis - corr_mat.x_axis,
        state.virial_tensor.y_axis - corr_mat.y_axis,
        state.virial_tensor.z_axis - corr_mat.z_axis,
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::math::{Vacuum, Vec3};
    use gromos_core::pairlist::StandardPairlistAlgorithm;
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
        topo.exclusions = vec![Vec::new(), Vec::new()];

        // LJ parameter matrix (1x1 for single atom type)
        let lj = LJParameters { c6: 6.2647e-3, c12: 9.847e-6, cs6: 6.2647e-3, cs12: 9.847e-6 };
        let lj_params = vec![vec![lj]];

        // CRF params (epsilon=1, rf_epsilon=1 => pure Coulomb cutoff, but charges=0)
        let crf_params = CRFParameters::new(1.4, 1.0, 1.0, 0.0);

        let periodicity = Periodicity::Vacuum(Vacuum);

        // Pairlist: manually set to contain the single pair
        let mut pairlist = PairlistContainer::new(1.4, 1.4, 0.0);
        pairlist.solute_short = vec![(0, 1)];
        pairlist.update_frequency = 1000; // Don't auto-update in test

        let pairlist_algorithm = PairlistAlgorithm::Standard(StandardPairlistAlgorithm::new(false));

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
