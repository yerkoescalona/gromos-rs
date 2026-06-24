//! Nonbonded interactions: Lennard-Jones and Coulomb (CRF)

pub mod params;
pub mod excluded;
pub mod innerloops;
pub mod perturbed;

pub use params::*;
pub use excluded::*;
pub use innerloops::*;
pub use perturbed::*;

/// MPI-parallel nonbonded force calculation
///
/// This module implements the MPI master-slave pattern from
/// md++/src/interaction/nonbonded/interaction/mpi_nonbonded_master.cc
#[cfg(feature = "use-mpi")]
pub mod mpi_nonbonded {
    use super::*;
    use crate::mpi::{MpiControl, MpiNonbondedMaster};

    /// Calculate nonbonded interactions with MPI parallelization
    ///
    /// Based on `MPI_Nonbonded_Master::calculate_interactions` in mpi_nonbonded_master.cc
    ///
    /// # Workflow
    ///
    /// 1. **Master broadcasts data** to all processes:
    ///    - Atomic positions
    ///    - Charges (for electrostatics)
    ///    - Box dimensions (for PBC)
    ///    - Lambda (for free energy calculations)
    ///
    /// 2. **All processes calculate** their assigned subset of pair interactions
    ///
    /// 3. **Master reduces results** from all processes:
    ///    - Forces (summed via MPI_Reduce)
    ///    - Energies (summed via MPI_Reduce)
    ///    - Virial tensor (summed via MPI_Reduce)
    ///
    /// # Arguments
    ///
    /// * `world` - MPI world communicator
    /// * `mpi_control` - MPI control structure with rank/size info
    /// * `positions` - Atomic positions (modified by broadcast on slaves)
    /// * `charges` - Atomic charges (modified by broadcast on slaves)
    /// * `iac` - Integer atom codes (atom types)
    /// * `pairlist` - List of atom pairs to calculate
    /// * `lj_params` - Lennard-Jones parameter matrix
    /// * `crf` - Coulomb Reaction Field parameters
    /// * `periodicity` - Periodic boundary conditions
    /// * `box_matrix` - Box dimensions (modified by broadcast on slaves)
    /// * `lambda` - Free energy lambda parameter (modified by broadcast on slaves)
    ///
    /// # Returns
    ///
    /// Force storage with forces, energies, and virial (only valid on master process)
    pub fn calculate_mpi<BC: gromos_core::math::BoundaryCondition>(
        world: &mpi::topology::SystemCommunicator,
        mpi_control: &MpiControl,
        positions: &mut [gromos_core::math::Vec3],
        charges: &mut [f64],
        iac: &[u32],
        pairlist: &[(u32, u32)],
        lj_params: &LJParamMatrix,
        crf: &CRFParameters,
        periodicity: &BC,
        box_matrix: &mut [f64; 9],
        lambda: &mut f64,
    ) -> ForceStorage {
        let mpi_master = MpiNonbondedMaster::new(mpi_control.clone());
        let n_atoms = positions.len();

        // ===================================================================
        // Phase 1: BROADCAST - Master sends data to all processes
        // ===================================================================
        // (Equivalent to mpi_nonbonded_master.cc:127-150)

        // Broadcast positions
        mpi_master.broadcast_positions(world, positions);

        // Broadcast charges (needed for QM/MM and dynamic charges)
        mpi_master.broadcast_charges(world, charges);

        // Broadcast box dimensions (needed for PBC)
        mpi_master.broadcast_box(world, box_matrix);

        // Broadcast lambda (for free energy perturbation)
        mpi_master.broadcast_lambda(world, lambda);

        // ===================================================================
        // Phase 2: COMPUTE - Each process calculates its subset of pairs
        // ===================================================================
        // (Equivalent to mpi_nonbonded_master.cc:152-163)

        // Divide pairlist among processes
        let (pair_start, pair_end) = mpi_master.calculate_pair_range(pairlist.len());
        let local_pairlist = &pairlist[pair_start..pair_end];

        // Calculate local forces using serial innerloop
        let mut local_storage = ForceStorage::new(n_atoms);
        lj_crf_innerloop(
            positions,
            charges,
            iac,
            local_pairlist,
            lj_params,
            crf,
            periodicity,
            &mut local_storage,
        );

        // ===================================================================
        // Phase 3: REDUCE - Master collects results from all processes
        // ===================================================================
        // (Equivalent to mpi_nonbonded_master.cc:165-240)

        let mut final_storage = ForceStorage::new(n_atoms);

        // Reduce forces (MPI_SUM)
        mpi_master.reduce_forces(world, &local_storage.forces, &mut final_storage.forces);

        // Reduce energies (MPI_SUM)
        let local_energies = [local_storage.e_lj, local_storage.e_crf];
        let mut total_energies = [0.0, 0.0];
        mpi_master.reduce_energies(world, &local_energies, &mut total_energies);

        final_storage.e_lj = total_energies[0];
        final_storage.e_crf = total_energies[1];

        // Reduce virial tensor (MPI_SUM)
        let local_virial_flat: [f64; 9] = [
            local_storage.virial[0][0],
            local_storage.virial[0][1],
            local_storage.virial[0][2],
            local_storage.virial[1][0],
            local_storage.virial[1][1],
            local_storage.virial[1][2],
            local_storage.virial[2][0],
            local_storage.virial[2][1],
            local_storage.virial[2][2],
        ];
        let mut total_virial_flat = [0.0; 9];
        mpi_master.reduce_virial(world, &local_virial_flat, &mut total_virial_flat);

        final_storage.virial = [
            [
                total_virial_flat[0],
                total_virial_flat[1],
                total_virial_flat[2],
            ],
            [
                total_virial_flat[3],
                total_virial_flat[4],
                total_virial_flat[5],
            ],
            [
                total_virial_flat[6],
                total_virial_flat[7],
                total_virial_flat[8],
            ],
        ];

        // Only the master process has valid results
        // Slaves have garbage in final_storage (not written to by reduce)
        final_storage
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_pair_range_calculation() {
            let control = MpiControl {
                enabled: true,
                rank: 0,
                size: 4,
                master_id: 0,
                simulation_id: 0,
            };

            let master = MpiNonbondedMaster::new(control);

            // Test with 100 pairs
            let (start, end) = master.calculate_pair_range(100);
            assert_eq!(start, 0);
            assert_eq!(end, 25);

            // Second process
            let control2 = MpiControl {
                enabled: true,
                rank: 1,
                size: 4,
                master_id: 0,
                simulation_id: 0,
            };
            let master2 = MpiNonbondedMaster::new(control2);
            let (start2, end2) = master2.calculate_pair_range(100);
            assert_eq!(start2, 25);
            assert_eq!(end2, 50);
        }
    }
}
