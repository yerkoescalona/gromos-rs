//! MPI support for distributed-memory parallelism
//!
//! This module provides MPI functionality for multi-node molecular dynamics simulations,
//! based on the implementation in md++/src/simulation/mpiControl.h and
//! md++/src/interaction/nonbonded/interaction/mpi_nonbonded_master.cc
//!
//! # Architecture
//!
//! The implementation follows a master-slave pattern:
//! - **Master process** (rank 0): Coordinates simulation, broadcasts data, collects results
//! - **Slave processes** (rank 1..N-1): Receive data, compute subset of interactions, send results
//!
//! # Workflow
//!
//! 1. **Initialization**: MPI_Init, determine rank/size, setup communicators
//! 2. **Distribution**: Master broadcasts positions, charges, box, parameters
//! 3. **Computation**: Each process calculates its assigned subset of nonbonded interactions
//! 4. **Collection**: Master reduces forces and energies from all processes
//!
//! # Example
//!
//! ```rust,no_run
//! use gromos_rs::mpi::*;
//!
//! // Initialize MPI
//! let universe = mpi::initialize().unwrap();
//! let world = universe.world();
//! let mpi_control = MpiControl::new(&world);
//!
//! if mpi_control.is_master() {
//!     println!("Master process with {} total processes", mpi_control.size);
//! }
//!
//! // ... run simulation with MPI support ...
//!
//! // MPI automatically finalizes when universe goes out of scope
//! ```

#[cfg(feature = "use-mpi")]
pub use mpi_impl::*;

#[cfg(feature = "use-mpi")]
mod mpi_impl {
    use gromos_core::math::Vec3;
    use mpi::traits::*;

    /// MPI control structure
    ///
    /// Equivalent to `simulation::MpiControl` in md++/src/simulation/mpiControl.h
    ///
    /// This structure holds all MPI-related information for a simulation:
    /// - Process rank and total number of processes
    /// - Master process ID
    /// - Communicator for this simulation
    #[derive(Debug, Clone)]
    pub struct MpiControl {
        pub enabled: bool,
        pub rank: i32,
        pub size: i32,
        pub master_id: i32,
        pub simulation_id: i32,
    }

    impl MpiControl {
        /// Create MPI control from world communicator
        pub fn new(world: &mpi::topology::SystemCommunicator) -> Self {
            let rank = world.rank();
            let size = world.size();

            Self {
                enabled: true,
                rank,
                size,
                master_id: 0,
                simulation_id: 0,
            }
        }

        /// Create disabled MPI control (serial execution)
        pub fn serial() -> Self {
            Self {
                enabled: false,
                rank: 0,
                size: 1,
                master_id: 0,
                simulation_id: 0,
            }
        }

        /// Check if this process is the master
        #[inline]
        pub fn is_master(&self) -> bool {
            self.rank == self.master_id
        }

        /// Check if this process is a slave
        #[inline]
        pub fn is_slave(&self) -> bool {
            self.rank != self.master_id
        }

        /// Get the number of slave processes
        #[inline]
        pub fn num_slaves(&self) -> i32 {
            if self.size > 0 {
                self.size - 1
            } else {
                0
            }
        }
    }

    /// MPI nonbonded master calculator
    ///
    /// Equivalent to `interaction::MPI_Nonbonded_Master` in
    /// md++/src/interaction/nonbonded/interaction/mpi_nonbonded_master.cc
    ///
    /// This handles the distribution of nonbonded force calculations across MPI processes:
    ///
    /// 1. **Broadcast phase**: Master sends positions, charges, box to all processes
    /// 2. **Compute phase**: Each process calculates its subset of pair interactions
    /// 3. **Reduce phase**: Master collects forces and energies from all processes
    pub struct MpiNonbondedMaster {
        pub control: MpiControl,
    }

    impl MpiNonbondedMaster {
        /// Create new MPI nonbonded master
        pub fn new(control: MpiControl) -> Self {
            Self { control }
        }

        /// Broadcast positions to all processes
        ///
        /// Equivalent to MPI_Bcast in mpi_nonbonded_master.cc:127-130
        pub fn broadcast_positions(
            &self,
            world: &mpi::topology::SystemCommunicator,
            positions: &mut [Vec3],
        ) {
            let root_process = world.process_at_rank(self.control.master_id);

            // Convert Vec3 slice to flat f32 slice for MPI
            let pos_flat: &mut [f32] = unsafe {
                std::slice::from_raw_parts_mut(
                    positions.as_mut_ptr() as *mut f32,
                    positions.len() * 3,
                )
            };

            root_process.broadcast_into(pos_flat);
        }

        /// Broadcast charges to all processes
        ///
        /// Equivalent to MPI_Bcast in mpi_nonbonded_master.cc:133-136
        pub fn broadcast_charges(
            &self,
            world: &mpi::topology::SystemCommunicator,
            charges: &mut [f64],
        ) {
            let root_process = world.process_at_rank(self.control.master_id);
            root_process.broadcast_into(charges);
        }

        /// Broadcast box dimensions to all processes
        ///
        /// Equivalent to MPI_Bcast in mpi_nonbonded_master.cc:142-144
        pub fn broadcast_box(
            &self,
            world: &mpi::topology::SystemCommunicator,
            box_matrix: &mut [f32; 9],
        ) {
            let root_process = world.process_at_rank(self.control.master_id);
            root_process.broadcast_into(box_matrix);
        }

        /// Broadcast lambda parameter (for free energy calculations)
        ///
        /// Equivalent to MPI_Bcast in mpi_nonbonded_master.cc:147-150
        pub fn broadcast_lambda(
            &self,
            world: &mpi::topology::SystemCommunicator,
            lambda: &mut f64,
        ) {
            let root_process = world.process_at_rank(self.control.master_id);
            root_process.broadcast_into(std::slice::from_mut(lambda));
        }

        /// Reduce forces from all processes to master
        ///
        /// Equivalent to MPI_Reduce in mpi_nonbonded_master.cc:168-173
        ///
        /// Each process has calculated a subset of forces. This sums them all
        /// and stores the result on the master process.
        pub fn reduce_forces(
            &self,
            world: &mpi::topology::SystemCommunicator,
            local_forces: &[Vec3],
            total_forces: &mut [Vec3],
        ) {
            let root_process = world.process_at_rank(self.control.master_id);

            // Convert Vec3 slices to flat f32 slices for MPI
            let local_flat: &[f32] = unsafe {
                std::slice::from_raw_parts(
                    local_forces.as_ptr() as *const f32,
                    local_forces.len() * 3,
                )
            };

            let total_flat: &mut [f32] = unsafe {
                std::slice::from_raw_parts_mut(
                    total_forces.as_mut_ptr() as *mut f32,
                    total_forces.len() * 3,
                )
            };

            // Sum forces from all processes
            root_process.reduce_into_root(
                local_flat,
                total_flat,
                mpi::collective::SystemOperation::sum(),
            );
        }

        /// Reduce energies from all processes to master
        ///
        /// Similar to force reduction but for scalar energy values
        pub fn reduce_energies(
            &self,
            world: &mpi::topology::SystemCommunicator,
            local_energies: &[f64],
            total_energies: &mut [f64],
        ) {
            let root_process = world.process_at_rank(self.control.master_id);
            root_process.reduce_into_root(
                local_energies,
                total_energies,
                mpi::collective::SystemOperation::sum(),
            );
        }

        /// Reduce virial tensor from all processes to master
        pub fn reduce_virial(
            &self,
            world: &mpi::topology::SystemCommunicator,
            local_virial: &[f64; 9],
            total_virial: &mut [f64; 9],
        ) {
            let root_process = world.process_at_rank(self.control.master_id);
            root_process.reduce_into_root(
                &local_virial[..],
                &mut total_virial[..],
                mpi::collective::SystemOperation::sum(),
            );
        }

        /// Calculate pair range for this process
        ///
        /// Divides the pairlist among processes. Each process gets approximately
        /// equal number of pairs to calculate.
        ///
        /// # Arguments
        /// * `total_pairs` - Total number of atom pairs to calculate
        ///
        /// # Returns
        /// * `(start, end)` - Range of pairs for this process [start, end)
        pub fn calculate_pair_range(&self, total_pairs: usize) -> (usize, usize) {
            let pairs_per_process = total_pairs / self.control.size as usize;
            let remainder = total_pairs % self.control.size as usize;

            let start = self.control.rank as usize * pairs_per_process
                + std::cmp::min(self.control.rank as usize, remainder);
            let end = start
                + pairs_per_process
                + if (self.control.rank as usize) < remainder {
                    1
                } else {
                    0
                };

            (start, end)
        }
    }

    /// Helper function to check if MPI is initialized
    pub fn is_initialized() -> bool {
        mpi::is_initialized()
    }

    /// Print MPI information (for debugging)
    pub fn print_info(control: &MpiControl) {
        if control.is_master() {
            println!("=== MPI Configuration ===");
            println!("  Total processes: {}", control.size);
            println!("  Master rank: {}", control.master_id);
            println!("  Slave processes: {}", control.num_slaves());
            println!("========================");
        } else {
            println!("Slave process {} ready", control.rank);
        }
    }
}

#[cfg(not(feature = "use-mpi"))]
mod mpi_impl {
    //! Dummy MPI implementation when MPI feature is disabled
    //!
    //! Provides the same API but only supports serial execution

    /// Dummy MPI control for serial execution
    #[derive(Debug, Clone)]
    pub struct MpiControl {
        pub enabled: bool,
        pub rank: i32,
        pub size: i32,
        pub master_id: i32,
        pub simulation_id: i32,
    }

    impl MpiControl {
        /// Create serial-only control
        pub fn serial() -> Self {
            Self {
                enabled: false,
                rank: 0,
                size: 1,
                master_id: 0,
                simulation_id: 0,
            }
        }

        #[inline]
        pub fn is_master(&self) -> bool {
            true
        }

        #[inline]
        pub fn is_slave(&self) -> bool {
            false
        }

        #[inline]
        pub fn num_slaves(&self) -> i32 {
            0
        }
    }

    /// Dummy nonbonded master for serial execution
    pub struct MpiNonbondedMaster {
        pub control: MpiControl,
    }

    impl MpiNonbondedMaster {
        pub fn new(control: MpiControl) -> Self {
            Self { control }
        }
    }

    pub fn is_initialized() -> bool {
        false
    }

    pub fn print_info(control: &MpiControl) {
        println!("MPI not enabled (serial execution)");
        println!("  Rank: {}", control.rank);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "use-mpi")]
    #[test]
    fn test_serial_control() {
        let control = MpiControl::serial();
        assert!(!control.enabled);
        assert_eq!(control.rank, 0);
        assert_eq!(control.size, 1);
        assert!(control.is_master());
        assert!(!control.is_slave());
    }

    #[cfg(feature = "use-mpi")]
    #[test]
    fn test_pair_range_division() {
        let control = MpiControl {
            enabled: true,
            rank: 0,
            size: 4,
            master_id: 0,
            simulation_id: 0,
        };

        let master = MpiNonbondedMaster::new(control.clone());

        // Test with 100 pairs divided among 4 processes
        let (start, end) = master.calculate_pair_range(100);
        assert_eq!(start, 0);
        assert_eq!(end, 25);

        // Test with 103 pairs (not evenly divisible)
        let (start, end) = master.calculate_pair_range(103);
        assert_eq!(start, 0);
        assert_eq!(end, 26); // First process gets 26 (extra from remainder)
    }
}
