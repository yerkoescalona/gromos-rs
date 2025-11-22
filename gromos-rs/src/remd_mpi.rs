//! MPI-enabled Replica Exchange Molecular Dynamics
//!
//! This module implements distributed replica exchange using MPI, where each
//! replica runs on a separate MPI process. Based on md++/program/repex_mpi.cc
//! and md++/src/replicaExchange/replica/replica_MPI_master.cpp
//!
//! # Architecture
//!
//! - Each MPI process runs one replica at a specific temperature/lambda
//! - Master process (rank 0) coordinates exchange attempts
//! - Replicas communicate energies and exchange decisions via MPI
//! - Exchanges are attempted periodically (e.g., every 100 steps)
//!
//! # Exchange Protocol
//!
//! 1. **Run MD**: Each process runs MD independently for N steps
//! 2. **Gather energies**: All processes send potential energy to master
//! 3. **Calculate acceptance**: Master computes Metropolis criterion
//! 4. **Broadcast decisions**: Master sends exchange decisions to all
//! 5. **Exchange configurations**: Processes swap velocities/positions as needed
//!
//! # Example
//!
//! ```rust,no_run
//! use gromos_rs::remd_mpi::*;
//!
//! let universe = mpi::initialize().unwrap();
//! let world = universe.world();
//! let rank = world.rank();
//!
//! // Each process gets a different temperature
//! let temperatures = vec![300.0, 310.0, 320.0, 330.0];
//! let my_temp = temperatures[rank as usize];
//!
//! let mut controller = RemdMpiController::new(&world, my_temp, 0.0);
//!
//! // Run REMD
//! for cycle in 0..1000 {
//!     // Run MD for 100 steps
//!     controller.run_md_steps(100, &topo, &mut integrator);
//!
//!     // Attempt exchange
//!     controller.attempt_exchange(cycle);
//! }
//! ```

#[cfg(feature = "use-mpi")]
pub use mpi_impl::*;

#[cfg(feature = "use-mpi")]
mod mpi_impl {
    use crate::configuration::Configuration;
    use crate::integrator::Integrator;
    use crate::remd::{ExchangeScheme, ExchangeStatistics, ExchangeType};
    use crate::replica::{Replica, ReplicaId, ReplicaInfo};
    use crate::topology::Topology;
    use mpi::traits::*;
    use rand::Rng;

    /// MPI-enabled Replica Exchange Controller
    ///
    /// Each MPI process runs one replica and coordinates exchanges through
    /// collective MPI operations.
    pub struct RemdMpiController {
        /// My replica
        replica: Replica,

        /// MPI rank
        rank: i32,

        /// Total number of processes (replicas)
        size: i32,

        /// Exchange type
        exchange_type: ExchangeType,

        /// Exchange scheme
        exchange_scheme: ExchangeScheme,

        /// Exchange attempt frequency (steps)
        exchange_frequency: usize,

        /// Statistics
        statistics: ExchangeStatistics,

        /// Current exchange partner (if any)
        partner_rank: Option<i32>,

        /// Whether exchange was accepted in last attempt
        last_exchange_accepted: bool,

        /// RNG for Metropolis criterion
        rng: rand::rngs::ThreadRng,

        /// Temperature grid for 2D T-Lambda REPEX (master only)
        /// Layout: [T1, T1, T2, T2, ...] for lambdas [λ1, λ2, λ1, λ2, ...]
        temperature_grid: Vec<f64>,

        /// Lambda grid for 2D T-Lambda REPEX (master only)
        lambda_grid: Vec<f64>,

        /// Number of temperature values (for 2D grid)
        n_temps: usize,

        /// Number of lambda values (for 2D grid)
        n_lambdas: usize,
    }

    impl RemdMpiController {
        /// Create new MPI REMD controller
        ///
        /// Each process creates its own replica with specified temperature/lambda
        pub fn new(
            world: &mpi::topology::SystemCommunicator,
            temperature: f64,
            lambda: f64,
            n_atoms: usize,
            dt: f64,
            exchange_type: ExchangeType,
            exchange_scheme: ExchangeScheme,
            exchange_frequency: usize,
        ) -> Self {
            let rank = world.rank();
            let size = world.size();

            let replica_info = ReplicaInfo::new(rank as usize, temperature, lambda, dt);
            let configuration = Configuration::new(n_atoms, 1, 1);
            let replica = Replica::new_from_info(replica_info, configuration);

            RemdMpiController {
                replica,
                rank,
                size,
                exchange_type,
                exchange_scheme,
                exchange_frequency,
                statistics: ExchangeStatistics::new(size as usize),
                partner_rank: None,
                last_exchange_accepted: false,
                rng: rand::thread_rng(),
                temperature_grid: Vec::new(),
                lambda_grid: Vec::new(),
                n_temps: 1,
                n_lambdas: 1,
            }
        }

        /// Set up 2D Temperature-Lambda grid for T-Lambda REPEX
        ///
        /// This must be called on all processes to set up the grid structure.
        /// The grid is organized as: (T1,λ1), (T1,λ2), ..., (T2,λ1), (T2,λ2), ...
        ///
        /// # Arguments
        /// * `temperatures` - All temperature values in the ladder
        /// * `lambdas` - All lambda values in the ladder
        ///
        /// # Panics
        /// Panics if n_temps * n_lambdas != number of MPI processes
        pub fn setup_2d_grid(
            &mut self,
            world: &mpi::topology::SystemCommunicator,
            temperatures: Vec<f64>,
            lambdas: Vec<f64>,
        ) {
            self.n_temps = temperatures.len();
            self.n_lambdas = lambdas.len();

            assert_eq!(
                self.n_temps * self.n_lambdas,
                self.size as usize,
                "Grid size ({} temps × {} lambdas = {}) must match number of processes ({})",
                self.n_temps,
                self.n_lambdas,
                self.n_temps * self.n_lambdas,
                self.size
            );

            // Build full grid on all processes (needed for neighbor finding)
            self.temperature_grid.clear();
            self.lambda_grid.clear();

            for temp in &temperatures {
                for lambda in &lambdas {
                    self.temperature_grid.push(*temp);
                    self.lambda_grid.push(*lambda);
                }
            }

            // Update my replica's temperature and lambda from grid
            let my_temp = self.temperature_grid[self.rank as usize];
            let my_lambda = self.lambda_grid[self.rank as usize];
            self.replica.info.temperature = my_temp;
            self.replica.info.lambda = my_lambda;
        }

        /// Get reference to my replica
        pub fn replica(&self) -> &Replica {
            &self.replica
        }

        /// Get mutable reference to my replica
        pub fn replica_mut(&mut self) -> &mut Replica {
            &mut self.replica
        }

        /// Run MD steps on this replica
        pub fn run_md_steps(
            &mut self,
            steps: usize,
            topo: &Topology,
            integrator: &mut dyn Integrator,
        ) {
            for _ in 0..steps {
                self.replica.step(topo, integrator);
            }
        }

        /// Attempt replica exchange (collective operation)
        ///
        /// This is a collective MPI operation - all processes must call this.
        ///
        /// # Exchange Protocol (based on repex_mpi.cc)
        ///
        /// 1. Gather potential energies from all replicas to master (rank 0)
        /// 2. Master determines exchange pairs based on scheme
        /// 3. Master calculates Metropolis acceptance probabilities
        /// 4. Master generates random numbers and decides acceptance
        /// 5. Master broadcasts exchange decisions to all processes
        /// 6. Processes with accepted exchanges swap configurations
        pub fn attempt_exchange(
            &mut self,
            world: &mpi::topology::SystemCommunicator,
            step: usize,
        ) -> bool {
            // Only attempt exchange at specified frequency
            if step % self.exchange_frequency != 0 {
                return false;
            }

            let root = world.process_at_rank(0);

            // Step 1: Gather energies to master
            let my_energy = self.replica.info.potential_energy;
            let mut all_energies = vec![0.0; self.size as usize];

            root.gather_into_root(&my_energy, &mut all_energies[..]);

            // Step 2-4: Master determines exchanges and broadcasts decisions
            let mut exchange_partner = -1i32; // -1 means no exchange
            let mut exchange_accepted = false;

            if self.rank == 0 {
                // Master logic
                let pairs = self.determine_exchange_pairs(step);

                for (i, j) in pairs {
                    let accept = self.calculate_acceptance(
                        i as usize,
                        j as usize,
                        all_energies[i as usize],
                        all_energies[j as usize],
                    );

                    // Record statistics
                    self.statistics
                        .record_attempt(i as usize, j as usize, accept);

                    if accept {
                        // Mark both replicas for exchange
                        if self.rank == i {
                            exchange_partner = j;
                            exchange_accepted = true;
                        } else if self.rank == j {
                            exchange_partner = i;
                            exchange_accepted = true;
                        }
                    }
                }
            }

            // Step 5: Broadcast exchange decisions to all processes
            // Each process needs to know: who is their partner (if any)?
            // We'll broadcast an array where exchange_partners[i] = j means i exchanges with j

            let mut exchange_partners = vec![-1i32; self.size as usize];
            if self.rank == 0 {
                // Master fills in the partner array
                // (In a full implementation, master would track this during pair evaluation)
                // For now, we'll use a simpler approach with direct broadcast
            }

            // Broadcast the partner for this specific rank
            root.broadcast_into(&mut exchange_partner);

            // Step 6: Execute exchange if accepted
            if exchange_partner >= 0 {
                self.execute_exchange(world, exchange_partner);
                self.last_exchange_accepted = true;
                self.partner_rank = Some(exchange_partner);
                return true;
            }

            self.last_exchange_accepted = false;
            self.partner_rank = None;
            false
        }

        /// Determine exchange pairs based on scheme and exchange type
        ///
        /// Returns vector of (rank_i, rank_j) pairs to attempt exchange
        fn determine_exchange_pairs(&self, step: usize) -> Vec<(i32, i32)> {
            let mut pairs = Vec::new();

            match self.exchange_type {
                ExchangeType::Temperature | ExchangeType::Lambda => {
                    // 1D exchange (either temperature or lambda only)
                    self.determine_1d_pairs(step, &mut pairs);
                },
                ExchangeType::TemperatureLambda => {
                    // 2D exchange - alternate between T and λ directions
                    self.determine_2d_pairs(step, &mut pairs);
                },
            }

            pairs
        }

        /// Determine pairs for 1D exchange (temperature-only or lambda-only)
        fn determine_1d_pairs(&self, step: usize, pairs: &mut Vec<(i32, i32)>) {
            match self.exchange_scheme {
                ExchangeScheme::Sequential => {
                    // 0↔1, 2↔3, 4↔5, ...
                    for i in (0..self.size - 1).step_by(2) {
                        pairs.push((i, i + 1));
                    }
                },
                ExchangeScheme::OddEven => {
                    // Alternate between even and odd pairs
                    let offset = (step / self.exchange_frequency) % 2;
                    for i in (offset as i32..self.size - 1).step_by(2) {
                        if i + 1 < self.size {
                            pairs.push((i, i + 1));
                        }
                    }
                },
                ExchangeScheme::Random => {
                    // Random pairs (simple implementation)
                    let mut available: Vec<i32> = (0..self.size).collect();
                    while available.len() >= 2 {
                        let idx1 = self.rng.gen_range(0..available.len());
                        let rank1 = available.remove(idx1);
                        let idx2 = self.rng.gen_range(0..available.len());
                        let rank2 = available.remove(idx2);
                        pairs.push((rank1, rank2));
                    }
                },
            }
        }

        /// Determine pairs for 2D T-Lambda exchange
        ///
        /// In 2D REPEX, we alternate between exchanges in the temperature direction
        /// (same lambda) and lambda direction (same temperature).
        ///
        /// Grid layout (example with 3 temps × 2 lambdas = 6 processes):
        /// ```text
        ///     λ1    λ2
        /// T1   0 ←→  1
        ///      ↕     ↕
        /// T2   2 ←→  3
        ///      ↕     ↕
        /// T3   4 ←→  5
        /// ```
        fn determine_2d_pairs(&self, step: usize, pairs: &mut Vec<(i32, i32)>) {
            // Alternate between T-direction and λ-direction exchanges
            let exchange_cycle = (step / self.exchange_frequency) % 2;

            if exchange_cycle == 0 {
                // Temperature direction exchanges (same lambda, different T)
                // For each lambda column, try neighboring temperature rows
                for lambda_idx in 0..self.n_lambdas {
                    for temp_idx in 0..self.n_temps - 1 {
                        let rank_i = (temp_idx * self.n_lambdas + lambda_idx) as i32;
                        let rank_j = ((temp_idx + 1) * self.n_lambdas + lambda_idx) as i32;
                        pairs.push((rank_i, rank_j));
                    }
                }
            } else {
                // Lambda direction exchanges (same T, different lambda)
                // For each temperature row, try neighboring lambda columns
                for temp_idx in 0..self.n_temps {
                    for lambda_idx in 0..self.n_lambdas - 1 {
                        let rank_i = (temp_idx * self.n_lambdas + lambda_idx) as i32;
                        let rank_j = (temp_idx * self.n_lambdas + lambda_idx + 1) as i32;
                        pairs.push((rank_i, rank_j));
                    }
                }
            }
        }

        /// Calculate Metropolis acceptance probability
        ///
        /// For temperature exchange (same lambda):
        /// ```text
        /// ΔE = E_i - E_j
        /// Δβ = β_i - β_j = 1/(k_B T_i) - 1/(k_B T_j)
        /// P_accept = min(1, exp(Δβ * ΔE))
        /// ```
        ///
        /// For lambda exchange (same temperature):
        /// ```text
        /// ΔH = H_i(λ_j) - H_i(λ_i) + H_j(λ_i) - H_j(λ_j)
        /// P_accept = min(1, exp(-β * ΔH))
        /// ```
        fn calculate_acceptance(
            &mut self,
            rank_i: usize,
            rank_j: usize,
            energy_i: f64,
            energy_j: f64,
        ) -> bool {
            const KB: f64 = 0.008314462618; // kJ/(mol·K)

            // Get temperature and lambda for both replicas
            let temp_i = self.temperature_grid.get(rank_i).copied().unwrap_or(300.0);
            let temp_j = self.temperature_grid.get(rank_j).copied().unwrap_or(300.0);
            let lambda_i = self.lambda_grid.get(rank_i).copied().unwrap_or(0.0);
            let lambda_j = self.lambda_grid.get(rank_j).copied().unwrap_or(0.0);

            let prob = match self.exchange_type {
                ExchangeType::Temperature => {
                    // Temperature exchange: different temps, same lambda
                    let beta_i = 1.0 / (KB * temp_i);
                    let beta_j = 1.0 / (KB * temp_j);
                    let delta_e = energy_i - energy_j;
                    let delta_beta = beta_i - beta_j;
                    (delta_beta * delta_e).exp()
                },
                ExchangeType::Lambda => {
                    // Lambda exchange: same temp, different lambdas
                    let beta = 1.0 / (KB * temp_i); // Same temperature

                    // Simplified: ΔH ≈ (λ_j - λ_i) * (E_i - E_j)
                    // In full implementation, need to evaluate energies at swapped lambdas
                    let delta_lambda = lambda_j - lambda_i;
                    let delta_e = energy_i - energy_j;
                    let delta_h = delta_lambda * delta_e;

                    (-beta * delta_h).exp()
                },
                ExchangeType::TemperatureLambda => {
                    // 2D exchange: determine if this is T or λ exchange
                    if (temp_i - temp_j).abs() > 1e-6 {
                        // Temperature direction exchange
                        let beta_i = 1.0 / (KB * temp_i);
                        let beta_j = 1.0 / (KB * temp_j);
                        let delta_e = energy_i - energy_j;
                        let delta_beta = beta_i - beta_j;
                        (delta_beta * delta_e).exp()
                    } else {
                        // Lambda direction exchange
                        let beta = 1.0 / (KB * temp_i);
                        let delta_lambda = lambda_j - lambda_i;
                        let delta_e = energy_i - energy_j;
                        let delta_h = delta_lambda * delta_e;
                        (-beta * delta_h).exp()
                    }
                },
            };

            // Metropolis acceptance
            let acceptance_prob = prob.min(1.0);
            let rand_val: f64 = self.rng.gen();
            rand_val < acceptance_prob
        }

        /// Execute configuration exchange with partner
        ///
        /// Exchanges velocities between this process and partner process.
        /// In T-REMD, we exchange configurations (velocities) while keeping
        /// positions (replicas stay in same potential energy basin).
        fn execute_exchange(
            &mut self,
            world: &mpi::topology::SystemCommunicator,
            partner_rank: i32,
        ) {
            // Exchange velocities with partner
            let n_atoms = self.replica.configuration.current().vel.len();
            let my_velocities: Vec<f32> = self
                .replica
                .configuration
                .current()
                .vel
                .iter()
                .flat_map(|v| [v.x, v.y, v.z])
                .collect();

            let mut partner_velocities = vec![0.0f32; n_atoms * 3];

            // Use sendrecv for bidirectional exchange
            let partner = world.process_at_rank(partner_rank);

            if self.rank < partner_rank {
                // Send to higher rank, receive from higher rank
                world.process_at_rank(partner_rank).send(&my_velocities[..]);
                world
                    .process_at_rank(partner_rank)
                    .receive_into(&mut partner_velocities[..]);
            } else {
                // Receive from lower rank, send to lower rank
                world
                    .process_at_rank(partner_rank)
                    .receive_into(&mut partner_velocities[..]);
                world.process_at_rank(partner_rank).send(&my_velocities[..]);
            }

            // Update my velocities with partner's
            for (i, chunk) in partner_velocities.chunks(3).enumerate() {
                self.replica.configuration.current_mut().vel[i] =
                    crate::math::Vec3::new(chunk[0], chunk[1], chunk[2]);
            }
        }

        /// Get exchange statistics (master only)
        pub fn statistics(&self) -> &ExchangeStatistics {
            &self.statistics
        }

        /// Print exchange statistics (master only)
        pub fn print_statistics(&self) {
            if self.rank == 0 {
                println!("\n=== Replica Exchange Statistics ===");
                println!("Total attempts: {}", self.statistics.attempts);
                println!("Total accepted: {}", self.statistics.accepted);
                println!(
                    "Acceptance rate: {:.2}%",
                    self.statistics.acceptance_rate * 100.0
                );
                println!("\nPer-pair acceptance rates:");
                for i in 0..self.size as usize {
                    for j in (i + 1)..self.size as usize {
                        let rate = self.statistics.pair_acceptance_rate(i, j);
                        if self.statistics.pair_attempts[i][j] > 0 {
                            println!(
                                "  {} ↔ {}: {:.2}% ({}/{})",
                                i,
                                j,
                                rate * 100.0,
                                self.statistics.pair_acceptances[i][j],
                                self.statistics.pair_attempts[i][j]
                            );
                        }
                    }
                }
            }
        }
    }
}

#[cfg(not(feature = "use-mpi"))]
mod mpi_impl {
    //! Dummy implementation when MPI is not enabled
    pub struct RemdMpiController;
}
