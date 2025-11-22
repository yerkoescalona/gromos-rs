//! GROMOS-RS: High-performance Rust molecular dynamics engine
//!
//! A complete Rust translation of GROMOS (GROningen MOlecular Simulation),
//! providing 2-3x performance improvements through modern optimization techniques.
//!
//! ## Features
//!
//! - **Pure Rust MD engine**: Complete translation, not just FFI wrappers
//! - **SIMD vectorization**: Automatic vectorization with glam + wide
//! - **Fearless concurrency**: Parallel algorithms with Rayon
//! - **Memory safety**: Zero-cost abstractions without runtime overhead
//! - **Modular design**: Easy to extend and customize
//!
//! ## Architecture
//!
//! ```text
//! gromos-rs
//! ├── math          - SIMD-accelerated math primitives
//! ├── topology      - Molecular structure and force field parameters
//! ├── configuration - System state management
//! ├── pairlist      - Neighbor list generation
//! ├── interaction   - Force calculations (bonded + nonbonded)
//! ├── integrator    - Time integration algorithms
//! └── engine        - Main MD simulation driver
//! ```
//!
//! ## Example
//!
//! ```rust,no_run
//! use gromos_rs::*;
//!
//! // Create system
//! let topo = topology::Topology::new();
//! let conf = configuration::Configuration::new(100, 1, 1);
//!
//! // Setup MD
//! let integrator = integrator::LeapFrog::new();
//! let forcefield = interaction::ForceField::new();
//!
//! // Run simulation
//! // engine::run_md(topo, conf, integrator, forcefield, 1000);
//! ```

pub mod algorithm;
pub mod configuration;
pub mod fep;
pub mod integrator;
pub mod interaction;
pub mod math;
pub mod pairlist;
pub mod selection;
pub mod topology;
// pub mod ffi;  // Temporarily disabled - needs update for refactored API
pub mod eds;
pub mod gamd;
pub mod gpu;
pub mod io;
pub mod logging;
pub mod mpi;
pub mod remd;
pub mod remd_mpi;
pub mod replica;
pub mod validation;

// Re-export main types for convenience
pub use configuration::{Configuration, Energy, State};
pub use math::{BoundaryCondition, Mat3, Periodicity, Vec3};
pub use topology::{LJParameters, Topology};
// pub use interaction::nonbonded::{lj_crf_interaction, ForceStorage};  // ForceStorage removed
pub use eds::{AEDSParameters, EDSForm, EDSParameters, EDSRunner, EDSState};
pub use gamd::{BoostForm, GamdParameters, GamdRunner, GamdStatistics, SearchMode, ThresholdType};
pub use integrator::{Integrator, LeapFrog, SteepestDescent, StochasticDynamics, VelocityVerlet};
pub use remd::{ExchangeScheme, ExchangeStatistics, ExchangeType, ReplicaController};
pub use replica::{Replica, ReplicaId, ReplicaInfo};

#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
