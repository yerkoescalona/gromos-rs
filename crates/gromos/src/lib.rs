//! GROMOS-RS: High-performance Rust molecular dynamics engine
//!
//! This is the facade crate that re-exports everything from the sub-crates.
//! For most users, this is the only crate you need to depend on.
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
//! gromos
//! ├── core          - Math types, topology, configuration
//! ├── forces        - Force calculations (bonded + nonbonded)
//! ├── integrators   - Time integration algorithms
//! ├── io            - File I/O
//! └── analysis      - Trajectory analysis
//! ```

// Re-export sub-crates as modules (Polars style)
pub use gromos_core as core;
pub use gromos_forces as forces;
pub use gromos_integrators as integrators;
pub use gromos_io as io;
pub use gromos_analysis as analysis;

// Also re-export with common names for compatibility with existing code
pub mod math {
    pub use gromos_core::math::*;
    pub use gromos_core::{Vec3, Mat3};
}

pub mod topology {
    pub use gromos_core::topology::*;
    pub use gromos_core::Topology;
}

pub mod configuration {
    pub use gromos_core::configuration::*;
    pub use gromos_core::{Configuration, State, Energy};
}

pub mod selection {
    pub use gromos_core::selection::*;
}

pub mod pairlist {
    pub use gromos_core::pairlist::*;
}

// Compatibility aliases for old module names
pub mod algorithm {
    pub use gromos_integrators::constraints;
    pub use gromos_integrators::constraints::*;
    pub use gromos_integrators::thermostats;
    pub use gromos_integrators::thermostats::*;
    pub use gromos_integrators::barostats;
    pub use gromos_integrators::barostats::*;
    pub use gromos_integrators::virtual_atoms;
    pub use gromos_integrators::virtual_atoms::*;
}

pub mod interaction {
    pub use gromos_forces::bonded;
    pub use gromos_forces::bonded::*;
    pub use gromos_forces::nonbonded;
    pub use gromos_forces::nonbonded::*;
    pub use gromos_forces::electrostatics;
    pub use gromos_forces::restraints;
    pub use gromos_forces::pme;
}

pub mod integrator {
    pub use gromos_integrators::integrator::*;
    pub use gromos_integrators::{Integrator, LeapFrog, VelocityVerlet, StochasticDynamics, SteepestDescent};
}

// Enhanced sampling re-exports
pub use gromos_integrators::eds;
pub use gromos_integrators::gamd;
pub use gromos_integrators::remd;
pub use gromos_integrators::replica;
pub use gromos_integrators::fep;

// Logging (re-export the macros)
pub use gromos_core::logging;
pub use gromos_core::{log_debug, log_info, log_warn, log_error};

// Validation
pub use gromos_core::validation;

// Re-export main types at crate root for convenience
pub use gromos_core::{
    Configuration, Energy, State,
    Vec3, Mat3, BoundaryCondition, Periodicity,
    Topology, LJParameters,
};

pub use gromos_integrators::{
    Integrator, LeapFrog, SteepestDescent, StochasticDynamics, VelocityVerlet,
};

#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
