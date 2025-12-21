//! GROMOS-RS Analysis: Trajectory analysis tools
//!
//! This crate provides analysis functions:
//! - RDF (radial distribution function)
//! - RMSD, RMSF
//! - Hydrogen bonds
//! - Gyration radius
//! - Diffusion coefficients
//! - Clustering
//! - And more...

// Re-export gromos-core types for convenience within this crate
pub use gromos_core::{
    math, topology, configuration,
    Vec3, Mat3,
    Configuration, State,
    Topology,
};

pub mod rdf;
pub mod rmsd;
pub mod hbond;
pub mod gyration;
pub mod diffusion;

// Re-export main functions
pub use rdf::calculate_rdf;
pub use rmsd::calculate_rmsd;
pub use gyration::calculate_gyration_radius;
pub use diffusion::{calculate_msd, calculate_diffusion_coefficient};

