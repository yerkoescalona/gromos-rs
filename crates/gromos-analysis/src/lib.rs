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
pub use gromos_core::{configuration, math, topology, Configuration, Mat3, State, Topology, Vec3};

pub mod diffusion;
pub mod fit;
pub mod gyration;
pub mod hbond;
pub mod rdf;
pub mod rmsd;

// Re-export main functions
pub use diffusion::{calculate_diffusion_coefficient, calculate_msd};
pub use gyration::calculate_gyration_radius;
pub use rdf::calculate_rdf;
pub use rmsd::calculate_rmsd;
