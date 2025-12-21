//! GROMOS-RS Forces: All force field calculations
//!
//! This crate provides force calculation routines:
//! - Bonded forces (bonds, angles, dihedrals, impropers)
//! - Nonbonded forces (Lennard-Jones, Coulomb)
//! - Electrostatics (PME, reaction field)
//! - Restraints (position, distance, dihedral, NMR)
//! - Enhanced sampling (local elevation, QM/MM, polarization)
//!
//! All calculations use f64 for molecular dynamics precision.

pub mod bonded;
pub mod nonbonded;
pub mod electrostatics;
pub mod restraints;
pub mod pme;

// Advanced force modules
pub mod local_elevation;
pub mod polarization;
pub mod qmmm;

// GPU acceleration (feature-gated)
#[cfg(feature = "gpu")]
pub mod gpu;

// MPI-enabled PME (feature-gated)
#[cfg(feature = "mpi")]
pub mod pme_mpi;

// Re-export main types
pub use bonded::{
    ForceEnergy, ForceEnergyLambda, LambdaController,
    calculate_bonded_forces, calculate_angle_forces, 
    calculate_dihedral_forces, calculate_improper_dihedral_forces,
};
pub use nonbonded::{
    LJParameters, CRFParameters, ForceStorage, lj_crf_interaction,
};
