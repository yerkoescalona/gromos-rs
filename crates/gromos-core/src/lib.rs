//! GROMOS-RS Core: Math types, topology, and fundamental data structures
//!
//! This crate provides the foundational types used throughout GROMOS-RS:
//! - Floating-point precision trait (`Float`)
//! - Generic vector and matrix types (`Vector3<F>`, `Matrix3<F>`)
//! - Boundary conditions and periodicity
//! - Atom and topology definitions
//! - System configuration and state
//! - Pairlist generation
//! - Selection language
//!
//! # Precision
//!
//! All types are generic over the floating-point precision. Use `f64` for
//! CPU-based MD simulations (better energy conservation) or `f32` for
//! GPU acceleration and reduced memory bandwidth.
//!
//! ```rust
//! use gromos_core::{Vector3, Vec3f, Vec3d};
//!
//! // Type aliases for convenience
//! let pos_f32: Vec3f = Vector3::new(1.0, 2.0, 3.0);
//! let pos_f64: Vec3d = Vector3::new(1.0, 2.0, 3.0);
//!
//! // Convert between precisions
//! let converted = pos_f32.to_f64();
//! ```

pub mod float;
pub mod vector;
pub mod matrix;
pub mod math;
pub mod topology;
pub mod configuration;
pub mod pairlist;
pub mod selection;
pub mod logging;
pub mod validation;

// FFI for C bindings (feature-gated)
#[cfg(feature = "ffi")]
pub mod ffi;

// Re-export Float trait and precision types
pub use float::{Float, F32, F64};

// Re-export generic vector and matrix types
pub use vector::{Vector3, Vec3f, Vec3d, VEC3F_ZERO, VEC3D_ZERO};
pub use matrix::{Matrix3, Mat3f, Mat3d};

// Legacy compatibility - keep the old Vec3/Mat3 names pointing to the glam-based types
// These will be phased out in favor of the generic types
pub use math::{BoundaryCondition, Mat3, Periodicity, Rectangular, Triclinic, Vacuum, Vec3};

// Re-export configuration types
pub use configuration::{Box, BoxType, Configuration, Energy, State, StochasticVariables};

// Re-export topology types
pub use topology::{LJParameters, Topology};

// Re-export pairlist types
pub use pairlist::{Pairlist, PairlistContainer, StandardPairlistAlgorithm};

// Re-export selection
pub use selection::AtomSelection;

