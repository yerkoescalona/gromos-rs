//! GROMOS-RS Integrators: Time integration algorithms
//!
//! This crate provides integration algorithms:
//! - Velocity Verlet
//! - Leap-frog
//! - Stochastic dynamics (Langevin)
//! - Steepest descent minimization
//! - Constraints (SHAKE, LINCS)
//! - Thermostats and barostats
//! - Virtual atoms
//! - Enhanced sampling: REMD, EDS, GAMD, FEP

pub mod integrator;
pub mod constraints;
pub mod thermostats;
pub mod barostats;
pub mod virtual_atoms;

// Enhanced sampling modules
pub mod eds;
pub mod fep;
pub mod gamd;
pub mod remd;
pub mod replica;

// MPI support (optional)
#[cfg(feature = "mpi")]
pub mod mpi;
#[cfg(feature = "mpi")]
pub mod remd_mpi;

// Re-export main types
pub use integrator::{Integrator, LeapFrog, VelocityVerlet, StochasticDynamics, SteepestDescent};
pub use constraints::{ShakeParameters, LincsParameters, ConstraintResult, shake, lincs, settle};
pub use thermostats::{BerendsenThermostatParameters, NoseHooverThermostatParameters, AndersenThermostatParameters};
pub use thermostats::{berendsen_thermostat, nose_hoover_thermostat, andersen_thermostat};
pub use barostats::{BerendsenBarostatParameters, ParrinelloRahmanBarostatParameters};
pub use barostats::{berendsen_barostat, parrinello_rahman_barostat, calculate_virial};
pub use virtual_atoms::{VirtualAtom, VirtualAtomsManager};
