//! Algorithm implementations for the MD sequence.
//!
//! Each struct implements `gromos_core::algorithm::Algorithm` and can be
//! added to an `AlgorithmSequence` to build the MD step.
//!
//! The standard GROMOS leap-frog MD sequence is:
//! ```text
//!   1. Forcefield (nonbonded + bonded forces)
//!   2. LeapFrogVelocity (exchange_state + velocity update)
//!   3. Thermostat (optional: Berendsen, Nosé-Hoover)
//!   4. LeapFrogPosition (position update)
//!   5. Constraints (optional: SHAKE)
//!   6. TemperatureCalculation (kinetic energy)
//!   7. PressureCalculation (optional)
//!   8. Barostat (optional: Berendsen)
//!   9. EnergyCalculation (total energy)
//! ```

mod leap_frog;
mod temperature;
mod energy;
mod forcefield;
mod shake_algorithm;

pub use leap_frog::{LeapFrogVelocity, LeapFrogPosition};
pub use temperature::TemperatureCalculation;
pub use energy::EnergyCalculation;
pub use forcefield::Forcefield;
pub use shake_algorithm::ShakeAlgorithm;
