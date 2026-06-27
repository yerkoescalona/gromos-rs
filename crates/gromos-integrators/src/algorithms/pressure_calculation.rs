//! Pressure tensor calculation algorithm.
//!
//! Equivalent to GROMOS `algorithm::Pressure_Calculation`.
//! Computes the pressure tensor from the kinetic energy tensor and virial tensor:
//!   P_tensor = (KE_tensor + 0.5 * virial_tensor) * (2.0 / V)
//!
//! The virial is stored as the raw outer product r⊗f (without -0.5 prefactor).
//! The 0.5 factor is applied here.
//!
//! KE tensor and virial correction are now computed inside the Forcefield
//! (prepare_virial + atomic_to_molecular_virial), matching GROMOS.
//! After exchange_state, they are in conf.old().
//!
//! Source: md++/src/algorithm/pressure/pressure_calculation.cc

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::math::Mat3;
use gromos_core::topology::Topology;

/// Virial type: how the virial and KE tensor are computed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VirialType {
    None,
    Atomic,
    Molecular,
}

/// Calculates the pressure tensor from KE tensor and virial tensor.
///
/// GROMOS sequence position: after TemperatureCalculation, before BerendsenBarostat.
/// Reads from conf.old() (after exchange_state).
///
/// The KE tensor and molecular virial correction have already been computed
/// by the Forcefield algorithm (on conf.current(), which is now conf.old()).
///
/// Formula: P_tensor = (KE_tensor + 0.5 * virial_tensor) * (2.0 / V)
#[derive(Debug, Clone)]
pub struct PressureCalculation {
    pub virial_type: VirialType,
}

impl PressureCalculation {
    pub fn new(virial_type: VirialType) -> Self {
        Self { virial_type }
    }
}

impl Algorithm for PressureCalculation {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        if self.virial_type == VirialType::None {
            return Ok(());
        }

        let volume = conf.old().box_config.volume();

        if volume < 1e-30 {
            return Ok(());
        }

        // KE tensor and virial tensor are already in conf.old()
        // (computed by Forcefield on current(), then swapped by exchange_state)
        let ke_tensor = conf.old().kinetic_energy_tensor;
        let virial_tensor = conf.old().virial_tensor;

        // P_tensor = (KE_tensor + 0.5 * virial_tensor) * (2.0 / V)
        // tot_cg_factor = 1.0 for non-coarse-grained systems
        let factor = 2.0 / volume;
        let half_virial = Mat3::from_cols(
            virial_tensor.x_axis * 0.5,
            virial_tensor.y_axis * 0.5,
            virial_tensor.z_axis * 0.5,
        );
        let sum = Mat3::from_cols(
            ke_tensor.x_axis + half_virial.x_axis,
            ke_tensor.y_axis + half_virial.y_axis,
            ke_tensor.z_axis + half_virial.z_axis,
        );
        let p_tensor = Mat3::from_cols(
            sum.x_axis * factor,
            sum.y_axis * factor,
            sum.z_axis * factor,
        );

        // Store pressure tensor in old() (GROMOS convention)
        conf.old_mut().pressure_tensor = p_tensor;

        Ok(())
    }

    fn name(&self) -> &str {
        "Pressure_Calculation"
    }
}
