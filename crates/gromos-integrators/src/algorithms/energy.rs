//! Energy calculation algorithm.
//!
//! Equivalent to gromosXX `algorithm::Energy_Calculation`.
//! Finalizes total energy = E_kin + E_pot.

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

/// Finalizes total energy computation.
///
/// This is the last algorithm in the standard sequence.
/// It ensures E_total = E_kinetic + E_potential is up to date.
#[derive(Debug, Clone)]
pub struct EnergyCalculation;

impl EnergyCalculation {
    pub fn new() -> Self {
        Self
    }
}

impl Default for EnergyCalculation {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for EnergyCalculation {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        // gromosXX convention: energies are stored in old() after exchange_state.
        // Forcefield wrote potential energies to what is now old(),
        // TemperatureCalculation wrote kinetic energy to old().
        // Finalize the potential total here.
        let state = conf.old_mut();
        state.energies.update_potential_total();
        log::debug!("  E_pot={:.10e}  E_kin={:.10e}  E_tot={:.10e}",
            state.energies.potential_total, state.energies.kinetic_total, state.energies.total());
        Ok(())
    }

    fn name(&self) -> &str {
        "Energy_Calculation"
    }
}
