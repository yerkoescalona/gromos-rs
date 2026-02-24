//! Temperature / kinetic energy calculation algorithm.
//!
//! Equivalent to gromosXX `algorithm::Temperature_Calculation`.
//! Computes kinetic energy using the gromosXX convention:
//! E_kin = 0.5 * sum_i( m_i * (|v_new_i|^2 + |v_old_i|^2) / 2 )
//! where v_new = conf.current().vel, v_old = conf.old().vel
//!
//! Source: md++/src/configuration/state_properties.cc `molecular_translational_ekin()`

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

/// Calculates kinetic energy and temperature using gromosXX's averaged formula.
///
/// In gromosXX, E_kin = 0.5 * sum_i( m_i * (|v_new|^2 + |v_old|^2) / 2 )
/// This averages the kinetic energy between the current and old velocities.
/// The result is stored in conf.old().energies (gromosXX convention).
#[derive(Debug, Clone)]
pub struct TemperatureCalculation;

impl TemperatureCalculation {
    pub fn new() -> Self {
        Self
    }
}

impl Default for TemperatureCalculation {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for TemperatureCalculation {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        // gromosXX convention: E_kin = 0.5 * sum_i( m_i * (|v_new|^2 + |v_old|^2) / 2 )
        // After Leap_Frog_Velocity's exchange_state():
        //   current().vel = new velocities v(t+dt/2)
        //   old().vel     = previous velocities v(t-dt/2)
        // The average gives the kinetic energy at time t.
        let n_atoms = topo.inverse_mass.len();
        let mut e_kin = 0.0;
        for i in 0..n_atoms {
            let v_new = conf.current().vel[i];
            let v_old = conf.old().vel[i];
            let m = topo.mass[i];
            e_kin += 0.5 * m * (v_new.length_squared() + v_old.length_squared()) / 2.0;
        }
        conf.old_mut().energies.kinetic_total = e_kin;
        Ok(())
    }

    fn name(&self) -> &str {
        "Temperature_Calculation"
    }
}
