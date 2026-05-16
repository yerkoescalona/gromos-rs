//! Temperature / kinetic energy calculation algorithm.
//!
//! Equivalent to gromosXX `algorithm::Temperature_Calculation`.
//! Computes kinetic energy using the gromosXX convention:
//! E_kin = 0.5 * sum_i( m_i * (|v_new_i|^2 + |v_old_i|^2) / 2 )
//! where v_new = conf.current().vel, v_old = conf.old().vel
//!
//! **Init behavior**: `init()` calls `apply()` once before the MD loop to
//! pre-compute `kinetic_energy_new` from the initial velocities. This is
//! required so that the Berendsen thermostat has a valid E_kin for its very
//! first scaling step. This matches gromosXX `Temperature_Calculation::init()`
//! (see md++/src/algorithm/temperature/temperature_calculation.cc line ~248-260).
//!
//! Note: this is NOT a gromosXX bug — it is physically correct. The thermostat
//! needs to know the current kinetic energy to compute the first scaling factor.
//! Without this init, E_kin_new=0 at step 0, the thermostat falls back to
//! T_free=T0 and scale=1.0, producing an unscaled first step that cascades
//! into diverging trajectories.
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
    /// Initialize by running apply() once to pre-compute kinetic_energy_new.
    ///
    /// gromosXX: Temperature_Calculation::init() calls apply() so that
    /// multibath.bath.ekin is populated before the first MD step. Without
    /// this, the Berendsen thermostat would see E_kin=0 at step 0 and
    /// skip scaling entirely.
    fn init(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        log::debug!("TemperatureCalculation::init() — pre-computing E_kin from initial velocities");
        self.apply(topo, conf, sim)
    }

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
        let mut e_kin_new = 0.0;
        for i in 0..n_atoms {
            let v_new = conf.current().vel[i];
            let v_old = conf.old().vel[i];
            let m = topo.mass[i];
            e_kin += 0.5 * m * (v_new.length_squared() + v_old.length_squared()) / 2.0;
            e_kin_new += 0.5 * m * v_new.length_squared();
        }
        conf.old_mut().energies.kinetic_total = e_kin;
        // Store "new" E_kin for thermostat scaling (gromosXX: multibath.bath.ekin)
        conf.old_mut().energies.kinetic_energy_new = e_kin_new;
        log::debug!("  E_kin={:.10e}  E_kin_new={:.10e}", e_kin, e_kin_new);
        Ok(())
    }

    fn name(&self) -> &str {
        "Temperature_Calculation"
    }
}
