//! Temperature / kinetic energy calculation algorithm.
//!
//! Equivalent to GROMOS `algorithm::Temperature_Calculation`.
//! Computes kinetic energy using the GROMOS convention:
//! E_kin = 0.5 * sum_i( m_i * (|v_new_i|^2 + |v_old_i|^2) / 2 )
//! where v_new = conf.current().vel, v_old = conf.old().vel
//!
//! **Init behavior**: `init()` calls `apply()` once before the MD loop to
//! pre-compute `kinetic_energy_new` from the initial velocities. This is
//! required so that the Berendsen thermostat has a valid E_kin for its very
//! first scaling step. This matches GROMOS `Temperature_Calculation::init()`
//! (see md++/src/algorithm/temperature/temperature_calculation.cc line ~248-260).
//!
//! Note: this is NOT a GROMOS bug — it is physically correct. The thermostat
//! needs to know the current kinetic energy to compute the first scaling factor.
//! Without this init, E_kin_new=0 at step 0, the thermostat falls back to
//! T_free=T0 and scale=1.0, producing an unscaled first step that cascades
//! into diverging trajectories.
//!
//! Source: md++/src/configuration/state_properties.cc `molecular_translational_ekin()`

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

/// Calculates kinetic energy and temperature using GROMOS's averaged formula.
///
/// In GROMOS, E_kin = 0.5 * sum_i( m_i * (|v_new|^2 + |v_old|^2) / 2 )
/// This averages the kinetic energy between the current and old velocities.
/// The result is stored in conf.old().energies (GROMOS convention).
#[derive(Debug, Clone, Default)]
pub struct TemperatureCalculation {
    /// Last atom (0-based, inclusive) of each temperature bath, in bath order. Empty = one bath
    /// over the whole system, where only the totals are written.
    bath_last_atom: Vec<usize>,
}

impl TemperatureCalculation {
    pub fn new() -> Self {
        Self {
            bath_last_atom: Vec::new(),
        }
    }

    /// The kinetic energy is accumulated per bath as well as in total, as gromosXX does
    /// (`Multibath::calculate_kinetic_energy`): each bath's thermostat scales on its own energy.
    pub fn with_baths(bath_last_atom: &[usize]) -> Self {
        Self {
            bath_last_atom: bath_last_atom.to_vec(),
        }
    }
}

impl Algorithm for TemperatureCalculation {
    /// Initialize by running apply() once to pre-compute kinetic_energy_new.
    ///
    /// GROMOS: Temperature_Calculation::init() calls apply() so that
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
        // GROMOS convention: E_kin = 0.5 * sum_i( m_i * (|v_new|^2 + |v_old|^2) / 2 )
        // After Leap_Frog_Velocity's exchange_state():
        //   current().vel = new velocities v(t+dt/2)
        //   old().vel     = previous velocities v(t-dt/2)
        // The average gives the kinetic energy at time t.
        let n_atoms = topo.inverse_mass.len();
        let n_baths = self.bath_last_atom.len();
        let mut e_kin = 0.0;
        let mut e_kin_new = 0.0;
        let mut per_bath = vec![0.0; n_baths];
        let mut per_bath_new = vec![0.0; n_baths];
        let mut bath = 0;
        for i in 0..n_atoms {
            let v_new = conf.current().vel[i];
            let v_old = conf.old().vel[i];
            let m = topo.mass[i];
            let avg = 0.5 * m * (v_new.length_squared() + v_old.length_squared()) / 2.0;
            let new = 0.5 * m * v_new.length_squared();
            e_kin += avg;
            e_kin_new += new;
            if n_baths > 0 {
                while bath + 1 < n_baths && i > self.bath_last_atom[bath] {
                    bath += 1;
                }
                per_bath[bath] += avg;
                per_bath_new[bath] += new;
            }
        }
        conf.old_mut().energies.kinetic_total = e_kin;
        // Store "new" E_kin for thermostat scaling (GROMOS: multibath.bath.ekin)
        conf.old_mut().energies.kinetic_energy_new = e_kin_new;
        if n_baths > 0 {
            conf.old_mut().energies.kinetic_energy = per_bath;
            conf.old_mut().energies.kinetic_energy_new_bath = per_bath_new;
        }
        log::debug!("  E_kin={:.10e}  E_kin_new={:.10e}", e_kin, e_kin_new);
        Ok(())
    }

    fn name(&self) -> &str {
        "Temperature_Calculation"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::math::Vec3;

    #[test]
    fn kinetic_energy_is_the_average_of_the_old_and_new_half_steps() {
        // gromosXX: E_kin = ½ m (v_old² + v_new²) / 2 — the leap-frog average
        let mut topo = Topology::new();
        topo.mass = vec![2.0];
        topo.inverse_mass = vec![0.5];
        let mut conf = Configuration::new(1, 1, 1);
        conf.old_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().vel[0] = Vec3::new(3.0, 0.0, 0.0);
        let sim = SimulationState::new(0.002, 1);
        TemperatureCalculation::new()
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        // ½·2·(1 + 9)/2 = 5 ; new-only: ½·2·9 = 9
        assert_eq!(conf.old().energies.kinetic_total, 5.0);
        assert_eq!(conf.old().energies.kinetic_energy_new, 9.0);
    }

    #[test]
    fn init_computes_the_same_thing_as_apply() {
        let mut topo = Topology::new();
        topo.mass = vec![1.0, 1.0];
        topo.inverse_mass = vec![1.0, 1.0];
        let mut conf = Configuration::new(2, 1, 1);
        conf.old_mut().vel[1] = Vec3::new(0.0, 2.0, 0.0);
        conf.current_mut().vel[1] = Vec3::new(0.0, 2.0, 0.0);
        let sim = SimulationState::new(0.002, 1);
        TemperatureCalculation::new()
            .init(&topo, &mut conf, &sim)
            .unwrap();
        assert_eq!(conf.old().energies.kinetic_total, 2.0);
    }
}
