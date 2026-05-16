//! Berendsen weak-coupling thermostat algorithm.
//!
//! Equivalent to gromosXX `algorithm::Berendsen_Thermostat`.
//! Rescales velocities to weakly couple the system to a heat bath.
//!
//! Placed between LeapFrogVelocity and LeapFrogPosition in the sequence.
//! Uses the "new" kinetic energy from the previous step's TemperatureCalculation
//! (stored in energies.kinetic_energy_new, equivalent to gromosXX multibath.bath.ekin).
//!
//! Source: md++/src/algorithm/temperature/berendsen_thermostat.cc
//!         md++/src/algorithm/temperature/thermostat.cc (scale method)

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

/// Boltzmann constant in kJ/(mol·K)
const K_BOLTZMANN: f64 = 0.00831441;

/// Berendsen thermostat parameters for a single temperature bath.
#[derive(Debug, Clone)]
pub struct BerendsenThermostatParams {
    /// Reference temperature T0 (K)
    pub temperature: f64,
    /// Coupling time constant τ (ps). τ < 0 → no coupling, τ = 0 → instantaneous
    pub tau: f64,
    /// Degrees of freedom for this bath
    pub dof: f64,
}

/// Berendsen weak-coupling thermostat.
///
/// gromosXX sequence position: after Leap_Frog_Velocity, before Leap_Frog_Position.
///
/// The scaling factor is:
///   λ = sqrt(1 + dt/τ · (T₀/T_free - 1))
/// where T_free = 2·E_kin / (DOF · k_B) using E_kin from the previous step.
#[derive(Debug, Clone)]
pub struct BerendsenThermostat {
    pub params: Vec<BerendsenThermostatParams>,
    /// Atom range per bath: (first_atom, last_atom) inclusive
    pub bath_ranges: Vec<(usize, usize)>,
}

impl BerendsenThermostat {
    /// Create a new Berendsen thermostat with a single bath covering all atoms.
    pub fn new_single_bath(temperature: f64, tau: f64, dof: f64, n_atoms: usize) -> Self {
        Self {
            params: vec![BerendsenThermostatParams {
                temperature,
                tau,
                dof,
            }],
            bath_ranges: vec![(0, n_atoms - 1)],
        }
    }
}

impl Algorithm for BerendsenThermostat {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        // After exchange_state in LeapFrogVelocity:
        //   current().energies.kinetic_energy_new = E_kin_new from previous step
        // (set by TemperatureCalculation into old().energies, then swapped by exchange_state)

        for (bath_idx, params) in self.params.iter().enumerate() {
            // Skip if no coupling (tau < 0)
            if params.tau < 0.0 {
                continue;
            }

            // Read E_kin from previous step (stored by TemperatureCalculation)
            let ekin = conf.current().energies.kinetic_energy_new;

            // Compute free temperature: T_free = 2·E_kin / (DOF · k_B)
            let mut free_temp = if params.dof > 0.0 {
                2.0 * ekin / (params.dof * K_BOLTZMANN)
            } else {
                0.0
            };

            // Fallback: if T_free ≈ 0, use reference temperature (gromosXX convention)
            if free_temp < 1e-10 {
                free_temp = params.temperature;
            }

            // Second fallback: if reference T is also 0, no scaling
            let scale = if free_temp < 1e-10 {
                1.0
            } else if params.tau == 0.0 {
                // Instantaneous scaling (tau = 0)
                (params.temperature / free_temp).sqrt()
            } else {
                // Weak coupling: λ = sqrt(1 + dt/τ · (T₀/T - 1))
                (1.0 + sim.dt / params.tau * (params.temperature / free_temp - 1.0)).sqrt()
            };

            log::debug!(
                "  Berendsen bath {}: E_kin_new={:.6e}, T_free={:.2}, T0={:.2}, scale={:.10}",
                bath_idx, ekin, free_temp, params.temperature, scale
            );

            if (scale - 1.0).abs() < 1e-15 {
                continue; // No scaling needed
            }

            // Scale velocities for atoms in this bath range
            let (first, last) = self.bath_ranges[bath_idx];
            for i in first..=last {
                conf.current_mut().vel[i] *= scale;
            }
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "Berendsen_Thermostat"
    }
}
