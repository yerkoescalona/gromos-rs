//! SHAKE constraint algorithm wrapper.

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

use crate::constraints::{shake, shake_positions, shake_velocities, ShakeParameters};

/// SHAKE constraint algorithm for the MD sequence.
///
/// Applied after position update to enforce bond length constraints.
/// Equivalent to gromosXX's constraint algorithm in the MD sequence.
pub struct ShakeAlgorithm {
    params: ShakeParameters,
    /// Whether to shake initial positions on init (gromosXX: sim.param().start.shake_pos)
    pub shake_initial_positions: bool,
    /// Whether to shake initial velocities on init (gromosXX: sim.param().start.shake_vel)
    pub shake_initial_velocities: bool,
}

impl ShakeAlgorithm {
    pub fn new(params: ShakeParameters) -> Self {
        Self {
            params,
            shake_initial_positions: false,
            shake_initial_velocities: false,
        }
    }
}

impl Algorithm for ShakeAlgorithm {
    fn init(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        if self.shake_initial_positions {
            log::info!("SHAKE: shaking initial positions");
            let result = shake_positions(topo, conf, sim.dt, &self.params);
            if !result.converged {
                return Err(format!(
                    "SHAKE failed to converge shaking initial positions after {} iterations",
                    result.iterations
                ));
            }

            if self.shake_initial_velocities {
                log::info!("SHAKE: shaking initial velocities");
                let result = shake_velocities(topo, conf, sim.dt, &self.params);
                if !result.converged {
                    return Err(format!(
                        "SHAKE failed to converge shaking initial velocities after {} iterations",
                        result.iterations
                    ));
                }
            }
        }
        Ok(())
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let result = shake(topo, conf, sim.dt, &self.params);
        if result.converged {
            Ok(())
        } else {
            Err(format!(
                "SHAKE failed to converge after {} iterations",
                result.iterations
            ))
        }
    }

    fn name(&self) -> &str {
        "ShakeConstraints"
    }
}
