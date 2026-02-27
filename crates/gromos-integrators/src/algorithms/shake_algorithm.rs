//! SHAKE constraint algorithm wrapper.

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

use crate::constraints::{shake, ShakeParameters};

/// SHAKE constraint algorithm for the MD sequence.
///
/// Applied after position update to enforce bond length constraints.
/// Equivalent to gromosXX's constraint algorithm in the MD sequence.
pub struct ShakeAlgorithm {
    params: ShakeParameters,
}

impl ShakeAlgorithm {
    pub fn new(params: ShakeParameters) -> Self {
        Self { params }
    }
}

impl Algorithm for ShakeAlgorithm {
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
