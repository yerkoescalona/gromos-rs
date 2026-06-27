//! SHAKE constraint algorithm wrapper.

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

use crate::constraints::{shake, shake_buffered, shake_positions, shake_velocities, ShakeBuffers, ShakeParameters};

/// SHAKE constraint algorithm for the MD sequence.
///
/// Applied after position update to enforce bond length constraints.
/// Equivalent to GROMOS's constraint algorithm in the MD sequence.
pub struct ShakeAlgorithm {
    params: ShakeParameters,
    /// Whether to shake initial positions on init (GROMOS: sim.param().start.shake_pos)
    pub shake_initial_positions: bool,
    /// Whether to shake initial velocities on init (GROMOS: sim.param().start.shake_vel)
    pub shake_initial_velocities: bool,
    /// Whether SHAKE should also constrain the solvent (false when NTCS selects
    /// a different solvent algorithm, e.g. SETTLE/LINCS)
    pub include_solvent: bool,
    /// Precomputed constraint data and reusable buffers (initialized in init())
    buffers: Option<ShakeBuffers>,
}

impl ShakeAlgorithm {
    pub fn new(params: ShakeParameters) -> Self {
        Self {
            params,
            shake_initial_positions: false,
            shake_initial_velocities: false,
            include_solvent: true,
            buffers: None,
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
        // Precompute constraint lists and allocate reusable buffers
        self.buffers = Some(ShakeBuffers::new(topo, self.params.ntc, self.include_solvent));

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
        let result = if let Some(ref mut buffers) = self.buffers {
            shake_buffered(topo, conf, sim.dt, &self.params, buffers)
        } else {
            // Fallback if init() wasn't called (shouldn't happen in normal flow)
            shake(topo, conf, sim.dt, &self.params)
        };
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
