//! SETTLE constraint algorithm wrapper.
//!
//! Equivalent to gromosXX `algorithm::Settle` — analytical 3-site rigid water
//! constraint solver (Miyamoto & Kollman, 1992). Solvent-only.
//!
//! Source: md++/src/algorithm/constraints/settle.cc

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

use crate::constraints::settle;

/// SETTLE constraint algorithm for the MD sequence.
///
/// Applied after the position update to analytically enforce the rigid
/// 3-site water geometry. Equivalent to gromosXX's `Settle::apply`.
pub struct SettleAlgorithm;

impl SettleAlgorithm {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SettleAlgorithm {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for SettleAlgorithm {
    fn init(
        &mut self,
        topo: &Topology,
        _conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        // gromosXX Settle::init validation (settle.cc:53-118)
        if topo.num_solvent_molecules() == 0 {
            return Err("SETTLE does only work if 1 solvent.".to_string());
        }
        if topo.atoms_per_solvent() != 3 {
            return Err("SETTLE does only work with water like molecules (3 atoms).".to_string());
        }
        let n_solute = topo.num_solute_atoms();
        if topo.mass[n_solute + 1] != topo.mass[n_solute + 2] {
            return Err(
                "SETTLE does only work with water like molecules (wrong masses).".to_string(),
            );
        }
        if topo.solvent_constraint_template.len() != 3 {
            return Err(
                "SETTLE does only work with water like molecules (3 distance constraints)."
                    .to_string(),
            );
        }
        if (topo.solvent_constraint_template[0].length - topo.solvent_constraint_template[1].length)
            .abs()
            > 1e-12
        {
            return Err(
                "SETTLE does only work with water like molecules (distance constraints wrong)."
                    .to_string(),
            );
        }
        Ok(())
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let result = settle(topo, conf, sim.dt);
        if result.converged {
            Ok(())
        } else {
            Err("SETTLE failed: sin(angle) > 1.0 in analytical solve".to_string())
        }
    }

    fn name(&self) -> &str {
        "SettleConstraints"
    }
}
