//! LINCS constraint algorithm wrapper.
//!
//! Equivalent to GROMOS `algorithm::Lincs` — linear constraint solver
//! (Hess et al. 1997), selectable for the solute (NTCP) and/or solvent (NTCS).
//!
//! Source: md++/src/algorithm/constraints/lincs.cc

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

use crate::constraints::{lincs_buffered, LincsBuffers, NtcMode};

/// LINCS constraint algorithm for the MD sequence.
///
/// `init()` precomputes the coupling matrices for the solute and/or solvent
/// (whichever is selected via `include_solute`/`include_solvent`), `apply()`
/// runs LINCS for both groups every step (GROMOS `Lincs::apply`,
/// `lincs.cc:238-284`).
pub struct LincsAlgorithm {
    ntc: NtcMode,
    solute_order: usize,
    solvent_order: usize,
    include_solute: bool,
    include_solvent: bool,
    buffers: Option<LincsBuffers>,
}

impl LincsAlgorithm {
    pub fn new(
        ntc: NtcMode,
        solute_order: usize,
        solvent_order: usize,
        include_solute: bool,
        include_solvent: bool,
    ) -> Self {
        Self {
            ntc,
            solute_order,
            solvent_order,
            include_solute,
            include_solvent,
            buffers: None,
        }
    }
}

impl Algorithm for LincsAlgorithm {
    fn init(
        &mut self,
        topo: &Topology,
        _conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        self.buffers = Some(LincsBuffers::new(
            topo,
            self.ntc,
            self.solute_order,
            self.solvent_order,
            self.include_solute,
            self.include_solvent,
        ));
        Ok(())
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        if let Some(ref buffers) = self.buffers {
            // GROMOS `Lincs::apply` (lincs.cc:238-288) always returns success —
            // LINCS is a fixed-order analytical expansion, not an iterative solver
            // that can fail to converge, so we don't gate on `max_error` here
            // (unlike SHAKE, whose iteration can genuinely fail).
            let _ = lincs_buffered(topo, conf, sim.dt, buffers);
            Ok(())
        } else {
            Err("LincsAlgorithm::apply called before init".to_string())
        }
    }

    fn name(&self) -> &str {
        "LincsConstraints"
    }
}
