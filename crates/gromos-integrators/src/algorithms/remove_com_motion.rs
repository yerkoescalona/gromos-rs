//! Centre of mass motion removal algorithm.
//!
//! Equivalent to gromosXX `algorithm::Remove_COM_Motion`.
//! Removes translational centre of mass motion from velocities.
//!
//! At step 0 (initial): controlled by NTICOM from INITIALISE block.
//! At step > 0 (periodic): every NSCM steps from COMTRANSROT block.
//!
//! Source: md++/src/algorithm/constraints/remove_com_motion.cc

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::math::Vec3;
use gromos_core::topology::Topology;

/// Removes centre of mass translational motion from velocities.
///
/// gromosXX convention: only modifies `conf.current().vel`.
/// The COM velocity is computed as `COM_v = Σ(m_i·v_i) / Σ(m_i)`
/// and subtracted from all atom velocities.
#[derive(Debug, Clone)]
pub struct RemoveCOMMotion {
    /// NTICOM from INITIALISE block: 0=off, 1=remove translation, 2=remove translation+rotation
    pub nticom: i32,
    /// NSCM from COMTRANSROT block: periodic removal every nscm steps (0=off)
    pub nscm: usize,
}

impl RemoveCOMMotion {
    pub fn new(nticom: i32, nscm: usize) -> Self {
        Self { nticom, nscm }
    }
}

/// Calculate and remove translational COM motion.
///
/// Returns the translational COM kinetic energy.
fn remove_com_translation(
    topo: &Topology,
    conf: &mut Configuration,
    remove: bool,
) -> f64 {
    let n_atoms = topo.num_atoms();
    let mut com_v = Vec3::ZERO;
    let mut com_mass = 0.0;

    for i in 0..n_atoms {
        com_mass += topo.mass[i];
        com_v += topo.mass[i] * conf.current().vel[i];
    }

    com_v /= com_mass;
    let ekin_trans = 0.5 * com_mass * com_v.length_squared();

    if remove {
        for i in 0..n_atoms {
            conf.current_mut().vel[i] -= com_v;
        }
        log::debug!(
            "  COM translation removed: v=({:.6e}, {:.6e}, {:.6e}), E_kin_trans={:.6e}",
            com_v.x, com_v.y, com_v.z, ekin_trans
        );
    }

    ekin_trans
}

impl Algorithm for RemoveCOMMotion {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let remove_trans = if sim.step == 0 {
            // Initial removal: controlled by NTICOM
            self.nticom >= 1
        } else {
            // Periodic removal: every NSCM steps
            self.nscm > 0 && (sim.step % self.nscm) == 0
        };

        if !remove_trans {
            return Ok(());
        }

        remove_com_translation(topo, conf, true);

        Ok(())
    }

    fn name(&self) -> &str {
        "Remove_COM_Motion"
    }
}
