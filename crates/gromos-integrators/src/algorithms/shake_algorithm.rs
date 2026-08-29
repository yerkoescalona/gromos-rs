//! SHAKE constraint algorithm wrapper.

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

use crate::constraints::{
    shake, shake_buffered, shake_positions, shake_velocities, ShakeBuffers, ShakeParameters,
};

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
        self.buffers = Some(ShakeBuffers::new(
            topo,
            self.params.ntc,
            self.include_solvent,
        ));

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::NtcMode;
    use gromos_core::math::Vec3;
    use gromos_core::topology::{Bond, BondParameters, MolTypeAtom};

    fn two_atom_bond_topology() -> Topology {
        let mut topo = Topology::new();
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0];
        topo.moltypes[0].atoms = vec![
            MolTypeAtom {
                name: "C1".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
            MolTypeAtom {
                name: "C2".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
        ];
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 0.0,
            r0: 0.15,
        });
        topo
    }

    /// Old positions at the constraint length, new ones stretched by 0.01 nm.
    fn stretched_pair() -> Configuration {
        let mut conf = Configuration::new(2, 1, 1);
        conf.old_mut().pos[1] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.16, 0.0, 0.0);
        conf
    }

    fn params() -> ShakeParameters {
        ShakeParameters {
            tolerance: 1e-8,
            max_iterations: 1000,
            ntc: NtcMode::AllBonds,
        }
    }

    #[test]
    fn apply_restores_the_bond_length() {
        let topo = two_atom_bond_topology();
        let mut conf = stretched_pair();
        let sim = SimulationState::new(0.002, 10);
        let mut alg = ShakeAlgorithm::new(params());
        alg.init(&topo, &mut conf, &sim).unwrap();
        alg.apply(&topo, &mut conf, &sim).unwrap();
        let d = (conf.current().pos[1] - conf.current().pos[0]).length();
        assert!((d - 0.15).abs() < 1e-7, "distance after SHAKE {d}");
        // equal masses: both atoms move towards each other by the same amount
        let (p0, p1) = (conf.current().pos[0], conf.current().pos[1]);
        assert!(
            (p0.x - 0.005).abs() < 1e-7 && (p1.x - 0.155).abs() < 1e-7,
            "p0 {p0:?} p1 {p1:?}"
        );
        assert_eq!(alg.name(), "ShakeConstraints");
    }

    #[test]
    fn init_can_shake_the_initial_positions() {
        let topo = two_atom_bond_topology();
        let mut conf = stretched_pair();
        conf.old_mut().pos[1] = Vec3::new(0.16, 0.0, 0.0);
        let sim = SimulationState::new(0.002, 10);
        assert_eq!(topo.num_atoms(), 2);
        let mut alg = ShakeAlgorithm::new(params());
        alg.shake_initial_positions = true;
        alg.init(&topo, &mut conf, &sim).unwrap();
        let d = (conf.current().pos[1] - conf.current().pos[0]).length();
        assert!((d - 0.15).abs() < 1e-7, "initial positions not shaken: {d}");
    }

    #[test]
    fn solvent_only_mode_leaves_a_solute_bond_alone() {
        let topo = two_atom_bond_topology();
        let mut conf = stretched_pair();
        let sim = SimulationState::new(0.002, 10);
        let mut alg = ShakeAlgorithm::new(ShakeParameters {
            ntc: NtcMode::SolventOnly,
            ..params()
        });
        alg.init(&topo, &mut conf, &sim).unwrap();
        alg.apply(&topo, &mut conf, &sim).unwrap();
        assert_eq!(conf.current().pos[1].x, 0.16);
    }
}
