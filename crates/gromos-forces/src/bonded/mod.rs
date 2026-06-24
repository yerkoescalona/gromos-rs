//! Bonded interaction calculations.
//!
//! Submodules by term type; this module re-exports everything and provides
//! the shared types and top-level combiners.

pub mod angles;
pub mod bonds;
pub mod dihedrals;
pub mod improper;
pub mod perturbed;

pub use angles::*;
pub use bonds::*;
pub use dihedrals::*;
pub use improper::*;
pub use perturbed::{
    calculate_perturbed_angle_forces,
    calculate_perturbed_bond_forces,
    calculate_perturbed_bonded_forces,
    calculate_perturbed_dihedral_forces,
    calculate_perturbed_improper_dihedral_forces,
};

use gromos_core::math::Vec3;
use gromos_core::topology::Topology;
use gromos_core::configuration::Configuration;

// ─── Shared result types ─────────────────────────────────────────────────────

/// Energy + per-atom forces + virial tensor from a bonded term.
#[derive(Debug, Clone)]
pub struct ForceEnergy {
    pub energy: f64,
    pub forces: Vec<Vec3>,
    /// virial[a][b] = Σ r[b] * f[a]  (gromosXX convention)
    pub virial: [[f64; 3]; 3],
}

/// Energy + per-atom forces + dH/dλ from a perturbed (FEP/TI) bonded term.
#[derive(Debug, Clone)]
pub struct ForceEnergyLambda {
    pub energy: f64,
    pub forces: Vec<Vec3>,
    pub lambda_derivative: f64,
}

impl ForceEnergy {
    pub fn new(num_atoms: usize) -> Self {
        Self { energy: 0.0, forces: vec![Vec3::ZERO; num_atoms], virial: [[0.0; 3]; 3] }
    }

    pub fn add(&mut self, other: &ForceEnergy) {
        self.energy += other.energy;
        for i in 0..self.forces.len().min(other.forces.len()) {
            self.forces[i] += other.forces[i];
        }
        for a in 0..3 {
            for b in 0..3 { self.virial[a][b] += other.virial[a][b]; }
        }
    }
}

impl ForceEnergyLambda {
    pub fn new(num_atoms: usize) -> Self {
        Self { energy: 0.0, forces: vec![Vec3::ZERO; num_atoms], lambda_derivative: 0.0 }
    }

    pub fn add(&mut self, other: &ForceEnergyLambda) {
        self.energy += other.energy;
        self.lambda_derivative += other.lambda_derivative;
        for i in 0..self.forces.len().min(other.forces.len()) {
            self.forces[i] += other.forces[i];
        }
    }
}

// ─── LambdaController (kept for back-compat with existing test code) ─────────

#[derive(Debug, Clone, Default)]
pub struct InteractionLambdas {
    pub bond: f64, pub angle: f64, pub dihedral: f64,
    pub improper: f64, pub lj: f64, pub coulomb: f64,
}

#[derive(Debug, Clone, Default)]
pub struct LambdaController {
    pub lambda: f64,
    pub interaction_lambdas: InteractionLambdas,
}

impl LambdaController {
    pub fn new() -> Self { Self::default() }

    pub fn with_lambda(mut self, lambda: f64) -> Self {
        self.lambda = lambda;
        self.interaction_lambdas = InteractionLambdas {
            bond: lambda, angle: lambda, dihedral: lambda,
            improper: lambda, lj: lambda, coulomb: lambda,
        };
        self
    }

    pub fn get_lambda(&self) -> f64 { self.lambda }
    pub fn lambda_derivative(&self) -> f64 { 1.0 }
}

// ─── Top-level combiners ─────────────────────────────────────────────────────

/// Calculate all bonded forces.
pub fn calculate_bonded_forces(
    topo: &Topology,
    conf: &Configuration,
    use_quartic_bonds: bool,
) -> ForceEnergy {
    calculate_bonded_forces_ntf(topo, conf, use_quartic_bonds, true, true, true, true)
}

/// Calculate bonded forces gated by NTF flags (FORCE block in gromosXX).
pub fn calculate_bonded_forces_ntf(
    topo: &Topology,
    conf: &Configuration,
    use_quartic_bonds: bool,
    ntf_bond: bool,
    ntf_angle: bool,
    ntf_dihedral: bool,
    ntf_improper: bool,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    if ntf_bond {
        let bf = if use_quartic_bonds {
            calculate_bond_forces_quartic(topo, conf)
        } else {
            calculate_bond_forces_harmonic(topo, conf)
        };
        log::debug!("  bonded: bond={:.6e}", bf.energy);
        result.add(&bf);
    }
    if ntf_angle {
        let af = calculate_angle_forces(topo, conf);
        log::debug!("  bonded: angle={:.6e}", af.energy);
        result.add(&af);
    }
    if ntf_dihedral {
        let df = calculate_dihedral_forces(topo, conf);
        log::debug!("  bonded: dihe={:.6e}", df.energy);
        result.add(&df);
    }
    if ntf_improper {
        let imf = calculate_improper_dihedral_forces(topo, conf);
        log::debug!("  bonded: impr={:.6e}", imf.energy);
        result.add(&imf);
    }

    log::debug!("  bonded total={:.6e}  max|f|={:.6e}",
        result.energy,
        result.forces.iter().map(|f| f.length()).fold(0.0_f64, f64::max));
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::configuration::Configuration;
    use gromos_core::topology::Topology;

    #[test]
    fn test_bonded_forces_complete() {
        let topo = Topology::new();
        let conf = Configuration::new(0, 1, 1);
        let result = calculate_bonded_forces(&topo, &conf, true);
        assert_eq!(result.energy, 0.0);
        assert!(result.forces.is_empty());
    }
}
