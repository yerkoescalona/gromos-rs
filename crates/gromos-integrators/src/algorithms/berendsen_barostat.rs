//! Berendsen barostat algorithm (weak pressure coupling).
//!
//! Equivalent to GROMOS `algorithm::Berendsen_Barostat`.
//! Isotropic scaling: mu = (1 - comp * dt / tau * (P0 - P))^(1/3)
//! Then scale box and positions by mu.
//!
//! GROMOS sequence position: after PressureCalculation, before EnergyCalculation.
//! Reads pressure_tensor from conf.old(), scales conf.current() positions and box.
//!
//! Source: md++/src/algorithm/pressure/berendsen_barostat.cc

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;

/// Berendsen barostat parameters.
#[derive(Debug, Clone)]
pub struct BerendsenBarostatParams {
    /// Reference pressure P0 (kJ/(mol·nm³)) — GROMOS uses internal units
    pub pressure0: f64,
    /// Isothermal compressibility κ (kJ/(mol·nm³))⁻¹
    pub compressibility: f64,
    /// Coupling time constant τ_P (ps)
    pub tau: f64,
}

/// Berendsen weak-coupling barostat (isotropic).
///
/// Scaling factor: μ = (1 - κ · dt/τ · (P₀ - P))^(1/3)
/// Scale box: box *= μ
/// Scale positions: pos[i] *= μ
#[derive(Debug, Clone)]
pub struct BerendsenBarostat {
    pub params: BerendsenBarostatParams,
}

impl BerendsenBarostat {
    pub fn new(params: BerendsenBarostatParams) -> Self {
        Self { params }
    }
}

impl Algorithm for BerendsenBarostat {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        // Read pressure tensor from PressureCalculation (stored in old())
        let pressure = conf.old().pressure_tensor;

        // Isotropic pressure: P = (P_xx + P_yy + P_zz) / 3
        let total_pressure = (pressure.x_axis.x + pressure.y_axis.y + pressure.z_axis.z) / 3.0;

        // Scaling factor: mu = (1 - comp * dt / tau * (P0 - P))^(1/3)
        let mu = (1.0
            - self.params.compressibility * sim.dt / self.params.tau
                * (self.params.pressure0 - total_pressure))
            .powf(1.0 / 3.0);

        // Scale the box
        {
            let box_vecs = conf.current().box_config.vectors;
            let scaled = gromos_core::math::Mat3::from_cols(
                box_vecs.x_axis * mu,
                box_vecs.y_axis * mu,
                box_vecs.z_axis * mu,
            );
            conf.current_mut().box_config.vectors = scaled;
            conf.current_mut().box_config.inv_vectors = scaled.inverse();
        }

        // Scale the positions
        let n = conf.current().pos.len();
        for i in 0..n {
            conf.current_mut().pos[i] *= mu;
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "Berendsen_Barostat"
    }
}
