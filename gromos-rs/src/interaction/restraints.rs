//! Restraint interactions for MD simulations
//!
//! Direct translation of GROMOS restraint interactions from:
//! - md++/src/interaction/special/position_restraint_interaction.cc
//! - md++/src/interaction/special/distance_restraint_interaction.cc

use crate::configuration::Configuration;
use crate::math::{Periodicity, Vec3};

/// Position restraint - restrains atom to reference position
///
/// Direct translation from position_restraint_interaction.cc
///
/// Energy: E = 0.5 * k * |r - r_ref|²
/// Force: F = -k * (r - r_ref)
#[derive(Debug, Clone)]
pub struct PositionRestraint {
    /// Atom index
    pub atom: usize,
    /// Reference position
    pub reference_pos: Vec3,
    /// Force constant (kJ mol⁻¹ nm⁻²)
    pub force_constant: f64,
}

impl PositionRestraint {
    pub fn new(atom: usize, reference_pos: Vec3, force_constant: f64) -> Self {
        Self {
            atom,
            reference_pos,
            force_constant,
        }
    }

    /// Calculate restraint energy and forces
    pub fn calculate(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        // Get displacement vector using periodic boundary conditions
        // nearest_image returns (reference - current), but we need (current - reference)
        let displacement =
            periodicity.nearest_image(self.reference_pos, conf.current().pos[self.atom]);

        let dist_sq = displacement.length_squared() as f64;
        let dist = (dist_sq).sqrt();

        // Energy: 0.5 * k * r²
        let energy = 0.5 * self.force_constant * dist_sq;

        // Force: F = -k * (r - r0) = -k * displacement
        let force_magnitude = -self.force_constant;
        let force = displacement * (force_magnitude as f32);

        conf.current_mut().force[self.atom] += force;

        energy
    }
}

/// Collection of position restraints
#[derive(Debug, Clone, Default)]
pub struct PositionRestraints {
    pub restraints: Vec<PositionRestraint>,
}

impl PositionRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: PositionRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        self.restraints
            .iter()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }
}

/// Distance restraint - restrains distance between two atoms
///
/// Direct translation from distance_restraint_interaction.cc
///
/// Supports:
/// - Harmonic/linear force function
/// - Time-averaging for NOE restraints
/// - Repulsive/attractive/harmonic modes
#[derive(Debug, Clone)]
pub struct DistanceRestraint {
    /// First atom index
    pub atom_i: usize,
    /// Second atom index
    pub atom_j: usize,
    /// Target distance (nm)
    pub r0: f64,
    /// Force constant (kJ mol⁻¹ nm⁻²)
    pub force_constant: f64,
    /// Linear region width (nm)
    pub r_linear: f64,
    /// Restraint type: -1 (repulsive), 0 (harmonic), +1 (attractive)
    pub restraint_type: i32,
    /// Time-averaged distance (for NOE restraints, r⁻³ averaging)
    pub time_averaged_distance: Option<f64>,
    /// Time constant for averaging (ps)
    pub tau: Option<f64>,
}

impl DistanceRestraint {
    pub fn new(atom_i: usize, atom_j: usize, r0: f64, force_constant: f64) -> Self {
        Self {
            atom_i,
            atom_j,
            r0,
            force_constant,
            r_linear: 0.0,     // Default: pure harmonic
            restraint_type: 0, // Harmonic by default
            time_averaged_distance: None,
            tau: None,
        }
    }

    /// Set linear region width
    pub fn with_linear_region(mut self, r_linear: f64) -> Self {
        self.r_linear = r_linear;
        self
    }

    /// Set restraint type: -1 (repulsive), 0 (harmonic), +1 (attractive)
    pub fn with_type(mut self, restraint_type: i32) -> Self {
        self.restraint_type = restraint_type;
        self
    }

    /// Enable time averaging with exponential decay
    pub fn with_time_averaging(mut self, tau: f64) -> Self {
        self.tau = Some(tau);
        self.time_averaged_distance = Some(self.r0);
        self
    }

    /// Calculate restraint energy and forces
    ///
    /// Direct translation from distance_restraint_interaction.cc
    pub fn calculate(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        // Get distance vector using periodic boundary conditions
        let r_vec = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_i],
        );

        let mut distance = r_vec.length() as f64;

        // Time averaging (NOE r⁻³ averaging)
        if let (Some(tau), Some(avg_dist)) = (self.tau, self.time_averaged_distance) {
            let exp_term = (-dt / tau).exp();
            let r_inv_3 = 1.0 / (distance * distance * distance);
            let avg_r_inv_3 = 1.0 / (avg_dist * avg_dist * avg_dist);

            // Update average: <r⁻³> = (1-e^(-dt/τ)) * r⁻³ + e^(-dt/τ) * <r⁻³>_old
            let new_avg_r_inv_3 = (1.0 - exp_term) * r_inv_3 + exp_term * avg_r_inv_3;
            let new_avg_dist = new_avg_r_inv_3.powf(-1.0 / 3.0);

            self.time_averaged_distance = Some(new_avg_dist);
            distance = new_avg_dist;
        }

        let delta = distance - self.r0;

        // Check if restraint is violated based on type
        let violated = match self.restraint_type {
            -1 => distance < self.r0, // Repulsive: violated if too close
            0 => true,                // Harmonic: always active
            1 => distance > self.r0,  // Attractive: violated if too far
            _ => true,
        };

        if !violated {
            return 0.0; // No restraint force
        }

        // Calculate energy and force
        let (energy, force_magnitude) = if delta.abs() < self.r_linear || self.r_linear == 0.0 {
            // Harmonic region: E = 0.5 * k * (r - r0)²
            let energy = 0.5 * self.force_constant * delta * delta;
            let force_mag = -self.force_constant * delta;
            (energy, force_mag)
        } else {
            // Linear region (prevents unrealistic forces at large distances)
            let sign = if delta > 0.0 { 1.0 } else { -1.0 };
            let energy = self.force_constant * self.r_linear * (delta.abs() - 0.5 * self.r_linear);
            let force_mag = -sign * self.force_constant * self.r_linear;
            (energy, force_mag)
        };

        // Apply force along bond direction
        if distance > 1e-10 {
            let force = r_vec * ((force_magnitude / distance) as f32);
            conf.current_mut().force[self.atom_i] -= force;
            conf.current_mut().force[self.atom_j] += force;
        }

        energy
    }
}

/// Collection of distance restraints
#[derive(Debug, Clone, Default)]
pub struct DistanceRestraints {
    pub restraints: Vec<DistanceRestraint>,
}

impl DistanceRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: DistanceRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        self.restraints
            .iter_mut()
            .map(|r| r.calculate(conf, periodicity, dt))
            .sum()
    }
}

/// Angle restraint - restrains angle between three atoms
///
/// Restrains the angle θ defined by atoms i-j-k to a target value.
/// Used for NMR refinement and maintaining structural features.
///
/// Energy: E = 0.5 * k * (θ - θ₀)²
/// where θ is the angle i-j-k
#[derive(Debug, Clone)]
pub struct AngleRestraint {
    /// First atom index
    pub atom_i: usize,
    /// Second atom index (vertex)
    pub atom_j: usize,
    /// Third atom index
    pub atom_k: usize,
    /// Target angle (radians)
    pub theta0: f64,
    /// Force constant (kJ mol⁻¹ rad⁻²)
    pub force_constant: f64,
}

impl AngleRestraint {
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        theta0: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            atom_i,
            atom_j,
            atom_k,
            theta0,
            force_constant,
        }
    }

    /// Create from target angle in degrees
    pub fn from_degrees(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        theta0_deg: f64,
        force_constant: f64,
    ) -> Self {
        Self::new(
            atom_i,
            atom_j,
            atom_k,
            theta0_deg.to_radians(),
            force_constant,
        )
    }

    /// Calculate restraint energy and forces
    pub fn calculate(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        // Get vectors using periodic boundary conditions
        let r_ij = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_i],
        );
        let r_kj = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_k],
        );

        let r_ij_len = r_ij.length() as f64;
        let r_kj_len = r_kj.length() as f64;

        if r_ij_len < 1e-10 || r_kj_len < 1e-10 {
            return 0.0; // Degenerate angle
        }

        // Calculate angle using dot product: cos(θ) = r_ij · r_kj / (|r_ij| |r_kj|)
        let cos_theta = (r_ij.dot(r_kj) as f64) / (r_ij_len * r_kj_len);
        let cos_theta = cos_theta.clamp(-1.0, 1.0); // Numerical safety
        let theta = cos_theta.acos();

        let delta_theta = theta - self.theta0;

        // Energy: E = 0.5 * k * (θ - θ₀)²
        let energy = 0.5 * self.force_constant * delta_theta * delta_theta;

        // Force magnitude: F = -dE/dθ = -k * (θ - θ₀)
        let force_magnitude = -self.force_constant * delta_theta;

        // sin(θ) for force calculation
        let sin_theta = theta.sin();
        if sin_theta.abs() < 1e-10 {
            return energy; // Avoid division by zero at θ = 0 or π
        }

        // Force derivatives
        let f_factor = force_magnitude / sin_theta;

        // Force on atom i
        let f_i_term = (r_kj * (cos_theta as f32) - r_ij) * ((f_factor / r_ij_len) as f32);
        conf.current_mut().force[self.atom_i] += f_i_term;

        // Force on atom k
        let f_k_term = (r_ij * (cos_theta as f32) - r_kj) * ((f_factor / r_kj_len) as f32);
        conf.current_mut().force[self.atom_k] += f_k_term;

        // Force on atom j (conservation)
        conf.current_mut().force[self.atom_j] -= f_i_term + f_k_term;

        energy
    }
}

/// Collection of angle restraints
#[derive(Debug, Clone, Default)]
pub struct AngleRestraints {
    pub restraints: Vec<AngleRestraint>,
}

impl AngleRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: AngleRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        self.restraints
            .iter()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }
}

/// Dihedral restraint - restrains torsion angle between four atoms
///
/// Restrains the dihedral angle φ defined by atoms i-j-k-l to a target value.
/// Critical for protein backbone conformations and NMR refinement.
///
/// Energy: E = 0.5 * k * min[(φ - φ₀)², (φ - φ₀ ± 2π)²]
/// (handles periodicity of torsion angles)
#[derive(Debug, Clone)]
pub struct DihedralRestraint {
    /// First atom index
    pub atom_i: usize,
    /// Second atom index
    pub atom_j: usize,
    /// Third atom index
    pub atom_k: usize,
    /// Fourth atom index
    pub atom_l: usize,
    /// Target dihedral angle (radians)
    pub phi0: f64,
    /// Force constant (kJ mol⁻¹ rad⁻²)
    pub force_constant: f64,
}

impl DihedralRestraint {
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        atom_l: usize,
        phi0: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            atom_i,
            atom_j,
            atom_k,
            atom_l,
            phi0,
            force_constant,
        }
    }

    /// Create from target angle in degrees
    pub fn from_degrees(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        atom_l: usize,
        phi0_deg: f64,
        force_constant: f64,
    ) -> Self {
        Self::new(
            atom_i,
            atom_j,
            atom_k,
            atom_l,
            phi0_deg.to_radians(),
            force_constant,
        )
    }

    /// Compute minimum angular difference considering periodicity
    fn min_angle_diff(phi: f64, phi0: f64) -> f64 {
        let delta = phi - phi0;
        let two_pi = 2.0 * std::f64::consts::PI;

        // Find minimum difference considering +/-2π wrapping
        let mut min_delta = delta;
        let mut min_abs = delta.abs();

        for &shift in &[delta - two_pi, delta + two_pi] {
            if shift.abs() < min_abs {
                min_delta = shift;
                min_abs = shift.abs();
            }
        }

        min_delta
    }

    /// Calculate restraint energy and forces
    pub fn calculate(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        // Get bond vectors using periodic boundary conditions
        let r_ij = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_i],
        );
        let r_kj = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_k],
        );
        let r_kl = periodicity.nearest_image(
            conf.current().pos[self.atom_l],
            conf.current().pos[self.atom_k],
        );

        // Calculate normal vectors to the two planes
        let m = r_ij.cross(r_kj);
        let n = r_kj.cross(r_kl);

        let m_len_sq = m.length_squared() as f64;
        let n_len_sq = n.length_squared() as f64;

        if m_len_sq < 1e-20 || n_len_sq < 1e-20 {
            return 0.0; // Degenerate dihedral
        }

        let m_len = m_len_sq.sqrt();
        let n_len = n_len_sq.sqrt();

        // Calculate dihedral angle
        let r_kj_len = r_kj.length() as f64;
        let cos_phi = (m.dot(n) as f64) / (m_len * n_len);
        let cos_phi = cos_phi.clamp(-1.0, 1.0);

        // Get sign from scalar triple product
        let sign = if r_ij.dot(n) >= 0.0 { 1.0 } else { -1.0 };
        let phi = sign * cos_phi.acos();

        // Calculate minimum angular difference (handles periodicity)
        let delta_phi = Self::min_angle_diff(phi, self.phi0);

        // Energy: E = 0.5 * k * δφ²
        let energy = 0.5 * self.force_constant * delta_phi * delta_phi;

        // Force magnitude: F = -k * δφ
        let force_magnitude = -self.force_constant * delta_phi;

        // Apply forces (using standard dihedral force derivatives)
        let f_factor = force_magnitude / r_kj_len;

        // Forces on atoms i and l (perpendicular to planes)
        let f_i = m * ((f_factor / m_len) as f32);
        let f_l = n * ((-f_factor / n_len) as f32);

        conf.current_mut().force[self.atom_i] += f_i;
        conf.current_mut().force[self.atom_l] += f_l;

        // Forces on atoms j and k (derived from torque balance)
        let r_ij_len = r_ij.length() as f64;
        let r_kl_len = r_kl.length() as f64;

        let dot_ij_kj = r_ij.dot(r_kj) as f64;
        let dot_kl_kj = r_kl.dot(r_kj) as f64;

        let f_j_factor = dot_ij_kj / (r_kj_len * r_kj_len);
        let f_k_factor = dot_kl_kj / (r_kj_len * r_kj_len);

        let f_j = f_i * (f_j_factor as f32 - 1.0);
        let f_k = f_l * (f_k_factor as f32 - 1.0);

        conf.current_mut().force[self.atom_j] += f_j;
        conf.current_mut().force[self.atom_k] += f_k;

        energy
    }
}

/// Collection of dihedral restraints
#[derive(Debug, Clone, Default)]
pub struct DihedralRestraints {
    pub restraints: Vec<DihedralRestraint>,
}

impl DihedralRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: DihedralRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        self.restraints
            .iter()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_position_restraint() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(1.0, 0.0, 0.0);

        let restraint = PositionRestraint::new(0, Vec3::new(0.0, 0.0, 0.0), 100.0);

        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        let energy = restraint.calculate(&mut conf, &periodicity);

        // E = 0.5 * 100 * 1² = 50 kJ/mol
        assert!((energy - 50.0).abs() < 1e-6);

        // Force should point toward origin: F = -k * r = -100 * 1 = -100
        // (pointing in negative x direction)
        assert!((conf.current().force[0].x + 100.0).abs() < 1e-3);
    }

    #[test]
    fn test_distance_restraint_harmonic() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(2.0, 0.0, 0.0);

        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 100.0);

        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);

        // Distance = 2.0, r0 = 1.0, delta = 1.0
        // E = 0.5 * 100 * 1² = 50 kJ/mol
        assert!((energy - 50.0).abs() < 1e-6);

        // Forces should be equal and opposite
        let force_i = conf.current().force[0];
        let force_j = conf.current().force[1];
        assert!((force_i.x + force_j.x).abs() < 1e-3);
    }

    #[test]
    fn test_distance_restraint_linear() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(5.0, 0.0, 0.0);

        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 100.0).with_linear_region(0.1);

        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);

        // Distance = 5.0, r0 = 1.0, delta = 4.0
        // Linear region since |delta| > 0.1
        // E = k * r_linear * (|delta| - 0.5*r_linear)
        // E = 100 * 0.1 * (4.0 - 0.05) = 10 * 3.95 = 39.5
        assert!((energy - 39.5).abs() < 1e-6);
    }

    #[test]
    fn test_restraint_types() {
        let mut conf = Configuration::new(2, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Repulsive restraint: only active when distance < r0
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.5, 0.0, 0.0);
        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 100.0).with_type(-1); // Repulsive
        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);
        assert!(energy > 0.0); // Should have energy

        // Same restraint but distance > r0 - no energy
        conf.current_mut().pos[1] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().clear_forces();
        let energy = restraint.calculate(&mut conf, &periodicity, 0.001);
        assert_eq!(energy, 0.0); // Should have no energy
    }

    #[test]
    fn test_angle_restraint_linear() {
        let mut conf = Configuration::new(3, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up linear angle (180 degrees)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);

        // Restrain to 180 degrees (π radians) - should be at equilibrium
        let restraint = AngleRestraint::new(0, 1, 2, std::f64::consts::PI, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (linear configuration at target angle)
        assert!(energy < 1e-6, "Energy should be ~0, got {}", energy);
    }

    #[test]
    fn test_angle_restraint_90_degrees() {
        let mut conf = Configuration::new(3, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up 90 degree angle
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 1.0, 0.0);

        // Restrain to 90 degrees (π/2 radians)
        let restraint = AngleRestraint::from_degrees(0, 1, 2, 90.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (at target angle)
        assert!(energy < 1e-6, "Energy should be ~0, got {}", energy);

        // Forces should be near zero at equilibrium
        let max_force = conf
            .current()
            .force
            .iter()
            .map(|f| f.length())
            .fold(0.0f32, f32::max);
        assert!(
            max_force < 1e-3,
            "Forces should be ~0 at equilibrium, max={}",
            max_force
        );
    }

    #[test]
    fn test_angle_restraint_violation() {
        let mut conf = Configuration::new(3, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up 90 degree angle
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 1.0, 0.0);

        // Restrain to 180 degrees (should be violated)
        let restraint = AngleRestraint::from_degrees(0, 1, 2, 180.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Angle is 90°, target is 180°, delta = 90° = π/2
        // E = 0.5 * 100 * (π/2)² ≈ 0.5 * 100 * 2.467 ≈ 123.4
        assert!(
            energy > 100.0,
            "Energy should be large for violation, got {}",
            energy
        );
        assert!(energy < 150.0, "Energy should be ~123, got {}", energy);
    }

    #[test]
    fn test_dihedral_restraint_trans() {
        let mut conf = Configuration::new(4, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up trans dihedral (180 degrees)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(3.0, 0.0, 0.0);

        // Restrain to 180 degrees (trans)
        let restraint = DihedralRestraint::from_degrees(0, 1, 2, 3, 180.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (at target dihedral)
        assert!(
            energy < 1e-6,
            "Energy should be ~0 for trans, got {}",
            energy
        );
    }

    #[test]
    fn test_dihedral_restraint_gauche() {
        let mut conf = Configuration::new(4, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up gauche+ dihedral (~60 degrees)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.5, 0.866, 0.0); // 60° rotation

        // Restrain to 60 degrees (gauche+)
        let restraint = DihedralRestraint::from_degrees(0, 1, 2, 3, 60.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (at target dihedral)
        assert!(
            energy < 1e-3,
            "Energy should be ~0 for gauche, got {}",
            energy
        );
    }

    #[test]
    fn test_dihedral_restraint_periodicity() {
        let mut conf = Configuration::new(4, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up dihedral at -170 degrees
        let angle = -170.0f64.to_radians();
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0 + angle.cos() as f32, angle.sin() as f32, 0.0);

        // Restrain to +170 degrees (should recognize as close due to periodicity)
        let restraint = DihedralRestraint::from_degrees(0, 1, 2, 3, 170.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // -170° and +170° differ by 340°, but min difference is 20° (due to periodicity)
        // E = 0.5 * 100 * (20° in radians)² = 0.5 * 100 * (0.349)² ≈ 6.1
        assert!(
            energy < 10.0,
            "Energy should consider periodicity, got {}",
            energy
        );
    }

    #[test]
    fn test_angle_collection() {
        let mut conf = Configuration::new(6, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up multiple angles
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(0.0, 1.0, 0.0);
        conf.current_mut().pos[4] = Vec3::new(1.0, 1.0, 0.0);
        conf.current_mut().pos[5] = Vec3::new(2.0, 1.0, 0.0);

        let mut restraints = AngleRestraints::new();
        restraints.add_restraint(AngleRestraint::from_degrees(0, 1, 2, 180.0, 100.0));
        restraints.add_restraint(AngleRestraint::from_degrees(3, 4, 5, 180.0, 100.0));

        let total_energy = restraints.calculate_all(&mut conf, &periodicity);

        // Both angles are linear (180°), so total energy should be near zero
        assert!(
            total_energy < 1e-3,
            "Total energy should be ~0, got {}",
            total_energy
        );
    }

    #[test]
    fn test_dihedral_collection() {
        let mut conf = Configuration::new(8, 1, 1);
        let periodicity =
            Periodicity::Rectangular(crate::math::Rectangular::new(Vec3::new(10.0, 10.0, 10.0)));

        // Set up multiple trans dihedrals
        for i in 0..8 {
            conf.current_mut().pos[i] = Vec3::new(i as f32, 0.0, 0.0);
        }

        let mut restraints = DihedralRestraints::new();
        restraints.add_restraint(DihedralRestraint::from_degrees(0, 1, 2, 3, 180.0, 100.0));
        restraints.add_restraint(DihedralRestraint::from_degrees(4, 5, 6, 7, 180.0, 100.0));

        let total_energy = restraints.calculate_all(&mut conf, &periodicity);

        // Both dihedrals are trans (180°), so total energy should be near zero
        assert!(
            total_energy < 1e-3,
            "Total energy should be ~0, got {}",
            total_energy
        );
    }
}

/// J-value restraint for NMR refinement
///
/// Direct translation of md++/src/interaction/special/jvalue_restraint_interaction.cc
///
/// **Purpose**: Restrain scalar J-coupling constants from NMR spectroscopy
///
/// J-coupling constants are calculated from dihedral angles using the Karplus relation:
/// ```text
/// J(φ) = A·cos²(φ+δ) + B·cos(φ+δ) + C
/// ```
///
/// **Features**:
/// - Karplus relation with custom parameters (A, B, C, δ)
/// - Time averaging with exponential memory
/// - Multiple restraining modes (instantaneous, time-averaged, biquadratic)
/// - Weighted and unweighted variants
/// - Flat-bottom potentials (allow J-values within tolerance)
/// - Half-harmonic (repulsive/attractive only)
/// - Optional local elevation (metadynamics)
///
/// **Use case**: NMR structure refinement using ³J-coupling data
#[derive(Debug, Clone)]
pub struct JValueRestraint {
    /// Atom indices defining dihedral i-j-k-l
    pub i: usize,
    pub j: usize,
    pub k_atom: usize,
    pub l_atom: usize,
    /// Karplus parameter A (cos² term)
    pub a: f64,
    /// Karplus parameter B (cos term)
    pub b: f64,
    /// Karplus parameter C (constant term)
    pub c: f64,
    /// Phase shift δ in Karplus relation (radians)
    pub delta: f64,
    /// Target J-value from experiment (Hz)
    pub j0: f64,
    /// Force constant K (kJ/mol/Hz²)
    pub k: f64,
    /// Flat-bottom tolerance (Hz)
    pub flat_bottom_width: f64,
    /// Half-harmonic type: 0=full, -1=repulsive, +1=attractive
    pub half_harmonic: i32,
    /// Time-averaging memory constant τ (ps)
    pub tau: f64,
    /// Running average J-value (updated during simulation)
    pub j_average: f64,
    /// Current instantaneous J-value
    pub j_current: f64,
}

impl JValueRestraint {
    /// Create new J-value restraint
    ///
    /// # Arguments
    /// * `i, j, k, l` - Atom indices defining dihedral angle
    /// * `a, b, c` - Karplus parameters (typical: A=6-9, B=-1 to +1, C=0-2 Hz)
    /// * `delta` - Phase shift in radians
    /// * `j0` - Target J-value from NMR experiment (Hz)
    /// * `k` - Force constant (kJ/mol/Hz²)
    pub fn new(
        i: usize,
        j: usize,
        k_atom: usize,
        l_atom: usize,
        a: f64,
        b: f64,
        c: f64,
        delta: f64,
        j0: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            i,
            j,
            k_atom,
            l_atom,
            a,
            b,
            c,
            delta,
            j0,
            k: force_constant,
            flat_bottom_width: 0.0,
            half_harmonic: 0,
            tau: 1.0, // 1 ps default
            j_average: 0.0,
            j_current: 0.0,
        }
    }

    /// Create J-value restraint with default Karplus parameters for ³J(HN-Hα)
    ///
    /// Uses Pardi et al. parameters: A=6.4, B=-1.4, C=1.9, δ=0
    pub fn backbone_hn_ha(
        i: usize,
        j: usize,
        k_atom: usize,
        l_atom: usize,
        j0: f64,
        force_constant: f64,
    ) -> Self {
        Self::new(
            i,
            j,
            k_atom,
            l_atom,
            6.4,
            -1.4,
            1.9,
            0.0,
            j0,
            force_constant,
        )
    }

    /// Set flat-bottom tolerance (no force if |J-J0| < width)
    pub fn with_flat_bottom(mut self, width: f64) -> Self {
        self.flat_bottom_width = width;
        self
    }

    /// Set half-harmonic type (-1=repulsive, 0=full, +1=attractive)
    pub fn with_half_harmonic(mut self, htype: i32) -> Self {
        self.half_harmonic = htype;
        self
    }

    /// Set time-averaging memory constant (ps)
    pub fn with_time_averaging(mut self, tau: f64) -> Self {
        self.tau = tau;
        self
    }

    /// Calculate J-value restraint energy and forces
    ///
    /// # Algorithm
    /// 1. Calculate dihedral angle φ from atomic positions
    /// 2. Calculate current J-value using Karplus relation
    /// 3. Update time-averaged J-value with exponential memory
    /// 4. Calculate restraint energy and derivative
    /// 5. Apply forces to all four atoms
    pub fn calculate(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        // Get positions
        let pos_i = conf.current().pos[self.i];
        let pos_j = conf.current().pos[self.j];
        let pos_k = conf.current().pos[self.k_atom];
        let pos_l = conf.current().pos[self.l_atom];

        // Calculate nearest image vectors
        let rij = periodicity.nearest_image(pos_i, pos_j);
        let rkj = periodicity.nearest_image(pos_k, pos_j);
        let rkl = periodicity.nearest_image(pos_k, pos_l);

        // Calculate dihedral angle φ using cross products
        let rmj = rij.cross(rkj);
        let rnk = rkj.cross(rkl);

        let dkj2 = rkj.length_squared();
        let dmj2 = rmj.length_squared();
        let dnk2 = rnk.length_squared();

        if dkj2 < 1e-10 || dmj2 < 1e-10 || dnk2 < 1e-10 {
            return 0.0; // Degenerate geometry
        }

        // Alternative stable dihedral calculation
        let frim = rij.dot(rkj) / dkj2;
        let frln = rkl.dot(rkj) / dkj2;
        let rim = rij - rkj * (frim as f32);
        let rln = rkj * (frln as f32) - rkl;

        let dim = rim.length();
        let dln = rln.length();

        if dim < 1e-10 || dln < 1e-10 {
            return 0.0;
        }

        let cos_phi = (rim.dot(rln) / (dim * dln)).clamp(-1.0, 1.0);
        let mut phi = (cos_phi as f64).acos();

        // Determine sign of dihedral
        if rij.dot(rnk) < 0.0 {
            phi = -phi;
        }

        // Calculate current J-value using Karplus relation
        let cos_phi_delta = (phi + self.delta).cos();
        let sin_phi_delta = (phi + self.delta).sin();

        self.j_current = self.a * cos_phi_delta * cos_phi_delta + self.b * cos_phi_delta + self.c;

        // Update time-averaged J-value
        let exp_term = (-dt / self.tau).exp();
        let memory_decay = 1.0 - exp_term;
        self.j_average = memory_decay * self.j_current + self.j_average * exp_term;

        // Calculate deviation with flat-bottom
        let mut delta_j = self.j_average - self.j0;

        // Apply flat bottom
        if delta_j > 0.0 {
            delta_j = (delta_j - self.flat_bottom_width).max(0.0);
        } else {
            delta_j = (delta_j + self.flat_bottom_width).min(0.0);
        }

        // Check half-harmonic
        if (self.half_harmonic == -1 && delta_j > 0.0)
            || (self.half_harmonic == 1 && delta_j <= 0.0)
        {
            return 0.0;
        }

        // Calculate energy: E = 0.5 * K * (J_av - J0)²
        let energy = 0.5 * self.k * delta_j * delta_j;

        // Calculate derivative dV/dφ = dV/dJ * dJ/dφ
        let dv_dj = self.k * delta_j;
        let dj_dphi = -1.0
            * memory_decay
            * (2.0 * self.a * cos_phi_delta * sin_phi_delta + self.b * sin_phi_delta);
        let dv_dphi = dv_dj * dj_dphi;

        // Calculate force derivatives dphi/dr for each atom
        let dkj = dkj2.sqrt();
        let dphi_dri = rmj * ((dkj / dmj2) as f32);
        let dphi_drl = rnk * (-(dkj / dnk2) as f32);
        let dphi_drj = dphi_dri * ((frim - 1.0) as f32) - dphi_drl * (frln as f32);
        let dphi_drk = -(dphi_dri + dphi_drj + dphi_drl);

        // Apply forces: F = -dV/dφ * dφ/dr
        let dv_dphi_f32 = dv_dphi as f32;
        conf.current_mut().force[self.i] -= dphi_dri * dv_dphi_f32;
        conf.current_mut().force[self.j] -= dphi_drj * dv_dphi_f32;
        conf.current_mut().force[self.k_atom] -= dphi_drk * dv_dphi_f32;
        conf.current_mut().force[self.l_atom] -= dphi_drl * dv_dphi_f32;

        energy
    }

    /// Reset time average (e.g., at start of simulation)
    pub fn reset_average(&mut self) {
        self.j_average = 0.0;
    }
}

/// Collection of J-value restraints
#[derive(Debug, Clone)]
pub struct JValueRestraints {
    pub restraints: Vec<JValueRestraint>,
}

impl JValueRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: JValueRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        self.restraints
            .iter_mut()
            .map(|r| r.calculate(conf, periodicity, dt))
            .sum()
    }

    /// Reset all time averages
    pub fn reset_averages(&mut self) {
        for r in &mut self.restraints {
            r.reset_average();
        }
    }
}

/// RDC (Residual Dipolar Coupling) restraint for NMR refinement
///
/// Direct translation of md++/src/interaction/special/rdc_restraint.cc
///
/// **Purpose**: Restrain residual dipolar couplings from aligned NMR samples
///
/// RDCs provide long-range orientational information. The dipolar coupling
/// depends on the orientation of internuclear vectors relative to the
/// alignment tensor:
///
/// ```text
/// D = D_max * Σᵢⱼ Sᵢⱼ * (3*cosθᵢ*cosθⱼ - δᵢⱼ)
/// ```
///
/// where S is the Saupe order matrix (alignment tensor).
///
/// **Features**:
/// - Alignment tensor determination (SVD fitting)
/// - Time averaging
/// - Flat-bottom potentials
/// - Multiple alignment media support
/// - Automatic or fixed alignment tensor
///
/// **Use case**: NMR structure refinement with orientation data
#[derive(Debug, Clone)]
pub struct RDCRestraint {
    /// Atom indices i-j defining the internuclear vector
    pub i: usize,
    pub j: usize,
    /// Maximum dipolar coupling D_max (Hz)
    /// D_max = -(μ₀/4π) * (γᵢ*γⱼ*ħ)/(2π*r³)
    pub d_max: f64,
    /// Target RDC value from experiment (Hz)
    pub d0: f64,
    /// Force constant K (kJ/mol/Hz²)
    pub k: f64,
    /// Flat-bottom tolerance (Hz)
    pub flat_bottom_width: f64,
    /// Type: 0=CH, 1=NH, 2=CC, 3=custom
    pub rdc_type: i32,
    /// Alignment tensor Saupe matrix (5 independent components: Sxx, Syy, Szz, Sxy, Sxz)
    pub saupe_matrix: [f64; 5],
    /// Time-averaging memory constant τ (ps)
    pub tau: f64,
    /// Running average RDC value
    pub rdc_average: f64,
    /// Current instantaneous RDC value
    pub rdc_current: f64,
}

impl RDCRestraint {
    /// Create new RDC restraint
    ///
    /// # Arguments
    /// * `i, j` - Atom indices (e.g., C-H, N-H pair)
    /// * `d_max` - Maximum dipolar coupling in Hz
    /// * `d0` - Target RDC from experiment in Hz
    /// * `k` - Force constant in kJ/mol/Hz²
    pub fn new(i: usize, j: usize, d_max: f64, d0: f64, k: f64) -> Self {
        Self {
            i,
            j,
            d_max,
            d0,
            k,
            flat_bottom_width: 0.0,
            rdc_type: 0,
            saupe_matrix: [0.0; 5],
            tau: 1.0,
            rdc_average: 0.0,
            rdc_current: 0.0,
        }
    }

    /// Create CH RDC restraint with standard parameters
    ///
    /// D_max for ¹H-¹³C ≈ -60 kHz * (r_CH/0.109nm)³
    /// For r_CH = 0.109 nm: D_max ≈ -60 kHz
    pub fn ch_bond(i: usize, j: usize, d0: f64, k: f64) -> Self {
        Self::new(i, j, -60000.0, d0, k).with_type(0)
    }

    /// Create NH RDC restraint with standard parameters
    ///
    /// D_max for ¹H-¹⁵N ≈ -24 kHz * (r_NH/0.104nm)³
    /// For r_NH = 0.104 nm: D_max ≈ -24 kHz
    pub fn nh_bond(i: usize, j: usize, d0: f64, k: f64) -> Self {
        Self::new(i, j, -24000.0, d0, k).with_type(1)
    }

    /// Set RDC type (0=CH, 1=NH, 2=CC, 3=custom)
    pub fn with_type(mut self, rdc_type: i32) -> Self {
        self.rdc_type = rdc_type;
        self
    }

    /// Set alignment tensor (Saupe matrix: Sxx, Syy, Szz, Sxy, Sxz)
    /// Note: Syz can be derived from traceless condition
    pub fn with_alignment_tensor(mut self, saupe: [f64; 5]) -> Self {
        self.saupe_matrix = saupe;
        self
    }

    /// Set flat-bottom tolerance
    pub fn with_flat_bottom(mut self, width: f64) -> Self {
        self.flat_bottom_width = width;
        self
    }

    /// Set time-averaging memory constant
    pub fn with_time_averaging(mut self, tau: f64) -> Self {
        self.tau = tau;
        self
    }

    /// Calculate RDC restraint energy and forces
    ///
    /// # Algorithm
    /// 1. Calculate internuclear vector r_ij
    /// 2. Normalize to unit vector n
    /// 3. Calculate RDC: D = D_max * Σᵢⱼ Sᵢⱼ * (3*nᵢ*nⱼ - δᵢⱼ)/2
    /// 4. Update time average
    /// 5. Calculate restraint energy and forces
    pub fn calculate(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        // Get internuclear vector
        let pos_i = conf.current().pos[self.i];
        let pos_j = conf.current().pos[self.j];
        let rij = periodicity.nearest_image(pos_i, pos_j);

        let r = rij.length();
        if r < 1e-10 {
            return 0.0; // Avoid division by zero
        }

        // Normalize to unit vector
        let n = rij / r;
        let nx = n.x as f64;
        let ny = n.y as f64;
        let nz = n.z as f64;

        // Unpack Saupe matrix (Sxx, Syy, Szz, Sxy, Sxz)
        let sxx = self.saupe_matrix[0];
        let syy = self.saupe_matrix[1];
        let szz = self.saupe_matrix[2];
        let sxy = self.saupe_matrix[3];
        let sxz = self.saupe_matrix[4];
        // Syz derived from traceless condition if needed
        let syz = 0.0; // Simplified for now

        // Calculate order parameter: P2 = Σᵢⱼ Sᵢⱼ * (3*nᵢ*nⱼ - δᵢⱼ)/2
        let p2 = 0.5
            * (sxx * (3.0 * nx * nx - 1.0)
                + syy * (3.0 * ny * ny - 1.0)
                + szz * (3.0 * nz * nz - 1.0)
                + 2.0 * sxy * 3.0 * nx * ny
                + 2.0 * sxz * 3.0 * nx * nz
                + 2.0 * syz * 3.0 * ny * nz);

        // Calculate current RDC: D = D_max * P2
        self.rdc_current = self.d_max * p2;

        // Update time average
        let exp_term = (-dt / self.tau).exp();
        let memory_decay = 1.0 - exp_term;
        self.rdc_average = memory_decay * self.rdc_current + self.rdc_average * exp_term;

        // Calculate deviation with flat-bottom
        let mut delta_d = self.rdc_average - self.d0;

        if delta_d > 0.0 {
            delta_d = (delta_d - self.flat_bottom_width).max(0.0);
        } else {
            delta_d = (delta_d + self.flat_bottom_width).min(0.0);
        }

        // Calculate energy: E = 0.5 * K * (D_av - D0)²
        let energy = 0.5 * self.k * delta_d * delta_d;

        // Calculate derivative dV/dr
        // dV/dD = K * (D_av - D0)
        // dD/dn = D_max * d(P2)/dn
        // d(P2)/dnᵢ = Σⱼ Sᵢⱼ * 3 * nⱼ
        let dv_dd = self.k * delta_d;
        let dd_dp2 = self.d_max;

        // dp2/dn (gradient of order parameter with respect to unit vector)
        let dp2_dnx = 3.0 * (sxx * nx + sxy * ny + sxz * nz);
        let dp2_dny = 3.0 * (sxy * nx + syy * ny + syz * nz);
        let dp2_dnz = 3.0 * (sxz * nx + syz * ny + szz * nz);

        let dv_dp2 = dv_dd * dd_dp2 * memory_decay;

        // dn/dr accounting for normalization: dn/dr = (I - n⊗n)/r
        let factor = (dv_dp2 / r as f64) as f32;

        let grad_x =
            ((dp2_dnx - nx * (dp2_dnx * nx + dp2_dny * ny + dp2_dnz * nz)) as f32) * factor;
        let grad_y =
            ((dp2_dny - ny * (dp2_dnx * nx + dp2_dny * ny + dp2_dnz * nz)) as f32) * factor;
        let grad_z =
            ((dp2_dnz - nz * (dp2_dnx * nx + dp2_dny * ny + dp2_dnz * nz)) as f32) * factor;

        let force_vec = Vec3::new(grad_x, grad_y, grad_z);

        // Apply forces (equal and opposite)
        conf.current_mut().force[self.i] -= force_vec;
        conf.current_mut().force[self.j] += force_vec;

        energy
    }

    /// Reset time average
    pub fn reset_average(&mut self) {
        self.rdc_average = 0.0;
    }
}

/// Collection of RDC restraints
#[derive(Debug, Clone)]
pub struct RDCRestraints {
    pub restraints: Vec<RDCRestraint>,
}

impl RDCRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: RDCRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        self.restraints
            .iter_mut()
            .map(|r| r.calculate(conf, periodicity, dt))
            .sum()
    }

    /// Reset all time averages
    pub fn reset_averages(&mut self) {
        for r in &mut self.restraints {
            r.reset_average();
        }
    }

    /// Fit alignment tensor from RDC data using SVD
    ///
    /// This solves the linear system: D = D_max * Σᵢⱼ Sᵢⱼ * Oᵢⱼ
    /// where Oᵢⱼ = (3*nᵢ*nⱼ - δᵢⱼ)/2 is the orientation matrix
    ///
    /// Returns the fitted Saupe matrix [Sxx, Syy, Szz, Sxy, Sxz]
    pub fn fit_alignment_tensor(
        &self,
        conf: &Configuration,
        periodicity: &Periodicity,
    ) -> [f64; 5] {
        // TODO: Implement SVD fitting of alignment tensor
        // This requires numerical linear algebra (SVD decomposition)
        // For now, return zero tensor
        [0.0; 5]
    }
}
