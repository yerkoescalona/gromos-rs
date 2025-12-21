//! Radius of gyration calculation

use gromos_core::Vec3;

/// Calculate radius of gyration
pub fn calculate_gyration_radius(positions: &[Vec3], masses: &[f64]) -> f64 {
    if positions.is_empty() || positions.len() != masses.len() {
        return 0.0;
    }

    let total_mass: f64 = masses.iter().sum();
    
    // Calculate center of mass
    let com = positions
        .iter()
        .zip(masses.iter())
        .fold(Vec3::ZERO, |acc, (p, m)| {
            Vec3::new(
                acc.x + (p.x as f64 * m) as f32,
                acc.y + (p.y as f64 * m) as f32,
                acc.z + (p.z as f64 * m) as f32,
            )
        });
    let com = Vec3::new(
        (com.x as f64 / total_mass) as f32,
        (com.y as f64 / total_mass) as f32,
        (com.z as f64 / total_mass) as f32,
    );

    // Calculate Rg^2
    let rg2: f64 = positions
        .iter()
        .zip(masses.iter())
        .map(|(p, m)| {
            let dx = p.x - com.x;
            let dy = p.y - com.y;
            let dz = p.z - com.z;
            m * (dx * dx + dy * dy + dz * dz) as f64
        })
        .sum::<f64>() / total_mass;

    rg2.sqrt()
}
