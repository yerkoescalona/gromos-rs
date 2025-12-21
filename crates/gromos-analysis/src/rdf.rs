//! Radial Distribution Function (RDF) calculation

use gromos_core::Vec3;

/// Calculate RDF between two groups of atoms
pub fn calculate_rdf(
    _positions1: &[Vec3],
    _positions2: &[Vec3],
    _box_size: &[f64; 3],
    _bins: usize,
    _rmax: f64,
) -> Vec<(f64, f64)> {
    // TODO: Implement RDF calculation
    Vec::new()
}
