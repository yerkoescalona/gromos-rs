//! LJ/CRF parameter types and force storage

use gromos_core::math::Vec3;

/// Lennard-Jones parameters for atom pair
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct LJParameters {
    pub c6: f64,  // LJ C6 coefficient (attraction)
    pub c12: f64, // LJ C12 coefficient (repulsion)
    pub cs6: f64,  // 1-4 C6 coefficient
    pub cs12: f64, // 1-4 C12 coefficient
}

impl From<&gromos_core::topology::LJParameters> for LJParameters {
    fn from(p: &gromos_core::topology::LJParameters) -> Self {
        Self { c6: p.c6, c12: p.c12, cs6: p.cs6, cs12: p.cs12 }
    }
}

/// Flat LJ parameter matrix for cache-friendly access.
///
/// Stores parameters in a contiguous array indexed as `[type_i * n_types + type_j]`.
/// Eliminates the double indirection of `Vec<Vec<LJParameters>>`.
#[derive(Debug, Clone)]
pub struct LJParamMatrix {
    data: Vec<LJParameters>,
    n_types: usize,
}

impl LJParamMatrix {
    /// Create from a nested Vec (converts from legacy format).
    pub fn from_nested(nested: &[Vec<LJParameters>]) -> Self {
        let n_types = nested.len();
        let mut data = Vec::with_capacity(n_types * n_types);
        for row in nested {
            assert_eq!(row.len(), n_types, "LJ parameter matrix must be square");
            data.extend_from_slice(row);
        }
        Self { data, n_types }
    }

    /// Get parameters for atom type pair (i, j).
    #[inline(always)]
    pub fn get(&self, type_i: usize, type_j: usize) -> LJParameters {
        unsafe { *self.data.get_unchecked(type_i * self.n_types + type_j) }
    }

    /// Number of atom types.
    pub fn n_types(&self) -> usize {
        self.n_types
    }
}

/// Coulomb constant: 1/(4*pi*eps0) in GROMOS units [kJ*nm/(mol*e^2)]
pub const FOUR_PI_EPS_I: f64 = 138.9354;

/// Coulomb Reaction Field parameters
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CRFParameters {
    pub crf_cut: f64,    // (1 - crf/2) / cutoff  (energy shift constant)
    pub crf_2cut3i: f64, // crf / (2 * cutoff^3)  (energy quadratic term)
    pub crf_cut3i: f64,  // crf / cutoff^3         (force correction term)
    pub cutoff_sq: f64,  // cutoff^2 for Heaviside truncation (GROMOS XXHEAVISIDE)
}

impl CRFParameters {
    /// Calculate CRF parameters from physical constants.
    ///
    /// # Arguments
    /// * `cutoff` - Nonbonded cutoff distance in nm
    /// * `epsilon` - Dielectric permittivity inside the cutoff
    /// * `rf_epsilon` - Dielectric permittivity outside the cutoff (reaction field)
    /// * `rf_kappa` - Debye-Hückel screening parameter in nm⁻¹
    pub fn new(cutoff: f64, epsilon: f64, rf_epsilon: f64, rf_kappa: f64) -> Self {
        let kappa_cut = rf_kappa * cutoff;
        let kappa_cut2 = kappa_cut * kappa_cut;

        let crf = (2.0 * (epsilon - rf_epsilon) * (1.0 + kappa_cut)
            - rf_epsilon * kappa_cut2)
            / ((epsilon + 2.0 * rf_epsilon) * (1.0 + kappa_cut)
                + rf_epsilon * kappa_cut2);

        let cut3i = 1.0 / (cutoff * cutoff * cutoff);
        let crf_cut3i = crf * cut3i;
        let crf_2cut3i = crf_cut3i / 2.0;
        let crf_cut = (1.0 - crf / 2.0) / cutoff;

        Self {
            crf_cut,
            crf_2cut3i,
            crf_cut3i,
            cutoff_sq: cutoff * cutoff,
        }
    }
}

/// Storage for forces and energies
#[repr(C)]
pub struct ForceStorage {
    pub forces: Vec<Vec3>,
    pub e_lj: f64,
    pub e_crf: f64,
    pub virial: [[f64; 3]; 3],
}

impl ForceStorage {
    pub fn new(n_atoms: usize) -> Self {
        Self {
            forces: vec![Vec3::ZERO; n_atoms],
            e_lj: 0.0,
            e_crf: 0.0,
            virial: [[0.0; 3]; 3],
        }
    }

    pub fn clear(&mut self) {
        self.forces.fill(Vec3::ZERO);
        self.e_lj = 0.0;
        self.e_crf = 0.0;
        self.virial = [[0.0; 3]; 3];
    }

    /// Accumulate another ForceStorage into this one (thread-local reduction).
    pub fn merge(&mut self, other: &ForceStorage) {
        for (dst, src) in self.forces.iter_mut().zip(other.forces.iter()) {
            *dst += *src;
        }
        self.e_lj += other.e_lj;
        self.e_crf += other.e_crf;
        for a in 0..3 {
            for b in 0..3 {
                self.virial[a][b] += other.virial[a][b];
            }
        }
    }
}

/// A group of atom pairs sharing the same PBC shift (same charge-group pair).
///
/// The `ref_atom_i` and `ref_atom_j` are used to compute nearest_image once;
/// all atom pairs in the group then use the resulting shift offset.
#[derive(Debug, Clone, Copy)]
pub struct CGPairGroup {
    /// Reference atom from CG i (typically first atom in the CG)
    pub ref_atom_i: u32,
    /// Reference atom from CG j (typically first atom in the CG)
    pub ref_atom_j: u32,
    /// Start index into the flat pair array
    pub start: u32,
    /// End index (exclusive) into the flat pair array
    pub end: u32,
}
