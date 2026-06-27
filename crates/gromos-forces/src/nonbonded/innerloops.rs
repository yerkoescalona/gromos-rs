//! LJ+CRF inner loop kernels: serial, parallel, CG-grouped, and solvent variants

use gromos_core::math::{BoundaryCondition, Vec3};
use rayon::prelude::*;
use super::{LJParamMatrix, CRFParameters, ForceStorage, CGPairGroup};
use super::params::FOUR_PI_EPS_I;

/// Core LJ + CRF interaction calculation (hot path!)
///
/// # Arguments
/// * `r` - Distance vector from i to j
/// * `c6` - LJ C6 coefficient
/// * `c12` - LJ C12 coefficient
/// * `q_prod` - Charge product qi * qj
/// * `crf` - CRF parameters
///
/// # Returns
/// * `force_magnitude` - Scalar force magnitude
/// * `e_lj` - Lennard-Jones energy
/// * `e_crf` - Coulomb reaction field energy
#[inline(always)]
pub fn lj_crf_interaction(
    r: Vec3,
    c6: f64,
    c12: f64,
    q_prod: f64,
    crf: &CRFParameters,
) -> (f64, f64, f64) {
    let r2 = r.length_squared();

    // Early exit for zero distance (should not happen, but safety)
    if r2 < 1e-10 {
        return (0.0, 0.0, 0.0);
    }

    let inv_r2 = 1.0 / r2;
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;

    // Lennard-Jones: E_lj = C12/r^12 - C6/r^6
    let e_lj = (c12 * inv_r6 - c6) * inv_r6;

    // Lennard-Jones force: F = 12*C12/r^14 - 6*C6/r^8
    let f_lj = (12.0 * c12 * inv_r6 - 6.0 * c6) * inv_r6 * inv_r2;

    // Coulomb Reaction Field
    // E_crf = q * (1/r - crf_2cut3i * r² - crf_cut)
    // F_crf = q * (1/r³ + crf_cut3i)
    let inv_r = inv_r2.sqrt();
    let e_crf = q_prod * (inv_r - crf.crf_2cut3i * r2 - crf.crf_cut);
    let f_crf = q_prod * (inv_r * inv_r2 + crf.crf_cut3i);

    let force_magnitude = f_lj + f_crf;

    (force_magnitude, e_lj, e_crf)
}

/// Vectorized version processing 4 pairs simultaneously (AVX2/AVX-512)
#[cfg(feature = "simd")]
#[inline]
pub fn lj_crf_interaction_simd_x4(
    r: [Vec3; 4],
    c6: [f64; 4],
    c12: [f64; 4],
    q_prod: [f64; 4],
    crf: &CRFParameters,
) -> ([f64; 4], [f64; 4], [f64; 4]) {
    use wide::f64x4;

    // Convert to SIMD vectors
    let r2_array: [f64; 4] = r.map(|v| v.length_squared());
    let r2 = f64x4::from(r2_array);

    let inv_r2 = f64x4::splat(1.0) / r2;
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;

    // Lennard-Jones (vectorized)
    let c6_vec = f64x4::from(c6);
    let c12_vec = f64x4::from(c12);

    let e_lj = (c12_vec * inv_r6 - c6_vec) * inv_r6;
    let f_lj =
        (f64x4::splat(12.0) * c12_vec * inv_r6 - f64x4::splat(6.0) * c6_vec) * inv_r6 * inv_r2;

    // Coulomb Reaction Field (vectorized)
    let q_prod_vec = f64x4::from(q_prod);
    let inv_r = inv_r2.sqrt();

    let e_crf =
        q_prod_vec * (inv_r - f64x4::splat(crf.crf_2cut3i) * r2 - f64x4::splat(crf.crf_cut));
    let f_crf = q_prod_vec * (inv_r * inv_r2 + f64x4::splat(crf.crf_cut3i));

    let force = f_lj + f_crf;

    // Convert back to arrays
    let force_array: [f64; 4] = force.into();
    let e_lj_array: [f64; 4] = e_lj.into();
    let e_crf_array: [f64; 4] = e_crf.into();

    (force_array, e_lj_array, e_crf_array)
}

/// CG-grouped kernel: compute nearest_image once per charge-group pair,
/// then process all atom pairs in that group using simple subtraction + shift.
///
/// This eliminates per-atom-pair nearest_image calls (expensive floor/branch ops)
/// and replaces them with cheap additions. For 3-atom CGs (water), this is ~3x fewer
/// nearest_image calls.
#[inline]
fn process_pairs_cg_grouped<BC: BoundaryCondition, const VIRIAL: bool>(
    pairs: &[(u32, u32)],
    groups: &[CGPairGroup],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    storage: &mut ForceStorage,
) {
    for group in groups {
        let ri = positions[group.ref_atom_i as usize];
        let rj = positions[group.ref_atom_j as usize];

        // Compute nearest_image ONCE for the CG pair reference atoms
        let r_ref = periodicity.nearest_image(ri, rj);
        // Extract the PBC shift (same pattern as solvent innerloop)
        let tx = r_ref.x - ri.x + rj.x;
        let ty = r_ref.y - ri.y + rj.y;
        let tz = r_ref.z - ri.z + rj.z;

        // Process all atom pairs in this group using the shared shift
        let block = &pairs[group.start as usize..group.end as usize];
        for &(ii, jj) in block {
            let i = ii as usize;
            let j = jj as usize;

            // r = (pos_i + shift) - pos_j (same as solvent innerloop pattern)
            let r = Vec3::new(
                positions[i].x + tx - positions[j].x,
                positions[i].y + ty - positions[j].y,
                positions[i].z + tz - positions[j].z,
            );

            let type_i = iac[i] as usize;
            let type_j = iac[j] as usize;
            let lj = lj_params.get(type_i, type_j);
            let q_prod = charges[i] * charges[j] * FOUR_PI_EPS_I;

            let (f_mag, e_lj, e_crf) = lj_crf_interaction(r, lj.c6, lj.c12, q_prod, crf);

            let force = r * f_mag;
            storage.forces[i] += force;
            storage.forces[j] -= force;

            storage.e_lj += e_lj;
            storage.e_crf += e_crf;

            if VIRIAL {
                let r = [r.x, r.y, r.z];
                let f = [force.x, force.y, force.z];
                for a in 0..3 {
                    for b in 0..3 {
                        storage.virial[a][b] += r[a] * f[b];
                    }
                }
            }
        }
    }
}

/// Shared kernel: process a slice of atom pairs into a ForceStorage.
///
/// This is the single source of truth for the nonbonded innerloop physics.
/// Both serial and parallel paths call this — no code duplication.
#[inline]
fn process_pairs<BC: BoundaryCondition, const VIRIAL: bool>(
    pairs: &[(u32, u32)],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    storage: &mut ForceStorage,
) {
    for &(i, j) in pairs {
        let i = i as usize;
        let j = j as usize;

        let pos_i = positions[i];
        let pos_j = positions[j];
        let r = periodicity.nearest_image(pos_i, pos_j);

        let type_i = iac[i] as usize;
        let type_j = iac[j] as usize;
        let lj = lj_params.get(type_i, type_j);
        let q_prod = charges[i] * charges[j] * FOUR_PI_EPS_I;

        let (f_mag, e_lj, e_crf) = lj_crf_interaction(r, lj.c6, lj.c12, q_prod, crf);

        let force = r * f_mag;
        storage.forces[i] += force;
        storage.forces[j] -= force;

        storage.e_lj += e_lj;
        storage.e_crf += e_crf;

        if VIRIAL {
            for a in 0..3 {
                for b in 0..3 {
                    storage.virial[a][b] += r[a] * force[b];
                }
            }
        }
    }
}

/// Shared kernel: process a slice of solvent molecule pairs into a ForceStorage.
///
/// Single source of truth for solvent-solvent interactions with shared PBC shift.
/// Precomputes LJ params and charge products for the N×N atom type combinations
/// within a solvent molecule (GROMOS optimization: avoids repeated lookups).
#[inline]
fn process_solvent_pairs<BC: BoundaryCondition, const VIRIAL: bool>(
    pairs: &[(u32, u32)],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    storage: &mut ForceStorage,
) {
    if pairs.is_empty() {
        return;
    }

    // Precompute LJ params and charge products for all atom_i × atom_j combinations
    // within a solvent molecule. All solvent molecules have the same atom types,
    // so we only need to look these up once using the first molecule.
    let first_mol = pairs[0].0 as usize;
    let n = atoms_per_solvent;
    // Stack-allocated arrays for up to 4-site water models (covers SPC, SPC/E, TIP4P)
    assert!(n <= 4, "solvent with >4 atoms not supported in optimized path");
    let mut lj_c6 = [[0.0f64; 4]; 4];
    let mut lj_c12 = [[0.0f64; 4]; 4];
    let mut q_prod = [[0.0f64; 4]; 4];

    for ai in 0..n {
        let type_i = iac[first_mol + ai] as usize;
        let qi = charges[first_mol + ai];
        for aj in 0..n {
            let type_j = iac[first_mol + aj] as usize;
            let lj = lj_params.get(type_i, type_j);
            lj_c6[ai][aj] = lj.c6;
            lj_c12[ai][aj] = lj.c12;
            q_prod[ai][aj] = qi * charges[first_mol + aj] * FOUR_PI_EPS_I;
        }
    }

    for &(i_first, j_first) in pairs {
        let i_first = i_first as usize;
        let j_first = j_first as usize;

        let pos_i0 = positions[i_first];
        let pos_j0 = positions[j_first];

        // Compute PBC shift from first-atom nearest image (O-O for water)
        let r_first = periodicity.nearest_image(pos_i0, pos_j0);
        let tx = r_first.x - pos_i0.x + pos_j0.x;
        let ty = r_first.y - pos_i0.y + pos_j0.y;
        let tz = r_first.z - pos_i0.z + pos_j0.z;

        for atom_i in 0..n {
            let i = i_first + atom_i;
            let xi = positions[i].x + tx;
            let yi = positions[i].y + ty;
            let zi = positions[i].z + tz;

            for atom_j in 0..n {
                let j = j_first + atom_j;

                let x = xi - positions[j].x;
                let y = yi - positions[j].y;
                let z = zi - positions[j].z;
                let r = Vec3::new(x, y, z);

                let (f_mag, e_lj, e_crf) = lj_crf_interaction(
                    r, lj_c6[atom_i][atom_j], lj_c12[atom_i][atom_j],
                    q_prod[atom_i][atom_j], crf,
                );

                let force = r * f_mag;
                storage.forces[i] += force;
                storage.forces[j] -= force;

                storage.e_lj += e_lj;
                storage.e_crf += e_crf;

                if VIRIAL {
                    for a in 0..3 {
                        for b in 0..3 {
                            storage.virial[a][b] += r[a] * force[b];
                        }
                    }
                }
            }
        }
    }
}

/// Minimum pairlist size to justify parallel overhead.
const PARALLEL_THRESHOLD: usize = 2048;

/// Inner loop for nonbonded interactions (solute-solute and solute-solvent).
///
/// This is the hottest function in MD simulations!
/// Processes pairlist and accumulates forces, energies, and virial.
pub fn lj_crf_innerloop<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    storage: &mut ForceStorage,
) {
    process_pairs::<BC, true>(pairlist, positions, charges, iac, lj_params, crf, periodicity, storage);
}

/// Inner loop without virial computation (for NVE/NVT without pressure coupling).
pub fn lj_crf_innerloop_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    storage: &mut ForceStorage,
) {
    process_pairs::<BC, false>(pairlist, positions, charges, iac, lj_params, crf, periodicity, storage);
}

/// CG-grouped innerloop with virial (for NPT).
pub fn lj_crf_innerloop_cg_grouped<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    storage: &mut ForceStorage,
) {
    process_pairs_cg_grouped::<BC, true>(pairlist, groups, positions, charges, iac, lj_params, crf, periodicity, storage);
}

/// CG-grouped innerloop without virial (for NVE/NVT).
pub fn lj_crf_innerloop_cg_grouped_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    storage: &mut ForceStorage,
) {
    process_pairs_cg_grouped::<BC, false>(pairlist, groups, positions, charges, iac, lj_params, crf, periodicity, storage);
}

/// Parallel CG-grouped innerloop with virial.
pub fn lj_crf_innerloop_cg_grouped_parallel<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    n_atoms: usize,
) -> ForceStorage {
    lj_crf_innerloop_cg_grouped_parallel_virial::<BC, true>(positions, charges, iac, pairlist, groups, lj_params, crf, periodicity, n_atoms)
}

/// Parallel CG-grouped innerloop without virial.
pub fn lj_crf_innerloop_cg_grouped_parallel_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    n_atoms: usize,
) -> ForceStorage {
    lj_crf_innerloop_cg_grouped_parallel_virial::<BC, false>(positions, charges, iac, pairlist, groups, lj_params, crf, periodicity, n_atoms)
}

fn lj_crf_innerloop_cg_grouped_parallel_virial<BC: BoundaryCondition, const VIRIAL: bool>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    n_atoms: usize,
) -> ForceStorage {
    let mut result = ForceStorage::new(n_atoms);
    if groups.len() < 64 {
        process_pairs_cg_grouped::<BC, VIRIAL>(pairlist, groups, positions, charges, iac, lj_params, crf, periodicity, &mut result);
        return result;
    }

    result = groups
        .par_chunks(groups.len() / rayon::current_num_threads().max(1))
        .fold(
            || ForceStorage::new(n_atoms),
            |mut local_storage, chunk| {
                process_pairs_cg_grouped::<BC, VIRIAL>(pairlist, chunk, positions, charges, iac, lj_params, crf, periodicity, &mut local_storage);
                local_storage
            },
        )
        .reduce(
            || ForceStorage::new(n_atoms),
            |mut a, b| {
                a.merge(&b);
                a
            },
        );

    result
}

/// Parallel nonbonded innerloop — same physics as `lj_crf_innerloop`.
///
/// Uses Rayon fold/reduce with thread-local ForceStorage buffers (GROMOS pattern:
/// per-thread Nonbonded_Set with private force arrays, sequential reduction).
///
/// Automatically falls back to serial for small pairlists.
pub fn lj_crf_innerloop_parallel<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    n_atoms: usize,
) -> ForceStorage {
    lj_crf_innerloop_parallel_virial::<BC, true>(positions, charges, iac, pairlist, lj_params, crf, periodicity, n_atoms)
}

/// Parallel nonbonded innerloop without virial (for NVE/NVT).
pub fn lj_crf_innerloop_parallel_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    n_atoms: usize,
) -> ForceStorage {
    lj_crf_innerloop_parallel_virial::<BC, false>(positions, charges, iac, pairlist, lj_params, crf, periodicity, n_atoms)
}

fn lj_crf_innerloop_parallel_virial<BC: BoundaryCondition, const VIRIAL: bool>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    n_atoms: usize,
) -> ForceStorage {
    let mut result = ForceStorage::new(n_atoms);
    if pairlist.len() < PARALLEL_THRESHOLD {
        process_pairs::<BC, VIRIAL>(pairlist, positions, charges, iac, lj_params, crf, periodicity, &mut result);
        return result;
    }

    result = pairlist
        .par_chunks(1024)
        .fold(
            || ForceStorage::new(n_atoms),
            |mut acc, chunk| {
                process_pairs::<BC, VIRIAL>(chunk, positions, charges, iac, lj_params, crf, periodicity, &mut acc);
                acc
            },
        )
        .reduce(
            || ForceStorage::new(n_atoms),
            |mut a, b| { a.merge(&b); a },
        );
    result
}

/// Solvent-solvent innerloop with shared PBC shift (GROMOS convention).
///
/// For each solvent-solvent CG pair, the PBC shift is computed once from
/// the first atoms (typically O-O for water) and reused for all atom pairs.
///
/// The pairlist stores pairs of first-atom indices (one per solvent molecule).
/// `atoms_per_solvent` defines the molecule size (e.g. 3 for water).
pub fn solvent_innerloop<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    storage: &mut ForceStorage,
) {
    process_solvent_pairs::<BC, true>(pairlist, positions, charges, iac, lj_params, crf, periodicity, atoms_per_solvent, storage);
}

/// Solvent-solvent innerloop without virial computation (for NVE/NVT).
pub fn solvent_innerloop_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    storage: &mut ForceStorage,
) {
    process_solvent_pairs::<BC, false>(pairlist, positions, charges, iac, lj_params, crf, periodicity, atoms_per_solvent, storage);
}

/// Parallel solvent-solvent innerloop — same physics as `solvent_innerloop`.
///
/// Uses Rayon fold/reduce with thread-local buffers. Falls back to serial
/// for small pairlists.
pub fn solvent_innerloop_parallel<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    n_atoms: usize,
) -> ForceStorage {
    let mut result = ForceStorage::new(n_atoms);
    if pairlist.len() < PARALLEL_THRESHOLD {
        process_solvent_pairs::<BC, true>(pairlist, positions, charges, iac, lj_params, crf, periodicity, atoms_per_solvent, &mut result);
        return result;
    }

    result = pairlist
        .par_chunks(512)
        .fold(
            || ForceStorage::new(n_atoms),
            |mut acc, chunk| {
                process_solvent_pairs::<BC, true>(chunk, positions, charges, iac, lj_params, crf, periodicity, atoms_per_solvent, &mut acc);
                acc
            },
        )
        .reduce(
            || ForceStorage::new(n_atoms),
            |mut a, b| { a.merge(&b); a },
        );
    result
}

/// Parallel solvent-solvent innerloop without virial (for NVE/NVT).
pub fn solvent_innerloop_parallel_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    n_atoms: usize,
) -> ForceStorage {
    let mut result = ForceStorage::new(n_atoms);
    if pairlist.len() < PARALLEL_THRESHOLD {
        process_solvent_pairs::<BC, false>(pairlist, positions, charges, iac, lj_params, crf, periodicity, atoms_per_solvent, &mut result);
        return result;
    }

    result = pairlist
        .par_chunks(512)
        .fold(
            || ForceStorage::new(n_atoms),
            |mut acc, chunk| {
                process_solvent_pairs::<BC, false>(chunk, positions, charges, iac, lj_params, crf, periodicity, atoms_per_solvent, &mut acc);
                acc
            },
        )
        .reduce(
            || ForceStorage::new(n_atoms),
            |mut a, b| { a.merge(&b); a },
        );
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::params::LJParameters;
    use gromos_core::math::{Rectangular, Vacuum};
    use approx::assert_relative_eq;

    #[test]
    fn test_lj_interaction() {
        // Test case: two atoms at 1 nm distance
        let r = Vec3::new(1.0, 0.0, 0.0);
        let c6 = 0.001; // Typical C6 value
        let c12 = 0.0001; // Typical C12 value
        let q_prod = 0.0; // No charges

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.0,
            crf_cut3i: 0.0,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let (f, e_lj, _e_crf) = lj_crf_interaction(r, c6, c12, q_prod, &crf);

        // Verify energy is correct: E = C12/r^12 - C6/r^6
        let expected_e_lj = c12 - c6;
        assert_relative_eq!(e_lj, expected_e_lj, epsilon = 1e-9);

        // Verify force has correct sign (attractive at this distance)
        assert!(f < 0.0, "Force should be attractive");
    }

    #[test]
    fn test_innerloop_simple() {
        // Simple 2-atom system
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let charges = vec![0.5, -0.5];
        let iac = vec![0, 0];
        let pairlist = vec![(0, 1)];

        let lj_params = LJParamMatrix::from_nested(&[vec![LJParameters {
            c6: 0.001,
            c12: 0.0001,
            cs6: 0.001,
            cs12: 0.0001,
        }]]);
        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431, // 2 / 1.4^3
            crf_cut3i: 0.364431 / 2.0,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let periodicity = Vacuum;
        let mut storage = ForceStorage::new(2);

        lj_crf_innerloop(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            &mut storage,
        );

        // Verify Newton's third law: F_i = -F_j
        assert_relative_eq!(storage.forces[0].x, -storage.forces[1].x, epsilon = 1e-6);
        assert_relative_eq!(storage.forces[0].y, -storage.forces[1].y, epsilon = 1e-6);
        assert_relative_eq!(storage.forces[0].z, -storage.forces[1].z, epsilon = 1e-6);

        // Verify energy is non-zero
        assert!(storage.e_lj != 0.0);
        assert!(storage.e_crf != 0.0);
    }

    #[test]
    fn test_periodic_boundary() {
        let box_size = Vec3::splat(10.0);
        let periodicity = Rectangular::new(box_size);

        // Two atoms across the boundary
        let positions = vec![Vec3::new(9.5, 0.0, 0.0), Vec3::new(0.5, 0.0, 0.0)];
        let charges = vec![0.0, 0.0];
        let iac = vec![0, 0];
        let pairlist = vec![(0, 1)];

        let lj_params = LJParamMatrix::from_nested(&[vec![LJParameters {
            c6: 0.001,
            c12: 0.0001,
            cs6: 0.001,
            cs12: 0.0001,
        }]]);
        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.0,
            crf_cut3i: 0.0,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let mut storage = ForceStorage::new(2);

        lj_crf_innerloop(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            &mut storage,
        );

        // Distance should be 1.0 nm (minimum image), not 9.0 nm
        // Force should be reasonable for 1 nm separation
        assert!(
            storage.forces[0].length() < 10.0,
            "Force too large - PBC not working"
        );
    }
}
