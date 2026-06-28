//! PBC molecule gathering (unwrapping).
//!
//! Faithful port of gromos-rs `bound/Boundary.cc` gather methods.
//! Operates on a flat position slice using `Periodicity::nearest_image`.
//!
//! # Convention
//! `Periodicity::nearest_image(ri, rj)` returns `ri − rj` under the minimum-image
//! convention (a *displacement* vector, not a position).  To place atom `ri` in the
//! image nearest to anchor `rj`:
//! ```text
//! gathered = rj + nearest_image(ri, rj)
//! ```
//!
//! # Methods
//! - [`gather_chain`]  — sequential: atom 0 → reference, then k → k-1.
//!   Fast; correct for linear/ring molecules whose atoms appear in bond order.
//! - [`gather_bond`]   — bond-connectivity BFS: propagates from the first gathered
//!   atom through the bond graph.  Correct for branched or badly-ordered topologies.
//! - [`gather_molecules`] — apply `gather_chain` to every molecule range, updating
//!   the reference after each molecule.  Mirrors gromos-rs `Boundary::gather()`.

use crate::math::{Periodicity, Vec3};

/// Chain-gather a single molecule's atoms in-place.
///
/// Atom at `mol_atoms[0]` is placed nearest to `reference`.
/// Each subsequent atom is placed nearest to the previous atom (chain).
pub fn gather_chain(
    positions: &mut [Vec3],
    mol_atoms: &[usize],
    periodicity: &Periodicity,
    reference: Vec3,
) {
    if mol_atoms.is_empty() {
        return;
    }

    let p0 = positions[mol_atoms[0]];
    positions[mol_atoms[0]] = reference + periodicity.nearest_image(p0, reference);

    for w in mol_atoms.windows(2) {
        let anchor = positions[w[0]];
        let pi = positions[w[1]];
        positions[w[1]] = anchor + periodicity.nearest_image(pi, anchor);
    }
}

/// Bond-connectivity gathering for a single molecule.
///
/// Starts with `mol_atoms[0]` gathered w.r.t. `reference`, then propagates
/// through the bond graph (BFS-like: iterates all bonds repeatedly until every
/// atom is gathered).  Correct for branched or disordered atom orderings.
///
/// `bonds` contains pairs of **local** (0-based within `mol_atoms`) atom indices.
/// Runs at most `mol_atoms.len()` passes — O(N·B).
pub fn gather_bond(
    positions: &mut [Vec3],
    mol_atoms: &[usize],
    bonds: &[(usize, usize)],
    periodicity: &Periodicity,
    reference: Vec3,
) {
    if mol_atoms.is_empty() {
        return;
    }
    let n = mol_atoms.len();
    let mut gathered = vec![false; n];

    let p0 = positions[mol_atoms[0]];
    positions[mol_atoms[0]] = reference + periodicity.nearest_image(p0, reference);
    gathered[0] = true;

    let mut remaining = n - 1;
    while remaining > 0 {
        let mut progress = false;
        for &(li, lj) in bonds {
            if li >= n || lj >= n {
                continue;
            }
            if gathered[li] && !gathered[lj] {
                let anchor = positions[mol_atoms[li]];
                let pi = positions[mol_atoms[lj]];
                positions[mol_atoms[lj]] = anchor + periodicity.nearest_image(pi, anchor);
                gathered[lj] = true;
                remaining -= 1;
                progress = true;
            } else if gathered[lj] && !gathered[li] {
                let anchor = positions[mol_atoms[lj]];
                let pi = positions[mol_atoms[li]];
                positions[mol_atoms[li]] = anchor + periodicity.nearest_image(pi, anchor);
                gathered[li] = true;
                remaining -= 1;
                progress = true;
            }
        }
        if !progress {
            break;
        }
    }
}

/// Gather all molecules using the chain method.
///
/// Each molecule's first atom is gathered w.r.t. the centre-of-geometry (COG) of the
/// previously gathered atoms.  This matches gromos-rs `Boundary::gather()` where
/// subsequent molecules are gathered towards the COG of the already-placed system.
///
/// `mol_ranges` — slice of `Topology::molecules` ranges (0-based global atom indices).
/// `reference`  — initial reference point (e.g. `Vec3::ZERO` or first-atom position).
pub fn gather_molecules(
    positions: &mut [Vec3],
    mol_ranges: &[std::ops::Range<usize>],
    periodicity: &Periodicity,
    reference: Vec3,
) {
    if mol_ranges.is_empty() {
        return;
    }
    let mut ref_pt = reference;
    let mut total_gathered = 0usize;
    let mut cog = Vec3::ZERO;

    for mol_range in mol_ranges {
        let atoms: Vec<usize> = mol_range.clone().collect();
        gather_chain(positions, &atoms, periodicity, ref_pt);

        // Update COG for next molecule's reference
        for &idx in &atoms {
            cog += positions[idx];
        }
        total_gathered += atoms.len();
        if total_gathered > 0 {
            ref_pt = cog / total_gathered as f64;
        }
    }
}

/// Centre of geometry of a set of atoms.
pub fn centre_of_geometry(positions: &[Vec3], indices: &[usize]) -> Vec3 {
    if indices.is_empty() {
        return Vec3::ZERO;
    }
    indices
        .iter()
        .map(|&i| positions[i])
        .fold(Vec3::ZERO, |a, b| a + b)
        / indices.len() as f64
}

/// Centre of mass of a set of atoms.
pub fn centre_of_mass(positions: &[Vec3], masses: &[f64], indices: &[usize]) -> Vec3 {
    if indices.is_empty() {
        return Vec3::ZERO;
    }
    let (sum, total_mass) = indices.iter().fold((Vec3::ZERO, 0.0f64), |(s, m), &i| {
        (s + positions[i] * masses[i], m + masses[i])
    });
    if total_mass > 0.0 {
        sum / total_mass
    } else {
        Vec3::ZERO
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Periodicity, Rectangular, Vec3};

    fn rect(l: f64) -> Periodicity {
        Periodicity::Rectangular(Rectangular::new(Vec3::new(l, l, l)))
    }

    #[test]
    fn gather_chain_wraps_split_molecule() {
        // Two atoms split across box boundary
        let mut pos = vec![Vec3::new(0.1, 0.0, 0.0), Vec3::new(0.9, 0.0, 0.0)];
        gather_chain(&mut pos, &[0, 1], &rect(1.0), Vec3::ZERO);
        // atom 1 should move to -0.1 (nearest image of 0.9 w.r.t. 0.1)
        assert!((pos[1].x - (-0.1)).abs() < 1e-10, "pos[1].x = {}", pos[1].x);
    }

    #[test]
    fn gather_chain_no_wrap_needed() {
        let mut pos = vec![Vec3::new(0.3, 0.0, 0.0), Vec3::new(0.5, 0.0, 0.0)];
        let orig = pos.clone();
        gather_chain(&mut pos, &[0, 1], &rect(1.0), Vec3::ZERO);
        assert!((pos[0] - orig[0]).length() < 1e-10);
        assert!((pos[1] - orig[1]).length() < 1e-10);
    }

    #[test]
    fn gather_bond_branched_molecule() {
        // 3-atom branched: bonds (0-1) and (0-2); atom 2 is across boundary from atom 0
        let mut pos = vec![
            Vec3::new(0.4, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.95, 0.0, 0.0),
        ];
        let bonds = vec![(0, 1), (0, 2)];
        gather_bond(&mut pos, &[0, 1, 2], &bonds, &rect(1.0), Vec3::ZERO);
        // atom 2 nearest to atom 0 (0.4): 0.95 - 1.0 = -0.05
        assert!(
            (pos[2].x - (-0.05)).abs() < 1e-10,
            "pos[2].x = {}",
            pos[2].x
        );
        // atom 1 stays (0.5 is nearest to 0.4)
        assert!((pos[1].x - 0.5).abs() < 1e-10);
    }

    #[test]
    fn gather_molecules_two_mols() {
        // mol0: atoms 0-1 near origin; mol1: atoms 2-3 split across boundary
        let mut pos = vec![
            Vec3::new(0.2, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(0.8, 0.0, 0.0),
            Vec3::new(0.05, 0.0, 0.0), // split: should gather near 0.8
        ];
        gather_molecules(&mut pos, &[0..2, 2..4], &rect(1.0), Vec3::ZERO);
        // mol1: atom 2 stays ~0.8; atom 3 gathers near atom 2 → ~1.05 - 1.0 = 0.05 + 1.0?
        // nearest_image(0.05, 0.8) = 0.05 - 0.8 - round(-0.75)*1 = -0.75 + 1.0 = 0.25 → gathered = 0.8+0.25=1.05
        let diff = (pos[3] - pos[2]).length();
        assert!(diff < 0.5, "gathered atoms within half-box: diff={diff:.4}");
    }

    #[test]
    fn cog_simple() {
        let pos = vec![Vec3::new(1.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0)];
        let cog = centre_of_geometry(&pos, &[0, 1]);
        assert!((cog.x - 2.0).abs() < 1e-10);
    }
}
