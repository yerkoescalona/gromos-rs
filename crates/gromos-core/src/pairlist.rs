//! Pairlist (neighbor list) generation.
//!
//! A pairlist pre-computes which atom pairs are within the nonbonded cutoff
//! so the inner force loop only evaluates close pairs, not all N² combinations.
//!
//! # GROMOS charge-group cutoff convention
//!
//! gromos-rs uses **charge groups** as the cutoff unit, matching gromosXX.
//! A charge group is a set of covalently-bonded atoms whose partial charges sum
//! to zero. The cutoff test is applied to charge-group centers-of-geometry (COG),
//! not individual atoms. This keeps neutral charge groups wholly inside or outside
//! the cutoff sphere, which is required for reaction-field electrostatics to be
//! correct: splitting a charge group across the cutoff boundary leaves a net-charge
//! imbalance that the RF formula is not designed to handle.
//!
//! # Twin-range cutoffs
//!
//! Two radii are supported: `short_range_cutoff` (`RCUTP`) and `long_range_cutoff`
//! (`RCUTL`). Pairs within the short range are evaluated every step; pairs in
//! `(RCUTP, RCUTL]` are evaluated only on pairlist-update steps and their forces
//! are reused between updates. When `RCUTP == RCUTL` the long lists are always empty.
//!
//! # Algorithms
//!
//! Two implementations are available. Both produce **identical pair sets**; the
//! difference is runtime complexity:
//!
//! - [`StandardPairlistAlgorithm`] — O(N²) over charge-group pairs; always correct;
//!   used as the reference oracle in tests.
//! - [`CellListPairlistAlgorithm`] — O(N) for rectangular boxes via a spatial grid;
//!   falls back to [`StandardPairlistAlgorithm`] for triclinic/vacuum boxes.
//!
//! See the [pairlist algorithm chapter](https://gromos-rs.github.io/algorithms/pairlist.html)
//! in the gromos-rs book for the scientific rationale.
//!
//! # Source references
//!
//! - `md++/src/interaction/nonbonded/pairlist/standard_pairlist_algorithm.cc`
//! - `md++/src/interaction/nonbonded/pairlist/grid_cell_pairlist.cc`

use crate::configuration::{BoxType, Configuration};
use crate::math::{BoundaryCondition, Vec3};
use crate::topology::Topology;
use rayon::prelude::*;
use std::collections::HashSet;

/// Ordered list of atom-index pairs `(i, j)` with `i < j` within the cutoff.
pub type Pairlist = Vec<(usize, usize)>;

/// Holds all nonbonded atom pairs, split by range and solute/solvent role.
///
/// The four lists correspond to the twin-range + solute/solvent classification
/// used by gromosXX. When `short_range_cutoff == long_range_cutoff` the `*_long`
/// lists are always empty and all forces are evaluated every step.
///
/// The `skin` extra distance is added to both cutoffs when building the lists so
/// that pairs which enter the true cutoff before the next rebuild are already
/// present. Rebuild is triggered every `update_frequency` steps (IMD `NSNB`).
#[derive(Debug, Clone)]
pub struct PairlistContainer {
    /// Solute–solute pairs within `short_range_cutoff` (1-4 and excluded pairs removed).
    pub solute_short: Pairlist,
    /// Solute–solute pairs in `(short_range_cutoff, long_range_cutoff]` (twin-range).
    pub solute_long: Pairlist,
    /// Solute–solvent and solvent–solvent pairs within `short_range_cutoff`.
    /// For solvent–solvent, stores one sentinel pair per molecule pair; the
    /// `solvent_innerloop` kernel expands it to all atom pairs within the molecule type.
    pub solvent_short: Pairlist,
    /// Solute–solvent and solvent–solvent pairs in `(short_range_cutoff, long_range_cutoff]`.
    pub solvent_long: Pairlist,

    /// Short-range cutoff radius in nm (IMD `RCUTP`).
    pub short_range_cutoff: f64,
    /// Long-range cutoff radius in nm (IMD `RCUTL`). Equal to `short_range_cutoff`
    /// when twin-range is disabled.
    pub long_range_cutoff: f64,
    /// Extra distance added to both cutoffs during list construction (nm).
    /// Ensures pairs that enter the true cutoff before the next rebuild are
    /// already present. See `needs_update` for the rebuild criterion.
    pub skin: f64,

    /// Steps elapsed since the last rebuild.
    pub steps_since_update: usize,
    /// Rebuild the pairlist every this many steps (IMD `NSNB`).
    pub update_frequency: usize,
}

impl PairlistContainer {
    /// Create a new container with the given cutoffs and skin distance (all in nm).
    pub fn new(short_cutoff: f64, long_cutoff: f64, skin: f64) -> Self {
        Self {
            solute_short: Vec::new(),
            solute_long: Vec::new(),
            solvent_short: Vec::new(),
            solvent_long: Vec::new(),
            short_range_cutoff: short_cutoff,
            long_range_cutoff: long_cutoff,
            skin,
            steps_since_update: 0,
            update_frequency: 5, // Default: update every 5 steps
        }
    }

    /// Returns `true` if the pairlist should be rebuilt this step.
    pub fn needs_update(&self) -> bool {
        self.steps_since_update >= self.update_frequency
    }

    /// Advance the step counter. Call once per MD step.
    pub fn step(&mut self) {
        self.steps_since_update += 1;
    }

    /// Reset the step counter after a rebuild. Called automatically by each algorithm's `update()`.
    pub fn reset_counter(&mut self) {
        self.steps_since_update = 0;
    }

    /// Total number of pairs across all four lists.
    pub fn total_pairs(&self) -> usize {
        self.solute_short.len()
            + self.solute_long.len()
            + self.solvent_short.len()
            + self.solvent_long.len()
    }
}

/// O(N²) pairlist algorithm — tests every charge-group pair.
///
/// Faithful translation of `standard_pairlist_algorithm.cc` from gromosXX.
/// Correct for all box types (rectangular, triclinic, vacuum) and all system
/// sizes. Kept as the reference oracle for validating [`CellListPairlistAlgorithm`].
///
/// For production runs with more than ~500 atoms in a rectangular box, prefer
/// [`CellListPairlistAlgorithm`] which scales O(N) instead of O(N²).
pub struct StandardPairlistAlgorithm {
    /// When `true`, uses charge-group COGs as the cutoff unit (GROMOS default).
    /// When `false`, uses individual atom positions (faster for small vacuum systems).
    pub use_chargegroups: bool,
}

impl StandardPairlistAlgorithm {
    /// Create a new algorithm instance.
    ///
    /// Pass `use_chargegroups: true` for all periodic-box simulations (required for
    /// correct reaction-field electrostatics). `false` is only appropriate for small
    /// vacuum systems without solvent.
    pub fn new(use_chargegroups: bool) -> Self {
        Self { use_chargegroups }
    }

    /// Update pairlist based on current configuration
    ///
    /// Translated from standard_pairlist_algorithm.cc:168-250 (chargegroup version)
    pub fn update<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        pairlist.solute_short.clear();
        pairlist.solute_long.clear();
        pairlist.solvent_short.clear();
        pairlist.solvent_long.clear();

        if self.use_chargegroups {
            self.update_chargegroup_based(topo, conf, pairlist, periodicity);
        } else {
            self.update_atom_based(topo, conf, pairlist, periodicity);
        }

        pairlist.reset_counter();
    }

    /// Chargegroup-based pairlist generation (GROMOS default)
    ///
    /// Matches GROMOS standard_pairlist_algorithm.cc:
    /// - Solute CGs: center-of-geometry distance, exclusion checks
    /// - Solvent CGs: first-atom position distance, no exclusion checks
    /// - Intra-CG non-excluded pairs added for solute CGs
    /// - Solvent-solvent pairs go to solvent_short/long, rest to solute_short/long
    fn update_chargegroup_based<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        let n_chargegroups = topo.chargegroups.len();
        let n_solute_cg = topo.num_solute_chargegroups;
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);

        // Precompute solute CG centers-of-geometry (GROMOS only computes COG for solute CGs)
        let cg_cog: Vec<Vec3> = (0..n_solute_cg)
            .map(|i| topo.chargegroups[i].center_of_geometry(&conf.current().pos))
            .collect();

        log::debug!("Pairlist update: {} CGs total, {} solute CGs, cutoff_short²={:.6}, cutoff_long²={:.6}",
            n_chargegroups, n_solute_cg, cutoff2_short, cutoff2_long);
        for i in 0..n_solute_cg {
            log::debug!("  Solute CG {}: atoms={:?}, COG=({:.6},{:.6},{:.6})",
                i, topo.chargegroups[i].atoms, cg_cog[i].x, cg_cog[i].y, cg_cog[i].z);
        }
        for i in n_solute_cg..n_chargegroups {
            let fa = topo.chargegroups[i].atoms[0];
            log::debug!("  Solvent CG {}: first_atom={}, pos=({:.6},{:.6},{:.6}), atoms={:?}",
                i, fa, conf.current().pos[fa].x, conf.current().pos[fa].y, conf.current().pos[fa].z,
                topo.chargegroups[i].atoms);
        }

        // === Solute CGs ===
        let mut n_solute_solute_cg_pairs = 0usize;
        let mut n_solute_solute_skip = 0usize;
        let mut n_solute_solvent_cg_pairs = 0usize;
        let mut n_solute_solvent_skip = 0usize;
        let mut n_intra_cg_pairs = 0usize;

        for cg1 in 0..n_solute_cg {
            // Intra-CG non-excluded pairs (GROMOS lines ~185-200)
            let cg1_atoms = &topo.chargegroups[cg1].atoms;
            for ai in 0..cg1_atoms.len() {
                let a1 = cg1_atoms[ai];
                for aj in (ai + 1)..cg1_atoms.len() {
                    let a2 = cg1_atoms[aj];
                    if !topo.is_excluded_or_14(a1, a2) {
                        pairlist.solute_short.push((a1, a2));
                        n_intra_cg_pairs += 1;
                    }
                }
            }

            // Solute-solute pairs
            for cg2 in (cg1 + 1)..n_solute_cg {
                let r = periodicity.nearest_image(cg_cog[cg1], cg_cog[cg2]);
                let dist2 = r.length_squared();

                if dist2 > cutoff2_long {
                    log::debug!("  Solute-Solute CG({},{}) SKIP dist={:.6} > long_cut={:.6}",
                        cg1, cg2, dist2.sqrt(), cutoff2_long.sqrt());
                    n_solute_solute_skip += 1;
                    continue;
                }

                n_solute_solute_cg_pairs += 1;
                if dist2 > cutoff2_short {
                    log::debug!("  Solute-Solute CG({},{}) LONG dist={:.6}", cg1, cg2, dist2.sqrt());
                    // Long-range (no exclusion check in GROMOS for long-range solute-solute)
                    for &a1 in &topo.chargegroups[cg1].atoms {
                        for &a2 in &topo.chargegroups[cg2].atoms {
                            pairlist.solute_long.push((a1, a2));
                        }
                    }
                    continue;
                }

                log::debug!("  Solute-Solute CG({},{}) SHORT dist={:.6}", cg1, cg2, dist2.sqrt());
                // Short-range with exclusion check
                for &a1 in &topo.chargegroups[cg1].atoms {
                    for &a2 in &topo.chargegroups[cg2].atoms {
                        if !topo.is_excluded_or_14(a1, a2) {
                            pairlist.solute_short.push((a1, a2));
                        } else {
                            log::debug!("    Excluded pair ({},{})", a1, a2);
                        }
                    }
                }
            }

            // Solute-solvent pairs (solvent CG distance uses first-atom position)
            for cg2 in n_solute_cg..n_chargegroups {
                let solvent_first_atom = topo.chargegroups[cg2].atoms[0];
                let r = periodicity.nearest_image(cg_cog[cg1], conf.current().pos[solvent_first_atom]);
                let dist2 = r.length_squared();

                if dist2 > cutoff2_long {
                    log::debug!("  Solute({})-Solvent({}) SKIP dist={:.6} > long_cut={:.6}",
                        cg1, cg2, dist2.sqrt(), cutoff2_long.sqrt());
                    n_solute_solvent_skip += 1;
                    continue;
                }

                n_solute_solvent_cg_pairs += 1;
                if dist2 > cutoff2_short {
                    log::debug!("  Solute({})-Solvent({}) LONG dist={:.6}", cg1, cg2, dist2.sqrt());
                    for &a1 in &topo.chargegroups[cg1].atoms {
                        for &a2 in &topo.chargegroups[cg2].atoms {
                            pairlist.solute_long.push((a1, a2));
                        }
                    }
                    continue;
                }

                log::debug!("  Solute({})-Solvent({}) SHORT dist={:.6}", cg1, cg2, dist2.sqrt());
                // Short-range: no exclusion check for solute-solvent (GROMOS convention)
                for &a1 in &topo.chargegroups[cg1].atoms {
                    for &a2 in &topo.chargegroups[cg2].atoms {
                        pairlist.solute_short.push((a1, a2));
                    }
                }
            }
        }

        // === Solvent-solvent CGs ===
        // Distance uses first-atom positions (GROMOS _solvent_solvent)
        // Store only first-atom pairs; solvent_innerloop expands all atom pairs internally
        let mut n_solvent_solvent_cg_pairs = 0usize;
        let mut n_solvent_solvent_skip = 0usize;

        for cg1 in n_solute_cg..n_chargegroups {
            let first_atom_1 = topo.chargegroups[cg1].atoms[0];

            for cg2 in (cg1 + 1)..n_chargegroups {
                let first_atom_2 = topo.chargegroups[cg2].atoms[0];
                let r = periodicity.nearest_image(
                    conf.current().pos[first_atom_1],
                    conf.current().pos[first_atom_2],
                );
                let dist2 = r.length_squared();

                if dist2 > cutoff2_long {
                    log::debug!("  Solvent({},a{})-Solvent({},a{}) SKIP dist={:.6} > long_cut={:.6}",
                        cg1, first_atom_1, cg2, first_atom_2, dist2.sqrt(), cutoff2_long.sqrt());
                    n_solvent_solvent_skip += 1;
                    continue;
                }

                n_solvent_solvent_cg_pairs += 1;
                if dist2 > cutoff2_short {
                    log::debug!("  Solvent({},a{})-Solvent({},a{}) LONG dist={:.6}",
                        cg1, first_atom_1, cg2, first_atom_2, dist2.sqrt());
                    // Long-range: expand to all atom pairs (GROMOS standard
                    // mode uses per-atom nearest_image, not shared PBC shift)
                    for &a1 in &topo.chargegroups[cg1].atoms {
                        for &a2 in &topo.chargegroups[cg2].atoms {
                            pairlist.solvent_long.push((a1, a2));
                        }
                    }
                    continue;
                }

                log::debug!("  Solvent({},a{})-Solvent({},a{}) SHORT dist={:.6}",
                    cg1, first_atom_1, cg2, first_atom_2, dist2.sqrt());
                // Short-range solvent-solvent: store first-atom pair only
                pairlist.solvent_short.push((first_atom_1, first_atom_2));
            }
        }

        log::debug!("CG pairlist summary:");
        log::debug!("  Solute intra-CG pairs: {}", n_intra_cg_pairs);
        log::debug!("  Solute-Solute CG pairs: {} included, {} skipped", n_solute_solute_cg_pairs, n_solute_solute_skip);
        log::debug!("  Solute-Solvent CG pairs: {} included, {} skipped", n_solute_solvent_cg_pairs, n_solute_solvent_skip);
        log::debug!("  Solvent-Solvent CG pairs: {} included, {} skipped", n_solvent_solvent_cg_pairs, n_solvent_solvent_skip);
        log::debug!("  solute_short={}, solute_long={}, solvent_short={}, solvent_long={}, cutoff={:.4}",
            pairlist.solute_short.len(), pairlist.solute_long.len(),
            pairlist.solvent_short.len(), pairlist.solvent_long.len(),
            pairlist.short_range_cutoff);
    }

    /// Atom-based pairlist generation (simpler, no chargegroups)
    fn update_atom_based<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        let n_atoms = topo.num_atoms();
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);

        for i in 0..n_atoms {
            for j in (i + 1)..n_atoms {
                // Check exclusions
                if topo.is_excluded_or_14(i, j) {
                    continue;
                }

                // Calculate distance
                let r = periodicity.nearest_image(conf.current().pos[i], conf.current().pos[j]);
                let dist2 = r.length_squared();

                // Add to appropriate pairlist
                if dist2 < cutoff2_short {
                    pairlist.solute_short.push((i, j));
                } else if dist2 < cutoff2_long {
                    pairlist.solute_long.push((i, j));
                }
            }
        }
    }
}

/// Charge-group-aware linked-cell pairlist algorithm — an O(N) alternative
/// to [`StandardPairlistAlgorithm`] for rectangular boxes (FUTURE.md Dim 9a).
///
/// Bins chargegroup reference positions into a grid that exactly tiles the
/// box, with each cell at least `long_range_cutoff + skin` wide along every
/// axis, then only tests CG pairs that share a cell or occupy a
/// periodic-wrapped neighboring cell. The reference position and the
/// pairwise distance/exclusion logic for each of the three cases
/// (solute-solute, solute-solvent, solvent-solvent) mirror
/// [`StandardPairlistAlgorithm`]'s chargegroup-based algorithm exactly, so
/// the resulting pair *set* is identical (see
/// `tests::test_cell_list_matches_standard_*`), just in a different
/// (cell-traversal) order.
///
/// For non-rectangular boxes (vacuum, triclinic, truncated octahedron) an
/// axis-aligned grid can't be made periodicity-safe without extra
/// machinery, so this falls back to `StandardPairlistAlgorithm`'s O(N²)
/// path — still correct, just not accelerated.
#[derive(Debug, Clone, Copy, Default)]
pub struct CellListPairlistAlgorithm;

impl CellListPairlistAlgorithm {
    /// Create a new cell-list algorithm instance.
    pub fn new() -> Self {
        Self
    }

    /// Update pairlist based on current configuration.
    ///
    /// Produces the same `solute_short`/`solute_long`/`solvent_short`/
    /// `solvent_long` pair sets as `StandardPairlistAlgorithm::new(true)`.
    pub fn update<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        pairlist.solute_short.clear();
        pairlist.solute_long.clear();
        pairlist.solvent_short.clear();
        pairlist.solvent_long.clear();

        if conf.current().box_config.box_type == BoxType::Rectangular {
            self.update_cell_list(topo, conf, pairlist, periodicity);
        } else {
            // Axis-aligned cells aren't periodicity-safe for vacuum/triclinic
            // boxes; fall back to the always-correct O(N^2) path.
            StandardPairlistAlgorithm::new(true)
                .update_chargegroup_based(topo, conf, pairlist, periodicity);
        }

        pairlist.reset_counter();
    }

    fn update_cell_list<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        let n_chargegroups = topo.chargegroups.len();
        if n_chargegroups == 0 {
            return;
        }
        let n_solute_cg = topo.num_solute_chargegroups;
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);
        let min_cell_size = pairlist.long_range_cutoff + pairlist.skin;

        let positions = &conf.current().pos;

        // Solute CG centers-of-geometry (StandardPairlistAlgorithm convention:
        // only solute CGs get a COG).
        let cg_cog: Vec<Vec3> = (0..n_solute_cg)
            .map(|i| topo.chargegroups[i].center_of_geometry(positions))
            .collect();

        // Intra-CG non-excluded pairs don't depend on the spatial grid.
        for cg in 0..n_solute_cg {
            let atoms = &topo.chargegroups[cg].atoms;
            for ai in 0..atoms.len() {
                let a1 = atoms[ai];
                for &a2 in &atoms[ai + 1..] {
                    if !topo.is_excluded_or_14(a1, a2) {
                        pairlist.solute_short.push((a1, a2));
                    }
                }
            }
        }

        // Reference position for binning AND for the CG-CG distance test:
        // solute CGs use their COG, solvent CGs use their first atom
        // (matches StandardPairlistAlgorithm).
        let ref_pos: Vec<Vec3> = (0..n_chargegroups)
            .map(|cg| {
                if cg < n_solute_cg {
                    cg_cog[cg]
                } else {
                    positions[topo.chargegroups[cg].atoms[0]]
                }
            })
            .collect();

        // A grid that exactly tiles the box, with cells >= min_cell_size
        // along each axis: any pair within cutoff is in the same or an
        // adjacent cell, modulo periodic wrap of cell indices.
        let box_dims = conf.current().box_config.dimensions();
        let grid_dim = [
            grid_dim_for_axis(box_dims.x, min_cell_size),
            grid_dim_for_axis(box_dims.y, min_cell_size),
            grid_dim_for_axis(box_dims.z, min_cell_size),
        ];
        let cell_size = [
            box_dims.x / grid_dim[0] as f64,
            box_dims.y / grid_dim[1] as f64,
            box_dims.z / grid_dim[2] as f64,
        ];

        let total_cells = grid_dim[0] * grid_dim[1] * grid_dim[2];
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); total_cells];
        for (cg, &pos) in ref_pos.iter().enumerate() {
            cells[cell_index(pos, &cell_size, &grid_dim)].push(cg);
        }

        for cell_a in 0..total_cells {
            if cells[cell_a].is_empty() {
                continue;
            }
            for cell_b in neighbor_cells(cell_a, &grid_dim) {
                if cell_b < cell_a || cells[cell_b].is_empty() {
                    continue;
                }
                if cell_a == cell_b {
                    let members = &cells[cell_a];
                    for i in 0..members.len() {
                        for &cg2 in &members[i + 1..] {
                            process_cg_pair(
                                members[i], cg2, topo, positions, periodicity, pairlist,
                                n_solute_cg, &cg_cog, cutoff2_short, cutoff2_long,
                            );
                        }
                    }
                } else {
                    for &cg1 in &cells[cell_a] {
                        for &cg2 in &cells[cell_b] {
                            process_cg_pair(
                                cg1, cg2, topo, positions, periodicity, pairlist,
                                n_solute_cg, &cg_cog, cutoff2_short, cutoff2_long,
                            );
                        }
                    }
                }
            }
        }
    }
}

/// Number of cells along one axis: at least 1, each cell at least `min_cell_size`.
fn grid_dim_for_axis(box_len: f64, min_cell_size: f64) -> usize {
    if min_cell_size <= 0.0 || box_len <= 0.0 {
        return 1;
    }
    ((box_len / min_cell_size).floor() as usize).max(1)
}

/// Cell index of a position, wrapping into `[0, grid_dim)` per axis.
fn cell_index(pos: Vec3, cell_size: &[f64; 3], grid_dim: &[usize; 3]) -> usize {
    let cx = ((pos.x / cell_size[0]).floor() as i64).rem_euclid(grid_dim[0] as i64) as usize;
    let cy = ((pos.y / cell_size[1]).floor() as i64).rem_euclid(grid_dim[1] as i64) as usize;
    let cz = ((pos.z / cell_size[2]).floor() as i64).rem_euclid(grid_dim[2] as i64) as usize;
    cx + cy * grid_dim[0] + cz * grid_dim[0] * grid_dim[1]
}

/// All cells (including `cell` itself) within one periodic-wrapped step of
/// `cell`, de-duplicated (matters when `grid_dim` is 1 or 2 along some axis).
fn neighbor_cells(cell: usize, grid_dim: &[usize; 3]) -> Vec<usize> {
    let cx = (cell % grid_dim[0]) as i64;
    let cy = ((cell / grid_dim[0]) % grid_dim[1]) as i64;
    let cz = (cell / (grid_dim[0] * grid_dim[1])) as i64;

    let mut neighbors = HashSet::with_capacity(27);
    for dx in -1..=1i64 {
        for dy in -1..=1i64 {
            for dz in -1..=1i64 {
                let nx = (cx + dx).rem_euclid(grid_dim[0] as i64) as usize;
                let ny = (cy + dy).rem_euclid(grid_dim[1] as i64) as usize;
                let nz = (cz + dz).rem_euclid(grid_dim[2] as i64) as usize;
                neighbors.insert(nx + ny * grid_dim[0] + nz * grid_dim[0] * grid_dim[1]);
            }
        }
    }
    neighbors.into_iter().collect()
}

/// Classify and test a single chargegroup pair, pushing atom pairs into
/// `pairlist` exactly as `StandardPairlistAlgorithm`'s chargegroup-based
/// algorithm would for the same pair.
#[allow(clippy::too_many_arguments)]
fn process_cg_pair<BC: BoundaryCondition>(
    cg_a: usize,
    cg_b: usize,
    topo: &Topology,
    positions: &[Vec3],
    periodicity: &BC,
    pairlist: &mut PairlistContainer,
    n_solute_cg: usize,
    cg_cog: &[Vec3],
    cutoff2_short: f64,
    cutoff2_long: f64,
) {
    // StandardPairlistAlgorithm's loops always visit CG pairs with cg1 < cg2.
    let cg1 = cg_a.min(cg_b);
    let cg2 = cg_a.max(cg_b);

    if cg2 < n_solute_cg {
        // Solute-solute: CG centers-of-geometry, exclusion check for short-range.
        let r = periodicity.nearest_image(cg_cog[cg1], cg_cog[cg2]);
        let dist2 = r.length_squared();
        if dist2 > cutoff2_long {
            return;
        }
        if dist2 > cutoff2_short {
            for &a1 in &topo.chargegroups[cg1].atoms {
                for &a2 in &topo.chargegroups[cg2].atoms {
                    pairlist.solute_long.push((a1, a2));
                }
            }
        } else {
            for &a1 in &topo.chargegroups[cg1].atoms {
                for &a2 in &topo.chargegroups[cg2].atoms {
                    if !topo.is_excluded_or_14(a1, a2) {
                        pairlist.solute_short.push((a1, a2));
                    }
                }
            }
        }
    } else if cg1 < n_solute_cg {
        // Solute-solvent: solute CG COG vs solvent CG first atom, no exclusion check.
        let solvent_first_atom = topo.chargegroups[cg2].atoms[0];
        let r = periodicity.nearest_image(cg_cog[cg1], positions[solvent_first_atom]);
        let dist2 = r.length_squared();
        if dist2 > cutoff2_long {
            return;
        }
        if dist2 > cutoff2_short {
            for &a1 in &topo.chargegroups[cg1].atoms {
                for &a2 in &topo.chargegroups[cg2].atoms {
                    pairlist.solute_long.push((a1, a2));
                }
            }
        } else {
            for &a1 in &topo.chargegroups[cg1].atoms {
                for &a2 in &topo.chargegroups[cg2].atoms {
                    pairlist.solute_short.push((a1, a2));
                }
            }
        }
    } else {
        // Solvent-solvent: first-atom positions; short-range stores the
        // first-atom pair only (solvent_innerloop expands all atom pairs).
        let first_atom_1 = topo.chargegroups[cg1].atoms[0];
        let first_atom_2 = topo.chargegroups[cg2].atoms[0];
        let r = periodicity.nearest_image(positions[first_atom_1], positions[first_atom_2]);
        let dist2 = r.length_squared();
        if dist2 > cutoff2_long {
            return;
        }
        if dist2 > cutoff2_short {
            for &a1 in &topo.chargegroups[cg1].atoms {
                for &a2 in &topo.chargegroups[cg2].atoms {
                    pairlist.solvent_long.push((a1, a2));
                }
            }
        } else {
            pairlist.solvent_short.push((first_atom_1, first_atom_2));
        }
    }
}

/// Dispatch enum for pairlist construction algorithms.
///
/// Owned by `Forcefield` and other integrators. The hot-path `update()` call
/// inlines to a direct dispatch with zero heap allocation (no trait objects).
///
/// # Selection
///
/// Use [`PairlistAlgorithm::from_imd`] to pick the algorithm from IMD
/// parameters. Manual construction is also supported for tests and advanced
/// use cases.
pub enum PairlistAlgorithm {
    /// O(N²) reference algorithm — always correct, used for all non-rectangular
    /// boxes and as the default for small systems.
    Standard(StandardPairlistAlgorithm),
    /// O(N) cell-list algorithm — rectangular boxes only; falls back to
    /// `Standard` internally for vacuum/triclinic input.
    CellList(CellListPairlistAlgorithm),
}

impl PairlistAlgorithm {
    /// Choose an algorithm from IMD PAIRLIST block parameters.
    ///
    /// Slot values match gromosXX `in_parameter.cc:1419-1422` exactly.
    /// Every value is an **explicit instruction** — matching gromosXX where
    /// `ALGORITHM standard` always means Standard, with no auto-heuristic.
    ///
    /// | `algorithm` | gromosXX class | gromos-rs result |
    /// |-------------|----------------|-----------------|
    /// | `0` / `"standard"` | `Standard_Pairlist_Algorithm` | `Standard` (forced) |
    /// | `1` / `"grid"`     | `Extended_Grid_Pairlist_Algorithm` | `Standard` (not yet ported; safe fallback) |
    /// | `2` / `"grid_cell"`| `Grid_Cell_Pairlist` (Heinz & Hünenberger 2004) | `CellList` (forced) |
    /// | any other          | — | `Standard` (safe fallback) |
    ///
    /// 9a-1 verified: CellList produces bit-identical energies to Standard on
    /// `water_216_box` (648 atoms, 100 steps, margin = 0.000e0). CellList is
    /// only activated when the input explicitly requests `ALGORITHM grid_cell`.
    pub fn from_imd(
        algorithm: i32,
        _n_atoms: usize,
        _box_type: BoxType,
        has_chargegroups: bool,
    ) -> Self {
        match algorithm {
            2 => Self::CellList(CellListPairlistAlgorithm::new()),
            _ => Self::Standard(StandardPairlistAlgorithm::new(has_chargegroups)),
        }
    }

    /// Update the pairlist, delegating to the selected algorithm.
    pub fn update<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        match self {
            Self::Standard(a) => a.update(topo, conf, pairlist, periodicity),
            Self::CellList(a) => a.update(topo, conf, pairlist, periodicity),
        }
    }
}

/// Parallel pairlist generation using Rayon
pub struct ParallelPairlistAlgorithm {
    base_algorithm: StandardPairlistAlgorithm,
}

impl ParallelPairlistAlgorithm {
    pub fn new(use_chargegroups: bool) -> Self {
        Self {
            base_algorithm: StandardPairlistAlgorithm::new(use_chargegroups),
        }
    }

    /// Update pairlist in parallel
    pub fn update<BC: BoundaryCondition + Sync>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) where
        BC: Clone,
    {
        let n_atoms = topo.num_atoms();
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);

        // Parallel generation of pairs
        let pairs: Vec<(usize, usize, bool)> = (0..n_atoms)
            .into_par_iter()
            .flat_map(|i| {
                let mut local_pairs = Vec::new();
                let periodicity = periodicity.clone();

                for j in (i + 1)..n_atoms {
                    if topo.is_excluded_or_14(i, j) {
                        continue;
                    }

                    let r = periodicity.nearest_image(conf.current().pos[i], conf.current().pos[j]);
                    let dist2 = r.length_squared();

                    if dist2 < cutoff2_short {
                        local_pairs.push((i, j, true)); // true = short range
                    } else if dist2 < cutoff2_long {
                        local_pairs.push((i, j, false)); // false = long range
                    }
                }

                local_pairs
            })
            .collect();

        // Separate into short and long range
        pairlist.solute_short.clear();
        pairlist.solute_long.clear();

        for (i, j, is_short) in pairs {
            if is_short {
                pairlist.solute_short.push((i, j));
            } else {
                pairlist.solute_long.push((i, j));
            }
        }

        pairlist.reset_counter();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::configuration::BoxType;
    use crate::math::Rectangular;
    use crate::topology::ChargeGroup;

    // --- PairlistAlgorithm::from_imd dispatch tests (pure logic, no floats) ---

    fn is_standard(a: &PairlistAlgorithm) -> bool {
        matches!(a, PairlistAlgorithm::Standard(_))
    }

    fn is_cell_list(a: &PairlistAlgorithm) -> bool {
        matches!(a, PairlistAlgorithm::CellList(_))
    }

    #[test]
    fn test_from_imd_explicit_standard() {
        // algorithm=0 "standard" → Standard always
        assert!(is_standard(&PairlistAlgorithm::from_imd(0, 6000, BoxType::Rectangular, true)));
        assert!(is_standard(&PairlistAlgorithm::from_imd(0, 10, BoxType::Vacuum, false)));
    }

    #[test]
    fn test_from_imd_grid_falls_back_to_standard() {
        // algorithm=1 "grid" → ExtendedGrid not yet ported; must give Standard (safe fallback),
        // never CellList — preserving faithful gromosXX slot semantics.
        assert!(is_standard(&PairlistAlgorithm::from_imd(1, 10, BoxType::Vacuum, false)));
        assert!(is_standard(&PairlistAlgorithm::from_imd(1, 10, BoxType::Rectangular, true)));
        assert!(is_standard(&PairlistAlgorithm::from_imd(1, 6000, BoxType::Rectangular, true)));
    }

    #[test]
    fn test_from_imd_grid_cell_gives_cell_list() {
        // algorithm=2 "grid_cell" → CellList (Heinz & Hünenberger 2004), regardless of size/box
        assert!(is_cell_list(&PairlistAlgorithm::from_imd(2, 10, BoxType::Vacuum, false)));
        assert!(is_cell_list(&PairlistAlgorithm::from_imd(2, 10, BoxType::Rectangular, true)));
        assert!(is_cell_list(&PairlistAlgorithm::from_imd(2, 6000, BoxType::Rectangular, true)));
    }

    #[test]
    fn test_from_imd_standard_is_always_standard() {
        // algorithm=0 "standard" forces Standard regardless of size or box type —
        // ALGORITHM standard means Standard, no auto-heuristic overrides it.
        assert!(is_standard(&PairlistAlgorithm::from_imd(0, 648, BoxType::Rectangular, true)));
        assert!(is_standard(&PairlistAlgorithm::from_imd(0, 2998, BoxType::Rectangular, true)));
        assert!(is_standard(&PairlistAlgorithm::from_imd(0, 10000, BoxType::Rectangular, true)));
    }

    #[test]
    fn test_from_imd_unknown_algorithm_falls_back_to_standard() {
        // Unrecognised algorithm values → Standard (safe fallback, no auto-heuristic)
        assert!(is_standard(&PairlistAlgorithm::from_imd(99, 10000, BoxType::Rectangular, true)));
        assert!(is_standard(&PairlistAlgorithm::from_imd(-1, 648, BoxType::Rectangular, true)));
    }

    #[test]
    fn test_pairlist_container() {
        let pairlist = PairlistContainer::new(1.0, 1.4, 0.2);

        assert_eq!(pairlist.short_range_cutoff, 1.0);
        assert_eq!(pairlist.long_range_cutoff, 1.4);
        assert!(!pairlist.needs_update());
    }

    #[test]
    fn test_pairlist_update_frequency() {
        let mut pairlist = PairlistContainer::new(1.0, 1.4, 0.2);
        pairlist.update_frequency = 3;

        pairlist.step();
        assert!(!pairlist.needs_update());

        pairlist.step();
        assert!(!pairlist.needs_update());

        pairlist.step();
        assert!(pairlist.needs_update());

        pairlist.reset_counter();
        assert!(!pairlist.needs_update());
    }

    #[test]
    fn test_standard_pairlist() {
        let mut topo = Topology::new();
        
        // Set solute atoms FIRST - populate the atoms vector
        use crate::topology::Atom;
        for i in 0..4 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("A{}", i + 1),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 1.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        
        // Now resize arrays (uses num_atoms() which counts solute.atoms)
        topo.iac = vec![0; 4];
        topo.mass = vec![1.0; 4];
        topo.charge = vec![0.0; 4];
        topo.resize_atom_arrays();
        topo.compute_inverse_masses();
        
        // Initialize LJ parameter matrix (needed for pairlist)
        topo.lj_parameters = vec![vec![crate::topology::LJParameters { c6: 0.001, c12: 0.0001, cs6: 0.0005, cs12: 0.00005 }; 1]; 1];

        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.5, 0.0, 0.0); // Close
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0); // Far
        conf.current_mut().pos[3] = Vec3::new(0.0, 0.5, 0.0); // Close

        let mut pairlist = PairlistContainer::new(1.0, 2.5, 0.2);
        let algorithm = StandardPairlistAlgorithm::new(false);
        let periodicity = crate::math::Vacuum;

        algorithm.update(&topo, &conf, &mut pairlist, &periodicity);

        // Should have pairs (0,1), (0,3), and (1,3) in short range
        // and (0,2) in long range
        assert!(pairlist.solute_short.len() >= 2);
        assert!(pairlist.solute_long.len() >= 1);
    }

    /// Sort + min/max-normalize a pair list so two pairlists can be compared
    /// for set-equality regardless of traversal order.
    fn normalize_pairs(pairs: &Pairlist) -> Vec<(usize, usize)> {
        let mut v: Vec<(usize, usize)> = pairs.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
        v.sort_unstable();
        v
    }

    /// Assert that two pairlist containers contain the same pair *sets* in
    /// all four categories (FUTURE.md Dim 9a validation requirement).
    fn assert_pairlists_match(standard: &PairlistContainer, cell: &PairlistContainer) {
        assert_eq!(
            normalize_pairs(&standard.solute_short),
            normalize_pairs(&cell.solute_short),
            "solute_short mismatch"
        );
        assert_eq!(
            normalize_pairs(&standard.solute_long),
            normalize_pairs(&cell.solute_long),
            "solute_long mismatch"
        );
        assert_eq!(
            normalize_pairs(&standard.solvent_short),
            normalize_pairs(&cell.solvent_short),
            "solvent_short mismatch"
        );
        assert_eq!(
            normalize_pairs(&standard.solvent_long),
            normalize_pairs(&cell.solvent_long),
            "solvent_long mismatch"
        );
    }

    /// Build a pure-solvent topology of `n_molecules` single-atom "molecules"
    /// (each its own chargegroup), scattered across a cubic box of side
    /// `box_len` (including near the box edges, to exercise periodic wrap).
    fn build_solvent_topology(n_molecules: usize, box_len: f64) -> (Topology, Configuration) {
        use crate::topology::SolventAtomTemplate;

        let mut topo = Topology::new();
        topo.solvent_atom_template.push(SolventAtomTemplate {
            iac: 0,
            name: "OW".to_string(),
            mass: 16.0,
            charge: 0.0,
        });
        topo.solvate(n_molecules);

        let mut conf = Configuration::new(n_molecules, 1, 1);
        conf.current_mut().box_config = crate::configuration::Box::rectangular(box_len, box_len, box_len);
        for (i, pos) in conf.current_mut().pos.iter_mut().enumerate() {
            let f = i as f64;
            *pos = Vec3::new(
                (f * 1.37).rem_euclid(box_len),
                (f * 0.91 + 0.3).rem_euclid(box_len),
                (f * 2.03 + 0.6).rem_euclid(box_len),
            );
        }

        (topo, conf)
    }

    #[test]
    fn test_cell_list_matches_standard_solvent_only() {
        let (topo, conf) = build_solvent_topology(40, 6.0);
        let periodicity = Rectangular::new(Vec3::new(6.0, 6.0, 6.0));

        let mut standard_pl = PairlistContainer::new(1.0, 1.4, 0.0);
        StandardPairlistAlgorithm::new(true).update(&topo, &conf, &mut standard_pl, &periodicity);

        let mut cell_pl = PairlistContainer::new(1.0, 1.4, 0.0);
        CellListPairlistAlgorithm::new().update(&topo, &conf, &mut cell_pl, &periodicity);

        assert!(standard_pl.total_pairs() > 0);
        assert_pairlists_match(&standard_pl, &cell_pl);
    }

    #[test]
    fn test_cell_list_matches_standard_single_cell() {
        // box_len < long_range_cutoff + skin -> grid_dim == [1,1,1]: exercises
        // the single-cell / 27-offset-deduplication fallback path.
        let (topo, conf) = build_solvent_topology(15, 1.0);
        let periodicity = Rectangular::new(Vec3::new(1.0, 1.0, 1.0));

        let mut standard_pl = PairlistContainer::new(1.0, 1.4, 0.0);
        StandardPairlistAlgorithm::new(true).update(&topo, &conf, &mut standard_pl, &periodicity);

        let mut cell_pl = PairlistContainer::new(1.0, 1.4, 0.0);
        CellListPairlistAlgorithm::new().update(&topo, &conf, &mut cell_pl, &periodicity);

        assert!(standard_pl.total_pairs() > 0);
        assert_pairlists_match(&standard_pl, &cell_pl);
    }

    #[test]
    fn test_cell_list_matches_standard_with_solute() {
        use crate::topology::{Atom, Bond, SolventAtomTemplate};

        let mut topo = Topology::new();
        for i in 0..4 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("A{}", i + 1),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 1.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        // Bonds: (0,1) and (2,3) within chargegroups, (1,2) across them
        // (an excluded cross-CG pair, to exercise the exclusion check).
        topo.moltypes[0].bonds.push(Bond { i: 0, j: 1, bond_type: 0 });
        topo.moltypes[0].bonds.push(Bond { i: 2, j: 3, bond_type: 0 });
        topo.moltypes[0].bonds.push(Bond { i: 1, j: 2, bond_type: 0 });
        topo.iac = vec![0; 4];
        topo.mass = vec![1.0; 4];
        topo.charge = vec![0.0; 4];
        topo.resize_atom_arrays();
        topo.build_exclusions(false);

        topo.chargegroups.push(ChargeGroup { atoms: vec![0, 1] });
        topo.chargegroups.push(ChargeGroup { atoms: vec![2, 3] });
        topo.num_solute_chargegroups = 2;

        topo.solvent_atom_template.push(SolventAtomTemplate {
            iac: 0,
            name: "OW".to_string(),
            mass: 16.0,
            charge: 0.0,
        });
        topo.solvate(20);
        topo.compute_inverse_masses();

        let n_total = topo.num_atoms();
        let mut conf = Configuration::new(n_total, 1, 1);
        conf.current_mut().box_config = crate::configuration::Box::rectangular(3.0, 3.0, 3.0);

        // CG0 = {0,1} and CG1 = {2,3} are close enough for a solute-solute pair.
        conf.current_mut().pos[0] = Vec3::new(0.1, 0.1, 0.1);
        conf.current_mut().pos[1] = Vec3::new(0.3, 0.1, 0.1);
        conf.current_mut().pos[2] = Vec3::new(0.5, 0.5, 0.5);
        conf.current_mut().pos[3] = Vec3::new(0.7, 0.5, 0.5);

        for i in 4..n_total {
            let f = i as f64;
            conf.current_mut().pos[i] = Vec3::new(
                (f * 1.13).rem_euclid(3.0),
                (f * 0.77 + 0.2).rem_euclid(3.0),
                (f * 1.91 + 0.5).rem_euclid(3.0),
            );
        }

        // box_len / (long_range_cutoff + skin) = 3.0 / 1.4 -> grid_dim == 2,
        // the trickiest case for the neighbor-cell deduplication.
        let periodicity = Rectangular::new(Vec3::new(3.0, 3.0, 3.0));

        let mut standard_pl = PairlistContainer::new(1.0, 1.4, 0.0);
        StandardPairlistAlgorithm::new(true).update(&topo, &conf, &mut standard_pl, &periodicity);

        let mut cell_pl = PairlistContainer::new(1.0, 1.4, 0.0);
        CellListPairlistAlgorithm::new().update(&topo, &conf, &mut cell_pl, &periodicity);

        assert!(!standard_pl.solute_short.is_empty());
        assert!(standard_pl.total_pairs() > 0);
        assert_pairlists_match(&standard_pl, &cell_pl);
    }
}
