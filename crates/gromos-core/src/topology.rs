//! Topology module - molecular structure and force field parameters
//!
//! This is a direct Rust translation of GROMOS topology structures from:
//! - md++/src/topology/topology.h
//! - md++/src/topology/solute.h
//! - md++/src/topology/solvent.h
//!
//! # Dimension 10 — instancing model (Phase 1)
//!
//! The legacy `Solute` / `Vec<Solvent>` fields remain the authoritative source
//! for bonded interactions during the transition.  The new fields alongside them
//! implement the *representation* half of Dim 10:
//!
//! - `MoleculeType` stores per-atom parameters **once** per unique molecule type.
//! - `MoleculeInstance` places one copy at an `atom_offset` in the flat arrays
//!   and carries a `Role` attribute (Solute / Solvent).
//! - `Topology::moltypes` + `Topology::instances` are populated by
//!   `init_solute_moltype()` and `solvate()`.
//!
//! **Phase 1 invariants:**
//! - `instances[k]` corresponds 1-to-1 with `molecules[k]`.
//! - `role_of_atom(i)` → `Some(Role::Solvent)` iff atom `i` belongs to a
//!   solvent instance.  Returns `None` when the instance model is not yet
//!   initialised (pre-`solvate()` test topologies).
//! - `s:` in `AtomSelection` routes through `role_of_atom` when instances are
//!   available; falls back to `molecules[1..]` otherwise.
//!
//! Phase 2 (future): move bonds/angles/dihedrals into `MoleculeType`; make the
//! flat `iac`/`mass`/`charge` arrays into rebuildable caches; remove `Solute`
//! and `Vec<Solvent>`; add `promote(instance)` with CG/exclusion renumbering.

use crate::math::Vec3;

// ─── Dim 10: instancing model ────────────────────────────────────────────────

/// Semantic role of a molecule instance in the simulation.
///
/// Routes SHAKE/SETTLE dispatch, the solvent fast-path in the nonbonded
/// innerloop, and the `s:`/`m:` AtomSpecifier grammar.
/// In Phase 1 this is set at `solvate()` time; Phase 2 adds `promote()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Role {
    Solute,
    Solvent,
}

/// Per-atom parameters stored *once* inside a `MoleculeType`.
///
/// Equivalent to `Atom` but without the boolean flags (those belong to the
/// instance or to future attribute layers).
#[derive(Debug, Clone)]
pub struct MolTypeAtom {
    pub name: String,
    pub residue_nr: usize,
    pub residue_name: String,
    pub iac: usize,
    pub mass: f64,
    pub charge: f64,
}

/// Topology stored once for a unique molecule type.
///
/// Phase 1: only atom parameters.  Phase 2: bonds/angles/dihedrals move here.
#[derive(Debug, Clone)]
pub struct MoleculeType {
    /// Human-readable name ("SOLUTE", "SPC", …).
    pub name: String,
    /// One entry per atom in the template, in order.
    pub atoms: Vec<MolTypeAtom>,
}

impl MoleculeType {
    pub fn num_atoms(&self) -> usize { self.atoms.len() }
}

/// One copy of a `MoleculeType` placed in the system.
#[derive(Debug, Clone)]
pub struct MoleculeInstance {
    /// Index into `Topology::moltypes`.
    pub moltype_id: usize,
    /// First atom's 0-based global index in the flat position/force arrays.
    pub atom_offset: usize,
    /// Semantic role — drives dispatch and selection grammar.
    pub role: Role,
}

/// Atom properties
#[derive(Debug, Clone)]
pub struct Atom {
    pub name: String,
    pub residue_nr: usize,
    pub residue_name: String,
    pub iac: usize, // Integer atom code (atom type)
    pub mass: f64,
    pub charge: f64,
    pub is_perturbed: bool,
    pub is_polarisable: bool,
    pub is_coarse_grained: bool,
}

/// Two-body bonded term (bonds, harmonic constraints)
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub bond_type: usize,
}

/// Three-body bonded term (angles)
#[derive(Debug, Clone, Copy)]
pub struct Angle {
    pub i: usize,
    pub j: usize, // Central atom
    pub k: usize,
    pub angle_type: usize,
}

/// Four-body bonded term (proper and improper dihedrals)
#[derive(Debug, Clone, Copy)]
pub struct Dihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub dihedral_type: usize,
}

/// Cross-dihedral term (8 atoms)
#[derive(Debug, Clone, Copy)]
pub struct CrossDihedral {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub d: usize,
    pub e: usize,
    pub f: usize,
    pub g: usize,
    pub h: usize,
    pub cross_dihedral_type: usize,
}

/// Perturbed atom — state A and B LJ+charge parameters (PERTATOMPARAM block).
/// `lj_soft` and `crf_soft` are the per-atom soft-core alphas from the file;
/// they are scaled by ALPHLJ/ALPHC from the imd PERTURBATION block at run time.
#[derive(Debug, Clone)]
pub struct PerturbedAtom {
    pub seq: usize,      // 0-indexed atom number in topology
    pub a_iac: usize,    // State A IAC (0-indexed)
    pub a_mass: f64,
    pub a_charge: f64,
    pub b_iac: usize,    // State B IAC (0-indexed)
    pub b_mass: f64,
    pub b_charge: f64,
    pub lj_soft: f64,    // per-atom alpha_LJ soft-core
    pub crf_soft: f64,   // per-atom alpha_CRF soft-core
}

/// Perturbed atom pair — special LJ interaction between two perturbed atoms
/// (PERTATOMPAIR block). These pairs are excluded from the regular pairlist
/// and handled with their own A/B LJ types.
///
/// **INVARIANT: all index fields are 0-indexed.** GROMOS files store `A_type`
/// and `B_type` as 1-indexed integers; `gromos-io/ptp.rs` subtracts 1 at the
/// parse boundary.  Never store a raw file value here.
///
/// `B_type` is `Option` because the file value `0` is a sentinel meaning "the
/// atom pair disappears in state B" (no interaction).  The parse boundary maps
/// file 0 → `None` and file N → `Some(N-1)`.
///
/// Interaction type encoding (same for A and B):
///   `0` (file 1) → normal full LJ (c6/c12)
///   `1` (file 2) → 1-4 LJ (cs6/cs12)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedAtomPair {
    pub i: usize,               // 0-indexed, i < j
    pub j: usize,
    pub a_type: usize,          // State A interaction type (0-indexed)
    pub b_type: Option<usize>,  // State B interaction type (0-indexed); None = absent in state B
}

/// Perturbed bond (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedBond {
    pub i: usize,
    pub j: usize,
    pub a_type: usize, // State A bond type
    pub b_type: usize, // State B bond type
}

/// Perturbed angle (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedAngle {
    pub i: usize,
    pub j: usize, // Central atom
    pub k: usize,
    pub a_type: usize, // State A angle type
    pub b_type: usize, // State B angle type
}

/// Perturbed dihedral (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedDihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub a_type: usize, // State A dihedral type
    pub b_type: usize, // State B dihedral type
}

/// Bond force field parameters (GROMOS format)
#[derive(Debug, Clone, Copy)]
pub struct BondParameters {
    pub k_quartic: f64,  // Quartic force constant
    pub k_harmonic: f64, // Harmonic force constant
    pub r0: f64,         // Equilibrium bond length
}

/// Angle force field parameters (GROMOS format)
#[derive(Debug, Clone, Copy)]
pub struct AngleParameters {
    pub k_cosine: f64,   // Force constant (harmonic in cosine)
    pub k_harmonic: f64, // Force constant (harmonic in angle)
    pub theta0: f64,     // Equilibrium angle in radians
}

/// Dihedral force field parameters
#[derive(Debug, Clone, Copy)]
pub struct DihedralParameters {
    pub k: f64,     // Force constant
    pub pd: f64,    // Phase shift
    pub cospd: f64, // Cos(phase)
    pub m: i32,     // Multiplicity
}

/// Improper dihedral force field parameters
#[derive(Debug, Clone, Copy)]
pub struct ImproperDihedralParameters {
    pub k: f64,  // Force constant
    pub q0: f64, // Equilibrium improper dihedral angle
}

/// Lennard-Jones parameters
#[derive(Debug, Clone, Copy)]
pub struct LJParameters {
    pub c6: f64,   // C6 coefficient (attractive: -C6/r^6)
    pub c12: f64,  // C12 coefficient (repulsive: C12/r^12)
    pub cs6: f64,  // Softcore C6 (for free energy calculations)
    pub cs12: f64, // Softcore C12
}

impl LJParameters {
    /// Create standard LJ parameters (without softcore)
    pub fn new(c6: f64, c12: f64) -> Self {
        Self {
            c6,
            c12,
            cs6: c6,
            cs12: c12,
        }
    }

    /// Create LJ parameters with explicit 1-4 (cs6, cs12) values
    pub fn new_with_14(c6: f64, c12: f64, cs6: f64, cs12: f64) -> Self {
        Self { c6, c12, cs6, cs12 }
    }

    /// Calculate sigma (size parameter) from C6 and C12
    /// sigma = (C12/C6)^(1/6)
    pub fn sigma(&self) -> f64 {
        if self.c6 > 0.0 {
            (self.c12 / self.c6).powf(1.0 / 6.0)
        } else {
            0.0
        }
    }

    /// Calculate epsilon (well depth) from C6 and C12
    /// epsilon = C6^2 / (4 * C12)
    pub fn epsilon(&self) -> f64 {
        if self.c12 > 0.0 {
            self.c6 * self.c6 / (4.0 * self.c12)
        } else {
            0.0
        }
    }
}

impl Default for LJParameters {
    fn default() -> Self {
        Self {
            c6: 0.0,
            c12: 0.0,
            cs6: 0.0,
            cs12: 0.0,
        }
    }
}

/// Exclusion list for an atom (atoms that don't interact via nonbonded forces)
///
/// Uses a sorted Vec for fast binary-search lookups. Matches the GROMOS C++
/// convention (sorted vector with binary search). For typical exclusion counts
/// (2-6 per atom), this is much faster than HashSet due to cache locality.
pub type Exclusions = Vec<usize>;

/// Chargegroup - atoms whose charge sum is zero (or small)
///
/// Chargegroups are used for cutoff optimizations:
/// - Distance calculated between chargegroup centers
/// - All atoms in both groups included if distance < cutoff
#[derive(Debug, Clone)]
pub struct ChargeGroup {
    pub atoms: Vec<usize>,
}

impl ChargeGroup {
    /// Calculate center of geometry for this chargegroup
    pub fn center_of_geometry(&self, positions: &[Vec3]) -> Vec3 {
        if self.atoms.is_empty() {
            return Vec3::ZERO;
        }

        let sum: Vec3 = self.atoms.iter().map(|&i| positions[i]).sum();

        sum / self.atoms.len() as f64
    }
}

/// Solute molecule structure
#[derive(Debug, Clone)]
pub struct Solute {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub proper_dihedrals: Vec<Dihedral>,
    pub improper_dihedrals: Vec<Dihedral>,
    pub cross_dihedrals: Vec<CrossDihedral>,
}

impl Solute {
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            proper_dihedrals: Vec::new(),
            improper_dihedrals: Vec::new(),
            cross_dihedrals: Vec::new(),
        }
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
}

/// Soft-core perturbed harmonic bond (PERTBONDSOFT).
/// Used when a bond is absent in one state (K=0), requiring soft-core treatment.
/// `b_type = None` means the bond is absent in state B (K_B = 0, r0_B = r0_A).
#[derive(Debug, Clone)]
pub struct SoftBond {
    pub i: usize,
    pub j: usize,
    pub a_type: usize,          // index into bond_parameters (uses k_harmonic)
    pub b_type: Option<usize>,  // None → K_B=0, r0_B=r0_A
    pub alpha: f64,
}

/// Soft-core perturbed cos-harmonic angle (PERTANGLESOFT).
/// `b_type = None` means the angle is absent in state B (K_B = 0, cos0_B = cos0_A).
#[derive(Debug, Clone)]
pub struct SoftAngle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub a_type: usize,
    pub b_type: Option<usize>,
    pub alpha: f64,
}

/// Soft-core perturbed improper dihedral (PERTIMPROPERDIHSOFT).
/// `b_type = None` means the improper is absent in state B (K_B = 0, q0_B = q0_A).
#[derive(Debug, Clone)]
pub struct SoftImproper {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub a_type: usize,
    pub b_type: Option<usize>,
    pub alpha: f64,
}

/// Perturbed solute for FEP calculations
///
/// Contains dual-topology (A/B state) bonded terms for free energy perturbation
#[derive(Debug, Clone)]
pub struct PerturbedSolute {
    /// Perturbed atoms: state A/B IAC, mass, charge, soft-core alphas.
    /// Indexed sparsely — only atoms listed in PERTATOMPARAM.
    pub atoms: Vec<PerturbedAtom>,
    /// Perturbed atom pairs (PERTATOMPAIR): special LJ between two perturbed atoms.
    pub atom_pairs: Vec<PerturbedAtomPair>,
    pub bonds: Vec<PerturbedBond>,
    pub angles: Vec<PerturbedAngle>,
    pub proper_dihedrals: Vec<PerturbedDihedral>,
    pub improper_dihedrals: Vec<PerturbedDihedral>,
    /// Soft-core perturbed bonds (PERTBONDSOFT): bond absent in one state.
    pub soft_bonds: Vec<SoftBond>,
    /// Soft-core perturbed angles (PERTANGLESOFT): angle absent in one state.
    pub soft_angles: Vec<SoftAngle>,
    /// Soft-core perturbed impropers (PERTIMPROPERDIHSOFT): improper absent in one state.
    pub soft_impropers: Vec<SoftImproper>,
}

impl PerturbedSolute {
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            atom_pairs: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            proper_dihedrals: Vec::new(),
            improper_dihedrals: Vec::new(),
            soft_bonds: Vec::new(),
            soft_angles: Vec::new(),
            soft_impropers: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty() && self.bonds.is_empty()
            && self.angles.is_empty() && self.proper_dihedrals.is_empty()
            && self.improper_dihedrals.is_empty()
    }
}

impl Default for PerturbedSolute {
    fn default() -> Self {
        Self::new()
    }
}

/// Solvent molecule structure (typically water)
#[derive(Debug, Clone)]
pub struct Solvent {
    pub name: String,
    pub atoms: Vec<Atom>,
    pub num_molecules: usize,
}

impl Solvent {
    pub fn new(name: String) -> Self {
        Self {
            name,
            atoms: Vec::new(),
            num_molecules: 0,
        }
    }

    pub fn atoms_per_molecule(&self) -> usize {
        self.atoms.len()
    }

    pub fn total_atoms(&self) -> usize {
        self.atoms.len() * self.num_molecules
    }
}

/// Solvent atom template (one molecule's worth of atom properties)
#[derive(Debug, Clone)]
pub struct SolventAtomTemplate {
    pub iac: usize,
    pub name: String,
    pub mass: f64,
    pub charge: f64,
}

/// Solvent constraint template (within one molecule)
#[derive(Debug, Clone)]
pub struct SolventConstraintTemplate {
    pub i: usize, // 0-indexed within solvent molecule
    pub j: usize,
    pub length: f64,
}

/// Distance restraint specification (virtual atom type 0 only)
#[derive(Debug, Clone)]
pub struct DistanceRestraintSpec {
    pub atom1: usize,
    pub atom2: usize,
    pub r0: f64,
    pub w0: f64,
    pub rah: i32,
}

/// Perturbed distance restraint specification (virtual atom type 0 only)
#[derive(Debug, Clone)]
pub struct PerturbedDistanceRestraintSpec {
    pub atom1: usize,
    pub atom2: usize,
    pub n: i32,
    pub m: i32,
    pub a_r0: f64,
    pub b_r0: f64,
    pub a_w0: f64,
    pub b_w0: f64,
    pub rah: i32,
}

/// Main topology structure containing all molecular information
#[derive(Debug, Clone)]
pub struct Topology {
    // Solute and solvent
    pub solute: Solute,
    pub perturbed_solute: PerturbedSolute, // FEP dual-topology terms
    pub solvents: Vec<Solvent>,

    // Per-atom properties (flat arrays for all atoms)

    /// Integer atom codes (atom types).
    ///
    /// **INVARIANT: 0-indexed** (0..N\_types-1), matching LJ matrix rows/cols.
    /// GROMOS files store IAC as 1-indexed; `gromos-io/topology.rs` subtracts 1
    /// at the parse boundary.  `ptp.rs` does the same for `PerturbedAtom.a_iac`
    /// and `b_iac`.  Never store a raw 1-indexed file value here.
    pub iac: Vec<usize>,
    pub mass: Vec<f64>,         // Atomic masses
    pub inverse_mass: Vec<f64>, // Precomputed 1/mass for efficiency
    pub charge: Vec<f64>,       // Atomic charges

    // Exclusions (atoms that don't interact via nonbonded)
    pub exclusions: Vec<Exclusions>, // exclusions[i] = set of atoms excluded from i
    pub one_four_pairs: Vec<Vec<usize>>, // 1-4 pairs (special scaling)

    // Chargegroups
    pub chargegroups: Vec<ChargeGroup>,
    pub atom_to_chargegroup: Vec<usize>, // atom -> chargegroup index
    pub num_solute_chargegroups: usize,  // boundary: CGs [0..n) are solute, [n..) are solvent

    // Temperature and pressure coupling groups
    pub temperature_groups: Vec<Vec<usize>>, // Atoms in each T-coupling group
    pub pressure_groups: Vec<std::ops::Range<usize>>, // Pressure groups as atom ranges (gromosXX: PRESSUREGROUPS)
    pub atom_to_temperature_group: Vec<usize>,
    pub atom_to_pressure_group: Vec<usize>,

    // Energy groups (for energy monitoring)
    pub energy_groups: Vec<Vec<usize>>,
    pub atom_to_energy_group: Vec<usize>,

    // Molecule boundaries (for molecule-based operations)
    pub molecules: Vec<std::ops::Range<usize>>, // molecule_start..molecule_end

    // ── Dim 10 instancing model ──────────────────────────────────────────────
    // Populated by init_solute_moltype() + solvate().
    // instances[k] corresponds 1-to-1 with molecules[k].
    // Empty until those methods are called (e.g. in raw test topologies).
    /// Registry: one entry per unique molecule type (parameters stored once).
    pub moltypes: Vec<MoleculeType>,
    /// All instances in system order; instances[k].atom_offset == molecules[k].start.
    pub instances: Vec<MoleculeInstance>,

    // Force field parameters
    pub bond_parameters: Vec<BondParameters>,
    pub angle_parameters: Vec<AngleParameters>,
    pub dihedral_parameters: Vec<DihedralParameters>,
    pub improper_dihedral_parameters: Vec<ImproperDihedralParameters>,
    pub lj_parameters: Vec<Vec<LJParameters>>, // [type_i][type_j] matrix

    // Solvent template (read from topology, expanded by solvate())
    pub solvent_atom_template: Vec<SolventAtomTemplate>,
    pub solvent_constraint_template: Vec<SolventConstraintTemplate>,

    // Chargegroup codes (CGC) from topology for solute atoms
    pub chargegroup_codes: Vec<usize>,

    // Distance restraints (loaded from .distres file)
    pub distance_restraints: Vec<DistanceRestraintSpec>,
    pub perturbed_distance_restraints: Vec<PerturbedDistanceRestraintSpec>,

    // FEP/TI: per-solute-atom flag, true if listed in PERTATOMPARAM
    pub is_perturbed: Vec<bool>,
}

impl Topology {
    pub fn new() -> Self {
        Self {
            solute: Solute::new(),
            perturbed_solute: PerturbedSolute::new(),
            solvents: Vec::new(),
            iac: Vec::new(),
            mass: Vec::new(),
            inverse_mass: Vec::new(),
            charge: Vec::new(),
            exclusions: Vec::new(),
            one_four_pairs: Vec::new(),
            chargegroups: Vec::new(),
            atom_to_chargegroup: Vec::new(),
            num_solute_chargegroups: 0,
            temperature_groups: Vec::new(),
            pressure_groups: Vec::new(),  // populated from PRESSUREGROUPS topology block
            atom_to_temperature_group: Vec::new(),
            atom_to_pressure_group: Vec::new(),
            energy_groups: Vec::new(),
            atom_to_energy_group: Vec::new(),
            molecules: Vec::new(),
            moltypes: Vec::new(),
            instances: Vec::new(),
            bond_parameters: Vec::new(),
            angle_parameters: Vec::new(),
            dihedral_parameters: Vec::new(),
            improper_dihedral_parameters: Vec::new(),
            lj_parameters: Vec::new(),
            solvent_atom_template: Vec::new(),
            solvent_constraint_template: Vec::new(),
            chargegroup_codes: Vec::new(),
            distance_restraints: Vec::new(),
            perturbed_distance_restraints: Vec::new(),
            is_perturbed: Vec::new(),
        }
    }

    /// Total number of atoms in the system
    pub fn num_atoms(&self) -> usize {
        self.solute.num_atoms() + self.solvents.iter().map(|s| s.total_atoms()).sum::<usize>()
    }

    /// Number of solute atoms
    pub fn num_solute_atoms(&self) -> usize {
        self.solute.num_atoms()
    }

    /// Number of atom types
    pub fn num_atom_types(&self) -> usize {
        self.lj_parameters.len()
    }

    /// Get LJ parameters for atom pair (i, j)
    #[inline]
    pub fn lj_parameter(&self, type_i: usize, type_j: usize) -> &LJParameters {
        &self.lj_parameters[type_i][type_j]
    }

    /// Check if atoms i and j are excluded from nonbonded interactions
    /// (1-2 and 1-3 bonded pairs only — used for RF excluded interactions)
    #[inline]
    pub fn is_excluded(&self, i: usize, j: usize) -> bool {
        self.exclusions[i].binary_search(&j).is_ok() || self.exclusions[j].binary_search(&i).is_ok()
    }

    /// Check if atoms i and j are excluded OR are a 1-4 pair.
    /// Used for pairlist exclusion (gromosXX: all_exclusion).
    #[inline]
    pub fn is_excluded_or_14(&self, i: usize, j: usize) -> bool {
        if self.exclusions[i].binary_search(&j).is_ok() || self.exclusions[j].binary_search(&i).is_ok() {
            return true;
        }
        if !self.one_four_pairs.is_empty() {
            if self.one_four_pairs[i].contains(&j) || self.one_four_pairs[j].contains(&i) {
                return true;
            }
        }
        false
    }

    /// Get chargegroup index for atom i
    #[inline]
    pub fn chargegroup(&self, i: usize) -> usize {
        self.atom_to_chargegroup[i]
    }

    /// Get temperature group index for atom i
    #[inline]
    pub fn temperature_group(&self, i: usize) -> usize {
        self.atom_to_temperature_group[i]
    }

    /// Get energy group index for atom i
    #[inline]
    pub fn energy_group(&self, i: usize) -> usize {
        self.atom_to_energy_group[i]
    }

    // ── Dim 10 accessors ─────────────────────────────────────────────────────

    /// Internal: look up the MolTypeAtom for global atom index `i` via instances.
    fn moltype_atom(&self, i: usize) -> Option<&MolTypeAtom> {
        let mol = self.molecule_nr(i)?;
        let inst = self.instances.get(mol)?;
        let mt = self.moltypes.get(inst.moltype_id)?;
        mt.atoms.get(i - inst.atom_offset)
    }

    /// Role of the molecule instance that owns atom `i`.
    /// Returns `None` when the instance model hasn't been initialised yet.
    pub fn role_of_atom(&self, i: usize) -> Option<Role> {
        self.molecule_nr(i)
            .and_then(|m| self.instances.get(m))
            .map(|inst| inst.role)
    }

    /// `true` iff atom `i` belongs to a `Role::Solvent` instance.
    pub fn is_solvent_atom(&self, i: usize) -> bool {
        self.role_of_atom(i) == Some(Role::Solvent)
    }

    /// Change the role of molecule instance `mol_idx` to `Role::Solute`.
    ///
    /// This is the "promote one solvent molecule to solute" operation from Dim 10.
    /// Phase 2 will add CG/exclusion renumbering; for now it's a label flip.
    pub fn promote(&mut self, mol_idx: usize) -> Result<(), String> {
        match self.instances.get_mut(mol_idx) {
            Some(inst) => { inst.role = Role::Solute; Ok(()) }
            None => Err(format!("molecule index {mol_idx} out of range (have {})", self.instances.len())),
        }
    }

    /// Build the solute `MoleculeType` + `MoleculeInstance` from the existing
    /// `solute.atoms` and flat parameter arrays.
    ///
    /// Safe to call repeatedly (no-op if already initialised).
    /// Called by `solvate()` and by `gromos-io::build_topology()`.
    pub fn init_solute_moltype(&mut self) {
        if !self.instances.is_empty() { return; }
        let n_sol = self.solute.num_atoms();
        if n_sol == 0 { return; }

        let mt_atoms: Vec<MolTypeAtom> = self.solute.atoms.iter().enumerate()
            .map(|(i, a)| MolTypeAtom {
                name:         a.name.clone(),
                residue_nr:   a.residue_nr,
                residue_name: a.residue_name.clone(),
                iac:   self.iac.get(i).copied().unwrap_or(a.iac),
                mass:  self.mass.get(i).copied().unwrap_or(a.mass),
                charge: self.charge.get(i).copied().unwrap_or(a.charge),
            })
            .collect();

        self.moltypes.push(MoleculeType { name: "SOLUTE".to_string(), atoms: mt_atoms });
        self.instances.push(MoleculeInstance { moltype_id: 0, atom_offset: 0, role: Role::Solute });

        // Ensure molecules[0] is set (solvate() may not have been called yet)
        if self.molecules.is_empty() && n_sol > 0 {
            self.molecules.push(0..n_sol);
        }
    }

    /// Atom name for atom `i` (covers solute + all expanded solvent). Returns `None` if out of range.
    ///
    /// Prefers the Dim 10 moltype path when instances are populated; falls back to
    /// the legacy `solute.atoms` + `solvent_atom_template` approach.
    pub fn atom_name(&self, i: usize) -> Option<&str> {
        if !self.instances.is_empty() {
            return self.moltype_atom(i).map(|a| a.name.as_str());
        }
        // Legacy fallback (pre-init_solute_moltype)
        let n_sol = self.solute.num_atoms();
        if i < n_sol {
            return self.solute.atoms.get(i).map(|a| a.name.as_str());
        }
        let aps = self.solvent_atom_template.len();
        if aps > 0 {
            return self.solvent_atom_template.get((i - n_sol) % aps).map(|a| a.name.as_str());
        }
        None
    }

    /// Residue number (1-based) for atom `i`.
    pub fn residue_nr(&self, i: usize) -> Option<usize> {
        if !self.instances.is_empty() {
            return self.moltype_atom(i).map(|a| a.residue_nr);
        }
        let n_sol = self.solute.num_atoms();
        if i < n_sol {
            return self.solute.atoms.get(i).map(|a| a.residue_nr);
        }
        let aps = self.solvent_atom_template.len();
        if aps > 0 { return Some((i - n_sol) / aps + 1); }
        None
    }

    /// Residue name for atom `i` (covers solute + solvent).
    pub fn residue_name(&self, i: usize) -> Option<&str> {
        if !self.instances.is_empty() {
            return self.moltype_atom(i).map(|a| a.residue_name.as_str());
        }
        let n_sol = self.solute.num_atoms();
        if i < n_sol {
            return self.solute.atoms.get(i).map(|a| a.residue_name.as_str());
        }
        if !self.solvents.is_empty() { return Some(self.solvents[0].name.as_str()); }
        None
    }

    /// 0-based molecule index for atom `i` (uses the `molecules` ranges populated by `solvate()`).
    /// Returns `None` if `molecules` is empty or atom is out of all molecule ranges.
    pub fn molecule_nr(&self, i: usize) -> Option<usize> {
        // molecules is sorted by start; binary-search on start to find candidate
        let idx = self.molecules.partition_point(|r| r.start <= i);
        if idx > 0 {
            let mol = idx - 1;
            if self.molecules[mol].contains(&i) {
                return Some(mol);
            }
        }
        None
    }

    /// Build exclusions from bonded topology
    ///
    /// Excludes:
    /// - 1-2 bonded neighbors
    /// - 1-3 angle neighbors
    /// - Optionally 1-4 dihedral neighbors (or use special 1-4 scaling)
    pub fn build_exclusions(&mut self, exclude_14: bool) {
        let n_atoms = self.num_atoms();
        self.exclusions = vec![Vec::new(); n_atoms];

        // Exclude bonded neighbors (1-2)
        for bond in &self.solute.bonds {
            self.exclusions[bond.i].push(bond.j);
            self.exclusions[bond.j].push(bond.i);
        }

        // Exclude angle neighbors (1-3)
        for angle in &self.solute.angles {
            self.exclusions[angle.i].push(angle.k);
            self.exclusions[angle.k].push(angle.i);
        }

        // Optionally exclude or mark 1-4 pairs
        if exclude_14 {
            for dihedral in &self.solute.proper_dihedrals {
                self.exclusions[dihedral.i].push(dihedral.l);
                self.exclusions[dihedral.l].push(dihedral.i);
            }
        } else {
            // Store as special 1-4 pairs
            self.one_four_pairs = vec![Vec::new(); n_atoms];
            for dihedral in &self.solute.proper_dihedrals {
                self.one_four_pairs[dihedral.i].push(dihedral.l);
                self.one_four_pairs[dihedral.l].push(dihedral.i);
            }
        }

        // Sort and deduplicate for fast binary_search lookups
        for excl in &mut self.exclusions {
            excl.sort_unstable();
            excl.dedup();
        }
    }

    /// Initialize LJ parameter matrix from combining rules
    ///
    /// GROMOS uses geometric mean combining rules:
    /// - C6_ij = sqrt(C6_ii * C6_jj)
    /// - C12_ij = sqrt(C12_ii * C12_jj)
    pub fn build_lj_matrix(&mut self, lj_types: &[LJParameters]) {
        let n_types = lj_types.len();
        self.lj_parameters = vec![vec![LJParameters::new(0.0, 0.0); n_types]; n_types];

        for i in 0..n_types {
            for j in 0..n_types {
                let c6 = (lj_types[i].c6 * lj_types[j].c6).sqrt();
                let c12 = (lj_types[i].c12 * lj_types[j].c12).sqrt();
                self.lj_parameters[i][j] = LJParameters::new(c6, c12);
            }
        }
    }

    /// Resize all per-atom arrays to match number of atoms
    pub fn resize_atom_arrays(&mut self) {
        let n_atoms = self.num_atoms();
        self.iac.resize(n_atoms, 0);
        self.mass.resize(n_atoms, 0.0);
        self.inverse_mass.resize(n_atoms, 0.0);
        self.charge.resize(n_atoms, 0.0);
        self.exclusions.resize(n_atoms, Vec::new());
        self.one_four_pairs.resize(n_atoms, Vec::new());
        self.atom_to_chargegroup.resize(n_atoms, 0);
        self.atom_to_temperature_group.resize(n_atoms, 0);
        self.atom_to_pressure_group.resize(n_atoms, 0);
        self.atom_to_energy_group.resize(n_atoms, 0);
    }

    /// Compute inverse masses for all atoms (for efficient integration)
    pub fn compute_inverse_masses(&mut self) {
        self.inverse_mass = self
            .mass
            .iter()
            .map(|&m| if m > 0.0 { 1.0 / m } else { 0.0 })
            .collect();
    }

    /// Expand solvent molecules (gromosXX: topo.solvate(0, nsm))
    ///
    /// Uses the solvent template stored on the topology to create `nsm` copies.
    /// Expands flat arrays (mass, charge, iac), chargegroups, and exclusions.
    /// Must be called after build_topology() and before reading coordinates.
    pub fn solvate(&mut self, nsm: usize) {
        if nsm == 0 || self.solvent_atom_template.is_empty() {
            return;
        }

        let atoms_per_solvent = self.solvent_atom_template.len();
        let n_solute = self.solute.num_atoms();

        // Dim 10: ensure solute moltype + instance are initialised before we add solvent.
        self.init_solute_moltype();

        // Add solute as one molecule (if it has atoms)
        if n_solute > 0 && self.molecules.is_empty() {
            self.molecules.push(0..n_solute);
        }

        let n_solvent_atoms = atoms_per_solvent * nsm;
        let n_total = n_solute + n_solvent_atoms;

        // Create Solvent entry
        let mut solvent = Solvent::new("SOLV".to_string());
        for sa in &self.solvent_atom_template {
            solvent.atoms.push(Atom {
                name: sa.name.clone(),
                residue_nr: 0,
                residue_name: "SOLV".to_string(),
                iac: sa.iac,
                mass: sa.mass,
                charge: sa.charge,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        solvent.num_molecules = nsm;
        self.solvents.push(solvent);

        // Expand flat arrays for all solvent copies
        for _mol in 0..nsm {
            for sa in &self.solvent_atom_template {
                self.mass.push(sa.mass);
                self.charge.push(sa.charge);
                self.iac.push(sa.iac);
            }
        }

        self.compute_inverse_masses();

        // Each solvent molecule is its own chargegroup
        for mol in 0..nsm {
            let base = n_solute + mol * atoms_per_solvent;
            let atoms: Vec<usize> = (base..base + atoms_per_solvent).collect();
            self.chargegroups.push(ChargeGroup { atoms });
        }

        // Rebuild atom_to_chargegroup mapping
        self.atom_to_chargegroup = vec![0; n_total];
        for (cg_idx, cg) in self.chargegroups.iter().enumerate() {
            for &atom in &cg.atoms {
                if atom < n_total {
                    self.atom_to_chargegroup[atom] = cg_idx;
                }
            }
        }

        // Extend exclusions for solvent atoms
        self.exclusions.resize(n_total, Vec::new());

        // Solvent intra-molecular exclusions: all atoms within each molecule exclude each other
        for mol in 0..nsm {
            let base = n_solute + mol * atoms_per_solvent;
            // Populate molecules list
            self.molecules.push(base..base + atoms_per_solvent);
            // Each solvent molecule is its own pressure group (gromosXX convention)
            self.pressure_groups.push(base..base + atoms_per_solvent);
            for a in 0..atoms_per_solvent {
                for b in (a + 1)..atoms_per_solvent {
                    self.exclusions[base + a].push(base + b);
                    self.exclusions[base + b].push(base + a);
                }
            }
        }

        // Sort solvent exclusions for binary_search
        for i in n_solute..n_total {
            self.exclusions[i].sort_unstable();
        }

        // Rebuild LJ matrix if solvent introduces new IAC types
        let max_iac = self.iac.iter().max().copied().unwrap_or(0);
        let n_types = max_iac + 1;
        if n_types > self.lj_parameters.len() {
            let old_len = self.lj_parameters.len();
            self.lj_parameters.resize(n_types, vec![LJParameters::default(); n_types]);
            for row in self.lj_parameters.iter_mut() {
                row.resize(n_types, LJParameters::default());
            }
            log::debug!("LJ matrix expanded from {}x{} to {}x{} for solvent types",
                old_len, old_len, n_types, n_types);
        }

        // Dim 10: add solvent MoleculeType (stored once) + one Instance per molecule.
        {
            let sv_name = self.solvents.last()
                .map(|sv| sv.name.clone())
                .unwrap_or_else(|| "SOLV".to_string());

            let sv_moltype = MoleculeType {
                name: sv_name,
                atoms: self.solvent_atom_template.iter().enumerate().map(|(idx, sa)| {
                    MolTypeAtom {
                        name:         sa.name.clone(),
                        residue_nr:   idx + 1,
                        residue_name: "SOLV".to_string(),
                        iac:   sa.iac,
                        mass:  sa.mass,
                        charge: sa.charge,
                    }
                }).collect(),
            };
            let solvent_moltype_id = self.moltypes.len();
            self.moltypes.push(sv_moltype);

            for mol in 0..nsm {
                let atom_offset = n_solute + mol * atoms_per_solvent;
                self.instances.push(MoleculeInstance {
                    moltype_id: solvent_moltype_id,
                    atom_offset,
                    role: Role::Solvent,
                });
            }
        }

        log::debug!("Solvated: {} solute + {} solvent ({} molecules × {} atoms) = {} total, {} chargegroups",
            n_solute, n_solvent_atoms, nsm, atoms_per_solvent, n_total,
            self.chargegroups.len());
    }
}

impl Default for Topology {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lj_parameters() {
        // For sigma = 0.34 nm and epsilon = 0.65 kJ/mol (typical for methane)
        // c6 = 4 * epsilon * sigma^6 = 4 * 0.65 * 0.34^6 = 0.00398
        // c12 = 4 * epsilon * sigma^12 = 4 * 0.65 * 0.34^12 = 4.6e-7
        // Or simpler: sigma = (c12/c6)^(1/6), epsilon = c6^2/(4*c12)
        
        // Use values that give known sigma/epsilon
        let sigma_expected: f64 = 0.34;
        let epsilon_expected: f64 = 0.65;
        let c6 = 4.0 * epsilon_expected * sigma_expected.powi(6);
        let c12 = 4.0 * epsilon_expected * sigma_expected.powi(12);
        
        let lj = LJParameters::new(c6, c12);

        // Check sigma
        let sigma = lj.sigma();
        assert!((sigma - sigma_expected).abs() < 0.01, "Sigma mismatch: got {}", sigma);

        // Check epsilon
        let epsilon = lj.epsilon();
        assert!((epsilon - epsilon_expected).abs() < 0.01, "Epsilon mismatch: got {}", epsilon);
    }

    #[test]
    fn test_topology_creation() {
        let mut topo = Topology::new();

        // Add a simple water molecule
        let mut water = Solvent::new("SOL".to_string());
        water.atoms.push(Atom {
            name: "OW".to_string(),
            residue_nr: 1,
            residue_name: "SOL".to_string(),
            iac: 0,
            mass: 15.9994,
            charge: -0.82,
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: false,
        });
        water.num_molecules = 100;

        topo.solvents.push(water);

        assert_eq!(topo.num_atoms(), 100);
    }

    #[test]
    fn test_exclusions() {
        let mut topo = Topology::new();
        
        // Add 3 atoms first so num_atoms() returns 3
        for i in 0..3 {
            topo.solute.atoms.push(Atom {
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
        
        topo.solute.bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
        });
        topo.solute.bonds.push(Bond {
            i: 1,
            j: 2,
            bond_type: 0,
        });
        topo.solute.angles.push(Angle {
            i: 0,
            j: 1,
            k: 2,
            angle_type: 0,
        });

        topo.resize_atom_arrays();
        topo.build_exclusions(false);

        // 0-1 should be excluded (bonded)
        assert!(topo.is_excluded(0, 1));

        // 0-2 should be excluded (1-3 neighbors)
        assert!(topo.is_excluded(0, 2));
    }

    #[test]
    fn test_lj_matrix() {
        let mut topo = Topology::new();

        // Two atom types
        let lj_types = vec![
            LJParameters::new(0.001, 0.0001), // Type 0
            LJParameters::new(0.002, 0.0002), // Type 1
        ];

        topo.build_lj_matrix(&lj_types);

        // Check diagonal elements
        assert!((topo.lj_parameter(0, 0).c6 - 0.001).abs() < 1e-10);
        assert!((topo.lj_parameter(1, 1).c6 - 0.002).abs() < 1e-10);

        // Check off-diagonal (geometric mean)
        let c6_01 = (0.001_f64 * 0.002_f64).sqrt();
        assert!((topo.lj_parameter(0, 1).c6 - c6_01).abs() < 1e-10);
    }
}
