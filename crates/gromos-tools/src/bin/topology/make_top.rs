//! make_top - Build GROMOS topology from building blocks
//!
//! Usage: make_top @build <mtb> @param <ifp> @seq <residues...> [@solv <solvent>]
//!    or: make_top @f <argfile>
//!    or: make_top --build <mtb> --param <ifp> --seq <residues...> [--solv <solvent>]
//!
//! Assembles a complete molecular topology from MTB building blocks
//! and IFP force field parameters, following the GROMOS convention.
//! Supports both GROMOS @-prefix and standard --prefix argument styles.

use clap::Parser;
use gromos_io::gromos_args;
use gromos_io::ifp::{self, ForceFieldParameters};
use gromos_io::mtb::{self, BbEnd, BbSolute, BlockRef, BuildingBlocks};
use std::collections::HashSet;
use std::process;

/// Build GROMOS topology from building blocks and force field parameters.
///
/// Supports GROMOS @-prefix args: make_top @build 54a7.mtb @param 54a7.ifp @seq NH3+ ALA COO- @solv H2O
#[derive(Parser)]
#[command(name = "make_top", version, about)]
struct Args {
    /// Building block (MTB) file
    #[arg(long)]
    build: String,

    /// Force field parameter (IFP) file
    #[arg(long)]
    param: String,

    /// Residue sequence (end groups + building blocks)
    #[arg(long, num_args = 1..)]
    seq: Vec<String>,

    /// Solvent building block name (e.g. H2O)
    #[arg(long)]
    solv: Option<String>,
}

/// Assembled atom in the final topology.
struct TopAtom {
    name: String,
    residue_nr: usize,
    residue_name: String,
    iac: usize,
    mass: f64,
    charge: f64,
    chargegroup: usize,
    exclusions: Vec<usize>,   // 0-based global indices
    one_four_pairs: Vec<usize>, // 0-based global indices
}

/// Assembled bonded interaction (indices are 0-based global).
struct TopBond {
    i: usize,
    j: usize,
    type_code: usize,
}

struct TopAngle {
    i: usize,
    j: usize,
    k: usize,
    type_code: usize,
}

struct TopImproper {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    type_code: usize,
}

struct TopDihedral {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    type_code: i32,
}

fn main() {
    let args = Args::parse_from(gromos_args());

    let solvent_name = args.solv.as_deref();

    // Read building blocks
    let bb = match mtb::read_mtb_file(&args.build) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("Error reading MTB file '{}': {:?}", args.build, e);
            process::exit(1);
        }
    };

    // Read force field parameters
    let ffp = match ifp::read_ifp_file(&args.param) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Error reading IFP file '{}': {:?}", args.param, e);
            process::exit(1);
        }
    };

    eprintln!("Building topology for {} residues", args.seq.len());
    eprintln!("Force field: {}", ffp.force_field);

    // Assemble topology
    let mut atoms: Vec<TopAtom> = Vec::new();
    let mut bonds: Vec<TopBond> = Vec::new();
    let mut angles: Vec<TopAngle> = Vec::new();
    let mut impropers: Vec<TopImproper> = Vec::new();
    let mut dihedrals: Vec<TopDihedral> = Vec::new();
    let mut residue_names: Vec<String> = Vec::new();

    let mut rep: i32 = 0; // overlap atoms from previous residue

    for (seq_idx, res_name) in args.seq.iter().enumerate() {
        match bb.find_block(res_name) {
            Some(BlockRef::End(end_bb)) => {
                if end_bb.n_replace > 0 {
                    // N-terminal / beginning group
                    add_begin(&mut atoms, &mut bonds, &mut angles, &mut impropers,
                              &mut dihedrals, &mut residue_names, end_bb, &ffp, &bb);
                    rep = end_bb.n_replace;
                } else {
                    // C-terminal / ending group
                    add_end(&mut atoms, &mut bonds, &mut angles, &mut impropers,
                            &mut dihedrals, end_bb, &ffp, &bb);
                }
            }
            Some(BlockRef::Solute(sol_bb)) => {
                add_solute(&mut atoms, &mut bonds, &mut angles, &mut impropers,
                           &mut dihedrals, &mut residue_names, sol_bb, &ffp, &bb, rep);
                rep = sol_bb.num_preceding_exclusions as i32;
            }
            None => {
                eprintln!("Error: Building block '{}' not found in MTB file (seq pos {})", res_name, seq_idx + 1);
                process::exit(1);
            }
        }
    }

    // Build exclusion lists from bonded connectivity
    build_exclusions(&mut atoms, &bonds, bb.link_exclusions);

    // Build 1-4 pair lists
    build_14_pairs(&mut atoms, &bonds, &angles);

    // Write topology to stdout
    write_topology(&atoms, &bonds, &angles, &impropers, &dihedrals,
                   &residue_names, &bb, &ffp, solvent_name);
}

/// Add a beginning (N-terminal) end group. All atoms are added directly.
fn add_begin(
    atoms: &mut Vec<TopAtom>,
    bonds: &mut Vec<TopBond>,
    angles: &mut Vec<TopAngle>,
    impropers: &mut Vec<TopImproper>,
    dihedrals: &mut Vec<TopDihedral>,
    residue_names: &mut Vec<String>,
    end_bb: &BbEnd,
    ffp: &ForceFieldParameters,
    _bb: &BuildingBlocks,
) {
    let offset = atoms.len();
    let res_nr = residue_names.len() + 1;
    residue_names.push(end_bb.name.clone());

    for bb_atom in &end_bb.atoms {
        let mass = ffp.get_mass(bb_atom.mass_code).unwrap_or(0.0);
        atoms.push(TopAtom {
            name: bb_atom.name.clone(),
            residue_nr: res_nr,
            residue_name: end_bb.name.clone(),
            iac: bb_atom.iac as usize,
            mass,
            charge: bb_atom.charge,
            chargegroup: bb_atom.chargegroup,
            exclusions: Vec::new(),
            one_four_pairs: Vec::new(),
        });
    }

    add_bonded_from_end(bonds, angles, impropers, dihedrals, end_bb, offset);
}

/// Add a solute residue, handling overlap with previous residue.
fn add_solute(
    atoms: &mut Vec<TopAtom>,
    bonds: &mut Vec<TopBond>,
    angles: &mut Vec<TopAngle>,
    impropers: &mut Vec<TopImproper>,
    dihedrals: &mut Vec<TopDihedral>,
    residue_names: &mut Vec<String>,
    sol_bb: &BbSolute,
    ffp: &ForceFieldParameters,
    _bb: &BuildingBlocks,
    rep: i32,
) {
    let rep_u = rep.unsigned_abs() as usize;
    // offset: global index of BB atom 1 is at (atoms.len() - rep_u)
    let offset = atoms.len().wrapping_sub(rep_u);
    let res_nr = residue_names.len() + 1;
    residue_names.push(sol_bb.name.clone());

    // Skip the first `rep` atoms (they overlap with previous residue)
    for bb_atom in sol_bb.atoms.iter().skip(rep_u) {
        let mass = ffp.get_mass(bb_atom.mass_code).unwrap_or(0.0);
        atoms.push(TopAtom {
            name: bb_atom.name.clone(),
            residue_nr: res_nr,
            residue_name: sol_bb.name.clone(),
            iac: bb_atom.iac as usize,
            mass,
            charge: bb_atom.charge,
            chargegroup: bb_atom.chargegroup,
            exclusions: Vec::new(),
            one_four_pairs: Vec::new(),
        });
    }

    // Add bonds (BB atom indices are 1-based; negative = previous residue link handled by offset)
    for &(ai, aj, btype) in &sol_bb.bonds {
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        bonds.push(TopBond { i: gi, j: gj, type_code: btype });
    }

    for &(ai, aj, ak, atype) in &sol_bb.angles {
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        let gk = resolve_bb_index(ak, offset);
        if gi == gj || gi == gk || gj == gk {
            eprintln!("Warning: skipping invalid angle with duplicate atoms: {}-{}-{}", gi+1, gj+1, gk+1);
            continue;
        }
        angles.push(TopAngle { i: gi, j: gj, k: gk, type_code: atype });
    }

    for &(ai, aj, ak, al, itype) in &sol_bb.improper_dihedrals {
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        let gk = resolve_bb_index(ak, offset);
        let gl = resolve_bb_index(al, offset);
        impropers.push(TopImproper { i: gi, j: gj, k: gk, l: gl, type_code: itype });
    }

    for &(ai, aj, ak, al, dtype) in &sol_bb.proper_dihedrals {
        if dtype < 0 {
            // Negative type = deletion marker, skip
            continue;
        }
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        let gk = resolve_bb_index(ak, offset);
        let gl = resolve_bb_index(al, offset);
        dihedrals.push(TopDihedral { i: gi, j: gj, k: gk, l: gl, type_code: dtype });
    }
}

/// Add a C-terminal end group. Replaces last |n_replace| atoms.
fn add_end(
    atoms: &mut Vec<TopAtom>,
    bonds: &mut Vec<TopBond>,
    angles: &mut Vec<TopAngle>,
    impropers: &mut Vec<TopImproper>,
    dihedrals: &mut Vec<TopDihedral>,
    end_bb: &BbEnd,
    ffp: &ForceFieldParameters,
    _bb: &BuildingBlocks,
) {
    let nrep = end_bb.n_replace.unsigned_abs() as usize;
    // The first |nrep| atoms in the end group REPLACE the last |nrep| atoms of the molecule
    let replace_start = atoms.len() - nrep;

    // Replace atoms: update IAC, mass, charge for the overlapping atoms
    for (idx, bb_atom) in end_bb.atoms.iter().enumerate() {
        if idx < nrep {
            // Check for iac == -2: charge-only update (GROMOS convention)
            if bb_atom.iac == -2 {
                atoms[replace_start + idx].charge += bb_atom.charge;
            } else {
                let mass = ffp.get_mass(bb_atom.mass_code).unwrap_or(0.0);
                atoms[replace_start + idx].iac = bb_atom.iac as usize;
                atoms[replace_start + idx].mass = mass;
                atoms[replace_start + idx].charge = bb_atom.charge;
                atoms[replace_start + idx].chargegroup = bb_atom.chargegroup;
            }
        } else {
            // New atoms beyond the overlap
            let mass = ffp.get_mass(bb_atom.mass_code).unwrap_or(0.0);
            atoms.push(TopAtom {
                name: bb_atom.name.clone(),
                residue_nr: atoms.last().map_or(1, |a| a.residue_nr),
                residue_name: end_bb.name.clone(),
                iac: bb_atom.iac as usize,
                mass,
                charge: bb_atom.charge,
                chargegroup: bb_atom.chargegroup,
                exclusions: Vec::new(),
                one_four_pairs: Vec::new(),
            });
        }
    }

    // offset for end group: BB atom 1 maps to replace_start
    let offset = replace_start;
    // But BB indices are 1-based, and negative indices refer backward
    // For end groups, all indices should be positive and refer to atoms within the group
    add_bonded_from_end(bonds, angles, impropers, dihedrals, end_bb, offset);
}

fn add_bonded_from_end(
    bonds: &mut Vec<TopBond>,
    angles: &mut Vec<TopAngle>,
    impropers: &mut Vec<TopImproper>,
    dihedrals: &mut Vec<TopDihedral>,
    end_bb: &BbEnd,
    offset: usize,
) {
    for &(ai, aj, btype) in &end_bb.bonds {
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        bonds.push(TopBond { i: gi, j: gj, type_code: btype });
    }

    for &(ai, aj, ak, atype) in &end_bb.angles {
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        let gk = resolve_bb_index(ak, offset);
        if gi == gj || gi == gk || gj == gk {
            eprintln!("Warning: skipping invalid end-group angle with duplicate atoms: {}-{}-{}", gi+1, gj+1, gk+1);
            continue;
        }
        angles.push(TopAngle { i: gi, j: gj, k: gk, type_code: atype });
    }

    for &(ai, aj, ak, al, itype) in &end_bb.improper_dihedrals {
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        let gk = resolve_bb_index(ak, offset);
        let gl = resolve_bb_index(al, offset);
        impropers.push(TopImproper { i: gi, j: gj, k: gk, l: gl, type_code: itype });
    }

    for &(ai, aj, ak, al, dtype) in &end_bb.proper_dihedrals {
        if dtype < 0 {
            continue;
        }
        let gi = resolve_bb_index(ai, offset);
        let gj = resolve_bb_index(aj, offset);
        let gk = resolve_bb_index(ak, offset);
        let gl = resolve_bb_index(al, offset);
        dihedrals.push(TopDihedral { i: gi, j: gj, k: gk, l: gl, type_code: dtype });
    }
}

/// Convert a BB atom index (1-based, can be negative for previous residue) to 0-based global index.
fn resolve_bb_index(bb_idx: i32, offset: usize) -> usize {
    if bb_idx > 0 {
        offset + (bb_idx as usize) - 1
    } else {
        // Negative index: refers to atoms before this residue
        // bb_idx=-1 means the last atom before offset, etc.
        (offset as i32 + bb_idx) as usize
    }
}

/// Build exclusion lists from connectivity (nearest-neighbour exclusions).
fn build_exclusions(atoms: &mut Vec<TopAtom>, bonds: &[TopBond], link_exclusions: usize) {
    let n = atoms.len();

    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for bond in bonds {
        if bond.i < n && bond.j < n {
            adj[bond.i].push(bond.j);
            adj[bond.j].push(bond.i);
        }
    }

    // For each atom, find all atoms within `link_exclusions` bonds
    for i in 0..n {
        let mut excluded = HashSet::new();
        let mut frontier = vec![i];
        let mut visited = HashSet::new();
        visited.insert(i);

        for _depth in 0..link_exclusions {
            let mut next_frontier = Vec::new();
            for &node in &frontier {
                for &neighbor in &adj[node] {
                    if visited.insert(neighbor) {
                        next_frontier.push(neighbor);
                        if neighbor > i {
                            excluded.insert(neighbor);
                        }
                    }
                }
            }
            frontier = next_frontier;
        }

        let mut excl_vec: Vec<usize> = excluded.into_iter().collect();
        excl_vec.sort();
        atoms[i].exclusions = excl_vec;
    }
}

/// Build 1-4 pair lists (atoms separated by exactly 3 bonds, not already excluded).
fn build_14_pairs(atoms: &mut Vec<TopAtom>, bonds: &[TopBond], _angles: &[TopAngle]) {
    let n = atoms.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for bond in bonds {
        if bond.i < n && bond.j < n {
            adj[bond.i].push(bond.j);
            adj[bond.j].push(bond.i);
        }
    }

    for i in 0..n {
        let mut pairs_14 = HashSet::new();
        // Walk 3 bonds from i
        for &j in &adj[i] {
            for &k in &adj[j] {
                if k == i { continue; }
                for &l in &adj[k] {
                    if l == j || l == i { continue; }
                    if l > i && !atoms[i].exclusions.contains(&l) {
                        pairs_14.insert(l);
                    }
                }
            }
        }
        let mut pairs: Vec<usize> = pairs_14.into_iter().collect();
        pairs.sort();
        atoms[i].one_four_pairs = pairs;
    }
}

/// Write the complete topology to stdout.
fn write_topology(
    atoms: &[TopAtom],
    bonds: &[TopBond],
    angles: &[TopAngle],
    impropers: &[TopImproper],
    dihedrals: &[TopDihedral],
    residue_names: &[String],
    bb: &BuildingBlocks,
    ffp: &ForceFieldParameters,
    solvent_name: Option<&str>,
) {
    // TITLE
    println!("TITLE");
    println!("    Topology generated by make_top (gromos-rs)");
    println!("    Force field: {}", ffp.force_field);
    println!("END");

    // PHYSICALCONSTANTS
    println!("PHYSICALCONSTANTS");
    println!("# FPEPSI");
    println!("{:15.4}", bb.fpepsi);
    println!("# HBAR");
    println!("{:15.7}", bb.hbar);
    println!("# SPDL");
    println!("{:15.3}", bb.spdl);
    println!("# BOLTZ");
    println!("{:15.8}", bb.boltz);
    println!("END");

    // TOPVERSION
    println!("TOPVERSION");
    println!("2.0");
    println!("END");

    // ATOMTYPENAME
    println!("ATOMTYPENAME");
    println!("# NRATT: number of atom types");
    println!("{:5}", ffp.num_atom_types);
    println!("# TYPE: atom type names");
    for alj in &ffp.atom_lj {
        println!("{}", alj.name);
    }
    println!("END");

    // RESNAME
    println!("RESNAME");
    println!("# NRAA2: number of residue names");
    println!("{:5}", residue_names.len());
    println!("# AANM: residue names");
    for name in residue_names {
        println!("{}", name);
    }
    println!("END");

    // SOLUTEATOM
    println!("SOLUTEATOM");
    println!("#   NRP: number of solute atoms");
    println!("{:6}", atoms.len());
    println!("#  ATNM: atom number");
    println!("#  MRES: residue number");
    println!("#  PANM: atom name");
    println!("#   IAC: integer atom code");
    println!("#  MASS: mass");
    println!("#    CG: charge");
    println!("#   CGC: charge group code");
    println!("#   INE: number of excluded atoms");
    println!("# INE14: number of 1-4 interactions");
    println!("#  ATNM  MRES  PANM   IAC      MASS        CG  CGC   INE");
    println!("#                                                     INE14");

    for (i, atom) in atoms.iter().enumerate() {
        let ne = atom.exclusions.len();
        println!(
            "{:6}{:5} {:>5}{:5}{:10.5}{:12.5}{:3}{:4}",
            i + 1,
            atom.residue_nr,
            atom.name,
            atom.iac,
            atom.mass,
            atom.charge,
            atom.chargegroup,
            ne
        );
        // Print exclusion list
        if !atom.exclusions.is_empty() {
            for (j, &excl) in atom.exclusions.iter().enumerate() {
                print!("{:5}", excl + 1);
                if (j + 1) % 20 == 0 {
                    println!();
                }
            }
            if atom.exclusions.len() % 20 != 0 {
                println!();
            }
        }
        // Print 1-4 pairs
        let n14 = atom.one_four_pairs.len();
        println!("   {:4}", n14);
        if !atom.one_four_pairs.is_empty() {
            for (j, &pair) in atom.one_four_pairs.iter().enumerate() {
                print!("{:5}", pair + 1);
                if (j + 1) % 20 == 0 {
                    println!();
                }
            }
            if atom.one_four_pairs.len() % 20 != 0 {
                println!();
            }
        }
    }
    println!("END");

    // BONDSTRETCHTYPE
    println!("BONDSTRETCHTYPE");
    println!("# NBTY: number of bond types");
    println!("{:5}", ffp.bond_types.len());
    println!("#  CB     HB     B0");
    for bt in &ffp.bond_types {
        println!(
            "{:15.7e}{:15.7e}{:13.7e}",
            bt.k_quartic, bt.k_harmonic, bt.r0
        );
    }
    println!("END");

    // Separate bonds into BONDH (involving H) and BOND (not involving H)
    // Convention: if either atom has mass < 2.0 (hydrogen), it's BONDH
    let (bondh, bond): (Vec<&TopBond>, Vec<&TopBond>) = bonds
        .iter()
        .partition(|b| atoms[b.i].mass < 2.0 || atoms[b.j].mass < 2.0);

    println!("BONDH");
    println!("# NBAH: number of bonds involving H atoms");
    println!("{:5}", bondh.len());
    println!("#  IB   JB  ICB");
    for b in &bondh {
        println!("{:5}{:5}{:5}", b.i + 1, b.j + 1, b.type_code);
    }
    println!("END");

    println!("BOND");
    println!("# NBA: number of bonds not involving H atoms");
    println!("{:5}", bond.len());
    println!("#  IB   JB  ICB");
    for b in &bond {
        println!("{:5}{:5}{:5}", b.i + 1, b.j + 1, b.type_code);
    }
    println!("END");

    // BONDANGLEBENDTYPE
    println!("BONDANGLEBENDTYPE");
    println!("# NTTY: number of angle types");
    println!("{:5}", ffp.angle_types.len());
    println!("#   CT     CHT     T0");
    for at in &ffp.angle_types {
        println!(
            "{:15.7e}{:15.7e}{:13.7e}",
            at.k_non_harmonic, at.k_harmonic, at.theta0
        );
    }
    println!("END");

    // Separate angles into BONDANGLEH and BONDANGLE
    let (angleh, angle): (Vec<&TopAngle>, Vec<&TopAngle>) = angles
        .iter()
        .partition(|a| atoms[a.i].mass < 2.0 || atoms[a.j].mass < 2.0 || atoms[a.k].mass < 2.0);

    println!("BONDANGLEH");
    println!("# NTHEH: number of angles involving H atoms");
    println!("{:5}", angleh.len());
    println!("#  IT   JT   KT  ICT");
    for a in &angleh {
        println!("{:5}{:5}{:5}{:5}", a.i + 1, a.j + 1, a.k + 1, a.type_code);
    }
    println!("END");

    println!("BONDANGLE");
    println!("# NTHE: number of angles not involving H atoms");
    println!("{:5}", angle.len());
    println!("#  IT   JT   KT  ICT");
    for a in &angle {
        println!("{:5}{:5}{:5}{:5}", a.i + 1, a.j + 1, a.k + 1, a.type_code);
    }
    println!("END");

    // IMPDIHEDRALTYPE
    println!("IMPDIHEDRALTYPECODE");
    println!("# NQTY: number of improper dihedral types");
    println!("{:5}", ffp.improper_types.len());
    println!("#  CQ     Q0");
    for it in &ffp.improper_types {
        println!("{:15.7e}{:13.4}", it.k, it.xi0);
    }
    println!("END");

    // Separate impropers into H / no-H
    let (impdihedralh, impdihedral): (Vec<&TopImproper>, Vec<&TopImproper>) = impropers
        .iter()
        .partition(|d| {
            atoms[d.i].mass < 2.0
                || atoms[d.j].mass < 2.0
                || atoms[d.k].mass < 2.0
                || atoms[d.l].mass < 2.0
        });

    println!("IMPDIHEDRALH");
    println!("# NQHIH: number of improper dihedrals involving H");
    println!("{:5}", impdihedralh.len());
    println!("#  IQ   JQ   KQ   LQ  ICQ");
    for d in &impdihedralh {
        println!(
            "{:5}{:5}{:5}{:5}{:5}",
            d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.type_code
        );
    }
    println!("END");

    println!("IMPDIHEDRAL");
    println!("# NQHI: number of improper dihedrals not involving H");
    println!("{:5}", impdihedral.len());
    println!("#  IQ   JQ   KQ   LQ  ICQ");
    for d in &impdihedral {
        println!(
            "{:5}{:5}{:5}{:5}{:5}",
            d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.type_code
        );
    }
    println!("END");

    // TORSDIHEDRALTYPECODE
    println!("TORSDIHEDRALTYPECODE");
    println!("# NPTY: number of dihedral types");
    println!("{:5}", ffp.torsion_types.len());
    println!("#   CP     PD    NP");
    for tt in &ffp.torsion_types {
        println!("{:15.7e}{:10.1}{:5}", tt.k, tt.delta, tt.multiplicity);
    }
    println!("END");

    // Separate dihedrals into H / no-H
    let (dihedralh, dihedral): (Vec<&TopDihedral>, Vec<&TopDihedral>) = dihedrals
        .iter()
        .partition(|d| {
            atoms[d.i].mass < 2.0
                || atoms[d.j].mass < 2.0
                || atoms[d.k].mass < 2.0
                || atoms[d.l].mass < 2.0
        });

    println!("DIHEDRALH");
    println!("# NPHIH: number of dihedrals involving H");
    println!("{:5}", dihedralh.len());
    println!("#  IP   JP   KP   LP  ICP");
    for d in &dihedralh {
        println!(
            "{:5}{:5}{:5}{:5}{:5}",
            d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.type_code
        );
    }
    println!("END");

    println!("DIHEDRAL");
    println!("# NPHI: number of dihedrals not involving H");
    println!("{:5}", dihedral.len());
    println!("#  IP   JP   KP   LP  ICP");
    for d in &dihedral {
        println!(
            "{:5}{:5}{:5}{:5}{:5}",
            d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.type_code
        );
    }
    println!("END");

    // LJPARAMETERS
    println!("LJPARAMETERS");
    println!("# NRATT2: number of LJ parameter pairs");
    let n_types = ffp.num_atom_types;
    let n_pairs = n_types * (n_types + 1) / 2;
    println!("{:5}", n_pairs);
    println!("#  IAC  JAC    C12          C6           CS12         CS6");
    for i in 1..=n_types {
        for j in i..=n_types {
            let c6 = ffp.compute_c6(i, j);
            let c12 = ffp.compute_c12(i, j);
            let cs6 = ffp.compute_c6_14(i, j);
            let cs12 = ffp.compute_c12_14(i, j);
            println!(
                "{:5}{:5}{:15.7e}{:15.7e}{:15.7e}{:15.7e}",
                i, j, c12, c6, cs12, cs6
            );
        }
    }
    println!("END");

    // SOLUTEMOLECULES
    println!("SOLUTEMOLECULES");
    println!("# NSPM: number of separate molecules");
    println!("    1");
    println!("# NSP: atom sequence number of last atom");
    println!("    {}", atoms.len());
    println!("END");

    // TEMPERATUREGROUPS
    println!("TEMPERATUREGROUPS");
    println!("# NSTM: number of temperature atom groups");
    println!("    1");
    println!("# NST: atom sequence number of last atom");
    println!("    {}", atoms.len());
    println!("END");

    // PRESSUREGROUPS
    println!("PRESSUREGROUPS");
    println!("# NSTM: number of pressure atom groups");
    println!("    1");
    println!("# NST: atom sequence number of last atom");
    println!("    {}", atoms.len());
    println!("END");

    // SOLVENTATOM / SOLVENTCONSTR
    if let Some(sname) = solvent_name {
        if let Some(sol) = bb.find_solvent(sname) {
            println!("SOLVENTATOM");
            println!("# NRAM: number of atoms per solvent molecule");
            println!("{:5}", sol.atoms.len());
            println!("#     I  ANMS   IACS  MASS     CGS");
            for (idx, (aname, iac, mass_code, charge)) in sol.atoms.iter().enumerate() {
                let mass = ffp.get_mass(*mass_code).unwrap_or(0.0);
                println!(
                    "{:5} {:>5}{:5}{:10.5}{:12.5}",
                    idx + 1,
                    aname,
                    iac,
                    mass,
                    charge
                );
            }
            println!("END");

            println!("SOLVENTCONSTR");
            println!("# NCONS: number of constraints");
            println!("{:5}", sol.constraints.len());
            println!("#  IC   JC       CC");
            for (ic, jc, dist) in &sol.constraints {
                println!("{:5}{:5}{:15.7e}", ic, jc, dist);
            }
            println!("END");
        }
    }
}
