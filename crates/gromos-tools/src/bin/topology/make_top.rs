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
use gromos_io::mtb::{self, BbAtom, BbEnd, BbSolute, BlockRef, BuildingBlocks};
use std::collections::BTreeSet;
use std::process;

/// Build a GROMOS topology from building blocks and force-field parameters.
///
/// Supports GROMOS @-prefix args: make_top @build 54a7.mtb @param 54a7.ifp @seq NH3+ ALA COO- @solv H2O
#[derive(Parser)]
#[command(name = "make_top", version, about)]
struct Args {
    /// Building block (MTB) files; later files add blocks to the first
    #[arg(long, num_args = 1..)]
    build: Vec<String>,

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

// ── The linear topology gromos++ assembles (`gcore::LinearTopology`) ──────────────────────
//
// Atom indices are 0-based, as in gromos++ after `InBuildingBlock` subtracted one from every
// index in the MTB file (`--i`); negative values are the "previous residue" references of the
// building blocks. IAC −1 marks an atom that a begin group replaces (removed afterwards), IAC −2
// an end-group atom whose charge is only added to the last atom of that name.

#[derive(Clone, Debug)]
struct LtAtom {
    name: String,
    iac: i32,
    mass: f64,
    charge: f64,
    chargegroup: usize,
    exclusions: BTreeSet<i64>,
    one_four: BTreeSet<usize>,
    res: usize,
}

/// gromos++ `std::set<Bond>` order: (i, j, type).
type Bond = (i64, i64, i64);
/// `std::set<Angle>` order: (j, i, k, type) — stored in that key order.
type AngleKey = (i64, i64, i64, i64);
/// `std::set<Improper>` order: (i, j, k, l, type).
type Improper = (i64, i64, i64, i64, i64);
/// `std::set<Dihedral>` order: (j, k, i, l, type) — stored in that key order.
type DihedralKey = (i64, i64, i64, i64, i64);

fn angle_key(i: i64, j: i64, k: i64, t: i64) -> AngleKey {
    (j, i, k, t)
}
fn angle_atoms(a: &AngleKey) -> (i64, i64, i64, i64) {
    (a.1, a.0, a.2, a.3)
}
fn dihedral_key(i: i64, j: i64, k: i64, l: i64, t: i64) -> DihedralKey {
    (j, k, i, l, t)
}
fn dihedral_atoms(d: &DihedralKey) -> (i64, i64, i64, i64, i64) {
    (d.2, d.0, d.1, d.3, d.4)
}

#[derive(Default)]
struct Lt {
    atoms: Vec<LtAtom>,
    res_names: Vec<String>,
    bonds: BTreeSet<Bond>,
    angles: BTreeSet<AngleKey>,
    impropers: BTreeSet<Improper>,
    dihedrals: BTreeSet<DihedralKey>,
}

/// A building block as gromos++ sees it: 0-based indices, resolved masses.
struct Block {
    name: String,
    rep: i32,
    atoms: Vec<LtAtom>,
    /// per-atom exclusions, 0-based relative to the block
    excl: Vec<Vec<i64>>,
    /// preceding exclusions: (atoms before the block, 0-based relative) with their lists
    pexcl: Vec<Vec<i64>>,
    bonds: Vec<(i64, i64, i64)>,
    angles: Vec<(i64, i64, i64, i64)>,
    impropers: Vec<(i64, i64, i64, i64, i64)>,
    dihedrals: Vec<(i64, i64, i64, i64, i64)>,
}

fn bb_atom(a: &BbAtom, ffp: &ForceFieldParameters) -> LtAtom {
    LtAtom {
        name: a.name.clone(),
        iac: a.iac - 1,
        mass: ffp.get_mass(a.mass_code).unwrap_or(0.0),
        charge: a.charge,
        chargegroup: a.chargegroup,
        exclusions: BTreeSet::new(),
        one_four: BTreeSet::new(),
        res: 0,
    }
}

fn from_solute(bb: &BbSolute, ffp: &ForceFieldParameters) -> Block {
    Block {
        name: bb.name.clone(),
        rep: 0,
        atoms: bb.atoms.iter().map(|a| bb_atom(a, ffp)).collect(),
        excl: bb
            .atoms
            .iter()
            .map(|a| a.exclusions.iter().map(|&e| e as i64 - 1).collect())
            .collect(),
        pexcl: bb
            .preceding_exclusions
            .iter()
            .map(|p| p.excluded_atoms.iter().map(|&e| e as i64 - 1).collect())
            .collect(),
        bonds: bb
            .bonds
            .iter()
            .map(|&(i, j, t)| (i as i64 - 1, j as i64 - 1, t as i64))
            .collect(),
        angles: bb
            .angles
            .iter()
            .map(|&(i, j, k, t)| (i as i64 - 1, j as i64 - 1, k as i64 - 1, t as i64))
            .collect(),
        impropers: bb
            .improper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
        dihedrals: bb
            .proper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
    }
}

fn from_end(bb: &BbEnd, ffp: &ForceFieldParameters) -> Block {
    Block {
        name: bb.name.clone(),
        rep: bb.n_replace,
        atoms: bb.atoms.iter().map(|a| bb_atom(a, ffp)).collect(),
        excl: bb
            .atoms
            .iter()
            .map(|a| a.exclusions.iter().map(|&e| e as i64 - 1).collect())
            .collect(),
        pexcl: Vec::new(),
        bonds: bb
            .bonds
            .iter()
            .map(|&(i, j, t)| (i as i64 - 1, j as i64 - 1, t as i64))
            .collect(),
        angles: bb
            .angles
            .iter()
            .map(|&(i, j, k, t)| (i as i64 - 1, j as i64 - 1, k as i64 - 1, t as i64))
            .collect(),
        impropers: bb
            .improper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
        dihedrals: bb
            .proper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
    }
}

impl Lt {
    fn set_res_name(&mut self, resnum: usize, name: &str) {
        if self.res_names.len() <= resnum {
            self.res_names.resize(resnum + 1, String::new());
        }
        self.res_names[resnum] = name.to_string();
    }

    fn has_bond(&self, i: i64, j: i64) -> bool {
        self.bonds.iter().any(|b| b.0 == i && b.1 == j)
    }
    fn has_angle(&self, i: i64, j: i64, k: i64) -> bool {
        self.angles
            .iter()
            .any(|a| angle_atoms(a).0 == i && a.0 == j && a.2 == k)
    }
    fn has_improper(&self, i: i64, j: i64, k: i64, l: i64) -> bool {
        self.impropers
            .iter()
            .any(|p| p.0 == i && p.1 == j && p.2 == k && p.3 == l)
    }

    /// gromos++ `bondedAtoms`: atoms bonded (as the first partner) to one of `atoms`, at an
    /// index not above `offset` and not in `atoms` — the previous residue's link candidates.
    fn bonded_atoms(&self, atoms: &BTreeSet<i64>, offset: i64) -> BTreeSet<i64> {
        let min = atoms.iter().copied().fold(offset, i64::min);
        let mut out = BTreeSet::new();
        for &a in atoms {
            for b in &self.bonds {
                if b.1 == a && b.0 <= min && !atoms.contains(&b.0) {
                    out.insert(b.0);
                }
            }
        }
        out
    }

    /// gromos++ `LinearTopology::removeAtoms`: drop atoms with a negative IAC (replaced by a
    /// begin group) and renumber everything.
    fn remove_atoms(&mut self) {
        let removed: BTreeSet<usize> = self
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.iac < 0)
            .map(|(i, _)| i)
            .collect();
        if removed.is_empty() {
            return;
        }
        let mut ren = vec![-1i64; self.atoms.len()];
        let mut corr = 0usize;
        for (i, r) in ren.iter_mut().enumerate() {
            if removed.contains(&i) {
                corr += 1;
            } else {
                *r = (i - corr) as i64;
            }
        }
        let map = |x: i64| -> Option<i64> {
            if x < 0 || x as usize >= ren.len() {
                Some(x) // out of range: keep (later checks drop it)
            } else if ren[x as usize] < 0 {
                None
            } else {
                Some(ren[x as usize])
            }
        };
        let mut atoms = Vec::new();
        for (i, a) in self.atoms.iter().enumerate() {
            if removed.contains(&i) {
                continue;
            }
            let mut a = a.clone();
            a.exclusions = a.exclusions.iter().filter_map(|&e| map(e)).collect();
            atoms.push(a);
        }
        self.atoms = atoms;
        self.bonds = self
            .bonds
            .iter()
            .filter_map(|b| Some((map(b.0)?, map(b.1)?, b.2)))
            .collect();
        self.angles = self
            .angles
            .iter()
            .filter_map(|a| Some((map(a.0)?, map(a.1)?, map(a.2)?, a.3)))
            .collect();
        self.impropers = self
            .impropers
            .iter()
            .filter_map(|p| Some((map(p.0)?, map(p.1)?, map(p.2)?, map(p.3)?, p.4)))
            .collect();
        self.dihedrals = self
            .dihedrals
            .iter()
            .filter_map(|d| Some((map(d.0)?, map(d.1)?, map(d.2)?, map(d.3)?, d.4)))
            .collect();
    }

    /// gromos++ `LinearTopology::get14s`: atoms exactly three bonds apart (and not closer)
    /// that are not excluded are 1-4 partners.
    fn get14s(&mut self) {
        let n = self.atoms.len();
        let mut neigh: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
        for b in &self.bonds {
            if b.0 >= 0 && b.1 >= 0 && (b.0 as usize) < n && (b.1 as usize) < n {
                neigh[b.0 as usize].insert(b.1 as usize);
                neigh[b.1 as usize].insert(b.0 as usize);
            }
        }
        for i in 0..n {
            let first = neigh[i].clone();
            let mut second = BTreeSet::new();
            for &f in &first {
                second.extend(neigh[f].iter().copied());
            }
            let mut third = BTreeSet::new();
            for &s in &second {
                third.extend(neigh[s].iter().copied());
            }
            let mut e = BTreeSet::new();
            for &t in &third {
                if i < t
                    && !first.contains(&t)
                    && !second.contains(&t)
                    && !self.atoms[i].exclusions.contains(&(t as i64))
                {
                    e.insert(t);
                }
            }
            self.atoms[i].one_four = e;
        }
    }
}

/// gromos++ `addSolute`.
fn add_solute(
    lt: &mut Lt,
    bb: &Block,
    resnum: usize,
    resname: &str,
    rep: i64,
    nn: i64,
) -> Result<(), String> {
    let na = lt.atoms.len() as i64;
    let strt = na - rep;
    let beg = if strt < 0 { -strt } else { 0 };
    for i in beg..rep {
        let e: BTreeSet<i64> = bb.excl[i as usize].iter().map(|x| x + strt).collect();
        lt.atoms[(strt + i) as usize].exclusions = e;
    }
    if rep == 0 {
        let pexl = strt - bb.pexcl.len() as i64;
        if pexl < 0 {
            return Err(format!(
                "{resname}: preceding exclusions, but no preceding atoms"
            ));
        }
        for (i, list) in bb.pexcl.iter().enumerate() {
            let e: BTreeSet<i64> = list.iter().map(|x| x + strt).collect();
            lt.atoms[pexl as usize + i].exclusions = e;
        }
    }
    for i in rep as usize..bb.atoms.len() {
        let mut a = bb.atoms[i].clone();
        a.exclusions = bb.excl[i].iter().map(|x| x + strt).collect();
        a.res = resnum;
        lt.atoms.push(a);
    }
    lt.set_res_name(resnum, resname);
    let offset = strt;
    for &(i, j, t) in &bb.bonds {
        let b = (i + offset, j + offset, t);
        if !(rep > 0 && lt.has_bond(b.0, b.1)) {
            lt.bonds.insert(b);
        }
    }
    for &(i, j, k, t) in &bb.angles {
        let (i, j, k) = (i + offset, j + offset, k + offset);
        let found = rep > 0 && lt.has_angle(i, j, k);
        if (i >= na || j >= na || k >= na) && !found {
            lt.angles.insert(angle_key(i, j, k, t));
        }
    }
    for &(i, j, k, l, t) in &bb.impropers {
        let (i, j, k, l) = (i + offset, j + offset, k + offset, l + offset);
        let found = rep > 0 && lt.has_improper(i, j, k, l);
        if (i >= na || j >= na || k >= na || l >= na) && !found {
            lt.impropers.insert((i, j, k, l, t));
        }
    }
    if rep == 0 {
        for &(d0, d1, d2, d3, t) in &bb.dihedrals {
            let (mut corr0, corr1, corr2, mut corr3) = (offset, offset, offset, offset);
            if d3 == -2 {
                corr3 = 0;
            }
            if d0 == -3 {
                corr0 = 0;
                for b in &lt.bonds {
                    if b.1 == d1 + offset {
                        corr0 = b.0 + 3;
                    }
                }
            }
            lt.dihedrals.insert(dihedral_key(
                d0 + corr0,
                d1 + corr1,
                d2 + corr2,
                d3 + corr3,
                t,
            ));
        }
    } else if rep > 0 {
        let mut old: Vec<(i64, i64, i64, i64, i64)> =
            lt.dihedrals.iter().map(dihedral_atoms).collect();
        lt.dihedrals.clear();
        for &(d0, d1, d2, d3, t) in &bb.dihedrals {
            let mut corr = offset;
            if d0 == -3 {
                for b in &lt.bonds {
                    if b.1 == d1 + offset {
                        corr = b.0 + 3;
                    }
                }
            }
            if d0 == -2 {
                for b in &lt.bonds {
                    if b.1 == d1 + offset {
                        corr = b.0 + 2;
                    }
                }
            }
            let b = (d0 + corr, d1 + offset, d2 + offset, d3 + offset, t);
            let mut add = true;
            for o in old.iter_mut() {
                if o.0 == b.0 && o.1 == b.1 && o.2 == b.2 {
                    o.3 = b.3;
                    add = false;
                }
            }
            if b.0 < nn || b.1 < nn || b.2 < nn || b.3 < nn {
                add = false;
            }
            if b.1 == -6 + offset && b.2 == -5 + offset && b.3 == -4 + offset {
                add = true;
            }
            if b.0 == -5 + offset && b.1 == -6 + offset {
                add = true;
            }
            if b.3 == -6 + offset {
                add = true;
            }
            if add {
                lt.dihedrals.insert(dihedral_key(b.0, b.1, b.2, b.3, b.4));
            }
        }
        for o in old {
            lt.dihedrals.insert(dihedral_key(o.0, o.1, o.2, o.3, o.4));
        }
    }
    Ok(())
}

/// gromos++ `addBegin`: the atoms of a begin group, exclusions shifted by the atoms so far.
fn add_begin(lt: &mut Lt, bb: &Block, resnum: usize) -> i64 {
    let na = lt.atoms.len() as i64;
    for (i, a) in bb.atoms.iter().enumerate() {
        let mut a = a.clone();
        a.exclusions = bb.excl[i].iter().map(|x| x + na).collect();
        a.res = resnum;
        lt.atoms.push(a);
    }
    bb.rep as i64
}

/// gromos++ `addEnd`: replace the last |rep| atoms of the chain by the end group.
fn add_end(lt: &mut Lt, bb: &Block, resnum: usize) -> Result<(), String> {
    let strt = lt.atoms.len() as i64 + bb.rep as i64;
    let mut search = Vec::new();
    for (i, a) in bb.atoms.iter().enumerate() {
        if a.iac < -2 {
            return Err(format!("IAC < -1 found in end group {}", bb.name));
        }
        if a.iac == -2 {
            search.push(i);
            eprintln!(
                "# WARNING: for atom {} in MTBUILDBLEND group {} only the CHARGE is transferred to the last atom with this name in the chain.",
                a.name, bb.name
            );
        }
    }
    let n_pop = (-bb.rep) as usize - search.len();
    for _ in 0..n_pop {
        lt.atoms.pop();
    }
    for &s in &search {
        let mut j = lt.atoms.len() - 1;
        while lt.atoms[j].name != bb.atoms[s].name {
            j -= 1;
        }
        lt.atoms[j].charge += bb.atoms[s].charge;
    }
    for (i, a) in bb.atoms.iter().enumerate() {
        if a.iac >= -1 {
            let mut a = a.clone();
            a.exclusions = bb.excl[i].iter().map(|x| x + strt).collect();
            a.res = resnum;
            lt.atoms.push(a);
        }
    }
    Ok(())
}

/// gromos++ `addCovEnd`: the covalent terms of an end group; for a terminal group the
/// references to the previous residue are resolved through the bond graph.
fn add_cov_end(lt: &mut Lt, bb: &Block, offset: i64) -> Result<(), String> {
    let name = &bb.name;
    for &(i, j, t) in &bb.bonds {
        let mut b = (i + offset, j + offset, t);
        if bb.rep < 0 {
            if i < 0 || bb.atoms[i as usize].iac == -2 {
                let atoms: BTreeSet<i64> = [j + offset].into_iter().collect();
                let cand = lt.bonded_atoms(&atoms, offset);
                if cand.is_empty() {
                    return Err(format!(
                        "bond to the previous residue in {name} is not found"
                    ));
                }
                if cand.len() != 1 {
                    return Err(format!(
                        "bond to the previous residue in {name} is ambiguous"
                    ));
                }
                b.0 = *cand.iter().next().unwrap();
            }
            lt.bonds.retain(|x| !(x.0 == b.0 && x.1 == b.1));
            lt.bonds.insert(b);
        } else {
            lt.bonds.insert(b);
        }
    }
    for &(i, j, k, t) in &bb.angles {
        let mut a = (i + offset, j + offset, k + offset, t);
        if bb.rep < 0 && (i < 0 || bb.atoms[i as usize].iac == -2) {
            let atoms: BTreeSet<i64> = [j + offset].into_iter().collect();
            let cand = lt.bonded_atoms(&atoms, offset);
            if cand.is_empty() {
                return Err(format!(
                    "angle to the previous residue in {name} is not found"
                ));
            }
            if cand.len() != 1 {
                return Err(format!(
                    "angle to the previous residue in {name} is ambiguous"
                ));
            }
            a.0 = *cand.iter().next().unwrap();
        }
        lt.angles.retain(|x| {
            let (xi, xj, xk, _) = angle_atoms(x);
            !(xi == a.0 && xj == a.1 && xk == a.2)
        });
        lt.angles.insert(angle_key(a.0, a.1, a.2, a.3));
    }
    for &(i, j, k, l, t) in &bb.impropers {
        let mut list = [i + offset, j + offset, k + offset, l + offset];
        if bb.rep < 0 {
            let raw = [i, j, k, l];
            let mut atoms = BTreeSet::new();
            let mut negs = Vec::new();
            for (n, &r) in raw.iter().enumerate() {
                if r >= 0 && bb.atoms[r as usize].iac >= -1 {
                    atoms.insert(r + offset);
                } else {
                    negs.push(n);
                }
            }
            if atoms.len() < 4 {
                let cand = lt.bonded_atoms(&atoms, offset);
                if cand.is_empty() {
                    return Err(format!(
                        "improper to the previous residue in {name} is not found"
                    ));
                }
                if cand.len() != 1 {
                    return Err(format!(
                        "improper to the previous residue in {name} is ambiguous"
                    ));
                }
                if negs.len() > 1 {
                    return Err(format!("multiple negative values in {name} impropers"));
                }
                list[negs[0]] = *cand.iter().next().unwrap();
            }
        }
        let set: BTreeSet<i64> = list.iter().copied().collect();
        lt.impropers.retain(|x| {
            !(set.contains(&x.0) && set.contains(&x.1) && set.contains(&x.2) && set.contains(&x.3))
        });
        lt.impropers.insert((list[0], list[1], list[2], list[3], t));
    }
    let mut added: BTreeSet<DihedralKey> = BTreeSet::new();
    for &(i, j, k, l, t) in &bb.dihedrals {
        let mut list = [i + offset, j + offset, k + offset, l + offset];
        if bb.rep < 0 {
            let raw = [i, j, k, l];
            let mut atoms = BTreeSet::new();
            let mut negs = Vec::new();
            for (n, &r) in raw.iter().enumerate() {
                if r >= 0 && bb.atoms[r as usize].iac >= -1 {
                    atoms.insert(r + offset);
                } else {
                    negs.push(n);
                }
            }
            if atoms.len() < 4 {
                let cand = lt.bonded_atoms(&atoms, offset);
                if cand.is_empty() {
                    return Err(format!(
                        "dihedral to the previous residue in {name} is not found"
                    ));
                }
                if cand.len() != 1 {
                    eprintln!("Warning: dihedral to the previous residue in {name} is ambiguous");
                }
                if negs.len() > 1 {
                    return Err(format!("multiple negative values in {name} dihedrals"));
                }
                list[negs[0]] = *cand.iter().next().unwrap();
            }
        }
        let key = dihedral_key(list[0], list[1], list[2], list[3], t);
        lt.dihedrals.retain(|x| {
            let (xi, xj, xk, xl, _) = dihedral_atoms(x);
            !(xi == list[0] && xj == list[1] && xk == list[2] && xl == list[3]) || added.contains(x)
        });
        lt.dihedrals.insert(key);
        added.insert(key);
    }
    Ok(())
}

fn is_h(mass: f64) -> bool {
    (mass - 1.008).abs() < 1e-9
}

fn main() {
    let args = Args::parse_from(gromos_args());
    let solvent_name = args.solv.as_deref();
    let mut bb: Option<BuildingBlocks> = None;
    for path in &args.build {
        match mtb::read_mtb_file(path) {
            Ok(b) => match bb.as_mut() {
                None => bb = Some(b),
                Some(first) => {
                    first.solute_blocks.extend(b.solute_blocks);
                    first.end_blocks.extend(b.end_blocks);
                    first.solvent_blocks.extend(b.solvent_blocks);
                },
            },
            Err(e) => {
                eprintln!("Error reading MTB file '{}': {:?}", path, e);
                process::exit(1);
            },
        }
    }
    let Some(bb) = bb else {
        eprintln!("Error: at least one @build file is required");
        process::exit(1);
    };
    let ffp = match ifp::read_ifp_file(&args.param) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Error reading IFP file '{}': {:?}", args.param, e);
            process::exit(1);
        },
    };
    eprintln!("Building topology for {} residues", args.seq.len());
    eprintln!("Force field: {}", ffp.force_field);

    // gromos++ make_top.cc: the residue loop with its status machine
    //   0 solute, 1 begin group, 2 first solute after a begin group, 3 end group
    let mut lt = Lt::default();
    let mut status = 0;
    let mut resnum = 0usize;
    let mut repforward = 0i64;
    let mut first_atom = 0i64;
    let fail = |e: String| -> ! {
        eprintln!("Error: {e}");
        process::exit(1);
    };
    for (seq_idx, res_name) in args.seq.iter().enumerate() {
        let block = match bb.find_block(res_name) {
            Some(BlockRef::End(e)) => {
                status = if e.n_replace < 0 { 3 } else { 1 };
                from_end(e, &ffp)
            },
            Some(BlockRef::Solute(s)) => {
                status = if status == 1 { 2 } else { 0 };
                from_solute(s, &ffp)
            },
            None => fail(format!(
                "building block '{}' not found in the MTB file(s) (sequence position {})",
                res_name,
                seq_idx + 1
            )),
        };
        match status {
            0 => {
                add_solute(&mut lt, &block, resnum, res_name, 0, first_atom)
                    .unwrap_or_else(|e| fail(e));
                resnum += 1;
            },
            1 => {
                repforward = add_begin(&mut lt, &block, resnum);
                first_atom = lt.atoms.len() as i64 - block.atoms.len() as i64;
                add_cov_end(&mut lt, &block, first_atom).unwrap_or_else(|e| fail(e));
            },
            2 => {
                add_solute(&mut lt, &block, resnum, res_name, repforward, first_atom)
                    .unwrap_or_else(|e| fail(e));
                lt.remove_atoms();
                resnum += 1;
            },
            _ => {
                if resnum == 0 {
                    fail(format!(
                        "end group {res_name} cannot be at the start of the sequence"
                    ));
                }
                resnum -= 1;
                add_end(&mut lt, &block, resnum).unwrap_or_else(|e| fail(e));
                let off = lt.atoms.len() as i64 - block.atoms.len() as i64;
                add_cov_end(&mut lt, &block, off).unwrap_or_else(|e| fail(e));
                resnum += 1;
            },
        }
    }
    lt.get14s();
    let n = lt.atoms.len() as i64;
    for (i, a) in lt.atoms.iter_mut().enumerate() {
        let before = a.exclusions.len();
        a.exclusions.retain(|&e| e >= 0 && e < n);
        if a.exclusions.len() != before {
            eprintln!(
                "WARNING: exclusion skipped since it is not within the solute (atom {})",
                i + 1
            );
        }
    }
    write_topology(&lt, &bb, &ffp, solvent_name);
}

/// Molecules from bond connectivity (union–find), as gromos++ `parse_molecules` does.
fn molecules(lt: &Lt) -> Vec<usize> {
    let n = lt.atoms.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for b in &lt.bonds {
        if b.0 >= 0 && b.1 >= 0 && (b.0 as usize) < n && (b.1 as usize) < n {
            let (a, c) = (
                find(&mut parent, b.0 as usize),
                find(&mut parent, b.1 as usize),
            );
            if a != c {
                parent[a.max(c)] = a.min(c);
            }
        }
    }
    // last atom of each molecule, molecules being contiguous in a linear topology
    let mut last_atoms = Vec::new();
    for i in 0..n {
        let next_same = i + 1 < n && find(&mut parent, i + 1) == find(&mut parent, i);
        if !next_same {
            last_atoms.push(i + 1);
        }
    }
    last_atoms
}

fn print_list(prefix_width: usize, items: &[i64]) {
    // gromos++ OutTopology: six per line, continuation lines indented
    for (k, x) in items.iter().enumerate() {
        if k % 6 == 0 && k != 0 {
            println!();
            print!("{:width$}", "", width = prefix_width);
        }
        print!(" {:5}", x);
    }
    println!();
}

fn write_topology(
    lt: &Lt,
    bb: &BuildingBlocks,
    ffp: &ForceFieldParameters,
    solvent_name: Option<&str>,
) {
    println!("TITLE");
    println!("MAKE_TOP topology (gromos-rs), using:");
    println!("Force-field code: {}", ffp.force_field);
    println!("END");
    println!("PHYSICALCONSTANTS");
    println!("# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)");
    println!("{:15.4}", bb.fpepsi);
    println!("# HBAR: Planck's constant HBAR = H/(2* PI)");
    println!("{:15.7}", bb.hbar);
    println!("# SPDL: Speed of light (nm/ps)");
    println!("{:15.3}", bb.spdl);
    println!("# BOLTZ: Boltzmann's constant kB");
    println!("{:15.8}", bb.boltz);
    println!("END");
    println!("TOPVERSION");
    println!("2.0");
    println!("END");
    println!("ATOMTYPENAME");
    println!("# NRATT: number of van der Waals atom types");
    println!("{:5}", ffp.num_atom_types);
    println!("# TYPE: atom type names");
    for alj in &ffp.atom_lj {
        println!("{}", alj.name);
    }
    println!("END");
    println!("RESNAME");
    println!("# NRAA2: number of residues in a solute molecule");
    println!("{:5}", lt.res_names.len());
    println!("# AANM: residue names");
    for name in &lt.res_names {
        println!("{}", name);
    }
    println!("END");
    println!("SOLUTEATOM");
    println!("#   NRP: number of solute atoms");
    println!("{:5}", lt.atoms.len());
    println!("#  ATNM: atom number");
    println!("#  MRES: residue number");
    println!("#  PANM: atom name of solute atom");
    println!("#   IAC: integer (van der Waals) atom type code");
    println!("#  MASS: mass of solute atom");
    println!("#    CG: charge of solute atom");
    println!("#   CGC: charge group code (0 or 1)");
    println!("#   INE: number of excluded atoms");
    println!("# INE14: number of 1-4 interactions");
    println!("# ATNM MRES PANM IAC     MASS       CG  CGC INE");
    println!("#                                           INE14");
    for (i, a) in lt.atoms.iter().enumerate() {
        print!(
            "{:6} {:4} {:4} {:3} {:8.5} {:8.5} {:2} {:5}",
            i + 1,
            a.res + 1,
            a.name,
            a.iac + 1,
            a.mass,
            a.charge,
            a.chargegroup,
            a.exclusions.len()
        );
        let excl: Vec<i64> = a.exclusions.iter().map(|e| e + 1).collect();
        print_list(47, &excl);
        print!("{:46}{}", "", a.one_four.len());
        let one_four: Vec<i64> = a.one_four.iter().map(|&e| e as i64 + 1).collect();
        print_list(45, &one_four);
    }
    println!("END");
    println!("BONDSTRETCHTYPE");
    println!("#  NBTY: number of covalent bond types");
    println!("{:5}", ffp.bond_types.len());
    println!("#  CB: force constant");
    println!("#  HB: harmonic force constant");
    println!("#  B0: bond length at minimum energy");
    println!("#         CB         HB         B0");
    for bt in &ffp.bond_types {
        println!(
            "{:15.7e}{:15.7e}{:13.7e}",
            bt.k_quartic, bt.k_harmonic, bt.r0
        );
    }
    println!("END");
    let atoms = &lt.atoms;
    let h = |i: i64| is_h(atoms[i as usize].mass);
    let (bondh, bond): (Vec<&Bond>, Vec<&Bond>) = lt.bonds.iter().partition(|b| h(b.0) || h(b.1));
    println!("BONDH");
    println!("#  NBONH: number of bonds involving H atoms in solute");
    println!("{:5}", bondh.len());
    println!("#  IBH, JBH: atom sequence numbers of atoms forming a bond");
    println!("#  ICBH: bond type code");
    println!("#   IBH    JBH ICBH");
    for b in &bondh {
        println!("{:7}{:7}{:5}", b.0 + 1, b.1 + 1, b.2);
    }
    println!("END");
    println!("BOND");
    println!("#  NBON: number of bonds NOT involving H atoms in solute");
    println!("{:5}", bond.len());
    println!("#  IB, JB: atom sequence numbers of atoms forming a bond");
    println!("#  ICB: bond type code");
    println!("#    IB     JB  ICB");
    for b in &bond {
        println!("{:7}{:7}{:5}", b.0 + 1, b.1 + 1, b.2);
    }
    println!("END");
    println!("BONDANGLEBENDTYPE");
    println!("#  NTTY: number of bond angle types");
    println!("{:5}", ffp.angle_types.len());
    println!("#  CT: force constant (based on potential");
    println!("#      harmonic in the angle cosine)");
    println!("#  CHT: force constant (based on potential");
    println!("#      harmonic in the angle)");
    println!("#  T0: bond angle at minimum energy in degrees");
    println!("#         CT         CHT          T0");
    for at in &ffp.angle_types {
        println!(
            "{:15.7e}{:15.7e}{:13.7e}",
            at.k_non_harmonic, at.k_harmonic, at.theta0
        );
    }
    println!("END");
    let (angleh, angle): (Vec<&AngleKey>, Vec<&AngleKey>) =
        lt.angles.iter().partition(|a| h(a.0) || h(a.1) || h(a.2));
    println!("BONDANGLEH");
    println!("#  NTHEH: number of bond angles involving H atoms in solute");
    println!("{:5}", angleh.len());
    println!("#  ITH, JTH, KTH: atom sequence numbers");
    println!("#    of atoms forming a bond angle in solute");
    println!("#  ICTH: bond angle type code");
    println!("#   ITH    JTH    KTH ICTH");
    for a in &angleh {
        let (i, j, k, t) = angle_atoms(a);
        println!("{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, t);
    }
    println!("END");
    println!("BONDANGLE");
    println!("#  NTHE: number of bond angles NOT");
    println!("#        involving H atoms in solute");
    println!("{:5}", angle.len());
    println!("#  IT, JT, KT: atom sequence numbers of atoms");
    println!("#     forming a bond angle");
    println!("#  ICT: bond angle type code");
    println!("#    IT     JT     KT  ICT");
    for a in &angle {
        let (i, j, k, t) = angle_atoms(a);
        println!("{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, t);
    }
    println!("END");
    println!("IMPDIHEDRALTYPE");
    println!("#  NQTY: number of improper dihedrals");
    println!("{:5}", ffp.improper_types.len());
    println!("#  CQ: force constant of improper dihedral per degrees square");
    println!("#  Q0: improper dihedral angle at minimum energy in degrees");
    println!("#         CQ          Q0");
    for it in &ffp.improper_types {
        println!("{:15.7e}{:15.7e}", it.k, it.xi0);
    }
    println!("END");
    let (imph, imp): (Vec<&Improper>, Vec<&Improper>) = lt
        .impropers
        .iter()
        .partition(|p| h(p.0) || h(p.1) || h(p.2) || h(p.3));
    println!("IMPDIHEDRALH");
    println!("#  NQHIH: number of improper dihedrals");
    println!("#         involving H atoms in the solute");
    println!("{:5}", imph.len());
    println!("#  IQH,JQH,KQH,LQH: atom sequence numbers");
    println!("#     of atoms forming an improper dihedral");
    println!("#  ICQH: improper dihedral type code");
    println!("#   IQH    JQH    KQH    LQH ICQH");
    for p in &imph {
        println!(
            "{:7}{:7}{:7}{:7}{:5}",
            p.0 + 1,
            p.1 + 1,
            p.2 + 1,
            p.3 + 1,
            p.4
        );
    }
    println!("END");
    println!("IMPDIHEDRAL");
    println!("#  NQHI: number of improper dihedrals NOT");
    println!("#    involving H atoms in solute");
    println!("{:5}", imp.len());
    println!("#  IQ,JQ,KQ,LQ: atom sequence numbers of atoms");
    println!("#    forming an improper dihedral");
    println!("#  ICQ: improper dihedral type code");
    println!("#    IQ     JQ     KQ     LQ  ICQ");
    for p in &imp {
        println!(
            "{:7}{:7}{:7}{:7}{:5}",
            p.0 + 1,
            p.1 + 1,
            p.2 + 1,
            p.3 + 1,
            p.4
        );
    }
    println!("END");
    println!("TORSDIHEDRALTYPE");
    println!("#  NPTY: number of dihedral types");
    println!("{:5}", ffp.torsion_types.len());
    println!("#  CP: force constant");
    println!("#  PD: phase-shift angle");
    println!("#  NP: multiplicity");
    println!("#       CP        PD  NP");
    for tt in &ffp.torsion_types {
        println!("{:10.5}{:10.5}{:4}", tt.k, tt.delta, tt.multiplicity);
    }
    println!("END");
    let (dihh, dih): (Vec<&DihedralKey>, Vec<&DihedralKey>) = lt.dihedrals.iter().partition(|d| {
        let (i, j, k, l, _) = dihedral_atoms(d);
        h(i) || h(j) || h(k) || h(l)
    });
    println!("DIHEDRALH");
    println!("#  NPHIH: number of dihedrals involving H atoms in solute");
    println!("{:5}", dihh.len());
    println!("#  IPH, JPH, KPH, LPH: atom sequence numbers");
    println!("#    of atoms forming a dihedral");
    println!("#  ICPH: dihedral type code");
    println!("#   IPH    JPH    KPH    LPH ICPH");
    for d in &dihh {
        let (i, j, k, l, t) = dihedral_atoms(d);
        println!("{:7}{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, l + 1, t);
    }
    println!("END");
    println!("DIHEDRAL");
    println!("#  NPHI: number of dihedrals NOT involving H atoms in solute");
    println!("{:5}", dih.len());
    println!("#  IP, JP, KP, LP: atom sequence numbers");
    println!("#     of atoms forming a dihedral");
    println!("#  ICP: dihedral type code");
    println!("#    IP     JP     KP     LP  ICP");
    for d in &dih {
        let (i, j, k, l, t) = dihedral_atoms(d);
        println!("{:7}{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, l + 1, t);
    }
    println!("END");
    println!("CROSSDIHEDRALH");
    println!("#  NPHIH: number of cross dihedrals involving H atoms in solute");
    println!("    0");
    println!("END");
    println!("CROSSDIHEDRAL");
    println!("#  NPPC: number of cross dihedrals NOT involving H atoms in solute");
    println!("    0");
    println!("END");
    println!("LJPARAMETERS");
    println!("#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2");
    let n_types = ffp.num_atom_types;
    println!("{:5}", n_types * (n_types + 1) / 2);
    println!("#  IAC,JAC: integer (van der Waals) atom type code");
    println!("#  C12: r**(-12) term in nonbonded interactions");
    println!("#   C6: r**(-6) term in nonbonded interactions");
    println!("# CS12: r**(-12) term in 1-4 nonbonded interactions");
    println!("#  CS6: r**(-6) term in 1-4 nonbonded interactions");
    println!("# IAC  JAC           C12            C6          CS12           CS6");
    for j in 1..=n_types {
        for i in 1..=j {
            let c6 = ffp.compute_c6(i, j);
            let c12 = ffp.compute_c12(i, j);
            let cs6 = ffp.compute_c6_14(i, j);
            let cs12 = ffp.compute_c12_14(i, j);
            println!(
                "{:5}{:5}{:15.6e}{:15.6e}{:15.6e}{:15.6e}",
                i, j, c12, c6, cs12, cs6
            );
        }
    }
    println!("END");
    let mols = molecules(lt);
    println!("SOLUTEMOLECULES");
    println!("# NSPM: number of separate molecules in solute block");
    println!("# NSP[1...NSPM]: atom sequence number of last atom");
    println!("#                of the successive submolecules");
    println!("#      NSPM  NSP[1...NSPM]");
    println!("{:5}", mols.len());
    let mols_i: Vec<i64> = mols.iter().map(|&m| m as i64).collect();
    for chunk in mols_i.chunks(10) {
        println!(
            "{}",
            chunk.iter().map(|m| format!("{:6}", m)).collect::<String>()
        );
    }
    println!("END");
    for (block, comment) in [
        ("TEMPERATUREGROUPS", "temperature"),
        ("PRESSUREGROUPS", "pressure"),
    ] {
        println!("{block}");
        println!("# NSTM: number of {comment} atom groups");
        println!("# NST[1...NSTM]: atom sequence number of last atom");
        println!("#                of the successive {comment} atom groups");
        println!("#      NSTM  NST[1...NSTM]");
        println!("    1");
        println!("{:6}", lt.atoms.len());
        println!("END");
    }
    println!("LJEXCEPTIONS");
    println!("# This block defines special LJ-interactions based on atom numbers");
    println!("# This overrules the normal LJ-parameters (including 1-4s)");
    println!("# NEX: number of exceptions");
    println!("    0");
    println!("END");
    if let Some(sname) = solvent_name {
        let Some(sol) = bb.find_solvent(sname) else {
            eprintln!("Error: solvent building block '{}' not found", sname);
            process::exit(1);
        };
        println!("SOLVENTATOM");
        println!("#  NRAM: number of atoms per solvent molecule");
        println!("{:5}", sol.atoms.len());
        println!("#     I: solvent atom sequence number");
        println!("#  IACS: integer (van der Waals) atom type code");
        println!("#  ANMS: atom name of solvent atom");
        println!("#  MASS: mass of solvent atom");
        println!("#   CGS: charge of solvent atom");
        println!("#  I  ANMS IACS      MASS        CGS");
        for (idx, (aname, iac, mass_code, charge)) in sol.atoms.iter().enumerate() {
            let mass = ffp.get_mass(*mass_code).unwrap_or(0.0);
            println!(
                "{:5} {:>5}{:5}{:11.5}{:11.5}",
                idx + 1,
                aname,
                iac,
                mass,
                charge
            );
        }
        println!("END");
        println!("SOLVENTCONSTR");
        println!("#  NCONS: number of constraints");
        println!("{:5}", sol.constraints.len());
        println!("#  ICONS, JCONS: atom sequence numbers forming constraint");
        println!("#   CONS constraint length");
        println!("#ICONS JCONS         CONS");
        for (ic, jc, dist) in &sol.constraints {
            println!("{:5}{:5}{:15.7}", ic, jc, dist);
        }
        println!("END");
    }
}
