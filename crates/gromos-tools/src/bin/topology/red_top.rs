//! red_top — reduce a topology to a subset of its solute atoms (gromos++ `red_top`).
//!
//! Usage: red_top @topo <file> @atoms <first>-<last>   (1-based, inclusive)
//!    or: red_top @topo <file> @natoms <n>              (the first n atoms)
//!
//! Bonded terms, exclusions and 1-4 pairs are kept when every atom they reference stays;
//! residues, molecules and the temperature/pressure groups are renumbered; the parameter
//! blocks and the solvent are carried over unchanged. Output in gromos++ `OutTopology` layout.

use gromos_io::topology::{read_topology_file, write_parsed_topology, ParsedTopology};
use std::env;
use std::process;

fn usage() -> ! {
    eprintln!("Usage: red_top @topo <file> @atoms <first>-<last> | @natoms <n>");
    process::exit(1);
}

/// The atoms `first..=last` (1-based) of `topo` as a topology of their own.
pub fn reduce(topo: &ParsedTopology, first: usize, last: usize) -> ParsedTopology {
    let keep = |a: usize| (first - 1..last).contains(&a);
    let map = |a: usize| a + 1 - first;
    let range = (first - 1)..last;
    let mut out = topo.clone();
    out.n_atoms = range.len();
    out.atom_names = topo.atom_names[range.clone()].to_vec();
    out.masses = topo.masses[range.clone()].to_vec();
    out.charges = topo.charges[range.clone()].to_vec();
    out.iac = topo.iac[range.clone()].to_vec();
    out.chargegroup_codes = topo.chargegroup_codes[range.clone()].to_vec();
    // residues: renumbered from 1, names of the kept residues only
    let kept_res: Vec<usize> = topo.residue_numbers[range.clone()].to_vec();
    let mut res_map: Vec<usize> = Vec::new();
    for &r in &kept_res {
        if !res_map.contains(&r) {
            res_map.push(r);
        }
    }
    out.residue_numbers = kept_res
        .iter()
        .map(|r| res_map.iter().position(|x| x == r).unwrap() + 1)
        .collect();
    out.residue_names = res_map
        .iter()
        .map(|&r| topo.residue_names.get(r - 1).cloned().unwrap_or_default())
        .collect();
    out.exclusions = range
        .clone()
        .map(|a| {
            topo.exclusions
                .get(a)
                .map(|e| e.iter().filter(|&&x| keep(x)).map(|&x| map(x)).collect())
                .unwrap_or_default()
        })
        .collect();
    out.one_four_pairs = range
        .clone()
        .map(|a| {
            topo.one_four_pairs
                .get(a)
                .map(|e| e.iter().filter(|&&x| keep(x)).map(|&x| map(x)).collect())
                .unwrap_or_default()
        })
        .collect();
    out.bonds = topo
        .bonds
        .iter()
        .filter(|b| keep(b.0) && keep(b.1))
        .map(|b| (map(b.0), map(b.1), b.2))
        .collect();
    out.angles = topo
        .angles
        .iter()
        .filter(|a| keep(a.0) && keep(a.1) && keep(a.2))
        .map(|a| (map(a.0), map(a.1), map(a.2), a.3))
        .collect();
    out.improper_dihedrals = topo
        .improper_dihedrals
        .iter()
        .filter(|d| keep(d.0) && keep(d.1) && keep(d.2) && keep(d.3))
        .map(|d| (map(d.0), map(d.1), map(d.2), map(d.3), d.4))
        .collect();
    out.proper_dihedrals = topo
        .proper_dihedrals
        .iter()
        .filter(|d| keep(d.0) && keep(d.1) && keep(d.2) && keep(d.3))
        .map(|d| (map(d.0), map(d.1), map(d.2), map(d.3), d.4))
        .collect();
    // molecules: the kept part of each original molecule
    let mut mols = Vec::new();
    for &end in &topo.solute_molecules {
        if end >= first && end <= last {
            mols.push(end + 1 - first);
        }
    }
    if mols.last() != Some(&out.n_atoms) && out.n_atoms > 0 {
        mols.push(out.n_atoms);
    }
    out.solute_molecules = mols;
    out.temperature_groups = vec![out.n_atoms - 1];
    out.pressure_groups = vec![out.n_atoms];
    out
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut topo_file = None;
    let mut first = 1usize;
    let mut last = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = args.get(i).cloned();
            },
            "@natoms" => {
                i += 1;
                last = args.get(i).and_then(|s| s.parse().ok());
            },
            "@atoms" => {
                i += 1;
                let spec = args.get(i).cloned().unwrap_or_default();
                // gromos++ "1:3-9" (molecule:range) or plain "3-9"
                let spec = spec.rsplit(':').next().unwrap_or("").to_string();
                let (a, b) = spec.split_once('-').unwrap_or((&spec, &spec));
                first = a.trim().parse().unwrap_or_else(|_| usage());
                last = Some(b.trim().parse().unwrap_or_else(|_| usage()));
            },
            _ => {},
        }
        i += 1;
    }
    let (Some(topo_file), Some(last)) = (topo_file, last) else {
        usage()
    };
    let topo = match read_topology_file(&topo_file) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error reading topology '{}': {}", topo_file, e);
            process::exit(1);
        },
    };
    if first == 0 || first > last || last > topo.n_atoms {
        eprintln!("Error: atoms {first}-{last} outside 1-{}", topo.n_atoms);
        process::exit(1);
    }
    eprintln!(
        "Original topology: {} atoms; keeping {}-{}",
        topo.n_atoms, first, last
    );
    let reduced = reduce(&topo, first, last);
    print!(
        "{}",
        write_parsed_topology(
            &reduced,
            &format!("red_top: atoms {first}-{last} of {topo_file}")
        )
    );
}
